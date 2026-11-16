"""Issue #599: zero-resolution family labels stay evaluable holds, never falsified.

A payload-declared gene_family / reference_gene_name / proposed_stable_id that
matches nothing anywhere in the family index is an index-vocabulary gap
(literature-standard names are often stored under different BLAST-symbol
stems), so it must be recorded as ``not_evaluated_family_unresolved_in_index``
-- never as a ``falsified_*`` biological absence claim. The downstream wet
gate still fails closed.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.family_candidate_bridge import BRIDGE_COLUMNS
from candidate_qtl.schemas import FAMILY_CANDIDATE_BRIDGE_SCHEMA, validate


def _hit13(query: str, subject: str, identity: float, bitscore: float, title: str, *, qlen: int = 300) -> str:
    # outfmt 6 with a stitle (13 columns): ... evalue bitscore stitle
    return f"{query}\t{subject}\t{identity}\t{qlen}\t{qlen}\t{qlen}\t1\t{qlen}\t1\t{qlen}\t1e-60\t{bitscore}\t{title}"


PANEL_HEADER = (
    "case_id\tbait_id\thypothesis_id\tpathway_or_process\tgene_family\t"
    "reference_gene_name\treference_species\texpected_role\tsequence_length\t"
    "sequence_sha256\tsource_provider\tsource_url\tretrieved_at\tfunction_source_ids\t"
    "rank_contribution\tpanel_status"
)

HYPOTHESES_ROUTE_A = """hypothesis_id\tpathway_or_process\tgene_family\treference_gene_name\tproposed_stable_id
H_NFR\tpea nodulation signalling\tNFR1\tNFR1\tpanel_nfr1
"""

# Bridge-shaped hypotheses table (same columns as the #530/#532 fixtures).
HYPOTHESES_BRIDGE = """hypothesis_id\tgene_family\treference_gene_name\tproposed_stable_id\tmember_interval_genes
H_NFR\tNFR1\tNFR1\tpanel_nfr1\tG_A;G_B
"""

BAIT_PANEL = """hypothesis_id\tfunction_source_ids\tsource_provider
H_NFR\tPMID:14534591\tlocal_verified_panel
"""

# Hand-written status row in the exact shape route-a-interval-match writes for
# a zero-resolution hypothesis: no interval genes, evaluable hold.
FAMILY_STATUS = """hypothesis_id\tfamily_status\tinterval_gene_ids\tsource_ids\treference_gene_name\tproposed_stable_id
H_NFR\tnot_evaluated_family_unresolved_in_index\t\tPMID:14534591\tNFR1\tpanel_nfr1
"""

INTERVAL_GENES = """gene_id\tchrom\tstart\tend
G_A\tchr1\t1000\t2000
G_B\tchr1\t5000\t6000
"""


def _write(root: Path, name: str, text: str) -> Path:
    path = root / name
    path.write_text(text, encoding="utf-8")
    return path


class Issue599ZeroResolutionHoldTests(unittest.TestCase):
    def _build_route_a_inputs(self, root: Path) -> tuple[Path, Path, Path]:
        """An index that only knows CERK/LYK stems while the verified payload
        declares the literature-standard name NFR1 (the #599 laneG shape)."""
        hits = root / "hits.tsv"
        hits.write_text("\n".join([
            _hit13("G_A", "sp|Q9M2A0|CERK1_ARATH", 82.0, 520,
                   "Cerberus kinase OS=Arabidopsis thaliana GN=CERK1 PE=1 SV=1"),
            _hit13("G_B", "tr|B9DGA0|LYK5_MEDTR", 78.0, 480,
                   "LysM receptor kinase OS=Medicago truncatula GN=LYK5 PE=3 SV=1"),
        ]) + "\n", encoding="utf-8")
        index_out = root / "index"
        self.assertEqual(0, main(["family-index", "--blast-hits", str(hits),
                                  "--genes", str(_write(root, "genes.tsv", INTERVAL_GENES)),
                                  "--out-dir", str(index_out), "--reference-label", "TestRef"]))
        hypotheses = _write(root, "hypotheses.tsv", HYPOTHESES_ROUTE_A)
        rows = [PANEL_HEADER,
                "\t".join(["case599", "panel_nfr1", "H_NFR", "", "NFR1", "NFR1",
                           "Pisum sativum", "", "96", "sha256nfr", "local_record_cache", "",
                           "2026-01-01T00:00:00Z", "PMID:14534591", "none", "frozen"])]
        panel = _write(root, "verified_bait_panel.tsv", "\n".join(rows) + "\n")
        return index_out / "family_index.tsv", hypotheses, panel

    def _run_bridge(self, root: Path, out_name: str = "bridge") -> Path:
        out = root / out_name
        args = [
            "family-candidate-bridge", "--out-dir", str(out), "--case-id", "case599",
            "--route-a-hypotheses", str(_write(root, "hypotheses_bridge.tsv", HYPOTHESES_BRIDGE)),
            "--route-a-family-status", str(_write(root, "status.tsv", FAMILY_STATUS)),
            "--route-a-bait-panel", str(_write(root, "bait.tsv", BAIT_PANEL)),
            "--genes", str(_write(root, "genes.tsv", INTERVAL_GENES)),
        ]
        self.assertEqual(0, main(args))
        return out

    def test_zero_resolution_route_a_status_is_not_evaluated_never_falsified(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index, hypotheses, panel = self._build_route_a_inputs(root)
            out = root / "match"
            self.assertEqual(0, main(["route-a-interval-match", "--panel", str(panel),
                                      "--hypotheses", str(hypotheses), "--out-dir", str(out),
                                      "--family-index", str(index), "--interval", "chr1:1-2000000"]))
            frame = pd.read_csv(out / "route_a_interval_family_match.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(1, len(frame))
            row = frame.iloc[0]
            self.assertEqual("not_evaluated_family_unresolved_in_index", row["verdict"])
            self.assertEqual("not_evaluated_family_unresolved_in_index", row["family_status"])
            self.assertEqual("none", row["match_level"])
            self.assertEqual("True", row["evaluable"])
            self.assertEqual("True", row["stage_completed"])
            manifest = json.loads((out / "route_a_interval_family_match_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(manifest["evaluated"])
            self.assertIsNone(manifest["not_evaluated_reason"])
            self.assertEqual({"not_evaluated_family_unresolved_in_index": 1},
                             manifest["family_status_counts"])
            falsification = pd.read_csv(out / "hypothesis_interval_falsification.tsv", sep="\t",
                                        dtype=str).fillna("")
            coarse = falsification[falsification.level == "coarse"].iloc[0]
            fine = falsification[falsification.level == "fine"].iloc[0]
            self.assertEqual("not_evaluated_family_unresolved_in_index", coarse["verdict"])
            self.assertEqual("", coarse["resolved_by"])
            self.assertEqual("not_computed_coarse_family_unresolved_in_index", fine["verdict"])
            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            review = reviews[reviews.reason == "hypothesis_family_not_detected_genome_wide"].iloc[0]
            self.assertEqual("not_evaluated_family_unresolved_in_index", review.status)
            manifest_falsify = json.loads(
                (out / "pathway_interval_falsify_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(0, manifest_falsify["coarse_falsified"])
            self.assertEqual(1, manifest_falsify["coarse_family_unresolved_in_index"])

    def test_zero_resolution_bridge_rows_are_tier6_with_explicit_reason(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = self._run_bridge(root)
            bridge = pd.read_csv(out / "family_candidate_bridge.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(BRIDGE_COLUMNS, list(bridge.columns))
            self.assertEqual(["G_A", "G_B"], sorted(bridge["gene_id"].tolist()))
            for _, row in bridge.iterrows():
                self.assertEqual("tier6_not_evaluated", row["mapping_tier"])
                self.assertEqual("not_evaluated", row["ambiguity_status"])
                self.assertEqual("none", row["rank_contribution"])
                self.assertEqual("not_evaluated_family_unresolved_in_index", row["reason"])
            validate(bridge, FAMILY_CANDIDATE_BRIDGE_SCHEMA)

    def test_zero_resolution_never_becomes_supported(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bridge_out = self._run_bridge(root)
            union_out = root / "union"
            self.assertEqual(0, main([
                "supported-member-union", "--out-dir", str(union_out), "--case-id", "case599",
                "--route-a-bridge", str(bridge_out / "family_candidate_bridge.tsv"),
            ]))
            union = pd.read_csv(union_out / "supported_member_union.tsv", sep="\t", dtype=str).fillna("")
            if not union.empty:
                self.assertNotIn("supported", set(union["support_status"]))


if __name__ == "__main__":
    unittest.main()
