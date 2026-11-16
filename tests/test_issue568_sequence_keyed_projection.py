"""Issue #568: reach Route A interval members by sequence, not by family name.

Route A proposes a reference gene in whatever species has the literature
(usually Arabidopsis); the interval index names its subjects after whatever the
pinned reference proteome contains (usually a legume). BOP1 and NOOT1, TFL1 and
CEN, or a bare TrEMBL accession never meet in a name space, so a name-only
projection loses verified baits. The frozen blast-panel-search evidence links
the two sides without either name.

The projection stays positive-only: no hit is not absence, the genome-wide
paralog guard stays unevaluated, and every threshold-passing copy is kept.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main

# Route A speaks Arabidopsis (BOP1); the interval index speaks Medicago
# (NOOT1). Deliberately no shared token between the two sides.
HYPOTHESES = """hypothesis_id\tgene_family\treference_gene_name\tproposed_stable_id\tmember_interval_genes
H_BOP\tBTB-ankyrin transcriptional co-activators of the BOP clade\tBOP1\tQ9M1I7\t
"""
STATUS = """hypothesis_id\tfamily_status\tinterval_gene_ids\tsource_ids
H_BOP\tnot_evaluated_missing_reference_genome_family_index\t\tPMID:123
"""
BAIT = """hypothesis_id\tfunction_source_ids\tsource_provider\tbait_id
H_BOP\tPMID:123\tsource_verified\tQ9M1I7
"""
INDEX = """gene_id\tassignment_status\tindex_scope\tbest_subject\tbest_subject_symbol\tbest_subject_stem\tfamily_symbols\tfamily_stems\tdiscovery_bucket_id\tevolutionary_family_id\tbest_identity\tbest_query_coverage\tbest_evalue
Ps2\tassigned\tinterval_candidates\tsp|Q2HW56|NOOT1_MEDTR\tNOOT1\tNOOT\tNOOT1\tNOOT\tNOOT\t\t71.2\t0.93\t1e-120
Ps4\tassigned\tinterval_candidates\tsp|Q9|PAL_ARATH\tPAL1\tPAL\tPAL1\tPAL\tPAL\t\t72.0\t0.91\t1e-100
"""
GENES = "gene_id\nPs1\nPs2\nPs3\nPs4\n"
BLAST = ("panel_id\tprovider_status\treview_status\tambiguity_status\ttop_hit_target_id\t"
         "top_hit_identity\ttop_hit_query_coverage\ttop_hit_target_coverage\ttop_hit_evalue\t"
         "paralog_count\tdb_fingerprint\tinterval_proteins_sha256\tblastp_version\n"
         "Q9M1I7\tevaluable\tverified\t{ambiguity}\tPs2\t62.5\t0.88\t0.9\t1e-90\t{paralogs}\t"
         "db1\tsha1\t2.12.0+\n")


def _long_form(rows: list[tuple[str, str, str]]) -> str:
    header = ("panel_id\ttarget_id\tidentity\tquery_coverage\ttarget_coverage\tevalue\tbitscore\t"
              "review_status\tprovider_status\trank_contribution\teligibility\n")
    body = "".join(
        f"Q9M1I7\t{target}\t{identity}\t{coverage}\t0.9\t1e-90\t400\tverified\tevaluable\tnone\teligible\n"
        for target, identity, coverage in rows)
    return header + body


def _write(root: Path, name: str, text: str) -> Path:
    path = root / name
    path.write_text(text, encoding="utf-8")
    return path


class SequenceKeyedProjectionTests(unittest.TestCase):
    def _inputs(self, root: Path, *, long_form_rows=None, ambiguity="ambiguous_paralog", paralogs="1"):
        files = {
            "hypotheses": _write(root, "hypotheses.tsv", HYPOTHESES),
            "status": _write(root, "status.tsv", STATUS),
            "bait": _write(root, "bait.tsv", BAIT),
            "genes": _write(root, "genes.tsv", GENES),
            "blast": _write(root, "blast_evidence.tsv",
                            BLAST.format(ambiguity=ambiguity, paralogs=paralogs)),
        }
        if long_form_rows is not None:
            files["eligible"] = _write(root, "eligible_interval_hits.tsv", _long_form(long_form_rows))
        return files

    def _run(self, root: Path, files, out_name="out", index=None):
        out = root / out_name
        argv = [
            "family-candidate-bridge", "--out-dir", str(out), "--case-id", "case568",
            "--route-a-hypotheses", str(files["hypotheses"]),
            "--route-a-family-status", str(files["status"]),
            "--route-a-bait-panel", str(files["bait"]),
            "--genes", str(files["genes"]),
            "--route-a-blast", str(files["blast"]),
        ]
        if "eligible" in files:
            argv += ["--route-a-blast-eligible", str(files["eligible"])]
        if index is not None:
            argv += ["--route-b-family-index", str(index)]
        self.assertEqual(0, main(argv))
        bridge = pd.read_csv(out / "family_candidate_bridge.tsv", sep="\t", dtype=str).fillna("")
        manifest = json.loads((out / "family_candidate_bridge_manifest.json").read_text(encoding="utf-8"))
        return bridge, manifest, out

    def test_verified_bait_reaches_interval_members_without_a_shared_name(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            files = self._inputs(root, long_form_rows=[("Ps2", "62.5", "0.88"), ("Ps3", "48.0", "0.71")])
            bridge, manifest, out = self._run(root, files)

            route_a = bridge[bridge["route"] == "a"]
            self.assertEqual(["Ps2", "Ps3"], sorted(route_a["gene_id"]))
            self.assertEqual({"blastp"}, set(route_a["match_type"]))
            self.assertEqual({"none"}, set(route_a["rank_contribution"]))
            projection = manifest["route_a_sequence_family_projection"]
            self.assertEqual(2, projection["projection_count"])
            self.assertFalse(projection["paralog_guard_evaluated"])
            self.assertTrue(all(not row["paralog_guard_evaluated"] for row in projection["projections"]))
            self.assertEqual({"not_evaluated_missing_reference_genome_family_index"},
                             {row["family_status"] for row in projection["projections"]})
            # The falsified-family ledger is a different claim and stays empty.
            self.assertEqual(0, manifest["route_a_blast"]["falsified_family_projection_count"])

            union = root / "union"
            self.assertEqual(0, main([
                "supported-member-union", "--out-dir", str(union), "--case-id", "case568",
                "--route-a-bridge", str(out / "family_candidate_bridge.tsv"),
            ]))
            supported = pd.read_csv(union / "supported_member_union.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(["Ps2", "Ps3"], sorted(supported["candidate_gene_id"]))

    def test_every_threshold_passing_copy_is_retained(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            files = self._inputs(root, long_form_rows=[
                ("Ps2", "62.5", "0.88"), ("Ps3", "61.0", "0.85"), ("Ps4", "59.0", "0.80")])
            bridge, _, _ = self._run(root, files)
            route_a = bridge[bridge["route"] == "a"]
            self.assertEqual(["Ps2", "Ps3", "Ps4"], sorted(route_a["gene_id"]))
            self.assertEqual({"paralog_ambiguous"}, set(route_a["ambiguity_status"]))

    def test_no_hit_is_not_absence(self):
        """No eligible interval target: no promoted row, and no negative row."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            files = self._inputs(root, long_form_rows=[("Outside1", "62.5", "0.88")])
            bridge, manifest, _ = self._run(root, files)
            promoted = bridge[bridge["mapping_tier"].isin(
                ["tier3_sequence_match_unique", "tier4_sequence_match_ambiguous"])]
            self.assertEqual(0, len(promoted))
            self.assertEqual(0, manifest["route_a_sequence_family_projection"]["projection_count"])
            self.assertNotIn("falsified", " ".join(bridge["reason"].tolist()))

    def test_name_and_sequence_projections_do_not_duplicate_a_gene(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index = _write(root, "family_index.tsv", INDEX.replace("NOOT1\tNOOT\tNOOT1", "BOP1\tBOP\tBOP1"))
            files = self._inputs(root, long_form_rows=[("Ps2", "62.5", "0.88")])
            bridge, manifest, _ = self._run(root, files, index=index)
            route_a = bridge[bridge["route"] == "a"]
            self.assertEqual(["Ps2"], sorted(route_a["gene_id"]))
            self.assertEqual(1, len(route_a))
            self.assertEqual(
                ["Ps2"],
                manifest["route_a_interval_family_projection"]["deduplicated_into_sequence_projection"])

    def test_unverified_bait_projects_nothing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            files = self._inputs(root, long_form_rows=[("Ps2", "62.5", "0.88")])
            files["bait"] = _write(root, "bait.tsv", "hypothesis_id\tfunction_source_ids\tsource_provider\tbait_id\n")
            bridge, manifest, _ = self._run(root, files, out_name="unverified")
            self.assertEqual(0, manifest["route_a_sequence_family_projection"]["projection_count"])
            self.assertEqual(0, len(bridge[bridge["match_type"] == "blastp"]))


if __name__ == "__main__":
    unittest.main()
