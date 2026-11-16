"""Issue #530: a verified Route A BLAST projection survives family falsification."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.canonical_family_recovery import build_canonical_family_recovery
from candidate_qtl.family_candidate_bridge import BRIDGE_COLUMNS
from candidate_qtl.family_to_copy import family_to_copy
from candidate_qtl.schemas import FAMILY_CANDIDATE_BRIDGE_SCHEMA, validate
from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff


HYPOTHESES = """hypothesis_id\tgene_family\treference_gene_name\tproposed_stable_id\tmember_interval_genes
H_FALSE\tFAM_FALSE\tREF_FALSE\tACC_FALSE\tG_UNSUPPORTED;G_OTHER
"""

BAIT_PANEL = """hypothesis_id\tfunction_source_ids\tsource_provider
H_FALSE\tROOT_FALSE\tlocal_verified_panel
"""

FAMILY_STATUS = """hypothesis_id\tfamily_status\tinterval_gene_ids\tsource_ids\treference_gene_name\tproposed_stable_id
H_FALSE\tfalsified_family_not_detected_genome_wide\t\tROOT_FALSE\tREF_FALSE\tACC_FALSE
"""

BLAST_EVIDENCE = """panel_id\tprovider_status\ttop_hit_target_id\ttop_hit_identity\ttop_hit_query_coverage\ttop_hit_target_coverage\ttop_hit_evalue\tparalog_count\tambiguity_status\treview_status\trank_contribution
ACC_FALSE\tevaluable\tG_BLAST\t84.924\t0.982625\t0.986792\t0.0\t0\tunique\tverified\tnone
"""

INTERVAL_GENES = """gene_id\tchrom\tstart\tend
G_BLAST\tchr1\t1000\t3000
G_OTHER\tchr1\t5000\t7000
"""


def _write(root: Path, name: str, text: str) -> Path:
    path = root / name
    path.write_text(text, encoding="utf-8")
    return path


class Issue530RouteABlastProjectionTests(unittest.TestCase):
    def _run_bridge(self, root: Path, blast: str = BLAST_EVIDENCE,
                    include_genes: bool = True) -> tuple[Path, pd.DataFrame]:
        inputs = {
            "hypotheses": _write(root, "hypotheses.tsv", HYPOTHESES),
            "bait": _write(root, "bait.tsv", BAIT_PANEL),
            "status": _write(root, "status.tsv", FAMILY_STATUS),
            "blast": _write(root, "blast.tsv", blast),
        }
        if include_genes:
            inputs["genes"] = _write(root, "genes.tsv", INTERVAL_GENES)
        out = root / "bridge"
        args = [
            "family-candidate-bridge", "--out-dir", str(out), "--case-id", "case530",
            "--route-a-hypotheses", str(inputs["hypotheses"]),
            "--route-a-family-status", str(inputs["status"]),
            "--route-a-bait-panel", str(inputs["bait"]),
            "--route-a-blast", str(inputs["blast"]),
        ]
        if include_genes:
            args.extend(["--genes", str(inputs["genes"])])
        self.assertEqual(0, main(args))
        return out, pd.read_csv(out / "family_candidate_bridge.tsv", sep="\t", dtype=str).fillna("")

    def test_falsified_family_emits_only_actual_verified_blast_target_and_union_supports_it(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out, bridge = self._run_bridge(root)
            self.assertEqual(BRIDGE_COLUMNS, list(bridge.columns))
            self.assertEqual(["G_BLAST"], bridge["gene_id"].tolist())
            row = bridge.iloc[0]
            self.assertEqual("FAM_FALSE", row["family_id"])
            self.assertEqual("blastp", row["match_type"])
            self.assertEqual("tier3_sequence_match_unique", row["mapping_tier"])
            self.assertEqual("none", row["rank_contribution"])
            self.assertNotEqual("", row["input_fingerprint"])
            validate(bridge, FAMILY_CANDIDATE_BRIDGE_SCHEMA)

            manifest = json.loads((out / "family_candidate_bridge_manifest.json").read_text(encoding="utf-8"))
            projection = manifest["route_a_blast"]["falsified_family_projections"]
            self.assertEqual(1, manifest["route_a_blast"]["falsified_family_projection_count"])
            self.assertEqual("falsified_family_not_detected_genome_wide", projection[0]["family_status"])
            self.assertEqual("G_BLAST", projection[0]["blast_top_hit_target_id"])
            self.assertNotEqual("", projection[0]["family_status_input_fingerprint"])
            self.assertNotEqual("", projection[0]["blast_input_fingerprint"])
            self.assertTrue(projection[0]["interval_gene_universe_supplied"])
            self.assertTrue(projection[0]["target_in_interval_gene_universe"])
            self.assertNotEqual("", projection[0]["bridge_output_fingerprint"])

            union_out = root / "union"
            self.assertEqual(0, main([
                "supported-member-union", "--out-dir", str(union_out), "--case-id", "case530",
                "--route-a-bridge", str(out / "family_candidate_bridge.tsv"),
            ]))
            union = pd.read_csv(union_out / "supported_member_union.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(["G_BLAST"], union["candidate_gene_id"].tolist())
            self.assertEqual(["a"], union["route_support"].tolist())
            self.assertEqual(["supported"], union["support_status"].tolist())

            family_out = root / "canonical_family"
            family_result = build_canonical_family_recovery(
                str(union_out / "supported_member_union.tsv"), str(family_out),
                legacy_bridge_arg=str(out / "family_candidate_bridge.tsv"), case_id="case530",
            )
            self.assertEqual(1, family_result["family_count"])
            family_copy_out = root / "canonical_family_copy"
            family_copy_result = family_to_copy(
                family_result["family_recovery"], str(family_copy_out),
                atlas_arg=str(root / "genes.tsv"),
                bridge_arg=family_result["family_bridge"], case_id="case530",
            )
            wet = wet_candidate_handoff(
                str(root / "wet"), route_a_bridge_arg=str(out / "family_candidate_bridge.tsv"),
                family_recovery_arg=family_result["family_recovery"],
                family_copy_arg=family_copy_result["family_copy_discrimination"],
                family_copy_shortlist_arg=family_copy_result["family_copy_shortlist"],
                wet_budget=20, case_id="case530",
            )
            self.assertEqual("wet_reviewable", wet["handoff_state"])
            wet_candidates = pd.read_csv(wet["wet_candidates"], sep="\t", dtype=str).fillna("")
            self.assertEqual(["G_BLAST"], wet_candidates["candidate_gene_id"].tolist())

    def test_falsified_family_still_abstains_for_noneligible_blast_rows(self) -> None:
        cases = {
            "ambiguous": BLAST_EVIDENCE.replace("\tunique\tverified", "\tambiguous_paralog\tverified"),
            "below_threshold": BLAST_EVIDENCE.replace("\t84.924\t", "\t20.0\t"),
            "not_evaluated": BLAST_EVIDENCE.replace("\tevaluable\t", "\tnot_evaluated_provider_failure\t"),
            "unverified": BLAST_EVIDENCE.replace("\tunique\tverified", "\tunique\treview_required"),
            "outside_interval_universe": BLAST_EVIDENCE.replace("G_BLAST", "G_OUTSIDE"),
        }
        for label, blast in cases.items():
            with self.subTest(label=label), tempfile.TemporaryDirectory() as tmp:
                _, bridge = self._run_bridge(Path(tmp), blast)
                self.assertTrue(bridge.empty)

    def test_falsified_family_abstains_without_interval_gene_universe(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            _, bridge = self._run_bridge(Path(tmp), include_genes=False)
            self.assertTrue(bridge.empty)

    def test_direct_projection_does_not_disable_existing_competitor_guard(self) -> None:
        """An independently observed same-family relative still holds BLASTP."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bridge_dir, bridge = self._run_bridge(root)
            self.assertEqual(["G_BLAST"], bridge["gene_id"].tolist())
            # Simulate a separately verified interval competitor in the same
            # family. The Route A projection remains one row; the existing
            # family-to-copy guard must hold its unique BLAST layer because a
            # same-family interval relative is present outside the family's
            # own member set.
            competitor = bridge.iloc[0].copy()
            competitor["gene_id"] = "G_COMPETITOR"
            competitor["route"] = "a"
            competitor["mapping_tier"] = "tier3_sequence_match_unique"
            competitor["output_fingerprint"] = "competitor"
            competitor_path = root / "bridge_with_competitor.tsv"
            pd.concat([bridge, pd.DataFrame([competitor])], ignore_index=True)[BRIDGE_COLUMNS].to_csv(
                competitor_path, sep="\t", index=False, lineterminator="\n"
            )
            families = _write(root, "families.tsv", "hypothesis_id\tgene_family\tmember_interval_genes\nH_FALSE\tFAM_FALSE\tG_BLAST\n")
            atlas = _write(root, "atlas.tsv", "gene_id\tchrom\tstart\tend\nG_BLAST\tchr1\t1000\t3000\n")
            copy_out = root / "family_copy"
            family_to_copy(
                str(families), str(copy_out), atlas_arg=str(atlas), bridge_arg=str(competitor_path),
                route_a_blast_arg=str(root / "blast.tsv"), case_id="case530",
            )
            discrimination = pd.read_csv(copy_out / "family_copy_discrimination.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("interval_paralog_present", discrimination.iloc[0]["route_a_blast_guard"])
            self.assertEqual("", discrimination.iloc[0]["copy_set_evidence_source"])
            self.assertFalse((copy_out / "family_copy_shortlist.tsv").read_text(encoding="utf-8").strip().splitlines()[1:])


if __name__ == "__main__":
    unittest.main()
