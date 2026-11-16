"""Issue #532: preserve every eligible Route A hit before copy review."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main


def _write(root: Path, name: str, text: str) -> Path:
    path = root / name
    path.write_text(text, encoding="utf-8")
    return path


class Issue532RouteACompleteHitSetTests(unittest.TestCase):
    def test_falsified_family_keeps_all_long_form_threshold_passing_targets(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hypotheses = _write(root, "hypotheses.tsv", """hypothesis_id\tgene_family\treference_gene_name\tproposed_stable_id\tmember_interval_genes
H_FALSE\tFAM_FALSE\tREF_FALSE\tACC_FALSE\tG_TOP;G_TRUTH
""")
            bait = _write(root, "bait.tsv", """hypothesis_id\tfunction_source_ids\tsource_provider
H_FALSE\tROOT_FALSE\tlocal_verified_panel
""")
            status = _write(root, "status.tsv", """hypothesis_id\tfamily_status\tinterval_gene_ids\tsource_ids\treference_gene_name\tproposed_stable_id
H_FALSE\tfalsified_family_not_detected_genome_wide\t\tROOT_FALSE\tREF_FALSE\tACC_FALSE
""")
            blast = _write(root, "blast_evidence.tsv", """panel_id\tprovider_status\ttop_hit_target_id\ttop_hit_identity\ttop_hit_query_coverage\ttop_hit_target_coverage\ttop_hit_evalue\thit_count\tparalog_count\tambiguity_status\treview_status\trank_contribution
ACC_FALSE\tevaluable\tG_TOP\t84.924\t0.982625\t0.986792\t0.0\t2\t1\tambiguous_paralog\tverified\tnone
""")
            eligible = _write(root, "eligible_interval_hits.tsv", """panel_id\ttarget_id\treview_status\tprovider_status\tambiguity_status\tpident\tquery_coverage\ttarget_coverage\tevalue\tbitscore\teligibility\trank_contribution
ACC_FALSE\tG_TOP\tverified\tevaluable\tambiguous_paralog\t84.924\t0.982625\t0.986792\t0.0\t441\teligible\tnone
ACC_FALSE\tG_TRUTH\tverified\tevaluable\tambiguous_paralog\t81.791\t0.962625\t0.976792\t0.0\t392\teligible\tnone
""")
            genes = _write(root, "genes.tsv", """gene_id\tchrom\tstart\tend
G_TOP\tchr1\t1000\t3000
G_TRUTH\tchr1\t5000\t7000
""")
            out = root / "bridge"
            self.assertEqual(0, main([
                "family-candidate-bridge", "--out-dir", str(out), "--case-id", "case532",
                "--route-a-hypotheses", str(hypotheses), "--route-a-family-status", str(status),
                "--route-a-bait-panel", str(bait), "--route-a-blast", str(blast),
                "--route-a-blast-hits", str(eligible),
                "--genes", str(genes),
            ]))
            bridge = pd.read_csv(out / "family_candidate_bridge.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(["G_TOP", "G_TRUTH"], bridge["gene_id"].tolist())
            self.assertEqual({"tier4_sequence_match_ambiguous"}, set(bridge["mapping_tier"]))
            self.assertEqual({"paralog_ambiguous"}, set(bridge["ambiguity_status"]))
            self.assertEqual({"none"}, set(bridge["rank_contribution"]))
            manifest = json.loads((out / "family_candidate_bridge_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(manifest["route_a_blast"]["long_form_input_present"])
            self.assertEqual(2, manifest["route_a_blast"]["long_form_evidence_rows"])
            self.assertEqual(2, manifest["route_a_blast"]["falsified_family_projection_count"])
            self.assertTrue(manifest["inputs"]["route_a_blast_hits"]["present"])

            # Even when the eligible artifact contains only one retained
            # target, a compact ambiguous/non-unique verdict must not be
            # upgraded to resolved_single_copy.
            eligible.write_text(
                "panel_id\ttarget_id\treview_status\tprovider_status\tambiguity_status\t"
                "pident\tquery_coverage\ttarget_coverage\tevalue\tbitscore\teligibility\trank_contribution\n"
                "ACC_FALSE\tG_TOP\tverified\tevaluable\tambiguous_paralog\t84.924\t"
                "0.982625\t0.986792\t0.0\t441\teligible\tnone\n",
                encoding="utf-8",
            )
            single_out = root / "bridge_single_eligible"
            self.assertEqual(0, main([
                "family-candidate-bridge", "--out-dir", str(single_out), "--case-id", "case532",
                "--route-a-hypotheses", str(hypotheses), "--route-a-family-status", str(status),
                "--route-a-bait-panel", str(bait), "--route-a-blast", str(blast),
                "--route-a-blast-hits", str(eligible), "--genes", str(genes),
            ]))
            single_bridge = pd.read_csv(
                single_out / "family_candidate_bridge.tsv", sep="\t", dtype=str
            ).fillna("")
            self.assertEqual(["G_TOP"], single_bridge["gene_id"].tolist())
            self.assertEqual({"paralog_ambiguous"}, set(single_bridge["ambiguity_status"]))
            self.assertEqual({"tier4_sequence_match_ambiguous"}, set(single_bridge["mapping_tier"]))

            # A present long-form artifact is authoritative. If it carries no
            # admissible rows, the bridge must not fall back to the compact
            # top hit and recreate the #532 collapse.
            eligible.write_text(
                "panel_id\ttarget_id\treview_status\tprovider_status\tambiguity_status\t"
                "pident\tquery_coverage\ttarget_coverage\tevalue\tbitscore\teligibility\trank_contribution\n",
                encoding="utf-8",
            )
            empty_out = root / "bridge_empty"
            self.assertEqual(0, main([
                "family-candidate-bridge", "--out-dir", str(empty_out), "--case-id", "case532",
                "--route-a-hypotheses", str(hypotheses), "--route-a-family-status", str(status),
                "--route-a-bait-panel", str(bait), "--route-a-blast", str(blast),
                "--route-a-blast-hits", str(eligible), "--genes", str(genes),
            ]))
            empty_bridge = pd.read_csv(
                empty_out / "family_candidate_bridge.tsv", sep="\t", dtype=str
            ).fillna("")
            self.assertEqual(0, len(empty_bridge))


if __name__ == "__main__":
    unittest.main()
