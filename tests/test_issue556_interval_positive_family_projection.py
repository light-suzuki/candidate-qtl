"""Issue #556: Route A positive projection through the frozen B0 interval index."""

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


class IntervalPositiveFamilyProjectionTests(unittest.TestCase):
    def _inputs(self, root: Path, *, with_bait: bool = True) -> tuple[Path, Path, Path, Path]:
        hypotheses = _write(root, "hypotheses.tsv", """hypothesis_id\tgene_family\treference_gene_name\tproposed_stable_id\tmember_interval_genes
H_MYB\tR2R3-MYB26/MS35 anther developmental regulators\tMYB26\tQ9FJA2\t
""")
        status = _write(root, "status.tsv", """hypothesis_id\tfamily_status\tinterval_gene_ids\tsource_ids
H_MYB\tnot_evaluated_missing_reference_genome_family_index\t\tPMID:123
""")
        bait_text = "hypothesis_id\tfunction_source_ids\tsource_provider\tbait_id\n"
        if with_bait:
            bait_text += "H_MYB\tPMID:123\tsource_verified\tQ9FJA2\n"
        bait = _write(root, "bait.tsv", bait_text)
        index = _write(root, "family_index.tsv", """gene_id\tassignment_status\tindex_scope\tbest_subject\tbest_subject_symbol\tbest_subject_stem\tfamily_symbols\tfamily_stems\tdiscovery_bucket_id\tevolutionary_family_id\tbest_identity\tbest_query_coverage\tbest_evalue
Ps1\tassigned\tinterval_candidates\tsp|Q9LK95|MYB21_ARATH\tMYB21\tMYB\tMYB21;MYB26\tMYB\tMYB\t\t63.927\t0.8247\t7.03e-87
Ps2\tassigned\tinterval_candidates\ttr|A0A|MYB_LOTJA\tMYB14\tMYB\tMYB14\tMYB\tMYB\t\t58.2\t0.71\t2e-50
Ps3\tassigned\tinterval_candidates\tsp|Q9|PAL_ARATH\tPAL1\tPAL\tPAL1\tPAL\tPAL\t\t72.0\t0.91\t1e-100
""")
        return hypotheses, status, bait, index

    def test_verified_family_retains_every_matching_interval_copy(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hypotheses, status, bait, index = self._inputs(root)
            out = root / "out"
            self.assertEqual(0, main([
                "family-candidate-bridge", "--out-dir", str(out), "--case-id", "case556",
                "--route-a-hypotheses", str(hypotheses),
                "--route-a-family-status", str(status),
                "--route-a-bait-panel", str(bait),
                "--route-b-family-index", str(index),
            ]))
            bridge = pd.read_csv(out / "family_candidate_bridge.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(["Ps1", "Ps2"], bridge["gene_id"].tolist())
            self.assertEqual({"a"}, set(bridge["route"]))
            self.assertEqual({"paralog_ambiguous"}, set(bridge["ambiguity_status"]))
            self.assertEqual({"tier4_sequence_match_ambiguous"}, set(bridge["mapping_tier"]))
            self.assertEqual({"none"}, set(bridge["rank_contribution"]))
            manifest = json.loads((out / "family_candidate_bridge_manifest.json").read_text(encoding="utf-8"))
            projection = manifest["route_a_interval_family_projection"]
            self.assertEqual(2, projection["projection_count"])
            self.assertTrue(projection["input_present"])
            self.assertTrue(manifest["inputs"]["route_b_family_index"]["present"])
            self.assertTrue(all(row["positive_only"] for row in projection["projections"]))
            self.assertTrue(all(not row["absence_inference_allowed"] for row in projection["projections"]))

            union = root / "union"
            self.assertEqual(0, main([
                "supported-member-union", "--out-dir", str(union), "--case-id", "case556",
                "--route-a-bridge", str(out / "family_candidate_bridge.tsv"),
            ]))
            supported = pd.read_csv(
                union / "supported_member_union.tsv", sep="\t", dtype=str
            ).fillna("")
            self.assertEqual(["Ps1", "Ps2"], supported["candidate_gene_id"].tolist())
            self.assertEqual({"none"}, set(supported["rank_contribution"]))

    def test_unverified_or_nonmatching_family_is_not_promoted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hypotheses, status, bait, index = self._inputs(root, with_bait=False)
            out = root / "unverified"
            self.assertEqual(0, main([
                "family-candidate-bridge", "--out-dir", str(out), "--case-id", "case556",
                "--route-a-hypotheses", str(hypotheses),
                "--route-a-family-status", str(status),
                "--route-a-bait-panel", str(bait), "--route-b-family-index", str(index),
            ]))
            bridge = pd.read_csv(out / "family_candidate_bridge.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(0, len(bridge))

            # A verified bait whose family has no interval hit is equally
            # non-negative: no promoted row, no fabricated miss.
            hypotheses, status, bait, index = self._inputs(root, with_bait=True)
            hypotheses.write_text(
                hypotheses.read_text(encoding="utf-8")
                .replace("R2R3-MYB26/MS35 anther developmental regulators", "NAC transcription factors")
                .replace("MYB26", "NAC42")
                .replace("Q9FJA2", "Q9NAC42"),
                encoding="utf-8",
            )
            verified_out = root / "verified_no_hit"
            self.assertEqual(0, main([
                "family-candidate-bridge", "--out-dir", str(verified_out), "--case-id", "case556",
                "--route-a-hypotheses", str(hypotheses),
                "--route-a-family-status", str(status),
                "--route-a-bait-panel", str(bait), "--route-b-family-index", str(index),
            ]))
            verified_bridge = pd.read_csv(
                verified_out / "family_candidate_bridge.tsv", sep="\t", dtype=str
            ).fillna("")
            self.assertEqual(0, len(verified_bridge))
            verified_manifest = json.loads(
                (verified_out / "family_candidate_bridge_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(0, verified_manifest["route_a_interval_family_projection"]["projection_count"])

    def test_reference_genome_scope_is_not_reused_as_interval_positive_seam(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hypotheses, status, bait, index = self._inputs(root)
            text = index.read_text(encoding="utf-8").replace("interval_candidates", "reference_genome")
            index.write_text(text, encoding="utf-8")
            out = root / "wrong_scope"
            self.assertEqual(0, main([
                "family-candidate-bridge", "--out-dir", str(out), "--case-id", "case556",
                "--route-a-hypotheses", str(hypotheses),
                "--route-a-family-status", str(status),
                "--route-a-bait-panel", str(bait), "--route-b-family-index", str(index),
            ]))
            bridge = pd.read_csv(out / "family_candidate_bridge.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(0, len(bridge))

    def test_verified_pi_symbol_in_broad_mads_label_matches_pistillata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hypotheses = _write(root, "hypotheses.tsv", """hypothesis_id\tgene_family\treference_gene_name\treference_genes\tproposed_stable_id\tmember_interval_genes
H_PI\tMADS-box floral homeotic factors (AP1/AP3/PI/AG)\t\tAP3\tP35632\t
""")
            status = _write(root, "status.tsv", """hypothesis_id\tfamily_status\tinterval_gene_ids\tsource_ids
H_PI\tnot_evaluated_missing_reference_genome_family_index\t\tPMID:9243505
""")
            bait = _write(root, "bait.tsv", """hypothesis_id\tfunction_source_ids\tsource_provider\tbait_id
H_PI\tPMID:9243505\tsource_verified\tP35632
""")
            index = _write(root, "family_index.tsv", """gene_id\tassignment_status\tindex_scope\tbest_subject\tbest_subject_symbol\tbest_subject_stem\tfamily_symbols\tfamily_stems\tdiscovery_bucket_id\tevolutionary_family_id\tbest_identity\tbest_query_coverage\tbest_evalue
PsPI\tassigned\tinterval_candidates\ttr|B7TX67|B7TX67_MEDTR\tPISTILLATA\tPISTILLATA\tPISTILLATA\tPISTILLATA\tPISTILLATA\t\t96.685\t1.0\t2.9e-129
PsPAL\tassigned\tinterval_candidates\tsp|P45724|PAL1_PEA\tPAL1\tPAL\tPAL1\tPAL\tPAL\t\t75.0\t0.9\t1e-100
""")
            out = root / "out"
            self.assertEqual(0, main([
                "family-candidate-bridge", "--out-dir", str(out), "--case-id", "case562",
                "--route-a-hypotheses", str(hypotheses),
                "--route-a-family-status", str(status),
                "--route-a-bait-panel", str(bait),
                "--route-b-family-index", str(index),
            ]))
            bridge = pd.read_csv(out / "family_candidate_bridge.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(["PsPI"], bridge["gene_id"].tolist())
            self.assertEqual({"a"}, set(bridge["route"]))
            manifest = json.loads((out / "family_candidate_bridge_manifest.json").read_text(encoding="utf-8"))
            projection = manifest["route_a_interval_family_projection"]
            self.assertEqual(1, projection["projection_count"])
            self.assertEqual("PISTILLATA", projection["projections"][0]["interval_family_label"])
            self.assertTrue(projection["projections"][0]["positive_only"])
            self.assertFalse(projection["projections"][0]["absence_inference_allowed"])

    def test_ambiguous_pi_kinase_does_not_match_pistillata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hypotheses = _write(root, "hypotheses.tsv", """hypothesis_id\tgene_family\treference_gene_name\tproposed_stable_id\tmember_interval_genes
H_PI_KINASE\tPI kinase family\tPI4K\tQ9SWR3\t
""")
            status = _write(root, "status.tsv", """hypothesis_id\tfamily_status\tinterval_gene_ids\tsource_ids
H_PI_KINASE\tnot_evaluated_missing_reference_genome_family_index\t\tPMID:1
""")
            bait = _write(root, "bait.tsv", """hypothesis_id\tfunction_source_ids\tsource_provider\tbait_id
H_PI_KINASE\tPMID:1\tsource_verified\tQ9SWR3
""")
            index = _write(root, "family_index.tsv", """gene_id\tassignment_status\tindex_scope\tbest_subject\tbest_subject_symbol\tbest_subject_stem\tfamily_symbols\tfamily_stems\tdiscovery_bucket_id\tevolutionary_family_id\tbest_identity\tbest_query_coverage\tbest_evalue
PsPI\tassigned\tinterval_candidates\ttr|B7TX67|B7TX67_MEDTR\tPISTILLATA\tPISTILLATA\tPISTILLATA\tPISTILLATA\tPISTILLATA\t\t96.685\t1.0\t2.9e-129
""")
            out = root / "out"
            self.assertEqual(0, main([
                "family-candidate-bridge", "--out-dir", str(out), "--case-id", "case562neg",
                "--route-a-hypotheses", str(hypotheses),
                "--route-a-family-status", str(status),
                "--route-a-bait-panel", str(bait),
                "--route-b-family-index", str(index),
            ]))
            bridge = pd.read_csv(out / "family_candidate_bridge.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(0, len(bridge))

            hypotheses.write_text(
                hypotheses.read_text(encoding="utf-8").replace(
                    "PI kinase family", "MADS-box and PI kinase families"
                ),
                encoding="utf-8",
            )
            mixed_out = root / "mixed_out"
            self.assertEqual(0, main([
                "family-candidate-bridge", "--out-dir", str(mixed_out),
                "--case-id", "case562mixedneg",
                "--route-a-hypotheses", str(hypotheses),
                "--route-a-family-status", str(status),
                "--route-a-bait-panel", str(bait),
                "--route-b-family-index", str(index),
            ]))
            mixed_bridge = pd.read_csv(
                mixed_out / "family_candidate_bridge.tsv", sep="\t", dtype=str
            ).fillna("")
            self.assertEqual(0, len(mixed_bridge))

            for index_value, label in enumerate((
                "MADS-box factors (non-PI/AG)",
                "MADS-box factors (without PI/AG)",
                "not PI-class MADS factor",
                "MADS factor not PI-class",
                "not PISTILLATA",
                "excluding PISTILLATA",
            )):
                hypotheses.write_text(
                    "hypothesis_id\tgene_family\treference_gene_name\tproposed_stable_id\tmember_interval_genes\n"
                    f"H_PI_KINASE\t{label}\tPI4K\tQ9SWR3\t\n",
                    encoding="utf-8",
                )
                negated_out = root / f"negated_out_{index_value}"
                self.assertEqual(0, main([
                    "family-candidate-bridge", "--out-dir", str(negated_out),
                    "--case-id", f"case562negated{index_value}",
                    "--route-a-hypotheses", str(hypotheses),
                    "--route-a-family-status", str(status),
                    "--route-a-bait-panel", str(bait),
                    "--route-b-family-index", str(index),
                ]))
                negated_bridge = pd.read_csv(
                    negated_out / "family_candidate_bridge.tsv", sep="\t", dtype=str
                ).fillna("")
                self.assertEqual(0, len(negated_bridge), label)

    def test_guided_profile_wires_route_b0_index_after_freeze(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "plan"
            self.assertEqual(0, main([
                "guided-profile", "--profile", "triple-route", "--out-dir", str(out),
            ]))
            manifest = json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            bridge = next(stage for stage in manifest["stages"] if stage["command"] == "family-candidate-bridge")
            self.assertEqual(
                {"stage": "route-b-prepare", "result_key": "family_index"},
                bridge["inputs_from"]["route_b_family_index"],
            )


if __name__ == "__main__":
    unittest.main()
