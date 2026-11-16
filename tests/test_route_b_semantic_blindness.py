from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.route_b_semantic_bundle import build_route_b_semantic_bundle
from candidate_qtl.route_b_semantic_candidates import freeze_route_b_semantic_candidates
from candidate_qtl.semantic_claims import CLAIM_COLUMNS


class RouteBSemanticBlindnessTests(unittest.TestCase):
    def test_route_b_bundle_is_route_a_blind(self) -> None:
        """Route B bundle consumes only B0 atlas and optional reverse literature, never Route A hypotheses."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atlas = root / "atlas.tsv"
            pd.DataFrame([
                {"gene_id": "G1", "evolutionary_family_id": "EV1", "orthogroup": "OG1", "family_id": "F1"},
            ]).to_csv(atlas, sep="\t", index=False, lineterminator="\n")

            result = build_route_b_semantic_bundle(
                str(atlas), str(root / "bundle"), phenotype="fasciated stem", case_id="blind_test"
            )
            bundle = json.loads(Path(result["bundle"]).read_text(encoding="utf-8"))
            self.assertEqual("interval_first", bundle["route"])
            self.assertIn("Route A artifacts and benchmark truth are not inputs", bundle["policy"])
            prompts = [json.loads(line) for line in Path(result["prompts"]).read_text(encoding="utf-8").splitlines()]
            for prompt in prompts:
                self.assertNotIn("hypothesis_id", prompt["user"])
                self.assertNotIn("route_a", prompt["user"].lower())

    def test_route_b_candidate_freeze_is_truth_blind(self) -> None:
        """Route B candidate freeze accepts only semantic claims with interval_first route, never benchmark truth."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            claim = {column: "" for column in CLAIM_COLUMNS}
            claim.update({
                "claim_id": "c1", "route": "interval_first", "semantic_unit_id": "bsem-blind",
                "family_id": "EV1", "member_interval_genes": "G1", "claim_type": "mutant_phenotype",
                "biological_function": "meristem maintenance", "mechanism": "receptor signaling",
                "phenotype_relation": "fasciated shoot meristem", "reference_gene": "REF1",
                "reference_species": "Arabidopsis thaliana", "evidence_scope": "family",
                "source_id": "PMID:1234567", "provider": "fixture", "model": "fixture",
                "prompt_fingerprint": "p" * 64, "input_fingerprint": "i" * 64,
                "output_fingerprint": "o" * 64, "verification_status": "verified",
                "verification_reason": "group verified", "rank_contribution": "none",
            })
            claims_file = root / "claims.tsv"
            pd.DataFrame([claim], columns=CLAIM_COLUMNS).to_csv(claims_file, sep="\t", index=False, lineterminator="\n")

            result = freeze_route_b_semantic_candidates(str(claims_file), str(root / "out"), case_id="blind_test")
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual("interval_first", manifest["route"])
            self.assertIn("this artifact contains Route B only and is intended to freeze before A/B comparison", manifest["policy"])

            claim_wrong = dict(claim)
            claim_wrong["route"] = "phenotype_first"
            wrong_file = root / "wrong.tsv"
            pd.DataFrame([claim_wrong], columns=CLAIM_COLUMNS).to_csv(wrong_file, sep="\t", index=False, lineterminator="\n")
            with self.assertRaisesRegex(ValueError, "interval_first"):
                freeze_route_b_semantic_candidates(str(wrong_file), str(root / "out2"))


if __name__ == "__main__":
    unittest.main()
