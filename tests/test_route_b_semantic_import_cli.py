from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main


class RouteBSemanticImportCliTests(unittest.TestCase):
    def test_cli_imports_only_non_scoring_proposals(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atlas = root / "atlas.tsv"
            pd.DataFrame([{
                "gene_id": "G1", "evolutionary_family_id": "EV1", "orthogroup": "OG1",
                "family_id": "F1", "orthology_reference_gene": "REF1",
            }]).to_csv(atlas, sep="\t", index=False, lineterminator="\n")
            literature = root / "literature.tsv"
            pd.DataFrame([{
                "interval_gene_ids": "G1", "reference_gene": "REF1", "reference_function": "meristem regulator",
                "query_status": "retrieved", "article_rank": "1", "source_id": "PMID:123",
                "title": "Reference meristem paper", "doi": "", "journal": "Plant", "year": "2020",
                "matched_phenotype_terms": "fasciation", "claim_type": "literature_link_candidate",
            }]).to_csv(literature, sep="\t", index=False, lineterminator="\n")

            self.assertEqual(0, main([
                "route-b-semantic-bundle",
                "--atlas", str(atlas),
                "--literature", str(literature),
                "--phenotype", "fasciated meristem",
                "--out-dir", str(root / "bundle"),
                "--case-id", "cli-import",
            ]))
            prompts = root / "bundle" / "route_b_semantic_prompts.jsonl"
            prompt = json.loads(prompts.read_text(encoding="utf-8").splitlines()[0])
            claims = [{
                "claim_type": "mutant_phenotype",
                "biological_function": "meristem maintenance",
                "mechanism": "receptor signaling",
                "phenotype_relation": "may explain fasciation",
                "reference_gene": "REF1",
                "reference_species": "Arabidopsis thaliana",
                "evidence_scope": "reference_gene",
                "source_ids": ["PMID:123"],
                "alternative_hypothesis": "",
                "uncertainty_note": "relation still requires verification",
            }]
            responses = root / "responses.json"
            responses.write_text(json.dumps({"candidates": [{
                "prompt_id": prompt["prompt_id"], "gene_id": prompt["gene_id"],
                "claims": claims, "model_json": {"claims": claims},
                "raw_text": json.dumps({"claims": claims}), "parse_status": "json",
                "review_status": "not_reviewed",
            }]}, ensure_ascii=False) + "\n", encoding="utf-8")

            self.assertEqual(0, main([
                "route-b-semantic-import",
                "--bundle", str(root / "bundle" / "route_b_semantic_bundle.json"),
                "--prompts", str(prompts),
                "--responses", str(responses),
                "--out-dir", str(root / "import"),
                "--case-id", "cli-import",
            ]))
            proposals = pd.read_csv(
                root / "import" / "semantic_claim_proposals.tsv", sep="\t", dtype=str
            ).fillna("")
            self.assertEqual(["proposed"], proposals["verification_status"].tolist())
            self.assertEqual(["none"], proposals["rank_contribution"].tolist())
            manifest = json.loads(
                (root / "import" / "semantic_claim_proposal_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(0, manifest["verified_claim_row_count"])


if __name__ == "__main__":
    unittest.main()
