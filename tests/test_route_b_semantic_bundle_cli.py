from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.route_b_semantic_import import import_route_b_semantic_proposals


class RouteBSemanticBundleCliTests(unittest.TestCase):
    def test_familyless_reference_prompt_reaches_prepare_with_schema_valid_example(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atlas = root / "atlas.tsv"
            pd.DataFrame([{
                "gene_id": "G1", "family_id": "", "evolutionary_family_id": "", "orthogroup": "",
                "orthology_reference_gene": "REF1", "domain_architecture": "kinase",
            }]).to_csv(atlas, sep="\t", index=False, lineterminator="\n")
            original_atlas = atlas.read_bytes()
            out = root / "bundle"
            env = {key: value for key, value in os.environ.items() if key != "PYTHONPATH"}
            subprocess.run([
                sys.executable, "-m", "candidate_qtl", "route-b-semantic-bundle",
                "--atlas", str(atlas), "--phenotype", "altered organ development",
                "--out-dir", str(out), "--case-id", "scope-contract",
            ], check=True, capture_output=True, text=True, env=env)
            prompts = out / "route_b_semantic_prompts.jsonl"
            prompt = json.loads(prompts.read_text(encoding="utf-8").strip())
            system = prompt["system"]
            self.assertIn("family requires a non-empty family_id", system)
            self.assertIn("do not emit family-scoped claims or invent a family_id", system)
            example = json.loads(system.split("Return JSON only as ", 1)[1].splitlines()[0].removesuffix("."))
            claim = example["claims"][0]
            self.assertEqual("reference_gene", claim["evidence_scope"])
            self.assertTrue(claim["reference_gene"])
            # The example becomes importable with supplied endpoints and actual
            # claim fields; changing only its scope to family stays rejected.
            claim.update({"claim_type": "mechanistic_review_seed", "reference_gene": "REF1", "source_ids": [],
                          "paper_reports_json": []})
            bundle = json.loads((out / "route_b_semantic_bundle.json").read_text(encoding="utf-8"))
            self.assertEqual(["G1"], bundle["units"][0]["member_interval_genes"])
            self.assertEqual([], bundle["units"][0]["identity"]["family_ids"])
            self.assertEqual([], bundle["units"][0]["identity"]["evolutionary_family_ids"])
            for scope, imported, rejected in (("reference_gene", 1, 0), ("family", 0, 1)):
                scoped = {**claim, "evidence_scope": scope}
                response = root / f"{scope}.json"
                response.write_text(json.dumps({"candidates": [{
                    "prompt_id": prompt["prompt_id"], "gene_id": prompt["gene_id"],
                    "model_json": {"claims": [scoped]}, "raw_text": json.dumps({"claims": [scoped]}),
                    "claims": [scoped], "parse_status": "json", "review_status": "not_reviewed",
                }]}), encoding="utf-8")
                result = import_route_b_semantic_proposals(
                    str(out / "route_b_semantic_bundle.json"), str(prompts), str(response), str(root / scope),
                )
                self.assertEqual(imported, result["claim_row_count"])
                self.assertEqual(rejected, result["rejected_claim_row_count"])
            script = Path(__file__).resolve().parents[1] / "scripts" / "agentic_route_b.py"
            prepared = root / "prepared"
            subprocess.run([
                sys.executable, str(script), "prepare", "--prompts", str(prompts), "--out", str(prepared),
            ], check=True, capture_output=True, text=True, env=env)
            exported = json.loads((prepared / "b0000" / "prompts.jsonl").read_text(encoding="utf-8").strip())
            self.assertEqual(prompt, exported)
            self.assertEqual(original_atlas, atlas.read_bytes())

    def test_cli_builds_post_b0_bundle_without_changing_candidate_identity(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atlas = root / "atlas.tsv"
            pd.DataFrame([
                {"gene_id": "G1", "evolutionary_family_id": "EV1", "orthogroup": "OG1", "family_id": "F1"},
                {"gene_id": "G2", "evolutionary_family_id": "EV1", "orthogroup": "OG1", "family_id": "F1"},
                {"gene_id": "DARK1", "evolutionary_family_id": "", "orthogroup": "", "family_id": ""},
            ]).to_csv(atlas, sep="\t", index=False, lineterminator="\n")

            self.assertEqual(
                0,
                main([
                    "route-b-semantic-bundle",
                    "--atlas", str(atlas),
                    "--phenotype", "fasciated meristem",
                    "--out-dir", str(root / "out"),
                    "--case-id", "cli403bundle",
                ]),
            )
            manifest = json.loads(
                (root / "out" / "route_b_semantic_bundle_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual("cli403bundle", manifest["case_id"])
            self.assertEqual(3, manifest["interval_gene_count"])
            self.assertEqual(2, manifest["semantic_unit_count"])
            self.assertTrue(manifest["candidate_identity_preserved"])
            self.assertEqual(1, manifest["gene_level_unclassified_unit_count"])
            self.assertEqual("none", manifest["rank_contribution"])


if __name__ == "__main__":
    unittest.main()
