"""Offline initial-material delivery and accounting; no biological re-rating."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from scripts.agentic_route_b_initial import collect_material, summarize_initial
from scripts.agentic_route_b_research import plan_research_view
from tests.test_agentic_route_b_research import fixture

ROOT = Path(__file__).resolve().parents[1]


class InitialMaterialTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.bundle, self.prompts = fixture()
        self.source = self.root / "bundle.json"

    def save(self):
        self.source.write_text(json.dumps(self.bundle), encoding="utf-8")
        self.sha = hashlib.sha256(self.source.read_bytes()).hexdigest()
        return self.sha

    def answers(self, rows):
        path = self.root / "answers.json"
        path.write_text(json.dumps({"bundle_sha256": self.sha, "assessments": rows}), encoding="utf-8")
        return path

    @staticmethod
    def judgment(gid="g2", state="assessed_no_specific_bridge", pointer="/units/0/genes/1", fields=None):
        return {"unit_id": "u1", "gene_id": gid, "state": state,
                "reason": "Original judgment, not a new scientific interpretation.",
                "used_material_pointers": [pointer] if pointer else [],
                "assessed_fields": fields if fields is not None else ["reference_function"]}

    def test_same_function_different_context_is_not_collapsed(self):
        cards = self.bundle["units"][0]["reference_material_cards"]
        base = {"interval_gene_id": "g2", "function_description": "same kinase name"}
        cards[:] = [{**base, "reference_gene": "r1", "keywords": "developmental context " * 30},
                    {**base, "reference_gene": "r2", "keywords": "stress context", "tissue_specificity": "root"},
                    {**base, "reference_gene": "r3", "keywords": "stress context", "tissue_specificity": "root"}]
        files, _ = plan_research_view(str(self.source), self.save(), self.prompts, "fixture", 400)
        manifest = json.loads(files["research/manifest.json"])
        material = manifest["initial_material"]
        groups = [g for p in material["pages"] for g in json.loads(files[p["file"]])]
        kinase = [g for g in groups if g["context"].get("function_description") == "same kinase name"]
        self.assertEqual(2, len(kinase))
        self.assertEqual([1, 2], sorted(len(g["bindings"]) for g in kinase))
        self.assertIn("developmental context " * 30, [g["context"]["keywords"] for g in kinase])
        self.assertTrue(any(p["over_bound_bytes"] for p in material["pages"]))
        self.assertEqual(["g3", "g4"], material["genes_without_material"])
        self.assertEqual("not_assessed", material["interpretation_state"])

    def test_mixed_judgments_preserve_interpreted_and_unknown_parts(self):
        self.save()
        rows = [self.judgment(), self.judgment(state="assessed_insufficient_material",
                                              pointer="/units/0/reference_material_cards/0", fields=["function_description"])]
        result = summarize_initial(str(self.source), self.sha, str(self.answers(rows)), str(self.root / "out"))
        genes = [json.loads(line) for line in (self.root / "out/initial_assessment_genes.jsonl").read_text().splitlines()]
        self.assertEqual(4, result["gene_count"])
        self.assertEqual({"not_assessed": 3, "partially_interpreted": 1}, result["interpretation_counts"])
        gene = next(g for g in genes if g["gene_id"] == "g2")
        self.assertTrue(gene["has_no_bridge_judgment"])
        self.assertTrue(gene["has_unresolved_judgment"])
        self.assertFalse(gene["has_original_claim"])
        self.assertEqual(rows, gene["original_assessments"])
        self.assertEqual(self.sha, hashlib.sha256(self.source.read_bytes()).hexdigest())

    def test_partial_field_receipt_does_not_claim_other_fields_read(self):
        self.bundle["units"][0]["genes"][1]["keywords"] = "context outside the original reading scope"
        self.save()
        out = self.root / "out"
        summarize_initial(str(self.source), self.sha, str(self.answers([self.judgment()])), str(out))
        genes = [json.loads(line) for line in (out / "initial_assessment_genes.jsonl").read_text().splitlines()]
        gene = next(g for g in genes if g["gene_id"] == "g2")
        self.assertEqual(["keywords"], gene["fields_without_read_receipt"]["/units/0/genes/1"])
        self.assertEqual("interpreted", gene["interpretation_state"])
        self.assertIn("/units/0/reference_material_cards/0", gene["fields_without_read_receipt"])

    def test_missing_material_can_remain_insufficient_without_search(self):
        self.save()
        rows = [{**self.judgment("g4", "assessed_insufficient_material", None, []), "unit_id": "u3"}]
        result = summarize_initial(str(self.source), self.sha, str(self.answers(rows)), str(self.root / "out"))
        self.assertEqual({"not_assessed": 3, "insufficient_only": 1}, result["interpretation_counts"])
        self.assertEqual({"present": 2, "absent": 2}, result["material_counts"])

    def test_legacy_hit_only_absent_material_shared_assessment_replays(self):
        self.assertNotIn("initial_material_projection", self.bundle)
        self.assertTrue(self.bundle["units"][1]["genes"][0]["sequence_similarity_hits"][0]["hit_protein_name"])
        self.save()
        before = self.source.read_bytes()
        answers = self.root / "legacy-shared.json"
        answers.write_text(json.dumps({"bundle_sha256": self.sha, "shared_assessments": [{
            "assessment_id": "legacy-no-material", "material_id": "", "member_gene_ids": ["g3"],
            "state": "assessed_insufficient_material", "reason": "Original frozen absent-material accounting.",
            "assessed_fields": [],
        }]}), encoding="utf-8")
        result = summarize_initial(str(self.source), self.sha, str(answers), str(self.root / "legacy-out"))
        self.assertEqual({"present": 2, "absent": 2}, result["material_counts"])
        self.assertEqual({"not_assessed": 3, "insufficient_only": 1}, result["interpretation_counts"])
        self.assertEqual(before, self.source.read_bytes())

    def test_wrong_binding_hash_and_duplicate_fail_before_writes(self):
        self.save()
        invalid = [
            [self.judgment("g1", pointer="/units/0/reference_material_cards/0", fields=["function_description"])],
            [self.judgment(fields=["nonexistent_field"])],
            [self.judgment(), self.judgment()],
        ]
        for rows in invalid:
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                summarize_initial(str(self.source), self.sha, str(self.answers(rows)), str(self.root / "out"))
            self.assertFalse((self.root / "out").exists())
        with self.assertRaisesRegex(ValueError, "SHA-256"):
            summarize_initial(str(self.source), "bad", str(self.answers([])), str(self.root / "out"))

    def test_explicit_card_outside_unit_is_rejected(self):
        self.bundle["units"][0]["reference_material_cards"][0]["interval_gene_id"] = "g3"
        with self.assertRaisesRegex(ValueError, "outside unit"):
            collect_material(self.bundle)

    def test_old_judgment_without_pointer_is_preserved_as_unrecorded(self):
        self.save()
        result = summarize_initial(str(self.source), self.sha,
                                   str(self.answers([self.judgment(pointer=None, fields=[])])), str(self.root / "out"))
        self.assertEqual(1, result["genes_with_judgments_without_material_pointer"])
        self.assertEqual(1, result["interpretation_counts"]["interpreted"])

    def test_reference_gene_id_is_not_mistaken_for_interval_gene(self):
        card = self.bundle["units"][0]["reference_material_cards"][0]
        del card["interval_gene_id"]
        card["gene_id"] = "reference_locus"
        groups, _ = collect_material(self.bundle)
        binding = next(b for g in groups for b in g["bindings"] if b["source_pointer"] == "/units/0/reference_material_cards/0")
        self.assertEqual("unit_context", binding["binding_scope"])
        self.assertEqual("reference_locus", binding["reference"]["gene_id"])

    def test_real_helper_cli_reaches_material_and_summary(self):
        self.save()
        prompts = self.root / "prompts.jsonl"
        prompts.write_text("".join(json.dumps(p) + "\n" for p in self.prompts), encoding="utf-8")
        script = str(ROOT / "scripts/agentic_route_b.py")
        first = subprocess.run([sys.executable, script, "prepare", "--prompts", str(prompts),
                                "--research-bundle", str(self.source), "--research-bundle-sha256", self.sha,
                                "--out", str(self.root / "prepare")], capture_output=True, text=True)
        self.assertEqual(0, first.returncode, first.stderr)
        manifest = json.loads((self.root / "prepare/research/manifest.json").read_text(encoding="utf-8"))
        self.assertTrue(manifest["initial_material"]["pages"])
        second = subprocess.run([sys.executable, script, "summarize-initial", "--bundle", str(self.source),
                                 "--bundle-sha256", self.sha, "--assessments", str(self.answers([self.judgment()])),
                                 "--out", str(self.root / "summary")], capture_output=True, text=True)
        self.assertEqual(0, second.returncode, second.stderr)
        self.assertEqual(4, json.loads(second.stdout)["gene_count"])


if __name__ == "__main__":
    unittest.main()
