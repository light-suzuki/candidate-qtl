"""Offline research navigation checks, not model/Type3 performance tests."""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from scripts import agentic_route_b as route_b
from scripts.agentic_route_b_research import TASK, plan_research_view

ROOT = Path(__file__).resolve().parents[1]


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def fixture() -> tuple[dict, list[dict]]:
    units = [
        {"unit_id": "u1", "member_interval_genes": ["g1", "g2"],
         "identity": {"family_ids": ["f1"]}, "request_state": "requested",
         "source_material_class": "annotation_context",
         "genes": [{"gene_id": "g1", "reference_gene": "ref1", "reference_function": "fixture function"},
                   {"gene_id": "g2", "reference_gene": "ref2", "reference_function": "different fixture function"}],
         "reference_material_cards": [{"interval_gene_id": "g2", "reference_gene": "ref3", "function_description": "alternative record"}]},
        {"unit_id": "u2", "member_interval_genes": ["g3"], "identity": {},
         "request_state": "not_requested_no_reviewable_material", "source_material_class": "retrieval_handles_only",
         "genes": [{"gene_id": "g3", "sequence_similarity_hits": [
             {"hit_accession": "acc1", "hit_species": "plant1", "hit_protein_name": "retrieval name"},
             {"hit_accession": "acc2", "hit_species": "plant2"}]}]},
        {"unit_id": "u3", "member_interval_genes": ["g4"], "identity": {},
         "request_state": "not_requested_no_reviewable_material", "source_material_class": "no_reviewable_material",
         "genes": [{"gene_id": "g4"}]},
    ]
    prompts = [{"case_id": "fixture", "prompt_id": "route-b-semantic-u1", "route_b_unit_id": "u1",
                "system": "TEST ONLY", "user": "Phenotype: test\n\n" + json.dumps(units[0], ensure_ascii=False, sort_keys=True)}]
    bundle = {"case_id": "fixture", "route": "interval_first", "stage": "post_b0_semantic_bundle",
              "phenotype": "test", "units": units, "interval_gene_count": 4, "semantic_unit_count": 3}
    seal_prompts(bundle, prompts)
    return bundle, prompts


def seal_prompts(bundle: dict, prompts: list[dict]) -> None:
    bundle["execution_prompt_count"] = len(prompts)
    bundle["execution_prompt_fingerprints"] = [
        digest(json.dumps(p, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()) for p in prompts
    ]


class ResearchViewTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.bundle, self.prompts = fixture()

    def save(self):
        self.source = self.root / "bundle.json"
        self.source.write_text(json.dumps(self.bundle, ensure_ascii=False), encoding="utf-8")
        self.prompts_path = self.root / "prompts.jsonl"
        self.prompts_path.write_text("".join(json.dumps(p, ensure_ascii=False) + "\n" for p in self.prompts), encoding="utf-8")
        return digest(self.source.read_bytes())

    def plan(self, budget=262144):
        sha = self.save()
        return plan_research_view(str(self.source), sha, self.prompts, "fixture", budget)

    @staticmethod
    def rows(files):
        return [json.loads(line) for name, data in sorted(files.items())
                if name.endswith(".jsonl") for line in data.splitlines()]

    def test_all_units_and_copies_visible_not_only_requested(self):
        files, metadata = self.plan()
        rows = self.rows(files)
        self.assertEqual(["u1", "u2", "u3"], [r["unit_id"] for r in rows])
        self.assertEqual(["g1", "g2", "g3", "g4"], [g for r in rows for g in r["member_interval_genes"]])
        self.assertEqual(4, metadata["gene_count"])
        self.assertEqual("not_started", metadata["execution_state"])
        manifest = json.loads(files["research/manifest.json"])
        self.assertEqual(2, manifest["not_requested_unit_count"])
        self.assertEqual(["acc1", "acc2"], [h["hit_accession"] for h in rows[1]["genes"][0]["sequence_handles"]])
        self.assertEqual("retrieval name", rows[1]["genes"][0]["sequence_handles"][0]["hit_protein_name"])
        self.assertTrue(all(r["view_status"] == "navigation_only_not_evaluated" for r in rows))

    def test_no_formal_prompts_still_has_research_roster(self):
        self.prompts = []
        self.bundle["units"][0]["request_state"] = "not_requested_no_reviewable_material"
        seal_prompts(self.bundle, self.prompts)
        files, metadata = self.plan()
        self.assertEqual(3, len(self.rows(files)))
        self.assertEqual(4, metadata["gene_count"])

    def test_empty_interval(self):
        self.bundle.update(units=[], interval_gene_count=0, semantic_unit_count=0)
        self.prompts = []
        seal_prompts(self.bundle, self.prompts)
        files, metadata = self.plan()
        self.assertEqual([], self.rows(files))
        self.assertEqual(0, metadata["page_count"])
        self.assertIn("research/task.md", files)

    def test_preview_omissions_and_atlas_identity_are_explicit(self):
        gene = self.bundle["units"][0]["genes"][0]
        gene["reference_function"] = "長" * 600
        gene["atlas_material_context"] = {"reference_gene": "other", "reference_function": "atlas only"}
        before = copy.deepcopy(self.bundle)
        files, _ = self.plan()
        view = self.rows(files)[0]["genes"][0]
        self.assertEqual("ref1", view["reference_gene"])
        self.assertEqual("other", view["atlas_context_preview"]["reference_gene"])
        self.assertEqual(["reference_function"], view["truncated_fields"])
        self.assertEqual(240, len(view["cues"]["reference_function"]))
        self.assertEqual("/units/0/genes/0", view["gene_pointer"])
        self.assertEqual(before, json.loads(self.source.read_text()))

    def test_oversized_unit_is_preserved(self):
        files, _ = self.plan(budget=1)
        manifest = json.loads(files["research/manifest.json"])
        self.assertEqual(3, len(manifest["pages"]))
        self.assertEqual(4, sum(len(r["member_interval_genes"]) for r in self.rows(files)))
        for page in manifest["pages"]:
            self.assertEqual(1, page["unit_count"])
            self.assertGreater(page["over_bound_bytes"], 0)
            self.assertEqual(page["sha256"], digest(files[page["file"]]))

    def test_1000_gene_roster_has_no_top_k(self):
        self.bundle["units"] = [{"unit_id": f"u{i}", "member_interval_genes": [f"g{i}"],
                                 "genes": [{"gene_id": f"g{i}"}],
                                 "request_state": "not_requested_no_reviewable_material"} for i in range(1000)]
        self.bundle.update(interval_gene_count=1000, semantic_unit_count=1000)
        self.prompts = []
        seal_prompts(self.bundle, self.prompts)
        files, meta = self.plan(budget=4096)
        self.assertEqual(1000, len(self.rows(files)))
        self.assertEqual(1000, meta["gene_count"])
        self.assertTrue(all(p["bytes"] <= 4096 for p in json.loads(files["research/manifest.json"])["pages"]))

    def test_integrity_errors(self):
        mutations = {
            "case": lambda b: b.update(case_id="another"),
            "stage": lambda b: b.update(stage="other"),
            "phenotype": lambda b: b.update(phenotype=""),
            "fingerprint": lambda b: b.update(execution_prompt_fingerprints=["bad"]),
            "count": lambda b: b.update(interval_gene_count=5),
            "bool_count": lambda b: b.update(execution_prompt_count=True),
            "unit_count": lambda b: b.update(semantic_unit_count=99),
            "unit_duplicate": lambda b: b["units"][1].update(unit_id="u1"),
            "member_duplicate": lambda b: b["units"][1].update(member_interval_genes=["g1"]),
            "gene_mismatch": lambda b: b["units"][1]["genes"][0].update(gene_id="unbound"),
            "gene_type": lambda b: b["units"][0]["genes"][0].update(gene_id=None),
            "requested": lambda b: b["units"][1].update(request_state="requested"),
            "bad_hits": lambda b: b["units"][1]["genes"][0].update(sequence_similarity_hits=["invalid"]),
        }
        for name, mutate in mutations.items():
            with self.subTest(name=name):
                self.bundle, self.prompts = fixture()
                mutate(self.bundle)
                with self.assertRaises(ValueError):
                    self.plan()

    def test_tampered_source_or_prompts(self):
        sha = self.save()
        self.source.write_bytes(self.source.read_bytes() + b" ")
        with self.assertRaisesRegex(ValueError, "SHA-256"):
            plan_research_view(str(self.source), sha, self.prompts, "fixture", 4096)
        sha = self.save()
        self.prompts[0]["user"] += "another question"
        with self.assertRaisesRegex(ValueError, "execution prompts"):
            plan_research_view(str(self.source), sha, self.prompts, "fixture", 4096)

    def test_prepare_reaches_view_without_changing_canonical_artifacts(self):
        sha = self.save()
        before = self.source.read_bytes(), self.prompts_path.read_bytes()
        default, research = self.root / "default", self.root / "research"
        baseline = route_b.prepare(str(self.prompts_path), str(default), case_id="fixture")
        actual = route_b.prepare(str(self.prompts_path), str(research), case_id="fixture",
                                 research_bundle=str(self.source), research_bundle_sha256=sha)
        self.assertNotIn("research_view", baseline)
        self.assertEqual(baseline, {k: v for k, v in actual.items() if k != "research_view"})
        for file in default.rglob("*"):
            if file.is_file() and file.name != "prepare_manifest.json":
                self.assertEqual(file.read_bytes(), (research / file.relative_to(default)).read_bytes())
        self.assertEqual(before, (self.source.read_bytes(), self.prompts_path.read_bytes()))
        self.assertTrue((research / "research/task.md").is_file())

    def test_invalid_opt_in_does_not_create_output(self):
        sha = self.save()
        for kwargs in ({"research_bundle": str(self.source)}, {"research_bundle_sha256": sha},
                       {"research_bundle": str(self.source), "research_bundle_sha256": "bad"}):
            with self.subTest(kwargs=kwargs):
                out = self.root / "invalid"
                with self.assertRaises(ValueError):
                    route_b.prepare(str(self.prompts_path), str(out), case_id="fixture", **kwargs)
                self.assertFalse(out.exists())

    def test_prepare_cli(self):
        sha = self.save()
        out = self.root / "cli"
        result = subprocess.run([sys.executable, str(ROOT / "scripts/agentic_route_b.py"), "prepare",
                                 "--prompts", str(self.prompts_path), "--out", str(out), "--case-id", "fixture",
                                 "--research-bundle", str(self.source), "--research-bundle-sha256", sha],
                                capture_output=True, text=True, check=True)
        manifest = json.loads(result.stdout)
        self.assertEqual(4, manifest["research_view"]["gene_count"])
        self.assertTrue(manifest["reading_view"]["exact_roundtrip_verified"])

    def test_task_is_research_not_an_unverified_claim_import(self):
        for text in ("具体的機構を既に説明できることは検索の開始条件ではない", "Type3", "新しいprepare出力をfreeze",
                     "まず区間gene/copy", "中立に判断する", "候補化を前提に", "無理に候補化せず",
                     "探索メモをcanonical responses", "全geneの定型ミニレビュー", "A出力", "hashを再確認"):
            self.assertIn(text, TASK)


if __name__ == "__main__":
    unittest.main()
