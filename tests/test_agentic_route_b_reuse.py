from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from scripts.agentic_route_b import merge
from scripts.agentic_route_b_reuse import reuse_plan


class RouteBReusePlanTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.old = self.root / "old.jsonl"
        self.new = self.root / "new.jsonl"
        self.raw = self.root / "raw.json"
        self.response = self.root / "response.json"
        self.import_file = self.root / "import.json"
        self.rejected = self.root / "rejected.json"
        self.prompts = [{"prompt_id": f"p{i}", "route_b_unit_id": f"u{i}", "case_id": "case-x",
                         "system": "frozen rules", "user": "full material", "extra": {"sources": ["s1", "s2"]}} for i in range(3)]
        self.write_prompts(self.old, self.prompts)
        self.write_prompts(self.new, self.prompts)
        self.write(self.raw, {"case_id": "case-x", "executor": "actual-agent", "model": "actual-model",
                             "responses": [{"prompt_id": f"p{i}", "route_b_unit_id": f"u{i}", "status": "completed",
                                            "raw_text": '{"claims":[]}'} for i in range(3)]})
        self.checkpoint()

    def write(self, path, value):
        path.write_text(json.dumps(value), encoding="utf-8")

    def write_prompts(self, path, value):
        path.write_text("".join(json.dumps(row) + "\n" for row in value), encoding="utf-8")

    def sha(self, path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def checkpoint(self, rejected=None, legacy=False):
        self.response.unlink(missing_ok=True)
        self.response.with_name("response_manifest.json").unlink(missing_ok=True)
        merge(str(self.old), [str(self.raw)], str(self.response))
        self.write(self.rejected, rejected or [])
        self.write(self.import_file, {"case_id": "case-x", "stage": "semantic_proposal_import", "response_count": 3,
                                     "prompt_sha256": None if legacy else self.sha(self.old),
                                     "response_sha256": None if legacy else self.sha(self.response),
                                     "rejected_claim_row_count": len(rejected or []), "rejected_rows_sha256": self.sha(self.rejected),
                                     "outputs": {"rejected_rows": str(self.rejected)}})

    def plan(self, name="plan.json", audited=True):
        return reuse_plan(str(self.old), str(self.new), str(self.response), str(self.root / name),
                          import_manifest=str(self.import_file) if audited else "")

    def test_same_changed_new_removed_and_sources_untouched(self):
        changed = {**self.prompts[1], "system": "new rules"}
        added = {**self.prompts[2], "prompt_id": "p3", "route_b_unit_id": "u3"}
        self.write_prompts(self.new, [self.prompts[0], changed, added])
        before = {p: self.sha(p) for p in (self.old, self.new, self.raw, self.response, self.import_file, self.rejected)}
        plan = self.plan()
        self.assertEqual({"same": 1, "changed": 1, "new": 1, "removed": 1}, plan["comparison_counts"])
        self.assertEqual((1, 2, 0), (plan["reusable_count"], plan["new_execution_required_count"], plan["new_model_execution_count"]))
        self.assertEqual(before, {p: self.sha(p) for p in before})
        self.assertEqual("actual-agent", plan["rows"][0]["raw_response"]["executor"])
        self.assertEqual(self.sha(self.raw), plan["rows"][0]["raw_response"]["sha256"])

    def test_every_prompt_field_changes_digest(self):
        for index, change in enumerate(({"system": "different"}, {"user": "different"}, {"case_id": "other"},
                                        {"extra": {"sources": ["s2", "s1"]}}, {"novel_field": True}, {"route_b_unit_id": "other"})):
            self.write_prompts(self.new, [{**p, **change, **({"route_b_unit_id": f"other-{n}"} if "route_b_unit_id" in change else {})}
                                         for n, p in enumerate(self.prompts)])
            with self.subTest(change=change):
                self.assertEqual(3, self.plan(f"change-{index}.json")["comparison_counts"]["changed"])
        self.write_prompts(self.old, [{**p, "typed": True} for p in self.prompts])
        self.checkpoint()
        self.write_prompts(self.new, [{**p, "typed": 1} for p in self.prompts])
        self.assertEqual(0, self.plan("types.json")["reusable_count"])

    def test_key_order_does_not_change_input(self):
        self.write_prompts(self.new, [dict(reversed(list(p.items()))) for p in reversed(self.prompts)])
        self.assertEqual(3, self.plan()["reusable_count"])

    def test_nonstring_ids_and_missing_case_are_not_canonical(self):
        for field, value in (("prompt_id", True), ("route_b_unit_id", 1), ("case_id", True), ("case_id", "")):
            self.write_prompts(self.new, [{**self.prompts[0], field: value}])
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                self.plan()

    def test_missing_or_unbound_import_audit_is_not_eligible(self):
        self.assertEqual(0, self.plan(audited=False)["reusable_count"])
        self.checkpoint(legacy=True)
        plan = self.plan("legacy.json")
        self.assertEqual(0, plan["reusable_count"])
        self.assertEqual({"import_audit_input_not_bound"}, {row["reason"] for row in plan["rows"]})

    def test_known_reject_and_healthy_siblings(self):
        claim = {"claim_type": "mechanistic_review_seed", "evidence_scope": "family"}
        raw = json.loads(self.raw.read_text())
        raw["responses"][1]["raw_text"] = json.dumps({"claims": [claim]})
        self.write(self.raw, raw)
        rejected = [{"prompt_id": "p1", "unit_id": "u1", "claim_index": 1, "raw_claim": claim, "reason": "missing family"}]
        self.checkpoint(rejected)
        plan = self.plan()
        self.assertEqual(2, plan["reusable_count"])
        self.assertEqual("known_import_rejected", plan["rows"][1]["reason"])
        self.checkpoint(rejected, legacy=True)
        self.assertEqual("known_import_rejected", self.plan("legacy.json")["rows"][1]["reason"])

    def test_invalid_partial_and_failed_keep_healthy_response(self):
        raw = json.loads(self.raw.read_text())
        raw["responses"][1]["raw_text"] = '{"claims":[{"claim_type":"valid"},{}]}'
        raw["responses"][2] = {"prompt_id": "p2", "status": "failed", "reason": "provider unavailable"}
        self.write(self.raw, raw)
        self.checkpoint()
        plan = self.plan()
        self.assertEqual(1, plan["reusable_count"])
        self.assertEqual("response_incomplete_or_invalid", plan["rows"][1]["reason"])
        self.assertEqual("response_incomplete_or_invalid", plan["rows"][2]["reason"])

    def test_raw_tamper_rejected_before_plan_written(self):
        self.raw.write_text(self.raw.read_text() + " ", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "raw response hash mismatch"):
            self.plan()
        self.assertFalse((self.root / "plan.json").exists())

    def test_failed_null_raw_keeps_healthy_rows(self):
        raw = json.loads(self.raw.read_text())
        raw["responses"][1] = {"prompt_id": "p1", "status": "failed", "reason": "provider unavailable", "raw_text": None}
        self.write(self.raw, raw)
        self.checkpoint()
        plan = self.plan()
        self.assertEqual(2, plan["reusable_count"])
        self.assertEqual("", plan["rows"][1]["raw_text_sha256"])
        self.assertEqual(self.sha(self.raw), plan["rows"][1]["raw_response"]["sha256"])

    def test_checkpoint_tamper_even_with_resealed_hash_is_rejected(self):
        response = json.loads(self.response.read_text())
        response["candidates"][0]["attempt_model"] = "invented"
        self.write(self.response, response)
        path = self.response.with_name("response_manifest.json")
        manifest = json.loads(path.read_text())
        manifest["response_sha256"] = self.sha(self.response)
        self.write(path, manifest)
        with self.assertRaisesRegex(ValueError, "original raw response/provenance"):
            self.plan()

    def test_rejected_evidence_must_match_frozen_claim(self):
        self.checkpoint([{"prompt_id": "p0", "unit_id": "u0", "claim_index": 1, "raw_claim": {}, "reason": "bad"}])
        with self.assertRaisesRegex(ValueError, "match old raw claim"):
            self.plan()

    def test_cli_and_existing_output_refusal(self):
        process = subprocess.run([sys.executable, "scripts/agentic_route_b.py", "reuse-plan", "--old-prompts", str(self.old),
                                  "--new-prompts", str(self.new), "--old-response", str(self.response),
                                  "--import-manifest", str(self.import_file), "--out", str(self.root / "plan.json")], capture_output=True, text=True)
        self.assertEqual(0, process.returncode, process.stderr)
        self.assertEqual(3, json.loads(process.stdout)["reusable_count"])
        with self.assertRaisesRegex(ValueError, "already exists"):
            self.plan()


if __name__ == "__main__":
    unittest.main()
