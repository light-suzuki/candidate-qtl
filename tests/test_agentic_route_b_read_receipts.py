"""Read receipts let merge tell an unread no-claim unit from a reviewed one."""
import json
from pathlib import Path
import tempfile
import unittest

from scripts.agentic_route_b import merge, preflight, prepare

from test_agentic_route_b_material import semantic_prompt

CLAIM = {"claims": [{"claim_type": "mechanistic_review_seed", "biological_function": "f",
                     "mechanism": "m", "phenotype_relation": "r", "reference_gene": "ref-1",
                     "reference_species": "Species one", "evidence_scope": "orthogroup",
                     "source_ids": ["PMID:1"], "alternative_hypothesis": "a",
                     "uncertainty_note": "u", "member_interval_genes": ["gene-1"]}]}


class ReadReceiptTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        source = self.root / "source.jsonl"
        self.rows = [semantic_prompt(1), semantic_prompt(2)]
        source.write_text("".join(json.dumps(r) + "\n" for r in self.rows), encoding="utf-8")
        self.prompts = str(source)
        self.manifest = prepare(self.prompts, str(self.root / "prepared"), max_prompts=2)
        self.reading_manifest = str(self.root / "prepared" / "prepare_manifest.json")
        self.shard = self.manifest["reading_view"]["shards"][0]

    def _receipt_of(self, display_file: str) -> str:
        return json.loads((self.root / "prepared" / display_file).read_text(encoding="utf-8"))["read_receipt"]

    def _response(self, name: str, records: list[dict]) -> str:
        path = self.root / name
        path.write_text(json.dumps({"case_id": "case-x", "executor": "e", "model": "m",
                                    "responses": records}), encoding="utf-8")
        return str(path)

    def _record(self, pid: str, claims: dict, receipts: list[str] | None) -> dict:
        record = {"prompt_id": pid, "route_b_unit_id": pid.replace("prompt", "unit"),
                  "raw_text": json.dumps(claims)}
        if receipts is not None:
            record["read_receipts"] = receipts
        return record

    def test_receipt_manifest_must_bind_the_exact_prompt_roster(self):
        response = self._response("empty.json", [self._record("prompt-1", {"claims": []}, [])])
        path = Path(self.reading_manifest)
        original = path.read_text(encoding="utf-8")
        for change in ("case", "hash", "missing_prompt", "missing_token"):
            data = json.loads(original)
            if change == "case":
                data["case_id"] = "other-case"
            elif change == "hash":
                data["source_sha256"] = "0" * 64
            elif change == "missing_prompt":
                del data["reading_view"]["shards"][0]["read_receipts_by_prompt"]["prompt-1"]
            else:
                data["reading_view"]["shards"][0]["read_receipts_by_prompt"]["prompt-1"].pop()
            path.write_text(json.dumps(data), encoding="utf-8")
            with self.subTest(change=change), self.assertRaises(ValueError):
                merge(self.prompts, [response], str(self.root / f"bad-{change}.json"), reading_manifest=str(path))
        path.write_text(original, encoding="utf-8")

    def test_receipt_is_the_last_field_of_every_display_and_bound_per_prompt(self):
        view = self.shard["display_view"]
        for descriptor in [view["assignment"], *view["materials"]]:
            text = (self.root / "prepared" / descriptor["file"]).read_text(encoding="utf-8")
            self.assertIn(descriptor["read_receipt"], text.splitlines()[-2])
            self.assertEqual(descriptor["read_receipt"], self._receipt_of(descriptor["file"]))
        required = self.shard["read_receipts_by_prompt"]
        self.assertEqual({"prompt-1", "prompt-2"}, set(required))
        self.assertIn(view["assignment"]["read_receipt"], required["prompt-1"])

    def test_unread_no_claim_answer_stays_not_evaluated_and_claims_are_flagged(self):
        required = self.shard["read_receipts_by_prompt"]
        response = self._response("r.json", [
            self._record("prompt-1", {"claims": []}, required["prompt-1"]),
            self._record("prompt-2", {"claims": []}, required["prompt-2"][:1]),
        ])
        result = merge(self.prompts, [response], str(self.root / "merged.json"),
                       reading_manifest=self.reading_manifest)
        roster = {row["prompt_id"]: row for row in result["response_roster"]}
        self.assertEqual("completed", roster["prompt-1"]["completion_status"])
        self.assertEqual("not_evaluated", roster["prompt-2"]["completion_status"])
        self.assertEqual("read_receipt_incomplete", roster["prompt-2"]["not_evaluated_reason"])
        self.assertEqual(1, result["read_receipts"]["demoted_to_not_evaluated"])

        flagged = self._response("c.json", [self._record("prompt-1", CLAIM, []),
                                            self._record("prompt-2", {"claims": []}, required["prompt-2"])])
        claims = merge(self.prompts, [flagged], str(self.root / "claims.json"),
                       reading_manifest=self.reading_manifest)
        self.assertEqual(1, claims["read_receipts"]["incomplete_with_claims"])
        self.assertEqual(2, claims["completed_count"])

    def test_merge_without_reading_manifest_keeps_previous_behavior(self):
        response = self._response("r.json", [self._record("prompt-1", {"claims": []}, None),
                                             self._record("prompt-2", {"claims": []}, None)])
        result = merge(self.prompts, [response], str(self.root / "merged.json"))
        self.assertEqual(2, result["completed_count"])
        self.assertNotIn("read_receipts", result)

    def test_preflight_accepts_cli_provenance_and_reports_receipts(self):
        path = self.root / "noprov.json"
        path.write_text(json.dumps({"case_id": "case-x", "responses": [
            self._record("prompt-1", {"claims": []}, None)]}), encoding="utf-8")
        self.assertEqual("rejected", preflight(self.prompts, [str(path)])["status"])
        checked = preflight(self.prompts, [str(path)], executor="e", model="m",
                            reading_manifest=self.reading_manifest)
        self.assertEqual("accepted", checked["status"])
        self.assertEqual("incomplete", checked["responses"][0]["read_receipt_state"])
        self.assertFalse(checked["responses"][0]["completed"])


if __name__ == "__main__":
    unittest.main()
