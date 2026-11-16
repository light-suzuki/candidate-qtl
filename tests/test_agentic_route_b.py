from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from scripts.agentic_route_b import merge, prepare


class AgenticRouteBTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.prompts = self.root / "route_b_semantic_prompts.jsonl"
        self.rows = [
            {"prompt_id": "p1", "route_b_unit_id": "u1", "case_id": "case-x", "user": "unit 1"},
            {"prompt_id": "p2", "route_b_unit_id": "u2", "case_id": "case-x", "user": "unit 2"},
            {"prompt_id": "p3", "route_b_unit_id": "u3", "case_id": "case-x", "user": "unit 3"},
        ]
        self.prompts.write_text("".join(json.dumps(r) + "\n" for r in self.rows), encoding="utf-8")

    def tearDown(self):
        self.tmp.cleanup()

    def test_prepare_shards_keep_all_prompt_and_unit_ids(self):
        manifest = prepare(str(self.prompts), str(self.root / "shards"), max_prompts=2)
        self.assertEqual(3, manifest["prompt_count"])
        self.assertEqual([["p1", "p2"], ["p3"]], [s["prompt_ids"] for s in manifest["shards"]])
        recovered = []
        for shard in manifest["shards"]:
            recovered.extend(json.loads(line) for line in (self.root / "shards" / shard["file"]).read_text().splitlines())
        self.assertEqual(["p1", "p2", "p3"], [r["prompt_id"] for r in recovered])
        self.assertEqual(["u1", "u2", "u3"], [r["route_b_unit_id"] for r in recovered])

    def test_merge_fills_untouched_as_not_evaluated_and_keeps_raw(self):
        partial = self.root / "partial.json"
        raw = '{"claims":[]}'
        partial.write_text(json.dumps({"case_id": "case-x", "executor": "Codex", "model": "example", "response_id": "r1", "responses": [
            {"prompt_id": "p1", "raw_text": raw},
            {"prompt_id": "p2", "status": "failed", "reason": "agent stopped"},
        ]}), encoding="utf-8")
        result = merge(str(self.prompts), [str(partial)], str(self.root / "responses.json"))
        roster = json.loads((self.root / "responses.json").read_text())["candidates"]
        self.assertEqual(["json", "not_evaluated", "not_evaluated"], [r["parse_status"] for r in roster])
        self.assertEqual(raw, roster[0]["raw_text"])
        self.assertEqual("agent stopped", roster[1]["not_evaluated_reason"])
        self.assertEqual("no_partial_response", roster[2]["not_evaluated_reason"])
        self.assertEqual("partial", result["completion_status"])
        self.assertEqual("Codex", result["partial_response_provenance"][0]["executor"])

    def test_empty_prompt_roster_is_explicit_not_evaluated(self):
        empty = self.root / "empty-prompts.jsonl"
        empty.write_text("", encoding="utf-8")
        out = self.root / "empty-merged.json"
        manifest = merge(str(empty), [], str(out), case_id="case-empty")
        self.assertEqual("not_evaluated", manifest["completion_status"])
        self.assertEqual(0, manifest["prompt_count"])
        self.assertEqual([], json.loads(out.read_text())["candidates"])

    def test_merge_writes_execution_contract_and_rejects_overwrite(self):
        partial = self.root / "contract-partial.json"
        partial.write_text(json.dumps({
            "case_id": "case-x", "executor": "Codex", "model": "example",
            "responses": [{"prompt_id": "p1", "raw_text": '{"claims":[]}'}],
        }), encoding="utf-8")
        out = self.root / "contract-responses.json"
        manifest = merge(str(self.prompts), [str(partial)], str(out))
        self.assertEqual("agentic-route-b-execution-1.1.0", manifest["contract"])
        self.assertEqual(1, manifest["completed_count"])
        self.assertEqual(2, manifest["not_evaluated_count"])
        with self.assertRaisesRegex(ValueError, "output pair already exists"):
            merge(str(self.prompts), [str(partial)], str(out))
        out.unlink()
        with self.assertRaisesRegex(ValueError, "output pair already exists"):
            merge(str(self.prompts), [str(partial)], str(out))

    def test_partial_merge_can_resume_with_a_new_attempt(self):
        first = self.root / "first-response.json"
        first.write_text(json.dumps({
            "case_id": "case-x", "executor": "Codex", "model": "example",
            "responses": [
                {"prompt_id": "p1", "raw_text": '{"claims":[]}'},
                {"prompt_id": "p2", "raw_text": '{"claims":[]}'},
            ],
        }), encoding="utf-8")
        retry = self.root / "retry-response.json"
        retry.write_text(json.dumps({
            "case_id": "case-x", "executor": "Codex", "model": "example",
            "responses": [{"prompt_id": "p3", "raw_text": '{"claims":[]}'}],
        }), encoding="utf-8")
        first_out = self.root / "first-merged.json"
        first_manifest = merge(str(self.prompts), [str(first)], str(first_out))
        self.assertEqual("partial", first_manifest["completion_status"])
        resumed_out = self.root / "resumed-merged.json"
        resumed = merge(
            str(self.prompts), [str(first_out), str(retry)], str(resumed_out),
            executor="Codex", model="example",
        )
        self.assertEqual("complete", resumed["completion_status"])
        self.assertEqual(3, resumed["completed_count"])
        self.assertEqual(0, resumed["not_evaluated_count"])
        self.assertTrue(first_out.is_file())
        self.assertTrue(resumed_out.is_file())

    def test_explicit_failed_response_can_be_retried_without_losing_healthy_results(self):
        first = self.root / "failed-first.json"
        first.write_text(json.dumps({
            "case_id": "case-x", "executor": "Codex", "model": "first-model",
            "response_id": "attempt-1",
            "responses": [
                {"prompt_id": "p1", "raw_text": '{"claims":[]}'},
                {"prompt_id": "p2", "status": "failed", "reason": "shard timeout"},
            ],
        }), encoding="utf-8")
        first_out = self.root / "failed-merged.json"
        first_manifest = merge(str(self.prompts), [str(first)], str(first_out))
        self.assertEqual("partial", first_manifest["completion_status"])
        self.assertEqual(1, first_manifest["completed_count"])

        retry = self.root / "failed-retry.json"
        retry.write_text(json.dumps({
            "case_id": "case-x", "executor": "Codex-retry", "model": "retry-model",
            "response_id": "attempt-2",
            "responses": [
                {"prompt_id": "p2", "raw_text": '{"claims":[]}'},
                {"prompt_id": "p3", "raw_text": '{"claims":[]}'},
            ],
        }), encoding="utf-8")
        resumed_out = self.root / "failed-resumed.json"
        resumed = merge(
            str(self.prompts), [str(first_out), str(retry)], str(resumed_out),
        )
        roster = json.loads(resumed_out.read_text())["candidates"]
        self.assertEqual("complete", resumed["completion_status"])
        self.assertEqual(["json", "json", "json"], [row["parse_status"] for row in roster])
        self.assertEqual("Codex-retry", roster[1]["attempt_executor"])
        self.assertEqual(1, len(resumed["retry_history"]))
        self.assertEqual("p2", resumed["retry_history"][0]["prompt_id"])
        self.assertEqual("shard timeout", resumed["retry_history"][0]["replaced_reason"])
        self.assertEqual(str(first.resolve()), resumed["retry_history"][0]["replaced_response_file"])
        self.assertEqual(hashlib.sha256(first.read_bytes()).hexdigest(), resumed["retry_history"][0]["replaced_response_sha256"])
        self.assertEqual(str(retry.resolve()), resumed["retry_history"][0]["replacement_response_file"])
        self.assertEqual(hashlib.sha256(retry.read_bytes()).hexdigest(), resumed["retry_history"][0]["replacement_response_sha256"])
        provenance = resumed["partial_response_provenance"]
        self.assertEqual({str(first.resolve()), str(retry.resolve())}, {row["file"] for row in provenance})
        self.assertEqual([0, 1], sorted(row["attempt_index"] for row in provenance))
        response_ids = {row["file"]: row["response_id"] for row in provenance}
        self.assertEqual("attempt-1", response_ids[str(first.resolve())])
        self.assertEqual("attempt-2", response_ids[str(retry.resolve())])

    def test_retry_cannot_duplicate_within_one_attempt_or_downgrade_completed(self):
        duplicate = self.root / "same-attempt-duplicate.json"
        duplicate.write_text(json.dumps({
            "case_id": "case-x", "responses": [
                {"prompt_id": "p1", "status": "failed", "reason": "timeout"},
                {"prompt_id": "p1", "raw_text": '{"claims":[]}'},
            ],
        }), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "duplicate prompt_id response in one attempt"):
            merge(str(self.prompts), [str(duplicate)], str(self.root / "same-attempt-out.json"), executor="TEST", model="TEST")

        completed = self.root / "completed.json"
        completed.write_text(json.dumps({
            "case_id": "case-x", "executor": "completed-executor", "model": "completed-model",
            "responses": [{"prompt_id": "p1", "raw_text": '{"claims":[]}'}],
        }), encoding="utf-8")
        failed = self.root / "late-failure.json"
        failed.write_text(json.dumps({
            "case_id": "case-x", "executor": "failed-executor", "model": "failed-model",
            "responses": [{"prompt_id": "p1", "status": "failed", "reason": "late timeout"}],
        }), encoding="utf-8")
        forward = self.root / "downgrade-forward.json"
        reverse = self.root / "downgrade-reverse.json"
        forward_result = merge(str(self.prompts), [str(completed), str(failed)], str(forward))
        reverse_result = merge(str(self.prompts), [str(failed), str(completed)], str(reverse))
        self.assertEqual("partial", forward_result["completion_status"])
        self.assertEqual("partial", reverse_result["completion_status"])
        self.assertEqual("json", json.loads(forward.read_text())["candidates"][0]["parse_status"])
        self.assertEqual("json", json.loads(reverse.read_text())["candidates"][0]["parse_status"])
        self.assertEqual(
            json.loads(forward.read_text())["candidates"][0]["attempt_executor"],
            json.loads(reverse.read_text())["candidates"][0]["attempt_executor"],
        )
        self.assertEqual("completed-executor", json.loads(forward.read_text())["candidates"][0]["attempt_executor"])

    def test_merge_rejects_unknown_and_duplicate_ids(self):
        for records, message in [
            ([{"prompt_id": "unknown", "raw_text": '{"claims":[]}'}], "unknown prompt_id"),
            ([{"prompt_id": "p1", "raw_text": '{"claims":[]}'}, {"prompt_id": "p1", "raw_text": '{"claims":[]}'}], "duplicate prompt_id"),
        ]:
            partial = self.root / "bad.json"
            partial.write_text(json.dumps({"case_id": "case-x", "responses": records}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, message):
                merge(str(self.prompts), [str(partial)], str(self.root / "out.json"), executor="TEST", model="TEST")

    def test_merge_rejects_case_mismatch_and_contradictory_model_json(self):
        partial = self.root / "case.json"
        partial.write_text(json.dumps({"case_id": "other", "responses": []}), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "case_id mismatch"):
            merge(str(self.prompts), [str(partial)], str(self.root / "out.json"))
        partial.write_text(json.dumps({
            "case_id": "case-x",
            "responses": [{
                "prompt_id": "p1",
                "raw_text": '{"claims":[]}',
                "model_json": {"claims": [{"claim_type": "unbound"}]},
            }],
        }), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "differs from raw_text"):
            merge(str(self.prompts), [str(partial)], str(self.root / "mismatch-out.json"), executor="TEST", model="TEST")

    def test_malformed_model_response_is_local_and_retry_preserves_healthy_rows_and_raw(self):
        first = self.root / "malformed-first.json"
        malformed_raw = ' {"claims": [\n broken 日本語\n'
        original_bad = {"prompt_id": "p2", "status": "completed", "raw_text": malformed_raw, "model_json": {"claims": []}, "request_id": "request-bad"}
        healthy_raw = '{"claims":[{"claim_type":"function","claim_text":"original healthy claim"}]}'
        first.write_text(json.dumps({
            "case_id": "case-x", "executor": "first-executor", "model": "first-model", "response_id": "attempt-bad",
            "responses": [{"prompt_id": "p1", "raw_text": healthy_raw}, original_bad, {"prompt_id": "p3", "raw_text": '{"claims":[]}'}],
        }, ensure_ascii=False), encoding="utf-8")
        first_bytes = first.read_bytes()
        checkpoint = self.root / "malformed-merged.json"
        initial = merge(str(self.prompts), [str(first)], str(checkpoint))
        initial_rows = json.loads(checkpoint.read_text(encoding="utf-8"))["candidates"]
        self.assertEqual("partial", initial["completion_status"])
        self.assertEqual(2, initial["completed_count"])
        failed = initial_rows[1]
        self.assertEqual("not_evaluated", failed["parse_status"])
        self.assertEqual("not_evaluated", failed["review_status"])
        self.assertEqual("model_response_invalid_json", failed["not_evaluated_reason"])
        self.assertIsNone(failed["model_json"])
        self.assertEqual([], failed["claims"])
        self.assertEqual(malformed_raw, failed["raw_text"])
        self.assertEqual(original_bad, failed["agent_response_raw"])
        self.assertEqual(hashlib.sha256(malformed_raw.encode("utf-8")).hexdigest(), initial["response_roster"][1]["raw_text_sha256"])

        retry = self.root / "malformed-retry.json"
        retry.write_text(json.dumps({
            "case_id": "case-x", "executor": "retry-executor", "model": "retry-model", "response_id": "attempt-good",
            "responses": [{"prompt_id": "p2", "raw_text": '{"claims":[]}'}],
        }), encoding="utf-8")
        resumed_out = self.root / "malformed-resumed.json"
        resumed = merge(str(self.prompts), [str(checkpoint), str(retry)], str(resumed_out))
        resumed_rows = json.loads(resumed_out.read_text(encoding="utf-8"))["candidates"]
        self.assertEqual("complete", resumed["completion_status"])
        self.assertEqual(3, resumed["completed_count"])
        self.assertEqual(0, resumed["not_evaluated_count"])
        self.assertEqual(initial_rows[0], resumed_rows[0])
        self.assertEqual(initial_rows[2], resumed_rows[2])
        self.assertEqual("retry-executor", resumed_rows[1]["attempt_executor"])
        history = resumed["retry_history"]
        self.assertEqual(1, len(history))
        self.assertEqual("p2", history[0]["prompt_id"])
        self.assertEqual("model_response_invalid_json", history[0]["replaced_reason"])
        self.assertEqual(malformed_raw, history[0]["replaced_raw_text"])
        self.assertEqual(original_bad, history[0]["replaced_agent_response_raw"])
        self.assertEqual(hashlib.sha256(malformed_raw.encode("utf-8")).hexdigest(), history[0]["replaced_raw_text_sha256"])
        self.assertEqual(str(first.resolve()), history[0]["replaced_response_file"])
        self.assertEqual(hashlib.sha256(first_bytes).hexdigest(), history[0]["replaced_response_sha256"])
        self.assertEqual(str(retry.resolve()), history[0]["replacement_response_file"])
        self.assertEqual(first_bytes, first.read_bytes())
        self.assertEqual({"attempt-bad", "attempt-good"}, {row["response_id"] for row in resumed["partial_response_provenance"]})

    def test_bad_claims_envelope_and_members_are_local_format_failures(self):
        for index, (raw, reason) in enumerate([
            ('[]', "model_response_invalid_claims_envelope"),
            ('null', "model_response_invalid_claims_envelope"),
            ('{}', "model_response_invalid_claims_envelope"),
            ('{"claims":null}', "model_response_invalid_claims_envelope"),
            ('{"claims":"unparsed"}', "model_response_invalid_claims_envelope"),
            ('{"claims":["unparsed"]}', "model_response_invalid_claims_members"),
        ]):
            with self.subTest(raw=raw):
                partial = self.root / f"bad-envelope-{index}.json"
                partial.write_text(json.dumps({
                    "case_id": "case-x", "executor": "TEST", "model": "TEST",
                    "responses": [{"prompt_id": "p1", "raw_text": raw, "model_json": json.loads(raw), "claims": []}, {"prompt_id": "p2", "raw_text": '{"claims":[]}'}],
                }), encoding="utf-8")
                out = self.root / f"bad-envelope-out-{index}.json"
                result = merge(str(self.prompts), [str(partial)], str(out))
                rows = json.loads(out.read_text())["candidates"]
                self.assertEqual(1, result["completed_count"])
                self.assertEqual(reason, rows[0]["not_evaluated_reason"])
                self.assertEqual(raw, rows[0]["raw_text"])
                self.assertEqual("json", rows[1]["parse_status"])

    def test_malformed_response_does_not_bypass_integrity_guards(self):
        for index, (records, message) in enumerate([
            ([{"prompt_id": "unknown", "raw_text": "bad-json"}], "unknown prompt_id"),
            ([{"prompt_id": "p1", "raw_text": "bad-json"}, {"prompt_id": "p1", "raw_text": '{"claims":[]}'}], "duplicate prompt_id"),
            ([{"prompt_id": "p1", "case_id": "other", "raw_text": "bad-json"}], "case_id mismatch"),
            ([{"prompt_id": "p1", "route_b_unit_id": "other", "raw_text": "bad-json"}], "unit_id mismatch"),
            ([{"prompt_id": "p1", "route_b_unit_id": "", "raw_text": "bad-json"}], "unit_id mismatch"),
            ([{"prompt_id": "p1", "raw_text": '{}', "model_json": {"claims": []}}], "differs from raw_text"),
            ([{"prompt_id": "p1", "raw_text": 12}], "raw_text must be a string"),
            ([{"prompt_id": "p1", "model_json": {"claims": []}}], "requires raw_text"),
        ]):
            with self.subTest(message=message):
                partial = self.root / f"integrity-{index}.json"
                partial.write_text(json.dumps({"case_id": "case-x", "executor": "TEST", "model": "TEST", "responses": records}), encoding="utf-8")
                out = self.root / f"integrity-out-{index}.json"
                with self.assertRaisesRegex(ValueError, message):
                    merge(str(self.prompts), [str(partial)], str(out))
                self.assertFalse(out.exists())
                self.assertFalse(out.with_name(out.stem + "_manifest.json").exists())

    def test_format_failure_checkpoint_rejects_changed_frozen_prompts(self):
        partial = self.root / "frozen-partial.json"
        partial.write_text(json.dumps({"case_id": "case-x", "executor": "TEST", "model": "TEST", "responses": [{"prompt_id": "p1", "raw_text": "bad-json"}]}), encoding="utf-8")
        checkpoint = self.root / "frozen-checkpoint.json"
        merge(str(self.prompts), [str(partial)], str(checkpoint))
        changed_rows = [dict(row) for row in self.rows]
        changed_rows[0]["user"] = "changed prompt same IDs"
        self.prompts.write_text("".join(json.dumps(row) + "\n" for row in changed_rows), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "execution manifest input hash mismatch"):
            merge(str(self.prompts), [str(checkpoint)], str(self.root / "frozen-resumed.json"))

    def test_local_format_failure_does_not_allow_two_completed_attempts(self):
        inputs = []
        for index in range(2):
            partial = self.root / f"two-completed-{index}.json"
            partial.write_text(json.dumps({
                "case_id": "case-x", "executor": "TEST", "model": "TEST", "response_id": f"attempt-{index}",
                "responses": [{"prompt_id": "p1", "raw_text": '{"claims":[]}'}, {"prompt_id": "p2", "raw_text": "bad-json"}],
            }), encoding="utf-8")
            inputs.append(str(partial))
        with self.assertRaisesRegex(ValueError, "multiple completed responses for prompt_id"):
            merge(str(self.prompts), inputs, str(self.root / "two-completed-out.json"))

    def test_late_malformed_response_keeps_completed_result_in_both_input_orders(self):
        completed = self.root / "healthy-before-malformed.json"
        completed.write_text(json.dumps({"case_id": "case-x", "executor": "healthy", "model": "TEST", "responses": [{"prompt_id": "p1", "raw_text": '{"claims":[]}'}]}), encoding="utf-8")
        malformed = self.root / "late-malformed.json"
        malformed.write_text(json.dumps({"case_id": "case-x", "executor": "failed", "model": "TEST", "responses": [{"prompt_id": "p1", "raw_text": "bad-json"}]}), encoding="utf-8")
        for index, inputs in enumerate(([completed, malformed], [malformed, completed])):
            out = self.root / f"late-malformed-{index}.json"
            result = merge(str(self.prompts), [str(path) for path in inputs], str(out))
            row = json.loads(out.read_text())["candidates"][0]
            self.assertEqual(1, result["completed_count"])
            self.assertEqual("json", row["parse_status"])
            self.assertEqual("healthy", row["attempt_executor"])

    def test_claim_type_schema_failure_is_retryable_and_keeps_raw(self):
        for index, claim in enumerate(({}, {"claim_type": ""}, {"claim_type": "  "}, {"claim_type": None}, {"claim_type": 12})):
            with self.subTest(claim=claim):
                first = self.root / f"claim-type-first-{index}.json"
                raw = json.dumps({"claims": [claim]})
                record = {"prompt_id": "p2", "raw_text": raw}
                first.write_text(json.dumps({"case_id": "case-x", "executor": "TEST", "model": "TEST", "responses": [
                    {"prompt_id": "p1", "raw_text": '{"claims":[]}'}, record,
                    {"prompt_id": "p3", "raw_text": '{"claims":[]}'},
                ]}))
                checkpoint = self.root / f"claim-type-checkpoint-{index}.json"
                initial = merge(str(self.prompts), [str(first)], str(checkpoint))
                rows = json.loads(checkpoint.read_text())["candidates"]
                self.assertEqual(("partial", 2, 1), (initial["completion_status"], initial["completed_count"], initial["not_evaluated_count"]))
                self.assertEqual("not_evaluated", rows[1]["parse_status"])
                self.assertEqual("model_response_invalid_claim_schema", rows[1]["not_evaluated_reason"])
                self.assertEqual(raw, rows[1]["raw_text"])
                self.assertEqual(record, rows[1]["agent_response_raw"])
                self.assertEqual(claim, rows[1]["rejected_claims"][0]["raw_claim"])
                retry = self.root / f"claim-type-retry-{index}.json"
                retry.write_text(json.dumps({"case_id": "case-x", "executor": "TEST", "model": "TEST", "responses": [{"prompt_id": "p2", "raw_text": '{"claims":[]}'}]}))
                out = self.root / f"claim-type-resumed-{index}.json"
                resumed = merge(str(self.prompts), [str(checkpoint), str(retry)], str(out))
                self.assertEqual("complete", resumed["completion_status"])
                self.assertEqual(raw, resumed["retry_history"][0]["replaced_raw_text"])
                self.assertEqual(rows[0], json.loads(out.read_text())["candidates"][0])
                self.assertEqual(rows[2], json.loads(out.read_text())["candidates"][2])

    def test_mixed_claim_schema_failure_keeps_healthy_claim_over_whole_failure(self):
        good = {"claim_type": "function", "biological_function": "healthy function"}
        invalid = {"mechanism": "missing claim type"}
        mixed = self.root / "mixed-schema.json"
        mixed.write_text(json.dumps({"case_id": "case-x", "executor": "TEST", "model": "TEST", "responses": [{"prompt_id": "p1", "raw_text": json.dumps({"claims": [good, invalid]})}]}))
        failed = self.root / "mixed-timeout.json"
        failed.write_text(json.dumps({"case_id": "case-x", "executor": "TEST", "model": "TEST", "responses": [{"prompt_id": "p1", "status": "failed", "reason": "timeout"}]}))
        for index, paths in enumerate(([mixed, failed], [failed, mixed])):
            out = self.root / f"mixed-schema-out-{index}.json"
            result = merge(str(self.prompts), [str(path) for path in paths], str(out))
            row = json.loads(out.read_text())["candidates"][0]
            self.assertEqual(0, result["completed_count"])
            self.assertEqual(3, result["not_evaluated_count"])
            self.assertEqual("partial", row["response_completion_status"])
            self.assertEqual("not_evaluated", result["response_roster"][0]["completion_status"])
            self.assertEqual("model_response_invalid_claim_schema", result["response_roster"][0]["not_evaluated_reason"])
            self.assertEqual("json", row["parse_status"])
            self.assertEqual([good, invalid], row["model_json"]["claims"])
            self.assertEqual(invalid, row["rejected_claims"][0]["raw_claim"])

    def test_claim_schema_failure_does_not_bypass_integrity_or_completed_duplicate_guards(self):
        raw = '{"claims":[{}]}'
        for index, (record, message) in enumerate((
            ({"prompt_id": "unknown", "raw_text": raw}, "unknown prompt_id"),
            ({"prompt_id": "p1", "route_b_unit_id": "other", "raw_text": raw}, "unit_id mismatch"),
            ({"prompt_id": "p1", "case_id": "other", "raw_text": raw}, "case_id mismatch"),
        )):
            partial = self.root / f"schema-integrity-{index}.json"
            partial.write_text(json.dumps({"case_id": "case-x", "executor": "TEST", "model": "TEST", "responses": [record]}))
            with self.assertRaisesRegex(ValueError, message):
                merge(str(self.prompts), [str(partial)], str(self.root / f"schema-integrity-out-{index}.json"))
        paths = []
        for index in range(2):
            partial = self.root / f"schema-two-completed-{index}.json"
            partial.write_text(json.dumps({"case_id": "case-x", "executor": "TEST", "model": "TEST", "responses": [
                {"prompt_id": "p1", "raw_text": '{"claims":[]}', "response_completion_status": "partial"}, {"prompt_id": "p2", "raw_text": raw},
            ]}))
            paths.append(str(partial))
        with self.assertRaisesRegex(ValueError, "multiple completed responses"):
            merge(str(self.prompts), paths, str(self.root / "schema-two-completed-out.json"))


if __name__ == "__main__":
    unittest.main()
