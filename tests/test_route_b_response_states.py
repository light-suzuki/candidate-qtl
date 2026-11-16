"""Execution explanations come from raw replies, not a reason keyword guess."""
import json
import tempfile
import unittest
from pathlib import Path

from candidate_qtl.route_b_response_validation import route_b_execution_breakdown, route_b_response_state
from scripts.agentic_route_b import merge


class RouteBResponseStatesTests(unittest.TestCase):
    def test_mixed_replies_and_retry_preserve_healthy_claims(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            prompts = root / "prompts.jsonl"
            prompts.write_text("".join(json.dumps({"prompt_id": f"p{i}", "route_b_unit_id": f"u{i}", "case_id": "case"}) + "\n"
                                       for i in range(1, 9)), encoding="utf-8")
            good = {"claim_type": "reference_function", "biological_function": "specific function"}
            first = root / "first.json"
            first.write_text(json.dumps({"case_id": "case", "executor": "fixture", "model": "fixture", "responses": [
                {"prompt_id": "p1", "raw_text": json.dumps({"claims": [good]})},
                {"prompt_id": "p2", "raw_text": '{"claims":[]}'},
                {"prompt_id": "p3", "status": "not_evaluated", "reason": "no_partial_response"},
                {"prompt_id": "p4", "status": "failed", "reason": "unavailable"},
                {"prompt_id": "p5", "raw_text": "broken"},
                {"prompt_id": "p6", "raw_text": json.dumps({"claims": [good, {"claim_type": ""}]})},
                {"prompt_id": "p8", "raw_text": '{"claims":[false]}'},
            ]}), encoding="utf-8")
            original = first.read_bytes()
            out = root / "merged.json"
            initial = merge(str(prompts), [str(first)], str(out))
            self.assertEqual({"completed_claims": 1, "completed_no_claim": 1, "declared_not_evaluated": 1,
                              "declared_failed": 1, "invalid_response": 2, "partial_claims": 1,
                              "missing_response": 1, "unclassified_not_evaluated": 0}, initial["response_state_counts"])
            self.assertEqual((2, 6), (initial["completed_count"], initial["not_evaluated_count"]))
            rows = json.loads(out.read_text())["candidates"]
            self.assertEqual(initial["response_state_counts"], route_b_execution_breakdown(rows, manifest=initial))
            self.assertEqual(original, first.read_bytes())
            self.assertEqual(1, len(rows[5]["rejected_claims"]))
            retry = root / "retry.json"
            retry.write_text(json.dumps({"case_id": "case", "executor": "retry", "model": "fixture", "responses": [
                {"prompt_id": f"p{i}", "raw_text": '{"claims":[]}'} for i in range(3, 9)]}), encoding="utf-8")
            resumed = root / "resumed.json"
            final = merge(str(prompts), [str(out), str(retry)], str(resumed))
            self.assertEqual("complete", final["completion_status"])
            self.assertEqual((1, 7), (final["response_state_counts"]["completed_claims"], final["response_state_counts"]["completed_no_claim"]))
            self.assertEqual(8, sum(final["response_state_counts"].values()))
            self.assertEqual(rows[:2], json.loads(resumed.read_text())["candidates"][:2])

    def test_legacy_unknown_and_declared_status_do_not_guess_a_biological_reason(self):
        self.assertEqual("unclassified_not_evaluated", route_b_response_state({
            "parse_status": "not_evaluated", "not_evaluated_reason": "provider failed or no material"}))
        self.assertEqual("declared_not_evaluated", route_b_response_state({
            "parse_status": "not_evaluated", "agent_response_raw": {"status": "not_evaluated"},
            "not_evaluated_reason": "provider failed"}))

    def test_raw_partial_cannot_be_promoted_by_completion_metadata(self):
        payload = {"claims": [{"claim_type": "reference_function"}, False]}
        self.assertEqual("partial_claims", route_b_response_state({
            "parse_status": "json", "raw_text": json.dumps(payload), "model_json": payload,
            "response_completion_status": "completed"}))

    def test_declared_breakdown_cannot_forge_counts_or_use_boolean_integers(self):
        rows = [{"parse_status": "json", "raw_text": '{"claims":[]}', "model_json": {"claims": []}}]
        counts = route_b_execution_breakdown(rows, manifest={})
        for forged in ({**counts, "completed_no_claim": 0}, {**counts, "completed_no_claim": True}, None, {}):
            with self.subTest(counts=forged), self.assertRaisesRegex(ValueError, "state counts"):
                route_b_execution_breakdown(rows, manifest={"response_state_counts": forged})
