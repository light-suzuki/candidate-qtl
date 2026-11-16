"""Original B replies can be checked without merging or rewriting evidence."""
import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

from scripts.agentic_route_b import main, merge, preflight


class RouteBPreflightTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.prompts = self.root / "prompts.jsonl"
        self.prompts.write_text(json.dumps({"prompt_id": "p1", "route_b_unit_id": "u1", "case_id": "case"}) + "\n", encoding="utf-8")
        self.response = self.root / "response.json"

    def write(self, claims, **record_changes):
        record = {"prompt_id": "p1", "route_b_unit_id": "u1", "raw_text": json.dumps({"claims": claims}), **record_changes}
        self.response.write_text(json.dumps({"case_id": "case", "executor": "agent", "model": "test", "responses": [record]}), encoding="utf-8")

    def paper_claim(self, **changes):
        return {"claim_type": "reference_function", "source_ids": ["PMID:123"], "paper_reports_json": [{
            "source_id": "PMID:123", "reported_gene": "REF1", "reported_species": "Plant alpha",
            "initial_tier": "T4", "initial_tier_reason": "Reported function", **changes}]}

    def test_enriched_report_accepted_without_writing_or_rewriting(self):
        self.write([self.paper_claim(research_annotations={"read_scope": "abstract", "counterevidence": ["Unresolved transfer"]})])
        before = {p.name: p.read_bytes() for p in self.root.iterdir()}
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            status = main(["preflight", "--prompts", str(self.prompts), "--response", str(self.response)])
        self.assertEqual(0, status)
        report = json.loads(output.getvalue())
        self.assertEqual("accepted", report["status"])
        self.assertTrue(report["responses"][0]["completed"])
        self.assertEqual(before, {p.name: p.read_bytes() for p in self.root.iterdir()})

    def test_invalid_papers_are_claim_errors_and_cannot_hide_behind_function_claim(self):
        for changes in ({"reported_gene": ["REF1"]}, {"invented_scientific_fact": "yes"}):
            with self.subTest(changes=changes):
                self.write([self.paper_claim(**changes)])
                report = preflight(str(self.prompts), [str(self.response)])
                self.assertEqual("rejected", report["status"])
                self.assertFalse(report["responses"][0]["completed"])
                self.assertEqual(1, report["errors"][0]["claim_index"])
                self.assertEqual("$.responses[0].raw_text.claims[0]", report["errors"][0]["path"])
                output = io.StringIO()
                with contextlib.redirect_stdout(output):
                    self.assertEqual(1, main(["preflight", "--prompts", str(self.prompts), "--response", str(self.response)]))

    def test_raw_and_parsed_disagreement_is_rejected(self):
        self.write([], model_json={"claims": [{"claim_type": "other"}]})
        report = preflight(str(self.prompts), [str(self.response)])
        self.assertEqual("rejected", report["status"])
        self.assertIn("differs from raw_text", report["errors"][0]["reason"])

    def test_parse_and_roster_errors_are_machine_readable(self):
        self.write([], raw_text="broken")
        report = preflight(str(self.prompts), [str(self.response)])
        self.assertEqual("$.responses[0]", report["errors"][0]["path"])
        self.write([], route_b_unit_id="wrong")
        report = preflight(str(self.prompts), [str(self.response)])
        self.assertEqual("$.responses[0].route_b_unit_id", report["errors"][0]["path"])

    def test_multiple_completed_files_match_merge_rejection_even_when_identical(self):
        self.write([])
        other = self.root / "second.json"
        other.write_bytes(self.response.read_bytes())
        before = {p.name: p.read_bytes() for p in self.root.iterdir()}
        for paths in ([self.response, other], [other, self.response], [self.response, self.response]):
            with self.subTest(paths=paths):
                report = preflight(str(self.prompts), list(map(str, paths)))
                self.assertEqual("rejected", report["status"])
                self.assertIn("multiple completed responses", report["errors"][0]["reason"])
                self.assertEqual(before, {p.name: p.read_bytes() for p in self.root.iterdir()})
                with self.assertRaisesRegex(ValueError, "multiple completed responses"):
                    merge(str(self.prompts), list(map(str, paths)), str(self.root / "merged.json"))
        self.assertEqual(before, {p.name: p.read_bytes() for p in self.root.iterdir()})

    def test_failed_then_completed_retry_matches_merge_and_preserves_inputs(self):
        self.write([], status="failed", reason="provider unavailable")
        failed = self.root / "failed.json"
        failed.write_bytes(self.response.read_bytes())
        self.write([])
        before = {p.name: p.read_bytes() for p in self.root.iterdir()}
        for index, paths in enumerate(([failed, self.response], [self.response, failed])):
            with self.subTest(paths=paths):
                report = preflight(str(self.prompts), list(map(str, paths)))
                self.assertEqual("accepted", report["status"])
                self.assertEqual(1, sum(row["completed"] for row in report["responses"]))
                self.assertEqual(before, {name: (self.root / name).read_bytes() for name in before})
                result = merge(str(self.prompts), list(map(str, paths)), str(self.root / f"merged-{index}.json"))
                self.assertEqual(1, result["completed_count"])
        self.assertEqual(before, {name: (self.root / name).read_bytes() for name in before})


if __name__ == "__main__":
    unittest.main()
