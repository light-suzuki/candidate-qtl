"""Coverage observations do not promote source support or change membership."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from candidate_qtl.wet_candidate_handoff import (
    _initial_a_input_coverage, _input_coverage, _route_c_input_coverage,
)


class ACInputCoverageTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def table(self, name, rows):
        path = self.root / name
        pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
        return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}

    def producer(self, *, policy="initial_bibliography_v1", bibliography="resolved", panel="frozen"):
        hypotheses = self.table("hypotheses.tsv", [{"hypothesis_id": "H1",
            "candidate_admission_policy": policy, "bibliography_state": bibliography}])
        projection = self.table("panel.tsv", [{"hypothesis_id": "H1", "panel_status": panel,
            "sequence_sha256": "a" * 64, "projection_eligibility": "review_only"}])
        review = self.table("review.tsv", [{"hypothesis_id": "H1", "candidate_gene_id": "G1"}])
        return {"stage": "route_a_review_candidate_freeze", "output": review["path"],
            "inputs": {"route_a_hypotheses": hypotheses, "route_a_projection_panel": projection}}

    def test_initial_a_usable_work_is_not_support_or_execution_completion(self):
        result = _initial_a_input_coverage([self.producer()])
        self.assertEqual(1, result["usable_candidate_count"])
        self.assertEqual("unknown", result["execution_state"])
        state, detail = _input_coverage({"route_a_supported_rows": 0}, route_a_supplied=True,
            route_b_supplied=False, initial_a_coverage=result)
        self.assertEqual("execution_coverage_unknown", state)
        self.assertEqual(0, detail["route_a_supported_rows"])
        self.assertNotIn("degraded", detail)

    def test_a_identity_or_bibliography_unavailable_does_not_become_usable(self):
        for kwargs in ({"panel": "unavailable"}, {"bibliography": "source_absent"}):
            with self.subTest(kwargs=kwargs):
                result = _initial_a_input_coverage([self.producer(**kwargs)])
                state, _ = _input_coverage({"route_a_supported_rows": 0}, route_a_supplied=True,
                    route_b_supplied=False, initial_a_coverage=result)
                self.assertEqual(0, result["usable_candidate_count"])
                self.assertEqual("route_a_initial_candidate_evidence_unavailable", state)

    def test_bound_partial_status_remains_partial_despite_valid_review(self):
        producer = self.producer()
        producer["inputs"]["import_status"] = self.table("status.tsv", [{"status": "imported"}, {"status": "invalid"}])
        result = _initial_a_input_coverage([producer])
        self.assertEqual(1, result["usable_candidate_count"])
        state, _ = _input_coverage({"route_a_supported_rows": 0}, route_a_supplied=True,
            route_b_supplied=False, initial_a_coverage=result)
        self.assertEqual("route_a_initial_import_incomplete", state)

    def test_successful_status_rows_do_not_prove_assigned_roster_completion(self):
        producer = self.producer()
        producer["inputs"]["import_status"] = self.table("status.tsv", [{"status": "imported"}])
        result = _initial_a_input_coverage([producer])
        self.assertEqual("unknown", result["import_state"])
        self.assertEqual("unknown", result["execution_state"])

    def test_missing_policy_and_content_policy_retain_legacy_meaning(self):
        for policy in ("content_v1", ""):
            self.assertIsNone(_initial_a_input_coverage([self.producer(policy=policy)]))
        producer = self.producer()
        producer["inputs"]["route_a_hypotheses"] = self.table("no-policy.tsv", [{"hypothesis_id": "H1"}])
        self.assertIsNone(_initial_a_input_coverage([producer]))
        state, _ = _input_coverage({"route_a_supported_rows": 0}, route_a_supplied=True, route_b_supplied=False)
        self.assertEqual("route_a_declared_but_unsupported", state)

    def c_coverage(self, records, projections, accepted=()):
        frozen = self.table("frozen.tsv", records)
        projection = self.table("projection.tsv", projections)
        member = self.table("members.tsv", [{"record_id": "R1", "support_status": "ledger_only"}])
        manifest = self.root / "manifest.json"
        manifest.write_text(json.dumps({"output": Path(member["path"]).name,
            "inputs": {"frozen": frozen, "projection": projection}}))
        return _route_c_input_coverage(str(manifest), [{"route": "c", "eligible_record_ids": list(accepted)}])

    def test_verified_no_interval_match_is_evaluated_not_negative(self):
        result = self.c_coverage([{"record_id": "R1", "source_verification_state": "verified"}],
            [{"record_id": "R1", "projection_state": "no_interval_match"}])
        self.assertEqual("evaluated_supplied_evidence", result["state"])
        self.assertFalse(result["no_interval_match_is_biological_negative"])
        state, _ = _input_coverage({"route_c_verified_support_rows": 0}, route_a_supplied=False,
            route_b_supplied=False, route_c_supplied=True, route_c_coverage=result)
        self.assertEqual("single_route_by_design", state)

    def test_inside_projection_requires_valid_review_or_support(self):
        for accepted, expected in (((), 1), (("R1",), 0)):
            result = self.c_coverage([{"record_id": "R1", "source_verification_state": "verified"}],
                [{"record_id": "R1", "projection_state": "inside_current_qtl"}], accepted)
            self.assertEqual(expected, result["unresolved_record_count"])

    def test_pending_not_projected_and_failed_are_not_completed(self):
        result = self.c_coverage([
            {"record_id": "R1", "source_verification_state": "pending"},
            {"record_id": "R2", "source_verification_state": "verified"},
            {"record_id": "R3", "source_verification_state": "failed"}], [
            {"record_id": "R1", "projection_state": "no_interval_match"},
            {"record_id": "R2", "projection_state": "not_projected"},
            {"record_id": "R3", "projection_state": "no_interval_match"}])
        self.assertEqual(3, result["unresolved_record_count"])
        state, _ = _input_coverage({"route_c_verified_support_rows": 0}, route_a_supplied=False,
            route_b_supplied=False, route_c_supplied=True, route_c_coverage=result)
        self.assertEqual("route_c_partial_input_coverage", state)

    def test_tampered_ancestor_is_not_coverage_evidence(self):
        producer = self.producer()
        Path(producer["inputs"]["route_a_hypotheses"]["path"]).write_text("changed")
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            _initial_a_input_coverage([producer])

    def test_source_pending_cannot_be_cleared_by_accepted_record_id(self):
        result = self.c_coverage([{"record_id": "R1", "source_verification_state": "pending"}],
            [{"record_id": "R1", "projection_state": "inside_current_qtl"}], ("R1",))
        self.assertEqual({"source_pending_or_unverified": 1}, result["record_state_counts"])

    def test_missing_c_ancestry_stays_unknown(self):
        manifest = self.root / "manifest.json"
        manifest.write_text(json.dumps({"inputs": {}}))
        result = _route_c_input_coverage(str(manifest), [])
        state, _ = _input_coverage({"route_c_verified_support_rows": 0}, route_a_supplied=False,
            route_b_supplied=False, route_c_supplied=True, route_c_coverage=result)
        self.assertEqual("execution_coverage_unknown", state)
