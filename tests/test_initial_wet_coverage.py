"""Initial adoption and strict support have distinct count/coverage meanings."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

from tests import test_route_b_initial_reference_identity as identity_fixture
from candidate_qtl.wet_candidate_handoff import _input_coverage


class InitialWetCoverageTests(unittest.TestCase):
    def check_chain(self, root):
        _, _, frozen, wet = identity_fixture.RouteBInitialReferenceIdentityTests()._chain(root, [identity_fixture.ARABIDOPSIS])
        manifest = json.loads(Path(wet["manifest"]).read_text())
        self.assertEqual({"G1"}, set(pd.read_csv(wet["wet_candidates"], sep="\t").candidate_gene_id))
        self.assertEqual(0, frozen["candidate_gene_count"])
        self.assertEqual(0, manifest["candidate_count"])
        self.assertEqual(1, manifest["final_wet_candidate_count"])
        self.assertEqual({"supported_union": 0, "eligible_first_stage": 1, "active_first_stage": 1, "final_wet": 1,
                          "adopted": 1, "pending": 0}, manifest["candidate_count_populations"])
        self.assertEqual(manifest["candidate_count_populations"], manifest["report"]["candidate_count_populations"])
        self.assertIn('"supported_union": 0', Path(manifest["report"]["path"]).read_text())
        return manifest

    def test_completed_initial_b_review_is_usable_without_verified_support(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = self.check_chain(Path(tmp))
            self.assertEqual("wet_reviewable", manifest["handoff_state"])
            self.assertEqual("single_route_by_design", manifest["input_coverage_state"])
            self.assertEqual(1, manifest["input_coverage_detail"]["route_b_initial_candidate_count"])

    def test_partial_import_keeps_valid_candidate_and_explicit_degradation(self):
        original = identity_fixture.import_route_b_semantic_proposals

        def partial(*args, **kwargs):
            response = Path(args[2])
            payload = json.loads(response.read_text())
            payload["candidates"][0]["model_json"]["claims"].append({})
            response.write_text(json.dumps(payload))
            return original(*args, **kwargs)

        with tempfile.TemporaryDirectory() as tmp, patch.object(identity_fixture, "import_route_b_semantic_proposals", partial):
            manifest = self.check_chain(Path(tmp))
            self.assertEqual("wet_reviewable_degraded_coverage", manifest["handoff_state"])
            self.assertEqual("import_incomplete", manifest["input_coverage_detail"]["initial_b_import_and_source_observations"][0]["status"])

    def test_legacy_coverage_and_independent_failed_routes_remain_visible(self):
        union = {"route_a_supported_rows": 0, "route_b_verified_candidate_rows": 0}
        state, _ = _input_coverage(union, route_a_supplied=False, route_b_supplied=True)
        self.assertEqual("route_b_declared_but_unverified", state)
        state, _ = _input_coverage(union, route_a_supplied=True, route_b_supplied=True, initial_b_candidate_count=1)
        self.assertEqual("route_a_declared_but_unsupported", state)
