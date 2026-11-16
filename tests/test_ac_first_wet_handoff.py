"""A/C-first Wet handoff: no Route B, no annotation fill, Route-A rethink signal."""
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.core import sha256
from candidate_qtl.route_c_direct_members import build_route_c_direct_members
from candidate_qtl.wet_candidate_handoff import (
    DEFAULT_MIN_ROUTE_CANDIDATES,
    route_a_rethink_signal,
    wet_candidate_handoff,
)
from tests.test_issue722_c_review_retention import Fixture


class RouteARethinkSignalTests(unittest.TestCase):
    def test_empty_wet_requires_fresh_route_a(self):
        signal = route_a_rethink_signal(0, min_route_candidates=3)
        self.assertEqual("required", signal["state"])
        self.assertEqual("rerun_route_a_fresh_generation", signal["action"])
        self.assertEqual("none", signal["membership_effect"])

    def test_small_wet_recommends_and_enough_does_not(self):
        self.assertEqual("recommended", route_a_rethink_signal(2, min_route_candidates=3)["state"])
        enough = route_a_rethink_signal(3, min_route_candidates=3)
        self.assertEqual("not_needed", enough["state"])
        self.assertIsNone(enough["action"])

    def test_invalid_threshold_is_rejected(self):
        with self.assertRaises(ValueError):
            route_a_rethink_signal(5, min_route_candidates=0)


class ACFirstHandoffTests(unittest.TestCase):
    def _a_and_c_without_b(self, root, *, min_route_candidates=DEFAULT_MIN_ROUTE_CANDIDATES):
        fixture = Fixture()
        rows = [{"gene_id": f"G{i:04d}", "orthology_reference_gene": ""} for i in range(999)]
        rows.append({"gene_id": "Psat.A", "orthology_reference_gene": ""})
        projected, frozen, freeze, atlas = fixture._project(root, rows, None, include_binding=True)
        members = build_route_c_direct_members(
            str(frozen), str(freeze), projected["projection"], projected["manifest"],
            str(atlas), str(root / "members"), case_id="case1")
        review_path = root / "route_a_review.tsv"
        pd.DataFrame([{
            "candidate_gene_id": "G0007", "route": "a", "support_status": "review_required",
            "rank_contribution": "none", "claim_id": "fixture-a",
            "mechanism": "fixture stem development mechanism",
            "phenotype_relation": "fixture altered stem growth",
        }]).to_csv(review_path, sep="\t", index=False)
        review_manifest = root / "route_a_review_manifest.json"
        review_manifest.write_text(json.dumps({
            "contract": "route-review-candidates-1.0.0", "stage": "review_candidate_freeze",
            "case_id": "case1", "output": str(review_path), "artifact_sha256": sha256(review_path),
            "inputs": {},
        }), encoding="utf-8")
        result = wet_candidate_handoff(
            str(root / "wet"), case_id="case1",
            route_c_members_arg=members["members"], route_c_manifest_arg=members["manifest"],
            review_candidates_arg=[str(review_path)],
            review_candidates_manifest_arg=[str(review_manifest)],
            interval_atlas_arg=str(atlas), min_route_candidates=min_route_candidates)
        wet = pd.read_csv(result["wet_candidates"], sep="\t", dtype=str).fillna("")
        manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
        return wet, result, manifest

    def test_wet_holds_only_route_candidates_from_a_thousand_gene_interval(self):
        with tempfile.TemporaryDirectory() as temp:
            wet, result, manifest = self._a_and_c_without_b(Path(temp))
            self.assertEqual({"Psat.A", "G0007"}, set(wet["candidate_gene_id"]))
            self.assertNotIn("interval_review_coverage", set(wet["membership_basis"]))
            self.assertFalse(manifest["input_coverage_detail"]["route_b_supplied"])
            self.assertEqual("recommended", result["route_a_rethink"]["state"])
            self.assertEqual(result["route_a_rethink"], manifest["route_a_rethink"])

    def test_threshold_is_operational_and_never_changes_membership(self):
        with tempfile.TemporaryDirectory() as temp:
            wet, result, _ = self._a_and_c_without_b(Path(temp), min_route_candidates=2)
            self.assertEqual({"Psat.A", "G0007"}, set(wet["candidate_gene_id"]))
            self.assertEqual("not_needed", result["route_a_rethink"]["state"])

    def test_invalid_threshold_fails_before_writing_output(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaises(ValueError):
                wet_candidate_handoff(str(Path(temp) / "wet"), case_id="case1", min_route_candidates=0)
            self.assertFalse((Path(temp) / "wet").exists())


if __name__ == "__main__":
    unittest.main()
