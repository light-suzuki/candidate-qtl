"""Operational target and provenance boundaries for unverified Wet reviews."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from candidate_qtl.boundary_tier_review import partition_relevance_tiers
from candidate_qtl.relevance_tiers import T2, T3, T4
from tests._issue664_helpers import manifest_for


class ReviewRetentionBoundaryTests(unittest.TestCase):
    def test_protected_c_is_retained_without_displacing_ordinary_tiers(self):
        rows = [
            {"candidate_gene_id": "c-low", "relevance_tier": T4},
            *({"candidate_gene_id": f"ab-{index}", "relevance_tier": T2}
              for index in range(20)),
            {"candidate_gene_id": "ab-lower", "relevance_tier": T3},
        ]
        plan = partition_relevance_tiers(rows, target=20, upper=25,
                                         protected_gene_ids={"c-low"})
        self.assertIn("c-low", {row["candidate_gene_id"] for row in plan["frozen"]})
        self.assertEqual(20, len([row for row in plan["frozen"]
                                  if row["candidate_gene_id"].startswith("ab-")]))
        self.assertEqual([], plan["boundary"])

        same_tier = [
            *({"candidate_gene_id": f"c-{index}", "relevance_tier": T2} for index in range(2)),
            *({"candidate_gene_id": f"ab-{index}", "relevance_tier": T2} for index in range(18)),
        ]
        same_plan = partition_relevance_tiers(same_tier, target=20, upper=25,
                                              protected_gene_ids={"c-0", "c-1"})
        self.assertEqual(20, len(same_plan["frozen"]))

        c_only = [
            *({"candidate_gene_id": f"c-{index}", "relevance_tier": T2} for index in range(30)),
            *({"candidate_gene_id": f"dark-{index}", "relevance_tier": T4} for index in range(100)),
        ]
        c_plan = partition_relevance_tiers(c_only, target=20, upper=25,
                                           protected_gene_ids={f"c-{index}" for index in range(30)})
        self.assertEqual(30, len(c_plan["frozen"]))
        self.assertEqual([], c_plan["boundary"])

    def _inputs(self, root: Path, count: int) -> dict:
        atlas = root / "atlas.tsv"
        genes = [f"G{i:03}" for i in range(count)]
        pd.DataFrame({"gene_id": genes}).to_csv(atlas, sep="\t", index=False)
        supported = root / "supported.tsv"
        pd.DataFrame(columns=["candidate_gene_id", "route", "support_status"]).to_csv(
            supported, sep="\t", index=False)
        reviews = root / "reviews.tsv"
        pd.DataFrame({"candidate_gene_id": genes, "route": "b",
                      "support_status": "unverified_route_b_review", "rank_contribution": "none"}).to_csv(
            reviews, sep="\t", index=False)
        seal = root / "reviews_manifest.json"
        seal.write_text(json.dumps({"contract": "route-review-candidates-1.0.0",
            "stage": "route_b_review_candidate_freeze", "case_id": "case",
            "output": str(reviews), "artifact_sha256": hashlib.sha256(reviews.read_bytes()).hexdigest(),
            "inputs": {"atlas": {"path": str(atlas), "sha256": hashlib.sha256(atlas.read_bytes()).hexdigest()}}}))
        return dict(route_b_candidates_arg=str(supported),
            route_b_candidates_manifest_arg=str(manifest_for(supported, route="b", case_id="case")),
            review_candidates_arg=str(reviews), review_candidates_manifest_arg=str(seal),
            interval_atlas_arg=str(atlas), case_id="case", wet_budget=20)

    def test_target_preserves_all_reviews_above_twenty_without_padding_small_set(self):
        for count in (2, 31):
            with self.subTest(count=count), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                result = wet_candidate_handoff(str(root / "out"), **self._inputs(root, count))
                wet = pd.read_csv(result["wet_candidates"], sep="\t")
                self.assertEqual(count, len(wet))
                self.assertTrue(wet.membership_basis.eq("route_review_candidate").all())
                union = pd.read_csv(result["frozen_candidate_set"], sep="\t")
                self.assertEqual(0, len(union))

    def test_mutated_review_artifact_is_not_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = self._inputs(root, 2)
            review = Path(inputs["review_candidates_arg"])
            review.write_text(review.read_text().replace("G001", "G999"))
            with self.assertRaisesRegex(ValueError, "producer manifest"):
                wet_candidate_handoff(str(root / "out"), **inputs)

    def test_sealed_out_of_interval_id_is_not_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = self._inputs(root, 2)
            review = Path(inputs["review_candidates_arg"])
            review.write_text(review.read_text().replace("G001", "G999"))
            seal_path = Path(inputs["review_candidates_manifest_arg"])
            seal = json.loads(seal_path.read_text())
            seal["artifact_sha256"] = hashlib.sha256(review.read_bytes()).hexdigest()
            seal_path.write_text(json.dumps(seal))
            with self.assertRaisesRegex(ValueError, "outside"):
                wet_candidate_handoff(str(root / "out"), **inputs)
