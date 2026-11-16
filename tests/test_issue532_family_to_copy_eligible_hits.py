"""Focused regression tests for Issue #532 long-form BLAST handoff."""

from __future__ import annotations

import unittest

import pandas as pd

from candidate_qtl.family_to_copy import _normalize_route_a_blast_hits


class FamilyToCopyEligibleHitsTests(unittest.TestCase):
    def test_normalizes_every_verified_eligible_target(self) -> None:
        frame = pd.DataFrame(
            [
                {
                    "panel_id": "Q70KR7",
                    "target_id": "Psat_NFR1",
                    "pident": "91.2",
                    "query_coverage": "0.94",
                    "target_coverage": "0.88",
                    "evalue": "1e-70",
                    "eligibility": "eligible",
                    "review_status": "verified",
                    "provider_status": "evaluable",
                    "ambiguity_status": "ambiguous_paralog",
                    "rank_contribution": "none",
                },
                {
                    "panel_id": "Q70KR7",
                    "target_id": "Psat_NFR1_like",
                    "pident": "87.4",
                    "query_coverage": "0.91",
                    "target_coverage": "0.83",
                    "evalue": "2e-60",
                    "eligibility": "eligible",
                    "review_status": "verified",
                    "provider_status": "evaluable",
                    "ambiguity_status": "ambiguous_paralog",
                    "rank_contribution": "none",
                },
            ]
        )

        normalized = _normalize_route_a_blast_hits(frame)

        self.assertEqual(["Psat_NFR1", "Psat_NFR1_like"], normalized["top_hit_target_id"].tolist())
        self.assertEqual(["91.2", "87.4"], normalized["top_hit_identity"].tolist())
        self.assertEqual(2, len(normalized))


if __name__ == "__main__":
    unittest.main()
