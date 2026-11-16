"""Truth-blind provisional Wet policy on neutral interval and screen records."""

from __future__ import annotations

import random
import unittest

from candidate_qtl.boundary_tier_review import partition_relevance_tiers
from candidate_qtl.provisional_wet import assess_provisional_wet
from candidate_qtl.relevance_tiers import T3, T4, UNRESOLVED


def roster(concrete: int = 2, identity: int = 30, unresolved: int = 3) -> list[dict]:
    return [
        {"candidate_gene_id": f"gene_{index:03d}", "relevance_tier": tier,
         "relevance_evidence_json": "[]", "membership_basis": "interval_review_coverage",
         "route_a_supported": "False", "route_b_supported": "False", "route_c_supported": "False"}
        for index, tier in enumerate([T3] * concrete + [T4] * identity + [UNRESOLVED] * unresolved)
    ]


def screen(rows: list[dict]) -> list[dict]:
    return [
        {"case_id": "neutral", "unit_id": f"unit_{index:03d}",
         "member_interval_genes": row["candidate_gene_id"],
         "screen_state": "screened" if index % 2 else "not_evaluated_material_gap",
         "answer_sha256": "a" * 64 if index % 2 else ""}
        for index, row in enumerate(rows)
    ]


def assess(rows: list[dict], *, protected: tuple[str, ...] = ()) -> dict:
    plan = partition_relevance_tiers(rows, protected_gene_ids=protected)
    return assess_provisional_wet(rows, plan, screen(rows), case_id="neutral",
                                  protected_gene_ids=protected)


class ProvisionalWetTests(unittest.TestCase):
    def test_full_screen_coverage_survives_exact_agreed_deferrals(self) -> None:
        original = roster(concrete=30)
        active = original[:2] + original[30:]
        deferred = [row["candidate_gene_id"] for row in original[2:30]]
        result = assess_provisional_wet(active, partition_relevance_tiers(active), screen(original),
                                        case_id="neutral", coverage_inventory=original,
                                        agreed_deferred_gene_ids=deferred)
        self.assertTrue(result["eligible"])
        self.assertEqual(2, len(result["adopted_gene_ids"]))
        self.assertEqual(len(original) - 2, len(result["unresolved_gene_ids"]))
        for bad in (deferred[:-1], deferred + [deferred[0]]):
            with self.assertRaisesRegex(ValueError, "exactly account"):
                assess_provisional_wet(active, partition_relevance_tiers(active), screen(original),
                                       case_id="neutral", coverage_inventory=original,
                                       agreed_deferred_gene_ids=bad)
        changed = [dict(row) for row in active]
        changed[0]["membership_basis"] = "changed"
        with self.assertRaisesRegex(ValueError, "preserve original rows"):
            assess_provisional_wet(changed, partition_relevance_tiers(changed), screen(original),
                                   case_id="neutral", coverage_inventory=original,
                                   agreed_deferred_gene_ids=deferred)

    def test_sparse_concrete_and_large_identity_tier_retains_all_concrete(self) -> None:
        rows = roster()
        result = assess(rows)
        self.assertTrue(result["eligible"])
        self.assertEqual(["gene_000", "gene_001"], result["adopted_gene_ids"])
        self.assertEqual(33, len(result["unresolved_gene_ids"]))
        self.assertEqual(18, result["material_gap_unit_count"])
        shuffled = rows.copy()
        random.Random(7).shuffle(shuffled)
        shuffled_result = assess_provisional_wet(shuffled, partition_relevance_tiers(shuffled),
                                                  screen(rows), case_id="neutral")
        self.assertEqual(result, shuffled_result)

    def test_informative_boundary_keeps_whole_tier_review(self) -> None:
        rows = roster()
        rows[3]["relevance_evidence_json"] = '[{"claim_id":"independent_claim"}]'
        result = assess(rows)
        self.assertFalse(result["eligible"])
        self.assertEqual("identity_tier_has_comparison_evidence", result["reason"])

    def test_no_concrete_or_uncovered_gene_cannot_be_adopted(self) -> None:
        rows = roster(concrete=0)
        self.assertEqual("no_concrete_hypothesis", assess(rows)["reason"])
        rows = roster()
        plan = partition_relevance_tiers(rows)
        result = assess_provisional_wet(rows, plan, screen(rows)[:-1], case_id="neutral")
        self.assertEqual("screen_gene_roster_incomplete", result["reason"])

    def test_protected_c_and_over_upper_are_never_trimmed(self) -> None:
        rows = roster(concrete=2, identity=35)
        protected = ("gene_003", "gene_004")
        result = assess(rows, protected=protected)
        self.assertTrue(result["eligible"])
        self.assertTrue(set(protected) <= set(result["adopted_gene_ids"]))
        over = tuple(f"gene_{index:03d}" for index in range(2, 29))
        result = assess(rows, protected=over)
        self.assertFalse(result["eligible"])
        self.assertEqual("concrete_or_protected_over_upper", result["reason"])

    def test_material_bearing_unanswered_unit_blocks_provisional_completion(self) -> None:
        rows = roster()
        ledger = screen(rows)
        ledger[3]["screen_state"] = "not_evaluated"
        result = assess_provisional_wet(rows, partition_relevance_tiers(rows), ledger,
                                        case_id="neutral")
        self.assertEqual("screen_not_complete", result["reason"])


if __name__ == "__main__":
    unittest.main()
