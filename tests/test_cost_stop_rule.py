"""Tests for Issue #434: Cheap-first execution and wet-budget stop rule contract."""

from __future__ import annotations

import unittest

from candidate_qtl.budget_account import BudgetAccount, stage_telemetry
from candidate_qtl.candidate_budget import budget_decision
from candidate_qtl.budget_engine import summarize_case_cost_telemetry


class TestCostStopRule(unittest.TestCase):
    """Verifies that cheap-first filtering halts expensive stages when within wet budget."""

    def test_cheap_stage_within_wet_budget_halts_heavy_stages(self) -> None:
        """When initial cheap deterministic narrowing yields <= wet_budget, no heavy stage is run."""
        account = BudgetAccount({"sequence_query_limit": 50, "llm_call_limit": 10})
        wet_budget = 20
        initial_candidates = 80

        # Step 1: Cheap deterministic stage
        with stage_telemetry(account, "cheap_deterministic_filter"):
            narrowed_candidates = 15  # <= wet_budget

        decision = budget_decision(narrowed_candidates, wet_budget=wet_budget, adaptive_round=0)
        self.assertEqual("within_wet_budget", decision["status"])
        self.assertEqual("within_wet_budget", decision["decision_reason"])
        self.assertFalse(decision["adaptive_refinement_required"])

        # Stop rule check: adaptive refinement is False, so expensive LLM/structure stage is skipped
        if decision["adaptive_refinement_required"]:
            with stage_telemetry(account, "heavy_llm_literature_synthesis", is_heavy=True,
                                 invocation_reason=decision["decision_reason"]):
                pass

        self.assertEqual(0, account.freeze()["usage"]["heavy_stage_invocation_count"])
        self.assertEqual([], account.heavy_stages())

        summary = summarize_case_cost_telemetry(
            account,
            initial_candidate_count=initial_candidates,
            final_candidate_count=narrowed_candidates,
            wet_budget=wet_budget,
            case_id="case_pea_test",
            execution_profile="standard",
            status="completed",
        )
        self.assertTrue(summary["candidate_counts"]["within_wet_budget"])
        self.assertEqual(0.8125, summary["candidate_counts"]["reduction_ratio"])
        self.assertEqual(0, summary["compute"]["heavy_stage_invocation_count"])
        self.assertEqual([], summary["invocation_reasons"])

    def test_cheap_stage_exceeds_wet_budget_invokes_heavy_stage_with_reason(self) -> None:
        """When candidates exceed wet budget, heavy refinement is conditionally invoked with reason."""
        account = BudgetAccount({"sequence_query_limit": 50, "llm_call_limit": 10})
        wet_budget = 20
        initial_candidates = 120

        # Step 1: Cheap deterministic stage
        narrowed_round0 = 45  # > wet_budget
        decision_round0 = budget_decision(narrowed_round0, wet_budget=wet_budget, adaptive_round=0, max_adaptive_rounds=1,
                                          budget_policy="strict")
        self.assertEqual("adaptive_refinement_required", decision_round0["status"])
        self.assertEqual("candidate_count_exceeds_wet_budget", decision_round0["decision_reason"])
        self.assertTrue(decision_round0["adaptive_refinement_required"])

        # Heavy refinement invoked because adaptive_refinement_required is True
        if decision_round0["adaptive_refinement_required"]:
            with stage_telemetry(
                account,
                "heavy_semantic_falsification",
                is_heavy=True,
                invocation_reason=decision_round0["decision_reason"],
            ):
                account.consume("llm_call_budget", 2)
                account.consume("llm_input_tokens", 1500)
                account.consume("llm_output_tokens", 300)
                narrowed_round1 = 12  # Now <= wet_budget

        decision_round1 = budget_decision(narrowed_round1, wet_budget=wet_budget, adaptive_round=1, max_adaptive_rounds=1,
                                          budget_policy="strict")
        self.assertEqual("within_wet_budget", decision_round1["status"])
        self.assertFalse(decision_round1["adaptive_refinement_required"])

        self.assertEqual(1, account.freeze()["usage"]["heavy_stage_invocation_count"])
        self.assertEqual(1, len(account.heavy_stages()))
        self.assertEqual("candidate_count_exceeds_wet_budget", account.heavy_stages()[0]["invocation_reason"])
        self.assertEqual("heavy_semantic_falsification", account.heavy_stages()[0]["stage"])

        summary = summarize_case_cost_telemetry(
            account,
            initial_candidate_count=initial_candidates,
            final_candidate_count=narrowed_round1,
            wet_budget=wet_budget,
            case_id="case_fas_heavy",
            execution_profile="standard",
            status="completed",
        )
        self.assertTrue(summary["candidate_counts"]["within_wet_budget"])
        self.assertEqual(0.9, summary["candidate_counts"]["reduction_ratio"])
        self.assertEqual(1, summary["compute"]["heavy_stage_invocation_count"])
        self.assertEqual(2, summary["compute"]["llm_call_count"])
        self.assertEqual(1800, summary["compute"]["llm_total_tokens"])
        self.assertEqual(["candidate_count_exceeds_wet_budget"], summary["invocation_reasons"])

    def test_exhausted_rounds_halts_further_heavy_refinement(self) -> None:
        """When max adaptive rounds are exhausted, heavy refinement stops with review required."""
        account = BudgetAccount({"sequence_query_limit": 50, "llm_call_limit": 10})
        wet_budget = 20

        # Attempt round 1 when max_adaptive_rounds=1
        decision = budget_decision(35, wet_budget=wet_budget, adaptive_round=1, max_adaptive_rounds=1,
                                   budget_policy="strict")
        self.assertEqual("review_required_over_budget", decision["status"])
        self.assertEqual("adaptive_round_budget_exhausted", decision["decision_reason"])
        self.assertFalse(decision["adaptive_refinement_required"])

        # No heavy refinement should be executed
        if decision["adaptive_refinement_required"]:
            with stage_telemetry(account, "heavy_stage", is_heavy=True, invocation_reason=decision["decision_reason"]):
                pass

        self.assertEqual(0, account.freeze()["usage"]["heavy_stage_invocation_count"])


if __name__ == "__main__":
    unittest.main()
