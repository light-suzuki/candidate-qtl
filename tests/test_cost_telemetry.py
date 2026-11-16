"""Tests for Issue #434: Cost telemetry aggregation, stage wrappers, and profile reports."""

from __future__ import annotations

import unittest

from candidate_qtl.budget_account import (
    BUDGET_PROFILES,
    BudgetAccount,
    stage_telemetry,
)
from candidate_qtl.budget_engine import (
    build_cost_profile_comparison_report,
    summarize_case_cost_telemetry,
)
from candidate_qtl.llm_usage import apply_llm_usage_to_account


class TestCostTelemetry(unittest.TestCase):
    """Verifies stage telemetry measurement, LLM token logging, and profile comparison reporting."""

    def test_stage_telemetry_captures_wall_cpu_and_tracemalloc_memory(self) -> None:
        account = BudgetAccount({"peak_memory_limit_mb": 512})

        with stage_telemetry(account, "python_stage_allocation") as tel:
            # Allocate memory to trigger tracemalloc
            data = [i for i in range(100000)]
            self.assertEqual(100000, len(data))

        self.assertGreater(tel["wall_time_seconds"], 0.0)
        self.assertGreaterEqual(tel["cpu_time_seconds"], 0.0)
        self.assertGreater(tel["peak_memory_mb"], 0.0)
        self.assertFalse(tel["is_heavy"])

        frozen = account.freeze()
        self.assertGreater(frozen["usage"]["peak_memory_mb"], 0.0)
        self.assertGreaterEqual(frozen["usage"]["cpu_time_seconds"], 0.0)
        self.assertEqual(0, frozen["usage"]["heavy_stage_invocation_count"])

    def test_stage_telemetry_heavy_stage_logs_reason(self) -> None:
        account = BudgetAccount()
        with stage_telemetry(
            account,
            "complex_homology_refinement",
            is_heavy=True,
            invocation_reason="candidate_count_exceeds_wet_budget",
            metadata={"focal_species": "Pisum sativum"},
        ) as tel:
            self.assertTrue(tel["is_heavy"])

        heavy = account.heavy_stages()
        self.assertEqual(1, len(heavy))
        self.assertEqual("complex_homology_refinement", heavy[0]["stage"])
        self.assertEqual("candidate_count_exceeds_wet_budget", heavy[0]["invocation_reason"])
        self.assertEqual("Pisum sativum", heavy[0]["metadata"]["focal_species"])
        self.assertEqual(1, account.freeze()["usage"]["heavy_stage_invocation_count"])

    def test_apply_llm_usage_to_account_consumes_provider_tokens(self) -> None:
        account = BudgetAccount({"llm_call_limit": 5, "token_limit": 10000})
        provider_manifest = {
            "call_count": 2,
            "input_tokens": 1200,
            "output_tokens": 350,
            "total_tokens": 1550,
        }

        res = apply_llm_usage_to_account(account, provider_manifest)
        self.assertEqual(2, res["call_count"])
        self.assertEqual(1550, res["total_tokens"])

        frozen = account.freeze()
        self.assertEqual(2, frozen["usage"]["llm_call_count"])
        self.assertEqual(1200, frozen["usage"]["llm_input_tokens"])
        self.assertEqual(350, frozen["usage"]["llm_output_tokens"])
        self.assertEqual(3, frozen["remaining"]["llm_call_limit"])
        self.assertEqual(8450, frozen["remaining"]["token_limit"])

    def test_build_cost_profile_comparison_report(self) -> None:
        # Case 1: light profile (deterministic CPU-only, 0 heavy stages, 0 LLM)
        acc_light = BudgetAccount(BUDGET_PROFILES["light"])
        acc_light.consume("sequence_search_budget", 5)
        acc_light.add_usage("cache_hit_count", 4)
        acc_light.add_usage("cache_miss_count", 1)
        acc_light.consume("wall_time_budget", 1.5)
        acc_light.observe_peak("peak_memory_budget", 45.0)

        sum_light = summarize_case_cost_telemetry(
            acc_light,
            initial_candidate_count=50,
            final_candidate_count=18,
            wet_budget=20,
            case_id="case_001",
            execution_profile="light",
        )

        # Case 1: heavy profile (with LLM and heavy refinement)
        acc_heavy = BudgetAccount(BUDGET_PROFILES["heavy"])
        acc_heavy.consume("sequence_search_budget", 25)
        acc_heavy.consume("llm_call_budget", 4)
        acc_heavy.consume("llm_input_tokens", 3000)
        acc_heavy.consume("llm_output_tokens", 800)
        acc_heavy.consume("wall_time_budget", 12.0)
        acc_heavy.observe_peak("peak_memory_budget", 180.0)
        acc_heavy.record_heavy_stage("llm_cross_species_synthesis", "candidate_count_exceeds_wet_budget")

        sum_heavy = summarize_case_cost_telemetry(
            acc_heavy,
            initial_candidate_count=50,
            final_candidate_count=8,
            wet_budget=20,
            case_id="case_001",
            execution_profile="heavy",
        )

        report = build_cost_profile_comparison_report([sum_light, sum_heavy], case_roster=["case_001"])
        self.assertEqual(["case_001"], report["roster"])
        self.assertIn("light", report["profiles"])
        self.assertIn("heavy", report["profiles"])

        light_prof = report["profiles"]["light"]
        self.assertEqual(1, light_prof["case_count"])
        self.assertEqual(1.0, light_prof["within_wet_budget_rate"])
        self.assertEqual(0.0, light_prof["heavy_stage_invocation_rate"])
        self.assertEqual(0, light_prof["total_llm_tokens"])

        heavy_prof = report["profiles"]["heavy"]
        self.assertEqual(1, heavy_prof["case_count"])
        self.assertEqual(1.0, heavy_prof["within_wet_budget_rate"])
        self.assertEqual(1.0, heavy_prof["heavy_stage_invocation_rate"])
        self.assertEqual(3800, heavy_prof["total_llm_tokens"])
        self.assertEqual(4, heavy_prof["total_llm_calls"])


if __name__ == "__main__":
    unittest.main()
