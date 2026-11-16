from __future__ import annotations

import unittest

from candidate_qtl.budget_account import (
    BUDGET_PROFILES,
    CONCEPTUAL_TO_LEDGER,
    LEDGER_LIMIT_FIELDS,
    LEDGER_USAGE_FIELDS,
    BudgetAccount,
)


class BudgetAccountTests(unittest.TestCase):
    def test_declare_consume_remaining_exhausted_basics(self) -> None:
        account = BudgetAccount()
        account.declare({"llm_call_limit": 2})
        self.assertEqual(2, account.remaining("llm_call_budget"))
        self.assertTrue(account.consume("llm_call_budget", 1))
        self.assertEqual(1, account.remaining("llm_call_limit"))
        self.assertFalse(account.exhausted("llm_call_budget"))
        self.assertTrue(account.consume("llm_call_limit", 1))
        self.assertEqual(0, account.remaining("llm_call_budget"))
        self.assertTrue(account.exhausted("llm_call_budget"))
        self.assertFalse(account.consume("llm_call_budget", 1))
        self.assertEqual(0, account.remaining("llm_call_budget"))
        self.assertEqual(2, account.freeze()["usage"]["llm_call_count"])

    def test_unlimited_when_no_limits_declared(self) -> None:
        account = BudgetAccount()
        self.assertIsNone(account.remaining("llm_call_budget"))
        self.assertIsNone(account.remaining("token_limit"))
        self.assertTrue(account.consume("llm_call_budget", 1000))
        self.assertFalse(account.exhausted("llm_call_budget"))
        self.assertEqual([], account.exhausted_budgets())
        self.assertIsNone(account.freeze()["limits"]["sequence_query_limit"])

    def test_consume_false_does_not_decrement(self) -> None:
        account = BudgetAccount({"llm_call_limit": 1})
        self.assertFalse(account.consume("llm_call_budget", 2))
        self.assertEqual(1, account.remaining("llm_call_budget"))
        self.assertEqual(0, account.freeze()["usage"]["llm_call_count"])
        self.assertTrue(account.exhausted("llm_call_budget"))

    def test_unknown_kind_raises(self) -> None:
        account = BudgetAccount()
        with self.assertRaises(ValueError):
            account.consume("bogus_budget", 1)
        with self.assertRaises(ValueError):
            account.remaining("bogus_budget")
        with self.assertRaises(ValueError):
            account.exhausted("bogus_budget")

    def test_double_declare_raises(self) -> None:
        account = BudgetAccount({"llm_call_limit": 1})
        with self.assertRaises(ValueError):
            account.declare({"token_limit": 100})
        second = BudgetAccount()
        second.declare({"token_limit": 100})
        with self.assertRaises(ValueError):
            second.declare({"token_limit": 100})

    def test_negative_amount_raises(self) -> None:
        account = BudgetAccount({"llm_call_limit": 1})
        with self.assertRaises(ValueError):
            account.consume("llm_call_budget", -1)
        with self.assertRaises(ValueError):
            account.consume("llm_call_budget", 0.5)
        with self.assertRaises(ValueError):
            account.add_usage("searched_protein_count", -1)

    def test_token_limit_shared_between_input_and_output(self) -> None:
        account = BudgetAccount({"token_limit": 100})
        self.assertTrue(account.consume("llm_input_tokens", 60))
        self.assertEqual(40, account.remaining("token_budget"))
        self.assertFalse(account.consume("llm_output_tokens", 50))
        self.assertEqual(0, account.freeze()["usage"]["llm_output_tokens"])
        self.assertTrue(account.consume("llm_output_tokens", 40))
        self.assertEqual(0, account.remaining("token_limit"))
        self.assertTrue(account.exhausted("token_budget"))
        self.assertEqual(60, account.freeze()["usage"]["llm_input_tokens"])

    def test_sequence_search_budget_maps_to_sequence_query_count(self) -> None:
        account = BudgetAccount({"sequence_query_limit": 2})
        self.assertTrue(account.consume("sequence_search_budget", 1))
        self.assertTrue(account.consume("sequence_search_budget", 1))
        self.assertEqual(2, account.freeze()["usage"]["sequence_query_count"])
        self.assertFalse(account.consume("sequence_search_budget", 1))
        self.assertEqual(2, account.freeze()["usage"]["sequence_query_count"])

    def test_wall_time_budget_mapping(self) -> None:
        account = BudgetAccount({"wall_time_limit_seconds": 60})
        self.assertTrue(account.consume("wall_time_budget", 30.5))
        usage = account.freeze()["usage"]
        self.assertEqual(30.5, usage["wall_time_seconds"])
        self.assertEqual(29.5, account.remaining("wall_time_budget"))

    def test_consume_rejects_peak_memory_kind(self) -> None:
        account = BudgetAccount({"peak_memory_limit_mb": 1024})
        for kind in ("peak_memory_budget", "peak_memory_limit_mb", "peak_memory_mb"):
            with self.assertRaises(ValueError):
                account.consume(kind, 1)
        self.assertEqual(0, account.freeze()["usage"]["peak_memory_mb"])

    def test_observe_peak_records_max_not_sum_and_marks_exhausted_on_overrun(self) -> None:
        account = BudgetAccount({"peak_memory_limit_mb": 1024})
        self.assertTrue(account.observe_peak("peak_memory_budget", 512))
        self.assertFalse(account.exhausted("peak_memory_budget"))
        self.assertTrue(account.observe_peak("peak_memory_budget", 300))
        # two observations record the max (512), never the sum (812)
        self.assertEqual(512, account.freeze()["usage"]["peak_memory_mb"])
        self.assertFalse(account.exhausted("peak_memory_budget"))
        self.assertFalse(account.observe_peak("peak_memory_budget", 2000))
        # the ledger still records the true observed maximum, even over limit
        self.assertEqual(2000, account.freeze()["usage"]["peak_memory_mb"])
        self.assertTrue(account.exhausted("peak_memory_budget"))

    def test_observe_peak_remaining_never_negative(self) -> None:
        account = BudgetAccount({"peak_memory_limit_mb": 100})
        self.assertFalse(account.observe_peak("peak_memory_budget", 500))
        self.assertEqual(0, account.remaining("peak_memory_budget"))
        self.assertEqual(0, account.freeze()["remaining"]["peak_memory_limit_mb"])

    def test_observe_peak_rejects_non_finite_and_leaves_usage_unchanged(self) -> None:
        account = BudgetAccount({"peak_memory_limit_mb": 1024})
        with self.assertRaises(ValueError):
            account.observe_peak("peak_memory_budget", float("nan"))
        with self.assertRaises(ValueError):
            account.observe_peak("peak_memory_budget", float("inf"))
        with self.assertRaises(ValueError):
            account.observe_peak("peak_memory_budget", -1)
        self.assertEqual(0, account.freeze()["usage"]["peak_memory_mb"])

    def test_observe_peak_under_limit_leaves_exhausted_false(self) -> None:
        account = BudgetAccount({"peak_memory_limit_mb": 1024})
        self.assertTrue(account.observe_peak("peak_memory_budget", 1000))
        self.assertFalse(account.exhausted("peak_memory_budget"))
        self.assertEqual(24, account.remaining("peak_memory_budget"))

    def test_freeze_exact_keys_and_none_for_unlimited(self) -> None:
        account = BudgetAccount({"llm_call_limit": 3})
        frozen = account.freeze()
        self.assertEqual({"limits", "usage", "remaining", "exhausted"}, set(frozen))
        self.assertEqual(LEDGER_LIMIT_FIELDS, tuple(frozen["limits"]))
        self.assertEqual(LEDGER_USAGE_FIELDS, tuple(frozen["usage"]))
        self.assertEqual(LEDGER_LIMIT_FIELDS, tuple(frozen["remaining"]))
        self.assertEqual(LEDGER_LIMIT_FIELDS, tuple(frozen["exhausted"]))
        self.assertEqual(3, frozen["limits"]["llm_call_limit"])
        self.assertIsNone(frozen["limits"]["token_limit"])
        self.assertIsNone(frozen["remaining"]["token_limit"])
        self.assertFalse(frozen["exhausted"]["token_limit"])

    def test_budget_profiles_exact_values(self) -> None:
        self.assertEqual({"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 0,
                          "wall_time_limit_seconds": 0, "peak_memory_limit_mb": 0}, BUDGET_PROFILES["zero"])
        self.assertEqual({"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 50,
                          "wall_time_limit_seconds": 3600, "peak_memory_limit_mb": 8192}, BUDGET_PROFILES["student_cpu"])
        self.assertEqual({"llm_call_limit": None, "token_limit": None, "sequence_query_limit": None,
                          "wall_time_limit_seconds": None, "peak_memory_limit_mb": None}, BUDGET_PROFILES["exhaustive"])
        self.assertEqual({"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 20,
                          "wall_time_limit_seconds": 1800, "peak_memory_limit_mb": 4096}, BUDGET_PROFILES["light"])
        self.assertEqual(BUDGET_PROFILES["student_cpu"], BUDGET_PROFILES["standard"])
        self.assertEqual({"llm_call_limit": 50, "token_limit": 150000, "sequence_query_limit": 200,
                          "wall_time_limit_seconds": 7200, "peak_memory_limit_mb": 16384}, BUDGET_PROFILES["heavy"])
        self.assertEqual(BUDGET_PROFILES["exhaustive"], BUDGET_PROFILES["diagnostic"])

    def test_conceptual_to_ledger_mapping(self) -> None:
        self.assertEqual({"llm_call_budget": "llm_call_limit", "token_budget": "token_limit",
                          "sequence_search_budget": "sequence_query_limit",
                          "wall_time_budget": "wall_time_limit_seconds",
                          "peak_memory_budget": "peak_memory_limit_mb"}, CONCEPTUAL_TO_LEDGER)

    def test_declare_accepts_conceptual_keys_and_rejects_unknown(self) -> None:
        account = BudgetAccount({"sequence_search_budget": 3})
        self.assertEqual(3, account.remaining("sequence_query_limit"))
        with self.assertRaises(ValueError):
            BudgetAccount({"unknown_limit": 1})

    def test_declare_does_not_mutate_shared_profiles(self) -> None:
        account = BudgetAccount(dict(BUDGET_PROFILES["student_cpu"]))
        self.assertTrue(account.consume("sequence_search_budget", 1))
        self.assertEqual(50, BUDGET_PROFILES["student_cpu"]["sequence_query_limit"])
        self.assertEqual(1, account.freeze()["usage"]["sequence_query_count"])

    def test_zero_profile_blocks_everything(self) -> None:
        account = BudgetAccount(dict(BUDGET_PROFILES["zero"]))
        self.assertFalse(account.consume("llm_call_budget", 1))
        self.assertFalse(account.consume("token_budget", 1))
        self.assertFalse(account.consume("sequence_search_budget", 1))
        self.assertFalse(account.consume("wall_time_budget", 1))
        with self.assertRaises(ValueError):
            account.consume("peak_memory_budget", 1)
        self.assertFalse(account.observe_peak("peak_memory_budget", 1))
        self.assertEqual(5, len(account.exhausted_budgets()))

    def test_add_usage_rejects_non_integer_counts(self) -> None:
        account = BudgetAccount()
        with self.assertRaises(ValueError):
            account.add_usage("searched_protein_count", 1.5)
        with self.assertRaises(ValueError):
            account.add_usage("cache_hit_count", 0.5)
        with self.assertRaises(ValueError):
            account.add_usage("external_retrieval_count", 0.5)
        with self.assertRaises(ValueError):
            account.add_usage("heavy_stage_invocation_count", 0.5)
        account.add_usage("searched_protein_count", 3)
        self.assertEqual(3, account.freeze()["usage"]["searched_protein_count"])
        account.add_usage("sequence_search_cpu_seconds", 2.5)
        self.assertEqual(2.5, account.freeze()["usage"]["sequence_search_cpu_seconds"])
        account.add_usage("external_retrieval_count", 2)
        self.assertEqual(2, account.freeze()["usage"]["external_retrieval_count"])
        account.add_usage("cpu_time_seconds", 1.25)
        self.assertEqual(1.25, account.freeze()["usage"]["cpu_time_seconds"])

    def test_add_usage_rejects_budgeted_kinds_and_unknown_fields(self) -> None:
        account = BudgetAccount()
        for field in ("llm_call_count", "llm_input_tokens", "llm_output_tokens",
                      "sequence_query_count", "wall_time_seconds", "peak_memory_mb"):
            with self.assertRaises(ValueError):
                account.add_usage(field, 1)
        with self.assertRaises(ValueError):
            account.add_usage("no_such_field", 1)

    def test_record_heavy_stage_updates_ledger_and_manifest(self) -> None:
        account = BudgetAccount()
        record = account.record_heavy_stage(
            "large_llm_literature_review",
            "candidate_count_exceeds_wet_budget",
            wall_time_seconds=12.5,
            cpu_seconds=8.2,
            peak_memory_mb=256.0,
            metadata={"family_id": "FAM001"},
        )
        self.assertEqual("large_llm_literature_review", record["stage"])
        self.assertEqual("candidate_count_exceeds_wet_budget", record["invocation_reason"])
        self.assertEqual(1, account.freeze()["usage"]["heavy_stage_invocation_count"])
        self.assertEqual(1, len(account.heavy_stages()))
        self.assertIn("heavy_stages", account.freeze())
        self.assertEqual("large_llm_literature_review", account.freeze()["heavy_stages"][0]["stage"])

    def test_record_heavy_stage_requires_reason(self) -> None:
        account = BudgetAccount()
        with self.assertRaisesRegex(ValueError, "invocation_reason"):
            account.record_heavy_stage("structure_prediction", "")

    def test_exhausted_budgets_sorted(self) -> None:
        account = BudgetAccount({"llm_call_limit": 1, "token_limit": 2, "sequence_query_limit": 0,
                                 "wall_time_limit_seconds": 3600})
        self.assertTrue(account.consume("llm_call_budget", 1))
        self.assertTrue(account.consume("token_budget", 2))
        self.assertFalse(account.consume("llm_call_budget", 1))
        self.assertEqual(["llm_call_limit", "sequence_query_limit", "token_limit"], account.exhausted_budgets())
        self.assertNotIn("wall_time_limit_seconds", account.exhausted_budgets())


class BudgetAccountFiniteTests(unittest.TestCase):
    def test_consume_rejects_nan(self) -> None:
        account = BudgetAccount({"llm_call_limit": 10, "wall_time_limit_seconds": 60})
        with self.assertRaises(ValueError):
            account.consume("llm_call_budget", float("nan"))
        with self.assertRaises(ValueError):
            account.consume("wall_time_budget", float("nan"))
        self.assertEqual(0, account.freeze()["usage"]["llm_call_count"])
        self.assertEqual(0, account.freeze()["usage"]["wall_time_seconds"])
        self.assertEqual(10, account.remaining("llm_call_budget"))
        self.assertEqual(60, account.remaining("wall_time_budget"))

    def test_consume_rejects_infinite(self) -> None:
        account = BudgetAccount({"wall_time_limit_seconds": 60, "peak_memory_limit_mb": 1024})
        with self.assertRaises(ValueError):
            account.consume("wall_time_budget", float("inf"))
        with self.assertRaises(ValueError):
            account.consume("wall_time_budget", float("-inf"))
        with self.assertRaises(ValueError):
            account.observe_peak("peak_memory_budget", float("inf"))
        self.assertEqual(0, account.freeze()["usage"]["wall_time_seconds"])
        self.assertEqual(0, account.freeze()["usage"]["peak_memory_mb"])

    def test_add_usage_rejects_nan_for_float_fields(self) -> None:
        account = BudgetAccount()
        with self.assertRaises(ValueError):
            account.add_usage("sequence_search_cpu_seconds", float("nan"))
        with self.assertRaises(ValueError):
            account.add_usage("sequence_search_cpu_seconds", float("inf"))
        with self.assertRaises(ValueError):
            account.add_usage("sequence_search_cpu_seconds", float("-inf"))
        self.assertEqual(0, account.freeze()["usage"]["sequence_search_cpu_seconds"])
        with self.assertRaises(ValueError):
            account.add_usage("cache_hit_count", float("nan"))
        self.assertEqual(0, account.freeze()["usage"]["cache_hit_count"])


if __name__ == "__main__":
    unittest.main()
