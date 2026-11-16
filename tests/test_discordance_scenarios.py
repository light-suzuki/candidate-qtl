"""Tests for the issue #363 A/B discordance -> review-required scenario adapter."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.causal_scenarios import ASSESSMENT_STATES, build_scenario_ledger
from candidate_qtl.discordance_scenarios import (
    PROVIDER_MISSING_SCENARIO_TYPE,
    STATE_SCENARIO_TYPES,
    classify_comparison_state,
    discordance_review_rows,
)

EXPECTED_STATE_TYPES = {state: set(types) for state, types in STATE_SCENARIO_TYPES.items()}


def _dual_frame() -> pd.DataFrame:
    return pd.DataFrame([
        {
            "comparison_id": "cmp_000", "route_a_hypothesis_id": "ha1",
            "route_b_hypothesis_id": "hb1", "pathway_or_process": "pathway X",
            "gene_family": "fam1", "route_a_source_ids": "SRC1",
            "route_b_source_ids": "SRC1;SRC2", "shared_source_ids": "SRC1",
            "source_independence_status": "shared_source_chain",
        },
        {
            "comparison_id": "cmp_001", "route_a_hypothesis_id": "ha2",
            "route_b_hypothesis_id": "", "pathway_or_process": "pathway A",
            "gene_family": "fam2", "route_a_source_ids": "SRC3",
            "route_b_source_ids": "", "shared_source_ids": "",
            "source_independence_status": "not_applicable",
        },
        {
            "comparison_id": "cmp_002", "route_a_hypothesis_id": "",
            "route_b_hypothesis_id": "hb3", "pathway_or_process": "pathway B",
            "gene_family": "fam3", "route_a_source_ids": "",
            "route_b_source_ids": "SRC4", "shared_source_ids": "",
            "source_independence_status": "not_applicable",
        },
        {
            "comparison_id": "cmp_003", "route_a_hypothesis_id": "",
            "route_b_hypothesis_id": "", "pathway_or_process": "",
            "gene_family": "fam4", "route_a_source_ids": "",
            "route_b_source_ids": "", "shared_source_ids": "",
            "source_independence_status": "not_applicable",
        },
    ])


def _triple_frame() -> pd.DataFrame:
    return pd.DataFrame([
        {"comparison_id": "cmp_010", "four_state": "a_and_b_supported",
         "pathway_or_process": "pathway X", "gene_family": "fam1",
         "shared_source_ids": "SRC1", "source_dependence_status": "shared_source_chain"},
        {"comparison_id": "cmp_011", "four_state": "a_only_supported",
         "pathway_or_process": "pathway A", "gene_family": "fam2"},
        {"comparison_id": "cmp_012", "four_state": "b_only_supported",
         "pathway_or_process": "pathway B", "gene_family": "fam3"},
        {"comparison_id": "cmp_013", "four_state": "neither_supported",
         "pathway_or_process": "", "gene_family": "fam4"},
    ])


def _state_from_reason(reason: str) -> str:
    return reason.split(";", 1)[0].split("=", 1)[1]


class DiscordanceScenarioTests(unittest.TestCase):
    def test_four_states_classify_deterministically(self) -> None:
        for frame in (_dual_frame(), _triple_frame()):
            with self.subTest(schema=tuple(frame.columns)):
                rows = discordance_review_rows(
                    frame, case_id="case1", input_fingerprint="fp-123",
                )
                by_state: dict[str, list[dict[str, str]]] = {}
                for row in rows:
                    by_state.setdefault(_state_from_reason(row["reason"]), []).append(row)
                self.assertEqual(
                    {"A∩B", "A-only", "B-only", "neither"}, set(by_state),
                )
                for state, state_rows in by_state.items():
                    self.assertEqual(
                        len(state_rows),
                        len(STATE_SCENARIO_TYPES[state]),
                        state,
                    )
                    self.assertEqual(
                        EXPECTED_STATE_TYPES[state],
                        {row["scenario_type"] for row in state_rows},
                        state,
                    )
                self.assertEqual(rows, discordance_review_rows(
                    frame, case_id="case1", input_fingerprint="fp-123",
                ))

    def test_tsv_path_input_matches_dataframe_input(self) -> None:
        frame = _dual_frame()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "dual_route_hypothesis_comparison.tsv"
            frame.to_csv(path, sep="\t", index=False, lineterminator="\n")
            from_path = discordance_review_rows(
                path, case_id="case1", input_fingerprint="fp-123",
            )
        self.assertEqual(
            from_path,
            discordance_review_rows(frame, case_id="case1", input_fingerprint="fp-123"),
        )

    def test_all_route_rows_are_review_required_not_excluded_or_promoted(self) -> None:
        rows = discordance_review_rows(
            _dual_frame(), case_id="case1", input_fingerprint="fp-123",
        )
        self.assertTrue(rows)
        for row in rows:
            self.assertEqual("review_required", row["assessment_state"])
            self.assertIn("review_required", ASSESSMENT_STATES)
            self.assertNotIn("rank_contribution", row)
            self.assertNotIn("exclud", row["reason"].casefold())
            self.assertNotIn("promot", row["reason"].casefold())
        ledger = build_scenario_ledger(rows, case_id="case1")
        self.assertTrue(ledger.rank_contribution.eq("none").all())
        self.assertTrue(ledger.assessment_state.eq("review_required").all())

    def test_provider_missing_is_not_evaluated_not_neither(self) -> None:
        a_only = _dual_frame().iloc[[1]].reset_index(drop=True)
        rows = discordance_review_rows(
            a_only, case_id="case1", input_fingerprint="fp-123", route_a_evaluated=False,
        )
        self.assertEqual(1, len(rows))
        row = rows[0]
        self.assertEqual("not_evaluated", row["assessment_state"])
        self.assertEqual(PROVIDER_MISSING_SCENARIO_TYPE, row["scenario_type"])
        self.assertEqual("", row["supporting_evidence"])
        self.assertIn("provider missing", row["reason"])
        self.assertIn("route_a", row["reason"])
        self.assertEqual(
            "provide route A provider or record route A as intentionally absent",
            row["recommended_next_experiment"],
        )

        b_only = _dual_frame().iloc[[2]].reset_index(drop=True)
        row_b = discordance_review_rows(
            b_only, case_id="case1", input_fingerprint="fp-123", route_b_evaluated=False,
        )[0]
        self.assertEqual("not_evaluated", row_b["assessment_state"])
        self.assertIn("route_b provider missing", row_b["reason"])

        both = discordance_review_rows(
            _triple_frame(), case_id="case1", input_fingerprint="fp-123",
            route_a_evaluated=False, route_b_evaluated=False,
        )
        self.assertEqual(4, len(both))
        self.assertTrue(all(
            r["assessment_state"] == "not_evaluated" for r in both
        ))
        self.assertIn("route_a and route_b providers missing", both[0]["reason"])

    def test_same_root_non_independent_convergence_is_recorded_not_forced(self) -> None:
        frame = pd.DataFrame([{
            "comparison_id": "cmp_000", "route_a_hypothesis_id": "ha1",
            "route_b_hypothesis_id": "hb1", "pathway_or_process": "pathway X",
            "gene_family": "fam1", "route_a_source_ids": "SRC1",
            "route_b_source_ids": "SRC1;SRC2", "shared_source_ids": "SRC1;SRC2",
            "source_independence_status": "shared_source_chain",
        }])
        rows = discordance_review_rows(
            frame, case_id="case1", input_fingerprint="fp-123",
        )
        self.assertEqual(2, len(rows))
        self.assertEqual(
            {"coding_gene", "regulatory_element"},
            {row["scenario_type"] for row in rows},
        )
        for row in rows:
            self.assertEqual("review_required", row["assessment_state"])
            self.assertIn("shared_source_ids=SRC1;SRC2", row["source_ids"])
            self.assertIn("same-root convergence", row["reason"])
            self.assertIn("shared_source_chain", row["reason"])

    def test_source_ids_only_differences_keep_distinct_scenario_ids(self) -> None:
        frame = pd.DataFrame([
            {"comparison_id": "cmp_000", "route_a_hypothesis_id": "ha1",
             "route_b_hypothesis_id": "hb1", "pathway_or_process": "pathway X"},
            {"comparison_id": "cmp_001", "route_a_hypothesis_id": "ha2",
             "route_b_hypothesis_id": "hb2", "pathway_or_process": "pathway X"},
        ])
        rows = discordance_review_rows(
            frame, case_id="case1", input_fingerprint="fp-123", route_b_evaluated=False,
        )
        self.assertEqual(2, len(rows))
        self.assertTrue(all(row["scenario_type"] == PROVIDER_MISSING_SCENARIO_TYPE for row in rows))
        ledger = build_scenario_ledger(rows, case_id="case1")
        self.assertEqual(2, len(ledger))
        self.assertEqual(2, ledger.scenario_id.nunique())
        self.assertEqual("not_evaluated", ledger.iloc[0].assessment_state)
        duplicate = rows + [dict(rows[0])]
        deduped = build_scenario_ledger(duplicate, case_id="case1")
        self.assertEqual(2, len(deduped))

    def test_input_order_does_not_change_output(self) -> None:
        frame = _dual_frame()
        reversed_frame = frame.iloc[::-1].reset_index(drop=True)
        ledger = build_scenario_ledger(
            discordance_review_rows(frame, case_id="case1", input_fingerprint="fp-123"),
            case_id="case1",
        )
        reversed_ledger = build_scenario_ledger(
            discordance_review_rows(reversed_frame, case_id="case1", input_fingerprint="fp-123"),
            case_id="case1",
        )
        self.assertTrue(ledger.equals(reversed_ledger))
        self.assertEqual(ledger.to_csv(index=False), reversed_ledger.to_csv(index=False))

    def test_unknown_comparison_schema_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "not_comparison.tsv"
            path.write_text(
                "family_id\tconvergence_status\nfam1\tphenotype_only_supported\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "not a canonical comparison artifact"):
                discordance_review_rows(path, case_id="case1", input_fingerprint="")
        frame = pd.DataFrame([
            {"comparison_id": "cmp_000", "four_state": "unexpected_value"},
        ])
        with self.assertRaisesRegex(ValueError, "unsupported four_state"):
            discordance_review_rows(frame, case_id="case1", input_fingerprint="")
        with self.assertRaisesRegex(ValueError, "unsupported four_state"):
            classify_comparison_state({"four_state": "unexpected_value"})

    def test_empty_comparison_returns_empty(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "empty.tsv"
            path.write_text(
                "comparison_id\troute_a_hypothesis_id\troute_b_hypothesis_id\n",
                encoding="utf-8",
            )
            self.assertEqual(
                [], discordance_review_rows(path, case_id="case1", input_fingerprint=""),
            )
        self.assertEqual(
            [], discordance_review_rows(
                pd.DataFrame(columns=["four_state"]), case_id="case1", input_fingerprint="",
            ),
        )

    def test_fingerprint_and_comparison_id_are_recorded(self) -> None:
        rows = discordance_review_rows(
            _dual_frame(), case_id="case1", input_fingerprint="sha256-abc123",
        )
        self.assertTrue(rows)
        for row in rows:
            self.assertIn("comparison=cmp_0", row["source_ids"])
            self.assertIn("fingerprint=sha256-abc123", row["source_ids"])
        self.assertIn(
            "route_a_source_ids=SRC1",
            rows[0]["source_ids"],
        )


if __name__ == "__main__":
    unittest.main()
