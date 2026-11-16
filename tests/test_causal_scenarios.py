from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from candidate_qtl.causal_scenarios import (
    ASSESSMENT_STATES,
    SCENARIO_COLUMNS,
    build_scenario_ledger,
    normalize_scenario,
    render_scenario_report,
    scenario_summary,
    write_scenario_ledger,
)


class CausalScenarioLedgerTests(unittest.TestCase):
    def test_coding_sv_and_mapping_scenarios_can_coexist(self) -> None:
        rows = [
            {
                "scenario_type": "coding_gene",
                "supporting_evidence": "candidate coding variant",
                "assessment_state": "supported",
                "recommended_next_experiment": "sequence candidate CDS",
                "reason": "coding hypothesis remains plausible",
            },
            {
                "scenario_type": "sv_te_cnv_pav",
                "assessment_state": "unassessed",
                "recommended_next_experiment": "inspect long-read/SV evidence",
                "reason": "coding route does not test structural variation",
            },
            {
                "scenario_type": "mapping_or_liftover_error",
                "contradicting_evidence": "flanking marker order is concordant",
                "assessment_state": "contradicted",
                "recommended_next_experiment": "review marker pairs only if discordance appears",
                "reason": "mapping failure is currently disfavored",
            },
        ]
        ledger = build_scenario_ledger(rows, case_id="case1")
        self.assertEqual(3, len(ledger))
        self.assertEqual(
            {"coding_gene", "sv_te_cnv_pav", "mapping_or_liftover_error"},
            set(ledger.scenario_type),
        )
        self.assertTrue(ledger.rank_contribution.eq("none").all())
        self.assertEqual(3, ledger.scenario_id.nunique())

    def test_evidence_absence_remains_unassessed_not_negative(self) -> None:
        ledger = build_scenario_ledger([
            {
                "scenario_type": "regulatory_element",
                "recommended_next_experiment": "promoter review",
                "reason": "noncoding space not evaluated",
            }
        ], case_id="case1")
        row = ledger.iloc[0]
        self.assertEqual("unassessed", row.assessment_state)
        self.assertEqual("", row.supporting_evidence)
        self.assertEqual("", row.contradicting_evidence)
        with self.assertRaisesRegex(ValueError, "must remain unassessed"):
            build_scenario_ledger([
                {"scenario_type": "regulatory_element", "assessment_state": "contradicted"}
            ], case_id="case1")

    def test_summary_distinguishes_support_contradiction_and_unassessed(self) -> None:
        ledger = build_scenario_ledger([
            {"scenario_type": "coding_gene", "supporting_evidence": "E1", "assessment_state": "supported"},
            {"scenario_type": "coding_gene", "contradicting_evidence": "E2", "assessment_state": "contradicted"},
            {"scenario_type": "coding_gene"},
        ], case_id="case1")
        summary = scenario_summary(ledger).iloc[0]
        self.assertEqual(3, summary.scenario_count)
        self.assertEqual(1, summary.supported_rows)
        self.assertEqual(1, summary.contradicted_rows)
        self.assertEqual(1, summary.unassessed_rows)
        self.assertEqual(1, summary.support_present_rows)
        self.assertEqual(1, summary.contradiction_present_rows)

    def test_report_never_selects_a_winner(self) -> None:
        ledger = build_scenario_ledger([
            {"scenario_type": "coding_gene", "supporting_evidence": "E1", "assessment_state": "supported",
             "recommended_next_experiment": "CDS assay", "reason": "coding evidence"},
            {"scenario_type": "environment_specific", "assessment_state": "unassessed",
             "recommended_next_experiment": "repeat environments", "reason": "environment not tested"},
        ], case_id="case1")
        report = render_scenario_report(ledger)
        self.assertIn("`coding_gene` / `supported`", report)
        self.assertIn("`environment_specific` / `unassessed`", report)
        self.assertNotIn("winner", report.casefold())

    def test_write_uses_canonical_columns(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "locus_causal_scenarios.tsv"
            write_scenario_ledger([{"scenario_type": "sv_te_cnv_pav"}], path, case_id="case1")
            header = path.read_text(encoding="utf-8").splitlines()[0].split("\t")
        self.assertEqual(SCENARIO_COLUMNS, header)

    def test_invalid_scenario_type_fails_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "unsupported scenario_type"):
            build_scenario_ledger([{"scenario_type": "unknown_scenario"}], case_id="case1")

    def test_not_evaluated_state_is_accepted_without_evidence(self) -> None:
        self.assertIn("not_evaluated", ASSESSMENT_STATES)
        ledger = build_scenario_ledger([
            {"scenario_type": "coding_gene", "assessment_state": "not_evaluated",
             "recommended_next_experiment": "provide route A provider",
             "reason": "route_a provider missing; discordance cannot be classified"},
        ], case_id="case1")
        row = ledger.iloc[0]
        self.assertEqual("not_evaluated", row.assessment_state)
        self.assertEqual("", row.supporting_evidence)
        self.assertEqual("", row.contradicting_evidence)
        self.assertEqual("", normalize_scenario(
            {"scenario_type": "coding_gene", "assessment_state": "not_evaluated"},
            case_id="case1",
        )["supporting_evidence"])

    def test_evidence_free_supported_still_fails_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "must remain unassessed"):
            build_scenario_ledger([
                {"scenario_type": "coding_gene", "assessment_state": "supported"}
            ], case_id="case1")

    def test_issue_363_review_types_accepted_unknown_still_fails(self) -> None:
        for scenario_type in (
            "novel_or_unannotated_mechanism",
            "phenotype_context_review",
            "mapping_adequacy_review",
            "provider_missing",
        ):
            with self.subTest(scenario_type=scenario_type):
                ledger = build_scenario_ledger(
                    [{"scenario_type": scenario_type}], case_id="case1",
                )
                self.assertEqual(scenario_type, ledger.iloc[0].scenario_type)
                self.assertEqual("unassessed", ledger.iloc[0].assessment_state)
        with self.assertRaisesRegex(ValueError, "unsupported scenario_type"):
            build_scenario_ledger(
                [{"scenario_type": "unknown_scenario"}], case_id="case1",
            )

    def test_case_mismatch_fails_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "does not match"):
            build_scenario_ledger([{"case_id": "case2", "scenario_type": "coding_gene"}], case_id="case1")

    def test_reloaded_tsv_blanks_are_not_evidence(self) -> None:
        rows = [
            {"scenario_type": "coding_gene", "supporting_evidence": "E1", "assessment_state": "supported"},
            {"scenario_type": "regulatory_element", "assessment_state": "unassessed"},
        ]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "locus_causal_scenarios.tsv"
            write_scenario_ledger(rows, path, case_id="case1")
            reloaded = pd.read_csv(path, sep="\t")
        summary = scenario_summary(reloaded).set_index("scenario_type")
        self.assertEqual(1, summary.loc["coding_gene", "support_present_rows"])
        self.assertEqual(0, summary.loc["coding_gene", "contradiction_present_rows"])
        self.assertEqual(0, summary.loc["regulatory_element", "support_present_rows"])
        self.assertEqual(0, summary.loc["regulatory_element", "contradiction_present_rows"])
        self.assertEqual(1, summary.loc["regulatory_element", "unassessed_rows"])
        report = render_scenario_report(reloaded)
        self.assertNotIn("nan", report)
        self.assertNotIn("<NA>", report)

    def test_pandas_missing_sentinels_do_not_bypass_guard(self) -> None:
        sentinels = [pd.NA, pd.NaT, np.nan, float("nan"), None, "", np.array(np.nan)]
        for sentinel in sentinels:
            with self.subTest(sentinel=repr(sentinel)):
                normalized = normalize_scenario(
                    {"scenario_type": "coding_gene", "supporting_evidence": sentinel},
                    case_id="case1",
                )
                self.assertEqual("", normalized["supporting_evidence"])
        for sentinel in sentinels:
            with self.subTest(guard=sentinel):
                with self.assertRaisesRegex(ValueError, "must remain unassessed"):
                    build_scenario_ledger(
                        [{"scenario_type": "coding_gene", "supporting_evidence": sentinel,
                          "assessment_state": "supported"}],
                        case_id="case1",
                    )
        frame = pd.DataFrame(
            [{"scenario_type": "coding_gene", "supporting_evidence": pd.NA,
              "assessment_state": "supported"}],
            dtype=object,
        )
        with self.assertRaisesRegex(ValueError, "must remain unassessed"):
            build_scenario_ledger(frame, case_id="case1")

    def test_multi_case_report_shows_case_and_scenario_ids(self) -> None:
        rows = [
            {"case_id": "case1", "scenario_type": "coding_gene",
             "supporting_evidence": "E1", "assessment_state": "supported"},
            {"case_id": "case2", "scenario_type": "coding_gene",
             "supporting_evidence": "E2", "assessment_state": "supported"},
        ]
        ledger = build_scenario_ledger(rows)
        report = render_scenario_report(ledger)
        scenario_by_case = dict(zip(ledger.case_id, ledger.scenario_id, strict=True))
        bullets = [line for line in report.splitlines() if line.startswith("- ")]
        self.assertEqual(2, len(bullets))
        for bullet in bullets:
            case_id = bullet[2:].split("/", 1)[0]
            self.assertIn(case_id, scenario_by_case)
            scenario_id = scenario_by_case[case_id]
            self.assertIn(f"{case_id}/{scenario_id}", bullet)
            self.assertEqual(1, bullet.count(scenario_id))

    def test_write_read_write_roundtrip_is_stable(self) -> None:
        rows = [
            {"scenario_type": "coding_gene", "supporting_evidence": "E1",
             "assessment_state": "supported", "recommended_next_experiment": "CDS assay",
             "reason": "coding evidence"},
            {"scenario_type": "regulatory_element", "assessment_state": "unassessed"},
        ]
        with tempfile.TemporaryDirectory() as tmp:
            first = Path(tmp) / "first.tsv"
            second = Path(tmp) / "second.tsv"
            write_scenario_ledger(rows, first, case_id="case1")
            reloaded = pd.read_csv(first, sep="\t")
            write_scenario_ledger(reloaded, second, case_id="case1")
            self.assertEqual(first.read_bytes(), second.read_bytes())
            rebuilt = pd.read_csv(second, sep="\t")
            self.assertEqual(reloaded.scenario_id.tolist(), rebuilt.scenario_id.tolist())
            reg_original = reloaded[reloaded.scenario_type == "regulatory_element"].iloc[0]
            reg_rebuilt = rebuilt[rebuilt.scenario_type == "regulatory_element"].iloc[0]
            self.assertEqual("unassessed", reg_rebuilt.assessment_state)
            self.assertTrue(pd.isna(reg_original.supporting_evidence))
            self.assertTrue(pd.isna(reg_rebuilt.supporting_evidence))


if __name__ == "__main__":
    unittest.main()
