from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.qtl_effect_contract import (
    QTL_EFFECT_COLUMNS,
    annotate_sign_concordance,
    contract_frame,
    load_qtl_effect_contract,
    normalize_qtl_effect_contract,
)


class QtlEffectContractTests(unittest.TestCase):
    def test_missing_contract_is_backward_compatible_and_explicit_unknown(self) -> None:
        contract = load_qtl_effect_contract(None, case_id="legacy_case")
        self.assertEqual("legacy_case", contract["case_id"])
        for column in QTL_EFFECT_COLUMNS:
            if column != "case_id":
                self.assertEqual("unknown", contract[column])

    def test_yaml_nested_contract_preserves_qtl_context_and_source(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "case.yaml"
            path.write_text(
                "case_id: pea_height\n"
                "qtl_effect_contract:\n"
                "  favorable_parent: ParentA\n"
                "  effect_direction: increase\n"
                "  additive_effect: '3.2'\n"
                "  dominance_class: additive\n"
                "  effect_ci: '2.1..4.3'\n"
                "  environment_id: field_2026\n"
                "  phenotyping_protocol: mature stem height\n"
                "  stable_across_environments: true\n"
                "  source_id: PMID:fixture\n"
                "  source_detail: Table 2\n",
                encoding="utf-8",
            )
            contract = load_qtl_effect_contract(path)
        self.assertEqual("pea_height", contract["case_id"])
        self.assertEqual("ParentA", contract["favorable_parent"])
        self.assertEqual("increase", contract["effect_direction"])
        self.assertEqual("true", contract["stable_across_environments"])
        self.assertEqual("PMID:fixture", contract["source_id"])

    def test_tsv_selects_exact_case_and_rejects_ambiguous_multirow_input(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "effects.tsv"
            path.write_text(
                "case_id\teffect_direction\tsource_id\n"
                "positive_case\tincrease\tPMID:1\n"
                "negative_case\tdecrease\tPMID:2\n",
                encoding="utf-8",
            )
            selected = load_qtl_effect_contract(path, case_id="negative_case")
            self.assertEqual("decrease", selected["effect_direction"])
            with self.assertRaisesRegex(ValueError, "exactly one row"):
                load_qtl_effect_contract(path)

    def test_case_mismatch_fails_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "does not match"):
            normalize_qtl_effect_contract({"case_id": "case_a"}, case_id="case_b")

    def test_yaml_top_level_case_mismatch_fails_closed_before_normalizing(self) -> None:
        # [P1] The nested block omits its own id; only the top-level case_id is
        # declared. A caller requesting a different case must fail closed instead
        # of silently attaching case_a's QTL-effect claim to case_b.
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "case.yaml"
            path.write_text(
                "case_id: case_a\n"
                "qtl_effect_contract:\n"
                "  effect_direction: increase\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "does not match"):
                load_qtl_effect_contract(path, case_id="case_b")
            # Matching ids still load normally.
            contract = load_qtl_effect_contract(path, case_id="case_a")
            self.assertEqual("case_a", contract["case_id"])
            self.assertEqual("increase", contract["effect_direction"])

    def test_header_only_tsv_with_requested_case_fails_closed(self) -> None:
        # [P2] A header-only TSV resolves no row; a requested case is a
        # fail-closed mismatch, while no caller case keeps the all-unknown
        # backward-compatible behaviour.
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "effects.tsv"
            path.write_text("case_id\teffect_direction\tsource_id\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "no row for case_id"):
                load_qtl_effect_contract(path, case_id="case_b")
            contract = load_qtl_effect_contract(path)
            self.assertEqual("unknown", contract["effect_direction"])

    def test_falsy_non_mapping_yaml_roots_fail_closed(self) -> None:
        for content in ("[]\n", "false\n", "0\n", '""\n'):
            with self.subTest(content=content.strip()):
                with tempfile.TemporaryDirectory() as tmp:
                    path = Path(tmp) / "effect.yaml"
                    path.write_text(content, encoding="utf-8")
                    with self.assertRaisesRegex(ValueError, "must contain an object"):
                        load_qtl_effect_contract(path, case_id="case_a")

    def test_empty_yaml_document_remains_backward_compatible(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "effect.yaml"
            path.write_text("", encoding="utf-8")
            contract = load_qtl_effect_contract(path, case_id="case_a")
        self.assertEqual("case_a", contract["case_id"])
        self.assertEqual("unknown", contract["effect_direction"])

    def test_scalar_missing_values_normalize_to_unknown(self) -> None:
        # [P2] pd.NA / pd.NaT / nan / None are all scalar missing markers:
        # stable_across_environments must not raise and every field must
        # normalize to the literal "unknown".
        for missing in (pd.NA, pd.NaT, float("nan"), None):
            contract = normalize_qtl_effect_contract(
                {"stable_across_environments": missing, "additive_effect": missing, "effect_direction": missing}
            )
            self.assertEqual("unknown", contract["stable_across_environments"])
            self.assertEqual("unknown", contract["additive_effect"])
            self.assertEqual("unknown", contract["effect_direction"])

    def test_list_values_are_not_treated_as_missing(self) -> None:
        # [P2] Array-like values are never missing: a YAML list must be kept as
        # its text form, not silently rewritten to "unknown".
        contract = normalize_qtl_effect_contract({"effect_ci": ["2.1", "4.3"]})
        self.assertEqual("['2.1', '4.3']", contract["effect_ci"])
        self.assertEqual("unknown", contract["additive_effect"])

    def test_invalid_stability_value_fails_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "stable_across_environments"):
            normalize_qtl_effect_contract({"stable_across_environments": "sometimes"})

    def test_contract_frame_uses_canonical_column_order(self) -> None:
        frame = contract_frame({"case_id": "c1", "effect_direction": "gain"})
        self.assertEqual(QTL_EFFECT_COLUMNS, list(frame.columns))
        self.assertEqual("gain", frame.iloc[0]["effect_direction"])


class QtlEffectSignConcordanceTests(unittest.TestCase):
    def test_positive_incompatible_and_unknown_fixtures_are_audit_only(self) -> None:
        hypotheses = pd.DataFrame(
            [
                {"hypothesis_id": "H_positive", "expected_mechanism_sign": "increase", "dry_score": "7.5", "dry_rank": "1"},
                {"hypothesis_id": "H_incompatible", "expected_mechanism_sign": "decrease", "dry_score": "4.0", "dry_rank": "2"},
                {"hypothesis_id": "H_unknown", "expected_mechanism_sign": "unknown", "dry_score": "1.0", "dry_rank": "3"},
            ]
        )
        before_scores = hypotheses[["hypothesis_id", "dry_score", "dry_rank"]].copy()
        contract = normalize_qtl_effect_contract(
            {
                "case_id": "case1",
                "effect_direction": "increase",
                "environment_id": "field_2026",
                "source_id": "PMID:123",
                "source_detail": "supplementary table S4",
            }
        )
        result = annotate_sign_concordance(hypotheses, contract).set_index("hypothesis_id")

        self.assertEqual("concordant", result.loc["H_positive", "sign_concordance"])
        self.assertEqual("incompatible", result.loc["H_incompatible", "sign_concordance"])
        self.assertEqual("unknown", result.loc["H_unknown", "sign_concordance"])
        self.assertEqual("PMID:123", result.loc["H_incompatible", "qtl_effect_source_id"])
        self.assertEqual("field_2026", result.loc["H_incompatible", "qtl_effect_environment_id"])
        self.assertIn("qtl_effect_sign=positive", result.loc["H_incompatible", "sign_concordance_reason"])

        after_scores = result.reset_index()[["hypothesis_id", "dry_score", "dry_rank"]]
        pd.testing.assert_frame_equal(before_scores.reset_index(drop=True), after_scores.reset_index(drop=True))

    def test_missing_expected_mechanism_sign_remains_unknown(self) -> None:
        hypotheses = pd.DataFrame([{"hypothesis_id": "H1"}])
        result = annotate_sign_concordance(
            hypotheses,
            normalize_qtl_effect_contract({"effect_direction": "decrease", "source_id": "local:effect"}),
        )
        self.assertEqual("unknown", result.iloc[0]["expected_mechanism_sign"])
        self.assertEqual("unknown", result.iloc[0]["sign_concordance"])


if __name__ == "__main__":
    unittest.main()
