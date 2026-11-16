from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main


class ExpressionSpecificityTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.candidates = self.root / "candidates.tsv"
        self.candidates.write_text("gene_id\nG1\n", encoding="utf-8")

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def _build(self, expression_text: str):
        expression = self.root / "expr.tsv"
        expression.write_text(expression_text, encoding="utf-8")
        out = self.root / "ev"
        self.assertEqual(0, main(["evidence-build", "--candidates", str(self.candidates), "--expression", str(expression), "--out-dir", str(out)]))
        evidence = pd.read_csv(out / "evidence.tsv", sep="\t").fillna("")
        summary = pd.read_csv(out / "expression_summary.tsv", sep="\t").fillna("")
        return evidence, summary

    def test_target_only_does_not_claim_specificity(self) -> None:
        evidence, summary = self._build(
            "gene_id\ttissue\tvalue\tunit\tis_target\n"
            "G1\tmeristem\t20\tTPM\ttrue\n"
        )
        codes = set(evidence[evidence["gene_id"] == "G1"]["evidence_code"])
        self.assertIn("expression_detected", codes)
        self.assertNotIn("expression_target_specific", codes)
        self.assertIn("expression_comparison_unavailable", codes)
        # The context flag must not be scored as support.
        flag = evidence[evidence["evidence_code"] == "expression_comparison_unavailable"].iloc[0]
        self.assertEqual(flag["direction"], "context")
        row = summary[summary["gene_id"] == "G1"].iloc[0]
        self.assertEqual(row["specificity_class"], "comparison_unavailable")
        self.assertEqual(row["comparison_available"], "no")

    def test_target_plus_background_can_be_specific(self) -> None:
        evidence, summary = self._build(
            "gene_id\ttissue\tvalue\tunit\tis_target\n"
            "G1\tmeristem\t20\tTPM\ttrue\n"
            "G1\tleaf\t2\tTPM\tfalse\n"
        )
        codes = set(evidence[evidence["gene_id"] == "G1"]["evidence_code"])
        self.assertIn("expression_target_specific", codes)
        self.assertNotIn("expression_comparison_unavailable", codes)
        self.assertEqual(summary.iloc[0]["specificity_class"], "target_specific")

    def test_target_only_is_not_top_on_expression_alone(self) -> None:
        # A target-only expression signal must not build a scored support family.
        _, _ = self._build("gene_id\ttissue\tvalue\tunit\tis_target\nG1\tmeristem\t20\tTPM\ttrue\n")
        out = self.root / "ev"
        rank = self.root / "rank"
        self.assertEqual(0, main(["evidence-rank", "--candidates", str(self.candidates), "--evidence", str(out / "evidence.tsv"), "--out-dir", str(rank)]))
        ranked = pd.read_csv(rank / "ranked_candidates.tsv", sep="\t").fillna("")
        g1 = ranked[ranked["gene_id"] == "G1"].iloc[0]
        self.assertNotEqual(g1["dry_group"], "top")


if __name__ == "__main__":
    unittest.main()
