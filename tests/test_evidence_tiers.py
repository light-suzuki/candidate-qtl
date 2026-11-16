"""Biological evidence tiers stay independent of the Stage-A assay budget (issue #229)."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import build_parser, main
from candidate_qtl.commands import handler_for
from candidate_qtl.evidence_tiers import evidence_tier_budget_benchmark, evidence_tier_report


def _eligibility_table(pairs: list[tuple[str, str]]) -> pd.DataFrame:
    """Ranked rows: one (gene_id, top_eligibility) pair per screening rank."""
    return pd.DataFrame([
        {"gene_id": gene_id, "dry_score": f"{len(pairs) - rank + 1}.0",
         "screening_rank": str(rank), "top_eligibility": eligibility}
        for rank, (gene_id, eligibility) in enumerate(pairs, start=1)
    ])


class EvidenceTierReportTests(unittest.TestCase):
    def _write_ranked(self, root: Path, name: str, pairs: list[tuple[str, str]]) -> Path:
        path = root / name
        _eligibility_table(pairs).to_csv(path, sep="\t", index=False)
        return path

    def _run_report(self, root: Path, out_name: str, ranked: Path, budget: int,
                    coverage: Path | None = None) -> tuple[pd.DataFrame, dict]:
        out = root / out_name
        argv = ["evidence-tiers", "--ranked", str(ranked), "--out-dir", str(out),
                "--review-budget", str(budget)]
        if coverage is not None:
            argv += ["--knowledge-coverage", str(coverage)]
        self.assertEqual(0, main(argv))
        report = pd.read_csv(out / "evidence_tiers.tsv", sep="\t", dtype=str).fillna("")
        manifest = json.loads((out / "evidence_tiers_manifest.json").read_text(encoding="utf-8"))
        return report, manifest

    def test_tier_is_identical_when_review_budget_changes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pairs = (
                [(f"Strong{index:02d}", "strong_genetic_evidence") for index in range(1, 6)]
                + [(f"Independent{index:02d}", "independent_experimental_support") for index in range(1, 6)]
                + [(f"Multi{index:02d}", "multiple_non_redundant_channels") for index in range(1, 6)]
                + [(f"Insufficient{index:02d}", "insufficient_independent_evidence") for index in range(1, 6)]
                + [(f"Conflicting{index:02d}", "conflicting_evidence") for index in range(1, 6)]
            )
            ranked = self._write_ranked(root, "ranked.tsv", pairs)
            report20, manifest20 = self._run_report(root, "budget20", ranked, 20)
            report40, manifest40 = self._run_report(root, "budget40", ranked, 40)

            self.assertEqual(report20["evidence_tier"].tolist(), report40["evidence_tier"].tolist())
            changed = report20[report20["operational_state"] != report40["operational_state"]]
            self.assertEqual(5, len(changed))
            self.assertEqual({"Conflicting01", "Conflicting02", "Conflicting03", "Conflicting04", "Conflicting05"},
                             set(changed["gene_id"]))
            self.assertEqual({"not_in_current_assay_budget"}, set(changed["operational_state"]))
            self.assertEqual({"stage_a_review_now"}, set(report40.loc[changed.index, "operational_state"]))
            self.assertEqual(20, manifest20["review_budget"])
            self.assertEqual(40, manifest40["review_budget"])
            self.assertEqual(20, manifest20["boundary"]["screening_rank_at_cut"])
            self.assertEqual(25, manifest40["boundary"]["screening_rank_at_cut"])

    def test_boundary_pair_with_the_same_tier_is_not_a_biological_difference(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pairs = ([(f"Gene{index:02d}", "strong_genetic_evidence") for index in range(1, 20)]
                     + [("Gene20", "insufficient_independent_evidence"), ("Gene21", "insufficient_independent_evidence")]
                     + [(f"Gene{index:02d}", "strong_genetic_evidence") for index in range(22, 26)])
            ranked = self._write_ranked(root, "ranked.tsv", pairs)
            report, manifest = self._run_report(root, "out", ranked, 20)

            comparison = manifest["boundary_tier_comparison"]
            self.assertTrue(comparison["evaluated"])
            self.assertEqual(20, comparison["rank_at_budget"])
            self.assertEqual(21, comparison["rank_after_budget"])
            self.assertEqual("unresolved", comparison["rank_at_budget_tier"])
            self.assertEqual("unresolved", comparison["rank_after_budget_tier"])
            self.assertTrue(comparison["same_biological_tier"])
            self.assertEqual("unresolved", report.iloc[19]["evidence_tier"])
            self.assertEqual("unresolved", report.iloc[20]["evidence_tier"])
            just_outside = pd.read_csv(root / "out" / "just_outside_budget.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(["Gene21"], just_outside["gene_id"].tolist())

    def test_top_eligibility_values_map_to_the_five_tiers(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pairs = [
                ("GeneStrong1", "strong_genetic_evidence"),
                ("GeneStrong2", "independent_experimental_support"),
                ("GenePlausible", "multiple_non_redundant_channels"),
                ("GeneContradictory", "conflicting_evidence"),
                ("GeneUnresolved", "insufficient_independent_evidence"),
            ]
            ranked = self._write_ranked(root, "ranked.tsv", pairs)
            report, manifest = self._run_report(root, "out", ranked, 20)
            by_gene = report.set_index("gene_id")["evidence_tier"]
            self.assertEqual("strong", by_gene["GeneStrong1"])
            self.assertEqual("strong", by_gene["GeneStrong2"])
            self.assertEqual("plausible", by_gene["GenePlausible"])
            self.assertEqual("contradictory", by_gene["GeneContradictory"])
            self.assertEqual("unresolved", by_gene["GeneUnresolved"])
            self.assertEqual(0, manifest["tier_counts"]["knowledge_limited"])

    def test_knowledge_limited_and_dark_coverage_are_knowledge_limited_tier(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pairs = [
                ("GeneLimited", "insufficient_independent_evidence"),
                ("GeneDark", "insufficient_independent_evidence"),
                ("GeneStudied", "insufficient_independent_evidence"),
                ("GeneStrongDark", "strong_genetic_evidence"),
            ]
            ranked = self._write_ranked(root, "ranked.tsv", pairs)
            coverage_path = root / "coverage.tsv"
            pd.DataFrame([
                {"gene_id": "GeneLimited", "knowledge_coverage_tier": "knowledge_limited"},
                {"gene_id": "GeneDark", "knowledge_coverage_tier": "dark"},
                {"gene_id": "GeneStudied", "knowledge_coverage_tier": "well_characterised"},
                {"gene_id": "GeneStrongDark", "knowledge_coverage_tier": "dark"},
            ]).to_csv(coverage_path, sep="\t", index=False)
            report, manifest = self._run_report(root, "out", ranked, 20, coverage=coverage_path)
            by_gene = report.set_index("gene_id")["evidence_tier"]
            self.assertEqual("knowledge_limited", by_gene["GeneLimited"])
            self.assertEqual("knowledge_limited", by_gene["GeneDark"])
            self.assertEqual("unresolved", by_gene["GeneStudied"])
            self.assertEqual("strong", by_gene["GeneStrongDark"])
            self.assertTrue(manifest["knowledge_coverage_provided"])
            self.assertEqual(2, manifest["tier_counts"]["knowledge_limited"])

    def test_knowledge_limited_tier_never_appears_without_coverage_input(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pairs = [(f"Gene{index:02d}", "insufficient_independent_evidence") for index in range(1, 6)]
            ranked = self._write_ranked(root, "ranked.tsv", pairs)
            report, manifest = self._run_report(root, "out", ranked, 20)
            self.assertEqual(["unresolved"], sorted(set(report["evidence_tier"])))
            self.assertEqual(0, manifest["tier_counts"]["knowledge_limited"])
            self.assertFalse(manifest["knowledge_coverage_provided"])

    def test_ranked_table_without_top_eligibility_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ranked = root / "ranked.tsv"
            pd.DataFrame([
                {"gene_id": "Gene01", "dry_score": "1.0", "screening_rank": "1"},
            ]).to_csv(ranked, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "top_eligibility"):
                evidence_tier_report(str(ranked), str(root / "out"))

    def test_report_writes_all_four_output_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pairs = [(f"Gene{index:02d}", "strong_genetic_evidence") for index in range(1, 6)]
            ranked = self._write_ranked(root, "ranked.tsv", pairs)
            out = root / "out"
            self.assertEqual(0, main(["evidence-tiers", "--ranked", str(ranked), "--out-dir", str(out)]))
            for name in ("evidence_tiers.tsv", "evidence_tier_budget_crosstab.tsv",
                         "just_outside_budget.tsv", "evidence_tiers_manifest.json"):
                self.assertTrue((out / name).is_file(), name)
            manifest = json.loads((out / "evidence_tiers_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(5, manifest["candidate_count"])
            self.assertEqual(5, manifest["tier_counts"]["strong"])
            crosstab = pd.read_csv(out / "evidence_tier_budget_crosstab.tsv", sep="\t")
            strong = crosstab[crosstab["evidence_tier"] == "strong"].iloc[0]
            self.assertEqual(5, strong["tier_total"])
            self.assertEqual(5, strong["in_budget_count"])


class EvidenceTierBenchmarkTests(unittest.TestCase):
    def _write_ranked(self, root: Path, name: str, pairs: list[tuple[str, str]]) -> Path:
        path = root / name
        _eligibility_table(pairs).to_csv(path, sep="\t", index=False)
        return path

    @staticmethod
    def _many_pairs(count: int) -> list[tuple[str, str]]:
        return [(f"Gene{index:02d}", "strong_genetic_evidence") for index in range(1, count + 1)]

    def _write_manifest(self, root: Path, rows: list[dict[str, object]]) -> Path:
        path = root / "manifest.tsv"
        pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
        return path

    def test_recall_at_budget_and_tier_recall_for_answer_in_and_out_of_the_cut(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case_a = self._write_ranked(root, "case_a.tsv", self._many_pairs(25))
            case_b = self._write_ranked(root, "case_b.tsv", self._many_pairs(25))
            manifest = self._write_manifest(root, [
                {"case_id": "case_a", "ranked_tsv": str(case_a), "answer_gene_id": "Gene01"},
                {"case_id": "case_b", "ranked_tsv": str(case_b), "answer_gene_id": "Gene25"},
            ])
            out = root / "bench"
            self.assertEqual(0, main(["evidence-tier-benchmark", "--manifest", str(manifest), "--out-dir", str(out)]))

            cases = pd.read_csv(out / "evidence_tier_benchmark_cases.tsv", sep="\t", dtype=str).fillna("").set_index("case_id")
            self.assertEqual("True", cases.loc["case_a", "recall_at_budget"])
            self.assertEqual("True", cases.loc["case_a", "tier_is_strong_or_plausible"])
            self.assertEqual("stage_a_review_now", cases.loc["case_a", "answer_operational_state"])
            self.assertEqual("strong", cases.loc["case_a", "answer_evidence_tier"])
            self.assertEqual("False", cases.loc["case_b", "recall_at_budget"])
            self.assertEqual("True", cases.loc["case_b", "tier_is_strong_or_plausible"])
            self.assertEqual("not_in_current_assay_budget", cases.loc["case_b", "answer_operational_state"])

            result_manifest = json.loads((out / "evidence_tier_benchmark_manifest.json").read_text(encoding="utf-8"))
            recall = result_manifest["metrics"]["recall_at_budget"]
            tier = result_manifest["metrics"]["tier_is_strong_or_plausible"]
            self.assertEqual(2, recall["evaluated"])
            self.assertEqual(1, recall["hits"])
            self.assertEqual(0.5, recall["recall"])
            self.assertEqual(2, tier["evaluated"])
            self.assertEqual(2, tier["hits"])
            self.assertEqual(1.0, tier["recall"])
            for name in ("evidence_tier_benchmark_cases.tsv", "evidence_tier_benchmark_strata.tsv",
                         "evidence_tier_benchmark_manifest.json"):
                self.assertTrue((out / name).is_file(), name)

    def test_answer_absent_from_candidates_is_a_real_miss_for_budget_recall(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ranked = self._write_ranked(root, "case.tsv", self._many_pairs(25))
            manifest = self._write_manifest(root, [
                {"case_id": "case_x", "ranked_tsv": str(ranked), "answer_gene_id": "NotACandidate"},
            ])
            out = root / "bench"
            self.assertEqual(0, main(["evidence-tier-benchmark", "--manifest", str(manifest), "--out-dir", str(out)]))
            cases = pd.read_csv(out / "evidence_tier_benchmark_cases.tsv", sep="\t", dtype=str).fillna("").set_index("case_id")
            self.assertEqual("False", cases.loc["case_x", "recall_at_budget"])
            self.assertEqual("False", cases.loc["case_x", "answer_present_in_candidates"])
            self.assertEqual("", cases.loc["case_x", "answer_evidence_tier"])
            self.assertEqual("", cases.loc["case_x", "tier_is_strong_or_plausible"])
            result_manifest = json.loads((out / "evidence_tier_benchmark_manifest.json").read_text(encoding="utf-8"))
            recall = result_manifest["metrics"]["recall_at_budget"]
            tier = result_manifest["metrics"]["tier_is_strong_or_plausible"]
            self.assertEqual(1, recall["evaluated"])
            self.assertEqual(0, recall["hits"])
            self.assertEqual(0, tier["evaluated"])
            self.assertEqual(1, tier["not_evaluated"])

    def test_missing_answer_gene_id_is_not_evaluated_not_a_miss(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ranked = self._write_ranked(root, "case.tsv", self._many_pairs(25))
            manifest = self._write_manifest(root, [{"case_id": "case_y", "ranked_tsv": str(ranked)}])
            out = root / "bench"
            self.assertEqual(0, main(["evidence-tier-benchmark", "--manifest", str(manifest), "--out-dir", str(out)]))
            cases = pd.read_csv(out / "evidence_tier_benchmark_cases.tsv", sep="\t", dtype=str).fillna("").set_index("case_id")
            self.assertEqual("", cases.loc["case_y", "recall_at_budget"])
            self.assertEqual("", cases.loc["case_y", "tier_is_strong_or_plausible"])
            result_manifest = json.loads((out / "evidence_tier_benchmark_manifest.json").read_text(encoding="utf-8"))
            recall = result_manifest["metrics"]["recall_at_budget"]
            self.assertEqual(0, recall["evaluated"])
            self.assertEqual(1, recall["not_evaluated"])
            self.assertIsNone(recall["recall"])

    def test_manifest_contract_errors_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ranked = self._write_ranked(root, "case.tsv", self._many_pairs(5))
            missing_column = self._write_manifest(root, [{"case_id": "case_z"}])
            with self.assertRaisesRegex(ValueError, "missing columns"):
                evidence_tier_budget_benchmark(str(missing_column), str(root / "out1"))
            duplicate = self._write_manifest(root, [
                {"case_id": "dup", "ranked_tsv": str(ranked)},
                {"case_id": "dup", "ranked_tsv": str(ranked)},
            ])
            with self.assertRaisesRegex(ValueError, "unique"):
                evidence_tier_budget_benchmark(str(duplicate), str(root / "out2"))


class EvidenceTierCommandRegistrationTests(unittest.TestCase):
    def test_evidence_tiers_commands_are_registered(self) -> None:
        parser = build_parser()
        choices = [action for action in parser._subparsers._group_actions if hasattr(action, "choices")][0].choices
        self.assertIn("evidence-tiers", choices)
        self.assertIn("evidence-tier-benchmark", choices)
        self.assertIsNotNone(handler_for("evidence-tiers"))
        self.assertIsNotNone(handler_for("evidence-tier-benchmark"))


if __name__ == "__main__":
    unittest.main()
