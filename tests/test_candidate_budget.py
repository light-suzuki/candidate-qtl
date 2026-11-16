from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.candidate_budget import budget_decision, candidate_budget_gate
from candidate_qtl.cli import main


class CandidateBudgetTests(unittest.TestCase):
    def _write_candidates(self, path: Path, count: int) -> None:
        pd.DataFrame({
            "candidate_gene_id": [f"G{index:03d}" for index in range(count)],
            "family_id": [f"F{index // 2:03d}" for index in range(count)],
            "rank_contribution": ["none"] * count,
        }).to_csv(path, sep="\t", index=False, lineterminator="\n")

    def test_empty_supported_set_is_explicit_abstention_not_success(self) -> None:
        decision = budget_decision(0, wet_budget=20, adaptive_round=0, max_adaptive_rounds=1)
        self.assertEqual("no_supported_candidates", decision["status"])
        self.assertEqual("stop_review_required", decision["action"])
        self.assertFalse(decision["adaptive_refinement_required"])
        self.assertEqual("candidate_set_empty", decision["decision_reason"])
        self.assertTrue(decision["review_required"])

    def test_at_budget_stops_without_adaptive_work(self) -> None:
        decision = budget_decision(20, wet_budget=20, adaptive_round=0, max_adaptive_rounds=1)
        self.assertEqual("within_wet_budget", decision["status"])
        self.assertEqual("freeze_candidate_set", decision["action"])
        self.assertFalse(decision["adaptive_refinement_required"])
        self.assertEqual("within_wet_budget", decision["decision_reason"])

    def test_over_budget_requests_only_one_bounded_refinement_round(self) -> None:
        first = budget_decision(21, wet_budget=20, adaptive_round=0, max_adaptive_rounds=1, budget_policy="strict")
        self.assertEqual("adaptive_refinement_required", first["status"])
        self.assertTrue(first["adaptive_refinement_required"])
        self.assertEqual("candidate_count_exceeds_wet_budget", first["decision_reason"])

        after_round = budget_decision(21, wet_budget=20, adaptive_round=1, max_adaptive_rounds=1, budget_policy="strict")
        self.assertEqual("review_required_over_budget", after_round["status"])
        self.assertEqual("stop_review_required", after_round["action"])
        self.assertFalse(after_round["adaptive_refinement_required"])
        self.assertEqual("adaptive_round_budget_exhausted", after_round["decision_reason"])
        self.assertTrue(after_round["review_required"])

    def test_gate_never_truncates_to_wet_budget(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "supported_members.tsv"
            self._write_candidates(source, 27)
            result = candidate_budget_gate(str(source), str(root / "gate"), wet_budget=20, case_id="case426")
            frozen = pd.read_csv(result["candidate_set"], sep="\t", dtype=str).fillna("")
            self.assertEqual(27, len(frozen))
            self.assertEqual({f"G{index:03d}" for index in range(27)}, set(frozen["candidate_gene_id"]))
            self.assertEqual("above_wet_target", result["status"])
            self.assertFalse(result["adaptive_refinement_required"])

            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual(27, manifest["candidate_count"])
            self.assertEqual(20, manifest["wet_budget"])
            self.assertEqual("retain_supported_members_above_operational_target", manifest["decision_reason"])
            self.assertEqual(7, manifest["above_target_count"])
            self.assertEqual("target", manifest["budget_policy"])
            self.assertNotIn("adaptive_synthesis_invoked", manifest)
            self.assertEqual("none", manifest["rank_contribution"])

    def test_within_budget_freezes_all_supported_copies_unchanged(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "supported_members.tsv"
            self._write_candidates(source, 17)
            original = source.read_bytes()
            result = candidate_budget_gate(str(source), str(root / "gate"), wet_budget=20)
            self.assertEqual("within_wet_budget", result["status"])
            self.assertFalse(result["adaptive_refinement_required"])
            self.assertEqual(original, Path(result["candidate_set"]).read_bytes())

    def test_cli_is_reachable_and_preserves_candidate_set(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "supported_members.tsv"
            self._write_candidates(source, 17)
            self.assertEqual(
                0,
                main([
                    "candidate-budget-gate",
                    "--candidates", str(source),
                    "--out-dir", str(root / "out"),
                    "--wet-budget", "20",
                    "--case-id", "cli-case",
                ]),
            )
            manifest = json.loads((root / "out" / "candidate_budget_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("cli-case", manifest["case_id"])
            self.assertEqual("within_wet_budget", manifest["status"])
            self.assertEqual(17, manifest["candidate_count"])

    def test_duplicate_candidate_ids_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "duplicates.tsv"
            pd.DataFrame({"candidate_gene_id": ["G1", "G1"]}).to_csv(
                source, sep="\t", index=False, lineterminator="\n"
            )
            with self.assertRaisesRegex(ValueError, "duplicate"):
                candidate_budget_gate(str(source), str(root / "gate"))

    def test_invalid_round_configuration_fails_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "cannot exceed"):
            budget_decision(30, wet_budget=20, adaptive_round=2, max_adaptive_rounds=1)


if __name__ == "__main__":
    unittest.main()
