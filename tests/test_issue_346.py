"""Issue #346: an unwired evidence lane must not silently degenerate the rank.

When inputs.evidence_resources is absent or empty, case-run still completes
(REPORT_READY, exit 0) but every candidate scores 0.0 and screening_rank is a
full tie that falls back to coordinate order. The run-level manifests
(case_state.json / run_manifest.json) must now say so explicitly:
``rank_semantic_status`` / ``rank_is_positional_fallback`` /
``rank_metric_eligible`` / ``stage_a_is_budget_only`` plus a warning, and the
guided-profile plan/execute must surface the same state in its manifest.
"""

from __future__ import annotations

import json
import shutil
import tempfile
import unittest
from pathlib import Path

import pandas as pd
import yaml

from candidate_qtl.cli import main


FIXTURES = Path(__file__).parent / "fixtures" / "case_run"

NO_SCORING_WARNING = "rank lane is non-informative: no scoring evidence (positional fallback)"


class CaseRunRankLaneSemanticsTests(unittest.TestCase):
    def _work(self, root: Path, name: str = "case") -> Path:
        work = root / name
        shutil.copytree(FIXTURES, work)
        return work

    def _write_case(self, work: Path, *, without_evidence: bool = False,
                    evidence_rows: str | None = None) -> Path:
        case = yaml.safe_load((work / "case.yaml").read_text(encoding="utf-8"))
        if without_evidence:
            case["inputs"].pop("evidence_resources", None)
        if evidence_rows is not None:
            (work / "evidence.tsv").write_text(
                "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\n"
                + evidence_rows, encoding="utf-8")
        case_path = work / "case_run.yaml"
        case_path.write_text(yaml.safe_dump(case, sort_keys=False), encoding="utf-8")
        return case_path

    def _run_case(self, work: Path, out: Path) -> Path:
        self.assertEqual(0, main(["case-run", "--case", str(work / "case_run.yaml"),
                                  "--resources", str(work / "resources.yaml"), "--out-dir", str(out)]))
        return out

    def _state(self, out: Path) -> dict:
        return json.loads((out / "case_state.json").read_text(encoding="utf-8"))

    def _run_manifest(self, out: Path) -> dict:
        return json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))

    def test_case_run_without_evidence_resources_reports_positional_fallback(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            self._write_case(work, without_evidence=True)
            out = self._run_case(work, root / "out")
            # The run still completes with the same exit-code / status contract.
            self.assertEqual("REPORT_READY", self._state(out)["state"])
            # Every row is a zero-score full tie: the degenerate state the issue
            # describes. It must now be labelled, not silent.
            ranked = pd.read_csv(out / "rank" / "ranked_candidates.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue((pd.to_numeric(ranked["dry_score"]) == 0.0).all())
            self.assertEqual(1, ranked["screening_rank"].nunique())
            self.assertTrue((ranked["dry_group"] == "evidence_insufficient").all())
            for manifest in (self._state(out), self._run_manifest(out)):
                self.assertEqual("not_evaluated_no_scoring_evidence", manifest["rank_semantic_status"])
                self.assertTrue(manifest["rank_is_positional_fallback"])
                self.assertFalse(manifest["rank_metric_eligible"])
                self.assertTrue(manifest["stage_a_is_budget_only"])
                self.assertEqual(NO_SCORING_WARNING, manifest["rank_lane_warning"])
                # Issue #397: the penalty-only declaration is not added to this
                # lane, so the pre-397 manifest stays byte-identical.
                self.assertNotIn("rank_is_penalty_only", manifest)
            # Stage-A is still generated as a budget cut, never called biology.
            self.assertTrue((out / "stage_a" / "deterministic_shortlist.tsv").exists())
            self.assertIn(NO_SCORING_WARNING, (out / "report.md").read_text(encoding="utf-8"))
            # The per-row flags and the run-level state agree.
            self.assertIn("insufficient_independent_evidence",
                          set(ranked["top_eligibility"]))

    def test_case_run_with_evidence_is_evaluable(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            self._write_case(work)
            out = self._run_case(work, root / "out")
            for manifest in (self._state(out), self._run_manifest(out)):
                self.assertEqual("evaluable", manifest["rank_semantic_status"])
                self.assertFalse(manifest["rank_is_positional_fallback"])
                self.assertTrue(manifest["rank_metric_eligible"])
                self.assertFalse(manifest["stage_a_is_budget_only"])
                self.assertEqual("", manifest["rank_lane_warning"])
                self.assertNotIn("rank_is_penalty_only", manifest)

    def test_case_run_zero_score_evidence_is_penalty_only_not_evaluated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            # One included evidence row whose code is not in the weight registry
            # and whose direction is rescue: the evidence lane is wired and the
            # row exists, but there is no support-direction row at all, so the
            # run is penalty-only -- scored as a not_evaluated
            # bucket with rank_metric_eligible false, never read as biology.
            self._write_case(work, evidence_rows="GeneStrong\tannotation_incomplete\trescue\tmedium\tno homology\tannotation\tfixture:gff3\tshort ORF\n")
            out = self._run_case(work, root / "out")
            evidence = pd.read_csv(out / "evidence.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(1, len(evidence))
            for manifest in (self._state(out), self._run_manifest(out)):
                self.assertEqual("not_evaluated_no_supporting_evidence", manifest["rank_semantic_status"])
                self.assertTrue(manifest["rank_is_positional_fallback"])
                self.assertTrue(manifest["rank_is_penalty_only"])
                self.assertFalse(manifest["rank_metric_eligible"])
                self.assertTrue(manifest["stage_a_is_budget_only"])
                self.assertIn("no supporting evidence", manifest["rank_lane_warning"])

    def test_case_run_evidence_provider_failure_degrades_to_review_and_distinct_status(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            # An unreadable evidence file is a provider failure, not an unwired
            # lane: it gets its own status and a surfaced review item while the
            # connection-diagnostic run still completes (exit 0, REPORT_READY).
            (work / "evidence.tsv").write_bytes(b"\xff\xfe\x00\x01\x02garbage-bytes\x00\xff")
            self._write_case(work)
            out = self._run_case(work, root / "out")
            self.assertEqual("REPORT_READY", self._state(out)["state"])
            for manifest in (self._state(out), self._run_manifest(out)):
                self.assertEqual("not_evaluated_evidence_provider_failed", manifest["rank_semantic_status"])
                self.assertTrue(manifest["rank_is_positional_fallback"])
                self.assertFalse(manifest["rank_metric_eligible"])
                self.assertTrue(manifest["stage_a_is_budget_only"])
                self.assertIn("evidence provider failed", manifest["rank_lane_warning"])
                self.assertNotIn("rank_is_penalty_only", manifest)
            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(1, len(reviews[reviews.reason.eq("evidence_provider_failed")]))


class GuidedProfileRankLanePreflightTests(unittest.TestCase):
    def test_guided_profile_plan_detects_empty_evidence_resources(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = root / "case"
            shutil.copytree(FIXTURES, work)
            case = yaml.safe_load((work / "case.yaml").read_text(encoding="utf-8"))
            case["inputs"].pop("evidence_resources", None)
            case_no_evidence = work / "case_no_evidence.yaml"
            case_no_evidence.write_text(yaml.safe_dump(case, sort_keys=False), encoding="utf-8")
            self.assertEqual(0, main(["guided-profile", "--profile", "full-case", "--out-dir", str(root / "out"),
                                      "--case", str(case_no_evidence), "--resources", str(work / "resources.yaml")]))
            manifest = json.loads((root / "out" / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            preflight = manifest["rank_lane_preflight"]
            self.assertTrue(preflight["plan_time_detection"])
            self.assertFalse(preflight["evidence_resources_declared"])
            self.assertEqual("not_evaluated_no_scoring_evidence", preflight["rank_semantic_status"])
            self.assertTrue(preflight["rank_is_positional_fallback"])
            self.assertFalse(preflight["rank_metric_eligible"])
            self.assertTrue(preflight["stage_a_is_budget_only"])
            self.assertEqual(NO_SCORING_WARNING, preflight["warning"])
            # Issue #397: penalty-only is never guessed at plan time (the
            # evidence rows are not known yet), so the preflight stays byte-identical.
            self.assertNotIn("rank_is_penalty_only", preflight)
            # The same plan with evidence declared preflights as evaluable.
            self.assertEqual(0, main(["guided-profile", "--profile", "full-case", "--out-dir", str(root / "out2"),
                                      "--case", str(work / "case.yaml"), "--resources", str(work / "resources.yaml")]))
            manifest = json.loads((root / "out2" / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            preflight = manifest["rank_lane_preflight"]
            self.assertTrue(preflight["evidence_resources_declared"])
            self.assertEqual("evaluable", preflight["rank_semantic_status"])
            self.assertEqual("", preflight["warning"])
            self.assertNotIn("rank_is_penalty_only", preflight)

    def test_guided_profile_execution_replaces_preflight_with_authoritative_state(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = root / "case"
            shutil.copytree(FIXTURES, work)
            case = yaml.safe_load((work / "case.yaml").read_text(encoding="utf-8"))
            case["inputs"].pop("evidence_resources", None)
            case_no_evidence = work / "case_no_evidence.yaml"
            case_no_evidence.write_text(yaml.safe_dump(case, sort_keys=False), encoding="utf-8")
            self.assertEqual(0, main(["guided-profile", "--profile", "full-case", "--out-dir", str(root / "out"),
                                      "--execute", "--case", str(case_no_evidence),
                                      "--resources", str(work / "resources.yaml")]))
            manifest = json.loads((root / "out" / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("completed", manifest["overall_status"])
            preflight = manifest["rank_lane_preflight"]
            # The executed case-run result replaces the plan-time guess with
            # row-level truth; a plan-time guess is never reported as fact.
            self.assertFalse(preflight["plan_time_detection"])
            self.assertEqual("not_evaluated_no_scoring_evidence", preflight["rank_semantic_status"])
            self.assertTrue(preflight["rank_is_positional_fallback"])
            self.assertEqual(NO_SCORING_WARNING, preflight["warning"])
            self.assertNotIn("rank_is_penalty_only", preflight)


if __name__ == "__main__":
    unittest.main()
