"""Issue #397: a run whose evidence rows carry no support-direction row is a
penalty-only run and is declared as such in the case-run / guided-profile
manifests instead of being read as evaluated biology.

Three states are fixed here:
  1. evidence rows exist but ``direction == "support"`` count is 0
     (context/against/rescue only) -> penalty-only: rank_is_penalty_only=true,
     rank_semantic_status=not_evaluated_no_supporting_evidence,
     rank_metric_eligible=false. This holds even when the against rows score,
     so the lane is NOT a positional fallback: penalty-only is orthogonal.
  2. at least one support row -> evaluable exactly as before.
  3. zero evidence rows -> the pre-397 not_evaluated_no_scoring_evidence path,
     byte-identical (no rank_is_penalty_only key added).

The benchmark gate test fixes that the rank_metric_eligible column is
authoritative in the gene-rank lane: a scored hit/miss whose rank lane declared
no supporting evidence leaves the denominator like any not_evaluated case.
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
from candidate_qtl import schemas

FIXTURES = Path(__file__).parent / "fixtures" / "case_run"


class CaseRunPenaltyOnlyStateTests(unittest.TestCase):
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

    def test_support_zero_but_scored_against_rows_is_penalty_only_not_positional_fallback(self) -> None:
        """Issue #397 core case: evidence rows exist and against rows score, but
        support count is 0. The run is penalty-only (rank_metric_eligible=false,
        stage-a budget-only) and must NOT be claimed as a positional fallback,
        because the rank really was scored -- only by penalties."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            self._write_case(work, evidence_rows="\n".join([
                "GeneStrong\thomology_context\tcontext\tmedium\tno homology\tannotation\tfixture:gff3\tconserved",
                "GeneStrong\tparalog_ambiguous\tagainst\tmedium\tno homology\tcomparative\tfixture:paralog\tnear-tied paralog",
                "GeneDark\tannotation_incomplete\trescue\tmedium\tno homology\tannotation\tfixture:gff3\tshort ORF",
            ]) + "\n")
            out = self._run_case(work, root / "out")
            evidence = pd.read_csv(out / "evidence.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(3, len(evidence))
            self.assertEqual(0, int((evidence["direction"] == "support").sum()))
            # The against row really scores: the lane is scored (not a
            # positional fallback) yet still penalty-only.
            contributions = pd.read_csv(out / "rank" / "score_contributions.tsv", sep="\t", dtype=str).fillna("")
            self.assertGreater(len(contributions), 0)
            for manifest in (self._state(out), self._run_manifest(out)):
                self.assertTrue(manifest["rank_is_penalty_only"])
                self.assertEqual("not_evaluated_no_supporting_evidence", manifest["rank_semantic_status"])
                self.assertFalse(manifest["rank_metric_eligible"])
                self.assertFalse(manifest["rank_is_positional_fallback"])
                # Design decision #2: penalty-only is budget-only even though the
                # scored contributions exist.
                self.assertTrue(manifest["stage_a_is_budget_only"])
                self.assertIn("no supporting evidence", manifest["rank_lane_warning"])

    def test_one_support_row_is_evaluable_even_with_context_rows(self) -> None:
        """A single support row (regardless of which gene it is on) keeps the
        lane evaluable and rank-metric-eligible."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            self._write_case(work, evidence_rows="\n".join([
                "GeneStrong\thomology_context\tcontext\tmedium\tno homology\tannotation\tfixture:gff3\tconserved",
                "GeneDark\ttrait_function_match\tsupport\tmedium\tstarch trait term\tcomparative\tfixture:term\ttrait-term match",
            ]) + "\n")
            out = self._run_case(work, root / "out")
            evidence = pd.read_csv(out / "evidence.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(1, int((evidence["direction"] == "support").sum()))
            for manifest in (self._state(out), self._run_manifest(out)):
                self.assertEqual("evaluable", manifest["rank_semantic_status"])
                self.assertTrue(manifest["rank_metric_eligible"])
                self.assertNotIn("rank_is_penalty_only", manifest)
                self.assertFalse(manifest["stage_a_is_budget_only"])

    def test_zero_evidence_rows_stays_no_scoring_evidence_without_penalty_key(self) -> None:
        """The pre-397 no-scoring-evidence path is preserved byte-for-byte: no
        rank_is_penalty_only key is added and the status keeps its old value."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            self._write_case(work, without_evidence=True)
            out = self._run_case(work, root / "out")
            for manifest in (self._state(out), self._run_manifest(out)):
                self.assertEqual("not_evaluated_no_scoring_evidence", manifest["rank_semantic_status"])
                self.assertNotIn("rank_is_penalty_only", manifest)
                self.assertFalse(manifest["rank_metric_eligible"])
                self.assertTrue(manifest["stage_a_is_budget_only"])


class GuidedProfilePenaltyOnlyPropagationTests(unittest.TestCase):
    def test_executed_profile_propagates_penalty_only_state(self) -> None:
        """The executed case-run's penalty-only state replaces the plan-time
        preflight guess, exactly like the other rank-lane states (#346)."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = root / "case"
            shutil.copytree(FIXTURES, work)
            (work / "evidence.tsv").write_text(
                "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\n"
                "GeneStrong\tparalog_ambiguous\tagainst\tmedium\tno homology\tcomparative\tfixture:paralog\tnear-tied paralog\n",
                encoding="utf-8")
            self.assertEqual(0, main(["guided-profile", "--profile", "full-case", "--out-dir", str(root / "out"),
                                      "--execute", "--case", str(work / "case.yaml"),
                                      "--resources", str(work / "resources.yaml")]))
            manifest = json.loads((root / "out" / "guided_profile_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("completed", manifest["overall_status"])
            preflight = manifest["rank_lane_preflight"]
            self.assertFalse(preflight["plan_time_detection"])
            self.assertTrue(preflight["rank_is_penalty_only"])
            self.assertEqual("not_evaluated_no_supporting_evidence", preflight["rank_semantic_status"])
            self.assertFalse(preflight["rank_metric_eligible"])
            self.assertTrue(preflight["stage_a_is_budget_only"])
            self.assertIn("no supporting evidence", preflight["warning"])


@unittest.skip("public release: supervised-benchmark commands are not shipped")
class BenchmarkGeneRankEligibilityGateTests(unittest.TestCase):
    """The rank_metric_eligible column is authoritative in the gene-rank lane
    (issue #397): a scored hit/miss whose lane declared no supporting evidence
    is not_evaluated and leaves the hit+miss denominator."""

    def _build(self, root: Path) -> tuple[Path, Path, Path, Path, Path]:
        truth = root / "truth.tsv"
        truth.write_text("\n".join([
            "case_id\ttruth_unit_type\ttruth_gene_id\ttruth_family_id",
            "case_g001\tgene\tG_1\tfam_1",
            "case_g002\tgene\tG_2\tfam_2",
            "case_g003\tgene\tG_3\tfam_3",
        ]) + "\n", encoding="utf-8")

        chain = root / "chain_cases.tsv"
        chain.write_text("\n".join([
            "case_id\ta_only_hit\ta_lane_state\tb0_only_hit\tb0_lane_state\tb1_only_hit\tb1_lane_state\tunion_hit\tunion_lane_state\tunion_b_lane",
            "case_g001\tTrue\thit\tTrue\thit\t\tnot_scoreable\tTrue\thit\tb0",
            "case_g002\tTrue\thit\tTrue\thit\t\tnot_scoreable\tTrue\thit\tb0",
            "case_g003\tFalse\tmiss\tFalse\tmiss\t\tnot_scoreable\tFalse\tmiss\tb0",
        ]) + "\n", encoding="utf-8")

        audit = root / "case_states.tsv"
        audit.write_text("\n".join([
            "case_id\tcase_state\tcase_state_reason",
            "case_g001\tclean\tno_red_or_unresolved_source",
            "case_g002\tclean\tno_red_or_unresolved_source",
            "case_g003\tclean\tno_red_or_unresolved_source",
        ]) + "\n", encoding="utf-8")

        # case_g002 hits but its rank lane declared no supporting evidence
        # (rank_metric_eligible=false, penalty-only); case_g003 misses with the
        # same eligibility. Only case_g001 stays eligible.
        score = root / "score_cases.tsv"
        score.write_text("\n".join([
            "case_id\ttruth_unit_type\tcleanliness_state\tcase_state\tstate_reason\tfamily_scoreable\trank_metric_eligible",
            "case_g001\tgene\tclean\thit\trecovered_sole_top\tTrue\tTrue",
            "case_g002\tgene\tclean\thit\trecovered_sole_top\tTrue\tFalse",
            "case_g003\tgene\tclean\tmiss\tranking_miss\tTrue\tFalse",
        ]) + "\n", encoding="utf-8")

        c_cases = root / "c_cases.tsv"
        c_cases.write_text("\n".join([
            "case_id\tfamily_recall_20",
            "case_g001\tTrue",
            "case_g002\tTrue",
            "case_g003\tFalse",
        ]) + "\n", encoding="utf-8")
        return truth, chain, audit, score, c_cases

    def test_penalty_only_hit_and_miss_leave_the_gene_rank_denominator(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            truth, chain, audit, score, c_cases = self._build(root)
            out = root / "metrics"
            self.assertEqual(0, main(["supervised-benchmark-route-family-metrics",
                                      "--chain-cases", str(chain), "--truth", str(truth),
                                      "--case-audit", str(audit), "--score-cases", str(score),
                                      "--c-cases", str(c_cases), "--out-dir", str(out)]))
            manifest_out = json.loads((out / "route_family_metrics_manifest.json").read_text(encoding="utf-8"))
            gene_rank = manifest_out["metrics"]["gene_rank_recall"]
            # Only case_g001 contributes: hit=1, miss=0, the penalty-only
            # hit (g002) and miss (g003) are not_evaluated, never converted
            # into misses.
            self.assertEqual(1, gene_rank["denominator"])
            self.assertEqual(1, gene_rank["hit"])
            self.assertEqual(0, gene_rank["miss"])
            self.assertEqual(2, gene_rank["not_evaluated"])
            self.assertEqual(0, gene_rank["not_scoreable"])
            self.assertEqual(1.0, gene_rank["recall"])

            cases = pd.read_csv(out / "route_family_metrics_cases.tsv", sep="\t", dtype=str).fillna("")
            for case_id in ("case_g002", "case_g003"):
                row = cases[(cases["case_id"] == case_id) & (cases["metric"] == "gene_rank_recall")].iloc[0]
                self.assertEqual("False", row["rank_metric_eligible"])
                self.assertEqual("not_evaluated", row["state"])
                self.assertEqual("gene_rank_metric_not_eligible", row["state_reason"])
                self.assertEqual("False", row["miss"])
                self.assertEqual("False", row["hit"])


class RankSemanticStatusVocabularyTests(unittest.TestCase):
    def test_vocabulary_registers_three_existing_plus_penalty_only(self) -> None:
        self.assertEqual(
            frozenset({
                "evaluable",
                "not_evaluated_no_scoring_evidence",
                "not_evaluated_evidence_provider_failed",
                "not_evaluated_no_supporting_evidence",
            }),
            schemas.RANK_SEMANTIC_STATUSES)

    def test_gate_accepts_registered_values_and_rejects_unknown(self) -> None:
        for status in schemas.RANK_SEMANTIC_STATUSES:
            self.assertEqual(status, schemas.require_rank_semantic_status(status))
        with self.assertRaises(ValueError):
            schemas.require_rank_semantic_status("penalty_only")
        with self.assertRaises(ValueError):
            schemas.require_rank_semantic_status("evaluable_but_penalty_only")


if __name__ == "__main__":
    unittest.main()
