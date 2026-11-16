"""Issue #287: case-run chains candidate-elimination into next-experiment.

The existing candidate_elimination() is reused as an optional stage between the
hypothesis ledger and next-experiment; deterministic rank never changes and a
failed elimination degrades to review-required without reusing stale output.
"""

from __future__ import annotations

import json
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import yaml

from candidate_qtl.cli import main
from candidate_qtl.workflow_extensions import case_run


class CaseRunEliminationWiringTests(unittest.TestCase):
    def _work(self, root: Path) -> Path:
        fixture = Path(__file__).parent / "fixtures" / "case_run"
        work = root / "case"
        shutil.copytree(fixture, work)
        candidates = pd.read_csv(work / "candidates.tsv", sep="\t", dtype=str).fillna("")
        candidates["strand"] = "+"
        extra = pd.DataFrame([{
            "gene_id": "GeneExcl", "chrom": "chr1", "start": "800", "end": "980",
            "annotation_quality": "complete_orf", "protein_length": "50", "strand": "+",
        }])
        pd.concat([candidates, extra], ignore_index=True).to_csv(
            work / "candidates.tsv", sep="\t", index=False, lineterminator="\n")
        return work

    def _write(self, work: Path, name: str, content: str) -> Path:
        path = work / name
        path.write_text(content, encoding="utf-8")
        return path

    def _elimination_evidence(self, work: Path, rows: str) -> Path:
        return self._write(work, "elimination_evidence.tsv",
                           "gene_id\thypothesis_id\tevidence_code\tevidence_id\tquality\tcoverage\n" + rows)

    def _run_case(self, work: Path, out: Path, *, inputs: dict | None = None,
                  resources: dict | None = None) -> Path:
        case = yaml.safe_load((work / "case.yaml").read_text(encoding="utf-8"))
        if inputs:
            case["inputs"].update(inputs)
        case_path = work / "case_run.yaml"
        case_path.write_text(yaml.safe_dump(case, sort_keys=False), encoding="utf-8")
        registry = yaml.safe_load((work / "resources.yaml").read_text(encoding="utf-8"))
        if resources:
            registry["resources"].update(resources)
        (work / "resources.yaml").write_text(yaml.safe_dump(registry, sort_keys=False), encoding="utf-8")
        self.assertEqual(0, main(["case-run", "--case", str(case_path),
                                  "--resources", str(work / "resources.yaml"), "--out-dir", str(out)]))
        return out

    def _facts(self, out: Path, experiment: str) -> str:
        options = pd.read_csv(out / "next_experiment" / "next_experiment_options.tsv", sep="\t", dtype=str).fillna("")
        return options[options.experiment_id.eq(experiment)].iloc[0]["input_facts"]

    def test_hard_exclusion_feeds_the_remaining_set_into_next_experiment(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # Elimination run: GeneExcl is hard-excluded, GeneStrong/GeneDark remain.
            work = self._work(root)
            self._elimination_evidence(work, "GeneExcl\t\tinformative_recombinant_excludes\tREC:1\thigh\t1.0\n")
            out = self._run_case(work, root / "out", inputs={"elimination_evidence_resource": "elim_evidence"},
                                 resources={"elim_evidence": {"type": "canonical_evidence", "path": "elimination_evidence.tsv"}})
            remaining = pd.read_csv(out / "candidate_elimination" / "remaining_candidate_set.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"GeneStrong", "GeneDark"}, set(remaining["candidate_gene_id"]))
            state = json.loads((out / "case_state.json").read_text(encoding="utf-8"))
            self.assertEqual(str(out / "candidate_elimination" / "remaining_candidate_set.tsv"),
                             state["remaining_candidate_set"])
            self.assertIsNone(state["integrated_narrowing"])
            self.assertIsNone(state["narrowing_manifest"])
            manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("complete", manifest["steps"]["candidate_elimination"]["status"])
            facts = self._facts(out, "top_n_parental_sanger")
            self.assertIn("candidate_count=2", facts)
            self.assertIn("candidate_count_before_elimination=3", facts)
            self.assertIn("candidate_count_after_elimination=2", facts)

    def test_elimination_never_changes_the_deterministic_rank(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            elim_work = self._work(root / "elim")
            self._elimination_evidence(elim_work, "GeneExcl\t\tinformative_recombinant_excludes\tREC:1\thigh\t1.0\n")
            out = self._run_case(elim_work, root / "out", inputs={"elimination_evidence_resource": "elim_evidence"},
                                 resources={"elim_evidence": {"type": "canonical_evidence", "path": "elimination_evidence.tsv"}})
            baseline = self._run_case(self._work(root / "baseline"), root / "baseline_out")
            self.assertEqual((baseline / "rank" / "ranked_candidates.tsv").read_bytes(),
                             (out / "rank" / "ranked_candidates.tsv").read_bytes())
            self.assertEqual((baseline / "rank" / "score_contributions.tsv").read_bytes(),
                             (out / "rank" / "score_contributions.tsv").read_bytes())

    def test_coding_hypothesis_rejection_keeps_the_gene_in_the_remaining_set(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            self._write(work, "hypotheses.tsv",
                        "hypothesis_id\tcase_id\thypothesis_type\tmember_gene_ids\tproposal_origin\treview_status\tstatus\n"
                        "h_coding\tbundled_fas_like\tsingle_gene_coding\tGeneDark\tmanual\treviewed\tready_to_test\n")
            self._elimination_evidence(work, "GeneDark\th_coding\tfull_cds_no_parental_difference\tSANGER:1\thigh\t1.0\n")
            out = self._run_case(work, root / "out",
                                 inputs={"elimination_evidence_resource": "elim_evidence",
                                         "hypotheses_resource": "hypotheses"},
                                 resources={"elim_evidence": {"type": "canonical_evidence", "path": "elimination_evidence.tsv"},
                                            "hypotheses": {"type": "hypotheses", "path": "hypotheses.tsv"}})
            remaining = pd.read_csv(out / "candidate_elimination" / "remaining_candidate_set.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("GeneDark", set(remaining["candidate_gene_id"]))
            elimination_manifest = json.loads((out / "candidate_elimination" / "elimination_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(1, elimination_manifest["mechanism_only_rejected_count"])
            self.assertEqual(3, elimination_manifest["remaining_candidates"])

    def test_without_elimination_input_the_stage_is_not_run(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            out = self._run_case(work, root / "out")
            self.assertFalse((out / "candidate_elimination").exists())
            state = json.loads((out / "case_state.json").read_text(encoding="utf-8"))
            self.assertIsNone(state["remaining_candidate_set"])
            self.assertIsNone(state["integrated_narrowing"])
            self.assertIsNone(state["narrowing_manifest"])
            facts = self._facts(out, "top_n_parental_sanger")
            self.assertIn("candidate_count=3", facts)
            self.assertNotIn("candidate_count_before_elimination", facts)
            self.assertNotIn("candidate_count_after_elimination", facts)

    def test_elimination_failure_degrades_to_review_and_never_forwards_stale_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            # Missing evidence_code column: candidate_elimination() fails closed.
            self._write(work, "elimination_evidence.tsv", "gene_id\nGeneStrong\n")
            out = self._run_case(work, root / "out", inputs={"elimination_evidence_resource": "elim_evidence"},
                                 resources={"elim_evidence": {"type": "canonical_evidence", "path": "elimination_evidence.tsv"}})
            self.assertFalse((out / "candidate_elimination" / "remaining_candidate_set.tsv").exists())
            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            failed = reviews[reviews.reason.eq("candidate_elimination_failed")]
            self.assertEqual(1, len(failed))
            self.assertIn("stale remaining sets are never reused", failed.iloc[0]["required_next_action"])
            # next-experiment still runs on the ranked input alone.
            facts = self._facts(out, "top_n_parental_sanger")
            self.assertIn("candidate_count=3", facts)
            self.assertNotIn("candidate_count_after_elimination", facts)

    def test_a_later_failure_never_forwards_an_earlier_remaining_set(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            evidence = self._elimination_evidence(work, "GeneExcl\t\tinformative_recombinant_excludes\tREC:1\thigh\t1.0\n")
            inputs = {"elimination_evidence_resource": "elim_evidence"}
            resources = {"elim_evidence": {"type": "canonical_evidence", "path": "elimination_evidence.tsv"}}
            # First run succeeds: a remaining set exists on disk.
            out = self._run_case(work, root / "out", inputs=inputs, resources=resources)
            self.assertTrue((out / "candidate_elimination" / "remaining_candidate_set.tsv").exists())
            self.assertIn("candidate_count_after_elimination=2", self._facts(out, "top_n_parental_sanger"))
            # Second run into the SAME out-dir fails: the stale remaining set must
            # not reach next-experiment, and the failed stage lands in review.
            evidence.write_text("gene_id\nGeneStrong\n", encoding="utf-8")
            out = self._run_case(work, root / "out", inputs=inputs, resources=resources)
            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(1, len(reviews[reviews.reason.eq("candidate_elimination_failed")]))
            facts = self._facts(out, "top_n_parental_sanger")
            self.assertIn("candidate_count=3", facts)
            self.assertNotIn("candidate_count_after_elimination", facts)


if __name__ == "__main__":
    unittest.main()


class CaseRunIntegratedNarrowingWiringTests(unittest.TestCase):
    """Issue #303: case-run chains the optional integrated-narrowing stage.

    The stage runs after candidate-elimination and only when a family artifact
    (family_to_copy_input.tsv / family_copy_discrimination.tsv) is registered;
    it never changes the deterministic rank, its digest covers every input, a
    failure degrades to review-required, and stale narrowing outputs are never
    forwarded.
    """

    # The case-run harness helpers (_work / _write / _run_case / _facts) are
    # defined on CaseRunEliminationWiringTests; reuse them unbound so both
    # classes drive the same fixture without duplicating the harness.
    _work = CaseRunEliminationWiringTests._work
    _write = CaseRunEliminationWiringTests._write
    _run_case = CaseRunEliminationWiringTests._run_case
    _facts = CaseRunEliminationWiringTests._facts

    FAMILY_RECOVERY = (
        "hypothesis_id\tgene_family\tmember_interval_genes\n"
        "cmp_001\tFAS_LIKE\tGeneStrong\n"
    )
    CONVERGENCE = (
        "comparison_id\tfamily_id\tconvergence_status\tsource_dependence_status\tfour_state\n"
        "cmp_001\tFAS_LIKE\ta_b_convergent\tsource_independent\ta_and_b_supported\n"
    )
    COPY_DISCRIMINATION = (
        "family_id\tcandidate_gene_id\tdisambiguation_status\trank_contribution\n"
        "FAS_LIKE\tGeneStrong\tpreferred_copy\tnone\n"
    )
    COPY_SHORTLIST = (
        "family_id\tcandidate_gene_id\tdisambiguation_status\tfirst_wet_priority\trank_contribution\n"
        "FAS_LIKE\tGeneStrong\tpreferred_copy\tfirst\tnone\n"
    )

    def _narrowing_inputs(self, work: Path, *, bad_copy: bool = False) -> tuple[dict, dict]:
        self._write(work, "family_recovery.tsv", self.FAMILY_RECOVERY)
        self._write(work, "convergence.tsv", self.CONVERGENCE)
        if bad_copy:
            self._write(work, "family_copy.tsv", "family_id\tcandidate_gene_id\nFAS_LIKE\tGeneStrong\n")
        else:
            self._write(work, "family_copy.tsv", self.COPY_DISCRIMINATION)
        self._write(work, "family_copy_shortlist.tsv", self.COPY_SHORTLIST)
        inputs = {
            "narrowing_family_recovery_resource": "family_recovery",
            "narrowing_convergence_resource": "convergence",
            "narrowing_family_copy_resource": "family_copy",
            "narrowing_family_copy_shortlist_resource": "family_copy_shortlist",
        }
        resources = {
            "family_recovery": {"type": "hypotheses", "path": "family_recovery.tsv"},
            "convergence": {"type": "comparison", "path": "convergence.tsv"},
            "family_copy": {"type": "family_copy_discrimination", "path": "family_copy.tsv"},
            "family_copy_shortlist": {"type": "family_copy_shortlist", "path": "family_copy_shortlist.tsv"},
        }
        return inputs, resources

    def test_narrowing_stage_runs_after_elimination_with_family_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            inputs, resources = self._narrowing_inputs(work)
            out = self._run_case(work, root / "out", inputs=inputs, resources=resources)
            narrowed = pd.read_csv(out / "integrated_narrowing" / "integrated_narrowing.tsv",
                                   sep="\t", dtype=str).fillna("")
            row = narrowed[narrowed.gene_id.eq("GeneStrong")].iloc[0]
            self.assertEqual("N1", row["narrowing_tier"])
            self.assertEqual("1", row["narrowing_priority"])
            self.assertEqual("none", row["rank_contribution"])
            manifest = json.loads((out / "integrated_narrowing" / "narrowing_manifest.json").read_text(encoding="utf-8"))
            state = json.loads((out / "case_state.json").read_text(encoding="utf-8"))
            self.assertIsNone(state["remaining_candidate_set"])
            self.assertEqual(str(out / "integrated_narrowing" / "integrated_narrowing.tsv"),
                             state["integrated_narrowing"])
            self.assertEqual(str(out / "integrated_narrowing" / "narrowing_manifest.json"),
                             state["narrowing_manifest"])
            # GeneStrong -> N1 (preferred copy); GeneDark -> N3 (base top/reserve/
            # rescue, no family info); GeneExcl -> N4 (base rank only).
            self.assertEqual({"N1": 1, "N2": 0, "N3": 1, "N4": 1}, manifest["tier_counts"])
            run_manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("complete", run_manifest["steps"]["integrated_narrowing"]["status"])
            # The narrowing stage must not change the deterministic rank.
            self.assertIn("integrated_narrowing", str(run_manifest["steps"]))
            self.assertTrue((out / "stage_a" / "deterministic_shortlist.tsv").exists())
            narrowed_stage_a = json.loads(
                (out / "stage_a_narrowed" / "shortlist_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(20, narrowed_stage_a["target_size"])
            self.assertEqual("A_narrowed_wet_priority", narrowed_stage_a["stage"])
            self.assertEqual("narrowing_priority", narrowed_stage_a["priority_column"])
            self.assertIsNone(narrowed_stage_a["boundary_screening_rank"])
            self.assertIsNotNone(narrowed_stage_a["boundary_priority"])
            self.assertTrue(narrowed_stage_a["selection_policy"].startswith("narrowing_priority asc"))
            facts = self._facts(out, "top_n_parental_sanger")
            self.assertIn("integrated_narrowing_status=evaluated", facts)
            self.assertIn("narrowing_wet_priority_count=1", facts)

    def test_narrowing_never_changes_the_deterministic_rank(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root / "narrow")
            inputs, resources = self._narrowing_inputs(work)
            out = self._run_case(work, root / "out", inputs=inputs, resources=resources)
            baseline = self._run_case(self._work(root / "baseline"), root / "baseline_out")
            self.assertEqual((baseline / "rank" / "ranked_candidates.tsv").read_bytes(),
                             (out / "rank" / "ranked_candidates.tsv").read_bytes())
            self.assertEqual((baseline / "rank" / "score_contributions.tsv").read_bytes(),
                             (out / "rank" / "score_contributions.tsv").read_bytes())

    def test_narrowing_stage_is_not_run_without_family_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            out = self._run_case(work, root / "out")
            self.assertFalse((out / "integrated_narrowing").exists())
            run_manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertNotIn("integrated_narrowing", run_manifest["steps"])
            provider = json.loads((out / "narrowing_provider_state.json").read_text(encoding="utf-8"))
            self.assertEqual("not_evaluated", provider["status"])
            self.assertEqual("none", provider["rank_contribution"])
            self.assertTrue((out / "stage_a" / "deterministic_shortlist.tsv").exists())
            stage_a = json.loads((out / "stage_a" / "shortlist_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(20, stage_a["target_size"])

    def test_triple_route_bridge_rejects_a_different_candidate_universe(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            self._write(work, "family_recovery.tsv", self.FAMILY_RECOVERY)
            self._write(work, "convergence.tsv", self.CONVERGENCE)
            self._write(work, "family_copy.tsv", self.COPY_DISCRIMINATION)
            self._write(work, "family_copy_shortlist.tsv", self.COPY_SHORTLIST)
            bridge_candidates = self._write(
                work, "bridge_candidates.tsv", "gene_id\nGeneStrong\n"
            )
            out = root / "out"
            result = case_run(
                str(work / "case.yaml"), str(work / "resources.yaml"), str(out),
                narrowing_overrides={
                    "bridge_case_id": "bundled_fas_like",
                    "bridge_candidates": str(bridge_candidates),
                    "family_recovery": str(work / "family_recovery.tsv"),
                    "convergence": str(work / "convergence.tsv"),
                    "family_copy": str(work / "family_copy.tsv"),
                    "family_copy_shortlist": str(work / "family_copy_shortlist.tsv"),
                },
            )
            self.assertIsNone(result["integrated_narrowing"])
            provider = json.loads((out / "narrowing_provider_state.json").read_text(encoding="utf-8"))
            self.assertEqual("triple_route_case_identity_mismatch", provider["reason"])
            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(1, len(reviews[reviews.reason.eq("triple_route_case_identity_mismatch")]))

    def test_failed_stage_a_rerun_does_not_forward_a_stale_shortlist(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            out = self._run_case(work, root / "out")
            stale = out / "stage_a" / "deterministic_shortlist.tsv"
            self.assertTrue(stale.exists())
            with patch("candidate_qtl.workflow_extensions.screening_shortlist", side_effect=ValueError("forced")):
                result = case_run(str(work / "case_run.yaml"), str(work / "resources.yaml"), str(out))
            self.assertIsNone(result["stage_a_shortlist"])
            self.assertTrue(stale.exists())
            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(1, len(reviews[reviews.reason.eq("stage_a_shortlist_failed")]))

    def test_narrowing_failure_degrades_to_review_and_never_forwards_stale(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            inputs, resources = self._narrowing_inputs(work, bad_copy=True)
            out = self._run_case(work, root / "out", inputs=inputs, resources=resources)
            # disambiguation_status column missing: integrated_narrowing() fails closed.
            self.assertFalse((out / "integrated_narrowing" / "integrated_narrowing.tsv").exists())
            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            failed = reviews[reviews.reason.eq("integrated_narrowing_failed")]
            self.assertEqual(1, len(failed))
            self.assertIn("stale narrowing outputs are never reused", failed.iloc[0]["required_next_action"])
            # next-experiment still runs on the ranked input alone.
            self.assertTrue((out / "next_experiment" / "next_experiment_decision.json").exists())

    def test_narrowing_digest_invalidates_when_inputs_change(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            inputs, resources = self._narrowing_inputs(work)
            out = self._run_case(work, root / "out", inputs=inputs, resources=resources)
            first = (out / "integrated_narrowing" / "integrated_narrowing.tsv").read_bytes()
            # Change the convergence state and the copy status: the digest must
            # catch both and the stage must recompute instead of serving the
            # stale cache (b_only + unresolved copy -> N2, not N1).
            convergence_path = work / "convergence.tsv"
            convergence_path.write_text(
                "comparison_id\tfamily_id\tconvergence_status\tsource_dependence_status\tfour_state\n"
                "cmp_001\tFAS_LIKE\tb_only_supported\tnot_applicable\tb_only_supported\n",
                encoding="utf-8")
            (work / "family_copy.tsv").write_text(
                "family_id\tcandidate_gene_id\tdisambiguation_status\trank_contribution\n"
                "FAS_LIKE\tGeneStrong\tfamily_supported_copy_unresolved\tnone\n",
                encoding="utf-8")
            out = self._run_case(work, root / "out", inputs=inputs, resources=resources)
            second = (out / "integrated_narrowing" / "integrated_narrowing.tsv").read_bytes()
            self.assertNotEqual(first, second)
            narrowed = pd.read_csv(out / "integrated_narrowing" / "integrated_narrowing.tsv",
                                   sep="\t", dtype=str).fillna("")
            row = narrowed[narrowed.gene_id.eq("GeneStrong")].iloc[0]
            self.assertEqual("b_only", row["family_recovery_state"])
            self.assertEqual("N2", row["narrowing_tier"])
            self.assertEqual("validated_family_copy_unresolved", row["narrowing_reason"])
