from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import build_parser
from candidate_qtl.next_experiment import select_next_experiment


class NextExperimentTests(unittest.TestCase):
    def _ranked(self, root: Path, count: int = 12) -> Path:
        path = root / "ranked.tsv"
        pd.DataFrame({
            "gene_id": [f"Gene{index:02d}" for index in range(count)],
            "dry_group": ["top" if index < 6 else "reserve" for index in range(count)],
            "dry_score": [str(count - index) for index in range(count)],
        }).to_csv(path, sep="\t", index=False, lineterminator="\n")
        return path

    def _manifest(self, root: Path, **overrides) -> Path:
        support = {
            "budget_class": "high",
            "allowed_technologies": [
                "additional_grid_markers", "top_n_parental_sanger",
                "multiplex_amplicon_sequencing", "target_capture",
                "sv_te_cnv_pav_branch", "expression_assay",
                "annotation_rescue", "non_gene_follow_up",
            ],
            "available_samples": {
                "parents": "yes", "informative_recombinants": "yes",
                "segregating_population": "yes", "target_tissue": "yes",
                "reference_proteome": "yes", "comparative_annotation_resources": "yes",
                "association_cohort": "yes", "population_phenotypes": "yes",
            },
            "local_cost_classes": {
                "additional_grid_markers": "low", "top_n_parental_sanger": "low",
                "multiplex_amplicon_sequencing": "medium", "target_capture": "high",
                "sv_te_cnv_pav_branch": "high", "expression_assay": "medium",
                "annotation_rescue": "low", "non_gene_follow_up": "medium",
            },
        }
        support.update(overrides)
        path = root / "case.yaml"
        path.write_text("case_id: test\ndecision_support:\n" + "\n".join(self._yaml_lines(support, 2)) + "\n", encoding="utf-8")
        return path

    def _yaml_lines(self, value, indent: int):
        lines = []
        prefix = " " * indent
        for key, item in value.items():
            if isinstance(item, dict):
                lines.append(f"{prefix}{key}:")
                lines.extend(self._yaml_lines(item, indent + 2))
            elif isinstance(item, list):
                lines.append(f"{prefix}{key}:")
                lines.extend(f"{' ' * (indent + 2)}- {entry}" for entry in item)
            else:
                lines.append(f"{prefix}{key}: {item}")
        return lines

    def test_marker_first_uses_observed_recombinant_counts_and_preserves_rank(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            ranked = self._ranked(root)
            original = ranked.read_bytes()
            case = self._manifest(root)
            mapping = root / "mapping.json"
            mapping.write_text(json.dumps({"mapping_status": "insufficient_marker_density"}), encoding="utf-8")
            recombinants = root / "recombinants.tsv"
            recombinants.write_text("individual_id\tcandidates_before\tcandidates_after\nR1\t12\t4\n", encoding="utf-8")
            result = select_next_experiment(
                str(ranked), str(case), str(root / "out"),
                mapping_gate_arg=str(mapping), informative_recombinants_arg=str(recombinants),
            )
            options = pd.read_csv(result["options"], sep="\t", dtype=str).fillna("")
            marker = options[options.experiment_id.eq("additional_grid_markers")].iloc[0]
            self.assertEqual("additional_grid_markers", result["primary_recommendation"])
            self.assertEqual("8", marker.candidates_potentially_eliminated)
            self.assertEqual("feasible", marker.feasibility_status)
            self.assertEqual(original, ranked.read_bytes())

    def test_sv_first_when_tandem_and_partial_amplification_signals_exist(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            ranked = self._ranked(root)
            case = self._manifest(root)
            flags = root / "flags.tsv"
            flags.write_text("gene_id\ttandem_cluster\tnotes\nGene00\tcluster1\tparalog repeat region\n", encoding="utf-8")
            outcomes = root / "outcomes.tsv"
            outcomes.write_text(
                "experiment_id\tstate\toutcome\tfailure_reason\n"
                "top_n_parental_sanger\tfailed\ttechnical_failure\tparent non-amplification\n",
                encoding="utf-8",
            )
            result = select_next_experiment(
                str(ranked), str(case), str(root / "out"),
                candidate_flags_arg=str(flags), assay_outcomes_arg=str(outcomes),
            )
            options = pd.read_csv(result["options"], sep="\t", dtype=str).fillna("")
            sanger = options[options.experiment_id.eq("top_n_parental_sanger")].iloc[0]
            self.assertEqual("sv_te_cnv_pav_branch", result["primary_recommendation"])
            self.assertEqual("blocked", sanger.feasibility_status)
            self.assertIn("R_TECHNICAL_FAILURE_IS_NOT_NEGATIVE", sanger.rules_fired)

    def test_unknown_cost_is_not_imputed_and_outputs_are_byte_stable(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            ranked = self._ranked(root, 3)
            case = self._manifest(root, budget_class="unknown", local_cost_classes={})
            out = root / "out"
            result = select_next_experiment(str(ranked), str(case), str(out))
            paths = [Path(result[key]) for key in ("options", "decision", "dependencies", "review_required")]
            first = [path.read_bytes() for path in paths]
            select_next_experiment(str(ranked), str(case), str(out))
            self.assertEqual(first, [path.read_bytes() for path in paths])
            options = pd.read_csv(result["options"], sep="\t", dtype=str).fillna("")
            self.assertTrue(options.estimated_cost_class.eq("unknown").all())
            decision = json.loads(Path(result["decision"]).read_text(encoding="utf-8"))
            self.assertFalse(decision["unknown_values_imputed"])
            self.assertEqual("pending", decision["human_approval"])

    def test_partial_and_complete_informative_negative_are_audited_separately(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            ranked = self._ranked(root, 4)
            case = self._manifest(root)
            outcomes = root / "outcomes.tsv"
            outcomes.write_text(
                "experiment_id\tstate\toutcome\tfailure_reason\n"
                "multiplex_amplicon_sequencing\tpartial\tinconclusive\t\n"
                "expression_assay\tcomplete\tinformative_negative\t\n",
                encoding="utf-8",
            )
            result = select_next_experiment(
                str(ranked), str(case), str(root / "out"), assay_outcomes_arg=str(outcomes)
            )
            options = pd.read_csv(result["options"], sep="\t", dtype=str).fillna("").set_index("experiment_id")
            self.assertIn("R_PARTIAL_COVERAGE_REQUIRES_COMPLETION_OR_ALTERNATIVE", options.loc["multiplex_amplicon_sequencing", "rules_fired"])
            self.assertIn("R_HISTORY_COMPLETE", options.loc["expression_assay", "rules_fired"])
            self.assertIn("R_INFORMATIVE_NEGATIVE_ONLY_ADDRESSES_MEASURED_SCOPE", options.loc["expression_assay", "rules_fired"])

    def test_budget_class_blocks_only_more_expensive_supplied_local_classes(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            ranked = self._ranked(root, 4)
            case = self._manifest(root, budget_class="low")
            result = select_next_experiment(str(ranked), str(case), str(root / "out"))
            options = pd.read_csv(result["options"], sep="\t", dtype=str).fillna("").set_index("experiment_id")
            self.assertEqual("feasible", options.loc["top_n_parental_sanger", "feasibility_status"])
            self.assertEqual("blocked", options.loc["target_capture", "feasibility_status"])
            self.assertEqual("high", options.loc["target_capture", "estimated_cost_class"])

    def test_normalized_wet_state_passes_incomplete_scope_to_sanger_decision(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            ranked = self._ranked(root, 4)
            case = self._manifest(root)
            wet = root / "candidate_wet_state.tsv"
            wet.write_text(
                "gene_id\tcandidate_state\tregulatory_sv_state\n"
                "Gene00\tpartially_covered\tnot_tested\n"
                "Gene01\tcoding_region_cleared\tnot_tested_retained\n",
                encoding="utf-8",
            )
            result = select_next_experiment(
                str(ranked), str(case), str(root / "out"), wet_candidate_state_arg=str(wet)
            )
            options = pd.read_csv(result["options"], sep="\t", dtype=str).fillna("").set_index("experiment_id")
            sanger = options.loc["top_n_parental_sanger"]
            self.assertIn("R_WET_STATE_INCOMPLETE_SCOPE_REQUIRES_COMPLETION", sanger.rules_fired)
            self.assertIn("R_WET_CODING_NEGATIVE_IS_SCOPE_LIMITED", sanger.rules_fired)
            self.assertIn("wet_candidates_incomplete=1", sanger.input_facts)

    def test_llm_hypothesis_cannot_change_action_priority_but_reviewed_manual_sv_can(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            ranked = self._ranked(root, 4)
            case = self._manifest(root)
            hypotheses = root / "hypotheses.tsv"
            hypotheses.write_text(
                "hypothesis_id\thypothesis_type\tproposal_origin\treview_status\tstatus\trequired_experiment_ids\n"
                "H1\tsv_cnv_pav\tllm\tpending_review\tpending_review\tsv_te_cnv_pav_branch\n",
                encoding="utf-8",
            )
            llm_result = select_next_experiment(
                str(ranked), str(case), str(root / "llm"), hypotheses_arg=str(hypotheses)
            )
            self.assertNotEqual("sv_te_cnv_pav_branch", llm_result["primary_recommendation"])
            hypotheses.write_text(
                "hypothesis_id\thypothesis_type\tproposal_origin\treview_status\tstatus\trequired_experiment_ids\n"
                "H1\tsv_cnv_pav\tmanual\treviewed\tready_to_test\tsv_te_cnv_pav_branch\n",
                encoding="utf-8",
            )
            manual_result = select_next_experiment(
                str(ranked), str(case), str(root / "manual"), hypotheses_arg=str(hypotheses)
            )
            self.assertEqual("sv_te_cnv_pav_branch", manual_result["primary_recommendation"])

    def test_hypothesis_action_context_requires_fail_closed_origin_schema(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            hypotheses = root / "hypotheses.tsv"
            hypotheses.write_text("hypothesis_type\tstatus\nsv_cnv_pav\tready_to_test\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "fail-closed"):
                select_next_experiment(
                    str(self._ranked(root, 4)), str(self._manifest(root)), str(root / "out"),
                    hypotheses_arg=str(hypotheses),
                )

    def test_cli_exposes_separate_action_rank_command(self) -> None:
        args = build_parser().parse_args([
            "next-experiment", "--ranked", "ranked.tsv", "--case", "case.yaml", "--out-dir", "out",
        ])
        self.assertEqual("next-experiment", args.command)

    def test_annotation_rescue_is_ranked_first_on_reviewed_manual_annotation_gap(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            hypotheses = root / "hypotheses.tsv"
            hypotheses.write_text(
                "hypothesis_id\thypothesis_type\tproposal_origin\treview_status\tstatus\trequired_experiment_ids\n"
                "H1\tannotation_missing\tmanual\treviewed\tready_to_test\t\n",
                encoding="utf-8",
            )
            result = select_next_experiment(
                str(self._ranked(root, 4)), str(self._manifest(root)), str(root / "out"),
                hypotheses_arg=str(hypotheses),
            )
            options = pd.read_csv(result["options"], sep="\t", dtype=str).fillna("").set_index("experiment_id")
            rescue = options.loc["annotation_rescue"]
            self.assertEqual("annotation_rescue", result["primary_recommendation"])
            self.assertEqual("feasible", rescue.feasibility_status)
            self.assertIn("R_ANNOTATION_GAP_FAVORS_RESCUE", rescue.rules_fired)
            self.assertIn("annotation_gap_hypotheses=1", rescue.input_facts)
            self.assertIn("reference_proteome", rescue.sample_requirement)
            self.assertEqual("high", rescue.expected_information_class)

    def test_llm_annotation_gap_hypothesis_cannot_trigger_annotation_rescue(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            hypotheses = root / "hypotheses.tsv"
            hypotheses.write_text(
                "hypothesis_id\thypothesis_type\tproposal_origin\treview_status\tstatus\trequired_experiment_ids\n"
                "H1\tannotation_missing\tllm\tpending_review\tpending_review\t\n",
                encoding="utf-8",
            )
            result = select_next_experiment(
                str(self._ranked(root, 4)), str(self._manifest(root)), str(root / "out"),
                hypotheses_arg=str(hypotheses),
            )
            self.assertNotEqual("annotation_rescue", result["primary_recommendation"])
            options = pd.read_csv(result["options"], sep="\t", dtype=str).fillna("").set_index("experiment_id")
            self.assertNotIn("R_ANNOTATION_GAP_FAVORS_RESCUE", options.loc["annotation_rescue", "rules_fired"])

    def test_annotation_rescue_can_be_required_by_hypothesis_experiment_ids(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            hypotheses = root / "hypotheses.tsv"
            hypotheses.write_text(
                "hypothesis_id\thypothesis_type\tproposal_origin\treview_status\tstatus\trequired_experiment_ids\n"
                "H1\tunknown_locus_mechanism\tmanual\treviewed\tready_to_test\tannotation_rescue\n",
                encoding="utf-8",
            )
            result = select_next_experiment(
                str(self._ranked(root, 4)), str(self._manifest(root)), str(root / "out"),
                hypotheses_arg=str(hypotheses),
            )
            self.assertEqual("annotation_rescue", result["primary_recommendation"])
            options = pd.read_csv(result["options"], sep="\t", dtype=str).fillna("").set_index("experiment_id")
            self.assertIn("linked_hypotheses=1", options.loc["annotation_rescue", "input_facts"])

    def test_non_gene_follow_up_ranks_first_on_low_recombination_region(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            mapping = root / "mapping.json"
            mapping.write_text(json.dumps({"statuses": ["low_recombination_region"], "primary_status": "low_recombination_region"}), encoding="utf-8")
            result = select_next_experiment(
                str(self._ranked(root, 4)), str(self._manifest(root)), str(root / "out"),
                mapping_gate_arg=str(mapping),
            )
            options = pd.read_csv(result["options"], sep="\t", dtype=str).fillna("").set_index("experiment_id")
            follow_up = options.loc["non_gene_follow_up"]
            marker = options.loc["additional_grid_markers"]
            self.assertEqual("non_gene_follow_up", result["primary_recommendation"])
            self.assertIn("R_NON_GENE_LOCUS_FAVORS_FOLLOW_UP", follow_up.rules_fired)
            self.assertEqual("high", follow_up.expected_information_class)
            self.assertNotIn("R_MAPPING_INSUFFICIENT_OR_RECOMBINANTS_FAVOR_INTERVAL_REFINEMENT", marker.rules_fired)
            self.assertEqual("5", marker.action_rank)

    def test_non_gene_follow_up_on_reviewed_linked_qtl_hypothesis(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            hypotheses = root / "hypotheses.tsv"
            hypotheses.write_text(
                "hypothesis_id\thypothesis_type\tproposal_origin\treview_status\tstatus\trequired_experiment_ids\n"
                "H1\tlinked_qtl\tmanual\treviewed\tready_to_test\t\n",
                encoding="utf-8",
            )
            result = select_next_experiment(
                str(self._ranked(root, 4)), str(self._manifest(root)), str(root / "out"),
                hypotheses_arg=str(hypotheses),
            )
            self.assertEqual("non_gene_follow_up", result["primary_recommendation"])
            options = pd.read_csv(result["options"], sep="\t", dtype=str).fillna("").set_index("experiment_id")
            self.assertIn("non_gene_locus_hypotheses=1", options.loc["non_gene_follow_up", "input_facts"])

    def test_linked_qtl_suspected_favors_interval_refinement_before_non_gene_follow_up(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            mapping = root / "mapping.json"
            mapping.write_text(json.dumps({"statuses": ["linked_qtl_suspected"], "primary_status": "linked_qtl_suspected"}), encoding="utf-8")
            result = select_next_experiment(
                str(self._ranked(root, 4)), str(self._manifest(root)), str(root / "out"),
                mapping_gate_arg=str(mapping),
            )
            self.assertEqual("additional_grid_markers", result["primary_recommendation"])
            options = pd.read_csv(result["options"], sep="\t", dtype=str).fillna("").set_index("experiment_id")
            marker = options.loc["additional_grid_markers"]
            self.assertIn("R_MAPPING_INSUFFICIENT_OR_RECOMBINANTS_FAVOR_INTERVAL_REFINEMENT", marker.rules_fired)
            self.assertLess(int(marker.action_rank), int(options.loc["non_gene_follow_up", "action_rank"]))

    def test_fine_mapping_required_token_triggers_grid_markers(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            mapping = root / "mapping.json"
            mapping.write_text(json.dumps({"statuses": ["fine_mapping_required"], "primary_status": "fine_mapping_required"}), encoding="utf-8")
            result = select_next_experiment(
                str(self._ranked(root, 4)), str(self._manifest(root)), str(root / "out"),
                mapping_gate_arg=str(mapping),
            )
            self.assertEqual("additional_grid_markers", result["primary_recommendation"])
            options = pd.read_csv(result["options"], sep="\t", dtype=str).fillna("").set_index("experiment_id")
            self.assertIn("R_MAPPING_INSUFFICIENT_OR_RECOMBINANTS_FAVOR_INTERVAL_REFINEMENT", options.loc["additional_grid_markers", "rules_fired"])

    def test_new_action_classes_keep_ordinal_decision_guarantees(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            hypotheses = root / "hypotheses.tsv"
            hypotheses.write_text(
                "hypothesis_id\thypothesis_type\tproposal_origin\treview_status\tstatus\trequired_experiment_ids\n"
                "H1\tannotation_missing\tmanual\treviewed\tready_to_test\t\n"
                "H2\tlinked_qtl\tmanual\treviewed\tready_to_test\t\n",
                encoding="utf-8",
            )
            mapping = root / "mapping.json"
            mapping.write_text(json.dumps({"statuses": ["low_recombination_region"]}), encoding="utf-8")
            result = select_next_experiment(
                str(self._ranked(root, 4)), str(self._manifest(root)), str(root / "out"),
                mapping_gate_arg=str(mapping), hypotheses_arg=str(hypotheses),
            )
            decision = json.loads(Path(result["decision"]).read_text(encoding="utf-8"))
            self.assertFalse(decision["causal_probability"])
            self.assertFalse(decision["bayesian_information_gain"])
            self.assertFalse(decision["unknown_values_imputed"])
            self.assertEqual("pending", decision["human_approval"])
            self.assertEqual("annotation_rescue", decision["primary_recommendation"]["experiment_id"])
            options = pd.read_csv(result["options"], sep="\t", dtype=str).fillna("").set_index("experiment_id")
            alternative_ids = {row["experiment_id"] for row in decision["alternatives"]}
            for experiment in ("annotation_rescue", "non_gene_follow_up"):
                self.assertIn(experiment, options.index)
                self.assertNotEqual("", options.loc[experiment, "hypotheses_addressed"])
            self.assertIn("non_gene_follow_up", alternative_ids)


class RemainingCandidatesTests(NextExperimentTests):
    """Issue #285: next-experiment consumes the canonical post-elimination
    remaining set, keeping pre/post counts auditable and never re-ranking."""

    def _remaining(self, root: Path, gene_ids: list[str]) -> Path:
        path = root / "remaining_candidate_set.tsv"
        path.write_text(
            "candidate_gene_id\tcurrent_gene_state\tsurviving_mechanisms\tfirst_wet_priority\trank_contribution\tnotes\n"
            + "".join(f"{gene}\tsurvives\tcoding;regulatory\tmedium\tnone\tsurvival is a triage state\n" for gene in gene_ids),
            encoding="utf-8")
        return path

    def _run(self, root: Path, ranked: Path, case: Path, remaining: Path | None) -> dict:
        return select_next_experiment(str(ranked), str(case), str(root / "out"),
                                      remaining_candidates_arg=str(remaining) if remaining else None)

    def test_post_elimination_counts_replace_the_pre_elimination_scale(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            # 20 ranked candidates, elimination left 4 (all dry_group=top).
            ranked = self._ranked(root, 20)
            case = self._manifest(root)
            remaining = self._remaining(root, [f"Gene{i:02d}" for i in range(4)])
            original = ranked.read_bytes()
            result = self._run(root, ranked, case, remaining)
            options = pd.read_csv(result["options"], sep="\t", dtype=str).fillna("")
            for experiment in ("top_n_parental_sanger", "multiplex_amplicon_sequencing", "target_capture"):
                facts = options[options.experiment_id.eq(experiment)].iloc[0]["input_facts"]
                self.assertIn("candidate_count=4", facts)
                self.assertIn("candidate_count_before_elimination=20", facts)
                self.assertIn("candidate_count_after_elimination=4", facts)
                self.assertIn("top_or_reserve_count=4", facts)
            decision = json.loads(Path(result["decision"]).read_text(encoding="utf-8"))
            self.assertIn("remaining_candidates", decision["input_fingerprints_sha256"])
            # The biological ranked input is untouched.
            self.assertEqual(original, ranked.read_bytes())

    def test_remaining_intersection_drives_top_or_reserve_count(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            # Remaining set includes reserve genes only: 12 ranked (6 top + 6
            # reserve), survivors = Gene06..Gene11 (reserve) -> top_or_reserve=6.
            ranked = self._ranked(root, 12)
            case = self._manifest(root)
            remaining = self._remaining(root, [f"Gene{i:02d}" for i in range(6, 12)])
            result = self._run(root, ranked, case, remaining)
            options = pd.read_csv(result["options"], sep="\t", dtype=str).fillna("")
            facts = options[options.experiment_id.eq("top_n_parental_sanger")].iloc[0]["input_facts"]
            self.assertIn("candidate_count=6", facts)
            self.assertIn("candidate_count_before_elimination=12", facts)
            self.assertIn("candidate_count_after_elimination=6", facts)
            self.assertIn("top_or_reserve_count=6", facts)

    def test_without_remaining_input_behavior_is_unchanged(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            ranked = self._ranked(root, 12)
            case = self._manifest(root)
            result = self._run(root, ranked, case, None)
            options = pd.read_csv(result["options"], sep="\t", dtype=str).fillna("")
            facts = options[options.experiment_id.eq("top_n_parental_sanger")].iloc[0]["input_facts"]
            self.assertIn("candidate_count=12", facts)
            self.assertNotIn("candidate_count_before_elimination", facts)
            decision = json.loads(Path(result["decision"]).read_text(encoding="utf-8"))
            self.assertNotIn("remaining_candidates", decision["input_fingerprints_sha256"])

    def test_empty_remaining_set_is_review_required_and_never_resurrects(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            ranked = self._ranked(root, 12)
            case = self._manifest(root)
            remaining = self._remaining(root, [])
            result = self._run(root, ranked, case, remaining)
            reviews = pd.read_csv(result["review_required"], sep="\t", dtype=str).fillna("")
            empty = reviews[reviews.review_reason.eq("remaining_candidates_empty")]
            self.assertEqual(1, len(empty))
            self.assertIn("do not resurrect", empty.iloc[0]["decision_effect"])
            options = pd.read_csv(result["options"], sep="\t", dtype=str).fillna("")
            self.assertIn("candidate_count=0", options.iloc[0]["input_facts"])
            self.assertIn("candidate_count_after_elimination=0", options.iloc[0]["input_facts"])

    def test_unknown_remaining_gene_is_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            ranked = self._ranked(root, 12)
            case = self._manifest(root)
            remaining = self._remaining(root, ["Gene00", "Gene99"])
            with self.assertRaises(ValueError) as context:
                self._run(root, ranked, case, remaining)
            self.assertIn("absent from the ranked table", str(context.exception))
            self.assertIn("Gene99", str(context.exception))


if __name__ == "__main__":
    unittest.main()
