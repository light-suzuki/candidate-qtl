from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import build_parser
from candidate_qtl.wet_state import ASSAY_COLUMNS, COVERAGE_COLUMNS, RESULT_COLUMNS, update_wet_state


class WetStateTests(unittest.TestCase):
    def _write(self, path: Path, rows: list[dict], columns: list[str]) -> Path:
        pd.DataFrame(rows, columns=columns).fillna("").to_csv(path, sep="\t", index=False, lineterminator="\n")
        return path

    def _assay(self, assay_id: str, gene: str, state: str, region: str = "chr1:100-199", region_class: str = "cds", **extra) -> dict:
        row = {column: "" for column in ASSAY_COLUMNS}
        row.update({
            "assay_id": assay_id, "gene_id": gene, "transcript_id": f"{gene}.1", "assay_type": "sanger",
            "planned_region": region, "region_class": region_class, "assay_state": state,
            "source_ref": f"plan:{assay_id}", "review_status": "reviewed", "reviewer": "tester",
        })
        row.update(extra)
        return row

    def _coverage(self, assay_id: str, gene: str, fraction: str = "1", **extra) -> dict:
        row = {column: "" for column in COVERAGE_COLUMNS}
        row.update({
            "assay_id": assay_id, "gene_id": gene, "transcript_id": f"{gene}.1",
            "planned_region": "chr1:100-199", "actual_covered_region": "chr1:100-199",
            "coverage_fraction": fraction, "region_class": "cds", "required_region": "yes",
            "quality_status": "pass", "read_direction": "bidirectional", "off_target_status": "none",
            "multi_product_status": "single", "source_ref": f"chromatogram:{assay_id}",
        })
        row.update(extra)
        return row

    def _result(self, result_id: str, assay_id: str, gene: str, role: str, result_type: str, value: str, **extra) -> dict:
        row = {column: "" for column in RESULT_COLUMNS}
        row.update({
            "result_id": result_id, "assay_id": assay_id, "gene_id": gene, "sample_id": f"sample:{role}",
            "sample_role": role, "result_type": result_type, "result_value": value,
            "quality_status": "pass", "source_ref": f"source:{result_id}", "review_status": "reviewed", "reviewer": "tester",
        })
        row.update(extra)
        return row

    def test_parent_non_amplification_is_sv_warning_not_negative(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            assays = self._write(root / "assays.tsv", [self._assay("A1", "GeneA", "failed", failure_reason="parent B no amplification")], ASSAY_COLUMNS)
            results = self._write(root / "results.tsv", [self._result("R1", "A1", "GeneA", "parent_b", "amplification", "no amplification")], RESULT_COLUMNS)
            out = root / "out"
            update_wet_state(str(out), assay_registry_arg=str(assays), results_arg=str(results))
            state = pd.read_csv(out / "candidate_wet_state.tsv", sep="\t", dtype=str).fillna("").iloc[0]
            review = pd.read_csv(out / "wet_review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("inconclusive", state.candidate_state)
            self.assertEqual("sv_pav_or_primer_mismatch_warning", state.regulatory_sv_state)
            self.assertEqual("no", state.reserve_transition_allowed)
            self.assertIn("parent_non_amplification_not_negative", set(review.review_reason))

    def test_full_reviewed_bidirectional_coding_negative_is_scope_limited(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            assays = [
                self._assay("A1", "GeneA", "completed", "chr1:100-199", "cds"),
                self._assay("A2", "GeneA", "completed", "chr1:300-349", "splice_boundary"),
            ]
            coverage = [
                self._coverage("A1", "GeneA"),
                self._coverage("A2", "GeneA", planned_region="chr1:300-349", actual_covered_region="chr1:300-349", region_class="splice_boundary"),
            ]
            results = [
                self._result("R1", "A1", "GeneA", "parent_a", "parental_difference", "no_difference"),
                self._result("R2", "A1", "GeneA", "parent_b", "parental_difference", "no_difference"),
                self._result("R3", "A2", "GeneA", "parent_a", "parental_difference", "no_difference"),
                self._result("R4", "A2", "GeneA", "parent_b", "parental_difference", "no_difference"),
            ]
            out = root / "out"
            update_wet_state(
                str(out),
                assay_registry_arg=str(self._write(root / "assays.tsv", assays, ASSAY_COLUMNS)),
                coverage_arg=str(self._write(root / "coverage.tsv", coverage, COVERAGE_COLUMNS)),
                results_arg=str(self._write(root / "results.tsv", results, RESULT_COLUMNS)),
            )
            state = pd.read_csv(out / "candidate_wet_state.tsv", sep="\t", dtype=str).fillna("").iloc[0]
            self.assertEqual("coding_region_cleared", state.candidate_state)
            self.assertEqual("not_tested_retained", state.regulatory_sv_state)
            self.assertEqual("scope_limited_coding_only", state.reserve_transition_allowed)
            self.assertEqual("no", state.hypothesis_reject_allowed)

    def test_partial_negative_cannot_request_reject(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            assays = [self._assay("A1", "GeneA", "partial_result"), self._assay("A2", "GeneA", "wet_pending", "chr1:300-399")]
            coverage = [self._coverage("A1", "GeneA")]
            results = [self._result("R1", "A1", "GeneA", "parent_a", "parental_difference", "no_difference", requested_candidate_action="reject")]
            out = root / "out"
            update_wet_state(
                str(out), assay_registry_arg=str(self._write(root / "assays.tsv", assays, ASSAY_COLUMNS)),
                coverage_arg=str(self._write(root / "coverage.tsv", coverage, COVERAGE_COLUMNS)),
                results_arg=str(self._write(root / "results.tsv", results, RESULT_COLUMNS)),
            )
            state = pd.read_csv(out / "candidate_wet_state.tsv", sep="\t", dtype=str).fillna("").iloc[0]
            review = pd.read_csv(out / "wet_review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("no_parental_difference", state.candidate_state)
            self.assertEqual("fail", state.coverage_quality_gate)
            self.assertEqual("no", state.hypothesis_reject_allowed)
            self.assertIn("candidate_transition_requested_without_coverage_quality_gate", set(review.review_reason))

    def test_every_required_assay_needs_reviewed_negative_before_clearance(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            assays = [self._assay("A1", "GeneA", "completed"), self._assay("A2", "GeneA", "completed", "chr1:300-399")]
            coverage = [self._coverage("A1", "GeneA"), self._coverage("A2", "GeneA", planned_region="chr1:300-399", actual_covered_region="chr1:300-399")]
            results = [self._result("R1", "A1", "GeneA", "parent_a", "parental_difference", "no_difference")]
            out = root / "out"
            update_wet_state(
                str(out), assay_registry_arg=str(self._write(root / "assays.tsv", assays, ASSAY_COLUMNS)),
                coverage_arg=str(self._write(root / "coverage.tsv", coverage, COVERAGE_COLUMNS)),
                results_arg=str(self._write(root / "results.tsv", results, RESULT_COLUMNS)),
            )
            state = pd.read_csv(out / "candidate_wet_state.tsv", sep="\t", dtype=str).fillna("").iloc[0]
            self.assertEqual("no_parental_difference", state.candidate_state)
            self.assertNotEqual("coding_region_cleared", state.coding_region_state)

    def test_forbidden_assay_transition_is_history_event_and_state_is_retained(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            first = root / "first"
            planned = self._write(root / "planned.tsv", [self._assay("A1", "GeneA", "planned")], ASSAY_COLUMNS)
            update_wet_state(str(first), assay_registry_arg=str(planned))
            second = root / "second"
            completed = self._write(root / "completed.tsv", [self._assay("A1", "GeneA", "completed")], ASSAY_COLUMNS)
            update_wet_state(str(second), assay_registry_arg=str(completed), previous_history_arg=str(first / "wet_state_history.tsv"))
            registry = pd.read_csv(second / "assay_registry.tsv", sep="\t", dtype=str).fillna("")
            history = pd.read_csv(second / "wet_state_history.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("planned", registry.iloc[0].assay_state)
            self.assertIn("rejected_forbidden_transition", set(history.transition_status))
            self.assertEqual(2, len(history))  # assay initial plus rejected assay transition; candidate remained not_tested

    def test_cosegregation_failure_and_recombinant_exclusion_are_distinct(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            assays = [self._assay("A1", "GeneA", "completed"), self._assay("B1", "GeneB", "completed")]
            results = [
                self._result("R1", "A1", "GeneA", "f2", "validation", "discordant", validation_result="fails_cosegregation", individual_id="F2-001", generation="F2"),
                self._result("R2", "B1", "GeneB", "recombinant", "validation", "breakpoint excludes locus", validation_result="excluded_by_recombinant", individual_id="REC-001"),
            ]
            out = root / "out"
            update_wet_state(
                str(out), assay_registry_arg=str(self._write(root / "assays.tsv", assays, ASSAY_COLUMNS)),
                results_arg=str(self._write(root / "results.tsv", results, RESULT_COLUMNS)),
            )
            states = pd.read_csv(out / "candidate_wet_state.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            self.assertEqual("fails_cosegregation", states.loc["GeneA", "candidate_state"])
            self.assertEqual("excluded_by_recombinant", states.loc["GeneB", "candidate_state"])

    def test_individual_validation_requires_individual_id(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            assays = self._write(root / "assays.tsv", [self._assay("A1", "GeneA", "completed")], ASSAY_COLUMNS)
            results = self._write(
                root / "results.tsv",
                [self._result("R1", "A1", "GeneA", "f2", "validation", "supports", validation_result="supports_qtl")],
                RESULT_COLUMNS,
            )
            with self.assertRaisesRegex(ValueError, "individual_id"):
                update_wet_state(str(root / "out"), assay_registry_arg=str(assays), results_arg=str(results))

    def test_explicit_cosegregation_field_normalizes_validation_state(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            assays = self._write(root / "assays.tsv", [self._assay("A1", "GeneA", "completed")], ASSAY_COLUMNS)
            result = self._result(
                "R1", "A1", "GeneA", "f2", "validation", "discordant",
                cosegregation_result="fails_cosegregation", individual_id="F2-001", generation="F2",
            )
            results = self._write(root / "results.tsv", [result], RESULT_COLUMNS)
            out = root / "out"
            update_wet_state(str(out), assay_registry_arg=str(assays), results_arg=str(results))
            state = pd.read_csv(out / "candidate_wet_state.tsv", sep="\t", dtype=str).fillna("").iloc[0]
            self.assertEqual("fails_cosegregation", state.candidate_state)

    def test_unreviewed_confirmed_variant_is_not_promoted(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            assays = self._write(root / "assays.tsv", [self._assay("A1", "GeneA", "completed")], ASSAY_COLUMNS)
            result = self._result(
                "R1", "A1", "GeneA", "parent_a", "variant", "SNP", variant_confirmation_status="confirmed",
                review_status="unreviewed",
            )
            results = self._write(root / "results.tsv", [result], RESULT_COLUMNS)
            out = root / "out"
            update_wet_state(str(out), assay_registry_arg=str(assays), results_arg=str(results))
            state = pd.read_csv(out / "candidate_wet_state.tsv", sep="\t", dtype=str).fillna("").iloc[0]
            review = pd.read_csv(out / "wet_review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertNotEqual("candidate_variant_confirmed", state.candidate_state)
            self.assertIn("variant_confirmation_without_reviewed_quality_source", set(review.review_reason))

    def test_exon_design_coverage_never_clears_wet_coding_state(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            design = root / "coverage_windows.tsv"
            design.write_text("gene_id\ttarget_interval\tstatus\toutput_json\nGeneA\tchr1:100-199\tcovered\tdesign.json\n", encoding="utf-8")
            out = root / "out"
            update_wet_state(str(out), exon_coverage_arg=str(design))
            state = pd.read_csv(out / "candidate_wet_state.tsv", sep="\t", dtype=str).fillna("").iloc[0]
            coverage = pd.read_csv(out / "assay_target_coverage.tsv", sep="\t", dtype=str).fillna("").iloc[0]
            self.assertNotEqual("coding_region_cleared", state.candidate_state)
            self.assertEqual("fail", state.coverage_quality_gate)
            self.assertEqual("design_only", coverage.quality_status)
            self.assertEqual("", coverage.actual_covered_region)

    def test_primer_import_cannot_claim_wet_completion(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            primer = root / "primer_candidates.tsv"
            primer.write_text(
                "marker_id\tgene_id\tpair_specificity\tcoverage_status\twet_validation_status\n"
                "M1\tGeneA\tunique\tcovered\tcompleted\n",
                encoding="utf-8",
            )
            out = root / "out"
            update_wet_state(str(out), primer_results_arg=str(primer))
            assay = pd.read_csv(out / "assay_registry.tsv", sep="\t", dtype=str).fillna("").iloc[0]
            self.assertEqual("primer_designed", assay.assay_state)

    def test_cli_exposes_wet_state(self) -> None:
        args = build_parser().parse_args(["wet-state", "--out-dir", "out", "--assays", "assays.tsv"])
        self.assertEqual("wet-state", args.command)
        self.assertEqual(1.0, args.min_coverage)


if __name__ == "__main__":
    unittest.main()
