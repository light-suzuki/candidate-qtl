from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import build_parser
from candidate_qtl.hypothesis_ledger import HYPOTHESIS_COLUMNS, update_hypothesis_ledger


class HypothesisLedgerTests(unittest.TestCase):
    def _hypothesis(self, hypothesis_id: str, hypothesis_type: str, genes: str = "GeneA", **extra) -> dict[str, str]:
        row = {column: "" for column in HYPOTHESIS_COLUMNS}
        row.update({
            "hypothesis_id": hypothesis_id, "case_id": "case01", "hypothesis_type": hypothesis_type,
            "member_gene_ids": genes, "proposed_mechanism": f"test {hypothesis_type}", "status": "proposed",
            "created_by": "manual:researcher", "proposal_origin": "manual", "review_status": "reviewed",
            "source_id": f"local:proposal:{hypothesis_id}",
        })
        row.update(extra)
        return row

    def _write(self, path: Path, rows: list[dict], columns: list[str] | None = None) -> Path:
        pd.DataFrame(rows, columns=columns).fillna("").to_csv(path, sep="\t", index=False, lineterminator="\n")
        return path

    def test_multiple_mechanisms_multigene_and_intergenic_are_representable(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            rows = [
                self._hypothesis("H1", "single_gene_coding"),
                self._hypothesis("H2", "promoter_or_utr"),
                self._hypothesis("H3", "multi_gene_combination", "GeneA;GeneB"),
                self._hypothesis("H4", "intergenic_regulatory", "", interval_or_element_id="chr1:1000-2000"),
            ]
            result = update_hypothesis_ledger(
                str(root / "out"), hypotheses_arg=str(self._write(root / "hypotheses.tsv", rows, HYPOTHESIS_COLUMNS))
            )
            ledger = pd.read_csv(result["locus_hypotheses"], sep="\t", dtype=str).fillna("")
            self.assertEqual(4, len(ledger))
            self.assertEqual(3, len(ledger[ledger.member_gene_ids.str.contains("GeneA")]))
            self.assertIn("intergenic_regulatory", set(ledger.hypothesis_type))

    def test_llm_proposal_is_forced_to_pending_review(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            row = self._hypothesis(
                "H1", "sv_cnv_pav", status="supported", created_by="llm:model",
                proposal_origin="llm", review_status="reviewed", source_id="model_request:req1",
            )
            result = update_hypothesis_ledger(
                str(root / "out"), hypotheses_arg=str(self._write(root / "hypotheses.tsv", [row], HYPOTHESIS_COLUMNS))
            )
            hypothesis = pd.read_csv(result["locus_hypotheses"], sep="\t", dtype=str).fillna("").iloc[0]
            reviews = pd.read_csv(result["hypothesis_review_required"], sep="\t", dtype=str).fillna("")
            history = pd.read_csv(result["hypothesis_history"], sep="\t", dtype=str).fillna("")
            self.assertEqual("pending_review", hypothesis.status)
            self.assertEqual("supported", hypothesis.requested_status)
            self.assertEqual("pending_review", hypothesis.review_status)
            self.assertIn("llm_proposal_cannot_self_promote", set(reviews.review_reason))
            self.assertEqual("rejected_guardrail_transition", history.iloc[0].transition_status)
            self.assertEqual("supported", history.iloc[0].requested_status)

    def test_missing_proposal_origin_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            row = self._hypothesis("H1", "single_gene_coding", proposal_origin="", created_by="gpt-model")
            with self.assertRaisesRegex(ValueError, "proposal_origin is required"):
                update_hypothesis_ledger(
                    str(root / "out"), hypotheses_arg=str(self._write(root / "hypotheses.tsv", [row], HYPOTHESIS_COLUMNS))
                )

    def test_evidence_and_wet_results_are_linked_by_id_without_payload_copy(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            row = self._hypothesis(
                "H1", "single_gene_coding", supporting_evidence_ids="obs_1;wet_1;sv_1",
                counterevidence_ids="obs_2", required_experiment_ids="top_n_parental_sanger",
            )
            evidence = self._write(root / "evidence.tsv", [{"observation_id": "obs_1", "payload": "large"}, {"observation_id": "obs_2", "payload": "counter"}])
            wet = self._write(root / "wet.tsv", [{"result_id": "wet_1", "result_value": "SNP"}])
            sv = self._write(root / "sv.tsv", [{"sv_event_id": "sv_1", "sv_type": "DEL"}])
            result = update_hypothesis_ledger(
                str(root / "out"), hypotheses_arg=str(self._write(root / "hypotheses.tsv", [row], HYPOTHESIS_COLUMNS)),
                evidence_arg=str(evidence), wet_results_arg=str(wet), sv_events_arg=str(sv),
            )
            links = pd.read_csv(result["hypothesis_evidence_links"], sep="\t", dtype=str).fillna("")
            self.assertEqual({"obs_1", "obs_2", "wet_1", "sv_1"}, set(links.evidence_id))
            self.assertEqual({"canonical_evidence", "wet_result", "sv_event"}, set(links.evidence_kind))
            self.assertNotIn("payload", links.columns)

    def test_evidence_namespace_collision_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            row = self._hypothesis("H1", "single_gene_coding", supporting_evidence_ids="same_id")
            evidence = self._write(root / "evidence.tsv", [{"observation_id": "same_id"}])
            wet = self._write(root / "wet.tsv", [{"result_id": "same_id"}])
            with self.assertRaisesRegex(ValueError, "namespace collision"):
                update_hypothesis_ledger(
                    str(root / "out"), hypotheses_arg=str(self._write(root / "hypotheses.tsv", [row], HYPOTHESIS_COLUMNS)),
                    evidence_arg=str(evidence), wet_results_arg=str(wet),
                )

    def test_cross_gene_evidence_link_requires_review(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            row = self._hypothesis("H1", "single_gene_coding", supporting_evidence_ids="obs_b")
            evidence = self._write(root / "evidence.tsv", [{"observation_id": "obs_b", "gene_id": "GeneB"}])
            result = update_hypothesis_ledger(
                str(root / "out"), hypotheses_arg=str(self._write(root / "hypotheses.tsv", [row], HYPOTHESIS_COLUMNS)),
                evidence_arg=str(evidence),
            )
            links = pd.read_csv(result["hypothesis_evidence_links"], sep="\t", dtype=str).fillna("")
            reviews = pd.read_csv(result["hypothesis_review_required"], sep="\t", dtype=str).fillna("")
            self.assertEqual("review_required", links.iloc[0].link_status)
            self.assertIn("evidence_gene_mismatch", set(reviews.review_reason))

    def test_coding_negative_weakens_mechanism_and_branches_without_changing_rank(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            ranked = root / "ranked.tsv"
            ranked.write_text("gene_id\tdry_rank\nGeneA\t1\n", encoding="utf-8")
            original = ranked.read_bytes()
            coding = self._hypothesis("H1", "single_gene_coding", status="testing")
            wet_state = self._write(root / "wet_state.tsv", [{
                "gene_id": "GeneA", "candidate_state": "coding_region_cleared",
                "regulatory_sv_state": "not_tested_retained", "source_assay_ids": "A1;A2",
            }])
            wet_results = self._write(root / "wet_results.tsv", [{
                "result_id": "wet_negative_1", "gene_id": "GeneA", "result_type": "parental_difference",
                "result_value": "no_difference",
            }])
            result = update_hypothesis_ledger(
                str(root / "out"), hypotheses_arg=str(self._write(root / "hypotheses.tsv", [coding], HYPOTHESIS_COLUMNS)),
                wet_candidate_state_arg=str(wet_state), wet_results_arg=str(wet_results), ranked_arg=str(ranked),
            )
            ledger = pd.read_csv(result["locus_hypotheses"], sep="\t", dtype=str).fillna("").set_index("hypothesis_type")
            self.assertEqual("weakened", ledger.loc["single_gene_coding", "status"])
            self.assertEqual("retained", ledger.loc["single_gene_coding", "gene_disposition"])
            self.assertEqual("coding_mechanism_only", ledger.loc["single_gene_coding", "rejection_scope"])
            self.assertIn("promoter_or_utr", ledger.index)
            self.assertIn("sv_cnv_pav", ledger.index)
            links = pd.read_csv(result["hypothesis_evidence_links"], sep="\t", dtype=str).fillna("")
            self.assertEqual({"against", "context"}, set(links.direction))
            self.assertTrue(links.evidence_kind.eq("wet_result").all())
            self.assertEqual(original, ranked.read_bytes())

    def test_partial_coverage_does_not_weaken_coding_hypothesis(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            coding = self._hypothesis("H1", "single_gene_coding", status="testing")
            wet_state = self._write(root / "wet_state.tsv", [{
                "gene_id": "GeneA", "candidate_state": "partially_covered",
                "regulatory_sv_state": "not_tested", "source_assay_ids": "A1",
            }])
            result = update_hypothesis_ledger(
                str(root / "out"), hypotheses_arg=str(self._write(root / "hypotheses.tsv", [coding], HYPOTHESIS_COLUMNS)),
                wet_candidate_state_arg=str(wet_state),
            )
            ledger = pd.read_csv(result["locus_hypotheses"], sep="\t", dtype=str).fillna("")
            self.assertEqual("testing", ledger.iloc[0].status)

    def test_rejected_status_without_recombinant_gate_is_downgraded(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            rejected = self._hypothesis("H1", "single_gene_coding", status="rejected")
            result = update_hypothesis_ledger(
                str(root / "out"), hypotheses_arg=str(self._write(root / "hypotheses.tsv", [rejected], HYPOTHESIS_COLUMNS))
            )
            ledger = pd.read_csv(result["locus_hypotheses"], sep="\t", dtype=str).fillna("").iloc[0]
            reviews = pd.read_csv(result["hypothesis_review_required"], sep="\t", dtype=str).fillna("")
            self.assertEqual("weakened", ledger.status)
            self.assertEqual("retained", ledger.gene_disposition)
            self.assertEqual("rejected", ledger.requested_status)
            self.assertIn("hypothesis_rejection_without_genetic_exclusion_gate", set(reviews.review_reason))

    def test_coding_hypothesis_can_be_rejected_after_full_gate_without_rejecting_gene(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            rejected = self._hypothesis("H1", "single_gene_coding", status="rejected")
            wet_state = self._write(root / "wet_state.tsv", [{
                "gene_id": "GeneA", "candidate_state": "coding_region_cleared",
                "regulatory_sv_state": "not_tested_retained",
            }])
            result = update_hypothesis_ledger(
                str(root / "out"), hypotheses_arg=str(self._write(root / "hypotheses.tsv", [rejected], HYPOTHESIS_COLUMNS)),
                wet_candidate_state_arg=str(wet_state),
            )
            ledger = pd.read_csv(result["locus_hypotheses"], sep="\t", dtype=str).fillna("").set_index("hypothesis_type")
            self.assertEqual("rejected", ledger.loc["single_gene_coding", "status"])
            self.assertEqual("retained", ledger.loc["single_gene_coding", "gene_disposition"])
            self.assertEqual("coding_mechanism_only", ledger.loc["single_gene_coding", "rejection_scope"])

    def test_multiple_unique_history_ids_are_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            rows = [self._hypothesis("H1", "single_gene_coding"), self._hypothesis("H2", "promoter_or_utr")]
            first = update_hypothesis_ledger(
                str(root / "first"), hypotheses_arg=str(self._write(root / "hypotheses.tsv", rows, HYPOTHESIS_COLUMNS))
            )
            second = update_hypothesis_ledger(
                str(root / "second"), hypotheses_arg=str(root / "hypotheses.tsv"),
                previous_history_arg=first["hypothesis_history"],
            )
            history = pd.read_csv(second["hypothesis_history"], sep="\t", dtype=str).fillna("")
            self.assertEqual(2, len(history))

    def test_hypothesis_relationship_cycle_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            rows = [
                self._hypothesis("H1", "single_gene_coding", parent_hypothesis_id="H2"),
                self._hypothesis("H2", "promoter_or_utr", parent_hypothesis_id="H1"),
            ]
            with self.assertRaisesRegex(ValueError, "cycle"):
                update_hypothesis_ledger(
                    str(root / "out"), hypotheses_arg=str(self._write(root / "hypotheses.tsv", rows, HYPOTHESIS_COLUMNS))
                )

    def test_forbidden_status_transition_is_retained_in_append_only_history(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            superseded = self._hypothesis("H1", "single_gene_coding", status="superseded")
            first = update_hypothesis_ledger(
                str(root / "first"),
                hypotheses_arg=str(self._write(root / "superseded.tsv", [superseded], HYPOTHESIS_COLUMNS)),
            )
            proposed = self._hypothesis("H1", "single_gene_coding", status="proposed", status_reason="reopen")
            second = update_hypothesis_ledger(
                str(root / "second"),
                hypotheses_arg=str(self._write(root / "proposed.tsv", [proposed], HYPOTHESIS_COLUMNS)),
                previous_history_arg=first["hypothesis_history"],
            )
            ledger = pd.read_csv(second["locus_hypotheses"], sep="\t", dtype=str).fillna("").iloc[0]
            history = pd.read_csv(second["hypothesis_history"], sep="\t", dtype=str).fillna("")
            self.assertEqual("superseded", ledger.status)
            self.assertIn("rejected_forbidden_transition", set(history.transition_status))

    def test_cli_exposes_hypothesis_ledger(self) -> None:
        args = build_parser().parse_args(["hypothesis-ledger", "--out-dir", "out", "--case-id", "case01"])
        self.assertEqual("hypothesis-ledger", args.command)


if __name__ == "__main__":
    unittest.main()
