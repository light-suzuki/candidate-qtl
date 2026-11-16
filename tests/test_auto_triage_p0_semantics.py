from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from candidate_qtl.auto_triage import (
    _automatic_context_evidence,
    _functional_annotation_coverage,
    _rank_semantics,
    _reference_annotations,
    _supporting_source_summary,
    _trait_match_evidence,
)


class AutoTriageP0SemanticsTests(unittest.TestCase):
    def test_reference_corpus_schema_is_adapted_without_losing_rich_fields(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            annotation = root / "arabidopsis_reference_corpus.tsv"
            pd.DataFrame([
                {
                    "reference_gene_id": "Q9LZS3",
                    "description": "Example receptor kinase",
                    "function_description": "Controls shoot apical meristem stem-cell homeostasis.",
                    "go_terms": "GO:0010015;GO:0048507",
                    "mutant_phenotype": "Loss of function causes an enlarged meristem.",
                    "pathway": "CLAVATA signaling pathway.",
                    "tissue_specificity": "Expressed in shoot apical meristem.",
                    "developmental_stage": "Vegetative and floral meristems.",
                    "protein_family": "LRR receptor-like kinase family",
                    "source_id": "UniProt:Q9LZS3",
                }
            ]).to_csv(annotation, sep="\t", index=False)
            panel = root / "panel.yaml"
            panel.write_text("references: []\n", encoding="utf-8")

            records = _reference_annotations(
                {"label": "ath", "annotation": str(annotation)},
                panel,
                {"sp|Q9LZS3|EXAMPLE_ARATH": "sp|Q9LZS3|EXAMPLE_ARATH Example receptor kinase"},
            )
            row = records["sp|Q9LZS3|EXAMPLE_ARATH"]
            self.assertEqual(row["gene_description"], "Example receptor kinase")
            self.assertIn("stem-cell homeostasis", row["function_description"])
            self.assertIn("enlarged meristem", row["mutant_phenotype"])
            self.assertIn("CLAVATA", row["pathway"])
            self.assertIn("tissue_specificity:", row["expression"])
            self.assertIn("developmental_stage:", row["expression"])
            self.assertEqual(row["source_id"], "UniProt:Q9LZS3")
            self.assertEqual(row["retrieval_status"], "source_retrieved")

    def test_rich_annotation_without_resolvable_source_is_not_called_source_backed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            annotation = root / "unverified_annotation.tsv"
            pd.DataFrame([
                {
                    "reference_gene_id": "Q9LZS3",
                    "description": "Example receptor kinase",
                    "function_description": "meristem regulation",
                    "mutant_phenotype": "enlarged meristem",
                    "pathway": "CLV signaling",
                    "source_id": "",
                }
            ]).to_csv(annotation, sep="\t", index=False)
            panel = root / "panel.yaml"
            panel.write_text("references: []\n", encoding="utf-8")
            records = _reference_annotations(
                {"label": "ath", "annotation": str(annotation)},
                panel,
                {"sp|Q9LZS3|EXAMPLE_ARATH": "sp|Q9LZS3|EXAMPLE_ARATH Example receptor kinase"},
            )
            row = records["sp|Q9LZS3|EXAMPLE_ARATH"]
            self.assertEqual(row["retrieval_status"], "annotation_loaded_source_unverified")
            functions = pd.DataFrame([{**row, "candidate_gene_id": "g1", "reference_label": "ath", "review_status": "unreviewed"}])
            coverage = _functional_annotation_coverage(functions)
            self.assertEqual(coverage["source_backed_rich_rows"], 0)
            self.assertEqual(coverage["status"], "insufficient_source_backed_function")

    def test_annotation_coverage_distinguishes_fasta_only_from_source_backed_biology(self) -> None:
        functions = pd.DataFrame([
            {
                "gene_description": "FASTA-only name",
                "function_description": "",
                "go_terms": "",
                "mutant_phenotype": "",
                "pathway": "",
                "expression": "",
                "protein_family": "",
                "source_id": "local_fasta:ath",
                "retrieval_status": "fasta_description_only",
            },
            {
                "gene_description": "curated name",
                "function_description": "meristem regulation",
                "go_terms": "GO:0010015",
                "mutant_phenotype": "enlarged meristem",
                "pathway": "CLV signaling",
                "expression": "shoot apex",
                "protein_family": "LRR-RLK",
                "source_id": "UniProt:Q9LZS3",
                "retrieval_status": "source_retrieved",
            },
        ])
        coverage = _functional_annotation_coverage(functions)
        self.assertEqual(coverage["source_backed_rich_rows"], 1)
        self.assertEqual(coverage["status"], "source_backed_function_available")
        self.assertEqual(coverage["field_coverage"]["mutant_phenotype"]["non_empty"], 1)

    def test_generic_mapping_is_context_only_and_paralog_ambiguity_is_not_against(self) -> None:
        candidates = pd.DataFrame([{"gene_id": "g1"}, {"gene_id": "g2"}])
        orthology = pd.DataFrame([
            {
                "candidate_gene_id": "g1",
                "reference_label": "ath",
                "reference_gene_id": "AT1G00010",
                "mapping_status": "mapped",
                "relationship": "one_to_one_candidate",
                "orthology_confidence": "high",
                "rbh_status": "yes",
                "identity": "71.2",
                "query_coverage": "0.98",
                "subject_coverage": "0.97",
                "paralog_ambiguity": "none_detected",
                "top2_bitscore_ratio": "0.50",
            },
            {
                "candidate_gene_id": "g2",
                "reference_label": "ath",
                "reference_gene_id": "AT1G00020",
                "mapping_status": "mapped",
                "relationship": "one_to_many",
                "orthology_confidence": "low",
                "rbh_status": "no",
                "identity": "61.0",
                "query_coverage": "0.91",
                "subject_coverage": "0.90",
                "paralog_ambiguity": "top2_near_tie",
                "top2_bitscore_ratio": "0.96",
            },
        ])
        rows = _automatic_context_evidence(candidates, orthology, Path("candidates.tsv"))
        frame = pd.DataFrame(rows)
        self.assertTrue(frame["direction"].eq("context").all())
        self.assertEqual(int(frame["evidence_code"].eq("interval_membership").sum()), 2)
        self.assertEqual(int(frame["evidence_code"].eq("homology_context").sum()), 2)
        ambiguity = frame[frame["evidence_code"].eq("paralog_ambiguous")]
        self.assertEqual(len(ambiguity), 1)
        self.assertEqual(ambiguity.iloc[0]["direction"], "context")

    def test_source_backed_lexical_match_can_emit_low_confidence_support(self) -> None:
        row = SimpleNamespace(
            candidate_gene_id="g1",
            source_id="UniProt:P30924",
            retrieval_status="source_retrieved",
        )
        evidence = _trait_match_evidence(row, ["meristem"])
        self.assertEqual(evidence["evidence_code"], "automated_candidate_function_match")
        self.assertEqual(evidence["direction"], "support")
        self.assertEqual(evidence["confidence"], "low")

    def test_fasta_only_or_unverified_lexical_match_is_context_not_support(self) -> None:
        for source_id, retrieval_status in (
            ("local_fasta:ath", "fasta_description_only"),
            ("", "annotation_loaded_source_unverified"),
        ):
            with self.subTest(source_id=source_id, retrieval_status=retrieval_status):
                row = SimpleNamespace(
                    candidate_gene_id="g1",
                    source_id=source_id,
                    retrieval_status=retrieval_status,
                )
                evidence = _trait_match_evidence(row, ["meristem"])
                self.assertEqual(evidence["evidence_code"], "trait_context")
                self.assertEqual(evidence["direction"], "context")

    def test_supporting_source_summary_uses_support_rows_only(self) -> None:
        evidence = pd.DataFrame([
            {"gene_id": "g1", "direction": "context", "source_id": "local_fasta:ath"},
            {"gene_id": "g1", "direction": "support", "source_id": "UniProt:P30924"},
            {"gene_id": "g1", "direction": "support", "source_id": "UniProt:P30924"},
            {"gene_id": "g2", "direction": "context", "source_id": "UniProt:Q9SXK2"},
        ])
        summary = _supporting_source_summary(evidence)
        self.assertEqual(summary, {"g1": "UniProt:P30924"})

    def test_support_zero_matches_issue_397_penalty_only_contract(self) -> None:
        context_only = pd.DataFrame([
            {"direction": "context", "evidence_code": "interval_membership"},
            {"direction": "context", "evidence_code": "homology_context"},
        ])
        semantics = _rank_semantics(context_only)
        self.assertFalse(semantics["rank_metric_eligible"])
        self.assertTrue(semantics["stage_a_is_budget_only"])
        self.assertTrue(semantics["rank_is_penalty_only"])
        self.assertEqual(semantics["rank_semantic_status"], "not_evaluated_no_supporting_evidence")

        with_support = pd.DataFrame([
            {"direction": "context", "evidence_code": "homology_context"},
            {"direction": "support", "evidence_code": "automated_candidate_function_match"},
        ])
        semantics = _rank_semantics(with_support)
        self.assertTrue(semantics["rank_metric_eligible"])
        self.assertFalse(semantics["stage_a_is_budget_only"])
        self.assertEqual(semantics["rank_semantic_status"], "evaluable")
        self.assertNotIn("rank_is_penalty_only", semantics)

    def test_zero_evidence_uses_existing_no_scoring_state(self) -> None:
        semantics = _rank_semantics(pd.DataFrame(columns=["direction", "evidence_code"]))
        self.assertEqual(semantics["rank_semantic_status"], "not_evaluated_no_scoring_evidence")
        self.assertFalse(semantics["rank_metric_eligible"])
        self.assertTrue(semantics["stage_a_is_budget_only"])
        self.assertNotIn("rank_is_penalty_only", semantics)


if __name__ == "__main__":
    unittest.main()
