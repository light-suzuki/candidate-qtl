from __future__ import annotations

import unittest

import numpy as np
import pandas as pd

from candidate_qtl.evidence_dependency import (
    annotate_evidence_dependencies,
    dependency_cluster_summary,
    dependency_reporting_counts,
)


class EvidenceDependencyTests(unittest.TestCase):
    def test_same_root_collapses_to_one_cluster(self) -> None:
        evidence = pd.DataFrame([
            {"source_id": "routeA", "route": "A", "evidence_kind": "primary_observation", "root_source_id": "ROOT1"},
            {"source_id": "routeB", "route": "B", "evidence_kind": "orthology_transfer", "root_source_id": "ROOT1"},
            {"source_id": "routeC", "route": "C", "evidence_kind": "curated_assertion", "root_source_id": "ROOT2"},
        ])
        annotated = annotate_evidence_dependencies(evidence)
        self.assertEqual(1, annotated[annotated.root_source_id.eq("ROOT1")].dependency_cluster.nunique())
        counts = dependency_reporting_counts(evidence)
        self.assertEqual(3, counts["evidence_count"])
        self.assertEqual(2, counts["independent_evidence_cluster_count"])
        summary = dependency_cluster_summary(evidence)
        root1 = summary[summary.root_source_id.eq("ROOT1")].iloc[0]
        self.assertEqual(2, root1.evidence_count)
        self.assertEqual("A;B", root1.routes)

    def test_unknown_root_is_not_assumed_independent(self) -> None:
        evidence = pd.DataFrame([
            {"source_id": "SRC1", "evidence_kind": "curated_assertion"},
            {"source_id": "SRC2", "evidence_kind": "prediction", "root_source_id": ""},
            {"source_id": "SRC3", "evidence_kind": "primary_observation", "root_source_id": "ROOT3"},
        ])
        annotated = annotate_evidence_dependencies(evidence)
        unknown = annotated[annotated.dependency_status.eq("dependency_unknown")]
        self.assertEqual(2, len(unknown))
        self.assertTrue(unknown.dependency_cluster.eq("dependency_unknown").all())
        counts = dependency_reporting_counts(evidence)
        self.assertEqual(1, counts["independent_evidence_cluster_count"])
        self.assertEqual(2, counts["dependency_unknown_evidence_count"])
        summary = dependency_cluster_summary(evidence)
        self.assertEqual(2, summary[summary.dependency_status.eq("dependency_unknown")].iloc[0].evidence_count)

    def test_legacy_table_is_backward_compatible(self) -> None:
        legacy = pd.DataFrame([
            {"gene_id": "G1", "evidence_code": "E1", "source_id": "SRC1"},
            {"gene_id": "G2", "evidence_code": "E2", "source_id": "SRC2"},
        ])
        annotated = annotate_evidence_dependencies(legacy)
        self.assertTrue(annotated.evidence_kind.eq("unknown").all())
        self.assertTrue(annotated.dependency_status.eq("dependency_unknown").all())
        self.assertEqual(0, dependency_reporting_counts(legacy)["independent_evidence_cluster_count"])

    def test_source_id_is_not_promoted_to_root(self) -> None:
        evidence = pd.DataFrame([{"source_id": "SRC42", "evidence_kind": "primary_observation"}])
        annotated = annotate_evidence_dependencies(evidence)
        self.assertEqual("", annotated.iloc[0].root_source_id)
        self.assertEqual("dependency_unknown", annotated.iloc[0].dependency_status)

    def test_human_readable_kind_is_normalized(self) -> None:
        evidence = pd.DataFrame([{"evidence_kind": "sequence domain call", "root_source_id": "ROOT1"}])
        self.assertEqual("sequence_domain_call", annotate_evidence_dependencies(evidence).iloc[0].evidence_kind)

    def test_invalid_kind_fails_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "unsupported evidence_kind"):
            annotate_evidence_dependencies(pd.DataFrame([{"evidence_kind": "invalid_kind", "root_source_id": "ROOT1"}]))

    def test_reporting_preserves_score_columns(self) -> None:
        evidence = pd.DataFrame([
            {"gene_id": "G1", "dry_score": "4.5", "dry_rank": "1", "root_source_id": "ROOT1"},
            {"gene_id": "G2", "dry_score": "2.0", "dry_rank": "2", "root_source_id": "ROOT2"},
        ])
        before = evidence[["gene_id", "dry_score", "dry_rank"]].copy()
        annotated = annotate_evidence_dependencies(evidence)
        pd.testing.assert_frame_equal(before.reset_index(drop=True), annotated[["gene_id", "dry_score", "dry_rank"]].reset_index(drop=True))
        self.assertEqual("none", dependency_reporting_counts(evidence)["rank_contribution"])

    def test_all_pandas_scalar_missing_values_are_empty(self) -> None:
        evidence = pd.DataFrame([
            {"source_id": "S1", "evidence_kind": pd.NA, "root_source_id": pd.NA},
            {"source_id": "S2", "evidence_kind": pd.NaT, "root_source_id": pd.NaT},
            {"source_id": "S3", "evidence_kind": np.nan, "root_source_id": np.nan},
            {"source_id": "S4", "evidence_kind": None, "root_source_id": None},
            {"source_id": "S5", "evidence_kind": "curated_assertion", "root_source_id": "ROOT5"},
        ])
        annotated = annotate_evidence_dependencies(evidence)
        # No "<NA>"/"nan" leakage: every missing cell renders as empty.
        rendered = annotated.astype(str)
        self.assertFalse(rendered.apply(lambda column: column.str.contains("<NA>")).any().any())
        self.assertTrue(annotated.root_source_id.eq("").head(4).all())
        self.assertTrue(annotated.evidence_kind.eq("unknown").head(4).all())
        self.assertEqual(4, annotated.dependency_status.eq("dependency_unknown").sum())
        counts = dependency_reporting_counts(evidence)
        self.assertEqual(5, counts["evidence_count"])
        self.assertEqual(1, counts["independent_evidence_cluster_count"])
        self.assertEqual(4, counts["dependency_unknown_evidence_count"])
        self.assertTrue(counts["dependency_unknown_present"])
        summary = dependency_cluster_summary(evidence)
        unknown = summary[summary.dependency_status.eq("dependency_unknown")].iloc[0]
        self.assertEqual(4, unknown.evidence_count)

    def test_additional_scalar_edge_cases(self) -> None:
        # np.float64(np.nan) is a missing value, like np.nan / pd.NA.
        evidence = pd.DataFrame([{"source_id": "S1", "evidence_kind": np.float64(np.nan),
                                  "root_source_id": np.float64(np.nan)}])
        annotated = annotate_evidence_dependencies(evidence)
        self.assertEqual("unknown", annotated.iloc[0].evidence_kind)
        self.assertEqual("", annotated.iloc[0].root_source_id)
        self.assertEqual("dependency_unknown", annotated.iloc[0].dependency_status)
        # pd.Timestamp is NOT missing: it renders as its (non-empty) string
        # when used as free-text root_source_id.
        timestamp = pd.Timestamp("2024-01-02 03:04:05")
        evidence = pd.DataFrame([{"source_id": "S2", "evidence_kind": "primary_observation",
                                  "root_source_id": timestamp}])
        annotated = annotate_evidence_dependencies(evidence)
        self.assertEqual(str(timestamp), annotated.iloc[0].root_source_id)
        self.assertNotEqual("", annotated.iloc[0].root_source_id)
        self.assertEqual("root_source_known", annotated.iloc[0].dependency_status)
        # bool scalars are not missing either; they render as text.
        evidence = pd.DataFrame([{"source_id": "S3", "evidence_kind": "curated_assertion",
                                  "root_source_id": True}])
        annotated = annotate_evidence_dependencies(evidence)
        self.assertEqual("True", annotated.iloc[0].root_source_id)
        self.assertEqual("root_source_known", annotated.iloc[0].dependency_status)
        # Compound cells stay fail-closed on the kind column too (pd.Series /
        # numpy.ndarray -> "must be a scalar" ValueError).
        for column, value in (("evidence_kind", pd.Series([1, 2])), ("evidence_kind", np.array([1, 2])),
                              ("root_source_id", pd.Series([1, 2])), ("root_source_id", np.array([1, 2]))):
            with self.assertRaisesRegex(ValueError, "must be a scalar"):
                annotate_evidence_dependencies(pd.DataFrame([{"source_id": "S4", column: value}]))

    def test_compound_cell_values_fail_closed_with_scalar_error(self) -> None:
        compound_values = ([1, 2], ("a", "b"), {"k": "v"}, {"a", "b"}, np.array([1, 2]), pd.Series([1, 2]))
        for column in ("evidence_kind", "root_source_id"):
            for value in compound_values:
                evidence = pd.DataFrame([{"source_id": "S1", column: value}])
                with self.assertRaisesRegex(ValueError, "must be a scalar"):
                    annotate_evidence_dependencies(evidence)
                with self.assertRaisesRegex(ValueError, "must be a scalar"):
                    dependency_reporting_counts(evidence)
                with self.assertRaisesRegex(ValueError, "must be a scalar"):
                    dependency_cluster_summary(evidence)


if __name__ == "__main__":
    unittest.main()
