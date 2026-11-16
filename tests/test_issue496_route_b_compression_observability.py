from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from candidate_qtl.route_b_semantic_bundle import build_route_b_semantic_bundle


ATLAS_COLUMNS = [
    "gene_id",
    "evolutionary_family_id",
    "evolutionary_family_status",
    "orthogroup",
    "family_id",
    "domain_architecture",
    "orthology_reference_gene",
    "conservation_class",
    "evidence_completeness",
    "notes",
    "broad_function_bins",
]


def _write_atlas(path: Path, rows: list[dict[str, str]]) -> None:
    path.write_text(
        "\t".join(ATLAS_COLUMNS)
        + "\n"
        + "".join(
            "\t".join(str(row.get(column, "")) for column in ATLAS_COLUMNS) + "\n"
            for row in rows
        ),
        encoding="utf-8",
    )


def _run(rows: list[dict[str, str]]) -> dict[str, object]:
    tempdir = tempfile.TemporaryDirectory()
    root = Path(tempdir.name)
    atlas = root / "atlas.tsv"
    _write_atlas(atlas, rows)
    result = build_route_b_semantic_bundle(
        str(atlas), str(root / "out"), phenotype="test phenotype", case_id="issue496"
    )
    manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
    tempdir.cleanup()
    return manifest


class RouteBSemanticCompressionObservabilityTests(unittest.TestCase):
    def test_two_multi_member_units_report_twofold_compression(self) -> None:
        manifest = _run([
            {"gene_id": "g1", "evolutionary_family_id": "F1", "evolutionary_family_status": "assigned",
             "broad_function_bins": "receptor kinase regulates meristem development"},
            {"gene_id": "g2", "evolutionary_family_id": "F1", "evolutionary_family_status": "assigned"},
            {"gene_id": "g3", "orthogroup": "OG1",
             "broad_function_bins": "transcription factor controls floral transition"},
            {"gene_id": "g4", "orthogroup": "OG1"},
        ])
        metrics = manifest["compression_observability"]
        self.assertEqual("measured", metrics["status"])
        self.assertEqual(2.0, metrics["compression_gene_per_unit"])
        self.assertEqual(0.5, metrics["unit_to_gene_ratio"])
        self.assertEqual(0, metrics["single_member_unit_count"])
        self.assertEqual(2, metrics["multi_member_unit_count"])
        self.assertEqual(0.0, metrics["single_member_unit_fraction"])
        self.assertEqual(1.0, metrics["multi_member_unit_fraction"])
        self.assertEqual(2, metrics["expected_llm_request_count"])
        self.assertEqual(50.0, metrics["expected_requests_per_100_genes"])
        self.assertEqual("none", metrics["rank_contribution"])

    def test_four_singletons_report_no_compression(self) -> None:
        manifest = _run([
            {"gene_id": "g1"},
            {"gene_id": "g2"},
            {"gene_id": "g3"},
            {"gene_id": "g4"},
        ])
        metrics = manifest["compression_observability"]
        self.assertEqual(1.0, metrics["compression_gene_per_unit"])
        self.assertEqual(1.0, metrics["unit_to_gene_ratio"])
        self.assertEqual(4, metrics["single_member_unit_count"])
        self.assertEqual(0, metrics["multi_member_unit_count"])
        self.assertEqual(1.0, metrics["single_member_unit_fraction"])
        # Phase 2: these four synthetic genes carry no family, domain,
        # reference gene or literature, so none of them can answer a semantic
        # request. The units stay; the requests are not issued.
        self.assertEqual(4, metrics["semantic_unit_count"])
        self.assertEqual(0, metrics["requested_unit_count"])
        self.assertEqual(4, metrics["skipped_unit_count"])
        self.assertEqual(0, metrics["expected_llm_request_count"])
        self.assertEqual(0.0, metrics["expected_requests_per_100_genes"])

    def test_coverage_counts_are_reused_without_reclassification(self) -> None:
        manifest = _run([
            {"gene_id": "assigned", "family_id": "KINASE"},
            {"gene_id": "unassigned"},
            {"gene_id": "unresolved", "evolutionary_family_status": "unresolved"},
            {"gene_id": "failed", "evolutionary_family_status": "failed"},
        ])
        self.assertEqual(
            {"assigned": 1, "unassigned": 1, "unresolved": 1, "failed": 1},
            manifest["b0_coverage"]["counts"],
        )
        self.assertEqual(
            manifest["b0_coverage"]["counts"],
            manifest["compression_observability"]["coverage_gene_counts"],
        )

    def test_empty_interval_has_defined_zero_metrics(self) -> None:
        manifest = _run([])
        metrics = manifest["compression_observability"]
        self.assertEqual("not_applicable_empty_interval", metrics["status"])
        self.assertEqual(0.0, metrics["compression_gene_per_unit"])
        self.assertEqual(0.0, metrics["unit_to_gene_ratio"])
        self.assertEqual(0, metrics["single_member_unit_count"])
        self.assertEqual(0, metrics["multi_member_unit_count"])
        self.assertEqual(0, metrics["expected_llm_request_count"])
        self.assertEqual(0.0, metrics["expected_requests_per_100_genes"])


if __name__ == "__main__":
    unittest.main()
