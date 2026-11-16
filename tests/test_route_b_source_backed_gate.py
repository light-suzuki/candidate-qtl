"""Route B may not call a reverse lookup source-backed with nothing behind it (issue #476).

`route_b_readiness` / `process_inference_completeness` gate whether a Route B run
counts as an evaluated lane. The mask behind them used to ask only whether
``reference_source_id`` was a non-empty string, so an atlas whose reference rows
carry a product name and a pipeline label -- no FUNCTION, no PATHWAY, no
DISRUPTION PHENOTYPE -- reported ``source_backed_reverse_lookup``, the strongest
state, and was indistinguishable from a run backed by curated UniProt text.

The tests pin the distinction in both directions: content-free rows must not
claim it, and a focal-genome annotation that *does* carry function text keeps it
even though its ``local:annotation`` marker is not independently citable.
"""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from candidate_qtl.cli import main

ATLAS = "\n".join([
    "gene_id\tchrom\tstart\tend\tfamily_id\tfamily_label\torthogroup\torthology_relationship"
    "\torthology_reference_gene\tdomain_architecture\tbroad_function_bins\tevidence_completeness",
    "G_CLV1\tchr1\t1000\t3000\tCLV\tCLV1\tOG1\tone-to-one\tCLV1\tLRR-Pkinase\t\thigh_confidence_family_assignment",
    "G_PXY\tchr1\t20000\t22000\tPXY\tPXY1\tOG2\tone-to-one\tPXY1\tLRR-Pkinase\t\thigh_confidence_family_assignment",
]) + "\n"

# What an auto-triage panel produces with no annotation table: a product name from
# the FASTA defline, and the panel label as the "source".
PRODUCT_NAME_ONLY = "\n".join([
    "reference_gene_id\tfunction_description\tpathway\tmutant_phenotype\tsource_id",
    "CLV1\t\t\t\tlocal_fasta:arath",
    "PXY1\t\t\t\tlocal_fasta:arath",
]) + "\n"

CURATED = "\n".join([
    "reference_gene_id\tfunction_description\tpathway\tmutant_phenotype\tsource_id",
    "CLV1\treceptor kinase restricting the stem cell pool\tCLAVATA-WUSCHEL\tenlarged meristem\tUniProt:Q9XID0",
    "PXY1\treceptor kinase for procambial cell division\tCLE-receptor vascular signalling\t\tUniProt:Q9ZQK1",
]) + "\n"

# A focal-genome annotation: real function text, but the marker is not a citable
# accession. This is a supported Route B reference and must stay one.
FOCAL_ANNOTATION = "\n".join([
    "reference_gene_id\tfunction_description\tpathway\tmutant_phenotype\tsource_id",
    "CLV1\treceptor kinase restricting the stem cell pool\tCLAVATA-WUSCHEL\t\tlocal:annotation",
]) + "\n"


class RouteBSourceBackedGateTests(unittest.TestCase):
    def _run(self, functions: str) -> dict:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atlas = root / "atlas.tsv"
            atlas.write_text(ATLAS, encoding="utf-8")
            reference_function = root / "functions.tsv"
            reference_function.write_text(functions, encoding="utf-8")
            out = root / "route_b0"
            self.assertEqual(0, main(["interval-reverse-lookup", "--atlas", str(atlas),
                                      "--reference-function", str(reference_function),
                                      "--out-dir", str(out), "--case-id", "gate"]))
            return json.loads((out / "interval_first_status.json").read_text(encoding="utf-8"))

    def test_declared_source_without_function_text_is_not_source_backed(self) -> None:
        status = self._run(PRODUCT_NAME_ONLY)
        self.assertEqual("structural_providers_without_source_backed_function", status["reverse_lookup_status"])
        self.assertEqual("partial", status["process_inference_completeness"])
        coverage = status["provider_coverage"]
        self.assertEqual(0, coverage["source_backed_reference_function_gene_count"])
        # The rows exist and are reported; they are simply not counted as function.
        self.assertEqual(2, coverage["declared_reference_source_gene_count"])
        self.assertEqual(2, coverage["declared_source_without_function_text_gene_count"])
        self.assertEqual(0, coverage["reference_function_text_gene_count"])
        self.assertEqual(0, coverage["citable_reference_source_gene_count"])

    def test_curated_function_text_is_source_backed(self) -> None:
        status = self._run(CURATED)
        self.assertEqual("source_backed_reverse_lookup", status["reverse_lookup_status"])
        self.assertEqual("complete", status["process_inference_completeness"])
        coverage = status["provider_coverage"]
        self.assertEqual(2, coverage["source_backed_reference_function_gene_count"])
        self.assertEqual(2, coverage["reference_function_text_gene_count"])
        self.assertEqual(2, coverage["citable_reference_source_gene_count"])
        self.assertEqual(0, coverage["declared_source_without_function_text_gene_count"])

    def test_two_runs_differing_only_in_reference_function_do_not_share_a_verdict(self) -> None:
        # The regression this issue was filed for: same atlas, same providers, the
        # reference-function table swapped, and the old mask returned one verdict.
        self.assertNotEqual(
            self._run(PRODUCT_NAME_ONLY)["process_inference_completeness"],
            self._run(CURATED)["process_inference_completeness"],
        )

    def test_non_citable_focal_annotation_with_function_text_stays_source_backed(self) -> None:
        status = self._run(FOCAL_ANNOTATION)
        self.assertEqual("source_backed_reverse_lookup", status["reverse_lookup_status"])
        coverage = status["provider_coverage"]
        self.assertEqual(1, coverage["source_backed_reference_function_gene_count"])
        # Not citable, and reported as such -- without that costing it the lane.
        self.assertEqual(0, coverage["citable_reference_source_gene_count"])


if __name__ == "__main__":
    unittest.main()
