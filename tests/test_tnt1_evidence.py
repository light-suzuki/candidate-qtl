"""Regression tests for #524 Tnt1 evidence semantics."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.tnt1_evidence import normalize_tnt1_evidence


HEADER = (
    "reported_gene_id\tgenetic_background\tmutant_line\tfst_id\tinsertion_site\t"
    "insertion_observed\tpublished_mutant_phenotype\tsegregation_supported\t"
    "independent_alleles_supported\tcomplementation_or_rescue_supported\t"
    "phenotype_or_function\tgene_symbol\tontology_ids\tsource_ids\tsource_url\tnotes"
)

ROWS = [
    # FST existence alone: observable but never causal function transfer.
    "medtr.R108.gnmHiC_1.ann1.G1\tR108\tNF1\tFST1\tchr1:10\tsupported\tnot_evaluated\tnot_evaluated\tnot_evaluated\tnot_evaluated\t\t\t\t\thttps://medicago-mutant.okstate.edu/\tFST only",
    # Published phenotype plus segregation/genotype support: causal-ready on an
    # exact R108-native ID with a publication identifier.
    "medtr.R108.gnmHiC_1.ann1.G2\tR108\tNF2\tFST2\tchr1:20\tsupported\tsupported\tsupported\tnot_evaluated\tnot_evaluated\tshort root phenotype\tMtG2\tTO:0000001\tPMID:12345\thttps://example.org/paper\tsegregating line",
    # Multiple independent alleles are a stronger causal class.
    "medtr.R108.gnmHiC_1.ann1.G3\tR108\tNF3\tFST3\tchr1:30\tsupported\tsupported\tsupported\tsupported\tnot_evaluated\tnodulation defect\tMtG3\t\tDOI:10.0000/example\thttps://example.org/alleles\ttwo alleles",
    # Historical A17-mapped FST from an R108 line: never silently treated as an
    # R108-native gene even when causal evidence is otherwise strong.
    "medtr.A17_HM341.gnm4.ann2.G4\tR108\tNF4\tFST4\tchr1:40\tsupported\tsupported\tsupported\tsupported\tsupported\tflowering defect\tMtG4\t\tPMID:67890\thttps://example.org/a17-map\trequires bridge",
    # Invalid state must be audited rather than guessed.
    "medtr.R108.gnmHiC_1.ann1.G5\tR108\tNF5\tFST5\tchr1:50\tmaybe\tnot_evaluated\tnot_evaluated\tnot_evaluated\tnot_evaluated\t\t\t\t\thttps://example.org/fst\tinvalid token",
]


class Tnt1EvidenceTests(unittest.TestCase):
    def _input(self, root: Path) -> Path:
        path = root / "tnt1.tsv"
        path.write_text(HEADER + "\n" + "\n".join(ROWS) + "\n", encoding="utf-8")
        return path

    def test_evidence_classes_and_identity_bridge_boundary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = self._input(root)
            out = root / "out"
            result = normalize_tnt1_evidence(source, out, case_id="pea_k1")

            ledger = pd.read_csv(result["ledger"], sep="\t", dtype=str).fillna("")
            functional = pd.read_csv(result["normalized_functional_evidence"], sep="\t", dtype=str).fillna("")
            audit = pd.read_csv(result["audit"], sep="\t", dtype=str).fillna("")
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))

            fst_only = ledger[ledger.reported_gene_id.eq("medtr.R108.gnmHiC_1.ann1.G1")].iloc[0]
            self.assertEqual("insertion_observed_only", fst_only.evidence_class)
            self.assertEqual("insertion_only", fst_only.causal_gene_support_state)
            self.assertEqual("false", fst_only.route_b_function_transfer_ready)
            self.assertEqual("causal_gene_support_not_established", fst_only.route_b_block_reason)

            genotype = ledger[ledger.reported_gene_id.eq("medtr.R108.gnmHiC_1.ann1.G2")].iloc[0]
            self.assertEqual("genotype_supported_mutant", genotype.evidence_class)
            self.assertEqual("causal_supported", genotype.causal_gene_support_state)
            self.assertEqual("native_background_match", genotype.identity_relation_state)
            self.assertEqual("true", genotype.route_b_function_transfer_ready)

            strong = ledger[ledger.reported_gene_id.eq("medtr.R108.gnmHiC_1.ann1.G3")].iloc[0]
            self.assertEqual("strong_causal", strong.evidence_class)
            self.assertEqual("true", strong.route_b_function_transfer_ready)

            cross_accession = ledger[ledger.reported_gene_id.eq("medtr.A17_HM341.gnm4.ann2.G4")].iloc[0]
            self.assertEqual("strong_causal", cross_accession.evidence_class)
            self.assertEqual("requires_identity_bridge", cross_accession.identity_relation_state)
            self.assertEqual("false", cross_accession.route_b_function_transfer_ready)
            self.assertEqual(
                "reported_gene_id_requires_explicit_background_identity_bridge",
                cross_accession.route_b_block_reason,
            )

            self.assertEqual(
                {"medtr.R108.gnmHiC_1.ann1.G2", "medtr.R108.gnmHiC_1.ann1.G3"},
                set(functional.native_gene_id),
            )
            self.assertTrue(functional.evidence_type.eq("tnt1_curated_mutant").all())
            self.assertTrue(functional.mutant_or_perturbation_background.eq("R108").all())

            invalid = audit[audit.reported_gene_id.eq("medtr.R108.gnmHiC_1.ann1.G5")].iloc[0]
            self.assertEqual("invalid", invalid.status)
            self.assertIn("invalid evidence state", invalid.reason)

            self.assertEqual(5, manifest["counts"]["input_records"])
            self.assertEqual(4, manifest["counts"]["normalized_records"])
            self.assertEqual(1, manifest["counts"]["invalid_records"])
            self.assertEqual(1, manifest["counts"]["insertion_observed_only"])
            self.assertEqual(1, manifest["counts"]["genotype_supported_mutant"])
            self.assertEqual(2, manifest["counts"]["strong_causal"])
            self.assertEqual(1, manifest["counts"]["requires_identity_bridge"])
            self.assertEqual(2, manifest["counts"]["route_b_function_transfer_ready"])
            self.assertFalse(manifest["cross_accession_identity_inferred"])
            self.assertEqual("none", manifest["rank_contribution"])

    def test_causal_row_without_publication_identifier_stays_blocked(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "tnt1.tsv"
            row = (
                "medtr.R108.gnmHiC_1.ann1.G6\tR108\tNF6\tFST6\tchr1:60\t"
                "supported\tsupported\tsupported\tnot_evaluated\tnot_evaluated\t"
                "root phenotype\t\t\t\thttps://medicago-mutant.okstate.edu/\tno paper id"
            )
            path.write_text(HEADER + "\n" + row + "\n", encoding="utf-8")
            result = normalize_tnt1_evidence(path, root / "out")
            ledger = pd.read_csv(result["ledger"], sep="\t", dtype=str).fillna("")
            functional = pd.read_csv(result["normalized_functional_evidence"], sep="\t", dtype=str).fillna("")

            self.assertEqual("genotype_supported_mutant", ledger.iloc[0].evidence_class)
            self.assertEqual("false", ledger.iloc[0].route_b_function_transfer_ready)
            self.assertEqual("causal_evidence_missing_publication_identifier", ledger.iloc[0].route_b_block_reason)
            self.assertTrue(functional.empty)

    def test_non_publication_source_id_stays_blocked(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "tnt1.tsv"
            row = (
                "medtr.R108.gnmHiC_1.ann1.G7\tR108\tNF7\tFST7\tchr1:70\t"
                "supported\tsupported\tsupported\tnot_evaluated\tnot_evaluated\t"
                "root phenotype\t\t\tTO:0000001\thttps://example.org/paper\tontology only"
            )
            path.write_text(HEADER + "\n" + row + "\n", encoding="utf-8")
            result = normalize_tnt1_evidence(path, root / "out")
            ledger = pd.read_csv(result["ledger"], sep="\t", dtype=str).fillna("")
            functional = pd.read_csv(result["normalized_functional_evidence"], sep="\t", dtype=str).fillna("")

            self.assertEqual("genotype_supported_mutant", ledger.iloc[0].evidence_class)
            self.assertEqual("false", ledger.iloc[0].route_b_function_transfer_ready)
            self.assertEqual("causal_evidence_missing_publication_identifier", ledger.iloc[0].route_b_block_reason)
            self.assertTrue(functional.empty)

    def test_malformed_publication_identifier_stays_blocked(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "tnt1.tsv"
            row = (
                "medtr.R108.gnmHiC_1.ann1.G7b\tR108\tNF7b\tFST7b\tchr1:75\t"
                "supported\tsupported\tsupported\tnot_evaluated\tnot_evaluated\t"
                "root phenotype\t\t\tPMID:NOT_A_PAPER\thttps://example.org/paper\tmalformed PMID"
            )
            path.write_text(HEADER + "\n" + row + "\n", encoding="utf-8")
            result = normalize_tnt1_evidence(path, root / "out")
            ledger = pd.read_csv(result["ledger"], sep="\t", dtype=str).fillna("")
            functional = pd.read_csv(result["normalized_functional_evidence"], sep="\t", dtype=str).fillna("")

            self.assertEqual("false", ledger.iloc[0].route_b_function_transfer_ready)
            self.assertEqual("causal_evidence_missing_publication_identifier", ledger.iloc[0].route_b_block_reason)
            self.assertTrue(functional.empty)

    def test_publication_identifier_requires_a_complete_token(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "tnt1.tsv"
            rows = []
            for suffix, source_ids in (
                ("pmid-suffix", "PMID:12345-typo"),
                ("pmcid-suffix", "PMCID:PMC12345.xyz"),
                ("doi-in-url", "https://example.org/10.1234/not-a-doi"),
            ):
                rows.append(
                    f"medtr.R108.gnmHiC_1.ann1.G8{suffix}\tR108\tNF8\tFST8\tchr1:80\t"
                    "supported\tsupported\tsupported\tnot_evaluated\tnot_evaluated\t"
                    f"root phenotype\t\t\t{source_ids}\thttps://example.org/paper\tmalformed token"
                )
            rows.append(
                "medtr.R108.gnmHiC_1.ann1.G8pmcid\tR108\tNF8\tFST8\tchr1:81\t"
                "supported\tsupported\tsupported\tnot_evaluated\tnot_evaluated\t"
                "root phenotype\t\t\tPMCID:PMC8122578\thttps://example.org/paper\tvalid PMCID"
            )
            path.write_text(HEADER + "\n" + "\n".join(rows) + "\n", encoding="utf-8")
            result = normalize_tnt1_evidence(path, root / "out")
            ledger = pd.read_csv(result["ledger"], sep="\t", dtype=str).fillna("")

            self.assertEqual(["false", "false", "false", "true"], ledger.route_b_function_transfer_ready.tolist())

    def test_non_medicago_gene_id_is_audited(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "tnt1.tsv"
            row = (
                "pissa.R108.gnm1.ann1.G8\tR108\tNF8\tFST8\tchr1:80\t"
                "supported\tsupported\tsupported\tnot_evaluated\tnot_evaluated\t"
                "root phenotype\t\t\tPMID:12346\thttps://example.org/paper\tnon-Medicago"
            )
            path.write_text(HEADER + "\n" + row + "\n", encoding="utf-8")
            result = normalize_tnt1_evidence(path, root / "out")
            audit = pd.read_csv(result["audit"], sep="\t", dtype=str).fillna("")

            self.assertEqual("invalid", audit.iloc[0].status)
            self.assertEqual("reported_gene_id_is_not_medicago_truncatula", audit.iloc[0].reason)
            self.assertEqual(0, result["counts"]["normalized_records"])

    def test_cli_wiring(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = self._input(root)
            out = root / "cli_out"
            self.assertEqual(0, main([
                "tnt1-evidence-normalize",
                "--evidence", str(source),
                "--out-dir", str(out),
                "--case-id", "pea_k1",
            ]))
            manifest = json.loads((out / "tnt1_evidence_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("pea_k1", manifest["case_id"])
            self.assertEqual(2, manifest["counts"]["route_b_function_transfer_ready"])


if __name__ == "__main__":
    unittest.main()
