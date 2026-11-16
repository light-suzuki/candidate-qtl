from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main


class ExpressionUnitTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.candidates = self.root / "candidates.tsv"
        self.candidates.write_text("gene_id\nGENE1\nGENE2\n", encoding="utf-8")

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_units_are_not_cross_divided(self) -> None:
        expression = self.root / "expr.tsv"
        expression.write_text(
            "gene_id\ttissue\tstage\tvalue\tunit\tis_target\tsource_id\n"
            "GENE1\tmeristem\ts1\t10\tTPM\ttrue\tmatrixA\n"
            "GENE1\tleaf\ts1\t2\tTPM\tfalse\tmatrixA\n"
            "GENE1\tmeristem\ts1\t100\tcount\ttrue\tmatrixB\n"
            "GENE1\tleaf\ts1\t50\tcount\tfalse\tmatrixB\n"
            "GENE2\tmeristem\ts1\t0.1\tTPM\ttrue\tmatrixA\n"
            "GENE2\tleaf\ts1\t0.05\tTPM\tfalse\tmatrixA\n",
            encoding="utf-8",
        )
        out = self.root / "ev"
        self.assertEqual(
            0,
            main(["evidence-build", "--candidates", str(self.candidates), "--expression", str(expression), "--out-dir", str(out)]),
        )
        summary = pd.read_csv(out / "expression_summary.tsv", sep="\t").fillna("")
        gene1 = summary[summary["gene_id"] == "GENE1"]
        # One summary row per unit, so TPM and count stay separate.
        self.assertEqual(set(gene1["unit"]), {"TPM", "count"})
        tpm = gene1[gene1["unit"] == "TPM"].iloc[0]
        count = gene1[gene1["unit"] == "count"].iloc[0]
        self.assertAlmostEqual(float(tpm["target_to_background_fold"]), 5.0, places=3)
        self.assertAlmostEqual(float(count["target_to_background_fold"]), 2.0, places=3)

        gene2 = summary[summary["gene_id"] == "GENE2"].iloc[0]
        self.assertEqual(gene2["non_expression"], "yes")
        self.assertEqual(gene2["expressed"], "no")

        evidence = pd.read_csv(out / "evidence.tsv", sep="\t").fillna("")
        gene1_codes = set(evidence[evidence["gene_id"] == "GENE1"]["evidence_code"])
        self.assertIn("expression_target_specific", gene1_codes)
        self.assertIn("expression_absent", set(evidence[evidence["gene_id"] == "GENE2"]["evidence_code"]))
        # Normalization unit is retained in the ledger description.
        self.assertTrue(evidence["description"].str.contains("unit=").any())

    def test_target_condition_filter(self) -> None:
        expression = self.root / "expr_cond.tsv"
        expression.write_text(
            "gene_id\ttissue\tcondition\tvalue\tunit\tis_target\n"
            "GENE1\tmeristem\tcontrol\t10\tTPM\ttrue\n"
            "GENE1\tmeristem\ttreatment\t1\tTPM\ttrue\n"
            "GENE1\tleaf\tcontrol\t2\tTPM\tfalse\n",
            encoding="utf-8",
        )
        out = self.root / "ev_cond"
        self.assertEqual(
            0,
            main(
                [
                    "evidence-build", "--candidates", str(self.candidates), "--expression", str(expression),
                    "--target-condition", "control", "--out-dir", str(out),
                ]
            ),
        )
        summary = pd.read_csv(out / "expression_summary.tsv", sep="\t").fillna("")
        gene1 = summary[summary["gene_id"] == "GENE1"].iloc[0]
        # Only the control target row (value 10) counts, not the treatment row.
        self.assertAlmostEqual(float(gene1["target_mean"]), 10.0, places=3)


class TraitOntologyTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_local_mapping_records_match_type_and_approval(self) -> None:
        terms = self.root / "terms.tsv"
        terms.write_text(
            "term\tterm_type\tsource_id\tapproved\n"
            "fasciation\tdirect_phenotype\tpaper1\ttrue\n"
            "band widening\tdirect_phenotype\tpaper1\ttrue\n"
            "novel unmapped trait\tmechanism_hypothesis\tpaper2\tfalse\n",
            encoding="utf-8",
        )
        ontology = self.root / "map.tsv"
        ontology.write_text(
            "term\tontology\tontology_id\tontology_name\tmatch_type\tconfidence\tapproved\n"
            "fasciation\tPTO\tTO:0000208\tfasciation\texact\thigh\ttrue\n"
            "band widening\tPTO\tTO:0000208\tfasciation\tllm_suggested\tlow\tfalse\n",
            encoding="utf-8",
        )
        out = self.root / "trait_ontology.tsv"
        self.assertEqual(
            0,
            main(["trait-ontology-map", "--trait-terms", str(terms), "--ontology-map", str(ontology), "--out", str(out)]),
        )
        result = pd.read_csv(out, sep="\t").fillna("")
        by_term = {row["original_term"]: row for _, row in result.iterrows()}

        self.assertEqual(by_term["fasciation"]["ontology_id"], "TO:0000208")
        self.assertEqual(by_term["fasciation"]["match_type"], "exact")
        self.assertEqual(by_term["fasciation"]["approval_status"], "approved")

        # A model-suggested synonym stays pending even if the term was approved.
        self.assertEqual(by_term["band widening"]["match_type"], "llm_suggested")
        self.assertEqual(by_term["band widening"]["approval_status"], "pending_review")

        # An unmapped term is retained, not dropped.
        self.assertEqual(by_term["novel unmapped trait"]["match_type"], "unmapped")
        self.assertEqual(by_term["novel unmapped trait"]["ontology_id"], "")

    def test_shipped_template_maps_cleanly(self) -> None:
        repo = Path(__file__).resolve().parents[1]
        terms = self.root / "terms.tsv"
        terms.write_text("term\tterm_type\tapproved\nfasciation\tdirect_phenotype\ttrue\nmeristem\torgan_stage\ttrue\n", encoding="utf-8")
        out = self.root / "mapped.tsv"
        self.assertEqual(
            0,
            main(
                [
                    "trait-ontology-map",
                    "--trait-terms", str(terms),
                    "--ontology-map", str(repo / "templates" / "trait_ontology_map.tsv"),
                    "--out", str(out),
                ]
            ),
        )
        result = pd.read_csv(out, sep="\t").fillna("")
        self.assertIn("TO:0000208", set(result["ontology_id"]))
        self.assertIn("PO:0009013", set(result["ontology_id"]))


if __name__ == "__main__":
    unittest.main()
