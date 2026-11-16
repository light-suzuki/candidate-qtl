from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main


class ReferenceCandidateTableTests(unittest.TestCase):
    def test_keeps_mapped_and_unmapped_candidates_and_flags_ambiguity(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            candidates = root / "candidates.tsv"
            candidates.write_text("gene_id\tgene_name\tchrom\tstart\tend\nG1\tAlpha\tchr1\t10\t20\nG2\tBeta\tchr1\t30\t40\n", encoding="utf-8")
            panel = root / "panel.tsv"
            panel.write_text("reference_gene_id\treference_gene_name\treference_species\ttrait_relation\tevidence_level\tsource_id\tdiscovery_method\treview_status\nAT1\tKNOWN1\tArabidopsis\tmeristem phenotype\tdirect_model_genetic\tPMID:1\tllm_assisted\tpending_review\n", encoding="utf-8")
            one = root / "one.tsv"
            one.write_text("candidate_gene_id\treference_gene_id\trelationship\nG1\tAT1\tone_to_many\n", encoding="utf-8")
            two = root / "two.tsv"
            two.write_text("candidate_gene_id\treference_gene_id\trelationship\nG1\tAT1\tone_to_many\n", encoding="utf-8")
            annotation = root / "annotation.tsv"
            annotation.write_text("gene_id\tannotation\nG1\tLRR receptor kinase\nG2\tunknown protein\n", encoding="utf-8")
            out = root / "out"
            self.assertEqual(0, main(["reference-candidate-table", "--candidates", str(candidates), "--reference-genes", str(panel), "--orthology", f"orthofinder={one}", "--orthology", f"omamer={two}", "--annotation", str(annotation), "--out-dir", str(out)]))
            table = pd.read_csv(out / "reference_candidate_table.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual("KNOWN1", table.loc["G1", "reference_gene_name"])
            self.assertEqual("omamer;orthofinder", table.loc["G1", "mapping_providers"])
            self.assertEqual("2", table.loc["G1", "mapping_provider_count"])
            self.assertEqual("True", table.loc["G1", "review_required"])
            self.assertEqual("llm_assisted", table.loc["G1", "discovery_method"])
            self.assertEqual("pending_review", table.loc["G1", "review_status"])
            self.assertEqual("no_reference_panel_match", table.loc["G2", "mapping_status"])
            self.assertEqual("unknown protein", table.loc["G2", "annotation"])

    def test_rejects_reference_gene_without_source(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            (root / "candidates.tsv").write_text("gene_id\nG1\n", encoding="utf-8")
            (root / "panel.tsv").write_text("reference_gene_id\tsource_id\nAT1\t\n", encoding="utf-8")
            (root / "orth.tsv").write_text("candidate_gene_id\treference_gene_id\trelationship\nG1\tAT1\tone_to_one\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                main(["reference-candidate-table", "--candidates", str(root / "candidates.tsv"), "--reference-genes", str(root / "panel.tsv"), "--orthology", f"test={root / 'orth.tsv'}", "--out-dir", str(root / "out")])


if __name__ == "__main__":
    unittest.main()
