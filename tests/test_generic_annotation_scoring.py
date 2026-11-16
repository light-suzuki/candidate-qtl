from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main


class GenericAnnotationNotTraitEvidenceTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        (self.root / "candidates.tsv").write_text(
            "gene_id\tchrom\tstart\tend\tcds_length\tannotation_quality\trepresentative_transcript\n"
            "GeneGeneric\tchr1\t100\t200\t300\tcomplete_orf\tGeneGeneric\n"
            "GeneGenetic\tchr1\t300\t400\t300\tcomplete_orf\tGeneGenetic\n",
            encoding="utf-8",
        )
        # Generic best hit, generic domain, generic GO -- all trait-independent.
        (self.root / "blast.tsv").write_text(
            "GeneGeneric\tHitX\t80\t90\t100\t100\t1\t90\t1\t90\t1e-40\t180\tgeneric kinase\n"
            "GeneGenetic\tHitY\t80\t90\t100\t100\t1\t90\t1\t90\t1e-40\t180\tgeneric kinase\n",
            encoding="utf-8",
        )
        (self.root / "hmmer.tbl").write_text(
            "PF00069 - GeneGeneric - 1e-20 50 0 1e-20 50 0 1 1 1 1 1 1 1 1 protein_kinase\n",
            encoding="utf-8",
        )
        (self.root / "eggnog.tsv").write_text(
            "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\n"
            "GeneGeneric\tS\t1e-10\t60\t-\tViridiplantae\tS\tuncharacterized\t-\tGO:0003674\n",
            encoding="utf-8",
        )
        (self.root / "manual.tsv").write_text(
            "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\n"
            "GeneGenetic\tknown_qtg\tsupport\thigh\tqtg\tcuration\tdoi:10.0/x\tknown QTG in a related species\n",
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_generic_annotation_is_context_not_scored_support(self) -> None:
        evidence_dir = self.root / "evidence"
        self.assertEqual(
            0,
            main([
                "evidence-build", "--candidates", str(self.root / "candidates.tsv"),
                "--blast", str(self.root / "blast.tsv"), "--hmmer", str(self.root / "hmmer.tbl"),
                "--eggnog", str(self.root / "eggnog.tsv"), "--manual-evidence", str(self.root / "manual.tsv"),
                "--out-dir", str(evidence_dir),
            ]),
        )
        evidence = pd.read_csv(evidence_dir / "evidence.tsv", sep="\t").fillna("")
        generic = evidence[evidence["gene_id"] == "GeneGeneric"]
        # Separate, clearly-named context codes, none scored as support.
        codes = dict(zip(generic["evidence_code"], generic["direction"]))
        self.assertEqual(codes.get("homology_context"), "context")
        self.assertEqual(codes.get("domain_context"), "context")
        self.assertEqual(codes.get("go_annotation_context"), "context")
        self.assertNotIn("domain_supported", set(generic["evidence_code"]))
        self.assertNotIn("homology_supported", set(generic["evidence_code"]))

        rank_dir = self.root / "rank"
        self.assertEqual(
            0,
            main(["evidence-rank", "--candidates", str(self.root / "candidates.tsv"), "--evidence", str(evidence_dir / "evidence.tsv"), "--out-dir", str(rank_dir), "--top-n", "5"]),
        )
        ranked = pd.read_csv(rank_dir / "ranked_candidates.tsv", sep="\t").fillna("")
        by_gene = {row["gene_id"]: row for _, row in ranked.iterrows()}
        # A gene rich in generic annotation only scores nothing and is not Top.
        self.assertEqual(float(by_gene["GeneGeneric"]["dry_score"]), 0.0)
        self.assertNotEqual(by_gene["GeneGeneric"]["dry_group"], "top")
        # The gene with a trait-linked/genetic line is Top.
        self.assertEqual(by_gene["GeneGenetic"]["dry_group"], "top")

        contributions = pd.read_csv(rank_dir / "score_contributions.tsv", sep="\t").fillna("")
        self.assertTrue(contributions[contributions["gene_id"] == "GeneGeneric"].empty)


if __name__ == "__main__":
    unittest.main()
