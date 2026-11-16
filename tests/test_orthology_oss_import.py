from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main


class OrthoFinderImportTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        (self.root / "candidates.tsv").write_text("gene_id\nG1\nG2\nG3\nG4\nG5\n", encoding="utf-8")
        (self.root / "Orthogroups.tsv").write_text(
            "Orthogroup\tpea\tarabidopsis\trice\n"
            "OG0000001\tG1\tAT1\tOS1\n"
            "OG0000002\tG2\tAT2, AT3\tOS2\n"
            "OG0000003\tG3, G3b\tAT4\tOS3\n"
            "OG0000004\tG4, G4b\tAT5, AT6\t\n"
            "OG0000005\tG5\t\t\n",
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_relationships_and_paralogs(self) -> None:
        out = self.root / "of_evidence.tsv"
        self.assertEqual(
            0,
            main(["orthofinder-import", "--orthogroups", str(self.root / "Orthogroups.tsv"), "--candidates", str(self.root / "candidates.tsv"), "--focal-species", "pea", "--out", str(out)]),
        )
        summary = pd.read_csv(out.with_name("of_evidence_orthofinder_summary.tsv"), sep="\t").fillna("")
        rel = {row["gene_id"]: row["relationship"] for _, row in summary.iterrows()}
        self.assertEqual(rel["G1"], "one_to_one")
        self.assertEqual(rel["G2"], "one_to_many")
        self.assertEqual(rel["G3"], "many_to_one")
        self.assertEqual(rel["G4"], "many_to_many")
        self.assertEqual(rel["G5"], "lineage_specific_or_focal_only")

        evidence = pd.read_csv(out, sep="\t").fillna("")
        g1 = evidence[evidence["gene_id"] == "G1"]
        self.assertIn("orthogroup_supported", set(g1["evidence_code"]))
        self.assertEqual(g1[g1["evidence_code"] == "orthogroup_supported"].iloc[0]["confidence"], "high")
        # In-species duplicates raise paralog ambiguity.
        self.assertIn("paralog_ambiguous", set(evidence[evidence["gene_id"] == "G3"]["evidence_code"]))
        # A focal-only orthogroup yields no cross-species orthogroup support.
        self.assertTrue(evidence[evidence["gene_id"] == "G5"].empty)


class OmamerImportTests(unittest.TestCase):
    def test_hog_assignment(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            (root / "candidates.tsv").write_text("gene_id\nG1\nG2\n", encoding="utf-8")
            (root / "omamer.tsv").write_text("qseqid\thogid\tfamily_p\nG1\tHOG:0042\t0.9\nG2\tna\t0.0\n", encoding="utf-8")
            out = root / "omamer_evidence.tsv"
            self.assertEqual(0, main(["omamer-import", "--omamer", str(root / "omamer.tsv"), "--candidates", str(root / "candidates.tsv"), "--out", str(out)]))
            evidence = pd.read_csv(out, sep="\t").fillna("")
            g1 = evidence[evidence["gene_id"] == "G1"].iloc[0]
            self.assertEqual(g1["evidence_code"], "orthogroup_supported")
            self.assertEqual(g1["value"], "HOG:0042")
            self.assertTrue(evidence[evidence["gene_id"] == "G2"].empty)


class SyntenyImportTests(unittest.TestCase):
    def test_anchors_format(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            (root / "candidates.tsv").write_text("gene_id\nG1\nG2\n", encoding="utf-8")
            (root / "pairs.anchors").write_text("###\nG1\tSbGene1\t100\nG2\tSbGene2\t95\n###\n", encoding="utf-8")
            out = root / "syn_evidence.tsv"
            self.assertEqual(0, main(["synteny-import", "--anchors", str(root / "pairs.anchors"), "--candidates", str(root / "candidates.tsv"), "--out", str(out)]))
            evidence = pd.read_csv(out, sep="\t").fillna("")
            g1 = evidence[evidence["gene_id"] == "G1"].iloc[0]
            self.assertEqual(g1["evidence_code"], "synteny_supported")
            self.assertEqual(g1["value"], "SbGene1")

    def test_collinearity_format(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            (root / "candidates.tsv").write_text("gene_id\nG1\n", encoding="utf-8")
            (root / "x.collinearity").write_text("## Alignment 0:\n  0-  0:\tG1\tSbGene1\t  1e-50\n", encoding="utf-8")
            out = root / "syn2.tsv"
            self.assertEqual(0, main(["synteny-import", "--anchors", str(root / "x.collinearity"), "--candidates", str(root / "candidates.tsv"), "--out", str(out), "--format", "collinearity"]))
            evidence = pd.read_csv(out, sep="\t").fillna("")
            self.assertEqual(evidence.iloc[0]["value"], "SbGene1")


if __name__ == "__main__":
    unittest.main()
