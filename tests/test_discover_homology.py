from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.homology_qc import quality_class


GENOME = ">chr1\n" + "ATG" * 300 + "\n"
GFF = (
    "##gff-version 3\n"
    "chr1\tfixture\tgene\t100\t300\t.\t+\t.\tID=GeneH;Note=candidate protein\n"
    "chr1\tfixture\tmRNA\t100\t300\t.\t+\t.\tID=GeneH.1;Parent=GeneH\n"
    "chr1\tfixture\tCDS\t100\t300\t.\t+\t0\tParent=GeneH.1\n"
)


class HomologyQcClassTests(unittest.TestCase):
    def test_quality_classes(self) -> None:
        self.assertEqual(quality_class(80, 0.9, 0.9), "full_length_match")
        self.assertEqual(quality_class(60, 0.6, 0.4), "conserved_domain_only")
        self.assertEqual(quality_class(35, 0.55, 0.35), "short_or_partial_match")
        self.assertEqual(quality_class(80, float("nan"), float("nan")), "coverage_unknown")


class DiscoverHomologyTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        (self.root / "genome.fa").write_text(GENOME, encoding="utf-8")
        (self.root / "genes.gff3").write_text(GFF, encoding="utf-8")

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def _run(self, hits: dict) -> tuple[pd.DataFrame, Path]:
        out = self.root / f"out_{abs(hash(str(hits))) % 100000}"
        args = ["discover", "--fasta", str(self.root / "genome.fa"), "--gff3", str(self.root / "genes.gff3"),
                "--interval", "chr1:1-900", "--phenotype", "xyzzy", "--out-dir", str(out), "--no-auto-local-db", "--allow-external-api"]
        with patch("candidate_qtl.discover._best_ncbi_hits", return_value=(hits, "", True, False)):
            self.assertEqual(0, main(args))
        return pd.read_csv(out / "04_evidence.tsv", sep="\t").fillna(""), out

    def test_single_clear_hit_is_context_with_quality_class(self) -> None:
        evidence, out = self._run({"GeneH": [{"subject": "sp|A", "identity": 80.0, "query_coverage": 0.9, "subject_coverage": 0.9, "evalue": 1e-30, "bitscore": 200.0, "description": "clear ortholog"}]})
        homology = evidence[evidence.evidence_code == "homology_context"]
        self.assertEqual(len(homology), 1)
        self.assertEqual(homology.iloc[0]["direction"], "context")
        self.assertIn("full_length_match", homology.iloc[0]["description"])
        self.assertTrue(evidence[evidence.evidence_code == "paralog_ambiguous"].empty)
        # An audit table lists every hit, not just the best.
        audit = pd.read_csv(out / "03_annotation" / "homology_audit.tsv", sep="\t").fillna("")
        self.assertEqual(list(audit["subject"]), ["sp|A"])

    def test_near_tie_produces_paralog_ambiguous(self) -> None:
        evidence, out = self._run({"GeneH": [
            {"subject": "sp|A", "identity": 80.0, "query_coverage": 0.9, "subject_coverage": 0.9, "evalue": 1e-30, "bitscore": 200.0, "description": "member A"},
            {"subject": "sp|B", "identity": 78.0, "query_coverage": 0.9, "subject_coverage": 0.9, "evalue": 1e-29, "bitscore": 195.0, "description": "member B"},
        ]})
        para = evidence[evidence.evidence_code == "paralog_ambiguous"]
        self.assertEqual(len(para), 1)
        self.assertEqual(para.iloc[0]["direction"], "against")
        audit = pd.read_csv(out / "03_annotation" / "homology_audit.tsv", sep="\t").fillna("")
        self.assertEqual(len(audit), 2)

    def test_domain_only_hit_is_labeled(self) -> None:
        evidence, _ = self._run({"GeneH": [{"subject": "sp|D", "identity": 60.0, "query_coverage": 0.6, "subject_coverage": 0.4, "evalue": 1e-15, "bitscore": 120.0, "description": "shared domain"}]})
        homology = evidence[evidence.evidence_code == "homology_context"].iloc[0]
        self.assertIn("conserved_domain_only", homology["description"])


if __name__ == "__main__":
    unittest.main()
