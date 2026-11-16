"""Real BLAST+ integration smoke: the local_blastp path end to end.

Skipped unless BLAST+ (blastp + makeblastdb) is on PATH, so the normal unit
suite stays hermetic; the CI `integration-blast` job installs ncbi-blast+ and
runs it for real.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main


@unittest.skipUnless(shutil.which("blastp") and shutil.which("makeblastdb"), "BLAST+ (blastp/makeblastdb) not installed")
class BlastIntegrationTests(unittest.TestCase):
    def test_local_blastp_path_produces_homology_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # A varied 79-aa ORF (poly-residue runs get composition-filtered by BLAST).
            orf = ("ATGAAAACTGCTTATATTGCTAAACAACGTCAAATTTCTTTTGTTAAATCTCATTTTTCTCGTCAACTGGAAGAACGT"
                   "CTGGGTCTGATTGAAGTTCAAGCTCCGATTCTGTCTCGTGTTGGTGATGGTACTCAAGATAATCTGTCTGGTGCTGAA"
                   "AAAGCTGTTCAAGTTAAAGTTAAAGCTCTGCCGGATGCTCAATTTGAAGTTGTTCATTCTCTGGCTAAATGGAAACGTTAA")
            (root / "genome.fa").write_text(">chr1\n" + orf + "C" * 200 + "\n", encoding="utf-8")
            (root / "genes.gff3").write_text(
                "##gff-version 3\n"
                "chr1\tx\tgene\t1\t237\t.\t+\t.\tID=G1;Note=test protein\n"
                "chr1\tx\tmRNA\t1\t237\t.\t+\t.\tID=G1.1;Parent=G1\n"
                "chr1\tx\tCDS\t1\t237\t.\t+\t0\tParent=G1.1\n",
                encoding="utf-8",
            )
            extract_out = root / "extract"
            self.assertEqual(0, main(["extract", "--fasta", str(root / "genome.fa"), "--gff3", str(root / "genes.gff3"), "--interval", "chr1:1-237", "--out-dir", str(extract_out)]))
            proteins = extract_out / "candidate_proteins.fa"
            self.assertTrue(proteins.exists())

            db = root / "protein_db"
            subprocess.run(["makeblastdb", "-in", str(proteins), "-dbtype", "prot", "-out", str(db)], check=True, capture_output=True, text=True)

            discover_out = root / "discover"
            self.assertEqual(
                0,
                main(["discover", "--fasta", str(root / "genome.fa"), "--gff3", str(root / "genes.gff3"), "--interval", "chr1:1-237",
                      "--phenotype", "test trait", "--protein-db", str(db), "--out-dir", str(discover_out), "--no-auto-local-db"]),
            )
            evidence = pd.read_csv(discover_out / "04_evidence.tsv", sep="\t").fillna("")
            self.assertIn("homology_context", set(evidence["evidence_code"]))
            manifest = json.loads((discover_out / "discover_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["annotation_provider"], "local_blastp")
            # The audit table from the real search exists.
            self.assertTrue((discover_out / "03_annotation" / "homology_audit.tsv").exists())


if __name__ == "__main__":
    unittest.main()
