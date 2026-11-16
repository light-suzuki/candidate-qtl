from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl.cli import main


# A clean short complete ORF (M + K*10 + stop = 11 aa) at chr1:101-136.
SHORT_ORF = "ATG" + "AAA" * 10 + "TAA"
GENOME = "C" * 100 + SHORT_ORF + "C" * (400 - 100 - len(SHORT_ORF))
GFF = (
    "##gff-version 3\n"
    "chr1\tfixture\tgene\t101\t136\t.\t+\t.\tID=DarkGene;Note=uncharacterized protein\n"
    "chr1\tfixture\tmRNA\t101\t136\t.\t+\t.\tID=DarkGene.1;Parent=DarkGene\n"
    "chr1\tfixture\tCDS\t101\t136\t.\t+\t0\tParent=DarkGene.1\n"
)


class DiscoverDarkResolutionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.fasta = self.root / "genome.fa"; self.fasta.write_text(">chr1\n" + GENOME + "\n", encoding="utf-8")
        self.gff = self.root / "genes.gff3"; self.gff.write_text(GFF, encoding="utf-8")

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def _args(self, out: Path) -> list[str]:
        return ["discover", "--fasta", str(self.fasta), "--gff3", str(self.gff), "--interval", "chr1:91-350",
                "--phenotype", "meristem enlargement", "--out-dir", str(out)]

    def _scope(self, out: Path) -> tuple[str, str]:
        flags = pd.read_csv(out / "02_candidate_flags.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
        return flags.loc["DarkGene", "candidate_scope"], flags.loc["DarkGene", "dark_state"]

    def test_offline_search_not_run_is_provisional(self) -> None:
        out = self.root / "offline"
        self.assertEqual(0, main(self._args(out) + ["--offline", "--no-auto-local-db"]))
        scope, state = self._scope(out)
        self.assertEqual(scope, "dark_pending")
        self.assertEqual(state, "homology_not_searched")
        manifest = json.loads((out / "discover_manifest.json").read_text())
        self.assertFalse(manifest["homology_searched"])

    def test_searched_with_no_hit_is_confirmed_dark(self) -> None:
        out = self.root / "nohit"
        with patch("candidate_qtl.discover._best_ncbi_hits", return_value=({}, "", True, False)):
            self.assertEqual(0, main(self._args(out) + ["--no-auto-local-db", "--allow-external-api"]))
        scope, state = self._scope(out)
        self.assertEqual(scope, "dark")
        self.assertEqual(state, "confirmed_no_homology")

    def test_short_orf_with_strong_hit_leaves_dark(self) -> None:
        out = self.root / "hit"
        strong = {"DarkGene": [{"subject": "sp|KIN", "identity": 80.0, "query_coverage": 1.0, "subject_coverage": 0.9, "evalue": 1e-20, "bitscore": 90.0, "description": "protein kinase"}]}
        with patch("candidate_qtl.discover._best_ncbi_hits", return_value=(strong, "", True, False)):
            self.assertEqual(0, main(self._args(out) + ["--no-auto-local-db", "--allow-external-api"]))
        scope, state = self._scope(out)
        self.assertEqual(scope, "positional")
        self.assertEqual(state, "homology_found_reclassified")

    def test_search_failure_is_distinct_state(self) -> None:
        out = self.root / "failed"
        with patch("candidate_qtl.discover._best_local_hits", side_effect=RuntimeError("blastp exploded")):
            self.assertEqual(0, main(self._args(out) + ["--protein-db", str(self.root / "fake_db")]))
        scope, state = self._scope(out)
        self.assertEqual(scope, "dark_pending")
        self.assertEqual(state, "homology_search_failed")
        self.assertIn("homology_search_failed", (out / "review_required.tsv").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
