from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from candidate_qtl.cli import main


class AnnotationSourceTests(unittest.TestCase):
    def test_stitle_and_word_boundary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = root / "candidates.tsv"
            candidates.write_text("gene_id\nG1\nG2\n", encoding="utf-8")
            blast = root / "hits.tsv"
            blast.write_text(
                "G1\tS1\t99\t10\t10\t10\t1\t10\t1\t10\t1e-9\t50\t0\tCLE peptide protein\n"
                "G2\tS2\t99\t10\t10\t10\t1\t10\t1\t10\t1e-9\t50\t0\tnuclear protein\n",
                encoding="utf-8",
            )
            out = root / "evidence.tsv"
            self.assertEqual(0, main(["annotation-keywords", "--candidates", str(candidates), "--blast", str(blast), "--use-stitle", "--keyword", "CLE", "--out", str(out)]))
            with out.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(["G1"], [row["gene_id"] for row in rows])


if __name__ == "__main__":
    unittest.main()
