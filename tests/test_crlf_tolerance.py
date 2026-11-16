"""CRLF tolerance for TSV inputs produced on Windows (issue #302).

Windows-authored TSV files (e.g. a batch manifest written by an editor that
uses CRLF) used to break parsing: readers stripped only ``\\n`` so the last
field carried a trailing ``\\r``, and ``parse_interval`` rejected the value.
Readers now strip ``\\r\\n`` and ``parse_interval`` tolerates surrounding
whitespace, so the same file parses identically on either platform.
"""

import tempfile
import unittest
from pathlib import Path

from candidate_qtl.core import iter_gff_rows, parse_interval, parse_seqid_map


class ParseIntervalCrlfTests(unittest.TestCase):
    def test_interval_with_trailing_carriage_return(self):
        self.assertEqual(parse_interval("Chr5:100-200\r"), ("Chr5", 100, 200))

    def test_interval_with_surrounding_whitespace(self):
        self.assertEqual(parse_interval("  Chr5:100-200  "), ("Chr5", 100, 200))

    def test_invalid_interval_still_rejected(self):
        with self.assertRaises(ValueError):
            parse_interval("not-an-interval")


class TsvReaderCrlfTests(unittest.TestCase):
    def test_seqid_map_with_crlf_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "seqid_map.tsv"
            path.write_bytes(b"Chr5\tscaffold_5\r\nChr3\tscaffold_3\r\n")
            self.assertEqual(parse_seqid_map(str(path)), {"Chr5": "scaffold_5", "Chr3": "scaffold_3"})

    def test_gff_rows_with_crlf_lines(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "genes.gff3"
            path.write_bytes(
                b"##gff-version 3\r\n"
                b"Chr5\tsoftmasked\tgene\t100\t200\t.\t+\t.\tID=gene1\r\n"
                b"Chr3\tsoftmasked\tgene\t300\t400\t.\t+\t.\tID=gene2\r\n"
            )
            rows = list(iter_gff_rows(path, "Chr5"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["attrs"], {"ID": "gene1"})

    def test_crlf_gzip_input(self):
        import gzip

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "seqid_map.tsv.gz"
            with gzip.open(path, "wb") as handle:
                handle.write(b"Chr5\tscaffold_5\r\nChr3\tscaffold_3\r\n")
            self.assertEqual(parse_seqid_map(str(path)),
                             {"Chr5": "scaffold_5", "Chr3": "scaffold_3"})

    def test_crlf_final_line_without_trailing_newline(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "seqid_map.tsv"
            path.write_bytes(b"Chr5\tscaffold_5\r\nChr3\tscaffold_3\r")
            self.assertEqual(parse_seqid_map(str(path)),
                             {"Chr5": "scaffold_5", "Chr3": "scaffold_3"})

    def test_evidence_best_hits_reader_tolerates_crlf(self):
        from candidate_qtl.evidence_io import best_hits

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "best_hits.tsv"
            path.write_bytes(
                b"G_1\tQ9LZS3\t85.0\t100\t110\t90\t1\t90\t1\t90\t1e-50\t200.0\r\n"
                b"G_2\tP0A8L9\t80.0\t90\t100\t80\t1\t80\t1\t80\t1e-40\t150.0\r\n"
            )
            grouped = best_hits(path)
            self.assertEqual(sorted(grouped), ["G_1", "G_2"])
            self.assertEqual(grouped["G_1"][0]["subject"], "Q9LZS3")


if __name__ == "__main__":
    unittest.main()
