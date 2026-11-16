import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from candidate_qtl.variants import (
    alignment_metrics,
    allele_compare,
    classify_alignment_set,
    left_normalize_indel,
    parse_cs,
)


def _paf(qname, qlen, qs, qe, strand, tname, tlen, ts, te, nmatch, block, mapq, cs=None):
    fields = [qname, str(qlen), str(qs), str(qe), strand, tname, str(tlen), str(ts), str(te), str(nmatch), str(block), str(mapq)]
    if cs is not None:
        fields.append(f"cs:Z:{cs}")
    return fields


class AlignmentClassificationTests(unittest.TestCase):
    def test_single_simple_alignment_has_no_structural_event(self) -> None:
        aln = [_paf("alt", 100, 0, 100, "+", "ref", 100, 0, 100, 100, 100, 60, ":100")]
        result = classify_alignment_set(aln, 100, 100)
        self.assertEqual(result["structural_events"], [])
        self.assertFalse(result["ambiguous_primary_alignment"])

    def test_mixed_strand_is_inversion(self) -> None:
        aln = [
            _paf("alt", 100, 0, 100, "+", "ref", 100, 0, 100, 100, 100, 60, ":100"),
            _paf("alt", 100, 0, 100, "-", "ref", 100, 0, 100, 80, 100, 60),
        ]
        result = classify_alignment_set(aln, 100, 100)
        self.assertIn("inversion_suspected", result["structural_events"])

    def test_disjoint_query_segments_are_split(self) -> None:
        aln = [
            _paf("alt", 100, 0, 50, "+", "ref", 600, 0, 50, 50, 50, 60, ":50"),
            _paf("alt", 100, 50, 100, "+", "ref", 600, 550, 600, 50, 50, 60),
        ]
        result = classify_alignment_set(aln, 600, 100)
        self.assertIn("split_alignment", result["structural_events"])

    def test_same_query_two_loci_is_ambiguous_duplication(self) -> None:
        aln = [
            _paf("alt", 100, 0, 50, "+", "ref", 600, 0, 50, 50, 50, 60, ":50"),
            _paf("alt", 100, 0, 50, "+", "ref", 600, 550, 600, 50, 50, 60),
        ]
        result = classify_alignment_set(aln, 600, 100)
        self.assertTrue(result["ambiguous_primary_alignment"])
        self.assertIn("duplication_suspected", result["structural_events"])


class LeftNormalizeTests(unittest.TestCase):
    def test_homopolymer_deletion_shifts_left(self) -> None:
        pos, ref, alt = left_normalize_indel("AAAAAAAAAA", 5, "A", "-")
        self.assertEqual((pos, ref, alt), (1, "A", "-"))

    def test_polyA_insertion_shifts_left(self) -> None:
        pos, ref, alt = left_normalize_indel("TAAAAG", 5, "-", "A")
        self.assertEqual((pos, ref, alt), (2, "-", "A"))

    def test_snv_is_unchanged(self) -> None:
        self.assertEqual(left_normalize_indel("ACGT", 2, "C", "G"), (2, "C", "G"))


class VariantParsingTests(unittest.TestCase):
    def test_long_cs_parses_snv_and_inframe_insertion(self) -> None:
        rows = parse_cs("=ACGT*ag=TT+GGAACCTTT=CC", reference_start=10, query_start=20)
        self.assertEqual(rows[0]["variant_type"], "snv")
        self.assertEqual(rows[0]["reference_pos_1based"], 15)
        self.assertEqual(rows[1]["variant_class"], "inframe_insertion")
        self.assertEqual(rows[1]["length"], 9)

    def test_frameshift_deletion_is_classified(self) -> None:
        rows = parse_cs(":5-AC:3")
        self.assertEqual(rows[0]["variant_class"], "frameshift_deletion")
        self.assertEqual(rows[0]["reference_pos_1based"], 6)

    def test_alignment_metrics_reject_partial_comparison(self) -> None:
        fields = ["alt", "50", "0", "50", "+", "ref", "100", "0", "50", "50", "50", "60"]
        metrics = alignment_metrics(fields, reference_length=100, alternate_length=50)
        self.assertEqual(0.5, metrics["reference_coverage"])
        self.assertEqual(1.0, metrics["alternate_coverage"])
        self.assertFalse(metrics["comparison_complete"])

    def test_partial_alignment_suppresses_genomic_marker_export(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            reference = root / "reference.fa"
            alternate = root / "alternate.fa"
            candidates = root / "candidates.tsv"
            reference.write_text(">GeneA\n" + "A" * 100 + "\n", encoding="utf-8")
            alternate.write_text(">GeneB\n" + "A" * 50 + "\n", encoding="utf-8")
            candidates.write_text(
                "gene_id\tchrom\tstart\tend\tstrand\nGeneA\tchr1\t1001\t1100\t+\n",
                encoding="utf-8",
            )
            paf = "GeneB\t50\t0\t50\t+\tGeneA\t100\t0\t50\t50\t50\t60\tcs:Z::50\n"
            with patch("candidate_qtl.variants.shutil.which", return_value="/fake/minimap2"), patch(
                "candidate_qtl.core.subprocess.run",
                return_value=SimpleNamespace(stdout=paf),
            ):
                result = allele_compare(
                    str(reference),
                    str(alternate),
                    str(root / "out"),
                    reference_candidates_arg=str(candidates),
                )
            self.assertEqual("partial_alignment", result["comparison_status"])
            self.assertIsNone(result["marker_candidates"])
            self.assertIn("partial_alignment", result["marker_suppression_reasons"])
            saved = json.loads((root / "out" / "allele_compare.json").read_text(encoding="utf-8"))
            self.assertEqual(0.5, saved["reference_coverage"])

    def test_complete_genomic_alignment_exports_forward_coordinate(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            reference = root / "reference.fa"
            alternate = root / "alternate.fa"
            candidates = root / "candidates.tsv"
            reference.write_text(">GeneA\nAAAAAAAAAA\n", encoding="utf-8")
            alternate.write_text(">GeneB\nAAAACAAAAA\n", encoding="utf-8")
            candidates.write_text(
                "gene_id\tchrom\tstart\tend\tstrand\nGeneA\tchr1\t101\t110\t+\n",
                encoding="utf-8",
            )
            paf = "GeneB\t10\t0\t10\t+\tGeneA\t10\t0\t10\t9\t10\t60\tcs:Z::4*ac:5\n"
            with patch("candidate_qtl.variants.shutil.which", return_value="/fake/minimap2"), patch(
                "candidate_qtl.core.subprocess.run",
                return_value=SimpleNamespace(stdout=paf),
            ):
                result = allele_compare(
                    str(reference),
                    str(alternate),
                    str(root / "out"),
                    reference_candidates_arg=str(candidates),
                )
            self.assertEqual("complete_alignment", result["comparison_status"])
            marker_text = Path(result["marker_candidates"]).read_text(encoding="utf-8")
            self.assertIn("chr1\t105\tA\tC", marker_text)


class StructuralAlleleCompareTests(unittest.TestCase):
    def _run(self, root: Path, ref_seq: str, alt_seq: str, paf: str, strand: str = "+", **kwargs):
        reference = root / "reference.fa"
        alternate = root / "alternate.fa"
        candidates = root / "candidates.tsv"
        reference.write_text(f">GeneA\n{ref_seq}\n", encoding="utf-8")
        alternate.write_text(f">GeneB\n{alt_seq}\n", encoding="utf-8")
        candidates.write_text(f"gene_id\tchrom\tstart\tend\tstrand\nGeneA\tchr1\t101\t{100 + len(ref_seq)}\t{strand}\n", encoding="utf-8")
        with patch("candidate_qtl.variants.shutil.which", return_value="/fake/minimap2"), patch(
            "candidate_qtl.core.subprocess.run", return_value=SimpleNamespace(stdout=paf)
        ):
            return allele_compare(str(reference), str(alternate), str(root / "out"), reference_candidates_arg=str(candidates), **kwargs)

    def test_inversion_suppresses_markers(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            paf = (
                "GeneB\t10\t0\t10\t+\tGeneA\t10\t0\t10\t10\t10\t60\tcs:Z::10\n"
                "GeneB\t10\t0\t10\t-\tGeneA\t10\t0\t10\t8\t10\t60\n"
            )
            result = self._run(root, "AAAAAAAAAA", "AAAAAAAAAA", paf)
            self.assertIn("inversion_suspected", result["structural_events"])
            self.assertIsNone(result["marker_candidates"])
            self.assertIn("inversion_suspected", result["marker_suppression_reasons"])

    def test_minus_strand_gene_emits_forward_genomic_allele(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            paf = "GeneB\t10\t0\t10\t+\tGeneA\t10\t0\t10\t9\t10\t60\tcs:Z::4*ac:5\n"
            result = self._run(root, "AAAAAAAAAA", "AAAACAAAAA", paf, strand="-")
            marker_text = Path(result["marker_candidates"]).read_text(encoding="utf-8")
            # SNV at gene offset 5 on a minus-strand gene: genomic 110-5+1=106, alleles reverse-complemented.
            self.assertIn("chr1\t106\tT\tG", marker_text)

    def test_indel_marker_is_left_normalized(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            paf = "GeneB\t9\t0\t9\t+\tGeneA\t10\t0\t10\t9\t10\t60\tcs:Z::4-A:5\n"
            result = self._run(root, "AAAAAAAAAA", "AAAAAAAAA", paf)
            marker_text = Path(result["marker_candidates"]).read_text(encoding="utf-8")
            lines = [line.split("\t") for line in marker_text.strip().splitlines()]
            header, row = lines[0], lines[1]
            record = dict(zip(header, row))
            # The homopolymer deletion normalizes to the leftmost position (101).
            self.assertEqual(record["pos"], "101")
            self.assertEqual(record["left_normalized"], "yes")

    def test_sequence_type_cds_keeps_frame_labels(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            paf = "GeneB\t9\t0\t9\t+\tGeneA\t10\t0\t10\t9\t10\t60\tcs:Z::4-A:5\n"
            result = self._run(root, "AAAAAAAAAA", "AAAAAAAAA", paf, sequence_type="cds")
            self.assertEqual(result["sequence_context"], "cds")
            variants = Path(result["variants"]).read_text(encoding="utf-8")
            self.assertIn("frameshift_deletion", variants)


if __name__ == "__main__":
    unittest.main()
