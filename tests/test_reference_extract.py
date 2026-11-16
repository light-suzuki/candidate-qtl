from __future__ import annotations

import csv
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

from Bio import SeqIO

from candidate_qtl.cli import main
from candidate_qtl.core import write_fasta


class ReferencePhaseExtractionTests(unittest.TestCase):
    def test_duplicate_parent_ids_are_sorted_unique_and_exclude_blank_and_cds_ids(self):
        from candidate_qtl.reference_extract import _normalize_features

        rows = []
        for kind, identifier, parent in (
            ("gene", "z", ""), ("gene", "z", ""), ("gene", "z", ""),
            ("gene", "g", ""), ("mRNA", "a", "g"), ("transcript", "a", "g"),
            ("ncRNA_gene", "n", ""), ("ncRNA_gene", "n", ""),
            ("pseudogene", "p", ""), ("pseudogene", "p", ""),
            ("gene", "", ""), ("gene", "   ", ""),
            ("CDS", "segment", "a"), ("CDS", "segment", "a"),
        ):
            rows.append({"seqid": "chr1", "type": kind, "start": 1, "end": 9,
                         "strand": "+", "phase": "0" if kind == "CDS" else ".",
                         "attrs": {"ID": identifier, "Parent": parent}})
        normalized, qc = _normalize_features(rows)
        self.assertEqual(qc["duplicate_parent_ids"], ["a", "n", "p", "z"])
        self.assertEqual(normalized, rows)
        self.assertEqual(_normalize_features([])[1]["duplicate_parent_ids"], [])

    def readonly_inputs(self, root):
        import pysam
        fasta = root / "genome.fa"
        fasta.write_text(">chr1\nATGAAATAA\n", encoding="utf-8")
        gff = root / "genes.gff3"
        gff.write_text("chr1\ttest\tgene\t1\t9\t.\t+\t.\tID=G1\n"
                       "chr1\ttest\tCDS\t1\t9\t.\t+\t0\tParent=G1\n", encoding="utf-8")
        pysam.faidx(str(fasta))
        return fasta, gff, Path(str(fasta) + ".fai")

    def test_readonly_extraction_rejects_output_overlap_and_alias_before_writing(self):
        from candidate_qtl.reference_extract import extract_candidates
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            fasta, gff, index = self.readonly_inputs(root)
            before = {p: p.read_bytes() for p in (fasta, gff, index)}
            with self.assertRaisesRegex(ValueError, "separate"):
                extract_candidates(str(fasta), str(gff), "chr1:1-9", str(root), read_only_reference=True)
            out = root / "derived"
            out.mkdir()
            (out / "candidate_genomic.fa").hardlink_to(fasta)
            with self.assertRaisesRegex(ValueError, "fresh output"):
                extract_candidates(str(fasta), str(gff), "chr1:1-9", str(out), read_only_reference=True)
            self.assertEqual(before, {p: p.read_bytes() for p in before})
            self.assertFalse((out / "candidates.tsv").exists())

    def test_readonly_invalid_existing_index_does_not_rebuild(self):
        from candidate_qtl.reference_extract import extract_candidates
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            fasta, gff, index = self.readonly_inputs(root)
            index.write_text("not a valid index\n", encoding="utf-8")
            before = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in (fasta, gff, index)}
            with patch("pysam.faidx", side_effect=AssertionError("shared index rebuild forbidden")):
                with self.assertRaises((OSError, ValueError)):
                    extract_candidates(str(fasta), str(gff), "chr1:1-9", str(root / "out"), read_only_reference=True)
            self.assertEqual(before, {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in before})

    def test_readonly_reference_mutation_during_extraction_is_not_sealed(self):
        from candidate_qtl import reference_extract
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            fasta, gff, _ = self.readonly_inputs(root)
            fingerprint = reference_extract.file_fingerprint
            def mutate(path, full_hash):
                if path == gff:
                    path.write_text(path.read_text() + "# concurrent update\n", encoding="utf-8")
                return fingerprint(path, full_hash)
            with patch.object(reference_extract, "file_fingerprint", side_effect=mutate):
                with self.assertRaisesRegex(ValueError, "changed during extraction"):
                    reference_extract.extract_candidates(str(fasta), str(gff), "chr1:1-9", str(root / "out"),
                                                         read_only_reference=True, full_hash=True)
            self.assertFalse((root / "out/extract_manifest.json").exists())

    def test_readonly_stale_fai_rejected_before_output_and_without_repair(self):
        import pysam
        from candidate_qtl.reference_extract import extract_candidates
        for replacement in (b">chr2\nAAAAAAAAAA\n>chr1\nCCCCCCCCCC\n",
                            b">chr1\nAAAAA\nAAAAA\n>chr2\nCCCCCCCCCC\n",
                            b">chr1\nAAAAAAAAAAA\n>chr2\nCCCCCCCCCC\n"):
            with self.subTest(replacement=replacement), tempfile.TemporaryDirectory() as tempdir:
                root = Path(tempdir)
                fasta, gff, index = self.readonly_inputs(root)
                fasta.write_bytes(b">chr1\nAAAAAAAAAA\n>chr2\nCCCCCCCCCC\n")
                pysam.faidx(str(fasta))
                fasta.write_bytes(replacement)
                before = {p: (p.read_bytes(), p.stat().st_mtime_ns, p.stat().st_ctime_ns) for p in (fasta, gff, index)}
                out = root / "out"
                with self.assertRaisesRegex(ValueError, "FASTA/.fai layout mismatch"):
                    extract_candidates(str(fasta), str(gff), "chr1:1-9", str(out), read_only_reference=True)
                self.assertFalse(out.exists())
                self.assertEqual(before, {p: (p.read_bytes(), p.stat().st_mtime_ns, p.stat().st_ctime_ns) for p in before})

    def test_readonly_fai_validation_keeps_seqid_mapping_and_crlf_sequences(self):
        import pysam
        from candidate_qtl.reference_extract import extract_candidates
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            fasta, gff, index = self.readonly_inputs(root)
            fasta.write_bytes(b">other description\r\nCCCCCCCCC\r\n>assembly1\r\nATGAA\r\nATAA")
            pysam.faidx(str(fasta))
            before = {p: (p.read_bytes(), p.stat().st_mtime_ns, p.stat().st_ctime_ns) for p in (fasta, gff, index)}
            manifest = extract_candidates(str(fasta), str(gff), "chr1:1-9", str(root / "out"),
                                          seqid_map_arg="chr1=assembly1", full_hash=True, read_only_reference=True)
            self.assertEqual(">G1\nATGAAATAA\n", (root / "out/candidate_genomic.fa").read_text())
            self.assertEqual("independent_faidx_layout", manifest["fasta_index_validation"]["method"])
            self.assertEqual(2, manifest["fasta_index_validation"]["contig_count"])
            self.assertEqual(before, {p: (p.read_bytes(), p.stat().st_mtime_ns, p.stat().st_ctime_ns) for p in before})

    def test_fasta_writer_rejects_unstable_ids(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            target = Path(tempdir) / "out.fa"
            for record_id in ("", "gene with space", "gene\nnext", "gene\x7fnext", "gene>next"):
                with self.subTest(record_id=record_id), self.assertRaisesRegex(ValueError, "FASTA record ID"):
                    write_fasta(target, [(record_id, "MAAA")])
            with self.assertRaisesRegex(ValueError, "duplicate FASTA record ID"):
                write_fasta(target, [("gene:A", "MAAA"), ("gene:A", "MBBB")])

    def test_fasta_ids_round_trip_without_sanitization_collision(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            fasta = root / "genome.fa"
            fasta.write_text(">chr1\nATGAAATAA" + "N" * 10 + "ATGAAATAA\n", encoding="utf-8")
            gff = root / "genes.gff3"
            lines = ["##gff-version 3"]
            for gene, start, end in (("gene:A", 1, 9), ("gene_A", 20, 28)):
                lines.extend((
                    f"chr1\ttest\tgene\t{start}\t{end}\t.\t+\t.\tID={gene}",
                    f"chr1\ttest\tmRNA\t{start}\t{end}\t.\t+\t.\tID={gene}.t1;Parent={gene}",
                    f"chr1\ttest\tCDS\t{start}\t{end}\t.\t+\t0\tParent={gene}.t1",
                ))
            gff.write_text("\n".join(lines) + "\n", encoding="utf-8")
            out = root / "out"
            self.assertEqual(0, main(["extract", "--fasta", str(fasta), "--gff3", str(gff),
                                      "--interval", "chr1:1-28", "--out-dir", str(out)]))
            with (out / "candidates.tsv").open(encoding="utf-8", newline="") as handle:
                roster = {row["gene_id"] for row in csv.DictReader(handle, delimiter="\t")}
            with (out / "id_map.tsv").open(encoding="utf-8", newline="") as handle:
                mapped = {row["sequence_record_id"] for row in csv.DictReader(handle, delimiter="\t")}
            with (out / "candidate_proteins.fa").open(encoding="utf-8") as handle:
                protein_ids = {record.id for record in SeqIO.parse(handle, "fasta")}
            self.assertEqual({"gene:A", "gene_A"}, roster)
            self.assertEqual(roster, mapped)
            self.assertEqual(roster, protein_ids)

    def _run_extract(self, sequence: str, gff_text: str, gene_id: str) -> tuple[dict[str, str], dict[str, str], str]:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            fasta = root / "genome.fa"
            fasta.write_text(f">chr1\n{sequence}\n", encoding="utf-8")
            gff = root / "genes.gff3"
            gff.write_text(gff_text, encoding="utf-8")
            out = root / "out"
            self.assertEqual(
                0,
                main(
                    [
                        "extract",
                        "--fasta",
                        str(fasta),
                        "--gff3",
                        str(gff),
                        "--interval",
                        f"chr1:1-{len(sequence)}",
                        "--out-dir",
                        str(out),
                    ]
                ),
            )
            with (out / "candidates.tsv").open(encoding="utf-8") as handle:
                candidate = next(csv.DictReader(handle, delimiter="\t"))
            with (out / "transcripts.tsv").open(encoding="utf-8") as handle:
                transcript = next(csv.DictReader(handle, delimiter="\t"))
            protein_records = list(SeqIO.parse(out / "candidate_proteins.fa", "fasta"))
            protein = str(protein_records[0].seq) if protein_records else ""
            with (out / "id_map.tsv").open(encoding="utf-8") as handle:
                id_rows = list(csv.DictReader(handle, delimiter="\t"))
            representative = next(row for row in id_rows if row["is_representative"] == "True")
            self.assertEqual(gene_id if protein_records else "", representative["sequence_record_id"])
            self.assertEqual(representative["transcript_id"] if protein_records else "", representative["protein_id"])
            self.assertEqual(gene_id, candidate["gene_id"])
            return candidate, transcript, protein

    def test_plus_strand_internal_phase_is_validated_not_trimmed(self) -> None:
        candidate, transcript, protein = self._run_extract(
            "ATGAAAAAAATAG",
            "##gff-version 3\n"
            "chr1\tx\tgene\t1\t13\t.\t+\t.\tID=GenePlus\n"
            "chr1\tx\tmRNA\t1\t13\t.\t+\t.\tID=TxPlus;Parent=GenePlus\n"
            "chr1\tx\tCDS\t1\t5\t.\t+\t0\tParent=TxPlus\n"
            "chr1\tx\tCDS\t10\t13\t.\t+\t1\tParent=TxPlus\n",
            "GenePlus",
        )
        self.assertEqual("valid", transcript["phase_status"])
        self.assertEqual("0,1", transcript["reported_phases"])
        self.assertEqual("0,1", transcript["expected_phases"])
        self.assertEqual("complete_orf", candidate["annotation_quality"])
        self.assertEqual("MK", protein)

    def test_minus_strand_phase_uses_transcript_order(self) -> None:
        candidate, transcript, protein = self._run_extract(
            "CTATAAAAATTCAT",
            "##gff-version 3\n"
            "chr1\tx\tgene\t1\t14\t.\t-\t.\tID=GeneMinus\n"
            "chr1\tx\tmRNA\t1\t14\t.\t-\t.\tID=TxMinus;Parent=GeneMinus\n"
            "chr1\tx\tCDS\t10\t14\t.\t-\t0\tParent=TxMinus\n"
            "chr1\tx\tCDS\t1\t4\t.\t-\t1\tParent=TxMinus\n",
            "GeneMinus",
        )
        self.assertEqual("valid", transcript["phase_status"])
        self.assertEqual("0,1", transcript["reported_phases"])
        self.assertEqual("complete_orf", candidate["annotation_quality"])
        self.assertEqual("MK", protein)

    def test_inconsistent_phase_suppresses_protein_translation(self) -> None:
        candidate, transcript, protein = self._run_extract(
            "ATGAAAAAAATAG",
            "##gff-version 3\n"
            "chr1\tx\tgene\t1\t13\t.\t+\t.\tID=GeneBad\n"
            "chr1\tx\tmRNA\t1\t13\t.\t+\t.\tID=TxBad;Parent=GeneBad\n"
            "chr1\tx\tCDS\t1\t5\t.\t+\t0\tParent=TxBad\n"
            "chr1\tx\tCDS\t10\t13\t.\t+\t0\tParent=TxBad\n",
            "GeneBad",
        )
        self.assertEqual("inconsistent", transcript["phase_status"])
        self.assertIn("expected_1", transcript["phase_warnings"])
        self.assertEqual("invalid_cds_phase", candidate["annotation_quality"])
        self.assertEqual("", protein)


if __name__ == "__main__":
    unittest.main()
