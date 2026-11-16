from __future__ import annotations

import json
import shutil
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main


class GffPreflightTests(unittest.TestCase):
    def test_provider_missing_transcript_is_not_silently_dropped_from_review(self) -> None:
        tutorial = Path(__file__).parents[1] / "docs" / "tutorial"
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); empty = root / "empty.fa"; empty.write_text("", encoding="utf-8")
            out = root / "out"
            self.assertEqual(0, main(["gff-preflight", "--fasta", str(tutorial / "genome.fa"), "--gff3", str(tutorial / "genes.gff3"),
                                      "--interval", "chr1:1-900", "--agat-cds", str(empty), "--gffread-cds", str(empty), "--out-dir", str(out)]))
            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertGreater(len(reviews), 0)
            self.assertEqual({"provider_missing_transcript"}, set(reviews["reason"]))

    def test_plan_mode_without_provider_crosscheck_does_not_flag_every_transcript(self) -> None:
        tutorial = Path(__file__).parents[1] / "docs" / "tutorial"
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "out"
            self.assertEqual(0, main(["gff-preflight", "--fasta", str(tutorial / "genome.fa"), "--gff3", str(tutorial / "genes.gff3"),
                                      "--interval", "chr1:1-900", "--out-dir", str(out)]))
            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertNotIn("provider_missing_transcript", set(reviews["reason"]))
            manifest = json.loads((out / "gff_preflight_manifest.json").read_text(encoding="utf-8"))
            self.assertFalse(manifest["cross_check_requested"])
            self.assertEqual(0, manifest["provider_missing_transcripts"])

    def test_crosscheck_reports_parser_disagreement_without_replacing_raw_gff(self) -> None:
        tutorial = Path(__file__).parents[1] / "docs" / "tutorial"
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); normalized = root / "normalized.gff3"; shutil.copyfile(tutorial / "genes.gff3", normalized)
            sequence = "ATG" + "A" * 54 + "TAA"
            agat = root / "agat.fa"; gffread = root / "gffread.fa"
            agat.write_text("".join(f">{gene}.1\n{sequence}\n" for gene in ("GENE_CAUSAL", "DECOY1", "DECOY2")), encoding="utf-8")
            gffread.write_text(f">GENE_CAUSAL.1\n{sequence[:-1]}C\n>DECOY1.1\n{sequence}\n>DECOY2.1\n{sequence}\n", encoding="utf-8")
            out = root / "out"
            self.assertEqual(0, main(["gff-preflight", "--fasta", str(tutorial / "genome.fa"), "--gff3", str(tutorial / "genes.gff3"),
                                      "--interval", "chr1:1-900", "--normalized-gff", str(normalized), "--agat-cds", str(agat),
                                      "--gffread-cds", str(gffread), "--out-dir", str(out)]))
            rows = pd.read_csv(out / "gff_concordance.tsv", sep="\t", dtype=str).fillna("").set_index("transcript_id")
            self.assertEqual("True", rows.loc["GENE_CAUSAL.1", "review_required"])
            self.assertEqual("False", rows.loc["DECOY1.1", "review_required"])
            manifest = json.loads((out / "gff_preflight_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(str(tutorial / "genes.gff3"), manifest["candidate_extraction_gff"])
            self.assertTrue(manifest["cross_check_only"])
            self.assertEqual(1, manifest["discordant_transcripts"])
            self.assertIn("agat_convert_sp_gxf2gxf.pl", (out / "run_gff_crosscheck.sh").read_text(encoding="utf-8"))


if __name__ == "__main__": unittest.main()
