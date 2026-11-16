from __future__ import annotations

import csv
import subprocess
import tempfile
import unittest
from pathlib import Path

from candidate_qtl.exon_assays import exon_sequencing_plan


class ExonAssayPlanTests(unittest.TestCase):
    def test_every_cds_base_is_covered_by_target_windows(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            candidates = root / "candidates.tsv"
            candidates.write_text(
                "gene_id\tchrom\tstrand\trepresentative_transcript\tdry_group\tdry_rank\n"
                "GeneA\tchr1\t+\tTxA\ttop\t1\n",
                encoding="utf-8",
            )
            gff = root / "genes.gff3"
            gff.write_text(
                "##gff-version 3\n"
                "chr1\tx\tgene\t100\t1200\t.\t+\t.\tID=GeneA\n"
                "chr1\tx\tmRNA\t100\t1200\t.\t+\t.\tID=TxA;Parent=GeneA\n"
                "chr1\tx\tCDS\t100\t199\t.\t+\t0\tParent=TxA\n"
                "chr1\tx\tCDS\t300\t1099\t.\t+\t0\tParent=TxA\n",
                encoding="utf-8",
            )
            genome = root / "genome.fa"
            genome.write_text(">chr1\n" + "ACGT" * 400 + "\n", encoding="utf-8")
            out = root / "plan"
            result = exon_sequencing_plan(
                str(candidates),
                str(gff),
                str(out),
                str(root / "primerblast-oss"),
                str(genome),
                [str(root / "genome_db")],
                top_n=1,
                max_target_span=300,
                overlap=50,
            )
            self.assertGreater(result["planned_template_count"], 2)
            with Path(result["plan"]).open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            covered = set()
            for row in rows:
                covered.update(range(int(row["target_start"]), int(row["target_end"]) + 1))
                self.assertTrue(Path(row["template_fasta"]).exists())
            expected = set(range(100, 200)) | set(range(300, 1100))
            self.assertTrue(expected.issubset(covered))
            script_text = Path(result["script"]).read_text(encoding="utf-8")
            self.assertIn(" primerblast_oss design ", script_text)
            self.assertIn(" --target ", script_text)
            subprocess.run(["bash", "-n", result["script"]], check=True)
            self.assertEqual([], result["genes_without_cds"])

    def test_gene_without_representative_cds_is_reported(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            candidates = root / "candidates.tsv"
            candidates.write_text(
                "gene_id\tchrom\tstrand\trepresentative_transcript\nGeneMissing\tchr1\t+\tTxMissing\n",
                encoding="utf-8",
            )
            gff = root / "genes.gff3"
            gff.write_text("##gff-version 3\n", encoding="utf-8")
            genome = root / "genome.fa"
            genome.write_text(">chr1\n" + "A" * 1000 + "\n", encoding="utf-8")
            result = exon_sequencing_plan(
                str(candidates), str(gff), str(root / "out"), str(root / "pb"), str(genome), [str(root / "db")]
            )
            self.assertEqual(["GeneMissing"], result["genes_without_cds"])
            self.assertEqual(0, result["planned_template_count"])


if __name__ == "__main__":
    unittest.main()
