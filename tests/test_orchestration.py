from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path

from candidate_qtl.cli import main


class OrchestrationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.candidates = self.root / "dry_candidates.tsv"
        self.candidates.write_text(
            "gene_id\tdry_group\tdry_rank\n"
            "GeneA\ttop\t1\n"
            "GeneB\treserve\t2\n",
            encoding="utf-8",
        )
        self.genomic = self.root / "candidate_genomic.fa"
        self.genomic.write_text(
            ">GeneA\nAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT\n"
            ">GeneB\nTTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA\n",
            encoding="utf-8",
        )
        self.variants = self.root / "parental_variants.tsv"
        self.variants.write_text(
            "variant_id\tgene_id\tchrom\tpos\tref\talt\tvalidation_status\n"
            "GeneA_snp1\tGeneA\tchr1\t10\tA\tG\tconfirmed\n"
            "GeneA_indel1\tGeneA\tchr1\t20\tAT\tA\tconfirmed\n"
            "GeneA_pending\tGeneA\tchr1\t30\tC\tT\tpending\n",
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_primer_plans_and_orchestration(self) -> None:
        direct = self.root / "direct"
        self.assertEqual(
            0,
            main(
                [
                    "primer-plan",
                    "--stage",
                    "direct-sequencing",
                    "--out-dir",
                    str(direct),
                    "--primerblast-root",
                    str(self.root / "primerblast-oss"),
                    "--genome",
                    str(self.root / "genome.fa"),
                    "--db",
                    str(self.root / "genome-db"),
                    "--candidates",
                    str(self.candidates),
                    "--top-n",
                    "1",
                    "--gff3",
                    str(self.root / "genes.gff3"),
                    "--candidate-genomic-fasta",
                    str(self.genomic),
                ]
            ),
        )
        direct_script = direct / "run_direct_sequencing.sh"
        text = direct_script.read_text(encoding="utf-8")
        self.assertIn("primerblast_oss assay", text)
        self.assertIn("primerblast_oss tile", text)
        self.assertIn("GeneA", text)
        self.assertNotIn("GeneB.assay", text)
        subprocess.run(["bash", "-n", str(direct_script)], check=True)

        marker = self.root / "marker"
        self.assertEqual(
            0,
            main(
                [
                    "primer-plan",
                    "--stage",
                    "variant-marker",
                    "--out-dir",
                    str(marker),
                    "--primerblast-root",
                    str(self.root / "primerblast-oss"),
                    "--genome",
                    str(self.root / "genome.fa"),
                    "--db",
                    str(self.root / "genome-db"),
                    "--candidates",
                    str(self.candidates),
                    "--top-n",
                    "1",
                    "--variants",
                    str(self.variants),
                ]
            ),
        )
        marker_script = marker / "run_variant_marker.sh"
        marker_text = marker_script.read_text(encoding="utf-8")
        self.assertIn("--snp chr1:10 --alt G", marker_text)
        self.assertIn("--interval chr1:1-322", marker_text)
        self.assertIn("--format order", marker_text)
        marker_plan = (marker / "variant_marker_plan.tsv").read_text(encoding="utf-8")
        self.assertNotIn("GeneA_pending", marker_plan)
        marker_manifest = (marker / "variant_marker_manifest.json").read_text(encoding="utf-8")
        self.assertIn('"excluded_unconfirmed_count": 1', marker_manifest)
        subprocess.run(["bash", "-n", str(marker_script)], check=True)

        grid = self.root / "grid"
        self.assertEqual(
            0,
            main(
                [
                    "primer-plan",
                    "--stage",
                    "qtl-grid",
                    "--out-dir",
                    str(grid),
                    "--primerblast-root",
                    str(self.root / "primerblast-oss"),
                    "--genome",
                    str(self.root / "genome.fa"),
                    "--db",
                    str(self.root / "genome-db"),
                    "--interval",
                    "chr1:1-1000",
                    "--n-markers",
                    "6",
                ]
            ),
        )
        grid_script = grid / "run_qtl_grid.sh"
        self.assertIn("primerblast_oss markers", grid_script.read_text(encoding="utf-8"))
        subprocess.run(["bash", "-n", str(grid_script)], check=True)

        orchestration = self.root / "orchestration"
        self.assertEqual(
            0,
            main(
                [
                    "orchestrate",
                    "--candidates",
                    str(self.candidates),
                    "--out-dir",
                    str(orchestration),
                    "--primerblast-root",
                    str(self.root / "primerblast-oss"),
                    "--genome",
                    str(self.root / "genome.fa"),
                    "--gff3",
                    str(self.root / "genes.gff3"),
                    "--nucleotide-db",
                    str(self.root / "genome-db"),
                    "--candidate-genomic-fasta",
                    str(self.genomic),
                    "--interval",
                    "chr1:1-1000",
                    "--variants",
                    str(self.variants),
                    "--top-n",
                    "1",
                ]
            ),
        )
        run_all = orchestration / "run_available_steps.sh"
        subprocess.run(["bash", "-n", str(run_all)], check=True)
        manifest = (orchestration / "orchestration_manifest.json").read_text(encoding="utf-8")
        self.assertIn("primerblast-oss assay/tile", manifest)
        self.assertIn("selection-marker design", manifest)
        self.assertIn('"memory_target_gb": 16', manifest)


if __name__ == "__main__":
    unittest.main()
