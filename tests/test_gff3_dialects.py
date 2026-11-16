from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path


from candidate_qtl.cli import main


class Gff3DialectTests(unittest.TestCase):
    def _extract(self, sequence: str, gff_text: str):
        self.tempdir = tempfile.TemporaryDirectory()
        root = Path(self.tempdir.name)
        fasta = root / "genome.fa"
        fasta.write_text(f">chr1\n{sequence}\n", encoding="utf-8")
        gff = root / "genes.gff3"
        gff.write_text(gff_text, encoding="utf-8")
        out = root / "out"
        self.assertEqual(
            0,
            main(["extract", "--fasta", str(fasta), "--gff3", str(gff), "--interval", f"chr1:1-{len(sequence)}", "--out-dir", str(out)]),
        )
        with (out / "candidates.tsv").open(encoding="utf-8") as handle:
            candidates = list(csv.DictReader(handle, delimiter="\t"))
        manifest = json.loads((out / "extract_manifest.json").read_text(encoding="utf-8"))
        return candidates, manifest, out

    def test_transcript_without_gene_gets_placeholder(self) -> None:
        candidates, manifest, _ = self._extract(
            "ATGAAAAAAATAG",
            "##gff-version 3\n"
            "chr1\tx\tmRNA\t1\t13\t.\t+\t.\tID=TxOnly\n"
            "chr1\tx\tCDS\t1\t5\t.\t+\t0\tParent=TxOnly\n"
            "chr1\tx\tCDS\t10\t13\t.\t+\t1\tParent=TxOnly\n",
        )
        self.assertEqual(len(candidates), 1)
        row = candidates[0]
        self.assertEqual(row["gene_id"], "gene:TxOnly")
        self.assertEqual(row["is_placeholder_gene"], "True")
        self.assertEqual(row["annotation_quality"], "complete_orf")
        self.assertEqual(manifest["placeholder_gene_count"], 1)
        self.assertEqual(manifest["gff_dialect_qc"]["transcripts_without_gene_parent"], 1)

    def test_cds_directly_on_gene(self) -> None:
        candidates, manifest, _ = self._extract(
            "ATGAAAAAAATAG",
            "##gff-version 3\n"
            "chr1\tx\tgene\t1\t13\t.\t+\t.\tID=GeneG\n"
            "chr1\tx\tCDS\t1\t5\t.\t+\t0\tParent=GeneG\n"
            "chr1\tx\tCDS\t10\t13\t.\t+\t1\tParent=GeneG\n",
        )
        row = candidates[0]
        self.assertEqual(row["gene_id"], "GeneG")
        self.assertEqual(row["is_placeholder_gene"], "False")
        self.assertEqual(row["annotation_quality"], "complete_orf")

    def test_cds_only_with_missing_transcript_and_gene(self) -> None:
        candidates, manifest, out = self._extract(
            "ATGAAAAAAATAG",
            "##gff-version 3\n"
            "chr1\tx\tCDS\t1\t5\t.\t+\t0\tParent=orphanTx\n"
            "chr1\tx\tCDS\t10\t13\t.\t+\t1\tParent=orphanTx\n",
        )
        row = candidates[0]
        self.assertEqual(row["gene_id"], "gene:orphanTx")
        self.assertEqual(row["is_placeholder_gene"], "True")
        self.assertEqual(manifest["gff_dialect_qc"]["orphan_cds_features"], 2)
        # An ID mapping table is emitted linking gene/transcript/protein.
        with (out / "id_map.tsv").open(encoding="utf-8") as handle:
            id_rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(id_rows[0]["gene_id"], "gene:orphanTx")
        self.assertEqual(id_rows[0]["transcript_id"], "orphanTx")

    def test_pseudogene_is_typed_and_not_forced_to_protein(self) -> None:
        candidates, manifest, _ = self._extract(
            "ATGAAAAAAATAG",
            "##gff-version 3\n"
            "chr1\tx\tpseudogene\t1\t13\t.\t+\t.\tID=PseudoA\n",
        )
        row = candidates[0]
        self.assertEqual(row["gene_id"], "PseudoA")
        self.assertEqual(row["gene_biotype"], "pseudogene")
        self.assertEqual(row["annotation_quality"], "noncoding_or_pseudogene")

    def test_duplicate_parent_ids_are_reported(self) -> None:
        candidates, manifest, out = self._extract(
            "ATGAAAAAAATAGATGAAAAAAATAG",
            "##gff-version 3\n"
            "chr1\tx\tgene\t1\t13\t.\t+\t.\tID=DupGene\n"
            "chr1\tx\tmRNA\t1\t13\t.\t+\t.\tID=DupTx;Parent=DupGene\n"
            "chr1\tx\tCDS\t1\t5\t.\t+\t0\tParent=DupTx\n"
            "chr1\tx\tCDS\t10\t13\t.\t+\t1\tParent=DupTx\n"
            "chr1\tx\tgene\t14\t26\t.\t+\t.\tID=DupGene\n"
            "chr1\tx\tmRNA\t14\t26\t.\t+\t.\tID=DupTx2;Parent=DupGene\n"
            "chr1\tx\tCDS\t14\t18\t.\t+\t0\tParent=DupTx2\n"
            "chr1\tx\tCDS\t23\t26\t.\t+\t1\tParent=DupTx2\n",
        )
        self.assertIn("DupGene", manifest["gff_dialect_qc"]["duplicate_parent_ids"])
        self.assertEqual(1, len(candidates))
        self.assertEqual("DupGene", candidates[0]["gene_id"])
        self.assertEqual(1, (out / "candidate_genomic.fa").read_text(encoding="utf-8").count(">DupGene\n"))


if __name__ == "__main__":
    unittest.main()
