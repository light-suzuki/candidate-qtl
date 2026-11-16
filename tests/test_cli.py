from __future__ import annotations

import csv
import json
import subprocess
import tempfile
import unittest
from pathlib import Path

from candidate_qtl.cli import main


class CandidateQtlCliTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        (self.root / "tiny.fa").write_text(">chr1\nAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT\n", encoding="utf-8")
        (self.root / "tiny.gff3").write_text(
            "##gff-version 3\n"
            "chr1\tfixture\tgene\t5\t16\t.\t+\t.\tID=GeneA;Name=GeneA\n"
            "chr1\tfixture\tmRNA\t5\t16\t.\t+\t.\tID=TxA;Parent=GeneA\n"
            "chr1\tfixture\tCDS\t5\t8\t.\t+\t0\tParent=TxA\n"
            "chr1\tfixture\tCDS\t13\t16\t.\t+\t0\tParent=TxA\n",
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_extract_accepts_gff_to_fasta_seqid_map(self) -> None:
        fasta = self.root / "alias.fa"
        gff = self.root / "alias.gff3"
        fasta.write_text(">NC_test\nATGAAATAA\n", encoding="utf-8")
        gff.write_text(
            "chrAlias\tfixture\tgene\t1\t9\t.\t+\t.\tID=AliasGene;Note=starch%20enzyme\n"
            "chrAlias\tfixture\tmRNA\t1\t9\t.\t+\t.\tID=AliasTx;Parent=AliasGene\n"
            "chrAlias\tfixture\tCDS\t1\t9\t.\t+\t0\tParent=AliasTx\n",
            encoding="utf-8",
        )
        output = self.root / "alias_extract"
        self.assertEqual(
            0,
            main(
                [
                    "extract", "--fasta", str(fasta), "--gff3", str(gff),
                    "--interval", "chrAlias:1-9", "--seqid-map", "chrAlias=NC_test",
                    "--out-dir", str(output),
                ]
            ),
        )
        row = next(csv.DictReader((output / "candidates.tsv").open(encoding="utf-8"), delimiter="\t"))
        self.assertEqual(row["description"], "starch enzyme")

    def test_extract_rank_and_cosegregation(self) -> None:
        extracted = self.root / "extract"
        self.assertEqual(
            0,
            main(
                [
                    "extract",
                    "--fasta",
                    str(self.root / "tiny.fa"),
                    "--gff3",
                    str(self.root / "tiny.gff3"),
                    "--interval",
                    "chr1:1-20",
                    "--out-dir",
                    str(extracted),
                ]
            ),
        )
        with (extracted / "candidates.tsv").open(encoding="utf-8") as handle:
            candidates = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual("GeneA", candidates[0]["gene_id"])
        ledger = self.root / "ledger.tsv"
        ledger.write_text(
            "gene_id\tortholog_ids\tfunction_evidence_source\tparental_variant\tvariant_class\tf2_cosegregation_status\n"
            "GeneA\tAT1G00010\tPMID:1\tc.1delA\tframeshift\tstrong\n",
            encoding="utf-8",
        )
        ranked = self.root / "ranked.tsv"
        self.assertEqual(0, main(["rank", "--ledger", str(ledger), "--out", str(ranked)]))
        self.assertIn("deterministic_priority_score", ranked.read_text(encoding="utf-8"))
        proteins = self.root / "proteins.fa"
        proteins.write_text(">GeneA\nMKKKK\n", encoding="utf-8")
        annotation = self.root / "annotation"
        self.assertEqual(
            0,
            main(
                [
                    "annotate-plan",
                    "--candidates",
                    str(ranked),
                    "--out-dir",
                    str(annotation),
                    "--protein-fasta",
                    str(proteins),
                ]
            ),
        )
        self.assertIn("GeneA", (annotation / "candidate_proteins.fa").read_text(encoding="utf-8"))
        subprocess.run(["bash", "-n", str(annotation / "run_local_annotation.sh")], check=True)
        packet = self.root / "packet"
        self.assertEqual(
            0,
            main(
                [
                    "research-packet",
                    "--candidates",
                    str(ranked),
                    "--out-dir",
                    str(packet),
                    "--species",
                    "Example species",
                    "--trait",
                    "example trait",
                ]
            ),
        )
        self.assertIn("GeneA", (packet / "ai_research_prompts.md").read_text(encoding="utf-8"))
        research_summary = self.root / "research_summary.tsv"
        research_summary.write_text(
            "gene_id\teurope_pmc_result_count\nGeneA\t3\n",
            encoding="utf-8",
        )
        go_evidence = self.root / "go_evidence.tsv"
        go_evidence.write_text(
            "gene_id\tgo_id\tsource\nGeneA\tGO:0000001\tUniProt\n",
            encoding="utf-8",
        )
        reranked = self.root / "reranked.tsv"
        self.assertEqual(
            0,
            main(
                [
                    "rank",
                    "--ledger",
                    str(ledger),
                    "--out",
                    str(reranked),
                    "--research-summary",
                    str(research_summary),
                    "--go-evidence",
                    str(go_evidence),
                ]
            ),
        )
        self.assertIn("research_article_count", reranked.read_text(encoding="utf-8"))
        assay_plan = self.root / "assay-plan.tsv"
        self.assertEqual(
            0,
            main(
                [
                    "assay-plan",
                    "--candidates",
                    str(ranked),
                    "--out",
                    str(assay_plan),
                    "--genome",
                    str(self.root / "tiny.fa"),
                    "--gff3",
                    str(self.root / "tiny.gff3"),
                    "--db",
                    str(self.root / "tiny-db"),
                    "--primerblast-root",
                    str(self.root),
                ]
            ),
        )
        self.assertIn("primerblast_oss_hint", assay_plan.read_text(encoding="utf-8"))
        self.assertIn("primerblast_oss assay", (self.root / "run-primerblast-oss-assays.sh").read_text(encoding="utf-8"))
        subprocess.run(["bash", "-n", str(self.root / "run-primerblast-oss-assays.sh")], check=True)
        f2 = self.root / "f2.tsv"
        f2.write_text(
            "candidate_marker_genotype\tphenotype_raw\nAA\t1.0\nAA\t1.2\nBB\t4.0\nBB\t4.2\n",
            encoding="utf-8",
        )
        summary = self.root / "summary.tsv"
        self.assertEqual(0, main(["cosegregation", "--input", str(f2), "--out", str(summary)]))
        self.assertIn("AA", summary.read_text(encoding="utf-8"))

    def test_dry_priority_separates_rescue_and_rejects_wet_evidence(self) -> None:
        candidates = self.root / "dry.tsv"
        candidates.write_text(
            "gene_id\tortholog_confidence\tmodel_phenotype_evidence\texpression_relevance\t"
            "known_qtg_or_qtl_evidence\tparalog_redundancy_risk\trescue_reason\tkey_counterevidence\n"
            "GeneStrong\thigh: RBH\tknockout phenotype\ttarget tissue\trelated QTG\tlow\t\t\n"
            "GeneWeak\tambiguous\t\t\t\tmedium\t\t\n"
            "GeneUnknown\t\t\t\t\t\tnear QTL peak\t\n",
            encoding="utf-8",
        )
        output = self.root / "dry_ranked.tsv"
        self.assertEqual(
            0,
            main(["dry-prioritize", "--candidates", str(candidates), "--out", str(output), "--top-n", "1", "--reserve-n", "1"]),
        )
        with output.open(encoding="utf-8") as handle:
            rows = {row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual("top", rows["GeneStrong"]["dry_group"])
        self.assertEqual("reserve", rows["GeneWeak"]["dry_group"])
        self.assertEqual("rescue", rows["GeneUnknown"]["dry_group"])
        self.assertIn("model_phenotype_evidence:+3", rows["GeneStrong"]["score_breakdown"])

        # Every candidate carries an explicit selection/non-selection reason and
        # Top/Reserve genes gain recommended parental-sequencing regions.
        self.assertTrue(rows["GeneStrong"]["selection_reason"].startswith("top:"))
        self.assertTrue(rows["GeneUnknown"]["selection_reason"].startswith("rescue:"))
        self.assertIn("CDS", rows["GeneStrong"]["recommended_regions"])
        self.assertIn("exon-intron", rows["GeneWeak"]["recommended_regions"])
        self.assertEqual("", rows["GeneUnknown"]["recommended_regions"])

        leaked = self.root / "dry_with_wet.tsv"
        leaked.write_text("gene_id\tparental_variant\nGeneA\tc.1A>G\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Wet-only evidence"):
            main(["dry-prioritize", "--candidates", str(leaked), "--out", str(self.root / "bad.tsv")])

    def test_dry_priority_completes_without_wet_columns_and_labels_all(self) -> None:
        candidates = self.root / "dry_all.tsv"
        candidates.write_text(
            "gene_id\tortholog_confidence\tmodel_phenotype_evidence\n"
            "G1\thigh: RBH\tknockout phenotype\n"
            "G2\t\t\n"
            "G3\t\t\n",
            encoding="utf-8",
        )
        output = self.root / "dry_all_ranked.tsv"
        self.assertEqual(0, main(["dry-prioritize", "--candidates", str(candidates), "--out", str(output), "--top-n", "1", "--reserve-n", "1"]))
        with output.open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        # No candidate is left without a reason, even with no Wet evidence present.
        self.assertTrue(all(row["selection_reason"] for row in rows))
        groups = {row["gene_id"]: row["dry_group"] for row in rows}
        self.assertEqual(groups["G1"], "top")
        self.assertIn("not_selected_yet", set(groups.values()))

    def test_cosegregation_categorical(self) -> None:
        f2 = self.root / "cat.tsv"
        f2.write_text(
            "candidate_marker_genotype\tphenotype_raw\n"
            "AA\tnormal\nAA\tnormal\nAB\tnormal\nBB\tfasciated\nBB\tfasciated\n",
            encoding="utf-8",
        )
        out = self.root / "cat_summary.tsv"
        self.assertEqual(0, main(["cosegregation", "--input", str(f2), "--out", str(out)]))
        text = out.read_text(encoding="utf-8")
        self.assertIn("phenotype_category", text)
        self.assertIn("fasciated", text)
        result = json.loads(out.with_suffix(".json").read_text(encoding="utf-8"))
        self.assertEqual("categorical", result["mode"])

    def test_rank_excludes_cosegregation(self) -> None:
        header = "gene_id\tortholog_ids\tfunction_evidence_source\tparental_variant\tvariant_class\tf2_cosegregation_status\n"
        strong = self.root / "led_strong.tsv"
        strong.write_text(header + "GeneA\tAT1G00010\tPMID:1\tc.1delA\tframeshift\tstrong\n", encoding="utf-8")
        excluded = self.root / "led_excluded.tsv"
        excluded.write_text(header + "GeneA\tAT1G00010\tPMID:1\tc.1delA\tframeshift\texcluded\n", encoding="utf-8")
        out_strong = self.root / "rank_strong.tsv"
        out_excluded = self.root / "rank_excluded.tsv"
        self.assertEqual(0, main(["rank", "--ledger", str(strong), "--out", str(out_strong)]))
        self.assertEqual(0, main(["rank", "--ledger", str(excluded), "--out", str(out_excluded)]))
        with out_strong.open(encoding="utf-8") as handle:
            row_strong = next(csv.DictReader(handle, delimiter="\t"))
        with out_excluded.open(encoding="utf-8") as handle:
            row_excluded = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(row_strong["deterministic_priority_score"], row_excluded["deterministic_priority_score"])


if __name__ == "__main__":
    unittest.main()
