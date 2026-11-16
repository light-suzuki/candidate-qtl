from __future__ import annotations

import csv
import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from Bio import SeqIO

from candidate_qtl.cli import main


class EvidencePipelineTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_annotation_keyword_evidence_links_pipe_ids_to_gff(self) -> None:
        candidates = self.root / "keyword_candidates.tsv"
        candidates.write_text("gene_id\trepresentative_transcript\nGeneA\tTxA\nGeneB\tTxB\n", encoding="utf-8")
        blast = self.root / "keyword_blast.tsv"
        blast.write_text(
            "GeneA\tZW6|rna|Psat03G0000100-T1\t99\t100\t100\t100\t1\t100\t1\t100\t1e-50\t200\ttitle\n"
            "GeneB\tZW6|rna|Psat03G0000200-T1\t99\t100\t100\t100\t1\t100\t1\t100\t1e-50\t200\ttitle\n",
            encoding="utf-8",
        )
        gff = self.root / "keyword.gff3"
        gff.write_text(
            "chr3\tx\tgene\t1\t100\t.\t+\t.\tID=pissa.ZW6.gnm1.ann1.Psat03G0000100;Note=starch%20branching%20enzyme\n"
            "chr3\tx\tgene\t200\t300\t.\t+\t.\tID=pissa.ZW6.gnm1.ann1.Psat03G0000200;Note=protein%20kinase\n",
            encoding="utf-8",
        )
        output = self.root / "keyword_evidence.tsv"
        self.assertEqual(
            0,
            main(
                [
                    "annotation-keywords", "--candidates", str(candidates), "--blast", str(blast),
                    "--annotation-gff3", str(gff), "--keyword", "starch", "--out", str(output),
                ]
            ),
        )
        text = output.read_text(encoding="utf-8")
        self.assertIn("GeneA\ttrait_function_match", text)
        self.assertNotIn("GeneB\ttrait_function_match", text)

    def test_eggnog_v2_header_aliases_are_imported(self) -> None:
        candidates = self.root / "v2_candidates.tsv"
        candidates.write_text("gene_id\trepresentative_transcript\nGeneA\tTxA\n", encoding="utf-8")
        eggnog = self.root / "v2_eggnog.tsv"
        eggnog.write_text(
            "#query_name\tseed_ortholog_evalue\tPreferred_name\tGOs\teggNOG OGs\teggNOG free text desc.\n"
            "TxA\t1e-30\tSBE1\tGO:0003844\tOG0001@Viridiplantae\tstarch branching enzyme\n",
            encoding="utf-8",
        )
        out = self.root / "v2_evidence"
        self.assertEqual(
            0,
            main(["evidence-build", "--candidates", str(candidates), "--eggnog", str(eggnog), "--out-dir", str(out)]),
        )
        text = (out / "evidence.tsv").read_text(encoding="utf-8")
        self.assertIn("orthogroup_supported", text)
        self.assertIn("starch branching enzyme", text)

    def test_transcript_aware_extract_does_not_merge_isoforms(self) -> None:
        fasta = self.root / "isoforms.fa"
        fasta.write_text(">chr1\nATGAAATAAATGCCCTAA\n", encoding="utf-8")
        gff = self.root / "isoforms.gff3"
        gff.write_text(
            "##gff-version 3\n"
            "chr1\tfixture\tgene\t1\t18\t.\t+\t.\tID=GeneIso\n"
            "chr1\tfixture\tmRNA\t1\t9\t.\t+\t.\tID=TxCanonical;Parent=GeneIso;tag=canonical\n"
            "chr1\tfixture\tCDS\t1\t9\t.\t+\t0\tParent=TxCanonical\n"
            "chr1\tfixture\tmRNA\t1\t18\t.\t+\t.\tID=TxLong;Parent=GeneIso\n"
            "chr1\tfixture\tCDS\t1\t9\t.\t+\t0\tParent=TxLong\n"
            "chr1\tfixture\tCDS\t10\t18\t.\t+\t0\tParent=TxLong\n",
            encoding="utf-8",
        )
        out = self.root / "extract"
        self.assertEqual(0, main(["extract", "--fasta", str(fasta), "--gff3", str(gff), "--interval", "chr1:1-18", "--out-dir", str(out)]))
        with (out / "candidates.tsv").open(encoding="utf-8") as handle:
            candidate = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual("TxCanonical", candidate["representative_transcript"])
        self.assertEqual("2", candidate["transcript_count"])
        self.assertEqual("9", candidate["cds_length"])
        with (out / "candidate_transcript_cds.fa").open(encoding="utf-8") as handle:
            transcript_records = list(SeqIO.parse(handle, "fasta"))
        self.assertEqual({"TxCanonical", "TxLong"}, {record.id for record in transcript_records})
        representative = SeqIO.read(out / "candidate_cds.fa", "fasta")
        self.assertEqual(9, len(representative.seq))

    def test_structured_evidence_ranking_and_top_only_orthology(self) -> None:
        candidates = self.root / "candidates.tsv"
        candidates.write_text(
            "gene_id\tchrom\tstart\tend\tcds_length\tannotation_quality\trepresentative_transcript\tinterval\n"
            "GeneA\tchr1\t100\t200\t300\tcomplete_orf\tTxA\tchr1:1-1000\n"
            "GeneB\tchr1\t800\t900\t0\tno_cds\t\tchr1:1-1000\n",
            encoding="utf-8",
        )
        blast = self.root / "blast.tsv"
        blast.write_text(
            "GeneA\tRelA1\t60\t80\t100\t120\t1\t80\t1\t80\t1e-30\t100\tknown regulator\n"
            "GeneA\tRelA2\t58\t78\t100\t120\t1\t78\t1\t78\t1e-25\t95\tclose paralog\n"
            "GeneB\tRelB1\t45\t70\t100\t110\t1\t70\t1\t70\t1e-10\t80\tuncharacterized conserved protein\n",
            encoding="utf-8",
        )
        hmmer = self.root / "hmmer.tbl"
        hmmer.write_text("PF0001 - GeneA - 1e-20 50.0 0.0 1e-20 50.0 0.0 1 1 1 1 1 1 1 1 regulator_domain\n", encoding="utf-8")
        eggnog = self.root / "eggnog.annotations"
        eggnog.write_text(
            "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\n"
            "GeneA\tSeed1\t1e-20\t100\tOG0001@Viridiplantae\tViridiplantae\tT\tregulatory protein\tREG1\tGO:0006355\n",
            encoding="utf-8",
        )
        expression = self.root / "expression.tsv"
        expression.write_text(
            "gene_id\ttissue\tvalue\tis_target\n"
            "GeneA\tmeristem\t10\ttrue\nGeneA\tleaf\t1\tfalse\n"
            "GeneB\tmeristem\t0.1\ttrue\nGeneB\tleaf\t0.2\tfalse\n",
            encoding="utf-8",
        )
        manual = self.root / "manual.tsv"
        manual.write_text(
            "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\n"
            "GeneA\tdirect_model_genetic\tsupport\thigh\tmutant phenotype\tcuration\tPMID:1\tdirect genetic evidence\n",
            encoding="utf-8",
        )
        evidence_dir = self.root / "evidence"
        self.assertEqual(
            0,
            main(
                [
                    "evidence-build", "--candidates", str(candidates), "--out-dir", str(evidence_dir),
                    "--blast", str(blast), "--hmmer", str(hmmer), "--eggnog", str(eggnog), "--expression", str(expression),
                    "--manual-evidence", str(manual), "--peak", "chr1:150", "--interval", "chr1:1-1000",
                ]
            ),
        )
        evidence_text = (evidence_dir / "evidence.tsv").read_text(encoding="utf-8")
        self.assertIn("expression_target_specific", evidence_text)
        self.assertIn("orthogroup_supported", evidence_text)
        self.assertIn("annotation_incomplete", evidence_text)
        self.assertIn("paralog_ambiguous", evidence_text)

        rank_dir = self.root / "rank"
        self.assertEqual(
            0,
            main(
                [
                    "evidence-rank", "--candidates", str(candidates), "--evidence", str(evidence_dir / "evidence.tsv"),
                    "--out-dir", str(rank_dir), "--top-n", "1", "--reserve-n", "0",
                ]
            ),
        )
        with (rank_dir / "ranked_candidates.tsv").open(encoding="utf-8") as handle:
            ranked = {row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual("top", ranked["GeneA"]["dry_group"])
        self.assertEqual("rescue", ranked["GeneB"]["dry_group"])
        contributions = (rank_dir / "score_contributions.tsv").read_text(encoding="utf-8")
        # Generic domain/homology/GO are imported as context and must not be scored.
        self.assertNotIn("domain_supported", contributions)
        self.assertNotIn("domain_context", contributions)
        self.assertNotIn("homology_context", contributions)
        evidence_text = (evidence_dir / "evidence.tsv").read_text(encoding="utf-8")
        self.assertIn("domain_context", evidence_text)
        self.assertIn("go_annotation_context", evidence_text)
        decision = json.loads((rank_dir / "candidate_path_decision.json").read_text(encoding="utf-8"))
        self.assertEqual("test_coding_candidates_first", decision["mode"])
        self.assertIn("interval_contains_gene_models_without_cds", decision["alerts"])

        proteins = self.root / "candidate_proteins.fa"
        sequence_a = "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPQ"
        sequence_b = "MNNIRRVAILGCGNIGTDHLAARLAAAGADVVVVATGS"
        proteins.write_text(f">GeneA\n{sequence_a}\n>GeneB\n{sequence_b}\n", encoding="utf-8")
        target = self.root / "target.fa"
        target.write_text(f">GeneA\n{sequence_a}\n>GeneB\n{sequence_b}\n", encoding="utf-8")
        relative = self.root / "relative.fa"
        relative.write_text(f">RelA1\n{sequence_a}\n", encoding="utf-8")
        orthology_dir = self.root / "orthology"
        self.assertEqual(
            0,
            main(
                [
                    "orthology-plan", "--ranked", str(rank_dir / "ranked_candidates.tsv"),
                    "--protein-fasta", str(proteins), "--target-proteome", str(target),
                    "--relative-proteome", f"relative={relative}", "--out-dir", str(orthology_dir), "--top-n", "1",
                ]
            ),
        )
        subprocess.run(["bash", "-n", str(orthology_dir / "run_orthology_refinement.sh")], check=True)
        manifest = json.loads((orthology_dir / "orthology_plan.json").read_text(encoding="utf-8"))
        self.assertEqual(["GeneA"], manifest["selected_genes"])
        result_set = manifest["result_sets"][0]
        if all(shutil.which(command) for command in ("blastp", "makeblastdb", "blastdbcmd")):
            subprocess.run(["bash", str(orthology_dir / "run_orthology_refinement.sh")], check=True)
        else:
            Path(result_set["forward"]).parent.mkdir(parents=True, exist_ok=True)
            Path(result_set["forward"]).write_text("GeneA\tRelA1\t90\t40\t40\t40\t1\t40\t1\t40\t1e-20\t100\thit\n", encoding="utf-8")
            Path(result_set["reverse"]).write_text("RelA1\tGeneA\t90\t40\t40\t40\t1\t40\t1\t40\t1e-20\t100\thit\n", encoding="utf-8")
        refined = self.root / "refined.tsv"
        self.assertEqual(0, main(["orthology-import", "--manifest", str(orthology_dir / "orthology_plan.json"), "--out", str(refined)]))
        self.assertIn("rbh_ortholog", refined.read_text(encoding="utf-8"))

        pipeline_dir = self.root / "pipeline"
        self.assertEqual(
            0,
            main(
                [
                    "candidate-pipeline", "--candidates", str(candidates), "--out-dir", str(pipeline_dir),
                    "--blast", str(blast), "--hmmer", str(hmmer), "--eggnog", str(eggnog), "--expression", str(expression),
                    "--manual-evidence", str(manual), "--refined-evidence", str(refined),
                    "--peak", "chr1:150", "--interval", "chr1:1-1000", "--top-n", "1", "--reserve-n", "0",
                    "--protein-fasta", str(proteins), "--target-proteome", str(target),
                    "--relative-proteome", f"relative={relative}",
                ]
            ),
        )
        pipeline_manifest = json.loads((pipeline_dir / "candidate_pipeline.json").read_text(encoding="utf-8"))
        self.assertTrue(str(pipeline_manifest["final_ranked_candidates"]).endswith("05_final_rank/ranked_candidates.tsv"))
        self.assertIn("rbh_ortholog", (pipeline_dir / "04_refined_evidence" / "evidence.tsv").read_text(encoding="utf-8"))

    def test_expression_only_and_zero_evidence_do_not_become_top(self) -> None:
        candidates = self.root / "minimal_candidates.tsv"
        candidates.write_text(
            "gene_id\tannotation_quality\tcds_length\n"
            "GeneExpression\tcomplete_orf\t300\nGeneZero\tcomplete_orf\t300\n",
            encoding="utf-8",
        )
        evidence = self.root / "minimal_evidence.tsv"
        evidence.write_text(
            "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\n"
            "GeneExpression\texpression_target_specific\tsupport\thigh\t10\texpression\texpr.tsv\ttarget enriched\n"
            "GeneExpression\texpression_detected\tsupport\tmedium\t10\texpression\texpr.tsv\tdetected\n",
            encoding="utf-8",
        )
        out = self.root / "minimal_rank"
        self.assertEqual(0, main(["evidence-rank", "--candidates", str(candidates), "--evidence", str(evidence), "--out-dir", str(out)]))
        with (out / "ranked_candidates.tsv").open(encoding="utf-8") as handle:
            ranked = {row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertNotEqual("top", ranked["GeneExpression"]["dry_group"])
        self.assertEqual("evidence_insufficient", ranked["GeneZero"]["dry_group"])
        self.assertEqual("insufficient_independent_evidence", ranked["GeneExpression"]["top_eligibility"])


if __name__ == "__main__":
    unittest.main()
