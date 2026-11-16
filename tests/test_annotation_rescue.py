from __future__ import annotations

import json
import shutil
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pandas as pd
import yaml

from candidate_qtl.cli import main


class AnnotationRescueTests(unittest.TestCase):
    def _fixture(self, root: Path) -> tuple[Path, Path, Path, Path]:
        fasta = root / "genome.fa"
        fasta.write_text(">chr1\n" + "ACGT" * 100 + "\n", encoding="utf-8")
        candidates = root / "candidates.tsv"
        candidates.write_text(
            "gene_id\tchrom\tstart\tend\tstrand\tannotation_quality\n"
            "NativeGood\tchr1\t50\t130\t+\tcomplete_orf\n"
            "NativeBroken\tchr1\t220\t300\t+\tinvalid_cds_phase\n",
            encoding="utf-8",
        )
        queries = root / "queries.faa"
        queries.write_text(">Q_NATIVE\n" + "M" * 20 + "\n>Q_RESCUE\n" + "M" * 20 + "\n>Q_UNMAPPED\n" + "M" * 20 + "\n", encoding="utf-8")
        metadata = root / "query_metadata.tsv"
        metadata.write_text(
            "query_protein_id\tspecies\tassembly_id\tannotation_release\tsource_id\tquery_role\tselection_method\treview_status\n"
            "Q_NATIVE\tPisum sativum\tV1a\tr1\tUniProt:Q_NATIVE\tfunctional_homolog\tautomated_sequence_backed_search\treviewed\n"
            "Q_RESCUE\tPisum sativum\tZW6\tr1\tUniProt:Q_RESCUE\tfunctional_homolog\tautomated_sequence_backed_search\treviewed\n"
            "Q_UNMAPPED\tPisum sativum\tZW6\tr1\tUniProt:Q_UNMAPPED\tnegative_control\tpredeclared_negative_control\treviewed\n",
            encoding="utf-8",
        )
        return fasta, candidates, queries, metadata

    def test_complete_import_separates_native_rescue_and_unmapped_without_promotion(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta, candidates, queries, metadata = self._fixture(root)
            projection = root / "projection.gff3"
            projection.write_text(
                "##gff-version 3\n"
                "chr1__1_400\tminiprot\tmRNA\t55\t125\t.\t+\t.\tID=MP1;Rank=1;Identity=0.9;Positive=0.95;Frameshift=0;StopCodon=0;Target=Q_NATIVE 1 20\n"
                "chr1__1_400\tminiprot\tCDS\t55\t80\t.\t+\t0\tParent=MP1\n"
                "chr1__1_400\tminiprot\tCDS\t100\t125\t.\t+\t0\tParent=MP1\n"
                "chr1__1_400\tminiprot\tstop_codon\t126\t128\t.\t+\t.\tParent=MP1\n"
                "chr1__1_400\tminiprot\tmRNA\t145\t205\t.\t+\t.\tID=MP2;Rank=1;Identity=0.8;Positive=0.9;Frameshift=0;StopCodon=0;Target=Q_RESCUE 1 20\n"
                "chr1__1_400\tminiprot\tCDS\t145\t205\t.\t+\t0\tParent=MP2\n",
                encoding="utf-8",
            )
            out = root / "out"
            self.assertEqual(0, main(["annotation-rescue", "--fasta", str(fasta), "--candidates", str(candidates),
                                      "--query-proteins", str(queries), "--query-metadata", str(metadata), "--interval", "chr1:1-400",
                                      "--projection-gff", str(projection), "--projection-complete", "--out-dir", str(out)]))
            rows = pd.read_csv(out / "annotation_rescue.tsv", sep="\t", dtype=str).fillna("").set_index("query_protein_id")
            self.assertEqual("native", rows.loc["Q_NATIVE", "model_class"])
            self.assertEqual("reference_rescue", rows.loc["Q_RESCUE", "model_class"])
            self.assertEqual("unmapped", rows.loc["Q_UNMAPPED", "model_class"])
            self.assertEqual({"not_promoted"}, set(rows["promotion_status"]))
            self.assertEqual("UniProt:Q_RESCUE", rows.loc["Q_RESCUE", "source_id"])
            manifest = json.loads((out / "annotation_rescue_manifest.json").read_text(encoding="utf-8"))
            self.assertFalse(manifest["external_sequence_submission_this_run"])
            self.assertEqual("none", manifest["ranking_effect"])

    def test_incomplete_import_never_converts_absent_query_to_biological_unmapped(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta, candidates, queries, metadata = self._fixture(root)
            projection = root / "partial.gff3"
            projection.write_text(
                "chr1__1_400\tminiprot\tmRNA\t55\t125\t.\t+\t.\tID=MP1;Rank=1;Identity=0.9;Target=Q_NATIVE 1 20\n"
                "chr1__1_400\tminiprot\tCDS\t55\t125\t.\t+\t0\tParent=MP1\n",
                encoding="utf-8",
            )
            out = root / "out"
            self.assertEqual(0, main(["annotation-rescue", "--fasta", str(fasta), "--candidates", str(candidates),
                                      "--query-proteins", str(queries), "--query-metadata", str(metadata), "--interval", "chr1:1-400",
                                      "--projection-gff", str(projection), "--out-dir", str(out)]))
            rows = pd.read_csv(out / "annotation_rescue.tsv", sep="\t", dtype=str).fillna("").set_index("query_protein_id")
            self.assertEqual("not_evaluated", rows.loc["Q_UNMAPPED", "model_class"])
            self.assertEqual("imported_partial", rows.loc["Q_UNMAPPED", "provider_status"])

    def test_missing_local_tool_is_provider_failure_not_no_hit(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta, candidates, queries, metadata = self._fixture(root)
            out = root / "out"
            self.assertEqual(0, main(["annotation-rescue", "--fasta", str(fasta), "--candidates", str(candidates),
                                      "--query-proteins", str(queries), "--query-metadata", str(metadata), "--interval", "chr1:1-400",
                                      "--execute", "--miniprot", "definitely_missing_miniprot", "--out-dir", str(out)]))
            rows = pd.read_csv(out / "annotation_rescue.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"not_evaluated"}, set(rows["model_class"]))
            self.assertEqual({"unavailable"}, set(rows["provider_status"]))
            self.assertFalse((out / "miniprot_projection.gff3").exists())

    def test_failed_execution_discards_partial_gff_and_never_reports_unmapped(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta, candidates, queries, metadata = self._fixture(root)
            out = root / "out"

            def fake_run(command: list[str], **kwargs: object) -> SimpleNamespace:
                if "--version" in command:
                    return SimpleNamespace(returncode=0, stdout="miniprot-0.17\n", stderr="")
                kwargs["stdout"].write("partial output that must not survive\n")
                return SimpleNamespace(returncode=2, stdout="", stderr="simulated execution failure")

            with patch("candidate_qtl.annotation_rescue.shutil.which", return_value="/fake/miniprot"), patch("candidate_qtl.annotation_rescue.subprocess.run", side_effect=fake_run):
                self.assertEqual(0, main(["annotation-rescue", "--fasta", str(fasta), "--candidates", str(candidates),
                                          "--query-proteins", str(queries), "--query-metadata", str(metadata), "--interval", "chr1:1-400",
                                          "--execute", "--out-dir", str(out)]))
            rows = pd.read_csv(out / "annotation_rescue.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"not_evaluated"}, set(rows["model_class"]))
            self.assertEqual({"failed"}, set(rows["provider_status"]))
            self.assertFalse((out / "miniprot_projection.gff3").exists())
            self.assertFalse((out / "miniprot_projection.gff3.tmp").exists())

    def test_quality_tier_excludes_boundary_and_multihit_and_review_is_bounded(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta, candidates, queries, metadata = self._fixture(root)
            projection = root / "projection.gff3"
            projection.write_text(
                "chr1__1_400\tminiprot\tmRNA\t1\t60\t.\t+\t.\tID=M1;Rank=1;Identity=0.9;Target=Q_RESCUE 1 20\n"
                "chr1__1_400\tminiprot\tCDS\t1\t60\t.\t+\t0\tParent=M1\n"
                "chr1__1_400\tminiprot\tmRNA\t150\t209\t.\t+\t.\tID=M2;Rank=2;Identity=0.9;Target=Q_RESCUE 1 20\n"
                "chr1__1_400\tminiprot\tCDS\t150\t209\t.\t+\t0\tParent=M2\n",
                encoding="utf-8",
            )
            out = root / "out"
            self.assertEqual(0, main(["annotation-rescue", "--fasta", str(fasta), "--candidates", str(candidates),
                                      "--query-proteins", str(queries), "--query-metadata", str(metadata), "--interval", "chr1:1-400",
                                      "--projection-gff", str(projection), "--projection-complete", "--review-limit", "1", "--out-dir", str(out)]))
            rows = pd.read_csv(out / "annotation_rescue.tsv", sep="\t", dtype=str).fillna("")
            self.assertNotIn("high", set(rows[rows["query_protein_id"].eq("Q_RESCUE")]["quality_tier"]))
            self.assertLessEqual(len(pd.read_csv(out / "annotation_rescue_review.tsv", sep="\t")), 1)
            self.assertGreater(len(pd.read_csv(out / "annotation_rescue_review_overflow.tsv", sep="\t")), 0)

    def test_case_run_rescue_lane_keeps_rank_and_score_contributions_byte_identical(self) -> None:
        fixture = Path(__file__).parent / "fixtures" / "case_run"
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); work = root / "case"; shutil.copytree(fixture, work)
            candidates = pd.read_csv(work / "candidates.tsv", sep="\t", dtype=str).fillna("")
            candidates["strand"] = "+"
            candidates.to_csv(work / "candidates.tsv", sep="\t", index=False, lineterminator="\n")
            (work / "genome.fa").write_text(">chr1\n" + "ACGT" * 250 + "\n", encoding="utf-8")
            (work / "queries.faa").write_text(">Q1\n" + "M" * 20 + "\n", encoding="utf-8")
            (work / "query_metadata.tsv").write_text(
                "query_protein_id\tspecies\tassembly_id\tannotation_release\tsource_id\tquery_role\tselection_method\treview_status\n"
                "Q1\tPisum sativum\tZW6\tr1\tUniProt:Q1\tfunctional_homolog\tautomated_sequence_backed_search\treviewed\n", encoding="utf-8",
            )
            (work / "projection.gff3").write_text(
                "chr1__1_1000\tminiprot\tmRNA\t700\t760\t.\t+\t.\tID=MP1;Rank=1;Identity=0.9;Target=Q1 1 20\n"
                "chr1__1_1000\tminiprot\tCDS\t700\t760\t.\t+\t0\tParent=MP1\n", encoding="utf-8",
            )
            resources = yaml.safe_load((work / "resources.yaml").read_text(encoding="utf-8"))
            resources["resources"].update({
                "genome":{"type":"fasta", "path":"genome.fa"}, "rescue_queries":{"type":"protein_fasta", "path":"queries.faa"},
                "rescue_metadata":{"type":"query_metadata", "path":"query_metadata.tsv"}, "projection":{"type":"miniprot_gff3", "path":"projection.gff3"},
            })
            (work / "resources.yaml").write_text(yaml.safe_dump(resources, sort_keys=False), encoding="utf-8")
            base_case = yaml.safe_load((work / "case.yaml").read_text(encoding="utf-8"))
            base_case["target"]["fasta_resource"] = "genome"
            baseline_case = work / "baseline.yaml"; baseline_case.write_text(yaml.safe_dump(base_case, sort_keys=False), encoding="utf-8")
            assisted_case = json.loads(json.dumps(base_case))
            assisted_case["inputs"].update({"annotation_rescue_query_proteins_resource":"rescue_queries", "annotation_rescue_query_metadata_resource":"rescue_metadata", "annotation_rescue_projection_gff_resource":"projection"})
            assisted_case["annotation_rescue"] = {"projection_complete":True, "review_limit":20}
            rescue_case = work / "rescue.yaml"; rescue_case.write_text(yaml.safe_dump(assisted_case, sort_keys=False), encoding="utf-8")
            base_out, rescue_out = root / "baseline", root / "rescue"
            self.assertEqual(0, main(["case-run", "--case", str(baseline_case), "--resources", str(work / "resources.yaml"), "--out-dir", str(base_out)]))
            self.assertEqual(0, main(["case-run", "--case", str(rescue_case), "--resources", str(work / "resources.yaml"), "--out-dir", str(rescue_out)]))
            self.assertEqual((base_out / "rank" / "ranked_candidates.tsv").read_bytes(), (rescue_out / "rank" / "ranked_candidates.tsv").read_bytes())
            self.assertEqual((base_out / "rank" / "score_contributions.tsv").read_bytes(), (rescue_out / "rank" / "score_contributions.tsv").read_bytes())
            self.assertTrue((rescue_out / "annotation_rescue" / "annotation_rescue.tsv").exists())

    def test_case_run_rescue_failure_degrades_to_review_and_preserves_rank(self) -> None:
        fixture = Path(__file__).parent / "fixtures" / "case_run"
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); work = root / "case"; shutil.copytree(fixture, work)
            candidates = pd.read_csv(work / "candidates.tsv", sep="\t", dtype=str).fillna("")
            candidates["strand"] = "+"
            candidates.to_csv(work / "candidates.tsv", sep="\t", index=False, lineterminator="\n")
            (work / "genome.fa").write_text(">chr1\n" + "ACGT" * 250 + "\n", encoding="utf-8")
            (work / "queries.faa").write_text(">Q1\n" + "M" * 20 + "\n", encoding="utf-8")
            (work / "query_metadata.tsv").write_text(
                "query_protein_id\tspecies\tassembly_id\tannotation_release\tsource_id\tquery_role\tselection_method\treview_status\n"
                "Q1\tPisum sativum\tZW6\tr1\tUniProt:Q1\tfunctional_homolog\tautomated_sequence_backed_search\treviewed\n", encoding="utf-8",
            )
            # Wrong interval-local seqid: annotation_rescue must raise, and the
            # deterministic pipeline must still finish with an unchanged rank.
            (work / "projection.gff3").write_text(
                "chr1__1_999\tminiprot\tmRNA\t700\t760\t.\t+\t.\tID=MP1;Rank=1;Identity=0.9;Target=Q1 1 20\n"
                "chr1__1_999\tminiprot\tCDS\t700\t760\t.\t+\t0\tParent=MP1\n", encoding="utf-8",
            )
            resources = yaml.safe_load((work / "resources.yaml").read_text(encoding="utf-8"))
            resources["resources"].update({
                "genome":{"type":"fasta", "path":"genome.fa"}, "rescue_queries":{"type":"protein_fasta", "path":"queries.faa"},
                "rescue_metadata":{"type":"query_metadata", "path":"query_metadata.tsv"}, "projection":{"type":"miniprot_gff3", "path":"projection.gff3"},
            })
            (work / "resources.yaml").write_text(yaml.safe_dump(resources, sort_keys=False), encoding="utf-8")
            base_case = yaml.safe_load((work / "case.yaml").read_text(encoding="utf-8"))
            base_case["target"]["fasta_resource"] = "genome"
            baseline_case = work / "baseline.yaml"; baseline_case.write_text(yaml.safe_dump(base_case, sort_keys=False), encoding="utf-8")
            broken_case = json.loads(json.dumps(base_case))
            broken_case["inputs"].update({"annotation_rescue_query_proteins_resource":"rescue_queries", "annotation_rescue_query_metadata_resource":"rescue_metadata", "annotation_rescue_projection_gff_resource":"projection"})
            broken_case["annotation_rescue"] = {"projection_complete":True, "review_limit":20}
            rescue_case = work / "rescue.yaml"; rescue_case.write_text(yaml.safe_dump(broken_case, sort_keys=False), encoding="utf-8")
            base_out, rescue_out = root / "baseline", root / "rescue"
            self.assertEqual(0, main(["case-run", "--case", str(baseline_case), "--resources", str(work / "resources.yaml"), "--out-dir", str(base_out)]))
            self.assertEqual(0, main(["case-run", "--case", str(rescue_case), "--resources", str(work / "resources.yaml"), "--out-dir", str(rescue_out)]))
            self.assertEqual((base_out / "rank" / "ranked_candidates.tsv").read_bytes(), (rescue_out / "rank" / "ranked_candidates.tsv").read_bytes())
            self.assertEqual((base_out / "rank" / "score_contributions.tsv").read_bytes(), (rescue_out / "rank" / "score_contributions.tsv").read_bytes())
            self.assertFalse((rescue_out / "annotation_rescue" / "annotation_rescue.tsv").exists())
            reviews = pd.read_csv(rescue_out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("annotation_rescue_failed", set(reviews["reason"]))


if __name__ == "__main__":
    unittest.main()
