"""Workstation-local Route C provider and frozen reference indexes (issue #217)."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import build_parser, main
from candidate_qtl.core import sha256
from candidate_qtl.retrieval_local import INDEX_MANIFEST_NAME, MODEL_VERSION
from candidate_qtl.retrieval_route import validate_retrieval_frame

CORPUS_A = "\n".join([
    "reference_gene_id\treference_species\tfamily_id\torthogroup_id\tdescription\tmutant_phenotype\tpathway\tsource_id",
    "AtSBE1\tArabidopsis thaliana\tSBE\tOG0001\tstarch branching enzyme I\twrinkled seed with reduced starch\tstarch_biosynthesis\tPMID:1978859",
    "AtCLV1\tArabidopsis thaliana\tCLV\tOG0002\tCLAVATA1 receptor kinase\tfasciated stem and enlarged meristem\tmeristem_maintenance\tPMID:8402901",
    "AtNOFAM\tArabidopsis thaliana\t\t\tuncharacterised protein\t\t\tPMID:9999999",
    "AtNOSOURCE\tArabidopsis thaliana\tPHY\tOG0003\tphytochrome A photoreceptor\t\tlight_signalling\t",
]) + "\n"

CORPUS_B = "\n".join([
    "reference_gene_id\treference_species\tfamily_id\torthogroup_id\tdescription\tmutant_phenotype\tsource_id",
    "MtSBE1\tMedicago truncatula\tSBE\tOG0001\tstarch branching enzyme\twrinkled seed phenotype\tPMID:1978860",
    "MtWOX\tMedicago truncatula\tWOX\tOG0004\tWUSCHEL related homeobox\tmeristem termination\tPMID:1978861",
]) + "\n"


def _write(root: Path, name: str, content: str) -> Path:
    path = root / name
    path.write_text(content, encoding="utf-8")
    return path


def _build(root: Path, corpus: Path, version: str, out_name: str = "index") -> Path:
    out = root / out_name
    assert main(["retrieval-index-build", "--corpus", str(corpus), "--out-dir", str(out),
                 "--database-version", version]) == 0
    return out / version / INDEX_MANIFEST_NAME


class IndexBuildTests(unittest.TestCase):
    def test_index_is_frozen_versioned_and_source_backed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest_path = _build(root, _write(root, "corpus_a.tsv", CORPUS_A), "plant_db_2026-08")
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual("plant_db_2026-08", manifest["database_version"])
            self.assertEqual(MODEL_VERSION, manifest["model_version"])
            self.assertTrue(manifest["built_at_utc"].endswith("+00:00"))
            self.assertFalse(manifest["interval_information_in_index"])
            self.assertEqual(2, manifest["document_count"])
            self.assertEqual(2, manifest["excluded_row_count"])
            self.assertEqual({"missing_family_and_orthogroup": 1, "function_not_source_backed": 1},
                             manifest["exclusion_reasons"])
            for name, digest in manifest["artifact_sha256"].items():
                self.assertEqual(digest, sha256(manifest_path.parent / name), name)
            excluded = pd.read_csv(manifest_path.parent / "excluded_rows.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"AtNOFAM", "AtNOSOURCE"}, set(excluded["reference_gene_id"]))
            self.assertTrue(all(manifest["corpus_inputs"][key]["sha256"] for key in manifest["corpus_inputs"]))

    def test_rebuilding_the_same_corpus_is_byte_identical(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            corpus = _write(root, "corpus_a.tsv", CORPUS_A)
            first = json.loads(_build(root, corpus, "v1", "one").read_text(encoding="utf-8"))
            second = json.loads(_build(root, corpus, "v1", "two").read_text(encoding="utf-8"))
            self.assertEqual(first["artifact_sha256"], second["artifact_sha256"])

    def test_provider_availability_is_recorded_not_required(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = json.loads(_build(root, _write(root, "c.tsv", CORPUS_A), "v1").read_text(encoding="utf-8"))
            provider = manifest["provider"]
            self.assertEqual("local_lexical_tfidf", provider["active_provider"])
            self.assertFalse(provider["active_provider_requires_third_party"])
            self.assertEqual({"sentence_transformers", "transformers", "faiss", "chromadb"},
                             set(provider["optional_backends"]))


class RetrievalRunTests(unittest.TestCase):
    def test_phase_one_has_no_interval_argument(self) -> None:
        parser = build_parser()
        for command in ("retrieval-run", "retrieval-index-build"):
            action = [item for item in parser._subparsers._group_actions if hasattr(item, "choices")][0]
            options = set()
            for sub_action in action.choices[command]._actions:
                options.update(sub_action.option_strings)
            for forbidden in ("--interval", "--atlas", "--candidates", "--interval-genes", "--gff3"):
                self.assertNotIn(forbidden, options, f"{command} must not accept {forbidden}")

    def test_retrieval_output_satisfies_the_route_c_contract(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index = _build(root, _write(root, "corpus_a.tsv", CORPUS_A), "plant_db_2026-08")
            out = root / "run"
            self.assertEqual(0, main(["retrieval-run", "--index", str(index), "--phenotype",
                                      "wrinkled seed with reduced starch content", "--out-dir", str(out),
                                      "--case-id", "demo"]))
            retrieval = pd.read_csv(out / "route_c_local_retrieval.tsv", sep="\t", dtype=str).fillna("")
            validate_retrieval_frame(retrieval, "local retrieval")
            self.assertEqual("AtSBE1", retrieval.iloc[0]["reference_gene_id"])
            status = json.loads((out / "route_c_local_retrieval_status.json").read_text(encoding="utf-8"))
            self.assertFalse(status["interval_information_in_model_input"])
            self.assertEqual([], status["interval_arguments"])
            self.assertEqual("retrieved", status["retrieval_status"])
            self.assertTrue(status["index_verified_before_retrieval"])
            self.assertEqual("none", status["rank_contribution"])
            self.assertEqual("none", status["external_transmission"])
            self.assertEqual(status["retrieval_sha256"], sha256(out / "route_c_local_retrieval.tsv"))

    def test_frozen_retrieval_is_importable_and_then_matched_separately(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index = _build(root, _write(root, "corpus_a.tsv", CORPUS_A), "plant_db_2026-08")
            run_dir = root / "run"
            self.assertEqual(0, main(["retrieval-run", "--index", str(index), "--phenotype",
                                      "fasciated stem with enlarged meristem", "--out-dir", str(run_dir)]))
            freeze = root / "phase1"
            self.assertEqual(0, main(["retrieval-import", "--retrieval", str(run_dir / "route_c_local_retrieval.tsv"),
                                      "--out-dir", str(freeze), "--phenotype", "fasciated stem"]))
            frozen = pd.read_csv(freeze / "route_c_retrieval_frozen.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("AtCLV1", frozen.iloc[0]["reference_gene_id"])
            atlas = _write(root, "atlas.tsv", "\n".join([
                "gene_id\tfamily_id\torthogroup\torthology_reference_gene",
                "Psat_CLV1\tCLV\tOG0002\tAtCLV1",
            ]) + "\n")
            phase2 = root / "phase2"
            self.assertEqual(0, main(["retrieval-interval-match", "--retrieval", str(freeze / "route_c_retrieval_frozen.tsv"),
                                      "--atlas", str(atlas), "--out-dir", str(phase2)]))
            match_status = json.loads((phase2 / "route_c_interval_match_status.json").read_text(encoding="utf-8"))
            self.assertTrue(match_status["retrieval_frozen_before_interval"])
            self.assertEqual(["Psat_CLV1"], match_status["matched_interval_genes"])

    def test_no_matching_token_is_not_evaluated_rather_than_a_miss(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index = _build(root, _write(root, "corpus_a.tsv", CORPUS_A), "v1")
            out = root / "run"
            self.assertEqual(0, main(["retrieval-run", "--index", str(index), "--phenotype",
                                      "zzzqqq wwwxxx", "--out-dir", str(out)]))
            status = json.loads((out / "route_c_local_retrieval_status.json").read_text(encoding="utf-8"))
            self.assertEqual(0, status["retrieval_count"])
            self.assertEqual("not_evaluated_no_indexed_token_matched_the_phenotype", status["retrieval_status"])

    def test_an_empty_retrieval_imports_and_matches_as_not_evaluated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index = _build(root, _write(root, "corpus_a.tsv", CORPUS_A), "v1")
            run_dir = root / "run"
            self.assertEqual(0, main(["retrieval-run", "--index", str(index), "--phenotype",
                                      "zzzqqq wwwxxx", "--out-dir", str(run_dir)]))
            freeze = root / "phase1"
            self.assertEqual(0, main(["retrieval-import", "--retrieval", str(run_dir / "route_c_local_retrieval.tsv"),
                                      "--out-dir", str(freeze), "--phenotype", "zzzqqq wwwxxx"]))
            status = json.loads((freeze / "route_c_retrieval_status.json").read_text(encoding="utf-8"))
            self.assertEqual(0, status["retrieval_count"])
            self.assertEqual("not_evaluated_no_indexed_token_matched_the_phenotype", status["retrieval_status"])
            frozen = pd.read_csv(freeze / "route_c_retrieval_frozen.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(frozen.empty)
            atlas = _write(root, "atlas.tsv", "\n".join([
                "gene_id\tfamily_id\torthogroup\torthology_reference_gene",
                "Psat_CLV1\tCLV\tOG0002\tAtCLV1",
            ]) + "\n")
            phase2 = root / "phase2"
            self.assertEqual(0, main(["retrieval-interval-match", "--retrieval", str(freeze / "route_c_retrieval_frozen.tsv"),
                                      "--atlas", str(atlas), "--out-dir", str(phase2)]))
            match_status = json.loads((phase2 / "route_c_interval_match_status.json").read_text(encoding="utf-8"))
            self.assertTrue(match_status["retrieval_frozen_before_interval"])
            self.assertEqual(0, match_status["matched_interval_gene_count"])
            self.assertEqual("empty_frozen_retrieval", match_status["zero_match_reason"])

    def test_a_changed_index_artifact_aborts_retrieval(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index = _build(root, _write(root, "corpus_a.tsv", CORPUS_A), "v1")
            documents = index.parent / "documents.tsv"
            documents.write_text(documents.read_text(encoding="utf-8").replace("AtSBE1", "AtEDITED"), encoding="utf-8")
            with self.assertRaises(ValueError) as context:
                main(["retrieval-run", "--index", str(index), "--phenotype", "starch", "--out-dir", str(root / "run")])
            self.assertIn("do not match", str(context.exception))

    def test_two_independent_databases_merge_by_within_database_rank(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index_a = _build(root, _write(root, "corpus_a.tsv", CORPUS_A), "arabidopsis_2026-08", "index_a")
            index_b = _build(root, _write(root, "corpus_b.tsv", CORPUS_B), "medicago_2026-08", "index_b")
            runs = []
            for name, index in (("a", index_a), ("b", index_b)):
                out = root / f"run_{name}"
                self.assertEqual(0, main(["retrieval-run", "--index", str(index), "--phenotype",
                                          "wrinkled seed with reduced starch", "--out-dir", str(out)]))
                runs.append(str(out / "route_c_local_retrieval.tsv"))
            merged_dir = root / "merged"
            self.assertEqual(0, main(["retrieval-merge", "--retrieval", runs[0], "--retrieval", runs[1],
                                      "--out-dir", str(merged_dir), "--top-k", "10"]))
            merged = pd.read_csv(merged_dir / "route_c_retrieval_merged.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"arabidopsis_2026-08", "medicago_2026-08"}, set(merged["database_version"]))
            status = json.loads((merged_dir / "route_c_retrieval_merge_status.json").read_text(encoding="utf-8"))
            self.assertFalse(status["cross_database_raw_scores_compared"])


if __name__ == "__main__":
    unittest.main()
