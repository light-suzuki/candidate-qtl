"""Semantic / neural-adapter / hybrid-fusion Route C providers (issue #225)."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import build_parser, main
from candidate_qtl.core import sha256
from candidate_qtl.retrieval_local import MODEL_VERSION as LEXICAL_MODEL_VERSION
from candidate_qtl.retrieval_route import RETRIEVAL_COLUMNS, retrieval_import, validate_retrieval_frame
from candidate_qtl.retrieval_semantic import (
    LSA_RANK_K,
    MODALITY,
    MODEL_NAME,
    MODEL_VERSION,
    MODALITY_NEURAL,
    MODEL_NAME_NEURAL,
    MODEL_VERSION_NEURAL,
    NOT_EVALUATED_PROVIDER_UNAVAILABLE,
    RRF_K,
    neural_provider_available,
    retrieval_hybrid_merge,
    retrieval_neural_index_build,
    retrieval_neural_run,
)
from candidate_qtl.route_c_benchmark import route_c_provider_comparison

CORPUS = "\n".join([
    "reference_gene_id\treference_species\tfamily_id\torthogroup_id\tdescription\tmutant_phenotype\tpathway\tsource_id",
    "AtSBE1\tArabidopsis thaliana\tSBE\tOG0001\tstarch branching enzyme I\twrinkled seed with reduced starch\tstarch_biosynthesis\tPMID:1978859",
    "AtCLV1\tArabidopsis thaliana\tCLV\tOG0002\tCLAVATA1 receptor kinase\tfasciated stem and enlarged meristem\tmeristem_maintenance\tPMID:8402901",
    "MtSBE1\tMedicago truncatula\tSBE\tOG0001\tstarch branching enzyme\twrinkled seed phenotype low starch\tstarch_biosynthesis\tPMID:1978860",
    "AtWOX\tArabidopsis thaliana\tWOX\tOG0004\tWUSCHEL related homeobox\tmeristem termination defect\tmeristem_maintenance\tPMID:1978861",
]) + "\n"


def _write(root: Path, name: str, content: str) -> Path:
    path = root / name
    path.write_text(content, encoding="utf-8")
    return path


def _build_lsa(root: Path, corpus: Path, version: str, out_name: str = "lsa_index") -> Path:
    out = root / out_name
    assert main(["retrieval-index-build", "--provider", "lsa", "--corpus", str(corpus),
                 "--out-dir", str(out), "--database-version", version]) == 0
    from candidate_qtl.retrieval_semantic import INDEX_MANIFEST_NAME
    return out / version / INDEX_MANIFEST_NAME


class LsaIndexBuildTests(unittest.TestCase):
    def test_index_is_frozen_versioned_and_source_backed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest_path = _build_lsa(root, _write(root, "corpus.tsv", CORPUS), "plant_db_2026-08")
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual("plant_db_2026-08", manifest["database_version"])
            self.assertEqual(MODEL_VERSION, manifest["model_version"])
            self.assertEqual(MODALITY, manifest["modality"])
            self.assertNotEqual(MODEL_VERSION, LEXICAL_MODEL_VERSION)
            self.assertTrue(manifest["built_at_utc"].endswith("+00:00"))
            self.assertFalse(manifest["interval_information_in_index"])
            self.assertEqual(4, manifest["document_count"])
            for name, digest in manifest["artifact_sha256"].items():
                self.assertEqual(digest, sha256(manifest_path.parent / name), name)
            # Honesty: latent-semantic, explicitly not a neural embedding.
            self.assertIn("not a neural embedding", " ".join(manifest["safety"]))
            self.assertEqual(LSA_RANK_K, manifest["preprocessing"]["lsa_rank_k_configured"])
            self.assertLessEqual(manifest["preprocessing"]["lsa_rank_k_effective"], LSA_RANK_K)

    def test_rebuilding_the_same_corpus_is_byte_identical(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            corpus = _write(root, "corpus.tsv", CORPUS)
            first = json.loads(_build_lsa(root, corpus, "v1", "one").read_text(encoding="utf-8"))
            second = json.loads(_build_lsa(root, corpus, "v1", "two").read_text(encoding="utf-8"))
            self.assertEqual(first["artifact_sha256"], second["artifact_sha256"])

    def test_sign_canonicalisation_forces_the_peak_term_loading_positive(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest_path = _build_lsa(root, _write(root, "corpus.tsv", CORPUS), "v1")
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            latent_path = manifest_path.parent / "lsa_latent_index.json"
            latent = json.loads(latent_path.read_text(encoding="utf-8"))
            for component_row in latent["term_loadings"]:
                peak_index = max(range(len(component_row)), key=lambda i: abs(component_row[i]))
                self.assertGreaterEqual(component_row[peak_index], 0.0)
            self.assertIn("sign_canonicalisation_rule", manifest["preprocessing"])

    def test_provider_availability_records_numpy_dependency_honestly(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = json.loads(_build_lsa(root, _write(root, "c.tsv", CORPUS), "v1").read_text(encoding="utf-8"))
            provider = manifest["provider"]
            self.assertEqual("local_lsa_svd", provider["active_provider"])
            self.assertIn("numpy", provider["third_party_dependency"])
            self.assertIn("sentence_transformers", provider["optional_neural_backends"])


class LsaRetrievalRunTests(unittest.TestCase):
    def test_provider_flag_has_no_interval_argument(self) -> None:
        parser = build_parser()
        action = [item for item in parser._subparsers._group_actions if hasattr(item, "choices")][0]
        for command in ("retrieval-run", "retrieval-index-build"):
            options = set()
            for sub_action in action.choices[command]._actions:
                options.update(sub_action.option_strings)
            self.assertIn("--provider", options)
            for forbidden in ("--interval", "--atlas", "--candidates", "--interval-genes", "--gff3"):
                self.assertNotIn(forbidden, options, f"{command} must not accept {forbidden}")

    def test_retrieval_output_satisfies_the_route_c_contract(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index = _build_lsa(root, _write(root, "corpus.tsv", CORPUS), "plant_db_2026-08")
            out = root / "run"
            self.assertEqual(0, main(["retrieval-run", "--provider", "lsa", "--index", str(index), "--phenotype",
                                      "wrinkled seed with reduced starch content", "--out-dir", str(out),
                                      "--case-id", "demo"]))
            retrieval = pd.read_csv(out / "route_c_lsa_retrieval.tsv", sep="\t", dtype=str).fillna("")
            validate_retrieval_frame(retrieval, "lsa retrieval")
            self.assertEqual({MODEL_NAME}, set(retrieval["model_name"]))
            status = json.loads((out / "route_c_lsa_retrieval_status.json").read_text(encoding="utf-8"))
            self.assertFalse(status["interval_information_in_model_input"])
            self.assertEqual([], status["interval_arguments"])
            self.assertEqual("none", status["rank_contribution"])
            self.assertEqual("latent-semantic (TF-IDF + truncated SVD), not a neural embedding",
                             status["provider_type_declaration"])
            self.assertEqual(status["retrieval_sha256"], sha256(out / "route_c_lsa_retrieval.tsv"))

    def test_a_changed_index_artifact_aborts_retrieval(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index = _build_lsa(root, _write(root, "corpus.tsv", CORPUS), "v1")
            documents = index.parent / "lsa_documents.tsv"
            documents.write_text(documents.read_text(encoding="utf-8").replace("AtSBE1", "AtEDITED"), encoding="utf-8")
            with self.assertRaises(ValueError) as context:
                main(["retrieval-run", "--provider", "lsa", "--index", str(index), "--phenotype", "starch",
                      "--out-dir", str(root / "run")])
            self.assertIn("do not match", str(context.exception))


class NeuralAdapterTests(unittest.TestCase):
    """sentence_transformers is genuinely absent from this workstation's venv
    (confirmed by the same importlib.util.find_spec probe retrieval_local.py
    uses); these tests exercise exactly the path that is real here."""

    def test_backend_is_absent_on_this_workstation(self) -> None:
        self.assertFalse(neural_provider_available())

    def test_index_build_reports_unavailable_without_writing_a_fake_index(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = retrieval_neural_index_build([], str(root / "neural_index"), "v1")
            self.assertEqual(NOT_EVALUATED_PROVIDER_UNAVAILABLE, result["index_build_status"])
            self.assertEqual(0, result["document_count"])
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertFalse(manifest["provider"]["backend_importable"])
            self.assertIn("sentence_transformers", manifest["reason"])
            self.assertEqual(MODEL_NAME_NEURAL, manifest["model_name"])
            self.assertEqual(MODEL_VERSION_NEURAL, manifest["model_version"])
            self.assertEqual(MODALITY_NEURAL, manifest["modality"])

    def test_cli_run_never_crashes_and_never_fabricates_a_row(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index_out = root / "neural_index"
            self.assertEqual(0, main(["retrieval-index-build", "--provider", "neural", "--corpus",
                                      str(_write(root, "corpus.tsv", CORPUS)), "--out-dir", str(index_out),
                                      "--database-version", "v1"]))
            manifest_path = index_out / "v1" / "neural_index_manifest.json"
            run_out = root / "run"
            self.assertEqual(0, main(["retrieval-run", "--provider", "neural", "--index", str(manifest_path),
                                      "--phenotype", "fasciated stem", "--out-dir", str(run_out)]))
            status = json.loads((run_out / "route_c_neural_retrieval_status.json").read_text(encoding="utf-8"))
            self.assertEqual(NOT_EVALUATED_PROVIDER_UNAVAILABLE, status["retrieval_status"])
            self.assertEqual(0, status["retrieval_count"])
            self.assertEqual("none", status["rank_contribution"])
            retrieval = pd.read_csv(run_out / "route_c_neural_retrieval.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(list(RETRIEVAL_COLUMNS), list(retrieval.columns))
            self.assertEqual(0, len(retrieval))

    def test_never_silently_falls_back_to_another_provider(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = retrieval_neural_index_build([], str(root / "idx"), "v1")
            manifest_path = Path(result["manifest"])
            run_result = retrieval_neural_run(str(manifest_path), "some phenotype", str(root / "run"))
            self.assertEqual(NOT_EVALUATED_PROVIDER_UNAVAILABLE, run_result["retrieval_status"])
            self.assertNotEqual(MODEL_NAME, run_result.get("model_version"))
            self.assertEqual(MODEL_VERSION_NEURAL, run_result["model_version"])


def _row(gene: str, species: str, family: str, score: float, model_name: str, model_version: str,
         database_version: str, rank: int) -> dict[str, object]:
    return {
        "reference_gene_id": gene, "reference_species": species, "family_id": family, "orthogroup_id": "",
        "retrieval_score": score, "matched_modality": "text_lexical",
        "matched_description_or_phenotype": f"{gene} description", "model_name": model_name,
        "model_version": model_version, "database_version": database_version, "rank": rank,
    }


class HybridFusionTests(unittest.TestCase):
    def test_fused_order_depends_only_on_rank_not_on_incomparable_raw_scores(self) -> None:
        # Provider A: plausible cosine-like scores in [0, 1], best rank -> GeneX.
        provider_a = pd.DataFrame([
            _row("GeneX", "Arabidopsis thaliana", "FAM1", 0.9, "provider-a", "v1", "db1", 1),
            _row("GeneY", "Arabidopsis thaliana", "FAM1", 0.85, "provider-a", "v1", "db1", 2),
            _row("GeneZ", "Arabidopsis thaliana", "FAM1", 0.80, "provider-a", "v1", "db1", 3),
        ], columns=RETRIEVAL_COLUMNS)
        # Provider B: deliberately incomparable scale (large, negative, tiny),
        # with GeneZ ranked best by this provider.
        provider_b = pd.DataFrame([
            _row("GeneZ", "Arabidopsis thaliana", "FAM1", 1_000_000.0, "provider-b", "v1", "db1", 1),
            _row("GeneX", "Arabidopsis thaliana", "FAM1", -5.0, "provider-b", "v1", "db1", 2),
            _row("GeneY", "Arabidopsis thaliana", "FAM1", 0.0000001, "provider-b", "v1", "db1", 3),
        ], columns=RETRIEVAL_COLUMNS)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path_a = root / "a.tsv"
            path_b = root / "b.tsv"
            provider_a.to_csv(path_a, sep="\t", index=False, lineterminator="\n")
            provider_b.to_csv(path_b, sep="\t", index=False, lineterminator="\n")
            out = root / "hybrid"
            result = retrieval_hybrid_merge([str(path_a), str(path_b)], str(out), top_k=10)
            merged = pd.read_csv(out / "route_c_retrieval_hybrid.tsv", sep="\t", dtype=str).fillna("")

            # Naively summing the raw scores would rank GeneZ >> GeneY > GeneX
            # (dominated by provider B's huge scale). RRF over rank alone
            # ranks GeneX first (best combined rank 1+2), proving the fused
            # order is driven by rank, not by the incomparable score scales.
            self.assertEqual(["GeneX", "GeneZ", "GeneY"], list(merged["reference_gene_id"]))

            expected_x = 1 / (RRF_K + 1) + 1 / (RRF_K + 2)
            expected_z = 1 / (RRF_K + 3) + 1 / (RRF_K + 1)
            expected_y = 1 / (RRF_K + 2) + 1 / (RRF_K + 3)
            scores = {row["reference_gene_id"]: float(row["retrieval_score"]) for _, row in merged.iterrows()}
            self.assertAlmostEqual(expected_x, scores["GeneX"], places=8)
            self.assertAlmostEqual(expected_z, scores["GeneZ"], places=8)
            self.assertAlmostEqual(expected_y, scores["GeneY"], places=8)

            status = json.loads(Path(result["status"]).read_text(encoding="utf-8"))
            self.assertFalse(status["cross_provider_raw_scores_compared"])
            self.assertFalse(status["cross_provider_raw_scores_summed"])
            self.assertEqual("reciprocal_rank_fusion_by_within_provider_rank", status["fusion_policy"])
            self.assertEqual(RRF_K, status["rrf_constant_k"])
            self.assertEqual("none", status["rank_contribution"])

    def test_hybrid_output_spanning_two_databases_can_be_frozen_and_tampering_cannot(self) -> None:
        # The freeze gate refuses a multi-source retrieval unless a sidecar proves
        # it was combined by rank alone. retrieval-hybrid-merge fuses across
        # providers by within-provider rank, so its output must be freezable --
        # and must stop being freezable the moment the merged rows are edited.
        provider_a = pd.DataFrame([
            _row("GeneX", "Arabidopsis thaliana", "FAM1", 0.9, "provider-a", "v1", "db_lexical", 1),
        ], columns=RETRIEVAL_COLUMNS)
        provider_b = pd.DataFrame([
            _row("GeneY", "Arabidopsis thaliana", "FAM1", 1_000_000.0, "provider-b", "v1", "db_lsa", 1),
        ], columns=RETRIEVAL_COLUMNS)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path_a = root / "a.tsv"
            path_b = root / "b.tsv"
            provider_a.to_csv(path_a, sep="\t", index=False, lineterminator="\n")
            provider_b.to_csv(path_b, sep="\t", index=False, lineterminator="\n")
            out = root / "hybrid"
            result = retrieval_hybrid_merge([str(path_a), str(path_b)], str(out), top_k=10)
            merged_path = Path(result["merged"])
            frozen = retrieval_import(str(merged_path), str(root / "frozen"), phenotype="fasciation", case_id="caseH")
            frozen_frame = pd.read_csv(frozen["frozen"], sep="\t", dtype=str).fillna("")
            self.assertEqual({"db_lexical", "db_lsa"}, set(frozen_frame["database_version"]))

            tampered = pd.read_csv(merged_path, sep="\t", dtype=str).fillna("")
            tampered.loc[0, "retrieval_score"] = "0.123456"
            tampered.to_csv(merged_path, sep="\t", index=False, lineterminator="\n")
            with self.assertRaisesRegex(ValueError, "does not match the input"):
                retrieval_import(str(merged_path), str(root / "frozen2"), phenotype="fasciation", case_id="caseH")

    def test_multi_provider_retrieval_without_a_rank_only_sidecar_is_refused(self) -> None:
        # Two providers concatenated by hand: no sidecar, so nothing certifies the
        # ordering was produced without comparing incomparable raw scores.
        hand_merged = pd.DataFrame([
            _row("GeneY", "Arabidopsis thaliana", "FAM1", 1_000_000.0, "provider-b", "v1", "db_lsa", 1),
            _row("GeneX", "Arabidopsis thaliana", "FAM1", 0.9, "provider-a", "v1", "db_lexical", 2),
        ], columns=RETRIEVAL_COLUMNS)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "hand_merged.tsv"
            hand_merged.to_csv(path, sep="\t", index=False, lineterminator="\n")
            with self.assertRaisesRegex(ValueError, "retrieval-merge or retrieval-hybrid-merge"):
                retrieval_import(str(path), str(root / "frozen"), phenotype="fasciation", case_id="caseH")

    def test_requires_at_least_two_distinct_providers(self) -> None:
        provider_a = pd.DataFrame([
            _row("GeneX", "Arabidopsis thaliana", "FAM1", 0.9, "provider-a", "v1", "db1", 1),
        ], columns=RETRIEVAL_COLUMNS)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path_a = root / "a.tsv"
            provider_a.to_csv(path_a, sep="\t", index=False, lineterminator="\n")
            with self.assertRaises(ValueError):
                retrieval_hybrid_merge([str(path_a)], str(root / "out"))
            with self.assertRaises(ValueError):
                retrieval_hybrid_merge([str(path_a), str(path_a)], str(root / "out2"))

    def test_cli_hybrid_merge_command_round_trips(self) -> None:
        provider_a = pd.DataFrame([
            _row("GeneX", "Arabidopsis thaliana", "FAM1", 0.9, "provider-a", "v1", "db1", 1),
        ], columns=RETRIEVAL_COLUMNS)
        provider_b = pd.DataFrame([
            _row("GeneX", "Arabidopsis thaliana", "FAM1", 42.0, "provider-b", "v1", "db1", 1),
        ], columns=RETRIEVAL_COLUMNS)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path_a = root / "a.tsv"
            path_b = root / "b.tsv"
            provider_a.to_csv(path_a, sep="\t", index=False, lineterminator="\n")
            provider_b.to_csv(path_b, sep="\t", index=False, lineterminator="\n")
            out = root / "hybrid"
            self.assertEqual(0, main(["retrieval-hybrid-merge", "--retrieval", str(path_a),
                                      "--retrieval", str(path_b), "--out-dir", str(out)]))
            self.assertTrue((out / "route_c_retrieval_hybrid.tsv").exists())


class ProviderComparisonTests(unittest.TestCase):
    def _manifest(self, root: Path, name: str, retrieval_tsv: Path, truth_family: str) -> Path:
        manifest = root / name
        manifest.write_text(
            "case_id\tretrieval_tsv\ttruth_family_id\n"
            f"case1\t{retrieval_tsv}\t{truth_family}\n",
            encoding="utf-8",
        )
        return manifest

    def test_side_by_side_recall_never_pools_raw_scores(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # Filler family labels deliberately avoid a trailing-digit form
            # (family_labels._TRAILING_COPY treats e.g. "FAM1"/"FAM2" as the
            # same numbered-copy family) so they never accidentally alias
            # the truth family "FAM1".
            # Provider "good": truth family FAM1 at rank 1 -> hits F@5/10/20.
            good_rows = [_row("GeneA", "Arabidopsis thaliana", "FAM1", 0.9, "provider-good", "v1", "db1", 1)]
            good_rows += [_row(f"Filler{i}", "Arabidopsis thaliana", f"OTHERFAM{chr(65 + i)}", 0.1, "provider-good", "v1", "db1", i + 2)
                          for i in range(4)]
            # Provider "weak": truth family only reachable at rank 15 -> misses F@5/10, hits F@20.
            weak_rows = [_row(f"Filler{i}", "Arabidopsis thaliana", f"OTHERFAM{chr(65 + i)}", 0.5, "provider-weak", "v1", "db1", i + 1)
                         for i in range(14)]
            weak_rows.append(_row("GeneA", "Arabidopsis thaliana", "FAM1", 0.05, "provider-weak", "v1", "db1", 15))
            good_path = root / "good.tsv"
            weak_path = root / "weak.tsv"
            pd.DataFrame(good_rows, columns=RETRIEVAL_COLUMNS).to_csv(good_path, sep="\t", index=False, lineterminator="\n")
            pd.DataFrame(weak_rows, columns=RETRIEVAL_COLUMNS).to_csv(weak_path, sep="\t", index=False, lineterminator="\n")

            manifest_good = self._manifest(root, "manifest_good.tsv", good_path, "FAM1")
            manifest_weak = self._manifest(root, "manifest_weak.tsv", weak_path, "FAM1")
            out = root / "comparison"
            result = route_c_provider_comparison(
                {"good": str(manifest_good), "weak": str(manifest_weak)}, str(out),
            )
            comparison = pd.read_csv(result["comparison"], sep="\t", dtype=str).fillna("")
            self.assertEqual({"good", "weak"}, set(comparison["provider_label"]))
            good_f5 = comparison[(comparison["provider_label"] == "good") & (comparison["metric"] == "family_recall_5")]
            weak_f5 = comparison[(comparison["provider_label"] == "weak") & (comparison["metric"] == "family_recall_5")]
            self.assertEqual("1.0", good_f5.iloc[0]["recall"])
            self.assertEqual("0.0", weak_f5.iloc[0]["recall"])
            weak_f20 = comparison[(comparison["provider_label"] == "weak") & (comparison["metric"] == "family_recall_20")]
            self.assertEqual("1.0", weak_f20.iloc[0]["recall"])

            manifest_out = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertFalse(manifest_out["cross_provider_raw_scores_compared"])
            self.assertFalse(manifest_out["cross_provider_raw_scores_pooled"])
            self.assertEqual(["good", "weak"], manifest_out["providers"])

    def test_requires_at_least_two_providers(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows = [_row("GeneA", "Arabidopsis thaliana", "FAM1", 0.9, "provider-good", "v1", "db1", 1)]
            path = root / "one.tsv"
            pd.DataFrame(rows, columns=RETRIEVAL_COLUMNS).to_csv(path, sep="\t", index=False, lineterminator="\n")
            manifest = self._manifest(root, "manifest.tsv", path, "FAM1")
            with self.assertRaises(ValueError):
                route_c_provider_comparison({"only": str(manifest)}, str(root / "out"))

    def test_cli_route_c_provider_comparison_round_trips(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows_a = [_row("GeneA", "Arabidopsis thaliana", "FAM1", 0.9, "provider-a", "v1", "db1", 1)]
            rows_b = [_row("GeneA", "Arabidopsis thaliana", "FAM1", 0.5, "provider-b", "v1", "db1", 1)]
            path_a = root / "a.tsv"
            path_b = root / "b.tsv"
            pd.DataFrame(rows_a, columns=RETRIEVAL_COLUMNS).to_csv(path_a, sep="\t", index=False, lineterminator="\n")
            pd.DataFrame(rows_b, columns=RETRIEVAL_COLUMNS).to_csv(path_b, sep="\t", index=False, lineterminator="\n")
            manifest_a = self._manifest(root, "manifest_a.tsv", path_a, "FAM1")
            manifest_b = self._manifest(root, "manifest_b.tsv", path_b, "FAM1")
            out = root / "comparison"
            self.assertEqual(0, main(["route-c-provider-comparison",
                                      "--provider-manifest", f"a={manifest_a}",
                                      "--provider-manifest", f"b={manifest_b}",
                                      "--out-dir", str(out)]))
            self.assertTrue((out / "route_c_provider_comparison.tsv").exists())


if __name__ == "__main__":
    unittest.main()
