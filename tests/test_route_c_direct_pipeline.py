from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.core import sha256
from candidate_qtl.route_c_direct_pipeline import (
    DIRECT_BENCHMARK_CASE_COLUMNS,
    DIRECT_PROPOSAL_COLUMNS,
    SPECIES_SCOPE_COLUMNS,
    _literature_query,
    direct_literature_search,
    direct_rescue_benchmark,
    direct_rescue_freeze,
    direct_rescue_project,
)


def _write_tsv(path: Path, rows: list[dict[str, str]], columns: tuple[str, ...]) -> Path:
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False, lineterminator="\n")
    return path


class RouteCDirectPipelineTests(unittest.TestCase):
    def _citation_cache(self, root: Path) -> Path:
        cache = root / "citation_cache.json"
        cache.write_text(json.dumps({
            "pmid": "12345678",
            "id": "12345678",
            "title": "CLV1 loss of function causes a fasciated meristem",
            "abstractText": (
                "Loss-of-function mutations in CLV1 caused an enlarged fasciated shoot apical meristem. "
                "Complementation with CLV1 rescued the phenotype."
            ),
            "journalTitle": "Plant Genetics",
            "pubYear": "2020",
            "doi": "10.1000/example",
            "authorString": "Example A",
        }), encoding="utf-8")
        return cache

    def _proposal(self) -> dict[str, str]:
        row = {column: "" for column in DIRECT_PROPOSAL_COLUMNS}
        row.update({
            "record_id": "c0001",
            "query_phenotype": "fasciated stem and enlarged meristem",
            "focal_species": "Pisum sativum",
            "reported_phenotype": "fasciated shoot apical meristem",
            "reference_species": "Arabidopsis thaliana",
            "reported_gene_id": "AT1G75820",
            "reported_gene_symbol": "CLV1",
            "species_distance": "S5_model_or_other_plant",
            "evidence_directness": "D0_causal_perturbation_or_rescue",
            "evidence_type": "loss_of_function_and_complementation",
            "source_id": "PMID:12345678",
            "reference_sequence_id": "UniProt:P93003",
        })
        return row

    def _freeze(self, root: Path) -> dict[str, object]:
        proposals = _write_tsv(root / "proposals.tsv", [self._proposal()], DIRECT_PROPOSAL_COLUMNS)
        return direct_rescue_freeze(str(proposals), str(root / "freeze"), str(self._citation_cache(root)))

    def test_search_is_interval_blind_and_supports_species_specific_similar_phenotypes(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            scopes = [
                {"query_id": "pea", "reference_species": "Pisum sativum",
                 "species_distance": "S1_same_species", "phenotype_query": "flattened fasciated stem"},
                {"query_id": "arabidopsis", "reference_species": "Arabidopsis thaliana",
                 "species_distance": "S5_model_or_other_plant", "phenotype_query": "enlarged fasciated meristem"},
            ]
            scope_path = _write_tsv(root / "scope.tsv", scopes, SPECIES_SCOPE_COLUMNS)
            cache_dir = root / "search_cache"
            cache_dir.mkdir()
            for index, scope in enumerate(scopes):
                query = _literature_query(scope["phenotype_query"], scope["reference_species"])
                (cache_dir / f"q{index}.json").write_text(json.dumps({
                    "query": query,
                    "results": [{"pmid": str(100 + index), "title": f"Direct report {index}",
                                 "abstractText": "A mutant allele caused the phenotype.", "pubYear": "2020"}],
                }), encoding="utf-8")
            result = direct_literature_search(
                "fasciated stem", "Pisum sativum", str(scope_path), str(root / "search"), str(cache_dir)
            )
            candidates = pd.read_csv(result["candidates"], sep="\t", dtype=str).fillna("")
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual({"Pisum sativum", "Arabidopsis thaliana"}, set(candidates["reference_species"]))
            self.assertTrue(manifest["interval_blind"])
            self.assertEqual("none", manifest["rank_contribution"])
            self.assertNotIn("qtl", " ".join(candidates.columns).casefold())
            self.assertTrue(manifest["candidates_sha256"])

    def test_independent_c_keeps_locus_variant_and_duplicate_reports_without_qtl(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            original = self._proposal()
            original.update(record_id="c-locus-1", reference_species="Pisum sativum",
                            species_distance="S1_same_species", reported_gene_id="",
                            reported_gene_symbol="", reference_sequence_id="",
                            reported_locus="historical-locus-A", reported_variant="")
            duplicate = {**original, "record_id": "c-locus-2", "source_id": "PMID:12345678"}
            variant = {**original, "record_id": "c-variant-3", "reported_locus":"",
                       "reported_variant": "reported-variant-A"}
            rows = [original, duplicate, variant]
            proposals = _write_tsv(root / "free-reports.tsv", rows, DIRECT_PROPOSAL_COLUMNS)
            result = direct_rescue_freeze(str(proposals), str(root / "freeze"), str(self._citation_cache(root)))
            frozen = pd.read_csv(result["frozen"], sep="\t", dtype=str).fillna("")
            self.assertEqual(["c-locus-1", "c-locus-2", "c-variant-3"], frozen.record_id.tolist())
            self.assertEqual(3, result["record_count"])
            self.assertEqual(2, (frozen.reported_locus == "historical-locus-A").sum())
            self.assertEqual(["", "", "reported-variant-A"], frozen.reported_variant.tolist())
            self.assertFalse(set(["qtl_id", "focal_gene_id", "candidate_gene_id"]) & set(frozen.columns))
            # Source/identity proof may remain unresolved; freezing is not evidence promotion.
            self.assertEqual(0, result["eligible_count"])

    def test_same_species_distance_is_checked(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            scope = _write_tsv(root / "scope.tsv", [{
                "query_id": "bad", "reference_species": "Arabidopsis thaliana",
                "species_distance": "S1_same_species", "phenotype_query": "fasciated meristem",
            }], SPECIES_SCOPE_COLUMNS)
            with self.assertRaisesRegex(ValueError, "requires reference_species to equal focal_species"):
                direct_literature_search("fasciated", "Pisum sativum", str(scope), str(root / "out"))

    def test_freeze_reuses_source_verify_and_emits_eligible_record(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            result = self._freeze(root)
            frozen = pd.read_csv(result["frozen"], sep="\t", dtype=str).fillna("")
            row = frozen.iloc[0]
            self.assertEqual("verified", row["source_verification_state"])
            self.assertEqual("verified_gene", row["gene_identity_state"])
            self.assertEqual("direct_target_trait", row["trait_relation_state"])
            self.assertEqual("eligible_direct_record", row["direct_rescue_eligibility"])
            self.assertEqual("none", row["rank_contribution"])

    def test_projection_verifies_freeze_hash_and_preserves_one_to_many(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            frozen = self._freeze(root)
            atlas = root / "atlas.tsv"
            pd.DataFrame([
                {"gene_id": "Psat.g1", "orthology_reference_gene": "AT1G75820", "chrom": "1"},
                {"gene_id": "Psat.g2", "orthology_reference_gene": "AT1G75820", "chrom": "1"},
            ]).to_csv(atlas, sep="\t", index=False, lineterminator="\n")
            projected = direct_rescue_project(
                str(frozen["frozen"]), str(frozen["manifest"]), str(atlas), "qtl1", str(root / "project")
            )
            frame = pd.read_csv(projected["projection"], sep="\t", dtype=str).fillna("")
            self.assertEqual({"Psat.g1", "Psat.g2"}, set(frame["focal_gene_id"]))
            self.assertEqual({"one_to_many"}, set(frame["ambiguity_state"]))
            self.assertEqual({"none"}, set(frame["rank_contribution"]))

            Path(frozen["frozen"]).write_text("tampered\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "hash does not match"):
                direct_rescue_project(
                    str(frozen["frozen"]), str(frozen["manifest"]), str(atlas), "qtl1", str(root / "bad")
                )

    def test_no_interval_match_is_not_outside_qtl(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            frozen = self._freeze(root)
            atlas = root / "atlas.tsv"
            pd.DataFrame([{"gene_id": "Psat.other", "orthology_reference_gene": "OTHER"}]).to_csv(
                atlas, sep="\t", index=False, lineterminator="\n"
            )
            result = direct_rescue_project(
                str(frozen["frozen"]), str(frozen["manifest"]), str(atlas), "qtl1", str(root / "project")
            )
            row = pd.read_csv(result["projection"], sep="\t", dtype=str).fillna("").iloc[0]
            self.assertEqual("no_interval_match", row["projection_state"])
            self.assertEqual("", row["focal_gene_id"])

    def test_benchmark_scores_source_and_projection_separately(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            frozen = self._freeze(root)
            atlas = root / "atlas.tsv"
            pd.DataFrame([{"gene_id": "Psat.g1", "orthology_reference_gene": "AT1G75820"}]).to_csv(
                atlas, sep="\t", index=False, lineterminator="\n"
            )
            projected = direct_rescue_project(
                str(frozen["frozen"]), str(frozen["manifest"]), str(atlas), "qtl1", str(root / "project")
            )
            case = {column: "" for column in DIRECT_BENCHMARK_CASE_COLUMNS}
            case.update({
                "case_id": "public-smoke-1", "expected_qtl_id": "qtl1",
                "expected_atlas_sha256": sha256(atlas), "direct_rescue_tsv": str(frozen["frozen"]),
                "direct_freeze_manifest": str(frozen["manifest"]), "projection_tsv": str(projected["projection"]),
                "projection_manifest": str(projected["manifest"]), "truth_reference_gene_id": "CLV1",
                "truth_focal_gene_id": "Psat.g1", "contamination_status": "clean",
            })
            manifest = _write_tsv(root / "benchmark.tsv", [case], DIRECT_BENCHMARK_CASE_COLUMNS)
            result = direct_rescue_benchmark(str(manifest), str(root / "benchmark"))
            self.assertEqual(1.0, result["metrics"]["source_recovery"]["recall"])
            self.assertEqual(1.0, result["metrics"]["projection"]["recall"])

            contaminated = dict(case)
            contaminated["case_id"] = "public-smoke-contaminated"
            contaminated["contamination_status"] = "contaminated"
            manifest = _write_tsv(
                root / "benchmark-contaminated.tsv", [case, contaminated], DIRECT_BENCHMARK_CASE_COLUMNS
            )
            result = direct_rescue_benchmark(str(manifest), str(root / "benchmark-contaminated"))
            self.assertEqual(1, result["metrics"]["source_recovery"]["evaluated"])
            self.assertEqual(1, result["metrics"]["source_recovery"]["contaminated"])
            self.assertEqual(1, result["metrics"]["projection"]["evaluated"])
            self.assertEqual(1, result["metrics"]["projection"]["contaminated"])

            Path(projected["projection"]).write_text("tampered\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "projection artifact hash mismatch"):
                direct_rescue_benchmark(str(_write_tsv(
                    root / "benchmark-tampered.tsv", [case], DIRECT_BENCHMARK_CASE_COLUMNS
                )), str(root / "benchmark-tampered"))

    def test_benchmark_rejects_unknown_contamination_and_wrong_qtl(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            frozen = self._freeze(root)
            atlas = root / "atlas.tsv"
            pd.DataFrame([{"gene_id": "Psat.g1", "orthology_reference_gene": "AT1G75820"}]).to_csv(
                atlas, sep="\t", index=False, lineterminator="\n"
            )
            projected = direct_rescue_project(
                str(frozen["frozen"]), str(frozen["manifest"]), str(atlas), "qtl1", str(root / "project")
            )
            case = {column: "" for column in DIRECT_BENCHMARK_CASE_COLUMNS}
            case.update({
                "case_id": "case-1", "expected_qtl_id": "qtl1", "expected_atlas_sha256": sha256(atlas),
                "direct_rescue_tsv": str(frozen["frozen"]), "direct_freeze_manifest": str(frozen["manifest"]),
                "projection_tsv": str(projected["projection"]), "projection_manifest": str(projected["manifest"]),
                "truth_reference_gene_id": "CLV1", "truth_focal_gene_id": "Psat.g1",
                "contamination_status": "typo",
            })
            manifest = _write_tsv(root / "bad-status.tsv", [case], DIRECT_BENCHMARK_CASE_COLUMNS)
            with self.assertRaisesRegex(ValueError, "contamination_status"):
                direct_rescue_benchmark(str(manifest), str(root / "bad-status"))

            case["contamination_status"] = "clean"
            case["expected_qtl_id"] = "qtl-other"
            manifest = _write_tsv(root / "wrong-qtl.tsv", [case], DIRECT_BENCHMARK_CASE_COLUMNS)
            with self.assertRaisesRegex(ValueError, "qtl_id mismatch"):
                direct_rescue_benchmark(str(manifest), str(root / "wrong-qtl"))


if __name__ == "__main__":
    unittest.main()
