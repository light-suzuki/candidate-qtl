"""Local Route B provider stack and reverse-lookup readiness (issue #216)."""

from __future__ import annotations

import hashlib
import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl.cli import build_parser, main
from candidate_qtl.interval_atlas import GENE_MODEL_STATUSES
from candidate_qtl.dual_route import _function_lookup
from candidate_qtl.route_b_prepare import _normalise_subject, route_b_prepare
from candidate_qtl.route_b_semantic_bundle import build_route_b_semantic_bundle

CANDIDATES = "\n".join([
    "gene_id\tchrom\tstart\tend\tdescription",
    "Psat_G1\tchr5\t1000\t3000\treceptor kinase with a meristem role",
    "Psat_G2\tchr5\t9000\t11000\tstarch branching enzyme like protein",
    "Psat_G3\tchr5\t20000\t20400\thypothetical protein",
]) + "\n"

BARE_CANDIDATES = "\n".join([
    "gene_id\tchrom\tstart\tend\tdescription",
    "Psat_G1\tchr5\t1000\t3000\t",
    "Psat_G2\tchr5\t9000\t11000\t",
]) + "\n"

ORTHOLOGY = "\n".join([
    "gene_id\treference_gene_id\torthogroup\trelationship\tprovider",
    "Psat_G1\tAtCLV1\tOG0002\tbest_hit_homology\tlocal_blastp:arabidopsis_uniprot",
    "Psat_G2\tAtSBE1\tOG0001\tbest_hit_homology\tlocal_blastp:arabidopsis_uniprot",
]) + "\n"

REFERENCE_CORPUS = "\n".join([
    "reference_gene_id\treference_species\tfamily_id\torthogroup_id\tdescription\tmutant_phenotype\tpathway\tsource_id",
    "AtCLV1\tArabidopsis thaliana\tCLV\tOG0002\tCLAVATA1 receptor kinase\tfasciated stem\tmeristem_maintenance\tPMID:8402901",
    "AtSBE1\tArabidopsis thaliana\tSBE\tOG0001\tstarch branching enzyme I\twrinkled seed\tstarch_biosynthesis\tPMID:1978859",
]) + "\n"


def _write(root: Path, name: str, content: str) -> Path:
    path = root / name
    path.write_text(content, encoding="utf-8")
    return path


class ProviderRichTests(unittest.TestCase):
    def test_focal_fallback_never_merges_distinct_reference_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "functions.tsv"
            pd.DataFrame([
                {"gene_id": "G1", "reference_gene_id": "R1", "function_description": "function A", "source_id": "PMID:1"},
                {"gene_id": "G1", "reference_gene_id": "R2", "function_description": "function B", "source_id": "PMID:2"},
            ]).to_csv(path, sep="\t", index=False)
            by_reference, by_focal = _function_lookup(str(path))
            self.assertEqual("function A", by_reference["r1"]["function_description"])
            self.assertEqual("PMID:1", by_reference["r1"]["source_id"])
            self.assertNotIn("g1", by_focal)

    def test_local_blast_autogenerates_an_interval_scoped_family_index(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "route_b"
            candidates = _write(root, "candidates.tsv", CANDIDATES.replace(
                "receptor kinase with a meristem role", "").replace("starch branching enzyme like protein", ""))
            proteins = _write(root, "proteins.fa", ">Psat_G1\nMASTPRQTVSEK\n")
            database = _write(root, "local_db", "placeholder")

            def fake_blast(_proteins, _database, out_dir, _threads, _min_identity, _min_coverage, _max_evalue):
                hits = out_dir / "route_b_local_blastp.tsv"
                hits.write_text(
                    "Psat_G1\tLotjaGi6g1v0243000.1\t80\t12\t12\t12\t1\t12\t1\t12\t1e-20\t100\t"
                    "LotjaGi6g1v0243000.1\n",
                    encoding="utf-8",
                )
                rows = [{
                    "gene_id": "Psat_G1", "reference_gene_id": "P0CLV1", "orthogroup": "",
                    "relationship": "best_hit_homology", "provider": "local_blastp:local_db",
                    "percent_identity": 80.0, "query_coverage": 1.0, "evalue": 1e-20, "bitscore": 100.0,
                }]
                return rows, {"available": True, "hits_table": str(hits), "database": str(_database),
                              "genes_with_a_usable_hit": 1, "command": "fake blastp"}

            with patch("candidate_qtl.route_b_prepare._run_blastp", side_effect=fake_blast):
                self.assertEqual(0, main([
                    "route-b-prepare", "--candidates", str(candidates), "--proteins", str(proteins),
                    "--protein-db", str(database), "--reference-function",
                    str(_write(root, "corpus.tsv", REFERENCE_CORPUS)), "--out-dir", str(out),
                    "--case-id", "auto-family",
                ]))

            manifest = json.loads((out / "route_b_preparation_manifest.json").read_text(encoding="utf-8"))
            provider = manifest["providers"]["family_index"]
            self.assertEqual("generated_from_local_blastp", provider["resolution"])
            self.assertEqual("interval_candidates", provider["index_scope"])
            self.assertFalse(provider["reusable_as_genome_wide_index"])
            self.assertTrue(Path(provider["path"]).exists())
            auto_manifest = json.loads(Path(provider["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual("interval_candidates", auto_manifest["index_scope"])
            atlas = pd.read_csv(out / "interval_gene_atlas.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("family_index", atlas.set_index("gene_id").loc["Psat_G1", "family_source"])
            atlas_manifest = json.loads((out / "atlas" / "interval_atlas_manifest.json").read_text(encoding="utf-8"))
            coverage = atlas_manifest["lineage_claim"]["family_index_search_provenance"]
            self.assertEqual("interval_candidates", coverage["index_scope"])
            self.assertFalse(coverage["reusable_as_genome_wide_index"])

            hypotheses = _write(
                root, "hypotheses.tsv",
                "hypothesis_id\tpathway_or_process\tgene_family\treference_gene_name\tproposed_stable_id\n"
                "h1\tmeristem\tCLV\tCLV1\tP0CLV1\n",
            )
            with self.assertRaisesRegex(ValueError, "requires a reference_genome family index"):
                main([
                    "pathway-interval-falsify", "--family-index", provider["path"],
                    "--hypotheses", str(hypotheses), "--candidates", str(candidates),
                    "--out-dir", str(root / "unsafe_falsify"),
                ])

    def test_explicit_family_index_is_never_replaced_by_the_fallback(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "route_b"
            family = _write(root, "family_index.tsv", "\n".join([
                "gene_id\tbest_subject\tbest_subject_symbol\tbest_subject_stem\tassignment_status",
                "Psat_G1\tAtCLV1\tCLV1\tCLV\tassigned",
            ]) + "\n")
            proteins = _write(root, "proteins.fa", ">Psat_G1\nMASTPRQTVSEK\n")
            database = _write(root, "local_db", "placeholder")

            def fake_blast(_proteins, _database, out_dir, *_args):
                hits = out_dir / "route_b_local_blastp.tsv"
                hits.write_text(
                    "Psat_G1\tsp|P0CLV1|CLV1_ARATH\t80\t12\t12\t12\t1\t12\t1\t12\t1e-20\t100\tCLV1\n",
                    encoding="utf-8",
                )
                rows = [{"gene_id": "Psat_G1", "reference_gene_id": "P0CLV1", "orthogroup": "",
                         "relationship": "best_hit_homology", "provider": "local_blastp:local_db",
                         "percent_identity": 80.0, "query_coverage": 1.0, "evalue": 1e-20, "bitscore": 100.0}]
                return rows, {"available": True, "hits_table": str(hits), "database": str(_database),
                              "genes_with_a_usable_hit": 1, "command": "fake blastp"}

            with patch("candidate_qtl.route_b_prepare._run_blastp", side_effect=fake_blast):
                self.assertEqual(0, main([
                    "route-b-prepare", "--candidates", str(_write(root, "candidates.tsv", CANDIDATES)),
                    "--proteins", str(proteins), "--protein-db", str(database),
                    "--family-index", str(family), "--out-dir", str(out),
                ]))
            manifest = json.loads((out / "route_b_preparation_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("explicit_argument", manifest["providers"]["family_index"]["resolution"])
            self.assertFalse((out / "auto_family_index").exists())

    def test_source_backed_stack_reaches_a_complete_reverse_lookup(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "route_b"
            self.assertEqual(0, main([
                "route-b-prepare",
                "--candidates", str(_write(root, "candidates.tsv", CANDIDATES)),
                "--orthology", str(_write(root, "orthology.tsv", ORTHOLOGY)),
                "--reference-function", str(_write(root, "corpus.tsv", REFERENCE_CORPUS)),
                "--out-dir", str(out), "--case-id", "demo",
            ]))
            manifest = json.loads((out / "route_b_preparation_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("source_backed_reverse_lookup", manifest["route_b_readiness"])
            self.assertEqual("complete", manifest["process_inference_completeness"])
            self.assertFalse(manifest["phenotype_used"])
            self.assertEqual(2, manifest["providers"]["reference_function"]["source_backed_rows"])
            status = json.loads(Path(manifest["outputs"]["route_b0_status"]).read_text(encoding="utf-8"))
            self.assertEqual("source_backed_reverse_lookup", status["reverse_lookup_status"])
            self.assertFalse(status["phenotype_used"])
            self.assertGreater(status["source_backed_genes"], 0)
            # the unknown gene is kept, not dropped
            self.assertEqual(1, status["unclassified_genes"])

    def test_verified_families_are_handed_off_without_a_score(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "route_b"
            self.assertEqual(0, main([
                "route-b-prepare",
                "--candidates", str(_write(root, "candidates.tsv", CANDIDATES)),
                "--orthology", str(_write(root, "orthology.tsv", ORTHOLOGY)),
                "--reference-function", str(_write(root, "corpus.tsv", REFERENCE_CORPUS)),
                "--out-dir", str(out),
            ]))
            handoff = pd.read_csv(out / "route_b_family_handoff.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(len(handoff))
            self.assertEqual({"none"}, set(handoff["rank_contribution"]))
            self.assertTrue(all(value for value in handoff["source_ids"]))
            self.assertEqual({"wet_handoff_single_interval_copy"}, set(handoff["next_action"]))

    def test_a_local_best_hit_never_claims_one_to_one_orthology(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "route_b"
            self.assertEqual(0, main([
                "route-b-prepare",
                "--candidates", str(_write(root, "candidates.tsv", CANDIDATES)),
                "--orthology", str(_write(root, "orthology.tsv", ORTHOLOGY)),
                "--reference-function", str(_write(root, "corpus.tsv", REFERENCE_CORPUS)),
                "--out-dir", str(out),
            ]))
            lookup = pd.read_csv(out / "route_b0" / "interval_reverse_lookup.tsv", sep="\t", dtype=str).fillna("")
            self.assertNotIn("exact_ortholog", set(lookup["relationship_tier"]))

    def test_uniprot_style_subject_ids_are_reduced_to_the_accession(self) -> None:
        self.assertEqual("P12345", _normalise_subject("sp|P12345|CLV1_ARATH"))
        self.assertEqual("AT1G75820", _normalise_subject("AT1G75820"))


class ReturnContractTests(unittest.TestCase):
    """Issue #374 T3: the public return dict must expose the orthology path for
    the other lanes, and the manifest must aggregate gene-model completeness."""

    def _source_backed_run(self, root: Path) -> dict[str, object]:
        out = root / "route_b"
        return route_b_prepare(
            candidates_arg=str(_write(root, "candidates.tsv", CANDIDATES)),
            out_dir_arg=str(out),
            orthology_arg=str(_write(root, "orthology.tsv", ORTHOLOGY)),
            reference_function_args=[str(_write(root, "corpus.tsv", REFERENCE_CORPUS))],
            case_id="contract",
        )

    def test_return_dict_exposes_the_orthology_path(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = self._source_backed_run(root)
            orthology = _write(root, "orthology.tsv", ORTHOLOGY)
            self.assertEqual(str(orthology), result["orthology"])
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            # the same path the manifest records, so a lane reading either
            # source sees the identical file
            self.assertEqual(result["orthology"], manifest["outputs"]["orthology"])

    def test_return_dict_reports_orthology_none_when_no_provider_exists(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "route_b"
            result = route_b_prepare(
                candidates_arg=str(_write(root, "candidates.tsv", CANDIDATES)),
                out_dir_arg=str(out), auto_local_db=False,
            )
            self.assertIsNone(result["orthology"])
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertIsNone(manifest["outputs"]["orthology"])

    def test_manifest_aggregates_gene_model_completeness_from_the_atlas(self) -> None:
        """The candidates fixture carries no protein_length / cds_length /
        translation_status columns, so every gene is no_model; the aggregation
        still reports the full vocabulary with zero counts."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = self._source_backed_run(root)
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            counts = manifest["gene_model_completeness"]
            self.assertEqual({"full_model": 0, "partial_model": 0, "no_model": 3}, counts)
            self.assertEqual(list(GENE_MODEL_STATUSES), manifest["gene_model_completeness_vocabulary"])
            # matches what the atlas itself aggregated (the atlas manifest is
            # sparse like evidence_completeness_counts; the lane manifest
            # reports the full vocabulary with explicit zero counts)
            atlas_manifest = json.loads(Path(manifest["outputs"]["interval_atlas_manifest"]).read_text(encoding="utf-8"))
            self.assertEqual({"no_model": 3}, atlas_manifest["gene_model_completeness_counts"])
            self.assertEqual(atlas_manifest["gene_model_completeness_counts"],
                             {status: value for status, value in counts.items() if value})
            # every interval gene is counted exactly once
            self.assertEqual(result["gene_count"], sum(int(value) for value in counts.values()))
            # the canonical atlas actually carries the column
            atlas = pd.read_csv(Path(manifest["outputs"]["interval_gene_atlas"]), sep="\t", dtype=str).fillna("")
            self.assertEqual({"no_model"}, set(atlas["gene_model_completeness"]))


class DegradedProviderTests(unittest.TestCase):
    def test_broad_bins_alone_are_not_an_evaluated_reverse_lookup(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "route_b"
            self.assertEqual(0, main([
                "route-b-prepare",
                "--candidates", str(_write(root, "candidates.tsv", CANDIDATES)),
                "--out-dir", str(out), "--no-auto-local-db",
            ]))
            manifest = json.loads((out / "route_b_preparation_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("not_evaluated_for_reverse_lookup", manifest["route_b_readiness"])
            status = json.loads(Path(manifest["outputs"]["route_b0_status"]).read_text(encoding="utf-8"))
            self.assertEqual("broad_bin_only", status["reverse_lookup_status"])
            # the old "did it have anything to group by" flag still says yes, which is
            # exactly the reading the strict status exists to correct
            self.assertTrue(status["semantically_evaluable"])
            reviews = pd.read_csv(out / "route_b0" / "reverse_lookup_review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(any(reason.startswith("route_b_reverse_lookup_not_complete") for reason in reviews["reason"]))

    def test_a_missing_provider_is_recorded_not_treated_as_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "route_b"
            self.assertEqual(0, main([
                "route-b-prepare",
                "--candidates", str(_write(root, "candidates.tsv", BARE_CANDIDATES)),
                "--out-dir", str(out), "--no-auto-local-db",
            ]))
            manifest = json.loads((out / "route_b_preparation_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("not_evaluated_for_reverse_lookup", manifest["route_b_readiness"])
            self.assertFalse(manifest["providers"]["reference_function"]["available"])
            self.assertFalse(manifest["providers"]["homology_search"]["available"])
            self.assertTrue(manifest["providers"]["homology_search"]["missing_provider_is_not_negative_evidence"])
            status = json.loads(Path(manifest["outputs"]["route_b0_status"]).read_text(encoding="utf-8"))
            self.assertEqual("no_provider", status["reverse_lookup_status"])
            self.assertEqual(2, status["unclassified_genes"])

    def test_enumeration_precedes_filtering(self) -> None:
        """Route B may see the phenotype -- withholding it was protecting an
        independence claim between A and B that never held. What must stay
        phenotype-free is the enumeration: the atlas."""
        parser = build_parser()
        action = [item for item in parser._subparsers._group_actions if hasattr(item, "choices")][0]

        def options(command: str) -> set[str]:
            found: set[str] = set()
            for sub_action in action.choices[command]._actions:
                found.update(sub_action.option_strings)
            return found

        self.assertIn("--phenotype", options("route-b-prepare"))
        self.assertNotIn("--phenotype", options("interval-atlas"))

    def test_the_atlas_stays_phenotype_independent_even_when_b_sees_it(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "route_b"
            self.assertEqual(0, main([
                "route-b-prepare",
                "--candidates", str(_write(root, "candidates.tsv", CANDIDATES)),
                "--orthology", str(_write(root, "orthology.tsv", ORTHOLOGY)),
                "--reference-function", str(_write(root, "corpus.tsv", REFERENCE_CORPUS)),
                "--phenotype", "fasciated stem with enlarged meristem",
                "--out-dir", str(out),
            ]))
            manifest = json.loads((out / "route_b_preparation_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(manifest["phenotype_used"])
            self.assertTrue(manifest["atlas_phenotype_independent"])
            self.assertIn("procedural", manifest["independence_claim"])
            atlas = json.loads(Path(manifest["outputs"]["interval_atlas_manifest"]).read_text(encoding="utf-8"))
            self.assertFalse(atlas["phenotype_used"])
            # every interval gene is still enumerated, including the unclassified one
            status = json.loads(Path(manifest["outputs"]["route_b0_status"]).read_text(encoding="utf-8"))
            self.assertEqual(3, status["interval_gene_count"])
            self.assertEqual(1, status["unclassified_genes"])

    def test_canonical_atlas_handoff_is_exposed_at_the_out_dir_root(self) -> None:
        """Issue #247: the atlas the next Route B stage consumes is one canonical
        machine-readable path at <out-dir>/interval_gene_atlas.tsv, not a nested
        layout a wrapper would have to infer."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "route_b"
            self.assertEqual(0, main([
                "route-b-prepare",
                "--candidates", str(_write(root, "candidates.tsv", CANDIDATES)),
                "--orthology", str(_write(root, "orthology.tsv", ORTHOLOGY)),
                "--reference-function", str(_write(root, "corpus.tsv", REFERENCE_CORPUS)),
                "--out-dir", str(out),
            ]))
            manifest = json.loads((out / "route_b_preparation_manifest.json").read_text(encoding="utf-8"))
            canonical = Path(manifest["outputs"]["interval_gene_atlas"])
            # the canonical handoff lives at the out-dir root, not inside atlas/.
            self.assertEqual(out / "interval_gene_atlas.tsv", canonical)
            self.assertTrue(canonical.exists())
            # byte-identical to the atlas inside the interval-atlas working dir.
            nested = out / "atlas" / "interval_gene_atlas.tsv"
            self.assertTrue(nested.exists())
            self.assertEqual(nested.read_bytes(), canonical.read_bytes())
            # the frozen fingerprint covers the canonical handoff file (full hash).
            fingerprint = manifest["frozen_output_fingerprints"]["interval_gene_atlas"]
            self.assertIsInstance(fingerprint, dict)
            self.assertEqual("full_sha256", fingerprint.get("method"))
            # the working dir and the canonical handoff are documented separately.
            self.assertEqual(str(out / "atlas"), manifest["outputs"]["interval_atlas_working_dir"])
            self.assertIn("canonical_atlas_rule", manifest)

    def test_reference_corpus_provenance_reaches_b0_lookup(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            corpus = _write(root, "rich_corpus.tsv", "\n".join([
                "reference_gene_id\tfunction_description\tmutant_phenotype\tpathway\tsource_id\tsource_url\t"
                "literature_source_ids\tphenotype_source_ids\tverified_source_ids\tevidence_codes\tevidence_type",
                "AtCLV1\treceptor kinase\tfasciated meristem\tmeristem maintenance\tUniProt:P123\thttps://example.test/P123\t"
                "PMID:1;DOI:10.1000/x\tPMID:1\tPMID:1;PMCID:PMC2\tEXP\texperimental",
            ]) + "\n")
            result = route_b_prepare(
                candidates_arg=str(_write(root, "candidates.tsv", CANDIDATES)),
                out_dir_arg=str(root / "route_b"),
                orthology_arg=str(_write(root, "orthology.tsv", ORTHOLOGY)),
                reference_function_args=[str(corpus)],
            )
            lookup = pd.read_csv(result["interval_reverse_lookup"], sep="\t", dtype=str).fillna("").set_index("gene_id")
            row = lookup.loc["Psat_G1"]
            self.assertEqual("https://example.test/P123", row["reference_source_url"])
            self.assertEqual({"PMID:1", "DOI:10.1000/x"}, set(row["literature_source_ids"].split(";")))
            self.assertEqual("PMID:1", row["phenotype_source_ids"])
            self.assertEqual({"PMID:1", "PMCID:PMC2"}, set(row["verified_source_ids"].split(";")))
            self.assertEqual("EXP", row["evidence_codes"])


class ReferenceFallbackTests(unittest.TestCase):
    def _candidate_file(self, root: Path, genes: list[str]) -> Path:
        rows = ["gene_id\tchrom\tstart\tend\tdescription"]
        rows.extend(f"{gene}\tchr1\t{index * 1000}\t{index * 1000 + 500}\t" for index, gene in enumerate(genes, 1))
        return _write(root, "fallback_candidates.tsv", "\n".join(rows) + "\n")

    def _family_file(self, root: Path, genes: list[str]) -> Path:
        rows = ["gene_id\tbest_subject\tbest_subject_symbol\tbest_subject_stem\tassignment_status"]
        rows.extend(f"{gene}\tREF_{gene}\tFAM{index}\tFAM{index}\tassigned" for index, gene in enumerate(genes, 1))
        return _write(root, "fallback_family.tsv", "\n".join(rows) + "\n")

    def test_complete_supplied_material_does_not_auto_search(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = self._candidate_file(root, ["G1", "G2"])
            family = self._family_file(root, ["G1", "G2"])
            proteins = _write(root, "fallback.fa", ">G1\nMASTPRQTVSEK\n>G2\nMASTPRQTVSEK\n")
            database = _write(root, "fallback.db", "placeholder")
            with patch("candidate_qtl.route_b_prepare._run_blastp", side_effect=AssertionError("search must be skipped")):
                result = route_b_prepare(str(candidates), str(root / "out"), proteins_arg=str(proteins),
                                         protein_db_arg=str(database), family_index_arg=str(family),
                                         auto_local_db=False)
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual(0, manifest["providers"]["homology_search"]["material_coverage"]["unresolved_gene_count"])

    def test_empty_supplied_orthology_falls_back_to_local_search(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = self._candidate_file(root, ["G1"])
            orthology = _write(root, "empty_orthology.tsv", "gene_id\treference_gene_id\torthogroup\n")
            proteins = _write(root, "fallback.fa", ">G1\nMASTPRQTVSEK\n")
            database = _write(root, "fallback.db", "placeholder")
            with patch("candidate_qtl.route_b_prepare._run_blastp", return_value=([], {"available": True})) as search:
                result = route_b_prepare(str(candidates), str(root / "out"), proteins_arg=str(proteins),
                                         protein_db_arg=str(database), orthology_arg=str(orthology),
                                         auto_local_db=False)
            search.assert_called_once()
            self.assertEqual(str(orthology), result["orthology"])

    def test_partial_material_search_preserves_supplied_rows_and_inventory(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = self._candidate_file(root, ["G1", "G2", "G3", "G4"])
            candidates.write_text(candidates.read_text().replace(
                "G1\tchr1\t1000\t1500\t", "G1\tchr1\t1000\t1500\tDNA mismatch repair activity"))
            family = self._family_file(root, ["G1"])
            orthology = _write(root, "partial_orthology.tsv", "gene_id\treference_gene_id\torthogroup\nG2\tREF_G2\tOG2\n")
            source_bytes = b">G1 annotated\nMASTPRQTVSEK\n>G2 mapped\nMASTPRQTVSEK\n"
            query_bytes = b">G3 unknown protein metadata\r\nMAST\r\nPRQTVSEK\r\n"
            source_bytes += query_bytes + b">OUTSIDE noninterval\nMASTPRQTVSEK\n"
            proteins = root / "fallback.fa"
            proteins.write_bytes(source_bytes)
            family_bytes, orthology_bytes = family.read_bytes(), orthology.read_bytes()
            database = _write(root, "fallback.db", "placeholder")
            rows = [{"gene_id": "G3", "reference_gene_id": "REF_G3", "orthogroup": "",
                     "relationship": "best_hit_homology", "provider": "local_blastp:fallback",
                     "percent_identity": 80.0, "query_coverage": 1.0, "evalue": 1e-20, "bitscore": 100.0}]
            with patch("candidate_qtl.route_b_prepare._run_blastp", return_value=(rows, {"available": True})) as search:
                result = route_b_prepare(str(candidates), str(root / "out"), proteins_arg=str(proteins),
                                         protein_db_arg=str(database), orthology_arg=str(orthology),
                                         family_index_arg=str(family), auto_local_db=False)
            merged = pd.read_csv(result["orthology"], sep="\t", dtype=str).fillna("")
            self.assertEqual({"G2", "G3"}, set(merged["gene_id"]))
            self.assertEqual("REF_G2", merged.set_index("gene_id").loc["G2", "reference_gene_id"])
            self.assertEqual("OG2", merged.set_index("gene_id").loc["G2", "orthogroup"])
            self.assertEqual(family_bytes, family.read_bytes())
            self.assertEqual(orthology_bytes, orthology.read_bytes())
            self.assertEqual(str(family), result["family_index"])
            self.assertEqual(source_bytes, proteins.read_bytes())
            query_path = search.call_args.args[0]
            self.assertNotEqual(proteins, query_path)
            self.assertEqual(query_bytes, query_path.read_bytes())
            manifest = json.loads(Path(result["manifest"]).read_text())
            subset = manifest["providers"]["homology_search"]["query_subset"]
            self.assertEqual(["G1", "G2", "G3", "OUTSIDE"], subset["original_gene_ids"])
            self.assertEqual(["G3"], subset["query_gene_ids"])
            self.assertEqual(["G3"], subset["completed_query_gene_ids"])
            self.assertEqual(1, subset["completed_query_gene_count"])
            self.assertEqual(["G1", "G2", "OUTSIDE"], subset["skipped_gene_ids"])
            self.assertEqual(["G4"], subset["missing_gene_ids"])
            self.assertEqual(["OUTSIDE"], subset["out_of_interval_gene_ids"])
            for label, count in (("original", 4), ("query", 1), ("skipped", 3), ("missing", 1), ("out_of_interval", 1)):
                self.assertEqual(count, subset[f"{label}_gene_count"])
            self.assertEqual(hashlib.sha256(source_bytes).hexdigest(), subset["original"]["sha256"])
            self.assertEqual(hashlib.sha256(query_bytes).hexdigest(), subset["query"]["sha256"])
            self.assertEqual(subset["original"], manifest["inputs"]["proteins"])
            self.assertEqual(str(query_path), manifest["outputs"]["reference_query_proteins"])
            inventory = pd.read_csv(Path(result["manifest"]).parent / "route_b_identity_inventory.tsv", sep="\t", dtype=str)
            self.assertEqual({"G1", "G2", "G3", "G4"}, set(inventory["gene_id"]))
            states = inventory.set_index("gene_id")["identity_state"].to_dict()
            self.assertEqual("identity_hit", states["G3"])
            for gene in ("G1", "G2", "G4"):
                self.assertEqual("not_evaluated_no_protein_query", states[gene])
            self.assertEqual(4, manifest["identity_coverage"]["identity_coverage"]["denominator_gene_count"])

    def test_fully_unresolved_fasta_still_searches_every_query(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = self._candidate_file(root, ["G1", "G2"])
            proteins = _write(root, "all.fa", ">G2 second\nMASTPRQTVSEK\n>G1 first\nMASTPRQTVSEK\n")
            database = _write(root, "local.db", "placeholder")
            with patch("candidate_qtl.route_b_prepare._run_blastp", return_value=([], {"available": True})) as search:
                result = route_b_prepare(str(candidates), str(root / "out"), proteins_arg=str(proteins),
                                         protein_db_arg=str(database), auto_local_db=False)
            search.assert_called_once()
            self.assertEqual(proteins.read_bytes(), search.call_args.args[0].read_bytes())
            manifest = json.loads(Path(result["manifest"]).read_text())
            subset = manifest["providers"]["homology_search"]["query_subset"]
            self.assertEqual(["G1", "G2"], subset["query_gene_ids"])
            self.assertEqual([], subset["skipped_gene_ids"])
            self.assertEqual([], subset["missing_gene_ids"])
            inventory = pd.read_csv(result["identity_inventory"], sep="\t")
            self.assertEqual({"evaluated_no_threshold_hit"}, set(inventory.identity_state))

    def test_empty_unresolved_query_skips_blast_and_database_resolution(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = self._candidate_file(root, ["G1", "G2"])
            family = self._family_file(root, ["G1"])
            proteins = _write(root, "empty_subset.fa", ">G1 already mapped\nMASTPRQTVSEK\n>G2 empty sequence\n\n")
            with patch("candidate_qtl.route_b_prepare._auto_local_protein_db", side_effect=AssertionError("no DB needed")), \
                 patch("candidate_qtl.route_b_prepare._run_blastp", side_effect=AssertionError("empty query")):
                result = route_b_prepare(str(candidates), str(root / "out"), proteins_arg=str(proteins),
                                         family_index_arg=str(family))
            manifest = json.loads(Path(result["manifest"]).read_text())
            provider = manifest["providers"]["homology_search"]
            self.assertEqual("not_evaluated_no_protein_query", provider["evaluation_state"])
            self.assertFalse(provider["available"])
            self.assertEqual(["G2"], provider["query_subset"]["missing_gene_ids"])
            self.assertEqual(0, provider["query_subset"]["query_gene_count"])
            self.assertEqual(b"", Path(manifest["outputs"]["reference_query_proteins"]).read_bytes())
            inventory = pd.read_csv(result["identity_inventory"], sep="\t")
            self.assertEqual(["G1", "G2"], inventory.gene_id.tolist())
            self.assertEqual({"not_evaluated_no_protein_query"}, set(inventory.identity_state))

    def test_query_subset_preserves_records_with_lf_crlf_or_cr_line_endings(self) -> None:
        for newline in ("\n", "\r\n", "\r"):
            with self.subTest(newline=newline), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                candidates = self._candidate_file(root, ["G1", "G2"])
                family = self._family_file(root, ["G1"])
                query_bytes = newline.join([">G2 unknown protein", "MAST", "PRQTVSEK", ""]).encode()
                source_bytes = query_bytes + newline.join([">G1 mapped", "MASTPRQTVSEK", ""]).encode()
                proteins = root / "line_endings.fa"
                proteins.write_bytes(source_bytes)
                database = _write(root, "local.db", "placeholder")
                with patch("candidate_qtl.route_b_prepare._run_blastp", return_value=([], {"available": True})) as search:
                    result = route_b_prepare(str(candidates), str(root / "out"), proteins_arg=str(proteins),
                                             protein_db_arg=str(database), family_index_arg=str(family), auto_local_db=False)
                search.assert_called_once()
                self.assertEqual(query_bytes, search.call_args.args[0].read_bytes())
                self.assertEqual(source_bytes, proteins.read_bytes())
                inventory = pd.read_csv(result["identity_inventory"], sep="\t").set_index("gene_id")
                self.assertEqual("not_evaluated_no_protein_query", inventory.loc["G1", "identity_state"])
                self.assertEqual("evaluated_no_threshold_hit", inventory.loc["G2", "identity_state"])

    def test_invalid_fasta_ids_fail_even_when_supplied_material_is_complete(self) -> None:
        for content in (">G1\nMAST\n>G1 duplicate\nMAST\n", "> \nMAST\n"):
            with self.subTest(content=content), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                candidates = self._candidate_file(root, ["G1"])
                family = self._family_file(root, ["G1"])
                proteins = _write(root, "bad.fa", content)
                with self.assertRaisesRegex(ValueError, "non-empty unique query IDs"):
                    route_b_prepare(str(candidates), str(root / "out"), proteins_arg=str(proteins),
                                    family_index_arg=str(family), auto_local_db=False)

    def test_invalid_candidate_ids_fail_before_search(self) -> None:
        for genes in (["G1", "G1"], ["G1", ""]):
            with self.subTest(genes=genes), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                candidates = self._candidate_file(root, genes)
                with self.assertRaisesRegex(ValueError, "non-empty unique interval gene IDs"):
                    route_b_prepare(str(candidates), str(root / "out"), auto_local_db=False)

    def test_missing_database_preserves_not_evaluated_inventory_and_query_snapshot(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = self._candidate_file(root, ["G1", "G2"])
            proteins = _write(root, "query.fa", ">G1\nMASTPRQTVSEK\n")
            with patch("candidate_qtl.route_b_prepare._run_blastp", side_effect=AssertionError("no provider")):
                result = route_b_prepare(str(candidates), str(root / "out"), proteins_arg=str(proteins),
                                         auto_local_db=False)
            manifest = json.loads(Path(result["manifest"]).read_text())
            provider = manifest["providers"]["homology_search"]
            self.assertEqual("not_evaluated_provider_unavailable", provider["evaluation_state"])
            self.assertEqual(["G1"], provider["query_subset"]["query_gene_ids"])
            self.assertEqual([], provider["query_subset"]["completed_query_gene_ids"])
            self.assertEqual(0, provider["query_subset"]["completed_query_gene_count"])
            self.assertEqual(["G2"], provider["query_subset"]["missing_gene_ids"])
            self.assertTrue(provider["missing_provider_is_not_negative_evidence"])
            inventory = pd.read_csv(result["identity_inventory"], sep="\t")
            self.assertEqual(["G1", "G2"], inventory.gene_id.tolist())
            self.assertEqual({"not_evaluated_provider_unavailable"}, set(inventory.identity_state))

    def test_unavailable_or_failed_blast_does_not_record_completed_queries(self) -> None:
        for binary, state in ((None, "not_evaluated_provider_unavailable"), ("/fake/blastp", "provider_failure")):
            with self.subTest(state=state), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                candidates = self._candidate_file(root, ["G1", "G2"])
                proteins = _write(root, "query.fa", ">G1\nMASTPRQTVSEK\n")
                database = _write(root, "local.db", "placeholder")
                failure = subprocess.CompletedProcess(["blastp"], 1, stdout="", stderr="provider test failure")
                with patch("candidate_qtl.route_b_prepare.shutil.which", return_value=binary), \
                     patch("candidate_qtl.route_b_prepare.run_tool", return_value=failure) as run:
                    result = route_b_prepare(str(candidates), str(root / "out"), proteins_arg=str(proteins),
                                             protein_db_arg=str(database), auto_local_db=False)
                if binary is None:
                    run.assert_not_called()
                else:
                    run.assert_called_once()
                manifest = json.loads(Path(result["manifest"]).read_text())
                provider = manifest["providers"]["homology_search"]
                self.assertEqual(state, provider["evaluation_state"])
                self.assertEqual(["G1"], provider["query_subset"]["query_gene_ids"])
                self.assertEqual([], provider["query_subset"]["completed_query_gene_ids"])
                self.assertEqual(0, provider["query_subset"]["completed_query_gene_count"])
                inventory = pd.read_csv(result["identity_inventory"], sep="\t")
                self.assertEqual({state}, set(inventory.identity_state))

    def test_inline_annotation_survives_no_search_prepare_to_semantic_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = _write(root, "annotated_candidates.tsv",
                                "gene_id\tchrom\tstart\tend\tdescription\tgo_terms\n"
                                "G1\tchr1\t100\t500\treceptor kinase meristem role\tGO:0004672\n")
            proteins = _write(root, "annotated.fa", ">G1\nMASTPRQTVSEK\n")
            database = _write(root, "annotated.db", "placeholder")
            out = root / "out"
            with patch("candidate_qtl.route_b_prepare._run_blastp", side_effect=AssertionError("search must be skipped")):
                prepared = route_b_prepare(str(candidates), str(out), proteins_arg=str(proteins),
                                           protein_db_arg=str(database), auto_local_db=False)
            atlas = pd.read_csv(prepared["interval_gene_atlas"], sep="\t", dtype=str).fillna("")
            self.assertEqual("receptor kinase meristem role", atlas.loc[0, "local_annotation_description"])
            bundle = build_route_b_semantic_bundle(prepared["interval_gene_atlas"], str(root / "bundle"),
                                                   phenotype="meristem phenotype")
            payload = json.loads(Path(bundle["bundle"]).read_text(encoding="utf-8"))
            gene = payload["units"][0]["genes"][0]
            self.assertEqual("receptor kinase meristem role", gene["local_annotation_description"])
            self.assertEqual("GO:0004672", gene["local_annotation_go_terms"])
            self.assertEqual("requested", payload["units"][0]["request_state"])
            self.assertTrue(payload["units"][0]["investigation_material"]["semantic_review_recommended"])
            identity = pd.read_csv(out / "route_b_identity_inventory.tsv", sep="\t")
            self.assertEqual("available", identity.loc[0, "functional_annotation_state"])

    def test_unknown_candidate_label_does_not_hide_supplied_function(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = _write(root, "c.tsv", "gene_id\tchrom\tstart\tend\tdescription\n"
                                "G1\tchr1\t100\t500\thypothetical protein\n")
            annotations = _write(root, "annotations.tsv", "gene_id\tfunction_description\n"
                                 "G1\tDNA mismatch repair activity\n")
            with patch("candidate_qtl.route_b_prepare._auto_local_protein_db",
                       side_effect=AssertionError("DB resolution must be skipped")):
                prepared = route_b_prepare(str(candidates), str(root / "out"),
                                           annotations_arg=str(annotations))
            bundle = build_route_b_semantic_bundle(prepared["interval_gene_atlas"], str(root / "bundle"),
                                                   phenotype="synthetic trait")
            unit = json.loads(Path(bundle["bundle"]).read_text())["units"][0]
            self.assertEqual("hypothetical protein\nDNA mismatch repair activity",
                             unit["genes"][0]["local_annotation_description"])
            self.assertEqual("requested", unit["request_state"])
            self.assertEqual(0, unit["source_backed_function_record_count"])

    def test_placeholder_or_off_interval_mapping_does_not_suppress_fallback(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = _write(root, "c.tsv", "gene_id\tchrom\tstart\tend\tdescription\n"
                                "G1\tchr1\t100\t500\thypothetical protein\n")
            orthology = _write(root, "other.tsv", "gene_id\treference_gene_id\nOTHER\tREF_OTHER\n")
            proteins = _write(root, "query.fa", ">G1\nMASTPRQTVSEK\n")
            with patch("candidate_qtl.route_b_prepare._auto_local_protein_db", return_value=root / "db"), \
                 patch("candidate_qtl.route_b_prepare._run_blastp",
                       return_value=([], {"available": True})) as search:
                prepared = route_b_prepare(str(candidates), str(root / "out"), proteins_arg=str(proteins),
                                           orthology_arg=str(orthology))
            search.assert_called_once()
            identity = pd.read_csv(Path(prepared["manifest"]).parent / "route_b_identity_inventory.tsv", sep="\t")
            self.assertEqual(["G1"], identity.gene_id.tolist())
            self.assertEqual("not_available", identity.loc[0, "functional_annotation_state"])


if __name__ == "__main__":
    unittest.main()
