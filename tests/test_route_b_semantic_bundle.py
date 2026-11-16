from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.route_b_semantic_bundle import build_route_b_semantic_bundle


class RouteBSemanticBundleTests(unittest.TestCase):
    def _atlas(self, path: Path) -> None:
        pd.DataFrame([
            {
                "gene_id": "G1", "chrom": "chr1", "start": "100", "end": "200",
                "family_id": "legacy-clv", "evolutionary_family_id": "EV_CLV",
                "evolutionary_family_status": "assigned", "orthogroup": "OG_CLV",
                "orthology_relationship": "one_to_one", "orthology_reference_gene": "CLV1",
                "subfamily_id": "CLV1-like", "domain_ids": "PF00069",
                "domain_architecture": "kinase", "broad_function_bins": "signaling",
                "primary_broad_function": "signaling", "conservation_class": "conserved_cross_species",
                "evidence_completeness": "family;orthology;domain", "annotation_rescue_needed": "no",
                "gene_model_completeness": "full_model", "classification_basis": "family+orthology", "notes": "",
            },
            {
                "gene_id": "G2", "chrom": "chr1", "start": "300", "end": "400",
                "family_id": "legacy-clv", "evolutionary_family_id": "EV_CLV",
                "evolutionary_family_status": "assigned", "orthogroup": "OG_CLV",
                "orthology_relationship": "one_to_many", "orthology_reference_gene": "CLV1",
                "subfamily_id": "CLV1-like", "domain_ids": "PF00069",
                "domain_architecture": "kinase", "broad_function_bins": "signaling",
                "primary_broad_function": "signaling", "conservation_class": "conserved_cross_species",
                "evidence_completeness": "family;orthology;domain", "annotation_rescue_needed": "no",
                "gene_model_completeness": "full_model", "classification_basis": "family+orthology", "notes": "",
            },
            {
                "gene_id": "DARK1", "chrom": "chr1", "start": "500", "end": "600",
                "family_id": "", "evolutionary_family_id": "", "evolutionary_family_status": "unresolved",
                "orthogroup": "", "orthology_relationship": "", "orthology_reference_gene": "",
                "subfamily_id": "", "domain_ids": "", "domain_architecture": "",
                "broad_function_bins": "", "primary_broad_function": "",
                "conservation_class": "unresolved", "evidence_completeness": "none",
                "annotation_rescue_needed": "yes", "gene_model_completeness": "partial_model",
                "classification_basis": "unclassified", "notes": "retain explicit dark gene",
            },
        ]).to_csv(path, sep="\t", index=False, lineterminator="\n")

    def _literature(self, literature: Path, status: Path) -> None:
        pd.DataFrame([
            {
                "interval_gene_ids": "G1;G2", "reference_gene": "CLV1",
                "reference_function": "receptor kinase controlling meristem size",
                "query": "CLV1 fasciation", "query_status": "retrieved", "article_rank": "1",
                "source_id": "PMID:123", "title": "CLV1 controls shoot meristem size",
                "journal": "Plant Journal", "year": "2020", "doi": "10.1000/example",
                "first_author": "Example", "matched_phenotype_terms": "fasciation",
                "claim_type": "literature_link_candidate", "provider": "fixture",
                "retrieved_at": "2026-08-16T00:00:00+00:00", "rank_contribution": "none",
            }
        ]).to_csv(literature, sep="\t", index=False, lineterminator="\n")
        pd.DataFrame([
            {"reference_gene": "CLV1", "interval_gene_count": "2", "query_status": "retrieved", "article_count": "1", "note": ""}
        ]).to_csv(status, sep="\t", index=False, lineterminator="\n")

    def test_bundle_compresses_family_but_never_drops_dark_gene(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atlas = root / "atlas.tsv"
            literature = root / "literature.tsv"
            status = root / "status.tsv"
            self._atlas(atlas)
            self._literature(literature, status)

            result = build_route_b_semantic_bundle(
                str(atlas), str(root / "out"), phenotype="fasciated shoot apical meristem",
                literature_arg=str(literature), literature_status_arg=str(status), case_id="b403",
            )
            self.assertEqual(3, result["interval_gene_count"])
            self.assertEqual(2, result["semantic_unit_count"])
            # The dark gene keeps its unit and its place in the universe, but
            # carries nothing to bridge from, so it is not requested.
            self.assertEqual(1, result["prompt_count"])
            self.assertEqual(1, result["skipped_unit_count"])

            units = pd.read_csv(result["units"], sep="\t", dtype=str).fillna("")
            family = units[units["unit_anchor"].eq("EV_CLV")].iloc[0]
            self.assertEqual("2", family["member_count"])
            self.assertEqual("G1;G2", family["member_interval_genes"])
            self.assertEqual("1", family["literature_link_count"])
            dark = units[units["unit_anchor"].eq("DARK1")].iloc[0]
            self.assertEqual("interval_gene", dark["unit_type"])
            self.assertEqual("1", dark["unclassified_member_count"])

            bundle = json.loads(Path(result["bundle"]).read_text(encoding="utf-8"))
            observed = {gene for unit in bundle["units"] for gene in unit["member_interval_genes"]}
            self.assertEqual({"G1", "G2", "DARK1"}, observed)
            self.assertTrue(bundle["candidate_identity_preserved"])

    def test_all_identity_hit_alternatives_reach_semantic_prompt_as_handles_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atlas = root / "atlas.tsv"
            identity_hits = root / "route_b_identity_hits.tsv"
            self._atlas(atlas)
            pd.DataFrame([
                {"focal_gene_id": "G1", "hit_accession": "A1", "hit_species": "Medicago truncatula",
                 "hit_gene_symbol": "LYK2", "hit_protein_name": "kinase", "hit_description": "header A1",
                 "hit_metadata_state": "partial", "hit_metadata_records_json": "[]",
                 "raw_reference_gene_id": "sp|A1|LYK2_MEDTR", "percent_identity": "81.0",
                 "query_coverage": "0.99", "evalue": "1e-20", "bitscore": "100",
                 "provider": "local_blastp:test", "database": "legume-test", "identity_policy": "test_v1"},
                {"focal_gene_id": "G1", "hit_accession": "A2", "hit_species": "Medicago truncatula",
                 "hit_gene_symbol": "LYK3", "hit_protein_name": "LysM receptor kinase", "hit_description": "header A2",
                 "hit_metadata_state": "partial", "hit_metadata_records_json": "[]",
                 "raw_reference_gene_id": "sp|A2|LYK3_MEDTR", "percent_identity": "80.0",
                 "query_coverage": "0.99", "evalue": "1e-19", "bitscore": "98",
                 "provider": "local_blastp:test", "database": "legume-test", "identity_policy": "test_v1"},
            ]).to_csv(identity_hits, sep="\t", index=False, lineterminator="\n")

            result = build_route_b_semantic_bundle(
                str(atlas), str(root / "out"), phenotype="short stem",
                identity_hits_arg=str(identity_hits), case_id="alternatives",
            )
            bundle = json.loads(Path(result["bundle"]).read_text(encoding="utf-8"))
            prompts = [json.loads(line) for line in Path(result["prompts"]).read_text(encoding="utf-8").splitlines()]
            unit = next(unit for unit in bundle["units"] if "G1" in unit["member_interval_genes"])
            g1 = next(gene for gene in unit["genes"] if gene["gene_id"] == "G1")
            hits = g1["sequence_similarity_hits"]
            self.assertEqual({"A1", "A2"}, {hit["hit_accession"] for hit in hits})
            self.assertTrue(all(hit["interpretation"] == "sequence_similarity_handle_only" for hit in hits))
            self.assertTrue(all(hit["rank_contribution"] == "none" for hit in hits))
            self.assertTrue(all("hit_metadata_records_json" in hit for hit in hits))
            prompt = next(row["user"] for row in prompts if row["route_b_unit_id"] == unit["unit_id"])
            self.assertIn("LYK2", prompt)
            self.assertIn("LYK3", prompt)
            self.assertNotIn("hit_metadata_records_json", prompt)
            self.assertEqual([], unit["stable_source_ids"])
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertIsNotNone(manifest["inputs"]["route_b_identity_hits"])

    def test_prompts_are_unit_level_proposals_not_scores(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atlas = root / "atlas.tsv"
            self._atlas(atlas)
            result = build_route_b_semantic_bundle(
                str(atlas), str(root / "out"), phenotype="fasciated meristem", case_id="b403",
            )
            lines = [json.loads(line) for line in Path(result["prompts"]).read_text(encoding="utf-8").splitlines()]
            # One prompt for the family unit; the dark gene is not requested.
            self.assertEqual(1, len(lines))
            self.assertTrue(all(row["gene_id"] == "__route_b_semantic_unit__" for row in lines))
            self.assertTrue(all("Return JSON only" in row["system"] for row in lines))
            self.assertTrue(all("do not rank genes numerically" in row["system"] for row in lines))
            self.assertTrue(all("Route: interval_first; B0 inventory already frozen" in row["user"] for row in lines))

    def test_exported_prompt_separates_source_observation_from_target_inference(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atlas = root / "atlas.tsv"
            self._atlas(atlas)
            result = build_route_b_semantic_bundle(
                str(atlas), str(root / "out"), phenotype="root and shoot development", case_id="atomic",
            )
            prompts = [json.loads(line) for line in Path(result["prompts"]).read_text().splitlines()]
            self.assertEqual(1, len(prompts))
            instruction = prompts[0]["system"]
            for required in (
                "mechanistic_review_seed",
                "specific phenotype-directed mechanism worth preserving for later review",
                "adopted as a first-stage candidate",
                "If no specific bridge can be stated from the supplied unit, return an empty claims list",
                "one independently source-supported observation or process",
                "different experimental genes, organisms, or phenotype components",
                "Each source_ids entry must support that claim's stated observation",
                "conditional relevance to the focal phenotype",
                "uncertainty_note or alternative_hypothesis",
                "need not explain every component",
                "homology alone does not establish a functional experiment",
                "reference_gene must exactly match a supplied unit reference endpoint",
                "not a ban on upstream Web research",
                "from the initial assessment and throughout reasoning",
                "do not impose a search-free preliminary gate, one-paper limit or uniform per-gene search quota",
                "Supplied-material-only benchmark assignments are a separate explicit execution mode",
                "new frozen bundle before submitting that claim",
                "Model knowledge and search snippets are leads, not verified paper evidence",
                "offline verifier consumes the supplied source cache",
                "reported reference observation -> homolog/family correspondence -> proposed focal mechanism",
                "without an experiment in the focal species",
                "First-stage admission and verified trait support are separate outcomes",
                "Start with the supplied interval genes/copies and their reference endpoints",
                "Do not presume that every unit must produce a candidate claim",
                "assessed neutrally",
                "generic kinase or transcription-factor involvement is not a specific bridge",
                "A multi-step but specifically grounded cross-species mechanism remains eligible for review",
                "A reviewed unit may legitimately return an empty claims list",
                "not biological evidence against the gene",
                "primary task is a cross-species biological bridge",
                "Same-species literature is supplementary evidence",
                "without species-specific or example-specific selection rules",
            ):
                self.assertIn(required, instruction)
            self.assertLess(
                instruction.index("Start with the supplied interval genes/copies"),
                instruction.index("Your primary task is a cross-species biological bridge"),
            )
            # Guidance must not change family/copy grouping or dark-gene retention.
            self.assertEqual(3, result["interval_gene_count"])
            self.assertEqual(2, result["semantic_unit_count"])
            self.assertEqual(1, result["skipped_unit_count"])

    def test_source_backed_b0_function_reaches_frozen_semantic_prompt(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atlas = root / "atlas.tsv"
            lookup = root / "lookup.tsv"
            self._atlas(atlas)
            pd.DataFrame([{
                "gene_id": "G1", "reference_gene": "CLV1",
                "reference_function": "receptor kinase controlling meristem size",
                "reference_pathway": "stem cell homeostasis", "reference_pathway_id": "GO:0019827",
                "reference_mutant_phenotype": "fasciated meristem",
                "reference_source_id": "PMID:123", "lookup_status": "source_backed",
            }]).to_csv(lookup, sep="\t", index=False, lineterminator="\n")

            result = build_route_b_semantic_bundle(
                str(atlas), str(root / "out"), phenotype="fasciated meristem",
                lookup_arg=str(lookup), case_id="b403",
            )
            bundle = json.loads(Path(result["bundle"]).read_text(encoding="utf-8"))
            family = next(unit for unit in bundle["units"] if unit["unit_anchor"] == "EV_CLV")
            self.assertEqual(["PMID:123"], family["stable_source_ids"])
            self.assertEqual(1, family["source_backed_function_record_count"])
            g1 = next(row for row in family["genes"] if row["gene_id"] == "G1")
            self.assertEqual("receptor kinase controlling meristem size", g1["reference_function"])
            self.assertEqual("fasciated meristem", g1["reference_mutant_phenotype"])
            prompt = next(
                json.loads(line) for line in Path(result["prompts"]).read_text(encoding="utf-8").splitlines()
                if "receptor kinase controlling meristem size" in line
            )
            self.assertIn("PMID:123", prompt["user"])
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual(str(lookup), manifest["inputs"]["interval_reverse_lookup"]["path"])

    def test_source_roster_deduplicates_literature_axes_and_keeps_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atlas = root / "atlas.tsv"
            lookup = root / "lookup.tsv"
            self._atlas(atlas)
            pd.DataFrame([{
                "gene_id": "G1", "reference_gene": "CLV1", "reference_function": "receptor kinase",
                "reference_source_id": "UniProt:P12345", "reference_source_url": "https://example.test/P12345",
                "literature_source_ids": "PMID:123;DOI:10.1000/x", "phenotype_source_ids": "PMID:123",
                "verified_source_ids": "PMID:123;PMCID:PMC999", "evidence_codes": "EXP",
                "evidence_type": "experimental", "provenance_tier": "retrieved", "source_database": "UniProtKB",
            }]).to_csv(lookup, sep="\t", index=False, lineterminator="\n")
            result = build_route_b_semantic_bundle(
                str(atlas), str(root / "out"), phenotype="fasciated meristem", lookup_arg=str(lookup)
            )
            bundle = json.loads(Path(result["bundle"]).read_text(encoding="utf-8"))
            family = next(unit for unit in bundle["units"] if unit["unit_anchor"] == "EV_CLV")
            self.assertEqual(["DOI:10.1000/x", "PMCID:PMC999", "PMID:123", "UniProt:P12345"], family["stable_source_ids"])
            self.assertEqual("https://example.test/P12345", family["source_provenance"][0]["source_url"])
            self.assertEqual("EXP", family["source_provenance"][0]["evidence_codes"])

    def test_duplicate_lookup_gene_ids_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atlas = root / "atlas.tsv"
            lookup = root / "lookup.tsv"
            self._atlas(atlas)
            pd.DataFrame([{"gene_id": "G1"}, {"gene_id": "G1"}]).to_csv(
                lookup, sep="\t", index=False, lineterminator="\n"
            )
            with self.assertRaisesRegex(ValueError, "unique"):
                build_route_b_semantic_bundle(
                    str(atlas), str(root / "out"), phenotype="trait", lookup_arg=str(lookup)
                )

    def test_duplicate_atlas_gene_ids_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atlas = root / "atlas.tsv"
            pd.DataFrame([{"gene_id": "G1"}, {"gene_id": "G1"}]).to_csv(
                atlas, sep="\t", index=False, lineterminator="\n"
            )
            with self.assertRaisesRegex(ValueError, "unique"):
                build_route_b_semantic_bundle(str(atlas), str(root / "out"), phenotype="trait")

    def test_b0_coverage_manifest_sum_matches_total_genes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atlas = root / "atlas.tsv"
            pd.DataFrame([
                {"gene_id": "G1", "evolutionary_family_id": "EV1", "evolutionary_family_status": "assigned"},
                {"gene_id": "G2", "evolutionary_family_id": "", "evolutionary_family_status": "unassigned"},
                {"gene_id": "G3", "evolutionary_family_id": "", "evolutionary_family_status": "unresolved"},
                {"gene_id": "G4", "evolutionary_family_id": "", "evolutionary_family_status": "failed"},
            ]).to_csv(atlas, sep="\t", index=False, lineterminator="\n")

            result = build_route_b_semantic_bundle(
                str(atlas), str(root / "out"), phenotype="flowering time", case_id="b0_cov"
            )
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            b0_cov = manifest["b0_coverage"]
            self.assertTrue(b0_cov["sum_matches_total"])
            self.assertEqual(4, b0_cov["total"])
            self.assertEqual(1, b0_cov["counts"]["assigned"])
            self.assertEqual(1, b0_cov["counts"]["unassigned"])
            self.assertEqual(1, b0_cov["counts"]["unresolved"])
            self.assertEqual(1, b0_cov["counts"]["failed"])
            self.assertEqual(
                4,
                b0_cov["counts"]["assigned"]
                + b0_cov["counts"]["unassigned"]
                + b0_cov["counts"]["unresolved"]
                + b0_cov["counts"]["failed"],
            )

    def test_unit_table_shows_curated_function_count_beside_literature_count(self) -> None:
        """A canonical run has no B1 literature links; the curated count must still show.

        `literature_link_count` is filled only by the offline-by-default legacy
        B1 reverse-literature lane, so a unit that carries a real curated
        function reads as `0` there.  The g7 bench read that column alone and
        concluded every unit was empty, so the units table publishes both.
        """
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atlas = root / "atlas.tsv"
            lookup = root / "lookup.tsv"
            self._atlas(atlas)
            pd.DataFrame([{
                "gene_id": "G1", "reference_gene": "CLV1",
                "reference_function": "receptor kinase controlling meristem size",
                "reference_source_id": "PMID:123", "lookup_status": "source_backed",
            }]).to_csv(lookup, sep="\t", index=False, lineterminator="\n")

            result = build_route_b_semantic_bundle(
                str(atlas), str(root / "out"), phenotype="fasciated meristem",
                lookup_arg=str(lookup), case_id="b403",
            )
            units = pd.read_csv(result["units"], sep="\t", dtype=str).fillna("")
            family = units[units["unit_anchor"].eq("EV_CLV")].iloc[0]
            self.assertEqual("0", family["literature_link_count"])
            self.assertEqual("1", family["source_backed_function_record_count"])

    def test_identity_signal_counts_expose_families_hidden_behind_unresolved_status(self) -> None:
        """`assigned: 0` must not be readable as "no family was assigned".

        `_classify_gene_b0_coverage` short-circuits on an unresolved
        evolutionary status before it reads `family_id`, which is exactly the
        g7 atlas shape: every gene unresolved, most genes carrying a family.
        """
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atlas = root / "atlas.tsv"
            pd.DataFrame([
                {"gene_id": "G1", "family_id": "FB121", "evolutionary_family_id": "",
                 "evolutionary_family_status": "unresolved", "orthogroup": "",
                 "domain_architecture": "", "orthology_reference_gene": "A0A072UJJ9"},
                {"gene_id": "G2", "family_id": "FB121", "evolutionary_family_id": "",
                 "evolutionary_family_status": "unresolved", "orthogroup": "",
                 "domain_architecture": "", "orthology_reference_gene": ""},
                {"gene_id": "G3", "family_id": "", "evolutionary_family_id": "",
                 "evolutionary_family_status": "unresolved", "orthogroup": "",
                 "domain_architecture": "", "orthology_reference_gene": ""},
            ]).to_csv(atlas, sep="\t", index=False, lineterminator="\n")

            result = build_route_b_semantic_bundle(
                str(atlas), str(root / "out"), phenotype="stem length", case_id="b0_signal"
            )
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            b0_cov = manifest["b0_coverage"]
            self.assertEqual(0, b0_cov["counts"]["assigned"])
            self.assertEqual(3, b0_cov["counts"]["unresolved"])
            signals = b0_cov["identity_signal_counts"]
            self.assertEqual(2, signals["family_id"])
            self.assertEqual(1, signals["orthology_reference_gene"])
            self.assertEqual(0, signals["evolutionary_family_id"])
            self.assertEqual(0, signals["orthogroup"])
            self.assertEqual(0, signals["domain_architecture"])
            self.assertIn("assigned=0 does not mean no family was assigned", b0_cov["reading_note"])

            bundle = json.loads(Path(result["bundle"]).read_text(encoding="utf-8"))
            self.assertEqual(signals, bundle["b0_coverage"]["identity_signal_counts"])

    def test_material_class_controls_requests_without_losing_sparse_or_noisy_units(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atlas, lookup, literature, hits = (root / name for name in ("atlas.tsv", "lookup.tsv", "links.tsv", "hits.tsv"))
            truncated = "ATP-dependent transport [truncated]"
            rows = [
                {"gene_id": "G1", "family_id": "F001", "family_label": "F001", "orthogroup": "OG001", "orthology_reference_gene": "P12345"},
                {"gene_id": "G2", "family_id": "F001", "orthogroup": "OG001", "orthology_reference_gene": "P67890"},
                {"gene_id": "G3"},
                {"gene_id": "G4", "domain_ids": "PF00069"},
                {"gene_id": "G5", "family_id": "F005", "family_label": "transferase family"},
                {"gene_id": "G6", "local_annotation_description": truncated},
                {"gene_id": "G7", "local_annotation_description": "unknown;uncharacterized protein;NaN;[truncated]"},
                {"gene_id": "G8", "orthology_reference_gene": "P88888"},
                {"gene_id": "G9"},
                {"gene_id": "G10", "reference_function": "measured enzyme activity", "reference_gene": "P10101"},
            ]
            pd.DataFrame(rows).fillna("").to_csv(atlas, sep="\t", index=False)
            pd.DataFrame([
                {"gene_id": "G1", "reference_gene": "P12345", "reference_source_id": "DB:P12345", "reference_source_url": "https://example.test/P12345"},
                {"gene_id": "G8", "reference_gene": "P88888", "reference_function": "measured transport activity", "reference_source_id": "PMID:888", "evidence_codes": "EXP"},
            ]).fillna("").to_csv(lookup, sep="\t", index=False)
            pd.DataFrame([{"interval_gene_ids": "G3", "source_id": "PMID:333", "source_url": "https://example.test/333", "title": "Unretrieved article title"}]).to_csv(literature, sep="\t", index=False)
            pd.DataFrame([{"focal_gene_id": "G9", "hit_accession": "P99999", "hit_protein_name": "transporter", "raw_reference_gene_id": "sp|P99999|TRANSPORT", "hit_metadata_state": "partial", "hit_metadata_records_json": "[]"}]).to_csv(hits, sep="\t", index=False)
            schedules = []
            for index, phenotype in enumerate(("altered development", "different neutral phenotype")):
                result = build_route_b_semantic_bundle(str(atlas), str(root / f"out-{index}"), phenotype=phenotype,
                                                      lookup_arg=str(lookup), literature_arg=str(literature), identity_hits_arg=str(hits))
                bundle = json.loads(Path(result["bundle"]).read_text(encoding="utf-8"))
                by_gene = {gene: unit for unit in bundle["units"] for gene in unit["member_interval_genes"]}
                self.assertEqual({row["gene_id"] for row in rows}, set(by_gene))
                self.assertEqual(["G1", "G2"], by_gene["G1"]["member_interval_genes"])
                self.assertEqual(["F001"], by_gene["G1"]["identity"]["family_ids"])
                self.assertEqual(["OG001"], by_gene["G1"]["identity"]["orthogroup_ids"])
                # #10: a hit's original protein annotation is delivered as unverified
                # initial material, so G9 is now requested (legacy bundles keep the
                # handles-only reading; see test_route_b_hit_initial_material).
                self.assertEqual(6, result["prompt_count"])
                schedules.append({gene: unit["request_state"] for gene, unit in by_gene.items()})
                for gene in ("G1", "G2", "G3", "G7"):
                    unit = by_gene[gene]
                    self.assertEqual("not_requested_no_reviewable_material", unit["request_state"])
                    self.assertEqual("not_evaluated", unit["review_status"])
                    self.assertEqual("no_reviewable_material", unit["not_evaluated_reason"])
                    self.assertFalse(unit["investigation_material"]["semantic_review_recommended"])
                    self.assertEqual("needs_material_retrieval", unit["investigation_material"]["state"])
                self.assertEqual("annotation_context", by_gene["G9"]["source_material_class"])
                self.assertEqual("retrieval_handles_only", by_gene["G3"]["source_material_class"])
                self.assertEqual("no_reviewable_material", by_gene["G7"]["source_material_class"])
                self.assertEqual(["DB:P12345"], by_gene["G1"]["stable_source_ids"])
                self.assertEqual("https://example.test/P12345", by_gene["G1"]["source_provenance"][0]["source_url"])
                self.assertEqual("PMID:333", by_gene["G3"]["literature_link_candidates"][0]["source_id"])
                self.assertEqual("sp|P99999|TRANSPORT", by_gene["G9"]["genes"][0]["sequence_similarity_hits"][0]["raw_reference_gene_id"])
                self.assertEqual(truncated, by_gene["G6"]["genes"][0]["local_annotation_description"])
                for gene in ("G4", "G5", "G6", "G8", "G9", "G10"):
                    self.assertEqual("requested", by_gene[gene]["request_state"])
                    self.assertEqual("not_evaluated", by_gene[gene]["review_status"])
                    self.assertEqual("none", by_gene[gene]["rank_contribution"])
                self.assertEqual("annotation_context", by_gene["G4"]["source_material_class"])
                self.assertEqual("function_context", by_gene["G8"]["source_material_class"])
                self.assertEqual("measured enzyme activity", by_gene["G10"]["genes"][0]["reference_function"])
                self.assertEqual(1, by_gene["G8"]["source_backed_function_record_count"])
                self.assertFalse(any("verified_support" in unit for unit in bundle["units"]))
                skipped = pd.read_csv(result["units_not_requested"], sep="\t", dtype=str)
                self.assertTrue(skipped.review_status.eq("not_evaluated").all())
                self.assertTrue(skipped.not_evaluated_reason.eq("no_reviewable_material").all())
            self.assertEqual(schedules[0], schedules[1])

    def test_review_cost_is_deterministic_input_estimate_before_any_model_call(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atlas = root / "atlas.tsv"
            self._atlas(atlas)
            results = [build_route_b_semantic_bundle(str(atlas), str(root / f"out-{i}"), phenotype="neutral phenotype", case_id="cost") for i in range(2)]
            bundles = [json.loads(Path(result["bundle"]).read_text(encoding="utf-8")) for result in results]
            self.assertEqual(bundles[0]["units"], bundles[1]["units"])
            self.assertEqual(bundles[0]["execution_prompt_fingerprints"], bundles[1]["execution_prompt_fingerprints"])
            by_unit = {unit["unit_id"]: unit for unit in bundles[0]["units"]}
            prompts = [json.loads(line) for line in Path(results[0]["prompts"]).read_text(encoding="utf-8").splitlines()]
            for prompt in prompts:
                unit = by_unit[prompt["route_b_unit_id"]]
                estimate = unit["estimated_review_cost"]
                prefix, body = prompt["user"].split("\n\n", 1)
                payload = json.loads(body)
                payload.pop("estimated_review_cost")
                text = prefix + "\n\n" + json.dumps(payload, ensure_ascii=False, sort_keys=True)
                size = len(prompt["system"].encode("utf-8")) + len(text.encode("utf-8"))
                self.assertEqual(size, estimate["input_utf8_bytes"])
                self.assertEqual((size + 3) // 4, estimate["estimated_input_tokens"])
                self.assertEqual("estimate_not_actual_usage", estimate["status"])
                self.assertIsNone(estimate["actual_input_tokens"])
                self.assertEqual(1, estimate["estimated_model_call_count"])
            total = bundles[0]["review_cost_estimate"]
            self.assertEqual(len(prompts), total["estimated_model_call_count"])
            self.assertEqual(sum(unit["estimated_review_cost"]["estimated_input_tokens"] for unit in by_unit.values() if unit["request_state"] == "requested"), total["estimated_input_tokens"])
            self.assertIsNone(total["actual_input_tokens"])
            skipped = next(unit for unit in by_unit.values() if unit["request_state"] != "requested")
            self.assertEqual(0, skipped["estimated_review_cost"]["estimated_model_call_count"])
            manifest = json.loads(Path(results[0]["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual(total, manifest["review_cost_estimate"])


if __name__ == "__main__":
    unittest.main()
