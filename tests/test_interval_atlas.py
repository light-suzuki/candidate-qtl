"""Phenotype-independent interval atlas (issue #168)."""

from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl import interval_atlas as atlas_module
from candidate_qtl.cli import main


class BinMapPackagingTests(unittest.TestCase):
    """The bin map ships in the wheel and stays in sync, like the weight registry."""

    def test_default_bin_map_resolves_and_loads(self) -> None:
        path = atlas_module.default_bin_map_path()
        self.assertTrue(path.exists(), f"default bin map not found at {path}")
        bins, meta = atlas_module.load_bin_map()
        self.assertGreater(len(bins), 0)
        self.assertGreater(meta["bin_count"], 5)
        self.assertIn("sha256", meta)

    def test_vendored_copy_is_byte_identical_to_source_tree(self) -> None:
        repo_root = Path(__file__).resolve().parents[1]
        source_tree = repo_root / "templates" / "broad_function_bins.tsv"
        vendored = repo_root / "src" / "candidate_qtl" / "resources" / "broad_function_bins.tsv"
        self.assertTrue(vendored.exists(), "vendored bin map is missing from the package")
        self.assertEqual(source_tree.read_bytes(), vendored.read_bytes(),
                         "templates/broad_function_bins.tsv and the vendored copy have drifted; re-sync them")

    def test_unknown_match_type_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bins.tsv"
            path.write_text("bin_id\tbin_label\tmatch_type\tpattern\nx\tX\tastrology\tfoo\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                atlas_module.load_bin_map(str(path))


class IntervalAtlasTests(unittest.TestCase):
    def _fixture(self, root: Path) -> dict[str, Path]:
        paths = {name: root / name for name in (
            "candidates.tsv", "family_index.tsv", "domains.tsv", "orthology.tsv", "annotations.tsv", "flags.tsv")}
        paths["candidates.tsv"].write_text("\n".join([
            "gene_id\tchrom\tstart\tend\tprotein_length",
            "G_RLK1\tchr1\t1000\t3000\t980",
            "G_RLK2\tchr1\t5000\t7000\t960",
            "G_RLK3\tchr1\t400000\t402000\t940",
            "G_TF1\tchr1\t20000\t22000\t420",
            "G_TRANS1\tchr1\t30000\t31000\t510",
            "G_ORPHAN\tchr1\t40000\t40300\t85",
            "G_CONFLICT\tchr1\t50000\t52000\t600",
        ]) + "\n", encoding="utf-8")
        paths["family_index.tsv"].write_text("\n".join([
            "gene_id\tchrom\tstart\tend\tbest_subject\tbest_subject_symbol\tbest_subject_stem\tfamily_stems\tassignment_status",
            "G_RLK1\tchr1\t1000\t3000\tsp|P1|CLV1_ARATH\tCLV1\tCLV\tCLV\tassigned",
            "G_RLK2\tchr1\t5000\t7000\tsp|P2|CLV2_ARATH\tCLV2\tCLV\tCLV\tassigned",
            "G_RLK3\tchr1\t400000\t402000\tsp|P3|CLV3_ARATH\tCLV3\tCLV\tCLV\tassigned",
            "G_TF1\tchr1\t20000\t22000\tsp|P4|WOX4_ARATH\tWOX4\tWOX\tWOX\tassigned",
            "G_TRANS1\tchr1\t30000\t31000\tsp|P5|PIN1_ARATH\tPIN1\tPIN\tPIN\tassigned",
            "G_ORPHAN\t\t\t\t\t\t\t\tdark_no_qualifying_hit",
            "G_CONFLICT\tchr1\t50000\t52000\tsp|P6|BAM1_ARATH\tBAM1\tBAM\tBAM\tassigned",
        ]) + "\n", encoding="utf-8")
        paths["domains.tsv"].write_text("\n".join([
            "gene_id\tdomain_id\tdomain_name\tstart\tsource",
            "G_RLK1\tPF00560\tLRR\t50\tPfam",
            "G_RLK1\tPF00069\tPkinase\t600\tPfam",
            "G_RLK2\tPF00560\tLRR\t60\tPfam",
            "G_RLK2\tPF00069\tPkinase\t610\tPfam",
            "G_TF1\tPF00046\tHomeobox\t80\tPfam",
        ]) + "\n", encoding="utf-8")
        paths["orthology.tsv"].write_text("\n".join([
            "gene_id\torthogroup\trelationship\treference_gene_id\tprovider",
            "G_RLK1\tOG0000123\tone-to-one\tCLV1\torthofinder",
            "G_RLK1\tOG0000123\tone-to-one\tCLV1\tomamer",
            "G_RLK2\tOG0000123\tone-to-many\tCLV2\torthofinder",
            "G_TF1\tOG0000456\tone-to-one\tWOX4\torthofinder",
            "G_CONFLICT\tOG0000789\tone-to-one\tPXY1\torthofinder",
        ]) + "\n", encoding="utf-8")
        paths["annotations.tsv"].write_text("\n".join([
            "gene_id\tdescription\tgo_terms",
            "G_RLK1\tleucine-rich repeat receptor kinase, meristem maintenance\tGO:0010073;GO:0004672",
            "G_RLK2\treceptor-like kinase\t",
            "G_TF1\thomeobox transcription factor\tGO:0006355",
            "G_TRANS1\tauxin efflux carrier component\t",
            "G_CONFLICT\treceptor kinase, vascular development\tGO:0010087",
        ]) + "\n", encoding="utf-8")
        paths["flags.tsv"].write_text("gene_id\tdark_state\nG_ORPHAN\tdark_unresolved\n", encoding="utf-8")
        return paths

    def _run(self, root: Path, extra: list[str] | None = None) -> Path:
        fixture = self._fixture(root)
        out = root / "atlas"
        argv = ["interval-atlas", "--candidates", str(fixture["candidates.tsv"]),
                "--family-index", str(fixture["family_index.tsv"]),
                "--domains", str(fixture["domains.tsv"]),
                "--orthology", str(fixture["orthology.tsv"]),
                "--annotations", str(fixture["annotations.tsv"]),
                "--flags", str(fixture["flags.tsv"]),
                "--out-dir", str(out), "--case-id", "demo"] + (extra or [])
        self.assertEqual(0, main(argv))
        return out

    def test_every_interval_gene_appears_exactly_once(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            atlas = pd.read_csv(out / "interval_gene_atlas.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(7, len(atlas))
            self.assertEqual(7, atlas["gene_id"].nunique())
            self.assertTrue(atlas["rank_contribution"].eq("none").all())

    def test_family_domain_and_orthogroup_layers_are_independent_columns(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            atlas = pd.read_csv(out / "interval_gene_atlas.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            row = atlas.loc["G_RLK1"]
            self.assertEqual("CLV", row.family_id)
            self.assertEqual("family_index", row.family_source)
            self.assertEqual("LRR-Pkinase", row.domain_architecture)
            self.assertEqual("OG0000123", row.orthogroup)
            self.assertEqual("one-to-one", row.orthology_relationship)
            # Two providers reaching the same orthogroup are recorded, not summed.
            self.assertEqual("omamer;orthofinder", row.orthology_providers)
            self.assertIn("family_index", row.classification_basis)
            self.assertIn("orthology", row.classification_basis)

    def test_family_index_subject_endpoints_and_search_uncertainty_are_preserved(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._fixture(root)
            family = fixture["family_index.tsv"]
            family.write_text("\n".join([
                "gene_id\tbest_subject\tbest_subject_symbol\tbest_subject_stem\tfamily_subjects\tassignment_status\tbest_identity\tbest_query_coverage\tbest_evalue\tbest_bitscore\tsecond_bitscore\ttop2_bitscore_ratio\tsymbol_source_rank\tsymbol_source_bitscore_ratio",
                "G_RLK1\tsp|P12345|KINASE_A\tKINASE_A\tKINASE\tsp|P12345|KINASE_A;tr|Q8XYZ1|KINASE_B\tassigned\t88.2\t0.91\t2e-40\t240\t221\t0.92\t2\t0.97",
                "G_RLK2\t\t\t\t\tdark_no_qualifying_hit\t\t\t\t\t\t\t\t",
                "G_RLK3\t\t\t\t\tdark_no_qualifying_hit\t\t\t\t\t\t\t\t",
                "G_TF1\t\t\t\t\tdark_no_qualifying_hit\t\t\t\t\t\t\t\t",
                "G_TRANS1\t\t\t\t\tdark_no_qualifying_hit\t\t\t\t\t\t\t\t",
                "G_ORPHAN\t\t\t\t\tdark_no_qualifying_hit\t\t\t\t\t\t\t\t",
                "G_CONFLICT\t\t\t\t\tdark_no_qualifying_hit\t\t\t\t\t\t\t\t",
            ]) + "\n", encoding="utf-8")
            out = root / "atlas"
            self.assertEqual(0, main(["interval-atlas", "--candidates", str(fixture["candidates.tsv"]),
                                      "--family-index", str(family), "--out-dir", str(out)]))
            atlas = pd.read_csv(out / "interval_gene_atlas.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            row = atlas.loc["G_RLK1"]
            self.assertEqual("sp|P12345|KINASE_A", row.family_index_best_subject)
            self.assertEqual("sp|P12345|KINASE_A;tr|Q8XYZ1|KINASE_B", row.family_index_subjects)
            self.assertEqual("P12345;Q8XYZ1", row.family_index_subject_endpoints)
            self.assertEqual("assigned", row.family_index_assignment_status)
            self.assertEqual("88.2", row.family_index_best_identity)
            self.assertEqual("221", row.family_index_second_bitscore)
            self.assertEqual("0.92", row.family_index_top2_bitscore_ratio)
            self.assertEqual("2", row.family_index_symbol_source_rank)
            self.assertEqual("none", row.rank_contribution)

    def test_issue_184_family_layers_are_separated(self) -> None:
        """Discovery bucket, evolutionary family, subfamily and copy set stay apart."""
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            atlas = pd.read_csv(out / "interval_gene_atlas.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            # Old-format family index: the coarse value is a discovery bucket
            # (symbol stem), and no curated evolutionary family exists.
            row = atlas.loc["G_RLK1"]
            self.assertEqual("CLV", row.discovery_bucket_id)
            self.assertEqual("CLV1", row.discovery_bucket_label)
            self.assertEqual("family_index", row.discovery_bucket_method)
            self.assertEqual("", row.evolutionary_family_id)
            self.assertEqual("", row.evolutionary_family_source)
            self.assertEqual("unresolved", row.evolutionary_family_status)
            # Subfamily is the finest honest proxy: the orthogroup, or unresolved.
            self.assertEqual("OG0000123", row.subfamily_id)
            self.assertEqual("orthogroup", row.subfamily_source)
            self.assertEqual("unresolved", atlas.loc["G_TRANS1", "subfamily_source"])
            # Copy set: homology layer id when present, the bucket as a
            # provisional placeholder when only a naming bucket exists.
            self.assertEqual("OG0000123", row.copy_set_id)
            self.assertEqual("validated", row.copy_set_status)
            self.assertEqual("OG0000456", atlas.loc["G_TF1", "copy_set_id"])
            self.assertEqual("CLV", atlas.loc["G_RLK3", "copy_set_id"])
            self.assertEqual("provisional", atlas.loc["G_RLK3", "copy_set_status"])
            self.assertEqual("", atlas.loc["G_ORPHAN", "copy_set_id"])
            self.assertEqual("unresolved", atlas.loc["G_ORPHAN", "copy_set_status"])
            manifest = json.loads((out / "interval_atlas_manifest.json").read_text(encoding="utf-8"))
            semantics = manifest["family_semantics"]
            for key in ("discovery_bucket", "evolutionary_family", "orthogroup", "subfamily", "copy_set", "gate", "dark"):
                self.assertIn(key, semantics)
            self.assertTrue(any("naming bucket alone never defines a valid copy set" in item for item in manifest["safety"]))

    def test_broad_function_bins_record_pattern_and_source(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            atlas = pd.read_csv(out / "interval_gene_atlas.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            row = atlas.loc["G_RLK1"]
            self.assertIn("meristem_development", row.broad_function_bins)
            self.assertIn("go_term:GO:0010073", row.broad_function_evidence)
            # A GO-backed bin outranks a keyword-only one as the primary.
            self.assertEqual("meristem_development", row.primary_broad_function)
            keyword_only = atlas.loc["G_TRANS1"]
            self.assertIn("transport", keyword_only.broad_function_bins)
            self.assertEqual("keyword", keyword_only.broad_function_source)
            self.assertIn("keyword:", keyword_only.broad_function_evidence)

    def test_dark_conflicting_and_unknown_genes_stay_as_buckets(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            atlas = pd.read_csv(out / "interval_gene_atlas.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            self.assertEqual("dark_or_short_orf", atlas.loc["G_ORPHAN", "evidence_completeness"])
            self.assertEqual("True", atlas.loc["G_ORPHAN", "annotation_rescue_needed"])
            # family index says BAM, orthology reference says PXY: recorded, not resolved.
            self.assertEqual("conflicting_annotation", atlas.loc["G_CONFLICT", "evidence_completeness"])
            self.assertIn("conflicting", atlas.loc["G_CONFLICT", "notes"])
            unclassified = pd.read_csv(out / "unclassified_dark_candidates.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"G_ORPHAN", "G_CONFLICT"}, set(unclassified["gene_id"]))
            self.assertTrue(unclassified["required_next_action"].ne("").all())

    def test_family_summary_counts_members_and_tandem_blocks(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            summary = pd.read_csv(out / "family_summary.tsv", sep="\t", dtype=str).fillna("").set_index("family_id")
            self.assertEqual("3", summary.loc["CLV", "member_count"])
            # RLK1/RLK2 are 4 kb apart, RLK3 is 400 kb away: two blocks.
            self.assertEqual("2", summary.loc["CLV", "tandem_blocks"])
            self.assertEqual("1", summary.loc["CLV", "orthology_one_to_one"])
            atlas = pd.read_csv(out / "interval_gene_atlas.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            self.assertEqual("True", atlas.loc["G_RLK1", "tandem_expansion"])
            self.assertEqual("False", atlas.loc["G_RLK3", "tandem_expansion"])
            self.assertEqual("2", atlas.loc["G_RLK1", "focal_paralog_count"])

    def test_review_queries_are_family_level_and_phenotype_free(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            queries = pd.read_csv(out / "family_review_queries.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(len(queries))
            self.assertTrue(queries["phenotype_used"].eq("False").all())
            self.assertIn("family", set(queries["scope"]))
            self.assertIn("broad_function", set(queries["scope"]))
            # Unknown genes get their own probe instead of being dropped.
            unknown = queries[queries.scope == "unclassified"].iloc[0]
            self.assertIn("G_ORPHAN", unknown.member_gene_ids)
            # Queries are gene-family / process scoped, not "which gene is causal".
            self.assertTrue(all("causal" not in value for value in queries["suggested_query"]))

    def test_manifest_records_phenotype_independence_and_atlas_fingerprint(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            manifest = json.loads((out / "interval_atlas_manifest.json").read_text(encoding="utf-8"))
            self.assertFalse(manifest["phenotype_used"])
            self.assertEqual(7, manifest["gene_count"])
            self.assertEqual(2, manifest["unclassified_count"])
            digest = hashlib.sha256((out / "interval_gene_atlas.tsv").read_bytes()).hexdigest()
            self.assertEqual(digest, manifest["atlas_sha256"])
            self.assertIn("sha256", manifest["broad_function_bin_map"])
            self.assertTrue(any("never causal evidence" in item for item in manifest["safety"]))

    def test_atlas_never_emits_a_score_column(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            for path in sorted(out.glob("*.tsv")):
                frame = pd.read_csv(path, sep="\t", dtype=str)
                for forbidden in ("dry_score", "score", "dry_group", "top_eligibility", "evidence_code"):
                    self.assertNotIn(forbidden, frame.columns, f"{path.name}:{forbidden}")

    def test_atlas_runs_without_any_optional_table(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._fixture(root)
            out = root / "bare"
            self.assertEqual(0, main(["interval-atlas", "--candidates", str(fixture["candidates.tsv"]), "--out-dir", str(out)]))
            atlas = pd.read_csv(out / "interval_gene_atlas.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(7, len(atlas))
            # With no annotation at all, every gene lands in an explicit unknown bucket
            # rather than being dropped. The short ORF is still recognised from the
            # candidate table's protein_length alone.
            self.assertEqual({"no_informative_annotation", "dark_or_short_orf"},
                             set(atlas["evidence_completeness"]))
            self.assertEqual("dark_or_short_orf",
                             atlas.set_index("gene_id").loc["G_ORPHAN", "evidence_completeness"])
            unclassified = pd.read_csv(out / "unclassified_dark_candidates.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(7, len(unclassified))

    def test_atlas_is_stable_for_the_same_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = self._run(root)
            fixture = self._fixture(root)
            again = root / "atlas_again"
            self.assertEqual(0, main(["interval-atlas", "--candidates", str(fixture["candidates.tsv"]),
                                      "--family-index", str(fixture["family_index.tsv"]),
                                      "--domains", str(fixture["domains.tsv"]),
                                      "--orthology", str(fixture["orthology.tsv"]),
                                      "--annotations", str(fixture["annotations.tsv"]),
                                      "--flags", str(fixture["flags.tsv"]),
                                      "--out-dir", str(again), "--case-id", "demo"]))
            for name in ("interval_gene_atlas", "family_summary", "broad_function_summary", "family_review_queries"):
                self.assertEqual((first / f"{name}.tsv").read_bytes(), (again / f"{name}.tsv").read_bytes(), name)


class ConservationLineageTests(unittest.TestCase):
    """Issue #187: a failed homology search is not evidence of lineage restriction.

    The safe default for a gene with no family / orthology / domain signal is
    ``no_detectable_homology``. ``lineage_restricted_supported`` requires explicit
    multi-reference negative-search coverage, which the current atlas inputs
    cannot supply, so the atlas never emits it.
    """

    LINEAGE_CLAIMS = {"conserved_cross_species", "relative_enriched", "lineage_restricted_supported"}

    def _write_inputs(self, root: Path, family_index: bool = True, domains: bool = True,
                      orthology: bool = True) -> dict[str, Path]:
        paths = {name: root / name for name in ("candidates.tsv", "family_index.tsv", "domains.tsv", "orthology.tsv")}
        paths["candidates.tsv"].write_text("\n".join([
            "gene_id\tchrom\tstart\tend\tprotein_length",
            "G_CONSERVED\tchr1\t1000\t3000\t980",
            "G_ENRICHED\tchr1\t20000\t22000\t420",
            "G_FAMILY_ONLY\tchr1\t30000\t31000\t510",
            "G_DOMAIN_ONLY\tchr1\t35000\t36000\t600",
            "G_ORPHAN\tchr1\t40000\t40300\t85",
        ]) + "\n", encoding="utf-8")
        if family_index:
            paths["family_index.tsv"].write_text("\n".join([
                "gene_id\tchrom\tstart\tend\tbest_subject\tbest_subject_symbol\tbest_subject_stem\tfamily_stems\tassignment_status",
                "G_CONSERVED\tchr1\t1000\t3000\tsp|P1|CLV1_ARATH\tCLV1\tCLV\tCLV\tassigned",
                "G_FAMILY_ONLY\tchr1\t30000\t31000\tsp|P5|PIN1_ARATH\tPIN1\tPIN\tPIN\tassigned",
                "G_ORPHAN\t\t\t\t\t\t\t\tdark_no_qualifying_hit",
            ]) + "\n", encoding="utf-8")
        if domains:
            paths["domains.tsv"].write_text("\n".join([
                "gene_id\tdomain_id\tdomain_name\tstart\tsource",
                "G_DOMAIN_ONLY\tPF00069\tPkinase\t600\tPfam",
            ]) + "\n", encoding="utf-8")
        if orthology:
            paths["orthology.tsv"].write_text("\n".join([
                "gene_id\torthogroup\trelationship\treference_gene_id\tprovider",
                "G_CONSERVED\tOG0000123\tone-to-one\tCLV1\torthofinder",
                "G_ENRICHED\tOG0000456\tone-to-many\tWOX4\torthofinder",
            ]) + "\n", encoding="utf-8")
        return paths

    def _run(self, root: Path, **include: bool) -> Path:
        fixture = self._write_inputs(root, **include)
        out = root / "atlas"
        argv = ["interval-atlas", "--candidates", str(fixture["candidates.tsv"]), "--out-dir", str(out)]
        if include.get("family_index", True):
            argv += ["--family-index", str(fixture["family_index.tsv"])]
        if include.get("domains", True):
            argv += ["--domains", str(fixture["domains.tsv"])]
        if include.get("orthology", True):
            argv += ["--orthology", str(fixture["orthology.tsv"])]
        self.assertEqual(0, main(argv))
        return out

    def _atlas(self, root: Path) -> pd.DataFrame:
        return pd.read_csv(root / "atlas" / "interval_gene_atlas.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")

    def test_reference_inputs_absent_means_no_detectable_homology_not_lineage_restricted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            self._run(Path(tmp), family_index=False, domains=False, orthology=False)
            atlas = self._atlas(Path(tmp))
            self.assertEqual({"no_detectable_homology"}, set(atlas["conservation_class"]))
            self.assertTrue(atlas["conservation_class"].isin(self.LINEAGE_CLAIMS).eq(False).all())
            self.assertTrue(atlas["rank_contribution"].eq("none").all())

    def test_gene_with_no_hit_despite_orthology_input_is_no_detectable_homology(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            self._run(Path(tmp))
            atlas = self._atlas(Path(tmp))
            self.assertEqual("no_detectable_homology", atlas.loc["G_ORPHAN", "conservation_class"])

    def test_dark_short_orf_carries_no_lineage_claim(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            self._run(Path(tmp))
            atlas = self._atlas(Path(tmp))
            row = atlas.loc["G_ORPHAN"]
            self.assertEqual("dark_or_short_orf", row.evidence_completeness)
            self.assertEqual("no_detectable_homology", row.conservation_class)
            self.assertNotIn(row.conservation_class, self.LINEAGE_CLAIMS)

    def test_family_only_and_domain_only_are_unresolved_not_lineage_restricted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            self._run(Path(tmp))
            atlas = self._atlas(Path(tmp))
            self.assertEqual("unresolved", atlas.loc["G_FAMILY_ONLY", "conservation_class"])
            self.assertEqual("unresolved", atlas.loc["G_DOMAIN_ONLY", "conservation_class"])
            self.assertNotIn("lineage_restricted_supported", set(atlas["conservation_class"]))

    def test_positive_conservation_classes_are_still_emitted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            self._run(Path(tmp))
            atlas = self._atlas(Path(tmp))
            self.assertEqual("conserved_cross_species", atlas.loc["G_CONSERVED", "conservation_class"])
            self.assertEqual("relative_enriched", atlas.loc["G_ENRICHED", "conservation_class"])

    def test_manifest_records_lineage_search_coverage_and_claim_basis(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            self._run(Path(tmp))
            manifest = json.loads((Path(tmp) / "atlas" / "interval_atlas_manifest.json").read_text(encoding="utf-8"))
            claim = manifest["lineage_claim"]
            for key in ("claim_basis", "reference_species_searched", "reference_roles", "search_method",
                        "thresholds", "provider_coverage", "no_hit_count", "outgroup_coverage"):
                self.assertIn(key, claim)
            self.assertEqual(list(atlas_module._CONSERVATION), claim["vocabulary"])
            self.assertNotIn("focal_lineage_restricted", claim["vocabulary"])
            self.assertIn("lineage_restricted_supported", claim["vocabulary"])
            self.assertEqual({"family_index", "orthology", "domains"}, set(claim["reference_roles"]))
            self.assertEqual({"orthofinder"}, set(claim["provider_coverage"]["orthology_providers"]))
            self.assertEqual(1, claim["no_hit_count"]["family_index"])
            self.assertEqual(3, claim["no_hit_count"]["orthology"])
            self.assertEqual(1, claim["no_hit_count"]["no_detectable_homology_genes"])
            self.assertTrue(claim["outgroup_coverage"].startswith("not_supplied"))
            self.assertTrue(any("failed homology search is not evidence" in item for item in manifest["safety"]))
            self.assertEqual("1.4", manifest["schema_version"])

    def test_lineage_restricted_supported_is_reserved_and_never_emitted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            self._run(Path(tmp))
            atlas = self._atlas(Path(tmp))
            self.assertNotIn("lineage_restricted_supported", set(atlas["conservation_class"]))


class GeneModelCompletenessTests(unittest.TestCase):
    """Issue #374 T3: every interval gene carries a gene-model completeness
    status derived from the candidate table's own columns only. The phenotype
    never reaches the atlas, so the status backs the "nothing was filtered on
    phenotype before model completeness was assessed" claim."""

    def _write(self, root: Path, content: str) -> Path:
        path = root / "candidates.tsv"
        path.write_text(content, encoding="utf-8")
        return path

    def _run(self, root: Path, candidates: Path) -> pd.DataFrame:
        out = root / "atlas"
        self.assertEqual(0, main(["interval-atlas", "--candidates", str(candidates), "--out-dir", str(out)]))
        return pd.read_csv(out / "interval_gene_atlas.tsv", sep="\t", dtype=str).fillna("")

    def _manifest(self, root: Path) -> dict[str, object]:
        return json.loads((root / "atlas" / "interval_atlas_manifest.json").read_text(encoding="utf-8"))

    def test_status_column_sits_last_in_the_atlas(self) -> None:
        self.assertEqual("gene_model_completeness", atlas_module.ATLAS_COLUMNS[-1])

    def test_derivation_branches_cover_the_full_vocabulary(self) -> None:
        """protein_length >= 1 -> full_model; cds_length / translation_status
        with an undeterminable protein length -> partial_model; neither -> no_model."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = self._write(root, "\n".join([
                "gene_id\tchrom\tstart\tend\tprotein_length\tcds_length\ttranslation_status",
                "G_FULL\tchr1\t1000\t3000\t980\t2940\tcomplete",
                "G_PART_CDS\tchr1\t5000\t7000\t\t1500\t",
                "G_PART_TRANS\tchr1\t9000\t11000\t\t\tpredicted",
                "G_ZERO_LEN\tchr1\t13000\t14000\t0\t\t",
                "G_NONE\tchr1\t15000\t16000\t\t\t",
            ]) + "\n")
            atlas = self._run(root, candidates).set_index("gene_id")
            self.assertEqual("full_model", atlas.loc["G_FULL", "gene_model_completeness"])
            self.assertEqual("partial_model", atlas.loc["G_PART_CDS", "gene_model_completeness"])
            self.assertEqual("partial_model", atlas.loc["G_PART_TRANS", "gene_model_completeness"])
            self.assertEqual("no_model", atlas.loc["G_ZERO_LEN", "gene_model_completeness"])
            self.assertEqual("no_model", atlas.loc["G_NONE", "gene_model_completeness"])
            # every emitted value is on the documented vocabulary, and the
            # fixture exercises all three statuses
            self.assertTrue(set(atlas["gene_model_completeness"]).issubset(atlas_module.GENE_MODEL_STATUSES))
            self.assertEqual(set(atlas_module.GENE_MODEL_STATUSES), set(atlas["gene_model_completeness"]))

    def test_missing_columns_fall_back_to_the_binary_status(self) -> None:
        """With only protein_length present the status degrades to full/no; the
        atlas must not invent a partial_model when the columns do not exist."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = self._write(root, "\n".join([
                "gene_id\tchrom\tstart\tend\tprotein_length",
                "G_FULL\tchr1\t1000\t3000\t980",
                "G_NONE\tchr1\t5000\t6000\t",
            ]) + "\n")
            atlas = self._run(root, candidates).set_index("gene_id")
            self.assertEqual("full_model", atlas.loc["G_FULL", "gene_model_completeness"])
            self.assertEqual("no_model", atlas.loc["G_NONE", "gene_model_completeness"])
            self.assertEqual({"full_model", "no_model"}, set(atlas["gene_model_completeness"]))

    def test_existing_full_fixture_genes_are_all_full_models_and_untouched(self) -> None:
        """The main fixture's genes all carry a protein length, so the new
        column is full_model everywhere and no existing value changed."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = IntervalAtlasTests()._run(root)
            atlas = pd.read_csv(out / "interval_gene_atlas.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"full_model"}, set(atlas["gene_model_completeness"]))
            self.assertEqual(7, len(atlas))

    def test_non_numeric_protein_length_is_no_model(self) -> None:
        """protein_length placeholders like "NA" / "?" cannot confirm a model:
        they are treated exactly like a missing length (no_model unless
        cds_length / translation_status upgrade to partial_model), and they
        never crash the short-ORF / dark detection."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = self._write(root, "\n".join([
                "gene_id\tchrom\tstart\tend\tprotein_length\tcds_length\ttranslation_status",
                "G_NA\tchr1\t1000\t3000\tNA\t\t",
                "G_NA_CDS\tchr1\t5000\t7000\tNA\t1500\t",
                "G_NA_TRANS\tchr1\t9000\t11000\tNA\t\tpredicted",
                "G_QUESTION\tchr1\t13000\t14000\t?\t\t",
                "G_FULL\tchr1\t15000\t16000\t980\t\t",
            ]) + "\n")
            atlas = self._run(root, candidates).set_index("gene_id")
            self.assertEqual("no_model", atlas.loc["G_NA", "gene_model_completeness"])
            self.assertEqual("partial_model", atlas.loc["G_NA_CDS", "gene_model_completeness"])
            self.assertEqual("partial_model", atlas.loc["G_NA_TRANS", "gene_model_completeness"])
            self.assertEqual("no_model", atlas.loc["G_QUESTION", "gene_model_completeness"])
            self.assertEqual("full_model", atlas.loc["G_FULL", "gene_model_completeness"])

    def test_manifest_aggregates_gene_model_counts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = self._write(root, "\n".join([
                "gene_id\tchrom\tstart\tend\tprotein_length\tcds_length\ttranslation_status",
                "G_FULL\tchr1\t1000\t3000\t980\t2940\tcomplete",
                "G_PART_CDS\tchr1\t5000\t7000\t\t1500\t",
                "G_PART_TRANS\tchr1\t9000\t11000\t\t\tpredicted",
                "G_NONE\tchr1\t15000\t16000\t\t\t",
            ]) + "\n")
            self._run(root, candidates)
            manifest = self._manifest(root)
            counts = manifest["gene_model_completeness_counts"]
            self.assertEqual({"full_model": 1, "partial_model": 2, "no_model": 1}, counts)
            self.assertEqual(manifest["gene_count"], sum(int(value) for value in counts.values()))
            self.assertEqual(list(atlas_module.GENE_MODEL_STATUSES),
                             manifest["gene_model_completeness_vocabulary"])

    def test_no_cds_translation_status_with_zero_cds_length_is_no_model(self) -> None:
        """Issue #374 T3 F3: cds_length="0" and translation_status="no_cds" both
        positively state there is no CDS. Before the fix, either non-empty
        string alone made the row read as partial_model; both must now read
        as no_model, individually and combined."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = self._write(root, "\n".join([
                "gene_id\tchrom\tstart\tend\tprotein_length\tcds_length\ttranslation_status",
                "G_BOTH\tchr1\t1000\t2000\t\t0\tno_cds",
                "G_LEN_ONLY\tchr1\t3000\t4000\t\t0\t",
                "G_STATUS_ONLY\tchr1\t5000\t6000\t\t\tno_cds",
                "G_STATUS_CASE\tchr1\t7000\t8000\t\t\tNO_CDS",
            ]) + "\n")
            atlas = self._run(root, candidates).set_index("gene_id")
            self.assertEqual("no_model", atlas.loc["G_BOTH", "gene_model_completeness"])
            self.assertEqual("no_model", atlas.loc["G_LEN_ONLY", "gene_model_completeness"])
            self.assertEqual("no_model", atlas.loc["G_STATUS_ONLY", "gene_model_completeness"])
            # translation_status is normalized (casefolded) before comparison.
            self.assertEqual("no_model", atlas.loc["G_STATUS_CASE", "gene_model_completeness"])

    def test_legitimate_complete_and_partial_translation_statuses_classify_correctly(self) -> None:
        """The real translation_status vocabulary this repo produces
        (reference_extract.write_candidates): a complete model is full_model
        via protein_length, and every positive translation_status that only
        ever accompanies a real CDS (translated / translated_complete_orf /
        translated_partial_or_unverified_orf / not_translated_invalid_phase /
        not_translated_length_not_multiple_of_3) keeps a gene with an
        undeterminable protein length at partial_model -- none of them is
        misread as no_model."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = self._write(root, "\n".join([
                "gene_id\tchrom\tstart\tend\tprotein_length\tcds_length\ttranslation_status",
                "G_COMPLETE\tchr1\t1000\t3000\t320\t963\ttranslated_complete_orf",
                "G_PARTIAL_ORF\tchr1\t5000\t6000\t\t500\ttranslated_partial_or_unverified_orf",
                "G_TRANSLATED\tchr1\t7000\t8000\t\t600\ttranslated",
                "G_INVALID_PHASE\tchr1\t9000\t10000\t\t400\tnot_translated_invalid_phase",
                "G_BAD_LENGTH\tchr1\t11000\t12000\t\t401\tnot_translated_length_not_multiple_of_3",
            ]) + "\n")
            atlas = self._run(root, candidates).set_index("gene_id")
            self.assertEqual("full_model", atlas.loc["G_COMPLETE", "gene_model_completeness"])
            self.assertEqual("partial_model", atlas.loc["G_PARTIAL_ORF", "gene_model_completeness"])
            self.assertEqual("partial_model", atlas.loc["G_TRANSLATED", "gene_model_completeness"])
            self.assertEqual("partial_model", atlas.loc["G_INVALID_PHASE", "gene_model_completeness"])
            self.assertEqual("partial_model", atlas.loc["G_BAD_LENGTH", "gene_model_completeness"])

    def test_malformed_cds_length_and_unrecognized_translation_status_stay_partial_model(self) -> None:
        """Issue #374 T3 F3: a non-numeric cds_length or an unrecognized
        translation_status must not be silently read as zero/absent -- the
        row keeps partial_model so the ambiguity stays visible rather than
        being asserted as no_model."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = self._write(root, "\n".join([
                "gene_id\tchrom\tstart\tend\tprotein_length\tcds_length\ttranslation_status",
                "G_BAD_CDS\tchr1\t1000\t2000\t\tnot_a_number\t",
                "G_UNKNOWN_STATUS\tchr1\t3000\t4000\t\t\tsome_unrecognized_value",
                "G_BAD_CDS_AND_UNKNOWN\tchr1\t5000\t6000\t\tN/A\tsome_unrecognized_value",
            ]) + "\n")
            atlas = self._run(root, candidates).set_index("gene_id")
            self.assertEqual("partial_model", atlas.loc["G_BAD_CDS", "gene_model_completeness"])
            self.assertEqual("partial_model", atlas.loc["G_UNKNOWN_STATUS", "gene_model_completeness"])
            self.assertEqual("partial_model", atlas.loc["G_BAD_CDS_AND_UNKNOWN", "gene_model_completeness"])


if __name__ == "__main__":
    unittest.main()
