"""Family candidate bridge (deterministic Route A/B/C normalization).

The bridge is a reporting-only table: rank_contribution is always ``none``,
there are no score columns, and every row carries fingerprints and a strict
mapping tier so promotion rules can be tested exactly.
"""

from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.family_candidate_bridge import BRIDGE_COLUMNS, dependency_cluster_id
from candidate_qtl.schemas import FAMILY_CANDIDATE_BRIDGE_SCHEMA, SchemaValidationError, validate

HYPOTHESES = "\n".join([
    "hypothesis_id\tgene_family\treference_gene_name\tproposed_stable_id\tmember_interval_genes",
    "H1\tFAM1\tREF1\tREF1_1\tG1",
    "H2\tFAM2\tREF2\tREF2_1\tG2;G3",
    "H3\tFAM3\tREF3\tREF3_1\tG7",
    "H4\tFAM4\tREF4\tREF4_1\tG9",
]) + "\n"

BAIT_PANEL = "\n".join([
    "hypothesis_id\tfunction_source_ids\tsource_provider",
    "H1\tROOT1\tuniprot_2026",
    "H2\tROOT1\tuniprot_2026",
]) + "\n"

FAMILY_STATUS = "\n".join([
    "hypothesis_id\tfamily_status\tinterval_gene_ids\tsource_ids\treference_gene_name\tproposed_stable_id",
    "H1\tconcrete_interval_gene_matched\tG1\tROOT1\tREF1\tREF1_1",
    "H2\tfamily_present_copy_unresolved\tG2;G3\tROOT1\tREF2\tREF2_1",
    "H4\tsurviving_sequence_only_annotation_gap\tG9\t\tREF4\tREF4_1",
]) + "\n"

ORTHOLOGY = "\n".join([
    "gene_id\treference_gene_id\trelationship\torthogroup\tpercent_identity\tquery_coverage\ttarget_coverage\tevalue\tprovider\tparalog_count\tsource_id",
    "G5\tREF5\tone_to_one\tOG5\t95\t0.9\t\t1e-50\tlegume3\t\tROOT1",
    "G5\tREF5\tone_to_one\tOG5\t95\t0.9\t\t1e-50\tlegume3\t\tROOT1",
    "G6\tREF6\tbest_hit_homology\tOG6\t\t0.3\t\t1e-5\tlegume3\t\t",
    "G4\tREF4\tone_to_many\tOG4\t\t0.4\t\t1e-10\tlegume3\t\t",
    "G9\tREF9\t\t\t\t\t\t\tlegume3\t\t",
    "G7\tREF7\tone_to_one\tOG7\t0\t0.8\t\t1e-40\tlegume3\t\t",
]) + "\n"

ATLAS = "\n".join([
    "gene_id\tfamily_id\tfamily_label\tevolutionary_family_id\torthogroup\tdomain_architecture\tdomain_ids\tdomain_source\tfocal_paralog_count\ttandem_expansion\torthology_reference_gene",
    "G4\tF4\tF4_label\t\t\t\t\t\t1\t\tREF4",
    "G5\t\t\tOG5\t\t\t\t\t1\t\tREF5",
    "G6\t\t\t\tOG6\t\t\t\t\t\tREF6",
]) + "\n"

MATCH = "\n".join([
    "reference_gene_id\trank\tmatched_interval_gene\tmatch_basis\tfamily_id\torthogroup_id",
    "R8\t1\tG8\tsame_orthogroup\t\tOG8",
    "R8\t2\tG8;G9\tsame_orthogroup\t\tOG8",
    "RC1\t1\tG1\tsame_family\tFAM1\t",
    "R9\t1\tG10\tsame_family\t\t",
]) + "\n"

MERGED = "\n".join([
    "reference_gene_id\trank\tdatabase_version\tmodel_version",
    "R8\t1\tdb2026\tesm3",
    "R8\t2\tdb2026\tesm3",
    "RC1\t1\tdb2026\tesm3",
    "R9\t1\tdb2026\tesm3",
]) + "\n"

GENES = "gene_id\n" + "\n".join([f"G{index}" for index in range(1, 11)]) + "\n"

_ROOT1_CLUSTER = "dep_" + hashlib.sha256(b"ROOT1").hexdigest()[:20]

CLUSTERS = "\n".join([
    "dependency_cluster\tdependency_status\troot_source_id\tsource_ids",
    f"{_ROOT1_CLUSTER}\troot_source_known\tROOT1\tROOT1;route_a_source",
]) + "\n"

# blast_evidence.tsv (blast-panel-search) columns the bridge joins by
# reference_accession == panel_id: REF1_1 is a unique threshold hit, REF2_1 is
# an ambiguous paralog, REF3_1 evaluated but found no threshold hit.
BLAST_EVIDENCE = "\n".join([
    "panel_id\tprovider_status\ttop_hit_target_id\ttop_hit_identity\ttop_hit_query_coverage\ttop_hit_target_coverage\ttop_hit_evalue\tparalog_count\tambiguity_status\treview_status\trank_contribution",
    "REF1_1\tevaluable\tG1\t99.0\t1.0\t1.0\t0.0\t0\tunique\tverified\tnone",
    "REF2_1\tevaluable\tG2\t98.0\t1.0\t1.0\t1e-50\t1\tambiguous_paralog\tverified\tnone",
    "REF3_1\tevaluable\t\t\t\t\t\t0\tno_threshold_hit\tverified\tnone",
]) + "\n"


class FamilyCandidateBridgeTests(unittest.TestCase):
    def _fixture(self, root: Path, with_clusters: bool = True) -> dict[str, Path]:
        paths = {name: root / name for name in (
            "hypotheses.tsv", "bait_panel.tsv", "family_status.tsv", "orthology.tsv",
            "atlas.tsv", "match.tsv", "merged.tsv", "genes.tsv", "clusters.tsv")}
        paths["hypotheses.tsv"].write_text(HYPOTHESES, encoding="utf-8")
        paths["bait_panel.tsv"].write_text(BAIT_PANEL, encoding="utf-8")
        paths["family_status.tsv"].write_text(FAMILY_STATUS, encoding="utf-8")
        paths["orthology.tsv"].write_text(ORTHOLOGY, encoding="utf-8")
        paths["atlas.tsv"].write_text(ATLAS, encoding="utf-8")
        paths["match.tsv"].write_text(MATCH, encoding="utf-8")
        paths["merged.tsv"].write_text(MERGED, encoding="utf-8")
        paths["genes.tsv"].write_text(GENES, encoding="utf-8")
        if with_clusters:
            paths["clusters.tsv"].write_text(CLUSTERS, encoding="utf-8")
        return paths

    def _run(
        self,
        root: Path,
        out_name: str = "out",
        with_clusters: bool = True,
        route_a: bool = True,
        route_b: bool = True,
        route_c: bool = True,
        with_genes: bool = False,
        fixture: dict[str, Path] | None = None,
    ) -> Path:
        if fixture is None:
            fixture = self._fixture(root, with_clusters=with_clusters)
        argv = ["family-candidate-bridge", "--out-dir", str(root / out_name), "--case-id", "demo"]
        if route_a:
            argv += ["--route-a-hypotheses", str(fixture["hypotheses.tsv"]),
                     "--route-a-family-status", str(fixture["family_status.tsv"]),
                     "--route-a-bait-panel", str(fixture["bait_panel.tsv"])]
        if route_b:
            argv += ["--route-b-orthology", str(fixture["orthology.tsv"]),
                     "--route-b-atlas", str(fixture["atlas.tsv"])]
        if route_c:
            argv += ["--route-c-retrieval", str(fixture["merged.tsv"]),
                     "--route-c-match", str(fixture["match.tsv"])]
        if with_clusters:
            argv += ["--dependency-clusters", str(fixture["clusters.tsv"])]
        if with_genes:
            argv += ["--genes", str(fixture["genes.tsv"])]
        self.assertEqual(0, main(argv))
        return root / out_name

    def _frame(self, out: Path) -> pd.DataFrame:
        return pd.read_csv(out / "family_candidate_bridge.tsv", sep="\t", dtype=str).fillna("")

    def _manifest(self, out: Path) -> dict[str, object]:
        return json.loads((out / "family_candidate_bridge_manifest.json").read_text(encoding="utf-8"))

    def _rows(self, out: Path, gene_id: str) -> pd.DataFrame:
        frame = self._frame(out)
        return frame[frame["gene_id"].eq(gene_id)]

    def test_strong_unique_source_backed_mapping_is_tier1(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            g1 = self._rows(out, "G1")
            g1 = g1[g1["route"].eq("a")]
            self.assertEqual(1, len(g1))
            self.assertEqual("tier1_source_backed_unique", g1.iloc[0]["mapping_tier"])
            self.assertEqual("resolved_single_copy", g1.iloc[0]["ambiguity_status"])
            self.assertEqual("complete", g1.iloc[0]["provider_status"])
            g5 = self._rows(out, "G5")
            rbh = g5[g5["match_type"].eq("rbh")]
            self.assertEqual(1, len(rbh))
            self.assertEqual("tier1_source_backed_unique", rbh.iloc[0]["mapping_tier"])
            self.assertEqual("resolved_single_copy", rbh.iloc[0]["ambiguity_status"])
            self.assertEqual("95", rbh.iloc[0]["identity"])
            self.assertEqual("0.9", rbh.iloc[0]["query_coverage"])

    def test_ambiguous_paralog_never_promoted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            for gene_id in ("G2", "G3"):
                row = self._rows(out, gene_id).iloc[0]
                self.assertEqual("a", row["route"])
                self.assertEqual("paralog_ambiguous", row["ambiguity_status"])
                self.assertEqual("tier2_source_backed_ambiguous", row["mapping_tier"])
            multi = self._frame(out)
            multi = multi[multi["route"].eq("c") & multi["paralog_count"].eq("2")]
            self.assertTrue(multi["ambiguity_status"].eq("paralog_ambiguous").all())
            self.assertTrue(multi["mapping_tier"].eq("tier4_sequence_match_ambiguous").all())
            for gene_id in ("G8", "G9"):
                row = self._rows(out, gene_id)
                row = row[row["paralog_count"].eq("2")].iloc[0]
                self.assertEqual("tier4_sequence_match_ambiguous", row["mapping_tier"])

    def test_name_only_is_strictly_weakest_evidence_tier(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            g4 = self._rows(out, "G4")
            name_only = g4[g4["match_type"].eq("name_only")]
            self.assertEqual(1, len(name_only))
            self.assertEqual("F4", name_only.iloc[0]["family_id"])
            self.assertEqual("tier5_name_only", name_only.iloc[0]["mapping_tier"])
            frame = self._frame(out)
            bucket_rows = frame[frame["match_type"].eq("name_only")]
            self.assertFalse(bucket_rows.empty)
            # name_only is strictly weakest: never above tier5, no matter what
            # the other columns carry.
            self.assertTrue(bucket_rows["mapping_tier"].eq("tier5_name_only").all())
            # The same gene also carries a real orthology row that stays above
            # the naming bucket row.
            rbh = g4[g4["match_type"].eq("rbh")]
            self.assertEqual(1, len(rbh))

    def test_provider_missing_rows_are_not_evaluated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            row = self._rows(out, "G7")
            row = row[row["route"].eq("a")].iloc[0]
            self.assertEqual("tier6_not_evaluated", row["mapping_tier"])
            self.assertEqual("not_evaluated", row["ambiguity_status"])
            self.assertEqual("not_run", row["provider_status"])
            self.assertEqual("not_evaluated_no_verified_bait", row["reason"])
            self.assertEqual("not_evaluated", row["review_status"])
            self.assertEqual("", row["match_type"])

    def test_coverage_below_threshold_never_tier1(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            g6 = self._rows(out, "G6")
            blastp = g6[g6["match_type"].eq("blastp")]
            self.assertEqual(1, len(blastp))
            row = blastp.iloc[0]
            self.assertEqual("tier6_not_evaluated", row["mapping_tier"])
            self.assertEqual("coverage_insufficient", row["ambiguity_status"])
            self.assertEqual("not_evaluated_coverage_below_threshold", row["reason"])
            self.assertEqual("not_evaluated", row["review_status"])
            self.assertNotEqual("tier1_source_backed_unique", row["mapping_tier"])
            # A source-backed row below threshold keeps its evidence at tier2:
            # ambiguity coverage_insufficient, no reason, review required.
            g4 = self._rows(out, "G4")
            rbh = g4[g4["match_type"].eq("rbh")]
            self.assertEqual(1, len(rbh))
            row = rbh.iloc[0]
            self.assertEqual("tier2_source_backed_ambiguous", row["mapping_tier"])
            self.assertEqual("coverage_insufficient", row["ambiguity_status"])
            self.assertEqual("", row["reason"])
            self.assertEqual("review_required", row["review_status"])
            self.assertNotEqual("tier1_source_backed_unique", row["mapping_tier"])
            g7 = self._rows(out, "G7")
            rbh7 = g7[g7["match_type"].eq("rbh")]
            self.assertEqual(1, len(rbh7))
            self.assertEqual("tier2_source_backed_ambiguous", rbh7.iloc[0]["mapping_tier"])
            self.assertEqual("coverage_insufficient", rbh7.iloc[0]["ambiguity_status"])

    def test_unclassified_annotation_gap_rows_not_evaluated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            row = self._rows(out, "G9")
            row = row[row["route"].eq("a")].iloc[0]
            self.assertEqual("tier6_not_evaluated", row["mapping_tier"])
            self.assertEqual("unclassified_gene", row["ambiguity_status"])
            self.assertEqual("complete", row["provider_status"])
            self.assertEqual("not_evaluated_unclassified_gene", row["reason"])
            self.assertEqual("", row["match_type"])
            self.assertEqual("not_evaluated", row["review_status"])

    def test_route_agreement_is_never_a_score(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            g1 = self._rows(out, "G1")
            g1 = g1[g1["family_id"].eq("FAM1")]
            # One Route A row and one Route C row for the same gene+family;
            # agreement produces two independent rows, never a fused score.
            self.assertEqual(2, len(g1))
            self.assertEqual({"a", "c"}, set(g1["route"]))
            self.assertTrue(g1["rank_contribution"].eq("none").all())
            frame = self._frame(out)
            self.assertEqual(BRIDGE_COLUMNS, list(frame.columns))
            self.assertFalse({"dry_score", "score", "dry_group", "rank"} & set(frame.columns))

    def test_a_b_same_source_shares_dependency_cluster(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            frame = self._frame(out)
            g1_a = frame[(frame["gene_id"].eq("G1")) & (frame["route"].eq("a"))].iloc[0]
            g5_b = frame[(frame["gene_id"].eq("G5")) & (frame["route"].eq("b")) & (frame["match_type"].eq("rbh"))].iloc[0]
            self.assertEqual(_ROOT1_CLUSTER, g1_a["dependency_cluster_id"])
            self.assertEqual(_ROOT1_CLUSTER, g5_b["dependency_cluster_id"])
            manifest = self._manifest(out)
            entries = {entry["dependency_cluster_id"]: entry for entry in manifest["dependency_clusters"]}
            self.assertIn(_ROOT1_CLUSTER, entries)
            self.assertEqual("a;b", entries[_ROOT1_CLUSTER]["routes"])
            self.assertEqual("ROOT1", entries[_ROOT1_CLUSTER]["root_source_id"])

    def test_source_id_never_auto_promoted_to_root(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp), with_clusters=False)
            frame = self._frame(out)
            g5_b = frame[(frame["gene_id"].eq("G5")) & (frame["route"].eq("b")) & (frame["match_type"].eq("rbh"))].iloc[0]
            # The B row's source_id ROOT1 is recorded but never promoted to a
            # root without the clusters table.
            self.assertEqual("ROOT1", g5_b["source_id"])
            self.assertEqual("dependency_unknown", g5_b["dependency_cluster_id"])
            # Route A verified bait with exactly one function_source_ids value
            # is a declared root.
            g1_a = frame[(frame["gene_id"].eq("G1")) & (frame["route"].eq("a"))].iloc[0]
            self.assertEqual(_ROOT1_CLUSTER, g1_a["dependency_cluster_id"])

    def test_route_b_orthogroup_source_resolves_only_through_the_clusters_table(self) -> None:
        """Issue #374 review: a Route B row without an explicit orthology
        source_id carries the orthogroup/HOG id as its source token. It is
        never auto-promoted: without the clusters table it stays
        dependency_unknown, and a route-b-prepared cluster table for that
        orthogroup resolves both the orthology row and the atlas orthogroup
        axis row to the same cluster."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._fixture(root, with_clusters=False)
            out = self._run(root, with_clusters=False, fixture=fixture)
            frame = self._frame(out)
            g6_blastp = frame[(frame["gene_id"].eq("G6")) & (frame["match_type"].eq("blastp"))].iloc[0]
            self.assertEqual("orthogroup:OG6", g6_blastp["source_id"])
            self.assertEqual("dependency_unknown", g6_blastp["dependency_cluster_id"])

            cluster = dependency_cluster_id("orthogroup:OG6")
            fixture["clusters.tsv"].write_text(
                "dependency_cluster\tdependency_status\troot_source_id\tsource_ids\n"
                f"{cluster}\troot_source_known\torthogroup:OG6\torthogroup:OG6\n",
                encoding="utf-8",
            )
            out2 = self._run(root, out_name="out2", fixture=fixture)
            frame2 = self._frame(out2)
            g6_blastp2 = frame2[(frame2["gene_id"].eq("G6")) & (frame2["match_type"].eq("blastp"))].iloc[0]
            self.assertEqual(cluster, g6_blastp2["dependency_cluster_id"])
            g6_orthogroup = frame2[(frame2["gene_id"].eq("G6")) & (frame2["match_type"].eq("orthogroup"))].iloc[0]
            self.assertEqual(cluster, g6_orthogroup["dependency_cluster_id"])

    def test_verified_bait_root_survives_unrelated_cluster_table(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._fixture(root)
            fixture["clusters.tsv"].write_text(
                "dependency_cluster\tdependency_status\troot_source_id\tsource_ids\n"
                "dep_unrelated\troot_source_known\tOTHER\tOTHER\n",
                encoding="utf-8",
            )
            out = self._run(root, fixture=fixture)
            frame = self._frame(out)
            g1_a = frame[(frame["gene_id"].eq("G1")) & (frame["route"].eq("a"))].iloc[0]
            self.assertEqual(_ROOT1_CLUSTER, g1_a["dependency_cluster_id"])
            manifest = self._manifest(out)
            clusters = {row["dependency_cluster_id"]: row for row in manifest["dependency_clusters"]}
            self.assertEqual("ROOT1", clusters[_ROOT1_CLUSTER]["root_source_id"])

    def test_route_b_atlas_emits_every_available_evidence_axis(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._fixture(root)
            with fixture["atlas.tsv"].open("a", encoding="utf-8", newline="") as stream:
                stream.write("G11\tlabel_bucket\tlabel\tEV11\tOG11\tPF1-PF2\tPF1;PF2\tinterpro\t1\t\tREF11\n")
            out = self._run(root, fixture=fixture)
            rows = self._rows(out, "G11")
            self.assertEqual(3, len(rows))
            self.assertEqual({"curated", "domain", "orthogroup"}, set(rows["match_type"]))

    def test_c_only_rescue_row_preserved(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            g8 = self._rows(out, "G8")
            single = g8[g8["paralog_count"].eq("1")]
            self.assertEqual(1, len(single))
            row = single.iloc[0]
            self.assertEqual("c", row["route"])
            self.assertEqual("orthogroup", row["match_type"])
            self.assertEqual("tier3_sequence_match_unique", row["mapping_tier"])
            self.assertEqual("db2026", row["source_id"])
            self.assertEqual("esm3", row["provider_version"])

    def test_empty_string_and_zero_are_distinguishable(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            frame = self._frame(out)
            g7_b = frame[(frame["gene_id"].eq("G7")) & (frame["route"].eq("b"))].iloc[0]
            self.assertEqual("0", g7_b["identity"])
            g6_b = frame[(frame["gene_id"].eq("G6")) & (frame["match_type"].eq("blastp"))].iloc[0]
            self.assertEqual("", g6_b["identity"])

    def test_deterministic_byte_identical_across_out_dirs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fixture = self._fixture(Path(tmp))
            first = self._run(Path(tmp), out_name="out1", fixture=fixture)
            second = self._run(Path(tmp), out_name="out2", fixture=fixture)
            self.assertEqual((first / "family_candidate_bridge.tsv").read_bytes(),
                             (second / "family_candidate_bridge.tsv").read_bytes())
            self.assertEqual((first / "family_candidate_bridge_manifest.json").read_bytes(),
                             (second / "family_candidate_bridge_manifest.json").read_bytes())

    def test_fingerprints_are_stable_and_path_free(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            fixture = self._fixture(Path(tmp))
            first = self._run(Path(tmp), out_name="out1", fixture=fixture)
            second = self._run(Path(tmp), out_name="out2", fixture=fixture)
            left = self._frame(first)
            right = self._frame(second)
            self.assertEqual(list(left["input_fingerprint"]), list(right["input_fingerprint"]))
            self.assertEqual(list(left["output_fingerprint"]), list(right["output_fingerprint"]))
            # output_fingerprint equals a recomputed hash over the 21 other
            # columns, sorted by column name, tab-joined.
            columns = sorted(column for column in BRIDGE_COLUMNS if column != "output_fingerprint")
            for record in left.to_dict("records"):
                text = "\t".join(str(record[column]) for column in columns)
                self.assertEqual(hashlib.sha256(text.encode("utf-8")).hexdigest(), record["output_fingerprint"])
            tsv_text = (first / "family_candidate_bridge.tsv").read_text(encoding="utf-8")
            manifest_text = (first / "family_candidate_bridge_manifest.json").read_text(encoding="utf-8")
            for text in (tsv_text, manifest_text):
                self.assertNotIn("C:\\", text)
                self.assertNotIn("/mnt/", text)
            manifest = json.loads(manifest_text)
            for name, filename in manifest["outputs"].items():
                if name.endswith("_sha256"):
                    continue
                self.assertNotIn("/", str(filename))
                self.assertNotIn("\\", str(filename))

    def test_schema_validation_via_schemas_validate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            frame = self._frame(out)
            validate(frame, FAMILY_CANDIDATE_BRIDGE_SCHEMA)
            broken = frame.copy()
            broken.loc[0, "mapping_tier"] = "tier9"
            with self.assertRaises(SchemaValidationError):
                validate(broken, FAMILY_CANDIDATE_BRIDGE_SCHEMA)
            broken = frame.copy()
            broken.loc[0, "rank_contribution"] = "scored"
            with self.assertRaises(SchemaValidationError):
                validate(broken, FAMILY_CANDIDATE_BRIDGE_SCHEMA)
            broken = frame.copy()
            broken.loc[0, "route"] = "x"
            with self.assertRaises(SchemaValidationError):
                validate(broken, FAMILY_CANDIDATE_BRIDGE_SCHEMA)

    def test_no_ranked_or_score_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            frame = self._frame(out)
            self.assertEqual(BRIDGE_COLUMNS, list(frame.columns))
            self.assertFalse({"dry_score", "score", "dry_group", "rank"} & set(frame.columns))
            self.assertTrue(frame["rank_contribution"].eq("none").all())
            manifest = self._manifest(out)
            self.assertEqual("none", manifest["rank_contribution"])

    def test_missing_optional_route_inputs_are_not_failures(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp), route_a=False, route_c=False)
            manifest = self._manifest(out)
            self.assertFalse(manifest["inputs"]["route_a_hypotheses"]["present"])
            self.assertFalse(manifest["inputs"]["route_c_match"]["present"])
            self.assertTrue(manifest["inputs"]["route_b_orthology"]["present"])
            self.assertEqual({"b": 7}, manifest["row_counts"]["by_route"])

    def test_route_b_atlas_is_a_sufficient_standalone_input(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._fixture(root, with_clusters=False)
            out = root / "out"
            self.assertEqual(0, main([
                "family-candidate-bridge", "--route-b-atlas", str(fixture["atlas.tsv"]),
                "--out-dir", str(out), "--case-id", "demo",
            ]))
            frame = self._frame(out)
            self.assertEqual(3, len(frame))
            self.assertEqual({"b"}, set(frame["route"]))

    def test_supplied_input_with_missing_columns_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            bad = Path(tmp) / "bad_hypotheses.tsv"
            bad.write_text("hypothesis_id\tgene_family\nH1\tFAM1\n", encoding="utf-8")
            with self.assertRaises(ValueError) as caught:
                main(["family-candidate-bridge", "--route-a-hypotheses", str(bad),
                      "--out-dir", str(Path(tmp) / "out"), "--case-id", "demo"])
            self.assertIn("member_interval_genes", str(caught.exception))

    def test_empty_case_id_fails_before_writing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._fixture(root)
            out = root / "out"
            with self.assertRaises(ValueError) as caught:
                from candidate_qtl.family_candidate_bridge import family_candidate_bridge
                family_candidate_bridge(str(out), route_b_atlas_arg=str(fixture["atlas.tsv"]), case_id="  ")
            self.assertIn("non-empty --case-id", str(caught.exception))
            self.assertFalse(out.exists())

    def test_route_b_row_without_family_skipped_and_counted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            manifest = self._manifest(out)
            self.assertEqual(1, manifest["skipped"]["route_b_rows_skipped_no_family"])
            frame = self._frame(out)
            self.assertTrue(frame[frame["gene_id"].eq("G9") & frame["route"].eq("b")].empty)

    def test_genes_without_bridge_rows_reported(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp), with_genes=True)
            manifest = self._manifest(out)
            self.assertEqual(["G10"], manifest["genes_without_bridge_rows"])

    def test_duplicate_evidence_rows_deduplicated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            frame = self._frame(out)
            duplicates = frame[frame["gene_id"].eq("G5") & frame["match_type"].eq("rbh")]
            self.assertEqual(1, len(duplicates))
            self.assertEqual(len(frame), len(frame.drop_duplicates()))

    def test_cli_main_entry(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            self.assertTrue((out / "family_candidate_bridge.tsv").exists())
            self.assertTrue((out / "family_candidate_bridge_manifest.json").exists())
            manifest = self._manifest(out)
            self.assertEqual("1.0", manifest["schema_version"])
            self.assertEqual("family-candidate-bridge", manifest["command"])
            self.assertEqual("demo", manifest["case_id"])

    def test_route_c_rows_without_family_skipped_and_counted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            manifest = self._manifest(out)
            self.assertEqual(1, manifest["skipped"]["route_c_rows_skipped_no_family"])
            frame = self._frame(out)
            self.assertTrue(frame[frame["gene_id"].eq("G10")].empty)

    def test_route_c_row_with_empty_interval_genes_skipped_and_counted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._fixture(root)
            fixture["match.tsv"].write_text("\n".join([
                "reference_gene_id\trank\tmatched_interval_gene\tmatch_basis\tfamily_id\torthogroup_id",
                "R8\t1\tG8\tsame_orthogroup\t\tOG8",
                "R8\t3\t\tsame_orthogroup\t\tOG8",
            ]) + "\n", encoding="utf-8")
            out = self._run(root, fixture=fixture)
            manifest = self._manifest(out)
            self.assertEqual(1, manifest["skipped"]["route_c_rows_skipped_empty_interval_genes"])
            frame = self._frame(out)
            self.assertFalse(frame[frame["route"].eq("c") & frame["gene_id"].eq("G8")].empty)

    def _run_with_blast(self, root: Path, blast: str) -> pd.DataFrame:
        fixture = self._fixture(root)
        blast_path = root / "blast_evidence.tsv"
        blast_path.write_text(blast, encoding="utf-8")
        out = root / "out"
        self.assertEqual(0, main([
            "family-candidate-bridge", "--out-dir", str(out), "--case-id", "demo",
            "--route-a-hypotheses", str(fixture["hypotheses.tsv"]),
            "--route-a-family-status", str(fixture["family_status.tsv"]),
            "--route-a-bait-panel", str(fixture["bait_panel.tsv"]),
            "--route-a-blast", str(blast_path),
        ]))
        return self._frame(out)

    def test_route_a_blast_evidence_is_transcribed_into_matched_rows(self) -> None:
        """Issue #374 Criterion 3: blast-panel-search identity / coverage /
        E-value / paralog ambiguity are joined onto Route A matched rows by
        reference_accession == panel_id. A blast hit never promotes and never
        changes rank_contribution."""
        with tempfile.TemporaryDirectory() as tmp:
            frame = self._run_with_blast(Path(tmp), BLAST_EVIDENCE)
            g1 = frame[(frame["gene_id"].eq("G1")) & frame["route"].eq("a")].iloc[0]
            self.assertEqual("99.0", g1["identity"])
            self.assertEqual("1.0", g1["query_coverage"])
            self.assertEqual("1.0", g1["target_coverage"])
            self.assertEqual("0.0", g1["evalue"])
            self.assertEqual("0", g1["paralog_count"])
            self.assertEqual("resolved_single_copy", g1["ambiguity_status"])
            self.assertEqual("tier1_source_backed_unique", g1["mapping_tier"])
            # An ambiguous-paralog blast verdict surfaces as copy ambiguity on a
            # matched row and its paralog count is transcribed.
            g2 = frame[(frame["gene_id"].eq("G2")) & frame["route"].eq("a")].iloc[0]
            self.assertEqual("98.0", g2["identity"])
            self.assertEqual("1", g2["paralog_count"])
            self.assertEqual("paralog_ambiguous", g2["ambiguity_status"])
            self.assertEqual("tier2_source_backed_ambiguous", g2["mapping_tier"])
            self.assertEqual("review_required", g2["review_status"])
            self.assertTrue(frame["rank_contribution"].eq("none").all())
            manifest = json.loads((Path(tmp) / "out" / "family_candidate_bridge_manifest.json")
                                  .read_text(encoding="utf-8"))
            self.assertTrue(manifest["route_a_blast"]["input_present"])
            self.assertEqual(3, manifest["route_a_blast"]["evidence_rows"])
            self.assertGreaterEqual(manifest["route_a_blast"]["rows_joined"], 2)

    def test_not_evaluated_blast_evidence_leaves_route_a_values_unchanged(self) -> None:
        """A blast provider gap (no verified bait / BLAST+ missing) is not a
        miss: the blast layer contributes no identity/coverage/evalue and the
        family-status evidence keeps its original tier."""
        with tempfile.TemporaryDirectory() as tmp:
            frame = self._run_with_blast(Path(tmp), "\n".join([
                "panel_id\tprovider_status\ttop_hit_target_id\ttop_hit_identity\ttop_hit_query_coverage\ttop_hit_target_coverage\ttop_hit_evalue\tparalog_count\tambiguity_status\treview_status\trank_contribution",
                "REF1_1\tnot_evaluated_blast_unavailable\t\t\t\t\t\t0\tnot_evaluated\tnot_evaluated\tnone",
            ]) + "\n")
            g1 = frame[(frame["gene_id"].eq("G1")) & frame["route"].eq("a")].iloc[0]
            self.assertEqual("", g1["identity"])
            self.assertEqual("", g1["query_coverage"])
            self.assertEqual("", g1["evalue"])
            self.assertEqual("resolved_single_copy", g1["ambiguity_status"])
            self.assertEqual("tier1_source_backed_unique", g1["mapping_tier"])


if __name__ == "__main__":
    unittest.main()
