"""Route C: embedding retrieval as a third discovery route (issue #193)."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import build_parser, main
from candidate_qtl.retrieval_route import triple_route_compare

RETRIEVAL = "\n".join([
    "reference_gene_id\treference_species\tfamily_id\torthogroup_id\tretrieval_score\tmatched_modality\tmatched_description_or_phenotype\tmodel_name\tmodel_version\tdatabase_version\trank",
    "AtCLV1\tArabidopsis thaliana\tCLV\tOG0001\t0.912\tprotein_text\treceptor kinase restricting the stem cell pool; meristem_development\tprotst\t0.1.0\tplant_db_2026-07\t1",
    "AtCLV2\tArabidopsis thaliana\tCLV\tOG0001\t0.887\tprotein_text\treceptor kinase; CLAVATA signalling\tprotst\t0.1.0\tplant_db_2026-07\t2",
    "AtPXY\tArabidopsis thaliana\tPXY\tOG0002\t0.854\tprotein_text\treceptor kinase for procambial cell division\tprotst\t0.1.0\tplant_db_2026-07\t3",
    "AtFT\tArabidopsis thaliana\tFT\tOG0003\t0.721\tprotein_text\tflowering locus T\tprotst\t0.1.0\tplant_db_2026-07\t4",
]) + "\n"

ATLAS = "\n".join([
    "gene_id\tchrom\tstart\tend\tfamily_id\tfamily_label\torthogroup\torthology_relationship\torthology_reference_gene\tdomain_architecture\tbroad_function_bins",
    "G_CLV1\tchr1\t1000\t3000\tCLV\tCLV1\tOG0001\tone-to-one\tAtCLV1\tLRR-Pkinase\tmeristem_development;signaling_receptor",
    "G_CLV2\tchr1\t5000\t7000\tCLV\tCLV2\tOG0001\tone-to-many\tAtCLV2\tLRR-Pkinase\tsignaling_receptor",
    "G_PXY\tchr1\t20000\t22000\tPXY\tPXY1\tOG0002\tone-to-one\tAtPXY\tLRR-Pkinase\tvascular_cell_wall;signaling_receptor",
    "G_DARK\tchr1\t40000\t40300\t\t\t\t\t\t\t",
]) + "\n"


def _write(root: Path, name: str, content: str) -> Path:
    path = root / name
    path.write_text(content, encoding="utf-8")
    return path


class RetrievalImportTests(unittest.TestCase):
    def test_import_freezes_contract_and_marks_the_lane_blind(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            retrieval = _write(root, "retrieval.tsv", RETRIEVAL)
            out = root / "route_c"
            self.assertEqual(0, main(["retrieval-import", "--retrieval", str(retrieval), "--out-dir", str(out),
                                      "--phenotype", "stem fasciation", "--case-id", "demo"]))
            frozen = pd.read_csv(out / "route_c_retrieval_frozen.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(4, len(frozen))
            self.assertEqual(["AtCLV1", "AtCLV2", "AtPXY", "AtFT"], list(frozen["reference_gene_id"]))
            self.assertEqual(["1", "2", "3", "4"], list(frozen["rank"]))
            status = json.loads((out / "route_c_retrieval_status.json").read_text(encoding="utf-8"))
            self.assertEqual("embedding_retrieval", status["route"])
            self.assertEqual("phenotype_only", status["lane"])
            self.assertFalse(status["interval_information_in_model_input"])
            self.assertTrue(status["blind_discovery_claim"])
            self.assertEqual([], status["interval_arguments"])
            self.assertEqual("none", status["rank_contribution"])

    def test_import_command_has_no_interval_input_flag(self) -> None:
        parser = build_parser()
        # The command line itself must not accept any way to pass interval
        # information: the retrieval is generated before the interval exists.
        with self.assertRaises(SystemExit):
            parser.parse_args(["retrieval-import", "--interval", "chr1:1-2",
                               "--retrieval", "x.tsv", "--out-dir", "y"])
        with self.assertRaises(SystemExit):
            parser.parse_args(["retrieval-import", "--candidates", "x.tsv",
                               "--retrieval", "x.tsv", "--out-dir", "y"])
        with self.assertRaises(SystemExit):
            parser.parse_args(["retrieval-import", "--atlas", "x.tsv",
                               "--retrieval", "x.tsv", "--out-dir", "y"])

    def test_import_rejects_columns_outside_the_contract(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            leaked = root / "leaked.tsv"
            leaked.write_text(
                "reference_gene_id\treference_species\tfamily_id\torthogroup_id\tretrieval_score\tmatched_modality\tmatched_description_or_phenotype\tmodel_name\tmodel_version\tdatabase_version\trank\tinterval_gene_id\n"
                "AtCLV1\tArabidopsis thaliana\tCLV\tOG0001\t0.912\tprotein_text\tstem cell pool\tprotst\t0.1.0\tdb\t1\tG_CLV1\n", encoding="utf-8")
            out = root / "route_c"
            with self.assertRaises(ValueError) as context:
                main(["retrieval-import", "--retrieval", str(leaked), "--out-dir", str(out)])
            self.assertIn("outside the retrieval contract", str(context.exception))
            self.assertIn("interval_gene_id", str(context.exception))

    def test_import_rejects_invalid_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "route_c"
            header = ("reference_gene_id\treference_species\tfamily_id\torthogroup_id\tretrieval_score\tmatched_modality\t"
                      "matched_description_or_phenotype\tmodel_name\tmodel_version\tdatabase_version\trank\n")
            cases = {
                "missing_columns.tsv": "reference_gene_id\trank\nAtCLV1\t1\n",
                "no_family.tsv": header + "AtCLV1\tArabidopsis thaliana\t\t\t0.9\tprotein_text\tdesc\tprotst\t0.1.0\tdb\t1\n",
                "bad_score.tsv": header + "AtCLV1\tArabidopsis thaliana\tCLV\tOG0001\tnan\tprotein_text\tdesc\tprotst\t0.1.0\tdb\t1\n",
                "bad_rank.tsv": header + "AtCLV1\tArabidopsis thaliana\tCLV\tOG0001\t0.9\tprotein_text\tdesc\tprotst\t0.1.0\tdb\t0\n",
                "fractional_rank.tsv": header + "AtCLV1\tArabidopsis thaliana\tCLV\tOG0001\t0.9\tprotein_text\tdesc\tprotst\t0.1.0\tdb\t1.5\n",
                "no_model.tsv": header + "AtCLV1\tArabidopsis thaliana\tCLV\tOG0001\t0.9\tprotein_text\tdesc\t\t0.1.0\tdb\t1\n",
                "duplicate_rank.tsv": header + "AtCLV1\tArabidopsis thaliana\tCLV\tOG0001\t0.9\tprotein_text\tdesc\tprotst\t0.1.0\tdb\t1\nAtCLV2\tArabidopsis thaliana\tCLV\tOG0001\t0.8\tprotein_text\tdesc\tprotst\t0.1.0\tdb\t1\n",
            }
            for name, content in cases.items():
                path = _write(root, name, content)
                with self.assertRaises(ValueError):
                    main(["retrieval-import", "--retrieval", str(path), "--out-dir", str(out)])


class RetrievalMergeTests(unittest.TestCase):
    def test_multi_database_merge_uses_within_db_rank_and_excludes_unannotated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            header = RETRIEVAL.splitlines()[0]
            db_a = _write(root, "db_a.tsv", "\n".join([
                header,
                "A_unannotated\tSpecies A\t\t\t999.0\tprotein_text\thigh raw score\tmodel\t1\tdb_a\t1",
                "A_PHY\tSpecies A\tPHY\tOG_PHY\t0.2\tprotein_text\tphytochrome\tmodel\t1\tdb_a\t2",
            ]) + "\n")
            db_b = _write(root, "db_b.tsv", "\n".join([
                header,
                "B_MLO\tSpecies B\tMLO\tOG_MLO\t0.99\tprotein_text\tmildew\tmodel\t1\tdb_b\t1",
                "B_FT\tSpecies B\tFT\tOG_FT\t0.98\tprotein_text\tflowering\tmodel\t1\tdb_b\t2",
            ]) + "\n")
            out = root / "merge"
            self.assertEqual(0, main(["retrieval-merge", "--retrieval", str(db_a), "--retrieval", str(db_b),
                                      "--top-k", "2", "--out-dir", str(out), "--case-id", "demo"]))
            merged = pd.read_csv(out / "route_c_retrieval_merged.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(["A_PHY", "B_MLO"], merged["reference_gene_id"].tolist())
            self.assertEqual(["1", "2"], merged["rank"].tolist())
            excluded = pd.read_csv(out / "route_c_retrieval_merge_excluded.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(["A_unannotated"], excluded["reference_gene_id"].tolist())
            self.assertEqual(["missing_family_or_orthogroup"], excluded["exclusion_reason"].tolist())
            status = json.loads((out / "route_c_retrieval_merge_status.json").read_text(encoding="utf-8"))
            self.assertFalse(status["cross_database_raw_scores_compared"])
            self.assertFalse(status["unannotated_rows_consume_top_k_slots"])
            self.assertEqual(1, status["per_database"]["db_a"]["selected_rows"])
            self.assertEqual(1, status["per_database"]["db_b"]["selected_rows"])
            frozen = root / "frozen"
            self.assertEqual(0, main(["retrieval-import", "--retrieval", str(out / "route_c_retrieval_merged.tsv"),
                                      "--out-dir", str(frozen)]))

    def test_merge_requires_multiple_distinct_databases(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            retrieval = _write(root, "retrieval.tsv", RETRIEVAL)
            with self.assertRaises(ValueError):
                main(["retrieval-merge", "--retrieval", str(retrieval), "--out-dir", str(root / "out")])

    def test_import_rejects_uncontracted_multi_database_table(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            lines = RETRIEVAL.splitlines()
            lines[2] = lines[2].replace("plant_db_2026-07", "other_db")
            raw_mix = _write(root, "raw_mix.tsv", "\n".join(lines) + "\n")
            with self.assertRaisesRegex(ValueError, "must be produced by retrieval-merge"):
                main(["retrieval-import", "--retrieval", str(raw_mix), "--out-dir", str(root / "out")])


class RetrievalIntervalMatchTests(unittest.TestCase):
    def test_match_happens_only_after_freeze_and_never_modifies_retrieval(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            retrieval = _write(root, "retrieval.tsv", RETRIEVAL)
            frozen_dir = root / "route_c"
            self.assertEqual(0, main(["retrieval-import", "--retrieval", str(retrieval), "--out-dir", str(frozen_dir)]))
            frozen_path = frozen_dir / "route_c_retrieval_frozen.tsv"
            before = frozen_path.read_bytes()
            atlas = _write(root, "atlas.tsv", ATLAS)
            out = root / "match"
            self.assertEqual(0, main(["retrieval-interval-match", "--retrieval", str(frozen_path), "--atlas", str(atlas),
                                      "--out-dir", str(out), "--case-id", "demo"]))
            # The frozen artifact is untouched by the interval step.
            self.assertEqual(before, frozen_path.read_bytes())
            matched = pd.read_csv(out / "route_c_interval_match.tsv", sep="\t", dtype=str).fillna("").set_index("reference_gene_id")
            self.assertEqual("G_CLV1;G_CLV2", matched.loc["AtCLV1", "matched_interval_gene"])
            self.assertIn("same_family", matched.loc["AtCLV1", "match_basis"])
            self.assertEqual("G_PXY", matched.loc["AtPXY", "matched_interval_gene"])
            self.assertEqual("", matched.loc["AtFT", "matched_interval_gene"])
            status = json.loads((out / "route_c_interval_match_status.json").read_text(encoding="utf-8"))
            self.assertEqual("after_freeze_only", status["interval_seen"])
            self.assertTrue(status["retrieval_frozen_before_interval"])
            self.assertTrue(status["retrieval_input_blind"])
            self.assertEqual("none", status["rank_contribution"])
            self.assertIsNone(status["zero_match_reason"])

    def test_match_uses_casefolded_orthology_reference_gene_without_alias_expansion(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            header = RETRIEVAL.splitlines()[0]
            retrieval = _write(root, "retrieval.tsv", "\n".join([
                header,
                "p12345\tReference species\tUNRELATED\tOG_UNRELATED\t0.9\tprotein_text\tdesc\tmodel\t1\tdb\t1",
                "P12345_ALIAS\tReference species\tOTHER\tOG_OTHER\t0.8\tprotein_text\tdesc\tmodel\t1\tdb\t2",
            ]) + "\n")
            frozen_dir = root / "route_c"
            self.assertEqual(0, main(["retrieval-import", "--retrieval", str(retrieval), "--out-dir", str(frozen_dir)]))
            frozen_path = frozen_dir / "route_c_retrieval_frozen.tsv"
            before = frozen_path.read_bytes()
            atlas = _write(root, "atlas.tsv", "\n".join([
                ATLAS.splitlines()[0],
                "Psat_TARGET\tchr1\t100\t200\tG7L4L7\t\t\t\tQ99999;P12345\t\t",
            ]) + "\n")
            out = root / "match"
            self.assertEqual(0, main(["retrieval-interval-match", "--retrieval", str(frozen_path), "--atlas", str(atlas),
                                      "--out-dir", str(out)]))
            self.assertEqual(before, frozen_path.read_bytes())
            matched = pd.read_csv(out / "route_c_interval_match.tsv", sep="\t", dtype=str).fillna("").set_index("reference_gene_id")
            self.assertEqual("Psat_TARGET", matched.loc["p12345", "matched_interval_gene"])
            self.assertEqual("same_orthology_reference_gene", matched.loc["p12345", "match_basis"])
            self.assertEqual("", matched.loc["P12345_ALIAS", "matched_interval_gene"])
            status = json.loads((out / "route_c_interval_match_status.json").read_text(encoding="utf-8"))
            self.assertEqual("none", status["rank_contribution"])
            self.assertEqual(1, status["zero_match_evidence"]["same_orthology_reference_gene"])

    def test_zero_match_reports_namespace_mismatch_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            header = RETRIEVAL.splitlines()[0]
            retrieval = _write(root, "retrieval.tsv", "\n".join([
                header,
                "P12345\tReference species\tMTR_FAMILY\tOG_MTR\t0.9\tprotein_text\tdesc\tmodel\t1\tdb\t1",
            ]) + "\n")
            frozen_dir = root / "route_c"
            self.assertEqual(0, main(["retrieval-import", "--retrieval", str(retrieval), "--out-dir", str(frozen_dir)]))
            atlas = _write(root, "atlas.tsv", "\n".join([
                ATLAS.splitlines()[0],
                "Psat_TARGET\tchr1\t100\t200\tG7L4L7\t\tOG_G7L4L7\t\tQ99999\t\t",
            ]) + "\n")
            out = root / "match"
            self.assertEqual(0, main(["retrieval-interval-match", "--retrieval", str(frozen_dir / "route_c_retrieval_frozen.tsv"),
                                      "--atlas", str(atlas), "--out-dir", str(out)]))
            status = json.loads((out / "route_c_interval_match_status.json").read_text(encoding="utf-8"))
            self.assertEqual(0, status["matched_retrieval_rows"])
            self.assertEqual("namespace_mismatch_suspected", status["zero_match_reason"])
            self.assertGreater(status["zero_match_evidence"]["atlas_orthology_reference_gene_key_count"], 0)
            self.assertGreater(status["zero_match_evidence"]["retrieval_reference_gene_id_count"], 0)
            self.assertEqual(0, status["zero_match_evidence"]["same_orthology_reference_gene"])
            self.assertEqual("none", status["rank_contribution"])

    def test_interval_match_requires_the_frozen_contract(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atlas = _write(root, "atlas.tsv", ATLAS)
            out = root / "match"
            bad = _write(root, "bad.tsv", "gene_id\trank\nG_X\t1\n")
            with self.assertRaises(ValueError):
                main(["retrieval-interval-match", "--retrieval", str(bad), "--atlas", str(atlas), "--out-dir", str(out)])


class TripleRouteCompareTests(unittest.TestCase):
    def _routes(self, root: Path, *, route_a_family: str, route_b_family: str = "CLV", retrieval_families: tuple[str, ...] = ("CLV", "PXY", "FT")) -> tuple[Path, Path, Path]:
        route_a = root / "route_a.tsv"
        route_a.write_text(
            "hypothesis_id\tgene_family\tmember_interval_genes\tsource_backing\n"
            f"hyp_001\t{route_a_family}\t{'G_CLV1;G_CLV2' if route_a_family == 'CLV' else ''}\tverified_bait\n", encoding="utf-8")
        route_b = root / "route_b.tsv"
        route_b.write_text(
            "hypothesis_id\tgene_family\tmember_interval_genes\tsource_backing\n"
            f"rev_001\t{route_b_family}\t{'G_CLV1;G_CLV2' if route_b_family == 'CLV' else ''}\tsource_backed\n", encoding="utf-8")
        rows = [
            "reference_gene_id\treference_species\tfamily_id\torthogroup_id\tretrieval_score\tmatched_modality\tmatched_description_or_phenotype\tmodel_name\tmodel_version\tdatabase_version\trank",
        ]
        for index, family in enumerate(retrieval_families, start=1):
            rows.append(f"Ref{family}\tArabidopsis thaliana\t{family}\tOG_{family}\t0.9\tprotein_text\tdesc\tprotst\t0.1.0\tdb\t{index}")
        route_c = _write(root, "route_c.tsv", "\n".join(rows) + "\n")
        return route_a, route_b, route_c

    def test_a_c_convergent_with_b_weak(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a, route_b, route_c = self._routes(root, route_a_family="CLV", route_b_family="PXY")
            out = root / "compare"
            self.assertEqual(0, main(["triple-route-compare", "--route-a", str(route_a), "--route-b", str(route_b),
                                      "--route-c", str(route_c), "--out-dir", str(out), "--case-id", "demo"]))
            convergence = json.loads((out / "triple_route_convergence.json").read_text(encoding="utf-8"))
            self.assertEqual("a_c_convergent", convergence["overall_status"])
            comparison = pd.read_csv(out / "triple_route_hypothesis_comparison.tsv", sep="\t", dtype=str).fillna("").set_index("family_id")
            self.assertEqual("a_c_convergent", comparison.loc["CLV", "convergence_status"])
            self.assertEqual("b_c_convergent", comparison.loc["PXY", "convergence_status"])
            self.assertEqual("c_only_supported", comparison.loc["FT", "convergence_status"])
            self.assertTrue(comparison["rank_contribution"].eq("none").all())

    def test_c_only_supported(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a, route_b, route_c = self._routes(root, route_a_family="", route_b_family="")
            out = root / "compare"
            self.assertEqual(0, main(["triple-route-compare", "--route-a", str(route_a), "--route-b", str(route_b),
                                      "--route-c", str(route_c), "--out-dir", str(out)]))
            convergence = json.loads((out / "triple_route_convergence.json").read_text(encoding="utf-8"))
            self.assertEqual("c_only_supported", convergence["overall_status"])
            self.assertIn("c_only_supported", set(pd.read_csv(out / "triple_route_review_required.tsv", sep="\t", dtype=str).fillna("")["reason"]))

    def test_all_weak_when_no_route_reaches_a_family(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            empty = _write(root, "empty.tsv", "hypothesis_id\tgene_family\tmember_interval_genes\tsource_backing\n")
            empty_c = _write(root, "empty_c.tsv", RETRIEVAL.split("\n", 1)[0] + "\n")
            out = root / "compare"
            self.assertEqual(0, main(["triple-route-compare", "--route-a", str(empty), "--route-b", str(empty),
                                      "--route-c", str(empty_c), "--out-dir", str(out)]))
            convergence = json.loads((out / "triple_route_convergence.json").read_text(encoding="utf-8"))
            self.assertEqual("all_weak", convergence["overall_status"])
            self.assertIn("all_routes_weak", set(pd.read_csv(out / "triple_route_review_required.tsv", sep="\t", dtype=str).fillna("")["reason"]))

    def test_three_way_convergence_is_recorded_without_scoring(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a, route_b, route_c = self._routes(root, route_a_family="CLV", route_b_family="CLV")
            out = root / "compare"
            self.assertEqual(0, main(["triple-route-compare", "--route-a", str(route_a), "--route-b", str(route_b),
                                      "--route-c", str(route_c), "--out-dir", str(out)]))
            convergence = json.loads((out / "triple_route_convergence.json").read_text(encoding="utf-8"))
            self.assertEqual("a_b_c_convergent", convergence["overall_status"])
            self.assertEqual(1, convergence["status_counts"].get("a_b_c_convergent"))
            self.assertIn("never enter a causal score", " ".join(convergence["safety"]))

    def test_convergence_on_shared_source_is_dependent_but_never_re_scored(self) -> None:
        """Issue #231: all three routes may rest on the same source id. The
        convergence label and rank_contribution stay untouched; the dependence
        is reported through the parallel source_dependence_status field."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a = _write(root, "route_a.tsv",
                             "hypothesis_id\tgene_family\tmember_interval_genes\tsource_backing\n"
                             "hyp_001\tCLV\tG_CLV1\tPMID:2001\n")
            route_b = _write(root, "route_b.tsv",
                             "hypothesis_id\tgene_family\tmember_interval_genes\tsource_backing\n"
                             "rev_001\tCLV\tG_CLV1\tPMID:2001\n")
            header = ("reference_gene_id\treference_species\tfamily_id\torthogroup_id\tretrieval_score\t"
                      "matched_modality\tmatched_description_or_phenotype\tmodel_name\tmodel_version\tdatabase_version\trank")
            route_c = _write(root, "route_c.tsv", "\n".join([
                header,
                "AtCLV1\tArabidopsis thaliana\tCLV\tOG0001\t0.9\tprotein_text\tdesc\tprotst\t0.1.0\tPMID:2001\t1",
            ]) + "\n")
            out = root / "compare"
            self.assertEqual(0, main(["triple-route-compare", "--route-a", str(route_a), "--route-b", str(route_b),
                                      "--route-c", str(route_c), "--out-dir", str(out), "--case-id", "demo"]))
            comparison = pd.read_csv(out / "triple_route_hypothesis_comparison.tsv", sep="\t", dtype=str).fillna("").set_index("family_id")
            row = comparison.loc["CLV"]
            self.assertEqual("a_b_c_convergent", row["convergence_status"])
            self.assertEqual("shared_source_dependent", row["source_dependence_status"])
            self.assertEqual("PMID:2001", row["route_a_source_ids"])
            self.assertEqual("PMID:2001", row["route_b_source_ids"])
            self.assertEqual("PMID:2001", row["route_c_source_ids"])
            self.assertEqual("none", row["rank_contribution"])
            self.assertIn("shared-source-dependent", row["notes"])
            self.assertIn("not counted as independent evidence", row["notes"])
            convergence = json.loads((out / "triple_route_convergence.json").read_text(encoding="utf-8"))
            self.assertEqual("a_b_c_convergent", convergence["overall_status"])
            dependence = convergence["source_dependence"]
            self.assertEqual(1, dependence["dependent_convergence_count"])
            self.assertEqual(0, dependence["independent_convergence_count"])
            self.assertEqual(["CLV"], dependence["families_dependent_on_shared_sources"])
            self.assertEqual({"1": 1}, dependence["converged_families_by_independent_source_count"])
            account = dependence["convergence_source_accounting"][0]
            self.assertEqual(3, account["routes_converged"])
            self.assertEqual(1, account["independent_sources"])
            self.assertEqual(["PMID:2001"], account["shared_source_ids"])
            self.assertEqual("shared_source_dependent", account["source_dependence_status"])
            self.assertIn("not counted as independent evidence", " ".join(convergence["safety"]))

    def test_per_route_source_ids_are_recorded_and_disjoint_ids_stay_independent(self) -> None:
        """Issue #231: each route's real source ids land in its own comparison
        column, and a convergence on disjoint ids is source_independent."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a = _write(root, "route_a.tsv",
                             "hypothesis_id\tgene_family\tmember_interval_genes\tsource_backing\n"
                             "hyp_001\tCLV\tG_CLV1\tPMID:1001;UniProt:Q8GRU6\n")
            route_b = _write(root, "route_b.tsv",
                             "hypothesis_id\tgene_family\tmember_interval_genes\tsource_backing\n"
                             "rev_001\tCLV\tG_CLV1\tPMID:2002\n")
            header = ("reference_gene_id\treference_species\tfamily_id\torthogroup_id\tretrieval_score\t"
                      "matched_modality\tmatched_description_or_phenotype\tmodel_name\tmodel_version\tdatabase_version\trank")
            route_c = _write(root, "route_c.tsv", "\n".join([
                header,
                "AtCLV1\tArabidopsis thaliana\tCLV\tOG0001\t0.9\tprotein_text\tdesc\tprotst\t0.1.0\tplant_db_2026-07\t1",
                "AtPXY\tArabidopsis thaliana\tPXY\tOG0002\t0.8\tprotein_text\tdesc\tprotst\t0.1.0\tplant_db_2026-07\t2",
            ]) + "\n")
            out = root / "compare"
            self.assertEqual(0, main(["triple-route-compare", "--route-a", str(route_a), "--route-b", str(route_b),
                                      "--route-c", str(route_c), "--out-dir", str(out)]))
            comparison = pd.read_csv(out / "triple_route_hypothesis_comparison.tsv", sep="\t", dtype=str).fillna("").set_index("family_id")
            clv = comparison.loc["CLV"]
            self.assertEqual("PMID:1001;UniProt:Q8GRU6", clv["route_a_source_ids"])
            self.assertEqual("PMID:2002", clv["route_b_source_ids"])
            self.assertEqual("plant_db_2026-07", clv["route_c_source_ids"])
            self.assertEqual("source_independent", clv["source_dependence_status"])
            self.assertEqual("a_b_c_convergent", clv["convergence_status"])
            self.assertEqual("none", clv["rank_contribution"])
            # Route C alone reaches PXY: no multi-route agreement to judge.
            pxy = comparison.loc["PXY"]
            self.assertEqual("", pxy["route_a_source_ids"])
            self.assertEqual("plant_db_2026-07", pxy["route_c_source_ids"])
            self.assertEqual("not_applicable", pxy["source_dependence_status"])
            convergence = json.loads((out / "triple_route_convergence.json").read_text(encoding="utf-8"))
            dependence = convergence["source_dependence"]
            self.assertEqual(1, dependence["independent_convergence_count"])
            self.assertEqual(0, dependence["dependent_convergence_count"])
            self.assertEqual({"4": 1}, dependence["converged_families_by_independent_source_count"])
            self.assertEqual("4", str(dependence["convergence_source_accounting"][0]["independent_sources"]))

    def test_descriptive_route_a_family_matches_symbol_route_c(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a, route_b, route_c = self._routes(root, route_a_family="Phytochrome (PHY)",
                                                     route_b_family="PXY", retrieval_families=("PHYA",))
            out = root / "compare"
            self.assertEqual(0, main(["triple-route-compare", "--route-a", str(route_a), "--route-b", str(route_b),
                                      "--route-c", str(route_c), "--out-dir", str(out)]))
            convergence = json.loads((out / "triple_route_convergence.json").read_text(encoding="utf-8"))
            self.assertEqual("a_c_convergent", convergence["overall_status"])

    def test_family_to_copy_handoff_is_a_shape_adapter_over_comparison_families(self) -> None:
        """Issue #283: one row per family with interval members; convergence
        status never filters the handoff; the path lands in the convergence
        manifest outputs."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a, route_b, route_c = self._routes(root, route_a_family="CLV", route_b_family="PXY")
            out = root / "compare"
            self.assertEqual(0, main(["triple-route-compare", "--route-a", str(route_a), "--route-b", str(route_b),
                                      "--route-c", str(route_c), "--out-dir", str(out)]))
            handoff = pd.read_csv(out / "family_to_copy_input.tsv", sep="\t", dtype=str).fillna("")
            comparison = pd.read_csv(out / "triple_route_hypothesis_comparison.tsv", sep="\t", dtype=str).fillna("").set_index("family_id")
            self.assertEqual(["hypothesis_id", "family_id", "gene_family", "member_interval_genes"], list(handoff.columns))
            self.assertEqual(1, len(handoff))
            row = handoff.iloc[0]
            self.assertEqual("CLV", row["family_id"])
            self.assertEqual("CLV", row["gene_family"])
            self.assertEqual(comparison.loc["CLV", "comparison_id"], row["hypothesis_id"])
            self.assertEqual("G_CLV1;G_CLV2", row["member_interval_genes"])
            convergence = json.loads((out / "triple_route_convergence.json").read_text(encoding="utf-8"))
            self.assertEqual(str(out / "family_to_copy_input.tsv"), convergence["outputs"]["family_to_copy_input"])
            # The comparison artifact itself is untouched.
            self.assertNotIn("family_to_copy_input", comparison.columns)

    def test_handoff_keeps_a_only_and_b_only_families_with_real_members(self) -> None:
        """Issue #283: a_only / b_only families pass through when their interval
        members exist; families without members are excluded."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a = _write(root, "route_a.tsv",
                             "hypothesis_id\tgene_family\tmember_interval_genes\tsource_backing\n"
                             "hyp_001\tCLV\tG_CLV1\tverified_bait\n"
                             "hyp_002\tWUS\t\tverified_bait\n")
            route_b = _write(root, "route_b.tsv",
                             "hypothesis_id\tgene_family\tmember_interval_genes\tsource_backing\n"
                             "rev_001\tCLV\tG_CLV1\tsource_backed\n"
                             "rev_002\tPXY\tG_PXY\tsource_backed\n")
            header = ("reference_gene_id\treference_species\tfamily_id\torthogroup_id\tretrieval_score\t"
                      "matched_modality\tmatched_description_or_phenotype\tmodel_name\tmodel_version\tdatabase_version\trank")
            route_c = _write(root, "route_c.tsv", "\n".join([
                header,
                "AtCLV1\tArabidopsis thaliana\tCLV\tOG0001\t0.9\tprotein_text\tdesc\tprotst\t0.1.0\tdb\t1",
                "AtFT\tArabidopsis thaliana\tFT\tOG0003\t0.7\tprotein_text\tdesc\tprotst\t0.1.0\tdb\t2",
            ]) + "\n")
            out = root / "compare"
            result = triple_route_compare(str(route_a), str(route_b), str(route_c), str(out), case_id="demo")
            handoff = pd.read_csv(out / "family_to_copy_input.tsv", sep="\t", dtype=str).fillna("").set_index("gene_family")
            self.assertEqual({"CLV", "PXY"}, set(handoff.index))
            self.assertEqual("G_CLV1", handoff.loc["CLV", "member_interval_genes"])
            self.assertEqual("G_PXY", handoff.loc["PXY", "member_interval_genes"])
            self.assertEqual(str(out / "family_to_copy_input.tsv"), result["family_to_copy_input"])

    def test_handoff_c_only_family_with_interval_match_is_kept(self) -> None:
        """Issue #283: a C-only family whose retrieval reached interval genes is
        handed over; route agreement is not a filter."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a = _write(root, "route_a.tsv", "hypothesis_id\tgene_family\tmember_interval_genes\tsource_backing\n")
            route_b = _write(root, "route_b.tsv", "hypothesis_id\tgene_family\tmember_interval_genes\tsource_backing\n")
            header = ("reference_gene_id\treference_species\tfamily_id\torthogroup_id\tretrieval_score\t"
                      "matched_modality\tmatched_description_or_phenotype\tmodel_name\tmodel_version\tdatabase_version\trank")
            route_c = _write(root, "route_c.tsv", "\n".join([
                header,
                "AtFT\tArabidopsis thaliana\tFT\tOG0003\t0.7\tprotein_text\tdesc\tprotst\t0.1.0\tdb\t1",
            ]) + "\n")
            interval_match = _write(root, "match.tsv",
                                    "reference_gene_id\tfamily_id\torthogroup_id\tmatched_interval_gene\n"
                                    "AtFT\tFT\t\tG_FT1\n")
            out = root / "compare"
            result = triple_route_compare(str(route_a), str(route_b), str(route_c), str(out), case_id="demo",
                                          route_c_match_arg=str(interval_match))
            handoff = pd.read_csv(out / "family_to_copy_input.tsv", sep="\t", dtype=str).fillna("").set_index("gene_family")
            self.assertEqual({"FT"}, set(handoff.index))
            self.assertEqual("G_FT1", handoff.loc["FT", "member_interval_genes"])
            self.assertEqual(str(out / "family_to_copy_input.tsv"), result["family_to_copy_input"])

    def test_handoff_dedupes_and_sorts_member_ids(self) -> None:
        """Issue #283: duplicate member ids are removed and the list is sorted."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a = _write(root, "route_a.tsv",
                             "hypothesis_id\tgene_family\tmember_interval_genes\tsource_backing\n"
                             "hyp_001\tCLV\tG_CLV2;G_CLV1;G_CLV2\tverified_bait\n")
            route_b = _write(root, "route_b.tsv",
                             "hypothesis_id\tgene_family\tmember_interval_genes\tsource_backing\n"
                             "rev_001\tCLV\tG_CLV2\tsource_backed\n")
            route_c = _write(root, "route_c.tsv", RETRIEVAL)
            out = root / "compare"
            triple_route_compare(str(route_a), str(route_b), str(route_c), str(out), case_id="demo")
            handoff = pd.read_csv(out / "family_to_copy_input.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(1, len(handoff))
            self.assertEqual("G_CLV1;G_CLV2", handoff.iloc[0]["member_interval_genes"])

    def test_handoff_is_byte_stable(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a, route_b, route_c = self._routes(root, route_a_family="CLV", route_b_family="PXY")
            first = root / "compare1"
            second = root / "compare2"
            for out in (first, second):
                self.assertEqual(0, main(["triple-route-compare", "--route-a", str(route_a), "--route-b", str(route_b),
                                          "--route-c", str(route_c), "--out-dir", str(out)]))
            self.assertEqual((first / "family_to_copy_input.tsv").read_bytes(),
                             (second / "family_to_copy_input.tsv").read_bytes())

    def test_handoff_is_absent_when_no_family_has_interval_members(self) -> None:
        """Issue #283/#284: no interval members anywhere -> no handoff artifact,
        so a downstream stage sees an upstream gap, not an empty success."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a = _write(root, "route_a.tsv",
                             "hypothesis_id\tgene_family\tmember_interval_genes\tsource_backing\n"
                             "hyp_002\tWUS\t\tverified_bait\n")
            route_b = _write(root, "route_b.tsv",
                             "hypothesis_id\tgene_family\tmember_interval_genes\tsource_backing\n")
            route_c = _write(root, "route_c.tsv", RETRIEVAL)
            out = root / "compare"
            result = triple_route_compare(str(route_a), str(route_b), str(route_c), str(out), case_id="demo")
            self.assertFalse((out / "family_to_copy_input.tsv").exists())
            self.assertNotIn("family_to_copy_input", result)
            convergence = json.loads((out / "triple_route_convergence.json").read_text(encoding="utf-8"))
            self.assertNotIn("family_to_copy_input", convergence["outputs"])


class RouteCBenchmarkTests(unittest.TestCase):
    def test_family_and_gene_recall_with_route_comparison(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hit = _write(root, "hit.tsv", RETRIEVAL)
            a_clv = _write(root, "a_clv.tsv", "hypothesis_id\tgene_family\nhyp_001\tCLV\n")
            a_ft = _write(root, "a_ft.tsv", "hypothesis_id\tgene_family\nhyp_001\tFT\n")
            manifest = _write(root, "manifest.tsv", "\n".join([
                "case_id\tretrieval_tsv\ttruth_family_id\ttruth_gene_id\troute_a_tsv",
                f"case_hit\t{hit}\tCLV\tAtCLV1\t{a_clv}",
                f"case_miss\t{hit}\tZZZ\tAtZZZ\t{a_ft}",
                f"case_a_only\t{hit}\tCLV\t\t{a_clv}",
            ]) + "\n")
            out = root / "bench"
            self.assertEqual(0, main(["route-c-benchmark", "--manifest", str(manifest), "--out-dir", str(out)]))
            cases = pd.read_csv(out / "route_c_benchmark_cases.tsv", sep="\t", dtype=str).fillna("").set_index("case_id")
            self.assertEqual("True", cases.loc["case_hit", "family_recall_5"])
            self.assertEqual("True", cases.loc["case_hit", "family_recall_20"])
            self.assertEqual("True", cases.loc["case_hit", "gene_recall_20"])
            self.assertEqual("1", cases.loc["case_hit", "c_family_rank"])
            self.assertEqual("False", cases.loc["case_miss", "family_recall_20"])
            self.assertEqual("False", cases.loc["case_miss", "c_only_hit_20"])
            # A already had CLV: A-only hit, C-only miss; union still hits.
            self.assertEqual("True", cases.loc["case_a_only", "a_only_hit"])
            self.assertEqual("False", cases.loc["case_a_only", "c_only_hit_20"])
            self.assertEqual("True", cases.loc["case_a_only", "ac_union_hit_20"])
            manifest_out = json.loads((out / "route_c_benchmark_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(3, manifest_out["case_count"])
            self.assertIn("no interval gene list", manifest_out["interval_blind"]["declaration"])
            self.assertEqual([], manifest_out["interval_blind"]["interval_arguments"])
            self.assertIn("never as causal evidence", manifest_out["scoring_rule"])

    def test_manifest_validation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bad = _write(root, "bad.tsv", "case_id\ttruth_family_id\ncase_1\tCLV\n")
            out = root / "bench"
            with self.assertRaises(ValueError) as context:
                main(["route-c-benchmark", "--manifest", str(bad), "--out-dir", str(out)])
            self.assertIn("retrieval_tsv", str(context.exception))

    def test_family_recall_uses_shared_descriptive_aliases(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hit = _write(root, "hit.tsv", RETRIEVAL.replace("\tCLV\tOG0001", "\tPHYA\tOG0001", 1))
            route_a = _write(root, "a.tsv", "hypothesis_id\tgene_family\nhyp_1\tPhytochrome (PHY)\n")
            manifest = _write(root, "manifest.tsv", "\n".join([
                "case_id\tretrieval_tsv\ttruth_family_id\troute_a_tsv",
                f"phy\t{hit}\tPhytochrome (PHY)\t{route_a}",
            ]) + "\n")
            out = root / "bench"
            self.assertEqual(0, main(["route-c-benchmark", "--manifest", str(manifest), "--out-dir", str(out)]))
            case = pd.read_csv(out / "route_c_benchmark_cases.tsv", sep="\t", dtype=str).fillna("").iloc[0]
            self.assertEqual("True", case["family_recall_5"])
            self.assertEqual("True", case["a_only_hit"])


class RouteCDenominatorTests(unittest.TestCase):
    """Optional inputs that were never supplied are not_evaluated, not misses (#214)."""

    def _run(self, root: Path, manifest_rows: list[str]) -> dict:
        manifest = _write(root, "manifest.tsv", "\n".join(manifest_rows) + "\n")
        out = root / "bench"
        self.assertEqual(0, main(["route-c-benchmark", "--manifest", str(manifest), "--out-dir", str(out)]))
        return json.loads((out / "route_c_benchmark_manifest.json").read_text(encoding="utf-8"))

    def test_missing_truth_gene_and_route_a_are_not_zero_misses(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hit = _write(root, "hit.tsv", RETRIEVAL)
            manifest = self._run(root, [
                "case_id\tretrieval_tsv\ttruth_family_id",
                f"case1\t{hit}\tCLV",
            ])
            metrics = manifest["metrics"]
            self.assertEqual(1.0, metrics["family_recall_20"]["recall"])
            for name in ("gene_recall_20", "a_only_hit", "c_only_hit_20", "ac_union_hit_20"):
                self.assertIsNone(metrics[name]["recall"], name)
                self.assertEqual(0, metrics[name]["evaluated"], name)
                self.assertEqual(1, metrics[name]["not_evaluated"], name)
            report = (root / "bench" / "route_c_benchmark_report.md").read_text(encoding="utf-8")
            gene_line = next(line for line in report.splitlines() if line.startswith("- Gene Recall@20"))
            self.assertIn("not_evaluated", gene_line)
            self.assertNotIn("%", gene_line)

    def test_mixed_manifest_uses_per_metric_denominators(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hit = _write(root, "hit.tsv", RETRIEVAL)
            route_a = _write(root, "a_clv.tsv", "hypothesis_id\tgene_family\nhyp_001\tCLV\n")
            manifest = self._run(root, [
                "case_id\tretrieval_tsv\ttruth_family_id\ttruth_gene_id\troute_a_tsv",
                f"with_gene\t{hit}\tCLV\tAtCLV1\t",
                f"with_route_a\t{hit}\tCLV\t\t{route_a}",
                f"family_only\t{hit}\tCLV\t\t",
            ])
            metrics = manifest["metrics"]
            self.assertEqual({"evaluated": 3, "hits": 3}, {key: metrics["family_recall_20"][key] for key in ("evaluated", "hits")})
            self.assertEqual(1, metrics["gene_recall_20"]["evaluated"])
            self.assertEqual(2, metrics["gene_recall_20"]["not_evaluated"])
            self.assertEqual(1.0, metrics["gene_recall_20"]["recall"])
            self.assertEqual(1, metrics["a_only_hit"]["evaluated"])
            self.assertEqual(1.0, metrics["ac_union_hit_20"]["recall"])
            summary = pd.read_csv(root / "bench" / "route_c_benchmark_summary.tsv", sep="\t", dtype=str).fillna("")
            # Issue #367 additive coverage columns (clean_cases / contaminated_cases /
            # evaluation_coverage / wilson_95_low / wilson_95_high).
            self.assertEqual({"metric", "requires", "cases", "evaluated", "hits", "recall", "not_evaluated",
                              "clean_cases", "contaminated_cases", "evaluation_coverage",
                              "wilson_95_low", "wilson_95_high"},
                             set(summary.columns))
            self.assertEqual("truth_gene_id", summary.set_index("metric").loc["gene_recall_20", "requires"])


if __name__ == "__main__":
    unittest.main()
