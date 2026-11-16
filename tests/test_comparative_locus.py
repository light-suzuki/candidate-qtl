from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from candidate_qtl.cli import main


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


class ComparativeLocusReconstructionTests(unittest.TestCase):
    def _inputs(self, root: Path, *, discordant: bool = False) -> dict[str, Path]:
        paths = {
            "markers": _write(root / "markers.fa", ">left\nACGTACGT\n>right\nTTGGCCAA\n"),
            "panel": _write(root / "panel.tsv", "\n".join([
                "species\trole\tavailable\tcache_version\tphylogenetic_distance",
                "RefA\tclose_reference\ttrue\tcache-a-v1\t0.1",
                "RefB\tclose_reference\ttrue\tcache-b-v2\t0.2",
                "RefC\toutgroup\ttrue\tcache-c-v1\t0.5",
            ]) + "\n"),
            "hits": _write(root / "hits.tsv", "\n".join([
                "marker_id\tspecies\tseqid\tstart\torientation\thit_count",
                "left\tRefA\tchr1\t100\t+\t1", "right\tRefA\tchr1\t900\t+\t1",
                "left\tRefB\tchr2\t100\t+\t1", "right\tRefB\tchr2\t900\t+\t1",
                "left\tRefC\tchr3\t900\t-\t1", "right\tRefC\tchr3\t100\t-\t1",
            ]) + "\n"),
            "native": _write(root / "native.tsv", "gene_id\torthogroup\nNativeA\tOG_A\n"),
            "focal": _write(root / "focal.tsv", "\n".join([
                "orthogroup\tfocal_gene_id\tsource_class\tsource_id\tidentity\tquery_coverage",
                "OG_B\tScafGeneB\tscaffold_recovered\tblast:hit-b\t88\t0.95",
                "OG_C\tTxC\ttranscript_supported\ttranscript:tx-c\t91\t0.80",
            ]) + "\n"),
            "absence": _write(root / "absence.tsv", "orthogroup\tevidence_status\tsource_id\nOG_D\tcallable_absent\tparents:callable\n"),
        }
        ref_b_order = ["OG_A", "OG_C", "OG_B", "OG_D"] if discordant else ["OG_A", "OG_B", "OG_C", "OG_D"]
        lines = ["species\tgene_id\torthogroup\torder_index"]
        for species, order in (("RefA", ["OG_A", "OG_B", "OG_C"]), ("RefB", ref_b_order), ("RefC", ["OG_D", "OG_C", "OG_B", "OG_A"])):
            for index, group in enumerate(order, 1):
                lines.append(f"{species}\t{species}_{group}\t{group}\t{index}")
        paths["genes"] = _write(root / "genes.tsv", "\n".join(lines) + "\n")
        return paths

    def _run(self, root: Path, *, discordant: bool = False) -> Path:
        inputs = self._inputs(root, discordant=discordant)
        out = root / "out"
        self.assertEqual(0, main([
            "comparative-locus-reconstruct", "--markers", str(inputs["markers"]),
            "--reference-panel", str(inputs["panel"]), "--marker-projections", str(inputs["hits"]),
            "--reference-locus-genes", str(inputs["genes"]), "--native-genes", str(inputs["native"]),
            "--focal-recovery", str(inputs["focal"]), "--absence-evidence", str(inputs["absence"]),
            "--out-dir", str(out), "--case-id", "draft-qtl",
        ]))
        return out

    def test_marker_only_multi_species_reconstruction_separates_membership_from_causality(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            manifest = json.loads((out / "comparative_reconstruction_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("consensus_strong", manifest["locus_reconstruction_status"])
            self.assertTrue(manifest["marker_sequence_start"])
            self.assertFalse(manifest["native_chromosome_coordinate_required"])
            self.assertFalse(manifest["standard_lane_dependency"])
            self.assertEqual(3, len(manifest["reference_species_with_gene_order"]))
            self.assertIn("cache-b-v2", manifest["reference_cache_versions"].values())
            universe = pd.read_csv(out / "reconstructed_candidate_universe.tsv", sep="\t", dtype=str).fillna("")
            by_group = universe.set_index("orthogroup")
            self.assertEqual("native_annotated", by_group.loc["OG_A", "candidate_class"])
            self.assertEqual("scaffold_recovered", by_group.loc["OG_B", "candidate_class"])
            self.assertEqual("transcript_supported", by_group.loc["OG_C", "candidate_class"])
            self.assertEqual("presence_absence_uncertain", by_group.loc["OG_D", "candidate_class"])
            self.assertEqual("expected:OG_D", by_group.loc["OG_D", "gene_id"])
            self.assertFalse(universe["gene_id"].eq("").any())
            self.assertFalse(universe["gene_id"].duplicated().any())
            self.assertEqual({"none"}, set(universe["rank_contribution"]))
            self.assertEqual({"not_assessed_here"}, set(universe["causal_confidence"]))
            expected = pd.read_csv(out / "expected_orthogroups.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"conserved_core", "lineage_variable"}, set(expected["expected_class"]))
            atlas_out = out.parent / "atlas"
            self.assertEqual(0, main([
                "interval-atlas", "--candidates", str(out / "reconstructed_candidate_universe.tsv"),
                "--orthology", str(out / "reconstructed_candidate_universe.tsv"),
                "--out-dir", str(atlas_out), "--case-id", "draft-qtl",
            ]))
            atlas_manifest = json.loads((atlas_out / "interval_atlas_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(4, atlas_manifest["gene_count"])
            self.assertEqual(0, atlas_manifest["unclassified_count"])
            route_b_out = out.parent / "route-b"
            self.assertEqual(0, main([
                "interval-reverse-lookup", "--atlas", str(atlas_out / "interval_gene_atlas.tsv"),
                "--out-dir", str(route_b_out), "--case-id", "draft-qtl",
            ]))
            route_b_status = json.loads((route_b_out / "interval_first_status.json").read_text(encoding="utf-8"))
            self.assertTrue(route_b_status["semantically_evaluable"])
            self.assertEqual(4, route_b_status["hypothesis_count"])

    def test_gene_order_disagreement_is_review_required_not_majority_truth(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp), discordant=True)
            manifest = json.loads((out / "comparative_reconstruction_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("rearrangement_suspected", manifest["locus_reconstruction_status"])
            self.assertGreater(sum(manifest["gene_order_inversions_by_species"].values()), 0)
            review = pd.read_csv(out / "comparative_reconstruction_review.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("rearrangement_suspected", set(review["reason"]))

    def test_multi_hit_marker_species_is_excluded_and_support_can_be_insufficient(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = self._inputs(root)
            hits = pd.read_csv(inputs["hits"], sep="\t", dtype=str)
            hits.loc[hits["species"].isin(["RefB", "RefC"]), "hit_count"] = "2"
            hits.to_csv(inputs["hits"], sep="\t", index=False, lineterminator="\n")
            out = root / "out"
            self.assertEqual(0, main([
                "comparative-locus-reconstruct", "--markers", str(inputs["markers"]),
                "--reference-panel", str(inputs["panel"]), "--marker-projections", str(inputs["hits"]),
                "--reference-locus-genes", str(inputs["genes"]), "--out-dir", str(out),
            ]))
            manifest = json.loads((out / "comparative_reconstruction_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("insufficient_reference_support", manifest["locus_reconstruction_status"])
            projections = pd.read_csv(out / "comparative_marker_projection.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(2, int((projections["projection_status"] == "multi_hit").sum()))

    def test_requires_two_marker_sequences(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = self._inputs(root)
            _write(inputs["markers"], ">only\nACGT\n")
            with self.assertRaisesRegex(ValueError, "at least two"):
                main([
                    "comparative-locus-reconstruct", "--markers", str(inputs["markers"]),
                    "--reference-panel", str(inputs["panel"]), "--marker-projections", str(inputs["hits"]),
                    "--reference-locus-genes", str(inputs["genes"]), "--out-dir", str(root / "out"),
                ])

    def test_unresolved_reference_families_do_not_form_a_false_consensus(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = self._inputs(root)
            genes = pd.read_csv(inputs["genes"], sep="\t", dtype=str).fillna("")
            genes["orthogroup"] = ""
            genes.to_csv(inputs["genes"], sep="\t", index=False, lineterminator="\n")
            out = root / "out"
            self.assertEqual(0, main([
                "comparative-locus-reconstruct", "--markers", str(inputs["markers"]),
                "--reference-panel", str(inputs["panel"]), "--marker-projections", str(inputs["hits"]),
                "--reference-locus-genes", str(inputs["genes"]), "--out-dir", str(out),
            ]))
            manifest = json.loads((out / "comparative_reconstruction_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("insufficient_reference_support", manifest["locus_reconstruction_status"])
            self.assertEqual(0, manifest["expected_orthogroup_count"])
            review = pd.read_csv(out / "comparative_reconstruction_review.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("no_resolved_orthogroup_consensus", set(review["reason"]))


class ComparativeLocusBenchmarkTests(unittest.TestCase):
    @unittest.skip("public release: bundled benchmark fixture under benchmarks/ (not shipped)")
    def test_bundled_degraded_reference_fixture(self) -> None:
        fixture = Path(__file__).resolve().parents[1] / "benchmarks" / "comparative_locus_v1"
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reconstruction = root / "reconstruction"
            self.assertEqual(0, main([
                "comparative-locus-reconstruct", "--markers", str(fixture / "markers.fa"),
                "--reference-panel", str(fixture / "reference_panel.tsv"),
                "--marker-projections", str(fixture / "marker_projections.tsv"),
                "--reference-locus-genes", str(fixture / "reference_locus_genes.tsv"),
                "--native-genes", str(fixture / "native_genes.tsv"),
                "--focal-recovery", str(fixture / "focal_recovery.tsv"),
                "--out-dir", str(reconstruction), "--case-id", "comparative-v1",
            ]))
            manifest = _write(root / "benchmark.tsv", "\n".join([
                "case_id\tdegradation_type\ttruth_orthogroup\treconstructed_universe_tsv\treconstruction_manifest_json\tnative_orthogroups_tsv\texpected_truth_set_tsv\texpected_discordance",
                f"comparative-v1\tcausal_annotation_deleted\tOG_B\t{reconstruction / 'reconstructed_candidate_universe.tsv'}\t{reconstruction / 'comparative_reconstruction_manifest.json'}\t{fixture / 'native_orthogroups.tsv'}\t{fixture / 'expected_truth_set.tsv'}\tfalse",
            ]) + "\n")
            benchmark = root / "benchmark"
            self.assertEqual(0, main(["comparative-locus-benchmark", "--manifest", str(manifest), "--out-dir", str(benchmark)]))
            result = json.loads((benchmark / "comparative_locus_benchmark_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(1.0, result["causal_orthogroup_recovery_rate"])
            self.assertEqual(1.0, result["missing_gene_recovery_rate"])
            self.assertEqual(0.0, result["mean_false_expected_inflation"])
            self.assertEqual(3.0, result["mean_reconstructed_candidate_set_size"])

    def test_degraded_reference_benchmark_scores_recovery_inflation_and_resources(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            helper = ComparativeLocusReconstructionTests()
            helper._run(root)
            _write(root / "native_groups.tsv", "orthogroup\nOG_A\n")
            _write(root / "truth_set.tsv", "orthogroup\nOG_A\nOG_B\nOG_C\nOG_D\n")
            bench_manifest = _write(root / "benchmark.tsv", "\n".join([
                "case_id\tdegradation_type\ttruth_orthogroup\treconstructed_universe_tsv\treconstruction_manifest_json\tnative_orthogroups_tsv\texpected_truth_set_tsv\texpected_discordance",
                "degraded-1\tcausal_annotation_deleted\tOG_B\tout/reconstructed_candidate_universe.tsv\tout/comparative_reconstruction_manifest.json\tnative_groups.tsv\ttruth_set.tsv\tfalse",
            ]) + "\n")
            bench = root / "bench"
            self.assertEqual(0, main(["comparative-locus-benchmark", "--manifest", str(bench_manifest), "--out-dir", str(bench)]))
            result = json.loads((bench / "comparative_locus_benchmark_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(1.0, result["causal_orthogroup_recovery_rate"])
            self.assertEqual(1.0, result["missing_gene_recovery_rate"])
            self.assertEqual(0.0, result["mean_false_expected_inflation"])
            self.assertGreater(result["by_reference_species_count"][0]["mean_peak_traced_memory_mb"], 0)
            cases = pd.read_csv(bench / "comparative_locus_benchmark_cases.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("True", cases.loc[0, "truth_absent_from_native"])
            self.assertEqual("True", cases.loc[0, "missing_gene_recovered"])

    def test_degraded_reference_benchmark_scores_expected_discordance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = ComparativeLocusReconstructionTests()._run(root, discordant=True)
            manifest = _write(root / "benchmark.tsv", "\n".join([
                "case_id\tdegradation_type\ttruth_orthogroup\treconstructed_universe_tsv\treconstruction_manifest_json\texpected_discordance",
                f"discordant\tmarker_only_coordinate_removed\tOG_B\t{out / 'reconstructed_candidate_universe.tsv'}\t{out / 'comparative_reconstruction_manifest.json'}\ttrue",
            ]) + "\n")
            bench = root / "bench"
            self.assertEqual(0, main(["comparative-locus-benchmark", "--manifest", str(manifest), "--out-dir", str(bench)]))
            result = json.loads((bench / "comparative_locus_benchmark_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(1.0, result["synteny_discordance_detection_rate"])

    def test_unknown_degradation_type_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = _write(root / "benchmark.tsv", "case_id\tdegradation_type\ttruth_orthogroup\treconstructed_universe_tsv\treconstruction_manifest_json\nx\timaginary\tOG1\tx.tsv\tx.json\n")
            with self.assertRaisesRegex(ValueError, "unknown degradation_type"):
                main(["comparative-locus-benchmark", "--manifest", str(manifest), "--out-dir", str(root / "out")])


if __name__ == "__main__":
    unittest.main()
