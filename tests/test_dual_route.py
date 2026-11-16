"""Forward-lookup and reverse-lookup discovery routes (issue #169)."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main


ATLAS = "\n".join([
    "gene_id\tchrom\tstart\tend\tfamily_id\tfamily_label\torthogroup\torthology_relationship\torthology_reference_gene\tdomain_architecture\tbroad_function_bins\tevidence_completeness",
    "G_CLV1\tchr1\t1000\t3000\tCLV\tCLV1\tOG1\tone-to-one\tCLV1\tLRR-Pkinase\tmeristem_development;signaling_receptor\thigh_confidence_family_assignment",
    "G_CLV2\tchr1\t5000\t7000\tCLV\tCLV2\tOG1\tone-to-many\tCLV2\tLRR-Pkinase\tsignaling_receptor\tfamily_assignment",
    "G_PXY\tchr1\t20000\t22000\tPXY\tPXY1\tOG2\tone-to-one\tPXY1\tLRR-Pkinase\tvascular_cell_wall;signaling_receptor\thigh_confidence_family_assignment",
    "G_WOX\tchr1\t30000\t32000\tWOX\tWOX4\tOG3\tone-to-one\tWOX4\tHomeobox\ttranscription_chromatin\thigh_confidence_family_assignment",
    "G_DARK\tchr1\t40000\t40300\t\t\t\t\t\t\t\tno_informative_annotation",
]) + "\n"

REFERENCE_FUNCTION = "\n".join([
    "reference_gene_id\tfunction_description\tpathway\tmutant_phenotype\tsource_id",
    "CLV1\treceptor kinase restricting the stem cell pool\tCLAVATA-WUSCHEL\tenlarged meristem\tPMID:1001",
    "PXY1\treceptor kinase for procambial cell division\tCLE-receptor vascular signalling\tdisorganised vasculature\tPMID:1002",
]) + "\n"


class PhenotypeFirstStatusTests(unittest.TestCase):
    def _fixture(self, root: Path, *, present_in_interval: bool) -> dict[str, Path]:
        paths = {name: root / name for name in ("hypotheses.tsv", "validation.tsv", "panel.tsv", "falsification.tsv")}
        paths["hypotheses.tsv"].write_text(
            "hypothesis_id\tpathway_or_process\tgene_family\treference_gene_name\tproposed_stable_id\treview_status\n"
            "hyp_001\tCLAVATA-WUSCHEL\tCLV\tCLV1\tP0CLV1\tpending_verification\n", encoding="utf-8")
        paths["validation.tsv"].write_text(
            "hypothesis_id\tverification_status\nhyp_001\tverified\n", encoding="utf-8")
        paths["panel.tsv"].write_text("bait_id\tcase_id\nP0CLV1\tdemo\n", encoding="utf-8")
        verdict = "surviving_enriched" if present_in_interval else "falsified_no_member_in_interval"
        genes = "G_CLV1;G_CLV2" if present_in_interval else ""
        paths["falsification.tsv"].write_text(
            "hypothesis_id\tlevel\tverdict\tinterval_gene_ids\n"
            f"hyp_001\tcoarse\t{verdict}\t{genes}\n"
            f"hyp_001\tfine\t{verdict}\t{genes}\n", encoding="utf-8")
        return paths

    def test_supported_route_a_reports_no_signal(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._fixture(root, present_in_interval=True)
            out = root / "route_a"
            self.assertEqual(0, main(["phenotype-first-status", "--out-dir", str(out), "--case-id", "demo",
                                      "--hypotheses", str(fixture["hypotheses.tsv"]),
                                      "--validation", str(fixture["validation.tsv"]),
                                      "--panel", str(fixture["panel.tsv"]),
                                      "--falsification", str(fixture["falsification.tsv"])]))
            status = json.loads((out / "phenotype_first_status.json").read_text(encoding="utf-8"))
            self.assertEqual([], status["signal_names"])
            self.assertFalse(status["route_b_recommended"])
            self.assertEqual("low", status["risk"]["hallucination_risk"])
            hypotheses = pd.read_csv(out / "phenotype_first_hypotheses.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("interval_anchored", hypotheses.iloc[0]["anchor_status"])
            self.assertEqual("phenotype_first", hypotheses.iloc[0]["route"])

    def test_absent_family_raises_the_signal_that_triggers_route_b(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._fixture(root, present_in_interval=False)
            out = root / "route_a"
            self.assertEqual(0, main(["phenotype-first-status", "--out-dir", str(out),
                                      "--hypotheses", str(fixture["hypotheses.tsv"]),
                                      "--validation", str(fixture["validation.tsv"]),
                                      "--panel", str(fixture["panel.tsv"]),
                                      "--falsification", str(fixture["falsification.tsv"])]))
            status = json.loads((out / "phenotype_first_status.json").read_text(encoding="utf-8"))
            self.assertIn("expected_family_absent_from_interval", status["signal_names"])
            self.assertIn("pathway_supported_but_interval_mismatch", status["signal_names"])
            self.assertTrue(status["route_b_recommended"])
            # A verified bait with no interval member is exactly where invention starts.
            self.assertEqual(["hyp_001"], status["unanchored_hypotheses"])
            self.assertEqual("high", status["risk"]["hallucination_risk"])

    def test_no_hypothesis_and_no_bait_are_separate_signals(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "route_a"
            self.assertEqual(0, main(["phenotype-first-status", "--out-dir", str(out)]))
            status = json.loads((out / "phenotype_first_status.json").read_text(encoding="utf-8"))
            self.assertIn("no_supported_pathway", status["signal_names"])
            self.assertIn("no_verified_bait", status["signal_names"])
            self.assertTrue(status["route_b_recommended"])

    def test_stable_process_id_passes_through_route_a(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._fixture(root, present_in_interval=True)
            fixture["hypotheses.tsv"].write_text(
                "hypothesis_id\tpathway_or_process\tprocess_id\tgene_family\treference_gene_name\tproposed_stable_id\treview_status\n"
                "hyp_001\tCLAVATA-WUSCHEL\tGO:0001234\tCLV\tCLV1\tP0CLV1\tpending_verification\n", encoding="utf-8")
            out = root / "route_a"
            self.assertEqual(0, main(["phenotype-first-status", "--out-dir", str(out),
                                      "--hypotheses", str(fixture["hypotheses.tsv"]),
                                      "--validation", str(fixture["validation.tsv"]),
                                      "--panel", str(fixture["panel.tsv"]),
                                      "--falsification", str(fixture["falsification.tsv"])]))
            hypotheses = pd.read_csv(out / "phenotype_first_hypotheses.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("GO:0001234", hypotheses.iloc[0]["process_id"])


class IntervalReverseLookupTests(unittest.TestCase):
    def _fixture(self, root: Path) -> dict[str, Path]:
        paths = {name: root / name for name in ("atlas.tsv", "functions.tsv")}
        paths["atlas.tsv"].write_text(ATLAS, encoding="utf-8")
        paths["functions.tsv"].write_text(REFERENCE_FUNCTION, encoding="utf-8")
        return paths

    def _run(self, root: Path, extra: list[str] | None = None) -> Path:
        fixture = self._fixture(root)
        out = root / "route_b"
        self.assertEqual(0, main(["interval-reverse-lookup", "--atlas", str(fixture["atlas.tsv"]),
                                  "--reference-function", str(fixture["functions.tsv"]),
                                  "--out-dir", str(out), "--case-id", "demo"] + (extra or [])))
        return out

    def test_hypotheses_come_back_from_the_interval_with_relationship_tiers(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            hypotheses = pd.read_csv(out / "family_process_hypotheses.tsv", sep="\t", dtype=str).fillna("").set_index("hypothesis_id")
            self.assertIn("rev_fam_CLV", hypotheses.index)
            self.assertIn("rev_fam_PXY", hypotheses.index)
            self.assertEqual("CLAVATA-WUSCHEL", hypotheses.loc["rev_fam_CLV", "pathway_or_process"])
            # PXY is an adjacent mechanism, kept but with its own tier, not merged into CLV.
            self.assertEqual("CLE-receptor vascular signalling", hypotheses.loc["rev_fam_PXY", "pathway_or_process"])
            self.assertEqual("exact_ortholog", hypotheses.loc["rev_fam_PXY", "relationship_tier"])
            self.assertTrue(hypotheses["rank_contribution"].eq("none").all())
            self.assertTrue(hypotheses["anchor_status"].eq("interval_anchored").all())

    def test_reference_may_be_the_focal_genome_annotation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._fixture(root)
            # No cross-species homolog table at all: annotation keyed on the focal gene.
            fixture["functions.tsv"].write_text(
                "gene_id\tfunction_description\tpathway\tsource_id\n"
                "G_WOX\thomeodomain regulator of the stem cell niche\tWUSCHEL-related\tlocal:annotation\n", encoding="utf-8")
            out = root / "route_b"
            self.assertEqual(0, main(["interval-reverse-lookup", "--atlas", str(fixture["atlas.tsv"]),
                                      "--reference-function", str(fixture["functions.tsv"]), "--out-dir", str(out)]))
            lookup = pd.read_csv(out / "interval_reverse_lookup.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            self.assertEqual("focal_reference_annotation", lookup.loc["G_WOX", "annotation_origin"])
            self.assertEqual("local:annotation", lookup.loc["G_WOX", "reference_source_id"])
            hypotheses = pd.read_csv(out / "family_process_hypotheses.tsv", sep="\t", dtype=str).fillna("").set_index("hypothesis_id")
            self.assertEqual("WUSCHEL-related", hypotheses.loc["rev_fam_WOX", "pathway_or_process"])
            self.assertEqual("source_backed", hypotheses.loc["rev_fam_WOX", "verification_status"])

    def test_stable_pathway_id_is_carried_into_route_b_hypotheses(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._fixture(root)
            fixture["functions.tsv"].write_text(
                "reference_gene_id\tfunction_description\tpathway\tpathway_id\tmutant_phenotype\tsource_id\n"
                "CLV1\treceptor kinase restricting the stem cell pool\tCLAVATA-WUSCHEL\tGO:0001234\tenlarged meristem\tPMID:1001\n"
                "PXY1\treceptor kinase for procambial cell division\tCLE-receptor vascular signalling\tGO:0005678\tdisorganised vasculature\tPMID:1002\n",
                encoding="utf-8")
            out = root / "route_b"
            self.assertEqual(0, main(["interval-reverse-lookup", "--atlas", str(fixture["atlas.tsv"]),
                                      "--reference-function", str(fixture["functions.tsv"]), "--out-dir", str(out)]))
            lookup = pd.read_csv(out / "interval_reverse_lookup.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            self.assertEqual("GO:0001234", lookup.loc["G_CLV1", "reference_pathway_id"])
            hypotheses = pd.read_csv(out / "family_process_hypotheses.tsv", sep="\t", dtype=str).fillna("").set_index("hypothesis_id")
            self.assertEqual("GO:0001234", hypotheses.loc["rev_fam_CLV", "process_id"])
            self.assertEqual("GO:0005678", hypotheses.loc["rev_fam_PXY", "process_id"])

    def test_unclassified_genes_become_a_probe_and_are_never_dropped(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            lookup = pd.read_csv(out / "interval_reverse_lookup.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(5, len(lookup))
            self.assertEqual("unclassified", lookup.set_index("gene_id").loc["G_DARK", "lookup_status"])
            hypotheses = pd.read_csv(out / "family_process_hypotheses.tsv", sep="\t", dtype=str).fillna("").set_index("hypothesis_id")
            self.assertIn("G_DARK", hypotheses.loc["rev_unknown_probe", "member_interval_genes"])
            reasons = set(pd.read_csv(out / "reverse_lookup_review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("unclassified_interval_genes_carried_forward", reasons)

    def test_unanchored_llm_proposal_is_rejected_as_a_bait(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._fixture(root)
            proposals = root / "llm.tsv"
            proposals.write_text(
                "hypothesis_id\tpathway_or_process\tgene_family\treference_gene_name\tproposed_stable_id\tcounterevidence\n"
                "hyp_001\tCLAVATA-WUSCHEL\tCLV\tCLV1\tP0CLV1\t\n"
                "hyp_002\tflowering time\tFT\tFT\tP0FT001\t\n", encoding="utf-8")
            out = root / "route_b"
            self.assertEqual(0, main(["interval-reverse-lookup", "--atlas", str(fixture["atlas.tsv"]),
                                      "--reference-function", str(fixture["functions.tsv"]),
                                      "--llm-proposals", str(proposals), "--out-dir", str(out)]))
            hypotheses = pd.read_csv(out / "family_process_hypotheses.tsv", sep="\t", dtype=str).fillna("").set_index("hypothesis_id")
            # CLV is in the interval listing, FT is not.
            self.assertEqual("interval_anchored", hypotheses.loc["rev_llm_hyp_001", "anchor_status"])
            self.assertEqual("unanchored", hypotheses.loc["rev_llm_hyp_002", "anchor_status"])
            self.assertEqual("rejected_unanchored", hypotheses.loc["rev_llm_hyp_002", "verification_status"])
            reasons = set(pd.read_csv(out / "reverse_lookup_review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("llm_proposal_not_anchored_in_interval", reasons)
            status = json.loads((out / "interval_first_status.json").read_text(encoding="utf-8"))
            self.assertEqual(1, status["unanchored_proposals"])

    def test_prompt_lists_only_interval_genes_and_forbids_invention(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp), extra=["--phenotype", "stem fasciation"])
            prompt = (out / "reverse_lookup_prompt.txt").read_text(encoding="utf-8")
            self.assertIn("Do not introduce a gene, family or accession that is absent", prompt)
            self.assertIn("is a valid and useful answer", prompt)
            listing = prompt.split("present in this QTL interval")[1].split("Return JSON only")[0]
            self.assertIn("CLV1", listing)
            self.assertIn("PXY1", listing)
            self.assertIn("(no annotated identity)", listing)
            # Only interval entries are shown; the focal gene ids are not sent.
            self.assertNotIn("G_CLV1", listing)
            record = json.loads((out / "reverse_lookup_prompts.jsonl").read_text(encoding="utf-8").strip())
            self.assertEqual({"prompt_id", "gene_id", "system", "user"}, set(record))

    def test_status_marks_the_lane_as_interval_informed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp), extra=["--phenotype", "stem fasciation"])
            status = json.loads((out / "interval_first_status.json").read_text(encoding="utf-8"))
            self.assertEqual("interval_informed", status["lane"])
            self.assertTrue(status["interval_information_in_model_input"])
            self.assertFalse(status["blind_discovery_claim"])
            self.assertTrue(status["phenotype_used"])
            self.assertTrue(status["command_completed"])
            self.assertTrue(status["semantically_evaluable"])
            self.assertEqual("evaluable", status["semantic_status"])
            self.assertEqual("b1_phenotype_overlay", status["lane_variant"])
            self.assertGreater(status["provider_coverage"]["semantic_coverage_fraction"], 0)
            self.assertIn("family_process_hypotheses", status["frozen_output_fingerprints"])
            self.assertIn("post-hoc rationalisation", status["risk"]["failure_mode"])
            self.assertIn("reference_homolog", status["annotation_origin_counts"])
            self.assertEqual("available_source_backed", status["process_inference_status"])
            self.assertTrue(status["phenotype_overlay"]["non_scoring"])

    def test_b1_overlay_is_auditable_and_changes_no_tier_or_score(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "b0").mkdir()
            (root / "b1").mkdir()
            b0 = self._run(root / "b0")
            b1 = self._run(root / "b1", extra=["--phenotype", "enlarged meristem stem pool"])
            lookup0 = pd.read_csv(b0 / "interval_reverse_lookup.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            lookup1 = pd.read_csv(b1 / "interval_reverse_lookup.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            self.assertEqual("not_applied", lookup0.loc["G_CLV1", "phenotype_overlay_status"])
            self.assertEqual("phenotype_consistent", lookup1.loc["G_CLV1", "phenotype_overlay_status"])
            self.assertIn("meristem", lookup1.loc["G_CLV1", "phenotype_overlap_terms"])
            self.assertEqual(lookup0["relationship_tier"].to_dict(), lookup1["relationship_tier"].to_dict())
            hypotheses = pd.read_csv(b1 / "family_process_hypotheses.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(hypotheses["rank_contribution"].eq("none").all())

    def test_shared_route_overlay_recovers_semantic_terms_without_hiding_raw_false_positive(self) -> None:
        """Issue #332: B1 keeps raw lexical and shared-route lexical overlays separate.

        These pea_r / pea_i / pea_fa analogues reproduce the observed lexical
        gaps.  In pea_i, raw ``green`` matches the ferredoxin decoy but not the
        true ``stay-green`` annotation; only the shared route's chlorophyll /
        stay-green vocabulary marks the truth.  This is diagnostic only.
        """
        scenarios = (
            ("pea_r", "wrinkled seed", "G_R_TRUE", "SBEI", "SBEI",
             "starch branching enzyme", "starch biosynthesis", "wrinkled",
             "seed_starch_metabolism", "starch|glucan"),
            ("pea_i", "green cotyledon", "G_I_TRUE", "SGR", "SGR",
             "stay-green regulator of chlorophyll degradation", "chlorophyll degradation", "green",
             "chlorophyll_degradation", "chlorophyll|stay-green|degreening"),
            ("pea_fa", "fasciated stem Mendel trait", "G_FA_TRUE", "CIK", "CIK1",
             "CLAVATA receptor kinase regulating meristem size", "meristem maintenance", "fasciated",
             "meristem_maintenance", "meristem|CLAVATA|receptor"),
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route = root / "shared_route.tsv"
            route.write_text("\n".join([
                "phenotype_token\tpathway_family\tannotation_keywords",
                *[f"{row[7]}\t{row[8]}\t{row[9]}" for row in scenarios],
            ]) + "\n", encoding="utf-8")
            for case_id, phenotype, truth_gene, family, reference, description, pathway, _token, expected_family, _keywords in scenarios:
                atlas_rows = [
                    "gene_id\tfamily_id\tfamily_label\torthogroup\torthology_relationship\torthology_reference_gene",
                    f"{truth_gene}\t{family}\t{reference}\tOG_{family}\tone-to-one\t{reference}",
                ]
                function_rows = [
                    "reference_gene_id\tfunction_description\tpathway\tmutant_phenotype\tsource_id",
                    f"{reference}\t{description}\t{pathway}\t\tPMID:{case_id}",
                ]
                if case_id == "pea_i":
                    atlas_rows.append("G_I_GREEN_DECOY\tFDX\tFDX\tOG_FDX\tone-to-one\tFDX")
                    function_rows.append("FDX\tpale green leaves ferredoxin\tphotosynthesis\t\tPMID:decoy")
                atlas = root / f"{case_id}_atlas.tsv"
                atlas.write_text("\n".join(atlas_rows) + "\n", encoding="utf-8")
                functions = root / f"{case_id}_functions.tsv"
                functions.write_text("\n".join(function_rows) + "\n", encoding="utf-8")
                out = root / case_id
                self.assertEqual(0, main([
                    "interval-reverse-lookup", "--atlas", str(atlas), "--reference-function", str(functions),
                    "--phenotype", phenotype, "--shared-route", str(route), "--out-dir", str(out), "--case-id", case_id,
                ]))
                lookup = pd.read_csv(out / "interval_reverse_lookup.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
                self.assertEqual("phenotype_not_matched", lookup.loc[truth_gene, "phenotype_overlay_status"])
                self.assertEqual("phenotype_consistent", lookup.loc[truth_gene, "shared_route_overlay_status"])
                hypotheses = pd.read_csv(out / "family_process_hypotheses.tsv", sep="\t", dtype=str).fillna("")
                truth_hypothesis = hypotheses.set_index("gene_family").loc[family]
                self.assertEqual("phenotype_not_matched", truth_hypothesis["phenotype_overlay_status"])
                self.assertEqual("phenotype_consistent", truth_hypothesis["shared_route_overlay_status"])
                self.assertTrue(hypotheses["rank_contribution"].eq("none").all())
                status = json.loads((out / "interval_first_status.json").read_text(encoding="utf-8"))
                shared = status["phenotype_overlay"]["shared_route"]
                self.assertEqual("applied", shared["status"])
                self.assertIn(expected_family, shared["pathway_families"])
                self.assertIn("shared_route_overlay_status", lookup.columns)
                vocabulary = pd.read_csv(out / "phenotype_overlay_vocabulary.tsv", sep="\t", dtype=str).fillna("")
                self.assertIn("raw_phenotype_token", set(vocabulary["vocabulary_source"]))
                self.assertIn("shared_route_pathway_family", set(vocabulary["vocabulary_source"]))
                self.assertIn("shared_route_annotation_keyword", set(vocabulary["vocabulary_source"]))
                if case_id == "pea_i":
                    self.assertEqual("phenotype_consistent", lookup.loc["G_I_GREEN_DECOY", "phenotype_overlay_status"])
                    self.assertEqual("phenotype_not_matched", lookup.loc["G_I_GREEN_DECOY", "shared_route_overlay_status"])

    def test_missing_reference_function_has_explicit_process_status(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._fixture(root)
            out = root / "route_b"
            self.assertEqual(0, main(["interval-reverse-lookup", "--atlas", str(fixture["atlas.tsv"]), "--out-dir", str(out)]))
            status = json.loads((out / "interval_first_status.json").read_text(encoding="utf-8"))
            self.assertEqual("unavailable_no_reference_function", status["process_inference_status"])
            self.assertGreater(status["relationship_tier_resolution"]["orthogroup_gene_count"], 0)

    def test_all_unclassified_atlas_completes_but_is_not_semantically_evaluable(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atlas = root / "minimal.tsv"
            atlas.write_text(
                "gene_id\tfamily_id\torthogroup\torthology_reference_gene\n"
                "GeneP\t\t\t\nGeneQ\t\t\t\nGeneR\t\t\t\n",
                encoding="utf-8",
            )
            out = root / "route_b"
            self.assertEqual(0, main([
                "interval-reverse-lookup", "--atlas", str(atlas),
                "--out-dir", str(out), "--case-id", "minimal",
            ]))
            status = json.loads((out / "interval_first_status.json").read_text(encoding="utf-8"))
            self.assertTrue(status["command_completed"])
            self.assertFalse(status["semantically_evaluable"])
            self.assertEqual("not_evaluated_provider_missing", status["semantic_status"])
            self.assertEqual(0, status["provider_coverage"]["semantically_annotated_gene_count"])
            self.assertEqual(0.0, status["provider_coverage"]["semantic_coverage_fraction"])
            hypotheses = pd.read_csv(out / "family_process_hypotheses.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(["unknown_probe"], hypotheses["scope"].tolist())
            reviews = pd.read_csv(out / "reverse_lookup_review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("route_b_semantic_inputs_insufficient", set(reviews["reason"]))

    def test_minimum_semantic_coverage_can_hold_a_partially_annotated_lane(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp), extra=["--min-semantic-coverage", "0.9"])
            status = json.loads((out / "interval_first_status.json").read_text(encoding="utf-8"))
            self.assertTrue(status["command_completed"])
            self.assertFalse(status["semantically_evaluable"])
            self.assertEqual(0.8, status["provider_coverage"]["semantic_coverage_fraction"])
            self.assertEqual(0.9, status["provider_coverage"]["min_semantic_coverage"])

    def test_orthogroup_only_atlas_produces_an_anchored_hypothesis(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atlas = root / "orthogroups.tsv"
            atlas.write_text(
                "gene_id\torthogroup\nGeneA\tOG_A\nGeneB\tOG_B\n",
                encoding="utf-8",
            )
            out = root / "route_b"
            self.assertEqual(0, main([
                "interval-reverse-lookup", "--atlas", str(atlas),
                "--out-dir", str(out), "--case-id", "orthogroup-only",
            ]))
            status = json.loads((out / "interval_first_status.json").read_text(encoding="utf-8"))
            self.assertTrue(status["semantically_evaluable"])
            hypotheses = pd.read_csv(out / "family_process_hypotheses.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"OG_A", "OG_B"}, set(hypotheses["gene_family"]))
            self.assertEqual({"orthogroup"}, set(hypotheses["scope"]))
            self.assertEqual({"interval_anchored"}, set(hypotheses["anchor_status"]))


class SemanticQualityTests(unittest.TestCase):
    """Issue #305: family-label quality section of interval_first_status.json."""

    def _run_atlas(self, root: Path, atlas_text: str) -> dict:
        atlas = root / "atlas.tsv"
        atlas.write_text(atlas_text, encoding="utf-8")
        out = root / "route_b"
        self.assertEqual(0, main([
            "interval-reverse-lookup", "--atlas", str(atlas),
            "--out-dir", str(out), "--case-id", "quality",
        ]))
        return json.loads((out / "interval_first_status.json").read_text(encoding="utf-8"))

    def test_symbol_family_labels_are_healthy(self) -> None:
        # The standard fixture carries CLV1 / CLV2 / PXY1 / WOX4 symbols; the
        # dark gene has no family label and counts as unresolved.
        with tempfile.TemporaryDirectory() as tmp:
            out = IntervalReverseLookupTests()._run(Path(tmp))
            status = json.loads((out / "interval_first_status.json").read_text(encoding="utf-8"))
            quality = status["semantic_quality"]
            self.assertEqual("healthy", quality["family_label_quality"])
            self.assertEqual(1.0, quality["usable_symbol_fraction"])
            self.assertEqual(0.0, quality["accession_or_locus_tag_fraction"])
            self.assertEqual(0.2, quality["unresolved_family_fraction"])
            self.assertEqual("1.4", status["schema_version"])
            self.assertIn("family_label", quality["rule"])

    def test_accession_only_family_labels_are_degraded(self) -> None:
        # A legume TrEMBL-style atlas whose family labels are all UniProt
        # accessions still passes the semantic coverage gate (family_id is
        # non-empty) but its quality is degraded: no label is a usable symbol.
        with tempfile.TemporaryDirectory() as tmp:
            status = self._run_atlas(Path(tmp),
                "gene_id\tfamily_id\tfamily_label\torthogroup\torthology_reference_gene\n"
                "G_A\tA0A072VKZ1\tA0A072VKZ1\tOG_A\tA0A072VKZ1\n"
                "G_B\tA0A072VKZ1\tA0A072VKZ1\tOG_A\tA0A072VKZ1\n")
            self.assertTrue(status["semantically_evaluable"])
            self.assertEqual("evaluable", status["semantic_status"])
            quality = status["semantic_quality"]
            self.assertEqual("degraded", quality["family_label_quality"])
            self.assertEqual(0.0, quality["usable_symbol_fraction"])
            self.assertEqual(1.0, quality["accession_or_locus_tag_fraction"])
            self.assertEqual(1.0, quality["unresolved_family_fraction"])

    def test_partial_usable_fraction_below_half_is_degraded(self) -> None:
        # 1 usable symbol out of 3 family-annotated genes: 0.333 < 0.5.
        with tempfile.TemporaryDirectory() as tmp:
            status = self._run_atlas(Path(tmp),
                "gene_id\tfamily_id\tfamily_label\torthogroup\n"
                "G_A\tCLV\tCLV1\tOG_A\n"
                "G_B\tA0A072VKZ1\tA0A072VKZ1\tOG_A\n"
                "G_C\tA0A072VKZ1\tA0A072VKZ1\tOG_A\n")
            self.assertTrue(status["semantically_evaluable"])
            quality = status["semantic_quality"]
            self.assertEqual("degraded", quality["family_label_quality"])
            self.assertEqual(0.333333, quality["usable_symbol_fraction"])
            self.assertEqual(0.666667, quality["accession_or_locus_tag_fraction"])
            self.assertEqual(0.666667, quality["unresolved_family_fraction"])

    def test_no_family_labels_is_not_evaluable_but_still_semantically_evaluable(self) -> None:
        # An orthogroup-only atlas has semantic coverage but no family_label
        # axis at all: quality is not_evaluable (consistent with provider
        # missing), while semantically_evaluable stays True.
        with tempfile.TemporaryDirectory() as tmp:
            status = self._run_atlas(Path(tmp),
                "gene_id\torthogroup\nGeneA\tOG_A\nGeneB\tOG_B\n")
            self.assertTrue(status["semantically_evaluable"])
            quality = status["semantic_quality"]
            self.assertEqual("not_evaluable", quality["family_label_quality"])
            self.assertEqual(0.0, quality["usable_symbol_fraction"])

    def test_usable_fraction_exactly_half_is_healthy(self) -> None:
        # 0.5 is the frozen boundary: >= 0.5 is healthy, < 0.5 is degraded.
        with tempfile.TemporaryDirectory() as tmp:
            status = self._run_atlas(Path(tmp),
                "gene_id\tfamily_id\tfamily_label\torthogroup\n"
                "G_A\tCLV\tCLV1\tOG_A\n"
                "G_B\tA0A072VKZ1\tA0A072VKZ1\tOG_A\n")
            self.assertTrue(status["semantically_evaluable"])
            quality = status["semantic_quality"]
            self.assertEqual("healthy", quality["family_label_quality"])
            self.assertEqual(0.5, quality["usable_symbol_fraction"])
            self.assertEqual(0.5, quality["accession_or_locus_tag_fraction"])


class DualRouteCompareTests(unittest.TestCase):
    ROUTE_COLUMNS = ["hypothesis_id", "route", "pathway_or_process", "process_id", "gene_family",
                     "member_interval_genes", "anchor_status", "relationship_tier", "source_backing"]

    def _routes(self, root: Path, *, route_a_family: str) -> tuple[Path, Path]:
        anchored = route_a_family == "CLV"
        pathway = "CLAVATA-WUSCHEL" if anchored else "flowering time"
        route_a = root / "route_a.tsv"
        route_a.write_text(
            "hypothesis_id\troute\tpathway_or_process\tgene_family\tmember_interval_genes\tanchor_status\trelationship_tier\tsource_backing\n"
            f"hyp_001\tphenotype_first\t{pathway}\t{route_a_family}\t"
            f"{'G_CLV1;G_CLV2' if anchored else ''}\t"
            f"{'interval_anchored' if anchored else 'unanchored'}\texact_ortholog\tverified_bait\n",
            encoding="utf-8")
        fixture = IntervalReverseLookupTests()._fixture(root)
        out_b = root / "route_b"
        self.assertEqual(0, main(["interval-reverse-lookup", "--atlas", str(fixture["atlas.tsv"]),
                                  "--reference-function", str(fixture["functions.tsv"]), "--out-dir", str(out_b)]))
        return route_a, out_b / "family_process_hypotheses.tsv"

    def _write_route(self, path: Path, rows: list[dict[str, str]]) -> None:
        frame = pd.DataFrame(rows, columns=self.ROUTE_COLUMNS)
        frame.to_csv(path, sep="\t", index=False, lineterminator="\n")

    def _compare(self, root: Path, route_a: Path, route_b: Path, extra: list[str] | None = None) -> Path:
        out = root / "compare"
        self.assertEqual(0, main(["dual-route-compare", "--route-a", str(route_a), "--route-b", str(route_b),
                                  "--out-dir", str(out), "--case-id", "demo"] + (extra or [])))
        return out

    def test_shared_interval_genes_mechanistic_pair_is_strictly_convergent(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a, route_b = self._routes(root, route_a_family="CLV")
            out = self._compare(root, route_a, route_b)
            convergence = json.loads((out / "dual_route_convergence.json").read_text(encoding="utf-8"))
            # Both routes name the same curated concept (CLAVATA-WUSCHEL) and the
            # same family with source backing on both sides -> strict tier 1.
            self.assertEqual("strict_mechanistic_convergent", convergence["overall_status"])
            comparison = pd.read_csv(out / "dual_route_hypothesis_comparison.tsv", sep="\t", dtype=str).fillna("")
            mechanistic = comparison[comparison.convergence_status == "mechanistic_convergence"]
            self.assertTrue(len(mechanistic))
            self.assertEqual("tier1_mechanistic", mechanistic.iloc[0]["match_tier"])
            self.assertIn("shared_interval_gene", mechanistic.iloc[0]["match_basis"])
            self.assertIn("same_process_concept", mechanistic.iloc[0]["match_basis"])
            self.assertEqual("shared_interval_genes", mechanistic.iloc[0]["interval_overlap_status"])
            # Route B reached PXY on its own; it stays visible as a process_related
            # match (same signaling architecture via the curated vocabulary, tier 4),
            # never merged into the CLV mechanism.
            process_related = comparison[comparison.convergence_status == "process_related"]
            self.assertIn("rev_fam_PXY", set(process_related["route_b_hypothesis_id"]))
            self.assertEqual("tier4_process_related", process_related.iloc[0]["match_tier"])
            interval_only = comparison[comparison.convergence_status == "interval_only_supported"]
            self.assertIn("rev_fam_WOX", set(interval_only["route_b_hypothesis_id"]))
            self.assertTrue(comparison["rank_contribution"].eq("none").all())

    def test_mirror_image_risks_are_reported_for_both_routes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a, route_b = self._routes(root, route_a_family="FT")
            out = self._compare(root, route_a, route_b)
            convergence = json.loads((out / "dual_route_convergence.json").read_text(encoding="utf-8"))
            self.assertEqual("discordant", convergence["overall_status"])
            risk = convergence["route_risk"]
            self.assertEqual("hallucinating a gene that is not in the interval", risk["route_a_failure_mode"])
            self.assertEqual("rationalising a gene that is in the interval", risk["route_b_failure_mode"])
            self.assertEqual(1, risk["route_a_unanchored"])
            self.assertEqual(0, risk["route_b_unanchored"])
            reasons = set(pd.read_csv(out / "dual_route_review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("routes_discordant", reasons)
            self.assertIn("unanchored_hypotheses_present", reasons)

    def test_both_weak_routes_back_out_instead_of_forcing_a_candidate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            empty_a = root / "a.tsv"
            empty_a.write_text("hypothesis_id\tgene_family\tmember_interval_genes\tsource_backing\n", encoding="utf-8")
            empty_b = root / "b.tsv"
            empty_b.write_text("hypothesis_id\tgene_family\tmember_interval_genes\tsource_backing\n", encoding="utf-8")
            out = self._compare(root, empty_a, empty_b)
            convergence = json.loads((out / "dual_route_convergence.json").read_text(encoding="utf-8"))
            self.assertEqual("both_weak", convergence["overall_status"])
            reasons = set(pd.read_csv(out / "dual_route_review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("both_routes_weak", reasons)

    def test_comparison_never_pools_the_lanes_for_benchmarking(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a, route_b = self._routes(root, route_a_family="CLV")
            out = self._compare(root, route_a, route_b)
            convergence = json.loads((out / "dual_route_convergence.json").read_text(encoding="utf-8"))
            self.assertIn("never pooled with blind phenotype-only recall", convergence["benchmark_rule"])
            self.assertIn("adds no score", " ".join(convergence["safety"]))

    # --- issue #186: convergence decomposition ---------------------------------

    def test_same_gene_different_mechanism_is_gene_overlap_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a = root / "a.tsv"
            self._write_route(route_a, [{
                "hypothesis_id": "hyp_001", "route": "phenotype_first",
                "pathway_or_process": "hormone signalling", "process_id": "",
                "gene_family": "CLV", "member_interval_genes": "G_CLV1;G_CLV2",
                "anchor_status": "interval_anchored", "relationship_tier": "exact_ortholog",
                "source_backing": "verified_bait"}])
            route_b = root / "b.tsv"
            self._write_route(route_b, [{
                "hypothesis_id": "rev_001", "route": "interval_first",
                "pathway_or_process": "cell-wall development", "process_id": "",
                "gene_family": "", "member_interval_genes": "G_CLV1;G_CLV2",
                "anchor_status": "interval_anchored", "relationship_tier": "same_broad_process",
                "source_backing": "PMID:2002"}])
            out = self._compare(root, route_a, route_b)
            convergence = json.loads((out / "dual_route_convergence.json").read_text(encoding="utf-8"))
            self.assertEqual("gene_overlap_only", convergence["overall_status"])
            comparison = pd.read_csv(out / "dual_route_hypothesis_comparison.tsv", sep="\t", dtype=str).fillna("")
            row = comparison.iloc[0]
            self.assertEqual("gene_overlap_only", row["convergence_status"])
            self.assertEqual("tier3_gene_overlap", row["match_tier"])
            self.assertIn("shared_interval_gene", row["match_basis"])
            self.assertNotIn("same_process_concept", row["match_basis"])
            self.assertEqual("shared_interval_genes", row["interval_overlap_status"])
            reasons = set(pd.read_csv(out / "dual_route_review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("weak_convergence_only", reasons)

    def test_same_family_different_genes_is_family_convergence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a = root / "a.tsv"
            self._write_route(route_a, [{
                "hypothesis_id": "hyp_001", "route": "phenotype_first",
                "pathway_or_process": "CLAVATA-WUSCHEL", "process_id": "",
                "gene_family": "CLV", "member_interval_genes": "G_CLV1",
                "anchor_status": "interval_anchored", "relationship_tier": "exact_ortholog",
                "source_backing": "verified_bait"}])
            route_b = root / "b.tsv"
            self._write_route(route_b, [{
                "hypothesis_id": "rev_001", "route": "interval_first",
                "pathway_or_process": "cell-wall regulation", "process_id": "",
                "gene_family": "CLV", "member_interval_genes": "G_CLV2",
                "anchor_status": "interval_anchored", "relationship_tier": "same_subfamily",
                "source_backing": "PMID:1001"}])
            out = self._compare(root, route_a, route_b)
            convergence = json.loads((out / "dual_route_convergence.json").read_text(encoding="utf-8"))
            self.assertEqual("family_convergent_only", convergence["overall_status"])
            comparison = pd.read_csv(out / "dual_route_hypothesis_comparison.tsv", sep="\t", dtype=str).fillna("")
            row = comparison.iloc[0]
            self.assertEqual("family_or_subfamily_convergence", row["convergence_status"])
            self.assertEqual("tier2_family", row["match_tier"])
            self.assertIn("same_gene_family", row["match_basis"])
            self.assertEqual("no_shared_interval_genes", row["interval_overlap_status"])

    def test_descriptive_and_symbol_family_labels_converge(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a = root / "a.tsv"
            route_b = root / "b.tsv"
            self._write_route(route_a, [{
                "hypothesis_id": "hyp_phy", "route": "phenotype_first",
                "pathway_or_process": "light response", "process_id": "",
                "gene_family": "Phytochrome (PHY)", "member_interval_genes": "G_PHY",
                "anchor_status": "interval_anchored", "relationship_tier": "same_family",
                "source_backing": "PMID:1",
            }])
            self._write_route(route_b, [{
                "hypothesis_id": "rev_phy", "route": "interval_first",
                "pathway_or_process": "photomorphogenesis", "process_id": "",
                "gene_family": "PHYA", "member_interval_genes": "G_OTHER",
                "anchor_status": "interval_anchored", "relationship_tier": "same_family",
                "source_backing": "UniProt:P14712",
            }])
            out = self._compare(root, route_a, route_b)
            comparison = pd.read_csv(out / "dual_route_hypothesis_comparison.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("family_or_subfamily_convergence", comparison.iloc[0]["convergence_status"])
            self.assertIn("same_gene_family", comparison.iloc[0]["match_basis"])

    def test_same_mechanism_different_wording_is_mechanistic(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a = root / "a.tsv"
            self._write_route(route_a, [{
                "hypothesis_id": "hyp_001", "route": "phenotype_first",
                "pathway_or_process": "stem-cell maintenance via CLV signaling", "process_id": "GO:0001234",
                "gene_family": "CLV", "member_interval_genes": "G_CLV1;G_CLV2",
                "anchor_status": "interval_anchored", "relationship_tier": "exact_ortholog",
                "source_backing": "verified_bait"}])
            route_b = root / "b.tsv"
            self._write_route(route_b, [{
                "hypothesis_id": "rev_001", "route": "interval_first",
                "pathway_or_process": "CLAVATA-WUSCHEL pathway", "process_id": "GO:0001234",
                "gene_family": "CLV", "member_interval_genes": "G_CLV1;G_CLV2",
                "anchor_status": "interval_anchored", "relationship_tier": "exact_ortholog",
                "source_backing": "PMID:1001"}])
            out = self._compare(root, route_a, route_b)
            convergence = json.loads((out / "dual_route_convergence.json").read_text(encoding="utf-8"))
            self.assertEqual("strict_mechanistic_convergent", convergence["overall_status"])
            comparison = pd.read_csv(out / "dual_route_hypothesis_comparison.tsv", sep="\t", dtype=str).fillna("")
            row = comparison.iloc[0]
            self.assertEqual("mechanistic_convergence", row["convergence_status"])
            self.assertEqual("tier1_mechanistic", row["match_tier"])
            self.assertIn("same_process_concept", row["match_basis"])

    def test_curated_vocabulary_aliases_count_as_same_concept(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a = root / "a.tsv"
            self._write_route(route_a, [{
                "hypothesis_id": "hyp_001", "route": "phenotype_first",
                "pathway_or_process": "CLAVATA-WUSCHEL", "process_id": "",
                "gene_family": "CLV", "member_interval_genes": "G_CLV1",
                "anchor_status": "interval_anchored", "relationship_tier": "exact_ortholog",
                "source_backing": "verified_bait"}])
            route_b = root / "b.tsv"
            self._write_route(route_b, [{
                "hypothesis_id": "rev_001", "route": "interval_first",
                "pathway_or_process": "CLV-WUS pathway", "process_id": "",
                "gene_family": "CLV", "member_interval_genes": "G_CLV1",
                "anchor_status": "interval_anchored", "relationship_tier": "same_family",
                "source_backing": "PMID:1001"}])
            out = self._compare(root, route_a, route_b)
            convergence = json.loads((out / "dual_route_convergence.json").read_text(encoding="utf-8"))
            self.assertEqual("strict_mechanistic_convergent", convergence["overall_status"])
            comparison = pd.read_csv(out / "dual_route_hypothesis_comparison.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("mechanistic_convergence", comparison.iloc[0]["convergence_status"])

    def test_known_stable_id_resolves_to_the_same_concept_as_vocabulary_text(self) -> None:
        # Route A carries only the free text, Route B carries the stable ID for
        # the same curated concept: they must still match mechanistically
        # (canonical pipeline: A from the LLM table, B from the reference table).
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a = root / "a.tsv"
            self._write_route(route_a, [{
                "hypothesis_id": "hyp_001", "route": "phenotype_first",
                "pathway_or_process": "CLAVATA-WUSCHEL", "process_id": "",
                "gene_family": "CLV", "member_interval_genes": "G_CLV1;G_CLV2",
                "anchor_status": "interval_anchored", "relationship_tier": "exact_ortholog",
                "source_backing": "verified_bait"}])
            route_b = root / "b.tsv"
            self._write_route(route_b, [{
                "hypothesis_id": "rev_001", "route": "interval_first",
                "pathway_or_process": "CLAVATA-WUSCHEL", "process_id": "GO:0001234",
                "gene_family": "CLV", "member_interval_genes": "G_CLV1;G_CLV2",
                "anchor_status": "interval_anchored", "relationship_tier": "exact_ortholog",
                "source_backing": "PMID:1001"}])
            out = self._compare(root, route_a, route_b)
            convergence = json.loads((out / "dual_route_convergence.json").read_text(encoding="utf-8"))
            self.assertEqual("strict_mechanistic_convergent", convergence["overall_status"])
            comparison = pd.read_csv(out / "dual_route_hypothesis_comparison.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("mechanistic_convergence", comparison.iloc[0]["convergence_status"])
            self.assertIn("same_process_concept", comparison.iloc[0]["match_basis"])

    def test_same_text_unrelated_ids_is_text_only_match(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a = root / "a.tsv"
            self._write_route(route_a, [{
                "hypothesis_id": "hyp_001", "route": "phenotype_first",
                "pathway_or_process": "novel flowering regulator cascade", "process_id": "GO:1111111",
                "gene_family": "", "member_interval_genes": "G_CLV1",
                "anchor_status": "interval_anchored", "relationship_tier": "same_broad_process",
                "source_backing": "verified_bait"}])
            route_b = root / "b.tsv"
            self._write_route(route_b, [{
                "hypothesis_id": "rev_001", "route": "interval_first",
                "pathway_or_process": "novel flowering regulator cascade", "process_id": "GO:2222222",
                "gene_family": "", "member_interval_genes": "G_PXY",
                "anchor_status": "interval_anchored", "relationship_tier": "same_broad_process",
                "source_backing": "PMID:3003"}])
            out = self._compare(root, route_a, route_b)
            convergence = json.loads((out / "dual_route_convergence.json").read_text(encoding="utf-8"))
            self.assertEqual("text_only_match", convergence["overall_status"])
            comparison = pd.read_csv(out / "dual_route_hypothesis_comparison.tsv", sep="\t", dtype=str).fillna("")
            row = comparison.iloc[0]
            self.assertEqual("text_only_match", row["convergence_status"])
            self.assertEqual("tier5_text_only", row["match_tier"])
            self.assertIn("same_pathway_text", row["match_basis"])
            self.assertNotIn("same_process_concept", row["match_basis"])

    def test_shared_source_chain_is_recorded(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a = root / "a.tsv"
            self._write_route(route_a, [{
                "hypothesis_id": "hyp_001", "route": "phenotype_first",
                "pathway_or_process": "CLAVATA-WUSCHEL", "process_id": "",
                "gene_family": "CLV", "member_interval_genes": "G_CLV1",
                "anchor_status": "interval_anchored", "relationship_tier": "exact_ortholog",
                "source_backing": "UniProt:Q8GRU6"}])
            route_b = root / "b.tsv"
            self._write_route(route_b, [{
                "hypothesis_id": "rev_001", "route": "interval_first",
                "pathway_or_process": "CLAVATA-WUSCHEL", "process_id": "",
                "gene_family": "CLV", "member_interval_genes": "G_CLV1",
                "anchor_status": "interval_anchored", "relationship_tier": "exact_ortholog",
                "source_backing": "UniProt:Q8GRU6"}])
            out = self._compare(root, route_a, route_b)
            convergence = json.loads((out / "dual_route_convergence.json").read_text(encoding="utf-8"))
            # The convergence classification itself is kept (mechanistic).
            self.assertEqual("strict_mechanistic_convergent", convergence["overall_status"])
            independence = convergence["independence"]
            self.assertEqual("shared_source_chain", independence["source_independence_status"])
            self.assertEqual(["UniProt:Q8GRU6"], independence["shared_source_ids"])
            comparison = pd.read_csv(out / "dual_route_hypothesis_comparison.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("shared_source_chain", comparison.iloc[0]["source_independence_status"])
            self.assertEqual("UniProt:Q8GRU6", comparison.iloc[0]["shared_source_ids"])

    def test_independent_source_chains_are_distinguished(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a = root / "a.tsv"
            self._write_route(route_a, [{
                "hypothesis_id": "hyp_001", "route": "phenotype_first",
                "pathway_or_process": "CLAVATA-WUSCHEL", "process_id": "",
                "gene_family": "CLV", "member_interval_genes": "G_CLV1",
                "anchor_status": "interval_anchored", "relationship_tier": "exact_ortholog",
                "source_backing": "UniProt:Q8GRU6"}])
            route_b = root / "b.tsv"
            self._write_route(route_b, [{
                "hypothesis_id": "rev_001", "route": "interval_first",
                "pathway_or_process": "CLAVATA-WUSCHEL", "process_id": "",
                "gene_family": "CLV", "member_interval_genes": "G_CLV1",
                "anchor_status": "interval_anchored", "relationship_tier": "exact_ortholog",
                "source_backing": "PMID:1001"}])
            out = self._compare(root, route_a, route_b)
            convergence = json.loads((out / "dual_route_convergence.json").read_text(encoding="utf-8"))
            self.assertEqual("independent_source_chains", convergence["independence"]["source_independence_status"])
            self.assertEqual([], convergence["independence"]["shared_source_ids"])
            comparison = pd.read_csv(out / "dual_route_hypothesis_comparison.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("independent_source_chains", comparison.iloc[0]["source_independence_status"])

    def test_b0_phenotype_blind_vs_b1_overlay_is_recorded(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_a = root / "route_a.tsv"
            self._write_route(route_a, [{
                "hypothesis_id": "hyp_001", "route": "phenotype_first",
                "pathway_or_process": "CLAVATA-WUSCHEL", "process_id": "",
                "gene_family": "CLV", "member_interval_genes": "G_CLV1;G_CLV2",
                "anchor_status": "interval_anchored", "relationship_tier": "exact_ortholog",
                "source_backing": "verified_bait"}])
            fixture = IntervalReverseLookupTests()._fixture(root)
            out_b0 = root / "route_b0"
            self.assertEqual(0, main(["interval-reverse-lookup", "--atlas", str(fixture["atlas.tsv"]),
                                      "--reference-function", str(fixture["functions.tsv"]),
                                      "--out-dir", str(out_b0)]))
            out_b1 = root / "route_b1"
            self.assertEqual(0, main(["interval-reverse-lookup", "--atlas", str(fixture["atlas.tsv"]),
                                      "--reference-function", str(fixture["functions.tsv"]),
                                      "--phenotype", "stem fasciation", "--out-dir", str(out_b1)]))
            compare_b0 = self._compare(root, route_a, out_b0 / "family_process_hypotheses.tsv",
                                       extra=["--route-b-manifest", str(out_b0 / "interval_first_status.json"),
                                              "--benchmark-mode"])
            convergence_b0 = json.loads((compare_b0 / "dual_route_convergence.json").read_text(encoding="utf-8"))
            self.assertTrue(convergence_b0["benchmark_mode"])
            self.assertFalse(convergence_b0["independence"]["phenotype_seen_by_route_b"])
            compare_b1 = self._compare(root, route_a, out_b1 / "family_process_hypotheses.tsv",
                                       extra=["--route-b-manifest", str(out_b1 / "interval_first_status.json")])
            convergence_b1 = json.loads((compare_b1 / "dual_route_convergence.json").read_text(encoding="utf-8"))
            self.assertFalse(convergence_b1["benchmark_mode"])
            self.assertTrue(convergence_b1["independence"]["phenotype_seen_by_route_b"])


if __name__ == "__main__":
    unittest.main()
