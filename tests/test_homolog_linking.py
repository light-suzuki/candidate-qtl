"""Tests for agent-in-the-loop homolog linking (issue #400 step 2)."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.homolog_linking import (
    calculate_counterfactual_ranks,
    export_homolog_linking_prompts,
    score_linking_and_probe_results,
    validate_and_import_linking_proposals,
)


class HomologLinkingExportTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = Path(tempfile.mkdtemp(prefix="hl-export-test-"))
        self.evidence_path = self.tmp / "04_evidence.tsv"
        self.manifest_path = self.tmp / "discover_manifest.json"

        # Synthetic evidence TSV
        edf = pd.DataFrame([
            {
                "gene_id": "Psat.cameor.v2.7g00050",
                "evidence_code": "homology_context",
                "direction": "context",
                "confidence": "medium",
                "value": "tr|G7JMM6|G7JMM6_MEDTR",
                "source_type": "local_blastp",
                "source_id": "tr|G7JMM6|G7JMM6_MEDTR",
                "description": "generic best-hit conservation: 26S proteasome regulatory subunit RPN7",
            },
            {
                "gene_id": "Psat.cameor.v2.7g03000",
                "evidence_code": "homology_context",
                "direction": "context",
                "confidence": "high",
                "value": "tr|A0A0C3WQQ0|A0A0C3WQQ0_MEDTR",
                "source_type": "local_blastp",
                "source_id": "tr|A0A0C3WQQ0|A0A0C3WQQ0_MEDTR",
                "description": "generic best-hit conservation: Beta-amyrin synthase OS=Medicago GN=bas1",
            },
        ])
        edf.to_csv(self.evidence_path, sep="\t", index=False)

        manifest = {
            "inputs": {
                "phenotype": "saponin (soyasaponin) content",
                "intervals": ["chr7:1-7646762"],
            },
            "candidate_count": 2,
        }
        self.manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    def test_export_is_self_contained_and_free_of_truth_leakage(self) -> None:
        out_dir = self.tmp / "export_out"
        res = export_homolog_linking_prompts(
            out_dir=out_dir,
            case_id="case_test",
            phenotype="",
            evidence_path=self.evidence_path,
            discover_manifest_path=self.manifest_path,
        )
        self.assertTrue(Path(res["prompt_json"]).is_file())
        self.assertTrue(Path(res["prompt_md"]).is_file())
        self.assertTrue(Path(res["manifest"]).is_file())

        prompt_data = json.loads(Path(res["prompt_json"]).read_text(encoding="utf-8"))
        self.assertEqual("saponin (soyasaponin) content", prompt_data["phenotype"])
        self.assertEqual(2, len(prompt_data["candidates"]))

        # Mechanical audit for truth leakage terms:
        prompt_text = Path(res["prompt_json"]).read_text(encoding="utf-8")
        for forbidden in ("case_selection", "truth_gene", "benchmark_score", "causal_rank"):
            self.assertNotIn(forbidden, prompt_text)


class HomologLinkingValidationGateTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = Path(tempfile.mkdtemp(prefix="hl-gates-test-"))
        self.prompt_json = self.tmp / "homolog_linking_prompt.json"
        prompt_data = {
            "schema_version": "1.0",
            "case_id": "case_bas1",
            "phenotype": "saponin content",
            "candidates": [
                {
                    "gene_id": "Psat.cameor.v2.7g00050",
                    "homology_description": "26S proteasome regulatory subunit RPN7",
                },
                {
                    "gene_id": "Psat.cameor.v2.7g03000",
                    "homology_description": "Beta-amyrin synthase OS=Medicago GN=bas1",
                },
            ],
        }
        self.prompt_json.write_text(json.dumps(prompt_data), encoding="utf-8")

    def test_closed_world_rejection_for_unknown_gene_id(self) -> None:
        resp = self.tmp / "resp_unknown_gene.json"
        resp.write_text(json.dumps({
            "case_id": "case_bas1",
            "abstain": False,
            "proposals": [
                {
                    "gene_id": "Psat.invented.gene99999",
                    "evidence_quote": "Beta-amyrin synthase",
                    "confidence_tier": "direct_catalytic_enzyme",
                    "mechanism_rationale": "Made up gene",
                }
            ],
        }), encoding="utf-8")

        with self.assertRaisesRegex(ValueError, "Closed-world violation"):
            validate_and_import_linking_proposals(
                prompt_export_manifest_or_json=self.prompt_json,
                response_payload_path=resp,
                out_dir=self.tmp / "out_closed_world",
                reviewed_by="test_reviewer",
            )

    def test_grounding_substring_rejection_for_non_matching_quote(self) -> None:
        resp = self.tmp / "resp_fake_quote.json"
        resp.write_text(json.dumps({
            "case_id": "case_bas1",
            "abstain": False,
            "proposals": [
                {
                    "gene_id": "Psat.cameor.v2.7g03000",
                    "evidence_quote": "Invented soyasaponin biosynthetic enzyme",
                    "confidence_tier": "direct_catalytic_enzyme",
                    "mechanism_rationale": "Hallucinated quote",
                }
            ],
        }), encoding="utf-8")

        with self.assertRaisesRegex(ValueError, "Grounding violation"):
            validate_and_import_linking_proposals(
                prompt_export_manifest_or_json=self.prompt_json,
                response_payload_path=resp,
                out_dir=self.tmp / "out_grounding",
                reviewed_by="test_reviewer",
            )

    def test_budget_exceeded_rejection(self) -> None:
        resp = self.tmp / "resp_overbudget.json"
        resp.write_text(json.dumps({
            "case_id": "case_bas1",
            "abstain": False,
            "proposals": [
                {"gene_id": "Psat.cameor.v2.7g03000", "evidence_quote": "Beta-amyrin synthase", "confidence_tier": "direct_catalytic_enzyme", "mechanism_rationale": "r1"},
                {"gene_id": "Psat.cameor.v2.7g00050", "evidence_quote": "proteasome", "confidence_tier": "broad_family_member", "mechanism_rationale": "r2"},
            ],
        }), encoding="utf-8")

        # Set max_proposals = 1
        with self.assertRaisesRegex(ValueError, "Budget exceeded"):
            validate_and_import_linking_proposals(
                prompt_export_manifest_or_json=self.prompt_json,
                response_payload_path=resp,
                out_dir=self.tmp / "out_budget",
                reviewed_by="test_reviewer",
                max_proposals=1,
            )

    def test_abstention_recorded_as_not_evaluated_no_link(self) -> None:
        resp = self.tmp / "resp_abstain.json"
        resp.write_text(json.dumps({
            "case_id": "case_bas1",
            "abstain": True,
            "proposals": [],
            "rejected_alternatives": [{"gene_id": "Psat.cameor.v2.7g00050", "reason": "proteasome unrelated"}],
        }), encoding="utf-8")

        res = validate_and_import_linking_proposals(
            prompt_export_manifest_or_json=self.prompt_json,
            response_payload_path=resp,
            out_dir=self.tmp / "out_abstain",
            reviewed_by="test_reviewer",
        )
        self.assertEqual("not_evaluated_no_link", res["status"])
        self.assertEqual(0, res["proposal_count"])

        mf = json.loads(Path(res["manifest"]).read_text(encoding="utf-8"))
        self.assertEqual("not_evaluated_no_link", mf["status"])
        self.assertEqual("none", mf["rank_contribution"])

    def test_successful_proposal_validation_provenance_and_reproducibility(self) -> None:
        resp = self.tmp / "resp_valid.json"
        resp.write_text(json.dumps({
            "case_id": "case_bas1",
            "abstain": False,
            "proposals": [
                {
                    "gene_id": "Psat.cameor.v2.7g03000",
                    "evidence_quote": "Beta-amyrin synthase",
                    "confidence_tier": "direct_catalytic_enzyme",
                    "mechanism_rationale": "Beta-amyrin is the direct precursor of soyasaponin.",
                }
            ],
            "rejected_alternatives": [
                {"gene_id": "Psat.cameor.v2.7g00050", "reason": "proteasome subunit unrelated to saponin"}
            ],
        }), encoding="utf-8")

        out_first = self.tmp / "out_first"
        res1 = validate_and_import_linking_proposals(
            prompt_export_manifest_or_json=self.prompt_json,
            response_payload_path=resp,
            out_dir=out_first,
            reviewed_by="reviewer_alpha",
            payload_source="agent",
            provider="subagent_test",
            model="test_model",
        )
        self.assertEqual("completed", res1["status"])
        self.assertEqual(1, res1["proposal_count"])

        out_second = self.tmp / "out_second"
        res2 = validate_and_import_linking_proposals(
            prompt_export_manifest_or_json=self.prompt_json,
            response_payload_path=resp,
            out_dir=out_second,
            reviewed_by="reviewer_alpha",
            payload_source="agent",
            provider="subagent_test",
            model="test_model",
        )

        # Byte-for-byte reproducibility check
        tsv1 = Path(res1["proposals_tsv"]).read_bytes()
        tsv2 = Path(res2["proposals_tsv"]).read_bytes()
        self.assertEqual(tsv1, tsv2)

        m1 = json.loads(Path(res1["manifest"]).read_text(encoding="utf-8"))
        m2 = json.loads(Path(res2["manifest"]).read_text(encoding="utf-8"))
        self.assertEqual(m1["prompt_sha256"], m2["prompt_sha256"])
        self.assertEqual(m1["response_sha256"], m2["response_sha256"])
        self.assertEqual(m1["status"], m2["status"])
        self.assertEqual(m1["proposal_count"], m2["proposal_count"])
        self.assertEqual(m1["rank_contribution"], m2["rank_contribution"])


class CounterfactualAndProbeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = Path(tempfile.mkdtemp(prefix="hl-cf-probe-test-"))
        self.ranked_tsv = self.tmp / "ranked_candidates.tsv"
        self.proposals_tsv = self.tmp / "homolog_linking_proposals.tsv"

        # Synthetic ranked TSV (10 candidates)
        rows = []
        for i in range(1, 11):
            gid = f"Psat.test.{i:05d}"
            rows.append({
                "gene_id": gid,
                "dry_score": "-0.562" if i == 5 else "0.0",
                "screening_rank": str(i),
            })
        pd.DataFrame(rows).to_csv(self.ranked_tsv, sep="\t", index=False)

        # Linking proposal for truth gene Psat.test.00005
        pd.DataFrame([{
            "case_id": "case_test",
            "proposal_id": "hl_01",
            "gene_id": "Psat.test.00005",
            "confidence_tier": "direct_catalytic_enzyme",
            "evidence_quote": "quote",
            "mechanism_rationale": "rationale",
            "review_required": "yes",
            "rank_contribution": "none",
            "status": "validated_proposal",
        }]).to_csv(self.proposals_tsv, sep="\t", index=False)

    def test_counterfactual_ranking_leaves_original_file_untouched(self) -> None:
        orig_bytes = self.ranked_tsv.read_bytes()
        cf_out = self.tmp / "counterfactual.tsv"

        res = calculate_counterfactual_ranks(
            ranked_tsv_path=self.ranked_tsv,
            linking_proposals_tsv_path=self.proposals_tsv,
            truth_gene_ids="Psat.test.00005",
            out_tsv_path=cf_out,
            counterfactual_bonus=0.5,
        )

        self.assertEqual(orig_bytes, self.ranked_tsv.read_bytes())
        self.assertTrue(cf_out.is_file())

        self.assertTrue(res["truth_proposed"])
        self.assertEqual(5, res["actual_rank"])
        # Score changed from -0.562 to -0.062, so new rank is higher than 0.0 scores?
        # Score -0.062 < 0.0, so it remains below 0.0 scores, but if bonus is +1.0 it would be 0.438 (rank 1)
        res_large_bonus = calculate_counterfactual_ranks(
            ranked_tsv_path=self.ranked_tsv,
            linking_proposals_tsv_path=self.proposals_tsv,
            truth_gene_ids="Psat.test.00005",
            counterfactual_bonus=1.0,
        )
        self.assertEqual(1, res_large_bonus["counterfactual_rank"])

    def test_memorization_probe_and_stratification(self) -> None:
        linking_manifest = self.tmp / "homolog_linking_manifest.json"
        linking_manifest.write_text(json.dumps({
            "case_id": "case_test",
            "status": "completed",
            "proposals_tsv": str(self.proposals_tsv),
        }), encoding="utf-8")

        probe_resp = self.tmp / "probe_resp.json"
        probe_resp.write_text(json.dumps({
            "case_id": "case_test",
            "predicted_gene_symbols": ["BAS1"],
            "predicted_enzyme_or_family": "Beta-amyrin synthase",
        }), encoding="utf-8")

        # Case 1: Probe shares a known identity key with truth -> match /
        # recall_positive / not_evaluated_memorization.
        res_recalled = score_linking_and_probe_results(
            linking_manifest_path=linking_manifest,
            memorization_probe_response_path=probe_resp,
            truth_gene_ids="Psat.test.00005",
            truth_family_or_symbol="BAS1",
        )
        self.assertEqual("match", res_recalled["probe_match_status"])
        self.assertEqual("recall_positive", res_recalled["stratification"])
        self.assertEqual("not_evaluated_memorization", res_recalled["verdict"])

        # Case 2: NO probe response at all -> the probe relation cannot be
        # evaluated, so it must NOT fall to recall_negative / reasoned_link_success
        # (the old default misreported memory as reasoning).
        res_no_probe = score_linking_and_probe_results(
            linking_manifest_path=linking_manifest,
            memorization_probe_response_path=None,
            truth_gene_ids="Psat.test.00005",
            truth_family_or_symbol="BAS1",
        )
        self.assertEqual("not_evaluated_probe_ambiguous", res_no_probe["probe_match_status"])
        self.assertEqual("not_evaluated_probe_ambiguous", res_no_probe["stratification"])
        self.assertEqual("not_evaluated_probe_ambiguous", res_no_probe["verdict"])

    def test_probe_alias_matching_brc1_tb1_tcp18(self) -> None:
        """PsBRC1 / TB1 / TCP18 are the same shoot-branching TCP factor and
        must be judged a probe match with auditable provenance (issue #419)."""
        linking_manifest = self.tmp / "homolog_linking_manifest.json"
        linking_manifest.write_text(json.dumps({
            "case_id": "pea_brc1",
            "status": "completed",
            "proposals_tsv": str(self.proposals_tsv),
        }), encoding="utf-8")

        probe_resp = self.tmp / "probe_brc1.json"
        probe_resp.write_text(json.dumps({
            "case_id": "pea_brc1",
            "predicted_gene_symbols": ["PsBRC1", "BRC1", "TB1"],
            "predicted_enzyme_or_family": "TCP family transcription factor (BRANCHED1/TB1 clade)",
        }), encoding="utf-8")

        res = score_linking_and_probe_results(
            linking_manifest_path=linking_manifest,
            memorization_probe_response_path=probe_resp,
            truth_gene_ids="Psat.test.00005",
            truth_family_or_symbol="TCP18",
        )
        self.assertEqual("match", res["probe_match_status"])
        self.assertEqual("recall_positive", res["stratification"])
        self.assertEqual("not_evaluated_memorization", res["verdict"])

        evidence = res["probe_match_evidence"]
        self.assertTrue(evidence)
        first = evidence[0]
        self.assertEqual("PsBRC1", first["term"])
        self.assertEqual("alias:brc1", first["matched_key"])
        self.assertIn("A1YKT1", first["provenance"])
        self.assertIn("Q93WI2", first["provenance"])

        # The score artifact row carries the same evidence as a column.
        tsv_path = Path(res["probe_score_tsv"])
        self.assertTrue(tsv_path.is_file())
        score_df = pd.read_csv(tsv_path, sep="\t", dtype=str)
        self.assertEqual(1, len(score_df))
        self.assertEqual("match", score_df.loc[0, "probe_match_status"])
        self.assertIn("A1YKT1", score_df.loc[0, "probe_match_evidence"])

    def test_probe_alias_outside_registry_is_ambiguous_not_negative(self) -> None:
        """A concrete probe answer that shares NO known identity key with the
        truth family must not fall to recall_negative; it is excluded as
        not_evaluated_probe_ambiguous (issue #419)."""
        linking_manifest = self.tmp / "homolog_linking_manifest.json"
        linking_manifest.write_text(json.dumps({
            "case_id": "case_test",
            "status": "completed",
            "proposals_tsv": str(self.proposals_tsv),
        }), encoding="utf-8")

        probe_resp = self.tmp / "probe_offtarget.json"
        probe_resp.write_text(json.dumps({
            "case_id": "case_test",
            "predicted_gene_symbols": ["COMT"],
            "predicted_enzyme_or_family": "caffeic acid O-methyltransferase",
            "confidence": "high",
        }), encoding="utf-8")

        res = score_linking_and_probe_results(
            linking_manifest_path=linking_manifest,
            memorization_probe_response_path=probe_resp,
            truth_gene_ids="Psat.test.00005",
            truth_family_or_symbol="CHS",
        )
        self.assertEqual("not_evaluated_probe_ambiguous", res["probe_match_status"])
        self.assertEqual("not_evaluated_probe_ambiguous", res["stratification"])
        self.assertEqual("not_evaluated_probe_ambiguous", res["verdict"])
        self.assertEqual([], res["probe_match_evidence"])
        self.assertFalse(res["probe_abstained"])

    def test_probe_explicit_abstention_is_mismatch(self) -> None:
        """Only an explicit probe abstention is a genuine non-recall and may
        support recall_negative / reasoned_link_success."""
        linking_manifest = self.tmp / "homolog_linking_manifest.json"
        linking_manifest.write_text(json.dumps({
            "case_id": "case_test",
            "status": "completed",
            "proposals_tsv": str(self.proposals_tsv),
        }), encoding="utf-8")

        probe_resp = self.tmp / "probe_unknown.json"
        probe_resp.write_text(json.dumps({
            "case_id": "case_test",
            "predicted_gene_symbols": [],
            "predicted_enzyme_or_family": "unknown",
            "confidence": "unknown",
        }), encoding="utf-8")

        res = score_linking_and_probe_results(
            linking_manifest_path=linking_manifest,
            memorization_probe_response_path=probe_resp,
            truth_gene_ids="Psat.test.00005",
            truth_family_or_symbol="BAS1",
        )
        self.assertEqual("mismatch", res["probe_match_status"])
        self.assertEqual("recall_negative", res["stratification"])
        self.assertEqual("reasoned_link_success", res["verdict"])
        self.assertTrue(res["probe_abstained"])

    def test_probe_unparseable_response_is_ambiguous_not_mismatch(self) -> None:
        """A probe response that cannot be read as JSON is NOT an explicit
        abstention: it must stay not_evaluated_probe_ambiguous and can never
        become reasoned_link_success (an unreadable answer is not evidence
        that the model reasoned)."""
        linking_manifest = self.tmp / "homolog_linking_manifest.json"
        linking_manifest.write_text(json.dumps({
            "case_id": "case_test",
            "status": "completed",
            "proposals_tsv": str(self.proposals_tsv),
        }), encoding="utf-8")

        bad_resp = self.tmp / "probe_bad.json"
        bad_resp.write_text("this is not json { broken", encoding="utf-8")

        res = score_linking_and_probe_results(
            linking_manifest_path=linking_manifest,
            memorization_probe_response_path=bad_resp,
            truth_gene_ids="Psat.test.00005",
            truth_family_or_symbol="BAS1",
        )
        self.assertEqual("not_evaluated_probe_ambiguous", res["probe_match_status"])
        self.assertEqual("not_evaluated_probe_ambiguous", res["stratification"])
        self.assertEqual("not_evaluated_probe_ambiguous", res["verdict"])
        self.assertFalse(res["probe_abstained"])

        empty_resp = self.tmp / "probe_empty.json"
        empty_resp.write_text("", encoding="utf-8")
        res_empty = score_linking_and_probe_results(
            linking_manifest_path=linking_manifest,
            memorization_probe_response_path=empty_resp,
            truth_gene_ids="Psat.test.00005",
            truth_family_or_symbol="BAS1",
        )
        self.assertEqual("not_evaluated_probe_ambiguous", res_empty["probe_match_status"])
        self.assertEqual("not_evaluated_probe_ambiguous", res_empty["verdict"])

    def test_probe_descriptive_family_text_can_match_symbol_truth(self) -> None:
        """The probe's free-text family label is matched with the expanded
        (Route A) matcher, so a descriptive answer meets a symbol truth."""
        linking_manifest = self.tmp / "homolog_linking_manifest.json"
        linking_manifest.write_text(json.dumps({
            "case_id": "pea_brc1",
            "status": "completed",
            "proposals_tsv": str(self.proposals_tsv),
        }), encoding="utf-8")

        probe_resp = self.tmp / "probe_descriptive.json"
        probe_resp.write_text(json.dumps({
            "case_id": "pea_brc1",
            "predicted_gene_symbols": [],
            "predicted_enzyme_or_family": "TCP family transcription factor (BRANCHED1/TB1 clade)",
            "confidence": "high",
        }), encoding="utf-8")

        res = score_linking_and_probe_results(
            linking_manifest_path=linking_manifest,
            memorization_probe_response_path=probe_resp,
            truth_gene_ids="Psat.test.00005",
            truth_family_or_symbol="TCP18",
        )
        self.assertEqual("match", res["probe_match_status"])
        self.assertEqual("recall_positive", res["stratification"])
        self.assertTrue(res["probe_match_evidence"])


if __name__ == "__main__":
    unittest.main()
