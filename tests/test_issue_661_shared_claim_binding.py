from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.claim_binding import (
    BOUND_CLAIM_COLUMNS,
    load_bound_claims,
    route_a_authorizes,
    route_b_authorizes,
    route_c_authorizes,
    validate_bound_claim_frame,
)
from candidate_qtl.route_c_claim_binding import BINDING_COLUMNS
from candidate_qtl.llm_pathway import bait_panel_freeze
from tests.bound_claim_fixtures import bound_claim


CASES = (
    ("unrelated_phenotype", {"object_match_state": "bounded_related_trait", "binding_state": "review_required"}),
    ("different_experiment_gene", {
        "experiment_subject_gene_identity": "OTHER_GENE",
        "experiment_verification_state": "contradicted",
        "binding_state": "contradicted",
    }),
    ("gene_only_in_background", {"quote_binding_state": "background_only", "binding_state": "review_required"}),
    ("secondary_review_article", {"quote_binding_state": "review_secondary", "binding_state": "review_required"}),
    ("multiple_genes_one_perturbation", {
        "experiment_verification_state": "not_evaluated",
        "binding_state": "not_evaluated",
    }),
    ("wrong_relation_label", {
        "relation_verification_state": "contradicted",
        "binding_state": "contradicted",
    }),
    ("wrong_species_or_accession", {
        "organism_verification_state": "contradicted",
        "binding_state": "contradicted",
    }),
    ("curated_function_without_phenotype_relation", {
        "predicate_relation": "curated_function_associated",
        "object_match_state": "not_evaluated",
        "binding_state": "review_required",
    }),
)


def _generic(route: str, **overrides: str) -> dict[str, str]:
    return bound_claim(
        route=route,
        claim_id="claim-661",
        citation="PMID:661",
        subject_gene="REF661",
        organism="Arabidopsis thaliana",
        phenotype="target fasciation",
        **overrides,
    )


def _direct_c_binding(case: str) -> dict[str, str]:
    row = {column: "" for column in BINDING_COLUMNS}
    row.update({
        "record_id": "claim-661",
        "frozen_record_sha256": "f" * 64,
        "citation": "PMID:661",
        "citation_resolution_state": "verified",
        "source_organism_state": "verified_match",
        "source_accession_or_cultivar_state": "not_applicable",
        "reported_identity_type": "stable_gene_id",
        "reported_identity_value": "REF661",
        "reported_identity_verification_state": "verified_exact",
        "identity_mapping_basis": "exact_source_text",
        "experimental_target_gene_state": "verified_same_identity",
        "reported_variant_state": "not_applicable",
        "phenotype_match_state": "target_trait",
        "perturbation_to_phenotype_binding_state": "verified_supports",
        "source_scope": "abstract",
        "source_record_sha256": "s" * 64,
        "quoted_support": "The REF661 mutant produced target fasciation.",
        "locator": "abstract:sentence-3",
        "reviewed_by": "fixture-reviewer",
        "binding_status": "review_required",
        "binding_reason": case,
        "rank_contribution": "none",
    })
    if case == "different_experiment_gene":
        row["experimental_target_gene_state"] = "contradicted_different_gene"
        row["binding_status"] = "contradicted"
    elif case == "unrelated_phenotype":
        row["phenotype_match_state"] = "bounded_related_trait"
    elif case == "wrong_species_or_accession":
        row["source_organism_state"] = "contradicted"
        row["binding_status"] = "contradicted"
    elif case == "wrong_relation_label":
        row["perturbation_to_phenotype_binding_state"] = "contradicts"
        row["binding_status"] = "contradicted"
    elif case in {"gene_only_in_background", "secondary_review_article", "multiple_genes_one_perturbation"}:
        row["experimental_target_gene_state"] = "not_evaluated"
        row["perturbation_to_phenotype_binding_state"] = "not_evaluated"
    elif case == "curated_function_without_phenotype_relation":
        row["phenotype_match_state"] = "not_evaluated"
        row["perturbation_to_phenotype_binding_state"] = "not_evaluated"
    return row


class SharedClaimBindingTests(unittest.TestCase):
    def test_header_only_binding_artifact_means_no_evaluated_claims(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bound_claims.tsv"
            pd.DataFrame(columns=BOUND_CLAIM_COLUMNS).to_csv(
                path, sep="\t", index=False, lineterminator="\n"
            )
            self.assertEqual({}, load_bound_claims(str(path), "a"))

    def test_positive_controls_authorize_all_three_route_adapters(self) -> None:
        a = _generic("a")
        b = _generic("b")
        self.assertTrue(validate_bound_claim_frame(pd.DataFrame([a], columns=BOUND_CLAIM_COLUMNS)).shape[0])
        self.assertTrue(route_a_authorizes(
            a, hypothesis_id="claim-661", proposed_accession="REF661",
            verified_source_ids={"PMID:661"},
        )[0])
        self.assertTrue(route_b_authorizes(
            b, claim_id="claim-661", reference_gene="REF661", source_id="PMID:661",
        )[0])
        direct = _direct_c_binding("positive")
        direct["binding_status"] = "verified_complete"
        projection = {
            "projection_basis": "exact_focal_id",
            "source_reference_gene_id": "REF661",
            "focal_gene_id": "REF661",
        }
        self.assertTrue(route_c_authorizes(direct, projection)[0])

    def test_same_eight_adversarial_cases_fail_closed_in_a_b_and_c(self) -> None:
        projection = {
            "projection_basis": "exact_focal_id",
            "source_reference_gene_id": "REF661",
            "focal_gene_id": "REF661",
        }
        for name, overrides in CASES:
            with self.subTest(case=name, route="a"):
                row = _generic("a", **overrides)
                validate_bound_claim_frame(pd.DataFrame([row], columns=BOUND_CLAIM_COLUMNS))
                self.assertFalse(route_a_authorizes(
                    row, hypothesis_id="claim-661", proposed_accession="REF661",
                    verified_source_ids={"PMID:661"},
                )[0])
            with self.subTest(case=name, route="b"):
                row = _generic("b", **overrides)
                validate_bound_claim_frame(pd.DataFrame([row], columns=BOUND_CLAIM_COLUMNS))
                self.assertFalse(route_b_authorizes(
                    row, claim_id="claim-661", reference_gene="REF661", source_id="PMID:661",
                )[0])
            with self.subTest(case=name, route="c"):
                self.assertFalse(route_c_authorizes(_direct_c_binding(name), projection)[0])

    def test_missing_binding_is_not_negative_evidence(self) -> None:
        self.assertEqual("not_evaluated", route_a_authorizes(
            None, hypothesis_id="h", proposed_accession="P1", verified_source_ids={"PMID:1"},
        )[1])
        self.assertEqual("not_evaluated", route_b_authorizes(
            None, claim_id="c", reference_gene="G1", source_id="PMID:1",
        )[1])

    def test_route_a_verified_accession_without_binding_never_freezes_as_bait(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sequence = "MAAA"
            hypotheses = root / "hypotheses.tsv"
            pd.DataFrame([{
                "hypothesis_id": "h1",
                "reference_gene_name": "REF1",
                "reference_species": "Arabidopsis thaliana",
            }]).to_csv(hypotheses, sep="\t", index=False)
            validation = root / "validation.tsv"
            pd.DataFrame([{
                "hypothesis_id": "h1",
                "verification_status": "verified",
                "proposed_stable_id": "REF1",
                "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                "bound_claim_authorized": False,
            }]).to_csv(validation, sep="\t", index=False)
            sequences = root / "sequences.json"
            sequences.write_text(json.dumps({"h1": sequence}), encoding="utf-8")
            result = bait_panel_freeze(
                str(hypotheses), str(validation), str(sequences), str(root / "out")
            )
            panel = pd.read_csv(result["panel"], sep="\t", dtype=str).fillna("")
            review = pd.read_csv(result["review_required"], sep="\t", dtype=str).fillna("")
            self.assertTrue(panel.empty)
            self.assertIn(
                "canonical_support_blocked:source_claim_binding_not_verified",
                set(review["reason"]),
            )
            self.assertIn(
                "review_projection_blocked:identity_species_sequence_not_verified",
                set(review["reason"]),
            )
            projection = pd.read_csv(result["projection_panel"], sep="\t", dtype=str).fillna("")
            self.assertTrue(projection.empty)

    def test_unnamed_accession_only_identity_can_remain_review_only(self) -> None:
        """An accession-only resolver receipt has no gene-name claim to reject."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sequence = "MAAA"
            hypotheses = root / "hypotheses.tsv"
            pd.DataFrame([{"hypothesis_id": "h1", "reference_gene_name": "", "reference_species": "Arabidopsis thaliana"}]).to_csv(
                hypotheses, sep="\t", index=False
            )
            validation = root / "validation.tsv"
            pd.DataFrame([{
                "hypothesis_id": "h1", "verification_status": "review_required", "proposed_stable_id": "P12345",
                "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(), "bound_claim_authorized": False,
                "identity_sequence_verification_status": "verified", "retrieved_species": "Arabidopsis thaliana",
                "source_provider": "local_record_cache", "source_url": "cache.json", "retrieved_at": "2026-01-01",
                "function_source_ids": "PMID:1", "reference_identity_status": "accession_only",
            }]).to_csv(validation, sep="\t", index=False)
            sequences = root / "sequences.json"
            sequences.write_text(json.dumps({"h1": sequence}), encoding="utf-8")
            result = bait_panel_freeze(
                str(hypotheses), str(validation), str(sequences), str(root / "out"),
                reviewable_sequences_arg=str(sequences),
            )
            projection = pd.read_csv(result["projection_panel"], sep="\t", dtype=str).fillna("")
            self.assertEqual(["review_only"], projection["projection_eligibility"].tolist())


if __name__ == "__main__":
    unittest.main()
