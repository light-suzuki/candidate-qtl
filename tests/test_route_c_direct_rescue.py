from __future__ import annotations

import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.route_c_direct_rescue import (
    DIRECT_RESCUE_COLUMNS,
    PROJECTION_COLUMNS,
    direct_record_can_seed_projection,
    validate_direct_rescue_frame,
    validate_projection_frame,
)

REPO_ROOT = Path(__file__).resolve().parents[1]


def _direct_row(**overrides: str) -> dict[str, str]:
    row = {column: "" for column in DIRECT_RESCUE_COLUMNS}
    row.update(
        {
            "record_id": "c0001",
            "query_phenotype": "fasciated stem",
            "focal_species": "Pisum sativum",
            "focal_cultivar_or_accession": "example accession",
            "reported_phenotype": "fasciated shoot apical meristem",
            "reference_species": "Arabidopsis thaliana",
            "reported_gene_id": "AT1G00010",
            "reported_gene_symbol": "A",
            "species_distance": "S5_model_or_other_plant",
            "evidence_directness": "D0_causal_perturbation_or_rescue",
            "evidence_type": "loss_of_function_mutant",
            "source_id": "PMID:12345678",
            "pmid": "12345678",
            "source_title": "Example direct phenotype paper",
            "source_year": "2026",
            "source_verification_state": "verified",
            "gene_identity_state": "verified_gene",
            "trait_relation_state": "direct_target_trait",
            "reference_sequence_id": "UniProt:EXAMPLE",
            "direct_rescue_eligibility": "eligible_direct_record",
            "rank_contribution": "none",
        }
    )
    row.update(overrides)
    return row


def _projection_row(**overrides: str) -> dict[str, str]:
    row = {column: "" for column in PROJECTION_COLUMNS}
    row.update(
        {
            "record_id": "c0001",
            "qtl_id": "qtl_chr4",
            "focal_gene_id": "Psat.example.gene1",
            "projection_state": "inside_current_qtl",
            "projection_basis": "blastp",
            "source_reference_gene_id": "AT1G00010",
            "identity": "0.61",
            "query_coverage": "0.92",
            "target_coverage": "0.88",
            "ambiguity_state": "unique",
            "interval_universe_fingerprint": "sha256:example",
            "rank_contribution": "none",
        }
    )
    row.update(overrides)
    return row


class RouteCDirectRescueTests(unittest.TestCase):
    def test_direct_rescue_template_matches_schema(self) -> None:
        template = pd.read_csv(REPO_ROOT / "templates" / "route_c_direct_rescue.tsv", sep="\t", dtype=str)
        self.assertEqual(tuple(template.columns), DIRECT_RESCUE_COLUMNS)

    def test_direct_rescue_accepts_cross_species_d0_record(self) -> None:
        frame = pd.DataFrame([_direct_row()], columns=DIRECT_RESCUE_COLUMNS)
        validated = validate_direct_rescue_frame(frame)
        self.assertEqual(validated.loc[0, "focal_species"], "Pisum sativum")
        self.assertEqual(validated.loc[0, "species_distance"], "S5_model_or_other_plant")
        self.assertTrue(direct_record_can_seed_projection(validated.iloc[0].to_dict()))

    def test_direct_rescue_rejects_interval_information_before_freeze(self) -> None:
        frame = pd.DataFrame([_direct_row()], columns=DIRECT_RESCUE_COLUMNS)
        frame["qtl_id"] = "qtl_chr4"
        with self.assertRaisesRegex(ValueError, "outside the contract"):
            validate_direct_rescue_frame(frame)

    def test_direct_rescue_requires_focal_species(self) -> None:
        frame = pd.DataFrame([_direct_row(focal_species="")], columns=DIRECT_RESCUE_COLUMNS)
        with self.assertRaisesRegex(ValueError, "requires non-empty focal_species"):
            validate_direct_rescue_frame(frame)

    def test_direct_rescue_requires_explicit_species_distance(self) -> None:
        frame = pd.DataFrame([_direct_row(species_distance="")], columns=DIRECT_RESCUE_COLUMNS)
        with self.assertRaisesRegex(ValueError, "requires non-empty species_distance"):
            validate_direct_rescue_frame(frame)

    def test_direct_rescue_record_ids_are_unique_after_trimming(self) -> None:
        frame = pd.DataFrame(
            [_direct_row(record_id="c0001"), _direct_row(record_id=" c0001 ")],
            columns=DIRECT_RESCUE_COLUMNS,
        )
        with self.assertRaisesRegex(ValueError, "record_id must be unique"):
            validate_direct_rescue_frame(frame)

    def test_context_only_record_cannot_be_direct_rescue_eligible(self) -> None:
        frame = pd.DataFrame(
            [
                _direct_row(
                    evidence_directness="D3_context_only",
                    trait_relation_state="context_only",
                    direct_rescue_eligibility="eligible_direct_record",
                )
            ],
            columns=DIRECT_RESCUE_COLUMNS,
        )
        with self.assertRaisesRegex(ValueError, "D3_context_only"):
            validate_direct_rescue_frame(frame)

    def test_unverified_source_cannot_be_direct_rescue_eligible(self) -> None:
        frame = pd.DataFrame(
            [_direct_row(source_verification_state="pending")],
            columns=DIRECT_RESCUE_COLUMNS,
        )
        with self.assertRaisesRegex(ValueError, "requires source_verification_state=verified"):
            validate_direct_rescue_frame(frame)

    def test_projection_keeps_outside_qtl_as_evaluated_state(self) -> None:
        frame = pd.DataFrame(
            [_projection_row(projection_state="outside_current_qtl")],
            columns=PROJECTION_COLUMNS,
        )
        validated = validate_projection_frame(frame)
        self.assertEqual(validated.loc[0, "projection_state"], "outside_current_qtl")
        self.assertEqual(validated.loc[0, "rank_contribution"], "none")

    def test_no_interval_match_is_not_an_outside_qtl_claim(self) -> None:
        frame = pd.DataFrame(
            [
                _projection_row(
                    focal_gene_id="",
                    projection_state="no_interval_match",
                    projection_basis="blastp",
                    identity="",
                    query_coverage="",
                    target_coverage="",
                    ambiguity_state="no_threshold_hit",
                )
            ],
            columns=PROJECTION_COLUMNS,
        )
        validated = validate_projection_frame(frame)
        self.assertEqual(validated.loc[0, "projection_state"], "no_interval_match")
        self.assertEqual(validated.loc[0, "focal_gene_id"], "")

    def test_no_interval_match_requires_qtl_context(self) -> None:
        frame = pd.DataFrame(
            [
                _projection_row(
                    qtl_id="",
                    focal_gene_id="",
                    projection_state="no_interval_match",
                    projection_basis="blastp",
                    ambiguity_state="no_threshold_hit",
                )
            ],
            columns=PROJECTION_COLUMNS,
        )
        with self.assertRaisesRegex(ValueError, "no_interval_match requires qtl_id"):
            validate_projection_frame(frame)

    def test_inside_qtl_requires_focal_gene(self) -> None:
        frame = pd.DataFrame(
            [_projection_row(focal_gene_id="")],
            columns=PROJECTION_COLUMNS,
        )
        with self.assertRaisesRegex(ValueError, "requires focal_gene_id"):
            validate_projection_frame(frame)

    def test_not_projected_has_no_fake_projection_basis(self) -> None:
        frame = pd.DataFrame(
            [
                _projection_row(
                    qtl_id="",
                    focal_gene_id="",
                    projection_state="not_projected",
                    projection_basis="not_evaluated",
                    ambiguity_state="not_evaluated",
                )
            ],
            columns=PROJECTION_COLUMNS,
        )
        validated = validate_projection_frame(frame)
        self.assertEqual(validated.loc[0, "projection_state"], "not_projected")

    def test_projection_requires_explicit_state(self) -> None:
        frame = pd.DataFrame(
            [_projection_row(projection_state="")],
            columns=PROJECTION_COLUMNS,
        )
        with self.assertRaisesRegex(ValueError, "requires non-empty projection_state"):
            validate_projection_frame(frame)


if __name__ == "__main__":
    unittest.main()
