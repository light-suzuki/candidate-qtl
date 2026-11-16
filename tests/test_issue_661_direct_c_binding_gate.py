"""Adversarial Direct-C support tests for the issue #661 safety gate."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from candidate_qtl.core import sha256
from candidate_qtl.route_c_claim_binding import (
    BINDING_COLUMNS,
    BINDING_CONTRACT,
    canonical_frozen_record_sha256,
)
from candidate_qtl.route_c_direct_members import build_route_c_direct_members
from candidate_qtl.route_c_direct_rescue import DIRECT_RESCUE_COLUMNS, PROJECTION_COLUMNS
from candidate_qtl.supported_member_union import build_supported_member_union

FAKE_ID = "pissa.Cameor.gnm1.ann1.PSAT_FAKE_EXACT"


def _write_spoof_case(
    root: Path,
    *,
    binding_identity_type: str | None = None,
    binding_identity_value: str = "",
) -> tuple[Path, Path, Path, Path, Path]:
    record = {column: "" for column in DIRECT_RESCUE_COLUMNS}
    record.update({
        "record_id": "spoof-1",
        "query_phenotype": "fasciated pea shoot",
        "focal_species": "Pisum sativum",
        "reported_phenotype": "fasciated shoot apical meristem",
        "reference_species": "Pisum sativum",
        "reported_gene_id": FAKE_ID,
        "reported_gene_symbol": "CLV1",
        "species_distance": "S1_same_species",
        "evidence_directness": "D0_causal_perturbation_or_rescue",
        "evidence_type": "knockout",
        "source_id": "PMID:12345678",
        "pmid": "12345678",
        "source_title": "Arabidopsis CLV1 knockout paper",
        "source_year": "2020",
        "source_verification_state": "verified",
        "gene_identity_state": "verified_gene",
        "trait_relation_state": "direct_target_trait",
        "direct_rescue_eligibility": "eligible_direct_record",
        "rank_contribution": "none",
        "notes": "old loose term/cue gate would promote this row",
    })
    frozen = root / "frozen.tsv"
    pd.DataFrame([record], columns=DIRECT_RESCUE_COLUMNS).to_csv(
        frozen, sep="\t", index=False
    )
    freeze_payload: dict[str, object] = {
        "schema_version": "1.0",
        "command": "route-c-direct-freeze",
        "frozen_sha256": sha256(frozen),
        "case_id": "spoof-case",
    }
    if binding_identity_type:
        binding_row = {column: "" for column in BINDING_COLUMNS}
        binding_row.update({
            "record_id": "spoof-1",
            "frozen_record_sha256": canonical_frozen_record_sha256(record),
            "citation": "PMID:12345678",
            "citation_resolution_state": "verified",
            "source_organism_state": "verified_match",
            "source_accession_or_cultivar_state": "not_applicable",
            "reported_identity_type": binding_identity_type,
            "reported_identity_value": binding_identity_value,
            "reported_identity_verification_state": "verified_exact",
            "identity_mapping_basis": "exact_source_text",
            "experimental_target_gene_state": "verified_same_identity",
            "reported_variant_state": "not_applicable",
            "phenotype_match_state": "target_trait",
            "perturbation_to_phenotype_binding_state": "verified_supports",
            "source_scope": "abstract",
            "source_record_sha256": "a" * 64,
            "quoted_support": "CLV1 knockout produced a fasciated shoot apical meristem.",
            "locator": "abstract",
            "reviewed_by": "fixture-reviewer",
            "binding_status": "verified_complete",
            "binding_reason": "the paper binds the CLV1 symbol to the phenotype",
            "rank_contribution": "none",
        })
        binding = root / "binding.tsv"
        pd.DataFrame([binding_row], columns=BINDING_COLUMNS).to_csv(
            binding, sep="\t", index=False
        )
        freeze_payload["claim_binding"] = {
            "contract": BINDING_CONTRACT,
            "artifact": str(binding),
            "sha256": sha256(binding),
            "record_count": 1,
            "complete": True,
        }
    freeze_manifest = root / "freeze.json"
    freeze_manifest.write_text(json.dumps(freeze_payload))

    projection_row = {column: "" for column in PROJECTION_COLUMNS}
    projection_row.update({
        "record_id": "spoof-1",
        "qtl_id": "qtl-spoof",
        "focal_gene_id": FAKE_ID,
        "projection_state": "inside_current_qtl",
        "projection_basis": "exact_focal_id",
        "source_reference_gene_id": FAKE_ID,
        "ambiguity_state": "unique",
        "interval_universe_fingerprint": "sha256:atlas",
        "rank_contribution": "none",
        "notes": "exact match of the separately submitted stable ID",
    })
    projection = root / "projection.tsv"
    pd.DataFrame([projection_row], columns=PROJECTION_COLUMNS).to_csv(
        projection, sep="\t", index=False
    )
    projection_manifest = root / "projection.json"
    projection_manifest.write_text(json.dumps({
        "schema_version": "1.0",
        "command": "route-c-direct-project",
        "frozen_sha256": sha256(frozen),
        "projection_sha256": sha256(projection),
        "qtl_id": "qtl-spoof",
    }))
    atlas = root / "atlas.tsv"
    pd.DataFrame([{"gene_id": FAKE_ID, "family_id": "CLV"}]).to_csv(
        atlas, sep="\t", index=False
    )
    return frozen, freeze_manifest, projection, projection_manifest, atlas


def _run_members(root: Path, artifacts: tuple[Path, Path, Path, Path, Path]):
    frozen, freeze_manifest, projection, projection_manifest, atlas = artifacts
    result = build_route_c_direct_members(
        str(frozen),
        str(freeze_manifest),
        str(projection),
        str(projection_manifest),
        str(atlas),
        str(root / "members"),
        case_id="spoof-case",
    )
    frame = pd.read_csv(result["members"], sep="\t", dtype=str).fillna("")
    union = build_supported_member_union(
        str(root / "union"),
        route_c_members_arg=str(result["members"]),
        case_id="spoof-case",
    )
    union_frame = pd.read_csv(
        union["candidate_set"], sep="\t", dtype=str
    ).fillna("")
    return frame, union_frame


def test_real_citation_and_exact_atlas_id_without_binding_remain_review_only(
    tmp_path: Path,
) -> None:
    """The old loose source status is not self-authenticating."""
    members, union = _run_members(tmp_path, _write_spoof_case(tmp_path))
    assert members.iloc[0]["support_status"] == "reference_phenotype_review"
    assert members.iloc[0]["membership_effect"] == "review_only"
    assert members.iloc[0]["claim_binding_status"] == "not_evaluated"
    assert members.iloc[0]["claim_binding_projection_authorized"] == "False"
    assert union.empty


def test_verified_symbol_binding_cannot_authorize_a_different_stable_id_projection(
    tmp_path: Path,
) -> None:
    """Binding CLV1 does not verify the separately supplied PSAT_FAKE ID."""
    members, union = _run_members(
        tmp_path,
        _write_spoof_case(
            tmp_path,
            binding_identity_type="gene_symbol",
            binding_identity_value="CLV1",
        ),
    )
    assert members.iloc[0]["claim_binding_status"] == "verified_complete"
    assert members.iloc[0]["claim_binding_projection_authorized"] == "False"
    assert members.iloc[0]["support_status"] == "reference_phenotype_review"
    assert members.iloc[0]["membership_effect"] == "review_only"
    assert union.empty
