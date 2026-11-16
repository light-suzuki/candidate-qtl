"""Regression tests for issue #665: group links cannot verify a focal copy."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from candidate_qtl.route_b_semantic_transfer_verify import verify_route_b_semantic_transfer
from candidate_qtl.claim_binding import BOUND_CLAIM_COLUMNS
from candidate_qtl.semantic_claims import CLAIM_COLUMNS
from candidate_qtl.supported_member_union import build_supported_member_union
from tests.bound_claim_fixtures import bound_claim


def _claim(**overrides: str) -> dict[str, str]:
    row = {column: "" for column in CLAIM_COLUMNS}
    row.update({
        "claim_id": "claim-665",
        "route": "interval_first",
        "semantic_unit_id": "unit-665",
        "family_id": "EV1",
        "orthogroup_id": "OG1",
        "subfamily_or_clade": "SUB1",
        "member_interval_genes": "G1",
        "claim_type": "mutant_phenotype",
        "biological_function": "meristem maintenance",
        "mechanism": "receptor signaling",
        "phenotype_relation": "fasciated shoot meristem",
        "reference_gene": "REF1",
        "reference_species": "Arabidopsis thaliana",
        "evidence_scope": "interval_gene",
        "source_id": "PMID:1234567",
        "alternative_hypothesis": "related receptor copy",
        "uncertainty_note": "copy identity must be reference-bound",
        "provider": "fixture",
        "model": "fixture-model",
        "prompt_fingerprint": "p" * 64,
        "input_fingerprint": "i" * 64,
        "output_fingerprint": "o" * 64,
        "verification_status": "proposed",
        "verification_reason": "source relation pending",
        "rank_contribution": "none",
    })
    row.update(overrides)
    return row


def _write_inputs(
    root: Path,
    *,
    claim: dict[str, str],
    genes: list[dict[str, str]],
    literature_reference: str = "REF1",
    literature_source: str = "PMID:1234567",
) -> tuple[Path, Path, Path]:
    proposals = root / "proposals.tsv"
    pd.DataFrame([claim], columns=CLAIM_COLUMNS).to_csv(
        proposals, sep="\t", index=False, lineterminator="\n"
    )

    source = root / "source.tsv"
    pd.DataFrame([{
        "claim_id": claim["claim_id"],
        "source_relation_status": "supported_mechanistic",
    }]).to_csv(source, sep="\t", index=False, lineterminator="\n")

    bundle = root / "bundle.json"
    bundle.write_text(json.dumps({
        "schema_version": "1.0",
        "case_id": "case665",
        "route": "interval_first",
        "units": [{
            "unit_id": "unit-665",
            "unit_type": "family",
            "unit_anchor": "EV1",
            "member_interval_genes": [gene["gene_id"] for gene in genes],
            "identity": {
                "evolutionary_family_ids": ["EV1"],
                "family_ids": ["F1"],
                "orthogroup_ids": ["OG1"],
                "subfamily_ids": ["SUB1"],
            },
            "genes": genes,
            "literature_link_candidates": [{
                "reference_gene": literature_reference,
                "source_id": literature_source,
            }],
            "stable_source_ids": ["PMID:1234567"],
            "rank_contribution": "none",
        }],
    }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return proposals, source, bundle


def _run(
    root: Path,
    *,
    claim: dict[str, str],
    genes: list[dict[str, str]],
    literature_reference: str = "REF1",
    literature_source: str = "PMID:1234567",
) -> dict[str, object]:
    proposals, source, bundle = _write_inputs(
        root,
        claim=claim,
        genes=genes,
        literature_reference=literature_reference,
        literature_source=literature_source,
    )
    binding = root / "bound.tsv"
    pd.DataFrame([bound_claim(
        route="b",
        claim_id=claim["claim_id"],
        citation=claim["source_id"],
        subject_gene=claim["reference_gene"],
        organism=claim["reference_species"],
        phenotype=claim["phenotype_relation"],
    )], columns=BOUND_CLAIM_COLUMNS).to_csv(
        binding, sep="\t", index=False, lineterminator="\n"
    )
    return verify_route_b_semantic_transfer(
        str(proposals),
        str(source),
        str(bundle),
        str(root / "transfer"),
        bound_claims_arg=str(binding),
        case_id="case665",
    )


def test_group_level_literature_link_cannot_borrow_one_to_one_for_another_reference(
    tmp_path: Path,
) -> None:
    result = _run(
        tmp_path,
        claim=_claim(),
        genes=[{
            "gene_id": "G1",
            "orthology_reference_gene": "OTHER_REF",
            "orthology_relationship": "1:1",
        }],
    )

    claims = pd.read_csv(result["claims"], sep="\t", dtype=str).fillna("")
    transfer = pd.read_csv(result["transfer_verification"], sep="\t", dtype=str).fillna("")
    row = transfer.iloc[0]

    assert claims.iloc[0]["verification_status"] == "not_evaluated"
    assert row["transfer_status"] == "identity_not_evaluated"
    assert row["transfer_scope"] == "family_or_orthogroup"
    assert row["reference_link_basis"] == "literature_unit_membership"
    assert row["reference_link_scope"] == "group"
    assert json.loads(row["orthology_relationship_edges"]) == [{
        "focal_gene_id": "G1",
        "reference_gene": "OTHER_REF",
        "relationship": "1:1",
    }]

    union = build_supported_member_union(
        str(tmp_path / "union"),
        route_b_claims_arg=str(result["claims"]),
        case_id="case665",
    )
    union_frame = pd.read_csv(union["candidate_set"], sep="\t", dtype=str).fillna("")
    assert union_frame.empty


def test_true_reference_bound_one_to_one_still_verifies_copy_scope(tmp_path: Path) -> None:
    result = _run(
        tmp_path,
        claim=_claim(),
        genes=[{
            "gene_id": "G1",
            "orthology_reference_gene": "REF1",
            "orthology_relationship": "1:1",
        }],
    )

    claims = pd.read_csv(result["claims"], sep="\t", dtype=str).fillna("")
    transfer = pd.read_csv(result["transfer_verification"], sep="\t", dtype=str).fillna("")
    row = transfer.iloc[0]

    assert claims.iloc[0]["verification_status"] == "verified"
    assert row["transfer_status"] == "copy_supported"
    assert row["transfer_scope"] == "interval_gene"
    assert row["reference_link_basis"] == "direct_gene_orthology"
    assert row["reference_link_scope"] == "copy"


def test_group_level_literature_link_still_supports_a_family_claim(tmp_path: Path) -> None:
    result = _run(
        tmp_path,
        claim=_claim(evidence_scope="family"),
        genes=[{
            "gene_id": "G1",
            "orthology_reference_gene": "OTHER_REF",
            "orthology_relationship": "1:1",
        }],
    )

    claims = pd.read_csv(result["claims"], sep="\t", dtype=str).fillna("")
    transfer = pd.read_csv(result["transfer_verification"], sep="\t", dtype=str).fillna("")
    row = transfer.iloc[0]

    assert claims.iloc[0]["verification_status"] == "verified"
    assert row["transfer_status"] == "family_supported_copy_unresolved"
    assert row["transfer_scope"] == "family"
    assert row["supported_member_genes"] == "G1"
    assert row["reference_link_basis"] == "literature_unit_membership"
    assert row["reference_link_scope"] == "group"


def test_literature_link_requires_claim_source_provenance(tmp_path: Path) -> None:
    result = _run(
        tmp_path,
        claim=_claim(),
        genes=[{"gene_id": "G1", "orthology_reference_gene": "OTHER_REF", "orthology_relationship": "1:1"}],
        literature_source="PMID:OTHER",
    )
    claims = pd.read_csv(result["claims"], sep="\t", dtype=str).fillna("")
    transfer = pd.read_csv(result["transfer_verification"], sep="\t", dtype=str).fillna("")
    assert claims.iloc[0]["verification_status"] == "not_evaluated"
    assert transfer.iloc[0]["reference_link_basis"] == "none"
    assert "identity link" in transfer.iloc[0]["verification_reason"]


def test_endpoint_species_and_accession_namespace_are_checked_when_supplied(tmp_path: Path) -> None:
    result = _run(
        tmp_path,
        claim=_claim(reference_gene="sp|REF1|NAME", reference_species="Arabidopsis thaliana"),
        genes=[{
            "gene_id": "G1",
            "orthology_reference_gene": "sp|REF1|NAME",
            "orthology_reference_species": "Pisum sativum",
            "orthology_reference_accession_namespace": "sp",
            "orthology_relationship": "1:1",
        }],
    )
    claims = pd.read_csv(result["claims"], sep="\t", dtype=str).fillna("")
    assert claims.iloc[0]["verification_status"] == "not_evaluated"


def test_endpoint_species_and_namespace_match_allows_copy_verification(tmp_path: Path) -> None:
    result = _run(
        tmp_path,
        claim=_claim(reference_gene="sp|REF1|NAME", reference_species="Arabidopsis thaliana"),
        genes=[{
            "gene_id": "G1",
            "orthology_reference_gene": "sp|REF1|NAME",
            "orthology_reference_species": "Arabidopsis thaliana",
            "orthology_reference_accession_namespace": "sp",
            "orthology_relationship": "1:1",
        }],
    )
    claims = pd.read_csv(result["claims"], sep="\t", dtype=str).fillna("")
    assert claims.iloc[0]["verification_status"] == "verified"
