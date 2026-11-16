"""Structured source-claim binding required before Direct-C support promotion.

The generic literature verifier establishes citation resolution and performs
conservative text triage. It does not prove that the experimented gene, source
organism, perturbation and phenotype belong to one claim. This module defines
the stricter, Route-C-specific sidecar contract required before a direct report
may enter canonical supported membership (issue #661).

The first contract is deliberately fail closed. It supports exact/mapped gene
or locus identities for the target trait from abstract/full-text/curated source
records. Variant-only and merely-related-trait records remain reviewable until
a later version defines an explicit variant-to-gene and bounded-trait mapping.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

from .core import runtime_path, sha256

BINDING_CONTRACT = "route-c-source-claim-binding-1.0.0"
RANK_CONTRIBUTION = "none"

BINDING_COLUMNS = (
    "record_id",
    "frozen_record_sha256",
    "citation",
    "citation_resolution_state",
    "source_organism_state",
    "source_accession_or_cultivar_state",
    "reported_identity_type",
    "reported_identity_value",
    "reported_identity_verification_state",
    "identity_mapping_basis",
    "experimental_target_gene_state",
    "reported_variant_state",
    "phenotype_match_state",
    "perturbation_to_phenotype_binding_state",
    "source_scope",
    "source_record_sha256",
    "quoted_support",
    "locator",
    "reviewed_by",
    "binding_status",
    "binding_reason",
    "rank_contribution",
)

CITATION_RESOLUTION_STATES = frozenset({"verified", "rejected", "not_evaluated"})
SOURCE_ORGANISM_STATES = frozenset({"verified_match", "contradicted", "not_evaluated"})
SOURCE_ACCESSION_STATES = frozenset(
    {"verified_match", "not_applicable", "contradicted", "not_evaluated"}
)
REPORTED_IDENTITY_TYPES = frozenset({"stable_gene_id", "gene_symbol", "locus", "variant"})
IDENTITY_VERIFICATION_STATES = frozenset(
    {"verified_exact", "verified_mapped", "contradicted", "not_evaluated"}
)
IDENTITY_MAPPING_BASES = frozenset(
    {"exact_source_text", "curated_id_map", "source_record_namespace", "none"}
)
EXPERIMENTAL_TARGET_STATES = frozenset(
    {"verified_same_identity", "contradicted_different_gene", "not_evaluated"}
)
REPORTED_VARIANT_STATES = frozenset(
    {"verified", "not_applicable", "contradicted", "not_evaluated"}
)
PHENOTYPE_MATCH_STATES = frozenset(
    {"target_trait", "bounded_related_trait", "contradicted", "not_evaluated"}
)
PERTURBATION_BINDING_STATES = frozenset(
    {
        "verified_supports",
        "contradicts",
        "no_effect",
        "failed_rescue",
        "not_evaluated",
    }
)
SOURCE_SCOPES = frozenset({"title", "abstract", "full_text", "curated_db_record"})
BINDING_STATUSES = frozenset(
    {"verified_complete", "review_required", "contradicted", "not_evaluated"}
)

_ENUMS = {
    "citation_resolution_state": CITATION_RESOLUTION_STATES,
    "source_organism_state": SOURCE_ORGANISM_STATES,
    "source_accession_or_cultivar_state": SOURCE_ACCESSION_STATES,
    "reported_identity_type": REPORTED_IDENTITY_TYPES,
    "reported_identity_verification_state": IDENTITY_VERIFICATION_STATES,
    "identity_mapping_basis": IDENTITY_MAPPING_BASES,
    "experimental_target_gene_state": EXPERIMENTAL_TARGET_STATES,
    "reported_variant_state": REPORTED_VARIANT_STATES,
    "phenotype_match_state": PHENOTYPE_MATCH_STATES,
    "perturbation_to_phenotype_binding_state": PERTURBATION_BINDING_STATES,
    "source_scope": SOURCE_SCOPES,
    "binding_status": BINDING_STATUSES,
}


def _truthy(value: object) -> bool:
    return str(value).strip().casefold() in {"true", "1", "yes"}


def _identity_key(value: object) -> str:
    return " ".join(str(value or "").strip().casefold().split())


def _resolve_artifact_path(value: object, manifest_path: Path) -> Path:
    path = runtime_path(str(value or ""))
    if not path.is_absolute():
        path = manifest_path.parent / path
    return path


def binding_is_complete(row: dict[str, object] | None) -> bool:
    """Return whether one row meets the v1 source-claim binding contract."""
    if not row or str(row.get("binding_status", "")).strip() != "verified_complete":
        return False
    if str(row.get("citation_resolution_state", "")).strip() != "verified":
        return False
    if str(row.get("source_organism_state", "")).strip() != "verified_match":
        return False
    if str(row.get("source_accession_or_cultivar_state", "")).strip() not in {
        "verified_match",
        "not_applicable",
    }:
        return False
    if str(row.get("reported_identity_type", "")).strip() == "variant":
        # Variant-only projection is not implemented in the v1 Direct-C member
        # path. A future contract must preserve variant -> gene/locus mapping.
        return False
    if str(row.get("reported_identity_verification_state", "")).strip() not in {
        "verified_exact",
        "verified_mapped",
    }:
        return False
    if str(row.get("identity_mapping_basis", "")).strip() == "none":
        return False
    if str(row.get("experimental_target_gene_state", "")).strip() != "verified_same_identity":
        return False
    if str(row.get("reported_variant_state", "")).strip() not in {
        "verified",
        "not_applicable",
    }:
        return False
    # The first contract admits only the actual query trait. A merely related
    # trait needs a separately verified bounded relationship before promotion.
    if str(row.get("phenotype_match_state", "")).strip() != "target_trait":
        return False
    if (
        str(row.get("perturbation_to_phenotype_binding_state", "")).strip()
        != "verified_supports"
    ):
        return False
    if str(row.get("source_scope", "")).strip() not in {
        "abstract",
        "full_text",
        "curated_db_record",
    }:
        return False
    return bool(
        str(row.get("frozen_record_sha256", "")).strip()
        and str(row.get("source_record_sha256", "")).strip()
        and str(row.get("quoted_support", "")).strip()
        and str(row.get("reviewed_by", "")).strip()
    )


def canonical_frozen_record_sha256(row: dict[str, object]) -> str:
    """Hash the complete normalized frozen Direct-C row.

    The sidecar must be bound to the row it explains, not merely to its
    ``record_id``. Including every frozen column prevents swapping a claim
    between records or changing source/species/trait metadata after review.
    """
    payload = json.dumps(
        {str(key): str(value or "").strip() for key, value in sorted(row.items())},
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def binding_authorizes_projection(
    row: dict[str, object] | None,
    projection: dict[str, object],
) -> bool:
    """Prove that the identity verified in the source is the one projected.

    A source mention of a symbol must not authorize a separately supplied stable
    focal ID. Exact focal projection therefore requires that the binding itself
    verified that stable ID. Assembly/alias projection requires a verified
    mapped stable-ID/locus relation and preserves the submitted source endpoint.
    """
    if not binding_is_complete(row):
        return False
    identity_type = str(row.get("reported_identity_type", "")).strip()
    identity_value = _identity_key(row.get("reported_identity_value", ""))
    verification = str(row.get("reported_identity_verification_state", "")).strip()
    mapping_basis = str(row.get("identity_mapping_basis", "")).strip()
    projection_basis = str(projection.get("projection_basis", "")).strip()
    source_reference = _identity_key(projection.get("source_reference_gene_id", ""))
    focal_gene = _identity_key(projection.get("focal_gene_id", ""))

    if projection_basis == "exact_focal_id":
        return (
            identity_type == "stable_gene_id"
            and verification == "verified_exact"
            and bool(identity_value)
            and identity_value == focal_gene
        )
    if projection_basis == "assembly_or_alias_mapping":
        return (
            identity_type in {"stable_gene_id", "locus"}
            and verification == "verified_mapped"
            and mapping_basis in {"curated_id_map", "source_record_namespace"}
            and bool(identity_value)
            and identity_value == source_reference
        )
    return False


def validate_claim_binding_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate and canonical-sort one claim-binding sidecar."""
    missing = sorted(set(BINDING_COLUMNS) - set(frame.columns))
    extra = sorted(set(frame.columns) - set(BINDING_COLUMNS))
    if missing or extra:
        parts = []
        if missing:
            parts.append(f"missing: {', '.join(missing)}")
        if extra:
            parts.append(f"outside the contract: {', '.join(extra)}")
        raise ValueError(f"Route C claim binding has {'; '.join(parts)}")
    normalized = frame[list(BINDING_COLUMNS)].copy().fillna("")
    for column in BINDING_COLUMNS:
        normalized[column] = normalized[column].astype(str).str.strip()
    if normalized.empty:
        raise ValueError("Route C claim binding is empty")
    if normalized["record_id"].eq("").any() or normalized["record_id"].duplicated().any():
        raise ValueError("Route C claim binding requires unique non-empty record_id")
    for column in (
        "citation",
        "reported_identity_value",
        "binding_reason",
        "rank_contribution",
        "frozen_record_sha256",
    ):
        if normalized[column].eq("").any():
            raise ValueError(f"Route C claim binding requires non-empty {column}")
    for column, allowed in _ENUMS.items():
        invalid = sorted(set(normalized[column]) - allowed)
        if invalid:
            raise ValueError(
                f"Route C claim binding has invalid {column}: {', '.join(invalid)}; "
                f"allowed={', '.join(sorted(allowed))}"
            )
    if not normalized["rank_contribution"].eq(RANK_CONTRIBUTION).all():
        raise ValueError("Route C claim binding rank_contribution must always be 'none'")

    for row in normalized.to_dict("records"):
        complete = binding_is_complete(row)
        if row["binding_status"] == "verified_complete" and not complete:
            raise ValueError(
                f"Route C claim binding {row['record_id']} is marked verified_complete "
                "but one or more required gene/species/trait/perturbation axes are not verified"
            )
        contradicted_axis = any(
            str(row[column]).startswith("contradicted")
            or str(row[column]) == "contradicts"
            for column in (
                "source_organism_state",
                "source_accession_or_cultivar_state",
                "reported_identity_verification_state",
                "experimental_target_gene_state",
                "reported_variant_state",
                "phenotype_match_state",
                "perturbation_to_phenotype_binding_state",
            )
        )
        if contradicted_axis and row["binding_status"] != "contradicted":
            raise ValueError(
                f"Route C claim binding {row['record_id']} contains a contradicted axis "
                "but binding_status is not contradicted"
            )
    return normalized.sort_values("record_id", kind="stable").reset_index(drop=True)


def load_claim_binding_from_freeze_manifest(
    freeze_manifest: dict[str, object],
    freeze_manifest_path: Path,
    expected_record_ids: set[str],
    expected_frozen_rows: dict[str, dict[str, object]] | None = None,
) -> tuple[dict[str, dict[str, str]], dict[str, object]]:
    """Load a sealed binding sidecar, or report that an old generation lacks one.

    Absence is an epistemic missingness state and therefore returns an empty map;
    Direct-C members remain review-only. A supplied descriptor that is corrupt,
    has an unsupported contract, or does not match the frozen record roster is a
    lineage error and fails closed.
    """
    descriptor = freeze_manifest.get("claim_binding")
    if descriptor in (None, "", {}):
        return {}, {
            "generation_state": "not_evaluated_missing_binding",
            "contract": "",
            "artifact": "",
            "sha256": "",
            "complete": False,
            "record_count": 0,
        }
    if not isinstance(descriptor, dict):
        raise ValueError("Route C freeze manifest claim_binding must be an object")
    contract = str(descriptor.get("contract", "")).strip()
    if contract != BINDING_CONTRACT:
        raise ValueError(
            f"unsupported Route C claim-binding contract: {contract or '<empty>'}; "
            f"expected {BINDING_CONTRACT}"
        )
    artifact_value = str(descriptor.get("artifact", "")).strip()
    expected_sha = str(descriptor.get("sha256", "")).strip().removeprefix("sha256:")
    if not artifact_value or not expected_sha:
        raise ValueError("Route C claim-binding descriptor requires artifact and sha256")
    artifact = _resolve_artifact_path(artifact_value, freeze_manifest_path)
    if not artifact.exists():
        raise ValueError(f"Route C claim-binding artifact does not exist: {artifact}")
    actual_sha = sha256(artifact)
    if actual_sha != expected_sha:
        raise ValueError("Route C claim-binding artifact hash does not match freeze manifest")
    frame = validate_claim_binding_frame(
        pd.read_csv(artifact, sep="\t", dtype=str, keep_default_na=False).fillna("")
    )
    actual_ids = set(frame["record_id"].astype(str))
    if actual_ids != expected_record_ids:
        missing = sorted(expected_record_ids - actual_ids)
        extra = sorted(actual_ids - expected_record_ids)
        raise ValueError(
            f"Route C claim-binding/frozen record roster mismatch: missing={missing}, extra={extra}"
        )
    if expected_frozen_rows is None:
        raise ValueError(
            "Route C claim binding requires the corresponding frozen Direct-C rows"
        )
    for row in frame.to_dict("records"):
        frozen_row = expected_frozen_rows.get(str(row["record_id"]))
        if frozen_row is None:
            raise ValueError(
                f"Route C claim binding has no corresponding frozen row for record_id={row['record_id']}"
            )
        expected_hash = canonical_frozen_record_sha256(frozen_row)
        if row["frozen_record_sha256"].removeprefix("sha256:") != expected_hash:
            raise ValueError(
                f"Route C claim binding frozen-row hash mismatch for record_id={row['record_id']}"
            )
    try:
        declared_count = int(str(descriptor.get("record_count", "")))
    except ValueError as error:
        raise ValueError("Route C claim-binding record_count must be an integer") from error
    if declared_count != len(frame):
        raise ValueError("Route C claim-binding record_count does not match artifact")
    generation_complete = _truthy(descriptor.get("complete", False))
    all_rows_complete = all(
        binding_is_complete(row) for row in frame.to_dict("records")
    )
    if generation_complete and not all_rows_complete:
        raise ValueError(
            "Route C claim-binding generation is marked complete but contains incomplete rows"
        )
    return {
        str(row["record_id"]): {key: str(value) for key, value in row.items()}
        for row in frame.to_dict("records")
    }, {
        "generation_state": (
            "verified_complete"
            if generation_complete and all_rows_complete
            else "incomplete"
        ),
        "contract": contract,
        "artifact": str(artifact),
        "sha256": actual_sha,
        "complete": generation_complete and all_rows_complete,
        "record_count": int(len(frame)),
        "status_counts": {
            str(key): int(value)
            for key, value in frame["binding_status"].value_counts().sort_index().items()
        },
    }


def canonical_binding_json(row: dict[str, object]) -> str:
    """Return deterministic JSON for embedding one binding row in provenance."""
    return json.dumps(
        {column: str(row.get(column, "")) for column in BINDING_COLUMNS},
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
