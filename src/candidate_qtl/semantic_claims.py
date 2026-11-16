"""Structured semantic-claim boundary for Route B (issue #403).

This module is deliberately small and provider-agnostic. Its job is to make the
boundary between an LLM proposal and machine-verified biological evidence
explicit and auditable:

    machine retrieval bundle -> LLM proposal -> machine verification
                              -> semantic_claims.tsv

It does **not** call an LLM, retrieve literature, rank genes, or choose a causal
copy. Only claims recorded as ``verified`` contribute supported membership.
Grounded proposals may separately reach Wet as explicitly unverified review
candidates. LLM prose or confidence never becomes a numerical gene score here.
"""

from __future__ import annotations

import hashlib
import json

import pandas as pd

from .core import runtime_path

SCHEMA_VERSION = "1.3"

VERIFICATION_STATUSES = (
    "proposed",
    "verified",
    "unsupported",
    "contradicted",
    "not_evaluated",
    "provider_failed",
)

EVIDENCE_SCOPES = (
    "family",
    "subfamily",
    "orthogroup",
    "reference_gene",
    "interval_gene",
    "process",
    "unknown",
)

REVIEW_ONLY_CLAIM_TYPES = frozenset({"mechanistic_review_seed"})

CLAIM_COLUMNS = [
    "claim_id",
    "route",
    "semantic_unit_id",
    "family_id",
    "orthogroup_id",
    "subfamily_or_clade",
    "member_interval_genes",
    "claim_type",
    "biological_function",
    "mechanism",
    "phenotype_relation",
    "reference_gene",
    "reference_species",
    "evidence_scope",
    "source_id",
    "source_url",
    "alternative_hypothesis",
    "uncertainty_note",
    "provider",
    "model",
    "prompt_fingerprint",
    "input_fingerprint",
    "output_fingerprint",
    "verification_status",
    "verification_reason",
    "rank_contribution",
]

from .relevance_tiers import INITIAL_FIELDS
OPTIONAL_CLAIM_COLUMNS = ["source_title", "source_identifier_aliases", "candidate_admission_policy", "focal_species", *INITIAL_FIELDS]
CLAIM_COLUMNS += OPTIONAL_CLAIM_COLUMNS

# These fields define claim identity/content. ``verification_*`` is deliberately
# excluded so a proposal can move through verification without changing merely
# because its verification state changed. Source identity remains semantic
# content: two otherwise identical statements backed by different sources are
# distinct source-backed claims. semantic_unit_id is part of provenance and
# content identity because it fixes which frozen Route-B unit the model saw.
_CONTENT_FINGERPRINT_COLUMNS = [
    "route",
    "semantic_unit_id",
    "family_id",
    "orthogroup_id",
    "subfamily_or_clade",
    "member_interval_genes",
    "claim_type",
    "biological_function",
    "mechanism",
    "phenotype_relation",
    "reference_gene",
    "reference_species",
    "evidence_scope",
    "source_id",
    "source_url",
    "alternative_hypothesis",
    "uncertainty_note",
    "provider",
    "model",
    "prompt_fingerprint",
    "input_fingerprint",
    "output_fingerprint",
]


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _canonical_gene_list(value: object) -> str:
    genes = {
        token.strip()
        for token in str(value).replace(",", ";").split(";")
        if token.strip()
    }
    return ";".join(sorted(genes))


def claim_content_fingerprint(row: dict[str, object]) -> str:
    """Return a deterministic fingerprint for one claim's semantic content."""
    payload = {
        column: (
            _canonical_gene_list(row.get(column, ""))
            if column == "member_interval_genes"
            else str(row.get(column, "")).strip()
        )
        for column in _CONTENT_FINGERPRINT_COLUMNS
    }
    payload.update({key: str(row[key]).strip() for key in OPTIONAL_CLAIM_COLUMNS if str(row.get(key, "")).strip()})
    return _sha256_text(json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":")))


def normalize_semantic_claims(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate and canonicalize a semantic-claim table.

    Fail closed instead of silently repairing scientific meaning. A claim must
    be anchored to at least one stable biological entity (family, orthogroup,
    interval gene, or reference gene). This deliberately permits unclassified
    interval genes to remain first-class instead of requiring a synthetic
    family assignment.
    """
    missing = sorted(set(CLAIM_COLUMNS) - set(OPTIONAL_CLAIM_COLUMNS) - set(frame.columns))
    if missing:
        raise ValueError(f"semantic claims table is missing columns: {', '.join(missing)}")

    normalized = frame.reindex(columns=CLAIM_COLUMNS, fill_value="").copy().fillna("")
    for column in CLAIM_COLUMNS:
        normalized[column] = normalized[column].astype(str).str.strip()
    normalized["member_interval_genes"] = normalized["member_interval_genes"].map(_canonical_gene_list)

    if normalized.empty:
        return normalized
    if normalized["claim_id"].eq("").any():
        raise ValueError("semantic claim_id must be non-empty")
    if normalized["claim_id"].duplicated().any():
        duplicates = sorted(set(normalized.loc[normalized["claim_id"].duplicated(keep=False), "claim_id"]))
        raise ValueError(f"semantic claim_id values must be unique: {', '.join(duplicates)}")
    if normalized["route"].eq("").any():
        raise ValueError("semantic claim route must be non-empty")
    route_b = normalized["route"].eq("interval_first")
    if normalized.loc[route_b, "semantic_unit_id"].eq("").any():
        raise ValueError("interval_first semantic claims require semantic_unit_id from the frozen Route B bundle")
    if normalized["claim_type"].eq("").any():
        raise ValueError("semantic claim claim_type must be non-empty")

    anchored = (
        normalized["family_id"].ne("")
        | normalized["orthogroup_id"].ne("")
        | normalized["member_interval_genes"].ne("")
        | normalized["reference_gene"].ne("")
    )
    if (~anchored).any():
        bad_ids = ", ".join(normalized.loc[~anchored, "claim_id"].tolist())
        raise ValueError(
            "semantic claims require at least one biological anchor "
            f"(family_id, orthogroup_id, member_interval_genes, reference_gene): {bad_ids}"
        )

    bad_status = sorted(set(normalized["verification_status"]) - set(VERIFICATION_STATUSES))
    if bad_status:
        raise ValueError(
            "unknown semantic verification_status: "
            + ", ".join(bad_status)
            + f"; allowed {', '.join(VERIFICATION_STATUSES)}"
        )
    bad_scope = sorted(set(normalized["evidence_scope"]) - set(EVIDENCE_SCOPES))
    if bad_scope:
        raise ValueError(
            "unknown semantic evidence_scope: "
            + ", ".join(bad_scope)
            + f"; allowed {', '.join(EVIDENCE_SCOPES)}"
        )

    review_only = normalized["claim_type"].isin(REVIEW_ONLY_CLAIM_TYPES)
    if (review_only & normalized["verification_status"].eq("verified")).any():
        raise ValueError("review-only semantic claims cannot be verified support")
    for column in ("mechanism", "phenotype_relation", "uncertainty_note"):
        if normalized.loc[review_only, column].eq("").any():
            raise ValueError(f"review-only semantic claims require a non-empty {column}")

    family_scope = normalized["evidence_scope"].eq("family")
    if normalized.loc[family_scope, "family_id"].eq("").any():
        raise ValueError("family-scoped semantic claims require family_id")
    orthogroup_scope = normalized["evidence_scope"].eq("orthogroup")
    if normalized.loc[orthogroup_scope, "orthogroup_id"].eq("").any():
        raise ValueError("orthogroup-scoped semantic claims require orthogroup_id")
    interval_gene_scope = normalized["evidence_scope"].eq("interval_gene")
    if normalized.loc[interval_gene_scope, "member_interval_genes"].eq("").any():
        raise ValueError("interval_gene-scoped semantic claims require member_interval_genes")
    reference_gene_scope = normalized["evidence_scope"].eq("reference_gene")
    if normalized.loc[reference_gene_scope, "reference_gene"].eq("").any():
        raise ValueError("reference_gene-scoped semantic claims require reference_gene")

    if not normalized["rank_contribution"].eq("none").all():
        raise ValueError("semantic claims are non-scoring; rank_contribution must be 'none' on every row")

    verified = normalized["verification_status"].eq("verified")
    if normalized.loc[verified, "source_id"].eq("").any():
        raise ValueError("verified semantic claims require a non-empty source_id")
    if normalized.loc[verified, "verification_reason"].eq("").any():
        raise ValueError("verified semantic claims require a non-empty verification_reason")

    provider_failed = normalized["verification_status"].eq("provider_failed")
    if normalized.loc[provider_failed, "verification_reason"].eq("").any():
        raise ValueError("provider_failed semantic claims require a verification_reason")

    return normalized.sort_values("claim_id", kind="stable").reset_index(drop=True)


def narrowing_eligible_claims(frame: pd.DataFrame) -> pd.DataFrame:
    """Return only machine-verified claims that may inform candidate narrowing.

    This is a *set handoff*, not a score. Multiple families/copies may remain;
    callers must not turn row count or LLM confidence into numerical evidence.
    """
    normalized = normalize_semantic_claims(frame)
    return normalized[normalized["verification_status"].eq("verified")].copy().reset_index(drop=True)


def freeze_semantic_claims(
    frame: pd.DataFrame,
    out_dir_arg: str,
    *,
    case_id: str = "",
    route: str = "interval_first",
) -> dict[str, object]:
    """Write canonical claims + manifest and freeze their artifact hash."""
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    normalized = normalize_semantic_claims(frame)
    if not normalized.empty and not normalized["route"].eq(route).all():
        raise ValueError(f"semantic claim route mismatch: expected all rows to be {route!r}")

    claims_path = out / "semantic_claims.tsv"
    normalized.to_csv(claims_path, sep="\t", index=False, lineterminator="\n")
    artifact_sha256 = hashlib.sha256(claims_path.read_bytes()).hexdigest()
    status_counts = {
        status: int(normalized["verification_status"].eq(status).sum())
        for status in VERIFICATION_STATUSES
    }
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id,
        "route": route,
        "claim_count": int(len(normalized)),
        "semantic_unit_count": int(normalized["semantic_unit_id"].nunique()) if len(normalized) else 0,
        "verification_status_counts": status_counts,
        "verified_claim_count": status_counts["verified"],
        "narrowing_eligible_claim_count": status_counts["verified"],
        "artifact_sha256": artifact_sha256,
        "output": str(claims_path),
        "policy": [
            "LLM output is a proposal until deterministic verification records verification_status=verified",
            "Route B claims retain the semantic_unit_id of the exact frozen unit shown to the model",
            "verified claims require a source_id and verification_reason",
            "family, orthogroup, interval-gene and reference-gene identities remain separate anchors",
            "unclassified interval genes remain representable and are never silently dropped for lacking family assignment",
            "review-only claim types require a mechanism, phenotype relation and uncertainty note and can never become verified support",
            "unsupported/contradicted/not_evaluated/provider_failed remain distinct states",
            "provider failure and not_evaluated are never biological negatives",
            "rank_contribution is always none; semantic claims never write a numerical gene score",
            "narrowing consumes the verified claim set, not LLM confidence or route agreement",
        ],
    }
    manifest_path = out / "semantic_claims_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        "claims": str(claims_path),
        "manifest": str(manifest_path),
        "claim_count": int(len(normalized)),
        "verified_claim_count": status_counts["verified"],
        "artifact_sha256": artifact_sha256,
    }
