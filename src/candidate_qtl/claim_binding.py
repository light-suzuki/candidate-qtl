"""Shared source-claim binding contract for canonical Route A/B/C support.

Citation resolution and term/cue matching are useful triage, but they do not
prove that a source's gene, experiment and phenotype form one proposition.
This module owns that stricter, route-neutral contract.  Route adapters retain
the responsibility for identity and membership semantics.
"""

from __future__ import annotations

from pathlib import Path
from dataclasses import dataclass

import pandas as pd

from .core import runtime_path

BOUND_CLAIM_CONTRACT = "source-bound-claim-1.0.0"
BOUND_CLAIM_COLUMNS = (
    "bound_claim_id",
    "route",
    "claim_id",
    "citation",
    "source_record_id",
    "subject_gene_identity",
    "subject_gene_verification_state",
    "experiment_perturbation",
    "experiment_subject_gene_identity",
    "experiment_verification_state",
    "predicate_relation",
    "relation_verification_state",
    "object_phenotype_process",
    "object_match_state",
    "organism",
    "organism_verification_state",
    "accession",
    "accession_verification_state",
    "quoted_support",
    "locator",
    "quote_binding_state",
    "source_scope",
    "binding_state",
    "curation_method",
    "reviewed_by",
    "binding_reason",
    "rank_contribution",
)

SUPPORTING_PREDICATES = frozenset({
    "causal_perturbation_supports_target_trait",
    "mapped_gene_functionally_validated_for_target_trait",
})
_ROUTES = frozenset({"a", "b", "c"})
_CURATION_METHODS = frozenset({"deterministic", "curated_human", "independent_agent_review"})


@dataclass(frozen=True)
class _AgentReviewReceipt:
    """In-memory result of checked source/review inputs, never a submitted flag."""

    claim_sha256: str
_SOURCE_SCOPES = frozenset({"title", "abstract", "full_text", "curated_db_record"})
_IDENTITY_STATES = frozenset({"verified_exact", "verified_mapped", "contradicted", "not_evaluated"})
_MATCH_STATES = frozenset({"target_trait", "bounded_related_trait", "contradicted", "not_evaluated"})
_VERIFICATION_STATES = frozenset({"verified_match", "not_applicable", "contradicted", "not_evaluated"})
_BINDING_STATES = frozenset({"verified_complete", "review_required", "contradicted", "not_evaluated"})
_EXPERIMENT_STATES = frozenset({"verified_same_subject", "contradicted", "not_evaluated"})
_RELATION_STATES = frozenset({"verified", "contradicted", "not_evaluated"})
_QUOTE_BINDING_STATES = frozenset({
    "verified_same_claim", "background_only", "review_secondary", "contradicted", "not_evaluated"
})


def validate_bound_claim_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate and canonicalize a shared bound-claim TSV."""
    missing = sorted(set(BOUND_CLAIM_COLUMNS) - set(frame.columns))
    extra = sorted(set(frame.columns) - set(BOUND_CLAIM_COLUMNS))
    if missing or extra:
        details = []
        if missing:
            details.append(f"missing: {', '.join(missing)}")
        if extra:
            details.append(f"outside the contract: {', '.join(extra)}")
        raise ValueError(f"bound-claim table has {'; '.join(details)}")
    result = frame[list(BOUND_CLAIM_COLUMNS)].copy().fillna("")
    for column in BOUND_CLAIM_COLUMNS:
        result[column] = result[column].astype(str).str.strip()
    for column in ("bound_claim_id", "claim_id", "binding_reason", "rank_contribution"):
        if result[column].eq("").any():
            raise ValueError(f"bound-claim table requires non-empty {column}")
    if result["bound_claim_id"].duplicated().any():
        raise ValueError("bound_claim_id values must be unique")
    if result[["route", "claim_id"]].duplicated().any():
        raise ValueError("bound-claim table requires one row per route/claim_id")
    enums = {
        "route": _ROUTES,
        "subject_gene_verification_state": _IDENTITY_STATES,
        "experiment_verification_state": _EXPERIMENT_STATES,
        "relation_verification_state": _RELATION_STATES,
        "object_match_state": _MATCH_STATES,
        "organism_verification_state": _VERIFICATION_STATES,
        "accession_verification_state": _VERIFICATION_STATES,
        "source_scope": _SOURCE_SCOPES,
        "quote_binding_state": _QUOTE_BINDING_STATES,
        "binding_state": _BINDING_STATES,
        "curation_method": _CURATION_METHODS,
    }
    for column, allowed in enums.items():
        invalid = sorted(set(result[column]) - allowed)
        if invalid:
            raise ValueError(f"bound-claim table has invalid {column}: {', '.join(invalid)}")
    if not result["rank_contribution"].eq("none").all():
        raise ValueError("bound-claim rank_contribution must always be 'none'")
    for row in result.to_dict("records"):
        authorized, _state, reason = bound_claim_authorizes(row)
        if row["binding_state"] == "verified_complete" and not authorized:
            raise ValueError(
                f"bound claim {row['bound_claim_id']} is marked verified_complete but is incomplete: {reason}"
            )
    return result.sort_values(["route", "claim_id"], kind="stable").reset_index(drop=True)


def load_bound_claims(path_arg: str | None, route: str) -> dict[str, dict[str, str]]:
    """Load one route's claims. Absence means not evaluated, never support."""
    if not path_arg:
        return {}
    path: Path = runtime_path(path_arg)
    if not path.exists():
        raise ValueError(f"bound-claim table does not exist: {path}")
    frame = validate_bound_claim_frame(pd.read_csv(path, sep="\t", dtype=str).fillna(""))
    selected = frame[frame["route"] == route]
    return {
        str(row["claim_id"]): {key: str(value) for key, value in row.items()}
        for row in selected.to_dict("records")
    }


def bound_claim_authorizes(row: dict[str, object] | None) -> tuple[bool, str, str]:
    """Authorize only a fully bound gene-experiment-phenotype proposition."""
    if not row:
        return False, "not_evaluated", "no shared bound claim was supplied"
    state = str(row.get("binding_state", "")).strip()
    if state == "contradicted":
        return False, "contradicted", str(row.get("binding_reason", "")) or "bound claim is contradicted"
    if state != "verified_complete":
        return False, state or "not_evaluated", str(row.get("binding_reason", "")) or "bound claim is incomplete"
    required = (
        "citation", "source_record_id", "subject_gene_identity", "experiment_perturbation",
        "experiment_subject_gene_identity",
        "predicate_relation", "object_phenotype_process", "organism", "quoted_support",
        "locator", "reviewed_by",
    )
    missing = [column for column in required if not str(row.get(column, "")).strip()]
    if missing:
        return False, "not_evaluated", f"bound claim is missing: {', '.join(missing)}"
    if str(row.get("subject_gene_verification_state", "")) not in {"verified_exact", "verified_mapped"}:
        return False, "review_required", "subject gene identity is not verified"
    if (
        str(row.get("experiment_verification_state", "")) != "verified_same_subject"
        or str(row.get("experiment_subject_gene_identity", "")).casefold()
        != str(row.get("subject_gene_identity", "")).casefold()
    ):
        return False, "review_required", "experiment is not bound to the verified subject gene"
    if str(row.get("relation_verification_state", "")) != "verified":
        return False, "review_required", "submitted predicate/relation is not independently verified"
    if str(row.get("organism_verification_state", "")) != "verified_match":
        return False, "review_required", "source organism is not a verified match"
    if str(row.get("accession_verification_state", "")) not in {"verified_match", "not_applicable"}:
        return False, "review_required", "source accession/cultivar is not verified"
    if str(row.get("object_match_state", "")) != "target_trait":
        return False, "review_required", "object is not verified as the target trait"
    if str(row.get("predicate_relation", "")) not in SUPPORTING_PREDICATES:
        return False, "review_required", "predicate is not a verified supporting relation"
    if str(row.get("source_scope", "")) == "title":
        return False, "review_required", "title-only evidence cannot bind an experiment"
    if str(row.get("quote_binding_state", "")) != "verified_same_claim":
        return False, "review_required", "quote/span does not bind the same gene-experiment-phenotype claim"
    if str(row.get("curation_method", "")) not in _CURATION_METHODS:
        return False, "not_evaluated", "curation method is not accepted"
    if row.get("curation_method") == "independent_agent_review":
        from .source_review_capture import content_sha256
        receipt = row.get("_agent_review_receipt")
        claim = {key: value for key, value in row.items() if key != "_agent_review_receipt"}
        if not isinstance(receipt, _AgentReviewReceipt) or receipt.claim_sha256 != content_sha256(claim):
            return False, "not_evaluated", "agent review requires replayed machine checks, not a submitted method label"
    return True, "verified_complete", "gene, experiment, relation, phenotype, organism and source are bound"


def route_a_authorizes(
    row: dict[str, object] | None,
    *,
    hypothesis_id: str,
    proposed_accession: str,
    verified_source_ids: set[str],
) -> tuple[bool, str, str]:
    ok, state, reason = bound_claim_authorizes(row)
    if not ok:
        return ok, state, reason
    if str(row.get("route")) != "a" or str(row.get("claim_id")) != hypothesis_id:
        return False, "review_required", "bound claim does not identify this Route-A hypothesis"
    if str(row.get("subject_gene_identity")).casefold() != proposed_accession.casefold():
        return False, "review_required", "bound subject gene differs from the verified proposed accession"
    if str(row.get("citation")) not in verified_source_ids:
        return False, "review_required", "bound citation is absent from verified proposal sources"
    return True, state, reason


def route_b_authorizes(
    row: dict[str, object] | None,
    *,
    claim_id: str,
    reference_gene: str,
    source_id: str,
) -> tuple[bool, str, str]:
    ok, state, reason = bound_claim_authorizes(row)
    if not ok:
        return ok, state, reason
    if str(row.get("route")) != "b" or str(row.get("claim_id")) != claim_id:
        return False, "review_required", "bound claim does not identify this Route-B claim"
    if str(row.get("subject_gene_identity")).casefold() != reference_gene.casefold():
        return False, "review_required", "bound subject gene differs from the submitted reference gene"
    if str(row.get("citation")) != source_id:
        return False, "review_required", "bound citation differs from the submitted source"
    return True, state, reason


def route_c_authorizes(
    binding: dict[str, object] | None,
    projection: dict[str, object],
) -> tuple[bool, str, str]:
    """Adapt the sealed Direct-C sidecar to the shared relation contract."""
    from .route_c_claim_binding import binding_authorizes_projection, binding_is_complete

    if not binding_is_complete(binding):
        return False, "not_evaluated", "Direct-C claim binding is incomplete"
    if not binding_authorizes_projection(binding, projection):
        return False, "review_required", "bound source identity does not authorize this projection"
    return True, "verified_complete", "sealed Direct-C binding and projection authorize support"
