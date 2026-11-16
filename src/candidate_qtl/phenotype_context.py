"""Structured phenotype context alongside free-text trait descriptions.

This module keeps context preservation separate from ranking.  It never replaces
free text with an inferred ontology term and never changes a candidate score.
Instead it records the biological/experimental context that makes a phenotype
statement interpretable and emits a deterministic context fingerprint for audit
and comparison.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd
import yaml

PHENOTYPE_CONTEXT_FIELDS = [
    "case_id",
    "phenotype_text",
    "trait_ontology_id",
    "trait_ontology_term",
    "unmapped_term",
    "developmental_stage",
    "tissue",
    "environment_or_treatment",
    "measurement_method",
    "comparison_or_control",
    "effect_direction",
    "replicate_info",
    "source_id",
    "source_detail",
]
PHENOTYPE_CONTEXT_COLUMNS = PHENOTYPE_CONTEXT_FIELDS + ["term_usage", "context_fingerprint"]
TERM_USAGE = frozenset({"normalized_term", "unmapped_term", "raw_text", "unknown"})

_FINGERPRINT_FIELDS = [
    "phenotype_text",
    "trait_ontology_id",
    "trait_ontology_term",
    "unmapped_term",
    "developmental_stage",
    "tissue",
    "environment_or_treatment",
    "measurement_method",
    "comparison_or_control",
    "effect_direction",
    "replicate_info",
]


# Scalar types that may carry a pandas/numpy missing marker. Array-like values
# (lists, ndarrays, dicts, strings) are deliberately absent: they are never run
# through ``pd.isna``, whose result on an array would be ambiguous.
_MISSING_SCALAR_TYPES: tuple[type, ...] = (
    float, int, complex, bool,
    np.floating, np.integer, np.complexfloating, np.bool_,
    np.datetime64, np.timedelta64,
    pd.Timestamp, pd.Timedelta,
)


def _is_missing_scalar(value: object) -> bool:
    if isinstance(value, _MISSING_SCALAR_TYPES) or value is pd.NaT or isinstance(value, type(pd.NA)):
        try:
            return bool(pd.isna(value))
        except (TypeError, ValueError):
            return False
    return False


def _text(value: object) -> str:
    if value is None or _is_missing_scalar(value):
        return "unknown"
    text = str(value).strip()
    return text if text else "unknown"


def _usage(normalized: Mapping[str, str]) -> str:
    if normalized["trait_ontology_id"] != "unknown" or normalized["trait_ontology_term"] != "unknown":
        return "normalized_term"
    if normalized["unmapped_term"] != "unknown":
        return "unmapped_term"
    if normalized["phenotype_text"] != "unknown":
        return "raw_text"
    return "unknown"


def context_fingerprint(context: Mapping[str, object]) -> str:
    """Hash only phenotype semantics/context, excluding case/source provenance."""
    payload = {field: _text(context.get(field)) for field in _FINGERPRINT_FIELDS}
    canonical = json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def normalize_phenotype_context(
    raw: Mapping[str, object] | None,
    *,
    case_id: str = "",
    phenotype_text: str = "",
) -> dict[str, str]:
    """Return one canonical optional phenotype-context record.

    ``phenotype_text`` is a compatibility bridge for existing case manifests
    that only carry free text.  Structured fields remain explicit ``unknown``;
    absence is never interpreted as negative biological evidence.

    Only ``None`` and an empty mapping count as the intended empty context;
    any other non-mapping payload fails closed instead of being silently
    converted to ``{}``.
    """
    if raw is None:
        source: dict[str, object] = {}
    elif isinstance(raw, Mapping):
        source = dict(raw)
    else:
        raise ValueError("phenotype_context must be a mapping (object)")
    case_id = str(case_id).strip()
    if phenotype_text and not str(source.get("phenotype_text", "")).strip():
        source["phenotype_text"] = phenotype_text
    normalized = {field: _text(source.get(field)) for field in PHENOTYPE_CONTEXT_FIELDS}
    if case_id:
        incoming_case = normalized["case_id"]
        if incoming_case not in {"unknown", case_id}:
            raise ValueError(f"phenotype_context case_id {incoming_case!r} does not match {case_id!r}")
        normalized["case_id"] = case_id
    normalized["term_usage"] = _usage(normalized)
    normalized["context_fingerprint"] = context_fingerprint(normalized)
    return normalized


def load_phenotype_context(
    path_arg: str | Path | None,
    *,
    case_id: str = "",
    phenotype_text: str = "",
) -> dict[str, str]:
    """Load ``phenotype_context.yaml``; absence preserves legacy free text.

    Fail-closed contract: a missing file, a non-.yaml/.yml suffix, a falsy
    non-mapping YAML root (``[]`` / ``false`` / ``0`` / ``""``), and a
    top-level ``case_id`` that conflicts with the caller's case all raise
    ``ValueError`` before any normalization.  Only an empty document
    (``None`` from ``yaml.safe_load`` on an empty file) or ``{}`` is the
    intended empty context.
    """
    if path_arg is None:
        return normalize_phenotype_context({}, case_id=case_id, phenotype_text=phenotype_text)
    path = Path(path_arg)
    if not path.exists():
        raise ValueError(f"phenotype context file not found: {path}")
    if path.suffix.casefold() not in {".yaml", ".yml"}:
        raise ValueError("phenotype context must be a .yaml/.yml file")
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if payload is None:
        payload = {}
    if not isinstance(payload, dict):
        raise ValueError("phenotype-context YAML must contain an object")
    top_case = str(payload.get("case_id", "")).strip()
    caller_case = str(case_id).strip()
    if caller_case and top_case and caller_case != top_case:
        raise ValueError(f"phenotype_context case_id {top_case!r} does not match {caller_case!r}")
    nested = payload.get("phenotype_context", payload)
    if not isinstance(nested, dict):
        raise ValueError("phenotype_context must be an object")
    trait = payload.get("trait")
    if trait is not None and not isinstance(trait, dict):
        raise ValueError("trait must be an object when present")
    legacy_text = str(phenotype_text).strip()
    if not legacy_text:
        legacy_text = str(nested.get("phenotype_text", "")).strip()
    if not legacy_text and trait:
        legacy_text = str(trait.get("phenotype_text", "")).strip()
    return normalize_phenotype_context(nested, case_id=caller_case or top_case, phenotype_text=legacy_text)


def context_frame(context: Mapping[str, object]) -> pd.DataFrame:
    """Materialize one canonical context record for audit/reporting."""
    normalized = normalize_phenotype_context(
        context,
        case_id=str(context.get("case_id", "")).strip(),
        phenotype_text=str(context.get("phenotype_text", "")).strip(),
    )
    return pd.DataFrame([normalized], columns=PHENOTYPE_CONTEXT_COLUMNS)


def structured_context_status(context: Mapping[str, object]) -> dict[str, str]:
    """Summarize whether structured context is present without scoring it."""
    normalized = normalize_phenotype_context(
        context,
        case_id=str(context.get("case_id", "")).strip(),
        phenotype_text=str(context.get("phenotype_text", "")).strip(),
    )
    structured_fields = [
        "trait_ontology_id",
        "trait_ontology_term",
        "unmapped_term",
        "developmental_stage",
        "tissue",
        "environment_or_treatment",
        "measurement_method",
        "comparison_or_control",
        "effect_direction",
        "replicate_info",
    ]
    present = [field for field in structured_fields if normalized[field] != "unknown"]
    return {
        "status": "structured_context_used" if present else "structured_context_missing",
        "term_usage": normalized["term_usage"],
        "present_fields": ";".join(present) if present else "none",
        "context_fingerprint": normalized["context_fingerprint"],
        "rank_contribution": "none",
    }


def stored_only_usage(reason: str) -> dict[str, object]:
    """The honest ``query_usage`` for a stage that records context but runs no
    query over it.  A consumer stage may only claim another ``usage`` value
    when a structured field actually reached its downstream query."""
    return {"usage": "stored_only", "reason": reason}


def structured_context_block(
    path_arg: str | Path | None,
    *,
    case_id: str = "",
    phenotype_text: str = "",
    query_usage: Mapping[str, object] | None = None,
) -> dict[str, object]:
    """One manifest-ready ``structured_context`` block for a consumer stage.

    The block is always present, even without a context file: it then records
    ``structured_context_missing`` together with the free text's term usage,
    so manifests stay machine-comparable across runs and stages.  The
    fingerprint is included so a different stage/environment input cannot be
    confused with a stale cached one; it never enters any rank/evidence hash.
    """
    context = load_phenotype_context(path_arg, case_id=case_id, phenotype_text=phenotype_text)
    status = structured_context_status(context)
    return {
        "status": status["status"],
        "term_usage": status["term_usage"],
        "phenotype_text": context["phenotype_text"],
        "unmapped_term": context["unmapped_term"],
        "present_fields": status["present_fields"],
        "context_fingerprint": status["context_fingerprint"],
        "rank_contribution": status["rank_contribution"],
        "query_usage": dict(query_usage) if query_usage is not None else stored_only_usage(
            "structured fields are recorded for audit only; no downstream query consumes them"
        ),
    }
