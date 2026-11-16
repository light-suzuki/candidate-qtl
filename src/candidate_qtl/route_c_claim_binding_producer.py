"""Prepare and seal reviewed Route-C source-claim binding sidecars.

Generic citation verification and verbatim quote matching are triage layers.
They do not prove that one gene, organism, perturbation and phenotype form the
same claim. This producer turns an explicit reviewed packet into the sealed
sidecar required by ``route_c_direct_members`` (issue #661).

The producer never infers a positive biological relation. A row can become
``verified_complete`` only when an accepted source-verification extraction is
present and a named reviewer explicitly attests every load-bearing axis. The
reviewer may be a human (``curated_human``) or an independent agent
(``independent_agent_review``); both pass the same machine checks of quote,
citation, frozen identity and retrieved scope, and the method is recorded.
"""

from __future__ import annotations

from .citation_tokens import split_citations

import argparse
import copy
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

from .core import runtime_path, sha256
from .route_c_claim_binding import (
    BINDING_COLUMNS,
    BINDING_CONTRACT,
    EXPERIMENTAL_TARGET_STATES,
    IDENTITY_MAPPING_BASES,
    IDENTITY_VERIFICATION_STATES,
    PERTURBATION_BINDING_STATES,
    PHENOTYPE_MATCH_STATES,
    REPORTED_IDENTITY_TYPES,
    REPORTED_VARIANT_STATES,
    SOURCE_ACCESSION_STATES,
    SOURCE_ORGANISM_STATES,
    SOURCE_SCOPES,
    binding_is_complete,
    canonical_frozen_record_sha256,
    validate_claim_binding_frame,
)
from .route_c_direct_rescue import validate_direct_rescue_frame
from .source_verify import classify_citation

PRODUCER_CONTRACT = "route-c-claim-binding-producer-1.0.0"
# Both reviewer kinds pass the same machine checks before a binding can be
# complete: the quote must be an accepted verbatim extraction of the selected
# frozen source, the citation must resolve to a verified source record, the
# identity must equal the frozen record and the scope must have been retrieved.
# An agent reviewer is accepted because those checks, not the reviewer label,
# carry the authority; an unlabelled or anonymous review stays review_required.
REVIEWER_CURATION_METHODS = frozenset({"curated_human", "independent_agent_review"})
CURATION_METHODS = frozenset({"", *REVIEWER_CURATION_METHODS})

REVIEW_COLUMNS = (
    "record_id",
    "selected_citation",
    "available_citations",
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
    "quoted_support",
    "locator",
    "curation_method",
    "reviewed_by",
    "binding_reason",
)

_REVIEW_ENUMS = {
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
    "curation_method": CURATION_METHODS,
}

_IDENTITY_FIELD_BY_TYPE = {
    "stable_gene_id": "reported_gene_id",
    "gene_symbol": "reported_gene_symbol",
    "locus": "reported_locus",
    "variant": "reported_variant",
}


def _resolve(value: object, parent: Path) -> Path:
    path = runtime_path(str(value or ""))
    if not path.is_absolute():
        path = parent / path
    return path


def _normalise_space(value: object) -> str:
    return re.sub(r"\s+", " ", str(value or "")).strip()


def _canonical_json_sha256(value: object) -> str:
    payload = json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _citation_key(value: object) -> str:
    kind, identifier = classify_citation(str(value or ""))
    if kind == "doi":
        identifier = identifier.casefold()
    elif kind == "pmcid":
        identifier = identifier.upper()
    return f"{kind}:{identifier}"


def _split_citations(value: object) -> tuple[str, ...]:
    values = set(split_citations(value))
    return tuple(sorted(values, key=lambda item: (_citation_key(item), item)))


def _read_json(path: Path, label: str) -> dict[str, object]:
    if not path.exists():
        raise ValueError(f"{label} does not exist: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is not valid JSON: {path}") from error
    if not isinstance(payload, dict):
        raise ValueError(f"{label} must contain a JSON object")
    return payload


def _require_sha(path: Path, expected: object, label: str) -> None:
    digest = str(expected or "").removeprefix("sha256:")
    if len(digest) != 64 or sha256(path) != digest:
        raise ValueError(f"{label} hash does not match its manifest")


def _source_context(
    frozen_arg: str,
    freeze_manifest_arg: str,
) -> tuple[
    Path,
    Path,
    dict[str, object],
    pd.DataFrame,
    pd.DataFrame,
    dict[str, object],
    Path,
    Path,
    Path,
]:
    frozen_path = runtime_path(frozen_arg)
    freeze_manifest_path = runtime_path(freeze_manifest_arg)
    freeze_manifest = _read_json(freeze_manifest_path, "Route C freeze manifest")
    if freeze_manifest.get("command") != "route-c-direct-freeze":
        raise ValueError("claim binding requires a route-c-direct-freeze manifest")
    _require_sha(
        frozen_path,
        freeze_manifest.get("frozen_sha256", ""),
        "Route C frozen artifact",
    )
    direct = validate_direct_rescue_frame(
        pd.read_csv(frozen_path, sep="\t", dtype=str, keep_default_na=False).fillna("")
    )

    source_manifest_value = (
        freeze_manifest.get("source_verification_manifest")
        or (freeze_manifest.get("outputs") or {}).get("source_verification")
    )
    source_manifest_path = _resolve(
        source_manifest_value, freeze_manifest_path.parent
    )
    source_manifest = _read_json(
        source_manifest_path, "Route C source-verification manifest"
    )
    if source_manifest.get("command") != "source-verify":
        raise ValueError("Route C freeze does not point to a source-verify manifest")
    outputs = source_manifest.get("outputs")
    if not isinstance(outputs, dict):
        raise ValueError("source-verify manifest lacks outputs")
    validation_path = _resolve(
        outputs.get("source_validation", ""), source_manifest_path.parent
    )
    sources_path = _resolve(
        outputs.get("verified_sources", ""), source_manifest_path.parent
    )
    _require_sha(
        validation_path,
        freeze_manifest.get("source_validation_sha256", ""),
        "Route C source-validation artifact",
    )
    _require_sha(
        sources_path,
        freeze_manifest.get("verified_sources_sha256", ""),
        "Route C verified-source artifact",
    )
    validation = pd.read_csv(validation_path, sep="\t", dtype=str, keep_default_na=False).fillna("")
    required = {
        "hypothesis_id",
        "citation",
        "verification_status",
        "citation_verified",
        "retrieval_scope",
        "extraction_status",
        "quoted_support",
    }
    missing = sorted(required - set(validation.columns))
    if missing:
        raise ValueError(
            f"source-validation artifact is missing columns: {', '.join(missing)}"
        )
    verified_sources = _read_json(
        sources_path, "Route C verified-source artifact"
    )
    return (
        frozen_path,
        freeze_manifest_path,
        freeze_manifest,
        direct,
        validation,
        verified_sources,
        source_manifest_path,
        validation_path,
        sources_path,
    )


def _identity_for_record(record: dict[str, object]) -> tuple[str, str]:
    for identity_type, field in _IDENTITY_FIELD_BY_TYPE.items():
        value = str(record.get(field, "") or "").strip()
        if value:
            return identity_type, value
    raise ValueError(f"Route C record {record.get('record_id', '')} has no identity")


def _validation_rows(
    validation: pd.DataFrame,
    record_id: str,
    allowed_citations: tuple[str, ...],
) -> pd.DataFrame:
    allowed = {_citation_key(value) for value in allowed_citations}
    rows = validation[
        validation["hypothesis_id"].astype(str).eq(record_id)
        & validation["citation"].map(_citation_key).isin(allowed)
    ].copy()
    if rows.empty:
        raise ValueError(
            f"Route C record {record_id} has no matching source-validation row"
        )
    return rows.sort_values(
        ["extraction_status", "citation"],
        ascending=[True, True],
        kind="stable",
    )


def _available_scope(retrieval_scope: object) -> str:
    value = str(retrieval_scope or "").strip().casefold()
    if "full" in value:
        return "full_text"
    if "abstract" in value:
        return "abstract"
    return "title"


def _source_record_by_citation(
    verified_sources: dict[str, object],
) -> dict[str, object]:
    result: dict[str, object] = {}
    for citation, record in verified_sources.items():
        result.setdefault(_citation_key(citation), record)
    return result


def _captured_sources(
    freeze_manifest: dict[str, object], source_manifest_path: Path,
) -> dict[str, object] | None:
    source_manifest = _read_json(source_manifest_path, "Route C source-verification manifest")
    value = (source_manifest.get("outputs") or {}).get("source_review_records")
    expected = freeze_manifest.get("source_review_records_sha256")
    if not value and not expected:
        return None  # Legacy freezes did not seal captured-source availability.
    if not value or not expected:
        raise ValueError("Route C captured-source availability requires artifact and hash")
    path = _resolve(value, source_manifest_path.parent)
    _require_sha(path, expected, "Route C captured-source artifact")
    payload = _read_json(path, "Route C captured-source artifact")
    records = payload.get("records")
    if not isinstance(records, dict):
        raise ValueError("Route C captured-source records must be an object")
    return _source_record_by_citation(records)


def _review_available_scope(
    source_row: dict[str, object], captured_sources: dict[str, object] | None,
) -> str:
    if captured_sources is None:
        return _available_scope(source_row.get("retrieval_scope", ""))
    record = captured_sources.get(_citation_key(source_row.get("citation", "")))
    if record is None:
        return "title"  # Unretrieved sources cannot authorize content review.
    if not isinstance(record, dict) or not isinstance(record.get("sections"), dict):
        raise ValueError("Route C captured-source availability lacks source sections")
    available = record.get("available_text_scope")
    scopes = {"full_text_cached": "full_text", "title_abstract": "abstract",
              "title_only": "title", "none": "title"}
    if not isinstance(available, str) or available not in scopes:
        raise ValueError("Route C captured-source availability is missing or malformed")
    sections = record["sections"]
    actual = next((scope for field, scope in (
        ("fullText", "full_text_cached"), ("abstractText", "title_abstract"),
        ("title", "title_only"),
    ) if isinstance(sections.get(field), str) and sections[field].strip()), "none")
    if actual != available:
        raise ValueError("Route C captured-source availability disagrees with source sections")
    return scopes[available]


def prepare_claim_binding_review(
    frozen_arg: str,
    freeze_manifest_arg: str,
    out_dir_arg: str,
) -> dict[str, object]:
    (
        frozen_path,
        freeze_manifest_path,
        freeze_manifest,
        direct,
        validation,
        _,
        source_manifest_path,
        validation_path,
        _,
    ) = _source_context(frozen_arg, freeze_manifest_arg)
    captured_sources = _captured_sources(freeze_manifest, source_manifest_path)
    rows: list[dict[str, str]] = []
    for record in direct.sort_values("record_id", kind="stable").to_dict("records"):
        record_id = str(record["record_id"]).strip()
        citations = _split_citations(record.get("source_id", ""))
        checked = _validation_rows(validation, record_id, citations)
        accepted = checked[
            checked["citation_verified"]
            .astype(str)
            .str.casefold()
            .isin({"true", "1"})
            & checked["extraction_status"].astype(str).eq("accepted")
            & checked["quoted_support"].astype(str).str.strip().ne("")
        ]
        selected = accepted.iloc[0] if not accepted.empty else checked.iloc[0]
        identity_type, identity_value = _identity_for_record(record)
        rows.append(
            {
                "record_id": record_id,
                "selected_citation": str(selected["citation"]),
                "available_citations": ";".join(
                    sorted(
                        set(checked["citation"].astype(str)),
                        key=lambda value: (_citation_key(value), value),
                    )
                ),
                "source_organism_state": "not_evaluated",
                "source_accession_or_cultivar_state": (
                    "not_applicable"
                    if not str(
                        record.get("reference_cultivar_or_accession", "")
                    ).strip()
                    else "not_evaluated"
                ),
                "reported_identity_type": identity_type,
                "reported_identity_value": identity_value,
                "reported_identity_verification_state": "not_evaluated",
                "identity_mapping_basis": "none",
                "experimental_target_gene_state": "not_evaluated",
                "reported_variant_state": (
                    "not_applicable"
                    if not str(record.get("reported_variant", "")).strip()
                    else "not_evaluated"
                ),
                "phenotype_match_state": "not_evaluated",
                "perturbation_to_phenotype_binding_state": "not_evaluated",
                "source_scope": _review_available_scope(selected.to_dict(), captured_sources),
                "quoted_support": str(selected["quoted_support"]),
                "locator": str(selected["retrieval_scope"]),
                "curation_method": "",
                "reviewed_by": "",
                "binding_reason": "review not completed",
            }
        )
    review = pd.DataFrame(rows, columns=REVIEW_COLUMNS)
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    review_path = out / "route_c_claim_binding_review.tsv"
    review.to_csv(review_path, sep="\t", index=False, lineterminator="\n")
    manifest = {
        "contract": PRODUCER_CONTRACT,
        "stage": "prepare_claim_binding_review",
        "record_count": int(len(review)),
        "interval_blind": True,
        "positive_relation_inference": "none",
        "inputs": {
            "frozen": {
                "path": str(frozen_path),
                "sha256": sha256(frozen_path),
            },
            "freeze_manifest": {
                "path": str(freeze_manifest_path),
                "sha256": sha256(freeze_manifest_path),
            },
            "source_validation": {
                "path": str(validation_path),
                "sha256": sha256(validation_path),
            },
        },
        "outputs": {"review": str(review_path)},
    }
    manifest_path = out / "route_c_claim_binding_review_manifest.json"
    manifest["outputs"]["manifest"] = str(manifest_path)
    manifest["review_sha256"] = sha256(review_path)
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        "review": str(review_path),
        "manifest": str(manifest_path),
        "record_count": len(review),
    }


def _validate_review_frame(frame: pd.DataFrame) -> pd.DataFrame:
    missing = sorted(set(REVIEW_COLUMNS) - set(frame.columns))
    extra = sorted(set(frame.columns) - set(REVIEW_COLUMNS))
    if missing or extra:
        parts = []
        if missing:
            parts.append(f"missing: {', '.join(missing)}")
        if extra:
            parts.append(f"outside the contract: {', '.join(extra)}")
        raise ValueError(f"Route C claim-binding review has {'; '.join(parts)}")
    normalized = frame[list(REVIEW_COLUMNS)].copy().fillna("")
    for column in REVIEW_COLUMNS:
        normalized[column] = normalized[column].astype(str).str.strip()
    if normalized.empty:
        raise ValueError("Route C claim-binding review is empty")
    if (
        normalized["record_id"].eq("").any()
        or normalized["record_id"].duplicated().any()
    ):
        raise ValueError(
            "Route C claim-binding review requires unique non-empty record_id"
        )
    for column in (
        "selected_citation",
        "available_citations",
        "reported_identity_type",
        "reported_identity_value",
        "source_scope",
        "binding_reason",
    ):
        if normalized[column].eq("").any():
            raise ValueError(
                f"Route C claim-binding review requires non-empty {column}"
            )
    for column, allowed in _REVIEW_ENUMS.items():
        invalid = sorted(set(normalized[column]) - allowed)
        if invalid:
            raise ValueError(
                f"Route C claim-binding review has invalid {column}: "
                f"{', '.join(invalid)}"
            )
    return normalized.sort_values("record_id", kind="stable").reset_index(
        drop=True
    )


def _contradicted(review: dict[str, str]) -> bool:
    return any(
        str(review[column]).startswith("contradicted")
        or str(review[column]) == "contradicts"
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


def _source_row_for_review(
    checked: pd.DataFrame,
    selected_citation: str,
) -> dict[str, str]:
    selected_key = _citation_key(selected_citation)
    matching = checked[checked["citation"].map(_citation_key).eq(selected_key)]
    if matching.empty:
        raise ValueError(
            f"selected citation {selected_citation} is absent from source validation"
        )
    rows = [
        {str(key): str(value) for key, value in row.items()}
        for row in matching.to_dict("records")
    ]
    rows.sort(
        key=lambda row: (
            row.get("extraction_status") != "accepted",
            row.get("verification_status") != "verified",
            json.dumps(row, ensure_ascii=False, sort_keys=True),
        )
    )
    return rows[0]


def _scope_is_available(requested: str, available: str) -> bool:
    order = {
        "title": 0,
        "abstract": 1,
        "full_text": 2,
        "curated_db_record": 2,
    }
    return order[requested] <= order[available]


def _binding_status(
    candidate: dict[str, str],
    *,
    attested_review: bool,
    quote_verified: bool,
) -> str:
    if _contradicted(candidate):
        return "contradicted"
    probe = dict(candidate)
    probe["binding_status"] = "verified_complete"
    if attested_review and quote_verified and binding_is_complete(probe):
        return "verified_complete"
    load_bearing = (
        "source_organism_state",
        "reported_identity_verification_state",
        "experimental_target_gene_state",
        "phenotype_match_state",
        "perturbation_to_phenotype_binding_state",
    )
    if all(str(candidate[column]) == "not_evaluated" for column in load_bearing):
        return "not_evaluated"
    return "review_required"


def seal_claim_binding_review(
    frozen_arg: str,
    freeze_manifest_arg: str,
    review_arg: str,
    out_dir_arg: str,
) -> dict[str, object]:
    (
        frozen_path,
        freeze_manifest_path,
        freeze_manifest,
        direct,
        validation,
        verified_sources,
        source_manifest_path,
        validation_path,
        sources_path,
    ) = _source_context(frozen_arg, freeze_manifest_arg)
    review_path = runtime_path(review_arg)
    review = _validate_review_frame(
        pd.read_csv(review_path, sep="\t", dtype=str, keep_default_na=False).fillna("")
    )
    records = {
        str(row["record_id"]): row for row in direct.to_dict("records")
    }
    if set(review["record_id"].astype(str)) != set(records):
        missing = sorted(set(records) - set(review["record_id"].astype(str)))
        extra = sorted(set(review["record_id"].astype(str)) - set(records))
        raise ValueError(
            f"claim-binding review/frozen roster mismatch: "
            f"missing={missing}, extra={extra}"
        )
    source_records = _source_record_by_citation(verified_sources)
    captured_sources = _captured_sources(freeze_manifest, source_manifest_path)
    binding_rows: list[dict[str, str]] = []
    for reviewed in review.to_dict("records"):
        record_id = str(reviewed["record_id"])
        record = records[record_id]
        actual_citations = _split_citations(record.get("source_id", ""))
        declared_citations = _split_citations(reviewed["available_citations"])
        if {_citation_key(value) for value in declared_citations} != {
            _citation_key(value) for value in actual_citations
        }:
            raise ValueError(
                f"available citations do not match frozen record {record_id}"
            )
        selected_key = _citation_key(reviewed["selected_citation"])
        if selected_key not in {
            _citation_key(value) for value in actual_citations
        }:
            raise ValueError(
                f"selected citation does not belong to frozen record {record_id}"
            )
        checked = _validation_rows(validation, record_id, actual_citations)
        source_row = _source_row_for_review(
            checked, str(reviewed["selected_citation"])
        )
        available_scope = _review_available_scope(source_row, captured_sources)
        requested_scope = str(reviewed["source_scope"])
        if not _scope_is_available(requested_scope, available_scope):
            raise ValueError(
                f"review source_scope exceeds retrieved scope for record {record_id}"
            )
        field = _IDENTITY_FIELD_BY_TYPE[
            str(reviewed["reported_identity_type"])
        ]
        expected_identity = str(record.get(field, "") or "").strip()
        if not expected_identity or expected_identity.casefold() != str(
            reviewed["reported_identity_value"]
        ).casefold():
            raise ValueError(
                f"reviewed identity is not the frozen {field} for record {record_id}"
            )

        reviewed_quote = _normalise_space(reviewed.get("quoted_support", ""))
        source_quote = _normalise_space(source_row.get("quoted_support", ""))
        quote_verified = (
            source_row.get("extraction_status") == "accepted"
            and bool(reviewed_quote)
            and reviewed_quote == source_quote
        )
        if reviewed_quote and not quote_verified:
            raise ValueError(
                f"reviewed quote is not an accepted verbatim extraction for "
                f"record {record_id}"
            )
        citation_verified = str(
            source_row.get("citation_verified", "")
        ).casefold() in {"true", "1"}
        source_record = source_records.get(selected_key)
        citation_state = (
            "verified"
            if citation_verified and source_record is not None
            else (
                "rejected"
                if source_row.get("verification_status") == "rejected"
                else "not_evaluated"
            )
        )
        candidate = {
            "record_id": record_id,
            "frozen_record_sha256": canonical_frozen_record_sha256(record),
            "citation": str(reviewed["selected_citation"]),
            "citation_resolution_state": citation_state,
            "source_organism_state": str(reviewed["source_organism_state"]),
            "source_accession_or_cultivar_state": str(
                reviewed["source_accession_or_cultivar_state"]
            ),
            "reported_identity_type": str(
                reviewed["reported_identity_type"]
            ),
            "reported_identity_value": str(
                reviewed["reported_identity_value"]
            ),
            "reported_identity_verification_state": str(
                reviewed["reported_identity_verification_state"]
            ),
            "identity_mapping_basis": str(reviewed["identity_mapping_basis"]),
            "experimental_target_gene_state": str(
                reviewed["experimental_target_gene_state"]
            ),
            "reported_variant_state": str(
                reviewed["reported_variant_state"]
            ),
            "phenotype_match_state": str(reviewed["phenotype_match_state"]),
            "perturbation_to_phenotype_binding_state": str(
                reviewed["perturbation_to_phenotype_binding_state"]
            ),
            "source_scope": requested_scope,
            "source_record_sha256": (
                _canonical_json_sha256(source_record)
                if source_record is not None
                else ""
            ),
            "quoted_support": reviewed_quote,
            "locator": str(reviewed["locator"]),
            "reviewed_by": str(reviewed["reviewed_by"]),
            "binding_status": "review_required",
            "binding_reason": str(reviewed["binding_reason"]),
            "rank_contribution": "none",
        }
        candidate["binding_status"] = _binding_status(
            candidate,
            attested_review=(
                reviewed["curation_method"] in REVIEWER_CURATION_METHODS
                and bool(str(reviewed["reviewed_by"]).strip())
            ),
            quote_verified=quote_verified,
        )
        binding_rows.append(candidate)

    binding = validate_claim_binding_frame(
        pd.DataFrame(binding_rows, columns=BINDING_COLUMNS)
    )
    all_complete = all(
        binding_is_complete(row) for row in binding.to_dict("records")
    )
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    binding_path = out / "route_c_source_claim_binding.tsv"
    binding.to_csv(binding_path, sep="\t", index=False, lineterminator="\n")
    bound_manifest_path = (
        out / "route_c_direct_rescue_freeze_manifest.bound.json"
    )
    bound_manifest = copy.deepcopy(freeze_manifest)
    bound_manifest["claim_binding"] = {
        "contract": BINDING_CONTRACT,
        "artifact": str(binding_path),
        "sha256": sha256(binding_path),
        "record_count": int(len(binding)),
        "complete": all_complete,
    }
    bound_manifest["claim_binding_producer"] = {
        "contract": PRODUCER_CONTRACT,
        "curation_policy": (
            "accepted verbatim extraction plus an attested curated_human or "
            "independent_agent_review, both after the same machine checks; "
            "the producer does not infer a positive relation"
        ),
        "review_artifact": str(review_path),
        "review_sha256": sha256(review_path),
        "source_verification_manifest": str(source_manifest_path),
        "source_verification_manifest_sha256": sha256(source_manifest_path),
        "source_validation": str(validation_path),
        "source_validation_sha256": sha256(validation_path),
        "verified_sources": str(sources_path),
        "verified_sources_sha256": sha256(sources_path),
        "status_counts": {
            str(key): int(value)
            for key, value in binding["binding_status"]
            .value_counts()
            .sort_index()
            .items()
        },
        "curation_method_counts": {
            (str(key) or "unattested"): int(value)
            for key, value in review["curation_method"]
            .astype(str)
            .value_counts()
            .sort_index()
            .items()
        },
    }
    bound_manifest.setdefault("outputs", {})["manifest"] = str(
        bound_manifest_path
    )
    bound_manifest["outputs"]["claim_binding"] = str(binding_path)
    bound_manifest["parent_freeze_manifest"] = {
        "path": str(freeze_manifest_path),
        "sha256": sha256(freeze_manifest_path),
    }
    bound_manifest_path.write_text(
        json.dumps(
            bound_manifest, ensure_ascii=False, indent=2, sort_keys=True
        )
        + "\n",
        encoding="utf-8",
    )
    return {
        "binding": str(binding_path),
        "freeze_manifest": str(bound_manifest_path),
        "record_count": len(binding),
        "complete": all_complete,
        "status_counts": bound_manifest["claim_binding_producer"][
            "status_counts"
        ],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("--frozen", required=True)
    prepare.add_argument("--freeze-manifest", required=True)
    prepare.add_argument("--out-dir", required=True)
    seal = subparsers.add_parser("seal")
    seal.add_argument("--frozen", required=True)
    seal.add_argument("--freeze-manifest", required=True)
    seal.add_argument("--review", required=True)
    seal.add_argument("--out-dir", required=True)
    args = parser.parse_args(argv)
    if args.command == "prepare":
        result = prepare_claim_binding_review(
            args.frozen, args.freeze_manifest, args.out_dir
        )
    else:
        result = seal_claim_binding_review(
            args.frozen, args.freeze_manifest, args.review, args.out_dir
        )
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
