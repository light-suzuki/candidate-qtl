"""Deterministic, provenance-complete Route-C projection reduction.

The legacy projection implementation kept one mutable mapping record per focal
gene, so a later BLAST row could overwrite an exact or orthology match. This
module records every exact/orthology/BLAST observation in an append-only ledger
and reduces that ledger with a versioned, order-invariant policy (issue #663).
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Iterable

import pandas as pd

from . import __version__
from .core import file_fingerprint, runtime_path, sha256
from .route_c_identity_mapping import IDENTITY_MAPPING_CONTRACT, PAPER_MAPPING_CONTRACT, load_identity_mapping, mapping_links_binding
from .route_c_direct_rescue import (
    PROJECTION_COLUMNS,
    direct_record_can_seed_projection,
    validate_direct_rescue_frame,
    validate_projection_frame,
)

SCHEMA_VERSION = "1.0"
PROJECTION_EVIDENCE_CONTRACT = "route-c-direct-projection-evidence-1.0.0"
PROJECTION_SELECTION_POLICY = "route-c-direct-projection-reducer-1.0.0"

PROJECTION_EVIDENCE_COLUMNS = (
    "record_id",
    "qtl_id",
    "focal_gene_id",
    "mapping_method",
    "source_identity_field",
    "source_identity_value",
    "source_identity_verification_state",
    "method_evaluation_state",
    "method_evidence_id",
    "identity",
    "query_coverage",
    "target_coverage",
    "relationship",
    "conflict_state",
    "interval_universe_fingerprint",
    "evidence_sha256",
    "rank_contribution",
    "notes",
)

_MAPPING_METHODS = frozenset({"exact_focal_id", "orthology_reference_gene", "blastp"})
_METHOD_EVALUATION_STATES = frozenset({
    "evaluated_with_hit",
    "evaluated_no_hit",
    "not_evaluated_not_supplied",
    "not_evaluated_no_panel",
    "not_evaluated_missing_atlas_column",
    "not_evaluated_record_ineligible",
    "not_evaluated_reference_bridge",
})
_CONFLICT_STATES = frozenset({
    "none",
    "corroborating_methods",
    "one_to_many_consistent",
    "mapping_conflict",
})
_METHOD_PRECEDENCE = {
    "exact_focal_id": 0,
    "orthology_reference_gene": 1,
    "blastp": 2,
}
_SOURCE_FIELD_PRECEDENCE = {
    "reported_gene_id": 0,
    "reported_gene_symbol": 1,
    "reference_sequence_id": 2,
    "panel_id": 3,
    "": 9,
}
_SOURCE_IDENTITY_FIELDS = (
    "reported_gene_id",
    "reported_gene_symbol",
    "reference_sequence_id",
)


def _blast_source_identity(
    identities: tuple[tuple[str, str], ...],
    record: dict[str, object],
) -> tuple[str, str]:
    """Select the source endpoint recorded on BLAST projection evidence.

    Locus-only Direct-C records intentionally do not enter the exact-ID or
    orthology joins.  When such a record is projected by an eligible BLAST
    hit, retain its submitted locus as provenance instead of manufacturing the
    record ID as ``panel_id``.  The fallback is scoped to BLAST evidence so it
    cannot broaden exact/orthology matching.
    """
    sequence = next(
        ((field, value) for field, value in identities
         if field == "reference_sequence_id"),
        None,
    )
    if sequence:
        return sequence
    if identities:
        return identities[0]
    locus = str(record.get("reported_locus", "") or "").strip()
    if locus:
        return "reported_locus", locus
    return "panel_id", str(record.get("record_id", ""))


def _canonical_json(payload: object) -> str:
    return json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _text_sha256(payload: object) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _split_values(value: object) -> tuple[str, ...]:
    values = {
        token.strip()
        for token in str(value or "").replace(",", ";").replace("|", ";").split(";")
        if token.strip()
    }
    return tuple(sorted(values, key=lambda item: (item.casefold(), item)))


def _identity_values(record: dict[str, object]) -> tuple[tuple[str, str], ...]:
    rows: list[tuple[str, str]] = []
    for field in _SOURCE_IDENTITY_FIELDS:
        rows.extend((field, value) for value in _split_values(record.get(field, "")))
    return tuple(rows)


def _identity_state(record: dict[str, object], field: str) -> str:
    state = str(record.get("gene_identity_state", "") or "").strip()
    if field == "reference_sequence_id":
        # The frozen schema does not separately prove that a sequence identifier
        # is bound to the source claim. The downstream #661 sidecar remains the
        # support gate even when the sequence projection itself is useful.
        return "not_evaluated_specific_sequence_binding"
    return state or "not_evaluated"


def _canonical_frame_sha256(frame: pd.DataFrame, *, sort_columns: tuple[str, ...]) -> str:
    canonical = frame.copy().fillna("")
    canonical = canonical.sort_values(
        list(sort_columns),
        kind="stable",
        key=lambda series: series.astype(str).str.casefold(),
    ).reset_index(drop=True)
    payload = canonical.to_csv(sep="\t", index=False, lineterminator="\n").encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _number(value: object) -> float:
    try:
        number = float(str(value).strip())
    except (TypeError, ValueError):
        return -math.inf
    return number if math.isfinite(number) else -math.inf


def _evidence_id(prefix: str, payload: object, ordinal: int | None = None) -> str:
    digest = _text_sha256(payload)[:20]
    return f"{prefix}:{digest}" if ordinal is None else f"{prefix}:{digest}:{ordinal:04d}"


def _base_evidence_row(
    *,
    record_id: str,
    qtl_id: str,
    focal_gene_id: str,
    mapping_method: str,
    source_identity_field: str,
    source_identity_value: str,
    source_identity_verification_state: str,
    method_evaluation_state: str,
    method_evidence_id: str,
    identity: str = "",
    query_coverage: str = "",
    target_coverage: str = "",
    relationship: str = "",
    interval_universe_fingerprint: str,
    notes: str = "",
) -> dict[str, str]:
    return {
        "record_id": record_id,
        "qtl_id": qtl_id,
        "focal_gene_id": focal_gene_id,
        "mapping_method": mapping_method,
        "source_identity_field": source_identity_field,
        "source_identity_value": source_identity_value,
        "source_identity_verification_state": source_identity_verification_state,
        "method_evaluation_state": method_evaluation_state,
        "method_evidence_id": method_evidence_id,
        "identity": identity,
        "query_coverage": query_coverage,
        "target_coverage": target_coverage,
        "relationship": relationship,
        "conflict_state": "none",
        "interval_universe_fingerprint": interval_universe_fingerprint,
        "evidence_sha256": "",
        "rank_contribution": "none",
        "notes": notes,
    }


def _status_row(
    *,
    record_id: str,
    qtl_id: str,
    method: str,
    state: str,
    interval_hash: str,
    note: str,
) -> dict[str, str]:
    payload = {"record_id": record_id, "method": method, "state": state, "note": note}
    return _base_evidence_row(
        record_id=record_id,
        qtl_id=qtl_id,
        focal_gene_id="",
        mapping_method=method,
        source_identity_field="",
        source_identity_value="",
        source_identity_verification_state="not_evaluated",
        method_evaluation_state=state,
        method_evidence_id=_evidence_id("status", payload),
        relationship="not_evaluated" if state.startswith("not_evaluated") else "no_interval_match",
        interval_universe_fingerprint=interval_hash,
        notes=note,
    )


def _finalize_evidence_rows(rows: list[dict[str, str]]) -> pd.DataFrame:
    by_record: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_record.setdefault(str(row["record_id"]), []).append(row)

    for record_rows in by_record.values():
        method_targets: dict[str, set[str]] = {}
        for row in record_rows:
            focal = str(row.get("focal_gene_id", "")).strip()
            if focal and row.get("method_evaluation_state") == "evaluated_with_hit":
                method_targets.setdefault(str(row["mapping_method"]), set()).add(focal)
        all_targets = set().union(*method_targets.values()) if method_targets else set()
        nonempty_sets = [tuple(sorted(targets)) for targets in method_targets.values() if targets]
        if len(all_targets) == 1 and len(method_targets) > 1:
            conflict_state = "corroborating_methods"
        elif len(all_targets) > 1:
            conflict_state = (
                "one_to_many_consistent"
                if nonempty_sets and len(set(nonempty_sets)) == 1
                else "mapping_conflict"
            )
        else:
            conflict_state = "none"
        for row in record_rows:
            row["conflict_state"] = conflict_state
            digest_payload = {key: value for key, value in row.items() if key != "evidence_sha256"}
            row["evidence_sha256"] = _text_sha256(digest_payload)

    frame = pd.DataFrame(rows, columns=PROJECTION_EVIDENCE_COLUMNS).fillna("")
    frame = frame.sort_values(
        [
            "record_id",
            "focal_gene_id",
            "mapping_method",
            "source_identity_field",
            "source_identity_value",
            "method_evidence_id",
        ],
        kind="stable",
        key=lambda series: series.astype(str).str.casefold(),
    ).reset_index(drop=True)
    return validate_projection_evidence_frame(frame)


def validate_projection_evidence_frame(frame: pd.DataFrame) -> pd.DataFrame:
    missing = sorted(set(PROJECTION_EVIDENCE_COLUMNS) - set(frame.columns))
    extra = sorted(set(frame.columns) - set(PROJECTION_EVIDENCE_COLUMNS))
    if missing or extra:
        parts = []
        if missing:
            parts.append(f"missing: {', '.join(missing)}")
        if extra:
            parts.append(f"outside the contract: {', '.join(extra)}")
        raise ValueError(f"Route C projection evidence has {'; '.join(parts)}")
    normalized = frame[list(PROJECTION_EVIDENCE_COLUMNS)].copy().fillna("")
    for column in PROJECTION_EVIDENCE_COLUMNS:
        normalized[column] = normalized[column].astype(str).str.strip()
    if normalized.empty:
        raise ValueError("Route C projection evidence is empty")
    required = (
        "record_id",
        "mapping_method",
        "method_evaluation_state",
        "method_evidence_id",
        "conflict_state",
        "interval_universe_fingerprint",
        "evidence_sha256",
        "rank_contribution",
    )
    for column in required:
        if normalized[column].eq("").any():
            raise ValueError(f"Route C projection evidence requires non-empty {column}")
    invalid_methods = sorted(set(normalized["mapping_method"]) - _MAPPING_METHODS)
    if invalid_methods:
        raise ValueError(f"Route C projection evidence has invalid mapping_method: {', '.join(invalid_methods)}")
    invalid_states = sorted(set(normalized["method_evaluation_state"]) - _METHOD_EVALUATION_STATES)
    if invalid_states:
        raise ValueError(f"Route C projection evidence has invalid method_evaluation_state: {', '.join(invalid_states)}")
    invalid_conflicts = sorted(set(normalized["conflict_state"]) - _CONFLICT_STATES)
    if invalid_conflicts:
        raise ValueError(f"Route C projection evidence has invalid conflict_state: {', '.join(invalid_conflicts)}")
    if not normalized["rank_contribution"].eq("none").all():
        raise ValueError("Route C projection evidence rank_contribution must always be 'none'")
    if normalized["method_evidence_id"].duplicated().any():
        raise ValueError("Route C projection evidence method_evidence_id must be unique")
    for row in normalized.to_dict("records"):
        has_hit = row["method_evaluation_state"] == "evaluated_with_hit"
        if has_hit != bool(row["focal_gene_id"]):
            raise ValueError("evaluated_with_hit and focal_gene_id must occur together")
        expected = _text_sha256({key: value for key, value in row.items() if key != "evidence_sha256"})
        if row["evidence_sha256"].removeprefix("sha256:") != expected:
            raise ValueError("Route C projection evidence row hash mismatch")
    return normalized


def _selected_evidence(rows: Iterable[dict[str, str]]) -> dict[str, str]:
    def key(row: dict[str, str]) -> tuple[object, ...]:
        method = str(row["mapping_method"])
        if method == "blastp":
            metric_key = (
                -_number(row.get("identity", "")),
                -_number(row.get("query_coverage", "")),
                -_number(row.get("target_coverage", "")),
            )
        else:
            metric_key = (0.0, 0.0, 0.0)
        return (
            _METHOD_PRECEDENCE[method],
            *metric_key,
            _SOURCE_FIELD_PRECEDENCE.get(str(row.get("source_identity_field", "")), 8),
            str(row.get("source_identity_value", "")).casefold(),
            str(row["method_evidence_id"]),
        )

    return min(rows, key=key)


def _projection_rows_for_record(
    record: dict[str, object],
    evidence_rows: list[dict[str, str]],
    qtl_id: str,
    interval_hash: str,
) -> list[dict[str, str]]:
    record_id = str(record["record_id"]).strip()
    if not direct_record_can_seed_projection(record):
        return [{
            "record_id": record_id,
            "qtl_id": "",
            "focal_gene_id": "",
            "projection_state": "not_projected",
            "projection_basis": "not_evaluated",
            "source_reference_gene_id": str(record.get("reported_gene_id", "")),
            "identity": "",
            "query_coverage": "",
            "target_coverage": "",
            "ambiguity_state": "not_evaluated",
            "interval_universe_fingerprint": interval_hash,
            "rank_contribution": "none",
            "notes": "direct record is not projection-eligible",
        }]

    real = [
        row for row in evidence_rows
        if row["method_evaluation_state"] == "evaluated_with_hit" and row["focal_gene_id"]
    ]
    if not real:
        states = {row["mapping_method"]: row["method_evaluation_state"] for row in evidence_rows}
        bridge_unavailable = states.get("blastp") == "not_evaluated_reference_bridge"
        if states.get("blastp") == "evaluated_no_hit":
            basis = "blastp"
        elif states.get("orthology_reference_gene") == "evaluated_no_hit":
            basis = "orthology_reference_gene"
        else:
            basis = "exact_focal_id"
        status_text = ";".join(
            f"{method}={states.get(method, 'not_evaluated')}"
            for method in sorted(_MAPPING_METHODS, key=_METHOD_PRECEDENCE.get)
        )
        return [{
            "record_id": record_id,
            "qtl_id": qtl_id,
            "focal_gene_id": "",
            "projection_state": "not_projected" if bridge_unavailable else "no_interval_match",
            "projection_basis": "not_evaluated" if bridge_unavailable else basis,
            "source_reference_gene_id": str(record.get("reported_gene_id", "")),
            "identity": "",
            "query_coverage": "",
            "target_coverage": "",
            "ambiguity_state": "not_evaluated" if bridge_unavailable else "no_threshold_hit",
            "interval_universe_fingerprint": interval_hash,
            "rank_contribution": "none",
            "notes": (
                ("reference sequence retrieval/search was not evaluated; no biological absence inferred; " if bridge_unavailable else
                 "evaluated interval search found no eligible match; outside-QTL presence is not inferred; ") +
                f"method_states={status_text}; selection_policy={PROJECTION_SELECTION_POLICY}"
            ),
        }]

    targets = sorted({row["focal_gene_id"] for row in real}, key=lambda item: (item.casefold(), item))
    conflict_state = str(real[0]["conflict_state"])
    ambiguity = "unique" if len(targets) == 1 else "one_to_many"
    all_methods = sorted({row["mapping_method"] for row in real}, key=_METHOD_PRECEDENCE.get)
    rows: list[dict[str, str]] = []
    for target in targets:
        target_evidence = [row for row in real if row["focal_gene_id"] == target]
        selected = _selected_evidence(target_evidence)
        target_methods = sorted({row["mapping_method"] for row in target_evidence}, key=_METHOD_PRECEDENCE.get)
        evidence_ids = ";".join(sorted(row["method_evidence_id"] for row in target_evidence))
        rows.append({
            "record_id": record_id,
            "qtl_id": qtl_id,
            "focal_gene_id": target,
            "projection_state": "inside_current_qtl",
            "projection_basis": selected["mapping_method"],
            "source_reference_gene_id": str(record.get("reported_gene_id", "")),
            "identity": selected["identity"],
            "query_coverage": selected["query_coverage"],
            "target_coverage": selected["target_coverage"],
            "ambiguity_state": ambiguity,
            "interval_universe_fingerprint": interval_hash,
            "rank_contribution": "none",
            "notes": (
                "projection evidence only; cross-species trait support is not automatically transferred; "
                f"selected_projection_basis={selected['mapping_method']}; "
                f"target_projection_bases={';'.join(target_methods)}; "
                f"all_projection_bases={';'.join(all_methods)}; "
                f"conflict_state={conflict_state}; method_evidence_ids={evidence_ids}; "
                f"selection_policy={PROJECTION_SELECTION_POLICY}"
            ),
        })
    return rows


def _collect_record_evidence(
    record: dict[str, object],
    atlas_rows: list[dict[str, str]],
    atlas_gene_ids: set[str],
    atlas_has_orthology: bool,
    hits: pd.DataFrame | None,
    qtl_id: str,
    interval_hash: str,
) -> list[dict[str, str]]:
    record_id = str(record["record_id"]).strip()
    if not direct_record_can_seed_projection(record):
        return [
            _status_row(
                record_id=record_id,
                qtl_id=qtl_id,
                method=method,
                state="not_evaluated_record_ineligible",
                interval_hash=interval_hash,
                note="direct record is not projection-eligible",
            )
            for method in sorted(_MAPPING_METHODS, key=_METHOD_PRECEDENCE.get)
        ]

    identities = _identity_values(record)
    rows: list[dict[str, str]] = []

    exact_count = 0
    orthology_count = 0
    for atlas_row in atlas_rows:
        gene_id = str(atlas_row["gene_id"]).strip()
        gene_key = gene_id.casefold()
        for field, value in identities:
            if gene_key == value.casefold():
                payload = {"record_id": record_id, "gene_id": gene_id, "field": field, "value": value}
                rows.append(_base_evidence_row(
                    record_id=record_id,
                    qtl_id=qtl_id,
                    focal_gene_id=gene_id,
                    mapping_method="exact_focal_id",
                    source_identity_field=field,
                    source_identity_value=value,
                    source_identity_verification_state=_identity_state(record, field),
                    method_evaluation_state="evaluated_with_hit",
                    method_evidence_id=_evidence_id("exact", payload),
                    relationship="exact_identifier_match",
                    interval_universe_fingerprint=interval_hash,
                ))
                exact_count += 1
        if atlas_has_orthology:
            references = {
                value.casefold()
                for value in _split_values(atlas_row.get("orthology_reference_gene", ""))
            }
            for field, value in identities:
                if value.casefold() in references:
                    payload = {"record_id": record_id, "gene_id": gene_id, "field": field, "value": value}
                    rows.append(_base_evidence_row(
                        record_id=record_id,
                        qtl_id=qtl_id,
                        focal_gene_id=gene_id,
                        mapping_method="orthology_reference_gene",
                        source_identity_field=field,
                        source_identity_value=value,
                        source_identity_verification_state=_identity_state(record, field),
                        method_evaluation_state="evaluated_with_hit",
                        method_evidence_id=_evidence_id("orthology", payload),
                        relationship="atlas_orthology_reference_match",
                        interval_universe_fingerprint=interval_hash,
                    ))
                    orthology_count += 1

    if exact_count == 0:
        rows.append(_status_row(
            record_id=record_id,
            qtl_id=qtl_id,
            method="exact_focal_id",
            state="evaluated_no_hit",
            interval_hash=interval_hash,
            note="interval atlas gene identifiers were evaluated",
        ))
    if not atlas_has_orthology:
        rows.append(_status_row(
            record_id=record_id,
            qtl_id=qtl_id,
            method="orthology_reference_gene",
            state="not_evaluated_missing_atlas_column",
            interval_hash=interval_hash,
            note="interval atlas has no orthology_reference_gene column",
        ))
    elif orthology_count == 0:
        rows.append(_status_row(
            record_id=record_id,
            qtl_id=qtl_id,
            method="orthology_reference_gene",
            state="evaluated_no_hit",
            interval_hash=interval_hash,
            note="orthology_reference_gene values were evaluated",
        ))

    if hits is None:
        rows.append(_status_row(
            record_id=record_id,
            qtl_id=qtl_id,
            method="blastp",
            state="not_evaluated_not_supplied",
            interval_hash=interval_hash,
            note="eligible interval hits artifact was not supplied",
        ))
    else:
        panel = hits[hits["panel_id"].astype(str).eq(record_id)].copy()
        if panel.empty:
            rows.append(_status_row(
                record_id=record_id,
                qtl_id=qtl_id,
                method="blastp",
                state="not_evaluated_no_panel",
                interval_hash=interval_hash,
                note="hits artifact was supplied but this record has no panel rows",
            ))
        else:
            eligible = panel[
                panel["eligibility"].astype(str).isin(
                    ("eligible", "eligible_interval_hit")
                )
                & panel["target_id"].astype(str).isin(atlas_gene_ids)
            ].copy()
            if eligible.empty:
                unavailable = ("provider_status" in panel and panel["provider_status"].astype(str).str.startswith("not_evaluated").all())
                rows.append(_status_row(
                    record_id=record_id,
                    qtl_id=qtl_id,
                    method="blastp",
                    state="not_evaluated_reference_bridge" if unavailable else "evaluated_no_hit",
                    interval_hash=interval_hash,
                    note="reference bridge not evaluated; no biological absence inferred" if unavailable else "this record had BLAST panel rows but no eligible target in the current interval atlas",
                ))
            else:
                raw_rows = eligible.to_dict("records")
                raw_rows.sort(
                    key=lambda row: _canonical_json({str(key): str(value) for key, value in row.items()})
                )
                duplicate_counts: dict[str, int] = {}
                source_field, source_value = _blast_source_identity(identities, record)
                for hit in raw_rows:
                    normalized_hit = {str(key): str(value) for key, value in sorted(hit.items())}
                    hit_digest = _text_sha256(normalized_hit)
                    duplicate_counts[hit_digest] = duplicate_counts.get(hit_digest, 0) + 1
                    target = str(hit["target_id"]).strip()
                    rows.append(_base_evidence_row(
                        record_id=record_id,
                        qtl_id=qtl_id,
                        focal_gene_id=target,
                        mapping_method="blastp",
                        source_identity_field=source_field,
                        source_identity_value=source_value,
                        source_identity_verification_state=_identity_state(record, source_field),
                        method_evaluation_state="evaluated_with_hit",
                        method_evidence_id=_evidence_id(
                            "blastp", normalized_hit, duplicate_counts[hit_digest]
                        ),
                        identity=str(hit.get("pident", "")),
                        query_coverage=str(hit.get("query_coverage", "")),
                        target_coverage=str(hit.get("target_coverage", "")),
                        relationship="threshold_gated_sequence_similarity",
                        interval_universe_fingerprint=interval_hash,
                        notes=f"raw_hit_sha256={hit_digest}",
                    ))
    return rows


def direct_rescue_project(
    frozen_arg: str,
    freeze_manifest_arg: str,
    atlas_arg: str,
    qtl_id: str,
    out_dir_arg: str,
    eligible_hits_arg: str | None = None,
    full_hash: bool = False,
    identity_mapping_arg: str | None = None,
    reference_bridge_arg: str | None = None,
) -> dict[str, object]:
    """Project frozen Direct-C records with an append-only evidence ledger."""
    if not str(qtl_id).strip():
        raise ValueError("qtl_id is required")
    if reference_bridge_arg and eligible_hits_arg:
        raise ValueError("reference bridge and raw eligible hits require separate projection generations")
    frozen_path = runtime_path(frozen_arg)
    freeze_manifest_path = runtime_path(freeze_manifest_arg)
    freeze_manifest = json.loads(freeze_manifest_path.read_text(encoding="utf-8"))
    expected_frozen_sha = str(freeze_manifest.get("frozen_sha256", "")).removeprefix("sha256:")
    if not expected_frozen_sha or expected_frozen_sha != sha256(frozen_path):
        raise ValueError("Route C frozen artifact hash does not match freeze manifest")

    atlas_path = runtime_path(atlas_arg)
    direct = validate_direct_rescue_frame(
        pd.read_csv(frozen_path, sep="\t", dtype=str, keep_default_na=False)
    )
    atlas = pd.read_csv(atlas_path, sep="\t", dtype=str, keep_default_na=False)
    if "gene_id" not in atlas.columns:
        raise ValueError("Route C projection atlas requires gene_id")
    if atlas["gene_id"].astype(str).str.strip().eq("").any():
        raise ValueError("Route C projection atlas requires non-empty gene_id")
    if atlas["gene_id"].astype(str).str.strip().duplicated().any():
        raise ValueError("Route C projection atlas gene_id must be unique")
    atlas = atlas.sort_values(
        "gene_id",
        kind="stable",
        key=lambda series: series.astype(str).str.casefold(),
    ).reset_index(drop=True)
    atlas_rows = [
        {str(key): str(value) for key, value in row.items()}
        for row in atlas.to_dict("records")
    ]
    atlas_gene_ids = set(atlas["gene_id"].astype(str))
    atlas_has_orthology = "orthology_reference_gene" in atlas.columns
    hits_path: Path | None = runtime_path(eligible_hits_arg) if eligible_hits_arg else None
    mapping_path = runtime_path(identity_mapping_arg) if identity_mapping_arg else None
    identity_mapping = load_identity_mapping(
        mapping_path, {str(row["record_id"]): row for row in direct.to_dict("records")},
        freeze_manifest_path=freeze_manifest_path, reference_bridge_arg=reference_bridge_arg,
        eligible_hits=hits_path,
    )
    provider_edges = "provider_accession" in identity_mapping
    from .route_c_provider_bait import CONTRACT as PROVIDER_EDGE_CONTRACT, EDGE_STATE as PROVIDER_EDGE_STATE
    from .route_c_paper_identity_mapping import paper_mapping_records
    projection_records = paper_mapping_records(
        {str(row["record_id"]): row for row in direct.to_dict("records")}, identity_mapping)

    hits: pd.DataFrame | None = None
    if hits_path is not None:
        hits = pd.read_csv(hits_path, sep="\t", dtype=str, keep_default_na=False).fillna("")
        missing = sorted(
            {"panel_id", "target_id", "eligibility", "pident", "query_coverage", "target_coverage"}
            - set(hits.columns)
        )
        if missing:
            raise ValueError(f"eligible interval hits are missing columns: {', '.join(missing)}")

    bridge_endpoints = {}
    if reference_bridge_arg:
        from .reference_bridge_consumers import direct_bridge_hits
        bridge_hits, bridge_endpoints = direct_bridge_hits(reference_bridge_arg,
            {str(row["record_id"]): row for row in direct.to_dict("records")}, identity_mapping)
        # Exact/orthology observations remain separate in the existing reducer.
        hits = bridge_hits

    interval_hash = f"sha256:{_canonical_frame_sha256(atlas, sort_columns=('gene_id',))}"
    all_evidence_rows: list[dict[str, str]] = []
    for record in sorted(projection_records.values(), key=lambda row: str(row["record_id"]).casefold()):
        all_evidence_rows.extend(
            _collect_record_evidence(
                record,
                atlas_rows,
                atlas_gene_ids,
                atlas_has_orthology,
                hits,
                str(qtl_id),
                interval_hash,
            )
        )

    for observation in all_evidence_rows:
        if observation["mapping_method"] == "blastp" and observation["record_id"] in bridge_endpoints:
            # Use a structured CDS qualifier or the exact frozen accession;
            # neither creates a name alias or authorizes trait support.
            field, value = bridge_endpoints[observation["record_id"]]
            observation["source_identity_field"] = field
            observation["source_identity_value"] = value
            observation["source_identity_verification_state"] = "sequence_correspondence_review_only"
        if provider_edges:
            # The BLAST query was the provider entry that cites this record's
            # paper under its symbol (rechecked on load). Record that symbol as
            # the hit's source endpoint; this is review provenance only.
            edge = identity_mapping[identity_mapping["record_id"] == observation["record_id"]]
            if observation["mapping_method"] == "blastp" and not edge.empty:
                observation["source_identity_field"] = "reported_gene_symbol"
                observation["source_identity_value"] = str(edge.iloc[0]["source_identity_value"])
                observation["source_identity_verification_state"] = PROVIDER_EDGE_STATE
            continue
        if observation["source_identity_field"] == "reference_sequence_id":
            mapping_rows = identity_mapping[identity_mapping["record_id"] == observation["record_id"]]
            if any(mapping_links_binding(mapping_rows, observation["record_id"], row["source_identity_value"], observation)
                   for row in mapping_rows.to_dict("records")):
                observation["source_identity_verification_state"] = (
                    "reviewed_paper_sequence_mapping" if "evidence_type" in mapping_rows
                    else "verified_structured_sequence_mapping")
    evidence = _finalize_evidence_rows(all_evidence_rows)
    finalized_by_record = {
        str(record_id): group.to_dict("records")
        for record_id, group in evidence.groupby("record_id", sort=False)
    }
    projection_rows: list[dict[str, str]] = []
    for record in sorted(projection_records.values(), key=lambda row: str(row["record_id"]).casefold()):
        projection_rows.extend(
            _projection_rows_for_record(
                record,
                finalized_by_record[str(record["record_id"])],
                str(qtl_id),
                interval_hash,
            )
        )
    projection = validate_projection_frame(
        pd.DataFrame(projection_rows, columns=PROJECTION_COLUMNS)
    )
    projection = projection.sort_values(
        ["record_id", "focal_gene_id", "projection_basis"],
        kind="stable",
        key=lambda series: series.astype(str).str.casefold(),
    ).reset_index(drop=True)

    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    projection_path = out / "route_c_direct_projection.tsv"
    evidence_ledger_path = out / "route_c_direct_projection_evidence.tsv"
    evidence_chain_path = out / "route_c_direct_evidence_chain.tsv"
    projection.to_csv(projection_path, sep="\t", index=False, lineterminator="\n")
    evidence.to_csv(evidence_ledger_path, sep="\t", index=False, lineterminator="\n")
    direct_sorted = direct.sort_values(
        "record_id", kind="stable", key=lambda series: series.astype(str).str.casefold()
    )
    evidence_chain = direct_sorted.merge(
        projection, on="record_id", how="left", suffixes=("_direct", "_projection")
    )
    evidence_chain.to_csv(evidence_chain_path, sep="\t", index=False, lineterminator="\n")

    conflict_ids = sorted(
        set(
            evidence.loc[
                evidence["conflict_state"].eq("mapping_conflict"), "record_id"
            ].astype(str)
        ),
        key=str.casefold,
    )
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "command": "route-c-direct-project",
        "candidate_qtl_version": __version__,
        "retrieval_frozen_before_interval": True,
        "frozen_sha256": sha256(frozen_path),
        "qtl_id": str(qtl_id),
        "interval_universe_fingerprint": interval_hash,
        "projection_sha256": sha256(projection_path),
        "evidence_chain_sha256": sha256(evidence_chain_path),
        "projection_evidence_contract": PROJECTION_EVIDENCE_CONTRACT,
        "projection_evidence_sha256": sha256(evidence_ledger_path),
        "projection_selection_policy": PROJECTION_SELECTION_POLICY,
        "projection_state_counts": {
            str(key): int(value)
            for key, value in projection["projection_state"].value_counts().items()
        },
        "projection_basis_counts": {
            str(key): int(value)
            for key, value in projection["projection_basis"].value_counts().items()
        },
        "method_evaluation_state_counts": {
            str(key): int(value)
            for key, value in evidence["method_evaluation_state"].value_counts().items()
        },
        "conflict_state_counts": {
            str(key): int(value)
            for key, value in evidence.groupby("record_id", sort=False)["conflict_state"]
            .first()
            .value_counts()
            .items()
        },
        "mapping_conflict_record_ids": conflict_ids,
        "mapping_conflict_count": len(conflict_ids),
        "rank_contribution": "none",
        "trait_transfer_policy": "projection/homology is not automatic trait-specific support",
        "inputs": {
            "frozen": file_fingerprint(frozen_path, full_hash),
            "freeze_manifest": file_fingerprint(freeze_manifest_path, full_hash or "evidence_type" in identity_mapping),
            "atlas": file_fingerprint(atlas_path, full_hash),
            "eligible_hits": file_fingerprint(hits_path, full_hash) if hits_path else None,
            "identity_mapping": {**file_fingerprint(mapping_path, True),
                                 "contract": PAPER_MAPPING_CONTRACT if "evidence_type" in identity_mapping
                                     else PROVIDER_EDGE_CONTRACT if provider_edges
                                     else IDENTITY_MAPPING_CONTRACT} if mapping_path else None,
            "reference_bridge": file_fingerprint(runtime_path(reference_bridge_arg), True) if reference_bridge_arg else None,
        },
        "outputs": {
            "projection": str(projection_path),
            "projection_evidence": str(evidence_ledger_path),
            "evidence_chain": str(evidence_chain_path),
        },
    }
    manifest_path = out / "route_c_direct_projection_manifest.json"
    manifest["outputs"]["manifest"] = str(manifest_path)
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        "projection": str(projection_path),
        "projection_evidence": str(evidence_ledger_path),
        "evidence_chain": str(evidence_chain_path),
        "manifest": str(manifest_path),
        "inside_count": int(
            projection["projection_state"].eq("inside_current_qtl").sum()
        ),
        "mapping_conflict_count": len(conflict_ids),
    }
