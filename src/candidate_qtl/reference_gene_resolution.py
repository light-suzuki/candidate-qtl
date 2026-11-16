"""Reference-gene resolution with typed UniProt entry outcomes (#655).

The established resolver remains in ``reference_gene_resolution_legacy``. This
facade preserves its API and rewrites only the accession/error taxonomy that
depends on the provider's concrete HTTP outcome.
"""

from __future__ import annotations

import json as _json
from collections import Counter as _Counter
from dataclasses import asdict as _asdict
from pathlib import Path as _Path

import pandas as _pd

from . import reference_gene_resolution_legacy as _legacy
from .uniprot_error_policy import (
    UniProtEntryEvent as _UniProtEntryEvent,
    capture_uniprot_entry_events as _capture_uniprot_entry_events,
    fetch_uniprot_entry,
)

for _name in dir(_legacy):
    if not _name.startswith("__"):
        globals()[_name] = getattr(_legacy, _name)

_ORIGINAL_REFERENCE_GENE_RESOLVE = _legacy.reference_gene_resolve
_ORIGINAL_CANDIDATE_RESOLVER = _legacy._resolve_accession_candidates
_legacy._fetch_uniprot = fetch_uniprot_entry
_fetch_uniprot = fetch_uniprot_entry


def _replace_review_row(
    frame: _pd.DataFrame,
    hypothesis_id: str,
    reason: str,
    action: str,
    *,
    preserve_existing: bool = False,
) -> _pd.DataFrame:
    if not frame.empty and not preserve_existing:
        frame = frame[frame["target"].astype(str) != hypothesis_id]
    addition = _pd.DataFrame(
        [
            {
                "reason": f"reference_gene_resolution_required:{reason}",
                "target": hypothesis_id,
                "required_next_action": action,
                "allowed_choices": (
                    "correct proposed_stable_id manually | rely on reference_gene_name resolution | discard the proposal"
                ),
            }
        ]
    )
    return _pd.concat([frame, addition], ignore_index=True)


def _apply_entry_outcomes(
    result: dict[str, object],
    events: list[_UniProtEntryEvent],
) -> dict[str, object]:
    if not events:
        return result

    ledger_path = _Path(str(result["ledger"]))
    ledger = _pd.read_csv(ledger_path, sep="\t", dtype=str).fillna("")
    resolved_path = _Path(str(result["mechanisms_resolved"]))
    resolved = _pd.read_csv(resolved_path, sep="\t", dtype=str).fillna("")
    review_path = _Path(str(result["review_required"]))
    reviews = _pd.read_csv(review_path, sep="\t", dtype=str).fillna("")
    invalid_hypothesis_ids: set[str] = set()

    for index, row in ledger.iterrows():
        provenance = _json.loads(str(row.get("retrieval_provenance", "")) or "{}")
        attempts = provenance.get("uniprot_entry_attempts", [])
        if not attempts:
            # A cache/model row did not make this network attempt. An error
            # from another alias or case must never rewrite its healthy record.
            continue
        attempt = attempts[-1]
        event_index = attempt["event_index"]
        if (not isinstance(event_index, int) or not 0 <= event_index < len(events)
                or {key: value for key, value in attempt.items() if key != "event_index"} != _asdict(events[event_index])):
            raise ValueError("reference resolution entry attempt provenance mismatch")
        event = events[event_index]
        hypothesis_id = str(row.get("hypothesis_id", ""))
        case_id = str(row.get("case_id", ""))
        shared_id = ("case_id" in ledger.columns and
                     ledger.loc[ledger["hypothesis_id"].astype(str).eq(hypothesis_id), "case_id"].nunique() > 1)
        if event.state in {"invalid_syntax", "invalid_http_400"}:
            ledger.at[index, "accession_resolution_status"] = event.state
            invalid_hypothesis_ids.add(hypothesis_id)
            ledger.at[index, "resolution_origin"] = "not_resolved"
            ledger.at[index, "resolved_stable_id"] = ""
            ledger.at[index, "resolved_source_ids"] = ""
            ledger.at[index, "record_found"] = "False"
            ledger.at[index, "function_source_backed"] = "False"
            ledger.at[index, "evaluation_status"] = "evaluated"
            ledger.at[index, "provider_error"] = ""
            ledger.at[index, "review_reason"] = "invalid_submitted_accession"
            if "hypothesis_id" in resolved.columns and "proposed_stable_id" in resolved.columns:
                target = resolved["hypothesis_id"].astype(str) == hypothesis_id
                if "case_id" in resolved.columns:
                    target &= resolved["case_id"].astype(str) == str(row.get("case_id", ""))
                resolved.loc[
                    target,
                    "proposed_stable_id",
                ] = ""
            reviews = _replace_review_row(
                reviews,
                hypothesis_id,
                "invalid_submitted_accession",
                (
                    "the submitted accession was syntactically invalid or rejected with HTTP 400; it was cleared before hypothesis verification because this is an input error, not provider unavailability"
                    + (f"; case_id={case_id}" if shared_id else "")
                ),
                preserve_existing=shared_id,
            )
        elif event.state == "not_found_http_404":
            ledger.at[index, "accession_resolution_status"] = event.state
            ledger.at[index, "resolution_origin"] = "not_resolved"
            ledger.at[index, "resolved_stable_id"] = ""
            ledger.at[index, "resolved_source_ids"] = ""
            ledger.at[index, "record_found"] = "False"
            ledger.at[index, "function_source_backed"] = "False"
            ledger.at[index, "evaluation_status"] = "evaluated"
            ledger.at[index, "provider_error"] = ""
            ledger.at[index, "review_reason"] = "accession_not_found_in_authoritative_database"
        elif event.state in {"rate_limited", "provider_unavailable", "response_parse_error"}:
            ledger.at[index, "accession_resolution_status"] = event.state
            ledger.at[index, "resolution_origin"] = "not_resolved"
            ledger.at[index, "evaluation_status"] = "not_evaluated"
            ledger.at[index, "review_reason"] = "accession_provider_unavailable"
            if not str(ledger.at[index, "provider_error"]):
                ledger.at[index, "provider_error"] = event.detail

    ledger.to_csv(ledger_path, sep="\t", index=False, lineterminator="\n")
    resolved.to_csv(resolved_path, sep="\t", index=False, lineterminator="\n")
    if not reviews.empty:
        reviews = reviews.drop_duplicates().sort_values(["reason", "target"], kind="stable")
    reviews.to_csv(review_path, sep="\t", index=False, lineterminator="\n")

    manifest_path = _Path(str(result["manifest"]))
    manifest = _json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict):
        raise ValueError("reference-gene-resolution manifest must be a JSON object")
    manifest["resolution_origin_counts"] = {
        origin: int((ledger["resolution_origin"] == origin).sum())
        for origin in _legacy.RESOLUTION_ORIGINS
    }
    manifest["evaluation_status_counts"] = {
        status: int((ledger["evaluation_status"] == status).sum())
        for status in _legacy.EVALUATION_STATUSES
    }
    manifest["not_evaluated_hypothesis_ids"] = sorted(
        ledger.loc[ledger["evaluation_status"] == "not_evaluated", "hypothesis_id"].astype(str)
    )
    manifest["review_required_count"] = int(len(reviews))
    manifest["uniprot_entry_outcome_counts"] = dict(
        sorted(_Counter(event.state for event in events).items())
    )
    manifest["invalid_submitted_accession_hypothesis_ids"] = sorted(invalid_hypothesis_ids)
    gates = manifest.setdefault("gates", [])
    if isinstance(gates, list):
        gates.append(
            "invalid UniProt syntax and HTTP 400 are evaluated input errors and are cleared before downstream verification; HTTP 404 is definitive not-found; only transient/provider failures are not_evaluated"
        )
    manifest_path.write_text(
        _json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    counts = manifest["resolution_origin_counts"]
    evaluations = manifest["evaluation_status_counts"]
    result["provider_resolved_count"] = counts["provider_resolved"]
    result["model_supplied_count"] = counts["model_supplied"]
    result["not_resolved_count"] = counts["not_resolved"]
    result["not_evaluated_count"] = evaluations["not_evaluated"]
    result["review_required_count"] = int(len(reviews))
    return result


def reference_gene_resolve(*args: object, **kwargs: object) -> dict[str, object]:
    """Run the existing resolver and preserve typed provider/input outcomes."""
    # The facade historically exposed these helpers for offline monkeypatching.
    # Synchronise them only for the call, then restore the retained module so a
    # completed patch context cannot leak a mock into later resolutions.
    previous_fetch = _legacy._fetch_uniprot
    previous_resolve = _legacy._resolve_accession
    previous_candidates = _legacy._resolve_accession_candidates
    previous_resolve_one = _legacy._resolve_one
    _legacy._fetch_uniprot = globals()["_fetch_uniprot"]
    _legacy._resolve_accession = globals()["_resolve_accession"]
    # Preserve a directly patched legacy helper (the established test seam).
    if previous_candidates is _ORIGINAL_CANDIDATE_RESOLVER:
        _legacy._resolve_accession_candidates = globals()["_resolve_accession_candidates"]
    try:
        with _capture_uniprot_entry_events() as events:
            def resolve_one_with_attempts(*row_args: object, **row_kwargs: object) -> dict[str, object]:
                start = len(events)
                row = previous_resolve_one(*row_args, **row_kwargs)
                provenance = _json.loads(str(row.get("retrieval_provenance", "")) or "{}")
                provenance["uniprot_entry_attempts"] = [
                    {"event_index": event_index, **_asdict(events[event_index])}
                    for event_index in range(start, len(events))
                ]
                row["retrieval_provenance"] = _json.dumps(provenance, sort_keys=True, separators=(",", ":"))
                return row

            _legacy._resolve_one = resolve_one_with_attempts
            result = _ORIGINAL_REFERENCE_GENE_RESOLVE(*args, **kwargs)
    finally:
        _legacy._fetch_uniprot = previous_fetch
        _legacy._resolve_accession = previous_resolve
        _legacy._resolve_accession_candidates = previous_candidates
        _legacy._resolve_one = previous_resolve_one
    return _apply_entry_outcomes(result, events)


_legacy.reference_gene_resolve = reference_gene_resolve
