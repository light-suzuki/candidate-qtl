"""Route-A hypothesis utilities with typed UniProt entry outcomes (#655).

The historical implementation remains in ``llm_pathway_legacy``. This facade
preserves its public API while replacing only the UniProt entry fetch and the
small amount of output taxonomy that depends on HTTP/input outcomes.
"""

from __future__ import annotations

import json as _json
import hashlib as _hashlib
from collections import Counter as _Counter
from pathlib import Path as _Path

import pandas as _pd

from . import llm_pathway_legacy as _legacy
from .identifiers import canonical_accession as _canonical_accession
from .uniprot_error_policy import (
    UniProtEntryEvent as _UniProtEntryEvent,
    fetch_uniprot_entry,
)

# Preserve the complete compatibility surface, including the low-level helpers
# intentionally used by offline tests and reference-gene resolution.
for _name in dir(_legacy):
    if not _name.startswith("__"):
        globals()[_name] = getattr(_legacy, _name)

_ORIGINAL_HYPOTHESIS_VERIFY = _legacy.hypothesis_verify
_legacy._fetch_uniprot = fetch_uniprot_entry
_fetch_uniprot = fetch_uniprot_entry


def _event_key(value: object) -> str:
    raw = str(value or "").strip()
    canonical = _canonical_accession(raw)
    return (canonical or raw).casefold()


def _rewrite_review_ledger(
    path: _Path,
    *,
    invalid_hypothesis_ids: set[str],
) -> None:
    if not path.exists():
        return
    frame = _pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if frame.empty:
        columns = list(frame.columns)
    else:
        columns = list(frame.columns)
        frame = frame[
            ~(
                frame["target"].astype(str).isin(invalid_hypothesis_ids)
                & frame["reason"].astype(str).str.startswith("proposal_rejected:")
            )
        ]
    additions = [
        {
            "reason": "proposal_rejected:invalid_submitted_accession",
            "target": hypothesis_id,
            "required_next_action": (
                "correct the submitted UniProt accession; HTTP 400 or invalid syntax is input error, not provider unavailability"
            ),
            "allowed_choices": "correct accession | propose another reference gene | keep rejected",
        }
        for hypothesis_id in sorted(invalid_hypothesis_ids)
    ]
    if additions:
        frame = _pd.concat([frame, _pd.DataFrame(additions)], ignore_index=True)
    if columns:
        frame = frame.reindex(columns=columns)
    if not frame.empty:
        frame = frame.drop_duplicates().sort_values(["reason", "target"], kind="stable")
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


def _apply_entry_outcomes(
    result: dict[str, object],
    events: list[_UniProtEntryEvent],
    *,
    events_by_hypothesis: dict[tuple[str, str], list[_UniProtEntryEvent]] | None = None,
) -> dict[str, object]:
    if not events:
        return result
    validation_path = _Path(str(result["validation"]))
    validation = _pd.read_csv(validation_path, sep="\t", dtype=str).fillna("")
    invalid_hypothesis_ids: set[str] = set()

    for index, row in validation.iterrows():
        row_key = (str(row.get("case_id", "")), str(row.get("hypothesis_id", "")))
        row_events = (events_by_hypothesis or {}).get(row_key, [])
        if not row_events:
            # Cache rows and other hypotheses did not make this request.
            # Unbound events contribute audit counts only, never row decisions.
            continue
        event = row_events[-1]
        if _event_key(event.accession) != _event_key(row.get("proposed_stable_id", "")):
            raise ValueError("hypothesis verification entry attempt does not match its accession")
        if event.state in {"invalid_syntax", "invalid_http_400"}:
            validation.at[index, "accession_resolution_status"] = event.state
            hypothesis_id = str(row.get("hypothesis_id", ""))
            invalid_hypothesis_ids.add(hypothesis_id)
            validation.at[index, "verification_status"] = "rejected"
            validation.at[index, "rejection_reason"] = "invalid_submitted_accession"
            validation.at[index, "review_reason"] = ""
            validation.at[index, "accession_provider_error"] = ""
        elif event.state == "not_found_http_404":
            validation.at[index, "accession_resolution_status"] = event.state
            validation.at[index, "verification_status"] = "rejected"
            validation.at[index, "rejection_reason"] = "accession_not_found_in_authoritative_database"
            validation.at[index, "review_reason"] = ""
            validation.at[index, "accession_provider_error"] = ""
        elif event.state in {"rate_limited", "provider_unavailable", "response_parse_error"}:
            validation.at[index, "accession_resolution_status"] = event.state
            validation.at[index, "verification_status"] = "review_required"
            validation.at[index, "rejection_reason"] = ""
            validation.at[index, "review_reason"] = "accession_provider_unavailable"
            if not str(validation.at[index, "accession_provider_error"]):
                validation.at[index, "accession_provider_error"] = event.detail

    validation.to_csv(validation_path, sep="\t", index=False, lineterminator="\n")
    rejected = validation[validation["verification_status"] == "rejected"]
    rejected_path = _Path(str(result["rejected"]))
    rejected.to_csv(rejected_path, sep="\t", index=False, lineterminator="\n")
    review_validation = validation[validation["verification_status"] == "review_required"]
    _rewrite_review_ledger(
        _Path(str(result["review_required"])),
        invalid_hypothesis_ids=invalid_hypothesis_ids,
    )

    manifest_path = _Path(str(result["manifest"]))
    manifest = _json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict):
        raise ValueError("hypothesis-verification manifest must be a JSON object")
    manifest["verified_count"] = int((validation["verification_status"] == "verified").sum())
    if manifest.get("agent_claim_review"):
        manifest["agent_review_validation_sha256"] = _hashlib.sha256(validation_path.read_bytes()).hexdigest()
    manifest["rejected_count"] = int(len(rejected))
    manifest["review_required_count"] = int(len(review_validation))
    manifest["rejection_reasons"] = {
        str(key): int(value)
        for key, value in rejected["rejection_reason"].value_counts().sort_index().items()
    }
    manifest["uniprot_entry_outcome_counts"] = dict(
        sorted(_Counter(event.state for event in events).items())
    )
    manifest["invalid_submitted_accession_hypothesis_ids"] = sorted(invalid_hypothesis_ids)
    gates = manifest.setdefault("gates", [])
    if isinstance(gates, list):
        gates.append(
            "invalid UniProt syntax and HTTP 400 are submitted-input errors; HTTP 404 is definitive not-found; only 429/5xx/timeouts/connection failures are provider-unavailable"
        )
    manifest_path.write_text(
        _json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    result["verified_count"] = manifest["verified_count"]
    result["rejected_count"] = manifest["rejected_count"]
    result["review_required_count"] = manifest["review_required_count"]
    return result


def hypothesis_verify(*args: object, **kwargs: object) -> dict[str, object]:
    """Run the existing gate and attach typed accession/provider outcomes."""
    events_by_hypothesis: dict[tuple[str, str], list[_UniProtEntryEvent]] = {}
    result = _ORIGINAL_HYPOTHESIS_VERIFY(
        *args, **kwargs, _entry_events_by_hypothesis=events_by_hypothesis,
    )
    events = [event for row_events in events_by_hypothesis.values() for event in row_events]
    return _apply_entry_outcomes(result, events, events_by_hypothesis=events_by_hypothesis)


_legacy.hypothesis_verify = hypothesis_verify
