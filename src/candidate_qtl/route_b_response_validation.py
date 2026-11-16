"""Local model-claim shape checks; artifact identity is validated separately."""
from __future__ import annotations

import json
from collections.abc import Iterable


def route_b_claim_schema_error(claim: object) -> str:
    """Return a quarantine reason without repairing a model's claim meaning."""
    if not isinstance(claim, dict):
        return "semantic claim must be a JSON object"
    claim_type = claim.get("claim_type")
    if not isinstance(claim_type, str) or not claim_type.strip():
        return "semantic claim claim_type must be a non-empty string"
    from .paper_reports import imported_paper_fields
    try:
        imported_paper_fields(claim)
    except (ValueError, TypeError) as error:
        return str(error)
    return ""


def route_b_response_completed(record: dict) -> bool:
    """Derive execution coverage from raw claim shape, never claimed metadata."""
    if record.get("parse_status") != "json":
        return False
    if not isinstance(record.get("model_json"), dict):
        raise ValueError("Route B completed response requires model_json")
    raw_text = record.get("raw_text")
    if not isinstance(raw_text, str):
        raise ValueError("Route B execution response raw text must be a string")
    try:
        raw_json = json.loads(raw_text)
    except json.JSONDecodeError as exc:
        raise ValueError("Route B execution response raw text is not valid JSON") from exc
    if record["model_json"] != raw_json:
        raise ValueError("Route B execution response model JSON differs from raw text")
    claims = raw_json.get("claims")
    return isinstance(claims, list) and all(not route_b_claim_schema_error(claim) for claim in claims)


def route_b_response_state(record: dict) -> str:
    """Describe transport execution, without interpreting a researcher's reason."""
    if route_b_response_completed(record):
        return "completed_claims" if record["model_json"]["claims"] else "completed_no_claim"
    if record.get("parse_status") == "json":
        claims = record["model_json"].get("claims")
        if isinstance(claims, list) and any(not route_b_claim_schema_error(claim) for claim in claims):
            return "partial_claims"
        return "invalid_response"
    original = record.get("agent_response_raw")
    if isinstance(original, dict):
        if original.get("status") == "not_evaluated":
            return "declared_not_evaluated"
        if original.get("status") == "failed":
            return "declared_failed"
    if record.get("response_format_error") or record.get("rejected_claims"):
        return "invalid_response"
    if (original is None and not record.get("attempt_source_file")
            and not str(record.get("raw_text", "")).strip()
            and record.get("not_evaluated_reason") == "no_partial_response"):
        return "missing_response"
    return "unclassified_not_evaluated"


def validate_route_b_not_run(records: object) -> None:
    """A deliberate omission cannot carry successful, partial, or hidden claims."""
    if not isinstance(records, list):
        raise ValueError("Route B not-run response requires a candidates list")
    if any(not isinstance(row, dict) or route_b_response_state(row) not in
           {"declared_not_evaluated", "missing_response"}
           or row.get("model_json") not in (None, {})
           or row.get("claims") not in (None, [])
           or (isinstance(row.get("agent_response_raw"), dict)
               and row["agent_response_raw"].get("claims") not in (None, [])) for row in records):
        raise ValueError("Route B not-run responses must all be not evaluated and contain no claims")


def route_b_execution_breakdown(records: Iterable[dict], *, manifest: dict | None = None) -> dict[str, int]:
    """Recompute an additive execution view; old manifests need no new fields."""
    counts = dict.fromkeys(("completed_claims", "completed_no_claim", "declared_not_evaluated",
                            "declared_failed", "invalid_response", "partial_claims", "missing_response",
                            "unclassified_not_evaluated"), 0)
    for record in records:
        counts[route_b_response_state(record)] += 1
    if manifest is not None and "response_state_counts" in manifest:
        declared = manifest["response_state_counts"]
        if (not isinstance(declared, dict) or set(declared) != set(counts)
                or any(type(value) is not int or value < 0 for value in declared.values())
                or declared != counts):
            raise ValueError("Route B response state counts do not match response rows")
    return counts
