"""Adapter for the frozen legacy 7-column mechanism payloads (issue #380).

mechanism-import accepts one strict 29-column hypothesis schema; the frozen T11
benchmark payloads (case_id, mechanism_id, mechanism_name, evidence_class,
source_log_entry, status, selection_reason) predate it. This module detects the
payload schema from the raw header line (never pd.read_csv, which mangles
duplicate columns into x, x.1), converts the legacy shape into canonical
29-column rows with every absent field left empty -- never inferred -- and
records the out-of-band legacy fields in a sidecar adapter map instead of
smuggling them into notes or confidence.

Fail-closed by design: an unknown schema, duplicate header columns, a forbidden
leakage column, a legacy row with an empty case_id/mechanism_id or a duplicated
mechanism_id, a missing payload, or a freeze-hash mismatch raises ValueError
with the fixed message formats below; a schema mismatch never yields a silent
0-row success. Schema acceptance is not bait verification: offline verification
stays review_required and legacy rows without an accession are rejected
downstream exactly as before.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from .claim_evidence import EXTRACTION_COLUMNS
from .relevance_tiers import INITIAL_FIELDS
from .phenotype_investigation import (
    INVESTIGATION_PROMPT_VERSION,
    MECHANISM_COLUMNS,
    PARSER_VERSION,
)

LEGACY_V1_7COL_COLUMNS = ("case_id", "mechanism_id", "mechanism_name", "evidence_class", "source_log_entry", "status", "selection_reason")
TEMPLATE_24_COLUMNS = frozenset(column for column in MECHANISM_COLUMNS
                                if column not in {"hypothesis_id", "parser_version", "review_status", "rank_contribution", "notes",
                                                  "source_title", "source_identifier_aliases", *INITIAL_FIELDS})
ADAPTER_VERSION = "mechanism-payload-adapter-1.0.0"
OUTPUT_SCHEMA_VERSION = "mechanism_hypotheses_29col_v2"
FORBIDDEN_TOKENS = ("truth", "gold", "answer", "candidate", "rank", "score", "interval", "accession", "bait")
_LEAKAGE_CLASS_LABEL = "truth|candidate|rank|score|interval|accession|bait"
_ADAPTER_MAP_COLUMNS = ("hypothesis_id", "source_row", "mechanism_id", "mechanism_name", "evidence_class", "source_log_entry", "status")
_CANONICAL_COLUMNS = set(MECHANISM_COLUMNS)
_LEGACY_COLUMNS = set(LEGACY_V1_7COL_COLUMNS)
_JSON_WRAPPER_KEYS = {"hypotheses", "mechanisms", "proposals", "results", "claim_extractions", "response", "parsed", "output", "text", "content", "output_text"}
_EXTRACTION_REQUIRED_KEYS = frozenset(column for column in EXTRACTION_COLUMNS
                                      if column in {"citation", "evidence_relation", "quoted_support"})
_MISSING = object()


class _PaperPayloadError(ValueError):
    """Keep the original record location available to read-only preflight."""

    def __init__(self, location: str, error: Exception) -> None:
        self.location = location
        super().__init__(f"{location}: {error}")


def preflight_payload_file(path: str | Path, freeze_manifest_arg: str | Path | None,
                           case_id: str) -> dict[str, object]:
    """Check the original import transport without creating or freezing artifacts.

    This shares the importer's schema, raw-byte hash and paper normalizer. It
    does not assert investigation completion, source verification or support.
    """
    from .core import runtime_path
    from .paper_reports import imported_paper_fields
    from .phenotype_investigation_legacy import _payload_records

    payload = runtime_path(str(path))
    try:
        validation = validate_payload_file(payload, freeze_manifest_arg, case_id)
        # JSON validation already visits every nested paper report. TSV cells
        # reach the very same normalizer through the actual import parser.
        records = _payload_records(payload, str(validation["schema"]))
        for index, record in enumerate(records):
            try:
                imported_paper_fields(record)
            except (ValueError, TypeError) as error:
                raise _PaperPayloadError(f"rows[{index}].paper_reports_json", error) from error
    except (ValueError, TypeError, OSError) as error:
        return {
            "command": "mechanism-import", "preflight": True, "status": "invalid",
            "payload": str(payload), "exit_code": 1,
            "errors": [{"path": getattr(error, "location", "$"), "error": str(error)}],
        }
    return {
        "command": "mechanism-import", "preflight": True, "status": "valid",
        "payload": str(payload), "validation_scope": "payload_schema_and_paper_transport",
        **validation,
        "advisories": _reference_id_advisories(records),
    }


def _reference_id_advisories(records: list) -> list[dict[str, object]]:
    """Name-only reference targets that later resolution may not carry to projection.

    Advisory only: the payload stays valid and nothing is filled in. A gene name
    that the provider cannot resolve, or that matches several records, keeps the
    hypothesis out of interval projection (contract: no copy is chosen from a
    name alone). The original researcher may add the stable ID before freeze.
    """
    advisories = []
    for index, record in enumerate(records):
        if not isinstance(record, dict):
            continue
        name = str(record.get("reference_gene_name", "") or "").strip()
        stable = str(record.get("proposed_stable_id", "") or "").strip()
        if name and not stable:
            advisories.append({
                "path": f"rows[{index}]", "advisory": "reference_gene_without_stable_id",
                "reference_gene_name": name,
                "reference_species": str(record.get("reference_species", "") or "").strip(),
                "effect": "name-only targets may stay unresolved or ambiguous and then cannot be projected onto the interval",
            })
    return advisories


def normalize_source_ids(value: object, *, location: str = "source_ids") -> str:
    """Return the canonical semicolon form without changing list order/content.

    Route A JSON payloads may use either the canonical string form or a
    structured list[str].  Do not stringify arbitrary JSON values: a nested
    list/object or a scalar member would become a Python representation that
    the downstream citation parser cannot interpret losslessly.
    """
    if isinstance(value, str):
        return value
    if isinstance(value, list):
        invalid = next((item for item in value if not isinstance(item, str)), _MISSING)
        if invalid is not _MISSING:
            raise ValueError(
                f"{location} must be a string or list[str]; "
                f"list member has type {type(invalid).__name__}"
            )
        return ";".join(value)
    raise ValueError(f"{location} must be a string or list[str]; got {type(value).__name__}")


def detect_payload_schema(header: list[str]) -> str:
    """Classify a payload header: legacy_v1_7col | canonical_29_col | canonical_24_col_template | canonical_subset | unknown.

    Legacy and 29-column matches are order-insensitive; duplicate or empty
    header columns always report unknown.
    """
    names = [name.strip() for name in header]
    if not names or any(name == "" for name in names) or len(names) != len(set(names)):
        return "unknown"
    name_set = set(names)
    if name_set == _LEGACY_COLUMNS:
        return "legacy_v1_7col"
    if name_set == _CANONICAL_COLUMNS:
        return "canonical_29_col"
    if name_set == TEMPLATE_24_COLUMNS:
        return "canonical_24_col_template"
    if name_set <= _CANONICAL_COLUMNS:
        return "canonical_subset"
    return "unknown"


def forbidden_columns(header: list[str]) -> list[str]:
    """Columns outside the canonical schema whose names carry leakage tokens.

    Canonical column names are exempt (rank_contribution, interval_informed);
    the denylist check applies only to columns never seen in MECHANISM_COLUMNS.
    """
    return [name for name in header
            if name.strip() not in _CANONICAL_COLUMNS and any(token in name.casefold() for token in FORBIDDEN_TOKENS)]


def convert_legacy_v1_7col(
    records: list[dict[str, object]],
    *,
    case_id: str = "",
    payload_source: str = "llm",
    prompt_version: str = INVESTIGATION_PROMPT_VERSION,
    proposal_round: int = 1,
    interval_informed: bool = False,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Map legacy 7-column records to canonical rows plus the adapter-map sidecar.

    Only mechanism_name -> mechanism_statement, source_log_entry -> source_ids
    and selection_reason -> selection_rationale are copied verbatim; case_id
    comes from the caller with the payload's own case_id as fallback. Every
    other canonical field is empty by design: the legacy payload carries no
    gene family, accession, counter-hypothesis or model evidence, and inventing
    them is forbidden. mechanism_id, evidence_class, source_log_entry and
    status are out-of-band -- adapter map only, never into the canonical row.
    The caller loop assigns hypothesis_id (mech_{index:03d}); map rows carry the
    same deterministic ids because the conversion is strictly 1:1.
    """
    canonical_rows: list[dict[str, object]] = []
    adapter_map_rows: list[dict[str, object]] = []
    for index, record in enumerate(records, 1):
        def cell(name: str) -> str:
            value = record.get(name)
            return "" if value is None else str(value)
        payload_case_id = cell("case_id").strip()
        if case_id and payload_case_id and payload_case_id != case_id:
            raise ValueError(
                f"mechanism payload case mismatch: requested case_id '{case_id}' does not match "
                f"payload case_id '{payload_case_id}' at row {index}"
            )
        row = {column: "" for column in MECHANISM_COLUMNS}
        row.update({
            "case_id": payload_case_id or case_id,
            "mechanism_statement": cell("mechanism_name"),
            "source_ids": cell("source_log_entry"),
            "selection_rationale": cell("selection_reason"),
            "confidence": "unknown",
            "payload_source": payload_source,
            "model": "",
            "prompt_version": prompt_version,
            "parser_version": PARSER_VERSION,
            "retrieval_provenance": "",
            "proposal_round": proposal_round,
            "interval_informed": bool(interval_informed),
            "review_status": "pending_verification",
            "rank_contribution": "none",
        })
        canonical_rows.append(row)
        adapter_map_rows.append({
            "hypothesis_id": f"mech_{index:03d}",
            "source_row": index,
            "mechanism_id": cell("mechanism_id"),
            "mechanism_name": cell("mechanism_name"),
            "evidence_class": cell("evidence_class"),
            "source_log_entry": cell("source_log_entry"),
            "status": cell("status"),
        })
    return canonical_rows, adapter_map_rows


def write_adapter_map(rows: list[dict[str, object]], out: Path) -> Path:
    """Write the out-of-band legacy fields as mechanism_payload_adapter_map.tsv."""
    path = Path(out) / "mechanism_payload_adapter_map.tsv"
    lines = ["\t".join(_ADAPTER_MAP_COLUMNS)]
    lines.extend("\t".join(str(row.get(name, "")) for name in _ADAPTER_MAP_COLUMNS) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return path


def validate_payload_file(path: str | Path, freeze_manifest_arg: str | Path | None, case_id: str) -> dict[str, object]:
    """Fail-closed schema and freeze validation of a mechanism payload.

    Returns {"schema", "rows_read", "raw_payload_sha256", "freeze_verification",
    "errors"} and raises ValueError on every rejection with the message formats
    fixed by the issue #380 design: unknown/duplicate/forbidden columns, legacy
    rows with empty case_id/mechanism_id or duplicate mechanism_id, and freeze
    hash mismatches. The sha256 is computed over the raw payload bytes exactly
    as frozen in coordination/payload_freeze.tsv; row numbers in legacy errors
    are 1-based on data rows.

    TSV payloads are classified by their header line (never pd.read_csv, which
    mangles duplicate columns). JSON/JSONL payloads are accepted -- mechanism-
    import keeps its own tolerant parse, so here only provenance and the freeze
    check are recorded -- which also lets the guided-profile preflight record a
    JSON payload as accepted instead of falsely rejecting it as an unknown TSV
    schema.
    """
    payload = Path(path)
    if not payload.is_file():
        raise ValueError(f"mechanism payload schema rejected: {payload}; file not found or unreadable")
    try:
        raw = payload.read_bytes()
    except OSError:
        raise ValueError(f"mechanism payload schema rejected: {payload}; file not found or unreadable")
    if not raw:
        raise ValueError(f"mechanism payload schema rejected: {payload}; file is empty")
    lines = raw.decode("utf-8", errors="replace").splitlines()
    if not lines or not lines[0].strip():
        raise ValueError(f"mechanism payload schema rejected: {payload}; file is empty")
    digest = hashlib.sha256(raw).hexdigest()
    if lines[0].lstrip().startswith(("{", "[")):
        return _json_payload_result(lines, raw, digest, payload, freeze_manifest_arg, case_id)
    header = [cell.strip() for cell in lines[0].rstrip("\r").split("\t")]
    duplicates = {name for name in header if header.count(name) > 1}
    if duplicates:
        entries: list[str] = []
        for name in duplicates:
            entries.extend([name] * header.count(name))
        raise ValueError(f"mechanism payload schema rejected (schema=unknown): {payload}; duplicate column(s): [{', '.join(entries)}]")
    schema = detect_payload_schema(header)
    data_lines = [line for line in lines[1:] if line.strip()]
    if schema == "unknown":
        parts = [
            f"unknown column(s): [{', '.join(name for name in header if name not in _CANONICAL_COLUMNS)}]",
            f"missing canonical column(s): [{', '.join(name for name in MECHANISM_COLUMNS if name not in set(header))}]",
        ]
        for name in forbidden_columns(header):
            parts.append(f"forbidden column '{name}' (leakage class: {_LEAKAGE_CLASS_LABEL}); forbidden columns are never accepted")
        raise ValueError(f"mechanism payload schema rejected (schema=unknown): {payload}; " + "; ".join(parts))
    payload_case_id = _tsv_payload_case_id(header, data_lines, payload, schema, case_id)
    if schema == "legacy_v1_7col":
        _validate_legacy_rows(header, data_lines, payload)
    return {
        "schema": schema,
        "rows_read": len(data_lines),
        "payload_case_id": payload_case_id,
        "raw_payload_sha256": digest,
        "freeze_verification": _verify_freeze(payload, freeze_manifest_arg, digest, case_id or payload_case_id),
        "errors": [],
    }


def _json_payload_result(lines: list[str], raw: bytes, digest: str, payload: Path,
                         freeze_manifest_arg: str | Path | None, case_id: str) -> dict[str, object]:
    """Provenance record for a JSON/JSONL payload: schema, row count and freeze."""
    schema, records, roots = _json_schema_and_records(lines, raw, payload)
    _validate_json_tree(roots, payload)
    payload_case_id = _records_case_id(records, payload, case_id)
    return {
        "schema": schema,
        "rows_read": len(records),
        "payload_case_id": payload_case_id,
        "raw_payload_sha256": digest,
        "freeze_verification": _verify_freeze(payload, freeze_manifest_arg, digest, case_id or payload_case_id),
        "errors": [],
    }


def _json_schema_and_records(lines: list[str], raw: bytes, payload: Path) -> tuple[str, list[object], list[object]]:
    try:
        parsed = json.loads(raw.decode("utf-8", errors="replace"))
    except ValueError:
        data_lines = [line for line in lines if line.strip()]
        if not data_lines or not all(_parses_json(line) for line in data_lines):
            raise ValueError(f"mechanism payload schema rejected (schema=unknown): {payload}; "
                             "payload is not a supported schema (canonical/legacy TSV, JSON or JSONL)")
        parsed_lines = [json.loads(line) for line in data_lines]
        return "jsonl", _extract_json_records(parsed_lines), parsed_lines
    return "json", _extract_json_records([parsed]), [parsed]


def _extract_json_records(payloads: list[object]) -> list[object]:
    records: list[object] = []
    for item in payloads:
        if isinstance(item, list):
            records.extend(_extract_json_records(item))
        elif isinstance(item, dict):
            for key in ("hypotheses", "mechanisms", "proposals", "results"):
                value = item.get(key)
                if isinstance(value, list):
                    records.extend(value)
                    break
            else:
                nested_found = False
                for key in ("response", "parsed", "output", "text", "content", "output_text"):
                    value = item.get(key)
                    if isinstance(value, str):
                        try:
                            value = json.loads(value)
                        except ValueError:
                            continue
                    if isinstance(value, (dict, list)):
                        records.extend(_extract_json_records([value]))
                        nested_found = True
                if not nested_found:
                    records.append(item)
    return records


def _validate_json_tree(values: list[object], payload: Path, location: str = "$") -> None:
    for index, value in enumerate(values):
        item_location = f"{location}[{index}]" if len(values) > 1 else location
        if isinstance(value, list):
            for position, member in enumerate(value):
                _validate_json_tree([member], payload, f"{item_location}[{position}]")
            continue
        if not isinstance(value, dict):
            raise ValueError(
                f"mechanism payload schema rejected (schema=unknown): {payload}; "
                f"JSON value at {item_location} is not an object"
            )
        keys = [str(key) for key in value]
        forbidden = forbidden_columns(keys)
        if forbidden:
            name = forbidden[0]
            raise ValueError(
                f"mechanism payload schema rejected (schema=unknown): {payload}; forbidden JSON key '{name}' "
                f"(leakage class: {_LEAKAGE_CLASS_LABEL}); forbidden keys are never accepted"
            )
        is_wrapper = any(name in _JSON_WRAPPER_KEYS for name in keys)
        allowed = _JSON_WRAPPER_KEYS if is_wrapper else _CANONICAL_COLUMNS
        unknown = [name for name in keys if name not in allowed]
        if unknown:
            raise ValueError(
                f"mechanism payload schema rejected (schema=unknown): {payload}; "
                f"unknown JSON key(s) at {item_location}: [{', '.join(unknown)}]"
            )
        if "paper_reports_json" in value:
            from .paper_reports import imported_paper_fields
            try:
                imported_paper_fields(value)
            except (ValueError, TypeError) as error:
                raise _PaperPayloadError(f"{item_location}.paper_reports_json", error) from error
        if "source_ids" in value:
            try:
                normalize_source_ids(value["source_ids"], location=f"{item_location}.source_ids")
            except ValueError as error:
                raise ValueError(
                    f"mechanism payload schema rejected (schema=json): {payload}; {error}"
                ) from error
        if is_wrapper:
            extractions = value.get("claim_extractions")
            if isinstance(extractions, list):
                # Agent-extracted claim spans ride in the same payload but use
                # the extraction schema, not the mechanism one.
                for position, row in enumerate(extractions):
                    row_location = f"{item_location}.claim_extractions[{position}]"
                    if not isinstance(row, dict):
                        raise ValueError(
                            f"mechanism payload schema rejected (schema=json): {payload}; "
                            f"JSON value at {row_location} is not an object"
                        )
                    missing = sorted(column for column in _EXTRACTION_REQUIRED_KEYS
                                     if not str(row.get(column, "")).strip())
                    if missing:
                        raise ValueError(
                            f"mechanism payload schema rejected (schema=json): {payload}; "
                            f"claim extraction at {row_location} is missing field(s): [{', '.join(missing)}]"
                        )
            for key, nested in value.items():
                if key == "claim_extractions":
                    continue  # validated above against the extraction schema
                if isinstance(nested, str):
                    try:
                        nested = json.loads(nested)
                    except ValueError:
                        continue
                if isinstance(nested, (dict, list)):
                    _validate_json_tree([nested], payload, f"{item_location}.{key}")


def _records_case_id(records: list[object], payload: Path, requested_case_id: str) -> str:
    case_ids = {str(record.get("case_id", "")).strip() for record in records if isinstance(record, dict) and str(record.get("case_id", "")).strip()}
    return _resolve_payload_case_id(case_ids, payload, requested_case_id)


def _tsv_payload_case_id(header: list[str], data_lines: list[str], payload: Path, schema: str,
                         requested_case_id: str) -> str:
    if "case_id" not in header:
        return ""
    case_col = header.index("case_id")
    case_ids = {cells[case_col].strip() for line in data_lines
                if len(cells := line.split("\t")) > case_col and cells[case_col].strip()}
    return _resolve_payload_case_id(case_ids, payload, requested_case_id, schema)


def _resolve_payload_case_id(case_ids: set[str], payload: Path, requested_case_id: str,
                             schema: str = "json") -> str:
    if len(case_ids) > 1:
        raise ValueError(
            f"mechanism payload schema rejected (schema={schema}): {payload}; multiple case_id values: "
            f"[{', '.join(sorted(case_ids))}]"
        )
    payload_case_id = next(iter(case_ids), "")
    if requested_case_id and payload_case_id and requested_case_id != payload_case_id:
        raise ValueError(
            f"mechanism payload case mismatch: requested case_id '{requested_case_id}' does not match "
            f"payload case_id '{payload_case_id}'"
        )
    return payload_case_id


def _parses_json(line: str) -> bool:
    try:
        json.loads(line)
        return True
    except ValueError:
        return False


def _validate_legacy_rows(header: list[str], data_lines: list[str], payload: Path) -> None:
    first_seen: dict[str, int] = {}
    for index, line in enumerate(data_lines, 1):
        cells = line.split("\t")
        row = {name: (cells[i] if i < len(cells) else "") for i, name in enumerate(header)}
        if not row["case_id"].strip():
            raise ValueError(f"mechanism payload schema rejected (schema=legacy_v1_7col): {payload}; row {index}: case_id is empty")
        mechanism_id = row["mechanism_id"].strip()
        if not mechanism_id:
            raise ValueError(f"mechanism payload schema rejected (schema=legacy_v1_7col): {payload}; row {index}: mechanism_id is empty")
        if mechanism_id in first_seen:
            raise ValueError(f"mechanism payload schema rejected (schema=legacy_v1_7col): {payload}; row {index}: "
                             f"duplicate mechanism_id '{mechanism_id}' (first at row {first_seen[mechanism_id]})")
        first_seen[mechanism_id] = index


def _verify_freeze(payload: Path, freeze_manifest_arg: str | Path | None, digest: str, case_id: str) -> dict[str, object] | None:
    if not freeze_manifest_arg:
        return None
    freeze = Path(freeze_manifest_arg)
    try:
        lines = freeze.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        raise ValueError(f"mechanism payload freeze mismatch: {freeze}; file not found or unreadable")
    if not lines:
        raise ValueError(f"mechanism payload freeze mismatch: {freeze}; file is empty")
    header = [cell.strip() for cell in lines[0].rstrip("\r").split("\t")]
    if "case_id" not in header or "mechanism_payload_sha256" not in header:
        raise ValueError(f"mechanism payload freeze mismatch: {freeze}; missing case_id/mechanism_payload_sha256 column(s)")
    case_col, sha_col = header.index("case_id"), header.index("mechanism_payload_sha256")
    expected: str | None = None
    for line in lines[1:]:
        cells = line.split("\t")
        if len(cells) > max(case_col, sha_col) and cells[case_col].strip() == case_id:
            expected = cells[sha_col].strip()
            break
    if expected is None:
        raise ValueError(f"mechanism payload freeze mismatch: {freeze} has no row for case_id '{case_id}'")
    if expected != digest:
        raise ValueError(f"mechanism payload freeze mismatch: {payload}; expected sha256 {expected} from {freeze}; actual {digest}")
    return {"source": str(freeze), "expected_sha256": expected, "actual_sha256": digest, "matched": True}
