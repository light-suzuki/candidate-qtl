"""Locus-level causal scenarios kept separate from gene ranking.

The ledger represents plausible causal forms for a mapped locus without forcing
one winner. Missing evidence remains ``unassessed``; a provider gap that made
evidence impossible to produce is ``not_evaluated`` (never converted into
contradicting evidence, never rounded into unassessed). All outputs are
reporting/review artifacts with ``rank_contribution=none``.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np
import pandas as pd
import yaml

# Issue #363 review categories: `novel_or_unannotated_mechanism` reviews an
# unknown/unannotated mechanism (not a known multi-gene hypothesis),
# `phenotype_context_review` reviews phenotype-measurement validity, and
# `mapping_adequacy_review` reviews mapping adequacy without asserting a
# mapping/liftover error. They are kept separate from the pre-existing types
# (multi_gene_or_epistasis / environment_specific / mapping_or_liftover_error)
# because folding them in would discard the review intent and the
# named-mechanism distinction. `provider_missing` marks a row whose discordance
# cannot be classified because a route's provider was absent; it never names a
# causal form, so it is not folded into any review category either.
SCENARIO_TYPES = frozenset({
    "coding_gene",
    "regulatory_element",
    "sv_te_cnv_pav",
    "reference_parent_mismatch",
    "mapping_or_liftover_error",
    "multi_gene_or_epistasis",
    "environment_specific",
    "novel_or_unannotated_mechanism",
    "phenotype_context_review",
    "mapping_adequacy_review",
    "provider_missing",
})
# `not_evaluated` marks a row whose evidence could not be produced at all
# (provider missing / semantically not evaluable), distinct from `unassessed`
# where evaluation was simply never attempted. It is a separate state so the
# ledger never silently rounds a provider gap into unassessed.
ASSESSMENT_STATES = frozenset({
    "unassessed", "supported", "contradicted", "mixed", "review_required", "not_evaluated",
})
SCENARIO_COLUMNS = [
    "case_id",
    "scenario_id",
    "scenario_type",
    "supporting_evidence",
    "contradicting_evidence",
    "assessment_state",
    "recommended_next_experiment",
    "reason",
    "source_ids",
    "rank_contribution",
]


def _is_missing(value: object) -> bool:
    if value is None:
        return True
    if isinstance(value, (str, bytes)):
        return False
    try:
        result = pd.isna(value)
        if isinstance(result, np.ndarray):
            return bool(result) if result.ndim == 0 else False
        return bool(result)
    except (TypeError, ValueError):
        return False


def _text(value: object, default: str = "unknown") -> str:
    if _is_missing(value):
        return default
    text = str(value).strip()
    return text if text else default


def _evidence_text(value: object) -> str:
    return _text(value, "")


def _scenario_id(row: Mapping[str, object], case_id: str) -> str:
    parts = [
        case_id,
        _text(row.get("scenario_type")),
        _evidence_text(row.get("supporting_evidence")),
        _evidence_text(row.get("contradicting_evidence")),
        _text(row.get("recommended_next_experiment")),
        _text(row.get("reason")),
        _text(row.get("source_ids")),
    ]
    digest = hashlib.sha256("\x1f".join(parts).encode("utf-8")).hexdigest()[:16]
    return f"SCN-{digest}"


def normalize_scenario(row: Mapping[str, object], *, case_id: str = "") -> dict[str, str]:
    """Normalize and validate one scenario without inferring negative evidence."""
    incoming_case = _text(row.get("case_id"))
    resolved_case = case_id or (incoming_case if incoming_case != "unknown" else "")
    if not resolved_case:
        raise ValueError("causal scenario requires case_id")
    if case_id and incoming_case not in {"unknown", case_id}:
        raise ValueError(f"causal scenario case_id {incoming_case!r} does not match {case_id!r}")

    scenario_type = _text(row.get("scenario_type")).casefold()
    if scenario_type not in SCENARIO_TYPES:
        raise ValueError(f"unsupported scenario_type: {scenario_type!r}")
    support = _evidence_text(row.get("supporting_evidence"))
    contradiction = _evidence_text(row.get("contradicting_evidence"))
    requested_state = _text(row.get("assessment_state"), "unassessed").casefold()
    if requested_state not in ASSESSMENT_STATES:
        raise ValueError(f"unsupported assessment_state: {requested_state!r}")
    if not support and not contradiction and requested_state not in {"unassessed", "not_evaluated"}:
        raise ValueError("scenario without support or contradiction must remain unassessed or not_evaluated")

    normalized = {
        "case_id": resolved_case,
        "scenario_id": _text(row.get("scenario_id"), ""),
        "scenario_type": scenario_type,
        "supporting_evidence": support,
        "contradicting_evidence": contradiction,
        "assessment_state": requested_state,
        "recommended_next_experiment": _text(row.get("recommended_next_experiment")),
        "reason": _text(row.get("reason")),
        "source_ids": _text(row.get("source_ids")),
        "rank_contribution": "none",
    }
    if not normalized["scenario_id"]:
        normalized["scenario_id"] = _scenario_id(normalized, resolved_case)
    return normalized


def build_scenario_ledger(rows: Iterable[Mapping[str, object]], *, case_id: str = "") -> pd.DataFrame:
    """Build a deterministic multi-scenario ledger for one or more cases."""
    if isinstance(rows, pd.DataFrame):
        rows = rows.to_dict(orient="records")
    normalized = [normalize_scenario(row, case_id=case_id) for row in rows]
    if not normalized:
        return pd.DataFrame(columns=SCENARIO_COLUMNS)
    frame = pd.DataFrame(normalized, columns=SCENARIO_COLUMNS)
    duplicate_ids = frame[frame["scenario_id"].duplicated(keep=False)]
    if not duplicate_ids.empty:
        conflicting = duplicate_ids.drop_duplicates()
        if len(conflicting) != duplicate_ids["scenario_id"].nunique():
            raise ValueError("scenario_id collision across non-identical rows")
        frame = frame.drop_duplicates(subset=["scenario_id"], keep="first")
    return frame.sort_values(["case_id", "scenario_type", "scenario_id"], kind="stable").reset_index(drop=True)


def write_scenario_ledger(rows: Iterable[Mapping[str, object]], out_path: str | Path, *, case_id: str = "") -> Path:
    """Write canonical ``locus_causal_scenarios.tsv`` output."""
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame = build_scenario_ledger(rows, case_id=case_id)
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")
    return path


def scenario_summary(ledger: pd.DataFrame) -> pd.DataFrame:
    """Report support/contradiction/unassessed/not_evaluated state by scenario type."""
    required = {"scenario_type", "assessment_state", "supporting_evidence", "contradicting_evidence"}
    missing = required - set(ledger.columns)
    if missing:
        raise ValueError(f"scenario ledger missing columns: {', '.join(sorted(missing))}")
    if ledger.empty:
        return pd.DataFrame(columns=[
            "scenario_type", "scenario_count", "supported_rows", "contradicted_rows",
            "unassessed_rows", "not_evaluated_rows", "support_present_rows", "contradiction_present_rows",
        ])
    rows: list[dict[str, object]] = []
    for scenario_type, group in ledger.groupby("scenario_type", sort=True):
        states = group["assessment_state"].map(_text)
        support_evidence = group["supporting_evidence"].map(_evidence_text)
        contradiction_evidence = group["contradicting_evidence"].map(_evidence_text)
        rows.append({
            "scenario_type": scenario_type,
            "scenario_count": len(group),
            "supported_rows": int(states.eq("supported").sum()),
            "contradicted_rows": int(states.eq("contradicted").sum()),
            "unassessed_rows": int(states.eq("unassessed").sum()),
            "not_evaluated_rows": int(states.eq("not_evaluated").sum()),
            "support_present_rows": int(support_evidence.ne("").sum()),
            "contradiction_present_rows": int(contradiction_evidence.ne("").sum()),
        })
    return pd.DataFrame(rows)


def render_scenario_report(ledger: pd.DataFrame) -> str:
    """Render a compact Markdown review section without choosing a winner."""
    if ledger.empty:
        return "## Locus causal scenarios\n\nNo causal scenarios have been registered; status is unassessed.\n"
    lines = ["## Locus causal scenarios", ""]
    for row in ledger.sort_values(["case_id", "scenario_type", "scenario_id"], kind="stable").to_dict(orient="records"):
        support = _evidence_text(row.get("supporting_evidence", "")) or "none recorded"
        contradiction = _evidence_text(row.get("contradicting_evidence", "")) or "none recorded"
        next_experiment = _text(row.get("recommended_next_experiment"))
        reason = _text(row.get("reason"))
        lines.append(
            f"- {_text(row.get('case_id'))}/{_text(row.get('scenario_id'), '')} "
            f"`{_text(row.get('scenario_type'))}` / `{_text(row.get('assessment_state'))}` — "
            f"support: {support}; contradiction: {contradiction}; next: {next_experiment}; "
            f"reason: {reason}"
        )
    return "\n".join(lines) + "\n"


def load_scenario_rows(path: str | Path) -> list[dict[str, object]]:
    """Load scenario rows from TSV or YAML (list of mappings) for case-run wiring.

    TSV is read with string dtype and ``keep_default_na=False`` so empty cells
    stay empty strings instead of NaN; YAML must be a list of mappings (or an
    empty document). Parse problems raise ``ValueError``/``OSError`` so the
    caller can degrade without crashing the case.
    """
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(f"scenario input does not exist: {path}")
    if source.suffix.casefold() in {".yaml", ".yml"}:
        try:
            payload = yaml.safe_load(source.read_text(encoding="utf-8"))
        except yaml.YAMLError as error:
            raise ValueError(f"scenario rows YAML could not be parsed: {error}") from error
        if payload is None:
            return []
        if not isinstance(payload, list):
            raise ValueError("scenario YAML input must be a list of scenario mappings")
        rows: list[dict[str, object]] = []
        for entry in payload:
            if not isinstance(entry, Mapping):
                raise ValueError("scenario YAML entries must be mappings")
            rows.append(dict(entry))
        return rows
    frame = pd.read_csv(source, sep="\t", dtype=str, keep_default_na=False)
    return frame.to_dict(orient="records")


def load_qtl_effect_contract(path: str | Path) -> list[dict[str, object]]:
    """Load the optional #360 qtl_effect_contract resource defensively.

    YAML accepts either a top-level list of contract mappings or a
    ``qtl_effect_contract:`` key holding one; TSV is read with string dtype and
    no NaN rewriting. Parse problems raise ``ValueError``/``OSError`` so the
    caller records ``invalid`` instead of crashing.
    """
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(f"qtl_effect_contract does not exist: {path}")
    if source.suffix.casefold() in {".yaml", ".yml"}:
        try:
            payload = yaml.safe_load(source.read_text(encoding="utf-8"))
        except yaml.YAMLError as error:
            raise ValueError(f"qtl_effect_contract YAML could not be parsed: {error}") from error
        if payload is None:
            return []
        if isinstance(payload, Mapping):
            payload = payload.get("qtl_effect_contract")
        if isinstance(payload, Mapping):
            payload = [payload]
        if not isinstance(payload, list):
            raise ValueError("qtl_effect_contract must be a list of contract mappings (or YAML with a qtl_effect_contract list)")
        rows: list[dict[str, object]] = []
        for entry in payload:
            if not isinstance(entry, Mapping):
                raise ValueError("qtl_effect_contract entries must be mappings")
            rows.append(dict(entry))
        return rows
    frame = pd.read_csv(source, sep="\t", dtype=str, keep_default_na=False)
    return frame.to_dict(orient="records")


def sign_concordance_evidence(contract_rows: Iterable[Mapping[str, object]], case_id: str) -> tuple[str | None, str | None]:
    """Return ``(sign_concordance, evidence_text)`` for the case_id's contract row.

    Only ``incompatible`` and ``concordant`` rows bind evidence:
    ``sign discordance: favorable_parent/effect_direction`` /
    ``sign concordance: ...`` (missing parts are dropped). A missing row, an
    ``unknown`` value, or a row without a usable value returns ``(None, None)``
    so the caller records ``unknown`` and appends nothing. The first row
    matching ``case_id`` wins.
    """
    if isinstance(contract_rows, pd.DataFrame):
        contract_rows = contract_rows.to_dict(orient="records")
    target = str(case_id).strip()
    for row in contract_rows:
        if not isinstance(row, Mapping):
            continue
        if str(row.get("case_id", "")).strip() != target:
            continue
        concordance = str(row.get("sign_concordance", "")).strip().casefold()
        if concordance not in {"incompatible", "concordant"}:
            return None, None
        parts = "/".join(
            str(row.get(key, "")).strip()
            for key in ("favorable_parent", "effect_direction")
            if str(row.get(key, "")).strip()
        )
        label = "sign discordance" if concordance == "incompatible" else "sign concordance"
        text = f"{label}: {parts}" if parts else label
        return concordance, text
    return None, None
