"""Optional QTL-effect context contract and sign-concordance audit helpers.

The contract is deliberately orthogonal to ranking: it records what the mapped
QTL says about allele effect and lets a hypothesis declare an expected mechanism
sign.  The resulting concordance is an auditable review signal only; this module
never drops candidates or changes a score.

Both YAML and TSV inputs are supported. Missing values are normalized to the
literal ``unknown`` so old cases that do not provide an effect contract remain
valid and explicit rather than acquiring an implicit biological interpretation.
"""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

import pandas as pd
import yaml

QTL_EFFECT_COLUMNS = [
    "case_id",
    "favorable_parent",
    "effect_direction",
    "additive_effect",
    "dominance_class",
    "effect_ci",
    "environment_id",
    "phenotyping_protocol",
    "stable_across_environments",
    "source_id",
    "source_detail",
]

SIGN_CONCORDANCE = frozenset({"concordant", "incompatible", "unknown"})
_SIGN_POSITIVE = frozenset({"+", "+1", "positive", "increase", "increased", "gain", "up", "higher"})
_SIGN_NEGATIVE = frozenset({"-", "-1", "negative", "decrease", "decreased", "loss", "down", "lower"})
_UNKNOWN = frozenset({"", "unknown", "na", "n/a", "none", "null", "not_evaluated", "not evaluated"})


def _is_missing_scalar(value: object) -> bool:
    """True only for scalar missing markers (None, pd.NA, pd.NaT, nan).

    List/array-like values are never missing: a YAML list in e.g. ``effect_ci``
    is a legitimate text value and must not be swallowed by ``pd.isna``.
    """
    if value is None:
        return True
    if not pd.api.types.is_scalar(value):
        return False
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def _text(value: object) -> str:
    if _is_missing_scalar(value):
        return "unknown"
    text = str(value).strip()
    return text if text else "unknown"


def _stable_across_environments(value: object) -> str:
    text = _text(value).casefold()
    if text in {"true", "yes", "y", "1", "stable"}:
        return "true"
    if text in {"false", "no", "n", "0", "unstable"}:
        return "false"
    if text in _UNKNOWN:
        return "unknown"
    raise ValueError("stable_across_environments must be true, false, or unknown")


def normalize_qtl_effect_contract(raw: Mapping[str, object] | None, *, case_id: str = "") -> dict[str, str]:
    """Return one canonical, backward-compatible QTL-effect contract."""
    source = dict(raw or {})
    normalized = {column: _text(source.get(column)) for column in QTL_EFFECT_COLUMNS}
    if case_id:
        incoming_case = normalized["case_id"]
        if incoming_case not in {"unknown", case_id}:
            raise ValueError(f"qtl_effect_contract case_id {incoming_case!r} does not match {case_id!r}")
        normalized["case_id"] = case_id
    normalized["stable_across_environments"] = _stable_across_environments(source.get("stable_across_environments"))
    return normalized


def load_qtl_effect_contract(path_arg: str | Path | None, *, case_id: str = "") -> dict[str, str]:
    """Load a contract from YAML or TSV; absence yields an all-unknown contract."""
    if path_arg is None:
        return normalize_qtl_effect_contract({}, case_id=case_id)
    path = Path(path_arg)
    suffix = path.suffix.casefold()
    if suffix in {".yaml", ".yml"}:
        payload = yaml.safe_load(path.read_text(encoding="utf-8"))
        if payload is None:
            payload = {}
        if not isinstance(payload, dict):
            raise ValueError("QTL-effect YAML must contain an object")
        top_case = str(payload.get("case_id", "")).strip()
        nested = payload.get("qtl_effect_contract", payload)
        if not isinstance(nested, dict):
            raise ValueError("qtl_effect_contract must be an object")
        # A declared top-level case_id that disagrees with the caller must fail
        # closed before normalizing: otherwise case_b could silently inherit
        # case_a's QTL-effect claim when the nested block omits its own id.
        if case_id and top_case and top_case != case_id:
            raise ValueError(f"qtl_effect_contract top-level case_id {top_case!r} does not match {case_id!r}")
        requested_case = case_id or top_case
        return normalize_qtl_effect_contract(nested, case_id=requested_case)
    if suffix == ".tsv":
        frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        if frame.empty:
            # A header-only (or otherwise empty) TSV resolves no row: with a
            # requested case that is a fail-closed mismatch, not an unknown.
            if case_id:
                raise ValueError(f"QTL-effect TSV has no row for case_id {case_id!r}")
            return normalize_qtl_effect_contract({}, case_id=case_id)
        if case_id:
            if "case_id" not in frame.columns:
                if len(frame) != 1:
                    raise ValueError("multi-row QTL-effect TSV requires a case_id column")
            else:
                frame = frame[frame["case_id"].str.strip().eq(case_id)]
                if frame.empty:
                    raise ValueError(f"QTL-effect TSV has no row for case_id {case_id!r}")
        if len(frame) != 1:
            raise ValueError("QTL-effect TSV must resolve to exactly one row per case")
        return normalize_qtl_effect_contract(frame.iloc[0].to_dict(), case_id=case_id)
    raise ValueError("QTL-effect contract must be .yaml/.yml or .tsv")


def contract_frame(contract: Mapping[str, object]) -> pd.DataFrame:
    """Materialize the canonical contract as a one-row audit table."""
    normalized = normalize_qtl_effect_contract(contract, case_id=str(contract.get("case_id", "")).strip())
    return pd.DataFrame([normalized], columns=QTL_EFFECT_COLUMNS)


def _mechanism_sign(value: object) -> str:
    text = _text(value).casefold()
    if text in _SIGN_POSITIVE:
        return "positive"
    if text in _SIGN_NEGATIVE:
        return "negative"
    return "unknown"


def assess_sign_concordance(expected_mechanism_sign: object, effect_direction: object) -> tuple[str, str]:
    """Classify expected mechanism vs observed QTL direction without scoring it."""
    expected = _mechanism_sign(expected_mechanism_sign)
    observed = _mechanism_sign(effect_direction)
    if expected == "unknown" or observed == "unknown":
        return "unknown", f"expected_mechanism_sign={expected}; qtl_effect_sign={observed}"
    if expected == observed:
        return "concordant", f"expected_mechanism_sign={expected}; qtl_effect_sign={observed}"
    return "incompatible", f"expected_mechanism_sign={expected}; qtl_effect_sign={observed}"


def annotate_sign_concordance(hypotheses: pd.DataFrame, contract: Mapping[str, object]) -> pd.DataFrame:
    """Append audit-only sign fields to a hypothesis ledger copy.

    ``expected_mechanism_sign`` may already be present on each hypothesis. If it
    is absent or blank the decision is ``unknown``. Existing score/rank columns
    are copied byte-for-value and are never interpreted here.
    """
    result = hypotheses.copy()
    normalized = normalize_qtl_effect_contract(contract, case_id=str(contract.get("case_id", "")).strip())
    if "expected_mechanism_sign" not in result.columns:
        result["expected_mechanism_sign"] = "unknown"
    result["expected_mechanism_sign"] = result["expected_mechanism_sign"].map(_text)
    decisions = result["expected_mechanism_sign"].map(
        lambda expected: assess_sign_concordance(expected, normalized["effect_direction"])
    )
    result["sign_concordance"] = decisions.map(lambda item: item[0])
    result["sign_concordance_reason"] = decisions.map(lambda item: item[1])
    result["qtl_effect_source_id"] = normalized["source_id"]
    result["qtl_effect_source_detail"] = normalized["source_detail"]
    result["qtl_effect_environment_id"] = normalized["environment_id"]
    return result
