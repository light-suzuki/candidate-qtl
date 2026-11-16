"""Reporting-only benchmark coverage, uncertainty, and paired comparison helpers.

The central rule is that conditional performance and evaluation coverage are
separate quantities. ``not_evaluated`` is never converted to a biological miss,
and contaminated cases are never included in performance confidence intervals.

Issue #367 wiring rules shared by every consumer:

* a missing fingerprint (``None`` / NaN / ``pd.NA`` / ``pd.NaT`` / blank
  string) is *missing*, so two runs that both carry an explicit null
  fingerprint are not comparable (``<key>:missing``), never equal.
* case ids are canonicalized before any roster comparison or merge
  (``_canonical_case_id``): ``1`` and ``"1"`` are the same case across runs,
  but the same run carrying both raw forms is a fail-closed collision.
* compound cells (list / dict / ndarray / Series) are fail-closed everywhere
  (``ValueError``) instead of being coerced through ``str()``.
"""

from __future__ import annotations

import hashlib
import json
import math
import random
from collections import Counter
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

EVALUATION_STATUSES = frozenset({"clean", "evaluated", "not_evaluated", "contaminated"})
_PERFORMANCE_ELIGIBLE = frozenset({"clean", "evaluated"})
_COMPARABILITY_KEYS = ("roster_fingerprint", "truth_fingerprint", "metric_contract_fingerprint")

_COMPOUND_TYPES = (list, dict, tuple, set, frozenset)


def _is_compound(value: object) -> bool:
    if isinstance(value, _COMPOUND_TYPES):
        return True
    if isinstance(value, (np.ndarray, pd.Series)):
        return True
    return False


def _scalar_text(value: object) -> str:
    """Fail-closed ``str()`` for scalar cells; compound values are rejected."""
    if _is_compound(value):
        raise ValueError(f"expected a scalar value, got {type(value).__name__}")
    return str(value)


def _is_missing(value: object) -> bool:
    """True for None, NaN (float or numpy), pd.NA, pd.NaT, and blank strings."""
    if value is None:
        return True
    if isinstance(value, str):
        return not value.strip()
    if _is_compound(value):
        return False
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def _canonical_case_id(value: object) -> str:
    """Canonical scalar case-id string; compound and missing values are fail-closed."""
    if _is_compound(value):
        raise ValueError(f"case id must be a scalar value, got {type(value).__name__}")
    if _is_missing(value):
        raise ValueError("case id is missing")
    return str(value).strip()


def _fingerprint_value(value: object) -> str | None:
    """Canonical fingerprint string; None for missing; compound values fail-closed.

    The reserved ``"missing"`` sentinel (produced by the benchmark-compare
    fallback when no fingerprint can be derived) is treated as missing so two
    runs that both lack a fingerprint are never comparable.
    """
    if _is_compound(value):
        raise ValueError(f"fingerprint value must be a scalar, got {type(value).__name__}")
    if _is_missing(value):
        return None
    if value == "missing":
        return None
    return str(value).strip()


def roster_fingerprint(case_ids: Sequence[object]) -> str:
    """SHA-256 over the canonical sorted case-id list, one ``"<id>\\n"`` per line."""
    canonical = sorted({_canonical_case_id(value) for value in case_ids})
    return hashlib.sha256("".join(f"{value}\n" for value in canonical).encode("utf-8")).hexdigest()


def _status(value: object) -> str:
    if _is_compound(value):
        raise ValueError(f"evaluation status must be a scalar value, got {type(value).__name__}")
    text = str(value).strip().casefold()
    if text not in EVALUATION_STATUSES:
        raise ValueError(f"unsupported evaluation status: {text or '<empty>'}")
    return text


def _hit(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if _is_compound(value):
        raise ValueError(f"binary hit value must be a scalar, got {type(value).__name__}")
    text = str(value).strip().casefold()
    if text in {"1", "true", "yes", "hit", "supported"}:
        return True
    if text in {"0", "false", "no", "miss", "unsupported"}:
        return False
    raise ValueError(f"binary hit value is not recognized: {value!r}")


def wilson_interval(successes: int, total: int, *, z: float = 1.959963984540054) -> tuple[float | None, float | None]:
    """Wilson score 95% interval by default; undefined for a zero denominator."""
    if total < 0 or successes < 0 or successes > total:
        raise ValueError("require 0 <= successes <= total")
    if total == 0:
        return None, None
    p = successes / total
    z2 = z * z
    denominator = 1.0 + z2 / total
    center = (p + z2 / (2.0 * total)) / denominator
    margin = z * math.sqrt((p * (1.0 - p) + z2 / (4.0 * total)) / total) / denominator
    return max(0.0, center - margin), min(1.0, center + margin)


def coverage_summary(statuses: Sequence[object]) -> dict[str, object]:
    """Count coverage states without folding them into a performance score."""
    normalized = [_status(value) for value in statuses]
    counts = Counter(normalized)
    total = len(normalized)
    evaluated = counts["clean"] + counts["evaluated"]
    clean = counts["clean"]
    return {
        "total_cases": total,
        "evaluated_cases": evaluated,
        "clean_cases": clean,
        "not_evaluated_cases": counts["not_evaluated"],
        "contaminated_cases": counts["contaminated"],
        "evaluation_coverage": (evaluated / total) if total else None,
        "clean_coverage": (clean / total) if total else None,
        "rank_contribution": "none",
    }


def conditional_binary_performance(
    frame: pd.DataFrame,
    *,
    hit_column: str = "hit",
    status_column: str = "evaluation_status",
) -> dict[str, object]:
    """Compute a hit rate only on evaluable, non-contaminated cases."""
    if hit_column not in frame.columns or status_column not in frame.columns:
        raise ValueError(f"frame requires {hit_column!r} and {status_column!r}")
    statuses = frame[status_column].map(_status)
    eligible = frame[statuses.isin(_PERFORMANCE_ELIGIBLE)]
    total = len(eligible)
    if total == 0:
        return {
            "metric_status": "not_evaluated",
            "hits": 0,
            "evaluated_cases": 0,
            "conditional_hit_rate": None,
            "wilson_95_low": None,
            "wilson_95_high": None,
            **coverage_summary(statuses.tolist()),
        }
    hits = sum(_hit(value) for value in eligible[hit_column])
    low, high = wilson_interval(hits, total)
    return {
        "metric_status": "evaluated",
        "hits": hits,
        "evaluated_cases": total,
        "conditional_hit_rate": hits / total,
        "wilson_95_low": low,
        "wilson_95_high": high,
        **coverage_summary(statuses.tolist()),
    }


def evaluation_coverage_block(*, total_cases: int, evaluated_cases: int, hit_count: int) -> dict[str, object]:
    """``conditional_binary_performance``-shaped block for an already-computed
    denominator (clean = evaluated, contaminated = 0). Used by consumers whose
    per-case status vocabulary has no contamination concept."""
    low, high = wilson_interval(hit_count, evaluated_cases) if evaluated_cases else (None, None)
    return {
        "metric_status": "evaluated" if evaluated_cases else "not_evaluated",
        "hits": hit_count,
        "evaluated_cases": evaluated_cases,
        "conditional_hit_rate": (hit_count / evaluated_cases) if evaluated_cases else None,
        "wilson_95_low": low,
        "wilson_95_high": high,
        "total_cases": total_cases,
        "clean_cases": evaluated_cases,
        "not_evaluated_cases": total_cases - evaluated_cases,
        "contaminated_cases": 0,
        "evaluation_coverage": (evaluated_cases / total_cases) if total_cases else None,
        "clean_coverage": (evaluated_cases / total_cases) if total_cases else None,
        "rank_contribution": "none",
    }


def comparability_status(
    left: Mapping[str, object],
    right: Mapping[str, object],
    *,
    required_keys: Sequence[str] = _COMPARABILITY_KEYS,
) -> dict[str, object]:
    """Fail closed unless roster, truth, and metric contracts are all identical."""
    mismatches: list[str] = []
    for key in required_keys:
        a = _fingerprint_value(left.get(key, ""))
        b = _fingerprint_value(right.get(key, ""))
        if a is None or b is None:
            mismatches.append(f"{key}:missing")
        elif a != b:
            mismatches.append(f"{key}:mismatch")
    return {
        "status": "comparable" if not mismatches else "not_comparable",
        "reasons": ";".join(mismatches) if mismatches else "none",
    }


def _case_outcome(status: str, hit: object) -> str:
    if status == "contaminated":
        return "contaminated"
    if status == "not_evaluated":
        return "not_evaluated"
    return "hit" if _hit(hit) else "miss"


def _percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    if not ordered:
        raise ValueError("cannot take percentile of an empty sequence")
    index = min(len(ordered) - 1, max(0, round(probability * (len(ordered) - 1))))
    return ordered[index]


def continuous_bootstrap_ci(
    values: Sequence[float],
    *,
    iterations: int = 2000,
    seed: int,
) -> dict[str, float | None]:
    """Case bootstrap 95% CI for a continuous quantity (MRR etc.).

    ``random.Random(seed)`` makes the resampling deterministic for a fixed
    seed; an empty ``values`` sequence yields ``None`` point estimate and CI
    (undefined, never NaN).
    """
    if not values:
        return {"point_estimate": None, "ci_95_low": None, "ci_95_high": None}
    if iterations <= 0:
        raise ValueError("iterations must be positive")
    sample = [float(value) for value in values]
    point_estimate = sum(sample) / len(sample)
    rng = random.Random(seed)
    boot: list[float] = []
    for _ in range(iterations):
        resampled = [sample[rng.randrange(len(sample))] for _ in range(len(sample))]
        boot.append(sum(resampled) / len(resampled))
    return {
        "point_estimate": point_estimate,
        "ci_95_low": _percentile(boot, 0.025),
        "ci_95_high": _percentile(boot, 0.975),
    }


def paired_binary_delta_ci(
    pairs: Sequence[tuple[bool, bool]],
    *,
    iterations: int = 2000,
    seed: int = 0,
) -> tuple[float | None, float | None, float | None]:
    """Paired case-bootstrap CI for hit-rate delta (right minus left)."""
    if not pairs:
        return None, None, None
    if iterations <= 0:
        raise ValueError("iterations must be positive")
    delta = sum((1 if right else 0) - (1 if left else 0) for left, right in pairs) / len(pairs)
    rng = random.Random(seed)
    boot: list[float] = []
    for _ in range(iterations):
        sampled = [pairs[rng.randrange(len(pairs))] for _ in range(len(pairs))]
        boot.append(sum((1 if right else 0) - (1 if left else 0) for left, right in sampled) / len(sampled))
    return delta, _percentile(boot, 0.025), _percentile(boot, 0.975)


def paired_binary_comparison(
    left: pd.DataFrame,
    right: pd.DataFrame,
    left_manifest: Mapping[str, object],
    right_manifest: Mapping[str, object],
    *,
    case_id_column: str = "case_id",
    hit_column: str = "hit",
    status_column: str = "evaluation_status",
    bootstrap_iterations: int = 2000,
    bootstrap_seed: int = 0,
) -> dict[str, object]:
    """Build case-level transitions and paired hit-rate delta for comparable runs.

    Case ids are canonicalized (``1`` == ``"1"``) before the roster check and
    the merge; both use the canonical key, and the merge is
    ``validate="one_to_one"``. The same run carrying two raw forms of one
    canonical id (``1`` and ``"1"``) is a fail-closed ``ValueError``.
    """
    comparable = comparability_status(left_manifest, right_manifest)
    if comparable["status"] != "comparable":
        return {**comparable, "transition_counts": {}, "paired_case_count": 0, "hit_rate_delta": None,
                "delta_ci_95_low": None, "delta_ci_95_high": None, "case_transitions": []}
    required = {case_id_column, hit_column, status_column}
    normalized_frames: dict[str, pd.DataFrame] = {}
    for name, frame in (("left", left), ("right", right)):
        missing = required - set(frame.columns)
        if missing:
            raise ValueError(f"{name} run missing columns: {', '.join(sorted(missing))}")
        normalized = frame[[case_id_column, hit_column, status_column]].copy()
        keys = normalized[case_id_column].map(_canonical_case_id)
        if keys.duplicated().any():
            raise ValueError(
                f"{name} run contains duplicate case_id (case id canonicalization collision)"
            )
        normalized["_case_key"] = keys
        normalized_frames[name] = normalized
    left_normalized = normalized_frames["left"]
    right_normalized = normalized_frames["right"]
    left_ids = set(left_normalized["_case_key"])
    right_ids = set(right_normalized["_case_key"])
    if left_ids != right_ids:
        return {"status": "not_comparable", "reasons": "case_roster:mismatch", "transition_counts": {},
                "paired_case_count": 0, "hit_rate_delta": None, "delta_ci_95_low": None,
                "delta_ci_95_high": None, "case_transitions": []}

    merged = left_normalized[["_case_key", hit_column, status_column]].merge(
        right_normalized[["_case_key", hit_column, status_column]], on="_case_key", suffixes=("_left", "_right"),
        validate="one_to_one"
    )
    transitions: Counter[str] = Counter()
    pairs: list[tuple[bool, bool]] = []
    case_transitions: list[dict[str, object]] = []
    for row in merged.sort_values("_case_key").to_dict(orient="records"):
        left_status = _status(row[f"{status_column}_left"])
        right_status = _status(row[f"{status_column}_right"])
        left_outcome = _case_outcome(left_status, row[f"{hit_column}_left"])
        right_outcome = _case_outcome(right_status, row[f"{hit_column}_right"])
        transition = f"{left_outcome}->{right_outcome}"
        transitions[transition] += 1
        case_transitions.append({
            "case_id": row["_case_key"],
            "left_evaluation_status": left_status,
            "left_hit": _hit(row[f"{hit_column}_left"]) if left_status in _PERFORMANCE_ELIGIBLE else None,
            "right_evaluation_status": right_status,
            "right_hit": _hit(row[f"{hit_column}_right"]) if right_status in _PERFORMANCE_ELIGIBLE else None,
            "transition": transition,
        })
        if left_status in _PERFORMANCE_ELIGIBLE and right_status in _PERFORMANCE_ELIGIBLE:
            pairs.append((_hit(row[f"{hit_column}_left"]), _hit(row[f"{hit_column}_right"])))
    delta, low, high = paired_binary_delta_ci(pairs, iterations=bootstrap_iterations, seed=bootstrap_seed)
    return {
        "status": "comparable",
        "reasons": "none",
        "transition_counts": dict(sorted(transitions.items())),
        "paired_case_count": len(pairs),
        "hit_rate_delta": delta,
        "delta_ci_95_low": low,
        "delta_ci_95_high": high,
        "case_transitions": case_transitions,
        "rank_contribution": "none",
    }


def no_nan_for_json(value: object) -> object:
    """Recursively convert float NaN (python or numpy floating) to None."""
    if isinstance(value, (float, np.floating)):
        try:
            return None if math.isnan(value) else value
        except TypeError:
            return value
    if isinstance(value, dict):
        return {key: no_nan_for_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [no_nan_for_json(item) for item in value]
    return value


def contract_fingerprint(fields: Mapping[str, object]) -> str:
    """SHA-256 over ``json.dumps(fields, sort_keys=True)`` (metric contract)."""
    return hashlib.sha256(
        json.dumps(fields, sort_keys=True).encode("utf-8")
    ).hexdigest()


def format_use_case_strata_summary_table(metrics: Mapping[str, object]) -> pd.DataFrame:
    """Format the use-case strata metrics into a tabular DataFrame for display/reporting (issue #444)."""
    rows: list[dict[str, object]] = []

    # Primary contract summary row
    prim = metrics.get("primary_retention_contract")
    if isinstance(prim, dict):
        rows.append({
            "use_case_stratum": "Type 1-3 Primary Retention Contract",
            "total_cases": prim.get("total_cases", 0),
            "evaluable_cases": prim.get("evaluable_cases", 0),
            "not_evaluated_cases": prim.get("not_evaluated_cases", 0),
            "contaminated_cases": prim.get("contaminated_cases", 0),
            "top5_retention": prim.get("top5_retention_rate"),
            "top10_retention": prim.get("top10_retention_rate"),
            "top20_retention": prim.get("top20_retention_rate"),
            "median_reduction": prim.get("median_reduction_ratio"),
            "wrong_confident_promotion": prim.get("wrong_confident_promotion_rate"),
            "guarantee_scope": "Primary Guarantee (Target >= 0.90)",
        })

    strata = metrics.get("use_case_strata")
    if isinstance(strata, dict):
        for uct, udata in strata.items():
            if isinstance(udata, dict):
                rows.append({
                    "use_case_stratum": uct,
                    "total_cases": udata.get("case_count", 0),
                    "evaluable_cases": udata.get("evaluable_count", 0),
                    "not_evaluated_cases": udata.get("not_evaluated_count", 0),
                    "contaminated_cases": udata.get("contaminated_count", 0),
                    "top5_retention": udata.get("top5_retention"),
                    "top10_retention": udata.get("top10_retention"),
                    "top20_retention": udata.get("top20_retention"),
                    "median_reduction": udata.get("median_reduction_ratio"),
                    "wrong_confident_promotion": udata.get("wrong_confident_promotion_rate"),
                    "guarantee_scope": "Retention Scope" if udata.get("retention_scope_eligible") else "Out of Scope / Honest Abstention",
                })

    return pd.DataFrame(rows)


def render_use_case_strata_report(metrics: Mapping[str, object]) -> str:
    """Render a Markdown section for the 4 QTL use-case strata and Top-20 retention contract (issue #444)."""
    prim = metrics.get("primary_retention_contract")
    if not isinstance(prim, dict):
        return "No use-case strata metrics available.\n"

    lines = [
        "### QTL Use-Case Strata & Top-20 Retention Contract (Issue #444)",
        "",
        f"- **Primary Scope (Type 1–3)**: {prim.get('evaluable_cases', 0)} evaluable cases (of {prim.get('total_cases', 0)} total)",
        f"- **Top-20 Retention Rate**: {prim.get('top20_retention_rate') if prim.get('top20_retention_rate') is not None else 'N/A'} (target >= {prim.get('target_top20_retention_min', 0.90)})",
        f"- **Top-10 Retention Rate**: {prim.get('top10_retention_rate') if prim.get('top10_retention_rate') is not None else 'N/A'}",
        f"- **Top-5 Retention Rate**: {prim.get('top5_retention_rate') if prim.get('top5_retention_rate') is not None else 'N/A'}",
        f"- **Median Candidate Reduction**: {prim.get('median_reduction_ratio') if prim.get('median_reduction_ratio') is not None else 'N/A'} (target >= {prim.get('target_median_reduction_min', 0.80)})",
        f"- **Wrong-Confident Promotion**: {prim.get('wrong_confident_promotion_rate') if prim.get('wrong_confident_promotion_rate') is not None else 'N/A'} (target <= {prim.get('target_wrong_confident_max', 0.05)})",
        "- **Type 4 Status**: Strictly excluded from Top-20 recall denominator; evaluated for honest abstention/escalation.",
        "",
    ]
    return "\n".join(lines)
