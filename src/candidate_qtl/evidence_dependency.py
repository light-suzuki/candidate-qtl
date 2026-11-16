"""Audit evidence dependence separately from evidence row counts.

Multiple evidence rows may derive from the same primary publication, curated
assertion, orthology annotation, or reference corpus.  This module groups rows
only when a caller supplies the same ``root_source_id``. Missing provenance is
never treated as independent evidence; it remains ``dependency_unknown``.

The output is reporting-only and never contributes to candidate ranking.
"""

from __future__ import annotations

import hashlib

import numpy as np
import pandas as pd

from .schemas import EVIDENCE_DEPENDENCY_KINDS as EVIDENCE_KINDS

DEPENDENCY_STATUSES = frozenset({"root_source_known", "dependency_unknown"})

# Compound, out-of-schema cell values: pd.isna() would return an array/Series
# for these and an ``if`` on that raises an ambiguous truth-value error, so they
# are rejected explicitly (fail-closed) instead of being coerced to text.
_COMPOUND_VALUE_TYPES = (list, tuple, dict, set, np.ndarray, pd.Series)


def _text(value: object) -> str:
    if isinstance(value, _COMPOUND_VALUE_TYPES):
        raise ValueError(
            f"evidence dependency value must be a scalar, got {type(value).__name__}"
        )
    if pd.isna(value):
        return ""
    return str(value).strip()


def _cluster_id(root_source_id: str) -> str:
    digest = hashlib.sha256(root_source_id.encode("utf-8")).hexdigest()[:20]
    return f"dep_{digest}"


def normalize_evidence_kind_column(evidence: pd.DataFrame) -> pd.Series:
    """Return ``evidence_kind`` normalized to the canonical vocabulary.

    Blank/missing cells become ``unknown``; values are casefolded and
    space-separated tokens joined with underscores (human-readable input like
    ``sequence domain call`` is accepted). Any value still outside the
    canonical vocabulary raises ``ValueError`` (fail-closed). Callers outside
    this module (the case-run evidence step) use this exact rule so written
    evidence.tsv can never carry a value the rank step's schema rejects.
    """
    kinds = evidence["evidence_kind"].map(
        lambda value: _text(value).casefold().replace(" ", "_") or "unknown"
    )
    invalid = sorted(set(kinds) - EVIDENCE_KINDS)
    if invalid:
        raise ValueError(f"unsupported evidence_kind: {', '.join(invalid)}")
    return kinds


def annotate_evidence_dependencies(evidence: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with explicit evidence-kind/root/dependency audit columns.

    Legacy evidence tables without either new field remain valid. Their rows get
    ``evidence_kind=unknown`` and ``dependency_status=dependency_unknown``;
    ``source_id`` is deliberately *not* promoted to a root source automatically.
    """
    result = evidence.copy()
    if "evidence_kind" not in result.columns:
        result["evidence_kind"] = "unknown"
    if "root_source_id" not in result.columns:
        result["root_source_id"] = ""

    kinds = normalize_evidence_kind_column(result)
    result["evidence_kind"] = kinds
    result["root_source_id"] = result["root_source_id"].map(_text)
    known = result["root_source_id"].ne("")
    result["dependency_status"] = "dependency_unknown"
    result.loc[known, "dependency_status"] = "root_source_known"
    result["dependency_cluster"] = "dependency_unknown"
    result.loc[known, "dependency_cluster"] = result.loc[known, "root_source_id"].map(_cluster_id)
    return result


def dependency_cluster_summary(evidence: pd.DataFrame) -> pd.DataFrame:
    """Summarize evidence rows by known root source plus one unknown bucket."""
    frame = annotate_evidence_dependencies(evidence)
    if frame.empty:
        return pd.DataFrame(columns=[
            "dependency_cluster", "dependency_status", "root_source_id", "evidence_count",
            "evidence_kinds", "routes", "source_ids",
        ])

    rows: list[dict[str, object]] = []
    for cluster, group in frame.groupby("dependency_cluster", sort=True, dropna=False):
        status = str(group.iloc[0]["dependency_status"])
        roots = sorted({value for value in group["root_source_id"].map(_text) if value})
        kinds = sorted({value for value in group["evidence_kind"].map(_text) if value})
        route_values: set[str] = set()
        for column in ("route", "route_id", "source_route"):
            if column in group.columns:
                route_values.update(value for value in group[column].map(_text) if value)
        source_ids = sorted({value for value in group.get("source_id", pd.Series(dtype=str)).map(_text) if value})
        rows.append({
            "dependency_cluster": cluster,
            "dependency_status": status,
            "root_source_id": roots[0] if len(roots) == 1 else "",
            "evidence_count": len(group),
            "evidence_kinds": ";".join(kinds) if kinds else "unknown",
            "routes": ";".join(sorted(route_values)) if route_values else "unknown",
            "source_ids": ";".join(source_ids) if source_ids else "unknown",
        })
    return pd.DataFrame(rows).sort_values(["dependency_status", "dependency_cluster"], kind="stable").reset_index(drop=True)


def dependency_reporting_counts(evidence: pd.DataFrame) -> dict[str, object]:
    """Return row count, known independent-cluster count, and unknown coverage."""
    frame = annotate_evidence_dependencies(evidence)
    known = frame[frame["dependency_status"].eq("root_source_known")]
    unknown = frame[frame["dependency_status"].eq("dependency_unknown")]
    return {
        "evidence_count": int(len(frame)),
        "independent_evidence_cluster_count": int(known["dependency_cluster"].nunique()),
        "dependency_unknown_evidence_count": int(len(unknown)),
        "dependency_unknown_present": bool(len(unknown)),
        "rank_contribution": "none",
    }
