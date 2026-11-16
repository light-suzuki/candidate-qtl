"""Wet-candidate budget gate for #426 / #434.

The Dry workflow is successful when it safely compresses the interval to a
Wet-reviewable candidate set; it does not need to force a single causal gene.
This module therefore implements the stop/escalation decision *without ranking
or truncating candidates*:

    supported members
      -> count
         == 0          : abstain / review_required
         <= wet_budget : freeze and stop
         >  wet_budget : freeze all members and report the target excess

An explicitly requested ``strict`` resource budget retains the historical
bounded refinement/stop behavior. The default target is not biological evidence.

No LLM is called here.  The gate only decides whether an expensive/adaptive
stage is necessary, records why, and preserves the entire candidate set.  It
therefore records ``adaptive_refinement_required`` rather than pretending that
an adaptive/LLM stage has already been invoked.
"""

from __future__ import annotations

import hashlib
import json

import pandas as pd

from .core import runtime_path

SCHEMA_VERSION = "1.2"
DEFAULT_WET_BUDGET = 20
DEFAULT_MAX_ADAPTIVE_ROUNDS = 1

STATUSES = (
    "no_supported_candidates",
    "within_wet_budget",
    "above_wet_target",
    "adaptive_refinement_required",
    "review_required_over_budget",
)
ACTIONS = (
    "stop_review_required",
    "freeze_candidate_set",
    "request_adaptive_refinement",
)


def _candidate_id_column(frame: pd.DataFrame) -> str:
    for column in ("candidate_gene_id", "gene_id"):
        if column in frame.columns:
            return column
    raise ValueError("candidate budget input requires candidate_gene_id or gene_id")


def normalize_candidate_set(frame: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    """Validate a candidate set without changing membership or review order."""
    column = _candidate_id_column(frame)
    normalized = frame.copy().fillna("")
    normalized[column] = normalized[column].astype(str).str.strip()
    if normalized[column].eq("").any():
        raise ValueError(f"{column} values must be non-empty")
    if normalized[column].duplicated().any():
        duplicates = sorted(set(normalized.loc[normalized[column].duplicated(keep=False), column]))
        raise ValueError(f"candidate set contains duplicate ids: {', '.join(duplicates)}")
    return normalized, column


def budget_decision(
    candidate_count: int,
    *,
    wet_budget: int = DEFAULT_WET_BUDGET,
    adaptive_round: int = 0,
    max_adaptive_rounds: int = DEFAULT_MAX_ADAPTIVE_ROUNDS,
    budget_policy: str = "target",
) -> dict[str, object]:
    """Return the deterministic stop/escalation decision for one candidate set."""
    if candidate_count < 0:
        raise ValueError("candidate_count must be non-negative")
    if budget_policy not in {"target", "strict"}:
        raise ValueError("budget_policy must be target or strict")
    if wet_budget < 1:
        raise ValueError("wet_budget must be at least 1")
    if adaptive_round < 0:
        raise ValueError("adaptive_round must be non-negative")
    if max_adaptive_rounds < 0:
        raise ValueError("max_adaptive_rounds must be non-negative")
    if adaptive_round > max_adaptive_rounds:
        raise ValueError("adaptive_round cannot exceed max_adaptive_rounds")

    if candidate_count == 0:
        return {
            "status": "no_supported_candidates",
            "action": "stop_review_required",
            "adaptive_refinement_required": False,
            "decision_reason": "candidate_set_empty",
            "review_required": True,
        }
    if candidate_count <= wet_budget:
        return {
            "status": "within_wet_budget",
            "action": "freeze_candidate_set",
            "adaptive_refinement_required": False,
            "decision_reason": "within_wet_budget",
            "review_required": False,
        }
    if budget_policy == "target":
        return {
            "status": "above_wet_target",
            "action": "freeze_candidate_set",
            "adaptive_refinement_required": False,
            "decision_reason": "retain_supported_members_above_operational_target",
            "review_required": False,
        }
    if adaptive_round < max_adaptive_rounds:
        return {
            "status": "adaptive_refinement_required",
            "action": "request_adaptive_refinement",
            "adaptive_refinement_required": True,
            "decision_reason": "candidate_count_exceeds_wet_budget",
            "review_required": False,
        }
    return {
        "status": "review_required_over_budget",
        "action": "stop_review_required",
        "adaptive_refinement_required": False,
        "decision_reason": "adaptive_round_budget_exhausted",
        "review_required": True,
    }


def candidate_budget_gate(
    candidates_arg: str,
    out_dir_arg: str,
    *,
    wet_budget: int = DEFAULT_WET_BUDGET,
    adaptive_round: int = 0,
    max_adaptive_rounds: int = DEFAULT_MAX_ADAPTIVE_ROUNDS,
    case_id: str = "",
    budget_policy: str = "target",
) -> dict[str, object]:
    """Freeze the input candidate set and record whether refinement is needed.

    The output candidate TSV is a membership-preserving copy.  This function
    never takes the first N rows, never assigns a numerical score, and never
    removes an unresolved paralog merely to satisfy the budget.
    """
    source = runtime_path(candidates_arg)
    if not source.exists():
        raise ValueError(f"candidate budget input does not exist: {source}")
    frame = pd.read_csv(source, sep="\t", dtype=str).fillna("")
    normalized, id_column = normalize_candidate_set(frame)

    decision = budget_decision(
        int(len(normalized)),
        wet_budget=wet_budget,
        adaptive_round=adaptive_round,
        max_adaptive_rounds=max_adaptive_rounds,
        budget_policy=budget_policy,
    )

    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    frozen_path = out / "wet_candidate_set.tsv"
    normalized.to_csv(frozen_path, sep="\t", index=False, lineterminator="\n")
    source_sha256 = hashlib.sha256(source.read_bytes()).hexdigest()
    output_sha256 = hashlib.sha256(frozen_path.read_bytes()).hexdigest()

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id,
        "candidate_id_column": id_column,
        "candidate_count": int(len(normalized)),
        "wet_budget": wet_budget,
        "budget_policy": budget_policy,
        "above_target_count": max(0, int(len(normalized)) - wet_budget),
        "adaptive_round": adaptive_round,
        "max_adaptive_rounds": max_adaptive_rounds,
        **decision,
        "input_sha256": source_sha256,
        "output_sha256": output_sha256,
        "output": str(frozen_path),
        "rank_contribution": "none",
        "policy": [
            "an empty supported set is an explicit abstention/review-required state, not a successful within-budget result",
            "candidate_count <= wet_budget stops immediately; no adaptive/heavy refinement is requested",
            "the default target is operational guidance, never a scientific exclusion or count-only stop",
            "only an explicitly requested strict resource budget requires bounded refinement or an over-budget stop",
            "the gate preserves every input candidate and never ranks, scores, truncates, or removes a supported copy",
            "this component records whether refinement is required; the orchestrator must separately record whether an adaptive/LLM stage was actually invoked",
        ],
    }
    manifest_path = out / "candidate_budget_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        "candidate_set": str(frozen_path),
        "manifest": str(manifest_path),
        "candidate_count": int(len(normalized)),
        **decision,
    }
