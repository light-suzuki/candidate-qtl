from __future__ import annotations

import hashlib
import json

import pandas as pd

from .core import runtime_path


def _sha256(path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


SCREENING_ORDER_HELPER_COLUMNS = (
    "_screening_rank", "_dry_score", "_support_family_count", "_independent_source_count", "_rescue_priority",
)


def prepare_screening_order(frame: pd.DataFrame, priority_column: str = "screening_rank") -> pd.DataFrame:
    """Deterministic screening-priority order shared by the shortlist cut and evidence-tier reporting.

    Extracted so ``screening_shortlist`` (the frozen assay-capacity cut, issue
    #144) and ``evidence_tier_report`` (the target_size-independent biological
    tiers, issue #229) apply *one* selection rule instead of two that could
    silently drift apart. The rule itself only orders candidates by economical
    review priority; it never filters on ``dry_group`` or ``top_eligibility``,
    which answer a different question (is the evidence sufficient), not this one
    (which candidate is cheapest to look at next).
    """
    if "gene_id" not in frame.columns:
        raise ValueError("ranked candidate table must contain gene_id")
    if frame["gene_id"].eq("").any() or frame["gene_id"].duplicated().any():
        raise ValueError("gene_id values must be unique and non-empty")

    work = frame.copy()
    if priority_column in work.columns:
        work["_screening_rank"] = pd.to_numeric(work[priority_column], errors="raise").astype(int)
    elif "dry_score" in work.columns:
        scores = pd.to_numeric(work["dry_score"], errors="raise")
        work["_screening_rank"] = scores.rank(method="min", ascending=False).astype(int)
    else:
        raise ValueError("ranked candidate table must contain screening_rank or dry_score")

    if "dry_score" in work.columns:
        work["_dry_score"] = pd.to_numeric(work["dry_score"], errors="raise")
    else:
        work["_dry_score"] = -work["_screening_rank"].astype(float)
    for source, target in (
        ("support_family_count", "_support_family_count"),
        ("independent_source_count", "_independent_source_count"),
        ("rescue_priority", "_rescue_priority"),
    ):
        work[target] = pd.to_numeric(work[source], errors="coerce").fillna(0) if source in work.columns else 0

    return work.sort_values(
        ["_screening_rank", "_dry_score", "_support_family_count", "_independent_source_count", "gene_id"],
        ascending=[True, False, False, False, True],
        kind="stable",
    )


def screening_shortlist(
    ranked_arg: str,
    out_dir_arg: str,
    target_size: int = 20,
    rescue_size: int = 5,
    *,
    priority_column: str = "screening_rank",
    stage: str = "A_deterministic_retrieval",
) -> dict[str, object]:
    """Freeze a deterministic Stage-A shortlist before human or LLM review.

    Selection deliberately uses screening rank, not ``dry_group`` or
    ``top_eligibility``.  The latter answer whether evidence is sufficient to
    promote a biological claim; this function only asks which candidates are
    economical to review next.  Boundary ties are capped deterministically and
    exported in full so a gene-ID tie break cannot be mistaken for biology.
    """
    if target_size < 1:
        raise ValueError("target_size must be at least 1")
    if rescue_size < 0:
        raise ValueError("rescue_size must be zero or greater")

    ranked_path = runtime_path(ranked_arg)
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(ranked_path, sep="\t", dtype=str).fillna("")

    ordered = prepare_screening_order(frame, priority_column=priority_column)
    primary = ordered.head(min(target_size, len(ordered))).copy()
    primary_ids = set(primary["gene_id"])
    primary["shortlist_stage"] = "deterministic_primary"
    primary["selection_reason"] = f"{priority_column}_with_audited_boundary_tie"
    primary["deterministic_review_order"] = range(1, len(primary) + 1)

    if primary.empty:
        boundary = ordered.head(0).copy()
        boundary_rank = None
    else:
        boundary_rank = int(primary.iloc[-1]["_screening_rank"])
        boundary = ordered[ordered["_screening_rank"].eq(boundary_rank)].copy()
    boundary["selected_in_primary"] = boundary["gene_id"].isin(primary_ids)

    rescue_mask = ~ordered["gene_id"].isin(primary_ids)
    if "rescue_reason" in ordered.columns:
        rescue_mask &= ordered["rescue_reason"].ne("")
    elif "dry_group" in ordered.columns:
        rescue_mask &= ordered["dry_group"].str.casefold().eq("rescue")
    else:
        rescue_mask &= False
    rescue = ordered[rescue_mask].sort_values(
        ["_rescue_priority", "_screening_rank", "gene_id"],
        ascending=[False, True, True],
        kind="stable",
    ).head(rescue_size).copy()
    rescue["shortlist_stage"] = "deterministic_rescue_appendix"
    rescue["selection_reason"] = "rescue_reason_retained_outside_primary_20"
    rescue["deterministic_review_order"] = range(1, len(rescue) + 1)

    helper_columns = [column for column in primary.columns if column.startswith("_")]
    primary = primary.drop(columns=helper_columns)
    rescue = rescue.drop(columns=[column for column in rescue.columns if column.startswith("_")])
    boundary = boundary.drop(columns=[column for column in boundary.columns if column.startswith("_")])
    excluded = frame[~frame["gene_id"].isin(primary_ids | set(rescue["gene_id"]))].copy()
    excluded["exclusion_reason"] = "outside_deterministic_primary_and_rescue_appendix"

    primary_path = out_dir / "deterministic_shortlist.tsv"
    rescue_path = out_dir / "rescue_appendix.tsv"
    boundary_path = out_dir / "boundary_ties.tsv"
    excluded_path = out_dir / "shortlist_excluded.tsv"
    primary.to_csv(primary_path, sep="\t", index=False)
    rescue.to_csv(rescue_path, sep="\t", index=False)
    boundary.to_csv(boundary_path, sep="\t", index=False)
    excluded.to_csv(excluded_path, sep="\t", index=False)

    truncated_boundary_tie = bool(len(boundary) and (~boundary["selected_in_primary"]).any())
    manifest = {
        "schema_version": "1.0",
        "stage": stage,
        "input": {"path": str(ranked_path), "sha256": _sha256(ranked_path)},
        "target_size": target_size,
        "primary_count": len(primary),
        "rescue_target_size": rescue_size,
        "rescue_count": len(rescue),
        "priority_column": priority_column,
        "boundary_priority": boundary_rank,
        "boundary_screening_rank": boundary_rank if priority_column == "screening_rank" else None,
        "boundary_tie_count": len(boundary),
        "boundary_tie_truncated": truncated_boundary_tie,
        "selection_policy": (
            f"{priority_column} asc, dry_score desc, support_family_count desc, "
            "independent_source_count desc, gene_id asc; primary size is capped"
        ),
        "eligibility_policy": (
            "dry_group and top_eligibility are retained for interpretation but do not filter Stage A retrieval"
        ),
        "llm_used": False,
        "llm_may_overwrite": False,
        "warning": (
            "A capped boundary tie uses deterministic non-biological tie breaks. "
            "Review boundary_ties.tsv before excluding tied genes."
            if truncated_boundary_tie
            else "No boundary tie was truncated."
        ),
    }
    manifest_path = out_dir / "shortlist_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return {
        "shortlist": str(primary_path),
        "rescue_appendix": str(rescue_path),
        "boundary_ties": str(boundary_path),
        "excluded": str(excluded_path),
        "manifest": str(manifest_path),
        "primary_count": len(primary),
        "rescue_count": len(rescue),
        "boundary_tie_truncated": truncated_boundary_tie,
    }
