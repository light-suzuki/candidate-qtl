"""Biological evidence tiers, kept apart from the Stage-A assay budget (issue #229).

A fixed Stage-A size of 20 is operationally useful -- a bench can only run so
many primers this month -- but 20 candidates mean something completely
different in an 18-gene interval than in an 800-gene interval. Before this
module, ``screening-shortlist`` (issue #144) answered exactly one question
("which candidates are cheapest to review next, given a frozen budget?") and
nothing in the codebase answered the other one ("how strong is the biological
case for each candidate, independent of how many we can afford to look at
right now?"). Reports that only ever showed the Top-20 could not help
conflating the two: a rank-21 candidate reads as *worse* than rank 20 even when
their evidence is identical, purely because of where an arbitrary cut fell.

This module keeps the two questions in two columns:

``evidence_tier``
    ``strong`` / ``plausible`` / ``unresolved`` / ``contradictory`` /
    ``knowledge_limited``. Derived only from :mod:`hardening`'s
    ``top_eligibility`` (itself already built from score contributions,
    support/against direction, and citable-source independence -- see
    ``prioritize_evidence``) and, where supplied, the
    ``knowledge_coverage_tier`` from :mod:`knowledge_coverage`. It never reads
    ``target_size``, ``top_n``, ``reserve_n`` or any other operational cap, so
    the tier of a candidate cannot change just because the assay budget did.

``operational_state``
    ``stage_a_review_now`` / ``reserve`` / ``not_in_current_assay_budget``.
    Purely a function of screening-priority rank against a configurable
    ``review_budget`` (``assay_capacity``), using the exact same ordering rule
    as ``screening-shortlist`` (:func:`shortlist.prepare_screening_order`, so
    the two commands cannot silently drift into different cuts). Default
    ``review_budget`` stays 20 for backwards compatibility.

The two are reported side by side, never merged into one score, and the
manifest says outright: **the operational cut is never a biological
exclusion**. Candidate #20 and candidate #21 are not implied to differ
biologically merely because the assay budget can only afford one of them this
round -- the crosstab in this report is the artifact a reviewer needs to see
that for themselves (their tiers may be identical, and often are).

``evidence_tier_budget_benchmark`` closes the loop on the benchmark side:
Recall@budget (did the truth gene make the Stage-A cut?) is reported next to
tier-based recall (was the truth gene's *biological* tier strong/plausible at
all, budget aside?), and interval-size strata are handed to the existing
:func:`knowledge_coverage.benchmark_stratify` rather than a second
stratification implementation.
"""

from __future__ import annotations

import json

import pandas as pd

from . import __version__
from .core import file_fingerprint, runtime_path
from .knowledge_coverage import benchmark_stratify
from .shortlist import prepare_screening_order

SCHEMA_VERSION = "1.0"

EVIDENCE_TIERS = ("strong", "plausible", "unresolved", "contradictory", "knowledge_limited")
OPERATIONAL_STATES = ("stage_a_review_now", "reserve", "not_in_current_assay_budget")
KNOWLEDGE_LIMITED_COVERAGE_TIERS = ("knowledge_limited", "dark")
INTERVAL_SIZE_BANDS = ("<=50", "51-200", "201-500", ">500")

EVIDENCE_TIER_COLUMNS = [
    "gene_id", "evidence_tier", "tier_basis", "top_eligibility", "knowledge_coverage_tier",
    "dry_score", "screening_rank", "review_budget_rank", "operational_state", "boundary_tie",
    "rank_contribution", "notes",
]
CROSSTAB_COLUMNS = [
    "evidence_tier", "tier_total", "tier_fraction_of_interval",
    "in_budget_count", "in_budget_fraction_of_tier", "reserve_count", "out_of_budget_count",
]
CASE_COLUMNS = [
    "case_id", "review_budget", "candidate_count", "interval_size_band",
    "answer_gene_id", "answer_present_in_candidates", "answer_evidence_tier", "answer_operational_state",
    "recall_at_budget", "tier_is_strong_or_plausible", "notes",
]
STRATA_COLUMNS = ["case_id", "interval_size_band"]


def _tier(top_eligibility: str, knowledge_coverage_tier: str) -> tuple[str, str]:
    """One deterministic, target_size-free rule, applied identically everywhere.

    Reads only categorical fields that already exist upstream of any budget
    cut: ``top_eligibility`` (support/against direction and citable-source
    independence, already resolved by ``prioritize_evidence``) and, where
    supplied, the knowledge-coverage tier. No numeric threshold is introduced
    here; a candidate never moves tier because a review_budget changed.
    """
    value = top_eligibility.strip()
    if value == "conflicting_evidence":
        return "contradictory", "top_eligibility=conflicting_evidence (support and against both present and unresolved)"
    if value == "strong_genetic_evidence":
        return "strong", "top_eligibility=strong_genetic_evidence"
    if value == "independent_experimental_support":
        return "strong", "top_eligibility=independent_experimental_support (>=2 independently sourced channels)"
    if value == "multiple_non_redundant_channels":
        return "plausible", "top_eligibility=multiple_non_redundant_channels (non-redundant support, not yet independently sourced)"
    if knowledge_coverage_tier in KNOWLEDGE_LIMITED_COVERAGE_TIERS:
        return ("knowledge_limited",
                f"knowledge_coverage_tier={knowledge_coverage_tier}: too little external knowledge exists to judge this candidate further, not weak evidence")
    return "unresolved", "top_eligibility=insufficient_independent_evidence and no contradiction and no knowledge-coverage input to explain it"


def _interval_size_band(count: int) -> str:
    if count <= 50:
        return "<=50"
    if count <= 200:
        return "51-200"
    if count <= 500:
        return "201-500"
    return ">500"


def _tier_frame(
    ranked_path,
    review_budget: int,
    reserve_size: int,
    knowledge_coverage_path=None,
) -> pd.DataFrame:
    frame = pd.read_csv(ranked_path, sep="\t", dtype=str).fillna("")
    if "top_eligibility" not in frame.columns:
        raise ValueError("ranked candidate table must contain top_eligibility for tier derivation (run evidence-rank first)")

    coverage_by_gene: dict[str, str] = {}
    if knowledge_coverage_path is not None:
        coverage = pd.read_csv(knowledge_coverage_path, sep="\t", dtype=str).fillna("")
        if "gene_id" not in coverage.columns or "knowledge_coverage_tier" not in coverage.columns:
            raise ValueError("knowledge-coverage table must contain gene_id and knowledge_coverage_tier")
        coverage_by_gene = {str(row["gene_id"]).casefold(): str(row["knowledge_coverage_tier"])
                            for row in coverage.to_dict("records")}

    ordered = prepare_screening_order(frame)
    ordered = ordered.reset_index(drop=True)
    ordered["review_budget_rank"] = range(1, len(ordered) + 1)

    budget_cut = min(review_budget, len(ordered))
    primary_ids = set(ordered.head(budget_cut)["gene_id"])
    if budget_cut and budget_cut <= len(ordered):
        boundary_rank = int(ordered.iloc[budget_cut - 1]["_screening_rank"])
    else:
        boundary_rank = None
    boundary_tie_mask = ordered["_screening_rank"].eq(boundary_rank) if boundary_rank is not None else pd.Series(False, index=ordered.index)

    reserve_ids: set[str] = set()
    if reserve_size > 0 and budget_cut < len(ordered):
        reserve_ids = set(ordered.iloc[budget_cut:budget_cut + reserve_size]["gene_id"])

    def _operational_state(gene_id: str) -> str:
        if gene_id in primary_ids:
            return "stage_a_review_now"
        if gene_id in reserve_ids:
            return "reserve"
        return "not_in_current_assay_budget"

    tiers = ordered.apply(
        lambda row: _tier(str(row["top_eligibility"]), coverage_by_gene.get(str(row["gene_id"]).casefold(), "")),
        axis=1,
    )
    ordered["evidence_tier"] = [item[0] for item in tiers]
    ordered["tier_basis"] = [item[1] for item in tiers]
    ordered["knowledge_coverage_tier"] = ordered["gene_id"].map(
        lambda gene_id: coverage_by_gene.get(str(gene_id).casefold(), ""))
    ordered["operational_state"] = ordered["gene_id"].map(_operational_state)
    ordered["boundary_tie"] = boundary_tie_mask
    ordered.attrs["boundary_rank"] = boundary_rank
    ordered.attrs["budget_cut"] = budget_cut
    ordered.attrs["primary_ids"] = primary_ids
    return ordered


def evidence_tier_report(
    ranked_arg: str,
    out_dir_arg: str,
    review_budget: int = 20,
    reserve_size: int = 0,
    knowledge_coverage_arg: str | None = None,
    case_id: str = "",
    full_hash: bool = False,
) -> dict[str, object]:
    """Report biological evidence tiers and the operational assay-budget cut, side by side.

    ``review_budget`` (assay_capacity) defaults to 20, matching
    ``screening-shortlist``'s ``target_size`` default, so existing callers who
    do not pass it see the same Stage-A cut. Raising or lowering it changes
    only ``operational_state``; ``evidence_tier`` never moves.
    """
    if review_budget < 1:
        raise ValueError("review_budget must be at least 1")
    if reserve_size < 0:
        raise ValueError("reserve_size must be zero or greater")

    ranked_path = runtime_path(ranked_arg)
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    knowledge_coverage_path = runtime_path(knowledge_coverage_arg) if knowledge_coverage_arg else None

    ordered = _tier_frame(ranked_path, review_budget, reserve_size, knowledge_coverage_path)
    boundary_rank = ordered.attrs["boundary_rank"]
    budget_cut = ordered.attrs["budget_cut"]
    primary_ids = ordered.attrs["primary_ids"]

    ordered["dry_score"] = ordered["_dry_score"]
    ordered["screening_rank"] = ordered["_screening_rank"]
    ordered["rank_contribution"] = "none"
    ordered["notes"] = (
        "evidence_tier is target_size-independent; operational_state alone reflects review_budget "
        "and is never a biological exclusion"
    )
    report = ordered[EVIDENCE_TIER_COLUMNS].copy()
    report_path = out / "evidence_tiers.tsv"
    report.to_csv(report_path, sep="\t", index=False, lineterminator="\n")

    # Candidates immediately outside the cut: the truncated tie group if the
    # boundary score was split by the cut, otherwise simply the next-ranked row.
    if budget_cut >= len(ordered):
        just_outside = report.iloc[0:0]
    else:
        tie_excluded = report[ordered["boundary_tie"] & ~ordered["gene_id"].isin(primary_ids)]
        just_outside = tie_excluded if len(tie_excluded) else report.iloc[budget_cut:budget_cut + 1]
    just_outside_path = out / "just_outside_budget.tsv"
    just_outside.to_csv(just_outside_path, sep="\t", index=False, lineterminator="\n")

    total = len(report)
    crosstab_rows: list[dict[str, object]] = []
    for tier in EVIDENCE_TIERS:
        subset = report[report["evidence_tier"] == tier]
        tier_total = int(len(subset))
        in_budget = int((subset["operational_state"] == "stage_a_review_now").sum())
        reserve = int((subset["operational_state"] == "reserve").sum())
        crosstab_rows.append({
            "evidence_tier": tier, "tier_total": tier_total,
            "tier_fraction_of_interval": round(tier_total / total, 6) if total else "",
            "in_budget_count": in_budget,
            "in_budget_fraction_of_tier": round(in_budget / tier_total, 6) if tier_total else "",
            "reserve_count": reserve,
            "out_of_budget_count": tier_total - in_budget - reserve,
        })
    crosstab = pd.DataFrame(crosstab_rows, columns=CROSSTAB_COLUMNS)
    crosstab_path = out / "evidence_tier_budget_crosstab.tsv"
    crosstab.to_csv(crosstab_path, sep="\t", index=False, lineterminator="\n")

    # Concretely check, rather than merely assert, whether the row at the edge
    # of the budget and the row right after it differ biologically.
    boundary_tier_comparison: dict[str, object] = {"evaluated": False}
    if budget_cut and budget_cut < len(report):
        at_cut = report.iloc[budget_cut - 1]
        after_cut = report.iloc[budget_cut]
        boundary_tier_comparison = {
            "evaluated": True,
            "rank_at_budget": int(at_cut["review_budget_rank"]), "rank_at_budget_tier": str(at_cut["evidence_tier"]),
            "rank_after_budget": int(after_cut["review_budget_rank"]), "rank_after_budget_tier": str(after_cut["evidence_tier"]),
            "same_biological_tier": bool(at_cut["evidence_tier"] == after_cut["evidence_tier"]),
        }

    tier_counts = {tier: int((report["evidence_tier"] == tier).sum()) for tier in EVIDENCE_TIERS}
    operational_counts = {state: int((report["operational_state"] == state).sum()) for state in OPERATIONAL_STATES}
    manifest = {
        "schema_version": SCHEMA_VERSION, "command": "evidence-tiers", "candidate_qtl_version": __version__,
        "case_id": case_id,
        "candidate_count": total,
        "review_budget": review_budget, "reserve_size": reserve_size,
        "interval_coverage_fraction_in_budget": round(budget_cut / total, 6) if total else "",
        "evidence_tiers": list(EVIDENCE_TIERS), "operational_states": list(OPERATIONAL_STATES),
        "tier_counts": tier_counts, "operational_state_counts": operational_counts,
        "boundary": {
            "screening_rank_at_cut": boundary_rank,
            "boundary_tie_row_count": int(ordered["boundary_tie"].sum()),
            "boundary_tie_truncated": bool((ordered["boundary_tie"] & ~ordered["gene_id"].isin(primary_ids)).any()),
        },
        "boundary_tier_comparison": boundary_tier_comparison,
        "knowledge_coverage_provided": bool(knowledge_coverage_arg),
        "rank_contribution": "none",
        "rules": [
            "evidence_tier never reads target_size, review_budget, top_n or reserve_n: it is the same value at any budget",
            "operational_state is purely a rank-against-review_budget cut using the same ordering as screening-shortlist",
            "the operational cut is never a biological exclusion; a candidate outside the assay budget keeps its evidence_tier",
            "candidate #20 and candidate #21 are not implied to differ biologically merely because one is inside "
            "review_budget and the other is not -- see boundary_tier_comparison for the actual answer on this run",
            "knowledge_limited describes how little could be said, not how weak the candidate is; it is not a demotion",
            "contradictory means top_eligibility already detected a genuine support/against conflict, not merely low evidence",
        ],
        "inputs": {"ranked": file_fingerprint(ranked_path, full_hash),
                   **({"knowledge_coverage": file_fingerprint(knowledge_coverage_path, full_hash)} if knowledge_coverage_path else {})},
        "outputs": {"report": str(report_path), "just_outside_budget": str(just_outside_path),
                    "crosstab": str(crosstab_path)},
    }
    manifest_path = out / "evidence_tiers_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), "report": str(report_path),
            "just_outside_budget": str(just_outside_path), "crosstab": str(crosstab_path),
            "candidate_count": total, "tier_counts": tier_counts, "operational_state_counts": operational_counts,
            "rank_contribution": "none"}


def evidence_tier_budget_benchmark(
    manifest_arg: str,
    out_dir_arg: str,
    review_budget: int = 20,
    full_hash: bool = False,
) -> dict[str, object]:
    """Recall@budget next to tier-based recall, stratified by interval-size band.

    Reads a manifest of frozen ``evidence-rank``/``candidate-pipeline`` outputs
    (one row per case: ``case_id``, ``ranked_tsv``, optional ``answer_gene_id``
    and ``knowledge_coverage_tsv``) and never re-runs ranking. Stratified
    reporting is delegated to :func:`knowledge_coverage.benchmark_stratify`
    rather than reimplemented here, so the per-stratum denominator rule stays
    in one place.
    """
    manifest_path = runtime_path(manifest_arg)
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str).fillna("")
    required = {"case_id", "ranked_tsv"}
    missing = sorted(required - set(manifest.columns))
    if missing:
        raise ValueError(f"evidence-tier-benchmark manifest missing columns: {', '.join(missing)}")
    if manifest.empty:
        raise ValueError("evidence-tier-benchmark manifest is empty")
    if manifest["case_id"].duplicated().any():
        raise ValueError("evidence-tier-benchmark case_id must be unique")

    rows: list[dict[str, object]] = []
    for record in manifest.to_dict("records"):
        case_id = str(record["case_id"]).strip()
        ranked_path = runtime_path(str(record["ranked_tsv"]).strip())
        knowledge_coverage_value = str(record.get("knowledge_coverage_tsv", "")).strip()
        knowledge_coverage_path = runtime_path(knowledge_coverage_value) if knowledge_coverage_value else None
        answer_gene_id = str(record.get("answer_gene_id", "")).strip()

        ordered = _tier_frame(ranked_path, review_budget, reserve_size=0, knowledge_coverage_path=knowledge_coverage_path)
        candidate_count = int(len(ordered))
        by_gene = {str(row["gene_id"]).casefold(): row for row in ordered.to_dict("records")}
        answer_row = by_gene.get(answer_gene_id.casefold()) if answer_gene_id else None

        if not answer_gene_id:
            recall_at_budget = ""
            tier_hit = ""
            answer_present = ""
            answer_tier = ""
            answer_state = ""
        elif answer_row is None:
            # A named truth gene that is not among this case's candidates at
            # all is a real, informative miss for Recall@budget -- not a
            # missing-input case -- but its evidence_tier cannot be read from a
            # row that does not exist, so the tier metric is not_evaluated.
            recall_at_budget = "False"
            tier_hit = ""
            answer_present = "False"
            answer_tier = ""
            answer_state = ""
        else:
            answer_present = "True"
            answer_tier = str(answer_row["evidence_tier"])
            answer_state = str(answer_row["operational_state"])
            recall_at_budget = str(answer_state == "stage_a_review_now")
            tier_hit = str(answer_tier in ("strong", "plausible"))

        rows.append({
            "case_id": case_id, "review_budget": review_budget, "candidate_count": candidate_count,
            "interval_size_band": _interval_size_band(candidate_count),
            "answer_gene_id": answer_gene_id, "answer_present_in_candidates": answer_present,
            "answer_evidence_tier": answer_tier, "answer_operational_state": answer_state,
            "recall_at_budget": recall_at_budget, "tier_is_strong_or_plausible": tier_hit,
            "notes": "recall_at_budget divides by cases with a named answer_gene_id; a named gene absent from "
                     "the candidate table is a real miss, not_evaluated is reserved for a missing answer_gene_id",
        })

    cases = pd.DataFrame(rows, columns=CASE_COLUMNS)
    cases_path = out / "evidence_tier_benchmark_cases.tsv"
    cases.to_csv(cases_path, sep="\t", index=False, lineterminator="\n")
    strata_path = out / "evidence_tier_benchmark_strata.tsv"
    cases[STRATA_COLUMNS].to_csv(strata_path, sep="\t", index=False, lineterminator="\n")

    total = len(cases)

    def _stats(column: str) -> dict[str, object]:
        values = cases[column].astype(str).str.strip()
        evaluated = int(values.ne("").sum())
        hits = int(values.eq("True").sum())
        return {"evaluated": evaluated, "hits": hits,
                "recall": round(hits / evaluated, 6) if evaluated else None, "not_evaluated": total - evaluated}

    metrics = {"recall_at_budget": _stats("recall_at_budget"), "tier_is_strong_or_plausible": _stats("tier_is_strong_or_plausible")}

    stratified = benchmark_stratify(
        str(cases_path), str(strata_path), str(out),
        metric_args=["recall_at_budget", "tier_is_strong_or_plausible"],
        stratum_column_args=["interval_size_band"],
    )

    manifest_out = {
        "schema_version": SCHEMA_VERSION, "command": "evidence-tier-benchmark", "candidate_qtl_version": __version__,
        "case_count": total, "review_budget": review_budget,
        "metrics": metrics,
        "interval_size_bands": list(INTERVAL_SIZE_BANDS),
        "stratified": stratified,
        "denominator_policy": (
            "recall_at_budget divides only by cases with a named answer_gene_id; tier_is_strong_or_plausible "
            "additionally requires the answer_gene_id to be present among that case's candidates"
        ),
        "rules": [
            "Recall@budget answers whether the truth gene made the frozen Stage-A cut, an operational question",
            "tier_is_strong_or_plausible answers whether the truth gene's biological evidence tier was strong or "
            "plausible at all, independent of any budget -- the two numbers are reported separately and never averaged",
            "interval-size strata are recomputed through knowledge_coverage.benchmark_stratify, not a second implementation",
        ],
        "inputs": {"manifest": file_fingerprint(manifest_path, full_hash)},
        "outputs": {"cases": str(cases_path), "strata": str(strata_path)},
    }
    manifest_path_out = out / "evidence_tier_benchmark_manifest.json"
    manifest_path_out.write_text(json.dumps(manifest_out, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path_out), "cases": str(cases_path), "strata": str(strata_path),
            "case_count": total, "metrics": metrics, "stratified": stratified}
