"""Tie-aware ranking-quality metrics for the candidate-first benchmarks.

The map-based Top-20 recall answers "given a legitimate cross-species pathway
homolog as bait, does the causal gene survive to the shortlist?". It does not,
on its own, say how good the *first* gene to sequence is, whether a confident
single call was actually the wrong paralog, or whether the system safely held
when the truth was missing. This module makes those separate, first-class
metrics (issue #154):

* Top-1 / Top-3 / Top-5 / Top-10 / Top-20 recall, reported both ``conservative``
  (the whole score-tie group must fit within k -- uses ``rank_max``) and
  ``optimistic`` (the answer's best possible rank -- uses ``rank_min``).
* reciprocal rank / MRR, again conservative and optimistic.
* tie-aware ``rank_min`` / ``rank_max`` / ``rank_percentile``.
* ``wrong_paralog_sole_top`` -- a single non-answer gene alone at the top score
  while the true gene is present in the interval (the pea_Tl failure mode).
* ``safe_hold`` / ``causal_tied_top`` -- the causal gene is not the sole top, but
  the system did not confidently commit to a single wrong gene either.
* for truth-missing / hard-negative cases: ``false_confident_sole_top`` vs
  ``abstained``.

Every value is computed from ``dry_score`` (and ``dry_group`` for the confidence
gate) only. Nothing here depends on gene-id order or any arbitrary tie break, so
the metrics are stable under row shuffling. This module never changes weights or
ranks; metric improvement must come from new signal, never from fitting weights
to a benchmark answer.

Issue #350 adds the **denominator engine** for family / route recall:
:func:`aggregate_denominator_metrics` builds each metric's denominator from
per-case states (``hit`` / ``miss`` / ``not_evaluated`` / ``not_scoreable`` /
``review_required`` / ``contaminated``). Only ``hit`` and ``miss`` form the
eligible denominator; every other state is excluded from it and reported
explicitly -- a case that was not_evaluated or not_scoreable is *never*
converted into a miss. Denominators are never hard-coded: each metric carries
its own ``eligible`` / ``hit`` / ``miss`` / ``not_evaluated`` / ``not_scoreable``
/ ``review_required`` / ``contaminated`` counts.
"""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

from .benchmark_reporting import continuous_bootstrap_ci

TOP_K_LEVELS = (1, 3, 5, 10, 20)

# Coarse, mutually exclusive outcome per case. Kept small so the summary reads
# cleanly; the individual booleans below carry the nuance.
RESOLUTION_OUTCOMES = (
    "causal_sole_top",       # answer present, unique top -> confident and correct
    "causal_tied_top",       # answer present, tied at the top score with others
    "wrong_sole_top",        # answer present, a non-answer is alone at the top
    "safe_hold",             # answer present, not top, but no single wrong gene committed
    "false_confident_sole_top",  # truth missing, a single gene confidently promoted
    "abstained",             # truth missing, no confident single call
)


def _answer_ids(value: object) -> list[str]:
    return [item.strip() for item in str(value).split("|") if item.strip()]


def resolution_metrics(ranked: pd.DataFrame, accepted: object) -> dict[str, object]:
    """Compute tie-aware ranking-quality metrics for one case.

    ``ranked`` needs a ``gene_id`` and ``dry_score`` column; ``dry_group`` is used
    when present to tell a confident promotion (``dry_group == 'top'``) from a
    gene that merely leads on score. ``accepted`` is the pipe-joined answer id
    string (or a list) of the genes that count as the causal answer; an empty /
    absent answer marks a truth-missing case.
    """
    accepted_set = set(accepted if isinstance(accepted, (list, set, tuple)) else _answer_ids(accepted))
    genes = ranked["gene_id"].astype(str).reset_index(drop=True)
    score = pd.to_numeric(ranked["dry_score"], errors="coerce").fillna(0.0).reset_index(drop=True)
    group = (
        ranked["dry_group"].astype(str).reset_index(drop=True)
        if "dry_group" in ranked.columns
        else pd.Series([""] * len(genes))
    )
    n = int(len(genes))

    present_mask = genes.isin(accepted_set)
    answer_present = bool(present_mask.any())
    max_score = float(score.max()) if n else 0.0
    # A gene with no positive evidence is not "recovered" -- it is merely present
    # in the interval. When nothing scores above zero the route has no signal, so
    # there is no top group and every candidate is an abstention, not a tie at top.
    has_top_signal = bool(n) and max_score > 0
    top_mask = (score == max_score) if has_top_signal else pd.Series([False] * n)
    tie_size_at_top = int(top_mask.sum())
    sole_top = tie_size_at_top == 1
    sole_top_gene_id = str(genes[top_mask].iloc[0]) if sole_top else ""
    # A confident single call: alone at the top score AND actually promoted to the
    # Top group. A sole-by-score gene that the gate left in reserve / insufficient
    # is a hold, not a confident call. Without a dry_group column, treat a sole
    # top as a confident call.
    if sole_top:
        top_group = str(group[top_mask].iloc[0])
        confident_sole_top = (top_group == "top") if "dry_group" in ranked.columns else True
    else:
        confident_sole_top = False
    causal_in_top_group = bool((top_mask & present_mask).any())

    if answer_present:
        answer_score = float(score[present_mask].max())
        rank_min = int((score > answer_score).sum()) + 1
        rank_max = int((score >= answer_score).sum())
        tie_size = int((score == answer_score).sum())
        # deterministic label: the best-scoring accepted gene, tie-broken by id so
        # it never flips on row order (metrics below never use this label)
        matched = sorted(str(g) for g in genes[present_mask & (score == answer_score)])
        matched_answer_gene_id = matched[0] if matched else ""
        has_positive_answer = answer_score > 0
        reciprocal_rank_conservative = 1.0 / rank_max if has_positive_answer else 0.0
        reciprocal_rank_optimistic = 1.0 / rank_min if has_positive_answer else 0.0
        rank_percentile = (n - rank_max + 1) / n if (n and has_positive_answer) else 0.0
    else:
        answer_score = None
        rank_min = rank_max = None
        tie_size = 0
        matched_answer_gene_id = ""
        has_positive_answer = False
        reciprocal_rank_conservative = reciprocal_rank_optimistic = 0.0
        rank_percentile = 0.0

    sole_top_is_causal = sole_top and sole_top_gene_id in accepted_set
    wrong_paralog_sole_top = bool(answer_present and sole_top and not sole_top_is_causal)
    tied_top_excludes_answer = bool(answer_present and tie_size_at_top > 1 and not causal_in_top_group)
    false_confident_sole_top = bool(not answer_present and confident_sole_top)
    abstained = bool(not answer_present and not confident_sole_top)

    # Confident wrong sole top is the failure to catch: a non-answer promoted to a
    # single Top call. A wrong paralog that merely leads on score but was left out
    # of the Top group (the gate held) is a safe hold, not a confident error --
    # detecting and holding is the intended behaviour, not tuning its weight down.
    confident_wrong_sole_top = bool(wrong_paralog_sole_top and confident_sole_top)
    if not answer_present:
        outcome = "false_confident_sole_top" if false_confident_sole_top else "abstained"
        safe_hold = False
    elif sole_top_is_causal:
        outcome = "causal_sole_top"
        safe_hold = False
    elif causal_in_top_group:
        outcome = "causal_tied_top"
        safe_hold = False
    elif confident_wrong_sole_top:
        outcome = "wrong_sole_top"
        safe_hold = False
    else:
        # answer present but not top, and no single wrong gene was confidently
        # promoted (a wrong paralog led on score but was held, or the top is a tie).
        outcome = "safe_hold"
        safe_hold = True

    metrics: dict[str, object] = {
        "answer_present": answer_present,
        "answer_present_count": int(present_mask.sum()),
        "matched_answer_gene_id": matched_answer_gene_id,
        "answer_score": answer_score,
        "candidate_count": n,
        "rank_min": rank_min,
        "rank_max": rank_max,
        "tie_size": tie_size,
        "is_tied_at_answer_rank": bool(answer_present and tie_size > 1),
        "tie_size_at_top": tie_size_at_top,
        "sole_top": sole_top,
        "sole_top_gene_id": sole_top_gene_id,
        "confident_sole_top": confident_sole_top,
        "causal_in_top_group": causal_in_top_group,
        "reciprocal_rank_conservative": round(reciprocal_rank_conservative, 6),
        "reciprocal_rank_optimistic": round(reciprocal_rank_optimistic, 6),
        "rank_percentile": round(rank_percentile, 6),
        "wrong_paralog_sole_top": wrong_paralog_sole_top,
        "confident_wrong_sole_top": confident_wrong_sole_top,
        "tied_top_excludes_answer": tied_top_excludes_answer,
        "false_confident_sole_top": false_confident_sole_top,
        "abstained": abstained,
        "safe_hold": safe_hold,
        "resolution_outcome": outcome,
    }
    for k in TOP_K_LEVELS:
        metrics[f"top{k}_conservative"] = bool(rank_max is not None and rank_max <= k and has_positive_answer)
        metrics[f"top{k}_optimistic"] = bool(rank_min is not None and rank_min <= k and has_positive_answer)
    return metrics


def aggregate_resolution_metrics(
    rows: list[dict[str, object]],
    *,
    mrr_bootstrap_iterations: int = 2000,
    mrr_bootstrap_seed: int = 20260812,
) -> dict[str, object]:
    """Aggregate per-case metrics into a benchmark summary.

    Recall / MRR are computed over the answer-present cases; abstention and
    false-confidence are computed over the truth-missing cases, so a run's
    Top-1 quality, its Top-20 recall, and its behaviour when the truth is absent
    stay separate numbers instead of one headline.

    Issue #367: MRR gains a deterministic case-bootstrap 95% CI
    (``continuous_bootstrap_ci`` over the answer-present rows' reciprocal
    ranks, the same population as the headline MRR). Defaults preserve
    existing callers.
    """
    present = [row for row in rows if row.get("answer_present")]
    missing = [row for row in rows if not row.get("answer_present")]
    n_present = len(present)
    n_missing = len(missing)

    def _rate(subset: list[dict[str, object]], key: str) -> float:
        return round(sum(1 for row in subset if row.get(key)) / len(subset), 6) if subset else 0.0

    def _mean(subset: list[dict[str, object]], key: str) -> float:
        return round(sum(float(row.get(key, 0.0) or 0.0) for row in subset) / len(subset), 6) if subset else 0.0

    summary: dict[str, object] = {
        "n_cases": len(rows),
        "n_answer_present": n_present,
        "n_truth_missing": n_missing,
        "mrr_conservative": _mean(present, "reciprocal_rank_conservative"),
        "mrr_optimistic": _mean(present, "reciprocal_rank_optimistic"),
    }
    mrr_conservative_ci = continuous_bootstrap_ci(
        [float(row.get("reciprocal_rank_conservative", 0.0) or 0.0) for row in present],
        iterations=mrr_bootstrap_iterations, seed=mrr_bootstrap_seed,
    )
    mrr_optimistic_ci = continuous_bootstrap_ci(
        [float(row.get("reciprocal_rank_optimistic", 0.0) or 0.0) for row in present],
        iterations=mrr_bootstrap_iterations, seed=mrr_bootstrap_seed,
    )
    summary["mrr_conservative_ci_95_low"] = mrr_conservative_ci["ci_95_low"]
    summary["mrr_conservative_ci_95_high"] = mrr_conservative_ci["ci_95_high"]
    summary["mrr_optimistic_ci_95_low"] = mrr_optimistic_ci["ci_95_low"]
    summary["mrr_optimistic_ci_95_high"] = mrr_optimistic_ci["ci_95_high"]
    summary.update({
        # raw = a non-answer alone at the top score (archive-continuous "wrong sole
        # top"); confident = the subset actually promoted to Top (the safety failure).
        "wrong_sole_top_count": sum(1 for row in present if row.get("wrong_paralog_sole_top")),
        "wrong_sole_top_rate": _rate(present, "wrong_paralog_sole_top"),
        "confident_wrong_sole_top_count": sum(1 for row in present if row.get("confident_wrong_sole_top")),
        "confident_wrong_sole_top_rate": _rate(present, "confident_wrong_sole_top"),
        "safe_hold_rate": _rate(present, "safe_hold"),
        "causal_sole_top_rate": round(sum(1 for row in present if row.get("resolution_outcome") == "causal_sole_top") / n_present, 6) if n_present else 0.0,
        "causal_tied_top_rate": round(sum(1 for row in present if row.get("resolution_outcome") == "causal_tied_top") / n_present, 6) if n_present else 0.0,
        "tied_top_excludes_answer_count": sum(1 for row in present if row.get("tied_top_excludes_answer")),
        "abstention_rate": _rate(missing, "abstained"),
        "false_confident_rate": _rate(missing, "false_confident_sole_top"),
        "median_candidate_count": float(pd.Series([row.get("candidate_count", 0) for row in rows]).median()) if rows else 0.0,
    })
    for k in TOP_K_LEVELS:
        summary[f"top{k}_recall_conservative"] = _rate(present, f"top{k}_conservative")
        summary[f"top{k}_recall_optimistic"] = _rate(present, f"top{k}_optimistic")
    return summary


# Issue #350 / #412: the denominator vocabulary shared by every route-family metric.
# Only hit and miss are eligible -- every other state is excluded from the
# denominator and reported in its own bucket; a case that was not_evaluated or
# not_scoreable or not_evaluated_not_in_index is never converted into a miss.
DENOMINATOR_STATES = ("hit", "miss", "not_evaluated", "not_evaluated_not_in_index",
                      "not_scoreable", "review_required", "contaminated")


def aggregate_denominator_metrics(
    rows: list[dict[str, object]], state_key: str = "state"
) -> dict[str, object]:
    """Build one metric's denominator from per-case states (issue #350, #412).

    ``rows`` are per-case records whose ``state_key`` value is one of
    :data:`DENOMINATOR_STATES`. The eligible denominator is ``hit + miss``;
    ``recall = hit / (hit + miss)`` and is ``None`` when the denominator is
    zero. ``not_evaluated`` / ``not_evaluated_not_in_index`` / ``not_scoreable`` /
    ``review_required`` / ``contaminated`` cases are counted in their own buckets
    and are never folded into the miss count or the denominator.

    The denominator is never hard-coded: every caller supplies the per-case
    states, and the requirements that produced them (family_scoreable, route
    evaluability, cleanliness, rank eligibility) stay in the per-case table
    that feeds this function.
    """
    state_counts = {state: sum(1 for row in rows if row.get(state_key) == state)
                    for state in DENOMINATOR_STATES}
    denominator = state_counts["hit"] + state_counts["miss"]
    return {
        "eligible": denominator,
        "hit": state_counts["hit"],
        "miss": state_counts["miss"],
        "not_evaluated": state_counts["not_evaluated"],
        "not_evaluated_not_in_index": state_counts["not_evaluated_not_in_index"],
        "not_scoreable": state_counts["not_scoreable"],
        "review_required": state_counts["review_required"],
        "contaminated": state_counts["contaminated"],
        "denominator": denominator,
        "recall": round(state_counts["hit"] / denominator, 6) if denominator else None,
    }


# --------------------------------------------------------------------------- #
# Issue #404: blind end-to-end benchmark contract (lane3 scope: the contract
# and schema are fixed here; no benchmark run is executed by this lane).
#
# The headline claim is "blind end-to-end Top-20 recall >= 0.90 while achieving
# substantial candidate reduction and a low wrong-confident-promotion rate".
# The contract makes the prohibitions of issue #404 structurally impossible
# instead of merely documented:
#   - Top-20-only reporting is impossible: top5/top10 columns are mandatory
#     whenever any Top-k column is reported;
#   - a reduction-less Top-20 headline is impossible: the interval/narrowed
#     counts and reduction_ratio are mandatory alongside every Top-k column;
#   - the four metric families are separate and never collapsed into one
#     aggregate score;
#   - preferred_copy stays an optional diagnostic (never auto-miss);
#   - the blind sequence (prediction -> freeze -> truth open -> score) is
#     mechanically verified through a manifest / hash contract;
#   - the target thresholds are fixed in the manifest before any result is
#     seen and may not be tuned after results.
# --------------------------------------------------------------------------- #

BENCHMARK_CONTRACT_SCHEMA_VERSION = "1.0"

# Issue #404 metric families: each answers a different question and is reported
# separately. Collapsing these into one aggregate score is a contract violation.
BENCHMARK_METRIC_FAMILIES = (
    ("phenotype_to_family", "phenotype -> proposed mechanism/function/family recall"),
    ("source_verified_family", "source/accession-verified family recall"),
    ("family_to_wet_candidate", "family -> Wet candidate set retention"),
    ("end_to_end_topk", "blind end-to-end truth gene retention in Top-5/10/20"),
)

# The mandatory report columns of a #404-contract benchmark. ``not_evaluated``
# is reason-broken (the not_evaluated_reasons column carries per-case reason
# counts); ``coverage`` is the evaluated/total read-out that must accompany
# every headline.
BENCHMARK_REQUIRED_REPORT_COLUMNS = (
    "top5_recall", "top10_recall", "top20_recall",
    "interval_candidate_count", "narrowed_candidate_count", "reduction_ratio",
    "family_recall", "wrong_confident_promotion",
    "not_evaluated_count", "not_evaluated_reasons", "coverage",
)

# Top-k columns must always be reported together; a lone top20 column is
# structurally impossible (issue #404: "Top-20 のみを報告することは禁止").
BENCHMARK_TOP_K_LEVELS = ("top5_recall", "top10_recall", "top20_recall")

# Issue #404: target thresholds fixed in the benchmark manifest BEFORE any
# result is seen. Tuning these after a result is a contract violation
# (benchmark discipline, AGENTS.md §4).
BENCHMARK_TARGET_THRESHOLDS = {
    "top20_recall_min": 0.90,
    "median_reduction_ratio_min": 0.80,
    "wrong_confident_promotion_max": 0.05,
}

# The blind sequence (AGENTS.md §5). A manifest carrying this phase ledger
# proves mechanically that prediction artifacts were frozen before the truth
# was opened and that scoring ran only after the truth was opened.
BLIND_SEQUENCE_PHASES = ("prediction", "freeze", "truth_open", "score")


def validate_benchmark_report_columns(columns: Sequence[str]) -> list[str]:
    """Issue #404: structurally enforce the report column contract.

    ``columns`` is the column list of a per-case or summary report table.
    Returns a list of human-readable violations (empty when the table may be
    published). Top-20-only reporting, reduction-less Top-20 headlines and a
    report with no Top-k column at all are rejected here, so the contract is
    enforced by schema, not by prose.
    """
    violations: list[str] = []
    present = set(columns)
    missing_top_k = [level for level in BENCHMARK_TOP_K_LEVELS if level not in present]
    if missing_top_k:
        violations.append(
            "Top-k columns are mandatory; missing: " + ", ".join(missing_top_k)
        )
    for required in ("interval_candidate_count", "narrowed_candidate_count",
                     "reduction_ratio", "family_recall", "wrong_confident_promotion",
                     "not_evaluated_count", "not_evaluated_reasons", "coverage"):
        if required not in present:
            violations.append(f"required report column missing: {required}")
    return violations


def validate_benchmark_thresholds(thresholds: dict[str, object]) -> list[str]:
    """Issue #404: the manifest thresholds must be exactly the fixed targets.

    The targets are fixed before any result is seen (they are part of the
    frozen manifest); a manifest carrying different values is a contract
    violation, not a new decision.
    """
    violations: list[str] = []
    for key, expected in BENCHMARK_TARGET_THRESHOLDS.items():
        actual = thresholds.get(key)
        if actual != expected:
            violations.append(
                f"target threshold {key}: expected {expected}, manifest has {actual!r}"
            )
    return violations


def verify_blind_sequence(manifest: dict[str, object]) -> list[str]:
    """Issue #404: verify the blind prediction sequence ledger.

    The manifest must carry a ``blind_sequence`` ledger covering every phase of
    :data:`BLIND_SEQUENCE_PHASES` (``prediction`` / ``freeze`` / ``truth_open``
    / ``score``), and the ``freeze`` phase must record the prediction artifact
    hashes that were fixed before the truth was opened. The ledger's own
    timestamps are attestations supplied by the run; the verifier checks
    presence and the freeze-hash requirement mechanically, not the ordering
    implied by free-text timestamps (a JSON object has no key order, so a
    dict-based ledger cannot itself prove temporal order).
    """
    violations: list[str] = []
    ledger = manifest.get("blind_sequence")
    if not isinstance(ledger, dict):
        return ["manifest has no blind_sequence ledger; the blind run order is unverified"]
    missing_phases = [phase for phase in BLIND_SEQUENCE_PHASES if phase not in ledger]
    if missing_phases:
        violations.append("blind_sequence missing phase(s): " + ", ".join(missing_phases))
        return violations
    freeze = ledger.get("freeze")
    if not isinstance(freeze, dict) or not freeze.get("hashes"):
        violations.append("blind_sequence.freeze must record prediction artifact hashes")
    return violations


def aggregate_contract_report(
    per_case_rows: list[dict[str, object]],
    threshold_rows: dict[str, object] | None = None,
) -> dict[str, object]:
    """Issue #404: assemble a contract report from per-case rows.

    ``per_case_rows`` must already carry the mandatory report columns
    (validate via :func:`validate_benchmark_report_columns`). The report keeps
    the four metric families separate -- it never collapses phenotype->family,
    source-verified family, family->Wet candidate and end-to-end Top-k into a
    single aggregate. ``preferred_copy`` results are optional diagnostics and
    are never turned into an automatic miss (a case without a preferred copy
    stays evaluable on the family and Top-k lanes).
    """
    rows = [row for row in per_case_rows if isinstance(row, dict)]
    n = len(rows)
    if n == 0:
        return {
            "case_count": 0,
            "top5_recall": None, "top10_recall": None, "top20_recall": None,
            "median_interval_candidate_count": None, "median_narrowed_candidate_count": None,
            "median_reduction_ratio": None, "family_recall": None,
            "wrong_confident_promotion_rate": None, "not_evaluated_count": 0,
            "coverage": None, "metric_families": {},
            "preferred_copy_policy": (
                "preferred_copy is an optional diagnostic; an undecided preferred copy is never "
                "an automatic benchmark miss (issue #404)"
            ),
        }
    top_k_rates: dict[str, float | None] = {}
    for level in BENCHMARK_TOP_K_LEVELS:
        eligible = [row for row in rows if row.get(level) is not None]
        top_k_rates[level] = (
            round(sum(1 for row in eligible if row.get(level)) / len(eligible), 6)
            if eligible else None
        )
    intervals = [float(row.get("interval_candidate_count", 0.0) or 0.0) for row in rows]
    narrowed = [float(row.get("narrowed_candidate_count", 0.0) or 0.0) for row in rows]
    reduction = [float(row.get("reduction_ratio", 0.0) or 0.0) for row in rows]
    family_eligible = [row for row in rows if row.get("family_recall") is not None]
    family_recall = (round(sum(1 for row in family_eligible if row.get("family_recall"))
                           / len(family_eligible), 6) if family_eligible else None)
    wrong_eligible = [row for row in rows if row.get("wrong_confident_promotion") is not None]
    wrong_rate = (round(sum(1 for row in wrong_eligible if row.get("wrong_confident_promotion"))
                        / len(wrong_eligible), 6) if wrong_eligible else None)
    evaluated = sum(1 for row in rows if row.get("coverage", 0) in (1, "1", True, "True"))
    return {
        "case_count": n,
        "top5_recall": top_k_rates["top5_recall"],
        "top10_recall": top_k_rates["top10_recall"],
        "top20_recall": top_k_rates["top20_recall"],
        "median_interval_candidate_count": float(pd.Series(intervals).median()) if n else 0.0,
        "median_narrowed_candidate_count": float(pd.Series(narrowed).median()) if n else 0.0,
        "median_reduction_ratio": float(pd.Series(reduction).median()) if n else 0.0,
        "family_recall": family_recall,
        "wrong_confident_promotion_rate": wrong_rate,
        "not_evaluated_count": sum(1 for row in rows if row.get("not_evaluated_count", 0)),
        "coverage": round(evaluated / n, 6) if n else 0.0,
        "metric_families": {name: label for name, label in BENCHMARK_METRIC_FAMILIES},
        "preferred_copy_policy": (
            "preferred_copy is an optional diagnostic; an undecided preferred copy is never "
            "an automatic benchmark miss (issue #404)"
        ),
    }


# --------------------------------------------------------------------------- #
# Issue #444: 4 QTL Use-case strata metrics and Top-20 retention contract
# --------------------------------------------------------------------------- #

from .schemas import (
    RETENTION_SCOPE_USE_CASE_TYPES,
    USE_CASE_TYPE_DESCRIPTIONS,
    USE_CASE_TYPES,
)


def _safe_rate(hits: int, total: int) -> float | None:
    return round(hits / total, 6) if total > 0 else None


def _is_truth_hit(row: dict[str, object], key: str) -> bool:
    val = row.get(key)
    if isinstance(val, bool):
        return val
    if isinstance(val, (int, float)):
        return bool(val > 0)
    text = str(val or "").strip().lower()
    return text in {"true", "1", "yes", "hit", "supported"}


def aggregate_use_case_strata_metrics(
    rows: list[dict[str, object]],
    *,
    use_case_key: str = "use_case_type",
) -> dict[str, object]:
    """Aggregate benchmark metrics stratified by the 4 QTL use-case types (issue #444).

    Primary contract:
    - Type 1–3 evaluable cases form the Primary Top-20 retention contract universe.
      Successful retention means the causal gene OR truth-containing family/copy-set
      is retained in the final <=20 Wet candidate set.
    - Type 4 (species_specific_novel) is NEVER mixed into the Type 1–3 recall denominator;
      it is evaluated separately for honest abstention/escalation (no wrong confident promotion).
    - not_evaluated and contaminated cases are never converted to misses and are
      excluded from the evaluable denominator.
    - non-gene truth units (regulatory_region, sv, pav, multi_gene, locus) and gene/family
      truths are segregated and tracked.
    """
    valid_rows = [row for row in rows if isinstance(row, dict)]

    # Normalized records
    normalized_records: list[dict[str, object]] = []
    for r in valid_rows:
        raw_uct = str(r.get(use_case_key, "") or "").strip().lower() or "unclassified"
        uct = raw_uct if raw_uct in USE_CASE_TYPES else "unclassified"
        normalized_records.append({**r, "_use_case_type": uct})

    per_type_metrics: dict[str, dict[str, object]] = {}

    for uct in sorted(USE_CASE_TYPES):
        type_rows = [r for r in normalized_records if r["_use_case_type"] == uct]
        n_total = len(type_rows)
        if n_total == 0:
            per_type_metrics[uct] = {
                "use_case_type": uct,
                "case_count": 0,
                "evaluable_count": 0,
                "not_evaluated_count": 0,
                "contaminated_count": 0,
                "top5_retention": None,
                "top10_retention": None,
                "top20_retention": None,
                "phenotype_to_family_recall": None,
                "family_to_member_retention": None,
                "median_reduction_ratio": None,
                "wrong_confident_promotion_count": 0,
                "wrong_confident_promotion_rate": None,
                "correct_abstention_or_escalation_rate": None,
                "retention_scope_eligible": uct in RETENTION_SCOPE_USE_CASE_TYPES,
                "interpretation": USE_CASE_TYPE_DESCRIPTIONS.get(uct, ""),
            }
            continue

        # Status counts
        contaminated_rows = [r for r in type_rows if str(r.get("cleanliness_state", r.get("evaluation_status", ""))).lower() == "contaminated"]
        not_eval_rows = [
            r for r in type_rows
            if str(r.get("cleanliness_state", r.get("evaluation_status", ""))).lower() in {
                "not_evaluated", "not_evaluated_no_scoring_evidence",
                "not_evaluated_evidence_provider_failed", "not_evaluated_no_supporting_evidence",
            } or r.get("not_evaluated", False) is True
        ]

        # Evaluable rows (clean / evaluated, excluding contaminated and not_evaluated)
        evaluable_rows = [r for r in type_rows if r not in contaminated_rows and r not in not_eval_rows]
        n_eval = len(evaluable_rows)

        # Top-k retention on evaluable cases
        # Retention accepts causal gene topk OR family_retention / copy_retention
        top5_hits = sum(1 for r in evaluable_rows if _is_truth_hit(r, "top5_recall") or _is_truth_hit(r, "top5_retention") or _is_truth_hit(r, "top5_conservative"))
        top10_hits = sum(1 for r in evaluable_rows if _is_truth_hit(r, "top10_recall") or _is_truth_hit(r, "top10_retention") or _is_truth_hit(r, "top10_conservative"))
        top20_hits = sum(1 for r in evaluable_rows if _is_truth_hit(r, "top20_recall") or _is_truth_hit(r, "top20_retention") or _is_truth_hit(r, "top20_conservative"))

        # Family recall and member retention
        fam_eligible = [r for r in evaluable_rows if r.get("family_recall") is not None or r.get("phenotype_to_family_recall") is not None]
        fam_hits = sum(1 for r in fam_eligible if _is_truth_hit(r, "family_recall") or _is_truth_hit(r, "phenotype_to_family_recall"))

        member_eligible = [r for r in evaluable_rows if r.get("family_to_member_retention") is not None or r.get("family_to_wet_candidate") is not None]
        member_hits = sum(1 for r in member_eligible if _is_truth_hit(r, "family_to_member_retention") or _is_truth_hit(r, "family_to_wet_candidate"))

        # Reduction ratio
        reductions = [float(r["reduction_ratio"]) for r in evaluable_rows if r.get("reduction_ratio") is not None]
        median_red = float(pd.Series(reductions).median()) if reductions else None

        # Wrong confident promotion
        wrong_conf = sum(1 for r in evaluable_rows if _is_truth_hit(r, "wrong_confident_promotion") or _is_truth_hit(r, "confident_wrong_sole_top") or _is_truth_hit(r, "false_confident_sole_top"))

        # Abstention / escalation (for Type 4 or non-gene)
        abstentions = sum(1 for r in evaluable_rows if _is_truth_hit(r, "abstained") or _is_truth_hit(r, "safe_hold") or _is_truth_hit(r, "escalated_to_mapping_or_wet"))

        per_type_metrics[uct] = {
            "use_case_type": uct,
            "case_count": n_total,
            "evaluable_count": n_eval,
            "not_evaluated_count": len(not_eval_rows),
            "contaminated_count": len(contaminated_rows),
            "top5_retention": _safe_rate(top5_hits, n_eval),
            "top10_retention": _safe_rate(top10_hits, n_eval),
            "top20_retention": _safe_rate(top20_hits, n_eval),
            "phenotype_to_family_recall": _safe_rate(fam_hits, len(fam_eligible)),
            "family_to_member_retention": _safe_rate(member_hits, len(member_eligible)),
            "median_reduction_ratio": median_red,
            "wrong_confident_promotion_count": wrong_conf,
            "wrong_confident_promotion_rate": _safe_rate(wrong_conf, n_eval),
            "correct_abstention_or_escalation_rate": _safe_rate(abstentions, n_eval) if uct == "type4_species_specific_novel" else None,
            "retention_scope_eligible": uct in RETENTION_SCOPE_USE_CASE_TYPES,
            "interpretation": USE_CASE_TYPE_DESCRIPTIONS.get(uct, ""),
        }

    # Primary retention contract summary (Type 1–3 evaluable cases ONLY)
    retention_rows = [r for r in normalized_records if r["_use_case_type"] in RETENTION_SCOPE_USE_CASE_TYPES]
    ret_contaminated = [r for r in retention_rows if str(r.get("cleanliness_state", r.get("evaluation_status", ""))).lower() == "contaminated"]
    ret_not_eval = [
        r for r in retention_rows
        if str(r.get("cleanliness_state", r.get("evaluation_status", ""))).lower() in {
            "not_evaluated", "not_evaluated_no_scoring_evidence",
            "not_evaluated_evidence_provider_failed", "not_evaluated_no_supporting_evidence",
        } or r.get("not_evaluated", False) is True
    ]
    ret_evaluable = [r for r in retention_rows if r not in ret_contaminated and r not in ret_not_eval]

    n_ret_total = len(retention_rows)
    n_ret_eval = len(ret_evaluable)

    ret_top5 = sum(1 for r in ret_evaluable if _is_truth_hit(r, "top5_recall") or _is_truth_hit(r, "top5_retention") or _is_truth_hit(r, "top5_conservative"))
    ret_top10 = sum(1 for r in ret_evaluable if _is_truth_hit(r, "top10_recall") or _is_truth_hit(r, "top10_retention") or _is_truth_hit(r, "top10_conservative"))
    ret_top20 = sum(1 for r in ret_evaluable if _is_truth_hit(r, "top20_recall") or _is_truth_hit(r, "top20_retention") or _is_truth_hit(r, "top20_conservative"))
    ret_reductions = [float(r["reduction_ratio"]) for r in ret_evaluable if r.get("reduction_ratio") is not None]
    ret_wrong_conf = sum(1 for r in ret_evaluable if _is_truth_hit(r, "wrong_confident_promotion") or _is_truth_hit(r, "confident_wrong_sole_top"))

    # Truth unit breakdown
    truth_units: dict[str, int] = {}
    for r in valid_rows:
        unit = str(r.get("truth_unit_type", r.get("truth_unit", "gene")) or "gene").lower()
        truth_units[unit] = truth_units.get(unit, 0) + 1

    return {
        "primary_retention_contract": {
            "scope": "type1_to_3_evaluable_cases",
            "total_cases": n_ret_total,
            "evaluable_cases": n_ret_eval,
            "not_evaluated_cases": len(ret_not_eval),
            "contaminated_cases": len(ret_contaminated),
            "top5_retention_rate": _safe_rate(ret_top5, n_ret_eval),
            "top10_retention_rate": _safe_rate(ret_top10, n_ret_eval),
            "top20_retention_rate": _safe_rate(ret_top20, n_ret_eval),
            "median_reduction_ratio": float(pd.Series(ret_reductions).median()) if ret_reductions else None,
            "wrong_confident_promotion_count": ret_wrong_conf,
            "wrong_confident_promotion_rate": _safe_rate(ret_wrong_conf, n_ret_eval),
            "target_top20_retention_min": 0.90,
            "target_median_reduction_min": 0.80,
            "target_wrong_confident_max": 0.05,
            "type4_strictly_excluded_from_denominator": True,
        },
        "use_case_strata": per_type_metrics,
        "truth_unit_distribution": truth_units,
        "total_benchmark_cases": len(valid_rows),
    }
