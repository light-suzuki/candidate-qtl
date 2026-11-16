from __future__ import annotations

import hashlib
import json
import random

import pandas as pd

from .core import runtime_path
from .evidence_io import CONFIDENCE_MULTIPLIERS
from .hardening import independence_group
from .ranking_policy import INDEPENDENCE_FAMILIES, load_weight_registry


BENCHMARK_VERSION = "1.0.0"
WET_BUDGETS = (10, 25, 50)


def _truth_ids(value: object) -> list[str]:
    return [part.strip() for part in str(value).split("|") if part.strip()]


def _numeric_profiles(registry: pd.DataFrame) -> dict[str, dict[str, float]]:
    weighted = registry[registry["base_weight_or_tier"].str.strip().ne("")].copy()
    current = {
        str(row.evidence_code): float(row.base_weight_or_tier)
        for row in weighted.itertuples()
    }
    uniform = {code: 1.0 for code in current}

    def tier(value: float) -> float:
        if value >= 3.0:
            return 3.0
        if value >= 1.0:
            return 2.0
        if value > 0.0:
            return 1.0
        return 0.0

    return {
        "current_registry": current,
        "uniform_family": uniform,
        "tier_only": {code: tier(value) for code, value in current.items()},
    }


def _eligibility(evidence: pd.DataFrame, gene_ids: list[str]) -> dict[str, str]:
    support = evidence[evidence["direction"].eq("support")].copy()
    explicit = support.get("independence_group", pd.Series("", index=support.index))
    support["_source_group"] = [
        independence_group(source, group)
        for source, group in zip(support["source_id"], explicit)
    ]
    codes = support.groupby("gene_id")["evidence_code"].apply(set)
    sources = support[support["_source_group"].ne("")].groupby("gene_id")["_source_group"].apply(set)
    citable_codes = support[support["_source_group"].ne("")].groupby("gene_id")["evidence_code"].apply(set)
    conflict = set()
    for (gene_id, _), rows in evidence.groupby(["gene_id", "evidence_code"]):
        if {"support", "against"}.issubset(set(rows["direction"])):
            conflict.add(str(gene_id))
    strong = {"direct_model_genetic", "ortholog_model_genetic", "known_qtg"}
    result: dict[str, str] = {}
    for gene_id in gene_ids:
        gene_codes = codes.get(gene_id, set())
        families = {INDEPENDENCE_FAMILIES.get(code, code) for code in gene_codes}
        result[gene_id] = (
            "conflicting_evidence" if gene_id in conflict
            else "strong_genetic_evidence" if citable_codes.get(gene_id, set()).intersection(strong)
            else "independent_experimental_support" if len(sources.get(gene_id, set())) >= 2
            else "multiple_non_redundant_channels" if len(families) >= 2
            else "insufficient_independent_evidence"
        )
    return result


def _scores(
    evidence: pd.DataFrame,
    gene_ids: list[str],
    weights: dict[str, float],
    *,
    family_cap: bool,
    confidence_multiplier: bool,
) -> tuple[dict[str, float], set[str]]:
    rows: list[dict[str, object]] = []
    observed_codes: set[str] = set()
    for (gene_id, code, direction), group in evidence.groupby(
        ["gene_id", "evidence_code", "direction"], sort=False
    ):
        if code not in weights or direction not in {"support", "against"}:
            continue
        observed_codes.add(str(code))
        sign = 1.0 if direction == "support" else -1.0
        candidates = []
        for item in group.itertuples():
            multiplier = CONFIDENCE_MULTIPLIERS.get(str(item.confidence), 0.0) if confidence_multiplier else 1.0
            candidates.append(sign * weights[str(code)] * multiplier)
        contribution = max(candidates, key=abs)
        rows.append({
            "gene_id": str(gene_id),
            "family": INDEPENDENCE_FAMILIES.get(str(code), str(code)),
            "direction": str(direction),
            "contribution": contribution,
        })
    contributions = pd.DataFrame(rows)
    totals = {gene_id: 0.0 for gene_id in gene_ids}
    if contributions.empty:
        return totals, observed_codes
    if family_cap:
        used = []
        for _, group in contributions.groupby(["gene_id", "family", "direction"], sort=False):
            used.append(group.loc[group["contribution"].abs().idxmax()])
        contributions = pd.DataFrame(used)
    summed = contributions.groupby("gene_id")["contribution"].sum()
    totals.update({str(gene_id): float(value) for gene_id, value in summed.items()})
    return totals, observed_codes


def _evaluate_scores(
    scores: dict[str, float], truth_ids: list[str], candidate_count: int
) -> dict[str, object]:
    present = [gene_id for gene_id in truth_ids if gene_id in scores]
    if not present:
        return {
            "covered": False, "rank_min": None, "rank_max": None, "tie_size": None,
            "mrr_conservative": None, "top1": False, "top5": False,
            "top10": False, "top20": False, "assays_to_truth": None,
            "compression_at_20": None,
        }
    best_score = max(scores[gene_id] for gene_id in present)
    greater = sum(value > best_score for value in scores.values())
    equal = sum(value == best_score for value in scores.values())
    rank_min = greater + 1
    rank_max = greater + equal
    return {
        "covered": True,
        "rank_min": rank_min,
        "rank_max": rank_max,
        "tie_size": equal,
        "mrr_conservative": 1.0 / rank_max,
        "top1": rank_max <= 1,
        "top5": rank_max <= 5,
        "top10": rank_max <= 10,
        "top20": rank_max <= 20,
        "assays_to_truth": rank_max,
        "compression_at_20": candidate_count / min(20, candidate_count),
    }


def _bootstrap(values: list[float], iterations: int, seed: int) -> tuple[float | None, float | None]:
    if not values:
        return None, None
    rng = random.Random(seed)
    estimates = []
    for _ in range(iterations):
        sample = [values[rng.randrange(len(values))] for _ in values]
        estimates.append(sum(sample) / len(sample))
    estimates.sort()
    low = estimates[int(0.025 * (len(estimates) - 1))]
    high = estimates[int(0.975 * (len(estimates) - 1))]
    return low, high


def weight_policy_benchmark(
    manifest_arg: str,
    out_dir_arg: str,
    blocked_terms: list[str] | None = None,
    bootstrap_iterations: int = 1000,
    seed: int = 20260718,
    crop: str = "Pisum sativum",
) -> dict[str, object]:
    """Compare numeric weight profiles and structural gates on frozen cases.

    This is an ordinal assay-prioritisation benchmark.  It does not fit weights
    and does not estimate causal probabilities.
    """
    if bootstrap_iterations < 100:
        raise ValueError("bootstrap_iterations must be at least 100")
    manifest_path = runtime_path(manifest_arg)
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    source_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    registry, registry_meta = load_weight_registry()
    profiles = _numeric_profiles(registry)
    summary_path = manifest_path.parent / "benchmark_summary.tsv"
    summary = pd.read_csv(summary_path, sep="\t", dtype=str).fillna("")
    truth_by_case = {
        str(row.case_id): _truth_ids(row.answer_gene_id)
        for row in summary.itertuples()
    }
    rows: list[dict[str, object]] = []
    observed_by_code: dict[str, set[str]] = {}
    structural_profiles = {
        "current_registry": (profiles["current_registry"], True, True),
        "uniform_family": (profiles["uniform_family"], True, True),
        "tier_only": (profiles["tier_only"], True, True),
        "current_no_family_cap": (profiles["current_registry"], False, True),
        "current_no_confidence_multiplier": (profiles["current_registry"], True, False),
    }
    for case in source_manifest.get("cases", []):
        case_id = str(case["case_id"])
        ranked_path = runtime_path(case["rank"]["ranked_candidates"])
        evidence_path = runtime_path(case["evidence"]["evidence"])
        ranked = pd.read_csv(ranked_path, sep="\t", dtype=str).fillna("")
        evidence = pd.read_csv(evidence_path, sep="\t", dtype=str).fillna("")
        gene_ids = [str(value) for value in ranked["gene_id"]]
        truths = truth_by_case.get(case_id, [])
        common = {
            "case_id": case_id,
            "locus_id": str(case.get("locus_id", "")),
            "crop": str(case.get("crop", crop)),
            "trait_profile": str(case.get("term_profile", "")),
            "truth_tier": str(case.get("truth_tier", "")),
            "gene_class": str(case.get("gene_class", "")),
            "candidate_count": len(gene_ids),
            "truth_gene_ids": "|".join(truths),
        }
        for profile_name, (weights, family_cap, confidence_multiplier) in structural_profiles.items():
            scores, observed = _scores(
                evidence, gene_ids, weights,
                family_cap=family_cap,
                confidence_multiplier=confidence_multiplier,
            )
            for code in observed:
                observed_by_code.setdefault(code, set()).add(case_id)
            evaluation = _evaluate_scores(scores, truths, len(gene_ids))
            rows.append({
                **common,
                "policy_profile": profile_name,
                "family_cap": family_cap,
                "confidence_multiplier": confidence_multiplier,
                "eligibility_only": False,
                **evaluation,
            })
        eligibility = _eligibility(evidence, gene_ids)
        eligibility_score = {
            "strong_genetic_evidence": 4.0,
            "independent_experimental_support": 3.0,
            "multiple_non_redundant_channels": 2.0,
            "insufficient_independent_evidence": 1.0,
            "conflicting_evidence": 0.0,
        }
        evaluation = _evaluate_scores(
            {gene_id: eligibility_score[eligibility[gene_id]] for gene_id in gene_ids},
            truths,
            len(gene_ids),
        )
        rows.append({
            **common,
            "policy_profile": "eligibility_only",
            "family_cap": True,
            "confidence_multiplier": False,
            "eligibility_only": True,
            **evaluation,
        })

    detail = pd.DataFrame(rows)
    detail_path = out_dir / "weight_policy_benchmark.tsv"
    detail.to_csv(detail_path, sep="\t", index=False)
    comparison_rows = []
    uncertainty_rows = []
    wet_rows = []
    for profile, group in detail.groupby("policy_profile", sort=False):
        covered = group[group["covered"] == True]  # noqa: E712
        metrics = {
            "policy_profile": profile,
            "case_count": len(group),
            "covered_count": len(covered),
            "uncovered_count": len(group) - len(covered),
            "top1_recall_covered": float(covered["top1"].mean()) if len(covered) else None,
            "top5_recall_covered": float(covered["top5"].mean()) if len(covered) else None,
            "top10_recall_covered": float(covered["top10"].mean()) if len(covered) else None,
            "top20_recall_covered": float(covered["top20"].mean()) if len(covered) else None,
            "mrr_conservative": float(covered["mrr_conservative"].mean()) if len(covered) else None,
            "median_rank_max": float(covered["rank_max"].median()) if len(covered) else None,
            "median_compression_at_20": float(covered["compression_at_20"].median()) if len(covered) else None,
        }
        comparison_rows.append(metrics)
        for metric in ("top5", "top20", "mrr_conservative"):
            values = [float(value) for value in covered[metric].tolist()]
            low, high = _bootstrap(values, bootstrap_iterations, seed + len(uncertainty_rows))
            uncertainty_rows.append({
                "policy_profile": profile, "metric": metric,
                "estimate": sum(values) / len(values) if values else None,
                "ci95_low": low, "ci95_high": high,
                "bootstrap_iterations": bootstrap_iterations, "seed": seed,
            })
        for budget in WET_BUDGETS:
            hits = covered["rank_max"].astype(float).le(budget)
            wet_rows.append({
                "policy_profile": profile, "wet_budget": budget,
                "covered_count": len(covered),
                "recall_covered": float(hits.mean()) if len(covered) else None,
            })
    comparison = pd.DataFrame(comparison_rows)
    comparison_path = out_dir / "weight_profile_comparison.tsv"
    comparison.to_csv(comparison_path, sep="\t", index=False)
    uncertainty_path = out_dir / "weight_uncertainty.tsv"
    pd.DataFrame(uncertainty_rows).to_csv(uncertainty_path, sep="\t", index=False)
    wet_path = out_dir / "wet_budget_recall.tsv"
    pd.DataFrame(wet_rows).to_csv(wet_path, sep="\t", index=False)
    ablation = comparison[comparison["policy_profile"].isin({
        "current_registry", "current_no_family_cap", "current_no_confidence_multiplier", "eligibility_only"
    })].copy()
    ablation_path = out_dir / "policy_ablation.tsv"
    ablation.to_csv(ablation_path, sep="\t", index=False)

    validity_rows = []
    for item in registry.itertuples():
        count = len(observed_by_code.get(str(item.evidence_code), set()))
        if count == 0:
            status = "insufficient_benchmark_evidence"
        elif str(item.base_weight_or_tier).strip() == "":
            status = "supported_as_is"
        else:
            status = "direction_supported_magnitude_uncertain"
        validity_rows.append({
            "evidence_code": item.evidence_code,
            "benchmark_case_count": count,
            "validity_status": status,
            "rationale": (
                "not observed in this benchmark" if count == 0
                else "context/rescue remained outside numeric scoring" if status == "supported_as_is"
                else "observed, but this single-crop calibration cannot identify a precise universal magnitude"
            ),
        })
    validity_path = out_dir / "evidence_code_validity.tsv"
    pd.DataFrame(validity_rows).to_csv(validity_path, sep="\t", index=False)

    crop_groups = sorted(set(detail["crop"]))
    holdout_rows = []
    for (profile, crop_name), group in detail.groupby(["policy_profile", "crop"], sort=False):
        covered = group[group["covered"] == True]  # noqa: E712
        holdout_rows.append({
            "policy_profile": profile,
            "held_out_crop": crop_name,
            "case_count": len(group),
            "top20_recall_covered": float(covered["top20"].mean()) if len(covered) else None,
            "cross_crop_holdout_estimable": len(crop_groups) > 1,
        })
    holdout_path = out_dir / "group_holdout.tsv"
    pd.DataFrame(holdout_rows).to_csv(holdout_path, sep="\t", index=False)

    current = comparison.set_index("policy_profile").loc["current_registry"]
    # Issue #367: when the source benchmark-reference manifest carries the MRR
    # bootstrap CI (mrr_conservative_ci_95_low/high in its resolution_metrics),
    # print it next to this module's own conservative MRR read-out.
    source_resolution = source_manifest.get("resolution_metrics", {}) or {}
    source_mrr_low = source_resolution.get("mrr_conservative_ci_95_low")
    source_mrr_high = source_resolution.get("mrr_conservative_ci_95_high")
    mrr_ci_suffix = (
        f" (95% CI {float(source_mrr_low):.3f}-{float(source_mrr_high):.3f})"
        if source_mrr_low is not None and source_mrr_high is not None else ""
    )
    report_lines = [
        "# Native weight policy validity benchmark",
        "",
        "This is an ordinal pre-Wet assay-prioritisation benchmark, not a causal probability model or statistical calibration.",
        "",
        "## Result",
        "",
        f"- cases: {int(current.case_count)}; uncovered: {int(current.uncovered_count)}",
        f"- current strict tie-aware Recall@20 (covered): {float(current.top20_recall_covered):.3f}",
        f"- current conservative MRR: {float(current.mrr_conservative):.3f}{mrr_ci_suffix}",
        "- no weight change is recommended from a single-crop, internally curated calibration set.",
        "",
        "## Interpretation",
        "",
        "Family caps, confidence multipliers, and eligibility-only ordering are reported as structural ablations; uniform and tier-only are numeric baselines. Boundary ties count as failures when the full tie exceeds the budget. Uncovered truth genes are excluded from covered-case recall and reported separately.",
        "",
        "## Remaining external validation",
        "",
        "Cross-crop hold-out is not estimable from this pea-only set. A prospective or independently curated crop/reference panel is required before changing the default registry.",
    ]
    report_path = out_dir / "weight_validity_report.md"
    report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    catalog_blocked: set[str] = set(blocked_terms or [])
    for column in ("gene_name", "locus_id", "answer_gene_id"):
        if column in summary.columns:
            for value in summary[column]:
                catalog_blocked.update(_truth_ids(value))
    blocked = sorted(term for term in catalog_blocked if term)
    leaked_terms: list[dict[str, str]] = []
    blocked_folded = {term.casefold(): term for term in blocked}
    for case in source_manifest.get("cases", []):
        for term in case.get("trait_terms", []):
            if str(term).strip().casefold() in blocked_folded:
                leaked_terms.append({
                    "case_id": str(case.get("case_id", "")),
                    "term": str(term),
                    "blocked_as": blocked_folded[str(term).strip().casefold()],
                })
    output_manifest = {
        "schema_version": "1.0",
        "benchmark_version": BENCHMARK_VERSION,
        "source_benchmark_manifest": str(manifest_path),
        "source_benchmark_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        "weight_registry": registry_meta,
        "blocked_terms": blocked,
        "blocked_terms_sha256": hashlib.sha256("\n".join(blocked).encode()).hexdigest(),
        "blocked_term_audit": {
            "policy": "exact case-insensitive match against gene, locus, and answer identifiers",
            "leak_count": len(leaked_terms),
            "leaks": leaked_terms,
        },
        "policy_profiles": list(comparison["policy_profile"]),
        "tie_policy": "strict: success only when rank_max <= budget",
        "coverage_policy": "uncovered is reported separately and never assigned the lowest rank",
        "probability_claim": False,
        "llm_used": False,
        "crop_groups": crop_groups,
        "cross_crop_holdout_estimable": len(crop_groups) > 1,
        "bootstrap": {"iterations": bootstrap_iterations, "seed": seed},
        "recommended_weight_registry_emitted": False,
        "recommendation": "retain current registry pending independent cross-crop validation",
    }
    output_manifest_path = out_dir / "weight_policy_benchmark_manifest.json"
    output_manifest_path.write_text(json.dumps(output_manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return {
        "detail": str(detail_path),
        "comparison": str(comparison_path),
        "ablation": str(ablation_path),
        "uncertainty": str(uncertainty_path),
        "wet_budget": str(wet_path),
        "validity": str(validity_path),
        "holdout": str(holdout_path),
        "report": str(report_path),
        "manifest": str(output_manifest_path),
        "case_count": int(current.case_count),
        "current_top20_recall": float(current.top20_recall_covered),
    }
