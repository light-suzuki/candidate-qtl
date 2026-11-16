from __future__ import annotations

import hashlib
import json
import tempfile
from pathlib import Path

import pandas as pd

from .core import runtime_path
from .hardening import prioritize_evidence, support_group
from .ranking_policy import INDEPENDENCE_FAMILIES, load_weight_registry, load_weights


# Rank positions at or beyond this absolute shift count as sensitivity even when
# the candidate never crosses the top_k boundary. Diagnostic threshold only;
# it never feeds back into scoring or the default rank.
RANK_SHIFT_THRESHOLD = 3


def _truth_ids(value: object) -> list[str]:
    return [part.strip() for part in str(value).split("|") if part.strip()]


def _rank_values(frame: pd.DataFrame) -> pd.Series:
    return pd.to_numeric(frame.set_index("gene_id")["screening_rank"], errors="raise").astype(int)


def _rank_max_values(frame: pd.DataFrame) -> pd.Series:
    return pd.to_numeric(frame.set_index("gene_id")["screening_rank_max"], errors="raise").astype(int)


def _kendall(a: pd.Series, b: pd.Series) -> float:
    pairs = list(zip(a.tolist(), b.tolist()))
    concordant = discordant = 0
    for i in range(len(pairs)):
        for j in range(i + 1, len(pairs)):
            product = (pairs[i][0] - pairs[j][0]) * (pairs[i][1] - pairs[j][1])
            if product > 0:
                concordant += 1
            elif product < 0:
                discordant += 1
    total = concordant + discordant
    return (concordant - discordant) / total if total else 1.0


def _spearman(a: pd.Series, b: pd.Series) -> float:
    ranked_a, ranked_b = a.rank(), b.rank()
    if ranked_a.nunique() == 1 or ranked_b.nunique() == 1:
        return 1.0 if ranked_a.equals(ranked_b) else 0.0
    return float(ranked_a.corr(ranked_b))


def _evidence_group_series(evidence: pd.DataFrame) -> pd.Series:
    """Resolve the scoring-layer independence group of every evidence row.

    Uses ``support_group`` -- the same resolver ``prioritize_evidence`` uses --
    so a source-dropout scenario removes exactly the rows the scorer would treat
    as one independent group (explicit group, citable source, or the shared
    comparative chain). Rows that resolve to "" carry no independence claim and
    are never a dropout target.
    """
    explicit = evidence["independence_group"] if "independence_group" in evidence.columns else pd.Series("", index=evidence.index)
    return pd.Series(
        [support_group(source, group, code) for source, group, code in zip(evidence["source_id"], explicit, evidence["evidence_code"])],
        index=evidence.index, name="independence_group_resolved",
    )


def _run_stability_scenarios(candidates_path: Path, evidence: pd.DataFrame, weights: dict[str, float], top_k: int,
                             perturbations: list[float], work: Path,
                             source_dropout_groups: list[str] | None = None,
                             provider_dropout_types: list[str] | None = None) -> dict[str, object]:
    """Run the deterministic default plus every adversarial scenario once.

    Scenarios (names are also the ranking_sensitivity.tsv scenario labels):

    * ``default_weighted`` -- the canonical registry ranking.
    * ``eligibility_first`` -- same score, ordered by eligibility group first.
    * ``leave_out:<family>`` -- rows of one INDEPENDENCE_FAMILIES family removed.
    * ``perturb:<code>:<minus|plus><pct>`` -- one evidence code's weight moved by
      a fraction; never tuned to truth, always centered on the default registry.
    * ``confidence:flat|conservative`` -- confidence re-graded.
    * ``provider:<source_type>`` -- keep only one provider (existing behavior).
    * ``source_dropout:<group>`` -- all rows whose resolved independence group
      equals ``group`` removed (groups collapsed via ``support_group``).
    * ``provider_dropout:<source_type>`` -- all rows of that source_type removed.

    The default rank is never overwritten by any scenario; scenario ranks are
    diagnostics only.
    """
    perturbations = perturbations or [0.10, 0.25, 0.50]
    scenarios: dict[str, pd.Series] = {}
    scenarios_max: dict[str, pd.Series] = {}
    metadata: list[dict[str, str]] = []

    def run(name: str, category: str, chain: str, evidence_frame: pd.DataFrame, scenario_weights: dict[str, float]) -> pd.DataFrame:
        safe = name.replace("/", "_").replace(":", "_")
        ev = work / f"{safe}.evidence.tsv"; wt = work / f"{safe}.weights.tsv"; rank_dir = work / safe
        evidence_frame.to_csv(ev, sep="\t", index=False)
        pd.DataFrame(sorted(scenario_weights.items()), columns=["evidence_code", "weight"]).to_csv(wt, sep="\t", index=False)
        prioritize_evidence(str(candidates_path), str(ev), str(rank_dir), str(wt), top_k, top_k, 0.0, top_k)
        ranked = pd.read_csv(rank_dir / "ranked_candidates.tsv", sep="\t", dtype=str).fillna("")
        scenarios[name] = _rank_values(ranked)
        scenarios_max[name] = _rank_max_values(ranked)
        metadata.append({"scenario": name, "category": category, "analysis_chain": chain})
        return ranked

    default_ranked = run("default_weighted", "default", "native", evidence, weights)
    hierarchy = default_ranked.copy()
    group_order = {"top": 0, "reserve": 1, "rescue": 2, "evidence_insufficient": 3}
    hierarchy["_group"] = hierarchy["dry_group"].map(group_order).fillna(9)
    hierarchy = hierarchy.sort_values(["_group", "screening_rank", "gene_id"], kind="stable")
    hierarchy["screening_rank"] = range(1, len(hierarchy) + 1)
    scenarios["eligibility_first"] = _rank_values(hierarchy)
    scenarios_max["eligibility_first"] = _rank_values(hierarchy)
    metadata.append({"scenario": "eligibility_first", "category": "eligibility", "analysis_chain": "native"})

    families = sorted({INDEPENDENCE_FAMILIES.get(code, code) for code in evidence["evidence_code"].unique()})
    for family in families:
        kept = evidence[evidence["evidence_code"].map(lambda code: INDEPENDENCE_FAMILIES.get(code, code) != family)]
        run(f"leave_out:{family}", "family_dropout", f"family:{family}", kept, weights)

    for code in sorted(set(evidence["evidence_code"]).intersection(weights)):
        for fraction in perturbations:
            for sign, label in ((-1, "minus"), (1, "plus")):
                changed = dict(weights); changed[code] = weights[code] * (1 + sign * fraction)
                run(f"perturb:{code}:{label}{int(fraction * 100)}", "weight_perturbation", f"weight:{code}", evidence, changed)

    high_confidence = evidence.copy(); high_confidence.loc[high_confidence["direction"].isin({"support", "against"}), "confidence"] = "high"
    run("confidence:flat", "confidence", "confidence", high_confidence, weights)
    conservative = evidence.copy(); conservative.loc[conservative["confidence"].eq("medium"), "confidence"] = "low"
    run("confidence:conservative", "confidence", "confidence", conservative, weights)

    if "source_type" in evidence.columns:
        for provider in sorted(value for value in evidence["source_type"].unique() if value):
            run(f"provider:{provider}", "provider", f"provider:{provider}", evidence[evidence["source_type"].eq(provider)], weights)

    resolved = _evidence_group_series(evidence)
    groups = sorted(set(resolved[resolved.ne("")]))
    if source_dropout_groups is not None:
        groups = sorted(set(groups) & set(source_dropout_groups))
    for group in groups:
        run(f"source_dropout:{group}", "source_dropout", f"source_group:{group}", evidence[resolved.ne(group)], weights)

    if "source_type" in evidence.columns:
        types = sorted({value for value in evidence["source_type"].unique() if value})
        if provider_dropout_types is not None:
            types = sorted(set(types) & set(provider_dropout_types))
        for stype in types:
            run(f"provider_dropout:{stype}", "provider_dropout", f"source_type:{stype}", evidence[evidence["source_type"].ne(stype)], weights)

    default = scenarios["default_weighted"]
    dropout_rows: list[dict[str, object]] = []
    for name in scenarios:
        entry = next(item for item in metadata if item["scenario"] == name)
        if entry["category"] not in {"source_dropout", "provider_dropout"}:
            continue
        dropped = entry["analysis_chain"].split(":", 1)[1]
        for gene_id in scenarios[name].index:
            default_rank = int(default[gene_id]); scenario_rank = int(scenarios[name][gene_id])
            dropout_rows.append({"gene_id": gene_id, "scenario": name, "scenario_type": entry["category"],
                                 "dropped_label": dropped, "default_rank": default_rank, "scenario_rank": scenario_rank,
                                 "rank_shift": scenario_rank - default_rank,
                                 "crossed_top_k": (default_rank <= top_k) != (scenario_rank <= top_k)})
    return {"scenarios": scenarios, "scenarios_max": scenarios_max, "metadata": metadata,
            "default_ranked": default_ranked, "dropout_rows": dropout_rows}


def _stability_class(gene_id: str, default: pd.Series, rank_matrix: pd.DataFrame, metadata_frame: pd.DataFrame,
                     tie_sizes: pd.Series, top_k: int) -> str:
    """Classify one candidate's rank stability (issue #228).

    Priority order: an underdetermined tie at the default rank is the most
    fundamental uncertainty, then source-group sensitivity, then provider
    sensitivity, then weight/structure sensitivity, else stable. A category is
    sensitive when one of its scenarios crosses the top_k boundary; for
    candidates inside top_k by default, an absolute shift of at least
    ``RANK_SHIFT_THRESHOLD`` also counts (deep-reserve candidates only ever
    move below the boundary, which never changes their review status).
    """
    if int(tie_sizes[gene_id]) > 1:
        return "underdetermined_ties"
    row = rank_matrix.loc[gene_id]
    default_rank = int(default[gene_id])
    sensitive: dict[str, bool] = {}
    for category in ("source_dropout", "provider_dropout", "weight"):
        columns = [column for column in rank_matrix if metadata_frame.loc[column, "category"] in
                   (("weight_perturbation", "confidence", "eligibility", "family_dropout") if category == "weight" else (category,))]
        if not columns:
            sensitive[category] = False
            continue
        crosses = any((int(row[column]) <= top_k) != (default_rank <= top_k) for column in columns)
        shift_trip = default_rank <= top_k and max(abs(int(row[column]) - default_rank) for column in columns) >= RANK_SHIFT_THRESHOLD
        sensitive[category] = crosses or shift_trip
    if sensitive["source_dropout"]:
        return "source_sensitive"
    if sensitive["provider_dropout"]:
        return "provider_sensitive"
    if sensitive["weight"]:
        return "weight_sensitive"
    return "stable"


def ranking_qa(candidates_arg: str, evidence_arg: str, out_dir_arg: str, weights_arg: str | None = None,
               top_k: int = 5, perturbations: list[float] | None = None,
               benchmark_mode: bool = False,
               source_dropout_groups: list[str] | None = None,
               provider_dropout_types: list[str] | None = None) -> dict[str, object]:
    if top_k < 1:
        raise ValueError("top_k must be at least 1")
    perturbations = perturbations or [0.10, 0.25, 0.50]
    if any(value <= 0 or value >= 1 for value in perturbations):
        raise ValueError("perturbations must be fractions between 0 and 1")
    candidates_path, evidence_path = runtime_path(candidates_arg), runtime_path(evidence_arg)
    evidence = pd.read_csv(evidence_path, sep="\t", dtype=str).fillna("")
    weights = load_weights(weights_arg)
    registry, registry_meta = load_weight_registry(weights_arg) if weights_arg and "base_weight_or_tier" in pd.read_csv(runtime_path(weights_arg), sep="\t", nrows=0).columns else load_weight_registry()
    out = runtime_path(out_dir_arg); out.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="candidate-qtl-ranking-qa-") as tmp:
        result = _run_stability_scenarios(candidates_path, evidence, weights, top_k, perturbations, Path(tmp),
                                          source_dropout_groups, provider_dropout_types)
    scenarios = result["scenarios"]
    metadata = result["metadata"]
    default_ranked = result["default_ranked"]

    rank_matrix = pd.DataFrame(scenarios).sort_index()
    default = rank_matrix["default_weighted"]
    concordance = []
    default_top = set(default.nsmallest(min(top_k, len(default))).index)
    for name in rank_matrix.columns:
        current = rank_matrix[name]
        current_top = set(current.nsmallest(min(top_k, len(current))).index)
        union = default_top | current_top
        concordance.append({"scenario": name, "spearman": round(_spearman(default, current), 6),
                            "kendall": round(_kendall(default, current), 6),
                            "top_k_jaccard": round(len(default_top & current_top) / len(union), 6) if union else 1.0})

    distributions = []
    metadata_frame = pd.DataFrame(metadata).set_index("scenario")
    dropout_columns = [column for column in rank_matrix if metadata_frame.loc[column, "category"] == "family_dropout"]
    perturb_columns = [column for column in rank_matrix if metadata_frame.loc[column, "category"] == "weight_perturbation"]
    provider_columns = [column for column in rank_matrix if metadata_frame.loc[column, "category"] == "provider"]
    source_dropout_columns = [column for column in rank_matrix if metadata_frame.loc[column, "category"] == "source_dropout"]
    provider_dropout_columns = [column for column in rank_matrix if metadata_frame.loc[column, "category"] == "provider_dropout"]
    # Collapse each correlated category before consensus so dozens of perturbations
    # cannot masquerade as independent votes.
    category_ranks = pd.DataFrame(index=rank_matrix.index)
    for category in ("default", "eligibility", "family_dropout", "confidence", "provider", "source_dropout", "provider_dropout"):
        cols = [column for column in rank_matrix if metadata_frame.loc[column, "category"] == category]
        if cols:
            category_ranks[category] = rank_matrix[cols].median(axis=1)
    consensus_rank = category_ranks.median(axis=1).rank(method="min").astype(int)
    tie_columns = {"screening_tie_size", "screening_rank_max"}
    default_tie = default_ranked.set_index("gene_id")[list(tie_columns & set(default_ranked.columns))]
    tie_sizes = default_tie["screening_tie_size"].astype(int)
    tie_maxes = default_tie["screening_rank_max"].astype(int)
    for gene_id in rank_matrix.index:
        row = rank_matrix.loc[gene_id]
        dropout_shift = max((abs(int(row[column]) - int(default[gene_id])) for column in dropout_columns), default=0)
        perturb_range = (int(row[perturb_columns].max() - row[perturb_columns].min())) if perturb_columns else 0
        top_rate = float((row <= top_k).mean())
        provider_hits = int((row[provider_columns] <= top_k).sum()) if provider_columns else 0
        source_crosses = any((int(row[column]) <= top_k) != (int(default[gene_id]) <= top_k) for column in source_dropout_columns)
        provider_dropout_crosses = any((int(row[column]) <= top_k) != (int(default[gene_id]) <= top_k) for column in provider_dropout_columns)
        scenario_shift_max = max((abs(int(row[column]) - int(default[gene_id])) for column in rank_matrix if column != "default_weighted"), default=0)
        reasons = []
        if dropout_shift > 0 and int(default[gene_id]) <= top_k and any(int(row[column]) > top_k for column in dropout_columns): reasons.append("single_family_dependent")
        if perturb_range > 0 and int(default[gene_id]) <= top_k and any(int(row[column]) > top_k for column in perturb_columns): reasons.append("weight_sensitive")
        if source_crosses: reasons.append("source_group_dependent")
        if provider_dropout_crosses: reasons.append("provider_dependent")
        if provider_columns and len(set(int(row[column]) for column in provider_columns)) > 1: reasons.append("provider_discordant")
        consensus_class = "stable_top" if int(default[gene_id]) <= top_k and top_rate >= 0.8 else ("stable_reserve" if top_rate >= 0.5 else (reasons[0] if reasons else "unstable"))
        distributions.append({"gene_id": gene_id, "default_rank": int(default[gene_id]), "consensus_rank": int(consensus_rank[gene_id]),
                              "median_rank": float(row.median()), "min_rank": int(row.min()), "max_rank": int(row.max()),
                              "top_k_rate": round(top_rate, 6), "family_dropout_max_shift": dropout_shift,
                              "weight_perturbation_rank_range": perturb_range, "provider_top_k_count": provider_hits,
                              "consensus_class": consensus_class, "instability_reason": ";".join(reasons),
                              "stability_class": _stability_class(gene_id, default, rank_matrix, metadata_frame, tie_sizes, top_k),
                              "screening_tie_size": int(tie_sizes[gene_id]), "screening_rank_max": int(tie_maxes[gene_id]),
                              "scenario_shift_max": scenario_shift_max})
    distribution = pd.DataFrame(distributions).sort_values(["consensus_rank", "default_rank", "gene_id"])
    sensitivity = rank_matrix.reset_index(names="gene_id").melt(id_vars="gene_id", var_name="scenario", value_name="rank").merge(pd.DataFrame(metadata), on="scenario", how="left")
    consensus = distribution[distribution["consensus_class"].isin({"stable_top", "stable_reserve"})]
    discordance = distribution[distribution["instability_reason"].ne("")]

    stability_rows = []
    for gene_id in rank_matrix.index:
        row = rank_matrix.loc[gene_id]
        min_rank, max_rank = int(row.min()), int(row.max())
        scenario_shift_max = max((abs(int(row[column]) - int(default[gene_id])) for column in rank_matrix if column != "default_weighted"), default=0)
        stability_rows.append({"gene_id": gene_id, "default_rank": int(default[gene_id]), "rank_min": min_rank,
                               "rank_max": max_rank, "rank_range": max_rank - min_rank,
                               "tie_size": int(tie_sizes[gene_id]),
                               "stability_class": _stability_class(gene_id, default, rank_matrix, metadata_frame, tie_sizes, top_k),
                               "instability_reason": distribution.set_index("gene_id").loc[gene_id, "instability_reason"],
                               "top_k_rate": round(float((row <= top_k).mean()), 6), "scenario_shift_max": scenario_shift_max})
    rank_stability = pd.DataFrame(stability_rows).sort_values(["default_rank", "gene_id"])

    perturbation_summary_rows = []
    for name in perturb_columns:
        body = name[len("perturb:"):]
        code, suffix = body.rsplit(":", 1)
        if suffix.startswith("minus"):
            direction, fraction = "minus", int(suffix[len("minus"):]) / 100
        else:
            direction, fraction = "plus", int(suffix[len("plus"):]) / 100
        shifts = (rank_matrix[name] - default).abs()
        worst_gene = str(shifts.idxmax()) if len(shifts) else ""
        perturbation_summary_rows.append({"evidence_code": code, "perturbation_fraction": fraction, "direction": direction,
                                          "max_rank_shift": int(shifts.max()) if len(shifts) else 0,
                                          "median_rank_shift": float(shifts.median()) if len(shifts) else 0.0,
                                          "candidates_affected": int((shifts > 0).sum()),
                                          "worst_affected_gene": worst_gene})
    weight_perturbation_summary = pd.DataFrame(perturbation_summary_rows)

    dropout_stability = pd.DataFrame(result["dropout_rows"])
    pd.DataFrame(concordance).to_csv(out / "rank_concordance.tsv", sep="\t", index=False)
    sensitivity.to_csv(out / "ranking_sensitivity.tsv", sep="\t", index=False)
    distribution.to_csv(out / "candidate_rank_distribution.tsv", sep="\t", index=False)
    consensus.to_csv(out / "consensus_shortlist.tsv", sep="\t", index=False)
    discordance.to_csv(out / "ranking_discordance_review.tsv", sep="\t", index=False)
    rank_stability.to_csv(out / "rank_stability.tsv", sep="\t", index=False)
    dropout_stability.to_csv(out / "source_dropout_stability.tsv", sep="\t", index=False)
    weight_perturbation_summary.to_csv(out / "weight_perturbation_summary.tsv", sep="\t", index=False)
    manifest = {"schema_version": "1.0", "benchmark_mode": benchmark_mode, "top_k": top_k,
                "perturbations": perturbations, "weight_registry": registry_meta, "scenario_count": len(scenarios),
                "scenario_categories": sorted({entry["category"] for entry in metadata}),
                "source_dropout_groups": source_dropout_groups, "provider_dropout_types": provider_dropout_types,
                "rank_shift_threshold": RANK_SHIFT_THRESHOLD,
                "consensus_rule": "Median of category-level ranks; correlated scenarios are collapsed before aggregation.",
                "stability_class_policy": (
                    "underdetermined_ties (default-rank tie, screening_tie_size > 1) > source_sensitive "
                    "(source_dropout crosses top_k; candidates inside top_k additionally on shift >= 3) "
                    "> provider_sensitive (provider_dropout, same rule) > weight_sensitive (weight perturbation, "
                    "confidence, eligibility, family dropout, same rule) > stable. "
                    "Classes are diagnostics and never feed back into scoring."
                ),
                "interpretation": "Rank + uncertainty: consensus reports robustness and never overwrites the deterministic default rank. "
                                  "Score differences are ordinal priorities, not probabilities; stability diagnostics are advisory and emit "
                                  "no recommended weights or hidden composite score."}
    (out / "ranking_qa_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return {"out_dir": str(out), "scenario_count": len(scenarios), "candidate_count": len(distribution),
            "stable_top_count": int((distribution["consensus_class"] == "stable_top").sum()),
            "review_count": len(discordance), "manifest": str(out / "ranking_qa_manifest.json"),
            "rank_stability": str(out / "rank_stability.tsv"),
            "source_dropout_stability": str(out / "source_dropout_stability.tsv"),
            "weight_perturbation_summary": str(out / "weight_perturbation_summary.tsv")}


def ranking_stability_benchmark(manifest_arg: str, out_dir_arg: str, top_k: int = 20) -> dict[str, object]:
    """Post-evaluate frozen prediction artifacts (issue #228).

    Each manifest case (``case_id``, ``candidates``, ``evidence``, ``ranked_tsv``,
    ``answer_gene_id``) is a frozen ranked artifact whose known-truth gene(s) are
    re-tested under weight perturbation and source/provider dropout. The frozen
    ``ranked_tsv`` supplies the default rank and is never rewritten; scenario
    ranks come from fresh runs against the same evidence. Cases with missing
    inputs or a truth gene absent from the frozen candidates are reported as
    ``not_evaluated`` and excluded from the retention denominator.

    This is evaluation-after-freeze for diagnostics, not a change to the
    audit-freeze contract (the frozen artifact stays byte-identical).
    """
    if top_k < 1:
        raise ValueError("top_k must be at least 1")
    manifest_path = runtime_path(manifest_arg)
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    source_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    registry, registry_meta = load_weight_registry()
    rows: list[dict[str, object]] = []
    not_evaluated: list[dict[str, object]] = []
    weights = load_weights(None)
    with tempfile.TemporaryDirectory(prefix="candidate-qtl-stability-benchmark-") as tmp:
        work = Path(tmp)
        for case in source_manifest.get("cases", []):
            case_id = str(case.get("case_id", ""))
            candidates_arg = str(case.get("candidates", ""))
            evidence_arg = str(case.get("evidence", ""))
            ranked_arg = str(case.get("ranked_tsv", ""))
            truths = _truth_ids(case.get("answer_gene_id"))
            missing = [name for name, value in (("candidates", candidates_arg), ("evidence", evidence_arg), ("ranked_tsv", ranked_arg))
                       if not value or not runtime_path(value).exists()]
            if missing:
                not_evaluated.append({"case_id": case_id, "reason": "missing_input:" + ",".join(missing)})
                continue
            ranked = pd.read_csv(runtime_path(ranked_arg), sep="\t", dtype=str).fillna("")
            truth_genes = [gene for gene in truths if gene in set(ranked["gene_id"])]
            if not truth_genes:
                not_evaluated.append({"case_id": case_id, "reason": "truth_not_in_candidates"})
                continue
            evidence = pd.read_csv(runtime_path(evidence_arg), sep="\t", dtype=str).fillna("")
            result = _run_stability_scenarios(runtime_path(candidates_arg), evidence, weights, top_k, None, work)
            # The scenario ranks come from a fresh run against `candidates`;
            # a truth gene absent there cannot be scored and must stay
            # not_evaluated instead of crashing the whole benchmark.
            scenario_index = next(iter(result["scenarios"].values())).index if result["scenarios"] else pd.Index([])
            absent_from_scenarios = [gene for gene in truth_genes if gene not in scenario_index]
            if absent_from_scenarios:
                not_evaluated.append({"case_id": case_id,
                                      "reason": "truth_not_in_candidates:" + ",".join(absent_from_scenarios)})
                continue
            ranked_index = ranked.set_index("gene_id")
            default_min = int(ranked_index.loc[truth_genes, "screening_rank"].astype(int).min())
            default_max = int(ranked_index.loc[truth_genes, "screening_rank_max"].astype(int).max())
            rows.append({"case_id": case_id, "scenario": "frozen_artifact", "category": "frozen",
                         "truth_gene_ids": "|".join(truth_genes),
                         "default_rank_min": default_min, "default_rank_max": default_max,
                         "scenario_rank_min": default_min, "scenario_rank_max": default_max,
                         "rank_shift_min": 0, "truth_retained_in_top_k": default_max <= top_k})
            for name in result["scenarios"]:
                scenario_series = result["scenarios"][name]
                scenario_max = result["scenarios_max"][name]
                entry = next(item for item in result["metadata"] if item["scenario"] == name)
                scenario_min_rank = int(scenario_series.loc[truth_genes].min())
                scenario_max_rank = int(scenario_max.loc[truth_genes].max())
                rows.append({"case_id": case_id, "scenario": name, "category": entry["category"],
                             "truth_gene_ids": "|".join(truth_genes),
                             "default_rank_min": default_min, "default_rank_max": default_max,
                             "scenario_rank_min": scenario_min_rank, "scenario_rank_max": scenario_max_rank,
                             "rank_shift_min": scenario_min_rank - default_min,
                             "truth_retained_in_top_k": scenario_max_rank <= top_k})

    detail = pd.DataFrame(rows)
    not_evaluated_frame = pd.DataFrame(not_evaluated)
    detail_path = out_dir / "stability_benchmark_cases.tsv"
    detail.to_csv(detail_path, sep="\t", index=False)
    not_evaluated_path = out_dir / "stability_benchmark_not_evaluated.tsv"
    if not_evaluated_frame.empty:
        not_evaluated_frame = pd.DataFrame(columns=["case_id", "reason"])
    not_evaluated_frame.to_csv(not_evaluated_path, sep="\t", index=False)

    summary_rows = []
    scenario_names = list(dict.fromkeys(detail["scenario"]))
    for scenario in scenario_names:
        group = detail[detail["scenario"].eq(scenario)]
        retained = group["truth_retained_in_top_k"].eq(True).sum()
        summary_rows.append({"scenario": scenario, "category": str(group["category"].iloc[0]),
                             "evaluated_case_count": int(len(group)),
                             "truth_retained_count": int(retained),
                             "truth_retained_rate": round(float(retained / len(group)), 6) if len(group) else None,
                             "not_evaluated_count": len(not_evaluated), "top_k": top_k})
    summary = pd.DataFrame(summary_rows)
    summary_path = out_dir / "stability_benchmark_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    evaluated = detail[detail["truth_retained_in_top_k"].notna()]
    default_row = summary[summary["scenario"].eq("frozen_artifact")]
    rate_rows = summary.dropna(subset=["truth_retained_rate"])
    worst_row = rate_rows.loc[rate_rows["truth_retained_rate"].idxmin()] if not rate_rows.empty else None
    evaluated_case_count = int(evaluated["case_id"].nunique()) if len(evaluated) else 0
    report_lines = [
        "# Ranking stability benchmark (frozen cases)",
        "",
        "Post-evaluation of frozen prediction artifacts (issue #228): the deterministic default rank is preserved, "
        "and this report only measures whether known-truth genes stay within Top-K under weight perturbation and "
        "source/provider dropout. It is not a causal probability and recommends no weight changes.",
        "",
        "## Result",
        "",
        f"- top_k: {top_k}; cases: {evaluated_case_count + len(not_evaluated)}; evaluated: {evaluated_case_count}; not_evaluated: {len(not_evaluated)}",
        f"- frozen-artifact truth retention (strict tie-aware, rank_max <= top_k): {float(default_row['truth_retained_rate'].iloc[0]):.3f}" if not default_row.empty else "- frozen-artifact truth retention: not evaluable",
        f"- worst scenario truth retention: {worst_row['scenario']} {float(worst_row['truth_retained_rate']):.3f}" if worst_row is not None else "- worst scenario truth retention: not evaluable",
        "",
        "## Interpretation",
        "",
        "Retention is strict: a truth gene counts as retained only when its worst tied rank (rank_max) is within top_k. "
        "not_evaluated cases (missing inputs, truth absent from the frozen candidate list) are excluded from every denominator "
        "and reported in stability_benchmark_not_evaluated.tsv.",
        "",
        "## Inputs",
        "",
        f"- source manifest: {manifest_path}",
        "- truth is read from answer_gene_id in the source manifest; frozen ranked artifacts are never rewritten.",
    ]
    report_path = out_dir / "stability_benchmark_report.md"
    report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    output_manifest = {
        "schema_version": "1.0",
        "benchmark": "ranking_stability_benchmark",
        "purpose": "Post-evaluation of frozen prediction artifacts: reports whether known-truth genes remain within Top-K under "
                   "weight perturbation, source-group dropout, and provider dropout. Evaluation-after-freeze; the frozen ranked "
                   "artifact is never rewritten and the audit-freeze contract is unchanged.",
        "source_benchmark_manifest": str(manifest_path),
        "source_benchmark_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        "weight_registry": registry_meta,
        "top_k": top_k,
        "scenarios": scenario_names,
        "evaluated_case_count": int(evaluated["case_id"].nunique()),
        "not_evaluated_case_count": len(not_evaluated),
        "not_evaluated": not_evaluated,
        "tie_policy": "strict: a truth gene is retained only when its worst rank (rank_max) is <= top_k",
        "denominator_policy": "not_evaluated cases are excluded from every retention denominator and reported separately",
        "probability_claim": False,
        "llm_used": False,
        "recommended_weight_registry_emitted": False,
        "recommendation": "diagnostic only; no weight changes are recommended from these outputs",
    }
    output_manifest_path = out_dir / "stability_benchmark_manifest.json"
    output_manifest_path.write_text(json.dumps(output_manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return {"cases": str(detail_path), "not_evaluated": str(not_evaluated_path), "summary": str(summary_path),
            "report": str(report_path), "manifest": str(output_manifest_path),
            "evaluated_case_count": evaluated_case_count, "not_evaluated_case_count": len(not_evaluated)}
