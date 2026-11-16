"""Route C retrieval benchmark (issue #193, Phase 1).

Measures the phenotype-only, interval-blind retrieval power of Route C against
frozen known-gene cases:

- **Family Recall@5 / @10 / @20**: is the truth family / orthogroup among the
  top-K distinct families of the frozen retrieval?
- **Gene Recall@20**: is the truth reference gene among the top-20 retrieval
  rows?
- **A-only / C-only / A+C union**: how much does the embedding retrieval (C)
  recover on its own, and how much does it add on top of the symbolic / literature
  route (A)?

The benchmark is interval-blind by construction: the manifest never carries an
interval gene list, and the retrieval tables are the frozen Route C artifacts
produced by ``retrieval-import`` (which rejects any extra column). The
``route_a_tsv`` input is the Route A hypothesis table (``gene_family`` column);
it is optional, and A-only/C-only/union columns are left empty when it is
absent.

**Per-metric denominators (issue #214).** ``truth_gene_id`` and ``route_a_tsv``
are documented as optional, so a manifest legitimately mixes family-only truth
with gene truth and cases with/without a Route A artifact. A metric that could
not be evaluated because its optional input is absent is *not* a miss: it is
``not_evaluated``. Every metric therefore reports its own ``evaluated``
denominator, ``hits``, ``recall`` (``None`` when nothing was evaluable) and
``not_evaluated`` count. Absent provider/input data is never converted into
negative biological evidence, which is the same rule Route B and Route C apply
to a missing provider.
"""

from __future__ import annotations

import json

import pandas as pd

from .benchmark_reporting import conditional_binary_performance, wilson_interval
from .core import file_fingerprint, runtime_path, safe_id
from .family_labels import canonical_family_key, families_match, family_set_contains
from .retrieval_route import validate_retrieval_frame

BENCHMARK_VERSION = "1.2.0"
MANIFEST_COLUMNS = ["case_id", "retrieval_tsv", "truth_family_id"]
OPTIONAL_COLUMNS = ("truth_gene_id", "route_a_tsv", "answer_derived")

# Which manifest input a metric needs before a case counts towards it. A case
# missing that input stays empty in the case table and is counted as
# ``not_evaluated``, never as a zero-valued miss. ``c_only_hit_20`` is a
# comparative claim ("C reached it and A did not"), so it needs Route A too.
METRIC_REQUIREMENTS = {
    "family_recall_5": "truth_family_id",
    "family_recall_10": "truth_family_id",
    "family_recall_20": "truth_family_id",
    "gene_recall_20": "truth_gene_id",
    "a_only_hit": "route_a_tsv",
    "c_only_hit_20": "route_a_tsv",
    "ac_union_hit_5": "route_a_tsv",
    "ac_union_hit_10": "route_a_tsv",
    "ac_union_hit_20": "route_a_tsv",
}
SUMMARY_COLUMNS = ["metric", "requires", "cases", "evaluated", "hits", "recall", "not_evaluated",
                   "clean_cases", "contaminated_cases", "evaluation_coverage", "wilson_95_low", "wilson_95_high"]


def _top_families(frame: pd.DataFrame, k: int) -> set[str]:
    """Distinct families/orthogroups ranked by their best row rank."""
    frame = frame.sort_values("rank", kind="stable")
    families: list[str] = []
    seen: set[str] = set()
    for row in frame.to_dict("records"):
        for value in (str(row.get("family_id", "")).strip(), str(row.get("orthogroup_id", "")).strip()):
            if value:
                key = canonical_family_key(value)
                if key not in seen:
                    seen.add(key)
                    families.append(value)
        if len(families) >= k:
            break
    return set(families)


def _top_genes(frame: pd.DataFrame, k: int) -> set[str]:
    return set(frame.sort_values("rank", kind="stable")["reference_gene_id"].head(k).astype(str).str.strip())


def _c_family_rank(frame: pd.DataFrame, truth_family: str) -> int | None:
    best: int | None = None
    for row in frame.to_dict("records"):
        for value in (str(row.get("family_id", "")).strip(), str(row.get("orthogroup_id", "")).strip()):
            if value and families_match(value, truth_family):
                best = int(row["rank"]) if best is None else min(best, int(row["rank"]))
    return best


def route_c_benchmark(manifest_arg: str, out_arg: str, full_hash: bool = False) -> dict[str, object]:
    """Score frozen Route C retrieval artifacts on known-gene cases."""
    manifest_path = runtime_path(manifest_arg)
    out = runtime_path(out_arg)
    out.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str).fillna("")
    missing = [column for column in MANIFEST_COLUMNS if column not in manifest.columns]
    if missing:
        raise ValueError(f"route C benchmark manifest missing columns: {', '.join(missing)}")
    if manifest.empty:
        raise ValueError("route C benchmark manifest is empty")
    if manifest["case_id"].duplicated().any():
        raise ValueError("route C benchmark case_id must be unique")
    if manifest["retrieval_tsv"].astype(str).str.strip().eq("").any():
        raise ValueError("route C benchmark requires a retrieval_tsv per case")
    if manifest["truth_family_id"].astype(str).str.strip().eq("").any():
        raise ValueError("route C benchmark requires a truth_family_id per case")

    config_dir = manifest_path.parent

    def resolve(value: str) -> str | None:
        text = str(value).strip()
        if not text:
            return None
        path = runtime_path(text)
        if not path.is_absolute():
            path = config_dir / path
        if not path.exists():
            raise ValueError(f"route C benchmark input does not exist: {path}")
        return str(path)

    rows: list[dict[str, object]] = []
    for record in manifest.to_dict("records"):
        case_id = str(record["case_id"]).strip()
        retrieval_path = resolve(record["retrieval_tsv"])
        truth_family = str(record["truth_family_id"]).strip()
        truth_gene = str(record.get("truth_gene_id", "")).strip()
        route_a_path = resolve(record.get("route_a_tsv", ""))
        retrieval = validate_retrieval_frame(pd.read_csv(runtime_path(retrieval_path), sep="\t", dtype=str).fillna(""), f"retrieval for case {case_id}")
        retrieval["rank"] = pd.to_numeric(retrieval["rank"], errors="coerce").astype(int)
        a_families: set[str] = set()
        if route_a_path:
            route_a = pd.read_csv(runtime_path(route_a_path), sep="\t", dtype=str).fillna("")
            if "gene_family" not in route_a.columns:
                raise ValueError(f"route A table for case {case_id} is missing gene_family")
            a_families = {str(value).strip() for value in route_a["gene_family"] if str(value).strip()}

        family_rank = _c_family_rank(retrieval, truth_family)
        top = {k: _top_families(retrieval, k) for k in (5, 10, 20)}
        hit = {k: family_set_contains(top[k], truth_family) for k in (5, 10, 20)}
        gene_hit = truth_gene and truth_gene in _top_genes(retrieval, 20)
        a_hit = family_set_contains(a_families, truth_family)
        rows.append({
            "case_id": case_id, "truth_family_id": truth_family, "truth_gene_id": truth_gene,
            "retrieval_rows": int(len(retrieval)),
            "c_family_rank": "" if family_rank is None else family_rank,
            "family_recall_5": hit[5], "family_recall_10": hit[10], "family_recall_20": hit[20],
            "gene_recall_20": bool(gene_hit) if truth_gene else "",
            "a_only_hit": bool(a_hit) if route_a_path else "",
            "c_only_hit_20": bool(hit[20] and not a_hit) if route_a_path else "",
            "ac_union_hit_5": bool(a_hit or hit[5]) if route_a_path else "",
            "ac_union_hit_10": bool(a_hit or hit[10]) if route_a_path else "",
            "ac_union_hit_20": bool(a_hit or hit[20]) if route_a_path else "",
        })
    cases = pd.DataFrame(rows).sort_values("case_id", kind="stable")
    cases_path = out / "route_c_benchmark_cases.tsv"
    cases.to_csv(cases_path, sep="\t", index=False, lineterminator="\n")

    n = len(cases)

    def stats(column: str) -> dict[str, object]:
        """Per-metric denominator: only cases whose required optional input was
        present are evaluated. An empty per-case value stays ``not_evaluated``
        instead of being coerced to a zero-valued miss (issue #214)."""
        values = cases[column].astype(str).str.strip()
        evaluated = int(values.ne("").sum())
        hits = int(values.eq("True").sum())
        low, high = wilson_interval(hits, evaluated) if evaluated else (None, None)
        return {
            "requires": METRIC_REQUIREMENTS[column],
            "evaluated": evaluated,
            "hits": hits,
            "recall": round(hits / evaluated, 6) if evaluated else None,
            "not_evaluated": n - evaluated,
            # Issue #367: coverage columns. These commands have no contamination
            # concept, so clean = evaluated and contaminated = 0.
            "clean_cases": evaluated,
            "contaminated_cases": 0,
            "evaluation_coverage": (evaluated / n) if n else None,
            "wilson_95_low": low,
            "wilson_95_high": high,
        }

    metrics = {column: stats(column) for column in METRIC_REQUIREMENTS}

    def coverage_block(column: str) -> dict[str, object]:
        """conditional_binary_performance-style counts; evaluated -> clean."""
        statuses = cases[column].astype(str).str.strip()
        frame = pd.DataFrame({
            "case_id": cases["case_id"],
            "evaluation_status": statuses.map(lambda value: "clean" if value else "not_evaluated"),
            "hit": statuses.eq("True"),
        })
        return conditional_binary_performance(frame)

    summary_path = out / "route_c_benchmark_summary.tsv"
    pd.DataFrame(
        [{"metric": column, "requires": record["requires"], "cases": n, "evaluated": record["evaluated"],
          "hits": record["hits"], "recall": "" if record["recall"] is None else record["recall"],
          "not_evaluated": record["not_evaluated"],
          "clean_cases": record["clean_cases"], "contaminated_cases": record["contaminated_cases"],
          "evaluation_coverage": "" if record["evaluation_coverage"] is None else record["evaluation_coverage"],
          "wilson_95_low": "" if record["wilson_95_low"] is None else record["wilson_95_low"],
          "wilson_95_high": "" if record["wilson_95_high"] is None else record["wilson_95_high"]}
         for column, record in metrics.items()],
        columns=SUMMARY_COLUMNS,
    ).to_csv(summary_path, sep="\t", index=False, lineterminator="\n")

    manifest_out = {
        "schema_version": "1.1", "benchmark_version": BENCHMARK_VERSION,
        "benchmark_type": "Route C embedding retrieval, phenotype-only and interval-blind",
        "case_count": n,
        "metrics": metrics,
        "evaluation_coverage": {column: coverage_block(column) for column in METRIC_REQUIREMENTS},
        "denominator_policy": (
            "each metric divides by the cases that carried its required optional input; a case without "
            "truth_gene_id or without a usable route_a_tsv is not_evaluated for that metric and is never "
            "counted as a 0% miss. Missing provider/input data is not negative biological evidence"
        ),
        "not_evaluated_totals": {column: record["not_evaluated"] for column, record in metrics.items()},
        "interval_blind": {
            "declaration": "the manifest carries no interval gene list; retrieval inputs are the frozen Route C artifacts, which reject any column outside the retrieval contract",
            "interval_arguments": [],
        },
        "scoring_rule": "retrieval similarity is scored as a discovery route, never as causal evidence; scores never enter the deterministic rank",
        "route_a_input": "optional; when absent the A-only / C-only / union columns stay empty and those metrics are not_evaluated",
        "manifest_fingerprint": file_fingerprint(manifest_path, full_hash),
        "outputs": {"cases": str(cases_path), "summary": str(summary_path),
                    "manifest": str(out / "route_c_benchmark_manifest.json"), "report": str(out / "route_c_benchmark_report.md")},
    }
    manifest_json = out / "route_c_benchmark_manifest.json"
    manifest_json.write_text(json.dumps(manifest_out, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write_report(cases, manifest_out, out)
    return {"cases": str(cases_path), "summary": str(summary_path), "manifest": str(manifest_json),
            "report": str(out / "route_c_benchmark_report.md"), "case_count": n, "metrics": metrics}


def _rate(metrics: dict[str, object], name: str) -> str:
    """Render a metric with its own denominator. A metric that was never
    evaluable prints ``not_evaluated``, never 0% (issue #214)."""
    record = metrics[name]
    assert isinstance(record, dict)
    if record["recall"] is None:
        return f"not_evaluated (no case supplied {record['requires']})"
    return f"{float(record['recall']):.1%} ({record['hits']}/{record['evaluated']})"


def _write_report(cases: pd.DataFrame, manifest: dict[str, object], out) -> None:
    metrics = manifest["metrics"]
    assert isinstance(metrics, dict)
    lines = [
        "# Route C retrieval benchmark (phenotype-only, interval-blind)",
        "",
        "Route C scores the embedding retrieval as an independent discovery route.",
        "Similarity is never causal evidence and never enters a gene score.",
        "",
        "Every rate carries its own evaluated denominator. A metric whose optional input",
        "(truth_gene_id / route_a_tsv) was absent is reported as not_evaluated, never as a",
        "0% miss: missing input is not negative biological evidence.",
        "",
        f"- cases: {manifest['case_count']}",
        f"- Family Recall@5 / @10 / @20: {_rate(metrics, 'family_recall_5')} / "
        f"{_rate(metrics, 'family_recall_10')} / {_rate(metrics, 'family_recall_20')}",
        f"- Gene Recall@20: {_rate(metrics, 'gene_recall_20')}",
        f"- A-only / C-only (top-20) / A+C union (top-20): {_rate(metrics, 'a_only_hit')} / "
        f"{_rate(metrics, 'c_only_hit_20')} / {_rate(metrics, 'ac_union_hit_20')}",
        "", "## Per case", "",
        "An empty cell means the case did not supply that metric's optional input.", "",
        "| case | truth family | C rank | F@5 | F@10 | F@20 | gene@20 | A-only | C-only | A+C@20 |", "|---|---|---:|---|---|---|---|---|---|---:|"]
    for row in cases.to_dict("records"):
        lines.append(f"| {row['case_id']} | {row['truth_family_id']} | {row['c_family_rank']} | "
                     f"{row['family_recall_5']} | {row['family_recall_10']} | {row['family_recall_20']} | "
                     f"{row['gene_recall_20']} | {row['a_only_hit']} | {row['c_only_hit_20']} | {row['ac_union_hit_20']} |")
    lines += ["", "The value of Route C is judged by whether it recovers the truth family in a failure mode "
              "where A and B are weak, not by pooled numbers."]
    (out / "route_c_benchmark_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


# --------------------------------------------------------------------------- #
# provider comparison (issue #225, #5): the smallest possible layer on top of
# the unmodified route_c_benchmark() above -- it runs that function once per
# provider over the same frozen cases and places the per-metric denominators
# side by side. It never restructures route_c_benchmark()'s own scoring, and
# it never pools a provider's raw retrieval_score with another provider's:
# only each provider's own evaluated / hits / recall / not_evaluated counts
# (booleans and their denominators, never a raw similarity number) are
# compared.
# --------------------------------------------------------------------------- #
PROVIDER_COMPARISON_SUMMARY_COLUMNS = [
    "provider_label", "metric", "requires", "cases", "evaluated", "hits", "recall", "not_evaluated",
]


def route_c_provider_comparison(
    provider_manifests: dict[str, str],
    out_arg: str,
    full_hash: bool = False,
) -> dict[str, object]:
    """Run ``route_c_benchmark`` once per provider (e.g. lexical / semantic /
    hybrid) over the same frozen known-gene cases, and report a side-by-side
    Family Recall@K / Gene Recall@K comparison.

    Each ``provider_manifests[label]`` is an ordinary ``route_c_benchmark``
    manifest (``case_id``/``retrieval_tsv``/``truth_family_id`` [+optional
    ``truth_gene_id``/``route_a_tsv``]) whose ``retrieval_tsv`` values are
    that provider's own frozen Route C retrieval for the same cases. Nothing
    about ``route_c_benchmark``'s per-case scoring changes here: this layer
    only calls it once per label and juxtaposes each provider's own
    denominators (``evaluated``/``hits``/``recall``/``not_evaluated``) --
    booleans and counts, never a provider's raw ``retrieval_score`` -- so
    providers are never pooled into one number.
    """
    if len(provider_manifests) < 2:
        raise ValueError("route-c-provider-comparison requires at least two --provider-manifest inputs")
    out = runtime_path(out_arg)
    out.mkdir(parents=True, exist_ok=True)

    per_provider: dict[str, dict[str, object]] = {}
    rows: list[dict[str, object]] = []
    for label in sorted(provider_manifests):
        manifest_arg = provider_manifests[label]
        result = route_c_benchmark(manifest_arg, str(out / safe_id(label)), full_hash)
        per_provider[label] = result
        metrics = result["metrics"]
        assert isinstance(metrics, dict)
        for metric, record in metrics.items():
            assert isinstance(record, dict)
            rows.append({
                "provider_label": label, "metric": metric, "requires": METRIC_REQUIREMENTS[metric],
                "cases": result["case_count"], "evaluated": record["evaluated"], "hits": record["hits"],
                "recall": "" if record["recall"] is None else record["recall"],
                "not_evaluated": record["not_evaluated"],
            })

    comparison = pd.DataFrame(rows, columns=PROVIDER_COMPARISON_SUMMARY_COLUMNS)
    comparison_path = out / "route_c_provider_comparison.tsv"
    comparison.to_csv(comparison_path, sep="\t", index=False, lineterminator="\n")

    manifest_out = {
        "schema_version": "1.0",
        "comparison_type": "side-by-side Family/Gene Recall@K across Route C providers on the same frozen cases",
        "providers": sorted(provider_manifests),
        "cross_provider_raw_scores_compared": False,
        "cross_provider_raw_scores_pooled": False,
        "policy": (
            "each provider's benchmark is run independently by the unmodified route_c_benchmark(); only each "
            "provider's own evaluated/hits/recall/not_evaluated denominators are placed side by side, never a "
            "merged score, and no provider's raw retrieval_score is read by this comparison layer"
        ),
        "per_provider_outputs": {label: result["manifest"] for label, result in per_provider.items()},
        "outputs": {"comparison": str(comparison_path), "manifest": str(out / "route_c_provider_comparison_manifest.json")},
    }
    manifest_path = out / "route_c_provider_comparison_manifest.json"
    manifest_path.write_text(json.dumps(manifest_out, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"comparison": str(comparison_path), "manifest": str(manifest_path),
            "providers": sorted(provider_manifests), "per_provider": {label: result["case_count"] for label, result in per_provider.items()}}
