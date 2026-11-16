"""Interpretation coverage, kept apart from candidate support (issue #227).

A candidate with a rich Arabidopsis literature is easier for every route to
interpret than an orphan, a lineage-specific gene, a recent duplicate, or a
locus whose annotation is thin. Routes A, B and C all read external knowledge,
so they are all easier on well-studied genes -- and if that is not measured, the
knowledge-availability prior reads exactly like biological confidence.

This module answers a different question from "is this candidate good": **how
much could anything have been said about it?** For every candidate it records,
per dimension, one of three states that must never be collapsed:

``available``
    The resource exists for this candidate and was read.
``knowledge_unavailable``
    The resource was consulted and holds nothing for this candidate. No
    literature, no reference function, no homolog -- *unknown*, not *unlikely*.
``not_evaluated``
    The provider was not supplied at all, so nothing was consulted. Absence of a
    provider is not absence of knowledge.

Neither of the last two is negative evidence, and neither may demote a
candidate. A dark gene stays visible; what changes is that the report can say
*why* it is dark, and a reviewer can see that a low-interpretation candidate is
an escalation target rather than a rejected one.

``benchmark-stratify`` closes the loop: given any per-case benchmark table and a
per-case stratum table (coverage tier, lineage specificity, interval size band),
it recomputes each metric inside each stratum, with the same per-metric
denominators. Recall over well-studied cases and recall over poorly annotated
cases are two numbers, and reporting only their average hides the whole problem.
"""

from __future__ import annotations

import json

import pandas as pd

from . import __version__
from .core import file_fingerprint, runtime_path

SCHEMA_VERSION = "1.0"

COVERAGE_STATES = ("available", "knowledge_unavailable", "not_evaluated")
COVERAGE_TIERS = ("well_characterised", "partially_characterised", "knowledge_limited", "dark")

# One dimension per external resource a route can read. The provider argument
# that supplies it decides between knowledge_unavailable and not_evaluated.
DIMENSIONS = ("annotation", "family", "homolog_orthology", "reference_function",
              "literature", "ontology", "expression", "variant")

COVERAGE_COLUMNS = [
    "gene_id", "family_id",
    *[f"{dimension}_coverage" for dimension in DIMENSIONS],
    "available_dimension_count", "evaluated_dimension_count", "knowledge_coverage_tier",
    "limiting_dimensions", "lineage_specificity", "review_state", "rank_contribution", "notes",
]
FAMILY_COLUMNS = ["family_id", "gene_count", "available_dimension_mean", "knowledge_coverage_tier",
                  "genes_knowledge_limited", "genes_dark"]
STRATUM_COLUMNS = ["metric", "stratum_column", "stratum_value", "cases", "evaluated", "hits",
                   "recall", "not_evaluated"]

LINEAGE_STATES = ("has_cross_species_homolog", "focal_lineage_specific_candidate", "not_evaluated")


def _read(path_arg: str | None, required: tuple[str, ...] = ()) -> pd.DataFrame | None:
    if not path_arg:
        return None
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{path_arg} is missing columns: {', '.join(missing)}")
    return frame


def _state(provider: pd.DataFrame | None, has_value: bool) -> str:
    """Three states, never two: no provider is not an empty provider."""
    if provider is None:
        return "not_evaluated"
    return "available" if has_value else "knowledge_unavailable"


def _nonempty(row: dict[str, object], *columns: str) -> bool:
    return any(str(row.get(column, "")).strip() for column in columns)


def knowledge_coverage(
    candidates_arg: str,
    out_dir_arg: str,
    atlas_arg: str | None = None,
    orthology_arg: str | None = None,
    reference_function_arg: str | None = None,
    literature_arg: str | None = None,
    expression_arg: str | None = None,
    variants_arg: str | None = None,
    focal_species: str = "",
    full_hash: bool = False,
) -> dict[str, object]:
    """Per-candidate interpretation coverage, per dimension, in three states."""
    candidates_path = runtime_path(candidates_arg)
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    candidates = pd.read_csv(candidates_path, sep="\t", dtype=str).fillna("")
    if "gene_id" not in candidates.columns:
        raise ValueError("knowledge-coverage requires a candidates table with gene_id")

    atlas = _read(atlas_arg, ("gene_id",))
    orthology = _read(orthology_arg, ("gene_id",))
    reference = _read(reference_function_arg)
    literature = _read(literature_arg)
    expression = _read(expression_arg, ("gene_id",))
    variants = _read(variants_arg, ("gene_id",))

    atlas_rows = {str(row["gene_id"]): row for row in atlas.to_dict("records")} if atlas is not None else {}
    orthology_by_gene: dict[str, list[dict[str, object]]] = {}
    if orthology is not None:
        for row in orthology.to_dict("records"):
            orthology_by_gene.setdefault(str(row["gene_id"]), []).append(row)
    reference_by_id: dict[str, dict[str, object]] = {}
    if reference is not None:
        for row in reference.to_dict("records"):
            for key in ("reference_gene_id", "gene_id"):
                value = str(row.get(key, "")).strip()
                if value:
                    reference_by_id.setdefault(value.casefold(), row)
    literature_genes: set[str] = set()
    literature_references: set[str] = set()
    if literature is not None:
        for row in literature.to_dict("records"):
            if str(row.get("query_status", "")).strip() not in ("", "retrieved"):
                continue
            if not str(row.get("source_id", "")).strip():
                continue
            for value in str(row.get("interval_gene_ids", "")).split(";"):
                if value.strip():
                    literature_genes.add(value.strip().casefold())
            reference_gene = str(row.get("reference_gene", "")).strip()
            if reference_gene:
                literature_references.add(reference_gene.casefold())
    expression_genes = {str(row["gene_id"]).casefold() for row in expression.to_dict("records")} \
        if expression is not None else set()
    variant_genes = {str(row["gene_id"]).casefold() for row in variants.to_dict("records")} \
        if variants is not None else set()

    rows: list[dict[str, object]] = []
    for record in candidates.to_dict("records"):
        gene_id = str(record["gene_id"]).strip()
        atlas_row = atlas_rows.get(gene_id, {})
        homologs = orthology_by_gene.get(gene_id, [])
        reference_ids = [str(item.get("reference_gene_id", "")).strip() for item in homologs
                         if str(item.get("reference_gene_id", "")).strip()]
        if not reference_ids and atlas_row:
            value = str(atlas_row.get("orthology_reference_gene", "")).strip()
            reference_ids = [value] if value else []
        reference_rows = [reference_by_id[key.casefold()] for key in reference_ids
                          if key.casefold() in reference_by_id]
        # The candidate's own id may itself be a curated reference (focal-species
        # annotation), which is a legitimate reference-function source.
        if gene_id.casefold() in reference_by_id:
            reference_rows.append(reference_by_id[gene_id.casefold()])

        coverage = {
            "annotation": _state(candidates, _nonempty(record, "description", "product", "annotation")),
            "family": _state(atlas if atlas is not None else None,
                             _nonempty(atlas_row, "family_id", "evolutionary_family_id", "orthogroup")),
            "homolog_orthology": _state(orthology if orthology is not None else (atlas if atlas is not None else None),
                                        bool(reference_ids)),
            "reference_function": _state(reference, any(
                _nonempty(row, "function_description", "mutant_phenotype", "pathway") for row in reference_rows)),
            "literature": _state(
                literature if literature is not None else reference,
                gene_id.casefold() in literature_genes
                or any(key.casefold() in literature_references for key in reference_ids)
                or any(_nonempty(row, "phenotype_source_ids", "literature_source_ids") for row in reference_rows)),
            "ontology": _state(reference if reference is not None else (atlas if atlas is not None else None),
                               any(_nonempty(row, "go_terms", "pathway_id") for row in reference_rows)
                               or _nonempty(atlas_row, "domain_ids", "broad_function_bins")),
            "expression": _state(expression, gene_id.casefold() in expression_genes),
            "variant": _state(variants, gene_id.casefold() in variant_genes),
        }
        available = sum(1 for state in coverage.values() if state == "available")
        evaluated = sum(1 for state in coverage.values() if state != "not_evaluated")
        unavailable = [dimension for dimension, state in coverage.items() if state == "knowledge_unavailable"]

        # Tiers describe how much could be said, never how good the candidate is.
        interpretable = sum(1 for dimension in ("family", "homolog_orthology", "reference_function", "literature")
                            if coverage[dimension] == "available")
        if interpretable >= 3:
            tier = "well_characterised"
        elif interpretable == 2:
            tier = "partially_characterised"
        elif interpretable == 1:
            tier = "knowledge_limited"
        else:
            tier = "dark"

        if coverage["homolog_orthology"] == "available":
            lineage = "has_cross_species_homolog"
        elif coverage["homolog_orthology"] == "knowledge_unavailable":
            lineage = "focal_lineage_specific_candidate"
        else:
            lineage = "not_evaluated"

        review_state = ""
        if tier in ("knowledge_limited", "dark") and evaluated:
            # The limitation is coverage, not contradiction: this is an
            # escalation target (annotation rescue, orthology, wet assay), and
            # explicitly not a demotion.
            review_state = "escalate_knowledge_limited_candidate"
        rows.append({
            "gene_id": gene_id,
            "family_id": str(atlas_row.get("family_id", "")) or str(record.get("family_id", "")),
            **{f"{dimension}_coverage": coverage[dimension] for dimension in DIMENSIONS},
            "available_dimension_count": available, "evaluated_dimension_count": evaluated,
            "knowledge_coverage_tier": tier, "limiting_dimensions": ";".join(sorted(unavailable)),
            "lineage_specificity": lineage, "review_state": review_state,
            "rank_contribution": "none",
            "notes": "coverage describes what could be interpreted, never how likely the candidate is; "
                     "knowledge_unavailable and not_evaluated are both non-evidence and never demote a candidate",
        })

    coverage_frame = pd.DataFrame(rows, columns=COVERAGE_COLUMNS)
    coverage_path = out / "knowledge_coverage.tsv"
    coverage_frame.to_csv(coverage_path, sep="\t", index=False, lineterminator="\n")

    family_rows: list[dict[str, object]] = []
    for family, group in coverage_frame[coverage_frame["family_id"].astype(str) != ""].groupby("family_id"):
        mean_available = round(float(group["available_dimension_count"].astype(float).mean()), 6)
        limited = int((group["knowledge_coverage_tier"] == "knowledge_limited").sum())
        dark = int((group["knowledge_coverage_tier"] == "dark").sum())
        family_rows.append({
            "family_id": family, "gene_count": int(len(group)), "available_dimension_mean": mean_available,
            "knowledge_coverage_tier": group["knowledge_coverage_tier"].mode().iloc[0],
            "genes_knowledge_limited": limited, "genes_dark": dark,
        })
    family_frame = pd.DataFrame(family_rows, columns=FAMILY_COLUMNS).sort_values("family_id", kind="stable")
    family_path = out / "knowledge_coverage_by_family.tsv"
    family_frame.to_csv(family_path, sep="\t", index=False, lineterminator="\n")

    tier_counts = {tier: int((coverage_frame["knowledge_coverage_tier"] == tier).sum()) for tier in COVERAGE_TIERS}
    manifest = {
        "schema_version": SCHEMA_VERSION, "command": "knowledge-coverage",
        "candidate_qtl_version": __version__,
        "focal_species": focal_species,
        "candidate_count": int(len(coverage_frame)),
        "dimensions": list(DIMENSIONS),
        "coverage_states": list(COVERAGE_STATES),
        "tier_counts": tier_counts,
        "lineage_counts": {state: int((coverage_frame["lineage_specificity"] == state).sum())
                           for state in LINEAGE_STATES},
        "dimension_state_counts": {
            dimension: {state: int((coverage_frame[f"{dimension}_coverage"] == state).sum())
                        for state in COVERAGE_STATES}
            for dimension in DIMENSIONS},
        "providers": {
            "atlas": bool(atlas_arg), "orthology": bool(orthology_arg),
            "reference_function": bool(reference_function_arg), "literature": bool(literature_arg),
            "expression": bool(expression_arg), "variants": bool(variants_arg),
        },
        "rank_contribution": "none",
        "rules": [
            "knowledge_unavailable means the resource was consulted and holds nothing: unknown, not unlikely",
            "not_evaluated means no provider was supplied: absence of a provider is not absence of knowledge",
            "neither state is negative evidence and neither demotes a candidate",
            "a knowledge-limited or dark candidate gets an escalation review state, never a lower rank",
            "coverage tiers describe interpretability, never candidate quality",
        ],
        "inputs": {"candidates": file_fingerprint(candidates_path, full_hash)},
        "outputs": {"coverage": str(coverage_path), "by_family": str(family_path)},
    }
    manifest_path = out / "knowledge_coverage_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), "coverage": str(coverage_path), "by_family": str(family_path),
            "candidate_count": manifest["candidate_count"], "tier_counts": tier_counts,
            "escalation_count": int((coverage_frame["review_state"] != "").sum()),
            "rank_contribution": "none"}


# --------------------------------------------------------------------------- #
# stratified benchmark reporting
# --------------------------------------------------------------------------- #

def benchmark_stratify(
    cases_arg: str,
    strata_arg: str,
    out_dir_arg: str,
    metric_args: list[str] | None = None,
    stratum_column_args: list[str] | None = None,
    case_column: str = "case_id",
    full_hash: bool = False,
) -> dict[str, object]:
    """Recompute per-case benchmark metrics inside each stratum.

    A single recall number over a mixed case set answers the wrong question when
    half the cases are well-studied Arabidopsis biology and half are orphans.
    The per-metric denominator rule is unchanged: an empty per-case value is
    ``not_evaluated`` and leaves the denominator rather than counting as a miss.
    """
    cases_path = runtime_path(cases_arg)
    strata_path = runtime_path(strata_arg)
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    cases = pd.read_csv(cases_path, sep="\t", dtype=str).fillna("")
    strata = pd.read_csv(strata_path, sep="\t", dtype=str).fillna("")
    for name, frame in (("cases", cases), ("strata", strata)):
        if case_column not in frame.columns:
            raise ValueError(f"{name} table is missing the case column {case_column!r}")
    if strata[case_column].duplicated().any():
        raise ValueError("the strata table must hold one row per case")

    boolean_like = [column for column in cases.columns
                    if column != case_column
                    and set(cases[column].astype(str).str.strip().unique()) <= {"True", "False", ""}
                    and cases[column].astype(str).str.strip().ne("").any()]
    metrics = [metric.strip() for metric in (metric_args or []) if metric.strip()] or boolean_like
    unknown = [metric for metric in metrics if metric not in cases.columns]
    if unknown:
        raise ValueError(f"metric column(s) not in the case table: {', '.join(unknown)}")
    if not metrics:
        raise ValueError("no boolean metric column was found; pass --metric explicitly")
    stratum_columns = [column.strip() for column in (stratum_column_args or []) if column.strip()] \
        or [column for column in strata.columns if column != case_column]
    unknown_strata = [column for column in stratum_columns if column not in strata.columns]
    if unknown_strata:
        raise ValueError(f"stratum column(s) not in the strata table: {', '.join(unknown_strata)}")

    # The case table keeps its own columns when the strata table repeats one
    # (e.g. evidence-tier-benchmark strata re-derive interval_size_band from
    # the cases they stratify), so the merge cannot produce _x/_y suffixed
    # duplicates that later column lookups would miss.
    merged = cases.merge(strata, on=case_column, how="left", suffixes=("", "_strata")).fillna("")
    unmatched = int((merged[stratum_columns[0]].astype(str).str.strip() == "").sum()) if stratum_columns else 0

    rows: list[dict[str, object]] = []
    for metric in metrics:
        for column in stratum_columns:
            values = merged[column].astype(str).str.strip()
            for value in sorted({item for item in values}):
                subset = merged[values == (value or "")]
                metric_values = subset[metric].astype(str).str.strip()
                evaluated = int(metric_values.ne("").sum())
                hits = int(metric_values.eq("True").sum())
                rows.append({
                    "metric": metric, "stratum_column": column,
                    "stratum_value": value or "(unstratified)",
                    "cases": int(len(subset)), "evaluated": evaluated, "hits": hits,
                    "recall": round(hits / evaluated, 6) if evaluated else "",
                    "not_evaluated": int(len(subset)) - evaluated,
                })
    summary = pd.DataFrame(rows, columns=STRATUM_COLUMNS)
    summary_path = out / "benchmark_strata_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False, lineterminator="\n")

    manifest = {
        "schema_version": SCHEMA_VERSION, "command": "benchmark-stratify",
        "candidate_qtl_version": __version__,
        "case_count": int(len(cases)), "cases_without_a_stratum": unmatched,
        "metrics": metrics, "stratum_columns": stratum_columns,
        "denominator_policy": (
            "inside every stratum a metric divides by the cases that could be evaluated; an empty per-case value "
            "is not_evaluated and leaves the denominator instead of counting as a miss"
        ),
        "reporting_rule": (
            "strata are reported side by side and never averaged into one headline number: recall over "
            "well-studied cases and recall over knowledge-limited cases are two different measurements"
        ),
        "inputs": {"cases": file_fingerprint(cases_path, full_hash),
                   "strata": file_fingerprint(strata_path, full_hash)},
        "outputs": {"summary": str(summary_path)},
    }
    manifest_path = out / "benchmark_strata_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), "summary": str(summary_path),
            "case_count": manifest["case_count"], "metrics": metrics,
            "stratum_columns": stratum_columns, "cases_without_a_stratum": unmatched}
