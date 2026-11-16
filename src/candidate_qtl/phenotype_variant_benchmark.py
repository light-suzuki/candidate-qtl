"""Phenotype wording robustness benchmark (issue #235, Phase 1).

Measures how much the discovery routes' outcomes depend on the exact wording
and description granularity of the phenotype: the same truth case is re-scored
across **frozen phenotype formulations** (variants) that differ only in how the
phenotype is phrased (original / short / synonymized / mechanism terms removed /
irrelevant detail added). Nothing is re-run here and no prediction logic is
changed: the benchmark reads the frozen per-variant artifacts that external
lane execution produced for each formulation and measures stability.

Per case the manifest names:

- ``variants_tsv`` -- the frozen formulation table (variant_id / variant_kind /
  phenotype_text). This table *is* the freeze; the benchmark never generates a
  variant and records the table's SHA-256.
- per-variant artifact patterns: ``route_a_families_tsv`` (Route A proposal
  families), ``route_c_retrieval_tsv`` (Route C frozen retrieval),
  ``route_b_overlay_status_json`` (the b1 ``interval_first_status.json`` of
  ``interval-reverse-lookup --phenotype``), ``route_b_atlas_manifest_json``
  (the interval-atlas manifest whose ``atlas_sha256`` must be identical for
  every variant: the phenotype-free atlas never changes),
  ``route_b_families_tsv`` (Route B family proposals) and ``ranked_tsv``
  (the deterministic Stage-A shortlist). A pattern containing ``{variant_id}``
  is substituted; otherwise ``<stem>_<variant_id>.<ext>`` is derived.

Measured per variant and per case:

- **Route A**: truth-family hit and family-set Jaccard across formulations.
- **Route C**: family Recall@5/@10/@20, gene Recall@20 and retrieval family-set
  Jaccard. A retrieval that produced no rows (the Route C
  ``not_evaluated_no_indexed_token_matched_the_phenotype`` state) is
  ``not_evaluated``, never a 0% miss.
- **Route B overlay**: per-variant phenotype_consistent state from the b1
  status JSON, overlay stability across formulations, and the atlas
  fingerprint check (any mismatch across variants is an immediate error).
- **A/B four-state** (family and gene level): ``a_and_b_supported`` /
  ``a_only_supported`` / ``b_only_supported`` / ``neither_supported``.
- **Stage-A stability**: is the truth gene still in the deterministic top-20
  shortlist per formulation.

Per case the variants are classified so that harmless movement is not confused
with a real failure: ``stable`` (nothing moved), ``stable_with_rank_movement``
(truth family and truth gene maintained; only Stage-A rank moved),
``catastrophic_family_dropout`` (the truth family left the A u B proposal
surface in at least one formulation), ``stage_a_dropout`` (family maintained,
truth gene left the Stage-A top-20), ``family_stable_gene_not_evaluated``
(gene level could not be scored) and ``insufficient_variants``.

**Leak audit.** Every variant's phenotype text is checked against the answer
terms (truth family / truth gene tokens) before any scoring. A variant that
leaks the answer is reported as ``excluded_answer_derived`` in the case table
and the manifest -- never silently dropped, and never counted in any metric
denominator. The fixture rule is the same one the benchmark enforces
mechanically: answer terms are never injected into a variant to make a route
succeed.
"""

from __future__ import annotations

import json
import re
from collections import Counter

import pandas as pd

from .core import file_fingerprint, runtime_path
from .family_labels import canonical_family_key, family_set_contains
from .retrieval_route import validate_retrieval_frame

BENCHMARK_VERSION = "1.0.0"
MANIFEST_COLUMNS = ["case_id", "truth_family_id", "variants_tsv"]
OPTIONAL_COLUMNS = (
    "truth_gene_id", "base_phenotype", "route_a_families_tsv",
    "route_c_retrieval_tsv", "route_b_overlay_status_json",
    "route_b_atlas_manifest_json", "route_b_families_tsv", "ranked_tsv",
)
VARIANT_COLUMNS = ["variant_id", "variant_kind", "phenotype_text"]

# Self-defined four-state vocabulary; the existing dual-route /
# retrieval-route pair vocabularies are per-family convergence states and are
# not reused here because this benchmark compares formulation pairs, not routes.
FOUR_STATES = (
    "a_and_b_supported", "a_only_supported", "b_only_supported", "neither_supported",
)
DROPOUT_CLASSES = (
    "stable", "stable_with_rank_movement", "catastrophic_family_dropout",
    "stage_a_dropout", "family_stable_gene_not_evaluated", "insufficient_variants",
)
RECALL_LEVELS = (5, 10, 20)
STAGE_A_TOP = 20
FAMILY_COLUMNS = ("gene_family", "family_id")
OVERLAY_CONSISTENT = "phenotype_consistent"
OVERLAY_NOT_MATCHED = "phenotype_not_matched"
OVERLAY_NOT_APPLIED = "not_applied"

_OUTPUTS = {
    "cases": "phenotype_variant_cases.tsv",
    "summary": "phenotype_variant_summary.tsv",
    "manifest": "phenotype_variant_manifest.json",
    "report": "phenotype_variant_report.md",
}


def _answer_terms(record: dict[str, object]) -> set[str]:
    """Answer tokens from the truth columns (blind_benchmark._answer_terms equivalent)."""
    terms: set[str] = set()
    for column in ("truth_gene_id", "truth_family_id", "truth_gene_name", "answer_locus"):
        for token in re.split(r"[|;,\s]+", str(record.get(column, "")).casefold()):
            token = token.strip()
            if len(token) >= 2:
                terms.add(token)
    return terms


def _leaked_terms(terms: set[str], text: str) -> list[str]:
    """Answer terms present in the phenotype text (word-boundary aware)."""
    haystack = str(text or "").casefold()
    leaked: list[str] = []
    for term in sorted(terms):
        if re.search(rf"(?<![a-z0-9]){re.escape(term)}(?![a-z0-9])", haystack):
            leaked.append(term)
    return leaked


def _four_state(in_a: bool, in_b: bool) -> str:
    if in_a and in_b:
        return "a_and_b_supported"
    if in_a:
        return "a_only_supported"
    if in_b:
        return "b_only_supported"
    return "neither_supported"


def _family_values(frame: pd.DataFrame) -> set[str]:
    """Distinct non-empty family values (gene_family preferred, family_id fallback)."""
    column = next((column for column in FAMILY_COLUMNS if column in frame.columns), None)
    if column is None:
        raise ValueError("family table must contain a gene_family (or family_id) column")
    return {str(value).strip() for value in frame[column] if str(value).strip()}


def _family_keys(values: set[str]) -> set[str]:
    return {key for value in values for key in (canonical_family_key(value),) if key}


def _gene_values(frame: pd.DataFrame) -> set[str]:
    """Gene ids a proposal artifact names (member_interval_genes and/or gene_id)."""
    genes: set[str] = set()
    if "member_interval_genes" in frame.columns:
        for value in frame["member_interval_genes"]:
            genes.update(item.strip() for item in str(value).replace(",", ";").split(";") if item.strip())
    if "gene_id" in frame.columns:
        genes.update(str(value).strip() for value in frame["gene_id"] if str(value).strip())
    return genes


def _resolve_input(record: dict[str, object], column: str, config_dir, variant_id: str) -> object:
    """Resolve one per-variant artifact path; None when the column is empty."""
    value = str(record.get(column, "")).strip()
    if not value:
        return None
    if "{variant_id}" in value:
        text = value.replace("{variant_id}", variant_id)
    else:
        path = runtime_path(value)
        text = str(path.parent / f"{path.stem}_{variant_id}{path.suffix}")
    path = runtime_path(text)
    if not path.is_absolute():
        path = config_dir / path
    if not path.exists():
        raise ValueError(f"phenotype variant benchmark input does not exist: {path} (case {record['case_id']}, variant {variant_id})")
    return path


def _read_tsv(path: object) -> pd.DataFrame:
    return pd.read_csv(runtime_path(path), sep="\t", dtype=str).fillna("")


def _resolve_atlas_manifest(record: dict[str, object], config_dir, variant_id: str) -> object:
    """Resolve the route B atlas manifest; a per-variant file, or one shared
    case-level file when no per-variant copy exists (the atlas is phenotype-free
    and must be the same object for every formulation)."""
    value = str(record.get("route_b_atlas_manifest_json", "")).strip()
    if not value:
        return None
    if "{variant_id}" in value:
        text = value.replace("{variant_id}", variant_id)
    else:
        path = runtime_path(value)
        text = str(path.parent / f"{path.stem}_{variant_id}{path.suffix}")
    candidate = runtime_path(text)
    if not candidate.is_absolute():
        candidate = config_dir / candidate
    if candidate.exists():
        return candidate
    shared = runtime_path(value)
    if not shared.is_absolute():
        shared = config_dir / shared
    if shared.exists():
        return shared
    raise ValueError(f"phenotype variant benchmark atlas manifest does not exist: {candidate} (case {record['case_id']}, variant {variant_id})")


def _c_family_rank(frame: pd.DataFrame, truth_family: str) -> int | None:
    best: int | None = None
    for row in frame.to_dict("records"):
        for value in (str(row.get("family_id", "")).strip(), str(row.get("orthogroup_id", "")).strip()):
            if value and family_set_contains({value}, truth_family):
                rank = int(row["rank"])
                best = rank if best is None else min(best, rank)
    return best


def _top_families(frame: pd.DataFrame, k: int) -> set[str]:
    """Distinct families/orthogroups ranked by their best row rank (route_c_benchmark equivalent)."""
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


def _retrieval_family_keys(frame: pd.DataFrame) -> set[str]:
    keys: set[str] = set()
    for row in frame.to_dict("records"):
        for value in (str(row.get("family_id", "")).strip(), str(row.get("orthogroup_id", "")).strip()):
            if value:
                key = canonical_family_key(value)
                if key:
                    keys.add(key)
    return keys


def _stage_a_rank(frame: pd.DataFrame, truth_gene: str) -> tuple[bool, int | None]:
    """Whether the truth gene is inside the deterministic top-20 shortlist, and its 1-based row rank.

    Order comes from ``deterministic_review_order`` (the actual Stage-A artifact
    column), else ``rank`` / ``dry_rank``, else the file's own row order.
    """
    if "gene_id" not in frame.columns:
        raise ValueError("ranked table must contain gene_id")
    ordered = frame.copy()
    order_column = next((column for column in ("deterministic_review_order", "rank", "dry_rank", "program_rank", "final_rank")
                         if column in ordered.columns), None)
    if order_column:
        ordered["_order"] = pd.to_numeric(ordered[order_column], errors="coerce").fillna(1 << 30)
        ordered = ordered.sort_values(["_order", "gene_id"], kind="stable")
    truth_key = truth_gene.casefold()
    names = [str(value).strip().casefold() for value in ordered["gene_id"] if str(value).strip()]
    position = next((index for index, name in enumerate(names, start=1) if name == truth_key), None)
    return bool(position is not None and position <= STAGE_A_TOP), position


def _mean_pairwise_jaccard(sets: list[set[str]]) -> float:
    """Mean Jaccard over all pairs. Pairs whose union is empty score 0.0 (two empty
    proposal sets are not stability; they are a shared failure)."""
    pairs = 0
    total = 0.0
    for i in range(len(sets)):
        for j in range(i + 1, len(sets)):
            left, right = sets[i], sets[j]
            union = left | right
            total += len(left & right) / len(union) if union else 0.0
            pairs += 1
    return round(total / pairs, 6) if pairs else 0.0


def _common_keys(sets: list[set[str]]) -> list[str]:
    if not sets:
        return []
    return sorted(set.intersection(*sets))


def _status_bool(value: object, default: bool = False) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value).strip().casefold()
    if not text:
        return default
    return text in {"true", "1", "yes"}


def _overlay_state(status_data: dict[str, object]) -> tuple[bool, str, bool]:
    """(phenotype_consistent, status_counts_text, overlay_evaluable) from a b1 status JSON."""
    counts = status_data.get("phenotype_overlay", {})
    if isinstance(counts, dict):
        counts = counts.get("status_counts", {}) if "status_counts" in counts else counts
    if not isinstance(counts, dict):
        counts = {}
    counts_text = ";".join(f"{key}={counts[key]}" for key in sorted(counts))
    phenotype_used = _status_bool(status_data.get("phenotype_used", False), False)
    if not phenotype_used:
        return False, counts_text, False
    consistent = int(counts.get(OVERLAY_CONSISTENT, 0)) > 0
    return consistent, counts_text, True


def _atlas_sha256(atlas_manifest_path: object, case_id: str) -> str:
    data = json.loads(runtime_path(atlas_manifest_path).read_text(encoding="utf-8"))
    if "atlas_sha256" not in data:
        raise ValueError(f"route B atlas manifest for case {case_id} is missing atlas_sha256: {atlas_manifest_path}")
    return str(data["atlas_sha256"])


def _evaluate_case(record: dict[str, object], config_dir) -> tuple[list[dict[str, object]], dict[str, object]]:
    case_id = str(record["case_id"]).strip()
    truth_family = str(record["truth_family_id"]).strip()
    truth_gene = str(record.get("truth_gene_id", "")).strip()
    base_phenotype = str(record.get("base_phenotype", "")).strip()
    answer_terms = _answer_terms(record)

    variants_path = runtime_path(str(record["variants_tsv"]).strip())
    if not variants_path.is_absolute():
        variants_path = config_dir / variants_path
    if not variants_path.exists():
        raise ValueError(f"variants table does not exist: {variants_path} (case {case_id})")
    variants = pd.read_csv(variants_path, sep="\t", dtype=str).fillna("")
    missing = [column for column in VARIANT_COLUMNS if column not in variants.columns]
    if missing:
        raise ValueError(f"variants table for case {case_id} is missing columns: {', '.join(missing)}")
    if variants["variant_id"].astype(str).str.strip().eq("").any():
        raise ValueError(f"variants table for case {case_id} requires a variant_id per row")
    if variants["variant_id"].astype(str).str.strip().duplicated().any():
        raise ValueError(f"variants table for case {case_id} has duplicated variant_id values")
    if variants["phenotype_text"].astype(str).str.strip().eq("").any():
        raise ValueError(f"variants table for case {case_id} requires a phenotype_text per row")

    rows: list[dict[str, object]] = []
    atlas_shas: set[str] = set()
    a_sets: list[set[str]] = []
    c_sets: list[set[str]] = []
    eligible_ids: list[str] = []

    for variant in variants.to_dict("records"):
        variant_id = str(variant["variant_id"]).strip()
        variant_kind = str(variant.get("variant_kind", "")).strip()
        phenotype_text = str(variant.get("phenotype_text", "")).strip()

        leaked = _leaked_terms(answer_terms, phenotype_text)
        row: dict[str, object] = {
            "case_id": case_id, "variant_id": variant_id, "variant_kind": variant_kind,
            "phenotype_text": phenotype_text,
            "leak_status": "excluded_answer_derived" if leaked else "eligible",
            "leaked_terms": ";".join(leaked),
        }
        if leaked:
            rows.append(row)
            continue

        a_path = _resolve_input(record, "route_a_families_tsv", config_dir, variant_id)
        c_path = _resolve_input(record, "route_c_retrieval_tsv", config_dir, variant_id)
        b_status_path = _resolve_input(record, "route_b_overlay_status_json", config_dir, variant_id)
        atlas_manifest_path = _resolve_atlas_manifest(record, config_dir, variant_id)
        b_families_path = _resolve_input(record, "route_b_families_tsv", config_dir, variant_id)
        ranked_path = _resolve_input(record, "ranked_tsv", config_dir, variant_id)

        a_families: set[str] = set()
        if a_path is not None:
            a_families = _family_values(_read_tsv(a_path))
        a_truth_hit = family_set_contains(a_families, truth_family) if a_path is not None else None

        c_rows = 0
        c_family_rank: int | None = None
        c_family_hits = {k: None for k in RECALL_LEVELS}
        c_gene_hit: bool | None = None
        if c_path is not None:
            raw_retrieval = _read_tsv(c_path)
            c_rows = int(len(raw_retrieval))
            # An empty retrieval is the Route C
            # not_evaluated_no_indexed_token_matched_the_phenotype state: there
            # is no ranking to score, so the variant is not_evaluated, never a
            # 0% miss. The contract validator only applies to non-empty frames.
            if c_rows:
                retrieval = validate_retrieval_frame(raw_retrieval, f"retrieval for case {case_id} variant {variant_id}")
                if "rank" in retrieval.columns:
                    retrieval["rank"] = pd.to_numeric(retrieval["rank"], errors="coerce").astype(int)
                top = {k: _top_families(retrieval, k) for k in RECALL_LEVELS}
                c_family_hits = {k: family_set_contains(top[k], truth_family) for k in RECALL_LEVELS}
                c_family_rank = _c_family_rank(retrieval, truth_family)
                if truth_gene:
                    top_genes = {str(value).strip() for value in
                                 retrieval.sort_values("rank", kind="stable")["reference_gene_id"].head(20)}
                    c_gene_hit = truth_gene in top_genes

        b_overlay_consistent: bool | None = None
        b_overlay_counts = ""
        if b_status_path is not None:
            status_data = json.loads(runtime_path(b_status_path).read_text(encoding="utf-8"))
            if str(status_data.get("route", "")) != "interval_first":
                raise ValueError(f"Route B status has an unexpected route: {b_status_path}")
            b_overlay_consistent, b_overlay_counts, overlay_evaluable = _overlay_state(status_data)
            if not overlay_evaluable:
                b_overlay_consistent = None
        if atlas_manifest_path is not None:
            atlas_shas.add(_atlas_sha256(atlas_manifest_path, case_id))

        b_families: set[str] = set()
        if b_families_path is not None:
            b_families = _family_values(_read_tsv(b_families_path))
        b_truth_hit = family_set_contains(b_families, truth_family) if b_families_path is not None else None

        a_genes: set[str] = set()
        b_genes: set[str] = set()
        if a_path is not None:
            a_genes = _gene_values(_read_tsv(a_path))
        if b_families_path is not None:
            b_genes = _gene_values(_read_tsv(b_families_path))
        gene_four_evaluable = bool(truth_gene) and a_path is not None and b_families_path is not None
        in_a_gene = truth_gene.casefold() in {value.casefold() for value in a_genes} if gene_four_evaluable else None
        in_b_gene = truth_gene.casefold() in {value.casefold() for value in b_genes} if gene_four_evaluable else None

        stage_a_present = ranked_path is not None
        stage_a_rank: int | None = None
        stage_a_top20_hit: bool | None = None
        if ranked_path is not None and truth_gene:
            stage_a_top20_hit, stage_a_rank = _stage_a_rank(_read_tsv(ranked_path), truth_gene)

        row.update({
            "a_families_count": len(a_families) if a_path is not None else "",
            "a_truth_family_hit": a_truth_hit,
            "c_retrieval_rows": c_rows if c_path is not None else "",
            "c_evaluable": bool(c_rows) if c_rows else "",
            "c_family_rank": c_family_rank,
            **{f"family_recall_{k}": c_family_hits[k] for k in RECALL_LEVELS},
            "gene_recall_20": c_gene_hit,
            "b_overlay_consistent": b_overlay_consistent,
            "b_overlay_status_counts": b_overlay_counts,
            "b_truth_family_hit": b_truth_hit,
            "four_state": _four_state(bool(a_truth_hit), bool(b_truth_hit)) if a_path is not None and b_families_path is not None else "",
            "four_state_gene": (_four_state(bool(in_a_gene), bool(in_b_gene)) if gene_four_evaluable else ""),
            "stage_a_present": bool(stage_a_present) if (truth_gene and ranked_path is not None) else "",
            "stage_a_rank": stage_a_rank,
            "stage_a_top20_hit": stage_a_top20_hit,
        })
        if a_path is not None:
            a_sets.append(_family_keys(a_families))
        if c_path is not None and c_rows:
            c_sets.append(_retrieval_family_keys(retrieval))
        eligible_ids.append(variant_id)
        rows.append(row)

    eligible = [row for row in rows if row["leak_status"] == "eligible"]
    def first_value(column: str) -> object:
        return eligible[0].get(column, "") if eligible else ""

    four_states = sorted({str(row["four_state"]) for row in eligible if str(row["four_state"])})
    gene_four_states = sorted({str(row["four_state_gene"]) for row in eligible if str(row["four_state_gene"])})
    a_original = first_value("a_truth_family_hit")
    b_original = first_value("b_truth_family_hit")
    a_b_union_hit_values = [
        bool(row["a_truth_family_hit"] or row["b_truth_family_hit"])
        for row in eligible
        if row["a_truth_family_hit"] is not None or row["b_truth_family_hit"] is not None
    ]
    family_union_sustained = bool(a_b_union_hit_values) and all(a_b_union_hit_values) if a_b_union_hit_values else None
    family_route_lost: list[str] = []
    for row in eligible:
        if row["a_truth_family_hit"] is True and row["b_truth_family_hit"] is True:
            continue
        if row["a_truth_family_hit"] is False or row["b_truth_family_hit"] is False:
            if not (row["a_truth_family_hit"] or row["b_truth_family_hit"]):
                family_route_lost.append(f"{row['variant_id']} (both A and B)")
            elif row["a_truth_family_hit"] is False:
                family_route_lost.append(f"{row['variant_id']} (A)")
            elif row["b_truth_family_hit"] is False:
                family_route_lost.append(f"{row['variant_id']} (B)")

    stage_a_ranks = [int(row["stage_a_rank"]) for row in eligible if row.get("stage_a_rank") is not None]
    stage_a_hits = [row["stage_a_top20_hit"] for row in eligible if row["stage_a_top20_hit"] is not None]
    gene_sustained: bool | None = bool(stage_a_hits) and all(stage_a_hits) if stage_a_hits else None
    gene_evaluable = bool(truth_gene) and any(row["stage_a_present"] is True for row in eligible)

    four_state_regressed = bool(
        any(str(row["four_state"]) == "a_and_b_supported" for row in eligible)
        and any(str(row["four_state"]) not in ("", "a_and_b_supported") for row in eligible)
    )
    gene_four_regressed = bool(
        any(str(row["four_state_gene"]) == "a_and_b_supported" for row in eligible)
        and any(str(row["four_state_gene"]) not in ("", "a_and_b_supported") for row in eligible)
    )

    dropout_class = ""
    dropout_reason = ""
    if len(eligible) < 2:
        dropout_class = "insufficient_variants"
        dropout_reason = f"only {len(eligible)} eligible variant(s); dropout classification needs at least 2"
    elif family_union_sustained is False:
        dropout_class = "catastrophic_family_dropout"
        dropout_reason = "the truth family left the A u B proposal surface in: " + "; ".join(family_route_lost)
    elif family_union_sustained is None:
        dropout_class = "insufficient_variants"
        dropout_reason = "no variant carried both a Route A and a Route B artifact"
    elif not gene_evaluable:
        dropout_class = "family_stable_gene_not_evaluated"
        dropout_reason = "the truth gene was not scoreable in Stage-A (no truth_gene_id and/or no ranked_tsv)"
    elif gene_sustained is False:
        dropout_class = "stage_a_dropout"
        dropout_reason = "the truth family stayed in A u B but the truth gene left the Stage-A top-20 in: " + "; ".join(
            str(row["variant_id"]) for row in eligible if row["stage_a_top20_hit"] is False
        )
    elif stage_a_ranks and max(stage_a_ranks) != min(stage_a_ranks):
        dropout_class = "stable_with_rank_movement"
        dropout_reason = f"truth family and truth gene maintained in every variant; Stage-A rank moved between {min(stage_a_ranks)} and {max(stage_a_ranks)}"
    else:
        dropout_class = "stable"
        dropout_reason = "truth family and truth gene maintained with an unchanged Stage-A rank across all variants"

    overlay_consistents = [row["b_overlay_consistent"] for row in eligible if row["b_overlay_consistent"] is not None]
    b_stability = round(sum(1 for value in overlay_consistents if value) / len(overlay_consistents), 6) if overlay_consistents else None
    b_original = first_value("b_overlay_consistent")

    summary = {
        "case_id": case_id, "base_phenotype": base_phenotype,
        "truth_family_id": truth_family, "truth_gene_id": truth_gene,
        "variant_count": len(rows), "eligible_count": len(eligible),
        "excluded_answer_derived": len(rows) - len(eligible),
        "a_truth_family_hit_original": a_original,
        "b_truth_family_hit_original": b_original,
        "a_common_families": ";".join(_common_keys(a_sets)),
        # A missing artifact is not a measured 0.0 overlap: leave the cell empty
        # so the corpus aggregate and readers can tell not_evaluated apart from
        # a genuinely empty proposal pair.
        "a_mean_pairwise_jaccard": _mean_pairwise_jaccard(a_sets) if a_sets else "",
        "c_mean_pairwise_jaccard": _mean_pairwise_jaccard(c_sets) if c_sets else "",
        "b_overlay_consistent_original": b_original,
        "b_overlay_stability": b_stability,
        "atlas_sha256_values": sorted(atlas_shas),
        "atlas_consistent": bool(atlas_shas and len(atlas_shas) == 1) if atlas_shas else None,
        "four_state_original": first_value("four_state"),
        "four_states_observed": ";".join(four_states),
        "four_state_gene_original": first_value("four_state_gene"),
        "four_state_gene_states_observed": ";".join(gene_four_states),
        "four_state_regressed": four_state_regressed,
        "gene_four_state_regressed": gene_four_regressed,
        "stage_a_original_rank": first_value("stage_a_rank"),
        "stage_a_rank_min": min(stage_a_ranks) if stage_a_ranks else "",
        "stage_a_rank_max": max(stage_a_ranks) if stage_a_ranks else "",
        "stage_a_top20_hit_all": gene_sustained,
        "family_union_sustained": family_union_sustained,
        "family_route_lost": "; ".join(family_route_lost),
        "gene_sustained": gene_sustained,
        "gene_evaluable": gene_evaluable,
        "dropout_class": dropout_class,
        "dropout_reason": dropout_reason,
    }
    return rows, summary


def phenotype_variant_benchmark(manifest_arg: str, out_dir_arg: str, full_hash: bool = False) -> dict[str, object]:
    """Score frozen per-formulation artifacts for phenotype-wording robustness."""
    manifest_path = runtime_path(manifest_arg)
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str).fillna("")
    missing = [column for column in MANIFEST_COLUMNS if column not in manifest.columns]
    if missing:
        raise ValueError(f"phenotype variant benchmark manifest missing columns: {', '.join(missing)}")
    if manifest.empty:
        raise ValueError("phenotype variant benchmark manifest is empty")
    if manifest["case_id"].astype(str).str.strip().duplicated().any():
        raise ValueError("phenotype variant benchmark case_id must be unique")
    if manifest["case_id"].astype(str).str.strip().eq("").any():
        raise ValueError("phenotype variant benchmark requires a case_id per case")
    if manifest["truth_family_id"].astype(str).str.strip().eq("").any():
        raise ValueError("phenotype variant benchmark requires a truth_family_id per case")
    if manifest["variants_tsv"].astype(str).str.strip().eq("").any():
        raise ValueError("phenotype variant benchmark requires a variants_tsv per case")

    config_dir = manifest_path.parent
    case_rows: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []
    fingerprint_targets: list = [manifest_path]
    for record in manifest.to_dict("records"):
        rows, summary = _evaluate_case(record, config_dir)
        case_rows.extend(rows)
        summaries.append(summary)
        variants_path = runtime_path(str(record["variants_tsv"]).strip())
        if not variants_path.is_absolute():
            variants_path = config_dir / variants_path
        fingerprint_targets.append(variants_path)

    cases = pd.DataFrame(case_rows).fillna("").sort_values(["case_id", "variant_id"], kind="stable")
    cases_path = out / _OUTPUTS["cases"]
    cases.to_csv(cases_path, sep="\t", index=False, lineterminator="\n")
    summary = pd.DataFrame(summaries).fillna("").sort_values("case_id", kind="stable")
    summary_path = out / _OUTPUTS["summary"]
    summary.to_csv(summary_path, sep="\t", index=False, lineterminator="\n")

    n = len(summary)
    eligible_total = int(cases["leak_status"].eq("eligible").sum())
    excluded_rows = cases[cases["leak_status"].ne("eligible")]
    excluded_variants = sorted(f"{row['case_id']}:{row['variant_id']}" for row in excluded_rows.to_dict("records"))
    leak_counts = {key: count for key, count in Counter(
        str(row["leaked_terms"]) for row in excluded_rows.to_dict("records")
    ).most_common()}

    def rate(column: str) -> float:
        values = [row[column] for row in summaries if row.get(column) is not None]
        return round(sum(1 for value in values if value is True) / len(values), 6) if values else None

    def per_variant_stats(column: str) -> dict[str, object]:
        values = cases[column].astype(str).str.strip()
        evaluated = int(values.ne("").sum())
        hits = int(values.eq("True").sum())
        return {"evaluated": evaluated, "hits": hits,
                "recall": round(hits / evaluated, 6) if evaluated else None,
                "not_evaluated": int(len(cases)) - evaluated}

    dropout_counts = {value: int((summary["dropout_class"] == value).sum()) for value in DROPOUT_CLASSES}
    classified = summary[summary["dropout_class"].isin(("stable", "stable_with_rank_movement", "catastrophic_family_dropout", "stage_a_dropout"))]
    a_jaccards = [float(row["a_mean_pairwise_jaccard"]) for row in summaries if row["a_mean_pairwise_jaccard"] and isinstance(row["a_mean_pairwise_jaccard"], float)]
    c_jaccards = [float(row["c_mean_pairwise_jaccard"]) for row in summaries if row["c_mean_pairwise_jaccard"] and isinstance(row["c_mean_pairwise_jaccard"], float)]
    overlays = [float(row["b_overlay_stability"]) for row in summaries if row["b_overlay_stability"] is not None]
    atlas_checks = [row["atlas_consistent"] for row in summaries if row["atlas_consistent"] is not None]
    if atlas_checks and not all(atlas_checks):
        offenders = sorted(str(row["case_id"]) for row in summaries if row["atlas_consistent"] is False)
        raise ValueError("Route B atlas fingerprint differs across variants; the phenotype-free atlas must never change: " + ", ".join(offenders))

    manifest_out = {
        "schema_version": "1.0",
        "benchmark_version": BENCHMARK_VERSION,
        "benchmark_type": "phenotype wording / description-granularity robustness (frozen per-formulation artifacts, issue #235)",
        "case_count": n,
        "variant_total": int(len(cases)),
        "eligible_variant_total": eligible_total,
        "excluded_answer_derived_count": int(len(excluded_rows)),
        "excluded_answer_derived_variants": excluded_variants,
        "excluded_leaked_terms": leak_counts,
        "freeze": {
            "declaration": "variants.tsv is the frozen formulation freeze; the benchmark reads frozen artifacts only and never generates or optimizes a variant",
            "variants_tsv_fingerprints": {str(path): file_fingerprint(path, full_hash) for path in fingerprint_targets[1:]},
            "mode": "frozen_artifacts_only",
        },
        "leak_audit": {
            "declaration": "each variant's phenotype_text is checked against the answer terms (truth family / truth gene tokens) before scoring; a leaking variant is reported as excluded_answer_derived and is never silently dropped",
            "rule": "answer terms are never injected into a variant to make a route succeed; variants are formulated before the truth is used for anything but scoring",
            "excluded_count": int(len(excluded_rows)),
        },
        "route_a": {
            "a_truth_family_hit_original_rate": rate("a_truth_family_hit_original"),
            "a_mean_pairwise_family_jaccard": round(sum(a_jaccards) / len(a_jaccards), 6) if a_jaccards else None,
            "note": "mean over cases of the per-case mean pairwise Jaccard between formulations' Route A family sets (canonical family keys)",
        },
        "route_c": {
            "family_recall": {f"family_recall_{k}": per_variant_stats(f"family_recall_{k}") for k in RECALL_LEVELS},
            "gene_recall_20": per_variant_stats("gene_recall_20"),
            "retrieval_family_jaccard": round(sum(c_jaccards) / len(c_jaccards), 6) if c_jaccards else None,
            "not_evaluated_rule": "a variant whose retrieval produced no rows (Route C no-indexed-token state) or whose artifact is absent is not_evaluated, never a 0% miss",
        },
        "route_b": {
            "overlay_stability_mean": round(sum(overlays) / len(overlays), 6) if overlays else None,
            "overlay_stability_rule": "fraction of eligible variants per case whose b1 overlay is phenotype_consistent; the phenotype-free atlas is never re-built",
            "atlas_fingerprint_checked": bool(atlas_checks),
            "atlas_fingerprint_consistent": bool(atlas_checks and all(atlas_checks)),
        },
        "four_state_vocabulary": list(FOUR_STATES),
        "four_state": {
            "definition": "family (or gene) level membership: in_a = truth in Route A proposals, in_b = truth in Route B proposals",
            "per_case_four_states_observed": {str(row["case_id"]): str(row["four_states_observed"]) for row in summaries},
        },
        "stage_a": {
            "rule": f"truth gene inside the deterministic top-{STAGE_A_TOP} shortlist per variant; absent truth/input is not_evaluated",
            "stage_a_top20_hit_all_rate": rate("stage_a_top20_hit_all"),
        },
        "dropout_classification": {
            "vocabulary": list(DROPOUT_CLASSES),
            "distinction": "stable / stable_with_rank_movement tolerate harmless movement; catastrophic_family_dropout is the truth family leaving the A u B proposal surface; stage_a_dropout is the truth gene leaving the Stage-A top-20 only",
            "counts": dropout_counts,
            "classified_cases": int(len(classified)),
            "stable_rate": round(int((summary["dropout_class"] == "stable").sum()) / len(classified), 6) if len(classified) else None,
            "stable_with_rank_movement_rate": round(int((summary["dropout_class"] == "stable_with_rank_movement").sum()) / len(classified), 6) if len(classified) else None,
            "catastrophic_family_dropout_rate": round(int((summary["dropout_class"] == "catastrophic_family_dropout").sum()) / len(classified), 6) if len(classified) else None,
            "stage_a_dropout_rate": round(int((summary["dropout_class"] == "stage_a_dropout").sum()) / len(classified), 6) if len(classified) else None,
            "cases": {str(row["case_id"]): str(row["dropout_class"]) for row in summaries},
        },
        "policy": [
            "no prediction logic is changed; this benchmark measures the frozen artifacts of the existing routes",
            "the phenotype-free atlas is invariant across formulations; a differing atlas fingerprint aborts the run",
            "Route B overlay is read from the b1 status JSONs; it never changes tier, rank, or score",
            "answer terms are never injected into a phenotype variant; leaking variants are excluded and reported",
            "formulations are frozen before truth scoring; the benchmark reads only frozen artifacts",
        ],
        "manifest_fingerprint": file_fingerprint(manifest_path, full_hash),
        "input_fingerprints": {str(path): file_fingerprint(path, full_hash) for path in fingerprint_targets[1:]},
        "fingerprint_mode": "strict_full_hash" if full_hash else "size_tiered",
        "outputs": {name: str(out / path) for name, path in _OUTPUTS.items()},
        "scientific_limits": [
            "artifacts are external frozen lane outputs; the harness verifies their column contract and fingerprints, never their generation process",
            "Jaccard pairs whose union is empty score 0.0 (two empty proposal sets are a shared failure, not stability)",
            "dropout classification is per case over eligible variants; variants excluded by the leak audit are reported but never counted",
            "measurement only: a stable number says the existing routes are wording-robust, not that they are correct",
        ],
    }
    manifest_json = out / _OUTPUTS["manifest"]
    manifest_json.write_text(json.dumps(manifest_out, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write_report(cases, summary, manifest_out, out)
    return {
        "cases": str(cases_path), "summary": str(summary_path), "manifest": str(manifest_json),
        "report": str(out / _OUTPUTS["report"]), "case_count": n,
        "excluded_answer_derived": int(len(excluded_rows)),
        "dropout_classification": manifest_out["dropout_classification"],
        "atlas_fingerprint_consistent": manifest_out["route_b"]["atlas_fingerprint_consistent"],
    }


def _rate(value: object) -> str:
    return "not_evaluated" if value is None else f"{float(value):.1%}"


def _write_report(cases: pd.DataFrame, summary: pd.DataFrame, manifest: dict[str, object], out) -> None:
    dropout = manifest["dropout_classification"]
    route_c = manifest["route_c"]
    lines = [
        "# Phenotype wording robustness benchmark (issue #235)",
        "",
        "The same truth case is scored across frozen phenotype formulations",
        "(original / short / synonymized / mechanism terms removed / irrelevant",
        "detail added). The benchmark reads frozen per-formulation lane artifacts;",
        "it never re-runs a lane and never changes a score.",
        "",
        f"- cases: {manifest['case_count']}; variants: {manifest['variant_total']}; "
        f"eligible variants: {manifest['eligible_variant_total']}",
        f"- leak audit: {manifest['excluded_answer_derived_count']} variant(s) excluded as answer-derived"
        + (f": {', '.join(manifest['excluded_answer_derived_variants'])}" if manifest["excluded_answer_derived_variants"] else ""),
        f"- Route A truth-family hit (original variants): {_rate(manifest['route_a']['a_truth_family_hit_original_rate'])}; "
        f"mean pairwise family Jaccard across formulations: {_rate(manifest['route_a']['a_mean_pairwise_family_jaccard'])}",
        f"- Route C family Recall@20: {_rate(route_c['family_recall']['family_recall_20']['recall'])}; "
        f"retrieval family Jaccard: {_rate(route_c['retrieval_family_jaccard'])}",
        f"- Route B overlay stability: {_rate(manifest['route_b']['overlay_stability_mean'])}; "
        f"atlas fingerprint consistent across variants: {manifest['route_b']['atlas_fingerprint_consistent']}",
        f"- Stage-A top-{STAGE_A_TOP} hit in every variant: {_rate(manifest['stage_a']['stage_a_top20_hit_all_rate'])}",
        f"- stable: {_rate(dropout['stable_rate'])}; stable with rank movement: {_rate(dropout['stable_with_rank_movement_rate'])}; "
        f"catastrophic family dropout: {_rate(dropout['catastrophic_family_dropout_rate'])}; stage-A dropout: {_rate(dropout['stage_a_dropout_rate'])}",
        "",
        "## Per case",
        "",
        "| case | truth family | eligible | four-state (orig) | A jaccard | C jaccard | B overlay | stage-A rank | dropout class |",
        "|---|---|---:|---|---|---|---:|---:|---|",
    ]
    for row in summary.to_dict("records"):
        lines.append(f"| {row['case_id']} | {row['truth_family_id']} | {row['eligible_count']} | "
                     f"{row['four_state_original']} | {row['a_mean_pairwise_jaccard']} | {row['c_mean_pairwise_jaccard']} | "
                     f"{row['b_overlay_stability']} | {row['stage_a_rank_min']}-{row['stage_a_rank_max']} | {row['dropout_class']} |")
    lines += [
        "",
        "Distinction: harmless rank movement (stable_with_rank_movement) is not a failure;",
        "catastrophic_family_dropout is the truth family leaving the A u B proposal surface,",
        "stage_a_dropout is the truth gene leaving the Stage-A top-20 only. Excluded answer-derived",
        "variants are reported, never silently dropped.",
    ]
    (out / _OUTPUTS["report"]).write_text("\n".join(lines) + "\n", encoding="utf-8")
