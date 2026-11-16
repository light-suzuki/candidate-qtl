"""Phenotype-first blind benchmark (the autonomous floor).

The map-based benchmark (``benchmark-reference`` / ``discovery-benchmark``)
measures recall *after* a human hands the pipeline a legitimate cross-species
pathway homolog to search the interval with. That is a real capability, but it
is not the ability to decide **which pathway / family to investigate from the
phenotype alone**. This benchmark measures exactly that, blind:

Per case the only admissible inputs are the QTL interval, the phenotype in
natural text (plus direction and inheritance), the reference interval genes with
their *generic* native annotation, and a single **shared, pre-fixed** phenotype
-> pathway/family route that is identical for every case. Barred: a human-chosen
homolog / query, the causal gene name / locus / answer-derived sequence, terms
that were tuned per case to describe the answer, and any weight tuning after
seeing results.

An eligibility gate rejects, *before* ranking, any case whose phenotype-derived
terms leak the answer (the causal gene id / name / locus), and flags a route
keyword that coincides with the answer (a case-tuned term). The deterministic,
no-LLM route is the mandatory baseline; an optional per-case LLM ranking over the
*same* frozen interval-gene set is scored separately and never merged in.

Evaluation is layered so a miss is attributable:

1. phenotype -> pathway/family hypothesis recovery,
2. pathway/family -> interval-gene recovery,
3. final gene rank (tie-aware Top-N / MRR / wrong-sole-top / abstention).

The initial autonomous floor is recorded as-is (including 0%). Improvement must
come from new, provenance-recorded signal; this module changes no weight and
holds no hidden rank logic, so a better number cannot come from fitting.
"""

from __future__ import annotations

import json
import re

import pandas as pd

from .benchmark_metrics import aggregate_resolution_metrics, resolution_metrics
from .core import file_fingerprint, runtime_path

BENCHMARK_VERSION = "1.0.0"
MANIFEST_COLUMNS = ["case_id", "phenotype_text", "interval", "interval_genes", "answer_gene_id"]
ROUTE_COLUMNS = ["phenotype_token", "pathway_family", "annotation_keywords"]
# Deliberately generic: broad phenotype words carry no pathway on their own and
# must not, alone, drive a hypothesis. Kept small and shared, never per case.
_STOPWORDS = frozenset(
    "the a an and or of to in on with without at is are be increased decreased reduced "
    "increase decrease reduction higher lower more less altered abnormal plant plants "
    "phenotype trait leaf leaves stem seed seeds size number shape color colour".split()
)
_TOKEN = re.compile(r"[a-z][a-z0-9]{2,}")


def phenotype_tokens(text: str) -> list[str]:
    """Generic phenotype tokens, stop-worded and de-duplicated.

    Public because the answerless prediction lane (``audit-predict``, issue #213)
    must derive its route through exactly the same tokens as this benchmark; a
    second implementation would mean the frozen prediction and the scored
    prediction were not the same object.
    """
    tokens = [token for token in _TOKEN.findall(str(text).casefold()) if token not in _STOPWORDS]
    # stable de-dup preserving first occurrence
    return list(dict.fromkeys(tokens))


def derive_route(tokens: list[str], route: pd.DataFrame) -> tuple[set[str], set[str]]:
    """Map phenotype tokens onto the shared route's pathways / annotation keywords.

    The route is a single shared file applied identically to every case, so a
    keyword cannot be tuned per case by construction.
    """
    pathways: set[str] = set()
    keywords: set[str] = set()
    for row in route.to_dict("records"):
        token = str(row["phenotype_token"]).strip().casefold()
        if token and token in tokens:
            pathways.add(str(row["pathway_family"]).strip())
            for keyword in re.split(r"[|;,]+", str(row["annotation_keywords"])):
                keyword = keyword.strip()
                if keyword:
                    keywords.add(keyword)
    return pathways, keywords


def _answer_terms(record: dict[str, object]) -> set[str]:
    terms: set[str] = set()
    for column in (
        "answer_gene_id", "answer_gene_name", "answer_locus",
        "known_answer_aliases", "causal_publication_ids", "truth_accession", "truth_gene_id", "truth_coordinate",
    ):
        raw = str(record.get(column, "")).casefold()
        for token in re.split(r"[|;,\s]+", raw):
            token = token.strip()
            if len(token) >= 2:
                terms.add(token)
        for token in _TOKEN.findall(raw):
            if len(token) >= 2 and token not in _STOPWORDS:
                terms.add(token)
    return terms


def _keyword_hit(annotation: str, keyword: str) -> bool:
    return re.search(rf"(?<![a-z0-9]){re.escape(keyword.casefold())}(?![a-z0-9])", annotation.casefold()) is not None


def _annotation_text(interval_genes: pd.DataFrame) -> pd.Series:
    """Resolve each row's functional-annotation text for keyword scoring.

    Column-name agnostic (issue #317): ``annotation``, the historical
    blind-benchmark fixture contract, takes priority; otherwise ``description``
    and ``ontology_terms``, the extract candidate-table contract
    (``reference_extract.write_candidates``), are joined with a single space so
    whole-word keyword matching is preserved. Empty / NaN cells are skipped; a
    row with none of these columns degrades to the historical empty string.
    """
    if "annotation" in interval_genes.columns:
        return interval_genes["annotation"].astype(str)
    parts: list[pd.Series] = []
    for column in ("description", "ontology_terms"):
        if column in interval_genes.columns:
            raw = interval_genes[column]
            parts.append(raw.astype(str).str.strip().mask(raw.isna(), ""))
    if not parts:
        return pd.Series([""] * len(interval_genes), index=interval_genes.index)
    text = parts[0]
    for part in parts[1:]:
        text = (text + " " + part).str.strip()
    return text


def rank_interval(interval_genes: pd.DataFrame, keywords: set[str]) -> pd.DataFrame:
    """Score each interval gene by how many pathway keywords its generic annotation matches.

    The annotation text is resolved column-name agnostically by
    :func:`_annotation_text`: ``annotation`` when present, otherwise
    ``description`` + ``ontology_terms`` (the extract output contract, issue
    #317). Public for the same reason as :func:`phenotype_tokens`: the
    answerless prediction lane freezes this exact ranking (issue #213).
    """
    ranked = interval_genes.copy()
    annotations = _annotation_text(ranked)
    scores = [float(sum(1 for keyword in keywords if _keyword_hit(annotation, keyword))) for annotation in annotations]
    ranked["dry_score"] = scores
    max_score = max(scores) if scores else 0.0
    # Confident call only when a single gene uniquely leads with a positive score;
    # otherwise the route abstains (the honest floor when nothing matches).
    sole_leader = max_score > 0 and scores.count(max_score) == 1
    ranked["dry_group"] = [
        "top" if (sole_leader and score == max_score) else ("reserve" if score > 0 else "evidence_insufficient")
        for score in scores
    ]
    return ranked[["gene_id", "dry_score", "dry_group"]]


def _evaluate_case(record: dict[str, object], route: pd.DataFrame, config_dir) -> dict[str, object]:
    case_id = str(record["case_id"]).strip()
    tokens = phenotype_tokens(record["phenotype_text"])
    answer = str(record["answer_gene_id"]).strip()
    answer_ids = [token.strip() for token in answer.split("|") if token.strip()]
    answer_terms = _answer_terms(record)

    # --- eligibility gate (runs before ranking) ---
    leaked_tokens = sorted(token for token in tokens if token in answer_terms)

    # Hypotheses from the shared route: a phenotype token selects a pathway/family
    # and its generic annotation keywords. The route is a single shared file applied
    # identically to every case (the manifest records its fingerprint) -- the only
    # per-case leakage vector is the phenotype text naming the answer, checked above.
    derived_pathways, keywords = derive_route(tokens, route)

    if leaked_tokens:
        return {
            "case_id": case_id, "eligibility": "excluded_answer_derived", "exclusion_reason": "answer_derived_phenotype_term",
            "leaked_terms": ";".join(leaked_tokens),
            "answer_gene_id": answer, "derived_pathways": ";".join(sorted(derived_pathways)),
        }

    # --- deterministic no-LLM route ---
    interval_path = runtime_path(str(record["interval_genes"]).strip())
    if not interval_path.is_absolute():
        interval_path = config_dir / interval_path
    interval_genes = pd.read_csv(interval_path, sep="\t", dtype=str).fillna("")
    ranked = rank_interval(interval_genes, keywords)

    metrics = resolution_metrics(ranked, answer_ids)

    # layered attribution
    answer_pathway = str(record.get("answer_pathway", "")).strip()
    pathway_recovered = bool(answer_pathway) and answer_pathway in derived_pathways
    # layer 2: is the answer among genes the derived keywords actually matched?
    matched_genes = set(ranked.loc[ranked["dry_score"].astype(float) > 0, "gene_id"].astype(str))
    gene_in_pathway_matches = bool(set(answer_ids) & matched_genes)
    if metrics["top1_conservative"]:
        layer = "recovered_sole_top"
    elif not answer_pathway:
        layer = "no_pathway_label"
    elif not pathway_recovered:
        layer = "missed_pathway"
    elif not gene_in_pathway_matches:
        layer = "pathway_hit_gene_missed"
    else:
        layer = "pathway_and_gene_hit_not_top"

    return {
        "case_id": case_id, "eligibility": "eligible", "exclusion_reason": "",
        "answer_gene_id": answer, "answer_pathway": answer_pathway,
        "derived_pathways": ";".join(sorted(derived_pathways)),
        "pathway_recovered": pathway_recovered, "gene_in_pathway_matches": gene_in_pathway_matches,
        "failure_layer": layer,
        **metrics,
    }


def _route_lane(records: list[dict[str, object]], lane: str) -> dict[str, object]:
    eligible = [row for row in records if row.get("eligibility") == "eligible"]
    resolution = aggregate_resolution_metrics(eligible) if eligible else {}
    pathway_present = [row for row in eligible if str(row.get("answer_pathway", ""))]
    return {
        "lane": lane,
        "eligible_cases": len(eligible),
        "pathway_recall": round(sum(1 for row in pathway_present if row.get("pathway_recovered")) / len(pathway_present), 6) if pathway_present else 0.0,
        "gene_in_pathway_recall": round(sum(1 for row in eligible if row.get("gene_in_pathway_matches")) / len(eligible), 6) if eligible else 0.0,
        "resolution_metrics": resolution,
    }


def blind_benchmark(
    manifest_arg: str,
    route_arg: str,
    out_arg: str,
    llm_rankings_arg: str | None = None,
    full_hash: bool = False,
) -> dict[str, object]:
    """Run the phenotype-first blind (autonomous-floor) benchmark."""
    manifest_path = runtime_path(manifest_arg)
    route_path = runtime_path(route_arg)
    out = runtime_path(out_arg)
    out.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str).fillna("")
    missing = [column for column in MANIFEST_COLUMNS if column not in manifest.columns]
    if missing:
        raise ValueError(f"blind benchmark manifest missing columns: {', '.join(missing)}")
    if manifest.empty:
        raise ValueError("blind benchmark manifest is empty")
    if manifest["case_id"].duplicated().any():
        raise ValueError("blind benchmark case_id must be unique")
    route = pd.read_csv(route_path, sep="\t", dtype=str).fillna("")
    route_missing = [column for column in ROUTE_COLUMNS if column not in route.columns]
    if route_missing:
        raise ValueError(f"phenotype route missing columns: {', '.join(route_missing)}")

    config_dir = manifest_path.parent
    records = [_evaluate_case(record, route, config_dir) for record in manifest.to_dict("records")]
    cases = pd.DataFrame(records)
    cases_path = out / "blind_cases.tsv"
    cases.to_csv(cases_path, sep="\t", index=False, lineterminator="\n")

    excluded = cases[cases["eligibility"] != "eligible"]
    deterministic = _route_lane(records, "deterministic_no_llm")

    # optional LLM lane over the SAME frozen interval-gene set, scored separately.
    llm_lane = None
    if llm_rankings_arg:
        # The LLM-rankings directory is a top-level argument (resolved like the
        # manifest itself), not a manifest-relative path like interval_genes.
        llm_dir = runtime_path(llm_rankings_arg)
        llm_records: list[dict[str, object]] = []
        for record in records:
            if record.get("eligibility") != "eligible":
                llm_records.append(record)
                continue
            ranking_path = llm_dir / f"{record['case_id']}.tsv"
            if not ranking_path.exists():
                llm_records.append({**record, "eligibility": "eligible", "answer_pathway": record.get("answer_pathway", "")})
                continue
            ranked = pd.read_csv(ranking_path, sep="\t", dtype=str).fillna("")
            answer_ids = [token.strip() for token in str(record["answer_gene_id"]).split("|") if token.strip()]
            llm_records.append({**record, **resolution_metrics(ranked, answer_ids)})
        llm_lane = _route_lane(llm_records, "optional_llm")

    manifest_out = {
        "schema_version": "1.0",
        "benchmark_version": BENCHMARK_VERSION,
        "benchmark_type": "phenotype-first blind autonomous-floor discovery",
        "distinct_from": "map-based homolog-informed benchmark (benchmark-reference / discovery-benchmark): a human-chosen pathway homolog is NOT provided here",
        "case_count": len(cases),
        "eligible_count": int((cases["eligibility"] == "eligible").sum()),
        "excluded_answer_derived": int((cases["eligibility"] == "excluded_answer_derived").sum()),
        "excluded_case_ids": sorted(excluded["case_id"]),
        "deterministic_no_llm": deterministic,
        "optional_llm": llm_lane,
        "manifest_fingerprint": file_fingerprint(manifest_path, full_hash),
        "route_fingerprint": file_fingerprint(route_path, full_hash),
        "policy": [
            "the phenotype->pathway route is a single shared file applied identically to every case (no per-case tuning)",
            "cases whose phenotype terms leak the answer, or whose route keyword coincides with the answer, are excluded before ranking",
            "the no-LLM deterministic lane is the mandatory baseline; the optional LLM lane is scored separately on the same frozen inputs",
            "initial results are recorded as-is, including 0%; improvement must come from new provenance-recorded signal, not weight tuning",
        ],
    }
    manifest_json = out / "blind_benchmark_manifest.json"
    manifest_json.write_text(json.dumps(manifest_out, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write_report(cases, manifest_out, out)
    return {
        "cases": str(cases_path),
        "manifest": str(manifest_json),
        "report": str(out / "blind_report.md"),
        "eligible_count": manifest_out["eligible_count"],
        "excluded_answer_derived": manifest_out["excluded_answer_derived"],
        "deterministic_no_llm": deterministic,
        "optional_llm": llm_lane,
    }


def _mrr_ci_suffix(resolution: dict) -> str:
    """Conservative MRR bootstrap CI suffix; empty when undefined."""
    low = resolution.get("mrr_conservative_ci_95_low")
    high = resolution.get("mrr_conservative_ci_95_high")
    if low is None or high is None:
        return ""
    return f" (95% CI {float(low):.3f}-{float(high):.3f})"


def _write_report(cases: pd.DataFrame, manifest: dict[str, object], out) -> None:
    det = manifest["deterministic_no_llm"]
    resolution = det.get("resolution_metrics", {}) or {}
    lines = [
        "# Phenotype-first blind benchmark (autonomous floor)",
        "",
        "Distinct from the map-based homolog-informed benchmark: no human-chosen",
        "pathway homolog is provided. This measures deriving the pathway/family",
        "from the phenotype alone. Numbers are recorded as-is, including 0%.",
        "",
        f"- cases: {manifest['case_count']}; eligible: {manifest['eligible_count']}; "
        f"excluded (answer-derived): {manifest['excluded_answer_derived']}",
        f"- deterministic no-LLM lane -- pathway recall: {det.get('pathway_recall', 0):.1%}; "
        f"gene-in-pathway recall: {det.get('gene_in_pathway_recall', 0):.1%}",
        f"- deterministic Top-1 / Top-20 recall (conservative): "
        f"{resolution.get('top1_recall_conservative', 0):.1%} / {resolution.get('top20_recall_conservative', 0):.1%}; "
        f"MRR: {resolution.get('mrr_conservative', 0):.3f}{_mrr_ci_suffix(resolution)}",
        f"- wrong sole top: {resolution.get('wrong_sole_top_count', 0)} case(s); "
        f"abstention (truth-missing): {resolution.get('abstention_rate', 0):.1%}",
    ]
    if manifest.get("optional_llm"):
        llm = manifest["optional_llm"]
        llm_res = llm.get("resolution_metrics", {}) or {}
        lines.append(
            f"- optional LLM lane (separate) -- pathway recall: {llm.get('pathway_recall', 0):.1%}; "
            f"Top-1 recall: {llm_res.get('top1_recall_conservative', 0):.1%}"
        )
    lines.extend(["", "## Failure attribution (eligible cases)", "", "| case | layer | pathway recovered | gene in pathway | rank_max |", "|---|---|---|---|---:|"])
    for row in cases[cases["eligibility"] == "eligible"].to_dict("records"):
        lines.append(
            f"| {row['case_id']} | {row.get('failure_layer', '')} | {row.get('pathway_recovered', '')} | "
            f"{row.get('gene_in_pathway_matches', '')} | {row.get('rank_max', '')} |"
        )
    (out / "blind_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
