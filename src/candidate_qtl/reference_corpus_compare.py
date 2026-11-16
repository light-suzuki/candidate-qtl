"""Compare reference corpora against each other (issue #232).

Enriching the corpus is a claim, and a claim needs a measurement. This command
answers two questions about two or more frozen corpora, without touching a single
ranking weight:

1. **What biology does each corpus actually contain?** Row count is the wrong
   number -- the description-heavy corpus that failed the fasciation query had
   39,254 rows. Coverage per field (mutant phenotype, pathway, GO, tissue,
   evidence) and per species is the number that moves when a corpus gets better.
2. **Does the extra biology change what the phenotype retrieves?** Given frozen
   Route C retrieval artifacts produced from each corpus for the same cases, the
   same Family Recall@5/10/20 and Gene Recall@20 are recomputed per corpus and
   reported as a delta against a baseline label.

Rules kept from the benchmarks this borrows from:

* A metric divides by the cases that could be evaluated. A case with no
  ``truth_gene_id`` is ``not_evaluated`` for gene recall, never a miss.
* Retrieval scores are never compared across corpora -- only ranks and hits are,
  because each corpus is its own database with its own scale.
* Nothing here writes a score onto a candidate. The comparison decides which
  reference resource to build on, not which gene is causal.
"""

from __future__ import annotations

import json

import pandas as pd

from . import __version__
from .core import file_fingerprint, runtime_path
from .family_labels import family_set_contains
from .reference_sources import ALL_SPECIES, COVERAGE_FIELDS, corpus_coverage
from .retrieval_route import validate_retrieval_frame
from .route_c_benchmark import _c_family_rank, _top_families, _top_genes

SCHEMA_VERSION = "1.0"
COMPARE_VERSION = "1.0.0"

COVERAGE_COLUMNS = ["field", "scope", "reference_species", "corpus_label", "rows", "non_empty",
                    "coverage", "baseline_coverage", "coverage_delta"]
CORPUS_SUMMARY_COLUMNS = ["corpus_label", "corpus_path", "rows", "species_count", "phenotype_rows",
                          "pathway_rows", "experimental_evidence_rows", "shared_accessions_with_baseline",
                          "accessions_only_here", "accessions_only_in_baseline"]
RETRIEVAL_CASE_COLUMNS = ["case_id", "corpus_label", "truth_family_id", "truth_gene_id", "retrieval_rows",
                          "family_rank", "family_recall_5", "family_recall_10", "family_recall_20", "gene_recall_20"]
RETRIEVAL_SUMMARY_COLUMNS = ["metric", "corpus_label", "cases", "evaluated", "hits", "recall",
                             "not_evaluated", "baseline_recall", "recall_delta"]
RETRIEVAL_MANIFEST_COLUMNS = ["corpus_label", "case_id", "retrieval_tsv", "truth_family_id"]
METRICS = ("family_recall_5", "family_recall_10", "family_recall_20", "gene_recall_20")


def _labelled(values: list[str] | None, flag: str) -> dict[str, str]:
    """``LABEL=PATH`` pairs, order preserved."""
    parsed: dict[str, str] = {}
    for value in values or []:
        label, separator, path = str(value).partition("=")
        if not separator or not label.strip() or not path.strip():
            raise ValueError(f"{flag} expects LABEL=PATH, got {value!r}")
        if label.strip() in parsed:
            raise ValueError(f"{flag} label {label.strip()!r} was given twice")
        parsed[label.strip()] = path.strip()
    return parsed


def reference_corpus_compare(
    corpus_args: list[str],
    out_dir_arg: str,
    baseline_label: str = "",
    retrieval_manifest_arg: str | None = None,
    full_hash: bool = False,
) -> dict[str, object]:
    """Coverage and (optionally) retrieval comparison across frozen corpora."""
    corpora = _labelled(corpus_args, "--corpus")
    if len(corpora) < 2:
        raise ValueError("reference-corpus-compare requires at least two --corpus LABEL=PATH inputs")
    baseline = str(baseline_label).strip() or next(iter(corpora))
    if baseline not in corpora:
        raise ValueError(f"--baseline-label {baseline!r} is not one of the supplied corpus labels: {', '.join(corpora)}")
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)

    frames: dict[str, pd.DataFrame] = {}
    coverage: dict[str, pd.DataFrame] = {}
    inputs: dict[str, object] = {}
    for label, value in corpora.items():
        path = runtime_path(value)
        frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        if "reference_gene_id" not in frame.columns:
            raise ValueError(f"corpus {path} has no reference_gene_id column")
        frames[label] = frame
        coverage[label] = corpus_coverage(frame)
        inputs[label] = file_fingerprint(path, full_hash)

    baseline_accessions = {str(value).strip().casefold() for value in frames[baseline]["reference_gene_id"]}
    coverage_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for label, frame in frames.items():
        baseline_map = {(str(row["scope"]), str(row["reference_species"]), str(row["field"])): float(row["coverage"])
                        for row in coverage[baseline].to_dict("records")}
        for row in coverage[label].to_dict("records"):
            key = (str(row["scope"]), str(row["reference_species"]), str(row["field"]))
            base = baseline_map.get(key)
            coverage_rows.append({
                "field": row["field"], "scope": row["scope"], "reference_species": row["reference_species"],
                "corpus_label": label, "rows": row["rows"], "non_empty": row["non_empty"],
                "coverage": row["coverage"],
                "baseline_coverage": "" if base is None else base,
                # Empty, not 0.0: a species absent from the baseline has no
                # comparison, which is a different statement from "no coverage".
                "coverage_delta": "" if base is None else round(float(row["coverage"]) - base, 6),
            })
        accessions = {str(value).strip().casefold() for value in frame["reference_gene_id"]}

        def filled(column: str) -> int:
            return int((frame[column].astype(str).str.strip() != "").sum()) if column in frame.columns else 0

        summary_rows.append({
            "corpus_label": label, "corpus_path": corpora[label], "rows": int(len(frame)),
            "species_count": len({str(value).strip() for value in frame.get("reference_species", pd.Series(dtype=str))
                                  if str(value).strip()}),
            "phenotype_rows": filled("mutant_phenotype"), "pathway_rows": filled("pathway"),
            "experimental_evidence_rows": int((frame["evidence_type"].astype(str).str.strip() == "experimental").sum())
            if "evidence_type" in frame.columns else 0,
            "shared_accessions_with_baseline": len(accessions & baseline_accessions),
            "accessions_only_here": len(accessions - baseline_accessions),
            "accessions_only_in_baseline": len(baseline_accessions - accessions),
        })

    coverage_frame = pd.DataFrame(coverage_rows, columns=COVERAGE_COLUMNS).sort_values(
        ["field", "scope", "reference_species", "corpus_label"], kind="stable")
    coverage_path = out / "corpus_coverage_comparison.tsv"
    coverage_frame.to_csv(coverage_path, sep="\t", index=False, lineterminator="\n")
    summary_frame = pd.DataFrame(summary_rows, columns=CORPUS_SUMMARY_COLUMNS)
    summary_path = out / "corpus_summary.tsv"
    summary_frame.to_csv(summary_path, sep="\t", index=False, lineterminator="\n")

    retrieval_cases = pd.DataFrame(columns=RETRIEVAL_CASE_COLUMNS)
    retrieval_summary = pd.DataFrame(columns=RETRIEVAL_SUMMARY_COLUMNS)
    retrieval_status = "not_evaluated_no_retrieval_manifest"
    if retrieval_manifest_arg:
        retrieval_cases, retrieval_summary = _retrieval_comparison(
            retrieval_manifest_arg, set(frames), baseline)
        retrieval_status = "evaluated" if len(retrieval_cases) else "not_evaluated_manifest_had_no_usable_case"
    retrieval_cases_path = out / "corpus_retrieval_cases.tsv"
    retrieval_cases.to_csv(retrieval_cases_path, sep="\t", index=False, lineterminator="\n")
    retrieval_summary_path = out / "corpus_retrieval_summary.tsv"
    retrieval_summary.to_csv(retrieval_summary_path, sep="\t", index=False, lineterminator="\n")

    overall = {
        row["corpus_label"]: {field: float(value) for field, value in
                              coverage_frame[(coverage_frame["corpus_label"] == row["corpus_label"])
                                             & (coverage_frame["scope"] == "all")]
                              .set_index("field")["coverage"].items()}
        for row in summary_rows
    }
    manifest = {
        "schema_version": SCHEMA_VERSION, "compare_version": COMPARE_VERSION,
        "command": "reference-corpus-compare", "candidate_qtl_version": __version__,
        "baseline_label": baseline,
        "corpus_labels": list(frames),
        "corpus_inputs": inputs,
        "coverage_fields": list(COVERAGE_FIELDS),
        "overall_coverage": overall,
        "corpus_summary": summary_rows,
        "retrieval_status": retrieval_status,
        "retrieval_metrics": retrieval_summary.to_dict("records") if len(retrieval_summary) else [],
        "rules": [
            "row count is not coverage: a large corpus with no mutant phenotype is the failure this compares against",
            "retrieval scores are never compared across corpora; only within-corpus ranks and hits are",
            "a metric divides by the cases that could be evaluated; a missing truth_gene_id is not_evaluated, not a miss",
            "ranking weights are untouched: this compares reference resources, never candidates",
            "an absent field is knowledge_unavailable for the routes downstream, never negative biological evidence",
        ],
        "outputs": {"coverage_comparison": str(coverage_path), "corpus_summary": str(summary_path),
                    "retrieval_cases": str(retrieval_cases_path), "retrieval_summary": str(retrieval_summary_path)},
    }
    manifest_path = out / "reference_corpus_comparison_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), "coverage_comparison": str(coverage_path),
            "corpus_summary": str(summary_path), "retrieval_cases": str(retrieval_cases_path),
            "retrieval_summary": str(retrieval_summary_path), "baseline_label": baseline,
            "corpus_labels": list(frames), "retrieval_status": retrieval_status,
            "overall_coverage": overall}


def _retrieval_comparison(
    manifest_arg: str,
    known_labels: set[str],
    baseline: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Recompute Route C recall per corpus label from frozen retrieval artifacts."""
    manifest_path = runtime_path(manifest_arg)
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str).fillna("")
    missing = [column for column in RETRIEVAL_MANIFEST_COLUMNS if column not in manifest.columns]
    if missing:
        raise ValueError(f"retrieval manifest is missing columns: {', '.join(missing)}")
    unknown = sorted({str(value).strip() for value in manifest["corpus_label"]} - known_labels)
    if unknown:
        raise ValueError(f"retrieval manifest names corpus labels that were not supplied as --corpus: {', '.join(unknown)}")
    config_dir = manifest_path.parent

    rows: list[dict[str, object]] = []
    for record in manifest.to_dict("records"):
        label = str(record["corpus_label"]).strip()
        case_id = str(record["case_id"]).strip()
        truth_family = str(record["truth_family_id"]).strip()
        truth_gene = str(record.get("truth_gene_id", "")).strip()
        if not (label and case_id and truth_family):
            raise ValueError("retrieval manifest rows need corpus_label, case_id and truth_family_id")
        path = runtime_path(str(record["retrieval_tsv"]).strip())
        if not path.is_absolute():
            path = config_dir / path
        if not path.exists():
            raise ValueError(f"retrieval artifact does not exist: {path}")
        retrieval = validate_retrieval_frame(pd.read_csv(path, sep="\t", dtype=str).fillna(""),
                                             f"retrieval for case {case_id} / corpus {label}")
        retrieval["rank"] = pd.to_numeric(retrieval["rank"], errors="coerce").astype(int)
        hits = {k: family_set_contains(_top_families(retrieval, k), truth_family) for k in (5, 10, 20)}
        family_rank = _c_family_rank(retrieval, truth_family)
        rows.append({
            "case_id": case_id, "corpus_label": label, "truth_family_id": truth_family,
            "truth_gene_id": truth_gene, "retrieval_rows": int(len(retrieval)),
            "family_rank": "" if family_rank is None else family_rank,
            "family_recall_5": bool(hits[5]), "family_recall_10": bool(hits[10]),
            "family_recall_20": bool(hits[20]),
            "gene_recall_20": bool(truth_gene in _top_genes(retrieval, 20)) if truth_gene else "",
        })
    cases = pd.DataFrame(rows, columns=RETRIEVAL_CASE_COLUMNS).sort_values(["case_id", "corpus_label"], kind="stable")

    summary_rows: list[dict[str, object]] = []
    baseline_recall: dict[str, float | None] = {}
    for metric in METRICS:
        for label in sorted({str(value) for value in cases["corpus_label"]}):
            subset = cases[cases["corpus_label"] == label]
            values = subset[metric].astype(str).str.strip()
            evaluated = int(values.ne("").sum())
            hit_count = int(values.eq("True").sum())
            recall = round(hit_count / evaluated, 6) if evaluated else None
            if label == baseline:
                baseline_recall[metric] = recall
            summary_rows.append({"metric": metric, "corpus_label": label, "cases": int(len(subset)),
                                 "evaluated": evaluated, "hits": hit_count,
                                 "recall": "" if recall is None else recall,
                                 "not_evaluated": int(len(subset)) - evaluated,
                                 "baseline_recall": "", "recall_delta": ""})
    for row in summary_rows:
        base = baseline_recall.get(str(row["metric"]))
        if base is None or row["recall"] == "":
            continue
        row["baseline_recall"] = base
        row["recall_delta"] = round(float(row["recall"]) - base, 6)
    summary = pd.DataFrame(summary_rows, columns=RETRIEVAL_SUMMARY_COLUMNS)
    return cases, summary


__all__ = ["reference_corpus_compare", "ALL_SPECIES"]
