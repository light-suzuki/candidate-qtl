"""Answerless prediction / freeze, then separate post-freeze scoring (issue #213).

``blind-benchmark`` measures the autonomous floor, but it takes the answer and
the prediction in one call: the manifest requires ``answer_gene_id`` *before*
anything is ranked, and the per-case ranking it computes is never written out.
A monitored audit therefore had to invent a sentinel answer column just to run
the prediction phase, and could not prove that the prediction had been fixed
before the truth was opened.

This module splits that into the two phases an audit actually needs:

``audit-predict`` (prediction only)
    Takes a manifest with **no answer-shaped column at all** -- one is actively
    rejected, so truth cannot ride in on a spare column -- derives the route from
    the phenotype text through the same shared route table, and writes one
    deterministic ranking per case plus a route-status table. It then freezes
    everything: SHA-256 per artifact, UTC timestamp, input fingerprints, command
    and version metadata, and per-lane evaluability. Truth-dependent lanes are
    recorded as ``not_evaluated`` rather than failing the run: an unavailable
    lane and a broken lane are different states.

``audit-score`` (post-freeze)
    Opens the truth for the first time, re-hashes every frozen artifact and
    refuses to score if one changed, and only then computes recall / MRR /
    failure attribution. It writes exclusively into its own output directory, so
    scoring cannot mutate a frozen ranking even by accident.

The answer-derived eligibility gate still runs, but it moves to the scoring
phase: at prediction time the answer is unknown by construction, so a case whose
phenotype text leaks the answer can only be detected once the truth is opened.
That ordering is stricter than doing it up front, not weaker -- the prediction
never saw the answer either way.

Nothing here ranks differently from ``blind_benchmark``: the token, route and
scoring functions are imported from it rather than reimplemented, so the frozen
prediction is the same object the blind benchmark would have produced.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from . import __version__
from .benchmark_metrics import aggregate_resolution_metrics, resolution_metrics
from .blind_benchmark import ROUTE_COLUMNS, derive_route, phenotype_tokens, rank_interval
from .core import file_fingerprint, runtime_path, safe_id, sha256

FREEZE_VERSION = "1.0.0"
SCHEMA_VERSION = "1.0"

PREDICT_MANIFEST_COLUMNS = ["case_id", "phenotype_text", "interval", "interval_genes"]

# Any column that could carry the answer into the prediction phase. The check is
# on the column *name*: a prediction manifest that carries one is rejected, so a
# sentinel/placeholder answer column is not needed and not accepted either.
from .prediction_contract import ANSWER_SHAPED_COLUMNS

TRUTH_COLUMNS = ["case_id", "answer_gene_id"]
TRUTH_OPTIONAL_COLUMNS = ("answer_pathway", "answer_gene_name", "answer_locus")

RANKING_COLUMNS = ["gene_id", "dry_score", "dry_group", "rank_min", "rank_max"]
ROUTE_STATUS_COLUMNS = [
    "case_id", "interval", "interval_gene_count", "phenotype_tokens", "phenotype_token_count",
    "derived_pathways", "derived_pathway_count", "keyword_count", "keyword_matched_gene_count",
    "route_status", "ranking_status", "ranking_tsv", "ranking_sha256",
]
SCORE_CASE_COLUMNS = [
    "case_id", "eligibility", "exclusion_reason", "leaked_terms", "answer_gene_id", "answer_pathway",
    "derived_pathways", "pathway_recovered", "gene_in_pathway_matches", "failure_layer",
]

# Lanes that need the truth. In prediction-only mode they are unavailable, which
# is a state, not a failure.
TRUTH_DEPENDENT_LANES = ("robustness_benchmark", "answer_recall", "answer_derived_eligibility_gate")

ROUTE_STATUS_NAME = "route_status.tsv"
FREEZE_NAME = "freeze_manifest.json"
PREDICTIONS_DIR = "predictions"


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def _answer_terms(record: dict[str, object]) -> set[str]:
    terms: set[str] = set()
    for column in ("answer_gene_id", "answer_gene_name", "answer_locus"):
        for token in str(record.get(column, "")).casefold().replace("|", " ").replace(";", " ").replace(",", " ").split():
            if len(token) >= 2:
                terms.add(token)
    return terms


def _tie_aware_ranks(scores: list[float]) -> tuple[list[int], list[int]]:
    """rank_min / rank_max for every row, computed only from the score column so
    the frozen ranking carries no order-dependent tie break."""
    rank_min = [int(sum(1 for other in scores if other > score)) + 1 for score in scores]
    rank_max = [int(sum(1 for other in scores if other >= score)) for score in scores]
    return rank_min, rank_max


# --------------------------------------------------------------------------- #
# phase 1: answerless prediction + freeze
# --------------------------------------------------------------------------- #

def audit_predict(
    manifest_arg: str,
    route_arg: str,
    out_arg: str,
    case_id: str = "",
    full_hash: bool = False,
) -> dict[str, object]:
    """Run the deterministic blind route with no truth input, and freeze it."""
    manifest_path = runtime_path(manifest_arg)
    route_path = runtime_path(route_arg)
    out = runtime_path(out_arg)
    predictions_dir = out / PREDICTIONS_DIR
    predictions_dir.mkdir(parents=True, exist_ok=True)

    manifest = _read_tsv(manifest_path)
    missing = [column for column in PREDICT_MANIFEST_COLUMNS if column not in manifest.columns]
    if missing:
        raise ValueError(f"prediction manifest missing columns: {', '.join(missing)}")
    if manifest.empty:
        raise ValueError("prediction manifest is empty")
    if manifest["case_id"].duplicated().any():
        raise ValueError("prediction manifest case_id must be unique")
    answer_shaped = sorted(set(manifest.columns) & ANSWER_SHAPED_COLUMNS)
    if answer_shaped:
        raise ValueError(
            "prediction manifest carries answer/truth column(s): " + ", ".join(answer_shaped)
            + ". The prediction phase must not receive truth in any form; supply them to `audit-score` instead"
        )

    route = _read_tsv(route_path)
    route_missing = [column for column in ROUTE_COLUMNS if column not in route.columns]
    if route_missing:
        raise ValueError(f"phenotype route missing columns: {', '.join(route_missing)}")

    config_dir = manifest_path.parent
    rows: list[dict[str, object]] = []
    interval_inputs: dict[str, object] = {}
    for record in manifest.to_dict("records"):
        case = str(record["case_id"]).strip()
        if not case:
            raise ValueError("prediction manifest contains an empty case_id")
        tokens = phenotype_tokens(record["phenotype_text"])
        pathways, keywords = derive_route(tokens, route)

        interval_path = runtime_path(str(record["interval_genes"]).strip())
        if not interval_path.is_absolute():
            interval_path = config_dir / interval_path
        if not interval_path.exists():
            raise ValueError(f"interval gene table for case {case} does not exist: {interval_path}")
        interval_genes = _read_tsv(interval_path)
        ranked = rank_interval(interval_genes, keywords).copy()
        scores = [float(value) for value in ranked["dry_score"]]
        rank_min, rank_max = _tie_aware_ranks(scores)
        ranked["rank_min"], ranked["rank_max"] = rank_min, rank_max
        ranked = ranked[RANKING_COLUMNS].sort_values(["rank_max", "gene_id"], kind="stable")
        ranking_path = predictions_dir / f"{safe_id(case)}.tsv"
        ranked.to_csv(ranking_path, sep="\t", index=False, lineterminator="\n")

        matched = int(sum(1 for score in scores if score > 0))
        groups = set(ranked["dry_group"].astype(str))
        ranking_status = ("sole_top" if "top" in groups
                          else ("tied_or_reserve_only" if "reserve" in groups else "abstained_no_keyword_match"))
        interval_inputs[case] = file_fingerprint(interval_path, full_hash)
        rows.append({
            "case_id": case, "interval": str(record["interval"]).strip(),
            "interval_gene_count": int(len(ranked)),
            "phenotype_tokens": ";".join(tokens), "phenotype_token_count": len(tokens),
            "derived_pathways": ";".join(sorted(pathways)), "derived_pathway_count": len(pathways),
            "keyword_count": len(keywords), "keyword_matched_gene_count": matched,
            "route_status": "route_hypothesis_derived" if pathways else "no_route_hypothesis",
            "ranking_status": ranking_status,
            "ranking_tsv": str(ranking_path), "ranking_sha256": sha256(ranking_path),
        })

    status = pd.DataFrame(rows, columns=ROUTE_STATUS_COLUMNS).sort_values("case_id", kind="stable")
    status_path = out / ROUTE_STATUS_NAME
    status.to_csv(status_path, sep="\t", index=False, lineterminator="\n")

    freeze = {
        "schema_version": SCHEMA_VERSION, "freeze_version": FREEZE_VERSION,
        "phase": "prediction_only",
        "command": "audit-predict",
        "candidate_qtl_version": __version__,
        "case_id": case_id,
        "frozen_at_utc": _now(),
        "case_count": int(len(status)),
        "truth_inputs_present": False,
        "answer_columns_rejected": sorted(ANSWER_SHAPED_COLUMNS),
        "prediction_contract": (
            "no answer/causal/truth column is accepted by this phase; the ranking is derived from the "
            "phenotype text and the shared route only, and is frozen before any truth file is opened"
        ),
        "ranking_provenance": (
            "phenotype tokens and the interval ranking come from candidate_qtl.blind_benchmark "
            "(phenotype_tokens / derive_route / rank_interval); this phase adds no ranking rule of its own"
        ),
        "route_status_sha256": sha256(status_path),
        "cases": {
            str(row["case_id"]): {
                "ranking_tsv": str(row["ranking_tsv"]), "sha256": str(row["ranking_sha256"]),
                "route_status": str(row["route_status"]), "ranking_status": str(row["ranking_status"]),
                "interval_gene_count": int(row["interval_gene_count"]),
            } for row in status.to_dict("records")
        },
        "route_evaluability": {
            "cases_with_route_hypothesis": int((status["route_status"] == "route_hypothesis_derived").sum()),
            "cases_without_route_hypothesis": int((status["route_status"] == "no_route_hypothesis").sum()),
            "cases_with_sole_top": int((status["ranking_status"] == "sole_top").sum()),
            "cases_abstained": int((status["ranking_status"] == "abstained_no_keyword_match").sum()),
            "rule": "a case whose phenotype matched no route token produced no hypothesis; that is a stage-1 miss, recorded here rather than hidden inside a rank",
        },
        "lanes": {
            "deterministic_ranking": {"status": "evaluated", "reason": "the deterministic no-LLM route ran for every case"},
            **{lane: {"status": "not_evaluated",
                      "reason": "truth-dependent lane; unavailable during prediction-only mode and not a stage failure"}
               for lane in TRUTH_DEPENDENT_LANES},
        },
        "inputs": {
            "manifest": file_fingerprint(manifest_path, full_hash),
            "route": file_fingerprint(route_path, full_hash),
            "interval_genes": interval_inputs,
        },
        "outputs": {"route_status": str(status_path), "predictions_dir": str(predictions_dir),
                    "freeze_manifest": str(out / FREEZE_NAME)},
        "policy": [
            "prediction is frozen with SHA-256 before any truth is read; scoring is a separate command",
            "an answer-shaped manifest column is rejected, so a sentinel answer is neither required nor accepted",
            "a truth-dependent lane reports not_evaluated instead of failing the profile",
            "no weight, rank rule, or keyword is changed here; improvement must come from new signal",
        ],
    }
    freeze_path = out / FREEZE_NAME
    freeze_path.write_text(json.dumps(freeze, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "freeze_manifest": str(freeze_path), "route_status": str(status_path),
        "predictions_dir": str(predictions_dir), "case_count": int(len(status)),
        "cases_with_route_hypothesis": freeze["route_evaluability"]["cases_with_route_hypothesis"],
        "cases_abstained": freeze["route_evaluability"]["cases_abstained"],
        "truth_inputs_present": False,
    }


# --------------------------------------------------------------------------- #
# phase 2: post-freeze scoring
# --------------------------------------------------------------------------- #

def _verify_freeze(freeze: dict[str, object], freeze_path: Path) -> tuple[pd.DataFrame, list[str]]:
    """Re-hash every frozen artifact; scoring refuses to run on a mutated freeze."""
    status_path = runtime_path(str(freeze["outputs"]["route_status"]))  # type: ignore[index]
    if not status_path.exists():
        status_path = freeze_path.parent / ROUTE_STATUS_NAME
    mutated: list[str] = []
    if not status_path.exists():
        raise ValueError(f"frozen route status table is missing: {status_path}")
    if sha256(status_path) != str(freeze.get("route_status_sha256", "")):
        mutated.append(str(status_path))
    status = _read_tsv(status_path)
    cases = freeze.get("cases", {})
    assert isinstance(cases, dict)
    for case, record in sorted(cases.items()):
        ranking_path = runtime_path(str(record["ranking_tsv"]))
        if not ranking_path.exists():
            ranking_path = freeze_path.parent / PREDICTIONS_DIR / f"{safe_id(case)}.tsv"
        if not ranking_path.exists():
            mutated.append(f"{case}:missing")
            continue
        if sha256(ranking_path) != str(record["sha256"]):
            mutated.append(f"{case}:{ranking_path}")
    return status, mutated


def verify_freeze_inventory(
    freeze_manifest: dict[str, object],
    roots: list[Path],
    max_workers: int = 8,
) -> dict[str, object]:
    """Re-verify a T8-style ``files[]`` freeze inventory against one or more roots (issue #348).

    The inventory-shaped manifest records one entry per artifact file as
    ``{"path", "size", "sha256"}`` (the audit-predict manifest instead records
    per-case ranking records; see ``_verify_freeze`` for that shape). Every
    entry is re-checked against the first root that contains its relative path,
    roots tried in order, and classified:

    - ``matched``: size and SHA-256 both equal the frozen values
    - ``mismatched``: the file exists but size or SHA-256 differ (expected and
      actual values are recorded so the caller can report the delta)
    - ``missing``: no root contains the path
    - ``invalid``: the entry is unusable (not a dict, or path is absolute or
      escapes a root with ``..``); it is reported, never silently skipped

    This never raises on content drift and never claims verification for a root
    that was not checked: the caller sees exactly which root_index each entry
    was resolved from (``None`` for missing/invalid) and decides what the
    availability means for its fixture layout.
    """
    import concurrent.futures

    files = freeze_manifest.get("files", [])
    if not isinstance(files, list):
        raise ValueError("freeze manifest has no files[] inventory")
    if not roots:
        raise ValueError("verify_freeze_inventory needs at least one root")

    def check(entry: object) -> dict[str, object]:
        base: dict[str, object] = {"status": "invalid", "root_index": None}
        if not isinstance(entry, dict):
            return base
        rel_text = str(entry.get("path", "")).strip()
        if not rel_text or rel_text.startswith(("/", "\\")) or ".." in rel_text.split("/"):
            base["path"] = rel_text or "<missing path>"
            return base
        rel = Path(rel_text)
        base.update({
            "path": rel_text,
            "expected_size": entry.get("size"),
            "expected_sha256": entry.get("sha256"),
        })
        for index, root in enumerate(roots):
            candidate = root / rel
            try:
                if not candidate.is_file():
                    continue
            except OSError:
                continue
            actual_size = candidate.stat().st_size
            actual_sha = sha256(candidate)
            base["root_index"] = index
            if actual_size == base["expected_size"] and actual_sha == base["expected_sha256"]:
                base["status"] = "matched"
            else:
                base.update({"status": "mismatched",
                             "actual_size": actual_size, "actual_sha256": actual_sha})
            return base
        base["status"] = "missing"
        return base

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as pool:
        entries = list(pool.map(check, files))
    counts = {status: sum(1 for entry in entries if entry["status"] == status)
              for status in ("matched", "mismatched", "missing", "invalid")}
    counts["total"] = len(entries)
    return {"counts": counts, "entries": entries}


def audit_score(
    freeze_arg: str,
    truth_arg: str,
    out_arg: str,
    full_hash: bool = False,
    non_gene_truth_case_reasons: dict[str, str] | None = None,
) -> dict[str, object]:
    """Score a frozen prediction against the truth, without touching the freeze."""
    freeze_path = runtime_path(freeze_arg)
    truth_path = runtime_path(truth_arg)
    out = runtime_path(out_arg)
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    if str(freeze.get("phase")) != "prediction_only":
        raise ValueError("audit-score expects a prediction-only freeze manifest from `audit-predict`")
    if out.resolve() == freeze_path.parent.resolve():
        raise ValueError("audit-score must write to a directory outside the freeze; scoring never writes into a frozen prediction")
    out.mkdir(parents=True, exist_ok=True)
    non_gene_truth_case_reasons = non_gene_truth_case_reasons or {}

    status, mutated = _verify_freeze(freeze, freeze_path)
    if mutated:
        raise ValueError("frozen prediction artifact(s) changed after the freeze; refusing to score: " + ", ".join(mutated))

    truth = _read_tsv(truth_path)
    truth_missing = [column for column in TRUTH_COLUMNS if column not in truth.columns]
    if truth_missing:
        raise ValueError(f"truth table missing columns: {', '.join(truth_missing)}")
    if truth["case_id"].duplicated().any():
        raise ValueError("truth table case_id must be unique")
    truth_records = {str(row["case_id"]).strip(): row for row in truth.to_dict("records")}

    rows: list[dict[str, object]] = []
    for record in status.to_dict("records"):
        case = str(record["case_id"]).strip()
        derived = str(record["derived_pathways"])
        base: dict[str, object] = {
            "case_id": case, "derived_pathways": derived,
            "route_status": str(record["route_status"]), "ranking_status": str(record["ranking_status"]),
        }
        answer_row = truth_records.get(case)
        if answer_row is None:
            non_gene_reason = non_gene_truth_case_reasons.get(case)
            if non_gene_reason:
                rows.append({**base, "eligibility": "not_scored_non_gene_truth_unit", "exclusion_reason": non_gene_reason,
                             "leaked_terms": "", "answer_gene_id": "", "answer_pathway": "",
                             "pathway_recovered": "", "gene_in_pathway_matches": "", "failure_layer": "not_scored"})
                continue
            rows.append({**base, "eligibility": "not_scored_no_truth_row", "exclusion_reason": "no truth row for this frozen case",
                         "leaked_terms": "", "answer_gene_id": "", "answer_pathway": "",
                         "pathway_recovered": "", "gene_in_pathway_matches": "", "failure_layer": "not_scored"})
            continue
        answer = str(answer_row.get("answer_gene_id", "")).strip()
        answer_ids = [token.strip() for token in answer.split("|") if token.strip()]
        answer_pathway = str(answer_row.get("answer_pathway", "")).strip()
        # The answer-derived eligibility gate, run now that the truth is open. The
        # prediction was blind either way; this only decides whether the case may
        # count towards the reported recall.
        tokens = [token for token in str(record["phenotype_tokens"]).split(";") if token]
        leaked = sorted(token for token in tokens if token in _answer_terms(answer_row))
        if leaked:
            rows.append({**base, "eligibility": "excluded_answer_derived",
                         "exclusion_reason": "answer_derived_phenotype_term", "leaked_terms": ";".join(leaked),
                         "answer_gene_id": answer, "answer_pathway": answer_pathway,
                         "pathway_recovered": "", "gene_in_pathway_matches": "", "failure_layer": "excluded"})
            continue

        ranking_path = runtime_path(str(record["ranking_tsv"]))
        if not ranking_path.exists():
            ranking_path = freeze_path.parent / PREDICTIONS_DIR / f"{safe_id(case)}.tsv"
        ranked = _read_tsv(ranking_path)
        metrics = resolution_metrics(ranked, answer_ids)
        derived_pathways = {value for value in derived.split(";") if value}
        pathway_recovered = bool(answer_pathway) and answer_pathway in derived_pathways
        matched_genes = set(ranked.loc[pd.to_numeric(ranked["dry_score"], errors="coerce").fillna(0.0) > 0, "gene_id"].astype(str))
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
        rows.append({**base, "eligibility": "eligible", "exclusion_reason": "", "leaked_terms": "",
                     "answer_gene_id": answer, "answer_pathway": answer_pathway,
                     "pathway_recovered": pathway_recovered, "gene_in_pathway_matches": gene_in_pathway_matches,
                     "failure_layer": layer, **metrics})

    cases = pd.DataFrame(rows).sort_values("case_id", kind="stable")
    cases_path = out / "audit_score_cases.tsv"
    cases.to_csv(cases_path, sep="\t", index=False, lineterminator="\n")

    eligible = [row for row in rows if row.get("eligibility") == "eligible"]
    pathway_present = [row for row in eligible if str(row.get("answer_pathway", ""))]
    summary = {
        "scored_cases": len(eligible),
        "excluded_answer_derived": sum(1 for row in rows if row.get("eligibility") == "excluded_answer_derived"),
        "not_scored_no_truth_row": sum(1 for row in rows if row.get("eligibility") == "not_scored_no_truth_row"),
        "not_scored_non_gene_truth_unit": sum(1 for row in rows if row.get("eligibility") == "not_scored_non_gene_truth_unit"),
        "pathway_recall": round(sum(1 for row in pathway_present if row.get("pathway_recovered")) / len(pathway_present), 6) if pathway_present else None,
        "pathway_recall_evaluated": len(pathway_present),
        "gene_in_pathway_recall": round(sum(1 for row in eligible if row.get("gene_in_pathway_matches")) / len(eligible), 6) if eligible else None,
        "resolution_metrics": aggregate_resolution_metrics(eligible) if eligible else {},
        "failure_layers": {layer: sum(1 for row in eligible if row.get("failure_layer") == layer)
                           for layer in sorted({str(row.get("failure_layer", "")) for row in eligible})},
    }

    manifest = {
        "schema_version": SCHEMA_VERSION, "phase": "post_freeze_scoring", "command": "audit-score",
        "candidate_qtl_version": __version__, "scored_at_utc": _now(),
        "frozen_prediction": {
            "freeze_manifest": str(freeze_path), "frozen_at_utc": freeze.get("frozen_at_utc"),
            "candidate_qtl_version_at_freeze": freeze.get("candidate_qtl_version"),
            "case_count": freeze.get("case_count"),
        },
        "freeze_verification": {
            "verified": True, "mutated_artifacts": [],
            "rule": "every frozen ranking and the route-status table are re-hashed before scoring; a changed artifact aborts the run",
        },
        "scoring_writes_into": str(out),
        "frozen_artifacts_written_by_scoring": [],
        "truth": file_fingerprint(truth_path, full_hash),
        "summary": summary,
        "eligibility_gate": (
            "the answer-derived phenotype-term gate runs here, not at prediction time: the prediction phase "
            "never saw the answer, so a leak can only be detected once the truth is opened"
        ),
        "policy": [
            "the prediction was frozen before this command ran; scoring cannot change it",
            "a case without a truth row is not_scored, never a miss",
            "a Supervisor-provided non-gene truth-unit reason is distinct from a missing truth row; the gene-level audit scorer does not score that unit",
            "truth-dependent lanes that were not_evaluated at prediction stay not_evaluated here unless their input is supplied",
            "no weight or rank rule is changed by scoring",
        ],
        "outputs": {"cases": str(cases_path), "manifest": str(out / "audit_score_manifest.json"),
                    "report": str(out / "audit_score_report.md")},
    }
    manifest_path = out / "audit_score_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write_report(cases, manifest, out)
    return {"cases": str(cases_path), "manifest": str(manifest_path), "report": str(out / "audit_score_report.md"),
            "scored_cases": summary["scored_cases"], "excluded_answer_derived": summary["excluded_answer_derived"],
            "not_scored_no_truth_row": summary["not_scored_no_truth_row"],
            "not_scored_non_gene_truth_unit": summary["not_scored_non_gene_truth_unit"], "summary": summary}


def _mrr_ci_suffix(resolution: dict) -> str:
    """Conservative MRR bootstrap CI suffix; empty when undefined."""
    low = resolution.get("mrr_conservative_ci_95_low")
    high = resolution.get("mrr_conservative_ci_95_high")
    if low is None or high is None:
        return ""
    return f" (95% CI {float(low):.3f}-{float(high):.3f})"


def _write_report(cases: pd.DataFrame, manifest: dict[str, object], out: Path) -> None:
    summary = manifest["summary"]
    assert isinstance(summary, dict)
    resolution = summary.get("resolution_metrics", {}) or {}
    frozen = manifest["frozen_prediction"]
    assert isinstance(frozen, dict)

    def rate(value: object) -> str:
        return "not_evaluated" if value is None else f"{float(value):.1%}"

    lines = [
        "# Post-freeze audit score",
        "",
        f"Prediction frozen at {frozen.get('frozen_at_utc')} by candidate-qtl "
        f"{frozen.get('candidate_qtl_version_at_freeze')}; scored at {manifest.get('scored_at_utc')}.",
        "Every frozen ranking was re-hashed before scoring, so these numbers are for the",
        "prediction as it stood before the truth was opened.",
        "",
        f"- scored cases: {summary['scored_cases']}; excluded (answer-derived phenotype term): "
        f"{summary['excluded_answer_derived']}; not scored (no truth row): {summary['not_scored_no_truth_row']}; "
        f"not scored (non-gene truth unit): {summary['not_scored_non_gene_truth_unit']}",
        f"- pathway recall: {rate(summary['pathway_recall'])} "
        f"(evaluated on {summary['pathway_recall_evaluated']} case(s) with an answer_pathway label)",
        f"- gene-in-pathway recall: {rate(summary['gene_in_pathway_recall'])}",
        f"- Top-1 / Top-20 recall (conservative): {resolution.get('top1_recall_conservative', 0):.1%} / "
        f"{resolution.get('top20_recall_conservative', 0):.1%}; MRR: {resolution.get('mrr_conservative', 0):.3f}"
        f"{_mrr_ci_suffix(resolution)}",
        f"- wrong sole top: {resolution.get('wrong_sole_top_count', 0)} case(s)",
        "", "## Failure attribution", "",
        "| case | eligibility | route | ranking | layer | pathway recovered | gene in pathway | rank_max |",
        "|---|---|---|---|---|---|---|---:|",
    ]
    for row in cases.to_dict("records"):
        lines.append(
            f"| {row.get('case_id', '')} | {row.get('eligibility', '')} | {row.get('route_status', '')} | "
            f"{row.get('ranking_status', '')} | {row.get('failure_layer', '')} | {row.get('pathway_recovered', '')} | "
            f"{row.get('gene_in_pathway_matches', '')} | {row.get('rank_max', '')} |"
        )
    lines += ["", "A `no_route_hypothesis` route status means the phenotype matched no shared-route token: "
              "the miss happened before any gene was ranked, and improving the rank cannot fix it."]
    (out / "audit_score_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
