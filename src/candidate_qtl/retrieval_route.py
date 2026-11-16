"""Route C: phenotype/protein embedding retrieval as a third discovery route (issue #193).

Two independent ways to reach a pathway hypothesis already exist:

* **Route A - phenotype-first (順引き).** Pin down the gene or family from the
  phenotype and known biology (symbolic / literature-informed).
* **Route B - interval-first reverse lookup (逆引き).** Start from the genes that
  are actually in the interval and ask what they are known to be.

This module adds the third, deliberately different route:

* **Route C - embedding-first / latent retrieval.** Phenotype text is embedded
  and compared against a frozen reference protein / gene / family embedding
  database by similarity. Top-K reference genes / orthogroups / families are
  frozen, and only afterwards - in a separate command - matched against the QTL
  interval.

The separation rule that keeps Route C honest is structural:

* `retrieval-import` (this module's `retrieval_import`) accepts **no interval
  input at all**: no interval, no interval-gene list, no atlas. The external
  retrieval TSV must contain exactly the documented contract columns and nothing
  else, so interval-derived columns cannot ride along.
* Interval matching happens only in `retrieval_interval_match`, against the
  **frozen** retrieval artifact, so the retrieval itself can never be tuned on
  interval contents.

Route C matches are recorded as convergence state in `triple_route_compare`;
they are **never added to a causal score**. Every output carries
`rank_contribution=none`, matching the A/B routes.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from .core import file_fingerprint, runtime_path
from .family_labels import (
    alias_canonical,
    canonical_family_key,
    families_match_expanded,
    split_respecting_brackets,
)
from .phenotype_context import structured_context_block

SCHEMA_VERSION = "1.1"
ROUTE_C = "embedding_retrieval"

# The documented Phase-1 contract for an external retrieval tool. `family_id` /
# `orthogroup_id` are alternatives: at least one must be present. A row-level
# extra column (e.g. an interval gene id) is rejected so interval information
# cannot enter the retrieval input.
RETRIEVAL_COLUMNS = [
    "reference_gene_id", "reference_species", "family_id", "orthogroup_id",
    "retrieval_score", "matched_modality", "matched_description_or_phenotype",
    "model_name", "model_version", "database_version", "rank",
]
FAMILY_COLUMNS = ("family_id", "orthogroup_id")
META_COLUMNS = ("model_name", "model_version", "database_version")

FROZEN_NAME = "route_c_retrieval_frozen.tsv"
MATCH_NAME = "route_c_interval_match.tsv"
STATUS_NAME = "route_c_retrieval_status.json"
MATCH_STATUS_NAME = "route_c_interval_match_status.json"
MERGED_NAME = "route_c_retrieval_merged.tsv"
MERGE_EXCLUDED_NAME = "route_c_retrieval_merge_excluded.tsv"
MERGE_STATUS_NAME = "route_c_retrieval_merge_status.json"

# A/B/C convergence states. The pair-wise A/B vocabulary from
# dual_route is kept so a Route C-absent run reads the same as before.
# Expanded family-label matches added by issue #500 are deliberately separate
# non-strict states: they recover a semantic A/B correspondence without
# pretending that the strict family identity keys were equal.
CONVERGENCE_STATES = (
    "a_b_c_convergent", "a_c_convergent", "b_c_convergent", "a_b_convergent",
    "a_b_expanded_c_supported", "a_b_expanded_match",
    "c_only_supported", "phenotype_only_supported", "interval_only_supported",
    "all_discordant", "all_weak",
)
# Labels that name a real multi-route agreement (used for source-dependence
# accounting; ``c_only_supported`` / single-route states and expanded matches
# are excluded). Expanded matches remain a_and_b_supported downstream but can
# never become strict_convergent merely because their source ids are disjoint.
CONVERGENT_STATES = (
    "a_b_c_convergent", "a_c_convergent", "b_c_convergent", "a_b_convergent",
)
# Source-dependence vocabulary for the parallel ``source_dependence_status``
# field. Routes that converge while citing the same source id
# provide one independent source, not N -- a shared reference is quoted, not
# corroborated. The native score is never touched: the state is reported
# alongside the retained convergence label.
SOURCE_DEPENDENCE_STATUSES = (
    "shared_source_dependent", "source_independent",
    "source_ids_incomplete", "not_applicable",
)
COMPARISON_COLUMNS = [
    "comparison_id", "family_id", "route_a_hypothesis_ids", "route_b_hypothesis_ids",
    "route_c_reference_genes", "route_c_retrieval_ranks", "matched_interval_genes",
    "convergence_status", "family_match_basis",
    # Per-family failure-mode accounting: which route's characteristic error was
    # actually checked off for THIS family, and which is still open.
    "applicable_failure_modes", "neutralised_failure_modes", "residual_failure_modes",
    "failure_mode_coverage", "interval_anchored", "interval_anchor_basis",
    "interval_blind_support", "mechanism_source_backed", "mechanism_source_ids",
    # A and B read the same literature, so agreement is cheap and disagreement is
    # the informative event. Both are reported per family. `four_state` is a
    # diagnostic label on the A/B lanes (which knowledge routes support this
    # family) and is not a score; `literature_link_candidate_count` counts the
    # reverse-literature rows whose interval genes intersect this family's
    # members -- candidates waiting for judgment, never evidence.
    "agreement_state", "agreement_strength", "shared_source_ids",
    # Per-route real source ids. A/B ids are reference-record ids
    # from ``source_backing`` (status strings are dropped); Route C's ids are
    # the ``database_version`` its hit rows were retrieved from.
    # ``source_dependence_status`` flags a convergence that rests on a source id
    # shared between routes -- one reference quoted by N routes is one
    # independent source, not N. ``convergence_status`` is retained as-is.
    "route_a_source_ids", "route_b_source_ids", "route_c_source_ids",
    "source_dependence_status",
    "disagreement_reason",
    "four_state", "literature_link_candidate_count",
    "rank_contribution", "notes",
]
DISAGREEMENT_COLUMNS = ["family_id", "agreement_state", "disagreement_reason", "routes_supporting",
                        "matched_interval_genes", "required_next_action"]
REVIEW_COLUMNS = ["reason", "target", "required_next_action", "allowed_choices"]
# Canonical family-to-copy handoff: pass family_id as join key
# and retain display gene_family for human review.
FAMILY_TO_COPY_INPUT_COLUMNS = ("hypothesis_id", "family_id", "gene_family", "member_interval_genes")

# Route C is not yet load-bearing: a lexical/embedding match may inform review but
# must not move a candidate's position.
ROUTE_C_STATUS = "experimental"

# Why a family that one knowledge route reached was missed by the other. For two
# routes over the same literature, this is where the information is.
DISAGREEMENT_REASONS = {
    "a_only_not_in_interval": "Route A proposed it from the literature but no interval gene carries it",
    "a_only_interval_annotation_thin": "Route A proposed it and the interval has members, but Route B could not classify them",
    "b_only_no_phenotype_link_in_literature": "the interval carries it, but no literature links this family to the phenotype -- a novelty candidate",
    "c_only_similarity_without_knowledge": "only the experimental retrieval reached it; no knowledge route did",
    "agreed": "",
}

# The three routes are worth running together only because they fail differently.
# Each mode is neutralised by a specific, checkable property -- not by agreement
# in general, and never by counting votes.
FAILURE_MODES = {
    "route_a_hallucinated_gene_absent_from_interval": {
        "route": "a",
        "neutralised_by": "interval_anchored",
        "why": "a family with a real interval member cannot be a gene invented outside the interval",
    },
    "route_b_post_hoc_rationalisation": {
        "route": "b",
        "neutralised_by": "interval_blind_support",
        "why": "a route that never saw the interval reached the same family, so the story was not written after looking",
    },
    "route_c_similarity_without_mechanism": {
        "route": "c",
        "neutralised_by": "mechanism_source_backed",
        "why": "a source-backed mechanism exists for this family, so the match is not similarity alone",
    },
}
COVERAGE_STATES = ("all_applicable_modes_neutralised", "partially_neutralised", "none_neutralised", "no_route_support")


def _read(path_arg: str, required: set[str], label: str) -> pd.DataFrame:
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing columns: {', '.join(missing)}")
    return frame


def _family_values(row: dict[str, object]) -> list[str]:
    values = [str(row.get(column, "")).strip() for column in FAMILY_COLUMNS if str(row.get(column, "")).strip()]
    return list(dict.fromkeys(values))


def validate_retrieval_frame(
    frame: pd.DataFrame,
    label: str = "retrieval table",
    *,
    require_family: bool = True,
) -> pd.DataFrame:
    """Validate an external retrieval table against the Phase-1 contract.

    The contract is strict on purpose: unknown columns are rejected so that
    interval-derived information cannot be smuggled into the retrieval input.
    """
    missing = sorted(set(RETRIEVAL_COLUMNS) - set(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing columns: {', '.join(missing)}")
    extra = sorted(set(frame.columns) - set(RETRIEVAL_COLUMNS))
    if extra:
        raise ValueError(f"{label} has columns outside the retrieval contract (rejecting input that could carry interval information): {', '.join(extra)}")
    if frame.empty:
        raise ValueError(f"{label} is empty")
    if frame["reference_gene_id"].astype(str).str.strip().eq("").any():
        raise ValueError(f"{label} contains an empty reference_gene_id")
    if require_family and not frame[list(FAMILY_COLUMNS)].astype(str).apply(lambda series: series.str.strip().ne("").any(), axis=1).all():
        raise ValueError(f"{label} requires family_id and/or orthogroup_id on every row")
    scores = pd.to_numeric(frame["retrieval_score"], errors="coerce")
    if scores.isna().any():
        raise ValueError(f"{label} has a non-numeric retrieval_score")
    ranks = pd.to_numeric(frame["rank"], errors="coerce")
    if ranks.isna().any() or (ranks < 1).any() or (ranks % 1 != 0).any():
        raise ValueError(f"{label} rank must be a positive integer")
    for column in META_COLUMNS:
        if frame[column].astype(str).str.strip().eq("").any():
            raise ValueError(f"{label} requires a non-empty {column} on every row (retrieval provenance)")
    if ranks.duplicated().any():
        raise ValueError(f"{label} rank values must be unique")
    return frame


def retrieval_merge(
    retrieval_args: list[str],
    out_dir_arg: str,
    top_k: int = 50,
    case_id: str = "",
) -> dict[str, object]:
    """Merge per-database retrievals without comparing cross-database scores.

    Rows lacking family/orthogroup annotation are retained in an excluded ledger
    and never consume a merged Top-K slot. Annotated rows are selected in a
    deterministic round-robin by within-database rank.
    """
    if len(retrieval_args) < 2:
        raise ValueError("retrieval-merge requires at least two --retrieval inputs")
    if top_k < 1:
        raise ValueError("top-k must be at least 1")
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    queues: list[tuple[str, Path, list[dict[str, object]]]] = []
    excluded_rows: list[dict[str, object]] = []
    database_counts: dict[str, dict[str, object]] = {}
    seen_databases: set[str] = set()
    for value in retrieval_args:
        path = runtime_path(value)
        frame = validate_retrieval_frame(
            pd.read_csv(path, sep="\t", dtype=str).fillna(""),
            f"retrieval input {path}", require_family=False,
        ).copy()
        database_versions = sorted({str(item).strip() for item in frame["database_version"] if str(item).strip()})
        if len(database_versions) != 1:
            raise ValueError(f"retrieval input {path} must contain exactly one database_version")
        database = database_versions[0]
        if database in seen_databases:
            raise ValueError(f"retrieval-merge database_version must be unique per input: {database}")
        seen_databases.add(database)
        frame["rank"] = pd.to_numeric(frame["rank"], errors="coerce").astype(int)
        frame["retrieval_score"] = pd.to_numeric(frame["retrieval_score"], errors="coerce")
        frame = frame.sort_values("rank", kind="stable")
        annotated_mask = frame[list(FAMILY_COLUMNS)].astype(str).apply(
            lambda series: series.str.strip().ne(""), axis=1,
        ).any(axis=1)
        annotated = frame[annotated_mask].to_dict("records")
        for row in frame[~annotated_mask].to_dict("records"):
            excluded_rows.append({
                **{column: row.get(column, "") for column in RETRIEVAL_COLUMNS},
                "source_database_rank": row.get("rank", ""),
                "source_input": str(path),
                "exclusion_reason": "missing_family_or_orthogroup",
            })
        queues.append((database, path, annotated))
        database_counts[database] = {
            "input_rows": int(len(frame)),
            "annotated_rows": int(annotated_mask.sum()),
            "excluded_unannotated_rows": int((~annotated_mask).sum()),
            "selected_rows": 0,
        }

    selected: list[dict[str, object]] = []
    positions = [0] * len(queues)
    seen_rows: set[tuple[str, str, str]] = set()
    while len(selected) < top_k:
        progressed = False
        for index, (database, _path, records) in enumerate(queues):
            while positions[index] < len(records):
                row = records[positions[index]]
                positions[index] += 1
                key = (
                    str(row.get("reference_gene_id", "")).casefold(),
                    canonical_family_key(row.get("family_id", "")),
                    canonical_family_key(row.get("orthogroup_id", "")),
                )
                if key in seen_rows:
                    continue
                seen_rows.add(key)
                selected.append(dict(row))
                database_counts[database]["selected_rows"] = int(database_counts[database]["selected_rows"]) + 1
                progressed = True
                break
            if len(selected) >= top_k:
                break
        if not progressed:
            break

    for rank, row in enumerate(selected, start=1):
        row["rank"] = rank
    merged = pd.DataFrame(selected, columns=RETRIEVAL_COLUMNS)
    if merged.empty:
        raise ValueError("retrieval-merge produced no annotated rows")
    validate_retrieval_frame(merged, "merged retrieval")
    merged_path = out / MERGED_NAME
    merged.to_csv(merged_path, sep="\t", index=False, lineterminator="\n")
    excluded_path = out / MERGE_EXCLUDED_NAME
    pd.DataFrame(
        excluded_rows,
        columns=RETRIEVAL_COLUMNS + ["source_database_rank", "source_input", "exclusion_reason"],
    ).to_csv(excluded_path, sep="\t", index=False, lineterminator="\n")
    status_path = out / MERGE_STATUS_NAME
    status = {
        "schema_version": SCHEMA_VERSION,
        "route": ROUTE_C,
        "case_id": case_id,
        "lane": "phenotype_only",
        "interval_information_in_model_input": False,
        "interval_arguments": [],
        "merge_policy": "deterministic_round_robin_by_within_database_rank",
        "cross_database_raw_scores_compared": False,
        "unannotated_rows_consume_top_k_slots": False,
        "top_k": top_k,
        "database_order": [database for database, _path, _records in queues],
        "selected_rows": int(len(merged)),
        "excluded_rows": int(len(excluded_rows)),
        "per_database": database_counts,
        "rank_contribution": "none",
        "inputs": {str(path): file_fingerprint(path) for _database, path, _records in queues},
        "output_fingerprints": {
            "merged": file_fingerprint(merged_path, full_hash=True),
            "excluded": file_fingerprint(excluded_path, full_hash=True),
        },
        "outputs": {"merged": str(merged_path), "excluded": str(excluded_path), "status": str(status_path)},
    }
    status_path.write_text(json.dumps(status, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"merged": str(merged_path), "excluded": str(excluded_path), "status": str(status_path),
            "selected_rows": int(len(merged)), "excluded_rows": int(len(excluded_rows))}


def _frozen_path(out: Path) -> Path:
    return out / FROZEN_NAME


_RANK_ONLY_SIDECARS = (
    (MERGE_STATUS_NAME, "merge_policy",
     "deterministic_round_robin_by_within_database_rank", "cross_database_raw_scores_compared"),
    ("route_c_retrieval_hybrid_status.json", "fusion_policy",
     "reciprocal_rank_fusion_by_within_provider_rank", "cross_provider_raw_scores_compared"),
)


def _require_rank_only_combination(retrieval_path: Path) -> None:
    """Refuse a multi-source retrieval unless a sidecar proves it was combined by rank."""
    actual = file_fingerprint(retrieval_path, full_hash=True).get("sha256")
    for status_name, policy_key, policy_value, raw_score_key in _RANK_ONLY_SIDECARS:
        status_path = retrieval_path.parent / status_name
        if not status_path.exists():
            continue
        status = json.loads(status_path.read_text(encoding="utf-8"))
        expected = status.get("output_fingerprints", {}).get("merged", {}).get("sha256")
        if (status.get(policy_key) != policy_value
                or status.get(raw_score_key) is not False
                or not expected or expected != actual):
            raise ValueError(
                f"multi-source retrieval status/fingerprint in {status_name} is missing or does not match the input"
            )
        return
    raise ValueError(
        "multi-source retrieval must be produced by retrieval-merge or retrieval-hybrid-merge; "
        "raw scores/ranks are not comparable across databases or providers"
    )


def retrieval_import(
    retrieval_arg: str,
    out_dir_arg: str,
    phenotype: str = "",
    case_id: str = "",
    phenotype_context: str | None = None,
) -> dict[str, object]:
    """Import an external retrieval TSV, validate it, and freeze it as Route C."""
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    retrieval_path = runtime_path(retrieval_arg)
    raw_frame = pd.read_csv(retrieval_path, sep="\t", dtype=str).fillna("")
    if raw_frame.empty:
        missing = sorted(set(RETRIEVAL_COLUMNS) - set(raw_frame.columns))
        if missing:
            raise ValueError(f"retrieval table is missing columns: {', '.join(missing)}")
        extra = sorted(set(raw_frame.columns) - set(RETRIEVAL_COLUMNS))
        if extra:
            raise ValueError(f"retrieval table has columns outside the retrieval contract (rejecting input that could carry interval information): {', '.join(extra)}")
        frozen = _frozen_path(out)
        raw_frame[RETRIEVAL_COLUMNS].to_csv(frozen, sep="\t", index=False, lineterminator="\n")
        status = {
            "schema_version": SCHEMA_VERSION, "route": ROUTE_C, "case_id": case_id,
            "lane": "phenotype_only",
            "interval_information_in_model_input": False,
            "blind_discovery_claim": True,
            "interval_arguments": [],
            "separation_rule": "this lane never sees the QTL interval; interval matching runs only against this frozen artifact in a separate command",
            "phenotype_used": bool(phenotype),
            "retrieval_count": 0,
            "family_count": 0,
            "retrieval_status": "not_evaluated_no_indexed_token_matched_the_phenotype",
            "models": [],
            "model_versions": [],
            "database_versions": [],
            "top10_reference_genes": [],
            "top_families": [],
            "rank_contribution": "none",
            "causal_claim": "independent discovery route only; never causal evidence",
            "structured_context": structured_context_block(
                phenotype_context, case_id=case_id, phenotype_text=phenotype,
                query_usage={
                    "usage": "stored_only",
                    "reason": "retrieval-import freezes an existing retrieval artifact and runs no query",
                },
            ),
            "inputs": {"retrieval": file_fingerprint(retrieval_path)},
            "outputs": {"frozen": str(frozen), "status": str(out / STATUS_NAME)},
        }
        (out / STATUS_NAME).write_text(json.dumps(status, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return {"frozen": str(frozen), "status": str(out / STATUS_NAME), "retrieval_count": 0,
                "family_count": 0, "rank_contribution": status["rank_contribution"]}
    frame = validate_retrieval_frame(raw_frame)
    database_versions = sorted({str(value).strip() for value in frame["database_version"] if str(value).strip()})
    if len(database_versions) > 1:
        _require_rank_only_combination(retrieval_path)
    frame["retrieval_score"] = pd.to_numeric(frame["retrieval_score"], errors="coerce")
    frame["rank"] = pd.to_numeric(frame["rank"], errors="coerce").astype(int)
    frame = frame.sort_values("rank", kind="stable").reset_index(drop=True)
    frozen = _frozen_path(out)
    frame[RETRIEVAL_COLUMNS].to_csv(frozen, sep="\t", index=False, lineterminator="\n")

    family_ids = sorted({value for row in frame.to_dict("records") for value in _family_values(row)})
    status = {
        "schema_version": SCHEMA_VERSION, "route": ROUTE_C, "case_id": case_id,
        "lane": "phenotype_only",
        "interval_information_in_model_input": False,
        "blind_discovery_claim": True,
        "interval_arguments": [],
        "separation_rule": "this lane never sees the QTL interval; interval matching runs only against this frozen artifact in a separate command",
        "phenotype_used": bool(phenotype),
        "retrieval_count": int(len(frame)),
        "family_count": len(family_ids),
        "models": sorted({str(value) for value in frame["model_name"].astype(str).str.strip() if value}),
        "model_versions": sorted({f"{row['model_name']}={row['model_version']}" for row in frame.to_dict("records")}),
        "database_versions": sorted({str(value) for value in frame["database_version"].astype(str).str.strip() if value}),
        "top10_reference_genes": [str(value) for value in frame["reference_gene_id"].head(10)],
        "top_families": family_ids[:20],
        "rank_contribution": "none",
        "causal_claim": "independent discovery route only; never causal evidence",
        "structured_context": structured_context_block(
            phenotype_context, case_id=case_id, phenotype_text=phenotype,
            query_usage={
                "usage": "stored_only",
                "reason": "retrieval-import freezes an existing retrieval artifact and runs no query",
            },
        ),
        "inputs": {"retrieval": file_fingerprint(retrieval_path)},
        "outputs": {"frozen": str(frozen), "status": str(out / STATUS_NAME)},
    }
    (out / STATUS_NAME).write_text(json.dumps(status, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"frozen": str(frozen), "status": str(out / STATUS_NAME), "retrieval_count": int(len(frame)),
            "family_count": len(family_ids), "rank_contribution": status["rank_contribution"]}


def retrieval_interval_match(
    retrieval_arg: str,
    atlas_arg: str,
    out_dir_arg: str,
    case_id: str = "",
) -> dict[str, object]:
    """Match the frozen Route C retrieval against the QTL interval."""
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    retrieval = _read(retrieval_arg, set(RETRIEVAL_COLUMNS), "frozen retrieval")
    retrieval["rank"] = pd.to_numeric(retrieval["rank"], errors="coerce").astype(int)
    atlas = _read(atlas_arg, {"gene_id"}, "interval atlas")
    if atlas.empty:
        raise ValueError("interval atlas is empty; run interval-atlas first")
    if retrieval.empty:
        match_path = out / MATCH_NAME
        pd.DataFrame(columns=["reference_gene_id", "reference_species", "family_id", "orthogroup_id",
                              "retrieval_score", "rank", "matched_interval_gene", "match_basis"]) \
            .to_csv(match_path, sep="\t", index=False, lineterminator="\n")
        status = {
            "schema_version": SCHEMA_VERSION, "route": ROUTE_C, "case_id": case_id,
            "interval_seen": "after_freeze_only",
            "retrieval_input_blind": True,
            "retrieval_frozen_before_interval": True,
            "retrieval_count": 0,
            "interval_gene_count": int(len(atlas)),
            "retrieval_rows": 0,
            "matched_retrieval_rows": 0,
            "matched_interval_genes": [],
            "matched_interval_gene_count": 0,
            "zero_match_reason": "empty_frozen_retrieval",
            "zero_match_evidence": {},
            "retrieval_status": "not_evaluated_no_indexed_token_matched_the_phenotype",
            "rank_contribution": "none",
            "causal_claim": "a similarity-ranked family with an interval member is a convergence candidate, not causal evidence",
            "inputs": {"retrieval": file_fingerprint(runtime_path(retrieval_arg)),
                       "interval_atlas": file_fingerprint(runtime_path(atlas_arg))},
            "outputs": {"interval_match": str(match_path), "status": str(out / MATCH_STATUS_NAME)},
        }
        (out / MATCH_STATUS_NAME).write_text(json.dumps(status, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return {"interval_match": str(match_path), "status": str(out / MATCH_STATUS_NAME),
                "matched_retrieval_rows": 0, "matched_interval_gene_count": 0,
                "zero_match_reason": "empty_frozen_retrieval"}
    retrieval = validate_retrieval_frame(retrieval, "frozen retrieval")
    interval_family: dict[str, set[str]] = {}
    interval_gene_family: dict[str, list[str]] = {}
    interval_orthology_reference_gene: dict[str, set[str]] = {}
    for row in atlas.to_dict("records"):
        gene_id = str(row["gene_id"])
        family = str(row.get("family_id", "")).strip()
        orthogroup = str(row.get("orthogroup", "")).strip()
        members = interval_family.setdefault(canonical_family_key(family), set()) if family else None
        if members is not None:
            members.add(gene_id)
        if orthogroup:
            interval_family.setdefault(canonical_family_key(orthogroup), set()).add(gene_id)
        interval_gene_family.setdefault(gene_id.casefold(), []).extend(value for value in (family, orthogroup) if value)
        for reference in str(row.get("orthology_reference_gene", "")).split(";"):
            reference = reference.strip()
            if reference:
                interval_orthology_reference_gene.setdefault(reference.casefold(), set()).add(gene_id)

    rows: list[dict[str, object]] = []
    match_basis_counts = {"same_family": 0, "same_orthogroup": 0, "same_reference_gene": 0,
                          "same_orthology_reference_gene": 0}
    for record in retrieval.to_dict("records"):
        family = str(record.get("family_id", "")).strip()
        orthogroup = str(record.get("orthogroup_id", "")).strip()
        matched: set[str] = set()
        basis: list[str] = []
        for key, label in ((family, "same_family"), (orthogroup, "same_orthogroup")):
            if key:
                members = interval_family.get(canonical_family_key(key), set())
                if members:
                    matched.update(members)
                    basis.append(label)
        reference = str(record.get("reference_gene_id", "")).strip()
        if reference.casefold() in interval_gene_family:
            matched.add(reference)
            basis.append("same_reference_gene")
        orthology_members = interval_orthology_reference_gene.get(reference.casefold(), set())
        if orthology_members:
            matched.update(orthology_members)
            basis.append("same_orthology_reference_gene")
        for label in set(basis):
            match_basis_counts[label] += 1
        rows.append({
            "reference_gene_id": reference,
            "reference_species": str(record.get("reference_species", "")),
            "family_id": family, "orthogroup_id": orthogroup,
            "retrieval_score": record.get("retrieval_score", ""),
            "rank": record["rank"],
            "matched_interval_gene": ";".join(sorted(matched)),
            "match_basis": ";".join(sorted(set(basis))),
        })
    matched_frame = pd.DataFrame(rows)
    match_path = out / MATCH_NAME
    matched_frame.to_csv(match_path, sep="\t", index=False, lineterminator="\n")
    matched = matched_frame["matched_interval_gene"].ne("")
    matched_family_rows = int(matched.sum())
    matched_genes = sorted({gene for value in matched_frame.loc[matched, "matched_interval_gene"] for gene in str(value).split(";")})
    zero_match_evidence = {
        "atlas_family_or_orthogroup_key_count": len(interval_family),
        "atlas_orthology_reference_gene_key_count": len(interval_orthology_reference_gene),
        "retrieval_family_or_orthogroup_key_count": len({canonical_family_key(str(record.get(column, "")).strip())
                                                           for record in retrieval.to_dict("records")
                                                           for column in ("family_id", "orthogroup_id")
                                                           if str(record.get(column, "")).strip()}),
        "retrieval_reference_gene_id_count": int(retrieval["reference_gene_id"].astype(str).str.strip().ne("").sum()),
        **match_basis_counts,
    }
    zero_match_reason = None
    if matched_family_rows == 0:
        zero_match_reason = (
            "interval_atlas_join_keys_unavailable"
            if not interval_family and not interval_orthology_reference_gene
            else "namespace_mismatch_suspected"
        )

    status = {
        "schema_version": SCHEMA_VERSION, "route": ROUTE_C, "case_id": case_id,
        "interval_seen": "after_freeze_only",
        "retrieval_input_blind": True,
        "retrieval_frozen_before_interval": True,
        "interval_gene_count": int(len(atlas)),
        "retrieval_rows": int(len(matched_frame)),
        "matched_retrieval_rows": matched_family_rows,
        "matched_interval_genes": matched_genes,
        "matched_interval_gene_count": len(matched_genes),
        "zero_match_reason": zero_match_reason,
        "zero_match_evidence": zero_match_evidence,
        "rank_contribution": "none",
        "causal_claim": "a similarity-ranked family with an interval member is a convergence candidate, not causal evidence",
        "inputs": {"retrieval": file_fingerprint(runtime_path(retrieval_arg)), "interval_atlas": file_fingerprint(runtime_path(atlas_arg))},
        "outputs": {"interval_match": str(match_path), "status": str(out / MATCH_STATUS_NAME)},
    }
    (out / MATCH_STATUS_NAME).write_text(json.dumps(status, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"interval_match": str(match_path), "status": str(out / MATCH_STATUS_NAME),
            "matched_retrieval_rows": matched_family_rows, "matched_interval_gene_count": len(matched_genes),
            "zero_match_reason": zero_match_reason}


def _route_families(frame: pd.DataFrame, column: str) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    for row in frame.to_dict("records"):
        for family in split_respecting_brackets(str(row.get(column, "")), separator=";"):
            family = family.strip()
            if family:
                result.setdefault(canonical_family_key(family), []).append(str(row["hypothesis_id"]))
    return result


def _c_families(frame: pd.DataFrame) -> dict[str, dict[str, object]]:
    """Route C families keyed casefolded; family_id and orthogroup_id are both kept."""
    result: dict[str, dict[str, object]] = {}
    for row in frame.to_dict("records"):
        for column in FAMILY_COLUMNS:
            value = str(row.get(column, "")).strip()
            if not value:
                continue
            record = result.setdefault(canonical_family_key(value), {
                "display": value, "reference_genes": [], "ranks": [], "sources": set(),
                "family_id": str(row.get("family_id", "") or row.get("orthogroup_id", "")).strip(),
            })
            record["display"] = value
            record["reference_genes"].append(str(row["reference_gene_id"]))
            record["ranks"].append(int(row["rank"]))
            database = str(row.get("database_version", "")).strip()
            if database:
                record["sources"].add(database)  # type: ignore[union-attr]
    return result


def _route_attributes(frame: pd.DataFrame) -> dict[str, dict[str, object]]:
    """Per family: name, interval members, source ids and mechanical family id."""
    result: dict[str, dict[str, object]] = {}
    for row in frame.to_dict("records"):
        raw_family_id = str(row.get("family_id", "")).strip()
        for value in split_respecting_brackets(str(row.get("gene_family", "")), separator=";"):
            if not value.strip():
                continue
            key = canonical_family_key(value)
            family_id_val = raw_family_id or alias_canonical(value) or value.strip()
            record = result.setdefault(key,
                                       {"display": value.strip(), "members": set(), "sources": set(), "family_id": family_id_val})
            record["members"].update(item for item in split_respecting_brackets(str(row.get("member_interval_genes", "")), separator=";") if item.strip())  # type: ignore[union-attr]
            record["sources"].update(item for item in split_respecting_brackets(str(row.get("source_backing", "")), separator=";") if item.strip())  # type: ignore[union-attr]
    return result


def _expanded_ab_aliases(
    families_a: dict[str, list[str]],
    families_b: dict[str, list[str]],
    attributes_a: dict[str, dict[str, object]],
    attributes_b: dict[str, dict[str, object]],
) -> dict[str, str]:
    """Return route-B strict keys that can safely alias one route-A strict key.

    Expanded matching is a bounded fallback for issue #500, not a replacement
    for strict family identity. Only strict-unmatched A/B keys are considered,
    and an alias is accepted only when the expanded relation is mutually unique
    (one A key matches exactly one B key and that B key matches exactly one A
    key). This prevents a broad token such as a paralog-family stem from silently
    collapsing several independently named families.
    """
    unmatched_a = [key for key in families_a if key not in families_b]
    unmatched_b = [key for key in families_b if key not in families_a]
    a_to_b: dict[str, set[str]] = {key: set() for key in unmatched_a}
    b_to_a: dict[str, set[str]] = {key: set() for key in unmatched_b}
    for a_key in unmatched_a:
        a_display = str(attributes_a.get(a_key, {}).get("display", "")).strip()
        if not a_display:
            continue
        for b_key in unmatched_b:
            b_display = str(attributes_b.get(b_key, {}).get("display", "")).strip()
            if not b_display:
                continue
            if (families_match_expanded(a_display, b_display)
                    or families_match_expanded(b_display, a_display)):
                a_to_b[a_key].add(b_key)
                b_to_a[b_key].add(a_key)
    return {
        next(iter(b_matches)): a_key
        for a_key, b_matches in a_to_b.items()
        if len(b_matches) == 1 and len(b_to_a[next(iter(b_matches))]) == 1
    }


def _remap_list_values(mapping: dict[str, list[str]], aliases: dict[str, str]) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    for key, values in mapping.items():
        target = aliases.get(key, key)
        result.setdefault(target, []).extend(values)
    return result


def _remap_attribute_records(
    mapping: dict[str, dict[str, object]], aliases: dict[str, str]
) -> dict[str, dict[str, object]]:
    result: dict[str, dict[str, object]] = {}
    for key, record in mapping.items():
        target = aliases.get(key, key)
        if target not in result:
            result[target] = {
                "display": record.get("display", ""),
                "members": set(record.get("members", set())),
                "sources": set(record.get("sources", set())),
                "family_id": record.get("family_id", ""),
            }
        else:
            result[target]["members"].update(record.get("members", set()))  # type: ignore[union-attr]
            result[target]["sources"].update(record.get("sources", set()))  # type: ignore[union-attr]
            if not str(result[target].get("family_id", "")).strip():
                result[target]["family_id"] = record.get("family_id", "")
    return result


def _remap_c_records(
    mapping: dict[str, dict[str, object]], aliases: dict[str, str]
) -> dict[str, dict[str, object]]:
    result: dict[str, dict[str, object]] = {}
    for key, record in mapping.items():
        target = aliases.get(key, key)
        if target not in result:
            result[target] = {
                "display": record.get("display", ""),
                "reference_genes": list(record.get("reference_genes", [])),
                "ranks": list(record.get("ranks", [])),
                "sources": set(record.get("sources", set())),
                "family_id": record.get("family_id", ""),
            }
        else:
            result[target]["reference_genes"].extend(record.get("reference_genes", []))  # type: ignore[union-attr]
            result[target]["ranks"].extend(record.get("ranks", []))  # type: ignore[union-attr]
            result[target]["sources"].update(record.get("sources", set()))  # type: ignore[union-attr]
    return result


def _route_source_ids(attributes_record: dict[str, object]) -> set[str]:
    """Real reference-record ids behind one knowledge route."""
    return {token for token in attributes_record.get("sources", set()) if ":" in token}


def _interval_match_families(path_arg: str | None) -> dict[str, set[str]]:
    """Route C families that reached an interval gene, from the Phase-2 artifact."""
    if not path_arg:
        return {}
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    if "matched_interval_gene" not in frame.columns:
        raise ValueError("route C interval-match table is missing matched_interval_gene; use route_c_interval_match.tsv")
    result: dict[str, set[str]] = {}
    for row in frame.to_dict("records"):
        genes = {value for value in str(row.get("matched_interval_gene", "")).split(";") if value.strip()}
        if not genes:
            continue
        for column in FAMILY_COLUMNS:
            value = str(row.get(column, "")).strip()
            if value:
                result.setdefault(canonical_family_key(value), set()).update(genes)
    return result


def _failure_mode_account(
    in_a: bool, in_b: bool, in_c: bool,
    interval_members: set[str], anchor_basis: list[str],
    source_ids: set[str],
) -> dict[str, object]:
    """Which of the three characteristic errors this family actually escapes."""
    interval_blind_support = in_a or in_c
    checks = {
        "interval_anchored": bool(interval_members),
        "interval_blind_support": bool(interval_blind_support and in_b),
        "mechanism_source_backed": bool(source_ids),
    }
    applicable = [name for name, mode in FAILURE_MODES.items()
                  if {"a": in_a, "b": in_b, "c": in_c}[str(mode["route"])]]
    neutralised = [name for name in applicable if checks[str(FAILURE_MODES[name]["neutralised_by"])]]
    residual = [name for name in applicable if name not in neutralised]
    if not applicable:
        coverage = "no_route_support"
    elif not residual:
        coverage = "all_applicable_modes_neutralised"
    elif neutralised:
        coverage = "partially_neutralised"
    else:
        coverage = "none_neutralised"
    return {
        "applicable_failure_modes": ";".join(applicable),
        "neutralised_failure_modes": ";".join(neutralised),
        "residual_failure_modes": ";".join(residual),
        "failure_mode_coverage": coverage,
        "interval_anchored": checks["interval_anchored"],
        "interval_anchor_basis": ";".join(sorted(set(anchor_basis))),
        "interval_blind_support": checks["interval_blind_support"],
        "mechanism_source_backed": checks["mechanism_source_backed"],
        "mechanism_source_ids": ";".join(sorted(source_ids)),
    }


def _literature_rows(literature_arg: str | None) -> dict[str, set[int]]:
    if not literature_arg:
        return {}
    frame = pd.read_csv(runtime_path(literature_arg), sep="\t", dtype=str).fillna("")
    if "interval_gene_ids" not in frame.columns:
        raise ValueError("reverse-literature table is missing interval_gene_ids; use reverse_literature.tsv")
    by_gene: dict[str, set[int]] = {}
    for index, row in frame.iterrows():
        for gene in str(row["interval_gene_ids"]).split(";"):
            if gene.strip():
                by_gene.setdefault(gene.strip().casefold(), set()).add(index)
    return by_gene


def _literature_status_counts(literature_status_arg: str | None) -> dict[str, int]:
    if not literature_status_arg:
        return {}
    frame = pd.read_csv(runtime_path(literature_status_arg), sep="\t", dtype=str).fillna("")
    if "query_status" not in frame.columns:
        raise ValueError("reverse-literature status table is missing query_status; use reverse_literature_status.tsv")
    values = frame["query_status"].astype(str)
    return {
        "queried": int((values == "retrieved").sum()),
        "no_hit": int((values == "no_hit").sum()),
        "not_evaluated_offline": int((values == "not_evaluated_offline").sum()),
    }


def triple_route_compare(
    route_a_arg: str,
    route_b_arg: str,
    route_c_arg: str,
    out_dir_arg: str,
    case_id: str = "",
    route_c_match_arg: str | None = None,
    literature_arg: str | None = None,
    literature_status_arg: str | None = None,
) -> dict[str, object]:
    """Compare all three independently generated routes; convergence adds no score."""
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    route_a = _read(route_a_arg, {"hypothesis_id"}, "route A hypothesis")
    route_b = _read(route_b_arg, {"hypothesis_id"}, "route B hypothesis")
    route_c = _read(route_c_arg, set(RETRIEVAL_COLUMNS), "route C retrieval")
    if not route_c.empty:
        route_c = validate_retrieval_frame(route_c, "route C retrieval")
        route_c["rank"] = pd.to_numeric(route_c["rank"], errors="coerce").astype(int)

    literature_by_gene = _literature_rows(literature_arg)
    literature_status_counts = _literature_status_counts(literature_status_arg)

    families_a = _route_families(route_a, "gene_family")
    families_b = _route_families(route_b, "gene_family")
    attributes_a = _route_attributes(route_a)
    attributes_b = _route_attributes(route_b)

    # Issue #500: strict identity remains the first tier. Only strict-unmatched
    # A/B keys may be joined by the existing expanded matcher, and only when the
    # relation is mutually unique. The route-B key is then aliased to the route-A
    # key for this comparison run. The alias is recorded and never treated as a
    # strict convergence.
    expanded_aliases = _expanded_ab_aliases(families_a, families_b, attributes_a, attributes_b)
    families_b = _remap_list_values(families_b, expanded_aliases)
    attributes_b = _remap_attribute_records(attributes_b, expanded_aliases)
    family_match_basis = {target: "expanded_symbol_token" for target in expanded_aliases.values()}

    families_c = _remap_c_records(_c_families(route_c), expanded_aliases)
    c_interval_raw = _interval_match_families(route_c_match_arg)
    c_interval: dict[str, set[str]] = {}
    for key, genes in c_interval_raw.items():
        c_interval.setdefault(expanded_aliases.get(key, key), set()).update(genes)

    interval_members = {canonical_family_key(value): str(row.get("member_interval_genes", ""))
                        for row in route_a.to_dict("records") for value in split_respecting_brackets(str(row.get("gene_family", "")), separator=";") if value.strip()}
    for row in route_b.to_dict("records"):
        for value in split_respecting_brackets(str(row.get("gene_family", "")), separator=";"):
            if value.strip():
                strict_key = canonical_family_key(value)
                key = expanded_aliases.get(strict_key, strict_key)
                interval_members.setdefault(key, str(row.get("member_interval_genes", "")))

    all_families = sorted(set(families_a) | set(families_b) | set(families_c), key=str.casefold)
    rows: list[dict[str, object]] = []
    comparison_family_ids: dict[str, str] = {}
    for family in all_families:
        in_a = family in families_a
        in_b = family in families_b
        in_c = family in families_c
        match_basis = family_match_basis.get(family, "strict_key")
        if in_a and in_b and match_basis == "expanded_symbol_token":
            status = "a_b_expanded_c_supported" if in_c else "a_b_expanded_match"
        elif in_a and in_b and in_c:
            status = "a_b_c_convergent"
        elif in_a and in_c:
            status = "a_c_convergent"
        elif in_b and in_c:
            status = "b_c_convergent"
        elif in_a and in_b:
            status = "a_b_convergent"
        elif in_c:
            status = "c_only_supported"
        elif in_a:
            status = "phenotype_only_supported"
        else:
            status = "interval_only_supported"
        record_c = families_c.get(family, {})
        members: set[str] = set()
        anchor_basis: list[str] = []
        for label, source in (("route_a", attributes_a), ("route_b", attributes_b)):
            found = source.get(family, {}).get("members", set())
            if found:
                members.update(found)  # type: ignore[arg-type]
                anchor_basis.append(label)
        if c_interval.get(family):
            members.update(c_interval[family])
            anchor_basis.append("route_c_interval_match")
        for value in split_respecting_brackets(str(interval_members.get(family, "")), separator=";"):
            if value.strip():
                members.add(value.strip())
        sources_a = set(attributes_a.get(family, {}).get("sources", set()))  # type: ignore[arg-type]
        sources_b = set(attributes_b.get(family, {}).get("sources", set()))  # type: ignore[arg-type]
        sources = sources_a | sources_b
        account = _failure_mode_account(in_a, in_b, in_c, members, anchor_basis, sources)
        shared_sources = sorted(_route_source_ids(attributes_a.get(family, {}))
                                & _route_source_ids(attributes_b.get(family, {})))
        if in_a and in_b:
            agreement_state = "agreed_a_b"
            if match_basis == "expanded_symbol_token":
                agreement_strength = "expanded_symbol_match"
            else:
                agreement_strength = "shared_sources" if shared_sources else "disjoint_sources"
            disagreement_reason = ""
        elif in_a:
            agreement_state = "a_only"
            disagreement_reason = ("a_only_interval_annotation_thin" if members
                                   else "a_only_not_in_interval")
            agreement_strength = "single_route"
        elif in_b:
            agreement_state = "b_only"
            disagreement_reason = "b_only_no_phenotype_link_in_literature"
            agreement_strength = "single_route"
        else:
            agreement_state = "c_only"
            disagreement_reason = "c_only_similarity_without_knowledge"
            agreement_strength = "experimental_route_only"
        account.update({
            "agreement_state": agreement_state, "agreement_strength": agreement_strength,
            "shared_source_ids": ";".join(shared_sources), "disagreement_reason": disagreement_reason,
        })
        source_ids_a = _route_source_ids(attributes_a.get(family, {}))
        source_ids_b = _route_source_ids(attributes_b.get(family, {}))
        source_ids_c = set(record_c.get("sources", set()))  # type: ignore[arg-type]
        source_sets = {"a": source_ids_a, "b": source_ids_b, "c": source_ids_c}
        supporting = [route for route in ("a", "b", "c")
                      if {"a": in_a, "b": in_b, "c": in_c}[route]]
        shared_across_routes = any(
            source_sets[left] & source_sets[right]
            for left, right in (("a", "b"), ("a", "c"), ("b", "c"))
        )
        if len(supporting) < 2:
            source_dependence = "not_applicable"
        elif shared_across_routes:
            source_dependence = "shared_source_dependent"
        elif all(source_sets[route] for route in supporting):
            source_dependence = "source_independent"
        else:
            source_dependence = "source_ids_incomplete"
        independent_source_count = len({token for route in supporting for token in source_sets[route]})
        notes = ("convergence is agreement between independent routes, not additional causal evidence; "
                 "a neutralised failure mode lowers that route's characteristic risk without proving the family")
        if match_basis == "expanded_symbol_token":
            notes += ("; A/B family identity used the expanded symbol-token matcher after strict-key mismatch; "
                      "this recovers a semantic correspondence but is non-strict and cannot become strict_convergent")
        if status in CONVERGENT_STATES and source_dependence == "shared_source_dependent":
            notes += (f"; shared-source-dependent: {len(supporting)} routes converged on "
                      f"{independent_source_count} independent source(s) -- shared source ids are not "
                      f"counted as independent evidence")

        if in_a and in_b:
            four_state = "a_and_b_supported"
        elif in_a:
            four_state = "a_only_supported"
        elif in_b:
            four_state = "b_only_supported"
        else:
            four_state = "neither_supported"
        if literature_arg:
            literature_row_indices: set[int] = set()
            for gene in members:
                literature_row_indices |= literature_by_gene.get(gene.casefold(), set())
            literature_count = len(literature_row_indices)
        else:
            literature_count = ""
        display = str(record_c.get("display", "")
                      or attributes_a.get(family, {}).get("display", "")
                      or attributes_b.get(family, {}).get("display", "")
                      or family)
        display_a = str(attributes_a.get(family, {}).get("display", "")).strip()
        mech_family_id = (
            str(attributes_b.get(family, {}).get("family_id", "")).strip()
            or str(record_c.get("family_id", "") or record_c.get("orthogroup_id", "")).strip()
            or (alias_canonical(display_a) or "")
            or (display_a if display_a and " " not in display_a and "(" not in display_a and ";" not in display_a else "")
        )
        cmp_id = f"cmp_{len(rows):03d}"
        comparison_family_ids[cmp_id] = mech_family_id
        rows.append({
            "comparison_id": cmp_id,
            "family_id": display,
            "route_a_hypothesis_ids": ";".join(families_a.get(family, [])),
            "route_b_hypothesis_ids": ";".join(families_b.get(family, [])),
            "route_c_reference_genes": ";".join(record_c.get("reference_genes", [])),
            "route_c_retrieval_ranks": ";".join(str(rank) for rank in sorted(record_c.get("ranks", []))),
            "matched_interval_genes": ";".join(sorted(members)),
            "convergence_status": status,
            "family_match_basis": match_basis,
            **account,
            "route_a_source_ids": ";".join(sorted(source_ids_a)),
            "route_b_source_ids": ";".join(sorted(source_ids_b)),
            "route_c_source_ids": ";".join(sorted(source_ids_c)),
            "source_dependence_status": source_dependence,
            "rank_contribution": "none",
            "four_state": four_state, "literature_link_candidate_count": literature_count,
            "notes": notes,
        })
    comparison = pd.DataFrame(rows, columns=COMPARISON_COLUMNS)

    counts = {state: int((comparison["convergence_status"] == state).sum()) for state in CONVERGENCE_STATES} if len(comparison) else {state: 0 for state in CONVERGENCE_STATES}
    n_a, n_b, n_c = len(families_a), len(families_b), len(families_c)
    if not (n_a or n_b or n_c):
        overall = "all_weak"
    elif counts["a_b_c_convergent"]:
        overall = "a_b_c_convergent"
    elif counts["a_c_convergent"]:
        overall = "a_c_convergent"
    elif counts["b_c_convergent"]:
        overall = "b_c_convergent"
    elif counts["a_b_convergent"]:
        overall = "a_b_convergent"
    elif counts["a_b_expanded_c_supported"]:
        overall = "a_b_expanded_c_supported"
    elif counts["a_b_expanded_match"]:
        overall = "a_b_expanded_match"
    elif n_c and not (n_a or n_b):
        overall = "c_only_supported"
    elif n_c and (n_a or n_b):
        overall = "all_discordant"
    elif n_a and n_b:
        overall = "all_discordant"
    elif n_a:
        overall = "phenotype_only_supported"
    else:
        overall = "interval_only_supported"

    review_rows: list[dict[str, str]] = []
    if overall == "all_discordant":
        review_rows.append({"reason": "routes_discordant", "target": case_id or "case",
                            "required_next_action": "3経路が別のmechanismを指している。どれも捨てず、識別できる実験を選ぶ",
                            "allowed_choices": "両方を保持して次実験で識別 | 追加文献検証 | retrieval再実行"})
    if overall == "all_weak":
        review_rows.append({"reason": "all_routes_weak", "target": case_id or "case",
                            "required_next_action": "どのrouteも仮説を出せていない。追加mapping、annotation rescue、またはretrieval DB/queryの見直しへ戻る",
                            "allowed_choices": "追加mapping | annotation rescue | retrieval DB/query見直し | abstain"})
    if counts["c_only_supported"]:
        review_rows.append({"reason": "c_only_supported", "target": case_id or "case",
                            "required_next_action": "retrievalが単独でfamilyに到達した。causal evidenceではなく、区間照合と実験で確認する",
                            "allowed_choices": "retrieval-interval-match | 文献検証 | 識別実験"})
    reviews = pd.DataFrame(review_rows, columns=REVIEW_COLUMNS) if review_rows else pd.DataFrame(columns=REVIEW_COLUMNS)

    disagreement_rows = [{
        "family_id": str(row["family_id"]), "agreement_state": str(row["agreement_state"]),
        "disagreement_reason": str(row["disagreement_reason"]),
        "routes_supporting": ";".join(name for name, present in (
            ("a", bool(row["route_a_hypothesis_ids"])), ("b", bool(row["route_b_hypothesis_ids"])),
            ("c", bool(row["route_c_reference_genes"]))) if present),
        "matched_interval_genes": str(row["matched_interval_genes"]),
        "required_next_action": {
            "a_only_not_in_interval": "annotation-rescue で配列レベルの不在を確認する。未確認のまま候補を捨てない",
            "a_only_interval_annotation_thin": "区間側の注釈が薄い。orthology/domain provider を足す",
            "b_only_no_phenotype_link_in_literature": "文献にこの表現型との結びつきが無い。新規性の候補として実験で識別する",
            "c_only_similarity_without_knowledge": "実験的routeのみ。順位には効かせず、知識routeで裏が取れるか確認する",
        }.get(str(row["disagreement_reason"]), ""),
    } for row in comparison.to_dict("records") if str(row["disagreement_reason"])] if len(comparison) else []
    disagreements = pd.DataFrame(disagreement_rows, columns=DISAGREEMENT_COLUMNS)

    paths = {"triple_route_hypothesis_comparison": out / "triple_route_hypothesis_comparison.tsv",
             "triple_route_disagreements": out / "triple_route_disagreements.tsv",
             "triple_route_review_required": out / "triple_route_review_required.tsv"}
    comparison.to_csv(paths["triple_route_hypothesis_comparison"], sep="\t", index=False, lineterminator="\n")
    disagreements.to_csv(paths["triple_route_disagreements"], sep="\t", index=False, lineterminator="\n")
    reviews.to_csv(paths["triple_route_review_required"], sep="\t", index=False, lineterminator="\n")

    handoff_rows = [{
        "hypothesis_id": str(row["comparison_id"]),
        "family_id": str(comparison_family_ids.get(str(row["comparison_id"]), "")),
        "gene_family": str(row["family_id"]),
        "member_interval_genes": str(row["matched_interval_genes"]),
    } for row in comparison.to_dict("records") if str(row["matched_interval_genes"])]
    if handoff_rows:
        handoff = pd.DataFrame(handoff_rows, columns=FAMILY_TO_COPY_INPUT_COLUMNS)
        paths["family_to_copy_input"] = out / "family_to_copy_input.tsv"
        handoff.to_csv(paths["family_to_copy_input"], sep="\t", index=False, lineterminator="\n")

    if len(comparison):
        convergent_rows = comparison[comparison["convergence_status"].isin(CONVERGENT_STATES)].to_dict("records")
    else:
        convergent_rows = []
    convergence_accounts: list[dict[str, object]] = []
    for row in convergent_rows:
        route_sets = [
            {token for token in str(row[f"route_{route}_source_ids"]).split(";") if token}
            for route in ("a", "b", "c")
        ]
        distinct = sorted(set().union(*route_sets))
        shared = [token for token in distinct if sum(token in ids for ids in route_sets) >= 2]
        convergence_accounts.append({
            "family_id": str(row["family_id"]),
            "routes_converged": int(sum(bool(str(row[column])) for column in (
                "route_a_hypothesis_ids", "route_b_hypothesis_ids", "route_c_reference_genes"))),
            "independent_sources": len(distinct),
            "shared_source_ids": shared,
            "source_dependence_status": str(row["source_dependence_status"]),
        })
    dependent_convergences = [str(account["family_id"]) for account in convergence_accounts
                              if account["source_dependence_status"] == "shared_source_dependent"]
    independent_convergences = [str(account["family_id"]) for account in convergence_accounts
                                if account["source_dependence_status"] == "source_independent"]

    expanded_pairs = [{
        "route_b_strict_key": b_key,
        "comparison_key": a_key,
        "route_a_display": str(attributes_a.get(a_key, {}).get("display", "")),
        "route_b_display": str(_route_attributes(route_b).get(b_key, {}).get("display", "")),
        "match_basis": "expanded_symbol_token",
    } for b_key, a_key in sorted(expanded_aliases.items())]

    convergence = {
        "schema_version": SCHEMA_VERSION, "case_id": case_id, "overall_status": overall,
        "route_a_families": n_a, "route_b_families": n_b, "route_c_families": n_c,
        "status_counts": {key: value for key, value in counts.items() if value},
        "family_matching": {
            "basis_counts": {
                basis: int((comparison["family_match_basis"] == basis).sum())
                for basis in sorted({str(value) for value in comparison["family_match_basis"]})
            } if len(comparison) else {},
            "expanded_match_count": len(expanded_pairs),
            "expanded_pairs": expanded_pairs,
            "rule": (
                "strict canonical-family keys are tried first. Only strict-unmatched A/B families may use "
                "the existing expanded matcher, and only a mutually unique expanded match is joined. "
                "expanded_symbol_token is recorded separately, adds no score, and is never a strict convergence"
            ),
        },
        "route_risk": {
            "route_a_failure_mode": "hallucinating a gene that is not in the interval",
            "route_b_failure_mode": "rationalising a gene that is in the interval",
            "route_c_failure_mode": "retrieval ranks by representation similarity, not mechanism: a plausible-sounding family can rank high yet be unrelated to the interval",
            "why_three_routes": "A and B fail in mirror-image ways; C fails differently (similarity is not mechanism), so a family reached independently by several routes has each route's distinct failure mode checked off",
            "residual_risk": "convergence lowers risk but is not proof and adds no score",
        },
        "failure_mode_cancellation": {
            "modes": FAILURE_MODES,
            "coverage_vocabulary": list(COVERAGE_STATES),
            "coverage_counts": {state: int((comparison["failure_mode_coverage"] == state).sum())
                                for state in COVERAGE_STATES} if len(comparison) else {},
            "fully_neutralised_families": sorted(
                str(row["family_id"]) for row in comparison.to_dict("records")
                if str(row["failure_mode_coverage"]) == "all_applicable_modes_neutralised"
            ) if len(comparison) else [],
            "route_c_interval_match_supplied": bool(route_c_match_arg),
            "rule": (
                "only the failure modes of the routes that support a family are applicable, and each is "
                "neutralised by its own specific check (interval anchoring / interval-blind support / "
                "source-backed mechanism), never by counting agreeing routes"
            ),
            "residual": (
                "a fully neutralised family has each supporting route's characteristic error checked off; "
                "that lowers risk, is not proof, and adds no score. Convergence on a shared source id is "
                "recorded as source-dependence: N routes citing the same reference provide one independent "
                "source, not N"
            ),
        },
        "disagreement": {
            "reasons": DISAGREEMENT_REASONS,
            "counts": {reason: int((comparison["disagreement_reason"] == reason).sum())
                       for reason in sorted({str(value) for value in comparison["disagreement_reason"] if str(value)})}
            if len(comparison) else {},
            "agreement_strength_counts": {state: int((comparison["agreement_strength"] == state).sum())
                                          for state in sorted({str(value) for value in comparison["agreement_strength"]})}
            if len(comparison) else {},
            "families_agreed_on_shared_sources": sorted(
                str(row["family_id"]) for row in comparison.to_dict("records")
                if str(row["agreement_strength"]) == "shared_sources") if len(comparison) else [],
            "rule": (
                "Routes A and B traverse the same literature, so agreement is the expected case and carries little "
                "information; the disagreements are the primary output. An agreement whose two routes cite the same "
                "source ids is downgraded to shared_sources: one source quoted twice is not two sources"
            ),
        },
        "source_dependence": {
            "status_counts": {state: int((comparison["source_dependence_status"] == state).sum())
                              for state in SOURCE_DEPENDENCE_STATUSES} if len(comparison) else {},
            "dependent_convergence_count": len(dependent_convergences),
            "independent_convergence_count": len(independent_convergences),
            "families_dependent_on_shared_sources": sorted(dependent_convergences),
            "converged_families_by_independent_source_count": {
                str(count): int(sum(1 for account in convergence_accounts
                                    if account["independent_sources"] == count))
                for count in sorted({account["independent_sources"] for account in convergence_accounts})
            },
            "convergence_source_accounting": convergence_accounts,
            "rule": (
                "N routes that converged on the same source id provide one independent source, not N: shared "
                "source ids are not counted as independent evidence. A convergent family whose routes cite "
                "disjoint real record ids is source_independent; one whose routes share a common id is "
                "shared_source_dependent. Provider or literature overlap downgrades the *reported independence* "
                "of a convergence only -- it never changes convergence_status for strict matches and never changes "
                "any native score. Expanded family matches are a separate non-strict state regardless of source independence"
            ),
        },
        "route_c_status": ROUTE_C_STATUS,
        "route_c_rule": "Route C is experimental: it may inform review, and never affects candidate ordering",
        "independence": (
            "procedural: the three lanes were generated without reading each other, and this comparison is the only "
            "place they meet. A and B are not informationally independent -- they read the same literature. Shared "
            "source ids are counted once: N routes that cite the same reference record provide one independent "
            "source, not N"
        ),
        "benchmark_rule": "Route C retrieval is phenotype-only and interval-blind; interval-first results are interval-informed and must be scored in their own lane, never pooled with blind phenotype-only recall",
        "safety": [
            "convergence is agreement between routes, not extra causal evidence, and adds no score",
            "Route C matches are recorded as convergence state only; they never enter a causal score",
            "expanded family-label matching is a separate non-strict tier and never upgrades to strict convergence",
            "a discordant result keeps all mechanisms alive instead of picking one",
            "all_weak routes back to additional mapping / annotation rescue / retrieval review rather than forcing a candidate",
            "a source id shared by several routes is one independent source, not several: convergence on the same reference record is not independent corroboration, and shared source ids are not counted as independent evidence",
        ],
        "inputs": {
            "route_a": file_fingerprint(runtime_path(route_a_arg)),
            "route_b": file_fingerprint(runtime_path(route_b_arg)),
            "route_c": file_fingerprint(runtime_path(route_c_arg)),
            "route_c_interval_match": file_fingerprint(runtime_path(route_c_match_arg)) if route_c_match_arg else None,
        },
        "outputs": {name: str(path) for name, path in sorted(paths.items())},
    }
    if literature_arg:
        convergence["route_b_literature"] = {
            "literature": str(runtime_path(literature_arg)),
            "status": str(runtime_path(literature_status_arg)) if literature_status_arg else None,
            "queried": literature_status_counts.get("queried", 0),
            "no_hit": literature_status_counts.get("no_hit", 0),
            "not_evaluated_offline": literature_status_counts.get("not_evaluated_offline", 0),
            "note": "a returned article is a link candidate, not evidence: retrieval is not relevance, "
                    "and these counts are never added to a score",
        }
    convergence_path = out / "triple_route_convergence.json"
    convergence_path.write_text(json.dumps(convergence, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    cancellation = convergence["failure_mode_cancellation"]
    assert isinstance(cancellation, dict)
    return {"convergence": str(convergence_path), "overall_status": overall,
            "comparison_rows": int(len(comparison)),
            "fully_neutralised_families": cancellation["fully_neutralised_families"],
            "coverage_counts": cancellation["coverage_counts"],
            **{name: str(path) for name, path in paths.items()}}
