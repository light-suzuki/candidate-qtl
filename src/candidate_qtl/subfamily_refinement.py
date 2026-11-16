"""Deterministic machine-verified subfamily refinement for an over-budget set.

The canonical Wet budget contract says an over-budget supported-member union
must go through *bounded evidence refinement* and re-gate, never a forced
Top-20 truncation. The adaptive-synthesis seam already exists, but its
verification table had no machine producer, so an automated run could only stop
as ``review_required_over_budget``. This module is that producer for one narrow,
fully deterministic evidence kind:

* the Route A hypothesis projected a verified reference target with an expected
  subfamily into the interval;
* ``family-to-copy`` recorded a member's curated subfamily assignment
  (``subfamily_source=curated_map``) and scored it target-relative;
* the member's curated subfamily differs from every expected target subfamily.

Only then does the member's Route-A-derived support count as falsified. An
orthogroup-proxy mismatch is never sufficient -- orthogroup granularity cannot
separate subfamilies -- and a pending-review target never falsifies anything.
Route-B-supported members are out of scope entirely: their support is
independent of the Route A projection, so a cross-subfamily verdict against the
Route A target says nothing about it.

The output is a verification table in exactly the vocabulary that
:func:`adaptive_synthesis.apply_verified_refinement` consumes; this module never
removes a candidate itself.
"""

from __future__ import annotations

import hashlib
import json

import pandas as pd

from .core import runtime_path

SCHEMA_VERSION = "1.0"
CONTRACT = "subfamily-refinement-verification-1.0.0"

VERIFICATION_TABLE_COLUMNS = ["candidate_gene_id", "verdict", "verification_source"]

# Same trust boundary family_to_copy applies to targets: a
# pending_review / LLM-only target may inform axes but never drives selection.
VERIFIED_TARGET_STATUSES = frozenset({"verified", "source_backed", "curated"})

CURATED_SOURCE = "curated_map"


def _sha256_path(path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _candidate_id_column(frame: pd.DataFrame) -> str:
    for column in ("candidate_gene_id", "gene_id"):
        if column in frame.columns:
            return column
    raise ValueError("candidate input requires candidate_gene_id or gene_id")


def _route_a_only(row: dict[str, str]) -> bool:
    """A candidate whose only support is the Route A family projection."""
    return str(row.get("route_support", "")).strip().casefold() == "a"


def build_subfamily_refinement_verification(
    candidates_arg: str,
    discrimination_arg: str,
    relationships_arg: str | None,
    out_dir_arg: str,
    *,
    case_id: str = "",
) -> dict[str, object]:
    """Emit per-candidate refinement verdicts from frozen family-to-copy artifacts.

    Every candidate in the input union receives exactly one verdict row:
    ``support_falsified``, ``support_confirmed``, or ``not_evaluated``.
    Candidates absent from the discrimination table stay ``not_evaluated``
    rather than being guessed about.
    """
    source = runtime_path(candidates_arg)
    if not source.exists():
        raise ValueError(f"refinement candidate input does not exist: {source}")
    candidates = pd.read_csv(source, sep="\t", dtype=str).fillna("")
    id_column = _candidate_id_column(candidates)
    if candidates[id_column].astype(str).str.strip().duplicated().any():
        raise ValueError("refinement candidate input contains duplicate ids")

    discrimination_path = runtime_path(discrimination_arg)
    if not discrimination_path.exists():
        raise ValueError(f"family copy discrimination table does not exist: {discrimination_path}")
    discrimination = pd.read_csv(discrimination_path, sep="\t", dtype=str).fillna("")
    for required in ("candidate_gene_id", "subfamily_vs_target"):
        if required not in discrimination.columns:
            raise ValueError(f"family copy discrimination table requires {required}")
    if discrimination["candidate_gene_id"].astype(str).str.strip().duplicated().any():
        raise ValueError("family copy discrimination table contains duplicate candidate_gene_id values")
    by_gene = {
        str(record["candidate_gene_id"]).strip(): record
        for record in discrimination.to_dict("records")
    }

    curated_sources: dict[str, bool] = {}
    if relationships_arg:
        relationships_path = runtime_path(relationships_arg)
        if not relationships_path.exists():
            raise ValueError(f"family copy relationships table does not exist: {relationships_path}")
        relationships = pd.read_csv(relationships_path, sep="\t", dtype=str).fillna("")
        for required in ("candidate_gene_id", "subfamily_source"):
            if required not in relationships.columns:
                raise ValueError(f"family copy relationships table requires {required}")
        curated_sources = {
            str(record["candidate_gene_id"]).strip():
                str(record.get("subfamily_source", "")).strip() == CURATED_SOURCE
            for record in relationships.to_dict("records")
        }

    rows: list[dict[str, str]] = []
    verdict_counts = {"support_falsified": 0, "support_confirmed": 0, "not_evaluated": 0}
    falsified: list[str] = []
    unmatched_candidates = 0

    for record in candidates.to_dict("records"):
        gene_id = str(record[id_column]).strip()
        row = by_gene.get(gene_id)
        if row is None:
            unmatched_candidates += 1
            verdict, evidence = "not_evaluated", ""
        else:
            relation = str(row.get("subfamily_vs_target", "")).strip()
            target_status = str(row.get("target_verification_status", "")).strip()
            subfamily = str(row.get("subfamily", "")).strip()
            family_id = str(row.get("family_id", "")).strip()
            target_ref = str(row.get("target_reference_gene", "")).strip() \
                or str(row.get("target_id", "")).strip()
            provenance = (
                f"family_copy_discrimination.tsv:{family_id or '-'}"
                f":subfamily={subfamily or '-'}:target={target_ref or '-'}"
            )
            if relation == "target_match":
                verdict = "support_confirmed"
                evidence = provenance
            elif (
                relation == "target_mismatch"
                and _route_a_only(record)
                and curated_sources.get(gene_id, False)
                and target_status in VERIFIED_TARGET_STATUSES
            ):
                verdict = "support_falsified"
                evidence = provenance + ":curated_cross_subfamily_vs_verified_target"
                falsified.append(gene_id)
            else:
                verdict = "not_evaluated"
                evidence = ""
        verdict_counts[verdict] += 1
        rows.append({
            "candidate_gene_id": gene_id,
            "verdict": verdict,
            "verification_source": evidence,
        })

    frame = pd.DataFrame(rows, columns=VERIFICATION_TABLE_COLUMNS)
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    table_path = out / "subfamily_refinement_verification.tsv"
    frame.to_csv(table_path, sep="\t", index=False, lineterminator="\n")

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "contract": CONTRACT,
        "case_id": case_id,
        "stage": "subfamily_refinement_verification",
        "input_candidate_count": int(len(candidates)),
        "verdict_counts": verdict_counts,
        "removed_candidate_count": len(falsified),
        "falsified_candidate_genes": sorted(set(falsified)),
        "candidates_without_discrimination_row": unmatched_candidates,
        "candidate_set_sha256": _sha256_path(source),
        "discrimination_sha256": _sha256_path(discrimination_path),
        "relationships_sha256": _sha256_path(runtime_path(relationships_arg)) if relationships_arg else None,
        "output_sha256": _sha256_path(table_path),
        "rank_contribution": "none",
        "policy": [
            "only a curated_map subfamily assignment that mismatches every expected subfamily "
            "of a verified/source-backed/curated target falsifies support",
            "an orthogroup-proxy subfamily mismatch stays not_evaluated: orthogroup granularity "
            "cannot separate subfamilies",
            "a pending_review target never falsifies anything",
            "candidates with any Route B support are out of scope: their support is independent "
            "of the Route A target projection",
            "a missing discrimination row stays not_evaluated rather than being guessed about",
            "this module emits verdicts only; apply_verified_refinement performs removal",
        ],
    }
    manifest_path = out / "subfamily_refinement_verification_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        "verification_table": str(table_path),
        "manifest": str(manifest_path),
        "input_candidate_count": int(len(candidates)),
        "verdict_counts": verdict_counts,
        "removed_candidate_count": len(falsified),
    }
