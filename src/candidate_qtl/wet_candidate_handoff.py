"""One canonical path from frozen Route A/B artifacts to a Wet-reviewable set (issue #426).

The pieces already existed -- ``supported_member_union``, ``candidate_budget``,
``narrowing`` -- but nothing joined them, so reaching a final Wet-reviewable set
from verified A/B output meant chaining commands by hand and re-deriving the stop
rule each time. This module is that join, and nothing more: it calls the existing
components in order and records what happened.

    A freeze + B freeze
      -> supported member union            (machine, non-scoring)
      -> candidate count vs Wet budget
         <= budget : freeze wet_candidates.tsv and stop membership changes
         >  target : retain all supported candidates and report the excess
         strict resource budget : bounded refinement, then re-gate
      -> reported upper pool when present; first overflow >25 compares that tier only
      -> later tiers stay in the full ledger without another comparison
      -> optional ordering inside the frozen Wet membership only

Two properties are load-bearing:

* ``wet_candidates.tsv`` is the only canonical production membership artifact.
  The historical ``stage_a/`` remains a baseline/review-budget artifact and
  ``stage_a_narrowed`` is ordering-only: neither is allowed to add or remove a
  Wet candidate.
* completion is judged on the final stage, provider state and coverage, not on
  an exit code. ``handoff_state`` carries that verdict.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import pandas as pd

from .adaptive_synthesis import (
    apply_verified_refinement,
    build_adaptive_synthesis_request,
    verify_adaptive_refinement,
)
from .candidate_budget import candidate_budget_gate
from .boundary_tier_review import EVIDENCE_SCOPED_POLICIES, UPPER_EVIDENCE_SCOPED_POLICY, REPORTED_UPPER_TIERS
from .core import runtime_path
from .narrowing import integrated_narrowing
from .subfamily_refinement import build_subfamily_refinement_verification
from .supported_member_union import build_supported_member_union
from .supported_member_union import MANIFEST_CONTRACT as SUPPORTED_UNION_CONTRACT

SCHEMA_VERSION = "1.2"
SECOND_STAGE_THRESHOLD = 25

# What the caller must do next. These are workflow states, not quality grades:
# `review_required` is a legitimate, contract-satisfying outcome.
HANDOFF_STATES = (
    "wet_reviewable",
    "wet_reviewable_degraded_coverage",
    "provisional_adopted_wet",
    "adaptive_refinement_requested",
    "review_required_over_budget",
    "review_required_no_supported_candidates",
)


def _sha256_path(path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _display_input(value):
    if isinstance(value, list):
        return [str(runtime_path(item)) for item in value]
    return str(runtime_path(value)) if value else None


def _candidate_id_column(frame: pd.DataFrame) -> str:
    for column in ("candidate_gene_id", "gene_id"):
        if column in frame.columns:
            return column
    raise ValueError("candidate set requires candidate_gene_id or gene_id")


def _has_direct_c_review_route(value: object) -> bool:
    """Read route tokens joined from already validated review records."""
    return "c" in {token.strip().casefold() for token in str(value).split("|")}


# What a Wet reviewer is being handed, per candidate. These are review handles,
# not evidence grades and not a ranking: tier 2 and tier 3 say what is known
# about a gene, never that the gene is supported.
TIER_VERIFIED = "tier1_verified_support"
TIER_REFERENCE_PHENOTYPE = "tier2_reference_phenotype"
TIER_FAMILY_OR_DOMAIN = "tier3_family_or_domain"
TIER_DARK = "tier4_dark"
FILL_TIERS = (TIER_REFERENCE_PHENOTYPE, TIER_FAMILY_OR_DOMAIN)
WET_FILL_MODES = ("tiered_interval_fill", "supported_only")
# Operational default for the A/C-first workflow: fewer final Wet genes than
# this asks the parent to rethink Route A from a fresh generation. It is a
# workload signal, not a scientific cutoff, and never pads or trims membership.
DEFAULT_MIN_ROUTE_CANDIDATES = 3


def route_a_rethink_signal(final_wet_count: int | None, *, min_route_candidates: int) -> dict[str, object]:
    """Say whether the parent should rethink Route A from scratch.

    The A/C-first workflow does not fill a short Wet set from interval
    annotation. A missing or empty Wet set requires a fresh Route-A
    generation; a set below ``min_route_candidates`` recommends one. The
    signal never changes the frozen membership.
    """
    if min_route_candidates < 1:
        raise ValueError("min_route_candidates must be at least 1")
    count = int(final_wet_count or 0)
    if count == 0:
        state, reason = "required", "no Wet membership was frozen from Route A/C candidates"
    elif count < min_route_candidates:
        state, reason = "recommended", f"final Wet set has {count} genes, below {min_route_candidates}"
    else:
        state, reason = "not_needed", f"final Wet set has {count} genes"
    return {
        "state": state,
        "final_wet_candidate_count": count,
        "min_route_candidates": min_route_candidates,
        "action": "rerun_route_a_fresh_generation" if state != "not_needed" else None,
        "reason": reason,
        "membership_effect": "none",
    }
_ASSIGNMENT_COLUMNS = (
    "family_id", "evolutionary_family_id", "orthogroup", "subfamily_id",
    "domain_ids", "domain_architecture", "orthology_reference_gene",
)


def _read_optional(path_arg: str | None) -> pd.DataFrame:
    if not path_arg:
        return pd.DataFrame()
    path = runtime_path(path_arg)
    if not path.exists():
        raise ValueError(f"explicitly supplied optional input does not exist: {path}")
    # Literal native/reference identities such as NA are valid IDs. Preserve
    # those fields while keeping the existing missing-value semantics for
    # phenotype/function/assignment metadata, booleans and numeric fields.
    columns = pd.read_csv(path, sep="\t", nrows=0).columns
    identity_columns = {
        "gene_id", "candidate_gene_id", "focal_gene_id", "reference_gene_id",
        "source_reference_gene_id", "reported_gene_id", "reported_gene_symbol",
        "reported_locus", "reference_sequence_id", "record_id", "claim_id", "hypothesis_id",
    }
    return pd.read_csv(path, sep="\t",
        dtype={column: str for column in columns if column not in identity_columns},
        converters={column: str for column in columns if column in identity_columns},
    ).fillna("")


def _read_one_review_candidates(path_arg: str, manifest_arg: str | None,
                            *, case_id: str, atlas_frame: pd.DataFrame,
                            supported_ids: set[str]) -> tuple[pd.DataFrame, dict[str, object] | None]:
    """Load a sealed A/B review lane and prove every ID is in the interval."""
    if not path_arg:
        return pd.DataFrame(), None
    path = runtime_path(path_arg)
    if not path.exists():
        raise ValueError(f"review candidates artifact does not exist: {path}")
    manifest_path = runtime_path(manifest_arg) if manifest_arg else path.with_name(f"{path.stem}_manifest.json")
    if not manifest_path.exists():
        raise ValueError("review candidates require a producer manifest")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("review candidates producer manifest is invalid JSON") from error
    if manifest.get("contract") != "route-review-candidates-1.0.0":
        raise ValueError("review candidates producer manifest has an unsupported contract")
    if manifest.get("stage") not in {"route_a_review_candidate_freeze", "route_b_review_candidate_freeze", "review_candidate_freeze"}:
        raise ValueError("review candidates producer manifest has an unsupported stage")
    if case_id and str(manifest.get("case_id", "")).strip() != case_id:
        raise ValueError("review candidates manifest case_id does not match Wet handoff")
    declared = Path(str(manifest.get("output", "")))
    if not declared.is_absolute():
        declared = manifest_path.parent / declared
    expected = str(manifest.get("artifact_sha256", "")).removeprefix("sha256:")
    if declared.resolve() != path.resolve() or len(expected) != 64 or _sha256_path(path) != expected:
        raise ValueError("review candidates artifact does not match its producer manifest")
    inputs = manifest.get("inputs")
    if not isinstance(inputs, dict):
        raise ValueError("review candidates manifest lacks input provenance")
    for name, entry in inputs.items():
        if entry is None:
            continue
        if not isinstance(entry, dict):
            raise ValueError(f"review candidates input {name} is invalid")
        input_path = runtime_path(str(entry.get("path", "")))
        input_hash = str(entry.get("sha256", "")).removeprefix("sha256:")
        if not input_path.exists() or len(input_hash) != 64 or _sha256_path(input_path) != input_hash:
            raise ValueError(f"review candidates input {name} hash mismatch or missing")
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    required = {"candidate_gene_id", "route", "support_status", "rank_contribution"}
    if not required <= set(frame.columns):
        raise ValueError("review candidates artifact lacks the common review schema")
    ids = frame["candidate_gene_id"].astype(str).str.strip()
    if ids.eq("").any():
        raise ValueError("review candidates require non-empty candidate_gene_id")
    if not ids.isin(set(atlas_frame["gene_id"].astype(str).str.strip()) if "gene_id" in atlas_frame.columns else set()).all():
        raise ValueError("review candidates contain genes outside the supplied interval atlas")
    if not frame["route"].astype(str).str.strip().str.casefold().isin({"a", "b"}).all():
        raise ValueError("review candidates route must be a or b")
    if not frame["support_status"].astype(str).str.strip().str.contains("unverified|review", case=False, regex=True).all():
        raise ValueError("review candidates must remain explicitly review-only")
    if not frame["rank_contribution"].astype(str).str.strip().eq("none").all():
        raise ValueError("review candidates must have rank_contribution=none")
    frame["candidate_gene_id"] = ids
    # Older B review artifacts omit the biological bridge columns. Their
    # already hash-checked ancestor still contains the original bound claims.
    # Rehydrate by BOTH claim ID and member ID, never by a family/name guess.
    claims_entry = inputs.get("claims")
    if manifest.get("stage") == "route_b_review_candidate_freeze" and claims_entry and "claim_id" in frame:
        original = pd.read_csv(runtime_path(str(claims_entry["path"])), sep="\t", dtype=str).fillna("")
        if {"claim_id", "member_interval_genes"} <= set(original.columns):
            bound = {}
            for record in original.to_dict("records"):
                for member in str(record["member_interval_genes"]).replace(",", ";").split(";"):
                    bound.setdefault((str(record["claim_id"]), member.strip()), []).append(record)
            if "relevance_evidence_json" not in frame:
                frame["relevance_evidence_json"] = ""
            for index, row in frame.iterrows():
                evidence = bound.get((str(row["claim_id"]), str(row["candidate_gene_id"])))
                if not evidence:
                    raise ValueError("Route B review claim is not bound to its original member")
                saved = json.loads(row["relevance_evidence_json"]) if row["relevance_evidence_json"] else []
                if not isinstance(saved, list) or any(not isinstance(item, dict) for item in saved):
                    raise ValueError("relevance evidence must be an object list")
                frame.at[index, "relevance_evidence_json"] = json.dumps(
                    saved + evidence, ensure_ascii=False, sort_keys=True)
    # Reconstruct first-stage homology from hash-checked producer evidence,
    # including historical freezes. Original paper answers stay unchanged.
    if "relevance_evidence_json" in frame:
        from .family_candidate_bridge import route_a_homolog_projection
        from .route_b_semantic_candidates import route_b_homolog_projection
        bundle_entry = inputs.get("bundle")
        units = {}
        review_taxonomy = {}
        if manifest["stage"] == "route_b_review_candidate_freeze" and bundle_entry:
            bundle = json.loads(runtime_path(bundle_entry["path"]).read_text(encoding="utf-8"))
            units = {unit["unit_id"]: unit for unit in bundle.get("units", [])}
            if manifest.get("review_membership_policy") == "initial_bibliography_v1" and inputs.get("claims_manifest"):
                from .reference_taxonomy import load_bound_taxonomy_artifacts
                transfer = json.loads(runtime_path(inputs["claims_manifest"]["path"]).read_text(encoding="utf-8"))
                if transfer.get("b_review_membership_policy") != "initial_bibliography_v1":
                    raise ValueError("Route B taxonomy projection policy mismatch")
                review_taxonomy = load_bound_taxonomy_artifacts(
                    transfer.get("reference_species_display_audit", {}).get("taxonomy_cache_artifacts", []))
        records = frame.to_dict("records")
        for index, row in frame.iterrows():
            saved = json.loads(row["relevance_evidence_json"] or "[]")
            for item in saved:
                projection = ""
                if manifest["stage"] == "route_a_review_candidate_freeze":
                    members = [r["candidate_gene_id"] for r in records if
                               (r.get("hypothesis_id"), r.get("reference_accession")) ==
                               (row.get("hypothesis_id"), row.get("reference_accession"))]
                    projection = route_a_homolog_projection(row.to_dict(), members)
                elif manifest["stage"] == "route_b_review_candidate_freeze":
                    unit = units.get(item.get("semantic_unit_id"), {})
                    members = {r["candidate_gene_id"] for r in records if r.get("claim_id") == row.get("claim_id")}
                    projection = route_b_homolog_projection(item, unit, row["candidate_gene_id"], members,
                                                           taxonomy_records=review_taxonomy)
                item["homolog_projection_json"] = projection
            frame.at[index, "relevance_evidence_json"] = json.dumps(saved, ensure_ascii=False, sort_keys=True)
    manifest = dict(manifest, output=str(path.resolve()))
    return frame, manifest


def _read_review_candidates(path_arg: str | list[str] | None,
                            manifest_arg: str | list[str] | None, *, case_id: str,
                            atlas_frame: pd.DataFrame,
                            supported_ids: set[str], route_c_members_arg: str | None = None,
                            route_c_manifest_arg: str | None = None) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    paths = [path_arg] if isinstance(path_arg, str) else list(path_arg or [])
    manifests = [manifest_arg] if isinstance(manifest_arg, str) else list(manifest_arg or [])
    if manifests and len(manifests) != len(paths):
        raise ValueError("review candidates and manifests must have matching lengths")
    frames = []
    producer_manifests = []
    if route_c_members_arg:
        from .route_c_review_candidates import read_direct_c_review_candidates
        c_frame, c_manifest = read_direct_c_review_candidates(
            route_c_members_arg, route_c_manifest_arg, case_id=case_id,
            atlas_frame=atlas_frame, supported_ids=set())
        c_manifest["eligible_record_ids"] = sorted(set(c_frame.get("record_id", [])))
        if len(c_frame):
            frames.append(c_frame)
        producer_manifests.append(c_manifest)
    for index, path in enumerate(paths):
        frame, producer = _read_one_review_candidates(
            str(path), manifests[index] if manifests else None,
            case_id=case_id, atlas_frame=atlas_frame, supported_ids=supported_ids,
        )
        frames.append(frame)
        if producer is not None:
            producer_manifests.append(producer)
    if not frames:
        return pd.DataFrame(), producer_manifests
    # Each route may legitimately use a different optional provenance vocabulary
    # (A review_status/reason versus B verification_status/reason). Concatenating
    # disjoint schemas introduces NaN; normalize before string aggregation.
    merged = pd.concat(frames, ignore_index=True).fillna("")
    rows = []
    from .relevance_tiers import annotate_relevance
    for candidate_id, group in merged.groupby("candidate_gene_id", sort=True):
        row = group.iloc[0].to_dict()
        for column in group.columns:
            if column in {"candidate_gene_id", "rank_contribution"}:
                continue
            values = []
            for value in group[column].astype(str):
                value = value.strip()
                if value and value not in values:
                    values.append(value)
            row[column] = " | ".join(values)
        row["candidate_gene_id"] = candidate_id
        row["rank_contribution"] = "none"
        # Display joins are not biological claims. Preserve and classify each
        # original record before grouping, never stitch fields across sources.
        originals = group.to_dict("records")
        annotated = annotate_relevance(
            [{"candidate_gene_id": candidate_id}], originals)
        row.update(annotated[0])
        rows.append(row)
    result = pd.DataFrame(rows)
    if len(result) and supported_ids:
        result = result[~result["candidate_gene_id"].isin(supported_ids)].copy()
    return result, producer_manifests


def _assert_optional_case(frame: pd.DataFrame, *, label: str, case_id: str) -> None:
    """Reject an explicitly case-labelled fill artifact from another case."""
    if not case_id or frame.empty or "case_id" not in frame.columns:
        return
    observed = {str(value).strip() for value in frame["case_id"] if str(value).strip()}
    if observed and observed != {case_id}:
        raise ValueError(f"explicit {label} case_id does not match the Wet handoff case")


def _reference_phenotype_genes(reference_function_frame: pd.DataFrame,
                               atlas_frame: pd.DataFrame) -> set[str]:
    """Interval genes whose reference gene has a reported mutant phenotype.

    This is the cross-species lane: the focal species has no report, but the
    reference gene the interval protein resolves to does. It is a reason to
    look at the gene in the lab, not a claim that the phenotype transfers --
    homology is not identity of function, so the tier label travels with the
    candidate and nothing about the candidate's support changes.
    """
    if reference_function_frame.empty or atlas_frame.empty:
        return set()
    if "reference_gene_id" not in reference_function_frame.columns:
        return set()
    phenotype_column = next(
        (column for column in ("mutant_phenotype", "reference_mutant_phenotype")
         if column in reference_function_frame.columns), None)
    if phenotype_column is None:
        return set()
    reported = reference_function_frame[
        reference_function_frame[phenotype_column].astype(str).str.strip().ne("")]
    reference_ids = {str(value).strip() for value in reported["reference_gene_id"] if str(value).strip()}
    if not reference_ids:
        return set()
    genes: set[str] = set()
    # The corpus may already carry the focal mapping; otherwise join through the
    # atlas' own reference gene.
    if "gene_id" in reported.columns:
        genes |= {str(value).strip() for value in reported["gene_id"] if str(value).strip()}
    if "orthology_reference_gene" in atlas_frame.columns:
        for record in atlas_frame.to_dict("records"):
            reference = str(record.get("orthology_reference_gene", "")).strip()
            if reference and reference in reference_ids:
                genes.add(str(record.get("gene_id", "")).strip())
    return {gene for gene in genes if gene}


def _direct_review_genes(route_c_members_frame: pd.DataFrame) -> set[str]:
    """Return interval genes carried by Direct-C as review-only evidence."""
    if route_c_members_frame.empty:
        return set()
    gene_column = next(
        (column for column in ("candidate_gene_id", "focal_gene_id", "gene_id")
         if column in route_c_members_frame.columns),
        None,
    )
    if gene_column is None or "support_status" not in route_c_members_frame.columns:
        return set()
    review = route_c_members_frame[
        route_c_members_frame["support_status"].astype(str).eq("reference_phenotype_review")
    ]
    return {str(value).strip() for value in review[gene_column] if str(value).strip()}


def _tier_interval_candidates(atlas_frame: pd.DataFrame,
                              reference_function_frame: pd.DataFrame,
                              route_c_members_frame: pd.DataFrame | None = None) -> dict[str, str]:
    """gene_id -> review tier, for every gene in the frozen interval inventory."""
    if atlas_frame.empty or "gene_id" not in atlas_frame.columns:
        return {}
    phenotype_genes = _reference_phenotype_genes(reference_function_frame, atlas_frame)
    if route_c_members_frame is not None:
        phenotype_genes |= _direct_review_genes(route_c_members_frame)
    tiers: dict[str, str] = {}
    for record in atlas_frame.to_dict("records"):
        gene_id = str(record.get("gene_id", "")).strip()
        if not gene_id:
            continue
        if gene_id in phenotype_genes:
            tiers[gene_id] = TIER_REFERENCE_PHENOTYPE
            continue
        assigned = any(str(record.get(column, "")).strip() for column in _ASSIGNMENT_COLUMNS)
        tiers[gene_id] = TIER_FAMILY_OR_DOMAIN if assigned else TIER_DARK
    return tiers


def _membership_basis(route_support: object) -> str:
    routes = {token.strip().casefold() for token in str(route_support or "").replace("+", ";").split(";") if token.strip()}
    if not routes or not routes <= {"a", "b", "c"}:
        raise ValueError(f"unsupported canonical route_support value: {route_support!r}")
    suffix = "_".join(route for route in ("a", "b", "c") if route in routes)
    return f"verified_route_{suffix}_support"


def _freeze_canonical_wet_membership(
    candidate_set_arg: str,
    out_dir,
    *,
    wet_budget: int,
    case_id: str,
    atlas_arg: str | None = None,
    reference_function_arg: str | None = None,
    route_c_members_arg: str | None = None,
    route_c_manifest_arg: str | None = None,
    review_candidates_arg: str | list[str] | None = None,
    review_candidates_manifest_arg: str | list[str] | None = None,
    fill_mode: str = "tiered_interval_fill",
    relevance_selection: bool = False,
    family_copy_arg: str | None = None,
    provisional_screen_ledger_arg: str | None = None,
    provisional_screen_ledger_sha256: str | None = None,
    identity_assessment_arg: str | None = None,
    identity_target_proteins_arg: str | None = None,
) -> tuple[str, str]:
    """Freeze canonical membership from bound evidence and interval context.

    Ordinary target handoff prefers eligible reported upper evidence when
    present and compares only the first overflowing tier. Later tiers remain
    in the ledger. Identity-only inventory does not fill a short set.
    Explicit strict/supported-only paths retain their legacy fill semantics.
    Support, mapping and copy attributes never become a relevance score.
    """
    source = runtime_path(candidate_set_arg)
    frame = pd.read_csv(source, sep="\t", dtype=str).fillna("")
    id_column = _candidate_id_column(frame)
    ids = frame[id_column].astype(str).str.strip()
    if ids.eq("").any() or ids.duplicated().any():
        raise ValueError("canonical Wet membership requires non-empty unique candidate ids")

    if fill_mode not in WET_FILL_MODES:
        raise ValueError(f"unknown Wet fill mode: {fill_mode}")
    frame["support_tier"] = TIER_VERIFIED
    if "route_support" not in frame.columns:
        raise ValueError("canonical Wet membership requires route_support provenance")
    # Keep review provenance visible in the final Wet artifact.  The supported
    # union has a deliberately narrower schema, so these fields must be added
    # before review rows are projected onto the membership columns.
    review_output_columns = (
        "review_routes", "review_claim_ids", "review_hypothesis_ids",
        "review_verification_statuses", "review_reasons", "review_source_ids",
    )
    for column in review_output_columns:
        if column not in frame.columns:
            frame[column] = ""
    frame["membership_basis"] = frame["route_support"].map(_membership_basis)

    atlas_frame = _read_optional(atlas_arg)
    reference_function_frame = _read_optional(reference_function_arg)
    route_c_members_frame = _read_optional(route_c_members_arg)
    if review_candidates_arg and not atlas_arg:
        raise ValueError("review candidates require an explicit interval atlas")
    _assert_optional_case(atlas_frame, label="interval atlas", case_id=case_id)
    _assert_optional_case(reference_function_frame, label="reference-function input", case_id=case_id)
    _assert_optional_case(route_c_members_frame, label="Route C members input", case_id=case_id)
    if atlas_arg and "gene_id" not in atlas_frame.columns:
        raise ValueError("explicit interval atlas requires gene_id")
    if atlas_arg:
        atlas_ids = atlas_frame["gene_id"].astype(str).str.strip()
        if atlas_ids.eq("").any() or atlas_ids.duplicated().any():
            raise ValueError("interval atlas requires non-empty unique gene IDs")
        atlas_frame["gene_id"] = atlas_ids
        if relevance_selection:
            outside_ids = set(ids) - set(atlas_ids)
            if outside_ids:
                raise ValueError(
                    "supported candidates are outside the supplied interval atlas: "
                    + ", ".join(sorted(outside_ids))
                )
    if reference_function_arg:
        phenotype_columns = {"mutant_phenotype", "reference_mutant_phenotype"}
        if "reference_gene_id" not in reference_function_frame.columns or not (
            phenotype_columns & set(reference_function_frame.columns)
        ):
            raise ValueError(
                "explicit reference-function input requires reference_gene_id and a phenotype column"
            )
    if route_c_members_arg:
        gene_columns = {"candidate_gene_id", "focal_gene_id", "gene_id"}
        if "support_status" not in route_c_members_frame.columns or not (
            gene_columns & set(route_c_members_frame.columns)
        ):
            raise ValueError("explicit Route C members input requires a gene id and support_status")
    tiers = _tier_interval_candidates(atlas_frame, reference_function_frame, route_c_members_frame)
    supported_ids = set(ids)
    review_frame, review_manifest = _read_review_candidates(
        review_candidates_arg, review_candidates_manifest_arg,
        case_id=case_id, atlas_frame=atlas_frame, supported_ids=set(),
        route_c_members_arg=route_c_members_arg if atlas_arg else None,
        route_c_manifest_arg=route_c_manifest_arg,
    )
    all_review_frame = review_frame
    if len(review_frame):
        review_frame = review_frame[~review_frame["candidate_gene_id"].isin(supported_ids)].copy()
    review_ids = set(review_frame["candidate_gene_id"]) if len(review_frame) else set()
    direct_c_review_ids = {
        str(row["candidate_gene_id"]).strip()
        for row in all_review_frame.to_dict("records")
        if _has_direct_c_review_route(row.get("route", ""))
    }
    direct_c_supported_ids = {
        str(row[id_column]).strip()
        for row in frame.to_dict("records")
        if "c" in {token.strip().casefold() for token in str(row.get("route_support", "")).replace("+", ";").split(";") if token.strip()}
    }
    protected_c_ids = direct_c_review_ids | direct_c_supported_ids
    fill_rows: list[dict[str, str]] = []
    tier_counts = {tier: 0 for tier in (TIER_VERIFIED, *FILL_TIERS)}
    tier_counts[TIER_VERIFIED] = int(len(frame))
    fill_state = "not_requested" if fill_mode == "supported_only" else "no_interval_inventory"
    if len(review_frame) and fill_mode == "tiered_interval_fill":
        for record in review_frame.to_dict("records"):
            row = {column: "" for column in frame.columns}
            for column, value in record.items():
                if column in row:
                    row[column] = str(value)
            row[id_column] = str(record["candidate_gene_id"])
            row["support_tier"] = "unverified_route_review"
            row["membership_basis"] = "route_review_candidate"
            row["route_support"] = ""
            row["review_routes"] = str(record.get("route", ""))
            row["review_claim_ids"] = str(record.get("claim_id", ""))
            row["review_hypothesis_ids"] = str(record.get("hypothesis_id", ""))
            statuses = [str(record.get(column, "")).strip() for column in ("verification_status", "review_status") if str(record.get(column, "")).strip()]
            reasons = [str(record.get(column, "")).strip() for column in ("verification_reason", "reason", "review_reason") if str(record.get(column, "")).strip()]
            row["review_verification_statuses"] = " | ".join(dict.fromkeys(statuses))
            row["review_reasons"] = " | ".join(dict.fromkeys(reasons))
            row["review_source_ids"] = str(record.get("source_id", ""))
            for column in ("route_a_supported", "route_b_supported", "route_c_supported"):
                if column in row:
                    row[column] = "False"
            fill_rows.append(row)
            tier_counts.setdefault("unverified_route_review", 0)
            tier_counts["unverified_route_review"] += 1
    if fill_mode == "tiered_interval_fill" and tiers:
        fill_state = "applied"
        available = sorted(
            (gene, tier) for gene, tier in tiers.items()
            if gene not in supported_ids and gene not in review_ids
            and (relevance_selection or tier in FILL_TIERS)
        )
        available.sort(key=lambda item: ((*FILL_TIERS, TIER_DARK).index(item[1]), item[0]))
        for gene, tier in available:
            if not relevance_selection and len(frame) + len(fill_rows) >= wet_budget:
                break
            row = {column: "" for column in frame.columns}
            row[id_column] = gene
            row["support_tier"] = tier
            row["membership_basis"] = "interval_review_coverage"
            for column in ("route_a_supported", "route_b_supported", "route_c_supported"):
                if column in row:
                    row[column] = "False"
            if "route_support" in row:
                row["route_support"] = ""
            fill_rows.append(row)
            tier_counts[tier] = tier_counts.get(tier, 0) + 1
    if frame.empty and not fill_rows:
        raise ValueError(
            "canonical Wet membership cannot be frozen from an empty supported set "
            "without reviewable interval candidates"
        )
    if fill_rows:
        fill_frame = pd.DataFrame(fill_rows, columns=frame.columns)
        membership = fill_frame if frame.empty else pd.concat([frame, fill_frame], ignore_index=True)
    else:
        membership = frame

    from .relevance_tiers import annotate_relevance
    # Rehydrate sealed source claims, including review evidence on a supported
    # gene. Membership/support policy remains independent of these annotations.
    claims = []
    from .relevance_tiers import INITIAL_TIER_POLICY
    tier_policy = INITIAL_TIER_POLICY if any(
        producer.get("review_membership_policy") == "initial_bibliography_v1" for producer in review_manifest
    ) else ""
    for record in all_review_frame.to_dict("records"):
        claims.extend(json.loads(record["relevance_evidence_json"]))
    identity_audit = None
    if identity_assessment_arg:
        if not relevance_selection or not atlas_arg:
            raise ValueError("post-freeze identity assessment requires ordinary target+atlas relevance selection")
        from .post_freeze_identity import load_identity_assessment
        identity_claims, identity_audit = load_identity_assessment(
            identity_assessment_arg, case_id=case_id, atlas_arg=atlas_arg, claims=claims,
            canonical_proteins_arg=identity_target_proteins_arg,
        )
        claims.extend(identity_claims)
    atlas_by_id = {str(row["gene_id"]): row for row in atlas_frame.to_dict("records")}
    copy_frame = _read_optional(family_copy_arg)
    copy_by_id = {}
    if family_copy_arg:
        _assert_optional_case(copy_frame, label="family/copy context", case_id=case_id)
        copy_id = _candidate_id_column(copy_frame)
        for gene, group in copy_frame.groupby(copy_id, sort=True):
            records = sorted(group.to_dict("records"), key=lambda row: json.dumps(row, sort_keys=True))
            copy_by_id[str(gene)] = {
                "copy_evidence_json": json.dumps(records, ensure_ascii=False, sort_keys=True),
                "copy_state": ";".join(sorted({str(row.get("disambiguation_status", "")).strip()
                                              for row in records if str(row.get("disambiguation_status", "")).strip()})),
            }
    annotated_rows = annotate_relevance([
        {**atlas_by_id.get(str(row[id_column]), {}), **row, **copy_by_id.get(str(row[id_column]), {})}
        for row in membership.to_dict("records")
    ], claims, policy=tier_policy)
    annotated_by_id = {str(row[id_column]): row for row in annotated_rows}
    membership = pd.DataFrame([annotated_by_id[str(gene)] for gene in membership[id_column]])
    inventory_rows = []
    for gene, row in atlas_by_id.items():
        inventory_rows.append(annotated_by_id.get(gene, {**row, "candidate_gene_id": gene}))
    inventory_rows.extend(row for gene, row in annotated_by_id.items() if gene not in atlas_by_id)
    # B eligibility is observation, never a veto over independently eligible A/C.
    b_eligibility: dict[str, set[str]] = {}
    for producer in review_manifest:
        entry = producer.get("eligibility_ledger")
        if not entry:
            continue
        ledger_path = runtime_path(str(entry.get("path", "")))
        if not ledger_path.is_file() or _sha256_path(ledger_path) != entry.get("sha256"):
            raise ValueError("Route B eligibility ledger hash mismatch")
        ledger = pd.read_csv(ledger_path, sep="\t", dtype=str).fillna("")
        if not set(ledger["candidate_gene_id"]) <= set(atlas_by_id):
            raise ValueError("Route B eligibility ledger contains genes outside interval")
        for record in ledger.to_dict("records"):
            b_eligibility.setdefault(record["candidate_gene_id"], set()).add(record["b_review_eligibility"])
    for row in inventory_rows:
        row["b_review_eligibility"] = ";".join(sorted(b_eligibility.get(str(row.get("candidate_gene_id", row.get("gene_id", ""))), set())))
    inventory_path = out_dir / "wet_relevance_inventory.tsv"
    inventory = pd.DataFrame(annotate_relevance(inventory_rows, claims, policy=tier_policy)).fillna("")
    inventory.to_csv(
        inventory_path, sep="\t", index=False, lineterminator="\n")
    boundary_scope = None
    provisional = None
    screen_snapshot = None
    if relevance_selection:
        from .boundary_tier_review import (
            DEFAULT_SELECTION_POLICY, LEGACY_POLICY,
            freeze_boundary_scope, partition_relevance_tiers,
        )
        from .relevance_tiers import T1
        eligible_gene_ids = supported_ids | review_ids
        # The explicit full-screen provisional workflow retains its historical
        # whole-tier accounting. Ordinary target handoff uses evidence-scoped
        # candidates and never expands identity-only inventory to fill a short set.
        selection_policy = (LEGACY_POLICY if provisional_screen_ledger_arg
                            else DEFAULT_SELECTION_POLICY)
        # All candidates carry their complete rows before the first membership
        # decision. Serialization order never decides inclusion within a tier.
        plan = partition_relevance_tiers(inventory.to_dict("records"), target=wet_budget,
                                         upper=max(wet_budget, SECOND_STAGE_THRESHOLD),
                                         protected_gene_ids=protected_c_ids,
                                         policy_version=selection_policy,
                                         eligible_gene_ids=(eligible_gene_ids if selection_policy in EVIDENCE_SCOPED_POLICIES else ()))
        if provisional_screen_ledger_arg or provisional_screen_ledger_sha256:
            if not provisional_screen_ledger_arg or not provisional_screen_ledger_sha256:
                raise ValueError("provisional screen ledger requires path and SHA-256")
            screen_source = runtime_path(provisional_screen_ledger_arg).resolve()
            if _sha256_path(screen_source) != provisional_screen_ledger_sha256.removeprefix("sha256:"):
                raise ValueError("provisional screen ledger SHA-256 mismatch")
            screen_snapshot = out_dir / "screen_coverage_ledger.tsv"
            shutil.copyfile(screen_source, screen_snapshot)
            if _sha256_path(screen_snapshot) != provisional_screen_ledger_sha256.removeprefix("sha256:"):
                raise ValueError("provisional screen ledger changed while copying")
            from .provisional_wet import assess_provisional_wet
            screen_rows = pd.read_csv(screen_snapshot, sep="\t", dtype=str).fillna("")
            required_screen_columns = {"case_id", "unit_id", "member_interval_genes", "screen_state", "answer_sha256"}
            if not required_screen_columns.issubset(screen_rows.columns):
                raise ValueError("provisional screen ledger lacks required columns")
            provisional = assess_provisional_wet(
                inventory.to_dict("records"), plan,
                screen_rows.loc[screen_rows["case_id"].eq(case_id)].to_dict("records"),
                case_id=case_id, protected_gene_ids=protected_c_ids,
                upper=max(wet_budget, SECOND_STAGE_THRESHOLD),
            )
        if provisional and provisional["eligible"]:
            adopted_ids = set(provisional["adopted_gene_ids"])
            membership = inventory.loc[inventory[id_column].isin(adopted_ids)].copy()
            fill_state = "provisional_evidence_limited"
        else:
            membership = pd.DataFrame(plan["frozen"] + plan["boundary"], columns=inventory.columns)
            boundary_scope = freeze_boundary_scope(
                str(inventory_path), str(out_dir / "boundary_scope"),
                inventory_sha256=_sha256_path(inventory_path), case_id=case_id,
                target=wet_budget, upper=max(wet_budget, SECOND_STAGE_THRESHOLD),
                protected_gene_ids=protected_c_ids,
                policy_version=selection_policy,
                eligible_gene_ids=(eligible_gene_ids if selection_policy in EVIDENCE_SCOPED_POLICIES else ()))
            fill_state = ("evidence_scoped_candidates" if selection_policy in EVIDENCE_SCOPED_POLICIES
                          else "whole_relevance_tiers")
        tier_counts = membership["support_tier"].value_counts().to_dict()
    # This ledger reports observable handoff events, not model competence or
    # causal truth. Missing upstream observations remain explicitly unassessed.
    retained_ids = set(membership[id_column].astype(str))
    boundary_ids = {str(row[id_column]) for row in plan["boundary"]} if relevance_selection else set()
    trace_rows = []
    for row in inventory.to_dict("records"):
        gene = str(row.get(id_column, row.get("gene_id", "")))
        trace_rows.append({
            "candidate_gene_id": gene, "relevance_tier": row["relevance_tier"],
            "candidate_admission_policy": row.get("candidate_admission_policy", ""),
            "initial_tier": row.get("initial_tier", ""),
            "initial_tier_reason": row.get("initial_tier_reason", ""),
            "paper_evidence_summary_json": row.get("paper_evidence_summary_json", ""),
            "support_state": row["support_state"],
            "supported_union_member": gene in supported_ids,
            "explicit_review_member": gene in set(all_review_frame.get("candidate_gene_id", [])),
            "bound_claim_count": len(json.loads(row["relevance_evidence_json"])),
            "adopted_first_stage": gene in retained_ids,
            "selection_disposition": ("provisional_adopted" if provisional and provisional["eligible"] and gene in retained_ids else
                                      "not_evaluated_identity_only" if provisional and provisional["eligible"] else
                                      "retained_pending_boundary_review" if gene in boundary_ids else
                                      "retained" if gene in retained_ids else
                                      "not_selected_after_first_boundary" if relevance_selection
                                      and plan.get("policy_version") == UPPER_EVIDENCE_SCOPED_POLICY
                                      and gene in set(plan.get("eligible_gene_ids", []))
                                      and gene not in set(plan.get("active_gene_ids", []))
                                      and (not plan.get("upper_priority_applied") or row["relevance_tier"] in REPORTED_UPPER_TIERS) else
                                      "not_selected_upper_tier_priority" if relevance_selection
                                      and plan.get("upper_priority_applied")
                                      and gene in set(plan.get("eligible_gene_ids", []))
                                      and gene not in set(plan.get("active_gene_ids", [])) else
                                      "not_selected_no_upstream_candidate" if relevance_selection
                                      and plan.get("policy_version") in EVIDENCE_SCOPED_POLICIES
                                      and gene not in set(plan.get("eligible_gene_ids", []))
                                      and row["relevance_tier"] != T1 else
                                      "not_reached_lower_tier" if relevance_selection else "legacy_not_selected"),
            "material_observation": (provisional["material_by_gene"].get(gene, "")
                                     if provisional and provisional["eligible"] else ""),
            "upstream_llm_assessment": "not_evaluated_by_handoff",
            "biological_exclusion": "not_established",
            "b_review_eligibility": row.get("b_review_eligibility", ""),
        })
    trace_path = out_dir / "wet_candidate_trace.tsv"
    pd.DataFrame(trace_rows).to_csv(trace_path, sep="\t", index=False, lineterminator="\n")
    unresolved_path = out_dir / "interval_unresolved.tsv" if provisional and provisional["eligible"] else None
    if unresolved_path:
        pd.DataFrame(trace_rows).to_csv(unresolved_path, sep="\t", index=False, lineterminator="\n")
    wet_path = out_dir / "wet_candidates.tsv"
    membership.to_csv(wet_path, sep="\t", index=False, lineterminator="\n")
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id,
        "producer_class": "canonical_wet_membership",
        "candidate_count": int(len(membership)),
        "supported_candidate_count": len(retained_ids & supported_ids),
        "review_coverage_candidate_count": len(retained_ids - supported_ids),
        "route_review_candidate_count": len(retained_ids & review_ids),
        "wet_budget": wet_budget,
        "above_target_count": max(0, int(len(membership)) - wet_budget),
        "wet_fill_mode": fill_mode,
        "wet_fill_state": fill_state,
        "support_tier_counts": tier_counts,
        "unfilled_budget_slots": max(0, wet_budget - int(len(membership))),
        "interval_inventory": str(runtime_path(atlas_arg)) if atlas_arg else None,
        "reference_function_input": str(runtime_path(reference_function_arg)) if reference_function_arg else None,
        "route_c_members_input": str(runtime_path(route_c_members_arg)) if route_c_members_arg else None,
        "review_candidates_input": _display_input(review_candidates_arg),
        "review_candidates_manifest": _display_input(review_candidates_manifest_arg),
        "post_freeze_identity_assessment": identity_audit,
        "candidate_id_column": id_column,
        "source_candidate_set": str(source),
        "source_sha256": _sha256_path(source),
        "wet_candidates_sha256": _sha256_path(wet_path),
        "relevance_selection_applied": relevance_selection,
        "first_stage_selection_policy_version": plan.get("policy_version") if relevance_selection else None,
        "eligible_first_stage_candidate_count": len(plan.get("eligible_gene_ids", [])) if relevance_selection else None,
        **({"active_first_stage_candidate_count": len(plan["active_gene_ids"]),
            "upper_priority_applied": plan["upper_priority_applied"]}
           if relevance_selection and plan.get("policy_version") == UPPER_EVIDENCE_SCOPED_POLICY else {}),
        "boundary_scope_manifest": boundary_scope["manifest"] if boundary_scope else None,
        "boundary_scope_sha256": boundary_scope["manifest_sha256"] if boundary_scope else None,
        "relevance_selection_state": ("provisional_evidence_limited" if provisional and provisional["eligible"] else
                                      boundary_scope["state"] if boundary_scope else "legacy_selection"),
        "candidate_trace": str(trace_path),
        "candidate_trace_sha256": _sha256_path(trace_path),
        "interval_unresolved": str(unresolved_path) if unresolved_path else None,
        "interval_unresolved_sha256": _sha256_path(unresolved_path) if unresolved_path else None,
        "provisional_wet": provisional,
        "screen_coverage_ledger": str(screen_snapshot) if screen_snapshot else None,
        "screen_coverage_ledger_sha256": _sha256_path(screen_snapshot) if screen_snapshot else None,
        "relevance_inventory": str(inventory_path),
        "relevance_inventory_sha256": _sha256_path(inventory_path),
        "family_copy_context": str(runtime_path(family_copy_arg)) if family_copy_arg else None,
        "family_copy_context_sha256": _sha256_path(runtime_path(family_copy_arg)) if family_copy_arg else None,
        "membership_rule": ("provisional adoption of every concrete T1-T3 and protected C candidate after "
            "complete screen accounting; all identity-only rows remain not_evaluated in a hashed interval ledger"
            if provisional and provisional["eligible"] else
            "upper_evidence_scoped_v3 retains eligible reported T1-T3 and protected C when any upper evidence exists; "
            "lower eligible candidates remain in the inventory without comparison or biological exclusion; "
            "if no upper evidence exists use the original eligible scope; over 25 compare only the first overflowing tier; "
            "higher tiers remain frozen and later tiers stay in the ledger without within-tier ranking or forced reduction"
            if relevance_selection and plan.get("policy_version") == UPPER_EVIDENCE_SCOPED_POLICY else
            "evidence-scoped selection retains all T1 and existing eligible union/review candidates; if their cumulative tiers exceed 25, "
            "the first overflowing tier and every lower-tier eligible candidate enter one comparison scope; "
            "identity-only inventory does not fill unused target slots"
            if relevance_selection and plan.get("policy_version") in EVIDENCE_SCOPED_POLICIES else
            "whole relevance tiers; preserve T1, review only the first overflowing tier; "
            "T4/U remain rescue material, no numerical score or within-tier ID truncation"
            if relevance_selection else
            "every supported candidate is retained as tier1 even above the operational target; "
            "remaining review slots are filled from the frozen interval inventory in tier order "
            "(tier2 reference phenotype, then tier3 family/domain) so the bench receives a budget's "
            "worth of reviewable material; dark genes are never used as padding, a filled candidate "
            "carries no support and no score, and an over-budget supported set is never truncated"
        ),
        "tier_semantics": {
            TIER_VERIFIED: "verified Route A, Route B and/or Route C support; membership_basis records the exact route set",
            "route_c": "verified Direct-C support is retained as provenance and never contributes a rank score",
            TIER_REFERENCE_PHENOTYPE: (
                "the reference gene this interval protein resolves to has a reported mutant "
                "phenotype in another species; homology is not identity of function, so this is a "
                "reason to look, not transferred support"
            ),
            TIER_FAMILY_OR_DOMAIN: "family, orthogroup, subfamily or domain assignment only",
            TIER_DARK: ("identity unresolved; excluded unless already an eligible route candidate" if relevance_selection
                        else "no assignment; never used to fill a legacy review slot"),
        },
        "ordering_rule": (
            "membership is frozen before any integrated-narrowing ordering; "
            "ordering may only reorder this exact set"
        ),
        "rank_contribution": "none",
    }
    union_manifest = out_dir / "supported_member_union" / "supported_member_union_manifest.json"
    manifest["supported_member_union_manifest"] = str(union_manifest) if union_manifest.exists() else None
    manifest["supported_member_union_manifest_sha256"] = _sha256_path(union_manifest) if union_manifest.exists() else None
    manifest["fill_input_sha256"] = {
        "interval_inventory": _sha256_path(runtime_path(atlas_arg)) if atlas_arg else None,
        "reference_function": _sha256_path(runtime_path(reference_function_arg)) if reference_function_arg else None,
        "route_c_members": _sha256_path(runtime_path(route_c_members_arg)) if route_c_members_arg else None,
    }
    manifest_path = out_dir / "wet_candidates_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return str(wet_path), str(manifest_path)


def _wet_membership_ids(wet_candidates_arg: str) -> set[str]:
    wet = pd.read_csv(runtime_path(wet_candidates_arg), sep="\t", dtype=str).fillna("")
    id_column = _candidate_id_column(wet)
    return {value.strip() for value in wet[id_column].astype(str) if value.strip()}


def _restrict_ranked_to_wet_membership(
    ranked_arg: str,
    wet_candidates_arg: str,
    out_dir,
) -> str:
    """Create an ordering-only ranked view containing exactly the Wet membership.

    This is intentionally not a candidate selector. It preserves the existing
    deterministic row order and columns, but drops every row outside the already
    frozen ``wet_candidates.tsv`` set. Therefore eliminated / unsupported genes
    can never consume a Wet priority slot (#503).
    """
    ranked_path = runtime_path(ranked_arg)
    ranked = pd.read_csv(ranked_path, sep="\t", dtype=str).fillna("")
    if "gene_id" not in ranked.columns:
        raise ValueError("ranked candidates require gene_id for Wet ordering")
    wet_ids = _wet_membership_ids(wet_candidates_arg)
    ranked_normalized = ranked["gene_id"].astype(str).str.strip()
    ranked_ids = set(ranked_normalized)
    missing = sorted(wet_ids - ranked_ids)
    if missing:
        raise ValueError(
            "canonical Wet membership contains candidates absent from ranked input: "
            + ", ".join(missing)
        )
    restricted = ranked[ranked_normalized.isin(wet_ids)].copy()
    if len(restricted) != len(wet_ids):
        raise ValueError("Wet ordering input failed to preserve the exact canonical membership")
    ordering_input_dir = out_dir / "ordering_input"
    ordering_input_dir.mkdir(parents=True, exist_ok=True)
    restricted_path = ordering_input_dir / "wet_membership_ranked.tsv"
    restricted.to_csv(restricted_path, sep="\t", index=False, lineterminator="\n")
    return str(restricted_path)


def _restrict_narrowing_context_to_wet_membership(
    wet_candidates_arg: str,
    out_dir,
    *,
    family_recovery_arg: str | None,
    family_copy_arg: str | None,
    family_copy_shortlist_arg: str | None,
) -> dict[str, str | None]:
    """Project gene-bearing narrowing inputs onto the frozen Wet membership.

    ``integrated_narrowing`` intentionally fails closed when an optional input
    names a gene absent from its ranked universe.  Once the ranked universe is
    restricted to ``wet_candidates.tsv``, the gene-bearing context tables must
    be projected to the same already-frozen set as well.  This projection cannot
    select membership: it only removes context rows/tokens for genes that are
    already outside canonical membership (#503).
    """
    wet_ids = _wet_membership_ids(wet_candidates_arg)
    ordering_dir = out_dir / "ordering_input"
    ordering_dir.mkdir(parents=True, exist_ok=True)
    result: dict[str, str | None] = {
        "family_recovery": None,
        "family_copy": None,
        "family_copy_shortlist": None,
    }

    if family_recovery_arg:
        source = pd.read_csv(runtime_path(family_recovery_arg), sep="\t", dtype=str).fillna("")
        if "member_interval_genes" not in source.columns:
            raise ValueError("family recovery input requires member_interval_genes")
        rows: list[dict[str, object]] = []
        for row in source.to_dict("records"):
            members = [
                token.strip() for token in str(row.get("member_interval_genes", "")).split(";")
                if token.strip() in wet_ids
            ]
            if not members:
                continue
            projected = dict(row)
            projected["member_interval_genes"] = ";".join(sorted(set(members)))
            rows.append(projected)
        projected_frame = pd.DataFrame(rows, columns=source.columns)
        path = ordering_dir / "wet_membership_family_recovery.tsv"
        projected_frame.to_csv(path, sep="\t", index=False, lineterminator="\n")
        result["family_recovery"] = str(path)

    for key, value, filename in (
        ("family_copy", family_copy_arg, "wet_membership_family_copy.tsv"),
        ("family_copy_shortlist", family_copy_shortlist_arg, "wet_membership_family_copy_shortlist.tsv"),
    ):
        if not value:
            continue
        source = pd.read_csv(runtime_path(value), sep="\t", dtype=str).fillna("")
        if "candidate_gene_id" not in source.columns:
            raise ValueError(f"{key} input requires candidate_gene_id")
        mask = source["candidate_gene_id"].astype(str).str.strip().isin(wet_ids)
        projected_frame = source[mask].copy()
        path = ordering_dir / filename
        projected_frame.to_csv(path, sep="\t", index=False, lineterminator="\n")
        result[key] = str(path)

    return result


# Supplying one route is a study-design choice, not degradation. Only a route
# that was supplied and produced nothing verified degrades the verdict.
DEGRADED_COVERAGE_STATES = frozenset({
    "route_a_declared_but_unsupported",
    "route_b_declared_but_unverified",
    "route_c_declared_but_unverified",
    "both_routes_degraded",
    "multiple_routes_degraded",
    "route_a_initial_candidate_evidence_unavailable",
    "route_a_initial_import_incomplete",
    "route_c_partial_input_coverage",
})


def _input_coverage(
    union_manifest: dict,
    *,
    route_a_supplied: bool,
    route_b_supplied: bool,
    route_c_supplied: bool = False,
    initial_b_candidate_count: int = 0,
    initial_a_coverage: dict | None = None,
    route_c_coverage: dict | None = None,
) -> tuple[str, dict[str, object]]:
    """Classify supplied evidence, preserving legacy strict-support meanings.

    A route that was never supplied is a study design choice, not degradation.
    Valid initial candidate evidence need not produce strict support. Bound
    failures remain visible; unavailable execution receipts remain unknown.
    """
    if not union_manifest:
        return "not_evaluated_reentry", {"reason": "re-entry after refinement reuses an already-built union"}

    a_supported = int(union_manifest.get("route_a_supported_rows", 0) or 0)
    b_verified = int(
        (union_manifest.get("route_b_verified_candidate_rows", 0) or 0)
        or (union_manifest.get("route_b_verified_claim_rows", 0) or 0)
    )
    b_declared = int(
        (union_manifest.get("route_b_candidate_artifact_rows", 0) or 0)
        or (union_manifest.get("route_b_claim_rows", 0) or 0)
    )
    c_verified = int(union_manifest.get("route_c_verified_support_rows", 0) or 0)
    c_declared = int(union_manifest.get("route_c_rows", 0) or 0)
    detail = {
        "route_a_supplied": route_a_supplied,
        "route_a_supported_rows": a_supported,
        "route_b_supplied": route_b_supplied,
        "route_b_declared_rows": b_declared,
        "route_b_verified_rows": b_verified,
        "route_b_initial_candidate_count": initial_b_candidate_count,
        "route_c_supplied": route_c_supplied,
        "route_c_declared_rows": c_declared,
        "route_c_verified_rows": c_verified,
    }
    if initial_a_coverage is not None:
        detail["route_a_initial_coverage"] = initial_a_coverage
    if route_c_coverage is not None:
        detail["route_c_input_observations"] = route_c_coverage
    degraded = []
    if route_a_supplied:
        if initial_a_coverage is None:
            if a_supported == 0:
                degraded.append("route_a_declared_but_unsupported")
        elif initial_a_coverage.get("import_state") == "incomplete":
            degraded.append("route_a_initial_import_incomplete")
        elif not initial_a_coverage.get("usable_candidate_count"):
            degraded.append("route_a_initial_candidate_evidence_unavailable")
    if route_b_supplied and b_verified == 0 and initial_b_candidate_count == 0:
        degraded.append("route_b_declared_but_unverified")
    if route_c_supplied:
        if route_c_coverage is None:
            if c_verified == 0:
                degraded.append("route_c_declared_but_unverified")
        elif route_c_coverage.get("unresolved_record_count", 0):
            degraded.append("route_c_partial_input_coverage")
    if not degraded:
        if (initial_a_coverage is not None and initial_a_coverage.get("execution_state") == "unknown"
                or route_c_coverage is not None and route_c_coverage.get("state", "").startswith("unknown")):
            return "execution_coverage_unknown", detail
        supplied_count = sum((route_a_supplied, route_b_supplied, route_c_supplied))
        return ("complete" if supplied_count > 1 else "single_route_by_design"), detail
    detail["degraded"] = degraded
    return ("multiple_routes_degraded" if len(degraded) > 1 else degraded[0]), detail


def _coverage_input(entry: dict, *, base: Path | None = None) -> pd.DataFrame:
    path = runtime_path(str(entry["path"]))
    if not path.is_absolute() and base is not None:
        path = base / path
    if not path.is_file() or _sha256_path(path) != str(entry.get("sha256", "")).removeprefix("sha256:"):
        raise ValueError("coverage input ancestor hash mismatch")
    return _read_optional(str(path))


def _initial_a_input_coverage(producers: list[dict]) -> dict | None:
    """Observe sealed initial-candidate work, never infer an execution receipt."""
    genes, observations = set(), []
    for producer in producers:
        if producer.get("stage") != "route_a_review_candidate_freeze":
            continue
        inputs = producer.get("inputs", {})
        if not inputs.get("route_a_hypotheses"):
            continue
        hypotheses = _coverage_input(inputs["route_a_hypotheses"])
        policy = hypotheses.get("candidate_admission_policy", pd.Series(dtype=str))
        if hypotheses.empty or len(policy) != len(hypotheses) or not policy.eq("initial_bibliography_v1").all():
            continue  # Legacy/content generations retain their original coverage meaning.
        bibliography = hypotheses.get("bibliography_state", pd.Series(dtype=str))
        observation = {"observed_hypothesis_count": len(hypotheses),
            "bibliography_state_counts": bibliography.value_counts().to_dict(),
            "import_state": "unknown", "execution_state": "unknown",
            "reason": "producer has no bound import completion/assigned-roster receipt; support is a separate attribute"}
        # Bound failures are observable, but successful rows alone do not prove
        # assigned-roster completeness. Never discover arbitrary siblings.
        if inputs.get("import_status"):
            status = _coverage_input(inputs["import_status"])
            successful_rows_only = bool(len(status) and "status" in status and status.status.isin({"imported", "proposals_imported", "no_claims"}).all())
            observation.update(import_state="unknown" if successful_rows_only else "incomplete",
                               execution_state="unknown" if successful_rows_only else "partial")
            observation["import_status_counts"] = status.get("status", pd.Series(dtype=str)).value_counts().to_dict()
        valid_ids = set(hypotheses.loc[bibliography.eq("resolved"), "hypothesis_id"])
        if inputs.get("route_a_projection_panel"):
            panel = _coverage_input(inputs["route_a_projection_panel"])
            valid = (panel.get("panel_status", pd.Series(index=panel.index, dtype=str)).eq("frozen")
                     & panel.get("sequence_sha256", pd.Series(index=panel.index, dtype=str)).str.fullmatch(r"[0-9a-f]{64}").fillna(False)
                     & panel.get("projection_eligibility", pd.Series(index=panel.index, dtype=str)).isin({"review_only", "verified_support"}))
            valid_ids &= set(panel.loc[valid, "hypothesis_id"])
            observation["projection_hypothesis_count"] = len(set(panel.get("hypothesis_id", [])))
            observation["panel_status_counts"] = panel.get("panel_status", pd.Series(dtype=str)).value_counts().to_dict()
        else:
            valid_ids = set()
            observation["projection_hypothesis_count"] = 0
        review = _read_optional(str(producer["output"]))
        if "hypothesis_id" in review:
            usable = review[review.hypothesis_id.isin(valid_ids)]
            genes.update(usable.get("candidate_gene_id", []))
        observations.append(observation)
    if not observations:
        return None
    return {"usable_candidate_count": len(genes), "observations": observations,
        "import_state": "incomplete" if any(item["import_state"] == "incomplete" for item in observations) else "unknown",
        "execution_state": "partial" if any(item["execution_state"] == "partial" for item in observations) else "unknown",
        "verified_support_required": False}


def _route_c_input_coverage(manifest_arg: str, producers: list[dict]) -> dict:
    """Describe supplied C source/projection coverage, not biological absence."""
    path = runtime_path(manifest_arg)
    manifest = json.loads(path.read_text(encoding="utf-8"))
    inputs = manifest.get("inputs", {})
    if not all(isinstance(inputs.get(key), dict) for key in ("frozen", "projection")):
        return {"state": "unknown_missing_projection_ancestry", "execution_state": "unknown"}
    frozen = _coverage_input(inputs["frozen"], base=path.parent)
    projection = _coverage_input(inputs["projection"], base=path.parent)
    accepted = set().union(*(set(item.get("eligible_record_ids", [])) for item in producers if item.get("route") == "c"))
    member_path = runtime_path(str(manifest.get("output", "")))
    if not member_path.is_absolute():
        member_path = path.parent / member_path
    member = _read_optional(str(member_path)) if manifest.get("output") else pd.DataFrame()
    if not member.empty and {"record_id", "support_status"} <= set(member):
        accepted |= set(member.loc[member.support_status.eq("verified_direct_support"), "record_id"])
    counts: dict[str, int] = {}
    for record in frozen.to_dict("records"):
        rows = projection[projection.record_id.eq(record["record_id"])]
        states = set(rows.get("projection_state", []))
        if record.get("source_verification_state") != "verified":
            state = "source_pending_or_unverified"
        elif record.get("record_id") in accepted and "inside_current_qtl" in states:
            state = "valid_inside_support_or_review"
        elif states == {"no_interval_match"}:
            state = "verified_source_no_interval_match"
        else:
            state = "source_verified_projection_unresolved"
        counts[state] = counts.get(state, 0) + 1
    unresolved = sum(value for key, value in counts.items() if key not in {"valid_inside_support_or_review", "verified_source_no_interval_match"})
    return {"state": "partial_input_coverage" if unresolved else "evaluated_supplied_evidence" if len(frozen) else "unknown_empty_ledger",
        "execution_state": "not_assessed_by_handoff", "record_count": len(frozen), "record_state_counts": counts,
        "source_state_counts": frozen.get("source_verification_state", pd.Series(dtype=str)).value_counts().to_dict(),
        "projection_state_counts": projection.get("projection_state", pd.Series(dtype=str)).value_counts().to_dict(),
        "unresolved_record_count": unresolved, "no_interval_match_is_biological_negative": False}


def _initial_b_coverage(producers: list[dict]) -> tuple[int, list[dict]]:
    """Read completed import/metadata work, separately from strict support.

    Producer inputs were validated by _read_review_candidates. Missing status
    material never becomes a completed execution merely because a gene survived.
    """
    genes, observations = set(), []
    for producer in producers:
        if producer.get("review_membership_policy") != "initial_bibliography_v1":
            continue
        transfer_spec = producer["inputs"]["claims_manifest"]
        transfer = json.loads(runtime_path(transfer_spec["path"]).read_text(encoding="utf-8"))
        ancestors = transfer["ancestors"]
        source_spec, proposal_spec = ancestors["source_verification"], ancestors["proposals"]
        for spec in (source_spec, proposal_spec):
            if _sha256_path(runtime_path(spec["path"])) != spec["sha256"]:
                raise ValueError("initial coverage ancestor hash mismatch")
        source = _read_optional(source_spec["path"])
        status_spec = ancestors.get("proposal_status")
        complete = False
        observation = {"status": "not_assessed_missing_import_status", "verified_support_required": False}
        if status_spec:
            status_path = runtime_path(status_spec["path"])
            digest = _sha256_path(status_path)
            if digest != status_spec["sha256"]:
                raise ValueError("initial coverage status hash mismatch")
            status = _read_optional(str(status_path))
            complete = bool(len(status) and "status" in status and status["status"].isin({"proposals_imported", "no_claims"}).all())
            observation.update(status="import_complete" if complete else "import_incomplete",
                               path=str(status_path), sha256=digest)
            if _sha256_path(status_path) != digest:
                raise ValueError("initial coverage status changed during handoff")
        unavailable = source.get("cache_investigation_state", pd.Series(dtype=str)).eq("provider_failure").any()
        unavailable |= source.get("bibliography_state", pd.Series(dtype=str)).isin({"provider_failure", "not_attempted"}).any()
        observation["source_unavailable"] = bool(unavailable)
        observations.append(observation)
        if complete and not unavailable:
            frame = _read_optional(str(producer["output"]))
            genes.update(frame.get("candidate_gene_id", []))
    complete = all(item["status"] == "import_complete" and not item["source_unavailable"] for item in observations)
    return len(genes) if complete else 0, observations


def wet_candidate_handoff(
    out_dir_arg: str,
    *,
    route_a_bridge_arg: str | None = None,
    route_b_claims_arg: str | None = None,
    route_b_candidates_arg: str | None = None,
    route_c_members_arg: str | None = None,
    route_a_manifest_arg: str | None = None,
    route_b_claims_manifest_arg: str | None = None,
    route_b_candidates_manifest_arg: str | None = None,
    route_c_manifest_arg: str | None = None,
    ranked_arg: str | None = None,
    family_recovery_arg: str | None = None,
    convergence_arg: str | None = None,
    family_copy_arg: str | None = None,
    family_copy_shortlist_arg: str | None = None,
    refined_candidates_arg: str | None = None,
    refinement_manifest_arg: str | None = None,
    refinement_original_candidates_arg: str | None = None,
    refinement_original_union_manifest_arg: str | None = None,
    wet_budget: int = 20,
    budget_policy: str = "target",
    adaptive_round: int = 0,
    max_adaptive_rounds: int = 1,
    ambiguity_arg: str | None = None,
    interval_atlas_arg: str | None = None,
    reference_function_arg: str | None = None,
    review_candidates_arg: str | list[str] | None = None,
    review_candidates_manifest_arg: str | list[str] | None = None,
    wet_fill_mode: str = "tiered_interval_fill",
    deterministic_subfamily_refinement: bool = False,
    provisional_screen_ledger_arg: str | None = None,
    provisional_screen_ledger_sha256: str | None = None,
    case_id: str = "",
    identity_assessment_arg: str | None = None,
    identity_target_proteins_arg: str | None = None,
    boundary_route_b_bundle_arg: str | None = None,
    min_route_candidates: int = DEFAULT_MIN_ROUTE_CANDIDATES,
) -> dict[str, object]:
    """Run the canonical A/B/C -> Wet-candidate handoff.

    ``refined_candidates_arg`` is how a completed adaptive round re-enters. It
    must be accompanied by the adaptive-refinement manifest and the original
    supported union so lineage, hashes, membership and sourced falsifications
    can be verified before the smaller set is trusted.

    ``deterministic_subfamily_refinement`` runs the one bounded refinement pass
    in-line instead of stopping at a request: when the gate says the supported
    union exceeds the Wet budget, curated cross-subfamily verdicts against the
    verified Route A target (see ``subfamily_refinement``) are applied through
    the existing machine-verified refinement contract and the gate re-runs once.
    It needs the frozen ``family_copy_discrimination.tsv`` as ``family_copy_arg``
    and never truncates: if the sourced verdicts do not bring the set within
    budget, the handoff stops as review-required exactly as before.
    """
    # Strict capacity/refinement and supported-only are explicit legacy modes.
    # The ordinary A/B handoff uses the full interval roster for relevance.
    relevance_selection = bool(interval_atlas_arg and wet_fill_mode == "tiered_interval_fill"
                               and budget_policy == "target" and not refined_candidates_arg)
    if (provisional_screen_ledger_arg or provisional_screen_ledger_sha256) and not relevance_selection:
        raise ValueError("provisional Wet requires ordinary target+atlas relevance selection")
    if relevance_selection and not case_id.strip():
        raise ValueError("ordinary target+atlas handoff requires a non-empty case_id for boundary lineage")
    if identity_assessment_arg and not relevance_selection:
        raise ValueError("post-freeze identity assessment requires ordinary target+atlas relevance selection")
    if bool(identity_assessment_arg) != bool(identity_target_proteins_arg):
        raise ValueError("identity assessment requires its canonical target proteins input")
    if relevance_selection and deterministic_subfamily_refinement:
        raise ValueError("relevance selection uses boundary review; legacy deterministic refinement requires strict capacity or supported_only")
    if min_route_candidates < 1:
        raise ValueError("min_route_candidates must be at least 1")
    out = runtime_path(out_dir_arg)
    # Every handoff is a new generation. Some downstream producers overwrite
    # their files before the immutable boundary scope rejects a reused path;
    # reject here so a retry cannot invalidate the original freeze or resume.
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError(
            "Wet handoff output directory must be new or empty; "
            "use a new output directory and preserve the previous generation"
        )
    out.mkdir(parents=True, exist_ok=True)
    stages: list[dict[str, object]] = []

    # --- 1. machine-first supported member union -------------------------
    if refined_candidates_arg:
        # Re-entry after a verified refinement round: the union was already
        # built, and rebuilding it would undo the refinement.
        union_path = runtime_path(refined_candidates_arg)
        if not union_path.exists():
            raise ValueError(f"refined candidate set does not exist: {union_path}")
        if adaptive_round < 1:
            raise ValueError("refined candidate re-entry requires adaptive_round >= 1")
        if (
            not refinement_manifest_arg
            or not refinement_original_candidates_arg
            or not refinement_original_union_manifest_arg
        ):
            raise ValueError(
                "refined candidate re-entry requires refinement_manifest_arg, "
                "refinement_original_candidates_arg and refinement_original_union_manifest_arg"
            )
        original_union_manifest_path = runtime_path(refinement_original_union_manifest_arg)
        if not original_union_manifest_path.exists():
            raise ValueError(f"original supported-union manifest does not exist: {original_union_manifest_path}")
        original_union_manifest = json.loads(original_union_manifest_path.read_text(encoding="utf-8"))
        if original_union_manifest.get("contract") != SUPPORTED_UNION_CONTRACT:
            raise ValueError("original candidate set is not bound to the canonical supported-union contract")
        if original_union_manifest.get("stage") != "supported_member_union":
            raise ValueError("original supported-union manifest has the wrong stage")
        original_inputs = original_union_manifest.get("inputs")
        expected_generation_digest = hashlib.sha256(
            json.dumps(original_inputs, sort_keys=True, ensure_ascii=False).encode()
        ).hexdigest() if isinstance(original_inputs, list) else ""
        if not expected_generation_digest or original_union_manifest.get("generation_digest") != expected_generation_digest:
            raise ValueError("original supported-union manifest generation digest mismatch")
        original_candidates_path = runtime_path(refinement_original_candidates_arg)
        if original_union_manifest.get("output_sha256") != _sha256_path(original_candidates_path):
            raise ValueError("original supported-union manifest hash mismatch")
        if case_id and str(original_union_manifest.get("case_id", "")) != case_id:
            raise ValueError("original supported-union manifest case_id does not match the Wet handoff case")
        refinement_lineage = verify_adaptive_refinement(
            str(union_path),
            refinement_manifest_arg,
            refinement_original_candidates_arg,
            case_id=case_id,
        )
        if int(refinement_lineage["input_candidate_count"]) <= wet_budget:
            raise ValueError("adaptive refinement re-entry requires an original set above the Wet budget")
        refinement_lineage["original_supported_union_manifest"] = str(original_union_manifest_path)
        refinement_lineage["original_supported_union_manifest_sha256"] = _sha256_path(
            original_union_manifest_path
        )
        union_result: dict[str, object] = {
            "candidate_set": str(union_path),
            "reentry": True,
            "refinement_lineage": refinement_lineage,
        }
        union_manifest = {}
        stages.append({
            "stage": "supported_member_union",
            "state": "reused_verified_refined_set",
            "refinement_lineage": refinement_lineage,
        })
    else:
        union_dir = out / "supported_member_union"
        union_result = build_supported_member_union(
            str(union_dir),
            route_a_bridge_arg=route_a_bridge_arg,
            route_b_claims_arg=route_b_claims_arg,
            route_b_candidates_arg=route_b_candidates_arg,
            route_c_members_arg=route_c_members_arg,
            route_a_manifest_arg=route_a_manifest_arg,
            route_b_claims_manifest_arg=route_b_claims_manifest_arg,
            route_b_candidates_manifest_arg=route_b_candidates_manifest_arg,
            route_c_manifest_arg=route_c_manifest_arg,
            case_id=case_id,
        )
        union_path = runtime_path(str(union_result.get("candidate_set")))
        union_manifest = json.loads(
            runtime_path(str(union_result["manifest"])).read_text(encoding="utf-8")
        )
        if union_manifest.get("contract") != SUPPORTED_UNION_CONTRACT or union_manifest.get("stage") != "supported_member_union":
            raise ValueError("Wet handoff requires the canonical supported-union manifest contract")
        if union_manifest.get("output_sha256") != _sha256_path(union_path):
            raise ValueError("supported-union manifest output hash mismatch")
        if case_id and str(union_manifest.get("case_id", "")).strip() != case_id:
            raise ValueError("supported-union manifest case_id does not match Wet handoff case")
        if not str(union_manifest.get("generation_digest", "")).strip():
            raise ValueError("supported-union manifest lacks generation digest")
        stages.append({"stage": "supported_member_union", "state": "built"})

    coverage_state, coverage_detail = _input_coverage(
        union_manifest,
        route_a_supplied=bool(route_a_bridge_arg),
        route_b_supplied=bool(route_b_claims_arg or route_b_candidates_arg),
        route_c_supplied=bool(route_c_members_arg),
    )

    union_frame = pd.read_csv(union_path, sep="\t", dtype=str).fillna("")
    id_column = _candidate_id_column(union_frame) if len(union_frame.columns) else "candidate_gene_id"
    candidate_count = int(len(union_frame))
    direct_c_supported_count = int(sum(
        "c" in {token.strip().casefold() for token in str(row.get("route_support", "")).replace("+", ";").split(";") if token.strip()}
        for row in union_frame.to_dict("records")
    ))

    # --- 2. budget gate ---------------------------------------------------
    if budget_policy not in {"target", "strict"}:
        raise ValueError("budget_policy must be target or strict")
    gate_dir = out / "candidate_budget"
    review_candidate_count = 0
    direct_c_review_count = 0
    versioned_review_input = False
    initial_b_candidate_count, initial_b_coverage = 0, []
    initial_a_coverage = None
    route_c_coverage = None
    if (review_candidates_arg or (route_c_members_arg and interval_atlas_arg)) and wet_fill_mode == "tiered_interval_fill":
        if not interval_atlas_arg:
            raise ValueError("review candidates require an explicit interval atlas")
        review_frame_for_gate, review_producers_for_gate = _read_review_candidates(
            review_candidates_arg, review_candidates_manifest_arg,
            case_id=case_id, atlas_frame=_read_optional(interval_atlas_arg),
            supported_ids=set(union_frame[id_column].astype(str).str.strip()),
            route_c_members_arg=route_c_members_arg, route_c_manifest_arg=route_c_manifest_arg,
        )
        versioned_review_input = any(producer.get("review_membership_policy") in {"source_function_v1", "initial_bibliography_v1"}
                                           for producer in review_producers_for_gate)
        if relevance_selection:
            initial_b_candidate_count, initial_b_coverage = _initial_b_coverage(review_producers_for_gate)
            initial_a_coverage = _initial_a_input_coverage(review_producers_for_gate)
            if route_c_manifest_arg:
                route_c_coverage = _route_c_input_coverage(route_c_manifest_arg, review_producers_for_gate)
            coverage_state, coverage_detail = _input_coverage(
                union_manifest, route_a_supplied=bool(route_a_bridge_arg),
                route_b_supplied=bool(route_b_claims_arg or route_b_candidates_arg),
                route_c_supplied=bool(route_c_members_arg), initial_b_candidate_count=initial_b_candidate_count,
                initial_a_coverage=initial_a_coverage, route_c_coverage=route_c_coverage)
            coverage_detail["initial_b_import_and_source_observations"] = initial_b_coverage
        review_candidate_count = int(len(review_frame_for_gate))
        direct_c_review_count = sum(
            _has_direct_c_review_route(route)
            for route in review_frame_for_gate.get("route", pd.Series(dtype=str))
        )
    reviewable_fill_count = 0
    if candidate_count == 0 and wet_fill_mode == "tiered_interval_fill":
        if not relevance_selection:
            reviewable_tiers = _tier_interval_candidates(
                _read_optional(interval_atlas_arg),
                _read_optional(reference_function_arg),
                _read_optional(route_c_members_arg),
            )
            reviewable_fill_count = sum(tier in FILL_TIERS for tier in reviewable_tiers.values())
        if review_candidates_arg or (route_c_members_arg and interval_atlas_arg):
            review_input = _read_review_candidates(
                review_candidates_arg, review_candidates_manifest_arg,
                case_id=case_id, atlas_frame=_read_optional(interval_atlas_arg), supported_ids=set(),
                route_c_members_arg=route_c_members_arg, route_c_manifest_arg=route_c_manifest_arg,
            )[0]
            reviewable_fill_count += int(len(review_input))
    if candidate_count == 0:
        # candidate_budget_gate rejects an empty id column. Under the tiered
        # review-budget contract, zero verified support may still produce
        # a useful Wet membership when the frozen interval contains tier2/tier3
        # review material. No support is fabricated: the rows remain explicitly
        # interval_review_coverage and the handoff is marked degraded.
        gate_dir.mkdir(parents=True, exist_ok=True)
        if reviewable_fill_count:
            gate_result = {
                "candidate_count": 0,
                "status": "within_wet_budget",
                "action": "freeze_review_coverage_membership",
                "adaptive_refinement_required": False,
                "review_required": False,
            }
            decision_reason = "zero_supported_candidates_with_reviewable_interval_fill"
        else:
            gate_result = {
                "candidate_count": 0,
                "status": "no_supported_candidates",
                "action": "stop_review_required",
                "adaptive_refinement_required": False,
                "review_required": True,
            }
            decision_reason = "candidate_set_empty_without_reviewable_interval_fill"
        gate_manifest_path = gate_dir / "candidate_budget_manifest.json"
        gate_manifest_path.write_text(
            json.dumps(
                {
                    "schema_version": SCHEMA_VERSION,
                    "case_id": case_id,
                    "candidate_count": 0,
                    "wet_budget": wet_budget,
                    "budget_policy": budget_policy,
                    **gate_result,
                    "decision_reason": decision_reason,
                    "reviewable_interval_fill_candidate_count": reviewable_fill_count,
                    "rank_contribution": "none",
                },
                ensure_ascii=False, indent=2, sort_keys=True,
            ) + "\n",
            encoding="utf-8",
        )
        gate_result["manifest"] = str(gate_manifest_path)
        gate_result["candidate_set"] = str(union_path)
    else:
        gate_result = candidate_budget_gate(
            str(union_path),
            str(gate_dir),
            wet_budget=wet_budget,
            adaptive_round=adaptive_round,
            max_adaptive_rounds=max_adaptive_rounds,
            budget_policy=budget_policy,
            case_id=case_id,
        )
    # Strict capacity applies to the complete Wet handoff, including explicit
    # review candidates. Never truncate them silently; escalate through the
    # existing refinement/review-required states instead.
    total_reviewable_count = candidate_count + review_candidate_count
    if (budget_policy == "strict" and review_candidate_count > 0
            and total_reviewable_count > wet_budget):
        gate_result.update({
            "status": "review_required_over_budget",
            "action": "stop_review_required",
            "adaptive_refinement_required": False,
            "review_required": True,
            "decision_reason": "supported_and_route_review_candidates_exceed_strict_wet_budget",
        })
        gate_manifest_path = runtime_path(str(gate_result.get("manifest", "")))
        if gate_manifest_path.exists():
            gate_payload = json.loads(gate_manifest_path.read_text(encoding="utf-8"))
            gate_payload.update({key: gate_result[key] for key in (
                "status", "action", "adaptive_refinement_required", "review_required", "decision_reason")})
            gate_payload["candidate_count_including_route_review"] = total_reviewable_count
            gate_manifest_path.write_text(json.dumps(gate_payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    status = str(gate_result.get("status", ""))
    # Sealed Direct-C review rows are an explicit rescue lane. They may be
    # appended to an over-target supported handoff, while the over-target
    # state remains visible and A/B refinement is not allowed to remove them.
    if ((direct_c_review_count > 0 or direct_c_supported_count > 0)
            and status in {"adaptive_refinement_required", "review_required_over_budget"}):
        gate_result.update({
            "status": "above_wet_target",
            "action": "freeze_candidate_set_with_protected_direct_c_review",
            "adaptive_refinement_required": False,
            "review_required": False,
            "decision_reason": "retain_supported_and_protected_direct_c_members_above_operational_target",
        })
        status = "above_wet_target"
        gate_manifest_path = runtime_path(str(gate_result["manifest"]))
        gate_payload = json.loads(gate_manifest_path.read_text(encoding="utf-8"))
        gate_payload.update({key: gate_result[key] for key in (
            "status", "action", "adaptive_refinement_required", "review_required", "decision_reason")})
        gate_payload.update(capacity_exceeded=True,
                            candidate_count_including_route_review=total_reviewable_count,
                            direct_c_review_count=direct_c_review_count,
                            direct_c_supported_count=direct_c_supported_count)
        gate_manifest_path.write_text(json.dumps(gate_payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
                                      encoding="utf-8")
    stages.append({"stage": "candidate_budget_gate", "state": str(gate_result.get("status"))})
    frozen_candidate_set = str(gate_result.get("candidate_set") or union_path)

    # --- 3. bounded refinement, only when the gate asked for it ----------
    adaptive_result: dict[str, object] | None = None
    heavy_stage_invoked = False
    heavy_stage_reason = "not_invoked_within_wet_budget"
    deterministic_refinement_summary: dict[str, object] | None = None
    if status == "adaptive_refinement_required" or (
        status == "above_wet_target" and deterministic_subfamily_refinement and family_copy_arg
        and adaptive_round < max_adaptive_rounds
    ):
        relationships_path = runtime_path(family_copy_arg).with_name(
            "family_copy_relationships.tsv"
        ) if family_copy_arg else None
        if deterministic_subfamily_refinement and family_copy_arg:
            # One bounded, fully machine-verified pass: curated cross-subfamily
            # verdicts against the verified Route A target are applied through
            # the existing refinement contract, then the gate re-runs once.
            verification_dir = out / "subfamily_refinement_verification"
            verification = build_subfamily_refinement_verification(
                frozen_candidate_set,
                str(family_copy_arg),
                str(relationships_path) if relationships_path and relationships_path.exists() else None,
                str(verification_dir),
                case_id=case_id,
            )
            refine_dir = out / "adaptive_refinement"
            refined = apply_verified_refinement(
                frozen_candidate_set,
                str(verification["verification_table"]),
                str(refine_dir),
                case_id=case_id,
            )
            lineage = verify_adaptive_refinement(
                str(refined["refined_candidate_set"]),
                str(refined["manifest"]),
                frozen_candidate_set,
                case_id=case_id,
            )
            regate_dir = out / "candidate_budget_round1"
            gate_result = candidate_budget_gate(
                str(refined["refined_candidate_set"]),
                str(regate_dir),
                wet_budget=wet_budget,
                adaptive_round=adaptive_round + 1,
                max_adaptive_rounds=max_adaptive_rounds,
                budget_policy=budget_policy,
                case_id=case_id,
            )
            frozen_candidate_set = str(gate_result.get("candidate_set") or refined["refined_candidate_set"])
            status = str(gate_result.get("status", ""))
            candidate_count = int(refined["output_candidate_count"])
            deterministic_refinement_summary = {
                "verification_manifest": str(verification["manifest"]),
                "verification_table": str(verification["verification_table"]),
                "verdict_counts": verification["verdict_counts"],
                "removed_candidate_count": refined["removed_candidate_count"],
                "input_candidate_count": refined["input_candidate_count"],
                "output_candidate_count": refined["output_candidate_count"],
                "regate_status": status,
                "lineage": lineage,
            }
            heavy_stage_reason = "not_invoked_deterministic_subfamily_refinement"
            stages.append({
                "stage": "subfamily_refinement",
                "state": "applied",
                "removed_candidate_count": refined["removed_candidate_count"],
            })
            stages.append({"stage": "candidate_budget_gate_round1", "state": status})
        else:
            if deterministic_subfamily_refinement:
                heavy_stage_reason = (
                    "deterministic_subfamily_refinement_unavailable_no_family_copy_input"
                )
            adaptive_dir = out / "adaptive_synthesis"
            adaptive_result = build_adaptive_synthesis_request(
                frozen_candidate_set,
                str(gate_result["manifest"]),
                str(adaptive_dir),
                ambiguity_arg=ambiguity_arg,
                case_id=case_id,
            )
            heavy_stage_invoked = True
            if not deterministic_subfamily_refinement:
                heavy_stage_reason = str(
                    adaptive_result.get("invocation_reason", "candidate_count_exceeds_wet_budget")
                )
            stages.append({"stage": "adaptive_synthesis_request", "state": "requested"})
    elif status == "review_required_over_budget":
        heavy_stage_reason = "not_invoked_adaptive_round_budget_exhausted"
    elif status == "no_supported_candidates":
        heavy_stage_reason = "not_invoked_no_supported_candidates"
    elif status == "above_wet_target":
        heavy_stage_reason = "not_invoked_target_is_advisory"

    # --- 4. freeze canonical membership before any ordering ---------------
    wet_candidates: str | None = None
    wet_candidates_manifest: str | None = None
    membership_counts: dict[str, object] = {}
    membership_state = "not_emitted"
    emit_empty_scoped_inventory = (status == "no_supported_candidates" and relevance_selection
                                   and versioned_review_input and not _read_optional(interval_atlas_arg).empty)
    if status in {"within_wet_budget", "above_wet_target"} or emit_empty_scoped_inventory:
        wet_candidates, wet_candidates_manifest = _freeze_canonical_wet_membership(
            frozen_candidate_set,
            out,
            wet_budget=wet_budget,
            case_id=case_id,
            atlas_arg=interval_atlas_arg,
            reference_function_arg=reference_function_arg,
            route_c_members_arg=route_c_members_arg,
            route_c_manifest_arg=route_c_manifest_arg,
            review_candidates_arg=review_candidates_arg,
            review_candidates_manifest_arg=review_candidates_manifest_arg,
            fill_mode=wet_fill_mode,
            relevance_selection=relevance_selection,
            family_copy_arg=family_copy_arg,
            provisional_screen_ledger_arg=provisional_screen_ledger_arg,
            provisional_screen_ledger_sha256=provisional_screen_ledger_sha256,
            identity_assessment_arg=identity_assessment_arg,
            identity_target_proteins_arg=identity_target_proteins_arg,
        )
        membership_state = "frozen"
        membership_counts = json.loads(
            runtime_path(wet_candidates_manifest).read_text(encoding="utf-8"))
        stages.append({"stage": "canonical_wet_membership", "state": "frozen"})
    elif status == "adaptive_refinement_required":
        membership_state = "not_emitted_over_budget_refinement_required"
    elif status == "review_required_over_budget":
        membership_state = "not_emitted_over_budget_review_required"
    elif status == "no_supported_candidates":
        membership_state = "not_emitted_no_supported_candidates"

    # Count the actual Wet set, including review candidates, after membership
    # freezes. The second-stage trigger is independent of the ~20 fill target.
    second_stage_state = "not_evaluated_no_wet_membership"
    second_stage_gate = None
    boundary_scope_manifest = membership_counts.get("boundary_scope_manifest")
    boundary_scope_sha256 = membership_counts.get("boundary_scope_sha256")
    if wet_candidates and boundary_scope_manifest:
        from .boundary_tier_review import build_boundary_review_request
        boundary_request = build_boundary_review_request(
            str(boundary_scope_manifest), str(out / "second_stage_request"),
            scope_sha256=str(boundary_scope_sha256), case_id=case_id,
            adaptive_round=adaptive_round, max_adaptive_rounds=max_adaptive_rounds,
            comparison_bundle_arg=boundary_route_b_bundle_arg)
        second_stage_state = str(boundary_request["state"])
        if boundary_request.get("request"):
            adaptive_result = boundary_request
            heavy_stage_invoked = True
            heavy_stage_reason = "boundary_relevance_tier_requires_comparison"
        stages.append({"stage": "second_stage_review", "state": second_stage_state})
    elif wet_candidates and (direct_c_review_count or direct_c_supported_count):
        # Legacy capacity modes have no protected boundary scope. Do not send
        # the complete rescued set to a count-only reducer which could defer C.
        # Capacity excess remains explicit in the budget manifest above.
        second_stage_state = "skipped_protected_direct_c_retained"
        stages.append({"stage": "second_stage_review", "state": second_stage_state})
    elif wet_candidates and membership_counts.get("candidate_count", 0):
        second_stage_gate = candidate_budget_gate(
            wet_candidates, str(out / "second_stage_gate"),
            wet_budget=SECOND_STAGE_THRESHOLD, budget_policy="strict",
            adaptive_round=adaptive_round, max_adaptive_rounds=max_adaptive_rounds,
            case_id=case_id,
        )
        if second_stage_gate["status"] == "adaptive_refinement_required":
            adaptive_result = build_adaptive_synthesis_request(
                wet_candidates, str(second_stage_gate["manifest"]),
                str(out / "second_stage_request"), case_id=case_id,
            )
            second_stage_state = "requested"
            heavy_stage_invoked = True
            heavy_stage_reason = "wet_candidate_count_exceeds_25"
        elif second_stage_gate["status"] == "review_required_over_budget":
            second_stage_state = "review_required_round_budget_exhausted"
        else:
            second_stage_state = "skipped_at_or_below_25"
        stages.append({"stage": "second_stage_review", "state": second_stage_state})

    # --- 5. ordering, only inside the already frozen Wet membership -------
    narrowing_result: dict[str, object] | None = None
    narrowing_state = "not_run"
    ordering_ranked_input: str | None = None
    ordering_context_inputs: dict[str, str | None] = {
        "family_recovery": None,
        "family_copy": None,
        "family_copy_shortlist": None,
    }
    ordering_skip_reason = ""
    if wet_candidates and ranked_arg:
        try:
            ordering_ranked_input = _restrict_ranked_to_wet_membership(
                ranked_arg,
                wet_candidates,
                out,
            )
        except ValueError as error:
            # Ordering may only reorder the exact frozen membership. A ranked
            # table that predates the membership (for example one built before
            # review slots were filled, #570) cannot do that, so ordering is
            # declared not evaluated instead of being run on a partial set or
            # allowed to change membership.
            ordering_ranked_input = None
            ordering_skip_reason = str(error)
    if ordering_ranked_input:
        ordering_context_inputs = _restrict_narrowing_context_to_wet_membership(
            wet_candidates,
            out,
            family_recovery_arg=family_recovery_arg,
            family_copy_arg=family_copy_arg,
            family_copy_shortlist_arg=family_copy_shortlist_arg,
        )
        narrowing_dir = out / "integrated_narrowing"
        narrowing_result = integrated_narrowing(
            ordering_ranked_input,
            str(narrowing_dir),
            family_recovery_arg=ordering_context_inputs["family_recovery"],
            convergence_arg=convergence_arg,
            family_copy_arg=ordering_context_inputs["family_copy"],
            family_copy_shortlist_arg=ordering_context_inputs["family_copy_shortlist"],
            remaining_arg=None,
            case_id=case_id,
        )
        narrowed = pd.read_csv(
            runtime_path(str(narrowing_result["integrated_narrowing"])),
            sep="\t", dtype=str,
        ).fillna("")
        narrowed_ids = set(narrowed["gene_id"].astype(str).str.strip())
        wet_ids = _wet_membership_ids(wet_candidates)
        if narrowed_ids != wet_ids:
            raise ValueError("ordering layer changed canonical Wet membership (fail-closed)")
        narrowing_state = "ordering_only"
        stages.append({"stage": "integrated_narrowing", "state": "ordering_only"})
    elif wet_candidates and ordering_skip_reason:
        narrowing_state = "not_evaluated_membership_not_covered_by_ranked_input"
    elif wet_candidates:
        narrowing_state = "not_evaluated_no_ranked_input"

    # --- 6. one verdict for the whole chain -------------------------------
    if status in {"within_wet_budget", "above_wet_target"}:
        # Supplied evidence coverage and strict support are separate. Generic
        # interval fill retains its historical degraded label; validated initial
        # candidate evidence need not establish strict support.
        handoff_state = (
            "wet_reviewable_degraded_coverage"
            if (candidate_count == 0 and initial_b_candidate_count == 0
                and not (initial_a_coverage or {}).get("usable_candidate_count")
                and not (route_c_coverage or {}).get("record_state_counts", {}).get("valid_inside_support_or_review"))
                or coverage_state in DEGRADED_COVERAGE_STATES
            else "wet_reviewable"
        )
    elif status == "adaptive_refinement_required":
        handoff_state = "adaptive_refinement_requested"
    elif status == "no_supported_candidates":
        handoff_state = "review_required_no_supported_candidates"
    else:
        handoff_state = "review_required_over_budget"
    if relevance_selection and second_stage_state in {"requested", "review_required_round_budget_exhausted"}:
        handoff_state = ("boundary_review_requested" if second_stage_state == "requested"
                         else "review_required_round_budget_exhausted")
    if (membership_counts.get("provisional_wet") or {}).get("eligible"):
        handoff_state = "provisional_adopted_wet"
    report = None
    count_populations = {"supported_union": int(len(union_frame)),
                         "eligible_first_stage": membership_counts.get("eligible_first_stage_candidate_count"),
                         "final_wet": membership_counts.get("candidate_count") if wet_candidates else None,
                         "adopted": None, "pending": None}
    if "active_first_stage_candidate_count" in membership_counts:
        count_populations["active_first_stage"] = membership_counts["active_first_stage_candidate_count"]
    route_a_rethink = route_a_rethink_signal(
        membership_counts.get("candidate_count") if wet_candidates else 0,
        min_route_candidates=min_route_candidates,
    )
    if wet_candidates:
        from .boundary_tier_review import write_wet_review_report
        pending_ids = set()
        if boundary_scope_manifest:
            scope = json.loads(runtime_path(str(boundary_scope_manifest)).read_text(encoding="utf-8"))
            pending = pd.read_csv(scope["artifacts"]["boundary"]["path"], sep="\t", dtype=str).fillna("")
            pending_ids = set(pending.get("candidate_gene_id", []))
        report = write_wet_review_report(
            wet_candidates, str(out / "wet_review_report.md"), case_id, handoff_state, "first_stage",
            adopted_gene_ids=_wet_membership_ids(wet_candidates) - pending_ids, pending_gene_ids=pending_ids,
            route_a_rethink=route_a_rethink,
            selection_policy_version=str(membership_counts.get("first_stage_selection_policy_version") or ""))
        count_populations.update(adopted=report["adopted_candidate_count"], pending=report["pending_candidate_count"])
        report_path = runtime_path(report["path"])
        with report_path.open("a", encoding="utf-8") as handle:
            handle.write("\nCandidate count populations (unique genes): " + json.dumps(count_populations, sort_keys=True) + "\n")
            handle.write("Input evidence coverage: " + coverage_state + "; execution outside the supplied import/source artifacts is not assessed by handoff.\n")
        report["sha256"] = _sha256_path(report_path)
        report["candidate_count_populations"] = count_populations

    manifest_policy_is_evidence_scoped = relevance_selection and (
        membership_counts.get("first_stage_selection_policy_version") in EVIDENCE_SCOPED_POLICIES
        or (membership_counts.get("first_stage_selection_policy_version") is None
            and not provisional_screen_ledger_arg)
    )
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id,
        "handoff_state": handoff_state,
        "report": report,
        "producer_class": "canonical_wet_membership_gate",
        "candidate_count": candidate_count,
        "candidate_count_basis": "supported_union_after_legacy_refinement",
        "candidate_count_populations": count_populations,
        "wet_budget": wet_budget,
        "budget_policy": budget_policy,
        "above_target_count": max(0, candidate_count - wet_budget),
        "adaptive_round": adaptive_round,
        "max_adaptive_rounds": max_adaptive_rounds,
        "budget_status": status,
        "second_stage_threshold": SECOND_STAGE_THRESHOLD,
        "second_stage_state": second_stage_state,
        "second_stage_gate_manifest": second_stage_gate["manifest"] if second_stage_gate else None,
        "candidate_id_column": id_column,
        "heavy_stage_invoked": heavy_stage_invoked,
        "heavy_stage_invocation_reason": heavy_stage_reason,
        "membership_state": membership_state,
        "final_wet_candidate_set": wet_candidates,
        "final_wet_candidate_manifest": wet_candidates_manifest,
        # The frozen membership can be larger than the supported union when
        # review slots were filled, so both numbers are reported.
        "final_wet_candidate_count": membership_counts.get("candidate_count") if wet_candidates else None,
        "final_supported_candidate_count": membership_counts.get("supported_candidate_count") if wet_candidates else None,
        "final_review_coverage_count": membership_counts.get("review_coverage_candidate_count") if wet_candidates else None,
        "final_support_tier_counts": membership_counts.get("support_tier_counts") if wet_candidates else None,
        "route_a_rethink": route_a_rethink,
        "wet_fill_mode": wet_fill_mode,
        "review_candidates_input": _display_input(review_candidates_arg),
        "review_candidates_manifest_input": _display_input(review_candidates_manifest_arg),
        "narrowing_state": narrowing_state,
        "stage_a_narrowed_producer_class": "ordering_only" if narrowing_result else None,
        "ordering_ranked_input": ordering_ranked_input,
        "ordering_skip_reason": ordering_skip_reason or None,
        "ordering_context_inputs": ordering_context_inputs,
        "ordering_membership_preserved": bool(narrowing_result),
        "ordering_policy": (
            "integrated narrowing may order only the already frozen wet_candidates.tsv set; "
            "the ranked table and every gene-bearing narrowing context input are first projected "
            "onto that same frozen set. This projection cannot select membership; it only removes "
            "context for genes already outside production membership. Ordering cannot create, add, "
            "remove, pad, or truncate production membership. first_wet_priority remains an auditable "
            "pass-through field; current ordering is narrowing tier then existing deterministic/base rank tie-break."
        ),
        "input_coverage_state": coverage_state,
        "input_coverage_detail": coverage_detail,
        "reduction_ratio": (
            round(1.0 - (candidate_count / float(len(union_frame))), 4) if len(union_frame) else None
        ),
        "candidate_set": str(union_path),
        "supported_member_union_sha256": _sha256_path(union_path),
        "supported_member_union_manifest": (
            str(out / "supported_member_union" / "supported_member_union_manifest.json")
            if (out / "supported_member_union" / "supported_member_union_manifest.json").exists() else None
        ),
        "supported_member_union_manifest_sha256": (
            _sha256_path(out / "supported_member_union" / "supported_member_union_manifest.json")
            if (out / "supported_member_union" / "supported_member_union_manifest.json").exists() else None
        ),
        "frozen_candidate_set": frozen_candidate_set,
        "refinement_lineage": union_result.get("refinement_lineage"),
        "candidate_budget_manifest": str(gate_result.get("manifest", "")),
        "adaptive_synthesis_request": str(adaptive_result.get("request")) if adaptive_result else None,
        "boundary_scope_manifest": boundary_scope_manifest,
        "boundary_scope_sha256": boundary_scope_sha256,
        "provisional_wet": membership_counts.get("provisional_wet"),
        "interval_unresolved": membership_counts.get("interval_unresolved"),
        "interval_unresolved_sha256": membership_counts.get("interval_unresolved_sha256"),
        "relevance_selection_applied": relevance_selection,
        "first_stage_selection_policy_version": membership_counts.get("first_stage_selection_policy_version"),
        "eligible_first_stage_candidate_count": membership_counts.get("eligible_first_stage_candidate_count"),
        **({"active_first_stage_candidate_count": membership_counts["active_first_stage_candidate_count"],
            "upper_priority_applied": membership_counts["upper_priority_applied"]}
           if "active_first_stage_candidate_count" in membership_counts else {}),
        "deterministic_subfamily_refinement": deterministic_refinement_summary,
        "stage_a_narrowed": str(narrowing_result.get("integrated_narrowing")) if narrowing_result else None,
        "stages": stages,
        "rank_contribution": "none",
        "policy": [
            "wet_candidates.tsv is the sole canonical production membership artifact",
            ("explicit provisional Wet adopts every concrete T1-T3 and protected C member after complete screen accounting; identity-only rows remain not_evaluated"
             if (membership_counts.get("provisional_wet") or {}).get("eligible") else
             "upper_evidence_scoped_v3 uses reported T1-T3 and protected C when upper evidence exists; compare only the first overflowing tier and preserve higher tiers; later eligibility remains in the ledger without padding, within-tier ranking, forced reduction or automatic continuation"
             if membership_counts.get("first_stage_selection_policy_version") == UPPER_EVIDENCE_SCOPED_POLICY else
             "ordinary target handoff retains T1 and existing eligible union/review candidates; if over 25, it compares the first overflowing tier and all lower eligible seeds together, excluding identity-only inventory"
             if manifest_policy_is_evidence_scoped else
             "ordinary target+atlas handoff selects whole relevance tiers; T1 is retained, only the first overflowing tier is compared, and T4/U are rescue pools"
             if relevance_selection else "legacy strict/supported-only handoff retains its verified refinement and tier2/tier3 fill contract"),
            "review-coverage fill never creates verified support, score, or route evidence",
            "sealed Route A/B review and supported claims reach relevance separately; review membership never upgrades support",
            "an over-budget set is never truncated to satisfy a Top-N file contract",
            "stage_a/ is a baseline/review-budget artifact, not canonical production membership",
            "stage_a_narrowed/integrated narrowing is ordering-only and receives exactly the frozen Wet membership",
            "gene-bearing narrowing context is projected to frozen membership before fail-closed ordering validation",
            "an adaptive/heavy stage runs only when the gate requested it",
            ("v3 compares only the first overflowing tier without within-tier ranking or forced reduction; higher tiers are retained and later tiers remain in the ledger"
             if membership_counts.get("first_stage_selection_policy_version") == UPPER_EVIDENCE_SCOPED_POLICY else
             "the evidence-scoped boundary compares eligible candidates after the first overflow; request generation is not model execution, and original Wet membership is retained"
             if manifest_policy_is_evidence_scoped else
             "the frozen Wet total including review candidates triggers second-stage review only above 25; request generation is not model execution, and membership is retained"),
            "deterministic subfamily refinement removes only curated, target-relative, sourced cross-subfamily support and re-gates once; it never truncates to fit the budget",
            "exit code 0 is not completion; handoff_state carries the verdict",
        ],
    }
    manifest_path = out / "wet_candidate_handoff_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    return {
        "manifest": str(manifest_path),
        "handoff_state": handoff_state,
        "report": report,
        "route_a_rethink": route_a_rethink,
        "budget_status": status,
        "second_stage_threshold": SECOND_STAGE_THRESHOLD,
        "second_stage_state": second_stage_state,
        "candidate_count": candidate_count,
        "heavy_stage_invoked": heavy_stage_invoked,
        "heavy_stage_invocation_reason": heavy_stage_reason,
        "input_coverage_state": coverage_state,
        "frozen_candidate_set": frozen_candidate_set,
        "wet_candidates": wet_candidates,
        "wet_candidates_manifest": wet_candidates_manifest,
        "membership_state": membership_state,
        "adaptive_synthesis_request": str(adaptive_result.get("request")) if adaptive_result else None,
        "boundary_scope_manifest": boundary_scope_manifest,
        "boundary_scope_sha256": boundary_scope_sha256,
        "provisional_wet": membership_counts.get("provisional_wet"),
        "interval_unresolved": membership_counts.get("interval_unresolved"),
        "interval_unresolved_sha256": membership_counts.get("interval_unresolved_sha256"),
        "adaptive_round": adaptive_round,
        "max_adaptive_rounds": max_adaptive_rounds,
        "relevance_selection_applied": relevance_selection,
        "first_stage_selection_policy_version": membership_counts.get("first_stage_selection_policy_version"),
        "stage_a_narrowed": str(narrowing_result.get("integrated_narrowing")) if narrowing_result else None,
        "narrowing_state": narrowing_state,
        "ordering_skip_reason": ordering_skip_reason or None,
    }
