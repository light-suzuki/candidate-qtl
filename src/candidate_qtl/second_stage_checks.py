"""Non-ranking review views over fresh alignments and native family context.

Numerical profiles are adjustable inspection heuristics, NOT validated gene
exclusion rules. Native family/target verdicts are reported, not re-authorized.
No reference, isoform, candidate or first-stage evidence is removed here.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

CONTRACT = "second-stage-review-checks-1.0"
DEFAULT_PROFILES = (
    {"profile_id": "broad_alignment_review", "min_identity": 30.0,
     "min_query_coverage": 0.5, "min_subject_coverage": 0.5, "max_evalue": 1e-5,
     "basis": "illustrative_review_heuristic_not_calibrated"},
    {"profile_id": "strict_full_length_review", "min_identity": 50.0,
     "min_query_coverage": 0.8, "min_subject_coverage": 0.8, "max_evalue": 1e-10,
     "basis": "illustrative_review_heuristic_not_calibrated"},
)
PROFILE_FIELDS = frozenset(DEFAULT_PROFILES[0])
ATLAS_FIELDS = ("gene_id", "family_id", "family_source", "discovery_bucket_id",
                "discovery_bucket_method", "evolutionary_family_id", "evolutionary_family_source",
                "evolutionary_family_status", "subfamily_id", "subfamily_source", "domain_ids",
                "domain_architecture", "domain_source", "orthogroup", "orthology_relationship",
                "orthology_providers", "copy_set_id", "copy_set_status", "focal_paralog_count",
                "gene_model_completeness", "classification_basis")
COPY_FIELDS = ("candidate_gene_id", "family_id", "copy_set_id", "copy_set_status", "subfamily",
               "orthology_class", "paralog_ambiguity", "disambiguation_status", "review_reason",
               "target_id", "target_reference_gene", "target_verification_status", "target_subfamily",
               "target_domain_architecture", "subfamily_vs_target", "domain_vs_target",
               "domain_vs_family_modal", "domain_architecture_match")
# Same reported-status vocabulary as family_to_copy; not a substitute for its
# source binding/verification. Every derived state below is explicitly reported.
REPORTED_TARGET_STATUSES = frozenset({"verified", "source_backed", "curated"})


def load_profiles(path: Path | None) -> list[dict[str, Any]]:
    profiles = json.loads(path.read_text(encoding="utf-8")) if path else [dict(p) for p in DEFAULT_PROFILES]
    if not isinstance(profiles, list) or not profiles:
        raise ValueError("review profiles must be a non-empty JSON list")
    seen = set()
    for p in profiles:
        if not isinstance(p, dict) or set(p) != PROFILE_FIELDS:
            raise ValueError("review profile requires exactly profile_id/basis and four numerical conditions")
        for key in ("profile_id", "basis"):
            if not isinstance(p[key], str) or not p[key].strip():
                raise ValueError(f"review profile requires a non-empty {key}")
        if p["profile_id"] in seen:
            raise ValueError("review profile IDs must be unique")
        seen.add(p["profile_id"])
        for key, upper in (("min_identity", 100), ("min_query_coverage", 1),
                           ("min_subject_coverage", 1), ("max_evalue", math.inf)):
            v = p[key]
            if type(v) not in (int, float) or not math.isfinite(v) or not 0 <= v <= upper:
                raise ValueError(f"invalid review profile {key}")
    return sorted(profiles, key=lambda p: p["profile_id"])


def alignment_checks(row: dict[str, Any], profiles: list[dict[str, Any]],
                     search_evalue: float) -> list[dict[str, Any]]:
    results = []
    for p in profiles:
        values = (("query_residue_identity", row["query_residue_identity"] >= p["min_identity"]),
                  ("subject_residue_identity", row["subject_residue_identity"] >= p["min_identity"]),
                  ("query_coverage", row["query_coverage"] >= p["min_query_coverage"]),
                  ("subject_coverage", row["subject_coverage"] >= p["min_subject_coverage"]),
                  ("evalue", row["evalue"] <= p["max_evalue"]))
        unmet = [key for key, passes in values if not passes]
        results.append({"profile_id": p["profile_id"],
                        "state": "meets_review_conditions" if not unmet else "below_review_conditions",
                        "unmet_conditions": unmet,
                        "search_more_restrictive_than_review_evalue": search_evalue < p["max_evalue"],
                        "candidate_membership_effect": "none"})
    return results


def load_context(path: Path | None, role: str, wanted: set[str]) -> tuple[dict[str, list[dict[str, str]]], str]:
    if path is None:
        return {}, "not_evaluated_not_supplied"
    if not path.is_file():
        return {}, "not_evaluated_missing_file"
    fields = ATLAS_FIELDS if role == "interval_atlas" else COPY_FIELDS
    key = fields[0]
    records: dict[str, list[dict[str, str]]] = {}
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if key not in (reader.fieldnames or []):
            raise ValueError(f"{role} requires native {key}")
        for row in reader:
            gene = (row.get(key) or "").strip()
            if not gene:
                raise ValueError(f"{role} contains an empty gene ID")
            if gene in wanted:
                # Keep every distinct annotation/target, including disagreements;
                # no gene-only last-writer join, symbol inference, or rank fields.
                item = {name: (row.get(name) or "").strip() for name in fields}
                records.setdefault(gene, []).append(item)
    for gene, rows in records.items():
        unique = {json.dumps(row, sort_keys=True, ensure_ascii=False): row for row in rows}
        records[gene] = [unique[k] for k in sorted(unique)]
    return records, "available"


def family_checks(gene: str, atlas: list[dict[str, str]], copies: list[dict[str, str]],
                  input_states: dict[str, str]) -> dict[str, Any]:
    flags: set[str] = set()
    if not atlas and not copies:
        flags.add("family_not_evaluated")
    if len(atlas) > 1:
        flags.add("multiple_atlas_assignments_require_review")
    for row in atlas:
        evolutionary = bool(row.get("evolutionary_family_id"))
        if not evolutionary and (row.get("family_id") or row.get("discovery_bucket_id")):
            flags.add("coarse_bucket_is_not_evolutionary_family_verification")
        if not evolutionary and row.get("domain_ids"):
            flags.add("domain_only_is_not_functional_family_identity")
        if evolutionary and not row.get("evolutionary_family_source"):
            flags.add("family_assignment_source_missing")
        if not row.get("subfamily_id"):
            flags.add("subfamily_unresolved")
        if row.get("subfamily_source") in {"orthogroup_proxy", "symbol_stem"}:
            flags.add("subfamily_proxy_requires_independent_check")
        if row.get("gene_model_completeness") in {"partial_model", "no_model"}:
            flags.add("model_incompleteness_may_explain_alignment_or_domain_difference")
        if "conflict" in row.get("evolutionary_family_status", "").lower():
            flags.add("family_assignment_conflict_requires_review")
    target_checks = []
    for row in copies:
        target_bound = (row.get("target_verification_status") in REPORTED_TARGET_STATUSES
                        and bool(row.get("target_id") or row.get("target_reference_gene")))
        axes = {}
        for axis in ("subfamily_vs_target", "domain_vs_target"):
            value = row.get(axis, "")
            axes[axis] = ("not_evaluated_target_not_bound" if not target_bound else
                          "reported_" + value if value in {"target_match", "target_mismatch"} else
                          "not_evaluated")
        if "reported_target_mismatch" in axes.values():
            flags.add("target_relative_mismatch_requires_review_not_automatic_exclusion")
        if row.get("domain_architecture_match") in {"architecture_mismatch", "architecture_partial"}:
            flags.add("family_modal_domain_difference_is_qc_not_selection")
        if row.get("orthology_class") in {"one_to_many", "many_to_many", "close_paralog", "unmapped_ambiguous"}:
            flags.add("copy_or_orthology_ambiguity_must_remain")
        target_checks.append({"family_id": row.get("family_id", ""),
                              "target_id": row.get("target_id", ""),
                              "target_reference_gene": row.get("target_reference_gene", ""),
                              "reported_target_verification_status": row.get("target_verification_status", ""),
                              **axes})
    return {"gene_id": gene, "state": "context_available" if atlas or copies else "not_evaluated",
            "input_states": input_states, "atlas_context": atlas, "copy_context": copies,
            "target_relative_checks": target_checks, "review_flags": sorted(flags),
            "independent_source_verification_performed": False,
            "new_family_assignment_performed": False,
            "candidate_membership_effect": "none"}


def review_checks(ids: list[str], reference_rows: list[dict[str, Any]], cohort_rows: list[dict[str, Any]],
                  profiles: list[dict[str, Any]], *, search_evalue: float,
                  atlas_path: Path | None = None, copy_path: Path | None = None) -> dict[str, Any]:
    """Annotate already-measured hits and native context without selecting genes."""
    atlas, astate = load_context(atlas_path, "interval_atlas", set(ids))
    copies, cstate = load_context(copy_path, "family_copy", set(ids))
    conditions = []
    for scope, rows in (("reference", reference_rows), ("cohort", cohort_rows)):
        for row in rows:
            row["review_conditions"] = alignment_checks(row, profiles, search_evalue)
            for check in row["review_conditions"]:
                conditions.append({"scope": scope, "gene_id": row["gene_id"], "protein_id": row["protein_id"],
                                   "subject_id": row["subject_id"], **check})
    states = {"interval_atlas": astate, "family_copy": cstate}
    families = {g: family_checks(g, atlas.get(g, []), copies.get(g, []), states) for g in sorted(ids)}
    return {"contract": CONTRACT, "profiles": profiles, "alignment_conditions": conditions,
            "families": families, "context_states": states,
            "policy": "review_conditions_not_candidate_filters; family_context_not_reverified_support"}
