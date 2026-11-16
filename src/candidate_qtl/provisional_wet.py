"""Truth-blind eligibility for an evidence-limited provisional Wet decision."""

from __future__ import annotations

import json
from collections import defaultdict
from collections.abc import Iterable, Mapping
from typing import Any

from .relevance_tiers import T1, T2, T3, T4, UNRESOLVED

CONCRETE_TIERS = frozenset((T1, T2, T3))
IDENTITY_TIERS = frozenset((T4, UNRESOLVED))
EXPLICIT_GAP_STATES = frozenset(("not_evaluated_material_gap", "not_evaluated_provider_failure"))


def assess_provisional_wet(
    inventory: Iterable[Mapping[str, Any]], plan: Mapping[str, Any],
    screen_rows: Iterable[Mapping[str, Any]], *, case_id: str,
    protected_gene_ids: Iterable[str] = (), upper: int = 25,
    coverage_inventory: Iterable[Mapping[str, Any]] | None = None,
    agreed_deferred_gene_ids: Iterable[str] = (),
) -> dict[str, Any]:
    """Retain concrete hypotheses when a fully accounted identity tier is unreadable.

    The screen ledger proves observation coverage, not negative biology. No
    ranking or within-tier selection is performed here.
    """
    rows = [dict(row) for row in inventory]
    ids = [str(row.get("candidate_gene_id") or row.get("gene_id") or "").strip() for row in rows]
    if not ids or any(not gene for gene in ids) or len(ids) != len(set(ids)):
        raise ValueError("provisional Wet requires a unique complete interval roster")
    by_id = dict(zip(ids, rows))
    coverage_rows = rows if coverage_inventory is None else [dict(row) for row in coverage_inventory]
    coverage_ids = [str(row.get("candidate_gene_id") or row.get("gene_id") or "").strip() for row in coverage_rows]
    if (not coverage_ids or any(not gene for gene in coverage_ids)
            or len(coverage_ids) != len(set(coverage_ids))):
        raise ValueError("provisional Wet requires a unique original coverage roster")
    original = dict(zip(coverage_ids, coverage_rows))
    deferred = list(agreed_deferred_gene_ids)
    if (len(deferred) != len(set(deferred)) or set(coverage_ids) - set(ids) != set(deferred)
            or not set(ids) <= set(coverage_ids) or any(original[gene] != row for gene, row in by_id.items())):
        raise ValueError("active roster must preserve original rows and exactly account agreed deferrals")
    protected = set(protected_gene_ids)
    if not protected <= set(ids):
        raise ValueError("protected Wet candidates are outside the interval roster")
    adopted = {gene for gene, row in by_id.items() if row.get("relevance_tier") in CONCRETE_TIERS} | protected
    if not adopted:
        return {"eligible": False, "reason": "no_concrete_hypothesis"}
    if len(adopted) > upper:
        return {"eligible": False, "reason": "concrete_or_protected_over_upper"}
    if plan.get("boundary_tier") not in IDENTITY_TIERS or not plan.get("boundary"):
        return {"eligible": False, "reason": "no_identity_only_boundary"}
    if any(row.get("relevance_tier") not in IDENTITY_TIERS for row in plan["boundary"]):
        return {"eligible": False, "reason": "informative_boundary"}
    for gene, row in by_id.items():
        if gene in adopted or row.get("relevance_tier") not in IDENTITY_TIERS:
            continue
        try:
            evidence = json.loads(str(row.get("relevance_evidence_json") or "[]"))
        except json.JSONDecodeError as error:
            raise ValueError(f"invalid relevance evidence for {gene}") from error
        if not isinstance(evidence, list):
            raise ValueError(f"relevance evidence must be a list for {gene}")
        if (evidence or any(str(row.get(field, "")).strip() for field in
                            ("review_claim_ids", "review_hypothesis_ids", "review_source_ids"))
                or any(str(row.get(field, "")).strip().casefold() in {"true", "1", "yes"}
                       for field in ("route_a_supported", "route_b_supported", "route_c_supported"))):
            return {"eligible": False, "reason": "identity_tier_has_comparison_evidence", "gene_id": gene}

    covered: dict[str, set[str]] = defaultdict(set)
    seen_units: set[str] = set()
    screen = [dict(row) for row in screen_rows]
    if not screen:
        return {"eligible": False, "reason": "screen_coverage_missing"}
    for row in screen:
        unit = str(row.get("unit_id", "")).strip()
        if str(row.get("case_id", "")).strip() != case_id or not unit or unit in seen_units:
            raise ValueError("screen ledger case or unit roster mismatch")
        seen_units.add(unit)
        members = [gene.strip() for gene in str(row.get("member_interval_genes", "")).split(";") if gene.strip()]
        if not members or len(members) != len(set(members)) or not set(members) <= set(coverage_ids):
            raise ValueError(f"screen unit members are outside the interval roster: {unit}")
        state = str(row.get("screen_state", "")).strip()
        if state not in {"screened", *EXPLICIT_GAP_STATES}:
            return {"eligible": False, "reason": "screen_not_complete", "unit_id": unit}
        if state == "screened" and not str(row.get("answer_sha256", "")).strip():
            return {"eligible": False, "reason": "screen_answer_not_frozen", "unit_id": unit}
        for gene in members:
            covered[gene].add(state)
    if set(covered) != set(coverage_ids):
        return {"eligible": False, "reason": "screen_gene_roster_incomplete",
                "missing_gene_count": len(set(coverage_ids) - set(covered))}
    material_by_gene = {
        gene: ("screened_with_gap" if "screened" in states and len(states) > 1 else
               "screened" if "screened" in states else "not_evaluated_material_gap")
        for gene, states in covered.items()
    }
    return {
        "eligible": True, "reason": "identity_boundary_has_no_comparison_evidence",
        "adopted_gene_ids": sorted(adopted), "unresolved_gene_ids": sorted(set(coverage_ids) - adopted),
        "material_by_gene": material_by_gene,
        "screen_unit_count": len(screen),
        "screened_unit_count": sum(row.get("screen_state") == "screened" for row in screen),
        "material_gap_unit_count": sum(row.get("screen_state") in EXPLICIT_GAP_STATES for row in screen),
    }
