"""Plan optional Route-B literature follow-up from frozen initial assessments.

This is a research scheduling ledger, NOT a Tier, candidate or Wet selector.
A light screen records what the original researcher checked, not that other
phenotype papers do not exist. The canonical prompt/claim roster is untouched.
"""
from __future__ import annotations

from collections import Counter
import json
from pathlib import Path
from typing import Any

try:
    from .agentic_route_b_initial import CONTRACT as INITIAL_CONTRACT, digest, encoded
except ImportError:  # direct script invocation
    from agentic_route_b_initial import CONTRACT as INITIAL_CONTRACT, digest, encoded

CONTRACT = "route-b-light-tier-screen-1.0.0"
SCREEN_STATES = {
    "no_concrete_upper_lead", "concrete_upper_lead", "not_assessed", "retrieval_failed"
}
SPECIES_SCOPES = {"focal_species", "other_plants"}
UPPER_TIERS = {"T1", "T2", "T3"}
LOWER_TIERS = {"T4", "T5"}
TIERS = UPPER_TIERS | LOWER_TIERS | {"T6", "U"}
# Missing context with potential developmental or phenotype specificity must
# not be treated as already interpreted just because the function title was.
T4_SUPPORT_FIELDS = {
    "reference_function", "function_description", "reference_mutant_phenotype",
    "mutant_phenotype", "reference_pathway", "pathway", "reference_pathway_id",
    "pathway_id",
}
T5_SUPPORT_FIELDS = {
    "go_terms", "keywords", "protein_family", "protein_name",
    "reference_domain_description", "domain_description",
    "domain_architecture", "domain_ids", "reference_domain_ids",
    "local_annotation_description", "local_annotation_go_terms",
}
UPPER_RELEVANT_FIELDS = (T4_SUPPORT_FIELDS | T5_SUPPORT_FIELDS | {
    "tissue_specificity", "developmental_stage",
})


def _read_json(path: Path) -> dict:
    data = json.loads(path.read_bytes())
    if not isinstance(data, dict):
        raise ValueError("expected a JSON object")
    return data


def _text(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be nonempty text")
    return value.strip()


def _strings(value: Any, label: str, *, required: bool = False) -> list[str]:
    if not isinstance(value, list) or any(not isinstance(v, str) or not v.strip() for v in value):
        raise ValueError(f"{label} must be an array of nonempty strings")
    result = [v.strip() for v in value]
    if len(set(result)) != len(result) or (required and not result):
        raise ValueError(f"{label} must be unique and nonempty when required")
    return result


def _bound(base: Path, descriptor: dict) -> Path:
    if not isinstance(descriptor, dict):
        raise ValueError("missing input descriptor")
    path = Path(_text(descriptor.get("path"), "input path"))
    if not path.is_absolute():
        path = base / path
    path = path.resolve()
    if not path.is_file() or digest(path.read_bytes()) != descriptor.get("sha256"):
        raise ValueError("input file missing or SHA-256 mismatch")
    return path


def _gene_reference_pairs(bundle: dict) -> tuple[dict[str, str], dict[str, set[tuple[str, str]]]]:
    """Use explicit gene/reference/species bindings, never a unit-level label."""
    roster: dict[str, str] = {}
    pairs: dict[str, set[tuple[str, str]]] = {}
    for unit in bundle.get("units", []):
        uid = _text(unit.get("unit_id"), "unit_id")
        members = _strings(unit.get("member_interval_genes"), "unit members", required=True)
        genes = unit.get("genes", [])
        if not isinstance(genes, list) or {g.get("gene_id") for g in genes} != set(members):
            raise ValueError("bundle unit gene roster mismatch")
        for gene in genes:
            gid = _text(gene.get("gene_id"), "gene_id")
            if gid in roster:
                raise ValueError("duplicate bundle gene")
            roster[gid] = uid
            pairs[gid] = set()

            for ref_key, species_key in (
                ("reference_gene", "reference_species"),
                ("reference_gene_id", "reference_species"),
                ("orthology_reference_gene", "orthology_reference_species"),
            ):
                species = str(gene.get(species_key) or "").strip()
                refs = str(gene.get(ref_key) or "").split(";")
                if species and ";" not in species:
                    pairs[gid].update((ref.strip(), species.casefold()) for ref in refs if ref.strip())
        for card in unit.get("reference_material_cards", []):
            gid = str(card.get("interval_gene_id") or "")
            # Only an explicitly supplied focal gene may bind a card.
            if not gid and card.get("gene_id") in members:
                gid = str(card["gene_id"])
            # An unbound card or a provider's reference gene_id is NOT a focal member.
            if gid not in members:
                continue
            ref = str(card.get("reference_gene") or card.get("reference_gene_id") or "").strip()
            species = str(card.get("reference_species") or "").strip()
            if ref and species and ";" not in species:
                pairs[gid].update((part.strip(), species.casefold()) for part in ref.split(";") if part.strip())
    return roster, pairs


def _unscreened_literature(unit: dict, gene_pairs: set[tuple[str, str]],
                          screened: dict[tuple[str, str], set[str]]) -> tuple[list[str], list[str]]:
    """Link candidates are retrieval leads, not absent phenotypes or verified facts.

    Require the originating researcher's screen to account for the cited handle
    before calling a candidate a no-lead screen. Explicitly unrelated reference
    genes need not block another gene; unscoped links remain unresolved.
    """
    pending: set[str] = set()
    for link in unit.get("literature_link_candidates", []):
        if not isinstance(link, dict):
            raise ValueError("invalid frozen literature link")
        reference = str(link.get("reference_gene") or "").strip()
        species = str(link.get("reference_species") or "").strip().casefold()
        if reference and any(pair[0] == reference for pair in gene_pairs):
            if species and (reference, species) not in gene_pairs:
                pending.add("reference_species_conflict:" + reference)
                continue
        elif reference and gene_pairs:
            # This record concerns an explicitly different reference in this
            # mixed unit; it is not a gene-specific paper for this copy.
            continue
        ids = [
            str(link.get(key) or "").strip() for key in ("source_id", "pmid", "doi", "pmcid", "source_url")
        ]
        handle = next((value for value in ids if value), "")
        relevant_pairs = {pair for pair in gene_pairs
                          if (not reference or pair[0] == reference) and (not species or pair[1] == species)}
        if not handle:
            pending.add("unidentified_link")
        elif not relevant_pairs or not all(handle in screened.get(pair, set()) for pair in relevant_pairs):
            pending.add(handle)
    uncertain_status: set[str] = set()
    for status in unit.get("literature_query_status", []):
        if not isinstance(status, dict):
            raise ValueError("invalid frozen literature query status")
        state = next((str(status[key]).strip().casefold()
                      for key in ("retrieval_status", "query_status", "status", "state")
                      if key in status and str(status[key]).strip()), "")
        if state not in {"completed", "success", "retrieved", "no_hit", "no_hits", "no_results"}:
            uncertain_status.add(state or "not_assessed")
    return sorted(pending), sorted(uncertain_status)


def load_followup_plan(path_arg: str, expected_sha256: str, bundle: dict, bundle_sha256: str) -> dict:
    """Consume the coordinator-pinned plan and verify all original ledger anchors."""
    path = Path(path_arg).resolve()
    if not expected_sha256 or digest(path.read_bytes()) != expected_sha256:
        raise ValueError("follow-up plan SHA-256 mismatch")
    plan = _read_json(path)
    if plan.get("contract") != CONTRACT or plan.get("case_id") != bundle["case_id"]:
        raise ValueError("follow-up plan case/contract mismatch")
    inputs = plan.get("input", {})
    for key in ("initial_summary", "screening", "original_bundle", "original_assessments", "initial_gene_ledger"):
        _bound(path.parent, inputs.get(key))
    if inputs["original_bundle"]["sha256"] != bundle_sha256:
        raise ValueError("follow-up plan belongs to a different bundle generation")
    roster, _ = _gene_reference_pairs(bundle)
    rows = plan.get("genes", [])
    actions = {"light_stop_eligible", "retain_upper_evidence", "targeted_followup",
               "needs_light_screen", "needs_grounded_assessment", "not_evaluated"}
    if (not isinstance(rows, list) or len(rows) != len(roster)
            or {r.get("gene_id") for r in rows} != set(roster)
            or any(r.get("unit_id") != roster[r["gene_id"]] or r.get("research_action") not in actions for r in rows)):
        raise ValueError("follow-up plan full gene/action roster mismatch")
    return plan


def plan_light_followup(summary_arg: str, screen_arg: str, out_arg: str) -> dict:
    """Derive an auditable *additional research* plan, never canonical claims.

    All inputs and output locations are bound. A gene can be light-stop eligible
    only on positive, interpreted T4/T5 evidence, an explicit checked upper-Tier
    search with no concrete lead, and no unread relevant input or unresolved
    original observation. Unknown/not-assessed is never a negative screen.
    """
    summary_path = Path(summary_arg).resolve()
    out = Path(out_arg).resolve()
    if out.exists():
        raise ValueError("follow-up output already exists")
    summary_raw = summary_path.read_bytes()
    summary = json.loads(summary_raw)
    if not isinstance(summary, dict) or summary.get("contract") != INITIAL_CONTRACT:
        raise ValueError("expected a sealed Route-B initial summary")
    case = _text(summary.get("case_id"), "case_id")
    assessments_path = _bound(summary_path.parent, summary.get("input_assessments"))
    screening_path = Path(screen_arg).resolve() if screen_arg else assessments_path
    screen_raw = screening_path.read_bytes()
    screening = json.loads(screen_raw)
    if not screen_arg:
        embedded = screening.get("research_screening")
        if not isinstance(embedded, dict) or set(embedded) - {"reference_screens", "gene_assessments"}:
            raise ValueError("initial ledger requires original research_screening records")
        screening = {**embedded, "contract": CONTRACT, "case_id": case,
                     "initial_summary_sha256": digest(summary_raw)}
    if not isinstance(screening, dict) or screening.get("contract") != CONTRACT:
        raise ValueError("invalid light-screen contract")
    if screening.get("case_id") != case or screening.get("initial_summary_sha256") != digest(summary_raw):
        raise ValueError("case or initial summary fingerprint mismatch")
    bundle_path = _bound(summary_path.parent, summary.get("input_bundle"))
    bundle = _read_json(bundle_path)
    if bundle.get("case_id") != case or bundle.get("stage") != "post_b0_semantic_bundle" or bundle.get("route") != "interval_first":
        raise ValueError("input bundle mismatch")
    roster, reference_pairs = _gene_reference_pairs(bundle)
    units_by_id = {unit["unit_id"]: unit for unit in bundle["units"]}
    if len(units_by_id) != len(bundle["units"]):
        raise ValueError("duplicate frozen unit")

    ledger_descriptor = summary.get("genes")
    if not isinstance(ledger_descriptor, dict):
        raise ValueError("initial gene ledger missing")
    ledger_path = _bound(summary_path.parent, {"path": ledger_descriptor.get("file"),
                                               "sha256": ledger_descriptor.get("sha256")})
    original_rows = [json.loads(line) for line in ledger_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    by_gene = {row["gene_id"]: row for row in original_rows}
    if (len(by_gene) != len(original_rows) or set(by_gene) != set(roster)
            or summary.get("gene_count") != len(roster)
            or any(row["unit_id"] != roster[gid] for gid, row in by_gene.items())):
        raise ValueError("initial gene ledger roster mismatch")

    screens: dict[str, dict] = {}
    for screen in screening.get("reference_screens", []):
        if not isinstance(screen, dict):
            raise ValueError("screen record must be an object")
        sid = _text(screen.get("screen_id"), "screen_id")
        if sid in screens:
            raise ValueError("duplicate screen_id")
        ref, species = _text(screen.get("reference_gene"), "reference_gene"), _text(screen.get("reference_species"), "reference_species")
        members = _strings(screen.get("member_gene_ids"), "screen member_gene_ids", required=True)
        status = screen.get("state")
        if status not in SCREEN_STATES:
            raise ValueError("invalid screen state")
        queries = _strings(screen.get("queries", []), "screen queries")
        sources = _strings(screen.get("screened_source_ids", []), "screen source_ids")
        scopes = _strings(screen.get("checked_species_scopes", []), "checked species scopes")
        if set(scopes) - SPECIES_SCOPES:
            raise ValueError("invalid checked species scope")
        reason = _text(screen.get("reason"), "screen reason")
        if status == "no_concrete_upper_lead":
            if not queries or screen.get("retrieval_state") != "completed":
                raise ValueError("no-lead screen requires completed retrieval and query scope")
        if status == "concrete_upper_lead" and not str(screen.get("upper_lead") or "").strip():
            raise ValueError("concrete lead requires an explicit description")
        for gid in members:
            if gid not in roster or (ref, species.casefold()) not in reference_pairs[gid]:
                raise ValueError("screen may not borrow another gene/reference/species")
        screens[sid] = {**screen, "screen_id": sid, "member_gene_ids": members,
                        "reference_gene": ref, "reference_species": species,
                        "queries": queries, "screened_source_ids": sources,
                        "checked_species_scopes": scopes,
                        "reason": reason}

    decisions: dict[str, dict] = {}
    authored_decisions = list(screening.get("gene_assessments", []))
    explicit = {row.get("gene_id") for row in authored_decisions if isinstance(row, dict)}
    conflicting = set()
    original_upper: dict[str, set[str]] = {}
    for gid, original in by_gene.items():
        candidates = {}
        for record in original["original_assessments"]:
            planning = record.get("research_planning")
            if planning is None:
                continue
            if not isinstance(planning, dict):
                raise ValueError("invalid original research planning")
            decision = {"gene_id": gid, "unit_id": roster[gid],
                        "provisional_tier": planning.get("provisional_tier"),
                        "tier_reason": planning.get("tier_reason", record["reason"]),
                        "screen_ids": planning.get("screen_ids", [])}
            if decision["provisional_tier"] not in TIERS:
                raise ValueError("invalid original provisional Tier")
            if decision["provisional_tier"] in UPPER_TIERS:
                original_upper.setdefault(gid, set()).add(decision["provisional_tier"])
            candidates[encoded(decision)] = decision
        if gid not in explicit:
            if len(candidates) == 1:
                authored_decisions.extend(candidates.values())
            elif candidates:
                conflicting.add(gid)
    for decision in authored_decisions:
        if not isinstance(decision, dict):
            raise ValueError("gene assessment must be an object")
        gid = _text(decision.get("gene_id"), "gene assessment gene_id")
        if gid not in roster or decision.get("unit_id") != roster[gid] or gid in decisions:
            raise ValueError("unknown, misplaced or duplicate gene assessment")
        tier = decision.get("provisional_tier")
        if tier not in TIERS:
            raise ValueError("invalid provisional Tier")
        reason = _text(decision.get("tier_reason"), "provisional Tier reason")
        ids = _strings(decision.get("screen_ids", []), "gene screen_ids")
        if any(sid not in screens or gid not in screens[sid]["member_gene_ids"] for sid in ids):
            raise ValueError("gene assessment screen binding mismatch")
        decisions[gid] = {"tier": tier, "reason": reason, "screen_ids": ids}

    result_rows = []
    for gid, uid in roster.items():
        original = by_gene[gid]
        declared = decisions.get(gid)
        connected = [screens[sid] for sid in declared["screen_ids"]] if declared else []
        bound = [screen for screen in screens.values() if gid in screen["member_gene_ids"]]
        unread = {pointer: [f for f in fields if f in UPPER_RELEVANT_FIELDS]
                  for pointer, fields in original.get("fields_without_read_receipt", {}).items()}
        unread = {pointer: fields for pointer, fields in unread.items() if fields}
        checked_scopes: dict[tuple[str, str], set[str]] = {}
        screened_sources: dict[tuple[str, str], set[str]] = {}
        for screen in connected:
            if screen["state"] == "no_concrete_upper_lead":
                pair = (screen["reference_gene"], screen["reference_species"].casefold())
                screened_sources.setdefault(pair, set()).update(screen["screened_source_ids"])
                checked_scopes.setdefault(pair, set()).update(screen["checked_species_scopes"])
        unscreened_refs = reference_pairs[gid] - {
            pair for pair, scopes in checked_scopes.items() if scopes == SPECIES_SCOPES}
        unscreened_literature, unresolved_queries = _unscreened_literature(
            units_by_id[uid], reference_pairs[gid], screened_sources)
        support_fields = (T4_SUPPORT_FIELDS if declared and declared["tier"] == "T4"
                          else T5_SUPPORT_FIELDS)
        grounded = any(
            record.get("state") == "assessed_claim" and record.get("used_material_pointers")
            and set(record.get("assessed_fields", [])) & support_fields
            for record in original.get("original_assessments", [])
        )
        # No decision or an unreported lighter screen cannot authorize a stop.
        if original_upper.get(gid):
            action, why = "retain_upper_evidence", "original_upper_evidence_present"
        elif gid in conflicting:
            action, why = "needs_grounded_assessment", "conflicting_original_planning"
        elif not declared:
            action, why = "not_evaluated", "no_gene_triage"
        elif declared["tier"] in UPPER_TIERS:
            action, why = "retain_upper_evidence", "provisional_upper_evidence_present"
        elif any(s["state"] == "concrete_upper_lead" or str(s.get("upper_lead") or "").strip()
                 for s in bound) or original["has_concrete_lead"]:
            action, why = "targeted_followup", "specific_upper_tier_or_mechanistic_lead"
        elif any(s["state"] == "retrieval_failed" for s in bound):
            action, why = "not_evaluated", "upper_literature_retrieval_failed"
        elif (not connected or any(s["state"] != "no_concrete_upper_lead" for s in bound)
              or unscreened_refs or unscreened_literature or unresolved_queries):
            action, why = "needs_light_screen", "upper_literature_not_fully_screened"
        elif declared["tier"] not in LOWER_TIERS:
            action, why = "not_evaluated", "no_grounded_T4_T5_stop_basis"
        elif not grounded or original.get("judgments_without_material_pointer"):
            action, why = "needs_grounded_assessment", "no_relevant_source_bound_T4_T5_claim"
        elif unread or original["has_unresolved_judgment"]:
            action, why = "targeted_followup", "unread_relevant_or_unresolved_material"
        else:
            action, why = "light_stop_eligible", "grounded_T4_T5_and_no_concrete_upper_lead"
        result_rows.append({
            "gene_id": gid, "unit_id": uid,
            "provisional_tier": declared["tier"] if declared else "not_assessed",
            "provisional_tier_reason": declared["reason"] if declared else "",
            "research_action": action, "reason": why,
            "screen_ids": declared["screen_ids"] if declared else [],
            "unlinked_screen_ids": [s["screen_id"] for s in bound if s not in connected],
            "original_upper_tiers": sorted(original_upper.get(gid, set())),
            "has_original_claim": original["has_original_claim"],
            "has_concrete_lead": original["has_concrete_lead"],
            "unread_relevant_material": unread,
            "unscreened_reference_pairs": [list(pair) for pair in sorted(unscreened_refs)],
            "unscreened_literature_handles": unscreened_literature,
            "unresolved_literature_query_states": unresolved_queries,
            "has_unresolved_judgment": original["has_unresolved_judgment"],
            "original_interpretation_state": original["interpretation_state"],
        })

    grouped: dict[str, list[dict]] = {}
    for row in result_rows:
        grouped.setdefault(row["unit_id"], []).append(row)
    units = []
    for uid, members in grouped.items():
        actions = {row["research_action"] for row in members}
        action = ("light_stop_eligible" if actions == {"light_stop_eligible"}
                  else "targeted_followup" if "targeted_followup" in actions
                  else "retain_upper_evidence" if "retain_upper_evidence" in actions
                  else "requires_assessment")
        units.append({"unit_id": uid, "member_gene_ids": [x["gene_id"] for x in members],
                      "research_action": action, "member_actions": {x["gene_id"]: x["research_action"] for x in members}})
    report = {
        "contract": CONTRACT, "case_id": case,
        "input": {
            "initial_summary": {"path": str(summary_path), "sha256": digest(summary_raw)},
            "screening": {"path": str(screening_path), "sha256": digest(screen_raw)},
            "original_bundle": summary["input_bundle"], "original_assessments": summary["input_assessments"],
            "initial_gene_ledger": {"path": str(ledger_path), "sha256": digest(ledger_path.read_bytes())},
        },
        "reference_screens": list(screens.values()), "genes": result_rows, "units": units,
        "gene_actions": dict(Counter(x["research_action"] for x in result_rows)),
        "unit_actions": dict(Counter(x["research_action"] for x in units)),
        "limitations": (
            "Scheduling advice for ADDITIONAL literature deep reading only. "
            "No gene is excluded, no T1-T6/U is assigned, and no canonical B prompt, "
            "claim, source, support or Wet selection is modified. A checked no-lead "
            "screen is not proof that an upper-Tier paper does not exist. "
            "Study-first model/latency/token impact has not been evaluated."
        ),
    }
    out.mkdir(parents=True)
    target = out / "followup_plan.json"
    target.write_bytes(encoded(report))
    return {"contract": CONTRACT, "case_id": case,
            "plan": {"path": str(target), "sha256": digest(target.read_bytes())},
            "gene_count": len(result_rows), "unit_count": len(units),
            "gene_actions": report["gene_actions"], "unit_actions": report["unit_actions"]}
