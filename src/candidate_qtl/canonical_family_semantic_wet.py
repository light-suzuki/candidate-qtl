"""Canonical family-semantic ordering wrapper (#505/#514).

Runs only after the supported-union-derived deterministic family/copy stage.
LLM/agent output is proposal-only; machine verification and competitor scan are
required before ordering changes. Membership is never changed.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from .canonical_route_b_wet import run_canonical_route_b_wet_extension
from .core import runtime_path
from .family_semantic_competitor_gate import enforce_semantic_competitor_gate
from .family_semantic_resolution import (
    build_family_semantic_request,
    import_family_semantic_responses,
    verify_family_semantic_proposals,
)
from .wet_candidate_handoff import wet_candidate_handoff

SCHEMA_VERSION = "1.0"
FAMILY_SEMANTIC_ROOT = Path("b04_canonical_semantic") / "10_family_semantic"
RESPONSES_NAME = "llm_candidate_responses.json"
LLM_MANIFEST_NAME = "llm_run_manifest.json"


def _existing(path: Path) -> str | None:
    return str(path) if path.exists() else None


def _strict_json_bool(value: object, *, default: bool) -> bool:
    """Accept only a JSON boolean for provenance gates; fail closed otherwise."""
    return value if isinstance(value, bool) else default


def _strict_json_string(value: object) -> str:
    """Do not coerce malformed producer/source provenance into valid labels."""
    return value if isinstance(value, str) else ""


def _membership_ids(path_arg: str | None) -> set[str]:
    if not path_arg:
        return set()
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    for column in ("candidate_gene_id", "gene_id"):
        if column in frame.columns:
            return {v.strip() for v in frame[column].astype(str) if v.strip()}
    raise ValueError("membership artifact requires candidate_gene_id or gene_id")


def _family_semantic_fallback(
    root: Path,
    *,
    atlas_arg: str,
    bridge_arg: str | None,
    family_copy_arg: str | None,
    family_copy_shortlist_arg: str | None,
    case_id: str,
) -> tuple[str | None, str | None, dict[str, object]]:
    semantic_root = root / FAMILY_SEMANTIC_ROOT
    if not family_copy_arg:
        return family_copy_arg, family_copy_shortlist_arg, {
            "state": "not_evaluated_no_family_copy_artifact",
            "rank_contribution": "none", "membership_effect": "none",
        }

    request = build_family_semantic_request(
        family_copy_arg, atlas_arg, str(semantic_root / "01_request"),
        bridge_arg=bridge_arg, case_id=case_id,
    )
    if int(request.get("prompt_count", 0) or 0) == 0:
        return family_copy_arg, family_copy_shortlist_arg, {
            "state": "not_required_deterministic_or_not_evaluable",
            "request": request.get("request"), "prompts": request.get("prompts"),
            "rank_contribution": "none", "membership_effect": "none",
        }

    run_dir = semantic_root / "02_llm_run"
    response = run_dir / RESPONSES_NAME
    llm_manifest = run_dir / LLM_MANIFEST_NAME
    if not response.exists():
        return family_copy_arg, family_copy_shortlist_arg, {
            "state": "semantic_review_available_no_response",
            "request": request.get("request"), "prompts": request.get("prompts"),
            "expected_semantic_run_dir": str(run_dir), "response": None,
            "rank_contribution": "none", "membership_effect": "none",
        }

    imported = import_family_semantic_responses(
        str(request["request"]), str(request["prompts"]), str(response),
        str(semantic_root / "03_import"),
        llm_manifest_arg=str(llm_manifest) if llm_manifest.exists() else None,
        case_id=case_id,
    )
    preliminary = verify_family_semantic_proposals(
        str(imported["proposals"]), family_copy_arg, atlas_arg,
        str(semantic_root / "04_preliminary_verify"),
        base_shortlist_arg=family_copy_shortlist_arg, case_id=case_id,
    )
    final = enforce_semantic_competitor_gate(
        str(preliminary["family_copy"]), family_copy_arg, atlas_arg,
        str(semantic_root / "05_competitor_gate"),
        base_shortlist_arg=family_copy_shortlist_arg, case_id=case_id,
    )
    promoted = int(final.get("promoted_family_count", 0) or 0)
    return str(final["family_copy"]), str(final["shortlist"]), {
        "state": "verified_semantic_copy_resolution" if promoted else "evaluated_no_verified_copy_mapping",
        "request": request.get("request"), "prompts": request.get("prompts"),
        "response": str(response), "proposal_manifest": imported.get("manifest"),
        "preliminary_verification_manifest": preliminary.get("manifest"),
        "competitor_gate_manifest": final.get("manifest"),
        "promoted_family_count": promoted,
        "promoted_copy_count": int(final.get("promoted_copy_count", 0) or 0),
        "competitor_expansion_count": int(final.get("competitor_expansion_count", 0) or 0),
        "rank_contribution": "none", "membership_effect": "none",
    }


def _write_manifest(root: Path, payload: dict[str, object]) -> str:
    out = root / FAMILY_SEMANTIC_ROOT
    out.mkdir(parents=True, exist_ok=True)
    path = out / "canonical_family_semantic_wet_manifest.json"
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return str(path)


def run_canonical_route_b_family_semantic_wet_extension(
    profile_out_dir_arg: str,
    *, phenotype: str,
    responses_arg: str | None = None,
    llm_manifest_arg: str | None = None,
    source_cache_arg: str | None = None,
    bound_claims_arg: str | None = None,
    agent_claim_review_arg: str | None = None,
    source_review_records_arg: str | None = None,
    record_cache_arg: str | None = None,
    allow_external_api: bool = False,
    wet_budget: int = 20,
    deterministic_subfamily_refinement: bool = False,
    route_c_members_arg: str | None = None,
    copy_evidence_args: dict[str, str | None] | None = None,
    case_id: str = "",
) -> dict[str, object]:
    root = runtime_path(profile_out_dir_arg)
    base = run_canonical_route_b_wet_extension(
        profile_out_dir_arg, phenotype=phenotype, responses_arg=responses_arg,
        llm_manifest_arg=llm_manifest_arg, source_cache_arg=source_cache_arg,
        bound_claims_arg=bound_claims_arg,
        agent_claim_review_arg=agent_claim_review_arg,
        source_review_records_arg=source_review_records_arg,
        record_cache_arg=record_cache_arg,
        allow_external_api=allow_external_api, wet_budget=wet_budget,
        deterministic_subfamily_refinement=deterministic_subfamily_refinement, case_id=case_id,
        route_c_members_arg=route_c_members_arg,
        **({"copy_evidence_args": copy_evidence_args} if copy_evidence_args is not None else {}),
    )
    base_manifest_path = runtime_path(str(base["manifest"]))
    base_manifest = json.loads(base_manifest_path.read_text(encoding="utf-8"))
    copy_provenance = ({"copy_evidence": base_manifest["copy_evidence"]}
                       if "copy_evidence" in base_manifest else {})

    if not base.get("wet_candidates"):
        semantic = {"state": "not_run_without_frozen_wet_membership", "rank_contribution": "none", "membership_effect": "none"}
        wrapper = {**copy_provenance, "schema_version": SCHEMA_VERSION, "case_id": case_id, "state": base.get("state"),
                   "base_canonical_manifest": str(base_manifest_path), "family_semantic_resolution": semantic,
                   "final_wet_candidates": base.get("wet_candidates"), "final_wet_handoff_manifest": base.get("wet_handoff_manifest"),
                   "legacy_family_context_consumed": False, "rank_contribution": "none"}
        path = _write_manifest(root, wrapper)
        return {**base, "manifest": path, "base_manifest": str(base_manifest_path), "family_semantic_resolution": semantic}

    producer = _strict_json_string(base_manifest.get("family_copy_producer_class", ""))
    source = _strict_json_string(base_manifest.get("family_recovery_source", ""))
    legacy = _strict_json_bool(
        base_manifest.get("legacy_triple_route_family_context_consumed", True),
        default=True,
    )
    if producer != "canonical_supported_family_copy" or source != "supported_member_union" or legacy:
        semantic = {"state": "not_evaluated_noncanonical_family_copy_provenance",
                    "producer_class": producer, "family_recovery_source": source,
                    "legacy_triple_route_family_context_consumed": legacy,
                    "rank_contribution": "none", "membership_effect": "none"}
        wrapper = {**copy_provenance, "schema_version": SCHEMA_VERSION, "case_id": case_id, "state": base.get("state"),
                   "base_canonical_manifest": str(base_manifest_path), "family_semantic_resolution": semantic,
                   "final_wet_candidates": base.get("wet_candidates"), "final_wet_handoff_manifest": base.get("wet_handoff_manifest"),
                   "legacy_family_context_consumed": legacy, "rank_contribution": "none"}
        path = _write_manifest(root, wrapper)
        return {**base, "manifest": path, "base_manifest": str(base_manifest_path), "family_semantic_resolution": semantic}

    atlas = root / "b01_route_b" / "interval_gene_atlas.tsv"
    deterministic_copy = str(base_manifest.get("family_copy") or "") or None
    deterministic_shortlist = str(base_manifest.get("family_copy_shortlist") or "") or None
    canonical_bridge = str(base_manifest.get("canonical_family_bridge") or "") or None
    effective_copy, effective_shortlist, semantic = _family_semantic_fallback(
        root, atlas_arg=str(atlas), bridge_arg=canonical_bridge,
        family_copy_arg=deterministic_copy, family_copy_shortlist_arg=deterministic_shortlist,
        case_id=case_id,
    )

    if effective_copy == deterministic_copy and effective_shortlist == deterministic_shortlist:
        wrapper = {**copy_provenance, "schema_version": SCHEMA_VERSION, "case_id": case_id, "state": base.get("state"),
                   "base_canonical_manifest": str(base_manifest_path), "family_semantic_resolution": semantic,
                   "effective_family_copy": effective_copy, "effective_family_copy_shortlist": effective_shortlist,
                   "final_wet_candidates": base.get("wet_candidates"), "final_wet_handoff_manifest": base.get("wet_handoff_manifest"),
                   "membership_identity_verified": True, "legacy_family_context_consumed": False, "rank_contribution": "none"}
        path = _write_manifest(root, wrapper)
        return {**base, "manifest": path, "base_manifest": str(base_manifest_path),
                "family_semantic_resolution": semantic, "effective_family_copy": effective_copy,
                "effective_family_copy_shortlist": effective_shortlist}

    route_a_bridge = _existing(root / "z01b_family_candidate_bridge" / "family_candidate_bridge.tsv")
    ranked = _existing(root / "z03_case_run" / "rank" / "ranked_candidates.tsv")
    route_b = str(base_manifest.get("route_b_verified_candidates") or "") or None
    recovery = str(base_manifest.get("canonical_family_recovery") or "") or None
    convergence = str(base_manifest.get("canonical_family_comparison") or "") or None
    if not route_b or not recovery or not convergence:
        raise RuntimeError("canonical semantic layer is missing #517 provenance artifacts")

    final_handoff = wet_candidate_handoff(
        str(root / "z05_canonical_family_semantic_wet_handoff"),
        **base_manifest.get("route_review_inputs", {}),
        route_a_bridge_arg=route_a_bridge, route_b_candidates_arg=route_b,
        route_c_members_arg=route_c_members_arg,
        route_a_manifest_arg=_existing(root / "z01b_family_candidate_bridge" / "family_candidate_bridge_manifest.json"),
        route_b_candidates_manifest_arg=_existing(Path(route_b).with_name("route_b_verified_candidates_manifest.json")),
        route_c_manifest_arg=_existing(Path(route_c_members_arg).with_name("route_c_direct_members_manifest.json")) if route_c_members_arg else None,
        ranked_arg=ranked, family_recovery_arg=recovery, convergence_arg=convergence,
        family_copy_arg=effective_copy, family_copy_shortlist_arg=effective_shortlist,
        wet_budget=wet_budget,
        interval_atlas_arg=str(atlas),
        reference_function_arg=_existing(root / "b01_route_b" / "route_b_reference_function.tsv"),
        deterministic_subfamily_refinement=deterministic_subfamily_refinement,
        case_id=case_id,
    )
    before = _membership_ids(str(base.get("wet_candidates") or ""))
    after = _membership_ids(str(final_handoff.get("wet_candidates") or ""))
    if before != after:
        raise RuntimeError("semantic family/copy ordering changed canonical Wet membership")

    final_manifest = json.loads(runtime_path(str(final_handoff["manifest"])).read_text(encoding="utf-8"))
    expected_sha = str(base_manifest.get("handoff_supported_union_sha256", ""))
    final_sha = str(final_manifest.get("supported_member_union_sha256", ""))
    if not expected_sha or expected_sha != final_sha:
        raise RuntimeError("semantic ordering rebuilt a different supported-member union")

    wrapper = {**copy_provenance, "schema_version": SCHEMA_VERSION, "case_id": case_id, "state": final_handoff.get("handoff_state"),
               "base_canonical_manifest": str(base_manifest_path), "family_semantic_resolution": semantic,
               "effective_family_copy": effective_copy, "effective_family_copy_shortlist": effective_shortlist,
               "final_wet_candidates": final_handoff.get("wet_candidates"), "final_wet_handoff_manifest": final_handoff.get("manifest"),
               "membership_identity_verified": True, "supported_union_identity_verified": True,
               "supported_union_sha256": final_sha, "legacy_family_context_consumed": False, "rank_contribution": "none"}
    path = _write_manifest(root, wrapper)
    return {**base, "state": wrapper["state"], "manifest": path, "base_manifest": str(base_manifest_path),
            "family_semantic_resolution": semantic, "effective_family_copy": effective_copy,
            "effective_family_copy_shortlist": effective_shortlist, "wet_handoff_manifest": final_handoff.get("manifest"),
            "wet_candidates": final_handoff.get("wet_candidates"), "stage_a_narrowed": final_handoff.get("stage_a_narrowed"),
            "membership_state": final_handoff.get("membership_state"), "input_coverage_state": final_handoff.get("input_coverage_state")}
