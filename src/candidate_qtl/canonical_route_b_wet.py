"""Canonical post-B0 Route-B semantic -> family/copy -> Wet extension.

The historical ``triple-route`` profile freezes Route A, Route B0 and Route C
artifacts, but its phenotype-aware reverse-lookup B1 comparison is diagnostic
only.  Production continuation is rebuilt from the frozen B0 atlas plus final
source+transfer-verified Route-B members:

    frozen B0 atlas
      -> semantic bundle/prompts
      -> agent/LLM responses
      -> proposal import
      -> source relation verification
      -> transfer verification
      -> frozen Route-B verified candidates
      -> supported A/B member union
      -> canonical supported-family recovery (#514)
      -> existing deterministic family-to-copy
      -> canonical Wet handoff / ordering

The first executed profile run emits prompts and stops at
``semantic_review_required``. Run the existing ``llm-run`` command with
``--out-dir <profile>/b04_canonical_semantic/02_llm_run`` (or place an
agent-authored artifact with the same contract there), then re-run the same
profile command. The extension discovers that artifact by convention, so no new
CLI surface is needed. Missing upstream/provider evidence is explicit
``not_runnable``/``not_evaluated``, never biological absence.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from .canonical_family_recovery import build_canonical_family_recovery
from .core import runtime_path
from .copy_evidence_inputs import bind_copy_evidence, verify_copy_evidence
from .family_to_copy import family_to_copy
from .route_b_semantic_bundle import build_route_b_semantic_bundle
from .route_b_semantic_candidates import freeze_route_b_semantic_candidates
from .route_b_semantic_import import (
    import_route_b_semantic_proposals,
    validate_route_b_semantic_artifacts,
)
from .route_b_semantic_source_verify import verify_route_b_semantic_sources
from .route_b_semantic_transfer_verify import verify_route_b_semantic_transfer
from .supported_member_union import build_supported_member_union
from .wet_candidate_handoff import wet_candidate_handoff

SCHEMA_VERSION = "1.1"
SEMANTIC_RUN_DIR = "02_llm_run"
SEMANTIC_RESPONSES_NAME = "llm_candidate_responses.json"
SEMANTIC_LLM_MANIFEST_NAME = "llm_run_manifest.json"


def _existing(path: Path) -> str | None:
    return str(path) if path.exists() else None


def _sha256(path_arg: str) -> str:
    return hashlib.sha256(runtime_path(path_arg).read_bytes()).hexdigest()


def _route_reviews(root: Path, frozen_b: dict[str, object] | None = None) -> dict[str, list[str]]:
    """Pass sealed review artifacts together without changing supported membership."""
    paths: list[str] = []
    manifests: list[str] = []
    a_path = root / "z01b_family_candidate_bridge" / "route_a_review_candidates.tsv"
    if a_path.exists():
        paths.append(str(a_path))
        manifests.append(str(a_path.with_name("route_a_review_candidates_manifest.json")))
    if frozen_b and frozen_b.get("review_candidates"):
        paths.append(str(frozen_b["review_candidates"]))
        manifests.append(str(frozen_b["review_manifest"]))
    return {"review_candidates_arg": paths, "review_candidates_manifest_arg": manifests}


def _write_terminal(
    semantic_root: Path,
    *,
    case_id: str,
    state: str,
    provider_state: str,
    reason: str,
    copy_provenance: dict[str, object] | None = None,
) -> dict[str, object]:
    semantic_root.mkdir(parents=True, exist_ok=True)
    manifest_path = semantic_root / "canonical_route_b_extension_manifest.json"
    manifest = {
        **(copy_provenance or {}),
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id,
        "state": state,
        "canonical_agent_path": "route_b_semantic_bundle_to_verified_candidates",
        "legacy_interval_reverse_lookup_role": "diagnostic_only",
        "legacy_b1_consumed": False,
        "legacy_triple_route_family_context_consumed": False,
        "provider_state": provider_state,
        "reason": reason,
        "route_b_verified_candidates": None,
        "canonical_family_recovery": None,
        "wet_candidates": None,
        "rank_contribution": "none",
        "policy": [
            "legacy interval-reverse-lookup B1 is never a production membership input",
            "legacy triple-route family/copy artifacts are never a production ordering input",
            "missing upstream/provider evidence is not_evaluated/not_runnable, never negative evidence",
        ],
    }
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        "state": state,
        "manifest": str(manifest_path),
        "route_b_verified_candidates": None,
        "wet_candidates": None,
        "legacy_b1_consumed": False,
        "legacy_triple_route_family_context_consumed": False,
    }


def _discover_semantic_response(
    semantic_root: Path,
    responses_arg: str | None,
    llm_manifest_arg: str | None,
) -> tuple[str | None, str | None, str]:
    """Resolve an explicit response or the standard post-prompt artifact path."""
    if responses_arg:
        response_path = runtime_path(responses_arg)
        manifest_path = runtime_path(llm_manifest_arg) if llm_manifest_arg else (
            response_path.parent / SEMANTIC_LLM_MANIFEST_NAME
        )
        return str(response_path), str(manifest_path), "explicit"
    run_dir = semantic_root / SEMANTIC_RUN_DIR
    response = run_dir / SEMANTIC_RESPONSES_NAME
    if not response.exists():
        manifest = run_dir / SEMANTIC_LLM_MANIFEST_NAME
        if manifest.exists() or llm_manifest_arg:
            return None, str(runtime_path(llm_manifest_arg)) if llm_manifest_arg else str(manifest), "manifest_without_response"
        return None, None, "not_available"
    llm_manifest = run_dir / SEMANTIC_LLM_MANIFEST_NAME
    manifest_path = runtime_path(llm_manifest_arg) if llm_manifest_arg else llm_manifest
    return str(response), str(manifest_path), "conventional_path"


def run_canonical_route_b_wet_extension(
    profile_out_dir_arg: str,
    *,
    phenotype: str,
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
    """Continue an executed triple-route profile through verified B -> Wet.

    ``responses_arg`` exists for direct module/tests. Normal guided-profile usage
    discovers ``b04_canonical_semantic/02_llm_run/llm_candidate_responses.json``
    automatically. Provider choice remains outside the deterministic verification
    boundary because this extension never invokes a model itself.
    """
    root = runtime_path(profile_out_dir_arg)
    if not root.exists():
        raise ValueError(f"guided-profile output directory does not exist: {root}")
    if wet_budget < 1:
        raise ValueError("wet_budget must be >= 1")

    copy_kwargs, copy_provenance = bind_copy_evidence(root, copy_evidence_args)
    semantic_root = root / "b04_canonical_semantic"
    if not str(phenotype).strip():
        return _write_terminal(
            semantic_root,
            copy_provenance=copy_provenance,
            case_id=case_id,
            state="not_runnable_missing_phenotype",
            provider_state="not_evaluated_missing_phenotype",
            reason="post-B0 semantic review requires a phenotype; no biological negative was inferred",
        )

    atlas = root / "b01_route_b" / "interval_gene_atlas.tsv"
    if not atlas.exists():
        return _write_terminal(
            semantic_root,
            copy_provenance=copy_provenance,
            case_id=case_id,
            state="not_runnable_upstream_b0_missing",
            provider_state="not_evaluated_upstream_b0_missing",
            reason="frozen Route B0 atlas was not produced by the profile",
        )

    extension_manifest_path = semantic_root / "canonical_route_b_extension_manifest.json"
    responses_arg, llm_manifest_arg, response_discovery = _discover_semantic_response(
        semantic_root, responses_arg, llm_manifest_arg
    )

    responses_path = runtime_path(responses_arg) if responses_arg else None
    llm_manifest_path = runtime_path(llm_manifest_arg) if llm_manifest_arg else None
    if llm_manifest_path and not responses_path:
        raise ValueError(
            "LLM run manifest exists without a semantic response; refusing to regenerate frozen artifacts"
        )
    if responses_path and not responses_path.exists():
        raise ValueError(f"Route B semantic responses do not exist: {responses_path}")

    artifact_hashes: dict[str, str] | None = None
    if responses_path:
        bundle_path = semantic_root / "01_bundle" / "route_b_semantic_bundle.json"
        prompts_path = semantic_root / "01_bundle" / "route_b_semantic_prompts.jsonl"
        if not bundle_path.exists() or not prompts_path.exists():
            raise ValueError(
                "frozen Route B semantic response exists but 01_bundle artifacts are missing; "
                "refusing to regenerate them"
            )
        artifact_hashes = validate_route_b_semantic_artifacts(
            str(bundle_path), str(prompts_path), str(responses_path),
            str(llm_manifest_path) if llm_manifest_path else None,
            require_manifest=True,
        )
        bundle_result = {
            "bundle": str(bundle_path),
            "prompts": str(prompts_path),
            "artifact_hashes": artifact_hashes,
        }
    else:
        literature = _existing(root / "b03_reverse_literature" / "reverse_literature.tsv")
        literature_status = _existing(root / "b03_reverse_literature" / "reverse_literature_status.tsv")
        lookup = _existing(root / "b01_route_b" / "route_b0" / "interval_reverse_lookup.tsv")
        identity_hits = _existing(root / "b01_route_b" / "route_b_identity_hits.tsv")
        bundle_result = build_route_b_semantic_bundle(
            str(atlas),
            str(semantic_root / "01_bundle"),
            phenotype=str(phenotype),
            literature_arg=literature,
            literature_status_arg=literature_status,
            lookup_arg=lookup,
            identity_hits_arg=identity_hits,
            case_id=case_id,
        )

    if not responses_arg and route_c_members_arg:
        # Direct C is an independent optional rescue lane.  A missing Route-B
        # semantic response must stay visible, but must not block a verified
        # same-species direct rescue or cross-species Wet review coverage.
        route_a_bridge = _existing(root / "z01b_family_candidate_bridge" / "family_candidate_bridge.tsv")
        ranked = _existing(root / "z03_case_run" / "rank" / "ranked_candidates.tsv")
        pre_union = build_supported_member_union(
            str(semantic_root / "07_supported_union_for_family"),
            route_a_bridge_arg=route_a_bridge,
            route_c_members_arg=route_c_members_arg,
            route_a_manifest_arg=_existing(root / "z01b_family_candidate_bridge" / "family_candidate_bridge_manifest.json"),
            route_c_manifest_arg=_existing(Path(route_c_members_arg).with_name("route_c_direct_members_manifest.json")) if route_c_members_arg else None,
            case_id=case_id,
        )
        canonical_family = build_canonical_family_recovery(
            str(pre_union["candidate_set"]),
            str(semantic_root / "08_canonical_family_recovery"),
            legacy_bridge_arg=route_a_bridge,
            case_id=case_id,
        )
        canonical_copy_root = semantic_root / "09_canonical_family_to_copy"
        family_copy = family_copy_shortlist = family_copy_manifest = None
        if int(canonical_family.get("family_count", 0) or 0) > 0:
            verify_copy_evidence(copy_kwargs, copy_provenance)
            family_to_copy(
                **copy_kwargs,
                families_arg=str(canonical_family["family_recovery"]),
                out_dir_arg=str(canonical_copy_root), atlas_arg=str(atlas),
                case_id=case_id, bridge_arg=str(canonical_family["family_bridge"]),
            )
            verify_copy_evidence(copy_kwargs, copy_provenance)
            if copy_provenance:
                copy_provenance["copy_evidence"]["consumer_state"] = "passed_to_family_to_copy"
            family_copy = _existing(canonical_copy_root / "family_copy_discrimination.tsv")
            family_copy_shortlist = _existing(canonical_copy_root / "family_copy_shortlist.tsv")
            family_copy_manifest = _existing(canonical_copy_root / "family_copy_manifest.json")
        handoff = wet_candidate_handoff(
            str(root / "z04_canonical_wet_handoff"),
            **_route_reviews(root),
            route_a_bridge_arg=route_a_bridge, route_c_members_arg=route_c_members_arg,
            route_a_manifest_arg=_existing(root / "z01b_family_candidate_bridge" / "family_candidate_bridge_manifest.json"),
            route_c_manifest_arg=_existing(Path(route_c_members_arg).with_name("route_c_direct_members_manifest.json")) if route_c_members_arg else None,
            ranked_arg=ranked, family_recovery_arg=str(canonical_family["family_recovery"]),
            convergence_arg=str(canonical_family["comparison"]), family_copy_arg=family_copy,
            family_copy_shortlist_arg=family_copy_shortlist, wet_budget=wet_budget,
            interval_atlas_arg=str(atlas),
            reference_function_arg=_existing(root / "b01_route_b" / "route_b_reference_function.tsv"),
            deterministic_subfamily_refinement=deterministic_subfamily_refinement, case_id=case_id,
        )
        handoff_payload = json.loads(runtime_path(str(handoff["manifest"])).read_text(encoding="utf-8"))
        pre_union_sha256 = _sha256(str(pre_union["candidate_set"]))
        if pre_union_sha256 != str(handoff_payload.get("supported_member_union_sha256", "")):
            raise RuntimeError("Direct-C family recovery and Wet handoff used different supported-member unions")
        manifest = {
            **copy_provenance,
            "schema_version": SCHEMA_VERSION, "case_id": case_id,
            "state": str(handoff["handoff_state"]),
            "provider_state": "route_b_not_evaluated_no_semantic_response",
            "semantic_response_discovery": response_discovery,
            "prompts": str(bundle_result["prompts"]), "bundle": str(bundle_result["bundle"]),
            "route_b_verified_candidates": None, "route_c_members": str(runtime_path(route_c_members_arg)),
            "route_review_inputs": _route_reviews(root),
            "pre_family_supported_union": str(pre_union["candidate_set"]),
            "pre_family_supported_union_sha256": pre_union_sha256,
            "canonical_family_recovery": str(canonical_family["family_recovery"]),
            "canonical_family_comparison": str(canonical_family["comparison"]),
            "canonical_family_bridge": str(canonical_family["family_bridge"]),
            "family_copy": family_copy, "family_copy_shortlist": family_copy_shortlist,
            "family_copy_manifest": family_copy_manifest,
            "family_copy_producer_class": "canonical_supported_family_copy" if family_copy else "not_evaluated_no_family_rows",
            "family_recovery_source": "supported_member_union",
            "handoff_supported_union_sha256": pre_union_sha256,
            "supported_union_identity_verified": True,
            "wet_handoff_manifest": str(handoff["manifest"]), "wet_candidates": handoff.get("wet_candidates"),
            "membership_state": handoff.get("membership_state"),
            "input_coverage_state": handoff.get("input_coverage_state"),
            "legacy_b1_consumed": False, "legacy_triple_route_family_context_consumed": False,
            "rank_contribution": "none",
            "policy": [
                "Direct C may continue while Route B semantic review is not_evaluated",
                "only verified_direct_support enters supported membership",
                "cross-species Direct C remains tier2 review coverage",
            ],
        }
        semantic_root.mkdir(parents=True, exist_ok=True)
        extension_manifest_path.write_text(
            json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        return {
            "state": manifest["state"], "manifest": str(extension_manifest_path),
            "prompts": str(bundle_result["prompts"]), "route_b_verified_candidates": None,
            "route_c_members": str(runtime_path(route_c_members_arg)),
            "canonical_family_recovery": str(canonical_family["family_recovery"]),
            "family_copy": family_copy, "wet_handoff_manifest": str(handoff["manifest"]),
            "wet_candidates": handoff.get("wet_candidates"), "membership_state": handoff.get("membership_state"),
            "input_coverage_state": handoff.get("input_coverage_state"),
            "legacy_b1_consumed": False, "legacy_triple_route_family_context_consumed": False,
        }

    if not responses_arg:
        manifest = {
            **copy_provenance,
            "schema_version": SCHEMA_VERSION,
            "case_id": case_id,
            "state": "semantic_review_required",
            "canonical_agent_path": "route_b_semantic_bundle_to_verified_candidates",
            "legacy_interval_reverse_lookup_role": "diagnostic_only",
            "legacy_b1_consumed": False,
            "legacy_triple_route_family_context_consumed": False,
            "provider_state": "not_evaluated_no_semantic_response",
            "semantic_response_discovery": response_discovery,
            "expected_semantic_run_dir": str(semantic_root / SEMANTIC_RUN_DIR),
            "prompts": str(bundle_result["prompts"]),
            "bundle": str(bundle_result["bundle"]),
            "route_b_verified_candidates": None,
            "canonical_family_recovery": None,
            "wet_candidates": None,
            "wet_handoff_manifest": None,
            "rank_contribution": "none",
            "policy": [
                "the canonical Route B agent/LLM path starts from the frozen B0 atlas",
                "legacy interval-reverse-lookup B1 output is never consumed by this extension",
                "legacy triple-route family/copy artifacts are not canonical production context",
                "no semantic response is review-required/not-evaluated, never negative evidence",
                "run llm-run into the conventional 02_llm_run directory or place an equivalent agent artifact there",
                "a response must be accompanied by a hash-bearing llm_run_manifest before import",
                "only source+transfer-verified semantic claims may reach production membership",
            ],
        }
        semantic_root.mkdir(parents=True, exist_ok=True)
        extension_manifest_path.write_text(
            json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return {
            "state": manifest["state"],
            "manifest": str(extension_manifest_path),
            "prompts": str(bundle_result["prompts"]),
            "bundle": str(bundle_result["bundle"]),
            "expected_semantic_run_dir": manifest["expected_semantic_run_dir"],
            "route_b_verified_candidates": None,
            "wet_candidates": None,
            "legacy_b1_consumed": False,
            "legacy_triple_route_family_context_consumed": False,
        }

    imported = import_route_b_semantic_proposals(
        str(bundle_result["bundle"]),
        str(bundle_result["prompts"]),
        str(responses_path),
        str(semantic_root / "03_import"),
        llm_manifest_arg=llm_manifest_arg,
        case_id=case_id,
    )
    source_verified = verify_route_b_semantic_sources(
        str(imported["proposals"]),
        str(semantic_root / "04_source_verify"),
        bundle_arg=str(bundle_result["bundle"]),
        cache_arg=source_cache_arg,
        allow_external_api=allow_external_api,
        case_id=case_id,
    )
    transfer_verified = verify_route_b_semantic_transfer(
        str(imported["proposals"]),
        str(source_verified["source_verification"]),
        str(bundle_result["bundle"]),
        str(semantic_root / "05_transfer_verify"),
        bound_claims_arg=bound_claims_arg,
        agent_claim_review_arg=agent_claim_review_arg,
        source_review_records_arg=source_review_records_arg,
        record_cache_arg=record_cache_arg,
        case_id=case_id,
    )
    frozen_b = freeze_route_b_semantic_candidates(
        str(transfer_verified["claims"]),
        str(semantic_root / "06_candidate_freeze"),
        case_id=case_id,
    )

    # Production membership is frozen by the same A/B artifacts later supplied
    # to wet_candidate_handoff.  Build it once here so family recovery can be a
    # pure projection of the canonical supported set.  The handoff intentionally
    # rebuilds the union with the same inputs; SHA equality below detects any
    # future drift between the two calls fail-closed without expanding the
    # wet_candidate_handoff API.
    route_a_bridge = _existing(root / "z01b_family_candidate_bridge" / "family_candidate_bridge.tsv")
    ranked = _existing(root / "z03_case_run" / "rank" / "ranked_candidates.tsv")
    pre_union = build_supported_member_union(
        str(semantic_root / "07_supported_union_for_family"),
        route_a_bridge_arg=route_a_bridge,
        route_b_candidates_arg=str(frozen_b["candidates"]),
        route_c_members_arg=route_c_members_arg,
        route_a_manifest_arg=_existing(root / "z01b_family_candidate_bridge" / "family_candidate_bridge_manifest.json"),
        route_b_candidates_manifest_arg=str(frozen_b["manifest"]),
        route_c_manifest_arg=_existing(Path(route_c_members_arg).with_name("route_c_direct_members_manifest.json")) if route_c_members_arg else None,
        case_id=case_id,
    )

    canonical_family = build_canonical_family_recovery(
        str(pre_union["candidate_set"]),
        str(semantic_root / "08_canonical_family_recovery"),
        legacy_bridge_arg=route_a_bridge,
        case_id=case_id,
    )

    canonical_copy_root = semantic_root / "09_canonical_family_to_copy"
    family_copy: str | None = None
    family_copy_shortlist: str | None = None
    family_copy_manifest: str | None = None
    if int(canonical_family.get("family_count", 0) or 0) > 0:
        verify_copy_evidence(copy_kwargs, copy_provenance)
        family_to_copy(
            **copy_kwargs,
            families_arg=str(canonical_family["family_recovery"]),
            out_dir_arg=str(canonical_copy_root),
            atlas_arg=str(atlas),
            case_id=case_id,
            bridge_arg=str(canonical_family["family_bridge"]),
        )
        verify_copy_evidence(copy_kwargs, copy_provenance)
        if copy_provenance:
            copy_provenance["copy_evidence"]["consumer_state"] = "passed_to_family_to_copy"
        family_copy = _existing(canonical_copy_root / "family_copy_discrimination.tsv")
        family_copy_shortlist = _existing(canonical_copy_root / "family_copy_shortlist.tsv")
        family_copy_manifest = _existing(canonical_copy_root / "family_copy_manifest.json")

    handoff = wet_candidate_handoff(
        str(root / "z04_canonical_wet_handoff"),
        **_route_reviews(root, frozen_b),
        route_a_bridge_arg=route_a_bridge,
        route_b_candidates_arg=str(frozen_b["candidates"]),
        route_c_members_arg=route_c_members_arg,
        ranked_arg=ranked,
        family_recovery_arg=str(canonical_family["family_recovery"]),
        convergence_arg=str(canonical_family["comparison"]),
        family_copy_arg=family_copy,
        family_copy_shortlist_arg=family_copy_shortlist,
        wet_budget=wet_budget,
        # The same frozen B0 inventory the family layer used. Under budget,
        # its unsupported genes fill the remaining review slots in tier
        # order rather than leaving the bench with a handful.
        interval_atlas_arg=str(atlas) if atlas else None,
        reference_function_arg=_existing(root / "b01_route_b" / "route_b_reference_function.tsv"),
        deterministic_subfamily_refinement=deterministic_subfamily_refinement,
        case_id=case_id,
    )

    # The family layer and Wet membership must see the same frozen supported
    # universe. Any mismatch means one producer changed semantics and the run is
    # invalid rather than something to reconcile heuristically.
    handoff_manifest_payload = json.loads(
        runtime_path(str(handoff["manifest"])).read_text(encoding="utf-8")
    )
    pre_union_sha256 = _sha256(str(pre_union["candidate_set"]))
    handoff_union_sha256 = str(handoff_manifest_payload.get("supported_member_union_sha256", ""))
    if not handoff_union_sha256 or pre_union_sha256 != handoff_union_sha256:
        raise RuntimeError(
            "canonical family recovery and Wet handoff used different supported-member unions"
        )

    manifest = {
        **copy_provenance,
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id,
        "state": str(handoff["handoff_state"]),
        "canonical_agent_path": "route_b_semantic_bundle_to_verified_candidates",
        "legacy_interval_reverse_lookup_role": "diagnostic_only",
        "legacy_b1_consumed": False,
        "legacy_triple_route_family_context_consumed": False,
        "provider_state": (
            "verified_semantic_candidates_available"
            if int(frozen_b.get("candidate_gene_count", 0) or 0) > 0
            else "evaluated_no_verified_semantic_candidates"
        ),
        "semantic_response_discovery": response_discovery,
        "semantic_response": str(responses_path),
        "llm_manifest": str(llm_manifest_path) if llm_manifest_path else None,
        "semantic_artifact_hashes": artifact_hashes,
        "bundle": str(bundle_result["bundle"]),
        "prompts": str(bundle_result["prompts"]),
        "proposal_manifest": str(imported["manifest"]),
        "source_verification_manifest": str(source_verified["manifest"]),
        "transfer_verification_manifest": str(transfer_verified["manifest"]),
        "route_b_verified_candidates": str(frozen_b["candidates"]),
        "route_b_verified_candidates_manifest": str(frozen_b["manifest"]),
        "route_b_review_candidates": frozen_b.get("review_candidates"),
        "route_b_review_candidates_manifest": frozen_b.get("review_manifest"),
        "route_review_inputs": _route_reviews(root, frozen_b),
        "pre_family_supported_union": str(pre_union["candidate_set"]),
        "pre_family_supported_union_sha256": pre_union_sha256,
        "canonical_family_recovery": str(canonical_family["family_recovery"]),
        "canonical_family_comparison": str(canonical_family["comparison"]),
        "canonical_family_bridge": str(canonical_family["family_bridge"]),
        "canonical_family_recovery_manifest": str(canonical_family["manifest"]),
        "family_copy": family_copy,
        "family_copy_shortlist": family_copy_shortlist,
        "family_copy_manifest": family_copy_manifest,
        "family_copy_producer_class": "canonical_supported_family_copy" if family_copy else "not_evaluated_no_family_rows",
        "family_recovery_source": "supported_member_union",
        "handoff_supported_union_sha256": handoff_union_sha256,
        "supported_union_identity_verified": True,
        "wet_handoff_manifest": str(handoff["manifest"]),
        "wet_candidates": handoff.get("wet_candidates"),
        "membership_state": handoff.get("membership_state"),
        "input_coverage_state": handoff.get("input_coverage_state"),
        "rank_contribution": "none",
        "policy": [
            "the canonical Route B semantic/agent path starts from the frozen B0 atlas",
            "only final source+transfer-verified Route B candidates enter supported_member_union",
            "canonical family recovery is projected only from supported_member_union",
            "legacy interval-reverse-lookup/triple-route family artifacts are diagnostic only and are not consumed",
            "family-to-copy reuses the existing deterministic implementation; no new copy score is introduced",
            "Route-B semantic verification is not translated into historical bridge mapping tiers",
            "pre-family and Wet-handoff supported unions must be byte-identical or the run fails closed",
            "the frozen bundle/prompt/response bytes must match the hash-bearing LLM manifest",
            "provider failure/not-evaluated remains explicit and never becomes biological absence",
            "the Wet handoff owns production membership and never forces an arbitrary Top-20",
        ],
    }
    semantic_root.mkdir(parents=True, exist_ok=True)
    extension_manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        "state": manifest["state"],
        "manifest": str(extension_manifest_path),
        "prompts": str(bundle_result["prompts"]),
        "route_b_verified_candidates": str(frozen_b["candidates"]),
        "canonical_family_recovery": str(canonical_family["family_recovery"]),
        "route_b_review_candidates": frozen_b.get("review_candidates"),
        "family_copy": family_copy,
        "wet_handoff_manifest": str(handoff["manifest"]),
        "wet_candidates": handoff.get("wet_candidates"),
        "membership_state": handoff.get("membership_state"),
        "input_coverage_state": handoff.get("input_coverage_state"),
        "legacy_b1_consumed": False,
        "legacy_triple_route_family_context_consumed": False,
    }
