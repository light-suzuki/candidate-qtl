"""CLI commands for deterministic compute/candidate budgets."""

from __future__ import annotations

from argparse import _SubParsersAction
from pathlib import Path
from typing import Any

from ..adaptive_synthesis import (
    apply_verified_refinement,
    build_adaptive_synthesis_request,
    import_adaptive_synthesis_proposals,
)
from ..budget_account import BUDGET_PROFILES
from ..budget_engine import run_budget_engine
from ..candidate_budget import candidate_budget_gate
from ..core import runtime_path
from ..llm_usage import summarize_llm_usage
from ..supported_member_union import build_supported_member_union
from ..tier0_inventory import load_provider_table
from ..wet_candidate_handoff import wet_candidate_handoff
from ._shared import add_json_flag

# Flag name -> provider kind consumed by tier0_inventory.load_provider_table.
# --pav-evidence is intentionally not exposed: Tier 0 supports the kind, but
# issue #375's CLI contract lists only expression/domains/orthogroup/family/
# variant/sv/promoter.
_TABLE_FLAGS = (
    ("expression", "expression", "expression_provider"),
    ("domains", "domains", "domains_provider"),
    ("orthogroup", "orthogroup", "orthogroup_provider"),
    ("family", "family_table", "family_provider"),
    ("variant", "variant_evidence", "variant_provider"),
    ("sv", "sv_evidence", "sv_provider"),
    ("promoter", "promoter_evidence", "promoter_provider"),
)


def _provider_spec(value: str) -> tuple[str, str]:
    """Split a NAME@VERSION provider label; malformed labels are rejected."""
    if "@" not in value:
        raise ValueError(f"provider must be NAME@VERSION, got {value!r}")
    name, version = value.split("@", 1)
    if not name or not version:
        raise ValueError(f"provider must be NAME@VERSION, got {value!r}")
    return name, version


def register(subparsers: _SubParsersAction) -> None:
    sub = subparsers.add_parser("budget-run", help="Deterministic budgeted engine: tier0 inventory, tier1 clusters, tier2 plan, handoff, resource-usage ledger (issue #375).")
    sub.add_argument("--case-id", default="case_001", help="Case identifier recorded in every manifest.")
    sub.add_argument("--extract-dir", required=True, help="Extract directory produced by `extract` (candidates.tsv etc.).")
    sub.add_argument("--out-dir", default=None, help="Output directory; default <extract-dir>/budget_engine.")
    for _kind, _flag, _provider_flag in _TABLE_FLAGS:
        sub.add_argument(f"--{_flag}", help=f"{_kind} provider table TSV.")
        sub.add_argument(f"--{_provider_flag}", default="unknown@0", help=f"{_kind} provider label NAME@VERSION.")
    sub.add_argument("--budget-profile", choices=sorted(BUDGET_PROFILES), default="student_cpu", help="Budget profile preset.")
    sub.add_argument("--cache-dir", default=None, help="Content-addressed cache directory; default <out-dir>/budget_cache.")
    sub.add_argument("--no-cache", action="store_true", help="Disable the content cache (recompute every query).")
    sub.add_argument("--protein-db", default=None, help="Proteome database file or directory searched by blastp; absent means provider missing.")
    sub.add_argument("--blastp-bin", default="blastp", help="blastp executable name or path.")
    sub.add_argument("--flag-gene", action="append", help="Force the Tier 2 expansion trigger for a gene; repeatable.")
    sub.add_argument("--exhaustive", action="store_true", help="Also run the unlimited exhaustive arm and the gating comparison.")
    sub.add_argument("--bridge", default=None, help="family_candidate_bridge.tsv (from `family-candidate-bridge`); supplies the route each candidate rode in on so the portfolio reducer can allocate its route-unique rescue slots. Without it the reducer records no_route_information_available instead of assuming a route.")
    sub.add_argument("--portfolio-size", type=int, default=20, help="Top-N portfolio size for the deterministic/LLM-assisted reducer (issue #375 T5).")
    add_json_flag(sub)

    sub = subparsers.add_parser(
        "supported-member-union",
        help="Build the non-scoring union of supported Route-A bridge members and final verified Route-B semantic members (#426).",
    )
    sub.add_argument("--route-a-bridge", help="Frozen family_candidate_bridge.tsv; only route=a rows that are not tier6/not_evaluated contribute.")
    sub.add_argument("--route-b-claims", help="Final semantic claims TSV; only verification_status=verified rows contribute all member_interval_genes.")
    sub.add_argument("--route-b-candidates", help="Frozen route_b_verified_candidates.tsv; only route=b and support_status=verified_semantic_support rows contribute.")
    sub.add_argument("--route-c-members", help="Frozen Direct-C member ledger; only support_status=verified_direct_support rows contribute.")
    sub.add_argument("--route-a-manifest", help="Producer manifest for Route A bridge (required for canonical input).")
    sub.add_argument("--route-b-claims-manifest", help="Producer manifest for Route B claims.")
    sub.add_argument("--route-b-candidates-manifest", help="Producer manifest for Route B frozen candidates.")
    sub.add_argument("--route-c-manifest", help="route_c_direct_members_manifest.json.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)

    sub = subparsers.add_parser(
        "candidate-budget-gate",
        help="Freeze every supported candidate; report the Wet target excess without truncation or a count-only stop.",
    )
    sub.add_argument("--candidates", required=True, help="Supported member TSV with candidate_gene_id or gene_id. Membership is never truncated.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--wet-budget", type=int, default=20, help="Operational review target, not a scientific cutoff.")
    sub.add_argument("--budget-policy", choices=("target", "strict"), default="target", help="Use strict only for an explicitly fixed resource capacity.")
    sub.add_argument("--adaptive-round", type=int, default=0)
    sub.add_argument("--max-adaptive-rounds", type=int, default=1)
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)

    sub = subparsers.add_parser(
        "wet-candidate-handoff",
        help="Canonical path from frozen Route A/B artifacts to a Wet-reviewable candidate set and stage_a_narrowed (#426).",
    )
    sub.add_argument("--route-a-bridge", help="Frozen family_candidate_bridge.tsv.")
    sub.add_argument("--route-b-claims", help="Final semantic claims TSV.")
    sub.add_argument("--route-b-candidates", help="Frozen route_b_verified_candidates.tsv.")
    sub.add_argument("--route-c-members", help="Frozen Direct-C member ledger. Verified direct support enters the union; review-only rows may fill tier2 Wet slots.")
    sub.add_argument("--route-a-manifest", help="Producer manifest for Route A bridge.")
    sub.add_argument("--route-b-claims-manifest", help="Producer manifest for Route B claims.")
    sub.add_argument("--route-b-candidates-manifest", help="Producer manifest for Route B frozen candidates.")
    sub.add_argument("--route-c-manifest", help="route_c_direct_members_manifest.json.")
    sub.add_argument("--refined-candidates", help="Machine-verified refined candidate set; re-entry requires adaptive-round >=1 and an original canonical union above the Wet budget.")
    sub.add_argument("--refinement-manifest", help="adaptive_refinement_manifest.json binding the verification and refined output.")
    sub.add_argument("--refinement-original-candidates", help="Original supported member union consumed by adaptive-refinement-apply.")
    sub.add_argument("--refinement-original-union-manifest", help="Canonical supported_member_union_manifest.json binding the original candidate set.")
    sub.add_argument("--ranked", help="Deterministic ranked_candidates.tsv. Read-only ordering inside the complete frozen Wet membership.")
    sub.add_argument("--family-recovery", help="Optional family recovery handoff for narrowing.")
    sub.add_argument("--convergence", help="Optional dual/triple-route convergence table for narrowing.")
    sub.add_argument("--family-copy", help="Optional family-to-copy table for narrowing.")
    sub.add_argument("--family-copy-shortlist", help="Optional family-copy shortlist for narrowing.")
    sub.add_argument("--ambiguity-context", help="Optional per-candidate ambiguity context handed to the adaptive synthesis request.")
    sub.add_argument("--interval-atlas", help="Frozen interval_gene_atlas.tsv. Enables tiered review-slot filling when the supported set is under the Wet budget.")
    sub.add_argument("--review-candidates", action="append", help="Sealed A/B review candidates TSV; repeat for multiple routes. Retained ahead of generic fill without verified support.")
    sub.add_argument("--review-candidates-manifest", action="append", help="Matching review producer manifest; repeat in the same order as --review-candidates.")
    sub.add_argument("--identity-assessment", help="Explicit hash-bound post-freeze parent/independent identity review; derives member-bound Tier evidence without changing support or frozen papers.")
    sub.add_argument("--identity-target-proteins", help="Canonical interval candidate protein FASTA; required with --identity-assessment to verify the exact target sequence.")
    sub.add_argument("--reference-function", help="Optional reference function/mutant-phenotype corpus; lets a filled candidate be tiered as tier2_reference_phenotype.")
    sub.add_argument("--wet-fill-mode", choices=("tiered_interval_fill", "supported_only"), default="tiered_interval_fill", help="How unused Wet review slots are filled. Filled candidates carry no support and no score.")
    sub.add_argument("--deterministic-subfamily-refinement", action="store_true",
                     help="When the gate reports over-budget, apply one machine-verified curated cross-subfamily refinement pass (requires --family-copy) and re-gate once instead of only emitting an adaptive synthesis request.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--wet-budget", type=int, default=20, help="Operational review target, not a scientific cutoff.")
    sub.add_argument("--budget-policy", choices=("target", "strict"), default="target", help="Use strict only for an explicitly fixed resource capacity.")
    sub.add_argument("--adaptive-round", type=int, default=0)
    sub.add_argument("--max-adaptive-rounds", type=int, default=1)
    sub.add_argument("--min-route-candidates", type=int, default=3,
                     help="Operational A/C-first signal: a final Wet set below this asks the parent to rethink Route A in a fresh generation. Never pads or trims membership.")
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)

    sub = subparsers.add_parser(
        "adaptive-synthesis-request",
        help="Emit the bounded ambiguity-resolution request for an over-budget candidate set; refuses to run within budget (#426/#434).",
    )
    sub.add_argument("--candidates", required=True)
    sub.add_argument("--budget-manifest", required=True, help="candidate_budget_manifest.json; must report adaptive_refinement_required.")
    sub.add_argument("--ambiguity-context", help="Optional per-candidate ambiguity context.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)

    sub = subparsers.add_parser(
        "adaptive-synthesis-import",
        help="Normalize reviewer ambiguity-resolution proposals; rank_contribution stays none and no proposal narrows the set (#426).",
    )
    sub.add_argument("--proposals", help="Reviewer proposals TSV. Omit to record not_evaluated_no_provider.")
    sub.add_argument("--request", required=True, help="adaptive_synthesis_request.json this response answers.")
    sub.add_argument("--provider", default="")
    sub.add_argument("--model", default="")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)

    sub = subparsers.add_parser(
        "adaptive-refinement-apply",
        help="Narrow a candidate set using machine-verified verdicts only; unsourced or not_evaluated verdicts never remove a candidate (#426).",
    )
    sub.add_argument("--candidates", required=True)
    sub.add_argument("--verification", help="Machine verification TSV with candidate_gene_id and verdict. Omit to leave the set untouched.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)

    sub = subparsers.add_parser(
        "llm-usage",
        help="Aggregate provider-reported token usage/call counts from llm-run raw responses without estimating missing tokens or prices (#434).",
    )
    sub.add_argument("--raw-responses", required=True, help="llm_raw_responses.jsonl from llm-run.")
    sub.add_argument("--llm-manifest", help="Optional llm_run_manifest.json for provider/model/prompt-count validation.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)


def _run_budget_run(args: Any) -> dict:
    providers: dict[str, dict[str, object]] = {}
    for kind, flag, provider_flag in _TABLE_FLAGS:
        table_path = getattr(args, flag)
        if table_path:
            provider, version = _provider_spec(getattr(args, provider_flag))
            providers[kind] = load_provider_table(kind, table_path, provider, version)
    out_dir = args.out_dir or str(Path(runtime_path(args.extract_dir)) / "budget_engine")
    manifest = run_budget_engine(
        args.extract_dir,
        out_dir,
        args.case_id,
        BUDGET_PROFILES[args.budget_profile],
        providers,
        cache_dir=args.cache_dir,
        no_cache=args.no_cache,
        protein_db=args.protein_db,
        blastp_bin=args.blastp_bin,
        flag_genes=args.flag_gene,
        exhaustive=args.exhaustive,
        bridge=args.bridge,
        portfolio_size=args.portfolio_size,
    )
    return {
        "schema_version": "1.0",
        "status": manifest["status"],
        "producer_class": "diagnostic_portfolio",
        "production_membership": False,
        "outputs": manifest["outputs"],
        "exit_code": 0,
    }


def _run_supported_member_union(args: Any) -> dict:
    return build_supported_member_union(
        args.out_dir,
        route_a_bridge_arg=args.route_a_bridge,
        route_b_claims_arg=args.route_b_claims,
        route_b_candidates_arg=args.route_b_candidates,
        route_c_members_arg=args.route_c_members,
        route_a_manifest_arg=args.route_a_manifest,
        route_b_claims_manifest_arg=args.route_b_claims_manifest,
        route_b_candidates_manifest_arg=args.route_b_candidates_manifest,
        route_c_manifest_arg=args.route_c_manifest,
        case_id=args.case_id,
    )


def _run_candidate_budget_gate(args: Any) -> dict:
    return candidate_budget_gate(
        args.candidates,
        args.out_dir,
        wet_budget=args.wet_budget,
        budget_policy=args.budget_policy,
        adaptive_round=args.adaptive_round,
        max_adaptive_rounds=args.max_adaptive_rounds,
        case_id=args.case_id,
    )


def _run_wet_candidate_handoff(args: Any) -> dict:
    return wet_candidate_handoff(
        args.out_dir,
        route_a_bridge_arg=args.route_a_bridge,
        route_b_claims_arg=args.route_b_claims,
        route_b_candidates_arg=args.route_b_candidates,
        route_c_members_arg=args.route_c_members,
        route_a_manifest_arg=args.route_a_manifest,
        route_b_claims_manifest_arg=args.route_b_claims_manifest,
        route_b_candidates_manifest_arg=args.route_b_candidates_manifest,
        route_c_manifest_arg=args.route_c_manifest,
        ranked_arg=args.ranked,
        family_recovery_arg=args.family_recovery,
        convergence_arg=args.convergence,
        family_copy_arg=args.family_copy,
        family_copy_shortlist_arg=args.family_copy_shortlist,
        refined_candidates_arg=args.refined_candidates,
        refinement_manifest_arg=args.refinement_manifest,
        refinement_original_candidates_arg=args.refinement_original_candidates,
        refinement_original_union_manifest_arg=args.refinement_original_union_manifest,
        wet_budget=args.wet_budget,
        budget_policy=args.budget_policy,
        adaptive_round=args.adaptive_round,
        max_adaptive_rounds=args.max_adaptive_rounds,
        ambiguity_arg=args.ambiguity_context,
        interval_atlas_arg=args.interval_atlas,
        review_candidates_arg=args.review_candidates,
        review_candidates_manifest_arg=args.review_candidates_manifest,
        identity_assessment_arg=getattr(args, "identity_assessment", None),
        identity_target_proteins_arg=getattr(args, "identity_target_proteins", None),
        min_route_candidates=getattr(args, "min_route_candidates", 3),
        reference_function_arg=args.reference_function,
        wet_fill_mode=args.wet_fill_mode,
        deterministic_subfamily_refinement=bool(getattr(args, "deterministic_subfamily_refinement", False)),
        case_id=args.case_id,
    )


def _run_adaptive_synthesis_request(args: Any) -> dict:
    return build_adaptive_synthesis_request(
        args.candidates,
        args.budget_manifest,
        args.out_dir,
        ambiguity_arg=args.ambiguity_context,
        case_id=args.case_id,
    )


def _run_adaptive_synthesis_import(args: Any) -> dict:
    return import_adaptive_synthesis_proposals(
        args.proposals,
        args.request,
        args.out_dir,
        provider=args.provider,
        model=args.model,
        case_id=args.case_id,
    )


def _run_adaptive_refinement_apply(args: Any) -> dict:
    return apply_verified_refinement(
        args.candidates,
        args.verification,
        args.out_dir,
        case_id=args.case_id,
    )


def _run_llm_usage(args: Any) -> dict:
    return summarize_llm_usage(
        args.raw_responses,
        args.out_dir,
        llm_manifest_arg=args.llm_manifest,
        case_id=args.case_id,
    )


HANDLERS = {
    "budget-run": _run_budget_run,
    "supported-member-union": _run_supported_member_union,
    "candidate-budget-gate": _run_candidate_budget_gate,
    "wet-candidate-handoff": _run_wet_candidate_handoff,
    "adaptive-synthesis-request": _run_adaptive_synthesis_request,
    "adaptive-synthesis-import": _run_adaptive_synthesis_import,
    "adaptive-refinement-apply": _run_adaptive_refinement_apply,
    "llm-usage": _run_llm_usage,
}
