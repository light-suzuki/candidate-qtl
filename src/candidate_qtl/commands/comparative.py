"""Comparative-genomics commands and optional heavy rescue lanes.

Registration and dispatch for comparative reference resources, external rank
ledgers, deterministic term propagation, Route-B semantic bundles/proposals/
source+transfer verification/claims, and comparative-locus reconstruction live
here.
"""

from __future__ import annotations

from argparse import _SubParsersAction
from typing import Any

import pandas as pd

from ..comparative_locus import comparative_locus_benchmark, comparative_locus_reconstruct
from ..core import runtime_path
from ..external_rank_ledger import external_rank_ledger
from ..reference_resource_panel import reference_resource_panel
from ..route_b_semantic_bundle import build_route_b_semantic_bundle
from ..route_b_semantic_candidates import freeze_route_b_semantic_candidates
from ..route_b_semantic_import import import_route_b_semantic_proposals
from ..route_b_semantic_source_verify import verify_route_b_semantic_sources
from ..route_b_semantic_transfer_verify import verify_route_b_semantic_transfer
from ..semantic_claims import freeze_semantic_claims
from ..term_propagation import propagate_terms
from ._shared import add_json_flag


def register(subparsers: _SubParsersAction) -> None:
    sub = subparsers.add_parser("reference-panel", help="Validate and availability-audit a comparative reference-resource panel (native required, others optional with machine-readable skip reasons).")
    sub.add_argument("--panel", required=True, help="TSV declaring resource_id/role/resource_type/taxon/source/uri_or_path per resource.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--full-hash", action="store_true")
    add_json_flag(sub)

    sub = subparsers.add_parser("external-rank", help="Import external comparative ranks into a coverage-aware ledger (RRF consensus) without touching the native dry_score.")
    sub.add_argument("--native-rank", required=True, help="Native ranked_candidates.tsv (source of the candidate universe and native rank).")
    sub.add_argument("--external-ranks", required=True, help="Long-format TSV: one row per (ranker_id, gene_id) with coverage.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--top-k", type=int, default=20, help="Top-K threshold used to compare native and consensus prioritization.")
    sub.add_argument("--rrf-k", type=int, default=60, help="Reciprocal-rank-fusion constant.")
    sub.add_argument("--discordance-spread", type=float, default=0.5, help="Normalized cross-family rank spread that flags disagreement.")
    sub.add_argument("--full-hash", action="store_true")
    add_json_flag(sub)

    sub = subparsers.add_parser("propagate-terms", help="Deterministically propagate reference-ortholog function terms onto QTL interval candidates (no LLM); native-only input degrades gracefully.")
    sub.add_argument("--candidates", required=True, help="Interval candidates.tsv (gene_id universe).")
    sub.add_argument("--orthology", help="TSV mapping candidate_gene_id -> reference_gene_id (+ optional relationship/source_genome/reference_role/reference_species/confidence).")
    sub.add_argument("--reference-annotation", help="TSV of reference_gene_id -> term (+ optional ontology_id/ontology_name/source).")
    sub.add_argument("--trait-terms", help="Optional analyst trait terms (one per line or a term column) for trait_relevance tagging.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--full-hash", action="store_true")
    add_json_flag(sub)

    sub = subparsers.add_parser(
        "route-b-semantic-bundle",
        help="After B0 freeze, compress the complete interval atlas into identity-preserving family/orthogroup/gene units and emit provider-neutral semantic prompts (#403).",
    )
    sub.add_argument("--atlas", help="Frozen interval_gene_atlas.tsv; required for an initial bundle.")
    sub.add_argument("--phenotype", help="Phenotype introduced only after B0; required for an initial bundle.")
    sub.add_argument("--extend-bundle", help="Frozen parent bundle to extend with newly retrieved sources; never edits the parent.")
    sub.add_argument("--parent-bundle-sha256", help="Coordinator-pinned parent hash for --extend-bundle.")
    sub.add_argument("--research-discoveries", help="Worker JSON with case_id, parent_bundle_sha256 and discovered_sources; not claim responses.")
    sub.add_argument("--focal-species", default="")
    sub.add_argument("--b-review-membership-policy", choices=["initial_bibliography_v1", "source_function_v1"], default="initial_bibliography_v1")
    sub.add_argument("--literature", help="Optional reverse_literature.tsv link candidates; retrieval is not evidence.")
    sub.add_argument("--literature-status", help="Optional reverse_literature_status.tsv preserving no-hit/not-evaluated states.")
    sub.add_argument("--lookup", help="Optional interval_reverse_lookup.tsv carrying source-backed B0 function/pathway/source fields into the frozen semantic unit.")
    sub.add_argument("--identity-hits", help="Optional route_b_identity_hits.tsv preserving all B0 sequence-similarity retrieval handles for each interval gene; these are not orthology or trait support.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--case-id", default="")
    sub.add_argument("--full-hash", action="store_true")
    add_json_flag(sub)

    sub = subparsers.add_parser(
        "route-b-semantic-import",
        help="Import normalized llm-run responses for frozen Route-B semantic prompts. Unknown sources/reference genes are unsupported; all other rows remain proposals pending relation verification (#403).",
    )
    sub.add_argument("--bundle", required=True, help="route_b_semantic_bundle.json")
    sub.add_argument("--prompts", required=True, help="route_b_semantic_prompts.jsonl used by llm-run")
    sub.add_argument("--responses", required=True, help="llm_candidate_responses.json from llm-run")
    sub.add_argument("--llm-manifest", help="Optional llm_run_manifest.json supplying provider/model metadata.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)

    sub = subparsers.add_parser(
        "route-b-semantic-source-verify",
        help="Verify proposal citations/source relations with the existing source-verify engine while keeping focal-family transfer as a separate unresolved axis (#403).",
    )
    sub.add_argument("--proposals", required=True, help="semantic_claim_proposals.tsv from route-b-semantic-import.")
    sub.add_argument("--bundle", help="Optional frozen route_b_semantic_bundle.json. Enables source-backed curated-function verification for database accessions grounded in the exact B0 unit; this is weaker than mechanistic literature and does not verify transfer.")
    sub.add_argument("--cache", help="Optional local literature-record cache used by source-verify.")
    sub.add_argument("--allow-external-api", action="store_true", help="Allow source-verify to resolve citation identifiers through Europe PMC with an official NCBI PubMed/PMC E-utilities fallback. Only citation identifiers are transmitted.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--case-id", default="")
    sub.add_argument("--full-hash", action="store_true")
    add_json_flag(sub)

    sub = subparsers.add_parser(
        "route-b-semantic-transfer-verify",
        help="Combine source-level mechanistic support with frozen family/orthogroup/copy identity. Group-level claims may verify without choosing one copy; interval_gene claims require one-to-one/RBH-like transfer (#403).",
    )
    sub.add_argument("--proposals", required=True, help="semantic_claim_proposals.tsv")
    sub.add_argument("--source-verification", required=True, help="semantic_claim_source_verification.tsv")
    sub.add_argument("--bundle", required=True, help="route_b_semantic_bundle.json")
    sub.add_argument("--bound-claims", help="Shared bound-claim TSV. Without a complete matching row, source/cue evidence cannot authorize canonical membership.")
    sub.add_argument("--agent-claim-review", help="Independent Route-B review packet for the exact semantic proposal/source snapshot.")
    sub.add_argument("--source-review-records", help="source_review_records.json matching the exact Route-B proposal citations.")
    sub.add_argument("--record-cache", "--records", dest="record_cache", help="Frozen UniProt-shaped provider records for exact accession/species checks.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)

    sub = subparsers.add_parser(
        "route-b-semantic-candidates",
        help="Freeze Route-B semantic candidate members from verified claims before comparison (#403).",
    )
    sub.add_argument("--claims", required=True, help="semantic_claims.tsv (or semantic_claims_frozen.tsv).")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)

    sub = subparsers.add_parser(
        "semantic-claims-freeze",
        help="Validate and freeze structured Route-B semantic claims; only verified source-backed rows are narrowing-eligible (#403).",
    )
    sub.add_argument("--claims", required=True, help="Structured semantic claims TSV following candidate_qtl.semantic_claims.CLAIM_COLUMNS.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--case-id", default="")
    sub.add_argument("--route", default="interval_first", help="Expected route label for every claim (default interval_first).")
    add_json_flag(sub)

    sub = subparsers.add_parser("comparative-locus-reconstruct", help="Reconstruct an optional candidate universe from flanking markers and precomputed multi-species locus caches; never adds causal score.")
    sub.add_argument("--markers", required=True, help="FASTA with at least two flanking marker sequences; native chromosome coordinates are not required.")
    sub.add_argument("--reference-panel", required=True, help="TSV: species/role/available/cache_version (+ optional phylogenetic_distance and cache paths).")
    sub.add_argument("--marker-projections", required=True, help="Precomputed marker hits: marker_id/species/seqid/start (+ optional orientation/hit_count).")
    sub.add_argument("--reference-locus-genes", required=True, help="Cached projected-block genes: species/gene_id/orthogroup/order_index.")
    sub.add_argument("--focal-recovery", help="Precomputed focal assembly/scaffold/transcriptome hits: orthogroup/focal_gene_id/source_class.")
    sub.add_argument("--native-genes", help="Optional native annotation genes: gene_id/orthogroup; retained as the primary candidate class.")
    sub.add_argument("--absence-evidence", help="Optional callable absence evidence: orthogroup/evidence_status/source_id.")
    sub.add_argument("--min-reference-species", type=int, default=2)
    sub.add_argument("--core-support-fraction", type=float, default=0.67)
    sub.add_argument("--case-id", default="")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--full-hash", action="store_true")
    add_json_flag(sub)

    sub = subparsers.add_parser("comparative-locus-benchmark", help="Evaluate frozen comparative reconstructions on synthetic/degraded-reference cases.")
    sub.add_argument("--manifest", required=True, help="TSV: case_id/degradation_type/truth_orthogroup/reconstructed_universe_tsv/reconstruction_manifest_json.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--full-hash", action="store_true")
    add_json_flag(sub)


def _run_reference_panel(args: Any) -> dict:
    return reference_resource_panel(args.panel, args.out_dir, args.full_hash)


def _run_external_rank(args: Any) -> dict:
    return external_rank_ledger(args.native_rank, args.external_ranks, args.out_dir, args.top_k, args.rrf_k, args.discordance_spread, args.full_hash)


def _run_propagate_terms(args: Any) -> dict:
    return propagate_terms(args.candidates, args.orthology, args.reference_annotation, args.out_dir, args.trait_terms, args.full_hash)


def _run_route_b_semantic_bundle(args: Any) -> dict:
    if args.extend_bundle or args.research_discoveries or args.parent_bundle_sha256:
        if not (args.extend_bundle and args.research_discoveries and args.parent_bundle_sha256):
            raise ValueError("source re-freeze requires --extend-bundle, --parent-bundle-sha256 and --research-discoveries")
        if any((args.atlas, args.phenotype, args.literature, args.literature_status, args.lookup, args.identity_hits, args.case_id, args.focal_species)):
            raise ValueError("source re-freeze inherits frozen case/phenotype/identity; do not supply initial-bundle inputs")
        from ..route_b_semantic_bundle import refreeze_route_b_sources
        return refreeze_route_b_sources(args.extend_bundle, args.research_discoveries, args.out_dir, parent_sha256=args.parent_bundle_sha256)
    if not args.atlas or not args.phenotype:
        raise ValueError("initial semantic bundle requires --atlas and --phenotype")
    return build_route_b_semantic_bundle(
        args.atlas,
        args.out_dir,
        phenotype=args.phenotype,
        focal_species=args.focal_species,
        b_review_membership_policy=args.b_review_membership_policy,
        literature_arg=args.literature,
        literature_status_arg=args.literature_status,
        lookup_arg=args.lookup,
        identity_hits_arg=args.identity_hits,
        case_id=args.case_id,
        full_hash=args.full_hash,
    )


def _run_route_b_semantic_import(args: Any) -> dict:
    return import_route_b_semantic_proposals(
        args.bundle,
        args.prompts,
        args.responses,
        args.out_dir,
        llm_manifest_arg=args.llm_manifest,
        case_id=args.case_id,
    )


def _run_route_b_semantic_source_verify(args: Any) -> dict:
    return verify_route_b_semantic_sources(
        args.proposals,
        args.out_dir,
        bundle_arg=args.bundle,
        cache_arg=args.cache,
        allow_external_api=args.allow_external_api,
        full_hash=args.full_hash,
        case_id=args.case_id,
    )


def _run_route_b_semantic_transfer_verify(args: Any) -> dict:
    return verify_route_b_semantic_transfer(
        args.proposals,
        args.source_verification,
        args.bundle,
        args.out_dir,
        bound_claims_arg=args.bound_claims,
        agent_claim_review_arg=args.agent_claim_review,
        source_review_records_arg=args.source_review_records,
        record_cache_arg=args.record_cache,
        case_id=args.case_id,
    )


def _run_route_b_semantic_candidates(args: Any) -> dict:
    return freeze_route_b_semantic_candidates(
        args.claims,
        args.out_dir,
        case_id=args.case_id,
    )


def _run_semantic_claims_freeze(args: Any) -> dict:
    claims_path = runtime_path(args.claims)
    if not claims_path.exists():
        raise ValueError(f"semantic claims input does not exist: {claims_path}")
    frame = pd.read_csv(claims_path, sep="\t", dtype=str).fillna("")
    return freeze_semantic_claims(
        frame,
        args.out_dir,
        case_id=args.case_id,
        route=args.route,
    )


def _run_comparative_locus_reconstruct(args: Any) -> dict:
    return comparative_locus_reconstruct(
        args.markers, args.reference_panel, args.marker_projections, args.reference_locus_genes,
        args.out_dir, args.focal_recovery, args.native_genes, args.absence_evidence,
        args.min_reference_species, args.core_support_fraction, args.case_id, args.full_hash,
    )


def _run_comparative_locus_benchmark(args: Any) -> dict:
    return comparative_locus_benchmark(args.manifest, args.out_dir, args.full_hash)


HANDLERS = {
    "reference-panel": _run_reference_panel,
    "external-rank": _run_external_rank,
    "propagate-terms": _run_propagate_terms,
    "route-b-semantic-bundle": _run_route_b_semantic_bundle,
    "route-b-semantic-import": _run_route_b_semantic_import,
    "route-b-semantic-source-verify": _run_route_b_semantic_source_verify,
    "route-b-semantic-transfer-verify": _run_route_b_semantic_transfer_verify,
    "route-b-semantic-candidates": _run_route_b_semantic_candidates,
    "semantic-claims-freeze": _run_semantic_claims_freeze,
    "comparative-locus-reconstruct": _run_comparative_locus_reconstruct,
    "comparative-locus-benchmark": _run_comparative_locus_benchmark,
}
