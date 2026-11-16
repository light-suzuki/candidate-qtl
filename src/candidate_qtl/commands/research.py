"""CLI commands for the research feature group (issue #144)."""

from __future__ import annotations

from ..core import candidate_research, europe_pmc_search
from ..phenotype_investigation import (
    MIN_COMPETING_MECHANISMS, PAYLOAD_SOURCES, investigation_attribution, mechanism_import, phenotype_investigate,
)
from ..reference_gene_resolution import DEFAULT_TIMEOUT as REFERENCE_GENE_RESOLUTION_TIMEOUT
from ..route_a_coverage import (
    mechanism_coverage_audit,
    phenotype_localization_freeze,
    phenotype_localization_packet,
)
from ..source_verify import claim_verification_benchmark, source_verify
from argparse import _SubParsersAction
from typing import Any
from ._shared import add_json_flag


def register(subparsers: _SubParsersAction) -> None:
    sub = subparsers.add_parser("phenotype-investigate", help="Route A stage 1: build the interval-blind phenotype research packet (direction/inheritance/tissue/stage/onset/penetrance/synonyms/ontology/uncertainty) plus a prompt package and a blank payload template an LLM or a human fills identically.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--case-id", required=True)
    sub.add_argument("--focal-species", required=True)
    sub.add_argument("--phenotype", required=True)
    sub.add_argument("--direction", default="unknown")
    sub.add_argument("--inheritance", default="unknown")
    sub.add_argument("--tissue", default="", help="Affected tissue/organ.")
    sub.add_argument("--stage", default="", help="Developmental stage.")
    sub.add_argument("--onset", default="")
    sub.add_argument("--penetrance", default="")
    sub.add_argument("--synonym", action="append", help="Phenotype synonym; repeatable.")
    sub.add_argument("--ontology-term", action="append", help="Trait/phenotype ontology term; repeatable.")
    sub.add_argument("--exclusion", action="append", help="Known exclusion; repeatable.")
    sub.add_argument("--constraint", action="append", help="Researcher constraint; repeatable.")
    sub.add_argument("--uncertainty", default="", help="What is genuinely uncertain about the phenotype description.")
    sub.add_argument("--candidates", help="Leak check only: the packet is rejected if it names an interval candidate id. Never inserted into the packet.")
    sub.add_argument("--benchmark-mode", action="store_true", help="Apply strict benchmark reference policy: reference targets must be non-focal.")
    add_json_flag(sub)

    sub = subparsers.add_parser("phenotype-localization-packet", help="Route A coverage stage A: create a phenotype-only prompt for a fresh causal-localization assessment before mechanism generation.")
    sub.add_argument("--investigation", required=True, help="phenotype_investigation.tsv from phenotype-investigate.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--candidates", help="Leak check only; candidate ids are never inserted into the prompt.")
    add_json_flag(sub)

    sub = subparsers.add_parser("phenotype-localization-freeze", help="Route A coverage stage B: strictly validate and hash-freeze the phenotype-only localization response before mechanism generation.")
    sub.add_argument("--payload", required=True, help="JSON with causal_localization/localization_confidence/supporting_phrase.")
    sub.add_argument("--investigation", required=True)
    sub.add_argument("--out-dir", required=True)
    add_json_flag(sub)

    sub = subparsers.add_parser("mechanism-import", help="Route A stage 2: import an LLM or human mechanism payload through one strict schema. Fewer than the required competing mechanisms, or a hypothesis without a counter-hypothesis, discriminating observation, source id or model-species evidence, is review-required.")
    sub.add_argument("--payload", required=True, help="Model output (.json/.jsonl/.txt) or a hand-authored .tsv following mechanism_payload_template.tsv.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--investigation", help="phenotype_investigation.tsv from phenotype-investigate.")
    sub.add_argument("--payload-source", choices=list(PAYLOAD_SOURCES), default="llm")
    sub.add_argument("--model", default="", help="Model identifier when the payload came from one.")
    sub.add_argument("--case-id", default="")
    sub.add_argument("--retrieval-provenance", default="", help="How the sources were retrieved (database, version, date).")
    sub.add_argument("--candidates", help="Leak check only: a payload naming an interval candidate id is rejected as not blind.")
    sub.add_argument("--min-mechanisms", type=int, default=MIN_COMPETING_MECHANISMS)
    sub.add_argument("--proposal-round", type=int, default=1, help="Generation round. A phenotype-only gap fill may be round 2 and remain blind.")
    sub.add_argument("--interval-informed", action="store_true", help="Mark that the payload saw an interval verdict. Never inferred from proposal-round.")
    sub.add_argument("--mechanism-freeze-manifest", help="TSV with case_id/mechanism_payload_sha256 (coordination/payload_freeze.tsv format); the payload bytes must match the frozen hash or the import fails closed.")
    sub.add_argument("--localization-freeze-manifest", help="Optional phenotype_localization_freeze_manifest.json; binds this import to the phenotype-only localization frozen before generation.")
    sub.add_argument("--benchmark-mode", action="store_true", help="Apply strict benchmark reference policy: reference targets must be non-focal.")
    add_json_flag(sub)

    sub = subparsers.add_parser("mechanism-coverage-audit", help="Route A coverage stage C: machine-audit 12 initial labeled hypotheses, emit a blind gap-fill prompt for at most 3 rows, or merge one imported gap-fill table and finalize at no more than 15.")
    sub.add_argument("--mechanisms", action="append", required=True, help="mechanism_hypotheses.tsv. Pass initial first; optionally repeat once for an imported gap-fill payload.")
    sub.add_argument("--mechanism-import-manifest", action="append", required=True, help="Matching mechanism_import_manifest.json; repeat in the same order as --mechanisms.")
    sub.add_argument("--localization-freeze-manifest", required=True)
    sub.add_argument("--investigation", required=True)
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--candidates", help="Leak check only; candidate ids are never inserted into the prompt.")
    sub.add_argument("--gap-request-audit", help="Required with a second --mechanisms table; binds the gap fill to the prior deterministic request.")
    add_json_flag(sub)

    sub = subparsers.add_parser("reference-gene-resolve", help="Route A stage 2.5: retain provider-linked reference-gene accession alternatives before source-verify, with parent-hypothesis lineage and per-accession citations. Gene-name retrieval includes unreviewed records and has no candidate-count cap. Existing supplied sources and explicit identifiers retain their verification path; provider failures remain not_evaluated.")
    sub.add_argument("--mechanisms", help="mechanism_hypotheses.tsv from mechanism-import; required for blind name resolution.")
    sub.add_argument("--out", required=True, help="Output directory.")
    sub.add_argument("--allow-external-api", action="store_true", help="Permit read-only UniProt REST lookups.")
    sub.add_argument("--records", help="Local UniProt-shaped JSON record file or directory for offline resolution (same format as hypothesis-verify's --record-cache).")
    sub.add_argument("--timeout", type=int, default=REFERENCE_GENE_RESOLUTION_TIMEOUT)
    sub.add_argument("--focal-species", default="", help="Optional blinding diagnostic only: flags a malformed proposed_stable_id that itself names the focal species as its own distinct review status. Not a verification gate -- that gate is hypothesis-verify's --focal-species.")
    sub.add_argument("--sequence-requests", help="Post-freeze INSDC accession request TSV; mutually exclusive with --mechanisms.")
    sub.add_argument("--post-freeze", action="store_true", help="Explicitly run sequence correspondence after research freeze; never return interval results to blind A.")
    sub.add_argument("--proteome", help="Whole supplied reference protein FASTA, including unplaced scaffolds.")
    sub.add_argument("--protein-gene-map", help="TSV with exact protein_id and gene_id columns.")
    sub.add_argument("--runner", choices=["auto", "native", "wsl"], default="auto")
    sub.add_argument("--threads", type=int, default=4)
    add_json_flag(sub)

    sub = subparsers.add_parser("source-verify", help="Verify literature IDs/titles with Europe PMC and an official NCBI PubMed/PMC E-utilities fallback. The default initial_bibliography_v1 checks metadata only; explicit content_v1 also classifies what the source supports. An unresolvable PMID/DOI is rejected only after a completed lookup; citation_verified and claim_support_status remain separate, and provider failure is not_evaluated rather than negative evidence.")
    sub.add_argument("--mechanisms", required=True, help="mechanism_hypotheses.tsv from mechanism-import, or any TSV with an id column and a source-id column.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--cache", help="Directory (or file) of cached literature result objects for offline verification.")
    sub.add_argument("--allow-external-api", action="store_true", help="Permit read-only Europe PMC and NCBI E-utilities lookups. Only citation identifiers are transmitted.")
    sub.add_argument("--id-column", default="hypothesis_id")
    sub.add_argument("--source-column", default="source_ids")
    sub.add_argument("--verification-policy", choices=["content_v1", "initial_bibliography_v1"], default="initial_bibliography_v1",
                     help="First-stage initial_bibliography_v1 checks IDs/titles without reading abstracts; content_v1 preserves strict legacy verification.")
    sub.add_argument("--focal-species", default="")
    sub.add_argument("--reference-bridge", help="Hash-verified sequence correspondence for post-freeze parent review only; no automatic Tier/support changes.")
    sub.add_argument("--claim-extraction", help="Optional agent/human claim extraction (.tsv/.json: citation/claim_type/evidence_relation/quoted_support/locator). Each quoted span must occur in the retrieved record or the row is rejected; an extraction never promotes an unresolved citation and never writes a score.")
    sub.add_argument("--full-hash", action="store_true")
    add_json_flag(sub)

    sub = subparsers.add_parser("claim-verification-benchmark", help="Score the claim-support classifier against labelled examples (mechanistic support, associative support, irrelevant mention, fabricated citation): per-class accuracy and a confusion table, with cases whose record could not be retrieved reported not_evaluated instead of wrong.")
    sub.add_argument("--cases", required=True, help="TSV: case_id/citation/expected_claim_support_status + optional reference_gene_name/gene_family/pathway_or_process/expected_citation_verified.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--cache", help="Directory (or file) of cached Europe PMC records; without it every resolvable case is not_evaluated rather than wrong.")
    sub.add_argument("--cache-is-complete", action="store_true", help="Assert that the cache holds every record that exists for these citations, so a cache miss is a genuine non-resolution. Required to evaluate fabricated-citation cases offline.")
    sub.add_argument("--allow-external-api", action="store_true")
    sub.add_argument("--full-hash", action="store_true")
    add_json_flag(sub)

    sub = subparsers.add_parser("investigation-attribution", help="Attribute a Route A outcome to the stage where it stopped: phenotype research, verification, bait freeze, interval match, or candidate recovery. Truth-free, so it runs inside a prediction-only audit.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--mechanisms", help="mechanism_hypotheses.tsv from mechanism-import.")
    sub.add_argument("--validation", help="proposal_validation.tsv from hypothesis-verify.")
    sub.add_argument("--panel", help="verified_bait_panel.tsv from bait-panel-freeze.")
    sub.add_argument("--interval-match", help="route_c_interval_match.tsv or family_process_hypotheses.tsv carrying matched interval genes.")
    sub.add_argument("--ranked", help="Ranked candidate table for the candidate-recovery stage.")
    sub.add_argument("--case-id", default="")
    sub.add_argument("--min-mechanisms", type=int, default=MIN_COMPETING_MECHANISMS)
    add_json_flag(sub)

    sub = subparsers.add_parser("literature", help="Run one read-only Europe PMC literature query and save normalized JSON.")
    sub.add_argument("--query", required=True, help="Europe PMC search query.")
    sub.add_argument("--out", required=True, help="Output JSON path.")
    sub.add_argument("--limit", type=int, default=10)
    add_json_flag(sub)

    sub = subparsers.add_parser("research", help="Research selected candidates through read-only APIs without using an AI.")
    sub.add_argument("--candidates", required=True, help="Candidate or ranked TSV.")
    sub.add_argument("--out-dir", required=True, help="Output directory.")
    sub.add_argument("--species", required=True)
    sub.add_argument("--trait", required=True)
    sub.add_argument("--top-n", type=int, default=10)
    sub.add_argument("--literature-limit", type=int, default=10)
    add_json_flag(sub)


def _run_literature(args: Any) -> dict:
    result = europe_pmc_search(args.query, args.out, args.limit)
    return result


def _run_research(args: Any) -> dict:
    result = candidate_research(args.candidates, args.out_dir, args.species, args.trait, args.top_n, args.literature_limit)
    return result


def _run_phenotype_investigate(args: Any) -> dict:
    return phenotype_investigate(
        args.out_dir, args.case_id, args.focal_species, args.phenotype, args.direction, args.inheritance,
        args.tissue, args.stage, args.onset, args.penetrance, args.synonym, args.ontology_term,
        args.exclusion, args.constraint, args.uncertainty, args.candidates,
        benchmark_mode=args.benchmark_mode,
    )


def _run_mechanism_import(args: Any) -> dict:
    return mechanism_import(
        payload_arg=args.payload, out_dir_arg=args.out_dir, investigation_arg=args.investigation,
        payload_source=args.payload_source, model=args.model, case_id=args.case_id,
        retrieval_provenance=args.retrieval_provenance, candidates_arg=args.candidates,
        min_mechanisms=args.min_mechanisms, proposal_round=args.proposal_round,
        interval_informed=args.interval_informed, mechanism_freeze_manifest=args.mechanism_freeze_manifest,
        localization_freeze_manifest=args.localization_freeze_manifest, benchmark_mode=args.benchmark_mode,
    )


def _run_phenotype_localization_packet(args: Any) -> dict:
    return phenotype_localization_packet(args.investigation, args.out_dir, args.candidates)


def _run_phenotype_localization_freeze(args: Any) -> dict:
    return phenotype_localization_freeze(args.payload, args.investigation, args.out_dir)


def _run_mechanism_coverage_audit(args: Any) -> dict:
    return mechanism_coverage_audit(
        args.mechanisms, args.localization_freeze_manifest, args.investigation,
        args.out_dir, args.candidates, args.gap_request_audit, args.mechanism_import_manifest,
    )


def _run_investigation_attribution(args: Any) -> dict:
    return investigation_attribution(
        args.out_dir, args.mechanisms, args.validation, args.panel, args.interval_match,
        args.ranked, args.case_id, args.min_mechanisms,
    )


def _run_reference_gene_resolve(args: Any) -> dict:
    from .production_research import run_reference_gene_resolve
    return run_reference_gene_resolve(args)


def _run_source_verify(args: Any) -> dict:
    return source_verify(
        args.mechanisms, args.out_dir, args.cache, args.allow_external_api,
        args.id_column, args.source_column, args.full_hash, args.claim_extraction,
        verification_policy=args.verification_policy, focal_species=args.focal_species, reference_bridge_arg=args.reference_bridge,
    )


def _run_claim_verification_benchmark(args: Any) -> dict:
    return claim_verification_benchmark(args.cases, args.out_dir, args.cache,
                                        args.allow_external_api, args.full_hash, args.cache_is_complete)


HANDLERS = {
    "reference-gene-resolve": _run_reference_gene_resolve,
    "source-verify": _run_source_verify,
    "claim-verification-benchmark": _run_claim_verification_benchmark,
    "literature": _run_literature,
    "research": _run_research,
    "phenotype-investigate": _run_phenotype_investigate,
    "mechanism-import": _run_mechanism_import,
    "phenotype-localization-packet": _run_phenotype_localization_packet,
    "phenotype-localization-freeze": _run_phenotype_localization_freeze,
    "mechanism-coverage-audit": _run_mechanism_coverage_audit,
    "investigation-attribution": _run_investigation_attribution,
}
