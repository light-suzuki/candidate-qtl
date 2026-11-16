"""CLI commands for the llm feature group (issue #144)."""

from __future__ import annotations

from ..hardening import llm_evidence_import, research_packet
from ..llm import benchmark_shortlist_reviews, export_llm_prompts, export_shortlist_review_prompt, import_shortlist_review, run_llm_prompts
from ..llm_pathway import (
    EXECUTION_MODES, bait_panel_freeze, hypothesis_verify, llm_pathway_import,
    llm_pathway_prompt, phenotype_pathway_benchmark,
)
from argparse import _SubParsersAction
from typing import Any
from ._shared import add_json_flag


def register(subparsers: _SubParsersAction) -> None:
    sub = subparsers.add_parser("research-packet", help="Create AI-optional candidate research packets and prompts.")
    sub.add_argument("--candidates", required=True, help="Candidate or ranked TSV.")
    sub.add_argument("--out-dir", required=True, help="Output directory.")
    sub.add_argument("--species", required=True)
    sub.add_argument("--trait", required=True)
    sub.add_argument("--top-n", type=int, default=10)
    sub.add_argument("--fetch-literature", action="store_true", help="Read Europe PMC for each selected candidate and cache JSON results.")
    sub.add_argument("--literature-limit", type=int, default=10)
    sub.add_argument("--evidence", help="Canonical evidence TSV to embed.")
    sub.add_argument("--contributions", help="Score contributions TSV to embed.")
    add_json_flag(sub)

    sub = subparsers.add_parser("llm-prompt-export", help="Export provider-neutral, reproducible candidate research prompts as JSONL.")
    sub.add_argument("--packet", required=True, help="research_packet.json")
    sub.add_argument("--out", required=True, help="Prompt JSONL output.")
    add_json_flag(sub)

    sub = subparsers.add_parser("llm-shortlist-prompt", help="Export one optional Stage-B prompt that ranks the exact frozen Stage-A set.")
    sub.add_argument("--packet", required=True, help="research_packet.json for the deterministic shortlist.")
    sub.add_argument("--out", required=True, help="Output prompt JSONL.")
    add_json_flag(sub)

    sub = subparsers.add_parser("llm-shortlist-import", help="Validate an LLM batch ranking without changing Stage A or scoring evidence.")
    sub.add_argument("--shortlist", required=True)
    sub.add_argument("--responses", required=True, help="llm_candidate_responses.json from llm-run.")
    sub.add_argument("--out-dir", required=True)
    add_json_flag(sub)

    sub = subparsers.add_parser("llm-review-benchmark", help="Compare deterministic and optional LLM review ranks on identical frozen shortlists.")
    sub.add_argument("--config", required=True, help="TSV with case, truth, shortlist, and two review ranking paths.")
    sub.add_argument("--out-dir", required=True)
    add_json_flag(sub)

    sub = subparsers.add_parser("llm-run", help="Run prompt JSONL through an API/LM Studio provider or normalize an offline agent artifact.")
    sub.add_argument("--prompts", required=True)
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--provider", required=True, choices=["openai-responses", "openai-compatible", "lmstudio-cli", "codex-agent"])
    sub.add_argument("--model", required=True)
    sub.add_argument("--endpoint", help="HTTP endpoint, or lms executable path for lmstudio-cli.")
    sub.add_argument("--agent-file", help="Offline Codex/OpenCode JSON/JSONL artifact; use with --provider codex-agent.")
    sub.add_argument("--api-key-env", default="OPENAI_API_KEY")
    sub.add_argument("--allow-external-api", action="store_true")
    sub.add_argument("--max-output-tokens", type=int, default=2500)
    sub.add_argument("--timeout", type=int, default=180)
    add_json_flag(sub)

    sub = subparsers.add_parser("llm-evidence-import", help="Import only human-reviewed, source-backed model claims. Each claim needs reviewed=true plus reviewed_by; unresolved source IDs are marked and cannot count as independent.")
    sub.add_argument("--candidates", required=True)
    sub.add_argument("--input", required=True)
    sub.add_argument("--out", required=True)
    sub.add_argument(
        "--resolve-against",
        nargs="+",
        metavar="PATH",
        help=(
            "Local corpus files (reference proteome FASTA, annotation TSV, citation list) a source_id must be found in. "
            "Offline and format-agnostic. Without this, accepted rows are marked unverified and share one independence "
            "group, so a well-formed but non-existent accession cannot reach independent_experimental_support."
        ),
    )
    add_json_flag(sub)

    sub = subparsers.add_parser("llm-pathway-prompt", help="Export a versioned phenotype-only prompt package asking for competing pathway/gene-family hypotheses and literature-supported reference baits. The model never sees or ranks QTL interval candidates.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--case-id", required=True)
    sub.add_argument("--focal-species", required=True)
    sub.add_argument("--phenotype", required=True)
    sub.add_argument("--direction", choices=["increase", "decrease", "gain", "loss", "unknown"], default="unknown")
    sub.add_argument("--inheritance", choices=["recessive", "dominant", "semidominant", "unknown"], default="unknown")
    sub.add_argument("--organ", default="", help="Affected organ or tissue.")
    sub.add_argument("--stage", default="", help="Developmental stage.")
    sub.add_argument("--mode", choices=list(EXECUTION_MODES), default="blind", help="blind: phenotype only. assisted: researcher context allowed, recorded separately in the manifest.")
    sub.add_argument("--constraint", action="append", dest="constraints", help="Researcher constraint; repeatable (assisted or blind).")
    sub.add_argument("--exclusion", action="append", dest="exclusions", help="Known exclusion; repeatable.")
    sub.add_argument("--assisted-context", help="Assisted mode only: file with researcher hypotheses / verified literature / ontology terms.")
    sub.add_argument("--candidates", help="Interval candidates TSV used ONLY to audit that no candidate id leaked into a blind prompt.")
    add_json_flag(sub)

    sub = subparsers.add_parser("llm-pathway-import", help="Normalize raw model output into the strict pathway-proposal schema, keeping the raw payload verbatim. Every proposal stays pending_verification.")
    sub.add_argument("--raw", required=True, help="Raw model output (.json/.jsonl/.txt), or a manually written proposal file.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--input-manifest", help="llm_input_manifest.json from llm-pathway-prompt.")
    sub.add_argument("--model", default="", help="Model/provider version that produced the output.")
    sub.add_argument("--proposal-source", default="", help="Proposal source, e.g. a provider name or `manual`.")
    sub.add_argument("--case-id", default="")
    sub.add_argument("--execution-mode", choices=list(EXECUTION_MODES), default="")
    add_json_flag(sub)

    sub = subparsers.add_parser("hypothesis-verify", help="Verify proposed accessions against an authoritative database: existence, species match, retrievable sequence, source-backed function. Unverified proposals never become baits.")
    sub.add_argument("--hypotheses", required=True, help="llm_pathway_hypotheses.tsv")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--focal-species", required=True, help="Focal species; a reference gene from it is rejected.")
    sub.add_argument("--record-cache", help="Local UniProt-shaped JSON record file or directory for offline verification.")
    sub.add_argument("--source-records", help="verified_sources.json from source-verify; required to validate publication annotations in canonical Route A.")
    sub.add_argument("--bound-claims", help="Shared bound-claim TSV joining citation, gene, experiment and phenotype. Without it, source/cue evidence remains review-only.")
    sub.add_argument("--agent-claim-review", help="Independent agent review packet with frozen raw response, prompt and separate executor provenance; alternative to --bound-claims.")
    sub.add_argument("--source-review-records", help="source_review_records.json from source-verify; required with --agent-claim-review for exact text/ID/species checks.")
    sub.add_argument("--allow-external-api", action="store_true", help="Allow accession lookups against the UniProt REST API.")
    sub.add_argument("--benchmark-mode", action="store_true", help="Additionally reject focal-lineage / answer-backmapped markers.")
    sub.add_argument("--min-panel-size", type=int, default=2, help="Panels below this size are review-required, not rejected.")
    add_json_flag(sub)

    sub = subparsers.add_parser("bait-panel-freeze", help="Freeze verified proposals into an immutable bait panel, FASTA, and a query manifest consumable by query-eligibility.")
    sub.add_argument("--hypotheses", required=True)
    sub.add_argument("--validation", required=True, help="proposal_validation.tsv from hypothesis-verify.")
    sub.add_argument("--sequences", required=True, help="verified_sequences.json from hypothesis-verify.")
    sub.add_argument("--reviewable-sequences", help="Inclusive identity/species/sequence-confirmed sequence store from hypothesis-verify, for review-only projection.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--case-id", default="")
    sub.add_argument("--min-panel-size", type=int, default=2)
    add_json_flag(sub)

    sub = subparsers.add_parser("phenotype-pathway-benchmark", help="Score the phenotype-only proposal lane stage by stage (pathway/family recall, accession existence, verification pass, downstream causal recovery). Never pooled with the map-based lane.")
    sub.add_argument("--config", required=True, help="TSV: case_id/expected_pathway/expected_gene_family[/causal_gene_id/hypotheses/validation/panel/ranked].")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--top-n", type=int, default=20)
    add_json_flag(sub)


def _run_research_packet(args: Any) -> dict:
    result = research_packet(
        args.candidates,
        args.out_dir,
        args.species,
        args.trait,
        args.top_n,
        args.fetch_literature,
        args.literature_limit,
        args.evidence,
        args.contributions,
    )
    return result


def _run_llm_prompt_export(args: Any) -> dict:
    result = export_llm_prompts(args.packet, args.out)
    return result


def _run_llm_shortlist_prompt(args: Any) -> dict:
    result = export_shortlist_review_prompt(args.packet, args.out)
    return result


def _run_llm_shortlist_import(args: Any) -> dict:
    result = import_shortlist_review(args.shortlist, args.responses, args.out_dir)
    return result


def _run_llm_review_benchmark(args: Any) -> dict:
    result = benchmark_shortlist_reviews(args.config, args.out_dir)
    return result


def _run_llm_run(args: Any) -> dict:
    result = run_llm_prompts(
        args.prompts,
        args.out_dir,
        args.provider,
        args.model,
        args.endpoint,
        args.api_key_env,
        args.allow_external_api,
        args.max_output_tokens,
        args.timeout,
        args.agent_file,
    )
    return result


def _run_llm_evidence_import(args: Any) -> dict:
    result = llm_evidence_import(args.candidates, args.input, args.out, args.resolve_against)
    return result


def _run_llm_pathway_prompt(args: Any) -> dict:
    return llm_pathway_prompt(
        args.out_dir, args.case_id, args.focal_species, args.phenotype, args.direction, args.inheritance,
        args.organ, args.stage, args.mode, args.constraints, args.exclusions, args.assisted_context, args.candidates,
    )


def _run_llm_pathway_import(args: Any) -> dict:
    return llm_pathway_import(
        args.raw, args.out_dir, args.input_manifest, args.model, args.proposal_source,
        args.case_id, args.execution_mode,
    )


def _run_hypothesis_verify(args: Any) -> dict:
    return hypothesis_verify(
        args.hypotheses, args.out_dir, args.focal_species, args.record_cache,
        args.allow_external_api, args.benchmark_mode, args.min_panel_size,
        source_records_arg=args.source_records,
        bound_claims_arg=args.bound_claims,
        agent_claim_review_arg=args.agent_claim_review,
        source_review_records_arg=args.source_review_records,
    )


def _run_bait_panel_freeze(args: Any) -> dict:
    return bait_panel_freeze(args.hypotheses, args.validation, args.sequences, args.out_dir, args.case_id, args.min_panel_size,
                             reviewable_sequences_arg=args.reviewable_sequences)


def _run_phenotype_pathway_benchmark(args: Any) -> dict:
    return phenotype_pathway_benchmark(args.config, args.out_dir, args.top_n)


HANDLERS = {
    "research-packet": _run_research_packet,
    "llm-pathway-prompt": _run_llm_pathway_prompt,
    "llm-pathway-import": _run_llm_pathway_import,
    "hypothesis-verify": _run_hypothesis_verify,
    "bait-panel-freeze": _run_bait_panel_freeze,
    "phenotype-pathway-benchmark": _run_phenotype_pathway_benchmark,
    "llm-prompt-export": _run_llm_prompt_export,
    "llm-shortlist-prompt": _run_llm_shortlist_prompt,
    "llm-shortlist-import": _run_llm_shortlist_import,
    "llm-review-benchmark": _run_llm_review_benchmark,
    "llm-run": _run_llm_run,
    "llm-evidence-import": _run_llm_evidence_import,
}
