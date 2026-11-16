"""CLI commands for agent-in-the-loop homolog linking (issue #400 step 2)."""

from __future__ import annotations

from argparse import _SubParsersAction
from typing import Any

from candidate_qtl.homolog_linking import (
    calculate_counterfactual_ranks,
    export_homolog_linking_prompts,
    export_memorization_probe_prompt,
    score_linking_and_probe_results,
    validate_and_import_linking_proposals,
)
from ._shared import add_json_flag


def register(subparsers: _SubParsersAction) -> None:
    sub = subparsers.add_parser(
        "homolog-linking-export",
        help="Export provider-neutral prompt package for agent-in-the-loop homolog linking.",
    )
    sub.add_argument("--out-dir", required=True, help="Output directory for prompt package.")
    sub.add_argument("--case-id", required=True, help="Case identifier.")
    sub.add_argument("--evidence", required=True, help="04_evidence.tsv from discovery.")
    sub.add_argument("--phenotype", default="", help="Phenotype description text.")
    sub.add_argument("--discover-manifest", default=None, help="discover_manifest.json from discovery.")
    add_json_flag(sub)

    sub = subparsers.add_parser(
        "homolog-linking-import",
        help="Mechanically validate and import agent linking proposals with strict closed-world and grounding substring gates.",
    )
    sub.add_argument("--prompt", required=True, help="homolog_linking_prompt.json or export manifest.")
    sub.add_argument("--response", required=True, help="Agent linking response JSON file.")
    sub.add_argument("--out-dir", required=True, help="Output directory for proposal artifacts.")
    sub.add_argument("--reviewed-by", required=True, help="Reviewer identity (separate from model proposer).")
    sub.add_argument("--payload-source", default="agent", help="Proposer category (default: agent).")
    sub.add_argument("--provider", default="subagent", help="Provider / engine identifier.")
    sub.add_argument("--model", default="agent-in-the-loop", help="Model name / version.")
    sub.add_argument("--max-proposals", type=int, default=3, help="Maximum number of candidate proposals allowed (budget).")
    add_json_flag(sub)

    sub = subparsers.add_parser(
        "homolog-linking-counterfactual",
        help="Evaluate what truth Stage-A rank would be if linking proposals acted as support, without altering actual rank.",
    )
    sub.add_argument("--ranked", required=True, help="ranked_candidates.tsv from discovery.")
    sub.add_argument("--proposals", required=True, help="homolog_linking_proposals.tsv.")
    sub.add_argument("--truth-genes", required=True, help="Semicolon-separated truth gene IDs.")
    sub.add_argument("--out-tsv", default=None, help="Optional output path for counterfactual ranking TSV.")
    sub.add_argument("--bonus", type=float, default=0.5, help="Counterfactual support score bonus (default: 0.5).")
    add_json_flag(sub)

    sub = subparsers.add_parser(
        "homolog-linking-probe-export",
        help="Export control memorization probe prompt with phenotype ONLY (no interval, no candidates).",
    )
    sub.add_argument("--out-dir", required=True, help="Output directory.")
    sub.add_argument("--case-id", required=True, help="Case identifier.")
    sub.add_argument("--phenotype", required=True, help="Phenotype description text.")
    add_json_flag(sub)

    sub = subparsers.add_parser(
        "homolog-linking-score",
        help="Score linking results stratified by memorization probe. Probe match is three-valued (match / mismatch / not_evaluated_probe_ambiguous): recall_positive only on a shared family_labels alias/symbol key, recall_negative only on explicit probe abstention, and un-catalogued non-matching relations are excluded from both (memory is never reported as reasoning). The alias used and its provenance are recorded in the memorization_probe_score.tsv artifact.",
    )
    sub.add_argument("--linking-manifest", required=True, help="homolog_linking_manifest.json.")
    sub.add_argument("--probe-response", default=None, help="Memorization probe response JSON (optional).")
    sub.add_argument("--truth-genes", required=True, help="Semicolon-separated truth gene IDs.")
    sub.add_argument("--truth-family", required=True, help="Truth family or gene symbol for probe matching.")
    add_json_flag(sub)


def _run_homolog_linking_export(args: Any) -> dict[str, Any]:
    return export_homolog_linking_prompts(
        out_dir=args.out_dir,
        case_id=args.case_id,
        phenotype=args.phenotype,
        evidence_path=args.evidence,
        discover_manifest_path=args.discover_manifest,
    )


def _run_homolog_linking_import(args: Any) -> dict[str, Any]:
    return validate_and_import_linking_proposals(
        prompt_export_manifest_or_json=args.prompt,
        response_payload_path=args.response,
        out_dir=args.out_dir,
        reviewed_by=args.reviewed_by,
        payload_source=args.payload_source,
        provider=args.provider,
        model=args.model,
        max_proposals=args.max_proposals,
    )


def _run_homolog_linking_counterfactual(args: Any) -> dict[str, Any]:
    return calculate_counterfactual_ranks(
        ranked_tsv_path=args.ranked,
        linking_proposals_tsv_path=args.proposals,
        truth_gene_ids=args.truth_genes,
        out_tsv_path=args.out_tsv,
        counterfactual_bonus=args.bonus,
    )


def _run_homolog_linking_probe_export(args: Any) -> dict[str, Any]:
    return export_memorization_probe_prompt(
        out_dir=args.out_dir,
        case_id=args.case_id,
        phenotype=args.phenotype,
    )


def _run_homolog_linking_score(args: Any) -> dict[str, Any]:
    return score_linking_and_probe_results(
        linking_manifest_path=args.linking_manifest,
        memorization_probe_response_path=args.probe_response,
        truth_gene_ids=args.truth_genes,
        truth_family_or_symbol=args.truth_family,
    )


HANDLERS = {
    "homolog-linking-export": _run_homolog_linking_export,
    "homolog-linking-import": _run_homolog_linking_import,
    "homolog-linking-counterfactual": _run_homolog_linking_counterfactual,
    "homolog-linking-probe-export": _run_homolog_linking_probe_export,
    "homolog-linking-score": _run_homolog_linking_score,
}
