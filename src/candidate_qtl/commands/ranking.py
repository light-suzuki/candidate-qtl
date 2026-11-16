"""CLI commands for the ranking feature group (issue #144)."""

from __future__ import annotations

from ..core import dry_prioritize, rank_ledger
from ..evidence_tiers import evidence_tier_budget_benchmark, evidence_tier_report
from ..hardening import candidate_pipeline, prioritize_evidence
from ..mapping_gate import mapping_adequacy_gate
from ..next_experiment import select_next_experiment
from ..ranking_qa import ranking_qa, ranking_stability_benchmark
from ..shortlist import screening_shortlist
from argparse import _SubParsersAction
from typing import Any
from ._shared import add_json_flag


def register(subparsers: _SubParsersAction) -> None:
    sub = subparsers.add_parser("rank", help="Rank a candidate ledger with deterministic, inspectable rules.")
    sub.add_argument("--ledger", required=True, help="TSV based on templates/candidate_ledger.tsv.")
    sub.add_argument("--out", required=True, help="Output ranked TSV.")
    sub.add_argument("--research-summary", help="Optional research_summary.tsv from candidate-qtl research.")
    sub.add_argument("--go-evidence", help="Optional gene-level GO evidence TSV from annotate-plan.")
    add_json_flag(sub)

    sub = subparsers.add_parser("dry-prioritize", help="Rank parental-sequencing candidates using pre-Wet evidence only.")
    sub.add_argument("--candidates", required=True, help="TSV based on templates/dry_candidate_prioritization.tsv.")
    sub.add_argument("--out", required=True, help="Output TSV with Dry score, group, rank, and score breakdown.")
    sub.add_argument("--top-n", type=int, default=7, help="Number of non-rescue candidates in the first sequencing batch.")
    sub.add_argument("--reserve-n", type=int, default=3, help="Number of non-rescue candidates retained as the next batch.")
    add_json_flag(sub)

    sub = subparsers.add_parser("evidence-rank", help="Rank candidates from the canonical evidence ledger with capped, inspectable contributions.")
    sub.add_argument("--candidates", required=True)
    sub.add_argument("--evidence", required=True)
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--weights", help="Optional TSV with evidence_code and weight.")
    sub.add_argument("--top-n", type=int, default=7)
    sub.add_argument("--reserve-n", type=int, default=3)
    sub.add_argument("--rescue-n", type=int, default=3, help="Maximum Rescue-group candidates retained.")
    sub.add_argument("--min-top-score", type=float, default=2.0)
    sub.add_argument("--case-id", default="", help="Optional case identifier recorded in the decision artifact.")
    add_json_flag(sub)

    sub = subparsers.add_parser("ranking-qa", help="Audit ranking stability across eligibility, family dropout, weight perturbation, confidence, provider, source-group dropout, and provider-dropout paths.")
    sub.add_argument("--candidates", required=True)
    sub.add_argument("--evidence", required=True)
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--weights")
    sub.add_argument("--top-k", type=int, default=5)
    sub.add_argument("--perturbation", action="append", type=float, help="Fractional weight perturbation; repeatable (defaults: 0.10, 0.25, 0.50).")
    sub.add_argument("--source-dropout", action="append", help="Resolved independence group to drop (e.g. comparative_annotation_chain, PMID:12345); repeatable. Defaults to every non-empty group found in the evidence.")
    sub.add_argument("--provider-dropout", action="append", help="source_type to drop (e.g. blastp, eggnog, qtl_position); repeatable. Defaults to every non-empty source_type found in the evidence.")
    sub.add_argument("--benchmark-mode", action="store_true", help="Mark this run as calibration/benchmark rather than a normal case.")
    add_json_flag(sub)

    sub = subparsers.add_parser("ranking-stability-benchmark", help="Report whether frozen-case truth genes stay within Top-K under weight perturbation and source/provider dropout.")
    sub.add_argument("--manifest", required=True, help="JSON manifest with cases: case_id, candidates, evidence, ranked_tsv, answer_gene_id.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--top-k", type=int, default=20)
    add_json_flag(sub)

    sub = subparsers.add_parser("screening-shortlist", help="Freeze the non-LLM Stage-A screening shortlist before human or LLM review.")
    sub.add_argument("--ranked", required=True, help="ranked_candidates.tsv from evidence-rank or candidate-pipeline.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--target-size", type=int, default=20, help="Primary deterministic review budget.")
    sub.add_argument("--rescue-size", type=int, default=5, help="Separate Rescue appendix budget; not counted in the primary shortlist.")
    add_json_flag(sub)

    sub = subparsers.add_parser("evidence-tiers", help="Issue #229: report target_size-independent biological evidence tiers (strong/plausible/unresolved/contradictory/knowledge_limited) side by side with the configurable review_budget/assay_capacity operational cut. The operational cut never implies a biological exclusion.")
    sub.add_argument("--ranked", required=True, help="ranked_candidates.tsv from evidence-rank or candidate-pipeline (needs top_eligibility).")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--review-budget", type=int, default=20, help="Assay-capacity cut (same meaning/default as screening-shortlist's --target-size); changes operational_state only, never evidence_tier.")
    sub.add_argument("--reserve-size", type=int, default=0, help="Extra rank window right after review_budget marked operational_state=reserve; 0 disables it.")
    sub.add_argument("--knowledge-coverage", help="Optional knowledge_coverage.tsv from knowledge-coverage; enables the knowledge_limited tier.")
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)

    sub = subparsers.add_parser("evidence-tier-benchmark", help="Issue #229: Recall@budget alongside evidence-tier recall, stratified by interval-size band via knowledge_coverage.benchmark_stratify. Reads frozen evidence-rank outputs only; never re-ranks.")
    sub.add_argument("--manifest", required=True, help="TSV: one row per case (case_id/ranked_tsv[/answer_gene_id/knowledge_coverage_tsv]).")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--review-budget", type=int, default=20)
    sub.add_argument("--full-hash", action="store_true")
    add_json_flag(sub)

    sub = subparsers.add_parser("weight-policy-benchmark", help="Compare numeric weight baselines and structural ranking gates on a frozen known-gene benchmark.")
    sub.add_argument("--manifest", required=True, help="benchmark_manifest.json from benchmark-reference.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--blocked-term", action="append", default=[], help="Known answer/locus term excluded from upstream benchmark vocabulary; repeatable and recorded.")
    sub.add_argument("--bootstrap-iterations", type=int, default=1000)
    sub.add_argument("--seed", type=int, default=20260718)
    sub.add_argument("--crop", default="Pisum sativum", help="Crop label when cases do not carry one.")
    add_json_flag(sub)

    sub = subparsers.add_parser("mapping-gate", help="Assess whether QTL mapping is adequate before candidate triage.")
    sub.add_argument("--metrics", required=True, help="TSV with metric, value, and source columns.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--thresholds", help="Optional TSV with threshold and value columns.")
    add_json_flag(sub)

    sub = subparsers.add_parser("next-experiment", help="Rank next experimental actions separately from biological candidate rank.")
    sub.add_argument("--ranked", required=True, help="Deterministic ranked candidate TSV; read-only biological priority input.")
    sub.add_argument("--case", required=True, help="Case manifest with decision_support constraints.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--mapping-gate")
    sub.add_argument("--candidate-path-decision")
    sub.add_argument("--candidate-flags")
    sub.add_argument("--assay-outcomes")
    sub.add_argument("--exon-coverage")
    sub.add_argument("--informative-recombinants")
    sub.add_argument("--wet-candidate-state", help="Normalized candidate_wet_state.tsv from wet-state.")
    sub.add_argument("--hypotheses", help="Locus hypothesis ledger; read-only action context, never biological rank input.")
    sub.add_argument("--remaining-candidates",
                    help="Canonical remaining_candidate_set.tsv from candidate-elimination; "
                         "post-elimination candidate counts are used when supplied, never ranked reordering.")
    add_json_flag(sub)

    sub = subparsers.add_parser("candidate-pipeline", help="Run consistent primary evidence/ranking and optionally plan or integrate Top-candidate orthology refinement.")
    sub.add_argument("--candidates", required=True)
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--blast", action="append")
    sub.add_argument("--diamond", action="append")
    sub.add_argument("--hmmer", action="append")
    sub.add_argument("--eggnog", action="append")
    sub.add_argument("--expression")
    sub.add_argument("--target-tissue", action="append")
    sub.add_argument("--target-stage", action="append")
    sub.add_argument("--min-expression", type=float, default=1.0)
    sub.add_argument("--target-fold", type=float, default=2.0)
    sub.add_argument("--manual-evidence")
    sub.add_argument("--refined-evidence", help="Canonical evidence produced by orthology-import.")
    sub.add_argument("--peak")
    sub.add_argument("--interval")
    sub.add_argument("--peak-window", type=int, default=0)
    sub.add_argument("--boundary-window", type=int, default=0)
    sub.add_argument("--weights")
    sub.add_argument("--top-n", type=int, default=7)
    sub.add_argument("--reserve-n", type=int, default=3)
    sub.add_argument("--rescue-n", type=int, default=3)
    sub.add_argument("--min-top-score", type=float, default=2.0)
    sub.add_argument("--protein-fasta", help="Candidate proteins normalized to gene_id.")
    sub.add_argument("--target-proteome")
    sub.add_argument("--relative-proteome", action="append", help="LABEL=FASTA; repeatable.")
    sub.add_argument("--num-threads", type=int, default=4)
    add_json_flag(sub)


def _run_rank(args: Any) -> dict:
    result = rank_ledger(args.ledger, args.out, args.research_summary, args.go_evidence)
    return result


def _run_dry_prioritize(args: Any) -> dict:
    result = dry_prioritize(args.candidates, args.out, args.top_n, args.reserve_n)
    return result


def _run_evidence_rank(args: Any) -> dict:
    result = prioritize_evidence(
        args.candidates,
        args.evidence,
        args.out_dir,
        args.weights,
        args.top_n,
        args.reserve_n,
        args.min_top_score,
        args.rescue_n,
        args.case_id,
    )
    return result


def _run_ranking_qa(args: Any) -> dict:
    result = ranking_qa(args.candidates, args.evidence, args.out_dir, args.weights, args.top_k, args.perturbation,
                        args.benchmark_mode, args.source_dropout, args.provider_dropout)
    return result


def _run_ranking_stability_benchmark(args: Any) -> dict:
    result = ranking_stability_benchmark(args.manifest, args.out_dir, args.top_k)
    return result


def _run_screening_shortlist(args: Any) -> dict:
    result = screening_shortlist(args.ranked, args.out_dir, args.target_size, args.rescue_size)
    return result


def _run_evidence_tiers(args: Any) -> dict:
    result = evidence_tier_report(
        args.ranked,
        args.out_dir,
        args.review_budget,
        args.reserve_size,
        args.knowledge_coverage,
        args.case_id,
    )
    return result


def _run_evidence_tier_benchmark(args: Any) -> dict:
    result = evidence_tier_budget_benchmark(args.manifest, args.out_dir, args.review_budget, args.full_hash)
    return result


def _run_weight_policy_benchmark(args: Any) -> dict:
    from ..weight_policy_benchmark import weight_policy_benchmark
    result = weight_policy_benchmark(
        args.manifest,
        args.out_dir,
        args.blocked_term,
        args.bootstrap_iterations,
        args.seed,
        args.crop,
    )
    return result


def _run_mapping_gate(args: Any) -> dict:
    result = mapping_adequacy_gate(args.metrics, args.out_dir, args.thresholds)
    return result


def _run_next_experiment(args: Any) -> dict:
    result = select_next_experiment(
        args.ranked,
        args.case,
        args.out_dir,
        mapping_gate_arg=args.mapping_gate,
        candidate_path_decision_arg=args.candidate_path_decision,
        candidate_flags_arg=args.candidate_flags,
        assay_outcomes_arg=args.assay_outcomes,
        exon_coverage_arg=args.exon_coverage,
        informative_recombinants_arg=args.informative_recombinants,
        wet_candidate_state_arg=args.wet_candidate_state,
        hypotheses_arg=args.hypotheses,
        remaining_candidates_arg=args.remaining_candidates,
    )
    return result


def _run_candidate_pipeline(args: Any) -> dict:
    result = candidate_pipeline(
        args.candidates,
        args.out_dir,
        args.blast,
        args.diamond,
        args.hmmer,
        args.eggnog,
        args.expression,
        args.manual_evidence,
        args.refined_evidence,
        args.target_tissue,
        args.target_stage,
        args.min_expression,
        args.target_fold,
        args.peak,
        args.interval,
        args.peak_window,
        args.boundary_window,
        args.weights,
        args.top_n,
        args.reserve_n,
        args.min_top_score,
        args.rescue_n,
        args.protein_fasta,
        args.target_proteome,
        args.relative_proteome,
        args.num_threads,
    )
    return result


HANDLERS = {
    "rank": _run_rank,
    "dry-prioritize": _run_dry_prioritize,
    "evidence-rank": _run_evidence_rank,
    "ranking-qa": _run_ranking_qa,
    "ranking-stability-benchmark": _run_ranking_stability_benchmark,
    "screening-shortlist": _run_screening_shortlist,
    "evidence-tiers": _run_evidence_tiers,
    "evidence-tier-benchmark": _run_evidence_tier_benchmark,
    "weight-policy-benchmark": _run_weight_policy_benchmark,
    "mapping-gate": _run_mapping_gate,
    "next-experiment": _run_next_experiment,
    "candidate-pipeline": _run_candidate_pipeline,
}
