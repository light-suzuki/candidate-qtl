"""CLI commands for the evidence feature group (issue #144)."""

from __future__ import annotations

from ..controlled_terms import promote_controlled_terms
from ..elimination import candidate_elimination
from ..evidence_policy import evidence_build
from ..peak_intervals import peak_context_import
from ..qtl_source import qtl_source_init
from ..trait_ontology import map_trait_terms
from ..trait_terms import build_trait_terms
from ..workflow_extensions import direction_evidence
from argparse import _SubParsersAction
from typing import Any
from ._shared import add_json_flag


def register(subparsers: _SubParsersAction) -> None:
    sub = subparsers.add_parser("candidate-elimination", help="Manage candidate and mechanism survival state. Only hard positional/genetic evidence may exclude a gene; mechanism-specific negatives reject one mechanism; assay failure, partial coverage and unmeasured data never eliminate anything.")
    sub.add_argument("--candidates", required=True, help="Candidate/ranked/shortlist TSV with gene_id.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--evidence", help="Elimination evidence TSV (gene_id/evidence_code[/hypothesis_id/evidence_id/quality/coverage]).")
    sub.add_argument("--hypotheses", help="locus_hypotheses.tsv, for mechanism-level state.")
    sub.add_argument("--wet-state", help="candidate_wet_state.tsv, for coverage and quality context.")
    sub.add_argument("--rules", help="Override the versioned elimination rule registry (templates/elimination_rules.tsv).")
    sub.add_argument("--previous-history", help="elimination_events.tsv from an earlier run; history is append-only.")
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)

    sub = subparsers.add_parser("direction-evidence", help="Mark phenotype-direction matches and counterevidence deterministically.")
    sub.add_argument("--evidence", required=True)
    sub.add_argument("--expected", required=True, choices=["increase", "decrease", "gain", "loss"])
    sub.add_argument("--out", required=True)
    add_json_flag(sub)

    sub = subparsers.add_parser("evidence-build", help="Normalize annotation, expression, position, and manual evidence into one source-backed ledger.")
    sub.add_argument("--candidates", required=True, help="Transcript-aware candidates.tsv from extract.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--blast", action="append", help="BLASTP outfmt 6 file; repeatable.")
    sub.add_argument("--diamond", action="append", help="DIAMOND outfmt 6 file; repeatable.")
    sub.add_argument("--hmmer", action="append", help="HMMER --tblout file; repeatable.")
    sub.add_argument("--eggnog", action="append", help="eggNOG-mapper annotations file; repeatable.")
    sub.add_argument("--expression", help="Long-form TSV with gene_id, tissue, value, and optional is_target.")
    sub.add_argument("--target-tissue", action="append", help="Target tissue name; repeatable unless expression TSV has is_target.")
    sub.add_argument("--target-stage", action="append", help="Target stage; combined with --target-tissue unless is_target is supplied.")
    sub.add_argument("--target-condition", action="append", help="Restrict target expression rows to these condition values; repeatable.")
    sub.add_argument("--min-expression", type=float, default=1.0)
    sub.add_argument("--target-fold", type=float, default=2.0)
    sub.add_argument("--manual-evidence", help="Source-backed TSV using the canonical evidence schema.")
    sub.add_argument("--existing-evidence", help="Existing canonical evidence TSV, for example orthology refinement results.")
    sub.add_argument("--peak", help="QTL peak as CHROM:POS or POS.")
    sub.add_argument("--interval", help="QTL interval; inferred from candidate table when unique.")
    sub.add_argument(
        "--peak-window",
        type=int,
        default=0,
        help="Analyst-supplied peak uncertainty half-width for scored peak_proximal support; 0 scores exact overlap only and reports the 5%%-interval heuristic as unscored context.",
    )
    sub.add_argument("--boundary-window", type=int, default=0)
    add_json_flag(sub)

    sub = subparsers.add_parser("evidence-term-promote", help="Promote approved annotation terms into low-confidence deterministic screening evidence.")
    sub.add_argument("--evidence", required=True)
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--term", action="append", required=True, help="Approved exact token/phrase; repeatable.")
    sub.add_argument("--demote-unapproved-existing", action="store_true", help="Convert pre-existing trait matches outside the approved vocabulary to unscored context.")
    add_json_flag(sub)

    sub = subparsers.add_parser("peak-context", help="Import QTL peak, LOD support interval, and recombinant core into auditable position evidence.")
    sub.add_argument("--candidates", help="Candidates.tsv with gene_id/chrom/start/end; position evidence is emitted per gene when supplied.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--lod-profile", help="TSV with chrom, position, and lod columns for a LOD-drop support interval.")
    sub.add_argument("--lod-drop", type=float, default=1.5, help="LOD drop for the support interval (default 1.5).")
    sub.add_argument("--peak-metadata", help="Single-QTL TSV following templates/qtl_peak_metadata.tsv.")
    sub.add_argument("--peak", help="QTL peak as CHROM:POS or POS when no profile/metadata is available.")
    sub.add_argument("--interval", help="QTL interval CHROM:START-END; used only as context, never as positive evidence.")
    sub.add_argument("--peak-window", type=int, default=0, help="Analyst-supplied uncertainty half-width for peak_proximal support.")
    add_json_flag(sub)

    sub = subparsers.add_parser("trait-terms", help="Build a controlled, approval-gated trait-term set for keyword evidence.")
    trait_sub = sub.add_subparsers(dest="trait_command", required=True)
    build = trait_sub.add_parser("build", help="Classify, approve, and export trait terms; reject answer-leakage in benchmark mode.")
    build.add_argument("--trait-terms", required=True, help="Raw TSV following templates/trait_terms.tsv.")
    build.add_argument("--out", required=True, help="Approved trait-term TSV consumed by annotation-keywords.")
    build.add_argument("--ontology-map", help="Optional local ontology table to attach ontology IDs.")
    build.add_argument("--manual-synonyms", help="Optional human synonym TSV with term and synonym columns.")
    build.add_argument("--llm-proposals", help="Optional JSON list of model-proposed synonyms; never auto-approved.")
    build.add_argument("--blocked-term", action="append", help="Known causal gene/locus name to reject as answer leakage; repeatable.")
    build.add_argument("--benchmark", action="store_true", help="Reject terms matching blocked answer names.")
    add_json_flag(build)

    sub = subparsers.add_parser("trait-ontology-map", help="Map free trait terms to ontology IDs using a local synonym table only.")
    sub.add_argument("--trait-terms", required=True, help="TSV following templates/trait_terms.tsv.")
    sub.add_argument("--ontology-map", required=True, help="Local TSV with term, ontology_id, and optional ontology_name/match_type/confidence.")
    sub.add_argument("--out", required=True, help="Normalized trait-ontology TSV.")
    add_json_flag(sub)

    sub = subparsers.add_parser("qtl-source-init", help="Validate coordinate provenance and choose a safe liftover route.")
    sub.add_argument("--manifest", required=True, help="YAML/JSON QTL source manifest.")
    sub.add_argument("--out-dir", required=True)
    add_json_flag(sub)


def _run_direction_evidence(args: Any) -> dict:
    result = direction_evidence(args.evidence, args.expected, args.out)
    return result


def _run_evidence_build(args: Any) -> dict:
    result = evidence_build(
        args.candidates,
        args.out_dir,
        args.blast,
        args.diamond,
        args.hmmer,
        args.eggnog,
        args.expression,
        args.manual_evidence,
        args.existing_evidence,
        args.target_tissue,
        args.target_stage,
        args.min_expression,
        args.target_fold,
        args.peak,
        args.interval,
        args.peak_window,
        args.boundary_window,
        args.target_condition,
    )
    return result


def _run_evidence_term_promote(args: Any) -> dict:
    result = promote_controlled_terms(args.evidence, args.out_dir, args.term, args.demote_unapproved_existing)
    return result


def _run_peak_context(args: Any) -> dict:
    result = peak_context_import(
        args.candidates,
        args.out_dir,
        args.lod_profile,
        args.peak_metadata,
        args.peak,
        args.interval,
        args.lod_drop,
        args.peak_window,
    )
    return result


def _run_trait_terms(args: Any) -> dict:
    result = build_trait_terms(
        args.trait_terms,
        args.out,
        args.ontology_map,
        args.manual_synonyms,
        args.llm_proposals,
        args.blocked_term,
        args.benchmark,
    )
    return result


def _run_trait_ontology_map(args: Any) -> dict:
    result = map_trait_terms(args.trait_terms, args.ontology_map, args.out)
    return result


def _run_qtl_source_init(args: Any) -> dict:
    result = qtl_source_init(args.manifest, args.out_dir)
    return result


def _run_candidate_elimination(args: Any) -> dict:
    return candidate_elimination(
        args.candidates, args.out_dir, args.evidence, args.hypotheses, args.wet_state,
        args.rules, args.previous_history, args.case_id,
    )


HANDLERS = {
    "candidate-elimination": _run_candidate_elimination,
    "direction-evidence": _run_direction_evidence,
    "evidence-build": _run_evidence_build,
    "evidence-term-promote": _run_evidence_term_promote,
    "peak-context": _run_peak_context,
    "trait-terms": _run_trait_terms,
    "trait-ontology-map": _run_trait_ontology_map,
    "qtl-source-init": _run_qtl_source_init,
}
