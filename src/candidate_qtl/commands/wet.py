"""CLI commands for the wet feature group (issue #144)."""

from __future__ import annotations

from ..hypothesis_ledger import update_hypothesis_ledger
from ..wet_state import update_wet_state
from argparse import _SubParsersAction
from typing import Any
from ._shared import add_json_flag


def register(subparsers: _SubParsersAction) -> None:
    sub = subparsers.add_parser("wet-state", help="Normalize assay coverage/results and derive append-only candidate Wet states.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--assays", help="Canonical assay registry TSV.")
    sub.add_argument("--coverage", help="Actual Wet target coverage TSV; in-silico design coverage does not pass the Wet gate.")
    sub.add_argument("--results", help="Parent/F1/individual-level Wet result TSV.")
    sub.add_argument("--previous-history", help="Previous wet_state_history.tsv to append without overwriting.")
    sub.add_argument("--primer-results", help="PrimerBlast import TSV; imported as design state only.")
    sub.add_argument("--exon-coverage", help="exon coverage_windows.tsv; imported as design coverage only.")
    sub.add_argument("--min-coverage", type=float, default=1.0)
    sub.add_argument("--required-read-direction", choices=["bidirectional", "forward", "reverse", "any"], default="bidirectional")
    add_json_flag(sub)

    sub = subparsers.add_parser("hypothesis-ledger", help="Track locus/mechanism hypotheses separately from gene rank.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--hypotheses", help="Hypothesis proposal/current-state TSV.")
    sub.add_argument("--case-id", default="")
    sub.add_argument("--evidence", help="Canonical evidence table with observation_id or evidence_id.")
    sub.add_argument("--wet-results", help="Normalized wet_results.tsv with result_id.")
    sub.add_argument("--sv-events", help="Canonical SV table with sv_event_id, call_id, or cluster_id.")
    sub.add_argument("--wet-candidate-state", help="Normalized candidate_wet_state.tsv.")
    sub.add_argument("--annotation-rescue", help="annotation-rescue audit TSV; native_discordant/reference_rescue/ambiguous_overlap rows seed annotation_missing hypotheses.")
    sub.add_argument("--dual-route-comparison", help="dual-route-compare comparison TSV; A/B disagreement seeds linked_qtl or unknown_locus_mechanism hypotheses.")
    sub.add_argument("--previous-history", help="Previous hypothesis_history.tsv for append-only transitions.")
    sub.add_argument("--ranked", help="Optional deterministic rank input; fingerprinted and never modified.")
    add_json_flag(sub)


def _run_wet_state(args: Any) -> dict:
    result = update_wet_state(
        args.out_dir,
        assay_registry_arg=args.assays,
        coverage_arg=args.coverage,
        results_arg=args.results,
        previous_history_arg=args.previous_history,
        primer_results_arg=args.primer_results,
        exon_coverage_arg=args.exon_coverage,
        min_coverage=args.min_coverage,
        required_read_direction=args.required_read_direction,
    )
    return result


def _run_hypothesis_ledger(args: Any) -> dict:
    result = update_hypothesis_ledger(
        args.out_dir,
        hypotheses_arg=args.hypotheses,
        case_id=args.case_id,
        evidence_arg=args.evidence,
        wet_results_arg=args.wet_results,
        sv_events_arg=args.sv_events,
        wet_candidate_state_arg=args.wet_candidate_state,
        annotation_rescue_arg=args.annotation_rescue,
        dual_route_comparison_arg=args.dual_route_comparison,
        previous_history_arg=args.previous_history,
        ranked_arg=args.ranked,
    )
    return result


HANDLERS = {
    "wet-state": _run_wet_state,
    "hypothesis-ledger": _run_hypothesis_ledger,
}
