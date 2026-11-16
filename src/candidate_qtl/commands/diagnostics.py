"""CLI commands for the diagnostics feature group (issue #144)."""

from __future__ import annotations

from ..core import doctor
from ..provenance_ledger import provenance_ledger, render_ledger_report, shared_derivation_paths
from argparse import _SubParsersAction
from typing import Any
from ._shared import add_json_flag


def register(subparsers: _SubParsersAction) -> None:
    sub = subparsers.add_parser("doctor", help="Inspect required and optional WSL tools.")
    add_json_flag(sub)

    sub = subparsers.add_parser(
        "provenance-ledger",
        help="Build and audit the normalized provenance DAG: nodes, edges, and shared derivation paths (issue #231).",
    )
    sub.add_argument("--entries", required=True, help="TSV with node_kind, source_id, label, optional parent_node; optional node_id/source_ref/edge_kind columns.")
    sub.add_argument("--ledger-dir", required=True, help="Directory receiving normalized nodes.tsv / edges.tsv.")
    sub.add_argument("--out-dir", required=True, help="Directory receiving provenance_graph.tsv and the manifest JSON.")
    sub.add_argument("--case-id", default="", help="Case identifier recorded into the manifest.")
    add_json_flag(sub)


def _run_doctor(args: Any) -> dict:
    result = doctor()
    return result


def _run_provenance_ledger(args: Any) -> dict:
    result = provenance_ledger(args.entries, args.ledger_dir, args.out_dir, args.case_id)
    result["shared_derivation_paths"] = shared_derivation_paths(args.ledger_dir)
    result["report"] = render_ledger_report(args.ledger_dir)
    return result


HANDLERS = {
    "doctor": _run_doctor,
    "provenance-ledger": _run_provenance_ledger,
}
