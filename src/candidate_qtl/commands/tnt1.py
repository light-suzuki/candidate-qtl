"""CLI command for curated Medicago Tnt1 evidence (Issue #524)."""

from __future__ import annotations

from argparse import _SubParsersAction
from typing import Any

from ..tnt1_evidence import normalize_tnt1_evidence
from ._shared import add_json_flag


def register(subparsers: _SubParsersAction) -> None:
    sub = subparsers.add_parser(
        "tnt1-evidence-normalize",
        help=(
            "Normalize curated Medicago Tnt1 evidence while keeping FST/insertion "
            "observations distinct from causal mutant evidence and preserving "
            "cross-accession identity-bridge requirements."
        ),
    )
    sub.add_argument("--evidence", required=True, help="Curated Tnt1 TSV following templates/tnt1_evidence.tsv.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument(
        "--registry",
        help=(
            "Optional reference-resource registry used only to label exact native "
            "namespaces as registered/unregistered; it does not infer identity."
        ),
    )
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)


def _run_tnt1_evidence_normalize(args: Any) -> dict:
    return normalize_tnt1_evidence(
        args.evidence,
        args.out_dir,
        registry_arg=args.registry,
        case_id=args.case_id,
    )


HANDLERS = {"tnt1-evidence-normalize": _run_tnt1_evidence_normalize}
