"""CLI commands for the deterministic blast-panel flow (issue #374 / #355).

Three commands wrap the closed, auditable blast layer:

* ``blast-panel-freeze``: the panel TSV is the only input; interval,
  candidate, truth and score inputs are structurally absent.
* ``blast-panel-verify``: re-hashes the frozen artifacts and refuses drift.
* ``blast-panel-search``: BLASTs only verified frozen baits against the
  interval proteins; missing tools or provider failure degrade to
  ``not_evaluated_*`` with exit code 0, never a miss.
"""

from __future__ import annotations

from ..blast_evidence import run_panel_blastp
from ..blast_panel import freeze_blast_panel, verify_blast_panel_freeze
from argparse import _SubParsersAction
from typing import Any
from ._shared import add_json_flag


def register(subparsers: _SubParsersAction) -> None:
    sub = subparsers.add_parser("blast-panel-freeze", help="Freeze a panel TSV into an immutable bait panel plus a manifest. Only the panel is read; interval, candidate, truth and score inputs are structurally absent.")
    sub.add_argument("--panel", required=True, help="Panel TSV in frozen_bait_panel.tsv format.")
    sub.add_argument("--out-dir", required=True, help="Output directory for frozen_bait_panel.tsv and blast_panel_freeze_manifest.json.")
    sub.add_argument("--case-id", default="")
    sub.add_argument("--min-panel-size", type=int, default=2, help="Panels with fewer verified rows than this are review_required_small_panel, not rejected.")
    add_json_flag(sub)

    sub = subparsers.add_parser("blast-panel-verify", help="Re-hash a blast-panel freeze manifest against its frozen artifacts and refuse any drift.")
    sub.add_argument("--freeze", required=True, help="blast_panel_freeze_manifest.json.")
    add_json_flag(sub)

    sub = subparsers.add_parser("blast-panel-search", help="BLAST the verified frozen baits against the interval proteins and emit auditable blast evidence. Missing tools or a provider failure degrade to not_evaluated_*, never a miss.")
    sub.add_argument("--freeze", required=True, help="blast_panel_freeze_manifest.json; verified before any BLAST tool runs.")
    sub.add_argument("--proteins", required=True, help="Interval protein FASTA; only opened after the freeze verifies.")
    sub.add_argument("--out-dir", required=True, help="Output directory for blast_evidence.tsv and blast_panel_search_manifest.json.")
    sub.add_argument("--runner", choices=["auto", "native", "wsl"], default="auto", help="auto: native when BLAST+ is on PATH, otherwise WSL.")
    sub.add_argument("--threads", type=int, default=4, help="blastp -num_threads.")
    add_json_flag(sub)


def _run_blast_panel_freeze(args: Any) -> dict:
    return freeze_blast_panel(args.panel, args.out_dir, args.case_id, min_panel_size=args.min_panel_size)


def _run_blast_panel_verify(args: Any) -> dict:
    return verify_blast_panel_freeze(args.freeze)


def _run_blast_panel_search(args: Any) -> dict:
    return run_panel_blastp(args.freeze, args.proteins, args.out_dir, args.threads, args.runner)


HANDLERS = {
    "blast-panel-freeze": _run_blast_panel_freeze,
    "blast-panel-verify": _run_blast_panel_verify,
    "blast-panel-search": _run_blast_panel_search,
}
