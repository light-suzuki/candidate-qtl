"""Locus causal scenario ledger commands (issue #362).

``scenario-create`` builds the canonical ``locus_causal_scenarios.tsv`` from TSV
or YAML scenario rows; ``scenario-report`` renders the ledger as a review
summary. The ledger is separate from gene ranking and never selects a winner:
every row keeps ``rank_contribution=none`` and evidence-free rows stay
``unassessed``.
"""

from __future__ import annotations

from argparse import _SubParsersAction
from pathlib import Path
from typing import Any

import pandas as pd

from ..causal_scenarios import (
    ASSESSMENT_STATES,
    build_scenario_ledger,
    load_scenario_rows,
    render_scenario_report,
    scenario_summary,
    write_scenario_ledger,
)
from ._shared import add_json_flag


def register(subparsers: _SubParsersAction) -> None:
    sub = subparsers.add_parser("scenario-create", help="Build the canonical locus causal scenario ledger (locus_causal_scenarios.tsv) from scenario rows; evidence-free rows stay unassessed and every row keeps rank_contribution=none.")
    sub.add_argument("--input", required=True, help="Scenario rows as TSV or YAML (list of mappings).")
    sub.add_argument("--case", default="", help="Optional case_id overriding/completing every row's case_id.")
    sub.add_argument("--out-dir", required=True, help="Directory receiving locus_causal_scenarios.tsv.")
    add_json_flag(sub)

    sub = subparsers.add_parser("scenario-report", help="Render the locus causal scenario ledger as a review summary without choosing a winner.")
    sub.add_argument("--ledger", required=True, help="locus_causal_scenarios.tsv from scenario-create.")
    add_json_flag(sub)


def _run_scenario_create(args: Any) -> dict:
    rows = load_scenario_rows(args.input)
    ledger = build_scenario_ledger(rows, case_id=args.case)
    out_path = write_scenario_ledger(ledger, Path(args.out_dir) / "locus_causal_scenarios.tsv", case_id=args.case)
    return {
        "out_path": str(out_path),
        "scenario_count": int(len(ledger)),
        "case_ids": sorted(ledger["case_id"].unique().tolist()),
        "assessment_breakdown": {
            state: int((ledger["assessment_state"] == state).sum())
            for state in sorted(ASSESSMENT_STATES)
        },
    }


def _run_scenario_report(args: Any) -> dict:
    ledger = pd.read_csv(args.ledger, sep="\t", dtype=str, keep_default_na=False)
    summary = scenario_summary(ledger)
    if not args.json:
        print(summary.to_string(index=False))
        print()
        print(render_scenario_report(ledger), end="")
        return {}
    return {"scenario_count": int(len(ledger)), "summary": summary.to_dict(orient="records")}


HANDLERS = {
    "scenario-create": _run_scenario_create,
    "scenario-report": _run_scenario_report,
}
