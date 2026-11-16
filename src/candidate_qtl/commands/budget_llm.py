"""CLI command for the batched LLM cluster review + critic + portfolio reducer (issue #375 T5).

Kept as a separate ``budget-`` prefixed subcommand rather than folding into the
older Stage-B ``llm-*`` commands in ``commands/llm.py``: this command runs
against an EXISTING ``budget-run`` output directory (tier0/tier1/tier2/handoff
already on disk) and rewrites only the LLM-review and Top-20 artifacts, never
touching rank/weights/shortlist state.
"""

from __future__ import annotations

from argparse import _SubParsersAction
from typing import Any

import pandas as pd

from ..budget_account import BUDGET_PROFILES, BudgetAccount
from ..budget_cache import ContentCache
from ..core import runtime_path
from ..llm_portfolio import FileReviewBackend, audit_proposer_diversity, build_critic, build_review_batches, reduce_portfolio, run_cluster_review
from ._shared import add_json_flag


def register(subparsers: _SubParsersAction) -> None:
    sub = subparsers.add_parser(
        "budget-llm-review",
        help="Batched, budget-gated LLM cluster review + critic + portfolio reducer over an existing budget-run output directory (issue #375 T5).",
    )
    sub.add_argument("--out-dir", required=True, help="Existing budget-run output directory (tier0_inventory.tsv/tier1_clusters.tsv/tier2_plan.tsv/budget_handoff.tsv already present).")
    sub.add_argument("--case-id", default="case_001", help="Case identifier recorded in every manifest.")
    sub.add_argument("--responses", required=True, help="Canned review responses JSON for the offline FileReviewBackend; never a network call.")
    sub.add_argument("--mechanisms", default=None, help="mechanism_hypotheses.tsv from `mechanism-import`; enables the proposer diversity audit and the mechanism-id allowlist.")
    sub.add_argument("--allowed-sources", default=None, help="TSV with a source_id column of verified source ids; enables the rationale source-id allowlist.")
    sub.add_argument("--bridge", default=None, help="family_candidate_bridge.tsv (from `family-candidate-bridge`); supplies the route each candidate rode in on so the portfolio reducer can allocate its route-unique rescue slots. Without it the reducer records no_route_information_available instead of assuming a route.")
    sub.add_argument("--budget-profile", choices=sorted(BUDGET_PROFILES), default="student-budget", help="Budget profile preset (default: student-budget).")
    sub.add_argument("--batch-size", type=int, default=5, help="Maximum Tier1 clusters packed per review batch.")
    sub.add_argument("--cache-dir", default=None, help="Content-addressed cache directory; default <out-dir>/budget_cache.")
    sub.add_argument("--no-cache", action="store_true", help="Disable the content cache (recompute every batch).")
    sub.add_argument("--portfolio-size", type=int, default=20, help="Top-N portfolio size for the reducer.")
    add_json_flag(sub)


def _read_tsv(path) -> list[dict[str, object]]:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("").to_dict("records")


def _run_budget_llm_review(args: Any) -> dict:
    out = runtime_path(args.out_dir)
    handoff_rows = _read_tsv(out / "budget_handoff.tsv")
    plan_rows = _read_tsv(out / "tier2_plan.tsv")
    cluster_rows = _read_tsv(out / "tier1_clusters.tsv")

    account = BudgetAccount()
    account.declare(BUDGET_PROFILES[args.budget_profile])
    cache = ContentCache(runtime_path(args.cache_dir) if args.cache_dir else out / "budget_cache", enabled=not args.no_cache)

    mechanism_ids: set[str] | None = None
    diversity_result: dict[str, object] | None = None
    if args.mechanisms:
        mechanism_rows = _read_tsv(runtime_path(args.mechanisms))
        mechanism_ids = {str(row.get("hypothesis_id", "")).strip() for row in mechanism_rows if str(row.get("hypothesis_id", "")).strip()}
        diversity_result = audit_proposer_diversity(mechanism_rows, out, args.case_id)

    allowed_source_ids: set[str] | None = None
    if args.allowed_sources:
        source_rows = _read_tsv(runtime_path(args.allowed_sources))
        allowed_source_ids = {str(row.get("source_id", "")).strip() for row in source_rows if str(row.get("source_id", "")).strip()}

    backend = FileReviewBackend(args.responses)
    batches = build_review_batches(handoff_rows, plan_rows, cluster_rows, batch_size=args.batch_size)
    review_manifest = run_cluster_review(
        batches, account, out, args.case_id, backend,
        mechanisms=mechanism_ids, allowed_source_ids=allowed_source_ids, cache=cache,
    )
    review_rows = _read_tsv(out / "llm_review.tsv")
    critic_result = build_critic(review_rows, batches, handoff_rows, cluster_rows, account, out, args.case_id)
    bridge_rows = _read_tsv(runtime_path(args.bridge)) if args.bridge else None
    portfolio_result = reduce_portfolio(
        handoff_rows, plan_rows, cluster_rows, out, args.case_id,
        review_rows=review_rows, bridge_rows=bridge_rows, portfolio_size=args.portfolio_size,
    )

    outputs = {
        "llm_review_tsv": "llm_review.tsv",
        "llm_review_manifest_json": "llm_review_manifest.json",
        "llm_review_critic_tsv": "llm_review_critic.tsv",
        "deterministic_top20_tsv": "deterministic_top20.tsv",
        "llm_assisted_top20_tsv": "llm_assisted_top20.tsv",
        "portfolio_manifest_json": "portfolio_manifest.json",
    }
    if diversity_result is not None:
        outputs["proposer_diversity_audit_tsv"] = "proposer_diversity_audit.tsv"

    return {
        "schema_version": "1.0",
        "case_id": args.case_id,
        "review": review_manifest,
        "critic_finding_count": len(critic_result["findings"]),
        "portfolio": portfolio_result["manifest"],
        "outputs": outputs,
        "exit_code": 0,
    }


HANDLERS = {
    "budget-llm-review": _run_budget_llm_review,
}
