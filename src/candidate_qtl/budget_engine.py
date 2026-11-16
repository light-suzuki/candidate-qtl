"""Deterministic budget engine (issue #375): tier0 -> tier1 -> tier2 -> handoff.

:func:`run_budget_engine` runs one fully deterministic, budget-gated pass over
an extract directory: Tier 0 inventory, Tier 1 clustering, the Tier 2
sequence-search plan (budgeted execution through the content-addressed cache),
a per-gene review handoff, the resource-usage ledger, and optionally a second
exhaustive (unlimited) arm plus the gating comparison. Every artifact is a
pure function of the inputs plus the ledger: no wall-clock values inside
manifests except the documented runtime carve-out (see below).

Status rule: ``budget_exhausted`` iff ``queries_deferred_budget > 0`` AND
(``queries_executed > 0`` OR a *positive* declared budget was found exhausted
in the frozen ledger -- see :func:`_status_for`). A run whose declared
limits are all zero or all unlimited never counts as exhausted on its own
(a zero/unset budget is by-design work-free), so it completes with a full
inventory, plan, and handoff even though every planned query is deferred.
But a run that starts over a positive declared limit (e.g. peak memory
already over budget before the first query) and therefore defers every
query without executing any is now correctly ``budget_exhausted`` too --
previously it was mislabelled ``completed`` because the rule only looked at
``queries_executed > 0``.

Runtime carve-out (kept out of ``tier0/tier1/tier2_plan/budget_handoff`` so
those artifacts stay byte-identical across reruns): measured wall time, peak
RSS, and sequence-search CPU seconds are recorded in
``benchmark_resource_usage.json`` ``usage.wall_time_seconds`` /
``usage.peak_memory_mb`` / ``usage.sequence_search_cpu_seconds``. The Tier 2
execution loop charges the declared wall-time budget per query (the measured
elapsed time via :meth:`BudgetAccount.consume`) and records the declared
peak-memory limit as an *observed maximum* (via
:meth:`BudgetAccount.observe_peak` -- never summed, never consumable, so a
truncation-causing overrun and the ledger always agree). The process
high-water mark measured at the end of the run (``resource.getrusage``) is
folded into that same observed maximum instead of being written over the
report, so ``budgets.peak_memory_limit_mb`` in ``budget_manifest.json`` and
``usage.peak_memory_mb`` in ``benchmark_resource_usage.json`` are always the
same quantity -- a run stopped for memory can never report a peak under its
own limit. Consequently ``budget_manifest.json`` is no longer byte-identical
across reruns either:
its ``budgets.wall_time_limit_seconds`` and ``budgets.peak_memory_limit_mb``
consumed/remaining/exhausted values are runtime measurements. Byte-identity
therefore holds for tier0/tier1/tier2_plan/handoff artifacts only;
``budget_manifest.json`` and ``benchmark_resource_usage.json`` are compared
structurally with the wall/memory consumption zeroed. ``sequence_search_cpu_seconds``
is charged but deterministic: a cache hit re-adds the originally measured
value stored in the payload, so warm reruns record the same cpu cost as the
cold run. Warm reruns intentionally change ``cache`` stats and
``sequence_query_count`` in ``budget_manifest.json``.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import time
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from .budget_account import BUDGET_PROFILES, LEDGER_LIMIT_FIELDS, LEDGER_USAGE_FIELDS, BudgetAccount
from .budget_cache import ContentCache
from .core import runtime_path, run_tool
from .llm_portfolio import audit_proposer_diversity, build_critic, build_review_batches, reduce_portfolio, run_cluster_review
from .tier0_inventory import build_tier0_inventory
from .tier1_cluster import build_tier1_clusters
from .tier2_planner import plan_sequence_search
# Imported as a module (not a bound name) so the RSS instrument is the one
# tests patch on tier2_planner, keeping one definition of "current RSS".
from . import tier2_planner as _tier2

HANDOFF_COLUMNS = (
    "gene_id",
    "cluster_id",
    "allocation",
    "rescue_bucket",
    "plan_step",
    "query_role",
    "query_status",
    "expansion_condition",
    "search_result_status",
    "review_state",
    "handoff_note",
)

GATING_COLUMNS = (
    "gene_id",
    "cluster_id",
    "hierarchical_query_status",
    "exhaustive_query_status",
    "hierarchical_hit_presence",
    "exhaustive_hit_presence",
    "gating_gap",
)

GATING_HEADER_COMMENT = "# measurement hook: reports gating differences; not a claim of causality"


def _resolve_tool(bin: str) -> tuple[str | None, str]:
    """Resolve an executable and capture its version.

    ``bin`` may be a bare command resolved through PATH or a direct file path;
    the resolved real path is returned. The version is the first non-empty
    line of ``[resolved, "-version"]`` (stdout+stderr, capped at 80 chars),
    run under a 10 s timeout; any failure or non-zero exit yields
    ``"unknown"``. Unresolvable input returns ``(None, "unknown")``.
    """
    if not bin:
        return None, "unknown"
    resolved = shutil.which(bin)
    if resolved is None and Path(bin).is_file():
        resolved = str(Path(bin).resolve())
    if resolved is None:
        return None, "unknown"
    resolved = os.path.realpath(resolved)
    version = "unknown"
    try:
        completed = subprocess.run([resolved, "-version"], capture_output=True, text=True, timeout=10, check=False)
    except (OSError, subprocess.TimeoutExpired):
        return resolved, version
    if completed.returncode != 0:
        return resolved, version
    for line in (completed.stdout + completed.stderr).splitlines():
        if line.strip():
            version = line.strip()[:80]
            break
    return resolved, version


class BlastpBackend:
    """Local BLASTP runner with a deterministic outfmt-6 parser.

    ``available`` is probed once at construction by resolving the actual
    binary (``shutil.which`` or a direct path) and running
    ``<bin> -version``; a missing binary is a provider gap (the engine plans
    ``not_evaluated`` and consumes nothing), so ``search`` raises
    ``RuntimeError("blastp not found")`` only as a defensive guard. ``search``
    runs under ``timeout`` (defaults to the engine's external-tool timeout);
    a stall raises :class:`candidate_qtl.core.ExternalDependencyError` with
    ``kind="timeout"``.
    """

    def __init__(self, cache: ContentCache | None = None, bin: str = "blastp") -> None:
        self.cache = cache
        self.resolved_bin, self.version = _resolve_tool(bin)
        self.available = self.resolved_bin is not None

    def search(self, query_fasta: Path, db: str, bin: str, timeout: float | None = None) -> tuple[list[dict[str, object]], float]:
        """Run ``blastp`` against ``db``; returns (hits sorted by (sseqid, qseqid), cpu_seconds)."""
        if not self.available:
            raise RuntimeError("blastp not found")
        command = [
            self.resolved_bin,
            "-query", str(query_fasta),
            "-db", str(db),
            "-outfmt", "6 qseqid sseqid pident evalue bitscore",
            "-max_target_seqs", "5",
            "-evalue", "10",
            "-num_threads", "1",
        ]
        start = time.perf_counter()
        result = run_tool(command, timeout=timeout, capture_output=True, text=True, check=False)
        elapsed = time.perf_counter() - start
        if result.returncode != 0:
            raise RuntimeError(f"blastp failed (exit {result.returncode}): {result.stderr.strip()[:200]}")
        hits: list[dict[str, object]] = []
        for line in result.stdout.splitlines():
            fields = line.split("\t")
            if len(fields) < 5:
                continue
            hits.append({"qseqid": fields[0], "sseqid": fields[1], "pident": fields[2], "evalue": fields[3], "bitscore": fields[4]})
        hits.sort(key=lambda hit: (str(hit["sseqid"]), str(hit["qseqid"])))
        return hits, elapsed


def _status_for(queries_deferred_budget: int, queries_executed: int, positive_budget_exhausted: bool = False) -> str:
    """Run status from the plan counters plus whether a positive budget is exhausted.

    ``positive_budget_exhausted`` is computed by the caller from the frozen
    ledger: any :data:`~candidate_qtl.budget_account.LEDGER_LIMIT_FIELDS` name
    whose declared limit is not ``None`` and ``> 0`` and whose ``exhausted``
    flag is True. It exists so a run that starts over a positive declared
    limit (e.g. peak memory already over budget) and defers every query
    before executing any is still reported ``budget_exhausted`` instead of
    ``completed``.
    """
    if queries_deferred_budget <= 0:
        return "completed"
    if queries_executed > 0:
        return "budget_exhausted"
    # executed == 0: an explicitly all-zero / zero-limit deterministic profile is
    # by-design work-free and stays "completed"; only the real overdraw or refusal
    # of a POSITIVE declared budget is an exhaustion.
    return "budget_exhausted" if positive_budget_exhausted else "completed"


def _hit_presence(search_result_status: str) -> str:
    return search_result_status if search_result_status in {"hit", "no_hit"} else "not_run"


def _peak_memory_mb() -> tuple[float, str | None]:
    try:
        import resource

        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0, None
    except (ImportError, AttributeError):
        return 0.0, "not_obtainable"


def _observe_run_peak(
    account: BudgetAccount,
    start_peak_mb: float | None,
    end_peak_mb: float,
    end_peak_note: str | None,
) -> None:
    """Fold only run-scoped memory measurements into the ledger's observed peak.

    Two sources are attributable to the run that just executed:

    * the end-of-run RSS from ``/proc/self/statm`` -- the same instrument the
      Tier 2 and review gates sample, so the ledger's peak keeps one meaning
      across the whole engine and is populated even when no gate ran (an
      undeclared peak limit means nothing samples RSS during the run);
    * ``ru_maxrss``, but ONLY when it grew past the baseline taken at entry.
      It is a whole-process lifetime high-water mark, so a value set by an
      earlier operation in a long-lived process says nothing about this run
      and must never mark this run's peak budget exhausted.

    ``observe_peak`` keeps a maximum, so nothing is summed and nothing already
    observed by the gates is lost.
    """
    current_rss, _note = _tier2._current_rss_mb()
    if current_rss is not None:
        account.observe_peak("peak_memory_budget", current_rss)
    if end_peak_note is None and start_peak_mb is not None and end_peak_mb > start_peak_mb:
        account.observe_peak("peak_memory_budget", end_peak_mb)


def _benchmark_resource_usage(
    case_id: str,
    arm: str,
    profile: str,
    status: str,
    limits: dict[str, object],
    usage: dict[str, object],
    heavy_stages: list[dict[str, object]] | None = None,
) -> dict[str, object]:
    result: dict[str, object] = {
        "schema_version": "1.0",
        "case_id": case_id,
        "benchmark_arm": arm,
        "execution_profile": profile,
        "status": status,
        "limits": {name: limits[name] for name in LEDGER_LIMIT_FIELDS},
        "usage": {name: usage[name] for name in LEDGER_USAGE_FIELDS},
    }
    if heavy_stages:
        result["heavy_stages"] = [dict(entry) for entry in heavy_stages]
    return result


def summarize_case_cost_telemetry(
    account: BudgetAccount,
    *,
    initial_candidate_count: int | None = None,
    final_candidate_count: int | None = None,
    wet_budget: int = 20,
    case_id: str = "",
    benchmark_arm: str = "",
    execution_profile: str = "",
    status: str = "completed",
) -> dict[str, object]:
    """Aggregate per-case compute cost and candidate reduction telemetry for Issue #434.

    Gathers:
    * wall time (seconds)
    * process CPU time (seconds)
    * peak memory (MB)
    * provider-reported LLM calls, input tokens, output tokens
    * cache hits and misses
    * BLAST / sequence queries and searched proteins
    * external retrieval queries
    * heavy-stage invocation count and list of invocation records with invocation_reason
    * candidate reduction metrics (initial, final, reduction ratio, wet_budget compliance)
    """
    frozen = account.freeze()
    usage = frozen["usage"]
    limits = frozen["limits"]

    initial = int(initial_candidate_count) if initial_candidate_count is not None else None
    final = int(final_candidate_count) if final_candidate_count is not None else None
    reduction_ratio: float | None = None
    if initial is not None and initial > 0 and final is not None:
        reduction_ratio = round(max(0.0, 1.0 - (final / float(initial))), 4)

    within_wet_budget: bool | None = None
    if final is not None:
        within_wet_budget = final <= wet_budget

    heavy_stage_list = account.heavy_stages()
    reasons = [str(entry.get("invocation_reason", "")) for entry in heavy_stage_list if entry.get("invocation_reason")]

    return {
        "schema_version": "1.0",
        "case_id": case_id,
        "benchmark_arm": benchmark_arm,
        "execution_profile": execution_profile,
        "status": status,
        "candidate_counts": {
            "initial": initial,
            "final": final,
            "wet_budget": int(wet_budget),
            "within_wet_budget": within_wet_budget,
            "reduction_ratio": reduction_ratio,
        },
        "compute": {
            "wall_time_seconds": usage["wall_time_seconds"],
            "cpu_time_seconds": usage["cpu_time_seconds"] or usage["sequence_search_cpu_seconds"],
            "peak_memory_mb": usage["peak_memory_mb"],
            "llm_call_count": usage["llm_call_count"],
            "llm_input_tokens": usage["llm_input_tokens"],
            "llm_output_tokens": usage["llm_output_tokens"],
            "llm_total_tokens": usage["llm_input_tokens"] + usage["llm_output_tokens"],
            "cache_hit_count": usage["cache_hit_count"],
            "cache_miss_count": usage["cache_miss_count"],
            "sequence_query_count": usage["sequence_query_count"],
            "searched_protein_count": usage["searched_protein_count"],
            "external_retrieval_count": usage["external_retrieval_count"],
            "heavy_stage_invocation_count": usage["heavy_stage_invocation_count"],
        },
        "heavy_stages": heavy_stage_list,
        "invocation_reasons": sorted(set(reasons)),
        "limits": {name: limits[name] for name in LEDGER_LIMIT_FIELDS},
    }


def build_cost_profile_comparison_report(
    profile_summaries: list[dict[str, object]],
    *,
    case_roster: list[str] | None = None,
) -> dict[str, object]:
    """Compare multiple profile runs (e.g. light, standard, heavy) on the same roster.

    Computes:
    * median candidate reduction
    * within-budget completion rate
    * heavy-stage invocation rate
    * total & max resource measurements across runs
    """
    profiles: dict[str, dict[str, object]] = {}
    for summary in profile_summaries:
        pname = str(summary.get("execution_profile", "unknown"))
        if pname not in profiles:
            profiles[pname] = {
                "case_count": 0,
                "total_wall_time_seconds": 0.0,
                "max_peak_memory_mb": 0.0,
                "total_llm_calls": 0,
                "total_llm_tokens": 0,
                "total_heavy_stage_invocations": 0,
                "within_wet_budget_count": 0,
                "reductions": [],
                "cases": [],
            }
        pdata = profiles[pname]
        pdata["case_count"] += 1
        comp = summary.get("compute", {})
        cand = summary.get("candidate_counts", {})
        pdata["total_wall_time_seconds"] = round(pdata["total_wall_time_seconds"] + float(comp.get("wall_time_seconds", 0.0)), 4)
        pdata["max_peak_memory_mb"] = max(pdata["max_peak_memory_mb"], float(comp.get("peak_memory_mb", 0.0)))
        pdata["total_llm_calls"] += int(comp.get("llm_call_count", 0))
        pdata["total_llm_tokens"] += int(comp.get("llm_total_tokens", 0))
        pdata["total_heavy_stage_invocations"] += int(comp.get("heavy_stage_invocation_count", 0))
        if cand.get("within_wet_budget"):
            pdata["within_wet_budget_count"] += 1
        if cand.get("reduction_ratio") is not None:
            pdata["reductions"].append(float(cand["reduction_ratio"]))
        pdata["cases"].append(str(summary.get("case_id", "")))

    aggregated: dict[str, object] = {}
    for pname, pdata in profiles.items():
        n = pdata["case_count"] or 1
        reds = sorted(pdata["reductions"])
        median_red = reds[len(reds) // 2] if reds else 0.0
        aggregated[pname] = {
            "case_count": pdata["case_count"],
            "median_reduction_ratio": median_red,
            "within_wet_budget_rate": round(pdata["within_wet_budget_count"] / n, 4),
            "heavy_stage_invocation_rate": round(pdata["total_heavy_stage_invocations"] / n, 4),
            "total_wall_time_seconds": pdata["total_wall_time_seconds"],
            "max_peak_memory_mb": pdata["max_peak_memory_mb"],
            "total_llm_tokens": pdata["total_llm_tokens"],
            "total_llm_calls": pdata["total_llm_calls"],
        }

    return {
        "schema_version": "1.0",
        "roster": sorted(set(case_roster or [])) if case_roster else None,
        "profiles": aggregated,
    }


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")


def run_budget_engine(
    extract_dir: str | Path,
    out_dir: str | Path,
    case_id: str,
    budget_limits: dict[str, object],
    providers: dict[str, dict[str, object]],
    cache_dir: str | Path | None = None,
    no_cache: bool = False,
    protein_db: str | Path | None = None,
    blastp_bin: str = "blastp",
    flag_genes: list[str] | None = None,
    exhaustive: bool = False,
    search_backend: object | None = None,
    review_backend: object | None = None,
    llm_responses: str | Path | None = None,
    mechanisms: str | Path | None = None,
    bridge: str | Path | None = None,
    portfolio_size: int = 20,
) -> dict[str, object]:
    """Run the full budgeted engine; returns the budget manifest (see step 8).

    ``providers`` maps provider kind to the dict returned by
    :func:`candidate_qtl.tier0_inventory.load_provider_table` (absent kinds are
    simply not provided). ``budget_limits`` are ledger or conceptual limit
    names (see :class:`candidate_qtl.budget_account.BudgetAccount`); when the
    dict equals one of the :data:`BUDGET_PROFILES` presets the profile name is
    recorded, otherwise ``"custom"``. ``flag_genes`` forces the Tier 2
    expansion trigger for the listed genes. ``exhaustive`` runs a second
    unlimited arm under ``<out-dir>/exhaustive`` (reusing the cache) and the
    gating comparison at ``<out-dir>/budget_gating_comparison.tsv``.

    ``outputs`` maps artifact names to file names relative to ``out_dir``
    (never absolute paths), so the deterministic artifacts
    (``tier0/tier1/tier2_plan/budget_handoff``) are byte-identical across
    reruns into different output directories. ``budget_manifest.json`` and
    ``benchmark_resource_usage.json`` embed runtime wall/memory measurements
    and are byte-identical only when compared with those values zeroed.
    """
    extract = runtime_path(extract_dir)
    out = runtime_path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    account = BudgetAccount()
    account.declare(budget_limits)
    start_time = time.perf_counter()
    # ``ru_maxrss`` is the high-water mark of the whole PROCESS, not of this
    # run, so it is only attributable to this run when the run itself raises
    # it. The baseline taken here is what makes that comparison possible; a
    # mark set by an earlier operation in a long-lived process must never
    # charge this run's peak budget.
    start_peak_mb, start_peak_note = _peak_memory_mb()

    # Tier 0 inventory, then Tier 1 clustering over the written rows.
    build_tier0_inventory(extract, out, providers, case_id)
    inventory_frame = pd.read_csv(out / "tier0_inventory.tsv", sep="\t", dtype=str).fillna("")
    if "cds_length" not in inventory_frame.columns:
        inventory_frame["cds_length"] = "0"
    inventory_rows = inventory_frame.to_dict("records")
    build_tier1_clusters(inventory_rows, out, case_id, providers)
    cluster_frame = pd.read_csv(out / "tier1_clusters.tsv", sep="\t", dtype=str).fillna("")
    cluster_rows = cluster_frame.to_dict("records")

    # Cache + backend; an injected backend always gets the engine cache attached
    # when it has none, so fake backends exercise the cache path identically.
    cache = ContentCache(runtime_path(cache_dir) if cache_dir else out / "budget_cache", enabled=not no_cache)
    if search_backend is None:
        search_backend = BlastpBackend(cache=cache, bin=blastp_bin)
    elif getattr(search_backend, "cache", None) is None:
        search_backend.cache = cache
    backend = search_backend

    proteins: dict[str, str] = {}
    proteins_path = extract / "candidate_proteins.fa"
    if proteins_path.exists():
        with proteins_path.open(encoding="utf-8") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                proteins[record.id] = str(record.seq)

    profile_name = next((name for name, preset in BUDGET_PROFILES.items() if preset == budget_limits), "custom")
    plan = plan_sequence_search(
        inventory_rows,
        cluster_rows,
        account,
        out,
        case_id,
        user_flagged=set(flag_genes or []),
        protein_db=protein_db,
        blastp_bin=blastp_bin,
        search_backend=backend,
        proteins=proteins,
    )
    plan_rows = pd.read_csv(out / "tier2_plan.tsv", sep="\t", dtype=str).fillna("").to_dict("records")
    inventory_by_gene = {str(row.get("gene_id", "")): row for row in inventory_rows}
    cluster_by_gene: dict[str, dict[str, object]] = {}
    for cluster in cluster_rows:
        for member in str(cluster.get("member_gene_ids", "")).split(";"):
            member = member.strip()
            if member:
                cluster_by_gene[member] = cluster

    handoff_rows = []
    for row in plan_rows:
        gene_id = str(row.get("gene_id", ""))
        cluster = cluster_by_gene.get(gene_id, {})
        query_status = str(row.get("query_status", ""))
        allocation = str(cluster.get("allocation", ""))
        search_result_status = str(row.get("search_result_status", ""))
        # A provider failure is checked before the allocation branches: a
        # rescue-allocation row that failed against the provider must never
        # be labelled rescue_pending_review (that would hide the failure
        # behind ordinary rescue triage) and is never a reviewed candidate.
        if query_status == "provider_failed":
            review_state = "not_evaluated"
        elif query_status == "executed":
            review_state = "priority_review"
        elif query_status in {"deferred_budget", "deferred_plan", "not_planned", "not_queried"} and allocation != "rescue":
            review_state = "deferred"
        elif allocation == "rescue":
            review_state = "rescue_pending_review"
        else:
            review_state = "not_evaluated"
        if query_status == "provider_failed":
            handoff_note = "provider_timeout"
        elif search_result_status == "not_evaluated":
            handoff_note = "not_evaluated_provider_missing"
        elif search_result_status == "provider_timeout":
            handoff_note = "provider_timeout"
        elif query_status == "executed":
            handoff_note = "searched"
        elif query_status == "deferred_budget":
            handoff_note = "budget_exhausted"
        elif query_status == "not_planned":
            handoff_note = str(row.get("deferred_reason", ""))
        else:
            handoff_note = "not_searched"
        handoff_rows.append(
            {
                "gene_id": gene_id,
                "cluster_id": str(row.get("cluster_id", "")),
                "allocation": allocation,
                "rescue_bucket": str(cluster.get("rescue_bucket", "")),
                "plan_step": str(row.get("plan_step", "")),
                "query_role": str(row.get("query_role", "")),
                "query_status": query_status,
                "expansion_condition": str(row.get("expansion_condition", "")),
                "search_result_status": search_result_status,
                "review_state": review_state,
                "handoff_note": handoff_note,
            }
        )
    handoff_rows.sort(
        key=lambda entry: (str(inventory_by_gene.get(entry["gene_id"], {}).get("chrom", "")), int(str(inventory_by_gene.get(entry["gene_id"], {}).get("start", "")) or 0), entry["gene_id"])
    )
    pd.DataFrame(handoff_rows, columns=HANDOFF_COLUMNS).to_csv(out / "budget_handoff.tsv", sep="\t", index=False, lineterminator="\n")
    reviewed_candidate_count = sum(1 for entry in handoff_rows if entry["review_state"] == "priority_review")
    unreviewed_candidate_count = len(handoff_rows) - reviewed_candidate_count
    account.add_usage("reviewed_candidate_count", reviewed_candidate_count)
    account.add_usage("unreviewed_candidate_count", unreviewed_candidate_count)

    # T5 LLM review + portfolio reduction. A positive declared
    # llm_call_limit plus a supplied backend (or --llm-responses file) runs
    # the batched, budget-gated cluster review and critic; otherwise it is
    # skipped and every candidate stays unreviewed. The portfolio reducer
    # always runs so a canonical run always carries both Top-20 artifacts.
    llm_review_manifest: dict[str, object] | None = None
    review_rows: list[dict[str, object]] = []
    diversity_result: dict[str, object] | None = None
    declared_llm_call_limit = account.freeze()["limits"]["llm_call_limit"]
    has_positive_llm_budget = declared_llm_call_limit is None or declared_llm_call_limit > 0

    resolved_review_backend = review_backend
    if resolved_review_backend is None and llm_responses is not None:
        from .llm_portfolio import FileReviewBackend

        resolved_review_backend = FileReviewBackend(llm_responses)

    # The canonical family_candidate_bridge supplies the route (a/b/c) each
    # candidate rode in on; without it the reducer's route-unique rescue slot
    # has nothing to allocate and says so instead of pretending (T5 audit, P2).
    bridge_rows_loaded: list[dict[str, object]] | None = None
    if bridge is not None:
        bridge_rows_loaded = pd.read_csv(runtime_path(bridge), sep="\t", dtype=str).fillna("").to_dict("records")

    mechanism_ids: set[str] | None = None
    if mechanisms is not None:
        mechanism_rows = pd.read_csv(runtime_path(mechanisms), sep="\t", dtype=str).fillna("").to_dict("records")
        mechanism_ids = {str(row.get("hypothesis_id", "")).strip() for row in mechanism_rows if str(row.get("hypothesis_id", "")).strip()}
        diversity_result = audit_proposer_diversity(mechanism_rows, out, case_id)

    if not has_positive_llm_budget:
        llm_review_status = "not_evaluated_no_budget"
    elif resolved_review_backend is None:
        llm_review_status = "not_run_no_backend"
    else:
        review_batches = build_review_batches(handoff_rows, plan_rows, cluster_rows)
        llm_review_manifest = run_cluster_review(
            review_batches, account, out, case_id, resolved_review_backend,
            mechanisms=mechanism_ids, cache=cache,
        )
        review_rows = pd.read_csv(out / "llm_review.tsv", sep="\t", dtype=str).fillna("").to_dict("records")
        build_critic(review_rows, review_batches, handoff_rows, cluster_rows, account, out, case_id)
        llm_review_status = "budget_exhausted" if llm_review_manifest["counts"]["deferred"] > 0 else "completed"

    reduce_portfolio(
        handoff_rows, plan_rows, cluster_rows, out, case_id,
        review_rows=review_rows or None, bridge_rows=bridge_rows_loaded, portfolio_size=portfolio_size,
    )

    # Wall time and peak RSS are sampled only AFTER the LLM review stage, so a
    # slow or memory-heavy review is inside the measured run rather than
    # invisible to the resource report (T5 audit, P1). run_cluster_review
    # additionally charges the declared wall budget per batch and observes the
    # peak, so the hard limits stop the review as well, not just Tier 2.
    elapsed = time.perf_counter() - start_time
    peak_memory_mb, peak_memory_note = _peak_memory_mb()
    # The peak that reaches the ledger is folded in, not written past it:
    # ``budgets.peak_memory_limit_mb`` and the resource report's
    # ``usage.peak_memory_mb`` have to be the same quantity, or a run that was
    # stopped for memory reports a peak comfortably under its own limit
    # (issue #375 T4 finding P1 / invariant 2 -- peaks mean the same thing in
    # every manifest). ``observe_peak`` keeps a maximum, so folding never
    # discards the per-query RSS the gates observed and never sums the two.
    #
    # Only run-scoped measurements are folded: the end-of-run RSS (the same
    # /proc instrument the Tier 2 and review gates sample) always, and the
    # process high-water mark ONLY when this run raised it above the baseline
    # taken at entry. A lifetime mark left by an earlier operation in a
    # long-lived process is not this run's peak and must not exhaust its
    # budget.
    _observe_run_peak(account, start_peak_mb if start_peak_note is None else None, peak_memory_mb, peak_memory_note)

    queries_deferred_budget = int(plan["queries_deferred_budget"])
    queries_executed = int(plan["queries_executed"])

    frozen = account.freeze()
    limits = frozen["limits"]
    # A positive declared budget that the frozen ledger marks exhausted makes
    # even a zero-executed run budget_exhausted (e.g. peak memory already
    # over limit before the first query); a genuinely zero/unset budget
    # still completes -- see _status_for.
    positive_budget_exhausted = any(
        limits[name] is not None and limits[name] > 0 and frozen["exhausted"][name]
        for name in LEDGER_LIMIT_FIELDS
    )
    status = _status_for(queries_deferred_budget, queries_executed, positive_budget_exhausted)
    zero_budget = bool(limits) and all(value is not None and value == 0 for value in limits.values())
    plan_mode = str(plan["plan_mode"])
    benchmark_arm = "deterministic" if (plan_mode == "bulk" or zero_budget) else "budgeted_hierarchical"

    usage = dict(frozen["usage"])
    usage["wall_time_seconds"] = round(elapsed, 6)
    # peak_memory_mb comes straight from the ledger (which already absorbed
    # ru_maxrss above), so the report and budget_manifest.json never disagree.
    benchmark = _benchmark_resource_usage(
        case_id, benchmark_arm, profile_name, status, limits, usage,
        heavy_stages=account.heavy_stages() if account.heavy_stages() else None,
    )
    _write_json(out / "benchmark_resource_usage.json", benchmark)

    gating_added_false_negative_count: int | None = None
    if exhaustive:
        exhaustive_dir = out / "exhaustive"
        exhaustive_account = BudgetAccount()
        exhaustive_account.declare(BUDGET_PROFILES["exhaustive"])
        exhaustive_start = time.perf_counter()
        exhaustive_start_peak_mb, exhaustive_start_peak_note = _peak_memory_mb()
        exhaustive_plan = plan_sequence_search(
            inventory_rows,
            cluster_rows,
            exhaustive_account,
            exhaustive_dir,
            case_id,
            user_flagged=set(flag_genes or []),
            protein_db=protein_db,
            blastp_bin=blastp_bin,
            search_backend=backend,
            proteins=proteins,
            exhaustive_override=True,
        )
        exhaustive_elapsed = time.perf_counter() - exhaustive_start
        exhaustive_peak_mb, exhaustive_peak_note = _peak_memory_mb()
        _observe_run_peak(
            exhaustive_account,
            exhaustive_start_peak_mb if exhaustive_start_peak_note is None else None,
            exhaustive_peak_mb,
            exhaustive_peak_note,
        )
        exhaustive_frozen = exhaustive_account.freeze()
        exhaustive_usage = dict(exhaustive_frozen["usage"])
        exhaustive_usage["wall_time_seconds"] = round(exhaustive_elapsed, 6)
        exhaustive_positive_budget_exhausted = any(
            exhaustive_frozen["limits"][name] is not None and exhaustive_frozen["limits"][name] > 0 and exhaustive_frozen["exhausted"][name]
            for name in LEDGER_LIMIT_FIELDS
        )
        exhaustive_benchmark = _benchmark_resource_usage(
            case_id,
            "exhaustive",
            "exhaustive",
            _status_for(
                int(exhaustive_plan["queries_deferred_budget"]),
                int(exhaustive_plan["queries_executed"]),
                exhaustive_positive_budget_exhausted,
            ),
            exhaustive_frozen["limits"],
            exhaustive_usage,
            heavy_stages=exhaustive_account.heavy_stages() if exhaustive_account.heavy_stages() else None,
        )
        _write_json(exhaustive_dir / "benchmark_resource_usage.json", exhaustive_benchmark)

        exhaustive_rows = pd.read_csv(exhaustive_dir / "tier2_plan.tsv", sep="\t", dtype=str).fillna("").to_dict("records")
        hierarchical_by_gene = {str(row.get("gene_id", "")): row for row in plan_rows}
        exhaustive_by_gene = {str(row.get("gene_id", "")): row for row in exhaustive_rows}
        gating_rows = []
        for gene_id in sorted(set(hierarchical_by_gene) | set(exhaustive_by_gene)):
            hrow = hierarchical_by_gene.get(gene_id, {})
            erow = exhaustive_by_gene.get(gene_id, {})
            hstatus = str(hrow.get("query_status", "not_planned"))
            estatus = str(erow.get("query_status", "not_planned"))
            hhit = _hit_presence(str(hrow.get("search_result_status", "not_run")))
            ehit = _hit_presence(str(erow.get("search_result_status", "not_run")))
            # A gap is a real missed candidate only when the exhaustive arm both
            # searched the gene AND found a hit; an exhaustive no_hit search is
            # not evidence that budget gating hid anything.
            gap = "yes" if estatus == "executed" and ehit == "hit" and hstatus != "executed" else "no"
            gating_rows.append(
                {
                    "gene_id": gene_id,
                    "cluster_id": str(hrow.get("cluster_id", erow.get("cluster_id", ""))),
                    "hierarchical_query_status": hstatus,
                    "exhaustive_query_status": estatus,
                    "hierarchical_hit_presence": hhit,
                    "exhaustive_hit_presence": ehit,
                    "gating_gap": gap,
                }
            )
        gating_rows.sort(
            key=lambda entry: (str(inventory_by_gene.get(entry["gene_id"], {}).get("chrom", "")), int(str(inventory_by_gene.get(entry["gene_id"], {}).get("start", "")) or 0), entry["gene_id"])
        )
        gating_text = GATING_HEADER_COMMENT + "\n" + pd.DataFrame(gating_rows, columns=GATING_COLUMNS).to_csv(sep="\t", index=False, lineterminator="\n")
        (out / "budget_gating_comparison.tsv").write_text(gating_text, encoding="utf-8", newline="\n")
        gating_added_false_negative_count = sum(1 for entry in gating_rows if entry["gating_gap"] == "yes")

    budgets = {}
    for name in LEDGER_LIMIT_FIELDS:
        declared = limits[name]
        remaining = frozen["remaining"][name]
        consumed = None if declared is None else declared - remaining
        budgets[name] = {"declared": declared, "consumed": consumed, "remaining": remaining, "exhausted": frozen["exhausted"][name]}

    outputs = {
        "tier0_inventory_tsv": "tier0_inventory.tsv",
        "tier0_manifest_json": "tier0_manifest.json",
        "tier1_clusters_tsv": "tier1_clusters.tsv",
        "tier1_manifest_json": "tier1_manifest.json",
        "tier2_plan_tsv": "tier2_plan.tsv",
        "tier2_manifest_json": "tier2_manifest.json",
        "query_file_map_tsv": "queries/query_file_map.tsv",
        "budget_handoff_tsv": "budget_handoff.tsv",
        "budget_manifest_json": "budget_manifest.json",
        "benchmark_resource_usage_json": "benchmark_resource_usage.json",
        "deterministic_top20_tsv": "deterministic_top20.tsv",
        "llm_assisted_top20_tsv": "llm_assisted_top20.tsv",
        "portfolio_manifest_json": "portfolio_manifest.json",
    }
    if llm_review_manifest is not None:
        outputs["llm_review_tsv"] = "llm_review.tsv"
        outputs["llm_review_manifest_json"] = "llm_review_manifest.json"
        outputs["llm_review_critic_tsv"] = "llm_review_critic.tsv"
    if diversity_result is not None:
        outputs["proposer_diversity_audit_tsv"] = "proposer_diversity_audit.tsv"
    if exhaustive:
        outputs["budget_gating_comparison_tsv"] = "budget_gating_comparison.tsv"
        outputs["exhaustive_benchmark_resource_usage_json"] = "exhaustive/benchmark_resource_usage.json"

    manifest: dict[str, object] = {
        "schema_version": "1.0",
        "case_id": case_id,
        "status": status,
        "llm_review_status": llm_review_status,
        "budget_profile": profile_name,
        "budgets": budgets,
        "exhausted_budgets": account.exhausted_budgets(),
        "unprocessed": {
            "deferred_query_count": queries_deferred_budget,
            "unreviewed_candidate_count": unreviewed_candidate_count,
            "truncation_points": plan["truncation_points"],
        },
        "cache": cache.stats(),
        "outputs": outputs,
        "exit_code": 0,
    }
    if account.heavy_stages():
        manifest["heavy_stages"] = account.heavy_stages()
        manifest["heavy_stage_invocation_count"] = len(account.heavy_stages())
    if peak_memory_note is not None:
        manifest["peak_memory_note"] = peak_memory_note
    if exhaustive:
        manifest["gating_added_false_negative_count"] = gating_added_false_negative_count
        if protein_db is None:
            manifest["gating_note"] = "exhaustive_arm_provider_missing"
    _write_json(out / "budget_manifest.json", manifest)
    return manifest
