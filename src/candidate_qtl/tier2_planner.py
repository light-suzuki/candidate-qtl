"""Deterministic Tier 2 sequence-search plan (issue #375, deterministic budget engine).

Tier 2 turns the Tier 0 inventory row dicts and the Tier 1 cluster rows into
one row per candidate gene: which queries are planned (representative first,
then condition-triggered expansion, defer-cluster representatives, rescue
escalations), in which deterministic order they execute, and what the budget
account allowed. The plan is a pure function of its inputs plus the budget
ledger: no wall-clock values, no RNG, no timestamps.

Bulk mode: when the hierarchical plan would query every gene anyway
(``queries_bulk <= queries_hierarchical``, i.e. the expansion triggers bought
nothing), the planner emits the simpler bulk plan where every gene is a
representative with exactly one query. The original design text wrote
``queries_bulk < queries_hierarchical`` here, but that comparison can never
hold -- hierarchical queries are always a subset of the gene set, so bulk
could never be selected. The equality case is included so bulk mode is
reachable exactly when it gains nothing over hierarchical.

Execution is cache-first: a content-addressed cache hit reuses the stored hit
rows and their original measured CPU cost (so the ledger stays deterministic
across cached reruns) and bypasses every budget gate -- a warm cache makes a
zero-budget rerun fully reproducible. A miss falls through to the fail-closed
gates, in order: sequence budget, declared wall-time budget, declared peak
memory (pre-run RSS), a per-query subprocess timeout derived from the wall
budget, wall-time charge-back after the run, peak-memory re-measurement, and
finally the sequence-budget charge. The first gate failure marks
``budget_stopped`` and every remaining planned row becomes ``deferred_budget``
(never silently truncated); a finished query whose wall charge is refused
keeps its executed result. The wall-time charge (both the post-run charge and
the charge on a provider timeout) is always the *measured elapsed* wall time
(``perf_counter`` around the call), never the backend's self-reported
``cpu_seconds`` -- a custom/test backend that misreports its cost cannot
bypass the wall hard budget this way; ``cpu_seconds`` itself is still
recorded separately as the non-budgeted ``sequence_search_cpu_seconds``
usage. Peak memory is a hard limit, not a consumable budget: every RSS
measurement (pre-query, post-query, and the final fallback measurement when
no query ran) is recorded via :meth:`~candidate_qtl.budget_account.
BudgetAccount.observe_peak`, which tracks the observed *maximum* (never
summed, never consumed) and marks ``peak_memory_limit_mb`` exhausted on an
overrun, so ``usage.peak_memory_mb`` and ``exhausted.peak_memory_limit_mb``
in the ledger now agree with the ``peak_memory`` truncation points recorded
in ``truncation_points`` (previously the ledger stayed at 0/false while the
run had actually stopped for memory). A provider timeout
(:class:`candidate_qtl.core.ExternalDependencyError`) records the row's
``search_result_status`` as ``provider_timeout`` and its terminal
``query_status`` as ``provider_failed`` (with ``deferred_reason
="provider_timeout"``), and, when a wall-time limit is declared, charges the
measured elapsed wall time (so provider hangs cannot bypass the wall budget);
without a declared wall-time limit no charge is made and the loop continues.
``deferred_plan`` is reserved vocabulary: no current path assigns it (P3
defer-cluster representatives are executed; unselected members are
``not_planned`` and budget-stop rows are ``deferred_budget``). When no search
provider is available (``protein_db`` missing or the backend unavailable)
planned rows are ``not_evaluated`` and nothing is consumed.
"""

from __future__ import annotations

import hashlib
import json
import os
from collections import Counter
from pathlib import Path
from time import perf_counter

import pandas as pd

from .budget_account import BudgetAccount
from .core import (
    DEFAULT_EXTERNAL_TOOL_TIMEOUT,
    TOOL_TIMEOUT_ENV,
    ExternalDependencyError,
    _env_timeout,
    runtime_path,
    safe_id,
    sha256,
)

TIER2_COLUMNS = (
    "gene_id",
    "cluster_id",
    "plan_step",
    "query_role",
    "expansion_condition",
    "query_status",
    "deferred_reason",
    "queries_planned",
    "cost_estimate",
    "ambiguity_status",
    "rescue_escalated",
    "a_b_discordance",
    "user_flagged",
    "search_result_status",
)

PLAN_MODES = frozenset({"hierarchical", "bulk"})
QUERY_ROLES = frozenset({"representative", "expansion_member", "not_queried"})
QUERY_STATUSES = frozenset({"planned", "executed", "deferred_budget", "deferred_plan", "not_planned", "provider_failed"})
# ``deferred_plan`` is reserved vocabulary: no current path assigns it (P3
# defer-cluster representatives are executed; unselected members are
# ``not_planned`` and budget-stop rows are ``deferred_budget``).
# ``provider_failed`` is terminal: a row whose search raised
# ExternalDependencyError (search_result_status "provider_timeout") is
# counted here instead of being left dangling in the non-terminal "planned"
# state. Every planned row that reached the provider therefore lands in
# exactly one of executed / deferred_budget / deferred_plan / not_planned /
# provider_failed. Rows that never reached the provider at all (no protein
# db, backend unavailable, no protein sequence for the gene, or a backend
# RuntimeError) are still left at "planned" and are identified by
# ``search_result_status == "not_evaluated"``; giving those a terminal
# status too is deliberately out of scope here (it would change the
# vocabulary for every provider-missing run) and is tracked as follow-up.
EXPANSION_CONDITIONS = frozenset({"representative_first", "ambiguity_status", "rescue_escalation", "a_b_discordance", "user_flagged"})
BULK_MAX_GENES = 20
PER_QUERY_COST = 1.0

_ALLOCATION_RANK = {"priority": 0, "defer": 1, "rescue": 2}
_CACHE_TOOL = "blastp"
_CACHE_PARAMS = {"max_target_seqs": 5, "evalue": 10}

# Windows reserved device names: a file named e.g. ``con.fa`` cannot be created
# on win32 (OSError), so a query stem whose base name (case-insensitive, with
# trailing dots/spaces stripped) matches one of these falls back to "gene".
_WIN32_RESERVED_STEMS = frozenset({"con", "nul", "aux", "prn"}) | {f"com{i}" for i in range(1, 10)} | {f"lpt{i}" for i in range(1, 10)}


def db_fingerprint(protein_db: str | Path | None) -> str | None:
    """Content fingerprint of a proteome database for cache keys.

    ``protein_db`` may be a single existing file (fingerprinted alone), a
    directory (every direct file child sorted by name), or a BLAST prefix
    whose prefix file is absent but whose ``name.*`` siblings exist (the glob
    sorted by name). Each entry records name, size, and the full-content
    SHA-256 of the file; the digest of the sorted JSON is the fingerprint. The
    paths themselves never enter the fingerprint, so moving the database does
    not invalidate cached hits. Nothing is fingerprinted (None is returned)
    when the path does not exist, the directory is empty, or no ``name.*``
    siblings exist; the caller treats that as a provider gap. Reading a
    multi-GB database in full on every run is the required cost of never
    reusing another database's cache entries; persisting DB hashes between
    runs is out of scope.
    """
    if protein_db is None:
        return None
    path = runtime_path(protein_db)
    if path.is_file():
        files = [path]
    elif path.is_dir():
        files = sorted((entry for entry in path.iterdir() if entry.is_file()), key=lambda entry: entry.name)
        if not files:
            return None
    else:
        files = sorted((entry for entry in path.parent.glob(path.name + ".*") if entry.is_file()), key=lambda entry: entry.name)
        if not files:
            return None
    entries = []
    for entry in files:
        entries.append(
            {
                "name": entry.name,
                "size_bytes": entry.stat().st_size,
                "sha256": sha256(entry),
            }
        )
    payload = json.dumps(entries, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _query_stem(gene_id: str, used: set[str]) -> str:
    """Deterministic filesystem-safe stem for a query FASTA; unique per gene_id.

    ``used`` tracks stems already handed out so colliding gene_ids (after
    ``core.safe_id``) get a deterministic 8-hex suffix. A stem whose base name
    is a win32 reserved device name (``con``/``nul``/``aux``/``prn``/``com1..9``/
    ``lpt1..9``, case-insensitive, trailing dots/spaces stripped) falls back to
    ``"gene"`` so query files can be created on Windows too; the usual
    collision-suffix logic then applies. Deterministic only when called in
    sorted gene_id order.
    """
    stem = safe_id(gene_id)
    if stem in {"", ".", ".."}:
        stem = "gene"
    if stem.rstrip(" .").lower() in _WIN32_RESERVED_STEMS:
        stem = "gene"
    if stem in used:
        stem = f"{stem}-{hashlib.sha256(gene_id.encode('utf-8')).hexdigest()[:8]}"
    used.add(stem)
    return stem


def _safe_query_paths(out: Path, gene_ids: list[str]) -> dict[str, Path]:
    """Map every gene_id to its query FASTA path inside ``<out>/queries``."""
    used: set[str] = set()
    stems = {gene_id: _query_stem(gene_id, used) for gene_id in sorted(gene_ids)}
    queries_dir = (out / "queries").resolve()
    paths: dict[str, Path] = {}
    for gene_id in gene_ids:
        path = out / "queries" / f"{stems[gene_id]}.fa"
        if not path.resolve().is_relative_to(queries_dir):
            raise ValueError(f"query path escapes the queries directory: {path}")
        paths[gene_id] = path
    return paths


def _write_query_file_map(out: Path, stems: dict[str, str]) -> Path:
    """Write ``queries/query_file_map.tsv``; payload_file is the search payload for the stem."""
    rows = [
        {
            "gene_id": gene_id,
            "stem": stems[gene_id],
            "query_file": f"queries/{stems[gene_id]}.fa",
            "payload_file": f"queries/{stems[gene_id]}.search_payload.json",
        }
        for gene_id in sorted(stems)
    ]
    path = out / "queries" / "query_file_map.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=("gene_id", "stem", "query_file", "payload_file")).to_csv(path, sep="\t", index=False, lineterminator="\n")
    return path


def _current_rss_mb() -> tuple[float | None, str | None]:
    """Current process RSS in MiB via /proc/self/statm; (None, "not_obtainable") on failure."""
    try:
        statm = Path("/proc/self/statm").read_text(encoding="ascii").split()
        rss_pages = int(statm[1])
        return round(rss_pages * os.sysconf("SC_PAGE_SIZE") / 1048576.0, 3), None
    except (OSError, ValueError, IndexError):
        return None, "not_obtainable"


def _effective_query_timeout(wall_remaining: float | int | None) -> float:
    """Per-query subprocess timeout: the default tool timeout, capped by the wall budget."""
    default = _env_timeout(TOOL_TIMEOUT_ENV, DEFAULT_EXTERNAL_TOOL_TIMEOUT)
    if wall_remaining is None:
        return default
    return min(default, float(wall_remaining))


def _cell_tokens(value: object) -> list[str]:
    return [token.strip() for token in str(value).split(";") if token.strip()]


def _ambiguity_status(row: dict[str, object]) -> str:
    conflict_flags = str(row.get("conflict_flags", "")).strip()
    domain_ids = {token for token in _cell_tokens(row.get("domain_hits", ""))}
    orthogroup_ids = {token for token in _cell_tokens(row.get("orthogroup", ""))}
    if conflict_flags or len(domain_ids) >= 2:
        return "ambiguous_annotation"
    if len(orthogroup_ids) >= 2:
        return "ambiguous_homology"
    return "not_ambiguous"


def _rescue_escalated(row: dict[str, object], allocation: str, ambiguity: str) -> str:
    if allocation != "rescue":
        return "no"
    if ambiguity != "not_ambiguous":
        return "yes"
    for cell in ("sv_evidence", "pav_evidence", "promoter_evidence"):
        if str(row.get(cell, "")).strip():
            return "yes"
    return "no"


def _hit_rows_from_payload(payload_path: Path) -> tuple[list[dict[str, object]], float]:
    payload = json.loads(payload_path.read_text(encoding="utf-8"))
    return list(payload.get("hits", [])), float(payload.get("cpu_seconds", 0.0))


def _payload_record(out: Path, stem: str, hits: list[dict[str, object]], cpu_seconds: float) -> Path:
    record = out / "queries" / f"{stem}.search_payload.json"
    record.write_text(
        json.dumps({"hits": hits, "cpu_seconds": cpu_seconds}, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return record


def plan_sequence_search(
    inventory_rows: list[dict[str, object]],
    cluster_rows: list[dict[str, object]],
    account: BudgetAccount,
    out_dir: str | Path,
    case_id: str,
    user_flagged: frozenset[str] | set[str] = frozenset(),
    protein_db: str | Path | None = None,
    blastp_bin: str = "blastp",
    search_backend: object | None = None,
    exhaustive_override: bool = False,
    proteins: dict[str, str] | None = None,
) -> dict[str, object]:
    """Plan and execute the budgeted sequence searches for one case.

    ``inventory_rows`` are the Tier 0 row dicts, ``cluster_rows`` the Tier 1
    cluster row dicts. ``account`` is the already-declared
    :class:`~candidate_qtl.budget_account.BudgetAccount`. ``search_backend``
    follows the :class:`~candidate_qtl.budget_engine.BlastpBackend` protocol
    (``available`` bool, optional ``cache`` ContentCache, and ``search(query_fasta,
    db, bin, timeout=None) -> (hits, cpu_seconds)``); when omitted a
    ``BlastpBackend`` is constructed lazily. ``proteins`` maps gene_id to the
    protein sequence used for the query FASTA (one file per gene under
    ``<out-dir>/queries``, named by the deterministic safe stem).

    Writes ``tier2_plan.tsv``, ``queries/query_file_map.tsv``, and
    ``tier2_manifest.json``. The returned manifest additionally carries
    ``truncation_points`` (the budget-stop ledger consumed by
    :func:`candidate_qtl.budget_engine.run_budget_engine`) and
    ``provider_errors`` (per-query provider failures); neither key is written
    into ``tier2_manifest.json``, which records only the deterministic
    ``provider_error_count``.
    """
    out = runtime_path(out_dir)
    flag_set = set(user_flagged)
    cluster_by_gene: dict[str, dict[str, str]] = {}
    for cluster in cluster_rows:
        allocation = str(cluster.get("allocation", ""))
        for member in _cell_tokens(cluster.get("member_gene_ids", "")):
            cluster_by_gene[member] = {
                "cluster_id": str(cluster.get("cluster_id", "")),
                "allocation": allocation,
                "rescue_bucket": str(cluster.get("rescue_bucket", "")),
                "representative_gene_id": str(cluster.get("representative_gene_id", "")),
            }

    gene_rows: list[dict[str, object]] = []
    for inventory in inventory_rows:
        gene_id = str(inventory.get("gene_id", ""))
        cluster = cluster_by_gene.get(gene_id, {})
        allocation = str(cluster.get("allocation", ""))
        ambiguity = _ambiguity_status(inventory)
        gene_rows.append(
            {
                "gene_id": gene_id,
                "cluster_id": cluster.get("cluster_id", ""),
                "allocation": allocation,
                "allocation_rank": _ALLOCATION_RANK.get(allocation, 3),
                "rescue_bucket": cluster.get("rescue_bucket", ""),
                "is_representative": gene_id == cluster.get("representative_gene_id", ""),
                "ambiguity_status": ambiguity,
                "rescue_escalated": _rescue_escalated(inventory, allocation, ambiguity),
                "a_b_discordance": "not_evaluated",
                "user_flagged": "yes" if gene_id in flag_set else "no",
            }
        )

    plan_rows: list[dict[str, object]] = []

    def add_plan(step: int, row: dict[str, object], role: str, condition: str, planned: int, deferred_reason: str = "") -> None:
        plan_rows.append(
            {
                "gene_id": row["gene_id"],
                "cluster_id": row["cluster_id"],
                "plan_step": step,
                "query_role": role,
                "expansion_condition": condition,
                "query_status": "planned" if planned else "not_planned",
                "deferred_reason": deferred_reason,
                "queries_planned": planned,
                "cost_estimate": str(planned * PER_QUERY_COST),
                "ambiguity_status": row["ambiguity_status"],
                "rescue_escalated": row["rescue_escalated"],
                "a_b_discordance": row["a_b_discordance"],
                "user_flagged": row["user_flagged"],
                "search_result_status": "not_run",
                "allocation_rank": row["allocation_rank"],
            }
        )

    order = sorted(gene_rows, key=lambda row: (row["allocation_rank"], str(row["cluster_id"]), str(row["gene_id"])))
    for row in order:
        allocation = str(row["allocation"])
        if allocation == "priority":
            if row["is_representative"]:
                add_plan(1, row, "representative", "representative_first", 1)
            elif str(row["ambiguity_status"]) != "not_ambiguous" or row["user_flagged"] == "yes":
                condition = "ambiguity_status" if str(row["ambiguity_status"]) != "not_ambiguous" else "user_flagged"
                add_plan(2, row, "expansion_member", condition, 1)
            else:
                add_plan(2, row, "not_queried", "", 0, "expansion_not_triggered")
        elif allocation == "defer":
            if row["is_representative"]:
                add_plan(3, row, "representative", "representative_first", 1)
        elif allocation == "rescue":
            if row["rescue_escalated"] == "yes":
                add_plan(4, row, "representative", "rescue_escalation", 1)
            else:
                add_plan(4, row, "not_queried", "", 0, "rescue_pending_review")
    covered = {str(row["gene_id"]) for row in plan_rows}
    for row in order:
        if str(row["gene_id"]) not in covered:
            add_plan(5, row, "not_queried", "", 0, "not_selected")

    queries_hierarchical = sum(1 for row in plan_rows if row["queries_planned"] == 1)
    queries_bulk = len(plan_rows)
    cost_estimate_hierarchical = queries_hierarchical * PER_QUERY_COST
    cost_estimate_bulk = queries_bulk * PER_QUERY_COST
    sequence_query_limit = account.remaining("sequence_search_budget")
    bulk_selected = (
        queries_bulk <= BULK_MAX_GENES
        and queries_bulk <= queries_hierarchical
        and (sequence_query_limit is None or queries_bulk <= sequence_query_limit)
    )

    plan_mode = "bulk" if bulk_selected else "hierarchical"
    if exhaustive_override:
        plan_mode = "exhaustive"
        cost_estimate_hierarchical = queries_bulk * PER_QUERY_COST
        cost_estimate_bulk = queries_bulk * PER_QUERY_COST
        bulk_selected = False
        for row in plan_rows:
            row["query_role"] = "representative"
            row["expansion_condition"] = "representative_first"
            row["query_status"] = "planned"
            row["deferred_reason"] = ""
            row["queries_planned"] = 1
            row["cost_estimate"] = str(PER_QUERY_COST)
    elif bulk_selected:
        for row in plan_rows:
            row["query_role"] = "representative"
            row["expansion_condition"] = "representative_first"
            row["query_status"] = "planned"
            row["deferred_reason"] = ""
            row["queries_planned"] = 1
            row["cost_estimate"] = str(PER_QUERY_COST)

    if search_backend is None:
        from .budget_engine import BlastpBackend

        search_backend = BlastpBackend(cache=None, bin=blastp_bin)
    cache = getattr(search_backend, "cache", None)
    backend_available = bool(getattr(search_backend, "available", True))
    provider_present = protein_db is not None and backend_available
    search_provider = "blastp" if provider_present else "not_evaluated_provider_missing"

    planned_ids = sorted(str(row["gene_id"]) for row in plan_rows if row["queries_planned"] == 1)
    used: set[str] = set()
    stems = {gene_id: _query_stem(gene_id, used) for gene_id in planned_ids}
    query_paths = _safe_query_paths(out, planned_ids)
    _write_query_file_map(out, stems)

    truncation_points: list[dict[str, object]] = []
    provider_errors: list[dict[str, object]] = []
    peak_memory_enforcement = "not_declared"
    peak_rss_obtained = False
    if not provider_present:
        for row in plan_rows:
            if row["queries_planned"] == 1:
                row["search_result_status"] = "not_evaluated"
    else:
        db_hash = db_fingerprint(protein_db)
        if db_hash is None:
            search_provider = "not_evaluated_provider_missing"
            for row in plan_rows:
                if row["queries_planned"] == 1:
                    row["search_result_status"] = "not_evaluated"
        else:
            db_hashes = {"db": db_hash}
            provider_identity = getattr(search_backend, "resolved_bin", None) or "blastp"
            provider_version = getattr(search_backend, "version", None) or "unknown"
            peak_memory_limit = account.remaining("peak_memory_budget")
            cache_enabled = cache is not None and getattr(cache, "_enabled", True)
            budget_stopped = False
            for row in sorted(plan_rows, key=lambda r: (r["plan_step"], r["allocation_rank"], str(r["cluster_id"]), str(r["gene_id"]))):
                if row["queries_planned"] != 1:
                    continue
                if budget_stopped:
                    row["query_status"] = "deferred_budget"
                    row["deferred_reason"] = "budget_exhausted"
                    continue
                gene_id = str(row["gene_id"])
                sequence = (proteins or {}).get(gene_id)
                if not sequence:
                    row["search_result_status"] = "not_evaluated"
                    continue
                query_fasta = query_paths[gene_id]
                query_fasta.parent.mkdir(parents=True, exist_ok=True)
                query_fasta.write_text(f">{gene_id}\n{sequence}\n", encoding="utf-8", newline="\n")
                input_hashes = {"query": sha256(query_fasta), **db_hashes}
                hits: list[dict[str, object]] | None = None
                cpu_seconds = 0.0
                # [CACHE-FIRST] a cache hit bypasses every budget gate (mission D).
                if cache_enabled:
                    hit, payload_path = cache.get(input_hashes, _CACHE_TOOL, provider_identity, provider_version, _CACHE_PARAMS)
                    if hit:
                        hits, cpu_seconds = _hit_rows_from_payload(payload_path)
                        account.add_usage("cache_hit_count", 1)
                        # Reuse the originally measured CPU cost (stored in the
                        # payload) so the ledger records the same cpu_seconds on a
                        # warm rerun as on the cold run that populated the cache.
                        account.add_usage("sequence_search_cpu_seconds", cpu_seconds)
                        row["query_status"] = "executed"
                        row["search_result_status"] = "hit" if hits else "no_hit"
                        continue
                # [BUDGET GATE] sequence budget first, then declared wall time.
                seq_remaining = account.remaining("sequence_search_budget")
                wall_remaining = account.remaining("wall_time_budget")
                if seq_remaining is not None and seq_remaining < 1:
                    budget_stopped = True
                    truncation_points.append({"kind": "sequence_query", "plan_step": int(row["plan_step"]), "remaining_before": round(float(seq_remaining), 3)})
                    row["query_status"] = "deferred_budget"
                    row["deferred_reason"] = "budget_exhausted"
                    continue
                if wall_remaining is not None and wall_remaining <= 0:
                    budget_stopped = True
                    truncation_points.append({"kind": "wall_time", "plan_step": int(row["plan_step"]), "remaining_before": round(float(wall_remaining), 3)})
                    row["query_status"] = "deferred_budget"
                    row["deferred_reason"] = "budget_exhausted"
                    continue
                # [MEMORY PRE] peak memory is a hard limit: a process already
                # over the declared limit must not run the next query (the
                # limit is never consumed; an overrun is recorded as a
                # truncation point).
                if peak_memory_limit is not None:
                    rss, _note = _current_rss_mb()
                    if rss is not None:
                        peak_rss_obtained = True
                        account.observe_peak("peak_memory_budget", rss)
                        if rss > peak_memory_limit:
                            budget_stopped = True
                            truncation_points.append({"kind": "peak_memory", "plan_step": int(row["plan_step"]), "remaining_before": round(rss, 3)})
                            row["query_status"] = "deferred_budget"
                            row["deferred_reason"] = "budget_exhausted"
                            continue
                # [TIMEOUT] per-query subprocess timeout derived from the wall budget.
                try:
                    search_started = perf_counter()
                    hits, cpu_seconds = search_backend.search(query_fasta, str(protein_db), blastp_bin, timeout=_effective_query_timeout(account.remaining("wall_time_budget")))
                except ExternalDependencyError as error:
                    row["search_result_status"] = "provider_timeout"
                    row["query_status"] = "provider_failed"
                    row["deferred_reason"] = "provider_timeout"
                    provider_errors.append({"gene_id": gene_id, "kind": error.kind, "message": str(error)})
                    # a provider hang must not bypass a declared wall budget:
                    # charge the actual elapsed time; a refused charge stops
                    # the loop and defers the remaining rows.
                    if wall_remaining is not None:
                        remaining_before_charge = account.remaining("wall_time_budget")
                        if not account.consume("wall_time_budget", perf_counter() - search_started):
                            budget_stopped = True
                            truncation_points.append({"kind": "wall_time", "plan_step": int(row["plan_step"]), "remaining_before": round(float(remaining_before_charge), 3)})
                    continue
                except RuntimeError:
                    row["search_result_status"] = "not_evaluated"
                    continue
                account.add_usage("cache_miss_count", 1)
                # [WALL POST] charge the measured elapsed wall time (never the
                # backend's self-reported cpu_seconds -- a custom/test backend
                # that misreports its cost must not bypass the wall hard
                # budget); a refused charge stops the loop but keeps the
                # finished row executed. remaining_before is captured before
                # the charge is attempted so the truncation point and the
                # exhausted flag agree with each other.
                elapsed_wall = perf_counter() - search_started
                if wall_remaining is not None:
                    remaining_before_wall_charge = account.remaining("wall_time_budget")
                    if not account.consume("wall_time_budget", elapsed_wall):
                        budget_stopped = True
                        truncation_points.append({"kind": "wall_time", "plan_step": int(row["plan_step"]), "remaining_before": round(float(remaining_before_wall_charge), 3)})
                # [MEMORY POST] re-measure after the run against the declared limit.
                if peak_memory_limit is not None:
                    rss_post, _note = _current_rss_mb()
                    if rss_post is not None:
                        peak_rss_obtained = True
                        account.observe_peak("peak_memory_budget", rss_post)
                        if rss_post > peak_memory_limit:
                            budget_stopped = True
                            truncation_points.append({"kind": "peak_memory", "plan_step": int(row["plan_step"]), "remaining_before": round(rss_post, 3)})
                if not account.consume("sequence_search_budget", 1):
                    budget_stopped = True
                    row["query_status"] = "deferred_budget"
                    row["deferred_reason"] = "budget_exhausted"
                    row["search_result_status"] = "not_run"
                    continue
                account.add_usage("searched_protein_count", 1)
                account.add_usage("sequence_search_cpu_seconds", cpu_seconds)
                row["query_status"] = "executed"
                row["search_result_status"] = "hit" if hits else "no_hit"
                if cache_enabled:
                    cache.put(input_hashes, _CACHE_TOOL, provider_identity, provider_version, _CACHE_PARAMS, _payload_record(out, stems[gene_id], hits, cpu_seconds))
            # peak_memory_enforcement derives from the declaration plus whether
            # RSS was obtainable anywhere, never from which loop path ran; the
            # final measurement below keeps all-hit (cache-warm) runs accurate.
            if peak_memory_limit is not None and not peak_rss_obtained:
                rss_final, _note = _current_rss_mb()
                peak_rss_obtained = rss_final is not None
                if rss_final is not None:
                    account.observe_peak("peak_memory_budget", rss_final)
            peak_memory_enforcement = "not_declared" if peak_memory_limit is None else ("enforced" if peak_rss_obtained else "not_obtainable")

    plan_rows.sort(key=lambda r: (r["plan_step"], r["allocation_rank"], str(r["cluster_id"]), str(r["gene_id"])))
    out.mkdir(parents=True, exist_ok=True)
    tsv_rows = [{column: row[column] for column in TIER2_COLUMNS} for row in plan_rows]
    pd.DataFrame(tsv_rows, columns=TIER2_COLUMNS).to_csv(out / "tier2_plan.tsv", sep="\t", index=False, lineterminator="\n")

    queries_planned_count = sum(1 for row in plan_rows if row["queries_planned"] == 1)
    queries_executed = sum(1 for row in plan_rows if row["query_status"] == "executed")
    queries_deferred_budget = sum(1 for row in plan_rows if row["query_status"] == "deferred_budget")
    queries_deferred_plan = sum(1 for row in plan_rows if row["query_status"] == "deferred_plan")
    queries_not_planned = sum(1 for row in plan_rows if row["query_status"] == "not_planned")
    queries_provider_failed = sum(1 for row in plan_rows if row["query_status"] == "provider_failed")
    ambiguity_counts = Counter(str(row["ambiguity_status"]) for row in plan_rows)
    expansion_counts = {
        "ambiguity_status": dict(sorted(ambiguity_counts.items())),
        "rescue_escalation": sum(1 for row in plan_rows if row["rescue_escalated"] == "yes"),
        "user_flagged": sum(1 for row in plan_rows if row["user_flagged"] == "yes"),
        "a_b_discordance": sum(1 for row in plan_rows if row["a_b_discordance"] == "yes"),
    }
    manifest: dict[str, object] = {
        "schema_version": "1.0",
        "case_id": case_id,
        "plan_mode": plan_mode,
        "gene_count": len(plan_rows),
        "cluster_count": len(cluster_rows),
        "queries_planned": queries_planned_count,
        "queries_executed": queries_executed,
        "queries_deferred_budget": queries_deferred_budget,
        "queries_deferred_plan": queries_deferred_plan,
        "queries_not_planned": queries_not_planned,
        "queries_provider_failed": queries_provider_failed,
        "cost_estimate_hierarchical": cost_estimate_hierarchical,
        "cost_estimate_bulk": cost_estimate_bulk,
        "bulk_selected": bulk_selected,
        "expansion_counts": expansion_counts,
        "search_provider": search_provider,
        "peak_memory_enforcement": peak_memory_enforcement,
        "truncation_points": truncation_points,
        "provider_errors": provider_errors,
    }
    written = {
        key: value
        for key, value in manifest.items()
        if key not in {"truncation_points", "provider_errors"}
    }
    written["provider_error_count"] = len(provider_errors)
    (out / "tier2_manifest.json").write_text(json.dumps(written, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest
