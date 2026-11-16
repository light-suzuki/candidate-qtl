"""Batched, budget-gated LLM cluster review + critic + recall-first portfolio (issue #375 T5).

This module is the only place an injected LLM-review backend is called. It is
strictly downstream of the Tier 0/1/2 budget engine (:mod:`candidate_qtl.
budget_engine`) and of the already-shipped ``mechanism-import`` proposer
(:mod:`candidate_qtl.mechanism_payload_adapter`); it never re-derives Tier 0-2
selection and never reimplements mechanism proposal.

Contract (mirrors :mod:`candidate_qtl.tier2_planner`):

* :func:`build_review_batches` selects ONLY Tier1/2 already-selected candidates
  (``review_state`` priority/rescue, or ``ambiguity_status`` ambiguous, or
  ``rescue_escalated``) and packs whole Tier1 clusters -- never single genes --
  into deterministic, content-hashed batches.
* :func:`run_cluster_review` is the budget-gated execution loop: a
  content-addressed cache hit is looked up FIRST, before every budget gate,
  and reuses the stored review rows and their recorded token usage at zero
  charge (mirroring :func:`candidate_qtl.tier2_planner.plan_sequence_search`'s
  cache-first execution). A miss then runs the fail-closed gates in order --
  ``llm_call_budget`` (>= 1 call per *batch*, never per gene), the backend
  call, then ``llm_input_tokens``/``llm_output_tokens`` -- and a refused
  charge stops the loop: every remaining batch becomes an explicit
  ``deferred_budget`` row (or ``not_evaluated_no_budget`` when the declared
  call budget was already non-positive before the very first batch was
  attempted) plus a ``truncation_points`` entry. Nothing is ever silently
  truncated.
* Hallucination rejection is mandatory and unconditional: a
  ``supported_mechanism_id``/``contradicted_mechanism_id`` absent from the
  caller-supplied mechanism id set, or a ``rationale_source_id`` absent from
  the caller-supplied verified source id set, is rejected -- never counted as
  support -- and recorded in ``rejected_claims``. The review artifact
  (``llm_review.tsv``) carries no score/confidence/probability column;
  ``rank_contribution`` is always ``"none"``.
* :func:`build_critic` audits the completed review for fixation, coverage
  holes, contradictions, unreviewed scope, budget exhaustion and
  single-cluster (family-proxy) saturation; a clean run still writes
  ``llm_review_critic.tsv`` with a header and zero rows.
* :func:`audit_proposer_diversity` is a read-only audit over the EXISTING
  ``mechanism_hypotheses.tsv`` produced by ``mechanism-import`` -- it never
  reimplements or mutates the proposer.
* :func:`reduce_portfolio` is the recall-first Top-N reducer. It writes TWO
  artifacts with an identical schema: ``deterministic_top20.tsv`` (a pure
  function of the budget-engine outputs -- no LLM input ever changes it, no
  score/probability/weight column, ``rank_contribution="none"``) and
  ``llm_assisted_top20.tsv`` (the same gene set; LLM review may only reorder
  within the already-selected slate through an explicit recorded
  ``selection_basis`` -- it can never add, remove, or score a gene). With zero
  LLM budget or no review at all, the two files are identical row-for-row.

Every writer in this module follows the house manifest/TSV conventions: JSON
via ``json.dumps(..., ensure_ascii=False, indent=2, sort_keys=True,
allow_nan=False) + "\\n"``, TSV via ``DataFrame.to_csv(sep="\\t", index=False,
lineterminator="\\n")``, multi-value cells as sorted ``;``-joined tokens.
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from time import perf_counter

import pandas as pd

from .budget_account import BudgetAccount
from .budget_cache import ContentCache
from .core import runtime_path
# The RSS measurement is shared with the Tier 2 loop rather than duplicated, so
# "current RSS" means exactly one thing across the engine. Imported as a module
# (not a bound name) so a test patching tier2_planner._current_rss_mb also
# controls this stage.
from . import tier2_planner as _tier2

# ---------------------------------------------------------------------------
# Shared column contracts
# ---------------------------------------------------------------------------
REVIEW_COLUMNS = (
    "batch_id", "cluster_id", "gene_id", "review_status",
    "supported_mechanism_ids", "contradicted_mechanism_ids", "rationale_source_ids",
    "dependency_cluster_id", "uncertainty", "missing_evidence", "next_experiment",
    "provider", "model", "prompt_sha256", "input_sha256", "output_sha256",
    # The complete request identity behind the row: prompt + provider/model/
    # version + the validation allowlists in force. A resume only reuses a row
    # whose fingerprint still matches, so a changed responses file, a changed
    # provider or a newly restrictive allowlist can never be answered with a
    # stale row (issue #375 T5 audit, P1).
    "request_fingerprint",
    "rank_contribution",
)
REVIEW_STATUSES = frozenset({
    "reviewed", "rejected_unverified", "deferred_budget", "not_evaluated_no_budget",
    "empty_response", "answered_over_budget",
})
REJECTION_REASONS = frozenset({"unknown_mechanism_id", "unverified_source_id"})

CRITIC_COLUMNS = ("finding_kind", "scope", "detail", "affected_ids", "severity", "rank_contribution")
CRITIC_FINDING_KINDS = frozenset({
    "fixation", "coverage_hole", "contradiction", "unreviewed_scope",
    "budget_exhaustion", "single_family_saturation",
})

DIVERSITY_AUDIT_COLUMNS = ("metric", "value", "threshold", "verdict", "detail")
# Named thresholds (module constants, per issue #375 T5): deliberately low bars
# -- this audit exists to catch a proposer that degenerates to one gene/family,
# not to second-guess a genuinely narrow but well-evidenced pathway.
PROPOSER_DIVERSITY_THRESHOLDS = {
    "distinct_mechanism_count_min": 3,
    "distinct_gene_family_count_min": 2,
    "distinct_reference_gene_count_min": 2,
    "distinct_source_id_count_min": 2,
    "top_reference_gene_share_max": 0.5,
}

TOP20_COLUMNS = (
    "slot_index", "gene_id", "cluster_id", "allocation", "rescue_bucket", "slot_kind", "slot_reason",
    "diversity_key", "dependency_cluster_id", "paralog_ambiguity", "route_rescue_state", "llm_status",
    "selection_basis", "rank_contribution",
)
LLM_STATUSES = frozenset({
    "reviewed", "rejected_unverified", "deferred_budget", "not_evaluated_no_budget",
    "empty_response", "answered_over_budget",
})
MAX_PER_DIVERSITY_KEY = 3
# ``single_family_saturation`` needs enough distinct families to be a real
# concentration signal; below this it is arithmetic (with two clusters one of
# them always holds more than half of an odd number of rows).
MIN_CLUSTERS_FOR_SATURATION = 3
RESCUE_BUCKET_ORDER = ("unclassified", "weak_annotation", "missing_model", "sv_pav_regulatory")
ROUTE_ORDER = ("a", "b", "c")

_REVIEW_STATE_RANK = {"priority_review": 0, "rescue_pending_review": 1, "deferred": 2, "not_evaluated": 3}
_REVIEW_TOOL = "llm_cluster_review"
# Statuses that mean the batch already consumed its call and must not be
# re-run (and re-charged) on a resume. ``answered_over_budget`` belongs here
# too: the provider was called and paid for, the token budget simply could not
# afford the answer, so the claims are refused -- but re-calling would pay a
# second time for a response already obtained (issue #375 T5 audit, P2).
_ANSWERED_REVIEW_STATUSES = frozenset({
    "reviewed", "rejected_unverified", "empty_response", "answered_over_budget",
})


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------
def _as_list(value: object) -> list[str]:
    """Normalize a cell/JSON field into a list of stripped, non-empty tokens."""
    if value is None:
        return []
    if isinstance(value, str):
        return [token.strip() for token in value.split(";") if token.strip()]
    if isinstance(value, (list, tuple, set, frozenset)):
        tokens: list[str] = []
        for item in value:
            text = str(item).strip()
            if text:
                tokens.append(text)
        return tokens
    text = str(value).strip()
    return [text] if text else []


def _joined(values) -> str:
    return ";".join(sorted(set(_as_list(values))))


def _write_manifest(path: Path, payload: dict[str, object]) -> None:
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _content_hash(payload: object) -> str:
    serialized = json.dumps(payload, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(serialized).hexdigest()


# ---------------------------------------------------------------------------
# A1: offline, deterministic review backend
# ---------------------------------------------------------------------------
class FileReviewBackend:
    """Offline, deterministic review "provider" reading canned responses from a JSON file.

    There is no network path anywhere in this class. The responses file may be
    ``{batch_id: {"items": [...], "input_tokens": N, "output_tokens": N}}``, a
    ``{"responses": [{"batch_id": ..., "items": [...], ...}, ...]}`` wrapper, or
    a bare list of the same per-batch objects. ``review(batch)`` looks the
    batch up by its ``batch_id`` and returns ``{"items": [...],
    "input_tokens": int, "output_tokens": int}``; a missing ``batch_id`` (or a
    malformed entry) yields an empty-response result -- never an exception,
    never an invented review.
    """

    def __init__(self, responses_path: str | Path, provider: str = "offline-canned", model: str = "canned-review-v1", version: str = "1.0") -> None:
        self.provider = provider
        self.model = model
        self.version = version
        self.cache: ContentCache | None = None
        self._responses = self._load(responses_path)
        self.available = True

    @staticmethod
    def _load(path: str | Path) -> dict[str, dict[str, object]]:
        resolved = runtime_path(path)
        try:
            raw = resolved.read_text(encoding="utf-8")
        except OSError as error:
            raise ValueError(f"llm responses file not found or unreadable: {resolved}: {error}") from error
        try:
            data = json.loads(raw)
        except ValueError as error:
            raise ValueError(f"llm responses file is not valid JSON: {resolved}: {error}") from error
        entries: list[object]
        if isinstance(data, dict) and isinstance(data.get("responses"), list):
            entries = data["responses"]
        elif isinstance(data, list):
            entries = data
        elif isinstance(data, dict):
            responses: dict[str, dict[str, object]] = {}
            for key, value in data.items():
                if isinstance(value, dict):
                    responses[str(key)] = value
            return responses
        else:
            raise ValueError(f"llm responses file must be a JSON object or array: {resolved}")
        responses = {}
        for entry in entries:
            if isinstance(entry, dict) and str(entry.get("batch_id", "")).strip():
                responses[str(entry["batch_id"])] = entry
        return responses

    def review(self, batch: dict[str, object]) -> dict[str, object]:
        """Look up a batch by ``batch_id``; a miss returns an empty-response result."""
        entry = self._responses.get(str(batch.get("batch_id", "")))
        if not isinstance(entry, dict):
            return {"items": [], "input_tokens": 0, "output_tokens": 0}
        items = entry.get("items", [])
        if not isinstance(items, list):
            items = []
        return {
            "items": items,
            "input_tokens": max(0, _safe_int(entry.get("input_tokens", 0))),
            "output_tokens": max(0, _safe_int(entry.get("output_tokens", 0))),
        }


def _safe_int(value: object) -> int:
    try:
        return int(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return 0


# ---------------------------------------------------------------------------
# A2: post-freeze batch selection
# ---------------------------------------------------------------------------
def build_review_batches(
    handoff_rows: list[dict[str, object]],
    plan_rows: list[dict[str, object]],
    cluster_rows: list[dict[str, object]],
    batch_size: int = 5,
) -> list[dict[str, object]]:
    """Select eligible Tier1/2 candidates and pack whole clusters into batches.

    A gene is eligible when its handoff ``review_state`` is
    ``priority_review`` or ``rescue_pending_review``, OR its Tier2 row has
    ``ambiguity_status != "not_ambiguous"``, OR ``rescue_escalated == "yes"``.
    Eligible genes are grouped by ``cluster_id``; whole clusters (never single
    genes) are then packed into batches of at most ``batch_size`` clusters
    each, in deterministic ``cluster_id`` order. Each batch's ``batch_id`` is
    a content hash (``sha256`` hex, first 20 chars, prefixed ``batch_``) of
    its normalized ``{cluster_ids, members}`` payload, so identical input
    always yields identical batch ids.
    """
    plan_by_gene = {str(row.get("gene_id", "")): row for row in plan_rows}
    cluster_ids_by_gene = {str(row.get("gene_id", "")): str(row.get("cluster_id", "")) for row in cluster_rows_flat(cluster_rows)}

    eligible_by_cluster: dict[str, set[str]] = defaultdict(set)
    for row in handoff_rows:
        gene_id = str(row.get("gene_id", ""))
        review_state = str(row.get("review_state", ""))
        plan_row = plan_by_gene.get(gene_id, {})
        ambiguity_status = str(plan_row.get("ambiguity_status", ""))
        rescue_escalated = str(plan_row.get("rescue_escalated", ""))
        eligible = (
            review_state in {"priority_review", "rescue_pending_review"}
            or (ambiguity_status not in {"", "not_ambiguous"})
            or rescue_escalated == "yes"
        )
        if not eligible:
            continue
        cluster_id = str(row.get("cluster_id", "")) or cluster_ids_by_gene.get(gene_id, "")
        eligible_by_cluster[cluster_id].add(gene_id)

    cluster_ids_sorted = sorted(eligible_by_cluster)
    size = max(1, int(batch_size))
    grouped: list[list[str]] = [cluster_ids_sorted[i:i + size] for i in range(0, len(cluster_ids_sorted), size)]

    batches: list[dict[str, object]] = []
    for cluster_ids in grouped:
        members = [
            {"cluster_id": cluster_id, "gene_id": gene_id}
            for cluster_id in cluster_ids
            for gene_id in sorted(eligible_by_cluster[cluster_id])
        ]
        payload = {"cluster_ids": cluster_ids, "members": members}
        batch_id = "batch_" + _content_hash(payload)[:20]
        batches.append({"batch_id": batch_id, "cluster_ids": list(cluster_ids), "members": members})
    return batches


def cluster_rows_flat(cluster_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    """Expand Tier1 cluster rows into one ``{gene_id, cluster_id}`` row per member."""
    flat: list[dict[str, object]] = []
    for cluster in cluster_rows:
        cluster_id = str(cluster.get("cluster_id", ""))
        for gene_id in _as_list(cluster.get("member_gene_ids", "")):
            flat.append({"gene_id": gene_id, "cluster_id": cluster_id})
    return flat


# ---------------------------------------------------------------------------
# A3: budget-gated review execution loop
# ---------------------------------------------------------------------------
def _deferred_row(batch_id: str, member: dict[str, object], status: str, provider: str, model: str, prompt_sha256: str, request_fingerprint: str = "") -> dict[str, object]:
    return {
        "batch_id": batch_id,
        "cluster_id": member["cluster_id"],
        "gene_id": member["gene_id"],
        "review_status": status,
        "supported_mechanism_ids": "",
        "contradicted_mechanism_ids": "",
        "rationale_source_ids": "",
        "dependency_cluster_id": "",
        "uncertainty": "",
        "missing_evidence": "",
        "next_experiment": "",
        "provider": provider,
        "model": model,
        "prompt_sha256": prompt_sha256,
        "input_sha256": "",
        "output_sha256": "",
        "request_fingerprint": request_fingerprint,
        "rank_contribution": "none",
    }


def _reviewed_row(
    batch_id: str, member: dict[str, object], status: str, provider: str, model: str,
    prompt_sha256: str, input_sha256: str, output_sha256: str, item: dict[str, object] | None,
    supported: list[str], contradicted: list[str], rationale: list[str],
    request_fingerprint: str = "",
) -> dict[str, object]:
    item = item or {}
    return {
        "batch_id": batch_id,
        "cluster_id": member["cluster_id"],
        "gene_id": member["gene_id"],
        "review_status": status,
        "supported_mechanism_ids": ";".join(sorted(set(supported))),
        "contradicted_mechanism_ids": ";".join(sorted(set(contradicted))),
        "rationale_source_ids": ";".join(sorted(set(rationale))),
        "dependency_cluster_id": str(item.get("dependency_cluster_id", "")).strip(),
        "uncertainty": str(item.get("uncertainty", "")).strip(),
        "missing_evidence": str(item.get("missing_evidence", "")).strip(),
        "next_experiment": str(item.get("next_experiment", "")).strip(),
        "provider": provider,
        "model": model,
        "prompt_sha256": prompt_sha256,
        "input_sha256": input_sha256,
        "output_sha256": output_sha256,
        "request_fingerprint": request_fingerprint,
        "rank_contribution": "none",
    }


def run_cluster_review(
    batches: list[dict[str, object]],
    account: BudgetAccount,
    out: str | Path,
    case_id: str,
    review_backend: object,
    mechanisms: set[str] | frozenset[str] | None = None,
    allowed_source_ids: set[str] | frozenset[str] | None = None,
    cache: ContentCache | None = None,
    resume: bool = True,
) -> dict[str, object]:
    """Run the budget-gated batched review loop; writes ``llm_review.tsv`` + manifest.

    Per batch: a cache hit (looked up FIRST, before every budget gate) reuses
    the stored rows at zero charge. A miss then runs the fail-closed gates in
    order -- ``llm_call_budget`` (1 call per *batch*), the backend call, then
    ``llm_input_tokens``/``llm_output_tokens`` -- and a refused charge stops
    the loop: every remaining batch (this one included) becomes an explicit
    row (``not_evaluated_no_budget`` when the declared call budget was
    already <= 0 before the very first batch was attempted, else
    ``deferred_budget``) plus one ``truncation_points`` entry. A returned
    ``supported_mechanism_id``/``contradicted_mechanism_id`` outside
    ``mechanisms`` (when supplied), or a ``rationale_source_id`` outside
    ``allowed_source_ids`` (when supplied), is rejected -- ``review_status``
    becomes ``rejected_unverified`` and the claim is recorded in
    ``rejected_claims``; it never becomes support.

    ``resume=True`` skips (no backend call, no charge) any batch that is both
    already answered -- every row on ``reviewed`` / ``rejected_unverified`` /
    ``empty_response`` / ``answered_over_budget`` -- AND whose stored
    ``request_fingerprint`` still equals the current one. The fingerprint
    covers the prompt, the provider/model/version and the validation
    allowlists, so a changed responses file, a changed provider or a newly
    restrictive allowlist re-runs the batch instead of silently reusing a
    stale answer. A batch the budget never reached (``deferred_budget`` /
    ``not_evaluated_no_budget``) is retried; a batch that was answered but
    whose output tokens could not be afforded is recorded
    ``answered_over_budget`` (the claims are refused, but the call was already
    paid for, so resume never buys the same response twice).

    The wall-time and peak-memory hard budgets cover this stage too: a batch
    is not attempted when the declared wall budget is already spent or the
    measured RSS is already over the declared peak limit, and the measured
    elapsed time of each call is charged to the wall budget exactly as the
    Tier 2 loop charges its searches.

    ``counts`` in the manifest are per review ROW (one row per eligible gene)
    and include rows reused from the cache or from a resume; ``batch_counts``
    are per batch (planned / cache_hit / called / resumed).
    """
    out = runtime_path(out)
    out.mkdir(parents=True, exist_ok=True)
    mechanism_ids = set(mechanisms) if mechanisms is not None else None
    source_ids = set(allowed_source_ids) if allowed_source_ids is not None else None

    provider = str(getattr(review_backend, "provider", ""))
    model = str(getattr(review_backend, "model", ""))
    version = str(getattr(review_backend, "version", ""))
    backend_cache = cache if cache is not None else getattr(review_backend, "cache", None)

    # The validation identity travels with every cache key and every stored
    # row: a cached or resumed answer was validated against THESE allowlists,
    # so reusing it under different ones would let a claim that a later,
    # stricter allowlist rejects survive as `reviewed` (T5 audit, P1).
    validation_fingerprint = _content_hash({
        "mechanisms": sorted(mechanism_ids) if mechanism_ids is not None else None,
        "allowed_source_ids": sorted(source_ids) if source_ids is not None else None,
    })

    def request_fingerprint_for(prompt_sha256: str) -> str:
        return _content_hash({
            "prompt": prompt_sha256, "provider": provider, "model": model,
            "version": version, "validation": validation_fingerprint,
        })

    counts = {
        "reviewed": 0, "deferred": 0, "rejected": 0, "empty": 0,
        "not_evaluated_no_budget": 0, "answered_over_budget": 0,
    }
    _COUNT_KEY_FOR_STATUS = {
        "reviewed": "reviewed",
        "rejected_unverified": "rejected",
        "empty_response": "empty",
        "deferred_budget": "deferred",
        "not_evaluated_no_budget": "not_evaluated_no_budget",
        "answered_over_budget": "answered_over_budget",
    }

    def tally(reused_rows: list[dict[str, object]]) -> None:
        """Count reused (cache-hit or resumed) rows exactly like freshly built ones."""
        for row in reused_rows:
            key = _COUNT_KEY_FOR_STATUS.get(str(row.get("review_status", "")))
            if key:
                counts[key] += 1

    review_path = out / "llm_review.tsv"
    rows_by_batch: dict[str, list[dict[str, object]]] = {}
    resumed_batch_ids: set[str] = set()
    if resume and review_path.is_file():
        existing = pd.read_csv(review_path, sep="\t", dtype=str).fillna("").to_dict("records")
        grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
        for row in existing:
            grouped[str(row.get("batch_id", ""))].append(row)
        expected_fingerprints = {
            str(batch["batch_id"]): request_fingerprint_for(
                _content_hash({"batch_id": str(batch["batch_id"]), "cluster_ids": batch["cluster_ids"], "members": batch["members"]})
            )
            for batch in batches
        }
        # A batch is resumable when it already got an answer AND the request
        # behind that answer is still the current one. Only batches the budget
        # never reached are retried, so a resume never re-charges a call that
        # was already paid for -- but a changed responses file, provider or
        # allowlist changes the fingerprint and forces a real re-run.
        for batch_id, rows in grouped.items():
            if not rows or any(row.get("review_status") not in _ANSWERED_REVIEW_STATUSES for row in rows):
                continue
            expected = expected_fingerprints.get(batch_id)
            if expected is None or any(str(row.get("request_fingerprint", "")) != expected for row in rows):
                continue
            rows_by_batch[batch_id] = rows
            resumed_batch_ids.add(batch_id)

    initial_remaining = account.remaining("llm_call_budget")
    zero_from_start = initial_remaining is not None and initial_remaining <= 0

    rows: list[dict[str, object]] = []
    rejected_claims: list[dict[str, object]] = []
    truncation_points: list[dict[str, object]] = []
    batches_cache_hit = 0
    batches_called = 0
    budget_stopped = False
    peak_memory_limit = account.remaining("peak_memory_budget")

    ordered_batches = sorted(batches, key=lambda batch: str(batch.get("batch_id", "")))
    for index, batch in enumerate(ordered_batches):
        batch_id = str(batch["batch_id"])
        members = batch["members"]
        if batch_id in resumed_batch_ids:
            rows.extend(rows_by_batch[batch_id])
            tally(rows_by_batch[batch_id])
            continue

        prompt_payload = {"batch_id": batch_id, "cluster_ids": batch["cluster_ids"], "members": members}
        prompt_sha256 = _content_hash(prompt_payload)
        fingerprint = request_fingerprint_for(prompt_sha256)

        def defer(status: str, reason_key: str) -> None:
            for member in members:
                rows.append(_deferred_row(batch_id, member, status, provider, model, prompt_sha256, fingerprint))
                counts[reason_key] += 1

        deferred_status = "not_evaluated_no_budget" if zero_from_start else "deferred_budget"
        deferred_key = "not_evaluated_no_budget" if zero_from_start else "deferred"

        # [CACHE-FIRST] a cache hit bypasses every budget gate, always -- even
        # after a prior batch already stopped the loop (invariant 6). The
        # validation fingerprint is part of the key, so a hit was validated
        # against exactly the allowlists in force now.
        if backend_cache is not None:
            input_hashes = {"prompt": prompt_sha256}
            params = {"provider": provider, "model": model, "validation": validation_fingerprint}
            hit, payload_path = backend_cache.get(input_hashes, _REVIEW_TOOL, provider, model, params)
            if hit:
                stored = json.loads(Path(payload_path).read_text(encoding="utf-8"))
                stored_rows = stored.get("rows", [])
                account.add_usage("cache_hit_count", 1)
                rows.extend(stored_rows)
                batches_cache_hit += 1
                tally(stored_rows)
                continue

        if budget_stopped:
            defer(deferred_status, deferred_key)
            continue

        # [BUDGET GATE 1] one llm call per batch, never per gene.
        remaining_before_call = account.remaining("llm_call_budget")
        if not account.consume("llm_call_budget", 1):
            budget_stopped = True
            truncation_points.append({"kind": "llm_call", "batch_index": index, "remaining_before": remaining_before_call})
            defer(deferred_status, deferred_key)
            continue

        # [WALL PRE] the declared wall budget covers this stage too: a spent
        # wall budget stops the review exactly as it stops a Tier 2 search.
        wall_remaining = account.remaining("wall_time_budget")
        if wall_remaining is not None and wall_remaining <= 0:
            budget_stopped = True
            truncation_points.append({"kind": "wall_time", "batch_index": index, "remaining_before": round(float(wall_remaining), 3)})
            defer(deferred_status, deferred_key)
            continue

        # [MEMORY PRE] peak memory is a hard limit, never consumed: a process
        # already over the declared limit must not run the next batch.
        if peak_memory_limit is not None:
            rss, _note = _tier2._current_rss_mb()
            if rss is not None:
                account.observe_peak("peak_memory_budget", rss)
                if rss > peak_memory_limit:
                    budget_stopped = True
                    truncation_points.append({"kind": "peak_memory", "batch_index": index, "remaining_before": round(rss, 3)})
                    defer(deferred_status, deferred_key)
                    continue

        call_started = perf_counter()
        response = review_backend.review(prompt_payload)
        elapsed_wall = perf_counter() - call_started
        batches_called += 1
        # [WALL POST] charge the measured elapsed, never a self-reported cost.
        if wall_remaining is not None:
            remaining_before_wall = account.remaining("wall_time_budget")
            if not account.consume("wall_time_budget", elapsed_wall):
                budget_stopped = True
                truncation_points.append({"kind": "wall_time", "batch_index": index, "remaining_before": round(float(remaining_before_wall), 3)})
        # [MEMORY POST] re-measure against the declared limit.
        if peak_memory_limit is not None:
            rss_post, _note = _tier2._current_rss_mb()
            if rss_post is not None:
                account.observe_peak("peak_memory_budget", rss_post)
                if rss_post > peak_memory_limit:
                    budget_stopped = True
                    truncation_points.append({"kind": "peak_memory", "batch_index": index, "remaining_before": round(rss_post, 3)})
        items = response.get("items", []) if isinstance(response, dict) else []
        input_tokens = _safe_int(response.get("input_tokens", 0)) if isinstance(response, dict) else 0
        output_tokens = _safe_int(response.get("output_tokens", 0)) if isinstance(response, dict) else 0

        # [BUDGET GATE 2/3] input then output tokens against the shared
        # token_limit. The call has already happened and been paid for, so a
        # refused token charge does NOT produce a retryable deferred row: the
        # claims are refused (never admitted unpaid) but the batch is recorded
        # ``answered_over_budget`` so a resume does not buy the same response
        # a second time (T5 audit, P2).
        remaining_before_input = account.remaining("token_budget")
        if not account.consume("llm_input_tokens", input_tokens):
            budget_stopped = True
            truncation_points.append({"kind": "token", "batch_index": index, "remaining_before": remaining_before_input})
            defer("answered_over_budget", "answered_over_budget")
            continue
        remaining_before_output = account.remaining("token_budget")
        if not account.consume("llm_output_tokens", output_tokens):
            budget_stopped = True
            truncation_points.append({"kind": "token", "batch_index": index, "remaining_before": remaining_before_output})
            defer("answered_over_budget", "answered_over_budget")
            continue

        output_sha256 = _content_hash(items)
        items_by_gene: dict[str, dict[str, object]] = {}
        for item in items:
            if isinstance(item, dict) and str(item.get("gene_id", "")).strip():
                items_by_gene[str(item["gene_id"]).strip()] = item

        batch_rows: list[dict[str, object]] = []
        for member in members:
            gene_id = str(member["gene_id"])
            item = items_by_gene.get(gene_id)
            if item is None:
                row = _reviewed_row(batch_id, member, "empty_response", provider, model, prompt_sha256, prompt_sha256, output_sha256, None, [], [], [], fingerprint)
                counts["empty"] += 1
                batch_rows.append(row)
                continue
            supported = _as_list(item.get("supported_mechanism_ids", item.get("supported_mechanism_id")))
            contradicted = _as_list(item.get("contradicted_mechanism_ids", item.get("contradicted_mechanism_id")))
            rationale = _as_list(item.get("rationale_source_ids"))
            reasons: list[str] = []
            if mechanism_ids is not None:
                bad = [m for m in (*supported, *contradicted) if m not in mechanism_ids]
                if bad:
                    reasons.append("unknown_mechanism_id")
            if source_ids is not None:
                bad_sources = [s for s in rationale if s not in source_ids]
                if bad_sources:
                    reasons.append("unverified_source_id")
            if reasons:
                for reason in sorted(set(reasons)):
                    rejected_claims.append({"batch_id": batch_id, "gene_id": gene_id, "reason": reason})
                row = _reviewed_row(batch_id, member, "rejected_unverified", provider, model, prompt_sha256, prompt_sha256, output_sha256, item, [], [], [], fingerprint)
                counts["rejected"] += 1
            elif not supported and not contradicted and not any(
                str(item.get(field, "")).strip() for field in ("uncertainty", "missing_evidence", "next_experiment")
            ):
                row = _reviewed_row(batch_id, member, "empty_response", provider, model, prompt_sha256, prompt_sha256, output_sha256, item, [], [], [], fingerprint)
                counts["empty"] += 1
            else:
                row = _reviewed_row(batch_id, member, "reviewed", provider, model, prompt_sha256, prompt_sha256, output_sha256, item, supported, contradicted, rationale, fingerprint)
                counts["reviewed"] += 1
            batch_rows.append(row)

        rows.extend(batch_rows)
        if backend_cache is not None:
            import tempfile

            fd, tmp_path = tempfile.mkstemp(prefix="llm_review_payload_", suffix=".json", dir=out)
            try:
                Path(tmp_path).write_text(json.dumps({"rows": batch_rows}, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
                backend_cache.put(input_hashes, _REVIEW_TOOL, provider, model, params, tmp_path)
            finally:
                try:
                    Path(tmp_path).unlink()
                except FileNotFoundError:
                    pass
                try:
                    import os

                    os.close(fd)
                except OSError:
                    pass

    pd.DataFrame(rows, columns=REVIEW_COLUMNS).to_csv(review_path, sep="\t", index=False, lineterminator="\n")

    total_planned = sum(len(batch["members"]) for batch in batches)
    manifest = {
        "schema_version": "1.0",
        "case_id": case_id,
        "provider": provider,
        "model": model,
        "version": str(getattr(review_backend, "version", "")),
        # Row counts cover every row written to llm_review.tsv, including rows
        # reused from the cache or from a resume, so the manifest and the TSV
        # never disagree (T5 audit, P2).
        "counts": {
            "planned": total_planned,
            "reviewed": counts["reviewed"],
            "deferred": counts["deferred"],
            "rejected": counts["rejected"],
            "empty": counts["empty"],
            "not_evaluated_no_budget": counts["not_evaluated_no_budget"],
            "answered_over_budget": counts["answered_over_budget"],
        },
        "batch_counts": {
            "planned": len(batches),
            "cache_hit": batches_cache_hit,
            "called": batches_called,
            "resumed": len(resumed_batch_ids),
        },
        "resume": {"resume_enabled": bool(resume), "resumed_batch_count": len(resumed_batch_ids)},
        "rejected_claims": rejected_claims,
        "truncation_points": truncation_points,
        "rank_contribution": "none",
    }
    _write_manifest(out / "llm_review_manifest.json", manifest)
    return manifest


# ---------------------------------------------------------------------------
# A4: critic over the completed review
# ---------------------------------------------------------------------------
def build_critic(
    review_rows: list[dict[str, object]],
    batches: list[dict[str, object]],
    handoff_rows: list[dict[str, object]],
    cluster_rows: list[dict[str, object]],
    account: BudgetAccount,
    out: str | Path,
    case_id: str,
) -> dict[str, object]:
    """Audit a completed review; writes ``llm_review_critic.tsv`` and merges counts into the review manifest.

    Always writes the file, header-only when clean. Findings never carry a
    score; ``rank_contribution`` is always ``"none"``.
    """
    out = runtime_path(out)
    out.mkdir(parents=True, exist_ok=True)

    reviewed = [row for row in review_rows if str(row.get("review_status", "")) == "reviewed"]
    total_reviewed = len(reviewed)
    findings: list[dict[str, object]] = []

    if total_reviewed > 0:
        support_counter: Counter[str] = Counter()
        for row in reviewed:
            for mech in _as_list(row.get("supported_mechanism_ids", "")):
                support_counter[mech] += 1
        for mech, count in sorted(support_counter.items()):
            if count > total_reviewed / 2:
                findings.append({
                    "finding_kind": "fixation", "scope": "mechanism",
                    "detail": f"mechanism {mech} is supported in {count}/{total_reviewed} reviewed row(s)",
                    "affected_ids": mech, "severity": "high", "rank_contribution": "none",
                })

    eligible_cluster_ids = sorted({str(member["cluster_id"]) for batch in batches for member in batch["members"]})
    reviewed_clusters = {str(row.get("cluster_id", "")) for row in reviewed}
    for cluster_id in eligible_cluster_ids:
        if cluster_id not in reviewed_clusters:
            findings.append({
                "finding_kind": "coverage_hole", "scope": "cluster",
                "detail": f"eligible cluster {cluster_id} has zero reviewed rows",
                "affected_ids": cluster_id, "severity": "medium", "rank_contribution": "none",
            })

    rescue_bucket_by_cluster = {str(row.get("cluster_id", "")): str(row.get("rescue_bucket", "")) for row in cluster_rows}
    eligible_buckets = sorted({
        rescue_bucket_by_cluster.get(cluster_id, "")
        for cluster_id in eligible_cluster_ids
        if rescue_bucket_by_cluster.get(cluster_id, "")
    })
    reviewed_buckets = {rescue_bucket_by_cluster.get(cluster_id, "") for cluster_id in reviewed_clusters}
    for bucket in eligible_buckets:
        if bucket not in reviewed_buckets:
            findings.append({
                "finding_kind": "coverage_hole", "scope": "rescue_bucket",
                "detail": f"rescue_bucket {bucket} has zero reviewed rows",
                "affected_ids": bucket, "severity": "medium", "rank_contribution": "none",
            })

    by_gene: dict[str, dict[str, set[str]]] = defaultdict(lambda: {"supported": set(), "contradicted": set()})
    for row in reviewed:
        gene_id = str(row.get("gene_id", ""))
        by_gene[gene_id]["supported"] |= set(_as_list(row.get("supported_mechanism_ids", "")))
        by_gene[gene_id]["contradicted"] |= set(_as_list(row.get("contradicted_mechanism_ids", "")))
    for gene_id in sorted(by_gene):
        overlap = sorted(by_gene[gene_id]["supported"] & by_gene[gene_id]["contradicted"])
        for mech in overlap:
            findings.append({
                "finding_kind": "contradiction", "scope": "gene",
                "detail": f"gene {gene_id} both supports and contradicts mechanism {mech}",
                "affected_ids": f"{gene_id}:{mech}", "severity": "high", "rank_contribution": "none",
            })

    eligible_genes = sorted({str(member["gene_id"]) for batch in batches for member in batch["members"]})
    reviewed_genes = {str(row.get("gene_id", "")) for row in reviewed}
    unreviewed = [gene_id for gene_id in eligible_genes if gene_id not in reviewed_genes]
    if unreviewed:
        findings.append({
            "finding_kind": "unreviewed_scope", "scope": "portfolio",
            "detail": f"{len(unreviewed)} eligible gene(s) were never reviewed",
            "affected_ids": ";".join(unreviewed), "severity": "medium", "rank_contribution": "none",
        })

    frozen = account.freeze()
    exhausted_budgets = sorted(
        name for name in ("llm_call_limit", "token_limit")
        if frozen["limits"].get(name) is not None and frozen["exhausted"].get(name)
    )
    if exhausted_budgets:
        findings.append({
            "finding_kind": "budget_exhaustion", "scope": "portfolio",
            "detail": f"budget(s) exhausted: {', '.join(exhausted_budgets)}",
            "affected_ids": ";".join(exhausted_budgets), "severity": "medium", "rank_contribution": "none",
        })

    # Saturation is a statement about *concentration across the families that
    # were actually available*: with fewer than MIN_CLUSTERS_FOR_SATURATION
    # eligible clusters, one cluster holding more than half of the reviewed
    # rows is arithmetic, not a finding, so the check would fire on every
    # small run and train reviewers to ignore it.
    if total_reviewed > 0 and len(eligible_cluster_ids) >= MIN_CLUSTERS_FOR_SATURATION:
        cluster_counter = Counter(str(row.get("cluster_id", "")) for row in reviewed)
        for cluster_id, count in sorted(cluster_counter.items()):
            if count > total_reviewed / 2:
                findings.append({
                    "finding_kind": "single_family_saturation", "scope": "cluster",
                    "detail": f"cluster {cluster_id} holds {count}/{total_reviewed} reviewed row(s)",
                    "affected_ids": cluster_id, "severity": "medium", "rank_contribution": "none",
                })

    findings.sort(key=lambda finding: (str(finding["finding_kind"]), str(finding["scope"]), str(finding["affected_ids"])))
    pd.DataFrame(findings, columns=CRITIC_COLUMNS).to_csv(out / "llm_review_critic.tsv", sep="\t", index=False, lineterminator="\n")

    manifest_path = out / "llm_review_manifest.json"
    manifest: dict[str, object] = {}
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["case_id"] = case_id
    manifest["critic"] = {
        "finding_count": len(findings),
        "finding_kind_counts": dict(sorted(Counter(str(f["finding_kind"]) for f in findings).items())),
    }
    manifest["rank_contribution"] = "none"
    _write_manifest(manifest_path, manifest)
    return {"findings": findings, "manifest": manifest}


# ---------------------------------------------------------------------------
# A5: pre-interval proposer diversity audit (read-only over mechanism-import output)
# ---------------------------------------------------------------------------
def audit_proposer_diversity(mechanism_rows: list[dict[str, object]], out: str | Path, case_id: str) -> dict[str, object]:
    """Diversity audit over ``mechanism_hypotheses.tsv`` rows; never mutates the input.

    Covers distinct mechanism (``hypothesis_id``), gene-family
    (``gene_family``), reference-gene (``reference_gene_name``) and source-id
    (``source_ids``) counts, plus the share of the single most frequent
    reference gene. Always writes ``proposer_diversity_audit.tsv``.
    """
    out = runtime_path(out)
    out.mkdir(parents=True, exist_ok=True)
    rows = list(mechanism_rows)

    mechanism_ids = sorted({str(row.get("hypothesis_id", "")).strip() for row in rows if str(row.get("hypothesis_id", "")).strip()})
    gene_families = sorted({str(row.get("gene_family", "")).strip() for row in rows if str(row.get("gene_family", "")).strip()})
    reference_genes = [str(row.get("reference_gene_name", "")).strip() for row in rows if str(row.get("reference_gene_name", "")).strip()]
    distinct_reference_genes = sorted(set(reference_genes))
    source_ids: set[str] = set()
    for row in rows:
        source_ids |= set(_as_list(row.get("source_ids", "")))

    reference_counter = Counter(reference_genes)
    top_gene, top_count = reference_counter.most_common(1)[0] if reference_counter else ("", 0)
    top_share = (top_count / len(reference_genes)) if reference_genes else 0.0

    def _row(metric: str, value: object, threshold: object, verdict: str, detail: str) -> dict[str, object]:
        return {"metric": metric, "value": str(value), "threshold": str(threshold), "verdict": verdict, "detail": detail}

    entries: list[dict[str, object]] = []
    insufficient = False

    mech_min = PROPOSER_DIVERSITY_THRESHOLDS["distinct_mechanism_count_min"]
    mech_pass = len(mechanism_ids) >= mech_min
    insufficient = insufficient or not mech_pass
    entries.append(_row("distinct_mechanism_count", len(mechanism_ids), mech_min, "pass" if mech_pass else "insufficient_diversity", f"{len(mechanism_ids)} distinct hypothesis_id value(s)"))

    fam_min = PROPOSER_DIVERSITY_THRESHOLDS["distinct_gene_family_count_min"]
    fam_pass = len(gene_families) >= fam_min
    insufficient = insufficient or not fam_pass
    entries.append(_row("distinct_gene_family_count", len(gene_families), fam_min, "pass" if fam_pass else "insufficient_diversity", f"{len(gene_families)} distinct gene_family value(s)"))

    ref_min = PROPOSER_DIVERSITY_THRESHOLDS["distinct_reference_gene_count_min"]
    ref_pass = len(distinct_reference_genes) >= ref_min
    insufficient = insufficient or not ref_pass
    entries.append(_row("distinct_reference_gene_count", len(distinct_reference_genes), ref_min, "pass" if ref_pass else "insufficient_diversity", f"{len(distinct_reference_genes)} distinct reference_gene_name value(s)"))

    src_min = PROPOSER_DIVERSITY_THRESHOLDS["distinct_source_id_count_min"]
    src_pass = len(source_ids) >= src_min
    insufficient = insufficient or not src_pass
    entries.append(_row("distinct_source_id_count", len(source_ids), src_min, "pass" if src_pass else "insufficient_diversity", f"{len(source_ids)} distinct source_ids token(s)"))

    share_max = PROPOSER_DIVERSITY_THRESHOLDS["top_reference_gene_share_max"]
    fixation = top_share > share_max
    detail = (
        f"most frequent reference gene {top_gene!r} appears {top_count}/{len(reference_genes)} time(s) (share={top_share:.4f})"
        if reference_genes else "no reference_gene_name values proposed"
    )
    entries.append(_row("top_reference_gene_share", f"{top_share:.4f}", share_max, "single_famous_gene_fixation" if fixation else "pass", detail))

    overall_verdict = "single_famous_gene_fixation" if fixation else ("insufficient_diversity" if insufficient else "pass")

    pd.DataFrame(entries, columns=DIVERSITY_AUDIT_COLUMNS).to_csv(out / "proposer_diversity_audit.tsv", sep="\t", index=False, lineterminator="\n")
    return {"rows": entries, "overall_verdict": overall_verdict, "case_id": case_id}


# ---------------------------------------------------------------------------
# A6: recall-first Top-N portfolio reducer
# ---------------------------------------------------------------------------
def _select_slots(
    order: list[str],
    by_gene: dict[str, dict[str, object]],
    bridge_by_gene: dict[str, list[dict[str, object]]],
    portfolio_size: int,
) -> tuple[dict[str, dict[str, str]], list[dict[str, object]]]:
    selected: dict[str, dict[str, str]] = {}
    unfilled_slots: list[dict[str, object]] = []

    def admit(gene_id: str, slot_kind: str, slot_reason: str) -> bool:
        if gene_id in selected or len(selected) >= portfolio_size:
            return False
        selected[gene_id] = {"slot_kind": slot_kind, "slot_reason": slot_reason}
        return True

    # 1. rescue_slot: at least one slot per non-empty rescue_bucket.
    for bucket in RESCUE_BUCKET_ORDER:
        bucket_candidates = [gene_id for gene_id in order if str(by_gene[gene_id]["rescue_bucket"]) == bucket]
        if not bucket_candidates:
            continue
        if not any(admit(gene_id, "rescue_slot", f"rescue_bucket:{bucket}") for gene_id in bucket_candidates):
            unfilled_slots.append({"slot_kind": "rescue_slot", "reason": "portfolio_full", "detail": bucket})

    # 2. route_unique_rescue: at least one slot per route that uniquely contributes a candidate.
    if bridge_by_gene:
        route_to_genes: dict[str, list[str]] = defaultdict(list)
        for gene_id in order:
            for bridge_row in bridge_by_gene.get(gene_id, []):
                route = str(bridge_row.get("route", ""))
                if route:
                    route_to_genes[route].append(gene_id)
        for route in ROUTE_ORDER:
            genes_for_route = route_to_genes.get(route, [])
            if not genes_for_route:
                continue
            unique_genes = [
                gene_id for gene_id in genes_for_route
                if {str(row.get("route", "")) for row in bridge_by_gene.get(gene_id, [])} == {route}
            ]
            candidates_for_route = unique_genes or genes_for_route
            if not any(admit(gene_id, "route_unique_rescue", f"route:{route}") for gene_id in candidates_for_route):
                unfilled_slots.append({"slot_kind": "route_unique_rescue", "reason": "portfolio_full", "detail": route})
    else:
        unfilled_slots.append({"slot_kind": "route_unique_rescue", "reason": "no_route_information_available", "detail": ""})

    # 3. dependency_dedup: no second member of the same dependency_cluster_id
    #    is admitted until every distinct dependency cluster has had a chance.
    dependency_groups: dict[str, list[str]] = defaultdict(list)
    for gene_id in order:
        dep = str(by_gene[gene_id]["dependency_cluster_id"])
        if dep:
            dependency_groups[dep].append(gene_id)
    represented = {str(by_gene[gene_id]["dependency_cluster_id"]) for gene_id in selected if str(by_gene[gene_id]["dependency_cluster_id"])}
    for dep in sorted(dependency_groups):
        if dep in represented:
            continue
        if any(gene_id in selected for gene_id in dependency_groups[dep]):
            represented.add(dep)
            continue
        for gene_id in dependency_groups[dep]:
            if admit(gene_id, "dependency_dedup", f"dependency_cluster:{dep}"):
                represented.add(dep)
                break

    # 4. family_diversity: at most MAX_PER_DIVERSITY_KEY per diversity_key (cluster_id)
    #    while other keys still have candidates.
    diversity_counts: dict[str, int] = defaultdict(int)
    for gene_id in selected:
        diversity_counts[str(by_gene[gene_id]["cluster_id"])] += 1
    pool_by_key: dict[str, list[str]] = defaultdict(list)
    for gene_id in order:
        if gene_id not in selected:
            pool_by_key[str(by_gene[gene_id]["cluster_id"])].append(gene_id)
    progressed = True
    while progressed and len(selected) < portfolio_size:
        progressed = False
        for key in sorted(pool_by_key):
            queue = pool_by_key[key]
            if not queue or diversity_counts[key] >= MAX_PER_DIVERSITY_KEY:
                continue
            gene_id = queue.pop(0)
            if admit(gene_id, "family_diversity", f"diversity_key:{key}"):
                diversity_counts[key] += 1
                progressed = True
            if len(selected) >= portfolio_size:
                break
    if len(selected) < portfolio_size:
        # every key hit its cap but slots remain and no key still has
        # "spare" supply below the cap -- overflow deterministically rather
        # than leave recall-relevant candidates unslotted.
        for key in sorted(pool_by_key):
            for gene_id in list(pool_by_key[key]):
                if len(selected) >= portfolio_size:
                    break
                if admit(gene_id, "family_diversity", f"diversity_key_overflow:{key}"):
                    diversity_counts[key] += 1

    # 5. deterministic_fill: remaining slots in base order.
    for gene_id in order:
        if len(selected) >= portfolio_size:
            break
        admit(gene_id, "deterministic_fill", "base_order")

    return selected, unfilled_slots


def _apply_paralog_guard(
    selected_order: list[str],
    by_gene: dict[str, dict[str, object]],
    is_paralog_ambiguous,
) -> dict[str, str]:
    diversity_all: dict[str, list[str]] = defaultdict(list)
    for gene_id, entry in by_gene.items():
        diversity_all[str(entry["cluster_id"])].append(gene_id)

    reasons: dict[str, str] = {}
    for gene_id in selected_order:
        if not is_paralog_ambiguous(gene_id):
            continue
        key = str(by_gene[gene_id]["cluster_id"])
        co_candidates = [other for other in diversity_all[key] if other != gene_id]
        reasons[gene_id] = "paralog_ambiguous_co_listed" if co_candidates else "paralog_ambiguous_no_independent_support"

    if selected_order and reasons.get(selected_order[0]) == "paralog_ambiguous_no_independent_support" and len(selected_order) > 1:
        selected_order[0], selected_order[1] = selected_order[1], selected_order[0]
    return reasons


def _llm_assisted_order(det_order: list[str], review_by_gene: dict[str, dict[str, object]]) -> tuple[list[str], dict[str, str]]:
    """Reorder only through an explicit, recorded ``selection_basis``; never adds/removes/scores."""
    selection_basis: dict[str, str] = {}
    if not review_by_gene:
        for gene_id in det_order:
            selection_basis[gene_id] = "deterministic_order_unchanged"
        return list(det_order), selection_basis
    promoted: list[str] = []
    rest: list[str] = []
    for gene_id in det_order:
        review = review_by_gene.get(gene_id)
        if review and str(review.get("review_status", "")) == "reviewed" and str(review.get("supported_mechanism_ids", "")).strip():
            promoted.append(gene_id)
            selection_basis[gene_id] = "llm_supported_promotion"
        else:
            rest.append(gene_id)
            selection_basis[gene_id] = "deterministic_order_unchanged"
    return promoted + rest, selection_basis


def reduce_portfolio(
    handoff_rows: list[dict[str, object]],
    plan_rows: list[dict[str, object]],
    cluster_rows: list[dict[str, object]],
    out: str | Path,
    case_id: str,
    review_rows: list[dict[str, object]] | None = None,
    bridge_rows: list[dict[str, object]] | None = None,
    portfolio_size: int = 20,
) -> dict[str, object]:
    """Recall-first Top-N reducer; writes ``deterministic_top20.tsv`` + ``llm_assisted_top20.tsv``.

    ``deterministic_top20.tsv`` is a pure function of ``handoff_rows`` /
    ``plan_rows`` / ``cluster_rows`` / ``bridge_rows`` only -- ``review_rows``
    never changes it, so it never carries a score/probability/weight column
    and ``llm_status``/``selection_basis`` are fixed placeholders on every
    row. ``llm_assisted_top20.tsv`` shares the exact same gene set and only
    reorders within it via :func:`_llm_assisted_order`; with no
    ``review_rows`` the two files are identical row-for-row (every
    ``llm_status`` is ``not_evaluated_no_budget``).
    """
    out = runtime_path(out)
    out.mkdir(parents=True, exist_ok=True)

    plan_by_gene = {str(row.get("gene_id", "")): row for row in plan_rows}
    cluster_by_id = {str(row.get("cluster_id", "")): row for row in cluster_rows}
    bridge_by_gene: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in (bridge_rows or []):
        bridge_by_gene[str(row.get("gene_id", ""))].append(row)
    review_by_gene = {str(row.get("gene_id", "")): row for row in (review_rows or [])}

    base_order = sorted(
        handoff_rows,
        key=lambda row: (_REVIEW_STATE_RANK.get(str(row.get("review_state", "")), 4), str(row.get("gene_id", ""))),
    )

    by_gene: dict[str, dict[str, object]] = {}
    order: list[str] = []
    for row in base_order:
        gene_id = str(row.get("gene_id", ""))
        if gene_id in by_gene:
            continue
        cluster_id = str(row.get("cluster_id", ""))
        cluster = cluster_by_id.get(cluster_id, {})
        plan_row = plan_by_gene.get(gene_id, {})
        dependency_ids = _as_list(cluster.get("dependency_cluster_ids", ""))
        by_gene[gene_id] = {
            "gene_id": gene_id,
            "cluster_id": cluster_id,
            "allocation": str(row.get("allocation", "")),
            "rescue_bucket": str(row.get("rescue_bucket", "")),
            "ambiguity_status": str(plan_row.get("ambiguity_status", "")),
            "dependency_cluster_id": dependency_ids[0] if dependency_ids else "",
        }
        order.append(gene_id)

    def is_paralog_ambiguous(gene_id: str) -> bool:
        entry = by_gene[gene_id]
        if str(entry["ambiguity_status"]) == "ambiguous_homology":
            return True
        return any(str(row.get("ambiguity_status", "")) == "paralog_ambiguous" for row in bridge_by_gene.get(gene_id, []))

    def route_state(gene_id: str) -> str:
        rows_for_gene = bridge_by_gene.get(gene_id)
        if not rows_for_gene:
            return "not_evaluated"
        routes = sorted({str(row.get("route", "")) for row in rows_for_gene if str(row.get("route", ""))})
        return ";".join(routes) if routes else "not_evaluated"

    portfolio_size = max(0, int(portfolio_size))
    selected, unfilled_slots = _select_slots(order, by_gene, bridge_by_gene, portfolio_size)
    selected_order = list(selected.keys())
    paralog_reasons = _apply_paralog_guard(selected_order, by_gene, is_paralog_ambiguous)

    def build_rows(order_list: list[str], llm_status_of, selection_basis_of) -> list[dict[str, object]]:
        rows: list[dict[str, object]] = []
        for index, gene_id in enumerate(order_list, start=1):
            entry = by_gene[gene_id]
            slot_info = selected[gene_id]
            slot_reason = paralog_reasons.get(gene_id, slot_info["slot_reason"])
            rows.append({
                "slot_index": index,
                "gene_id": gene_id,
                "cluster_id": entry["cluster_id"],
                "allocation": entry["allocation"],
                "rescue_bucket": entry["rescue_bucket"],
                "slot_kind": slot_info["slot_kind"],
                "slot_reason": slot_reason,
                "diversity_key": entry["cluster_id"],
                "dependency_cluster_id": entry["dependency_cluster_id"],
                "paralog_ambiguity": "yes" if is_paralog_ambiguous(gene_id) else "no",
                "route_rescue_state": route_state(gene_id),
                "llm_status": llm_status_of(gene_id),
                "selection_basis": selection_basis_of(gene_id),
                "rank_contribution": "none",
            })
        return rows

    deterministic_rows = build_rows(selected_order, lambda _gene_id: "not_evaluated_no_budget", lambda _gene_id: "deterministic_only")
    llm_order, selection_basis_map = _llm_assisted_order(selected_order, review_by_gene)
    llm_rows = build_rows(
        llm_order,
        lambda gene_id: str(review_by_gene.get(gene_id, {}).get("review_status", "") or "not_evaluated_no_budget"),
        lambda gene_id: selection_basis_map.get(gene_id, "deterministic_order_unchanged"),
    )

    pd.DataFrame(deterministic_rows, columns=TOP20_COLUMNS).to_csv(out / "deterministic_top20.tsv", sep="\t", index=False, lineterminator="\n")
    pd.DataFrame(llm_rows, columns=TOP20_COLUMNS).to_csv(out / "llm_assisted_top20.tsv", sep="\t", index=False, lineterminator="\n")

    manifest = {
        "schema_version": "1.0",
        "case_id": case_id,
        "portfolio_size": portfolio_size,
        "candidate_count": len(order),
        "deterministic_count": len(deterministic_rows),
        "llm_assisted_count": len(llm_rows),
        "unfilled_slots": unfilled_slots,
        "rank_contribution": "none",
    }
    _write_manifest(out / "portfolio_manifest.json", manifest)
    return {"deterministic_rows": deterministic_rows, "llm_assisted_rows": llm_rows, "manifest": manifest}
