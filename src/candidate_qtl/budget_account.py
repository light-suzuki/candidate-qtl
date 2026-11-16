"""Deterministic budget engine: per-run spend ledger.

A :class:`BudgetAccount` is the single source of truth for what a deterministic
run may spend: LLM calls, tokens, sequence-search queries, and wall time are
*consumable* budgets charged via :meth:`consume` (fail-closed: it never
overdraws a declared limit, a refused consume charges nothing and marks the
budget exhausted). Peak memory is not consumable -- it is a hard limit
observed via :meth:`observe_peak`, which records the highest amount seen
(a max, never a sum) and marks the limit exhausted when an observation
exceeds it; ``consume`` rejects peak-memory kinds outright so a caller cannot
accidentally sum peaks. :meth:`BudgetAccount.freeze` produces a plain
JSON-serializable snapshot for manifests so two runs with the same inputs and
limits freeze identically (peak-memory ``remaining`` is clamped at 0, never
negative, since an overrun is not a debt).
"""

from __future__ import annotations

from contextlib import contextmanager
import math
import time
import tracemalloc

LEDGER_LIMIT_FIELDS = (
    "llm_call_limit",
    "token_limit",
    "sequence_query_limit",
    "wall_time_limit_seconds",
    "peak_memory_limit_mb",
)

LEDGER_USAGE_FIELDS = (
    "llm_call_count",
    "llm_input_tokens",
    "llm_output_tokens",
    "sequence_query_count",
    "searched_protein_count",
    "sequence_search_cpu_seconds",
    "wall_time_seconds",
    "peak_memory_mb",
    "cache_hit_count",
    "cache_miss_count",
    "reviewed_candidate_count",
    "unreviewed_candidate_count",
    "external_retrieval_count",
    "heavy_stage_invocation_count",
    "cpu_time_seconds",
)

LEDGER_INTEGER_FIELDS = frozenset({
    "llm_call_count", "llm_input_tokens", "llm_output_tokens", "sequence_query_count",
    "searched_protein_count", "cache_hit_count", "cache_miss_count",
    "reviewed_candidate_count", "unreviewed_candidate_count",
    "external_retrieval_count", "heavy_stage_invocation_count",
})

LEDGER_STATUSES = frozenset({"completed", "budget_exhausted", "failed"})

CONCEPTUAL_TO_LEDGER = {
    "llm_call_budget": "llm_call_limit",
    "token_budget": "token_limit",
    "sequence_search_budget": "sequence_query_limit",
    "wall_time_budget": "wall_time_limit_seconds",
    "peak_memory_budget": "peak_memory_limit_mb",
}

BUDGET_PROFILES = {
    "zero": {"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 0, "wall_time_limit_seconds": 0, "peak_memory_limit_mb": 0},
    "student_cpu": {"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 50, "wall_time_limit_seconds": 3600, "peak_memory_limit_mb": 8192},
    "exhaustive": {"llm_call_limit": None, "token_limit": None, "sequence_query_limit": None, "wall_time_limit_seconds": None, "peak_memory_limit_mb": None},
    # issue #375's own wording for this profile name uses a hyphen (unlike
    # every other key here, which uses an underscore); kept as-is so
    # ``--budget-profile student-budget`` matches the issue text literally.
    "student-budget": {"llm_call_limit": 20, "token_limit": 60000, "sequence_query_limit": 50, "wall_time_limit_seconds": 3600, "peak_memory_limit_mb": 8192},
    # Issue #434 profiles
    "light": {"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 20, "wall_time_limit_seconds": 1800, "peak_memory_limit_mb": 4096},
    "standard": {"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 50, "wall_time_limit_seconds": 3600, "peak_memory_limit_mb": 8192},
    "heavy": {"llm_call_limit": 50, "token_limit": 150000, "sequence_query_limit": 200, "wall_time_limit_seconds": 7200, "peak_memory_limit_mb": 16384},
    "diagnostic": {"llm_call_limit": None, "token_limit": None, "sequence_query_limit": None, "wall_time_limit_seconds": None, "peak_memory_limit_mb": None},
}

# Every accepted ``consume`` kind maps to (ledger limit name, usage field).
# ``token_limit`` is special: its charged usage is the *combined* input+output
# tokens, so its usage field is None and both token fields share the same limit.
_BUDGET_KIND_MAP = {
    "llm_call_budget": ("llm_call_limit", "llm_call_count"),
    "llm_call_limit": ("llm_call_limit", "llm_call_count"),
    "llm_call_count": ("llm_call_limit", "llm_call_count"),
    "token_budget": ("token_limit", None),
    "token_limit": ("token_limit", None),
    "llm_input_tokens": ("token_limit", "llm_input_tokens"),
    "llm_output_tokens": ("token_limit", "llm_output_tokens"),
    "sequence_search_budget": ("sequence_query_limit", "sequence_query_count"),
    "sequence_query_limit": ("sequence_query_limit", "sequence_query_count"),
    "sequence_query_count": ("sequence_query_limit", "sequence_query_count"),
    "wall_time_budget": ("wall_time_limit_seconds", "wall_time_seconds"),
    "wall_time_limit_seconds": ("wall_time_limit_seconds", "wall_time_seconds"),
    "wall_time_seconds": ("wall_time_limit_seconds", "wall_time_seconds"),
    "peak_memory_budget": ("peak_memory_limit_mb", "peak_memory_mb"),
    "peak_memory_limit_mb": ("peak_memory_limit_mb", "peak_memory_mb"),
    "peak_memory_mb": ("peak_memory_limit_mb", "peak_memory_mb"),
}

_BUDGETED_USAGE_FIELDS = frozenset({
    "llm_call_count", "llm_input_tokens", "llm_output_tokens", "sequence_query_count",
    "wall_time_seconds", "peak_memory_mb",
})


class BudgetAccount:
    """Tracks declared limits, charged usage, and remaining budget."""

    def __init__(self, limits: dict[str, object] | None = None) -> None:
        self._limits: dict[str, int | None] = {}
        self._declared: bool = False
        self._exhausted_flags: set[str] = set()
        self._usage: dict[str, int | float] = {
            field: 0.0 if field not in LEDGER_INTEGER_FIELDS else 0 for field in LEDGER_USAGE_FIELDS
        }
        self._heavy_stages: list[dict[str, object]] = []
        if limits is not None:
            self.declare(limits)

    def declare(self, limits: dict[str, object]) -> None:
        """Declare limits; each key is a ledger or conceptual name, value int or None (unlimited)."""
        if self._declared:
            raise ValueError("budget already declared")
        resolved: dict[str, int | None] = {}
        for raw_key, raw_value in limits.items():
            if raw_key in CONCEPTUAL_TO_LEDGER:
                key = CONCEPTUAL_TO_LEDGER[raw_key]
            elif raw_key in LEDGER_LIMIT_FIELDS:
                key = raw_key
            else:
                raise ValueError(f"unknown budget limit key {raw_key!r}")
            if raw_value is None:
                resolved[key] = None
            elif isinstance(raw_value, bool) or not isinstance(raw_value, int) or raw_value < 0:
                raise ValueError(f"budget limit {key} must be a non-negative int or None, got {raw_value!r}")
            else:
                resolved[key] = raw_value
        self._limits.update(resolved)
        self._declared = True

    def _resolve(self, kind: str) -> tuple[str, str | None]:
        if kind in _BUDGET_KIND_MAP:
            return _BUDGET_KIND_MAP[kind]
        raise ValueError(f"unknown budget kind {kind!r}")

    @staticmethod
    def _check_amount(amount: object, field: str | None) -> int | float:
        if isinstance(amount, bool) or not isinstance(amount, (int, float)):
            raise ValueError(f"budget amount must be a number, got {amount!r}")
        if not math.isfinite(amount):
            raise ValueError(f"budget amount must be finite, got {amount!r}")
        if amount < 0:
            raise ValueError(f"budget amount must be non-negative, got {amount!r}")
        if field in LEDGER_INTEGER_FIELDS and not isinstance(amount, int):
            raise ValueError(f"amount for integer usage field {field} must be an int, got {amount!r}")
        return amount

    def _usage_for(self, limit_name: str) -> int | float:
        if limit_name == "token_limit":
            return self._usage["llm_input_tokens"] + self._usage["llm_output_tokens"]
        return self._usage[_BUDGET_KIND_MAP[limit_name][1]]

    def _remaining_for(self, limit_name: str) -> int | float | None:
        limit = self._limits.get(limit_name)
        if limit is None:
            return None
        remaining = limit - self._usage_for(limit_name)
        # Peaks are hard limits, never consumed, so an observed overrun must
        # not surface as a negative "remaining"; every other limit is
        # fail-closed via consume() and never overdraws in the first place.
        if limit_name == "peak_memory_limit_mb" and remaining < 0:
            return 0
        return remaining

    def consume(self, kind: str, amount: object) -> bool:
        """Charge ``amount`` to a budget; returns False (no charge) when it would overdraw."""
        limit_name, usage_field = self._resolve(kind)
        if limit_name == "peak_memory_limit_mb":
            raise ValueError(f"{kind!r} is a peak, not a consumable budget; use observe_peak() instead")
        target = "llm_input_tokens" if limit_name == "token_limit" and usage_field is None else usage_field
        checked = self._check_amount(amount, target)
        remaining = self._remaining_for(limit_name)
        if remaining is not None and remaining < checked:
            self._exhausted_flags.add(limit_name)
            return False
        self._usage[target] += checked
        return True

    def observe_peak(self, kind: str, amount: object) -> bool:
        """Record a peak (max) observation; never consumed, never summed.

        Only peak-memory kinds (``peak_memory_budget`` / ``peak_memory_limit_mb``
        / ``peak_memory_mb``) are meaningful here. Unlike :meth:`consume`,
        the usage field is set to ``max(current usage, amount)`` -- repeated
        observations record the highest value seen, not a running total --
        so the ledger reflects the true observed peak even when it exceeds
        the declared limit. Returns False and marks the limit exhausted when
        the declared limit is exceeded (``amount > limit``); True otherwise
        (including when no limit is declared).
        """
        limit_name, usage_field = self._resolve(kind)
        checked = self._check_amount(amount, usage_field)
        self._usage[usage_field] = max(self._usage[usage_field], checked)
        limit = self._limits.get(limit_name)
        if limit is not None and checked > limit:
            self._exhausted_flags.add(limit_name)
            return False
        return True

    def add_usage(self, field: str, amount: object) -> None:
        """Record non-budgeted usage; budgeted fields must be charged via consume()."""
        if field not in LEDGER_USAGE_FIELDS:
            raise ValueError(f"unknown usage field {field!r}")
        if field in _BUDGETED_USAGE_FIELDS:
            raise ValueError(f"usage field {field!r} is budgeted; charge it via consume() instead")
        checked = self._check_amount(amount, field)
        self._usage[field] += checked

    def remaining(self, kind: str) -> int | float | None:
        """Remaining budget for a kind; None means unlimited."""
        limit_name, _usage_field = self._resolve(kind)
        return self._remaining_for(limit_name)

    def exhausted(self, kind: str) -> bool:
        """True when a budget is spent (limit reached or an overdraw was refused)."""
        limit_name, _usage_field = self._resolve(kind)
        limit = self._limits.get(limit_name)
        if limit is None:
            return False
        remaining = self._remaining_for(limit_name)
        return limit_name in self._exhausted_flags or remaining <= 0

    def exhausted_budgets(self) -> list[str]:
        """Sorted ledger-limit names that are fully spent."""
        return sorted(
            name for name in LEDGER_LIMIT_FIELDS
            if self._limits.get(name) is not None
            and (name in self._exhausted_flags or self._remaining_for(name) <= 0)
        )

    def record_heavy_stage(
        self,
        stage: str,
        invocation_reason: str,
        *,
        wall_time_seconds: float = 0.0,
        cpu_seconds: float = 0.0,
        peak_memory_mb: float = 0.0,
        metadata: dict[str, object] | None = None,
    ) -> dict[str, object]:
        """Record an explicit heavy-stage invocation event with its decision/invocation reason."""
        stage_name = str(stage).strip()
        if not stage_name:
            raise ValueError("stage name must not be empty")
        reason = str(invocation_reason).strip()
        if not reason:
            raise ValueError(f"heavy stage {stage_name!r} requires a non-empty invocation_reason")
        record = {
            "stage": stage_name,
            "invocation_reason": reason,
            "wall_time_seconds": round(float(wall_time_seconds), 6),
            "cpu_seconds": round(float(cpu_seconds), 6),
            "peak_memory_mb": round(float(peak_memory_mb), 3),
            "metadata": metadata or {},
        }
        self._heavy_stages.append(record)
        self._usage["heavy_stage_invocation_count"] += 1
        return record

    def heavy_stages(self) -> list[dict[str, object]]:
        """Return recorded heavy-stage invocation events."""
        return [dict(entry) for entry in self._heavy_stages]

    def freeze(self) -> dict[str, object]:
        """Plain JSON-serializable snapshot for manifests; limits None where unlimited."""
        frozen: dict[str, object] = {
            "limits": {name: self._limits.get(name) for name in LEDGER_LIMIT_FIELDS},
            "usage": {name: self._usage[name] for name in LEDGER_USAGE_FIELDS},
            "remaining": {name: self._remaining_for(name) for name in LEDGER_LIMIT_FIELDS},
            "exhausted": {name: self.exhausted(name) for name in LEDGER_LIMIT_FIELDS},
        }
        if self._heavy_stages:
            frozen["heavy_stages"] = [dict(entry) for entry in self._heavy_stages]
        return frozen


@contextmanager
def stage_telemetry(
    account: BudgetAccount,
    stage_name: str,
    *,
    invocation_reason: str = "",
    is_heavy: bool = False,
    metadata: dict[str, object] | None = None,
):
    """Context manager for portable stage runtime telemetry.

    Measures:
    * wall-clock time via :func:`time.perf_counter`;
    * process CPU time via :func:`time.process_time`;
    * Python memory allocation peak via :mod:`tracemalloc`.

    Peak memory is observed and folded into ``account.observe_peak("peak_memory_mb", ...)``.
    When ``is_heavy=True``, records the heavy-stage invocation event with ``invocation_reason``.
    """
    tracemalloc_started_here = False
    if not tracemalloc.is_tracing():
        tracemalloc.start()
        tracemalloc_started_here = True
    tracemalloc.reset_peak()

    start_wall = time.perf_counter()
    start_cpu = time.process_time()

    telemetry_result: dict[str, object] = {
        "stage": stage_name,
        "wall_time_seconds": 0.0,
        "cpu_time_seconds": 0.0,
        "peak_memory_mb": 0.0,
        "invocation_reason": invocation_reason,
        "is_heavy": is_heavy,
    }

    try:
        yield telemetry_result
    finally:
        elapsed_wall = time.perf_counter() - start_wall
        elapsed_cpu = time.process_time() - start_cpu
        _current_bytes, peak_bytes = tracemalloc.get_traced_memory()
        peak_mb = peak_bytes / (1024.0 * 1024.0)

        if tracemalloc_started_here:
            tracemalloc.stop()

        telemetry_result["wall_time_seconds"] = round(elapsed_wall, 6)
        telemetry_result["cpu_time_seconds"] = round(elapsed_cpu, 6)
        telemetry_result["peak_memory_mb"] = round(peak_mb, 3)

        account.observe_peak("peak_memory_mb", peak_mb)
        account.add_usage("cpu_time_seconds", round(elapsed_cpu, 6))

        if is_heavy:
            account.record_heavy_stage(
                stage_name,
                invocation_reason,
                wall_time_seconds=elapsed_wall,
                cpu_seconds=elapsed_cpu,
                peak_memory_mb=peak_mb,
                metadata=metadata,
            )
