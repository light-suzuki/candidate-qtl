"""Post-hoc provider-reported LLM usage telemetry for issue #434.

``llm-run`` already preserves each raw provider payload.  This module aggregates
only usage values explicitly returned by the provider; it does not estimate
missing token counts with a tokenizer and does not attach a price table.

Supported common shapes:

* OpenAI Responses: ``usage.input_tokens`` / ``usage.output_tokens``;
* OpenAI-compatible chat completions: ``usage.prompt_tokens`` /
  ``usage.completion_tokens``.

Providers such as LM Studio CLI that do not return structured usage remain
``not_available`` rather than receiving invented estimates.
"""

from __future__ import annotations

import json

import pandas as pd

from .core import runtime_path

SCHEMA_VERSION = "1.0"
USAGE_COLUMNS = [
    "prompt_id",
    "gene_id",
    "provider",
    "input_tokens",
    "output_tokens",
    "total_tokens",
    "usage_status",
    "usage_shape",
]


def _int_or_none(value: object) -> int | None:
    if value is None or value == "":
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _usage_from_payload(payload: dict[str, object]) -> tuple[int | None, int | None, int | None, str]:
    usage = payload.get("usage")
    if not isinstance(usage, dict):
        return None, None, None, "none"

    input_tokens = _int_or_none(usage.get("input_tokens"))
    output_tokens = _int_or_none(usage.get("output_tokens"))
    shape = "responses" if input_tokens is not None or output_tokens is not None else ""
    if input_tokens is None and output_tokens is None:
        input_tokens = _int_or_none(usage.get("prompt_tokens"))
        output_tokens = _int_or_none(usage.get("completion_tokens"))
        shape = "chat_completions" if input_tokens is not None or output_tokens is not None else "unknown"
    total_tokens = _int_or_none(usage.get("total_tokens"))
    if total_tokens is None and input_tokens is not None and output_tokens is not None:
        total_tokens = input_tokens + output_tokens
    return input_tokens, output_tokens, total_tokens, shape


def summarize_llm_usage(
    raw_responses_arg: str,
    out_dir_arg: str,
    *,
    llm_manifest_arg: str | None = None,
    case_id: str = "",
) -> dict[str, object]:
    """Aggregate provider-reported usage from ``llm_raw_responses.jsonl``."""
    raw_path = runtime_path(raw_responses_arg)
    if not raw_path.exists():
        raise ValueError(f"LLM raw responses do not exist: {raw_path}")

    provider = ""
    model = ""
    expected_prompt_count: int | None = None
    if llm_manifest_arg:
        manifest_path = runtime_path(llm_manifest_arg)
        if not manifest_path.exists():
            raise ValueError(f"LLM manifest does not exist: {manifest_path}")
        llm_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if not isinstance(llm_manifest, dict):
            raise ValueError("LLM manifest must be a JSON object")
        provider = str(llm_manifest.get("provider", ""))
        model = str(llm_manifest.get("model", ""))
        expected_prompt_count = _int_or_none(llm_manifest.get("prompt_count"))

    rows: list[dict[str, object]] = []
    seen_prompts: set[str] = set()
    for line_number, line in enumerate(raw_path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        record = json.loads(line)
        if not isinstance(record, dict) or not isinstance(record.get("response"), dict):
            raise ValueError(f"invalid LLM raw response at line {line_number}")
        prompt_id = str(record.get("prompt_id", "")).strip()
        if not prompt_id:
            raise ValueError(f"LLM raw response is missing prompt_id at line {line_number}")
        if prompt_id in seen_prompts:
            raise ValueError(f"duplicate LLM raw response prompt_id: {prompt_id}")
        seen_prompts.add(prompt_id)
        response = record["response"]
        payload_provider = str(response.get("provider", "")).strip()
        row_provider = provider or payload_provider
        input_tokens, output_tokens, total_tokens, usage_shape = _usage_from_payload(response)
        values = (input_tokens, output_tokens, total_tokens)
        if input_tokens is not None and output_tokens is not None:
            usage_status = "complete"
        elif any(value is not None for value in values):
            usage_status = "partial"
        else:
            usage_status = "not_available"
        rows.append({
            "prompt_id": prompt_id,
            "gene_id": str(record.get("gene_id", "")),
            "provider": row_provider,
            "input_tokens": "" if input_tokens is None else input_tokens,
            "output_tokens": "" if output_tokens is None else output_tokens,
            "total_tokens": "" if total_tokens is None else total_tokens,
            "usage_status": usage_status,
            "usage_shape": usage_shape,
        })

    frame = pd.DataFrame(rows, columns=USAGE_COLUMNS).sort_values("prompt_id", kind="stable") if rows \
        else pd.DataFrame(columns=USAGE_COLUMNS)
    call_count = int(len(frame))
    if expected_prompt_count is not None and call_count != expected_prompt_count:
        raise ValueError(
            f"LLM usage roster mismatch: raw call count {call_count} != manifest prompt_count {expected_prompt_count}"
        )

    def numeric_sum(column: str) -> int:
        if frame.empty:
            return 0
        return int(pd.to_numeric(frame[column], errors="coerce").fillna(0).sum())

    complete_count = int(frame["usage_status"].eq("complete").sum()) if len(frame) else 0
    partial_count = int(frame["usage_status"].eq("partial").sum()) if len(frame) else 0
    unavailable_count = int(frame["usage_status"].eq("not_available").sum()) if len(frame) else 0
    if call_count == 0:
        aggregate_status = "no_calls"
    elif complete_count == call_count:
        aggregate_status = "complete"
    elif unavailable_count == call_count:
        aggregate_status = "not_available"
    else:
        aggregate_status = "partial"

    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    usage_path = out / "llm_usage.tsv"
    frame.to_csv(usage_path, sep="\t", index=False, lineterminator="\n")
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id,
        "provider": provider,
        "model": model,
        "call_count": call_count,
        "input_tokens": numeric_sum("input_tokens"),
        "output_tokens": numeric_sum("output_tokens"),
        "total_tokens": numeric_sum("total_tokens"),
        "usage_complete_call_count": complete_count,
        "usage_partial_call_count": partial_count,
        "usage_unavailable_call_count": unavailable_count,
        "usage_status": aggregate_status,
        "token_counts_are_provider_reported": True,
        "token_estimation_performed": False,
        "price_estimation_performed": False,
        "outputs": {
            "usage": str(usage_path),
            "manifest": str(out / "llm_usage_manifest.json"),
        },
        "policy": [
            "token counts are aggregated only when explicitly returned by the provider",
            "missing token usage is not estimated with a tokenizer",
            "no API price table is applied here; cost should be computed separately from a versioned price source if needed",
            "call_count is auditable from the raw response roster",
        ],
    }
    manifest_path = out / "llm_usage_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "usage": str(usage_path),
        "manifest": str(manifest_path),
        "call_count": call_count,
        "input_tokens": manifest["input_tokens"],
        "output_tokens": manifest["output_tokens"],
        "total_tokens": manifest["total_tokens"],
        "usage_status": aggregate_status,
    }


def apply_llm_usage_to_account(
    account: object,
    summary_or_manifest: dict[str, object],
) -> dict[str, object]:
    """Charge provider-reported usage onto a BudgetAccount ledger."""
    call_count = int(summary_or_manifest.get("call_count", 0))
    input_tokens = int(summary_or_manifest.get("input_tokens", 0))
    output_tokens = int(summary_or_manifest.get("output_tokens", 0))
    if call_count > 0 and hasattr(account, "consume"):
        account.consume("llm_call_budget", call_count)
    if input_tokens > 0 and hasattr(account, "consume"):
        account.consume("llm_input_tokens", input_tokens)
    if output_tokens > 0 and hasattr(account, "consume"):
        account.consume("llm_output_tokens", output_tokens)
    return {
        "call_count": call_count,
        "input_tokens": input_tokens,
        "output_tokens": output_tokens,
        "total_tokens": input_tokens + output_tokens,
    }
