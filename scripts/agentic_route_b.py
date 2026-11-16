"""Prepare bounded Route B prompt shards and merge partial agent JSON replies.

This is orchestration glue: it preserves frozen prompt IDs and unit IDs and
never decides which biological units or claims matter.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

if __package__:
    from .agentic_route_b_material import CONTRACT as MATERIAL_CONTRACT, DEFAULT_MAX_MATERIAL_BYTES, plan_reading_views, restore_reading_shard
else:
    from agentic_route_b_material import CONTRACT as MATERIAL_CONTRACT, DEFAULT_MAX_MATERIAL_BYTES, plan_reading_views, restore_reading_shard


EXECUTION_CONTRACT = "agentic-route-b-execution-1.1.0"


class _ModelResponseFormatError(ValueError):
    """Untrusted model syntax/schema failure, distinct from artifact integrity."""

    def __init__(self, reason: str, message: str) -> None:
        super().__init__(message)
        self.reason = reason


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _prompts(path: Path, case_id: str = "") -> tuple[list[dict[str, Any]], str]:
    rows = []
    seen: set[str] = set()
    seen_units: set[str] = set()
    case_ids: set[str] = set()
    for n, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        row = json.loads(line)
        if not isinstance(row, dict) or not row.get("prompt_id") or not row.get("route_b_unit_id"):
            raise ValueError(f"invalid prompt record at line {n}: prompt_id and route_b_unit_id required")
        pid = str(row["prompt_id"])
        if pid in seen:
            raise ValueError(f"duplicate prompt_id: {pid}")
        uid = str(row["route_b_unit_id"])
        if uid in seen_units:
            raise ValueError(f"duplicate route_b_unit_id: {uid}")
        seen.add(pid)
        seen_units.add(uid)
        case_ids.add(str(row.get("case_id", "")))
        rows.append(row)
    if rows:
        if len(case_ids) != 1 or not next(iter(case_ids)):
            raise ValueError("prompts must have one non-empty case_id")
        row_case_id = next(iter(case_ids))
        if case_id and row_case_id != case_id:
            raise ValueError("prompt roster case_id mismatch")
        return rows, row_case_id
    if not case_id:
        for manifest_name in ("route_b_semantic_bundle_manifest.json", "prepare_manifest.json"):
            manifest_path = path.parent / manifest_name
            if manifest_path.is_file():
                manifest = _read_json(manifest_path)
                if isinstance(manifest, dict) and str(manifest.get("case_id", "")).strip():
                    case_id = str(manifest["case_id"]).strip()
                    break
    if not case_id:
        raise ValueError("empty prompt roster requires an explicit case_id")
    return rows, case_id


def _research_ledgers(paths: list[str], digests: list[str]) -> list[tuple[dict[str, str], bytes]]:
    """Pin optional research ledgers for a fresh-context transcription worker.

    Each path must match its coordinator-supplied SHA-256 exactly; a run that
    cannot prove which ledger bytes it hands over is rejected before any output.
    The verified bytes are copied under the prepare output root so the worker
    keeps its read root and the task text does not depend on the host path.
    """
    if len(paths) != len(digests):
        raise ValueError("research ledger path and SHA-256 must be supplied together")
    ledgers: list[tuple[dict[str, str], bytes]] = []
    for raw, expected in zip(paths, digests):
        path = Path(raw).resolve()
        if not path.is_file():
            raise ValueError(f"research ledger is not a readable file: {raw}")
        data = path.read_bytes()
        digest = str(expected).strip().lower()
        if not digest or hashlib.sha256(data).hexdigest() != digest:
            raise ValueError(f"research ledger SHA-256 mismatch: {path}")
        name = f"research-ledger/{len(ledgers):02d}-{digest[:16]}{path.suffix}"
        ledgers.append(({"file": name, "sha256": digest, "source_path": str(path)}, data))
    return ledgers


def prepare(prompts_path: str, out_dir: str, max_prompts: int = 8, *, case_id: str = "",
            max_material_bytes: int = DEFAULT_MAX_MATERIAL_BYTES,
            prefer_shared_material: bool = False,
            research_bundle: str = "", research_bundle_sha256: str = "",
            followup_plan: str = "", followup_plan_sha256: str = "",
            research_ledger: list[str] | None = None,
            research_ledger_sha256: list[str] | None = None,
            research_batch_units: int = 0) -> dict[str, Any]:
    if max_prompts < 1:
        raise ValueError("max_prompts must be positive")
    research_ledgers = _research_ledgers(list(research_ledger or []), list(research_ledger_sha256 or []))
    src, out = Path(prompts_path), Path(out_dir)
    if out.exists():
        raise ValueError(f"output already exists: {out}")
    prompts, case_id = _prompts(src, case_id)
    source_sha256 = _sha(src)
    research_files, research_view = {}, None
    if bool(research_bundle) != bool(research_bundle_sha256):
        raise ValueError("research bundle path and SHA-256 must be supplied together")
    if bool(followup_plan) != bool(followup_plan_sha256) or (followup_plan and not research_bundle):
        raise ValueError("follow-up requires a pinned plan and research bundle")
    if type(research_batch_units) is not int or research_batch_units < 0:
        raise ValueError("research batch size must be a non-negative integer")
    if research_batch_units and not research_bundle:
        raise ValueError("--research-batch-units requires --research-bundle")
    if research_batch_units and followup_plan:
        raise ValueError("research batching is limited to the initial phase and cannot be combined with a follow-up plan")
    if research_bundle:
        if __package__:
            from .agentic_route_b_research import plan_research_view
        else:
            from agentic_route_b_research import plan_research_view
        research_files, research_view = plan_research_view(
            research_bundle, research_bundle_sha256, prompts, case_id, max_material_bytes,
            followup_plan=followup_plan, followup_plan_sha256=followup_plan_sha256,
            batch_units=research_batch_units,
        )
    parent_bundle = src.parent / "route_b_semantic_bundle.json"
    parent_hash = ""
    if parent_bundle.is_file():
        parent = _read_json(parent_bundle)
        fingerprints = [hashlib.sha256(json.dumps(prompt, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()).hexdigest() for prompt in prompts]
        if parent.get("case_id") != case_id or parent.get("execution_prompt_fingerprints") != fingerprints:
            raise ValueError("research parent bundle does not bind these prompts")
        parent_hash = _sha(parent_bundle)
    view_files, reading_shards = plan_reading_views(
        prompts, case_id, source_sha256, max_prompts, max_material_bytes,
        prefer_shared_material=prefer_shared_material,
        parent_bundle_sha256=parent_hash,
        research_ledgers=tuple(row for row, _ in research_ledgers),
    )
    shards = []
    for start in range(0, len(prompts), max_prompts):
        sid = f"b{start // max_prompts:04d}"
        items = prompts[start:start + max_prompts]
        target = out / sid / "prompts.jsonl"
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("".join(json.dumps(p, ensure_ascii=False, separators=(",", ":")) + "\n" for p in items), encoding="utf-8")
        shards.append({"shard_id": sid, "file": str(target.relative_to(out)), "prompt_ids": [p["prompt_id"] for p in items], "unit_ids": [p["route_b_unit_id"] for p in items], "sha256": _sha(target)})
    for name, data in [*view_files.items(), *((row["file"], data) for row, data in research_ledgers)]:
        target = out / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
    manifest = {"schema_version": 1, "case_id": case_id, "source_prompts": str(src.resolve()), "source_sha256": source_sha256, "prompt_count": len(prompts), "max_prompts": max_prompts, "shards": shards,
                "reading_view": {"contract": MATERIAL_CONTRACT, "max_material_bytes": max_material_bytes, "shards": reading_shards,
                                 "oversized_unit_count": sum(s["reading_state"] == "oversized_full_material" for s in reading_shards),
                                 "shared_material_count": sum(name.startswith("shared-material/") for name in view_files),
                                 "stored_bytes": sum(len(data) for data in view_files.values()), "exact_roundtrip_verified": False}}
    if prefer_shared_material:
        reading_ids = [pid for shard in reading_shards for pid in shard["prompt_ids"]]
        manifest["reading_view"]["packing_strategy"] = (
            "exact_shared_reference_material_affinity_opt_in"
            if reading_ids != [p["prompt_id"] for p in prompts]
            else "canonical_order_no_safe_affinity_gain"
        )
        manifest["reading_view"]["total_reading_bytes"] = sum(
            shard["total_reading_bytes"] for shard in reading_shards
        )
    _write_json(out / "prepare_manifest.json", manifest)
    recovered = [prompt for shard in reading_shards for prompt in restore_reading_shard(out, shard["file"])]
    if prefer_shared_material:
        restored_exactly = (len(recovered) == len(prompts) and
                            {row["prompt_id"]: row for row in recovered} ==
                            {row["prompt_id"]: row for row in prompts})
    else:
        restored_exactly = recovered == prompts
    if not restored_exactly:
        raise ValueError("reading view does not restore the complete frozen prompt roster")
    manifest["reading_view"]["exact_roundtrip_verified"] = True
    if research_view is not None:
        for name, data in research_files.items():
            target = out / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
        manifest["research_view"] = research_view
    manifest["research_source_parent"] = {"path": str(parent_bundle.resolve()), "sha256": parent_hash} if parent_hash else None
    if research_ledgers:
        manifest["research_ledgers"] = [row for row, _ in research_ledgers]
    _write_json(out / "prepare_manifest.json", manifest)
    return manifest


def _parsed_model_response(record: dict[str, Any]) -> dict[str, Any]:
    """Normalize one completed response only when raw and parsed JSON agree."""
    from candidate_qtl.route_b_response_validation import route_b_claim_schema_error

    raw_text = record.get("raw_text")
    model_json = record.get("model_json")
    if raw_text is None:
        raise ValueError(f"response {record.get('prompt_id')} requires raw_text")
    if not isinstance(raw_text, str):
        raise ValueError(f"response {record.get('prompt_id')} raw_text must be a string")
    if isinstance(raw_text, str):
        try:
            raw_json = json.loads(raw_text)
        except json.JSONDecodeError as exc:
            raise _ModelResponseFormatError("model_response_invalid_json", f"response {record.get('prompt_id')} raw_text is not valid JSON") from exc
        if model_json is not None and model_json != raw_json:
            raise ValueError(f"response {record.get('prompt_id')} model_json differs from raw_text")
        model_json = raw_json
    if not isinstance(model_json, dict) or not isinstance(model_json.get("claims"), list):
        raise _ModelResponseFormatError("model_response_invalid_claims_envelope", f"response {record.get('prompt_id')} must contain a JSON object with claims list")
    rejected = [
        {"claim_index": index, "raw_claim": claim, "reason": reason}
        for index, claim in enumerate(model_json["claims"], 1)
        if (reason := route_b_claim_schema_error(claim))
    ]
    normalized = {**record, "model_json": model_json, "raw_text": raw_text, "claims": model_json["claims"], "parse_status": "json"}
    # Completion is computed from the raw response, never from caller-supplied
    # bookkeeping that could bypass the multiple-completed-attempt guard.
    normalized.pop("rejected_claims", None)
    normalized.pop("response_completion_status", None)
    if rejected:
        normalized["rejected_claims"] = rejected
        normalized["response_completion_status"] = "partial"
        normalized["review_status"] = "not_evaluated"
        normalized["not_evaluated_reason"] = (
            "model_response_invalid_claims_members"
            if all(not isinstance(row["raw_claim"], dict) for row in rejected)
            else "model_response_invalid_claim_schema"
        )
        normalized["response_format_error"] = "; ".join(f"claim {row['claim_index']}: {row['reason']}" for row in rejected)
        if len(rejected) == len(model_json["claims"]):
            normalized.update(parse_status="not_evaluated", model_json=None, claims=[])
    return normalized


def _response_completed(record: dict[str, Any]) -> bool:
    from candidate_qtl.route_b_response_validation import route_b_response_completed

    return route_b_response_completed(record)


def preflight(prompts_path: str, response_paths: list[str], *, case_id: str = "", executor: str = "",
              model: str = "", reading_manifest: str = "") -> dict[str, Any]:
    """Check original replies with the merge parser without writing artifacts.

    This checks transport/schema readiness, not frozen source membership or
    biological support. Failed/unprocessed replies remain explicitly incomplete.
    Executor/model given here stand in for reply fields, exactly as in merge.
    """
    prompts, case_id = _prompts(Path(prompts_path), case_id)
    expected = {str(row["prompt_id"]): row for row in prompts}
    required_receipts = _required_read_receipts(reading_manifest, Path(prompts_path), prompts, case_id)
    provenance = {"executor": executor, "model": model}
    errors, responses = [], []
    completed_prompts: set[str] = set()
    for arg in response_paths:
        path = Path(arg)
        def error(location: str, reason: str, **details: Any) -> None:
            errors.append({"file": str(path), "path": location, "reason": reason, **details})
        try:
            payload = _read_json(path)
        except (OSError, ValueError) as exc:
            error("$", str(exc))
            continue
        if not isinstance(payload, dict):
            error("$", "partial response must be a JSON object")
            continue
        if payload.get("case_id") != case_id:
            error("$.case_id", "case_id mismatch")
        if payload.get("discovered_sources"):
            error("$.discovered_sources", "new research sources require re-freeze before claim merge")
        reused = {}
        if "reuse" in payload:
            try:
                reused = _validate_reuse(payload, prompts_path)
            except (OSError, ValueError, KeyError) as exc:
                error("$.reuse", str(exc))
                continue
        for key in ("executor", "model"):
            if not str(payload.get(key) or provenance[key]).strip():
                error(f"$.{key}", "response requires executor/model provenance (pass --executor/--model as merge does)")
        records = payload.get("responses")
        if not isinstance(records, list):
            error("$.responses", "preflight requires original responses list; checkpoints use merge manifest validation")
            continue
        seen = set()
        for index, record in enumerate(records):
            location = f"$.responses[{index}]"
            if not isinstance(record, dict) or not record.get("prompt_id"):
                error(location, "response requires prompt_id", response_index=index + 1)
                continue
            pid = str(record["prompt_id"])
            details = {"response_index": index + 1, "prompt_id": pid}
            if pid not in expected:
                error(location + ".prompt_id", "unknown prompt_id", **details)
                continue
            if pid in seen:
                error(location + ".prompt_id", "duplicate prompt_id response in one attempt", **details)
            seen.add(pid)
            if str(record.get("case_id", case_id)) != case_id:
                error(location + ".case_id", "case_id mismatch for response", **details)
            if "route_b_unit_id" in record and str(record["route_b_unit_id"]) != str(expected[pid]["route_b_unit_id"]):
                error(location + ".route_b_unit_id", "unit_id mismatch for response", **details)
            state = str(record.get("status", "not_evaluated" if record.get("parse_status") == "not_evaluated" else "completed"))
            if state in {"not_evaluated", "failed"}:
                reason = str(record.get("reason") or record.get("not_evaluated_reason", ""))
                if not reason:
                    error(location + ".reason", "not_evaluated response requires reason", **details)
                responses.append({"file": str(path), **details, "completed": False,
                                  "not_evaluated_reason": reason})
                continue
            if state != "completed":
                error(location + ".status", f"invalid response status: {state}", **details)
                continue
            try:
                parsed = _parsed_model_response(record)
            except (ValueError, TypeError) as exc:
                error(location, str(exc), **details)
                continue
            for rejected in parsed.get("rejected_claims", []):
                claim_index = rejected["claim_index"]
                error(location + f".raw_text.claims[{claim_index - 1}]", rejected["reason"],
                      **details, claim_index=claim_index)
            if reused:
                _apply_reused_receipt(parsed, reused[pid])
            elif required_receipts:
                _apply_read_receipt(parsed, record, required_receipts.get(pid, set()))
                details["read_receipt_state"] = parsed["read_receipt_state"]
            completed = _response_completed(parsed)
            if completed:
                if pid in completed_prompts:
                    error(location + ".prompt_id", f"multiple completed responses for prompt_id: {pid}", **details)
                completed_prompts.add(pid)
            responses.append({"file": str(path), **details, "completed": completed})
    return {"contract": EXECUTION_CONTRACT, "case_id": case_id,
            "status": "rejected" if errors else "accepted", "response_count": len(responses),
            "validation_scope": "original_response_transport_and_claim_schema; no_support_verification",
            "responses": responses, "errors": errors}


def _required_read_receipts(reading_manifest: str, prompts_file: Path,
                            prompts: list[dict[str, Any]], case_id: str) -> dict[str, set[str]]:
    """Per-prompt read receipts sealed by ``prepare``; empty for old generations."""
    if not reading_manifest:
        return {}
    path = Path(reading_manifest)
    manifest = _read_json(path)
    if manifest.get("case_id") != case_id or manifest.get("source_sha256") != _sha(prompts_file):
        raise ValueError("reading manifest case/prompt hash mismatch")
    required: dict[str, set[str]] = {}
    recovered = []
    for shard in (manifest.get("reading_view") or {}).get("shards", []):
        recovered.extend(restore_reading_shard(path.parent, shard["file"]))
        receipts_by_prompt = shard.get("read_receipts_by_prompt")
        if receipts_by_prompt is None:
            raise ValueError("reading manifest has no read receipts; it predates read-receipt prepare")
        if set(receipts_by_prompt) != set(shard["prompt_ids"]):
            raise ValueError("reading receipt prompt roster mismatch")
        assignment = _read_json(path.parent / shard["file"])
        view = shard["display_view"]
        descriptors = {item["source_file"]: item for item in view["materials"]}
        material_paths = {item["material_id"]: item["file"] for item in assignment["materials"]}
        for pid, receipts in receipts_by_prompt.items():
            if pid in required or not isinstance(receipts, list) or not receipts:
                raise ValueError("invalid or duplicate reading receipt roster")
            # Reconstruct the sealed mapping from the canonical assignment;
            # deleting a required token from manifest metadata must not pass.
            prompt_view = next(item for item in assignment["prompts"] if item["prompt"]["prompt_id"] == pid)
            references = {prompt_view[key] for key in ("system_material_id", "user_prefix_material_id") if key in prompt_view}
            references.update(row["material_id"] for rows in prompt_view.get("unit_tables", {}).values() for row in rows)
            if __package__:
                from .agentic_route_b_material import read_receipt_token
            else:
                from agentic_route_b_material import read_receipt_token
            expected = {read_receipt_token(view["assignment"]["source_sha256"])}
            expected.update(read_receipt_token(descriptors[material_paths[mid]]["source_sha256"]) for mid in references)
            if set(receipts) != expected:
                raise ValueError("reading receipt material roster mismatch")
            required[pid] = expected
    if not required and prompts:
        raise ValueError("reading manifest has no read receipts; it predates read-receipt prepare")
    if {p["prompt_id"]: p for p in recovered} != {p["prompt_id"]: p for p in prompts} or len(recovered) != len(prompts):
        raise ValueError("reading manifest does not restore the full prompt roster")
    return required


def _validate_reuse(payload: dict[str, Any], prompts_path: str) -> dict[str, dict[str, Any]]:
    if __package__:
        from .agentic_route_b_extend import validate_reuse
    else:
        from agentic_route_b_extend import validate_reuse
    return validate_reuse(payload, prompts_path)


def _apply_reused_receipt(candidate: dict[str, Any], prior: dict[str, Any]) -> None:
    candidate["read_receipt_state"] = "reused_prior_assessment"
    candidate["prior_read_receipt_state"] = prior.get("read_receipt_state", "legacy_not_recorded")


def _apply_read_receipt(candidate: dict[str, Any], record: dict[str, Any], required: set[str] | None) -> None:
    """Mark a completed response's receipt coverage; demote unread no-claim answers.

    A receipt is transport evidence that the worker opened each linked file to
    its end, not proof of understanding. A no-claim answer without it cannot be
    distinguished from an unread unit, so it stays not_evaluated (never negative).
    Claims are kept for review but flagged.
    """
    if required is None:
        return
    supplied = record.get("read_receipts")
    supplied = {str(item) for item in supplied} if isinstance(supplied, list) else set()
    missing = sorted(required - supplied)
    candidate["read_receipt_state"] = "complete" if not missing else "incomplete"
    if missing:
        candidate["read_receipt_missing_count"] = len(missing)
        if not candidate.get("claims"):
            candidate.update(parse_status="not_evaluated", model_json=None, claims=[],
                             review_status="not_evaluated",
                             not_evaluated_reason="read_receipt_incomplete")


def merge(prompts_path: str, response_paths: list[str], out_path: str, *, executor: str = "", model: str = "",
          case_id: str = "", reading_manifest: str = "") -> dict[str, Any]:
    prompts_file, out = Path(prompts_path), Path(out_path)
    manifest_path = out.with_name(out.stem + "_manifest.json")
    output_tmp = out.with_name(out.name + ".tmp")
    manifest_tmp = manifest_path.with_name(manifest_path.name + ".tmp")
    if out.exists() or manifest_path.exists():
        raise ValueError(f"output pair already exists: {out}")
    if output_tmp.exists() or manifest_tmp.exists():
        raise ValueError(f"temporary output already exists: {output_tmp}")
    prompts, case_id = _prompts(prompts_file, case_id)
    required_receipts = _required_read_receipts(reading_manifest, prompts_file, prompts, case_id)
    expected = {str(p["prompt_id"]): p for p in prompts}
    gathered: dict[str, dict[str, Any]] = {}
    attempts_by_pid: dict[str, list[dict[str, Any]]] = {}
    provenance = []
    provenance_keys: set[tuple[str, str]] = set()
    retry_history: list[dict[str, Any]] = []
    ordered_inputs = []
    for arg in response_paths:
        candidate_path = Path(arg)
        if not candidate_path.is_file():
            raise ValueError(f"response path does not exist: {candidate_path}")
        ordered_inputs.append((_sha(candidate_path), str(candidate_path.resolve()), arg))
    for attempt_index, (_, _, arg) in enumerate(sorted(ordered_inputs)):
        path = Path(arg)
        payload = _read_json(path)
        if not isinstance(payload, dict):
            raise ValueError(f"partial response must be a JSON object: {path}")
        if not str(payload.get("case_id", "")).strip():
            raise ValueError(f"partial response requires case_id: {path}")
        if str(payload["case_id"]) != case_id:
            raise ValueError(f"case_id mismatch in {path}")
        if payload.get("discovered_sources"):
            raise ValueError("new research sources require re-freeze before claim merge; use route-b-semantic-bundle --extend-bundle with this discovery artifact")
        reused = _validate_reuse(payload, prompts_path) if "reuse" in payload else {}
        is_aggregate = "candidates" in payload and path.with_name(path.stem + "_manifest.json").is_file()
        agent = str(payload.get("executor") or executor).strip()
        agent_model = str(payload.get("model") or model).strip()
        if not is_aggregate and (not agent or not agent_model):
            raise ValueError(f"response requires executor/model provenance: {path}")
        records = payload.get("responses", payload.get("candidates"))
        if not is_aggregate and "responses" not in payload and "candidates" in payload:
            raise ValueError(f"aggregate response requires sibling execution manifest: {path.with_name(path.stem + '_manifest.json')}")
        if not isinstance(records, list):
            raise ValueError(f"partial response requires responses list: {path}")
        if is_aggregate:
            checkpoint_manifest_path = path.with_name(path.stem + "_manifest.json")
            if not checkpoint_manifest_path.is_file():
                raise ValueError(f"aggregate response requires sibling execution manifest: {checkpoint_manifest_path}")
            checkpoint = _read_json(checkpoint_manifest_path)
            if not isinstance(checkpoint, dict) or checkpoint.get("contract") != EXECUTION_CONTRACT:
                raise ValueError("aggregate response has an unsupported execution manifest")
            if checkpoint.get("case_id") != case_id:
                raise ValueError("aggregate response execution manifest case_id mismatch")
            if checkpoint.get("prompt_sha256") != _sha(prompts_file) or checkpoint.get("response_sha256") != _sha(path):
                raise ValueError("aggregate response execution manifest input hash mismatch")
            prior_provenance = checkpoint.get("partial_response_provenance")
            if not isinstance(prior_provenance, list) or not prior_provenance:
                raise ValueError("aggregate response execution manifest lacks attempt provenance")
            for prior in prior_provenance:
                if not isinstance(prior, dict) or type(prior.get("attempt_index")) is not int or prior["attempt_index"] < 0:
                    raise ValueError("aggregate response attempt provenance is invalid")
                key = (str(prior.get("file", "")), str(prior.get("sha256", "")))
                if key in provenance_keys:
                    continue
                provenance.append(prior)
                provenance_keys.add(key)
            prior_history = checkpoint.get("retry_history")
            if not isinstance(prior_history, list) or any(not isinstance(row, dict) for row in prior_history):
                raise ValueError("aggregate response retry history is invalid")
            retry_history.extend(dict(row) for row in prior_history)
            for prior in prior_provenance:
                if "reuse" not in prior.get("metadata", {}):
                    continue
                prior_file = Path(prior["file"])
                if _sha(prior_file) != prior["sha256"]:
                    raise ValueError("reused raw response hash mismatch")
                prior_payload = _read_json(prior_file)
                validated = _validate_reuse(prior_payload, prompts_path)
                original_rows = {r["prompt_id"]: r for r in prior_payload["responses"]}
                for record in records:
                    pid = record.get("prompt_id")
                    if record.get("attempt_source_file") != str(prior_file.resolve()) or pid not in validated:
                        continue
                    if (record.get("agent_response_raw") != original_rows[pid]
                            or record.get("raw_text") != original_rows[pid].get("raw_text")
                            or record.get("attempt_executor") != prior_payload["executor"]
                            or record.get("attempt_model") != prior_payload["model"]):
                        raise ValueError("checkpoint differs from reused original response")
                    reused[pid] = validated[pid]
        else:
            provenance_row = {"attempt_index": attempt_index, "file": str(path.resolve()), "sha256": _sha(path), "executor": agent, "model": agent_model, "response_id": payload.get("response_id", ""), "metadata": {k: v for k, v in payload.items() if k not in {"responses", "candidates"}}}
            provenance_key = (provenance_row["file"], provenance_row["sha256"])
            if provenance_key not in provenance_keys:
                provenance.append(provenance_row)
                provenance_keys.add(provenance_key)
        seen_in_file: set[str] = set()
        for record in records:
            if not isinstance(record, dict) or not record.get("prompt_id"):
                raise ValueError(f"response in {path} requires prompt_id")
            pid = str(record["prompt_id"])
            if pid not in expected:
                raise ValueError(f"unknown prompt_id: {pid}")
            if pid in seen_in_file:
                raise ValueError(f"duplicate prompt_id response in one attempt: {pid}")
            seen_in_file.add(pid)
            prompt_case = str(record.get("case_id", case_id))
            if prompt_case != case_id:
                raise ValueError(f"case_id mismatch for response {pid}")
            if "route_b_unit_id" in record and str(record["route_b_unit_id"]) != str(expected[pid]["route_b_unit_id"]):
                raise ValueError(f"unit_id mismatch for response {pid}")
            if (
                is_aggregate
                and record.get("parse_status") == "not_evaluated"
                and not str(record.get("raw_text", "")).strip()
                and not record.get("agent_response_raw")
                and record.get("not_evaluated_reason") == "no_partial_response"
            ):
                continue
            state = str(record.get("status", "not_evaluated" if record.get("parse_status") == "not_evaluated" else "completed"))
            if is_aggregate:
                attempt_source_file = str(record.get("attempt_source_file") or path.resolve())
                attempt_executor = str(record.get("attempt_executor") or agent).strip()
                attempt_model = str(record.get("attempt_model") or agent_model).strip()
            else:
                attempt_source_file = str(path.resolve())
                attempt_executor = agent
                attempt_model = agent_model
            if not attempt_executor or not attempt_model:
                raise ValueError(f"response attempt requires executor/model provenance: {path}")
            original_response = record.get("agent_response_raw", record) if is_aggregate else record
            format_error = None
            if state == "completed":
                try:
                    normalized = _parsed_model_response(record)
                except _ModelResponseFormatError as exc:
                    format_error = exc
                    state = "not_evaluated"
            if state in {"not_evaluated", "failed"}:
                reason = format_error.reason if format_error else str(record.get("reason") or record.get("not_evaluated_reason", ""))
                if not reason:
                    raise ValueError(f"not_evaluated response {pid} requires reason")
                candidate = {"prompt_id": pid, "route_b_unit_id": expected[pid]["route_b_unit_id"], "parse_status": "not_evaluated", "model_json": None, "claims": [], "raw_text": record.get("raw_text", ""), "review_status": "not_evaluated", "not_evaluated_reason": reason, "agent_response_raw": original_response, "attempt_source_file": attempt_source_file, "attempt_executor": attempt_executor, "attempt_model": attempt_model}
                if format_error or record.get("response_format_error"):
                    candidate["response_format_error"] = str(format_error) if format_error else record["response_format_error"]
                if record.get("rejected_claims"):
                    candidate["rejected_claims"] = record["rejected_claims"]
                    candidate["response_completion_status"] = "partial"
            elif state == "completed":
                normalized.setdefault("route_b_unit_id", expected[pid]["route_b_unit_id"])
                if str(normalized["route_b_unit_id"]) != str(expected[pid]["route_b_unit_id"]):
                    raise ValueError(f"unit_id mismatch for response {pid}")
                normalized["agent_response_raw"] = original_response
                normalized["attempt_source_file"] = attempt_source_file
                normalized["attempt_executor"] = attempt_executor
                normalized["attempt_model"] = attempt_model
                candidate = normalized
                if pid in reused:
                    _apply_reused_receipt(candidate, reused[pid])
                elif required_receipts:
                    _apply_read_receipt(candidate, record, required_receipts.get(pid, set()))
            else:
                raise ValueError(f"invalid response status for {pid}: {state}")
            attempts_by_pid.setdefault(pid, []).append({
                "candidate": candidate,
                "source_file": attempt_source_file,
                "executor": attempt_executor,
                "model": attempt_model,
            })
            previous = gathered.get(pid)
            if previous is not None:
                previous_state = "completed" if _response_completed(previous) else "not_evaluated"
                current_state = "completed" if _response_completed(candidate) else "not_evaluated"
                if previous_state == "completed":
                    if current_state == "completed":
                        raise ValueError(f"multiple completed responses for prompt_id: {pid}")
                    # A late failure must not replace a completed response. Inputs
                    # are sorted by content hash, so this is order-independent.
                    continue
                if current_state != "completed":
                    # Prefer an incomplete response with usable claims over a
                    # wholly failed attempt, independently of input order.
                    if previous.get("parse_status") == "json" or candidate.get("parse_status") != "json":
                        continue
            gathered[pid] = candidate
    for pid, attempts in attempts_by_pid.items():
        final = gathered.get(pid)
        if final is None or not _response_completed(final):
            continue
        failed_attempt = next((attempt for attempt in attempts if not _response_completed(attempt["candidate"])), None)
        if failed_attempt is None:
            continue
        replaced_file = Path(str(failed_attempt["source_file"]))
        replacement_file = Path(str(final["attempt_source_file"]))
        retry_history.append({
            "prompt_id": pid,
            "route_b_unit_id": expected[pid]["route_b_unit_id"],
            "replaced_state": "not_evaluated",
            "replaced_reason": failed_attempt["candidate"].get("not_evaluated_reason", ""),
            "replaced_response_file": str(replaced_file),
            "replaced_response_sha256": _sha(replaced_file) if replaced_file.is_file() else "",
            "replacement_response_file": str(replacement_file),
            "replacement_response_sha256": _sha(replacement_file) if replacement_file.is_file() else "",
            "replacement_executor": final["attempt_executor"],
            "replacement_model": final["attempt_model"],
            **({
                "replaced_raw_text": failed_attempt["candidate"]["raw_text"],
                "replaced_raw_text_sha256": hashlib.sha256(failed_attempt["candidate"]["raw_text"].encode("utf-8")).hexdigest(),
                "replaced_agent_response_raw": failed_attempt["candidate"]["agent_response_raw"],
                "replaced_response_format_error": failed_attempt["candidate"]["response_format_error"],
            } if failed_attempt["candidate"].get("response_format_error") else {}),
        })
    candidates = []
    for prompt in prompts:
        pid = str(prompt["prompt_id"])
        if pid in gathered:
            row = gathered[pid]
        else:
            row = {"prompt_id": pid, "route_b_unit_id": prompt["route_b_unit_id"], "parse_status": "not_evaluated", "model_json": None, "claims": [], "raw_text": "", "review_status": "not_evaluated", "not_evaluated_reason": "no_partial_response"}
        row.setdefault("gene_id", "")
        row.setdefault("review_status", "not_reviewed")
        candidates.append(row)
    retry_history.sort(key=lambda row: (
        str(row.get("prompt_id", "")), str(row.get("replaced_response_file", "")),
        str(row.get("replacement_response_file", "")),
    ))
    provenance.sort(key=lambda row: (str(row.get("file", "")), str(row.get("sha256", ""))))
    for index, row in enumerate(provenance):
        row["attempt_index"] = index
    payload = {"case_id": case_id, "candidates": candidates}
    out.parent.mkdir(parents=True, exist_ok=True)
    _write_json(output_tmp, payload)
    response_digest = _sha(output_tmp)
    manifest = {
        "contract": EXECUTION_CONTRACT,
        "schema_version": 1,
        "case_id": case_id,
        "source_prompts": str(prompts_file.resolve()),
        "output": str(out.resolve()),
        "completion_status": "not_evaluated" if not prompts else ("complete" if len(gathered) == len(expected) and all(_response_completed(r) for r in gathered.values()) else "partial"),
        "prompt_count": len(prompts),
        "completed_count": sum(_response_completed(r) for r in gathered.values()),
        "not_evaluated_count": len(prompts) - sum(_response_completed(r) for r in gathered.values()),
        "prompt_sha256": _sha(prompts_file),
        "response_sha256": response_digest,
        "partial_response_provenance": provenance,
        "retry_history": retry_history,
        "response_roster": [
            {
                "prompt_id": r["prompt_id"],
                "route_b_unit_id": r["route_b_unit_id"],
                "parse_status": r["parse_status"],
                "completion_status": "completed" if _response_completed(r) else "not_evaluated",
                "not_evaluated_reason": r.get("not_evaluated_reason", ""),
                "raw_text_sha256": hashlib.sha256(str(r.get("raw_text", "")).encode()).hexdigest(),
            }
            for r in candidates
        ],
    }
    from candidate_qtl.route_b_response_validation import route_b_execution_breakdown

    manifest["response_state_counts"] = route_b_execution_breakdown(candidates)
    if required_receipts:
        states = [str(r.get("read_receipt_state", "")) for r in candidates]
        manifest["read_receipts"] = {
            "reading_manifest": str(Path(reading_manifest).resolve()),
            "reading_manifest_sha256": _sha(Path(reading_manifest)),
            "complete": states.count("complete"),
            "reused_prior_assessment": states.count("reused_prior_assessment"),
            "incomplete_with_claims": sum(r.get("read_receipt_state") == "incomplete" and bool(r.get("claims"))
                                          for r in candidates),
            "demoted_to_not_evaluated": sum(r.get("not_evaluated_reason") == "read_receipt_incomplete"
                                            for r in candidates),
            "scope": "transport receipt that each linked display file was opened to its end; not proof of understanding",
        }
    _write_json(manifest_tmp, manifest)
    try:
        output_tmp.replace(out)
        manifest_tmp.replace(manifest_path)
    except Exception:
        if output_tmp.exists():
            output_tmp.unlink()
        if manifest_tmp.exists():
            manifest_tmp.unlink()
        raise
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare"); p.add_argument("--prompts", required=True); p.add_argument("--out", required=True); p.add_argument("--max-prompts", type=int, default=8); p.add_argument("--max-material-bytes", type=int, default=DEFAULT_MAX_MATERIAL_BYTES, help="reading assignment plus linked material byte budget; oversized single units retain full material separately"); p.add_argument("--case-id", default=""); p.add_argument("--prefer-shared-material", action="store_true", help="opt-in source-affinity reading assignment; canonical prompts and decisions are unchanged")
    p.add_argument("--research-bundle", default="", help="opt-in ordinary research navigation over ALL frozen B units; not a model run or claim response")
    p.add_argument("--research-bundle-sha256", default="", help="parent-pinned SHA-256 of the full research bundle; required with --research-bundle")
    p.add_argument("--followup-plan", default="", help="consume a shared-ledger research plan; deliver only unread follow-up context")
    p.add_argument("--followup-plan-sha256", default="", help="coordinator-pinned follow-up plan SHA-256")
    p.add_argument("--research-ledger", action="append", default=[], help="opt-in original research ledger to transcribe from a fresh context; repeat in the same order as --research-ledger-sha256")
    p.add_argument("--research-ledger-sha256", action="append", default=[], help="coordinator-pinned SHA-256 of each --research-ledger; mismatch is rejected")
    p.add_argument("--research-batch-units", type=int, default=0, help="opt-in fresh-context research batches of at most N requested units; requires --research-bundle and no follow-up plan")
    p = sub.add_parser("merge"); p.add_argument("--prompts", required=True); p.add_argument("--response", action="append", default=[]); p.add_argument("--out", required=True); p.add_argument("--executor", default=""); p.add_argument("--model", default=""); p.add_argument("--case-id", default="")
    p.add_argument("--reading-manifest", default="", help="prepare_manifest.json; verify per-unit read receipts and keep unread no-claim answers not_evaluated")
    p = sub.add_parser("preflight", help="read-only transport/schema check of original replies before merge")
    p.add_argument("--prompts", required=True); p.add_argument("--response", action="append", required=True); p.add_argument("--case-id", default="")
    p.add_argument("--executor", default="", help="executor provenance merge will record, when replies omit it")
    p.add_argument("--model", default="", help="model provenance merge will record, when replies omit it")
    p.add_argument("--reading-manifest", default="", help="also report replies whose read receipts are incomplete")
    p = sub.add_parser("reuse-plan", help="compare full frozen prompts and audit old raw responses without merging or executing models")
    p.add_argument("--old-prompts", required=True); p.add_argument("--new-prompts", required=True); p.add_argument("--old-response", required=True)
    p.add_argument("--import-manifest", default=""); p.add_argument("--out", required=True); p.add_argument("--case-id", default="")
    p = sub.add_parser("extend", help="combine discoveries, re-freeze, prepare and export hash-eligible reusable replies; no model execution")
    p.add_argument("--bundle", required=True); p.add_argument("--bundle-sha256", required=True)
    p.add_argument("--old-prompts", required=True); p.add_argument("--old-response", required=True)
    p.add_argument("--discovery", action="append", required=True); p.add_argument("--import-manifest", default="")
    p.add_argument("--out", required=True); p.add_argument("--max-material-bytes", type=int, default=DEFAULT_MAX_MATERIAL_BYTES)
    p = sub.add_parser("summarize-initial", help="account for original initial judgments without changing claims or Tier")
    p.add_argument("--bundle", required=True); p.add_argument("--bundle-sha256", required=True)
    p.add_argument("--assessments", required=True); p.add_argument("--out", required=True)
    p = sub.add_parser("merge-initial", help="concatenate original fresh-context batch ledgers in batch order; no new judgment")
    p.add_argument("--bundle", required=True); p.add_argument("--bundle-sha256", required=True)
    p.add_argument("--batches", required=True, help="research/batches/index.json from prepare --research-batch-units")
    p.add_argument("--batches-sha256", required=True, help="coordinator-pinned SHA-256 of the batch index")
    p.add_argument("--out", required=True, help="new directory for the merged ledger and receipt")
    p = sub.add_parser("plan-light-followup", help="opt-in shared reference screen and provisional T4/T5 research plan; never a Tier/Wet selector")
    p.add_argument("--summary", required=True, help="sealed summarize-initial/initial_assessment_summary.json")
    p.add_argument("--screening", default="", help="optional legacy screen file; default reuses research_screening in the original initial ledger")
    p.add_argument("--out", required=True)
    p = sub.add_parser("light-prepare", help="opt-in light Route B batches: frozen material + roster + task per bounded unit batch; no model call")
    p.add_argument("--bundle", required=True); p.add_argument("--bundle-sha256", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--batch-units", type=int, default=20, help="maximum requested units per light batch")
    p.add_argument("--output-contract", choices=["legacy-v1", "science-v2"], default="legacy-v1")
    p.add_argument("--packing", choices=["original", "reference-affinity"], default="original")
    p.add_argument("--units-file", default="", help="optional newline list of requested unit IDs for a small trial")
    p = sub.add_parser("light-convert", help="check a light Route B run and convert it into discoveries plus an llm-run agent file; no judgment rewrite")
    p.add_argument("--bundle", required=True); p.add_argument("--bundle-sha256", required=True)
    p.add_argument("--light-index", required=True); p.add_argument("--light-index-sha256", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--cache", default="", help="optional JSON metadata cache keyed by paper id")
    p.add_argument("--offline", action="store_true", help="use only the supplied cache; a missing entry is not_attempted")
    for command in ("light-preflight", "light-apply", "light-read-material"):
        p = sub.add_parser(command, help="local hash-bound Light research helper; no model or network call")
        p.add_argument("--bundle", required=True); p.add_argument("--bundle-sha256", required=True)
        p.add_argument("--light-index", required=True); p.add_argument("--light-index-sha256", required=True)
        p.add_argument("--batch-id", required=True)
        if command == "light-read-material":
            p.add_argument("--unit-id", action="append", default=[])
            p.add_argument("--gene-id", action="append", default=[])
            p.add_argument("--field", action="append", default=[])
            p.add_argument("--max-chars", type=int, default=12000)
        else:
            p.add_argument("--results", required=True)
            if command == "light-preflight":
                p.add_argument("--out", required=True)
            else:
                p.add_argument("--patch", required=True); p.add_argument("--patch-sha256", required=True)
                p.add_argument("--out-results", required=True); p.add_argument("--receipt", required=True)
    p = sub.add_parser("light-read-source", help="read selected saved-source ranges with exact input SHA; no retrieval or summary")
    p.add_argument("--source", required=True); p.add_argument("--source-sha256", required=True)
    p.add_argument("--pointer", action="append", default=[])
    p.add_argument("--line-start", type=int); p.add_argument("--line-end", type=int)
    p.add_argument("--max-chars", type=int, default=12000)
    p = sub.add_parser("light-review-batch", help="collect local findings and selected saved facts once; detailed receipts stay on disk")
    p.add_argument("--request", required=True); p.add_argument("--request-sha256", required=True)
    p.add_argument("--out", required=True)
    p = sub.add_parser("light-apply-batch", help="apply already-authored decisions with machine-built hashes and patches")
    p.add_argument("--review", required=True); p.add_argument("--review-sha256", required=True)
    p.add_argument("--decisions", required=True); p.add_argument("--decisions-sha256", required=True)
    p.add_argument("--out", required=True)
    p = sub.add_parser("light-postprocess", help="check and convert saved original Light answers with compact artifact references; no model call")
    p.add_argument("--request", required=True); p.add_argument("--request-sha256", required=True); p.add_argument("--out", required=True)
    p.add_argument("--decisions", default=""); p.add_argument("--decisions-sha256", default="")
    p.add_argument("--cache", default=""); p.add_argument("--offline", action="store_true")
    args = parser.parse_args(argv)
    if args.command == "preflight":
        try:
            result = preflight(args.prompts, args.response, case_id=args.case_id, executor=args.executor,
                               model=args.model, reading_manifest=args.reading_manifest)
        except (OSError, ValueError) as exc:
            result = {"status": "rejected", "errors": [{"file": args.prompts, "path": "$", "reason": str(exc)}]}
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return int(bool(result["errors"]))
    if args.command == "summarize-initial":
        if __package__:
            from .agentic_route_b_initial import summarize_initial
        else:
            from agentic_route_b_initial import summarize_initial
        result = summarize_initial(args.bundle, args.bundle_sha256, args.assessments, args.out)
    elif args.command == "merge-initial":
        if __package__:
            from .agentic_route_b_research import merge_initial_batches
        else:
            from agentic_route_b_research import merge_initial_batches
        result = merge_initial_batches(args.bundle, args.bundle_sha256, args.batches,
                                       args.batches_sha256, args.out)
    elif args.command == "plan-light-followup":
        if __package__:
            from .agentic_route_b_triage import plan_light_followup
        else:
            from agentic_route_b_triage import plan_light_followup
        result = plan_light_followup(args.summary, args.screening, args.out)
    elif args.command == "light-prepare":
        if __package__:
            from .agentic_route_b_light import light_prepare
        else:
            from agentic_route_b_light import light_prepare
        result = light_prepare(args.bundle, args.bundle_sha256, args.out,
                               batch_units=args.batch_units, units_file=args.units_file,
                               output_contract=args.output_contract, packing=args.packing)
    elif args.command == "light-convert":
        if __package__:
            from .agentic_route_b_light import light_convert
        else:
            from agentic_route_b_light import light_convert
        result = light_convert(args.bundle, args.bundle_sha256, args.light_index,
                               args.light_index_sha256, args.out,
                               cache_arg=args.cache, offline=args.offline)
    elif args.command == "light-postprocess":
        if __package__:
            from .agentic_route_b_light_local import light_postprocess
        else:
            from agentic_route_b_light_local import light_postprocess
        result = light_postprocess(args.request, args.request_sha256, args.out, decisions_arg=args.decisions,
                                   decisions_sha256=args.decisions_sha256, cache_arg=args.cache, offline=args.offline)
    elif args.command in {"light-review-batch", "light-apply-batch"}:
        if __package__:
            from .agentic_route_b_light_local import light_review_batch, light_apply_batch
        else:
            from agentic_route_b_light_local import light_review_batch, light_apply_batch
        if args.command == "light-review-batch":
            result = light_review_batch(args.request, args.request_sha256, args.out)
        else:
            result = light_apply_batch(args.review, args.review_sha256, args.decisions,
                                       args.decisions_sha256, args.out)
    elif args.command in {"light-preflight", "light-apply"}:
        if __package__:
            from .agentic_route_b_light_local import light_preflight, light_apply
        else:
            from agentic_route_b_light_local import light_preflight, light_apply
        bound = (args.bundle, args.bundle_sha256, args.light_index, args.light_index_sha256,
                 args.batch_id, args.results)
        if args.command == "light-preflight":
            result = light_preflight(*bound, args.out)
        else:
            result = light_apply(*bound, args.patch, args.patch_sha256, args.out_results, args.receipt)
    elif args.command in {"light-read-material", "light-read-source"}:
        if __package__:
            from .agentic_route_b_material_read import material_read, source_read
        else:
            from agentic_route_b_material_read import material_read, source_read
        if args.command == "light-read-material":
            result = material_read(args.bundle, args.bundle_sha256, args.light_index, args.light_index_sha256,
                                   batch_id=args.batch_id, unit_ids=args.unit_id, gene_ids=args.gene_id,
                                   fields=args.field, max_chars=args.max_chars)
        else:
            result = source_read(args.source, args.source_sha256, pointers=args.pointer or None,
                                 line_start=args.line_start, line_end=args.line_end, max_chars=args.max_chars)
    elif args.command == "extend":
        if __package__:
            from .agentic_route_b_extend import extend
        else:
            from agentic_route_b_extend import extend
        result = extend(args.bundle, args.bundle_sha256, args.old_prompts, args.discovery,
                        args.old_response, args.out, import_manifest=args.import_manifest,
                        max_material_bytes=args.max_material_bytes)
    elif args.command == "reuse-plan":
        if __package__:
            from .agentic_route_b_reuse import reuse_plan
        else:
            from agentic_route_b_reuse import reuse_plan
        result = reuse_plan(args.old_prompts, args.new_prompts, args.old_response, args.out, import_manifest=args.import_manifest, case_id=args.case_id)
    else:
        result = prepare(args.prompts, args.out, args.max_prompts, case_id=args.case_id, max_material_bytes=args.max_material_bytes, prefer_shared_material=args.prefer_shared_material, research_bundle=args.research_bundle, research_bundle_sha256=args.research_bundle_sha256, followup_plan=args.followup_plan, followup_plan_sha256=args.followup_plan_sha256, research_ledger=args.research_ledger, research_ledger_sha256=args.research_ledger_sha256, research_batch_units=args.research_batch_units) if args.command == "prepare" else merge(args.prompts, args.response, args.out, executor=args.executor, model=args.model, case_id=args.case_id, reading_manifest=args.reading_manifest)
    summary = {k: v for k, v in result.items() if k not in {"shards", "partial_response_provenance", "retry_history", "response_roster", "rows"}}
    if "reading_view" in summary:
        summary["reading_view"] = {k: v for k, v in summary["reading_view"].items() if k != "shards"}
        summary["reading_view"]["shard_count"] = len(result["reading_view"]["shards"])
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
