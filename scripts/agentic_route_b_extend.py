"""Compose the existing B discovery re-freeze, prepare and hash-only reuse plan."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

if __package__:
    from .agentic_route_b import _prompts, _read_json, _sha, _write_json, prepare
    from .agentic_route_b_reuse import prompt_digest, reuse_plan
    from .agentic_route_b_material import DEFAULT_MAX_MATERIAL_BYTES
else:
    from agentic_route_b import _prompts, _read_json, _sha, _write_json, prepare
    from agentic_route_b_reuse import prompt_digest, reuse_plan
    from agentic_route_b_material import DEFAULT_MAX_MATERIAL_BYTES

CONTRACT = "agentic-route-b-extension-1.0.0"


def extend(bundle_path: str, parent_sha256: str, old_prompts: str,
           discovery_paths: list[str], old_response: str, out_dir: str, *,
           import_manifest: str = "", max_material_bytes: int = DEFAULT_MAX_MATERIAL_BYTES) -> dict[str, Any]:
    from candidate_qtl.route_b_semantic_bundle import refreeze_route_b_sources

    parent_file, out = Path(bundle_path).resolve(), Path(out_dir).resolve()
    if out.exists():
        raise ValueError(f"output already exists: {out}")
    if _sha(parent_file) != parent_sha256:
        raise ValueError("research parent bundle SHA-256 mismatch")
    parent = _read_json(parent_file)
    prompts, case_id = _prompts(Path(old_prompts), parent["case_id"])
    if parent.get("execution_prompt_fingerprints") != [prompt_digest(p) for p in prompts]:
        raise ValueError("research parent bundle does not bind old prompts")
    if out.is_relative_to(parent_file.parent):
        raise ValueError("extension output must be separate from the frozen bundle directory")
    discoveries, receipts, normalizations, seen = [], [], [], set()
    for path in sorted({p.resolve() for p in map(Path, discovery_paths)}, key=str):
        payload = _read_json(path)
        if payload.get("case_id") != case_id or payload.get("parent_bundle_sha256") != parent_sha256:
            raise ValueError("research discovery case/parent hash mismatch")
        found = payload.get("discovered_sources")
        if not isinstance(found, list) or not found:
            raise ValueError("research requires non-empty discovered_sources")
        receipts.append({"file": str(path), "sha256": _sha(path)})
        for index, original in enumerate(found):
            if not isinstance(original, dict):
                raise ValueError("research discovery must be an object")
            row = dict(original)
            aliases = row.get("source_identifier_aliases", "")
            if isinstance(aliases, list):
                if any(not isinstance(value, str) or not value.strip() for value in aliases):
                    raise ValueError("research source aliases must contain nonempty strings")
                row["source_identifier_aliases"] = ";".join(aliases)
                normalizations.append({"file": str(path), "sha256": receipts[-1]["sha256"],
                                       "record_index": index, "field": "source_identifier_aliases",
                                       "operation": "string_array_to_semicolon_text"})
            key = json.dumps(row, ensure_ascii=False, sort_keys=True, allow_nan=False)
            if key not in seen:
                seen.add(key)
                discoveries.append(row)
    if not discoveries:
        raise ValueError("at least one discovery artifact is required")
    # Audit the old raw/checkpoint/import chain before writing a generation.
    reuse_plan(old_prompts, old_prompts, old_response, None, import_manifest=import_manifest, case_id=case_id)
    if not import_manifest:
        from candidate_qtl.route_b_semantic_import import import_route_b_semantic_proposals
        audit_input = out / "prior_import_inputs.json"
        _write_json(audit_input, {"bundle_sha256": parent_sha256,
                                  "prompt_sha256": _sha(Path(old_prompts)),
                                  "response_sha256": _sha(Path(old_response)),
                                  "artifact_role": "existing_checkpoint_hashes_for_import_audit_not_new_model_execution"})
        audited = import_route_b_semantic_proposals(str(parent_file), old_prompts, old_response,
                                                   str(out / "prior_import_audit"),
                                                   llm_manifest_arg=str(audit_input), case_id=case_id)
        import_manifest = str(audited["manifest"])
    merged_path = out / "discoveries.json"
    _write_json(merged_path, {"case_id": case_id, "parent_bundle_sha256": parent_sha256,
                              "discovered_sources": discoveries, "discovery_provenance": receipts,
                              "transport_normalizations": normalizations})
    bundle = refreeze_route_b_sources(str(parent_file), str(merged_path), str(out / "bundle"),
                                      parent_sha256=parent_sha256)
    prepared_dir = out / "prepared"
    # One complete unit per reading assignment lets retry scope exclude all
    # reusable answers without giving workers conflicting partial assignments.
    prepared = prepare(bundle["prompts"], str(prepared_dir), max_prompts=1,
                       max_material_bytes=max_material_bytes, case_id=case_id)
    plan_path = out / "reuse_plan.json"
    plan = reuse_plan(old_prompts, bundle["prompts"], old_response, str(plan_path),
                      import_manifest=import_manifest, case_id=case_id)
    reusable = {row["prompt_id"]: row for row in plan["rows"] if row["reuse_eligible"]}
    groups: dict[str, list[dict[str, Any]]] = {}
    for row in reusable.values():
        groups.setdefault(row["raw_response"]["file"], []).append(row)
    reuse_files = []
    for index, (source, rows) in enumerate(sorted(groups.items())):
        original = _read_json(Path(source))
        ids = {row["prompt_id"] for row in rows}
        receipt = rows[0]["raw_response"]
        target = out / "reusable" / f"responses-{index:04d}.json"
        _write_json(target, {"case_id": case_id, "executor": receipt["executor"], "model": receipt["model"],
                            "responses": [r for r in original["responses"] if r["prompt_id"] in ids],
                            "reuse": {"contract": CONTRACT, "plan": str(plan_path), "plan_sha256": _sha(plan_path),
                                      "original_response": receipt}})
        reuse_files.append({"file": str(target), "sha256": _sha(target), "prompt_ids": sorted(ids),
                            "original_response": receipt})
    pending = {row["prompt_id"] for row in plan["rows"] if row["comparison"] != "removed" and not row["reuse_eligible"]}
    dispatch = [{**shard, "prepare_directory": str(prepared_dir)}
                for shard in prepared["reading_view"]["shards"] if pending.intersection(shard["prompt_ids"])]
    _write_json(out / "redispatch.json", {"case_id": case_id, "prompt_ids": sorted(pending),
                                         "reading_shards": dispatch})
    result = {"contract": CONTRACT, "case_id": case_id,
              "parent_bundle": str(parent_file), "parent_bundle_sha256": parent_sha256,
              "discoveries": str(merged_path), "discovery_provenance": receipts,
              "transport_normalizations": normalizations,
              "bundle": bundle["bundle"], "prompts": bundle["prompts"],
              "reading_manifest": str(prepared_dir / "prepare_manifest.json"),
              "reuse_plan": str(plan_path), "redispatch": str(out / "redispatch.json"),
              "reusable_responses": reuse_files, "reusable_count": len(reusable),
              "new_execution_required_count": len(pending), "new_model_execution_count": 0,
              "next_action": "execute only redispatch reading_shards, then merge reusable responses and new raw replies against prompts and reading_manifest"}
    result["files"] = {str(p.relative_to(out)): _sha(p) for p in sorted(out.rglob("*")) if p.is_file()}
    _write_json(out / "extension_manifest.json", result)
    return result


def validate_reuse(payload: dict[str, Any], prompts_path: str) -> dict[str, dict[str, Any]]:
    """Re-run the existing reuse gate before consuming a transported old answer."""
    receipt = payload.get("reuse")
    if not isinstance(receipt, dict) or receipt.get("contract") != CONTRACT:
        raise ValueError("unsupported response reuse contract")
    path = Path(receipt["plan"])
    if _sha(path) != receipt.get("plan_sha256"):
        raise ValueError("reuse plan hash mismatch")
    plan = _read_json(path)
    if _sha(Path(prompts_path)) != plan["new_prompts_sha256"]:
        raise ValueError("reuse destination prompt hash mismatch")
    refreshed = reuse_plan(plan["old_prompts"], plan["new_prompts"], plan["old_response"], None,
                           import_manifest=(plan.get("import_audit") or {}).get("file", ""))
    if refreshed != plan:
        raise ValueError("reuse source artifacts no longer match the frozen plan")
    eligible = {row["prompt_id"]: row for row in plan["rows"] if row["reuse_eligible"]}
    source = receipt["original_response"]
    original = _read_json(Path(source["file"]))
    originals = {r["prompt_id"]: r for r in original["responses"]}
    records = payload.get("responses")
    if not isinstance(records, list) or not records:
        raise ValueError("reused responses require a non-empty roster")
    for record in records:
        pid = record.get("prompt_id")
        if (pid not in eligible or eligible[pid]["raw_response"] != source or record != originals.get(pid)
                or payload.get("executor") != source["executor"] or payload.get("model") != source["model"]):
            raise ValueError("reused response differs from eligible original raw/provenance")
    requested = {row["prompt_id"] for row in records}
    return {row["prompt_id"]: row for row in _read_json(Path(plan["old_response"]))["candidates"]
            if row["prompt_id"] in requested}
