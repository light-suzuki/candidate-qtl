"""Read-only comparison of frozen Route B inputs and existing response receipts."""
from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any


def prompt_digest(prompt: dict[str, Any]) -> str:
    # JSON types remain distinct (True is not 1); list order and user/system
    # strings remain exact. Object key order alone is not prompt content.
    return hashlib.sha256(json.dumps(prompt, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")).hexdigest()


def reuse_plan(old_prompts: str, new_prompts: str, old_response: str, out_path: str | None, *, import_manifest: str = "", case_id: str = "") -> dict[str, Any]:
    if __package__:
        from .agentic_route_b import EXECUTION_CONTRACT, _ModelResponseFormatError, _parsed_model_response, _prompts, _read_json, _response_completed, _sha
    else:
        from agentic_route_b import EXECUTION_CONTRACT, _ModelResponseFormatError, _parsed_model_response, _prompts, _read_json, _response_completed, _sha

    out = Path(out_path) if out_path is not None else None
    if out is not None and out.exists():
        raise ValueError(f"output already exists: {out}")
    old_file, new_file, response_file = map(Path, (old_prompts, new_prompts, old_response))
    old_rows, old_case = _prompts(old_file, case_id)
    new_rows, new_case = _prompts(new_file, case_id)
    if any(not isinstance(p.get(field), str) or not p[field].strip()
           for p in old_rows + new_rows for field in ("prompt_id", "route_b_unit_id", "case_id")):
        raise ValueError("canonical prompt case/prompt/unit IDs must be non-empty strings")
    old = {str(p["prompt_id"]): p for p in old_rows}
    new = {str(p["prompt_id"]): p for p in new_rows}
    manifest_file = response_file.with_name(response_file.stem + "_manifest.json")
    manifest, response = _read_json(manifest_file), _read_json(response_file)
    if (manifest.get("contract") != EXECUTION_CONTRACT or manifest.get("case_id") != old_case
            or response.get("case_id") != old_case or manifest.get("prompt_sha256") != _sha(old_file)
            or manifest.get("response_sha256") != _sha(response_file)):
        raise ValueError("old checkpoint case/input/response hash mismatch")
    candidates = response.get("candidates")
    if not isinstance(candidates, list) or any(not isinstance(r, dict) or not isinstance(r.get("prompt_id"), str) for r in candidates):
        raise ValueError("old checkpoint requires candidates list")
    records = {str(r.get("prompt_id", "")): r for r in candidates}
    if len(records) != len(candidates) or set(records) != set(old):
        raise ValueError("old checkpoint response roster mismatch")
    raw_by_file = {}
    provenance = manifest.get("partial_response_provenance")
    if not isinstance(provenance, list):
        raise ValueError("old checkpoint lacks raw provenance")
    for receipt in provenance:
        path = Path(receipt["file"])
        if _sha(path) != receipt.get("sha256"):
            raise ValueError(f"raw response hash mismatch: {path}")
        raw = _read_json(path)
        executor, model = str(receipt.get("executor", "")).strip(), str(receipt.get("model", "")).strip()
        if not executor or not model or raw.get("case_id") != old_case:
            raise ValueError("raw response lacks executor/model or has a case mismatch")
        if ((raw.get("executor") and raw["executor"] != executor)
                or (raw.get("model") and raw["model"] != model)):
            raise ValueError("raw response executor/model mismatch")
        raw_records = raw.get("responses")
        if not isinstance(raw_records, list) or any(not isinstance(r, dict) or not r.get("prompt_id") for r in raw_records):
            raise ValueError("raw response requires response roster")
        indexed = {str(r["prompt_id"]): r for r in raw_records}
        if len(indexed) != len(raw_records) or not set(indexed) <= set(old) or any(not isinstance(r["prompt_id"], str) for r in raw_records):
            raise ValueError("raw response has duplicate/unknown prompt IDs")
        key = str(path.resolve())
        if key in raw_by_file:
            raise ValueError("duplicate raw provenance file")
        raw_by_file[key] = (indexed, {"file": key, "sha256": _sha(path), "executor": executor, "model": model})
    completed_raw = set()
    for pid, record in records.items():
        if record.get("route_b_unit_id", "") != old[pid]["route_b_unit_id"]:
            raise ValueError("old checkpoint unit mismatch")
        # Missing responses have no source; their state remains not reusable.
        if not record.get("attempt_source_file") and not _response_completed(record):
            continue
        key = str(Path(record.get("attempt_source_file", "")).resolve())
        if key not in raw_by_file or pid not in raw_by_file[key][0]:
            raise ValueError("checkpoint response lacks original raw attempt")
        raw_record, receipt = raw_by_file[key][0][pid], raw_by_file[key][1]
        if (record.get("agent_response_raw") != raw_record or record.get("raw_text", "") != raw_record.get("raw_text", "")
                or record.get("attempt_executor") != receipt["executor"] or record.get("attempt_model") != receipt["model"]
                or raw_record.get("route_b_unit_id", old[pid]["route_b_unit_id"]) != old[pid]["route_b_unit_id"]
                or raw_record.get("case_id", old_case) != old_case):
            raise ValueError("checkpoint differs from original raw response/provenance")
        if raw_record.get("status", "completed") == "completed":
            try:
                if _response_completed(_parsed_model_response(raw_record)):
                    completed_raw.add(pid)
            except _ModelResponseFormatError:
                pass

    rejected: dict[str, list[str]] = {}
    import_receipt = None
    if import_manifest:
        import_file = Path(import_manifest)
        imported = _read_json(import_file)
        if (imported.get("stage") != "semantic_proposal_import" or imported.get("case_id") != old_case
                or imported.get("response_count") != len(old)):
            raise ValueError("import receipt case/roster mismatch")
        for field, path in (("prompt_sha256", old_file), ("response_sha256", response_file)):
            if imported.get(field) and imported[field] != _sha(path):
                raise ValueError("import receipt input hash mismatch")
        rejected_file = Path(imported["outputs"]["rejected_rows"])
        if _sha(rejected_file) != imported.get("rejected_rows_sha256"):
            raise ValueError("import rejected rows hash mismatch")
        rejected_rows = _read_json(rejected_file)
        if not isinstance(rejected_rows, list) or len(rejected_rows) != imported.get("rejected_claim_row_count"):
            raise ValueError("import rejected row count mismatch")
        for row in rejected_rows:
            pid = str(row.get("prompt_id", ""))
            index = row.get("claim_index")
            if pid not in old or row.get("unit_id") != old[pid]["route_b_unit_id"] or type(index) is not int or index < 1:
                raise ValueError("import rejected row prompt/unit mismatch")
            raw_claims = json.loads(records[pid]["raw_text"]).get("claims", [])
            if index > len(raw_claims) or raw_claims[index - 1] != row.get("raw_claim"):
                raise ValueError("import rejection does not match old raw claim")
            rejected.setdefault(pid, []).append(str(row.get("reason", "import_rejected")))
        import_receipt = {"file": str(import_file.resolve()), "sha256": _sha(import_file),
                          "rejected_rows_file": str(rejected_file.resolve()), "rejected_rows_sha256": _sha(rejected_file),
                          "input_binding": "hash_bound" if imported.get("prompt_sha256") and imported.get("response_sha256") else "legacy_case_roster_raw_claim_match"}
    rows = []
    for pid in sorted(set(old) | set(new)):
        old_digest = prompt_digest(old[pid]) if pid in old else ""
        new_digest = prompt_digest(new[pid]) if pid in new else ""
        comparison = "new" if pid not in old else "removed" if pid not in new else "same" if old_digest == new_digest else "changed"
        record = records.get(pid)
        reason = ""
        if comparison != "same":
            reason = f"prompt_{comparison}"
        elif not record or not _response_completed(record) or pid not in completed_raw:
            reason = "response_incomplete_or_invalid"
        elif pid in rejected:
            reason = "known_import_rejected"
        elif import_receipt is None:
            reason = "import_audit_not_recorded"
        elif import_receipt["input_binding"] != "hash_bound":
            reason = "import_audit_input_not_bound"
        key = str(Path(record["attempt_source_file"]).resolve()) if record and record.get("attempt_source_file") else ""
        raw_receipt = raw_by_file[key][1] if key else None
        rows.append({"prompt_id": pid, "old_unit_id": old.get(pid, {}).get("route_b_unit_id", ""),
                     "new_unit_id": new.get(pid, {}).get("route_b_unit_id", ""), "comparison": comparison,
                     "old_prompt_digest": old_digest, "new_prompt_digest": new_digest,
                     "reuse_eligible": not reason, "reason": reason, "import_rejection_reasons": rejected.get(pid, []),
                     "raw_response": raw_receipt,
                     "raw_text_sha256": hashlib.sha256(record["raw_text"].encode("utf-8")).hexdigest()
                     if record and isinstance(record.get("raw_text"), str) else ""})
    counts = Counter(row["comparison"] for row in rows)
    plan = {"schema_version": 1, "mode": "read_only_reuse_plan", "old_case_id": old_case, "new_case_id": new_case,
            "old_prompts": str(old_file.resolve()), "old_prompts_sha256": _sha(old_file),
            "new_prompts": str(new_file.resolve()), "new_prompts_sha256": _sha(new_file),
            "old_response": str(response_file.resolve()), "old_response_sha256": _sha(response_file),
            "execution_manifest": str(manifest_file.resolve()), "execution_manifest_sha256": _sha(manifest_file),
            "import_audit": import_receipt, "comparison_counts": {state: counts[state] for state in ("same", "changed", "new", "removed")},
            "reusable_count": sum(row["reuse_eligible"] for row in rows),
            "new_execution_required_count": sum(row["comparison"] != "removed" and not row["reuse_eligible"] for row in rows),
            "new_model_execution_count": 0, "rows": rows}
    if out is not None:
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("x", encoding="utf-8") as stream:
            stream.write(json.dumps(plan, ensure_ascii=False, indent=2) + "\n")
    return plan
