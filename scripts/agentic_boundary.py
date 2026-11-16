"""Agent-facing shard/merge harness for an immutable Wet boundary review request."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd
try:
    from scripts.agentic_boundary_material import (
        DEFAULT_MAX_MATERIAL_BYTES, build_comparison_index, canonical_json,
        paper_evidence_for_candidates, partition_candidate_ids, validate_reading_projections, validate_material_artifacts,
        write_material_artifacts, _compact_row,
    )
except ModuleNotFoundError:
    from agentic_boundary_material import (
        DEFAULT_MAX_MATERIAL_BYTES, build_comparison_index, canonical_json,
        paper_evidence_for_candidates, partition_candidate_ids, validate_reading_projections, validate_material_artifacts,
        write_material_artifacts, _compact_row,
    )

DECISIONS = {"adopt", "defer", "not_evaluated"}


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _write(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _comparison_context(path: Path, case_id: str, ids: list[str]) -> tuple[dict, dict]:
    raw = path.read_bytes()
    context = json.loads(raw)
    if not isinstance(context, dict) or context.get("case_id") != case_id:
        raise ValueError("comparison context case_id mismatch")
    if not isinstance(context.get("phenotype"), str) or not context["phenotype"].strip():
        raise ValueError("comparison context requires the recorded phenotype")
    material = context.get("candidate_material")
    if (not isinstance(material, dict) or set(material) != set(ids)
            or any(not isinstance(row, dict) for row in material.values())):
        raise ValueError("comparison context roster mismatch")
    return context, {"path": str(path), "sha256": hashlib.sha256(raw).hexdigest()}


def prepare(request_path: str, out_dir: str, max_candidates: int = 8, *,
            comparison_context_path: str = "", max_material_bytes: int = DEFAULT_MAX_MATERIAL_BYTES) -> dict[str, Any]:
    request_file, out = Path(request_path).resolve(), Path(out_dir).resolve()
    if max_candidates < 1:
        raise ValueError("max_candidates must be positive")
    if max_material_bytes < 1:
        raise ValueError("max_material_bytes must be positive")
    if out.exists():
        raise ValueError("output already exists")
    request = _read(request_file)
    if request.get("stage") != "boundary_review_request" or request.get("state") != "awaiting_review":
        raise ValueError("request is not an active boundary_review_request")
    case_id = str(request.get("case_id", ""))
    if not case_id:
        raise ValueError("boundary request requires case_id")
    source = Path(str(request.get("candidate_input", {}).get("path", ""))).resolve()
    expected_hash = str(request.get("source_wet_sha256", ""))
    if not source.is_file() or not expected_hash or _sha(source) != expected_hash or request["candidate_input"].get("sha256") != expected_hash:
        raise ValueError("candidate_input/source_wet_sha256 mismatch")
    rows = pd.read_csv(source, sep="\t", dtype=str).fillna("").to_dict("records")
    ids = [str(row.get("candidate_gene_id", "")).strip() for row in rows]
    if not ids or any(not value for value in ids) or len(ids) != len(set(ids)):
        raise ValueError("boundary candidate_input requires non-empty unique candidate_gene_id")
    if ids != request.get("candidate_gene_ids") or len(ids) != request.get("candidate_count"):
        raise ValueError("candidate_input roster does not match request")
    higher = request.get("higher_tier_context", {})
    if higher.get("state") != "read_only" or not isinstance(higher.get("candidate_gene_ids"), list):
        raise ValueError("higher-tier context must be explicitly read_only")
    higher_path = Path(str(higher.get("candidate_input", {}).get("path", ""))).resolve()
    higher_hash = str(higher.get("candidate_input", {}).get("sha256", ""))
    if not higher_path.is_file() or _sha(higher_path) != higher_hash:
        raise ValueError("higher-tier context hash mismatch")
    higher_rows = pd.read_csv(higher_path, sep="\t", dtype=str).fillna("").to_dict("records")
    higher_ids = [str(r.get("candidate_gene_id", "")).strip() for r in higher_rows]
    if any(not value for value in higher_ids) or len(higher_ids) != len(set(higher_ids)):
        raise ValueError("higher-tier context requires non-empty unique candidate_gene_id")
    if set(ids) & set(higher_ids):
        raise ValueError("higher-tier context must be disjoint from boundary candidates")
    if higher_ids != higher.get("candidate_gene_ids"):
        raise ValueError("higher-tier context roster mismatch")

    context = None
    context_descriptor = None
    derived_context = False
    recorded_context = request.get("comparison_context")
    if recorded_context is not None:
        if not isinstance(recorded_context, dict) or not recorded_context.get("path"):
            raise ValueError("recorded comparison context descriptor is invalid")
        context_file = Path(recorded_context["path"]).resolve()
        if not context_file.is_file():
            raise ValueError("recorded comparison context SHA mismatch")
        context, context_descriptor = _comparison_context(context_file, case_id, ids + higher_ids)
        if context_descriptor["sha256"] != recorded_context.get("sha256"):
            raise ValueError("recorded comparison context SHA mismatch")
    if comparison_context_path:
        context_file = Path(comparison_context_path).resolve()
        additional, additional_descriptor = _comparison_context(context_file, case_id, ids + higher_ids)
        if context is not None and additional_descriptor["sha256"] != context_descriptor["sha256"]:
            if additional["phenotype"] != context["phenotype"]:
                raise ValueError("comparison context recorded phenotype mismatch")
            context = {**context, "source_contexts": [context_descriptor, additional_descriptor],
                       "candidate_material": {cid: {
                           "recorded_source_context": context["candidate_material"][cid],
                           "additional_source_context": additional["candidate_material"][cid],
                       } for cid in ids + higher_ids}}
            derived_context = True
        elif context is None:
            context, context_descriptor = additional, additional_descriptor

    row_by_id = dict(zip(ids, rows))
    partition_rows = row_by_id
    if "primary_paper_evidence" in request:
        partition_rows = {
            cid: {**row, "primary_paper_evidence": paper_evidence_for_candidates(request["primary_paper_evidence"], [cid])}
            for cid, row in row_by_id.items()
        }
    higher_papers = (paper_evidence_for_candidates(higher["primary_paper_evidence"], higher_ids)
                     if "primary_paper_evidence" in higher else None)
    comparison_index = build_comparison_index(ids, higher_ids, row_by_id, higher_rows, context)
    for entry in comparison_index:
        entry["canonical_shard"] = {"file": "boundary-0000/shard.json", "sha256": "0" * 64}
    fixed_read_bytes = (len(canonical_json(comparison_index))
                        + len(canonical_json({"rows": [_compact_row(row) for row in higher_rows]}))
                        + len(canonical_json(request.get("policy", [])))
                        + (len(canonical_json(context["phenotype"])) if context is not None else 4)
                        + 512)
    if higher_papers is not None:
        fixed_read_bytes += len(canonical_json({"primary_paper_evidence": higher_papers}))
    if context is not None:
        higher_material_sizes = {hashlib.sha256(canonical_json(context["candidate_material"][cid])).hexdigest():
                                 len(canonical_json(context["candidate_material"][cid])) for cid in higher_ids}
        fixed_read_bytes += sum(higher_material_sizes.values())
    partition_materials = context["candidate_material"] if context is not None else {cid: {} for cid in ids}
    groups, oversized = partition_candidate_ids(ids, partition_materials, max_material_bytes, max_candidates,
                                                candidate_rows=partition_rows, fixed_bytes=fixed_read_bytes)
    out.mkdir(parents=True)
    if derived_context:
        context_file = out / "comparison_context.json"
        _write(context_file, context)
        context_descriptor = {"path": str(context_file), "sha256": _sha(context_file)}
    shards = []
    for index, selected_ids in enumerate(groups):
        sid = f"boundary-{index:04d}"
        selected = [row_by_id[cid] for cid in selected_ids]
        item = {"shard_id": sid, "case_id": case_id, "source_wet_sha256": expected_hash,
                "candidate_gene_ids": selected_ids,
                "candidate_input": selected,
                "higher_tier_context": {"state": "read_only", "candidate_gene_ids": higher_ids, "rows": higher_rows}}
        item["request_policy"] = request.get("policy", [])
        if request.get("assessment_mode"):
            item["assessment_mode"] = request["assessment_mode"]
        if "primary_paper_evidence" in request:
            item["primary_paper_evidence"] = paper_evidence_for_candidates(request["primary_paper_evidence"], selected_ids)
        if higher_papers is not None:
            item["higher_tier_context"]["primary_paper_evidence"] = higher_papers
        item["request_path"] = str(request_file)
        item["request_sha256"] = _sha(request_file)
        if context is not None:
            item["comparison_context"] = context_descriptor
            item["phenotype"] = context["phenotype"]
            item["candidate_comparison_material"] = {
                cid: context["candidate_material"][cid]
                for cid in item["candidate_gene_ids"] + higher_ids
            }
        path = out / sid / "shard.json"
        _write(path, item)
        shards.append({"shard_id": sid, "file": str(path.relative_to(out)), "candidate_gene_ids": item["candidate_gene_ids"], "sha256": _sha(path)})
    manifest = {"schema_version": 1, "case_id": case_id, "request_path": str(request_file),
                "request_sha256": _sha(request_file), "source_wet_path": str(source),
                "source_wet_sha256": expected_hash, "higher_tier_state": "read_only",
                "higher_tier_source_sha256": higher_hash, "candidate_count": len(ids), "shards": shards}
    if request.get("assessment_mode"):
        manifest["assessment_mode"] = request["assessment_mode"]
    material_meta = write_material_artifacts(out, shards, context=context, all_ids=ids, higher_ids=higher_ids,
                                             max_material_bytes=max_material_bytes, oversized=oversized,
                                             candidate_rows=row_by_id, higher_rows=higher_rows,
                                             phenotype=context["phenotype"] if context is not None else "",
                                             request_policy=request.get("policy", []),
                                             shared_context_bytes=fixed_read_bytes)
    manifest.update(material_meta)
    try:
        from scripts.agentic_boundary_reading import CONTRACT as READING_DELIVERY_CONTRACT, write_reading_pages
    except ModuleNotFoundError:
        from agentic_boundary_reading import CONTRACT as READING_DELIVERY_CONTRACT, write_reading_pages
    if context_descriptor is not None:
        manifest["comparison_context"] = context_descriptor
    manifest["reading_pages"] = write_reading_pages(out, manifest)
    manifest["reading_delivery_contract"] = READING_DELIVERY_CONTRACT
    _write(out / "manifest.json", manifest)
    return manifest


def merge(prepared_dir: str, response_paths: list[str], out_path: str, *, role: str,
          executor_id: str = "", executor_model: str = "", reviewer_id: str = "", reviewer_model: str = "",
          proposal_path: str = "") -> dict[str, Any]:
    root, out = Path(prepared_dir).resolve(), Path(out_path).resolve()
    manifest_path = root / "manifest.json"
    manifest = _read(manifest_path)
    validate_material_artifacts(root, manifest)
    if manifest.get("material_catalog"):
        catalog_data = _read(root / manifest["material_catalog"]["file"]).get("materials", {})
        validate_reading_projections(root, manifest, catalog_data)
    request_path = Path(manifest["request_path"])
    source_path = Path(manifest["source_wet_path"])
    if _sha(request_path) != manifest["request_sha256"]:
        raise ValueError("boundary request SHA mismatch")
    request = _read(request_path)
    if request.get("assessment_mode") == "original_assessment_v1" and role != "proposal":
        raise ValueError("original assessment uses proposal role without independent review")
    if request.get("case_id") != manifest["case_id"]:
        raise ValueError("case_id mismatch")
    if _sha(source_path) != manifest["source_wet_sha256"] or request.get("source_wet_sha256") != manifest["source_wet_sha256"]:
        raise ValueError("source SHA mismatch")
    context_descriptor = manifest.get("comparison_context")
    if context_descriptor and _sha(Path(context_descriptor["path"])) != context_descriptor["sha256"]:
        raise ValueError("comparison context SHA mismatch")
    if role not in {"proposal", "review"}:
        raise ValueError("role must be proposal or review")
    if role == "proposal" and (not executor_id.strip() or not executor_model.strip()):
        raise ValueError("proposal requires executor_id and executor_model")
    if role == "review" and (not reviewer_id.strip() or not reviewer_model.strip() or not proposal_path):
        raise ValueError("review requires reviewer identity/model and --proposal")
    proposal_sha = _sha(Path(proposal_path).resolve()) if proposal_path else ""
    if role == "review" and context_descriptor:
        proposal_manifest_path = Path(proposal_path).resolve().with_suffix(Path(proposal_path).suffix + ".manifest.json")
        if not proposal_manifest_path.is_file():
            raise ValueError("comparison review requires the merged proposal manifest")
        proposal_manifest = _read(proposal_manifest_path)
        if (proposal_manifest.get("role") != "proposal"
                or proposal_manifest.get("output_sha256") != proposal_sha
                or proposal_manifest.get("case_id") != manifest["case_id"]
                or proposal_manifest.get("source_wet_sha256") != manifest["source_wet_sha256"]
                or proposal_manifest.get("comparison_context", {}).get("sha256") != context_descriptor["sha256"]):
            raise ValueError("proposal comparison context or lineage mismatch")
    assigned: dict[str, str] = {}
    for shard in manifest["shards"]:
        path = (root / shard["file"]).resolve()
        if (not path.is_relative_to(root)
                or (not manifest.get("material_catalog") and _sha(path) != shard["sha256"])):
            raise ValueError("prepared shard hash mismatch")
        for cid in shard["candidate_gene_ids"]:
            if cid in assigned:
                raise ValueError("duplicate prepared candidate ID")
            assigned[cid] = shard["shard_id"]
    observed: dict[str, dict[str, Any]] = {}
    provenance = []
    for arg in response_paths:
        path = Path(arg).resolve()
        payload = _read(path)
        if not isinstance(payload, dict) or str(payload.get("case_id", "")) != manifest["case_id"]:
            raise ValueError(f"case_id mismatch in response: {path}")
        if payload.get("source_wet_sha256") != manifest["source_wet_sha256"]:
            raise ValueError(f"source SHA mismatch in response: {path}")
        if context_descriptor and payload.get("comparison_context_sha256") != context_descriptor["sha256"]:
            raise ValueError(f"comparison context SHA mismatch in response: {path}")
        if role == "review" and payload.get("proposal_sha256") != proposal_sha:
            raise ValueError(f"proposal SHA mismatch in independent response: {path}")
        if role == "proposal" and (payload.get("executor_id") != executor_id or payload.get("executor_model") != executor_model):
            raise ValueError(f"proposal executor identity mismatch in response: {path}")
        if role == "review" and (payload.get("reviewer_id") != reviewer_id or payload.get("reviewer_model") != reviewer_model):
            raise ValueError(f"reviewer identity mismatch in response: {path}")
        sid = str(payload.get("shard_id", ""))
        shard = next((s for s in manifest["shards"] if s["shard_id"] == sid), None)
        if shard is None:
            raise ValueError(f"unknown shard_id: {sid}")
        provenance.append({"path": str(path), "sha256": _sha(path), "shard_id": sid,
                          "executor_id": payload.get("executor_id", ""), "executor_model": payload.get("executor_model", ""),
                          "reviewer_id": payload.get("reviewer_id", ""), "reviewer_model": payload.get("reviewer_model", ""),
                          "response_id": payload.get("response_id", ""), "status": payload.get("status", "completed")})
        status = payload.get("status", "completed")
        records = payload.get("responses", [])
        if status not in {"completed", "failed", "not_evaluated"}:
            raise ValueError("invalid response status")
        if status == "completed" and not isinstance(records, list):
            raise ValueError("completed response requires responses list")
        if status in {"failed", "not_evaluated"}:
            reason = str(payload.get("reason", "")).strip()
            if not reason:
                raise ValueError("failed/not_evaluated shard requires reason")
            records = [{"candidate_gene_id": cid, "decision": "not_evaluated", "reason": reason} for cid in shard["candidate_gene_ids"]]
        for row in records:
            if not isinstance(row, dict) or not row.get("candidate_gene_id"):
                raise ValueError("every response requires candidate_gene_id")
            cid = str(row["candidate_gene_id"])
            if cid not in assigned:
                raise ValueError(f"unknown candidate ID: {cid}")
            if assigned[cid] != sid:
                raise ValueError(f"candidate assigned to another shard: {cid}")
            if cid in observed:
                raise ValueError(f"duplicate candidate response ID: {cid}")
            decision = str(row.get("decision", ""))
            reason = str(row.get("reason", "")).strip()
            if decision not in DECISIONS:
                raise ValueError(f"invalid decision for {cid}: {decision}")
            if not reason:
                raise ValueError(f"decision requires reason for {cid}")
            observed[cid] = {"decision": decision, "reason": reason, "raw_response": row}
    columns = ["candidate_gene_id", "decision", "reason", "source_wet_sha256"]
    if role == "proposal":
        columns += ["executor_id", "executor_model"]
        identity = (executor_id, executor_model)
    else:
        columns += ["proposal_sha256", "reviewer_id", "reviewer_model"]
        identity = (proposal_sha, reviewer_id, reviewer_model)
    source_sha = manifest["source_wet_sha256"]
    rows_out = []
    for cid in assigned:
        record = observed.get(cid, {"decision": "not_evaluated", "reason": "missing_or_failed_response", "raw_response": None})
        base = {"candidate_gene_id": cid, "decision": record["decision"], "reason": record["reason"], "source_wet_sha256": source_sha}
        if role == "proposal":
            base.update(executor_id=identity[0], executor_model=identity[1])
        else:
            base.update(proposal_sha256=identity[0], reviewer_id=identity[1], reviewer_model=identity[2])
        rows_out.append(base)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows_out, columns=columns).to_csv(out, sep="\t", index=False, lineterminator="\n")
    provenance.sort(key=lambda row: (str(row.get("path", "")), str(row.get("sha256", ""))))
    sidecar = {"schema_version": 1, "role": role, "case_id": manifest["case_id"],
               "source_wet_sha256": source_sha, "output_sha256": _sha(out), "candidate_count": len(rows_out),
               "completion_status": "complete" if all(row["decision"] != "not_evaluated" for row in rows_out) else "partial",
               "state_counts": {state: sum(row["decision"] == state for row in rows_out) for state in sorted(DECISIONS)},
               "proposal_sha256": proposal_sha or None, "partial_response_provenance": provenance,
               "response_records": [{"candidate_gene_id": cid, "raw_response": observed.get(cid, {}).get("raw_response")} for cid in assigned]}
    if request.get("assessment_mode"):
        sidecar["assessment_mode"] = request["assessment_mode"]
    if context_descriptor is not None:
        sidecar["comparison_context"] = context_descriptor
    _write(out.with_suffix(out.suffix + ".manifest.json"), sidecar)
    return sidecar


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare"); p.add_argument("--request", required=True); p.add_argument("--out", required=True); p.add_argument("--max-candidates", type=int, default=8); p.add_argument("--max-material-bytes", type=int, default=DEFAULT_MAX_MATERIAL_BYTES)
    p.add_argument("--comparison-context", default="", help="Hash-bound recorded phenotype and candidate comparison material; does not change support or selection")
    p = sub.add_parser("merge"); p.add_argument("--prepared", required=True); p.add_argument("--response", action="append", default=[]); p.add_argument("--role", choices=["proposal", "review"], required=True); p.add_argument("--out", required=True); p.add_argument("--executor-id", default=""); p.add_argument("--executor-model", default=""); p.add_argument("--reviewer-id", default=""); p.add_argument("--reviewer-model", default=""); p.add_argument("--proposal", default="")
    args = parser.parse_args(argv)
    if args.command == "prepare":
        result = prepare(args.request, args.out, args.max_candidates, comparison_context_path=args.comparison_context, max_material_bytes=args.max_material_bytes)
    else:
        result = merge(args.prepared, args.response, args.out, role=args.role, executor_id=args.executor_id, executor_model=args.executor_model, reviewer_id=args.reviewer_id, reviewer_model=args.reviewer_model, proposal_path=args.proposal)
    print(json.dumps({k: v for k, v in result.items() if k not in {"shards", "response_records", "partial_response_provenance"}}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
