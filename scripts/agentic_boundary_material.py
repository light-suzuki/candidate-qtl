"""Compact, hash-bound reading views for canonical boundary shards.

The shard JSON files remain the audit source. These derived files intern repeated
material and keep per-gene/copy bindings explicit for bounded reading.
"""
from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any


DEFAULT_MAX_MATERIAL_BYTES = 262_144
STRUCTURED_READING = "structured_json_catalog_refs_v1"


def canonical_json(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")


def sha_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _write(path: Path, value: Any) -> str:
    data = json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(data, encoding="utf-8")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def material_id(material: Any) -> str:
    return sha_bytes(canonical_json(material))


def _unique_json_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate JSON key cannot be losslessly projected")
        result[key] = value
    return result


def _invalid_json_constant(value: str) -> None:
    raise ValueError(f"non-JSON constant: {value}")


def paper_reading_projection(value: Any, path: str = "", *, structured_json: bool = False) -> Any:
    """Read all distinct-paper summaries; keep lengthy raw ledgers in shards."""
    if isinstance(value, list):
        return [paper_reading_projection(child, f"{path}[{index}]", structured_json=structured_json) for index, child in enumerate(value)]
    if not isinstance(value, dict):
        return value
    result = {}
    for key, child in value.items():
        if structured_json and key.endswith("_json") and isinstance(child, str):
            try:
                parsed = json.loads(child, object_pairs_hook=_unique_json_object, parse_constant=_invalid_json_constant)
                if isinstance(parsed, (dict, list)):
                    child = parsed
            except ValueError:
                pass  # Invalid/non-container JSON remains the original literal.
        result[key] = paper_reading_projection(child, f"{path}.{key}" if path else key, structured_json=structured_json)
    if (value.get("tier_aggregation") == "strongest_paper_v1"
            and value.get("relevance_tier") in {"T1_reported_focal_trait", "T2_reported_similar_phenomenon", "T3_reported_shared_process"}
            and value.get("tier_evidence_json") not in (None, "", "{}")
            and value.get("paper_evidence_summary_json") not in (None, "", "[]")):
        ledger = value.get("relevance_evidence_json", "[]")
        if isinstance(ledger, str):
            try:
                records = json.loads(ledger)
            except ValueError:
                return result  # Unparseable evidence is not safe to summarize.
            cautions = set()
            for record in records if isinstance(records, list) else []:
                if not isinstance(record, dict):
                    continue
                for key in ("counterevidence", "contradiction_reason"):
                    if record.get(key):
                        cautions.add(str(record[key]))
                if record.get("verification_status") == "contradicted":
                    cautions.add(str(record.get("verification_reason") or "Recorded contradicted evidence; inspect audit ledger"))
            if cautions:
                result["known_counterevidence"] = sorted(cautions)
        for key in ("relevance_evidence_json", "paper_reports_json", "source_ids"):
            if key not in value:
                continue
            if isinstance(value[key], dict) and value[key].get("audit_only"):
                continue
            result[key] = {"audit_only": True, "sha256": sha_bytes(canonical_json(value[key])),
                           "full_value_in": f"{path}.{key}" if path else key,
                           "read_when": "specific unresolved decision or known conflict"}
        result["paper_reading_policy"] = "all_distinct_paper_summaries; counts are not independent replications; full records in canonical shard"
    return result


def candidate_material_reading_projection(materials: dict[str, Any], *, structured_json: bool = False) -> dict[str, Any]:
    """Project each candidate at its own root, as in the original reading views.

    ``full_value_in`` is relative to the candidate's canonical material, not
    the ID-keyed envelope. Keep this root identical when writing and checking
    views; projecting the envelope would silently prefix audit paths with IDs.
    """
    return {cid: paper_reading_projection(material, structured_json=structured_json) for cid, material in materials.items()}


def expand_reading_view(view: dict[str, Any], catalog: dict[str, Any]) -> dict[str, Any]:
    """Resolve existing catalog bindings, without interpreting scientific markers."""
    result = dict(view)
    result["candidate_material"] = {cid: _expand_tree(binding, catalog) for cid, binding in view.get("candidate_material", {}).items()}
    if view.get("shared_field_bindings"):
        result["candidate_input"] = {cid: _expand_tree(binding, catalog) for cid, binding in view["candidate_input"].items()}
        for key in ("higher_tier_context", "primary_paper_evidence", "request_policy", "comparison_scope"):
            if key in view:
                result[key] = _expand_tree(view[key], catalog)
    return result


def _intern_tree(value: Any, catalog: dict[str, Any], counts: Counter | None = None) -> dict[str, Any]:
    """Intern exact substantial subtrees; inline small values to avoid hash overhead."""
    if isinstance(value, str) and len(canonical_json(value)) >= 128:
        node = {"kind": "scalar", "value": value}
        key = sha_bytes(canonical_json(node))
        catalog.setdefault(key, node)
        if counts is not None:
            count = counts[material_id(value)]
            if count * len(canonical_json({"$inline": value})) <= (
                    len(canonical_json(node)) + len(key) + 4 + count * len(canonical_json({"$ref": key}))):
                return {"$inline": value}
        return {"$ref": key}
    if not isinstance(value, (dict, list)) or len(canonical_json(value)) < 128:
        return {"$inline": value}
    if isinstance(value, dict):
        node = {"kind": "dict", "value": {key: _intern_tree(child, catalog, counts)
                                             for key, child in sorted(value.items())}}
    else:
        node = {"kind": "list", "value": [_intern_tree(child, catalog, counts) for child in value]}
    key = sha_bytes(canonical_json(node))
    catalog.setdefault(key, node)
    count = counts[material_id(value)] if counts is not None else 0
    if counts is not None and count * len(canonical_json({"$inline": value})) <= (
            len(canonical_json(node)) + len(key) + 4 + count * len(canonical_json({"$ref": key}))):
        return {"$inline": value}
    return {"$ref": key}


def _count_subtrees(value: Any, counts: Counter) -> None:
    if isinstance(value, (dict, list, str)) and len(canonical_json(value)) >= 128:
        key = material_id(value)
        counts[key] += 1
        if counts[key] > 1:
            return  # Reused parents do not duplicate their catalog child edges.
    if isinstance(value, dict):
        for child in value.values():
            _count_subtrees(child, counts)
    elif isinstance(value, list):
        for child in value:
            _count_subtrees(child, counts)


def candidate_material_bytes(material: Any) -> int:
    """UTF-8 size of full material payload, excluding JSON envelope overhead."""
    return len(canonical_json(material))


def _expand_tree(binding: dict[str, Any], catalog: dict[str, Any]) -> Any:
    if "$inline" in binding:
        return binding["$inline"]
    node = catalog[binding["$ref"]]
    kind, value = node["kind"], node["value"]
    if kind == "scalar":
        return value
    if kind == "dict":
        return {name: _expand_tree(ref, catalog) for name, ref in value.items()}
    if kind == "list":
        return [_expand_tree(ref, catalog) for ref in value]
    raise ValueError("invalid shared material node")


def _reachable_keys(binding: dict[str, Any], catalog: dict[str, Any], found: set[str]) -> None:
    key = binding.get("$ref")
    if key is None or key in found:
        return
    found.add(key)
    node = catalog[key]
    value = node["value"]
    if node["kind"] == "dict":
        for ref in value.values():
            _reachable_keys(ref, catalog, found)
    elif node["kind"] == "list":
        for ref in value:
            _reachable_keys(ref, catalog, found)


def paper_evidence_for_candidates(panels: list[dict[str, Any]], candidate_ids: list[str]) -> list[dict[str, Any]]:
    """Project selected papers onto assigned copies without ranking or reselection."""
    assigned = set(candidate_ids)
    result = []
    for panel in panels:
        panel_ids = [cid for cid in panel["candidate_gene_ids"] if cid in assigned]
        if not panel_ids:
            continue
        papers = []
        for paper in panel["papers"]:
            paper_ids = [cid for cid in paper["candidate_gene_ids"] if cid in panel_ids]
            if not paper_ids:
                continue
            reports = []
            for report in paper["reports"]:
                report_ids = [cid for cid in report["candidate_gene_ids"] if cid in paper_ids]
                if report_ids:
                    reports.append({**report, "candidate_gene_ids": report_ids})
            papers.append({**paper, "candidate_gene_ids": paper_ids, "reports": reports,
                           "record_state": "original_reports" if reports else "original_report_unavailable"})
        if papers:
            result.append({**panel, "candidate_gene_ids": panel_ids, "papers": papers})
    return result


def partition_candidate_ids(ids: list[str], materials: dict[str, Any], max_material_bytes: int,
                            max_candidates: int, *, candidate_rows: dict[str, dict[str, Any]] | None = None,
                            fixed_bytes: int = 0) -> tuple[list[list[str]], list[str]]:
    if max_material_bytes < 1:
        raise ValueError("max_material_bytes must be positive")
    groups: list[list[str]] = []
    oversized: list[str] = []
    group: list[str] = []
    group_keys: set[str] = set()
    # If the invariant shared index/context alone exceeds the reading budget,
    # keep normal candidate groups and report the common overage separately.
    fixed_for_partition = fixed_bytes if fixed_bytes <= max_material_bytes else 0
    used = fixed_for_partition
    for cid in ids:
        value = materials[cid]
        key = material_id(value)
        material_size = candidate_material_bytes(value)
        row_size = len(canonical_json(candidate_rows[cid])) if candidate_rows else 0
        incremental = (0 if key in group_keys else material_size) + row_size
        if group and (len(group) >= max_candidates or used + incremental > max_material_bytes):
            groups.append(group)
            group, group_keys, used = [], set(), fixed_for_partition
            incremental = material_size + row_size
        group.append(cid)
        if key not in group_keys:
            used += material_size
            group_keys.add(key)
        used += row_size
        if used > max_material_bytes:
            oversized.append(cid)
            groups.append(group)
            group, group_keys, used = [], set(), fixed_for_partition
    if group:
        groups.append(group)
    return groups, oversized


def write_material_artifacts(root: Path, shards: list[dict[str, Any]], *, context: dict[str, Any] | None,
                             all_ids: list[str], higher_ids: list[str], max_material_bytes: int,
                             oversized: list[str], candidate_rows: dict[str, dict[str, Any]],
                             higher_rows: list[dict[str, Any]], phenotype: str, request_policy: Any,
                             shared_context_bytes: int) -> dict[str, Any]:
    """Write interned catalog and per-shard reading views; return manifest metadata."""
    catalog: dict[str, Any] = {}
    bindings: dict[str, str] = {}
    owner_map = {}
    first_shard = shards[0] if shards else None
    for shard in shards:
        for cid in shard["candidate_gene_ids"]:
            owner_map[cid] = {"file": shard["file"], "sha256": shard["sha256"]}
    if first_shard:
        for cid in higher_ids:
            owner_map[cid] = {"file": first_shard["file"], "sha256": first_shard["sha256"]}
    index = build_comparison_index(all_ids, higher_ids, candidate_rows, higher_rows, context, owner_map, structured_json=True)
    index_path = root / "comparison_index.json"
    index_sha = _write(index_path, {"schema_version": 1,
                                    "projection_policy": "Selected relevance, family, mechanism, support, copy, source, evidence, status, and function fields; follow canonical_shard file/SHA for all omitted or full values.",
                                    "comparison_index": index})
    index_file_bytes = index_path.stat().st_size
    counts: Counter = Counter()
    projected_material = (candidate_material_reading_projection(context["candidate_material"], structured_json=True)
                          if context is not None else {})
    if context is not None:
        for shard in shards:
            for cid in shard["candidate_gene_ids"] + higher_ids:
                _count_subtrees(projected_material[cid], counts)
    for shard in shards:
        item = json.loads((root / shard["file"]).read_text(encoding="utf-8"))
        for row in item["candidate_input"]:
            _count_subtrees(paper_reading_projection(row, structured_json=True), counts)
        for value in (item["higher_tier_context"], request_policy,
                      {"boundary_candidate_gene_ids": all_ids, "higher_tier_read_only_gene_ids": higher_ids}):
            _count_subtrees(paper_reading_projection(value, structured_json=True), counts)
        if "primary_paper_evidence" in item:
            _count_subtrees(paper_reading_projection(item["primary_paper_evidence"], structured_json=True), counts)
    if context is not None:
        for cid in all_ids + higher_ids:
            bindings[cid] = _intern_tree(projected_material[cid], catalog, counts)
    oversized_higher = ([cid for cid in higher_ids
                         if candidate_material_bytes(context["candidate_material"][cid]) > max_material_bytes]
                        if context is not None else [])

    # All scientific reading fields use the same catalog as candidate material.
    # Previously shared higher-tier papers/context were copied into every view.
    shared_fields = {}
    for shard in shards:
        canonical = json.loads((root / shard["file"]).read_text(encoding="utf-8"))
        higher_context = paper_reading_projection(canonical["higher_tier_context"], structured_json=True)
        fields = {"candidate_input": {cid: _intern_tree(paper_reading_projection(candidate_rows[cid], structured_json=True), catalog, counts)
                                        for cid in canonical["candidate_gene_ids"]},
                  "higher_tier_context": _intern_tree(higher_context, catalog, counts),
                  "request_policy": _intern_tree(request_policy, catalog, counts),
                  "comparison_scope": _intern_tree({"boundary_candidate_gene_ids": all_ids,
                                                     "higher_tier_read_only_gene_ids": higher_ids}, catalog, counts)}
        if "primary_paper_evidence" in canonical:
            fields["primary_paper_evidence"] = _intern_tree(
                paper_reading_projection(canonical["primary_paper_evidence"], structured_json=True), catalog, counts)
        shared_fields[shard["shard_id"]] = fields
    reachable: set[str] = set()
    for binding in bindings.values():
        _reachable_keys(binding, catalog, reachable)
    for fields in shared_fields.values():
        for binding in fields["candidate_input"].values():
            _reachable_keys(binding, catalog, reachable)
        for key, binding in fields.items():
            if key != "candidate_input":
                _reachable_keys(binding, catalog, reachable)
    catalog = {key: catalog[key] for key in sorted(reachable)}
    catalog_path = root / "material_catalog.json"
    catalog_sha = _write(catalog_path, {"schema_version": 1, "materials": catalog})
    catalog_bytes = catalog_path.stat().st_size
    view_entries = []
    expanded_read_bytes = 0
    for shard in shards:
        shard_path = root / shard["file"]
        canonical = json.loads(shard_path.read_text(encoding="utf-8"))
        ids = canonical["candidate_gene_ids"]
        shard_bindings = {cid: bindings[cid] for cid in ids + higher_ids} if context is not None else {}
        shard_material_ids = set()
        for binding in shard_bindings.values():
            _reachable_keys(binding, catalog, shard_material_ids)
        view = {"schema_version": 1, "shard_id": shard["shard_id"],
                "canonical_shard": {"path": shard["file"], "sha256": shard["sha256"]},
                "catalog": {"path": "material_catalog.json", "sha256": catalog_sha},
                "request_policy": request_policy,
                "phenotype": phenotype if context is not None else None,
                "shared_field_bindings": True,
                **shared_fields[shard["shard_id"]],
                "candidate_material": shard_bindings,
                "comparison_index": {"path": "comparison_index.json", "sha256": index_sha},
                "reader_instruction": "Resolve candidate_material inline/ref bindings through catalog.materials recursively. Load the shared comparison_index.json once and compare every candidate/copy distinctly. Read paper_evidence_summary_json across Tiers and paper_tier_counts_json with tier_evidence_json; do not treat counts as independent replications or stop at a single T1. Raw audit_only fields remain in the hash-bound canonical shard for citation/content checks when needed. Other full fields are available through each index row canonical_shard file/SHA descriptor. Read oversized material or shared index in chunks without truncating or excluding candidates."}
        view["reader_instruction"] += (
            " Resolve candidate_input rows, higher_tier_context, primary_paper_evidence, request_policy and comparison_scope "
            "through the same catalog bindings as candidate_material. Structured *_json fields contain the complete decoded "
            "JSON value; exact original strings remain in the canonical shard. Read each shared catalog node once, preserving "
            "every candidate-specific binding. "
            " Follow request_policy when comparing experiment priority. Keep the original Tier separate from "
            "source-to-copy correspondence; inspect candidate-specific function, aligned regions, available domains "
            "and conflicting identities in candidate_material. Source and unit-wide context are not copy-specific "
            "functional proof. A reasoned defer postpones an experiment, without proving noncausality; explain the "
            "comparative reason and the question that would change it. Missing material or a partial alignment "
            "alone does not justify automatic deferral. Do not invent observations, force a quota, change Tier, "
            "or add a Tier reviewer."
        )
        if "primary_paper_evidence" in canonical or "primary_paper_evidence" in canonical["higher_tier_context"]:
            view["reader_instruction"] = view["reader_instruction"].replace(
                "Read paper_evidence_summary_json across Tiers and paper_tier_counts_json with tier_evidence_json; do not treat counts as independent replications or stop at a single T1.",
                "Read every retained primary_paper_evidence report first, with its original Tier, reason, and candidate binding. "
                "The stored ranks identify the existing handoff set; do not replace it with a newly ranked five. "
                "Read the other paper summaries for a concrete unresolved question and retain known counterevidence; "
                "counts are not independent replications or a priority bonus."
            )
        oversized_here = [cid for cid in ids if cid in oversized] + [cid for cid in higher_ids if cid in oversized_higher]
        if oversized_here:
            view["oversized_material_gene_ids"] = oversized_here
            view["chunk_reading_required"] = True
        relative = str((Path(shard["file"]).parent / "reading_view.json").as_posix())
        view_path = root / relative
        view_sha = _write(view_path, view)
        expanded = len(canonical_json(expand_reading_view(view, catalog))) + index_file_bytes
        if expanded > max_material_bytes:
            view["chunk_reading_required"] = True
            view["reading_budget_exceeded_bytes"] = expanded - max_material_bytes
            view_sha = _write(view_path, view)
            entry_oversized = True
        else:
            entry_oversized = False
        expanded_read_bytes += expanded
        entry = {"shard_id": shard["shard_id"], "file": relative, "sha256": view_sha,
                 "expanded_read_bytes": expanded,
                 "reading_budget_exceeded": entry_oversized,
                 "oversized_material_gene_ids": [cid for cid in ids if cid in oversized]}
        view_entries.append(entry)

    canonical_bytes = sum((root / shard["file"]).stat().st_size for shard in shards)
    unique_shared_bytes = sum(len(canonical_json(value)) for value in catalog.values())
    index_bytes = index_file_bytes
    expanded_total = sum(item["expanded_read_bytes"] for item in view_entries)
    return {"material_reading_projection": STRUCTURED_READING,
            "material_catalog": {"file": "material_catalog.json", "sha256": catalog_sha},
            "comparison_index_artifact": {"file": "comparison_index.json", "sha256": index_sha,
                                          "bytes": index_file_bytes},
            "reading_views": view_entries,
            "material_delivery_bytes": {"canonical_shard_bytes": canonical_bytes,
                                        "unique_shared_material_bytes": unique_shared_bytes,
                                        "material_catalog_bytes": catalog_bytes,
                                        "comparison_index_bytes": index_file_bytes,
                                        "shard_reading_view_bytes": sum((root / item["file"]).stat().st_size for item in view_entries),
                                        "expanded_read_bytes": expanded_total,
                                        "comparison_index_bytes_per_shard": index_bytes,
                                        "expanded_read_bytes_basis": "compact JSON reading view with shared material references expanded",
                                        "max_material_bytes": max_material_bytes,
                                        "shared_context_bytes_per_shard": shared_context_bytes,
                                        "shared_context_exceeds_budget": shared_context_bytes > max_material_bytes,
                                        "oversized_gene_ids": oversized,
                                        "oversized_higher_tier_gene_ids": oversized_higher,
                                        "oversized_shard_ids": [item["shard_id"] for item in view_entries if item["reading_budget_exceeded"]]}}


_INDEX_TERMS = ("relevance", "tier", "mechanism", "family", "copy", "source", "support", "evidence", "status", "function")


def _compact_value(value: Any, path: str) -> Any:
    if isinstance(value, str) and len(value) > 160:
        return {"preview": value[:157] + "...", "full_value_sha256": sha_bytes(value.encode("utf-8")),
                "full_value_in": f"canonical_shards[candidate_gene_id].candidate_input_or_candidate_comparison_material[{path}]"}
    if isinstance(value, dict):
        return {key: _compact_value(child, f"{path}.{key}") for key, child in value.items()}
    if isinstance(value, list):
        return [_compact_value(child, f"{path}[{index}]") for index, child in enumerate(value)]
    return value


def _compact_row(row: dict[str, Any]) -> dict[str, Any]:
    return {key: _compact_value(value, key) for key, value in paper_reading_projection(row).items()}


def _compact_index_value(value: Any, path: str) -> Any:
    if isinstance(value, str) and len(value) > 64:
        return {"preview": value[:23] + "...", "sha256": sha_bytes(value.encode("utf-8")),
                "field": path, "truncated": True}
    if isinstance(value, dict):
        nested = {}
        for key, child in value.items():
            nested.update(_collect_index_fields(child, f"{path}.{key}", key))
        return nested or {path: {"sha256": sha_bytes(canonical_json(value)), "item_count": len(value),
                                 "truncated": True, "field": path}}
    if isinstance(value, list):
        return {path: {"sha256": sha_bytes(canonical_json(value)), "item_count": len(value),
                       "truncated": True, "field": path}}
    return value


def _collect_index_fields(value: Any, path: str, key: str) -> dict[str, Any]:
    if isinstance(value, dict):
        result = {}
        for child_key, child in value.items():
            result.update(_collect_index_fields(child, f"{path}.{child_key}", child_key))
        return result
    if isinstance(value, list):
        if any(term in key.lower() for term in _INDEX_TERMS):
            return {path: {"sha256": sha_bytes(canonical_json(value)), "item_count": len(value),
                           "truncated": True, "field": path}}
        return {}
    if key == "candidate_gene_id" or any(term in key.lower() for term in _INDEX_TERMS):
        return {path: _compact_index_value(value, path)}
    return {}


def build_comparison_index(all_ids: list[str], higher_ids: list[str], candidate_rows: dict[str, dict[str, Any]],
                           higher_rows: list[dict[str, Any]], context: dict[str, Any] | None,
                           canonical_owners: dict[str, dict[str, str]] | None = None, *, structured_json: bool = False) -> list[dict[str, Any]]:
    higher_by_id = {str(row.get("candidate_gene_id", "")): row for row in higher_rows}
    result = []
    for cid, role in [(cid, "boundary") for cid in all_ids] + [(cid, "higher_read_only") for cid in higher_ids]:
        sources = [candidate_rows.get(cid, {}) if role == "boundary" else higher_by_id.get(cid, {})]
        if context is not None and isinstance(context["candidate_material"].get(cid), dict):
            sources.append(context["candidate_material"][cid])
        fields: dict[str, Any] = {}
        for source in sources:
            for key, value in paper_reading_projection(source, structured_json=structured_json).items():
                if key != "candidate_gene_id" and any(term in key.lower() for term in _INDEX_TERMS):
                    fields.update(_collect_index_fields(value, key, key))
        entry = {"candidate_gene_id": cid, "boundary_role": role, "fields": fields,
                 "fields_are_projection": True}
        if canonical_owners and cid in canonical_owners:
            entry["canonical_shard"] = canonical_owners[cid]
        result.append(entry)
    return result


def validate_reading_projections(root: Path, manifest: dict[str, Any], catalog: dict[str, Any] | None = None) -> dict[str, dict[str, Any]]:
    """Check each view against canonical shards and their complete cross-shard index."""
    if not manifest.get("material_catalog"):
        return {}
    projection = manifest.get("material_reading_projection")
    if projection not in (None, STRUCTURED_READING):
        raise ValueError("unknown material reading projection")
    structured = projection == STRUCTURED_READING
    canonical = {item["shard_id"]: json.loads((root / item["file"]).read_text(encoding="utf-8"))
                 for item in manifest["shards"]}
    if catalog is None:
        catalog = json.loads((root / manifest["material_catalog"]["file"]).read_text(encoding="utf-8"))["materials"]
    ids = [cid for item in manifest["shards"] for cid in canonical[item["shard_id"]]["candidate_gene_ids"]]
    candidate_rows = {row["candidate_gene_id"]: row for item in manifest["shards"]
                      for row in canonical[item["shard_id"]]["candidate_input"]}
    first = canonical[manifest["shards"][0]["shard_id"]]
    higher = first["higher_tier_context"]
    higher_ids, higher_rows = higher["candidate_gene_ids"], higher["rows"]
    materials: dict[str, Any] = {}
    for cid in ids + higher_ids:
        owner = next((canonical[item["shard_id"]] for item in manifest["shards"]
                      if cid in item["candidate_gene_ids"]), first)
        material = owner.get("candidate_comparison_material", {}).get(cid)
        if material is not None:
            materials[cid] = material
    context = {"candidate_material": materials} if materials else None
    owner_map = {}
    first_entry = manifest["shards"][0]
    for shard in manifest["shards"]:
        owner = {"file": shard["file"], "sha256": shard["sha256"]}
        for cid in shard["candidate_gene_ids"]:
            owner_map[cid] = owner
    for cid in higher_ids:
        owner_map[cid] = {"file": first_entry["file"], "sha256": first_entry["sha256"]}
    expected_index = build_comparison_index(ids, higher_ids, candidate_rows, higher_rows, context, owner_map, structured_json=structured)
    index_descriptor = manifest["comparison_index_artifact"]
    index_value = json.loads((root / index_descriptor["file"]).read_text(encoding="utf-8"))
    if index_value.get("comparison_index") != expected_index:
        raise ValueError("shared comparison index does not match canonical shard projection")
    expected_policy = first.get("request_policy")
    expected_phenotype = first.get("phenotype")
    views = {item["shard_id"]: json.loads((root / item["file"]).read_text(encoding="utf-8"))
             for item in manifest.get("reading_views", [])}
    if set(views) != set(canonical):
        raise ValueError("reading view roster mismatch")
    for shard in manifest["shards"]:
        item = canonical[shard["shard_id"]]
        view = views[shard["shard_id"]]
        if bool(view.get("shared_field_bindings")) != structured:
            raise ValueError("reading view projection/catalog mode mismatch")
        expanded_view = expand_reading_view(view, catalog)
        if item.get("candidate_gene_ids") != shard.get("candidate_gene_ids"):
            raise ValueError("canonical shard roster mismatch")
        expected_rows = {row["candidate_gene_id"]: paper_reading_projection(row, structured_json=structured) for row in item["candidate_input"]}
        expected_material = candidate_material_reading_projection(item.get("candidate_comparison_material", {}), structured_json=structured)
        restored_material = expanded_view.get("candidate_material", {})
        expected_higher = (paper_reading_projection(higher, structured_json=True) if structured else
                           {"state": "read_only", "candidate_gene_ids": higher_ids, "rows": [_compact_row(row) for row in higher_rows]})
        if not structured and "primary_paper_evidence" in higher:
            expected_higher["primary_paper_evidence"] = higher["primary_paper_evidence"]
        if (expanded_view.get("candidate_input") != expected_rows
                or expanded_view.get("request_policy") != expected_policy
                or view.get("phenotype") != expected_phenotype
                or expanded_view.get("higher_tier_context") != expected_higher
                or ("primary_paper_evidence" in view) != ("primary_paper_evidence" in item)
                or expanded_view.get("primary_paper_evidence") != paper_reading_projection(item.get("primary_paper_evidence"), structured_json=structured)
                or restored_material != expected_material
                or view.get("comparison_index") != {"path": index_descriptor["file"], "sha256": index_descriptor["sha256"]}
                or expanded_view.get("comparison_scope") != {"boundary_candidate_gene_ids": ids,
                                                    "higher_tier_read_only_gene_ids": higher_ids}):
            raise ValueError("reading view does not restore canonical shard fields")
    try:
        from scripts.agentic_boundary_reading import validate_reading_pages
    except ModuleNotFoundError:
        from agentic_boundary_reading import validate_reading_pages
    validate_reading_pages(root, manifest)
    return canonical


def restore_reading_shard(root: Path, manifest: dict[str, Any], shard_id: str) -> dict[str, Any]:
    """Return the exact hash-verified canonical shard after checking all view projections."""
    validate_material_artifacts(root, manifest)
    restored = validate_reading_projections(root, manifest)
    if shard_id not in restored:
        raise ValueError("unknown shard_id or missing material reading view")
    return restored[shard_id]


def validate_material_artifacts(root: Path, manifest: dict[str, Any]) -> None:
    """Reject changed canonical shards, shared catalog, or derived reading views."""
    try:
        from scripts.agentic_boundary_reading import validate_reading_page_generation
    except ModuleNotFoundError:
        from agentic_boundary_reading import validate_reading_page_generation
    validate_reading_page_generation(root, manifest)
    catalog = manifest.get("material_catalog")
    if not catalog:
        return
    catalog_path = (root / catalog["file"]).resolve()
    if not catalog_path.is_relative_to(root.resolve()) or not catalog_path.is_file() or hashlib.sha256(catalog_path.read_bytes()).hexdigest() != catalog["sha256"]:
        raise ValueError("prepared material catalog hash mismatch")
    index = manifest.get("comparison_index_artifact")
    if not index:
        raise ValueError("prepared comparison index descriptor missing")
    index_path = (root / index["file"]).resolve()
    if not index_path.is_relative_to(root.resolve()) or not index_path.is_file() or hashlib.sha256(index_path.read_bytes()).hexdigest() != index["sha256"]:
        raise ValueError("prepared comparison index hash mismatch")
    for shard in manifest["shards"]:
        shard_path = (root / shard["file"]).resolve()
        if not shard_path.is_relative_to(root.resolve()) or not shard_path.is_file() or hashlib.sha256(shard_path.read_bytes()).hexdigest() != shard["sha256"]:
            raise ValueError("prepared canonical shard hash mismatch")
    for view in manifest.get("reading_views", []):
        path = (root / view["file"]).resolve()
        if not path.is_relative_to(root.resolve()) or not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != view["sha256"]:
            raise ValueError("prepared reading view hash mismatch")
        content = json.loads(path.read_text(encoding="utf-8"))
        if content.get("catalog", {}).get("sha256") != catalog["sha256"]:
            raise ValueError("prepared reading view catalog binding mismatch")
        if content.get("comparison_index", {}).get("sha256") != index["sha256"]:
            raise ValueError("prepared reading view comparison index binding mismatch")
        shard = next(item for item in manifest["shards"] if item["shard_id"] == view["shard_id"])
        if content.get("canonical_shard", {}).get("sha256") != shard["sha256"]:
            raise ValueError("prepared reading view canonical binding mismatch")
