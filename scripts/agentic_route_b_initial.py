"""Lossless initial reading context and accounting of the original B judgments.

These views never infer a biological judgment or change canonical responses.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any

from candidate_qtl.route_b_semantic_bundle import INITIAL_MATERIAL_PROJECTION, similarity_annotation_context

CONTRACT = "route-b-initial-material-1.0.0"
CONTEXT_FIELDS = (
    "reference_function", "function_description", "reference_mutant_phenotype", "mutant_phenotype",
    "reference_pathway", "pathway", "reference_pathway_id", "pathway_id", "go_terms", "keywords",
    "tissue_specificity", "developmental_stage", "protein_family", "protein_name",
    "reference_domain_description", "domain_description", "domain_architecture", "domain_ids",
    "reference_domain_ids", "local_annotation_description", "local_annotation_go_terms",
)
IDENTITY_FIELDS = (
    "gene_id", "reference_gene", "reference_gene_id", "reference_species", "reference_source_id", "source_id",
    "source_ids", "source_url", "reference_source_url", "literature_source_ids", "phenotype_source_ids",
    "evidence_type", "evidence_codes", "provenance_tier", "provenance_detail", "source_database",
    "source_release", "retrieved_at", "mapping_state", "mapping_relationship", "function_transfer_scope",
)
STATES = {"assessed_claim", "assessed_concrete_lead", "assessed_no_specific_bridge", "assessed_insufficient_material", "not_assessed"}
INTERPRETED = STATES - {"assessed_insufficient_material", "not_assessed"}


def encoded(value: Any) -> bytes:
    return (json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def collect_material(bundle: dict) -> tuple[list[dict], dict[str, str]]:
    """Share only identical full context; keep each source/member binding separate."""
    grouped: dict[str, dict] = {}
    roster: dict[str, str] = {}
    unit_ids = set()
    projection = bundle.get("initial_material_projection")
    if projection not in (None, INITIAL_MATERIAL_PROJECTION):
        raise ValueError("unsupported initial material projection")
    for i, unit in enumerate(bundle["units"]):
        uid, members, genes = unit["unit_id"], unit["member_interval_genes"], unit["genes"]
        if (not isinstance(uid, str) or not uid or uid in unit_ids or not members
                or any(not isinstance(g, str) or not g or g in roster for g in members)
                or len(set(members)) != len(members)
                or sorted(g["gene_id"] for g in genes) != sorted(members)):
            raise ValueError("invalid initial material roster")
        unit_ids.add(uid)
        roster.update({g: uid for g in members})
        records = []
        for j, gene in enumerate(genes):
            pointer = f"/units/{i}/genes/{j}"
            records.append((gene, [gene["gene_id"]], pointer, "gene"))
            if isinstance(gene.get("atlas_material_context"), dict):
                records.append((gene["atlas_material_context"], [gene["gene_id"]], pointer + "/atlas_material_context", "gene"))
            if projection == INITIAL_MATERIAL_PROJECTION:
                for k, hit in enumerate(gene.get("sequence_similarity_hits", [])):
                    records.append((hit, [gene["gene_id"]], pointer + f"/sequence_similarity_hits/{k}", "sequence_similarity_hit"))
        for j, card in enumerate(unit.get("reference_material_cards", [])):
            gid = card.get("interval_gene_id")
            if not gid and card.get("gene_id") in members:
                gid = card["gene_id"]
            if gid and gid not in members:
                raise ValueError("initial material member outside unit")
            records.append((card, [gid] if gid else members, f"/units/{i}/reference_material_cards/{j}", "gene" if gid else "unit_context"))
        for record, bound, pointer, scope in records:
            is_hit = scope == "sequence_similarity_hit"
            context = (similarity_annotation_context(record) if is_hit else
                       {k: record[k] for k in CONTEXT_FIELDS if record.get(k) not in (None, "", [])})
            if not context:
                continue
            mid = "material-" + digest(encoded(context))
            group = grouped.setdefault(mid, {"material_id": mid, "context": context, "bindings": []})
            if is_hit:
                group["material_status"] = "UNVERIFIED_sequence_similarity_annotation_context"
            group["bindings"].append({
                "unit_id": uid, "member_interval_genes": list(bound), "source_pointer": pointer,
                "binding_scope": scope,
                "reference": dict(record) if is_hit else {k: record[k] for k in IDENTITY_FIELDS if record.get(k) not in (None, "", [])},
            })
    return list(grouped.values()), roster


def plan_initial_material(bundle: dict, max_bytes: int, *, required_pointers: set[str] | None = None) -> tuple[dict[str, bytes], dict]:
    groups, roster = collect_material(bundle)
    full_count = len(groups)
    full_covered = {gid for group in groups for binding in group["bindings"] for gid in binding["member_interval_genes"]}
    if required_pointers is not None:
        groups = [{**g, "bindings": [b for b in g["bindings"] if b["source_pointer"] in required_pointers]}
                  for g in groups if any(b["source_pointer"] in required_pointers for b in g["bindings"])]
    files, pages, current = {}, [], []

    def flush():
        data = encoded(current)
        name = f"research/initial-material-{len(pages):04d}.json"
        files[name] = data
        pages.append({"file": name, "sha256": digest(data), "bytes": len(data),
                      "material_ids": [g["material_id"] for g in current],
                      "over_bound_bytes": max(0, len(data) - max_bytes)})

    for group in groups:
        if current and len(encoded(current + [group])) > max_bytes:
            flush()
            current = []
        current.append(group)
    if current:
        flush()
    covered = {g for m in groups for b in m["bindings"] for g in b["member_interval_genes"]}
    return files, {"contract": CONTRACT, "pages": pages, "material_group_count": len(groups),
                   "full_material_group_count": full_count,
                   "delivery_scope": "all_initial_context" if required_pointers is None else "unread_followup_context",
                   "gene_count": len(roster), "genes_with_material": len(full_covered),
                   "genes_without_material": sorted(set(roster) - full_covered),
                   "delivered_gene_count": len(covered),
                   "interpretation_state": "not_assessed", "full_context_fields": list(CONTEXT_FIELDS) + (
                       ["hit_protein_name", "hit_description"] if bundle.get("initial_material_projection") == INITIAL_MATERIAL_PROJECTION else []),
                   "meaning": "Available context, not characterized genes, verified support, or completed reading. Unit-context binding is not gene-specific evidence."}


def expand_assessments(payload: dict, groups: list[dict], roster: dict[str, str]) -> list[dict]:
    """Expand explicitly authored applications of identical context, never infer them."""
    records = payload.get("assessments", [])
    shared = payload.get("shared_assessments", [])
    if "assessments" not in payload and "shared_assessments" not in payload:
        raise ValueError("initial assessments are required")
    if not isinstance(records, list) or not isinstance(shared, list):
        raise ValueError("initial assessments must be lists")
    material = {g["material_id"]: g for g in groups}
    present = {gid for g in groups for b in g["bindings"] for gid in b["member_interval_genes"]}
    result, seen = list(records), set()
    for record in shared:
        if not isinstance(record, dict):
            raise ValueError("invalid shared assessment")
        aid, mid = record.get("assessment_id"), record.get("material_id", "")
        members = record.get("member_gene_ids")
        if not isinstance(aid, str) or not aid.strip() or aid in seen:
            raise ValueError("invalid or duplicate shared assessment ID")
        seen.add(aid)
        if (not isinstance(members, list) or not members
                or any(not isinstance(gid, str) or gid not in roster for gid in members)
                or len(set(members)) != len(members)):
            raise ValueError("shared assessment requires explicit unique interval members")
        if mid:
            if mid not in material:
                raise ValueError("unknown shared material ID")
            bindings = material[mid]["bindings"]
        else:
            if (any(gid in present for gid in members) or record.get("assessed_fields")
                    or record.get("state") not in {"not_assessed", "assessed_insufficient_material"}):
                raise ValueError("absent-material sharing requires actual absent context and unresolved state")
            bindings = []
        pointers_by_gene: dict[str, set[str]] = {}
        for binding in bindings:
            for member in binding["member_interval_genes"]:
                pointers_by_gene.setdefault(member, set()).add(binding["source_pointer"])
        for gid in members:
            pointers = sorted(pointers_by_gene.get(gid, set()))
            if mid and not pointers:
                raise ValueError("shared assessment member does not bind this exact material")
            row = {k: v for k, v in record.items() if k not in {"assessment_id", "material_id", "member_gene_ids"}}
            if "used_material_pointers" in row or "gene_id" in row or "unit_id" in row:
                raise ValueError("shared application bindings must come from the frozen context")
            result.append({**row, "gene_id": gid, "unit_id": roster[gid],
                           "used_material_pointers": pointers, "shared_assessment_id": aid})
    return result


def summarize_initial(bundle_path: str, expected_sha256: str, assessments_path: str, out_dir: str) -> dict:
    source, answers, out = Path(bundle_path), Path(assessments_path), Path(out_dir)
    if out.exists():
        raise ValueError("initial summary output already exists")
    raw, answer_raw = source.read_bytes(), answers.read_bytes()
    if not expected_sha256 or digest(raw) != expected_sha256:
        raise ValueError("initial assessment bundle SHA-256 mismatch")
    bundle, payload = json.loads(raw), json.loads(answer_raw)
    if (bundle.get("stage") != "post_b0_semantic_bundle" or bundle.get("route") != "interval_first"
            or not bundle.get("case_id") or payload.get("bundle_sha256") != expected_sha256):
        raise ValueError("initial assessments must bind the frozen semantic bundle")
    groups, roster = collect_material(bundle)
    available: dict[str, dict[str, set[str]]] = {gid: {} for gid in roster}
    for group in groups:
        for binding in group["bindings"]:
            for gid in binding["member_interval_genes"]:
                available[gid][binding["source_pointer"]] = set(group["context"])
    assessments = expand_assessments(payload, groups, roster)
    by_gene: dict[str, list[dict]] = {gid: [] for gid in roster}
    seen = set()
    for record in assessments:
        if not isinstance(record, dict):
            raise ValueError("invalid initial assessment record")
        gid, uid, state = record.get("gene_id"), record.get("unit_id"), record.get("state")
        if gid not in roster or roster[gid] != uid or state not in STATES:
            raise ValueError("invalid initial assessment member or state")
        pointers, fields = record.get("used_material_pointers", []), record.get("assessed_fields", [])
        if (not isinstance(pointers, list) or not all(isinstance(p, str) for p in pointers)
                or len(set(pointers)) != len(pointers) or any(p not in available[gid] for p in pointers)
                or not isinstance(fields, list) or not all(isinstance(f, str) for f in fields)
                or len(set(fields)) != len(fields)
                # A field read from any one of the listed materials counts; it need not be in all.
                or any(not any(f in available[gid][p] for p in pointers) for f in fields)
                or (fields and not pointers)):
            raise ValueError("initial assessment material pointer/field mismatch")
        if not isinstance(record.get("reason"), str) or not record["reason"].strip():
            raise ValueError("initial assessment requires an original reason")
        key = encoded({**record, "used_material_pointers": sorted(pointers), "assessed_fields": sorted(fields)})
        if key in seen:
            raise ValueError("duplicate initial assessment")
        seen.add(key)
        by_gene[gid].append(record)
    rows = []
    for gid, uid in roster.items():
        records = by_gene[gid]
        counts = Counter(r["state"] for r in records)
        interpreted = sum(counts[s] for s in INTERPRETED)
        unresolved = counts["assessed_insufficient_material"] + counts["not_assessed"]
        used: dict[str, set[str]] = {}
        for record in records:
            if record["state"] != "not_assessed":
                for pointer in record.get("used_material_pointers", []):
                    used.setdefault(pointer, set()).update(record.get("assessed_fields", []))
        unread = {p: sorted(fs - used.get(p, set())) for p, fs in available[gid].items() if fs - used.get(p, set())}
        state = ("partially_interpreted" if interpreted and unresolved else "interpreted" if interpreted
                 else "insufficient_only" if counts["assessed_insufficient_material"] else "not_assessed")
        rows.append({"gene_id": gid, "unit_id": uid,
                     "material_availability": "present" if available[gid] else "absent",
                     "interpretation_state": state, "original_state_counts": dict(counts),
                     "has_original_claim": bool(counts["assessed_claim"]),
                     "has_concrete_lead": bool(counts["assessed_concrete_lead"]),
                     "has_no_bridge_judgment": bool(counts["assessed_no_specific_bridge"]),
                     "has_unresolved_judgment": bool(unresolved),
                     "judgments_without_material_pointer": sum(
                         r["state"] != "not_assessed" and not r.get("used_material_pointers") for r in records),
                     "used_material_pointers": sorted(used), "fields_without_read_receipt": unread,
                     "original_assessments": records})
    summary = {"contract": CONTRACT, "case_id": bundle["case_id"], "gene_count": len(rows),
               "input_bundle": {"path": str(source.resolve()), "sha256": expected_sha256},
               "input_assessments": {"path": str(answers.resolve()), "sha256": digest(answer_raw)},
               "assessment_encoding": {"shared_count": len(payload.get("shared_assessments", [])),
                                       "expanded_count": len(assessments), "original_bytes": len(answer_raw),
                                       "expanded_bytes": len(encoded(assessments))},
               "material_counts": dict(Counter(r["material_availability"] for r in rows)),
               "interpretation_counts": dict(Counter(r["interpretation_state"] for r in rows)),
               "genes_with_fields_without_read_receipt": sum(bool(r["fields_without_read_receipt"]) for r in rows),
               "genes_with_judgments_without_material_pointer": sum(bool(r["judgments_without_material_pointer"]) for r in rows),
               "limitations": "Counts of original judgments, not new biology or reading verification. Missing field receipts do not require exhaustive reading/search or veto candidates. Canonical responses, Tier and Wet are unchanged."}
    out.mkdir(parents=True)
    data = b"".join(encoded(r) for r in rows)
    (out / "initial_assessment_genes.jsonl").write_bytes(data)
    summary["genes"] = {"file": "initial_assessment_genes.jsonl", "sha256": digest(data)}
    (out / "initial_assessment_summary.json").write_bytes(encoded(summary))
    return summary
