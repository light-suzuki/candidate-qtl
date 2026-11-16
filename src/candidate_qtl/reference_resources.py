"""Accession-aware adapters for frozen legume reference resources.

This module normalizes upstream LIS Gene Function Registry cards and LIS gene
family assignment (``*.gfa.tsv``) files without erasing their native identifier
namespace.  The registry can also describe external native namespaces (for
example LotusBase Gifu) without pretending those identifiers are LIS IDs.

It is intentionally a resource/provenance layer: it does not rank candidates,
infer cross-accession identity, or choose a causal copy.
"""

from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Sequence

import pandas as pd
import yaml

from . import __version__
from .core import file_fingerprint, open_text, runtime_path

SCHEMA_VERSION = "1.1"
RESOURCE_ROLES = {"direct_target", "near_reference", "mutant_reference", "broad_reference"}
NAMESPACE_SCHEMES = {"lis_full_id", "external"}
RUNTIME_FILE_STATUSES = {"available", "unavailable", "not_evaluated"}
DATASET_KEY4_PATTERN = re.compile(r"^[A-Za-z0-9]{4}$")
SNAPSHOT_FREEZER_VERSION = "reference-resource-snapshot-freezer-1.1"
RESOURCE_REQUIRED_FIELDS = (
    "resource_id",
    "species",
    "accession_or_ecotype",
    "assembly_id",
    "annotation_id",
    "native_gene_id_namespace",
    "namespace_scheme",
    "role",
    "source_snapshot_date",
    "source_url_or_upstream_id",
    "license_or_redistribution_note",
)

FUNCTION_COLUMNS = [
    "source_record_id",
    "native_gene_id",
    "native_species",
    "native_species_code",
    "native_accession",
    "native_assembly",
    "native_annotation",
    "native_gene_id_namespace",
    "resource_match_state",
    "gene_symbol",
    "gene_model_pub_name",
    "classical_locus",
    "phenotype_or_function",
    "ontology_ids",
    "source_ids",
    "source_confidence",
    "evidence_type",
    "mutant_or_perturbation_background",
]

FAMILY_COLUMNS = [
    "native_gene_id",
    "native_species_code",
    "native_accession",
    "native_assembly",
    "native_annotation",
    "native_gene_id_namespace",
    "resource_match_state",
    "family_id",
    "protein_id",
    "evalue",
    "score",
    "best_domain_score",
    "pct_hmm_coverage",
]

AUDIT_COLUMNS = ["source_kind", "record_index", "native_gene_id", "status", "reason"]


def parse_lis_full_gene_id(value: object) -> dict[str, str] | None:
    """Parse ``gensp.Accession.gnm#.ann#.gene`` without guessing missing parts."""
    text = str(value or "").strip()
    parts = text.split(".", 4)
    if len(parts) != 5 or not all(parts):
        return None
    species_code, accession, assembly_id, annotation_id, local_gene_id = parts
    if not assembly_id.startswith("gnm") or not annotation_id.startswith("ann"):
        return None
    return {
        "full_id": text,
        "species_code": species_code,
        "accession": accession,
        "assembly_id": assembly_id,
        "annotation_id": annotation_id,
        "local_gene_id": local_gene_id,
        "namespace": ".".join(parts[:4]),
    }


def load_reference_resource_registry(registry_arg: str | Path) -> list[dict[str, str]]:
    """Load and validate the small tracked registry describing frozen resources.

    ``lis_full_id`` namespaces are structurally checked against the declared
    accession/assembly/annotation. ``external`` namespaces are kept opaque:
    their identity is defined by the upstream resource, not by LIS syntax.
    """
    path = runtime_path(registry_arg)
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    resources = payload.get("resources") if isinstance(payload, dict) else payload
    if not isinstance(resources, list):
        raise ValueError("reference resource registry must be a list or contain a resources list")
    normalized: list[dict[str, str]] = []
    seen_ids: set[str] = set()
    seen_namespaces: set[tuple[str, str]] = set()
    seen_dataset_key4: set[str] = set()
    for index, item in enumerate(resources, 1):
        if not isinstance(item, dict):
            raise ValueError(f"resource registry entry {index} must be a mapping")
        missing = [field for field in RESOURCE_REQUIRED_FIELDS if not str(item.get(field, "")).strip()]
        if missing:
            raise ValueError(f"resource registry entry {index} missing required fields: {', '.join(missing)}")
        row = {str(key): str(value).strip() for key, value in item.items() if value is not None}
        resource_id = row["resource_id"]
        namespace = row["native_gene_id_namespace"]
        scheme = row["namespace_scheme"]
        namespace_key = (scheme, namespace)
        if resource_id in seen_ids:
            raise ValueError(f"duplicate resource_id: {resource_id}")
        if namespace_key in seen_namespaces:
            raise ValueError(f"duplicate native namespace for scheme {scheme}: {namespace}")
        if row["role"] not in RESOURCE_ROLES:
            raise ValueError(f"invalid resource role {row['role']!r}; expected one of {sorted(RESOURCE_ROLES)}")
        if scheme not in NAMESPACE_SCHEMES:
            raise ValueError(f"invalid namespace_scheme {scheme!r}; expected one of {sorted(NAMESPACE_SCHEMES)}")
        if scheme == "lis_full_id":
            parsed_namespace = parse_lis_full_gene_id(namespace + ".__placeholder__")
            if parsed_namespace is None or parsed_namespace["namespace"] != namespace:
                raise ValueError(f"invalid LIS native_gene_id_namespace: {namespace}")
            if parsed_namespace["accession"] != row["accession_or_ecotype"]:
                raise ValueError(f"registry accession does not match namespace for {resource_id}")
            if parsed_namespace["assembly_id"] != row["assembly_id"] or parsed_namespace["annotation_id"] != row["annotation_id"]:
                raise ValueError(f"registry assembly/annotation does not match namespace for {resource_id}")
        dataset_key4 = row.get("dataset_key4", "").strip()
        if dataset_key4:
            if DATASET_KEY4_PATTERN.fullmatch(dataset_key4) is None:
                raise ValueError(f"invalid dataset_key4 for {resource_id}: {dataset_key4!r}")
            if dataset_key4 in seen_dataset_key4:
                raise ValueError(f"duplicate dataset_key4: {dataset_key4}")
            expected_resource_id = f"{namespace}.{dataset_key4}"
            if resource_id != expected_resource_id:
                raise ValueError(f"resource_id does not match dataset_key4 for {resource_id}; expected {expected_resource_id}")
            dataset_name = resource_id.split(".", 1)[1]
            required_text = row.get("required", "").casefold()
            if required_text not in {"true", "false"}:
                raise ValueError(f"dataset resource {resource_id} requires explicit boolean required: true|false")
            if row["role"] == "direct_target" and required_text != "true":
                raise ValueError(f"direct_target dataset resource must be required: {resource_id}")
            expected_runtime_names = {
                "manifest_path": f"MANIFEST.{dataset_name}.yml",
                "runtime_manifest_path": f"MANIFEST.{dataset_name}.yml",
                "checksum_path": f"CHECKSUM.{dataset_name}.md5",
                "runtime_checksum_path": f"CHECKSUM.{dataset_name}.md5",
            }
            for field, expected_name in expected_runtime_names.items():
                value = row.get(field, "")
                if value:
                    parts = value.replace("\\", "/").rstrip("/").split("/")
                    if parts[-1] != expected_name or len(parts) < 2 or parts[-2] != dataset_name:
                        raise ValueError(
                            f"{field} does not match dataset_key4 for {resource_id}; expected "
                            f"{dataset_name}/{expected_name}"
                        )
            seen_dataset_key4.add(dataset_key4)
        seen_ids.add(resource_id)
        seen_namespaces.add(namespace_key)
        normalized.append(row)
    return normalized


def _runtime_path_value(row: Mapping[str, str], keys: Sequence[str]) -> str:
    for key in keys:
        value = str(row.get(key, "")).strip()
        if value:
            return value
    return ""


def _declared_runtime_payloads(row: Mapping[str, str]) -> dict[str, str]:
    """Return each declared runtime payload once, preserving future ``*_path`` fields."""
    manifest_aliases = ("manifest_path", "runtime_manifest_path", "manifest", "MANIFEST")
    checksum_aliases = ("checksum_path", "runtime_checksum_path", "checksum", "CHECKSUM")
    payloads: dict[str, str] = {}
    manifest_value = _runtime_path_value(row, manifest_aliases)
    checksum_value = _runtime_path_value(row, checksum_aliases)
    if manifest_value:
        payloads["manifest"] = manifest_value
    if checksum_value:
        payloads["checksum"] = checksum_value

    alias_path_fields = {"manifest_path", "runtime_manifest_path", "checksum_path", "runtime_checksum_path"}
    for field in sorted(row):
        if not field.endswith("_path") or field in alias_path_fields:
            continue
        value = str(row.get(field, "")).strip()
        if value:
            payloads[field.removesuffix("_path")] = value
    return payloads


def _resolve_runtime_file(value: str, runtime_root: Path | None) -> tuple[Path | None, str]:
    """Resolve a declared local runtime payload path without fetching remote data."""
    if not value:
        return None, "path_not_declared"
    if "://" in value:
        return None, "remote_source_not_fetched"
    path = Path(value)
    if not path.is_absolute():
        if runtime_root is None:
            return None, "runtime_root_required_for_relative_path"
        path = runtime_root / path
    return path, ""


def _freeze_runtime_file(value: str, runtime_root: Path | None) -> dict[str, object]:
    path, reason = _resolve_runtime_file(value, runtime_root)
    if path is None:
        return {"status": "not_evaluated", "reason": reason, "path": value or None, "sha256": None}
    if not path.is_file():
        return {"status": "unavailable", "reason": "path_missing", "path": str(path), "sha256": None}
    fingerprint = file_fingerprint(path, full_hash=True)
    return {
        "status": "available",
        "reason": "frozen_full_sha256",
        "path": str(path),
        "size_bytes": fingerprint["size_bytes"],
        "sha256": fingerprint["sha256"],
        "method": "full_sha256",
    }


def freeze_reference_resource_snapshots(
    registry_arg: str | Path,
    out_dir_arg: str | Path,
    *,
    runtime_root_arg: str | Path | None = None,
) -> dict[str, object]:
    """Freeze every declared per-dataset runtime payload with full SHA-256.

    Legacy ``availability`` remains MANIFEST/CHECKSUM-based for compatibility.
    ``snapshot_completeness`` separately reports whether all declared payloads
    were frozen. Missing optional files stay explicit ``unavailable`` or
    ``not_evaluated`` records and never become negative biological evidence.
    This function never joins native IDs, symbols, accessions, or candidates.
    """
    registry = load_reference_resource_registry(registry_arg)
    registry_path = runtime_path(registry_arg)
    runtime_root = runtime_path(runtime_root_arg) if runtime_root_arg is not None else None
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    snapshots: list[dict[str, object]] = []
    for row in registry:
        declared_payloads = _declared_runtime_payloads(row)
        payloads = {
            name: _freeze_runtime_file(value, runtime_root)
            for name, value in declared_payloads.items()
        }
        manifest = payloads.get("manifest", _freeze_runtime_file("", runtime_root))
        checksum = payloads.get("checksum", _freeze_runtime_file("", runtime_root))
        file_statuses = {str(manifest["status"]), str(checksum["status"])}
        if file_statuses == {"available"}:
            availability = "available"
        elif "not_evaluated" in file_statuses:
            availability = "not_evaluated"
        else:
            availability = "unavailable"

        payload_statuses = [str(payload["status"]) for payload in payloads.values()]
        if not payload_statuses or all(status == "not_evaluated" for status in payload_statuses):
            snapshot_completeness = "not_evaluated"
        elif all(status == "available" for status in payload_statuses):
            snapshot_completeness = "complete"
        elif any(status == "available" for status in payload_statuses):
            snapshot_completeness = "partial"
        else:
            snapshot_completeness = "unavailable"

        required = str(row.get("required", "")).strip().lower() in {"1", "true", "yes", "required"}
        snapshots.append(
            {
                "resource_id": row["resource_id"],
                "dataset_key4": row.get("dataset_key4", ""),
                "native_gene_id_namespace": row["native_gene_id_namespace"],
                "accession_or_ecotype": row["accession_or_ecotype"],
                "assembly_id": row["assembly_id"],
                "annotation_id": row["annotation_id"],
                "role": row["role"],
                "source_snapshot_date": row["source_snapshot_date"],
                "source_url_or_upstream_id": row["source_url_or_upstream_id"],
                "required": required,
                "availability": availability,
                "snapshot_completeness": snapshot_completeness,
                "declared_payload_count": len(payloads),
                "frozen_payload_count": sum(payload["status"] == "available" for payload in payloads.values()),
                "payloads": payloads,
                "manifest": manifest,
                "checksum": checksum,
            }
        )
    snapshot_path = out_dir / "reference_resource_snapshot_manifest.json"
    summary = {
        "schema_version": SCHEMA_VERSION,
        "registry_schema_version": SCHEMA_VERSION,
        "parser_version": SNAPSHOT_FREEZER_VERSION,
        "candidate_qtl_version": __version__,
        "producer_class": "reference_resource_runtime_snapshot_freezer",
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "registry": file_fingerprint(registry_path, full_hash=True),
        "rank_contribution": "none",
        "cross_accession_identity_inferred": False,
        "resource_count": len(snapshots),
        "counts": {
            "available": sum(item["availability"] == "available" for item in snapshots),
            "unavailable": sum(item["availability"] == "unavailable" for item in snapshots),
            "not_evaluated": sum(item["availability"] == "not_evaluated" for item in snapshots),
            "required_blocked": sum(item["required"] and item["availability"] != "available" for item in snapshots),
            "optional_unavailable_or_not_evaluated": sum(not item["required"] and item["availability"] != "available" for item in snapshots),
            "full_snapshot_complete": sum(item["snapshot_completeness"] == "complete" for item in snapshots),
            "full_snapshot_partial": sum(item["snapshot_completeness"] == "partial" for item in snapshots),
            "full_snapshot_unavailable": sum(item["snapshot_completeness"] == "unavailable" for item in snapshots),
            "full_snapshot_not_evaluated": sum(item["snapshot_completeness"] == "not_evaluated" for item in snapshots),
            "declared_payloads": sum(int(item["declared_payload_count"]) for item in snapshots),
            "frozen_payloads": sum(int(item["frozen_payload_count"]) for item in snapshots),
        },
        "resources": snapshots,
        "policy": [
            "every declared local resource *_path payload is frozen independently with full SHA-256 when available",
            "legacy resource availability remains based on MANIFEST and CHECKSUM for backward compatibility",
            "full snapshot completeness is reported separately from resource availability",
            "optional unavailable resources remain explicit unavailable or not_evaluated and do not fail available resources",
            "native full-ID namespaces and dataset KEY4 identities are preserved",
            "cross-accession identity, symbol mapping, scoring, and candidate-universe filtering are not performed",
        ],
    }
    summary["state"] = (
        "blocked_required_resource_missing"
        if summary["counts"]["required_blocked"]
        else "ready_available_resources_frozen"
        if summary["counts"]["available"]
        else "not_ready_no_available_resources"
    )
    snapshot_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary | {"manifest": str(snapshot_path)}


def _source_record_id(document: dict[str, object]) -> str:
    encoded = json.dumps(document, ensure_ascii=False, sort_keys=True, default=str, separators=(",", ":"))
    return "gfr-" + hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:16]


def _string_list(value: object) -> list[str]:
    if isinstance(value, list):
        return [str(item).strip() for item in value if str(item).strip()]
    text = str(value or "").strip()
    return [text] if text else []


def _ontology_ids(document: dict[str, object]) -> list[str]:
    values: set[str] = set()
    for trait in document.get("traits", []) if isinstance(document.get("traits", []), list) else []:
        if not isinstance(trait, dict):
            continue
        for key in ("entity", "relation", "quality"):
            text = str(trait.get(key, "")).strip()
            if text:
                values.add(text)
    return sorted(values)


def _source_ids(document: dict[str, object]) -> list[str]:
    values: set[str] = set()
    for reference in document.get("references", []) if isinstance(document.get("references", []), list) else []:
        if not isinstance(reference, dict):
            continue
        pmid = str(reference.get("pmid", "")).strip()
        doi = str(reference.get("doi", "")).strip()
        if pmid:
            values.add(f"PMID:{pmid}")
        if doi:
            values.add(f"DOI:{doi}")
    return sorted(values)


def normalize_lis_gene_function_registry(
    traits_arg: str | Path,
    registered_namespaces: set[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Normalize LIS Gene Function Registry YAML documents without cross-ID inference."""
    path = runtime_path(traits_arg)
    rows: list[dict[str, object]] = []
    audit: list[dict[str, object]] = []
    namespaces = registered_namespaces or set()
    with open_text(path) as handle:
        documents = list(yaml.safe_load_all(handle))
    for index, document in enumerate(documents, 1):
        if document is None:
            continue
        if not isinstance(document, dict):
            audit.append({"source_kind": "gene_function", "record_index": index, "native_gene_id": "", "status": "invalid", "reason": "document_not_mapping"})
            continue
        full_id = document.get("gene_model_full_id", document.get("gene_model_full_ID", ""))
        parsed = parse_lis_full_gene_id(full_id)
        if parsed is None:
            audit.append({"source_kind": "gene_function", "record_index": index, "native_gene_id": str(full_id or ""), "status": "unresolved", "reason": "missing_or_invalid_full_lis_gene_id"})
            continue
        symbols = _string_list(document.get("gene_symbols", []))
        namespace = parsed["namespace"]
        row = {
            "source_record_id": _source_record_id(document),
            "native_gene_id": parsed["full_id"],
            "native_species": str(document.get("scientific_name", "")).strip(),
            "native_species_code": parsed["species_code"],
            "native_accession": parsed["accession"],
            "native_assembly": parsed["assembly_id"],
            "native_annotation": parsed["annotation_id"],
            "native_gene_id_namespace": namespace,
            "resource_match_state": "registered_namespace" if namespace in namespaces else "unregistered_namespace",
            "gene_symbol": symbols[0] if symbols else "",
            "gene_model_pub_name": str(document.get("gene_model_pub_name", document.get("gene_model_pub_names", ""))).strip(),
            "classical_locus": str(document.get("classical_locus", "")).strip(),
            "phenotype_or_function": str(document.get("phenotype_synopsis", "")).strip(),
            "ontology_ids": ";".join(_ontology_ids(document)),
            "source_ids": ";".join(_source_ids(document)),
            "source_confidence": str(document.get("confidence", "")).strip(),
            "evidence_type": "lis_gene_function_registry",
            "mutant_or_perturbation_background": str(document.get("mutant_background", "")).strip(),
        }
        rows.append(row)
        audit.append({"source_kind": "gene_function", "record_index": index, "native_gene_id": parsed["full_id"], "status": "normalized", "reason": row["resource_match_state"]})
    return pd.DataFrame(rows, columns=FUNCTION_COLUMNS), pd.DataFrame(audit, columns=AUDIT_COLUMNS)


def normalize_lis_gene_family_assignments(
    family_arg: str | Path,
    registered_namespaces: set[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Normalize LIS ``*.gfa.tsv`` family assignments while preserving quality fields."""
    path = runtime_path(family_arg)
    namespaces = registered_namespaces or set()
    rows: list[dict[str, object]] = []
    audit: list[dict[str, object]] = []
    with open_text(path) as handle:
        raw_lines = [line.rstrip("\r\n") for line in handle if line.strip()]
    if not raw_lines:
        return pd.DataFrame(columns=FAMILY_COLUMNS), pd.DataFrame(columns=AUDIT_COLUMNS)
    header = raw_lines[0].lstrip("#").split("\t")
    if len(header) == 1:
        header = raw_lines[0].lstrip("#").split()
    required = ["gene", "family", "protein", "e-value", "score"]
    missing = [column for column in required if column not in header]
    if missing:
        raise ValueError(f"LIS gene-family assignment missing required columns: {', '.join(missing)}")
    for index, raw in enumerate(raw_lines[1:], 1):
        if raw.startswith("#"):
            continue
        fields = raw.split("\t")
        if len(fields) != len(header):
            audit.append({"source_kind": "family_assignment", "record_index": index, "native_gene_id": fields[0] if fields else "", "status": "invalid", "reason": "column_count_mismatch"})
            continue
        record = dict(zip(header, fields, strict=True))
        parsed = parse_lis_full_gene_id(record.get("gene", ""))
        if parsed is None:
            audit.append({"source_kind": "family_assignment", "record_index": index, "native_gene_id": record.get("gene", ""), "status": "unresolved", "reason": "missing_or_invalid_full_lis_gene_id"})
            continue
        namespace = parsed["namespace"]
        rows.append({
            "native_gene_id": parsed["full_id"],
            "native_species_code": parsed["species_code"],
            "native_accession": parsed["accession"],
            "native_assembly": parsed["assembly_id"],
            "native_annotation": parsed["annotation_id"],
            "native_gene_id_namespace": namespace,
            "resource_match_state": "registered_namespace" if namespace in namespaces else "unregistered_namespace",
            "family_id": record.get("family", ""),
            "protein_id": record.get("protein", ""),
            "evalue": record.get("e-value", ""),
            "score": record.get("score", ""),
            "best_domain_score": record.get("best_domain_score", ""),
            "pct_hmm_coverage": record.get("pct_hmm_coverage", ""),
        })
        audit.append({"source_kind": "family_assignment", "record_index": index, "native_gene_id": parsed["full_id"], "status": "normalized", "reason": "registered_namespace" if namespace in namespaces else "unregistered_namespace"})
    return pd.DataFrame(rows, columns=FAMILY_COLUMNS), pd.DataFrame(audit, columns=AUDIT_COLUMNS)


def build_legume_reference_snapshot(
    registry_arg: str,
    traits_arg: str,
    family_arg: str,
    out_dir_arg: str,
    *,
    full_hash: bool = False,
    runtime_root_arg: str | Path | None = None,
) -> dict[str, object]:
    """Build deterministic normalized resource artifacts; no candidate scoring occurs."""
    registry_path = runtime_path(registry_arg)
    traits_path = runtime_path(traits_arg)
    family_path = runtime_path(family_arg)
    registry = load_reference_resource_registry(registry_path)
    all_namespaces = {row["native_gene_id_namespace"] for row in registry}
    lis_namespaces = {
        row["native_gene_id_namespace"]
        for row in registry
        if row["namespace_scheme"] == "lis_full_id"
    }
    functions, function_audit = normalize_lis_gene_function_registry(traits_path, lis_namespaces)
    families, family_audit = normalize_lis_gene_family_assignments(family_path, lis_namespaces)

    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    function_path = out_dir / "normalized_functional_evidence.tsv"
    family_path_out = out_dir / "family_assignment.tsv"
    audit_path = out_dir / "resource_ingest_audit.tsv"
    manifest_path = out_dir / "reference_resources_manifest.json"
    functions.to_csv(function_path, sep="\t", index=False, lineterminator="\n")
    families.to_csv(family_path_out, sep="\t", index=False, lineterminator="\n")
    audit = pd.concat([function_audit, family_audit], ignore_index=True)
    audit.to_csv(audit_path, sep="\t", index=False, lineterminator="\n")

    runtime_snapshot: dict[str, object] | None = None
    if runtime_root_arg is not None or any(
        _runtime_path_value(row, ("manifest_path", "runtime_manifest_path", "manifest", "MANIFEST", "checksum_path", "runtime_checksum_path", "checksum", "CHECKSUM"))
        for row in registry
    ):
        runtime_snapshot = freeze_reference_resource_snapshots(registry_path, out_dir, runtime_root_arg=runtime_root_arg)

    family_ids = set(families["native_gene_id"].astype(str)) if not families.empty else set()
    function_ids = set(functions["native_gene_id"].astype(str)) if not functions.empty else set()
    duplicate_function_records = int(functions["native_gene_id"].duplicated(keep=False).sum()) if not functions.empty else 0
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "producer_class": "legume_reference_resource_normalizer",
        "rank_contribution": "none",
        "cross_accession_identity_inferred": False,
        "resource_count": len(registry),
        "resources": registry,
        "registered_namespaces": sorted(all_namespaces),
        "registered_lis_namespaces": sorted(lis_namespaces),
        "counts": {
            "functional_evidence_records": int(len(functions)),
            "functional_duplicate_native_gene_record_count": duplicate_function_records,
            "family_assignment_records": int(len(families)),
            "functional_registered_namespace_records": int((functions["resource_match_state"] == "registered_namespace").sum()) if not functions.empty else 0,
            "functional_unregistered_namespace_records": int((functions["resource_match_state"] == "unregistered_namespace").sum()) if not functions.empty else 0,
            "family_registered_namespace_records": int((families["resource_match_state"] == "registered_namespace").sum()) if not families.empty else 0,
            "family_unregistered_namespace_records": int((families["resource_match_state"] == "unregistered_namespace").sum()) if not families.empty else 0,
            "function_genes_with_family_assignment": len(function_ids & family_ids),
            "audit_unresolved_or_invalid_records": int(audit["status"].isin(["unresolved", "invalid"]).sum()) if not audit.empty else 0,
        },
        "inputs": {
            "registry": file_fingerprint(registry_path, full_hash),
            "gene_function_registry": file_fingerprint(traits_path, full_hash),
            "gene_family_assignments": file_fingerprint(family_path, full_hash),
        },
        "outputs": {
            "functional_evidence": str(function_path),
            "family_assignments": str(family_path_out),
            "audit": str(audit_path),
            "manifest": str(manifest_path),
        },
        "policy": [
            "native species/accession/assembly/annotation namespaces are preserved",
            "LIS namespaces are structurally validated while external native namespaces remain explicit opaque identifiers",
            "external namespaces are never treated as LIS Gene Function Registry or GFA join keys without an explicit later mapping adapter",
            "bare gene symbols are never used to merge identities",
            "gene-function confidence is preserved as source metadata and is not converted to a ranking weight",
            "duplicate native-gene function cards remain separate source records and are counted, not silently collapsed",
            "family assignment quality fields are preserved without inventing thresholds",
            "cross-accession or cross-species identity is not inferred in this phase",
        ],
    }
    if runtime_snapshot is not None:
        manifest["runtime_snapshot"] = {
            "manifest": runtime_snapshot["manifest"],
            "counts": runtime_snapshot["counts"],
        }
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest | {"manifest": str(manifest_path)}
