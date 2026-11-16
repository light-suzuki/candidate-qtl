"""Comparative reference-resource panel and local-only fallback.

This panel is the resource-contract counterpart to the LLM-proposed known-gene
panel (``reference_gene_panel``): it declares the FASTA / GFF3 / proteome / DB /
API resources used for comparative genomics, grouped by generic role, and audits
their availability.  The native role alone is required; every other resource is
optional and, when absent, is recorded with a machine-readable skip reason
instead of silently dropped.  Roles, resource types, taxa, and identifiers are
kept generic so no crop name or gene-ID system is hard-coded.

Resource availability is provenance, not comparative evidence: a resource being
present does not establish homology, orthology, or trait relevance.
"""

from __future__ import annotations

import json

import pandas as pd

from . import schemas
from .core import file_fingerprint, runtime_path


PANEL_VERSION = "1.0.0"
INPUT_COLUMNS = ["resource_id", "role", "resource_type", "taxon", "source", "uri_or_path"]
OPTIONAL_COLUMNS = ["assembly_id", "annotation_version", "proteome_version", "db_version", "external_api", "review_status", "notes"]
REMOTE_SCHEMES = ("http://", "https://", "ftp://", "ftps://", "s3://", "gs://")
TRUE_TOKENS = frozenset({"1", "true", "yes", "y", "external", "remote", "api"})


def _is_truthy(value: str) -> bool:
    return str(value).strip().casefold() in TRUE_TOKENS


def _is_remote(uri: str, external_api: str) -> bool:
    return _is_truthy(external_api) or str(uri).strip().casefold().startswith(REMOTE_SCHEMES)


def reference_resource_panel(panel_arg: str, out_arg: str, full_hash: bool = False) -> dict[str, object]:
    """Validate and availability-audit a comparative reference-resource panel."""
    panel_path = runtime_path(panel_arg)
    out = runtime_path(out_arg)
    out.mkdir(parents=True, exist_ok=True)

    panel = pd.read_csv(panel_path, sep="\t", dtype=str).fillna("")
    missing_columns = [column for column in INPUT_COLUMNS if column not in panel.columns]
    if missing_columns:
        raise ValueError(f"resource panel missing columns: {', '.join(missing_columns)}")
    for column in OPTIONAL_COLUMNS:
        if column not in panel.columns:
            panel[column] = ""
    if panel.empty:
        raise ValueError("resource panel is empty")
    if panel["resource_id"].duplicated().any():
        raise ValueError("resource panel resource_id must be unique")
    for column in ("resource_id", "role", "resource_type", "taxon", "source", "uri_or_path"):
        if panel[column].str.strip().eq("").any():
            raise ValueError(f"resource panel {column} is mandatory for every resource")
    bad_roles = sorted(set(panel["role"]) - schemas.REFERENCE_RESOURCE_ROLES)
    if bad_roles:
        raise ValueError(f"unknown resource role(s): {', '.join(bad_roles)}; allowed {sorted(schemas.REFERENCE_RESOURCE_ROLES)}")
    bad_types = sorted(set(panel["resource_type"]) - schemas.REFERENCE_RESOURCE_TYPES)
    if bad_types:
        raise ValueError(f"unknown resource_type(s): {', '.join(bad_types)}; allowed {sorted(schemas.REFERENCE_RESOURCE_TYPES)}")
    if not (panel["role"] == "native").any():
        raise ValueError("resource panel must declare at least one native role resource")

    rows: list[dict[str, object]] = []
    fingerprints: dict[str, object] = {}
    external_destinations: list[str] = []
    for record in panel.to_dict("records"):
        resource_id = str(record["resource_id"]).strip()
        uri = str(record["uri_or_path"]).strip()
        external_api = str(record.get("external_api", "")).strip()
        remote = _is_remote(uri, external_api)
        availability, skip_reason, sha256_value, size_bytes = "", "", "", ""
        if remote:
            availability = "declared_external"
            skip_reason = "external_resource_not_fetched_local_first"
            fingerprints[resource_id] = {"status": "external", "uri": uri}
            external_destinations.append(uri)
        else:
            path = runtime_path(uri)
            if not path.is_absolute():
                path = panel_path.parent / path
            if path.is_file():
                fingerprint = file_fingerprint(path, full_hash=full_hash)
                fingerprints[resource_id] = fingerprint
                availability = "available"
                sha256_value = "" if fingerprint.get("sha256") in (None, "") else str(fingerprint["sha256"])
                size_bytes = str(fingerprint.get("size_bytes", ""))
            else:
                availability = "skipped"
                skip_reason = "local_resource_missing"
                fingerprints[resource_id] = {"status": "missing", "path": str(path)}
        if str(record["role"]) == "native" and availability != "available":
            raise ValueError(f"native resource {resource_id!r} must be a locally available file; got {availability} ({skip_reason or uri})")
        rows.append({
            "resource_id": resource_id, "role": record["role"], "resource_type": record["resource_type"], "taxon": record["taxon"],
            "assembly_id": record.get("assembly_id", ""), "annotation_version": record.get("annotation_version", ""),
            "proteome_version": record.get("proteome_version", ""), "db_version": record.get("db_version", ""),
            "source": record["source"], "uri_or_path": uri, "external_api": "true" if remote and _is_truthy(external_api) else external_api,
            "review_status": record.get("review_status", ""), "availability": availability, "skip_reason": skip_reason,
            "resource_sha256": sha256_value, "resource_size_bytes": size_bytes, "notes": record.get("notes", ""),
        })

    frame = pd.DataFrame(rows, columns=schemas.REFERENCE_RESOURCE_PANEL_COLUMNS).sort_values(["role", "resource_id"], kind="stable")
    schemas.validate(frame.fillna(""), schemas.REFERENCE_RESOURCE_PANEL_SCHEMA)
    panel_out = out / "reference_resource_panel.tsv"
    frame.to_csv(panel_out, sep="\t", index=False, lineterminator="\n")

    availability_counts = {str(key): int(value) for key, value in frame["availability"].value_counts().sort_index().items()}
    role_counts = {str(key): int(value) for key, value in frame["role"].value_counts().sort_index().items()}
    external_used = bool(external_destinations)
    manifest = {
        "schema_version": "1.0", "panel_version": PANEL_VERSION, "panel_fingerprint": file_fingerprint(panel_path, full_hash),
        "resource_count": len(frame), "role_counts": role_counts, "roles_present": sorted(role_counts),
        "availability_counts": availability_counts, "native_available": bool((frame["role"].eq("native") & frame["availability"].eq("available")).any()),
        "optional_roles_absent": sorted(schemas.REFERENCE_RESOURCE_ROLES - set(frame["role"]) - {"native"}),
        "skip_reasons": {str(row["resource_id"]): str(row["skip_reason"]) for row in rows if row["skip_reason"]},
        "external_api_used": external_used, "external_destinations": sorted(set(external_destinations)), "local_only": not external_used,
        "fingerprint_mode": "strict_full_hash" if full_hash else "size_tiered", "resource_fingerprints": fingerprints,
        "scientific_limits": [
            "resource availability is provenance, not homology/orthology/trait evidence",
            "external resources are declared but not fetched by this local-first contract",
            "this resource panel is distinct from the LLM-proposed known-gene panel",
        ],
    }
    manifest_out = out / "reference_resource_panel_manifest.json"
    manifest_out.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"panel": str(panel_out), "manifest": str(manifest_out), "resource_count": len(frame),
            "availability_counts": availability_counts, "native_available": manifest["native_available"], "local_only": manifest["local_only"]}
