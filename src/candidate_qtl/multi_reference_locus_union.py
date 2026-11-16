"""Multi-reference candidate universe and locus-entity union.

This is a reporting/identity layer.  It does not make an assembly canonical:
the primary table is an explicit run choice, while alternate-only and
ambiguous copies remain visible.  Only a verified 1:1 bridge may collapse two
native IDs into one locus entity.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

from .assembly_bridge import _adapter_eligible, _read
from .core import file_fingerprint, runtime_path

UNION_SCHEMA_VERSION = "1.0"
CANDIDATE_COLUMNS = [
    "candidate_key", "native_gene_id", "assembly_label", "assembly_role",
    "chrom", "start", "end", "protein_length", "input_status", "special_status", "provenance",
    "candidate_status", "locus_entity_id", "entity_relation", "recall_universe",
    "conditional_compression_input",
]
COMPRESSION_COLUMNS = [
    "gene_id", "locus_entity_id", "assembly_label", "assembly_role", "chrom",
    "start", "end", "protein_length", "candidate_status", "special_status",
    "provenance",
]
ENTITY_COLUMNS = [
    "locus_entity_id", "representative_candidate_key", "representative_native_gene_id",
    "member_candidate_keys", "member_native_gene_ids", "member_count", "assembly_labels",
    "entity_status", "recall_universe", "conditional_compression_input",
]
LEDGER_COLUMNS = [
    "source_full_id", "target_full_id", "source_candidate_key", "target_candidate_key",
    "mapping_state", "relationship", "source_family_id", "target_family_id",
    "synteny_state", "link_state", "collapse_allowed", "function_transfer_allowed",
    "function_transfer_scope", "reason",
]


def _safe(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "unknown"


def _first_column(frame: pd.DataFrame, names: tuple[str, ...], label: str) -> str:
    for name in names:
        if name in frame.columns:
            return name
    raise ValueError(f"{label} requires one of: {', '.join(names)}")


def _candidate_table(path_arg: str | Path, label: str, role: str) -> list[dict[str, str]]:
    frame = _read(path_arg, f"{role} candidate table")
    id_column = _first_column(frame, ("native_gene_id", "gene_id", "full_id"), f"{role} candidate table")
    assembly_column = next((c for c in ("assembly", "assembly_id", "native_gene_id_namespace", "namespace") if c in frame.columns), "")
    status_column = next((c for c in ("candidate_status", "status", "mapping_status") if c in frame.columns), "")
    provenance_column = next((c for c in ("provenance", "source", "source_id") if c in frame.columns), "")
    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    for row in frame.to_dict("records"):
        gene = str(row.get(id_column, "")).strip()
        if not gene:
            continue
        if gene in seen:
            raise ValueError(f"{role} candidate table contains duplicate native_gene_id: {gene}")
        seen.add(gene)
        status = str(row.get(status_column, "")).strip() if status_column else ""
        special = ";".join(sorted({token for token in ("annotation_rescue", "PAV") if token.casefold() in status.casefold()}))
        rows.append({
            "candidate_key": f"{label}:{gene}", "native_gene_id": gene,
            "assembly_label": str(row.get(assembly_column, "")).strip() if assembly_column else label,
            "assembly_role": role, "chrom": str(row.get("chrom", row.get("contig", ""))).strip(),
            "start": str(row.get("start", row.get("start_bp", ""))).strip(),
            "end": str(row.get("end", row.get("end_bp", ""))).strip(),
            "protein_length": str(row.get("protein_length", "")).strip(),
            "input_status": status, "special_status": special,
            "provenance": str(row.get(provenance_column, "")).strip() if provenance_column else f"candidate_table:{label}",
        })
    if not rows:
        raise ValueError(f"{role} candidate table has no usable native IDs: {path_arg}")
    return rows


def _parse_alternate(value: str) -> tuple[str, str]:
    if "=" not in value:
        raise ValueError("alternate candidate tables must use LABEL=PATH")
    label, path = value.split("=", 1)
    label, path = label.strip(), path.strip()
    if not label or not path:
        raise ValueError("alternate candidate tables must use non-empty LABEL=PATH")
    return label, path


class _DisjointSet:
    def __init__(self, keys: list[str]) -> None:
        self.parent = {key: key for key in keys}

    def find(self, key: str) -> str:
        while self.parent[key] != key:
            self.parent[key] = self.parent[self.parent[key]]
            key = self.parent[key]
        return key

    def union(self, left: str, right: str) -> None:
        a, b = self.find(left), self.find(right)
        if a != b:
            self.parent[b] = a


def build_multi_reference_locus_union(
    primary_candidates_arg: str | Path,
    alternate_candidates_args: list[str],
    mapping_arg: str | Path,
    out_dir_arg: str | Path,
    *,
    primary_label: str = "primary",
    case_id: str = "",
) -> dict[str, object]:
    """Build an explicit candidate universe and conservative locus entities."""
    primary_rows = _candidate_table(primary_candidates_arg, primary_label, "primary")
    rows = primary_rows[:]
    alternate_labels: list[str] = []
    for value in alternate_candidates_args:
        label, path = _parse_alternate(value)
        if label == primary_label or label in alternate_labels:
            raise ValueError(f"duplicate assembly label: {label}")
        alternate_labels.append(label)
        rows.extend(_candidate_table(path, label, "alternate"))
    by_id = {row["native_gene_id"]: row for row in rows}
    by_key = {row["candidate_key"]: row for row in rows}
    if len(by_id) != len(rows):
        # A full native ID can occur in multiple assembly tables only if the
        # caller has mislabeled the inputs; never silently merge such rows.
        duplicates = sorted({row["native_gene_id"] for row in rows if sum(other["native_gene_id"] == row["native_gene_id"] for other in rows) > 1})
        raise ValueError(f"native IDs occur in multiple candidate tables: {', '.join(duplicates)}")

    mapping = _read(mapping_arg, "assembly mapping")
    required = {"source_full_id", "target_full_id", "mapping_state", "relationship"}
    missing = sorted(required - set(mapping.columns))
    if missing:
        raise ValueError(f"assembly mapping is missing columns: {', '.join(missing)}")
    dsu = _DisjointSet(list(by_key))
    link_rows: list[dict[str, str]] = []
    ambiguous_keys: set[str] = set()
    syntenic_keys: set[str] = set()
    for record in mapping.to_dict("records"):
        source = str(record.get("source_full_id", "")).strip()
        target = str(record.get("target_full_id", "")).strip()
        source_row, target_row = by_id.get(source), by_id.get(target)
        state, relationship = str(record.get("mapping_state", "")), str(record.get("relationship", ""))
        adapter_ok, adapter_reason = _adapter_eligible(record)
        if source_row and target_row:
            source_key, target_key = source_row["candidate_key"], target_row["candidate_key"]
            if state == "verified" and relationship == "one_to_one":
                dsu.union(source_key, target_key)
                link_state, collapse = "verified_equivalent", "true"
            elif adapter_ok and relationship == "one_to_many":
                ambiguous_keys.update((source_key, target_key))
                link_state, collapse = "ambiguous_copy_set", "false"
            elif str(record.get("synteny_state", "")).casefold() in {"verified", "concordant", "consistent"}:
                syntenic_keys.update((source_key, target_key))
                link_state, collapse = "alternate_only_syntenic", "false"
            else:
                link_state, collapse = "not_collapsible", "false"
            reason = adapter_reason if not adapter_ok else str(record.get("mapping_reason", ""))
        else:
            source_key = source_row["candidate_key"] if source_row else ""
            target_key = target_row["candidate_key"] if target_row else ""
            link_state, collapse, reason = "endpoint_not_in_candidate_universe", "false", "mapping_endpoint_not_in_candidate_table"
        link_rows.append({
            "source_full_id": source, "target_full_id": target,
            "source_candidate_key": source_key, "target_candidate_key": target_key,
            "mapping_state": state, "relationship": relationship,
            "source_family_id": str(record.get("source_family_id", "")),
            "target_family_id": str(record.get("target_family_id", "")),
            "synteny_state": str(record.get("synteny_state", "")), "link_state": link_state,
            "collapse_allowed": collapse,
            "function_transfer_allowed": "true" if (
                adapter_ok and state == "verified" and relationship == "one_to_one" and source_row and target_row
            ) else "false",
            "function_transfer_scope": (
                "exact_gene" if adapter_ok and state == "verified" and relationship == "one_to_one"
                else "family_set" if adapter_ok and relationship == "one_to_many"
                else "none"
            ),
            "reason": reason,
        })

    groups: dict[str, list[str]] = {}
    for key in by_key:
        groups.setdefault(dsu.find(key), []).append(key)
    # Stable representative: primary native row first, then lexical key.
    representative: dict[str, str] = {}
    for root, members in groups.items():
        representative[root] = sorted(members, key=lambda key: (0 if by_key[key]["assembly_role"] == "primary" else 1, key))[0]
    entity_ids = {root: f"entity:{_safe(representative[root])}" for root in groups}
    candidate_output: list[dict[str, str]] = []
    for key, row in by_key.items():
        root = dsu.find(key)
        is_ambiguous = key in ambiguous_keys
        if row["assembly_role"] == "primary":
            status = "ambiguous_copy_set" if is_ambiguous else "primary_native"
        elif is_ambiguous:
            status = "ambiguous_copy_set"
        elif key in syntenic_keys:
            status = "alternate_only_syntenic"
        elif len(groups[root]) > 1:
            status = "verified_equivalent"
        else:
            status = "alternate_only"
        candidate_output.append({**row, "candidate_status": status,
                                 "locus_entity_id": entity_ids[root],
                                 "entity_relation": "representative" if key == representative[root] else ("verified_equivalent" if len(groups[root]) > 1 else "retained_copy"),
                                 "recall_universe": "true",
                                 "conditional_compression_input": "true" if key == representative[root] else "false"})

    entity_output: list[dict[str, str]] = []
    for root, members in sorted(groups.items(), key=lambda item: entity_ids[item[0]]):
        member_rows = [by_key[key] for key in sorted(members)]
        ambiguous = any(key in ambiguous_keys for key in members)
        entity_output.append({
            "locus_entity_id": entity_ids[root], "representative_candidate_key": representative[root],
            "representative_native_gene_id": by_key[representative[root]]["native_gene_id"],
            "member_candidate_keys": ";".join(sorted(members)),
            "member_native_gene_ids": ";".join(row["native_gene_id"] for row in member_rows),
            "member_count": str(len(members)), "assembly_labels": ";".join(sorted({row["assembly_label"] for row in member_rows})),
            "entity_status": "ambiguous_copy_set" if ambiguous else ("verified_equivalent" if len(members) > 1 else "single_native"),
            "recall_universe": "true", "conditional_compression_input": "true",
        })

    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    candidate_path = out / "candidate_universe.tsv"
    compression_path = out / "compression_candidates.tsv"
    entity_path = out / "locus_entities.tsv"
    ledger_path = out / "locus_link_ledger.tsv"
    pd.DataFrame(candidate_output, columns=CANDIDATE_COLUMNS).sort_values("candidate_key", kind="stable").to_csv(candidate_path, sep="\t", index=False, lineterminator="\n")
    compression_rows = [
        {
            "gene_id": row["native_gene_id"],
            **{column: row.get(column, "") for column in COMPRESSION_COLUMNS if column != "gene_id"},
        }
        for row in candidate_output
        if row["conditional_compression_input"] == "true"
    ]
    pd.DataFrame(compression_rows, columns=COMPRESSION_COLUMNS).sort_values("gene_id", kind="stable").to_csv(
        compression_path, sep="\t", index=False, lineterminator="\n"
    )
    pd.DataFrame(entity_output, columns=ENTITY_COLUMNS).to_csv(entity_path, sep="\t", index=False, lineterminator="\n")
    pd.DataFrame(link_rows, columns=LEDGER_COLUMNS).to_csv(ledger_path, sep="\t", index=False, lineterminator="\n")
    status_counts = pd.Series([row["candidate_status"] for row in candidate_output]).value_counts().sort_index().to_dict()
    manifest = {
        "schema_version": UNION_SCHEMA_VERSION, "command": "multi-reference-locus-union", "case_id": case_id,
        "primary_label": primary_label, "alternate_labels": alternate_labels,
        "inputs": {"primary_candidates": file_fingerprint(runtime_path(primary_candidates_arg), full_hash=True),
                    "mapping": file_fingerprint(runtime_path(mapping_arg), full_hash=True),
                    "alternate_candidates": {label: file_fingerprint(runtime_path(path), full_hash=True) for label, path in (_parse_alternate(value) for value in alternate_candidates_args)}},
        "outputs": {"candidate_universe": str(candidate_path), "compression_candidates": str(compression_path),
                    "locus_entities": str(entity_path), "link_ledger": str(ledger_path)},
        "counts": {"candidate_universe_count": len(candidate_output), "recall_universe_count": len(candidate_output),
                   "conditional_compression_input_count": sum(row["conditional_compression_input"] == "true" for row in candidate_output),
                   "locus_entity_count": len(entity_output), "link_count": len(link_rows),
                   "verified_equivalent_link_count": sum(row["link_state"] == "verified_equivalent" for row in link_rows),
                   "ambiguous_copy_set_link_count": sum(row["link_state"] == "ambiguous_copy_set" for row in link_rows),
                   "candidate_status_counts": status_counts},
        "policy": [
            "primary is an explicit run label, not a universal biological authority",
            "source-only and alternate-only native IDs are retained in the recall universe",
            "only verified one-to-one mappings collapse entity duplicates",
            "one-to-many mappings retain every copy as an ambiguous copy set",
            "provisional/family-mismatch/synteny-conflict/not-evaluated edges never collapse or transfer function",
            "candidate-universe recall and conditional-compression input are separate counts",
            "compression_candidates is the Route-B-facing entity representative table; verified one-to-one duplicates alone are removed",
        ],
        "fingerprints": {
            "candidate_universe": file_fingerprint(candidate_path, full_hash=True),
            "compression_candidates": file_fingerprint(compression_path, full_hash=True),
            "locus_entities": file_fingerprint(entity_path, full_hash=True),
            "link_ledger": file_fingerprint(ledger_path, full_hash=True),
        },
    }
    manifest_path = out / "multi_reference_locus_union_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"candidate_universe": str(candidate_path), "compression_candidates": str(compression_path),
            "locus_entities": str(entity_path),
            "link_ledger": str(ledger_path), "manifest": str(manifest_path), "counts": manifest["counts"]}
