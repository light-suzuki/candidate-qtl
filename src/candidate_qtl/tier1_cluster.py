"""Deterministic Tier 1 clustering (issue #375, deterministic budget engine).

Tier 1 consumes the exact inventory row dicts Tier 0 wrote (``TIER0_COLUMNS``
keys, optionally plus the ``cds_length``/``protein_length`` auxiliaries) and
assigns every gene to exactly one cluster by a fixed, exhaustive rule order:

1. ``family``     - remaining genes grouped by family id from the family
                     provider table; a gene with several family ids joins the
                     cluster of its lexicographically smallest family id
                     (priority). A single-member family cluster is still
                     ``priority`` by provider-derived evidence, symmetric with
                     the domain rule below.
2. ``annotation``  - genes sharing a non-empty gene_name with at least one
                     other remaining gene (priority); a unique gene_name never
                     forms a cluster by itself (audit finding G1)
3. ``domain``      - remaining genes, keyed by their lexicographically smallest
                     domain id from ``domain_hits`` (priority)
4. ``orthogroup``  - remaining genes, keyed by their smallest orthogroup id
                     (defer)
5. ``model_state_rescue`` - remaining genes with a missing/broken model,
                     missing annotation, sv/pav/promoter evidence, a provider
                     gap (provider tables provided, no resource status cell is
                     ``present``, and at least one cell is ``missing``), or an
                     unevaluated model; one cluster per gene (rescue). The
                     rescue bucket is first-match-wins:
                     sv_pav_regulatory (sv/pav/promoter evidence) >
                     missing_model (missing/broken model) > weak_annotation
                     (missing annotation) > provider_gap > unclassified
                     (unevaluated model); all five buckets are reachable.
                     ``not_evaluated`` cells never count as negative evidence
                     (audit finding G4).
6. ``singleton``   - every remaining gene alone (defer)

``cluster_id`` is content-addressed (``cls_`` + the first 20 hex chars of the
SHA-256 of the sorted member ids joined with the unit separator ``\\x1f``), so
identical inputs always produce identical output. ``dep_id`` uses the same
formula as the evidence-dependency module (``dep_`` + 20 hex chars of the
SHA-256 of the source id). No wall-clock values, no RNG, no timestamps.

Writes ``tier1_clusters.tsv`` (``TIER1_COLUMNS``) and ``tier1_manifest.json``
into ``out_dir`` and returns the manifest.
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path

import pandas as pd

from .core import runtime_path

TIER1_COLUMNS = (
    "cluster_id",
    "seed_reason",
    "representative_gene_id",
    "representative_reason",
    "member_gene_ids",
    "member_count",
    "rule_key",
    "cluster_key",
    "family_ids",
    "domain_ids",
    "orthogroup_ids",
    "allocation",
    "allocation_reason",
    "rescue_bucket",
    "dependency_cluster_ids",
)

ALLOCATIONS = frozenset({"priority", "defer", "rescue"})
RESCUE_BUCKETS = frozenset({"unclassified", "weak_annotation", "missing_model", "sv_pav_regulatory", "provider_gap"})
RULE_KEYS = frozenset({"family", "annotation", "domain", "orthogroup", "model_state_rescue", "singleton"})

# Model-state priority for representative selection; unknown values rank last.
_MODEL_STATUS_RANK = {"complete_orf": 0, "partial_or_unverified_orf": 1, "no_cds": 2, "not_evaluated": 3}
# Tier 0 evidence-status cell names (note: domains -> domain_status).
_RESOURCE_STATUS_CELLS = (
    "domain_status",
    "orthogroup_status",
    "expression_status",
    "variant_status",
    "sv_status",
    "pav_status",
    "promoter_status",
)
_RESCUE_TRIGGER_MODEL = frozenset({"no_cds", "partial_or_unverified_orf"})
_DEPENDENCY_SOURCE_CELLS = ("variant_source_id", "sv_source_id", "pav_source_id", "promoter_source_id")


def dep_id(source_id: str) -> str:
    """Content-addressed dependency-cluster id for one source id.

    Same formula as the evidence-dependency module: ``dep_`` followed by the
    first 20 hex chars of the SHA-256 digest of the UTF-8 encoded source id.
    """
    return "dep_" + hashlib.sha256(source_id.encode("utf-8")).hexdigest()[:20]


def _cell_tokens(value: object) -> list[str]:
    return [token.strip() for token in str(value).split(";") if token.strip()]


def _sort_key(row: dict[str, object]) -> tuple[str, int, str]:
    return (str(row.get("chrom", "")), int(str(row.get("start", "")) or 0), str(row.get("gene_id", "")))


def _representative_key(row: dict[str, object]) -> tuple[int, int, int, str]:
    cds_length = int(str(row.get("cds_length", "")) or 0)
    protein_length = int(str(row.get("protein_length", "")) or 0)
    return (_MODEL_STATUS_RANK.get(str(row.get("model_status", "")), 4), -cds_length, -protein_length, str(row.get("gene_id", "")))


def _provider_gap(row: dict[str, object]) -> bool:
    """True when providers are applicable and the gene has no positive
    evidence while at least one provided table has no row for it.

    A mix of ``missing`` (table provided, no row for the gene) and
    ``not_evaluated`` (kind not provided) cells still counts as a gap; a
    fully ``not_evaluated`` gene (no providers at all) never does.
    """
    statuses = [str(row.get(cell, "")) for cell in _RESOURCE_STATUS_CELLS]
    return (
        str(row.get("provider_applicability", "")) == "applicable"
        and "present" not in statuses
        and "missing" in statuses
    )


def _rescue_trigger(row: dict[str, object]) -> bool:
    return (
        str(row.get("model_status", "")) in _RESCUE_TRIGGER_MODEL
        or str(row.get("annotation_status", "")) == "missing"
        or any(str(row.get(f"{kind}_status", "")) == "present" for kind in ("sv", "pav", "promoter"))
        or _provider_gap(row)
        or str(row.get("model_status", "")) == "not_evaluated"
    )


def _rescue_bucket(row: dict[str, object]) -> str:
    if any(str(row.get(f"{kind}_status", "")) == "present" for kind in ("sv", "pav", "promoter")):
        return "sv_pav_regulatory"
    if str(row.get("model_status", "")) in _RESCUE_TRIGGER_MODEL:
        return "missing_model"
    if str(row.get("annotation_status", "")) == "missing":
        return "weak_annotation"
    if _provider_gap(row):
        return "provider_gap"
    return "unclassified"


def _family_value_column(family_table: dict[str, object] | None) -> str:
    if not family_table:
        return ""
    from .tier0_inventory import _value_column

    columns = tuple(map(str, family_table["rows"][0].keys())) if family_table["rows"] else ()
    return _value_column("family", columns)


def _family_ids(member_rows: list[dict[str, object]], family_table: dict[str, object] | None) -> str:
    value_column = _family_value_column(family_table)
    if not family_table or not value_column:
        return ""
    members = {str(row.get("gene_id", "")) for row in member_rows}
    values = {
        str(record[value_column]).strip()
        for record in family_table["rows"]
        if str(record.get("gene_id", "")) in members and str(record.get(value_column, "")).strip()
    }
    return ";".join(sorted(values)) if values else ""


def _dependency_ids(member_rows: list[dict[str, object]]) -> str:
    ids: set[str] = set()
    for row in member_rows:
        for column in _DEPENDENCY_SOURCE_CELLS:
            for token in _cell_tokens(row.get(column, "")):
                ids.add(dep_id(token))
    return ";".join(sorted(ids)) if ids else ""


def build_tier1_clusters(
    inventory_rows: list[dict[str, object]],
    out_dir: str | Path,
    case_id: str,
    providers: dict[str, dict[str, object]] | None = None,
) -> dict[str, object]:
    """Assign every inventory gene to exactly one Tier 1 cluster.

    ``inventory_rows`` are the row dicts written by
    ``tier0_inventory.build_tier0_inventory`` (hand-built dicts carrying at
    least the ``TIER0_COLUMNS`` keys are accepted; missing ``cds_length`` /
    ``protein_length`` default to 0). ``providers`` maps kind to the dict from
    ``load_provider_table``; only ``family`` is consumed here: its value
    column values are the R0 cluster keys and, for all member rows, the
    ``family_ids`` cell.
    """
    out = runtime_path(out_dir)
    rows = [dict(row) for row in inventory_rows]
    rows.sort(key=_sort_key)
    family_table = (providers or {}).get("family")

    clusters: list[dict[str, object]] = []

    def add_cluster(
        member_rows: list[dict[str, object]],
        rule_key: str,
        cluster_key: str,
        allocation: str,
        reason: str,
        rescue_bucket: str = "not_applicable",
    ) -> None:
        members = sorted({str(row.get("gene_id", "")) for row in member_rows})
        cluster_id = "cls_" + hashlib.sha256("\x1f".join(members).encode("utf-8")).hexdigest()[:20]
        representative = min(member_rows, key=_representative_key)
        model_status = str(representative.get("model_status", ""))
        cds_length = int(str(representative.get("cds_length", "")) or 0)
        protein_length = int(str(representative.get("protein_length", "")) or 0)
        representative_reason = f"model_status={model_status};cds_length={cds_length};protein_length={protein_length}"
        representative_conflicts = str(representative.get("conflict_flags", ""))
        if representative_conflicts:
            representative_reason += f";conflict_flags={representative_conflicts}"
        member_conflicts = sorted(str(row.get("gene_id", "")) for row in member_rows if str(row.get("conflict_flags", "")))
        if member_conflicts:
            representative_reason += ";member_conflicts=" + ";".join(member_conflicts)
        clusters.append(
            {
                "cluster_id": cluster_id,
                "seed_reason": f"{rule_key}:{cluster_key}",
                "representative_gene_id": str(representative.get("gene_id", "")),
                "representative_reason": representative_reason,
                "member_gene_ids": ";".join(members),
                "member_count": len(members),
                "rule_key": rule_key,
                "cluster_key": cluster_key,
                "family_ids": _family_ids(member_rows, family_table),
                "domain_ids": ";".join(sorted({token for row in member_rows for token in _cell_tokens(row.get("domain_hits", ""))})),
                "orthogroup_ids": ";".join(sorted({token for row in member_rows for token in _cell_tokens(row.get("orthogroup", ""))})),
                "allocation": allocation,
                "allocation_reason": reason,
                "rescue_bucket": rescue_bucket,
                "dependency_cluster_ids": _dependency_ids(member_rows),
            }
        )

    # R0: shared family id (smallest family id wins for multi-family genes).
    family_value_column = _family_value_column(family_table)
    family_groups: dict[str, list[dict[str, object]]] = {}
    family_members: set[str] = set()
    if family_value_column:
        for row in rows:
            gene_id = str(row.get("gene_id", ""))
            family_values = sorted(
                {
                    str(record[family_value_column])
                    for record in family_table["rows"]
                    if str(record.get("gene_id", "")) == gene_id and str(record.get(family_value_column, "")).strip()
                }
            )
            if family_values:
                family_groups.setdefault(family_values[0], []).append(row)
                family_members.add(gene_id)
    for family_id in sorted(family_groups):
        add_cluster(family_groups[family_id], "family", family_id, "priority", "shared_family")
    pool = [row for row in rows if str(row.get("gene_id", "")) not in family_members]

    # R1: shared non-empty gene_name; a unique gene_name never clusters alone.
    name_groups: dict[str, list[dict[str, object]]] = {}
    for row in pool:
        name = str(row.get("gene_name", "")).strip()
        if name:
            name_groups.setdefault(name, []).append(row)
    shared_names = {name for name, group in name_groups.items() if len(group) >= 2}
    for name in sorted(shared_names):
        add_cluster(name_groups[name], "annotation", name, "priority", "shared_annotation")
    pool = [row for row in pool if str(row.get("gene_name", "")).strip() not in shared_names]

    # R2: shared smallest domain id.
    domain_groups: dict[str, list[dict[str, object]]] = {}
    for row in pool:
        hits = sorted(_cell_tokens(row.get("domain_hits", "")))
        if hits:
            domain_groups.setdefault(hits[0], []).append(row)
    for domain_id in sorted(domain_groups):
        add_cluster(domain_groups[domain_id], "domain", domain_id, "priority", "shared_domain")
    pool = [row for row in pool if not _cell_tokens(row.get("domain_hits", ""))]

    # R3: shared smallest orthogroup id.
    orthogroup_groups: dict[str, list[dict[str, object]]] = {}
    for row in pool:
        orthogroups = sorted(_cell_tokens(row.get("orthogroup", "")))
        if orthogroups:
            orthogroup_groups.setdefault(orthogroups[0], []).append(row)
    for orthogroup_id in sorted(orthogroup_groups):
        add_cluster(orthogroup_groups[orthogroup_id], "orthogroup", orthogroup_id, "defer", "shared_orthogroup")
    pool = [row for row in pool if not _cell_tokens(row.get("orthogroup", ""))]

    # R4: model-state rescue (one cluster per gene).
    for row in sorted(pool, key=_sort_key):
        if not _rescue_trigger(row):
            continue
        bucket = _rescue_bucket(row)
        add_cluster([row], "model_state_rescue", str(row.get("gene_id", "")), "rescue", f"model_state_rescue:{bucket}", rescue_bucket=bucket)
    pool = [row for row in pool if not _rescue_trigger(row)]

    # R5: singleton default.
    for row in sorted(pool, key=_sort_key):
        add_cluster([row], "singleton", str(row.get("gene_id", "")), "defer", "singleton_default")

    clusters.sort(key=lambda cluster: str(cluster["cluster_id"]))
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(clusters, columns=TIER1_COLUMNS).to_csv(out / "tier1_clusters.tsv", sep="\t", index=False, lineterminator="\n")

    dependency_ids = sorted({token for cluster in clusters for token in _cell_tokens(cluster.get("dependency_cluster_ids", ""))})
    manifest = {
        "schema_version": "1.0",
        "case_id": case_id,
        "cluster_count": len(clusters),
        "gene_count": len(rows),
        "genes_in_cluster_count": sum(int(cluster["member_count"]) for cluster in clusters),
        "allocation_counts": dict(sorted(Counter(str(cluster["allocation"]) for cluster in clusters).items())),
        "rescue_bucket_counts": dict(sorted(Counter(str(cluster["rescue_bucket"]) for cluster in clusters).items())),
        "rule_key_counts": dict(sorted(Counter(str(cluster["rule_key"]) for cluster in clusters).items())),
        "dependency_cluster_count": len(dependency_ids),
        "empty_dependency_source_count": sum(1 for cluster in clusters if not _cell_tokens(cluster.get("dependency_cluster_ids", ""))),
    }
    (out / "tier1_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
    )
    return manifest
