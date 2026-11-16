"""Normalized provenance DAG ledger (issue #231: explicit route/source dependence).

Different routes (A/B/C), retrieval passes, and annotation imports each bring
their own candidate lists into the same case. When several of them are derived
from the same underlying source, the workflow must be able to say so explicitly
instead of silently counting the shared source once per route. This module
builds a small directed acyclic graph (DAG) of provenance nodes and edges and
materialises it as a normalized ledger:

- ``nodes.tsv``: one row per provenance node (experiment, annotation, review,
  literature, or route output).
- ``edges.tsv``: one row per child -> parent derivation link.
- ``provenance_graph.tsv``: the node+edge join view for humans.
- ``provenance_ledger_manifest.json``: counts, shared-derivation summary, and
  the scientific-limits statement that a shared source is never counted as
  independent evidence per route.

Every node belongs to one of the canonical kinds
(``primary_experiment`` / ``database_annotation`` / ``review`` /
``route_output`` / ``literature``) and every edge to one of
(``derived_from`` / ``annotates`` / ``reviewed_in`` / ``produced``). Sources
declared without a parent are auto-registered as ``primary_experiment`` root
nodes whose ``node_id`` is the ``source_id`` itself.
"""

from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import pandas as pd

from . import schemas
from .core import file_fingerprint, runtime_path

LEDGER_VERSION = "1.0.0"

# Default edge kind per child node kind, used when the entries TSV does not
# declare an explicit edge_kind column. The vocabulary mirrors the canonical
# relationships in PROVENANCE_EDGE_KINDS.
DEFAULT_EDGE_KIND_BY_CHILD = {
    "primary_experiment": "produced",
    "database_annotation": "annotates",
    "review": "reviewed_in",
    "literature": "reviewed_in",
    "route_output": "derived_from",
}

REQUIRED_ENTRY_COLUMNS = ["node_kind", "source_id", "label", "parent_node"]

# Manifest wording for the double-counting hazard the ledger exists to expose.
SHARED_SOURCE_LIMIT = (
    "a source shared by multiple routes is counted once, never as "
    "independent evidence per route"
)


def _node_id_for(node_kind: str, source_id: str) -> str:
    """Deterministic auto-generated node_id for a child node."""
    return f"{node_kind}:{source_id}"


def _load_entries(entries_path: Path) -> pd.DataFrame:
    frame = pd.read_csv(entries_path, sep="\t", dtype=str).fillna("")
    missing = [column for column in REQUIRED_ENTRY_COLUMNS if column not in frame.columns]
    if missing:
        raise ValueError(f"entries TSV requires columns: {', '.join(missing)}")
    for column in REQUIRED_ENTRY_COLUMNS:
        frame[column] = frame[column].astype(str).str.strip()
    return frame


def provenance_ledger(
    entries_tsv: str,
    ledger_dir: str,
    out_dir: str,
    case_id: str = "",
) -> dict[str, object]:
    """Build, validate, and write a normalized provenance DAG ledger.

    Parameters
    ----------
    entries_tsv:
        TSV with one row per node: ``node_kind``, ``source_id``, ``label`` and
        optional ``parent_node`` (the parent's ``node_id``; blank means root).
        Optional ``node_id`` (explicit id, validated for uniqueness) and
        ``source_ref`` columns are honoured when present.
    ledger_dir:
        Directory receiving the normalized ``nodes.tsv`` / ``edges.tsv``.
    out_dir:
        Directory receiving ``provenance_graph.tsv`` and the manifest JSON.
    case_id:
        Case identifier recorded into the manifest, if any.

    Raises
    ------
    ValueError
        On duplicate ``node_id``, an edge to a parent that does not exist, or a
        child that declares itself as its own parent.
    """
    entries_path = runtime_path(entries_tsv)
    ledger = runtime_path(ledger_dir)
    out = runtime_path(out_dir)
    ledger.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    entries = _load_entries(entries_path)

    node_ids: dict[str, str] = {}  # node_id -> node_kind
    for index, row in entries.iterrows():
        node_kind = str(row["node_kind"])
        source_id = str(row["source_id"])
        if node_kind not in schemas.PROVENANCE_NODE_KINDS:
            raise ValueError(f"unknown node_kind {node_kind!r}; expected one of {sorted(schemas.PROVENANCE_NODE_KINDS)}")
        if not source_id:
            raise ValueError("source_id must be non-empty")
        explicit = str(row.get("node_id", "")).strip()
        if explicit:
            node_id = explicit
        elif not str(row["parent_node"]):
            # Parent-less sources are primary-experiment-equivalent roots whose
            # node_id is the source_id itself.
            node_id = source_id
        else:
            node_id = _node_id_for(node_kind, source_id)
        if node_id in node_ids:
            raise ValueError(f"duplicate node_id {node_id!r} in entries TSV")
        node_ids[node_id] = node_kind

    # Second pass resolves edges, which may reference parents declared later.
    node_rows: list[dict[str, object]] = []
    edge_rows: list[dict[str, object]] = []
    declared_edge_kind = "edge_kind" in entries.columns
    for index, row in entries.iterrows():
        node_kind = str(row["node_kind"])
        source_id = str(row["source_id"])
        parent = str(row["parent_node"])
        explicit = str(row.get("node_id", "")).strip()
        node_id = explicit or (source_id if not parent else _node_id_for(node_kind, source_id))
        node_rows.append({
            "node_id": node_id,
            "node_kind": node_kind,
            "source_id": source_id,
            "label": str(row.get("label", "")).strip(),
            "source_ref": str(row.get("source_ref", "")).strip(),
        })
        if not parent:
            continue
        if parent not in node_ids:
            raise ValueError(f"edge from {node_id!r} references unknown parent {parent!r}")
        if parent == node_id:
            raise ValueError(f"self-loop edge: child node_id {node_id!r} must differ from its parent")
        edge_kind = str(row.get("edge_kind", "")).strip() if declared_edge_kind else ""
        if not edge_kind:
            edge_kind = DEFAULT_EDGE_KIND_BY_CHILD[node_kind]
        if edge_kind not in schemas.PROVENANCE_EDGE_KINDS:
            raise ValueError(f"unknown edge_kind {edge_kind!r}; expected one of {sorted(schemas.PROVENANCE_EDGE_KINDS)}")
        edge_rows.append({"child_node": node_id, "parent_node": parent, "edge_kind": edge_kind})

    nodes = pd.DataFrame(node_rows, columns=schemas.PROVENANCE_NODE_COLUMNS)
    edges = pd.DataFrame(edge_rows, columns=schemas.PROVENANCE_EDGE_COLUMNS)
    schemas.validate(nodes, schemas.PROVENANCE_NODE_SCHEMA)
    schemas.validate(edges, schemas.PROVENANCE_EDGE_SCHEMA)

    nodes_path = ledger / "nodes.tsv"
    edges_path = ledger / "edges.tsv"
    nodes.sort_values("node_id", kind="stable").to_csv(nodes_path, sep="\t", index=False, lineterminator="\n")
    edges.sort_values(["child_node", "parent_node"], kind="stable").to_csv(edges_path, sep="\t", index=False, lineterminator="\n")

    # Node+edge join view: parent_node/edge_kind stay blank for root nodes.
    graph = nodes.merge(
        edges.rename(columns={"child_node": "node_id"}), how="left", on="node_id"
    ).fillna("")
    graph = graph[schemas.PROVENANCE_GRAPH_COLUMNS]
    schemas.validate(graph, schemas.PROVENANCE_GRAPH_SCHEMA)
    graph_path = out / "provenance_graph.tsv"
    graph.sort_values("node_id", kind="stable").to_csv(graph_path, sep="\t", index=False, lineterminator="\n")

    shared = shared_derivation_paths(ledger)
    shared_sources = [entry for entry in shared if len(entry["route_outputs"]) >= 2]
    manifest = {
        "schema_version": "1.0",
        "provenance_ledger_version": LEDGER_VERSION,
        "case_id": case_id,
        "node_count": len(nodes),
        "edge_count": len(edges),
        "node_kind_counts": {str(key): int(value) for key, value in nodes["node_kind"].value_counts().sort_index().items()},
        "edge_kind_counts": {str(key): int(value) for key, value in edges["edge_kind"].value_counts().sort_index().items()},
        "shared_derivation_source_count": len(shared_sources),
        "shared_derivation_pair_count": sum(len(entry["pairs"]) for entry in shared_sources),
        "files": {
            "nodes": str(nodes_path),
            "edges": str(edges_path),
            "graph": str(graph_path),
            "manifest": str(out / "provenance_ledger_manifest.json"),
        },
        "entries": file_fingerprint(entries_path),
        "scientific_limits": [
            SHARED_SOURCE_LIMIT,
            "parent-less sources are auto-registered as primary_experiment roots",
            "the ledger records dependence explicitly; it does not score evidence",
        ],
    }
    manifest_path = out / "provenance_ledger_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "nodes": str(nodes_path),
        "edges": str(edges_path),
        "graph": str(graph_path),
        "manifest": str(manifest_path),
        "node_count": len(nodes),
        "edge_count": len(edges),
        "shared_derivation_source_count": len(shared_sources),
        "shared_derivation_pair_count": sum(len(entry["pairs"]) for entry in shared_sources),
    }


def shared_derivation_paths(ledger_dir: str) -> list[dict[str, object]]:
    """List pairs of ``route_output`` nodes that converge on the same source.

    Two route outputs share a derivation path when walking each output up its
    parent edges reaches the same root source. Each returned entry carries the
    root ``source_id``, the converging ``route_output`` node ids, and every
    unordered pair of those outputs. A shared source must not be counted once
    per route, so these pairs are exactly the structures issue #231 asks to
    surface ("A/B/C are all derived from the same source").
    """
    ledger = runtime_path(ledger_dir)
    nodes = pd.read_csv(ledger / "nodes.tsv", sep="\t", dtype=str).fillna("")
    edges = pd.read_csv(ledger / "edges.tsv", sep="\t", dtype=str).fillna("")
    source_by_id = dict(zip(nodes["node_id"], nodes["source_id"]))
    kind_by_id = dict(zip(nodes["node_id"], nodes["node_kind"]))
    parent_by_child: dict[str, str] = {}
    for _, edge in edges.iterrows():
        child, parent = str(edge["child_node"]), str(edge["parent_node"])
        if child in parent_by_child:
            raise ValueError(f"multiple parents for child {child!r} are not supported")
        parent_by_child[child] = parent

    def root_source(node_id: str) -> str:
        seen: set[str] = set()
        current = node_id
        while current in parent_by_child and current not in seen:
            seen.add(current)
            current = parent_by_child[current]
        return source_by_id.get(current, source_by_id.get(node_id, ""))

    groups: dict[str, list[str]] = {}
    for node_id, kind in kind_by_id.items():
        if kind != "route_output":
            continue
        groups.setdefault(root_source(node_id), []).append(node_id)

    result: list[dict[str, object]] = []
    for source_id, outputs in sorted(groups.items()):
        outputs = sorted(outputs)
        if len(outputs) < 2:
            continue
        result.append({
            "source_id": source_id,
            "route_outputs": outputs,
            "pairs": [list(pair) for pair in combinations(outputs, 2)],
        })
    return result


def render_ledger_report(ledger_dir: str) -> str:
    """Human-readable summary of the ledger's shared-derivation structure."""
    shared = shared_derivation_paths(ledger_dir)
    lines = ["provenance ledger", "-----------------"]
    if shared:
        for entry in shared:
            lines.append(f"- shared source {entry['source_id']!r} feeds route outputs "
                         f"{', '.join(entry['route_outputs'])}")
    else:
        lines.append("- no route outputs converge on a shared source")
    lines.append(f"- shared derivation sources: {len(shared)}")
    lines.append(f"- {SHARED_SOURCE_LIMIT}.")
    return "\n".join(lines)
