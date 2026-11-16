"""Shared evidence I/O, vocabulary, and provenance helpers.

This module is the single public surface for the small building blocks that
several evidence-producing modules need (candidate/evidence table loading, the
canonical evidence row, best-hit parsing, and the shared vocabulary constants).
Keeping them here -- rather than importing private ``_``-prefixed helpers out of
``evidence`` -- gives ranking, benchmark, discovery, and orthology code one
explicit dependency boundary instead of reaching into another module's internals
(issue #152).
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from . import schemas
from .core import runtime_path

# Canonical column order and enum vocabularies live in ``schemas`` (the single
# source of truth); re-exported here for the call sites that import them.
EVIDENCE_COLUMNS = schemas.EVIDENCE_COLUMNS
CONFIDENCE_MULTIPLIERS = schemas.CONFIDENCE_MULTIPLIERS
STRONG_SINGLE_CODES = {"direct_model_genetic", "ortholog_model_genetic", "known_qtg"}


def annotation_aliases(identifier: str) -> set[str]:
    aliases = {identifier}
    aliases.update(token for token in identifier.split("|") if token)
    for token in list(aliases):
        aliases.add(token.rsplit(".", 1)[-1])
        aliases.add(re.sub(r"-T\d+$", "", token))
        aliases.add(re.sub(r"\.\d+$", "", token))
    return {alias for alias in aliases if alias}


def candidate_frame(path_arg: str) -> pd.DataFrame:
    path = runtime_path(path_arg)
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    required = {"gene_id"}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"candidate table is missing required columns: {', '.join(missing)}")
    if frame["gene_id"].str.strip().eq("").any() or frame["gene_id"].duplicated().any():
        raise ValueError("candidate gene_id values must be non-empty and unique")
    return frame


def id_map(candidates: pd.DataFrame) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for _, row in candidates.iterrows():
        gene_id = str(row["gene_id"])
        mapping[gene_id] = gene_id
        transcript_id = str(row.get("representative_transcript", "")).strip()
        if transcript_id:
            mapping[transcript_id] = gene_id
    return mapping


def evidence_row(
    gene_id: str,
    code: str,
    direction: str,
    confidence: str,
    value: object,
    source_type: str,
    source_id: object,
    description: str,
    *,
    evidence_kind: str = "",
    root_source_id: str = "",
) -> dict[str, str]:
    row = {
        "gene_id": gene_id,
        "evidence_code": code,
        "direction": direction,
        "confidence": confidence,
        "value": str(value),
        "source_type": source_type,
        "source_id": str(source_id),
        "description": description,
    }
    if evidence_kind:
        row["evidence_kind"] = evidence_kind
    if root_source_id:
        row["root_source_id"] = root_source_id
    return row


def best_hits(path: Path) -> dict[str, list[dict[str, object]]]:
    grouped: dict[str, list[dict[str, object]]] = {}
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 12:
                continue
            grouped.setdefault(fields[0], []).append({"subject": fields[1], "evalue": float(fields[10]), "bitscore": float(fields[11])})
    for hits in grouped.values():
        hits.sort(key=lambda hit: (-float(hit["bitscore"]), float(hit["evalue"]), str(hit["subject"])))
    return grouped
