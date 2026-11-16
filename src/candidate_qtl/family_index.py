"""Genome-wide family index with bounded annotation-aware hit selection.

The historical implementation is retained byte-for-byte in
``family_index_legacy``. This facade replaces its broad symbol-only promotion
with one sequence-first equivalence rule and appends explicit selection
provenance to the generated index.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pandas as pd

from . import family_index_legacy as _legacy
from .core import runtime_path
from .reference_hit_selection import (
    MAX_IDENTITY_DEFICIT_PERCENTAGE_POINTS,
    MAX_QUERY_COVERAGE_DEFICIT,
    MAX_SUBJECT_COVERAGE_DEFICIT,
    MIN_EQUIVALENT_BITSCORE_RATIO,
    is_reviewed_subject,
    select_reference_hit,
)

# Preserve the complete public/private compatibility surface. Several historical
# tests and benchmark helpers intentionally import low-level family-index helpers.
for _export_name in dir(_legacy):
    if not _export_name.startswith("__"):
        globals()[_export_name] = getattr(_legacy, _export_name)

REFERENCE_SELECTION_COLUMNS = [
    "sequence_top_subject",
    "sequence_top_bitscore",
    "selected_reference_subject",
    "selected_reference_rank",
    "selected_bitscore_ratio",
    "reference_selection_reason",
    "equivalent_reference_hit_count",
    "reference_selection_policy",
]
_SELECTION_POLICY = "bounded_annotation_tiebreak_v1"


def _prepared_records(
    records: list[dict[str, object]],
    overrides: dict[str, dict[str, str]],
) -> list[dict[str, object]]:
    prepared: list[dict[str, object]] = []
    for record in records:
        copy = dict(record)
        subject = str(record.get("subject_id", ""))
        override = overrides.get(subject, {})
        symbol, symbol_source, _confidence = _legacy._best_symbol(
            subject,
            str(record.get("subject_title", "")),
            str(override.get("symbol", "")),
        )
        copy["reviewed"] = bool(record.get("reviewed_subject", False)) or is_reviewed_subject(subject)
        # Identifier-derived fallbacks such as a TrEMBL accession are not gene
        # symbols, even when they happen to pass the broad lexical helper.
        copy["usable_symbol"] = symbol_source != "none" and _legacy.is_usable_symbol(symbol)
        prepared.append(copy)
    return prepared


def _selection_decision(
    records: list[dict[str, object]],
    overrides: dict[str, dict[str, str]],
) -> tuple[dict[str, object], dict[str, object]]:
    """Apply one narrow sequence-equivalence rule to the production path.

    The previous selector could promote a usable-symbol record at half the top
    bitscore. That was safe against grossly distant hits, but it is broader than
    issue #647 permits for annotation-driven selection. All symbol, reviewed-
    record and informative-title preferences now operate only inside the same
    fixed homology-equivalence window.
    """
    prepared = _prepared_records(records, overrides)
    annotation = select_reference_hit(prepared)
    selected_index = int(annotation["selected_index"])
    sequence_top_index = int(annotation["sequence_top_index"])
    sequence_top = records[sequence_top_index]
    selected = records[selected_index]
    top_bitscore = float(sequence_top.get("bitscore", 0.0) or 0.0)
    selected_bitscore = float(selected.get("bitscore", 0.0) or 0.0)
    audit = {
        "sequence_top_index": sequence_top_index,
        "selected_index": selected_index,
        "sequence_top_subject": str(sequence_top.get("subject_id", "")),
        "sequence_top_bitscore": top_bitscore,
        "selected_reference_subject": str(selected.get("subject_id", "")),
        "selected_reference_rank": selected_index + 1,
        "selected_bitscore_ratio": (selected_bitscore / top_bitscore if top_bitscore else 1.0),
        "selection_reason": str(annotation["selection_reason"]),
        "equivalent_hit_count": int(annotation["equivalent_hit_count"]),
        "policy": _SELECTION_POLICY,
    }
    return selected, audit


def _pick_best_record(
    records: list[dict[str, object]],
    overrides: dict[str, dict[str, str]],
) -> dict[str, object]:
    return _selection_decision(records, overrides)[0]


# Functions defined in the retained implementation resolve globals in that
# module, so install the bounded selector there before exposing ``family_index``.
_legacy._pick_best_record = _pick_best_record


def _selection_audits(
    blast_hits_arg: str,
    *,
    id_map_arg: str | None,
    family_map_arg: str | None,
    family_hits: int,
    min_identity: float,
    min_query_coverage: float,
    max_evalue: float,
) -> dict[str, dict[str, object]]:
    hits = _legacy._read_blast_hits(runtime_path(blast_hits_arg))
    overrides = _legacy._family_overrides(family_map_arg)
    protein_gene = _legacy._protein_to_gene(id_map_arg)
    hits["gene_id"] = [protein_gene.get(value, value) for value in hits["query_id"].astype(str)]
    quality = (
        (hits["identity"] >= min_identity)
        & (hits["merged_query_coverage"] >= min_query_coverage)
        & (hits["evalue"] <= max_evalue)
    )
    hits = hits[quality].copy().sort_values(
        ["gene_id", "bitscore", "evalue", "reviewed_subject", "subject_id"],
        ascending=[True, False, True, False, True],
        kind="stable",
    )
    audits: dict[str, dict[str, object]] = {}
    for gene_id, group in hits.groupby("gene_id", sort=True):
        records = group.head(family_hits).to_dict("records")
        if records:
            _selected, audit = _selection_decision(records, overrides)
            audits[str(gene_id)] = audit
    return audits


def _append_selection_provenance(
    result: dict[str, object],
    audits: dict[str, dict[str, object]],
) -> None:
    index_path = Path(str(result["family_index"]))
    frame = pd.read_csv(index_path, sep="\t", dtype=str).fillna("")
    for column in REFERENCE_SELECTION_COLUMNS:
        frame[column] = ""
    for index, gene_id in frame["gene_id"].astype(str).items():
        audit = audits.get(gene_id)
        if audit is None:
            continue
        frame.at[index, "sequence_top_subject"] = audit["sequence_top_subject"]
        frame.at[index, "sequence_top_bitscore"] = str(audit["sequence_top_bitscore"])
        frame.at[index, "selected_reference_subject"] = audit["selected_reference_subject"]
        frame.at[index, "selected_reference_rank"] = str(audit["selected_reference_rank"])
        frame.at[index, "selected_bitscore_ratio"] = f"{float(audit['selected_bitscore_ratio']):.6f}"
        frame.at[index, "reference_selection_reason"] = audit["selection_reason"]
        frame.at[index, "equivalent_reference_hit_count"] = str(audit["equivalent_hit_count"])
        frame.at[index, "reference_selection_policy"] = audit["policy"]
    frame.to_csv(index_path, sep="\t", index=False, lineterminator="\n")

    manifest_path = Path(str(result["manifest"]))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict):
        raise ValueError("family-index manifest must be a JSON object")
    reasons = Counter(str(audit["selection_reason"]) for audit in audits.values())
    manifest["reference_hit_selection"] = {
        "policy": _SELECTION_POLICY,
        "sequence_first": True,
        "annotation_is_causal_evidence": False,
        "legacy_half_bitscore_symbol_promotion_active": False,
        "bounds": {
            "min_bitscore_ratio": MIN_EQUIVALENT_BITSCORE_RATIO,
            "max_identity_deficit_percentage_points": MAX_IDENTITY_DEFICIT_PERCENTAGE_POINTS,
            "max_query_coverage_deficit": MAX_QUERY_COVERAGE_DEFICIT,
            "max_subject_coverage_deficit": MAX_SUBJECT_COVERAGE_DEFICIT,
        },
        "generic_title_policy": (
            "only blank/uncharacterized/hypothetical/unknown titles are generic; "
            "putative names remain informative"
        ),
        "reason_counts": dict(sorted(reasons.items())),
        "selected_non_sequence_top_count": sum(
            int(audit["selected_index"] != audit["sequence_top_index"])
            for audit in audits.values()
        ),
        "columns": REFERENCE_SELECTION_COLUMNS,
    }
    family_definition = manifest.get("family_definition")
    if isinstance(family_definition, dict):
        family_definition["reference_record_selection"] = (
            "maximum bitscore remains the sequence anchor; usable-symbol, reviewed or informative records may "
            "replace it only inside the fixed #647 homology-equivalence window, with both records retained in the TSV"
        )
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def family_index(
    blast_hits_arg: str,
    out_dir_arg: str,
    genes_arg: str | None = None,
    id_map_arg: str | None = None,
    family_map_arg: str | None = None,
    reference_label: str = "",
    reference_species: str = "",
    family_hits: int = _legacy.DEFAULT_FAMILY_HITS,
    min_identity: float = 25.0,
    min_query_coverage: float = 0.5,
    max_evalue: float = 1e-5,
    full_hash: bool = False,
    case_id: str = "",
    index_scope: str = "reference_genome",
) -> dict[str, object]:
    result = _legacy.family_index(
        blast_hits_arg,
        out_dir_arg,
        genes_arg=genes_arg,
        id_map_arg=id_map_arg,
        family_map_arg=family_map_arg,
        reference_label=reference_label,
        reference_species=reference_species,
        family_hits=family_hits,
        min_identity=min_identity,
        min_query_coverage=min_query_coverage,
        max_evalue=max_evalue,
        full_hash=full_hash,
        case_id=case_id,
        index_scope=index_scope,
    )
    audits = _selection_audits(
        blast_hits_arg,
        id_map_arg=id_map_arg,
        family_map_arg=family_map_arg,
        family_hits=family_hits,
        min_identity=min_identity,
        min_query_coverage=min_query_coverage,
        max_evalue=max_evalue,
    )
    _append_selection_provenance(result, audits)
    return result


# Ensure direct imports and legacy internal calls both use the wrapped producer.
globals()["family_index"] = family_index
