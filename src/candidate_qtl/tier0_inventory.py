"""Deterministic Tier 0 inventory (issue #375, deterministic budget engine).

Tier 0 consumes the extract outputs (``candidates.tsv``, ``transcripts.tsv``,
``candidate_proteins.fa``, ``extract_manifest.json``) plus user-provided
provider tables (expression, domains, orthogroup, family, variant, sv, pav,
promoter) and produces one deterministic row per candidate gene.

Row contract: every row dict carries ``TIER0_COLUMNS`` plus the two auxiliary
keys ``cds_length`` and ``protein_length`` (used by Tier 1 representative
selection). ``tier0_inventory.tsv`` writes ``TIER0_COLUMNS`` with a trailing
``cds_length`` column; ``tier0_manifest.json`` is the returned manifest.

Determinism contract (issue #375): no wall-clock values, no RNG, no timestamps.
Multi-value cells and every manifest list are sorted; the manifest is written
with ``json.dumps(..., ensure_ascii=False, indent=2, sort_keys=True) + "\\n"``
and the TSV with ``pandas .to_csv(sep="\\t", index=False, lineterminator="\\n")``.

Defensive reading: extract outputs only a subset of the canonical candidate
columns, so every candidate field is fetched with ``.get(column, "")`` and the
provider value/source-id cells are derived from the columns actually present
(see ``PROVIDER_VALUE_COLUMNS`` for the deterministic first-match priority).

Status cells are drawn from closed vocabularies (``TIER0_MODEL_STATUSES``,
``TIER0_ANNOTATION_STATUSES``, ``TIER0_EVIDENCE_STATUSES``,
``TIER0_APPLICABILITY``); a ``*_status`` or ``provider_applicability`` cell is
never empty. A provider table that is absent for a kind yields
``not_evaluated`` with empty evidence/provider/source-id cells.

Evidence-status semantics (audit finding F): a kind's status is ``present``
only when the gene has at least one row whose evidence value cell is non-empty
after stripping; rows with empty/whitespace-only values make the status
``missing``. A ``*_source_id`` value alone is never positive evidence, and
``provider_status[kind]["genes_covered"]`` counts only genes with at least one
non-empty evidence value.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pandas as pd

from .core import file_fingerprint, fingerprint_is_weak, open_text, runtime_path

TIER0_COLUMNS = (
    "gene_id",
    "chrom",
    "start",
    "end",
    "strand",
    "transcript_count",
    "representative_transcript",
    "model_status",
    "protein_length",
    "translation_status",
    "gene_name",
    "gene_biotype",
    "symbol_status",
    "description",
    "dbxref",
    "ontology_terms",
    "annotation_status",
    "domain_hits",
    "domain_provider",
    "domain_status",
    "orthogroup",
    "orthogroup_provider",
    "orthogroup_status",
    "expression_evidence",
    "expression_provider",
    "expression_status",
    "variant_evidence",
    "variant_source_id",
    "variant_provider",
    "variant_status",
    "sv_evidence",
    "sv_source_id",
    "sv_provider",
    "sv_status",
    "pav_evidence",
    "pav_source_id",
    "pav_provider",
    "pav_status",
    "promoter_evidence",
    "promoter_source_id",
    "promoter_provider",
    "promoter_status",
    "provider_applicability",
    "conflict_flags",
)

TIER0_MODEL_STATUSES = frozenset({"complete_orf", "partial_or_unverified_orf", "no_cds", "not_evaluated"})
TIER0_ANNOTATION_STATUSES = frozenset({"present", "partial", "missing", "not_evaluated"})
TIER0_EVIDENCE_STATUSES = frozenset({"present", "missing", "not_evaluated"})
TIER0_APPLICABILITY = frozenset({"applicable", "not_evaluated"})

# Provider kinds accepted by load_provider_table. ``family`` is consumed by
# Tier 1 only; it still participates in fingerprints and provider_status.
PROVIDER_KINDS = frozenset({"expression", "domains", "orthogroup", "family", "variant", "sv", "pav", "promoter"})

# Evidence kinds that produce Tier 0 per-gene cells (family has no Tier 0 cell).
EVIDENCE_KINDS = ("domains", "orthogroup", "expression", "variant", "sv", "pav", "promoter")

# Deterministic first-match priority for the evidence value column of each kind.
PROVIDER_VALUE_COLUMNS: dict[str, tuple[str, ...]] = {
    "expression": ("expression", "value", "tpm", "fpkm"),
    "domains": ("domain_id", "domain", "domain_name"),
    "orthogroup": ("orthogroup",),
    "variant": ("effect", "consequence", "variant"),
    "sv": ("sv_type", "type", "sv"),
    "pav": ("pav_type", "type", "pav"),
    "promoter": ("motif", "evidence", "promoter"),
    "family": ("family_id", "family", "value"),
}

# Tier 0 cell names per evidence kind: (evidence, provider, status, source_id).
# Only variant/sv/pav/promoter expose a ``*_source_id`` cell in TIER0_COLUMNS.
_EVIDENCE_CELL_COLUMNS: dict[str, tuple[str, str, str, str | None]] = {
    "domains": ("domain_hits", "domain_provider", "domain_status", None),
    "orthogroup": ("orthogroup", "orthogroup_provider", "orthogroup_status", None),
    "expression": ("expression_evidence", "expression_provider", "expression_status", None),
    "variant": ("variant_evidence", "variant_provider", "variant_status", "variant_source_id"),
    "sv": ("sv_evidence", "sv_provider", "sv_status", "sv_source_id"),
    "pav": ("pav_evidence", "pav_provider", "pav_status", "pav_source_id"),
    "promoter": ("promoter_evidence", "promoter_provider", "promoter_status", "promoter_source_id"),
}

# annotation_quality -> model_status; any other (or absent) value is not_evaluated.
_MODEL_QUALITY_TO_STATUS = {
    "complete_orf": "complete_orf",
    "partial_or_unverified_orf": "partial_or_unverified_orf",
    "no_cds": "no_cds",
    "invalid_cds_phase": "partial_or_unverified_orf",
}


def _gene_column_value_sets(rows: list[dict[str, str]]) -> dict[str, dict[str, set[str]]]:
    """Map gene_id -> column -> set of distinct non-empty stripped values (table-local)."""
    grouped: dict[str, dict[str, set[str]]] = {}
    for row in rows:
        gene_id = str(row.get("gene_id", ""))
        columns = grouped.setdefault(gene_id, {})
        for column, value in row.items():
            if column == "gene_id" or str(value).strip() == "":
                continue
            columns.setdefault(column, set()).add(str(value).strip())
    return grouped


def _value_column(kind: str, columns: tuple[str, ...]) -> str:
    return next((column for column in PROVIDER_VALUE_COLUMNS[kind] if column in columns), "")


def load_provider_table(kind: str, path: str | Path, provider: str, version: str) -> dict[str, object]:
    """Read one user-provided provider TSV into a reusable table dict.

    Returns ``{"kind", "provider", "version", "path", "rows", "gene_id_set",
    "conflicts"}``. ``rows`` is the table verbatim as a list of string-valued
    dicts (file order); ``gene_id_set`` is the sorted distinct non-empty gene
    ids; ``conflicts`` counts the genes with >=2 distinct non-empty values in
    any non-gene_id column. Raises ``ValueError`` (fail-closed) for an
    unsupported kind, a missing file, a missing ``gene_id`` column, or when
    none of the kind's documented evidence columns is present.
    """
    if kind not in PROVIDER_KINDS:
        raise ValueError(f"unsupported provider kind {kind!r}; expected one of {sorted(PROVIDER_KINDS)}")
    table_path = runtime_path(path)
    if not table_path.exists():
        raise ValueError(f"provider table not found: {table_path}")
    with open_text(table_path) as handle:
        frame = pd.read_csv(handle, sep="\t", dtype=str).fillna("")
    if "gene_id" not in frame.columns:
        raise ValueError(f"{kind} provider table requires a gene_id column; observed columns: {sorted(map(str, frame.columns))}")
    value_column = _value_column(kind, tuple(map(str, frame.columns)))
    if not value_column:
        raise ValueError(
            f"{kind} provider table requires one of the evidence columns {PROVIDER_VALUE_COLUMNS[kind]}; "
            f"observed columns: {sorted(map(str, frame.columns))}"
        )
    rows: list[dict[str, str]] = [dict(record) for record in frame.to_dict("records")]
    grouped = _gene_column_value_sets(rows)
    conflicts = sum(1 for values_by_column in grouped.values() if any(len(values) >= 2 for values in values_by_column.values()))
    gene_id_set = sorted({str(row.get("gene_id", "")) for row in rows if str(row.get("gene_id", ""))})
    return {
        "kind": kind,
        "provider": provider,
        "version": version,
        "path": str(table_path),
        "rows": rows,
        "gene_id_set": gene_id_set,
        "conflicts": conflicts,
    }


def build_tier0_inventory(
    extract_dir: str | Path,
    out_dir: str | Path,
    providers: dict[str, dict[str, object]],
    case_id: str,
) -> dict[str, object]:
    """Build the Tier 0 per-gene inventory from an extract directory.

    ``providers`` maps a provider kind to the dict returned by
    ``load_provider_table`` (an absent or ``None`` value means "not provided").
    Writes ``tier0_inventory.tsv`` and ``tier0_manifest.json`` into ``out_dir``
    and returns the manifest. Rows are sorted by (chrom, start as int, gene_id)
    and one row is emitted per candidate gene; ``candidates.tsv`` must exist
    (``ValueError`` otherwise). transcripts.tsv / candidate_proteins.fa /
    extract_manifest.json are fingerprinted but never required.
    """
    extract = runtime_path(extract_dir)
    out = runtime_path(out_dir)
    candidates_path = extract / "candidates.tsv"
    if not candidates_path.exists():
        raise ValueError(f"extract candidates table missing: {candidates_path}")
    frame = pd.read_csv(candidates_path, sep="\t", dtype=str).fillna("")
    candidates_columns_observed = sorted(str(column) for column in frame.columns)
    rows = [dict(record) for record in frame.to_dict("records")]
    duplicate_gene_ids = sorted(
        {gene_id for gene_id, count in Counter(str(row.get("gene_id", "")) for row in rows).items() if count > 1}
    )
    if duplicate_gene_ids:
        raise ValueError(f"candidates.tsv contains duplicate gene_id values: {', '.join(duplicate_gene_ids)}")
    rows.sort(key=lambda row: (str(row.get("chrom", "")), int(str(row.get("start", "")) or 0), str(row.get("gene_id", ""))))
    provided_kinds = {kind for kind, table in providers.items() if table is not None}
    table_columns: dict[str, tuple[str, ...]] = {}
    value_columns: dict[str, str] = {}
    for kind in EVIDENCE_KINDS:
        table = providers.get(kind)
        if table is None:
            table_columns[kind] = ()
            value_columns[kind] = ""
            continue
        columns = tuple(map(str, table["rows"][0].keys())) if table["rows"] else ()
        table_columns[kind] = columns
        value_columns[kind] = _value_column(kind, columns)

    inventory: list[dict[str, str]] = []
    for row in rows:
        gene_id = str(row.get("gene_id", ""))
        gene_name = str(row.get("gene_name", ""))
        description = str(row.get("description", ""))
        dbxref = str(row.get("dbxref", ""))
        ontology_terms = str(row.get("ontology_terms", ""))
        gene_biotype = str(row.get("gene_biotype", ""))
        if description or dbxref or ontology_terms:
            annotation_status = "present"
        elif gene_name or gene_biotype:
            annotation_status = "partial"
        else:
            annotation_status = "missing"
        entry: dict[str, str] = {
            "gene_id": gene_id,
            "chrom": str(row.get("chrom", "")),
            "start": str(row.get("start", "")),
            "end": str(row.get("end", "")),
            "strand": str(row.get("strand", "")),
            "transcript_count": str(row.get("transcript_count", "")),
            "representative_transcript": str(row.get("representative_transcript", "")),
            "model_status": _MODEL_QUALITY_TO_STATUS.get(str(row.get("annotation_quality", "")), "not_evaluated"),
            "protein_length": str(row.get("protein_length", "0")),
            "translation_status": str(row.get("translation_status", "")),
            "gene_name": gene_name,
            "gene_biotype": gene_biotype,
            "symbol_status": "present" if gene_name else "missing",
            "description": description,
            "dbxref": dbxref,
            "ontology_terms": ontology_terms,
            "annotation_status": annotation_status,
            "provider_applicability": "applicable" if provided_kinds else "not_evaluated",
            "conflict_flags": "",
        }
        conflict_flags: list[str] = []
        for kind in EVIDENCE_KINDS:
            table = providers.get(kind)
            evidence_cell, provider_cell_name, status_cell, source_id_cell_name = _EVIDENCE_CELL_COLUMNS[kind]
            if table is None:
                entry[evidence_cell] = ""
                entry[provider_cell_name] = ""
                entry[status_cell] = "not_evaluated"
                if source_id_cell_name is not None:
                    entry[source_id_cell_name] = ""
                continue
            gene_rows = [record for record in table["rows"] if str(record.get("gene_id", "")) == gene_id]
            value_column = value_columns[kind]
            if value_column:
                evidence = ";".join(
                    sorted({str(record[value_column]).strip() for record in gene_rows if str(record.get(value_column, "")).strip()})
                )
            else:
                evidence = ""
            entry[evidence_cell] = evidence
            entry[status_cell] = "present" if evidence else "missing"
            entry[provider_cell_name] = f"{table['provider']}@{table['version']}"
            if source_id_cell_name is not None:
                if "source_id" in table_columns[kind]:
                    entry[source_id_cell_name] = ";".join(
                        sorted({str(record["source_id"]) for record in gene_rows if str(record.get("source_id", "")).strip()})
                    )
                else:
                    entry[source_id_cell_name] = ""
            for column, values in _gene_column_value_sets(gene_rows).get(gene_id, {}).items():
                if len(values) >= 2:
                    conflict_flags.append(f"{kind}:{gene_id}:{column}")
        entry["conflict_flags"] = ";".join(sorted(conflict_flags))
        entry["cds_length"] = str(row.get("cds_length", "0"))
        inventory.append(entry)

    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(inventory, columns=[*TIER0_COLUMNS, "cds_length"]).to_csv(
        out / "tier0_inventory.tsv", sep="\t", index=False, lineterminator="\n"
    )

    input_fingerprints: dict[str, object] = {}
    for key, filename in (
        ("candidates", "candidates.tsv"),
        ("transcripts", "transcripts.tsv"),
        ("candidate_proteins", "candidate_proteins.fa"),
        ("extract_manifest", "extract_manifest.json"),
    ):
        path = extract / filename
        if path.exists():
            input_fingerprints[key] = file_fingerprint(path, False)
    for kind in sorted(provided_kinds):
        input_fingerprints[kind] = file_fingerprint(runtime_path(providers[kind]["path"]), False)
    weak_input_fingerprints = sorted(
        str(fingerprint["path"]) for fingerprint in input_fingerprints.values() if fingerprint_is_weak(fingerprint)
    )
    provider_status: dict[str, object] = {}
    candidate_gene_ids = {str(row.get("gene_id", "")) for row in rows}
    for kind in sorted(provided_kinds):
        table = providers[kind]
        value_column = value_columns.get(kind) or (
            _value_column(kind, tuple(map(str, table["rows"][0].keys()))) if table["rows"] else ""
        )
        covered = sorted(
            {
                str(record.get("gene_id", ""))
                for record in table["rows"]
                if str(record.get("gene_id", "")) in candidate_gene_ids
                and value_column
                and str(record.get(value_column, "")).strip()
            }
        )
        provider_status[kind] = {
            "provided": True,
            "provider": str(table["provider"]),
            "version": str(table["version"]),
            "rows": len(table["rows"]),
            "genes_covered": len(covered),
            "conflicts": int(table["conflicts"]),
            # the actually-used evidence value column; "" when the table was
            # provided but had no rows (nothing consumed) or the fallback was
            # exhausted, making a secondary-column fallback visible.
            "value_column": value_column,
        }
    counts: dict[str, dict[str, int]] = {}
    for column in (
        "model_status",
        "annotation_status",
        "symbol_status",
        "domain_status",
        "orthogroup_status",
        "expression_status",
        "variant_status",
        "sv_status",
        "pav_status",
        "promoter_status",
        "provider_applicability",
    ):
        counts[column] = dict(sorted(Counter(entry[column] for entry in inventory).items()))
    manifest = {
        "schema_version": "1.0",
        "case_id": case_id,
        "gene_count": len(inventory),
        "candidates_columns_observed": candidates_columns_observed,
        "input_fingerprints": input_fingerprints,
        "weak_input_fingerprints": weak_input_fingerprints,
        "provider_status": provider_status,
        "counts": counts,
    }
    (out / "tier0_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
    )
    return manifest
