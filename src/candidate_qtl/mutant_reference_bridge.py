"""Project curated mutant evidence onto the Route B reference-function contract.

Issues #478 / #525.  The chain ``interval gene -> reference gene -> that
reference gene's reported mutant phenotype`` must stay traceable at accession
level.  The LORE1/Tnt1 normalizers already emit accession-native
``normalized_functional_evidence.tsv`` rows (``FUNCTION_COLUMNS``; Gifu/R108
backgrounds).  This module projects those rows onto the existing
``REFERENCE_FUNCTION_COLUMNS`` contract consumed by ``route-b-prepare
--reference-function``, applying an explicit cross-accession mapping table.

Rules (fail-closed unless stated otherwise):

- native identifiers are matched case-preservingly and never casefolded,
  because upstream namespaces differ by case (``LotjaGi1g1v0001500`` vs
  ``Lj2G3V3373110``);
- every crosswalk row must declare a machine-checkable ``mapping_basis``
  (explicit_id_map | pangene_family | synteny | reciprocal_sequence |
  composite); a bare symbol match is not a mapping and is rejected;
- only ``verified`` + ``one_to_one`` mappings may replace the native ID with
  the target ID as ``reference_gene_id``; every other case keeps the NATIVE
  full ID so no false identity is manufactured while the declared state and
  target multiplicity stay visible;
- evidence rows whose class does not assert published/genotyped mutant support
  are skipped and counted, never converted into negative evidence;
- genes without any crosswalk record stay ``not_evaluated``.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

from .route_b_prepare import REFERENCE_FUNCTION_COLUMNS

CROSSWALK_REQUIRED = [
    "source_full_id", "target_full_id", "mapping_basis",
    "relationship", "mapping_state",
]
CROSSWALK_OPTIONAL = [
    "mapping_record", "fam3_family_id", "target_species",
    "sequence_identity", "query_coverage", "subject_coverage",
    "synteny_state", "domain_concordance", "mapping_reason",
]
ALLOWED_MAPPING_BASIS = {
    "explicit_id_map", "pangene_family", "synteny",
    "reciprocal_sequence", "composite",
}
ALLOWED_RELATIONSHIP = {"one_to_one", "one_to_many", "many_to_many", "unresolved"}
ALLOWED_MAPPING_STATE = {"verified", "provisional", "ambiguous", "not_evaluated"}

# Evidence classes that carry a source-published mutant phenotype.  The LORE1 /
# Tnt1 normalizers put this class in the ``source_confidence`` column of
# normalized_functional_evidence.tsv.
MUTANT_EVIDENCE_CLASSES = {
    "strong_causal", "genotype_supported_mutant", "published_phenotype_only",
}

EVIDENCE_REQUIRED = ["native_gene_id", "phenotype_or_function", "source_confidence"]

OUTPUT_NAME = "mutant_reference_functions.tsv"
MANIFEST_NAME = "mutant_reference_bridge_manifest.json"


def _load_crosswalk(paths: list[Path]) -> dict[str, list[dict[str, str]]]:
    """Load explicit mapping rows keyed by exact-case source_full_id."""
    table: dict[str, list[dict[str, str]]] = {}
    for path in paths:
        frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        missing = [column for column in CROSSWALK_REQUIRED if column not in frame.columns]
        if missing:
            raise ValueError(f"cross-accession map {path} is missing columns: {missing}")
        for record in frame.to_dict("records"):
            basis = str(record["mapping_basis"]).strip()
            if basis not in ALLOWED_MAPPING_BASIS:
                raise ValueError(
                    f"cross-accession map {path} row for {record['source_full_id']!r} "
                    f"declares non-machine mapping_basis={basis!r}; symbol-only "
                    "identity transfer is forbidden (issues #478/#525)"
                )
            relationship = str(record["relationship"]).strip()
            if relationship not in ALLOWED_RELATIONSHIP:
                raise ValueError(
                    f"cross-accession map {path} row for {record['source_full_id']!r} "
                    f"has unknown relationship={relationship!r}"
                )
            state = str(record["mapping_state"]).strip()
            if state not in ALLOWED_MAPPING_STATE:
                raise ValueError(
                    f"cross-accession map {path} row for {record['source_full_id']!r} "
                    f"has unknown mapping_state={state!r}"
                )
            key = str(record["source_full_id"])
            table.setdefault(key, []).append({column: str(record.get(column, "")).strip()
                                              for column in CROSSWALK_REQUIRED + CROSSWALK_OPTIONAL})
    return table


def _resolve_mapping(native_id: str, table: dict[str, list[dict[str, str]]]) -> tuple[str, str, str, dict[str, object]]:
    """Return (reference_gene_id, mapping_state, mapping_relationship, detail).

    Only verified one_to_one rows promote the target ID; everything else keeps
    the native full ID and records why.
    """
    rows = table.get(native_id, [])
    if not rows:
        return native_id, "not_evaluated", "", {"targets": []}
    targets = [row["target_full_id"] for row in rows]
    unique_targets = sorted(set(targets))
    basis = "|".join(sorted({row["mapping_basis"] for row in rows}))
    reason = "; ".join(sorted({row["mapping_reason"] for row in rows if row["mapping_reason"]}))
    detail: dict[str, object] = {
        "basis": basis,
        "reason": reason,
        "records": [{k: v for k, v in row.items() if v} for row in rows],
        "targets": unique_targets,
    }
    if len(rows) == 1 and rows[0]["relationship"] == "one_to_one" \
            and rows[0]["mapping_state"] == "verified":
        state = "verified"
    else:
        # Declared ambiguity/provisional state wins; multiplicity alone also
        # degrades to ambiguous even when a row claims verified one-to-one.
        declared = {row["mapping_state"] for row in rows}
        if len(unique_targets) > 1 or len(rows) > 1:
            state = "ambiguous" if "verified" in declared or not declared - {"not_evaluated"} \
                else sorted(declared)[0]
        else:
            state = sorted(declared)[0]
    relationship = "|".join(sorted({row["relationship"] for row in rows}))
    if state == "verified":
        return rows[0]["target_full_id"], state, relationship, detail
    return native_id, state, relationship, detail


def build_mutant_reference_functions(
    evidence_paths: list[Path],
    crosswalk_paths: list[Path],
    out_dir_arg: str | Path,
    *,
    case_id: str = "",
) -> dict[str, object]:
    """Project normalized mutant evidence onto REFERENCE_FUNCTION_COLUMNS."""
    out_dir = Path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    table = _load_crosswalk(crosswalk_paths) if crosswalk_paths else {}
    rows: list[dict[str, object]] = []
    counts = {
        "input_rows": 0, "emitted_rows": 0, "rows_with_mutant_phenotype": 0,
        "skipped_not_mutant_evidence": 0, "unmapped_rows": 0,
        "mapping_state_counts": {"verified": 0, "provisional": 0,
                                 "ambiguous": 0, "not_evaluated": 0},
        "resolved_one_to_one_rows": 0,
    }
    inputs: dict[str, object] = {}
    for path in evidence_paths:
        frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        missing = [column for column in EVIDENCE_REQUIRED if column not in frame.columns]
        if missing:
            raise ValueError(f"mutant functional evidence {path} is missing columns: {missing}")
        inputs[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        for record in frame.to_dict("records"):
            counts["input_rows"] += 1
            evidence_class = str(record.get("source_confidence", "")).strip()
            if evidence_class not in MUTANT_EVIDENCE_CLASSES:
                counts["skipped_not_mutant_evidence"] += 1
                continue
            native_id = str(record["native_gene_id"]).strip()
            phenotype = str(record["phenotype_or_function"]).strip()
            reference_gene_id, mapping_state, mapping_relationship, detail = \
                _resolve_mapping(native_id, table)
            if mapping_state == "not_evaluated":
                counts["unmapped_rows"] += 1
            counts["mapping_state_counts"][mapping_state] = \
                counts["mapping_state_counts"].get(mapping_state, 0) + 1
            if mapping_state == "verified":
                counts["resolved_one_to_one_rows"] += 1
            provenance = {
                "native": {
                    "gene_id": native_id,
                    "species": str(record.get("native_species", "")).strip(),
                    "accession": str(record.get("native_accession", "")).strip(),
                    "assembly": str(record.get("native_assembly", "")).strip(),
                    "annotation": str(record.get("native_annotation", "")).strip(),
                    "namespace": str(record.get("native_gene_id_namespace", "")).strip(),
                    "background": str(record.get("mutant_or_perturbation_background", "")).strip(),
                },
                "mapping": dict(detail),
                "evidence_class": evidence_class,
                "rule": "only verified one_to_one mappings replace the native ID; "
                        "everything else keeps the accession-native identity",
            }
            source_ids = str(record.get("source_ids", "")).strip()
            rows.append({
                "reference_gene_id": reference_gene_id,
                "gene_id": "",
                "function_description": phenotype,
                "pathway": "", "pathway_id": "",
                "mutant_phenotype": phenotype,
                "source_id": source_ids,
                "source_url": "",
                "reference_species": (detail.get("records") or [{}])[0].get("target_species", "")
                if mapping_state == "verified" else str(record.get("native_species", "")).strip(),
                "literature_source_ids": source_ids,
                "phenotype_source_ids": source_ids,
                "verified_source_ids": "",
                "evidence_codes": "",
                "evidence_type": str(record.get("evidence_type", "")).strip(),
                "provenance_tier": evidence_class,
                "provenance_detail": json.dumps(provenance, ensure_ascii=False, sort_keys=True),
                "source_database": "curated_mutant_evidence_bridge",
                "source_release": str(record.get("native_annotation", "")).strip(),
                "retrieved_at": "",
                "corpus_version": "",
                "tissue_specificity": "", "developmental_stage": "",
                "protein_family": "", "go_terms": str(record.get("ontology_ids", "")).strip(),
                "keywords": "",
                "taxon_id": "",
                "gene_symbol": str(record.get("gene_symbol", "")).strip(),
                "family_id": ((detail.get("records") or [{}])[0].get("fam3_family_id", "")
                              if detail.get("records") else ""),
                "orthogroup_id": "",
                "mapping_state": mapping_state,
                "mapping_relationship": mapping_relationship,
                "function_transfer_scope": "",
            })
            counts["emitted_rows"] += 1
            if phenotype:
                counts["rows_with_mutant_phenotype"] += 1
    out_path = out_dir / OUTPUT_NAME
    frame = pd.DataFrame(rows, columns=REFERENCE_FUNCTION_COLUMNS)
    frame.to_csv(out_path, sep="\t", index=False, lineterminator="\n")
    manifest: dict[str, object] = {
        "available": True,
        "path": str(out_path),
        "sha256": hashlib.sha256(out_path.read_bytes()).hexdigest(),
        "case_sensitive_native_match": True,
        "rank_contribution": "none",
        "inputs": {"evidence": inputs,
                   "cross_accession_map": {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                                           for p in crosswalk_paths}},
        **counts,
    }
    if case_id:
        manifest["case_id"] = case_id
    (out_dir / MANIFEST_NAME).write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8")
    return manifest
