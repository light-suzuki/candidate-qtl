"""Phenotype-free identity metadata and all-gene coverage for Route B (#751).

Headers and database taxonomy are annotation handles, not orthology or trait
support. Conflicting headers are preserved as records rather than last-writer
wins. Coverage reports reference material availability, not causal evidence.
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Iterable, Mapping

UNKNOWN = "unknown"
METADATA_COLUMNS = [
    "hit_accession", "hit_species", "hit_taxon_id", "hit_gene_symbol",
    "hit_protein_name", "hit_description", "hit_metadata_state",
    "hit_metadata_records_json",
]
IDENTITY_COLUMNS = [
    "focal_gene_id", *METADATA_COLUMNS, "raw_reference_gene_id",
    "percent_identity", "query_coverage", "evalue", "bitscore",
    "provider", "database", "identity_policy", "rank_contribution",
]
INVENTORY_COLUMNS = [
    "gene_id", "identity_state", "identity_hit_count", "identity_metadata_state",
    "identity_accessions", "identity_reference_species", "identity_gene_symbols",
    "identity_protein_names", "functional_annotation_state", "literature_material_state",
    "unresolved_state",
]
_MISSING = {"", "unknown", "n/a", "na", "none", "nan", "not_evaluated", "unresolved", "-"}


def usable(value: object) -> bool:
    return str(value or "").strip().casefold() not in _MISSING


def _usable_taxon_id(value: object) -> bool:
    # BLAST staxids=0 denotes absent database taxonomy. Keep the raw value in
    # metadata records, but do not compare it with an explicit header OX value.
    return usable(value) and str(value).strip() != "0"


def annotation_material(value: object) -> bool:
    """Recognize available annotation text, without treating unknown labels as function."""
    placeholders = ("hypothetical protein", "uncharacterized protein", "uncharacterised protein",
                    "protein of unknown function")
    return any(usable(line) and line.strip().casefold() not in {
        "unknown protein", "predicted protein", "uncharacterized", "uncharacterised",
    } and not any(marker in line.casefold() for marker in placeholders)
               for line in str(value or "").splitlines())


def _json(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def header_metadata(title: str, taxon_ids: str = "", species: str = "") -> dict[str, str]:
    """Parse only explicit UniProt/NCBI fields; never expand a species suffix."""
    title = title.strip()
    result = {key: UNKNOWN for key in (
        "hit_species", "hit_taxon_id", "hit_gene_symbol", "hit_protein_name")}
    result["hit_description"] = title if usable(title) else UNKNOWN
    # Retain native fields so disagreements remain inspectable at record level.
    result["database_taxon_ids"] = taxon_ids.strip() or UNKNOWN
    result["database_species"] = species.strip() or UNKNOWN
    for field, tag in (("hit_species", "OS"), ("hit_taxon_id", "OX"), ("hit_gene_symbol", "GN")):
        match = re.search(rf"(?:^|\s){tag}=(.*?)(?=\s[A-Z]{{2}}=|$)", title)
        if match and usable(match.group(1)):
            result[field] = match.group(1).strip()
    # NCBI uses a terminal organism bracket; Ensembl's terminal Source bracket
    # is attribution instead. Keep the latter in hit_description for audit,
    # but do not promote it to species or conflict with database taxonomy.
    bracket = re.search(r"\[([^\[\]]+)\]\s*$", title)
    source_bracket = bool(bracket and re.match(r"\s*Source\s*:", bracket.group(1), re.IGNORECASE))
    if result["hit_species"] == UNKNOWN and bracket and not source_bracket:
        result["hit_species"] = bracket.group(1).strip()
    product = re.split(r"\sOS=", title, maxsplit=1)[0]
    if bracket:
        product = product[:product.rfind("[")].strip()
    # Remove a FASTA subject identifier only for recognized header layouts.
    if " OS=" in title or bracket:
        parts = product.split(None, 1)
        if parts and ("|" in parts[0] or re.fullmatch(r"[A-Z]{1,4}_?\d+(?:\.\d+)?", parts[0])):
            product = parts[1] if len(parts) == 2 else ""
        if usable(product):
            result["hit_protein_name"] = product
    return result


def merge_header_metadata(records: Iterable[tuple[str, str, str]]) -> dict[str, str]:
    """Deterministic, lossless header roster plus conservative scalar summary."""
    parsed = [header_metadata(*record) for record in sorted(set(records))]
    result: dict[str, str] = {}
    conflicts = False
    for field in ("hit_species", "hit_taxon_id", "hit_gene_symbol", "hit_protein_name", "hit_description"):
        field_usable = _usable_taxon_id if field == "hit_taxon_id" else usable
        values = {row[field] for row in parsed if field_usable(row[field])}
        # Database and header values are both retained; a mismatch is not resolved
        # by trusting whichever HSP happened to occur last.
        database_field = {"hit_species": "database_species", "hit_taxon_id": "database_taxon_ids"}.get(field)
        if database_field:
            values.update(row[database_field] for row in parsed if field_usable(row[database_field]))
        if len(values) > 1:
            conflicts = True
        result[field] = next(iter(values)) if len(values) == 1 else UNKNOWN
    fields = ("hit_species", "hit_taxon_id", "hit_gene_symbol", "hit_protein_name", "hit_description")
    present = sum(usable(result[field]) for field in fields)
    result["hit_metadata_state"] = ("conflicting" if conflicts else "complete" if present == len(fields)
                                    else "partial" if present else "missing")
    result["hit_metadata_records_json"] = _json(parsed)
    return result


def _write(path: Path, columns: list[str], rows: Iterable[Mapping[str, Any]]) -> dict[str, str]:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def protein_query_ids(path: Path | None) -> set[str]:
    """IDs are exact FASTA query IDs, never inferred gene/isoform aliases."""
    if path is None:
        return set()
    ids: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                value = line[1:].strip().split()
                if not value or value[0] in ids:
                    raise ValueError("candidate protein FASTA requires non-empty unique query IDs")
                ids.add(value[0])
    return ids


def write_identity_layer(
    atlas_records: Iterable[Mapping[str, Any]], hit_records: Iterable[Mapping[str, Any]],
    provider: Mapping[str, Any], reference_records: Iterable[Mapping[str, Any]],
    out_dir: Path, *, query_ids: set[str],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Write all-gene status and all-hit metadata, independently of literature.

    Supplied orthology/family annotations can supply identity when local BLAST
    was not run. Literature coverage means reference source-roster availability;
    it does not perform or claim source/trait verification.
    """
    atlas = [dict(row) for row in atlas_records]
    genes = [str(row.get("gene_id", "")).strip() for row in atlas]
    if any(not gene for gene in genes) or len(genes) != len(set(genes)):
        raise ValueError("identity inventory requires non-empty unique interval gene IDs")
    gene_set = set(genes)
    hits: dict[str, list[dict[str, Any]]] = {}
    rows: dict[str, dict[str, Any]] = {}
    unmatched_queries: set[str] = set()
    for item in hit_records:
        hit = dict(item)
        gene = str(hit.get("gene_id", "")).strip()
        if gene not in gene_set:
            unmatched_queries.add(gene)
            continue
        if str(hit.get("threshold_pass", "true")).lower() not in {"true", "1"}:
            continue
        row = {"focal_gene_id": gene, **{key: hit.get(key, UNKNOWN) for key in METADATA_COLUMNS},
               **{key: hit.get(key, "") for key in ("raw_reference_gene_id", "percent_identity",
                  "query_coverage", "evalue", "bitscore", "provider", "identity_policy")},
               "hit_accession": hit.get("hit_accession") or hit.get("reference_gene_id", UNKNOWN),
               "database": provider.get("database", UNKNOWN), "rank_contribution": "none"}
        rows[_json(row)] = row
    for row in [rows[key] for key in sorted(rows)]:
        hits.setdefault(row["focal_gene_id"], []).append(row)
    references: dict[str, list[Mapping[str, Any]]] = {}
    for row in reference_records:
        reference = str(row.get("reference_gene_id", "")).strip()
        if reference:
            references.setdefault(reference, []).append(row)
    available = bool(provider.get("available"))
    inventory: list[dict[str, Any]] = []
    for source in sorted(atlas, key=lambda row: str(row["gene_id"])):
        gene = str(source["gene_id"]).strip()
        gene_hits = hits.get(gene, [])
        accessions = {str(row["hit_accession"]) for row in gene_hits if usable(row["hit_accession"])}
        refs = {part.strip() for part in str(source.get("orthology_reference_gene", "")).split(";") if usable(part)}
        material = [row for ref in accessions | refs for row in references.get(ref, [])]
        annotation_fields = ("function_description", "description", "go_terms", "domain_ids", "domain_architecture")
        functional = (any(usable(source.get(key)) for key in annotation_fields)
                      or any(annotation_material(source.get(key)) for key in (
                          "local_annotation_description", "local_annotation_go_terms"))
                      or any(usable(row.get(key)) for row in material for key in annotation_fields)
                      or any(usable(row.get("hit_protein_name")) for row in gene_hits))
        literature = any(usable(row.get(key)) for row in material for key in (
            "source_id", "literature_source_ids", "phenotype_source_ids", "verified_source_ids"))
        annotated_identity = any(usable(source.get(key)) for key in (
            "family_id", "evolutionary_family_id", "orthogroup", "subfamily_id",
            "orthology_reference_gene", "domain_ids", "domain_architecture"))
        if gene_hits:
            state = "identity_hit"
        elif available and gene in query_ids:
            state = "evaluated_no_threshold_hit"
        elif available:
            state = "not_evaluated_no_protein_query"
        else:
            state = str(provider.get("evaluation_state") or "not_evaluated_provider_unavailable")
        metadata_states = {row["hit_metadata_state"] for row in gene_hits}
        row = {"gene_id": gene, "identity_state": state, "identity_hit_count": len(gene_hits),
               "identity_metadata_state": ";".join(sorted(metadata_states)) or "not_evaluated",
               "functional_annotation_state": "available" if functional else "not_available",
               "literature_material_state": "reference_source_roster_present" if literature else "not_available",
               "unresolved_state": "resolved_handle" if gene_hits or annotated_identity or functional else "unresolved"}
        for column, key in (("identity_accessions", "hit_accession"),
                            ("identity_reference_species", "hit_species"),
                            ("identity_gene_symbols", "hit_gene_symbol"),
                            ("identity_protein_names", "hit_protein_name")):
            row[column] = ";".join(sorted({str(hit[key]) for hit in gene_hits if usable(hit[key])}))
        inventory.append(row)
    out_dir.mkdir(parents=True, exist_ok=True)
    hit_output = _write(out_dir / "route_b_identity_hits.tsv", IDENTITY_COLUMNS, (rows[key] for key in sorted(rows)))
    inventory_output = _write(out_dir / "route_b_identity_inventory.tsv", INVENTORY_COLUMNS, inventory)
    total = len(inventory)
    def coverage(count: int) -> dict[str, Any]:
        return {"gene_count": count, "denominator_gene_count": total, "fraction": count / total if total else None}
    meta = {
        "contract": "route-b-identity-1.0.0", "phenotype_free": True,
        "identity_coverage": coverage(sum(row["unresolved_state"] == "resolved_handle" for row in inventory)),
        "blast_identity_coverage": coverage(sum(row["identity_hit_count"] > 0 for row in inventory)),
        "functional_annotation_coverage": coverage(sum(row["functional_annotation_state"] == "available" for row in inventory)),
        "source_backed_literature_coverage": coverage(sum(row["literature_material_state"] == "reference_source_roster_present" for row in inventory)),
        "dark_unresolved_coverage": coverage(sum(row["unresolved_state"] == "unresolved" for row in inventory)),
        "literature_coverage_basis": "reference source-roster availability, not source/trait verification",
        "provider_evaluation_state": provider.get("evaluation_state", "completed" if available else "not_evaluated"),
        "unmatched_hit_query_ids": sorted(unmatched_queries),
        "query_ids_outside_inventory": sorted(query_ids - gene_set),
        "identity_hits": hit_output, "identity_inventory": inventory_output,
        "rank_contribution": "none", "species_quota": None,
        "retention_scope": "all threshold-passing retrieved hits; upstream BLAST retrieval limit remains explicit",
        "relationship_claim": "identity/annotation only; no one-to-one orthology or trait support",
    }
    return inventory, meta
