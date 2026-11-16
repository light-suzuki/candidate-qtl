"""Offline, explicit taxonomy records for reference-species display audits."""
from __future__ import annotations

import hashlib
import json
import re

from .core import runtime_path


def species_key(value: object) -> str:
    return " ".join(str(value or "").split()).casefold()


def load_taxonomy_cache(path_arg: str | None) -> tuple[dict[str, dict], list[dict]]:
    """Read reserved taxonomy-ID files beside the existing UniProt gene cache.

    No network, name search, genus inference or raw-input repair is performed.
    Other gene-cache JSON files retain their existing loader and semantics.
    """
    if not path_arg:
        return {}, []
    path = runtime_path(path_arg)
    files = sorted(path.glob("taxonomy-*.json")) if path.is_dir() else [path]
    records, artifacts = {}, []
    for item in files:
        match = re.fullmatch(r"taxonomy-([1-9][0-9]*)\.json", item.name)
        if not match:
            continue
        raw = item.read_bytes()
        record = json.loads(raw)
        if (not isinstance(record, dict) or type(record.get("taxonId")) is not int
                or str(record["taxonId"]) != match[1]
                or not isinstance(record.get("scientificName"), str) or not record["scientificName"].strip()
                or ("commonName" in record and not isinstance(record["commonName"], str))
                or not isinstance(record.get("synonyms", []), list)
                or any(not isinstance(value, str) or not value.strip() for value in record.get("synonyms", []))):
            raise ValueError(f"invalid or mismatched reference taxonomy cache: {item}")
        records[match[1]] = record
        artifacts.append({"path": str(item), "sha256": hashlib.sha256(raw).hexdigest()})
    return records, artifacts


def resolve_species_display(display: object, expected_species: str, taxon_id: object,
                            records: dict[str, dict]) -> tuple[str, dict | None]:
    """Resolve only a display exactly reproducible from this endpoint's taxon.

    The claim must use the official scientific name, including any strain or
    subspecies. Synonyms are display components, never accepted claim aliases.
    """
    original = species_key(display)
    record = records.get(str(taxon_id).strip())
    if not record or str(record.get("taxonId")) != str(taxon_id).strip():
        return original, None
    scientific = record["scientificName"]
    if species_key(expected_species) != species_key(scientific):
        return original, None
    decorated = scientific
    if record.get("commonName"):
        decorated += " (" + record["commonName"] + ")"
    for synonym in record.get("synonyms", []):
        decorated += " (" + synonym + ")"
    if original != species_key(decorated) or original == species_key(scientific):
        return original, None
    return str(expected_species).strip().casefold(), {"taxon_id": str(record["taxonId"]), "scientific_name": scientific,
        "frozen_display": str(display), "basis": "same_endpoint_taxon_official_scientific_name_and_reproduced_display"}


def load_bound_taxonomy_artifacts(artifacts: list[dict]) -> dict[str, dict]:
    """Recheck only taxonomy files explicitly sealed by the transfer producer."""
    records = {}
    for artifact in artifacts:
        path = runtime_path(artifact["path"])
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != artifact["sha256"]:
            raise ValueError("Route B reference taxonomy cache hash mismatch")
        cached, descriptors = load_taxonomy_cache(str(path))
        if len(descriptors) != 1 or descriptors[0]["sha256"] != artifact["sha256"]:
            raise ValueError("Route B reference taxonomy cache artifact mismatch")
        if any(key in records and records[key] != value for key, value in cached.items()):
            raise ValueError("Route B reference taxonomy cache records conflict")
        records.update(cached)
    return records
