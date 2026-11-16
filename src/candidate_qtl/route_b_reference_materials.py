"""Read the explicitly manifest-bound Route-B reference corpus as context.

The canonical preparation manifest is a declared sidecar, not a resource
search. Standalone atlases without it retain their existing behavior. Cards
are attached to their own reference identity; they never authorize support.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _subject_endpoint(value: object) -> str:
    """Unwrap only a complete UniProt sp|ACCESSION|entry wrapper."""
    raw = str(value or "").strip()
    fields = raw.split("|")
    if len(fields) == 3 and fields[0] in {"sp", "tr"} and fields[1].strip() and fields[2].strip():
        return fields[1].strip()
    return raw


def _row_reference_endpoints(row: dict[str, Any]) -> dict[str, dict[str, set[str]]]:
    """Collect exact orthology and family-index endpoints with provenance."""
    endpoints: dict[str, dict[str, set[str]]] = {}

    def add(value: object, source: str, raw: object | None = None) -> None:
        endpoint = str(value or "").strip()
        if not endpoint:
            return
        item = endpoints.setdefault(endpoint, {"sources": set(), "raw": set()})
        item["sources"].add(source)
        item["raw"].add(str(raw if raw is not None else value).strip())

    for raw in str(row.get("orthology_reference_gene", "") or "").split(";"):
        add(raw, "orthology")

    # New atlases carry their syntactically normalized family-index endpoints.
    # Reparse raw subjects as a compatibility fallback for older frozen atlases.
    for endpoint in str(row.get("family_index_subject_endpoints", "") or "").split(";"):
        add(endpoint, "family_index", endpoint)
    raw_subjects = [str(row.get("family_index_best_subject", "") or "")]
    raw_subjects.extend(str(row.get("family_index_subjects", "") or "").split(";"))
    for raw in raw_subjects:
        if raw.strip():
            add(_subject_endpoint(raw), "family_index", raw)
    return endpoints


def load_reference_materials(
    atlas_path: Path, atlas_rows: list[dict[str, Any]],
) -> tuple[dict[str, list[dict[str, str]]], dict[str, Any]]:
    """Join every exact reference endpoint to its hash-verified source cards.

    A present but inconsistent sidecar is an error, not permission to guess a
    corpus. No alias inference, name-based orthology, or home-directory scan.
    """
    manifest_path = atlas_path.parent / "route_b_preparation_manifest.json"
    if not manifest_path.exists():
        return {}, {"status": "not_supplied", "reason": "no preparation sidecar"}
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("command") != "route-b-prepare":
        raise ValueError("reference material sidecar is not a Route-B preparation manifest")
    atlas_fingerprint = manifest.get("frozen_output_fingerprints", {}).get("interval_gene_atlas", {})
    if atlas_fingerprint.get("sha256") != _sha(atlas_path):
        raise ValueError("reference material atlas fingerprint mismatch")
    provider = manifest.get("providers", {}).get("reference_function", {})
    if provider.get("available") is not True:
        return {}, {"status": "not_supplied", "reason": "reference corpus unavailable",
                    "manifest_sha256": _sha(manifest_path)}
    corpus_value = manifest.get("outputs", {}).get("reference_function")
    if not corpus_value or corpus_value != provider.get("path"):
        raise ValueError("inconsistent declared reference corpus path")
    # route_b_prepare records runtime paths. No relocation or alternate-file guessing.
    corpus_path = Path(corpus_value)
    if not corpus_path.is_file() or provider.get("sha256") != _sha(corpus_path):
        raise ValueError("reference material corpus missing or fingerprint mismatch")
    by_reference: dict[str, dict[str, dict[str, str]]] = {}
    with corpus_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "reference_gene_id" not in (reader.fieldnames or []):
            raise ValueError("reference material corpus requires reference_gene_id")
        for raw in reader:
            if None in raw or any(value is None for value in raw.values()):
                raise ValueError("malformed reference material row")
            record = {key: value for key, value in raw.items()}
            ref = record["reference_gene_id"].strip()
            if not ref:
                continue
            record["reference_gene_id"] = ref
            encoded = json.dumps(record, ensure_ascii=False, sort_keys=True)
            by_reference.setdefault(ref, {})[encoded] = record
    cards: dict[str, list[dict[str, str]]] = {}
    for row in atlas_rows:
        gene = str(row["gene_id"])
        refs = _row_reference_endpoints(row)
        selected = []
        for ref in sorted(refs):
            for encoded, record in sorted(by_reference.get(ref, {}).items()):
                selected.append({**record, "interval_gene_id": gene,
                                 "reference_gene": ref, "rank_contribution": "none",
                                 "reference_endpoint_sources": ";".join(sorted(refs[ref]["sources"])),
                                 "reference_endpoint_raw_values": ";".join(sorted(refs[ref]["raw"])),
                                 "material_status": "investigation_context_not_verified_support"})
        if selected:
            cards[gene] = selected
    return cards, {
        "status": "loaded", "manifest_sha256": _sha(manifest_path),
        "corpus_sha256": _sha(corpus_path), "corpus_path": str(corpus_path),
        "gene_count": len(cards), "card_count": sum(map(len, cards.values())),
        "policy": "exact endpoint join; source identity kept per card; no support promotion",
    }
