"""Source-backed typed identity edges for optional Direct-C projection (#663).

Mapping review is a curator assertion, not a proof of trait causality. Raw
identities are never rewritten; consumers recheck the mapping and source bytes.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from .core import runtime_path, sha256
from .route_c_claim_binding import IDENTITY_MAPPING_BASES

IDENTITY_MAPPING_CONTRACT = "route-c-source-identity-mapping-1.0.0"
PAPER_MAPPING_CONTRACT = "route-c-paper-source-identity-mapping-1.0.0"
IDENTITY_MAPPING_COLUMNS = (
    "record_id", "source_identity_field", "source_identity_value",
    "target_identity_field", "target_identity_value", "reference_species",
    "source_namespace", "mapping_basis", "source_record_path", "source_record_sha256",
    "quoted_mapping", "reviewed_by", "rank_contribution",
)
PAPER_MAPPING_COLUMNS = (*IDENTITY_MAPPING_COLUMNS, "evidence_type", "citation", "source_scope",
                         "locator", "curation_method", "mapping_relation", "mapping_reason",
                         "taxonomy_record_path", "taxonomy_record_sha256")


def load_identity_mapping(path: Path | None, records: dict[str, dict], *,
                          freeze_manifest_path: Path | None = None,
                          reference_bridge_arg: str | None = None,
                          eligible_hits: Path | None = None) -> pd.DataFrame:
    """Authenticate explicitly reviewed edges; missing mapping is not negative."""
    if path is None:
        return pd.DataFrame(columns=IDENTITY_MAPPING_COLUMNS)
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False).fillna("")
    from .route_c_provider_bait import EDGE_COLUMNS, load_provider_edges
    if set(frame.columns) == set(EDGE_COLUMNS):
        return load_provider_edges(frame, path, records, eligible_hits)
    if set(frame.columns) == set(PAPER_MAPPING_COLUMNS):
        from .route_c_paper_identity_mapping import validate_paper_mapping
        return validate_paper_mapping(frame, path, records, freeze_manifest_path, reference_bridge_arg)
    if set(frame.columns) != set(IDENTITY_MAPPING_COLUMNS):
        raise ValueError("Route C identity mapping columns do not match the contract")
    frame = frame[list(IDENTITY_MAPPING_COLUMNS)].apply(lambda col: col.astype(str).str.strip())
    if frame.eq("").any().any():
        raise ValueError("Route C identity mapping requires non-empty fields")
    for row in frame.to_dict("records"):
        record = records.get(row["record_id"])
        if record is None:
            raise ValueError("Route C identity mapping unknown record")
        source_field, target_field = row["source_identity_field"], row["target_identity_field"]
        if source_field not in {"reported_gene_id", "reported_locus"}:
            raise ValueError("Route C identity mapping source must be a stable gene or locus")
        if row["source_identity_value"] != str(record.get(source_field, "")).strip():
            raise ValueError("Route C identity mapping source endpoint mismatch")
        if row["reference_species"] != str(record.get("reference_species", "")).strip():
            raise ValueError("Route C identity mapping species mismatch")
        if target_field == "reference_sequence_id":
            if row["target_identity_value"] != str(record.get(target_field, "")).strip():
                raise ValueError("Route C identity mapping sequence endpoint mismatch")
        else:
            raise ValueError("Route C identity mapping target must be a reference sequence")
        if row["mapping_basis"] not in IDENTITY_MAPPING_BASES - {"none", "exact_source_text"}:
            raise ValueError("Route C identity mapping requires a reviewed mapping basis")
        # These UniProt cross-reference databases explicitly identify genes.
        # Other xrefs (e.g. GO, PDB, transcript IDs) cannot authorize a gene edge.
        if row["source_namespace"] not in {"TAIR", "GeneID"}:
            raise ValueError("Route C identity mapping unsupported gene namespace")
        if row["rank_contribution"] != "none":
            raise ValueError("Route C identity mapping must not contribute a score")
        source = runtime_path(row["source_record_path"])
        if not source.is_absolute():
            source = path.parent / source
        if not source.is_file() or sha256(source) != row["source_record_sha256"].removeprefix("sha256:"):
            raise ValueError("Route C identity mapping source hash mismatch")
        quote = row["quoted_mapping"]
        if quote not in source.read_text(encoding="utf-8"):
            raise ValueError("Route C identity mapping quote is not in the source")
        # v1 consumes an authoritative structured UniProt cross-reference. ID
        # cooccurrence in a quote, gene name or free text does not establish a
        # gene -> protein relationship. Other provider formats remain missing
        # mapping rather than being guessed from headers or prose.
        payload = json.loads(source.read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            raise ValueError("Route C identity mapping provider record must be an object")
        accession = str(payload.get("primaryAccession", ""))
        if not accession or row["target_identity_value"] not in {accession, f"UniProt:{accession}"}:
            raise ValueError("Route C identity mapping provider sequence endpoint mismatch")
        organism = payload.get("organism", {})
        if not isinstance(organism, dict) or organism.get("scientificName") != row["reference_species"]:
            raise ValueError("Route C identity mapping provider species mismatch")
        xrefs = payload.get("uniProtKBCrossReferences", [])
        if not isinstance(xrefs, list) or not any(isinstance(item, dict)
                   and item.get("database") == row["source_namespace"]
                   and item.get("id") == row["source_identity_value"] for item in xrefs):
            raise ValueError("Route C identity mapping lacks a structured gene cross-reference")
    edges = frame.drop_duplicates()
    # One edge may have multiple corroborating sources, but a conflicting target
    # cannot be resolved by input order or by choosing an arbitrary alias.
    grouped = edges.groupby(["record_id", "source_identity_field", "source_identity_value",
                             "target_identity_field"], sort=False)
    if grouped["target_identity_value"].nunique().gt(1).any():
        raise ValueError("Route C identity mapping has conflicting endpoints")
    return edges.sort_values(list(IDENTITY_MAPPING_COLUMNS), kind="stable").reset_index(drop=True)


def load_projection_identity_mapping(manifest: dict, manifest_path: Path,
                                     records: dict[str, dict]) -> pd.DataFrame:
    descriptor = manifest.get("inputs", {}).get("identity_mapping")
    if descriptor is None:
        return load_identity_mapping(None, records)
    if not isinstance(descriptor, dict):
        raise ValueError("Route C identity mapping descriptor must be an object")
    from .route_c_provider_bait import CONTRACT as PROVIDER_CONTRACT, EDGE_COLUMNS
    if descriptor.get("contract") not in {IDENTITY_MAPPING_CONTRACT, PAPER_MAPPING_CONTRACT, PROVIDER_CONTRACT}:
        raise ValueError("Route C identity mapping unsupported contract")
    path = runtime_path(str(descriptor.get("path", "")))
    if not path.is_absolute():
        path = manifest_path.parent / path
    if not descriptor.get("sha256") or not path.is_file() or sha256(path) != descriptor["sha256"]:
        raise ValueError("Route C identity mapping artifact hash mismatch")
    freeze_descriptor = manifest.get("inputs", {}).get("freeze_manifest", {})
    freeze_path = runtime_path(str(freeze_descriptor.get("path", "")))
    if not freeze_path.is_absolute():
        freeze_path = manifest_path.parent / freeze_path
    bridge_descriptor = manifest.get("inputs", {}).get("reference_bridge")
    bridge_path = None
    if bridge_descriptor:
        bridge_path = runtime_path(str(bridge_descriptor.get("path", "")))
        if not bridge_path.is_absolute():
            bridge_path = manifest_path.parent / bridge_path
        if not bridge_path.is_file() or sha256(bridge_path) != bridge_descriptor.get("sha256"):
            raise ValueError("Route C identity mapping bridge hash mismatch")
    columns = set(pd.read_csv(path, sep="\t", nrows=0).columns)
    paper = columns == set(PAPER_MAPPING_COLUMNS)
    if paper != (descriptor.get("contract") == PAPER_MAPPING_CONTRACT):
        raise ValueError("Route C identity mapping schema and contract disagree")
    if (columns == set(EDGE_COLUMNS)) != (descriptor.get("contract") == PROVIDER_CONTRACT):
        raise ValueError("Route C identity mapping schema and contract disagree")
    hits_path = None
    hits_descriptor = manifest.get("inputs", {}).get("eligible_hits")
    if hits_descriptor:
        hits_path = runtime_path(str(hits_descriptor.get("path", "")))
        if not hits_path.is_absolute():
            hits_path = manifest_path.parent / hits_path
        if not hits_path.is_file() or sha256(hits_path) != hits_descriptor.get("sha256"):
            raise ValueError("Route C identity mapping eligible hits hash mismatch")
    if paper:
        if not freeze_path.is_file() or sha256(freeze_path) != freeze_descriptor.get("sha256"):
            raise ValueError("Route C identity mapping freeze manifest hash mismatch")
    return load_identity_mapping(path, records, freeze_manifest_path=freeze_path,
                                 reference_bridge_arg=str(bridge_path) if bridge_path else None,
                                 eligible_hits=hits_path)


def mapping_links_binding(frame: pd.DataFrame, record_id: str, identity: str,
                          observation: dict, identity_type: str | None = None) -> bool:
    """Match the source-bound gene to the original sequence endpoint on a hit."""
    if "target_identity_field" not in frame:
        return False  # provider bait edges are matched by route_c_provider_bait
    source_field = {"stable_gene_id": "reported_gene_id", "locus": "reported_locus",
                    "gene_symbol": "reported_gene_symbol"}.get(identity_type)
    eligible = frame if identity_type is None else frame[frame["source_identity_field"] == source_field]
    return bool(((eligible["record_id"] == record_id)
                 & (eligible["source_identity_value"] == identity)
                 & (eligible["target_identity_field"] == observation["source_identity_field"])
                 & (eligible["target_identity_value"] == observation["source_identity_value"])).any())
