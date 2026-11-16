"""Attested paper gene-symbol/sequence edges, authenticated for review only."""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

from .core import runtime_path, sha256
from .observation_tiers import is_unknown
from .paper_reports import source_tokens
from .reference_taxonomy import load_taxonomy_cache, species_key
from .route_c_identity_mapping import PAPER_MAPPING_COLUMNS

PAPER_EDGE = "reviewed_paper_symbol_sequence_v1"


def prepare_paper_mapping_review(frozen_arg: str, freeze_manifest_arg: str, out_arg: str) -> dict:
    """Write an unapproved paper-edge template; blank review fields fail closed."""
    from .route_c_claim_binding_producer import _resolve, _source_context
    _, freeze_path, freeze, records, _, _, source_manifest, _, _ = _source_context(frozen_arg, freeze_manifest_arg)
    source_outputs = json.loads(source_manifest.read_text(encoding="utf-8"))["outputs"]
    captured = _resolve(source_outputs.get("source_review_records", ""), source_manifest.parent)
    if not captured.is_file() or sha256(captured) != freeze.get("source_review_records_sha256"):
        raise ValueError("paper identity mapping requires frozen captured paper text")
    rows = []
    for record in records.to_dict("records"):
        if not record.get("reported_gene_symbol"):
            continue
        row = dict.fromkeys(PAPER_MAPPING_COLUMNS, "")
        row.update(record_id=record["record_id"], source_identity_field="reported_gene_symbol",
                   source_identity_value=record["reported_gene_symbol"], target_identity_field="reference_sequence_id",
                   target_identity_value=record["reference_sequence_id"], reference_species=record["reference_species"],
                   source_namespace="paper_reported_gene", mapping_basis="curated_id_map",
                   source_record_path=str(captured.resolve()), source_record_sha256=sha256(captured),
                   evidence_type=PAPER_EDGE, mapping_relation="reported_gene_sequence", rank_contribution="none")
        rows.append(row)
    out = runtime_path(out_arg)
    out.mkdir(parents=True, exist_ok=True)
    target = out / "paper_identity_mapping_review.tsv"
    pd.DataFrame(rows, columns=PAPER_MAPPING_COLUMNS).to_csv(target, sep="\t", index=False, lineterminator="\n")
    receipt = {"review": str(target), "review_sha256": sha256(target), "state": "not_reviewed",
               "frozen_sha256": sha256(runtime_path(frozen_arg)), "freeze_manifest": str(freeze_path),
               "freeze_manifest_sha256": sha256(freeze_path), "membership_effect": "none",
               "instructions": "Read the frozen paper and record an explicit symbol-to-accession relation, verbatim quote, "
                   "citation, scope, locator, named review, curation_method and reason. Preserve unknowns; blanks cannot be consumed. "
                   "The quote's endpoint occurrence alone does not establish the relation. Do not rewrite frozen records."}
    (out / "paper_identity_mapping_review_manifest.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return receipt


def _literal(text: str, value: str) -> bool:
    return bool(re.search(r"(?<![A-Za-z0-9_.])" + re.escape(value) + r"(?![A-Za-z0-9_]|\.[0-9])", text))


def validate_paper_mapping(frame: pd.DataFrame, path: Path, records: dict[str, dict],
                           freeze_path: Path | None, bridge_arg: str | None) -> pd.DataFrame:
    """Check frozen source bytes and explicit review, never infer a relation."""
    from .route_c_claim_binding import load_claim_binding_from_freeze_manifest
    from .route_c_claim_binding_producer import (
        REVIEWER_CURATION_METHODS, _captured_sources, _citation_key, _normalise_space,
        _resolve, _review_available_scope, _scope_is_available, _source_context,
        _source_row_for_review, _split_citations, _validation_rows,
    )
    from .reference_sequence_bridge import load_reference_sequence_bridge

    frame = frame[list(PAPER_MAPPING_COLUMNS)].apply(lambda col: col.astype(str).str.strip())
    required = [column for column in PAPER_MAPPING_COLUMNS if not column.startswith("taxonomy_record_")]
    if frame.empty or frame[required].eq("").any().any():
        raise ValueError("paper identity mapping requires non-empty fields")
    if any(is_unknown(value) for column in ("reviewed_by", "mapping_reason", "locator") for value in frame[column]):
        raise ValueError("paper identity mapping requires a named reviewer and assessed relation")
    if freeze_path is None or not bridge_arg:
        raise ValueError("paper identity mapping requires frozen source context and official sequence bridge")
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    frozen_path = _resolve((freeze.get("outputs") or {}).get("frozen", ""), freeze_path.parent)
    _, _, _, frozen, validation, _, source_manifest, _, _ = _source_context(str(frozen_path), str(freeze_path))
    frozen_records = {row["record_id"]: row for row in frozen.to_dict("records")}
    if frozen_records != records:
        raise ValueError("paper identity mapping frozen roster mismatch")
    captured = _captured_sources(freeze, source_manifest)
    if captured is None:
        raise ValueError("paper identity mapping requires frozen captured paper text")
    source_outputs = json.loads(source_manifest.read_text(encoding="utf-8"))["outputs"]
    captured_path = _resolve(source_outputs["source_review_records"], source_manifest.parent).resolve()
    bindings, _ = load_claim_binding_from_freeze_manifest(freeze, freeze_path, set(records), records)
    _, correspondence, _ = load_reference_sequence_bridge(bridge_arg)
    for row in frame.to_dict("records"):
        record = records.get(row["record_id"])
        if not record:
            raise ValueError("paper identity mapping unknown record")
        if (record.get("source_verification_state") != "verified"
                or record.get("trait_relation_state") not in {"direct_target_trait", "direct_related_trait"}
                or record.get("evidence_directness") not in {"D0_causal_perturbation_or_rescue",
                    "D1_mapped_gene_with_functional_validation", "D2_direct_trait_gene_association"}):
            raise ValueError("paper identity mapping cannot resolve source or trait eligibility")
        if (row["evidence_type"] != PAPER_EDGE or row["source_namespace"] != "paper_reported_gene"
                or row["mapping_basis"] != "curated_id_map"
                or row["mapping_relation"] != "reported_gene_sequence"
                or row["curation_method"] not in REVIEWER_CURATION_METHODS
                or row["rank_contribution"] != "none"):
            raise ValueError("paper identity mapping requires an explicitly attested paper relation")
        if (row["source_identity_field"] != "reported_gene_symbol"
                or row["source_identity_value"] != record.get("reported_gene_symbol")
                or row["target_identity_field"] != "reference_sequence_id"):
            raise ValueError("paper identity mapping source endpoint mismatch")
        target = row["target_identity_value"]
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*[0-9](?:\.[0-9]+)?", target):
            raise ValueError("paper identity mapping requires a public INSDC accession")
        if record.get("reference_sequence_id") and target != record["reference_sequence_id"]:
            raise ValueError("paper identity mapping contradicts frozen sequence endpoint")
        if (row["reference_species"] != record.get("reference_species")
                or record.get("species_distance") not in {"S0_same_cultivar_or_accession", "S1_same_species"}):
            raise ValueError("paper identity mapping species mismatch")
        binding = bindings.get(row["record_id"], {})
        if (binding.get("reported_identity_type") != "gene_symbol"
                or binding.get("reported_identity_value") != row["source_identity_value"]
                or binding.get("reported_identity_verification_state") not in {"verified_exact", "verified_mapped"}
                or binding.get("citation_resolution_state") != "verified"
                or binding.get("binding_status") in {"contradicted", "not_evaluated"}
                or binding.get("source_organism_state") != "verified_match"
                or binding.get("source_accession_or_cultivar_state") not in {"verified_match", "not_applicable"}
                or binding.get("reported_variant_state") not in {"verified", "not_applicable"}
                or _citation_key(binding.get("citation")) != _citation_key(row["citation"])):
            raise ValueError("paper identity mapping lacks verified source-bound symbol identity")
        citations = _split_citations(record.get("source_id", ""))
        if _citation_key(row["citation"]) not in {_citation_key(item) for item in citations}:
            raise ValueError("paper identity mapping citation mismatch")
        checked = _validation_rows(validation, row["record_id"], citations)
        source_row = _source_row_for_review(checked, row["citation"])
        if (source_row.get("extraction_status") != "accepted"
                or source_row.get("citation_verified", "").casefold() not in {"true", "1"}):
            raise ValueError("paper identity mapping requires accepted source extraction")
        if row["source_scope"] not in {"abstract", "full_text"} or not _scope_is_available(
                row["source_scope"], _review_available_scope(source_row, captured)):
            raise ValueError("paper identity mapping scope exceeds frozen retrieval")
        source = runtime_path(row["source_record_path"])
        if not source.is_absolute():
            source = path.parent / source
        if (source.resolve() != captured_path or not source.is_file()
                or sha256(source) != row["source_record_sha256"].removeprefix("sha256:")):
            raise ValueError("paper identity mapping source is not the frozen paper artifact")
        paper = captured.get(_citation_key(row["citation"]), {})
        sections = paper.get("sections", {})
        text = sections.get("fullText" if row["source_scope"] == "full_text" else "abstractText", "")
        quote = _normalise_space(row["quoted_mapping"])
        if not quote or quote not in _normalise_space(text):
            raise ValueError("paper identity mapping quote is not in frozen source scope")
        # Presence checks are necessary mechanical guards. The named review,
        # not token cooccurrence, attests the semantic relation.
        if not _literal(quote, target) or not _literal(quote, row["source_identity_value"]):
            raise ValueError("paper identity mapping quote lacks exact relation endpoints")
        bridge = correspondence[correspondence.request_id.eq(row["record_id"])]
        if bridge.empty or any(
                item["sequence_state"] != "same_species_official_translation"
                or not item["query_sequence_sha256"]
                or item["reference_taxid"] != item["focal_taxid"]
                or item["reported_gene"] != row["source_identity_value"]
                or target not in {item["reference_accession"], item["cds_protein_id"]}
                or not source_tokens(item["source_id"]) & source_tokens(row["citation"])
                for item in bridge.to_dict("records")):
            raise ValueError("paper identity mapping lacks matching official sequence bridge")
        taxonomy_path, taxonomy_sha = row["taxonomy_record_path"], row["taxonomy_record_sha256"]
        if bool(taxonomy_path) != bool(taxonomy_sha):
            raise ValueError("paper identity mapping taxonomy requires path and hash")
        taxonomy = {}
        if taxonomy_path:
            artifact = runtime_path(taxonomy_path)
            if not artifact.is_absolute():
                artifact = path.parent / artifact
            if not artifact.is_file() or sha256(artifact) != taxonomy_sha.removeprefix("sha256:"):
                raise ValueError("paper identity mapping taxonomy hash mismatch")
            taxonomy, descriptors = load_taxonomy_cache(str(artifact))
            if len(descriptors) != 1:
                raise ValueError("paper identity mapping requires one structured taxonomy artifact")
        for item in bridge.to_dict("records"):
            official = species_key(item["source_organism"])
            expected = species_key(row["reference_species"])
            if not taxonomy_path:
                if not official or official != expected or official != species_key(record["focal_species"]):
                    raise ValueError("paper identity mapping official organism differs; bound taxonomy required")
                continue
            taxon = taxonomy.get(item["reference_taxid"], {})
            names = {species_key(name) for name in [taxon.get("scientificName", ""), *taxon.get("synonyms", [])] if name}
            if not taxon or not {official, expected, species_key(record["focal_species"])} <= names:
                raise ValueError("paper identity mapping taxonomy does not bind source and case species to bridge taxid")
    if frame.groupby(["record_id", "source_identity_field", "source_identity_value"])["target_identity_value"].nunique().gt(1).any():
        raise ValueError("paper identity mapping has conflicting endpoints")
    return frame.drop_duplicates().sort_values(list(PAPER_MAPPING_COLUMNS), kind="stable").reset_index(drop=True)


def paper_mapping_records(records: dict[str, dict], mapping: pd.DataFrame) -> dict[str, dict]:
    """Enrich only an effective view; preserve raw frozen researcher records."""
    result = {key: dict(value) for key, value in records.items()}
    if "evidence_type" in mapping:
        for row in mapping[mapping.evidence_type.eq(PAPER_EDGE)].to_dict("records"):
            record = result[row["record_id"]]
            record["reference_sequence_id"] = row["target_identity_value"]
            if (record.get("source_verification_state") == "verified"
                    and record.get("trait_relation_state") in {"direct_target_trait", "direct_related_trait"}
                    and record.get("evidence_directness") in {"D0_causal_perturbation_or_rescue",
                        "D1_mapped_gene_with_functional_validation", "D2_direct_trait_gene_association"}
                    and record.get("gene_identity_state") == "unresolved"):
                record["gene_identity_state"] = "verified_gene"
                if record.get("direct_rescue_eligibility") == "pending_verification":
                    record["direct_rescue_eligibility"] = "eligible_direct_record"
    return result
