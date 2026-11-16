"""Read sealed Direct-C identities as review candidates, never as trait support."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from .core import runtime_path, sha256
from .route_c_claim_binding import load_claim_binding_from_freeze_manifest
from .route_c_direct_members import _verified_mapping_conflicts
from .route_c_direct_projection import _canonical_frame_sha256, validate_projection_evidence_frame
from .route_c_identity_mapping import load_projection_identity_mapping, mapping_links_binding
from .supported_member_union import _load_sealed_input, _verify_route_c_ancestors


def read_direct_c_review_candidates(members_arg: str, manifest_arg: str | None,
                                   *, case_id: str, atlas_frame: pd.DataFrame,
                                   supported_ids: set[str]) -> tuple[pd.DataFrame, dict]:
    from .route_c_provider_bait import provider_edge_links

    path, manifest, manifest_name = _load_sealed_input(
        members_arg, manifest_arg, label="Route C members", case_id=case_id)
    manifest_path = Path(manifest_name)
    _verify_route_c_ancestors(manifest, manifest_path, case_id)

    def ancestor(name: str) -> Path:
        item = manifest["inputs"][name]
        result = runtime_path(str(item["path"]))
        if not result.is_absolute():
            result = manifest_path.parent / result
        if not result.is_file() or sha256(result) != item["sha256"]:
            raise ValueError(f"Route C review ancestor {name} hash mismatch")
        return result

    conflicts, evidence = _verified_mapping_conflicts(str(ancestor("projection_manifest")))
    if evidence is None:
        return pd.DataFrame(), {"route": "c", "state": "not_evaluated_missing_projection_ledger"}
    original_atlas = pd.read_csv(ancestor("atlas"), sep="\t", dtype=str, keep_default_na=False).fillna("")
    if "gene_id" not in atlas_frame or "gene_id" not in original_atlas:
        raise ValueError("Route C review requires an interval atlas")
    atlas_ids = set(atlas_frame["gene_id"].astype(str).str.strip())
    if atlas_ids != set(original_atlas["gene_id"].astype(str).str.strip()):
        raise ValueError("Route C review interval roster mismatch")
    frozen = pd.read_csv(ancestor("frozen"), sep="\t", dtype=str, keep_default_na=False).fillna("")
    records = {str(row["record_id"]): row for row in frozen.to_dict("records")}
    projection_manifest_path = ancestor("projection_manifest")
    identity_mapping = load_projection_identity_mapping(
        json.loads(projection_manifest_path.read_text(encoding="utf-8")),
        projection_manifest_path, records,
    )
    freeze_path = ancestor("freeze_manifest")
    bindings, _ = load_claim_binding_from_freeze_manifest(
        json.loads(freeze_path.read_text(encoding="utf-8")), freeze_path, set(records), records)
    projections = pd.read_csv(ancestor("projection"), sep="\t", dtype=str, keep_default_na=False).fillna("")
    projected_pairs = {(str(row["record_id"]), str(row["focal_gene_id"]))
                       for row in projections.to_dict("records")
                       if row["projection_state"] == "inside_current_qtl"}
    observations = validate_projection_evidence_frame(pd.read_csv(
        runtime_path(str(evidence["path"])), sep="\t", dtype=str, keep_default_na=False).fillna(""))
    interval_fingerprint = f"sha256:{_canonical_frame_sha256(original_atlas, sort_columns=('gene_id',))}"
    if not observations["interval_universe_fingerprint"].eq(interval_fingerprint).all():
        raise ValueError("Route C review projection ledger interval fingerprint mismatch")
    members = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False).fillna("")
    rows = []
    for member in members.to_dict("records"):
        gene, record_id = str(member["candidate_gene_id"]), str(member["record_id"])
        if member["support_status"] != "reference_phenotype_review" or gene in supported_ids:
            continue
        if (record_id, gene) not in projected_pairs or gene not in atlas_ids:
            raise ValueError("Route C review member is not bound to its interval projection")
        record, binding = records.get(record_id, {}), bindings.get(record_id, {})
        # Verify source identity, not focal trait causality. Related traits and
        # unverified perturbation transfer can still be useful review material.
        if record_id in conflicts or not binding:
            continue
        # A cultivar/line detail is not gene identity. An unverified
        # (not_evaluated) cultivar keeps the row review-only; only a
        # contradiction excludes it. Support still requires the complete
        # binding contract. The reported variant stays verified-or-NA.
        if (binding.get("citation_resolution_state") != "verified"
                or binding.get("binding_status") == "contradicted"
                or binding.get("source_organism_state") != "verified_match"
                or binding.get("source_accession_or_cultivar_state") not in {
                    "verified_match", "not_applicable", "not_evaluated"}
                or binding.get("reported_identity_type") not in {"stable_gene_id", "locus", "gene_symbol"}
                or binding.get("reported_identity_verification_state") not in {"verified_exact", "verified_mapped"}
                or binding.get("identity_mapping_basis") in {None, "", "none"}
                or binding.get("reported_variant_state") not in {"verified", "not_applicable"}
                or binding.get("source_scope") not in {"abstract", "full_text", "curated_db_record"}
                or not all(str(binding.get(key, "")).strip() for key in
                           ("reported_identity_value", "quoted_support", "source_record_sha256", "reviewed_by"))
                or record.get("source_verification_state") != "verified"
                or binding.get("citation") != record.get("source_id")):
            continue
        matching = observations[(observations["record_id"] == record_id)
                                & (observations["focal_gene_id"] == gene)
                                & (observations["method_evaluation_state"] == "evaluated_with_hit")]
        identity = str(binding["reported_identity_value"]).strip()
        identity_type = str(binding.get("reported_identity_type"))
        # A provider entry that cites the same paper under the bound name and
        # was the frozen BLAST query links the source identity to this hit.
        provider_linked = any(provider_edge_links(identity_mapping, record, identity_type, identity, item)
                              for item in matching.to_dict("records"))
        if not provider_linked and identity_type == "gene_symbol":
            # Otherwise a symbol reaches a BLAST query only through a reviewed paper edge.
            matching = matching[matching["mapping_method"].eq("blastp")
                                & matching["source_identity_verification_state"].eq("reviewed_paper_sequence_mapping")]
            if not any(mapping_links_binding(identity_mapping, record_id, identity, item, "gene_symbol")
                       for item in matching.to_dict("records")):
                continue
        if not provider_linked and not any((str(item["source_identity_value"]).strip() == identity
                    and str(record.get(str(item["source_identity_field"]), "")).strip() == identity)
                   or mapping_links_binding(identity_mapping, record_id, identity, item,
                                            str(binding["reported_identity_type"]))
                   for item in matching.to_dict("records")):
            continue
        # Preserve the already validated original biological record. A sparse
        # identity-only adapter would otherwise demote a direct trait report
        # during relevance selection. These fields come from one record, not
        # display joins across sources, and do not upgrade its support state.
        rows.append({**record, "candidate_gene_id": gene, "route": "c",
                     "source_reference_gene_id": identity,
                     "projection_state": "inside_current_qtl",
                     "support_status": "unverified_route_review", "rank_contribution": "none",
                     "claim_id": record_id, "source_id": record["source_id"],
                     "verification_status": "review_required",
                     "verification_reason": "sealed source identity and interval projection; focal trait support remains unverified"})
    return pd.DataFrame(rows), {"route": "c", "state": "evaluated", "members": str(path),
                               "members_sha256": sha256(path), "manifest": str(manifest_path),
                               "manifest_sha256": sha256(manifest_path), "review_record_count": len(rows)}
