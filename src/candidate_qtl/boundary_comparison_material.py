"""Carry recorded B context to the sealed boundary without making judgments."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


_ALIGNMENT_FIELDS = {
    "hit_accession", "hit_gene_symbol", "hit_description", "hit_taxid", "database",
    "identity", "query_coverage", "target_coverage", "subject_coverage", "evalue", "bitscore",
    "query_start", "query_end", "subject_start", "subject_end", "qstart", "qend", "sstart", "send",
}


def build_comparison_context(bundle_arg: str, case_id: str, candidate_ids: list[str]) -> dict[str, Any]:
    """Project existing records, retaining per-copy bindings and original source pointers."""
    source = Path(bundle_arg).resolve()
    raw = source.read_bytes()
    bundle = json.loads(raw)
    if (bundle.get("stage") != "post_b0_semantic_bundle" or bundle.get("route") != "interval_first"
            or bundle.get("case_id") != case_id):
        raise ValueError("boundary comparison bundle case/stage mismatch")
    phenotype = bundle.get("phenotype")
    if not isinstance(phenotype, str) or not phenotype.strip():
        raise ValueError("boundary comparison bundle requires the recorded phenotype")
    if len(candidate_ids) != len(set(candidate_ids)) or any(not cid for cid in candidate_ids):
        raise ValueError("boundary comparison requires unique candidate IDs")
    materials: dict[str, Any] = {}
    source_descriptor = {"path": str(source), "sha256": hashlib.sha256(raw).hexdigest()}
    wanted = set(candidate_ids)
    seen: set[str] = set()
    unit_ids: set[str] = set()
    for index, unit in enumerate(bundle["units"]):
        uid, members, genes = unit["unit_id"], unit["member_interval_genes"], unit["genes"]
        if (not uid or uid in unit_ids or len(members) != len(set(members))
                or sorted(g["gene_id"] for g in genes) != sorted(members)
                or seen.intersection(members)):
            raise ValueError("boundary comparison bundle has an invalid member roster")
        unit_ids.add(uid)
        seen.update(members)
        bound_cards = []
        for card_index, card in enumerate(unit.get("reference_material_cards", [])):
            gid = card.get("interval_gene_id")
            if not gid and card.get("gene_id") in members:
                gid = card["gene_id"]
            if gid and gid not in members:
                raise ValueError("boundary comparison reference material is outside its unit")
            bound_cards.append((gid, {
                "source_pointer": f"/units/{index}/reference_material_cards/{card_index}",
                "binding_scope": "gene" if gid else "unit_context", "record": card,
            }))
        for gene_index, gene in enumerate(genes):
            cid = gene["gene_id"]
            if cid not in wanted:
                continue
            pointer = f"/units/{index}/genes/{gene_index}"
            # Raw hit metadata remains available in the hash-bound source bundle;
            # all hits keep their coordinates/coverage without repeating provider dumps.
            alignments = [
                {"source_pointer": f"{pointer}/sequence_similarity_hits/{hit_index}",
                 **{key: value for key, value in hit.items() if key in _ALIGNMENT_FIELDS}}
                for hit_index, hit in enumerate(gene.get("sequence_similarity_hits", []))
            ]
            materials[cid] = {
                "route_b_unit_id": uid, "source_pointer": pointer, "source_bundle": source_descriptor,
                "function_and_copy_context": {key: value for key, value in gene.items() if key != "sequence_similarity_hits"},
                "alignment_evidence": alignments,
                "reference_material_evidence": [record for gid, record in bound_cards if not gid or gid == cid],
            }
    return {
        "case_id": case_id, "phenotype": phenotype,
        "source_bundle": source_descriptor,
        "candidate_material": {cid: materials.get(cid, {"material_state": "not_available_in_route_b_bundle"})
                               for cid in candidate_ids},
        "material_scope": "Recorded context only, not verified function or a new assessment. Missing material is not negative evidence. Full hit metadata remains at the source pointer in source_bundle.",
    }
