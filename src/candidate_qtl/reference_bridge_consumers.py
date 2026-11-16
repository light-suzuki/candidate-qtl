"""Views/adapters of measured sequence correspondence, never trait judgments."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from .core import runtime_path, sha256
from .paper_reports import PAPER_FIELD, normalize_paper_reports, source_tokens


def direct_bridge_hits(bridge_arg: str, records: dict[str, dict],
                       identity_mapping: pd.DataFrame | None = None) -> tuple[pd.DataFrame, dict[str, tuple[str, str]]]:
    from .reference_sequence_bridge import load_reference_sequence_bridge
    from .route_c_paper_identity_mapping import paper_mapping_records
    if identity_mapping is not None:
        records = paper_mapping_records(records, identity_mapping)
    _, correspondence, hits = load_reference_sequence_bridge(bridge_arg)
    endpoints = {}
    for request_id, rows in correspondence.groupby("request_id", sort=True):
        record = records.get(str(request_id))
        if record is None:
            raise ValueError("reference bridge unknown Direct C record")
        for row in rows.to_dict("records"):
            if not source_tokens(row["source_id"]) & source_tokens(record.get("source_id", "")):
                raise ValueError("reference bridge source endpoint mismatch")
            sequence = str(record.get("reference_sequence_id", "")).strip()
            # A contradictory submitted accession cannot be ignored just because
            # a symbol matches. Neither display names nor unversioned prefixes
            # can substitute for the accession actually retrieved.
            if not sequence or sequence not in {row["reference_accession"], row.get("cds_protein_id", "")}:
                raise ValueError("reference bridge sequence endpoint mismatch")
            tokens = set(json.loads(row.get("source_gene_names_json") or "[]"))
            choices = [(field, str(record.get(field, "")).strip())
                       for field in ("reported_gene_id", "reported_gene_symbol", "reference_sequence_id")]
            matching = [(field, value) for field, value in choices if value and value in tokens]
            stable = str(record.get("reported_gene_id", "")).strip()
            if stable and tokens and stable not in tokens:
                raise ValueError("reference bridge reported gene endpoint mismatch")
            if matching:
                endpoints[str(request_id)] = matching[0]
            elif (row.get("sequence_state") == "same_species_official_translation"
                  and row.get("query_sequence_sha256")
                  and row.get("reported_gene") in {value for field, value in choices
                                                    if field != "reference_sequence_id" and value}):
                # The exact frozen accession authorizes sequence transport,
                # without declaring the paper label an official CDS alias.
                endpoints[str(request_id)] = ("reference_sequence_id", sequence)
            elif str(row.get("evaluation_state", "")).startswith("not_evaluated"):
                continue
            else:
                raise ValueError("reference bridge lacks a structured frozen identity endpoint")
            if (identity_mapping is not None and "evidence_type" in identity_mapping
                    and str(request_id) in set(identity_mapping["record_id"])):
                endpoints[str(request_id)] = ("reference_sequence_id", sequence)
    # The bridge stores all reference competitors. The existing projector
    # performs the interval-atlas join; no top-hit or scaffold preference here.
    # Keep unexecuted retrieval/search distinct from a completed zero-hit panel.
    absent = correspondence[~correspondence["request_id"].isin(set(hits["panel_id"]))]
    statuses = []
    for request_id, rows in absent.groupby("request_id", sort=True):
        unavailable = rows["evaluation_state"].astype(str).str.startswith("not_evaluated").all()
        statuses.append({"panel_id": request_id, "target_id": "", "eligibility": "not_evaluated" if unavailable else "no_hit",
                         "pident": "", "query_coverage": "", "target_coverage": "",
                         "provider_status": "not_evaluated_reference_bridge" if unavailable else "evaluable"})
    if statuses:
        hits = pd.concat([hits, pd.DataFrame(statuses)], ignore_index=True).fillna("")
    return hits, endpoints


def verify_projection_reference_bridge(manifest: dict, manifest_path: Path) -> set[tuple[str, str]]:
    descriptor = manifest.get("inputs", {}).get("reference_bridge")
    if descriptor is None:
        return set()
    path = runtime_path(str(descriptor.get("path", "")))
    if not path.is_absolute():
        path = manifest_path.parent / path
    if not path.is_file() or sha256(path) != descriptor.get("sha256"):
        raise ValueError("reference bridge manifest hash mismatch")
    frozen_descriptor = manifest["inputs"]["frozen"]
    frozen = runtime_path(str(frozen_descriptor["path"]))
    if not frozen.is_absolute():
        frozen = manifest_path.parent / frozen
    if not frozen.is_file() or sha256(frozen) != manifest.get("frozen_sha256"):
        raise ValueError("reference bridge frozen artifact hash mismatch")
    records = pd.read_csv(frozen, sep="\t", dtype=str, keep_default_na=False).to_dict("records")
    from .route_c_identity_mapping import load_projection_identity_mapping
    by_id = {row["record_id"]: row for row in records}
    mapping = load_projection_identity_mapping(manifest, manifest_path, by_id)
    hits, _ = direct_bridge_hits(str(path), by_id, mapping)
    # An unavailable or ordinary corroborating bridge is not counterevidence
    # against an independently supported exact mapping. Block only the models
    # for which the measured correspondence exposes a concrete model mismatch.
    affected = hits[hits["whole_annotation_model_identity"].isin(
        {"partial_model_correspondence", "multiple_reported_genes_in_model"})]
    return set(zip(affected["panel_id"].astype(str), affected["target_id"].astype(str)))


def attach_paper_sequence_correspondence(result: dict, mechanisms_arg: str, bridge_arg: str,
                                        *, id_column: str) -> dict:
    """Attach a parent-only post-freeze view without changing original reports."""
    from .reference_sequence_bridge import load_reference_sequence_bridge
    _, correspondence, _ = load_reference_sequence_bridge(bridge_arg)
    source = runtime_path(mechanisms_arg)
    records = pd.read_csv(source, sep="\t", dtype=str, keep_default_na=False).to_dict("records")
    paper_rows = [(record[id_column], report) for record in records
                  for report in json.loads(normalize_paper_reports(record.get(PAPER_FIELD)) or "[]")]
    joined = []
    for row in correspondence.to_dict("records"):
        matches = [(owner, report) for owner, report in paper_rows
                   if source_tokens(row["source_id"]) & source_tokens([report["source_id"],
                       *report.get("source_identifier_aliases", "").split(";")])
                   and row["reported_gene"] in {report.get("reported_gene"), report.get("reported_gene_id")}]
        if not matches:
            raise ValueError("reference bridge has no matching original paper/gene")
        for owner, report in matches:
            tokens = set(json.loads(row.get("source_gene_names_json") or "[]"))
            structured = bool(tokens & {report.get("reported_gene"), report.get("reported_gene_id")})
            joined.append({**row, "proposal_id": owner,
                           "paper_reported_species": report.get("reported_species", ""),
                           "paper_to_sequence_binding_state": "structured_identity_match_taxon_requires_assessment" if structured else "unresolved_reported_identity_to_sequence",
                           "original_same_gene_as_focal": json.loads(report.get("paper_checklist") or "{}").get("same_gene_as_focal", ""),
                           "assessment_state": "post_freeze_parent_review_required", "membership_effect": "none"})
    out = runtime_path(str(result["manifest"])).parent
    path = out / "paper_sequence_correspondence.tsv"
    pd.DataFrame(joined).to_csv(path, sep="\t", index=False, lineterminator="\n")
    descriptor = {"path": str(path), "sha256": sha256(path),
                  "reference_bridge": {"path": str(runtime_path(bridge_arg).resolve()), "sha256": sha256(runtime_path(bridge_arg))},
                  "source_sha256": sha256(source), "membership_effect": "none",
                  "audience": "post_freeze_parent_only",
                  "next_action": "A request's reported_gene is not a verified gene-to-sequence edge. "
                      "Resolve any unmatched official CDS identity and the paper's experimental taxon before using a sequence as same-gene evidence. "
                      "Assess locus identity separately from whole annotation model and paper phenotype. "
                      "Preserve competing copies and checked unknowns. Do not return interval results to blind A, "
                      "rewrite frozen reports, or infer T1/support from BLAST. Record any revised assessment in a new generation."}
    manifest_path = runtime_path(str(result["manifest"]))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["paper_sequence_correspondence"] = descriptor
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return {**result, "paper_sequence_correspondence": descriptor}
