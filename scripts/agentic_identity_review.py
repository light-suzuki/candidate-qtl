"""Prepare and transport explicit identity reviews to the existing Wet gate.

No scientific answers are generated here. Old paper answers, partial matches
and competing proteins remain visible; the canonical consumer owns acceptance.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pandas as pd

from candidate_qtl.core import sha256
from candidate_qtl.paper_reports import bound_paper_rows, source_tokens
from candidate_qtl.post_freeze_identity import (
    IDENTITY_CONTRACT, REVIEW_AXES, REVIEW_FIELDS, identity_material_bindings,
    load_identity_assessment, official_identity_record,
)
from candidate_qtl.reference_sequence_bridge import load_reference_sequence_bridge
from candidate_qtl.source_verify import _cache_key, _load_cache, classify_citation
from candidate_qtl.wet_candidate_handoff import _read_review_candidates

if __package__:
    from .agentic_identity_requests import _load_stage, _sealed_input
else:
    from agentic_identity_requests import _load_stage, _sealed_input

CONTRACT = "agentic-identity-review-1.0.0"


def _write(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def _descriptor(path, **extra):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": sha256(path), **extra}


def _quote(path):
    text = Path(path).read_bytes().decode("utf-8")
    return {"start": 0, "end": len(text), "text": text}


def _claims(path, manifest, case_id, atlas):
    frame, _ = _read_review_candidates(str(path), str(manifest), case_id=case_id,
                                      atlas_frame=pd.read_csv(atlas, sep="\t", dtype=str).fillna(""), supported_ids=set())
    return [claim for row in frame.to_dict("records")
            for claim in json.loads(row["relevance_evidence_json"])]


def prepare(prepared, expected_sha256, bridge_path, members, members_manifest,
            paper_records, out_dir, *, provider, model, run_id, reviewer_role="parent_agent"):
    anchor, out = Path(prepared).resolve(), Path(out_dir).resolve()
    if out.exists():
        raise ValueError("identity review requires a new output generation")
    if reviewer_role not in {"parent_agent", "independent_agent"} or not all(str(v).strip() for v in (provider, model, run_id)):
        raise ValueError("explicit reviewer provenance is required")
    stage = _load_stage(anchor, expected_sha256)
    atlas = _sealed_input(anchor, stage, stage["b0"]["interval_gene_atlas"])
    proteins = _sealed_input(anchor, stage, str(Path(stage["inputs"]["packet_b"]) / "candidate_proteins.fa"))
    papers = _sealed_input(anchor, stage, stage["a_sources"]["mechanisms_with_verified_sources"])
    claims = _claims(members, members_manifest, stage["case_id"], atlas)
    bridge, correspondence, _ = load_reference_sequence_bridge(bridge_path)
    cache = _load_cache(paper_records)
    data = {"contract": IDENTITY_CONTRACT, "case_id": stage["case_id"], "post_freeze": True,
            "support_effect": "none", "atlas": _descriptor(atlas),
            "candidate_proteins": _descriptor(proteins),
            "paper_input": _descriptor(papers, id_column="hypothesis_id"),
            "reference_bridge": _descriptor(bridge_path),
            "review": {"reviewer_role": reviewer_role, "provider": provider, "model": model, "run_id": run_id}}
    items, skipped, seen = [], [], set()
    for claim in claims:
        gene = str(claim.get("candidate_gene_id") or claim.get("gene_id") or "")
        for report in bound_paper_rows(claim):
            key = (claim.get("hypothesis_id"), report.get("source_id"), report.get("reported_gene"), gene)
            if key in seen:
                continue
            seen.add(key)
            reason = ""
            checks = json.loads(report.get("paper_checklist") or "{}")
            matches = correspondence[
                correspondence.reported_gene.eq(report.get("reported_gene", ""))
                & correspondence.source_id.map(lambda value: bool(source_tokens(value) & source_tokens(report.get("source_id"))))
                & correspondence.target_id.eq(gene)]
            if checks.get("focal_species") != "yes" or checks.get("same_gene_as_focal") == "no":
                reason = "original_report_not_focal_or_identity_negative"
            elif matches.empty:
                reason = "no_source_bound_sequence_correspondence_for_existing_member"
            kind, identifier = classify_citation(report.get("source_id", ""))
            record = cache.get(_cache_key(kind, identifier))
            if not reason and (not record or not (record.get("abstractText") or record.get("fullText"))):
                reason = "paper_text_not_available_in_supplied_cache"
            if reason:
                skipped.append({"proposal_id": key[0], "source_id": key[1], "reported_gene": key[2], "candidate_gene_id": gene, "reason": reason})
                continue
            if matches.query_sequence_sha256.nunique() != 1:
                skipped.append({"proposal_id": key[0], "source_id": key[1], "candidate_gene_id": gene,
                                "reason": "multiple_reference_sequences_require_explicit_selection",
                                "request_ids": sorted(set(matches.request_id))})
                continue
            # Several papers/accessions may request the same bytes. Keep all
            # competitors in the view, but only one review per paper/member.
            matches = matches[matches.request_id == sorted(set(matches.request_id))[0]]
            for request_id, group in matches.groupby("request_id", sort=True):
                row = group.iloc[0]
                cds_path, cds_sha = official_identity_record(bridge, row.reference_accession)
                if not cds_path:
                    skipped.append({"request_id": request_id, "reason": "official_record_unavailable_or_conflicting"})
                    continue
                source = out / "papers" / f"paper-{len(items):04d}.json"
                _write(source, record)
                item = {"proposal_id": key[0], "source_id": key[1], "reported_gene": key[2],
                        "reported_species": report["reported_species"], "candidate_gene_id": gene,
                        "request_id": request_id, "sequence_reported_gene": row.reported_gene,
                        **{name: row[name] for name in ("reference_accession", "cds_protein_id", "reference_taxid", "focal_taxid", "query_sequence_sha256")},
                        "decision": "not_evaluated", "reason": "",
                        "review_axes": dict.fromkeys(REVIEW_AXES, "unknown"),
                        "paper_source": _descriptor(source, source_id=key[1]),
                        "paper_quote": _quote(source), "cds_quote": _quote(cds_path)}
                if "focal_taxid" in data and data["focal_taxid"] != row.focal_taxid:
                    raise ValueError("mixed focal taxids in review packet")
                data["focal_taxid"] = row.focal_taxid
                item["material_bindings"] = identity_material_bindings(data, item, cds_sha)
                items.append(item)
    data["assessments"] = items
    _write(out / "assessment_materials.json", data)
    _write(out / "review_template.json", [{key: item[key] for key in REVIEW_FIELDS} for item in items])
    # Keep every supplied competitor, including outside-interval and partial hits.
    correspondence.to_csv(out / "all_correspondence.tsv", sep="\t", index=False, lineterminator="\n")
    (out / "REVIEW.md").write_text(
        "# Post-freeze identity review\n\nRead assessment_materials.json and every linked paper, official CDS, "
        "whole-proteome search and all_correspondence.tsv. Inspect the canonical proteins and frozen paper checklist. "
        "The template contains no assessment. Save actual answers as a separate JSON array; change only decision, reason and review_axes. "
        "Account for paper-to-accession, taxon, locus, whole model and all competing copies. "
        "An exact match alone is not a scientific identity judgement. Partial/ambiguous/unavailable evidence stays unknown/not_evaluated. "
        "Do not change original phenotype answers or use benchmark truth. Compile into the existing identity consumer, then finish in a new generation.\n",
        encoding="utf-8")
    manifest = {"contract": CONTRACT, "case_id": stage["case_id"], "prepared": _descriptor(anchor),
                "members": _descriptor(members), "members_manifest": _descriptor(members_manifest),
                "assessment_materials": str(out / "assessment_materials.json"), "review_count": len(items),
                "skipped": skipped, "state": "awaiting_explicit_identity_review" if items else "no_reviewable_correspondence",
                "tier_effect": "none", "support_effect": "none",
                "files": {str(p.relative_to(out)): sha256(p) for p in out.rglob("*") if p.is_file()}}
    _write(out / "prepare_manifest.json", manifest)
    return manifest


def compile_review(manifest_path, expected_sha256, review_path, out_dir):
    path, out = Path(manifest_path).resolve(), Path(out_dir).resolve()
    if out.exists():
        raise ValueError("identity assessment requires a new output generation")
    manifest = _load_stage(path, expected_sha256)
    if manifest.get("contract") != CONTRACT:
        raise ValueError("unsupported identity review packet")
    materials = _sealed_input(path, manifest, manifest["assessment_materials"])
    data = json.loads(materials.read_text(encoding="utf-8"))
    text = Path(review_path).read_bytes().decode("utf-8")
    reviews = json.loads(text)
    if not isinstance(reviews, list) or len(reviews) != len(data["assessments"]) or not reviews:
        raise ValueError("identity review must answer the exact nonempty packet roster")
    # Preserve exact raw JSON spans instead of rewriting the review as our answer.
    offset = text.index("[") + 1
    decoder = json.JSONDecoder()
    for item, review in zip(data["assessments"], reviews):
        while offset < len(text) and text[offset] in " \r\n\t,":
            offset += 1
        parsed, end = decoder.raw_decode(text, offset)
        if parsed != review or set(review) != set(REVIEW_FIELDS):
            raise ValueError("identity review record/schema mismatch")
        immutable = set(REVIEW_FIELDS) - {"decision", "reason", "review_axes"}
        if any(review[key] != item[key] for key in immutable):
            raise ValueError("identity review changed frozen material or endpoint")
        item.update(review)
        item["review_quote"] = {"start": offset, "end": end, "text": text[offset:end]}
        offset = end
    data["review"].update(_descriptor(review_path))
    claims = _claims(manifest["members"]["path"], manifest["members_manifest"]["path"], data["case_id"], data["atlas"]["path"])
    for descriptor in (manifest["prepared"], manifest["members"], manifest["members_manifest"]):
        if sha256(Path(descriptor["path"])) != descriptor["sha256"]:
            raise ValueError("identity review ancestor changed")
    output = out / "identity_assessment.json"
    _write(output, data)
    _, audit = load_identity_assessment(str(output), case_id=data["case_id"],
                                       atlas_arg=data["atlas"]["path"], canonical_proteins_arg=data["candidate_proteins"]["path"], claims=claims)
    result = {"state": "validated_identity_assessment", "assessment": str(output), **audit}
    _write(out / "compile_manifest.json", result)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare")
    for name in ("prepared", "sha256", "bridge", "members", "members-manifest", "paper-records", "out", "provider", "model", "run-id"):
        prep.add_argument("--" + name, required=True)
    prep.add_argument("--reviewer-role", default="parent_agent", choices=["parent_agent", "independent_agent"])
    comp = sub.add_parser("compile")
    for name in ("manifest", "sha256", "review", "out"):
        comp.add_argument("--" + name, required=True)
    args = parser.parse_args(argv)
    try:
        result = (prepare(args.prepared, args.sha256, args.bridge, args.members, args.members_manifest,
                          args.paper_records, args.out, provider=args.provider, model=args.model, run_id=args.run_id,
                          reviewer_role=args.reviewer_role) if args.command == "prepare" else
                  compile_review(args.manifest, args.sha256, args.review, args.out))
        print(json.dumps(result, ensure_ascii=False))
        return 0
    except (ValueError, OSError) as error:
        print(json.dumps({"error": str(error)}), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
