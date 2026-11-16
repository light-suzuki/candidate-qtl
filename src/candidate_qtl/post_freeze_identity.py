"""Validate explicit post-freeze locus reviews without promoting support."""
from __future__ import annotations

import io
import hashlib
import json
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from .core import runtime_path, sha256
from .observation_tiers import CHECKLIST_POLICY, POST_FREEZE_IDENTITY_BASIS, POST_FREEZE_IDENTITY_CONTRACT, normalize_checklist
from .paper_reports import PAPER_FIELD, PAPER_POLICY, REPORT_FIELDS, bound_paper_rows, normalize_paper_reports, source_tokens
from .reference_sequence_bridge import _extract_query, load_reference_sequence_bridge

IDENTITY_CONTRACT = POST_FREEZE_IDENTITY_CONTRACT
REVIEW_AXES = ("paper_to_accession", "experimental_taxon", "locus_identity",
               "whole_model_identity", "competitors_reviewed")
REVIEW_FIELDS = ("proposal_id", "source_id", "reported_gene", "reported_species", "request_id",
                "query_sequence_sha256", "candidate_gene_id", "reference_taxid", "focal_taxid",
                "reference_accession", "cds_protein_id", "sequence_reported_gene", "decision", "reason", "review_axes",
                "material_bindings")


def official_identity_record(bridge: dict, accession: str) -> tuple[str, str]:
    """Repeated requests may save identical official bytes under different paths.

    The bridge loader verifies every saved file. Uniqueness concerns the record
    contents, not how many papers requested them; differing contents still fail.
    """
    records = [record for record in bridge.get("retrieval", [])
               if record.get("reference_accession") == accession]
    hashes = {record.get("raw_sha256") for record in records}
    if not records or len(hashes) != 1 or not next(iter(hashes)) or any(not r.get("raw_path") for r in records):
        return "", ""
    return min(record["raw_path"] for record in records), next(iter(hashes))


def identity_material_bindings(data: dict, item: dict, official_record_sha256: str) -> dict:
    """Fingerprint the exact materials and locations read by the raw reviewer.

    The review file cannot include its own hash. Its parsed response instead
    commits to every scientific input, the quoted spans and executor identity.
    Paths may differ across machines; identical bytes remain the same material.
    """
    def location(value):
        if value is None:
            return None
        if (not isinstance(value, dict) or type(value.get("start")) is not int
                or type(value.get("end")) is not int or not isinstance(value.get("text"), str)):
            raise ValueError("identity assessment invalid material quote locator")
        return {"start": value["start"], "end": value["end"],
                "text_sha256": hashlib.sha256(value["text"].encode("utf-8")).hexdigest()}

    return {"case_id": data["case_id"],
            **{key: data[key]["sha256"] for key in ("atlas", "paper_input", "reference_bridge", "candidate_proteins")},
            "paper_source": item["paper_source"]["sha256"], "official_record": official_record_sha256,
            "paper_quote": location(item.get("paper_quote")), "cds_quote": location(item.get("cds_quote")),
            "reviewer": {key: data["review"][key] for key in ("reviewer_role", "provider", "model", "run_id")}}


def _file(descriptor: dict, base: Path) -> Path:
    if not isinstance(descriptor, dict) or not descriptor.get("path") or not descriptor.get("sha256"):
        raise ValueError("identity assessment requires a hash-bound file")
    path = runtime_path(str(descriptor["path"]))
    if not path.is_absolute():
        path = base / path
    path = path.resolve()
    if not path.is_file() or sha256(path) != descriptor["sha256"]:
        raise ValueError("identity assessment file hash mismatch")
    return path


def _quote(path: Path, quote: dict) -> str:
    # Decode the original bytes without newline conversion; offsets are Unicode
    # character offsets, end exclusive, into this exact UTF-8 source.
    text = path.read_bytes().decode("utf-8")
    if not isinstance(quote, dict):
        raise ValueError("identity assessment requires an exact quoted span")
    start, end, excerpt = quote.get("start"), quote.get("end"), quote.get("text")
    if (type(start) is not int or type(end) is not int or not isinstance(excerpt, str)
            or not excerpt.strip() or not 0 <= start < end <= len(text) or text[start:end] != excerpt):
        raise ValueError("identity assessment quoted span mismatch")
    return excerpt


def _members(row: dict) -> set[str]:
    return set(str(row.get("candidate_gene_id") or row.get("gene_id") or row.get("member_interval_genes") or "").split(";")) - {""}


def _paper_source(path: Path, item: dict, report: dict) -> None:
    """Reuse the existing provider ID/title contract on the actual quoted JSON.

    A labelled arbitrary text file is not a source record. The exact span must
    contain a complete abstract/full-text value belonging to the matched record,
    so another result or a title-only quote cannot supply the scientific text.
    """
    from .bibliography import identifier_matches, normalize_title
    from .source_verify import _cache_key, _load_cache, classify_citation
    kind, identifier = classify_citation(item["source_id"])
    record = _load_cache(str(path)).get(_cache_key(kind, identifier))
    if kind not in {"pmid", "pmcid", "doi"} or not record or not identifier_matches(kind, identifier, record):
        raise ValueError("identity assessment saved paper record identifier mismatch")
    if report.get("source_title") and normalize_title(report["source_title"]) != normalize_title(record.get("title")):
        raise ValueError("identity assessment saved paper record title mismatch")
    _quote(path, item.get("paper_quote"))
    quote = item["paper_quote"]
    raw = path.read_bytes().decode("utf-8")
    decoder = json.JSONDecoder()
    for match in re.finditer(r'"(abstractText|fullText)"\s*:\s*', raw):
        value, length = decoder.raw_decode(raw[match.end():])
        if (isinstance(value, str) and value.strip() and value == record.get(match[1])
                and quote["start"] <= match.end() and quote["end"] >= match.end() + length):
            return
    raise ValueError("identity assessment paper quote must contain matched record abstract/full text")


def load_identity_assessment(assessment_arg: str, *, case_id: str, atlas_arg: str, canonical_proteins_arg: str,
                             claims: list[dict]) -> tuple[list[dict], dict]:
    """Derive member-bound observations only from an existing frozen A/B paper.

    This intentionally narrow consumer accepts a reviewed unique, exact complete
    protein model. Similarity, symbols, unavailable searches and partial models
    alone never resolve locus identity. Other outcomes remain explicit audit.
    """
    path = runtime_path(assessment_arg).resolve()
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("contract") != IDENTITY_CONTRACT or not case_id or data.get("case_id") != case_id:
        raise ValueError("identity assessment contract/case mismatch")
    if data.get("post_freeze") is not True or data.get("support_effect") != "none":
        raise ValueError("identity assessment must be post-freeze and review-only")
    atlas_path = _file(data.get("atlas"), path.parent)
    if atlas_path != runtime_path(atlas_arg).resolve():
        raise ValueError("identity assessment interval atlas mismatch")
    atlas = pd.read_csv(atlas_path, sep="\t", dtype=str, keep_default_na=False)
    if "gene_id" not in atlas or atlas.gene_id.eq("").any() or atlas.gene_id.duplicated().any():
        raise ValueError("identity assessment requires a unique interval roster")
    paper_path = _file(data.get("paper_input"), path.parent)
    papers = pd.read_csv(paper_path, sep="\t", dtype=str, keep_default_na=False)
    id_column = data["paper_input"].get("id_column")
    if id_column not in papers or papers[id_column].eq("").any() or papers[id_column].duplicated().any():
        raise ValueError("identity assessment paper proposal roster mismatch")
    if "case_id" in papers and set(papers.case_id) - {"", case_id}:
        raise ValueError("identity assessment paper case mismatch")
    review_path = _file(data.get("review"), path.parent)
    review = data["review"]
    if (review.get("reviewer_role") not in {"parent_agent", "independent_agent"}
            or any(not str(review.get(key, "")).strip() for key in ("provider", "model", "run_id"))):
        raise ValueError("identity assessment lacks explicit agent review provenance")
    bridge_path = _file(data.get("reference_bridge"), path.parent)
    bridge, correspondence, _ = load_reference_sequence_bridge(str(bridge_path))
    canonical_path = _file(data.get("candidate_proteins"), path.parent)
    if not canonical_proteins_arg or canonical_path != runtime_path(canonical_proteins_arg).resolve():
        raise ValueError("identity assessment canonical candidate proteins mismatch")

    def proteins(path):
        result = {}
        for record in SeqIO.parse(str(path), "fasta"):
            sequence = str(record.seq).upper().removesuffix("*")
            if record.id in result or not sequence or not re.fullmatch("[A-Z]+", sequence):
                raise ValueError("identity assessment requires unique valid protein sequences")
            result[record.id] = sequence
        if not result:
            raise ValueError("identity assessment requires nonempty candidate proteins")
        return result

    canonical_sequences = proteins(canonical_path)
    target_sequences = proteins(_file(bridge["inputs"]["proteome"], bridge_path.parent))
    target_map = pd.read_csv(_file(bridge["inputs"]["protein_gene_map"], bridge_path.parent),
                             sep="\t", dtype=str, keep_default_na=False)
    if (not {"protein_id", "gene_id"} <= set(target_map)
            or target_map.protein_id.duplicated().any()
            or target_map[["protein_id", "gene_id"]].eq("").any().any()):
        raise ValueError("identity assessment bridge protein/gene mapping mismatch")
    target_genes = dict(zip(target_map.protein_id, target_map.gene_id))
    if not set(target_sequences) <= set(target_genes):
        raise ValueError("identity assessment whole proteome has unmapped proteins")
    atlas_genes = set(atlas.gene_id)
    sequences_by_gene, proteins_by_digest = {}, {}
    for protein, sequence in target_sequences.items():
        sequences_by_gene.setdefault(target_genes[protein], set()).add(sequence)
        digest = hashlib.sha256(sequence.encode("utf-8")).hexdigest()
        proteins_by_digest.setdefault(digest, set()).add(protein)
    canonical_by_gene = {}
    for endpoint, sequence in canonical_sequences.items():
        gene = endpoint if endpoint in atlas_genes else target_genes.get(endpoint)
        if (gene not in atlas_genes or sequence not in sequences_by_gene.get(gene, set())
                or (endpoint in atlas_genes and endpoint in target_genes and target_genes[endpoint] != gene)):
            raise ValueError("identity assessment canonical roster/protein differs from the whole-proteome mapping")
        canonical_by_gene.setdefault(gene, set()).add(sequence)
    if set(canonical_by_gene) != atlas_genes:
        raise ValueError("identity assessment canonical proteins do not cover the complete interval roster")
    assessments = data.get("assessments")
    if not isinstance(assessments, list) or not assessments:
        raise ValueError("identity assessment requires explicit assessed records")
    outputs, seen = [], set()
    for item in assessments:
        if not isinstance(item, dict) or any(key not in item for key in REVIEW_FIELDS):
            raise ValueError("identity assessment incomplete record")
        key = (item["proposal_id"], item["source_id"], item["reported_gene"], item["candidate_gene_id"])
        if key in seen:
            raise ValueError("identity assessment duplicate paper/member")
        seen.add(key)
        if item["decision"] not in {"same_gene", "unknown", "not_evaluated"} or not str(item["reason"]).strip():
            raise ValueError("identity assessment invalid decision/reason")
        axes = item["review_axes"]
        if (not isinstance(axes, dict) or set(axes) != set(REVIEW_AXES)
                or set(axes.values()) - {"yes", "no", "unknown"}):
            raise ValueError("identity assessment invalid review axes")
        reviewed = json.loads(_quote(review_path, item.get("review_quote")))
        if reviewed != {key: item[key] for key in REVIEW_FIELDS}:
            raise ValueError("identity assessment contradicts raw review")
        gene = item["candidate_gene_id"]
        if gene not in set(atlas.gene_id):
            raise ValueError("identity assessment member outside interval")
        owner = papers[papers[id_column] == item["proposal_id"]]
        if len(owner) != 1:
            raise ValueError("identity assessment unknown paper proposal")
        reports = json.loads(normalize_paper_reports(owner.iloc[0].get(PAPER_FIELD)) or "[]")
        matches = [report for report in reports if report.get("reported_gene") == item["reported_gene"]
                   and report.get("reported_species") == item["reported_species"]
                   and source_tokens(item["source_id"]) & source_tokens([report["source_id"], *report.get("source_identifier_aliases", "").split(";")])]
        if len(matches) != 1 or not matches[0].get("paper_checklist"):
            raise ValueError("identity assessment paper/gene/species binding mismatch")
        report = matches[0]
        paper_source = _file(item.get("paper_source"), path.parent)
        if not source_tokens(item["paper_source"].get("source_id")) & source_tokens(item["source_id"]):
            raise ValueError("identity assessment saved source identifier mismatch")
        _paper_source(paper_source, item, report)
        # Require this original observation on the same already sealed member.
        # The assessment cannot introduce a new paper or expand a proxy's copies.
        applicable = [part for claim in claims if gene in _members(claim)
                      and str(claim.get(id_column, "")) == item["proposal_id"]
                      for part in bound_paper_rows(claim)
                      if part.get("reported_gene") == report["reported_gene"]
                      and part.get("reported_species") == report["reported_species"]
                      and source_tokens(part.get("source_id")) & source_tokens(report["source_id"])
                      and all(part.get(field, "") == report.get(field, "") for field in REPORT_FIELDS)]
        if not applicable:
            raise ValueError("identity assessment is not bound to an existing frozen member/paper")
        rows = correspondence[correspondence.request_id == item["request_id"]]
        if rows.empty:
            raise ValueError("identity assessment unknown sequence request")
        first = rows.iloc[0]
        for field in ("reference_accession", "cds_protein_id", "reference_taxid", "focal_taxid", "query_sequence_sha256"):
            if set(rows[field]) != {item[field]} or (item["decision"] == "same_gene" and not item[field]):
                raise ValueError("identity assessment sequence endpoint mismatch")
        if (str(data.get("focal_taxid", "")) != item["focal_taxid"]
                or not source_tokens(first["source_id"]) & source_tokens(item["source_id"])
                or first["reported_gene"] != item["sequence_reported_gene"]):
            raise ValueError("identity assessment request source/gene/taxon mismatch")
        raw_name, official_digest = official_identity_record(bridge, item["reference_accession"])
        expected_bindings = identity_material_bindings(data, item, official_digest or "")
        if item["material_bindings"] != expected_bindings:
            raise ValueError("identity assessment raw review material bindings mismatch")
        if item["decision"] != "same_gene" and not item["query_sequence_sha256"]:
            # An unavailable source/CDS/search is a recorded unassessed outcome,
            # not a negative gene judgment or a reason to alter the original Tier.
            continue
        if not raw_name:
            raise ValueError("identity assessment requires a unique saved official record")
        cds_path = _file({"path": raw_name, "sha256": official_digest}, bridge_path.parent)
        _quote(cds_path, item.get("cds_quote"))
        record = None
        for fmt in ("embl", "genbank"):
            try:
                record = SeqIO.read(io.StringIO(cds_path.read_bytes().decode("utf-8")), fmt)
                break
            except ValueError:
                continue
        if record is None or record.id != first["resolved_accession"]:
            raise ValueError("identity assessment official accession mismatch")
        _, extracted = _extract_query(record, item)
        if extracted.get("query_sequence_sha256") != item["query_sequence_sha256"]:
            raise ValueError("identity assessment official CDS/taxon/query mismatch")
        if item["decision"] != "same_gene":
            continue
        if any(axes[axis] != "yes" for axis in REVIEW_AXES):
            raise ValueError("same-gene identity requires completed explicit review of every axis")
        exact = rows[(pd.to_numeric(rows.pident, errors="coerce") == 100)
                     & (pd.to_numeric(rows.query_coverage, errors="coerce") == 1)
                     & (pd.to_numeric(rows.target_coverage, errors="coerce") == 1)]
        if (bridge.get("search", {}).get("state") != "evaluated" or item["reference_taxid"] != item["focal_taxid"]
                or set(rows.sequence_state) != {"same_species_official_translation"}
                or set(exact.target_protein_id) == {""} or exact.target_protein_id.nunique() != 1
                or set(exact.target_id) != {gene}
                or set(exact.eligibility) != {"eligible"}
                or set(exact.whole_annotation_model_identity) != {"model_coverage_passes_review_threshold"}):
            raise ValueError("same-gene identity requires a unique exact whole-proteome model; partial/ambiguous/unassessed remains unknown")
        protein_id = exact.iloc[0]["target_protein_id"]
        actual = canonical_by_gene[gene]
        if (target_genes.get(protein_id) != gene or protein_id not in target_sequences
                or target_sequences[protein_id] not in actual
                or hashlib.sha256(target_sequences[protein_id].encode("utf-8")).hexdigest() != item["query_sequence_sha256"]):
            raise ValueError("identity assessment bridge target sequence differs from the canonical candidate protein")
        # Audit the actual sealed search universe as well as its BLAST ledger.
        # Exact matching is an integrity check here, not a replacement for the
        # original standard search or its retained non-exact competitors.
        if (proteins_by_digest.get(item["query_sequence_sha256"]) != {protein_id}
                or {member for member, sequences in canonical_by_gene.items()
                    if target_sequences[protein_id] in sequences} != {gene}):
            raise ValueError("identity assessment actual proteome/canonical proteins contain competing exact models")
        checks = json.loads(report["paper_checklist"])
        if checks["same_gene_as_focal"] == "no" or checks["focal_species"] != "yes":
            raise ValueError("identity assessment cannot override negative identity or unresolved paper species")
        checks["same_gene_as_focal"] = "yes"
        # The reviewed exact match is the member correspondence that paper-side
        # T1 is conditional on (#7). Record it as an explicit projection so the
        # Tier gate consumes it instead of a same-gene answer standing in for it.
        projection = {"basis": POST_FREEZE_IDENTITY_BASIS, "route": "post_freeze_identity",
                      "candidate_gene_id": gene, "member_genes": [gene], "copy_count": 1,
                      "reference_gene_id": report.get("reference_gene_id", ""),
                      "evidence": {"target_protein_id": protein_id, "pident": 100.0,
                                   "query_coverage": 1.0, "target_coverage": 1.0,
                                   "query_sequence_sha256": item["query_sequence_sha256"],
                                   "identity_assessment_sha256": sha256(path)}}
        # This is a new gene-bound observation; the blind report stays verbatim
        # in the audit and in the original claim. No support flags are added.
        outputs.append({**report, "candidate_gene_id": gene, "paper_checklist": normalize_checklist(checks),
                        "homolog_projection_json": json.dumps(projection, sort_keys=True),
                        "tier_classification_policy": CHECKLIST_POLICY, "tier_gene_binding": "matched",
                        "tier_aggregation": PAPER_POLICY,
                        "initial_tier_reason": item["reason"], "post_freeze_identity_contract": IDENTITY_CONTRACT,
                        "identity_assessment_sha256": sha256(path), "identity_assessment": item,
                        "original_paper_report": report, "support_effect": "none"})
    audit = {"path": str(path), "sha256": sha256(path), "contract": IDENTITY_CONTRACT,
             "assessment_count": len(assessments), "resolved_evidence_count": len(outputs), "support_effect": "none"}
    return outputs, audit
