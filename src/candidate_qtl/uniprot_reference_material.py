"""Lossless investigation material from local UniProt JSON; never verified claims."""
from __future__ import annotations

import json
from pathlib import Path

from .core import sha256


def _join(values) -> str:
    return "; ".join(dict.fromkeys(str(v).strip() for v in values if str(v).strip()))


def reference_material(record: dict, *, path: Path) -> dict[str, str]:
    accession = str(record.get("primaryAccession", "")).strip()
    if not accession:
        raise ValueError(f"UniProt record lacks primaryAccession: {path}")
    comments = record.get("comments", [])
    def text(kind):
        return _join(t.get("value", "") for c in comments if c.get("commentType") == kind for t in c.get("texts", []))
    def citations(kind):
        return _join("PMID:" + e["id"] for c in comments if c.get("commentType") == kind
                     for t in c.get("texts", []) for e in t.get("evidences", [])
                     if e.get("source") == "PubMed" and e.get("id"))
    description = record.get("proteinDescription", {})
    names = [description.get("recommendedName", {})] + description.get("alternativeNames", []) + description.get("submissionNames", [])
    protein_name = _join(n.get("fullName", {}).get("value", "") for n in names)
    xrefs = record.get("uniProtKBCrossReferences", [])
    domains = [x for x in xrefs if x.get("database") in {"InterPro", "Pfam", "PROSITE"}]
    go = [x for x in xrefs if x.get("database") == "GO"]
    literature = _join("PMID:" + x["id"] for r in record.get("references", [])
                       for x in r.get("citation", {}).get("citationCrossReferences", [])
                       if x.get("database") == "PubMed" and x.get("id"))
    gene_names = _join(g.get("geneName", {}).get("value", "") for g in record.get("genes", []))
    # Raw record hash preserves detailed evidence without repeating an entire
    # provider record for every interval copy in a 1,000-gene packet.
    return {
        "reference_gene_id": accession, "gene_id": "", "protein_name": protein_name,
        "function_description": text("FUNCTION"), "pathway": text("PATHWAY"),
        "mutant_phenotype": text("DISRUPTION PHENOTYPE"),
        "protein_family": text("SIMILARITY"),
        "domain_ids": _join(x.get("database", "") + ":" + x.get("id", "") for x in domains),
        "domain_description": _join(p.get("value", "") for x in domains for p in x.get("properties", []) if p.get("key") == "EntryName"),
        "go_terms": _join(x.get("id", "") + " " + _join(p.get("value", "") for p in x.get("properties", []) if p.get("key") == "GoTerm") for x in go),
        "keywords": _join(k.get("name", "") for k in record.get("keywords", [])),
        "tissue_specificity": text("TISSUE SPECIFICITY"), "developmental_stage": text("DEVELOPMENTAL STAGE"),
        "reference_species": str(record.get("organism", {}).get("scientificName", "")),
        "taxon_id": str(record.get("organism", {}).get("taxonId", "")), "gene_symbol": gene_names,
        "source_id": "UniProt:" + accession, "source_url": f"https://rest.uniprot.org/uniprotkb/{accession}.json",
        "literature_source_ids": literature, "phenotype_source_ids": citations("DISRUPTION PHENOTYPE"),
        "verified_source_ids": "", "source_database": "UniProtKB",
        "source_release": str(record.get("entryAudit", {}).get("lastAnnotationUpdateDate", "")),
        "evidence_type": "reference_annotation_for_investigation", "provenance_tier": "unverified_annotation",
        "provenance_detail": json.dumps({"record_path": str(path.resolve()), "record_sha256": sha256(path),
                                          "scope": "annotation_and_citation_links_not_verified_claims"}, ensure_ascii=False, sort_keys=True),
        "mapping_state": "exact_uniprot_accession", "mapping_relationship": "reference_accession_record",
        "function_transfer_scope": "reference_record_only",
    }


def read_reference_directory(directory: Path) -> tuple[list[dict[str, str]], dict]:
    rows, files, seen = [], {}, set()
    for path in sorted(directory.glob("*.json")):
        record = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(record, dict) or not record.get("primaryAccession"):
            continue  # failure ledgers are not biological records
        accession = record["primaryAccession"]
        if accession in seen:
            raise ValueError(f"duplicate UniProt accession in record directory: {accession}")
        seen.add(accession)
        rows.append(reference_material(record, path=path))
        files[str(path)] = sha256(path)
    return rows, {"kind": "local_uniprot_json_directory", "files_sha256": files, "record_count": len(rows)}
