"""Cheap reported-gene metadata audit, never a paper-content or Tier reviewer.

Only explicit public database IDs are retrieved. Name/display differences and
provider gaps remain unconfirmed; this audit never vetoes candidate membership.
"""
from __future__ import annotations

import json
import re
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

from .core import runtime_path
from .reference_name_pairs import reference_handle, record_reference_handles


def display_key(value: object) -> str:
    value = unicodedata.normalize("NFKC", str(value or "")).casefold()
    return " ".join(re.sub(r"\s*\([^()]*\)\s*$", "", value).split())


def _json(url: str) -> dict:
    request = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "candidate-qtl gene-metadata"})
    with urllib.request.urlopen(request, timeout=20) as response:
        return json.load(response)


def uniprot_gene_metadata(record: dict) -> dict:
    names = []
    for gene in record.get("genes", []):
        names.append(gene.get("geneName", {}).get("value", ""))
        names.extend(item.get("value", "") for field in ("synonyms", "orderedLocusNames", "orfNames") for item in gene.get(field, []))
    organism = record.get("organism", {})
    return {"ids": [f"{kind}:{value}" for kind, value in sorted(record_reference_handles(record))],
            "names": sorted({name for name in names if name}),
            "species": [organism.get("scientificName", ""), organism.get("commonName", ""), *organism.get("synonyms", [])]}


def fetch_gene_metadata(kind: str, identifier: str) -> dict:
    """Direct records for UniProt/NCBI; explicit UniProt xrefs for plant DBs.

    An empty cross-reference search is not evidence that the native ID is fake.
    No sequences, phenotype text or interval identifiers are sent.
    """
    quoted = urllib.parse.quote(identifier, safe="")
    if kind == "uniprot":
        try:
            record = _json(f"https://rest.uniprot.org/uniprotkb/{quoted}.json")
        except urllib.error.HTTPError as error:
            if error.code == 404:
                return {"status": "not_found", "provider": "uniprot"}
            raise
        return {**uniprot_gene_metadata(record), "provider": "uniprot"}
    if kind == "geneid":
        payload = _json("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" +
                        urllib.parse.urlencode({"db": "gene", "id": identifier, "retmode": "json"}))
        record = payload.get("result", {}).get(identifier)
        if not isinstance(record, dict) or record.get("error"):
            return {"status": "unconfirmed", "provider": "ncbi_gene", "reason": "record_unavailable"}
        organism = record.get("organism", {})
        return {"ids": [f"geneid:{record.get('uid', '')}"], "provider": "ncbi_gene",
                "names": [record.get("name", ""), record.get("nomenclaturesymbol", ""),
                          *str(record.get("otheraliases", "")).split(", ")],
                "species": [organism.get("scientificname", ""), organism.get("commonname", "")]}
    if kind in {"embl", "refseq"}:
        protein = bool(re.match(r"(?:[NXYW]P_|[A-Z]{3}\d)", identifier))
        query = urllib.parse.urlencode({"db": "protein" if protein else "nuccore", "id": identifier,
                                       "rettype": "gb", "retmode": "xml"})
        request = urllib.request.Request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + query,
                                         headers={"User-Agent": "candidate-qtl gene-metadata"})
        with urllib.request.urlopen(request, timeout=20) as response:
            payload = ET.fromstring(response.read())
        records = payload.findall(".//GBSeq")
        if not records:
            return {"status": "unconfirmed", "provider": "ncbi_sequence", "reason": "record_unavailable"}
        names, ids, species = set(), set(), set()
        for record in records:
            # Do not strip version suffixes to force an identity match.
            for field in ("GBSeq_primary-accession", "GBSeq_accession-version"):
                if record.findtext(field):
                    ids.add(f"{kind}:{record.findtext(field)}")
            species.add(record.findtext("GBSeq_organism", ""))
            for qualifier in record.findall(".//GBQualifier"):
                if qualifier.findtext("GBQualifier_name") in {"gene", "gene_synonym", "locus_tag"}:
                    names.update(re.split(r";\s*", qualifier.findtext("GBQualifier_value", "")))
        return {"ids": sorted(ids), "names": sorted(names), "species": sorted(species),
                "provider": "ncbi_sequence", "entity_scope": "sequence_record_gene_annotations"}
    query = urllib.parse.urlencode({"query": f"xref:{kind}-{identifier}", "format": "json", "size": 25})
    payload = _json("https://rest.uniprot.org/uniprotkb/search?" + query)
    records = [uniprot_gene_metadata(record) for record in payload.get("results", [])
               if (kind, identifier) in record_reference_handles(record)]
    if len(records) != 1 or payload.get("nextPage"):
        return {"status": "unconfirmed", "provider": "uniprot_xref", "reason": "xref_missing_or_ambiguous"}
    return {**records[0], "provider": "uniprot_xref"}


class GeneMetadataChecks:
    """Cache once per unique ID, reusing an optional cache's gene_records map."""
    def __init__(self, cache_arg: str | None, allow_external_api: bool):
        self.records: dict[str, dict] = {}
        self.allow_external_api = allow_external_api
        self.lookup_count = 0
        if cache_arg:
            path = runtime_path(cache_arg)
            for item in sorted(path.glob("*.json")) if path.is_dir() else [path]:
                payload = json.loads(item.read_text(encoding="utf-8"))
                if isinstance(payload, dict) and isinstance(payload.get("gene_records"), dict):
                    self.records.update(payload["gene_records"])

    def seed_record(self, record: dict | None) -> None:
        """Reuse metadata already obtained by Route A's reference retrieval."""
        if isinstance(record, dict):
            for value in record.get("ids", []):
                handle = reference_handle(value, allow_untyped=True)
                if handle:
                    self.records.setdefault(":".join(handle), record)

    def check(self, report: dict) -> dict:
        gene_id = report.get("reported_gene_id") or report.get("reported_gene", "")
        handle = reference_handle(gene_id, allow_untyped=True)
        result = {"reported_gene": report.get("reported_gene", ""), "reported_gene_id": gene_id,
                  "reported_species": report.get("reported_species", ""), "source_id": report.get("source_id", ""),
                  "state": "unconfirmed", "reason": "no_supported_database_id", "membership_effect": "none"}
        if not handle:
            return result
        kind, identifier = handle
        key = f"{kind}:{identifier}"
        if key not in self.records:
            if not self.allow_external_api:
                self.records[key] = {"status": "unconfirmed", "reason": "not_attempted_offline"}
            else:
                self.lookup_count += 1
                try:
                    self.records[key] = fetch_gene_metadata(kind, identifier)
                except (OSError, ValueError, TypeError, KeyError, AttributeError, ET.ParseError) as error:
                    self.records[key] = {"status": "unconfirmed", "reason": "provider_failure", "error": str(error)}
        record = self.records[key]
        result.update(provider=record.get("provider", "local_record_cache"), database_record=record)
        if record.get("status") == "not_found":
            return {**result, "state": "mismatch", "reason": "id_not_found"}
        if record.get("status") == "unconfirmed":
            return {**result, "reason": record.get("reason", "record_unavailable")}
        actual_ids = {reference_handle(value, allow_untyped=True) for value in record.get("ids", [])}
        if handle not in actual_ids:
            return {**result, "state": "mismatch" if actual_ids - {None} else "unconfirmed", "reason": "returned_identifier_mismatch"}
        name = display_key(report.get("reported_gene"))
        named = bool(name and name in {display_key(value) for value in record.get("names", [])})
        # No arbitrary name-length/edit-distance cutoff. Missing aliases or
        # different long names are unresolved metadata, never a hard rejection.
        if named:
            return {**result, "state": "matched", "reason": "id_name_match"}
        # Species labels are saved as context, never compared to decide this
        # lightweight name/ID audit. Do not build another taxonomy resolver.
        return {**result, "state": "id_resolved", "reason": "gene_name_not_confirmed", "id_exists": True,
                "name_match": False}
