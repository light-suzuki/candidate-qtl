"""Hash-bound, gene-specific one-to-many reference transport; no trait judgment."""
from __future__ import annotations

import json
from pathlib import Path

from Bio import SeqIO

from .core import runtime_path, sha256
from .reference_name_pairs import record_reference_handles, record_species_matches, reference_handle


def non_unique_insdc_handles(cache_arg: str, handles) -> set[tuple[str, str]]:
    """An available multi-CDS record cannot masquerade as one fetched owner."""
    cache = runtime_path(cache_arg) if cache_arg else None
    if cache is None or not cache.is_dir():
        return set()
    blocked = set()
    for handle in handles:
        if handle[0] != "embl":
            continue
        for suffix, fmt in ((".embl", "embl"), (".gb", "genbank"), (".gbk", "genbank")):
            path = cache / (handle[1] + suffix)
            if not path.is_file():
                continue
            try:
                record = SeqIO.read(path, fmt)
                if (sum(feature.type == "CDS" for feature in record.features) != 1
                        or sum(feature.type == "gene" for feature in record.features) != 1):
                    blocked.add(handle)
            except (ValueError, TypeError):
                blocked.add(handle)
            break
    return blocked


def _records(path: Path) -> dict[str, dict]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    rows = (payload.get("results", [payload] if payload.get("primaryAccession") else list(payload.values()))
            if isinstance(payload, dict) else payload)
    return {row["primaryAccession"]: row for row in rows if isinstance(row, dict) and row.get("primaryAccession")}


def _checked(proof: dict) -> set[str]:
    try:
        return _checked_data(proof)
    except (KeyError, TypeError, AttributeError, IndexError) as error:
        raise ValueError("shared reference proof has malformed structured evidence") from error


def _checked_data(proof: dict) -> set[str]:
    """Recompute membership from saved official records, never from names alone."""
    handle = tuple(proof["handle"])
    if handle[0] not in {"geneid", "embl"}:
        raise ValueError("shared reference proof requires a gene-specific identity")
    descriptors = proof["source_records"]
    if sorted(item["accession"] for item in descriptors) != proof["candidate_roster"]:
        raise ValueError("shared reference proof has an incomplete sibling roster")
    records = {}
    for descriptor in descriptors:
        source = runtime_path(descriptor["path"])
        if not source.is_file() or sha256(source) != descriptor["sha256"]:
            raise ValueError("shared reference proof provider hash mismatch")
        record = _records(source).get(descriptor["accession"])
        if not record or not record_species_matches(record, proof["reference_species"]):
            raise ValueError("shared reference proof provider identity/species mismatch")
        records[descriptor["accession"]] = record
    declared_taxids = [record.get("organism", {}).get("taxonId") for record in records.values()]
    if any(type(taxon) is not int or taxon <= 0 for taxon in declared_taxids) or len(set(declared_taxids)) != 1:
        raise ValueError("shared reference proof requires one explicit taxon")
    taxids = set(declared_taxids)
    owners = {accession for accession, record in records.items() if handle in record_reference_handles(record)}
    if len(owners) < 2 or sorted(owners) != proof["owners"]:
        raise ValueError("shared reference proof owner mismatch")
    if handle[0] == "geneid":
        if proof["basis"] != "official_geneid_one_to_many":
            raise ValueError("shared reference proof basis mismatch")
        if any({value for kind, value in record_reference_handles(records[owner]) if kind == "geneid"} != {handle[1]}
               for owner in owners):
            raise ValueError("shared reference proof contains multiple gene IDs")
        return owners
    if proof["basis"] != "official_insdc_single_gene_alternative_cds":
        raise ValueError("shared reference proof basis mismatch")
    descriptor = proof["insdc_record"]
    source = runtime_path(descriptor["path"])
    if not source.is_file() or sha256(source) != descriptor["sha256"]:
        raise ValueError("shared reference proof INSDC hash mismatch")
    record = SeqIO.read(source, descriptor["format"])
    if not (record.id == handle[1] if "." in handle[1] else record.id.split(".")[0] == handle[1]):
        raise ValueError("shared reference proof INSDC accession mismatch")
    source_taxids = {value.partition(":")[2] for feature in record.features if feature.type == "source"
                    for value in feature.qualifiers.get("db_xref", []) if value.startswith("taxon:")}
    if source_taxids != {str(next(iter(taxids)))}:
        raise ValueError("shared reference proof INSDC taxon mismatch")
    genes = [feature for feature in record.features if feature.type == "gene"]
    cds = [feature for feature in record.features if feature.type == "CDS"]
    if len(genes) != 1 or len(cds) < 2:
        raise ValueError("shared reference proof requires a single gene with alternative CDSs")
    gene = genes[0]
    if gene.location is None or gene.location.strand not in {-1, 1}:
        raise ValueError("shared reference proof lacks an explicit gene location/strand")
    identities = {(key, value) for key in ("gene", "locus_tag") for value in gene.qualifiers.get(key, [])}
    if not identities or any(len(set(gene.qualifiers.get(key, []))) > 1 for key in ("gene", "locus_tag")):
        raise ValueError("shared reference proof lacks an annotated gene identity")
    from .reference_sequence_bridge import _extract_query
    proteins = {}
    for feature in cds:
        if feature.location is None:
            raise ValueError("shared reference proof lacks an explicit CDS location")
        values = feature.qualifiers.get("protein_id", [])
        feature_ids = {(key, value) for key in ("gene", "locus_tag") for value in feature.qualifiers.get(key, [])}
        if (len(values) != 1 or values[0] in proteins or not identities & feature_ids
                or any(value not in identities for value in feature_ids)
                or feature.location.strand != gene.location.strand
                or int(feature.location.start) < int(gene.location.start)
                or int(feature.location.end) > int(gene.location.end)
                or not any("alternativ" in note.casefold() and "splic" in note.casefold()
                           for note in feature.qualifiers.get("note", []))):
            raise ValueError("shared reference proof CDS is not an explicit isoform of this gene")
        translation, _ = _extract_query(record, {"reference_taxid": str(next(iter(taxids))),
            "focal_taxid": str(next(iter(taxids))), "cds_protein_id": values[0]})
        if translation is None:
            raise ValueError("shared reference proof lacks an official CDS translation")
        proteins[values[0]] = translation
    matched = set()
    for owner in owners:
        source_proteins = {prop.get("value") for xref in records[owner].get("uniProtKBCrossReferences", [])
                           if xref.get("database") == "EMBL" and xref.get("id") == handle[1]
                           for prop in xref.get("properties", []) if prop.get("key") == "ProteinId"}
        if not source_proteins or not source_proteins <= proteins.keys():
            raise ValueError("shared reference proof protein endpoint mismatch")
        if any(proteins[protein] != "".join(records[owner].get("sequence", {}).get("value", "").split()).upper().removesuffix("*")
               for protein in source_proteins):
            raise ValueError("shared reference proof official translation/UniProt sequence mismatch")
        matched.update(source_proteins)
    if matched != set(proteins):
        raise ValueError("shared reference proof lacks a sibling CDS/protein owner")
    return owners


def shared_gene_proofs(cache_arg: str, candidates: set[str], species: str,
                       handle_owners: dict[tuple[str, str], set[str]], *, complete: bool) -> dict[tuple[str, str], dict]:
    """Authorize only frozen complete cache records; missing evidence stays unknown."""
    if not complete or not cache_arg:
        return {}
    cache = runtime_path(cache_arg)
    sources = []
    for accession in sorted(candidates):
        path = cache / f"{accession}.json" if cache.is_dir() else cache
        if not path.is_file() or accession not in _records(path):
            return {}
        sources.append({"accession": accession, "path": str(path.resolve()), "sha256": sha256(path)})
    result = {}
    for handle, owners in sorted(handle_owners.items()):
        if len(owners) < 2 or handle[0] not in {"geneid", "embl"}:
            continue
        proof = {"handle": list(handle), "reference_species": species, "candidate_roster": sorted(candidates),
                 "source_records": sources, "owners": sorted(owners),
                 "basis": "official_geneid_one_to_many" if handle[0] == "geneid" else "official_insdc_single_gene_alternative_cds"}
        if handle[0] == "embl":
            if not cache.is_dir():
                continue
            for suffix, fmt in ((".embl", "embl"), (".gb", "genbank"), (".gbk", "genbank")):
                path = cache / (handle[1] + suffix)
                if path.is_file():
                    proof["insdc_record"] = {"path": str(path.resolve()), "sha256": sha256(path), "format": fmt}
                    break
            if "insdc_record" not in proof:
                continue
        try:
            _checked(proof)
        except (ValueError, KeyError, TypeError):
            continue
        result[handle] = proof
    return result


def validate_shared_reference_transports(row) -> None:
    """Final paper consumers recheck the exact original-to-provider mapping."""
    encoded = str(row.get("retrieval_provenance") or "{}")
    if "provider_record_shared_gene_" not in encoded and "paper_shared_reference_proofs" not in encoded:
        return
    provenance = json.loads(encoded)
    if not isinstance(provenance, dict):
        raise ValueError("shared reference transport provenance must be an object")
    proofs = provenance.get("paper_shared_reference_proofs", [])
    receipts = [item for item in provenance.get("paper_reference_rebindings", [])
                if item.get("basis", "").startswith("provider_record_shared_gene_")]
    if not receipts or not proofs:
        raise ValueError("shared reference transport requires matching proofs and receipts")
    try:
        proof_handles = [tuple(proof["handle"]) for proof in proofs]
        receipt_handles = {reference_handle(item["input_reference_gene_id"], allow_untyped=True) for item in receipts}
    except (KeyError, TypeError) as error:
        raise ValueError("shared reference transport has malformed proof/receipt identities") from error
    if len(set(proof_handles)) != len(proof_handles) or set(proof_handles) != receipt_handles:
        raise ValueError("shared reference transport has orphan or duplicate proofs")
    if receipts and (provenance.get("lookup_search_state") != "complete" or provenance.get("lookup_failures")):
        raise ValueError("shared reference transport requires a complete lookup")
    for receipt in receipts:
        handle = reference_handle(receipt["input_reference_gene_id"], allow_untyped=True)
        if not handle or receipt["basis"] != "provider_record_shared_gene_" + handle[0]:
            raise ValueError("shared reference transport receipt basis mismatch")
        matching = [proof for proof in proofs if tuple(proof["handle"]) == handle]
        if len(matching) != 1 or matching[0]["reference_species"] != row.get("reference_species"):
            raise ValueError("shared reference transport lacks its exact proof")
        owners = _checked(matching[0])
        if receipt["reference_gene_id"] != row.get("proposed_stable_id") or row.get("proposed_stable_id") not in owners:
            raise ValueError("shared reference transport provider endpoint mismatch")
