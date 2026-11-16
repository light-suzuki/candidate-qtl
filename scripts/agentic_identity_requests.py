"""Build post-freeze same-species identity requests from a frozen practical stage.

After A freezes, a focal-species paper report whose reference is a UniProt entry
of that same species can be checked against the focal genome: does the reported
gene's own CDS land on an interval gene? The existing post-freeze resolver
(``reference-gene-resolve --sequence-requests --post-freeze``) needs INSDC
accessions, while A records UniProt accessions. This helper reads the sealed
stage's paper evidence table and record cache, follows each record's EMBL
cross-references to INSDC accession/protein_id pairs, and writes the request
table the resolver consumes. It never judges identity, edits a freeze, changes a
Tier or contacts a provider; skipped rows are recorded with a reason.

  python scripts/agentic_identity_requests.py requests --prepared STAGE --sha256 SHA --out DIR
  python scripts/agentic_identity_requests.py protein-map --proteome FASTA --out MAP.tsv
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import sys
from typing import Any

CONTRACT = "post-freeze-identity-requests-1.0.0"
REQUEST_COLUMNS = ("request_id", "source_id", "reported_gene", "reference_accession",
                   "reference_taxid", "focal_taxid", "cds_protein_id")
# UniProt's documented accession pattern.
_UNIPROT = re.compile(r"^(?:UniProt(?:KB)?:)?([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")


def _sha(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _norm(text: str) -> str:
    return " ".join(str(text or "").casefold().split())


def _load_stage(anchor: Path, expected: str) -> dict[str, Any]:
    if _sha(anchor) != expected:
        raise ValueError("stage anchor mismatch")
    stage = json.loads(anchor.read_text(encoding="utf-8"))
    for name, sha in stage.get("files", {}).items():
        path = (anchor.parent / name).resolve()
        if not path.is_relative_to(anchor.parent.resolve()) or _sha(path) != sha:
            raise ValueError(f"stage artifact mismatch: {name}")
    return stage


def _records(record_cache: Path) -> dict[str, dict[str, Any]]:
    files = sorted(record_cache.glob("*.json")) if record_cache.is_dir() else [record_cache]
    found: dict[str, dict[str, Any]] = {}
    for path in files:
        if not path.is_file():
            continue
        payload = json.loads(path.read_text(encoding="utf-8"))
        for entry in (payload if isinstance(payload, list) else [payload]):
            if isinstance(entry, dict) and entry.get("primaryAccession"):
                for accession in [entry["primaryAccession"], *entry.get("secondaryAccessions", [])]:
                    key = str(accession)
                    if key in found and found[key] != entry:
                        raise ValueError(f"conflicting sealed UniProt accession: {key}")
                    found[key] = entry
    return found


def _sealed_input(anchor: Path, stage: dict[str, Any], value: str, *, directory: bool = False) -> Path:
    path = Path(value)
    path = path.resolve() if path.is_absolute() else (anchor.parent / path).resolve()
    if not path.is_relative_to(anchor.parent) or not path.exists():
        raise ValueError("identity input is not inside the sealed stage")
    inputs = sorted(path.glob("*.json")) if directory and path.is_dir() else [path]
    sealed = {(anchor.parent / name).resolve(): sha for name, sha in stage.get("files", {}).items()}
    if not inputs or any(p.resolve() not in sealed or _sha(p) != sealed[p.resolve()] for p in inputs):
        raise ValueError("identity input is not a sealed stage artifact")
    return path


def _embl_pairs(record: dict[str, Any]) -> list[tuple[str, str]]:
    pairs = []
    for xref in record.get("uniProtKBCrossReferences", []):
        if xref.get("database") != "EMBL":
            continue
        props = {p.get("key"): p.get("value") for p in xref.get("properties", [])}
        protein = str(props.get("ProteinId") or "").strip()
        if protein and protein != "-":
            pairs.append((str(xref["id"]), protein))
    return sorted(set(pairs))


def build_requests(prepared: str, expected_sha256: str, out_dir: str) -> dict[str, Any]:
    anchor = Path(prepared).resolve()
    stage = _load_stage(anchor, expected_sha256)
    focal = _norm(stage.get("focal_species", ""))
    if not focal:
        raise ValueError("stage lacks focal_species")
    table = _sealed_input(anchor, stage, stage["a_sources"]["paper_evidence_table"])
    record_cache = _sealed_input(anchor, stage, stage["inputs"]["record_cache"], directory=True)
    out = Path(out_dir)
    if out.exists():
        raise ValueError(f"output already exists: {out}")
    records = _records(record_cache)
    with table.open(encoding="utf-8", newline="") as stream:
        papers = list(csv.DictReader(stream, delimiter="\t"))
    requests, provenance, skipped, seen = [], [], [], set()
    for row in papers:
        if row.get("record_kind") != "paper_report":
            continue
        reason = ""
        species = _norm(row.get("reported_species", ""))
        if species != focal and row.get("focal_species") != "yes":
            reason = "reported_gene_not_focal_species"
        elif row.get("gene_correspondence") in {"different_gene", "homolog"}:
            reason = "reference_is_not_the_reported_gene"
        match = _UNIPROT.match(str(row.get("reference_gene_id", "")).strip())
        record = records.get(match.group(1)) if match else None
        if not reason and match is None:
            reason = "reference_not_uniprot_accession"
        elif not reason and record is None:
            reason = "reference_record_not_in_sealed_cache"
        if not reason and _norm(record.get("organism", {}).get("scientificName", "")) != focal:
            reason = "reference_record_is_not_focal_species"
        pairs = _embl_pairs(record) if not reason else []
        if not reason and not pairs:
            reason = "no_embl_protein_cross_reference"
        if reason:
            skipped.append({"record_id": row.get("record_id", ""), "source_id": row.get("source_id", ""),
                            "reference_gene_id": row.get("reference_gene_id", ""), "reason": reason})
            continue
        taxid = str(record.get("organism", {}).get("taxonId", ""))
        if not taxid.isdigit() or int(taxid) <= 0:
            skipped.append({"record_id": row.get("record_id", ""), "source_id": row.get("source_id", ""),
                            "reference_gene_id": row.get("reference_gene_id", ""), "reason": "reference_taxid_missing"})
            continue
        for accession, protein in pairs:
            key = (row["source_id"], row.get("reported_gene", ""), accession, protein)
            if key in seen:
                continue
            seen.add(key)
            requests.append({
                "request_id": f"idreq-{len(requests) + 1:03d}", "source_id": row["source_id"],
                "reported_gene": row.get("reported_gene", "") or match.group(1),
                "reference_accession": accession, "reference_taxid": taxid, "focal_taxid": taxid,
                "cds_protein_id": protein,
            })
            provenance.append({"request_id": requests[-1]["request_id"], "record_id": row.get("record_id", ""),
                               "reference_gene_id": row["reference_gene_id"],
                               "official_primary_accession": record["primaryAccession"],
                               "accession_binding": "primary" if match.group(1) == record["primaryAccession"] else "official_secondary"})
    out.mkdir(parents=True)
    target = out / "sequence_requests.tsv"
    with target.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=REQUEST_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(requests)
    manifest = {
        "contract": CONTRACT, "case_id": stage.get("case_id", ""), "focal_species": stage.get("focal_species", ""),
        "prepared": str(anchor), "prepared_sha256": expected_sha256,
        "paper_evidence_table": str(table), "paper_evidence_table_sha256": _sha(table),
        "requests": str(target), "requests_sha256": _sha(target), "request_count": len(requests),
        "request_provenance": provenance, "skipped": skipped,
        "scope": ("post_freeze_parent_only request table for the existing resolver; "
                  "not an identity judgement, Tier change or support promotion"),
    }
    (out / "identity_requests_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def protein_map(proteome: str, out_path: str) -> dict[str, Any]:
    """Map each protein FASTA header's first token to its explicit ``gene=`` token."""
    source, out = Path(proteome), Path(out_path)
    if out.exists():
        raise ValueError(f"output already exists: {out}")
    rows, missing = [], 0
    with source.open(encoding="utf-8") as stream:
        for line in stream:
            if not line.startswith(">"):
                continue
            match = re.search(r"(?:^|\s)gene=(\S+)", line)
            if not match:
                missing += 1
                continue
            rows.append({"protein_id": line[1:].split()[0], "gene_id": match.group(1)})
    if missing:
        raise ValueError(f"{missing} FASTA headers lack an explicit gene= token; supply a curated map")
    if len(rows) != len({row["protein_id"] for row in rows}):
        raise ValueError("duplicate protein IDs in proteome")
    with out.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=("protein_id", "gene_id"), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return {"proteome": str(source.resolve()), "proteome_sha256": _sha(source), "map": str(out.resolve()),
            "map_sha256": _sha(out), "protein_count": len(rows),
            "gene_count": len({row["gene_id"] for row in rows}),
            "derivation": "FASTA header first token -> explicit gene= token; no inference"}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("requests", help="build sequence_requests.tsv from a sealed practical stage")
    p.add_argument("--prepared", required=True)
    p.add_argument("--sha256", required=True)
    p.add_argument("--out", required=True)
    p = sub.add_parser("protein-map", help="derive protein_id -> gene_id from explicit gene= FASTA headers")
    p.add_argument("--proteome", required=True)
    p.add_argument("--out", required=True)
    args = parser.parse_args(argv)
    if args.command == "requests":
        result = build_requests(args.prepared, args.sha256, args.out)
    else:
        result = protein_map(args.proteome, args.out)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
