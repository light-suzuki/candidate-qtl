"""Query eligibility gate: forbid Pisum-derived homolog queries.

The Cameor V2 known-gene benchmark claims *non-Pisum, literature-informed,
map-based* candidate prioritization: use a QTL interval on Cameor V2, its native
annotation, and functional homologs from other plants (Medicago, Lotus,
Arabidopsis, ...). A gene id differing between assemblies is not proof of
non-leakage: BLASTing the same causal gene's sequence from another Pisum
reference (Cameor v1a, ZW6, JI lines) -- or a projection / renamed record /
answer-backmapped sequence -- into V2 is searching the answer itself.

This gate validates a homolog-query manifest *before* the benchmark runs and
admits only non-Pisum functional homologs with confirmable provenance. It is a
pre-execution filter; ``discovery-benchmark``'s real-answer-source leave-out
remains a separate, additional audit, not a substitute for this gate.
"""

from __future__ import annotations

import hashlib
import json
import re

import pandas as pd

from . import schemas
from .core import file_fingerprint, runtime_path


GATE_VERSION = "2.0.0"
REQUIRED_COLUMNS = {"case_id", "homolog_query_id", "query_species", "sequence_origin", "selection_basis"}
OPTIONAL_COLUMNS = ["homolog_source_id", "query_genus", "query_taxon_id", "source_assembly", "homolog_source_type"]
# Any of these in a genus / species / assembly / id field means the query is a
# Pisum sequence (including Cameor v1a/v2, ZW6, JI lines and renamed records).
PISUM_MARKERS = ("pisum", "p. sativum", "psat", "cameor", "v1a", "zw6", "ji2822", "ji1794", "ji281")
PISUM_SEQUENCE_ORIGINS = frozenset({"pisum_reference_proteome", "pisum_derived_projection"})
ANSWER_DERIVED_BASES = frozenset({"answer_backmapped", "pisum_sequence_based", "answer_derived", "pisum_causal_paper"})
_OS = re.compile(r"OS=([A-Z][a-z]+ [a-z]+)")
_OX = re.compile(r"OX=(\d+)")
_BRACKET = re.compile(r"\[([A-Z][a-z]+ [a-z]+)\]")


def _looks_pisum(*values: str) -> bool:
    blob = " ".join(str(v) for v in values).casefold()
    return any(marker in blob for marker in PISUM_MARKERS)


def _split_answer_terms(value: object) -> set[str]:
    return {token.strip().casefold() for token in re.split(r"[;,|\s]+", str(value or "")) if token.strip()}


def _parse_known_answer_sources(sources_arg: str | list[str] | None) -> set[str]:
    """Extract normalized answer identity / publication / accession tokens from args or files."""
    if not sources_arg:
        return set()
    terms: set[str] = set()
    args_list = sources_arg if isinstance(sources_arg, list) else [sources_arg]
    for item in args_list:
        raw_str = str(item).strip()
        if not raw_str:
            continue
        path = runtime_path(raw_str)
        if path.exists() and path.is_file():
            try:
                table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
                columns_to_scan = [
                    c for c in table.columns
                    if c in {
                        "causal_publication_ids", "known_answer_derived_sources", "known_answer_aliases",
                        "truth_accession", "truth_gene_id", "truth_coordinate", "answer_gene_id",
                        "answer_source_id", "homolog_source_id", "truth_primary_source",
                    }
                ]
                if columns_to_scan:
                    for col in columns_to_scan:
                        for val in table[col]:
                            terms |= _split_answer_terms(val)
                else:
                    for col in table.columns:
                        for val in table[col]:
                            terms |= _split_answer_terms(val)
            except Exception:
                for line in path.read_text(encoding="utf-8").splitlines():
                    terms |= _split_answer_terms(line)
        else:
            terms |= _split_answer_terms(raw_str)
    return {t for t in terms if t}


def _audit_fasta(path_arg: str | None) -> dict[str, dict[str, str]]:
    """Map query id -> {species, genus, sha256} parsed from FASTA headers."""
    if not path_arg:
        return {}
    from Bio import SeqIO

    audited: dict[str, dict[str, str]] = {}
    with runtime_path(path_arg).open(encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            header = record.description
            species = ""
            match = _OS.search(header) or _BRACKET.search(header)
            if match:
                species = match.group(1)
            audited[record.id] = {
                "species": species, "genus": species.split(" ")[0] if species else "",
                "taxon": (_OX.search(header).group(1) if _OX.search(header) else ""),
                "sha256": hashlib.sha256(str(record.seq).upper().encode("ascii", "ignore")).hexdigest(),
            }
    return audited


def query_eligibility(
    manifest_arg: str,
    out_arg: str,
    query_fasta_arg: str | None = None,
    strict: bool = False,
    full_hash: bool = False,
    known_answer_sources_arg: str | list[str] | None = None,
) -> dict[str, object]:
    """Gate a homolog-query manifest to non-Pisum functional homologs only."""
    manifest_path = runtime_path(manifest_arg)
    out = runtime_path(out_arg)
    out.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str).fillna("")
    missing = REQUIRED_COLUMNS - set(manifest.columns)
    if missing:
        raise ValueError(f"query manifest missing columns: {', '.join(sorted(missing))}")
    for column in OPTIONAL_COLUMNS:
        if column not in manifest.columns:
            manifest[column] = ""
    if manifest.empty:
        raise ValueError("query manifest is empty")
    bad_basis = sorted(set(manifest["selection_basis"]) - schemas.QUERY_SELECTION_BASIS - {""})
    if bad_basis:
        raise ValueError(f"unknown selection_basis: {', '.join(bad_basis)}; allowed {sorted(schemas.QUERY_SELECTION_BASIS)}")
    bad_origin = sorted(set(manifest["sequence_origin"]) - schemas.QUERY_SEQUENCE_ORIGIN - {""})
    if bad_origin:
        raise ValueError(f"unknown sequence_origin: {', '.join(bad_origin)}; allowed {sorted(schemas.QUERY_SEQUENCE_ORIGIN)}")
    audited = _audit_fasta(query_fasta_arg)
    known_answer_terms = _parse_known_answer_sources(known_answer_sources_arg)

    rows: list[dict[str, object]] = []
    for record in manifest.to_dict("records"):
        query_id = str(record["homolog_query_id"]).strip()
        species = str(record["query_species"]).strip()
        declared_genus = str(record.get("query_genus", "")).strip() or (species.split(" ")[0] if species else "")
        assembly = str(record.get("source_assembly", "")).strip()
        source_id = str(record.get("homolog_source_id", "")).strip()
        origin = str(record["sequence_origin"]).strip()
        basis = str(record["selection_basis"]).strip()
        audit = audited.get(query_id, {})
        audited_species = audit.get("species", "")
        audited_genus = audit.get("genus", "")
        query_sha = audit.get("sha256", "")
        provenance_method = "fasta_header" if audit and audited_species else ("manifest_self_report" if (species or declared_genus or str(record.get("query_taxon_id", "")).strip()) else "none")

        # Order matters: Pisum detection (incl. FASTA override) beats self-report.
        pisum_by_source = _looks_pisum(declared_genus, species, assembly, source_id, query_id) or origin in PISUM_SEQUENCE_ORIGINS
        pisum_by_audit = bool(audited_genus) and audited_genus.casefold() == "pisum"

        # Check answer laundering / known answer derived sources
        matched_answer_term = ""
        if known_answer_terms:
            probe_tokens = {
                str(val).strip().casefold()
                for val in (query_id, source_id, assembly, basis, origin, record.get("homolog_source_type", ""))
                if str(val).strip()
            }
            # Direct token match
            direct_matches = probe_tokens & known_answer_terms
            if direct_matches:
                matched_answer_term = sorted(direct_matches)[0]
            else:
                # Substring match for publication IDs or complex accession identifiers
                for term in known_answer_terms:
                    if len(term) >= 4 and any(term in p for p in probe_tokens):
                        matched_answer_term = term
                        break

        if pisum_by_audit or pisum_by_source:
            eligibility, reason = "excluded_pisum", ("fasta_taxon_is_pisum" if pisum_by_audit else "pisum_derived_query")
        elif basis in ANSWER_DERIVED_BASES:
            eligibility, reason = "excluded_answer_derived", basis
        elif matched_answer_term:
            eligibility, reason = "excluded_answer_derived", f"matches_known_answer_source:{matched_answer_term}"
        elif provenance_method == "none" or origin == "unknown":
            eligibility, reason = "excluded_unknown_provenance", "no_confirmable_taxon_provenance"
        elif audited_genus and species and audited_genus.casefold() != declared_genus.casefold():
            eligibility, reason = "excluded_unknown_provenance", "declared_vs_audited_species_mismatch"
        elif strict and query_fasta_arg and not query_sha:
            eligibility, reason = "excluded_unknown_provenance", "query_missing_from_fasta"
        else:
            eligibility, reason = "eligible", ""
        rows.append({
            "case_id": str(record["case_id"]).strip(), "homolog_query_id": query_id, "homolog_source_id": source_id,
            "query_species": species, "query_genus": declared_genus, "query_taxon_id": str(record.get("query_taxon_id", "")).strip(),
            "source_assembly": assembly, "selection_basis": basis, "sequence_origin": origin,
            "audited_species": audited_species, "audited_genus": audited_genus, "query_sha256": query_sha,
            "provenance_method": provenance_method, "benchmark_eligibility": eligibility, "exclusion_reason": reason,
        })

    frame = pd.DataFrame(rows, columns=schemas.QUERY_ELIGIBILITY_COLUMNS).sort_values(["case_id", "homolog_query_id"], kind="stable")
    schemas.validate(frame.fillna(""), schemas.QUERY_ELIGIBILITY_SCHEMA)
    audit_path = out / "query_provenance_audit.tsv"
    frame.to_csv(audit_path, sep="\t", index=False, lineterminator="\n")
    eligible = frame[frame["benchmark_eligibility"].eq("eligible")]
    excluded = frame[~frame["benchmark_eligibility"].eq("eligible")]
    eligible.to_csv(out / "eligible_queries.tsv", sep="\t", index=False, lineterminator="\n")
    excluded.to_csv(out / "excluded_queries.tsv", sep="\t", index=False, lineterminator="\n")
    summary = frame.groupby(["benchmark_eligibility", "exclusion_reason"], as_index=False).size().rename(columns={"size": "count"})
    summary.sort_values(["benchmark_eligibility", "exclusion_reason"], kind="stable").to_csv(out / "eligibility_summary.tsv", sep="\t", index=False, lineterminator="\n")

    pisum_present = int(frame["benchmark_eligibility"].eq("excluded_pisum").sum())
    unknown_present = int(frame["benchmark_eligibility"].eq("excluded_unknown_provenance").sum())
    answer_present = int(frame["benchmark_eligibility"].eq("excluded_answer_derived").sum())
    fingerprints_verified = int(frame["query_sha256"].astype(str).str.strip().ne("").sum())

    query_ledger = [
        {
            "case_id": str(r["case_id"]),
            "query_id": str(r["homolog_query_id"]),
            "source_id": str(r["homolog_source_id"]),
            "species": str(r["query_species"]),
            "genus": str(r["query_genus"]),
            "assembly": str(r["source_assembly"]),
            "origin": str(r["sequence_origin"]),
            "selection_basis": str(r["selection_basis"]),
            "sha256": str(r["query_sha256"]),
            "benchmark_eligibility": str(r["benchmark_eligibility"]),
            "exclusion_reason": str(r["exclusion_reason"]),
        }
        for r in frame.to_dict("records")
    ]

    manifest_out = {
        "schema_version": "1.0", "gate_version": GATE_VERSION, "strict": strict,
        "benchmark_claim": "non-Pisum literature-informed map-based prioritization on Cameor V2",
        "query_count": len(frame), "eligible_count": int(len(eligible)), "excluded_count": int(len(excluded)),
        "pisum_queries_present": pisum_present, "pisum_queries_executed": False,
        "unknown_provenance_queries": unknown_present, "answer_derived_queries": answer_present,
        "fingerprints_verified": fingerprints_verified,
        "queries": query_ledger,
        "exclusion_reasons": {str(k): int(v) for k, v in excluded["exclusion_reason"].value_counts().sort_index().items()},
        "taxon_validation_method": "fasta_header_os_ox_or_brackets" if audited else "manifest_self_report_only",
        "query_manifest": file_fingerprint(manifest_path, full_hash),
        "query_fasta": file_fingerprint(runtime_path(query_fasta_arg), full_hash) if query_fasta_arg else {"status": "not_provided"},
        "known_answer_sources_checked": sorted(known_answer_terms) if known_answer_terms else [],
        "guardrails": [
            "a differing gene id is not accepted as proof of non-leakage",
            "any Pisum-genus / Pisum-assembly / Pisum-derived sequence is barred, including renamed records",
            "answer-backmapped or Pisum-sequence-based selection is barred",
            "queries without confirmable non-Pisum taxon provenance are excluded from the main benchmark",
            "non-Pisum homologs are kept as ordinary map-based comparative-function evidence, not causal proof",
        ],
    }
    manifest_json = out / "query_eligibility_manifest.json"
    manifest_json.write_text(json.dumps(manifest_out, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    if strict and len(excluded):
        reasons = ", ".join(f"{r}={c}" for r, c in manifest_out["exclusion_reasons"].items())
        raise ValueError(f"strict query eligibility gate failed: {len(excluded)} ineligible query(ies) ({reasons}). See {audit_path}")
    return {"audit": str(audit_path), "eligible_queries": str(out / "eligible_queries.tsv"), "excluded_queries": str(out / "excluded_queries.tsv"),
            "manifest": str(manifest_json), "eligible_count": int(len(eligible)), "excluded_count": int(len(excluded)),
            "pisum_queries_executed": False, "fingerprints_verified": fingerprints_verified}

