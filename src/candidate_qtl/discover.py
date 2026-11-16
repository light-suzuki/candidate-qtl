from __future__ import annotations

import hashlib
import json
import os
import re
import sys
import threading
import time
import urllib.error
from pathlib import Path

import pandas as pd
import yaml
from Bio import SeqIO

from . import homology_qc
from .core import ExternalDependencyError, _env_timeout, runtime_path, run_tool
from .evidence_io import EVIDENCE_COLUMNS, evidence_row as _row
from .hardening import prioritize_evidence
from .workflow_extensions import candidate_flags, extract_multi


# Words that carry no trait specificity on their own. Broad anatomy/direction
# terms still match (recorded as context) but must not, alone, become scored
# trait evidence; stopwords are dropped entirely.
_STOPWORDS = {"the", "and", "of", "in", "a", "an", "to", "with", "for", "on", "at", "by",
              "increase", "increased", "decrease", "decreased", "gain", "loss", "high", "low",
              "more", "less", "reduced", "enhanced", "altered", "abnormal", "defect", "defective"}
_GENERIC_ANATOMY = {"stem", "leaf", "leaves", "flower", "flowers", "seed", "seeds", "root", "roots",
                    "shoot", "shoots", "organ", "tissue", "plant", "growth", "size", "number", "length",
                    "width", "height", "development", "developmental", "phenotype", "trait", "cell", "cells"}

# Dry ranking defaults for the discovery lane. Reserve and Rescue take the same
# slot count as Top, and the score gate is lenient (1.0) so a weak keyword match
# stays visible; the case-run lane defaults to top_n=3/reserve_n=3/rescue_n=3
# and min_top_score=2.0 unless the case manifest declares otherwise.
DEFAULT_TOP_N = 10
DEFAULT_MIN_TOP_SCORE = 1.0


def _build_terms(phenotype: str, extra: list[str] | None, approved: list[tuple[str, str]] | None = None) -> dict[str, object]:
    """Build a controlled trait vocabulary instead of splitting into all tokens.

    With --approved-terms the vocabulary is the human-approved set (weight_class
    ``context`` becomes broad/unscored). Otherwise the phenotype phrase is kept
    intact, stopwords are dropped, and generic anatomy/direction words are set
    aside as broad (context-only) rather than scored trait terms.
    """
    if approved is not None:
        scored = sorted({term for term, weight_class in approved if weight_class != "context"})
        broad = sorted({term for term, weight_class in approved if weight_class == "context"})
        return {"phrases": [], "scored": scored, "broad": broad, "source": "approved_terms"}
    values = [phenotype, *(extra or [])]
    phrases = sorted({value.strip().casefold() for value in values if value.strip()})
    tokens: set[str] = set()
    for value in values:
        tokens.update(token for token in re.findall(r"[\w-]+", value.casefold()) if len(token) >= 3)
    tokens -= _STOPWORDS
    broad = sorted(tokens & _GENERIC_ANATOMY)
    scored = sorted(tokens - _GENERIC_ANATOMY)
    return {"phrases": phrases, "scored": scored, "broad": broad, "source": "auto_split"}


def _match(text: str, words: list[str]) -> list[str]:
    lowered = text.casefold()
    return [word for word in words if word and re.search(rf"(?<![0-9a-z]){re.escape(word)}(?![0-9a-z])", lowered)]


# Trait vocabulary normalization (issue #400 step 1). Substring and shared-stem
# matching are additive to the exact whole-word boundary match in ``_match`` and
# only fire when no scored term matched exactly. Minimum word lengths guard
# against spurious hits from short shared prefixes ("colo" in colour/color) or
# tiny contained tokens ("red" inside a longer word). These are morphological
# rules only; the emitted row keeps direction=support at low confidence so it
# stays visible but never outranks an exact token/phrase match.
_MIN_SUBSTRING_LEN = 6
_MIN_STEM_LEN = 6

# Process / entity-category words (プロセス語・実体カテゴリ語) name a biological
# process or a broad category of entities rather than a specific trait or
# function, so on their own they carry no discriminative power that maps a
# phenotype to a function. Substring/stem normalization must therefore not fire
# on them alone; exact phrase/token matches in ``_match`` are unaffected. The
# suffix rule extends the class to compound trailing forms (photosynthesis,
# glycoprotein, subfamily) whose final morpheme is itself a process word.
_PROCESS_TERMS = {
    "synthesis", "syntheses", "synthetic",
    "biosynthesis", "biosynthetic",
    "metabolism", "metabolic", "metabolite", "metabolites",
    "activity", "activities",
    "protein", "proteins",
    "family", "families",
}
_PROCESS_SUFFIXES = ("synthesis", "metabolism", "activity", "protein", "family")


def _is_process_term(word: str) -> bool:
    """True for process/entity-category words that must not fire normalization alone."""
    return word in _PROCESS_TERMS or word.endswith(_PROCESS_SUFFIXES)


def _common_prefix_len(left: str, right: str) -> int:
    limit = min(len(left), len(right))
    offset = 0
    while offset < limit and left[offset] == right[offset]:
        offset += 1
    return offset


def _match_normalized(text: str, words: list[str]) -> list[tuple[str, str]]:
    """Return (word, method) pairs that match ``text`` via normalization rules.

    Substring containment is checked in both directions (the term inside a
    description word, or a description word inside the term) with both sides
    required to clear ``_MIN_SUBSTRING_LEN``. Stem matches require a shared
    prefix of at least ``_MIN_STEM_LEN`` characters (flavonoid/flavonone share
    "flavon" = 6 and match; colour/color share only "colo" = 4 and must not).
    Process/entity-category terms (``_is_process_term``) never fire these rules;
    they carry no discriminative power on their own, so normalization must not
    invent a trait-function link from them alone (issue #400 step 2).
    """
    lowered = text.casefold()
    # Hyphens separate words (consistent with the exact boundary regex in
    # ``_match``, where "-" is not a [0-9a-z] boundary), so split on \w+ to keep
    # e.g. "chalcone-flavonone" as two tokens.
    tokens = [tok for tok in re.findall(r"\w+", lowered) if len(tok) >= 3]
    matched: list[tuple[str, str]] = []
    for word in words:
        if not word or _is_process_term(word):
            continue
        hit: str | None = None
        if len(word) >= _MIN_SUBSTRING_LEN:
            contained = [tok for tok in tokens if len(tok) >= _MIN_SUBSTRING_LEN and tok != word]
            if any(word in tok for tok in contained) or any(tok in word for tok in contained):
                hit = "substring"
        if hit is None and len(word) >= _MIN_STEM_LEN:
            if any(len(tok) >= _MIN_STEM_LEN and _common_prefix_len(word, tok) >= _MIN_STEM_LEN
                   for tok in tokens if tok != word):
                hit = "stem"
        if hit is not None:
            matched.append((word, hit))
    return matched


def _trait_evidence_rows(terms: dict[str, object], gene_id: str, description: str, source_type: str, source_id: str) -> list[dict[str, str]]:
    """Emit a scored trait_function_match for phrase/specific matches, or an
    unscored trait_context row when only generic anatomy/direction words match.

    Exact phrase/token matches keep their existing high/medium confidence. When
    no scored term matches exactly, the normalization layer (``_match_normalized``)
    may still match scored terms by substring or shared stem at low confidence.
    Broad generic-anatomy terms are never normalized, so context behavior is
    unchanged (issue #400 step 1). Process/entity-category terms are likewise
    never normalized, so a process word alone cannot fabricate a
    trait_function_match (issue #400 step 2).
    """
    phrase = _match(description, list(terms["phrases"]))
    scored = _match(description, list(terms["scored"]))
    broad = _match(description, list(terms["broad"]))
    strong = sorted(set(phrase + scored))
    if strong:
        confidence = "high" if phrase else "medium"
        kind = "exact phrase" if phrase else "trait token"
        return [_row(gene_id, "trait_function_match", "support", confidence, ";".join(strong), source_type, source_id,
                     f"{kind} match ({', '.join(strong)}): {description}")]
    normalized = _match_normalized(description, list(terms["scored"]))
    if normalized:
        words = ";".join(sorted({word for word, _ in normalized}))
        methods = ";".join(sorted({method for _, method in normalized}))
        return [_row(gene_id, "trait_function_match", "support", "low", words, source_type, source_id,
                     f"normalized trait match ({methods}): {description}")]
    if broad:
        return [_row(gene_id, "trait_context", "context", "low", ";".join(sorted(set(broad))), source_type, source_id,
                     f"Only generic anatomy/direction term(s) matched; context, not scored: {description}")]
    return []


def _load_approved_terms(path_arg: str) -> list[tuple[str, str]]:
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    if "term" not in frame.columns:
        raise ValueError("approved terms file requires a term column (from `trait-terms build`)")
    weight_col = "weight_class" if "weight_class" in frame.columns else None
    return [(str(row["term"]).strip().casefold(), str(row[weight_col]).strip() if weight_col else "direct")
            for _, row in frame.iterrows() if str(row["term"]).strip()]


def _sorted_top(hits: dict[str, list[dict[str, object]]], keep: int = 5) -> dict[str, list[dict[str, object]]]:
    for query in hits:
        hits[query].sort(key=lambda hit: (-float(hit["bitscore"]), float(hit["evalue"]), str(hit["subject"])))
        hits[query] = hits[query][:keep]
    return hits


def _best_local_hits(proteins: Path, database: Path, out_dir: Path, threads: int) -> dict[str, list[dict[str, object]]]:
    output = out_dir / "local_blastp.tsv"
    command = ["blastp", "-query", str(proteins), "-db", str(database), "-evalue", "1e-5", "-max_target_seqs", "5",
               "-num_threads", str(threads), "-outfmt", "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore staxids stitle", "-out", str(output)]
    completed = run_tool(command, text=True, capture_output=True)
    if completed.returncode: raise RuntimeError(f"local blastp failed: {completed.stderr.strip()}")
    hits: dict[str, list[dict[str, object]]] = {}
    with output.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t", 13)
            if len(fields) < 14: continue
            record = {"subject": fields[1], "identity": float(fields[2]),
                      "query_coverage": int(fields[3]) / max(1, int(fields[4])),
                      "subject_coverage": int(fields[3]) / max(1, int(fields[5])),
                      "evalue": float(fields[10]), "bitscore": float(fields[11]), "description": fields[13]}
            hits.setdefault(fields[0], []).append(record)
    return _sorted_top(hits)


def _auto_local_protein_db() -> Path | None:
    candidates = [os.environ.get("CANDIDATE_QTL_PROTEIN_DB"),
                  str(Path.home() / ".codex" / "blast_databases" / "arabidopsis_uniprot" / "arabidopsis_uniprot_pep")]
    for value in candidates:
        if not value: continue
        path = runtime_path(value)
        if any(Path(str(path) + suffix).exists() for suffix in (".pin", ".pdb")): return path
    return None


DEFAULT_NCBI_QBLAST_TIMEOUT = 60
NCBI_QBLAST_TIMEOUT_ENV = "CANDIDATE_QTL_NCBI_TIMEOUT"
_NCBI_ATTEMPTS = 2
_NCBI_RETRY_BACKOFF = 2.0


def _qblast_error_kind(error: Exception) -> str:
    if isinstance(error, ExternalDependencyError):
        return error.kind
    if isinstance(error, TimeoutError):
        return "timeout"
    return "network"


def _fetch_ncbi_xml(payload: str, email: str | None, timeout: float) -> str:
    """Run NCBIWWW.qblast in a daemon worker bounded by ``timeout`` seconds.

    The worker is abandoned on timeout; daemon threads do not block interpreter
    exit the way ThreadPoolExecutor workers do, so a hung network call cannot
    stall the process at shutdown.
    """
    from Bio.Blast import NCBIWWW

    box: dict[str, object] = {}

    def fetch() -> None:
        try:
            if email:
                NCBIWWW.email = email
            handle = NCBIWWW.qblast("blastp", "swissprot", payload, hitlist_size=5, expect=1e-5, format_type="XML")
            try:
                box["result"] = handle.read()
            finally:
                handle.close()
        except Exception as error:
            box["error"] = error

    thread = threading.Thread(target=fetch, daemon=True)
    thread.start()
    thread.join(timeout)
    if thread.is_alive():
        raise ExternalDependencyError("timeout", f"NCBI qblast for swissprot did not respond within {timeout:g}s")
    if "error" in box:
        raise box["error"]
    return box["result"]


def _fetch_ncbi_xml_with_retry(payload: str, email: str | None) -> str:
    timeout = _env_timeout(NCBI_QBLAST_TIMEOUT_ENV, DEFAULT_NCBI_QBLAST_TIMEOUT)
    last_kind = "network"
    last_error: Exception | None = None
    for attempt in range(_NCBI_ATTEMPTS):
        if attempt:
            time.sleep(_NCBI_RETRY_BACKOFF)
        try:
            return _fetch_ncbi_xml(payload, email, timeout)
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError, ExternalDependencyError) as error:
            last_kind = _qblast_error_kind(error)
            last_error = error
    detail = str(last_error) if last_error is not None else last_kind
    raise ExternalDependencyError(last_kind, f"NCBI qblast for swissprot failed after {_NCBI_ATTEMPTS} attempts: {detail}")


def _write_ncbi_failure_manifest(out: Path, inputs: dict[str, object], extract_settings: dict[str, object], error: ExternalDependencyError) -> None:
    manifest = {
        "schema_version": "1.0",
        "inputs": inputs,
        "annotation_provider": "ncbi_blast_unavailable",
        "external_api_failure": {"kind": error.kind, "detail": str(error)},
        "external_results_used": False,
        "external_sequence_submission_this_run": None,
        "external_cache_used": None,
        "homology_searched": False,
        "workflow_stage": "annotation_failed",
        "review_ready": False,
        "extraction": extract_settings,
        "unrun_stages": ["human_evidence_review", "case_run_source_gate", "manifest_resume", "assay_planning", "wet_validation"],
    }
    (out / "discover_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _best_ncbi_hits(proteins: Path, out_dir: Path, max_candidates: int, email: str | None) -> tuple[dict[str, list[dict[str, object]]], str, bool, bool]:
    from Bio.Blast import NCBIXML
    records = list(SeqIO.parse(proteins, "fasta"))
    if not records:
        return {}, "external_annotation_skipped:no translated candidate proteins", False, False
    if len(records) > max_candidates:
        return {}, f"external_annotation_skipped:{len(records)} candidates exceeds --external-max-candidates {max_candidates}", False, False
    payload = "".join(f">{record.id}\n{record.seq}\n" for record in records)
    query_hash = hashlib.sha256(payload.encode()).hexdigest()
    cache = out_dir / "external_cache" / f"ncbi_swissprot_{query_hash}.xml"; cache.parent.mkdir(parents=True, exist_ok=True)
    cache_used = cache.exists(); submitted = False
    if not cache_used:
        xml_data = _fetch_ncbi_xml_with_retry(payload, email)
        cache.write_text(xml_data, encoding="utf-8")
        submitted = True
    hits: dict[str, list[dict[str, object]]] = {}
    with cache.open(encoding="utf-8") as handle:
        for record in NCBIXML.parse(handle):
            query_id = str(record.query).split()[0]
            for alignment in record.alignments[:5]:
                hsp = alignment.hsps[0]
                hits.setdefault(query_id, []).append({
                    "subject": alignment.hit_id, "identity": 100 * hsp.identities / max(1, hsp.align_length),
                    "query_coverage": hsp.align_length / max(1, record.query_length),
                    "subject_coverage": hsp.align_length / max(1, alignment.length),
                    "evalue": hsp.expect, "bitscore": hsp.bits, "description": alignment.hit_def})
    return _sorted_top(hits), "", submitted, cache_used


def discover_candidates(fasta_arg: str, gff3_arg: str, intervals: list[str], phenotype: str, out_dir_arg: str,
                        direction: str = "unknown", inheritance: str = "unknown", phenotype_terms: list[str] | None = None,
                        protein_db_arg: str | None = None, allow_external_api: bool = False,
                        external_provider: str = "ncbi-blast", external_max_candidates: int = 100,
                        email: str | None = None, top_n: int = DEFAULT_TOP_N, threads: int = 4, flank: int = 500,
                        auto_local_db: bool = True, approved_terms_arg: str | None = None,
                        seqid_map_arg: str | None = None, gff_qc_report_arg: str | None = None,
                        provenance_arg: str | None = None, full_hash: bool = False,
                        evidence_args: list[str] | None = None, dry_run: bool = False) -> dict[str, object]:
    out = runtime_path(out_dir_arg); out.mkdir(parents=True, exist_ok=True)
    extract_dir = out / "01_extract"
    extraction = extract_multi(fasta_arg, gff3_arg, intervals, str(extract_dir), seqid_map_arg, full_hash, gff_qc_report_arg, flank, provenance_arg)
    # Carry the extract manifest's QC/provenance forward so discover records the
    # same seqid mapping, GFF QC, fingerprint strength, and provenance as extract.
    extract_manifest: dict[str, object] = {}
    extract_manifest_path = extract_dir / "extract_manifest.json"
    for candidate_manifest in (extract_dir / "extract_manifest.json", extract_dir / "block_01" / "extract_manifest.json"):
        if candidate_manifest.exists():
            extract_manifest_path = candidate_manifest
            extract_manifest = json.loads(candidate_manifest.read_text(encoding="utf-8"))
            break
    extract_settings = {
        "seqid_map": seqid_map_arg, "gff_qc_report": gff_qc_report_arg,
        "qtl_source_manifest": provenance_arg, "full_hash": full_hash, "flank": flank,
        "extract_manifest": str(extract_manifest_path),
        "seqid_mapping_applied": extract_manifest.get("seqid_mapping_applied"),
        "weak_input_fingerprints": extract_manifest.get("weak_input_fingerprints"),
        "gff_dialect_qc": extract_manifest.get("gff_dialect_qc"),
        "external_gff_qc": extract_manifest.get("external_gff_qc"),
    }
    candidates_path = runtime_path(str(extraction.get("candidates_tsv") or extract_dir / "candidates.tsv"))
    proteins = runtime_path(str(extraction.get("candidate_proteins_fasta") or extract_dir / "candidate_proteins.fa"))
    flagged_path = out / "02_candidate_flags.tsv"; candidate_flags(str(candidates_path), str(flagged_path))
    candidates = pd.read_csv(flagged_path, sep="\t", dtype=str).fillna("")
    approved = _load_approved_terms(approved_terms_arg) if approved_terms_arg else None
    terms = _build_terms(phenotype, phenotype_terms, approved)
    interval_source = ";".join(intervals)

    # Homology search runs before the Dark call is finalized, so classification is
    # decided by results, not processing order.
    annotation_dir = out / "03_annotation"; annotation_dir.mkdir(parents=True, exist_ok=True)
    hits: dict[str, list[dict[str, object]]] = {}; external_warning = ""; external_submitted = False; external_cache_used = False
    provider = "gff3_only"; homology_error = ""
    external_transmission_warning = False
    dry_run_plan: dict[str, object] | None = None
    selected_local_db = runtime_path(protein_db_arg) if protein_db_arg else (_auto_local_protein_db() if auto_local_db else None)
    if dry_run:
        if selected_local_db:
            dry_run_plan = {"planned_provider": "local_blastp", "external_transmission_would_occur": False,
                            "reason": "local protein database selected; no external transmission"}
        elif allow_external_api:
            if external_provider != "ncbi-blast": raise ValueError("supported external provider: ncbi-blast")
            protein_count = sum(1 for _ in SeqIO.parse(proteins, "fasta"))
            would_occur = 0 < protein_count <= external_max_candidates
            skip_reason = "no candidate proteins" if protein_count == 0 else f"{protein_count} candidates exceed --external-max-candidates {external_max_candidates}"
            dry_run_plan = {"planned_provider": "ncbi_blast_swissprot", "external_transmission_would_occur": would_occur,
                            "reason": (f"{protein_count} candidate protein sequences would be sent to the NCBI BLAST API (swissprot)" if would_occur
                                       else f"external search would be skipped: {skip_reason}")}
        else:
            dry_run_plan = {"planned_provider": "gff3_only", "external_transmission_would_occur": False,
                            "reason": "external API not allowed; no external transmission"}
        transmission_state = "WOULD transmit" if dry_run_plan["external_transmission_would_occur"] else "would NOT transmit"
        print(f"discover: DRY RUN: no sequences submitted and no BLAST search ran; a real run {transmission_state} candidate protein sequences externally ({dry_run_plan['reason']}). Rerun without --dry-run to execute.", file=sys.stderr)
    elif selected_local_db:
        try:
            hits = _best_local_hits(proteins, selected_local_db, annotation_dir, threads); provider = "local_blastp"
        except RuntimeError as error:
            provider = "homology_search_failed"; homology_error = str(error)
    elif allow_external_api:
        if external_provider != "ncbi-blast": raise ValueError("supported external provider: ncbi-blast")
        external_transmission_warning = True
        print("discover: WARNING: no local protein database selected, so candidate protein sequences may be transmitted to the external NCBI BLAST API (swissprot). Use --protein-db or --offline (or omit --allow-external-api) to keep sequence processing local.", file=sys.stderr)
        try:
            hits, external_warning, external_submitted, external_cache_used = _best_ncbi_hits(proteins, annotation_dir, external_max_candidates, email)
        except ExternalDependencyError as error:
            _write_ncbi_failure_manifest(out, {"fasta": str(runtime_path(fasta_arg)), "gff3": str(runtime_path(gff3_arg)), "intervals": intervals, "phenotype": phenotype, "direction": direction, "inheritance": inheritance}, extract_settings, error)
            raise
        # A completed search with no hits still counts as searched; only an actual
        # skip (too many candidates, no proteins) leaves Dark unresolved.
        provider = "external_skipped" if external_warning.startswith("external_annotation_skipped") else "ncbi_blast_swissprot"

    qualifying = {gene_id for gene_id, gene_hits in hits.items()
                  if gene_hits and homology_qc.passes_quality(gene_hits[0]["identity"], gene_hits[0]["query_coverage"], gene_hits[0]["evalue"])}
    homology_searched = provider in {"local_blastp", "ncbi_blast_swissprot"}

    def _resolve_dark(scope: str, gene_id: str) -> tuple[str, str]:
        if scope not in {"dark", "dark_pending"}:
            return scope, ""
        if gene_id in qualifying:
            return "positional", "homology_found_reclassified"
        if homology_searched:
            return "dark", "confirmed_no_homology"
        if provider == "homology_search_failed":
            return "dark_pending", "homology_search_failed"
        if provider == "external_skipped":
            return "dark_pending", "homology_search_skipped"
        return "dark_pending", "homology_not_searched"
    resolved = [_resolve_dark(str(getattr(row, "candidate_scope", "")), str(row.gene_id)) for row in candidates.itertuples()]
    candidates["candidate_scope"] = [item[0] for item in resolved]
    candidates["dark_state"] = [item[1] for item in resolved]
    candidates.to_csv(flagged_path, sep="\t", index=False, lineterminator="\n")

    evidence_rows: list[dict[str, str]] = []
    for row in candidates.itertuples():
        # Interval membership defines the candidate set; it is context, not a
        # scored, candidate-distinguishing line of evidence.
        evidence_rows.append(_row(row.gene_id, "interval_membership", "context", "low", row.interval if hasattr(row,"interval") else interval_source,
                                  "coordinate", interval_source, "Candidate overlaps an explicitly supplied target interval (context, not scored)."))
        description = str(getattr(row, "description", ""))
        evidence_rows.extend(_trait_evidence_rows(terms, str(row.gene_id), description, "gff3_annotation", str(runtime_path(gff3_arg))))
        scope = str(getattr(row, "candidate_scope", "")); dark_state = str(getattr(row, "dark_state", ""))
        if scope in {"dark", "dark_pending"}:
            evidence_rows.append(_row(row.gene_id, "annotation_incomplete", "rescue", "medium", dark_state or scope,
                                      "reference_annotation", str(runtime_path(gff3_arg)), f"Short complete ORF retained in the Dark bucket ({dark_state or scope})."))

    audit_rows: list[dict[str, object]] = []
    concordance_rows: list[dict[str, object]] = []

    # Process all candidates for concordance tracking. Concordance is only
    # recorded when a homology search actually ran (hits present); a search that
    # never ran leaves no rows, matching homology_audit.tsv.
    candidate_ids = list(candidates["gene_id"].astype(str)) if hits else []
    for gene_id in candidate_ids:
        gene_hits = hits.get(gene_id, [])
        for hit_rank, hit in enumerate(gene_hits, 1):
            pname = homology_qc.extract_protein_name(hit.get("description", ""), hit.get("subject", ""))
            species = homology_qc.extract_species(hit.get("description", "") or hit.get("subject", ""))
            acc = homology_qc.extract_accession(hit.get("subject", ""))
            reviewed = homology_qc.is_reviewed_entry(hit.get("subject", ""))
            is_info = homology_qc.is_informative_annotation(pname, hit.get("subject", ""))
            audit_rows.append({
                "gene_id": gene_id, "hit_rank": hit_rank, "subject": hit["subject"],
                "identity": round(float(hit["identity"]), 2), "query_coverage": round(float(hit["query_coverage"]), 3),
                "subject_coverage": round(float(hit.get("subject_coverage", float("nan"))), 3) if hit.get("subject_coverage") is not None else "",
                "evalue": hit["evalue"], "bitscore": hit["bitscore"],
                "quality_class": homology_qc.quality_class(hit["identity"], hit["query_coverage"], hit.get("subject_coverage")),
                "provider": provider, "description": hit["description"],
                "species": species, "accession": acc, "reviewed": reviewed,
                "protein_name": pname, "is_informative": is_info,
            })

        gene_qualifying = [h for h in gene_hits[:5] if homology_qc.passes_quality(h["identity"], h["query_coverage"], h["evalue"])]
        state, state_details = homology_qc.assess_annotation_concordance(gene_hits)
        info_hits = [h for h in gene_qualifying if homology_qc.is_informative_annotation(
            homology_qc.extract_protein_name(h.get("description", ""), h.get("subject", "")), h.get("subject", ""))]
        species_set = [homology_qc.extract_species(h.get("description", "") or h.get("subject", "")) for h in gene_qualifying]
        acc_set = [homology_qc.extract_accession(h.get("subject", "")) for h in gene_qualifying]
        concordance_rows.append({
            "gene_id": gene_id,
            "qualifying_hits_count": len(gene_qualifying),
            "informative_hits_count": len(info_hits),
            "concordance_state": state,
            "concordance_details": state_details,
            "species_list": ";".join([s for s in species_set if s]),
            "accessions": ";".join([a for a in acc_set if a]),
        })

        # Keep the original best-hit gate so ranking, dark reclassification and
        # trait evidence stay byte-identical to the pre-#422 behavior; the
        # multi-hit lift below only adds weight-0 context rows.
        if gene_id not in qualifying:
            continue

        # Lift all qualifying hits (up to top 5) for each gene into 04_evidence.tsv as homology_context rows
        for q_rank, q_hit in enumerate(gene_qualifying, 1):
            quality = homology_qc.quality_class(q_hit["identity"], q_hit["query_coverage"], q_hit.get("subject_coverage"))
            conf = homology_qc.confidence(q_hit["identity"], q_hit["query_coverage"], q_hit.get("subject_coverage"), q_hit["evalue"])
            scov = q_hit.get("subject_coverage")
            scov_text = "unknown" if scov is None else f"{float(scov):.3f}"
            if q_rank == 1:
                desc = f"generic best-hit conservation ({quality}; qcov={float(q_hit['query_coverage']):.3f}; scov={scov_text}): {q_hit['description']}"
            else:
                desc = f"generic homology hit #{q_rank} conservation ({quality}; qcov={float(q_hit['query_coverage']):.3f}; scov={scov_text}): {q_hit['description']}"
            # Generic homology is conservation context, not scored trait support (weight 0)
            evidence_rows.append(_row(gene_id, "homology_context", "context", conf, q_hit["subject"], provider, q_hit["subject"], desc))

        best = gene_hits[0]
        if len(gene_hits) > 1 and homology_qc.is_paralog_ambiguous(best["bitscore"], gene_hits[1]["bitscore"]):
            evidence_rows.append(_row(gene_id, "paralog_ambiguous", "against", "medium",
                                      f"top2_bitscore_ratio={homology_qc.bitscore_ratio(best['bitscore'], gene_hits[1]['bitscore']):.3f}",
                                      provider, best["subject"], "Top two similarity hits are too close for a simple one-hit ortholog assumption."))
        evidence_rows.extend(_trait_evidence_rows(terms, gene_id, str(best["description"]), provider, best["subject"]))

    if audit_rows:
        pd.DataFrame(audit_rows).to_csv(annotation_dir / "homology_audit.tsv", sep="\t", index=False, lineterminator="\n")
    if concordance_rows:
        pd.DataFrame(concordance_rows).to_csv(annotation_dir / "annotation_concordance.tsv", sep="\t", index=False, lineterminator="\n")

    evidence = pd.DataFrame(evidence_rows, columns=EVIDENCE_COLUMNS)
    external_frames = [pd.read_csv(runtime_path(path), sep="\t", dtype=str).fillna("") for path in (evidence_args or [])]
    excluded = pd.DataFrame()
    if external_frames:
        external_evidence = pd.concat(external_frames, ignore_index=True, sort=False).fillna("")
        for column in EVIDENCE_COLUMNS:
            if column not in external_evidence.columns:
                external_evidence[column] = ""
        external_evidence = external_evidence[external_evidence["gene_id"].isin(set(candidates["gene_id"]))].copy()
        excluded = external_evidence[external_evidence["source_id"].str.strip().eq("")].copy()
        external_evidence = external_evidence[external_evidence["source_id"].str.strip().ne("")].copy()
        evidence = pd.concat([evidence, external_evidence], ignore_index=True, sort=False).fillna("")

    direction_matches = direction_opposites = inheritance_matches = inheritance_mismatches = 0
    if direction != "unknown" and "phenotype_direction" in evidence.columns:
        opposite = {"increase": "decrease", "decrease": "increase", "gain": "loss", "loss": "gain"}
        values = evidence["phenotype_direction"].str.strip().str.casefold()
        evidence["direction_match"] = "unknown"
        match_mask = values.eq(direction); opposite_mask = values.eq(opposite.get(direction, "__none__"))
        evidence.loc[match_mask, "direction_match"] = "match"
        evidence.loc[opposite_mask, "direction_match"] = "opposite"
        evidence.loc[opposite_mask, "direction"] = "against"
        direction_matches, direction_opposites = int(match_mask.sum()), int(opposite_mask.sum())
    if inheritance != "unknown" and {"inheritance_model", "molecular_effect"}.issubset(evidence.columns):
        models = evidence["inheritance_model"].str.strip().str.casefold()
        mechanism_backed = evidence["molecular_effect"].str.strip().ne("") & evidence["source_id"].str.strip().ne("") & models.ne("")
        evidence["inheritance_match"] = "unknown"
        match_mask = mechanism_backed & models.eq(inheritance); mismatch_mask = mechanism_backed & ~models.eq(inheritance)
        evidence.loc[match_mask, "inheritance_match"] = "match"
        evidence.loc[mismatch_mask, "inheritance_match"] = "mismatch"
        evidence.loc[mismatch_mask, "direction"] = "against"
        inheritance_matches, inheritance_mismatches = int(match_mask.sum()), int(mismatch_mask.sum())
    evidence_path = out / "04_evidence.tsv"
    evidence.sort_values(["gene_id","evidence_code","source_id"], kind="stable").to_csv(evidence_path, sep="\t", index=False, lineterminator="\n")
    excluded_path = out / "evidence_excluded.tsv"
    excluded.to_csv(excluded_path, sep="\t", index=False, lineterminator="\n")
    rank_dir = out / "05_rank"; ranking = prioritize_evidence(str(flagged_path), str(evidence_path), str(rank_dir), None, top_n, top_n, DEFAULT_MIN_TOP_SCORE, top_n)
    ranked = pd.read_csv(rank_dir / "ranked_candidates.tsv", sep="\t", dtype=str).fillna("")
    shortlist = ranked[ranked.dry_group.isin({"top","reserve","rescue"})].copy()
    shortlist_path = out / "candidate_shortlist.tsv"; shortlist.to_csv(shortlist_path, sep="\t", index=False, lineterminator="\n")
    review_rows = []
    if provider == "gff3_only": review_rows.append({"reason":"protein_annotation_not_configured","required_next_action":"register a local protein BLAST DB or explicitly allow the NCBI API","allowed_choices":"local --protein-db | --allow-external-api | review GFF3-only candidates"})
    if provider == "homology_search_failed": review_rows.append({"reason":f"homology_search_failed: {homology_error}","required_next_action":"fix the BLAST database/tool and rerun; Dark calls stay provisional until a search completes","allowed_choices":"repair local database | narrower interval"})
    if external_warning: review_rows.append({"reason":external_warning,"required_next_action":"use a local BLAST DB for this candidate count","allowed_choices":"local database | narrower interval"})
    if not excluded.empty: review_rows.append({"reason":f"evidence_without_source_id:{len(excluded)}","required_next_action":"add a stable publication, database, or local-analysis source_id before ranking","allowed_choices":"supply source_id | keep quarantined"})
    review_path = out / "review_required.tsv"; pd.DataFrame(review_rows, columns=["reason","required_next_action","allowed_choices"]).to_csv(review_path, sep="\t", index=False, lineterminator="\n")
    report = out / "report.md"; counts = ranked.dry_group.value_counts().to_dict()
    direction_used = direction_matches + direction_opposites > 0
    inheritance_used = inheritance_matches + inheritance_mismatches > 0
    lines = [f"# Candidate discovery: {phenotype}", "", "- Workflow stage: `preliminary_screening`", f"- Intervals: `{interval_source}`", f"- Direction: `{direction}` (used for ranking: `{str(direction_used).lower()}`)", f"- Inheritance: `{inheritance}` (used for ranking: `{str(inheritance_used).lower()}`)", f"- Annotation provider: `{provider}`", "", "## Groups", ""]
    lines += [f"- {key}: {counts[key]}" for key in sorted(counts)]
    lines += ["", "## Candidates", "", "| gene_id | group | score | reason |", "|---|---:|---:|---|"]
    for row in ranked.head(max(20, top_n)).itertuples(): lines.append(f"| {row.gene_id} | {row.dry_group} | {row.dry_score} | {row.score_breakdown} |")

    # Top Candidate Homology Evidence section
    top_candidates = ranked.head(max(20, top_n))
    lines += [
        "",
        "## Top Candidate Homology Evidence",
        "",
        f"> {homology_qc.HOMOLOGY_DISCLAIMER}",
        "",
        "| gene_id | group | concordance | homology_rationale |",
        "|---|---|---|---|",
    ]
    concordance_by_gene = {row["gene_id"]: row for row in concordance_rows} if concordance_rows else {}
    for row in top_candidates.itertuples():
        gene_id = str(row.gene_id)
        gene_concordance = concordance_by_gene.get(gene_id)
        if gene_concordance is None:
            c_state = "not_searched" if not hits else "no_hits"
        else:
            c_state = gene_concordance["concordance_state"]
        gene_hits_list = hits.get(gene_id, [])
        rationale = homology_qc.format_grounded_rationale(hits=gene_hits_list, include_disclaimer=True)
        lines.append(f"| {gene_id} | {row.dry_group} | {c_state} | {rationale} |")

    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    manifest = {"schema_version":"1.0", "inputs":{"fasta":str(runtime_path(fasta_arg)),"gff3":str(runtime_path(gff3_arg)),"intervals":intervals,"phenotype":phenotype,"direction":direction,"inheritance":inheritance},
                "annotation_provider":provider,"external_api_allowed":allow_external_api,"external_transmission_warning":external_transmission_warning,"external_results_used":provider.startswith("ncbi_"),
                "external_sequence_submission_this_run":external_submitted,"external_cache_used":external_cache_used,
                "local_protein_db":str(selected_local_db) if selected_local_db else None,"phenotype_terms":terms,
                "extraction":extract_settings,
                "homology_searched":homology_searched,"dark_state_counts":{str(k):int(v) for k,v in candidates["dark_state"].value_counts().items() if str(k)},
                "homology_audit":str(annotation_dir / "homology_audit.tsv") if (annotation_dir / "homology_audit.tsv").exists() else None,
                "annotation_concordance":str(annotation_dir / "annotation_concordance.tsv") if (annotation_dir / "annotation_concordance.tsv").exists() else None,
                "annotation_concordance_counts":{str(k):int(v) for k,v in pd.DataFrame(concordance_rows)["concordance_state"].value_counts().items()} if concordance_rows else {},
                "homology_qc_policy":{"min_identity":homology_qc.MIN_IDENTITY,"min_query_coverage":homology_qc.MIN_QUERY_COVERAGE,"max_evalue":homology_qc.MAX_EVALUE,"paralog_bitscore_ratio":homology_qc.PARALOG_BITSCORE_RATIO},
                "workflow_stage":"preliminary_screening","review_ready":False,
                "dry_run":dry_run,"dry_run_plan":dry_run_plan,
                "evidence_constraints":{"direction":{"value":direction,"used_for_ranking":direction_used,"matches":direction_matches,"opposites":direction_opposites,"unused_reason":None if direction_used else "no source-backed evidence row supplied phenotype_direction"},
                                        "inheritance":{"value":inheritance,"used_for_ranking":inheritance_used,"matches":inheritance_matches,"mismatches":inheritance_mismatches,"unused_reason":None if inheritance_used else "no source-backed evidence row supplied both inheritance_model and molecular_effect"}},
                "unrun_stages":["human_evidence_review","case_run_source_gate","manifest_resume","assay_planning","wet_validation"],
                "candidate_count":len(candidates),"shortlist_count":len(shortlist),"groups":{str(k):int(v) for k,v in counts.items()},
                "ranking_defaults":{"min_top_score":DEFAULT_MIN_TOP_SCORE,"top_n":top_n,"reserve_n":top_n,"rescue_n":top_n},
                "outputs":{"shortlist":str(shortlist_path),"ranked":ranking["ranked_candidates"],"evidence":str(evidence_path),"evidence_excluded":str(excluded_path),"report":str(report),"review_required":str(review_path)}}
    manifest_path = out / "discover_manifest.json"; manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    return manifest | {"manifest":str(manifest_path)}


def promote_discover(discover_dir_arg: str, out_dir_arg: str, case_id: str | None = None) -> dict[str, object]:
    discover_dir, out = runtime_path(discover_dir_arg), runtime_path(out_dir_arg)
    manifest_path = discover_dir / "discover_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    candidates = discover_dir / "01_extract" / "candidates.tsv"
    evidence = discover_dir / "04_evidence.tsv"
    if not candidates.exists() or not evidence.exists():
        raise ValueError("discover output requires 01_extract/candidates.tsv and 04_evidence.tsv")
    out.mkdir(parents=True, exist_ok=True)
    identifier = case_id or re.sub(r"[^A-Za-z0-9_.-]+", "_", str(manifest["inputs"]["phenotype"])).strip("_") or "promoted_case"
    intervals = manifest["inputs"]["intervals"]
    # Carry the discover run's effective ranking defaults into the case manifest
    # so case-run ranks the same way the Dry shortlist was built; without
    # this the case-run lane falls back to its own stricter defaults.
    ranking_defaults = manifest.get("ranking_defaults") or {}
    min_top_score = float(ranking_defaults.get("min_top_score", DEFAULT_MIN_TOP_SCORE))
    top_n = int(ranking_defaults.get("top_n", DEFAULT_TOP_N))
    reserve_n = int(ranking_defaults.get("reserve_n", top_n))
    rescue_n = int(ranking_defaults.get("rescue_n", top_n))
    case = {"case_id": identifier,
            "trait":{"name":manifest["inputs"]["phenotype"],"phenotype_direction":manifest["inputs"].get("direction", "unknown"),"inheritance":manifest["inputs"].get("inheritance", "unknown")},
            "source":{"source_id":"REVIEW_REQUIRED","target_intervals":intervals,"coordinate_basis":"physical","literature_search_status":"incomplete"},
            "target":{"assembly":"REVIEW_REQUIRED"},
            "inputs":{"candidates_resource":"discover_candidates","evidence_resources":["discover_evidence"]},
            "ranking_policy":"deterministic_registry","min_top_score":min_top_score,
            "assay":{"top_n":top_n,"reserve_n":reserve_n,"rescue_n":rescue_n}}
    resources = {"resources":{"discover_candidates":{"type":"candidates_tsv","path":str(candidates)},"discover_evidence":{"type":"canonical_evidence","path":str(evidence)}}}
    case_path, resources_path = out / "case.yaml", out / "resources.yaml"
    case_path.write_text(yaml.safe_dump(case, sort_keys=False, allow_unicode=True), encoding="utf-8")
    resources_path.write_text(yaml.safe_dump(resources, sort_keys=False, allow_unicode=True), encoding="utf-8")
    review = out / "PROMOTION_REVIEW.md"
    review.write_text("# Discover promotion review\n\nReplace `REVIEW_REQUIRED` source/assembly fields, review quarantined evidence, then run:\n\n```text\ncandidate-qtl case-run --case case.yaml --resources resources.yaml --out-dir case_run\n```\n", encoding="utf-8")
    return {"case":str(case_path),"resources":str(resources_path),"review":str(review),"review_ready":False}
