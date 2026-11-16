from __future__ import annotations

import gzip
import hashlib
import importlib.metadata
import json
import os
import platform
import re
import shutil
import shlex
import statistics
import subprocess
import sys
import unicodedata
import urllib.error
from datetime import UTC, datetime
from collections import defaultdict
from pathlib import Path
from typing import Iterable
from urllib.parse import unquote, urlencode
from urllib.request import Request, urlopen

import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO

from .optional_deps import PYSAM_DEPENDENT_COMMANDS, pysam_available, require_pysam


WINDOWS_PATH = re.compile(r"^([A-Za-z]):[\\/](.*)$")
INTERVAL = re.compile(r"^([^:]+):(\d+)-(\d+)$")


class ExternalDependencyError(Exception):
    """A temporary environment failure: network outage, service timeout, or an external tool exceeding its timeout."""

    def __init__(self, kind: str, message: str) -> None:
        super().__init__(f"external dependency failure (kind={kind}): {message}")
        self.kind = kind


DEFAULT_EXTERNAL_TOOL_TIMEOUT = 300
TOOL_TIMEOUT_ENV = "CANDIDATE_QTL_TOOL_TIMEOUT"


def _env_timeout(env_name: str, default: float) -> float:
    raw = os.environ.get(env_name, "").strip()
    try:
        return float(raw) if raw else default
    except ValueError:
        return default


def run_tool(command, timeout=None, **kwargs):
    """Run an external tool (blastp, minimap2, Rscript, ...) under a default timeout.

    A tool that exceeds its timeout is an environment failure
    (ExternalDependencyError kind=timeout), not a defect in this package.
    """
    if timeout is None:
        timeout = _env_timeout(TOOL_TIMEOUT_ENV, DEFAULT_EXTERNAL_TOOL_TIMEOUT)
    try:
        return subprocess.run(command, timeout=timeout, **kwargs)
    except subprocess.TimeoutExpired as error:
        raise ExternalDependencyError("timeout", f"external tool {command!r} exceeded the {timeout:g}s timeout") from error


def runtime_path(value: str | Path) -> Path:
    text = str(value)
    match = WINDOWS_PATH.match(text)
    if match and os.name != "nt":
        return Path("/mnt") / match.group(1).lower() / match.group(2).replace("\\", "/")
    if os.name != "nt" and "\\" in text:
        return Path(text.replace("\\", "/"))
    return Path(text)


def parse_interval(value: str) -> tuple[str, int, int]:
    value = str(value).strip()
    match = INTERVAL.match(value)
    if not match:
        raise ValueError("interval must be CHROM:START-END with 1-based inclusive coordinates")
    chrom, start, end = match.groups()
    if int(start) > int(end):
        raise ValueError("interval start must be less than or equal to end")
    return chrom, int(start), int(end)


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open("r", encoding="utf-8")


def attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[unquote(key)] = unquote(value)
    return result


def parse_seqid_map(value: str | None) -> dict[str, str]:
    """Parse GFF_SEQID=FASTA_SEQID pairs or a two-column TSV file."""
    if not value:
        return {}
    path = runtime_path(value)
    if path.exists():
        mapping: dict[str, str] = {}
        with open_text(path) as handle:
            for line_number, raw in enumerate(handle, 1):
                if not raw.strip() or raw.startswith("#"):
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) != 2 or not all(fields):
                    raise ValueError(f"{path}:{line_number}: seqid map must have two non-empty TSV columns")
                mapping[fields[0]] = fields[1]
        return mapping
    mapping = {}
    for item in value.split(","):
        if "=" not in item:
            raise ValueError("seqid map must be GFF_SEQID=FASTA_SEQID or a two-column TSV file")
        source, target = item.split("=", 1)
        if not source or not target:
            raise ValueError("seqid map entries must have non-empty GFF and FASTA sequence IDs")
        mapping[source] = target
    return mapping


def iter_gff_rows(path: Path, chrom: str):
    with open_text(path) as handle:
        for raw in handle:
            if not raw or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != 9 or fields[0] != chrom:
                continue
            yield {
                "seqid": fields[0],
                "type": fields[2],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
                "phase": fields[7],
                "attrs": attributes(fields[8]),
            }


def safe_id(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


FULL_HASH_MAX_BYTES = 512 * 1024 * 1024

# How strongly a fingerprint proves content identity, weakest to strongest.
FINGERPRINT_STRENGTH = {
    "full_sha256": "strong",
    "metadata_plus_index_hash": "index",
    "metadata_only": "weak",
}


def file_fingerprint(path: Path, full_hash: bool = False, full_hash_max_bytes: int = FULL_HASH_MAX_BYTES) -> dict[str, object]:
    """Fingerprint large references without rereading multi-gigabyte FASTA by default.

    The ``method`` names a fixed enum (full_sha256 / metadata_plus_index_hash /
    metadata_only) and ``strength`` grades how strongly content identity is
    proven (strong / index / weak) so a manifest reader can see, per input,
    whether byte-identity was actually verified.
    """
    stat = path.stat()
    result: dict[str, object] = {"path": str(path), "size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    if full_hash or stat.st_size <= full_hash_max_bytes:
        result["sha256"] = sha256(path)
        result["method"] = "full_sha256"
        result["strength"] = FINGERPRINT_STRENGTH["full_sha256"]
        return result
    index = Path(str(path) + ".fai")
    result["method"] = "metadata_plus_index_hash" if index.exists() else "metadata_only"
    result["strength"] = FINGERPRINT_STRENGTH[result["method"]]
    result["sha256"] = None
    if index.exists():
        result["fai_path"] = str(index)
        result["fai_sha256"] = sha256(index)
    return result


def fingerprint_is_weak(fingerprint: dict[str, object]) -> bool:
    """True when a fingerprint did not verify full byte content."""
    return str(fingerprint.get("strength", "weak")) != "strong"


def write_fasta(path: Path, records: Iterable[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        seen_ids: set[str] = set()
        for record_id, sequence in records:
            record_id = str(record_id)
            if not record_id or ">" in record_id or any(char.isspace() or unicodedata.category(char) == "Cc" for char in record_id):
                raise ValueError("FASTA record ID must be non-empty and contain no whitespace, control characters, or '>'")
            if record_id in seen_ids:
                raise ValueError(f"duplicate FASTA record ID: {record_id}")
            seen_ids.add(record_id)
            handle.write(f">{record_id}\n")
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset : offset + 80] + "\n")


def extract_candidates(
    fasta_arg: str,
    gff_arg: str,
    interval_arg: str,
    out_arg: str,
    seqid_map_arg: str | None = None,
    full_hash: bool = False,
) -> dict[str, object]:
    fasta_path, gff_path, out_dir = runtime_path(fasta_arg), runtime_path(gff_arg), runtime_path(out_arg)
    chrom, interval_start, interval_end = parse_interval(interval_arg)
    fasta_chrom = parse_seqid_map(seqid_map_arg).get(chrom, chrom)
    out_dir.mkdir(parents=True, exist_ok=True)
    feature_rows = list(iter_gff_rows(gff_path, chrom))
    genes = [
        row
        for row in feature_rows
        if row["type"] == "gene" and not (int(row["end"]) < interval_start or int(row["start"]) > interval_end)
    ]
    gene_ids = {
        str(row["attrs"].get("ID") or row["attrs"].get("Name") or f"{chrom}_{row['start']}_{row['end']}")
        for row in genes
    }
    gene_features = {
        str(row["attrs"].get("ID") or row["attrs"].get("Name") or f"{chrom}_{row['start']}_{row['end']}"): row
        for row in genes
    }
    transcripts = [
        row
        for row in feature_rows
        if row["type"] in {"mRNA", "transcript"}
        and any(parent in gene_ids for parent in str(row["attrs"].get("Parent", "")).split(","))
    ]
    transcript_to_gene: dict[str, str] = {}
    transcript_features: dict[str, dict[str, object]] = {}
    for row in transcripts:
        transcript_id = str(row["attrs"].get("ID", ""))
        parents = [parent for parent in str(row["attrs"].get("Parent", "")).split(",") if parent in gene_ids]
        if transcript_id and parents:
            transcript_to_gene[transcript_id] = parents[0]
            transcript_features[transcript_id] = row
    transcript_ids = set(transcript_to_gene)
    transcript_cds: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in feature_rows:
        if row["type"] != "CDS":
            continue
        for parent in str(row["attrs"].get("Parent", "")).split(","):
            if parent in transcript_ids:
                transcript_cds[parent].append(row)
            elif parent in gene_ids:
                transcript_to_gene.setdefault(parent, parent)
                transcript_features.setdefault(parent, gene_features[parent])
                transcript_cds[parent].append(row)
    pysam = require_pysam("extract")
    try:
        fasta = pysam.FastaFile(str(fasta_path))
    except OSError:
        pysam.faidx(str(fasta_path))
        fasta = pysam.FastaFile(str(fasta_path))
    if fasta_chrom not in fasta.references:
        available = ", ".join(fasta.references[:7])
        fasta.close()
        raise ValueError(
            f"GFF sequence ID {chrom!r} maps to FASTA sequence ID {fasta_chrom!r}, which is absent; "
            f"available examples: {available}. Supply --seqid-map GFF_SEQID=FASTA_SEQID."
        )
    table_rows: list[dict[str, object]] = []
    transcript_rows: list[dict[str, object]] = []
    genomic_records: list[tuple[str, str]] = []
    representative_cds_records: list[tuple[str, str]] = []
    transcript_cds_records: list[tuple[str, str]] = []
    representative_protein_records: list[tuple[str, str]] = []
    for gene in sorted(genes, key=lambda item: (int(item["start"]), int(item["end"]))):
        attrs = gene["attrs"]
        gene_id = str(attrs.get("ID") or attrs.get("Name") or f"{chrom}_{gene['start']}_{gene['end']}")
        sequence = fasta.fetch(fasta_chrom, int(gene["start"]) - 1, int(gene["end"]))
        if gene["strand"] == "-":
            sequence = str(Seq(sequence).reverse_complement())
        genomic_records.append((gene_id, sequence))
        isoforms: list[dict[str, object]] = []
        for transcript_id, parent_gene in transcript_to_gene.items():
            if parent_gene != gene_id:
                continue
            feature = transcript_features.get(transcript_id, gene)
            cds_parts = sorted(transcript_cds.get(transcript_id, []), key=lambda item: int(item["start"]))
            cds_sequence = "".join(fasta.fetch(fasta_chrom, int(part["start"]) - 1, int(part["end"])) for part in cds_parts).upper()
            if gene["strand"] == "-" and cds_sequence:
                cds_sequence = str(Seq(cds_sequence).reverse_complement())
            protein = str(Seq(cds_sequence).translate()) if cds_sequence and len(cds_sequence) % 3 == 0 else ""
            complete_orf = bool(protein and cds_sequence.startswith("ATG") and protein.endswith("*") and "*" not in protein[:-1])
            tag_text = " ".join(str(feature["attrs"].get(key, "")) for key in ("tag", "biotype", "transcript_type")).lower()
            canonical = any(term in tag_text for term in ("canonical", "representative", "mane_select"))
            isoforms.append(
                {
                    "gene_id": gene_id,
                    "transcript_id": transcript_id,
                    "transcript_start": feature["start"],
                    "transcript_end": feature["end"],
                    "strand": gene["strand"],
                    "cds_length": len(cds_sequence),
                    "cds_sequence": cds_sequence,
                    "protein_sequence": protein[:-1] if protein.endswith("*") else protein,
                    "complete_orf": complete_orf,
                    "canonical_tag": canonical,
                }
            )
            if cds_sequence:
                transcript_cds_records.append((transcript_id, cds_sequence))
        isoforms.sort(
            key=lambda row: (
                -int(bool(row["canonical_tag"])),
                -int(bool(row["complete_orf"])),
                -int(row["cds_length"]),
                str(row["transcript_id"]),
            )
        )
        representative = isoforms[0] if isoforms else None
        if representative and representative["cds_sequence"]:
            representative_cds_records.append((gene_id, str(representative["cds_sequence"])))
        if representative and representative["protein_sequence"]:
            representative_protein_records.append((gene_id, str(representative["protein_sequence"])))
        for isoform in isoforms:
            is_representative = bool(representative and isoform["transcript_id"] == representative["transcript_id"])
            reason = ""
            if is_representative:
                if isoform["canonical_tag"]:
                    reason = "canonical annotation tag"
                elif isoform["complete_orf"]:
                    reason = "longest complete ORF"
                else:
                    reason = "longest annotated CDS"
            transcript_rows.append(
                {
                    key: value
                    for key, value in isoform.items()
                    if key not in {"cds_sequence", "protein_sequence"}
                }
                | {"is_representative": is_representative, "selection_reason": reason}
            )
        representative_length = int(representative["cds_length"]) if representative else 0
        annotation_quality = "no_cds"
        if representative:
            annotation_quality = "complete_orf" if representative["complete_orf"] else "partial_or_unverified_orf"
        table_rows.append(
            {
                "gene_id": gene_id,
                "gene_name": attrs.get("Name", ""),
                "description": attrs.get("Note", attrs.get("description", "")),
                "dbxref": attrs.get("Dbxref", ""),
                "ontology_terms": attrs.get("Ontology_term", ""),
                "interval": interval_arg,
                "chrom": chrom,
                "start": gene["start"],
                "end": gene["end"],
                "strand": gene["strand"],
                "gene_length": int(gene["end"]) - int(gene["start"]) + 1,
                "transcript_count": len(isoforms),
                "representative_transcript": representative["transcript_id"] if representative else "",
                "representative_reason": next((row["selection_reason"] for row in transcript_rows if row["gene_id"] == gene_id and row["is_representative"]), ""),
                "cds_length": representative_length,
                "annotation_quality": annotation_quality,
                "overlaps_qtl_boundary": int(gene["start"]) < interval_start or int(gene["end"]) > interval_end,
            }
        )
    fasta.close()
    table_path = out_dir / "candidates.tsv"
    transcript_path = out_dir / "transcripts.tsv"
    candidate_columns = [
        "gene_id", "gene_name", "description", "dbxref", "ontology_terms", "interval", "chrom", "start", "end", "strand", "gene_length",
        "transcript_count", "representative_transcript", "representative_reason", "cds_length",
        "annotation_quality", "overlaps_qtl_boundary",
    ]
    transcript_columns = [
        "gene_id", "transcript_id", "transcript_start", "transcript_end", "strand", "cds_length",
        "complete_orf", "canonical_tag", "is_representative", "selection_reason",
    ]
    pd.DataFrame(table_rows, columns=candidate_columns).to_csv(table_path, sep="\t", index=False)
    pd.DataFrame(transcript_rows, columns=transcript_columns).to_csv(transcript_path, sep="\t", index=False)
    write_fasta(out_dir / "candidate_genomic.fa", genomic_records)
    write_fasta(out_dir / "candidate_cds.fa", representative_cds_records)
    write_fasta(out_dir / "candidate_transcript_cds.fa", transcript_cds_records)
    write_fasta(out_dir / "candidate_proteins.fa", representative_protein_records)
    manifest = {
        "interval": interval_arg,
        "gff_seqid": chrom,
        "fasta_seqid": fasta_chrom,
        "seqid_mapping_applied": chrom != fasta_chrom,
        "fasta": str(fasta_path),
        "fasta_fingerprint": file_fingerprint(fasta_path, full_hash),
        "gff3": str(gff_path),
        "gff3_fingerprint": file_fingerprint(gff_path, full_hash),
        "candidate_count": len(table_rows),
        "with_cds_count": len(representative_cds_records),
        "transcript_count": len(transcript_rows),
        "outputs": {
            "candidates_tsv": str(table_path),
            "transcripts_tsv": str(transcript_path),
            "candidate_genomic_fasta": str(out_dir / "candidate_genomic.fa"),
            "candidate_cds_fasta": str(out_dir / "candidate_cds.fa"),
            "candidate_transcript_cds_fasta": str(out_dir / "candidate_transcript_cds.fa"),
            "candidate_proteins_fasta": str(out_dir / "candidate_proteins.fa"),
        },
    }
    (out_dir / "extract_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


VARIANT_POINTS = {
    "stop_gain": 5,
    "frameshift": 5,
    "splice": 5,
    "splice_site": 5,
    "large_indel": 4,
    "sv_suspected": 4,
    "missense": 3,
    "promoter_indel": 2,
    "utr": 1,
    "synonymous": 0,
}


def nonempty(value: object) -> bool:
    return str(value).strip().lower() not in {"", "na", "nan", "none", "required", "pending"}


DRY_EVIDENCE_WEIGHTS = {
    "model_phenotype_evidence": 3,
    "domain_pathway_evidence": 1,
    "expression_relevance": 2,
    "trait_ontology_match": 1,
    "known_qtg_or_qtl_evidence": 2,
    "synteny_evidence": 1,
}

WET_ONLY_COLUMNS = {
    "parental_variant",
    "parental_variant_evidence",
    "variant_class",
    "variant_validation",
    "marker_validation_status",
    "f2_cosegregation_status",
    "f2_or_f2_3_evidence",
    "recombinant_evidence",
    "functional_validation_evidence",
}


def dry_prioritize(candidate_arg: str, out_arg: str, top_n: int = 7, reserve_n: int = 3) -> dict[str, object]:
    """Prioritize Wet assays from pre-Wet evidence only.

    The score intentionally rewards evidence categories, not publication or GO
    term counts. Rescue candidates are carried in a separate group so a weak
    annotation cannot silently remove them.
    """
    candidate_path, out_path = runtime_path(candidate_arg), runtime_path(out_arg)
    if top_n < 0 or reserve_n < 0:
        raise ValueError("top_n and reserve_n must be zero or greater")
    frame = pd.read_csv(candidate_path, sep="\t", dtype=str).fillna("")
    if "gene_id" not in frame.columns:
        raise ValueError("candidate table must contain gene_id")
    if frame["gene_id"].str.strip().eq("").any():
        raise ValueError("candidate table contains an empty gene_id")
    duplicated = sorted(frame.loc[frame["gene_id"].duplicated(keep=False), "gene_id"].unique())
    if duplicated:
        raise ValueError(f"candidate table contains duplicate gene_id values: {', '.join(duplicated)}")

    leaked = [column for column in sorted(WET_ONLY_COLUMNS) if column in frame.columns and frame[column].map(nonempty).any()]
    if leaked:
        raise ValueError(
            "Dry priority input contains Wet-only evidence; remove these values before ranking: "
            + ", ".join(leaked)
        )

    def score_row(row: pd.Series) -> tuple[int, str]:
        points: list[tuple[str, int]] = []
        orthology = str(row.get("ortholog_confidence", "")).strip().lower()
        if nonempty(orthology):
            if any(term in orthology for term in ("high", "probable", "rbh")):
                points.append(("orthology", 2))
            elif "ambiguous" in orthology or "homology_only" in orthology:
                points.append(("orthology_ambiguous", -1))
            else:
                points.append(("orthology", 1))
        for column, weight in DRY_EVIDENCE_WEIGHTS.items():
            if nonempty(row.get(column, "")):
                points.append((column, weight))
        paralog_risk = str(row.get("paralog_redundancy_risk", "")).strip().lower()
        if nonempty(paralog_risk):
            if "high" in paralog_risk:
                points.append(("high_paralog_risk", -2))
            elif "medium" in paralog_risk:
                points.append(("medium_paralog_risk", -1))
        annotation = str(row.get("annotation_quality", "")).strip().lower()
        if any(term in annotation for term in ("complete", "high")):
            points.append(("annotation_quality", 1))
        if nonempty(row.get("key_counterevidence", "")):
            points.append(("counterevidence", -2))
        total = sum(weight for _, weight in points)
        breakdown = ";".join(f"{label}:{weight:+d}" for label, weight in points) or "no_scored_evidence:0"
        return total, breakdown

    scored = frame.apply(score_row, axis=1)
    frame["knowledge_prior_score"] = [item[0] for item in scored]
    frame["score_breakdown"] = [item[1] for item in scored]
    frame["_rescue"] = frame.get("rescue_reason", pd.Series("", index=frame.index)).map(nonempty)
    frame = frame.sort_values(["knowledge_prior_score", "gene_id"], ascending=[False, True], kind="stable")
    frame["dry_group"] = "not_selected_yet"
    eligible = frame.index[~frame["_rescue"]]
    frame.loc[eligible[:top_n], "dry_group"] = "top"
    frame.loc[eligible[top_n:top_n + reserve_n], "dry_group"] = "reserve"
    frame.loc[frame["_rescue"], "dry_group"] = "rescue"
    group_order = {"top": 0, "reserve": 1, "rescue": 2, "not_selected_yet": 3}
    frame["_group_order"] = frame["dry_group"].map(group_order)
    frame = frame.sort_values(["_group_order", "knowledge_prior_score", "gene_id"], ascending=[True, False, True], kind="stable")
    frame["dry_rank"] = frame.groupby("dry_group", sort=False).cumcount() + 1
    frame = frame.drop(columns=["_rescue", "_group_order"])

    # Every candidate keeps an explicit selection or non-selection reason, and
    # Top/Reserve genes carry recommended parental-sequencing regions. These are
    # ordering decisions, not causal claims; not_selected_yet is deliberately not
    # an exclusion.
    default_regions = "all CDS; exon-intron junctions; UTRs; proximal promoter if mechanistically relevant"

    def _selection_reason(row: pd.Series) -> str:
        group = row["dry_group"]
        if group == "top":
            return "top: highest pre-Wet knowledge prior; sequence parents first"
        if group == "reserve":
            return "reserve: sequence if Top candidates are excluded"
        if group == "rescue":
            reason = str(row.get("rescue_reason", "")).strip()
            return f"rescue: retained independent of knowledge score ({reason})" if reason else "rescue: retained independent of knowledge score"
        return "not_selected_yet: not scheduled for Wet; not excluded as a causal candidate"

    frame["selection_reason"] = frame.apply(_selection_reason, axis=1)
    if "recommended_regions" not in frame.columns:
        frame["recommended_regions"] = ""
    sequencing_mask = frame["dry_group"].isin({"top", "reserve"}) & ~frame["recommended_regions"].map(nonempty)
    frame.loc[sequencing_mask, "recommended_regions"] = default_regions

    out_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out_path, sep="\t", index=False)
    counts = frame["dry_group"].value_counts().reindex(group_order, fill_value=0)
    return {
        "input": str(candidate_path),
        "output": str(out_path),
        "candidate_count": int(len(frame)),
        "groups": {group: int(count) for group, count in counts.items()},
        "warning": "Dry evidence orders parental sequencing; it does not establish causality. Rescue candidates remain separate from the knowledge score.",
    }


def rank_ledger(ledger_arg: str, out_arg: str, research_arg: str | None = None, go_arg: str | None = None) -> dict[str, object]:
    ledger_path, out_path = runtime_path(ledger_arg), runtime_path(out_arg)
    frame = pd.read_csv(ledger_path, sep="\t", dtype=str).fillna("")
    required = {"gene_id", "ortholog_ids", "function_evidence_source", "parental_variant", "variant_class", "f2_cosegregation_status"}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"ledger is missing required columns: {', '.join(missing)}")
    frame["research_article_count"] = 0
    frame["go_annotation_count"] = 0
    if research_arg:
        research = pd.read_csv(runtime_path(research_arg), sep="\t", dtype=str).fillna("")
        if not {"gene_id", "europe_pmc_result_count"}.issubset(research.columns):
            raise ValueError("research summary must contain gene_id and europe_pmc_result_count")
        research = research[["gene_id", "europe_pmc_result_count"]].copy()
        research["research_article_count"] = pd.to_numeric(research["europe_pmc_result_count"], errors="coerce").fillna(0).astype(int)
        frame = frame.merge(research[["gene_id", "research_article_count"]], on="gene_id", how="left", suffixes=("", "_api"))
        frame["research_article_count"] = pd.to_numeric(frame["research_article_count_api"], errors="coerce").fillna(frame["research_article_count"]).astype(int)
        frame = frame.drop(columns=["research_article_count_api"])
    if go_arg:
        go_evidence = pd.read_csv(runtime_path(go_arg), sep="\t", dtype=str).fillna("")
        if not {"gene_id", "go_id", "source"}.issubset(go_evidence.columns):
            raise ValueError("GO evidence must contain gene_id, go_id, and source")
        accepted = go_evidence[(go_evidence["go_id"].str.strip() != "") & (go_evidence["source"].str.strip() != "")]
        counts = accepted.groupby("gene_id").size().rename("go_annotation_count").reset_index()
        frame = frame.drop(columns=["go_annotation_count"]).merge(counts, on="gene_id", how="left")
        frame["go_annotation_count"] = frame["go_annotation_count"].fillna(0).astype(int)

    def score(row: pd.Series) -> int:
        total = 0
        if nonempty(row["ortholog_ids"]):
            total += 1
        if nonempty(row["function_evidence_source"]) and "required" not in row["function_evidence_source"].lower():
            total += 2
        total += VARIANT_POINTS.get(row["variant_class"].strip().lower(), 0)
        if nonempty(row["parental_variant"]) and "pending" not in row["parental_variant"].lower():
            total += 1
        if int(row["research_article_count"]) > 0:
            total += 1
        if int(row["go_annotation_count"]) > 0:
            total += 1
        return total

    frame["deterministic_priority_score"] = frame.apply(score, axis=1)
    frame["program_rank"] = frame["deterministic_priority_score"].rank(method="dense", ascending=False).astype(int)
    frame = frame.sort_values(["program_rank", "gene_id"], kind="stable")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out_path, sep="\t", index=False)
    return {
        "input": str(ledger_path),
        "output": str(out_path),
        "candidate_count": int(len(frame)),
        "warning": "F2 cosegregation is intentionally excluded from this priority score and is kept as a separate Wet-validation column, not a ranking input. Literature and GO presence add at most one point each. This score prioritizes assays; it does not establish causality.",
    }


def cosegregation(input_arg: str, out_arg: str, genotype_column: str, phenotype_column: str) -> dict[str, object]:
    input_path, out_path = runtime_path(input_arg), runtime_path(out_arg)
    frame = pd.read_csv(input_path, sep="\t", dtype=str).fillna("")
    missing = [column for column in [genotype_column, phenotype_column] if column not in frame.columns]
    if missing:
        raise ValueError(f"input is missing columns: {', '.join(missing)}")
    work = pd.DataFrame(
        {
            "_geno": frame[genotype_column].astype(str).str.strip(),
            "_pheno": frame[phenotype_column].astype(str).str.strip(),
        }
    )
    present = work[(work["_geno"] != "") & (work["_pheno"] != "")].copy()
    numeric = pd.to_numeric(present["_pheno"], errors="coerce")
    numeric_fraction = float(numeric.notna().mean()) if len(present) else 0.0
    quantitative = numeric_fraction >= 0.8
    rows: list[dict[str, object]] = []
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if quantitative:
        valid = present.assign(_value=numeric)
        valid = valid[valid["_value"].notna()]
        for genotype, group in valid.groupby("_geno", sort=True):
            numbers = [float(value) for value in group["_value"]]
            rows.append(
                {
                    "genotype": genotype,
                    "n": len(numbers),
                    "mean": round(statistics.mean(numbers), 6),
                    "median": round(statistics.median(numbers), 6),
                    "sd": round(statistics.stdev(numbers), 6) if len(numbers) > 1 else "",
                    "min": round(min(numbers), 6),
                    "max": round(max(numbers), 6),
                }
            )
        pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)
        means = [float(row["mean"]) for row in rows]
        result = {
            "input": str(input_path),
            "output": str(out_path),
            "mode": "quantitative",
            "genotype_groups": len(rows),
            "phenotyped_individuals": int(len(valid)),
            "mean_range": round(max(means) - min(means), 6) if means else None,
            "warning": "Descriptive result only. Use F2:3, blocks, environment, and recombinant genotypes before causal interpretation.",
        }
    else:
        valid = present
        for genotype, group in valid.groupby("_geno", sort=True):
            n_genotype = len(group)
            for category, count in sorted(group["_pheno"].value_counts().items()):
                rows.append(
                    {
                        "genotype": genotype,
                        "phenotype_category": category,
                        "n_in_category": int(count),
                        "n_genotype": n_genotype,
                        "proportion_within_genotype": round(int(count) / n_genotype, 6),
                    }
                )
        pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)
        result = {
            "input": str(input_path),
            "output": str(out_path),
            "mode": "categorical",
            "genotype_groups": int(valid["_geno"].nunique()),
            "phenotype_categories": int(valid["_pheno"].nunique()),
            "phenotyped_individuals": int(len(valid)),
            "warning": "Descriptive categorical counts only; no association test is computed. Use recombinant genotypes, F2:3, and independent evidence before causal interpretation.",
        }
    (out_path.with_suffix(".json")).write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def assay_plan(
    candidate_arg: str,
    out_arg: str,
    top_n: int,
    genome_arg: str | None = None,
    gff_arg: str | None = None,
    db_args: list[str] | None = None,
    primerblast_root_arg: str | None = None,
) -> dict[str, object]:
    candidate_path, out_path = runtime_path(candidate_arg), runtime_path(out_arg)
    frame = pd.read_csv(candidate_path, sep="\t", dtype=str).fillna("")
    if "gene_id" not in frame.columns:
        raise ValueError("candidate table must contain gene_id")
    if "program_rank" in frame.columns:
        frame["program_rank_numeric"] = pd.to_numeric(frame["program_rank"], errors="coerce").fillna(999999)
        frame = frame.sort_values(["program_rank_numeric", "gene_id"], kind="stable")
    selected = frame.head(top_n).copy()
    selected["assay_intent"] = "parental PCR/Sanger then marker conversion"
    selected["primerblast_oss_hint"] = "Use candidate_genomic.fa or candidate_cds.fa, then require parent/F1 PCR validation before F2 scale-out."
    out_path.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(out_path, sep="\t", index=False)
    result: dict[str, object] = {"input": str(candidate_path), "output": str(out_path), "selected_count": int(len(selected))}
    if genome_arg and gff_arg and db_args and primerblast_root_arg:
        genome_path = runtime_path(genome_arg)
        gff_path = runtime_path(gff_arg)
        primerblast_root = runtime_path(primerblast_root_arg)
        db_paths = [str(runtime_path(value)) for value in db_args]
        command_file = out_path.parent / "run-primerblast-oss-assays.sh"
        output_dir = out_path.parent / "primerblast_oss"
        db_flags = " ".join(f"--db {shlex.quote(value)}" for value in db_paths)
        lines = [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            f"mkdir -p {shlex.quote(str(output_dir))}",
        ]
        for gene_id in selected["gene_id"]:
            output = output_dir / f"{safe_id(str(gene_id))}.json"
            lines.append(
                "PYTHONPATH={root} python3 -m primerblast_oss assay --gene {gene} --genome {genome} --gff3 {gff3} {dbs} --format json --out {out}".format(
                    root=shlex.quote(str(primerblast_root)),
                    gene=shlex.quote(str(gene_id)),
                    genome=shlex.quote(str(genome_path)),
                    gff3=shlex.quote(str(gff_path)),
                    dbs=db_flags,
                    out=shlex.quote(str(output)),
                )
            )
        command_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
        command_file.chmod(0o755)
        result["primerblast_oss_script"] = str(command_file)
    else:
        result["primerblast_oss_script"] = None
        result["next_step"] = "Pass --genome, --gff3, --db, and --primerblast-root to generate executable primerblast-oss assay commands."
    return result


def _tool_version(path: str) -> str:
    """Best-effort external tool version string (first non-empty output line)."""
    for flag in ("-version", "--version", "-h"):
        try:
            completed = subprocess.run([path, flag], capture_output=True, text=True, timeout=10, check=False)
        except (OSError, subprocess.TimeoutExpired):
            continue
        for line in (completed.stdout + completed.stderr).splitlines():
            if line.strip():
                return line.strip()
    return "unknown"


def _package_versions(distributions: list[str]) -> dict[str, str | None]:
    versions: dict[str, str | None] = {}
    for dist in distributions:
        try:
            versions[dist] = importlib.metadata.version(dist)
        except importlib.metadata.PackageNotFoundError:
            versions[dist] = None
    return versions


def doctor() -> dict[str, object]:
    command_names = ["python3", "blastn", "makeblastdb", "primer3_core", "minimap2", "hmmscan", "emapper.py", "Rscript", "bcftools", "snakemake"]
    commands = {name: shutil.which(name) for name in command_names}
    modules = {}
    for module in ["pysam", "Bio", "pandas", "yaml"]:
        try:
            __import__(module)
            modules[module] = True
        except ImportError:
            modules[module] = False
    primerblast_root = os.environ.get("PRIMERBLAST_OSS_ROOT", "")
    primerblast_module = Path(primerblast_root) / "primerblast_oss" / "__init__.py" if primerblast_root else None
    r_packages = {"qtl2": False}
    if commands.get("Rscript"):
        check = run_tool(
            [commands["Rscript"], "-e", "quit(status=ifelse(requireNamespace('qtl2', quietly=TRUE), 0, 1))"],
            capture_output=True,
            check=False,
        )
        r_packages["qtl2"] = check.returncode == 0
    versioned_tools = ["blastn", "makeblastdb", "primer3_core", "minimap2", "hmmscan", "bcftools"]
    tool_versions = {name: _tool_version(commands[name]) for name in versioned_tools if commands.get(name)}
    return {
        "platform": os.uname().sysname if hasattr(os, "uname") else os.name,
        "python_version": platform.python_version(),
        "python_executable": sys.executable,
        "package_versions": _package_versions(["biopython", "pandas", "pysam", "PyYAML", "candidate-qtl"]),
        "tool_versions": tool_versions,
        "commands": commands,
        "python_modules": modules,
        "r_packages": r_packages,
        "ready_for_mvp": all(commands[name] for name in ["blastn", "makeblastdb", "primer3_core"]) and all(modules.values()),
        # pysam has no native Windows build. Without it the table-only commands
        # still run; the sequence-I/O ones do not. Report the split here so it is
        # discoverable before a command fails.
        "pysam_available": pysam_available(),
        "pysam_dependent_commands": list(PYSAM_DEPENDENT_COMMANDS),
        "optional_missing": [name for name in ["bcftools", "snakemake"] if not commands[name]],
        "primerblast_oss_root": primerblast_root or None,
        "primerblast_oss_ready": bool(primerblast_module and primerblast_module.exists()),
    }


def subset_proteins(candidate_arg: str, protein_arg: str, out_arg: str) -> dict[str, object]:
    candidate_path, protein_path, out_path = runtime_path(candidate_arg), runtime_path(protein_arg), runtime_path(out_arg)
    candidates = pd.read_csv(candidate_path, sep="\t", dtype=str).fillna("")
    if "gene_id" not in candidates.columns:
        raise ValueError("candidate table must contain gene_id")
    candidate_ids = set(candidates["gene_id"])
    transcript_to_gene = {
        str(row["representative_transcript"]): str(row["gene_id"])
        for _, row in candidates.iterrows()
        if nonempty(row.get("representative_transcript", ""))
    }
    selected_by_gene: dict[str, object] = {}
    with protein_path.open("r", encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            aliases = [record.id]
            aliases.extend(token for token in record.id.split("|") if token)
            gene_id = next((alias for alias in aliases if alias in candidate_ids), None)
            if not gene_id:
                gene_id = next((transcript_to_gene[alias] for alias in aliases if alias in transcript_to_gene), None)
            if gene_id and gene_id not in selected_by_gene:
                normalized = record[:]
                normalized.id = gene_id
                normalized.name = gene_id
                normalized.description = f"source_id={record.id}"
                selected_by_gene[gene_id] = normalized
    records = [selected_by_gene[gene_id] for gene_id in candidates["gene_id"] if gene_id in selected_by_gene]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, str(out_path), "fasta")
    missing = sorted(candidate_ids.difference(selected_by_gene))
    return {
        "input_candidates": str(candidate_path),
        "input_proteins": str(protein_path),
        "output": str(out_path),
        "selected_proteins": len(records),
        "unmatched_candidate_ids": missing,
        "id_normalization": "Exact gene/transcript IDs and pipe-delimited FASTA ID tokens are normalized to gene_id.",
    }


def annotation_plan(
    candidate_arg: str,
    out_arg: str,
    protein_arg: str | None = None,
    blastp_db_args: list[str] | None = None,
    diamond_db_args: list[str] | None = None,
    hmmer_db_arg: str | None = None,
    eggnog_data_dir_arg: str | None = None,
    eggnog_bin: str = "emapper.py",
    num_threads: int = 4,
) -> dict[str, object]:
    candidate_path, out_dir = runtime_path(candidate_arg), runtime_path(out_arg)
    if num_threads < 1:
        raise ValueError("num_threads must be at least 1")
    out_dir.mkdir(parents=True, exist_ok=True)
    candidates = pd.read_csv(candidate_path, sep="\t", dtype=str).fillna("")
    if "gene_id" not in candidates.columns:
        raise ValueError("candidate table must contain gene_id")
    gene_ids = list(candidates["gene_id"])
    go_template = out_dir / "go_evidence.tsv"
    pd.DataFrame(
        {
            "gene_id": gene_ids,
            "go_id": "",
            "go_term": "",
            "source": "",
            "evidence_type": "",
            "ortholog_or_hit": "",
            "source_url_or_accession": "",
            "note": "",
        }
    ).to_csv(go_template, sep="\t", index=False)
    result: dict[str, object] = {
        "candidate_count": len(gene_ids),
        "go_evidence_template": str(go_template),
        "script": None,
        "warning": "Use GO terms as candidate evidence, not as proof. Do not claim enrichment from one narrow QTL interval.",
    }
    if not protein_arg:
        result["next_step"] = "Provide --protein-fasta to create a candidate-only protein FASTA and local homology/domain script."
        return result
    protein_subset = out_dir / "candidate_proteins.fa"
    subset = subset_proteins(candidate_arg, protein_arg, str(protein_subset))
    lines = ["#!/usr/bin/env bash", "set -euo pipefail", f"OUT={shlex.quote(str(out_dir))}", "mkdir -p \"$OUT\""]
    for index, db in enumerate(blastp_db_args or [], start=1):
        lines.append(
            "blastp -query {query} -db {db} -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle' -max_target_seqs 25 -num_threads {threads} -out \"$OUT/blastp_hits_{index}.tsv\"".format(
                query=shlex.quote(str(protein_subset)),
                db=shlex.quote(str(runtime_path(db))),
                index=index,
                threads=num_threads,
            )
        )
    for index, db in enumerate(diamond_db_args or [], start=1):
        lines.append(
            "diamond blastp --query {query} --db {db} --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle --max-target-seqs 25 --threads {threads} --out \"$OUT/diamond_hits_{index}.tsv\"".format(
                query=shlex.quote(str(protein_subset)),
                db=shlex.quote(str(runtime_path(db))),
                index=index,
                threads=num_threads,
            )
        )
    if hmmer_db_arg:
        lines.append(
            "hmmscan --cpu {threads} --tblout \"$OUT/hmmer_hits.tbl\" {db} {query}".format(
                db=shlex.quote(str(runtime_path(hmmer_db_arg))),
                query=shlex.quote(str(protein_subset)),
                threads=num_threads,
            )
        )
    if eggnog_data_dir_arg:
        lines.append(
            "{bin} -i {query} -m diamond --data_dir {data_dir} --output candidate_eggnog --output_dir \"$OUT\" --cpu {threads}".format(
                bin=shlex.quote(eggnog_bin),
                query=shlex.quote(str(protein_subset)),
                data_dir=shlex.quote(str(runtime_path(eggnog_data_dir_arg))),
                threads=num_threads,
            )
        )
    if len(lines) == 4:
        lines.append("echo 'No homology or HMM database was supplied; candidate_proteins.fa is ready for selected local or remote annotation.'")
    script_path = out_dir / "run_local_annotation.sh"
    script_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    script_path.chmod(0o755)
    result.update({"script": str(script_path), "protein_subset": subset})
    return result


def europe_pmc_search(query: str, out_arg: str, limit: int) -> dict[str, object]:
    out_path = runtime_path(out_arg)
    if limit < 1 or limit > 100:
        raise ValueError("limit must be between 1 and 100")
    params = urlencode({"query": query, "format": "json", "pageSize": limit})
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?{params}"
    request = Request(url, headers={"User-Agent": "candidate-qtl/0.1 (read-only literature search)"})
    try:
        with urlopen(request, timeout=30) as response:
            payload = json.load(response)
    except (urllib.error.URLError, TimeoutError, OSError) as error:
        kind = "timeout" if isinstance(error, TimeoutError) else "network"
        raise ExternalDependencyError(kind, f"Europe PMC search failed: {error}") from error
    records = []
    for item in payload.get("resultList", {}).get("result", []):
        records.append(
            {
                "id": item.get("id", ""),
                "source": item.get("source", ""),
                "title": item.get("title", ""),
                "authors": item.get("authorString", ""),
                "journal": item.get("journalTitle", ""),
                "year": item.get("pubYear", ""),
                "doi": item.get("doi", ""),
                "pmid": item.get("pmid", ""),
                "cited_by_count": item.get("citedByCount", ""),
            }
        )
    result = {
        "provider": "Europe PMC",
        "query": query,
        "url": url,
        "retrieved_at": datetime.now(UTC).isoformat(),
        "result_count": len(records),
        "records": records,
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    return {"output": str(out_path), "query": query, "result_count": len(records), "url": url}


def research_packet(
    candidate_arg: str,
    out_arg: str,
    species: str,
    trait: str,
    top_n: int,
    fetch_literature: bool,
    literature_limit: int,
) -> dict[str, object]:
    candidate_path, out_dir = runtime_path(candidate_arg), runtime_path(out_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(candidate_path, sep="\t", dtype=str).fillna("")
    if "gene_id" not in frame.columns:
        raise ValueError("candidate table must contain gene_id")
    if "program_rank" in frame.columns:
        frame["_rank"] = pd.to_numeric(frame["program_rank"], errors="coerce").fillna(999999)
        frame = frame.sort_values(["_rank", "gene_id"], kind="stable")
    selected = frame.head(top_n).copy()
    packet_rows = []
    prompts = ["# Candidate-gene research prompts", "", "Each section is independent. Require direct URLs, DOI, PMID, or database accessions. Do not claim causality from functional plausibility or co-segregation alone."]
    all_literature = []
    for _, row in selected.iterrows():
        gene_id = row["gene_id"]
        orthologs = row.get("ortholog_ids", "")
        query_terms = " ".join(term for term in [gene_id, orthologs.replace(";", " "), species, trait] if term)
        packet_rows.append(
            {
                "gene_id": gene_id,
                "species": species,
                "trait": trait,
                "literature_query": query_terms,
                "parental_variant": row.get("parental_variant", ""),
                "variant_class": row.get("variant_class", ""),
                "function_evidence_source": row.get("function_evidence_source", ""),
                "uniprot_search_url": "https://rest.uniprot.org/uniprotkb/search?" + urlencode({"query": query_terms, "format": "json", "size": 10}),
            }
        )
        prompts.extend(
            [
                "",
                f"## {gene_id}",
                f"Species: {species}",
                f"Trait: {trait}",
                f"Ortholog IDs: {orthologs or 'not supplied'}",
                f"Observed parental variant: {row.get('parental_variant', '') or 'not supplied'} ({row.get('variant_class', '') or 'not classified'})",
                "Research task: return direct evidence, ortholog evidence, and contradictory evidence separately. For every claim provide a DOI, PMID, or stable database URL. State whether evidence is from this species or inferred from an ortholog. List GO/domain evidence separately from phenotype evidence. Do not state that this gene is causal without independent genetic or functional evidence.",
            ]
        )
        if fetch_literature:
            literature_path = out_dir / f"{safe_id(gene_id)}_europe_pmc.json"
            europe_pmc_search(query_terms, str(literature_path), literature_limit)
            all_literature.append(str(literature_path))
    packet_path = out_dir / "research_packet.tsv"
    selected_prompts = out_dir / "ai_research_prompts.md"
    pd.DataFrame(packet_rows).to_csv(packet_path, sep="\t", index=False)
    selected_prompts.write_text("\n".join(prompts) + "\n", encoding="utf-8")
    return {
        "packet": str(packet_path),
        "prompts": str(selected_prompts),
        "candidate_count": len(packet_rows),
        "literature_files": all_literature,
        "ai_optional": True,
    }


def candidate_research(
    candidate_arg: str,
    out_arg: str,
    species: str,
    trait: str,
    top_n: int,
    literature_limit: int,
) -> dict[str, object]:
    candidate_path, out_dir = runtime_path(candidate_arg), runtime_path(out_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(candidate_path, sep="\t", dtype=str).fillna("")
    if "gene_id" not in frame.columns:
        raise ValueError("candidate table must contain gene_id")
    if "program_rank" in frame.columns:
        frame["_rank"] = pd.to_numeric(frame["program_rank"], errors="coerce").fillna(999999)
        frame = frame.sort_values(["_rank", "gene_id"], kind="stable")
    selected = frame.head(top_n)
    summary_rows = []
    evidence_rows = []
    for _, row in selected.iterrows():
        gene_id = row["gene_id"]
        orthologs = row.get("ortholog_ids", "")
        query_terms = " ".join(term for term in [gene_id, orthologs.replace(";", " "), species, trait] if term)
        article_path = out_dir / "literature" / f"{safe_id(gene_id)}_europe_pmc.json"
        article_result = europe_pmc_search(query_terms, str(article_path), literature_limit)
        payload = json.loads(article_path.read_text(encoding="utf-8"))
        summary_rows.append(
            {
                "gene_id": gene_id,
                "species": species,
                "trait": trait,
                "literature_query": query_terms,
                "europe_pmc_result_count": article_result["result_count"],
                "literature_json": str(article_path),
                "uniprot_search_url": "https://rest.uniprot.org/uniprotkb/search?" + urlencode({"query": query_terms, "format": "json", "size": 10}),
            }
        )
        for article in payload["records"]:
            evidence_rows.append({"gene_id": gene_id, "literature_query": query_terms, **article})
    summary_path = out_dir / "research_summary.tsv"
    evidence_path = out_dir / "literature_evidence.tsv"
    pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)
    pd.DataFrame(evidence_rows).to_csv(evidence_path, sep="\t", index=False)
    return {
        "summary": str(summary_path),
        "evidence": str(evidence_path),
        "candidate_count": len(summary_rows),
        "api_provider": "Europe PMC",
        "ai_used": False,
    }
