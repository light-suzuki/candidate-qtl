from __future__ import annotations

import json
import tempfile
from time import perf_counter
from collections import Counter, defaultdict
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd
from Bio.Seq import Seq

from . import schemas
from .core import (
    file_fingerprint,
    fingerprint_is_weak,
    iter_gff_rows,
    parse_interval,
    parse_seqid_map,
    runtime_path,
    write_fasta,
)
from .optional_deps import require_pysam


def _fai_layout(path: Path) -> list[tuple[str, int, int, int, int]]:
    """Read the five uncompressed FASTA index fields without trusting offsets."""
    layout = []
    names = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) != 5 or not fields[0] or fields[0] in names:
            raise ValueError("invalid shared FASTA index layout")
        try:
            values = tuple(int(value) for value in fields[1:])
        except ValueError as exc:
            raise ValueError("invalid shared FASTA index layout") from exc
        if any(value < 0 for value in values) or values[3] < values[2]:
            raise ValueError("invalid shared FASTA index layout")
        names.add(fields[0])
        layout.append((fields[0], *values))
    if not layout:
        raise ValueError("empty shared FASTA index layout")
    return layout


def _validate_readonly_fai(fasta_path: Path, index_path: Path, pysam) -> dict[str, object]:
    """Independently scan the FASTA; never create or repair its shared sidecar."""
    started = perf_counter()
    supplied = _fai_layout(index_path)
    with tempfile.TemporaryDirectory(prefix="cqtl-fai-check-") as temporary:
        derived_path = Path(temporary) / "reference.fai"
        pysam.faidx("--fai-idx", str(derived_path), str(fasta_path))
        if supplied != _fai_layout(derived_path):
            raise ValueError("shared FASTA/.fai layout mismatch; supply an index for this exact FASTA")
    return {"method": "independent_faidx_layout", "contig_count": len(supplied),
            "elapsed_seconds": perf_counter() - started}


if TYPE_CHECKING:  # annotation only; pysam is imported at the point of use
    import pysam


# GFF3 dialects vary in which feature types wrap coding models. These sets keep
# gene/transcript recognition in one place so placeholder synthesis and the main
# extractor agree.
GENE_TYPES = {
    "gene",
    "protein_coding_gene",
    "pseudogene",
    "ncRNA_gene",
    "transposable_element_gene",
}
TRANSCRIPT_TYPES = {
    "mRNA",
    "transcript",
    "primary_transcript",
    "pseudogenic_transcript",
    "mRNA_TE_gene",
}


def _is_transcript_type(feature_type: str) -> bool:
    return feature_type in TRANSCRIPT_TYPES or feature_type.endswith(("RNA", "_primary_transcript"))


def _is_transcript_feature(row: dict[str, object]) -> bool:
    return _is_transcript_type(row["type"]) and not row.get("_rna_product", False)


def _first_attr(attrs: dict[str, object], *keys: str) -> str:
    for key in keys:
        value = str(attrs.get(key, "")).strip()
        if value:
            return value
    return ""


def _synthetic_row(feature_type: str, seqid: str, start: int, end: int, strand: str, attrs: dict[str, object]) -> dict[str, object]:
    return {
        "seqid": seqid,
        "type": feature_type,
        "start": int(start),
        "end": int(end),
        "strand": strand,
        "phase": ".",
        "attrs": attrs,
    }


def _normalize_features(feature_rows: list[dict[str, object]]) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Synthesize placeholder gene/transcript features for gene-less dialects.

    Some GFF3 files attach CDS/exon directly to a gene, name a Parent gene that
    has no feature line, or ship transcripts with no gene at all. Rather than
    silently dropping those loci, this step manufactures the missing parent
    features so the phase-aware extractor sees a consistent gene -> transcript ->
    CDS hierarchy. It also collects QC counters (duplicate IDs, orphan features,
    unresolved parents, placeholder counts).
    """

    gene_ids = {
        _first_attr(row["attrs"], "ID")
        for row in feature_rows
        if row["type"] in GENE_TYPES and _first_attr(row["attrs"], "ID")
    }
    by_id: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in feature_rows:
        if identifier := _first_attr(row["attrs"], "ID"):
            by_id[identifier].append(row)
    cross_type_ids = {identifier for identifier, rows in by_id.items()
                      if len({row["type"] for row in rows}) > 1}
    transcript_like_ids = {identifier for identifier, rows in by_id.items()
                           if any(_is_transcript_type(row["type"]) for row in rows)}

    def relations(row):
        return {parent for key in ("Parent", "Derives_from")
                for parent in str(row["attrs"].get(key, "")).split(",") if parent}

    # A primary transcript is a transcript; an RNA derived from it is a product.
    # Keep the original IDs and relationships; products never invent gene parents.
    for row in feature_rows:
        kind = row["type"]
        if kind.endswith("RNA") and kind not in TRANSCRIPT_TYPES:
            row["_rna_product"] = bool(row["attrs"].get("Derives_from")) or bool(relations(row) & transcript_like_ids)
    transcript_rows = [row for row in feature_rows if _is_transcript_feature(row)]
    transcript_ids = {_first_attr(row["attrs"], "ID") for row in transcript_rows if _first_attr(row["attrs"], "ID")}

    # Duplicate gene/transcript IDs (multi-line CDS legitimately share an ID and
    # are excluded from this check).
    parent_level_ids: list[str] = []
    for row in feature_rows:
        if row["type"] in GENE_TYPES or _is_transcript_feature(row):
            identifier = _first_attr(row["attrs"], "ID")
            if identifier:
                parent_level_ids.append(identifier)
    parent_id_counts = Counter(parent_level_ids)
    duplicate_ids = sorted(identifier for identifier, count in parent_id_counts.items() if count > 1)

    synthetic_rows: list[dict[str, object]] = []
    placeholder_genes: dict[str, dict[str, object]] = {}
    unresolved_parents = 0
    cyclic_ids: set[str] = set()
    ambiguous_ids: set[str] = set()
    missing_parent_ids: set[str] = set()
    inferred_gene_by_transcript: dict[str, str] = {}

    def ancestors(identifier: str, visiting: frozenset[str] = frozenset()) -> tuple[set[str], bool]:
        if identifier in cross_type_ids:
            ambiguous_ids.add(identifier)
            return set(), True
        if identifier in gene_ids or identifier in placeholder_genes:
            return {identifier}, False
        if identifier in inferred_gene_by_transcript:
            return {inferred_gene_by_transcript[identifier]}, False
        if identifier in visiting:
            cyclic_ids.add(identifier)
            return set(), True
        rows = by_id.get(identifier, [])
        if not rows:
            missing_parent_ids.add(identifier)
            return set(), True
        parents = set().union(*(relations(row) for row in rows))
        genes: set[str] = set()
        unresolved = not parents
        for parent in sorted(parents):
            found, failed = ancestors(parent, visiting | {identifier})
            genes.update(found)
            unresolved |= failed
        if len(genes) > 1:
            ambiguous_ids.add(identifier)
            unresolved = True
        return genes, unresolved

    def touch_gene(gene_id: str, seqid: str, start: int, end: int, strand: str) -> None:
        span = placeholder_genes.get(gene_id)
        if span is None:
            placeholder_genes[gene_id] = {"seqid": seqid, "start": int(start), "end": int(end), "strand": strand}
        else:
            span["start"] = min(span["start"], int(start))
            span["end"] = max(span["end"], int(end))

    # Only genuinely missing gene features receive placeholders. A known
    # transcript/product/other feature ID must never be recycled as a gene ID.
    initial_bindings = [ancestors(_first_attr(row["attrs"], "ID")) for row in transcript_rows]
    for row, (genes, failed) in zip(transcript_rows, initial_bindings):
        parents = relations(row)
        transcript_id = _first_attr(row["attrs"], "ID") or f"tx:{row['seqid']}:{row['start']}:{row['end']}"
        if len(genes) == 1 and not failed:
            continue
        unresolved_parents += 1
        if parents and (len(parents) != 1 or any(parent in by_id for parent in parents)):
            continue
        gene_id = next(iter(parents)) if parents else f"gene:{transcript_id}"
        if gene_id in by_id or transcript_id in cross_type_ids:
            ambiguous_ids.add(transcript_id)
            continue
        touch_gene(gene_id, str(row["seqid"]), int(row["start"]), int(row["end"]), str(row["strand"]))
        inferred_gene_by_transcript[transcript_id] = gene_id

    # Resolve after all placeholders exist, so child-before-parent order is safe.
    for row in transcript_rows:
        genes, failed = ancestors(_first_attr(row["attrs"], "ID"))
        row["_resolved_gene_parents"] = sorted(genes) if len(genes) == 1 and not failed else []
        if len(genes) == 1 and not failed and (gene_id := next(iter(genes))) in placeholder_genes:
            touch_gene(gene_id, str(row["seqid"]), int(row["start"]), int(row["end"]), str(row["strand"]))
    unresolved_products = []
    for row in feature_rows:
        if row.get("_rna_product"):
            identifier = _first_attr(row["attrs"], "ID")
            genes, failed = ancestors(identifier)
            if len(genes) != 1 or failed:
                unresolved_products.append(identifier)

    # CDS whose parent is neither a known transcript nor a known/placeholder gene
    # get both a placeholder transcript and gene, grouped by the CDS Parent.
    orphan_cds = 0
    cds_groups: dict[str, list[dict[str, object]]] = {}
    unresolved_cds_parent_ids: set[str] = set()
    for row in feature_rows:
        if row["type"] != "CDS":
            continue
        parents = [parent for parent in str(row["attrs"].get("Parent", "")).split(",") if parent]
        resolvable = any(
            parent in transcript_ids or parent in gene_ids or parent in placeholder_genes
            for parent in parents
        )
        if resolvable:
            continue
        orphan_cds += 1
        if len(parents) > 1 or any(parent in by_id for parent in parents):
            unresolved_cds_parent_ids.update(parents)
            continue
        key = parents[0] if parents else f"cds:{row['seqid']}:{row['start']}"
        if f"gene:{key}" in by_id:
            unresolved_cds_parent_ids.add(key)
            continue
        cds_groups.setdefault(key, []).append(row)

    for parent_id, rows in cds_groups.items():
        seqid = str(rows[0]["seqid"])
        strand = str(rows[0]["strand"])
        start = min(int(r["start"]) for r in rows)
        end = max(int(r["end"]) for r in rows)
        gene_id = f"gene:{parent_id}"
        synthetic_rows.append(
            _synthetic_row("mRNA", seqid, start, end, strand, {"ID": parent_id, "Parent": gene_id, "_placeholder": "cds_only"})
        )
        touch_gene(gene_id, seqid, start, end, strand)

    for gene_id, span in placeholder_genes.items():
        synthetic_rows.append(
            _synthetic_row(
                "gene",
                str(span["seqid"]),
                int(span["start"]),
                int(span["end"]),
                str(span["strand"]),
                {"ID": gene_id, "_placeholder": "synthesized_gene", "biotype": "inferred"},
            )
        )

    total_parent_features = len(transcript_rows) + sum(1 for row in feature_rows if row["type"] == "CDS")
    qc = {
        "duplicate_parent_ids": duplicate_ids,
        "cross_type_id_collisions": sorted(cross_type_ids),
        "cyclic_parent_ids": sorted(cyclic_ids),
        "ambiguous_hierarchy_ids": sorted(ambiguous_ids),
        "missing_parent_ids": sorted(missing_parent_ids),
        "rna_product_ids": sorted(_first_attr(row["attrs"], "ID") for row in feature_rows if row.get("_rna_product")),
        "unresolved_rna_product_ids": sorted(set(unresolved_products)),
        "unresolved_cds_parent_ids": sorted(unresolved_cds_parent_ids),
        "orphan_cds_features": orphan_cds,
        "transcripts_without_gene_parent": unresolved_parents,
        "placeholder_gene_count": len(placeholder_genes),
        "unmatched_parent_rate": round((unresolved_parents + orphan_cds) / total_parent_features, 4) if total_parent_features else 0.0,
    }
    return feature_rows + synthetic_rows, qc


def _load_validator_report(path_arg: str | None) -> dict[str, object] | None:
    """Attach an external GFF3 validator report (GenomeTools/AGAT) to the manifest.

    The report is stored by reference with a fingerprint and a short excerpt; it
    is not re-run here, only recorded so the extraction remains reproducible.
    """
    if not path_arg:
        return None
    path = runtime_path(path_arg)
    if not path.exists():
        raise ValueError(f"GFF QC report not found: {path}")
    excerpt = path.read_text(encoding="utf-8", errors="replace").splitlines()[:20]
    return {
        "path": str(path),
        "fingerprint": file_fingerprint(path, False),
        "excerpt": excerpt,
    }


def _phase_value(raw: object) -> int | None:
    text = str(raw).strip()
    if text in {"0", "1", "2"}:
        return int(text)
    return None


def _assemble_transcript_cds(
    fasta: pysam.FastaFile,
    fasta_seqid: str,
    parts: list[dict[str, object]],
    strand: str,
) -> dict[str, object]:
    """Assemble one CDS in transcript 5'->3' order and validate GFF3 phase.

    GFF3 phase is not an instruction to trim every CDS segment. It states where
    the next codon begins relative to the 5' end of that segment. Internal
    phases are therefore validated against cumulative CDS length; only a
    non-zero phase on the first, 5'-partial CDS is removed before translation.
    """
    ordered = sorted(
        parts,
        key=lambda row: (int(row["start"]), int(row["end"])),
        reverse=strand == "-",
    )
    raw_sequences: list[str] = []
    reported_phases: list[str] = []
    expected_phases: list[str] = []
    warnings: list[str] = []
    cumulative_translation_bases = 0
    first_phase = 0
    missing_phase = False

    for index, part in enumerate(ordered):
        sequence = fasta.fetch(
            fasta_seqid,
            int(part["start"]) - 1,
            int(part["end"]),
        ).upper()
        if strand == "-":
            sequence = str(Seq(sequence).reverse_complement())
        raw_sequences.append(sequence)

        reported = _phase_value(part.get("phase", "."))
        reported_phases.append("." if reported is None else str(reported))
        if reported is None:
            missing_phase = True
            expected_phases.append("unknown")
            warnings.append(f"segment_{index + 1}:missing_phase")
            if index == 0:
                first_phase = 0
            cumulative_translation_bases += len(sequence)
            continue

        if index == 0:
            first_phase = reported
            expected = reported
            cumulative_translation_bases += max(0, len(sequence) - first_phase)
        else:
            expected = (3 - (cumulative_translation_bases % 3)) % 3
            if reported != expected:
                warnings.append(
                    f"segment_{index + 1}:reported_phase_{reported}_expected_{expected}"
                )
            cumulative_translation_bases += len(sequence)
        expected_phases.append(str(expected))

    raw_cds = "".join(raw_sequences)
    cds_sequence = raw_cds[first_phase:] if raw_cds else ""
    inconsistent = any("expected_" in warning for warning in warnings)
    phase_usable = bool(raw_cds) and not missing_phase and not inconsistent
    if not ordered:
        phase_status = "no_cds"
    elif missing_phase:
        phase_status = "missing"
    elif inconsistent:
        phase_status = "inconsistent"
    elif first_phase:
        phase_status = "valid_partial_5prime"
    else:
        phase_status = "valid"

    protein = ""
    translation_status = "no_cds"
    if cds_sequence:
        if not phase_usable:
            translation_status = "not_translated_invalid_phase"
        elif len(cds_sequence) % 3 != 0:
            translation_status = "not_translated_length_not_multiple_of_3"
            warnings.append(f"translated_length_not_multiple_of_3:{len(cds_sequence)}")
        else:
            protein = str(Seq(cds_sequence).translate())
            translation_status = "translated"

    complete_orf = bool(
        protein
        and first_phase == 0
        and cds_sequence.startswith("ATG")
        and protein.endswith("*")
        and "*" not in protein[:-1]
    )
    if complete_orf:
        translation_status = "translated_complete_orf"
    elif protein:
        translation_status = "translated_partial_or_unverified_orf"

    return {
        "cds_sequence": cds_sequence,
        "protein_sequence": protein[:-1] if protein.endswith("*") else protein,
        "complete_orf": complete_orf,
        "phase_usable": phase_usable,
        "phase_status": phase_status,
        "first_cds_phase": first_phase,
        "reported_phases": ",".join(reported_phases),
        "expected_phases": ",".join(expected_phases),
        "phase_warnings": ";".join(warnings),
        "translation_status": translation_status,
        "cds_segment_count": len(ordered),
    }


def extract_candidates(
    fasta_arg: str,
    gff_arg: str,
    interval_arg: str,
    out_arg: str,
    seqid_map_arg: str | None = None,
    full_hash: bool = False,
    gff_qc_report_arg: str | None = None,
    *,
    read_only_reference: bool = False,
) -> dict[str, object]:
    """Extract reference-anchored QTL candidates with phase-aware CDS models."""
    fasta_path = runtime_path(fasta_arg)
    gff_path = runtime_path(gff_arg)
    out_dir = runtime_path(out_arg)
    chrom, interval_start, interval_end = parse_interval(interval_arg)
    seqid_map = parse_seqid_map(seqid_map_arg)
    fasta_chrom = seqid_map.get(chrom, chrom)
    # Check the inputs here rather than letting pysam surface a missing genome as
    # a SamtoolsError about a failed .fai build, which names the index instead of
    # the argument the caller got wrong.
    for label, path in (("--fasta", fasta_path), ("--gff3", gff_path)):
        if not path.exists():
            raise FileNotFoundError(f"{label} does not exist: {path}")
        if path.is_dir():
            raise IsADirectoryError(f"{label} must be a file, not a directory: {path}")
    index_path = Path(str(fasta_path) + ".fai")
    source_stats = {}
    index_fingerprint = None
    index_validation = None
    if read_only_reference:
        if any(path.resolve().is_relative_to(out_dir.resolve())
               for path in (fasta_path, gff_path, index_path)):
            raise ValueError("read-only extraction output must be separate from shared inputs")
        output_names = ("candidates.tsv", "transcripts.tsv", "id_map.tsv", "extract_manifest.json",
                        "candidate_genomic.fa", "candidate_cds.fa", "candidate_transcript_cds.fa",
                        "candidate_proteins.fa")
        if any((out_dir / name).exists() or (out_dir / name).is_symlink() for name in output_names):
            raise ValueError("read-only extraction requires fresh output artifacts")
        # Supplying an explicit existing index prevents pysam from creating a
        # sidecar beside a shared reference. Compressed references need a second
        # index and are not supported by this opt-in mode yet.
        with fasta_path.open("rb") as stream:
            if stream.read(2) == b"\x1f\x8b":
                raise ValueError("read-only reference requires an uncompressed FASTA")
        for path in (fasta_path, gff_path, index_path):
            if path.is_symlink() or not path.is_file():
                raise ValueError(f"read-only reference requires regular inputs and an existing .fai: {path}")
            stat = path.stat()
            source_stats[path] = (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
        index_fingerprint = file_fingerprint(index_path, True)
        index_validation = _validate_readonly_fai(fasta_path, index_path, require_pysam('extract'))
    out_dir.mkdir(parents=True, exist_ok=True)

    feature_rows, feature_qc = _normalize_features(list(iter_gff_rows(gff_path, chrom)))
    genes = [
        row
        for row in feature_rows
        if row["type"] in GENE_TYPES
        and not (int(row["end"]) < interval_start or int(row["start"]) > interval_end)
    ]
    gene_ids = {
        str(
            row["attrs"].get("ID")
            or row["attrs"].get("Name")
            or f"{chrom}_{row['start']}_{row['end']}"
        )
        for row in genes
    }
    gene_features = {
        str(
            row["attrs"].get("ID")
            or row["attrs"].get("Name")
            or f"{chrom}_{row['start']}_{row['end']}"
        ): row
        for row in genes
    }

    transcript_to_gene: dict[str, str] = {}
    transcript_features: dict[str, dict[str, object]] = {}
    for row in feature_rows:
        if not _is_transcript_feature(row):
            continue
        transcript_id = str(row["attrs"].get("ID", ""))
        parents = [
            parent
            for parent in row.get("_resolved_gene_parents", str(row["attrs"].get("Parent", "")).split(","))
            if parent in gene_ids
        ]
        if transcript_id and parents:
            transcript_to_gene[transcript_id] = parents[0]
            transcript_features[transcript_id] = row

    transcript_cds: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in feature_rows:
        if row["type"] != "CDS":
            continue
        for parent in str(row["attrs"].get("Parent", "")).split(","):
            if parent in transcript_to_gene:
                transcript_cds[parent].append(row)
            elif parent in gene_ids:
                transcript_to_gene.setdefault(parent, parent)
                transcript_features.setdefault(parent, gene_features[parent])
                transcript_cds[parent].append(row)

    pysam = require_pysam('extract')
    if read_only_reference:
        fasta = pysam.FastaFile(str(fasta_path), filepath_index=str(index_path))
    else:
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
        gene_id = str(
            attrs.get("ID")
            or attrs.get("Name")
            or f"{chrom}_{gene['start']}_{gene['end']}"
        )
        genomic_sequence = fasta.fetch(
            fasta_chrom,
            int(gene["start"]) - 1,
            int(gene["end"]),
        ).upper()
        if gene["strand"] == "-":
            genomic_sequence = str(Seq(genomic_sequence).reverse_complement())
        genomic_records.append((gene_id, genomic_sequence))

        isoforms: list[dict[str, object]] = []
        for transcript_id, parent_gene in transcript_to_gene.items():
            if parent_gene != gene_id:
                continue
            feature = transcript_features.get(transcript_id, gene)
            assembled = _assemble_transcript_cds(
                fasta,
                fasta_chrom,
                transcript_cds.get(transcript_id, []),
                str(gene["strand"]),
            )
            tag_text = " ".join(
                str(feature["attrs"].get(key, ""))
                for key in ("tag", "biotype", "transcript_type")
            ).casefold()
            canonical = any(
                term in tag_text
                for term in ("canonical", "representative", "mane_select")
            )
            isoform = {
                "gene_id": gene_id,
                "transcript_id": transcript_id,
                "transcript_start": feature["start"],
                "transcript_end": feature["end"],
                "strand": gene["strand"],
                "cds_length": len(str(assembled["cds_sequence"])),
                "protein_length": len(str(assembled["protein_sequence"])),
                "canonical_tag": canonical,
                **assembled,
            }
            isoforms.append(isoform)
            if assembled["cds_sequence"]:
                transcript_cds_records.append(
                    (transcript_id, str(assembled["cds_sequence"]))
                )

        isoforms.sort(
            key=lambda row: (
                -int(bool(row["phase_usable"])),
                -int(bool(row["canonical_tag"])),
                -int(bool(row["complete_orf"])),
                -int(row["cds_length"]),
                str(row["transcript_id"]),
            )
        )
        representative = isoforms[0] if isoforms else None
        if representative and representative["cds_sequence"]:
            representative_cds_records.append(
                (gene_id, str(representative["cds_sequence"]))
            )
        if representative and representative["protein_sequence"]:
            representative_protein_records.append(
                (gene_id, str(representative["protein_sequence"]))
            )

        for isoform in isoforms:
            is_representative = bool(
                representative
                and isoform["transcript_id"] == representative["transcript_id"]
            )
            reason = ""
            if is_representative:
                if not isoform["cds_segment_count"]:
                    reason = "best available transcript; no annotated CDS"
                elif not isoform["phase_usable"]:
                    reason = "best available transcript; CDS phase invalid or missing"
                elif isoform["canonical_tag"]:
                    reason = "phase-valid canonical annotation tag"
                elif isoform["complete_orf"]:
                    reason = "phase-valid longest complete ORF"
                else:
                    reason = "phase-valid longest annotated CDS"
            transcript_rows.append(
                {
                    key: value
                    for key, value in isoform.items()
                    if key not in {"cds_sequence", "protein_sequence"}
                }
                | {
                    "is_representative": is_representative,
                    "selection_reason": reason,
                }
            )

        annotation_quality = "no_cds"
        if representative and representative["cds_segment_count"]:
            if not representative["phase_usable"]:
                annotation_quality = "invalid_cds_phase"
            elif representative["complete_orf"]:
                annotation_quality = "complete_orf"
            else:
                annotation_quality = "partial_or_unverified_orf"
        gene_biotype = _first_attr(attrs, "biotype", "gene_biotype", "gene_type") or (
            gene["type"] if gene["type"] != "gene" else "protein_coding"
        )
        is_placeholder = bool(_first_attr(attrs, "_placeholder"))
        if annotation_quality == "no_cds" and gene_biotype not in {"protein_coding", "inferred"}:
            annotation_quality = "noncoding_or_pseudogene"
        table_rows.append(
            {
                "gene_id": gene_id,
                "gene_name": attrs.get("Name", ""),
                "gene_biotype": gene_biotype,
                "is_placeholder_gene": is_placeholder,
                "description": attrs.get(
                    "Note",
                    attrs.get("description", ""),
                ),
                "dbxref": attrs.get("Dbxref", ""),
                "ontology_terms": attrs.get("Ontology_term", ""),
                "interval": interval_arg,
                "coordinate_basis": "reference_gff3",
                "gff_seqid": chrom,
                "fasta_seqid": fasta_chrom,
                "chrom": chrom,
                "start": gene["start"],
                "end": gene["end"],
                "strand": gene["strand"],
                "gene_length": int(gene["end"]) - int(gene["start"]) + 1,
                "transcript_count": len(isoforms),
                "representative_transcript": (
                    representative["transcript_id"] if representative else ""
                ),
                "representative_reason": next(
                    (
                        row["selection_reason"]
                        for row in transcript_rows
                        if row["gene_id"] == gene_id and row["is_representative"]
                    ),
                    "",
                ),
                "cds_length": int(representative["cds_length"])
                if representative
                else 0,
                "protein_length": int(representative["protein_length"])
                if representative
                else 0,
                "cds_phase_status": representative["phase_status"]
                if representative
                else "no_cds",
                "translation_status": representative["translation_status"]
                if representative
                else "no_cds",
                "annotation_quality": annotation_quality,
                "overlaps_qtl_boundary": int(gene["start"]) < interval_start
                or int(gene["end"]) > interval_end,
            }
        )

    fasta.close()
    table_path = out_dir / "candidates.tsv"
    transcript_path = out_dir / "transcripts.tsv"
    candidate_columns = [
        "gene_id",
        "gene_name",
        "gene_biotype",
        "is_placeholder_gene",
        "description",
        "dbxref",
        "ontology_terms",
        "interval",
        "coordinate_basis",
        "gff_seqid",
        "fasta_seqid",
        "chrom",
        "start",
        "end",
        "strand",
        "gene_length",
        "transcript_count",
        "representative_transcript",
        "representative_reason",
        "cds_length",
        "protein_length",
        "cds_phase_status",
        "translation_status",
        "annotation_quality",
        "overlaps_qtl_boundary",
    ]
    transcript_columns = [
        "gene_id",
        "transcript_id",
        "transcript_start",
        "transcript_end",
        "strand",
        "cds_segment_count",
        "cds_length",
        "protein_length",
        "complete_orf",
        "canonical_tag",
        "phase_usable",
        "phase_status",
        "first_cds_phase",
        "reported_phases",
        "expected_phases",
        "phase_warnings",
        "translation_status",
        "is_representative",
        "selection_reason",
    ]
    candidate_frame = pd.DataFrame(table_rows, columns=candidate_columns)
    # A malformed GFF3 can repeat a gene ID (reported in gff_dialect_qc as
    # duplicate_parent_ids). Emitting both rows would let gene_id-keyed evidence
    # double-join, so keep the first occurrence and rely on the QC report to
    # surface the duplication. This makes candidates.tsv honour the gene_id
    # uniqueness contract enforced by schemas.CANDIDATE_SCHEMA.
    candidate_frame = candidate_frame.drop_duplicates(subset="gene_id", keep="first")
    schemas.validate(candidate_frame, schemas.CANDIDATE_SCHEMA)
    candidate_frame.to_csv(
        table_path,
        sep="\t",
        index=False,
    )
    pd.DataFrame(transcript_rows, columns=transcript_columns).to_csv(
        transcript_path,
        sep="\t",
        index=False,
    )
    # The candidate table keeps the first repeated GFF gene ID and reports the
    # malformed parent in QC. Keep FASTA joins aligned with that same choice.
    def first_records(records: list[tuple[str, str]]) -> list[tuple[str, str]]:
        seen: set[str] = set()
        result: list[tuple[str, str]] = []
        for record_id, sequence in records:
            if record_id not in seen:
                result.append((record_id, sequence))
                seen.add(record_id)
        return result

    genomic_records = first_records(genomic_records)
    representative_cds_records = first_records(representative_cds_records)
    transcript_cds_records = first_records(transcript_cds_records)
    representative_protein_records = first_records(representative_protein_records)
    write_fasta(out_dir / "candidate_genomic.fa", genomic_records)
    write_fasta(out_dir / "candidate_cds.fa", representative_cds_records)
    write_fasta(out_dir / "candidate_transcript_cds.fa", transcript_cds_records)
    write_fasta(out_dir / "candidate_proteins.fa", representative_protein_records)

    # Biological protein IDs follow the transcript model, while representative
    # protein FASTA records remain keyed by gene ID for candidate-table joins.
    # Keep both identities explicit instead of silently treating them as equal.
    id_map_rows = [
        {
            "gene_id": str(row["gene_id"]),
            "transcript_id": str(row["transcript_id"]),
            "protein_id": str(row["transcript_id"]) if int(row["protein_length"]) > 0 else "",
            "sequence_record_id": str(row["gene_id"]) if bool(row["is_representative"]) and int(row["protein_length"]) > 0 else "",
            "is_representative": bool(row["is_representative"]),
        }
        for row in transcript_rows
    ]
    id_map_path = out_dir / "id_map.tsv"
    pd.DataFrame(id_map_rows, columns=["gene_id", "transcript_id", "protein_id", "sequence_record_id", "is_representative"]).to_csv(
        id_map_path, sep="\t", index=False
    )

    phase_status_counts = (
        pd.Series(
            [row["phase_status"] for row in transcript_rows],
            dtype=str,
        )
        .value_counts()
        .to_dict()
        if transcript_rows
        else {}
    )
    fasta_fingerprint = file_fingerprint(fasta_path, full_hash)
    gff3_fingerprint = file_fingerprint(gff_path, full_hash)
    if read_only_reference:
        for path, expected in source_stats.items():
            stat = path.stat()
            actual = (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
            if path.is_symlink() or actual != expected:
                raise ValueError(f"shared reference changed during extraction: {path}")
    weak_input_fingerprints = [
        fingerprint["path"]
        for fingerprint in (fasta_fingerprint, gff3_fingerprint)
        if fingerprint_is_weak(fingerprint)
    ]
    manifest = {
        "interval": interval_arg,
        "coordinate_contract": (
            "QTL interval, candidate gene coordinates, CDS models, and downstream IDs "
            "are anchored to this FASTA/GFF3 reference pair."
        ),
        "gff_seqid": chrom,
        "fasta_seqid": fasta_chrom,
        "seqid_mapping_applied": chrom != fasta_chrom,
        "fasta": str(fasta_path),
        "read_only_reference": read_only_reference,
        **({"fasta_index_fingerprint": index_fingerprint,
            "fasta_index_validation": index_validation} if read_only_reference else {}),
        "fasta_fingerprint": fasta_fingerprint,
        "gff3": str(gff_path),
        "gff3_fingerprint": gff3_fingerprint,
        "weak_input_fingerprints": weak_input_fingerprints,
        "reproducibility_note": (
            "Inputs above 512 MB are fingerprinted by metadata (+FAI) unless --full-hash is given; "
            "weak_input_fingerprints lists inputs whose byte content was not fully verified."
        ),
        "candidate_count": len(table_rows),
        "with_cds_count": len(representative_cds_records),
        "with_protein_count": len(representative_protein_records),
        "transcript_count": len(transcript_rows),
        "placeholder_gene_count": sum(1 for row in table_rows if row["is_placeholder_gene"]),
        "cds_phase_status_counts": phase_status_counts,
        "gene_biotype_counts": (
            pd.Series([row["gene_biotype"] for row in table_rows], dtype=str).value_counts().to_dict()
            if table_rows
            else {}
        ),
        "gff_dialect_qc": feature_qc,
        "representative_transcript_rule": (
            "phase-valid > canonical/representative tag > complete ORF > longest annotated CDS"
        ),
        "phase_rule": (
            "Internal CDS phases are validated, not trimmed; only the first "
            "non-zero phase of a 5'-partial CDS is removed before translation."
        ),
        "external_gff_qc": _load_validator_report(gff_qc_report_arg),
        "outputs": {
            "candidates_tsv": str(table_path),
            "transcripts_tsv": str(transcript_path),
            "id_map_tsv": str(id_map_path),
            "candidate_genomic_fasta": str(out_dir / "candidate_genomic.fa"),
            "candidate_cds_fasta": str(out_dir / "candidate_cds.fa"),
            "candidate_transcript_cds_fasta": str(
                out_dir / "candidate_transcript_cds.fa"
            ),
            "candidate_proteins_fasta": str(out_dir / "candidate_proteins.fa"),
        },
    }
    manifest_path = out_dir / "extract_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest | {"manifest": str(manifest_path)}
