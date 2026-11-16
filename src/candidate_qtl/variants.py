from __future__ import annotations

import json
import re
import shutil
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from .core import runtime_path, run_tool, safe_id, sha256


CS_TOKEN = re.compile(r"(:\d+|=[A-Za-z]+|\*[A-Za-z][A-Za-z]|[+-][A-Za-z]+|~[A-Za-z]{2}\d+[A-Za-z]{2})")
MIN_COMPLETE_COVERAGE = 0.90


def parse_cs(cs: str, reference_start: int = 0, query_start: int = 0) -> list[dict[str, object]]:
    """Convert a minimap2 long-form cs tag to a compact, 1-based variant table."""
    reference_position = reference_start
    query_position = query_start
    variants: list[dict[str, object]] = []
    tokens = CS_TOKEN.findall(cs)
    if "".join(tokens) != cs:
        raise ValueError("unsupported or malformed minimap2 cs tag")
    for token in tokens:
        operator = token[0]
        if operator == ":":
            length = int(token[1:])
            reference_position += length
            query_position += length
        elif operator == "=":
            length = len(token) - 1
            reference_position += length
            query_position += length
        elif operator == "*":
            variants.append(
                {
                    "variant_type": "snv",
                    "reference_pos_1based": reference_position + 1,
                    "query_pos_1based": query_position + 1,
                    "reference_allele": token[1].upper(),
                    "alternate_allele": token[2].upper(),
                    "length": 1,
                    "variant_class": "snv",
                }
            )
            reference_position += 1
            query_position += 1
        elif operator == "+":
            sequence = token[1:].upper()
            variants.append(
                {
                    "variant_type": "insertion",
                    "reference_pos_1based": reference_position + 1,
                    "query_pos_1based": query_position + 1,
                    "reference_allele": "-",
                    "alternate_allele": sequence,
                    "length": len(sequence),
                    "variant_class": "inframe_insertion" if len(sequence) % 3 == 0 else "frameshift_insertion",
                }
            )
            query_position += len(sequence)
        elif operator == "-":
            sequence = token[1:].upper()
            variants.append(
                {
                    "variant_type": "deletion",
                    "reference_pos_1based": reference_position + 1,
                    "query_pos_1based": query_position + 1,
                    "reference_allele": sequence,
                    "alternate_allele": "-",
                    "length": len(sequence),
                    "variant_class": "inframe_deletion" if len(sequence) % 3 == 0 else "frameshift_deletion",
                }
            )
            reference_position += len(sequence)
        elif operator == "~":
            match = re.fullmatch(r"~([A-Za-z]{2})(\d+)([A-Za-z]{2})", token)
            if not match:
                raise ValueError(f"malformed cs intron token: {token}")
            length = int(match.group(2))
            variants.append(
                {
                    "variant_type": "long_gap",
                    "reference_pos_1based": reference_position + 1,
                    "query_pos_1based": query_position + 1,
                    "reference_allele": f"<{length}bp>",
                    "alternate_allele": "-",
                    "length": length,
                    "variant_class": "large_indel",
                }
            )
            reference_position += length
    return variants


def alignment_metrics(
    fields: list[str],
    reference_length: int,
    alternate_length: int,
    min_coverage: float = MIN_COMPLETE_COVERAGE,
) -> dict[str, float | int | str | bool]:
    """Summarize one PAF alignment and decide whether absence of variants is interpretable."""
    if len(fields) < 12:
        raise ValueError("PAF alignment requires at least 12 columns")
    if reference_length < 1 or alternate_length < 1:
        raise ValueError("sequence lengths must be positive")
    reference_aligned = max(0, int(fields[8]) - int(fields[7]))
    alternate_aligned = max(0, int(fields[3]) - int(fields[2]))
    block_length = max(1, int(fields[10]))
    reference_coverage = reference_aligned / reference_length
    alternate_coverage = alternate_aligned / alternate_length
    identity = int(fields[9]) / block_length
    return {
        "alignment_strand": fields[4],
        "reference_aligned_bases": reference_aligned,
        "alternate_aligned_bases": alternate_aligned,
        "reference_coverage": reference_coverage,
        "alternate_coverage": alternate_coverage,
        "alignment_identity": identity,
        "mapping_quality": int(fields[11]),
        "comparison_complete": reference_coverage >= min_coverage and alternate_coverage >= min_coverage,
    }


def _interval_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def _alignment_score(fields: list[str]) -> tuple[int, int, int]:
    return (int(fields[11]), int(fields[9]), int(fields[10]))


def classify_alignment_set(
    alignments: list[list[str]],
    reference_length: int,
    alternate_length: int,
    min_coverage: float = MIN_COMPLETE_COVERAGE,
) -> dict[str, object]:
    """Choose a primary alignment and classify multi-alignment structure.

    A single simple full-length alignment yields no structural events. Multiple
    alignments are inspected for inversions (mixed strands), split mappings
    (disjoint query segments), duplications/ambiguity (one query region placed at
    several reference loci with comparable scores), and presence/absence or large
    indels (one sequence only partially aligned). These are returned so callers
    can suppress naive SNV/InDel markers instead of misclassifying the event.
    """
    if not alignments:
        raise ValueError("no alignments to classify")
    ordered = sorted(alignments, key=_alignment_score, reverse=True)
    best = ordered[0]
    events: list[str] = []
    ambiguous = False

    strands = {fields[4] for fields in ordered}
    if len(strands) >= 2:
        events.append("inversion_suspected")

    if len(ordered) >= 2:
        second = ordered[1]
        best_q = (int(best[2]), int(best[3]))
        second_q = (int(second[2]), int(second[3]))
        best_r = (int(best[7]), int(best[8]))
        second_r = (int(second[7]), int(second[8]))
        query_overlap = _interval_overlap(*best_q, *second_q)
        ref_overlap = _interval_overlap(*best_r, *second_r)
        shorter_q = max(1, min(best_q[1] - best_q[0], second_q[1] - second_q[0]))

        # Comparable score on the same query region placed at a different
        # reference locus: ambiguous primary and a duplication signal.
        comparable = _alignment_score(second)[0] == _alignment_score(best)[0] and int(second[9]) >= 0.9 * int(best[9])
        if query_overlap >= 0.5 * shorter_q and comparable:
            ambiguous = True
            if ref_overlap < 0.5 * shorter_q:
                events.append("duplication_suspected")
        elif query_overlap < 0.5 * shorter_q:
            # Different parts of the query align separately.
            events.append("split_alignment")

    metrics = alignment_metrics(best, reference_length, alternate_length, min_coverage)
    if not metrics["comparison_complete"]:
        # One side only partially aligned: presence/absence or a large indel,
        # not a set of small variants.
        events.append("presence_absence_or_large_indel")

    return {
        "primary": best,
        "primary_metrics": metrics,
        "structural_events": sorted(dict.fromkeys(events)),
        "ambiguous_primary_alignment": ambiguous,
        "alignment_count": len(alignments),
    }


def left_normalize_indel(reference_sequence: str, position_1based: int, ref_allele: str, alt_allele: str) -> tuple[int, str, str]:
    """Left-align a pure insertion or deletion against the reference sequence.

    ``position_1based`` is the 1-based reference coordinate of the first affected
    base (deletion) or the base after which sequence is inserted (insertion).
    SNVs and complex alleles are returned unchanged. Coordinates outside the
    supplied reference are returned unchanged so callers never fabricate context.
    """
    is_deletion = alt_allele == "-" and ref_allele != "-"
    is_insertion = ref_allele == "-" and alt_allele != "-"
    if not (is_deletion or is_insertion):
        return position_1based, ref_allele, alt_allele

    seq = reference_sequence.upper()
    pos = position_1based
    if is_deletion:
        allele = ref_allele.upper()
        length = len(allele)
        # Shift left while the base before the deletion equals the last base of it.
        while pos > 1 and (pos - 2) < len(seq) and (pos - 2 + length) <= len(seq) and seq[pos - 2] == seq[pos - 2 + length]:
            pos -= 1
            allele = seq[pos - 1] + allele[:-1]
        return pos, allele, "-"

    allele = alt_allele.upper()
    length = len(allele)
    # Insertion sits after position pos-1 (1-based). Shift left while the base
    # before the insertion equals the last inserted base.
    while pos > 1 and (pos - 2) < len(seq) and seq[pos - 2] == allele[-1]:
        allele = seq[pos - 2] + allele[:-1]
        pos -= 1
    return pos, "-", allele


def _select_record(path: Path, requested_id: str | None):
    records = list(SeqIO.parse(str(path), "fasta"))
    if not records:
        raise ValueError(f"no FASTA records found in {path}")
    if requested_id:
        for record in records:
            aliases = {record.id, *record.id.split("|")}
            if requested_id in aliases:
                return record
        raise ValueError(f"FASTA record {requested_id!r} was not found in {path}")
    if len(records) != 1:
        raise ValueError(f"{path} contains {len(records)} records; select one with an ID option")
    return records[0]


SEQUENCE_TYPES = {"genomic", "cds", "transcript"}


def allele_compare(
    reference_arg: str,
    alternate_arg: str,
    out_dir_arg: str,
    reference_id: str | None = None,
    alternate_id: str | None = None,
    minimap2_bin: str = "minimap2",
    reference_candidates_arg: str | None = None,
    min_coverage: float = MIN_COMPLETE_COVERAGE,
    sequence_type: str | None = None,
) -> dict[str, object]:
    if not 0.0 < min_coverage <= 1.0:
        raise ValueError("min_coverage must be in the interval (0, 1]")
    if sequence_type is not None and sequence_type not in SEQUENCE_TYPES:
        raise ValueError(f"sequence_type must be one of {sorted(SEQUENCE_TYPES)}")
    reference_path = runtime_path(reference_arg)
    alternate_path = runtime_path(alternate_arg)
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    reference = _select_record(reference_path, reference_id)
    alternate = _select_record(alternate_path, alternate_id)
    selected_reference = out_dir / "reference.fa"
    selected_alternate = out_dir / "alternate.fa"
    SeqIO.write([reference], str(selected_reference), "fasta")
    SeqIO.write([alternate], str(selected_alternate), "fasta")
    executable = shutil.which(minimap2_bin) if not Path(minimap2_bin).exists() else minimap2_bin
    if not executable:
        raise ValueError(f"minimap2 executable not found: {minimap2_bin}")
    completed = run_tool(
        [str(executable), "-x", "asm5", "-c", "--cs=long", str(selected_reference), str(selected_alternate)],
        check=True,
        capture_output=True,
        text=True,
    )
    paf_path = out_dir / "alignment.paf"
    paf_path.write_text(completed.stdout, encoding="utf-8")
    alignments = [line.split("\t") for line in completed.stdout.splitlines() if line.strip()]
    if not alignments:
        raise ValueError("minimap2 produced no alignment")
    classification = classify_alignment_set(alignments, len(reference.seq), len(alternate.seq), min_coverage)
    best = classification["primary"]
    metrics = classification["primary_metrics"]
    structural_events = list(classification["structural_events"])
    ambiguous = bool(classification["ambiguous_primary_alignment"])
    tags = {field.split(":", 2)[0]: field.split(":", 2)[2] for field in best[12:] if field.count(":") >= 2}
    if "cs" not in tags:
        raise ValueError("minimap2 alignment is missing the cs tag")
    variants = parse_cs(tags["cs"], int(best[7]), int(best[2]))

    # Modulo-3 labels are meaningful for CDS/transcript comparisons only. For
    # genomic sequence the frame is unknown, so retain neutral genomic classes.
    if sequence_type is not None:
        sequence_context = sequence_type
    elif reference_candidates_arg:
        sequence_context = "genomic"
    else:
        sequence_context = "unspecified_sequence"
    genomic_context = sequence_context in {"genomic", "unspecified_sequence"} and (reference_candidates_arg is not None or sequence_type == "genomic")
    if genomic_context:
        for row in variants:
            if row["variant_type"] == "insertion":
                row["variant_class"] = "genomic_insertion"
            elif row["variant_type"] == "deletion":
                row["variant_class"] = "genomic_deletion"

    for row in variants:
        row["reference_id"] = reference.id
        row["alternate_id"] = alternate.id
        row["alignment_strand"] = best[4]
    columns = [
        "reference_id", "alternate_id", "alignment_strand", "variant_type", "variant_class",
        "reference_pos_1based", "query_pos_1based", "reference_allele", "alternate_allele", "length",
    ]
    variant_path = out_dir / "variants.tsv"
    pd.DataFrame(variants, columns=columns).to_csv(variant_path, sep="\t", index=False)

    marker_path = None
    marker_suppression_reasons: list[str] = []
    if reference_candidates_arg and not bool(metrics["comparison_complete"]):
        marker_suppression_reasons.append("partial_alignment")
    if reference_candidates_arg and best[4] != "+":
        marker_suppression_reasons.append("reverse_alignment_requires_manual_allele_orientation_review")
    if reference_candidates_arg and ambiguous:
        marker_suppression_reasons.append("ambiguous_primary_alignment")
    if reference_candidates_arg:
        marker_suppression_reasons.extend(structural_events)
    if reference_candidates_arg and not marker_suppression_reasons:
        candidates = pd.read_csv(runtime_path(reference_candidates_arg), sep="\t", dtype=str).fillna("")
        required = {"gene_id", "chrom", "start", "end", "strand"}
        if not required.issubset(candidates.columns):
            raise ValueError("--reference-candidates requires gene_id, chrom, start, end, and strand columns")
        matches = candidates[candidates["gene_id"] == reference.id]
        if len(matches) != 1:
            raise ValueError("--reference-candidates must contain exactly one row matching the selected reference ID")
        candidate = matches.iloc[0]
        marker_rows = []
        reference_sequence = str(reference.seq)
        for index, row in enumerate(variants, 1):
            norm_pos, norm_ref, norm_alt = left_normalize_indel(
                reference_sequence,
                int(row["reference_pos_1based"]),
                str(row["reference_allele"]),
                str(row["alternate_allele"]),
            )
            normalized = norm_pos != int(row["reference_pos_1based"])
            offset = norm_pos
            if candidate["strand"] == "-":
                genomic_position = int(candidate["end"]) - offset + 1
                reference_allele = str(Seq(norm_ref).reverse_complement()) if norm_ref != "-" else "-"
                alternate_allele = str(Seq(norm_alt).reverse_complement()) if norm_alt != "-" else "-"
            else:
                genomic_position = int(candidate["start"]) + offset - 1
                reference_allele = norm_ref
                alternate_allele = norm_alt
            marker_rows.append(
                {
                    "variant_id": f"{safe_id(reference.id)}_{index}",
                    "gene_id": reference.id,
                    "chrom": candidate["chrom"],
                    "pos": genomic_position,
                    "ref": reference_allele,
                    "alt": alternate_allele,
                    "variant_type": row["variant_type"],
                    "left_normalized": "yes" if normalized else "no",
                    "validation_status": "unvalidated",
                    "parent_a_allele": reference_allele,
                    "parent_b_allele": alternate_allele,
                    "sanger_trace_or_vcf_source": str(variant_path),
                    "notes": "Assembly comparison only; confirm in the actual parents before marker design.",
                }
            )
        marker_path = out_dir / "marker_candidates.tsv"
        pd.DataFrame(marker_rows).to_csv(marker_path, sep="\t", index=False)

    comparison_status = "complete_alignment" if metrics["comparison_complete"] else "partial_alignment"
    warnings = [
        "Assembly-to-assembly differences are candidate parental variants, not causal proof; verify by parent reads and F2/F3 cosegregation."
    ]
    if not metrics["comparison_complete"]:
        warnings.append(
            "Less than 90% of one or both sequences aligned; variant absence outside the aligned block is unknown and marker export is suppressed."
        )
    if best[4] != "+":
        warnings.append("The best alignment is reverse-strand; marker allele orientation requires manual review and marker export is suppressed.")
    if structural_events:
        warnings.append(
            "Structural signal detected (" + ", ".join(structural_events) + "); these are not simple SNV/InDel events and marker export is suppressed."
        )
    if ambiguous:
        warnings.append("The query region aligns equally well at more than one reference locus; primary placement is ambiguous and marker export is suppressed.")
    summary = {
        "reference_id": reference.id,
        "alternate_id": alternate.id,
        "sequence_context": sequence_context,
        "sequence_type": sequence_type or "inferred",
        "min_coverage": min_coverage,
        "structural_events": structural_events,
        "ambiguous_primary_alignment": ambiguous,
        "reference_length": len(reference.seq),
        "alternate_length": len(alternate.seq),
        "alignment_strand": best[4],
        "reference_aligned_start_1based": int(best[7]) + 1,
        "reference_aligned_end_1based": int(best[8]),
        "alternate_aligned_start_1based": int(best[2]) + 1,
        "alternate_aligned_end_1based": int(best[3]),
        "reference_coverage": round(float(metrics["reference_coverage"]), 6),
        "alternate_coverage": round(float(metrics["alternate_coverage"]), 6),
        "alignment_identity": round(float(metrics["alignment_identity"]), 6),
        "comparison_status": comparison_status,
        "alignment_count": len(alignments),
        "matching_bases": int(best[9]),
        "alignment_block_length": int(best[10]),
        "mapping_quality": int(best[11]),
        "variant_count": len(variants),
        "variant_types": pd.Series([row["variant_class"] for row in variants], dtype=str).value_counts().to_dict(),
        "variants": str(variant_path),
        "alignment_paf": str(paf_path),
        "marker_candidates": str(marker_path) if marker_path else None,
        "marker_suppression_reasons": marker_suppression_reasons,
        "input_sha256": {"reference": sha256(reference_path), "alternate": sha256(alternate_path)},
        "warnings": warnings,
    }
    summary_path = out_dir / "allele_compare.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return summary | {"summary": str(summary_path)}
