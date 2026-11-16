"""QTL-limited protein-to-genome annotation rescue.

The projection is deliberately a review channel, not ranking evidence.  A
protein alignment can reveal a missing or damaged reference model, but it does
not establish that the projected locus causes the trait.
"""

from __future__ import annotations

import json
import hashlib
import shlex
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from . import schemas
from .core import attributes, file_fingerprint, parse_interval, runtime_path, run_tool, safe_id
from .optional_deps import require_pysam


RESCUE_COLUMNS = schemas.ANNOTATION_RESCUE_COLUMNS


def _float_attr(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _int_attr(value: str) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def _query_from_target(target: str) -> tuple[str, int | None, int | None]:
    fields = str(target).split()
    if not fields:
        return "", None, None
    start = _int_attr(fields[1]) if len(fields) > 2 else None
    end = _int_attr(fields[2]) if len(fields) > 2 else None
    return fields[0], start, end


def _parse_miniprot_gff(path: Path) -> list[dict[str, object]]:
    """Parse the documented miniprot mRNA/CDS/stop_codon GFF3 subset."""
    models: dict[str, dict[str, object]] = {}
    cds: dict[str, list[tuple[int, int]]] = defaultdict(list)
    noncanonical_splice: dict[str, int] = defaultdict(int)
    stop_features: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{line_number}: projection GFF requires 9 columns")
            attrs = attributes(fields[8])
            feature_type = fields[2]
            if feature_type == "mRNA":
                projection_id = attrs.get("ID", "").strip()
                query_id, query_start, query_end = _query_from_target(attrs.get("Target", ""))
                if not projection_id or not query_id:
                    raise ValueError(f"{path}:{line_number}: miniprot mRNA requires ID and Target query coordinates")
                if projection_id in models:
                    raise ValueError(f"{path}:{line_number}: duplicate projection ID {projection_id}")
                models[projection_id] = {
                    "query_protein_id": query_id, "projection_id": projection_id,
                    "projection_rank": _int_attr(attrs.get("Rank", "")) or 1,
                    "seqid": fields[0], "start": int(fields[3]), "end": int(fields[4]), "strand": fields[6],
                    "identity": _float_attr(attrs.get("Identity", "")),
                    "positive_fraction": _float_attr(attrs.get("Positive", "")),
                    "frameshift_count": _int_attr(attrs.get("Frameshift", "")),
                    "inframe_stop_count": _int_attr(attrs.get("StopCodon", "")),
                    "query_start": query_start, "query_end": query_end,
                }
            elif feature_type in {"CDS", "stop_codon"}:
                parents = [value for value in attrs.get("Parent", "").split(",") if value]
                for parent in parents:
                    if feature_type == "CDS":
                        cds[parent].append((int(fields[3]), int(fields[4])))
                        if attrs.get("Donor") or attrs.get("Acceptor"):
                            noncanonical_splice[parent] += 1
                    else:
                        stop_features.add(parent)
    for projection_id, model in models.items():
        model["cds"] = cds.get(projection_id, [])
        model["terminal_stop_feature"] = projection_id in stop_features
        model["noncanonical_splice_count"] = noncanonical_splice.get(projection_id, 0)
    dangling = sorted((set(cds) | stop_features) - set(models))
    if dangling:
        raise ValueError(f"{path}: CDS/stop_codon parent has no mRNA: {', '.join(dangling[:5])}")
    return list(models.values())


def _overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def _classify(model: dict[str, object], native: pd.DataFrame) -> tuple[str, list[str], str, float, float]:
    start, end = int(model["start"]), int(model["end"])
    same_chrom = native[
        native["chrom"].astype(str).eq(str(model["chrom"]))
        & native["strand"].astype(str).eq(str(model["strand"]))
    ]
    overlaps: list[tuple[pd.Series, int]] = []
    for _, row in same_chrom.iterrows():
        overlap = _overlap(start, end, int(row["start"]), int(row["end"]))
        if overlap:
            overlaps.append((row, overlap))
    if not overlaps:
        return "reference_rescue", [], "", 0.0, 0.0
    projection_length = end - start + 1
    overlap_total = min(projection_length, sum(value for _, value in overlaps))
    native_ids = sorted({str(row["gene_id"]) for row, _ in overlaps})
    qualities = sorted({str(row.get("annotation_quality", "")) for row, _ in overlaps if str(row.get("annotation_quality", ""))})
    best_native_fraction = max(value / (int(row["end"]) - int(row["start"]) + 1) for row, value in overlaps)
    projection_fraction = overlap_total / projection_length
    bad_quality = bool(set(qualities) & {"no_cds", "invalid_cds_phase", "partial_or_unverified_orf"})
    if len(native_ids) > 1:
        model_class = "ambiguous_overlap"
    elif bad_quality:
        model_class = "reference_rescue"
    elif projection_fraction >= 0.5 and best_native_fraction >= 0.5:
        model_class = "native"
    else:
        model_class = "native_discordant"
    return model_class, native_ids, ";".join(qualities), round(best_native_fraction, 6), round(projection_fraction, 6)


def annotation_rescue(
    fasta_arg: str,
    candidates_arg: str,
    query_proteins_arg: str,
    query_metadata_arg: str,
    interval_arg: str,
    out_arg: str,
    projection_gff_arg: str | None = None,
    projection_complete: bool = False,
    execute: bool = False,
    miniprot_executable: str = "miniprot",
    threads: int = 4,
    full_hash: bool = False,
    review_limit: int = 20,
) -> dict[str, object]:
    """Plan/run a local miniprot projection and emit a non-promoting rescue audit."""
    fasta_path = runtime_path(fasta_arg)
    candidate_path = runtime_path(candidates_arg)
    query_path = runtime_path(query_proteins_arg)
    query_metadata_path = runtime_path(query_metadata_arg)
    out = runtime_path(out_arg)
    out.mkdir(parents=True, exist_ok=True)
    chrom, interval_start, interval_end = parse_interval(interval_arg)
    if threads < 1:
        raise ValueError("threads must be at least 1")
    if review_limit < 0:
        raise ValueError("review_limit must be zero or greater")
    if projection_complete and not projection_gff_arg:
        raise ValueError("projection_complete requires a precomputed projection GFF")
    if projection_gff_arg and execute:
        raise ValueError("choose either a precomputed projection GFF or local execution")

    candidates = pd.read_csv(candidate_path, sep="\t", dtype=str).fillna("")
    schemas.validate(candidates, schemas.CANDIDATE_SCHEMA)
    query_lengths: dict[str, int] = {}
    query_hashes: dict[str, str] = {}
    with query_path.open(encoding="utf-8") as query_handle:
        for record in SeqIO.parse(query_handle, "fasta"):
            if record.id in query_lengths:
                raise ValueError(f"duplicate query protein ID: {record.id}")
            query_lengths[record.id] = len(record.seq)
            query_hashes[record.id] = hashlib.sha256(str(record.seq).upper().encode("ascii")).hexdigest()
    if not query_lengths:
        raise ValueError("query protein FASTA is empty")
    query_metadata = pd.read_csv(query_metadata_path, sep="\t", dtype=str).fillna("")
    metadata_columns = ["query_protein_id", "species", "assembly_id", "annotation_release", "source_id", "query_role", "selection_method", "review_status"]
    missing_metadata = set(metadata_columns) - set(query_metadata.columns)
    if missing_metadata:
        raise ValueError(f"query metadata missing columns: {', '.join(sorted(missing_metadata))}")
    if query_metadata["query_protein_id"].duplicated().any() or set(query_metadata["query_protein_id"]) != set(query_lengths):
        raise ValueError("query metadata must contain each FASTA query ID exactly once")
    for mandatory_column in ("source_id", "query_role", "selection_method"):
        if query_metadata[mandatory_column].str.strip().eq("").any():
            raise ValueError(f"query metadata {mandatory_column} is mandatory")
    if query_metadata["review_status"].str.casefold().isin({"", "rejected"}).any():
        raise ValueError("query metadata review_status must be non-empty and not rejected")
    metadata_by_query = query_metadata.set_index("query_protein_id").to_dict("index")

    pysam = require_pysam('annotation-rescue')
    try:
        fasta = pysam.FastaFile(str(fasta_path))
    except OSError:
        pysam.faidx(str(fasta_path))
        fasta = pysam.FastaFile(str(fasta_path))
    if chrom not in fasta.references:
        fasta.close()
        raise ValueError(f"interval sequence ID {chrom!r} is absent from FASTA")
    interval_end = min(interval_end, fasta.get_reference_length(chrom))
    interval_sequence = fasta.fetch(chrom, interval_start - 1, interval_end).upper()
    fasta.close()
    interval_record_id = f"{safe_id(chrom)}__{interval_start}_{interval_end}"
    interval_fasta = out / "qtl_interval.fa"
    interval_fasta.write_text(f">{interval_record_id} source={chrom}:{interval_start}-{interval_end}\n{interval_sequence}\n", encoding="utf-8", newline="\n")

    generated_gff = out / "miniprot_projection.gff3"
    command = [miniprot_executable, f"-t{threads}", "--gff", "--gff-only", "-I", str(interval_fasta), str(query_path)]
    script = out / "run_local_annotation_rescue.sh"
    script_tmp = out / "miniprot_projection.gff3.tmp"
    script.write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\n"
        + " ".join(shlex.quote(value) for value in command)
        + f" > {shlex.quote(str(script_tmp))}\n"
        + f"mv -- {shlex.quote(str(script_tmp))} {shlex.quote(str(generated_gff))}\n",
        encoding="utf-8", newline="\n",
    )

    tool_status = "not_run"
    tool_error = ""
    projection_path = runtime_path(projection_gff_arg) if projection_gff_arg else generated_gff
    if projection_gff_arg:
        if not projection_path.exists():
            raise ValueError(f"projection GFF not found: {projection_path}")
        tool_status = "imported_complete" if projection_complete else "imported_partial"
    tool_version = ""
    if execute:
        executable = shutil.which(miniprot_executable)
        if not executable:
            tool_status = "unavailable"
            tool_error = f"{miniprot_executable} not found; run the generated local script where miniprot is installed"
        else:
            command[0] = executable
            version_result = run_tool([executable, "--version"], capture_output=True, text=True, check=False)
            tool_version = (version_result.stdout or version_result.stderr).strip().splitlines()[0] if (version_result.stdout or version_result.stderr).strip() else "unknown"
            temporary_gff = out / "miniprot_projection.gff3.tmp"
            with temporary_gff.open("w", encoding="utf-8", newline="\n") as handle:
                completed = run_tool(command, stdout=handle, stderr=subprocess.PIPE, text=True, check=False)
            if completed.returncode:
                tool_status = "failed"
                tool_error = completed.stderr.strip()[:2000]
                temporary_gff.unlink(missing_ok=True)
            else:
                try:
                    _parse_miniprot_gff(temporary_gff)
                except ValueError as error:
                    tool_status = "failed"
                    tool_error = f"invalid miniprot GFF: {error}"
                    temporary_gff.unlink(missing_ok=True)
                else:
                    temporary_gff.replace(generated_gff)
                    tool_status = "complete"

    projection_path = generated_gff if tool_status == "complete" else projection_path
    models = _parse_miniprot_gff(projection_path) if tool_status in {"complete", "imported_complete", "imported_partial"} else []
    rows: list[dict[str, object]] = []
    mapped_queries: set[str] = set()
    hit_counts = pd.Series([str(model["query_protein_id"]) for model in models]).value_counts().to_dict() if models else {}
    for model in models:
        query_id = str(model["query_protein_id"])
        if query_id not in query_lengths:
            raise ValueError(f"projection query {query_id!r} is absent from query protein FASTA")
        mapped_queries.add(query_id)
        # The generated interval FASTA resets coordinates to 1.  Imported GFFs
        # must use that same interval-local coordinate basis.
        local_start, local_end = int(model["start"]), int(model["end"])
        if str(model["seqid"]) != interval_record_id or local_start < 1 or local_end > len(interval_sequence):
            raise ValueError(
                "projection GFF must use coordinates from the generated qtl_interval.fa "
                f"record {interval_record_id!r}"
            )
        model["chrom"] = chrom
        model["start"] = interval_start + local_start - 1
        model["end"] = interval_start + local_end - 1
        qstart, qend = model.get("query_start"), model.get("query_end")
        if qstart is not None and qend is not None and not (1 <= int(qstart) <= int(qend) <= query_lengths[query_id]):
            raise ValueError(f"projection Target coordinates are outside query {query_id!r}")
        query_coverage = ((int(qend) - int(qstart) + 1) / query_lengths[query_id]) if qstart is not None and qend is not None else None
        model_class, native_ids, native_quality, native_fraction, projection_fraction = _classify(model, candidates)
        cds_parts = list(model["cds"])
        cds_length = sum(end - start + 1 for start, end in cds_parts)
        frameshifts, internal_stops = int(model["frameshift_count"]), int(model["inframe_stop_count"])
        if frameshifts or internal_stops:
            orf_audit = "disrupted_projection"
        elif query_coverage is not None and query_coverage >= 0.95:
            orf_audit = "near_complete_query_projection_no_reported_disruption"
        elif query_coverage is not None:
            orf_audit = "partial_query_projection"
        else:
            orf_audit = "coverage_unknown"
        splice_audit = "multi_exon_projection" if len(cds_parts) > 1 else "single_exon_projection" if cds_parts else "cds_missing"
        boundary_truncation = local_start == 1 or local_end == len(interval_sequence)
        identity = model["identity"]
        positive = model["positive_fraction"]
        if identity is not None and not 0 <= identity <= 1:
            raise ValueError(f"projection Identity must be a fraction from 0 to 1 for {query_id!r}")
        if positive is not None and not 0 <= positive <= 1:
            raise ValueError(f"projection Positive must be a fraction from 0 to 1 for {query_id!r}")
        metadata = metadata_by_query[query_id]
        strict_projection = (
            query_coverage is not None and query_coverage >= 0.8 and identity is not None and identity >= 0.5
            and not frameshifts and not internal_stops and len(cds_parts) > 0 and not boundary_truncation
            and int(model["noncanonical_splice_count"]) == 0 and hit_counts[query_id] == 1 and int(model["projection_rank"]) == 1
            and str(metadata["review_status"]).casefold() in {"reviewed", "approved"}
        )
        quality_tier = "high" if strict_projection else "medium" if query_coverage is not None and query_coverage >= 0.5 and identity is not None and identity >= 0.3 else "low"
        review_reasons = []
        if model_class != "native": review_reasons.append(model_class)
        if orf_audit != "near_complete_query_projection_no_reported_disruption": review_reasons.append(orf_audit)
        if query_coverage is not None and query_coverage < 0.5: review_reasons.append("low_query_coverage")
        if boundary_truncation: review_reasons.append("interval_boundary_truncation")
        if int(model["noncanonical_splice_count"]): review_reasons.append("noncanonical_splice")
        rows.append({
            "query_protein_id": query_id, "query_sequence_sha256":query_hashes[query_id], "query_source_id":metadata["source_id"], "query_species":metadata["species"], "query_assembly_id":metadata["assembly_id"],
            "query_annotation_release":metadata["annotation_release"], "query_role":metadata["query_role"], "query_selection_method":metadata["selection_method"], "query_review_status":metadata["review_status"],
            "projection_id": model["projection_id"], "projection_rank": model["projection_rank"], "hit_count":hit_counts[query_id],
            "near_tie":hit_counts[query_id] > 1 and int(model["projection_rank"]) <= 2,
            "model_class": model_class, "quality_tier":quality_tier, "provider_status":tool_status,
            "chrom": chrom, "start": model["start"], "end": model["end"], "strand": model["strand"],
            "native_gene_ids": ";".join(native_ids), "native_annotation_quality": native_quality,
            "query_coverage": "" if query_coverage is None else round(query_coverage, 6),
            "identity": "" if identity is None else identity,
            "positive_fraction": "" if positive is None else positive,
            "projected_cds_count": len(cds_parts), "projected_cds_length": cds_length,
            "frameshift_count": frameshifts, "inframe_stop_count": internal_stops, "orf_audit": orf_audit,
            "splice_audit": splice_audit, "noncanonical_splice_count":model["noncanonical_splice_count"], "terminal_stop_feature":model["terminal_stop_feature"],
            "interval_boundary_truncation":boundary_truncation, "native_overlap_fraction": native_fraction, "projection_overlap_fraction": projection_fraction,
            "review_status": "pending" if review_reasons else "not_required", "promotion_status": "not_promoted",
            "source_type": "local_protein_to_genome_projection", "source_id": metadata["source_id"], "notes": ";".join(review_reasons),
        })
    for query_id in sorted(set(query_lengths) - mapped_queries):
        metadata = metadata_by_query[query_id]
        evaluated = tool_status in {"complete", "imported_complete"}
        rows.append({
            "query_protein_id": query_id, "query_sequence_sha256":query_hashes[query_id], "query_source_id":metadata["source_id"], "query_species":metadata["species"], "query_assembly_id":metadata["assembly_id"],
            "query_annotation_release":metadata["annotation_release"], "query_role":metadata["query_role"], "query_selection_method":metadata["selection_method"], "query_review_status":metadata["review_status"],
            "projection_id": "", "projection_rank": "", "hit_count":0, "near_tie":False, "model_class": "unmapped" if evaluated else "not_evaluated",
            "quality_tier":"not_assessed", "provider_status":tool_status,
            "chrom": chrom, "start": "", "end": "", "strand": "", "native_gene_ids": "", "native_annotation_quality": "",
            "query_coverage": 0, "identity": "", "positive_fraction": "", "projected_cds_count": 0, "projected_cds_length": 0,
            "frameshift_count": 0, "inframe_stop_count": 0, "orf_audit": "unmapped", "splice_audit": "unmapped",
            "noncanonical_splice_count":0, "terminal_stop_feature":False, "interval_boundary_truncation":False,
            "native_overlap_fraction": 0, "projection_overlap_fraction": 0, "review_status": "pending", "promotion_status": "not_promoted",
            "source_type": "local_protein_to_genome_projection", "source_id": metadata["source_id"],
            "notes": "no_projection_in_qtl_interval" if evaluated else f"provider_{tool_status}",
        })
    frame = pd.DataFrame(rows, columns=RESCUE_COLUMNS).sort_values(["model_class", "query_protein_id", "projection_rank"], kind="stable")
    schemas.validate(frame.fillna(""), schemas.ANNOTATION_RESCUE_SCHEMA)
    rescue_path = out / "annotation_rescue.tsv"
    frame.to_csv(rescue_path, sep="\t", index=False, lineterminator="\n")
    pending_review = frame[frame["review_status"].eq("pending")].copy()
    pending_review["_quality_order"] = pending_review["quality_tier"].map({"high":0, "medium":1, "low":2, "not_assessed":3}).fillna(4)
    pending_review = pending_review.sort_values(
        ["_quality_order", "query_protein_id", "projection_rank"], ascending=[True, True, True], kind="stable"
    )
    review = pending_review.head(review_limit)[
        ["query_protein_id", "projection_id", "model_class", "native_gene_ids", "orf_audit", "splice_audit", "notes"]
    ].copy()
    review["required_next_action"] = "inspect projection/native model; validate with RNA evidence, parental sequence, or assembly/SV evidence"
    review["allowed_actions"] = "retain rescue | repair annotation | branch SV/PAV | reject projection; never auto-promote"
    review_path = out / "annotation_rescue_review.tsv"
    review.to_csv(review_path, sep="\t", index=False, lineterminator="\n")
    overflow_path = out / "annotation_rescue_review_overflow.tsv"
    pending_review.iloc[review_limit:].drop(columns=["_quality_order"]).to_csv(overflow_path, sep="\t", index=False, lineterminator="\n")

    manifest = {
        "schema_version": "1.0", "interval": f"{chrom}:{interval_start}-{interval_end}", "coordinate_basis": "1-based inclusive",
        "projection_coordinate_basis": "interval_local_converted_to_reference", "projection_provider": "miniprot",
        "projection_command": command, "tool_status": tool_status, "tool_version":tool_version, "tool_error": tool_error,
        "external_sequence_submission_this_run": False, "external_sequence_destinations": [], "local_only": True,
        "causal_promotion_allowed": False, "ranking_effect": "none", "model_count": len(models),
        "threshold_profile": {"high_query_coverage_min":0.8, "high_identity_min":0.5, "medium_query_coverage_min":0.5, "medium_identity_min":0.3,
                              "high_requires_unique_rank1":True, "high_requires_no_boundary_or_reported_disruption":True,
                              "thresholds_are_screening_policy_not_biological_truth":True},
        "class_counts": {str(key): int(value) for key, value in frame["model_class"].value_counts().sort_index().items()},
        "review_count": len(review), "review_overflow_count":max(0, len(pending_review) - review_limit), "review_limit":review_limit,
        "inputs": {
            "fasta": file_fingerprint(fasta_path, full_hash), "candidates": file_fingerprint(candidate_path, full_hash),
            "query_proteins": file_fingerprint(query_path, full_hash),
            "query_metadata": file_fingerprint(query_metadata_path, full_hash),
            "projection_gff": file_fingerprint(projection_path, full_hash) if projection_path.exists() else {"path": str(projection_path), "status": "not_available"},
        },
        "scientific_limits": [
            "protein projection is annotation evidence, not trait-causality evidence",
            "near-complete projection means query coverage without reported disruption, not experimentally proven ORF function",
            "unmapped does not prove gene absence; divergence, sequence gaps, thresholds, and query choice remain alternatives",
        ],
    }
    manifest_path = out / "annotation_rescue_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"rescue_table": str(rescue_path), "review_table": str(review_path), "review_overflow":str(overflow_path), "manifest": str(manifest_path), "script": str(script),
            "class_counts": manifest["class_counts"], "review_count": len(review), "tool_status": tool_status}
