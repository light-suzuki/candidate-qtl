from __future__ import annotations

import gzip
import json
import shlex
from pathlib import Path
from urllib.parse import unquote

import pandas as pd
from Bio.Seq import Seq

from .core import parse_interval, parse_seqid_map, runtime_path, run_tool, safe_id
from .optional_deps import require_pysam


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def _attrs(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[unquote(key)] = unquote(value)
    return result


def _candidate_order(frame: pd.DataFrame) -> pd.DataFrame:
    ordered = frame.copy()
    if "dry_group" in ordered.columns:
        group_order = {"top": 0, "reserve": 1, "rescue": 2, "evidence_insufficient": 3, "not_selected_yet": 3}
        ordered["_group_order"] = ordered["dry_group"].str.casefold().map(group_order).fillna(9)
    else:
        ordered["_group_order"] = 0
    for column in ("dry_rank", "program_rank", "rank"):
        if column in ordered.columns:
            ordered["_rank_order"] = pd.to_numeric(ordered[column], errors="coerce").fillna(999999)
            break
    else:
        ordered["_rank_order"] = range(1, len(ordered) + 1)
    return ordered.sort_values(["_group_order", "_rank_order", "gene_id"], kind="stable")


def _cds_by_parent(gff_path: Path) -> dict[str, list[dict[str, object]]]:
    grouped: dict[str, list[dict[str, object]]] = {}
    with _open_text(gff_path) as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            attrs = _attrs(fields[8])
            parents = [value for value in attrs.get("Parent", "").split(",") if value]
            for parent in parents:
                grouped.setdefault(parent, []).append(
                    {
                        "chrom": fields[0],
                        "start": int(fields[3]),
                        "end": int(fields[4]),
                        "strand": fields[6],
                        "phase": fields[7],
                        "line_number": line_number,
                    }
                )
    for segments in grouped.values():
        segments.sort(key=lambda item: (int(item["start"]), int(item["end"])))
    return grouped


def _tile_segment(start: int, end: int, max_target_span: int, overlap: int) -> list[tuple[int, int]]:
    if start > end:
        raise ValueError("CDS segment start must not exceed end")
    if max_target_span < 50:
        raise ValueError("max_target_span must be at least 50 bp")
    if overlap < 0 or overlap >= max_target_span:
        raise ValueError("overlap must be zero or greater and smaller than max_target_span")
    if end - start + 1 <= max_target_span:
        return [(start, end)]
    step = max_target_span - overlap
    windows: list[tuple[int, int]] = []
    cursor = start
    while cursor <= end:
        window_end = min(end, cursor + max_target_span - 1)
        windows.append((cursor, window_end))
        if window_end == end:
            break
        cursor += step
    return windows


def _write_template(path: Path, record_id: str, sequence: str) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f">{record_id}\n")
        for offset in range(0, len(sequence), 80):
            handle.write(sequence[offset : offset + 80] + "\n")


def exon_sequencing_plan(
    candidates_arg: str,
    gff_arg: str,
    out_dir_arg: str,
    primerblast_root_arg: str,
    genome_arg: str,
    db_args: list[str],
    top_n: int = 7,
    flank: int = 200,
    max_target_span: int = 500,
    overlap: int = 80,
    num_threads: int = 4,
    seqid_map_arg: str | None = None,
) -> dict[str, object]:
    """Extract CDS windows plus flank and design one target-spanning amplicon per window."""
    if top_n < 1:
        raise ValueError("top_n must be at least 1")
    if flank < 0 or num_threads < 1:
        raise ValueError("flank must be non-negative and num_threads must be at least 1")
    if not db_args:
        raise ValueError("at least one nucleotide BLAST --db is required")

    candidates_path = runtime_path(candidates_arg)
    gff_path = runtime_path(gff_arg)
    genome = runtime_path(genome_arg)
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    results_dir = out_dir / "results"
    template_dir = out_dir / "templates"
    results_dir.mkdir(parents=True, exist_ok=True)
    template_dir.mkdir(parents=True, exist_ok=True)

    candidates = pd.read_csv(candidates_path, sep="\t", dtype=str).fillna("")
    required = {"gene_id", "chrom", "strand"}
    missing = sorted(required.difference(candidates.columns))
    if missing:
        raise ValueError(f"candidate table is missing required columns: {', '.join(missing)}")
    selected = _candidate_order(candidates).head(top_n)
    cds_map = _cds_by_parent(gff_path)
    seqid_map = parse_seqid_map(seqid_map_arg)

    pysam = require_pysam('exon-sequencing-plan')
    try:
        fasta = pysam.FastaFile(str(genome))
    except OSError:
        pysam.faidx(str(genome))
        fasta = pysam.FastaFile(str(genome))

    root = runtime_path(primerblast_root_arg)
    prefix = f"PYTHONPATH={shlex.quote(str(root))} python3 -m primerblast_oss"
    db_flags = " ".join(f"--db {shlex.quote(str(runtime_path(value)))}" for value in db_args)
    lines = ["#!/usr/bin/env bash", "set -euo pipefail", f"mkdir -p {shlex.quote(str(results_dir))}"]
    rows: list[dict[str, object]] = []
    missing_cds: list[str] = []

    try:
        for _, candidate in selected.iterrows():
            gene_id = str(candidate["gene_id"])
            transcript_id = str(candidate.get("representative_transcript", "")).strip() or gene_id
            segments = cds_map.get(transcript_id) or cds_map.get(gene_id) or []
            if not segments:
                missing_cds.append(gene_id)
                continue
            for exon_index, segment in enumerate(segments, 1):
                if str(segment["chrom"]) != str(candidate["chrom"]):
                    raise ValueError(
                        f"{gene_id}: representative CDS chromosome {segment['chrom']} does not match candidate chromosome {candidate['chrom']}"
                    )
                fasta_chrom = seqid_map.get(str(segment["chrom"]), str(segment["chrom"]))
                if fasta_chrom not in fasta.references:
                    raise ValueError(
                        f"{gene_id}: GFF sequence ID {segment['chrom']!r} maps to absent FASTA ID {fasta_chrom!r}; use --seqid-map"
                    )
                chrom_length = fasta.get_reference_length(fasta_chrom)
                windows = _tile_segment(int(segment["start"]), int(segment["end"]), max_target_span, overlap)
                for window_index, (start, end) in enumerate(windows, 1):
                    assay_id = safe_id(f"{gene_id}.{transcript_id}.cds{exon_index}.w{window_index}")
                    template_start = max(1, start - flank)
                    template_end = min(chrom_length, end + flank)
                    sequence = fasta.fetch(fasta_chrom, template_start - 1, template_end).upper()
                    if str(segment["strand"]) == "-":
                        sequence = str(Seq(sequence).reverse_complement())
                        local_target_start = template_end - end
                        local_target_end = template_end - start
                    else:
                        local_target_start = start - template_start
                        local_target_end = end - template_start
                    template_path = template_dir / f"{assay_id}.fa"
                    _write_template(template_path, assay_id, sequence)
                    out_json = results_dir / f"{assay_id}.json"
                    target_length = local_target_end - local_target_start + 1
                    product_max = min(len(sequence), 850)
                    product_min = min(product_max, max(100, target_length + 60))
                    lines.append(
                        f"{prefix} design --template-fasta {shlex.quote(str(template_path))} {db_flags} "
                        f"--genome-fasta {shlex.quote(str(genome))} "
                        f"--target {local_target_start},{target_length} --product-size {product_min}-{product_max} "
                        f"--num-return 10 --num-threads {num_threads} "
                        f"--format json --out {shlex.quote(str(out_json))}"
                    )
                    rows.append(
                        {
                            "gene_id": gene_id,
                            "representative_transcript": transcript_id,
                            "strand": candidate["strand"],
                            "cds_segment": exon_index,
                            "window": window_index,
                            "target_interval": f"{segment['chrom']}:{start}-{end}",
                            "target_start": start,
                            "target_end": end,
                            "target_length": end - start + 1,
                            "template_interval": f"{fasta_chrom}:{template_start}-{template_end}",
                            "template_length": len(sequence),
                            "local_target_start": local_target_start,
                            "local_target_end": local_target_end,
                            "template_fasta": str(template_path),
                            "gff_phase": segment["phase"],
                            "output_json": str(out_json),
                            "purpose": "design one amplicon spanning the full CDS window plus splice-boundary flank for parental PCR/Sanger",
                        }
                    )
    finally:
        fasta.close()

    plan_path = out_dir / "exon_sequencing_plan.tsv"
    columns = [
        "gene_id", "representative_transcript", "strand", "cds_segment", "window", "target_interval",
        "target_start", "target_end", "target_length", "template_interval", "template_length", "template_fasta",
        "local_target_start", "local_target_end", "gff_phase", "output_json", "purpose",
    ]
    pd.DataFrame(rows, columns=columns).to_csv(plan_path, sep="\t", index=False)
    script_path = out_dir / "run_exon_sequencing_plan.sh"
    script_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    script_path.chmod(0o755)
    manifest = {
        "input_candidates": str(candidates_path),
        "gff3": str(gff_path),
        "genome": str(genome),
        "selected_gene_count": len(selected),
        "planned_template_count": len(rows),
        "genes_without_cds": missing_cds,
        "max_target_span": max_target_span,
        "target_window_overlap": overlap,
        "flank": flank,
        "plan": str(plan_path),
        "script": str(script_path),
        "interpretation": (
            "Every representative-transcript CDS base is included in a Primer3 target window. PrimerBlast OSS design must return "
            "an amplicon spanning that full target; collect-sanger-markers rejects partial coverage. Parent/F1 PCR validation remains mandatory."
        ),
    }
    manifest_path = out_dir / "exon_sequencing_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest | {"manifest": str(manifest_path)}


def _window_is_covered(payload: dict, target_start: int, target_end: int) -> bool:
    if payload.get("mode") == "tile":
        for tile in payload.get("tiles", []):
            cover = tile.get("covers", [])
            if len(cover) == 2 and int(cover[0]) <= target_start and int(cover[1]) >= target_end:
                return True
        return False
    for pair in payload.get("pairs", []):
        if int(pair.get("left_start", target_start + 1)) <= target_start and int(pair.get("right_start", target_end - 1)) >= target_end:
            return True
    return False


def exon_coverage_check(
    plan_arg: str,
    out_dir_arg: str,
    relax_flank: int = 200,
    relax_product_extra: int = 300,
) -> dict[str, object]:
    """Detect CDS windows the PrimerBlast OSS results do not span and plan retries.

    Reads an exon_sequencing_plan.tsv and the PrimerBlast OSS JSON each window
    points at, then reports covered vs gap windows per gene and writes a relaxed
    retry plan (wider flank/product) for only the failed windows. Missing result
    files are treated as not-yet-run, distinct from a design that failed to span.
    """
    plan_path = runtime_path(plan_arg)
    plan = pd.read_csv(plan_path, sep="\t", dtype=str).fillna("")
    required = {"gene_id", "target_interval", "target_length", "local_target_start", "local_target_end", "output_json"}
    missing = sorted(required.difference(plan.columns))
    if missing:
        raise ValueError(f"exon sequencing plan is missing required columns: {', '.join(missing)}")
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)

    window_rows: list[dict[str, object]] = []
    retry_rows: list[dict[str, object]] = []
    for _, assay in plan.iterrows():
        target_length = int(assay["target_length"])
        result_path = runtime_path(str(assay["output_json"]))
        if not result_path.exists():
            status = "result_missing"
        else:
            payload = json.loads(result_path.read_text(encoding="utf-8"))
            covered = _window_is_covered(payload, int(assay["local_target_start"]), int(assay["local_target_end"]))
            status = "covered" if covered else "coverage_gap"
        window_rows.append(
            {
                "gene_id": assay["gene_id"],
                "target_interval": assay["target_interval"],
                "target_length": target_length,
                "status": status,
                "output_json": str(result_path),
            }
        )
        if status != "covered":
            retry_rows.append(
                {
                    "gene_id": assay["gene_id"],
                    "target_interval": assay["target_interval"],
                    "target_length": target_length,
                    "reason": status,
                    "suggested_flank": relax_flank,
                    "suggested_product_min": target_length + 60,
                    "suggested_product_max": target_length + 60 + relax_product_extra,
                    "note": "Re-run PrimerBlast OSS design with a wider flank/product for this window only.",
                }
            )

    windows = pd.DataFrame(window_rows, columns=["gene_id", "target_interval", "target_length", "status", "output_json"])
    windows_path = out_dir / "coverage_windows.tsv"
    windows.to_csv(windows_path, sep="\t", index=False)

    by_gene = []
    for gene_id, group in windows.groupby("gene_id", sort=True):
        total = len(group)
        covered = int((group["status"] == "covered").sum())
        gaps = int((group["status"] == "coverage_gap").sum())
        missing_results = int((group["status"] == "result_missing").sum())
        covered_bp = int(group.loc[group["status"] == "covered", "target_length"].astype(int).sum())
        total_bp = int(group["target_length"].astype(int).sum())
        by_gene.append(
            {
                "gene_id": gene_id,
                "cds_windows": total,
                "covered_windows": covered,
                "coverage_gap_windows": gaps,
                "result_missing_windows": missing_results,
                "covered_fraction": round(covered / total, 4) if total else 0.0,
                "covered_bp_fraction": round(covered_bp / total_bp, 4) if total_bp else 0.0,
                "fully_covered": covered == total,
            }
        )
    by_gene_frame = pd.DataFrame(by_gene, columns=["gene_id", "cds_windows", "covered_windows", "coverage_gap_windows", "result_missing_windows", "covered_fraction", "covered_bp_fraction", "fully_covered"])
    by_gene_path = out_dir / "coverage_by_gene.tsv"
    by_gene_frame.to_csv(by_gene_path, sep="\t", index=False)

    retry_path = out_dir / "coverage_retry_plan.tsv"
    pd.DataFrame(retry_rows, columns=["gene_id", "target_interval", "target_length", "reason", "suggested_flank", "suggested_product_min", "suggested_product_max", "note"]).to_csv(retry_path, sep="\t", index=False)

    manifest = {
        "plan": str(plan_path),
        "coverage_windows": str(windows_path),
        "coverage_by_gene": str(by_gene_path),
        "retry_plan": str(retry_path),
        "window_count": int(len(windows)),
        "covered_windows": int((windows["status"] == "covered").sum()),
        "coverage_gap_windows": int((windows["status"] == "coverage_gap").sum()),
        "result_missing_windows": int((windows["status"] == "result_missing").sum()),
        "fully_covered_genes": int(by_gene_frame["fully_covered"].sum()) if not by_gene_frame.empty else 0,
        "complete": bool((windows["status"] == "covered").all()) if len(windows) else True,
        "principle": "Every representative-transcript CDS window must be spanned; gaps get a relaxed retry plan, they are not silently dropped.",
    }
    manifest_path = out_dir / "coverage_check_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest | {"manifest": str(manifest_path)}


def collect_sanger_markers(plan_arg: str, out_arg: str, run_script: bool = False) -> dict[str, object]:
    """Collect PrimerBlast OSS tile JSON into a wet-lab-ready Sanger primer table."""
    plan_path = runtime_path(plan_arg)
    output = runtime_path(out_arg)
    plan = pd.read_csv(plan_path, sep="\t", dtype=str).fillna("")
    required = {
        "gene_id", "representative_transcript", "target_interval", "local_target_start", "local_target_end",
        "template_fasta", "output_json",
    }
    missing = sorted(required.difference(plan.columns))
    if missing:
        raise ValueError(f"exon sequencing plan is missing required columns: {', '.join(missing)}")
    if run_script:
        script = plan_path.parent / "run_exon_sequencing_plan.sh"
        if not script.exists():
            raise FileNotFoundError(f"sequencing script not found: {script}")
        completed = run_tool(["bash", str(script)], text=True, capture_output=True, check=False)
        (plan_path.parent / "primerblast_stdout.log").write_text(completed.stdout, encoding="utf-8")
        (plan_path.parent / "primerblast_stderr.log").write_text(completed.stderr, encoding="utf-8")
        if completed.returncode:
            raise RuntimeError(f"PrimerBlast OSS failed (exit {completed.returncode}); see primerblast_stderr.log")

    rows: list[dict[str, object]] = []
    missing_results: list[str] = []
    failed_windows: list[str] = []
    for _, assay in plan.iterrows():
        result_path = runtime_path(str(assay["output_json"]))
        if not result_path.exists():
            missing_results.append(str(result_path))
            continue
        payload = json.loads(result_path.read_text(encoding="utf-8"))
        target_start = int(assay["local_target_start"])
        target_end = int(assay["local_target_end"])
        candidates: list[dict[str, object]] = []
        if payload.get("mode") == "tile":
            for tile in payload.get("tiles", []):
                cover = tile.get("covers", [])
                if len(cover) == 2 and int(cover[0]) <= target_start and int(cover[1]) >= target_end:
                    candidates.append(tile)
        else:
            for pair in payload.get("pairs", []):
                if int(pair.get("left_start", target_start + 1)) <= target_start and int(pair.get("right_start", target_end - 1)) >= target_end:
                    candidates.append(pair)
        if not candidates:
            failed_windows.append(str(assay["target_interval"]))
            continue
        best = candidates[0]
        specificity = best.get("specificity") or {}
        if payload.get("mode") == "tile":
            local_amplicon_start, local_amplicon_end = [int(value) for value in best.get("covers", [0, 0])]
        else:
            local_amplicon_start = int(best.get("left_start", 0))
            local_amplicon_end = int(best.get("right_start", 0))
        genomic_amplicon = ""
        if assay.get("template_interval"):
            chrom, template_start, template_end = parse_interval(str(assay["template_interval"]))
            if str(assay.get("strand", "+")) == "-":
                genomic_start = template_end - local_amplicon_end
                genomic_end = template_end - local_amplicon_start
            else:
                genomic_start = template_start + local_amplicon_start
                genomic_end = template_start + local_amplicon_end
            genomic_amplicon = f"{chrom}:{genomic_start}-{genomic_end}"
        rows.append(
            {
                "gene_id": assay["gene_id"],
                "representative_transcript": assay["representative_transcript"],
                "target_interval": assay["target_interval"],
                "tile_index": best.get("index", ""),
                "forward_primer": best.get("forward", ""),
                "reverse_primer": best.get("reverse", ""),
                "product_size": best.get("product_size", ""),
                "reference_amplicon_interval": genomic_amplicon,
                "tm_forward": best.get("tm_f", ""),
                "tm_reverse": best.get("tm_r", ""),
                "specificity_rank": specificity.get("rank", ""),
                "specific_all_databases": specificity.get("specific_all_db", ""),
                "gel_distinguishable": specificity.get("gel_distinguishable", ""),
                "predicted_off_targets": specificity.get("total_off_target", ""),
                "target_fully_spanned": True,
                "template_fasta": assay["template_fasta"],
                "source_json": str(result_path),
                "wet_status": "design_only_validate_in_parents_and_F1",
            }
        )
    columns = [
        "gene_id", "representative_transcript", "target_interval", "tile_index", "forward_primer", "reverse_primer",
        "product_size", "reference_amplicon_interval", "tm_forward", "tm_reverse", "specificity_rank", "specific_all_databases",
        "gel_distinguishable", "predicted_off_targets", "target_fully_spanned", "template_fasta", "source_json", "wet_status",
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=columns).to_csv(output, sep="\t", index=False)
    manifest = {
        "plan": str(plan_path),
        "output": str(output),
        "marker_count": len(rows),
        "gene_count": len({str(row["gene_id"]) for row in rows}),
        "missing_result_count": len(missing_results),
        "missing_results": missing_results,
        "failed_windows": failed_windows,
        "complete": not missing_results and not failed_windows,
        "warning": "Primer sequences are computational designs; validate amplification, specificity, and parental polymorphism before F2 genotyping.",
    }
    manifest_path = output.with_suffix(".json")
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest | {"manifest": str(manifest_path)}
