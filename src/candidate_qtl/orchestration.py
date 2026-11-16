from __future__ import annotations

import json
import shlex
from pathlib import Path
from typing import Iterable

import pandas as pd
from Bio import SeqIO

from .core import annotation_plan, runtime_path, safe_id


def _sorted_candidates(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    if "dry_group" in frame.columns:
        order = {"top": 0, "reserve": 1, "rescue": 2, "not_selected_yet": 3}
        frame["_dry_group_order"] = frame["dry_group"].str.strip().str.lower().map(order).fillna(9)
    else:
        frame["_dry_group_order"] = 0
    for column in ("dry_rank", "program_rank", "rank"):
        if column in frame.columns:
            frame["_rank_numeric"] = pd.to_numeric(frame[column], errors="coerce").fillna(999999)
            break
    else:
        frame["_rank_numeric"] = range(1, len(frame) + 1)
    return frame.sort_values(["_dry_group_order", "_rank_numeric", "gene_id"], kind="stable")


def _select_candidates(candidate_arg: str, top_n: int) -> tuple[Path, pd.DataFrame]:
    candidate_path = runtime_path(candidate_arg)
    frame = pd.read_csv(candidate_path, sep="\t", dtype=str).fillna("")
    if "gene_id" not in frame.columns:
        raise ValueError("candidate table must contain gene_id")
    if top_n < 1:
        raise ValueError("top_n must be at least 1")
    return candidate_path, _sorted_candidates(frame).head(top_n).copy()


def _db_flags(db_args: Iterable[str]) -> str:
    return " ".join(f"--db {shlex.quote(str(runtime_path(value)))}" for value in db_args)


def _primerblast_prefix(primerblast_root_arg: str) -> str:
    root = runtime_path(primerblast_root_arg)
    return f"PYTHONPATH={shlex.quote(str(root))} python3 -m primerblast_oss"


def _write_script(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    path.chmod(0o755)


def _write_selected_fastas(
    selected: pd.DataFrame,
    fasta_arg: str,
    output_dir: Path,
) -> tuple[list[Path], list[str]]:
    fasta_path = runtime_path(fasta_arg)
    wanted = set(selected["gene_id"].astype(str))
    records_by_id = {}
    with fasta_path.open("r", encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in wanted:
                records_by_id[record.id] = record
    template_dir = output_dir / "templates"
    template_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for gene_id in selected["gene_id"].astype(str):
        record = records_by_id.get(gene_id)
        if record is None:
            continue
        path = template_dir / f"{safe_id(gene_id)}.fa"
        SeqIO.write([record], str(path), "fasta")
        paths.append(path)
    missing = sorted(wanted.difference(records_by_id))
    return paths, missing


def primer_plan(
    stage: str,
    out_dir_arg: str,
    primerblast_root_arg: str,
    genome_arg: str,
    db_args: list[str],
    candidate_arg: str | None = None,
    top_n: int = 10,
    gff_arg: str | None = None,
    candidate_genomic_fasta_arg: str | None = None,
    variants_arg: str | None = None,
    vcf_arg: str | None = None,
    interval: str | None = None,
    n_markers: int = 0,
    spacing: int = 0,
    num_threads: int = 4,
) -> dict[str, object]:
    if stage not in {"direct-sequencing", "variant-marker", "qtl-grid"}:
        raise ValueError("stage must be direct-sequencing, variant-marker, or qtl-grid")
    if not db_args:
        raise ValueError("at least one --db is required")
    if num_threads < 1:
        raise ValueError("num_threads must be at least 1")

    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    genome = runtime_path(genome_arg)
    prefix = _primerblast_prefix(primerblast_root_arg)
    db_flags = _db_flags(db_args)
    common = f"--genome {shlex.quote(str(genome))} {db_flags} --num-threads {num_threads}"
    script_path = out_dir / f"run_{stage.replace('-', '_')}.sh"
    rows: list[dict[str, object]] = []
    lines = ["#!/usr/bin/env bash", "set -euo pipefail", f"mkdir -p {shlex.quote(str(out_dir / 'results'))}"]
    result: dict[str, object] = {"stage": stage, "output_dir": str(out_dir), "script": str(script_path)}

    if stage == "direct-sequencing":
        if not candidate_arg or not gff_arg:
            raise ValueError("direct-sequencing requires --candidates and --gff3")
        candidate_path, selected = _select_candidates(candidate_arg, top_n)
        gff = runtime_path(gff_arg)
        for gene_id in selected["gene_id"].astype(str):
            out_json = out_dir / "results" / f"{safe_id(gene_id)}.assay.json"
            rows.append(
                {
                    "gene_id": gene_id,
                    "stage": stage,
                    "assay_type": "gene_cds_boundary_assay",
                    "purpose": "parental PCR/Sanger before variant-marker conversion",
                    "output": str(out_json),
                }
            )
            lines.append(
                f"{prefix} assay --gene {shlex.quote(gene_id)} --gene-feature cds "
                f"--gff3 {shlex.quote(str(gff))} {common} --flank 250 "
                f"--product-size 350-900 --num-return 10 --format json "
                f"--out {shlex.quote(str(out_json))}"
            )
        missing_templates: list[str] = []
        if candidate_genomic_fasta_arg:
            template_paths, missing_templates = _write_selected_fastas(selected, candidate_genomic_fasta_arg, out_dir)
            for template_path in template_paths:
                gene_id = template_path.stem
                out_json = out_dir / "results" / f"{gene_id}.tiling.json"
                lines.append(
                    f"{prefix} tile --template-fasta {shlex.quote(str(template_path))} "
                    f"{db_flags} --genome-fasta {shlex.quote(str(genome))} "
                    f"--amplicon-min 450 --amplicon-max 850 --overlap 80 "
                    f"--candidates-per-tile 6 --num-threads {num_threads} "
                    f"--format json --out {shlex.quote(str(out_json))}"
                )
        plan_path = out_dir / "direct_sequencing_plan.tsv"
        pd.DataFrame(rows).to_csv(plan_path, sep="\t", index=False)
        result.update(
            {
                "input_candidates": str(candidate_path),
                "selected_count": int(len(selected)),
                "plan": str(plan_path),
                "missing_genomic_templates": missing_templates,
                "note": "Gene assays target CDS spans; optional genomic tiling covers long loci. Validate parent/F1 PCR before sequencing all candidates.",
            }
        )

    elif stage == "variant-marker":
        if not variants_arg:
            raise ValueError("variant-marker requires --variants")
        variant_path = runtime_path(variants_arg)
        variants = pd.read_csv(variant_path, sep="\t", dtype=str).fillna("")
        required = {"variant_id", "chrom", "pos", "ref", "alt", "validation_status"}
        missing = sorted(required.difference(variants.columns))
        if missing:
            raise ValueError(f"variant table is missing required columns: {', '.join(missing)}")
        confirmed_statuses = {"confirmed", "validated"}
        status = variants["validation_status"].str.strip().str.lower()
        excluded_unconfirmed = int((~status.isin(confirmed_statuses)).sum())
        variants = variants[status.isin(confirmed_statuses)].copy()
        if candidate_arg and "gene_id" in variants.columns:
            _, selected = _select_candidates(candidate_arg, top_n)
            variants = variants[variants["gene_id"].isin(set(selected["gene_id"]))].copy()
        vcf_flag = f" --vcf {shlex.quote(str(runtime_path(vcf_arg)))}" if vcf_arg else ""
        for _, row in variants.iterrows():
            variant_id = str(row["variant_id"])
            chrom = str(row["chrom"])
            pos = int(row["pos"])
            ref = str(row["ref"])
            alt = str(row["alt"])
            basename = safe_id(variant_id)
            out_json = out_dir / "results" / f"{basename}.marker.json"
            out_order = out_dir / "results" / f"{basename}.order.tsv"
            if len(ref) == 1 and len(alt) == 1:
                target = f"--snp {shlex.quote(f'{chrom}:{pos}')} --alt {shlex.quote(alt)}"
                marker_class = "SNP/CAPS candidate"
            else:
                start = max(1, pos - 300)
                end = pos + max(len(ref), len(alt)) + 300
                target = f"--interval {shlex.quote(f'{chrom}:{start}-{end}')} --name {shlex.quote(variant_id)}"
                marker_class = "InDel/sequence-length candidate"
            base_command = (
                f"{prefix} assay {target} {common}{vcf_flag} --flank 300 "
                f"--product-size 100-650 --num-return 12"
            )
            lines.append(f"{base_command} --format json --out {shlex.quote(str(out_json))}")
            lines.append(f"{base_command} --format order --out {shlex.quote(str(out_order))}")
            rows.append(
                {
                    "variant_id": variant_id,
                    "gene_id": row.get("gene_id", ""),
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "marker_class": marker_class,
                    "json_output": str(out_json),
                    "order_output": str(out_order),
                }
            )
        plan_path = out_dir / "variant_marker_plan.tsv"
        pd.DataFrame(
            rows,
            columns=["variant_id", "gene_id", "chrom", "pos", "ref", "alt", "marker_class", "json_output", "order_output"],
        ).to_csv(plan_path, sep="\t", index=False)
        result.update(
            {
                "input_variants": str(variant_path),
                "variant_count": int(len(rows)),
                "excluded_unconfirmed_count": excluded_unconfirmed,
                "plan": str(plan_path),
                "note": "Only validation_status=confirmed/validated variants are planned. SNPs are scanned for CAPS-compatible assays; InDels use interval assays.",
            }
        )

    else:
        if not interval:
            raise ValueError("qtl-grid requires --interval")
        if n_markers <= 0 and spacing <= 0:
            n_markers = 8
        out_json = out_dir / "results" / "qtl_grid_markers.json"
        density = f"--n-markers {n_markers}" if n_markers > 0 else f"--spacing {spacing}"
        lines.append(
            f"{prefix} markers --interval {shlex.quote(interval)} {common} "
            f"{density} --marker-flank 300 --product-size 100-500 "
            f"--num-return 8 --format json --out {shlex.quote(str(out_json))}"
        )
        rows.append(
            {
                "stage": stage,
                "interval": interval,
                "n_markers": n_markers,
                "spacing": spacing,
                "output": str(out_json),
                "purpose": "coarse interval reduction or recombinant screening",
            }
        )
        plan_path = out_dir / "qtl_grid_plan.tsv"
        pd.DataFrame(rows).to_csv(plan_path, sep="\t", index=False)
        result.update({"plan": str(plan_path), "interval": interval, "note": "Grid markers are for interval reduction, not causal inference."})

    _write_script(script_path, lines)
    manifest_path = out_dir / f"{stage.replace('-', '_')}_manifest.json"
    manifest_path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    result["manifest"] = str(manifest_path)
    return result


def orchestrate_plan(
    candidates_arg: str,
    out_dir_arg: str,
    protein_fasta_arg: str | None = None,
    blastp_db_args: list[str] | None = None,
    diamond_db_args: list[str] | None = None,
    hmmer_db_arg: str | None = None,
    eggnog_data_dir_arg: str | None = None,
    primerblast_root_arg: str | None = None,
    genome_arg: str | None = None,
    gff_arg: str | None = None,
    nucleotide_db_args: list[str] | None = None,
    candidate_genomic_fasta_arg: str | None = None,
    interval: str | None = None,
    variants_arg: str | None = None,
    vcf_arg: str | None = None,
    top_n: int = 10,
    num_threads: int = 4,
) -> dict[str, object]:
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    steps: list[dict[str, object]] = []
    runnable_scripts: list[Path] = []

    annotation_dir = out_dir / "01_annotation"
    annotation = annotation_plan(
        candidates_arg,
        str(annotation_dir),
        protein_fasta_arg,
        blastp_db_args,
        diamond_db_args,
        hmmer_db_arg,
        eggnog_data_dir_arg,
        "emapper.py",
        num_threads,
    )
    steps.append({"step": 1, "name": "candidate annotation", "tooling": "BLAST/DIAMOND, HMMER/Pfam, eggNOG", "result": annotation})
    if annotation.get("script"):
        runnable_scripts.append(Path(str(annotation["script"])))

    primer_ready = bool(primerblast_root_arg and genome_arg and nucleotide_db_args)
    if primer_ready and gff_arg:
        direct = primer_plan(
            "direct-sequencing",
            str(out_dir / "02_direct_sequencing"),
            str(primerblast_root_arg),
            str(genome_arg),
            list(nucleotide_db_args or []),
            candidate_arg=candidates_arg,
            top_n=top_n,
            gff_arg=gff_arg,
            candidate_genomic_fasta_arg=candidate_genomic_fasta_arg,
            num_threads=num_threads,
        )
        steps.append({"step": 2, "name": "parental direct-sequencing primer design", "tooling": "primerblast-oss assay/tile", "result": direct})
        runnable_scripts.append(Path(str(direct["script"])))
    else:
        steps.append(
            {
                "step": 2,
                "name": "parental direct-sequencing primer design",
                "status": "not generated",
                "requires": ["--primerblast-root", "--genome", "--gff3", "--nucleotide-db"],
            }
        )

    if primer_ready and interval:
        grid = primer_plan(
            "qtl-grid",
            str(out_dir / "00_optional_qtl_grid"),
            str(primerblast_root_arg),
            str(genome_arg),
            list(nucleotide_db_args or []),
            interval=interval,
            n_markers=8,
            num_threads=num_threads,
        )
        steps.append({"step": 0, "name": "optional coarse QTL marker grid", "tooling": "primerblast-oss markers", "result": grid})
        runnable_scripts.insert(0, Path(str(grid["script"])))

    if primer_ready and variants_arg:
        marker = primer_plan(
            "variant-marker",
            str(out_dir / "03_variant_markers"),
            str(primerblast_root_arg),
            str(genome_arg),
            list(nucleotide_db_args or []),
            candidate_arg=candidates_arg,
            top_n=top_n,
            variants_arg=variants_arg,
            vcf_arg=vcf_arg,
            num_threads=num_threads,
        )
        steps.append({"step": 3, "name": "selection-marker design from parental variants", "tooling": "primerblast-oss assay/CAPS/InDel", "result": marker})
        runnable_scripts.append(Path(str(marker["script"])))
    else:
        steps.append(
            {
                "step": 3,
                "name": "selection-marker design from parental variants",
                "status": "deferred until variants are available",
                "requires": ["--variants", "--primerblast-root", "--genome", "--nucleotide-db"],
            }
        )

    run_all = out_dir / "run_available_steps.sh"
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "# Generated orchestration executes only currently available steps.",
        "# Wet gates remain manual: review Dry ranking before sequencing and validate variants before F2 scale-out.",
    ]
    for script in runnable_scripts:
        lines.append(f"bash {shlex.quote(str(script))}")
    _write_script(run_all, lines)

    manifest = {
        "candidates": str(runtime_path(candidates_arg)),
        "output_dir": str(out_dir),
        "laptop_profile": {"threads": num_threads, "memory_target_gb": 16, "execution": "sequential scripts"},
        "orchestration_principle": "Use existing OSS for computation; candidate-qtl coordinates inputs, evidence, gates, and outputs.",
        "steps": sorted(steps, key=lambda item: int(item["step"])),
        "run_available_steps": str(run_all),
        "manual_gates": [
            "Review annotation and optional LLM evidence before choosing Top candidates.",
            "Validate direct-sequencing primers in parents and F1 before scaling.",
            "Only variants confirmed between parents proceed to marker design.",
            "Validate marker calls before F2/F2:3 association.",
        ],
    }
    manifest_path = out_dir / "orchestration_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), "run_available_steps": str(run_all), "step_count": len(steps)}
