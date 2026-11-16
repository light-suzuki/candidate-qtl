from __future__ import annotations

import json
import shutil

import pandas as pd

from .core import runtime_path, run_tool


VALID_GENOTYPES = {"AA", "AB", "BA", "BB", "-", "NA", "N", ""}


def qtl2_plan(
    genotypes_arg: str,
    phenotypes_arg: str,
    markers_arg: str,
    out_dir_arg: str,
    trait: str,
    covariates: list[str] | None = None,
    n_perm: int = 1000,
    alpha: float = 0.05,
    error_prob: float = 0.002,
    execute: bool = False,
    rscript: str = "Rscript",
) -> dict[str, object]:
    """Validate F2 tables, make an R/qtl2 bundle, and optionally execute it."""
    if n_perm < 0:
        raise ValueError("n_perm must be zero or greater")
    if not 0 < alpha < 1:
        raise ValueError("alpha must be between 0 and 1")
    if not 0 <= error_prob < 0.5:
        raise ValueError("error_prob must be between 0 and 0.5")

    genotypes_path = runtime_path(genotypes_arg)
    phenotypes_path = runtime_path(phenotypes_arg)
    markers_path = runtime_path(markers_arg)
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)

    geno = pd.read_csv(genotypes_path, sep="\t", dtype=str).fillna("")
    pheno = pd.read_csv(phenotypes_path, sep="\t", dtype=str).fillna("")
    markers = pd.read_csv(markers_path, sep="\t", dtype=str).fillna("")
    if geno.columns.empty or geno.columns[0] != "individual_id":
        raise ValueError("genotypes must be wide TSV with individual_id as the first column")
    if "individual_id" not in pheno.columns or trait not in pheno.columns:
        raise ValueError(f"phenotypes require individual_id and trait column {trait!r}")
    required_marker_columns = {"marker_id", "chrom", "pos_bp"}
    missing = sorted(required_marker_columns.difference(markers.columns))
    if missing:
        raise ValueError(f"markers table is missing required columns: {', '.join(missing)}")
    marker_ids = markers["marker_id"].astype(str).tolist()
    if len(marker_ids) != len(set(marker_ids)):
        raise ValueError("marker_id values must be unique")
    if set(geno.columns[1:]) != set(marker_ids):
        absent = sorted(set(marker_ids).difference(geno.columns[1:]))
        extra = sorted(set(geno.columns[1:]).difference(marker_ids))
        raise ValueError(f"genotype marker columns do not match markers table; absent={absent}, extra={extra}")
    if geno["individual_id"].duplicated().any() or pheno["individual_id"].duplicated().any():
        raise ValueError("individual_id values must be unique in both tables")
    common_ids = [value for value in geno["individual_id"] if value in set(pheno["individual_id"])]
    if len(common_ids) < 20:
        raise ValueError("at least 20 F2 individuals with both genotype and phenotype are required")
    geno = geno.set_index("individual_id").loc[common_ids, marker_ids].reset_index()
    pheno = pheno.set_index("individual_id").loc[common_ids].reset_index()
    covariates = covariates or []
    unknown_covariates = [value for value in covariates if value not in pheno.columns]
    if unknown_covariates:
        raise ValueError(f"unknown phenotype covariates: {', '.join(unknown_covariates)}")

    invalid: set[str] = set()
    for marker_id in marker_ids:
        values = {str(value).strip().upper() for value in geno[marker_id]}
        invalid.update(values.difference(VALID_GENOTYPES))
        geno[marker_id] = geno[marker_id].astype(str).str.strip().str.upper().replace({"BA": "AB", "": "-", "NA": "-", "N": "-"})
    if invalid:
        raise ValueError(f"unsupported genotype codes: {', '.join(sorted(invalid))}; use AA/AB/BB or missing '-'")
    numeric_trait = pd.to_numeric(pheno[trait], errors="coerce")
    if numeric_trait.notna().sum() < 20:
        raise ValueError("R/qtl2 scan requires at least 20 numeric phenotypes for the selected trait")
    pheno[trait] = numeric_trait
    for covariate in covariates:
        numeric_covar = pd.to_numeric(pheno[covariate], errors="coerce")
        if numeric_covar.isna().any():
            raise ValueError(f"covariate {covariate!r} must be numeric and complete; encode factors as numeric columns")
        pheno[covariate] = numeric_covar

    markers["pos_bp"] = pd.to_numeric(markers["pos_bp"], errors="raise").astype(int)
    if "pos_cm" in markers.columns and pd.to_numeric(markers["pos_cm"], errors="coerce").notna().all():
        markers["pos_cm"] = pd.to_numeric(markers["pos_cm"], errors="raise")
        map_source = "provided_pos_cm"
    else:
        markers = markers.sort_values(["chrom", "pos_bp", "marker_id"], kind="stable")
        markers["pos_cm"] = markers.groupby("chrom", sort=False).cumcount().astype(float)
        map_source = "provisional_1cM_marker_order"

    geno_csv = out_dir / "geno.csv"
    pheno_csv = out_dir / "pheno.csv"
    covar_csv = out_dir / "covar.csv"
    gmap_csv = out_dir / "gmap.csv"
    pmap_csv = out_dir / "pmap.csv"
    control = out_dir / "qtl2_control.yaml"
    script = out_dir / "run_qtl2.R"
    geno.to_csv(geno_csv, index=False)
    pheno[["individual_id", trait]].to_csv(pheno_csv, index=False)
    if covariates:
        pheno[["individual_id"] + covariates].to_csv(covar_csv, index=False)
    markers[["marker_id", "chrom", "pos_cm"]].to_csv(gmap_csv, index=False, header=["marker", "chr", "pos"])
    physical_map = markers[["marker_id", "chrom", "pos_bp"]].copy()
    physical_map["pos_bp"] = physical_map["pos_bp"] / 1_000_000.0
    physical_map.to_csv(pmap_csv, index=False, header=["marker", "chr", "pos"])
    covar_control = "covar: covar.csv\n" if covariates else ""
    control.write_text(
        "crosstype: f2\n"
        "geno: geno.csv\npheno: pheno.csv\n" + covar_control + "gmap: gmap.csv\npmap: pmap.csv\n"
        "alleles: [A, B]\n"
        "genotypes:\n  AA: 1\n  AB: 2\n  BB: 3\n"
        "na.strings: ['-', 'NA', '']\n",
        encoding="utf-8",
    )

    covar_vector = "NULL" if not covariates else f"as.matrix(cross$covar[, c({', '.join(json.dumps(x) for x in covariates)}), drop=FALSE])"
    provisional_r = "TRUE" if map_source.startswith("provisional") else "FALSE"
    script.write_text(
        "args <- commandArgs(trailingOnly=TRUE)\n"
        "outdir <- if (length(args)) normalizePath(args[[1]], mustWork=TRUE) else getwd()\n"
        "suppressPackageStartupMessages(library(qtl2))\n"
        "cross <- read_cross2(file.path(outdir, 'qtl2_control.yaml'))\n"
        f"trait <- {json.dumps(trait)}\n"
        f"addcovar <- {covar_vector}\n"
        "if (!is.null(addcovar)) storage.mode(addcovar) <- 'numeric'\n"
        f"if ({provisional_r}) cross$gmap <- est_map(cross, error_prob={error_prob})\n"
        f"pr <- calc_genoprob(cross, error_prob={error_prob})\n"
        "scan <- scan1(pr, cross$pheno[, trait, drop=FALSE], addcovar=addcovar)\n"
        "scan_df <- cbind(marker=rownames(scan), as.data.frame(scan), stringsAsFactors=FALSE)\n"
        "pmap_df <- data.frame(marker=names(unlist(cross$pmap)), chr=rep(names(cross$pmap), lengths(cross$pmap)), pos_Mbp=as.numeric(unlist(cross$pmap)))\n"
        "scan_df <- merge(scan_df, pmap_df, by='marker', all.x=TRUE)\n"
        "write.table(scan_df, file.path(outdir, 'lod_profile.tsv'), sep='\\t', quote=FALSE, row.names=FALSE)\n"
        f"n_perm <- {n_perm}\nalpha <- {alpha}\n"
        "threshold <- NA_real_\n"
        "if (n_perm > 0) {\n"
        "  set.seed(20260711)\n"
        "  perm <- scan1perm(pr, cross$pheno[, trait, drop=FALSE], addcovar=addcovar, n_perm=n_perm)\n"
        "  threshold <- as.numeric(summary(perm, alpha=alpha)[1])\n"
        "  write.table(data.frame(alpha=alpha, n_perm=n_perm, lod_threshold=threshold), file.path(outdir, 'permutation_threshold.tsv'), sep='\\t', quote=FALSE, row.names=FALSE)\n"
        "}\n"
        "peaks <- find_peaks(scan, cross$pmap, threshold=ifelse(is.na(threshold), 0, threshold), drop=1.5)\n"
        "write.table(peaks, file.path(outdir, 'qtl_peaks_mbp.tsv'), sep='\\t', quote=FALSE, row.names=FALSE)\n"
        "peaks_bp <- peaks\n"
        "for (nm in intersect(c('pos','ci_lo','ci_hi'), names(peaks_bp))) peaks_bp[[nm]] <- round(peaks_bp[[nm]] * 1e6)\n"
        "names(peaks_bp)[names(peaks_bp) == 'pos'] <- 'pos_bp'\n"
        "names(peaks_bp)[names(peaks_bp) == 'ci_lo'] <- 'ci_lo_bp'\n"
        "names(peaks_bp)[names(peaks_bp) == 'ci_hi'] <- 'ci_hi_bp'\n"
        "write.table(peaks_bp, file.path(outdir, 'qtl_peaks.tsv'), sep='\\t', quote=FALSE, row.names=FALSE)\n"
        "intervals <- data.frame(chr=character(), start=integer(), end=integer(), interval=character())\n"
        "if (all(c('chr','ci_lo_bp','ci_hi_bp') %in% names(peaks_bp))) {\n"
        "  intervals <- data.frame(chr=peaks_bp$chr, start=peaks_bp$ci_lo_bp, end=peaks_bp$ci_hi_bp, interval=paste0(peaks_bp$chr, ':', peaks_bp$ci_lo_bp, '-', peaks_bp$ci_hi_bp))\n"
        "}\n"
        "write.table(intervals, file.path(outdir, 'candidate_qtl_intervals.tsv'), sep='\\t', quote=FALSE, row.names=FALSE)\n"
        "saveRDS(calc_errorlod(cross, pr), file.path(outdir, 'genotyping_error_lod.rds'))\n"
        f"vit <- viterbi(cross, error_prob={error_prob})\n"
        "xo <- locate_xo(vit, cross$gmap)\n"
        "saveRDS(xo, file.path(outdir, 'crossovers.rds'))\n"
        "xo_count <- as.data.frame(count_xo(vit))\n"
        "xo_count <- cbind(individual_id=rownames(xo_count), xo_count)\n"
        "write.table(xo_count, file.path(outdir, 'crossover_counts.tsv'), sep='\\t', quote=FALSE, row.names=FALSE)\n"
        "writeLines(capture.output(print(summary(cross))), file.path(outdir, 'qtl2_summary.txt'))\n",
        encoding="utf-8",
    )

    run = None
    if execute:
        executable = shutil.which(rscript)
        if not executable:
            raise RuntimeError(f"Rscript not found: {rscript}")
        run = run_tool([executable, str(script), str(out_dir)], text=True, capture_output=True, check=False)
        (out_dir / "qtl2_stdout.log").write_text(run.stdout, encoding="utf-8")
        (out_dir / "qtl2_stderr.log").write_text(run.stderr, encoding="utf-8")
        if run.returncode:
            detail = run.stderr.strip()[-1500:] if run.stderr.strip() else "no stderr"
            raise RuntimeError(f"R/qtl2 failed (exit {run.returncode}): {detail}; see qtl2_stderr.log")

    manifest = {
        "cross_type": "f2",
        "individual_count": len(common_ids),
        "marker_count": len(marker_ids),
        "trait": trait,
        "covariates": covariates,
        "map_source": map_source,
        "map_estimated_on_execute": execute and map_source.startswith("provisional"),
        "control": str(control),
        "script": str(script),
        "executed": execute,
        "outputs": {
            "lod_profile": str(out_dir / "lod_profile.tsv"),
            "permutation_threshold": str(out_dir / "permutation_threshold.tsv"),
            "qtl_peaks": str(out_dir / "qtl_peaks.tsv"),
            "qtl_peaks_mbp": str(out_dir / "qtl_peaks_mbp.tsv"),
            "candidate_intervals": str(out_dir / "candidate_qtl_intervals.tsv"),
            "crossovers": str(out_dir / "crossovers.rds"),
            "crossover_counts": str(out_dir / "crossover_counts.tsv"),
            "genotyping_error_lod": str(out_dir / "genotyping_error_lod.rds"),
        },
        "warning": (
            "The provisional 1 cM marker-order map is re-estimated from F2 recombination when --execute is used; "
            "the generated-but-not-run bundle is not an inferred genetic map."
            if map_source.startswith("provisional") else ""
        ),
    }
    manifest_path = out_dir / "qtl2_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest | {"manifest": str(manifest_path)}
