"""CLI commands for the mapping feature group (issue #144)."""

from __future__ import annotations

from ..core import cosegregation
from ..qtl_mapping import qtl2_plan
from ..scaffold_rescue import scaffold_rescue
from argparse import _SubParsersAction
from typing import Any
from ._shared import add_json_flag


def register(subparsers: _SubParsersAction) -> None:
    sub = subparsers.add_parser("scaffold-rescue", help="Conservatively rescue unplaced / unanchored scaffolds that are plausibly linked to the QTL locus, as a separate candidate class. Being unplaced is not itself a reason; family agreement alone is never high confidence.")
    sub.add_argument("--scaffolds", required=True, help="Scaffold table (scaffold_id[/length]).")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--marker-hits", help="Marker/flanking alignment (marker_id/scaffold_id[/pair_id/identity]); paired markers beat lone hits and multi-mapping markers are rejected.")
    sub.add_argument("--synteny", help="Synteny bridge (scaffold_id[/block_id/anchor_species]).")
    sub.add_argument("--orthogroup-context", help="Gene-content context (scaffold_id[/orthogroup]); weak support only.")
    sub.add_argument("--placement-support", help="Hi-C / linkage map / optical map / assembler placement (scaffold_id/support_type[/source]).")
    sub.add_argument("--scaffold-genes", help="Genes on scaffolds (gene_id/scaffold_id).")
    sub.add_argument("--max-marker-scaffolds", type=int, default=1, help="A marker hitting more scaffolds than this is treated as a repeat.")
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)

    sub = subparsers.add_parser("cosegregation", help="Summarize phenotype distribution by a candidate marker genotype.")
    sub.add_argument("--input", required=True, help="TSV based on templates/f2_genotype_phenotype.tsv.")
    sub.add_argument("--out", required=True, help="Output summary TSV.")
    sub.add_argument("--genotype-column", default="candidate_marker_genotype")
    sub.add_argument("--phenotype-column", default="phenotype_raw")
    add_json_flag(sub)

    sub = subparsers.add_parser("qtl2-plan", help="Validate F2 tables and generate or run an R/qtl2 interval scan.")
    sub.add_argument("--genotypes", required=True, help="Wide TSV: individual_id then AA/AB/BB marker columns.")
    sub.add_argument("--phenotypes", required=True, help="TSV with individual_id, trait, and optional covariates.")
    sub.add_argument("--markers", required=True, help="TSV with marker_id, chrom, pos_bp, and optional pos_cm.")
    sub.add_argument("--trait", required=True)
    sub.add_argument("--covariate", action="append", help="Phenotype covariate; repeatable.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--n-perm", type=int, default=1000)
    sub.add_argument("--alpha", type=float, default=0.05)
    sub.add_argument("--error-prob", type=float, default=0.002)
    sub.add_argument("--execute", action="store_true", help="Run the generated script; requires R/qtl2.")
    sub.add_argument("--rscript", default="Rscript")
    add_json_flag(sub)


def _run_cosegregation(args: Any) -> dict:
    result = cosegregation(args.input, args.out, args.genotype_column, args.phenotype_column)
    return result


def _run_qtl2_plan(args: Any) -> dict:
    result = qtl2_plan(
        args.genotypes,
        args.phenotypes,
        args.markers,
        args.out_dir,
        args.trait,
        args.covariate,
        args.n_perm,
        args.alpha,
        args.error_prob,
        args.execute,
        args.rscript,
    )
    return result


def _run_scaffold_rescue(args: Any) -> dict:
    return scaffold_rescue(
        args.scaffolds, args.out_dir, args.marker_hits, args.synteny, args.orthogroup_context,
        args.placement_support, args.scaffold_genes, args.max_marker_scaffolds, args.case_id,
    )


HANDLERS = {
    "scaffold-rescue": _run_scaffold_rescue,
    "cosegregation": _run_cosegregation,
    "qtl2-plan": _run_qtl2_plan,
}
