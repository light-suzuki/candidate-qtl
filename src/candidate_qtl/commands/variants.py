"""CLI commands for the variants feature group (issue #144)."""

from __future__ import annotations

from ..multigenome_private import multigenome_private
from ..promoter_review import promoter_review
from ..structural_variants import sv_plan
from ..variant_evidence import variant_evidence
from ..variants import allele_compare
from argparse import _SubParsersAction
from typing import Any
from ._shared import add_json_flag


def register(subparsers: _SubParsersAction) -> None:
    sub = subparsers.add_parser("variant-evidence", help="Normalize a parental VCF plan and compare SnpEff/VEP transcript consequences without treating predictions as Wet evidence.")
    sub.add_argument("--vcf", required=True)
    sub.add_argument("--fasta", required=True)
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--interval")
    sub.add_argument("--snpeff-vcf")
    sub.add_argument("--vep-vcf")
    sub.add_argument("--execute-norm", action="store_true")
    add_json_flag(sub)

    sub = subparsers.add_parser("multigenome-private", help="Classify focal events against multiple control accessions as private_confirmed / private_possible / shared_confirmed / not_evaluated. Keeps control callability and projection denominators explicit; never adds score and never promotes a candidate.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--focal-vcf", required=True, help="Normalized focal accession VCF/VCF.GZ.")
    sub.add_argument("--gff3", required=True, help="Reference GFF3 used for gene assignment and CDS denominators.")
    sub.add_argument("--case", help="Case manifest containing private_variant_model (focal/controls/callability/loci).")
    sub.add_argument("--control-vcf", action="append", help="ACCESSION=path to a control VCF; repeatable.")
    sub.add_argument("--callable-bed", action="append", help="ACCESSION=path to a callable-regions BED; repeatable.")
    sub.add_argument("--projection-ledger", help="Per reference gene projection state TSV (gene_id/projection_state[/projection_run_complete]).")
    sub.add_argument("--native-genes", help="Optional native annotation gene table, reported as a denominator.")
    sub.add_argument("--reference-consequence", help="Reference-frame consequence TSV (chrom/pos/ref/alt/consequence[/gene_id]).")
    sub.add_argument("--native-consequence", help="Focal native-frame consequence TSV with the same columns.")
    sub.add_argument("--focal-accession", help="Focal accession; overrides the case manifest.")
    sub.add_argument("--phenotype-control", action="append", dest="phenotype_controls", help="Required phenotype control accession; repeatable.")
    sub.add_argument("--comparison-accession", action="append", dest="comparison_accessions", help="Additional comparison accession kept in separate columns; repeatable.")
    sub.add_argument("--case-id", default="")
    sub.add_argument("--full-hash", action="store_true", help="SHA-256 all inputs, including large VCFs.")
    add_json_flag(sub)

    sub = subparsers.add_parser("sv-plan", help="Plan assembly SV/PAV analysis and reconcile SyRI/read-based calls as unscored context.")
    sub.add_argument("--reference", required=True)
    sub.add_argument("--query", required=True)
    sub.add_argument("--qtl-interval", required=True)
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--candidates")
    sub.add_argument("--syri-out")
    sub.add_argument("--call-table", action="append", help="PROVIDER=canonical.tsv; repeatable")
    sub.add_argument("--breakpoint-tolerance", type=int, default=100)
    sub.add_argument("--decision", help="candidate_path_decision.json; must request broaden_to_regulatory_sv_te_cnv")
    add_json_flag(sub)

    sub = subparsers.add_parser("promoter-review", help="Limited non-promoting promoter/5'UTR/3'UTR/SV review route after the coding hypothesis is examined. Emits strand-aware regions and intersects parental variants / SV-PAV-TE as unscored, review-required context; never adds to a gene score.")
    sub.add_argument("--genes", required=True, help="Target genes TSV (gene_id + optional representative_transcript).")
    sub.add_argument("--gff3", required=True, help="Reference GFF3 with gene/mRNA/exon/UTR features.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--parental-variants", help="Parental SNV/indel TSV: chrom/start/end (+ optional flags/feature_id/source_id).")
    sub.add_argument("--sv-calls", help="SV/PAV/TE call TSV: chrom/start/end (+ optional observation_class/flags/feature_id/source_id).")
    sub.add_argument("--coverage", help="Coding exon coverage / Wet-state TSV; gates route start and records unexamined scope.")
    sub.add_argument("--motif-hits", help="Optional reference-only motif hit TSV; recorded as context, never causal.")
    sub.add_argument("--hypotheses", help="Optional locus_hypotheses.tsv; matching typed observations are linked by observation ID.")
    sub.add_argument("--promoter-window", type=int, action="append", help="Promoter window in bp upstream of TSS; repeatable for multiple scopes (default 2000).")
    sub.add_argument("--case-id", default="", help="Locus/case id recorded in the manifest.")
    sub.add_argument("--force", action="store_true", help="Run the review scope even if coding coverage has not been supplied/examined.")
    add_json_flag(sub)

    sub = subparsers.add_parser("allele-compare", help="Call candidate variants from a pair of homologous FASTA sequences with minimap2.")
    sub.add_argument("--reference-fasta", required=True)
    sub.add_argument("--alternate-fasta", required=True)
    sub.add_argument("--reference-id", help="Record ID when reference FASTA contains multiple records.")
    sub.add_argument("--alternate-id", help="Record ID when alternate FASTA contains multiple records.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--minimap2", default="minimap2", help="minimap2 executable name or path.")
    sub.add_argument("--reference-candidates", help="Optional candidates.tsv used to emit genomic marker_candidates.tsv.")
    sub.add_argument("--min-coverage", type=float, default=0.90, help="Minimum reference/alternate coverage for a complete comparison (0-1].")
    sub.add_argument("--sequence-type", choices=["genomic", "cds", "transcript"], help="Declare the compared sequence type; controls frame-aware vs genomic variant classes.")
    add_json_flag(sub)


def _run_variant_evidence(args: Any) -> dict:
    result = variant_evidence(args.vcf, args.fasta, args.out_dir, args.interval, args.snpeff_vcf, args.vep_vcf, args.execute_norm)
    return result


def _run_multigenome_private(args: Any) -> dict:
    return multigenome_private(
        args.out_dir, args.focal_vcf, args.gff3,
        case_arg=args.case, control_vcf_args=args.control_vcf, callable_bed_args=args.callable_bed,
        projection_ledger_arg=args.projection_ledger, native_genes_arg=args.native_genes,
        reference_consequence_arg=args.reference_consequence, native_consequence_arg=args.native_consequence,
        focal_accession=args.focal_accession, phenotype_controls=args.phenotype_controls,
        comparison_accessions=args.comparison_accessions, case_id=args.case_id, full_hash=args.full_hash,
    )


def _run_sv_plan(args: Any) -> dict:
    result = sv_plan(args.reference, args.query, args.qtl_interval, args.out_dir, args.candidates, args.syri_out, args.call_table, args.breakpoint_tolerance, args.decision)
    return result


def _run_promoter_review(args: Any) -> dict:
    return promoter_review(
        args.genes, args.gff3, args.out_dir, args.parental_variants, args.sv_calls,
        args.coverage, args.motif_hits, args.hypotheses, args.promoter_window, args.case_id, args.force,
    )


def _run_allele_compare(args: Any) -> dict:
    result = allele_compare(
        args.reference_fasta,
        args.alternate_fasta,
        args.out_dir,
        args.reference_id,
        args.alternate_id,
        args.minimap2,
        args.reference_candidates,
        args.min_coverage,
        args.sequence_type,
    )
    return result


HANDLERS = {
    "variant-evidence": _run_variant_evidence,
    "multigenome-private": _run_multigenome_private,
    "sv-plan": _run_sv_plan,
    "promoter-review": _run_promoter_review,
    "allele-compare": _run_allele_compare,
}
