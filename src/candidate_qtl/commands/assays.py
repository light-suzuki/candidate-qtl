"""CLI commands for the assays feature group (issue #144)."""

from __future__ import annotations

from ..core import assay_plan
from ..exon_assays import collect_sanger_markers, exon_coverage_check, exon_sequencing_plan
from ..orchestration import primer_plan
from ..primerblast_import import primerblast_import
from argparse import _SubParsersAction
from typing import Any
from ._shared import add_json_flag, add_primerblast_common


def register(subparsers: _SubParsersAction) -> None:
    sub = subparsers.add_parser("exon-coverage-check", help="Detect CDS windows PrimerBlast OSS results do not span and plan relaxed retries.")
    sub.add_argument("--plan", required=True, help="exon_sequencing_plan.tsv from exon-sequencing-plan.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--relax-flank", type=int, default=200, help="Suggested wider flank for retry windows.")
    sub.add_argument("--relax-product-extra", type=int, default=300, help="Extra product-size span suggested for retries.")
    add_json_flag(sub)

    sub = subparsers.add_parser(
        "benchmark-reference",
        help="Run known-gene virtual-QTL recovery cases on one reference without using answer IDs for scoring.",
    )
    sub.add_argument("--config", required=True, help="TSV with case_id, trait, interval, answer_gene_id, trait_terms, and answer_source_id.")
    sub.add_argument("--fasta", required=True)
    sub.add_argument("--gff3", required=True)
    sub.add_argument("--protein-fasta", required=True, help="Official reference protein FASTA.")
    sub.add_argument("--eggnog", required=True, help="Official eggNOG-mapper annotation table for the protein FASTA.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--seqid-map")
    sub.add_argument("--annotation-source-id", default="doi:10.57745/5MDE99#CAMEOR_V2_ANNOTATION_1_prot_eggnog.emapper.annotations")
    sub.add_argument(
        "--interval-size",
        type=int,
        default=0,
        help="Build an answer-centred virtual interval of this size when the config interval is blank; position is not scored.",
    )
    sub.add_argument("--homolog-proteins", action="append", help="FASTA containing curated homolog queries; repeatable.")
    sub.add_argument("--homolog-map", help="TSV keyed by case_id with homolog_query_id and source provenance.")
    sub.add_argument("--query-eligibility-manifest", help="Optional query_eligibility_manifest.json to enforce fail-closed query validation.")
    sub.add_argument("--blastp-runner", choices=("auto", "native", "wsl"), default="auto")
    sub.add_argument("--num-threads", type=int, default=4)
    sub.add_argument("--resume", action="store_true", help="Reuse complete shared interval files and BLAST tables in out-dir.")
    add_json_flag(sub)

    sub = subparsers.add_parser("primerblast-import", help="Normalize PrimerBlast OSS JSON reports into an orderable marker ledger.")
    sub.add_argument("--results", action="append", required=True, help="PrimerBlast OSS JSON report (design or CAPS/InDel/assay); repeatable.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--candidates", help="Optional candidate ledger to join gene IDs by gene_id/transcript/name.")
    sub.add_argument("--export", action="append", choices=["kasp", "asp", "hrm"], help="Emit an allele-specific ordering scaffold; repeatable.")
    add_json_flag(sub)

    sub = subparsers.add_parser("assay-plan", help="Legacy top-candidate primerblast-oss assay plan.")
    sub.add_argument("--candidates", required=True, help="Ranked candidate TSV.")
    sub.add_argument("--out", required=True, help="Output assay-plan TSV.")
    sub.add_argument("--top-n", type=int, default=10)
    sub.add_argument("--genome", help="Reference FASTA for primerblast-oss assay commands.")
    sub.add_argument("--gff3", help="Reference GFF3 for primerblast-oss gene assays.")
    sub.add_argument("--db", action="append", help="BLAST nucleotide DB; repeat for multi-reference screening.")
    sub.add_argument("--primerblast-root", help="primerblast-oss checkout root.")
    add_json_flag(sub)

    sub = subparsers.add_parser(
        "primer-plan",
        help="Generate primerblast-oss plans for direct sequencing, parental-variant markers, or a coarse QTL marker grid.",
    )
    sub.add_argument("--stage", required=True, choices=["direct-sequencing", "variant-marker", "qtl-grid"])
    sub.add_argument("--out-dir", required=True)
    add_primerblast_common(sub)
    sub.add_argument("--candidates", help="Dry-ranked candidate TSV; required for direct-sequencing.")
    sub.add_argument("--top-n", type=int, default=10)
    sub.add_argument("--gff3", help="Reference GFF3; required for direct-sequencing.")
    sub.add_argument("--candidate-genomic-fasta", help="Optional candidate genomic FASTA for tiled Sanger amplicons.")
    sub.add_argument("--variants", help="Confirmed parental variant TSV for variant-marker stage.")
    sub.add_argument("--vcf", help="Optional parental VCF supplied to primerblast-oss.")
    sub.add_argument("--interval", help="QTL interval for qtl-grid, CHROM:START-END.")
    sub.add_argument("--n-markers", type=int, default=0)
    sub.add_argument("--spacing", type=int, default=0)
    add_json_flag(sub)

    sub = subparsers.add_parser("exon-sequencing-plan", help="Tile every Top-candidate representative CDS for parental PCR/Sanger.")
    sub.add_argument("--candidates", required=True)
    sub.add_argument("--gff3", required=True)
    sub.add_argument("--out-dir", required=True)
    add_primerblast_common(sub)
    sub.add_argument("--top-n", type=int, default=7)
    sub.add_argument("--flank", type=int, default=200)
    sub.add_argument("--max-target-span", type=int, default=500)
    sub.add_argument("--overlap", type=int, default=80)
    sub.add_argument("--seqid-map")
    add_json_flag(sub)

    sub = subparsers.add_parser("collect-sanger-markers", help="Collect PrimerBlast OSS tile results into an ordering table.")
    sub.add_argument("--plan", required=True, help="exon_sequencing_plan.tsv")
    sub.add_argument("--out", required=True, help="Output candidate_sanger_markers.tsv")
    sub.add_argument("--run", action="store_true", help="Run the generated PrimerBlast OSS shell script first.")
    add_json_flag(sub)


def _run_exon_coverage_check(args: Any) -> dict:
    result = exon_coverage_check(args.plan, args.out_dir, args.relax_flank, args.relax_product_extra)
    return result


def _run_benchmark_reference(args: Any) -> dict:
    from ..benchmarks import reference_benchmark
    result = reference_benchmark(
        args.config,
        args.fasta,
        args.gff3,
        args.protein_fasta,
        args.eggnog,
        args.out_dir,
        args.seqid_map,
        args.annotation_source_id,
        args.interval_size,
        args.homolog_proteins,
        args.blastp_runner,
        args.num_threads,
        args.homolog_map,
        args.resume,
        query_eligibility_manifest_arg=getattr(args, "query_eligibility_manifest", None),
    )
    return result


def _run_primerblast_import(args: Any) -> dict:
    result = primerblast_import(args.results, args.out_dir, args.candidates, args.export)
    return result


def _run_assay_plan(args: Any) -> dict:
    result = assay_plan(args.candidates, args.out, args.top_n, args.genome, args.gff3, args.db, args.primerblast_root)
    return result


def _run_primer_plan(args: Any) -> dict:
    result = primer_plan(
        args.stage,
        args.out_dir,
        args.primerblast_root,
        args.genome,
        args.db,
        candidate_arg=args.candidates,
        top_n=args.top_n,
        gff_arg=args.gff3,
        candidate_genomic_fasta_arg=args.candidate_genomic_fasta,
        variants_arg=args.variants,
        vcf_arg=args.vcf,
        interval=args.interval,
        n_markers=args.n_markers,
        spacing=args.spacing,
        num_threads=args.num_threads,
    )
    return result


def _run_exon_sequencing_plan(args: Any) -> dict:
    result = exon_sequencing_plan(
        args.candidates,
        args.gff3,
        args.out_dir,
        args.primerblast_root,
        args.genome,
        args.db,
        args.top_n,
        args.flank,
        args.max_target_span,
        args.overlap,
        args.num_threads,
        args.seqid_map,
    )
    return result


def _run_collect_sanger_markers(args: Any) -> dict:
    result = collect_sanger_markers(args.plan, args.out, args.run)
    return result


HANDLERS = {
    "exon-coverage-check": _run_exon_coverage_check,
    "benchmark-reference": _run_benchmark_reference,
    "primerblast-import": _run_primerblast_import,
    "assay-plan": _run_assay_plan,
    "primer-plan": _run_primer_plan,
    "exon-sequencing-plan": _run_exon_sequencing_plan,
    "collect-sanger-markers": _run_collect_sanger_markers,
}
