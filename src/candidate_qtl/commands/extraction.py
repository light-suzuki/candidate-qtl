"""CLI commands for the extraction feature group (issue #144)."""

from __future__ import annotations

from ..annotation_rescue import annotation_rescue
from ..core import annotation_plan
from ..gff_preflight import gff_preflight
from ..hardening import annotation_keyword_evidence
from ..workflow_extensions import candidate_flags, dark_candidates, extract_multi, liftover_markers
from argparse import _SubParsersAction
from typing import Any
from ._shared import add_json_flag


def register(subparsers: _SubParsersAction) -> None:
    sub = subparsers.add_parser("extract", help="Extract genes and sequences from one QTL interval.")
    sub.add_argument("--fasta", required=True, help="Reference FASTA; Windows paths are accepted through the PowerShell wrapper.")
    sub.add_argument("--gff3", required=True, help="Reference GFF3 or GFF3.GZ.")
    sub.add_argument("--interval", action="append", required=True, help="CHROM:START-END; repeat for split/non-collinear blocks.")
    sub.add_argument("--out-dir", required=True, help="Output directory.")
    sub.add_argument(
        "--seqid-map",
        help="Optional GFF_SEQID=FASTA_SEQID pairs (comma-separated) or a two-column TSV map.",
    )
    sub.add_argument("--full-hash", action="store_true", help="Read and SHA-256 all inputs, including multi-gigabyte FASTA files.")
    sub.add_argument("--gff-qc-report", help="Optional GenomeTools/AGAT validator report to record in the manifest.")
    sub.add_argument("--flank", type=int, default=0, help="Add gene flanks to candidate_genomic.fa for primer design.")
    sub.add_argument("--qtl-source-manifest", help="qtl_source_preflight.json provenance gate.")
    add_json_flag(sub)

    sub = subparsers.add_parser("gff-preflight", help="Cross-check candidate-qtl GFF interpretation against optional AGAT/gffread outputs.")
    sub.add_argument("--fasta", required=True)
    sub.add_argument("--gff3", required=True)
    sub.add_argument("--interval", action="append", required=True)
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--seqid-map")
    sub.add_argument("--normalized-gff")
    sub.add_argument("--agat-cds")
    sub.add_argument("--gffread-cds")
    sub.add_argument("--execute", action="store_true")
    add_json_flag(sub)

    sub = subparsers.add_parser("annotation-rescue", help="Project known proteins only into a QTL interval and audit missing/broken native models without changing rank.")
    sub.add_argument("--fasta", required=True, help="Target reference genome FASTA.")
    sub.add_argument("--candidates", required=True, help="Native candidates.tsv for the same interval.")
    sub.add_argument("--query-proteins", required=True, help="Source-backed homolog protein FASTA to project.")
    sub.add_argument("--query-metadata", required=True, help="TSV provenance for every query protein; source_id is mandatory.")
    sub.add_argument("--interval", required=True, help="CHROM:START-END; projection is limited to this sequence.")
    sub.add_argument("--out-dir", required=True)
    rescue_source = sub.add_mutually_exclusive_group()
    rescue_source.add_argument("--projection-gff", help="Precomputed miniprot GFF3 against the generated interval-local FASTA.")
    rescue_source.add_argument("--execute", action="store_true", help="Run local miniprot now when available; otherwise emit a runnable script and review rows.")
    sub.add_argument("--projection-complete", action="store_true", help="Assert that an imported GFF is a complete successful run; otherwise absent queries remain not_evaluated.")
    sub.add_argument("--miniprot", default="miniprot", help="miniprot executable name or path.")
    sub.add_argument("--threads", type=int, default=4)
    sub.add_argument("--review-limit", type=int, default=20, help="Maximum review appendix rows; the full ledger is always retained.")
    sub.add_argument("--full-hash", action="store_true")
    add_json_flag(sub)

    sub = subparsers.add_parser("liftover-markers", help="Resolve pair-level marker hits into safe target interval blocks.")
    source = sub.add_mutually_exclusive_group(required=True)
    source.add_argument("--hits", help="Pair table or 10-column BLAST TSV.")
    source.add_argument("--marker-fasta", help="Primer FASTA with IDs ending _F/_R; runs blastn-short.")
    sub.add_argument("--db", help="Nucleotide BLAST DB for --marker-fasta.")
    sub.add_argument("--threads", type=int, default=4)
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--max-product", type=int, default=5000)
    sub.add_argument("--max-block-gap", type=int, default=1000000)
    add_json_flag(sub)

    sub = subparsers.add_parser("dark-candidates", help="Export short complete no-homology ORFs as an explicit dark bucket.")
    sub.add_argument("--candidates", required=True)
    sub.add_argument("--out", required=True)
    sub.add_argument("--max-aa", type=int, default=200)
    add_json_flag(sub)

    sub = subparsers.add_parser("candidate-flags", help="Flag tandem redundancy and dark short-ORF candidates.")
    sub.add_argument("--candidates", required=True)
    sub.add_argument("--out", required=True)
    sub.add_argument("--tandem-distance", type=int, default=100000)
    sub.add_argument("--max-aa", type=int, default=200)
    add_json_flag(sub)

    sub = subparsers.add_parser("annotation-keywords", help="Link best protein hits to GFF descriptions and emit auditable trait-keyword evidence.")
    sub.add_argument("--candidates", required=True)
    sub.add_argument("--blast", required=True, help="Candidate protein BLASTP outfmt 6 with candidate-qtl columns.")
    sources = sub.add_mutually_exclusive_group(required=True)
    sources.add_argument("--annotation-gff3", help="GFF3/GFF3.GZ used to annotate BLAST subjects.")
    sources.add_argument("--subject-annotations", help="Two-column subject_id/description TSV.")
    sources.add_argument("--use-stitle", action="store_true", help="Read stitle from column 14 after the 13-column quality-gated BLAST fields.")
    sub.add_argument("--keyword", action="append", required=True, help="Trait/function keyword or comma-separated keywords; repeatable.")
    sub.add_argument("--out", required=True, help="Canonical evidence TSV output.")
    add_json_flag(sub)

    sub = subparsers.add_parser("annotate-plan", help="Create candidate-only protein, GO evidence template, and local homology/domain commands.")
    sub.add_argument("--candidates", required=True, help="Candidate or ranked TSV containing gene_id.")
    sub.add_argument("--out-dir", required=True, help="Output directory.")
    sub.add_argument("--protein-fasta", help="Protein FASTA whose record IDs match gene_id.")
    sub.add_argument("--blastp-db", action="append", help="BLASTP database; repeatable.")
    sub.add_argument("--diamond-db", action="append", help="DIAMOND database; repeatable.")
    sub.add_argument("--hmmer-db", help="HMMER profile database, for example a local Pfam database.")
    sub.add_argument("--eggnog-data-dir", help="Optional eggNOG-mapper data directory for candidate-only GO/orthology annotation.")
    sub.add_argument("--eggnog-bin", default="emapper.py")
    sub.add_argument("--num-threads", type=int, default=4)
    add_json_flag(sub)


def _run_extract(args: Any) -> dict:
    result = extract_multi(args.fasta, args.gff3, args.interval, args.out_dir, args.seqid_map, args.full_hash, args.gff_qc_report, args.flank, args.qtl_source_manifest)
    return result


def _run_gff_preflight(args: Any) -> dict:
    result = gff_preflight(args.fasta, args.gff3, args.interval, args.out_dir, args.seqid_map, args.normalized_gff, args.agat_cds, args.gffread_cds, args.execute)
    return result


def _run_annotation_rescue(args: Any) -> dict:
    result = annotation_rescue(
        args.fasta,
        args.candidates,
        args.query_proteins,
        args.query_metadata,
        args.interval,
        args.out_dir,
        args.projection_gff,
        args.projection_complete,
        args.execute,
        args.miniprot,
        args.threads,
        args.full_hash,
        args.review_limit,
    )
    return result


def _run_liftover_markers(args: Any) -> dict:
    result = liftover_markers(args.hits, args.out_dir, args.max_product, args.marker_fasta, args.db, args.threads, args.max_block_gap)
    return result


def _run_dark_candidates(args: Any) -> dict:
    result = dark_candidates(args.candidates, args.out, args.max_aa)
    return result


def _run_candidate_flags(args: Any) -> dict:
    result = candidate_flags(args.candidates, args.out, args.tandem_distance, args.max_aa)
    return result


def _run_annotation_keywords(args: Any) -> dict:
    result = annotation_keyword_evidence(
        args.candidates,
        args.blast,
        args.annotation_gff3,
        args.keyword,
        args.out,
        subject_annotations_arg=args.subject_annotations,
        use_stitle=args.use_stitle,
    )
    return result


def _run_annotate_plan(args: Any) -> dict:
    result = annotation_plan(
        args.candidates,
        args.out_dir,
        args.protein_fasta,
        args.blastp_db,
        args.diamond_db,
        args.hmmer_db,
        args.eggnog_data_dir,
        args.eggnog_bin,
        args.num_threads,
    )
    return result


HANDLERS = {
    "extract": _run_extract,
    "gff-preflight": _run_gff_preflight,
    "annotation-rescue": _run_annotation_rescue,
    "liftover-markers": _run_liftover_markers,
    "dark-candidates": _run_dark_candidates,
    "candidate-flags": _run_candidate_flags,
    "annotation-keywords": _run_annotation_keywords,
    "annotate-plan": _run_annotate_plan,
}
