#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from candidate_qtl.exon_assays import exon_sequencing_plan  # noqa: E402
from candidate_qtl.hardening import llm_evidence_import, research_packet  # noqa: E402


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Candidate-QTL extension commands kept separate from the stable v0.3 CLI.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    sub = subparsers.add_parser("exon-plan", help="Cover every representative-transcript CDS segment with PrimerBlast OSS tiling assays.")
    sub.add_argument("--candidates", required=True)
    sub.add_argument("--gff3", required=True)
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--primerblast-root", required=True)
    sub.add_argument("--genome", required=True)
    sub.add_argument("--db", action="append", required=True)
    sub.add_argument("--seqid-map", help="Optional GFF_SEQID=FASTA_SEQID pairs or a two-column TSV map.")
    sub.add_argument("--top-n", type=int, default=7)
    sub.add_argument("--flank", type=int, default=200)
    sub.add_argument("--max-target-span", type=int, default=500)
    sub.add_argument("--overlap", type=int, default=80)
    sub.add_argument("--num-threads", type=int, default=4)

    sub = subparsers.add_parser("research-packet", help="Create source-backed OSS-to-LLM JSON packets.")
    sub.add_argument("--candidates", required=True)
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--species", required=True)
    sub.add_argument("--trait", required=True)
    sub.add_argument("--top-n", type=int, default=10)
    sub.add_argument("--evidence")
    sub.add_argument("--contributions")
    sub.add_argument("--fetch-literature", action="store_true")
    sub.add_argument("--literature-limit", type=int, default=10)

    sub = subparsers.add_parser("llm-evidence-import", help="Import only human-reviewed LLM claims with stable source IDs.")
    sub.add_argument("--candidates", required=True)
    sub.add_argument("--input", required=True)
    sub.add_argument("--out", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "exon-plan":
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
        elif args.command == "research-packet":
            result = research_packet(
                args.candidates,
                args.out_dir,
                args.species,
                args.trait,
                args.top_n,
                args.fetch_literature,
                args.literature_limit,
                args.evidence,
                args.contributions,
            )
        else:
            result = llm_evidence_import(args.candidates, args.input, args.out)
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 0
    except Exception as error:
        print(json.dumps({"error": str(error), "command": args.command}, ensure_ascii=False), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
