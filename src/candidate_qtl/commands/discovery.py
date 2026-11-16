"""CLI commands for the discovery feature group (issue #144)."""

from __future__ import annotations

from ..auto_triage import auto_triage
from ..discover import discover_candidates
from ..homolog_panel import automated_homolog_match, automated_homolog_panel, benchmark_automated_homolog_panels
from argparse import _SubParsersAction
from typing import Any
from ._shared import add_json_flag


def register(subparsers: _SubParsersAction) -> None:
    sub = subparsers.add_parser("discover", help="Run preliminary candidate screening from genome, interval, and phenotype inputs; promote reviewed outputs to case-run for the full workflow.")
    sub.add_argument("--fasta", required=True)
    sub.add_argument("--gff3", required=True)
    sub.add_argument("--interval", action="append", required=True, help="CHROM:START-END; repeat for split blocks.")
    sub.add_argument("--phenotype", required=True)
    sub.add_argument("--phenotype-term", action="append")
    sub.add_argument("--approved-terms", help="Approved trait-term TSV from `trait-terms build`; replaces auto token splitting.")
    sub.add_argument("--direction", choices=["increase", "decrease", "gain", "loss", "unknown"], default="unknown")
    sub.add_argument(
        "--inheritance",
        choices=["recessive", "dominant", "semidominant", "unknown"],
        default="unknown",
        help="Expected segregation model. It affects evidence only when a source-backed row supplies inheritance_model and molecular_effect.",
    )
    sub.add_argument("--evidence", action="append", help="Optional canonical source-backed evidence TSV; repeatable. Rows without source_id are quarantined.")
    sub.add_argument("--protein-db", help="Preferred local protein BLAST database.")
    sub.add_argument("--no-auto-local-db", dest="auto_local_db", action="store_false", help="Disable discovery of the configured local protein DB.")
    sub.set_defaults(auto_local_db=True)
    external = sub.add_mutually_exclusive_group()
    external.add_argument("--offline", dest="allow_external_api", action="store_false", help="Prohibit all external sequence submission; use GFF3/local databases only.")
    external.add_argument("--allow-external-api", dest="allow_external_api", action="store_true", help="Allow candidate protein sequences to be transmitted to the external NCBI BLAST API when no local protein database is selected.")
    sub.set_defaults(allow_external_api=False)
    sub.add_argument("--dry-run", action="store_true", help="Plan only: report whether a real run would transmit candidate sequences externally and write the plan to the manifest, without running BLAST or submitting anything.")
    sub.add_argument("--external-provider", choices=["ncbi-blast"], default="ncbi-blast")
    sub.add_argument("--external-max-candidates", type=int, default=100)
    sub.add_argument("--email")
    sub.add_argument("--top-n", type=int, default=10, help="Dry ranking batch size; also caps Reserve and Rescue slots (default: 10). The Dry score gate is min_top_score=1.0.")
    sub.add_argument("--threads", type=int, default=4)
    sub.add_argument("--flank", type=int, default=500)
    sub.add_argument("--seqid-map", help="GFF_SEQID=FASTA_SEQID pairs or a two-column TSV, for references whose GFF and FASTA seqids differ.")
    sub.add_argument("--gff-qc-report", help="Optional GenomeTools/AGAT validator report to record in the manifest.")
    sub.add_argument("--qtl-source-manifest", help="qtl_source_preflight.json coordinate-provenance gate (same gate as extract).")
    sub.add_argument("--full-hash", action="store_true", help="SHA-256 all inputs, including multi-gigabyte FASTA files.")
    sub.add_argument("--out-dir", required=True)
    add_json_flag(sub)

    sub = subparsers.add_parser("discover-promote", help="Create a reviewable case-run manifest from preliminary discover outputs; writes the discover run's ranking defaults (min_top_score, top_n, reserve_n, rescue_n) into case.yaml.")
    sub.add_argument("--discover-dir", required=True)
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--case-id")
    add_json_flag(sub)

    sub = subparsers.add_parser("homolog-panel-auto", help="Generate a sequence-backed, pending-review homolog panel from phenotype text without inspecting the QTL interval.")
    sub.add_argument("--phenotype", required=True, help="Free-text phenotype supplied before candidate inspection.")
    sub.add_argument("--species", required=True, help="Focal species name, used to exclude exact-species reference entries.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--direction", default="unknown")
    sub.add_argument("--inheritance", default="unknown")
    sub.add_argument("--term", action="append", help="Additional user-supplied phenotype phrase; repeatable.")
    sub.add_argument("--llm-provider", choices=["openai-responses", "openai-compatible", "lmstudio-cli"], help="Optional provider used only to propose search hypotheses.")
    sub.add_argument("--model", help="Model name required with --llm-provider.")
    sub.add_argument("--endpoint", help="Optional LLM HTTP endpoint, or lms executable path for lmstudio-cli.")
    sub.add_argument("--allow-external-api", action="store_true", help="Allow phenotype text to be sent to the UniProt REST API.")
    sub.add_argument("--allow-external-llm", action="store_true", help="Allow phenotype text to be sent to a non-loopback LLM endpoint.")
    sub.add_argument("--panel-size", type=int, default=20)
    sub.add_argument("--per-query", type=int, default=10)
    add_json_flag(sub)

    sub = subparsers.add_parser("homolog-panel-match", help="Match an automated homolog panel only against interval proteins and optionally freeze the deterministic Stage-A 20.")
    sub.add_argument("--candidates", required=True)
    sub.add_argument("--candidate-proteins", required=True)
    sub.add_argument("--panel", required=True)
    sub.add_argument("--panel-fasta", required=True)
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--base-evidence", help="Existing canonical evidence TSV; when supplied, rank and freeze Stage A.")
    sub.add_argument("--target-size", type=int, default=20)
    sub.add_argument("--runner", choices=["auto", "native", "wsl"], default="auto")
    sub.add_argument("--threads", type=int, default=4)
    add_json_flag(sub)

    sub = subparsers.add_parser("homolog-panel-benchmark", help="Blindly compare automated phenotype-to-homolog retrieval on frozen known-gene cases.")
    sub.add_argument("--reference-manifest", required=True, help="benchmark-reference manifest containing frozen phenotype-only cases.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--species", default="Pisum sativum")
    sub.add_argument("--case", action="append", help="Locus ID to benchmark; repeatable. Default: all phenotype-only cases.")
    sub.add_argument("--llm-provider", choices=["openai-responses", "openai-compatible", "lmstudio-cli"])
    sub.add_argument("--model")
    sub.add_argument("--endpoint")
    sub.add_argument("--allow-external-api", action="store_true")
    sub.add_argument("--allow-external-llm", action="store_true")
    sub.add_argument("--panel-size", type=int, default=20)
    sub.add_argument("--per-query", type=int, default=10)
    sub.add_argument("--target-size", type=int, default=20)
    sub.add_argument("--runner", choices=["auto", "native", "wsl"], default="auto")
    sub.add_argument("--threads", type=int, default=4)
    sub.add_argument("--resume", action="store_true")
    sub.add_argument("--benchmark-curated-terms", action="store_true", help="Ablation only: include pre-curated phenotype terms. Primary benchmark uses raw trait text only.")
    add_json_flag(sub)

    sub = subparsers.add_parser("auto-triage", help="Process every interval candidate against role-labelled reference proteomes without a curated reference-gene panel.")
    sub.add_argument("--phenotype", required=True)
    sub.add_argument("--resource-panel", required=True, help="YAML with native and role-labelled local reference proteomes.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--candidates", help="Existing interval candidates TSV; requires --candidate-proteins.")
    sub.add_argument("--candidate-proteins", help="Existing interval candidate protein FASTA; requires --candidates.")
    sub.add_argument("--id-map", help="Optional native gene/transcript/protein/FASTA-record identity table from extract.")
    sub.add_argument("--fasta", help="Native genome FASTA; use with --gff3 and --interval when candidates are not pre-extracted.")
    sub.add_argument("--gff3")
    sub.add_argument("--interval")
    sub.add_argument("--base-evidence", help="Optional canonical evidence TSV to retain alongside automatic screening evidence.")
    sub.add_argument("--parent-vcf", help="Optional parental VCF/VCF.GZ; absence does not block the route.")
    sub.add_argument("--trait-term", action="append", help="Additional input phenotype phrase; repeatable.")
    sub.add_argument("--llm-provider", choices=["openai-responses", "openai-compatible", "lmstudio-cli"], help="Optional phenotype-query expansion only.")
    sub.add_argument("--model")
    sub.add_argument("--endpoint")
    sub.add_argument("--allow-external-llm", action="store_true")
    sub.add_argument("--runner", choices=["auto", "native", "wsl"], default="auto")
    sub.add_argument("--threads", type=int, default=4)
    sub.add_argument("--target-size", type=int, default=20)
    sub.add_argument("--resume", action="store_true")
    add_json_flag(sub)


def _run_discover(args: Any) -> dict:
    result = discover_candidates(
        args.fasta,
        args.gff3,
        args.interval,
        args.phenotype,
        args.out_dir,
        args.direction,
        args.inheritance,
        args.phenotype_term,
        args.protein_db,
        args.allow_external_api,
        args.external_provider,
        args.external_max_candidates,
        args.email,
        args.top_n,
        args.threads,
        args.flank,
        args.auto_local_db,
        args.approved_terms,
        args.seqid_map,
        args.gff_qc_report,
        args.qtl_source_manifest,
        args.full_hash,
        args.evidence,
        args.dry_run,
    )
    return result


def _run_discover_promote(args: Any) -> dict:
    from ..discover import promote_discover

    result = promote_discover(args.discover_dir, args.out_dir, args.case_id)
    return result


def _run_homolog_panel_auto(args: Any) -> dict:
    result = automated_homolog_panel(
        args.phenotype,
        args.species,
        args.out_dir,
        args.direction,
        args.inheritance,
        args.term,
        args.llm_provider,
        args.model,
        args.endpoint,
        args.allow_external_api,
        args.allow_external_llm,
        args.panel_size,
        args.per_query,
    )
    return result


def _run_homolog_panel_match(args: Any) -> dict:
    result = automated_homolog_match(
        args.candidates,
        args.candidate_proteins,
        args.panel,
        args.panel_fasta,
        args.out_dir,
        args.base_evidence,
        args.target_size,
        args.runner,
        args.threads,
    )
    return result


def _run_homolog_panel_benchmark(args: Any) -> dict:
    result = benchmark_automated_homolog_panels(
        args.reference_manifest,
        args.out_dir,
        args.species,
        args.case,
        args.llm_provider,
        args.model,
        args.endpoint,
        args.allow_external_api,
        args.allow_external_llm,
        args.panel_size,
        args.per_query,
        args.target_size,
        args.runner,
        args.threads,
        args.resume,
        args.benchmark_curated_terms,
    )
    return result


def _run_auto_triage(args: Any) -> dict:
    result = auto_triage(
        args.phenotype,
        args.resource_panel,
        args.out_dir,
        args.candidates,
        args.candidate_proteins,
        args.fasta,
        args.gff3,
        args.interval,
        args.base_evidence,
        args.parent_vcf,
        args.trait_term,
        args.llm_provider,
        args.model,
        args.endpoint,
        args.allow_external_llm,
        args.runner,
        args.threads,
        args.target_size,
        args.resume,
        args.id_map,
    )
    return result


HANDLERS = {
    "discover": _run_discover,
    "discover-promote": _run_discover_promote,
    "homolog-panel-auto": _run_homolog_panel_auto,
    "homolog-panel-match": _run_homolog_panel_match,
    "homolog-panel-benchmark": _run_homolog_panel_benchmark,
    "auto-triage": _run_auto_triage,
}
