"""CLI commands for Route C: embedding retrieval discovery (issue #193).

Route C is a third, independent discovery route: phenotype text is embedded and
retrieved against a frozen reference embedding database; the QTL interval is
matched only after the retrieval is frozen. The import command deliberately has
no interval argument, so interval information cannot enter the retrieval input.
"""

from __future__ import annotations

import argparse
from typing import Any

from ..corpus_masking import TIERS as MASK_TIERS
from ..corpus_masking import corpus_mask
from ..reference_corpus import reference_corpus_build
from ..reference_corpus_compare import reference_corpus_compare
from ..reference_sources import uniprot_fetch
from ..retrieval_local import retrieval_index_build, retrieval_run
from ..retrieval_route import retrieval_import, retrieval_interval_match, retrieval_merge, triple_route_compare
from ..retrieval_semantic import (
    PROVIDER_LEXICAL,
    PROVIDER_LSA,
    PROVIDER_NEURAL,
    PROVIDERS,
    retrieval_hybrid_merge,
    retrieval_neural_index_build,
    retrieval_neural_run,
    retrieval_semantic_index_build,
    retrieval_semantic_run,
)
from ..route_c_benchmark import route_c_benchmark, route_c_provider_comparison
from ..route_c_direct_pipeline import (
    direct_literature_search,
    direct_rescue_benchmark,
    direct_rescue_freeze,
    direct_rescue_project,
)
from ..route_c_direct_members import build_route_c_direct_members
from ..phenotype_variant_benchmark import phenotype_variant_benchmark
from ._shared import add_json_flag

# Provider dispatch for the shared --provider flag on retrieval-index-build /
# retrieval-run (issue #225 #3). "lexical" (issue #217's deterministic TF-IDF
# baseline) stays the default so existing behaviour is unchanged; "lsa" is
# the latent-semantic provider (issue #225 #1); "neural" is the
# sentence_transformers adapter, correctly gated to
# not_evaluated_provider_unavailable on a workstation that lacks the backend
# (issue #225 #2). None of the three functions behind these ids accept an
# interval/atlas/candidate argument.
_INDEX_BUILD_PROVIDERS = {
    PROVIDER_LEXICAL: retrieval_index_build,
    PROVIDER_LSA: retrieval_semantic_index_build,
    PROVIDER_NEURAL: retrieval_neural_index_build,
}
_RUN_PROVIDERS = {
    PROVIDER_LEXICAL: retrieval_run,
    PROVIDER_LSA: retrieval_semantic_run,
    PROVIDER_NEURAL: retrieval_neural_run,
}


def register(subparsers: argparse._SubParsersAction) -> None:
    sub = subparsers.add_parser("reference-corpus-fetch", help="Download source-backed reference records (UniProtKB: function, disruption phenotype, pathway, tissue, GO, evidence codes and the PubMed ids each block cites) for one or more taxa. Requires --allow-external-api; there is no free-text query argument, so a case phenotype cannot shape the corpus.")
    sub.add_argument("--taxon", action="append", help="NCBI taxon id (3702 Arabidopsis, 3880 Medicago, 3847 soybean, 4530 rice, 4577 maize, 4081 tomato ...). Repeatable; may be combined with --preset (union; an explicit --taxon wins over the same id from the preset).")
    sub.add_argument("--preset", help="Versioned corpus preset from the packaged refcorpus_presets.json (e.g. legume_v1): expands to its taxa and per-taxon reviewed flags, so a taxon that needs TrEMBL entries is fetched with them. May be combined with --taxon.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--allow-external-api", action="store_true", help="Required. Only the taxon id and the field list are transmitted.")
    sub.add_argument("--include-unreviewed", action="store_true", help="Include unreviewed (TrEMBL) entries; reviewed-only by default because curated comment blocks are what the routes need.")
    sub.add_argument("--max-rows", type=int, default=0, help="Stop after this many rows per taxon (0 = all).")
    sub.add_argument("--page-size", type=int, default=500)
    add_json_flag(sub)

    sub = subparsers.add_parser("reference-corpus-build", help="Build the versioned, phenotype-blind reference corpus that Routes B and C both consume. Retrieved rows come from local authoritative sources; a model-proposed row is promoted only when its accession appears in a retrieved source and a citation resolves. Unverified proposals stay in an exclusion ledger.")
    sub.add_argument("--from-uniprot", action="append", dest="uniprot", help="Downloaded UniProtKB table from reference-corpus-fetch (or an equivalent export): curated function, disruption phenotype, pathway, tissue/stage, GO and evidence codes. Repeatable.")
    sub.add_argument("--from-blastdb", action="append", dest="blastdb", help="Local protein BLAST database prefix; every entry is read with blastdbcmd and its UniProt-style defline parsed. Repeatable.")
    sub.add_argument("--from-dump", action="append", dest="dump", help="Exported reference table (TSV with reference_gene_id ...). Repeatable.")
    sub.add_argument("--from-payload", action="append", dest="payload", help="Model/agent-proposed rows (.tsv/.json/.jsonl) adding function and mutant-phenotype text. Repeatable.")
    sub.add_argument("--orthogroup-map", action="append", dest="orthogroup_map", help="TSV mapping reference_gene_id to orthogroup_id/family_id from a source-backed orthology resource. Repeatable.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--corpus-version", required=True, help="Version label; the corpus is frozen under it and referenced by that name downstream.")
    sub.add_argument("--citation-cache", help="Cached Europe PMC records used to resolve the citations of proposed rows offline.")
    sub.add_argument("--benchmark-mode", action="store_true", help="Bar focal-lineage / answer-shaped rows and record that the corpus was built for a benchmark.")
    sub.add_argument("--preset-acceptance", help="Preset id from the packaged refcorpus_presets.json: the build verifies the preset's required gene symbols are present in the corpus as usable gene symbols (casefold; fail-closed: one missing required term aborts before anything is written, so a previously frozen corpus cannot be silently replaced) and records the per-term result in the manifest. Optional terms are recorded only.")
    sub.add_argument("--full-hash", action="store_true")
    add_json_flag(sub)

    sub = subparsers.add_parser("reference-corpus-compare", help="Compare two or more frozen corpora: per-field and per-species coverage (mutant phenotype, pathway, GO, evidence) against a baseline label, plus optional Route C Recall@K recomputed from frozen retrievals per corpus. Ranking weights are untouched and scores are never compared across corpora.")
    sub.add_argument("--corpus", action="append", required=True, help="LABEL=PATH for a reference corpus TSV. Repeat at least twice.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--baseline-label", default="", help="Label the deltas are measured against; defaults to the first --corpus.")
    sub.add_argument("--retrieval-manifest", help="Optional TSV (corpus_label/case_id/retrieval_tsv/truth_family_id[/truth_gene_id]) of frozen Route C retrievals produced from each corpus for the same cases.")
    sub.add_argument("--full-hash", action="store_true")
    add_json_flag(sub)

    sub = subparsers.add_parser("corpus-mask", help="Issue #237: remove the truth gene / its family / its close orthogroup from a frozen reference corpus BEFORE prediction, so a benchmark can measure recovery of genuinely unfamiliar biology. Mask rules come only from a benchmark-construction mask spec and are fingerprinted in the manifest; no prediction command ever sees them, only the masked corpus TSV. Refuses to write a masked corpus that still names the masked gene/family (self-check).")
    sub.add_argument("--corpus", required=True, help="Frozen reference corpus TSV from reference-corpus-build (or a prior corpus-mask output, for stacking tiers).")
    sub.add_argument("--mask-spec", required=True, help="TSV: case_id/truth_gene_id (+ optional truth_family_id, truth_orthogroup_id, aliases), generated from truth in the benchmark-construction phase. May hold every case; --case-id selects one row.")
    sub.add_argument("--tier", required=True, choices=list(MASK_TIERS), help="leave_gene_out removes the exact reference gene; leave_family_out also removes its family/subfamily and clearly-equivalent aliases (requires truth_family_id); leave_close_orthogroup_out additionally removes the whole orthogroup (requires truth_orthogroup_id too). Each tier is a superset of the one before it.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--case-id", default="", help="Select one row of --mask-spec; optional when the spec has exactly one row.")
    sub.add_argument("--full-hash", action="store_true")
    add_json_flag(sub)

    sub = subparsers.add_parser("retrieval-index-build", help="Route C: build an immutable, versioned local reference index from source-backed plant records. Rows without a family/orthogroup or without a source id are excluded before Top-K and kept in a ledger; every artifact is SHA-256 fingerprinted.")
    sub.add_argument("--corpus", action="append", required=True, help="Reference corpus TSV (reference_gene_id/reference_species/source_id + family_id and/or orthogroup_id + description/function_description/mutant_phenotype/pathway/go_terms/keywords). Repeatable.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--database-version", required=True, help="Version label for this index; retrieval scores are only comparable within one database_version.")
    sub.add_argument("--corpus-label", default="", help="Human-readable corpus name recorded in the manifest.")
    sub.add_argument("--provider", choices=list(PROVIDERS), default=PROVIDER_LEXICAL, help="Retrieval provider (issue #225): 'lexical' (default, deterministic TF-IDF, issue #217), 'lsa' (latent-semantic: TF-IDF + truncated SVD, no new dependency), or 'neural' (sentence_transformers adapter; reports not_evaluated_provider_unavailable when the backend is not installed, as it is not on this workstation).")
    sub.add_argument("--full-hash", action="store_true")
    add_json_flag(sub)

    sub = subparsers.add_parser("retrieval-run", help="Route C: run a phenotype-only local retrieval against a frozen index and emit the retrieval contract. The command has no interval/atlas/candidate argument, so Phase 1 cannot see the QTL interval.")
    sub.add_argument("--index", required=True, help="index_manifest.json (or the index directory) from retrieval-index-build. Its artifacts are re-hashed before retrieval.")
    sub.add_argument("--phenotype", required=True, help="Phenotype text. This is the only query input; it never leaves the machine.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--top-k", type=int, default=50)
    sub.add_argument("--case-id", default="")
    sub.add_argument("--phenotype-context", help="Optional phenotype_context.yaml (case_id, developmental_stage, tissue, environment_or_treatment, ontology term ...). Recorded for audit in the retrieval status; the retrieval query stays the phenotype free text and no structured field is appended to it. Supported for the default lexical provider.")
    sub.add_argument("--provider", choices=list(PROVIDERS), default=PROVIDER_LEXICAL, help="Must match the provider the --index was built with ('lexical' default / 'lsa' / 'neural').")
    add_json_flag(sub)

    sub = subparsers.add_parser("retrieval-merge", help="Route C: merge two or more per-database retrieval TSVs by within-database rank. Raw scores are never compared across databases; unannotated rows are excluded before Top-K selection.")
    sub.add_argument("--retrieval", action="append", required=True, help="Per-database retrieval TSV. Repeat at least twice; each input must carry one unique database_version.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--top-k", type=int, default=50)
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)

    sub = subparsers.add_parser("retrieval-hybrid-merge", help="Route C (issue #225): fuse two or more per-provider retrieval TSVs (e.g. lexical + lsa) by reciprocal-rank fusion over each provider's within-provider rank. Raw retrieval_score values are never compared or summed across providers -- only rank is read.")
    sub.add_argument("--retrieval", action="append", required=True, help="Per-provider retrieval TSV. Repeat at least twice; each input must carry exactly one model_name/model_version pair (one provider).")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--top-k", type=int, default=50)
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)

    sub = subparsers.add_parser("retrieval-import", help="Route C: import an external phenotype/protein retrieval TSV, validate it against the documented contract, and freeze it. No interval input exists on this command: the retrieval stays phenotype-only and interval-blind until the separate interval-match command.")
    sub.add_argument("--retrieval", required=True, help="External retrieval TSV following templates/route_c_retrieval.tsv (reference_gene_id/reference_species/family_id|orthogroup_id/retrieval_score/matched_modality/matched_description_or_phenotype/model_name/model_version/database_version/rank). Columns outside the contract are rejected.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--phenotype", default="", help="Phenotype text recorded in the status; not sent to any model by this command.")
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)

    sub = subparsers.add_parser("retrieval-interval-match", help="Route C: match the frozen retrieval artifact against the QTL interval (interval atlas). Runs only after retrieval-import; the frozen retrieval is never modified here.")
    sub.add_argument("--retrieval", required=True, help="route_c_retrieval_frozen.tsv produced by retrieval-import.")
    sub.add_argument("--atlas", required=True, help="interval_gene_atlas.tsv (gene_id/family_id/orthogroup/orthology_reference_gene): the canonical <out-dir>/interval_gene_atlas.tsv from route-b-prepare or the root-level interval-atlas file.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)

    sub = subparsers.add_parser("triple-route-compare", help="Compare all three independently generated routes (A: phenotype-first, B: interval-first, C: embedding retrieval). Convergence is recorded as state and adds no score.")
    sub.add_argument("--route-a", required=True, help="phenotype_first_hypotheses.tsv (route A).")
    sub.add_argument("--route-b", required=True, help="family_process_hypotheses.tsv (route B).")
    sub.add_argument("--route-c", required=True, help="route_c_retrieval_frozen.tsv (route C).")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--case-id", default="")
    sub.add_argument("--route-c-interval-match", help="Optional route_c_interval_match.tsv (Phase 2). Without it a Route C-only family cannot claim interval anchoring, which stays the honest state rather than an assumed one.")
    add_json_flag(sub)

    sub = subparsers.add_parser("route-c-benchmark", help="Score frozen Route C retrieval artifacts on known-gene cases: Family Recall@5/@10/@20, Gene Recall@20, and A-only / C-only / A+C union. Interval-blind by construction.")
    sub.add_argument("--manifest", required=True, help="TSV: one row per case (case_id/retrieval_tsv/truth_family_id[/truth_gene_id/route_a_tsv]).")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--full-hash", action="store_true")
    add_json_flag(sub)

    sub = subparsers.add_parser("route-c-provider-comparison", help="Route C (issue #225): run route-c-benchmark once per provider (e.g. lexical/lsa/hybrid) over the same frozen known-gene cases and report Family/Gene Recall@K side by side. Never restructures route-c-benchmark's own scoring and never pools a provider's raw retrieval_score.")
    sub.add_argument("--provider-manifest", action="append", required=True, help="LABEL=PATH for a route-c-benchmark manifest whose retrieval_tsv values are that provider's own retrieval for the same cases. Repeat at least twice.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--full-hash", action="store_true")
    add_json_flag(sub)

    sub = subparsers.add_parser("route-c-direct-search", help="Route C direct rescue (#640): interval-blind Europe PMC search for directly reported trait-gene/locus/variant relations across explicit same/related/model-plant phenotype scopes. Retrieval rows are candidates only, never support.")
    sub.add_argument("--phenotype", required=True)
    sub.add_argument("--focal-species", required=True)
    sub.add_argument("--species-scope", required=True, help="TSV: query_id/reference_species/species_distance/phenotype_query. Similar phenotype wording may differ by species; no interval fields are permitted.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--cache", help="Optional query-keyed Europe PMC cache for offline reproducibility.")
    sub.add_argument("--allow-external-api", action="store_true", help="Send phenotype query wording and reference species to Europe PMC. No interval/candidate/sequence is sent.")
    sub.add_argument("--limit", type=int, default=10)
    sub.add_argument("--case-id", default="")
    sub.add_argument("--full-hash", action="store_true")
    add_json_flag(sub)

    sub = subparsers.add_parser("route-c-direct-freeze", help="Route C direct rescue (#640): machine-verify agent/researcher extracted direct reports with the generic citation/claim gate, then freeze an interval-blind direct-report artifact.")
    sub.add_argument("--proposals", required=True, help="TSV following templates/route_c_direct_proposals.tsv. QTL/interval/candidate/truth columns are rejected.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--citation-cache", help="Local Europe PMC-shaped records for offline source verification.")
    sub.add_argument("--allow-external-api", action="store_true", help="Resolve citation identifiers through Europe PMC/NCBI. Only identifiers are sent by this stage.")
    sub.add_argument("--claim-extraction", help="Optional source-quote extraction artifact; every quote is checked against retrieved text.")
    sub.add_argument("--full-hash", action="store_true")
    add_json_flag(sub)

    sub = subparsers.add_parser("route-c-direct-project", help="Route C direct rescue (#640): after freeze/hash verification, project eligible direct records into an interval atlas and optional existing Route-A eligible BLAST-hit artifact. Homology is projection evidence, not automatic trait transfer.")
    sub.add_argument("--frozen", required=True, help="route_c_direct_rescue_frozen.tsv from route-c-direct-freeze.")
    sub.add_argument("--freeze-manifest", required=True, help="route_c_direct_rescue_freeze_manifest.json; the frozen TSV hash is rechecked before projection.")
    sub.add_argument("--atlas", required=True, help="interval_gene_atlas.tsv.")
    sub.add_argument("--qtl-id", required=True)
    sub.add_argument("--eligible-hits", help="Existing blast-panel-search eligible_interval_hits.tsv; panel_id must equal Route C record_id.")
    sub.add_argument("--identity-mapping", help="Optional source-backed typed gene-to-protein mapping; authorizes review only.")
    sub.add_argument("--reference-bridge", help="Hash-verified post-freeze sequence bridge, kept review-only; mutually exclusive with --eligible-hits.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--full-hash", action="store_true")
    add_json_flag(sub)

    sub = subparsers.add_parser("route-c-direct-benchmark", help="Score frozen Route C direct source recovery and post-freeze projection with separate denominators. Contaminated/not-evaluated cases are excluded, never misses.")
    sub.add_argument("--manifest", required=True, help="TSV following templates/route_c_direct_benchmark.tsv; truth is opened only after freeze.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--full-hash", action="store_true")
    add_json_flag(sub)

    sub = subparsers.add_parser("route-c-direct-members", help="Verify Direct-C freeze/projection lineage and split same-species direct support from cross-species Wet review evidence.")
    sub.add_argument("--frozen", required=True)
    sub.add_argument("--freeze-manifest", required=True)
    sub.add_argument("--projection", required=True)
    sub.add_argument("--projection-manifest", required=True)
    sub.add_argument("--atlas", required=True)
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)

    sub = subparsers.add_parser("phenotype-variant-benchmark", help="Score frozen per-formulation lane artifacts for phenotype wording / description-granularity robustness: Route A family overlap, Route C Recall@K and retrieval overlap, Route B overlay stability with an atlas-fingerprint check, A/B four-states, Stage-A stability, and stable / rank-movement / catastrophic-family-dropout / stage-A-dropout classification. Never re-runs a lane and never changes a score.")
    sub.add_argument("--manifest", required=True, help="TSV: one row per case (case_id/truth_family_id/variants_tsv + optional truth_gene_id/base_phenotype and per-variant artifact patterns route_a_families_tsv/route_c_retrieval_tsv/route_b_overlay_status_json/route_b_atlas_manifest_json/route_b_families_tsv/ranked_tsv; {variant_id} is substituted, else <stem>_<variant_id>.<ext>). variants.tsv rows: variant_id/variant_kind/phenotype_text.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--full-hash", action="store_true")
    add_json_flag(sub)


def _run_retrieval_import(args: Any) -> dict:
    return retrieval_import(args.retrieval, args.out_dir, args.phenotype, args.case_id)


def _run_retrieval_merge(args: Any) -> dict:
    return retrieval_merge(args.retrieval, args.out_dir, args.top_k, args.case_id)


def _run_retrieval_interval_match(args: Any) -> dict:
    return retrieval_interval_match(args.retrieval, args.atlas, args.out_dir, args.case_id)


def _run_triple_route_compare(args: Any) -> dict:
    return triple_route_compare(args.route_a, args.route_b, args.route_c, args.out_dir, args.case_id,
                                args.route_c_interval_match)


def _run_route_c_benchmark(args: Any) -> dict:
    return route_c_benchmark(args.manifest, args.out_dir, args.full_hash)


def _parse_label_path_pairs(values: list[str], flag: str) -> dict[str, str]:
    pairs: dict[str, str] = {}
    for value in values:
        label, separator, path = str(value).partition("=")
        if not separator or not label.strip() or not path.strip():
            raise ValueError(f"{flag} must be LABEL=PATH, got: {value}")
        if label in pairs:
            raise ValueError(f"{flag} label must be unique: {label}")
        pairs[label] = path
    return pairs


def _run_route_c_provider_comparison(args: Any) -> dict:
    provider_manifests = _parse_label_path_pairs(args.provider_manifest, "--provider-manifest")
    return route_c_provider_comparison(provider_manifests, args.out_dir, args.full_hash)


def _run_route_c_direct_search(args: Any) -> dict:
    return direct_literature_search(args.phenotype, args.focal_species, args.species_scope, args.out_dir,
                                    args.cache, args.allow_external_api, args.limit, args.case_id, args.full_hash)


def _run_route_c_direct_freeze(args: Any) -> dict:
    return direct_rescue_freeze(args.proposals, args.out_dir, args.citation_cache,
                                args.allow_external_api, args.claim_extraction, args.full_hash)


def _run_route_c_direct_project(args: Any) -> dict:
    extra = {"reference_bridge_arg": args.reference_bridge} if getattr(args, "reference_bridge", None) else {}
    return direct_rescue_project(args.frozen, args.freeze_manifest, args.atlas, args.qtl_id,
                                 args.out_dir, args.eligible_hits, args.full_hash, args.identity_mapping,
                                 **extra)


def _run_route_c_direct_benchmark(args: Any) -> dict:
    return direct_rescue_benchmark(args.manifest, args.out_dir, args.full_hash)


def _run_route_c_direct_members(args: Any) -> dict:
    return build_route_c_direct_members(
        args.frozen, args.freeze_manifest, args.projection, args.projection_manifest,
        args.atlas, args.out_dir, case_id=args.case_id,
    )


def _run_phenotype_variant_benchmark(args: Any) -> dict:
    return phenotype_variant_benchmark(args.manifest, args.out_dir, args.full_hash)


def _run_corpus_mask(args: Any) -> dict:
    return corpus_mask(args.corpus, args.mask_spec, args.tier, args.out_dir, args.case_id, args.full_hash)


def _run_retrieval_index_build(args: Any) -> dict:
    build = _INDEX_BUILD_PROVIDERS[args.provider]
    return build(args.corpus, args.out_dir, args.database_version, args.corpus_label, args.full_hash)


def _run_retrieval_run(args: Any) -> dict:
    run = _RUN_PROVIDERS[args.provider]
    if args.phenotype_context and args.provider != PROVIDER_LEXICAL:
        # Fail closed rather than silently dropping the audit record: the
        # structured_context block is written by the lexical provider only.
        raise ValueError("--phenotype-context is recorded by the lexical retrieval provider only; "
                         "use --provider lexical (default) or omit the flag")
    if args.provider == PROVIDER_LEXICAL:
        return run(args.index, args.phenotype, args.out_dir, args.top_k, args.case_id,
                   phenotype_context=args.phenotype_context)
    return run(args.index, args.phenotype, args.out_dir, args.top_k, args.case_id)


def _run_retrieval_hybrid_merge(args: Any) -> dict:
    return retrieval_hybrid_merge(args.retrieval, args.out_dir, args.top_k, args.case_id)


def _run_reference_corpus_build(args: Any) -> dict:
    return reference_corpus_build(args.out_dir, args.corpus_version, args.blastdb, args.dump, args.payload,
                                  args.citation_cache, args.benchmark_mode, args.full_hash,
                                  args.uniprot, args.orthogroup_map, args.preset_acceptance)


def _run_reference_corpus_fetch(args: Any) -> dict:
    return uniprot_fetch(args.taxon, args.out_dir, args.allow_external_api, args.include_unreviewed,
                         args.max_rows, args.page_size, args.preset)


def _run_reference_corpus_compare(args: Any) -> dict:
    return reference_corpus_compare(args.corpus, args.out_dir, args.baseline_label,
                                    args.retrieval_manifest, args.full_hash)


HANDLERS = {
    "reference-corpus-fetch": _run_reference_corpus_fetch,
    "reference-corpus-build": _run_reference_corpus_build,
    "reference-corpus-compare": _run_reference_corpus_compare,
    "corpus-mask": _run_corpus_mask,
    "retrieval-index-build": _run_retrieval_index_build,
    "retrieval-run": _run_retrieval_run,
    "retrieval-merge": _run_retrieval_merge,
    "retrieval-hybrid-merge": _run_retrieval_hybrid_merge,
    "retrieval-import": _run_retrieval_import,
    "retrieval-interval-match": _run_retrieval_interval_match,
    "triple-route-compare": _run_triple_route_compare,
    "route-c-benchmark": _run_route_c_benchmark,
    "route-c-provider-comparison": _run_route_c_provider_comparison,
    "route-c-direct-search": _run_route_c_direct_search,
    "route-c-direct-freeze": _run_route_c_direct_freeze,
    "route-c-direct-project": _run_route_c_direct_project,
    "route-c-direct-benchmark": _run_route_c_direct_benchmark,
    "route-c-direct-members": _run_route_c_direct_members,
    "phenotype-variant-benchmark": _run_phenotype_variant_benchmark,
}
