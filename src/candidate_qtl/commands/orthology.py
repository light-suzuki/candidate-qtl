"""Orthology / homology commands (issue #144 phase 3).

Provider concordance, the reference-candidate table, reciprocal-BLAST planning
and import, and OrthoFinder / OMAmer / synteny / focal-paralog imports.
"""

from __future__ import annotations

from argparse import _SubParsersAction
from typing import Any

from ..evidence import orthology_import, orthology_plan
from ..core import runtime_path
from ..family_index import (
    DEFAULT_FAMILY_HITS, DEFAULT_MAX_Q, DEFAULT_TANDEM_DISTANCE,
    family_index, pathway_interval_falsify, route_a_interval_family_match,
)
from ..orthology_staged import (
    focal_paralog_import,
    omamer_import,
    orthofinder_import,
    orthology_tier1_plan,
    synteny_import,
)
from ..family_to_copy import DEFAULT_TANDEM_DISTANCE as COPY_TANDEM_DISTANCE, family_to_copy
from ..dual_route import dual_route_compare, interval_reverse_lookup, phenotype_first_status
from ..narrowing import integrated_narrowing
from ..homology_concordance import homology_concordance
from ..interval_atlas import interval_atlas
from ..reference_candidates import reference_candidate_table
from ..reverse_literature import reverse_literature
from ..mutant_reference_bridge import OUTPUT_NAME, build_mutant_reference_functions
from ..route_b_prepare import route_b_prepare
from ._shared import add_json_flag


def register(subparsers: _SubParsersAction) -> None:
    sub = subparsers.add_parser("family-index", help="Build the cached genome-wide family index for one reference genome from a BLAST/DIAMOND hit table. Run once per reference; per-case falsification is then table lookup only.")
    sub.add_argument("--blast-hits", required=True, help="Reference proteome vs model/relative proteome, 12/13-column outfmt 6.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--genes", help="Gene coordinate table (gene_id/chrom/start/end); genes without a qualifying hit are kept as an explicit blind spot.")
    sub.add_argument("--id-map", help="Protein/transcript to gene id map from extract.")
    sub.add_argument("--family-map", help="Optional curated family labels (subject_id/family_id[/symbol]); overrides the symbol-stem heuristic.")
    sub.add_argument("--reference-label", default="", help="Focal reference genome label recorded in the manifest.")
    sub.add_argument("--reference-species", default="")
    sub.add_argument("--family-hits", type=int, default=DEFAULT_FAMILY_HITS, help="Top hits kept per gene for the coarse family bucket.")
    sub.add_argument("--min-identity", type=float, default=25.0)
    sub.add_argument("--min-query-coverage", type=float, default=0.5)
    sub.add_argument("--max-evalue", type=float, default=1e-5)
    sub.add_argument("--case-id", default="", help="Optional case identifier recorded in the manifest.")
    sub.add_argument("--full-hash", action="store_true")
    add_json_flag(sub)

    sub = subparsers.add_parser("interval-atlas", help="Structure every QTL-interval gene by family / domain / orthogroup / broad function without using any phenotype. An exploration index for forming hypotheses from what the interval actually contains; it never scores or promotes a gene.")
    sub.add_argument("--candidates", required=True, help="Interval candidates TSV (gene_id + optional chrom/start/end/protein_length).")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--family-index", help="family_index.tsv; reused rather than recomputed.")
    sub.add_argument("--domains", help="Domain table (gene_id/domain_id[/domain_name/start/source]).")
    sub.add_argument("--orthology", help="Orthology table (gene_id[/orthogroup/relationship/reference_gene_id/provider]).")
    sub.add_argument("--annotations", help="Reference annotation table (gene_id[/description/function_description/go_terms/protein_family]).")
    sub.add_argument("--flags", help="candidate_flags.tsv, for dark_state and tandem context.")
    sub.add_argument("--function-bins", help="Override the versioned broad-function bin map (templates/broad_function_bins.tsv).")
    sub.add_argument("--tandem-distance", type=int, default=100000)
    sub.add_argument("--dark-max-aa", type=int, default=200)
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)

    sub = subparsers.add_parser("family-to-copy", help="Compare the members of a recovered gene family against each other and record what does and does not discriminate between the copies. A first-Wet priority layer, never a causal claim.")
    sub.add_argument("--families", required=True, help="family_process_hypotheses.tsv (or any TSV with hypothesis_id/gene_family/member_interval_genes).")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--atlas", help="interval_gene_atlas.tsv, for coordinates and domain architecture (canonical <out-dir>/interval_gene_atlas.tsv from route-b-prepare or the root-level interval-atlas file).")
    sub.add_argument("--orthology", help="Orthology table (gene_id[/relationship/provider/evidence_chain/reference_gene_id/orthogroup/paralog_count]).")
    sub.add_argument("--synteny", help="Synteny table (gene_id/synteny_status[/provider]).")
    sub.add_argument("--expression", help="Expression table (gene_id/relevant_expression_status); unmeasured is never a negative.")
    sub.add_argument("--position", help="Positional evidence table (gene_id/position_class); broad interval membership does not discriminate.")
    sub.add_argument("--subfamily-map", help="Curated subfamily labels (subfamily + gene_id/reference_gene_id/family_id); otherwise an orthogroup is used, or the answer stays subfamily_unresolved.")
    sub.add_argument("--family-target", help="Verified target specification (target_id + optional hypothesis_id/expected_family_id/expected_orthogroup/reference_gene_id/expected_subfamily/expected_domain_architecture/expected_synteny_anchor/expected_expression_context/source_ids/verification_status). Copy selection is scored against this target; without it the family-internal modal subfamily/domain are descriptive statistics and outlier QC only and never pick the preferred copy. verification_status other than verified/source_backed/curated is treated as pending_review and can never be the sole basis of a preferred copy.")
    sub.add_argument("--bridge", help="family_candidate_bridge.tsv (canonical 22-column bridge, e.g. from the family-candidate-bridge command). Looked up per (gene_id, family_id) pair, never by gene_id alone, so a gene's bridge row under a different family never leaks into this family's gate. Supplies mapping tier / dependency cluster / ambiguity / route support: members with a tier5/tier6 bridge tier are never promoted into the shortlist (recorded as bridge_promotion_gate=blocked_tier5/blocked_tier6, disambiguation status unchanged), and family-level route rescue (route_rescue_state: rescue_a_only/rescue_b_only/rescue_c_only/multi_route_supported/not_evaluated, the last when every bridge row a family has is tier6_not_evaluated) and bridge-backed copy-selection loss (copy_selection_loss) are recorded as separate columns.")
    sub.add_argument("--route-a-blast", help="blast_evidence.tsv from blast-panel-search. Only evaluable verified rows with ambiguity_status=unique and threshold-passing top hits become a route_a_blastp homology layer for the hit member gene; provenance is recorded in the discrimination column copy_set_evidence_source (present only when this input is supplied).")
    sub.add_argument("--route-a-blast-hits", "--route-a-blast-eligible", dest="route_a_blast_hits",
                     help="eligible_interval_hits.tsv from blast-panel-search; all threshold-passing interval targets retained for copy review (#532).")
    sub.add_argument("--tandem-distance", type=int, default=COPY_TANDEM_DISTANCE)
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)

    sub = subparsers.add_parser("integrated-narrowing", help="Issue #303: compress a ranked candidate set into an ordinal Wet-review order by combining family recovery (convergence), copy discrimination (family-to-copy), and survival state (candidate-elimination). Every ranked column is passed through verbatim; narrowing_priority is the global 1..N review order (tier-major), ordinal tiers N1-N4 only, never a numeric score, rank_contribution stays none, and no candidate is ever excluded here (elimination alone may exclude).")
    sub.add_argument("--ranked", required=True, help="ranked_candidates.tsv from evidence-rank or candidate-pipeline; all columns are passed through verbatim.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--family-recovery", help="family_to_copy_input.tsv (canonical triple-route handoff: hypothesis_id/gene_family/member_interval_genes) resolving gene -> family.")
    sub.add_argument("--convergence", help="triple_route_hypothesis_comparison.tsv (full per-family recovery states) or triple_route_convergence.json (degraded: convergent families only).")
    sub.add_argument("--family-copy", help="family_copy_discrimination.tsv; disambiguation_status is passed through as copy_state.")
    sub.add_argument("--family-copy-shortlist", help="family_copy_shortlist.tsv; first_wet_priority is passed through when present.")
    sub.add_argument("--remaining", help="remaining_candidate_set.tsv from candidate-elimination; survival-state columns only.")
    sub.add_argument("--case-id", default="")
    sub.add_argument("--family-copy-manifest", help="family_copy_manifest.json from family-to-copy.")
    add_json_flag(sub)

    sub = subparsers.add_parser("phenotype-first-status", help="Route A (forward lookup): normalize the phenotype-first artifacts and detect the failure signals that should trigger the interval-first route.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--hypotheses", help="llm_pathway_hypotheses.tsv from the phenotype-first proposal route.")
    sub.add_argument("--validation", help="proposal_validation.tsv from hypothesis-verify.")
    sub.add_argument("--panel", help="verified_bait_panel.tsv from bait-panel-freeze.")
    sub.add_argument("--falsification", help="hypothesis_interval_falsification.tsv from pathway-interval-falsify.")
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)

    sub = subparsers.add_parser("interval-reverse-lookup", help="Route B (reverse lookup): start from the interval's genes, go out to their model-species homologs, and derive family/process hypotheses from what is actually there. Interval-informed, so never scored as blind discovery.")
    sub.add_argument("--atlas", required=True, help="interval_gene_atlas.tsv: the canonical <out-dir>/interval_gene_atlas.tsv from route-b-prepare, or the root-level file from interval-atlas.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--reference-function", help="Reference gene function table (reference_gene_id[/function_description/pathway/mutant_phenotype/source_id]).")
    sub.add_argument("--llm-proposals", help="llm_pathway_hypotheses.tsv imported from an LLM run of the emitted reverse-lookup prompt.")
    sub.add_argument("--phenotype", default="", help="Optional phenotype overlay; recorded in the status, and it makes the lane non-blind.")
    sub.add_argument("--shared-route", help="Optional pre-fixed phenotype->pathway route TSV (phenotype_token/pathway_family/annotation_keywords). Its pathway family and keywords form a separately reported, non-scoring B1 overlay vocabulary.")
    sub.add_argument("--min-members", type=int, default=1)
    sub.add_argument("--min-semantic-coverage", type=float, default=0.0, help="Minimum fraction (0-1) of interval genes with family/orthogroup/domain/broad-function/source-backed-function annotation. At least one annotated gene is always required for semantic evaluability.")
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)

    sub = subparsers.add_parser("route-b-prepare", help="Build the local Route B provider stack (homology search -> orthology table, cached source-backed reference functions, interval atlas) and run the phenotype-free B0 reverse lookup. Reports readiness: a broad-bin-only run is not_evaluated_for_reverse_lookup, never an equivalent Route B result. The canonical atlas handoff is <out-dir>/interval_gene_atlas.tsv.")
    sub.add_argument("--candidates", required=True, help="Interval candidate table (gene_id ...), e.g. candidates.tsv from extract.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--proteins", help="Candidate protein FASTA (candidate_proteins.fa) for the local homology search. When --family-index is absent, the same standard 13-column BLAST output generates an explicitly interval-scoped fallback family index.")
    sub.add_argument("--protein-db", help="Local BLASTP database used when interval reference/annotation material is incomplete; defaults to CANDIDATE_QTL_PROTEIN_DB. Complete usable material skips the search.")
    sub.add_argument("--reference-function", action="append", help="Local reference TSV or directory of UniProt JSON records; repeatable. Annotation and citation links remain investigation material, not verified claims. Defaults to CANDIDATE_QTL_REFERENCE_FUNCTION.")
    sub.add_argument("--mutant-functional-evidence", action="append", help="normalized_functional_evidence.tsv from tnt1-evidence-normalize or the LORE1 adapter (accession-native FUNCTION_COLUMNS rows). Projected onto the reference-function contract via --cross-accession-map; repeatable. Rows without published/genotyped mutant support are skipped and counted, never turned into negatives.")
    sub.add_argument("--cross-accession-map", action="append", help="Explicit cross-accession mapping table (source_full_id/target_full_id/mapping_basis/relationship/mapping_state ...). mapping_basis must be machine-checkable (explicit_id_map/pangene_family/synteny/reciprocal_sequence/composite); symbol-only identity transfer is rejected. Repeatable; without it every emitted row keeps its native ID with mapping_state=not_evaluated.")
    sub.add_argument("--orthology", help="Pre-computed orthology table. Usable interval mappings are retained; local homology search fills missing material when a protein query and database are available. Empty/off-interval rows do not suppress fallback.")
    sub.add_argument("--domains", help="Domain table (gene_id/domain_id[/domain_name/source/start]).")
    sub.add_argument("--family-index", help="Reusable genome-wide family_index.tsv from family-index. Preferred over the interval-scoped fallback generated from --proteins + --protein-db.")
    sub.add_argument("--annotations", help="Per-gene annotation table (gene_id/description/go_terms...).")
    sub.add_argument("--case-id", default="")
    sub.add_argument("--threads", type=int, default=4)
    sub.add_argument("--min-identity", type=float, default=30.0, help="Minimum percent identity for a usable local homology hit.")
    sub.add_argument("--min-coverage", type=float, default=0.5, help="Minimum query coverage (0-1) for a usable local homology hit.")
    sub.add_argument("--max-evalue", type=float, default=1e-5)
    sub.add_argument("--no-auto-local-db", dest="auto_local_db", action="store_false", help="Disable discovery of the configured local protein DB.")
    sub.add_argument("--phenotype", default="", help="Marks which entries of the (phenotype-independent) interval inventory look related. The atlas is always built without it, so enumeration precedes filtering.")
    sub.set_defaults(auto_local_db=True)
    add_json_flag(sub)

    sub = subparsers.add_parser("reverse-literature", help="Route B stage 2: for every interval gene's reference homolog, query the literature focused by the phenotype. A returned article is a link candidate, not evidence; an offline query is not_evaluated, never 'no literature exists'.")
    sub.add_argument("--lookup", required=True, help="interval_reverse_lookup.tsv from route-b-prepare (gene_id/reference_gene[/reference_function]).")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--phenotype", default="", help="Focuses each query. Only generic phenotype words are transmitted.")
    sub.add_argument("--cache", help="Query-keyed cached Europe PMC results, so a rerun is byte-identical offline.")
    sub.add_argument("--allow-external-api", action="store_true", help="Permit read-only Europe PMC queries.")
    sub.add_argument("--limit", type=int, default=5, help="Articles kept per reference gene.")
    sub.add_argument("--case-id", default="")
    sub.add_argument("--full-hash", action="store_true")
    add_json_flag(sub)

    sub = subparsers.add_parser("dual-route-compare", help="Compare the independently generated forward and reverse routes and record the convergence tier: strict mechanistic / family / gene overlap / process-related / text-only, plus phenotype-only / interval-only / discordant / both-weak. Convergence adds no score.")
    sub.add_argument("--route-a", required=True, help="phenotype_first_hypotheses.tsv")
    sub.add_argument("--route-b", required=True, help="family_process_hypotheses.tsv")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--case-id", default="")
    sub.add_argument("--route-b-manifest", help="interval_first_status.json from interval-reverse-lookup; records whether Route B saw the phenotype (B0 phenotype-blind = false, B1 phenotype overlay = true).")
    sub.add_argument("--benchmark-mode", action="store_true", help="Mark the artifact as a benchmark-mode comparison (B0 phenotype-blind Route B).")
    add_json_flag(sub)

    sub = subparsers.add_parser("pathway-interval-falsify", help="Falsify pathway/gene-family hypotheses against the QTL interval, coarse bucket first and the fine ortholog level only for survivors. Enrichment filters hypotheses; it never scores a gene.")
    sub.add_argument("--family-index", required=True, help="family_index.tsv from family-index.")
    sub.add_argument("--hypotheses", required=True, help="llm_pathway_hypotheses.tsv, or any TSV with hypothesis_id + gene_family/reference_gene_name/proposed_stable_id.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--annotation-rescue", help="annotation_rescue.tsv (miniprot projection of the reference proteins onto the interval sequence). Without it an interval miss stays absent_annotation_only -- a hold -- because 'not in the annotated gene list' is not 'absent'.")
    sub.add_argument("--interval", action="append", help="CHROM:START-END resolved against the index coordinates; repeatable.")
    sub.add_argument("--candidates", help="Interval candidates TSV, as an alternative to --interval.")
    sub.add_argument("--tandem-distance", type=int, default=DEFAULT_TANDEM_DISTANCE, help="Same-family genes closer than this collapse into one block.")
    sub.add_argument("--max-q", type=float, default=DEFAULT_MAX_Q, help="BH-adjusted threshold separating enriched from uninformative survivors.")
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)

    sub = subparsers.add_parser("route-a-interval-match", help="Route A stage (issue #355): forward-match only the verified/source-backed baits frozen by bait-panel-freeze to concrete genes inside the QTL interval via a reference_genome-scoped family index. Wires pathway-interval-falsify with verified targets only; zero verified baits, a missing index or an unresolved interval degrade to not_evaluated, never to a miss.")
    sub.add_argument("--panel", required=True, help="verified_bait_panel.tsv from bait-panel-freeze; only its hypothesis_ids are matched.")
    sub.add_argument("--hypotheses", required=True, help="llm_pathway_hypotheses.tsv from mechanism-import.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--family-index", help="reference_genome-scoped family_index.tsv; absent -> not_evaluated_missing_reference_genome_family_index.")
    sub.add_argument("--interval", action="append", help="CHROM:START-END resolved against the index coordinates; repeatable.")
    sub.add_argument("--candidates", help="Interval candidates TSV, as an alternative to --interval.")
    sub.add_argument("--annotation-rescue", help="annotation_rescue.tsv; a rescue provider that did not complete holds absence rather than falsifying it.")
    sub.add_argument("--tandem-distance", type=int, default=DEFAULT_TANDEM_DISTANCE)
    sub.add_argument("--max-q", type=float, default=DEFAULT_MAX_Q)
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)

    sub = subparsers.add_parser("homology-concordance", help="Compare BLAST/DIAMOND/MMseqs2, domain, and orthology providers without vote inflation.")
    sub.add_argument("--hit-table",action="append",required=True,help="PROVIDER=12/13-column outfmt TSV; repeatable"); sub.add_argument("--out-dir",required=True); sub.add_argument("--top-n",type=int,default=5)
    sub.add_argument("--domain-table",action="append",help="PROVIDER=TSV with gene_id/domain_id/start/end"); sub.add_argument("--orthology-table",action="append",help="PROVIDER=TSV with gene_id/orthogroup/relationship")
    sub.add_argument("--candidates"); sub.add_argument("--provider-metadata",action="append",help="PROVIDER=JSON with tool/database version and command"); add_json_flag(sub)

    sub = subparsers.add_parser("reference-candidate-table", help="Create a review table linking QTL candidates to source-backed reference genes through orthology mappings.")
    sub.add_argument("--candidates", required=True, help="Candidate or ranked TSV with gene_id.")
    sub.add_argument("--reference-genes", required=True, help="Source-backed TSV following templates/reference_gene_panel.tsv.")
    sub.add_argument("--orthology", action="append", required=True, help="PROVIDER=TSV following templates/candidate_reference_orthology.tsv; repeatable.")
    sub.add_argument("--annotation", help="Optional TSV with gene_id and annotation columns.")
    sub.add_argument("--out-dir", required=True)
    add_json_flag(sub)

    sub = subparsers.add_parser("orthology-plan", help="Generate reciprocal-BLAST refinement only for Top/Rescue candidates.")
    sub.add_argument("--ranked", required=True)
    sub.add_argument("--protein-fasta", required=True, help="Candidate proteins normalized to gene_id.")
    sub.add_argument("--target-proteome", required=True, help="Complete target-species proteome with gene-compatible IDs.")
    sub.add_argument("--relative-proteome", action="append", required=True, help="LABEL=FASTA; repeatable for near relatives.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--top-n", type=int, default=7)
    sub.add_argument("--num-threads", type=int, default=4)
    add_json_flag(sub)

    sub = subparsers.add_parser("orthology-import", help="Convert reciprocal-BLAST refinement results into canonical evidence rows.")
    sub.add_argument("--manifest", required=True, help="orthology_plan.json after running its shell script.")
    sub.add_argument("--out", required=True)
    add_json_flag(sub)

    sub = subparsers.add_parser("orthology-tier1-plan", help="Profile-aware Tier-1 plan: candidate query vs focal self search and a reference panel.")
    sub.add_argument("--candidates", required=True, help="Candidate table; only its gene_id proteins are used as query.")
    sub.add_argument("--protein-fasta", required=True, help="Candidate proteins normalized to gene_id.")
    sub.add_argument("--focal-proteome", required=True, help="Whole focal-species proteome (search target, reused across QTLs).")
    sub.add_argument("--reference-proteome", action="append", default=[], help="LABEL=FASTA reference panel proteome; repeatable.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--profile", choices=["accuracy64", "laptop16"], default="laptop16")
    sub.add_argument("--profiles-table", help="Optional runtime_profiles.tsv overriding the built-in defaults.")
    sub.add_argument("--aligner", choices=["diamond", "blastp"], default="diamond")
    sub.add_argument("--threads", type=int, help="Override the profile thread count.")
    sub.add_argument("--max-hits", type=int, help="Override the profile per-query hit cap.")
    add_json_flag(sub)

    sub = subparsers.add_parser("focal-paralog-import", help="Flag out-of-interval paralog redundancy from a focal self-search.")
    sub.add_argument("--focal-blast", required=True, help="Candidate-vs-focal-proteome BLAST/DIAMOND outfmt 6 result.")
    sub.add_argument("--candidates", required=True, help="Candidate table defining the QTL-interval gene set.")
    sub.add_argument("--out", required=True, help="Canonical evidence TSV output.")
    sub.add_argument("--close-ratio", type=float, default=0.4, help="Paralog bitscore fraction of the self hit to count as close.")
    add_json_flag(sub)

    sub = subparsers.add_parser("orthofinder-import", help="Import an OrthoFinder Orthogroups.tsv and classify 1:1/1:many/many:many.")
    sub.add_argument("--orthogroups", required=True, help="OrthoFinder Orthogroups.tsv.")
    sub.add_argument("--candidates", required=True)
    sub.add_argument("--focal-species", required=True, help="Column name for the focal species in Orthogroups.tsv.")
    sub.add_argument("--out", required=True, help="Canonical evidence TSV output.")
    add_json_flag(sub)

    sub = subparsers.add_parser("omamer-import", help="Import an OMAmer HOG placement table into orthogroup evidence.")
    sub.add_argument("--omamer", required=True, help="OMAmer search output TSV.")
    sub.add_argument("--candidates", required=True)
    sub.add_argument("--out", required=True)
    add_json_flag(sub)

    sub = subparsers.add_parser("synteny-import", help="Import MCScanX/jcvi synteny anchors into synteny evidence.")
    sub.add_argument("--anchors", required=True, help="jcvi .anchors or MCScanX .collinearity file.")
    sub.add_argument("--candidates", required=True)
    sub.add_argument("--out", required=True)
    sub.add_argument("--format", choices=["anchors", "collinearity"], default="anchors")
    add_json_flag(sub)


def _run_homology_concordance(args: Any) -> dict:
    return homology_concordance(args.hit_table, args.out_dir, args.top_n, args.domain_table, args.orthology_table, args.candidates, args.provider_metadata)


def _run_reference_candidate_table(args: Any) -> dict:
    return reference_candidate_table(args.candidates, args.reference_genes, args.orthology, args.out_dir, args.annotation)


def _run_orthology_plan(args: Any) -> dict:
    return orthology_plan(args.ranked, args.protein_fasta, args.target_proteome, args.relative_proteome, args.out_dir, args.top_n, args.num_threads)


def _run_orthology_import(args: Any) -> dict:
    return orthology_import(args.manifest, args.out)


def _run_orthology_tier1_plan(args: Any) -> dict:
    return orthology_tier1_plan(args.candidates, args.protein_fasta, args.focal_proteome, args.reference_proteome, args.out_dir,
                                args.profile, args.profiles_table, args.aligner, args.threads, args.max_hits)


def _run_focal_paralog_import(args: Any) -> dict:
    return focal_paralog_import(args.focal_blast, args.candidates, args.out, args.close_ratio)


def _run_orthofinder_import(args: Any) -> dict:
    return orthofinder_import(args.orthogroups, args.candidates, args.focal_species, args.out)


def _run_omamer_import(args: Any) -> dict:
    return omamer_import(args.omamer, args.candidates, args.out)


def _run_synteny_import(args: Any) -> dict:
    return synteny_import(args.anchors, args.candidates, args.out, args.format)


def _run_family_index(args: Any) -> dict:
    return family_index(
        args.blast_hits, args.out_dir, args.genes, args.id_map, args.family_map,
        args.reference_label, args.reference_species, args.family_hits,
        args.min_identity, args.min_query_coverage, args.max_evalue, args.full_hash, args.case_id,
    )


def _run_pathway_interval_falsify(args: Any) -> dict:
    return pathway_interval_falsify(
        args.family_index, args.hypotheses, args.out_dir, args.interval, args.candidates,
        args.tandem_distance, args.max_q, args.case_id, args.annotation_rescue,
    )


def _run_route_a_interval_match(args: Any) -> dict:
    return route_a_interval_family_match(
        args.panel, args.hypotheses, args.out_dir, args.family_index, args.interval,
        args.candidates, args.annotation_rescue, args.tandem_distance, args.max_q, args.case_id,
    )


def _run_interval_atlas(args: Any) -> dict:
    return interval_atlas(
        args.candidates, args.out_dir, args.family_index, args.domains, args.orthology,
        args.annotations, args.flags, args.function_bins, args.tandem_distance,
        args.dark_max_aa, args.case_id,
    )


def _run_phenotype_first_status(args: Any) -> dict:
    return phenotype_first_status(args.out_dir, args.hypotheses, args.validation, args.panel,
                                  args.falsification, args.case_id)


def _run_interval_reverse_lookup(args: Any) -> dict:
    return interval_reverse_lookup(args.atlas, args.out_dir, args.reference_function, args.llm_proposals,
                                   args.phenotype, args.case_id, args.min_members,
                                   args.min_semantic_coverage, args.shared_route)


def _run_route_b_prepare(args: Any) -> dict:
    reference_function = list(args.reference_function or [])
    bridge_manifest: dict | None = None
    if getattr(args, "mutant_functional_evidence", None):
        out_dir = runtime_path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        bridge_manifest = build_mutant_reference_functions(
            [runtime_path(p) for p in args.mutant_functional_evidence],
            [runtime_path(p) for p in (args.cross_accession_map or [])],
            out_dir,
            case_id=args.case_id,
        )
        reference_function.append(str(out_dir / OUTPUT_NAME))
    result = route_b_prepare(
        args.candidates, args.out_dir, args.proteins, args.protein_db, reference_function,
        args.orthology, args.domains, args.family_index, args.annotations, args.case_id,
        args.threads, args.min_identity, args.min_coverage, args.max_evalue, args.auto_local_db,
        args.phenotype,
    )
    if bridge_manifest is not None:
        result["mutant_reference_bridge"] = bridge_manifest
    return result


def _run_reverse_literature(args: Any) -> dict:
    return reverse_literature(args.lookup, args.out_dir, args.phenotype, args.cache,
                              args.allow_external_api, args.limit, args.case_id, args.full_hash)


def _run_dual_route_compare(args: Any) -> dict:
    return dual_route_compare(args.route_a, args.route_b, args.out_dir, args.case_id,
                              args.route_b_manifest, args.benchmark_mode)


def _run_family_to_copy(args: Any) -> dict:
    return family_to_copy(
        args.families, args.out_dir, args.atlas, args.orthology, args.synteny,
        args.expression, args.position, args.subfamily_map, args.family_target,
        args.tandem_distance, args.case_id, args.bridge, args.route_a_blast,
        args.route_a_blast_hits,
    )


def _run_integrated_narrowing(args: Any) -> dict:
    return integrated_narrowing(
        args.ranked, args.out_dir, args.family_recovery, args.convergence,
        args.family_copy, args.family_copy_shortlist, args.remaining, args.case_id,
        args.family_copy_manifest,
    )


HANDLERS = {
    "integrated-narrowing": _run_integrated_narrowing,
    "family-to-copy": _run_family_to_copy,
    "phenotype-first-status": _run_phenotype_first_status,
    "interval-reverse-lookup": _run_interval_reverse_lookup,
    "route-b-prepare": _run_route_b_prepare,
    "reverse-literature": _run_reverse_literature,
    "dual-route-compare": _run_dual_route_compare,
    "interval-atlas": _run_interval_atlas,
    "family-index": _run_family_index,
    "pathway-interval-falsify": _run_pathway_interval_falsify,
    "route-a-interval-match": _run_route_a_interval_match,
    "homology-concordance": _run_homology_concordance,
    "reference-candidate-table": _run_reference_candidate_table,
    "orthology-plan": _run_orthology_plan,
    "orthology-import": _run_orthology_import,
    "orthology-tier1-plan": _run_orthology_tier1_plan,
    "focal-paralog-import": _run_focal_paralog_import,
    "orthofinder-import": _run_orthofinder_import,
    "omamer-import": _run_omamer_import,
    "synteny-import": _run_synteny_import,
}
