"""Family candidate bridge CLI (deterministic Route A/B/C normalization).

The bridge merges the independently produced route artifacts into one
reporting table with a strict mapping tier. Reporting-only: every row carries
``rank_contribution=none`` and no score column is ever emitted.
"""

from __future__ import annotations

from argparse import _SubParsersAction
from typing import Any

from ..family_candidate_bridge import family_candidate_bridge
from ._shared import add_json_flag


def register(subparsers: _SubParsersAction) -> None:
    sub = subparsers.add_parser(
        "family-candidate-bridge",
        help="Normalize Route A / B / C family evidence (phenotype-first hypotheses, orthology/atlas, frozen retrieval match) into one reporting table with a strict mapping tier (tier1_source_backed_unique .. tier6_not_evaluated). Reporting-only: rank_contribution is always none, there are no score columns, and the family_to_copy / narrowing wiring is untouched. At least one route input is required; a Route B atlas is independently sufficient.",
    )
    sub.add_argument("--route-a-hypotheses", help="phenotype_first_hypotheses.tsv (hypothesis_id + member_interval_genes[/gene_family/proposed_stable_id]).")
    sub.add_argument("--route-a-family-status", help="route_a_interval_family_match.tsv (hypothesis_id + family_status[/interval_gene_ids/source_ids]); per-hypothesis family status from the Route A interval match stage.")
    sub.add_argument("--route-a-bait-panel", help="verified_bait_panel.tsv from bait-panel-freeze (hypothesis_id + function_source_ids[/source_provider]); hypotheses outside the panel are unverified bait, never a match.")
    sub.add_argument("--route-a-projection-panel", help="projection_bait_panel.tsv for identity-confirmed review-only BLAST members; never contributes verified support.")
    sub.add_argument("--route-a-blast", help="blast_evidence.tsv from blast-panel-search (panel_id + top_hit_identity/top_hit_query_coverage/top_hit_target_coverage/top_hit_evalue/paralog_count/ambiguity_status/provider_status); joined to Route A matched rows by reference_accession=panel_id (the verified bait accession). A blast hit is similarity evidence only for an existing matched family row; an eligible hit may create one sequence-match row only when family-index status explicitly falsifies that family (#530).")
    sub.add_argument("--route-a-blast-hits", "--route-a-blast-eligible", dest="route_a_blast_hits",
                     help="eligible_interval_hits.tsv from blast-panel-search; all long-form threshold-passing interval hits consumed by Route A family/copy review (#532).")
    sub.add_argument("--route-b-orthology", help="Orthology table (gene_id/reference_gene_id + relationship/orthogroup/percent_identity/query_coverage/evalue/provider[/source_id]); best_hit_homology is a best hit, not a reciprocal one.")
    sub.add_argument("--route-b-atlas", help="interval_gene_atlas.tsv (gene_id + evolutionary_family_id/orthogroup/domain_architecture/domain_ids/domain_source/family_id/family_label/focal_paralog_count/tandem_expansion); also emits one row per atlas evidence axis.")
    sub.add_argument("--route-b-family-index", help="Route B0 family_index.tsv. An interval_candidates index is consumed only for positive Route A family projection after both lanes freeze: every matching interval copy is retained; no-hit never proves absence.")
    sub.add_argument("--route-c-retrieval", help="route_c_retrieval_merged.tsv (reference_gene_id/rank + database_version/model_version); joined to the interval match on (reference_gene_id, rank).")
    sub.add_argument("--route-c-match", help="route_c_interval_match.tsv (reference_gene_id/rank/matched_interval_gene + match_basis/family_id/orthogroup_id); rows without any family id are skipped and counted.")
    sub.add_argument("--dependency-clusters", help="evidence_dependency_clusters.tsv (dependency_cluster/dependency_status/root_source_id/source_ids); overrides the route-default dependency roots.")
    sub.add_argument("--genes", help="One-column gene_id TSV; genes with no bridge row are reported in the manifest.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--case-id", required=True, help="Non-empty case identifier recorded in every bridge row and the manifest.")
    add_json_flag(sub)


def _run_family_candidate_bridge(args: Any) -> dict:
    return family_candidate_bridge(
        out_dir_arg=args.out_dir,
        route_a_hypotheses_arg=args.route_a_hypotheses,
        route_a_family_status_arg=args.route_a_family_status,
        route_a_bait_panel_arg=args.route_a_bait_panel,
        route_a_projection_panel_arg=args.route_a_projection_panel,
        route_b_orthology_arg=args.route_b_orthology,
        route_b_atlas_arg=args.route_b_atlas,
        route_c_retrieval_arg=args.route_c_retrieval,
        route_c_match_arg=args.route_c_match,
        dependency_clusters_arg=args.dependency_clusters,
        genes_arg=args.genes,
        case_id=args.case_id,
        route_a_blast_arg=args.route_a_blast,
        route_a_blast_eligible_arg=args.route_a_blast_hits,
        route_b_family_index_arg=args.route_b_family_index,
    )


HANDLERS = {"family-candidate-bridge": _run_family_candidate_bridge}
