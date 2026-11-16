"""Same-species assembly bridge commands (Issue #527 / #524)."""

from __future__ import annotations

from argparse import _SubParsersAction
from typing import Any

from ..assembly_bridge import build_assembly_bridge
from ..multi_reference_locus_union import build_multi_reference_locus_union
from ._shared import add_json_flag


def register(subparsers: _SubParsersAction) -> None:
    sub = subparsers.add_parser(
        "assembly-bridge",
        help="Bridge native same-species assembly IDs using optional LIS pangene candidate correspondences, fam3, and explicit pairwise reciprocal sequence evidence; emit mapping and existing Route B TSV adapters without symbol-only matching.",
    )
    sub.add_argument("--functional-evidence", "--source-functional-evidence", "--source-functions", dest="functional_evidence", required=True,
                     help="normalized_functional_evidence.tsv from the native resource snapshot.")
    sub.add_argument("--source-fam3", "--source-family-assignments", dest="source_fam3", required=True,
                     help="Source native-ID fam3 assignment TSV (native_gene_id/gene + family_id/fam3).")
    sub.add_argument("--target-fam3", "--target-family-assignments", dest="target_fam3", required=True,
                     help="Target native-ID fam3 assignment TSV (native_gene_id/gene + family_id/fam3).")
    sub.add_argument("--sequence-evidence", "--pairwise-sequence", "--sequences", dest="sequence_evidence", required=True,
                     help="Explicit source_full_id/target_full_id pairwise sequence TSV; reciprocal state is required for verification.")
    sub.add_argument("--pangene-hash", help="Optional frozen LIS pangene *.hsh.tsv[.gz] table (set ID first, gene second; named gene/representative headers also accepted). Membership only proposes full-ID candidate relations; it never verifies one-copy identity.")
    sub.add_argument("--pangene-dataset", default="", help="Exact pangene dataset generation recorded with --pangene-hash, e.g. Medicago.pan2.451K. Required when --pangene-hash is supplied.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--source-namespace", default="", help="Optional declared source namespace; otherwise inferred from each full ID.")
    sub.add_argument("--target-namespace", default="", help="Optional declared target namespace; otherwise inferred from each full ID.")
    sub.add_argument(
        "--mapping-scope",
        default="same_species_cross_assembly",
        choices=["same_species_cross_assembly", "same_accession_annotation_generation", "cross_accession"],
        help=(
            "Mapping scope recorded in every mapping row and the manifest. "
            "cross_accession is for explicit same-species accession-to-accession bridges "
            "such as current Medicago A17 -> R108; it does not weaken identity verification."
        ),
    )
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)

    sub = subparsers.add_parser(
        "multi-reference-locus-union",
        help="Build a full-ID-preserving multi-reference candidate universe and conservative locus-entity union; only verified 1:1 mappings collapse duplicates.",
    )
    sub.add_argument("--primary-candidates", "--primary", dest="primary_candidates", required=True,
                     help="Primary/run-target candidate table.")
    sub.add_argument("--alternate-candidates", "--alternate", dest="alternate_candidates", action="append", default=[], metavar="LABEL=PATH",
                     help="Alternate candidate table; repeat LABEL=PATH. Alternate-only IDs are retained.")
    sub.add_argument("--mapping", "--assembly-mapping", dest="mapping", required=True,
                     help="Full assembly_bridge_mapping.tsv; all states are audited, only safe links collapse.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--primary-label", default="primary")
    sub.add_argument("--case-id", default="")
    add_json_flag(sub)


def _run_assembly_bridge(args: Any) -> dict:
    return build_assembly_bridge(
        args.functional_evidence, args.source_fam3, args.target_fam3,
        args.sequence_evidence, args.out_dir,
        source_namespace=args.source_namespace, target_namespace=args.target_namespace,
        mapping_scope=args.mapping_scope, case_id=args.case_id,
        pangene_hash_arg=args.pangene_hash, pangene_dataset=args.pangene_dataset,
    )


def _run_multi_reference_locus_union(args: Any) -> dict:
    return build_multi_reference_locus_union(
        args.primary_candidates, args.alternate_candidates, args.mapping, args.out_dir,
        primary_label=args.primary_label, case_id=args.case_id,
    )


HANDLERS = {
    "assembly-bridge": _run_assembly_bridge,
    "multi-reference-locus-union": _run_multi_reference_locus_union,
}
