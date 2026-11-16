"""Adapt frozen Route-C direct records into the supported-member ledger.

This is deliberately a ledger stage, rather than a scoring stage. Literature
and projection provenance are retained for review. Canonical support additionally
requires a sealed Route-C-specific claim-binding artifact proving that source
organism, reported identity, experimented gene, perturbation and target phenotype
belong to one claim; generic citation/cue matching alone remains review-only.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from .core import runtime_path, sha256
from .claim_binding import route_c_authorizes
from .route_c_claim_binding import (
    BINDING_CONTRACT,
    binding_authorizes_projection,
    binding_is_complete,
    load_claim_binding_from_freeze_manifest,
)
from .route_c_direct_rescue import (
    validate_direct_rescue_frame,
    validate_projection_frame,
)

SCHEMA_VERSION = "1.1"
CONTRACT = "route-c-direct-members-1.1.0"

OUTPUT_COLUMNS = (
    "candidate_gene_id",
    "record_id",
    "source_id",
    "query_phenotype",
    "reported_phenotype",
    "focal_species",
    "reference_species",
    "species_distance",
    "evidence_directness",
    "evidence_type",
    "source_verification_state",
    "gene_identity_state",
    "trait_relation_state",
    "direct_rescue_eligibility",
    "projection_state",
    "projection_basis",
    "projection_ambiguity_state",
    "source_reference_gene_id",
    "identity",
    "query_coverage",
    "target_coverage",
    "family_ids",
    "orthogroup_ids",
    "claim_binding_status",
    "claim_binding_reason",
    "claim_binding_source_organism_state",
    "claim_binding_identity_state",
    "claim_binding_phenotype_state",
    "claim_binding_perturbation_state",
    "claim_binding_projection_authorized",
    "support_status",
    "membership_effect",
    "rank_contribution",
    "notes",
)

SUPPORTED_STATUS = "verified_direct_support"
REVIEW_STATUS = "reference_phenotype_review"
LEDGER_STATUS = "ledger_only"


def _read_manifest(path_arg: str, label: str) -> tuple[Path, dict[str, object]]:
    path = runtime_path(path_arg)
    if not path.exists():
        raise ValueError(f"{label} does not exist: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is not valid JSON: {path}") from error
    if not isinstance(payload, dict):
        raise ValueError(f"{label} must contain a JSON object")
    return path, payload


def _manifest_sha(manifest: dict[str, object], *keys: str) -> str:
    value: object = manifest
    for key in keys:
        if not isinstance(value, dict):
            return ""
        value = value.get(key, "")
    return str(value or "").removeprefix("sha256:")


def _atlas_values(row: dict[str, object], names: tuple[str, ...]) -> str:
    values: set[str] = set()
    for name in names:
        value = str(row.get(name, "") or "")
        values.update(
            token.strip()
            for token in value.replace(",", ";").split(";")
            if token.strip()
        )
    return ";".join(sorted(values))


def _base_direct_support_eligible(
    record: dict[str, object], projection: dict[str, object]
) -> bool:
    return (
        str(record.get("species_distance", "")).strip()
        in {"S0_same_cultivar_or_accession", "S1_same_species"}
        and str(record.get("evidence_directness", "")).strip()
        in {
            "D0_causal_perturbation_or_rescue",
            "D1_mapped_gene_with_functional_validation",
        }
        and str(record.get("source_verification_state", "")).strip() == "verified"
        and str(record.get("direct_rescue_eligibility", "")).strip()
        == "eligible_direct_record"
        and str(record.get("gene_identity_state", "")).strip()
        in {"verified_gene", "verified_locus_only"}
        and str(record.get("trait_relation_state", "")).strip()
        in {"direct_target_trait", "direct_related_trait"}
        and str(projection.get("projection_state", "")).strip()
        == "inside_current_qtl"
        and str(projection.get("projection_basis", "")).strip()
        in {"exact_focal_id", "assembly_or_alias_mapping"}
        and bool(str(projection.get("focal_gene_id", "")).strip())
    )


def _is_verified_direct(
    record: dict[str, object],
    projection: dict[str, object],
    binding: dict[str, object] | None,
    binding_generation_complete: bool,
) -> bool:
    bound_ok, _bound_state, _bound_reason = route_c_authorizes(binding, projection)
    return (
        binding_generation_complete
        and binding_is_complete(binding)
        and binding_authorizes_projection(binding, projection)
        and bound_ok
        and _base_direct_support_eligible(record, projection)
    )


def build_route_c_direct_members(
    frozen_arg: str,
    freeze_manifest_arg: str,
    projection_arg: str,
    projection_manifest_arg: str,
    atlas_arg: str,
    out_dir_arg: str,
    *,
    case_id: str = "",
) -> dict[str, object]:
    """Build a lineage-checked Direct-C member/review ledger.

    The frozen and projection artifacts are validated independently and their
    hashes are checked against their manifests. The atlas is used only to attach
    family/orthogroup provenance. A missing claim-binding generation is explicit
    ``not_evaluated`` and blocks support promotion without turning into negative
    biological evidence; a corrupt supplied binding artifact is an error.
    """
    freeze_manifest_path, freeze_manifest = _read_manifest(
        freeze_manifest_arg, "Route C freeze manifest"
    )
    projection_manifest_path, projection_manifest = _read_manifest(
        projection_manifest_arg, "Route C projection manifest"
    )
    frozen_tsv = runtime_path(frozen_arg)
    projection_tsv = runtime_path(projection_arg)
    atlas_path = runtime_path(atlas_arg)
    if not frozen_tsv.exists() or not projection_tsv.exists() or not atlas_path.exists():
        raise ValueError("Route C direct-members input artifact is missing")
    expected_frozen = _manifest_sha(freeze_manifest, "frozen_sha256")
    if not expected_frozen or expected_frozen != sha256(frozen_tsv):
        raise ValueError("Route C frozen artifact hash does not match freeze manifest")
    expected_projection = _manifest_sha(projection_manifest, "projection_sha256")
    if not expected_projection or expected_projection != sha256(projection_tsv):
        raise ValueError("Route C projection artifact hash does not match projection manifest")
    projection_qtl = str(projection_manifest.get("qtl_id", "")).strip()
    direct = validate_direct_rescue_frame(
        pd.read_csv(frozen_tsv, sep="\t", dtype=str, keep_default_na=False)
    )
    binding_by_record, binding_meta = load_claim_binding_from_freeze_manifest(
        freeze_manifest,
        freeze_manifest_path,
        set(direct["record_id"].astype(str)),
        {
            str(row["record_id"]): row
            for row in direct.to_dict("records")
        },
    )
    projection = validate_projection_frame(
        pd.read_csv(projection_tsv, sep="\t", dtype=str, keep_default_na=False)
    )
    if case_id and str(freeze_manifest.get("case_id", case_id)).strip() not in {
        "",
        case_id,
    }:
        raise ValueError("Route C freeze manifest case_id does not match requested case_id")
    if projection_qtl and set(projection["qtl_id"].astype(str)) - {"", projection_qtl}:
        raise ValueError("Route C projection rows contain qtl_id different from projection manifest")
    atlas = pd.read_csv(atlas_path, sep="\t", dtype=str, keep_default_na=False)
    if "gene_id" not in atlas.columns:
        raise ValueError("Route C direct-members atlas requires gene_id")
    atlas_by_gene = {
        str(row["gene_id"]).strip(): row for row in atlas.to_dict("records")
    }
    records = {
        str(row["record_id"]): row for row in direct.to_dict("records")
    }
    if set(projection["record_id"].astype(str)) - set(records):
        raise ValueError("Route C projection contains record_id absent from frozen direct table")

    rows: list[dict[str, object]] = []
    blocked_incomplete_binding = 0
    binding_generation_complete = bool(binding_meta.get("complete", False))
    for projected in projection.to_dict("records"):
        record_id = str(projected["record_id"])
        record = records[record_id]
        binding = binding_by_record.get(record_id, {})
        projection_authorized = binding_authorizes_projection(binding, projected)
        state = str(projected.get("projection_state", "")).strip()
        inside = state == "inside_current_qtl" and str(
            projected.get("focal_gene_id", "")
        ).strip()
        focal_gene = (
            str(projected.get("focal_gene_id", "")).strip() if inside else ""
        )
        atlas_row = atlas_by_gene.get(focal_gene, {})
        base_eligible = _base_direct_support_eligible(record, projected)
        if base_eligible and not (
            binding_generation_complete
            and binding_is_complete(binding)
            and projection_authorized
        ):
            blocked_incomplete_binding += 1
        if _is_verified_direct(
            record, projected, binding, binding_generation_complete
        ):
            support_status, membership = SUPPORTED_STATUS, "supported_union"
        elif (
            inside
            and str(record.get("species_distance", "")).strip()
            in {
                "S0_same_cultivar_or_accession",
                "S1_same_species",
                "S2_same_genus",
                "S3_close_clade",
                "S4_broader_legume",
                "S5_model_or_other_plant",
            }
            and str(record.get("evidence_directness", "")).strip()
            != "D3_context_only"
        ):
            support_status, membership = REVIEW_STATUS, "review_only"
        else:
            support_status, membership = LEDGER_STATUS, "ledger_only"
        family_ids = _atlas_values(
            atlas_row, ("family_id", "family_ids", "gene_family")
        )
        if inside and not family_ids:
            source_family = str(
                projected.get("source_reference_gene_id", "")
            ).strip()
            family_ids = f"direct:{source_family or focal_gene}"
        rows.append(
            {
                "candidate_gene_id": focal_gene,
                "record_id": record_id,
                "source_id": record["source_id"],
                "query_phenotype": record["query_phenotype"],
                "reported_phenotype": record["reported_phenotype"],
                "focal_species": record["focal_species"],
                "reference_species": record["reference_species"],
                "species_distance": record["species_distance"],
                "evidence_directness": record["evidence_directness"],
                "evidence_type": record["evidence_type"],
                "source_verification_state": record["source_verification_state"],
                "gene_identity_state": record["gene_identity_state"],
                "trait_relation_state": record["trait_relation_state"],
                "direct_rescue_eligibility": record["direct_rescue_eligibility"],
                "projection_state": state,
                "projection_basis": projected["projection_basis"],
                "projection_ambiguity_state": projected["ambiguity_state"],
                "source_reference_gene_id": projected[
                    "source_reference_gene_id"
                ],
                "identity": projected["identity"],
                "query_coverage": projected["query_coverage"],
                "target_coverage": projected["target_coverage"],
                "family_ids": family_ids,
                "orthogroup_ids": _atlas_values(
                    atlas_row,
                    ("orthogroup_id", "orthogroup_ids", "orthogroup"),
                ),
                "claim_binding_status": str(
                    binding.get("binding_status", "not_evaluated")
                ),
                "claim_binding_reason": str(
                    binding.get(
                        "binding_reason",
                        "structured Route C claim binding was not supplied",
                    )
                ),
                "claim_binding_source_organism_state": str(
                    binding.get("source_organism_state", "not_evaluated")
                ),
                "claim_binding_identity_state": str(
                    binding.get(
                        "reported_identity_verification_state", "not_evaluated"
                    )
                ),
                "claim_binding_phenotype_state": str(
                    binding.get("phenotype_match_state", "not_evaluated")
                ),
                "claim_binding_perturbation_state": str(
                    binding.get(
                        "perturbation_to_phenotype_binding_state",
                        "not_evaluated",
                    )
                ),
                "claim_binding_projection_authorized": projection_authorized,
                "support_status": support_status,
                "membership_effect": membership,
                "rank_contribution": "none",
                "notes": record.get("notes", ""),
            }
        )
    frame = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    output_path = out / "route_c_direct_members.tsv"
    frame.to_csv(output_path, sep="\t", index=False, lineterminator="\n")
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "contract": CONTRACT,
        "stage": "route_c_direct_members",
        "case_id": case_id,
        "qtl_id": projection_qtl,
        "record_count": int(len(frame)),
        "supported_count": int(
            frame["support_status"].eq(SUPPORTED_STATUS).sum()
        ),
        "review_count": int(frame["support_status"].eq(REVIEW_STATUS).sum()),
        "ledger_only_count": int(
            frame["support_status"].eq(LEDGER_STATUS).sum()
        ),
        "one_to_many_count": int(
            frame["projection_ambiguity_state"].eq("one_to_many").sum()
        ),
        "claim_binding": binding_meta,
        "required_claim_binding_contract": BINDING_CONTRACT,
        "support_blocked_incomplete_claim_binding_count": blocked_incomplete_binding,
        "output_sha256": hashlib.sha256(output_path.read_bytes()).hexdigest(),
        "rank_contribution": "none",
        "inputs": {
            "frozen": {
                "path": str(frozen_tsv),
                "sha256": sha256(frozen_tsv),
            },
            "freeze_manifest": {
                "path": str(freeze_manifest_path),
                "sha256": sha256(freeze_manifest_path),
            },
            "projection": {
                "path": str(projection_tsv),
                "sha256": sha256(projection_tsv),
            },
            "projection_manifest": {
                "path": str(projection_manifest_path),
                "sha256": sha256(projection_manifest_path),
            },
            "atlas": {
                "path": str(atlas_path),
                "sha256": sha256(atlas_path),
            },
        },
        "policy": [
            "generic citation resolution and cue/term matching are triage, not a gene/species/trait binding proof",
            "canonical Direct-C support requires a sealed route-c-source-claim-binding-1.0.0 sidecar with every load-bearing axis verified",
            "the source-bound identity must be the same stable ID/locus used by exact or assembly/alias projection; a symbol cannot authorize another submitted ID",
            "missing claim binding is not_evaluated and leaves an otherwise eligible inside-QTL record review-only",
            "a corrupt or roster-mismatched supplied claim-binding artifact fails closed",
            "only S0/S1, D0/D1, verified, eligible, inside-QTL exact or assembly/alias mappings may support after binding passes",
            "S0/S1 D2 and all cross-species records remain reference_phenotype_review",
            "one-to-many projections retain every mapped focal copy and never add a rank score",
            "outside, no-match, and unresolved projections remain ledger-only",
        ],
        "output": str(output_path),
    }
    manifest_path = out / "route_c_direct_members_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        "members": str(output_path),
        "manifest": str(manifest_path),
        "record_count": len(frame),
        "supported_count": manifest["supported_count"],
        "review_count": manifest["review_count"],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frozen", required=True)
    parser.add_argument("--freeze-manifest", required=True)
    parser.add_argument("--projection", required=True)
    parser.add_argument("--projection-manifest", required=True)
    parser.add_argument("--atlas", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--case-id", default="")
    args = parser.parse_args(argv)
    result = build_route_c_direct_members(
        args.frozen,
        args.freeze_manifest,
        args.projection,
        args.projection_manifest,
        args.atlas,
        args.out_dir,
        case_id=args.case_id,
    )
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
