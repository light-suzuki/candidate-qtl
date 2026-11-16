"""Freeze Route-B semantic candidate members before A/B comparison (#403).

Only final ``verification_status=verified`` semantic claims may contribute to
the supported set. Grounded unresolved claims are sealed separately as
explicitly unverified review candidates; they never become support.
"""

from __future__ import annotations

import hashlib
import json

import pandas as pd

from .core import runtime_path
from .semantic_claims import normalize_semantic_claims
from .relevance_tiers import CONTRACT as RELEVANCE_CONTRACT, RELEVANCE_COLUMNS, annotate_relevance

SCHEMA_VERSION = "1.0"
CANDIDATE_COLUMNS = list(dict.fromkeys([
    "candidate_gene_id",
    "family_ids",
    "orthogroup_ids",
    "subfamilies",
    "semantic_unit_ids",
    "claim_ids",
    "evidence_scopes",
    "reference_genes",
    "source_ids",
    "route",
    "support_status",
    "rank_contribution",
    *RELEVANCE_COLUMNS,
]))
REVIEW_COLUMNS = list(dict.fromkeys([
    "candidate_admission_policy", "candidate_admission_state",
    "candidate_gene_id", "claim_id", "semantic_unit_id", "family_id",
    "orthogroup_id", "subfamily_or_clade", "evidence_scope", "reference_gene",
    "source_id", "verification_status", "verification_reason", "route",
    "support_status", "rank_contribution", *RELEVANCE_COLUMNS,
    "reference_function_question_id", "reference_function_review_reason",
]))


def _split(value: object) -> set[str]:
    return {
        token.strip()
        for token in str(value or "").replace(",", ";").split(";")
        if token.strip()
    }


def _join(values: set[str]) -> str:
    return ";".join(sorted(values))


def route_b_homolog_projection(claim: dict, unit: dict, gene_id: str, members: set[str], *,
                               taxonomy_records: dict | None = None) -> str:
    """Use this reference's frozen gene relationship, never a unit/domain label."""
    if gene_id not in members or not any(gene.get("gene_id") == gene_id for gene in unit.get("genes", [])):
        return ""
    from math import isfinite
    from .family_candidate_bridge import _coverage_pass
    from .route_b_semantic_transfer_verify import _endpoint_namespace, _reference_endpoint_matches, _initial_reference_conflicts, _initial_reference_endpoints, _initial_reference_metadata
    from .reference_taxonomy import resolve_species_display
    reference = str(claim.get("reference_gene", ""))
    species = str(claim.get("reference_species", ""))
    gene = next((g for g in unit.get("genes", []) if g.get("gene_id") == gene_id), {})
    cards = unit.get("reference_material_cards", []) if taxonomy_records else []
    if gene_id not in members or _initial_reference_conflicts(_split(reference), _initial_reference_metadata(species),
            _initial_reference_endpoints(gene, cards, expected_species=species, taxonomy_records=taxonomy_records)):
        return ""
    evidence, basis = [], ""
    comparison_gene = gene
    if taxonomy_records:
        orthology = str(gene.get("orthology_reference_gene", "")).strip()
        same_reference = str(gene.get("reference_gene", "")).strip() in {"", orthology}
        taxon = gene.get("orthology_reference_taxon_id", "") or (gene.get("reference_taxon_id", "") if same_reference else "")
        endpoint_species = gene.get("orthology_reference_species", "") or gene.get("orthology_species", "") or (
            gene.get("reference_species", "") if same_reference else "")
        _, audit = resolve_species_display(endpoint_species, species, taxon, taxonomy_records)
        if audit is not None:
            comparison_gene = {**gene, "orthology_reference_species": species}
    if _reference_endpoint_matches(reference, species, comparison_gene):
        evidence = [{key: gene[key] for key in gene if key.startswith("orthology_")}]
        basis = "frozen_gene_orthology"
    else:
        for hit in gene.get("sequence_similarity_hits", []):
            if hit.get("hit_accession") not in _split(reference):
                continue
            requested_species = _initial_reference_metadata(species)
            hit_species = _initial_reference_metadata(hit.get("hit_species", ""))
            if taxonomy_records:
                requested_namespace = _initial_reference_metadata(_endpoint_namespace({}, reference))
                hit_namespace = _initial_reference_metadata(_endpoint_namespace(hit, str(hit.get("hit_accession", ""))))
                if requested_namespace and hit_namespace and requested_namespace != hit_namespace:
                    continue
                hit_species, _ = resolve_species_display(hit.get("hit_species", ""), species,
                    hit.get("hit_taxon_id", ""), taxonomy_records)
            if requested_species and hit_species and requested_species != hit_species:
                continue
            values = (hit.get("percent_identity", ""), hit.get("query_coverage", ""))
            try:
                usable = all(str(value).strip() and isfinite(float(value)) for value in values)
            except (ValueError, TypeError):
                usable = False
            if usable and _coverage_pass(*values, "")[0]:
                evidence.append(hit)
        if evidence:
            basis = "qualified_sequence_hit"
    return json.dumps({"basis": basis, "route": "b", "candidate_gene_id": gene_id,
                       "reference_gene_id": reference, "reference_species": species,
                       "member_genes": sorted(members), "copy_count": len(members),
                       "ambiguity": "multiple_interval_copies" if len(members) > 1 else "not_assessed",
                       "semantic_unit_id": unit.get("unit_id", ""), "evidence": evidence},
                      ensure_ascii=False, sort_keys=True) if basis else ""


def freeze_route_b_semantic_candidates(
    claims_arg: str,
    out_dir_arg: str,
    *,
    case_id: str = "",
) -> dict[str, object]:
    claims_path = runtime_path(claims_arg)
    if not claims_path.exists():
        raise ValueError(f"Route B semantic claims do not exist: {claims_path}")
    claims = normalize_semantic_claims(
        pd.read_csv(claims_path, sep="\t", dtype=str).fillna("")
    )
    if not claims.empty and not claims["route"].eq("interval_first").all():
        raise ValueError("Route B semantic candidate freeze accepts interval_first claims only")

    claims_manifest_path = claims_path.with_name("semantic_claim_transfer_verification_manifest.json")
    claims_manifest = None
    if claims_manifest_path.exists():
        try:
            claims_manifest = json.loads(claims_manifest_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as error:
            raise ValueError("Route B transfer-verification manifest is invalid JSON") from error
        if claims_manifest.get("stage") != "route_b_semantic_transfer_verification":
            raise ValueError("Route B claims were not produced by transfer verification")
        if claims_manifest.get("contract") != "route-b-semantic-transfer-verification-1.0.0":
            raise ValueError("Route B transfer-verification contract is unsupported")
        if str(claims_manifest.get("case_id", "")).strip() != str(case_id).strip():
            raise ValueError("Route B transfer-verification case_id mismatch")

    support_by_claim = {}
    if claims_manifest is not None:
        outputs = claims_manifest.get("outputs", {})
        transfer_path = runtime_path(str(outputs.get("transfer_verification", "")))
        if (claims_manifest.get("artifact_sha256") != hashlib.sha256(claims_path.read_bytes()).hexdigest()
                or not transfer_path.is_file()
                or outputs.get("transfer_verification_sha256") != hashlib.sha256(transfer_path.read_bytes()).hexdigest()):
            raise ValueError("Route B transfer-verification input hash mismatch")
        assessed_support = pd.read_csv(transfer_path, sep="\t", dtype=str).fillna("")
        if not {"claim_id", "supported_member_genes"} <= set(assessed_support.columns):
            raise ValueError("Route B transfer table lacks supported member columns")
        if assessed_support["claim_id"].duplicated().any() or set(assessed_support["claim_id"]) != set(claims["claim_id"]):
            raise ValueError("Route B transfer claim roster mismatch")
        support_by_claim = {row["claim_id"]: _split(row["supported_member_genes"])
                            for row in assessed_support.to_dict("records")}

    projection_units = {}
    bundle_input = ((claims_manifest or {}).get("ancestors") or {}).get("bundle")
    if isinstance(bundle_input, dict):
        projection_bundle = runtime_path(str(bundle_input.get("path", "")))
        if not projection_bundle.is_file() or hashlib.sha256(projection_bundle.read_bytes()).hexdigest() != str(bundle_input.get("sha256", "")).removeprefix("sha256:"):
            raise ValueError("Route B transfer-verification bundle hash mismatch")
        projection_units = {unit["unit_id"]: unit for unit in json.loads(projection_bundle.read_text(encoding="utf-8")).get("units", [])}
    verified = claims[claims["verification_status"].eq("verified")].copy()
    records: dict[str, dict[str, set[str]]] = {}
    supported_relevance_claims = []
    for row in verified.to_dict("records"):
        members = _split(row.get("member_interval_genes", ""))
        if not members:
            raise ValueError(f"verified semantic claim {row.get('claim_id', '')} has no member_interval_genes")
        if support_by_claim:
            supported_members = support_by_claim[str(row["claim_id"])]
            if not supported_members or not supported_members <= members:
                raise ValueError("Route B verified support members are not a non-empty claim subset")
            members = supported_members
        for member in sorted(members):
            supported_relevance_claims.append({**row, "candidate_gene_id": member, "member_interval_genes": _join(members),
                "homolog_projection_json": route_b_homolog_projection(row, projection_units.get(row.get("semantic_unit_id"), {}), member, members)})
        for gene_id in members:
            item = records.setdefault(gene_id, {
                "family_ids": set(),
                "orthogroup_ids": set(),
                "subfamilies": set(),
                "semantic_unit_ids": set(),
                "claim_ids": set(),
                "evidence_scopes": set(),
                "reference_genes": set(),
                "source_ids": set(),
            })
            for value, key in (
                (row.get("family_id", ""), "family_ids"),
                (row.get("orthogroup_id", ""), "orthogroup_ids"),
                (row.get("subfamily_or_clade", ""), "subfamilies"),
                (row.get("semantic_unit_id", ""), "semantic_unit_ids"),
                (row.get("claim_id", ""), "claim_ids"),
                (row.get("evidence_scope", ""), "evidence_scopes"),
                (row.get("reference_gene", ""), "reference_genes"),
                (row.get("source_id", ""), "source_ids"),
            ):
                if str(value).strip():
                    item[key].add(str(value).strip())

    rows: list[dict[str, object]] = []
    for gene_id in sorted(records):
        item = records[gene_id]
        rows.append({
            "candidate_gene_id": gene_id,
            "family_ids": _join(item["family_ids"]),
            "orthogroup_ids": _join(item["orthogroup_ids"]),
            "subfamilies": _join(item["subfamilies"]),
            "semantic_unit_ids": _join(item["semantic_unit_ids"]),
            "claim_ids": _join(item["claim_ids"]),
            "evidence_scopes": _join(item["evidence_scopes"]),
            "reference_genes": _join(item["reference_genes"]),
            "source_ids": _join(item["source_ids"]),
            "route": "b",
            "support_status": "verified_semantic_support",
            "rank_contribution": "none",
        })

    rows = annotate_relevance(rows, supported_relevance_claims)
    frame = pd.DataFrame(rows, columns=CANDIDATE_COLUMNS)
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    candidate_path = out / "route_b_verified_candidates.tsv"
    frame.to_csv(candidate_path, sep="\t", index=False, lineterminator="\n")
    artifact_sha256 = hashlib.sha256(candidate_path.read_bytes()).hexdigest()
    # Keep grounded model-selected members whose verification is incomplete in
    # a separate sealed review lane. The transfer manifest's frozen bundle is
    # authoritative for both the member roster and identity anchors.
    review_statuses = {"proposed", "not_evaluated", "provider_failed", "unsupported"}
    bundle_path = None
    review_policy = "legacy_grounded_review"
    eligibility_by_claim: dict[str, dict[str, str]] = {}
    bundle_units: dict[str, dict[str, object]] = {}
    review_taxonomy = {}
    if claims_manifest is not None:
        entry = (claims_manifest.get("ancestors") or {}).get("bundle")
        if isinstance(entry, dict):
            bundle_path = runtime_path(str(entry.get("path", "")))
            expected = str(entry.get("sha256", "")).removeprefix("sha256:")
            if not bundle_path.exists() or hashlib.sha256(bundle_path.read_bytes()).hexdigest() != expected:
                raise ValueError("Route B transfer-verification bundle hash mismatch")
            bundle = json.loads(bundle_path.read_text(encoding="utf-8"))
            if not isinstance(bundle, dict) or bundle.get("route") != "interval_first":
                raise ValueError("Route B transfer-verification bundle is invalid")
            units = bundle.get("units")
            if not isinstance(units, list):
                raise ValueError("Route B frozen bundle lacks units")
            bundle_units = {str(unit.get("unit_id", "")): unit for unit in units if isinstance(unit, dict)}
            review_policy = str(bundle.get("b_review_membership_policy", "legacy_grounded_review"))
            if review_policy not in {"legacy_grounded_review", "source_function_v1", "initial_bibliography_v1"}:
                raise ValueError("unsupported Route B review membership policy")
            if review_policy == "initial_bibliography_v1":
                from .reference_taxonomy import load_bound_taxonomy_artifacts
                review_taxonomy = load_bound_taxonomy_artifacts(
                    claims_manifest.get("reference_species_display_audit", {}).get("taxonomy_cache_artifacts", []))
            if review_policy in {"source_function_v1", "initial_bibliography_v1"}:
                if claims_manifest.get("b_review_membership_policy") != review_policy:
                    raise ValueError("Route B eligibility policy differs from frozen bundle")
                outputs = claims_manifest.get("outputs", {})
                transfer_path = runtime_path(str(outputs.get("transfer_verification", "")))
                if (claims_manifest.get("artifact_sha256") != hashlib.sha256(claims_path.read_bytes()).hexdigest()
                        or not transfer_path.is_file()
                        or outputs.get("transfer_verification_sha256") != hashlib.sha256(transfer_path.read_bytes()).hexdigest()):
                    raise ValueError("Route B reference-function eligibility input hash mismatch")
                assessed = pd.read_csv(transfer_path, sep="\t", dtype=str).fillna("")
                if not {"claim_id", "review_eligibility", "review_member_genes"} <= set(assessed.columns):
                    raise ValueError("Route B eligibility table lacks assessment columns")
                if assessed["claim_id"].duplicated().any() or set(assessed["claim_id"]) != set(claims["claim_id"]):
                    raise ValueError("Route B eligibility claim roster mismatch")
                eligibility_by_claim = {row["claim_id"]: row for row in assessed.to_dict("records")}
    review_rows: list[dict[str, str]] = []
    for row in claims[claims["verification_status"].isin(review_statuses)].to_dict("records"):
        unit_id = str(row.get("semantic_unit_id", "")).strip()
        unit = bundle_units.get(unit_id)
        members = _split(row.get("member_interval_genes", ""))
        frozen_members = _split(";".join(str(x) for x in (unit or {}).get("member_interval_genes", [])))
        if not unit or not members or not members <= frozen_members:
            continue
        assessment = eligibility_by_claim.get(str(row.get("claim_id", "")), {})
        versioned = review_policy in {"source_function_v1", "initial_bibliography_v1"}
        review_members = _split(assessment.get("review_member_genes", "")) if versioned else members
        expected_eligibility = "eligible_initial_candidate" if review_policy == "initial_bibliography_v1" else "eligible_reference_function_review"
        if versioned and (
            assessment.get("review_eligibility") != expected_eligibility
            or not review_members or not review_members <= members
        ):
            continue
        identity = unit.get("identity", {}) if isinstance(unit.get("identity"), dict) else {}
        families = {str(x).strip() for key in ("evolutionary_family_ids", "family_ids") for x in identity.get(key, []) if str(x).strip()}
        orthogroups = {str(x).strip() for x in identity.get("orthogroup_ids", []) if str(x).strip()}
        subfamilies = {str(x).strip() for x in identity.get("subfamily_ids", []) if str(x).strip()}
        if str(row.get("family_id", "")).strip() and str(row["family_id"]).strip() not in families:
            continue
        if str(row.get("orthogroup_id", "")).strip() and str(row["orthogroup_id"]).strip() not in orthogroups:
            continue
        if str(row.get("subfamily_or_clade", "")).strip() and str(row["subfamily_or_clade"]).strip() not in subfamilies:
            continue
        for gene_id in sorted(review_members):
            review_row = {
                "candidate_admission_policy": review_policy,
                "candidate_admission_state": "adopted_first_stage" if review_policy == "initial_bibliography_v1" else "legacy_review_candidate",
                "candidate_gene_id": gene_id, "claim_id": str(row.get("claim_id", "")),
                "semantic_unit_id": unit_id, "family_id": str(row.get("family_id", "")),
                "orthogroup_id": str(row.get("orthogroup_id", "")),
                "subfamily_or_clade": str(row.get("subfamily_or_clade", "")),
                "evidence_scope": str(row.get("evidence_scope", "")),
                "reference_gene": str(row.get("reference_gene", "")),
                "source_id": str(row.get("source_id", "")),
                "verification_status": str(row.get("verification_status", "")),
                "verification_reason": str(row.get("verification_reason", "")),
                "route": "b", "support_status": "unverified_route_b_review", "rank_contribution": "none",
                "reference_function_question_id": assessment.get("reference_function_question_id", ""),
                "reference_function_review_reason": assessment.get("reference_function_review_reason", ""),
            }
            # Classify this exact claim/member pair before any gene aggregation.
            application = ({**row, "candidate_gene_id": gene_id}
                           if review_policy == "initial_bibliography_v1" else row)
            application = {**application, "homolog_projection_json": route_b_homolog_projection(
                row, unit, gene_id, review_members, taxonomy_records=review_taxonomy)}
            review_rows.extend(annotate_relevance([review_row], [application]))
    review_frame = pd.DataFrame(review_rows, columns=REVIEW_COLUMNS).drop_duplicates(
        subset=["candidate_gene_id", "claim_id"], keep="first"
    )
    review_path = out / "route_b_review_candidates.tsv"
    review_frame.to_csv(review_path, sep="\t", index=False, lineterminator="\n")
    review_sha256 = hashlib.sha256(review_path.read_bytes()).hexdigest()
    ledger_rows = []
    for unit_id, unit in sorted(bundle_units.items()):
        unit_claims = claims[claims["semantic_unit_id"].eq(unit_id)].to_dict("records")
        for gene in sorted(str(x) for x in unit.get("member_interval_genes", [])):
            for claim in unit_claims or [{}]:
                assessment = eligibility_by_claim.get(str(claim.get("claim_id", "")), {})
                unlinked_copy = (bool(claim) and (
                    gene not in _split(claim.get("member_interval_genes", ""))
                    or (assessment.get("review_eligibility") in {"eligible_reference_function_review", "eligible_initial_candidate"}
                        and gene not in _split(assessment.get("review_member_genes", "")))))
                ledger_rows.append({"candidate_gene_id": gene, "semantic_unit_id": unit_id,
                                    "claim_id": claim.get("claim_id", ""), "review_membership_policy": review_policy,
                                    "b_review_eligibility": "reference_correspondence_not_evaluated" if unlinked_copy else assessment.get("review_eligibility", "not_evaluated"),
                                    "b_review_eligibility_reason": "specific reference function is not linked to this frozen copy" if unlinked_copy else assessment.get("review_eligibility_reason", "no reference-function assessment recorded"),
                                    "verification_status": claim.get("verification_status", "not_evaluated"),
                                    "request_state": unit.get("request_state", "not_recorded"),
                                    "biological_exclusion": "not_asserted"})
    ledger_path = out / "route_b_review_eligibility.tsv"
    pd.DataFrame(ledger_rows, columns=["candidate_gene_id", "semantic_unit_id", "claim_id", "review_membership_policy",
                                      "b_review_eligibility", "b_review_eligibility_reason", "verification_status",
                                      "request_state", "biological_exclusion"]).to_csv(ledger_path, sep="\t", index=False, lineterminator="\n")
    review_manifest_path = out / "route_b_review_candidates_manifest.json"
    review_manifest = {
        "schema_version": SCHEMA_VERSION, "contract": "route-review-candidates-1.0.0",
        "relevance_contract": RELEVANCE_CONTRACT,
        "stage": "route_b_review_candidate_freeze", "case_id": case_id, "route": "interval_first",
        "review_membership_policy": review_policy,
        "eligibility_ledger": {"path": str(ledger_path), "sha256": hashlib.sha256(ledger_path.read_bytes()).hexdigest()},
        "candidate_gene_count": int(review_frame["candidate_gene_id"].nunique()) if len(review_frame) else 0,
        "review_row_count": int(len(review_frame)), "artifact_sha256": review_sha256, "output": str(review_path),
        "inputs": {
            "claims": {"path": str(claims_path), "sha256": hashlib.sha256(claims_path.read_bytes()).hexdigest()},
            "claims_manifest": {"path": str(claims_manifest_path), "sha256": hashlib.sha256(claims_manifest_path.read_bytes()).hexdigest()} if claims_manifest is not None else None,
            "bundle": {"path": str(bundle_path), "sha256": hashlib.sha256(bundle_path.read_bytes()).hexdigest()} if bundle_path is not None else None,
        },
        "policy": ["grounded members of the imported frozen Route-B bundle only",
                   "unverified_route_b_review is not support or scientific proof",
                   "verification status and reason are preserved; contradicted claims are excluded",
                   "rank_contribution is none"],
    }
    review_manifest_path.write_text(json.dumps(review_manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "contract": "route-b-verified-candidates-1.0.0",
        "relevance_contract": RELEVANCE_CONTRACT,
        "stage": "route_b_verified_candidate_freeze",
        "case_id": case_id,
        "route": "interval_first",
        "input_claim_count": int(len(claims)),
        "verified_claim_count": int(len(verified)),
        "candidate_gene_count": int(len(frame)),
        "semantic_unit_count": int(verified["semantic_unit_id"].nunique()) if len(verified) else 0,
        "artifact_sha256": artifact_sha256,
        "review_candidates": str(review_path),
        "review_manifest": str(review_manifest_path),
        "review_candidate_gene_count": int(review_frame["candidate_gene_id"].nunique()) if len(review_frame) else 0,
        "inputs": {
            "claims": {"path": str(claims_path), "sha256": hashlib.sha256(claims_path.read_bytes()).hexdigest()},
            "claims_manifest": (
                {"path": str(claims_manifest_path), "sha256": hashlib.sha256(claims_manifest_path.read_bytes()).hexdigest()}
                if claims_manifest is not None else None
            ),
        },
        "rank_contribution": "none",
        "output": str(candidate_path),
        "policy": [
            "this artifact contains Route B only and is intended to freeze before A/B comparison",
            "only final verification_status=verified semantic claims contribute members",
            "all supported members of verified group-level claims are retained",
            "multiple copies are not collapsed and no numerical rank is assigned",
            "route agreement cannot influence this pre-comparison Route B artifact",
        ],
    }
    manifest_path = out / "route_b_verified_candidates_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        "candidates": str(candidate_path),
        "manifest": str(manifest_path),
        "review_candidates": str(review_path),
        "review_manifest": str(review_manifest_path),
        "review_candidate_gene_count": int(review_frame["candidate_gene_id"].nunique()) if len(review_frame) else 0,
        "candidate_gene_count": int(len(frame)),
        "verified_claim_count": int(len(verified)),
        "artifact_sha256": artifact_sha256,
    }
