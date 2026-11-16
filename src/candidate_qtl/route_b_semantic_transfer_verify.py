"""Transfer-scope verification for Route-B semantic claims (#403).

This is the second deterministic verification axis after source verification.
It asks how far a source-backed reference-gene claim may be transferred into the
focal interval:

* family / orthogroup / reference-gene relevance may be verified while keeping
  every supported focal member;
* copy-level (``interval_gene``) transfer requires an explicit one-to-one / RBH-
  like relation to the submitted reference identity and exactly one focal member;
* one-to-many / many-to-many / unresolved relations never force one copy.

The output is a final semantic-claim table suitable for ``semantic-claims-freeze``.
No numerical gene score is produced.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import re

import pandas as pd

from .claim_binding import load_bound_claims, route_b_authorizes
from .agent_claim_review import AgentClaimReview, reference_function_review_authorizes
from .llm_pathway_legacy import _record_source_ids, _load_record_cache
from .core import runtime_path
from .semantic_claims import CLAIM_COLUMNS, REVIEW_ONLY_CLAIM_TYPES, normalize_semantic_claims
from .reference_taxonomy import load_taxonomy_cache, resolve_species_display

SCHEMA_VERSION = "1.3"

_MEMBERSHIP_SOURCE_STATUSES = {"supported_mechanistic", "supported_curated_function"}
TRANSFER_STATUSES = (
    "copy_supported",
    "family_supported_copy_unresolved",
    "identity_not_evaluated",
    "identity_mismatch",
    "source_support_insufficient",
    "source_contradicted",
)

TRANSFER_COLUMNS = [
    "claim_id",
    "semantic_unit_id",
    "evidence_scope",
    "source_relation_status",
    "bound_claim_state",
    "bound_claim_reason",
    "transfer_status",
    "transfer_scope",
    "supported_member_genes",
    "reference_gene",
    "reference_link_basis",
    "reference_link_scope",
    "orthology_relationships",
    "orthology_relationship_edges",
    "final_verification_status",
    "verification_reason",
    "rank_contribution",
    "review_eligibility", "review_eligibility_reason", "review_member_genes",
    "reference_function_question_id", "reference_function_review_reason",
]

_STRONG_COPY_RELATIONSHIPS = {
    "1:1",
    "one:one",
    "one to one",
    "one-to-one",
    "one_to_one",
    "reciprocal best hit",
    "reciprocal_best_hit",
    "rbh",
}


@dataclass(frozen=True)
class _ReferenceLink:
    """How the submitted reference identity is linked to a frozen unit.

    ``direct_gene_orthology`` may support copy scope because each relationship
    endpoint is explicitly the submitted reference gene. A
    ``literature_unit_membership`` link proves only that the reference belongs
    with the semantic unit; it never authorizes reuse of copy relationships
    whose endpoint is another reference gene.
    """

    genes: tuple[dict[str, object], ...]
    basis: str
    max_transfer_scope: str


def _norm_relation(value: object) -> str:
    return re.sub(r"\s+", " ", str(value or "").strip().casefold())


def _tokens(value: object) -> set[str]:
    return {
        token.strip()
        for token in str(value or "").replace(",", ";").split(";")
        if token.strip()
    }


def _reference_identity_tokens(value: object) -> set[str]:
    """Tokenise a possibly joined reference endpoint without broad aliases."""
    return _tokens(value)


def _source_identity_tokens(value: object) -> set[str]:
    return _tokens(value)


def _endpoint_species(gene: dict[str, object]) -> str:
    for key in (
        "orthology_reference_species", "reference_species", "orthology_species",
        "species",
    ):
        value = str(gene.get(key, "")).strip()
        if value:
            return value.casefold()
    return ""


def _endpoint_namespace(gene: dict[str, object], reference_gene: str) -> str:
    for key in (
        "orthology_reference_accession_namespace", "reference_accession_namespace",
        "accession_namespace", "orthology_reference_namespace", "reference_namespace",
    ):
        value = str(gene.get(key, "")).strip()
        if value:
            return value.casefold()
    token = next(iter(sorted(_reference_identity_tokens(reference_gene))), "")
    # Common accession displays use ``namespace|accession|name``.  Do not
    # invent a namespace for plain gene symbols/accessions.
    if "|" in token:
        return token.split("|", 1)[0].strip().casefold()
    return ""


def _reference_endpoint_matches(
    reference_gene: str,
    reference_species: str,
    gene: dict[str, object],
) -> bool:
    endpoint = _reference_identity_tokens(gene.get("orthology_reference_gene", ""))
    requested = _reference_identity_tokens(reference_gene)
    if not endpoint or not requested or not endpoint.intersection(requested):
        return False
    endpoint_species = _endpoint_species(gene)
    if reference_species and endpoint_species and reference_species.casefold().strip() != endpoint_species:
        return False
    requested_namespace = _endpoint_namespace({}, reference_gene)
    endpoint_namespace = _endpoint_namespace(gene, str(gene.get("orthology_reference_gene", "")))
    if requested_namespace and endpoint_namespace and requested_namespace != endpoint_namespace:
        return False
    return True


def _literature_link_matches(
    reference_gene: str,
    reference_species: str,
    source_id: str,
    record: dict[str, object],
    *,
    initial_policy: bool = False,
) -> bool:
    candidate_ref = str(record.get("reference_gene", "")).strip()
    candidate_source = str(record.get("source_id", "")).strip()
    if not (_reference_identity_tokens(reference_gene) & _reference_identity_tokens(candidate_ref)):
        return False
    if not source_id or not (_source_identity_tokens(source_id) & _source_identity_tokens(candidate_source)):
        return False
    requested_species = reference_species.casefold().strip()
    candidate_species = _endpoint_species(record)
    candidate_namespace = _endpoint_namespace(record, candidate_ref)
    requested_namespace = _endpoint_namespace({}, reference_gene)
    if initial_policy:
        requested_species = _initial_reference_metadata(requested_species)
        candidate_species = _initial_reference_metadata(candidate_species)
        requested_namespace = _initial_reference_metadata(requested_namespace)
        candidate_namespace = _initial_reference_metadata(candidate_namespace)
    if requested_species and candidate_species and requested_species != candidate_species:
        return False
    return not (requested_namespace and candidate_namespace and requested_namespace != candidate_namespace)


def _unit_maps(bundle: dict[str, object]) -> dict[str, dict[str, object]]:
    units_raw = bundle.get("units", [])
    if not isinstance(units_raw, list):
        raise ValueError("Route B semantic bundle requires units list")
    units = {str(unit.get("unit_id", "")): unit for unit in units_raw if isinstance(unit, dict)}
    if (units and "" in units) or len(units) != len(units_raw):
        raise ValueError("Route B semantic bundle requires non-empty unique unit_id values")
    return units


def _identity_check(claim: pd.Series, unit: dict[str, object]) -> tuple[bool, str, list[dict[str, object]]]:
    members = _tokens(claim.get("member_interval_genes", ""))
    unit_members = {str(value).strip() for value in unit.get("member_interval_genes", []) if str(value).strip()}
    if not members or not members <= unit_members:
        return False, "claim member set is not a non-empty subset of frozen semantic unit", []

    identity = unit.get("identity", {}) if isinstance(unit.get("identity"), dict) else {}
    family_ids = {
        str(value).strip()
        for key in ("evolutionary_family_ids", "family_ids")
        for value in identity.get(key, [])
        if str(value).strip()
    }
    orthogroups = {str(value).strip() for value in identity.get("orthogroup_ids", []) if str(value).strip()}
    subfamilies = {str(value).strip() for value in identity.get("subfamily_ids", []) if str(value).strip()}
    if str(claim.get("family_id", "")).strip() and str(claim["family_id"]).strip() not in family_ids:
        return False, "claim family_id differs from frozen semantic unit", []
    if str(claim.get("orthogroup_id", "")).strip() and str(claim["orthogroup_id"]).strip() not in orthogroups:
        return False, "claim orthogroup_id differs from frozen semantic unit", []
    if str(claim.get("subfamily_or_clade", "")).strip() and str(claim["subfamily_or_clade"]).strip() not in subfamilies:
        return False, "claim subfamily differs from frozen semantic unit", []

    # Claim applicability cannot erase frozen competing copy evidence.
    genes = [record for record in unit.get("genes", []) if isinstance(record, dict)]
    return True, "identity anchors match frozen semantic unit", genes


def _reference_link(
    reference_gene: str,
    reference_species: str,
    source_id: str,
    unit: dict[str, object],
    genes: list[dict[str, object]],
) -> _ReferenceLink:
    """Return the linked genes together with the strongest justified scope."""
    if not reference_gene:
        return _ReferenceLink((), "none", "none")

    direct = tuple(
        gene
        for gene in genes
        if _reference_endpoint_matches(reference_gene, reference_species, gene)
    )
    if direct:
        return _ReferenceLink(direct, "direct_gene_orthology", "copy")

    # reverse_literature rows are linked to this semantic unit through
    # interval_gene_ids before bundling. That establishes group membership
    # only: relationships on the unit genes may point to entirely different
    # reference identities and therefore cannot be borrowed for copy transfer.
    literature_refs = [
        record
        for record in unit.get("literature_link_candidates", [])
        if isinstance(record, dict)
        and str(record.get("reference_gene", "")).strip()
        and str(record.get("source_id", "")).strip()
    ]
    if any(_literature_link_matches(reference_gene, reference_species, source_id, record) for record in literature_refs):
        return _ReferenceLink(tuple(genes), "literature_unit_membership", "group")
    return _ReferenceLink((), "none", "none")


def _copy_relationships(linked_genes: tuple[dict[str, object], ...]) -> list[str]:
    return sorted({
        str(gene.get("orthology_relationship", "")).strip()
        for gene in linked_genes
        if str(gene.get("orthology_relationship", "")).strip()
    })


def _initial_reference_metadata(value: object) -> str:
    text = str(value or "").strip().casefold()
    return "" if text in {"", "unknown", "n/a", "na", "none", "null", "not_evaluated", "unresolved", "-",
                           "mixed_or_unspecified", "mixed_or_unspecified_protein_database"} else text


def _initial_reference_endpoints(gene: dict[str, object], cards: list[dict[str, object]], *,
                                 expected_species: str = "", taxonomy_records: dict | None = None,
                                 display_audit: list | None = None) -> dict[str, list[tuple[str, str]]]:
    """Keep species/namespace on the exact endpoint that supplied them."""
    endpoints: dict[str, list[tuple[str, str]]] = {}

    def add(references: object, species: object = "", namespace: object = "", taxon_id: object = "") -> None:
        for reference in _tokens(references):
            species_value = _initial_reference_metadata(species)
            if taxonomy_records:
                resolved, audit = resolve_species_display(species, expected_species, taxon_id, taxonomy_records)
                if audit is not None:
                    species_value = resolved
                    if display_audit is not None:
                        display_audit.append({"reference_gene": reference, "interval_gene_id": str(gene.get("gene_id", "")), **audit})
            endpoints.setdefault(reference, []).append((
                species_value,
                _initial_reference_metadata(namespace),
            ))

    orthology = str(gene.get("orthology_reference_gene", "")).strip()
    # A lookup for another reference cannot lend its species to this orthology
    # endpoint. Existing explicit orthology metadata remains authoritative.
    same_reference = str(gene.get("reference_gene", "")).strip() in {"", orthology}
    add(orthology, gene.get("orthology_reference_species", "") or gene.get("orthology_species", "")
        or (gene.get("reference_species", "") if same_reference else ""),
        gene.get("orthology_reference_accession_namespace", "") or gene.get("orthology_reference_namespace", "")
        or (gene.get("reference_accession_namespace", "") or gene.get("accession_namespace", "") if same_reference else ""),
        gene.get("orthology_reference_taxon_id", "") or (gene.get("reference_taxon_id", "") if same_reference else ""))
    add(gene.get("reference_gene", ""), gene.get("reference_species", ""),
        gene.get("reference_accession_namespace", "") or gene.get("accession_namespace", ""), gene.get("reference_taxon_id", ""))
    add(gene.get("family_index_subject_endpoints", ""))
    for hit in gene.get("sequence_similarity_hits", []):
        if isinstance(hit, dict):
            add(hit.get("hit_accession", ""), hit.get("hit_species", ""), hit.get("reference_accession_namespace", ""), hit.get("hit_taxon_id", ""))
    for card in cards:
        if str(card.get("interval_gene_id", "")).strip() == str(gene.get("gene_id", "")).strip():
            add(card.get("reference_gene", ""), card.get("reference_species", ""), card.get("reference_accession_namespace", ""), card.get("taxon_id", ""))
    return endpoints


def _initial_reference_conflicts(references: set[str], species: str, endpoints: dict[str, list[tuple[str, str]]]) -> bool:
    """An unspecified group link cannot override explicit endpoint conflicts.

    Alternative records may give a valid correspondence for this same copy.
    Require one such record to resolve every conflicting axis together; never
    combine species from one endpoint record with namespace from another.
    """
    conflicted = False
    for reference in references & endpoints.keys():
        expected = (species, _initial_reference_metadata(_endpoint_namespace({}, reference)))
        records = endpoints[reference]
        axes = {i for record in records for i in range(2) if expected[i] and record[i] and expected[i] != record[i]}
        if not axes or any(
            all(record[i] == expected[i] for i in axes)
            and all(not expected[i] or not record[i] or record[i] == expected[i] for i in range(2))
            for record in records
        ):
            return False
        conflicted = True
    return conflicted


def _initial_reference_members(claim: pd.Series, unit: dict, link: _ReferenceLink, *,
                               taxonomy_records: dict | None = None, display_audit: list | None = None) -> set[str]:
    """Retain initial-policy copies without overriding known endpoint identity."""
    references = _tokens(claim.get("reference_gene", ""))
    species = _initial_reference_metadata(claim.get("reference_species", ""))
    members = {str(gene.get("gene_id", "")).strip() for gene in link.genes if str(gene.get("gene_id", "")).strip()}
    # Unknown literature metadata is missing identity information, not a known
    # mismatch. Restore only the initial group mapping; strict transfer keeps
    # its original link and the frozen source records retain their raw values.
    if not members and any(
        _literature_link_matches(str(claim.get("reference_gene", "")), str(claim.get("reference_species", "")),
                                 str(claim.get("source_id", "")), record, initial_policy=True)
        for record in unit.get("literature_link_candidates", []) if isinstance(record, dict)
    ):
        members.update(str(gene.get("gene_id", "")).strip() for gene in unit.get("genes", []) if isinstance(gene, dict))
    cards_by_gene: dict[str, list[dict[str, object]]] = {}
    for card in unit.get("reference_material_cards", []):
        if isinstance(card, dict):
            cards_by_gene.setdefault(str(card.get("interval_gene_id", "")).strip(), []).append(card)
    for gene in unit.get("genes", []):
        if not isinstance(gene, dict):
            continue
        gene_id = str(gene.get("gene_id", "")).strip()
        endpoints = _initial_reference_endpoints(gene, cards_by_gene.get(gene_id, []),
            expected_species=str(claim.get("reference_species", "")), taxonomy_records=taxonomy_records,
            display_audit=display_audit)
        if references & endpoints.keys():
            members.add(gene_id)
        if _initial_reference_conflicts(references, species, endpoints):
            members.discard(gene_id)
    return members & _tokens(claim.get("member_interval_genes", ""))


def _relationship_edges(linked_genes: tuple[dict[str, object], ...]) -> list[dict[str, str]]:
    """Keep relationship type and both identity endpoints in one record."""
    edges = []
    for gene in linked_genes:
        gene_id = str(gene.get("gene_id", "")).strip()
        reference_gene = str(gene.get("orthology_reference_gene", "")).strip()
        relationship = str(gene.get("orthology_relationship", "")).strip()
        if not any((gene_id, reference_gene, relationship)):
            continue
        edges.append({
            "focal_gene_id": gene_id,
            "reference_gene": reference_gene,
            "relationship": relationship,
        })
    return sorted(
        edges,
        key=lambda edge: (
            edge["focal_gene_id"],
            edge["reference_gene"],
            edge["relationship"],
        ),
    )


def verify_route_b_semantic_transfer(
    proposals_arg: str,
    source_verification_arg: str,
    bundle_arg: str,
    out_dir_arg: str,
    *,
    bound_claims_arg: str | None = None,
    agent_claim_review_arg: str | None = None,
    source_review_records_arg: str | None = None,
    record_cache_arg: str | None = None,
    case_id: str = "",
) -> dict[str, object]:
    """Combine source-level support with frozen identity/transfer scope."""
    proposals_path = runtime_path(proposals_arg)
    source_path = runtime_path(source_verification_arg)
    bundle_path = runtime_path(bundle_arg)
    for label, path in (("proposals", proposals_path), ("source verification", source_path), ("bundle", bundle_path)):
        if not path.exists():
            raise ValueError(f"Route B semantic {label} input does not exist: {path}")

    proposals = normalize_semantic_claims(pd.read_csv(proposals_path, sep="\t", dtype=str).fillna(""))
    source = pd.read_csv(source_path, sep="\t", dtype=str).fillna("")
    required_source = {"claim_id", "source_relation_status"}
    missing_source = sorted(required_source - set(source.columns))
    if missing_source:
        raise ValueError(f"source verification table is missing columns: {', '.join(missing_source)}")
    if source["claim_id"].duplicated().any():
        raise ValueError("source verification claim_id values must be unique")
    if set(source["claim_id"]) != set(proposals["claim_id"]):
        missing = sorted(set(proposals["claim_id"]) - set(source["claim_id"]))
        extra = sorted(set(source["claim_id"]) - set(proposals["claim_id"]))
        raise ValueError(f"proposal/source verification roster mismatch: missing={missing}, extra={extra}")
    source_by_claim = source.set_index("claim_id")
    if bound_claims_arg and agent_claim_review_arg:
        raise ValueError("Route B accepts either --bound-claims or --agent-claim-review, not both")
    bound_claims = load_bound_claims(bound_claims_arg, "b")
    agent_review = None
    records: dict[str, dict[str, object]] = {}
    if agent_claim_review_arg:
        if not source_review_records_arg or not record_cache_arg:
            raise ValueError("Route B independent review requires --source-review-records and --record-cache")
        records = _load_record_cache(record_cache_arg)
        agent_review = AgentClaimReview(agent_claim_review_arg, proposals_arg,
                                        source_review_records_arg, route="b")

    bundle = json.loads(bundle_path.read_text(encoding="utf-8"))
    if not isinstance(bundle, dict) or bundle.get("route") != "interval_first":
        raise ValueError("invalid Route B semantic bundle")
    units = _unit_maps(bundle)
    if bundle.get("b_review_membership_policy") == "initial_bibliography_v1":
        manifest_path = source_path.with_name("semantic_claim_source_verification_manifest.json")
        if not manifest_path.is_file():
            raise ValueError("initial bibliography admission requires the source-verification manifest")
        sealed_source = json.loads(manifest_path.read_text(encoding="utf-8"))
        if (sealed_source.get("verification_policy") != "initial_bibliography_v1"
                or sealed_source.get("proposals_sha256") != hashlib.sha256(proposals_path.read_bytes()).hexdigest()
                or sealed_source.get("bundle_sha256") != hashlib.sha256(bundle_path.read_bytes()).hexdigest()
                or sealed_source.get("artifact_sha256") != hashlib.sha256(source_path.read_bytes()).hexdigest()):
            raise ValueError("initial bibliography source/proposal/bundle hash mismatch")

    taxonomy_records, taxonomy_artifacts = (load_taxonomy_cache(record_cache_arg)
        if bundle.get("b_review_membership_policy") == "initial_bibliography_v1" else ({}, []))
    species_display_audit = []

    final_rows: list[dict[str, str]] = []
    transfer_rows: list[dict[str, object]] = []
    for _, claim in proposals.iterrows():
        claim_id = str(claim["claim_id"])
        unit_id = str(claim["semantic_unit_id"])
        source_row = source_by_claim.loc[claim_id]
        source_status = str(source_row.get("source_relation_status", "")).strip()
        reference_gene = str(claim.get("reference_gene", "")).strip()
        reference_species = str(claim.get("reference_species", "")).strip()
        source_id = str(claim.get("source_id", "")).strip()
        if agent_review is not None:
            accession = reference_gene
            review_claim = claim.to_dict()
            # Some B proposals omit species even though the frozen unit carries
            # an authoritative reference endpoint. Supply only that exact
            # endpoint species; never infer a genus or alias.
            record = records.get(accession)
            if not str(review_claim.get("reference_species", "")).strip():
                unit_for_review = units.get(unit_id, {})
                endpoints = [gene for gene in unit_for_review.get("genes", [])
                             if isinstance(gene, dict) and accession in (
                                 _tokens(gene.get("orthology_reference_gene", ""))
                                 | _tokens(gene.get("reference_gene", "")))]
                species = {str(gene.get(key, "")).strip() for gene in endpoints
                           for key in ("orthology_reference_species", "reference_species", "orthology_species")
                           if str(gene.get(key, "")).strip()}
                if len(species) == 1:
                    review_claim["reference_species"] = next(iter(species))
                elif not species and record:
                    organism = record.get("organism", {})
                    taxon = str(organism.get("taxonId", "")).strip() if isinstance(organism, dict) else ""
                    taxa = {str(gene.get("reference_taxon_id", "")).strip() for gene in endpoints
                            if str(gene.get("reference_taxon_id", "")).strip()}
                    if taxon and taxa == {taxon}:
                        review_claim["reference_species"] = str(organism.get("scientificName", "")).strip()
                        review_claim["species_resolution_basis"] = "frozen_reference_endpoint_and_exact_provider_taxon"
            bound_claims[claim_id] = agent_review.bind(
                review_claim, accession, record, _record_source_ids(record) if record else [],
            )
            agent_review.audit[-1]["reference_species_resolution"] = {
                "input": reference_species, "resolved": review_claim.get("reference_species", ""),
                "basis": review_claim.get("species_resolution_basis", "declared_or_frozen_endpoint_species"),
            }
        bound_ok, bound_state, bound_reason = route_b_authorizes(
            bound_claims.get(claim_id),
            claim_id=claim_id,
            reference_gene=reference_gene,
            source_id=source_id,
        )
        reference_link = _ReferenceLink((), "none", "none")
        relationships: list[str] = []
        relationship_edges: list[dict[str, str]] = []

        if unit_id not in units:
            transfer_status = "identity_mismatch"
            transfer_scope = "none"
            supported_members: list[str] = []
            final_status = "unsupported"
            reason = "semantic_unit_id is absent from frozen Route B bundle"
        else:
            unit = units[unit_id]
            identity_ok, identity_reason, genes = _identity_check(claim, unit)
            supported_members = sorted(_tokens(claim.get("member_interval_genes", "")))
            reference_link = _reference_link(reference_gene, reference_species, source_id, unit, genes)
            relationships = _copy_relationships(reference_link.genes)
            relationship_edges = _relationship_edges(reference_link.genes)
            scope = str(claim.get("evidence_scope", "unknown")).strip()

            if not identity_ok:
                transfer_status = "identity_mismatch"
                transfer_scope = "none"
                final_status = "unsupported"
                reason = identity_reason
            elif source_status == "contradicted":
                transfer_status = "source_contradicted"
                transfer_scope = "none"
                final_status = "contradicted"
                reason = "source verification determined that evidence contradicts the proposed mechanism"
            elif source_status == "unsupported":
                transfer_status = "source_support_insufficient"
                transfer_scope = "none"
                final_status = "unsupported"
                reason = "source verification marked the claim unsupported"
            elif source_status not in _MEMBERSHIP_SOURCE_STATUSES and not (
                bound_ok and source_status == "review_required"
            ):
                transfer_status = "source_support_insufficient"
                transfer_scope = "none"
                final_status = "not_evaluated"
                reason = (
                    f"source relation is {source_status or 'not_evaluated'}, not source-level mechanistic "
                    "or curated functional support"
                )
            elif not bound_ok:
                transfer_status = (
                    "source_contradicted"
                    if bound_state == "contradicted"
                    else "source_support_insufficient"
                )
                transfer_scope = "none"
                final_status = "contradicted" if bound_state == "contradicted" else "not_evaluated"
                reason = (
                    "source triage cannot authorize Route-B membership without a complete "
                    f"gene/experiment/phenotype bound claim: {bound_reason}"
                )
            elif not reference_gene or not reference_link.genes:
                transfer_status = "identity_not_evaluated"
                transfer_scope = "none"
                final_status = "not_evaluated"
                reason = f"{source_status} exists, but no frozen reference-gene-to-unit identity link was available"
            elif scope == "interval_gene" and reference_link.max_transfer_scope != "copy":
                transfer_status = "identity_not_evaluated"
                transfer_scope = "family_or_orthogroup"
                final_status = "not_evaluated"
                reason = (
                    "the submitted reference is linked to the frozen semantic unit only at group scope; "
                    "copy relationships whose endpoint is another reference identity cannot verify an interval_gene claim"
                )
            elif scope == "interval_gene":
                claim_members = sorted(_tokens(claim.get("member_interval_genes", "")))
                strong = [
                    gene
                    for gene in reference_link.genes
                    if _reference_endpoint_matches(reference_gene, reference_species, gene)
                    and _norm_relation(gene.get("orthology_relationship", "")) in _STRONG_COPY_RELATIONSHIPS
                ]
                strong_ids = sorted({
                    str(gene.get("gene_id", "")).strip()
                    for gene in strong
                    if str(gene.get("gene_id", "")).strip()
                })
                if len(claim_members) == 1 and strong_ids == claim_members:
                    transfer_status = "copy_supported"
                    transfer_scope = "interval_gene"
                    supported_members = strong_ids
                    final_status = "verified"
                    reason = (
                        f"{source_status} plus explicit one-to-one/RBH-like focal copy relation "
                        "whose endpoint is the submitted reference gene"
                    )
                else:
                    transfer_status = "family_supported_copy_unresolved"
                    transfer_scope = "family_or_orthogroup"
                    final_status = "not_evaluated"
                    reason = "family/reference relation is supported but copy-level claim is unresolved; retain all supported members"
            elif scope in {"family", "subfamily", "orthogroup", "reference_gene", "process"}:
                if scope == "family" and not str(claim.get("family_id", "")).strip():
                    transfer_status = "identity_not_evaluated"
                    transfer_scope = "none"
                    final_status = "not_evaluated"
                    reason = "family-scoped claim lacks family identity"
                elif scope == "subfamily" and not str(claim.get("subfamily_or_clade", "")).strip():
                    transfer_status = "identity_not_evaluated"
                    transfer_scope = "none"
                    final_status = "not_evaluated"
                    reason = "subfamily-scoped claim lacks subfamily identity"
                elif scope == "orthogroup" and not str(claim.get("orthogroup_id", "")).strip():
                    transfer_status = "identity_not_evaluated"
                    transfer_scope = "none"
                    final_status = "not_evaluated"
                    reason = "orthogroup-scoped claim lacks orthogroup identity"
                else:
                    transfer_status = "family_supported_copy_unresolved"
                    transfer_scope = scope
                    final_status = "verified"
                    reason = (
                        f"{source_status} plus frozen reference-gene/group identity; "
                        "claim verified only at group/reference scope and all focal copies remain"
                    )
            else:
                transfer_status = "identity_not_evaluated"
                transfer_scope = "none"
                final_status = "not_evaluated"
                reason = "unknown evidence scope cannot be transferred automatically"

        if final_status == "verified" and str(claim.get("claim_type", "")) in REVIEW_ONLY_CLAIM_TYPES:
            # A bound source may support a seed's reference component without
            # turning the complete mechanistic hypothesis into verified support.
            # Keep the original claim and independent-review audit for review.
            final_status = "not_evaluated"
            transfer_status = "source_support_insufficient"
            transfer_scope = "none"
            supported_members = []
            reason = (
                f"review_only_claim_type:{claim['claim_type']}; semantic support promotion prohibited; "
                f"source/reference evidence retained for review; {reason}"
            )

        final = {column: str(claim.get(column, "")) for column in CLAIM_COLUMNS}
        final["verification_status"] = final_status
        final["verification_reason"] = reason
        final["rank_contribution"] = "none"
        if bundle.get("b_review_membership_policy") == "initial_bibliography_v1":
            final["candidate_admission_policy"] = "initial_bibliography_v1"
            final["focal_species"] = str(bundle.get("focal_species", ""))
        final_rows.append(final)
        # #901 checks reference function and conditional inference independently
        # of focal trait-support promotion. Neither a citation nor triage suffices.
        function_ok = reference_function_review_authorizes(bound_claims.get(claim_id))
        # The frozen unit/claim roster is the interval anchor. Reference-source
        # identity must be checked, but first-stage Type3 family/paralog reasoning
        # does not require exact focal-copy endpoints or one-to-one orthology.
        identity_ok = unit_id in units and _identity_check(claim, units[unit_id])[0]
        review_members = _tokens(claim.get("member_interval_genes", "")) if identity_ok else set()
        review_ok = (identity_ok and (function_ok or bound_ok)
                     and final_status != "contradicted" and source_status != "contradicted"
                     and (source_status != "unsupported" or (function_ok and claim.get("verification_status") == "proposed")))
        review_eligibility = "eligible_reference_function_review" if review_ok else "not_evaluated"
        review_reason = "checked reference-source function, frozen correspondence and independent semantic bridge" if review_ok else bound_reason
        if not source_id:
            review_eligibility, review_reason = "source_absent", "no external reference-function source supplied; ledger only"
        elif not review_ok and source_status in {"unsupported", "contradicted"}:
            review_eligibility, review_reason = "source_not_supporting", "source/import did not support this claim; not gene-level biological exclusion"
        elif not identity_ok:
            review_reason = "claim identity/member anchors differ from frozen unit"
        elif not review_ok and bound_state == "review_required":
            review_eligibility = "review_not_authorized"
        from .bibliography import INITIAL_POLICY, RESOLVED_STATES, concrete_bridge
        if bundle.get("b_review_membership_policy") == INITIAL_POLICY:
            # Original proposal owns the biological interpretation. Metadata is
            # neither an independent function attestation nor verified support.
            bibliography_ok = (str(source_row.get("bibliography_state", "")) in RESOLVED_STATES
                               or source_status == "supported_curated_function")
            claim_display_audit = []
            review_members = _initial_reference_members(claim, units.get(unit_id, {}), reference_link,
                taxonomy_records=taxonomy_records, display_audit=claim_display_audit)
            species_display_audit.extend({"claim_id": claim_id, **item} for item in claim_display_audit)
            review_ok = bool(identity_ok and review_members and bibliography_ok
                             and concrete_bridge(claim.to_dict())
                             and claim.get("verification_status") == "proposed"
                             and final_status != "contradicted")
            review_eligibility = "eligible_initial_candidate" if review_ok else "not_evaluated"
            review_reason = ("original A/B functional judgment with resolved bibliography and frozen reference correspondence; "
                             "content not independently checked") if review_ok else (
                             "bibliography unresolved, missing specific original bridge, or reference correspondence unresolved; ledger only")
        transfer_rows.append({
            "claim_id": claim_id,
            "semantic_unit_id": unit_id,
            "evidence_scope": str(claim.get("evidence_scope", "")),
            "source_relation_status": source_status,
            "bound_claim_state": bound_state,
            "bound_claim_reason": bound_reason,
            "transfer_status": transfer_status,
            "transfer_scope": transfer_scope,
            "supported_member_genes": ";".join(supported_members),
            "reference_gene": reference_gene,
            "reference_link_basis": reference_link.basis,
            "reference_link_scope": reference_link.max_transfer_scope,
            "orthology_relationships": ";".join(relationships),
            "orthology_relationship_edges": json.dumps(
                relationship_edges,
                ensure_ascii=False,
                sort_keys=True,
                separators=(",", ":"),
            ),
            "final_verification_status": final_status,
            "verification_reason": reason,
            "rank_contribution": "none",
            "review_eligibility": review_eligibility, "review_eligibility_reason": review_reason,
            "review_member_genes": ";".join(sorted(review_members)) if review_ok else "",
            "reference_function_question_id": agent_review.shared_question_by_claim.get(claim_id, "") if agent_review else "",
            "reference_function_review_reason": bound_reason if review_ok else "",
        })

    final_frame = normalize_semantic_claims(pd.DataFrame(final_rows, columns=CLAIM_COLUMNS)) if final_rows \
        else pd.DataFrame(columns=CLAIM_COLUMNS)
    transfer_frame = pd.DataFrame(transfer_rows, columns=TRANSFER_COLUMNS).sort_values("claim_id", kind="stable") \
        if transfer_rows else pd.DataFrame(columns=TRANSFER_COLUMNS)

    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    agent_audit_path = agent_review.write_audit(out) if agent_review is not None else None
    final_path = out / "semantic_claims_verified.tsv"
    transfer_path = out / "semantic_claim_transfer_verification.tsv"
    final_frame.to_csv(final_path, sep="\t", index=False, lineterminator="\n")
    transfer_frame.to_csv(transfer_path, sep="\t", index=False, lineterminator="\n")

    transfer_counts = {status: int(transfer_frame["transfer_status"].eq(status).sum()) for status in TRANSFER_STATUSES} \
        if len(transfer_frame) else {status: 0 for status in TRANSFER_STATUSES}
    final_counts = {
        status: int(final_frame["verification_status"].eq(status).sum())
        for status in ("verified", "unsupported", "not_evaluated", "proposed", "contradicted", "provider_failed")
    } if len(final_frame) else {}
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "contract": "route-b-semantic-transfer-verification-1.0.0",
        "stage": "route_b_semantic_transfer_verification",
        "case_id": case_id or str(bundle.get("case_id", "")),
        "claim_count": int(len(final_frame)),
        "transfer_status_counts": transfer_counts,
        "final_verification_status_counts": final_counts,
        "verified_claim_count": int(final_frame["verification_status"].eq("verified").sum()) if len(final_frame) else 0,
        "narrowing_eligible_claim_count": int(final_frame["verification_status"].eq("verified").sum()) if len(final_frame) else 0,
        "rank_contribution": "none",
        "artifact_sha256": hashlib.sha256(final_path.read_bytes()).hexdigest(),
        "b_review_membership_policy": bundle.get("b_review_membership_policy", "legacy_grounded_review"),
        **({"reference_species_display_audit": {"taxonomy_cache_artifacts": taxonomy_artifacts,
            "resolutions": species_display_audit, "raw_inputs_preserved": True,
            "scope": "initial candidate identity only; strict support transfer unchanged"}} if taxonomy_artifacts else {}),
        "ancestors": {
            "proposal_status": ({"path": str(proposals_path.with_name("semantic_claim_proposal_status.tsv")),
                "sha256": hashlib.sha256(proposals_path.with_name("semantic_claim_proposal_status.tsv").read_bytes()).hexdigest()}
                if bundle.get("b_review_membership_policy") == "initial_bibliography_v1"
                and proposals_path.with_name("semantic_claim_proposal_status.tsv").is_file() else None),
            "source_verification": {"path": str(source_path), "sha256": hashlib.sha256(source_path.read_bytes()).hexdigest()},
            "proposals": {"path": str(proposals_path), "sha256": hashlib.sha256(proposals_path.read_bytes()).hexdigest()},
            "bundle": {"path": str(bundle_path), "sha256": hashlib.sha256(bundle_path.read_bytes()).hexdigest()},
            "bound_claims": (
                {
                    "path": str(runtime_path(bound_claims_arg)),
                    "sha256": hashlib.sha256(runtime_path(bound_claims_arg).read_bytes()).hexdigest(),
                }
                if bound_claims_arg
                else None
            ),
        },
        "outputs": {
            "claims": str(final_path),
            "claims_sha256": hashlib.sha256(final_path.read_bytes()).hexdigest(),
            "transfer_verification": str(transfer_path),
            "transfer_verification_sha256": hashlib.sha256(transfer_path.read_bytes()).hexdigest(),
            "manifest": str(out / "semantic_claim_transfer_verification_manifest.json"),
            "agent_claim_review_audit": str(agent_audit_path) if agent_audit_path else None,
        },
        "agent_claim_review": {
            "packet": str(runtime_path(agent_claim_review_arg)) if agent_claim_review_arg else None,
            "source_review_records": str(runtime_path(source_review_records_arg)) if source_review_records_arg else None,
            "record_cache": str(runtime_path(record_cache_arg)) if record_cache_arg else None,
            "audit": str(agent_audit_path) if agent_audit_path else None,
            "audit_sha256": hashlib.sha256(agent_audit_path.read_bytes()).hexdigest() if agent_audit_path else None,
        } if agent_review is not None else None,
        "policy": [
            "review-only claim types remain unverified review hypotheses even when a bound source supports a reference component",
            "final membership verification requires mechanistic literature or frozen curated functional support "
            "plus a frozen focal identity/transfer link",
            "curated functional support is a weaker membership tier and is never relabeled as mechanistic or causal proof",
            "family/group claims may be verified while every supported focal copy remains in the candidate set",
            "literature-unit membership establishes group scope only and cannot satisfy interval_gene copy identity",
            "copy-level interval_gene claims require exactly one focal member plus an explicit one-to-one/RBH-like "
            "relation whose endpoint is the submitted reference identity",
            "relationship type and reference endpoint remain paired in orthology_relationship_edges",
            "one-to-many/many-to-many/unresolved identity never forces a preferred copy",
            "associative/background/review-required source support is kept explicit but is insufficient for automatic semantic verification",
            "cue/term-only source support cannot authorize membership without a shared complete bound claim",
            "not_evaluated is never biological counterevidence and no numerical gene score is written",
        ],
    }
    manifest_path = out / "semantic_claim_transfer_verification_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "claims": str(final_path),
        "transfer_verification": str(transfer_path),
        "manifest": str(manifest_path),
        "claim_count": int(len(final_frame)),
        "verified_claim_count": manifest["verified_claim_count"],
        "narrowing_eligible_claim_count": manifest["narrowing_eligible_claim_count"],
    }
