"""Source/relation verification adapter for Route-B semantic proposals (#403).

This module deliberately keeps two scientific questions separate:

1. Does the cited PMID/PMCID/DOI resolve, and what does its retrievable text
   support about the reference gene / process?  Existing ``source_verify`` owns
   this axis.
2. Is that supported function transferable from the reference gene to the focal
   interval family/subfamily/copy?  This adapter does *not* answer that question.

Accordingly, even strong mechanistic source support does not upgrade a semantic
proposal to final ``verification_status=verified``.  It produces a separate
source-relation verification artifact for a later transfer-scope gate.
"""

from __future__ import annotations

import hashlib
import json

import pandas as pd

from .core import runtime_path
from .semantic_claims import normalize_semantic_claims
from .source_verify import source_verify
from .source_review_capture import content_sha256

SCHEMA_VERSION = "1.3"

SOURCE_RELATION_STATUSES = (
    "supported_mechanistic",
    "supported_curated_function",
    "supported_associative",
    "background_only",
    "contradicted",
    "unsupported",
    "review_required",
    "not_evaluated",
)

# Cache-investigation state (#403 5.6) is a bookkeeping axis, deliberately kept
# separate from ``source_relation_status``: whether a citation lookup was
# answered from cache, freshly investigated, or could not be attempted at all
# says nothing about the biology and must never be folded into a support/
# evidence judgement.
#   cache_hit        -- resolved from a locally supplied record cache.
#   new_investigation -- resolved (or attempted) through a live provider call
#                        this run.
#   provider_failure  -- no cache entry and no provider was available/allowed;
#                        the claim was not investigated, which is not evidence
#                        of anything about the claim.
#   not_applicable    -- no literature lookup applies to this row at all
#                        (no source_id, a non-literature/database accession,
#                        an unparsable citation, or a claim the frozen-bundle
#                        import gate already rejected before any lookup ran).
CACHE_INVESTIGATION_STATES = ("cache_hit", "new_investigation", "provider_failure", "not_applicable")

OUTPUT_COLUMNS = [
    "gene_identity_checks_json",
    "bibliography_state",
    "claim_id",
    "source_id",
    "citation_verified",
    "citation_verification_status",
    "claim_support_status",
    "evidence_relation",
    "retrieval_scope",
    "source_relation_status",
    "source_relation_reason",
    "resolved_id",
    "title",
    "journal",
    "year",
    "doi",
    "first_author",
    "claim_terms_found",
    "evidence_cues",
    "source_provider",
    "source_url",
    "retrieved_at",
    "cache_investigation_state",
    "final_semantic_verification_status",
    "transfer_verification_required",
    "rank_contribution",
]

_MECHANISTIC = {
    "supported_direct_genetic_test",
    "supported_perturbation_phenotype",
}
_ASSOCIATIVE = {"associative_support_only"}
_BACKGROUND = {"background_or_review_only"}
_UNSUPPORTED = {"unsupported_mention"}
_NOT_EVALUATED = {"not_evaluated"}

_GENERIC_FUNCTION_VALUES = {
    "", "unknown", "unknown function", "uncharacterized", "uncharacterized protein",
    "hypothetical", "hypothetical protein", "predicted protein", "protein",
    "putative protein", "unspecified",
}

_CURATED_MEMBERSHIP_CLAIM_TYPES = {
    "direct_functional_bridge",
    "mechanistic_process_bridge",
    "mutant_phenotype_bridge",
    "mechanistic_review_seed",
}


def _tokens(value: object) -> set[str]:
    return {
        token.strip()
        for token in str(value or "").replace(",", ";").split(";")
        if token.strip()
    }


def _meaningful_function(value: object) -> bool:
    text = " ".join(str(value or "").strip().casefold().split())
    return text not in _GENERIC_FUNCTION_VALUES and len(text) >= 4


def _bundle_units(bundle_arg: str | None) -> dict[str, dict[str, object]]:
    if not bundle_arg:
        return {}
    path = runtime_path(bundle_arg)
    if not path.exists():
        raise ValueError(f"Route B semantic bundle does not exist: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or payload.get("route") != "interval_first":
        raise ValueError("invalid Route B semantic bundle for source verification")
    units = payload.get("units", [])
    if not isinstance(units, list):
        raise ValueError("Route B semantic bundle units must be a list")
    result: dict[str, dict[str, object]] = {}
    for unit in units:
        if not isinstance(unit, dict) or not str(unit.get("unit_id", "")).strip():
            raise ValueError("Route B semantic bundle contains an invalid unit")
        unit_id = str(unit["unit_id"])
        if unit_id in result:
            raise ValueError(f"duplicate Route B semantic unit_id: {unit_id}")
        result[unit_id] = unit
    return result


def _curated_function_relation(
    proposal: dict[str, object],
    unit: dict[str, object] | None,
    source_id: str,
) -> tuple[bool, str, dict[str, str]]:
    """Verify a non-literature functional record frozen into the exact B0 unit.

    Accession existence alone is deliberately insufficient. The proposal's
    source and reference gene must be present in the same frozen unit, the unit
    must retain traceable provider provenance, and the linked gene must carry a
    non-generic function/process/phenotype annotation. This verifies a
    source-backed functional bridge, not causal or mechanistic literature.
    """
    if unit is None:
        return False, "database accession has no frozen Route B semantic unit", {}
    claim_type = str(proposal.get("claim_type", "")).strip()
    if claim_type not in _CURATED_MEMBERSHIP_CLAIM_TYPES:
        return False, f"claim_type {claim_type or '<empty>'} is not an explicit curated membership bridge", {}
    stable_sources = {str(value).strip() for value in unit.get("stable_source_ids", []) if str(value).strip()}
    if source_id not in stable_sources:
        return False, "database accession is absent from frozen unit stable_source_ids", {}

    reference_gene = str(proposal.get("reference_gene", "")).strip()
    if not reference_gene:
        return False, "database functional record has no grounded reference_gene", {}

    matching_provenance: dict[str, object] | None = None
    for record in unit.get("source_provenance", []):
        if not isinstance(record, dict):
            continue
        record_sources = {str(value).strip() for value in record.get("source_ids", []) if str(value).strip()}
        if source_id not in record_sources:
            continue
        if str(record.get("reference_gene", "")).strip() != reference_gene:
            continue
        if not any(str(record.get(field, "")).strip() for field in ("source_database", "provenance_detail", "source_url")):
            continue
        matching_provenance = record
        break
    if matching_provenance is None:
        return False, "database accession lacks matching frozen reference-gene provenance", {}

    functional_fields = (
        "reference_mutant_phenotype", "reference_function", "reference_pathway",
        "protein_family", "go_terms", "keywords",
    )
    matching_gene: dict[str, object] | None = None
    matched_fields: list[str] = []
    for gene in unit.get("genes", []):
        if not isinstance(gene, dict):
            continue
        gene_references = {
            str(gene.get("reference_gene", "")).strip(),
            str(gene.get("orthology_reference_gene", "")).strip(),
        }
        if reference_gene not in gene_references:
            continue
        if source_id not in _tokens(gene.get("reference_source_id", "")):
            continue
        if str(gene.get("lookup_status", "")).strip() != "source_backed":
            continue
        current_fields = [field for field in functional_fields if _meaningful_function(gene.get(field, ""))]
        if not current_fields:
            continue
        matching_gene = gene
        matched_fields = current_fields
        break
    matching_card: dict[str, object] | None = None
    if matching_gene is None:
        # B0 family/hit genes can carry only IDs while the actual reference
        # function is retained in the existing frozen reference-material card.
        # Preserve that content for independent review; ID/provenance alone
        # still cannot establish functional or focal-trait support.
        for card in unit.get("reference_material_cards", []):
            if not isinstance(card, dict):
                continue
            if (str(card.get("source_id", "")).strip() != source_id
                    or str(card.get("reference_gene", "")).strip() != reference_gene):
                continue
            if not any(str(card.get(field, "")).strip() for field in ("source_database", "provenance_detail", "source_url")):
                continue
            fields = {"reference_function": card.get("function_description", ""),
                      "reference_mutant_phenotype": card.get("mutant_phenotype", ""),
                      "reference_pathway": card.get("pathway", "")}
            matched_fields = [field for field, value in fields.items() if _meaningful_function(value)]
            if not matched_fields:
                continue
            matching_card = card
            matching_gene = {"reference_gene": reference_gene, "reference_source_id": source_id,
                             "reference_species": card.get("reference_species", ""), **fields,
                             "source_database": card.get("source_database", ""),
                             "source_url": card.get("source_url", ""), "retrieved_at": card.get("retrieved_at", "")}
            break
    if matching_gene is None:
        return False, "database accession is grounded but has no non-generic frozen functional annotation", {}

    database = str(matching_provenance.get("source_database", "") or matching_gene.get("source_database", "")).strip()
    metadata = {
        "claim_support_status": "source_backed_curated_function",
        "evidence_relation": "curated_function_annotation",
        "retrieval_scope": "frozen_route_b_bundle",
        "resolved_id": source_id,
        "source_provider": f"frozen_bundle:{database}" if database else "frozen_route_b_bundle",
        "source_url": str(matching_provenance.get("source_url", "") or matching_gene.get("source_url", "")).strip(),
        "retrieved_at": str(matching_provenance.get("retrieved_at", "") or matching_gene.get("retrieved_at", "")).strip(),
        # Exact frozen records for source reading; this is not a positive review.
        "curated_record": json.dumps({
            "gene": {key: value for key, value in matching_gene.items() if key.startswith("reference_")
                     or key in {"source_database", "source_url", "retrieved_at", "lookup_status"}},
            "provenance": {key: value for key, value in matching_provenance.items() if key != "gene_id"},
            **({"reference_material_card": matching_card} if matching_card is not None else {}),
        }, ensure_ascii=False, sort_keys=True),
    }
    return True, f"source-backed curated function in frozen unit fields: {','.join(matched_fields)}", metadata


def _bool_text(value: object) -> bool:
    return str(value).strip().casefold() in {"true", "1", "yes"}


def _source_relation(row: pd.Series | None, source_id: str) -> tuple[str, str]:
    if not source_id:
        return "not_evaluated", "semantic proposal supplied no literature source_id"
    if row is None:
        return "not_evaluated", "source verifier produced no citation row"

    verification = str(row.get("verification_status", "")).strip()
    support = str(row.get("claim_support_status", "")).strip()
    review_reason = str(row.get("review_reason", "")).strip()
    rejection_reason = str(row.get("rejection_reason", "")).strip()

    if verification == "rejected":
        if "contradict" in support.casefold() or "contradict" in rejection_reason.casefold() or "counterevidence" in rejection_reason.casefold():
            return "contradicted", rejection_reason or "citation evidence contradicts proposed mechanism"
        return "unsupported", rejection_reason or "citation was rejected by source verification"
    if "contradict" in support.casefold() or "contradict" in review_reason.casefold():
        return "contradicted", review_reason or support
    if verification == "review_required" and review_reason in {
        "verification_not_performed_offline",
        "literature_provider_unavailable",
    }:
        return "not_evaluated", review_reason
    if support in _MECHANISTIC and _bool_text(row.get("citation_verified", "")):
        return "supported_mechanistic", support
    if support in _ASSOCIATIVE and _bool_text(row.get("citation_verified", "")):
        return "supported_associative", support
    if support in _BACKGROUND and _bool_text(row.get("citation_verified", "")):
        return "background_only", support
    if support in _UNSUPPORTED:
        return "unsupported", support
    if support in _NOT_EVALUATED:
        return "not_evaluated", review_reason or support
    return "review_required", review_reason or support or "source relation requires review"


def _cache_investigation_state(row: pd.Series | None, source_id: str, claim_verification_status: str) -> str:
    """Classify how (or whether) this claim's citation lookup actually ran.

    Derived only from fields ``source_verify`` already writes per citation
    (``source_provider`` / ``verification_status`` / ``review_reason``), so no
    behaviour of the shared ``source_verify`` engine changes. This function
    must not be used to adjust ``source_relation_status``: a provider that
    could not be reached is not evidence, so it never becomes ``unsupported``
    -- ``_source_relation`` already resolves that case to ``not_evaluated``.
    """
    if claim_verification_status == "unsupported":
        # Already rejected by the frozen-bundle import gate; no lookup ran.
        return "not_applicable"
    if not source_id or row is None:
        return "not_applicable"

    provider = str(row.get("source_provider", "")).strip()
    verification = str(row.get("verification_status", "")).strip()
    review_reason = str(row.get("review_reason", "")).strip()

    if verification == "not_literature_source":
        # A database accession, not a literature citation: out of this axis.
        return "not_applicable"
    if provider == "local_record_cache":
        return "cache_hit"
    if str(row.get("bibliography_state", "")) == "provider_failure":
        return "provider_failure"
    if verification == "review_required" and review_reason in {
        "verification_not_performed_offline",
        "literature_provider_unavailable",
    }:
        return "provider_failure"
    provider_parts = {part for part in provider.split("+") if part}
    if provider_parts & {"europe_pmc_rest", "europe_pmc_lite", "ncbi_eutils_pubmed", "ncbi_eutils_pmc"}:
        # A live provider call was made this run, whether or not it found a
        # matching record -- "not found" is still a completed investigation.
        return "new_investigation"
    # Citation text itself was unparsable (rejected with no provider ever
    # consulted); this is a citation-format problem, not a provider outage.
    return "not_applicable"


def _adapter_table(proposals: pd.DataFrame) -> pd.DataFrame:
    """Create the minimal source_verify input without changing claim meaning."""
    rows: list[dict[str, str]] = []
    for row in proposals.to_dict("records"):
        process_parts = [
            str(row.get("biological_function", "")).strip(),
            str(row.get("mechanism", "")).strip(),
            str(row.get("phenotype_relation", "")).strip(),
        ]
        process = "; ".join(part for part in process_parts if part)
        family_parts = [
            str(row.get("family_id", "")).strip(),
            str(row.get("subfamily_or_clade", "")).strip(),
            str(row.get("orthogroup_id", "")).strip(),
        ]
        family = "; ".join(part for part in family_parts if part)
        rows.append({
            "claim_id": str(row["claim_id"]),
            "pathway_or_process": process,
            "gene_family": family,
            "reference_gene_name": str(row.get("reference_gene", "")).strip(),
            "source_ids": str(row.get("source_id", "")).strip(),
            "source_title": str(row.get("source_title", "")).strip(),
            "source_identifier_aliases": str(row.get("source_identifier_aliases", "")).strip(),
            "paper_reports_json": str(row.get("paper_reports_json", "")).strip(),
        })
    return pd.DataFrame(
        rows,
        columns=["claim_id", "pathway_or_process", "gene_family", "reference_gene_name", "source_ids", "source_title", "source_identifier_aliases", "paper_reports_json"],
    )


def _initial_bibliography_budget(proposed: pd.DataFrame) -> tuple[pd.DataFrame, set[str], list[dict]]:
    """Consume the imported paper handoff; unmarked inputs keep the old budget.

    Shared scope selection precedes per-paper/copy expansion. The selector reads
    machine-stamped primary papers without choosing another five for metadata.
    """
    from .paper_reports import PAPER_FIELD, PAPER_HANDOFF_POLICY, normalize_paper_reports, select_bibliography_papers, source_tokens
    from .source_verify import classify_citation
    keep, skipped, panels = [], set(), []
    for key, group in proposed.groupby(["semantic_unit_id", "reference_gene", "reference_species"], sort=False):
        papers = [paper for value in group[PAPER_FIELD] for paper in json.loads(normalize_paper_reports(value) or "[]")]
        selection = select_bibliography_papers({PAPER_FIELD: papers if papers else ""}, group.source_id.tolist())
        selected = source_tokens(selection["selected_source_ids"])
        panel = {"semantic_unit_id": key[0], "reference_gene": key[1], "reference_species": key[2], **selection}
        if selection.get("policy") == PAPER_HANDOFF_POLICY:
            panel["audit_only_claim_ids"] = []
        panels.append(panel)
        for index, row in group.iterrows():
            if classify_citation(row.source_id)[0] not in {"doi", "pmid", "pmcid"} or source_tokens([row.source_id]) & selected:
                keep.append(index)
            else:
                skipped.add(str(row.claim_id))
                if "audit_only_claim_ids" in panel:
                    panel["audit_only_claim_ids"].append(str(row.claim_id))
    return proposed.loc[keep].copy(), skipped, panels


def verify_route_b_semantic_sources(
    proposals_arg: str,
    out_dir_arg: str,
    *,
    bundle_arg: str | None = None,
    cache_arg: str | None = None,
    allow_external_api: bool = False,
    full_hash: bool = False,
    case_id: str = "",
) -> dict[str, object]:
    """Verify citation existence/source relation while leaving transfer unresolved."""
    proposals_path = runtime_path(proposals_arg)
    if not proposals_path.exists():
        raise ValueError(f"semantic proposal table does not exist: {proposals_path}")
    proposals = pd.read_csv(proposals_path, sep="\t", dtype=str).fillna("")
    proposals = normalize_semantic_claims(proposals)
    units = _bundle_units(bundle_arg)
    from .bibliography import INITIAL_POLICY
    policy = json.loads(runtime_path(bundle_arg).read_text(encoding="utf-8")).get("b_review_membership_policy") if bundle_arg else None
    initial = policy == INITIAL_POLICY
    if proposals.empty:
        out = runtime_path(out_dir_arg)
        out.mkdir(parents=True, exist_ok=True)
        result_path = out / "semantic_claim_source_verification.tsv"
        pd.DataFrame(columns=OUTPUT_COLUMNS).to_csv(result_path, sep="\t", index=False, lineterminator="\n")
        empty_sources = out / "source_review_records.json"
        empty_sources.write_text(json.dumps({"artifact_role": "source_text_for_review_only", "records": {}}, sort_keys=True) + "\n", encoding="utf-8")
        manifest = {
            "verification_policy": INITIAL_POLICY if initial else "content_v1",
            "proposals_sha256": hashlib.sha256(proposals_path.read_bytes()).hexdigest(),
            "bundle_sha256": hashlib.sha256(runtime_path(bundle_arg).read_bytes()).hexdigest() if bundle_arg else None,
            "artifact_sha256": hashlib.sha256(result_path.read_bytes()).hexdigest(),
            "schema_version": SCHEMA_VERSION,
            "case_id": case_id,
            "claim_count": 0,
            "source_relation_status_counts": {status: 0 for status in SOURCE_RELATION_STATUSES},
            "final_verified_claim_count": 0,
            "transfer_verification_required_count": 0,
            "cache_investigation_state_counts": {state: 0 for state in CACHE_INVESTIGATION_STATES},
            "cache_accounting": {
                "cache_hit_count": 0,
                "new_investigation_count": 0,
                "provider_failure_count": 0,
                "not_applicable_count": 0,
                "provider": "none",
            },
            "rank_contribution": "none",
            "outputs": {"source_verification": str(result_path), "source_review_records": str(empty_sources)},
            "policy": ["no proposals is a valid abstention; no source verification work is fabricated"],
        }
        manifest_path = out / "semantic_claim_source_verification_manifest.json"
        manifest["outputs"]["manifest"] = str(manifest_path)
        manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return {"source_verification": str(result_path), "source_review_records": str(empty_sources), "manifest": str(manifest_path), "claim_count": 0}

    # Only model proposals are eligible for source verification. Rows already
    # rejected as unsupported by the frozen-bundle import gate remain unsupported
    # and do not need an external lookup to become more rejected.
    proposed = proposals[proposals["verification_status"].eq("proposed")].copy()
    other = proposals[~proposals["verification_status"].isin({"proposed", "unsupported"})].copy()
    if not other.empty:
        statuses = ", ".join(sorted(set(other["verification_status"])))
        raise ValueError(
            "route-b semantic source verification expects proposal-import output "
            f"(proposed/unsupported only); found: {statuses}"
        )

    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    validation = pd.DataFrame()
    source_manifest: dict[str, object] = {}
    lookup_proposals, budget_skipped, panels = proposed, set(), []
    selection_path = None
    if initial:
        lookup_proposals, budget_skipped, panels = _initial_bibliography_budget(proposed)
        selection_path = out / "bibliography_selection.json"
        selection_path.write_text(json.dumps({"panels": panels}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    handoff_skipped = {claim for panel in panels for claim in panel.get("audit_only_claim_ids", [])}
    if not proposed.empty:
        adapter_path = out / "semantic_source_verify_input.tsv"
        _adapter_table(lookup_proposals).to_csv(adapter_path, sep="\t", index=False, lineterminator="\n")
        source_out = out / "source_verify"
        source_result = source_verify(
            str(adapter_path),
            str(source_out),
            cache_arg=cache_arg,
            allow_external_api=allow_external_api,
            id_column="claim_id",
            source_column="source_ids",
            full_hash=full_hash,
            verification_policy=INITIAL_POLICY if initial else "content_v1",
        )
        validation_path = source_out / "source_validation.tsv"
        validation = pd.read_csv(validation_path, sep="\t", dtype=str).fillna("")
        source_manifest_path = runtime_path(str(source_result["manifest"]))
        source_manifest = json.loads(source_manifest_path.read_text(encoding="utf-8"))

    validation_by_claim: dict[str, pd.Series] = {}
    source_snapshot = {"artifact_role": "source_text_for_review_only", "records": {}}
    if not proposed.empty:
        source_snapshot = json.loads(runtime_path(str(source_result["source_review_records"])).read_text(encoding="utf-8"))
    if not validation.empty:
        # Proposal import intentionally emits one semantic-claim row per source,
        # so each claim_id maps to at most one literature validation row.
        for claim_id, group in validation.groupby("hypothesis_id", sort=False):
            if len(group) > 1 and not initial:
                raise ValueError(f"semantic claim {claim_id} unexpectedly produced multiple source validation rows")
            chosen = group.iloc[0].copy()
            if initial:
                failed = group[~group["citation_verified"].astype(str).str.casefold().eq("true")]
                if len(failed):
                    chosen = failed.iloc[0].copy()
            validation_by_claim[str(claim_id)] = chosen

    rows: list[dict[str, object]] = []
    for proposal in proposals.to_dict("records"):
        claim_id = str(proposal["claim_id"])
        source_id = str(proposal.get("source_id", "")).strip()
        if str(proposal["verification_status"]) == "unsupported":
            relation_status = "unsupported"
            relation_reason = str(proposal.get("verification_reason", "")) or "proposal import rejected frozen-bundle identity/source"
            validation_row = None
        else:
            validation_row = validation_by_claim.get(claim_id)
            relation_status, relation_reason = _source_relation(validation_row, source_id)
            if claim_id in budget_skipped:
                relation_reason = ("outside_primary_paper_handoff; original paper observation retained"
                                   if claim_id in handoff_skipped else
                                   "metadata_budget_not_selected; original paper observation retained")

        curated_metadata: dict[str, str] = {}
        if (
            relation_status == "not_evaluated"
            and validation_row is not None
            and (str(validation_row.get("verification_status", "")).strip() == "not_literature_source"
                 or (initial and str(validation_row.get("citation_type", "")) == "database_accession"))
        ):
            curated_ok, curated_reason, curated_metadata = _curated_function_relation(
                proposal,
                units.get(str(proposal.get("semantic_unit_id", ""))),
                source_id,
            )
            if curated_ok:
                relation_status = "supported_curated_function"
                relation_reason = curated_reason
                sections = {"curatedRecord": curated_metadata["curated_record"]}
                captured = source_snapshot["records"].setdefault(source_id, {
                    "citation": source_id, "source_record_id": source_id,
                    "source_identifiers": {"accession": str(proposal.get("reference_gene", "")), "source_id": source_id},
                    "sections": sections, "sections_sha256": content_sha256(sections),
                    "hypothesis_ids": [], "semantic_binding_state": "not_evaluated",
                    "available_text_scope": "curated_db_record", "classification_scope": "frozen_route_b_bundle",
                })
                snapshots = captured.setdefault("curated_record_snapshots", [captured["sections"]["curatedRecord"]])
                if sections["curatedRecord"] not in snapshots:
                    snapshots.append(sections["curatedRecord"])
                if claim_id not in captured["hypothesis_ids"]:
                    captured["hypothesis_ids"].append(claim_id)

        def field(name: str) -> str:
            if name in curated_metadata:
                return curated_metadata[name]
            return "" if validation_row is None else str(validation_row.get(name, ""))

        # Crucial boundary: source/relation support alone never verifies transfer
        # to the focal interval member.  Final status remains proposed (or
        # unsupported) until a separate family/subfamily/orthology transfer gate.
        final_semantic_status = "unsupported" if relation_status == "unsupported" else "proposed"
        transfer_required = relation_status in {
            "supported_mechanistic", "supported_curated_function", "supported_associative",
            "background_only", "review_required"
        }
        cache_state = _cache_investigation_state(validation_row, source_id, str(proposal["verification_status"]))
        rows.append({
            "gene_identity_checks_json": field("gene_identity_checks_json"),
            "bibliography_state": ("paper_handoff_audit_only" if claim_id in handoff_skipped else
                                   "metadata_budget_not_selected" if claim_id in budget_skipped else field("bibliography_state")),
            "claim_id": claim_id,
            "source_id": source_id,
            "citation_verified": field("citation_verified"),
            "citation_verification_status": field("verification_status"),
            "claim_support_status": field("claim_support_status"),
            "evidence_relation": field("evidence_relation"),
            "retrieval_scope": field("retrieval_scope"),
            "source_relation_status": relation_status,
            "source_relation_reason": relation_reason,
            "resolved_id": field("resolved_id"),
            "title": field("title"),
            "journal": field("journal"),
            "year": field("year"),
            "doi": field("doi"),
            "first_author": field("first_author"),
            "claim_terms_found": field("claim_terms_found"),
            "evidence_cues": field("evidence_cues"),
            "source_provider": field("source_provider"),
            "source_url": field("source_url"),
            "retrieved_at": field("retrieved_at"),
            "cache_investigation_state": cache_state,
            "final_semantic_verification_status": final_semantic_status,
            "transfer_verification_required": transfer_required,
            "rank_contribution": "none",
        })

    result = pd.DataFrame(rows, columns=OUTPUT_COLUMNS).sort_values("claim_id", kind="stable")
    result_path = out / "semantic_claim_source_verification.tsv"
    result.to_csv(result_path, sep="\t", index=False, lineterminator="\n")
    source_snapshot_path = out / "source_review_records.json"
    for record in source_snapshot["records"].values():
        if "curated_record_snapshots" in record:
            record["sections"] = {"curatedRecord": "\n".join(sorted(record.pop("curated_record_snapshots")))}
            record["sections_sha256"] = content_sha256(record["sections"])
    source_snapshot_path.write_text(json.dumps(source_snapshot, ensure_ascii=False, sort_keys=True, indent=2) + "\n", encoding="utf-8")

    counts = {status: int(result["source_relation_status"].eq(status).sum()) for status in SOURCE_RELATION_STATUSES}
    cache_state_counts = {
        state: int(result["cache_investigation_state"].eq(state).sum()) for state in CACHE_INVESTIGATION_STATES
    }
    if sum(cache_state_counts.values()) != len(result):
        # Fail closed rather than silently drop a row from the accounting: every
        # claim must land in exactly one cache-investigation state.
        raise RuntimeError(
            "cache_investigation_state_counts do not sum to claim_count "
            f"({sum(cache_state_counts.values())} != {len(result)})"
        )
    manifest = {
        "verification_policy": INITIAL_POLICY if initial else "content_v1",
        "proposals_sha256": hashlib.sha256(proposals_path.read_bytes()).hexdigest(),
        "bundle_sha256": hashlib.sha256(runtime_path(bundle_arg).read_bytes()).hexdigest() if bundle_arg else None,
        "artifact_sha256": hashlib.sha256(result_path.read_bytes()).hexdigest(),
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id,
        "claim_count": int(len(result)),
        "source_relation_status_counts": counts,
        "final_verified_claim_count": 0,
        "transfer_verification_required_count": int(result["transfer_verification_required"].astype(bool).sum()),
        "source_verify_provider": str(source_manifest.get("provider", "none")),
        "paper_selection": {"path": str(selection_path), "sha256": hashlib.sha256(selection_path.read_bytes()).hexdigest(),
                            "max_papers_per_unit_reference": 5, "unselected_claim_count": len(budget_skipped),
                            **({"handoff_consumer": "audit_imported_primary_papers_without_reselection",
                                "handoff_audit_only_claim_count": len(handoff_skipped)}
                               if any("audit_only_claim_ids" in panel for panel in panels) else {})} if selection_path else None,
        "cache_investigation_state_counts": cache_state_counts,
        "cache_accounting": {
            "cache_hit_count": cache_state_counts["cache_hit"],
            "new_investigation_count": cache_state_counts["new_investigation"],
            "provider_failure_count": cache_state_counts["provider_failure"],
            "not_applicable_count": cache_state_counts["not_applicable"],
            "provider": str(source_manifest.get("provider", "none")),
        },
        "external_api_allowed": allow_external_api,
        "rank_contribution": "none",
        "outputs": {
            "gene_identity_checks": source_manifest.get("gene_identity_checks"),
            "source_verification": str(result_path),
            "source_review_records": str(source_snapshot_path),
            "source_verify_dir": str(out / "source_verify") if not proposed.empty else "",
            "manifest": str(out / "semantic_claim_source_verification_manifest.json"),
        },
        "policy": [
            "citation resolution and claim support are separate axes inherited from source_verify",
            "fabricated/unrecognised citations are unsupported; provider unavailability is not_evaluated",
            "direct genetic/perturbation text can establish source-level mechanistic support only",
            "a database accession can establish supported_curated_function only when the exact frozen unit also "
            "contains matching reference-gene provenance and non-generic source-backed functional annotation",
            "supported_curated_function is not mechanistic or causal proof and still requires transfer verification",
            "associative/background sources remain weaker explicit states rather than mechanistic support",
            "source-level support never verifies functional transfer to the focal family/subfamily/copy",
            "final semantic verified count is therefore zero until a separate transfer-scope gate passes",
            "nothing in this stage writes a numerical gene score",
            "cache_investigation_state (cache_hit/new_investigation/provider_failure/not_applicable) is a "
            "bookkeeping axis only; provider_failure is never converted into unsupported or treated as evidence "
            "-- the offline-without-provider case already resolves source_relation_status to not_evaluated",
        ],
    }
    manifest_path = out / "semantic_claim_source_verification_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "source_verification": str(result_path),
        "source_review_records": str(source_snapshot_path),
        "manifest": str(manifest_path),
        "claim_count": int(len(result)),
        "source_relation_status_counts": counts,
        "final_verified_claim_count": 0,
        "cache_investigation_state_counts": cache_state_counts,
    }
