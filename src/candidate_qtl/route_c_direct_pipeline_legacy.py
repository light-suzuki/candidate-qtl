"""Executable Route C direct trait-gene rescue stages (issue #640).

The lane is deliberately split into interval-blind discovery/freeze and a
post-freeze projection stage.  Literature retrieval proposes links; only the
generic source-verification gate can make a direct record projection-eligible.
Projection never writes a score and never transfers a cross-species phenotype
to a focal gene by homology alone.
"""

from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from . import __version__
from .core import file_fingerprint, runtime_path, safe_id, sha256
from .reverse_literature import _fetch as europe_pmc_fetch
from .reverse_literature import _load_cache as load_query_cache
from .route_c_direct_rescue import (
    DIRECT_RESCUE_COLUMNS,
    PROJECTION_COLUMNS,
    SPECIES_DISTANCE,
    direct_record_can_seed_projection,
    validate_direct_rescue_frame,
    validate_projection_frame,
)
from .source_verify import source_verify
from .source_review_capture import content_sha256

SCHEMA_VERSION = "1.0"

SPECIES_SCOPE_COLUMNS = (
    "query_id",
    "reference_species",
    "species_distance",
    "phenotype_query",
)

SEARCH_COLUMNS = (
    "query_id",
    "query_phenotype",
    "phenotype_query",
    "focal_species",
    "reference_species",
    "species_distance",
    "query",
    "query_status",
    "article_rank",
    "source_id",
    "pmid",
    "pmcid",
    "doi",
    "title",
    "abstract_text",
    "journal",
    "year",
    "first_author",
    "provider",
    "retrieved_at",
    "rank_contribution",
)

DIRECT_PROPOSAL_COLUMNS = (
    "record_id",
    "query_phenotype",
    "focal_species",
    "focal_cultivar_or_accession",
    "reported_phenotype",
    "reference_species",
    "reference_cultivar_or_accession",
    "reported_gene_id",
    "reported_gene_symbol",
    "reported_locus",
    "reported_variant",
    "species_distance",
    "evidence_directness",
    "evidence_type",
    "source_id",
    "reference_sequence_id",
    "notes",
)

DIRECT_BENCHMARK_CASE_COLUMNS = (
    "case_id",
    "expected_qtl_id",
    "expected_atlas_sha256",
    "direct_rescue_tsv",
    "direct_freeze_manifest",
    "projection_tsv",
    "projection_manifest",
    "truth_reference_gene_id",
    "truth_focal_gene_id",
    "contamination_status",
)

BENCHMARK_CONTAMINATION_STATES = frozenset({"clean", "contaminated", "not_evaluated"})

_DIRECT_SUPPORT = frozenset({
    "supported_direct_genetic_test",
    "supported_perturbation_phenotype",
})
_ASSOCIATIVE_SUPPORT = "associative_support_only"
_DIRECT_CUES = (
    "mutant OR mutation OR allele OR knockout OR knock-out OR loss-of-function OR "
    "gene-editing OR CRISPR OR complementation OR rescue OR transgenic OR fine-mapped"
)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _require_exact(frame: pd.DataFrame, columns: tuple[str, ...], label: str) -> pd.DataFrame:
    missing = sorted(set(columns) - set(frame.columns))
    extra = sorted(set(frame.columns) - set(columns))
    if missing or extra:
        detail = []
        if missing:
            detail.append(f"missing: {', '.join(missing)}")
        if extra:
            detail.append(f"outside the contract: {', '.join(extra)}")
        raise ValueError(f"{label} has {'; '.join(detail)}")
    return frame.fillna("")


def _literature_query(phenotype_query: str, species: str) -> str:
    phenotype = re.sub(r"\s+", " ", str(phenotype_query)).strip()
    species_value = re.sub(r"\s+", " ", str(species)).strip()
    return f'"{phenotype}" AND "{species_value}" AND ({_DIRECT_CUES})'


def direct_literature_search(
    phenotype: str,
    focal_species: str,
    species_scope_arg: str,
    out_dir_arg: str,
    cache_arg: str | None = None,
    allow_external_api: bool = False,
    limit: int = 10,
    case_id: str = "",
    full_hash: bool = False,
) -> dict[str, object]:
    """Search direct trait-gene reports across explicit species/phenotype scopes.

    The scope table is the auditable place to supply close phenotype wording for
    each species.  No QTL, interval atlas, candidate id, or sequence argument is
    accepted by this function.
    """
    if not str(phenotype).strip() or not str(focal_species).strip():
        raise ValueError("phenotype and focal_species are required")
    if limit < 1:
        raise ValueError("limit must be at least 1")
    scope_path = runtime_path(species_scope_arg)
    scope = _require_exact(
        pd.read_csv(scope_path, sep="\t", dtype=str), SPECIES_SCOPE_COLUMNS, "Route C species scope"
    )
    if scope.empty:
        raise ValueError("Route C species scope is empty")
    for column in SPECIES_SCOPE_COLUMNS:
        if scope[column].astype(str).str.strip().eq("").any():
            raise ValueError(f"Route C species scope requires non-empty {column}")
    invalid_distance = sorted(set(scope["species_distance"].astype(str)) - SPECIES_DISTANCE)
    if invalid_distance:
        raise ValueError(f"Route C species scope has invalid species_distance: {', '.join(invalid_distance)}")
    if scope["query_id"].astype(str).str.strip().duplicated().any():
        raise ValueError("Route C species scope query_id must be unique")
    focal_key = str(focal_species).strip().casefold()
    for item in scope.to_dict("records"):
        reference_key = str(item["reference_species"]).strip().casefold()
        distance = str(item["species_distance"]).strip()
        if distance in {"S0_same_cultivar_or_accession", "S1_same_species"} and reference_key != focal_key:
            raise ValueError(f"{distance} requires reference_species to equal focal_species")
        if reference_key == focal_key and distance not in {"S0_same_cultivar_or_accession", "S1_same_species"}:
            raise ValueError("same-species scope must use S0 or S1 species_distance")

    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    cache = load_query_cache(cache_arg)
    provider = "europe_pmc_rest" if allow_external_api else ("local_query_cache" if cache else "none")
    rows: list[dict[str, object]] = []
    status_rows: list[dict[str, object]] = []
    for item in scope.to_dict("records"):
        query_id = str(item["query_id"]).strip()
        reference_species = str(item["reference_species"]).strip()
        phenotype_query = str(item["phenotype_query"]).strip()
        query = _literature_query(phenotype_query, reference_species)
        key = safe_id(query)
        provider_error = ""
        if key in cache:
            articles = cache[key]
            query_status = "retrieved" if articles else "no_hit"
            retrieved_at = _now()
        elif allow_external_api:
            try:
                articles = europe_pmc_fetch(query, limit)
                query_status = "retrieved" if articles else "no_hit"
            except RuntimeError as error:
                articles = []
                query_status = "not_evaluated_provider_unavailable"
                provider_error = str(error)
            retrieved_at = _now()
        else:
            articles = []
            query_status = "not_evaluated_offline"
            retrieved_at = ""
        status_rows.append({
            "query_id": query_id,
            "reference_species": reference_species,
            "species_distance": str(item["species_distance"]),
            "phenotype_query": phenotype_query,
            "query": query,
            "query_status": query_status,
            "article_count": len(articles),
            "provider_error": provider_error,
        })
        for rank, article in enumerate(articles[:limit], start=1):
            pmid = str(article.get("pmid", "") or article.get("id", ""))
            pmcid = str(article.get("pmcid", ""))
            doi = str(article.get("doi", ""))
            source_id = f"PMID:{pmid}" if pmid else (f"PMCID:{pmcid}" if pmcid else (f"DOI:{doi}" if doi else ""))
            rows.append({
                "query_id": query_id,
                "query_phenotype": str(phenotype).strip(),
                "phenotype_query": phenotype_query,
                "focal_species": str(focal_species).strip(),
                "reference_species": reference_species,
                "species_distance": str(item["species_distance"]),
                "query": query,
                "query_status": query_status,
                "article_rank": rank,
                "source_id": source_id,
                "pmid": pmid,
                "pmcid": pmcid,
                "doi": doi,
                "title": re.sub(r"\s+", " ", str(article.get("title", ""))).strip(),
                "abstract_text": re.sub(r"\s+", " ", str(article.get("abstractText", ""))).strip(),
                "journal": str(article.get("journalTitle", "")),
                "year": str(article.get("pubYear", "")),
                "first_author": str(article.get("authorString", "")).split(",")[0],
                "provider": provider,
                "retrieved_at": retrieved_at,
                "rank_contribution": "none",
            })
    candidates = pd.DataFrame(rows, columns=SEARCH_COLUMNS)
    statuses = pd.DataFrame(status_rows)
    candidates_path = out / "route_c_direct_search_candidates.tsv"
    status_path = out / "route_c_direct_search_status.tsv"
    candidates.to_csv(candidates_path, sep="\t", index=False, lineterminator="\n")
    statuses.to_csv(status_path, sep="\t", index=False, lineterminator="\n")
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "command": "route-c-direct-search",
        "candidate_qtl_version": __version__,
        "case_id": case_id,
        "route": "C_direct_known_report_rescue",
        "interval_blind": True,
        "provider": provider,
        "external_api_allowed": allow_external_api,
        "external_transmission": (
            "phenotype query wording and each explicit reference species; no interval, candidate id, or sequence"
            if allow_external_api else "none"
        ),
        "query_count": int(len(scope)),
        "candidate_article_count": int(len(candidates)),
        "candidates_sha256": sha256(candidates_path),
        "status_sha256": sha256(status_path),
        "query_status_counts": {str(k): int(v) for k, v in statuses["query_status"].value_counts().items()},
        "retrieval_is_not_claim_verification": True,
        "rank_contribution": "none",
        "inputs": {"species_scope": file_fingerprint(scope_path, full_hash)},
        "outputs": {"candidates": str(candidates_path), "status": str(status_path)},
    }
    manifest_path = out / "route_c_direct_search_manifest.json"
    manifest["outputs"]["manifest"] = str(manifest_path)
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"candidates": str(candidates_path), "status": str(status_path), "manifest": str(manifest_path),
            "query_count": len(scope), "candidate_article_count": len(candidates)}


def _proposal_claim_table(proposals: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in proposals.to_dict("records"):
        identity = ";".join(filter(None, [
            str(row.get("reported_gene_id", "")).strip(),
            str(row.get("reported_gene_symbol", "")).strip(),
            str(row.get("reported_locus", "")).strip(),
        ]))
        rows.append({
            "record_id": str(row["record_id"]).strip(),
            "source_ids": str(row["source_id"]).strip(),
            "reference_gene_name": identity,
            "gene_family": "",
            "pathway_or_process": str(row["reported_phenotype"]).strip(),
        })
    return pd.DataFrame(rows)


def _tokens(value: object) -> set[str]:
    return {token.casefold() for token in re.split(r"[;,|\s]+", str(value)) if token.strip()}


def _short_identity_in_source(proposal: dict[str, str], supported_row: dict[str, str], source_records: dict[str, object]) -> bool:
    """Check explicit short identifiers in this accepted citation's frozen text.

    The generic claim lexicon intentionally omits two-letter words. This exact
    case-sensitive identity check neither expands that lexicon nor infers an
    experimental relation; its caller supplies only already-supported rows.
    """
    identities = {
        token.strip() for field in ("reported_gene_id", "reported_gene_symbol", "reported_locus")
        for token in re.split(r"[;,|]+", str(proposal.get(field, "")))
        if re.fullmatch(r"[A-Za-z]{2}", token.strip())
    }
    if not identities:
        return False
    citation = str(supported_row.get("citation", ""))
    record = source_records.get(citation)
    if not isinstance(record, dict) or record.get("citation") != citation:
        return False
    if str(record.get("source_record_id", "")) != str(supported_row.get("resolved_id", "")) or str(proposal["record_id"]) not in record.get("hypothesis_ids", []):
        return False
    sections = record.get("sections")
    if not isinstance(sections, dict) or content_sha256(sections) != record.get("sections_sha256"):
        raise ValueError("Route C source review sections fingerprint mismatch")
    return any(
        re.search(r"(?<!\w)" + re.escape(identity) + r"(?!\w)", sections.get(section, "")) is not None
        for identity in identities for section in ("title", "abstractText")
        if isinstance(sections.get(section), str)
    )


def direct_rescue_freeze(
    proposals_arg: str,
    out_dir_arg: str,
    citation_cache_arg: str | None = None,
    allow_external_api: bool = False,
    claim_extraction_arg: str | None = None,
    full_hash: bool = False,
) -> dict[str, object]:
    """Verify agent/researcher-extracted direct reports and freeze Route C."""
    proposals_path = runtime_path(proposals_arg)
    proposals = _require_exact(
        pd.read_csv(proposals_path, sep="\t", dtype=str, keep_default_na=False), DIRECT_PROPOSAL_COLUMNS, "Route C direct proposal table"
    )
    if proposals.empty:
        raise ValueError("Route C direct proposal table is empty")
    for column in ("record_id", "query_phenotype", "focal_species", "reported_phenotype",
                   "reference_species", "species_distance", "evidence_directness", "evidence_type", "source_id"):
        if proposals[column].astype(str).str.strip().eq("").any():
            raise ValueError(f"Route C direct proposal table requires non-empty {column}")
    if proposals["record_id"].astype(str).str.strip().duplicated().any():
        raise ValueError("Route C direct proposal record_id must be unique")
    identity_present = proposals[["reported_gene_id", "reported_gene_symbol", "reported_locus", "reported_variant"]] \
        .astype(str).apply(lambda row: row.str.strip().ne("").any(), axis=1)
    if not identity_present.all():
        raise ValueError("Route C direct proposal requires a reported gene/locus/variant identifier")

    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    adapter_path = out / "source_verify_input.tsv"
    _proposal_claim_table(proposals).to_csv(adapter_path, sep="\t", index=False, lineterminator="\n")
    verification_dir = out / "source_verification"
    verification = source_verify(
        str(adapter_path), str(verification_dir), citation_cache_arg, allow_external_api,
        id_column="record_id", source_column="source_ids", full_hash=full_hash,
        claim_extraction_arg=claim_extraction_arg,
    )
    validation = pd.read_csv(verification["source_validation"], sep="\t", dtype=str).fillna("")
    source_records_path = runtime_path(verification["source_review_records"])
    source_records = json.loads(source_records_path.read_text(encoding="utf-8")).get("records", {})
    if not isinstance(source_records, dict):
        raise ValueError("Route C source review records must be an object")
    grouped = {str(key): group for key, group in validation.groupby("hypothesis_id", sort=False)}
    direct_rows: list[dict[str, str]] = []
    for proposal in proposals.to_dict("records"):
        record_id = str(proposal["record_id"]).strip()
        checked = grouped.get(record_id, pd.DataFrame())
        verified_rows = checked[checked.get("citation_verified", pd.Series(dtype=str)).astype(str).eq("True")] \
            if not checked.empty else pd.DataFrame()
        support_states = set(verified_rows.get("claim_support_status", pd.Series(dtype=str)).astype(str))
        directness = str(proposal["evidence_directness"]).strip()
        accepted_states = set(_DIRECT_SUPPORT)
        if directness == "D2_direct_trait_gene_association":
            accepted_states.add(_ASSOCIATIVE_SUPPORT)
        supported = verified_rows[verified_rows.get("claim_support_status", pd.Series(dtype=str)).isin(accepted_states)] \
            if not verified_rows.empty else pd.DataFrame()
        identity_terms = _tokens(";".join(filter(None, [
            str(proposal.get("reported_gene_id", "")), str(proposal.get("reported_gene_symbol", "")),
            str(proposal.get("reported_locus", "")),
        ])))
        # A citation-verified source whose abstract matched no direct-test cue
        # is undecided, not negative. Its identity can still be checked
        # against the same verified source, so the record keeps a review-only
        # projection path; support still requires a supported claim below. The
        # short-symbol rescue stays limited to sources with an accepted relation.
        cue_supported = not supported.empty
        identity_basis = supported if cue_supported else verified_rows
        identity_rows = identity_basis[[
            bool(identity_terms & _tokens(row.get("claim_terms_found", "")))
            or (cue_supported and _short_identity_in_source(proposal, row, source_records))
            for row in identity_basis.to_dict("records")
        ]] if not identity_basis.empty else pd.DataFrame()
        identity_verified = not identity_rows.empty
        if not supported.empty:
            source_state = "verified"
        elif not checked.empty and checked["verification_status"].eq("rejected").all():
            source_state = "failed"
        elif not checked.empty:
            source_state = "pending"
        else:
            source_state = "not_evaluated"
        if identity_verified and (str(proposal.get("reported_gene_id", "")).strip() or
                                  str(proposal.get("reported_gene_symbol", "")).strip()):
            gene_state = "verified_gene"
        elif identity_verified and str(proposal.get("reported_locus", "")).strip():
            gene_state = "verified_locus_only"
        else:
            gene_state = "unresolved" if source_state != "not_evaluated" else "not_evaluated"
        if support_states & _DIRECT_SUPPORT:
            trait_state = "direct_target_trait"
        elif _ASSOCIATIVE_SUPPORT in support_states:
            trait_state = "direct_related_trait"
        elif directness == "D3_context_only":
            trait_state = "context_only"
        else:
            trait_state = "unresolved"
        if directness == "D3_context_only":
            eligibility = "context_only"
        elif source_state == "verified" and gene_state in {"verified_gene", "verified_locus_only"} and \
                trait_state in {"direct_target_trait", "direct_related_trait"}:
            eligibility = "eligible_direct_record"
        elif source_state == "failed":
            eligibility = "rejected"
        else:
            eligibility = "pending_verification"
        source_row = identity_rows.iloc[0].to_dict() if not identity_rows.empty else (supported.iloc[0].to_dict() if not supported.empty else (
            verified_rows.iloc[0].to_dict() if not verified_rows.empty else (checked.iloc[0].to_dict() if not checked.empty else {})
        ))
        citation = str(source_row.get("citation", proposal["source_id"]))
        citation_kind = str(source_row.get("citation_type", ""))
        identifier = re.sub(r"^(PMID|PMCID|DOI):", "", citation, flags=re.IGNORECASE)
        direct_rows.append({
            **{column: str(proposal.get(column, "")) for column in DIRECT_PROPOSAL_COLUMNS if column in DIRECT_RESCUE_COLUMNS},
            "record_id": record_id,
            "source_id": str(proposal["source_id"]),
            "pmid": identifier if citation_kind == "pmid" else "",
            "pmcid": identifier if citation_kind == "pmcid" else "",
            "doi": str(source_row.get("doi", "")) or (identifier if citation_kind == "doi" else ""),
            "source_title": str(source_row.get("title", "")),
            "source_year": str(source_row.get("year", "")),
            "source_verification_state": source_state,
            "gene_identity_state": gene_state,
            "trait_relation_state": trait_state,
            "direct_rescue_eligibility": eligibility,
            "rank_contribution": "none",
            "notes": "; ".join(filter(None, [str(proposal.get("notes", "")),
                f"claim_support={','.join(sorted(support_states))}" if support_states else "claim_support=not_evaluated"])),
        })
    direct = pd.DataFrame(direct_rows, columns=DIRECT_RESCUE_COLUMNS)
    direct = validate_direct_rescue_frame(direct)
    frozen_path = out / "route_c_direct_rescue_frozen.tsv"
    direct.to_csv(frozen_path, sep="\t", index=False, lineterminator="\n")
    frozen_sha = sha256(frozen_path)
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "command": "route-c-direct-freeze",
        "candidate_qtl_version": __version__,
        "interval_blind": True,
        "source_verification_manifest": verification["manifest"],
        "source_validation_sha256": sha256(runtime_path(verification["source_validation"])),
        "verified_sources_sha256": sha256(runtime_path(verification["verified_sources"])),
        "source_review_records_sha256": sha256(source_records_path),
        "record_count": int(len(direct)),
        "eligibility_counts": {str(k): int(v) for k, v in direct["direct_rescue_eligibility"].value_counts().items()},
        "frozen_sha256": frozen_sha,
        "rank_contribution": "none",
        "inputs": {"proposals": file_fingerprint(proposals_path, full_hash)},
        "outputs": {"frozen": str(frozen_path), "source_verification": verification["manifest"]},
        "gates": [
            "citation resolution and direct claim support are both required",
            "the claimed gene/locus identifier must occur in retrieved source text",
            "explicit two-letter identities require literal-case word-boundary occurrence in the same accepted citation's frozen title/abstract",
            "D3 context never becomes a direct rescue record",
            "the frozen artifact contains no QTL interval or focal candidate fields",
        ],
    }
    manifest_path = out / "route_c_direct_rescue_freeze_manifest.json"
    manifest["outputs"]["manifest"] = str(manifest_path)
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"frozen": str(frozen_path), "manifest": str(manifest_path), "frozen_sha256": frozen_sha,
            "record_count": len(direct), "eligible_count": int(direct_rescue_eligible_count(direct))}


def direct_rescue_eligible_count(frame: pd.DataFrame) -> int:
    return int(frame["direct_rescue_eligibility"].eq("eligible_direct_record").sum())


def _split_identifiers(value: object) -> set[str]:
    return {token.strip().casefold() for token in re.split(r"[;,|]+", str(value)) if token.strip()}


def direct_rescue_project(
    frozen_arg: str,
    freeze_manifest_arg: str,
    atlas_arg: str,
    qtl_id: str,
    out_dir_arg: str,
    eligible_hits_arg: str | None = None,
    full_hash: bool = False,
) -> dict[str, object]:
    """Project frozen direct records into an interval atlas/eligible BLAST hits."""
    if not str(qtl_id).strip():
        raise ValueError("qtl_id is required")
    frozen_path = runtime_path(frozen_arg)
    freeze_manifest_path = runtime_path(freeze_manifest_arg)
    freeze_manifest = json.loads(freeze_manifest_path.read_text(encoding="utf-8"))
    expected_frozen_sha = str(freeze_manifest.get("frozen_sha256", ""))
    if not expected_frozen_sha or expected_frozen_sha != sha256(frozen_path):
        raise ValueError("Route C frozen artifact hash does not match freeze manifest")
    atlas_path = runtime_path(atlas_arg)
    direct = validate_direct_rescue_frame(pd.read_csv(frozen_path, sep="\t", dtype=str).fillna(""))
    atlas = pd.read_csv(atlas_path, sep="\t", dtype=str).fillna("")
    if "gene_id" not in atlas.columns:
        raise ValueError("Route C projection atlas requires gene_id")
    hits = pd.DataFrame()
    hits_path: Path | None = None
    if eligible_hits_arg:
        hits_path = runtime_path(eligible_hits_arg)
        hits = pd.read_csv(hits_path, sep="\t", dtype=str).fillna("")
        missing = sorted({"panel_id", "target_id", "eligibility", "pident", "query_coverage", "target_coverage"} - set(hits.columns))
        if missing:
            raise ValueError(f"eligible interval hits are missing columns: {', '.join(missing)}")
        hits = hits[hits["eligibility"].eq("eligible_interval_hit")]
    interval_hash = f"sha256:{sha256(atlas_path)}"
    projection_rows: list[dict[str, str]] = []
    for record in direct.to_dict("records"):
        record_id = str(record["record_id"])
        if not direct_record_can_seed_projection(record):
            projection_rows.append({
                "record_id": record_id, "qtl_id": "", "focal_gene_id": "",
                "projection_state": "not_projected", "projection_basis": "not_evaluated",
                "source_reference_gene_id": str(record.get("reported_gene_id", "")),
                "identity": "", "query_coverage": "", "target_coverage": "",
                "ambiguity_state": "not_evaluated", "interval_universe_fingerprint": interval_hash,
                "rank_contribution": "none", "notes": "direct record is not projection-eligible",
            })
            continue
        keys = _split_identifiers(";".join(filter(None, [
            str(record.get("reported_gene_id", "")), str(record.get("reported_gene_symbol", "")),
            str(record.get("reference_sequence_id", "")),
        ])))
        matched: dict[str, dict[str, str]] = {}
        for atlas_row in atlas.to_dict("records"):
            gene_id = str(atlas_row["gene_id"]).strip()
            if gene_id.casefold() in keys:
                matched[gene_id] = {"basis": "exact_focal_id", "identity": "", "query_coverage": "", "target_coverage": ""}
                continue
            references = _split_identifiers(atlas_row.get("orthology_reference_gene", ""))
            if keys & references:
                matched[gene_id] = {"basis": "orthology_reference_gene", "identity": "", "query_coverage": "", "target_coverage": ""}
        if not hits.empty:
            for hit in hits[hits["panel_id"].astype(str).eq(record_id)].to_dict("records"):
                target = str(hit["target_id"]).strip()
                if target in set(atlas["gene_id"].astype(str)):
                    matched[target] = {
                        "basis": "blastp", "identity": str(hit["pident"]),
                        "query_coverage": str(hit["query_coverage"]), "target_coverage": str(hit["target_coverage"]),
                    }
        if not matched:
            basis = "blastp" if hits_path else ("orthology_reference_gene" if "orthology_reference_gene" in atlas.columns else "exact_focal_id")
            projection_rows.append({
                "record_id": record_id, "qtl_id": str(qtl_id), "focal_gene_id": "",
                "projection_state": "no_interval_match", "projection_basis": basis,
                "source_reference_gene_id": str(record.get("reported_gene_id", "")),
                "identity": "", "query_coverage": "", "target_coverage": "",
                "ambiguity_state": "no_threshold_hit", "interval_universe_fingerprint": interval_hash,
                "rank_contribution": "none", "notes": "evaluated interval search found no eligible match; outside-QTL presence is not inferred",
            })
            continue
        ambiguity = "unique" if len(matched) == 1 else "one_to_many"
        for gene_id in sorted(matched):
            evidence = matched[gene_id]
            projection_rows.append({
                "record_id": record_id, "qtl_id": str(qtl_id), "focal_gene_id": gene_id,
                "projection_state": "inside_current_qtl", "projection_basis": evidence["basis"],
                "source_reference_gene_id": str(record.get("reported_gene_id", "")),
                "identity": evidence["identity"], "query_coverage": evidence["query_coverage"],
                "target_coverage": evidence["target_coverage"], "ambiguity_state": ambiguity,
                "interval_universe_fingerprint": interval_hash, "rank_contribution": "none",
                "notes": "projection evidence only; cross-species trait support is not automatically transferred",
            })
    projection = validate_projection_frame(pd.DataFrame(projection_rows, columns=PROJECTION_COLUMNS))
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    projection_path = out / "route_c_direct_projection.tsv"
    projection.to_csv(projection_path, sep="\t", index=False, lineterminator="\n")
    evidence = direct.merge(projection, on="record_id", how="left", suffixes=("_direct", "_projection"))
    evidence_path = out / "route_c_direct_evidence_chain.tsv"
    evidence.to_csv(evidence_path, sep="\t", index=False, lineterminator="\n")
    projection_sha = sha256(projection_path)
    evidence_chain_sha = sha256(evidence_path)
    manifest = {
        "schema_version": SCHEMA_VERSION, "command": "route-c-direct-project",
        "candidate_qtl_version": __version__, "retrieval_frozen_before_interval": True,
        "frozen_sha256": sha256(frozen_path), "qtl_id": str(qtl_id),
        "projection_sha256": projection_sha,
        "evidence_chain_sha256": evidence_chain_sha,
        "projection_state_counts": {str(k): int(v) for k, v in projection["projection_state"].value_counts().items()},
        "rank_contribution": "none",
        "trait_transfer_policy": "projection/homology is not automatic trait-specific support",
        "inputs": {"frozen": file_fingerprint(frozen_path, full_hash),
                   "freeze_manifest": file_fingerprint(freeze_manifest_path, full_hash),
                   "atlas": file_fingerprint(atlas_path, full_hash),
                   "eligible_hits": file_fingerprint(hits_path, full_hash) if hits_path else None},
        "outputs": {"projection": str(projection_path), "evidence_chain": str(evidence_path)},
    }
    manifest_path = out / "route_c_direct_projection_manifest.json"
    manifest["outputs"]["manifest"] = str(manifest_path)
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"projection": str(projection_path), "evidence_chain": str(evidence_path), "manifest": str(manifest_path),
            "inside_count": int(projection["projection_state"].eq("inside_current_qtl").sum())}


def direct_rescue_benchmark(manifest_arg: str, out_dir_arg: str, full_hash: bool = False) -> dict[str, object]:
    """Score frozen source recovery and post-freeze projection separately."""
    manifest_path = runtime_path(manifest_arg)
    cases = _require_exact(
        pd.read_csv(manifest_path, sep="\t", dtype=str), DIRECT_BENCHMARK_CASE_COLUMNS,
        "Route C direct benchmark manifest",
    )
    rows: list[dict[str, object]] = []
    for case in cases.to_dict("records"):
        contamination = str(case["contamination_status"]).strip()
        if contamination not in BENCHMARK_CONTAMINATION_STATES:
            raise ValueError(
                "Route C direct benchmark contamination_status must be one of "
                f"{', '.join(sorted(BENCHMARK_CONTAMINATION_STATES))}; got {contamination or '<empty>'}"
            )
        direct_path = runtime_path(case["direct_rescue_tsv"])
        freeze_manifest_path = runtime_path(case["direct_freeze_manifest"])
        freeze_manifest = json.loads(freeze_manifest_path.read_text(encoding="utf-8"))
        freeze_valid = str(freeze_manifest.get("frozen_sha256", "")) == sha256(direct_path)
        if not freeze_valid:
            raise ValueError(f"Route C benchmark frozen artifact hash mismatch for case {case['case_id']}")
        direct = validate_direct_rescue_frame(pd.read_csv(direct_path, sep="\t", dtype=str).fillna(""))
        eligible = direct[direct["direct_rescue_eligibility"].eq("eligible_direct_record")]
        truth_reference = str(case["truth_reference_gene_id"]).strip().casefold()
        clean = contamination == "clean"
        source_evaluated = bool(truth_reference and clean)
        source_hit = bool(source_evaluated and any(
            truth_reference in _split_identifiers(";".join([str(row.get("reported_gene_id", "")), str(row.get("reported_gene_symbol", "")), str(row.get("reference_sequence_id", ""))]))
            for row in eligible.to_dict("records")
        ))
        truth_focal = str(case["truth_focal_gene_id"]).strip()
        projection_evaluated = bool(truth_focal and str(case["projection_tsv"]).strip() and clean)
        projection_hit = False
        if projection_evaluated:
            projection_path = runtime_path(case["projection_tsv"])
            projection_manifest_path = runtime_path(case["projection_manifest"])
            projection_manifest = json.loads(projection_manifest_path.read_text(encoding="utf-8"))
            if str(projection_manifest.get("frozen_sha256", "")) != sha256(direct_path):
                raise ValueError(f"Route C benchmark projection lineage mismatch for case {case['case_id']}")
            if str(projection_manifest.get("projection_sha256", "")) != sha256(projection_path):
                raise ValueError(f"Route C benchmark projection artifact hash mismatch for case {case['case_id']}")
            expected_qtl_id = str(case["expected_qtl_id"]).strip()
            if not expected_qtl_id or str(projection_manifest.get("qtl_id", "")) != expected_qtl_id:
                raise ValueError(f"Route C benchmark projection qtl_id mismatch for case {case['case_id']}")
            expected_atlas_sha = str(case["expected_atlas_sha256"]).strip().removeprefix("sha256:")
            manifest_atlas_sha = str(projection_manifest.get("inputs", {}).get("atlas", {}).get("sha256", ""))
            if not expected_atlas_sha or manifest_atlas_sha != expected_atlas_sha:
                raise ValueError(f"Route C benchmark projection atlas mismatch for case {case['case_id']}")
            projection = validate_projection_frame(pd.read_csv(projection_path, sep="\t", dtype=str).fillna(""))
            if set(projection["qtl_id"].astype(str)) - {"", expected_qtl_id}:
                raise ValueError(f"Route C benchmark projection rows contain another qtl_id for case {case['case_id']}")
            expected_fingerprint = f"sha256:{expected_atlas_sha}"
            if set(projection["interval_universe_fingerprint"].astype(str)) != {expected_fingerprint}:
                raise ValueError(f"Route C benchmark projection rows contain another atlas for case {case['case_id']}")
            projection_hit = bool(((projection["projection_state"] == "inside_current_qtl") &
                                   (projection["focal_gene_id"] == truth_focal)).any())
        rows.append({
            "case_id": str(case["case_id"]), "contamination_status": contamination,
            "source_evaluation_status": "clean" if source_evaluated else (contamination if not clean else "not_evaluated"),
            "source_recovery_hit": source_hit if source_evaluated else "", "eligible_direct_record_count": len(eligible),
            "projection_evaluation_status": "clean" if projection_evaluated else (contamination if not clean else "not_evaluated"),
            "projection_hit": projection_hit if projection_evaluated else "",
        })
    scored = pd.DataFrame(rows)
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    cases_path = out / "route_c_direct_benchmark_cases.tsv"
    scored.to_csv(cases_path, sep="\t", index=False, lineterminator="\n")

    def metric(status_column: str, hit_column: str) -> dict[str, object]:
        evaluated = scored[status_column].eq("clean")
        denominator = int(evaluated.sum())
        hits = int(scored.loc[evaluated, hit_column].astype(str).eq("True").sum())
        return {"evaluated": denominator, "hits": hits, "recall": hits / denominator if denominator else None,
                "not_evaluated": int(scored[status_column].eq("not_evaluated").sum()),
                "contaminated": int(scored[status_column].eq("contaminated").sum())}

    metrics = {
        "source_recovery": metric("source_evaluation_status", "source_recovery_hit"),
        "projection": metric("projection_evaluation_status", "projection_hit"),
    }
    summary = pd.DataFrame([{"stage": key, **value} for key, value in metrics.items()])
    summary_path = out / "route_c_direct_benchmark_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False, lineterminator="\n")
    output_manifest = {
        "schema_version": SCHEMA_VERSION, "command": "route-c-direct-benchmark",
        "case_count": len(scored), "metrics": metrics,
        "denominator_policy": "not_evaluated and contaminated cases are not misses; source recovery and projection have separate denominators",
        "scoring_order": "prediction/freeze -> truth open -> score",
        "inputs": {"manifest": file_fingerprint(manifest_path, full_hash)},
        "outputs": {"cases": str(cases_path), "summary": str(summary_path)},
    }
    output_manifest_path = out / "route_c_direct_benchmark_manifest.json"
    output_manifest["outputs"]["manifest"] = str(output_manifest_path)
    output_manifest_path.write_text(json.dumps(output_manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"cases": str(cases_path), "summary": str(summary_path), "manifest": str(output_manifest_path), "metrics": metrics}
