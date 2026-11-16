"""Direct-C interval-blind search backed by shared literature cassettes."""

from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from . import __version__
from .core import file_fingerprint, runtime_path, sha256
from .literature_cassette import (
    CASSETTE_CONTRACT,
    LiteratureRetrieval,
    retrieve_europe_pmc,
    write_query_provenance,
)
from .route_c_direct_pipeline_legacy import (
    SCHEMA_VERSION,
    SEARCH_COLUMNS,
    SPECIES_SCOPE_COLUMNS,
    _literature_query,
    _require_exact,
)
from .route_c_direct_rescue import SPECIES_DISTANCE

STATUS_COLUMNS = (
    "query_id",
    "reference_species",
    "species_distance",
    "phenotype_query",
    "query",
    "query_status",
    "article_count",
    "provider",
    "retrieval_source",
    "provider_retrieved_at",
    "request_key_sha256",
    "raw_response_sha256",
    "cassette_path",
    "provider_error",
)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _query_status(retrieval: LiteratureRetrieval) -> str:
    return {
        "retrieved": "retrieved",
        "no_hit": "no_hit",
        "not_evaluated_offline": "not_evaluated_offline",
        "rate_limited": "not_evaluated_rate_limited",
        "provider_unavailable": "not_evaluated_provider_unavailable",
        "response_parse_error": "not_evaluated_response_parse_error",
        "invalid_request": "invalid_request",
    }[retrieval.retrieval_status]


def _provider_label(retrieval: LiteratureRetrieval) -> str:
    return f"{retrieval.provider}:{retrieval.retrieval_source}"


def _generation_sha256(*digests: str) -> str:
    return hashlib.sha256(
        json.dumps(
            list(digests), ensure_ascii=False, separators=(",", ":")
        ).encode("utf-8")
    ).hexdigest()


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
    """Search direct trait-gene reports with immutable per-query provenance."""
    if not str(phenotype).strip() or not str(focal_species).strip():
        raise ValueError("phenotype and focal_species are required")
    if limit < 1:
        raise ValueError("limit must be at least 1")
    scope_path = runtime_path(species_scope_arg)
    scope = _require_exact(
        pd.read_csv(scope_path, sep="\t", dtype=str),
        SPECIES_SCOPE_COLUMNS,
        "Route C species scope",
    )
    if scope.empty:
        raise ValueError("Route C species scope is empty")
    for column in SPECIES_SCOPE_COLUMNS:
        if scope[column].astype(str).str.strip().eq("").any():
            raise ValueError(
                f"Route C species scope requires non-empty {column}"
            )
    invalid_distance = sorted(
        set(scope["species_distance"].astype(str)) - SPECIES_DISTANCE
    )
    if invalid_distance:
        raise ValueError(
            "Route C species scope has invalid species_distance: "
            + ", ".join(invalid_distance)
        )
    if scope["query_id"].astype(str).str.strip().duplicated().any():
        raise ValueError("Route C species scope query_id must be unique")
    focal_key = str(focal_species).strip().casefold()
    for item in scope.to_dict("records"):
        reference_key = str(item["reference_species"]).strip().casefold()
        distance = str(item["species_distance"]).strip()
        if (
            distance in {"S0_same_cultivar_or_accession", "S1_same_species"}
            and reference_key != focal_key
        ):
            raise ValueError(
                f"{distance} requires reference_species to equal focal_species"
            )
        if reference_key == focal_key and distance not in {
            "S0_same_cultivar_or_accession",
            "S1_same_species",
        }:
            raise ValueError(
                "same-species scope must use S0 or S1 species_distance"
            )

    scope = scope.sort_values(
        "query_id",
        kind="stable",
        key=lambda series: series.astype(str).str.casefold(),
    ).reset_index(drop=True)
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    cassette_dir = out / "literature_cassettes"
    rows: list[dict[str, object]] = []
    status_rows: list[dict[str, object]] = []
    retrievals: list[LiteratureRetrieval] = []
    for item in scope.to_dict("records"):
        query_id = str(item["query_id"]).strip()
        reference_species = str(item["reference_species"]).strip()
        phenotype_query = str(item["phenotype_query"]).strip()
        query = _literature_query(phenotype_query, reference_species)
        retrieval = retrieve_europe_pmc(
            query,
            limit,
            cache_arg=cache_arg,
            write_dir_arg=cassette_dir,
            allow_external_api=allow_external_api,
            user_agent=(
                f"candidate-qtl/{__version__} route-c-direct-search "
                "(read-only literature lookup)"
            ),
        )
        retrievals.append(retrieval)
        query_status = _query_status(retrieval)
        provider_label = _provider_label(retrieval)
        status_rows.append(
            {
                "query_id": query_id,
                "reference_species": reference_species,
                "species_distance": str(item["species_distance"]),
                "phenotype_query": phenotype_query,
                "query": query,
                "query_status": query_status,
                "article_count": len(retrieval.articles),
                "provider": retrieval.provider,
                "retrieval_source": retrieval.retrieval_source,
                "provider_retrieved_at": retrieval.provider_retrieved_at,
                "request_key_sha256": retrieval.request_key_sha256,
                "raw_response_sha256": retrieval.raw_response_sha256,
                "cassette_path": (
                    Path(retrieval.cassette_path).name
                    if retrieval.cassette_path
                    else ""
                ),
                "provider_error": retrieval.error,
            }
        )
        for rank, article in enumerate(retrieval.articles[:limit], start=1):
            pmid = str(article.get("pmid", "") or article.get("id", ""))
            pmcid = str(article.get("pmcid", ""))
            doi = str(article.get("doi", ""))
            source_id = (
                f"PMID:{pmid}"
                if pmid
                else (
                    f"PMCID:{pmcid}"
                    if pmcid
                    else (f"DOI:{doi}" if doi else "")
                )
            )
            rows.append(
                {
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
                    "title": re.sub(
                        r"\s+", " ", str(article.get("title", ""))
                    ).strip(),
                    "abstract_text": re.sub(
                        r"\s+", " ", str(article.get("abstractText", ""))
                    ).strip(),
                    "journal": str(article.get("journalTitle", "")),
                    "year": str(article.get("pubYear", "")),
                    "first_author": str(article.get("authorString", "")).split(",")[0],
                    "provider": provider_label,
                    "retrieved_at": retrieval.provider_retrieved_at,
                    "rank_contribution": "none",
                }
            )

    candidates = pd.DataFrame(rows, columns=SEARCH_COLUMNS)
    statuses = pd.DataFrame(status_rows, columns=STATUS_COLUMNS)
    candidates_path = out / "route_c_direct_search_candidates.tsv"
    status_path = out / "route_c_direct_search_status.tsv"
    provenance_path = out / "route_c_direct_search_query_provenance.tsv"
    candidates.to_csv(
        candidates_path, sep="\t", index=False, lineterminator="\n"
    )
    statuses.to_csv(status_path, sep="\t", index=False, lineterminator="\n")
    write_query_provenance(retrievals, provenance_path)
    candidates_sha = sha256(candidates_path)
    status_sha = sha256(status_path)
    provenance_sha = sha256(provenance_path)
    source_counts: dict[str, int] = {}
    for retrieval in retrievals:
        source_counts[retrieval.retrieval_source] = (
            source_counts.get(retrieval.retrieval_source, 0) + 1
        )
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "command": "route-c-direct-search",
        "candidate_qtl_version": __version__,
        "case_id": case_id,
        "route": "C_direct_known_report_rescue",
        "interval_blind": True,
        "provider": "per_query",
        "retrieval_source_counts": dict(sorted(source_counts.items())),
        "cassette_contract": CASSETTE_CONTRACT,
        "external_api_allowed": allow_external_api,
        "external_transmission": (
            "phenotype query wording and each explicit reference species; no "
            "interval, candidate id, or sequence"
            if allow_external_api
            else "none"
        ),
        "query_count": int(len(scope)),
        "candidate_article_count": int(len(candidates)),
        "candidates_sha256": candidates_sha,
        "status_sha256": status_sha,
        "query_provenance_sha256": provenance_sha,
        "search_generation_sha256": _generation_sha256(
            candidates_sha, status_sha, provenance_sha
        ),
        "query_status_counts": {
            str(key): int(value)
            for key, value in statuses["query_status"]
            .value_counts()
            .sort_index()
            .items()
        },
        "retrieval_is_not_claim_verification": True,
        "rank_contribution": "none",
        "inputs": {"species_scope": file_fingerprint(scope_path, full_hash)},
        "outputs": {
            "candidates": str(candidates_path),
            "status": str(status_path),
            "query_provenance": str(provenance_path),
            "cassette_directory": str(cassette_dir),
        },
    }
    manifest_path = out / "route_c_direct_search_manifest.json"
    manifest["outputs"]["manifest"] = str(manifest_path)
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    run_envelope_path = out / "route_c_direct_search_run_envelope.json"
    run_envelope_path.write_text(
        json.dumps(
            {
                "command": "route-c-direct-search",
                "run_started_at_utc": _now(),
                "scientific_manifest": str(manifest_path),
                "scientific_manifest_sha256": sha256(manifest_path),
                "query_execution_sources": [
                    {
                        "request_key_sha256": retrieval.request_key_sha256,
                        "execution_source": retrieval.execution_source,
                        "retrieval_status": retrieval.retrieval_status,
                    }
                    for retrieval in sorted(
                        retrievals, key=lambda item: item.request_key_sha256
                    )
                ],
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return {
        "candidates": str(candidates_path),
        "status": str(status_path),
        "query_provenance": str(provenance_path),
        "manifest": str(manifest_path),
        "run_envelope": str(run_envelope_path),
        "query_count": len(scope),
        "candidate_article_count": len(candidates),
    }
