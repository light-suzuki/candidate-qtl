"""Route-B literature retrieval backed by shared immutable cassettes."""

from __future__ import annotations

import json
import re
from datetime import datetime, timezone

import pandas as pd

from . import __version__
from .core import file_fingerprint, runtime_path, sha256
from .literature_cassette import (
    CASSETTE_CONTRACT,
    EUROPE_PMC_RESULT_TYPE,
    LiteratureRetrieval,
    retrieve_europe_pmc,
    write_query_provenance,
)
from .identifiers import is_usable_symbol
from .reverse_literature_legacy import (
    LITERATURE_COLUMNS,
    _phenotype_terms,
    _query_for,
    _searchable_reference,
)

SCHEMA_VERSION = "1.1"
QUERY_STATES = (
    "retrieved",
    "no_hit",
    "not_evaluated_offline",
    "not_evaluated_rate_limited",
    "not_evaluated_provider_unavailable",
    "not_evaluated_response_parse_error",
    "invalid_request",
    "skipped_no_reference_gene",
    "skipped_unusable_reference_identifier",
    "not_evaluated_conflicting_reference_gene_symbols",
)
CASSETTE_STATUS_COLUMNS = [
    "reference_gene", "reference_gene_symbol", "query_term", "identifier_source",
    "interval_gene_count", "query_status", "article_count", "note",
]


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


def reverse_literature(
    lookup_arg: str,
    out_dir_arg: str,
    phenotype: str = "",
    cache_arg: str | None = None,
    allow_external_api: bool = False,
    limit: int = 5,
    case_id: str = "",
    full_hash: bool = False,
) -> dict[str, object]:
    """Query each interval reference with per-query cassette provenance."""
    if limit < 1:
        raise ValueError("limit must be at least 1")
    lookup_path = runtime_path(lookup_arg)
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    cassette_dir = out / "literature_cassettes"
    lookup = pd.read_csv(lookup_path, sep="\t", dtype=str).fillna("")
    missing = sorted({"gene_id", "reference_gene"} - set(lookup.columns))
    if missing:
        raise ValueError(
            f"interval reverse-lookup table is missing columns: {', '.join(missing)}"
        )

    terms = _phenotype_terms(phenotype)
    by_reference: dict[str, dict[str, object]] = {}
    for row in lookup.to_dict("records"):
        reference = str(row.get("reference_gene", "")).strip()
        record = by_reference.setdefault(
            reference,
            {
                "genes": [],
                "function": str(row.get("reference_function", "")),
                "symbols": set(),
            },
        )
        record["genes"].append(str(row["gene_id"]))  # type: ignore[union-attr]
        if not record["function"]:
            record["function"] = str(row.get("reference_function", ""))
        symbol = str(row.get("reference_gene_symbol", "")).strip()
        if is_usable_symbol(symbol):
            record["symbols"].add(symbol)  # type: ignore[union-attr]

    rows: list[dict[str, object]] = []
    status_rows: list[dict[str, object]] = []
    retrievals: list[LiteratureRetrieval] = []
    for reference in sorted(by_reference, key=lambda value: (value.casefold(), value)):
        genes = sorted(by_reference[reference]["genes"])  # type: ignore[arg-type]
        function = str(by_reference[reference]["function"])
        symbols = sorted(by_reference[reference]["symbols"], key=lambda value: (value.casefold(), value))  # type: ignore[arg-type]
        if not reference:
            status_rows.append(
                {
                    "reference_gene": "",
                    "reference_gene_symbol": "",
                    "query_term": "",
                    "identifier_source": "none",
                    "interval_gene_count": len(genes),
                    "query_status": "skipped_no_reference_gene",
                    "article_count": 0,
                    "note": (
                        "these interval genes have no characterised relative; "
                        "they stay an unresolved bucket"
                    ),
                }
            )
            continue
        if len(symbols) > 1:
            status_rows.append(
                {
                    "reference_gene": reference,
                    "reference_gene_symbol": ";".join(symbols),
                    "query_term": "",
                    "identifier_source": "conflicting_reference_gene_symbols",
                    "interval_gene_count": len(genes),
                    "query_status": "not_evaluated_conflicting_reference_gene_symbols",
                    "article_count": 0,
                    "note": "multiple usable symbols map to this reference identifier; no arbitrary symbol was searched",
                }
            )
            continue
        if symbols:
            query_term, identifier_source = symbols[0], "reference_gene_symbol"
        else:
            query_term, identifier_source = _searchable_reference(reference)
        if not query_term:
            status_rows.append(
                {
                    "reference_gene": reference,
                    "reference_gene_symbol": symbols[0] if symbols else "",
                    "query_term": "",
                    "identifier_source": identifier_source,
                    "interval_gene_count": len(genes),
                    "query_status": "skipped_unusable_reference_identifier",
                    "article_count": 0,
                    "note": (
                        "no symbol or accession could be derived from this identifier; "
                        "no search was performed and absence of literature is not established"
                    ),
                }
            )
            continue
        query = _query_for(query_term, terms)
        retrieval = retrieve_europe_pmc(
            query,
            limit,
            cache_arg=cache_arg,
            write_dir_arg=cassette_dir,
            allow_external_api=allow_external_api,
            user_agent=(
                f"candidate-qtl/{__version__} reverse-literature "
                "(read-only literature lookup)"
            ),
        )
        retrievals.append(retrieval)
        query_status = _query_status(retrieval)
        status_rows.append(
            {
                "reference_gene": reference,
                "reference_gene_symbol": symbols[0] if symbols else "",
                "query_term": query_term,
                "identifier_source": identifier_source,
                "interval_gene_count": len(genes),
                "query_status": query_status,
                "article_count": len(retrieval.articles),
                "note": (
                    retrieval.error
                    or (
                        "no search was performed; absence of literature is not established"
                        if query_status == "not_evaluated_offline"
                        else ""
                    )
                ),
            }
        )
        for rank, article in enumerate(retrieval.articles[:limit], start=1):
            blob = " ".join(
                str(article.get(key, "")) for key in ("title", "abstractText")
            ).casefold()
            matched = [term for term in terms if term in blob]
            identifier = str(article.get("pmid", "") or article.get("id", ""))
            rows.append(
                {
                    "interval_gene_ids": ";".join(genes),
                    "reference_gene": reference,
                    "reference_function": function,
                    "query": query,
                    "query_status": query_status,
                    "article_rank": rank,
                    "source_id": f"PMID:{identifier}" if identifier else "",
                    "title": re.sub(
                        r"\s+", " ", str(article.get("title", ""))
                    ).strip(),
                    "journal": str(article.get("journalTitle", "")),
                    "year": str(article.get("pubYear", "")),
                    "doi": str(article.get("doi", "")),
                    "first_author": str(article.get("authorString", "")).split(",")[0],
                    "matched_phenotype_terms": ";".join(matched),
                    "abstract_available": str(
                        bool(str(article.get("abstractText", "")).strip())
                    ),
                    "claim_type": "literature_link_candidate",
                    "provider": _provider_label(retrieval),
                    "retrieved_at": retrieval.provider_retrieved_at,
                    "rank_contribution": "none",
                }
            )

    literature = pd.DataFrame(rows, columns=LITERATURE_COLUMNS)
    status = pd.DataFrame(status_rows, columns=CASSETTE_STATUS_COLUMNS)
    if not status.empty:
        status = status.sort_values(
            "reference_gene",
            kind="stable",
            key=lambda series: series.astype(str).str.casefold(),
        ).reset_index(drop=True)
    literature_path = out / "reverse_literature.tsv"
    status_path = out / "reverse_literature_status.tsv"
    provenance_path = out / "reverse_literature_query_provenance.tsv"
    literature.to_csv(literature_path, sep="\t", index=False, lineterminator="\n")
    status.to_csv(status_path, sep="\t", index=False, lineterminator="\n")
    write_query_provenance(retrievals, provenance_path)

    retrieved_at_values = sorted(
        {
            retrieval.provider_retrieved_at
            for retrieval in retrievals
            if retrieval.provider_retrieved_at
        }
    )
    source_counts: dict[str, int] = {}
    for retrieval in retrievals:
        source_counts[retrieval.retrieval_source] = (
            source_counts.get(retrieval.retrieval_source, 0) + 1
        )
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "command": "reverse-literature",
        "case_id": case_id,
        "candidate_qtl_version": __version__,
        "queried_at_utc": retrieved_at_values[-1] if retrieved_at_values else "",
        "provider_retrieved_at_values": retrieved_at_values,
        "route": "interval_first",
        "stage": "reference_gene_to_literature",
        "provider": "per_query",
        "retrieval_source_counts": dict(sorted(source_counts.items())),
        "external_api_allowed": allow_external_api,
        "result_type": EUROPE_PMC_RESULT_TYPE,
        "cassette_contract": CASSETTE_CONTRACT,
        "abstract_available_article_count": (
            int(literature["abstract_available"].eq("True").sum())
            if not literature.empty
            else 0
        ),
        "phenotype_terms": terms,
        "reference_gene_count": int(len(status[status["reference_gene"] != ""])),
        "queried": int((status["query_status"] == "retrieved").sum()),
        "query_attempted": int(
            status["query_status"]
            .isin(
                {
                    "retrieved",
                    "no_hit",
                    "not_evaluated_rate_limited",
                    "not_evaluated_provider_unavailable",
                    "not_evaluated_response_parse_error",
                    "invalid_request",
                }
            )
            .sum()
        ),
        "no_hit": int((status["query_status"] == "no_hit").sum()),
        "not_evaluated_offline": int(
            (status["query_status"] == "not_evaluated_offline").sum()
        ),
        "provider_failure_count": int(
            status["query_status"]
            .isin(
                {
                    "not_evaluated_rate_limited",
                    "not_evaluated_provider_unavailable",
                    "not_evaluated_response_parse_error",
                }
            )
            .sum()
        ),
        "skipped_unusable_reference_identifier": int(
            (
                status["query_status"]
                == "skipped_unusable_reference_identifier"
            ).sum()
        ),
        "not_evaluated_conflicting_reference_gene_symbols": int(
            (
                status["query_status"]
                == "not_evaluated_conflicting_reference_gene_symbols"
            ).sum()
        ),
        "identifier_source_counts": {
            str(key): int(value)
            for key, value in status["identifier_source"]
            .value_counts()
            .sort_index()
            .items()
        },
        "interval_genes_without_reference": (
            int(
                status[status["reference_gene"] == ""][
                    "interval_gene_count"
                ].sum()
            )
            if len(status)
            else 0
        ),
        "article_rows": int(len(literature)),
        "external_transmission": (
            "reference gene names and generic phenotype words only; interval "
            "coordinates, candidate ids and sequences are never transmitted"
            if allow_external_api
            else "none"
        ),
        "claims": [
            "a returned article is a link candidate, not evidence: retrieval is not relevance",
            "an offline or provider-failed query is not_evaluated, never 'no literature exists'",
            "interval genes without a characterised relative stay an explicit unresolved bucket",
            "nothing here writes a score",
        ],
        "next_step": (
            "judge each link candidate, then source-verify the citations that are kept"
        ),
        "inputs": {"reverse_lookup": file_fingerprint(lookup_path, full_hash)},
        "outputs": {
            "literature": str(literature_path),
            "status": str(status_path),
            "query_provenance": str(provenance_path),
            "cassette_directory": str(cassette_dir),
        },
        "output_sha256": {
            "reverse_literature": sha256(literature_path),
            "status": sha256(status_path),
            "query_provenance": sha256(provenance_path),
        },
    }
    manifest_path = out / "reverse_literature_manifest.json"
    manifest["outputs"]["manifest"] = str(manifest_path)
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    run_envelope_path = out / "reverse_literature_run_envelope.json"
    run_envelope_path.write_text(
        json.dumps(
            {
                "command": "reverse-literature",
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
        "literature": str(literature_path),
        "status": str(status_path),
        "query_provenance": str(provenance_path),
        "manifest": str(manifest_path),
        "run_envelope": str(run_envelope_path),
        "article_rows": len(literature),
    }
