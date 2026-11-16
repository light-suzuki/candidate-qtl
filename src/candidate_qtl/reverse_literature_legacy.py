"""Route B, second stage: go from the interval's genes out to the literature.

``route-b-prepare`` stops at *what* each interval gene's closest characterised
relative is called:

    Psat_A -> Q9SYQ8 -> "Receptor protein kinase CLAVATA1" (UniProt:Q9SYQ8)

That is an annotation join, not a literature search, and it is where the reverse
route was quietly incomplete: the point of going gene -> literature is to find
out **what happens when that gene is broken**, which a description never says.

This stage takes the reference gene of every interval gene and asks the
literature about it, focused by the phenotype. The work is bounded by the number
of genes in the interval, so it is a few hundred queries at worst and it caches.

What it does and does not claim:

* A returned article is a **link candidate**, not evidence. Retrieval is not
  relevance: a paper mentioning CLV1 and fasciation may still be about something
  else. Every row is ``literature_link_candidate`` and needs a human or an agent
  to judge it.
* Nothing here writes a score. The rows feed the same review path as every other
  proposal.
* Offline without a cache is ``not_evaluated``, never "no literature exists".
  A search that did not run is not an absence.
* Only the reference gene name and generic phenotype words are sent. The interval
  coordinates, the candidate ids and every sequence stay local.
"""

from __future__ import annotations

import json
import re
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone

import pandas as pd

from . import __version__
from .core import file_fingerprint, runtime_path, safe_id, sha256
from .identifiers import canonical_accession, is_uniprot_accession, is_usable_symbol

SCHEMA_VERSION = "1.0"
EUROPE_PMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
USER_AGENT = f"candidate-qtl/{__version__} reverse-literature (read-only literature lookup)"

LITERATURE_COLUMNS = [
    "interval_gene_ids", "reference_gene", "reference_function", "query", "query_status",
    "article_rank", "source_id", "title", "journal", "year", "doi", "first_author",
    "matched_phenotype_terms", "abstract_available", "claim_type", "provider", "retrieved_at", "rank_contribution",
]
STATUS_COLUMNS = ["reference_gene", "query_term", "identifier_source", "interval_gene_count",
                  "query_status", "article_count", "note"]

QUERY_STATES = ("retrieved", "no_hit", "not_evaluated_offline", "skipped_no_reference_gene",
                "skipped_unusable_reference_identifier")

# Generic words that would match half of plant biology on their own.
_STOPWORDS = frozenset(
    "the a an and or of to in on with without at is are be increased decreased reduced plant plants "
    "phenotype trait gene genes protein mutant abnormal altered".split()
)
_TERM = re.compile(r"[A-Za-z][A-Za-z0-9-]{2,}")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _phenotype_terms(phenotype: str) -> list[str]:
    return list(dict.fromkeys(
        token for token in _TERM.findall(str(phenotype).casefold()) if token not in _STOPWORDS))


def _searchable_reference(reference_gene: str) -> tuple[str, str]:
    """Return ``(query term, how it was derived)`` for a reference-gene identifier.

    The family index writes the BLAST subject verbatim, so a reference gene
    arrives as a UniProt FASTA defline: ``sp|P93003|TFL1_ARATH``. Europe PMC
    finds nothing for that string -- it returns 0 hits where ``TFL1`` returns
    thousands -- so searching on it silently turned every reference into
    ``no_hit`` and left the Route B semantic lane with no citation to verify
    (issue #480). Search on the entry symbol when the defline carries a usable
    one, fall back to the accession, and record which was used so a zero result
    can be read as an answer rather than an artefact of the identifier.
    """
    token = str(reference_gene or "").strip()
    if not token:
        return "", "none"
    if "|" not in token:
        return token, "as_supplied"
    parts = [part for part in token.split("|") if part]
    entry_name = parts[-1] if len(parts) >= 3 else ""
    symbol = entry_name.split("_")[0] if entry_name else ""
    if symbol and is_usable_symbol(symbol) and not is_uniprot_accession(symbol):
        return symbol, "uniprot_defline_symbol"
    accession = canonical_accession(token)
    if accession and accession != token:
        return accession, "uniprot_defline_accession"
    return "", "unusable_identifier"


def _query_for(reference_gene: str, terms: list[str]) -> str:
    focus = " OR ".join(f'"{term}"' for term in terms[:6])
    return f'"{reference_gene}" AND ({focus})' if focus else f'"{reference_gene}"'


def _load_cache(path_arg: str | None) -> dict[str, list[dict[str, object]]]:
    """Query-keyed cached Europe PMC results, so a rerun is byte-identical."""
    if not path_arg:
        return {}
    path = runtime_path(path_arg)
    cache: dict[str, list[dict[str, object]]] = {}
    for item in (sorted(path.glob("*.json")) if path.is_dir() else [path]):
        payload = json.loads(item.read_text(encoding="utf-8"))
        if isinstance(payload, dict) and "query" in payload:
            cache[safe_id(str(payload["query"]))] = [record for record in payload.get("results", [])
                                                     if isinstance(record, dict)]
        elif isinstance(payload, dict):
            for key, value in payload.items():
                if isinstance(value, list):
                    cache[safe_id(str(key))] = [record for record in value if isinstance(record, dict)]
    return cache


# Europe PMC returns an abstract only for resultType=core. The claim
# classifier reads `abstractText`, and its cue lexicon is built from phrases
# that live in abstracts -- "loss-of-function", "complementation", "T-DNA
# insertion". Asking for "lite" therefore handed the classifier titles alone,
# so mechanistic support could not be detected and every claim fell through to
# review_required. The cue lexicon is unchanged; only what gets
# fetched is.
EUROPE_PMC_RESULT_TYPE = "core"


def _fetch(query: str, limit: int, timeout: int = 30) -> list[dict[str, object]]:
    url = f"{EUROPE_PMC_SEARCH}?" + urllib.parse.urlencode(
        {"query": query, "format": "json", "pageSize": limit,
         "resultType": EUROPE_PMC_RESULT_TYPE})
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = json.load(response)
    except urllib.error.HTTPError as error:
        raise RuntimeError(f"Europe PMC HTTP {error.code} for query {query!r}") from error
    except urllib.error.URLError as error:
        raise RuntimeError(f"Europe PMC connection failed for query {query!r}: {error.reason}") from error
    results = payload.get("resultList", {}).get("result", []) if isinstance(payload, dict) else []
    return [record for record in results if isinstance(record, dict)]


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
    """Query the literature for every interval gene's reference homolog."""
    if limit < 1:
        raise ValueError("limit must be at least 1")
    lookup_path = runtime_path(lookup_arg)
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    lookup = pd.read_csv(lookup_path, sep="\t", dtype=str).fillna("")
    missing = sorted({"gene_id", "reference_gene"} - set(lookup.columns))
    if missing:
        raise ValueError(f"interval reverse-lookup table is missing columns: {', '.join(missing)}")

    terms = _phenotype_terms(phenotype)
    cache = _load_cache(cache_arg)
    provider = "europe_pmc_rest" if allow_external_api else ("local_query_cache" if cache else "none")

    by_reference: dict[str, dict[str, object]] = {}
    for row in lookup.to_dict("records"):
        reference = str(row.get("reference_gene", "")).strip()
        record = by_reference.setdefault(reference, {"genes": [], "function": str(row.get("reference_function", ""))})
        record["genes"].append(str(row["gene_id"]))  # type: ignore[union-attr]
        if not record["function"]:
            record["function"] = str(row.get("reference_function", ""))

    rows: list[dict[str, object]] = []
    status_rows: list[dict[str, object]] = []
    for reference in sorted(by_reference):
        genes = sorted(by_reference[reference]["genes"])  # type: ignore[arg-type]
        function = str(by_reference[reference]["function"])
        if not reference:
            status_rows.append({"reference_gene": "", "query_term": "", "identifier_source": "none",
                                "interval_gene_count": len(genes),
                                "query_status": "skipped_no_reference_gene", "article_count": 0,
                                "note": "these interval genes have no characterised relative; they stay an unresolved bucket"})
            continue
        query_term, identifier_source = _searchable_reference(reference)
        if not query_term:
            # Not searched is not "no literature". Keep the two apart so an
            # identifier this stage cannot use never reads as an absence of
            # published work about the reference gene.
            status_rows.append({"reference_gene": reference, "query_term": "",
                                "identifier_source": identifier_source,
                                "interval_gene_count": len(genes),
                                "query_status": "skipped_unusable_reference_identifier", "article_count": 0,
                                "note": "no symbol or accession could be derived from this identifier; "
                                        "no search was performed and absence of literature is not established"})
            continue
        query = _query_for(query_term, terms)
        key = safe_id(query)
        if key in cache:
            articles, query_status, retrieved_at = cache[key], "retrieved", _now()
        elif allow_external_api:
            articles, query_status, retrieved_at = _fetch(query, limit), "retrieved", _now()
        else:
            articles, query_status, retrieved_at = [], "not_evaluated_offline", ""
        if query_status == "retrieved" and not articles:
            query_status = "no_hit"
        status_rows.append({
            "reference_gene": reference, "query_term": query_term,
            "identifier_source": identifier_source, "interval_gene_count": len(genes),
            "query_status": query_status, "article_count": len(articles),
            "note": ("no search was performed; absence of literature is not established"
                     if query_status == "not_evaluated_offline" else ""),
        })
        for rank, article in enumerate(articles[:limit], start=1):
            blob = " ".join(str(article.get(key_, "")) for key_ in ("title", "abstractText")).casefold()
            matched = [term for term in terms if term in blob]
            identifier = str(article.get("pmid", "") or article.get("id", ""))
            rows.append({
                "interval_gene_ids": ";".join(genes), "reference_gene": reference,
                "reference_function": function, "query": query, "query_status": query_status,
                "article_rank": rank, "source_id": f"PMID:{identifier}" if identifier else "",
                "title": re.sub(r"\s+", " ", str(article.get("title", ""))).strip(),
                "journal": str(article.get("journalTitle", "")), "year": str(article.get("pubYear", "")),
                "doi": str(article.get("doi", "")),
                "first_author": str(article.get("authorString", "")).split(",")[0],
                "matched_phenotype_terms": ";".join(matched),
                "abstract_available": str(bool(str(article.get("abstractText", "")).strip())),
                # Retrieval is not relevance: a returned paper is a candidate link.
                "claim_type": "literature_link_candidate",
                "provider": provider, "retrieved_at": retrieved_at, "rank_contribution": "none",
            })

    literature = pd.DataFrame(rows, columns=LITERATURE_COLUMNS)
    status = pd.DataFrame(status_rows, columns=STATUS_COLUMNS).sort_values("reference_gene", kind="stable")
    literature_path = out / "reverse_literature.tsv"
    status_path = out / "reverse_literature_status.tsv"
    literature.to_csv(literature_path, sep="\t", index=False, lineterminator="\n")
    status.to_csv(status_path, sep="\t", index=False, lineterminator="\n")

    manifest = {
        "schema_version": SCHEMA_VERSION, "command": "reverse-literature", "case_id": case_id,
        "candidate_qtl_version": __version__, "queried_at_utc": _now(),
        "route": "interval_first", "stage": "reference_gene_to_literature",
        "provider": provider, "external_api_allowed": allow_external_api,
        # Recorded so a cached run taken before the abstract fix is not mistaken
        # for one that had abstracts to classify.
        "result_type": EUROPE_PMC_RESULT_TYPE,
        "abstract_available_article_count": int(
            literature["abstract_available"].eq("True").sum()
        ) if "abstract_available" in literature.columns and not literature.empty else 0,
        "phenotype_terms": terms,
        "reference_gene_count": int(len(status[status["reference_gene"] != ""])),
        # `queried` counts references that came back with at least one article;
        # `query_attempted` counts the searches actually issued, so a run where
        # every search returned nothing is distinguishable from one where no
        # search was issued at all.
        "queried": int((status["query_status"] == "retrieved").sum()),
        "query_attempted": int(status["query_status"].isin({"retrieved", "no_hit"}).sum()),
        "no_hit": int((status["query_status"] == "no_hit").sum()),
        "skipped_unusable_reference_identifier": int(
            (status["query_status"] == "skipped_unusable_reference_identifier").sum()
        ),
        "identifier_source_counts": {
            str(key): int(value)
            for key, value in status["identifier_source"].value_counts().sort_index().items()
        },
        "not_evaluated_offline": int((status["query_status"] == "not_evaluated_offline").sum()),
        "interval_genes_without_reference": int(status[status["reference_gene"] == ""]["interval_gene_count"].sum())
        if len(status) else 0,
        "article_rows": int(len(literature)),
        "external_transmission": (
            "reference gene names and generic phenotype words only; interval coordinates, candidate ids and "
            "sequences are never transmitted"
        ) if allow_external_api else "none",
        "claims": [
            "a returned article is a link candidate, not evidence: retrieval is not relevance",
            "an offline query is not_evaluated, never 'no literature exists'",
            "interval genes without a characterised relative stay an explicit unresolved bucket",
            "nothing here writes a score",
        ],
        "next_step": "judge each link candidate, then source-verify the citations that are kept",
        "inputs": {"reverse_lookup": file_fingerprint(lookup_path, full_hash)},
        "output_sha256": {"reverse_literature": sha256(literature_path),
                          "reverse_literature_status": sha256(status_path)},
        "outputs": {"literature": str(literature_path), "status": str(status_path),
                    "manifest": str(out / "reverse_literature_manifest.json")},
    }
    manifest_path = out / "reverse_literature_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), "literature": str(literature_path), "status": str(status_path),
            "reference_gene_count": manifest["reference_gene_count"], "article_rows": int(len(literature)),
            "not_evaluated_offline": manifest["not_evaluated_offline"], "provider": provider}
