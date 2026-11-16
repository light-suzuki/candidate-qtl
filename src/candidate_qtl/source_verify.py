"""Verify the literature a proposal cites, not just that it cited something.

``hypothesis-verify`` closed the accession hole: a proposed reference gene is
promoted only after its accession, species and sequence are re-retrieved from an
authoritative database and fingerprinted. Model memory is not a citable source.

The same rule was not applied to the literature. ``mechanism-import`` checks that
``source_ids`` is non-empty and nothing more, so a fabricated ``PMID:1234567``
passes exactly like a real one. That is the one remaining place where a model's
recollection is accepted as evidence, and it is the place an LLM is most likely
to be wrong: plausible-looking citation ids are cheap to generate.

This module resolves every literature citation against Europe PMC and records
what came back:

* a citation that does not resolve is ``rejected`` with
  ``citation_not_found_in_literature_database`` -- the anti-hallucination gate;
* a citation that resolves carries its title, journal, year, DOI and authors as
  provenance, so a later reader can check the claim without re-running anything;
* when the hypothesis names a gene / family / pathway, the retrieved title and
  abstract are searched for those terms. A citation that resolves but never
  mentions the claimed subject is **not** auto-rejected -- whether a paper
  supports a mechanism is a judgement, not a string match -- it is
  ``review_required`` with the terms that were and were not found.

Boundaries:

* Only citation identifiers leave this machine. The phenotype text, the interval,
  the candidate ids and the sequences are never transmitted.
* Offline without a cache is ``review_required``, never ``rejected``: an
  unavailable provider is not evidence that a citation is fake.
* A database accession (UniProt / TAIR / Ensembl) is not a literature citation;
  it is recorded as ``not_literature_source`` and left to ``hypothesis-verify``.
* Nothing here writes a score. Verification decides eligibility, not rank.
"""

from __future__ import annotations

import hashlib
import json
import re
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone

import pandas as pd

from . import __version__
from .claim_evidence import (
    CLAIM_SUPPORT_STATUSES, EXTRACTION_COLUMNS, classifiable_text, classify_claim, gate_extraction, scope_for,
)
from .core import file_fingerprint, runtime_path, safe_id
from .source_review_capture import capture_source

SCHEMA_VERSION = "1.1"
EUROPE_PMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
USER_AGENT = f"candidate-qtl/{__version__} source-verify (read-only literature lookup)"

VALIDATION_COLUMNS = [
    "hypothesis_id", "citation", "citation_type", "verification_status", "resolved_id",
    "title", "journal", "year", "doi", "first_author", "citation_query",
    "source_provider", "source_url", "retrieved_at",
    "claim_terms", "claim_terms_found", "claim_support",
    # issue #226: citation resolution and claim support are two different axes.
    "citation_verified", "claim_support_status", "evidence_relation", "evidence_cues",
    "retrieval_scope", "extraction_status", "quoted_support",
    "rejection_reason", "review_reason",
]
HYPOTHESIS_STATUS_COLUMNS = [
    "hypothesis_id", "citation_count", "literature_citation_count", "verified_citation_count",
    "rejected_citation_count", "review_required_citation_count", "database_accession_count",
    "verified_source_ids", "source_verification_status",
    "mechanistic_support_citation_count", "associative_support_citation_count",
    "unsupported_mention_citation_count", "claim_support_summary",
    "publication_source_provenance", "publication_records_sha256",
]
REVIEW_COLUMNS = ["reason", "target", "required_next_action", "allowed_choices"]

CITATION_TYPES = ("pmid", "pmcid", "doi", "database_accession", "unrecognised")
HYPOTHESIS_STATUSES = ("all_verified", "partly_verified", "none_verified", "no_literature_citation")

# Terms this generic in a title/abstract say nothing about whether the paper is
# about the claimed family, so they never count as claim support.
_GENERIC_CLAIM_TERMS = frozenset(
    "gene genes protein proteins family plant plants mutant mutants pathway process "
    "expression function functional analysis analyses study studies".split()
)
_TERM = re.compile(r"[A-Za-z][A-Za-z0-9_-]{2,}")

_PMID = re.compile(r"^(?:pubmed|pmid|medline|med)[:\s]*(\d{1,9})$", re.IGNORECASE)
_BARE_PMID = re.compile(r"^\d{5,9}$")
_PUBMED_URL = re.compile(r"^https?://pubmed\.ncbi\.nlm\.nih\.gov/(\d{1,9})/?(?:[?#].*)?$", re.IGNORECASE)
_PMC_URL = re.compile(r"^https?://(?:pmc\.ncbi\.nlm\.nih\.gov/articles/|(?:www\.)?ncbi\.nlm\.nih\.gov/pmc/articles/)(PMC\d{4,9})/?(?:[?#].*)?$", re.IGNORECASE)
_PMCID = re.compile(r"^(?:pmcid[:\s]*|pmc:\s*)?(PMC\d{4,9})$", re.IGNORECASE)
_DOI = re.compile(r"^(?:doi[:\s]*|https?://(?:dx\.)?doi\.org/)?(10\.\d{4,9}/\S+)$", re.IGNORECASE)
_ACCESSION_PREFIX = re.compile(
    r"^(uniprot|uniprotkb|tair|ensembl|ensemblplants|refseq|genbank|embl|geneid|ncbi|araport|phytozome|gramene)[:\s]",
    re.IGNORECASE,
)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def classify_citation(value: str) -> tuple[str, str]:
    """Return (citation_type, normalised id) for one citation string."""
    text = unicodedata.normalize("NFKC", str(value)).strip()
    if not text:
        return "unrecognised", ""
    if _ACCESSION_PREFIX.match(text):
        return "database_accession", text
    match = _PMID.match(text) or _PUBMED_URL.match(text)
    if match:
        return "pmid", match.group(1)
    if _BARE_PMID.match(text):
        return "pmid", text
    match = _PMCID.match(text) or _PMC_URL.match(text)
    if match:
        return "pmcid", match.group(1).upper()
    match = _DOI.match(text)
    if match:
        # DOI identifiers are case-insensitive.  Canonicalising here also
        # prevents duplicate live lookups for casing-only variants.
        return "doi", match.group(1).rstrip(".,;").casefold()
    return "unrecognised", text


def _query_for(citation_type: str, identifier: str) -> str:
    if citation_type == "pmid":
        return f"EXT_ID:{identifier} AND SRC:MED"
    if citation_type == "pmcid":
        return f"PMCID:{identifier}"
    return f'DOI:"{identifier}"'


def _cache_key(citation_type: str, identifier: str) -> str:
    token = str(identifier)
    kind, normalized = classify_citation(token)
    if kind != citation_type:
        kind, normalized = classify_citation(f"{citation_type}:{token}")
    if kind == citation_type:
        identifier = normalized
    if citation_type == "doi":
        identifier = str(identifier).casefold()
    elif citation_type == "pmcid":
        identifier = str(identifier).upper()
    return safe_id(f"{citation_type}_{identifier}")


def _load_cache(path_arg: str | None) -> dict[str, dict[str, object]]:
    """Locally cached Europe PMC result objects, keyed by cache key or by id.

    A cache lets the gate run with no network at all, which is the default for a
    benchmark rerun: the same citations must resolve the same way every time.
    """
    if not path_arg:
        return {}
    path = runtime_path(path_arg)
    records: dict[str, dict[str, object]] = {}
    files = sorted(path.glob("*.json")) if path.is_dir() else [path]
    for item in files:
        payload = json.loads(item.read_text(encoding="utf-8"))
        entries = payload if isinstance(payload, list) else [payload]
        if isinstance(payload, dict):
            nested = payload.get("resultList", {}).get("result") if isinstance(payload.get("resultList"), dict) else None
            if isinstance(nested, list):
                entries = nested
            elif isinstance(payload.get("records"), list):
                entries = payload["records"]
        for entry in entries:
            if not isinstance(entry, dict):
                continue
            for key in (str(entry.get("cache_key", "")), _cache_key("pmid", str(entry.get("pmid", "") or entry.get("id", ""))),
                        _cache_key("pmcid", str(entry.get("pmcid", ""))),
                        _cache_key("doi", str(entry.get("doi", "")))):
                if key and not key.endswith("_"):
                    records.setdefault(key, entry)
    return records


# Europe PMC returns an abstract only for resultType=core. The claim
# classifier reads `abstractText`, and its cue lexicon is built from phrases
# that live in abstracts -- "loss-of-function", "complementation", "T-DNA
# insertion". Asking for "lite" therefore handed the classifier titles alone,
# so mechanistic support could not be detected and every claim fell through to
# review_required. The cue lexicon is unchanged; only what gets
# fetched is.
EUROPE_PMC_RESULT_TYPE = "core"
NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _fetch_europe_pmc(query: str, timeout: int = 30) -> dict[str, object] | None:
    url = f"{EUROPE_PMC_SEARCH}?" + urllib.parse.urlencode(
        {"query": query, "format": "json", "pageSize": 1,
         "resultType": EUROPE_PMC_RESULT_TYPE})
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = json.load(response)
    except urllib.error.HTTPError as error:
        raise RuntimeError(f"Europe PMC HTTP {error.code} for query {query!r}") from error
    except urllib.error.URLError as error:
        raise RuntimeError(f"Europe PMC connection failed for query {query!r}: {error.reason}") from error
    except (TimeoutError, OSError) as error:
        raise RuntimeError(f"Europe PMC connection failed for query {query!r}: {error}") from error
    except ValueError as error:
        raise RuntimeError(f"Europe PMC returned invalid JSON for query {query!r}") from error
    if not isinstance(payload, dict):
        raise RuntimeError(f"Europe PMC returned an invalid payload shape for query {query!r}")
    result_list = payload.get("resultList")
    if (not isinstance(result_list, dict) or "result" not in result_list
            or not isinstance(result_list.get("result"), list)):
        raise RuntimeError(f"Europe PMC returned an invalid result list for query {query!r}")
    results = result_list.get("result", [])
    if results and not isinstance(results[0], dict):
        raise RuntimeError(f"Europe PMC returned an invalid result record for query {query!r}")
    return results[0] if results else None


def _ncbi_text(element: ET.Element | None) -> str:
    return re.sub(r"\s+", " ", "".join(element.itertext()) if element is not None else "").strip()


def _fetch_ncbi_eutils(citation_type: str, identifier: str, timeout: int = 30) -> dict[str, object] | None:
    """Retrieve title/abstract metadata from official NCBI E-utilities.

    Europe PMC remains the primary provider.  This function is a deliberately
    small fallback for an already-recognised citation identifier; it never
    searches by phenotype/gene text and therefore transmits only the citation
    id.  ``None`` means a clean no-hit, while transport/HTTP failures raise so
    the caller can record ``review_required`` rather than fabricate absence.
    """
    if citation_type not in {"pmid", "pmcid", "doi"} or not str(identifier).strip():
        return None
    if citation_type == "pmid":
        db, article_id = "pubmed", str(identifier).strip()
        search_url = ""
    else:
        db = "pmc" if citation_type == "pmcid" else "pubmed"
        term = f"{identifier}[accn]" if citation_type == "pmcid" else f'"{identifier}"[doi]'
        search_url = f"{NCBI_EUTILS_BASE}/esearch.fcgi?" + urllib.parse.urlencode(
            {"db": db, "term": term, "retmode": "json", "retmax": 1}
        )
        request = urllib.request.Request(search_url, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                payload = json.load(response)
        except urllib.error.HTTPError as error:
            raise RuntimeError(f"NCBI ESearch HTTP {error.code} for {citation_type}:{identifier}") from error
        except urllib.error.URLError as error:
            raise RuntimeError(f"NCBI ESearch connection failed for {citation_type}:{identifier}: {error.reason}") from error
        except (TimeoutError, OSError) as error:
            raise RuntimeError(f"NCBI ESearch connection failed for {citation_type}:{identifier}: {error}") from error
        except ValueError as error:
            raise RuntimeError(f"NCBI ESearch returned invalid JSON for {citation_type}:{identifier}") from error
        if not isinstance(payload, dict):
            raise RuntimeError(f"NCBI ESearch returned an invalid payload shape for {citation_type}:{identifier}")
        search_result = payload.get("esearchresult")
        if (not isinstance(search_result, dict) or "idlist" not in search_result
                or not isinstance(search_result.get("idlist"), list)):
            raise RuntimeError(f"NCBI ESearch returned an invalid result list for {citation_type}:{identifier}")
        ids = search_result.get("idlist", [])
        if not ids:
            return None
        article_id = str(ids[0])

    fetch_url = f"{NCBI_EUTILS_BASE}/efetch.fcgi?" + urllib.parse.urlencode(
        {"db": db, "id": article_id, "rettype": "abstract", "retmode": "xml"}
    )
    request = urllib.request.Request(fetch_url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            root = ET.fromstring(response.read())
    except urllib.error.HTTPError as error:
        raise RuntimeError(f"NCBI EFetch HTTP {error.code} for {citation_type}:{identifier}") from error
    except urllib.error.URLError as error:
        raise RuntimeError(f"NCBI EFetch connection failed for {citation_type}:{identifier}: {error.reason}") from error
    except (TimeoutError, OSError) as error:
        raise RuntimeError(f"NCBI EFetch connection failed for {citation_type}:{identifier}: {error}") from error
    except ET.ParseError as error:
        raise RuntimeError(f"NCBI EFetch returned invalid XML for {citation_type}:{identifier}") from error

    article = next(
        (node for node in root.iter()
         if str(node.tag).rsplit("}", 1)[-1] in {"PubmedArticle", "article"}),
        None,
    )
    if article is None:
        return None
    def descendants(*names: str) -> list[ET.Element]:
        wanted = {name.casefold() for name in names}
        return [node for node in article.iter()
                if str(node.tag).rsplit("}", 1)[-1].casefold() in wanted]

    title_nodes = descendants("ArticleTitle", "article-title")
    title = _ncbi_text(title_nodes[0]) if title_nodes else ""
    abstract_nodes = descendants("AbstractText", "abstract")
    abstract = " ".join(text for text in (_ncbi_text(node) for node in abstract_nodes) if text)
    journal_nodes = descendants("Title", "journal-title")
    journal = _ncbi_text(journal_nodes[0]) if journal_nodes else ""
    year_nodes = descendants("Year", "year")
    medline_nodes = descendants("MedlineDate")
    year = (_ncbi_text(year_nodes[0]) if year_nodes else "") or (_ncbi_text(medline_nodes[0]) if medline_nodes else "")
    authors = []
    for node in descendants("Author"):
        name_nodes = [child for child in node.iter()
                      if str(child.tag).rsplit("}", 1)[-1].casefold() in {"lastname", "surname", "collectivename"}]
        name = _ncbi_text(name_nodes[0]) if name_nodes else ""
        if name:
            authors.append(name)
    doi = ""
    id_nodes = descendants("ArticleId", "article-id")
    pmid = next((_ncbi_text(node) for node in id_nodes if str(node.attrib.get("IdType", node.attrib.get("pub-id-type", ""))).casefold() in {"pubmed", "pmid"}), article_id)
    pmcid = next((_ncbi_text(node) for node in id_nodes if str(node.attrib.get("IdType", node.attrib.get("pub-id-type", ""))).casefold() in {"pmc", "pmcid"}), "")
    for node in id_nodes:
        if str(node.attrib.get("IdType", node.attrib.get("pub-id-type", ""))).casefold() == "doi":
            doi = _ncbi_text(node)
            break
    return {
        "id": pmid, "pmid": pmid, "pmcid": pmcid, "title": title,
        "abstractText": abstract, "journalTitle": journal, "pubYear": year,
        "doi": doi, "authorString": ", ".join(authors),
        "_source_provider": "ncbi_eutils_pmc" if db == "pmc" else "ncbi_eutils_pubmed",
        "_source_url": fetch_url if not search_url else f"{search_url} -> {fetch_url}",
    }


def _claim_terms(row: dict[str, object]) -> list[str]:
    """Distinctive terms the citation is supposed to be about."""
    terms: list[str] = []
    for column in ("reference_gene_name", "gene_family", "pathway_or_process"):
        for token in _TERM.findall(str(row.get(column, ""))):
            if token.casefold() not in _GENERIC_CLAIM_TERMS:
                terms.append(token)
    return list(dict.fromkeys(terms))


def _claim_text(record: dict[str, object]) -> str:
    return " ".join(str(record.get(key, "")) for key in
                    ("title", "abstractText", "journalTitle", "keywordList")).casefold()


def _claim_support(
    record: dict[str, object], terms: list[str], *, normalised_text: str | None = None,
) -> tuple[list[str], str]:
    if not terms:
        return [], "not_checked_no_claim_term"
    blob = _claim_text(record) if normalised_text is None else normalised_text
    if not blob.strip():
        return [], "not_checked_no_retrievable_text"
    found = [term for term in terms if term.casefold() in blob]
    return found, "supported" if found else "unsupported"


def _record_fields(record: dict[str, object]) -> dict[str, str]:
    return {
        "resolved_id": str(record.get("id", "") or record.get("pmid", "")),
        "title": re.sub(r"\s+", " ", str(record.get("title", ""))).strip(),
        "journal": str(record.get("journalTitle", "")),
        "year": str(record.get("pubYear", "")),
        "doi": str(record.get("doi", "")),
        "first_author": str(record.get("authorString", "")).split(",")[0],
    }


BENCHMARK_CASE_COLUMNS = [
    "case_id", "citation", "citation_type", "expected_claim_support_status", "claim_support_status",
    "evidence_relation", "evidence_cues", "retrieval_scope", "citation_verified",
    "expected_citation_verified", "outcome", "not_evaluated_reason",
]
BENCHMARK_SUMMARY_COLUMNS = ["expected_claim_support_status", "cases", "evaluated", "correct",
                             "accuracy", "not_evaluated"]


def claim_verification_benchmark(
    cases_arg: str,
    out_dir_arg: str,
    cache_arg: str | None = None,
    allow_external_api: bool = False,
    full_hash: bool = False,
    cache_is_complete: bool = False,
) -> dict[str, object]:
    """Score the claim-support classifier against labelled examples.

    The four labels the design cares about are genuinely different failures:
    real mechanistic support, support that is only an association, a real paper
    that merely mentions the gene, and a citation that does not exist at all.
    A case whose record could not be retrieved is ``not_evaluated`` and leaves
    the accuracy denominator -- an unavailable provider is not a wrong answer.
    """
    cases_path = runtime_path(cases_arg)
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(cases_path, sep="\t", dtype=str).fillna("")
    missing = [column for column in ("case_id", "citation", "expected_claim_support_status")
               if column not in frame.columns]
    if missing:
        raise ValueError(f"claim benchmark cases are missing columns: {', '.join(missing)}")
    if frame["case_id"].duplicated().any():
        raise ValueError("claim benchmark case_id must be unique")

    cache = _load_cache(cache_arg)
    rows: list[dict[str, object]] = []
    providers_used: set[str] = set()
    live_lookup: dict[str, tuple[dict[str, object] | None, str, str]] = {}
    for record in frame.to_dict("records"):
        citation = str(record["citation"]).strip()
        expected = str(record["expected_claim_support_status"]).strip()
        expected_verified = str(record.get("expected_citation_verified", "")).strip()
        citation_type, identifier = classify_citation(citation)
        terms = _claim_terms(record)
        row: dict[str, object] = {
            "case_id": str(record["case_id"]).strip(), "citation": citation, "citation_type": citation_type,
            "expected_claim_support_status": expected, "claim_support_status": "", "evidence_relation": "",
            "evidence_cues": "", "retrieval_scope": "none", "citation_verified": "",
            "expected_citation_verified": expected_verified, "outcome": "", "not_evaluated_reason": "",
        }
        if citation_type in {"database_accession", "unrecognised"}:
            # Not a literature claim at all: the citation axis already answers it.
            row.update({"citation_verified": False, "claim_support_status": "not_evaluated",
                        "outcome": "not_evaluated",
                        "not_evaluated_reason": f"citation_type_{citation_type}_is_not_a_literature_claim"})
            rows.append(row)
            continue
        key = _cache_key(citation_type, identifier)
        result = cache.get(key)
        result_provider = ""
        provider_error = ""
        if result is not None:
            result_provider = "local_record_cache"
        if result is None and allow_external_api:
            if key not in live_lookup:
                attempted = ["europe_pmc_rest"]
                try:
                    result = _fetch_europe_pmc(_query_for(citation_type, identifier))
                except RuntimeError as error:
                    provider_error = str(error)
                if result is not None:
                    result_provider = "europe_pmc_rest"
                else:
                    attempted.append("ncbi_eutils_pubmed")
                    try:
                        result = _fetch_ncbi_eutils(citation_type, identifier)
                    except RuntimeError as error:
                        provider_error = f"{provider_error}; {error}".strip("; ")
                    else:
                        provider_error = ""
                        if result is not None:
                            result_provider = str(result.get("_source_provider", "ncbi_eutils_pubmed"))
                live_lookup[key] = (result, provider_error, result_provider or "+".join(attempted))
            result, provider_error, result_provider = live_lookup[key]
        if result_provider:
            providers_used.update(result_provider.split("+"))
        if result is None and not (allow_external_api or cache_is_complete):
            row.update({"outcome": "not_evaluated",
                        "not_evaluated_reason": "record_not_in_cache_and_no_provider_was_allowed"})
            rows.append(row)
            continue
        if result is None and provider_error:
            row.update({"outcome": "not_evaluated", "not_evaluated_reason": "literature_provider_unavailable"})
            rows.append(row)
            continue
        if result is None:
            # A declared-complete cache (or a live lookup) that returns nothing is
            # the fabricated-citation case, and it is evaluable.
            row.update({"citation_verified": False, "claim_support_status": "not_evaluated",
                        "evidence_relation": "not_evaluated_no_text",
                        "outcome": "correct" if expected in ("", "not_evaluated") or expected_verified == "False"
                        else "incorrect"})
            row["not_evaluated_reason"] = "citation_did_not_resolve"
            rows.append(row)
            continue
        scope = scope_for(result)
        found, _support = _claim_support(result, terms)
        relation, status, cues = classify_claim(classifiable_text(result), found, terms, scope)
        row.update({"citation_verified": True, "claim_support_status": status, "evidence_relation": relation,
                    "evidence_cues": ";".join(cues), "retrieval_scope": scope,
                    "outcome": "correct" if status == expected else "incorrect"})
        rows.append(row)

    cases = pd.DataFrame(rows, columns=BENCHMARK_CASE_COLUMNS).sort_values("case_id", kind="stable")
    cases_out = out / "claim_benchmark_cases.tsv"
    cases.to_csv(cases_out, sep="\t", index=False, lineterminator="\n")

    summary_rows: list[dict[str, object]] = []
    for expected in sorted({str(value) for value in cases["expected_claim_support_status"]}):
        subset = cases[cases["expected_claim_support_status"] == expected]
        evaluated = int((subset["outcome"] != "not_evaluated").sum())
        correct = int((subset["outcome"] == "correct").sum())
        summary_rows.append({"expected_claim_support_status": expected, "cases": int(len(subset)),
                             "evaluated": evaluated, "correct": correct,
                             "accuracy": round(correct / evaluated, 6) if evaluated else "",
                             "not_evaluated": int(len(subset)) - evaluated})
    summary = pd.DataFrame(summary_rows, columns=BENCHMARK_SUMMARY_COLUMNS)
    summary_out = out / "claim_benchmark_summary.tsv"
    summary.to_csv(summary_out, sep="\t", index=False, lineterminator="\n")
    confusion = (cases.groupby(["expected_claim_support_status", "claim_support_status"]).size()
                 .reset_index(name="cases") if len(cases)
                 else pd.DataFrame(columns=["expected_claim_support_status", "claim_support_status", "cases"]))
    confusion_out = out / "claim_benchmark_confusion.tsv"
    confusion.to_csv(confusion_out, sep="\t", index=False, lineterminator="\n")

    evaluated_total = int((cases["outcome"] != "not_evaluated").sum()) if len(cases) else 0
    correct_total = int((cases["outcome"] == "correct").sum()) if len(cases) else 0
    manifest = {
        "schema_version": SCHEMA_VERSION, "command": "claim-verification-benchmark",
        "candidate_qtl_version": __version__, "run_at_utc": _now(),
        "case_count": int(len(cases)), "evaluated": evaluated_total, "correct": correct_total,
        "accuracy": round(correct_total / evaluated_total, 6) if evaluated_total else None,
        "not_evaluated": int(len(cases)) - evaluated_total,
        "provider": "+".join(sorted(providers_used)) if providers_used else (
            "europe_pmc_rest+ncbi_eutils_pubmed" if allow_external_api else ("local_record_cache" if cache else "none")
        ),
        "providers_used": sorted(providers_used),
        "cache_is_complete": cache_is_complete,
        "denominator_policy": (
            "a case whose record could not be retrieved is not_evaluated and leaves the denominator; an "
            "unavailable provider is not a wrong classification. A citation that genuinely does not resolve is "
            "only evaluable with a live lookup or an explicitly complete cache"
        ),
        "labels": list(CLAIM_SUPPORT_STATUSES),
        "inputs": {"cases": file_fingerprint(cases_path, full_hash)},
        "outputs": {"cases": str(cases_out), "summary": str(summary_out), "confusion": str(confusion_out)},
    }
    manifest_path = out / "claim_benchmark_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), "cases": str(cases_out), "summary": str(summary_out),
            "confusion": str(confusion_out), "case_count": manifest["case_count"],
            "evaluated": evaluated_total, "correct": correct_total, "accuracy": manifest["accuracy"],
            "not_evaluated": manifest["not_evaluated"]}


def _support_bucket(status: str) -> str | None:
    """Which support_counts bucket a claim_support_status increments (or none)."""
    if status in ("supported_direct_genetic_test", "supported_perturbation_phenotype"):
        return "mechanistic"
    if status == "associative_support_only":
        return "associative"
    if status == "unsupported_mention":
        return "unsupported_mention"
    return None


def _claim_summary(support_counts: dict[str, int], support_states: list[str]) -> str:
    """Roll per-citation evidence judgements up to one hypothesis-level summary."""
    if support_counts["mechanistic"]:
        return "mechanistic_support"
    if support_counts["associative"]:
        return "associative_support_only"
    if any(state == "review_required" for state in support_states):
        return "review_required"
    if support_counts["unsupported_mention"]:
        return "unsupported_mention_only"
    if support_states:
        return "background_or_review_only"
    return "not_evaluated"


def _extraction_rows(path_arg: str | None) -> list[dict[str, object]]:
    """Agent/LLM claim extractions, read through one schema."""
    if not path_arg:
        return []
    path = runtime_path(path_arg)
    if path.suffix.casefold() == ".tsv":
        frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        missing = [column for column in ("citation", "evidence_relation", "quoted_support")
                   if column not in frame.columns]
        if missing:
            raise ValueError(f"claim extraction {path} is missing columns: {', '.join(missing)}")
        return frame.to_dict("records")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(payload, dict):
        for key in ("records", "rows", "claims", "extractions"):
            if isinstance(payload.get(key), list):
                return [item for item in payload[key] if isinstance(item, dict)]
        return [payload]
    return [item for item in payload if isinstance(item, dict)]


def source_verify(
    mechanisms_arg: str, out_dir_arg: str, cache_arg: str | None = None,
    allow_external_api: bool = False, id_column: str = "hypothesis_id",
    source_column: str = "source_ids", full_hash: bool = False,
    claim_extraction_arg: str | None = None, verification_policy: str = "content_v1",
    focal_species: str = "", reference_bridge_arg: str | None = None,
) -> dict[str, object]:
    if reference_bridge_arg:
        from .reference_sequence_bridge import load_reference_sequence_bridge
        load_reference_sequence_bridge(reference_bridge_arg)
    result = _source_verify(mechanisms_arg, out_dir_arg, cache_arg, allow_external_api,
                           id_column, source_column, full_hash, claim_extraction_arg,
                           verification_policy, focal_species)
    if reference_bridge_arg:
        from .reference_bridge_consumers import attach_paper_sequence_correspondence
        result = attach_paper_sequence_correspondence(result, mechanisms_arg, reference_bridge_arg, id_column=id_column)
    return result


def _source_verify(
    mechanisms_arg: str,
    out_dir_arg: str,
    cache_arg: str | None = None,
    allow_external_api: bool = False,
    id_column: str = "hypothesis_id",
    source_column: str = "source_ids",
    full_hash: bool = False,
    claim_extraction_arg: str | None = None,
    verification_policy: str = "content_v1",
    focal_species: str = "",
) -> dict[str, object]:
    """Resolve every literature citation before it may count as a source."""
    from .citation_tokens import split_citations
    from .bibliography import CONTENT_POLICY, INITIAL_POLICY, verify_bibliography
    if verification_policy == INITIAL_POLICY:
        if claim_extraction_arg:
            raise ValueError("bibliography-only verification cannot accept content extraction")
        return verify_bibliography(mechanisms_arg, out_dir_arg, cache_arg=cache_arg,
            allow_external_api=allow_external_api, id_column=id_column, source_column=source_column,
            focal_species=focal_species)
    if verification_policy != CONTENT_POLICY:
        raise ValueError("unsupported source verification policy")
    mechanisms_path = runtime_path(mechanisms_arg)
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(mechanisms_path, sep="\t", dtype=str).fillna("")
    missing = sorted({id_column, source_column} - set(frame.columns))
    if missing:
        raise ValueError(f"hypothesis table is missing columns: {', '.join(missing)}")
    if frame[id_column].duplicated().any():
        raise ValueError(f"hypothesis table requires unique {id_column}")

    cache = _load_cache(cache_arg)
    provider = "europe_pmc_rest" if allow_external_api else ("local_record_cache" if cache else "none")
    providers_used: set[str] = set()
    rows: list[dict[str, object]] = []
    hypothesis_rows: list[dict[str, object]] = []
    review_rows: list[dict[str, str]] = []
    verified_records: dict[str, dict[str, object]] = {}
    review_source_records: dict[str, dict[str, object]] = {}
    retrieved_text: dict[str, str] = {}
    # Provider/cache records are immutable within this invocation. Reuse their
    # text preparation and cue scans, but still check each hypothesis's terms.
    source_texts: dict[str, tuple[str, str]] = {}
    cue_cache: dict[str, tuple[str, str, tuple[str, ...]]] = {}
    # Per-hypothesis aggregates stay reachable after the citation loop so an
    # accepted agent extraction can re-open and update the frozen roll-up
    # instead of annotating only the validation rows.
    hypothesis_aggregates: dict[str, dict[str, object]] = {}
    # One frozen claim set often cites the same paper for several interval
    # members.  Resolve each citation identifier once per command invocation;
    # repeated claims reuse the same immutable provider result and provenance.
    live_lookup: dict[str, tuple[dict[str, object] | None, str, str, str]] = {}

    for record in frame.to_dict("records"):
        hypothesis_id = str(record[id_column])
        citations = split_citations(record[source_column])
        terms = _claim_terms(record)
        counts = {"verified": 0, "rejected": 0, "review_required": 0, "not_literature_source": 0}
        support_counts = {"mechanistic": 0, "associative": 0, "unsupported_mention": 0}
        support_states: list[str] = []
        verified_ids: list[str] = []
        verified_provenance: list[dict[str, str]] = []

        for citation in citations:
            citation_type, identifier = classify_citation(citation)
            row: dict[str, object] = {
                "hypothesis_id": hypothesis_id, "citation": citation, "citation_type": citation_type,
                "resolved_id": "", "title": "", "journal": "", "year": "", "doi": "", "first_author": "",
                "citation_query": "", "source_provider": "", "source_url": "", "retrieved_at": "",
                "claim_terms": ";".join(terms), "claim_terms_found": "", "claim_support": "not_checked",
                "citation_verified": False, "claim_support_status": "not_evaluated",
                "evidence_relation": "", "evidence_cues": "", "retrieval_scope": "none",
                "extraction_status": "", "quoted_support": "",
                "rejection_reason": "", "review_reason": "",
            }
            if citation_type == "database_accession":
                # Not a literature claim: hypothesis-verify owns accession identity.
                row.update({"verification_status": "not_literature_source",
                            "review_reason": "database accession; verified by hypothesis-verify, not here"})
                counts["not_literature_source"] += 1
                rows.append(row)
                continue
            if citation_type == "unrecognised":
                row.update({"verification_status": "rejected",
                            "rejection_reason": "citation_is_not_a_recognisable_pmid_pmcid_or_doi"})
                counts["rejected"] += 1
                review_rows.append({"reason": "citation_unparsable", "target": f"{hypothesis_id}:{citation}",
                                    "required_next_action": "supply a PMID, PMCID or DOI, or drop the claim",
                                    "allowed_choices": "fix the citation | drop the claim | keep as unsourced"})
                rows.append(row)
                continue

            query = _query_for(citation_type, identifier)
            row["citation_query"] = query
            key = _cache_key(citation_type, identifier)
            result: dict[str, object] | None = None
            provider_error = ""
            if key in cache:
                result = cache[key]
                providers_used.add("local_record_cache")
                row.update({"source_provider": "local_record_cache",
                            "source_url": str(runtime_path(cache_arg)) if cache_arg else "",
                            "retrieved_at": str(result.get("candidate_qtl_retrieved_at", "")) or _now()})
            elif allow_external_api:
                if key not in live_lookup:
                    # Europe PMC is primary. A clean no-hit and provider
                    # failure remain distinct before the NCBI fallback.
                    attempted = ["europe_pmc_rest"]
                    live_provider = ""
                    live_url = ""
                    try:
                        result = _fetch_europe_pmc(query)
                    except RuntimeError as error:
                        provider_error = str(error)
                    if result is not None:
                        live_provider = "europe_pmc_rest"
                        live_url = f"{EUROPE_PMC_SEARCH}?query={urllib.parse.quote(query)}&format=json"
                    else:
                        attempted.append("ncbi_eutils_pubmed")
                        try:
                            result = _fetch_ncbi_eutils(citation_type, identifier)
                        except RuntimeError as error:
                            provider_error = f"{provider_error}; {error}".strip("; ")
                        else:
                            # A clean NCBI response (including a clean no-hit)
                            # completes the fallback investigation.
                            provider_error = ""
                        if result is not None:
                            live_provider = str(result.get("_source_provider", "ncbi_eutils_pubmed"))
                            live_url = str(result.get("_source_url", ""))
                    # A clean no-hit is still a completed investigation. Keep
                    # every provider that was consulted so downstream semantic
                    # review can distinguish it from an offline/not-run lane.
                    if not live_provider:
                        live_provider = "+".join(attempted)
                    live_lookup[key] = (result, provider_error, live_provider, live_url)
                result, provider_error, live_provider, live_url = live_lookup[key]
                if live_provider:
                    providers_used.update(live_provider.split("+"))
                    row.update({"source_provider": live_provider, "source_url": live_url, "retrieved_at": _now()})
            else:
                # No provider is not evidence that the citation is fake.
                row.update({"verification_status": "review_required",
                            "review_reason": "verification_not_performed_offline"})
                counts["review_required"] += 1
                review_rows.append({"reason": "citation_not_verified_offline", "target": f"{hypothesis_id}:{citation}",
                                    "required_next_action": "supply --cache, or allow the Europe PMC lookup with --allow-external-api",
                                    "allowed_choices": "supply a cache | allow the API | keep unverified"})
                rows.append(row)
                continue

            if result is None and provider_error:
                # External access was explicitly allowed, but the providers
                # were unavailable.  This is not evidence that the citation is
                # fabricated and must not be counted as a rejection.
                row.update({"verification_status": "review_required",
                            "review_reason": "literature_provider_unavailable"})
                counts["review_required"] += 1
                review_rows.append({"reason": "citation_provider_unavailable", "target": f"{hypothesis_id}:{citation}",
                                    "required_next_action": "retry with a working Europe PMC/NCBI provider or supply --cache",
                                    "allowed_choices": "retry the API | supply a cache | keep unverified"})
                rows.append(row)
                continue
            if result is None:
                row.update({"verification_status": "rejected",
                            "rejection_reason": "citation_not_found_in_literature_database"})
                counts["rejected"] += 1
                review_rows.append({"reason": "citation_not_found", "target": f"{hypothesis_id}:{citation}",
                                    "required_next_action": "the citation did not resolve; treat the claim as unsourced until a real reference is supplied",
                                    "allowed_choices": "replace the citation | drop the claim | keep as unsourced"})
                rows.append(row)
                continue

            row.update(_record_fields(result))
            if key not in source_texts:
                source_texts[key] = _claim_text(result), classifiable_text(result)
            claim_text, source_text = source_texts[key]
            found, support = _claim_support(result, terms, normalised_text=claim_text)
            # Two axes, kept apart on purpose: the identifier resolved (a citation
            # fact) and what the retrieved text supports (an evidence judgement).
            scope = scope_for(result)
            if citation not in review_source_records:
                review_source_records[citation] = capture_source(
                    result, citation=citation, citation_type=citation_type,
                    resolved_id=str(row["resolved_id"]), classification_scope=scope,
                )
            captured = review_source_records[citation]
            if hypothesis_id not in captured["hypothesis_ids"]:
                captured["hypothesis_ids"].append(hypothesis_id)
            relation, support_status, cues = classify_claim(
                source_text, found, terms, scope, cue_cache=cue_cache)
            retrieved_text[citation] = source_text
            row.update({"claim_terms_found": ";".join(found), "claim_support": support,
                        "citation_verified": True, "claim_support_status": support_status,
                        "evidence_relation": relation, "evidence_cues": ";".join(cues),
                        "retrieval_scope": scope})
            support_states.append(support_status)
            if support_status in ("supported_direct_genetic_test", "supported_perturbation_phenotype"):
                support_counts["mechanistic"] += 1
            elif support_status == "associative_support_only":
                support_counts["associative"] += 1
            elif support_status == "unsupported_mention":
                support_counts["unsupported_mention"] += 1
            if support_status in ("associative_support_only", "background_or_review_only"):
                # Resolvable, on-topic, and weaker than the mechanistic claim it
                # was cited for. Surfaced for a reader rather than scored.
                review_rows.append({"reason": f"claim_support_is_{support_status}",
                                    "target": f"{hypothesis_id}:{citation}",
                                    "required_next_action": "the source supports an association or repeats background; "
                                                            "cite the primary genetic test or weaken the claim",
                                    "allowed_choices": "cite the primary test | weaken the claim | keep as associative"})
            if support == "unsupported":
                # It resolves, but never mentions what it was cited for. Whether a
                # paper supports a mechanism is a judgement, so this is a review,
                # not a rejection.
                row["verification_status"] = "review_required"
                row["review_reason"] = "citation_resolves_but_does_not_mention_the_claimed_subject"
                counts["review_required"] += 1
                review_rows.append({"reason": "citation_does_not_mention_claim", "target": f"{hypothesis_id}:{citation}",
                                    "required_next_action": "read the reference and confirm it supports this mechanism, or replace it",
                                    "allowed_choices": "confirm after reading | replace the citation | drop the claim"})
            else:
                row["verification_status"] = "verified"
                counts["verified"] += 1
                verified_ids.append(citation)
                provenance = {
                    "citation": citation, "resolved_id": str(row.get("resolved_id", "")),
                    "provider": str(row.get("source_provider", "")),
                    "url": str(row.get("source_url", "")),
                    "retrieved_at": str(row.get("retrieved_at", "")),
                    "claim_support_status": support_status,
                }
                source_record = verified_records.setdefault(
                    citation,
                    {**_record_fields(result), "citation": citation,
                     "citation_type": citation_type, "query": query,
                     "candidate_qtl_hypothesis_verifications": {}},
                )
                verifications = source_record["candidate_qtl_hypothesis_verifications"]
                if isinstance(verifications, dict):
                    verifications[hypothesis_id] = provenance
                verified_provenance.append(provenance)
            rows.append(row)

        literature_total = len(citations) - counts["not_literature_source"]
        if not literature_total:
            status = "no_literature_citation"
            review_rows.append({"reason": "hypothesis_without_literature_citation", "target": hypothesis_id,
                                "required_next_action": "a mechanism with no literature citation cannot be called source-backed",
                                "allowed_choices": "supply a citation | keep as an unsourced hypothesis | drop"})
        elif counts["verified"] == literature_total:
            status = "all_verified"
        elif counts["verified"]:
            status = "partly_verified"
        else:
            status = "none_verified"
            review_rows.append({"reason": "hypothesis_citations_all_unresolved", "target": hypothesis_id,
                                "required_next_action": "no citation for this mechanism resolved; it is not source-backed",
                                "allowed_choices": "supply real citations | keep as an unsourced hypothesis | drop"})
        claim_summary = _claim_summary(support_counts, support_states)
        hypothesis_aggregates[hypothesis_id] = {
            "support_counts": support_counts, "support_states": support_states,
            "verified_provenance": verified_provenance,
        }
        hypothesis_rows.append({
            "hypothesis_id": hypothesis_id, "citation_count": len(citations),
            "literature_citation_count": literature_total, "verified_citation_count": counts["verified"],
            "rejected_citation_count": counts["rejected"],
            "review_required_citation_count": counts["review_required"],
            "database_accession_count": counts["not_literature_source"],
            "verified_source_ids": ";".join(verified_ids), "source_verification_status": status,
            "mechanistic_support_citation_count": support_counts["mechanistic"],
            "associative_support_citation_count": support_counts["associative"],
            "unsupported_mention_citation_count": support_counts["unsupported_mention"],
            "claim_support_summary": claim_summary,
            "publication_source_provenance": json.dumps(
                verified_provenance, ensure_ascii=False, sort_keys=True, separators=(",", ":")),
            # Filled after every hypothesis has been evaluated. A citation record
            # can be shared by several hypotheses, so only the final frozen record
            # graph has a stable digest.
            "publication_records_sha256": "",
        })

    # An agent may extract what a paper says; it may not assert it. Every quoted
    # span must occur in the text that was actually retrieved, or the row is
    # rejected with its reason and annotates nothing.
    hypothesis_citations = {str(row["hypothesis_id"]): set() for row in hypothesis_rows}
    for row in rows:
        hypothesis_citations[str(row["hypothesis_id"])].add(str(row["citation"]))
    extraction_rows = gate_extraction(
        _extraction_rows(claim_extraction_arg), retrieved_text, hypothesis_citations,
    )
    extraction_by_target = {
        (str(row.get("hypothesis_id", "")), str(row["citation"])): row
        for row in extraction_rows if row["extraction_status"] == "accepted"
    }
    for row in rows:
        extracted = extraction_by_target.get((str(row["hypothesis_id"]), str(row["citation"])))
        if extracted is None:
            extracted = extraction_by_target.get(("", str(row["citation"])))
        if not extracted:
            continue
        row["extraction_status"] = "accepted"
        row["quoted_support"] = extracted["quoted_support"]
        # The extraction may sharpen the relation the cue lexicon guessed, but it
        # never turns an unresolved citation into a supported claim.
        if not row.get("citation_verified"):
            continue
        new_status = str(extracted["claim_support_status"])
        row["evidence_relation"] = extracted["evidence_relation"]
        old_status = str(row.get("claim_support_status", ""))
        if old_status == new_status:
            continue
        row["claim_support_status"] = new_status
        # Propagate the accepted-extraction judgement into the frozen hypothesis
        # aggregates.  The provenance dicts are the same objects stored inside
        # verified_records, so updating them here keeps both the TSV roll-up and
        # the frozen record graph (and therefore its sha256 digest) consistent.
        target = f"{row['hypothesis_id']}:{row['citation']}"
        review_rows[:] = [entry for entry in review_rows if not (
            entry["target"] == target and entry["reason"].startswith("claim_support_is_"))]
        aggregate = hypothesis_aggregates.get(str(row["hypothesis_id"]))
        if not aggregate:
            continue
        counts_by_bucket: dict[str, int] = aggregate["support_counts"]  # type: ignore[assignment]
        states: list[str] = aggregate["support_states"]  # type: ignore[assignment]
        old_bucket = _support_bucket(old_status)
        new_bucket = _support_bucket(new_status)
        if old_bucket:
            counts_by_bucket[old_bucket] -= 1
        if new_bucket:
            counts_by_bucket[new_bucket] += 1
        if old_status in states:
            states[states.index(old_status)] = new_status
        else:
            states.append(new_status)
        for provenance in aggregate["verified_provenance"]:  # type: ignore[union-attr]
            if provenance.get("citation") == row["citation"] and provenance.get("claim_support_status") == old_status:
                provenance["claim_support_status"] = new_status
                break
        if new_bucket == "associative":
            review_rows.append({"reason": "claim_support_is_associative_support_only",
                                "target": target,
                                "required_next_action": "the source supports an association or repeats background; "
                                                        "cite the primary genetic test or weaken the claim",
                                "allowed_choices": "cite the primary test | weaken the claim | keep as associative"})
    for row in extraction_rows:
        if row["extraction_status"] != "accepted":
            review_rows.append({"reason": f"claim_extraction_{row['extraction_status']}",
                                "target": str(row["citation"]) or "(missing citation)",
                                "required_next_action": "the extracted claim could not be located in the retrieved "
                                                        "record; re-extract from the real text or drop the claim",
                                "allowed_choices": "re-extract | drop the claim | keep as unverified"})

    # Accepted extractions may have moved per-citation judgements; re-open the
    # frozen roll-ups so the hypothesis table, the provenance it publishes and
    # the digest computed below all describe the same final state.
    for hypothesis_row in hypothesis_rows:
        aggregate = hypothesis_aggregates.get(str(hypothesis_row["hypothesis_id"]))
        if not aggregate:
            continue
        counts_by_bucket: dict[str, int] = aggregate["support_counts"]  # type: ignore[assignment]
        states: list[str] = aggregate["support_states"]  # type: ignore[assignment]
        hypothesis_row.update({
            "mechanistic_support_citation_count": counts_by_bucket["mechanistic"],
            "associative_support_citation_count": counts_by_bucket["associative"],
            "unsupported_mention_citation_count": counts_by_bucket["unsupported_mention"],
            "claim_support_summary": _claim_summary(counts_by_bucket, states),
            "publication_source_provenance": json.dumps(
                aggregate["verified_provenance"], ensure_ascii=False, sort_keys=True, separators=(",", ":")),
        })

    validation = pd.DataFrame(rows, columns=VALIDATION_COLUMNS) if rows else pd.DataFrame(columns=VALIDATION_COLUMNS)
    extraction = (pd.DataFrame(extraction_rows) if extraction_rows
                  else pd.DataFrame(columns=[*EXTRACTION_COLUMNS, "extraction_status", "claim_support_status"]))
    # The graph is now final, including accepted extractions. Shared citation
    # sets have identical digests even when their hypothesis-specific terms and
    # provenance differ. Preserve duplicate citations in the original hash input.
    publication_digests: dict[tuple[str, ...], str] = {}
    for hypothesis_row in hypothesis_rows:
        citation_ids = tuple(sorted(split_citations(hypothesis_row["verified_source_ids"])))
        if citation_ids:
            if citation_ids not in publication_digests:
                publication_digests[citation_ids] = hashlib.sha256(json.dumps(
                    [verified_records[value] for value in citation_ids],
                    ensure_ascii=False, sort_keys=True, separators=(",", ":"),
                ).encode("utf-8")).hexdigest()
            hypothesis_row["publication_records_sha256"] = publication_digests[citation_ids]

    hypothesis_status = pd.DataFrame(hypothesis_rows, columns=HYPOTHESIS_STATUS_COLUMNS) if hypothesis_rows \
        else pd.DataFrame(columns=HYPOTHESIS_STATUS_COLUMNS)
    reviews = (pd.DataFrame(review_rows, columns=REVIEW_COLUMNS).drop_duplicates()
               .sort_values(["reason", "target"], kind="stable") if review_rows else pd.DataFrame(columns=REVIEW_COLUMNS))

    # The downstream table: the same hypotheses, plus what survived verification.
    annotated = frame.merge(hypothesis_status, left_on=id_column, right_on="hypothesis_id", how="left") \
        if id_column != "hypothesis_id" else frame.merge(hypothesis_status, on="hypothesis_id", how="left")
    annotated = annotated.fillna("")

    paths = {
        "source_validation": out / "source_validation.tsv",
        "source_rejected": out / "source_rejected.tsv",
        "hypothesis_source_status": out / "hypothesis_source_status.tsv",
        "mechanisms_with_verified_sources": out / "mechanisms_with_verified_sources.tsv",
        "review_required": out / "review_required.tsv",
        "verified_sources": out / "verified_sources.json",
        "source_review_records": out / "source_review_records.json",
        "claim_extraction": out / "claim_extraction_gate.tsv",
    }
    extraction.to_csv(paths["claim_extraction"], sep="\t", index=False, lineterminator="\n")
    validation.to_csv(paths["source_validation"], sep="\t", index=False, lineterminator="\n")
    validation[validation["verification_status"] == "rejected"].to_csv(
        paths["source_rejected"], sep="\t", index=False, lineterminator="\n")
    hypothesis_status.to_csv(paths["hypothesis_source_status"], sep="\t", index=False, lineterminator="\n")
    annotated.to_csv(paths["mechanisms_with_verified_sources"], sep="\t", index=False, lineterminator="\n")
    reviews.to_csv(paths["review_required"], sep="\t", index=False, lineterminator="\n")
    paths["verified_sources"].write_text(
        json.dumps({key: verified_records[key] for key in sorted(verified_records)},
                   ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    review_only_material = []
    for record in frame.to_dict("records"):
        hypothesis_id = str(record[id_column])
        for citation in sorted(set(split_citations(record.get("review_only_source_ids", "")))):
            kind, identifier = classify_citation(citation)
            if kind not in {"pmid", "pmcid", "doi"}:
                continue
            result = cache.get(_cache_key(kind, identifier))
            item = {"hypothesis_id": hypothesis_id, "citation": citation,
                    "parent_hypothesis_id": str(record.get("parent_hypothesis_id", "")) or hypothesis_id,
                    "state": "not_evaluated_cache_missing"}
            review_only_material.append(item)
            if result is None:
                continue
            from .bibliography import identifier_matches
            if not identifier_matches(kind, identifier, result):
                item["state"] = "not_evaluated_cached_identifier_mismatch"
                continue
            if citation not in review_source_records:
                review_source_records[citation] = capture_source(
                    result, citation=citation, citation_type=kind,
                    resolved_id=_record_fields(result)["resolved_id"], classification_scope="review_only_not_classified",
                )
            captured = review_source_records[citation]
            if hypothesis_id not in captured["hypothesis_ids"]:
                captured["hypothesis_ids"].append(hypothesis_id)
            review_ids = captured.setdefault("review_only_hypothesis_ids", [])
            if hypothesis_id not in review_ids:
                review_ids.append(hypothesis_id)
            lineage = captured.setdefault("review_only_lineage", [])
            link = {"hypothesis_id": hypothesis_id, "parent_hypothesis_id": item["parent_hypothesis_id"]}
            if link not in lineage:
                lineage.append(link)
            item["state"] = "captured_review_only"
            providers_used.add("local_record_cache")
    for captured in review_source_records.values():
        captured["hypothesis_ids"].sort()
        if "review_only_hypothesis_ids" in captured:
            captured["review_only_hypothesis_ids"].sort()
            captured["review_only_lineage"].sort(key=lambda row: (row["hypothesis_id"], row["parent_hypothesis_id"]))
    paths["source_review_records"].write_text(json.dumps({
        "schema_version": "1.0", "artifact_role": "source_text_for_review_only",
        "records": {key: review_source_records[key] for key in sorted(review_source_records)},
        "membership_effect": "none",
    }, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    literature_citations = int((validation["citation_type"].isin(["pmid", "pmcid", "doi", "unrecognised"])).sum()) if len(validation) else 0
    rejected = int((validation["verification_status"] == "rejected").sum()) if len(validation) else 0
    manifest = {
        "schema_version": SCHEMA_VERSION, "command": "source-verify",
        "candidate_qtl_version": __version__, "verified_at_utc": _now(),
        "provider": "+".join(sorted(providers_used)) if providers_used else provider,
        "providers_used": sorted(providers_used), "external_api_allowed": allow_external_api,
        "cache": str(runtime_path(cache_arg)) if cache_arg else None,
        "review_only_material": {
            "policy": "cached original citations only; no retrieval, classification or support promotion",
            "associations": sorted(review_only_material, key=lambda row: (row["hypothesis_id"], row["citation"])),
            "captured_association_count": sum(item["state"] == "captured_review_only" for item in review_only_material),
        },
        "hypothesis_count": int(len(frame)),
        "citation_count": int(len(validation)),
        "literature_citation_count": literature_citations,
        "verified_citation_count": int((validation["verification_status"] == "verified").sum()) if len(validation) else 0,
        "rejected_citation_count": rejected,
        "review_required_citation_count": int((validation["verification_status"] == "review_required").sum()) if len(validation) else 0,
        "database_accession_count": int((validation["verification_status"] == "not_literature_source").sum()) if len(validation) else 0,
        "unresolvable_citation_rate": round(rejected / literature_citations, 6) if literature_citations else None,
        "hypothesis_status_counts": {status: int((hypothesis_status["source_verification_status"] == status).sum())
                                     for status in HYPOTHESIS_STATUSES} if len(hypothesis_status) else {},
        "claim_support_counts": {str(key): int(value) for key, value in
                                 validation["claim_support"].value_counts().sort_index().items()} if len(validation) else {},
        "citation_verified_count": int(validation["citation_verified"].astype(str).eq("True").sum()) if len(validation) else 0,
        "claim_support_status_counts": {str(key): int(value) for key, value in
                                        validation["claim_support_status"].value_counts().sort_index().items()} if len(validation) else {},
        "evidence_relation_counts": {str(key): int(value) for key, value in
                                     validation["evidence_relation"].value_counts().sort_index().items()} if len(validation) else {},
        "retrieval_scope_counts": {str(key): int(value) for key, value in
                                   validation["retrieval_scope"].value_counts().sort_index().items()} if len(validation) else {},
        "claim_extraction": {
            "supplied": bool(claim_extraction_arg),
            "path": str(runtime_path(claim_extraction_arg)) if claim_extraction_arg else None,
            "row_count": int(len(extraction)),
            "status_counts": {str(key): int(value) for key, value in
                              extraction["extraction_status"].value_counts().sort_index().items()} if len(extraction) else {},
            "gate": ("an extracted quote must occur verbatim in the retrieved record; an unlocatable quote is "
                     "rejected whoever produced it, and an extraction never promotes an unresolved citation"),
        },
        "claim_axes": {
            "citation_verified": "did the identifier resolve to a real record",
            "claim_support_status": "what the retrieved text supports about the claim it was cited for",
            "note": ("a real paper that mentions the gene without testing the mechanism is "
                     "citation_verified=True, claim_support_status=unsupported_mention"),
        },
        "external_transmission": (
            "citation identifiers only (PMID / PMCID / DOI). The phenotype text, the interval, the candidate ids "
            "and every sequence stay on this machine"
        ) if allow_external_api else "none",
        "gates": [
            "a citation that does not resolve in the literature database is rejected: model memory is not a citable source",
            "citation resolution and claim support are separate columns: a real paper may mention a gene without testing the mechanism",
            "an evidence relation is inferred from cue terms in the retrieved text; anything short of a clear cue is review_required, never a pass",
            "an agent extraction must quote text that occurs in the retrieved record, and never creates scoring evidence",
            "gene-name co-occurrence in a real paper is unsupported_mention, not mechanistic support",
            "a citation that resolves but never mentions the claimed gene/family/pathway is review_required, not rejected",
            "an unavailable provider is review_required, never a rejection: a missing lookup is not evidence of fabrication",
            "a database accession is not a literature citation and is left to hypothesis-verify",
            "rejected citations are retained with their reason; nothing is silently dropped",
            "verification decides eligibility, never rank: no score is written here",
        ],
        "inputs": {"mechanisms": file_fingerprint(mechanisms_path, full_hash)},
        "source_review_records_fingerprint": file_fingerprint(paths["source_review_records"], True),
        "outputs": {name: str(path) for name, path in sorted(paths.items())},
    }
    manifest_path = out / "source_verification_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "manifest": str(manifest_path),
        **{name: str(path) for name, path in paths.items()},
        "citation_count": manifest["citation_count"],
        "verified_citation_count": manifest["verified_citation_count"],
        "rejected_citation_count": manifest["rejected_citation_count"],
        "review_required_citation_count": manifest["review_required_citation_count"],
        "unresolvable_citation_rate": manifest["unresolvable_citation_rate"],
        "claim_support_status_counts": manifest["claim_support_status_counts"],
        "evidence_relation_counts": manifest["evidence_relation_counts"],
    }
