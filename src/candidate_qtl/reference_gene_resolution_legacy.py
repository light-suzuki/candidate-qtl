"""Deterministic citation resolution in front of hypothesis-verify (issue #650).

A real 80-case Route A benchmark run produced a 0-row verified bait panel in
80/80 cases. The root cause traced upstream of the verification gate itself:

* ``phenotype-investigate``'s prompt requires every mechanism hypothesis to
  name a non-focal-species reference gene "with a source id (PMID/DOI/UniProt/
  TAIR accession)".
* But that stage runs interval-blind from a phenotype sentence alone, has no
  retrieval tool, and its system prompt forbids inventing citations. The
  compliant model behaviour is therefore to leave ``source_ids`` empty.
* ``hypothesis_verify`` (a04)'s ``proposal_literature_verified`` gate requires a
  non-empty ``proposal_verified_source_ids`` on every code path. That gate is
  correct and is not touched here: an unsourced proposal must stay
  ``review_required``, never be promoted on trust.

The self-contradiction is in the a01 contract, not the a04 gate: a01 demands a
citation from a tool that cannot produce one. This module is the fix -- a new,
deterministic stage inserted between ``mechanism-import`` (a02) and
``source-verify`` (a03) that supplies missing citations from an authoritative
provider (UniProt) wherever that can be done without guessing:

    mechanism-import (a02)
      -> reference-gene-resolve   (this module)
      -> source-verify (a03) -> hypothesis-verify (a04) -> bait-panel-freeze (a05)

Resolution semantics, enforced here and covered by tests:

* A complete reference name/species pair is resolved from provider metadata,
  even when the model supplied a well-formed accession or a citation. Original
  IDs and citations remain in the audit ledger; replacement identities receive
  their own provider-linked publications rather than borrowed model citations.
  When one resolved accession confirms the supplied accession, original
  citations also remain available for independent review, separately labelled
  from provider-linked publications. This is not claim verification.
* For accession-only inputs without a complete named target, supplied literature
  or unknown sources are preserved. One unambiguous UniProt/TAIR database-only
  source may retrieve linked publications with provenance.
* A real 80-case benchmark run found 242/640 ``proposed_stable_id`` values
  were fabricated pseudo-accessions such as ``Pisum_sativum__NIN`` -- the
  model treated the field as "the focal-species gene I am hypothesising"
  instead of "the reference gene's stable accession". Because
  ``mechanism-import`` gates ``source_ids`` but never validates
  ``proposed_stable_id``'s format, these reached this stage unchallenged and,
  in the general benchmark analysis, would have gone on to hypothesis_verify
  as ``supplied_by_proposal`` and queried UniProt with a garbage string
  (HTTP 400). A ``proposed_stable_id`` that is NOT a recognisable UniProt
  accession or AGI locus is therefore treated as ABSENT for resolution
  purposes -- never queried -- and resolution falls back to
  ``reference_gene_name`` + ``reference_species`` exactly like the clean
  path. The original malformed string is preserved verbatim in the ledger's
  ``input_proposed_stable_id`` column with its own ``malformed_proposed_stable_id_status``
  (``malformed_proposed_stable_id_ignored``, or
  ``malformed_proposed_stable_id_names_focal_species`` when the fabricated
  string itself names the run's focal species -- a blinding-relevant event
  that must stay visible, not be buried inside an ordinary-malformation
  bucket). The written ``mechanism_hypotheses_resolved.tsv`` never carries a
  malformed string forward: a malformed ``proposed_stable_id`` cell is always
  cleaned to either a provider-resolved accession or empty.
* A gene name that resolves to multiple exact provider records is never
  tie-broken. Every accession is emitted as a deterministic child hypothesis
  and verified independently; ambiguous parentage remains visible in lineage.
* Offline (no ``--allow-external-api``), a bare gene name/symbol can still be
  resolved against a frozen ``--records`` cache through a gene-symbol index
  built the same way ``llm_pathway._cached_tair_index`` builds its TAIR
  index: only from the record's own ``genes[].geneName.value`` /
  ``genes[].synonyms[].value`` fields, matched to the record's own organism
  genus. "A cache is an authoritative record snapshot, not an alias guess" --
  every matching accession is retained so an ambiguous cache fails closed
  instead of picking one.
* A provider failure (network error, timeout, HTTP error) is recorded as
  ``not_evaluated`` -- never as a rejection or a negative result -- in both the
  ledger and the manifest, matching this repository's not_evaluated idiom
  (see ``blast_panel.REVIEW_STATUS_VOCAB`` and ``hypothesis_verify``'s
  ``accession_provider_unavailable`` review reason).
* Without ``--allow-external-api`` and without a cache hit (accession-keyed or
  gene-symbol-indexed), nothing is attempted at all: ``not_attempted_offline``,
  never a guess.
* A record with no source-backed function citation contributes nothing; the
  hypothesis stays unresolved and review-required rather than being filled
  with an uncited accession.
* Every filled cell carries a machine-readable ``retrieval_provenance`` JSON
  object (provider, exact query, resolution status, hit count, UTC retrieval
  timestamp, record version) so a model-supplied value and a
  provider-resolved value can never be confused downstream.

This module intentionally reuses ``llm_pathway``'s existing UniProt resolution
primitives (``_resolve_accession``, ``_resolve_tair_locus``, ``_fetch_uniprot``,
``_load_record_cache``, ``_cached_tair_index``, ``_searchable_symbol``, record
accessors) rather than reimplementing or moving them. ``hypothesis_verify``,
``source_verify`` and ``blast_panel`` are not imported for their logic and are
not modified; their behaviour on unchanged inputs is unaffected by this
module's existence.
"""

from __future__ import annotations

import hashlib
import inspect
import json
import re
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

from .core import file_fingerprint, runtime_path
from .identifiers import canonical_accession, is_agi_locus, is_uniprot_accession
from .llm_pathway import (
    UNIPROT_SEARCH_ENDPOINT,
    UNIPROT_ENTRY_ENDPOINT,
    _cached_tair_index,
    _fetch_uniprot,
    _load_record_cache,
    _now,
    _record_source_ids,
    _record_version,
    _resolve_accession,
    _resolve_accession_candidates,
    _resolve_tair_locus,
    _searchable_symbol,
)
from .phenotype_investigation import MECHANISM_COLUMNS
from .reference_name_pairs import (IDENTITY_DECLARED_ACCESSION, IDENTITY_EXACT, IDENTITY_MISMATCH,
                                   IDENTITY_REVIEWED_ALIAS, PROJECTABLE_IDENTITIES,
                                   atomic_reference_pairs, declared_record_identity,
                                   declared_uniprot_accessions, exact_record_candidates, load_identity_reviews,
                                   record_gene_names,
                                   reference_handle, reference_identity_names, record_reference_handles,
                                   reference_name_for_input,
                                   reference_lookup_handles,
                                   reference_lookup_names_for_inputs)
from .source_verify import classify_citation

SCHEMA_VERSION = "1.2"
DEFAULT_TIMEOUT = 60  # matches llm_pathway._resolve_accession / _resolve_tair_locus / _fetch_uniprot


def _resolve_reference_candidates(gene_name: str, species: str, raw_inputs: list[dict], timeout: int):
    """Search complete name forms and typed ID handles without selecting a record.

    Returned accessions are only lookup candidates. Callers must fetch the
    authoritative record and apply exact_record_candidates before adoption.
    """
    accessions: set[str] = set()
    primary_accessions: set[str] = set()
    queries: list[str] = []
    query_urls: list[str] = []
    failures: list[str] = []

    def search(query: str, citation_id: str = "") -> None:
        cursor = "*"
        seen: set[str] = set()
        while cursor:
            if cursor in seen:
                raise RuntimeError("UniProt reference search repeated pagination cursor")
            seen.add(cursor)
            params = {"query": query, "format": "json", "size": 500}
            if cursor != "*":
                params["cursor"] = cursor
            url = UNIPROT_SEARCH_ENDPOINT + "?" + urllib.parse.urlencode(params)
            if url not in query_urls:
                query_urls.append(url)
            request = urllib.request.Request(url, headers={"User-Agent": "candidate-qtl reference-gene-resolve"})
            try:
                with urllib.request.urlopen(request, timeout=timeout) as response:
                    payload = json.loads(response.read().decode("utf-8"))
                    links = response.headers.get("Link", "")
            except urllib.error.HTTPError as error:
                raise RuntimeError(f"UniProt API HTTP {error.code} for reference query {query!r}") from error
            except urllib.error.URLError as error:
                raise RuntimeError(f"UniProt API connection failed for reference query {query!r}: {error.reason}") from error
            if not isinstance(payload, dict) or not isinstance(payload.get("results"), list):
                raise ValueError(f"UniProt API returned an invalid payload for reference query {query!r}")
            for item in payload["results"]:
                if isinstance(item, dict) and item.get("primaryAccession"):
                    accession = str(item["primaryAccession"]).strip()
                    if citation_id:
                        # A citation finds records; it is neither a gene ID nor
                        # evidence of trait support. Require explicit record
                        # identity and the requested citation before forwarding.
                        if not exact_record_candidates(
                                {accession: item}, lookup_names, species):
                            continue
                        if f"PMID:{citation_id}" not in _record_source_ids(item):
                            continue
                    accessions.add(accession)
            next_link = re.search(r'<([^>]+)>;\s*rel="next"', links)
            cursor = ""
            if next_link:
                parsed = urllib.parse.urlparse(next_link.group(1))
                cursor = urllib.parse.parse_qs(parsed.query).get("cursor", [""])[0]
                if not cursor:
                    raise ValueError("UniProt reference search next-page link has no cursor; search is incomplete")

    # Keep the established symbol lookup and add complete phrases from the
    # original display name (including parenthetical formal-name alternatives).
    primary_trace: dict = {}
    try:
        # Keep 3-argument injected resolvers working; production also records
        # actual pagination URLs and exact partial hits through this trace.
        supports_trace = "query_trace" in inspect.signature(_resolve_accession_candidates).parameters
        if supports_trace:
            primary, status, query, _hits = _resolve_accession_candidates(
                gene_name, species, timeout, query_trace=primary_trace)
        else:
            primary, status, query, _hits = _resolve_accession_candidates(gene_name, species, timeout)
        queries.append(query or "")
        query_urls.extend(primary_trace.get("query_urls", []))
        if query and not primary_trace.get("query_urls"):
            query_urls.append(UNIPROT_SEARCH_ENDPOINT + "?" + urllib.parse.urlencode(
                {"query": query, "format": "json", "size": 500}))
        primary_accessions.update(primary)
        accessions.update(primary)
    except (RuntimeError, TimeoutError, OSError, ValueError) as error:
        status = "provider_unavailable"
        if primary_trace.get("query"):
            queries.append(primary_trace["query"])
        failures.append(f"symbol_lookup: {type(error).__name__}: {error}")
        query_urls.extend(primary_trace.get("query_urls", []))
        primary_accessions.update(primary_trace.get("accessions", []))
        accessions.update(primary_trace.get("accessions", []))
    lookup_names = reference_lookup_names_for_inputs(raw_inputs, gene_name, species)
    for name in lookup_names:
        if (name.casefold() == str(gene_name).strip().casefold()
                and not re.search(r"\s", name) and status != "no_gene_name_to_resolve"):
            continue
        formal_query = f'gene:"{name.replace(chr(34), "")}" AND organism_name:"{species.strip()}"'
        if formal_query not in queries:
            queries.append(formal_query)
            try:
                search(formal_query)
                if status == "no_gene_name_to_resolve":
                    status = "unresolved_gene_name"
            except (RuntimeError, TimeoutError, OSError, ValueError) as error:
                failures.append(f"{formal_query}: {type(error).__name__}: {error}")

    for kind, handle in reference_lookup_handles(raw_inputs):
        if kind == "uniprot":
            # Direct entry fetch follows; the record's primary accession,
            # declared name, species and sequence remain independently checked.
            accessions.add(handle)
            queries.append(f"handle:uniprot:{handle}")
        else:
            handle_query = f'xref:{kind}-{handle} AND organism_name:"{species.strip()}"'
            queries.append(handle_query)
            try:
                search(handle_query)
            except (RuntimeError, TimeoutError, OSError, ValueError) as error:
                failures.append(f"{handle_query}: {type(error).__name__}: {error}")

    # Only a completed, empty name/handle search falls back to the hypothesis's
    # own citations. Never broaden species, invent aliases, or mask a failed
    # provider query. Entry fetch and downstream sequence checks still follow.
    if not accessions and not failures:
        pmids = sorted({identifier
                        for row in raw_inputs
                        for token in str(row.get("source_ids", "")).split(";")
                        for kind, identifier in [classify_citation(token)]
                        if kind == "pmid"})
        for pmid in pmids:
            citation_query = f'lit_pubmed:{pmid} AND organism_name:"{species.strip()}"'
            queries.append(citation_query)
            try:
                search(citation_query, citation_id=pmid)
            except (RuntimeError, TimeoutError, OSError, ValueError) as error:
                failures.append(f"{citation_query}: {type(error).__name__}: {error}")

    # Deterministic query ledger text also records alternatives and partial
    # failures while preserving every discovered accession for record checks.
    query_text = " | ".join(dict.fromkeys(value for value in queries if value))
    if not query_text and failures:
        query_text = "reference_lookup:provider_unavailable"
    if failures:
        status = "provider_unavailable"
    elif accessions and accessions == primary_accessions:
        status = status or "resolved_from_gene_name"
    elif len(accessions) == 1:
        status = "resolved_from_reference_lookup"
    elif len(accessions) > 1:
        status = "ambiguous_gene_name"
    else:
        status = status or "unresolved_gene_name"
    return sorted(accessions), status, query_text, len(accessions), failures, query_urls

# How this stage arrived at the citations/accession recorded for a hypothesis.
# model_supplied: existing sources preserved without provider enrichment.
# provider_resolved: this stage filled proposed_stable_id and/or source_ids from
#   an authoritative record with a source-backed citation.
# not_resolved: nothing was filled (ambiguous, provider failure, offline without
#   a cache hit, no citation on the retrieved record, or nothing to resolve).
RESOLUTION_ORIGINS = ("model_supplied", "provider_resolved", "not_resolved")
# not_attempted: no provider call was made (already sourced, nothing to resolve,
#   or offline without a cache hit).
# evaluated: a provider call (cache or API) completed and produced a definitive
#   answer (a fill, an ambiguity, a not-found, or a no-citation record).
# not_evaluated: a provider call was attempted but failed (network/timeout/HTTP);
#   this is never treated as a negative result anywhere downstream.
EVALUATION_STATUSES = ("not_attempted", "evaluated", "not_evaluated")
# A malformed proposed_stable_id is never a rejection reason on its own; it is
# an input-quality event recorded alongside whatever the gene-name fallback
# resolution produces.
MALFORMED_STATUSES = ("", "malformed_proposed_stable_id_ignored", "malformed_proposed_stable_id_names_focal_species")

REVIEW_COLUMNS = ["reason", "target", "required_next_action", "allowed_choices"]
REFERENCE_GENE_RESOLUTION_COLUMNS = [
    "hypothesis_id", "case_id", "resolution_origin",
    "reference_gene_name", "reference_species",
    "input_proposed_stable_id", "input_source_ids",
    "malformed_proposed_stable_id_status",
    "resolved_stable_id", "resolved_source_ids",
    "accession_resolution_status", "accession_resolution_query", "accession_resolution_hit_count",
    "record_found", "function_source_backed",
    "evaluation_status",
    "source_provider", "source_url", "retrieved_at", "source_record_version",
    "provider_error", "review_reason", "retrieval_provenance",
    "reference_identity_status", "reference_identity_match_basis",
    "parent_hypothesis_id", "derived_accession",
    "input_reference_gene_name", "input_reference_species", "atomic_pair_status",
]
RESOLVED_MECHANISM_COLUMNS = MECHANISM_COLUMNS + [
    "review_only_source_ids",
    "reference_identity_status", "reference_identity_match_basis",
    "parent_hypothesis_id", "derived_hypothesis_index", "derived_hypothesis_count",
    "derived_accession",
    "input_reference_gene_name", "input_reference_species", "atomic_pair_status",
]


def _frame(rows: list[dict[str, object]], columns: list[str]) -> pd.DataFrame:
    if not rows:
        return pd.DataFrame(columns=columns)
    frame = pd.DataFrame(rows)
    for column in columns:
        if column not in frame.columns:
            frame[column] = ""
    return frame[columns]


def _write(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


def _provenance_json(provider: str, query: str, resolution_status: str, hit_count: int,
                      retrieved_at: str, record_version: str) -> str:
    payload = {
        "provider": provider, "query": query, "resolution_status": resolution_status,
        "hit_count": int(hit_count), "retrieved_at": retrieved_at, "record_version": record_version,
    }
    return json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _with_lookup_provenance(serialized: str, *, search_state: str, query_urls: list[str],
                            failures: list[str], record_identity_state: str,
                            record_fetch_failures: dict | None = None,
                            record_identity_failures: dict | None = None) -> str:
    provenance = json.loads(serialized or "{}")
    provenance.update(lookup_search_state=search_state, query_urls=list(query_urls),
                      lookup_failures=list(failures), record_identity_state=record_identity_state)
    if record_fetch_failures:
        provenance["record_fetch_failures"] = record_fetch_failures
    if record_identity_failures:
        provenance["record_identity_failures"] = record_identity_failures
    return json.dumps(provenance, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _cache_lookup_receipt(receipt: dict) -> tuple[str, list[str], list[str], dict, dict]:
    """Interpret current and older cache receipts without turning incomplete into no-hit."""
    query_urls = list(receipt.get("query_urls", []))
    lookup_failures = list(receipt.get("lookup_failures", []))
    record_failures = dict(receipt.get("record_fetch_failures", receipt.get("failed_accessions", {})))
    identity_failures = dict(receipt.get("record_identity_failures", {}))
    # The old name/species error followed a valid dictionary/primary check.
    # Its generic identity error also covered non-dict payloads, so that legacy
    # message remains a fetch failure unless an explicit identity field exists.
    for accession, reason in list(record_failures.items()):
        if reason == f"ValueError: provider record name/species mismatch: {accession}":
            identity_failures.setdefault(accession, record_failures.pop(accession))
        elif accession in identity_failures and reason in (
                identity_failures[accession], f"ValueError: provider record identity mismatch: {accession}"):
            record_failures.pop(accession)
    status = str(receipt.get("status", ""))
    reason = str(receipt.get("reason", ""))
    attempted = bool(receipt.get("logical_query_attempts", 0) or receipt.get("query") or query_urls)
    # A record limit can coexist with failed lookup alternatives. Preserve
    # both axes; a complete record roster is not a complete search receipt.
    if status == "not_evaluated" and reason == "record_budget_exceeded":
        record_failures.setdefault("cache", reason)
    if lookup_failures:
        search_state = "partial_not_evaluated"
    elif status == "not_evaluated" and not attempted:
        search_state = "not_evaluated"
    elif status == "not_evaluated" and (record_failures or identity_failures):
        search_state = "complete"
    elif status == "not_evaluated":
        search_state = "partial_not_evaluated"
        if reason:
            lookup_failures.append(reason)
    elif status in {"complete", "no_provider_match"}:
        search_state = "complete"
    else:
        search_state = "not_evaluated"
    return search_state, query_urls, lookup_failures, record_failures, identity_failures


def _is_well_formed_accession(value: str) -> bool:
    """Whether a proposed_stable_id is a recognisable UniProt accession or AGI locus.

    Mirrors the format acceptance of ``llm_pathway._declared_source_accessions``
    (``is_uniprot_accession(accession) or is_agi_locus(accession)`` after
    canonicalisation). A fabricated pseudo-accession such as
    ``Pisum_sativum__NIN`` never matches either grammar and is treated as
    absent rather than queried against the provider (issue #649 amendment 1).
    """
    canonical = canonical_accession(value)
    return bool(canonical) and (is_uniprot_accession(canonical) or is_agi_locus(canonical))


def _names_focal_species(value: str, focal_species: str) -> bool:
    """Whether a free-text value appears to name the run's focal species.

    Genus-token comparison only -- the same granularity ``hypothesis_verify``
    already uses to exclude a focal-species bait (``retrieved_genus ==
    focal_genus``) -- and parameterised entirely by the ``focal_species``
    argument threaded in from the caller. Never a hardcoded species string.
    """
    if not focal_species or not value:
        return False
    focal_genus = focal_species.split(" ")[0].casefold()
    if not focal_genus:
        return False
    tokens = [token.casefold() for token in re.split(r"[^A-Za-z0-9]+", value) if token]
    if focal_genus in tokens:
        return True
    normalized_value = " ".join(tokens)
    normalized_species = " ".join(focal_species.split()).casefold()
    return bool(normalized_species) and normalized_species in normalized_value


def _cached_gene_symbol_index(records: dict[str, dict[str, object]]) -> dict[tuple[str, str], list[str]]:
    """Index gene symbols declared in frozen UniProt records, keyed by (genus, symbol).

    A cache is an authoritative record snapshot, not an alias guess: only the
    record's own ``genes[].geneName.value`` and ``genes[].synonyms[].value``
    fields are indexed, matched to the record's own organism genus, and every
    matching accession is retained so an ambiguous cache fails closed instead
    of selecting one -- exactly the discipline
    ``llm_pathway._cached_tair_index`` already applies to TAIR cross-references.
    A free-text gene name is normalised through the existing
    ``_searchable_symbol`` before lookup, never a new parser.
    """
    index: dict[tuple[str, str], set[str]] = {}
    for accession, record in records.items():
        organism = record.get("organism", {})
        species = str(organism.get("scientificName", "")) if isinstance(organism, dict) else ""
        genus = species.split(" ")[0].casefold() if species else ""
        if not genus:
            continue
        genes = record.get("genes", [])
        names: list[str] = []
        for gene in genes if isinstance(genes, list) else []:
            if not isinstance(gene, dict):
                continue
            gene_name = gene.get("geneName", {})
            if isinstance(gene_name, dict) and gene_name.get("value"):
                names.append(str(gene_name["value"]))
            synonyms = gene.get("synonyms", [])
            for synonym in synonyms if isinstance(synonyms, list) else []:
                if isinstance(synonym, dict) and synonym.get("value"):
                    names.append(str(synonym["value"]))
        for name in names:
            symbol = _searchable_symbol(name)
            if symbol:
                index.setdefault((genus, symbol.casefold()), set()).add(accession)
    return {key: sorted(accessions) for key, accessions in index.items()}


def _cached_gene_symbol_candidates_exact(records: dict[str, dict[str, object]], gene_name: str,
                                         species: str) -> list[str]:
    return exact_record_candidates(records, gene_name, species)


def _resolve_one(
    row: dict[str, object],
    cache: dict[str, dict[str, object]],
    cache_path_str: str,
    cached_tair: dict[str, list[str]],
    cached_symbols: dict[tuple[str, str], list[str]],
    allow_external_api: bool,
    timeout: int,
    focal_species: str,
    candidate_override: dict[str, object] | None = None,
) -> dict[str, object]:
    """Resolve one mechanism hypothesis row into a ledger record.

    Never called for a row whose ``source_ids`` is already non-empty -- the
    caller short-circuits that case before this function runs, so a
    model-supplied citation is never re-touched.
    """
    hypothesis_id = str(row.get("hypothesis_id", ""))
    case_id = str(row.get("case_id", ""))
    proposed_stable_id = str(row.get("proposed_stable_id", "")).strip()
    reference_species = str(row.get("reference_species", "")).strip()
    reference_gene_name = str(row.get("reference_gene_name", "")).strip()

    # --- Phase 0: a malformed proposed_stable_id is never queried; it is
    # treated as absent and resolution falls back to the gene name exactly
    # like a hypothesis that never supplied an accession at all. ---
    supplied = canonical_accession(proposed_stable_id) if proposed_stable_id else ""
    well_formed = bool(proposed_stable_id) and _is_well_formed_accession(proposed_stable_id)
    malformed_status = ""
    effective_proposed_stable_id = proposed_stable_id
    if proposed_stable_id and not well_formed:
        malformed_status = ("malformed_proposed_stable_id_names_focal_species"
                            if _names_focal_species(proposed_stable_id, focal_species)
                            else "malformed_proposed_stable_id_ignored")
        effective_proposed_stable_id = ""
        supplied = ""

    # --- Phase 1: find an accession to fetch a record for. ---
    accession, resolution_status, query, hits, provider_error = "", "", "", 0, ""
    if candidate_override is not None:
        accession = str(candidate_override.get("accession", ""))
        resolution_status = str(candidate_override.get("status", "unresolved_gene_name"))
        query = str(candidate_override.get("query", ""))
        hits = int(candidate_override.get("hits", 0) or 0)
        provider_error = str(candidate_override.get("provider_error", ""))
    elif effective_proposed_stable_id and is_agi_locus(supplied):
        if allow_external_api:
            try:
                accession, resolution_status, query, hits = _resolve_tair_locus(supplied, reference_species, timeout)
            except (RuntimeError, TimeoutError, OSError, ValueError) as error:
                resolution_status, provider_error = "provider_unavailable", str(error)
        else:
            cached_hits = cached_tair.get(supplied, [])
            query, hits = f"record_cache:xref:tair-{supplied}", len(cached_hits)
            if len(cached_hits) == 1:
                accession, resolution_status = cached_hits[0], "resolved_from_tair_locus_cache"
            elif len(cached_hits) > 1:
                resolution_status = "ambiguous_tair_locus"
            else:
                accession, resolution_status = supplied, "not_attempted_offline"
    elif effective_proposed_stable_id:
        # Use the supplied, well-formed accession as-is; never re-derive it
        # from the gene name.
        accession, resolution_status = supplied, "supplied_by_proposal"
    elif reference_gene_name:
        if allow_external_api:
            try:
                accession, resolution_status, query, hits = _resolve_accession(reference_gene_name, reference_species, timeout)
            except (RuntimeError, TimeoutError, OSError, ValueError) as error:
                resolution_status, provider_error = "provider_unavailable", str(error)
        else:
            symbol = _searchable_symbol(reference_gene_name)
            genus = reference_species.split(" ")[0].casefold() if reference_species else ""
            if symbol and genus:
                query = f"record_cache:gene:{symbol} AND organism_genus:{genus}"
                cached_hits = cached_symbols.get((genus, symbol.casefold()), [])
                hits = len(cached_hits)
                if len(cached_hits) == 1:
                    accession, resolution_status = cached_hits[0], "resolved_from_gene_name_cache"
                elif len(cached_hits) > 1:
                    resolution_status = "ambiguous_gene_name"
                else:
                    resolution_status = "not_attempted_offline"
            else:
                resolution_status = "not_attempted_offline"
    else:
        resolution_status = "no_gene_name_to_resolve"

    # --- Phase 2: retrieve the record for a discovered accession. ---
    record: dict[str, object] | None = None
    provider, source_url, retrieved_at = "", "", ""
    if accession and not provider_error:
        if accession in cache:
            record = cache[accession]
            provider, source_url = "local_record_cache", cache_path_str
            retrieved_at = str(record.get("candidate_qtl_retrieved_at", "")) or _now()
        elif allow_external_api:
            try:
                record = _fetch_uniprot(accession, timeout)
                provider = "uniprot_rest"
                source_url = UNIPROT_ENTRY_ENDPOINT + accession + ".json"
                retrieved_at = _now()
            except (RuntimeError, TimeoutError, OSError, ValueError) as error:
                provider_error = str(error)
        elif resolution_status != "not_attempted_offline":
            # The accession is known (supplied directly, or from an offline
            # cache index) but there is neither a cache hit nor API
            # permission: nothing was actually attempted.
            resolution_status = "not_attempted_offline"

    identity_names = (candidate_override or {}).get("identity_names", [reference_gene_name])
    record_identity = ""
    if record and candidate_override is not None and reference_gene_name and reference_species:
        record_identity = declared_record_identity(
            record, accession, identity_names, reference_species,
            set(candidate_override.get("declared_accessions", ())),
            candidate_override.get("identity_reviews") or {})
    identity_mismatch = record_identity == IDENTITY_MISMATCH
    if identity_mismatch:
        record = None
        resolution_status = "provider_record_identity_mismatch"
    record_source_ids = _record_source_ids(record) if record else []
    publication_ids = [value for value in record_source_ids
                       if classify_citation(value)[0] in {"pmid", "pmcid", "doi"}]
    record_version = _record_version(record) if record else ""

    # --- Phase 3: classify the outcome. Never tie-break, never fabricate. ---
    origin, review_reason, evaluation_status = "not_resolved", "", "evaluated"
    if provider_error:
        review_reason, evaluation_status = "accession_provider_unavailable", "not_evaluated"
    elif resolution_status in ("ambiguous_gene_name", "ambiguous_tair_locus"):
        review_reason = ("declared_tair_locus_resolves_to_multiple_accessions"
                         if resolution_status == "ambiguous_tair_locus"
                         else "reference_gene_name_resolves_to_multiple_accessions")
    elif resolution_status == "not_attempted_offline":
        review_reason, evaluation_status = "not_attempted_offline", "not_attempted"
    elif resolution_status == "no_gene_name_to_resolve":
        review_reason, evaluation_status = "no_reference_gene_to_resolve", "not_attempted"
    elif resolution_status == "provider_record_identity_mismatch":
        review_reason = "provider_record_identity_mismatch"
    elif not accession:
        review_reason = "reference_gene_not_resolved"
    elif identity_mismatch:
        review_reason = "provider_record_identity_mismatch"
    elif record is None:
        review_reason = "accession_not_found_in_authoritative_database"
    elif not publication_ids:
        review_reason = "reference_gene_function_is_not_source_backed"
    else:
        origin = "provider_resolved"

    resolved_stable_id = accession if origin == "provider_resolved" else ""
    resolved_source_ids = ";".join(record_source_ids) if origin == "provider_resolved" else ""

    provenance = _provenance_json(provider, query, resolution_status, hits, retrieved_at, record_version)
    if record:
        from .gene_metadata import uniprot_gene_metadata
        provenance_data = json.loads(provenance)
        provenance_data["validated_reference_handles"] = [list(value) for value in sorted(record_reference_handles(record))]
        provenance_data["validated_reference_metadata"] = {**uniprot_gene_metadata(record), "provider": provider}
        provenance = json.dumps(provenance_data, sort_keys=True, separators=(",", ":"))
    result = {
        "hypothesis_id": hypothesis_id, "case_id": case_id, "resolution_origin": origin,
        "reference_gene_name": reference_gene_name, "reference_species": reference_species,
        "input_proposed_stable_id": proposed_stable_id, "input_source_ids": "",
        "malformed_proposed_stable_id_status": malformed_status,
        "resolved_stable_id": resolved_stable_id, "resolved_source_ids": resolved_source_ids,
        "accession_resolution_status": resolution_status, "accession_resolution_query": query,
        "accession_resolution_hit_count": hits,
        "record_found": bool(record), "function_source_backed": bool(publication_ids),
        "evaluation_status": evaluation_status,
        "source_provider": provider, "source_url": source_url, "retrieved_at": retrieved_at,
        "source_record_version": record_version,
        "provider_error": provider_error, "review_reason": review_reason,
        "retrieval_provenance": provenance,
    }
    if record_identity and record:
        result["reference_identity_status"] = record_identity
    return result


def _review_row(hypothesis_id: str, review_reason: str, hits: int) -> dict[str, str]:
    detail = {
        "not_attempted_offline": (
            "no --allow-external-api and this accession/gene name was not found in --records; "
            "supply an offline record cache, allow the API, or resolve the citation manually"),
        "no_reference_gene_to_resolve": (
            "the hypothesis has neither a proposed_stable_id nor a reference_gene_name to resolve; "
            "supply one or keep the hypothesis unverifiable"),
        "accession_provider_unavailable": (
            "the UniProt provider was unreachable, timed out, or returned an error; "
            "this is not evaluated as a negative result -- retry or supply an offline record cache"),
        "reference_gene_name_resolves_to_multiple_accessions": (
            f"{hits} reviewed UniProt accession(s) match this gene name for the declared species; "
            "never tie-broken -- supply the correct accession manually"),
        "declared_tair_locus_resolves_to_multiple_accessions": (
            f"{hits} reviewed UniProt accession(s) cross-reference this TAIR locus; "
            "never tie-broken -- supply the correct accession manually"),
        "reference_gene_not_resolved": (
            "no eligible UniProt record matches this gene name or locus for the declared species; "
            "supply an accession manually or keep the hypothesis unverifiable"),
        "accession_not_found_in_authoritative_database": (
            "the accession does not exist in the authoritative database; correct it or keep the hypothesis unverifiable"),
        "reference_gene_function_is_not_source_backed": (
            "the retrieved record carries no PMID-backed function citation; "
            "supply a literature source manually or keep the hypothesis unverifiable"),
    }.get(review_reason, "resolve the citation gap before this hypothesis can be treated as sourced")
    return {
        "reason": f"reference_gene_resolution_required:{review_reason}", "target": hypothesis_id,
        "required_next_action": detail,
        "allowed_choices": "supply a record cache | allow the database API | supply the citation manually | keep as unverifiable proposal",
    }


def _malformed_review_row(hypothesis_id: str, malformed_status: str, original_value: str) -> dict[str, str]:
    if malformed_status == "malformed_proposed_stable_id_names_focal_species":
        detail = (
            f"proposed_stable_id '{original_value}' is not a recognisable UniProt accession or AGI locus AND "
            "appears to name this run's focal species -- a blinding-relevant fabrication, not an ordinary "
            "formatting slip. It was never queried; resolution fell back to reference_gene_name + "
            "reference_species. The original string is preserved verbatim in input_proposed_stable_id."
        )
    else:
        detail = (
            f"proposed_stable_id '{original_value}' is not a recognisable UniProt accession or AGI locus and "
            "was never queried; resolution fell back to reference_gene_name + reference_species. The original "
            "string is preserved verbatim in input_proposed_stable_id."
        )
    return {
        "reason": f"reference_gene_resolution_required:{malformed_status}", "target": hypothesis_id,
        "required_next_action": detail,
        "allowed_choices": "correct proposed_stable_id manually | rely on the reference_gene_name fallback | discard the proposal",
    }


def _resolved_paper_anchor_basis(report: dict, row: dict, provenance: dict,
                                 *, unique_reference: bool, unique_handles: set[tuple[str, str]],
                                 supplied_id_conflict: bool, shared_handles: set[tuple[str, str]] | None = None) -> str:
    """Transport an existing judgment only across a verified reference endpoint.

    A display label is scoped to its original atomic pair, never registered as
    a global alias. Name-only anchors cannot choose among provider copies.
    """
    anchor = str(report.get("reference_gene_id", "")).strip()
    handles = {tuple(value) for value in provenance.get("validated_reference_handles", [])}
    handle = reference_handle(anchor, allow_untyped=True)
    if handle:
        primary = reference_handle(provenance.get("validated_reference_accession", ""), allow_untyped=True)
        if handle in handles and (handle == primary or handle in unique_handles):
            return "provider_record_cross_reference"
        if handle in handles and handle in (shared_handles or set()):
            return "provider_record_shared_gene_" + handle[0]
        return ""
    if not unique_reference or supplied_id_conflict:
        return ""
    normalize = lambda value: " ".join(str(value).split()).casefold()
    species = str(row.get("reference_species", ""))
    if (report.get("gene_correspondence") != "homolog" and report.get("reported_species")
            and normalize(report["reported_species"]) != normalize(species)):
        return ""
    names = reference_lookup_names_for_inputs([row], str(row.get("reference_gene_name", "")), species)
    validated_names = {normalize(value) for value in provenance.get("validated_reference_metadata", {}).get("names", [])}
    # Every name accepted here must also occur in the provider record. A label
    # supplied only by the paper report is not an authoritative alias.
    if normalize(anchor) in {normalize(value) for value in names} & validated_names:
        return "provider_record_declared_alias"
    raw_name = str(row.get("input_reference_gene_name") or row.get("reference_gene_name", ""))
    raw_species = str(row.get("input_reference_species") or species)
    pairs, _ = atomic_reference_pairs(raw_name, raw_species)
    anchor_pairs, _ = atomic_reference_pairs(anchor, species)
    current = (normalize(row.get("reference_gene_name", "")), normalize(species))
    # The complete original single-reference label may have a parenthetical
    # short name validated by the resolver. Composite/sibling labels cannot
    # borrow that validation, nor can a conflicting explicit ID in the label.
    if (len(pairs) == 1 and len(anchor_pairs) == 1 and normalize(anchor) == normalize(raw_name)
            and tuple(map(normalize, anchor_pairs[0])) == current
            and {normalize(value) for value in names} & validated_names):
        explicit = [reference_handle(value, allow_untyped=True) for value in re.findall(r"\(([^()]*)\)", anchor)]
        if all(value is None or value in handles for value in explicit):
            return "provider_record_original_atomic_label"
    return ""


def reference_gene_resolve(
    mechanisms_arg: str,
    out_dir_arg: str,
    allow_external_api: bool = False,
    record_cache_arg: str | None = None,
    timeout: int = DEFAULT_TIMEOUT,
    focal_species: str = "",
    identity_review_arg: str | None = None,
) -> dict[str, object]:
    """Deterministically fill missing mechanism-hypothesis citations from UniProt.

    Reads ``mechanism_hypotheses.tsv`` (a02 output, ``MECHANISM_COLUMNS`` shape)
    and writes, under ``out_dir_arg``:

    * ``mechanism_hypotheses_resolved.tsv`` -- the same canonical schema, with
      ``source_ids`` / ``proposed_stable_id`` / ``retrieval_provenance`` filled
      in only for rows this stage could resolve without guessing. A malformed
      ``proposed_stable_id`` is cleaned where resolution is attempted. Named
      targets retain every exact provider alternative; unresolved names never
      silently fall back to the model accession.
    * ``reference_gene_resolution.tsv`` -- one audit-ledger row per resolved alternative,
      including the verbatim malformed value and its own status when
      applicable.
    * ``review_required.tsv`` -- this repository's standard 4-column ledger.
    * ``reference_gene_resolution_manifest.json`` -- run manifest with counts,
      config, and an explicit not-evaluated list (never silently a rejection).

    ``focal_species`` is optional and used only to flag, as its own distinct
    ledger status, a malformed ``proposed_stable_id`` that names the focal
    species -- it is a blinding diagnostic, not a verification gate (that gate
    is ``hypothesis_verify``'s and is untouched here).
    """
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    mechanisms_path = runtime_path(mechanisms_arg)
    mechanisms = pd.read_csv(mechanisms_path, sep="\t", dtype=str).fillna("")
    if "causal_layer" not in mechanisms.columns:
        mechanisms.insert(2, "causal_layer", "")
    missing = sorted({"hypothesis_id", "case_id", "source_ids", "proposed_stable_id",
                      "reference_species", "reference_gene_name"} - set(mechanisms.columns))
    if missing:
        raise ValueError(f"mechanism table is missing columns: {', '.join(missing)}")

    original_count = len(mechanisms)
    original_ids = [(str(row["case_id"]), str(row["hypothesis_id"])) for row in mechanisms.to_dict("records")]
    if (any(not all(pair) or any(value != value.strip() for value in pair) for pair in original_ids)
            or len({tuple(value.casefold() for value in pair) for pair in original_ids}) != len(original_ids)):
        raise ValueError("mechanism IDs must be nonempty, unique, and case-exact per case")
    expanded = []
    used_ids = {(case.casefold(), identifier.casefold()) for case, identifier in original_ids}
    for row in mechanisms.to_dict("records"):
        original_name, species = str(row["reference_gene_name"]), str(row["reference_species"])
        name = reference_name_for_input(row)
        pairs, pair_status = atomic_reference_pairs(name, species)
        for pair in pairs or [(name, species)]:
            atomic = dict(row)
            if len(pairs) > 1:
                suffix = hashlib.sha256(json.dumps(pair, ensure_ascii=False).encode()).hexdigest()[:12]
                child_id = f"{row['hypothesis_id']}__pair_{suffix}"
                child_key = (str(row["case_id"]).casefold(), child_id.casefold())
                if child_key in used_ids:
                    raise ValueError(f"derived hypothesis id collides with input roster: {child_id}")
                used_ids.add(child_key)
                atomic["hypothesis_id"] = child_id
            atomic.update(reference_gene_name=pair[0], reference_species=pair[1],
                          input_reference_gene_name=original_name, input_reference_species=species,
                          atomic_pair_status=pair_status, _original_parent_id=row["hypothesis_id"])
            expanded.append(atomic)
    mechanisms = pd.DataFrame(expanded, columns=list(mechanisms.columns) + ["input_reference_gene_name", "input_reference_species",
                                                                          "atomic_pair_status", "_original_parent_id"])

    cache = _load_record_cache(record_cache_arg)
    cache_path_str = str(runtime_path(record_cache_arg)) if record_cache_arg else ""
    identity_reviews = load_identity_reviews(runtime_path(identity_review_arg)) if identity_review_arg else {}
    identity_review_requests: list[dict[str, str]] = []
    cached_tair = _cached_tair_index(cache)
    cached_symbols = _cached_gene_symbol_index(cache)
    cache_receipts = {}
    if record_cache_arg and runtime_path(record_cache_arg).is_dir():
        candidates = [runtime_path(record_cache_arg).parent / name for name in ("manifest.json", "record_cache_manifest.json")]
        manifests = [path for path in candidates if path.exists()]
        if len(manifests) > 1:
            raise ValueError("ambiguous reference record cache manifests")
        if manifests:
            cache_manifest = manifests[0]
            from .reference_record_cache import validate_reference_record_cache
            receipt_manifest = json.loads(cache_manifest.read_text(encoding="utf-8"))
            if receipt_manifest.get("query_ledger") or receipt_manifest.get("cache_complete") is False:
                _, receipts = validate_reference_record_cache(str(cache_manifest), mechanisms_arg=mechanisms_arg,
                                                              records_dir_arg=record_cache_arg)
                cache_receipts = {(item["reference_gene_name"], item["reference_species"]): item for item in receipts}

    resolved_rows: list[dict[str, object]] = []
    ledger_rows: list[dict[str, object]] = []
    review_rows: list[dict[str, str]] = []
    used_hypothesis_ids = {(str(row["case_id"]).casefold(), str(row["hypothesis_id"]).casefold())
                           for row in mechanisms.to_dict("records")}
    for record in mechanisms.to_dict("records"):
        hypothesis_id = str(record.get("hypothesis_id", ""))
        source_ids = str(record.get("source_ids", "")).strip()
        source_tokens = [value.strip() for value in re.split(r"[;,|]+", source_ids) if value.strip()]
        database_only = bool(source_tokens) and all(
            classify_citation(value)[0] == "database_accession" for value in source_tokens)
        source_accessions = {canonical_accession(value) for value in source_tokens
                             if re.match(r"^(?:uniprot(?:kb)?|tair)[:\s]", value, re.I)}
        source_accession = next(iter(source_accessions)) if len(source_accessions) == 1 else ""
        supplied_id = str(record.get("proposed_stable_id", "")).strip()
        supplied_id_malformed_status = ""
        if supplied_id and not (
            reference_handle(supplied_id, allow_untyped=True)
            or is_uniprot_accession(canonical_accession(supplied_id))
            or is_agi_locus(canonical_accession(supplied_id))
        ):
            supplied_id_malformed_status = (
                "malformed_proposed_stable_id_names_focal_species"
                if _names_focal_species(supplied_id, focal_species)
                else "malformed_proposed_stable_id_ignored"
            )
        named_reference_target = bool(
            str(record.get("reference_gene_name", "")).strip()
            and str(record.get("reference_species", "")).strip()
        )
        enrich_database_source = bool(
            database_only and len(source_accessions) == 1
            and all(canonical_accession(value) == source_accession for value in source_tokens)
            and (is_uniprot_accession(source_accession) or is_agi_locus(source_accession))
            and (not supplied_id or canonical_accession(supplied_id) == source_accession)
            and not named_reference_target)
        if source_ids and not enrich_database_source and not named_reference_target:
            resolved_rows.append(dict(record))  # byte-identical: never re-touch a sourced hypothesis
            ledger_rows.append({
                "hypothesis_id": hypothesis_id, "case_id": str(record.get("case_id", "")),
                "resolution_origin": "model_supplied",
                "reference_gene_name": str(record.get("reference_gene_name", "")),
                "reference_species": str(record.get("reference_species", "")),
                "input_proposed_stable_id": str(record.get("proposed_stable_id", "")).strip(),
                "input_source_ids": source_ids,
                "malformed_proposed_stable_id_status": "",
                "resolved_stable_id": str(record.get("proposed_stable_id", "")).strip(),
                "resolved_source_ids": source_ids,
                "accession_resolution_status": "model_supplied", "accession_resolution_query": "",
                "accession_resolution_hit_count": 0, "record_found": False,
                "function_source_backed": not database_only, "evaluation_status": "not_attempted",
                "source_provider": "", "source_url": "", "retrieved_at": "", "source_record_version": "",
                "provider_error": "", "review_reason": "",
                "retrieval_provenance": _provenance_json("", "", "model_supplied", 0, "", ""),
            })
            if database_only:
                ledger_rows[-1]["review_reason"] = "database_source_identity_unresolved"
                review_rows.append(_review_row(hypothesis_id, "database_source_identity_unresolved", len(source_accessions)))
            continue

        resolution_record = {**record, "proposed_stable_id": source_accession} if enrich_database_source else record

        # A gene-name query can return several exact provider-linked records.
        # Preserve every alternative as a derived row; do not select a winner.
        candidates: list[str] = []
        candidate_query = ""
        candidate_status = ""
        candidate_error = ""
        receipt = None
        candidate_hits = 0
        lookup_failures: list[str] = []
        query_urls: list[str] = []
        record_fetch_failures: dict = {}
        record_identity_failures: dict = {}
        lookup_search_state = "not_evaluated"
        # A named model hypothesis has two identity axes.  Resolve the exact
        # name/species pair first, even when the model also supplied a
        # well-formed-looking accession: the latter is an input claim, not
        # authoritative evidence that it names the declared gene.  Keep the
        # accession-only path below for rows without a complete name/species
        # pair, and retain the original value in the audit ledger.
        identity_status = "name_species_unresolved" if named_reference_target else "accession_only"
        identity_basis = "exact_name_species_lookup" if named_reference_target else ""
        if not enrich_database_source and named_reference_target:
            declared_name = str(record.get("reference_gene_name", "")).strip()
            locus = declared_name.upper()
            symbol = locus if is_agi_locus(locus) else _searchable_symbol(declared_name)
            genus = str(record.get("reference_species", "")).split(" ")[0].casefold()
            if record["atomic_pair_status"] not in {"atomic_name_species", "incomplete_name_species"}:
                candidate_query = "atomic_pairing:unresolved"
                candidate_status = str(record["atomic_pair_status"])
                candidate_error = candidate_status
            elif allow_external_api:
                try:
                    candidates, candidate_status, candidate_query, candidate_hits, failures, query_urls = _resolve_reference_candidates(
                        str(record.get("reference_gene_name", "")), str(record.get("reference_species", "")),
                        [record], timeout)
                    lookup_failures = failures
                    lookup_search_state = "partial_not_evaluated" if failures else "complete"
                    if failures:
                        if not candidates:
                            candidate_error = "; ".join(failures)
                except (RuntimeError, TimeoutError, OSError, ValueError) as error:
                    candidates = []
                    candidate_status = "provider_unavailable"
                    candidate_error = f"{type(error).__name__}: {error}"
                    lookup_failures = [candidate_error]
                    lookup_search_state = "not_evaluated"
            elif genus:
                # Offline exact record names also include ORF/locus identifiers
                # that intentionally are not admitted as general gene symbols.
                lookup_names = reference_identity_names(
                    [record], str(record.get("reference_gene_name", "")),
                    str(record.get("reference_species", "")))
                candidates = sorted({accession for name in lookup_names
                                     for accession in exact_record_candidates(
                                         cache, name, str(record.get("reference_species", "")))})
                candidate_hits = len(candidates)
                candidate_status = ("ambiguous_gene_name" if candidate_hits > 1 else
                                    "resolved_from_gene_name_cache" if candidate_hits == 1 else
                                    "not_attempted_offline")
                candidate_query = f"record_cache:gene:{symbol} AND organism_name:{str(record.get('reference_species', '')).strip()}"
                lookup_search_state = "complete"
                if not candidates:
                    # No exact name/species record: the worker's own declared
                    # accession, when cached, is graded record by record
                    # (provisional if only the name is unconfirmed).
                    declared_cached = sorted(accession for accession in declared_uniprot_accessions([record])
                                             if accession in cache)
                    if declared_cached:
                        candidates = declared_cached
                        candidate_hits = len(candidates)
                        candidate_status = "resolved_from_declared_accession_cache"
            receipt = cache_receipts.get((str(record["reference_gene_name"]), str(record["reference_species"])))
            if receipt and not allow_external_api:
                candidate_query = str(receipt.get("query", "")) or candidate_query
                (lookup_search_state, query_urls, lookup_failures,
                 record_fetch_failures, record_identity_failures) = _cache_lookup_receipt(receipt)
                if not candidates:
                    if lookup_failures:
                        candidate_status = "provider_unavailable"
                        candidate_error = "; ".join(lookup_failures)
                    elif record_fetch_failures:
                        candidate_status = "provider_unavailable"
                        candidate_error = json.dumps(record_fetch_failures, sort_keys=True)
                    elif record_identity_failures:
                        candidate_status = "provider_record_identity_mismatch"
                    elif receipt.get("status") == "not_evaluated":
                        candidate_status = "provider_unavailable"
                        candidate_error = str(receipt.get("reason", "cache_lookup_not_evaluated"))
                    else:
                        candidate_status = "unresolved_gene_name"
        if candidate_query:
            parent_id = hypothesis_id
            if not candidates:
                # A complete name/species declaration must not silently fall
                # back to an unrelated supplied accession.  Leave the input
                # accession visible only in the ledger and keep the resolved
                # hypothesis unresolved (or not_evaluated on provider error).
                parent_ledger = _resolve_one(
                    {**record, "proposed_stable_id": ""}, cache, cache_path_str,
                    cached_tair, cached_symbols, allow_external_api,
                    timeout, focal_species,
                    {"accession": "", "status": candidate_status or "unresolved_gene_name",
                     "query": candidate_query, "hits": candidate_hits,
                     "provider_error": candidate_error})
                parent_ledger["input_proposed_stable_id"] = supplied_id
                parent_ledger["input_source_ids"] = source_ids
                parent_ledger["malformed_proposed_stable_id_status"] = supplied_id_malformed_status
                parent_ledger["reference_identity_status"] = identity_status
                parent_ledger["reference_identity_match_basis"] = identity_basis
                parent_ledger["retrieval_provenance"] = _with_lookup_provenance(
                    parent_ledger["retrieval_provenance"], search_state=lookup_search_state,
                    query_urls=query_urls, failures=lookup_failures,
                    record_identity_state="mismatch" if record_identity_failures and not candidate_error else "unresolved",
                    record_fetch_failures=record_fetch_failures, record_identity_failures=record_identity_failures)
                ledger_rows.append(parent_ledger)
                unresolved = dict(record)
                unresolved["proposed_stable_id"] = ""
                unresolved["source_ids"] = ""
                unresolved["reference_identity_status"] = identity_status
                unresolved["reference_identity_match_basis"] = identity_basis
                unresolved["retrieval_provenance"] = parent_ledger["retrieval_provenance"]
                resolved_rows.append(unresolved)
                review_rows.append(_review_row(parent_id, str(parent_ledger["review_reason"]), candidate_hits))
                if supplied_id_malformed_status:
                    review_rows.append(_malformed_review_row(parent_id, supplied_id_malformed_status, supplied_id))
                continue
            # A declared accession is adopted provisionally only when no other
            # candidate is an exact name/species record for this target.
            declared_for_target = declared_uniprot_accessions([record])
            if declared_for_target & set(candidates) and set(candidates) - declared_for_target:
                exact_names = reference_identity_names([record], str(record.get("reference_gene_name", "")),
                                                       str(record.get("reference_species", "")))
                for other in sorted(set(candidates) - declared_for_target):
                    other_record = cache.get(other)
                    if other_record is None and allow_external_api:
                        try:
                            other_record = _fetch_uniprot(other, timeout)
                        except (RuntimeError, TimeoutError, OSError, ValueError):
                            other_record = None
                    if other_record and other in exact_record_candidates(
                            {other: other_record}, exact_names, str(record.get("reference_species", ""))):
                        declared_for_target = set()
                        break
            paper_binding_targets = []
            for derived_index, accession in enumerate(sorted(set(candidates)), 1):
                derived = dict(record)
                derived_id = (parent_id if len(set(candidates)) == 1 else
                              f"{parent_id}__ref_{canonical_accession(accession).casefold()}")
                derived_key = (str(record["case_id"]).casefold(), derived_id.casefold())
                if derived_key in used_hypothesis_ids and derived_id != parent_id:
                    raise ValueError(f"derived hypothesis id collides with input roster: {derived_id}")
                if derived_id != parent_id:
                    used_hypothesis_ids.add(derived_key)
                derived.update({"hypothesis_id": derived_id, "proposed_stable_id": accession,
                                "parent_hypothesis_id": parent_id,
                                "derived_hypothesis_index": str(derived_index),
                                "derived_hypothesis_count": str(len(set(candidates))),
                                "derived_accession": accession})
                child_ledger = _resolve_one(derived, cache, cache_path_str, cached_tair, cached_symbols,
                                            allow_external_api, timeout, focal_species,
                                            {"accession": accession, "status": candidate_status if candidate_hits == 1 else "resolved_from_gene_name",
                                             "query": candidate_query, "hits": candidate_hits,
                                             "identity_names": reference_identity_names(
                                                 [record], str(record.get("reference_gene_name", "")),
                                                 str(record.get("reference_species", ""))),
                                             "declared_accessions": sorted(declared_for_target),
                                             "identity_reviews": identity_reviews})
                graded_identity = str(child_ledger.get("reference_identity_status") or IDENTITY_EXACT)
                identity_basis_value = {
                    IDENTITY_EXACT: "exact_name_species_lookup",
                    IDENTITY_REVIEWED_ALIAS: "identity_review_ledger",
                    IDENTITY_DECLARED_ACCESSION: "declared_accession_species",
                }.get(graded_identity, "exact_name_species_lookup")
                child_ledger.update({"hypothesis_id": derived_id, "case_id": str(record.get("case_id", "")),
                                     "input_proposed_stable_id": supplied_id, "input_source_ids": source_ids,
                                     "malformed_proposed_stable_id_status": supplied_id_malformed_status,
                                     "accession_resolution_query": candidate_query,
                                     "accession_resolution_hit_count": candidate_hits,
                                     "parent_hypothesis_id": parent_id, "derived_accession": accession})
                child_ledger["reference_identity_status"] = graded_identity
                child_ledger["reference_identity_match_basis"] = identity_basis_value
                ledger_rows.append(child_ledger)
                resolved = dict(derived)
                resolved["parent_hypothesis_id"] = parent_id
                resolved["derived_hypothesis_index"] = str(derived_index)
                resolved["derived_hypothesis_count"] = str(len(set(candidates)))
                resolved["derived_accession"] = accession
                resolved["reference_identity_status"] = graded_identity
                resolved["reference_identity_match_basis"] = identity_basis_value
                if child_ledger["record_found"] and graded_identity == IDENTITY_DECLARED_ACCESSION:
                    # Provisional adoption: the worker's own accession and the
                    # species agree, the name does not. An independent identity
                    # reviewer may confirm it; until then it is review-only.
                    found_record = cache.get(accession) or {}
                    identity_review_requests.append({
                        "hypothesis_id": derived_id, "case_id": str(record.get("case_id", "")),
                        "reference_gene_name": str(record.get("reference_gene_name", "")),
                        "input_reference_gene_name": str(record.get("input_reference_gene_name", "")),
                        "reference_species": str(record.get("reference_species", "")),
                        "accession": accession,
                        "record_gene_names": ";".join(record_gene_names(found_record)) if found_record else "",
                        "record_organism": str((found_record.get("organism") or {}).get("scientificName", ""))
                        if found_record else "",
                        "source_ids": source_ids,
                        "question": "Does the paper's gene name refer to the gene on this UniProt record? "
                                    "Answer same_gene only with a citation that links the two names.",
                    })
                if not child_ledger["record_found"]:
                    child_ledger["reference_identity_status"] = resolved["reference_identity_status"] = "name_species_unresolved"
                    resolved["proposed_stable_id"] = ""
                provenance = json.loads(child_ledger["retrieval_provenance"] or "{}")
                provider_sources = [value for value in str(child_ledger["resolved_source_ids"]).split(";") if value]
                validated_handles = {tuple(value) for value in provenance.get("validated_reference_handles", [])}
                declared_handle = reference_handle(supplied_id, allow_untyped=True)
                same_record = bool(declared_handle and declared_handle in validated_handles)
                retain_original_sources = bool(
                    child_ledger["record_found"]
                    and resolved["reference_identity_status"] in PROJECTABLE_IDENTITIES
                    and supplied_id and not supplied_id_malformed_status
                    and (same_record or canonical_accession(supplied_id) == canonical_accession(accession))
                )
                # Reconfirming the same endpoint must not discard the original
                # claim's citations merely because the provider bibliography is
                # incomplete. Keep them as review inputs, never as provider
                # links. Changed/ambiguous endpoints retain audit provenance
                # only and cannot borrow the parent's source association.
                if retain_original_sources:
                    child_ledger["resolved_source_ids"] = ";".join(sorted(set(provider_sources) | set(source_tokens)))
                record_identity_state = (("resolved_exact_name_species"
                                          if graded_identity in {IDENTITY_EXACT, IDENTITY_REVIEWED_ALIAS}
                                          else "resolved_declared_accession_species")
                                         if child_ledger["record_found"] else
                                         "mismatch" if child_ledger["accession_resolution_status"] ==
                                         "provider_record_identity_mismatch" else "unresolved")
                candidate_record_failures = dict(record_fetch_failures)
                if child_ledger.get("provider_error"):
                    candidate_record_failures.setdefault(accession, str(child_ledger["provider_error"]))
                provenance.update(lookup_search_state=lookup_search_state, query_urls=query_urls,
                                 lookup_failures=lookup_failures,
                                 record_identity_state=record_identity_state)
                if candidate_record_failures:
                    provenance["record_fetch_failures"] = candidate_record_failures
                if record_identity_failures:
                    provenance["record_identity_failures"] = record_identity_failures
                provenance.update({
                    "parent_hypothesis_id": parent_id, "derived_accession": accession,
                    "candidate_query": candidate_query, "candidate_count": candidate_hits,
                    "child_hypothesis_id": derived_id,
                    "input_proposed_stable_id": supplied_id,
                    "original_source_ids": source_tokens,
                    "provider_linked_source_ids": provider_sources,
                    "validated_reference_handles": [list(value) for value in sorted(validated_handles)],
                    "validated_reference_accession": accession if child_ledger["record_found"] else "",
                    "original_source_retention_state": (
                        "same_record_review_material" if retain_original_sources
                        else "audit_only_identity_not_confirmed")})
                if child_ledger["record_found"] and resolved.get("paper_reports_json"):
                    paper_binding_targets.append((resolved, child_ledger, provenance, bool(declared_handle and not same_record)))
                resolved["retrieval_provenance"] = json.dumps(provenance, sort_keys=True, separators=(",", ":"))
                child_ledger["retrieval_provenance"] = resolved["retrieval_provenance"]
                if child_ledger["resolution_origin"] == "provider_resolved" or retain_original_sources:
                    resolved["source_ids"] = child_ledger["resolved_source_ids"]
                    resolved["retrieval_provenance"] = child_ledger["retrieval_provenance"]
                else:
                    resolved["source_ids"] = ""
                resolved_rows.append(resolved)
                if child_ledger["resolution_origin"] != "provider_resolved":
                    review_rows.append(_review_row(derived_id, str(child_ledger["review_reason"]),
                                                   int(child_ledger["accession_resolution_hit_count"])))
                if supplied_id_malformed_status:
                    review_rows.append(_malformed_review_row(derived_id, supplied_id_malformed_status, supplied_id))
            # An explicit external ID may distinguish same-name provider
            # records. Wait for every sibling record before transporting it;
            # a shared handle or an unavailable sibling remains ambiguous.
            handle_owners: dict[tuple[str, str], set[str]] = {}
            for target, _, target_provenance, _ in paper_binding_targets:
                for handle in target_provenance.get("validated_reference_handles", []):
                    handle_owners.setdefault(tuple(handle), set()).add(target["proposed_stable_id"])
            unique_handles = ({handle for handle, owners in handle_owners.items() if len(owners) == 1}
                              if len(paper_binding_targets) == len(set(candidates)) else set())
            from .reference_shared_gene import non_unique_insdc_handles, shared_gene_proofs
            unique_handles -= non_unique_insdc_handles(cache_path_str, handle_owners)
            shared_proofs = shared_gene_proofs(cache_path_str, set(candidates), str(record["reference_species"]),
                handle_owners, complete=(len(paper_binding_targets) == len(set(candidates))
                    and lookup_search_state == "complete" and not lookup_failures
                    and not record_fetch_failures and not record_identity_failures))
            for target, target_ledger, target_provenance, conflict in paper_binding_targets:
                reports = json.loads(target["paper_reports_json"])
                rebindings = []
                used_shared_proofs = {}
                for report in reports:
                    anchor = str(report.get("reference_gene_id", ""))
                    basis = _resolved_paper_anchor_basis(
                        report, record, target_provenance, unique_reference=len(set(candidates)) == 1,
                        unique_handles=unique_handles, supplied_id_conflict=conflict, shared_handles=set(shared_proofs))
                    if basis:
                        report["reference_gene_id"] = target["proposed_stable_id"]
                        rebindings.append({"input_reference_gene_id": anchor, "reference_gene_id": target["proposed_stable_id"],
                                           "basis": basis})
                        if basis.startswith("provider_record_shared_gene_"):
                            handle = reference_handle(anchor, allow_untyped=True)
                            used_shared_proofs[handle] = shared_proofs[handle]
                if rebindings:
                    target["paper_reports_json"] = json.dumps(reports, sort_keys=True, separators=(",", ":"))
                    target_provenance["paper_reference_rebindings"] = rebindings
                    if used_shared_proofs:
                        target_provenance["paper_shared_reference_proofs"] = [used_shared_proofs[key] for key in sorted(used_shared_proofs)]
                    target["retrieval_provenance"] = json.dumps(target_provenance, sort_keys=True, separators=(",", ":"))
                    target_ledger["retrieval_provenance"] = target["retrieval_provenance"]
            if receipt and not allow_external_api and receipt["status"] == "not_evaluated":
                pending_id = parent_id + "__pending"
                pending_key = (str(record["case_id"]).casefold(), pending_id.casefold())
                if pending_key in used_hypothesis_ids:
                    raise ValueError(f"derived hypothesis id collides with input roster: {pending_id}")
                used_hypothesis_ids.add(pending_key)
                pending = {**record, "hypothesis_id": pending_id, "proposed_stable_id": "", "source_ids": ""}
                identity_only = bool(record_identity_failures) and not (lookup_failures or record_fetch_failures)
                pending_errors = {f"lookup:{index}": failure for index, failure in enumerate(lookup_failures)}
                pending_errors.update(record_fetch_failures)
                pending_ledger = _resolve_one(pending, cache, cache_path_str, cached_tair, cached_symbols, False, timeout, focal_species,
                    {"accession": "", "status": "provider_record_identity_mismatch" if identity_only else "provider_unavailable",
                     "query": candidate_query, "hits": receipt["provider_hit_count"],
                     "provider_error": "" if identity_only else
                         json.dumps(pending_errors, sort_keys=True) if pending_errors else receipt["reason"]})
                pending_ledger.update(input_proposed_stable_id=supplied_id, input_source_ids=source_ids,
                                      reference_identity_status="name_species_unresolved",
                                      reference_identity_match_basis="exact_name_species_lookup")
                pending.update(reference_identity_status="name_species_unresolved", reference_identity_match_basis="exact_name_species_lookup")
                pending_ledger["retrieval_provenance"] = _with_lookup_provenance(
                    pending_ledger["retrieval_provenance"], search_state=lookup_search_state,
                    query_urls=query_urls, failures=lookup_failures,
                    record_identity_state="mismatch" if identity_only else "unresolved",
                    record_fetch_failures=record_fetch_failures, record_identity_failures=record_identity_failures)
                pending["retrieval_provenance"] = pending_ledger["retrieval_provenance"]
                ledger_rows.append(pending_ledger)
                resolved_rows.append(pending)
                review_rows.append(_review_row(pending_id, str(pending_ledger["review_reason"]), int(receipt["provider_hit_count"])))
            continue

        ledger = _resolve_one(resolution_record, cache, cache_path_str, cached_tair, cached_symbols,
                              allow_external_api, timeout, focal_species)
        ledger["reference_identity_status"] = identity_status
        ledger["reference_identity_match_basis"] = identity_basis
        if enrich_database_source:
            ledger["input_source_ids"] = source_ids
            ledger["input_proposed_stable_id"] = supplied_id
            if ledger["resolution_origin"] == "provider_resolved":
                ledger["resolved_source_ids"] = ";".join(dict.fromkeys(
                    [*source_tokens, *str(ledger["resolved_source_ids"]).split(";")]))
        ledger_rows.append(ledger)
        origin = str(ledger["resolution_origin"])
        malformed_status = str(ledger["malformed_proposed_stable_id_status"])

        if origin == "provider_resolved":
            augmented = dict(record)
            if not str(augmented.get("proposed_stable_id", "")).strip() or malformed_status:
                augmented["proposed_stable_id"] = ledger["resolved_stable_id"]
            augmented["source_ids"] = ledger["resolved_source_ids"]
            augmented["retrieval_provenance"] = ledger["retrieval_provenance"]
            augmented["reference_identity_status"] = identity_status
            augmented["reference_identity_match_basis"] = identity_basis
            resolved_rows.append(augmented)
        elif malformed_status:
            # Never forward a fabricated/malformed identifier downstream, even
            # when the gene-name fallback could not resolve it either.
            augmented = dict(record)
            augmented["proposed_stable_id"] = ""
            augmented["retrieval_provenance"] = ledger["retrieval_provenance"]
            augmented["reference_identity_status"] = identity_status
            augmented["reference_identity_match_basis"] = identity_basis
            resolved_rows.append(augmented)
        else:
            unchanged = dict(record)
            unchanged["reference_identity_status"] = identity_status
            unchanged["reference_identity_match_basis"] = identity_basis
            resolved_rows.append(unchanged)

        if malformed_status:
            review_rows.append(_malformed_review_row(hypothesis_id, malformed_status,
                                                      str(ledger["input_proposed_stable_id"])))
        if origin != "provider_resolved":
            review_rows.append(_review_row(hypothesis_id, str(ledger["review_reason"]),
                                           int(ledger["accession_resolution_hit_count"])))

    expanded_by_id = {(str(row["case_id"]), str(row["hypothesis_id"])): row for row in expanded}
    for output_row in ledger_rows + resolved_rows:
        row_id = str(output_row["hypothesis_id"])
        atomic_parent = str(output_row.get("parent_hypothesis_id", "")) or row_id.removesuffix("__pending")
        case = str(output_row["case_id"])
        original = expanded_by_id.get((case, atomic_parent)) or expanded_by_id.get((case, row_id))
        if original:
            current_sources = {value.strip() for value in re.split(r"[;,|]+", str(output_row.get("source_ids", ""))) if value.strip()}
            original_sources = {value.strip() for value in re.split(r"[;,|]+", str(original.get("source_ids", ""))) if value.strip()}
            output_row["review_only_source_ids"] = ";".join(sorted(original_sources - current_sources))
            for key in ("input_reference_gene_name", "input_reference_species", "atomic_pair_status"):
                output_row[key] = original[key]
            if str(original["_original_parent_id"]) != row_id or output_row.get("parent_hypothesis_id"):
                output_row["parent_hypothesis_id"] = original["_original_parent_id"]
            if original["atomic_pair_status"] == "incomplete_name_species":
                continue
            provenance = json.loads(str(output_row.get("retrieval_provenance", "")) or "{}")
            provenance.update({key: original[key] for key in ("input_reference_gene_name", "input_reference_species", "atomic_pair_status")})
            provenance.update(atomic_reference_gene_name=original["reference_gene_name"], atomic_reference_species=original["reference_species"],
                              parent_hypothesis_id=original["_original_parent_id"])
            output_row["retrieval_provenance"] = json.dumps(provenance, sort_keys=True, separators=(",", ":"))
    resolved = _frame(resolved_rows, RESOLVED_MECHANISM_COLUMNS)
    ledger_frame = _frame(ledger_rows, REFERENCE_GENE_RESOLUTION_COLUMNS)
    reviews = (_frame(review_rows, REVIEW_COLUMNS).drop_duplicates()
               .sort_values(["reason", "target"], kind="stable") if review_rows else _frame([], REVIEW_COLUMNS))

    paths = {
        "mechanisms_resolved": out / "mechanism_hypotheses_resolved.tsv",
        "ledger": out / "reference_gene_resolution.tsv",
        "review_required": out / "review_required.tsv",
    }
    _write(resolved, paths["mechanisms_resolved"])
    request_columns = ["hypothesis_id", "case_id", "reference_gene_name", "input_reference_gene_name",
                       "reference_species", "accession", "record_gene_names", "record_organism",
                       "source_ids", "question"]
    paths["identity_review_requests"] = out / "identity_review_requests.tsv"
    _write(_frame(identity_review_requests, request_columns), paths["identity_review_requests"])
    for request in identity_review_requests:
        review_rows_extra = {
            "reason": "reference_name_unconfirmed_on_declared_accession",
            "target": request["hypothesis_id"],
            "required_next_action": ("provisional: the declared accession and species match but the gene name is not "
                                     "on the record; an independent identity reviewer may confirm same_gene with a "
                                     "citation (identity_review_requests.tsv), otherwise it stays review-only"),
            "allowed_choices": "confirm same_gene with citation | keep provisional | replace the reference",
        }
        reviews = pd.concat([reviews, pd.DataFrame([review_rows_extra], columns=REVIEW_COLUMNS)], ignore_index=True)
    if identity_review_requests:
        reviews = reviews.drop_duplicates().sort_values(["reason", "target"], kind="stable")
    _write(ledger_frame, paths["ledger"])
    _write(reviews, paths["review_required"])

    origin_counts = {origin: int((ledger_frame["resolution_origin"] == origin).sum()) if len(ledger_frame) else 0
                     for origin in RESOLUTION_ORIGINS}
    evaluation_counts = {status: int((ledger_frame["evaluation_status"] == status).sum()) if len(ledger_frame) else 0
                         for status in EVALUATION_STATUSES}
    not_evaluated_ids = (sorted(ledger_frame.loc[ledger_frame["evaluation_status"] == "not_evaluated", "hypothesis_id"])
                         if len(ledger_frame) else [])
    malformed_counts = ({status: int((ledger_frame["malformed_proposed_stable_id_status"] == status).sum())
                        for status in MALFORMED_STATUSES if status}
                        if len(ledger_frame) else {status: 0 for status in MALFORMED_STATUSES if status})
    malformed_focal_species_hypothesis_ids = (
        sorted(ledger_frame.loc[ledger_frame["malformed_proposed_stable_id_status"]
                                == "malformed_proposed_stable_id_names_focal_species", "hypothesis_id"])
        if len(ledger_frame) else []
    )

    manifest = {
        "schema_version": SCHEMA_VERSION, "command": "reference-gene-resolve",
        "external_api_allowed": allow_external_api,
        "provider": "uniprot_rest" if allow_external_api else ("local_record_cache" if cache else "none"),
        "record_cache": cache_path_str or None,
        "timeout_seconds": timeout,
        "focal_species": focal_species or None,
        "mechanisms_input": file_fingerprint(mechanisms_path, False),
        "hypothesis_count": int(original_count), "atomic_hypothesis_count": int(len(mechanisms)),
        "resolution_origin_counts": origin_counts,
        "evaluation_status_counts": evaluation_counts,
        "malformed_proposed_stable_id_counts": malformed_counts,
        "malformed_proposed_stable_id_names_focal_species_hypothesis_ids": malformed_focal_species_hypothesis_ids,
        "review_required_count": int(len(reviews)),
        "identity_review_request_count": len(identity_review_requests),
        "identity_review_ledger": file_fingerprint(runtime_path(identity_review_arg)) if identity_review_arg else None,
        "not_evaluated_hypothesis_ids": not_evaluated_ids,
        "gates": [
            "model-supplied literature/unknown sources are preserved; database-only IDs may retrieve linked publications from one unambiguous matching identity without replacing original IDs",
            "a supplied, well-formed proposed_stable_id (UniProt accession or AGI locus) is used as-is to "
            "fetch citations; it is never re-derived from the gene name",
            "a supplied proposed_stable_id that is not a recognisable UniProt accession or AGI locus is never "
            "queried; it is treated as absent, resolution falls back to reference_gene_name + reference_species, "
            "and the original string is preserved verbatim in the ledger with its own status",
            "a malformed proposed_stable_id that names the focal species gets its own distinct ledger status, "
            "separate from ordinary malformation",
            "a gene name resolving to multiple exact provider accessions is fanned out deterministically and never tie-broken",
            "original citations absent from resolved source_ids remain in review_only_source_ids with parent lineage; cached text is review context only, never automatic support for an alternative or unresolved accession",
            "offline, a gene symbol is resolved only through an index built from the cache's own "
            "genes[].geneName.value / genes[].synonyms[].value fields matched to the record's scientific species for bare-name expansion -- "
            "a cache is an authoritative record snapshot, not an alias guess",
            "a provider failure (network/timeout/HTTP) is not_evaluated, never a rejection or a negative result",
            "without --allow-external-api and without a cache hit, nothing is attempted: not_attempted_offline",
            "a retrieved record without a source-backed function citation fills nothing",
            "mechanism_hypotheses_resolved.tsv never carries a malformed proposed_stable_id forward into "
            "hypothesis-verify: it is always cleaned to a provider-resolved accession or empty",
            "hypothesis_verify's accession/species/sequence/source verification gate is unchanged by this stage",
            "a record naming an infraspecific taxon (subsp./var./cv.) of the declared binomial matches that species; "
            "a declared subspecies is never confirmed by a species-level record",
            "a worker-declared UniProt accession whose record matches the species but not the gene name is adopted "
            "provisionally (declared_accession_species_name_unconfirmed): review-only, never verified support; "
            "an identity review row (same_gene, citation, reviewer, a gene name present on the record) confirms it",
        ],
        "outputs": {name: str(path) for name, path in sorted(paths.items())},
    }
    manifest_path = out / "reference_gene_resolution_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    return {
        "manifest": str(manifest_path),
        "mechanisms_resolved": str(paths["mechanisms_resolved"]),
        "ledger": str(paths["ledger"]),
        "review_required": str(paths["review_required"]),
        "hypothesis_count": int(len(mechanisms)),
        "provider_resolved_count": origin_counts["provider_resolved"],
        "model_supplied_count": origin_counts["model_supplied"],
        "not_resolved_count": origin_counts["not_resolved"],
        "not_evaluated_count": evaluation_counts["not_evaluated"],
        "malformed_proposed_stable_id_count": sum(malformed_counts.values()),
        "review_required_count": int(len(reviews)),
        "identity_review_requests": str(paths["identity_review_requests"]),
        "identity_review_request_count": len(identity_review_requests),
    }
