"""Cache-first citation identity checks; deliberately no claim/content inference."""
from __future__ import annotations

from .citation_tokens import split_citations

import hashlib
import html
import json
import re
import unicodedata
import urllib.parse
import urllib.request
from datetime import datetime, timezone

import pandas as pd

from .core import runtime_path
from .identifiers import canonical_accession, is_agi_locus, is_uniprot_accession
from .gene_metadata import GeneMetadataChecks
from .paper_reports import PAPER_FIELD, PAPER_HANDOFF_POLICY, normalize_paper_reports, source_tokens, select_bibliography_papers

INITIAL_POLICY = "initial_bibliography_v1"
CONTENT_POLICY = "content_v1"
RESOLVED_STATES = {"id_title_matched", "id_resolved_title_not_compared"}
PMC_ID_CONVERTER = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"


def normalize_title(value: object) -> str:
    text = html.unescape(str(value or ""))
    text = re.sub(r"</?[A-Za-z][^>]*>", "", text)
    text = unicodedata.normalize("NFKC", text).casefold()
    return " ".join("".join(c if c.isalnum() else " " for c in text).split())


def fetch_metadata(kind: str, identifier: str) -> dict | None:
    """Europe PMC lite only. Never silently fall back to abstract/full-text."""
    from .source_verify import EUROPE_PMC_SEARCH, USER_AGENT, _query_for
    query = urllib.parse.urlencode({"query": _query_for(kind, identifier),
                                   "format": "json", "pageSize": 2, "resultType": "lite"})
    request = urllib.request.Request(EUROPE_PMC_SEARCH + "?" + query,
                                     headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            payload = json.load(response)
        results = payload["resultList"]["result"]
        if not isinstance(results, list) or any(not isinstance(row, dict) for row in results):
            raise ValueError("invalid metadata response")
        # An exact identifier query must resolve uniquely; never trust first hit.
        matches = [row for row in results if identifier_matches(kind, identifier, row)]
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            raise ValueError("ambiguous metadata identifier")
        return results[0] if results else None
    except (OSError, ValueError, KeyError, TypeError) as error:
        raise RuntimeError(f"bibliography provider unavailable: {error}") from error


def identifier_matches(kind: str, identifier: str, result: dict) -> bool:
    from .source_verify import classify_citation
    actual = result.get(kind) or (result.get("id") if kind == "pmid" and result.get("source", "MED") == "MED" else "")
    token = str(actual) if classify_citation(str(actual))[0] == kind else f"{kind}:{actual}"
    return bool(actual) and source_tokens([f"{kind}:{identifier}"]) == source_tokens([token])


def fetch_pmc_mapping(identifier: str) -> tuple[dict | None, dict]:
    """Resolve a missing PMC metadata alias through official ID metadata only."""
    from .source_verify import USER_AGENT, classify_citation
    url = PMC_ID_CONVERTER + "?" + urllib.parse.urlencode({"ids": identifier, "idtype": "pmcid", "format": "json", "tool": "candidate-qtl"})
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": USER_AGENT}), timeout=30) as response:
            raw = response.read()
        payload = json.loads(raw)
        records = payload["records"]
        if payload.get("status") != "ok" or not isinstance(records, list) or any(not isinstance(row, dict) for row in records):
            raise ValueError("invalid PMC ID Converter response")
        provenance = {"provider": "ncbi_pmc_id_converter", "url": url,
                      "retrieved_at": datetime.now(timezone.utc).isoformat(),
                      "response_sha256": hashlib.sha256(raw).hexdigest(), "requested_pmcid": identifier}
        if not records or (len(records) == 1 and records[0].get("status") == "error"):
            return None, provenance
        if len(records) != 1 or not identifier_matches("pmcid", identifier, records[0]):
            raise ValueError("ambiguous or mismatched PMC ID Converter mapping")
        mapping = {"pmcid": identifier}
        for kind in ("pmid", "doi"):
            if records[0].get(kind):
                parsed = classify_citation(f"{kind}:{records[0][kind]}")
                if parsed[0] != kind:
                    raise ValueError("invalid PMC ID Converter alias")
                mapping[kind] = parsed[1]
        if len(mapping) == 1:
            return None, provenance
        return mapping, {**provenance, "mapping": mapping}
    except (OSError, ValueError, KeyError, TypeError) as error:
        raise RuntimeError(f"bibliography provider unavailable: ncbi_pmc_id_converter: {error}") from error


def concrete_bridge(row: dict) -> bool:
    """Shape/generic-label guard, not an independent biological adjudication.

    Original A/B owns scientific relevance. Text length or a real PMID is not a
    biological score. Keep specific functional observations and a stated bridge.
    """
    function = str(row.get("biological_function") or row.get("expected_role") or "").strip()
    bridge = str(row.get("phenotype_relation") or row.get("mechanistic_link") or row.get("mechanism_statement") or "").strip()
    generic = {"kinase", "protein kinase", "transcription factor", "tf", "unknown",
               "unknown function", "protein", "gene", "go", "housekeeping"}
    return bool(function and bridge and normalize_title(function) not in generic
                and normalize_title(bridge) not in {"unknown", "possibly relevant", "may affect phenotype"})


def _declared_input_accession(provenance: dict) -> str:
    """The researcher's declared accession, or "" when none was declared.

    Reference resolution already ignores a malformed proposed_stable_id (for
    example "unknown (accession not verified)") and records it only for audit.
    Reading that audit text back as a declared endpoint would make every
    resolved record look like a different accession.
    """
    from .reference_name_pairs import reference_handle
    supplied = str(provenance.get("input_proposed_stable_id", "")).strip()
    if not supplied:
        return ""
    bare = canonical_accession(supplied)
    if (reference_handle(supplied, allow_untyped=True)
            or is_uniprot_accession(bare) or is_agi_locus(bare)):
        return supplied
    return ""


def _rebound_papers(row: dict, provenance: dict, reports: list[dict]) -> list[dict]:
    """Consume resolved paper endpoints without restoring the parent bibliography."""
    from .reference_name_pairs import reference_handle
    current = reference_handle(row.get("proposed_stable_id", ""), allow_untyped=True)
    if (_declared_input_accession(provenance) or not current
            or provenance.get("record_identity_state") not in {"resolved_exact_name_species",
                                                               "resolved_declared_accession_species"}
            or current != reference_handle(provenance.get("validated_reference_accession", ""), allow_untyped=True)):
        return []
    receipts = provenance.get("paper_reference_rebindings", [])
    if not isinstance(receipts, list) or not any(
        isinstance(receipt, dict) and receipt.get("basis") == "provider_record_cross_reference"
        and reference_handle(receipt.get("reference_gene_id", ""), allow_untyped=True) == current
        and list(reference_handle(receipt.get("input_reference_gene_id", ""), allow_untyped=True) or ())
        in provenance.get("validated_reference_handles", []) for receipt in receipts
    ):
        return []
    normalize = lambda value: " ".join(str(value).split()).casefold()
    return [report for report in reports
            if reference_handle(report["reference_gene_id"], allow_untyped=True) == current
            and (report.get("gene_correspondence") == "homolog" or not report.get("reported_species")
                 or normalize(report["reported_species"]) == normalize(row.get("reference_species", "")))]


def verify_bibliography(mechanisms_arg: str, out_dir_arg: str, *, cache_arg: str | None,
                        allow_external_api: bool, id_column: str, source_column: str, focal_species: str = "") -> dict:
    from .source_verify import _cache_key, _load_cache, classify_citation
    frame = pd.read_csv(runtime_path(mechanisms_arg), sep="\t", dtype=str).fillna("")
    if not {id_column, source_column} <= set(frame) or frame[id_column].duplicated().any():
        raise ValueError("bibliography input requires unique IDs and source column")
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    cache = _load_cache(cache_arg)
    lookups: dict[tuple[str, str], tuple[dict | None, str]] = {}
    lookup_provenance: dict[tuple[str, str], list[dict]] = {}

    def lookup(kind: str, identifier: str) -> tuple[dict | None, str]:
        key = (kind, identifier)
        if key in lookups:
            return lookups[key]
        record = cache.get(_cache_key(kind, identifier))
        if record is not None:
            lookups[key] = (record, "local_record_cache")
        elif not allow_external_api:
            lookups[key] = (None, "not_attempted_offline")
        else:
            try:
                record, provider = fetch_metadata(kind, identifier), "europe_pmc_lite"
                if record is None and kind == "pmcid":
                    mapping, provenance = fetch_pmc_mapping(identifier)
                    lookup_provenance[key] = [provenance]
                    provider += "+ncbi_pmc_id_converter"
                    aliases = [alias for alias in ("pmid", "doi") if mapping and alias in mapping]
                    aliases.sort(key=lambda alias: not (
                        _cache_key(alias, mapping[alias]) in cache or lookups.get((alias, mapping[alias]), (None, ""))[0] is not None))
                    failures = []
                    for alias in aliases:
                        candidate, metadata_provider = lookup(alias, mapping[alias])
                        provider = "europe_pmc_lite+ncbi_pmc_id_converter+" + metadata_provider
                        if candidate is None:
                            if metadata_provider.startswith("bibliography provider unavailable"):
                                failures.append(metadata_provider)
                            continue
                        if not identifier_matches(alias, mapping[alias], candidate) or any(
                            candidate.get(field) and not identifier_matches(field, value, candidate)
                            for field, value in mapping.items()
                        ):
                            raise ValueError("PMC mapping conflicts with title metadata")
                        if str(candidate.get("title", "")).strip():
                            if failures:
                                provenance = {**provenance, "metadata_lookup_failures": failures}
                                lookup_provenance[key] = [provenance]
                            record = {**candidate, **mapping, "metadata_provenance": [*candidate.get("metadata_provenance", []), provenance]}
                            for field, value in mapping.items():
                                if lookups.get((field, value), (None, ""))[0] is None:
                                    lookups[(field, value)] = (record, provider)
                            break
                    if record is None and failures:
                        raise RuntimeError("; ".join(failures))
                lookups[key] = (record, provider)
            except (RuntimeError, ValueError) as failure:
                message = str(failure)
                lookups[key] = (None, message if message.startswith("bibliography provider unavailable") else
                                f"bibliography provider unavailable: {message}")
        return lookups[key]

    metadata: dict[str, dict] = {}
    gene_checks = GeneMetadataChecks(cache_arg, allow_external_api)
    gene_rows = []
    rows, annotated, selections = [], [], []
    for original in frame.to_dict("records"):
        paper_reports = normalize_paper_reports(original.get(PAPER_FIELD))
        delivery_reports = json.loads(paper_reports or "[]")
        paper_scoped = False
        citations = split_citations(original[source_column])
        # Resolution enriches source_ids with the provider's whole bibliography.
        # Those are retrieval leads, not additional sources for this hypothesis.
        provenance = original.get("retrieval_provenance", "")
        try:
            provenance = json.loads(provenance) if provenance else {}
        except (ValueError, TypeError):
            provenance = {}
        if isinstance(provenance, dict):
            gene_checks.seed_record(provenance.get("validated_reference_metadata"))
        endpoint_changed = False
        provider_sources = ""
        if isinstance(provenance, dict) and isinstance(provenance.get("original_source_ids"), list):
            provider_sources = ";".join(provenance.get("provider_linked_source_ids", []))
            declared_input = _declared_input_accession(provenance)
            declared = canonical_accession(declared_input)
            current = canonical_accession(str(original.get("proposed_stable_id", "")))
            same_endpoint = bool(current and (
                current == declared if declared else provenance.get("candidate_count") == 1))
            from .reference_name_pairs import reference_handle
            declared_handle = reference_handle(declared_input, allow_untyped=True) if declared_input else None
            same_endpoint = same_endpoint or bool(
                declared_handle and list(declared_handle) in provenance.get("validated_reference_handles", [])
                and current == provenance.get("validated_reference_accession")
            )
            identity_checked = provenance.get("record_identity_state") in {
                "resolved_exact_name_species", "resolved_declared_accession_species"}
            if same_endpoint and identity_checked:
                citations = [str(value) for value in provenance["original_source_ids"] if value]
            else:
                # A different accession must never borrow the declared target's
                # citations/title merely because it shares a gene-name lookup.
                delivery_reports = _rebound_papers(original, provenance, delivery_reports)
                paper_ids = set().union(*(source_tokens([report["source_id"]])
                                         | source_tokens(report["source_identifier_aliases"]) for report in delivery_reports))
                citations = [str(value) for value in provenance["original_source_ids"]
                             if source_tokens([value]) & paper_ids]
                paper_scoped = True
                endpoint_changed = not citations
        claim_citations = ";".join(citations)
        aliases_supplied = "" if endpoint_changed or paper_scoped else str(original.get("source_identifier_aliases", "")).strip()
        alias_citations = set(split_citations(aliases_supplied))
        if aliases_supplied:
            alias_citations.update(citations[:1])
        title_citations = alias_citations or set(citations[:1])
        citations += [s for s in sorted(alias_citations) if s not in citations]
        paper_titles, paper_alias_groups = {}, []
        if paper_reports and not endpoint_changed:
            # Each paper owns its title/aliases. Never attach the first paper's
            # title to another paper or restore a rejected reference endpoint.
            for report in delivery_reports:
                ids = source_tokens([report["source_id"]]) | source_tokens(report["source_identifier_aliases"])
                if not ids & source_tokens(claim_citations):
                    continue
                for identifier in ids:
                    paper_titles.setdefault(identifier, set()).add(report["source_title"])
                paper_alias_groups.append(ids)
                citations += [s for s in [report["source_id"], *split_citations(report["source_identifier_aliases"])]
                              if s.strip() and not source_tokens([s]) & source_tokens(";".join(citations))]
        selection = select_bibliography_papers({**original, PAPER_FIELD: delivery_reports}, citations)
        selections.append({"hypothesis_id": original[id_column], **selection})
        selected_ids = source_tokens(selection["selected_source_ids"])
        report_checks = [gene_checks.check(report) for report in delivery_reports
                         if (source_tokens([report["source_id"]]) | source_tokens(report["source_identifier_aliases"])) & selected_ids]
        gene_rows.extend({"hypothesis_id": original[id_column], **check} for check in report_checks)
        citations = [citation for citation in dict.fromkeys(citations)
                     if source_tokens([citation]) & selected_ids or classify_citation(citation)[0] == "database_accession"]
        group = []
        for citation in citations:
            kind, identifier = classify_citation(citation)
            state, result, error = "unrecognized_identifier", None, ""
            provider = "none"
            if kind in {"pmid", "pmcid", "doi"}:
                result, provider = lookup(kind, identifier)
                if result is None:
                    state = ("not_attempted" if provider == "not_attempted_offline" else
                             "provider_failure" if provider.startswith("bibliography provider unavailable") else "not_found")
                    error = provider if state == "provider_failure" else ""
                elif not identifier_matches(kind, identifier, result):
                    state = "identifier_mismatch"
                else:
                    titles = (set().union(*(paper_titles.get(key, set()) for key in source_tokens([citation]))) if paper_reports else
                              {(original.get("source_title") or original.get("proposed_title") or "")
                               if citation in title_citations else ""})
                    titles = {normalize_title(title) for title in titles if str(title).strip()}
                    state = ("id_title_matched" if titles == {normalize_title(result.get("title"))} else
                             "title_mismatch") if titles else "id_resolved_title_not_compared"
                    # Save metadata only, even when the supplied cache has text.
                    metadata[citation] = {k: result[k] for k in
                        ("id", "source", "pmid", "pmcid", "doi", "title", "pubYear", "journalTitle", "metadata_provenance") if k in result}
                    # Reuse one paper's aliases for subsequent claims in this call.
                    for alias in ("pmid", "pmcid", "doi"):
                        if result.get(alias):
                            raw_alias = str(result[alias])
                            parsed = classify_citation(raw_alias)
                            if parsed[0] != alias:
                                parsed = classify_citation(f"{alias}:{raw_alias}")
                            if parsed[0] == alias:
                                lookups.setdefault(parsed, (result, provider))
            record = {"hypothesis_id": original[id_column], "citation": citation,
                      "gene_identity_checks_json": json.dumps(report_checks, ensure_ascii=False, sort_keys=True),
                      "citation_type": kind, "bibliography_state": state,
                      "citation_verified": state in RESOLVED_STATES,
                      "verification_status": "verified" if state in RESOLVED_STATES else "review_required",
                      "claim_support_status": "not_evaluated", "claim_support": "not_checked",
                      "retrieval_scope": "bibliography_only", "source_provider": provider,
                      "metadata_provenance_json": json.dumps((result or {}).get("metadata_provenance", lookup_provenance.get((kind, identifier), [])), sort_keys=True),
                      "resolved_id": str((result or {}).get("pmid") or (result or {}).get("id") or ""),
                      "title": str((result or {}).get("title", "")), "doi": str((result or {}).get("doi", "")),
                      "review_reason": state, "rejection_reason": "", "provider_error": error}
            if kind == "database_accession":
                record["verification_status"] = "not_literature_source"
            group.append(record)
        # Only identifiers explicitly declared as aliases must be one paper.
        # Separate source_ids may cite different papers for the same observation.
        alias_groups = paper_alias_groups if paper_reports else ([source_tokens(list(alias_citations))] if aliases_supplied else [])
        for alias_group in alias_groups:
            resolved = [metadata[r["citation"]] for r in group
                        if r["citation_verified"] and source_tokens([r["citation"]]) & alias_group]
            if len(resolved) < 2:
                continue
            aliases = [{(kind, str(r.get(kind, r.get("id", "") if kind == "pmid" else "")).casefold())
                        for kind in ("pmid", "pmcid", "doi") if r.get(kind) or (kind == "pmid" and r.get("id"))} for r in resolved]
            if any(not aliases[0].intersection(other) for other in aliases[1:]):
                for record in group:
                    if not paper_reports or source_tokens([record["citation"]]) & alias_group:
                        record.update(bibliography_state="multiple_identifiers_mismatch", citation_verified=False,
                                      verification_status="review_required")
        literature_group = [r for r in group if r["citation_type"] != "database_accession"]
        good = bool(literature_group) and all(r["citation_verified"] for r in literature_group)
        rows.extend(group)
        annotated.append({**original, source_column: claim_citations,
                          "bibliography_selection_json": json.dumps(selection, ensure_ascii=False, sort_keys=True),
                          "gene_identity_checks_json": json.dumps(report_checks, ensure_ascii=False, sort_keys=True),
                          "reference_record_source_ids": provider_sources,
                          "candidate_admission_policy": INITIAL_POLICY,
                          "focal_species": focal_species or original.get("focal_species", ""),
                          "bibliography_state": "resolved" if good else ("reference_endpoint_changed" if endpoint_changed else
                              group[0]["bibliography_state"] if group else "source_absent"),
                          "verified_source_ids": ";".join(r["citation"] for r in group if r["citation_verified"]),
                          "source_verification_status": "bibliography_resolved" if good else "not_evaluated",
                          "claim_support_summary": "not_evaluated", "content_check_state": "not_independently_checked"})
    validation = out / "source_validation.tsv"
    pd.DataFrame(rows, columns=["hypothesis_id", "citation", "citation_type", "bibliography_state", "citation_verified",
        "verification_status", "claim_support_status", "claim_support", "retrieval_scope", "source_provider",
        "resolved_id", "title", "doi", "review_reason", "rejection_reason", "provider_error", "gene_identity_checks_json", "metadata_provenance_json"]).to_csv(validation, sep="\t", index=False)
    table = out / "mechanisms_with_verified_sources.tsv"
    pd.DataFrame(annotated, columns=list(dict.fromkeys([*frame.columns, "candidate_admission_policy", "bibliography_state",
        "verified_source_ids", "source_verification_status", "claim_support_summary", "content_check_state", "focal_species", "reference_record_source_ids", "gene_identity_checks_json", "bibliography_selection_json"]))).to_csv(table, sep="\t", index=False)
    selection_path = out / "bibliography_selection.json"
    selection_path.write_text(json.dumps({"panels": selections}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    gene_path = out / "gene_identity_checks.json"
    gene_path.write_text(json.dumps({"checks": gene_rows, "gene_records": gene_checks.records,
        "lookup_count": gene_checks.lookup_count, "membership_effect": "none", "paper_content_checked": False}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    metadata_path = out / "bibliography_records.json"
    metadata_path.write_text(json.dumps({"records": list(metadata.values())}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    sources = out / "source_review_records.json"
    sources.write_text(json.dumps({"artifact_role": "bibliography_only_not_content_support", "records": {},
                                  "bibliography_records": metadata, "membership_effect": "none"}, ensure_ascii=False) + "\n", encoding="utf-8")
    from .paper_evidence_table import write_paper_evidence_table
    paper_evidence = write_paper_evidence_table(annotated, out,
        route="interval_first" if id_column == "claim_id" else "phenotype_first",
        source_artifact=str(runtime_path(mechanisms_arg)), metadata_rows=rows)
    manifest = out / "source_verification_manifest.json"
    manifest.write_text(json.dumps({"schema_version": "2.0", "verification_policy": INITIAL_POLICY,
        "paper_evidence": paper_evidence,
        "provider": "+".join(sorted({p for _, p in lookups.values()})),
        "content_check_state": "not_independently_checked", "lookup_count": len(lookups),
        "paper_selection": {"path": str(selection_path), "sha256": hashlib.sha256(selection_path.read_bytes()).hexdigest(),
                            "max_papers_per_reference_proposal": 5, "selected_paper_count": sum(p["selected_paper_count"] for p in selections),
                            **({"handoff_consumer": "audit_imported_primary_papers_without_reselection"}
                               if any(p.get("policy") == PAPER_HANDOFF_POLICY for p in selections) else {})},
        "gene_identity_checks": {"path": str(gene_path), "sha256": hashlib.sha256(gene_path.read_bytes()).hexdigest(),
                                 "lookup_count": gene_checks.lookup_count},
        "input_sha256": hashlib.sha256(runtime_path(mechanisms_arg).read_bytes()).hexdigest(),
        "metadata_sha256": hashlib.sha256(metadata_path.read_bytes()).hexdigest()}, indent=2) + "\n", encoding="utf-8")
    return {"paper_evidence": paper_evidence, **paper_evidence["outputs"],
            "manifest": str(manifest), "source_validation": str(validation),
            "mechanisms_with_verified_sources": str(table), "source_review_records": str(sources),
            "bibliography_records": str(metadata_path), "gene_identity_checks": str(gene_path), "bibliography_selection": str(selection_path)}
