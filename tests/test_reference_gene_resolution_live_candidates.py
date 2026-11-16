"""Offline contract tests for exact, all-record reference-gene fanout.

These tests deliberately mock every provider call.  They pin the intended
pagination and no-fallback behavior of the new live candidate resolver without
using benchmark truth or real biological records.
"""

from __future__ import annotations

import json
import urllib.parse
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import pandas as pd

from candidate_qtl import llm_pathway_legacy as pathway
from candidate_qtl import reference_gene_resolution_legacy as resolution_legacy
from candidate_qtl.phenotype_investigation import MECHANISM_COLUMNS
from candidate_qtl.reference_gene_resolution import reference_gene_resolve


class _Response:
    def __init__(self, payload: object, *, next_url: str | None = None) -> None:
        self._body = json.dumps(payload).encode("utf-8")
        self.headers = {"Link": f'<{next_url}>; rel="next"'} if next_url else {}

    def read(self) -> bytes:
        return self._body

    def __enter__(self) -> "_Response":
        return self

    def __exit__(self, *args: object) -> None:
        return None


def _result(accession: str, species: str = "Arabidopsis thaliana", *, reviewed: bool = False) -> dict:
    return {
        "primaryAccession": accession,
        "entryType": "UniProtKB reviewed (Swiss-Prot)" if reviewed else "UniProtKB unreviewed (TrEMBL)",
        "organism": {"scientificName": species},
        "genes": [{"geneName": {"value": "FAN1"}, "synonyms": [{"value": "FAN-1"}]}],
    }


def _mechanism(path: Path, *, hypothesis_id: str = "h1") -> None:
    row = {column: "" for column in MECHANISM_COLUMNS}
    row.update({"hypothesis_id": hypothesis_id, "case_id": "live", "reference_species": "Arabidopsis thaliana",
                "reference_gene_name": "FAN1", "gene_family": "synthetic", "source_ids": ""})
    pd.DataFrame([row])[MECHANISM_COLUMNS].to_csv(path, sep="\t", index=False, lineterminator="\n")


class ReferenceGeneLiveCandidatesTests(unittest.TestCase):
    def test_next_page_without_cursor_never_returns_partial_candidates(self) -> None:
        response = _Response({"results": [_result("Q9ABC1")]},
                             next_url="https://rest.uniprot.org/uniprotkb/search?page=2")
        with mock.patch.object(pathway.urllib.request, "urlopen", return_value=response):
            with self.assertRaisesRegex(ValueError, "search is incomplete"):
                pathway._resolve_accession_candidates("FAN1", "Arabidopsis thaliana")

    def setUp(self) -> None:
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.tmp_path = Path(tmp.name)

    def test_candidates_use_link_header_all_pages_include_unreviewed_and_are_sorted_deduped(self) -> None:
        """UniProt pagination is Link-header based; payload.next is not required."""
        calls: list[str] = []
        page2 = "https://rest.uniprot.org/uniprotkb/search?cursor=page-2"

        def fake_open(request: object, timeout: int = 0) -> _Response:
            url = request.full_url  # type: ignore[attr-defined]
            calls.append(url)
            if "cursor=page-2" in url:
                return _Response({"results": [_result("Q9ABC2"), _result("Q9ABC1"), _result("Q9ABC2")]})
            return _Response({"results": [_result("Q9ABC3", reviewed=True), _result("Q9ABC2")]}, next_url=page2)

        with mock.patch.object(pathway.urllib.request, "urlopen", side_effect=fake_open):
            candidates, status, query, hits = pathway._resolve_accession_candidates("FAN1", "Arabidopsis thaliana")
        assert candidates == ["Q9ABC1", "Q9ABC2", "Q9ABC3"]
        assert status == "ambiguous_gene_name"
        assert hits == 3
        assert len(calls) == 2
        assert "reviewed%3Atrue" not in urllib.parse.urlparse(calls[0]).query
        assert query == 'gene:FAN1 AND organism_name:"Arabidopsis thaliana"'


    def test_malformed_provider_payload_fails_closed(self) -> None:
        for payload in ({"results": None}, {"results": [{"organism": []}]}):
            with self.subTest(payload=payload), mock.patch.object(pathway.urllib.request, "urlopen", return_value=_Response(payload)):
                with self.assertRaises(ValueError):
                    pathway._resolve_accession_candidates("FAN1", "Arabidopsis thaliana")


    def test_repeated_link_cursor_does_not_return_partial_candidates(self) -> None:
        loop_url = "https://rest.uniprot.org/uniprotkb/search?cursor=loop"

        def fake_open(request: object, timeout: int = 0) -> _Response:
            return _Response({"results": [_result("Q9ABC1")]}, next_url=loop_url)

        with mock.patch.object(pathway.urllib.request, "urlopen", side_effect=fake_open):
            with self.assertRaises((RuntimeError, ValueError)):
                pathway._resolve_accession_candidates("FAN1", "Arabidopsis thaliana")


    def test_single_exact_unreviewed_candidate_never_calls_old_reviewed_resolver(self) -> None:
        tmp_path = self.tmp_path
        mechanisms = tmp_path / "mechanisms.tsv"
        _mechanism(mechanisms)
        cache = tmp_path / "records.json"
        cache.write_text(json.dumps({"results": [{
            **_result("Q9ABC1", reviewed=False),
            "sequence": {"value": "M" * 24, "length": 24},
            "references": [{"citation": {"citationCrossReferences": [{"database": "PubMed", "id": "12345678"}]}}],
        }]}), encoding="utf-8")
        out = tmp_path / "out"
        with mock.patch("candidate_qtl.reference_gene_resolution._resolve_accession_candidates",
                        return_value=(['Q9ABC1'], "resolved_from_gene_name", "exact-query", 1)), \
             mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_accession_candidates",
                         return_value=(['Q9ABC1'], "resolved_from_gene_name", "exact-query", 1)), \
             mock.patch("candidate_qtl.reference_gene_resolution._resolve_accession",
                         side_effect=AssertionError("old reviewed-only resolver must not be called")), \
             mock.patch.object(pathway.urllib.request, "urlopen",
                               side_effect=AssertionError("unexpected provider call")):
            result = reference_gene_resolve(str(mechanisms), str(out), allow_external_api=True,
                                            record_cache_arg=str(cache))
        resolved = pd.read_csv(out / "mechanism_hypotheses_resolved.tsv", sep="\t", dtype=str).fillna("")
        assert resolved.loc[0, "proposed_stable_id"] == "Q9ABC1"
        assert result["provider_resolved_count"] == 1


    def test_successful_zero_exact_candidates_stays_unresolved_without_old_fallback(self) -> None:
        tmp_path = self.tmp_path
        mechanisms = tmp_path / "mechanisms.tsv"
        _mechanism(mechanisms)
        out = tmp_path / "out"
        with mock.patch("candidate_qtl.reference_gene_resolution._resolve_accession_candidates",
                        return_value=([], "unresolved_gene_name", "exact-query", 0)), \
             mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_accession_candidates",
                         return_value=([], "unresolved_gene_name", "exact-query", 0)), \
             mock.patch("candidate_qtl.reference_gene_resolution._resolve_accession",
                         side_effect=AssertionError("old reviewed-only resolver must not be called")):
            reference_gene_resolve(str(mechanisms), str(out), allow_external_api=True)
        ledger = pd.read_csv(out / "reference_gene_resolution.tsv", sep="\t", dtype=str).fillna("")
        assert ledger.loc[0, "resolution_origin"] == "not_resolved"
        assert ledger.loc[0, "resolved_stable_id"] == ""
        assert ledger.loc[0, "accession_resolution_status"] == "unresolved_gene_name"


    def test_zero_or_late_page_failure_does_not_fallback_to_old_resolver(self) -> None:
        tmp_path = self.tmp_path
        mechanisms = tmp_path / "mechanisms.tsv"
        _mechanism(mechanisms)
        out = tmp_path / "out"
        with mock.patch("candidate_qtl.reference_gene_resolution._resolve_accession_candidates",
                        side_effect=RuntimeError("late provider failure")), \
             mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_accession_candidates",
                         side_effect=RuntimeError("late provider failure")), \
             mock.patch("candidate_qtl.reference_gene_resolution._resolve_accession",
                         side_effect=AssertionError("old resolver fallback")):
            result = reference_gene_resolve(str(mechanisms), str(out), allow_external_api=True)
        ledger = pd.read_csv(out / "reference_gene_resolution.tsv", sep="\t", dtype=str).fillna("")
        assert ledger.loc[0, "resolution_origin"] == "not_resolved"
        assert ledger.loc[0, "evaluation_status"] == "not_evaluated"
        assert result["provider_resolved_count"] == 0

    def test_primary_late_page_failure_keeps_page_one_candidate_and_actual_urls(self) -> None:
        mechanisms = self.tmp_path / "mechanisms.tsv"
        _mechanism(mechanisms)
        out = self.tmp_path / "out"
        page2 = "https://rest.uniprot.org/uniprotkb/search?cursor=second"
        calls: list[str] = []

        def fake_open(request: object, timeout: int = 0) -> _Response:
            url = request.full_url  # type: ignore[attr-defined]
            calls.append(url)
            if len(calls) == 1:
                return _Response({"results": [_result("Q9ABC1")]}, next_url=page2)
            raise RuntimeError("page two unavailable")

        record = {
            **_result("Q9ABC1"),
            "sequence": {"value": "M" * 24, "length": 24},
            "references": [{"citation": {"citationCrossReferences": [{"database": "PubMed", "id": "12345678"}]}}],
        }
        with mock.patch.object(pathway.urllib.request, "urlopen", side_effect=fake_open), \
             mock.patch("candidate_qtl.reference_gene_resolution._fetch_uniprot", return_value=record):
            reference_gene_resolve(str(mechanisms), str(out), allow_external_api=True)

        ledger = pd.read_csv(out / "reference_gene_resolution.tsv", sep="\t", dtype=str).fillna("")
        child = ledger[ledger["derived_accession"].eq("Q9ABC1")].iloc[0]
        provenance = json.loads(child["retrieval_provenance"])
        assert calls == provenance["query_urls"]
        assert len(calls) == 2
        assert provenance["lookup_search_state"] == "partial_not_evaluated"
        assert any("page two unavailable" in failure for failure in provenance["lookup_failures"])
        assert child["resolution_origin"] == "provider_resolved"
        assert child["provider_error"] == ""

    def test_direct_entry_failure_is_record_fetch_failure_with_complete_search(self) -> None:
        mechanisms = self.tmp_path / "mechanisms.tsv"
        _mechanism(mechanisms)
        out = self.tmp_path / "out"
        with mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_accession_candidates",
                        return_value=(["Q9ABC1"], "resolved_from_gene_name", "exact-query", 1)), \
             mock.patch("candidate_qtl.reference_gene_resolution._fetch_uniprot",
                        side_effect=TimeoutError("entry fetch timed out")):
            reference_gene_resolve(str(mechanisms), str(out), allow_external_api=True)

        ledger = pd.read_csv(out / "reference_gene_resolution.tsv", sep="\t", dtype=str).fillna("")
        row = ledger[ledger["derived_accession"].eq("Q9ABC1")].iloc[0]
        provenance = json.loads(row["retrieval_provenance"])
        assert provenance["lookup_search_state"] == "complete"
        assert provenance["record_fetch_failures"] == {"Q9ABC1": "entry fetch timed out"}

    def test_reference_search_next_link_without_cursor_is_reported_partial(self) -> None:
        raw = [{"input_reference_gene_name": "TERMINAL FLOWER 1",
                "reference_gene_name": "TERMINAL FLOWER 1"}]
        response = _Response({"results": []}, next_url="https://rest.uniprot.org/uniprotkb/search?page=2")
        with mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_accession_candidates",
                        return_value=([], "unresolved_gene_name", "symbol query", 0)), \
             mock.patch.object(pathway.urllib.request, "urlopen", return_value=response):
            candidates, _status, _query, _hits, failures, urls = resolution_legacy._resolve_reference_candidates(
                "TERMINAL FLOWER 1", "Arabidopsis thaliana", raw, 1)
        assert candidates == []
        assert urls
        assert any("query=" in url for url in urls)
        assert any("search is incomplete" in failure for failure in failures)
