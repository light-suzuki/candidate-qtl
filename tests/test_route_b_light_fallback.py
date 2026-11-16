"""Exact DOI fallback must preserve the light transport's scientific gates."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import tempfile
import unittest
import urllib.error
from unittest import mock

from scripts import agentic_route_b_light as light
from tests.test_route_b_light import checklist, digest, fixture


DOI = "10.1234/example.article"
TITLE = "A homolog controls leaflet development"
ABSTRACT = "Loss of the homolog produced fused leaflets in Arabidopsis."


def crossref_payload(**changes):
    message = {
        "DOI": DOI, "title": ["A homolog controls <i>leaflet</i> development"],
        "abstract": "<jats:p>" + ABSTRACT + "</jats:p>",
        "container-title": ["Plant Journal"], "published": {"date-parts": [[2024, 1, 2]]},
        "author": [{"given": "Example", "family": "Researcher"}],
    }
    message.update(changes)
    return {"status": "ok", "message": message}


class CrossrefResolverTests(unittest.TestCase):
    def fetch(self, payload=None, primary=None):
        with mock.patch.object(light, "_europe_pmc", return_value=primary) as epmc, \
                mock.patch.object(light, "_http_json", return_value=payload or crossref_payload()) as http:
            record = light._fetch_paper("doi", DOI, timeout=17)
        return record, epmc, http

    def test_exact_doi_fallback_carries_plain_title_abstract_and_provenance(self):
        record, epmc, http = self.fetch()
        epmc.assert_called_once_with("doi", DOI, 17)
        http.assert_called_once_with(light.CROSSREF_WORKS + "10.1234%2Fexample.article", timeout=17)
        self.assertEqual((DOI, TITLE, ABSTRACT), (record["doi"], record["title"], record["abstract"]))
        self.assertEqual("crossref", record["_source_provider"])
        self.assertEqual("not_resolved", record["metadata_provenance"][0]["state"])
        self.assertEqual("retrieved", record["metadata_provenance"][1]["abstract_state"])
        self.assertEqual(64, len(record["metadata_provenance"][1]["payload_sha256"]))
        self.assertEqual(["doi:" + DOI], record["aliases"])

    def test_existing_europe_pmc_record_does_not_call_crossref(self):
        primary = {"doi": DOI, "title": TITLE, "abstract": ABSTRACT}
        record, _, http = self.fetch(primary=primary)
        self.assertEqual(TITLE, record["title"])
        http.assert_not_called()

    def test_fallback_never_borrows_another_doi_even_with_matching_title(self):
        record, _, _ = self.fetch(crossref_payload(DOI="10.1234/other"))
        self.assertIsNone(record)

    def test_missing_abstract_stays_missing_with_an_explicit_trace(self):
        payload = crossref_payload()
        del payload["message"]["abstract"]
        record, _, _ = self.fetch(payload)
        self.assertEqual("", record["abstract"])
        self.assertEqual("not_in_metadata", record["metadata_provenance"][1]["abstract_state"])

    def test_invalid_response_is_a_provider_failure(self):
        for payload in ({"status": "ok", "message": []}, crossref_payload(title="not a list"),
                        crossref_payload(abstract={"invented": "text"})):
            with self.subTest(payload=payload), mock.patch.object(light, "_http_json", return_value=payload):
                with self.assertRaises(RuntimeError):
                    light._crossref_doi(DOI, 30)

    def test_crossref_404_is_a_clean_no_record(self):
        error = urllib.error.HTTPError("https://api.crossref.org", 404, "missing", {}, None)
        with mock.patch.object(light, "_europe_pmc", return_value=None), \
                mock.patch.object(light, "_http_json", side_effect=error):
            self.assertIsNone(light._fetch_paper("doi", DOI))

    def test_primary_outage_can_resolve_but_is_preserved_in_provenance(self):
        with mock.patch.object(light, "_europe_pmc", side_effect=TimeoutError("primary timeout")), \
                mock.patch.object(light, "_http_json", return_value=crossref_payload()):
            record = light._fetch_paper("doi", DOI)
        self.assertEqual("provider_failure", record["metadata_provenance"][0]["state"])
        self.assertIn("primary timeout", record["metadata_provenance"][0]["detail"])

    def test_primary_outage_and_fallback_miss_are_not_paper_absence(self):
        with mock.patch.object(light, "_europe_pmc", side_effect=TimeoutError("primary timeout")), \
                mock.patch.object(light, "_crossref_doi", return_value=None):
            record, state, _ = light._paper_state({}, "DOI:" + DOI, TITLE, False, light._fetch_paper)
        self.assertIsNone(record)
        self.assertEqual("provider_failure", state)

    def test_malformed_primary_response_is_recorded_as_failure_before_fallback(self):
        with mock.patch.object(light, "_http_json", side_effect=[{}, crossref_payload()]):
            record = light._fetch_paper("doi", DOI)
        self.assertEqual("provider_failure", record["metadata_provenance"][0]["state"])
        self.assertIn("invalid result list", record["metadata_provenance"][0]["detail"])

    def test_cache_with_a_conflicting_explicit_doi_is_rejected(self):
        cache = {"doi:" + DOI: {"doi": "10.1234/other", "title": TITLE, "abstract": ABSTRACT}}
        fetch = mock.Mock()
        record, state, detail = light._paper_state(cache, "DOI:" + DOI, TITLE, True, fetch)
        fetch.assert_not_called()
        self.assertIsNone(record)
        self.assertEqual("not_resolved", state)
        self.assertIn("conflicts", detail)


class CrossrefLightTransportTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        bundle = fixture()
        bundle["units"] = bundle["units"][:1]
        bundle["interval_gene_count"] = bundle["semantic_unit_count"] = 1
        self.bundle = self.root / "bundle.json"
        self.bundle.write_text(json.dumps(bundle), encoding="utf-8")
        self.bundle_sha = digest(self.bundle.read_bytes())
        light.light_prepare(str(self.bundle), self.bundle_sha, str(self.root / "light"))
        self.index = self.root / "light/light/index.json"
        self.index_sha = digest(self.index.read_bytes())
        self.row = {
            "interval_gene_id": "g1", "reference_gene": "REF-A", "tier": "T2",
            "function": "leaflet boundary specification", "reason": "Original homolog phenotype bridge.",
            "papers": [{"paper_id": "DOI:" + DOI, "title": TITLE, "reported_gene": "HOM1",
                        "reported_species": "Arabidopsis thaliana", "gene_correspondence": "homolog",
                        "finding": "Original reported fused leaflets.", "checklist": checklist("T2")}],
        }

    def convert(self, payload=None, *, offline=False, cache="", out="converted"):
        results = self.root / "light/light/b01/results.jsonl"
        results.write_text(json.dumps(self.row) + "\n", encoding="utf-8")
        original = results.read_bytes()
        with mock.patch.object(light, "_europe_pmc", return_value=None), \
                mock.patch.object(light, "_http_json", return_value=payload or crossref_payload()) as http:
            light.light_convert(str(self.bundle), self.bundle_sha, str(self.index), self.index_sha,
                                str(self.root / out), offline=offline, cache_arg=cache)
        self.assertEqual(original, results.read_bytes())
        self.assertEqual(self.bundle_sha, digest(self.bundle.read_bytes()))
        return json.loads((self.root / out / "light_check.json").read_text()), http

    def test_convert_keeps_original_checklist_and_exports_reusable_provider_trace(self):
        check, _ = self.convert()
        self.assertEqual(1, check["totals"]["papers_checked"])
        responses = json.loads((self.root / "converted/agent_responses.json").read_text())["responses"]
        claim = json.loads(responses[0]["raw_text"])["claims"][0]
        self.assertEqual("T2", claim["initial_tier"])
        self.assertEqual(self.row["papers"][0]["checklist"], claim["paper_reports_json"][0]["paper_checklist"])
        self.assertEqual("REF-A", claim["reference_gene"])
        discovery = json.loads((self.root / "converted/discoveries.json").read_text())["discovered_sources"][0]
        self.assertEqual(ABSTRACT, discovery["source_text"])
        metadata = json.loads((self.root / "converted/metadata_cache.json").read_text())["doi:" + DOI]
        source = json.loads((self.root / "converted/source_verify_cache.json").read_text())["records"][0]
        self.assertEqual("crossref", metadata["_source_provider"])
        self.assertEqual(metadata["metadata_provenance"], source["metadata_provenance"])
        self.assertEqual(metadata["_source_url"], source["_source_url"])
        check, http = self.convert(offline=True, cache=str(self.root / "converted/metadata_cache.json"), out="reused")
        http.assert_not_called()
        self.assertEqual(1, check["totals"]["papers_checked"])
        reused = json.loads((self.root / "reused/metadata_cache.json").read_text())["doi:" + DOI]
        self.assertEqual(metadata, reused)

    def test_metadata_without_abstract_cannot_form_a_discovery_or_retain_t2(self):
        payload = crossref_payload()
        del payload["message"]["abstract"]
        check, _ = self.convert(payload)
        self.assertEqual("missing_abstract", check["paper_issues"][0]["state"])
        self.assertEqual(0, check["totals"]["claims"])
        self.assertFalse((self.root / "converted/discoveries.json").exists())
        responses = json.loads((self.root / "converted/agent_responses.json").read_text())["responses"]
        self.assertEqual("not_evaluated", responses[0]["status"])

    def test_existing_title_threshold_rejects_fallback_mismatch(self):
        check, _ = self.convert(crossref_payload(title=["Completely unrelated photosynthetic physiology"]))
        self.assertEqual("title_mismatch", check["paper_issues"][0]["state"])
        self.assertEqual(0, check["totals"]["claims"])

    def test_publisher_abstract_enrichment_keeps_its_distinct_source_offline(self):
        with mock.patch.object(light, "_http_json", return_value=crossref_payload(abstract="")):
            record = light._crossref_doi(DOI, 30)
        record["abstract"] = ABSTRACT
        record["abstract_provenance"] = [{
            "provider": "publisher", "url": "https://publisher.example/article/" + DOI,
            "retrieval_transport": "web_search_primary_page", "abstract_sha256": digest(ABSTRACT.encode()),
            "direct_http_status": 403,
        }]
        cache = self.root / "publisher-cache.json"
        cache.write_text(json.dumps({"doi:" + DOI: record}), encoding="utf-8")
        check, http = self.convert(offline=True, cache=str(cache))
        http.assert_not_called()
        self.assertEqual(1, check["totals"]["papers_checked"])
        metadata = json.loads((self.root / "converted/metadata_cache.json").read_text())["doi:" + DOI]
        source = json.loads((self.root / "converted/source_verify_cache.json").read_text())["records"][0]
        for entry in (metadata, source):
            self.assertEqual(record["abstract_provenance"], entry["abstract_provenance"])
            self.assertEqual("crossref", entry["_source_provider"])
            self.assertEqual("not_in_metadata", entry["metadata_provenance"][0]["abstract_state"])

    def test_offline_miss_never_fetches_and_is_not_a_biological_negative(self):
        check, http = self.convert(offline=True)
        http.assert_not_called()
        self.assertEqual("not_attempted", check["paper_issues"][0]["state"])
        self.assertEqual(0, check["totals"]["claims"])

    def test_frozen_endpoint_gate_still_prevents_lookup(self):
        self.row["reference_gene"] = "UNFROZEN"
        check, http = self.convert()
        http.assert_not_called()
        self.assertEqual("reference_not_in_frozen_endpoints", check["gene_issues"][0]["issue"])
        self.assertEqual(0, check["totals"]["claims"])

    def test_fallback_outage_remains_provider_failure_and_is_not_cached(self):
        cache = {}
        with mock.patch.object(light, "_europe_pmc", return_value=None), \
                mock.patch.object(light, "_http_json", side_effect=TimeoutError("Crossref timeout")):
            record, state, detail = light._paper_state(cache, "DOI:" + DOI, TITLE, False, light._fetch_paper)
        self.assertIsNone(record)
        self.assertEqual("provider_failure", state)
        self.assertIn("Crossref timeout", detail)
        self.assertEqual({}, cache)

    def test_missing_original_checklist_answer_is_not_repaired(self):
        original = copy.deepcopy(self.row["papers"][0]["checklist"])
        del self.row["papers"][0]["checklist"]["shared_specific_process"]
        check, _ = self.convert()
        self.assertEqual("paper_report_incomplete", check["paper_issues"][0]["state"])
        self.assertEqual(8, len(self.row["papers"][0]["checklist"]))
        self.assertEqual(9, len(original))


if __name__ == "__main__":
    unittest.main()
