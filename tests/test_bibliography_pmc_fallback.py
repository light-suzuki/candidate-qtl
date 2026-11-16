"""Metadata transport fallback; no paper content or biological judgment."""
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

from candidate_qtl.bibliography import INITIAL_POLICY, fetch_metadata, fetch_pmc_mapping
from candidate_qtl.source_verify import source_verify


class PmcMetadataFallbackTests(unittest.TestCase):
    def run_case(self, root, *, cache=None, external=True, mapping=None, fetched=None):
        table = root / "input.tsv"
        pd.DataFrame([{"hypothesis_id": str(i), "source_ids": "PMC1234567",
                       "source_identifier_aliases": "PMID:1234567", "source_title": "Synthetic paper"}
                      for i in range(2)]).to_csv(table, sep="\t", index=False)
        cache_path = root / "cache.json"
        cache_path.write_text(json.dumps(cache or []), encoding="utf-8")
        payload = {"status": "ok", "records": mapping if mapping is not None else [
            {"pmcid": "PMC1234567", "pmid": 1234567, "doi": "10.1234/paper"}]}
        def metadata_lookup(kind, identifier):
            value = (fetched or {}).get((kind, identifier))
            if isinstance(value, Exception):
                raise value
            return value
        with patch("candidate_qtl.bibliography.urllib.request.urlopen", return_value=io.BytesIO(json.dumps(payload).encode())) as converter, \
                patch("candidate_qtl.bibliography.fetch_metadata", side_effect=metadata_lookup) as lite, \
                patch("candidate_qtl.source_verify._claim_text", side_effect=AssertionError("content was read")):
            result = source_verify(str(table), str(root / "out"), cache_arg=str(cache_path),
                                   allow_external_api=external, verification_policy=INITIAL_POLICY)
        validation = pd.read_csv(result["source_validation"], sep="\t", dtype=str).fillna("")
        annotated = pd.read_csv(result["mechanisms_with_verified_sources"], sep="\t", dtype=str).fillna("")
        return result, validation, annotated, converter, lite

    def test_missing_pmc_uses_official_mapping_and_cached_title_once(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result, validation, annotated, converter, lite = self.run_case(root, cache=[
                {"pmid": "1234567", "doi": "10.1234/paper", "title": "Synthetic paper", "abstractText": "never retain"}])
            self.assertEqual({"resolved"}, set(annotated.bibliography_state))
            self.assertEqual(1, converter.call_count)
            self.assertEqual([("pmcid", "PMC1234567")], [call.args for call in lite.call_args_list])
            pmc = validation[validation.citation == "PMC1234567"]
            self.assertTrue(pmc.citation_verified.eq("True").all())
            self.assertTrue(pmc.source_provider.str.contains("ncbi_pmc_id_converter").all())
            metadata = Path(result["bibliography_records"]).read_text(encoding="utf-8")
            self.assertIn("metadata_provenance", metadata)
            self.assertNotIn("abstractText", metadata)
            self.assertEqual({"not_evaluated"}, set(validation.claim_support_status))
            replay = root / "replay"
            replay.mkdir()
            _, replay_validation, replay_annotated, converter, lite = self.run_case(
                replay, cache=json.loads(metadata)["records"], external=False)
            self.assertEqual({"resolved"}, set(replay_annotated.bibliography_state))
            self.assertTrue(replay_validation.metadata_provenance_json.str.contains("ncbi_pmc_id_converter").all())
            converter.assert_not_called()
            lite.assert_not_called()

    def test_live_title_and_doi_only_mapping_reuse_one_paper(self):
        for mapping, fetched, calls in [
            (None, {("pmid", "1234567"): {"pmid": "1234567", "doi": "10.1234/paper", "title": "Synthetic paper"}},
             [("pmcid", "PMC1234567"), ("pmid", "1234567")]),
            ([{"pmcid": "PMC1234567", "doi": "10.1234/paper"}],
             {("doi", "10.1234/paper"): {"pmid": "1234567", "doi": "10.1234/paper", "title": "Synthetic paper"}},
             [("pmcid", "PMC1234567"), ("doi", "10.1234/paper")]),
        ]:
            with self.subTest(mapping=mapping), tempfile.TemporaryDirectory() as tmp:
                _, _, annotated, converter, lite = self.run_case(Path(tmp), mapping=mapping, fetched=fetched)
                self.assertEqual({"resolved"}, set(annotated.bibliography_state))
                self.assertEqual(calls, [call.args for call in lite.call_args_list])
                self.assertEqual(1, converter.call_count)

    def test_cached_doi_precedes_unfetched_mapped_pmid(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, _, annotated, converter, lite = self.run_case(Path(tmp), cache=[
                {"doi": "10.1234/paper", "title": "Synthetic paper"}])
            self.assertEqual({"resolved"}, set(annotated.bibliography_state))
            self.assertEqual([("pmcid", "PMC1234567")], [call.args for call in lite.call_args_list])
            self.assertEqual(1, converter.call_count)

    def test_title_mismatch_and_nonmatching_primary_result_do_not_bypass_gate(self):
        for fetched, conversion_count in [
            ({("pmid", "1234567"): {"pmid": "1234567", "title": "Different title"}}, 1),
            ({("pmcid", "PMC1234567"): {"pmcid": "PMC9999999", "title": "Synthetic paper"}}, 0),
        ]:
            with self.subTest(fetched=fetched), tempfile.TemporaryDirectory() as tmp:
                _, _, annotated, converter, _ = self.run_case(Path(tmp), fetched=fetched)
                self.assertFalse(annotated.bibliography_state.eq("resolved").any())
                self.assertEqual(conversion_count, converter.call_count)

    def test_converter_http_failure_is_explicit_and_never_content_fallback(self):
        with patch("candidate_qtl.bibliography.urllib.request.urlopen", side_effect=OSError("HTTP 503")):
            with self.assertRaisesRegex(RuntimeError, "provider unavailable.*503"):
                fetch_pmc_mapping("PMC1234567")

    def test_primary_provider_failure_does_not_start_converter(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, validation, annotated, converter, _ = self.run_case(Path(tmp), fetched={
                ("pmcid", "PMC1234567"): RuntimeError("bibliography provider unavailable: HTTP 429")})
            converter.assert_not_called()
            self.assertFalse(annotated.bibliography_state.eq("resolved").any())
            self.assertEqual({"provider_failure"}, set(validation[validation.citation == "PMC1234567"].bibliography_state))

    def test_mapped_provider_failure_survives_later_no_hit_but_allows_later_success(self):
        for doi_result, expected in [(None, "provider_failure"),
                                     ({"doi": "10.1234/paper", "title": "Synthetic paper"}, "id_title_matched")]:
            with self.subTest(expected=expected), tempfile.TemporaryDirectory() as tmp:
                _, validation, annotated, converter, lite = self.run_case(Path(tmp), fetched={
                    ("pmid", "1234567"): RuntimeError("bibliography provider unavailable: HTTP 503"),
                    ("doi", "10.1234/paper"): doi_result})
                pmc = validation[validation.citation == "PMC1234567"]
                self.assertEqual({expected}, set(pmc.bibliography_state))
                self.assertEqual(1, converter.call_count)
                self.assertEqual(3, lite.call_count)
                if doi_result:
                    self.assertEqual({"resolved"}, set(annotated.bibliography_state))
                    self.assertTrue(pmc.metadata_provenance_json.str.contains("HTTP 503").all())
                else:
                    self.assertTrue(pmc.provider_error.str.contains("HTTP 503").all())

    def test_title_lookup_rejects_multiple_matching_records(self):
        payload = {"resultList": {"result": [{"pmid": "1234567", "title": "Synthetic paper"}] * 2}}
        with patch("candidate_qtl.bibliography.urllib.request.urlopen", return_value=io.BytesIO(json.dumps(payload).encode())):
            with self.assertRaisesRegex(RuntimeError, "ambiguous metadata identifier"):
                fetch_metadata("pmid", "1234567")

    def test_converter_request_is_scoped_to_the_missing_pmc_identifier(self):
        payload = {"status": "ok", "records": [{"pmcid": "PMC1234567", "pmid": 1234567}]}
        with patch("candidate_qtl.bibliography.urllib.request.urlopen", return_value=io.BytesIO(json.dumps(payload).encode())) as fetch:
            mapping, provenance = fetch_pmc_mapping("PMC1234567")
        self.assertEqual({"pmcid": "PMC1234567", "pmid": "1234567"}, mapping)
        request = fetch.call_args.args[0]
        self.assertTrue(request.full_url.startswith("https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?"))
        self.assertIn("ids=PMC1234567", request.full_url)
        self.assertIn("idtype=pmcid", request.full_url)
        self.assertEqual(64, len(provenance["response_sha256"]))

    def test_no_external_fallback_for_offline_cache_or_existing_pmc(self):
        for external, cache in [(False, []), (True, [{"pmcid": "PMC1234567", "pmid": "1234567", "title": "Synthetic paper"}])]:
            with self.subTest(external=external), tempfile.TemporaryDirectory() as tmp:
                _, _, _, converter, lite = self.run_case(Path(tmp), cache=cache, external=external)
                converter.assert_not_called()
                lite.assert_not_called()

    def test_bad_mapping_or_conflicting_metadata_never_resolves_all_aliases(self):
        for mapping, fetched in [
            ([{"pmcid": "PMC9999999", "pmid": 1234567}], {}),
            ([{"pmcid": "PMC1234567", "pmid": 1234567}] * 2, {}),
            ([], {}),
            ([{"pmcid": "PMC1234567", "pmid": 1234567}], {("pmid", "1234567"): {"pmid": "7654321", "title": "Synthetic paper"}}),
            ([{"pmcid": "PMC1234567", "pmid": 1234567}], {("pmid", "1234567"): {"pmid": "1234567", "pmcid": "PMC9999999", "title": "Synthetic paper"}}),
            ([{"pmcid": "PMC1234567", "pmid": 1234567}], {("pmid", "1234567"): {"pmid": "1234567", "title": ""}}),
        ]:
            with self.subTest(mapping=mapping, fetched=fetched), tempfile.TemporaryDirectory() as tmp:
                _, validation, annotated, converter, _ = self.run_case(Path(tmp), mapping=mapping, fetched=fetched)
                self.assertFalse(annotated.bibliography_state.eq("resolved").any())
                self.assertFalse(validation[validation.citation == "PMC1234567"].citation_verified.eq("True").any())
                self.assertEqual(1, converter.call_count)
