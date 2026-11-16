"""Offline contract tests for explicit Route A provider-record acquisition."""

from __future__ import annotations

import csv
import hashlib
import json
import tempfile
import unittest
import urllib.parse
from pathlib import Path
from unittest import mock

import pandas as pd

from candidate_qtl.reference_gene_resolution import reference_gene_resolve
from candidate_qtl.reference_record_cache import build_reference_record_cache, validate_reference_record_cache


def _mechanisms(path: Path) -> None:
    columns = ["case_id", "hypothesis_id", "reference_gene_name", "reference_species",
               "source_ids", "proposed_stable_id"]
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow({"case_id": "neutral_case", "hypothesis_id": "h001",
                         "reference_gene_name": "GENE1", "reference_species": "Species neutralis",
                         "source_ids": "PMID:12345678", "proposed_stable_id": ""})


def _record(accession: str, pmid: str) -> dict:
    return {"primaryAccession": accession, "organism": {"scientificName": "Species neutralis"},
            "genes": [{"geneName": {"value": "GENE1"}}],
            "sequence": {"value": "MKLLTAVWQ" * 12, "length": 108},
            "references": [{"citation": {"citationCrossReferences": [
                {"database": "PubMed", "id": pmid}]}}]}


def _search_response(records: list[dict]):
    response = mock.MagicMock()
    response.__enter__.return_value = response
    response.read.return_value = json.dumps({"results": records}).encode()
    response.headers = {}
    return response


class ReferenceRecordCacheTests(unittest.TestCase):
    def test_cited_record_fallback_keeps_exact_aliases_and_all_accessions(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _mechanisms(root / "in.tsv")
            records = [_record("Q00001", "12345678"), _record("Q00002", "12345678")]
            records[1]["genes"] = [{"geneName": {"value": "FORMAL1"},
                                   "synonyms": [{"value": "GENE1"}]}]
            by_id = {record["primaryAccession"]: record for record in records}
            with mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_accession_candidates",
                            return_value=([], "unresolved_gene_name", "gene query", 0)), \
                 mock.patch("candidate_qtl.reference_gene_resolution_legacy.urllib.request.urlopen",
                            return_value=_search_response(records)) as search, \
                 mock.patch("candidate_qtl.reference_record_cache._fetch_uniprot",
                            side_effect=lambda accession, timeout: by_id[accession]):
                result = build_reference_record_cache(str(root / "in.tsv"), str(root / "cache"),
                                                      allow_external_api=True)
            _, ledger = validate_reference_record_cache(result["manifest"])
            self.assertEqual(2, result["record_count"])
            self.assertEqual(["Q00001", "Q00002"], ledger[0]["cached_accessions"])
            url = search.call_args.args[0].full_url
            self.assertEqual(['lit_pubmed:12345678 AND organism_name:"Species neutralis"'],
                             urllib.parse.parse_qs(urllib.parse.urlparse(url).query)["query"])
            self.assertIn(url, ledger[0]["query_urls"])
            resolved = reference_gene_resolve(str(root / "in.tsv"), str(root / "resolved"),
                                              record_cache_arg=result["records_dir"])
            children = pd.read_csv(resolved["mechanisms_resolved"], sep="\t").fillna("")
            self.assertEqual(2, len(children))
            self.assertEqual({"Q00001", "Q00002"}, set(children["proposed_stable_id"]))

    def test_citation_fallback_rejects_wrong_name_species_and_citation(self):
        from candidate_qtl.reference_gene_resolution_legacy import _resolve_reference_candidates
        records = [_record("Q00001", "12345678"), _record("Q00002", "12345678"),
                   _record("Q00003", "999")]
        records[0]["genes"][0]["geneName"]["value"] = "OTHER"
        records[1]["organism"]["scientificName"] = "Species secunda"
        with mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_accession_candidates",
                        return_value=([], "unresolved_gene_name", "query", 0)), \
             mock.patch("candidate_qtl.reference_gene_resolution_legacy.urllib.request.urlopen",
                        return_value=_search_response(records)):
            candidates = _resolve_reference_candidates("GENE1", "Species neutralis",
                                                       [{"source_ids": "PMID:12345678"}], 10)
        self.assertEqual([], candidates[0])
        self.assertEqual([], candidates[4])

    def test_citation_fallback_does_not_mask_failure_or_repeat_success(self):
        from candidate_qtl.reference_gene_resolution_legacy import _resolve_reference_candidates
        for result in [(["Q00001"], "resolved", "query", 1), TimeoutError("provider failed")]:
            with self.subTest(result=result), \
                 mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_accession_candidates") as lookup, \
                 mock.patch("candidate_qtl.reference_gene_resolution_legacy.urllib.request.urlopen") as search:
                if isinstance(result, Exception):
                    lookup.side_effect = result
                else:
                    lookup.return_value = result
                candidates = _resolve_reference_candidates("GENE1", "Species neutralis",
                                                           [{"source_ids": "PMID:12345678"}], 10)
                search.assert_not_called()
                self.assertEqual(bool(isinstance(result, Exception)), bool(candidates[4]))

    def test_failed_citation_lookup_stays_provider_unavailable(self):
        from candidate_qtl.reference_gene_resolution_legacy import _resolve_reference_candidates
        with mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_accession_candidates",
                        return_value=([], "unresolved_gene_name", "query", 0)), \
             mock.patch("candidate_qtl.reference_gene_resolution_legacy.urllib.request.urlopen",
                        side_effect=TimeoutError("citation lookup failed")):
            candidates = _resolve_reference_candidates("GENE1", "Species neutralis",
                                                       [{"source_ids": "PMID:12345678"}], 10)
        self.assertEqual([], candidates[0])
        self.assertEqual("provider_unavailable", candidates[1])
        self.assertIn("citation lookup failed", candidates[4][0])

    def test_late_healthy_records_survive_an_earlier_accession_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _mechanisms(root / "in.tsv")
            def fetch(accession, timeout):
                if accession == "Q00002":
                    raise TimeoutError("local accession failure")
                return _record(accession, "123")
            with mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_accession_candidates",
                            return_value=(["Q00001", "Q00002", "Q00003"], "ambiguous_gene_name", "query", 3)), \
                 mock.patch("candidate_qtl.reference_record_cache._fetch_uniprot", side_effect=fetch):
                result = build_reference_record_cache(str(root / "in.tsv"), str(root / "cache"), allow_external_api=True)
            self.assertEqual({"Q00001", "Q00003"}, set(result["records_sha256"]))
            _, ledger = validate_reference_record_cache(result["manifest"])
            self.assertEqual({"Q00002"}, set(ledger[0]["failed_accessions"]))
            self.assertFalse(result["cache_complete"])

    def test_ledger_case_and_hash_negative_guards(self):
        import hashlib
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _mechanisms(root / "in.tsv")
            with mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_accession_candidates",
                            return_value=([], "unresolved_gene_name", "query", 0)), \
                 mock.patch("candidate_qtl.reference_gene_resolution_legacy.urllib.request.urlopen",
                            return_value=_search_response([])):
                result = build_reference_record_cache(str(root / "in.tsv"), str(root / "cache"), allow_external_api=True)
            ledger_path = Path(result["query_ledger"])
            ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
            ledger[0]["hypothesis_keys"] = ["Neutral_case/h001"]
            ledger_path.write_text(json.dumps(ledger), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "ledger hash"):
                validate_reference_record_cache(result["manifest"])
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            manifest["query_ledger_sha256"] = hashlib.sha256(ledger_path.read_bytes()).hexdigest()
            Path(result["manifest"]).write_text(json.dumps(manifest), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "case/hypothesis"):
                validate_reference_record_cache(result["manifest"])
            with (root / "in.tsv").open("a", encoding="utf-8") as stream:
                stream.write("other\th002\tGENE2\tSpecies neutralis\t\t\n")
            with self.assertRaisesRegex(ValueError, "mechanisms hash"):
                validate_reference_record_cache(result["manifest"])

    def test_retry_only_failed_query_preserves_frozen_checkpoint_and_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanisms.tsv"
            _mechanisms(mechanisms)
            with mechanisms.open("a", encoding="utf-8") as stream:
                stream.write("neutral_case\th002\tGENE2\tSpecies neutralis\t\t\n")
            with mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_accession_candidates",
                            side_effect=[(["Q00001"], "resolved", "first", 1), TimeoutError("failed second")]), \
                 mock.patch("candidate_qtl.reference_record_cache._fetch_uniprot", return_value=_record("Q00001", "123")):
                first = build_reference_record_cache(str(mechanisms), str(root / "first"), allow_external_api=True)
            frozen = {str(p): p.read_bytes() for p in (root / "first").rglob("*") if p.is_file()}
            record = _record("Q00002", "456")
            record["genes"][0]["geneName"]["value"] = "GENE2"
            with mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_accession_candidates",
                            return_value=(["Q00002"], "resolved", "second", 1)) as query, \
                 mock.patch("candidate_qtl.reference_record_cache._fetch_uniprot", return_value=record):
                second = build_reference_record_cache(str(mechanisms), str(root / "second"), allow_external_api=True,
                                                       resume_manifest_arg=first["manifest"])
            query.assert_called_once_with("GENE2", "Species neutralis", 60)
            self.assertTrue(second["cache_complete"])
            self.assertEqual(2, second["record_count"])
            self.assertEqual(1, second["logical_query_attempt_count"])
            for path, contents in frozen.items():
                self.assertEqual(contents, Path(path).read_bytes())
            validate_reference_record_cache(second["manifest"])
            (root / "second" / "records" / "Q00001.json").write_text("{}", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "roster/hash"):
                validate_reference_record_cache(second["manifest"])

    def test_composite_cache_queries_each_atomic_pair_and_preserves_raw_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanisms.tsv"
            _mechanisms(mechanisms)
            frame = pd.read_csv(mechanisms, sep="\t", dtype=str).fillna("")
            frame.loc[0, "reference_gene_name"] = "GENE1 (one); GENE2/ALT2 (two)"
            frame.loc[0, "reference_species"] = "Species neutralis; Species secunda"
            frame.to_csv(mechanisms, sep="\t", index=False)
            with mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_accession_candidates",
                            return_value=([], "unresolved_gene_name", "exact query", 0)) as query, \
                 mock.patch("candidate_qtl.reference_gene_resolution_legacy.urllib.request.urlopen",
                            return_value=_search_response([])):
                result = build_reference_record_cache(str(mechanisms), str(root / "cache"), allow_external_api=True)
            self.assertEqual({("GENE1", "Species neutralis", 60), ("GENE2", "Species secunda", 60),
                              ("ALT2", "Species secunda", 60)}, {call.args for call in query.call_args_list})
            ledger = json.loads(Path(result["query_ledger"]).read_text(encoding="utf-8"))
            self.assertTrue(all(item["raw_inputs"][0]["reference_gene_name"] == frame.loc[0, "reference_gene_name"]
                                for item in ledger))

    def test_formal_name_and_typed_handles_are_lookup_only_until_record_identity_matches(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanisms.tsv"
            _mechanisms(mechanisms)
            frame = pd.read_csv(mechanisms, sep="\t", dtype=str).fillna("")
            frame.loc[0, "reference_gene_name"] = "TERMINAL FLOWER 1"
            frame.loc[0, "reference_species"] = "Species neutralis"
            frame.loc[0, "proposed_stable_id"] = "UniProt:P93003"
            frame.loc[0, "source_ids"] = "TAIR:AT5G03840"
            frame.to_csv(mechanisms, sep="\t", index=False)

            def search_response(request, timeout):
                actual = urllib.parse.parse_qs(urllib.parse.urlparse(request.full_url).query)["query"][0]
                accession = "Q00001" if "TERMINAL FLOWER 1" in actual else "P93003"
                response = mock.MagicMock()
                response.read.return_value = json.dumps({"results": [{"primaryAccession": accession}]}).encode()
                response.headers.get.return_value = ""
                response.__enter__.return_value = response
                return response

            good = _record("Q00001", "11111111")
            good["genes"][0]["geneName"]["value"] = "TFL1"
            good["genes"][0]["synonyms"] = [{"value": "TERMINAL FLOWER 1"}]
            mismatch = _record("P93003", "22222222")
            mismatch["genes"][0]["geneName"]["value"] = "OTHER"
            with mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_accession_candidates",
                            side_effect=TimeoutError("symbol query timeout")), \
                 mock.patch("candidate_qtl.reference_gene_resolution_legacy.urllib.request.urlopen",
                            side_effect=search_response) as search, \
                 mock.patch("candidate_qtl.reference_record_cache._fetch_uniprot",
                            side_effect=lambda accession, timeout: {"Q00001": good, "P93003": mismatch}[accession]):
                result = build_reference_record_cache(str(mechanisms), str(root / "cache"), allow_external_api=True)

            ledger = json.loads(Path(result["query_ledger"]).read_text(encoding="utf-8"))[0]
            self.assertEqual({"Q00001"}, set(ledger["cached_accessions"]), ledger)
            self.assertEqual([("uniprot", "P93003"), ("tair", "AT5G03840")],
                             [tuple(item) for item in ledger["lookup_handles"]])
            self.assertIn('gene:"TERMINAL FLOWER 1"', ledger["query"])
            self.assertIn("xref:tair-AT5G03840", ledger["query"])
            self.assertEqual(2, len(ledger["query_urls"]))
            self.assertTrue(all("query=" in url for url in ledger["query_urls"]))
            self.assertEqual("", ledger["search_url"])
            self.assertIn("symbol query timeout", ledger["reason"])
            self.assertIn("P93003", ledger["failed_accessions"])
            self.assertFalse(result["cache_complete"])
            self.assertEqual(2, search.call_count)
            resumed = reference_gene_resolve(str(mechanisms), str(root / "resolved"),
                                             record_cache_arg=result["records_dir"])
            resolved_rows = pd.read_csv(resumed["ledger"], sep="\t", dtype=str).fillna("")
            healthy = resolved_rows[resolved_rows["resolved_stable_id"].eq("Q00001")].iloc[0]
            provenance = json.loads(healthy["retrieval_provenance"])
            self.assertEqual("partial_not_evaluated", provenance["lookup_search_state"])
            self.assertEqual(2, len(provenance["query_urls"]))
            self.assertTrue(any("symbol query timeout" in value for value in provenance["lookup_failures"]))
            self.assertNotIn("record_fetch_failures", provenance)
            self.assertEqual({"P93003"}, set(provenance["record_identity_failures"]))
            self.assertEqual("provider_resolved", healthy["resolution_origin"])
            self.assertEqual("", healthy["provider_error"])

    def test_ambiguous_pairing_is_explicit_not_evaluated_without_queries(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanisms.tsv"
            _mechanisms(mechanisms)
            frame = pd.read_csv(mechanisms, sep="\t", dtype=str).fillna("")
            frame.loc[0, "reference_species"] = "Species alpha; Species beta"
            frame.to_csv(mechanisms, sep="\t", index=False)
            with mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_accession_candidates") as query:
                result = build_reference_record_cache(str(mechanisms), str(root / "cache"), allow_external_api=True)
            query.assert_not_called()
            self.assertFalse(result["cache_complete"])
            self.assertEqual("ambiguous_name_species_pairing", result["unresolved_inputs"][0]["reason"])

    def test_one_exact_record_is_frozen_then_resolved_offline(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanisms.tsv"
            _mechanisms(mechanisms)
            with mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_accession_candidates",
                            return_value=(["Q00001"], "resolved_from_gene_name", "query", 1)), \
                 mock.patch("candidate_qtl.reference_record_cache._fetch_uniprot",
                            return_value=_record("Q00001", "99999999")):
                result = build_reference_record_cache(str(mechanisms), str(root / "frozen"),
                                                      allow_external_api=True)
            self.assertTrue(result["cache_usable_for_prepare"])
            self.assertEqual(1, result["record_count"])
            self.assertEqual(64, len(result["records_sha256"]["Q00001"]))
            ledger = json.loads((root / "frozen" / "query_ledger.json").read_text(encoding="utf-8"))
            self.assertTrue(ledger[0]["retrieved_at"])
            self.assertIn("uniprotkb/search", ledger[0]["search_url"])
            self.assertIn("Q00001.json", ledger[0]["record_urls"]["Q00001"])
            reference_gene_resolve(str(mechanisms), str(root / "resolved"),
                                   record_cache_arg=result["records_dir"])
            rows = pd.read_csv(root / "resolved" / "mechanism_hypotheses_resolved.tsv",
                               sep="\t", dtype=str).fillna("")
            self.assertEqual(["Q00001"], rows["proposed_stable_id"].tolist())
            self.assertIn("PMID:99999999", rows.loc[0, "source_ids"])
            self.assertNotIn("PMID:12345678", rows.loc[0, "source_ids"])

    def test_parenthetical_short_symbol_survives_cache_validation_and_offline_resolver(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanisms.tsv"
            _mechanisms(mechanisms)
            frame = pd.read_csv(mechanisms, sep="\t", dtype=str).fillna("")
            frame.loc[0, "reference_gene_name"] = "LONG  NODE (LN1)"
            frame.to_csv(mechanisms, sep="\t", index=False)
            record = _record("Q00001", "99999999")
            record["genes"][0]["geneName"]["value"] = "LN1"
            with mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_accession_candidates",
                            return_value=(["Q00001"], "resolved_from_gene_name", "query", 1)), \
                 mock.patch("candidate_qtl.reference_record_cache._fetch_uniprot", return_value=record):
                cache = build_reference_record_cache(str(mechanisms), str(root / "cache"),
                                                     allow_external_api=True)
            self.assertEqual(["Q00001"], json.loads(Path(cache["query_ledger"]).read_text(encoding="utf-8"))[0]["cached_accessions"])
            validate_reference_record_cache(cache["manifest"])
            resolved = reference_gene_resolve(str(mechanisms), str(root / "resolved"),
                                              record_cache_arg=cache["records_dir"])
            rows = pd.read_csv(resolved["ledger"], sep="\t", dtype=str).fillna("")
            self.assertEqual(["Q00001"], rows["resolved_stable_id"].tolist())
            self.assertEqual(["provider_resolved"], rows["resolution_origin"].tolist())

    def test_multiple_exact_records_remain_distinct_children(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanisms.tsv"
            _mechanisms(mechanisms)
            records = {"Q00001": _record("Q00001", "99999991"),
                       "Q00002": _record("Q00002", "99999992")}
            with mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_accession_candidates",
                            return_value=(["Q00002", "Q00001"], "ambiguous_gene_name", "query", 2)), \
                 mock.patch("candidate_qtl.reference_record_cache._fetch_uniprot",
                            side_effect=lambda accession, timeout: records[accession]):
                result = build_reference_record_cache(str(mechanisms), str(root / "frozen"),
                                                      allow_external_api=True)
            self.assertTrue(result["cache_usable_for_prepare"])
            reference_gene_resolve(str(mechanisms), str(root / "resolved"),
                                   record_cache_arg=result["records_dir"])
            rows = pd.read_csv(root / "resolved" / "mechanism_hypotheses_resolved.tsv",
                               sep="\t", dtype=str).fillna("")
            self.assertEqual({"Q00001", "Q00002"}, set(rows["proposed_stable_id"]))
            self.assertEqual({"h001__ref_q00001", "h001__ref_q00002"}, set(rows["hypothesis_id"]))
            self.assertEqual({"PMID:99999991", "PMID:99999992"}, set(rows["source_ids"]))
            self.assertNotIn("PMID:12345678", " ".join(rows["source_ids"]))

    def test_provider_failure_retains_not_evaluated_and_healthy_partial_records(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanisms.tsv"
            _mechanisms(mechanisms)
            with mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_accession_candidates",
                            return_value=(["Q00001", "Q00002"], "ambiguous_gene_name", "query", 2)), \
                 mock.patch("candidate_qtl.reference_record_cache._fetch_uniprot",
                            side_effect=[_record("Q00001", "99999991"), RuntimeError("provider down")]):
                result = build_reference_record_cache(str(mechanisms), str(root / "frozen"),
                                                      allow_external_api=True)
            self.assertTrue(result["cache_usable_for_prepare"])
            self.assertFalse(result["cache_complete"])
            self.assertTrue((root / "frozen" / "records" / "Q00001.json").exists())
            ledger = json.loads((root / "frozen" / "query_ledger.json").read_text(encoding="utf-8"))
            self.assertEqual("not_evaluated", ledger[0]["status"])
            # Simulate a receipt from before query URLs/failure details were added.
            for key in ("query_urls", "lookup_failures", "record_fetch_failures"):
                ledger[0].pop(key, None)
            ledger_path = root / "frozen" / "query_ledger.json"
            ledger_path.write_text(json.dumps(ledger), encoding="utf-8")
            manifest_path = root / "frozen" / "manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["query_ledger_sha256"] = hashlib.sha256(ledger_path.read_bytes()).hexdigest()
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            reference_gene_resolve(str(mechanisms), str(root / "resolved"), record_cache_arg=result["records_dir"])
            resolution = pd.read_csv(root / "resolved" / "reference_gene_resolution.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"evaluated", "not_evaluated"}, set(resolution["evaluation_status"]))
            self.assertEqual({"Q00001", ""}, set(resolution["resolved_stable_id"]))
            self.assertIn("provider down", resolution.loc[resolution["evaluation_status"] == "not_evaluated", "provider_error"].iloc[0])
            self.assertTrue(all(value == "query" for value in resolution["accession_resolution_query"]))
            healthy = resolution[resolution["resolved_stable_id"].eq("Q00001")].iloc[0]
            provenance = json.loads(healthy["retrieval_provenance"])
            self.assertEqual("provider_resolved", healthy["resolution_origin"])
            self.assertEqual("complete", provenance["lookup_search_state"])
            self.assertEqual({"Q00002"}, set(provenance["record_fetch_failures"]))

    def test_elapsed_budget_before_query_stays_not_evaluated_after_offline_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanisms.tsv"
            _mechanisms(mechanisms)
            with mock.patch("candidate_qtl.reference_record_cache.time.monotonic", side_effect=[10.0, 20.0]), \
                 mock.patch("candidate_qtl.reference_record_cache._resolve_reference_candidates") as query:
                cache = build_reference_record_cache(str(mechanisms), str(root / "cache"),
                                                     allow_external_api=True, max_elapsed_seconds=1)
            query.assert_not_called()
            receipt = json.loads(Path(cache["query_ledger"]).read_text(encoding="utf-8"))[0]
            self.assertEqual("not_evaluated", receipt["status"])
            self.assertEqual("elapsed_budget_exhausted", receipt["reason"])
            self.assertEqual(0, receipt["logical_query_attempts"])
            self.assertEqual([], receipt["lookup_failures"])
            resolved = reference_gene_resolve(str(mechanisms), str(root / "resolved"),
                                             record_cache_arg=cache["records_dir"])
            row = pd.read_csv(resolved["ledger"], sep="\t", dtype=str).fillna("").iloc[0]
            provenance = json.loads(row["retrieval_provenance"])
            self.assertEqual("not_evaluated", row["evaluation_status"])
            self.assertEqual("provider_unavailable", row["accession_resolution_status"])
            self.assertIn("elapsed_budget_exhausted", row["provider_error"])
            self.assertEqual("not_evaluated", provenance["lookup_search_state"])
            self.assertEqual([], provenance["query_urls"])

    def test_record_budget_preserves_partial_lookup_receipt_and_retry(self):
        for failures in ([], ["TimeoutError: alternate name query unavailable"]):
            with self.subTest(failures=failures), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                mechanisms = root / "mechanisms.tsv"
                _mechanisms(mechanisms)
                urls = ["https://rest.uniprot.org/uniprotkb/search?query=GENE1"]
                lookup = (["Q00001", "Q00002"], "ambiguous_gene_name", "fixture query", 2, failures, urls)
                with mock.patch("candidate_qtl.reference_record_cache._resolve_reference_candidates", return_value=lookup), \
                     mock.patch("candidate_qtl.reference_record_cache._fetch_uniprot") as fetch:
                    cache = build_reference_record_cache(str(mechanisms), str(root / "cache"),
                                                         allow_external_api=True, max_records=1)
                fetch.assert_not_called()
                _, ledger = validate_reference_record_cache(cache["manifest"])
                self.assertFalse(cache["cache_complete"])
                self.assertEqual("not_evaluated", ledger[0]["status"])
                self.assertEqual("record_budget_exceeded", ledger[0]["reason"])
                self.assertEqual(failures, ledger[0]["lookup_failures"])
                self.assertEqual(urls, ledger[0]["query_urls"])
                self.assertEqual(["Q00001", "Q00002"], ledger[0]["accessions"])
                frozen = {str(p): p.read_bytes() for p in (root / "cache").rglob("*") if p.is_file()}
                resolved = reference_gene_resolve(str(mechanisms), str(root / "resolved"), record_cache_arg=cache["records_dir"])
                row = pd.read_csv(resolved["ledger"], sep="\t", dtype=str).fillna("").iloc[0]
                provenance = json.loads(row["retrieval_provenance"])
                self.assertEqual("not_evaluated", row["evaluation_status"])
                self.assertEqual("provider_unavailable", row["accession_resolution_status"])
                self.assertEqual("partial_not_evaluated" if failures else "complete", provenance["lookup_search_state"])
                self.assertEqual(failures, provenance["lookup_failures"])
                self.assertEqual({"cache": "record_budget_exceeded"}, provenance["record_fetch_failures"])
                with mock.patch("candidate_qtl.reference_record_cache._resolve_reference_candidates",
                                return_value=(*lookup[:4], [], urls)) as query, \
                     mock.patch("candidate_qtl.reference_record_cache._fetch_uniprot",
                                side_effect=lambda accession, timeout: _record(accession, "123")):
                    retry = build_reference_record_cache(str(mechanisms), str(root / "retry"), allow_external_api=True,
                                                         max_records=2, resume_manifest_arg=cache["manifest"])
                query.assert_called_once()
                self.assertTrue(retry["cache_complete"])
                self.assertEqual(2, retry["record_count"])
                validate_reference_record_cache(retry["manifest"])
                for path, contents in frozen.items():
                    self.assertEqual(contents, Path(path).read_bytes())

    def test_no_hit_is_complete_but_identity_remains_unresolved(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanisms.tsv"
            _mechanisms(mechanisms)
            with mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_accession_candidates",
                            return_value=([], "unresolved_gene_name", "query", 0)), \
                 mock.patch("candidate_qtl.reference_gene_resolution_legacy.urllib.request.urlopen",
                            return_value=_search_response([])), \
                 mock.patch("candidate_qtl.reference_record_cache._fetch_uniprot") as fetch:
                result = build_reference_record_cache(str(mechanisms), str(root / "frozen"),
                                                      allow_external_api=True)
                fetch.assert_not_called()
            self.assertTrue(result["cache_usable_for_prepare"])
            self.assertEqual(0, result["record_count"])
            self.assertEqual(1, result["status_counts"]["no_provider_match"])
            reference_gene_resolve(str(mechanisms), str(root / "resolved"),
                                   record_cache_arg=result["records_dir"])
            rows = pd.read_csv(root / "resolved" / "mechanism_hypotheses_resolved.tsv",
                               sep="\t", dtype=str).fillna("")
            self.assertEqual("", rows.loc[0, "proposed_stable_id"])
            self.assertEqual("", rows.loc[0, "source_ids"])

    def test_alias_mismatch_and_timeout_do_not_yield_complete_cache(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanisms.tsv"
            _mechanisms(mechanisms)
            for label, response in (("mismatch", {**_record("Q00001", "99999991"),
                                                  "genes": [{"geneName": {"value": "OTHER"}}]}),
                                    ("timeout", TimeoutError("provider timeout"))):
                with self.subTest(label=label), \
                     mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_accession_candidates",
                                return_value=(["Q00001"], "resolved_from_gene_name", "query", 1)), \
                     mock.patch("candidate_qtl.reference_record_cache._fetch_uniprot",
                                side_effect=response if isinstance(response, Exception) else None,
                                return_value=response if isinstance(response, dict) else None):
                    result = build_reference_record_cache(str(mechanisms), str(root / label),
                                                          allow_external_api=True)
                    self.assertFalse(result["cache_complete"])
                    self.assertEqual(1, result["status_counts"]["not_evaluated"])

    def test_network_opt_in_and_query_budget_are_checked_before_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mechanisms = root / "mechanisms.tsv"
            _mechanisms(mechanisms)
            with mechanisms.open("a", encoding="utf-8", newline="") as stream:
                stream.write("neutral_case\th002\tGENE2\tSpecies neutralis\t\t\n")
            with mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_accession_candidates") as query:
                with self.assertRaisesRegex(ValueError, "--allow-external-api"):
                    build_reference_record_cache(str(mechanisms), str(root / "no_api"))
                with self.assertRaisesRegex(ValueError, "exceed max_queries"):
                    build_reference_record_cache(str(mechanisms), str(root / "no_budget"),
                                                 allow_external_api=True, max_queries=1)
                query.assert_not_called()
            self.assertFalse((root / "no_api").exists())
            self.assertFalse((root / "no_budget").exists())


if __name__ == "__main__":
    unittest.main()
