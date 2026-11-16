from __future__ import annotations

import json
import tempfile
import unittest
import urllib.error
from pathlib import Path
from unittest import mock

import pandas as pd

from candidate_qtl.reference_gene_resolution import reference_gene_resolve


def _input(path: Path, name: str, species: str, stable_id: str = "bad__model_id"):
    pd.DataFrame([{"case_id": "neutral", "hypothesis_id": "h001", "reference_gene_name": name,
                   "reference_species": species, "source_ids": "PMID:999", "proposed_stable_id": stable_id}]).to_csv(path, sep="\t", index=False)


def _record(accession: str, name: str, species: str):
    return {"primaryAccession": accession, "organism": {"scientificName": species},
            "genes": [{"geneName": {"value": name}}], "sequence": {"value": "M" * 30},
            "references": [{"citation": {"citationCrossReferences": [{"database": "PubMed", "id": "123"}]}}]}


class AtomicReferenceResolutionTests(unittest.TestCase):
    def test_partial_symbol_search_keeps_exact_handle_record_and_records_query_receipt(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _input(root / "in.tsv", "TERMINAL FLOWER 1", "Species alpha", "UniProt:P93003")
            mechanisms = pd.read_csv(root / "in.tsv", sep="\t", dtype=str).fillna("")
            mechanisms.loc[0, "source_ids"] = "PMID:999; TAIR:AT5G03840"
            mechanisms.to_csv(root / "in.tsv", sep="\t", index=False)
            good = _record("Q00001", "TFL1", "Species alpha")
            good["genes"][0]["synonyms"] = [{"value": "TERMINAL FLOWER 1"}]
            mismatch = _record("P93003", "OTHER", "Species alpha")
            response = mock.MagicMock()
            response.read.return_value = json.dumps({"results": [{"primaryAccession": "Q00001"}]}).encode()
            response.headers.get.return_value = ""
            response.__enter__.return_value = response
            def search(request, timeout):
                if "xref%3Atair" in request.full_url or "xref:tair" in request.full_url:
                    raise TimeoutError("TAIR query timeout")
                return response
            with mock.patch("candidate_qtl.reference_gene_resolution._resolve_accession_candidates",
                            side_effect=TimeoutError("symbol query timeout")), \
                 mock.patch("candidate_qtl.reference_gene_resolution_legacy.urllib.request.urlopen",
                            side_effect=search) as search_calls, \
                 mock.patch("candidate_qtl.reference_gene_resolution._fetch_uniprot",
                            side_effect=lambda accession, timeout: {"Q00001": good, "P93003": mismatch}[accession]):
                result = reference_gene_resolve(str(root / "in.tsv"), str(root / "out"), allow_external_api=True)

            ledger = pd.read_csv(result["ledger"], sep="\t", dtype=str).fillna("").set_index("resolved_stable_id")
            self.assertIn("Q00001", ledger.index)
            record = ledger.loc["Q00001"]
            self.assertEqual("provider_resolved", record["resolution_origin"])
            self.assertEqual("", record["provider_error"])
            provenance = json.loads(record["retrieval_provenance"])
            self.assertEqual("partial_not_evaluated", provenance["lookup_search_state"])
            self.assertIn("symbol query timeout", " ".join(provenance["lookup_failures"]))
            self.assertEqual(2, len(provenance["query_urls"]))
            self.assertEqual("resolved_exact_name_species", provenance["record_identity_state"])
            self.assertEqual(2, search_calls.call_count)
            self.assertIn("provider_record_identity_mismatch", set(ledger["accession_resolution_status"]))

    def test_partial_search_failure_remains_on_unresolved_parent_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _input(root / "in.tsv", "TERMINAL FLOWER 1", "Species alpha", "")
            with mock.patch("candidate_qtl.reference_gene_resolution._resolve_accession_candidates",
                            side_effect=TimeoutError("symbol query timeout")), \
                 mock.patch("candidate_qtl.reference_gene_resolution_legacy.urllib.request.urlopen",
                            side_effect=TimeoutError("formal query timeout")):
                result = reference_gene_resolve(str(root / "in.tsv"), str(root / "out"), allow_external_api=True)
            ledger = pd.read_csv(result["ledger"], sep="\t", dtype=str).fillna("").iloc[0]
            provenance = json.loads(ledger["retrieval_provenance"])
            self.assertEqual("partial_not_evaluated", provenance["lookup_search_state"])
            self.assertEqual(2, len(provenance["lookup_failures"]))
            self.assertEqual(1, len(provenance["query_urls"]))
            self.assertEqual("unresolved", provenance["record_identity_state"])
            self.assertEqual("not_evaluated", ledger["evaluation_status"])

    def test_model_row_without_entry_attempt_is_not_rewritten_by_other_case_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows = [
                {"case_id": "model-case", "hypothesis_id": "h001", "reference_gene_name": "",
                 "reference_species": "Species alpha", "source_ids": "PMID:222", "proposed_stable_id": "Q00001"},
                {"case_id": "lookup-case", "hypothesis_id": "h001", "reference_gene_name": "NODE1",
                 "reference_species": "Species alpha", "source_ids": "", "proposed_stable_id": "Q00001"},
            ]
            pd.DataFrame(rows).to_csv(root / "in.tsv", sep="\t", index=False)
            failure = urllib.error.HTTPError("https://rest.uniprot.org/uniprotkb/Q00001.json", 400, "bad request", None, None)
            with mock.patch("candidate_qtl.reference_gene_resolution._resolve_accession_candidates",
                            return_value=(["Q00001"], "resolved_from_gene_name", "query", 1)), \
                 mock.patch("candidate_qtl.uniprot_error_policy.urllib.request.urlopen", side_effect=failure) as provider:
                result = reference_gene_resolve(str(root / "in.tsv"), str(root / "out"), allow_external_api=True)
            self.assertEqual(1, provider.call_count)
            resolved = pd.read_csv(result["mechanisms_resolved"], sep="\t", dtype=str).fillna("").set_index("case_id")
            ledger = pd.read_csv(result["ledger"], sep="\t", dtype=str).fillna("").set_index("case_id")
            self.assertEqual("Q00001", resolved.loc["model-case", "proposed_stable_id"])
            self.assertEqual("PMID:222", resolved.loc["model-case", "source_ids"])
            self.assertEqual("model_supplied", ledger.loc["model-case", "resolution_origin"])
            self.assertEqual("", resolved.loc["lookup-case", "proposed_stable_id"])
            self.assertEqual("invalid_http_400", ledger.loc["lookup-case", "accession_resolution_status"])

    def test_same_accession_alias_attempt_failure_does_not_rewrite_healthy_attempt(self):
        failures = [
            (urllib.error.HTTPError("https://rest.uniprot.org/uniprotkb/Q00001.json", 400, "bad request", None, None), "invalid_http_400", "evaluated"),
            (urllib.error.HTTPError("https://rest.uniprot.org/uniprotkb/Q00001.json", 404, "missing", None, None), "not_found_http_404", "evaluated"),
            (TimeoutError("second alias timed out"), "provider_unavailable", "not_evaluated"),
        ]
        for failure, expected_status, expected_evaluation in failures:
            with self.subTest(status=expected_status), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                _input(root / "in.tsv", "NODE1/NODE2", "Species alpha", "")
                record = _record("Q00001", "NODE1", "Species alpha")
                record["genes"][0]["synonyms"] = [{"value": "NODE2"}]
                response = mock.MagicMock()
                response.read.return_value = json.dumps(record).encode()
                response.__enter__.return_value = response

                def query(name, species, timeout):
                    return (["Q00001"], "resolved_from_gene_name", f'gene:{name} AND organism_name:"{species}"', 1)

                with mock.patch("candidate_qtl.reference_gene_resolution._resolve_accession_candidates", side_effect=query), \
                     mock.patch("candidate_qtl.uniprot_error_policy.urllib.request.urlopen", side_effect=[response, failure]) as provider:
                    result = reference_gene_resolve(str(root / "in.tsv"), str(root / "out"), allow_external_api=True)
                self.assertEqual(2, provider.call_count)
                ledger = pd.read_csv(result["ledger"], sep="\t", dtype=str).fillna("").set_index("reference_gene_name")
                resolved = pd.read_csv(result["mechanisms_resolved"], sep="\t", dtype=str).fillna("").set_index("reference_gene_name")
                self.assertEqual("Q00001", resolved.loc["NODE1", "proposed_stable_id"])
                self.assertEqual("PMID:123", resolved.loc["NODE1", "source_ids"])
                self.assertEqual("True", ledger.loc["NODE1", "record_found"])
                self.assertEqual("provider_resolved", ledger.loc["NODE1", "resolution_origin"])
                self.assertEqual("", resolved.loc["NODE2", "proposed_stable_id"])
                self.assertEqual("", resolved.loc["NODE2", "source_ids"])
                self.assertEqual(expected_status, ledger.loc["NODE2", "accession_resolution_status"])
                self.assertEqual(expected_evaluation, ledger.loc["NODE2", "evaluation_status"])
                for name, event_state, attempt_index in (("NODE1", "retrieved", 0), ("NODE2", expected_status, 1)):
                    self.assertIn(f"gene:{name}", ledger.loc[name, "accession_resolution_query"])
                    provenance = json.loads(ledger.loc[name, "retrieval_provenance"])
                    self.assertEqual(name, provenance["atomic_reference_gene_name"])
                    attempt = provenance["uniprot_entry_attempts"][0]
                    self.assertEqual(attempt_index, attempt["event_index"])
                    self.assertEqual(event_state, attempt["state"])
                    self.assertEqual("Q00001", attempt["accession"])
                    self.assertEqual(provenance["uniprot_entry_attempts"], json.loads(resolved.loc[name, "retrieval_provenance"])["uniprot_entry_attempts"])
                manifest = json.loads(Path(result["manifest"]).read_text())
                self.assertEqual({"retrieved": 1, expected_status: 1}, manifest["uniprot_entry_outcome_counts"])

    def test_http_failure_does_not_clear_other_case_with_same_hypothesis_id(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows = [{"case_id": case, "hypothesis_id": "h001", "reference_gene_name": name,
                     "reference_species": "Species alpha", "source_ids": "", "proposed_stable_id": accession}
                    for case, name, accession in (("case-good", "NODE1", "Q00001"), ("case-bad", "NODE2", "Q00002"))]
            pd.DataFrame(rows).to_csv(root / "in.tsv", sep="\t", index=False)
            (root / "cache.json").write_text(json.dumps({"results": [_record("Q00001", "NODE1", "Species alpha")]}), encoding="utf-8")
            error = urllib.error.HTTPError("https://rest.uniprot.org/uniprotkb/Q00002.json", 400, "bad request", None, None)
            def query(name, species, timeout):
                return (["Q00001" if name == "NODE1" else "Q00002"], "resolved_from_gene_name", "exact", 1)

            with mock.patch("candidate_qtl.reference_gene_resolution._resolve_accession_candidates", side_effect=query), \
                 mock.patch("candidate_qtl.uniprot_error_policy.urllib.request.urlopen", side_effect=error):
                result = reference_gene_resolve(str(root / "in.tsv"), str(root / "out"),
                                               record_cache_arg=str(root / "cache.json"), allow_external_api=True)
            resolved = pd.read_csv(result["mechanisms_resolved"], sep="\t", dtype=str).fillna("").set_index("case_id")
            self.assertEqual("Q00001", resolved.loc["case-good", "proposed_stable_id"])
            self.assertEqual("", resolved.loc["case-bad", "proposed_stable_id"])
            ledger = pd.read_csv(result["ledger"], sep="\t", dtype=str).fillna("").set_index("case_id")
            self.assertEqual("True", ledger.loc["case-good", "record_found"])
            self.assertEqual("invalid_http_400", ledger.loc["case-bad", "accession_resolution_status"])

    def test_same_hypothesis_id_in_different_cases_keeps_distinct_provenance(self):
        for names in (("NODE1/NODE2", "NODE1/NODE2 (second)"), ("NODE1", "NODE2")):
            with self.subTest(names=names), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                rows = [{"case_id": case, "hypothesis_id": "h001", "reference_gene_name": name,
                         "reference_species": "Species alpha", "source_ids": "PMID:999", "proposed_stable_id": ""}
                        for case, name in zip(("case-one", "case-two"), names)]
                pd.DataFrame(rows).to_csv(root / "in.tsv", sep="\t", index=False)
                records = [_record("Q00001", "NODE1", "Species alpha"), _record("Q00002", "NODE2", "Species alpha")]
                (root / "cache.json").write_text(json.dumps({"results": records}), encoding="utf-8")
                result = reference_gene_resolve(str(root / "in.tsv"), str(root / "out"), record_cache_arg=str(root / "cache.json"))
                ledger = pd.read_csv(result["ledger"], sep="\t", dtype=str).fillna("")
                self.assertEqual(len(names[0].split("/")) * 2, len(ledger))
                for case, name in zip(("case-one", "case-two"), names):
                    part = ledger[ledger["case_id"] == case]
                    self.assertEqual({name}, set(part["input_reference_gene_name"]))
                    for provenance in part["retrieval_provenance"]:
                        self.assertEqual(name, json.loads(provenance)["input_reference_gene_name"])

    def test_http_outcome_preserves_other_case_review_reason(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows = [{"case_id": case, "hypothesis_id": "h001", "reference_gene_name": name,
                     "reference_species": "Species alpha", "source_ids": "", "proposed_stable_id": accession}
                    for case, name, accession in (("case-no-publication", "NODE1", "Q00001"), ("case-bad", "NODE2", "Q00002"))]
            pd.DataFrame(rows).to_csv(root / "in.tsv", sep="\t", index=False)
            record = _record("Q00001", "NODE1", "Species alpha")
            record["references"] = []
            (root / "cache.json").write_text(json.dumps({"results": [record]}), encoding="utf-8")
            error = urllib.error.HTTPError("https://rest.uniprot.org/uniprotkb/Q00002.json", 400, "bad request", None, None)
            def query(name, species, timeout):
                return (["Q00001" if name == "NODE1" else "Q00002"], "resolved_from_gene_name", "exact", 1)
            with mock.patch("candidate_qtl.reference_gene_resolution._resolve_accession_candidates", side_effect=query), \
                 mock.patch("candidate_qtl.uniprot_error_policy.urllib.request.urlopen", side_effect=error):
                result = reference_gene_resolve(str(root / "in.tsv"), str(root / "out"),
                                               record_cache_arg=str(root / "cache.json"), allow_external_api=True)
            reviews = pd.read_csv(result["review_required"], sep="\t", dtype=str).fillna("")
            self.assertIn("reference_gene_resolution_required:reference_gene_function_is_not_source_backed", set(reviews["reason"]))
            typed = reviews[reviews["reason"].eq("reference_gene_resolution_required:invalid_submitted_accession")]
            self.assertEqual(1, len(typed))
            self.assertIn("case_id=case-bad", typed.iloc[0]["required_next_action"])

    def test_composite_fans_out_all_exact_records_with_atomic_endpoints_and_original_parent(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw_name = "NODE1 (one); NODE2/ALT2 (two)"
            raw_species = "Species alpha; Species beta"
            _input(root / "in.tsv", raw_name, raw_species)
            records = [_record("Q00001", "NODE1", "Species alpha"), _record("Q00002", "NODE1", "Species alpha"),
                       _record("Q00003", "NODE2", "Species beta"), _record("Q00004", "ALT2", "Species beta")]
            (root / "cache.json").write_text(json.dumps({"results": records}), encoding="utf-8")
            result = reference_gene_resolve(str(root / "in.tsv"), str(root / "out"), record_cache_arg=str(root / "cache.json"))
            resolved = pd.read_csv(result["mechanisms_resolved"], sep="\t", dtype=str).fillna("")
            ledger = pd.read_csv(result["ledger"], sep="\t", dtype=str).fillna("")
            self.assertEqual({"Q00001", "Q00002", "Q00003", "Q00004"}, set(resolved["proposed_stable_id"]))
            self.assertEqual({"NODE1", "NODE2", "ALT2"}, set(resolved["reference_gene_name"]))
            self.assertEqual({"h001"}, set(resolved["parent_hypothesis_id"]))
            self.assertEqual(4, resolved["hypothesis_id"].nunique())
            self.assertEqual({"bad__model_id"}, set(ledger["input_proposed_stable_id"]))
            self.assertEqual({"malformed_proposed_stable_id_ignored"}, set(ledger["malformed_proposed_stable_id_status"]))
            self.assertEqual({raw_name}, set(ledger["input_reference_gene_name"]))
            for value in ledger["retrieval_provenance"]:
                provenance = json.loads(value)
                self.assertEqual(raw_species, provenance["input_reference_species"])
                self.assertIn(provenance["atomic_reference_gene_name"], {"NODE1", "NODE2", "ALT2"})

    def test_live_query_and_record_identity_accept_only_explicit_atomic_short_name(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _input(root / "in.tsv", "LONG NODE (LN1)", "Species alpha")
            response = mock.MagicMock()
            response.read.return_value = json.dumps({"results": [{"primaryAccession": "Q00001"}]}).encode()
            response.headers.get.return_value = ""
            response.__enter__.return_value = response
            alias_record = _record("Q00001", "LN1", "Species alpha")
            with mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_accession_candidates",
                            return_value=([], "unresolved_gene_name", "primary", 0)), \
                 mock.patch("candidate_qtl.reference_gene_resolution_legacy.urllib.request.urlopen",
                            return_value=response) as search, \
                 mock.patch("candidate_qtl.reference_gene_resolution._fetch_uniprot",
                            return_value=alias_record):
                result = reference_gene_resolve(str(root / "in.tsv"), str(root / "out"),
                                                allow_external_api=True)
            self.assertTrue(any('gene%3A%22LN1%22' in call.args[0].full_url for call in search.call_args_list))
            rows = pd.read_csv(result["ledger"], sep="\t", dtype=str).fillna("")
            self.assertEqual(["Q00001"], rows["resolved_stable_id"].tolist())
            self.assertEqual(["provider_resolved"], rows["resolution_origin"].tolist())

    def test_ambiguous_pairing_does_not_query_or_authorize_supplied_id(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _input(root / "in.tsv", "NODE1/NODE2", "Species alpha; Species beta", "Q00001")
            with mock.patch("candidate_qtl.reference_gene_resolution._resolve_accession_candidates") as query:
                result = reference_gene_resolve(str(root / "in.tsv"), str(root / "out"), allow_external_api=True)
            query.assert_not_called()
            ledger = pd.read_csv(result["ledger"], sep="\t", dtype=str).fillna("")
            self.assertEqual("not_evaluated", ledger.loc[0, "evaluation_status"])
            self.assertEqual("ambiguous_name_species_pairing", ledger.loc[0, "provider_error"])
            self.assertEqual("", ledger.loc[0, "resolved_stable_id"])

    def test_live_search_hit_does_not_bypass_fetched_record_identity(self):
        for record in [_record("Q00001", "NODE1 unrelated", "Species alpha"),
                       _record("Q00002", "NODE1", "Species alpha"),
                       _record("Q00001", "NODE1", "Species beta")]:
            with self.subTest(record=record), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                _input(root / "in.tsv", "NODE1", "Species alpha", "Q00009")
                with mock.patch("candidate_qtl.reference_gene_resolution._resolve_accession_candidates",
                                return_value=(["Q00001"], "resolved_from_gene_name", "exact", 1)), \
                     mock.patch("candidate_qtl.reference_gene_resolution._fetch_uniprot", return_value=record):
                    result = reference_gene_resolve(str(root / "in.tsv"), str(root / "out"), allow_external_api=True)
                ledger = pd.read_csv(result["ledger"], sep="\t", dtype=str).fillna("")
                resolved = pd.read_csv(result["mechanisms_resolved"], sep="\t", dtype=str).fillna("")
                self.assertEqual("provider_record_identity_mismatch", ledger.loc[0, "review_reason"])
                self.assertEqual("", resolved.loc[0, "proposed_stable_id"])
                self.assertEqual("", resolved.loc[0, "source_ids"])

    def test_duplicate_case_variant_input_ids_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _input(root / "in.tsv", "NODE1", "Species alpha")
            frame = pd.read_csv(root / "in.tsv", sep="\t", dtype=str).fillna("")
            duplicate = frame.copy()
            duplicate.loc[0, "case_id"] = "Neutral"
            pd.concat([frame, duplicate]).to_csv(root / "in.tsv", sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "case-exact"):
                reference_gene_resolve(str(root / "in.tsv"), str(root / "out"))
