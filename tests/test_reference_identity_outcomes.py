"""Cache-to-resolver identity/transport separation (workflow2 issue #15)."""
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl.reference_gene_resolution import reference_gene_resolve
from candidate_qtl.reference_gene_resolution_legacy import _cache_lookup_receipt
from candidate_qtl.reference_record_cache import build_reference_record_cache, validate_reference_record_cache
from tests.test_reference_record_cache import _mechanisms, _record, _search_response


class ReferenceIdentityOutcomesTests(unittest.TestCase):
    def _roundtrip(self, root, outcomes, name="GENE1"):
        mechanisms = root / "in.tsv"
        _mechanisms(mechanisms)
        frame = pd.read_csv(mechanisms, sep="\t", dtype=str).fillna("")
        frame.loc[0, "reference_gene_name"] = name
        frame.to_csv(mechanisms, sep="\t", index=False)

        def fetch(accession, timeout):
            value = outcomes[accession]
            if isinstance(value, Exception):
                raise value
            return value

        with patch("candidate_qtl.reference_record_cache._resolve_reference_candidates",
                   return_value=(sorted(outcomes), "ambiguous_gene_name", "synthetic query", len(outcomes), [], [])), \
             patch("candidate_qtl.reference_record_cache._fetch_uniprot", side_effect=fetch):
            cache = build_reference_record_cache(str(mechanisms), str(root / "cache"), allow_external_api=True)
        _, receipts = validate_reference_record_cache(cache["manifest"])
        resolved = reference_gene_resolve(str(mechanisms), str(root / "resolved"), record_cache_arg=cache["records_dir"])
        ledger = pd.read_csv(resolved["ledger"], sep="\t", dtype=str, keep_default_na=False)
        return cache, receipts[0], ledger

    def test_identity_only_is_evaluated_unresolved_without_provider_error(self):
        for mismatch in ("name", "species", "primary_accession"):
            with self.subTest(mismatch=mismatch), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                record = _record("Q00001", "12345678")
                if mismatch == "name":
                    record["genes"] = [{"geneName": {"value": "OTHER"}}]
                elif mismatch == "species":
                    record["organism"]["scientificName"] = "Species beta"
                else:
                    record["primaryAccession"] = "Q00002"
                cache, receipt, ledger = self._roundtrip(root, {"Q00001": record})
                self.assertFalse(cache["cache_complete"])
                self.assertEqual(0, cache["record_count"])
                self.assertEqual({}, receipt["record_fetch_failures"])
                self.assertEqual({"Q00001"}, set(receipt["record_identity_failures"]))
                self.assertEqual(["provider_record_identity_mismatch"], list(ledger["accession_resolution_status"]))
                self.assertEqual(["evaluated"], list(ledger["evaluation_status"]))
                self.assertEqual([""], list(ledger["provider_error"]))
                self.assertEqual(["name_species_unresolved"], list(ledger["reference_identity_status"]))
                self.assertEqual([""], list(ledger["resolved_stable_id"]))
                provenance = json.loads(ledger.loc[0, "retrieval_provenance"])
                self.assertEqual("complete", provenance["lookup_search_state"])
                self.assertEqual("mismatch", provenance["record_identity_state"])
                self.assertNotIn("record_fetch_failures", provenance)
                self.assertEqual({"Q00001"}, set(provenance["record_identity_failures"]))

    def test_transport_and_parse_failures_remain_not_evaluated(self):
        for failure in (TimeoutError("timeout"), RuntimeError("HTTP 503"), ValueError("invalid JSON")):
            with self.subTest(failure=failure), tempfile.TemporaryDirectory() as temp:
                cache, receipt, ledger = self._roundtrip(Path(temp), {"Q00001": failure})
                self.assertEqual(0, cache["record_count"])
                self.assertEqual({}, receipt["record_identity_failures"])
                self.assertEqual({"Q00001"}, set(receipt["record_fetch_failures"]))
                self.assertEqual(["provider_unavailable"], list(ledger["accession_resolution_status"]))
                self.assertEqual(["not_evaluated"], list(ledger["evaluation_status"]))
                self.assertTrue(ledger.loc[0, "provider_error"])

    def test_healthy_identity_transport_alternatives_keep_separate_provenance(self):
        for with_transport in (False, True):
            with self.subTest(with_transport=with_transport), tempfile.TemporaryDirectory() as temp:
                outcomes = {"Q00001": _record("Q00001", "12345678"),
                            "Q00002": {**_record("Q00002", "12345678"), "genes": []}}
                if with_transport:
                    outcomes["Q00003"] = TimeoutError("transport timeout")
                cache, receipt, ledger = self._roundtrip(Path(temp), outcomes)
                self.assertEqual(1, cache["record_count"])
                self.assertEqual(["Q00001"], receipt["cached_accessions"])
                self.assertEqual({"Q00002"}, set(receipt["record_identity_failures"]))
                self.assertEqual({"Q00003"} if with_transport else set(), set(receipt["record_fetch_failures"]))
                self.assertEqual({"h001", "h001__pending"}, set(ledger["hypothesis_id"]))
                healthy = ledger[ledger["resolved_stable_id"].eq("Q00001")].iloc[0]
                self.assertEqual("evaluated", healthy["evaluation_status"])
                self.assertEqual("", healthy["provider_error"])
                pending = ledger[ledger["hypothesis_id"].eq("h001__pending")].iloc[0]
                self.assertEqual("provider_unavailable" if with_transport else "provider_record_identity_mismatch",
                                 pending["accession_resolution_status"])
                self.assertEqual("not_evaluated" if with_transport else "evaluated", pending["evaluation_status"])
                self.assertNotIn("Q00002", pending["provider_error"])
                for row in (healthy, pending):
                    provenance = json.loads(row["retrieval_provenance"])
                    self.assertEqual({"Q00002"}, set(provenance["record_identity_failures"]))
                    self.assertEqual({"Q00003"} if with_transport else set(), set(provenance.get("record_fetch_failures", {})))

    def test_orf_only_records_roundtrip_retains_all_exact_alternatives(self):
        with tempfile.TemporaryDirectory() as temp:
            outcomes = {key: {**_record(key, "12345678"), "genes": [{"orfNames": [{"value": "LOC_001"}]}]}
                        for key in ("Q00001", "Q00002")}
            cache, _, ledger = self._roundtrip(Path(temp), outcomes, name="LOC_001")
            self.assertTrue(cache["cache_complete"])
            self.assertEqual({"Q00001", "Q00002"}, set(ledger["resolved_stable_id"]))
            self.assertTrue(ledger["reference_identity_status"].eq("exact_name_species").all())

    def test_orf_literal_query_runs_when_symbol_lookup_cannot_query_the_name(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            mechanisms = root / "in.tsv"
            _mechanisms(mechanisms)
            frame = pd.read_csv(mechanisms, sep="\t", dtype=str).fillna("")
            frame.loc[0, "reference_gene_name"] = "LOC_001"
            frame.to_csv(mechanisms, sep="\t", index=False)
            record = {**_record("Q00001", "12345678"), "genes": [{"orfNames": [{"value": "LOC_001"}]}]}
            with patch("candidate_qtl.reference_gene_resolution_legacy.urllib.request.urlopen",
                       return_value=_search_response([record])) as search, \
                 patch("candidate_qtl.reference_record_cache._fetch_uniprot", return_value=record):
                cache = build_reference_record_cache(str(mechanisms), str(root / "cache"), allow_external_api=True)
            self.assertEqual(1, search.call_count)
            self.assertEqual(1, cache["record_count"])
            receipt = json.loads(Path(cache["query_ledger"]).read_text(encoding="utf-8"))[0]
            self.assertEqual('gene:"LOC_001" AND organism_name:"Species neutralis"', receipt["query"])

    def test_legacy_receipt_identity_messages_are_separated_and_unknown_errors_retained(self):
        for failure_key in ("failed_accessions", "record_fetch_failures"):
            receipt = {"status": "not_evaluated", "query": "synthetic query", failure_key: {
                "Q00001": "ValueError: provider record name/species mismatch: Q00001",
                "Q00002": "TimeoutError: timeout",
                "Q00003": "ValueError: unknown failure",
                "Q00004": "ValueError: provider record identity mismatch: Q00004",
            }}
            state, _, lookups, fetches, identities = _cache_lookup_receipt(receipt)
            self.assertEqual("complete", state)
            self.assertEqual([], lookups)
            self.assertEqual({"Q00001"}, set(identities))
            self.assertEqual({"Q00002", "Q00003", "Q00004"}, set(fetches))
        ambiguous = {"status": "not_evaluated", "query": "synthetic query", "failed_accessions": {
            "Q00001": "ValueError: provider record identity mismatch: Q00001"}}
        _, _, _, fetches, identities = _cache_lookup_receipt(ambiguous)
        self.assertEqual({"Q00001"}, set(fetches))
        self.assertEqual({}, identities)
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            _, receipt, _ = self._roundtrip(root, {"Q00001": {**_record("Q00001", "12345678"), "genes": []}})
            receipt.pop("record_identity_failures")
            receipt.pop("record_fetch_failures")
            receipt["failed_accessions"] = {"Q00001": "ValueError: provider record name/species mismatch: Q00001"}
            path = root / "cache" / "query_ledger.json"
            path.write_text(json.dumps([receipt]), encoding="utf-8")
            manifest_path = root / "cache" / "manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["query_ledger_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            result = reference_gene_resolve(str(root / "in.tsv"), str(root / "legacy-resolved"),
                                            record_cache_arg=str(root / "cache" / "records"))
            ledger = pd.read_csv(result["ledger"], sep="\t", dtype=str, keep_default_na=False)
            self.assertEqual("provider_record_identity_mismatch", ledger.loc[0, "accession_resolution_status"])
            self.assertEqual("", ledger.loc[0, "provider_error"])
            receipt["failed_accessions"] = {"Q00001": "ValueError: provider record identity mismatch: Q00001"}
            path.write_text(json.dumps([receipt]), encoding="utf-8")
            manifest["query_ledger_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            result = reference_gene_resolve(str(root / "in.tsv"), str(root / "legacy-ambiguous-resolved"),
                                            record_cache_arg=str(root / "cache" / "records"))
            ledger = pd.read_csv(result["ledger"], sep="\t", dtype=str, keep_default_na=False)
            self.assertEqual("provider_unavailable", ledger.loc[0, "accession_resolution_status"])
            self.assertEqual("not_evaluated", ledger.loc[0, "evaluation_status"])
            self.assertIn("Q00001", ledger.loc[0, "provider_error"])


if __name__ == "__main__":
    unittest.main()
