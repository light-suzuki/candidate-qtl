"""Offline policy wiring fixtures; not scientific accuracy measurements."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl.bibliography import INITIAL_POLICY
from candidate_qtl.source_verify import source_verify
from candidate_qtl.route_b_semantic_source_verify import verify_route_b_semantic_sources
from candidate_qtl.route_b_semantic_transfer_verify import verify_route_b_semantic_transfer
from candidate_qtl.route_b_semantic_candidates import freeze_route_b_semantic_candidates
from tests.test_route_b_semantic_transfer_verify import _claim, _bundle


class InitialBibliographyTests(unittest.TestCase):
    def test_initial_production_prompt_has_no_pending_semantic_review_gate(self):
        from candidate_qtl.route_b_semantic_bundle import _prompt_record
        unit = {"unit_id": "u", "genes": []}
        initial = _prompt_record(unit, "flowering", "case")["system"]
        legacy = _prompt_record(unit, "flowering", "case", "source_function_v1")["system"]
        self.assertIn("adopted as a first-stage candidate", initial)
        self.assertIn("initial_tier_reason", initial)
        self.assertNotIn("until source/identity checks and semantic review", initial)
        self.assertNotIn("source_function_v1 bundles", initial)
        self.assertIn("until source/identity checks and semantic review", legacy)

    def test_id_title_aliases_missing_and_failures_are_separate_from_support(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cache = root / "cache.json"
            cache.write_text(json.dumps([{"id": "1234567", "pmid": "1234567", "doi": "10.1234/paper",
                                         "title": "A function: in plants", "abstractText": "must not be read"},
                                        {"id": "7654321", "title": "Another paper"}]))
            cases = [
                {"source_ids": "PMID:1234567", "source_title": "A FUNCTION—in plants."},
                {"source_ids": "PMID:1234567", "source_title": "Unrelated study"},
                {"source_ids": "PMID:1234567"},
                {"source_ids": "PMID:1234567", "source_identifier_aliases": "doi:10.1234/paper"},
                {"source_ids": "PMID:1234567", "source_identifier_aliases": "PMID:7654321"},
                {"source_ids": "PMID:9999999"}, {"source_ids": ""},
            ]
            table = root / "a.tsv"
            pd.DataFrame([{"hypothesis_id": str(i), **r} for i, r in enumerate(cases)]).fillna("").to_csv(table, sep="\t", index=False)
            with patch("candidate_qtl.source_verify._claim_text", side_effect=AssertionError("content was read")):
                result = source_verify(str(table), str(root / "out"), cache_arg=str(cache), verification_policy=INITIAL_POLICY)
            val = pd.read_csv(result["source_validation"], sep="\t", dtype=str).fillna("")
            self.assertEqual(["id_title_matched", "title_mismatch", "id_resolved_title_not_compared",
                              "id_resolved_title_not_compared", "multiple_identifiers_mismatch", "not_attempted"],
                             [val[val.hypothesis_id == str(i)].iloc[0].bibliography_state for i in range(6)])
            self.assertEqual({"not_evaluated"}, set(val.claim_support_status))
            self.assertNotIn("abstractText", Path(result["bibliography_records"]).read_text())
            self.assertEqual({}, json.loads(Path(result["source_review_records"]).read_text())["records"])
            annotated = pd.read_csv(result["mechanisms_with_verified_sources"], sep="\t", dtype=str).fillna("")
            self.assertEqual("source_absent", annotated.iloc[6].bibliography_state)

    def test_resolved_a_checks_original_papers_without_provider_bibliography_or_other_copy(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cache = root / "cache.json"
            cache.write_text(json.dumps([
                {"id": "1234567", "pmid": "1234567", "doi": "10.1234/first", "title": "First observation"},
                {"id": "7654321", "pmid": "7654321", "title": "Another study of the same process"}]))
            provenance = {"original_source_ids": ["PMID:1234567", "PMID:7654321"],
                          "provider_linked_source_ids": ["PMID:9999999"],
                          "input_proposed_stable_id": "Q96524",
                          "record_identity_state": "resolved_exact_name_species", "candidate_count": 2}
            base = {"source_ids": "PMID:1234567;PMID:7654321;PMID:9999999",
                    "source_title": "First observation",
                    "source_identifier_aliases": "PMID:1234567;DOI:10.1234/first",
                    "retrieval_provenance": json.dumps(provenance)}
            table = root / "a.tsv"
            pd.DataFrame([{**base, "hypothesis_id": "declared", "proposed_stable_id": "Q96524"},
                          {**base, "hypothesis_id": "other-copy", "proposed_stable_id": "Q96525"}]).to_csv(table, sep="\t", index=False)
            result = source_verify(str(table), str(root / "out"), cache_arg=str(cache), verification_policy=INITIAL_POLICY)
            annotated = pd.read_csv(result["mechanisms_with_verified_sources"], sep="\t").fillna("").set_index("hypothesis_id")
            self.assertEqual("resolved", annotated.loc["declared", "bibliography_state"])
            self.assertEqual("PMID:1234567;PMID:7654321", annotated.loc["declared", "source_ids"])
            self.assertEqual("PMID:9999999", annotated.loc["declared", "reference_record_source_ids"])
            self.assertEqual("reference_endpoint_changed", annotated.loc["other-copy", "bibliography_state"])
            self.assertEqual("", annotated.loc["other-copy", "verified_source_ids"])
            validation = pd.read_csv(result["source_validation"], sep="\t").set_index("citation")
            self.assertEqual("id_resolved_title_not_compared", validation.loc["PMID:7654321", "bibliography_state"])
            self.assertNotIn("PMID:9999999", validation.index)

    def test_malformed_declared_accession_is_undeclared_not_a_changed_endpoint(self):
        """Audit text such as "unknown (...)" is not an accession the researcher declared.

        With one resolved record the original papers stay checkable; with several
        records the paper-bound contract (#752) still withholds them.
        """
        for candidate_count, expected_state, expected_sources in [
                (1, "resolved", "PMID:1234567"), (2, "reference_endpoint_changed", "")]:
            with self.subTest(candidate_count=candidate_count), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                cache = root / "cache.json"
                cache.write_text(json.dumps([{"id": "1234567", "pmid": "1234567", "title": "First observation"}]))
                provenance = {"original_source_ids": ["PMID:1234567"],
                              "provider_linked_source_ids": ["PMID:9999999"],
                              "input_proposed_stable_id": "unknown (accession not verified)",
                              "record_identity_state": "resolved_exact_name_species",
                              "candidate_count": candidate_count}
                table = root / "a.tsv"
                pd.DataFrame([{"hypothesis_id": "h", "proposed_stable_id": "Q96524",
                               "source_ids": "PMID:1234567;PMID:9999999", "source_title": "First observation",
                               "retrieval_provenance": json.dumps(provenance)}]).to_csv(table, sep="\t", index=False)
                result = source_verify(str(table), str(root / "out"), cache_arg=str(cache),
                                       verification_policy=INITIAL_POLICY)
                annotated = pd.read_csv(result["mechanisms_with_verified_sources"], sep="\t").fillna("")
                self.assertEqual(expected_state, annotated.iloc[0].bibliography_state)
                self.assertEqual(expected_sources, annotated.iloc[0].verified_source_ids)

    def test_wrong_returned_id_no_hit_and_429(self):
        for returned, expected in [({"id": "7654321", "title": "Wrong hit"}, "identifier_mismatch"),
                                   (None, "not_found"), (RuntimeError("bibliography provider unavailable: HTTP 429"), "provider_failure")]:
            with self.subTest(expected=expected), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                table = root / "a.tsv"
                pd.DataFrame([{"hypothesis_id": "a", "source_ids": "PMID:1234567"}]).to_csv(table, sep="\t", index=False)
                with patch("candidate_qtl.bibliography.fetch_metadata", side_effect=returned if isinstance(returned, Exception) else None,
                           return_value=returned):
                    result = source_verify(str(table), str(root / "out"), allow_external_api=True, verification_policy=INITIAL_POLICY)
                self.assertEqual(expected, pd.read_csv(result["source_validation"], sep="\t").iloc[0].bibliography_state)

    def test_b_initial_candidate_reaches_freeze_without_content_review_and_preserves_copy_endpoint(self):
        for function, source, expected in [("meristem maintenance", "PMID:1234567", {"G1"}),
                                            ("kinase", "PMID:1234567", set()),
                                            ("meristem maintenance", "", set()),
                                            ("meristem maintenance", "PMID:9999999", set())]:
            with self.subTest(function=function, source=source), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                bundle = root / "bundle.json"
                _bundle(bundle, members=[("G1", "1:many"), ("G2", "1:many")])
                value = json.loads(bundle.read_text())
                value["b_review_membership_policy"] = INITIAL_POLICY
                value["units"][0]["genes"][1]["orthology_reference_gene"] = "OTHER_REF"
                bundle.write_text(json.dumps(value))
                proposals = root / "proposals.tsv"
                pd.DataFrame([_claim(source_id=source, biological_function=function)]).to_csv(proposals, sep="\t", index=False)
                cache = root / "cache.json"
                cache.write_text(json.dumps({"id": "1234567", "title": "Synthetic meristem study"}))
                sources = verify_route_b_semantic_sources(str(proposals), str(root / "source"), bundle_arg=str(bundle), cache_arg=str(cache))
                transfer = verify_route_b_semantic_transfer(str(proposals), sources["source_verification"], str(bundle), str(root / "transfer"), case_id="case403")
                frozen = freeze_route_b_semantic_candidates(transfer["claims"], str(root / "freeze"), case_id="case403")
                review = pd.read_csv(frozen["review_candidates"], sep="\t").fillna("")
                self.assertEqual(expected, set(review.candidate_gene_id))
                self.assertEqual(0, frozen["candidate_gene_count"])
                if expected:
                    self.assertEqual({"adopted_first_stage"}, set(review.candidate_admission_state))
                    self.assertEqual({"review"}, set(review.support_state))
                    source_file = Path(sources["source_verification"])
                    source_file.write_bytes(source_file.read_bytes() + b"\n")
                    with self.assertRaisesRegex(ValueError, "hash mismatch"):
                        verify_route_b_semantic_transfer(str(proposals), str(source_file), str(bundle), str(root / "tampered"))

    def test_curated_function_original_judgment_is_adopted_without_independent_review(self):
        from tests.test_route_b_semantic_source_verify import _bundle as curated_bundle, _proposal
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bundle = curated_bundle(root)
            value = json.loads(bundle.read_text())
            value["b_review_membership_policy"] = INITIAL_POLICY
            bundle.write_text(json.dumps(value))
            proposals = root / "proposals.tsv"
            pd.DataFrame([_proposal(source_id="UniProt:REF1", claim_type="direct_functional_bridge")]).to_csv(proposals, sep="\t", index=False)
            sources = verify_route_b_semantic_sources(str(proposals), str(root / "sources"), bundle_arg=str(bundle))
            transfer = verify_route_b_semantic_transfer(str(proposals), sources["source_verification"], str(bundle), str(root / "transfer"), case_id="case403")
            frozen = freeze_route_b_semantic_candidates(transfer["claims"], str(root / "freeze"), case_id="case403")
            self.assertEqual({"G1"}, set(pd.read_csv(frozen["review_candidates"], sep="\t").candidate_gene_id))
            self.assertEqual(0, frozen["candidate_gene_count"])


if __name__ == "__main__":
    unittest.main()
