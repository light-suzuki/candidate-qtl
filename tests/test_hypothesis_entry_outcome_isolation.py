"""Fixed Route A provider failures stay attached to the requesting hypothesis."""

from __future__ import annotations

import json
import tempfile
import unittest
import urllib.error
from io import BytesIO
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl.bibliography import INITIAL_POLICY
from candidate_qtl.llm_pathway import (
    _ORIGINAL_HYPOTHESIS_VERIFY,
    bait_panel_freeze,
    hypothesis_verify,
)
from candidate_qtl.phenotype_investigation import mechanism_import
from candidate_qtl.reference_gene_resolution import reference_gene_resolve
from candidate_qtl.source_verify import source_verify


def _entry() -> dict:
    return {
        "primaryAccession": "P12345",
        "organism": {"scientificName": "Arabidopsis thaliana"},
        "sequence": {"value": "MAAGSVDPLA" * 12},
        "genes": [{"geneName": {"value": "REF1"}}],
        "references": [{"citation": {"citationCrossReferences": [
            {"database": "PubMed", "id": "12345"},
        ]}}],
    }


def _hypothesis(hypothesis_id: str, accession: str = "P12345") -> dict:
    return {
        "hypothesis_id": hypothesis_id,
        "pathway_or_process": "meristem growth",
        "gene_family": "reference family",
        "reference_gene_name": "",
        "reference_species": "Arabidopsis thaliana",
        "proposed_stable_id": accession,
        "expected_role": "meristem proliferation",
        "mechanism_statement": "Reduced meristem proliferation reduces branching",
        "source_ids": "PMID:12345",
    }


def _prepare(root: Path, rows: list[dict]) -> str:
    payload = root / "payload.json"
    payload.write_text(json.dumps({"hypotheses": rows}), encoding="utf-8")
    imported = mechanism_import(str(payload), str(root / "import"), case_id="isolation")
    resolved = reference_gene_resolve(
        imported["mechanisms"], str(root / "resolve"), focal_species="Pisum sativum",
    )
    publications = root / "publications.json"
    publications.write_text(json.dumps({
        "id": "12345", "pmid": "12345", "title": "A synthetic meristem study",
    }), encoding="utf-8")
    sources = source_verify(
        resolved["mechanisms_resolved"], str(root / "sources"),
        cache_arg=str(publications), verification_policy=INITIAL_POLICY,
    )
    return sources["mechanisms_with_verified_sources"]


class HypothesisEntryOutcomeIsolationTests(unittest.TestCase):
    def _verify_and_freeze(self, root, hypotheses, replies, cache=None):
        with patch("candidate_qtl.uniprot_error_policy.urllib.request.urlopen", side_effect=replies) as fetch:
            verified = hypothesis_verify(
                hypotheses, str(root / "verified"), "Pisum sativum",
                record_cache_arg=str(cache) if cache else None, allow_external_api=True,
            )
        frozen = bait_panel_freeze(
            hypotheses, verified["validation"], verified["verified_sequences"], str(root / "frozen"),
            case_id="isolation", reviewable_sequences_arg=verified["reviewable_sequences"],
        )
        validation = pd.read_csv(verified["validation"], sep="\t", dtype=str).fillna("").set_index("hypothesis_id")
        projection = pd.read_csv(frozen["projection_panel"], sep="\t", dtype=str).fillna("")
        manifest = json.loads(Path(verified["manifest"]).read_text(encoding="utf-8"))
        self.assertEqual({"healthy"}, set(projection.hypothesis_id))
        self.assertEqual({"review_only"}, set(projection.projection_eligibility))
        self.assertEqual({"healthy"}, set(json.loads(Path(verified["reviewable_sequences"]).read_text())))
        self.assertEqual("verified", validation.loc["healthy", "identity_sequence_verification_status"])
        self.assertEqual("review_required", validation.loc["healthy", "verification_status"])
        self.assertEqual("", validation.loc["healthy", "rejection_reason"])
        self.assertEqual("False", validation.loc["healthy", "bound_claim_authorized"])
        self.assertEqual("not_verified", validation.loc["failed", "identity_sequence_verification_status"])
        return validation, manifest, fetch.call_count

    def test_same_accession_failure_never_rewrites_an_independent_healthy_hypothesis(self):
        failures = [
            (400, "invalid_http_400", "rejected", "invalid_submitted_accession"),
            (404, "not_found_http_404", "rejected", "accession_not_found_in_authoritative_database"),
            (429, "rate_limited", "review_required", ""),
            (500, "provider_unavailable", "review_required", ""),
            (None, "provider_unavailable", "review_required", ""),
        ]
        for code, state, status, rejection in failures:
            for healthy_first in (True, False):
                with self.subTest(code=code, healthy_first=healthy_first), tempfile.TemporaryDirectory() as tmp:
                    root = Path(tmp)
                    order = ["healthy", "failed"] if healthy_first else ["failed", "healthy"]
                    hypotheses = _prepare(root, [_hypothesis(value) for value in order])
                    failure = (urllib.error.HTTPError("https://rest.uniprot.org/uniprotkb/P12345.json",
                               code, "fixed error", None, None) if code else TimeoutError("fixed timeout"))
                    replies = [BytesIO(json.dumps(_entry()).encode()) if value == "healthy" else failure for value in order]
                    validation, manifest, calls = self._verify_and_freeze(root, hypotheses, replies)
                    self.assertEqual(2, calls)
                    self.assertEqual(state, validation.loc["failed", "accession_resolution_status"])
                    self.assertEqual(status, validation.loc["failed", "verification_status"])
                    self.assertEqual(rejection, validation.loc["failed", "rejection_reason"])
                    self.assertEqual({"retrieved": 1, state: 1}, manifest["uniprot_entry_outcome_counts"])
                    if status == "review_required":
                        self.assertEqual("accession_provider_unavailable", validation.loc["failed", "review_reason"])

    def test_network_failure_for_an_alias_does_not_reject_the_healthy_cache_row(self):
        for healthy_first in (True, False):
            with self.subTest(healthy_first=healthy_first), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                rows = [_hypothesis("healthy"), _hypothesis("failed", "UniProt:P12345")]
                hypotheses = _prepare(root, rows if healthy_first else list(reversed(rows)))
                cache = root / "records"
                cache.mkdir()
                (cache / "P12345.json").write_text(json.dumps(_entry()), encoding="utf-8")
                failure = urllib.error.HTTPError("https://rest.uniprot.org/uniprotkb/P12345.json", 404,
                                                 "fixed missing alias response", None, None)
                validation, manifest, calls = self._verify_and_freeze(root, hypotheses, [failure], cache)
                self.assertEqual(1, calls)
                self.assertEqual("local_record_cache", validation.loc["healthy", "source_provider"])
                self.assertEqual("not_found_http_404", validation.loc["failed", "accession_resolution_status"])
                self.assertEqual({"not_found_http_404": 1}, manifest["uniprot_entry_outcome_counts"])

    def test_unversioned_cache_only_validation_keeps_its_existing_shape_and_decisions(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hypotheses = root / "hypotheses.tsv"
            pd.DataFrame([_hypothesis("healthy")]).to_csv(hypotheses, sep="\t", index=False)
            cache = root / "cache"
            cache.mkdir()
            (cache / "record.json").write_text(json.dumps(_entry()), encoding="utf-8")
            with patch("candidate_qtl.llm_pathway_legacy._now", return_value="2026-10-19T00:00:00Z"):
                original = _ORIGINAL_HYPOTHESIS_VERIFY(str(hypotheses), str(root / "legacy"), "Pisum sativum",
                                                      record_cache_arg=str(cache))
                actual = hypothesis_verify(str(hypotheses), str(root / "actual"), "Pisum sativum",
                                           record_cache_arg=str(cache))
            for key in ("validation", "rejected", "review_required", "verified_sequences", "reviewable_sequences"):
                self.assertEqual(Path(original[key]).read_bytes(), Path(actual[key]).read_bytes())
            self.assertEqual(1, actual["verified_count"])
            manifest = json.loads(Path(actual["manifest"]).read_text(encoding="utf-8"))
            self.assertNotIn("uniprot_entry_outcome_counts", manifest)


if __name__ == "__main__":
    unittest.main()
