"""Curated-function alternative path for the bait literature gate (issue #611 sibling).

A proposal whose accession carries an authoritative, source-backed curated
function (UniProt REST returns a function-bearing record) used to die as
``proposal_literature_is_not_verified`` whenever the abstracts of its real,
resolved papers lacked mechanistic cue language. The gate now accepts the
two-point combination: authoritative curated function + one provider-resolved
citation. Mechanistic support alone keeps passing; curated grounding without
any resolved citation still does not.
"""
from __future__ import annotations

import hashlib
import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from candidate_qtl.llm_pathway import hypothesis_verify

import test_hypothesis_verify_accession_resolution as harness


class CuratedFunctionGateTests(unittest.TestCase):
    def setUp(self) -> None:
        self.base = harness.AccessionResolutionTests()

    def _row(self, **overrides):
        row = {
            "hypothesis_id": "mech_001",
            "case_id": "c",
            "gene_family": "CYCLOPS/IPD3",
            "reference_gene_name": "",
            "reference_species": "Medicago truncatula",
            "proposed_stable_id": "P93003",
        }
        row.update(overrides)
        return row

    def test_curated_function_plus_resolved_citation_passes(self) -> None:
        """UniProt curated function + one resolved paper: two-point pass."""
        result = self.base._verify([self._row(
            verified_source_ids="PMID:12345",
            verified_citation_count="1",
            claim_support_summary="associative_support",
            mechanistic_support_citation_count="0",
            publication_source_provenance='[{"citation":"PMID:12345"}]',
            publication_records_sha256="AUTO",
        )], harness._api([]), cache_records=[
            harness._entry("P93003", "Medicago truncatula")
        ]).iloc[0]
        self.assertEqual("verified", result["verification_status"])
        self.assertEqual("", result["review_reason"])
        self.assertEqual("True", result["function_source_backed"])

    def test_associative_paper_alone_without_reviewed_record_still_fails(self) -> None:
        """A TrEMBL-tier record carries no authoritative curated function, so an
        associative paper alone still cannot pass the two-point combination."""
        entry = harness._entry("P93003", "Medicago truncatula", include_references=False)
        entry["entryType"] = "UniProtKB unreviewed (TrEMBL)"
        result = self.base._verify([self._row(
            verified_source_ids="PMID:12345",
            verified_citation_count="1",
            claim_support_summary="associative_support",
            mechanistic_support_citation_count="0",
            publication_source_provenance='[{"citation":"PMID:12345"}]',
            publication_records_sha256="AUTO",
        )], harness._api([]), cache_records=[entry]).iloc[0]
        self.assertEqual("review_required", result["verification_status"])
        self.assertEqual("proposal_literature_is_not_verified", result["review_reason"])

    def test_mechanistic_path_is_unchanged(self) -> None:
        result = self.base._verify([self._row(
            verified_source_ids="PMID:12345",
            verified_citation_count="1",
            claim_support_summary="mechanistic_support",
            mechanistic_support_citation_count="1",
            publication_source_provenance='[{"citation":"PMID:12345"}]',
            publication_records_sha256="AUTO",
        )], harness._api([]), cache_records=[
            harness._entry("P93003", "Medicago truncatula")
        ]).iloc[0]
        self.assertEqual("verified", result["verification_status"])

    def test_unresolved_provenance_blocks_the_curated_path(self) -> None:
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            records = {"PMID:12345": {
                "citation": "PMID:12345",
                "candidate_qtl_hypothesis_verifications": {"mech_001": {
                    "citation": "PMID:12345", "resolved_id": "",
                    "provider": "local_record_cache",
                    "url": "frozen-literature-cache.json",
                    "retrieved_at": "2026-08-22T00:00:00+00:00",
                    "claim_support_status": "review_required",
                }},
            }}
            provenance = dict(
                records["PMID:12345"]["candidate_qtl_hypothesis_verifications"]["mech_001"])
            rows = [self._row(
                verified_source_ids="PMID:12345",
                claim_support_summary="associative_support",
                mechanistic_support_citation_count="0",
                publication_source_provenance=json.dumps(provenance, sort_keys=True,
                                                         separators=(",", ":")),
                publication_records_sha256=hashlib.sha256(json.dumps(
                    [records["PMID:12345"]], ensure_ascii=False, sort_keys=True,
                    separators=(",", ":")).encode("utf-8")).hexdigest(),
            )]
            hypotheses = root / "hypotheses.tsv"
            extra = sorted(set().union(*(r.keys() for r in rows)) - set(
                ["hypothesis_id", "case_id", "gene_family", "reference_gene_name",
                 "reference_species", "proposed_stable_id"]))
            pd.DataFrame(rows, columns=["hypothesis_id", "case_id", "gene_family",
                                        "reference_gene_name", "reference_species",
                                        "proposed_stable_id"] + extra).to_csv(
                hypotheses, sep="\t", index=False)
            source_records_path = root / "source_records.json"
            source_records_path.write_text(json.dumps(records), encoding="utf-8")
            cache = root / "cache"
            cache.mkdir()
            (cache / "record.json").write_text(json.dumps(
                harness._entry("P93003", "Medicago truncatula")), encoding="utf-8")
            verified_dir = root / "verified"
            hypothesis_verify(str(hypotheses), str(verified_dir), "Pisum sativum",
                              record_cache_arg=str(cache), allow_external_api=False,
                              source_records_arg=str(source_records_path))
            validation = pd.read_csv(verified_dir / "proposal_validation.tsv",
                                     sep="\t", dtype=str).fillna("")
            row = validation.iloc[0]
            self.assertEqual("review_required", row["verification_status"])
            self.assertEqual("proposal_literature_is_not_verified", row["review_reason"])


if __name__ == "__main__":
    unittest.main()
