"""Original Route-A citations are review material, never provider evidence.

Synthetic offline contract fixtures; no biological benchmark or model result.
"""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.llm_pathway import _record_source_ids
from candidate_qtl.reference_gene_resolution import reference_gene_resolve
from candidate_qtl.source_verify import source_verify
from tests.test_agent_claim_review import ReviewFixture, write_json
from tests.test_reference_gene_resolution import (
    _mechanism_row, _uniprot_record, _write_cache, _write_mechanisms,
)


class OriginalRouteASourceTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def resolve(self, accessions, supplied="P12345"):
        mechanisms = self.root / "mechanisms.tsv"
        _write_mechanisms(mechanisms, [_mechanism_row(
            "h1", source_ids="PMID:12345", proposed_stable_id=supplied,
            reference_gene_name="XYZ1", reference_species="Arabidopsis thaliana")])
        cache = self.root / "records.json"
        _write_cache(cache, [_uniprot_record(
            accession, "Arabidopsis thaliana", "MAAA", pmid=pmid,
            gene_name="XYZ1") for accession, pmid in accessions])
        original = mechanisms.read_bytes()
        result = reference_gene_resolve(str(mechanisms), str(self.root / "resolved"),
                                        record_cache_arg=str(cache))
        self.assertEqual(original, mechanisms.read_bytes())
        rows = pd.read_csv(result["mechanisms_resolved"], sep="\t", dtype=str).fillna("")
        ledger = pd.read_csv(result["ledger"], sep="\t", dtype=str).fillna("")
        return result, rows, ledger

    def test_same_accession_keeps_original_citation_with_separate_provenance(self):
        _, rows, ledger = self.resolve([("P12345", "54321")])
        row = rows.iloc[0]
        self.assertEqual({"PMID:12345", "PMID:54321", "UniProt:P12345"}, set(row.source_ids.split(";")))
        provenance = json.loads(row.retrieval_provenance)
        self.assertEqual(["PMID:12345"], provenance["original_source_ids"])
        self.assertEqual(["PMID:54321", "UniProt:P12345"], provenance["provider_linked_source_ids"])
        self.assertEqual("same_record_review_material", provenance["original_source_retention_state"])
        self.assertEqual("P12345", provenance["input_proposed_stable_id"])
        self.assertEqual(row.source_ids, ledger.iloc[0].resolved_source_ids)
        self.assertEqual("pending_verification", row.review_status)

    def test_changed_accession_does_not_inherit_original_citation(self):
        _, rows, ledger = self.resolve([("Q12345", "54321")])
        row = rows.iloc[0]
        self.assertEqual("Q12345", row.proposed_stable_id)
        self.assertEqual("PMID:54321;UniProt:Q12345", row.source_ids)
        self.assertEqual("PMID:12345", ledger.iloc[0].input_source_ids)
        provenance = json.loads(row.retrieval_provenance)
        self.assertEqual(["PMID:12345"], provenance["original_source_ids"])
        self.assertEqual("audit_only_identity_not_confirmed", provenance["original_source_retention_state"])

    def test_exact_same_accession_without_provider_bibliography_keeps_review_citation(self):
        _, rows, ledger = self.resolve([("P12345", "")])
        row = rows.iloc[0]
        self.assertEqual("P12345", row.proposed_stable_id)
        self.assertEqual("exact_name_species", row.reference_identity_status)
        self.assertEqual("PMID:12345", row.source_ids)
        provenance = json.loads(row.retrieval_provenance)
        self.assertEqual([], provenance["provider_linked_source_ids"])
        self.assertEqual("same_record_review_material", provenance["original_source_retention_state"])
        self.assertEqual(row.source_ids, ledger.iloc[0].resolved_source_ids)
        self.assertEqual("pending_verification", row.review_status)

    def test_changed_accession_without_provider_bibliography_does_not_borrow_citation(self):
        _, rows, ledger = self.resolve([("Q12345", "")])
        self.assertEqual("Q12345", rows.iloc[0].proposed_stable_id)
        self.assertEqual("", rows.iloc[0].source_ids)
        self.assertEqual("PMID:12345", ledger.iloc[0].input_source_ids)

    def test_only_the_confirmed_same_record_child_keeps_original_citation(self):
        _, rows, _ = self.resolve([("P12345", "54321"), ("Q12345", "67890")])
        self.assertEqual(2, len(rows))
        self.assertEqual({"PMID:12345;PMID:54321;UniProt:P12345", "PMID:67890;UniProt:Q12345"}, set(rows.source_ids))
        for row in rows.to_dict("records"):
            provenance = json.loads(row["retrieval_provenance"])
            self.assertEqual(["PMID:12345"], provenance["original_source_ids"])
            same_record = row["proposed_stable_id"] == "P12345"
            self.assertEqual("same_record_review_material" if same_record else "audit_only_identity_not_confirmed",
                             provenance["original_source_retention_state"])
            self.assertEqual(same_record, "PMID:12345" in row["source_ids"].split(";"))
            self.assertEqual(row["proposed_stable_id"], provenance["validated_reference_accession"])

    def test_missing_supplied_accession_does_not_bind_original_citation_to_resolved_name(self):
        _, rows, _ = self.resolve([("P12345", "54321")], supplied="")
        self.assertEqual("PMID:54321;UniProt:P12345", rows.iloc[0].source_ids)

    def test_retained_source_reaches_independent_review_without_bypassing_binding(self):
        self.check_retained_source_delivery(with_provider_bibliography=True)

    def test_source_less_provider_still_delivers_original_text_to_independent_review(self):
        self.check_retained_source_delivery(with_provider_bibliography=False)

    def check_retained_source_delivery(self, *, with_provider_bibliography):
        fixture = ReviewFixture(self.root)
        frame = pd.read_csv(fixture.hypotheses, sep="\t", dtype=str).fillna("")
        frame["case_id"] = "source-provenance-contract"
        frame.to_csv(fixture.hypotheses, sep="\t", index=False)
        if with_provider_bibliography:
            fixture.record["references"][0]["citation"]["citationCrossReferences"][0]["id"] = "54321"
        else:
            fixture.record["references"] = []
        write_json(fixture.records_path, fixture.record)
        result = reference_gene_resolve(str(fixture.hypotheses), str(self.root / "resolved"),
                                        record_cache_arg=str(fixture.records_path))
        source = source_verify(result["mechanisms_resolved"], str(self.root / "source-after-resolve"),
                               cache_arg=str(fixture.cache))
        fixture.hypotheses = Path(source["mechanisms_with_verified_sources"])
        fixture.sources_path = Path(source["source_review_records"])
        fixture.source_records = Path(source["verified_sources"])
        fixture.seal()
        prompt = json.loads(fixture.prompt_path.read_text())
        captured = prompt["source_review_records"]["records"]["PMID:12345"]
        self.assertEqual(["h1"], captured["hypothesis_ids"])
        self.assertEqual(fixture.text, captured["sections"]["abstractText"])
        self.assertIn("PMID:12345", prompt["hypotheses"][0]["source_ids"])
        bound = fixture.checker().bind({"hypothesis_id": "h1"}, "P12345", fixture.record,
                                       _record_source_ids(fixture.record))
        self.assertEqual("review_required", bound["binding_state"])
        self.assertEqual("symbol-to-accession binding lacks a provider citation link", bound["binding_reason"])
        verified = fixture.verify(out="verify-after-resolve")
        self.assertEqual(0, verified["verified_count"])


if __name__ == "__main__":
    unittest.main()
