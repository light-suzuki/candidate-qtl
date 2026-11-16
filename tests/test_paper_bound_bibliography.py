"""Offline citation delivery through the real resolver and bibliography consumer."""
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.bibliography import INITIAL_POLICY
from candidate_qtl.cli import main
from candidate_qtl.paper_reports import PAPER_FIELD, normalize_paper_reports
from candidate_qtl.reference_gene_resolution import reference_gene_resolve
from tests.test_resolved_paper_reference_binding import paper, record


class PaperBoundBibliographyTests(unittest.TestCase):
    def verify(self, reports, *, shared_id=False, stable="", provider_species="Species alpha"):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name)
        records = [record(species=provider_species), record("Q00002", species=provider_species)]
        if not shared_id:
            records[1]["uniProtKBCrossReferences"][0]["id"] = "99999"
        cache = root / "records.json"
        cache.write_text(json.dumps({"results": records}))
        source = root / "input.tsv"
        pd.DataFrame([{"case_id": "neutral", "hypothesis_id": "h1", "reference_gene_name": "LN1",
                       "reference_species": "Species alpha", "proposed_stable_id": stable,
                       "source_ids": ";".join(report["source_id"] for report in reports),
                       "source_identifier_aliases": "PMID:888", "source_title": "Unbound top-level title",
                       PAPER_FIELD: normalize_paper_reports(reports)}]).to_csv(source, sep="\t", index=False)
        original_bytes = source.read_bytes()
        resolved = reference_gene_resolve(str(source), str(root / "resolve"), record_cache_arg=str(cache))
        resolved_path = Path(resolved["mechanisms_resolved"])
        resolved_bytes = resolved_path.read_bytes()
        before = pd.read_csv(resolved_path, sep="\t", dtype=str).fillna("")
        bibliography = root / "bibliography.json"
        bibliography.write_text(json.dumps({"records": [
            {"id": "123", "title": "Bound observation", "doi": "10.1234/bound"},
            {"id": "456", "title": "Unbound observation"},
            {"id": "888", "title": "Unbound top-level title"}]}))
        out = root / "source"
        self.assertEqual(0, main(["source-verify", "--mechanisms", str(resolved_path), "--out-dir", str(out),
                                  "--cache", str(bibliography), "--verification-policy", INITIAL_POLICY]))
        after = pd.read_csv(out / "mechanisms_with_verified_sources.tsv", sep="\t", dtype=str).fillna("")
        self.assertEqual(original_bytes, source.read_bytes())
        self.assertEqual(resolved_bytes, resolved_path.read_bytes())
        for field in (PAPER_FIELD, "retrieval_provenance"):
            self.assertEqual(before[field].tolist(), after[field].tolist())
        return after.set_index("proposed_stable_id"), pd.read_csv(out / "source_validation.tsv", sep="\t").fillna("")

    def test_only_uniquely_bound_paper_and_its_alias_reach_bibliography(self):
        bound = {**paper("GeneID:12345"), "source_title": "Bound observation",
                 "source_identifier_aliases": "DOI:10.1234/bound"}
        unbound = {**paper(), "source_id": "PMID:456", "source_title": "Unbound observation"}
        rows, validation = self.verify([bound, unbound])
        self.assertEqual("PMID:123", rows.loc["Q00001", "source_ids"])
        self.assertEqual("resolved", rows.loc["Q00001", "bibliography_state"])
        self.assertEqual({"PMID:123", "DOI:10.1234/bound"}, set(validation.citation))
        self.assertEqual({"id_title_matched"}, set(validation.bibliography_state))
        self.assertEqual({"not_evaluated"}, set(validation.claim_support_status))
        self.assertEqual("", rows.loc["Q00002", "source_ids"])
        self.assertEqual("reference_endpoint_changed", rows.loc["Q00002", "bibliography_state"])
        receipt, = json.loads(rows.loc["Q00001", "retrieval_provenance"])["paper_reference_rebindings"]
        self.assertEqual("provider_record_cross_reference", receipt["basis"])
        reports = json.loads(rows.loc["Q00001", PAPER_FIELD])
        self.assertEqual(bound["paper_checklist"], json.loads(reports[0]["paper_checklist"]))
        self.assertEqual(unbound["paper_checklist"], json.loads(reports[1]["paper_checklist"]))

    def test_mixed_reports_do_not_union_unbound_citations_titles_or_aliases(self):
        bound = {**paper("GeneID:12345"), "source_title": "Bound observation",
                 "source_identifier_aliases": "DOI:10.1234/bound"}
        cases = [(paper(), {}), (paper("GeneID:99999"), {}), (paper("UniProt:Q00003"), {}),
                 (paper("GeneID:12345", species="Species beta"), {}),
                 (paper("GeneID:12345"), {"shared_id": True})]
        for other, options in cases:
            with self.subTest(other=other, options=options):
                bound_report = {**bound, "reference_gene_id": "UniProt:Q00001"} if options else bound
                other = {**other, "source_id": "PMID:456", "source_title": "Unbound observation",
                         "source_identifier_aliases": "PMID:888"}
                rows, validation = self.verify([bound_report, other], **options)
                target = rows.loc["Q00001"]
                self.assertEqual("PMID:123", target.source_ids)
                self.assertEqual("resolved", target.bibliography_state)
                citations = validation[validation.hypothesis_id == target.hypothesis_id]
                self.assertEqual({"PMID:123", "DOI:10.1234/bound"}, set(citations.citation))
                self.assertEqual({"id_title_matched"}, set(citations.bibliography_state))
                self.assertEqual(["PMID:123", "PMID:456"],
                                 json.loads(target.retrieval_provenance)["original_source_ids"])
                self.assertEqual(2, len(json.loads(target[PAPER_FIELD])))

    def test_homolog_report_keeps_its_distinct_reported_species(self):
        report = {**paper("GeneID:12345", species="Species beta"), "gene_correspondence": "homolog",
                  "source_title": "Bound observation"}
        rows, _ = self.verify([report])
        self.assertEqual("resolved", rows.loc["Q00001", "bibliography_state"])
        self.assertEqual("Species beta", json.loads(rows.loc["Q00001", PAPER_FIELD])[0]["reported_species"])

    def test_unresolved_or_conflicting_paper_cannot_restore_citations(self):
        cases = [(paper(), {}), (paper("GeneID:12345"), {"shared_id": True}),
                 (paper("GeneID:77777"), {}), (paper("UniProt:Q00003"), {}),
                 (paper("GeneID:12345", species="Species beta"), {}),
                 (paper("GeneID:12345"), {"provider_species": "Species beta"})]
        for report, options in cases:
            with self.subTest(report=report, options=options):
                rows, validation = self.verify([report], **options)
                self.assertEqual({""}, set(rows.source_ids))
                self.assertTrue(validation.empty)

    def test_explicit_parent_endpoint_cannot_be_overridden_by_paper_rebinding(self):
        rows, validation = self.verify([paper("GeneID:12345")], stable="Q00002")
        self.assertEqual("", rows.loc["Q00001", "source_ids"])
        self.assertEqual("reference_endpoint_changed", rows.loc["Q00001", "bibliography_state"])
        self.assertNotIn(rows.loc["Q00001", "hypothesis_id"], set(validation.hypothesis_id))


if __name__ == "__main__":
    unittest.main()
