"""Offline reference transport fixtures; no new paper or biological judgment."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from candidate_qtl.observation_tiers import classify_observation
from candidate_qtl.paper_reports import bound_paper_rows
from candidate_qtl.reference_gene_resolution import reference_gene_resolve


def record(accession="Q00001", species="Species alpha", name="LN1"):
    return {"primaryAccession": accession, "organism": {"scientificName": species},
            "genes": [{"geneName": {"value": name}}], "sequence": {"value": "M" * 30},
            "uniProtKBCrossReferences": [{"database": "GeneID", "id": "12345"}],
            "references": [{"citation": {"citationCrossReferences": [{"database": "PubMed", "id": "123"}]}}]}


def paper(anchor="LONG NODE (LN1)", species="Species alpha"):
    return {"source_id": "PMID:123", "reported_gene": "LN1", "reported_species": species,
            "reference_gene_id": anchor, "gene_correspondence": "same_gene",
            "initial_tier_reason": "Original report reason retained.",
            "paper_checklist": {"focal_species": "no", "same_gene_as_focal": "unknown",
                                "phenotype_reported": "yes", "equivalent_target_phenotype": "unknown",
                                "similar_cellular_or_developmental_phenotype": "yes", "shared_specific_process": "yes",
                                "function_supports_hypothesis": "yes", "annotation_supports_hypothesis": "unknown",
                                "identity_resolved": "yes"}}


class ResolvedPaperReferenceBindingTests(unittest.TestCase):
    def resolve(self, reports, records=None, name="LONG NODE (LN1)", species="Species alpha", stable="", allow_external_api=False):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name)
        source = root / "input.tsv"
        pd.DataFrame([{"case_id": "neutral", "hypothesis_id": "h001", "reference_gene_name": name,
                       "reference_species": species, "source_ids": "PMID:123", "proposed_stable_id": stable,
                       "paper_reports_json": json.dumps(reports)}]).to_csv(source, sep="\t", index=False)
        original = source.read_bytes()
        cache = root / "records.json"
        cache.write_text(json.dumps({"results": records if records is not None else [record()]}), encoding="utf-8")
        result = reference_gene_resolve(str(source), str(root / "resolved"), record_cache_arg=str(cache),
                                        allow_external_api=allow_external_api)
        self.assertEqual(original, source.read_bytes())
        return pd.read_csv(result["mechanisms_resolved"], sep="\t", dtype=str).fillna("").to_dict("records")

    def test_original_display_label_reaches_paper_tier_without_rewriting_answers(self):
        original = paper()
        row, = self.resolve([original])
        resolved, = json.loads(row["paper_reports_json"])
        self.assertEqual({**original, "reference_gene_id": "Q00001"}, resolved)
        bound, = bound_paper_rows(row)
        self.assertEqual("matched", bound["tier_gene_binding"])
        self.assertEqual("T2", classify_observation(bound)[0])
        self.assertEqual("unknown", resolved["paper_checklist"]["same_gene_as_focal"])
        receipt, = json.loads(row["retrieval_provenance"])["paper_reference_rebindings"]
        self.assertEqual("LONG NODE (LN1)", receipt["input_reference_gene_id"])
        self.assertEqual("provider_record_original_atomic_label", receipt["basis"])

    def test_declared_provider_alias_and_stable_cross_reference(self):
        for anchor, basis in [("LN1", "provider_record_declared_alias"),
                              ("GeneID:12345", "provider_record_cross_reference")]:
            with self.subTest(anchor=anchor):
                row, = self.resolve([paper(anchor)])
                self.assertEqual("Q00001", json.loads(row["paper_reports_json"])[0]["reference_gene_id"])
                self.assertEqual(basis, json.loads(row["retrieval_provenance"])["paper_reference_rebindings"][0]["basis"])

    def test_arbitrary_label_other_reference_and_conflicting_stable_ids_stay_unchanged(self):
        for anchor in ["INVENTED (LN1)", "OTHER", "GeneID:99999", "UniProt:Q00002"]:
            with self.subTest(anchor=anchor):
                original = paper(anchor)
                row, = self.resolve([original])
                self.assertEqual([original], json.loads(row["paper_reports_json"]))
                self.assertNotIn("paper_reference_rebindings", json.loads(row["retrieval_provenance"]))
        original = paper()
        row, = self.resolve([original], stable="Q00002")
        self.assertEqual([original], json.loads(row["paper_reports_json"]))

    def test_conflicting_id_in_original_display_label_is_not_a_valid_alias(self):
        original = paper("LONG NODE (UniProt:Q00002)")
        actual = record(name="LONG NODE")
        row, = self.resolve([original], records=[actual], name=original["reference_gene_id"])
        self.assertEqual([original], json.loads(row["paper_reports_json"]))

    def test_name_only_copy_ambiguity_stays_unknown_but_explicit_id_can_bind_one_copy(self):
        rows = self.resolve([paper(), paper("UniProt:Q00001"), paper("GeneID:12345")], records=[record(), record("Q00002")])
        self.assertEqual(2, len(rows))
        for row in rows:
            reports = json.loads(row["paper_reports_json"])
            self.assertEqual("LONG NODE (LN1)", reports[0]["reference_gene_id"])
            self.assertEqual("Q00001" if row["proposed_stable_id"] == "Q00001" else "UniProt:Q00001",
                             reports[1]["reference_gene_id"])
            self.assertEqual("GeneID:12345", reports[2]["reference_gene_id"])

    def test_same_name_in_other_species_does_not_receive_paper_binding(self):
        original = paper(species="Species beta")
        row, = self.resolve([original], records=[record(), record("Q00002", "Species beta")])
        self.assertEqual("Q00001", row["proposed_stable_id"])
        self.assertEqual([original], json.loads(row["paper_reports_json"]))

    def test_unique_external_id_distinguishes_same_name_provider_records(self):
        other = record("Q00002")
        other["uniProtKBCrossReferences"][0]["id"] = "99999"
        original = paper("GeneID:12345")
        rows = self.resolve([original], records=[record(), other])
        self.assertEqual(2, len(rows))
        for row in rows:
            expected = "Q00001" if row["proposed_stable_id"] == "Q00001" else "GeneID:12345"
            resolved, = json.loads(row["paper_reports_json"])
            self.assertEqual({**original, "reference_gene_id": expected}, resolved)
            self.assertEqual("unknown", resolved["paper_checklist"]["same_gene_as_focal"])
            bound, = bound_paper_rows(row)
            self.assertEqual("matched" if row["proposed_stable_id"] == "Q00001" else "unresolved", bound["tier_gene_binding"])

    def test_fetched_sibling_records_also_disambiguate_external_id_without_extra_lookup(self):
        records = {"Q00001": record(), "Q00002": record("Q00002")}
        records["Q00002"]["uniProtKBCrossReferences"][0]["id"] = "99999"
        with mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_reference_candidates",
                        return_value=(["Q00001", "Q00002"], "ambiguous_gene_name", "fixture", 2, [], [])), \
             mock.patch("candidate_qtl.reference_gene_resolution._fetch_uniprot",
                        side_effect=lambda accession, timeout: records[accession]) as fetch:
            rows = self.resolve([paper("GeneID:12345")], records=[], allow_external_api=True)
        self.assertEqual(2, fetch.call_count)
        anchors = {row["proposed_stable_id"]: json.loads(row["paper_reports_json"])[0]["reference_gene_id"] for row in rows}
        self.assertEqual({"Q00001": "Q00001", "Q00002": "GeneID:12345"}, anchors)

    def test_shared_supplied_external_id_does_not_choose_a_provider_copy(self):
        original = paper("GeneID:12345")
        rows = self.resolve([original], records=[record(), record("Q00002")], stable="GeneID:12345")
        self.assertEqual(2, len(rows))
        for row in rows:
            self.assertEqual([original], json.loads(row["paper_reports_json"]))

    def test_missing_sibling_record_leaves_external_id_ownership_unconfirmed(self):
        original = paper("GeneID:12345")
        with mock.patch("candidate_qtl.reference_gene_resolution_legacy._resolve_reference_candidates",
                        return_value=(["Q00001", "Q00002"], "ambiguous_gene_name", "fixture", 2, [], [])), \
             mock.patch("candidate_qtl.reference_gene_resolution._fetch_uniprot",
                        side_effect=lambda accession, timeout: record() if accession == "Q00001" else None):
            rows = self.resolve([original], records=[], allow_external_api=True)
        self.assertEqual(2, len(rows))
        for row in rows:
            self.assertEqual([original], json.loads(row["paper_reports_json"]))

    def test_missing_or_wrong_species_provider_record_cannot_authorize_transport(self):
        for records in [[], [record(species="Species beta")]]:
            with self.subTest(records=records):
                original = paper()
                row, = self.resolve([original], records=records)
                self.assertEqual("", row["proposed_stable_id"])
                self.assertEqual([original], json.loads(row["paper_reports_json"]))
                self.assertNotIn("paper_reference_rebindings", json.loads(row["retrieval_provenance"]))

    def test_homolog_reported_species_does_not_replace_or_veto_reference_species(self):
        original = {**paper(species="Species beta"), "reported_gene": "BETA1", "gene_correspondence": "homolog"}
        row, = self.resolve([original])
        self.assertEqual("Species alpha", row["reference_species"])
        resolved, = json.loads(row["paper_reports_json"])
        self.assertEqual({**original, "reference_gene_id": "Q00001"}, resolved)
        bound, = bound_paper_rows(row)
        self.assertEqual("matched", bound["tier_gene_binding"])
        self.assertEqual("T2", classify_observation(bound)[0])

    def test_composite_parent_label_does_not_bind_to_a_sibling_atomic_reference(self):
        label = "LONG NODE (LN1); OTHER NODE (ON1)"
        original = paper(label)
        rows = self.resolve([original], records=[record(), record("Q00002", name="ON1")], name=label)
        self.assertEqual(2, len(rows))
        for row in rows:
            self.assertEqual([original], json.loads(row["paper_reports_json"]))


if __name__ == "__main__":
    unittest.main()
