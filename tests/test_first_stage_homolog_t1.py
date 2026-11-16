"""#1002: initial T1 is a same-species reported-phenotype homolog candidate."""
import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from candidate_qtl.observation_tiers import CHECKLIST_FIELDS, CHECKLIST_POLICY, checklist_tier, classify_observation
from candidate_qtl.paper_reports import PAPER_FIELD, normalize_paper_reports
from candidate_qtl.relevance_tiers import NEW_TIERS
from candidate_qtl.route_b_semantic_candidates import route_b_homolog_projection
from tests import test_route_a_review_relevance as a_fixture
from tests import test_route_b_initial_reference_identity as b_fixture


def report(**changes):
    checks = dict.fromkeys(CHECKLIST_FIELDS, "unknown")
    checks.update(focal_species="yes", phenotype_reported="yes", equivalent_target_phenotype="yes")
    checks.update(changes.pop("checks", {}))
    return {"source_id": "PMID:1234567", "source_title": "Synthetic meristem observation",
            "reported_gene": "reported focal gene", "reported_species": "Solanum lycopersicum",
            "reference_gene_id": "REF1", "gene_correspondence": "homolog",
            "reported_phenotype": "altered meristem development",
            "initial_tier_reason": "Focal gene mutation alters meristem development; the supplied Arabidopsis REF1 is its sequence homolog proxy.",
            "paper_checklist": checks, **changes}


class FirstStageHomologT1Tests(unittest.TestCase):
    def test_same_locus_is_not_required_and_close_specific_phenotype_qualifies(self):
        for identity in ("unknown", "no", "yes"):
            for match in ("equivalent_target_phenotype", "similar_cellular_or_developmental_phenotype"):
                checks = {"focal_species": "yes", "phenotype_reported": "yes", match: "yes", "same_gene_as_focal": identity}
                self.assertEqual("T1", checklist_tier(checks))
                self.assertEqual("T2", checklist_tier({**checks, "focal_species": "no"}))

    def test_a_nonidentical_two_copy_projection_reaches_wet_without_identity_review(self):
        paper = report(reference_gene_id="bait1", reported_gene="old assembly gene", reported_species="Plant alpha")
        with tempfile.TemporaryDirectory() as tmp:
            wet, _, original, _ = a_fixture.RouteAReviewRelevanceTests()._run(Path(tmp),
                {"source_ids": "PMID:1234567", "reference_species": "Arabidopsis thaliana",
                 PAPER_FIELD: normalize_paper_reports([paper])}, alignment_span=90, review_only=True)
            rows = pd.read_csv(wet["wet_candidates"], sep="\t", dtype=str).fillna("")
            self.assertEqual({"g1", "g2"}, set(rows.candidate_gene_id))
            self.assertEqual({NEW_TIERS[0]}, set(rows.relevance_tier))
            self.assertEqual({"review"}, set(rows.support_state))
            for encoded in rows.tier_evidence_json:
                evidence = json.loads(encoded)
                self.assertEqual("unknown", json.loads(evidence["paper_checklist"])["same_gene_as_focal"])
                self.assertEqual("Arabidopsis thaliana", evidence["proxy_reference_species"])
                projection = json.loads(evidence["homolog_projection_json"])
                self.assertEqual(2, projection["copy_count"])
                self.assertEqual({"g1", "g2"}, set(projection["member_genes"]))
            self.assertIn("unknown", original.read_text())

    def test_a_short_local_projection_stays_t2(self):
        for identity in ("unknown", "no", "yes"):
            paper = report(reference_gene_id="bait1", checks={"same_gene_as_focal": identity})
            with self.subTest(identity=identity), tempfile.TemporaryDirectory() as tmp:
                wet, _, _, _ = a_fixture.RouteAReviewRelevanceTests()._run(Path(tmp),
                    {"source_ids": "PMID:1234567", PAPER_FIELD: normalize_paper_reports([paper])}, alignment_span=40)
                rows = pd.read_csv(wet["wet_candidates"], sep="\t", dtype=str)
                self.assertEqual({"g1", "g2"}, set(rows.candidate_gene_id))
                self.assertEqual({NEW_TIERS[1]}, set(rows.relevance_tier))
                for encoded in rows.tier_evidence_json:
                    evidence = json.loads(encoded)
                    self.assertEqual(paper["paper_checklist"], json.loads(evidence["paper_checklist"]))

    def test_a_verified_reference_still_retains_all_nonidentical_focal_report_copies(self):
        for identity in ("unknown", "no", "yes"):
            paper = report(reference_gene_id="bait1", checks={"same_gene_as_focal": identity})
            with self.subTest(identity=identity), tempfile.TemporaryDirectory() as tmp:
                wet, _, _, _ = a_fixture.RouteAReviewRelevanceTests()._run(Path(tmp),
                    {"source_ids": "PMID:1234567", PAPER_FIELD: normalize_paper_reports([paper])},
                    alignment_span=90)
                rows = pd.read_csv(wet["wet_candidates"], sep="\t", dtype=str)
                self.assertEqual({"g1", "g2"}, set(rows.candidate_gene_id))
                self.assertEqual({NEW_TIERS[0]}, set(rows.relevance_tier))
                for encoded in rows.tier_evidence_json:
                    evidence = json.loads(encoded)
                    self.assertEqual(paper["paper_checklist"], json.loads(evidence["paper_checklist"]))
                    projection = json.loads(evidence["homolog_projection_json"])
                    self.assertEqual(2, projection["copy_count"])
                    self.assertEqual({"g1", "g2"}, set(projection["member_genes"]))

    def test_paper_side_t1_remains_conditional_and_requires_a_reason_for_all_identity_answers(self):
        for identity in ("unknown", "no", "yes"):
            paper = report(checks={"same_gene_as_focal": identity})
            row = {**paper, "tier_classification_policy": CHECKLIST_POLICY, "tier_gene_binding": "matched"}
            with self.subTest(identity=identity):
                self.assertEqual("T1", classify_observation(row)[0])
                self.assertEqual("T2", classify_observation({**row, "initial_tier_reason": ""})[0])

    def test_b_sequence_endpoint_span_and_species_are_not_borrowed_from_other_hits(self):
        claim = {"reference_gene": "REF1", "reference_species": "Arabidopsis thaliana"}
        hit = {"hit_accession": "REF1", "hit_species": "Arabidopsis thaliana", "percent_identity": "75", "query_coverage": "0.8"}
        def projection(hits, gene="G1"):
            return route_b_homolog_projection(claim, {"unit_id": "unit", "genes": [
                {"gene_id": "G1", "sequence_similarity_hits": hits}]}, gene, {"G1", "outside"})
        self.assertEqual("qualified_sequence_hit", json.loads(projection([hit]))["basis"])
        for hits in ([{**hit, "hit_accession": "OTHER"}], [{**hit, "query_coverage": "0.1"}],
                     [{**hit, "query_coverage": ""}], [{**hit, "query_coverage": "nan"}],
                     [{**hit, "hit_species": "Medicago truncatula"}]):
            self.assertEqual("", projection(hits))
        self.assertEqual("", projection([hit], "outside"))
        self.assertEqual("", projection([{**hit, "query_coverage": "0.1"},
                                         {**hit, "hit_species": "Medicago truncatula"}]))

    def test_b_reasoned_foreign_proxy_keeps_one_three_and_26_copies(self):
        for copies in (1, 3, 26):
            with self.subTest(copies=copies), tempfile.TemporaryDirectory() as tmp:
                _, _, frozen, wet = b_fixture.RouteBInitialReferenceIdentityTests()._chain(
                    Path(tmp), [b_fixture.ARABIDOPSIS], paper_reports=[report()], copies=copies)
                trace = pd.read_csv(Path(wet["wet_candidates"]).with_name("wet_candidate_trace.tsv"), sep="\t", dtype=str)
                self.assertEqual({"retained"}, set(trace.selection_disposition))
                for path in (frozen["review_candidates"], wet["wet_candidates"]):
                    rows = pd.read_csv(path, sep="\t", dtype=str).fillna("")
                    self.assertEqual(copies, rows.candidate_gene_id.nunique())
                    self.assertEqual({NEW_TIERS[0]}, set(rows.relevance_tier))
                    self.assertEqual({"review"}, set(rows.support_state))
                    for encoded in rows.paper_evidence_summary_json:
                        item = json.loads(encoded)[0]
                        self.assertEqual("Solanum lycopersicum", item["reported_species"])
                        self.assertEqual("Arabidopsis thaliana", item["proxy_reference_species"])
                        self.assertEqual("unknown", json.loads(item["paper_checklist"])["same_gene_as_focal"])

    def test_b_wrong_species_distant_phenotype_unbound_proxy_and_domain_only_do_not_promote(self):
        cases = [({"checks": {"focal_species": "no"}, "reported_species": "Arabidopsis thaliana"}, True),
                 ({"checks": {"equivalent_target_phenotype": "no", "similar_cellular_or_developmental_phenotype": "no", "shared_specific_process": "yes"}}, True),
                 ({"gene_correspondence": "unknown"}, True), ({}, False)]
        for change, endpoint in cases:
            with self.subTest(change=change, endpoint=endpoint), tempfile.TemporaryDirectory() as tmp:
                _, _, frozen, _ = b_fixture.RouteBInitialReferenceIdentityTests()._chain(
                    Path(tmp), [b_fixture.ARABIDOPSIS], paper_reports=[report(**change)], gene_endpoint=endpoint)
                rows = pd.read_csv(frozen["review_candidates"], sep="\t", dtype=str).fillna("")
                self.assertNotIn(NEW_TIERS[0], set(rows.relevance_tier))


if __name__ == "__main__":
    unittest.main()
