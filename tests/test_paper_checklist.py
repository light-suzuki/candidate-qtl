"""Original research answers a factual checklist; machines apply a tiny table."""
import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.observation_tiers import CHECKLIST_POLICY, checklist_tier, normalize_checklist
from candidate_qtl.cli import main
from candidate_qtl.paper_reports import PAPER_FIELD, normalize_paper_reports
from candidate_qtl.relevance_tiers import NEW_TIERS, annotate_relevance, classify_relevance
from tests import test_paper_tiers as fixtures
from tests import test_route_a_review_relevance as route_a


def report(**answers):
    return {"source_id": "PMID:123", "reported_gene": "NCBI Gene:12345",
            "reference_gene_id": "REF1", "gene_correspondence": "same_gene",
            "initial_tier_reason": "Source reports the target phenotype for this mapped gene.",
            "paper_checklist": {"identity_resolved": "yes", **answers}}


class PaperChecklistTests(unittest.TestCase):
    def test_native_a_response_without_a_tier_imports_the_checklist(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = {**fixtures.MECHANISMS[0], "reference_gene_name": "REF1", "source_ids": "PMID:123",
                PAPER_FIELD: [report(phenotype_reported="yes", equivalent_target_phenotype="yes",
                                     focal_species="yes", same_gene_as_focal="yes")]}
            self.assertEqual(0, main(["mechanism-import", "--payload", str(fixtures._payload(root, [payload])), "--out-dir", str(root / "import")]))
            imported = pd.read_csv(root / "import/mechanism_hypotheses.tsv", sep="\t").fillna("")
            stored = json.loads(imported.iloc[0][PAPER_FIELD])[0]
            self.assertEqual("T1", stored["initial_tier"])
            self.assertEqual("yes", json.loads(stored["paper_checklist"])["same_gene_as_focal"])

    def test_decision_table_does_not_read_names_species_or_reason_as_classification(self):
        cases = [
            ({"phenotype_reported": "yes", "equivalent_target_phenotype": "yes", "focal_species": "yes", "same_gene_as_focal": "yes"}, "T1"),
            ({"phenotype_reported": "yes", "equivalent_target_phenotype": "yes", "focal_species": "no"}, "T2"),
            ({"phenotype_reported": "yes", "similar_cellular_or_developmental_phenotype": "yes"}, "T2"),
            ({"phenotype_reported": "yes", "shared_specific_process": "yes"}, "T3"),
            ({"function_supports_hypothesis": "yes"}, "T4"),
            ({"annotation_supports_hypothesis": "yes"}, "T5"),
            ({"identity_resolved": "yes"}, "T6"), ({}, "U"),
        ]
        for checks, expected in cases:
            with self.subTest(expected=expected):
                self.assertEqual(expected, checklist_tier(checks))
                payload = report(**checks)
                payload["paper_checklist"] = checks
                payload["initial_tier"] = "T1"  # Inconsistent free Tier cannot override facts.
                normalized = json.loads(normalize_paper_reports([payload]))[0]
                self.assertEqual(expected, normalized["initial_tier"])
                expected_label = NEW_TIERS[int(expected[1:]) - 1] if expected != "U" else NEW_TIERS[-1]
                self.assertEqual(expected_label, classify_relevance(fixtures.claim([payload], candidate_gene_id=""))[0])

    def test_same_locus_answer_alone_never_creates_trait_evidence(self):
        # same_gene_as_focal is retained for downstream member/projection checks,
        # but it cannot fabricate the phenotype evidence needed for T1-T3.
        for same_gene in ("yes", "no", "unknown"):
            with self.subTest(same_gene=same_gene):
                self.assertEqual("U", checklist_tier({"same_gene_as_focal": same_gene}))
                self.assertEqual("T6", checklist_tier({"same_gene_as_focal": same_gene, "identity_resolved": "yes"}))
                self.assertEqual("T1", checklist_tier({"same_gene_as_focal": same_gene,
                    "phenotype_reported": "yes", "equivalent_target_phenotype": "yes", "focal_species": "yes"}))

    def test_unknown_is_preserved_and_not_treated_as_a_yes(self):
        checks = {"phenotype_reported": "unknown", "equivalent_target_phenotype": "yes", "identity_resolved": "yes"}
        normalized = json.loads(normalize_checklist(checks))
        self.assertEqual("unknown", normalized["focal_species"])
        self.assertEqual("T6", checklist_tier(checks))
        for value in ({"phenotype_reported": "maybe"}, {"free_question": "yes"}):
            with self.assertRaises(ValueError):
                normalize_checklist(value)

    def test_original_answers_survive_a_to_wet_without_a_tier_declaration(self):
        payload = report(phenotype_reported="yes", equivalent_target_phenotype="yes",
                         focal_species="yes", same_gene_as_focal="yes")
        payload["reference_gene_id"] = "bait1"
        with tempfile.TemporaryDirectory() as tmp:
            wet, _, _, _ = route_a.RouteAReviewRelevanceTests()._run(Path(tmp),
                {"source_ids": "PMID:123", PAPER_FIELD: normalize_paper_reports([payload])}, alignment_span=90)
            rows = pd.read_csv(wet["wet_candidates"], sep="\t").fillna("")
            self.assertEqual({NEW_TIERS[0]}, set(rows.relevance_tier))
            self.assertEqual({CHECKLIST_POLICY}, set(rows.tier_classification_policy))
            for summary in rows.paper_evidence_summary_json:
                checks = json.loads(json.loads(summary)[0]["paper_checklist"])
                self.assertEqual("yes", checks["same_gene_as_focal"])

    def test_unrelated_reference_anchor_cannot_supply_a_paper_tier(self):
        payload = report(phenotype_reported="yes", equivalent_target_phenotype="yes",
                         focal_species="yes", same_gene_as_focal="yes")
        for changes in ({"reference_gene_id": "OTHER"}, {"gene_correspondence": "unknown"}):
            self.assertEqual(NEW_TIERS[-1], classify_relevance(fixtures.claim([{**payload, **changes}]))[0])

    def test_reported_gene_and_proxy_are_separate_for_t1_and_t2(self):
        payload = {**report(phenotype_reported="yes", equivalent_target_phenotype="yes",
                           focal_species="no", same_gene_as_focal="no"),
                   "reported_gene": "Brassica locus X", "reported_gene_id": "GenBank:AB123456",
                   "reported_species": "Brassica rapa", "gene_correspondence": "homolog",
                   "initial_tier_reason": "Brassica paper reports a similar phenotype; supplied At reference is the documented sequence proxy."}
        row = annotate_relevance([{"candidate_gene_id": "G1"}], [fixtures.claim([payload])])[0]
        self.assertEqual(NEW_TIERS[1], row["relevance_tier"])
        summary = json.loads(row["paper_evidence_summary_json"])[0]
        self.assertEqual("Brassica locus X", summary["reported_gene"])
        self.assertEqual("GenBank:AB123456", summary["reported_gene_id"])
        # Paper-side T1 remains conditional on the eventual member projection,
        # even when a different-species homolog supplies the sequence proxy.
        payload["paper_checklist"].update(focal_species="yes", same_gene_as_focal="yes")
        payload.update(reported_species="Pisum sativum", reported_gene="old focal locus",
                       initial_tier_reason="Old and new focal loci correspond; At is only the sequence proxy.")
        self.assertEqual(NEW_TIERS[0], classify_relevance(fixtures.claim([payload], candidate_gene_id=""))[0])

    def test_strongest_paper_uses_checklist_and_replays_without_count_bonus(self):
        strong = report(phenotype_reported="yes", similar_cellular_or_developmental_phenotype="yes")
        weak = {**report(function_supports_hypothesis="yes"), "source_id": "PMID:456"}
        row = annotate_relevance([{"candidate_gene_id": "G1"}], [fixtures.claim([strong, weak])])[0]
        self.assertEqual(NEW_TIERS[1], row["relevance_tier"])
        self.assertEqual({"T2": 1, "T4": 1}, json.loads(row["paper_tier_counts_json"]))
        self.assertEqual(row["relevance_tier"], annotate_relevance([row])[0]["relevance_tier"])

    def test_winning_checklist_metadata_does_not_inherit_an_unresolved_paper(self):
        cases = [
            dict(phenotype_reported="yes", equivalent_target_phenotype="yes",
                 focal_species="yes", same_gene_as_focal="yes"),
            dict(phenotype_reported="yes", similar_cellular_or_developmental_phenotype="yes"),
            dict(phenotype_reported="yes", shared_specific_process="yes"),
        ]
        unresolved = {**report(function_supports_hypothesis="yes"),
                      "source_id": "PMID:456", "gene_correspondence": "unknown"}
        for index, answers in enumerate(cases):
            with self.subTest(tier=index + 1):
                strong = report(**answers)
                for papers in ([strong, unresolved], [unresolved, strong]):
                    projection = json.dumps({"basis": "frozen_gene_orthology", "candidate_gene_id": "G1",
                        "reference_gene_id": "REF1", "member_genes": ["G1"],
                        "evidence": {"fixture": "frozen REF1 to G1 orthology"}})
                    row = annotate_relevance([{"candidate_gene_id": "G1"}],
                        [fixtures.claim(papers, homolog_projection_json=projection)])[0]
                    self.assertEqual(NEW_TIERS[index], row["relevance_tier"])
                    for current in (row, annotate_relevance([row])[0]):
                        self.assertEqual(CHECKLIST_POLICY, current["tier_classification_policy"])
                        self.assertEqual("classified", current["tier_classification_state"])
                        self.assertEqual("matched", current["tier_gene_binding"])
                        self.assertEqual([], json.loads(current["tier_missing_fields_json"]))
                        selected = json.loads(current["tier_evidence_json"])
                        self.assertEqual(current["tier_classification_state"], selected["tier_classification_state"])
                        self.assertEqual(current["tier_gene_binding"], selected["tier_gene_binding"])
                        evidence = json.loads(current["relevance_evidence_json"])
                        reports = [paper for item in evidence for paper in json.loads(item.get(PAPER_FIELD, "[]"))]
                        self.assertTrue(any(paper["source_id"] == "PMID:456"
                                            and paper["gene_correspondence"] == "unknown" for paper in reports))

    def test_winning_checklist_metadata_reaches_wet_cli_with_unresolved_paper_preserved(self):
        strong = report(phenotype_reported="yes", equivalent_target_phenotype="yes",
                        focal_species="yes", same_gene_as_focal="yes")
        strong["reference_gene_id"] = "bait1"
        unresolved = {**report(function_supports_hypothesis="yes"),
                      "source_id": "PMID:456", "gene_correspondence": "unknown"}
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            wet, bridge, _, _ = route_a.RouteAReviewRelevanceTests()._run(root,
                {"source_ids": "PMID:123;PMID:456", PAPER_FIELD: normalize_paper_reports([strong, unresolved])}, alignment_span=90)
            output = root / "wet-cli"
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(0, main([
                    "wet-candidate-handoff", "--out-dir", str(output), "--case-id", "route-a-review",
                    "--wet-budget", "2", "--route-a-bridge", bridge["bridge_tsv"],
                    "--route-a-manifest", bridge["manifest"],
                    "--review-candidates", bridge["review_candidates"],
                    "--review-candidates-manifest", bridge["review_manifest"],
                    "--interval-atlas", str(root / "genes.tsv"),
                ]))
            original = pd.read_csv(wet["wet_candidates"], sep="\t", dtype=str).fillna("")
            rows = pd.read_csv(output / "wet_candidates.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(set(original.candidate_gene_id), set(rows.candidate_gene_id))
            self.assertEqual({NEW_TIERS[0]}, set(rows.relevance_tier))
            self.assertEqual({CHECKLIST_POLICY}, set(rows.tier_classification_policy))
            self.assertEqual({"classified"}, set(rows.tier_classification_state))
            self.assertEqual({"matched"}, set(rows.tier_gene_binding))
            self.assertEqual({"[]"}, set(rows.tier_missing_fields_json))
            for summary in rows.paper_evidence_summary_json:
                papers = {paper["source_id"]: paper for paper in json.loads(summary)}
                self.assertEqual({"PMID:123", "PMID:456"}, set(papers))
                self.assertEqual("partial", papers["PMID:456"]["tier_classification_state"])
                self.assertEqual(["gene_correspondence"], json.loads(papers["PMID:456"]["tier_missing_fields_json"]))

    def test_b_native_response_to_1000_gene_wet_uses_the_same_checklist(self):
        fixtures.PaperTierTests()._assert_b_model_response_reaches_wet(checklist=True)
