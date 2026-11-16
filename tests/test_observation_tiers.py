"""Source observations drive Tier without a model-supplied label or rereview."""
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.paper_reports import PAPER_FIELD, normalize_paper_reports
from candidate_qtl.relevance_tiers import NEW_TIERS, annotate_relevance, classify_relevance
from candidate_qtl.observation_tiers import gene_binding
from tests import test_paper_tiers as fixtures
from tests import test_route_a_review_relevance as a_fixtures


def report(types=("phenotype",), **fields):
    value = fixtures.paper()
    del value["initial_tier"]
    del value["initial_tier_reason"]
    return {**value, "evidence_types": list(types), **fields}


def annotate(reports, **fields):
    return annotate_relevance([{"candidate_gene_id": "G1", "family_id": "F", "support_status": "unverified_route_b_review"}],
                              [fixtures.claim(reports, **fields)])[0]


class ObservationTierTests(unittest.TestCase):
    def test_a_native_json_observations_import_without_any_tier_label(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            observation = report(("phenotype", "function"), reported_function="meristem maintenance")
            payload = {**fixtures.MECHANISMS[0], "reference_gene_name": "REF1", "source_ids": "PMID:123",
                       PAPER_FIELD: [observation]}
            self.assertEqual(0, main(["mechanism-import", "--payload", str(fixtures._payload(root, [payload])),
                                      "--out-dir", str(root / "import")]))
            rows = pd.read_csv(root / "import/mechanism_hypotheses.tsv", sep="\t").fillna("")
            stored = json.loads(rows.iloc[0][PAPER_FIELD])[0]
            self.assertEqual(["function", "phenotype"], json.loads(stored["evidence_types"]))
            self.assertEqual("", stored["initial_tier"])
            self.assertEqual("meristem maintenance", stored["reported_function"])

    def test_rule_matrix_needs_no_declared_tier(self):
        examples = [(report(), NEW_TIERS[0]),
                    (report(phenotype_comparison="same_cellular_process"), NEW_TIERS[1]),
                    (report(phenotype_comparison="shared_process", reported_organ="root"), NEW_TIERS[2]),
                    (report(("function",), reported_function="meristem maintenance"), NEW_TIERS[3]),
                    (report(("annotation",), reported_annotation="specific receptor signaling domain"), NEW_TIERS[4]),
                    (report(()), NEW_TIERS[5]), (report(("unknown",)), NEW_TIERS[6])]
        for observation, expected in examples:
            with self.subTest(expected=expected):
                self.assertEqual(expected, classify_relevance(fixtures.claim([observation]))[0])

    def test_declared_tier_cannot_change_machine_classification_or_display(self):
        outputs = [annotate([report(initial_tier=tier, initial_tier_reason="ignored declaration")]) for tier in ("T1", "T3", "T6", "")]
        for row in outputs:
            self.assertEqual(NEW_TIERS[0], row["relevance_tier"])
            self.assertEqual("T1", json.loads(row["tier_evidence_json"])["initial_tier"])
            self.assertEqual("review", row["support_state"])
            self.assertEqual("observation_rules_v1", row["tier_classification_policy"])

    def test_other_species_same_trait_is_t2_not_t1(self):
        row = annotate([report(reported_species="Plant beta")], reference_species="Plant beta")
        self.assertEqual(NEW_TIERS[1], row["relevance_tier"])

    def test_unknown_comparison_retains_supported_function_with_partial_state(self):
        row = annotate([report(("phenotype", "function"), phenotype_comparison="unknown", reported_function="meristem maintenance")])
        self.assertEqual(NEW_TIERS[3], row["relevance_tier"])
        self.assertEqual("partial", row["tier_classification_state"])
        self.assertIn("phenotype_comparison", json.loads(row["tier_missing_fields_json"]))
        self.assertEqual("review", row["support_state"])
        self.assertEqual(1, len(json.loads(row["paper_evidence_summary_json"])))

    def test_unknown_species_does_not_become_a_negative_or_an_automatic_t1(self):
        row = annotate([report(reported_species="")])
        self.assertEqual(NEW_TIERS[1], row["relevance_tier"])
        self.assertEqual(["species_comparison"], json.loads(row["tier_missing_fields_json"]))

    def test_wrong_gene_or_species_cannot_borrow_the_papers_t1(self):
        for changes in ({"reported_gene": "OTHER_GENE"}, {"reported_species": "Plant beta"}):
            with self.subTest(changes=changes):
                observation = report(**changes)
                self.assertEqual(NEW_TIERS[6], classify_relevance(fixtures.claim([observation]))[0])
                row = annotate([observation])
                self.assertEqual(NEW_TIERS[5], row["relevance_tier"])
                self.assertIn("gene_correspondence", json.loads(row["tier_missing_fields_json"]))
                self.assertEqual("G1", row["candidate_gene_id"])

    def test_existing_endpoint_formats_and_names_are_used_without_lookups(self):
        self.assertEqual("matched", gene_binding({"reference_gene": "sp|P12345|NAME;sp|Q12345|OTHER"}, "UniProt:Q12345"))
        self.assertEqual("matched", gene_binding({"proposed_stable_id": "Q12345", "reference_gene_name": "REF1"}, "REF1"))
        self.assertEqual("matched", gene_binding({"reference_gene": "TAIR:AT1G12345"}, "AT1G12345"))
        self.assertEqual("unresolved", gene_binding({"reference_gene": "REF1"}, "REF2"))

    def test_function_is_not_invented_from_the_proposed_mechanism(self):
        row = annotate([report(("function",), reported_function="", initial_tier="T1")])
        self.assertEqual(NEW_TIERS[5], row["relevance_tier"])
        self.assertIn("reported_function", json.loads(row["tier_missing_fields_json"]))

    def test_annotation_without_paper_uses_original_annotation_and_no_tier_label(self):
        claim = {"candidate_gene_id": "G1", "reference_gene": "REF1", "source_id": "UniProt:P12345",
                 "evidence_types": ["annotation"], "reported_annotation": "specific meristem receptor domain",
                 "phenotype_relation": "could alter meristem size", PAPER_FIELD: "[]"}
        row = annotate_relevance([{"candidate_gene_id": "G1"}], [claim])[0]
        self.assertEqual(NEW_TIERS[4], row["relevance_tier"])
        self.assertEqual("classified", row["tier_classification_state"])

    def test_invalid_types_and_comparisons_are_not_silently_coerced(self):
        for observation in (report(evidence_types="phenotype"), report(evidence_types=[True]),
                            report(evidence_types=["inferred_phenotype"]), report(phenotype_comparison="probably_same")):
            with self.subTest(observation=observation), self.assertRaises(ValueError):
                normalize_paper_reports([observation])

    def test_missing_reported_gene_remains_partial_without_rejecting_the_candidate(self):
        row = annotate([report(reported_gene="")])
        self.assertEqual(NEW_TIERS[5], row["relevance_tier"])
        self.assertEqual("G1", row["candidate_gene_id"])
        self.assertEqual("partial", row["tier_classification_state"])
        self.assertIn("gene_correspondence", json.loads(row["tier_missing_fields_json"]))

    def test_complete_same_tier_report_is_displayed_without_losing_partial_report(self):
        row = annotate([report(reported_species=""),
                        report(source_id="PMID:456", phenotype_comparison="same_cellular_process")])
        self.assertEqual(NEW_TIERS[1], row["relevance_tier"])
        self.assertEqual("PMID:456", row["tier_source_id"])
        self.assertEqual("classified", row["tier_classification_state"])
        self.assertEqual(2, len(json.loads(row["paper_evidence_summary_json"])))
    def test_counts_and_replay_preserve_observations_and_missingness(self):
        first = report(("phenotype", "function"), phenotype_comparison="unknown", reported_function="meristem maintenance", publication_type="review")
        second = report(source_id="PMID:456", publication_type="original")
        row = annotate([first, second])
        replay = annotate_relevance([row])[0]
        for key in ("relevance_tier", "paper_tier_counts_json", "paper_evidence_summary_json", "tier_evidence_json", "tier_classification_state"):
            self.assertEqual(row[key], replay[key], key)
        self.assertEqual({"T1": 1, "T4": 1}, json.loads(row["paper_tier_counts_json"]))
        summaries = json.loads(row["paper_evidence_summary_json"])
        partial = next(item for item in summaries if item["tier"] == "T4")
        self.assertEqual("meristem maintenance", partial["reported_function"])
        self.assertIn("phenotype_comparison", json.loads(partial["tier_missing_fields_json"]))

    def test_route_b_without_tier_reaches_wet_and_counts_in_1000_gene_inventory(self):
        fixtures.PaperTierTests()._assert_b_model_response_reaches_wet(machine=True)

    def test_route_a_observations_reach_projection_and_wet_without_tier(self):
        with tempfile.TemporaryDirectory() as tmp:
            reports = [report(phenotype_comparison="same_cellular_process")]
            wet, _, _, _ = a_fixtures.RouteAReviewRelevanceTests()._run(Path(tmp),
                {"source_ids": "PMID:123", PAPER_FIELD: normalize_paper_reports(reports)})
            rows = pd.read_csv(wet["wet_candidates"], sep="\t").fillna("")
            self.assertEqual({"g1", "g2"}, set(rows.candidate_gene_id))
            self.assertEqual({NEW_TIERS[1]}, set(rows.relevance_tier))
            self.assertEqual({"observation_rules_v1"}, set(rows.tier_classification_policy))
