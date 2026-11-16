"""Original research judgments survive identifier namespaces and A/B transport."""
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.observation_tiers import JUDGMENT_POLICY, gene_binding
from candidate_qtl.paper_reports import PAPER_FIELD, imported_paper_fields, normalize_paper_reports
from candidate_qtl.relevance_tiers import INITIAL_TIER_PROMPT, NEW_TIERS, annotate_relevance, classify_relevance
from tests import test_paper_tiers as fixtures
from tests import test_route_a_review_relevance as a_fixtures


def judgment(**changes):
    return {**fixtures.paper(tier="T2", reported_gene="NCBI Gene:111111", reported_species="source common name",
                          reference_gene_id="REF1", gene_correspondence="same_gene",
                          evidence_types=["phenotype"], initial_tier_reason="Source gene corresponds to the supplied reference; similar cellular phenotype.",
                          ), **changes}


class PaperJudgmentTests(unittest.TestCase):
    def test_original_judgment_does_not_require_a_name_or_species_dictionary(self):
        for gene in ("NCBI Gene:111111", "TAIR:AT1G12345", "old-assembly.model.7", "a source gene name"):
            report = judgment(reported_gene=gene)
            record = fixtures.claim([report])
            self.assertEqual(NEW_TIERS[1], classify_relevance(record)[0])
            row = annotate_relevance([{"candidate_gene_id": "G1", "support_status": "unverified_route_b_review"}], [record])[0]
            self.assertEqual(NEW_TIERS[1], row["relevance_tier"])
            self.assertEqual(JUDGMENT_POLICY, row["tier_classification_policy"])
            self.assertEqual("review", row["support_state"])
            saved = json.loads(row["paper_evidence_summary_json"])[0]
            self.assertEqual(gene, saved["reported_gene"])
            self.assertEqual("REF1", saved["reference_gene_id"])

    def test_existing_reference_anchor_is_required_even_for_a_model_tier(self):
        for changes in ({"reference_gene_id": "OTHER"}, {"gene_correspondence": "different_gene"},
                        {"gene_correspondence": "unknown"}):
            with self.subTest(changes=changes):
                record = fixtures.claim([judgment(**changes)])
                self.assertEqual(NEW_TIERS[6], classify_relevance(record)[0])
                row = annotate_relevance([{"candidate_gene_id": "G1"}], [record])[0]
                self.assertEqual("partial", row["tier_classification_state"])
                self.assertIn("gene_correspondence", json.loads(row["tier_missing_fields_json"]))

    def test_homolog_relation_can_carry_t2_but_cannot_authorize_t1(self):
        self.assertEqual(NEW_TIERS[1], classify_relevance(fixtures.claim([judgment(gene_correspondence="homolog")]))[0])
        self.assertEqual(NEW_TIERS[6], classify_relevance(fixtures.claim([judgment(gene_correspondence="homolog", initial_tier="T1")]))[0])

    def test_judgment_requires_a_tier_reason_and_explicit_relation(self):
        for changes in ({"initial_tier": ""}, {"initial_tier_reason": ""}, {"reference_gene_id": ""},
                        {"gene_correspondence": "probably"}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                normalize_paper_reports([judgment(**changes)])

    def test_strongest_paper_and_replay_keep_semantic_judgment_and_no_count_bonus(self):
        strong = judgment(initial_tier="T2")
        weak = judgment(source_id="PMID:456", initial_tier="T4", evidence_types=["function"], reported_function="meristem maintenance")
        row = annotate_relevance([{"candidate_gene_id": "G1"}], [fixtures.claim([strong, weak])])[0]
        replay = annotate_relevance([row])[0]
        self.assertEqual(NEW_TIERS[1], row["relevance_tier"])
        self.assertEqual({"T2": 1, "T4": 1}, json.loads(row["paper_tier_counts_json"]))
        for field in ("relevance_tier", "paper_tier_counts_json", "paper_evidence_summary_json", "tier_classification_policy", "tier_evidence_json"):
            self.assertEqual(row[field], replay[field], field)

    def test_empty_paper_list_preserves_database_function_without_inheriting_phenotype_or_tier(self):
        record = {"source_ids": ["UniProt:P12345"], PAPER_FIELD: [], "evidence_types": ["annotation", "function"],
                  "reported_annotation": "specific receptor signaling domain", "reported_function": "not a paper-bound function",
                  "reported_phenotype": "must not leak", "initial_tier": "T1"}
        result = imported_paper_fields(record, "UniProt:P12345")
        self.assertEqual(["annotation", "function"], json.loads(result["evidence_types"]))
        self.assertEqual(record["reported_annotation"], result["reported_annotation"])
        self.assertEqual(record["reported_function"], result["reported_function"])
        for key in ("reported_phenotype", "initial_tier"):
            self.assertEqual("", result[key])
        self.assertEqual("", imported_paper_fields(record, "UniProt:Q12345")["evidence_types"])

    def test_compatibility_formats_preserve_machine_t2_without_false_gene_matches(self):
        report = fixtures.paper(initial_tier="", initial_tier_reason="", reported_gene="SOURCE (UniProt:P12345)",
                                reported_species="Plant beta", evidence_types=["phenotype"], phenotype_comparison="same_equivalent")
        record = fixtures.claim([report], reference_gene="P12345", reference_species="Plant beta (common name)")
        self.assertEqual(NEW_TIERS[1], classify_relevance(record)[0])
        self.assertEqual("unresolved", gene_binding({"reference_gene": "P12345", "reference_gene_name": "SOURCE"}, "SOURCE (Q12345)"))
        self.assertEqual("unresolved", gene_binding({"reference_gene": "P12345", "reference_species": "Plant beta"}, "P12345", "Plant gamma"))

    def test_native_a_import_preserves_original_gene_id_and_original_judgment(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = {**fixtures.MECHANISMS[0], "reference_gene_name": "REF1", "source_ids": "PMID:123", PAPER_FIELD: [judgment()]}
            self.assertEqual(0, main(["mechanism-import", "--payload", str(fixtures._payload(root, [payload])), "--out-dir", str(root / "import")]))
            imported = pd.read_csv(root / "import/mechanism_hypotheses.tsv", sep="\t").fillna("")
            report = json.loads(imported.iloc[0][PAPER_FIELD])[0]
            self.assertEqual("NCBI Gene:111111", report["reported_gene"])
            self.assertEqual("REF1", report["reference_gene_id"])
            self.assertEqual("T2", report["initial_tier"])

    def test_a_projection_and_wet_preserve_judgment_and_both_copies(self):
        with tempfile.TemporaryDirectory() as tmp:
            wet, _, _, _ = a_fixtures.RouteAReviewRelevanceTests()._run(Path(tmp),
                {"source_ids": "PMID:123", PAPER_FIELD: normalize_paper_reports([judgment()])})
            rows = pd.read_csv(wet["wet_candidates"], sep="\t").fillna("")
            self.assertEqual({"g1", "g2"}, set(rows.candidate_gene_id))
            self.assertEqual({NEW_TIERS[1]}, set(rows.relevance_tier))
            self.assertEqual({JUDGMENT_POLICY}, set(rows.tier_classification_policy))

    def test_b_import_source_transfer_and_1000_gene_wet_preserve_judgment(self):
        fixtures.PaperTierTests()._assert_b_model_response_reaches_wet(semantic=True)

    def test_prompt_assigns_semantics_in_original_research_without_extra_stage(self):
        for term in ("reference_gene_id", "gene_correspondence", "NCBI Gene", "T2:", "No Tier-only search"):
            self.assertIn(term, INITIAL_TIER_PROMPT)


if __name__ == "__main__":
    unittest.main()
