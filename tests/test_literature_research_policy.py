"""Ordinary research breadth and representation-independent paper Tier wiring."""
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl import phenotype_investigation as pi
from candidate_qtl.cli import main
from candidate_qtl.observation_tiers import checklist_tier, normalize_checklist
from candidate_qtl.paper_reports import PAPER_FIELD, normalize_paper_reports, select_bibliography_papers
from candidate_qtl.relevance_tiers import NEW_TIERS, annotate_relevance
from candidate_qtl.route_b_semantic_import import import_route_b_semantic_proposals
from tests import test_paper_tiers as papers
from tests import test_route_b_semantic_import as route_b
from tests import test_route_a_review_relevance as route_a
from tests.test_paper_checklist import report


class LiteratureResearchPolicyTests(unittest.TestCase):
    def test_real_a_and_b_prompts_export_twenty_candidates_and_five_reports(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            a = pi.phenotype_investigate(str(root / "a"), "case", "Plant alpha", "altered meristem")
            a_policy = json.loads(Path(a["manifest"]).read_text())["payload_contract"]["literature_research"]
            (root / "b").mkdir()
            b = route_b.RouteBSemanticImportTests()._bundle(root / "b")
            prompt = json.loads(Path(b["prompts"]).read_text().splitlines()[0])
            self.assertEqual(a_policy, prompt["literature_research"])
            self.assertEqual((20, 5), (a_policy["candidate_paper_target"], a_policy["evidence_paper_target"]))
            for text in (Path(a["prompt"]).read_text(), prompt["system"]):
                for instruction in ("screen about 20", "five", "Originals and reviews", "next",
                                    "100 annotation-only", "shortfall", "not 20 full-text"):
                    self.assertIn(instruction, text)

    def test_one_t2_review_precedes_one_hundred_t5_reports_without_count_bonus(self):
        weak = [dict(report(annotation_supports_hypothesis="yes"), source_id=f"PMID:{i}") for i in range(101, 201)]
        strong = dict(report(phenotype_reported="yes", similar_cellular_or_developmental_phenotype="yes"),
                      source_id="PMID:999", publication_type="review")
        records = [*weak, strong]
        selection = select_bibliography_papers({PAPER_FIELD: records}, [r["source_id"] for r in records])
        self.assertEqual(101, selection["candidate_paper_count"])
        self.assertEqual(5, selection["selected_paper_count"])
        self.assertEqual("PMID:999", selection["selected_source_ids"][0])
        result = annotate_relevance([{"candidate_gene_id": "G1"}], [papers.claim(records)])[0]
        self.assertEqual(NEW_TIERS[1], result["relevance_tier"])
        only_weak = annotate_relevance([{"candidate_gene_id": "G1"}], [papers.claim(weak)])[0]
        self.assertEqual(NEW_TIERS[4], only_weak["relevance_tier"])

    def test_boolean_case_and_alias_array_reach_t2_through_real_a_consumer(self):
        paper = dict(report(phenotype_reported=True, similar_cellular_or_developmental_phenotype=" YES ",
                            focal_species=False), source_identifier_aliases=["DOI:10.1234/one"])
        normalized = json.loads(normalize_paper_reports([paper]))[0]
        self.assertEqual("T2", normalized["initial_tier"])
        self.assertEqual("DOI:10.1234/one", normalized["source_identifier_aliases"])
        self.assertEqual("no", json.loads(normalized["paper_checklist"])["focal_species"])
        with tempfile.TemporaryDirectory() as tmp:
            wet, _, _, _ = route_a.RouteAReviewRelevanceTests()._run(Path(tmp),
                {"source_ids": "PMID:123", PAPER_FIELD: normalize_paper_reports([paper])})
            rows = pd.read_csv(wet["wet_candidates"], sep="\t").fillna("")
            self.assertTrue(len(rows))
            self.assertEqual({NEW_TIERS[1]}, set(rows.relevance_tier))

    def test_real_a_import_preserves_top_level_alias_array(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = {**papers.MECHANISMS[0], "source_identifier_aliases": ["DOI:10.1234/one", "PMID:123"]}
            self.assertEqual(0, main(["mechanism-import", "--payload", str(papers._payload(root, [payload])),
                                      "--out-dir", str(root / "import")]))
            rows = pd.read_csv(root / "import/mechanism_hypotheses.tsv", sep="\t").fillna("")
            self.assertEqual("DOI:10.1234/one;PMID:123", rows.iloc[0].source_identifier_aliases)

    def test_real_b_import_accepts_same_facts_without_manual_type_repair(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = route_b.RouteBSemanticImportTests()
            built = fixture._bundle(root)
            paper = dict(report(phenotype_reported=True, equivalent_target_phenotype=True, focal_species=False),
                         source_identifier_aliases=["DOI:10.1234/one"])
            claim = {"claim_type": "mutant_phenotype", "biological_function": "meristem maintenance",
                     "mechanism": "receptor signaling", "phenotype_relation": "could explain enlarged meristems",
                     "reference_gene": "REF1", "reference_species": "Plant beta", "evidence_scope": "reference_gene",
                     "source_ids": ["PMID:123"], PAPER_FIELD: [paper]}
            response = root / "response.json"
            fixture._response(Path(built["prompts"]), [claim], response)
            result = import_route_b_semantic_proposals(built["bundle"], built["prompts"], str(response), str(root / "import"))
            self.assertEqual(1, result["proposed_claim_row_count"])
            row = pd.read_csv(result["proposals"], sep="\t").fillna("").iloc[0]
            self.assertEqual("T2", row.initial_tier)
            self.assertEqual("DOI:10.1234/one", row.source_identifier_aliases)

    def test_missing_facts_and_unknown_question_never_become_positive_answers(self):
        checks = json.loads(normalize_checklist({"identity_resolved": True, "phenotype_reported": " UNKNOWN "}))
        self.assertEqual("unknown", checks["equivalent_target_phenotype"])
        self.assertEqual("T6", checklist_tier(checks))
        for value in ({"phenotype_reported": 1}, {"unexpected": True}):
            with self.assertRaises(ValueError):
                normalize_checklist(value)
        with self.assertRaises(ValueError):
            normalize_paper_reports([dict(report(), source_identifier_aliases=[123])])


if __name__ == "__main__":
    unittest.main()
