"""Exact shared function judgments with separate copy application, not a benchmark."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from candidate_qtl.agent_claim_review import AgentClaimReview, B_CONTRACT, FUNCTION_AXES, FUNCTION_PREDICATE, prepare_review_prompt
from candidate_qtl.core import sha256
from candidate_qtl.route_b_semantic_transfer_verify import verify_route_b_semantic_transfer
from candidate_qtl.route_b_semantic_candidates import freeze_route_b_semantic_candidates
from candidate_qtl.semantic_claims import CLAIM_COLUMNS
from tests.test_agent_claim_review import ReviewFixture, write_json, spec
from tests.test_route_b_semantic_transfer_verify import _claim, _bundle


class SharedFunctionReviewTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.fixture = ReviewFixture(self.root)
        self.rows = [_claim(claim_id=f"c{i}", semantic_unit_id=f"unit-{i}", member_interval_genes=f"G{i}",
                            reference_gene="P12345", source_id="PMID:12345", claim_type="mechanistic_review_seed",
                            biological_function="reference root growth", mechanism="conditional cell elongation in focal organs",
                            phenotype_relation="possible shorter focal stem organs") for i in range(6)]
        self.claims = self.root / "claims.tsv"
        self.prompt = self.root / "prompt-b.json"
        sources = json.loads(self.fixture.sources_path.read_text())
        sources["records"]["PMID:12345"]["hypothesis_ids"] = [row["claim_id"] for row in self.rows]
        write_json(self.fixture.sources_path, sources)
        self.prepare()

    def prepare(self):
        pd.DataFrame(self.rows, columns=CLAIM_COLUMNS).to_csv(self.claims, sep="\t", index=False)
        prepare_review_prompt(str(self.claims), str(self.fixture.sources_path), str(self.prompt), route="b")
        return json.loads(self.prompt.read_text())

    def response(self, omit=()):
        prompt = json.loads(self.prompt.read_text())
        base = copy.deepcopy(self.fixture.review)
        base.pop("hypothesis_id")
        base.update(predicate_relation=FUNCTION_PREDICATE, reference_function="short roots",
                    focal_bridge="conditional cell elongation in focal organs; focal causality untested",
                    axes=dict.fromkeys(FUNCTION_AXES, True))
        questions = prompt["reference_function_questions"]["questions"]
        shared = [{**base, "question_id": q["question_id"]} for q in questions]
        application = [{"claim_id": cid, "reference_function_question_id": q["question_id"],
                        "copy_application_state": "uncertain",
                        "copy_application_reason": "Apply the shared conditional hypothesis to this frozen endpoint; copy remains ambiguous."}
                       for q in questions for cid in q["claim_ids"] if cid not in omit]
        return {"shared_reference_function_reviews": shared, "reviews": application}

    def packet(self, response):
        raw = write_json(self.root / "review-b.json", response)
        return write_json(self.root / "packet-b.json", {
            "contract": B_CONTRACT, "hypotheses_sha256": sha256(self.claims),
            "source_review_records_sha256": sha256(self.fixture.sources_path),
            "proposal_run": {"executor_id": "proposer", "model": "fixture", "raw_response": spec(self.fixture.proposal_raw)},
            "review_run": {"executor_id": "independent-reviewer", "model": "fixture", "raw_response": spec(raw), "prompt": spec(self.prompt)},
        })

    def run_transfer(self, response):
        bundle = self.root / "bundle.json"
        _bundle(bundle, members=[("G0", "many-to-many")])
        payload = json.loads(bundle.read_text())
        template = payload["units"][0]
        template["genes"][0].update(orthology_reference_gene="P12345", reference_species="Arabidopsis thaliana")
        template["literature_link_candidates"][0].update(reference_gene="P12345", source_id="PMID:12345")
        payload["units"] = []
        for i, row in enumerate(self.rows):
            unit = copy.deepcopy(template)
            unit.update(unit_id=row["semantic_unit_id"], member_interval_genes=[f"G{i}"])
            unit["genes"][0]["gene_id"] = f"G{i}"
            payload["units"].append(unit)
        payload["b_review_membership_policy"] = "source_function_v1"
        write_json(bundle, payload)
        source = self.root / "source.tsv"
        pd.DataFrame([{"claim_id": row["claim_id"], "source_relation_status": "review_required"} for row in self.rows]).to_csv(source, sep="\t", index=False)
        result = verify_route_b_semantic_transfer(str(self.claims), str(source), str(bundle), str(self.root / "transfer"),
                    agent_claim_review_arg=str(self.packet(response)), source_review_records_arg=str(self.fixture.sources_path),
                    record_cache_arg=str(self.fixture.records_path), case_id="case403")
        audit = json.loads((self.root / "transfer/agent_claim_review_audit.json").read_text())
        candidates = freeze_route_b_semantic_candidates(result["claims"], str(self.root / "candidates"), case_id="case403")
        return audit, pd.read_csv(candidates["review_candidates"], sep="\t", dtype=str).fillna("")

    def test_six_copy_applications_share_one_function_judgment_through_real_consumer(self):
        response = self.response()
        self.assertEqual(1, len(response["shared_reference_function_reviews"]))
        audit, candidates = self.run_transfer(response)
        self.assertEqual(1, audit["shared_reference_function_judgment_count"])
        self.assertEqual(6, audit["copy_application_judgment_count"])
        self.assertEqual({f"G{i}" for i in range(6)}, set(candidates.candidate_gene_id))
        self.assertEqual(1, candidates.reference_function_question_id.nunique())
        self.assertIsNone(audit["actual_source_read_count"])
        self.assertIsNone(audit["actual_tokens"])

    def test_function_species_source_or_copy_difference_separates_questions(self):
        original = copy.deepcopy(self.rows)
        for field, value in (("reference_gene", "other"), ("reference_species", "different species"),
                             ("biological_function", "distinct function"), ("mechanism", "different mechanism"),
                             ("phenotype_relation", "different observation"), ("uncertainty_note", "copy-specific domain absence")):
            self.rows = copy.deepcopy(original)
            self.rows[-1][field] = value
            with self.subTest(field=field):
                self.assertEqual(2, len(self.prepare()["reference_function_questions"]["questions"]))
        self.rows = copy.deepcopy(original)
        self.rows[-1]["source_id"] = "unavailable source"
        questions = self.prepare()["reference_function_questions"]["questions"]
        self.assertNotIn("c5", questions[0]["claim_ids"])

    def test_missing_copy_application_remains_unread_not_a_negative(self):
        audit, candidates = self.run_transfer(self.response(omit={"c5"}))
        self.assertEqual(5, len(candidates))
        absent = next(row for row in audit["claims"] if row["claim"]["claim_id"] == "c5")
        self.assertEqual("not_evaluated", absent["claim"]["binding_state"])
        ledger = pd.read_csv(self.root / "candidates/route_b_review_eligibility.tsv", sep="\t", dtype=str)
        self.assertEqual(6, len(ledger))

    def test_a_different_question_cannot_borrow_shared_assessment(self):
        self.rows[-1]["uncertainty_note"] = "different copy domain question"
        self.prepare()
        response = self.response()
        target = next(row for row in response["reviews"] if row["claim_id"] == "c5")
        target["reference_function_question_id"] = next(row["question_id"] for row in response["shared_reference_function_reviews"]
                                                       if row["question_id"] != target["reference_function_question_id"])
        with self.assertRaisesRegex(ValueError, "claim/copy application"):
            AgentClaimReview(str(self.packet(response)), str(self.claims), str(self.fixture.sources_path), route="b")

    def test_changed_frozen_phenotype_question_invalidates_old_reuse(self):
        response = self.response()
        self.rows[0]["phenotype_relation"] = "new focal observation"
        self.prepare()
        with self.assertRaisesRegex(ValueError, "unknown reference-function question"):
            AgentClaimReview(str(self.packet(response)), str(self.claims), str(self.fixture.sources_path), route="b")

    def test_copy_specific_non_applicability_cannot_borrow_positive_function(self):
        response = self.response()
        next(row for row in response["reviews"] if row["claim_id"] == "c5")["copy_application_state"] = "not_applicable"
        audit, candidates = self.run_transfer(response)
        self.assertEqual(5, len(candidates))
        self.assertNotIn("G5", set(candidates.candidate_gene_id))
        absent = next(row for row in audit["claims"] if row["claim"]["claim_id"] == "c5")
        self.assertEqual("review_required", absent["claim"]["binding_state"])
