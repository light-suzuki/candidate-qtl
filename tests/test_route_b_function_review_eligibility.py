"""#901 neutral offline checks; no model predictions or scientific performance claim."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from candidate_qtl.agent_claim_review import B_CONTRACT, FUNCTION_AXES, FUNCTION_PREDICATE, prepare_review_prompt
from candidate_qtl.cli import main
from candidate_qtl.core import sha256
from candidate_qtl.route_b_semantic_candidates import freeze_route_b_semantic_candidates
from candidate_qtl.route_b_semantic_transfer_verify import verify_route_b_semantic_transfer
from candidate_qtl.route_b_semantic_source_verify import verify_route_b_semantic_sources
from candidate_qtl.semantic_claims import CLAIM_COLUMNS
from tests.test_agent_claim_review import ReviewFixture, write_json, spec
from tests.test_route_b_semantic_transfer_verify import _bundle, _claim


class FunctionReviewEligibilityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.fixture = ReviewFixture(self.root)
        self.claim = _claim(claim_id="h1", reference_gene="P12345", source_id="PMID:12345",
                            claim_type="mechanistic_review_seed", biological_function="reference root growth",
                            mechanism="conditional shared cell elongation process in focal stems",
                            phenotype_relation="could shorten focal stem organs; focal causality untested")
        self.review = copy.deepcopy(self.fixture.review)
        self.review.update(claim_id=self.review.pop("hypothesis_id"), predicate_relation=FUNCTION_PREDICATE,
                           reference_function="short roots", focal_bridge=self.claim["phenotype_relation"],
                           axes=dict.fromkeys(FUNCTION_AXES, True))
        self.bundle = self.root / "bundle.json"
        _bundle(self.bundle, members=[("G1", "many-to-many"), ("G2", "many-to-many")])
        self.payload = json.loads(self.bundle.read_text())
        self.payload["b_review_membership_policy"] = "source_function_v1"
        self.payload["phenotype"] = "shortened stems, not the reported reference roots"
        unit = self.payload["units"][0]
        unit["stable_source_ids"] = ["PMID:12345"]
        unit["literature_link_candidates"][0].update(reference_gene="P12345", source_id="PMID:12345")
        for gene in unit["genes"]:
            gene.update(orthology_reference_gene="P12345", reference_species="Arabidopsis thaliana")
        self.proposals = self.root / "proposals.tsv"

    def run_chain(self, *, review=True, source_status="review_required", cli=False):
        write_json(self.bundle, self.payload)
        pd.DataFrame([self.claim], columns=CLAIM_COLUMNS).to_csv(self.proposals, sep="\t", index=False)
        source = self.root / "source.tsv"
        pd.DataFrame([{"claim_id": "h1", "source_relation_status": source_status}]).to_csv(source, sep="\t", index=False)
        prompt = self.root / "review-prompt.json"
        prepare_review_prompt(str(self.proposals), str(self.fixture.sources_path), str(prompt), route="b")
        raw = write_json(self.root / "function-review.json", {"reviews": [self.review] if review else []})
        packet = write_json(self.root / "function-packet.json", {
            "contract": B_CONTRACT, "hypotheses_sha256": sha256(self.proposals),
            "source_review_records_sha256": sha256(self.fixture.sources_path),
            "proposal_run": {"executor_id": "fixture-proposer", "model": "fixture", "raw_response": spec(self.fixture.proposal_raw)},
            "review_run": {"executor_id": "fixture-reviewer", "model": "fixture", "raw_response": spec(raw), "prompt": spec(prompt)},
        })
        before = {p: p.read_bytes() for p in (self.proposals, self.bundle, self.fixture.sources_path, raw)}
        out = self.root / "transfer"
        if cli:
            self.assertEqual(0, main(["route-b-semantic-transfer-verify", "--proposals", str(self.proposals),
                                     "--source-verification", str(source), "--bundle", str(self.bundle),
                                     "--agent-claim-review", str(packet), "--source-review-records", str(self.fixture.sources_path),
                                     "--record-cache", str(self.fixture.records_path), "--out-dir", str(out), "--case-id", "case403"]))
            claims = str(out / "semantic_claims_verified.tsv")
        else:
            claims = verify_route_b_semantic_transfer(str(self.proposals), str(source), str(self.bundle), str(out),
                agent_claim_review_arg=str(packet), source_review_records_arg=str(self.fixture.sources_path),
                record_cache_arg=str(self.fixture.records_path), case_id="case403")["claims"]
        result = freeze_route_b_semantic_candidates(claims, str(self.root / "candidates"), case_id="case403")
        for path, original in before.items():
            self.assertEqual(original, path.read_bytes())
        self.result = result
        self.ledger = pd.read_csv(self.root / "candidates/route_b_review_eligibility.tsv", sep="\t", dtype=str).fillna("")
        return pd.read_csv(result["review_candidates"], sep="\t", dtype=str).fillna("")

    def test_type3_function_review_reaches_cli_and_retains_all_ambiguous_copies(self):
        rows = self.run_chain(cli=True)
        self.assertEqual({"G1", "G2"}, set(rows.candidate_gene_id))
        self.assertTrue(rows.support_status.eq("unverified_route_b_review").all())
        self.assertEqual(0, self.result["verified_claim_count"])
        self.assertTrue(self.ledger.b_review_eligibility.eq("eligible_reference_function_review").all())

    def test_citation_alone_and_unread_sources_are_ledger_only(self):
        self.assertTrue(self.run_chain(review=False).empty)
        self.assertTrue(self.ledger.b_review_eligibility.eq("not_evaluated").all())
        self.assertEqual({"G1", "G2"}, set(self.ledger.candidate_gene_id))
        self.assertTrue(self.ledger.biological_exclusion.eq("not_asserted").all())

    def test_function_review_does_not_require_focal_trait_match_in_source_triage(self):
        self.assertEqual({"G1", "G2"}, set(self.run_chain(source_status="unsupported").candidate_gene_id))
        self.assertEqual(0, self.result["verified_claim_count"])

    def test_wrong_reference_species_or_function_or_generic_bridge_cannot_enter(self):
        for change in ({"reference_function": "unrelated unquoted function"},
                       {"reported_gene_identity": "unrelated gene"}, {"organism": "wrong species"},
                       {"axes": {**dict.fromkeys(FUNCTION_AXES, True), "specific_focal_bridge": False}}):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as tmp:
                original_root, original_review = self.root, self.review
                self.root = Path(tmp)
                self.review = {**original_review, **change}
                self.assertTrue(self.run_chain().empty)
                self.root, self.review = original_root, original_review

    def test_source_absence_and_retrieval_failure_remain_distinct(self):
        self.claim["source_id"] = ""
        rows = self.run_chain(review=False, source_status="not_evaluated")
        self.assertTrue(rows.empty)
        self.assertTrue(self.ledger.b_review_eligibility.eq("source_absent").all())

    def test_zero_b_eligible_still_freezes_empty_wet_and_full_unresolved_inventory(self):
        from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff

        self.claim["source_id"] = ""
        self.run_chain(review=False, source_status="not_evaluated")
        atlas = self.root / "empty-eligible-atlas.tsv"
        pd.DataFrame([{"gene_id": g, "family_id": "EV1"} for g in ("G1", "G2")]).to_csv(atlas, sep="\t", index=False)
        result = wet_candidate_handoff(str(self.root / "wet-empty-eligible"),
            route_b_candidates_arg=self.result["candidates"], route_b_candidates_manifest_arg=self.result["manifest"],
            review_candidates_arg=self.result["review_candidates"], review_candidates_manifest_arg=self.result["review_manifest"],
            interval_atlas_arg=str(atlas), case_id="case403")
        self.assertEqual("review_required_no_supported_candidates", result["handoff_state"])
        self.assertTrue(pd.read_csv(result["wet_candidates"], sep="\t").empty)
        self.assertFalse(result["heavy_stage_invoked"])
        inventory = pd.read_csv(self.root / "wet-empty-eligible/wet_relevance_inventory.tsv", sep="\t", dtype=str)
        self.assertEqual({"G1", "G2"}, set(inventory.candidate_gene_id))
        self.assertTrue(inventory.b_review_eligibility.eq("source_absent").all())

    def test_paralog_group_review_does_not_require_exact_focal_copy_endpoint(self):
        for gene in self.payload["units"][0]["genes"]:
            gene["orthology_reference_gene"] = "other endpoint"
        self.assertEqual({"G1", "G2"}, set(self.run_chain().candidate_gene_id))
        self.assertEqual(0, self.result["verified_claim_count"])

    def test_exact_card_mapping_is_usable_without_strict_orthology(self):
        unit = self.payload["units"][0]
        unit["literature_link_candidates"] = []
        for gene in unit["genes"]:
            gene["orthology_reference_gene"] = "other endpoint"
        unit["reference_material_cards"] = [{"interval_gene_id": g, "reference_gene": "P12345",
                                             "reference_species": "Arabidopsis thaliana", "source_id": "PMID:12345"}
                                            for g in ("G1", "G2")]
        self.assertEqual({"G1", "G2"}, set(self.run_chain().candidate_gene_id))

    def test_n1000_generic_bridge_is_not_admitted_without_specific_semantic_review(self):
        unit = self.payload["units"][0]
        unit["member_interval_genes"] = [f"G{i}" for i in range(1, 1001)]
        unit["genes"].extend({"gene_id": f"G{i}", "orthology_reference_gene": f"other-ref-{i}",
                              "reference_species": "Arabidopsis thaliana", "domain_architecture": "generic kinase"}
                             for i in range(3, 1001))
        self.claim["member_interval_genes"] = ";".join(unit["member_interval_genes"])
        self.review["axes"]["specific_focal_bridge"] = False
        rows = self.run_chain(cli=True)
        self.assertTrue(rows.empty)
        self.assertEqual(1000, len(self.ledger))
        self.assertTrue(self.ledger.b_review_eligibility.eq("review_not_authorized").all())
        from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff

        atlas = self.root / "atlas-n1000.tsv"
        pd.DataFrame([{"gene_id": gene, "family_id": "EV1", "domain_architecture": "generic kinase"}
                      for gene in unit["member_interval_genes"]]).to_csv(atlas, sep="\t", index=False)
        wet_result = wet_candidate_handoff(str(self.root / "wet-n1000"),
            route_b_candidates_arg=self.result["candidates"], route_b_candidates_manifest_arg=self.result["manifest"],
            review_candidates_arg=self.result["review_candidates"], review_candidates_manifest_arg=self.result["review_manifest"],
            interval_atlas_arg=str(atlas), case_id="case403")
        wet = pd.read_csv(wet_result["wet_candidates"], sep="\t", dtype=str)
        self.assertTrue(wet.empty)
        inventory = pd.read_csv(self.root / "wet-n1000/wet_relevance_inventory.tsv", sep="\t", dtype=str)
        self.assertEqual(1000, len(inventory))

    def test_legacy_bundle_replays_original_grounded_review_policy(self):
        self.payload.pop("b_review_membership_policy")
        self.claim["source_id"] = ""
        self.assertEqual({"G1", "G2"}, set(self.run_chain(review=False).candidate_gene_id))

    def test_b_source_absence_has_no_veto_over_a_and_c_independent_support(self):
        from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
        from tests._issue664_helpers import manifest_for

        self.claim["source_id"] = ""
        self.run_chain(review=False)
        atlas = self.root / "atlas.tsv"
        pd.DataFrame([{"gene_id": g, "family_id": "EV1"} for g in ("G1", "G2", "G3")]).to_csv(atlas, sep="\t", index=False)
        a = self.root / "a.tsv"
        pd.DataFrame([{"route": "a", "gene_id": "G1", "family_id": "EV1", "hypothesis_id": "A1",
                       "mapping_tier": "tier3_sequence_match_unique", "ambiguity_status": "unique", "review_status": "supported"}]).to_csv(a, sep="\t", index=False)
        c = self.root / "c.tsv"
        pd.DataFrame([{"candidate_gene_id": "G2", "record_id": "C1", "source_id": "PMID:2", "family_ids": "EV1",
                       "species_distance": "S1_same_species", "evidence_directness": "D0_causal_perturbation_or_rescue",
                       "projection_state": "inside_current_qtl", "projection_basis": "exact_focal_id",
                       "support_status": "verified_direct_support", "rank_contribution": "none"}]).to_csv(c, sep="\t", index=False)
        result = wet_candidate_handoff(str(self.root / "wet"), route_a_bridge_arg=str(a),
            route_a_manifest_arg=str(manifest_for(a, route="a", case_id="case403")),
            route_c_members_arg=str(c), route_c_manifest_arg=str(manifest_for(c, route="c", case_id="case403")),
            interval_atlas_arg=str(atlas), review_candidates_arg=self.result["review_candidates"],
            review_candidates_manifest_arg=self.result["review_manifest"], case_id="case403")
        wet = pd.read_csv(result["wet_candidates"], sep="\t", dtype=str)
        self.assertEqual({"G1", "G2"}, set(wet.candidate_gene_id))
        inventory = pd.read_csv(self.root / "wet/wet_relevance_inventory.tsv", sep="\t", dtype=str).fillna("")
        self.assertEqual(3, len(inventory))
        self.assertTrue(inventory[inventory.candidate_gene_id.isin({"G1", "G2"})].b_review_eligibility.eq("source_absent").all())

    def test_curated_function_is_captured_and_independently_reviewable_without_paper(self):
        source_id = "UniProt:P12345"
        self.claim["source_id"] = source_id
        self.payload["units"][0]["stable_source_ids"] = [source_id]
        unit = self.payload["units"][0]
        unit["source_provenance"] = [{"gene_id": "G1", "reference_gene": "P12345", "source_ids": [source_id],
                                      "source_database": "UniProt", "source_url": "https://example.invalid/curated",
                                      "reference_species": "Arabidopsis thaliana"}]
        for gene in unit["genes"]:
            gene.update(reference_gene="P12345", reference_source_id=source_id, lookup_status="source_backed",
                        reference_function="P12345 controls cell elongation in Arabidopsis thaliana")
        write_json(self.bundle, self.payload)
        pd.DataFrame([self.claim], columns=CLAIM_COLUMNS).to_csv(self.proposals, sep="\t", index=False)
        source = verify_route_b_semantic_sources(str(self.proposals), str(self.root / "curated-sources"), bundle_arg=str(self.bundle))
        self.fixture.sources_path = Path(source["source_review_records"])
        captured = json.loads(self.fixture.sources_path.read_text())["records"][source_id]
        text = captured["sections"]["curatedRecord"]
        self.review.update(citation=source_id, source_record_id=source_id, section="curatedRecord", quoted_support=text,
                           quote_start=0, quote_end=len(text), sections_sha256=captured["sections_sha256"],
                           reported_gene_identity="P12345", organism_quote="Arabidopsis thaliana",
                           reference_function="controls cell elongation")
        self.assertEqual({"G1", "G2"}, set(self.run_chain(source_status="supported_curated_function").candidate_gene_id))

    def test_eligibility_table_tamper_is_rejected_by_candidate_consumer(self):
        self.run_chain()
        path = self.root / "transfer/semantic_claim_transfer_verification.tsv"
        path.write_bytes(path.read_bytes() + b" ")
        with self.assertRaisesRegex(ValueError, "Route B transfer-verification input hash mismatch"):
            freeze_route_b_semantic_candidates(str(self.root / "transfer/semantic_claims_verified.tsv"), str(self.root / "tamper"), case_id="case403")
        self.assertFalse((self.root / "tamper").exists(), "hash rejection must precede candidate output")

    def test_reference_card_function_reaches_independent_review_when_gene_row_has_only_ids(self):
        source_id = "UniProt:P12345"
        self.claim["source_id"] = source_id
        unit = self.payload["units"][0]
        unit["stable_source_ids"] = [source_id]
        unit["source_provenance"] = [{"reference_gene": "P12345", "source_ids": [source_id],
                                      "source_url": "https://example.invalid/curated"}]
        unit["reference_material_cards"] = [{"reference_gene": "P12345", "source_id": source_id,
            "reference_species": "Arabidopsis thaliana", "source_database": "UniProt",
            "source_url": "https://example.invalid/curated", "function_description": "controls cell elongation"}]
        write_json(self.bundle, self.payload)
        pd.DataFrame([self.claim], columns=CLAIM_COLUMNS).to_csv(self.proposals, sep="\t", index=False)
        source = verify_route_b_semantic_sources(str(self.proposals), str(self.root / "card-sources"), bundle_arg=str(self.bundle))
        self.fixture.sources_path = Path(source["source_review_records"])
        captured = json.loads(self.fixture.sources_path.read_text())["records"][source_id]
        text = captured["sections"]["curatedRecord"]
        self.assertEqual(unit["reference_material_cards"][0], json.loads(text)["reference_material_card"])
        self.review.update(citation=source_id, source_record_id=source_id, section="curatedRecord", quoted_support=text,
                           quote_start=0, quote_end=len(text), sections_sha256=captured["sections_sha256"],
                           reported_gene_identity="P12345", organism_quote="Arabidopsis thaliana",
                           reference_function="controls cell elongation")
        self.assertEqual({"G1", "G2"}, set(self.run_chain(source_status="supported_curated_function", cli=True).candidate_gene_id))

    def test_reference_card_id_or_uninformative_content_cannot_stand_in_for_function(self):
        for field, value in (("function_description", "unknown"), ("reference_gene", "OTHER"), ("source_id", "UniProt:OTHER")):
            with self.subTest(field=field):
                source_id = "UniProt:P12345"
                self.claim["source_id"] = source_id
                unit = self.payload["units"][0]
                unit["stable_source_ids"] = [source_id]
                unit["source_provenance"] = [{"reference_gene": "P12345", "source_ids": [source_id],
                                              "source_url": "https://example.invalid/curated"}]
                unit["reference_material_cards"] = [{"reference_gene": "P12345", "source_id": source_id,
                    "source_url": "https://example.invalid/curated", "function_description": "controls cell elongation", field: value}]
                write_json(self.bundle, self.payload)
                pd.DataFrame([self.claim], columns=CLAIM_COLUMNS).to_csv(self.proposals, sep="\t", index=False)
                result = verify_route_b_semantic_sources(str(self.proposals), str(self.root / field), bundle_arg=str(self.bundle))
                records = json.loads(Path(result["source_review_records"]).read_text())["records"]
                self.assertNotIn(source_id, records)
