"""#688 contract fixtures, NOT biological benchmark/model predictions."""

from __future__ import annotations

import copy
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.agent_claim_review import AgentClaimReview, AXES, CONTRACT, METHOD, prepare_review_prompt
from candidate_qtl.claim_binding import bound_claim_authorizes, validate_bound_claim_frame
from candidate_qtl.cli import main
from candidate_qtl.core import sha256
from candidate_qtl.guided_profiles import _plan, _execute_stage
from candidate_qtl.llm_pathway import hypothesis_verify, bait_panel_freeze
from candidate_qtl.source_verify import source_verify
from tests.bound_claim_fixtures import bound_claim


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, sort_keys=True), encoding="utf-8")
    return path


def spec(path):
    return {"path": str(path), "sha256": sha256(path)}


class ReviewFixture:
    def __init__(self, root, reference_name="XYZ1"):
        self.root = root
        self.hypotheses = root / "mechanisms.tsv"
        pd.DataFrame([{
            "hypothesis_id": "h1", "reference_gene_name": reference_name,
            "reference_species": "Arabidopsis thaliana", "proposed_stable_id": "P12345",
            "pathway_or_process": "organ elongation", "mechanistic_link": "loss causes short roots",
            "source_ids": "PMID:12345",
        }]).to_csv(self.hypotheses, sep="\t", index=False)
        self.text = "In Arabidopsis thaliana, loss of XYZ1 caused short roots. Complementation restored root length."
        self.cache = write_json(root / "source-cache.json", {
            "id": "12345", "pmid": "12345", "title": "XYZ1 root development",
            "abstractText": self.text,
        })
        source = source_verify(str(self.hypotheses), str(root / "source"), cache_arg=str(self.cache))
        self.hypotheses = Path(source["mechanisms_with_verified_sources"])
        self.source_records = Path(source["verified_sources"])
        self.sources_path = Path(source["source_review_records"])
        self.sources = json.loads(self.sources_path.read_text())
        self.record = {
            "primaryAccession": "P12345", "organism": {"scientificName": "Arabidopsis thaliana"},
            "sequence": {"value": "MAAA"}, "genes": [{"geneName": {"value": "XYZ1"}}],
            "references": [{"citation": {"citationCrossReferences": [{"database": "PubMed", "id": "12345"}]}}],
        }
        self.records_path = root / "records.json"
        write_json(self.records_path, self.record)
        self.review = {
            "hypothesis_id": "h1", "citation": "PMID:12345", "source_record_id": "12345",
            "sections_sha256": self.sources["records"]["PMID:12345"]["sections_sha256"],
            "subject_gene_identity": "P12345", "reported_gene_identity": "XYZ1",
            "organism": "Arabidopsis thaliana", "reported_organism": "Arabidopsis thaliana",
            "organism_quote": "In Arabidopsis thaliana", "section": "abstractText",
            "quoted_support": self.text, "quote_start": 0, "quote_end": len(self.text),
            "experiment_perturbation": "loss of XYZ1", "reported_phenotype": "short roots",
            "predicate_relation": "causal_perturbation_supports_target_trait",
            "decision": "supported", "reason": "Synthetic fixture: same gene loss and complementation bind the phenotype.",
            "axes": dict.fromkeys(AXES, True),
        }
        self.proposal_raw = write_json(root / "proposal.json", {"fixture": "proposal; not real model execution"})
        self.prompt_path = root / "prompt.json"
        self.raw_path = root / "review.json"
        self.packet_path = root / "packet.json"
        self.seal()

    def seal(self, reviews=None):
        self.sources = json.loads(self.sources_path.read_text())
        prepare_review_prompt(str(self.hypotheses), str(self.sources_path), str(self.prompt_path))
        write_json(self.raw_path, {"reviews": [self.review] if reviews is None else reviews})
        self.packet = {
            "contract": CONTRACT, "hypotheses_sha256": sha256(self.hypotheses),
            "source_review_records_sha256": sha256(self.sources_path),
            "proposal_run": {"executor_id": "fixture-proposer", "model": "fixture-not-model",
                             "raw_response": spec(self.proposal_raw)},
            "review_run": {"executor_id": "fixture-reviewer", "model": "fixture-not-model",
                           "raw_response": spec(self.raw_path), "prompt": spec(self.prompt_path)},
        }
        write_json(self.packet_path, self.packet)

    def checker(self):
        return AgentClaimReview(str(self.packet_path), str(self.hypotheses), str(self.sources_path))

    def bind(self, record=None):
        checker = self.checker()
        return checker.bind({"hypothesis_id": "h1"}, "P12345", self.record if record is None else record,
                            ["PMID:12345"])

    def verify(self, out="verify", **kwargs):
        return hypothesis_verify(str(self.hypotheses), str(self.root / out), "Pisum sativum",
                                 record_cache_arg=str(self.records_path),
                                 source_records_arg=str(self.source_records),
                                 agent_claim_review_arg=str(self.packet_path),
                                 source_review_records_arg=str(self.sources_path), **kwargs)


class AgentClaimReviewTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.fixture = ReviewFixture(Path(self.tmp.name))

    def test_benchmark_truth_columns_are_rejected_before_prompt_export(self):
        fixture = self.fixture
        rows = pd.read_csv(fixture.hypotheses, sep="\t", dtype=str).fillna("")
        rows["truth_gene_rank_le_20"] = "1"
        rows["answerGeneId"] = "BENCH_ONLY"
        rows.to_csv(fixture.hypotheses, sep="\t", index=False)
        with self.assertRaisesRegex(ValueError, "benchmark truth/answer fields"):
            prepare_review_prompt(
                str(fixture.hypotheses), str(fixture.sources_path),
                str(fixture.root / "leaked-prompt.json"),
            )
        self.assertFalse((fixture.root / "leaked-prompt.json").exists())

    def test_nested_benchmark_truth_is_rejected_before_prompt_export(self):
        fixture = self.fixture
        source = json.loads(fixture.sources_path.read_text(encoding="utf-8"))
        source["records"]["PMID:12345"]["truth_variant"] = "BENCH_ONLY"
        source["records"]["PMID:12345"]["truthVariant"] = "BENCH_ONLY"
        fixture.sources_path.write_text(json.dumps(source), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "benchmark truth/answer fields"):
            prepare_review_prompt(
                str(fixture.hypotheses), str(fixture.sources_path),
                str(fixture.root / "nested-leaked-prompt.json"),
            )
        self.assertFalse((fixture.root / "nested-leaked-prompt.json").exists())

    def test_positive_cli_to_bait_and_method_provenance(self):
        f = self.fixture
        out = f.root / "cli"
        self.assertEqual(0, main([
            "hypothesis-verify", "--hypotheses", str(f.hypotheses), "--out-dir", str(out),
            "--focal-species", "Pisum sativum", "--record-cache", str(f.records_path),
            "--source-records", str(f.source_records), "--agent-claim-review", str(f.packet_path),
            "--source-review-records", str(f.sources_path), "--json",
        ]))
        validation = pd.read_csv(out / "proposal_validation.tsv", sep="\t")
        self.assertTrue(validation.iloc[0].bound_claim_authorized)
        self.assertEqual(METHOD, validation.iloc[0].bound_claim_curation_method)
        freeze = bait_panel_freeze(str(f.hypotheses), str(out / "proposal_validation.tsv"),
                                  str(out / "verified_sequences.json"), str(f.root / "freeze"))
        panel = pd.read_csv(freeze["panel"], sep="\t")
        self.assertEqual(["P12345"], panel.bait_id.tolist())
        self.assertEqual([METHOD], panel.bound_claim_curation_method.tolist())
        self.assertEqual(["none"], panel.rank_contribution.tolist())
        audit = json.loads((out / "agent_claim_review_audit.json").read_text())
        self.assertEqual("independent_agent_attestation_not_deterministic_proof", audit["semantic_verification"])
        self.assertEqual(f.record, audit["claims"][0]["accession_record"])
        manifest = json.loads((out / "hypothesis_verification_manifest.json").read_text())
        self.assertEqual(sha256(out / "agent_claim_review_audit.json"), manifest["agent_claim_review_audit"]["sha256"])

    def test_each_missing_semantic_axis_stays_review_only(self):
        f = self.fixture
        for axis in AXES:
            with self.subTest(axis=axis):
                f.review["axes"] = dict.fromkeys(AXES, True)
                f.review["axes"][axis] = False
                f.seal()
                self.assertFalse(bound_claim_authorizes(f.bind())[0])

    def test_independent_review_can_resolve_lexical_triage_without_rewriting_it(self):
        directory = self.fixture.root / "lexical"
        directory.mkdir()
        f = ReviewFixture(directory, reference_name="P12345")
        before = f.source_records.read_bytes()
        source_status = pd.read_csv(f.hypotheses, sep="\t").iloc[0]
        self.assertNotEqual("mechanistic_support", source_status.claim_support_summary)
        result = f.verify()
        self.assertEqual(1, result["verified_count"])
        self.assertTrue(pd.read_csv(result["validation"], sep="\t").iloc[0].bound_claim_authorized)
        self.assertEqual(before, f.source_records.read_bytes())

    def test_missing_source_or_wrong_species_never_freezes(self):
        f = self.fixture
        for mutation in ({"citation": "PMID:99999"}, {"organism": "Oryza sativa"}):
            f.review.update(mutation)
            f.seal()
            result = f.verify()
            freeze = bait_panel_freeze(str(f.hypotheses), result["validation"], result["verified_sequences"],
                                      str(f.root / "freeze"))
            self.assertEqual(0, freeze["bait_count"])

    def test_wrong_text_gene_species_experiment_and_scope_fail_closed(self):
        f = self.fixture
        original = copy.deepcopy(f.review)
        for field, value in (
            ("subject_gene_identity", "Q99999"), ("reported_gene_identity", "OTHER"),
            ("organism", "Oryza sativa"), ("reported_organism", "Oryza sativa"),
            ("organism_quote", "not retrieved"), ("source_record_id", "54321"),
            ("sections_sha256", "0" * 64), ("quote_start", 1),
            ("quoted_support", "Fabricated source"), ("section", "title"),
            ("experiment_perturbation", "loss of OTHER"), ("reported_phenotype", "long roots"),
            ("predicate_relation", "curated_function_associated"),
        ):
            with self.subTest(field=field):
                f.review = {**original, field: value}
                f.seal()
                self.assertFalse(bound_claim_authorizes(f.bind())[0])

    def test_symbol_without_provider_publication_link_cannot_authorize_id(self):
        f = self.fixture
        claim = f.checker().bind({"hypothesis_id": "h1"}, "P12345", f.record, [])
        self.assertFalse(bound_claim_authorizes(claim)[0])
        self.assertIn("provider citation", claim["binding_reason"])

    def test_primary_accession_and_provider_alias_checked(self):
        for modified in ({**self.fixture.record, "primaryAccession": "Q99999"},
                         {**self.fixture.record, "genes": []}):
            self.assertFalse(bound_claim_authorizes(self.fixture.bind(modified))[0])

    def test_malformed_review_types_are_explicit_not_evaluated(self):
        f = self.fixture
        original = copy.deepcopy(f.review)
        for field, value in (("subject_gene_identity", None), ("organism", []), ("axes", True)):
            f.review = {**original, field: value}
            f.seal()
            self.assertEqual("not_evaluated", bound_claim_authorizes(f.bind())[1])

    def test_same_genus_different_declared_species_is_not_enough(self):
        f = self.fixture
        hypotheses = pd.read_csv(f.hypotheses, sep="\t", dtype=str).fillna("")
        hypotheses.loc[0, "reference_species"] = "Arabidopsis lyrata"
        hypotheses.to_csv(f.hypotheses, sep="\t", index=False)
        f.seal()
        self.assertFalse(bound_claim_authorizes(f.bind())[0])

    def test_species_elsewhere_in_section_is_not_the_claim_species(self):
        f = self.fixture
        f.review["quoted_support"] = "loss of XYZ1 caused short roots."
        f.review["quote_start"] = f.text.index(f.review["quoted_support"])
        f.review["quote_end"] = f.review["quote_start"] + len(f.review["quoted_support"])
        f.seal()
        self.assertFalse(bound_claim_authorizes(f.bind())[0])

    def test_freeze_rechecks_validation_and_audit_hashes(self):
        f = self.fixture
        result = f.verify()
        audit = Path(result["agent_claim_review_audit"])
        audit.write_text(audit.read_text() + "\n")
        with self.assertRaisesRegex(ValueError, "differ from the verified"):
            bait_panel_freeze(str(f.hypotheses), result["validation"], result["verified_sequences"],
                              str(f.root / "freeze"))

    def test_provider_outcome_facade_reseals_final_validation(self):
        from candidate_qtl.llm_pathway import _apply_entry_outcomes
        from candidate_qtl.uniprot_error_policy import UniProtEntryEvent
        f = self.fixture
        result = f.verify()
        result = _apply_entry_outcomes(result, [UniProtEntryEvent("P12345", "success", 200, True, "fixture")])
        freeze = bait_panel_freeze(str(f.hypotheses), result["validation"], result["verified_sequences"],
                                   str(f.root / "freeze"))
        self.assertEqual(1, freeze["bait_count"])

    def test_raw_method_label_and_mutated_receipt_cannot_authorize(self):
        row = bound_claim(route="a", claim_id="h1", citation="PMID:12345", subject_gene="P12345",
                          organism="Arabidopsis thaliana", curation_method=METHOD)
        with self.assertRaises(ValueError):
            validate_bound_claim_frame(pd.DataFrame([row]))
        row = self.fixture.bind()
        self.assertTrue(bound_claim_authorizes(row)[0])
        row["object_phenotype_process"] = "other phenotype"
        self.assertFalse(bound_claim_authorizes(row)[0])

    def test_same_executor_and_tampered_raw_response_rejected(self):
        f = self.fixture
        f.packet["review_run"]["executor_id"] = "fixture-proposer"
        write_json(f.packet_path, f.packet)
        with self.assertRaisesRegex(ValueError, "separate executor"):
            f.checker()
        f.seal()
        write_json(f.raw_path, {"reviews": []})
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            f.checker()

    def test_prompt_must_contain_exact_inputs(self):
        f = self.fixture
        prompt = json.loads(f.prompt_path.read_text())
        prompt["source_review_records"] = {}
        write_json(f.prompt_path, prompt)
        f.packet["review_run"]["prompt"] = spec(f.prompt_path)
        write_json(f.packet_path, f.packet)
        with self.assertRaisesRegex(ValueError, "exact frozen source"):
            f.checker()

    def test_missing_review_is_not_evaluated_and_never_freezes(self):
        f = self.fixture
        f.seal([])
        result = f.verify()
        val = pd.read_csv(result["validation"], sep="\t").iloc[0]
        self.assertEqual("not_evaluated", val.bound_claim_state)
        self.assertEqual("review_required", val.verification_status)
        panel = bait_panel_freeze(str(f.hypotheses), result["validation"], result["verified_sequences"],
                                 str(f.root / "freeze"))
        self.assertTrue(pd.read_csv(panel["panel"], sep="\t").empty)

    def test_missing_record_keeps_provider_gap_and_never_promotes(self):
        f = self.fixture
        write_json(f.records_path, {})
        result = f.verify()
        val = pd.read_csv(result["validation"], sep="\t").iloc[0]
        self.assertEqual("not_evaluated", val.bound_claim_state)
        self.assertEqual("review_required", val.verification_status)

    def test_duplicate_unknown_or_replayed_hypotheses_rejected(self):
        f = self.fixture
        for reviews in ([f.review, f.review], [{**f.review, "hypothesis_id": "other"}]):
            f.seal(reviews)
            with self.assertRaisesRegex(ValueError, "duplicate or unknown"):
                f.checker()
        f.seal()
        with f.hypotheses.open("a") as stream:
            stream.write("\n")
        with self.assertRaisesRegex(ValueError, "frozen hypotheses"):
            f.checker()

    def test_guided_chain_wires_sources_and_dispatches_existing_stage(self):
        f = self.fixture
        stages, _ = _plan("triple-route", {"agent_claim_review": str(f.packet_path)}, f.root / "plan")
        stage = next(row for row in stages if row["command"] == "hypothesis-verify")
        self.assertIn("source_review_records", stage["inputs_from"])
        result = _execute_stage({"command": "hypothesis-verify", "arguments": {
            "hypotheses": str(f.hypotheses), "out_dir": str(f.root / "guided"),
            "focal_species": "Pisum sativum", "record_cache": str(f.records_path),
            "source_records": str(f.source_records), "agent_claim_review": str(f.packet_path),
            "source_review_records": str(f.sources_path),
        }})
        self.assertEqual(1, result["verified_count"])


if __name__ == "__main__":
    unittest.main()
