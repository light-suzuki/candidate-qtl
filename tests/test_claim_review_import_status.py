"""Frozen B import coverage survives sessions without promoting scientific support."""
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.bibliography import INITIAL_POLICY
from candidate_qtl.claim_review_session import complete_session, prepare_session, sha256
from candidate_qtl.route_b_semantic_import import import_route_b_semantic_proposals
from candidate_qtl.route_b_semantic_source_verify import verify_route_b_semantic_sources
from candidate_qtl.wet_candidate_handoff import _initial_b_coverage, wet_candidate_handoff
from tests import test_route_b_semantic_import as import_tests


class ClaimReviewImportStatusTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def inputs(self, *, failed=False, skipped=False, no_claims=False, with_llm_manifest=False):
        fixture = import_tests.RouteBSemanticImportTests()
        built = fixture._bundle(self.root)
        if skipped:
            bundle = json.loads(Path(built["bundle"]).read_text())
            bundle["units"].append({"unit_id": "skipped", "genes": [],
                                    "request_state": "not_requested_no_reviewable_material"})
            Path(built["bundle"]).write_text(json.dumps(bundle))
        response = self.root / "response.json"
        claims = [{
            "claim_type": "mutant_phenotype", "reference_gene": "REF1",
            "reference_species": "Arabidopsis thaliana", "evidence_scope": "reference_gene",
            "source_ids": ["PMID:123"], "biological_function": "meristem maintenance",
            "mechanism": "receptor signaling", "phenotype_relation": "could explain fasciation",
        }]
        fixture._response(Path(built["prompts"]), [] if no_claims else claims, response)
        if failed:
            value = json.loads(response.read_text())
            value["candidates"][0].update(parse_status="not_evaluated", model_json=None, claims=[])
            response.write_text(json.dumps(value))
        llm = None
        if with_llm_manifest:
            llm = self.root / "llm_run_manifest.json"
            llm.write_text(json.dumps({"provider": "synthetic", "model": "fixture",
                "bundle_sha256": sha256(built["bundle"]), "prompt_sha256": sha256(built["prompts"]),
                "response_sha256": sha256(response), "execution_evidence": "synthetic_fixture",
                "counts_as_per_prompt_evaluation": False}))
        imported = import_route_b_semantic_proposals(built["bundle"], built["prompts"],
            str(response), str(self.root / "import"), case_id="case403",
            llm_manifest_arg=str(llm) if llm else None)
        cache = self.root / "cache.json"
        cache.write_text(json.dumps({"id": "123", "title": "Synthetic meristem study"}))
        source = verify_route_b_semantic_sources(imported["proposals"], str(self.root / "sources"),
            bundle_arg=built["bundle"], cache_arg=str(cache), case_id="case403")
        records = self.root / "records.json"
        records.write_text("{}")
        return {"hypotheses": imported["proposals"], "bundle": built["bundle"],
                "source_verification": source["source_verification"],
                "sources": source["source_review_records"], "proposal": str(response), "records": str(records)}

    def prepare(self, inputs):
        return prepare_session(str(self.root / "session"), route="b", inputs=inputs,
            proposer_id="synthetic-b", proposer_model="fixture", case_id="case403", review_policy=INITIAL_POLICY)

    def finish(self, inputs):
        prepared = self.prepare(inputs)
        result = complete_session(str(self.root / "session"), prepared["request_sha256"])
        producer = json.loads(Path(result["candidates"]["review_manifest"]).read_text())
        _, observations = _initial_b_coverage([producer])
        manifest = json.loads(Path(result["verification"]["manifest"]).read_text())
        return prepared, result, manifest, observations

    def test_status_survives_real_prepare_complete_and_coverage_consumer(self):
        inputs = self.inputs()
        original = Path(inputs["hypotheses"]).with_name("semantic_claim_proposal_status.tsv")
        prepared = self.prepare(inputs)
        request = json.loads((self.root / "session/request.json").read_text())
        snapshot = self.root / "session/inputs/semantic_claim_proposal_status.tsv"
        self.assertEqual(original.read_bytes(), snapshot.read_bytes())
        self.assertIn(str(snapshot.relative_to(self.root / "session")), request["files"])
        original.write_text("upstream removed after snapshot")
        result = complete_session(str(self.root / "session"), prepared["request_sha256"])
        transfer = json.loads(Path(result["verification"]["manifest"]).read_text())
        self.assertEqual({"path": str(snapshot), "sha256": sha256(snapshot)}, transfer["ancestors"]["proposal_status"])
        producer = json.loads(Path(result["candidates"]["review_manifest"]).read_text())
        _, observations = _initial_b_coverage([producer])
        self.assertEqual("import_complete", observations[0]["status"])
        self.assertEqual(0, result["candidates"]["verified_claim_count"])
        wet = wet_candidate_handoff(str(self.root / "wet"), case_id="case403",
            route_b_claims_arg=result["verification"]["claims"],
            route_b_claims_manifest_arg=result["verification"]["manifest"],
            review_candidates_arg=result["candidates"]["review_candidates"],
            review_candidates_manifest_arg=result["candidates"]["review_manifest"],
            interval_atlas_arg=str(self.root / "atlas.tsv"))
        handoff = json.loads(Path(wet["manifest"]).read_text())
        coverage = handoff["input_coverage_detail"]["initial_b_import_and_source_observations"]
        self.assertEqual("import_complete", coverage[0]["status"])
        self.assertEqual("single_route_by_design", wet["input_coverage_state"])
        self.assertEqual(2, len(pd.read_csv(wet["wet_candidates"], sep="\t")))

    def test_failed_and_unrequested_units_remain_incomplete(self):
        for failed, skipped in ((True, False), (False, True)):
            with self.subTest(failed=failed, skipped=skipped), tempfile.TemporaryDirectory() as temporary:
                self.root = Path(temporary)
                _, result, manifest, observations = self.finish(self.inputs(failed=failed, skipped=skipped))
                self.assertIsNotNone(manifest["ancestors"]["proposal_status"])
                self.assertEqual("import_incomplete", observations[0]["status"])
                self.assertEqual(0, result["candidates"]["verified_claim_count"])

    def test_nonempty_llm_provenance_survives_prepare_complete(self):
        inputs = self.inputs(with_llm_manifest=True)
        llm = self.root / "llm_run_manifest.json"
        original = llm.read_bytes()
        _, result, _, observations = self.finish(inputs)
        self.assertEqual(original, llm.read_bytes())
        self.assertEqual(original, (self.root / "session/inputs/llm_provenance/manifest.json").read_bytes())
        frame = pd.read_csv(self.root / "session/inputs/hypotheses.tsv", sep="\t").fillna("")
        self.assertEqual({"synthetic"}, set(frame["provider"]))
        self.assertEqual({"fixture"}, set(frame["model"]))
        self.assertEqual("import_complete", observations[0]["status"])
        self.assertEqual(0, result["candidates"]["verified_claim_count"])

    def test_llm_provenance_snapshot_is_self_contained(self):
        inputs = self.inputs(with_llm_manifest=True)
        prepared = self.prepare(inputs)
        (self.root / "llm_run_manifest.json").unlink()
        result = complete_session(str(self.root / "session"), prepared["request_sha256"])
        self.assertEqual(0, result["candidates"]["verified_claim_count"])

    def test_llm_manifest_basename_cannot_overwrite_session_inputs(self):
        for name in ("bundle.json", "sources.json", "semantic_claim_proposal_manifest.json"):
            with self.subTest(name=name), tempfile.TemporaryDirectory() as temporary:
                self.root = Path(temporary)
                inputs = self.inputs(with_llm_manifest=True)
                original = self.root / "llm_run_manifest.json"
                relocated = self.root / "provider" / name
                relocated.parent.mkdir()
                original.rename(relocated)
                import_manifest = Path(inputs["hypotheses"]).with_name("semantic_claim_proposal_manifest.json")
                value = json.loads(import_manifest.read_text())
                value["llm_manifest"] = str(relocated)
                import_manifest.write_text(json.dumps(value))
                prepared = self.prepare(inputs)
                snapshot = self.root / "session/inputs"
                self.assertEqual(Path(inputs["bundle"]).read_bytes(), (snapshot / "bundle.json").read_bytes())
                self.assertEqual(relocated.read_bytes(), (snapshot / "llm_provenance/manifest.json").read_bytes())
                relocated.unlink()
                result = complete_session(str(self.root / "session"), prepared["request_sha256"])
                self.assertEqual(0, result["candidates"]["verified_claim_count"])

    def test_original_llm_provenance_tampering_is_rejected(self):
        inputs = self.inputs(with_llm_manifest=True)
        llm = self.root / "llm_run_manifest.json"
        original = json.loads(llm.read_text())
        for key, value in (("provider", "other"), ("model", "other"),
                           ("bundle_sha256", "0" * 64), ("prompt_sha256", "0" * 64),
                           ("response_sha256", "0" * 64), ("execution_evidence", "model_executed")):
            with self.subTest(key=key):
                llm.write_text(json.dumps({**original, key: value}))
                with self.assertRaisesRegex(ValueError, "LLM.*mismatch"):
                    self.prepare(inputs)
        llm.write_text(json.dumps(original))

    def test_proposal_provenance_and_claim_changes_are_rejected_even_with_updated_source_hash(self):
        inputs = self.inputs(with_llm_manifest=True)
        proposals = Path(inputs["hypotheses"])
        original = pd.read_csv(proposals, sep="\t", dtype=str).fillna("")
        source_manifest = Path(inputs["source_verification"]).with_name("semantic_claim_source_verification_manifest.json")
        source = json.loads(source_manifest.read_text())
        for key, value in (("provider", "other"), ("model", "other"), ("mechanism", "altered claim")):
            with self.subTest(key=key):
                frame = original.copy()
                frame.loc[0, key] = value
                frame.to_csv(proposals, sep="\t", index=False)
                source_manifest.write_text(json.dumps({**source, "proposals_sha256": sha256(proposals)}))
                with self.assertRaisesRegex(ValueError, "frozen response/roster mismatch"):
                    self.prepare(inputs)

    def test_legacy_missing_status_remains_explicit_unknown(self):
        inputs = self.inputs()
        Path(inputs["hypotheses"]).with_name("semantic_claim_proposal_status.tsv").unlink()
        _, _, manifest, observations = self.finish(inputs)
        self.assertIsNone(manifest["ancestors"]["proposal_status"])
        self.assertEqual("not_assessed_missing_import_status", observations[0]["status"])

    def test_reviewed_empty_claim_roster_remains_complete_without_candidates(self):
        _, result, _, observations = self.finish(self.inputs(no_claims=True))
        self.assertEqual("import_complete", observations[0]["status"])
        self.assertEqual(0, result["candidates"]["candidate_gene_count"])
        self.assertEqual(0, result["candidates"]["review_candidate_gene_count"])

    def test_status_tampering_before_and_after_snapshot_is_rejected(self):
        inputs = self.inputs()
        status = Path(inputs["hypotheses"]).with_name("semantic_claim_proposal_status.tsv")
        original = status.read_bytes()
        frame = pd.read_csv(status, sep="\t")
        frame.loc[0, "status"] = "no_claims"
        frame.to_csv(status, sep="\t", index=False)
        with self.assertRaisesRegex(ValueError, "frozen response/roster mismatch"):
            self.prepare(inputs)
        status.write_bytes(original)
        prepared = self.prepare(inputs)
        snapshot = self.root / "session/inputs/semantic_claim_proposal_status.tsv"
        snapshot.write_bytes(snapshot.read_bytes() + b"\n")
        with self.assertRaisesRegex(ValueError, "artifact mismatch"):
            complete_session(str(self.root / "session"), prepared["request_sha256"])

    def test_unrelated_output_paths_and_source_binding_are_rejected(self):
        inputs = self.inputs()
        import_manifest = Path(inputs["hypotheses"]).with_name("semantic_claim_proposal_manifest.json")
        original = import_manifest.read_bytes()
        value = json.loads(original)
        value["outputs"]["proposals"] = str(self.root / "unrelated.tsv")
        import_manifest.write_text(json.dumps(value))
        with self.assertRaisesRegex(ValueError, "not bound"):
            self.prepare(inputs)
        value = json.loads(original)
        value["case_id"] = "another-case"
        import_manifest.write_text(json.dumps(value))
        with self.assertRaisesRegex(ValueError, "case_id mismatch"):
            self.prepare(inputs)
        import_manifest.unlink()
        with self.assertRaisesRegex(ValueError, "requires import and source manifests"):
            self.prepare(inputs)
        import_manifest.write_bytes(original)
        bundle = Path(inputs["bundle"])
        bundle.write_bytes(bundle.read_bytes() + b"\n")
        with self.assertRaisesRegex(ValueError, "source binding hash mismatch"):
            self.prepare(inputs)

    def test_duplicate_missing_and_unknown_roster_units_are_rejected(self):
        inputs = self.inputs()
        status = Path(inputs["hypotheses"]).with_name("semantic_claim_proposal_status.tsv")
        frame = pd.read_csv(status, sep="\t").fillna("")
        variants = [pd.concat([frame, frame]), frame.iloc[:0], frame.assign(unit_id="unrelated")]
        for variant in variants:
            with self.subTest(rows=len(variant)):
                variant.to_csv(status, sep="\t", index=False)
                with self.assertRaises(ValueError):
                    self.prepare(inputs)


if __name__ == "__main__":
    unittest.main()
