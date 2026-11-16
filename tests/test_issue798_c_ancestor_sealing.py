"""Bound C ancestor relocation through the canonical A/B/C handoff.

All evidence is synthetic; these tests do not measure biological performance.
"""
import json
import shutil
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.claim_review_session import complete_session, handoff_sessions, prepare_session, sha256
from candidate_qtl.route_c_direct_members import build_route_c_direct_members
from candidate_qtl.route_c_direct_pipeline import direct_rescue_project
from tests.test_agent_claim_review import ReviewFixture, write_json
from tests import test_route_b_agent_claim_review as b_tests
from tests import test_issue_663_route_c_projection_evidence as c_tests


@unittest.skipUnless(shutil.which("blastp") and shutil.which("makeblastdb"), "BLAST+ required")
class BoundCAncestorSealingTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        ar = self.root / "a-inputs"
        ar.mkdir()
        f = ReviewFixture(ar)
        f.record["sequence"]["value"] = "M" + "ACDEFGHIKLMNPQRSTVWY" * 8
        write_json(f.records_path, f.record)
        self.a = self.root / "a"
        prepared = prepare_session(str(self.a), route="a", inputs={
            "hypotheses": str(f.hypotheses), "sources": str(f.sources_path),
            "proposal": str(f.proposal_raw), "records": str(f.records_path)},
            proposer_id="fixture-a-proposer", proposer_model="fixture",
            focal_species="Pisum sativum", case_id="case403")
        complete_session(str(self.a), prepared["request_sha256"],
            response_arg=str(f.raw_path), reviewer_id="fixture-a-reviewer", reviewer_model="fixture")
        harness = b_tests.RouteBAgentClaimReviewTests()
        harness.setUp()
        self.addCleanup(harness.doCleanups)
        br = Path(harness.tmp.name)
        proposals, source, bundle = harness._cli_inputs(br)
        self.b = self.root / "b"
        prepared = prepare_session(str(self.b), route="b", inputs={
            "hypotheses": str(proposals), "sources": str(harness.sources),
            "proposal": str(harness.proposal), "records": str(br / "records.json"),
            "bundle": str(bundle), "source_verification": str(source)},
            proposer_id="fixture-b-proposer", proposer_model="fixture", case_id="case403")
        complete_session(str(self.b), prepared["request_sha256"],
            response_arg=str(harness.raw), reviewer_id="fixture-b-reviewer", reviewer_model="fixture")
        self.genes = self.root / "genes.tsv"
        pd.DataFrame([{"gene_id": gene, "family_id": "EV1", "seqid": "synthetic",
            "start": i * 100 + 1, "end": i * 100 + 100}
            for i, gene in enumerate(("G1", "G2", "Psat.A"))]).to_csv(self.genes, sep="\t", index=False)
        self.proteins = self.root / "proteins.fa"
        self.proteins.write_text(">G1\n" + f.record["sequence"]["value"] + "\n>G2\nM" + "V" * 160 + "\n>Psat.A\nM" + "W" * 160 + "\n")
        cr = self.root / "c-inputs"
        cr.mkdir()
        frozen, self.freeze = c_tests.RouteCProjectionEvidenceTests()._frozen(cr, include_binding=True)
        payload = json.loads(self.freeze.read_text())
        payload["case_id"] = "case403"
        write_json(self.freeze, payload)
        project = direct_rescue_project(str(frozen), str(self.freeze), str(self.genes),
            "q1", str(cr / "projection"))
        self.c = build_route_c_direct_members(str(frozen), str(self.freeze),
            project["projection"], project["manifest"], str(self.genes), str(cr / "members"), case_id="case403")
        self.manifest = Path(self.c["manifest"])
        self.original = {p: p.read_bytes() for p in cr.rglob("*") if p.is_file()}

    def run_handoff(self, name="handoff"):
        return handoff_sessions(str(self.a), str(self.b), str(self.root / name),
            genes_arg=str(self.genes), proteins_arg=str(self.proteins), atlas_arg=str(self.genes),
            route_c_members_arg=self.c["members"], route_c_manifest_arg=str(self.manifest),
            allow_unrecorded_route_b=True)

    def test_fresh_bound_handoff_creates_relocated_manifest_once(self):
        result = self.run_handoff()
        self.assertEqual("first_stage_frozen", result["state"])
        final = pd.read_csv(result["wet"]["wet_candidates"], sep="\t")
        self.assertIn("Psat.A", set(final.candidate_gene_id))
        self.assertTrue(final.loc[final.candidate_gene_id.eq("Psat.A"), "route_c_supported"].iloc[0])
        sealed = json.loads(Path(result["execution_context"]["route_c_execution"]["manifest_path"]).read_text())
        spec = sealed["inputs"]["freeze_manifest"]
        rebound = json.loads(Path(spec["path"]).read_text())
        self.assertEqual(sha256(spec["path"]), spec["sha256"])
        self.assertNotEqual(sha256(self.freeze), spec["sha256"])
        self.assertEqual(sealed["claim_binding"]["artifact"], rebound["claim_binding"]["artifact"])
        self.assertEqual(1, len(list((self.root / "handoff").glob("route_c_claim_binding_*.tsv"))))
        for path, content in self.original.items():
            self.assertEqual(content, path.read_bytes(), str(path))
        with self.assertRaises(FileExistsError):
            self.run_handoff()

    def test_repeated_binding_descriptor_must_match_its_own_digest(self):
        frozen = json.loads(self.freeze.read_text())
        frozen["claim_binding"]["sha256"] = "0" * 64
        write_json(self.freeze, frozen)
        manifest = json.loads(self.manifest.read_text())
        manifest["inputs"]["freeze_manifest"]["sha256"] = sha256(self.freeze)
        write_json(self.manifest, manifest)
        with self.assertRaisesRegex(ValueError, "claim-binding artifact hash mismatch"):
            self.run_handoff()

    def test_ancestor_only_binding_is_relocated_without_top_descriptor(self):
        payload = json.loads(self.manifest.read_text())
        payload.pop("claim_binding")
        write_json(self.manifest, payload)
        result = self.run_handoff()
        sealed = json.loads(Path(result["execution_context"]["route_c_execution"]["manifest_path"]).read_text())
        rebound = json.loads(Path(sealed["inputs"]["freeze_manifest"]["path"]).read_text())
        self.assertTrue(Path(rebound["claim_binding"]["artifact"]).is_relative_to(self.root / "handoff"))
        self.assertEqual(sha256(rebound["claim_binding"]["artifact"]), rebound["claim_binding"]["sha256"])
        self.assertEqual("first_stage_frozen", result["state"])

    def test_rebound_sidecar_roster_mismatch_still_rejected_by_consumer(self):
        payload = json.loads(self.manifest.read_text())
        binding = Path(payload["claim_binding"]["artifact"])
        frame = pd.read_csv(binding, sep="\t", dtype=str).fillna("")
        frame.loc[0, "record_id"] = "different-record"
        frame.to_csv(binding, sep="\t", index=False)
        payload["claim_binding"]["sha256"] = sha256(binding)
        frozen = json.loads(self.freeze.read_text())
        frozen["claim_binding"]["sha256"] = sha256(binding)
        write_json(self.freeze, frozen)
        payload["inputs"]["freeze_manifest"]["sha256"] = sha256(self.freeze)
        write_json(self.manifest, payload)
        with self.assertRaisesRegex(ValueError, "record|roster"):
            self.run_handoff()

    def test_corrupted_binding_or_ancestor_bytes_fail_closed(self):
        payload = json.loads(self.manifest.read_text())
        binding = Path(payload["claim_binding"]["artifact"])
        binding.write_text(binding.read_text() + "\n")
        with self.assertRaisesRegex(ValueError, "claim-binding artifact hash mismatch"):
            self.run_handoff("bad-binding")
        binding.write_bytes(self.original[binding])
        self.freeze.write_text(self.freeze.read_text() + "\n")
        with self.assertRaisesRegex(ValueError, "producer ancestor hash mismatch"):
            self.run_handoff("bad-ancestor")


if __name__ == "__main__":
    unittest.main()
