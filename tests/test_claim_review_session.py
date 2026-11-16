"""Session lifecycle tests use explicitly synthetic evidence, never performance data."""
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

import pandas as pd

from candidate_qtl.claim_review_session import (
    _files,
    _validate_route_b_execution,
    _validate_route_c_execution,
    complete_session,
    finish_wet_review,
    handoff_sessions,
    prepare_session,
    sha256,
)
from tests.test_agent_claim_review import ReviewFixture, write_json
from tests import test_route_b_agent_claim_review as b_tests
from scripts.agentic_no_source_review import build as build_no_source_review


class ClaimReviewSessionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.fixture = ReviewFixture(self.root)
        self.session = self.root / "session"
        f = self.fixture
        f.record["sequence"]["value"] = "M" + "ACDEFGHIKLMNPQRSTVWY" * 8
        write_json(f.records_path, f.record)
        self.prepared = prepare_session(str(self.session), route="a", inputs={
            "hypotheses": str(f.hypotheses), "sources": str(f.sources_path),
            "proposal": str(f.proposal_raw), "records": str(f.records_path)},
            proposer_id="fixture-proposer", proposer_model="fixture", focal_species="Pisum sativum", case_id="case403")

    def complete(self, **kwargs):
        return complete_session(str(self.session), self.prepared["request_sha256"],
                                response_arg=str(self.fixture.raw_path),
                                reviewer_id="fixture-reviewer", reviewer_model="fixture", **kwargs)

    def test_a_consumer_and_idempotent_resume_preserve_original_source(self):
        self.fixture.sources_path.write_text("{}")
        result = self.complete()
        self.assertEqual("completed", result["state"])
        self.assertEqual(1, result["verification"]["verified_count"])
        self.assertTrue(Path(result["candidates"]["panel"]).exists())
        self.assertEqual(result, self.complete())

    def test_strict_reference_policy_survives_session_snapshot(self):
        f = self.fixture
        policy = f.hypotheses.parent / "mechanism_import_manifest.json"
        write_json(policy, {"benchmark_mode": True})
        session = self.root / "strict-session"
        prepared = prepare_session(str(session), route="a", inputs={
            "hypotheses": str(f.hypotheses), "sources": str(f.sources_path),
            "proposal": str(f.proposal_raw), "records": str(f.records_path)},
            proposer_id="fixture-proposer", proposer_model="fixture",
            focal_species="Arabidopsis thaliana", case_id="case403")
        policy.unlink()
        request = json.loads((session / "request.json").read_text())
        self.assertIn("inputs/mechanism_import_manifest.json",
                      {name.replace("\\", "/") for name in request["files"]})
        result = complete_session(str(session), prepared["request_sha256"],
                                  response_arg=str(f.raw_path), reviewer_id="fixture-reviewer",
                                  reviewer_model="fixture")
        validation = pd.read_csv(result["verification"]["validation"], sep="\t").fillna("")
        self.assertEqual("reference_gene_is_from_the_focal_species", validation.iloc[0]["rejection_reason"])
        self.assertEqual("not_verified", validation.iloc[0]["identity_sequence_verification_status"])

    def test_partial_review_is_not_complete_and_never_promotes(self):
        write_json(self.fixture.raw_path, {"reviews": []})
        result = self.complete()
        self.assertEqual("completed_partial_review", result["state"])
        self.assertEqual(1, result["review_counts"]["missing"])
        self.assertEqual(0, result["verification"]["verified_count"])

    def test_no_source_machine_state_is_partial_not_independent_review(self):
        f = self.fixture
        source = json.loads(f.sources_path.read_text(encoding="utf-8"))
        source["records"] = {}
        write_json(f.sources_path, source)
        session = self.root / "no-source-session"
        prepared = prepare_session(str(session), route="a", inputs={
            "hypotheses": str(f.hypotheses), "sources": str(f.sources_path),
            "proposal": str(f.proposal_raw), "records": str(f.records_path)},
            proposer_id="fixture-proposer", proposer_model="fixture",
            focal_species="Pisum sativum", case_id="case403")
        response = self.root / "no-source.json"
        write_json(response, build_no_source_review(Path(prepared["prompt"]), "a"))
        result = complete_session(str(session), prepared["request_sha256"],
                                  response_arg=str(response), reviewer_id="no-source-machine-a",
                                  reviewer_model="not_invoked")
        self.assertEqual("completed_partial_review", result["state"])
        self.assertEqual("not_executed_no_source", result["review_execution_state"])
        self.assertEqual({"requested": 1, "returned": 0, "missing": 1}, result["review_counts"])
        self.assertEqual(0, result["verification"]["verified_count"])

    def test_execution_context_validation_is_fail_closed(self):
        prompts = self.root / "route-b-prompts.jsonl"
        prompts.write_text(json.dumps({"prompt_id": "p1", "route_b_unit_id": "u1", "case_id": "case403"}) + "\n", encoding="utf-8")
        prompt_row = {"prompt_id": "p1", "route_b_unit_id": "u1", "case_id": "case403"}
        bundle = self.root / "route-b-bundle.json"
        write_json(bundle, {
            "units": [{"unit_id": "u1"}],
            "semantic_unit_count": 1,
            "compression_observability": {
                "semantic_unit_count": 1,
                "requested_unit_count": 1,
                "skipped_unit_count": 0,
            },
            "execution_prompt_count": 1,
            "execution_prompt_fingerprints": [hashlib.sha256(
                json.dumps(prompt_row, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
            ).hexdigest()],
        })
        response = self.root / "route-b-response.json"
        write_json(response, {
            "case_id": "case403",
            "candidates": [{
                "prompt_id": "p1", "route_b_unit_id": "u1", "parse_status": "json",
                "raw_text": '{"claims":[]}', "model_json": {"claims": []},
            }],
        })
        attempt = response
        manifest = self.root / "route-b-execution.json"
        payload = {
            "contract": "agentic-route-b-execution-1.1.0",
            "state": "complete",
            "case_id": "case403",
            "source_prompts": str(prompts),
            "prompt_sha256": sha256(prompts),
            "output": str(response),
            "response_sha256": sha256(response),
            "prompt_count": 1,
            "completed_count": 1,
            "not_evaluated_count": 0,
            "response_roster": [{
                "prompt_id": "p1", "route_b_unit_id": "u1", "parse_status": "json",
                "raw_text_sha256": hashlib.sha256(b'{"claims":[]}').hexdigest(),
            }],
            "partial_response_provenance": [{
                "attempt_index": 0, "file": str(attempt), "sha256": sha256(attempt),
                "executor": "fixture-b", "model": "fixture",
            }],
            "retry_history": [],
        }
        write_json(manifest, payload)
        route_b = {
            **payload,
            "artifact_path": str(manifest),
            "artifact_sha256": sha256(manifest),
        }
        validated = _validate_route_b_execution(route_b, "case403", bundle_path=bundle)
        self.assertEqual("complete", validated["state"])
        broken_bundle = json.loads(bundle.read_text(encoding="utf-8"))
        broken_bundle["compression_observability"]["requested_unit_count"] = 0
        write_json(bundle, broken_bundle)
        with self.assertRaisesRegex(ValueError, "semantic/prompt accounting mismatch"):
            _validate_route_b_execution(route_b, "case403", bundle_path=bundle)
        write_json(bundle, {
            "units": [{"unit_id": "u1"}],
            "semantic_unit_count": 1,
            "compression_observability": {
                "semantic_unit_count": 1,
                "requested_unit_count": 1,
                "skipped_unit_count": 0,
            },
            "execution_prompt_count": 1,
            "execution_prompt_fingerprints": [hashlib.sha256(
                json.dumps(prompt_row, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
            ).hexdigest()],
        })
        with self.assertRaisesRegex(ValueError, "counts do not cover"):
            _validate_route_b_execution({**route_b, "completed_count": 999}, "case403")
        with self.assertRaisesRegex(ValueError, "case_id mismatch"):
            _validate_route_b_execution(route_b, "other-case")

        status = self.root / "route-c-status.json"
        write_json(status, {"case_id": "case403", "state": "provider_failure", "reason": "fixture unavailable"})
        route_c = _validate_route_c_execution({
            **json.loads(status.read_text()),
            "artifact_path": str(status), "artifact_sha256": sha256(status),
        }, "case403", None)
        self.assertEqual("provider_failure", route_c["state"])
        with self.assertRaisesRegex(ValueError, "requires a reason"):
            _validate_route_c_execution({"case_id": "case403", "state": "provider_failure"}, "case403", None)
        with self.assertRaisesRegex(ValueError, "requires a member artifact"):
            _validate_route_c_execution({
                "case_id": "case403", "state": "completed_members",
                "artifact_path": str(status), "artifact_sha256": sha256(status),
            }, "case403", None)
        with self.assertRaisesRegex(ValueError, "non-completed execution state"):
            _validate_route_c_execution({
                "case_id": "case403", "state": "provider_failure", "reason": "fixture unavailable",
                "artifact_path": str(status), "artifact_sha256": sha256(status),
            }, "case403", str(response))

        from tests._issue664_helpers import manifest_for
        empty_members = self.root / "empty-route-c-members.tsv"
        pd.DataFrame(columns=["candidate_gene_id", "support_status"]).to_csv(
            empty_members, sep="\t", index=False
        )
        empty_manifest = manifest_for(empty_members, route="c", case_id="case403")
        with self.assertRaisesRegex(ValueError, "explicit completed_empty status"):
            _validate_route_c_execution(
                None, "case403", str(empty_members), str(empty_manifest)
            )
        empty_status = self.root / "empty-route-c-status.json"
        write_json(empty_status, {
            "case_id": "case403", "state": "completed_empty",
            "reason": "fixture says search completed",
        })
        with self.assertRaisesRegex(ValueError, "sealed search-generation lineage"):
            _validate_route_c_execution({
                **json.loads(empty_status.read_text()),
                "status_path": str(empty_status),
                "status_sha256": sha256(empty_status),
            }, "case403", str(empty_members), str(empty_manifest))

    def test_input_tampering_and_same_executor_rejected(self):
        with self.assertRaisesRegex(ValueError, "separate actual"):
            complete_session(str(self.session), self.prepared["request_sha256"],
                             response_arg=str(self.fixture.raw_path), reviewer_id="fixture-proposer", reviewer_model="fixture")
        request = json.loads((self.session / "request.json").read_text())
        (self.session / request["inputs"]["sources"]).write_text("{}")
        with self.assertRaisesRegex(ValueError, "artifact mismatch"):
            self.complete()

    def test_completed_outputs_cannot_be_tampered_or_replaced(self):
        result = self.complete()
        Path(result["candidates"]["panel"]).write_text("tampered")
        with self.assertRaisesRegex(ValueError, "artifact mismatch"):
            self.complete()

    def test_failed_attempt_is_preserved_and_not_silently_retried(self):
        self.fixture.raw_path.write_text("not JSON")
        with self.assertRaises(ValueError):
            self.complete()
        failure = self.session / "completed" / "failure.json"
        self.assertEqual("failed", json.loads(failure.read_text())["state"])
        before = failure.read_bytes()
        with self.assertRaises(FileExistsError):
            self.complete()
        self.assertEqual(before, failure.read_bytes())
        self.assertFalse((self.session / "completion.json").exists())

    def test_a_and_b_support_both_reach_wet(self):
        harness = b_tests.RouteBAgentClaimReviewTests()
        harness.setUp()
        self.addCleanup(harness.doCleanups)
        root = Path(harness.tmp.name)
        proposals, source, bundle = harness._cli_inputs(root)
        prepared = prepare_session(str(root / "session"), route="b", inputs={
            "hypotheses": str(proposals), "sources": str(harness.sources),
            "proposal": str(harness.proposal), "records": str(root / "records.json"),
            "bundle": str(bundle), "source_verification": str(source)},
            proposer_id="fixture-b-proposer", proposer_model="fixture", case_id="case403")
        result = complete_session(str(root / "session"), prepared["request_sha256"],
                                  response_arg=str(harness.raw), reviewer_id="fixture-b-reviewer", reviewer_model="fixture")
        self.assertEqual("completed", result["state"])
        self.assertEqual(1, result["candidates"]["candidate_gene_count"])
        if shutil.which("blastp") and shutil.which("makeblastdb"):
            self.complete()
            genes = root / "genes.tsv"
            pd.DataFrame([{"gene_id": gene, "family_id": "EV1", "seqid": "synthetic",
                           "start": i * 100 + 1, "end": i * 100 + 100} for i, gene in enumerate(("G1", "G2"))]).to_csv(genes, sep="\t", index=False)
            proteins = root / "proteins.fa"
            proteins.write_text(">G1\n" + self.fixture.record["sequence"]["value"] + "\n>G2\n" + "M" + "V" * 160 + "\n")
            from tests._issue664_helpers import manifest_for
            route_c = root / "route_c_members.tsv"
            pd.DataFrame([
                {"candidate_gene_id": "G1", "support_status": "verified_direct_support",
                 "record_id": "RC1", "source_id": "PMID:1", "species_distance": "S1_same_species",
                 "evidence_directness": "D1_mapped_gene_with_functional_validation",
                 "projection_state": "inside_current_qtl", "projection_basis": "exact_focal_id"},
                {"candidate_gene_id": "G2", "support_status": "reference_phenotype_review",
                 "record_id": "RC2", "source_id": "PMID:2", "species_distance": "S5_model_or_other_plant",
                 "evidence_directness": "D0_causal_perturbation_or_rescue",
                 "projection_state": "inside_current_qtl", "projection_basis": "orthology_reference_gene"},
            ]).to_csv(route_c, sep="\t", index=False)
            route_c_manifest = manifest_for(route_c, route="c", case_id="case403")
            script = Path(__file__).parents[1] / "scripts" / "claim_review_session.py"
            forged_b = root / "forged-route-b.json"
            write_json(forged_b, {"state": "complete", "completed_count": 999})
            base_command = [
                sys.executable, str(script), "handoff",
                "--a-session", str(self.session), "--b-session", str(root / "session"),
                "--out", str(root / "forbidden-handoff"), "--genes", str(genes),
                "--proteins", str(proteins), "--atlas", str(genes),
                "--allow-unrecorded-route-b",
            ]
            rejected_b = subprocess.run(
                [*base_command, "--route-b-execution-manifest", str(forged_b)],
                cwd=Path(__file__).parents[1], capture_output=True, text=True,
            )
            self.assertNotEqual(0, rejected_b.returncode)
            self.assertIn("unsupported contract", rejected_b.stderr)
            forged_c = root / "forged-route-c.json"
            write_json(forged_c, {"case_id": "case403", "state": "completed_members"})
            rejected_c = subprocess.run(
                [*base_command, "--out", str(root / "forbidden-c-handoff"), "--route-c-status", str(forged_c)],
                cwd=Path(__file__).parents[1], capture_output=True, text=True,
            )
            self.assertNotEqual(0, rejected_c.returncode)
            self.assertIn("requires a member artifact", rejected_c.stderr)
            wet = handoff_sessions(str(self.session), str(root / "session"), str(root / "handoff"),
                genes_arg=str(genes), proteins_arg=str(proteins), atlas_arg=str(genes),
                route_c_members_arg=str(route_c), route_c_manifest_arg=str(route_c_manifest),
                allow_unrecorded_route_b=True)
            self.assertEqual("first_stage_frozen", wet["state"])
            bridge = pd.read_csv(wet["a"]["bridge_tsv"], sep="\t")
            self.assertIn("G1", bridge.loc[bridge.route.eq("a"), "gene_id"].tolist())
            self.assertEqual(2, wet["family"]["family_count"])
            final = pd.read_csv(wet["wet"]["wet_candidates"], sep="\t")
            self.assertEqual(["G1"], final["candidate_gene_id"].tolist())
            g1 = final.loc[final.candidate_gene_id.eq("G1")].iloc[0]
            self.assertTrue(g1["route_a_supported"])
            self.assertTrue(g1["route_b_supported"])
            self.assertTrue(g1["route_c_supported"])
            union_manifest = json.loads(Path(wet["union"]["manifest"]).read_text())
            self.assertEqual(1, union_manifest["route_c_verified_support_rows"])
            route_c_input = next(
                row for row in union_manifest["inputs"] if row.get("route") == "c"
            )
            sealed_c = wet["execution_context"]["route_c_execution"]
            self.assertEqual(
                Path(route_c_input["artifact"]).resolve(),
                Path(sealed_c["members_path"]).resolve(),
            )
            self.assertEqual(
                Path(route_c_input["producer_manifest"]).resolve(),
                Path(sealed_c["manifest_path"]).resolve(),
            )
            final = pd.read_csv(wet["wet"]["wet_candidates"], sep="\t", dtype=str).fillna("")
            wet_manifest = json.loads(Path(wet["wet"]["wet_candidates_manifest"]).read_text())
            trace = pd.read_csv(wet_manifest["candidate_trace"], sep="\t", dtype=str).fillna("")
            g2 = trace.set_index("candidate_gene_id").loc["G2"]
            # The hand-written C row above has no sealed projection/binding
            # ancestry and is not an eligible review seed at the Wet consumer.
            self.assertEqual("not_selected_no_upstream_candidate", g2.selection_disposition)

    def test_handoff_relocates_bound_c_freeze_without_rewriting_sources(self):
        from candidate_qtl.route_c_direct_members import build_route_c_direct_members
        from candidate_qtl.route_c_review_candidates import read_direct_c_review_candidates
        from tests.test_issue_663_route_c_projection_evidence import RouteCProjectionEvidenceTests

        harness = b_tests.RouteBAgentClaimReviewTests()
        harness.setUp()
        self.addCleanup(harness.doCleanups)
        root = Path(harness.tmp.name)
        proposals, source, bundle = harness._cli_inputs(root)
        b_session = root / "b-session"
        prepared = prepare_session(str(b_session), route="b", inputs={
            "hypotheses": str(proposals), "sources": str(harness.sources),
            "proposal": str(harness.proposal), "records": str(root / "records.json"),
            "bundle": str(bundle), "source_verification": str(source)},
            proposer_id="fixture-b-proposer", proposer_model="fixture", case_id="case403")
        complete_session(str(b_session), prepared["request_sha256"], response_arg=str(harness.raw),
                         reviewer_id="fixture-b-reviewer", reviewer_model="fixture")
        self.complete()

        c_root = root / "c-source"
        c_root.mkdir()
        fixture = RouteCProjectionEvidenceTests()
        original_binding = fixture._binding_row

        def review_only_binding(record):
            return {**original_binding(record), "binding_status": "review_required",
                    "perturbation_to_phenotype_binding_state": "not_evaluated"}

        fixture._binding_row = review_only_binding
        projected, frozen, freeze, atlas = fixture._project(
            c_root, [{"gene_id": "Psat.A", "chrom": "chr1", "start": "1", "end": "100"}],
            None, include_binding=True)
        freeze_payload = json.loads(freeze.read_text(encoding="utf-8"))
        freeze_payload["case_id"] = "case403"
        freeze_payload["claim_binding"]["complete"] = False
        write_json(freeze, freeze_payload)
        members = build_route_c_direct_members(
            str(frozen), str(freeze), projected["projection"], projected["manifest"],
            str(atlas), str(c_root / "members"), case_id="case403")
        self.assertEqual(0, members["supported_count"])
        originals = {path: path.read_bytes() for path in c_root.rglob("*") if path.is_file()}
        copied_root = root / "handoff-bound-c"
        proteins = root / "proteins.fa"
        proteins.write_text(">Psat.A\n" + self.fixture.record["sequence"]["value"] + "\n", encoding="utf-8")
        kwargs = dict(genes_arg=str(atlas), proteins_arg=str(proteins),
                      atlas_arg=str(atlas), route_c_members_arg=members["members"],
                      route_c_manifest_arg=members["manifest"], allow_unrecorded_route_b=True)
        # Stop after real handoff sealing, before native BLAST; the C consumer
        # below reads the produced ancestry and binding without mocking gates.
        with patch("candidate_qtl.blast_evidence.run_panel_blastp", side_effect=RuntimeError("sealing checkpoint")):
            with self.assertRaisesRegex(RuntimeError, "sealing checkpoint"):
                handoff_sessions(str(self.session), str(b_session), str(copied_root), **kwargs)
        copied_manifest_path = copied_root / "route_c_manifest.json"
        copied_manifest = json.loads(copied_manifest_path.read_text(encoding="utf-8"))
        freeze_spec = copied_manifest["inputs"]["freeze_manifest"]
        copied_freeze_path = Path(freeze_spec["path"])
        copied_freeze = json.loads(copied_freeze_path.read_text(encoding="utf-8"))
        binding_path = Path(copied_freeze["claim_binding"]["artifact"])
        self.assertEqual(copied_root.resolve(), binding_path.parent.resolve())
        self.assertEqual(sha256(binding_path), copied_freeze["claim_binding"]["sha256"])
        self.assertEqual(sha256(copied_freeze_path), freeze_spec["sha256"])
        self.assertNotEqual(sha256(freeze), freeze_spec["sha256"])
        self.assertEqual(copied_manifest["claim_binding"]["artifact"], str(binding_path))
        self.assertEqual(1, len(list(copied_root.glob("route_c_claim_binding_*.tsv"))))
        candidates, receipt = read_direct_c_review_candidates(
            str(copied_root / "route_c_members.tsv"), str(copied_manifest_path),
            case_id="case403", atlas_frame=pd.read_csv(atlas, sep="\t", dtype=str), supported_ids=set())
        self.assertEqual(1, receipt["review_record_count"])
        self.assertEqual("Psat.A", candidates.iloc[0]["candidate_gene_id"])
        self.assertEqual("unverified_route_review", candidates.iloc[0]["support_status"])
        copied_bytes = {path: path.read_bytes() for path in copied_root.rglob("*") if path.is_file()}
        with self.assertRaises(FileExistsError):
            handoff_sessions(str(self.session), str(b_session), str(copied_root), **kwargs)
        self.assertEqual(copied_bytes, {path: path.read_bytes() for path in copied_root.rglob("*") if path.is_file()})
        for path, contents in originals.items():
            self.assertEqual(contents, path.read_bytes(), str(path))

    def test_second_stage_import_preserves_members_without_verified_exclusions(self):
        from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
        from tests._issue664_helpers import manifest_for
        source = self.root / "supported.tsv"
        pd.DataFrame({"candidate_gene_id": [f"G{i}" for i in range(26)], "route": "b",
                      "support_status": "verified_semantic_support"}).to_csv(source, sep="\t", index=False)
        handoff = self.root / "first_stage"
        wet = wet_candidate_handoff(str(handoff / "wet"), route_b_candidates_arg=str(source),
            route_b_candidates_manifest_arg=str(manifest_for(source, route="b", case_id="case403")), case_id="case403")
        write_json(handoff / "completion.json", {"case_id": "case403", "wet": wet, "files": _files(handoff)})
        proposals = self.root / "proposals.tsv"
        pd.DataFrame([{"ambiguity_kind": "family_vs_subfamily", "ambiguity_target": "F1", "affected_candidate_genes": "G0",
                       "evidence_kind": "sequence_check", "discriminating_evidence_request": "Compare sequences"}]).to_csv(proposals, sep="\t", index=False)
        before = Path(wet["wet_candidates"]).read_bytes()
        result = finish_wet_review(str(handoff), sha256(handoff / "completion.json"), str(self.root / "second_stage"),
            proposals_arg=str(proposals), reviewer_id="fixture-reviewer", reviewer_model="fixture")
        self.assertEqual("completed_review_only", result["state"])
        self.assertEqual(26, result["refinement"]["output_candidate_count"])
        self.assertEqual(before, Path(wet["wet_candidates"]).read_bytes())

    def test_review_inclusive_finish_requires_bound_independent_review_and_preserves_freeze(self):
        from candidate_qtl.adaptive_synthesis import _sha256_path
        from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
        from tests._issue664_helpers import manifest_for

        source = self.root / "supported-selection.tsv"
        pd.DataFrame({"candidate_gene_id": [f"G{i}" for i in range(26)],
                      "support_status": "verified_semantic_support",
                      "route": "b"}).to_csv(source, sep="\t", index=False)
        handoff = self.root / "selection-first-stage"
        wet = wet_candidate_handoff(str(handoff / "wet"), route_b_candidates_arg=str(source),
            route_b_candidates_manifest_arg=str(manifest_for(source, route="b", case_id="case-selection")), case_id="case-selection")
        write_json(handoff / "completion.json", {"case_id": "case-selection", "wet": wet, "files": _files(handoff)})
        completion_hash = sha256(handoff / "completion.json")
        frozen_before = Path(wet["wet_candidates"]).read_bytes()
        wet_hash = _sha256_path(Path(wet["wet_candidates"]))
        proposal = self.root / "selection-proposal.tsv"
        rows = [{"candidate_gene_id": f"G{i}", "decision": "defer" if i == 0 else ("not_evaluated" if i == 1 else "adopt"),
                 "reason": "independent ambiguity remains" if i == 0 else f"reviewed candidate G{i}",
                 "source_wet_sha256": wet_hash, "executor_id": "luna-executor", "executor_model": "gpt-5.6-luna"}
                for i in range(26)]
        pd.DataFrame(rows).to_csv(proposal, sep="\t", index=False)
        proposal_hash = _sha256_path(proposal)
        independent = self.root / "selection-independent.tsv"
        pd.DataFrame([{**row, "proposal_sha256": proposal_hash, "reviewer_id": "luna-reviewer",
                       "reviewer_model": "gpt-5.6-luna"} for row in rows]).to_csv(independent, sep="\t", index=False)
        result = finish_wet_review(str(handoff), completion_hash, str(self.root / "selection-out"),
            proposals_arg=str(proposal), reviewer_id="luna-reviewer", reviewer_model="gpt-5.6-luna",
            independent_review_arg=str(independent), executor_id="luna-executor", executor_model="gpt-5.6-luna")
        self.assertEqual("completed_review_inclusive_selection", result["state"])
        selected = pd.read_csv(result["selection"]["wet_candidates"], sep="\t", dtype=str)
        self.assertEqual(25, len(selected))
        self.assertNotIn("G0", set(selected.candidate_gene_id))
        self.assertIn("G1", set(selected.candidate_gene_id))
        self.assertEqual("supported", selected.loc[selected.candidate_gene_id.eq("G1"), "support_status"].iloc[0])
        ledger = pd.read_csv(result["selection"]["deferral_ledger"], sep="\t", dtype=str)
        self.assertEqual(["G0"], ledger.candidate_gene_id.tolist())
        self.assertEqual(frozen_before, Path(wet["wet_candidates"]).read_bytes())
        verification = self.root / "selection-verification.tsv"
        verification.write_text("candidate_gene_id\tverdict\tverification_source\nG0\tnot_evaluated\tlocal\n")
        with self.assertRaisesRegex(ValueError, "cannot be combined"):
            finish_wet_review(str(handoff), completion_hash, str(self.root / "selection-and-refinement"),
                proposals_arg=str(proposal), reviewer_id="luna-reviewer", reviewer_model="gpt-5.6-luna",
                independent_review_arg=str(independent), executor_id="luna-executor",
                verification_arg=str(verification))

    def test_review_selection_rejects_empty_or_mismatched_artifact_bindings(self):
        from candidate_qtl.adaptive_synthesis import consume_review_inclusive_selection, _sha256_path
        source = self.root / "selection-input.tsv"
        pd.DataFrame({"candidate_gene_id": ["G0", "G1"], "support_status": ["verified", "not_evaluated"]}).to_csv(source, sep="\t", index=False)
        wet_hash = _sha256_path(source)
        proposal = self.root / "selection-bad.tsv"
        pd.DataFrame({"candidate_gene_id": ["G0", "G1"], "decision": ["defer", "adopt"],
                      "reason": ["defer reason", "adopt reason"], "source_wet_sha256": [wet_hash, ""],
                      "executor_id": ["luna-executor"] * 2, "executor_model": ["gpt-5.6-luna"] * 2}).to_csv(proposal, sep="\t", index=False)
        with self.assertRaisesRegex(ValueError, "frozen Wet hash"):
            consume_review_inclusive_selection(str(source), str(proposal), str(self.root / "bad-out"),
                reviewer_id="luna-reviewer", reviewer_model="gpt-5.6-luna", executor_id="luna-executor", executor_model="gpt-5.6-luna")

    def test_unresolved_independent_review_keeps_all_candidates_and_rejects_empty_reason(self):
        from candidate_qtl.adaptive_synthesis import consume_review_inclusive_selection, _sha256_path
        source = self.root / "disagreement-input.tsv"
        pd.DataFrame({"candidate_gene_id": ["G0", "G1"], "support_status": ["supported", "not_evaluated"]}).to_csv(source, sep="\t", index=False)
        wet_hash = _sha256_path(source)
        proposal = self.root / "disagreement-proposal.tsv"
        pd.DataFrame({"candidate_gene_id": ["G0", "G1"], "decision": ["defer", "adopt"],
                      "reason": ["keep ambiguity", "retain pending review"], "source_wet_sha256": [wet_hash, wet_hash],
                      "executor_id": ["luna-executor"] * 2, "executor_model": ["gpt-5.6-luna"] * 2}).to_csv(proposal, sep="\t", index=False)
        independent = self.root / "disagreement-independent.tsv"
        pd.DataFrame({"candidate_gene_id": ["G0", "G1"], "decision": ["adopt", "adopt"],
                      "reason": ["disagree", "independent review"], "source_wet_sha256": [wet_hash, wet_hash],
                      "proposal_sha256": [_sha256_path(proposal)] * 2, "reviewer_id": ["luna-reviewer"] * 2,
                      "reviewer_model": ["gpt-5.6-luna"] * 2}).to_csv(independent, sep="\t", index=False)
        result = consume_review_inclusive_selection(str(source), str(proposal), str(self.root / "disagreement-out"),
            reviewer_id="luna-reviewer", reviewer_model="gpt-5.6-luna", executor_id="luna-executor",
            executor_model="gpt-5.6-luna",
            independent_review_arg=str(independent))
        self.assertEqual(2, result["output_candidate_count"])
        empty_reason = proposal
        frame = pd.read_csv(empty_reason, sep="\t")
        frame.loc[0, "reason"] = ""
        frame.to_csv(empty_reason, sep="\t", index=False)
        independent_frame = pd.read_csv(independent, sep="\t")
        independent_frame["proposal_sha256"] = _sha256_path(empty_reason)
        independent_frame.loc[0, "decision"] = "defer"
        independent_frame.to_csv(independent, sep="\t", index=False)
        with self.assertRaisesRegex(ValueError, "selection proposal requires non-empty reason"):
            consume_review_inclusive_selection(str(source), str(empty_reason), str(self.root / "empty-reason-out"),
                reviewer_id="luna-reviewer", reviewer_model="gpt-5.6-luna", executor_id="luna-executor",
                executor_model="gpt-5.6-luna",
                independent_review_arg=str(independent))

    def test_adopt_without_independent_review_is_not_evaluated_and_keeps_wet(self):
        from candidate_qtl.adaptive_synthesis import consume_review_inclusive_selection, _sha256_path
        source = self.root / "no-independent-input.tsv"
        pd.DataFrame({"candidate_gene_id": ["G0", "G1"], "support_status": ["supported", "not_evaluated"]}).to_csv(source, sep="\t", index=False)
        wet_hash = _sha256_path(source)
        proposal = self.root / "no-independent-proposal.tsv"
        pd.DataFrame({"candidate_gene_id": ["G0"], "decision": ["adopt"], "reason": ["retain pending independent review"],
                      "source_wet_sha256": [wet_hash], "executor_id": ["luna-executor"], "executor_model": ["gpt-5.6-luna"]}).to_csv(proposal, sep="\t", index=False)
        result = consume_review_inclusive_selection(str(source), str(proposal), str(self.root / "no-independent-out"),
            reviewer_id="luna-reviewer", reviewer_model="gpt-5.6-luna", executor_id="luna-executor", executor_model="gpt-5.6-luna")
        self.assertEqual("not_evaluated_no_independent_review", result["state"])
        self.assertEqual("partial", result["review_coverage"])
        self.assertEqual(2, result["output_candidate_count"])

    def test_selection_applies_only_agreed_defers_and_records_candidate_disagreements(self):
        from candidate_qtl.adaptive_synthesis import consume_review_inclusive_selection, _sha256_path

        source = self.root / "per-candidate-input.tsv"
        pd.DataFrame({"candidate_gene_id": ["G0", "G1", "G2"],
                      "support_status": ["supported", "not_evaluated", "review"]}).to_csv(source, sep="\t", index=False)
        wet_hash = _sha256_path(source)
        proposal = self.root / "per-candidate-proposal.tsv"
        pd.DataFrame({"candidate_gene_id": ["G0", "G1", "G2"],
                      "decision": ["defer", "defer", "retain"],
                      "reason": ["copy ambiguity", "specific unresolved mechanism", "retain candidate"],
                      "source_wet_sha256": [wet_hash] * 3,
                      "executor_id": ["luna-executor"] * 3,
                      "executor_model": ["gpt-5.6-luna"] * 3}).to_csv(proposal, sep="\t", index=False)
        independent = self.root / "per-candidate-independent.tsv"
        pd.DataFrame({"candidate_gene_id": ["G0", "G1", "G2"],
                      "decision": ["adopt", "deferred", "adopt"],
                      "reason": ["retain pending adjudication", "same unresolved mechanism", "retain pending review"],
                      "source_wet_sha256": [wet_hash] * 3,
                      "proposal_sha256": [_sha256_path(proposal)] * 3,
                      "reviewer_id": ["luna-reviewer"] * 3,
                      "reviewer_model": ["gpt-5.6-luna"] * 3}).to_csv(independent, sep="\t", index=False)
        result = consume_review_inclusive_selection(
            str(source), str(proposal), str(self.root / "per-candidate-out"),
            reviewer_id="luna-reviewer", reviewer_model="gpt-5.6-luna",
            executor_id="luna-executor", executor_model="gpt-5.6-luna",
            independent_review_arg=str(independent))
        selected = pd.read_csv(result["wet_candidates"], sep="\t", dtype=str)
        self.assertEqual(["G0", "G2"], selected.candidate_gene_id.tolist())
        ledger = pd.read_csv(result["deferral_ledger"], sep="\t", dtype=str)
        self.assertEqual(["G1"], ledger.candidate_gene_id.tolist())
        decisions = pd.read_csv(result["review_decisions"], sep="\t", dtype=str)
        self.assertEqual(["unresolved_disagreement", "agreement", "agreement"], decisions.resolution_state.tolist())
        self.assertEqual("adopt", decisions.loc[decisions.candidate_gene_id.eq("G2"), "applied_decision"].iloc[0])
        manifest = json.loads(Path(result["manifest"]).read_text())
        self.assertEqual(_sha256_path(Path(result["review_decisions"])), manifest["review_decisions_sha256"])
        self.assertEqual(1, result["disagreement_candidate_count"])

    def test_independent_unknown_selection_decision_is_rejected(self):
        from candidate_qtl.adaptive_synthesis import consume_review_inclusive_selection, _sha256_path

        source = self.root / "unknown-decision-input.tsv"
        pd.DataFrame({"candidate_gene_id": ["G0"], "support_status": ["review"]}).to_csv(source, sep="\t", index=False)
        wet_hash = _sha256_path(source)
        proposal = self.root / "unknown-decision-proposal.tsv"
        pd.DataFrame({"candidate_gene_id": ["G0"], "decision": ["adopt"], "reason": ["retain"],
                      "source_wet_sha256": [wet_hash], "executor_id": ["luna-executor"],
                      "executor_model": ["gpt-5.6-luna"]}).to_csv(proposal, sep="\t", index=False)
        independent = self.root / "unknown-decision-independent.tsv"
        pd.DataFrame({"candidate_gene_id": ["G0"], "decision": ["maybe"], "reason": ["unclear"],
                      "source_wet_sha256": [wet_hash], "proposal_sha256": [_sha256_path(proposal)],
                      "reviewer_id": ["luna-reviewer"], "reviewer_model": ["gpt-5.6-luna"]}).to_csv(independent, sep="\t", index=False)
        with self.assertRaisesRegex(ValueError, "unsupported independent review selection decision"):
            consume_review_inclusive_selection(
                str(source), str(proposal), str(self.root / "unknown-decision-out"),
                reviewer_id="luna-reviewer", reviewer_model="gpt-5.6-luna",
                executor_id="luna-executor", executor_model="gpt-5.6-luna",
                independent_review_arg=str(independent))

    def test_proposal_executor_must_match_and_defer_requires_independent(self):
        from candidate_qtl.adaptive_synthesis import consume_review_inclusive_selection, _sha256_path
        source = self.root / "executor-input.tsv"
        pd.DataFrame({"candidate_gene_id": ["G0"], "support_status": ["supported"]}).to_csv(source, sep="\t", index=False)
        wet_hash = _sha256_path(source)
        proposal = self.root / "executor-proposal.tsv"
        pd.DataFrame({"candidate_gene_id": ["G0"], "decision": ["defer"], "reason": ["specific unresolved mechanism"],
                      "source_wet_sha256": [wet_hash], "executor_id": ["other-executor"], "executor_model": ["gpt-5.6-luna"]}).to_csv(proposal, sep="\t", index=False)
        with self.assertRaisesRegex(ValueError, "executor_id does not match"):
            consume_review_inclusive_selection(str(source), str(proposal), str(self.root / "executor-out"),
                reviewer_id="luna-reviewer", reviewer_model="gpt-5.6-luna", executor_id="luna-executor", executor_model="gpt-5.6-luna")
