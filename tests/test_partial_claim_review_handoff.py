"""Handoff orchestration tests; biology/BLAST consumers are mocked explicitly."""
import json
import os
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

import pandas as pd

from candidate_qtl import claim_review_session as session


class PartialClaimReviewHandoffTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.out = self.root / "handoff"
        self.resources = self.root / "resources.tsv"
        self.resources.write_text("gene_id\nG1\n", encoding="utf-8")
        self.states = {}
        for route in ("a", "b"):
            root = self.root / route
            root.mkdir()
            data = root / "input.tsv"
            data.write_text("synthetic frozen input\n", encoding="utf-8")
            request = {"route": route, "case_id": "case1", "inputs": {"hypotheses": "input.tsv"}}
            session._write(root / "request.json", request)
            completed = {
                "state": "completed", "review_counts": {"requested": 2, "returned": 2, "missing": 0},
                "identity": {"request_sha256": session.sha256(root / "request.json")},
                "verification": {"validation": str(data), "reviewable_sequences": str(data)},
                "candidates": {key: str(data) for key in (
                    "panel", "projection_panel", "candidates", "review_candidates", "review_manifest")},
                "files": session._files(root),
            }
            self.states[route] = completed
            session._write(root / "completion.json", completed)
        self.consumers = {}
        modules = {
            "blast_panel": {"freeze_blast_panel": lambda *a, **k: {"manifest": "synthetic-freeze"},
                            "panel_from_verified_bait_panel": lambda *a, **k: pd.DataFrame({"panel_id": ["bait1"]})},
            "blast_evidence": {"run_panel_blastp": lambda *a, **k: {"evidence": "evidence", "eligible_hits": "hits"}},
            "canonical_family_recovery": {"build_canonical_family_recovery": lambda *a, **k: {
                "family_count": 0, "family_recovery": "family", "comparison": "comparison"}},
            "dual_route": {"phenotype_first_status": lambda *a, **k: {"hypotheses": "hypotheses"}},
            "family_candidate_bridge": {"family_candidate_bridge": lambda *a, **k: {
                "bridge_tsv": "bridge", "review_candidates": "a-review", "review_manifest": "a-review-manifest"}},
            "family_index": {"route_a_interval_family_match": lambda *a, **k: {"status": "interval-status"}},
            "family_to_copy": {"family_to_copy": lambda *a, **k: None},
            "supported_member_union": {"build_supported_member_union": self.union},
            "wet_candidate_handoff": {"wet_candidate_handoff": self.wet},
        }
        replacements = {}
        for name, functions in modules.items():
            module = types.ModuleType("candidate_qtl." + name)
            for function_name, implementation in functions.items():
                function = Mock(side_effect=implementation)
                setattr(module, function_name, function)
                self.consumers[function_name] = function
            replacements[module.__name__] = module
        dependency_patch = patch.dict(sys.modules, replacements)
        dependency_patch.start()
        self.addCleanup(dependency_patch.stop)
        self.wet_state = "wet_reviewable"

    def change_state(self, route, **updates):
        self.states[route].update(updates)
        (self.root / route / "completion.json").write_text(json.dumps(self.states[route]), encoding="utf-8")

    def partial(self, route):
        self.change_state(route, state="completed_partial_review",
                          review_counts={"requested": 2, "returned": 1, "missing": 1})

    def union(self, out, **kwargs):
        path = Path(out) / "supported_member_union.tsv"
        path.parent.mkdir(parents=True)
        path.write_text("candidate_gene_id\nG1\n", encoding="utf-8")
        self.union_path = path
        return {"candidate_set": str(path)}

    def wet(self, out, **kwargs):
        root = Path(out)
        root.mkdir()
        candidates = root / "wet_candidates.tsv"
        candidates.write_text("candidate_gene_id\tsupport_status\nG1\tsupported\nG2\treview_only\n", encoding="utf-8")
        self.wet_bytes = candidates.read_bytes()
        manifest = root / "wet_candidate_handoff_manifest.json"
        session._write(manifest, {"handoff_state": self.wet_state,
                                 "supported_member_union_sha256": session.sha256(self.union_path)})
        return {"manifest": str(manifest), "wet_candidates": str(candidates), "handoff_state": self.wet_state,
                "second_stage_state": "skipped_at_or_below_25"}

    def run_handoff(self, **kwargs):
        return session.handoff_sessions(str(self.root / "a"), str(self.root / "b"), str(self.out),
            genes_arg=str(self.resources), proteins_arg=str(self.resources), atlas_arg=str(self.resources),
            allow_unrecorded_route_b=True, **kwargs)

    def check_partial(self, routes):
        for route in routes:
            self.partial(route)
        before = {route: session._files(self.root / route) for route in ("a", "b")}
        result = self.run_handoff()
        self.assertEqual("first_stage_frozen", result["state"])
        self.assertEqual("partial", result["claim_review"]["state"])
        self.assertEqual("wet_reviewable_degraded_coverage", result["wet"]["handoff_state"])
        wm = session._read(Path(result["wet"]["manifest"]))
        self.assertEqual(result["claim_review"], wm["claim_review"])
        for route in ("a", "b"):
            self.assertEqual(int(route in routes), wm["claim_review"]["routes"][route]["missing"])
            self.assertEqual(before[route], session._files(self.root / route))
        self.assertEqual(self.wet_bytes, Path(result["wet"]["wet_candidates"]).read_bytes())
        kwargs = self.consumers["wet_candidate_handoff"].call_args.kwargs
        self.assertEqual(["a-review", self.states["b"]["candidates"]["review_candidates"]], kwargs["review_candidates_arg"])
        session._verify(self.out, result["files"])
        return result

    def test_partial_a_reaches_wet_without_changing_candidates(self):
        self.check_partial({"a"})

    def test_partial_b_reaches_wet_without_changing_candidates(self):
        self.check_partial({"b"})

    def test_both_partial_reach_wet_with_separate_counts(self):
        self.check_partial({"a", "b"})

    def test_complete_reviews_keep_normal_wet_state(self):
        result = self.run_handoff()
        self.assertEqual("complete", result["claim_review"]["state"])
        self.assertEqual("wet_reviewable", result["wet"]["handoff_state"])

    def _family_index(self, *, scope="reference_genome"):
        directory = self.root / "index"
        directory.mkdir()
        index = directory / "family_index.tsv"
        index.write_text(
            "gene_id\tassignment_status\tbest_subject\tfamily_stems\tindex_scope\n"
            f"G1\tassigned\tREF1\tFAM1\t{scope}\n", encoding="utf-8")
        (directory / "family_index_manifest.json").write_text(
            json.dumps({"index_scope": scope, "reference_label": "synthetic_reference"}) + "\n",
            encoding="utf-8")
        return index

    def test_handoff_wires_explicit_family_index_and_hashes_manifest(self):
        index = self._family_index()
        result = self.run_handoff(family_index_arg=str(index))
        call = self.consumers["route_a_interval_family_match"].call_args.kwargs
        self.assertEqual(str(index.resolve()), call["family_index_arg"])
        resolution = result["execution_context"]["family_index_resolution"]
        self.assertEqual("explicit_cli", resolution["source"])
        self.assertEqual("valid", resolution["validation"]["status"])
        fingerprints = session._read(self.out / "inputs.json")["inputs"]
        for key, path in (("family_index", index),
                          ("family_index_manifest", index.parent / "family_index_manifest.json")):
            self.assertEqual(str(path.resolve()), fingerprints[key]["path"])
            self.assertEqual(session.sha256(path), fingerprints[key]["sha256"])
        self.assertEqual(result["execution_context"], session._read(
            Path(result["wet"]["manifest"]))["execution_context"])

    def test_auto_resolved_family_index_wires_from_environment(self):
        index = self._family_index()
        with patch.dict(os.environ, {"CANDIDATE_QTL_FAMILY_INDEX": str(index)}), \
                patch("candidate_qtl.family_index_resolution._default_candidates", return_value=[]):
            result = self.run_handoff()
        self.assertEqual(str(index.resolve()),
                         self.consumers["route_a_interval_family_match"].call_args.kwargs["family_index_arg"])
        self.assertEqual("environment", result["execution_context"]["family_index_resolution"]["source"])

    def test_invalid_auto_scope_remains_unevaluated_and_not_negative(self):
        index = self._family_index(scope="interval_candidates")
        with patch.dict(os.environ, {"CANDIDATE_QTL_FAMILY_INDEX": str(index)}), \
                patch("candidate_qtl.family_index_resolution._default_candidates", return_value=[]):
            result = self.run_handoff()
        self.assertIsNone(
            self.consumers["route_a_interval_family_match"].call_args.kwargs["family_index_arg"])
        resolution = result["execution_context"]["family_index_resolution"]
        self.assertEqual("not_found", resolution["status"])
        attempt = next(row for row in resolution["attempts"] if row["source"] == "environment")
        self.assertEqual("invalid_reusable_index", attempt["status"])
        self.assertEqual("not_reusable_reference_genome_index", attempt["validation"]["reason"])
        self.assertNotIn("family_index", session._read(self.out / "inputs.json")["inputs"])

    def test_explicit_bad_scope_fails_before_any_handoff_output(self):
        index = self._family_index(scope="interval_candidates")
        with self.assertRaisesRegex(ValueError, "reference_genome family index"):
            self.run_handoff(family_index_arg=str(index))
        self.assertFalse(self.out.exists())

    def test_changed_family_manifest_is_detected_before_wet_completion(self):
        index = self._family_index()

        def mutate_index_manifest(*_args, **_kwargs):
            (index.parent / "family_index_manifest.json").write_text(
                json.dumps({"index_scope": "reference_genome", "changed": True}) + "\n",
                encoding="utf-8")
            return {"status": "interval-status"}

        self.consumers["route_a_interval_family_match"].side_effect = mutate_index_manifest
        with self.assertRaisesRegex(ValueError, "projection input changed"):
            self.run_handoff(family_index_arg=str(index))
        self.assertEqual("failed", session._read(self.out / "failure.json")["state"])
        self.assertFalse((self.out / "completion.json").exists())

    def test_partial_does_not_replace_a_more_restrictive_wet_state(self):
        self.partial("a")
        self.wet_state = "review_required_no_supported_candidates"
        result = self.run_handoff()
        self.assertEqual(self.wet_state, result["wet"]["handoff_state"])
        self.assertEqual("partial", result["claim_review"]["state"])

    def test_failed_or_awaiting_session_is_rejected_before_consumers(self):
        for state in ("failed", "awaiting_independent_review", "unknown"):
            with self.subTest(state=state):
                self.change_state("a", state=state)
                with self.assertRaises(ValueError):
                    self.run_handoff()
        self.consumers["freeze_blast_panel"].assert_not_called()
        self.assertFalse(self.out.exists())

    def test_inconsistent_counts_are_rejected(self):
        for counts in ({"requested": 2, "returned": 1, "missing": 0},
                       {"requested": 2, "returned": 1, "missing": 1},
                       {"requested": 1, "returned": True, "missing": 0},
                       {"requested": 1, "returned": 2, "missing": -1}, {}):
            with self.subTest(counts=counts):
                self.change_state("a", state="completed", review_counts=counts)
                with self.assertRaises(ValueError):
                    self.run_handoff()
        self.assertFalse(self.out.exists())

    def test_partial_state_with_zero_missing_is_rejected(self):
        self.change_state("a", state="completed_partial_review")
        with self.assertRaises(ValueError):
            self.run_handoff()
        self.assertFalse(self.out.exists())

    def test_partial_request_hash_tampering_is_still_rejected(self):
        self.partial("a")
        (self.root / "a" / "request.json").write_text(json.dumps({
            "route": "a", "case_id": "changed", "inputs": {"hypotheses": "input.tsv"}}), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "request anchor mismatch"):
            self.run_handoff()
        self.assertFalse(self.out.exists())

    def test_partial_frozen_artifact_tampering_is_still_rejected(self):
        self.partial("b")
        (self.root / "b" / "input.tsv").write_text("changed", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "artifact mismatch"):
            self.run_handoff()
        self.assertFalse(self.out.exists())

    def test_partial_downstream_failure_remains_failure(self):
        self.partial("a")
        self.consumers["run_panel_blastp"].side_effect = RuntimeError("synthetic provider failure")
        with self.assertRaisesRegex(RuntimeError, "provider failure"):
            self.run_handoff()
        self.assertEqual("failed", session._read(self.out / "failure.json")["state"])
        self.assertFalse((self.out / "completion.json").exists())


if __name__ == "__main__":
    unittest.main()
