"""Production entry/forwarding regression, no benchmark truth or model calls.

Biology producers are mocked in the integration fixture. These tests establish
routing and provenance, not correctness of the existing copy-discrimination
policy or biological performance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from candidate_qtl import copy_evidence_inputs as evidence
from candidate_qtl import canonical_route_b_wet as canonical
from candidate_qtl import canonical_family_semantic_wet as semantic
from candidate_qtl.commands import orchestration as cli


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return str(path)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class CopyInputsTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.inputs = {axis: write(self.root / (axis + ".tsv"), "gene_id\tstatus\nx\tnot_evaluated\n")
                       for axis in evidence.COPY_EVIDENCE_AXES}

    def test_no_input_no_pin_or_hidden_discovery(self):
        self.assertEqual(({}, {}), evidence.bind_copy_evidence(self.root, None))
        self.assertFalse((self.root / "b04_canonical_semantic").exists())

    def test_all_six_native_names_and_full_hashes(self):
        kwargs, provenance = evidence.bind_copy_evidence(self.root, self.inputs)
        self.assertEqual({axis + "_arg" for axis in self.inputs}, set(kwargs))
        record = provenance["copy_evidence"]
        self.assertFalse(record["fingerprint_is_biological_verification"])
        self.assertFalse(record["missing_is_negative"])
        self.assertEqual([], record["omitted_axes"])
        for axis, path in self.inputs.items():
            self.assertEqual(sha(path), record["inputs"][axis]["sha256"])

    def test_omitted_axes_are_not_empty_negative_files(self):
        kwargs, p = evidence.bind_copy_evidence(self.root, {"expression": self.inputs["expression"]})
        self.assertEqual({"expression_arg"}, set(kwargs))
        self.assertEqual(5, len(p["copy_evidence"]["omitted_axes"]))

    def test_repeated_same_input_keeps_pin_bytes(self):
        k, p = evidence.bind_copy_evidence(self.root, self.inputs)
        path = Path(p["copy_evidence"]["input_manifest"])
        before = path.read_bytes(), path.stat().st_mtime_ns
        self.assertEqual((k, p), evidence.bind_copy_evidence(self.root, self.inputs))
        self.assertEqual(before, (path.read_bytes(), path.stat().st_mtime_ns))

    def test_resume_changed_content_refused(self):
        evidence.bind_copy_evidence(self.root, self.inputs)
        Path(self.inputs["synteny"]).write_text("changed", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "changed"):
            evidence.bind_copy_evidence(self.root, self.inputs)

    def test_resume_omitted_input_refused(self):
        evidence.bind_copy_evidence(self.root, self.inputs)
        with self.assertRaisesRegex(ValueError, "omitted"):
            evidence.bind_copy_evidence(self.root, None)

    def test_existing_unpinned_wet_not_overwritten_by_new_axes(self):
        wet = Path(write(self.root / "z04_canonical_wet_handoff/wet_candidates.tsv", "gene_id\nx\n"))
        with self.assertRaisesRegex(ValueError, "new output revision"):
            evidence.bind_copy_evidence(self.root, self.inputs)
        self.assertEqual("gene_id\nx\n", wet.read_text())

    def test_unknown_axis_cannot_override_family_universe(self):
        with self.assertRaises(ValueError):
            evidence.inspect_copy_evidence({"families": self.inputs["orthology"]})

    def test_wrong_container_type_refused(self):
        with self.assertRaises(ValueError):
            evidence.inspect_copy_evidence([])

    def test_missing_explicit_path_fails_not_silently_omitted(self):
        with self.assertRaisesRegex(ValueError, "does not exist"):
            evidence.inspect_copy_evidence({"expression": str(self.root / "absent")})

    def test_empty_explicit_path_refused(self):
        with self.assertRaises(ValueError):
            evidence.inspect_copy_evidence({"expression": " "})

    def test_directory_is_not_evidence_file(self):
        with self.assertRaises(ValueError):
            evidence.inspect_copy_evidence({"expression": str(self.root)})

    def test_native_missing_provider_status_is_not_reinterpreted(self):
        path = Path(self.inputs["expression"])
        before = path.read_bytes()
        k, p = evidence.bind_copy_evidence(self.root, {"expression": str(path)})
        evidence.verify_copy_evidence(k, p)
        self.assertEqual(before, path.read_bytes())
        self.assertEqual("not_evaluated_before_family_to_copy", p["copy_evidence"]["consumer_state"])

    def test_changed_during_consumer_execution_refused(self):
        k, p = evidence.bind_copy_evidence(self.root, self.inputs)
        Path(self.inputs["target"]).write_text("new target", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "during execution"):
            evidence.verify_copy_evidence(k, p)


class CanonicalCopyRoutingTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name) / "run"
        self.root.mkdir()
        self.paths = {a: write(self.root.parent / (a + ".tsv"), "gene_id\tstatus\ngA\tnot_evaluated\ngB\tone_to_many\n")
                      for a in evidence.COPY_EVIDENCE_AXES}
        # Real production input size deliberately differs from benchmark N=1000.
        self.atlas = write(self.root / "b01_route_b/interval_gene_atlas.tsv",
                           "gene_id\n" + "gA\ngB\n" + "".join(f"dark{i}\n" for i in range(1245)))
        self.sroot = self.root / "b04_canonical_semantic"
        self.bundle = write(self.sroot / "01_bundle/route_b_semantic_bundle.json", "{}")
        self.prompts = write(self.sroot / "01_bundle/route_b_semantic_prompts.jsonl", "{}\n")
        self.response = write(self.sroot / "02_llm_run/llm_candidate_responses.json", "[]")
        write(self.sroot / "02_llm_run/llm_run_manifest.json", "{}")
        self.union = write(self.sroot / "07_supported_union_for_family/supported.tsv", "gene_id\ngA\ngB\n")
        self.family = write(self.sroot / "08_canonical_family_recovery/family.tsv", "family\tmembers\nF\tgA;gB\n")
        self.bridge = write(self.sroot / "08_canonical_family_recovery/bridge.tsv", "gene_id\ngA\ngB\n")
        self.family_count = 1
        self.copy_kwargs = None
        self.handoff_kwargs = None
        self.fail_copy = False
        self.patch(canonical, "validate_route_b_semantic_artifacts", return_value={"test_fixture": "not_real_validation"})
        self.patch(canonical, "build_route_b_semantic_bundle", return_value={"bundle": self.bundle, "prompts": self.prompts})
        self.patch(canonical, "import_route_b_semantic_proposals", return_value={"proposals": "proposal", "manifest": "import"})
        self.patch(canonical, "verify_route_b_semantic_sources", return_value={"source_verification": "source", "manifest": "source-manifest"})
        self.patch(canonical, "verify_route_b_semantic_transfer", return_value={"claims": "claim", "manifest": "transfer"})
        self.patch(canonical, "freeze_route_b_semantic_candidates", return_value={"candidates": self.union, "manifest": "b-manifest", "candidate_gene_count": 2})
        self.patch(canonical, "build_supported_member_union", return_value={"candidate_set": self.union})
        self.patch(canonical, "build_canonical_family_recovery", side_effect=lambda *a, **k: {
            "family_count": self.family_count, "family_recovery": self.family, "family_bridge": self.bridge,
            "comparison": self.family, "manifest": "recovery-manifest"})
        self.copy_mock = self.patch(canonical, "family_to_copy", side_effect=self.copy_stage)
        self.patch(canonical, "wet_candidate_handoff", side_effect=self.handoff)
        self.patch(semantic, "build_family_semantic_request", return_value={"prompt_count": 0})
        self.patch(cli, "resolve_family_index", return_value=(None, {"attempted": False}))
        self.profile_mock = self.patch(cli, "guided_profile", side_effect=self.profile)

    def patch(self, obj, name, **kw):
        p = mock.patch.object(obj, name, **kw)
        result = p.start()
        self.addCleanup(p.stop)
        return result

    def copy_stage(self, **kwargs):
        self.copy_kwargs = kwargs
        if self.fail_copy:
            raise ValueError("native schema error")
        out = Path(kwargs["out_dir_arg"])
        write(out / "family_copy_discrimination.tsv", "gene_id\tstate\ngA\tco_preferred\ngB\tunresolved\n")
        write(out / "family_copy_shortlist.tsv", "gene_id\ngA\ngB\n")
        write(out / "family_copy_manifest.json", json.dumps({"inputs": {k: v for k, v in kwargs.items() if k.endswith("_arg")}}))
        return {}

    def handoff(self, out, **kwargs):
        self.handoff_kwargs = kwargs
        wet = write(Path(out) / "wet_candidates.tsv", "gene_id\ngA\ngB\n")
        manifest = write(Path(out) / "wet_manifest.json", json.dumps({"supported_member_union_sha256": sha(self.union)}))
        return {"wet_candidates": wet, "manifest": manifest, "handoff_state": "wet_reviewable"}

    def profile(self, *args, **kwargs):
        path = write(self.root / "guided_manifest.json", "{}")
        return {"manifest": path}

    def args(self, extra=()):
        parser = argparse.ArgumentParser()
        cli.register(parser.add_subparsers(dest="command", required=True))
        return parser.parse_args(["guided-profile", "--profile", "triple-route", "--out-dir", str(self.root),
                                  "--phenotype", "observed organ morphology", "--execute", *extra])

    def run_core(self, values=None):
        return canonical.run_canonical_route_b_wet_extension(str(self.root), phenotype="observed organ morphology",
                                                              copy_evidence_args=values)

    def test_real_cli_to_wrapper_to_core_to_native_consumer_all_six(self):
        flags = [v for axis, path in self.paths.items() for v in ("--copy-" + axis.replace("_", "-"), path)]
        before = Path(self.atlas).read_bytes()
        result = cli.HANDLERS["guided-profile"](self.args(flags))
        for axis, path in self.paths.items():
            self.assertEqual(path, self.copy_kwargs[axis + "_arg"])
        self.assertEqual(self.family, self.copy_kwargs["families_arg"])
        self.assertEqual(self.bridge, self.copy_kwargs["bridge_arg"])
        self.assertEqual(before, Path(self.atlas).read_bytes())
        self.assertEqual(1247, len(before.splitlines()) - 1)
        # The copy-only inputs never leak into phenotype/A or B0 generation.
        self.assertFalse(any(k.startswith("copy") for k in self.profile_mock.call_args.kwargs))
        manifests = [result["manifest"], result["canonical_route_b_wet"]["manifest"],
                     result["canonical_route_b_wet"]["base_manifest"]]
        for path in manifests:
            data = json.loads(Path(path).read_text())
            self.assertEqual(6, len(data["copy_evidence"]["inputs"]))
            self.assertEqual("passed_to_family_to_copy", data["copy_evidence"]["consumer_state"])
        self.assertEqual("gene_id\ngA\ngB\n", Path(result["canonical_route_b_wet"]["wet_candidates"]).read_text())

    def test_no_optional_arguments_preserve_native_call_shape(self):
        result = self.run_core()
        self.assertEqual({"families_arg", "out_dir_arg", "atlas_arg", "case_id", "bridge_arg"}, set(self.copy_kwargs))
        self.assertNotIn("copy_evidence", json.loads(Path(result["manifest"]).read_text()))

    def test_only_supplied_axis_passed_without_other_required_inputs(self):
        self.run_core({"synteny": self.paths["synteny"]})
        self.assertIn("synteny_arg", self.copy_kwargs)
        self.assertNotIn("target_arg", self.copy_kwargs)
        self.assertNotIn("expression_arg", self.copy_kwargs)

    def test_plain_b0_orthology_not_reused_as_copy_evidence(self):
        cli.HANDLERS["guided-profile"](self.args(["--orthology", self.paths["orthology"]]))
        self.assertNotIn("orthology_arg", self.copy_kwargs)

    def test_c_fallback_also_reaches_native_copy_consumer(self):
        Path(self.response).unlink()
        (self.sroot / "02_llm_run/llm_run_manifest.json").unlink()
        result = canonical.run_canonical_route_b_wet_extension(
            str(self.root), phenotype="observed", route_c_members_arg="c-members.tsv", copy_evidence_args=self.paths)
        self.assertIn("expression_arg", self.copy_kwargs)
        self.assertIsNone(result["route_b_verified_candidates"])
        self.assertIn("copy_evidence", json.loads(Path(result["manifest"]).read_text()))

    def test_zero_family_does_not_promote_evidence_rows_into_candidates(self):
        self.family_count = 0
        result = self.run_core(self.paths)
        self.copy_mock.assert_not_called()
        self.assertIsNone(result["family_copy"])
        data = json.loads(Path(result["manifest"]).read_text())
        self.assertEqual("not_evaluated_before_family_to_copy", data["copy_evidence"]["consumer_state"])

    def test_no_b_response_records_inputs_without_claiming_consumption(self):
        Path(self.response).unlink()
        (self.sroot / "02_llm_run/llm_run_manifest.json").unlink()
        result = self.run_core(self.paths)
        self.assertEqual("semantic_review_required", result["state"])
        self.copy_mock.assert_not_called()
        self.assertEqual("not_evaluated_before_family_to_copy", json.loads(Path(result["manifest"]).read_text())["copy_evidence"]["consumer_state"])

    def test_missing_atlas_not_negative_and_no_copy_call(self):
        Path(self.atlas).unlink()
        result = self.run_core(self.paths)
        self.assertEqual("not_runnable_upstream_b0_missing", result["state"])
        self.copy_mock.assert_not_called()
        self.assertIn("copy_evidence", json.loads(Path(result["manifest"]).read_text()))

    def test_native_schema_failure_is_not_swallowed(self):
        self.fail_copy = True
        with self.assertRaisesRegex(ValueError, "native schema"):
            self.run_core(self.paths)
        self.assertIsNone(self.handoff_kwargs)

    def test_bad_explicit_input_fails_before_upstream_profile(self):
        with self.assertRaisesRegex(ValueError, "does not exist"):
            cli.HANDLERS["guided-profile"](self.args(["--copy-expression", str(self.root / "missing")]))
        self.profile_mock.assert_not_called()

    def test_unsupported_profile_is_not_silently_accepted(self):
        args = self.args(["--copy-expression", self.paths["expression"]])
        args.profile = "screening"
        with self.assertRaisesRegex(ValueError, "require --profile"):
            cli.HANDLERS["guided-profile"](args)
        self.profile_mock.assert_not_called()

    def test_dry_run_reports_planned_not_consumed(self):
        args = self.args(["--copy-expression", str(self.root / "not-created-yet.tsv")])
        args.execute = False
        result = cli.HANDLERS["guided-profile"](args)
        self.assertEqual("planned_not_consumed", result["copy_evidence"]["state"])
        self.copy_mock.assert_not_called()

    def test_changed_copy_input_rejected_before_b_reprocessing(self):
        result = self.run_core(self.paths)
        wet_before = Path(result["wet_candidates"]).read_bytes()
        Path(self.paths["position"]).write_text("changed")
        self.copy_mock.reset_mock()
        with self.assertRaisesRegex(ValueError, "changed"):
            self.run_core(self.paths)
        self.copy_mock.assert_not_called()
        self.assertEqual(wet_before, Path(result["wet_candidates"]).read_bytes())

    def test_mutation_during_copy_stage_fails_before_handoff(self):
        original = self.copy_stage
        def mutate(**kwargs):
            result = original(**kwargs)
            Path(self.paths["expression"]).write_text("mutated")
            return result
        self.copy_mock.side_effect = mutate
        with self.assertRaisesRegex(ValueError, "during execution"):
            self.run_core(self.paths)
        self.assertIsNone(self.handoff_kwargs)

    def test_cli_changed_evidence_rejected_before_upstream_side_effects(self):
        self.run_core(self.paths)
        Path(self.paths["expression"]).write_text("changed")
        flags = [v for a, p in self.paths.items() for v in ("--copy-" + a.replace("_", "-"), p)]
        with self.assertRaisesRegex(ValueError, "changed"):
            cli.HANDLERS["guided-profile"](self.args(flags))
        self.profile_mock.assert_not_called()

    def test_cli_omitted_evidence_rejected_before_upstream_side_effects(self):
        self.run_core(self.paths)
        with self.assertRaisesRegex(ValueError, "omitted"):
            cli.HANDLERS["guided-profile"](self.args())
        self.profile_mock.assert_not_called()

    def test_corrupt_pin_refused_without_overwrite(self):
        evidence.bind_copy_evidence(self.root, self.paths)
        pin = self.sroot / "canonical_copy_inputs.json"
        pin.write_text("{")
        with self.assertRaises(ValueError):
            self.run_core(self.paths)
        self.assertEqual("{", pin.read_text())
        self.copy_mock.assert_not_called()


if __name__ == "__main__":
    unittest.main()
