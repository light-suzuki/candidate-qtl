"""Exercise the actual guided CLI parser/handler with unrelated upstream
producers stubbed. Fresh homology logic and file I/O run unmocked; the provider
is unavailable here. Native BLAST interface tests are in the companion file.
"""
from __future__ import annotations
import argparse
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest import mock

from candidate_qtl import second_stage_recheck as recheck
from candidate_qtl.commands import orchestration as installed_cli


class GuidedIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);self.output=self.root/'run';self.output.mkdir()
        self.request=self.output/'adaptive_synthesis_request.json'
        ids=[f'g{i:02}' for i in range(26)]
        self.request.write_text(json.dumps({'stage':'adaptive_synthesis_request','candidate_count':26,'candidate_gene_ids':ids,'case_id':'fixture'}))
        self.wet=self.output/'wet_candidates.tsv';self.wet.write_text('gene_id\n'+'\n'.join(ids)+'\n')
        self.handoff=self.output/'wet_candidate_handoff_manifest.json'
        self.handoff.write_text(json.dumps({'second_stage_state':'requested','adaptive_synthesis_request':str(self.request)}))
        self.outer=self.output/'guided_manifest.json';self.outer.write_text('{}')
        self.proteins=self.root/'proteins.fa';self.proteins.write_text('>g00\nMACDEFGHI\n')
        self.profile=mock.Mock(return_value={'manifest':str(self.outer)})
        self.canonical=mock.Mock(return_value={'wet_handoff_manifest':str(self.handoff),'wet_candidates':str(self.wet)})
        attrs={
            'canonical_family_semantic_wet':{'run_canonical_route_b_family_semantic_wet_extension':self.canonical},
            'family_index_resolution':{'resolve_family_index':lambda *a:(None,{'attempted':False})},
            'orchestration':{'orchestrate_plan':mock.Mock()},
            'guided_profiles':{'PROFILE_NAMES':('triple-route','screening'),'guided_profile':self.profile},
            'phenotype_investigation':{'PAYLOAD_SOURCES':('llm','agent','human')},
            'route_c_direct_pipeline':{'direct_rescue_project':mock.Mock()},
            'route_c_direct_members':{'build_route_c_direct_members':mock.Mock()},
            'workflow_extensions':{'case_run':mock.Mock()},
            'commands._shared':{'add_json_flag':lambda p:p.add_argument('--json',action='store_true')},
        }
        stubs={}
        for name,values in attrs.items():
            full='candidate_qtl.'+name;module=types.ModuleType(full);module.__dict__.update(values);stubs[full]=module
        cp=types.ModuleType('candidate_qtl.commands');cp.__path__=[];stubs['candidate_qtl.commands']=cp
        patch=mock.patch.dict(sys.modules,stubs);patch.start();self.addCleanup(patch.stop)
        source=Path(installed_cli.__file__)
        spec=importlib.util.spec_from_file_location('candidate_qtl.commands._test_orchestration',source)
        self.cli=importlib.util.module_from_spec(spec);spec.loader.exec_module(self.cli)
    def args(self,execute=True,profile='triple-route'):
        parser=argparse.ArgumentParser();self.cli.register(parser.add_subparsers(dest='command',required=True))
        return parser.parse_args(['guided-profile','--profile',profile,'--out-dir',str(self.output),
            '--phenotype','observed morphology','--proteins',str(self.proteins),'--protein-db',str(self.root/'db'),
            *(['--execute'] if execute else [])])
    def test_actual_cli_existing_arguments_reach_recheck_and_persist_link(self):
        before=self.wet.read_bytes(),self.request.read_bytes()
        with mock.patch.object(recheck.shutil,'which',return_value=None):
            result=self.cli.HANDLERS['guided-profile'](self.args())
        self.assertEqual(26,result['second_stage_recheck']['candidate_count'])
        self.assertEqual('not_evaluated_missing_blastp',result['second_stage_recheck']['reference_search']['state'])
        self.assertEqual(26,result['second_stage_recheck']['prompt_count'])
        self.assertEqual(before,(self.wet.read_bytes(),self.request.read_bytes()))
        outer=json.loads(self.outer.read_text());self.assertEqual(result['second_stage_recheck']['manifest'],outer['second_stage_recheck']['manifest'])
        self.assertEqual(0,outer['second_stage_recheck']['model_calls'])
    def test_existing_skip25_result_has_no_new_work(self):
        self.handoff.write_text(json.dumps({'second_stage_state':'skipped_at_or_below_25'}))
        with mock.patch.object(recheck,'fresh_recheck') as work:
            result=self.cli.HANDLERS['guided-profile'](self.args())
        work.assert_not_called();self.assertNotIn('second_stage_recheck',result)
    def test_dry_run_has_no_recheck_or_canonical_execution(self):
        with mock.patch.object(recheck,'fresh_recheck') as work:
            result=self.cli.HANDLERS['guided-profile'](self.args(execute=False))
        work.assert_not_called();self.canonical.assert_not_called()
        self.assertNotIn('second_stage_recheck',result)
    def test_screening_has_no_added_second_stage(self):
        with mock.patch.object(recheck,'fresh_recheck') as work:
            self.cli.HANDLERS['guided-profile'](self.args(profile='screening'))
        work.assert_not_called();self.canonical.assert_not_called()
    def test_parser_surface_is_unchanged(self):
        parser=argparse.ArgumentParser();self.cli.register(parser.add_subparsers(dest='command',required=True))
        args=self.args()
        self.assertNotIn('second_stage_recheck',vars(args))
        self.assertEqual('guided-profile',args.command)


class CurrentBoundaryGuidedTests(unittest.TestCase):
    def fixture(self, root, count):
        from tests.test_wet_relevance_selection import handoff

        wet = handoff(root / "producer", (0, 0, count, 0, 0))
        out = root / "guided"
        out.mkdir()
        outer = out / "guided_manifest.json"
        outer.write_text(json.dumps({"overall_status": "completed", "final": True}))
        return wet, out, outer

    def run_cli(self, wet, out, outer, canonical_override=None):
        from candidate_qtl.cli import run

        upstream = {"manifest": str(outer), "overall_status": "completed", "exit_code": 0}
        canonical = {"wet_handoff_manifest": wet["manifest"], "wet_candidates": wet["wet_candidates"],
                     "state": wet["handoff_state"]}
        canonical.update(canonical_override or {})
        canonical["manifest"] = str(out / "canonical.json")
        Path(canonical["manifest"]).write_text(json.dumps(canonical))
        output, errors = io.StringIO(), io.StringIO()
        with mock.patch.object(installed_cli, "guided_profile", return_value=upstream), \
             mock.patch.object(installed_cli, "run_canonical_route_b_family_semantic_wet_extension", return_value=canonical), \
             mock.patch.object(recheck, "fresh_recheck") as legacy, \
             contextlib.redirect_stdout(output), contextlib.redirect_stderr(errors):
            code = run(["guided-profile", "--profile", "triple-route", "--execute", "--case-id", "case1",
                        "--focal-species", "Fixture plant", "--phenotype", "altered leaf shape",
                        "--out-dir", str(out), "--json"])
        legacy.assert_not_called()
        return code, output.getvalue(), errors.getvalue()

    def test_real_25_and_26_wet_producers_use_current_boundary_consumer(self):
        for count in (25, 26):
            with self.subTest(count=count), tempfile.TemporaryDirectory() as tmp:
                wet, out, outer = self.fixture(Path(tmp), count)
                before = Path(wet["wet_candidates"]).read_bytes()
                code, output, errors = self.run_cli(wet, out, outer)
                self.assertEqual(0, code, errors)
                result = json.loads(output)
                self.assertEqual(before, Path(wet["wet_candidates"]).read_bytes())
                if count == 26:
                    link = result["boundary_review"]
                    self.assertEqual(26, link["candidate_count"])
                    self.assertEqual("awaiting_review", link["state"])
                    self.assertEqual("not_executed", link["review_state"])
                    self.assertEqual(0, link["model_calls"])
                    self.assertEqual("candidate_qtl.boundary_tier_review.consume_boundary_selection", link["next_consumer"])
                    self.assertEqual(link, json.loads(outer.read_text())["boundary_review"])
                    self.assertEqual("partial", result["overall_status"])
                    self.assertFalse(result["final"])
                    self.assertTrue(result["canonical_execution"]["required_work_remaining"])
                else:
                    self.assertNotIn("boundary_review", result)
                    self.assertEqual("completed", result["overall_status"])
                    self.assertTrue(result["final"])
                    self.assertEqual("first_stage", result["canonical_execution"]["completed_stage"])

    def test_semantic_pending_is_persisted_with_identity_and_resume_location(self):
        from candidate_qtl.core import sha256

        with tempfile.TemporaryDirectory() as tmp:
            wet, out, outer = self.fixture(Path(tmp), 1)
            material = out / "material.json"
            material.write_text('{"fixture":true}')
            code, output, errors = self.run_cli(wet, out, outer, {
                "state": "semantic_review_required", "wet_candidates": None, "wet_handoff_manifest": None,
                "expected_semantic_run_dir": str(out / "02_llm_run"),
                "prompts": str(material), "bundle": str(material),
            })
            self.assertEqual(0, code, errors)
            result, saved = json.loads(output), json.loads(outer.read_text())
            self.assertEqual("partial", saved["overall_status"])
            self.assertFalse(saved["final"])
            self.assertEqual("completed", saved["legacy_profile"]["overall_status"])
            self.assertTrue(saved["legacy_profile"]["final"])
            canonical = saved["canonical_execution"]
            self.assertEqual(canonical, result["canonical_execution"])
            self.assertEqual("semantic_review_required", canonical["state"])
            self.assertTrue(canonical["required_work_remaining"])
            self.assertIsNone(canonical["wet_candidates"])
            self.assertEqual(sha256(Path(canonical["manifest"]["path"])), canonical["manifest"]["sha256"])
            self.assertEqual(str(out / "02_llm_run"), canonical["expected_semantic_run_dir"])
            self.assertEqual(sha256(material), canonical["bundle"]["sha256"])

    def test_invalid_stage_hash_and_ids_still_reject_without_changing_wet(self):
        for damage in ("stage", "scope_sha", "ids"):
            with self.subTest(damage=damage), tempfile.TemporaryDirectory() as tmp:
                wet, out, outer = self.fixture(Path(tmp), 26)
                path = Path(wet["adaptive_synthesis_request"])
                request = json.loads(path.read_text())
                if damage == "stage":
                    request["stage"] = "unrecognized_request"
                elif damage == "scope_sha":
                    request["scope_manifest"]["sha256"] = "0" * 64
                else:
                    request["candidate_gene_ids"][0] = "not_in_scope"
                path.write_text(json.dumps(request))
                before = Path(wet["wet_candidates"]).read_bytes()
                code, _, _ = self.run_cli(wet, out, outer)
                self.assertNotEqual(0, code)
                self.assertEqual(before, Path(wet["wet_candidates"]).read_bytes())


if __name__=='__main__':unittest.main()
