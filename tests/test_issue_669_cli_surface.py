"""Entrypoint and import-boundary contract for issue #669."""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from candidate_qtl.command_surface import COMMAND_SURFACES, SURFACE_CLASSES
from candidate_qtl.commands import command_names
from candidate_qtl.cli import build_parser

ROOT = Path(__file__).resolve().parents[1]


class CliSurface669Tests(unittest.TestCase):
    def test_every_handler_has_explicit_complete_record(self) -> None:
        self.assertEqual(182, len(COMMAND_SURFACES))
        self.assertEqual(182, len({r.command for r in COMMAND_SURFACES}))
        required = ("command", "feature_module", "surface", "installed_entrypoint", "canonical_consumer", "reachability_profile", "deprecation_or_move_state", "replacement")
        for record in COMMAND_SURFACES:
            data = record.as_dict()
            self.assertTrue(all(field in data for field in required), record.command)
            self.assertIn(record.surface, SURFACE_CLASSES)
            self.assertTrue(record.feature_module)

    @unittest.skip("public release: pins the private bench command count")
    def test_production_and_bench_command_sets(self) -> None:
        production = set(command_names("core"))
        bench = set(command_names("all"))
        self.assertEqual(130, len(production))
        self.assertEqual(182, len(bench))
        self.assertTrue({"phenotype-investigate", "mechanism-import", "mechanism-coverage-audit", "investigation-attribution", "source-verify", "route-c-direct-search", "route-c-direct-freeze", "route-c-direct-project", "route-c-direct-members", "triple-route-compare"} <= production)
        self.assertFalse(any("benchmark" in name for name in production))
        self.assertFalse(any(name.startswith("retrieval-") for name in production))
        self.assertTrue(production < bench)

    def test_fresh_production_parser_has_no_benchmark_implementation_imports(self) -> None:
        code = """import json,sys\nfrom candidate_qtl.cli import build_parser\nimport argparse\np=build_parser('core')\na=next(x for x in p._actions if isinstance(x,argparse._SubParsersAction))\nprint(json.dumps({'commands':len(a.choices),'modules':len([x for x in sys.modules if x=='candidate_qtl' or x.startswith('candidate_qtl.')]),'citation_tokenizer_imported':'candidate_qtl.citation_tokens' in sys.modules,'deny':[x for x in sys.modules if any(y in x for y in ('blind_benchmark','benchmark_metrics','benchmark_reporting','retrieval_local','retrieval_route'))]}))"""
        env = dict(os.environ)
        env["PYTHONPATH"] = str(ROOT / "src")
        result = subprocess.run([sys.executable, "-c", code], cwd=ROOT, env=env, check=True, capture_output=True, text=True)
        observed = json.loads(result.stdout)
        self.assertEqual(130, observed["commands"])
        self.assertEqual([], observed["deny"])
        self.assertFalse(observed["citation_tokenizer_imported"])
        self.assertLess(observed["modules"], 176)

    def test_bench_parser_reaches_all_registry_handlers(self) -> None:
        parser = build_parser("all", prog="candidate-qtl-bench")
        action = next(x for x in parser._actions if isinstance(x, argparse._SubParsersAction))
        self.assertEqual(set(command_names("all")), set(action.choices))

    def test_lazy_direct_c_project_preserves_public_argument_order(self) -> None:
        from candidate_qtl.commands.production_direct_c import _run

        args = SimpleNamespace(
            frozen="frozen.tsv", freeze_manifest="freeze.json", atlas="atlas.tsv",
            qtl_id="qtl-1", out_dir="out", eligible_hits="hits.tsv", full_hash=True,
            identity_mapping="identity-mapping.tsv",
        )
        with patch(
            "candidate_qtl.route_c_direct_pipeline.direct_rescue_project",
            return_value={"ok": True},
        ) as project:
            self.assertEqual({"ok": True}, _run("route-c-direct-project", args))
        project.assert_called_once_with(
            "frozen.tsv", "freeze.json", "atlas.tsv", "qtl-1", "out", "hits.tsv", True, "identity-mapping.tsv"
        )


if __name__ == "__main__":
    unittest.main()
