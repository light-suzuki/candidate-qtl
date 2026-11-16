"""CLI wiring tests for build-executor-package (issue #365)."""

from __future__ import annotations

import contextlib
import io
import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from candidate_qtl.cli import build_parser, main, run
from candidate_qtl.core import runtime_path
from candidate_qtl.commands import COMMAND_NAMES

REPO = runtime_path(r"/path/to/home/AppData\Local\Temp\opencode\wt-cqtl-365-20260812162842")


def _git_init(root: Path) -> None:
    subprocess.run(["git", "init", "-q", "-b", "main", str(root)], check=True, capture_output=True)
    subprocess.run(["git", "-C", str(root), "config", "user.email", "test@example.com"], check=True, capture_output=True)
    subprocess.run(["git", "-C", str(root), "config", "user.name", "test"], check=True, capture_output=True)
    if not any(child for child in root.iterdir() if child.name != ".git"):
        (root / ".keep").write_text("", encoding="utf-8")
    subprocess.run(["git", "-C", str(root), "add", "-A"], check=True, capture_output=True)
    subprocess.run(["git", "-C", str(root), "commit", "-qm", "init"], check=True, capture_output=True)


class PackageCliTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)
        self.out = Path(tempfile.mkdtemp(prefix="cqtl-pkgcli-", dir=self.root.parent))
        self.addCleanup(lambda: shutil.rmtree(self.out, ignore_errors=True))
        _git_init(self.root)
        (self.root / "src/candidate_qtl").mkdir(parents=True)
        (self.root / "src/candidate_qtl/__init__.py").write_text("__version__ = '0.0.0'\n", encoding="utf-8")
        (self.root / "pyproject.toml").write_text('[project]\nname = "candidate-qtl"\nversion = "0.0.0"\n', encoding="utf-8")
        (self.root / "constraints.txt").write_text("pandas\n", encoding="utf-8")
        subprocess.run(["git", "-C", str(self.root), "add", "-A"], check=True, capture_output=True)
        subprocess.run(["git", "-C", str(self.root), "commit", "-qm", "fixture"], check=True, capture_output=True)

    def test_command_registered(self):
        self.assertIn("build-executor-package", COMMAND_NAMES)
        parser = build_parser()
        self.assertIn("build-executor-package", parser._subparsers._group_actions[0].choices)

    def test_build_success_exit_0(self):
        code = main(["build-executor-package", "--out-dir", str(self.out), "--repo-root", str(self.root)])
        self.assertEqual(code, 0)
        self.assertTrue((self.out / "package_manifest.json").is_file())
        self.assertTrue((self.out / "run_capsule.json").is_file())

    def test_build_json_output(self):
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            code = main(
                ["build-executor-package", "--out-dir", str(self.out), "--repo-root", str(self.root), "--json"]
            )
        self.assertEqual(code, 0)
        payload = json.loads(buf.getvalue())
        self.assertEqual(payload["command"], "build-executor-package")
        self.assertEqual(payload["exit_code"], 0)
        self.assertFalse(payload["truth_inputs_present"])
        self.assertEqual(payload["rejected_count"], 0)
        self.assertIn("file_count", payload)

    def test_scanner_rejection_exit_1_with_rule_ids(self):
        (self.root / "src/data").mkdir(parents=True)
        (self.root / "src/data/leak.tsv").write_text("case_id\tanswer_gene_id\nc1\tg1\n", encoding="utf-8")
        with contextlib.redirect_stderr(io.StringIO()) as err:
            code = run(["build-executor-package", "--out-dir", str(self.out), "--repo-root", str(self.root), "--include", "src/**"])
        self.assertEqual(code, 1)
        self.assertIn("rejected 1 file(s)", err.getvalue())
        self.assertIn("C1", err.getvalue())

    def test_missing_out_dir_required(self):
        with self.assertRaises(SystemExit):
            main(["build-executor-package", "--repo-root", str(self.root)])

    def test_out_dir_inside_repo_rejected(self):
        with self.assertRaises(ValueError):
            main(["build-executor-package", "--out-dir", str(self.root / "out"), "--repo-root", str(self.root)])

    def test_unmatched_include_rejected(self):
        with self.assertRaises(ValueError):
            main(
                [
                    "build-executor-package",
                    "--out-dir", str(self.out),
                    "--repo-root", str(self.root),
                    "--include", "no/such/**",
                ]
            )

    def test_run_command_recorded_in_capsule(self):
        code = main(
            [
                "build-executor-package", "--out-dir", str(self.out), "--repo-root", str(self.root),
                "--run-command", "candidate-qtl audit-predict --manifest m.tsv",
            ]
        )
        self.assertEqual(code, 0)
        capsule = json.loads((self.out / "run_capsule.json").read_text(encoding="utf-8"))
        self.assertEqual(
            capsule["build"]["declared_executor_command"],
            ["candidate-qtl", "audit-predict", "--manifest", "m.tsv"],
        )

    def test_help_lists_command(self):
        with contextlib.redirect_stdout(io.StringIO()) as buf, self.assertRaises(SystemExit):
            main(["--help"])
        self.assertIn("build-executor-package", buf.getvalue())

    def test_windows_style_paths_converted(self):
        if not REPO.is_dir():
            self.skipTest("candidate-qtl worktree not present in this environment")
        win_repo = r"/path/to/home/AppData\Local\Temp\opencode\wt-cqtl-365-20260812162842"
        win_out = win_repo + "-pkg-cli"
        cleanup = runtime_path(win_out)
        self.addCleanup(lambda: shutil.rmtree(cleanup, ignore_errors=True))
        code = main(["build-executor-package", "--out-dir", win_out, "--repo-root", win_repo])
        self.assertEqual(code, 0)
        manifest = json.loads((cleanup / "package_manifest.json").read_text(encoding="utf-8"))
        self.assertIn("src/candidate_qtl/package_build.py", manifest["files"])
        self.assertNotIn("\\", str(cleanup / "package_manifest.json").replace(os.sep, "/"))


if __name__ == "__main__":
    unittest.main()
