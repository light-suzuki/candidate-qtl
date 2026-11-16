"""The process entry point's error contract (#190).

`main` is the library entry point and must keep re-raising so callers and the
rest of the test-suite see the real exception. `run` is what the console script
and `python -m candidate_qtl` execute, and it must turn an expected input error
into one actionable line plus a non-zero exit -- never a traceback that buries
the message explaining how to fix the input.

Subprocess cases run the real installed entry points, because the defect this
guards against was a mismatch between `pyproject.toml`, `__main__.py` and the
function that actually handles errors.
"""

from __future__ import annotations

import io
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from candidate_qtl import cli
from candidate_qtl.cli import EXIT_EXPECTED_FAILURE, EXIT_INTERNAL_ERROR, EXIT_OK, main, run


REPO = Path(__file__).resolve().parents[1]
TUTORIAL = REPO / "docs" / "tutorial"


def _run_module(*arguments: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "candidate_qtl", *arguments],
        capture_output=True,
        text=True,
        cwd=REPO,
    )


class MainKeepsRaisingTests(unittest.TestCase):
    """`main` is the library contract: exceptions propagate unchanged."""

    def test_main_reraises_expected_errors(self) -> None:
        buffer = io.StringIO()
        with redirect_stderr(buffer), self.assertRaises(ValueError):
            main(["extract", "--fasta", "x.fa", "--gff3", "x.gff3", "--interval", "chr1:abc", "--out-dir", "out"])
        self.assertIn("error: interval must be CHROM:START-END", buffer.getvalue())


class RunExitCodeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def _run(self, arguments: list[str]) -> tuple[int, str, str]:
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            code = run(arguments)
        return code, out.getvalue(), err.getvalue()

    def test_malformed_interval_exits_one_without_a_traceback(self) -> None:
        code, _, err = self._run(
            ["extract", "--fasta", str(TUTORIAL / "genome.fa"), "--gff3", str(TUTORIAL / "genes.gff3"),
             "--interval", "chr1:abc", "--out-dir", str(self.root / "out")]
        )
        self.assertEqual(EXIT_EXPECTED_FAILURE, code)
        self.assertIn("error: interval must be CHROM:START-END with 1-based inclusive coordinates", err)
        self.assertNotIn("Traceback", err)

    def test_missing_fasta_is_an_expected_error_naming_the_argument(self) -> None:
        code, _, err = self._run(
            ["extract", "--fasta", str(self.root / "absent.fa"), "--gff3", str(TUTORIAL / "genes.gff3"),
             "--interval", "chr1:1-900", "--out-dir", str(self.root / "out")]
        )
        self.assertEqual(EXIT_EXPECTED_FAILURE, code)
        # Names --fasta, not the .fai index pysam failed to build.
        self.assertIn("--fasta does not exist", err)
        self.assertNotIn("Traceback", err)
        self.assertNotIn("SamtoolsError", err)

    def test_unknown_seqid_keeps_its_actionable_hint_visible(self) -> None:
        code, _, err = self._run(
            ["extract", "--fasta", str(TUTORIAL / "genome.fa"), "--gff3", str(TUTORIAL / "genes.gff3"),
             "--interval", "chrZ:1-900", "--out-dir", str(self.root / "out")]
        )
        self.assertEqual(EXIT_EXPECTED_FAILURE, code)
        self.assertIn("--seqid-map GFF_SEQID=FASTA_SEQID", err)
        self.assertNotIn("Traceback", err)

    def test_success_exits_zero(self) -> None:
        code, _, _ = self._run(
            ["extract", "--fasta", str(TUTORIAL / "genome.fa"), "--gff3", str(TUTORIAL / "genes.gff3"),
             "--interval", "chr1:1-900", "--out-dir", str(self.root / "ok")]
        )
        self.assertEqual(EXIT_OK, code)

    def test_json_mode_keeps_a_machine_readable_error(self) -> None:
        code, _, err = self._run(
            ["extract", "--fasta", str(TUTORIAL / "genome.fa"), "--gff3", str(TUTORIAL / "genes.gff3"),
             "--interval", "chr1:abc", "--out-dir", str(self.root / "out"), "--json"]
        )
        self.assertEqual(EXIT_EXPECTED_FAILURE, code)
        import json

        payload = json.loads(err.strip().splitlines()[0])
        self.assertEqual("extract", payload["command"])
        self.assertIn("CHROM:START-END", payload["error"])
        self.assertNotIn("Traceback", err)

    def test_traceback_flag_restores_the_traceback(self) -> None:
        with redirect_stderr(io.StringIO()), self.assertRaises(ValueError):
            run(["--traceback", "extract", "--fasta", str(TUTORIAL / "genome.fa"),
                 "--gff3", str(TUTORIAL / "genes.gff3"), "--interval", "chr1:abc",
                 "--out-dir", str(self.root / "out")])

    def test_traceback_env_var_restores_the_traceback(self) -> None:
        import os
        from unittest import mock

        with mock.patch.dict(os.environ, {cli.TRACEBACK_ENV: "1"}):
            with redirect_stderr(io.StringIO()), self.assertRaises(ValueError):
                run(["extract", "--fasta", str(TUTORIAL / "genome.fa"), "--gff3", str(TUTORIAL / "genes.gff3"),
                     "--interval", "chr1:abc", "--out-dir", str(self.root / "out")])

    def test_env_var_off_values_do_not_enable_tracebacks(self) -> None:
        import os
        from unittest import mock

        for value in ("", "0", "false", "no", "  ", "FALSE"):
            with mock.patch.dict(os.environ, {cli.TRACEBACK_ENV: value}):
                self.assertFalse(cli.traceback_requested([]), f"{value!r} must not enable tracebacks")

    def test_unexpected_error_keeps_its_traceback_and_exits_two(self) -> None:
        from unittest import mock

        def _boom(_args: object) -> dict[str, object]:
            raise KeyError("internal invariant")

        with mock.patch.object(cli, "handler_for", return_value=_boom):
            code, _, err = self._run(["doctor"])
        self.assertEqual(EXIT_INTERNAL_ERROR, code)
        self.assertIn("Traceback", err)
        self.assertIn("internal error", err)

    def test_environment_failure_exits_one_naming_the_dependency(self) -> None:
        from unittest import mock

        from candidate_qtl.core import ExternalDependencyError

        def _network_failure(_args: object) -> dict[str, object]:
            raise ExternalDependencyError("network", "simulated service outage")

        with mock.patch.object(cli, "handler_for", return_value=_network_failure):
            code, _, err = self._run(["doctor"])
        self.assertEqual(EXIT_EXPECTED_FAILURE, code)
        self.assertIn("error: external dependency failure (kind=network): simulated service outage", err)
        self.assertNotIn("Traceback", err)
        self.assertNotIn("defect in candidate-qtl", err)

    def test_raw_network_and_timeout_errors_exit_one_with_kind(self) -> None:
        import subprocess
        import urllib.error
        from unittest import mock

        cases = [
            (urllib.error.URLError("connection refused"), "network"),
            (TimeoutError("server slow"), "timeout"),
            (subprocess.TimeoutExpired(["blastp"], 300), "timeout"),
        ]
        for exception, kind in cases:

            def _fail(_args: object, exc: Exception = exception) -> dict[str, object]:
                raise exc

            with mock.patch.object(cli, "handler_for", return_value=_fail):
                code, _, err = self._run(["doctor"])
            self.assertEqual(EXIT_EXPECTED_FAILURE, code, kind)
            self.assertIn(f"external dependency failure (kind={kind})", err)
            self.assertNotIn("defect in candidate-qtl", err)


class InstalledEntryPointTests(unittest.TestCase):
    """Both shipped entry points must agree; they disagreed before #190."""

    def test_module_entry_point_reports_expected_failure_cleanly(self) -> None:
        completed = _run_module(
            "extract", "--fasta", str(TUTORIAL / "genome.fa"), "--gff3", str(TUTORIAL / "genes.gff3"),
            "--interval", "chr1:abc", "--out-dir", "/tmp/candidate_qtl_cli_error_surface",
        )
        self.assertEqual(EXIT_EXPECTED_FAILURE, completed.returncode)
        self.assertIn("error: interval must be CHROM:START-END", completed.stderr)
        self.assertNotIn("Traceback", completed.stderr)

    def test_console_script_target_is_the_error_handling_entry_point(self) -> None:
        # Guards the defect this issue reported: pyproject pointed at `main`, so
        # the console script re-raised and `_cli` was unreachable dead code.
        import tomllib

        config = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))
        self.assertEqual("candidate_qtl.cli:run", config["project"]["scripts"]["candidate-qtl"])

    def test_help_and_version_still_exit_zero(self) -> None:
        for arguments in (["--help"], ["--version"]):
            completed = _run_module(*arguments)
            self.assertEqual(EXIT_OK, completed.returncode, arguments)


if __name__ == "__main__":
    unittest.main()
