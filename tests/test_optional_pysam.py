"""The package must import and half the CLI must run without pysam (#192).

pysam has no native Windows build, and it used to be imported at the top of
`core`, which everything imports. That made even `--help` and the table-only
commands unavailable on native Windows.

pysam is installed in this environment, so absence is simulated in a subprocess
with a meta-path finder that refuses the module. Running out-of-process is the
point: a same-process patch would leave the already-imported `pysam` in
`sys.modules` and prove nothing about import time.
"""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

from candidate_qtl.optional_deps import PYSAM_DEPENDENT_COMMANDS, pysam_available, require_pysam


REPO = Path(__file__).resolve().parents[1]
TUTORIAL = REPO / "docs" / "tutorial"

BLOCK_PYSAM = textwrap.dedent(
    """
    import sys

    class _RefusePysam:
        def find_spec(self, name, path=None, target=None):
            if name == "pysam" or name.startswith("pysam."):
                raise ImportError("pysam is unavailable in this simulated environment")
            return None

    sys.meta_path.insert(0, _RefusePysam())
    assert "pysam" not in sys.modules, "fixture is only valid before pysam is imported"
    """
)


def _without_pysam(body: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-c", BLOCK_PYSAM + textwrap.dedent(body)],
        capture_output=True,
        text=True,
        cwd=REPO,
    )


class RequirePysamTests(unittest.TestCase):
    def test_available_here(self) -> None:
        # Guards the fixture itself: the simulated-absence tests below only mean
        # something if pysam really is importable in this environment.
        self.assertTrue(pysam_available())
        self.assertEqual("pysam", require_pysam("extract").__name__)

    def test_dependent_command_list_names_real_commands(self) -> None:
        from candidate_qtl.commands import COMMAND_NAMES

        unknown = sorted(set(PYSAM_DEPENDENT_COMMANDS).difference(COMMAND_NAMES))
        self.assertEqual([], unknown, "PYSAM_DEPENDENT_COMMANDS must reference real subcommands")


class ImportWithoutPysamTests(unittest.TestCase):
    def test_package_and_parser_import(self) -> None:
        from candidate_qtl.commands import COMMAND_NAMES

        completed = _without_pysam(
            """
            import candidate_qtl
            from candidate_qtl.cli import build_parser
            parser = build_parser()
            action = next(a for a in parser._actions if getattr(a, "choices", None))
            print("commands", len(action.choices))
            """
        )
        self.assertEqual(0, completed.returncode, completed.stderr)
        self.assertIn(f"commands {len(COMMAND_NAMES)}", completed.stdout)

    def test_table_only_modules_import(self) -> None:
        modules = [
            "candidate_qtl.core",
            "candidate_qtl.hardening",
            "candidate_qtl.ranking_policy",
            "candidate_qtl.wet_state",
            "candidate_qtl.hypothesis_ledger",
            "candidate_qtl.next_experiment",
            "candidate_qtl.elimination",
            "candidate_qtl.external_rank_ledger",
            "candidate_qtl.schemas",
            # Imported for their CLI handlers even when the command needs pysam.
            "candidate_qtl.reference_extract",
            "candidate_qtl.workflow_extensions",
            "candidate_qtl.auto_triage",
            "candidate_qtl.annotation_rescue",
            "candidate_qtl.exon_assays",
            "candidate_qtl.benchmarks",
            "candidate_qtl.retrieval_route",
            "candidate_qtl.retrieval_semantic",
            "candidate_qtl.route_c_benchmark",
        ]
        completed = _without_pysam(
            f"""
            import importlib
            for name in {modules!r}:
                importlib.import_module(name)
            print("imported", len({modules!r}))
            """
        )
        self.assertEqual(0, completed.returncode, completed.stderr)
        self.assertIn(f"imported {len(modules)}", completed.stdout)

    def test_help_version_and_doctor_still_work(self) -> None:
        for arguments in (["--help"], ["--version"], ["doctor", "--json"]):
            completed = _without_pysam(
                f"""
                from candidate_qtl.cli import run
                raise SystemExit(run({arguments!r}))
                """
            )
            self.assertEqual(0, completed.returncode, f"{arguments}: {completed.stderr}")

    def test_doctor_reports_the_split(self) -> None:
        completed = _without_pysam(
            """
            from candidate_qtl.core import doctor
            import json
            report = doctor()
            print(json.dumps({
                "pysam_available": report["pysam_available"],
                "dependent": report["pysam_dependent_commands"],
                "module_flag": report["python_modules"]["pysam"],
            }))
            """
        )
        self.assertEqual(0, completed.returncode, completed.stderr)
        report = json.loads(completed.stdout.strip().splitlines()[-1])
        self.assertFalse(report["pysam_available"])
        self.assertFalse(report["module_flag"])
        self.assertIn("extract", report["dependent"])


class RunWithoutPysamTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_table_only_command_completes(self) -> None:
        ledger = self.root / "ledger.tsv"
        ledger.write_text(
            "gene_id\tortholog_ids\tfunction_evidence_source\tparental_variant\tvariant_class\tf2_cosegregation_status\n"
            "GeneA\tAT1G00010\tPMID:1\tc.1delA\tframeshift\tstrong\n",
            encoding="utf-8",
        )
        out = self.root / "ranked.tsv"
        completed = _without_pysam(
            f"""
            from candidate_qtl.cli import run
            raise SystemExit(run(["rank", "--ledger", {str(ledger)!r}, "--out", {str(out)!r}]))
            """
        )
        self.assertEqual(0, completed.returncode, completed.stderr)
        self.assertTrue(out.exists())
        self.assertIn("GeneA", out.read_text(encoding="utf-8"))

    def test_sequence_command_fails_with_an_actionable_message(self) -> None:
        completed = _without_pysam(
            f"""
            from candidate_qtl.cli import run
            raise SystemExit(run([
                "extract", "--fasta", {str(TUTORIAL / "genome.fa")!r},
                "--gff3", {str(TUTORIAL / "genes.gff3")!r},
                "--interval", "chr1:1-900", "--out-dir", {str(self.root / "out")!r},
            ]))
            """
        )
        # Expected failure (1), not an internal defect (2), and not a bare
        # ModuleNotFoundError the caller has to interpret.
        self.assertEqual(1, completed.returncode, completed.stderr)
        self.assertIn("needs pysam", completed.stderr)
        self.assertIn("candidate-qtl.ps1", completed.stderr)
        self.assertNotIn("ModuleNotFoundError", completed.stderr)
        self.assertNotIn("Traceback", completed.stderr)


if __name__ == "__main__":
    unittest.main()
