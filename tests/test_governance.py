"""Guards #250: closed hardening issues must leave a machine-checkable
enforcement behind, not just a fix.

Each test pins the CI/config wiring that makes a historical hardening fix
stick: #191 (coverage and CLI-reference drift gates) and #266 (single-source
version). If a future change removes a gate, the corresponding test fails.
"""

from __future__ import annotations

import tomllib
import unittest
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]


def _pyproject() -> dict:
    return tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))


def _ci_steps() -> list[dict]:
    ci = yaml.safe_load((REPO / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8"))
    return ci["jobs"]["test"]["steps"]


class CoverageGateTests(unittest.TestCase):
    def test_pyproject_coverage_fail_under_gate(self) -> None:
        report = _pyproject()["tool"]["coverage"]["report"]
        self.assertIn("fail_under", report, "no [tool.coverage.report] fail_under; the #191 coverage gate is gone")
        self.assertGreaterEqual(report["fail_under"], 80, "fail_under was lowered below the #191 acceptance bar")

    def test_ci_runs_coverage_report(self) -> None:
        runs = [s.get("run", "") for s in _ci_steps()]
        self.assertTrue(
            any("python -m coverage report" in run for run in runs),
            "test job never runs `coverage report`, so fail_under cannot gate CI",
        )


class CliReferenceGateTests(unittest.TestCase):
    def test_ci_has_named_cli_reference_drift_step(self) -> None:
        steps = [s for s in _ci_steps() if s.get("name") == "Check the generated CLI reference is current"]
        self.assertEqual(len(steps), 1, "the named CLI-reference drift step is missing")
        self.assertIn(
            "generate_cli_reference.py --check",
            steps[0]["run"],
            "the drift step does not run the --check form of the generator",
        )


class VersionSingleSourceEnforcementTests(unittest.TestCase):
    def test_version_single_source_enforcement_present(self) -> None:
        config = _pyproject()
        project = config["project"]
        self.assertNotIn("version", project, "a static [project].version reintroduces the #266 drift")
        self.assertIn("version", project.get("dynamic", []), "[project].dynamic must declare version")
        self.assertEqual(
            {"attr": "candidate_qtl.__version__"},
            config["tool"]["setuptools"]["dynamic"]["version"],
            "the dynamic version must resolve from candidate_qtl.__version__",
        )
        self.assertTrue(
            (REPO / "tests" / "test_version_single_source.py").is_file(),
            "tests/test_version_single_source.py is missing; nothing guards the #266 wiring",
        )
        runs = [s.get("run", "") for s in _ci_steps()]
        self.assertTrue(
            any("run_test_partition.py --partition core" in run for run in runs),
            "test job does not run the core partition, so the #266 guard never executes in CI",
        )


if __name__ == "__main__":
    unittest.main()
