"""Guards #266: a manifest-stamped version must correspond to an installable
release.

`pyproject.toml` and `src/candidate_qtl/__init__.py` used to carry two
independent version literals (0.21.0 vs 0.22.0) that drifted apart. Manifests
stamp `__version__`, so a frozen benchmark artifact could record a version no
wheel built from the same tree would ever produce. The fix removes the
duplication instead of just re-syncing the numbers: `__init__.py` is the only
literal, and `pyproject.toml` resolves its `version` field dynamically from
that attribute via `[tool.setuptools.dynamic]`. These tests fail if either
half of that wiring regresses -- a hardcoded `version` reappearing in
`[project]`, or the dynamic declaration pointing somewhere else.
"""

from __future__ import annotations

import re
import tomllib
import unittest
from pathlib import Path

import candidate_qtl

REPO = Path(__file__).resolve().parents[1]


class VersionSingleSourceTests(unittest.TestCase):
    def test_init_version_is_a_release_looking_string(self) -> None:
        self.assertRegex(candidate_qtl.__version__, r"^\d+\.\d+\.\d+$")

    def test_pyproject_does_not_hardcode_a_second_version_literal(self) -> None:
        config = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))
        project = config["project"]
        self.assertNotIn(
            "version",
            project,
            "a static [project].version reintroduces the #266 drift; "
            "declare it in [project].dynamic instead",
        )
        self.assertIn("version", project.get("dynamic", []))

    def test_pyproject_dynamic_version_points_at_the_init_literal(self) -> None:
        config = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))
        dynamic_version = config["tool"]["setuptools"]["dynamic"]["version"]
        self.assertEqual({"attr": "candidate_qtl.__version__"}, dynamic_version)

    def test_setuptools_resolves_the_same_version_when_available(self) -> None:
        # Belt-and-braces: where setuptools is installed (build/dev environments),
        # actually run its resolver instead of just checking the declaration shape.
        try:
            from setuptools.config.pyprojecttoml import read_configuration
        except ImportError:
            self.skipTest("setuptools is not installed in this environment")
        config = read_configuration(str(REPO / "pyproject.toml"))
        self.assertEqual(candidate_qtl.__version__, config["project"]["version"])

    def test_readmes_report_the_same_version_as_init(self) -> None:
        # The README "Current status" line is prose, not a build input, so it
        # cannot be wired to the single source of truth mechanically -- but it
        # must not silently go stale either.
        pattern = re.compile(r"Version: `([^`]+)`")
        for name in ("README.md", "README.ja.md"):
            text = (REPO / name).read_text(encoding="utf-8")
            match = pattern.search(text)
            self.assertIsNotNone(match, f"{name} is missing a 'Version: `...`' status line")
            self.assertEqual(candidate_qtl.__version__, match.group(1), name)

    def test_citation_cff_reports_the_same_version_as_init(self) -> None:
        cff = (REPO / "CITATION.cff").read_text(encoding="utf-8")
        match = re.search(r"^version:\s*[\"']?([^\"'\s]+)[\"']?\s*$", cff, re.MULTILINE)
        self.assertIsNotNone(match, "CITATION.cff is missing a version field")
        self.assertEqual(candidate_qtl.__version__, match.group(1))


if __name__ == "__main__":
    unittest.main()
