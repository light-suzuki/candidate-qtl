from __future__ import annotations

import tomllib
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class PackagingBoundaryTests(unittest.TestCase):
    def test_pysam_is_alignment_extra_not_default_dependency(self) -> None:
        data = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
        dependencies = data["project"]["dependencies"]
        alignment = data["project"]["optional-dependencies"]["alignment"]
        self.assertFalse(any(item.lower().startswith("pysam") for item in dependencies))
        self.assertTrue(any(item.lower().startswith("pysam") for item in alignment))

    @unittest.skip("public release: candidate-qtl-bench entry point is not shipped")
    def test_production_and_benchmark_entrypoints_are_distinct(self) -> None:
        scripts = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]["scripts"]
        self.assertEqual("candidate_qtl.cli:run", scripts["candidate-qtl"])
        self.assertEqual("candidate_qtl.bench_cli:run", scripts["candidate-qtl-bench"])


if __name__ == "__main__":
    unittest.main()
