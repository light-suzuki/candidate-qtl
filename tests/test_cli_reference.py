from __future__ import annotations

from pathlib import Path
import unittest

from candidate_qtl.cli_reference import render_cli_reference


class CliReferenceTests(unittest.TestCase):
    def test_committed_references_match_each_command_surface(self) -> None:
        docs = Path(__file__).resolve().parents[1] / "docs"
        self.assertEqual((docs / "CLI-REFERENCE.md").read_text(encoding="utf-8"), render_cli_reference("core"))
        self.assertEqual((docs / "CLI-BENCH-REFERENCE.md").read_text(encoding="utf-8"), render_cli_reference("all"))
        former = (docs / "BENCHMARK-CLI-REFERENCE.md").read_text(encoding="utf-8")
        self.assertIn("(CLI-BENCH-REFERENCE.md)", former)
        self.assertNotIn("usage: candidate-qtl", former)


if __name__ == "__main__":
    unittest.main()
