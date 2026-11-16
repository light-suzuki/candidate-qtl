from __future__ import annotations

from pathlib import Path
import unittest

from candidate_qtl.cli import build_parser


def _registered_command_names() -> list[str]:
    parser = build_parser()
    subparsers = next(
        action for action in parser._subparsers._group_actions if hasattr(action, "choices")
    )
    return sorted(subparsers.choices)


def _first_cell(line: str) -> str | None:
    stripped = line.strip()
    if not stripped.startswith("|"):
        return None
    parts = [part.strip() for part in stripped.strip("|").split("|")]
    return parts[0] if parts else None


class CliDocsCompletenessTests(unittest.TestCase):
    DOCS_CLI_JA = Path(__file__).resolve().parents[1] / "docs" / "CLI.ja.md"

    def test_every_registered_command_appears_in_cli_docs(self) -> None:
        cells = {
            cell
            for line in self.DOCS_CLI_JA.read_text(encoding="utf-8").splitlines()
            if (cell := _first_cell(line)) is not None
        }
        missing = [cmd for cmd in _registered_command_names() if cmd not in cells]
        self.assertEqual(missing, [])


if __name__ == "__main__":
    unittest.main()
