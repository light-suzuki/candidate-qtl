"""Deterministic CLI reference generation for documentation drift checks."""

from __future__ import annotations

import os

from .cli import build_parser


def render_cli_reference(surface: str = "core") -> str:
    previous_columns = os.environ.get("COLUMNS")
    os.environ["COLUMNS"] = "120"
    try:
        prog = "candidate-qtl" if surface == "core" else "candidate-qtl-bench"
        parser = build_parser(surface=surface, prog=prog)
        subparser_action = next(action for action in parser._subparsers._group_actions if hasattr(action, "choices"))
        sections = [
            "# CLI reference",
            "",
            "Generated from the argparse registry. Run `python scripts/generate_cli_reference.py` after an intentional CLI change.",
            "",
            f"## {prog}",
            "",
            "```text",
            parser.format_help().rstrip(),
            "```",
        ]
        for command, command_parser in sorted(subparser_action.choices.items()):
            sections.extend(
                [
                    "",
                    f"## {prog} {command}",
                    "",
                    "```text",
                    command_parser.format_help().rstrip(),
                    "```",
                ]
            )
        return "\n".join(sections) + "\n"
    finally:
        if previous_columns is None:
            os.environ.pop("COLUMNS", None)
        else:
            os.environ["COLUMNS"] = previous_columns
