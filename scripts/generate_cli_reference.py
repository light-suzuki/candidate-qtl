#!/usr/bin/env python3
"""Generate or verify docs/CLI-REFERENCE.md from the command registry."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from candidate_qtl.cli_reference import render_cli_reference


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="Fail when the committed reference is stale.")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    targets = {
        root / "docs" / "CLI-REFERENCE.md": render_cli_reference("core"),
        root / "docs" / "CLI-BENCH-REFERENCE.md": render_cli_reference("all"),
        # Preserve the former filename for links without duplicating full help.
        root / "docs" / "BENCHMARK-CLI-REFERENCE.md": (
            "# Research and benchmark CLI reference\n\n"
            "The generated reference is [CLI-BENCH-REFERENCE.md](CLI-BENCH-REFERENCE.md).\n\n"
            "This former filename is retained for existing links. Run `python scripts/generate_cli_reference.py` "
            "to update the canonical references and this pointer.\n"
        ),
    }
    if args.check:
        stale = [target for target, rendered in targets.items() if not target.exists() or target.read_text(encoding="utf-8") != rendered]
        if stale:
            for target in stale:
                print(f"CLI reference is stale: {target}", file=sys.stderr)
            return 1
        for target in targets:
            print(f"CLI reference is current: {target}")
        return 0
    for target, rendered in targets.items():
        target.write_text(rendered, encoding="utf-8")
        print(target)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
