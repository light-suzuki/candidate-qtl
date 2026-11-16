"""Record an unexecuted claim review when a frozen request has no source text.

This is a machine coverage result, never an independent semantic review. When
source records exist the agent must run the actual review instead.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from candidate_qtl.claim_review_session import build_no_source_review as build


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prompt", type=Path, required=True)
    parser.add_argument("--route", choices=("a", "b"), required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = build(args.prompt, args.route)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
