#!/usr/bin/env python3
"""Run the canonical feature-reachability gate (issue #373).

Executes every integrated entry of the committed reachability registry
(tests/fixtures/reachability/registry.json) as a real subprocess
(``python -m candidate_qtl <command> ...``) and asserts the declared
artifacts, manifest keys, fixture immutability and no-rank-artifacts policy.

Modes:

* ``--mode source`` -- cwd is the repo root; the subprocess python resolves
  the package the way the editable CI install does.
* ``--mode installed`` -- first proves via assert_no_source_fallback that the
  given interpreter imports candidate_qtl from site-packages (not from the
  source checkout), then runs every entry with cwd outside the repo root.

The tooling layer is test infrastructure: this script may run under a python
without candidate_qtl installed (or with an editable install pointing at a
different checkout), so the reachability module itself is imported from the
source checkout when it is not already importable.  The integration proof is
the subprocess run of the *installed* CLI, which this import never shadows.
"""

from __future__ import annotations

import argparse
import importlib.util
import shutil
import sys
import tempfile
import traceback
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REGISTRY = REPO_ROOT / "tests" / "fixtures" / "reachability" / "registry.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the canonical feature-reachability gate.")
    parser.add_argument("--registry", default=str(DEFAULT_REGISTRY), help="Registry JSON (default: tests/fixtures/reachability/registry.json).")
    parser.add_argument("--repo-root", default=str(REPO_ROOT), help="Repository root; fixtures resolve against it (default: repository root).")
    parser.add_argument("--mode", required=True, choices=["source", "installed"],
                        help="source: run with cwd=repo root; installed: prove site-packages import, then run outside the repo.")
    parser.add_argument("--python", default=sys.executable, help="Interpreter that runs `python -m candidate_qtl` (installed mode: point at the clean venv, e.g. /tmp/cleanenv/bin/python).")
    parser.add_argument("--out-root", default="", help="Output directory for the gate runs (default: fresh temp dir, removed on success unless --keep-out).")
    parser.add_argument("--keep-out", action="store_true", help="Keep the auto-created output directory.")
    parser.add_argument("--timeout", type=int, default=300, help="Per-entry subprocess timeout in seconds (default: 300).")
    args = parser.parse_args()

    repo_root = Path(args.repo_root).resolve()
    # Import the reachability module from the source checkout only when it is
    # not already importable (see module docstring).  The proof of clean
    # installation is the SUBPROCESS import probe, not this import.
    if importlib.util.find_spec("candidate_qtl.reachability") is None:
        src = repo_root / "src"
        if str(src) not in sys.path:
            sys.path.insert(0, str(src))

    from candidate_qtl import reachability  # noqa: E402

    out_root = Path(args.out_root).resolve() if args.out_root else Path(tempfile.mkdtemp(prefix="reachability-gate-"))
    out_root.mkdir(parents=True, exist_ok=True)
    auto_created = not args.out_root

    try:
        if args.mode == "installed":
            try:
                resolved = reachability.assert_no_source_fallback(args.python, repo_root)
                print(f"clean-install import resolved to {resolved}")
            except reachability.ReachabilityError as error:
                print(f"FAIL <clean-install>: {error}", file=sys.stderr)
                return 1
        result = reachability.run_gate(
            args.registry, repo_root, args.mode, args.python,
            out_root=out_root, timeout=args.timeout,
        )
    except reachability.ReachabilityError as error:
        print(f"FAIL <registry>: {error}", file=sys.stderr)
        return 1
    except Exception:
        print("internal error: unexpected failure in the reachability gate", file=sys.stderr)
        traceback.print_exc()
        return 1

    failed = False
    for entry in result["entries"]:
        feature_id = entry["feature_id"]
        if entry.get("error"):
            failed = True
            print(f"FAIL {feature_id}: {entry['error']}")
        else:
            artifacts = len(entry["artifacts"])
            manifest_keys = len(entry["manifest_keys"])
            print(f"PASS {feature_id} (exit {entry['exit_code']}, {artifacts} artifacts, {manifest_keys} manifest keys)")
    print(f"summary: integrated {result['integrated_count']}, core-only {result['core_only_count']}, valid={result['valid']}")
    if result["core_only_count"]:
        print("core-only entries are not counted as integrated coverage")
    if result["valid"] and not failed:
        if auto_created and not args.keep_out:
            shutil.rmtree(out_root, ignore_errors=True)
        print(f"reachability gate passed ({result['integrated_count']} integrated entries)")
        return 0
    print(f"reachability gate failed; outputs kept in {result['out_root']}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
