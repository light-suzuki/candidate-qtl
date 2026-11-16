"""Executor package commands: physical Executor/Supervisor/run-capsule separation (issue #365)."""

from __future__ import annotations

import shlex
from argparse import _SubParsersAction
from typing import Any

from ..core import runtime_path
from ..package_build import build_executor_package, detect_repo_root
from ._shared import add_json_flag


def register(subparsers: _SubParsersAction) -> None:
    sub = subparsers.add_parser(
        "build-executor-package",
        help="Build a physically separated Executor package (code + executor manifests + declared reference metadata only) with an allowlist package manifest and a run_capsule.json. Fails closed when any truth/scoring artifact or link would enter the package.",
    )
    sub.add_argument("--out-dir", required=True, help="New output directory outside the repository root.")
    sub.add_argument("--repo-root", help="Repository root; default: detected from the candidate-qtl pyproject.toml.")
    sub.add_argument("--include", action="append", help="Repo-relative glob to include; repeatable. Default: the built-in allowlist patterns.")
    sub.add_argument("--reference-assets", help="Optional TSV (path_wsl required; resource_key/provider/version/retrieved_at_utc optional) of external reference assets fingerprinted into the capsule.")
    sub.add_argument("--run-command", help="Quoted Executor command line (e.g. \"candidate-qtl audit-predict --manifest m.tsv\") recorded in the capsule as the declared execution command.")
    sub.add_argument("--full-hash", action="store_true", help="Force full SHA-256 for large reference assets.")
    add_json_flag(sub)


def _run_build_executor_package(args: Any) -> dict:
    repo_root = runtime_path(args.repo_root) if args.repo_root else detect_repo_root()
    run_command = shlex.split(args.run_command) if args.run_command else None
    result = build_executor_package(
        repo_root=repo_root,
        out_dir=runtime_path(args.out_dir),
        include_patterns=args.include,
        reference_assets=runtime_path(args.reference_assets) if args.reference_assets else None,
        full_hash=args.full_hash,
        run_command=run_command,
    )
    return {
        "command": "build-executor-package",
        "out_dir": str(result.out_dir),
        "package_manifest": str(result.manifest_path),
        "run_capsule": str(result.capsule_path),
        "file_count": result.file_count,
        "total_bytes": result.total_bytes,
        "rejected_count": len(result.rejected),
        "review_required_count": len(result.review_required),
        "truth_inputs_present": result.truth_inputs_present,
        "exit_code": 0,
    }


HANDLERS = {
    "build-executor-package": _run_build_executor_package,
}
