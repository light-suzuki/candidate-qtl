from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import traceback
import urllib.error

from . import __version__
from .commands import handler_for, register_all
from .core import ExternalDependencyError


# Exit-code contract for the process entry point (`run`).
EXIT_OK = 0
EXIT_EXPECTED_FAILURE = 1
EXIT_INTERNAL_ERROR = 2

# Failures the caller can act on: rejected input, a missing or unreadable file, a
# refused gate, or an external tool that returned non-zero. The codebase raises
# ValueError for rejected input at ~400 sites and RuntimeError when a spawned
# tool fails, so these are the workflow's own "your input is wrong" channel and
# get one actionable line instead of a traceback. Anything else is a defect in
# this package and keeps its traceback.
EXPECTED_ERRORS = (
    ValueError,          # includes schemas.SchemaValidationError
    FileNotFoundError,
    IsADirectoryError,
    NotADirectoryError,
    PermissionError,
    RuntimeError,
)

# Failures the environment causes, not this package: a network outage, a hung
# external tool, or an expired timeout. They exit like expected input errors
# (one actionable line, no traceback) but the message names the dependency.
ENVIRONMENT_ERRORS = (
    ExternalDependencyError,
    urllib.error.URLError,   # includes urllib.error.HTTPError
    TimeoutError,
    subprocess.TimeoutExpired,
)

TRACEBACK_ENV = "CANDIDATE_QTL_TRACEBACK"


def build_parser(surface: str = "all", prog: str | None = None) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=prog or "candidate-qtl", description="Evidence-preserving plant QTL candidate triage and experiment decision-support workflow.")
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument(
        "--traceback",
        action="store_true",
        help=f"Show the Python traceback for expected input errors too (or set {TRACEBACK_ENV}=1).",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    register_all(subparsers, surface=surface)
    return parser


def emit(result: dict[str, object], as_json: bool) -> None:
    if as_json:
        print(json.dumps(result, ensure_ascii=False, indent=2))
    else:
        for key, value in result.items():
            print(f"{key}: {value}")


def traceback_requested(argv: list[str] | None = None) -> bool:
    """Report whether the caller asked to keep tracebacks on expected errors.

    Read from raw argv rather than the parsed namespace because the decision is
    made in ``run`` after ``main`` has already raised.
    """
    arguments = sys.argv[1:] if argv is None else argv
    if "--traceback" in arguments:
        return True
    return os.environ.get(TRACEBACK_ENV, "").strip().lower() not in {"", "0", "false", "no"}


def _environment_kind(error: Exception) -> str | None:
    if isinstance(error, (TimeoutError, subprocess.TimeoutExpired)):
        return "timeout"
    if isinstance(error, (urllib.error.URLError, urllib.error.HTTPError)):
        return "network"
    return None


def main(argv: list[str] | None = None, *, surface: str = "all", prog: str | None = None) -> int:
    """Parse and dispatch one command, re-raising on failure.

    This is the library entry point and the one the tests drive, so a failure
    still propagates as its original exception. ``run`` is the process entry
    point that turns those exceptions into exit codes.
    """
    parser = build_parser(surface=surface, prog=prog)
    args = parser.parse_args(argv)
    try:
        handler = handler_for(args.command, surface=surface)
        if handler is None:
            raise ValueError(f"unsupported command: {args.command}")
        result = handler(args)
        emit(result, args.json)
        # Commands that deliberately record a failure instead of raising (so their
        # audit manifest survives) still signal it through exit_code.
        return int(result.get("exit_code", 0)) if isinstance(result, dict) else 0
    except Exception as error:
        kind = _environment_kind(error)
        message = f"external dependency failure (kind={kind}): {error}" if kind else str(error)
        if getattr(args, "json", False):
            print(json.dumps({"error": message, "command": args.command}, ensure_ascii=False), file=sys.stderr)
        else:
            print(f"error: {message}", file=sys.stderr)
        raise


def run(argv: list[str] | None = None) -> int:
    """Process entry point: report expected failures without a traceback.

    ``main`` has already written the one-line ``error:`` summary (or the JSON
    error object) by the time an exception arrives here, so an expected failure
    needs no further output -- only the exit code. A traceback on a mistyped
    interval buries the message that says how to fix it.
    """
    try:
        return main(argv, surface="core", prog="candidate-qtl")
    except SystemExit:
        raise
    except ENVIRONMENT_ERRORS:
        if traceback_requested(argv):
            raise
        return EXIT_EXPECTED_FAILURE
    except EXPECTED_ERRORS:
        if traceback_requested(argv):
            raise
        return EXIT_EXPECTED_FAILURE
    except Exception:
        # Unexpected: keep the traceback, and separate it from the summary line
        # main already printed so the two are not read as one message.
        print(f"internal error: this is a defect in candidate-qtl {__version__}, not an input problem.", file=sys.stderr)
        traceback.print_exc()
        return EXIT_INTERNAL_ERROR


if __name__ == "__main__":
    raise SystemExit(run())
