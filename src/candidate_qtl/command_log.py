"""Canonical command-log schema and writer (issue #381).

The T11 benchmark run produced 98 ``command_log.tsv`` files in 12 distinct
legacy header shapes, three timestamp notations, and two missing-value
markers (``-`` and empty cells).  This module is the single source of truth
for the canonical form that new runs must write; the legacy reader/audit for
the historical files lives in :mod:`candidate_qtl.command_log_validation` so
the writer never has to interpret a legacy file to keep appending to one.

Canonical contract (``docs/COMMAND-LOG.ja.md``):

- one row per executed stage, columns fixed by :data:`CANONICAL_COLUMNS`
- ``argv_json`` is the structured argv (a JSON array), secrets redacted
- timestamps are UTC, aware, ``YYYY-MM-DDTHH:MM:SSZ``
- a stage that was not executed is recorded as a ``skipped`` row when the
  runner knows about it; a stage that is simply absent is an audit finding,
  never silently padded (the validator reports missingness, it does not
  fabricate completeness)
"""

from __future__ import annotations

import csv
import json
import os
import re
import shlex
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator

from . import __version__

# Canonical schema version. Mirrors the repo's JSON manifest convention of a
# dotted string (e.g. "1.0") so a future breaking schema change can be
# detected by the validator without guessing.
CANONICAL_SCHEMA_VERSION = "1.0"

# Column order is fixed; the validator requires an exact header match.
CANONICAL_COLUMNS: list[str] = [
    "schema_version",
    "case_id",
    "stage",
    "command",
    "argv_json",
    "status",
    "exit_code",
    "started_at_utc",
    "completed_at_utc",
    "wall_time_seconds",
    "working_directory",
    "tool_version",
    "notes",
]

# Canonical stage names. Legacy spellings map onto these
# (command_log_validation.STAGE_NORMALIZATION).
CANONICAL_STAGES = ("extract", "dry_run", "full_execute")

STATUS_VALUES = ("success", "failure", "skipped")

# Allowed |wall_time_seconds - (completed - started)| drift before the strict
# validator warns that the two measurements disagree.
WALL_TIME_TOLERANCE_SECONDS = 30.0

# Empty values in canonical rows: any column may be left empty by the
# producer (e.g. exit_code unknown, wall time not recorded, no notes), and the
# reader/validator preserves that missingness instead of inventing a value.
MISSING_MARKERS = ("", "-")

_ISO_UTC_RE = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(\.\d+)?(Z|[+-]\d{2}:\d{2})$")

# --------------------------------------------------------------------------- #
# time formatting
# --------------------------------------------------------------------------- #


def format_utc(dt_or_none: datetime | None) -> str:
    """Format an aware datetime as ``YYYY-MM-DDTHH:MM:SSZ`` (UTC).

    ``None`` formats as the empty missing marker. A naive datetime is a
    caller defect and raises (the canonical form must not record a zone-less
    wall-clock time as if it were UTC).
    """
    if dt_or_none is None:
        return ""
    if dt_or_none.tzinfo is None:
        raise ValueError("naive datetime cannot be written to a canonical command log; provide an aware UTC value")
    return dt_or_none.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def parse_utc(value: str | None) -> datetime | None:
    """Parse a canonical (ISO, UTC) timestamp into an aware datetime.

    Returns ``None`` for missing or unparseable input and for naive timestamps
    (a canonical row must be zone-qualified; zone-less values are reported by
    the validator instead of being assumed UTC).
    """
    if value is None:
        return None
    text = str(value).strip()
    if not text or text in MISSING_MARKERS:
        return None
    if not _ISO_UTC_RE.match(text):
        return None
    try:
        parsed = datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError:
        return None
    if parsed.tzinfo is None:
        return None
    return parsed.astimezone(timezone.utc)


# --------------------------------------------------------------------------- #
# secret redaction
# --------------------------------------------------------------------------- #

# Option flags whose next token is a secret value (``--api-key secret``).
_SECRET_FLAGS = frozenset({
    "api-key", "apikey", "api_key",
    "token", "auth-token", "auth_token", "access-token", "access_token",
    "password", "passwd",
    "secret", "client-secret", "client_secret", "client-id-secret",
    "bearer", "private-key", "private_key",
})

_SECRET_FLAG_RE = re.compile(r"^--([a-z][a-z0-9_-]*)$", re.IGNORECASE)

# Built-in single-token patterns (``--api-key=secret``, ``api_key=secret``,
# ``TOKEN=secret``).  Custom patterns supplied to :func:`redact_argv` are
# merged on top (custom wins on name conflict); secrets are never left
# unredacted because a caller forgot to opt in.
DEFAULT_REDACT_PATTERNS: dict[str, str] = {
    "api_key": r"(?i)(api[_-]?key)=(\S+)",
    "token": r"(?i)\b(token)=(\S+)",
    "password": r"(?i)\b(passw(?:or)?d)=(\S+)",
    "secret": r"(?i)\b(secret)=(\S+)",
}


def redact_argv(argv: list[str], patterns: dict[str, str] | None = None) -> list[str]:
    """Return a copy of ``argv`` with secret values replaced.

    Applies two layers:

    1. flag/value pairs: ``--api-key VALUE`` (and the other
       :data:`_SECRET_FLAGS` names) replace the value token with
       ``<redacted:NAME>``; ``--api-key=VALUE`` becomes ``--api-key=<redacted:api-key>``.
    2. ``key=value`` assignments inside a token. Custom ``patterns`` are
       applied before the built-in defaults (a provider-specific rule such as
       ``OPENAI_API_KEY=...`` wins over the generic ``api_key`` rule so the
       whole assignment is replaced, not a substring); each match is replaced
       by ``<redacted:NAME>``.
    """
    compiled = [(name, re.compile(regex)) for name, regex in (patterns or {}).items()]
    compiled.extend((name, re.compile(regex)) for name, regex in DEFAULT_REDACT_PATTERNS.items())

    def _redact_token(token: str) -> str:
        for name, regex in compiled:
            if regex.search(token):
                token = regex.sub(f"<redacted:{name}>", token)
        return token

    result: list[str] = []
    index = 0
    while index < len(argv):
        token = argv[index]
        match = _SECRET_FLAG_RE.match(token)
        if match and match.group(1).lower() in _SECRET_FLAGS:
            result.append(token)
            if index + 1 < len(argv):
                result.append(f"<redacted:{match.group(1).lower()}>")
                index += 2
                continue
        result.append(_redact_token(token))
        index += 1
    return result


# --------------------------------------------------------------------------- #
# writer
# --------------------------------------------------------------------------- #


class CommandLogWriter:
    """Append-only canonical command-log writer with header verification.

    - the header is written only when the file is created
    - an existing file must already carry exactly :data:`CANONICAL_COLUMNS`,
      otherwise appending is refused (a legacy/mixed file is never silently
      extended with canonical rows)
    - every append is flushed and fsync'ed before returning, so a later crash
      cannot make a recorded run look successful
    - writer failures propagate as exceptions (never a fake-success result)
    """

    def __init__(
        self,
        path: Path | str,
        *,
        tool_version: str = __version__,
        redact_patterns: dict[str, str] | None = None,
    ) -> None:
        self.path = Path(path)
        self.tool_version = tool_version
        self.redact_patterns = redact_patterns
        existing = self.path.exists() and self.path.stat().st_size > 0
        if existing:
            self._verify_existing_header()
        self._file = self.path.open("a", encoding="utf-8", newline="")
        self._writer = csv.writer(self._file, delimiter="\t", lineterminator="\n")
        if not existing:
            self._writer.writerow(CANONICAL_COLUMNS)
            self._flush()

    def _verify_existing_header(self) -> None:
        try:
            with self.path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
                first = handle.readline().rstrip("\r\n")
        except OSError as error:
            raise OSError(f"cannot read existing command log {self.path}: {error}") from error
        if first != "\t".join(CANONICAL_COLUMNS):
            raise ValueError(
                f"refusing to append to {self.path}: existing header does not match the canonical "
                f"command_log schema v{CANONICAL_SCHEMA_VERSION} (found {first!r})")

    def _flush(self) -> None:
        self._file.flush()
        os.fsync(self._file.fileno())

    def append(
        self,
        *,
        case_id: str,
        stage: str,
        command: str,
        argv: list[str],
        status: str,
        exit_code: int | None,
        started_at_utc: datetime | None,
        completed_at_utc: datetime | None,
        wall_time_seconds: float | None,
        working_directory: str,
        notes: str = "",
    ) -> None:
        """Append one canonical record. Values are written as given; empty
        cells are the missing marker and are preserved on read-back."""
        if not case_id:
            raise ValueError("case_id must not be empty in a canonical command log")
        if stage not in CANONICAL_STAGES:
            raise ValueError(f"unknown stage {stage!r}; canonical stages are {CANONICAL_STAGES}")
        if status not in STATUS_VALUES:
            raise ValueError(f"unknown status {status!r}; status values are {STATUS_VALUES}")
        # The human-readable command summary may echo secret flag values; it
        # goes through the same redaction as argv_json so no secret is
        # recorded in either column.
        try:
            command_tokens = shlex.split(command, posix=True)
        except ValueError:
            # A malformed shell summary must never bypass redaction.  The
            # structured argv is authoritative and already tokenized.
            command_tokens = list(argv)
        redacted_command = " ".join(redact_argv(command_tokens, self.redact_patterns))
        row = [
            CANONICAL_SCHEMA_VERSION,
            case_id,
            stage,
            redacted_command,
            json.dumps(redact_argv(argv, self.redact_patterns), ensure_ascii=False),
            status,
            "" if exit_code is None else str(exit_code),
            format_utc(started_at_utc),
            format_utc(completed_at_utc),
            "" if wall_time_seconds is None else repr(wall_time_seconds),
            working_directory,
            self.tool_version,
            notes,
        ]
        self._writer.writerow(row)
        self._flush()

    def close(self) -> None:
        if not self._file.closed:
            self._file.close()

    def __enter__(self) -> CommandLogWriter:
        return self

    def __exit__(self, *exc_info: object) -> None:
        self.close()


@contextmanager
def timing() -> Iterator[object]:
    """Context manager that yields a callable returning the elapsed seconds.

    The callable is live while the block runs and keeps its final value
    afterwards (3-decimal precision, monotonic clock):

        with timing() as elapsed:
            run_something()
        wall = elapsed()   # rounded seconds
    """
    started = time.monotonic()

    def elapsed() -> float:
        return round(time.monotonic() - started, 3)

    yield elapsed


# --------------------------------------------------------------------------- #
# reader
# --------------------------------------------------------------------------- #


def read_command_log(path: Path | str) -> tuple[list[str], list[dict[str, str]]]:
    """Read a canonical command log without normalizing any value.

    Returns ``(header, rows)`` where every row is a dict keyed by the
    canonical column names and missing cells stay empty strings (missingness
    is preserved, never padded). Raises ``ValueError`` when the file is not a
    canonical command log (wrong header), and ``OSError`` when it cannot be
    read.
    """
    log_path = Path(path)
    with log_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = [row for row in reader]
    if not rows:
        raise ValueError(f"empty command log: {log_path}")
    header = [cell.strip() for cell in rows[0]]
    if header != CANONICAL_COLUMNS:
        raise ValueError(
            f"{log_path} is not a canonical command log: header {header!r} "
            f"does not match {CANONICAL_COLUMNS!r}")
    records: list[dict[str, str]] = []
    for line_number, row in enumerate(rows[1:], 2):
        if len(row) > len(CANONICAL_COLUMNS):
            raise ValueError(
                f"{log_path} is not a canonical command log: row {line_number} has "
                f"{len(row)} cells; expected at most {len(CANONICAL_COLUMNS)}"
            )
        padded = row + [""] * (len(CANONICAL_COLUMNS) - len(row))
        records.append(dict(zip(CANONICAL_COLUMNS, padded)))
    return header, records
