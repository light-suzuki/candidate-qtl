"""Deterministic content-addressable cache for expensive stage outputs.

Entries are keyed by a sha256 of every input that can change the result
(input file hashes, tool, provider, provider_version, params), so a
deterministic engine reuses outputs only when the inputs match exactly.
Reads are fail-closed: schema, staleness, and payload corruption reject the
entry and count a rejection instead of raising, so the engine recomputes
instead of trusting a possibly-bad entry.
"""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from pathlib import Path

from .core import sha256

CACHE_SCHEMA_VERSION = "1.0"
_PAYLOAD_NAME = "payload"


def _unlink_quiet(path: str) -> None:
    try:
        os.unlink(path)
    except FileNotFoundError:
        pass


class ContentCache:
    """Content-addressed cache on disk: one metadata JSON plus one payload file per key."""

    def __init__(self, dir: str | Path, enabled: bool = True) -> None:
        self._dir = Path(dir)
        self._enabled = enabled
        self._hit_count = 0
        self._miss_count = 0
        self._corrupt_reject_count = 0
        self._stale_reject_count = 0

    def key(self, input_hashes: dict[str, str], tool: str, provider: str, provider_version: str, params: dict[str, object] | None) -> str:
        """Content key: sha256 of every input that can change the result."""
        serialized = json.dumps({
            "input_hashes": input_hashes, "tool": tool, "provider": provider,
            "provider_version": provider_version, "params": params,
        }, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
        return hashlib.sha256(serialized).hexdigest()

    def _meta_path(self, key: str) -> Path:
        return self._dir / f"{key}.json"

    def _payload_path(self, key: str) -> Path:
        return self._dir / f"{key}.{_PAYLOAD_NAME}"

    def get(self, input_hashes, tool, provider, provider_version, params) -> tuple[bool, Path | None]:
        """Return (True, payload path) on a valid hit, else (False, None); never raises."""
        if not self._enabled:
            return False, None
        key = self.key(input_hashes, tool, provider, provider_version, params)
        meta_path = self._meta_path(key)
        if not meta_path.is_file():
            self._miss_count += 1
            return False, None
        try:
            metadata = json.loads(meta_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            self._corrupt_reject_count += 1
            self._miss_count += 1
            return False, None
        if not isinstance(metadata, dict):
            self._corrupt_reject_count += 1
            self._miss_count += 1
            return False, None
        if metadata.get("schema_version") != CACHE_SCHEMA_VERSION:
            self._stale_reject_count += 1
            self._miss_count += 1
            return False, None
        if (metadata.get("tool") != tool
                or metadata.get("provider") != provider
                or metadata.get("provider_version") != provider_version
                or metadata.get("input_hashes") != input_hashes
                or metadata.get("params") != params):
            self._stale_reject_count += 1
            return False, None
        payload_path = self._payload_path(key)
        try:
            payload_digest = sha256(payload_path)
        except OSError:
            self._corrupt_reject_count += 1
            self._miss_count += 1
            return False, None
        if payload_digest != metadata.get("payload_sha256"):
            self._corrupt_reject_count += 1
            self._miss_count += 1
            return False, None
        self._hit_count += 1
        return True, payload_path

    def put(self, input_hashes, tool, provider, provider_version, params, payload_path) -> str:
        """Store a payload under its content key (payload replaced before metadata is published)."""
        key = self.key(input_hashes, tool, provider, provider_version, params)
        if not self._enabled:
            return key
        self._dir.mkdir(parents=True, exist_ok=True)
        source = Path(payload_path)
        payload_digest = sha256(source)
        payload_target = self._payload_path(key)

        payload_fd, payload_tmp = tempfile.mkstemp(dir=self._dir, prefix=f".{key}.payload.", suffix=".tmp")
        try:
            with os.fdopen(payload_fd, "wb") as handle:
                handle.write(source.read_bytes())
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(payload_tmp, payload_target)
        finally:
            _unlink_quiet(payload_tmp)

        metadata = {
            "schema_version": CACHE_SCHEMA_VERSION, "key": key, "tool": tool, "provider": provider,
            "provider_version": provider_version, "params": params, "input_hashes": input_hashes,
            "payload_sha256": payload_digest, "payload_name": _PAYLOAD_NAME,
        }
        meta_path = self._meta_path(key)
        meta_fd, meta_tmp = tempfile.mkstemp(dir=self._dir, prefix=f".{key}.meta.", suffix=".tmp")
        try:
            with os.fdopen(meta_fd, "w", encoding="utf-8") as handle:
                handle.write(json.dumps(metadata, ensure_ascii=False, indent=2, sort_keys=True) + "\n")
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(meta_tmp, meta_path)
        finally:
            _unlink_quiet(meta_tmp)
        return key

    def stats(self) -> dict[str, object]:
        """Plain JSON-serializable counters; cache_status enabled|disabled."""
        return {
            "hit_count": self._hit_count, "miss_count": self._miss_count,
            "corrupt_reject_count": self._corrupt_reject_count,
            "stale_reject_count": self._stale_reject_count,
            "cache_status": "enabled" if self._enabled else "disabled",
        }
