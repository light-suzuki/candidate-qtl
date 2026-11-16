"""Canonical feature-reachability gate (issue #373).

A repository-owned machine-readable registry (``tests/fixtures/reachability/
registry.json``) maps every feature that *claims* pipeline integration to its
canonical CLI entrypoint, committed fixtures, expected artifacts, manifest keys
and rank policy.  ``run_gate`` executes each integrated entry as a real
subprocess (``python -m candidate_qtl <command> ...``), asserts the declared
artifacts and manifest keys, proves fixture immutability, and proves that the
entry writes no rank artifact.  ``assert_no_source_fallback`` proves that a
clean (non-editable) install resolves ``candidate_qtl`` into site-packages,
not into the source checkout.

A feature declared ``integrated_closing_slice`` is only *called* integrated
when the gate passes for it.  A ``core_only_non_closing_slice`` declaration
records remaining canonical wiring and is never counted as integrated
coverage.

This module is deliberately dependency-free beyond the stdlib: the registry
schema is validated by hand (the repository pins no ``jsonschema``), and the
gate is hermetic (no network, no external binaries, outputs only into temp
directories).
"""

from __future__ import annotations

import fnmatch
import hashlib
import importlib
import json
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

SCHEMA_VERSION = 1

# A feature whose canonical wiring is complete: it must pass the gate before
# the issue that closes the slice may be marked closed.
STATUS_INTEGRATED = "integrated_closing_slice"
# A feature whose wiring is still incomplete: it records what remains and is
# never counted as integrated coverage by the gate.
STATUS_CORE_ONLY = "core_only_non_closing_slice"

# Every field an integrated entry must carry.  Presence alone is not enough;
# validate_registry checks each value shape separately.
REQUIRED_INTEGRATED_FIELDS = frozenset({
    "entrypoint", "fixtures", "expected_artifacts", "expected_manifest_keys",
    "rank_contribution", "rank_assertion", "clean_install_smoke",
    "cache_fingerprint_assertion", "owning_issue", "doc_reference",
})

# A registry entry is not allowed to turn an arbitrary helper into an
# "integrated" feature merely because that helper is presently exported by the
# CLI.  Adding an integration requires an explicit canonical ownership record
# here, including the JSON manifest which makes the wiring externally visible.
# Keep this list deliberately small and grow it together with the canonical
# orchestration that owns a new closing slice.
CANONICAL_INTEGRATIONS = {
    "budget-run": {
        "command": "budget-run",
        "manifest": "budget_manifest.json",
    },
    "discovery-chain-benchmark": {
        "command": "discovery-chain-benchmark",
        "manifest": "discovery_chain_manifest.json",
    },
    "family-candidate-bridge": {
        "command": "family-candidate-bridge",
        "manifest": "family_candidate_bridge_manifest.json",
    },
    "family-index": {
        "command": "family-index",
        "manifest": "family_index_manifest.json",
    },
    "mechanism-import": {
        "command": "mechanism-import",
        "manifest": "mechanism_import_manifest.json",
    },
    "mechanism-import-coverage": {
        "command": "mechanism-coverage-audit",
        "manifest": "mechanism_coverage_manifest.json",
    },
    "route-a-interval-match": {
        "command": "route-a-interval-match",
        "manifest": "route_a_interval_family_match_manifest.json",
    },
    "route-b-research-refreeze": {
        "command": "route-b-semantic-bundle",
        "manifest": "route_b_semantic_bundle_manifest.json",
    },
    "route-family-metrics-not-in-index": {
        "command": "supervised-benchmark-route-family-metrics",
        "manifest": "route_family_metrics_manifest.json",
    },
    "triple-route-e2e": {
        "command": "guided-profile",
        "manifest": "guided_profile_manifest.json",
    },
    "tier-bibliography-budget-a": {
        "command": "source-verify",
        "manifest": "source_verification_manifest.json",
    },
    "tier-bibliography-budget-b": {
        "command": "route-b-semantic-source-verify",
        "manifest": "semantic_claim_source_verification_manifest.json",
    },
}

# These are the rank-sensitive outputs of the repository's smallest canonical
# case-run fixture.  A reachability slice must not perturb them: compare a
# clean run before and after the full gate, rather than trusting a feature's
# claim that it never touches rank state.
CANONICAL_BASELINE_FIXTURES = (
    "tests/fixtures/case_run/case.yaml",
    "tests/fixtures/case_run/resources.yaml",
)
CANONICAL_BASELINE_ARTIFACTS = (
    "candidates.tsv",
    "rank/ranked_candidates.tsv",
    "rank/score_contributions.tsv",
    "stage_a/deterministic_shortlist.tsv",
)

# Strict whitelist of per-entry keys; an unknown key is a schema error so a
# typo cannot silently weaken a declaration.
ENTRY_KEY_WHITELIST = frozenset({
    "feature_id", "integration_status", "entrypoint", "fixtures",
    "expected_artifacts", "expected_manifest_keys", "expected_exit_codes",
    "rank_contribution", "rank_assertion", "clean_install_smoke",
    "cache_fingerprint_assertion", "owning_issue", "doc_reference",
    "remaining_wiring", "assert_values",
})

RANK_ASSERTION_TYPES = frozenset({
    "deterministic_and_no_rank_artifacts",
    "emits_canonical_rank_unchanged",
})

# Files an entry with rank_contribution == "none" must never write into its
# out-dir.  Matched with fnmatch against the relative path (and the basename)
# of every file under the out-dir.
RANK_IMPACTING_PATTERNS = (
    "rank/ranked_candidates.tsv",
    "score_contributions.tsv",
    "candidates.tsv",
    "flagged_candidates.tsv",
    "stage_a/**",
    "dry_rank*",
    "ranked_candidates*",
)

_FEATURE_ID_RE = re.compile(r"^[a-z0-9][a-z0-9-]*$")
# Windows drive paths, POSIX absolute paths, UNC shares and the //-prefixed
# paths argparse would otherwise accept are all forbidden in the registry.
_ABSOLUTE_PATH_RE = re.compile(r"^(?:[A-Za-z]:[\\/]|/|\\|//)")
_PLACEHOLDER_RE = re.compile(r"\{([^}]*)\}")

# candidate_qtl.commands is the source of truth for command membership, but
# importing it pulls in the pandas-based feature modules.  The gate keeps the
# import lazy so the module itself loads in a bare interpreter; the import
# fails closed (ImportError) when the feature modules' dependencies are
# absent, which is an environment error, not a registry error.
_COMMAND_NAMES: frozenset[str] | None = None


def _command_names() -> frozenset[str]:
    global _COMMAND_NAMES
    if _COMMAND_NAMES is None:
        commands = importlib.import_module("candidate_qtl.commands")
        _COMMAND_NAMES = frozenset(commands.COMMAND_NAMES)
    return _COMMAND_NAMES


class ReachabilityError(Exception):
    """A gate or validation failure bound to one feature.

    ``feature_id`` identifies the registry entry (or a gate-level pseudo-id
    such as ``<registry>`` / ``<clean-install>``); ``reason`` is the concrete
    failure.  The message carries both so log lines stay actionable.
    """

    def __init__(self, feature_id: str, reason: str) -> None:
        self.feature_id = feature_id
        self.reason = reason
        super().__init__(f"feature '{feature_id}': {reason}")


def _fixture_path(repo_root: Path, relpath: str) -> Path:
    return (repo_root / relpath).resolve()


def _is_absolute_fixture(path: str) -> bool:
    return bool(_ABSOLUTE_PATH_RE.match(path))


def _has_dotdot_component(path: str) -> bool:
    return ".." in [part for part in re.split(r"[\\/]+", path) if part]


def _output_path_error(path: object, *, allow_glob: bool) -> str | None:
    """Return the fail-closed reason for one out-dir-relative declaration."""
    if not isinstance(path, str) or not path:
        return "must be a non-empty string"
    if "\\" in path:
        return "must use forward slashes"
    if _is_absolute_fixture(path):
        return "must be out-dir-relative, not absolute"
    components = path.split("/")
    if any(component in ("", ".", "..") for component in components):
        return "must not contain empty, . or .. path components"
    if not allow_glob and any(character in path for character in "*?["):
        return "must name one declared JSON artifact, not a glob"
    return None


def _is_declared_json_artifact(artifact: str, artifacts: list[object]) -> bool:
    """Whether ``artifact`` is an exact declared JSON output, never a glob."""
    return (
        artifact.endswith(".json")
        and artifact in artifacts
        and not any(character in artifact for character in "*?[")
    )


def _dotpath_get(data: Any, dot_path: str) -> tuple[bool, Any]:
    """Resolve ``a.b.0.c`` into a JSON structure.

    Returns ``(found, value)``; list segments are resolved by integer index.
    """
    current: Any = data
    for part in dot_path.split("."):
        if isinstance(current, list):
            if not part.isdigit():
                return False, None
            index = int(part)
            if index >= len(current):
                return False, None
            current = current[index]
        elif isinstance(current, dict):
            if part not in current:
                return False, None
            current = current[part]
        else:
            return False, None
    return True, current


def _emitted_manifests(out_dir: Path) -> list[Path]:
    return sorted(path for path in out_dir.rglob("*.json") if path.is_file())


def validate_registry(data: dict, repo_root: Path | None = None) -> list[str]:
    """Validate the registry data structure; returns error strings (empty = valid).

    The rules are fail-closed: duplicates, path traversal, absolute fixture
    paths, unknown statuses, unknown keys and missing required fields are
    errors, never warnings.  When ``repo_root`` is given, every fixture path
    must additionally exist under it.
    """
    errors: list[str] = []
    if not isinstance(data, dict):
        return ["registry must be a JSON object"]
    if data.get("schema_version") != SCHEMA_VERSION:
        errors.append(f"schema_version must be {SCHEMA_VERSION}")
    entries = data.get("entries")
    if not isinstance(entries, list):
        errors.append("entries must be a list")
        return errors

    feature_ids: set[str] = set()
    for entry in entries:
        if not isinstance(entry, dict):
            errors.append("each entry must be a JSON object")
            continue
        feature_id = entry.get("feature_id")
        if not isinstance(feature_id, str) or not _FEATURE_ID_RE.fullmatch(feature_id):
            errors.append(f"entry {feature_id!r}: feature_id must match ^[a-z0-9][a-z0-9-]*$")
            continue
        if feature_id in feature_ids:
            errors.append(f"entry {feature_id!r}: duplicate feature_id")
        feature_ids.add(feature_id)

        unknown = sorted(set(entry) - ENTRY_KEY_WHITELIST)
        if unknown:
            errors.append(f"entry {feature_id!r}: unknown key(s): {', '.join(unknown)}")

        status = entry.get("integration_status")
        if status not in (STATUS_INTEGRATED, STATUS_CORE_ONLY):
            errors.append(f"entry {feature_id!r}: unknown integration_status {status!r}")
            continue
        if status == STATUS_CORE_ONLY:
            remaining = entry.get("remaining_wiring")
            doc_reference = entry.get("doc_reference")
            if not isinstance(remaining, str) or not remaining.strip():
                errors.append(f"entry {feature_id!r}: core_only entry requires a non-empty remaining_wiring")
            if not isinstance(doc_reference, str) or not doc_reference.strip():
                errors.append(f"entry {feature_id!r}: core_only entry requires a non-empty doc_reference")
            for forbidden in ("entrypoint", "expected_artifacts", "expected_manifest_keys"):
                if forbidden in entry:
                    errors.append(f"entry {feature_id!r}: core_only entry must not carry {forbidden}")
            continue

        for field in sorted(REQUIRED_INTEGRATED_FIELDS):
            if field not in entry:
                errors.append(f"entry {feature_id!r}: integrated entry missing required field {field!r}")

        entrypoint = entry.get("entrypoint", {})
        if isinstance(entrypoint, dict):
            command = entrypoint.get("command")
            surface = entrypoint.get("surface", "core")
            if surface not in ("core", "benchmark"):
                errors.append(f"entry {feature_id!r}: entrypoint.surface must be 'core' or 'benchmark'")
            if command not in _command_names():
                errors.append(f"entry {feature_id!r}: unknown command {command!r} (not in candidate_qtl.commands.COMMAND_NAMES)")
            canonical = CANONICAL_INTEGRATIONS.get(feature_id)
            if canonical is None:
                errors.append(
                    f"entry {feature_id!r}: command {command!r} is not wired into canonical orchestration"
                )
            elif command != canonical["command"]:
                errors.append(
                    f"entry {feature_id!r}: canonical orchestration requires command {canonical['command']!r}, got {command!r}"
                )
            argv = entrypoint.get("argv")
            if not isinstance(argv, list) or not argv:
                errors.append(f"entry {feature_id!r}: entrypoint.argv must be a non-empty list")
            else:
                for item in argv:
                    if not isinstance(item, str) or not item:
                        errors.append(f"entry {feature_id!r}: entrypoint.argv items must be non-empty strings")
                        continue
                    if "{" in item or "}" in item:
                        match = _PLACEHOLDER_RE.fullmatch(item)
                        if match is None:
                            errors.append(
                                f"entry {feature_id!r}: argv placeholder must be a full-match token, got {item!r}"
                            )
                            continue
                        placeholder = match.group(1)
                        if placeholder == "out":
                            continue
                        if not placeholder.startswith("fixture:"):
                            errors.append(f"entry {feature_id!r}: unknown argv placeholder {{{placeholder}}}")
                            continue
                        fixture_rel = placeholder[len("fixture:"):]
                        if fixture_rel not in entry.get("fixtures", []):
                            errors.append(
                                f"entry {feature_id!r}: argv placeholder mentions undeclared fixture {fixture_rel!r}"
                            )
        else:
            errors.append(f"entry {feature_id!r}: entrypoint must be an object with command and argv")

        fixtures = entry.get("fixtures")
        if not isinstance(fixtures, list) or not fixtures:
            errors.append(f"entry {feature_id!r}: fixtures must be a non-empty list")
        else:
            seen_fixtures: set[str] = set()
            for fixture in fixtures:
                if not isinstance(fixture, str) or not fixture:
                    errors.append(f"entry {feature_id!r}: fixture paths must be non-empty strings")
                    continue
                if "\\" in fixture:
                    errors.append(f"entry {feature_id!r}: fixture paths must use forward slashes, got {fixture!r}")
                if _is_absolute_fixture(fixture):
                    errors.append(f"entry {feature_id!r}: fixture path must be repo-relative, got {fixture!r}")
                if _has_dotdot_component(fixture):
                    errors.append(f"entry {feature_id!r}: fixture path must not contain a .. component, got {fixture!r}")
                if fixture in seen_fixtures:
                    errors.append(f"entry {feature_id!r}: duplicate fixture path {fixture!r}")
                seen_fixtures.add(fixture)
                if repo_root is not None:
                    candidate = _fixture_path(repo_root, fixture)
                    if not candidate.is_relative_to(repo_root.resolve()):
                        errors.append(f"entry {feature_id!r}: fixture escapes repo root: {fixture!r}")
                    elif not candidate.exists():
                        errors.append(f"entry {feature_id!r}: fixture does not exist: {fixture!r}")

        artifacts = entry.get("expected_artifacts")
        if not isinstance(artifacts, list) or not artifacts:
            errors.append(f"entry {feature_id!r}: expected_artifacts must be a non-empty list")
        else:
            seen_artifacts: set[str] = set()
            for artifact in artifacts:
                path_error = _output_path_error(artifact, allow_glob=True)
                if path_error is not None:
                    errors.append(f"entry {feature_id!r}: expected artifact {artifact!r} {path_error}")
                    continue
                if artifact in seen_artifacts:
                    errors.append(f"entry {feature_id!r}: duplicate expected artifact {artifact!r}")
                seen_artifacts.add(artifact)

        manifest_keys = entry.get("expected_manifest_keys")
        if not isinstance(manifest_keys, list) or not manifest_keys:
            errors.append(f"entry {feature_id!r}: expected_manifest_keys must be a non-empty list")
        else:
            for item in manifest_keys:
                if isinstance(item, dict):
                    artifact = item.get("artifact")
                    keys = item.get("keys")
                    unexpected = sorted(set(item) - {"artifact", "keys"})
                    if unexpected:
                        errors.append(f"entry {feature_id!r}: unexpected manifest-key entry field(s): {', '.join(unexpected)}")
                    if not isinstance(artifact, str) or not artifact:
                        errors.append(f"entry {feature_id!r}: manifest-key entry requires a non-empty artifact name")
                    elif not isinstance(artifacts, list) or not artifacts:
                        errors.append(f"entry {feature_id!r}: manifest-key artifact {artifact!r} cannot be checked without expected_artifacts")
                    elif not _is_declared_json_artifact(artifact, artifacts):
                        errors.append(
                            f"entry {feature_id!r}: manifest-key artifact {artifact!r} must be an exact declared JSON artifact"
                        )
                    elif _output_path_error(artifact, allow_glob=False) is not None:
                        errors.append(f"entry {feature_id!r}: manifest-key artifact {artifact!r} is not a safe out-dir-relative path")
                    if not isinstance(keys, list) or not keys or any(not isinstance(key, str) or not key for key in keys):
                        errors.append(f"entry {feature_id!r}: manifest-key entry requires a non-empty list of key paths")
                    continue
                errors.append(f"entry {feature_id!r}: manifest key entries must be {{artifact, keys}} objects bound to a declared JSON artifact")

        canonical = CANONICAL_INTEGRATIONS.get(feature_id)
        if canonical is not None:
            declared_manifests = {
                item.get("artifact") for item in manifest_keys
                if isinstance(item, dict) and isinstance(item.get("artifact"), str)
            } if isinstance(manifest_keys, list) else set()
            if canonical["manifest"] not in declared_manifests:
                errors.append(
                    f"entry {feature_id!r}: canonical manifest {canonical['manifest']!r} is not exposed by expected_manifest_keys"
                )

        rank_contribution = entry.get("rank_contribution")
        if rank_contribution not in ("none",) and not (
            isinstance(rank_contribution, str)
            and rank_contribution.startswith("policy:")
            and rank_contribution[len("policy:"):]
        ):
            errors.append(
                f"entry {feature_id!r}: rank_contribution must be 'none' or 'policy:<name>', got {rank_contribution!r}"
            )

        rank_assertion = entry.get("rank_assertion")
        if not isinstance(rank_assertion, dict):
            errors.append(f"entry {feature_id!r}: rank_assertion must be an object")
        else:
            assertion_type = rank_assertion.get("type")
            if assertion_type not in RANK_ASSERTION_TYPES:
                errors.append(
                    f"entry {feature_id!r}: rank_assertion.type must be one of {sorted(RANK_ASSERTION_TYPES)}, got {assertion_type!r}"
                )
            reason = rank_assertion.get("reason")
            if not isinstance(reason, str) or not reason.strip():
                errors.append(f"entry {feature_id!r}: rank_assertion.reason must be a non-empty string")

        for field in ("clean_install_smoke", "cache_fingerprint_assertion"):
            if field in entry and not isinstance(entry[field], bool):
                errors.append(f"entry {feature_id!r}: {field} must be a boolean")

        for field in ("owning_issue", "doc_reference"):
            value = entry.get(field)
            if not isinstance(value, str) or not value.strip():
                errors.append(f"entry {feature_id!r}: {field} must be a non-empty string")

        exit_codes = entry.get("expected_exit_codes", [0])
        if not isinstance(exit_codes, list) or not exit_codes or any(not isinstance(code, int) for code in exit_codes):
            errors.append(f"entry {feature_id!r}: expected_exit_codes must be a non-empty list of integers")

        assert_values = entry.get("assert_values")
        if assert_values is not None:
            if not isinstance(assert_values, dict) or set(assert_values) != {"artifact", "values"}:
                errors.append(f"entry {feature_id!r}: assert_values must be an object with artifact and values")
            else:
                artifact = assert_values.get("artifact")
                values = assert_values.get("values")
                if not isinstance(artifact, str) or not _is_declared_json_artifact(artifact, artifacts if isinstance(artifacts, list) else []):
                    errors.append(f"entry {feature_id!r}: assert_values.artifact must be an exact declared JSON artifact")
                if not isinstance(values, dict) or not values:
                    errors.append(f"entry {feature_id!r}: assert_values.values must be a non-empty object mapping dot-paths to literals")
                elif isinstance(artifact, str) and isinstance(manifest_keys, list):
                    bound_artifacts = {item.get("artifact") for item in manifest_keys if isinstance(item, dict)}
                    if artifact not in bound_artifacts:
                        errors.append(f"entry {feature_id!r}: assert_values.artifact {artifact!r} is not bound by expected_manifest_keys")
                for dot_path in values if isinstance(values, dict) else ():
                    if not isinstance(dot_path, str) or not dot_path.strip():
                        errors.append(f"entry {feature_id!r}: assert_values keys must be non-empty dot-paths")

    sorted_ids = [entry.get("feature_id") for entry in entries if isinstance(entry, dict)]
    if sorted_ids != sorted(sorted_ids, key=lambda value: value if isinstance(value, str) else ""):
        errors.append("entries must be sorted by feature_id ascending")
    return errors


def load_registry(registry_path: Path) -> dict:
    """Load and decode the registry JSON into a dict."""
    try:
        data = json.loads(Path(registry_path).read_text(encoding="utf-8"))
    except OSError as error:
        raise ReachabilityError("<registry>", f"cannot read registry {registry_path}: {error}") from error
    except json.JSONDecodeError as error:
        raise ReachabilityError("<registry>", f"registry {registry_path} is not valid JSON: {error}") from error
    if not isinstance(data, dict):
        raise ReachabilityError("<registry>", f"registry {registry_path} must be a JSON object")
    return data


def render_argv(entry: dict, repo_root: Path, out_dir: Path) -> list[str]:
    """Resolve one integrated entry's argv to concrete, absolute arguments.

    ``{out}`` is replaced with the entry's out-dir and ``{fixture:<relpath>}``
    with the repo-relative fixture's absolute path.  No shell is involved:
    the argv is a structured list passed straight to subprocess.
    """
    rendered: list[str] = []
    for item in entry["entrypoint"]["argv"]:
        if item == "{out}":
            rendered.append(str(out_dir))
            continue
        match = _PLACEHOLDER_RE.fullmatch(item)
        if match and match.group(1).startswith("fixture:"):
            relpath = match.group(1)[len("fixture:"):]
            rendered.append(str(_fixture_path(repo_root, relpath)))
            continue
        rendered.append(item)
    return rendered


def _feature_env(env: dict[str, str] | None) -> dict[str, str]:
    """Caller-supplied environment with PYTHONPATH removed.

    A leftover interactive ``PYTHONPATH`` would silently redirect the feature
    interpreter's import (weakening source mode) or defeat the clean-install
    proof (weakening installed mode), so it is stripped from every gate
    subprocess.
    """
    merged = os.environ.copy()
    if env is not None:
        merged.update(env)
    merged.pop("PYTHONPATH", None)
    return merged


def _expected_artifact_matches(out_dir: Path, name: str) -> list[Path]:
    path_error = _output_path_error(name, allow_glob=True)
    if path_error is not None:
        return []
    if any(character in name for character in "*?["):
        return sorted(
            path for path in out_dir.rglob("*")
            if path.is_file() and fnmatch.fnmatch(path.relative_to(out_dir).as_posix(), name)
        )
    path = out_dir / name
    return [path] if path.is_file() else []


def _check_manifest_key(feature_id: str, out_dir: Path, key_path: str) -> bool:
    manifests = _emitted_manifests(out_dir)
    for path in manifests:
        found, _ = _dotpath_get(json.loads(path.read_text(encoding="utf-8")), key_path)
        if found:
            return True
    raise ReachabilityError(feature_id, f"manifest key not found in any emitted manifest: {key_path!r}")


def run_entry(
    entry: dict,
    repo_root: Path,
    out_dir: Path,
    python_exe: str,
    timeout: int = 300,
    env: dict[str, str] | None = None,
    cwd: Path | None = None,
) -> dict:
    """Run one integrated entry as a real subprocess and assert its contract.

    Raises :class:`ReachabilityError` (carrying the entry's feature_id) on a
    non-expected exit code, a missing or empty expected artifact, a missing
    manifest key, or an ``assert_values`` mismatch.  The return dict carries
    the observed exit code, output tails, artifact sizes and manifest-key
    presence for the summary line.
    """
    feature_id = str(entry["feature_id"])
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    command = str(entry["entrypoint"]["command"])
    argv = render_argv(entry, repo_root, out_dir)
    expected_exit_codes = [int(code) for code in entry.get("expected_exit_codes", [0])]

    try:
        module = "candidate_qtl.bench_cli" if entry["entrypoint"].get("surface", "core") == "benchmark" else "candidate_qtl"
        completed = subprocess.run(
            [python_exe, "-m", module, command, *argv],
            capture_output=True,
            text=True,
            timeout=timeout,
            env=_feature_env(env),
            cwd=str(cwd) if cwd is not None else None,
        )
    except subprocess.TimeoutExpired as error:
        raise ReachabilityError(feature_id, f"command timed out after {timeout}s") from error
    except OSError as error:
        raise ReachabilityError(feature_id, f"cannot run {python_exe}: {error}") from error

    stdout_tail = (completed.stdout or "")[-2000:]
    stderr_tail = (completed.stderr or "")[-2000:]
    exit_code = int(completed.returncode)
    if exit_code not in expected_exit_codes:
        raise ReachabilityError(
            feature_id,
            f"command exited {exit_code}; expected {expected_exit_codes}; stderr tail: {stderr_tail.strip() or '(empty)'}",
        )

    artifacts: dict[str, dict[str, object]] = {}
    for name in entry["expected_artifacts"]:
        matches = _expected_artifact_matches(out_dir, name)
        if not matches:
            raise ReachabilityError(feature_id, f"expected artifact not produced: {name!r}")
        size = matches[0].stat().st_size
        if size == 0:
            raise ReachabilityError(feature_id, f"expected artifact is empty: {name!r}")
        artifacts[name] = {"exists": True, "size": size}

    manifest_keys: dict[str, bool] = {}
    try:
        for item in entry["expected_manifest_keys"]:
            artifact = item["artifact"]
            if not _is_declared_json_artifact(artifact, entry["expected_artifacts"]):
                raise ReachabilityError(feature_id, f"manifest artifact is not an exact declared JSON artifact: {artifact!r}")
            matches = _expected_artifact_matches(out_dir, artifact)
            if len(matches) != 1:
                raise ReachabilityError(feature_id, f"manifest artifact not produced exactly once: {artifact!r}")
            data = json.loads(matches[0].read_text(encoding="utf-8"))
            for key_path in item["keys"]:
                found, _ = _dotpath_get(data, key_path)
                manifest_keys[f"{artifact}:{key_path}"] = found
                if not found:
                    raise ReachabilityError(
                        feature_id,
                        f"manifest key not found in {artifact!r}: {key_path!r}",
                    )
        assertion = entry.get("assert_values")
        if assertion:
            artifact = assertion["artifact"]
            matches = _expected_artifact_matches(out_dir, artifact)
            if len(matches) != 1:
                raise ReachabilityError(feature_id, f"assert_values artifact not produced exactly once: {artifact!r}")
            data = json.loads(matches[0].read_text(encoding="utf-8"))
            for dot_path, expected in assertion["values"].items():
                found, actual = _dotpath_get(data, dot_path)
                if not found or actual != expected:
                    raise ReachabilityError(
                        feature_id,
                        f"assert_values mismatch in {artifact!r} for {dot_path!r}: expected {expected!r}, got {actual!r}",
                    )
    except ReachabilityError:
        raise
    except (OSError, json.JSONDecodeError) as error:
        raise ReachabilityError(feature_id, f"cannot read emitted manifest: {error}") from error

    return {
        "feature_id": feature_id,
        "command": command,
        "exit_code": exit_code,
        "stdout_tail": stdout_tail,
        "stderr_tail": stderr_tail,
        "artifacts": artifacts,
        "manifest_keys": manifest_keys,
    }


def _matches_rank_impacting(relative_path: str) -> bool:
    """Whether an out-dir file matches the rank-impacting artifact patterns."""
    relative = relative_path.replace("\\", "/")
    name = relative.rsplit("/", 1)[-1]
    return any(
        fnmatch.fnmatch(relative, pattern) or fnmatch.fnmatch(name, pattern)
        for pattern in RANK_IMPACTING_PATTERNS
    )


def _fixture_sha256(repo_root: Path, fixture: str) -> str:
    """Content hash of one declared fixture; a directory hashes its whole tree.

    A directory fixture (e.g. an ``extract`` folder handed to a subcommand via
    ``{fixture:...}``) is hashed deterministically: every file below it,
    sorted by relative path, contributes its relative path and bytes.
    """
    path = _fixture_path(repo_root, fixture)
    if path.is_dir():
        digest = hashlib.sha256()
        for child in sorted(path.rglob("*")):
            if child.is_file():
                digest.update(child.relative_to(path).as_posix().encode("utf-8"))
                digest.update(b"\x00")
                digest.update(child.read_bytes())
        return digest.hexdigest()
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _run_canonical_baseline(repo_root: Path, out_dir: Path, python_exe: str, timeout: int, cwd: Path) -> dict[str, str]:
    """Execute the smallest full case-run fixture and hash rank-sensitive outputs."""
    case, resources = (_fixture_path(repo_root, item) for item in CANONICAL_BASELINE_FIXTURES)
    try:
        completed = subprocess.run(
            [
                python_exe, "-m", "candidate_qtl", "case-run",
                "--case", str(case), "--resources", str(resources), "--out-dir", str(out_dir),
            ],
            capture_output=True,
            text=True,
            timeout=timeout,
            env=_feature_env(None),
            cwd=str(cwd),
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ReachabilityError("<canonical-baseline>", f"cannot run case-run baseline: {error}") from error
    if completed.returncode != 0:
        raise ReachabilityError(
            "<canonical-baseline>",
            f"case-run exited {completed.returncode}; stderr tail: {(completed.stderr or '')[-2000:].strip() or '(empty)'}",
        )
    hashes: dict[str, str] = {}
    for artifact in CANONICAL_BASELINE_ARTIFACTS:
        path = out_dir / artifact
        if not path.is_file() or path.stat().st_size == 0:
            raise ReachabilityError("<canonical-baseline>", f"baseline artifact not produced: {artifact!r}")
        hashes[artifact] = hashlib.sha256(path.read_bytes()).hexdigest()
    return hashes


def run_gate(
    registry_path: Path,
    repo_root: Path,
    mode: str,
    python_exe: str,
    out_root: Path | None = None,
    timeout: int = 300,
) -> dict:
    """Run the full reachability gate: validate, execute, assert, aggregate.

    ``mode="source"`` runs each entry with cwd = repo_root (the import is the
    source checkout); ``mode="installed"`` runs each entry with cwd = a fresh
    temp directory outside repo_root (any import must come from the installed
    package).  Every integrated entry runs even when a previous one failed;
    the failure of any entry, of fixture immutability, or of the no-rank-artifacts
    assertion makes ``valid`` False.  A registry with zero integrated entries
    is a hard failure.
    """
    registry_path = Path(registry_path)
    repo_root = Path(repo_root).resolve()
    data = load_registry(registry_path)
    errors = validate_registry(data, repo_root=repo_root)
    if errors:
        raise ReachabilityError(
            "<registry>",
            f"registry validation failed ({len(errors)} error(s)): {'; '.join(errors)}",
        )
    if mode not in ("source", "installed"):
        raise ReachabilityError("<registry>", f"unknown mode {mode!r}")

    integrated = sorted(
        (entry for entry in data["entries"] if entry["integration_status"] == STATUS_INTEGRATED),
        key=lambda entry: entry["feature_id"],
    )
    core_only_count = sum(1 for entry in data["entries"] if entry["integration_status"] == STATUS_CORE_ONLY)
    if not integrated:
        raise ReachabilityError("<registry>", "no integrated entries in the registry")

    out_root = Path(out_root) if out_root is not None else Path(tempfile.mkdtemp(prefix="reachability-gate-"))
    out_root.mkdir(parents=True, exist_ok=True)
    installed_cwd = Path(tempfile.mkdtemp(dir=out_root, prefix="installed-cwd-")) if mode == "installed" else None
    run_cwd = repo_root if mode == "source" else installed_cwd

    baseline_before = _run_canonical_baseline(
        repo_root, out_root / "canonical-baseline-before", python_exe, timeout, run_cwd
    )

    entries_result: list[dict[str, object]] = []
    valid = True
    for entry in integrated:
        feature_id = str(entry["feature_id"])
        out_dir = out_root / feature_id
        out_dir.mkdir(parents=True, exist_ok=True)
        entry_result: dict[str, object] = {
            "feature_id": feature_id,
            "command": str(entry["entrypoint"]["command"]),
            "out_dir": str(out_dir),
            "error": None,
            "fixture_hashes": {},
        }
        try:
            before = {fixture: _fixture_sha256(repo_root, fixture) for fixture in entry["fixtures"]}
            entry_result["fixture_hashes"] = before
            run_result = run_entry(
                entry,
                repo_root,
                out_dir,
                python_exe,
                timeout=timeout,
                cwd=run_cwd,
            )
            entry_result.update(run_result)
            if entry["rank_contribution"] == "none":
                offending = [
                    path.relative_to(out_dir).as_posix()
                    for path in out_dir.rglob("*")
                    if path.is_file() and _matches_rank_impacting(path.relative_to(out_dir).as_posix())
                ]
                if offending:
                    raise ReachabilityError(feature_id, f"rank-impacting artifact(s) written: {', '.join(sorted(offending))}")
            after = {fixture: _fixture_sha256(repo_root, fixture) for fixture in entry["fixtures"]}
            changed = sorted(fixture for fixture in before if before[fixture] != after[fixture])
            if changed:
                raise ReachabilityError(
                    feature_id, f"fixture(s) changed during the run: {', '.join(changed)}"
                )
        except ReachabilityError as error:
            valid = False
            entry_result["error"] = str(error)
        entries_result.append(entry_result)

    try:
        baseline_after = _run_canonical_baseline(
            repo_root, out_root / "canonical-baseline-after", python_exe, timeout, run_cwd
        )
        if baseline_before != baseline_after:
            valid = False
            changed = sorted(name for name in baseline_before if baseline_before[name] != baseline_after.get(name))
            baseline_error = f"rank-sensitive canonical baseline changed: {', '.join(changed)}"
        else:
            baseline_error = None
    except ReachabilityError as error:
        valid = False
        baseline_after = {}
        baseline_error = str(error)

    return {
        "valid": valid,
        "entries": entries_result,
        "integrated_count": len(integrated),
        "core_only_count": core_only_count,
        "out_root": str(out_root),
        "mode": mode,
        "canonical_baseline": {
            "artifacts": list(CANONICAL_BASELINE_ARTIFACTS),
            "before": baseline_before,
            "after": baseline_after,
            "error": baseline_error,
        },
    }


def assert_no_source_fallback(python_exe: str, repo_root: Path, cwd: Path | None = None) -> str:
    """Prove the given interpreter imports candidate_qtl from site-packages.

    Runs ``import candidate_qtl`` in a temp directory outside ``repo_root`` and
    fails when the resolved module path lives inside the source checkout or
    outside site-packages.  This is the clean-install proof: a bare wheel
    install must resolve into site-packages, never into the repo's ``src/``.
    """
    repo_root = Path(repo_root).resolve()
    created_cwd = cwd is None
    run_dir = Path(cwd) if cwd is not None else Path(tempfile.mkdtemp(prefix="reachability-cwd-"))
    try:
        completed = subprocess.run(
            [python_exe, "-c", "import candidate_qtl,sys;print(candidate_qtl.__file__)"],
            capture_output=True,
            text=True,
            timeout=60,
            env=_feature_env(None),
            cwd=str(run_dir),
        )
    except (subprocess.TimeoutExpired, OSError) as error:
        raise ReachabilityError("<clean-install>", f"cannot probe {python_exe}: {error}") from error
    finally:
        if created_cwd:
            shutil.rmtree(run_dir, ignore_errors=True)
    if completed.returncode != 0:
        raise ReachabilityError(
            "<clean-install>",
            f"candidate_qtl is not importable in {python_exe} (exit {completed.returncode}): "
            f"{(completed.stderr or '').strip()[-400:]}",
        )
    printed = [line.strip() for line in (completed.stdout or "").splitlines() if line.strip()]
    if not printed:
        raise ReachabilityError("<clean-install>", f"{python_exe} printed no import path")
    import_path = printed[-1]
    resolved = Path(import_path).resolve()
    if resolved.is_relative_to(repo_root):
        raise ReachabilityError(
            "<clean-install>",
            f"import falls back to the source checkout: {import_path} (repo root {repo_root})",
        )
    if "site-packages" not in import_path.lower():
        raise ReachabilityError(
            "<clean-install>",
            f"import path is not from site-packages (editable install?): {import_path}",
        )
    return import_path
