"""Resolve an existing reusable family index for guided orchestration.

Resolution is deliberately an orchestration concern: it only locates an existing
artifact and never infers families, changes evidence, or turns missing data into
negative evidence. Downstream family-index consumers retain their own schema and
scope validation.
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
from typing import Any

import yaml

from .core import runtime_path


FAMILY_INDEX_ENV = "CANDIDATE_QTL_FAMILY_INDEX"
REQUIRED_COLUMNS = {"gene_id", "assignment_status", "best_subject", "family_stems"}


def _default_candidates() -> list[Path]:
    """Return portable per-user cache locations only.

    Repository defaults must not depend on a maintainer-specific workstation
    layout. A built-once index stored elsewhere should be exposed through
    ``CANDIDATE_QTL_FAMILY_INDEX``, ``resources.yaml``, or a symlink/copy in one
    of these conventional cache locations. This keeps no-provider behavior and
    tests deterministic across machines (#643).
    """
    cache_root = runtime_path(os.environ.get("XDG_CACHE_HOME") or (Path.home() / ".cache"))
    return [
        cache_root / "candidate-qtl" / "family_index.tsv",
        Path.home() / ".candidate-qtl" / "family_index.tsv",
    ]


def _existing_file(path: Path) -> bool:
    try:
        return path.is_file()
    except OSError:
        return False


def _reusable_index_validation(path: Path) -> dict[str, object]:
    """Validate that an auto-discovered table is a reusable genome-wide index."""
    validation: dict[str, object] = {
        "status": "invalid",
        "index_scope": None,
        "manifest": None,
        "reference_label": None,
        "reference_species": None,
    }
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            columns = set(reader.fieldnames or ())
            missing = sorted(REQUIRED_COLUMNS - columns)
            if missing:
                validation.update({"reason": "missing_columns", "missing_columns": missing})
                return validation
            scopes = {str(row.get("index_scope") or "").strip() for row in reader}
            scopes.discard("")
    except (OSError, UnicodeDecodeError, csv.Error) as error:
        validation.update({"reason": "unreadable_index", "error": f"{type(error).__name__}: {error}"})
        return validation

    manifest_path = path.parent / "family_index_manifest.json"
    manifest_scope = ""
    if _existing_file(manifest_path):
        validation["manifest"] = str(manifest_path)
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
            validation.update({"reason": "unreadable_manifest", "error": f"{type(error).__name__}: {error}"})
            return validation
        if not isinstance(manifest, dict):
            validation["reason"] = "invalid_manifest"
            return validation
        manifest_scope = str(manifest.get("index_scope") or "").strip()
        validation["reference_label"] = str(manifest.get("reference_label") or "").strip() or None
        validation["reference_species"] = str(manifest.get("reference_species") or "").strip() or None

    observed_scopes = set(scopes)
    if manifest_scope:
        observed_scopes.add(manifest_scope)
    validation["index_scope"] = ";".join(sorted(observed_scopes)) or None
    if observed_scopes != {"reference_genome"}:
        validation["reason"] = "not_reusable_reference_genome_index"
        return validation
    validation.update({"status": "valid", "reason": None})
    return validation


def _validated_candidate(path: Path, attempt: dict[str, object]) -> bool:
    if not _existing_file(path):
        attempt["status"] = "file_missing"
        return False
    validation = _reusable_index_validation(path)
    attempt["validation"] = validation
    attempt["status"] = "resolved" if validation["status"] == "valid" else "invalid_reusable_index"
    return validation["status"] == "valid"


def _registry_candidate(resources_arg: str | None) -> tuple[Path | None, dict[str, object]]:
    attempt: dict[str, object] = {"source": "resources", "status": "not_configured", "path": None}
    if not resources_arg:
        return None, attempt
    resources_path = runtime_path(resources_arg)
    attempt["registry"] = str(resources_path)
    if not _existing_file(resources_path):
        attempt["status"] = "registry_unavailable"
        return None, attempt
    try:
        registry = yaml.safe_load(resources_path.read_text(encoding="utf-8")) or {}
    except (OSError, UnicodeDecodeError, yaml.YAMLError) as error:
        attempt.update({"status": "registry_unreadable", "error": f"{type(error).__name__}: {error}"})
        return None, attempt
    if not isinstance(registry, dict):
        attempt["status"] = "registry_invalid"
        return None, attempt
    resource_map = registry.get("resources", registry)
    if not isinstance(resource_map, dict):
        attempt["status"] = "registry_invalid"
        return None, attempt
    entry: Any = resource_map.get("family_index")
    if not entry:
        return None, attempt
    raw = entry.get("path") if isinstance(entry, dict) else entry
    if not raw:
        attempt["status"] = "path_missing"
        return None, attempt
    candidate = runtime_path(str(raw))
    if not candidate.is_absolute():
        candidate = resources_path.parent / candidate
    attempt["path"] = str(candidate)
    if not _validated_candidate(candidate, attempt):
        return None, attempt
    return candidate, attempt


def resolve_family_index(
    explicit_arg: str | None,
    resources_arg: str | None = None,
) -> tuple[str | None, dict[str, object]]:
    """Resolve a reusable family index without changing missing-data semantics.

    Priority is explicit CLI > resources registry > environment > conventional
    per-user cache paths. Only auto-discovered paths that exist as files are
    returned. An explicit CLI value is preserved as-is so downstream validation
    remains byte-for-byte compatible with the previous command behavior.
    """
    if explicit_arg:
        return explicit_arg, {
            "attempted": False,
            "status": "explicit",
            "source": "explicit_cli",
            "path": explicit_arg,
            "environment_variable": FAMILY_INDEX_ENV,
            "attempts": [],
            "rank_contribution": "none",
        }

    attempts: list[dict[str, object]] = []
    registry_path, registry_attempt = _registry_candidate(resources_arg)
    attempts.append(registry_attempt)
    if registry_path is not None:
        return str(registry_path), {
            "attempted": True,
            "status": "resolved",
            "source": "resources",
            "path": str(registry_path),
            "environment_variable": FAMILY_INDEX_ENV,
            "attempts": attempts,
            "rank_contribution": "none",
        }

    env_value = os.environ.get(FAMILY_INDEX_ENV, "").strip()
    env_attempt: dict[str, object] = {"source": "environment", "status": "not_configured", "path": None}
    if env_value:
        env_path = runtime_path(env_value)
        env_attempt["path"] = str(env_path)
        _validated_candidate(env_path, env_attempt)
        attempts.append(env_attempt)
        if env_attempt["status"] == "resolved":
            return str(env_path), {
                "attempted": True,
                "status": "resolved",
                "source": "environment",
                "path": str(env_path),
                "environment_variable": FAMILY_INDEX_ENV,
                "attempts": attempts,
                "rank_contribution": "none",
            }
    else:
        attempts.append(env_attempt)

    for candidate in _default_candidates():
        attempt = {
            "source": "default_cache",
            "status": "not_checked",
            "path": str(candidate),
        }
        _validated_candidate(candidate, attempt)
        attempts.append(attempt)
        if attempt["status"] == "resolved":
            return str(candidate), {
                "attempted": True,
                "status": "resolved",
                "source": "default_cache",
                "path": str(candidate),
                "environment_variable": FAMILY_INDEX_ENV,
                "attempts": attempts,
                "rank_contribution": "none",
            }

    return None, {
        "attempted": True,
        "status": "not_found",
        "source": None,
        "path": None,
        "environment_variable": FAMILY_INDEX_ENV,
        "attempts": attempts,
        "rank_contribution": "none",
    }
