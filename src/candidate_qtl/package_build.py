"""Executor package builder (issue #365).

Builds a physically separate Executor-visible package: code, executor
manifests and declared reference metadata only.  The allowlist is the
primary contract; the negative scanner is the fail-closed backstop for
every file that matches an include pattern.  A rejected file aborts the
build; an undeclared file is excluded and recorded in the manifest as
review_required instead of being silently packaged.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import platform
import re
import shutil
import stat
import subprocess
import sys
import tomllib
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from fnmatch import fnmatch
from pathlib import Path

from . import __version__
from .prediction_contract import ANSWER_SHAPED_COLUMNS
from .command_log import CANONICAL_COLUMNS
from .core import file_fingerprint, runtime_path, sha256

PACKAGE_MANIFEST_VERSION = "1.0"
PACKAGE_VERSION = "1.0.0"
CAPSULE_VERSION = "1.0"
POST_FREEZE_PHASE = "post_freeze_scoring"

REPARSE_POINT = getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0x400)


@dataclass(frozen=True)
class AllowlistRule:
    pattern: str
    kind: str
    reason: str


@dataclass(frozen=True)
class FileRecord:
    rel: str
    size: int
    sha256: str
    classification: str
    rules: tuple[str, ...]
    reason: str
    kind: str
    provenance: str
    allowlist_rule: str | None
    executor_input_schema: str | None = None


@dataclass
class ScanResult:
    files: dict[str, FileRecord]
    rejected: list[str]
    review_required: list[str]
    unmatched_includes: list[str]


@dataclass
class BuildResult:
    out_dir: Path
    manifest_path: Path
    capsule_path: Path
    file_count: int
    total_bytes: int
    rejected: list[str]
    review_required: list[str]
    truth_inputs_present: bool


DEFAULT_ALLOWLIST = (
    AllowlistRule("src/**", "source_code", "candidate-qtl source package"),
    AllowlistRule("pyproject.toml", "config", "package metadata"),
    AllowlistRule("constraints.txt", "config", "pinned dependency lock"),
    AllowlistRule("benchmarks/*/executor/**", "executor_data", "answerless executor case inputs"),
    AllowlistRule("benchmarks/practical_blind_v1/input/**", "executor_data", "answerless executor inputs"),
    AllowlistRule("templates/benchmark_case_manifest.tsv", "executor_data", "executor case manifest template"),
    AllowlistRule("templates/benchmark_source_audit.tsv", "executor_data", "executor source audit template"),
)

DENY_PATH_PATTERNS = (
    "benchmarks/*/truth/**",
    "benchmarks/*/source/**",
    "benchmarks/*/score*/**",
    "benchmarks/*/freeze/**",
    "benchmarks/*/run_manifest.json",
    "benchmarks/*negative_control*/**",
    "benchmarks/dual_route_v1/truth.tsv",
    "benchmarks/dual_route_v1/cases.tsv",
    "benchmarks/phenotype_variant_v1/manifest.tsv",
    "benchmarks/comparative_locus_v1/expected_truth_set.tsv",
    "benchmarks/comparative_locus_v1/benchmark_manifest.tsv",
    "benchmarks/comparative_locus_v1/focal_recovery.tsv",
    "templates/benchmark_truth_package.tsv",
    "templates/v2_*_benchmark.tsv",
    "templates/pea_truth_catalog_expansion_*.tsv",
    "templates/cross_crop_benchmark.tsv",
    "docs/tutorial/known_qtg_benchmark.tsv",
    "tests/fixtures/blind_benchmark/manifest.tsv",
)

TRUTH_COLUMN_NAMES = frozenset(
    ANSWER_SHAPED_COLUMNS
    | {
        "answer_gene_name",
        "answer_locus",
        "answer_source_id",
        "answer_mapping_source",
        "truth_accession",
        "truth_primary_source",
        "truth_mapping_source",
        "truth_provenance_tier",
        "truth_unit_type",
        "truth_confidence_notes",
        "truth_tier",
        "known_answer_aliases",
        "known_answer_derived_sources",
        "profile_disposition",
        "disposition_reason",
        "gene_class",
        "causal_variant",
        "causal_variant_class",
        "causal_publication_ids",
        "expected_recovery",
        "expected_action_class",
        "decoy_family",
        "control_type",
        "gene_label",
    }
)

TRUTH_JSON_KEYS = frozenset(
    {
        "truth",
        "truth_gene_id",
        "truth_family_id",
        "truth_unit_type",
        "truth_provenance_tier",
        "truth_accession",
        "truth_primary_source",
        "truth_mapping_source",
        "answer_gene_id",
        "answer_pathway",
        "answer_source_id",
        "answer_mapping_source",
        "known_answer_aliases",
        "known_answer_derived_sources",
        "causal_publication_ids",
        "profile_disposition",
        "disposition_reason",
        "score_key",
        "historical_answer",
        "ground_truth",
        "expected_gene_id",
        "expected_recovery",
        "decoy_family",
    }
)

_SHA256 = re.compile(r"^[0-9a-f]{64}$")

# Executor-visible data is deliberately a small, typed surface.  This is a
# physical boundary rather than a claim that a generic text search can infer a
# secret biological answer from arbitrary prose.  New input types require an
# explicit schema addition and review; an arbitrary TSV cannot be smuggled in
# merely by passing --include.
EXECUTOR_TSV_SCHEMAS: dict[str, tuple[frozenset[str], ...]] = {
    "case_manifest": (
        frozenset({"case_id", "benchmark_profile", "phenotype_text", "focal_species", "inheritance", "tissue",
                   "developmental_stage", "interval_id", "interval_source", "reference_build", "candidate_input",
                   "reference_corpus", "knowledge_cutoff", "case_selection_tier_public", "allowed_external_search",
                   "notes_for_executor"}),
        frozenset({"case_id", "benchmark_profile", "phenotype_text", "focal_species", "interval_id",
                   "interval_source", "reference_build", "candidate_input", "case_selection_tier_public",
                   "allowed_external_search"}),
        frozenset({"case_id", "benchmark_profile", "benchmark_stratum", "phenotype_text", "focal_species",
                   "inheritance", "tissue", "developmental_stage", "interval_id", "interval_source",
                   "reference_build", "interval", "candidate_input", "reference_corpus", "knowledge_cutoff",
                   "case_selection_tier_public", "allowed_external_search", "notes_for_executor"}),
        frozenset({"case_id", "phenotype_text", "interval", "coordinate_basis", "focal_species",
                   "inheritance_if_source_explicit", "candidate_input", "knowledge_policy", "input_provenance_id"}),
    ),
    "interval_candidates": (
        frozenset({"gene_id", "annotation"}),
        frozenset({"gene_id", "gene_name", "gene_biotype", "is_placeholder_gene", "description", "dbxref",
                   "ontology_terms", "interval", "coordinate_basis", "gff_seqid", "fasta_seqid", "chrom",
                   "start", "end", "strand", "gene_length", "transcript_count", "representative_transcript",
                   "representative_reason", "cds_length", "protein_length", "cds_phase_status",
                   "translation_status", "annotation_quality", "overlaps_qtl_boundary"}),
    ),
    "route": (frozenset({"phenotype_token", "pathway_family", "annotation_keywords"}),),
    "source_log": (frozenset({"case_id", "route", "query_or_lookup", "source_id", "source_type",
                               "selected_entity", "selection_reason", "retrieved_at", "actually_used_in_prediction"}),),
    "locus_hypotheses": (frozenset({"hypothesis_id", "case_id", "hypothesis_type", "interval_or_element_id",
                                     "member_gene_ids", "target_gene_ids", "proposed_mechanism",
                                     "predicted_molecular_effect", "predicted_phenotype_direction",
                                     "supporting_evidence_ids", "counterevidence_ids", "context_evidence_ids",
                                     "required_experiment_ids", "status", "requested_status", "status_guardrail",
                                     "status_reason", "created_by", "proposal_origin", "review_status", "source_id",
                                     "parent_hypothesis_id", "supersedes_hypothesis_id", "gene_disposition",
                                     "rejection_scope"}),),
    "reference_assets": (
        frozenset({"reference_build", "asset_kind", "path_wsl", "notes"}),
        frozenset({"reference_build", "asset_kind", "path_wsl", "asset_version", "source_fasta_sha256",
                   "sequence_count", "usage_scope", "pin_policy", "notes"}),
    ),
    "command_log": (
        frozenset({"case_id", "stage", "command", "exit_code", "start_utc", "end_utc", "notes"}),
        frozenset(CANONICAL_COLUMNS),
    ),
    "executor_manifest": (frozenset({"case_id", "phenotype_text", "interval", "coordinate_basis", "focal_species",
                                      "inheritance_if_source_explicit", "candidate_input", "knowledge_policy",
                                      "input_provenance_id"}),),
    "benchmark_source_audit_template": (frozenset({"case_id", "route", "query_or_lookup", "source_id",
                                                     "source_type", "selected_entity", "selection_reason",
                                                     "retrieved_at", "supervisor_class", "supervisor_action",
                                                     "decision_reason"}),),
}

_TRUTH_VALUE_MARKER = re.compile(
    r"(?:known[ _-]?(?:answer|gene)|ground[ _-]?truth|causal[ _-]?(?:gene|variant)|"
    r"(?:answer|truth)[ _-]?alias)\s*[:=]",
    re.IGNORECASE,
)

_SKIP_DIRNAMES = {".git", "__pycache__", ".venv"}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _write_manifest(path: Path, payload: dict[str, object]) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _format_record(record: FileRecord) -> str:
    return f"{record.rel} | {','.join(record.rules) or '-'} | {record.classification} | {record.reason}"


def _is_reparse_point(path: Path) -> bool:
    if os.name != "nt":
        return False
    try:
        attributes = os.lstat(path).st_file_attributes
    except OSError:
        return False
    return bool(attributes & REPARSE_POINT)


def _resolved_under(path: Path, root: Path) -> bool:
    try:
        return path.resolve().is_relative_to(root)
    except OSError:
        return False


def detect_repo_root(start: Path | None = None) -> Path:
    current = (start or Path(__file__).resolve()).resolve()
    if current.is_file():
        current = current.parent
    for candidate in (current, *current.parents):
        pyproject = candidate / "pyproject.toml"
        if not pyproject.is_file():
            continue
        try:
            with pyproject.open("rb") as handle:
                data = tomllib.load(handle)
        except (OSError, tomllib.TOMLDecodeError):
            continue
        if data.get("project", {}).get("name") == "candidate-qtl":
            return candidate
    raise ValueError("could not locate candidate-qtl repository root (pyproject.toml with project.name = candidate-qtl)")


def _tsv_header(path: Path) -> list[str]:
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            first = handle.readline().lstrip("\ufeff")
    except OSError:
        return []
    if not first.strip():
        return []
    try:
        row = next(csv.reader([first], delimiter="\t"))
    except csv.Error:
        return []
    return [column.strip() for column in row if column and column.strip()]


_TRUTH_PREFIXES = (
    "answer_",
    "truth_",
    "causal_",
    "known_gene",
    "known_answer",
    "expected_gene",
    "ground_truth",
    "solution_",
    "decoy_",
    "leak_",
)

_TRUTH_INFIXES = (
    "_truth",
    "_answer",
    "_causal",
    "_known_gene",
    "_known_answer",
    "_expected_gene",
    "_ground_truth",
    "_solution",
    "_decoy",
    "_leak",
)


def _truth_shaped_columns(columns: list[str]) -> list[str]:
    hits: list[str] = []
    for column in columns:
        lowered = column.strip().lower()
        if column in TRUTH_COLUMN_NAMES or lowered in TRUTH_COLUMN_NAMES:
            hits.append(column)
            continue
        if lowered.startswith(_TRUTH_PREFIXES) or any(infix in lowered for infix in _TRUTH_INFIXES):
            hits.append(column)
    return sorted(hits)


def _executor_input_schema(rel: str, path: Path) -> tuple[str | None, str | None]:
    """Return a typed Executor input classification or its rejection reason."""
    if path.suffix == ".sha256":
        try:
            content = path.read_text(encoding="ascii").strip()
        except OSError:
            return None, "input schema: cannot read SHA-256 sidecar"
        if _SHA256.fullmatch(content.lower()):
            return "sha256_sidecar", None
        return None, "input schema: SHA-256 sidecar must contain exactly one digest"
    if path.suffix != ".tsv":
        return None, "input schema: Executor data must be a UTF-8 TSV or .sha256 sidecar"

    columns = _tsv_header(path)
    column_set = frozenset(columns)
    if not columns or len(columns) != len(column_set):
        return None, "input schema: TSV header must be non-empty and have unique columns"
    for name, accepted in EXECUTOR_TSV_SCHEMAS.items():
        if column_set in accepted:
            return name, None
    return None, f"input schema: undeclared Executor TSV columns {sorted(column_set)}"


def _truth_value_marker(path: Path) -> str | None:
    """Find explicit answer-label/value markers in Executor input cells.

    This intentionally does not reject ordinary scientific prose containing the
    word "truth".  Typed schemas provide the primary boundary; this catches the
    common attempt to smuggle a known-answer alias as a labelled cell value.
    """
    if path.suffix != ".tsv":
        return None
    try:
        with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            next(reader, None)
            for line_no, row in enumerate(reader, start=2):
                for value in row:
                    if _TRUTH_VALUE_MARKER.search(value):
                        return f"value rules: labelled truth/answer value at TSV line {line_no}"
    except (OSError, csv.Error):
        return "value rules: could not scan Executor TSV values"
    return None


def _json_document(path: Path):
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            return json.load(handle)
    except (OSError, json.JSONDecodeError):
        return None


def _iter_json_keys(value: object):
    if isinstance(value, dict):
        for key, item in value.items():
            yield str(key)
            yield from _iter_json_keys(item)
    elif isinstance(value, list):
        for item in value:
            yield from _iter_json_keys(item)


def _json_has_post_freeze_phase(value: object) -> bool:
    if isinstance(value, dict):
        if value.get("phase") == POST_FREEZE_PHASE:
            return True
        return any(_json_has_post_freeze_phase(item) for item in value.values())
    if isinstance(value, list):
        return any(_json_has_post_freeze_phase(item) for item in value)
    return False


def _match_allowlist(rel: str, patterns: list[str]) -> AllowlistRule | None:
    for rule in DEFAULT_ALLOWLIST:
        if rule.pattern in patterns and fnmatch(rel, rule.pattern):
            return rule
    default_patterns = {rule.pattern for rule in DEFAULT_ALLOWLIST}
    for pattern in patterns:
        if pattern not in default_patterns and fnmatch(rel, pattern):
            return AllowlistRule(pattern, "declared", "user-declared include")
    return None


def _is_executor_input_rel(rel: str) -> bool:
    return "/executor/" in f"/{rel}" or rel.startswith("benchmarks/practical_blind_v1/input/")


def scan_package_source(source_root: Path, include_patterns: list[str]) -> ScanResult:
    root = source_root.resolve()
    if not root.is_dir():
        raise ValueError(f"source root is not a directory: {root}")
    for pattern in include_patterns:
        if pattern.startswith(("/", "\\")) or re.match(r"^[A-Za-z]:[\\/]", pattern) or ".." in Path(pattern).parts:
            raise ValueError(f"include pattern must be a repo-relative glob without traversal: {pattern!r}")

    files: dict[str, FileRecord] = {}
    rejected: list[str] = []
    review_required: list[str] = []
    matched_includes = {pattern: False for pattern in include_patterns}

    for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
        dir_path = Path(dirpath)
        safe_dirnames: list[str] = []
        for name in sorted(dirnames):
            if name in _SKIP_DIRNAMES:
                continue
            directory = dir_path / name
            if directory.is_symlink() or _is_reparse_point(directory):
                rel = directory.relative_to(root).as_posix()
                if any(fnmatch(rel, pattern) or fnmatch(f"{rel}/x", pattern) for pattern in include_patterns):
                    for pattern in include_patterns:
                        if fnmatch(rel, pattern) or fnmatch(f"{rel}/x", pattern):
                            matched_includes[pattern] = True
                    rule = "S1" if directory.is_symlink() else "S2"
                    files[rel] = FileRecord(
                        rel=rel, size=0, sha256="", classification="reject", rules=(rule,),
                        reason="path rules: link or reparse-point directory is not allowed inside a package",
                        kind="", provenance="repo", allowlist_rule=None,
                    )
                    rejected.append(rel)
                continue
            safe_dirnames.append(name)
        dirnames[:] = safe_dirnames
        for name in sorted(filenames):
            path = dir_path / name
            rel = path.relative_to(root).as_posix()
            for pattern in include_patterns:
                if fnmatch(rel, pattern):
                    matched_includes[pattern] = True

            classification = "reject"
            rules: tuple[str, ...] = ()
            reason = ""
            kind = ""
            allowlist_rule: str | None = None
            executor_input_schema: str | None = None
            size = 0
            digest = ""

            if not any(fnmatch(rel, pattern) for pattern in include_patterns):
                classification, rules, reason = "review_required", ("R1",), "allowlist rules: not declared by any include pattern"
            elif any(fnmatch(rel, deny) for deny in DENY_PATH_PATTERNS):
                rules, reason = ("P1",), "path rules: path matches a deny pattern (truth/scoring artifact)"
            elif path.is_symlink():
                rules, reason = ("S1",), "path rules: symlink is not allowed inside a package"
            elif _is_reparse_point(path):
                rules, reason = ("S2",), "path rules: reparse point is not allowed inside a package"
            elif not _resolved_under(path, root):
                rules, reason = ("S3",), "path rules: resolved target escapes the source root"
            else:
                rule = _match_allowlist(rel, include_patterns)
                if rule is None:
                    classification, rules, reason = "review_required", ("R1",), "allowlist rules: not declared by any include pattern"
                else:
                    classification = "allow"
                    kind = rule.kind
                    allowlist_rule = rule.pattern
                    if rel.endswith((".tsv", ".txt")):
                        hits = _truth_shaped_columns(_tsv_header(path))
                        if hits:
                            classification, rules, reason = (
                                "reject",
                                ("C1",),
                                f"column rules: truth-derived column(s) {hits} in TSV header",
                            )
                    elif rel.endswith(".json"):
                        document = _json_document(path)
                        if document is not None:
                            violations = sorted(set(_iter_json_keys(document)) & TRUTH_JSON_KEYS)
                            if violations:
                                classification, rules, reason = (
                                    "reject",
                                    ("C2",),
                                    f"content rules: truth/scoring JSON key(s) {violations}",
                                )
                            elif _json_has_post_freeze_phase(document):
                                classification, rules, reason = (
                                    "reject",
                                    ("C2",),
                                "content rules: JSON carries a post_freeze_scoring phase",
                            )
                    if classification == "allow" and (kind == "executor_data" or _is_executor_input_rel(rel)):
                        schema, schema_error = _executor_input_schema(rel, path)
                        if schema_error:
                            classification, rules, reason = "reject", ("C3",), schema_error
                        else:
                            value_error = _truth_value_marker(path)
                            if value_error:
                                classification, rules, reason = "reject", ("C4",), value_error
                            else:
                                executor_input_schema = schema
                    else:
                        executor_input_schema = None

            if classification == "allow":
                try:
                    size = path.stat().st_size
                except OSError:
                    size = 0
                digest = sha256(path)
            elif classification == "reject":
                rejected.append(rel)
            else:
                review_required.append(rel)

            files[rel] = FileRecord(
                rel=rel,
                size=size,
                sha256=digest,
                classification=classification,
                rules=rules,
                reason=reason,
                kind=kind,
                provenance="repo",
                allowlist_rule=allowlist_rule,
                executor_input_schema=executor_input_schema,
            )

    casefold_map: dict[str, list[str]] = {}
    for rel in files:
        casefold_map.setdefault(rel.casefold(), []).append(rel)
    for collisions in casefold_map.values():
        if len(collisions) <= 1:
            continue
        for rel in collisions:
            record = files[rel]
            if record.classification == "reject":
                continue
            others = sorted(other for other in collisions if other != rel)
            files[rel] = FileRecord(
                rel=rel,
                size=record.size,
                sha256="",
                classification="reject",
                rules=("S5",),
                reason=f"path rules: case-collision with {others}",
                kind="",
                provenance="repo",
                allowlist_rule=None,
                executor_input_schema=None,
            )
            if rel not in rejected:
                rejected.append(rel)

    unmatched = [pattern for pattern, hit in matched_includes.items() if not hit]
    return ScanResult(files=files, rejected=rejected, review_required=review_required, unmatched_includes=unmatched)


def _manifest_record(record: FileRecord) -> dict[str, object]:
    return {
        "size": record.size,
        "sha256": record.sha256,
        "classification": record.classification,
        "rules": list(record.rules),
        "reason": record.reason,
        "kind": record.kind,
        "provenance": record.provenance,
        "allowlist_rule": record.allowlist_rule,
        "executor_input_schema": record.executor_input_schema,
    }


def _package_tree_digest(files: dict[str, dict[str, object]]) -> str:
    """Digest the immutable Executor-visible file inventory, not host paths."""
    entries = [
        {"rel": rel, "size": record["size"], "sha256": record["sha256"]}
        for rel, record in sorted(files.items())
    ]
    encoded = json.dumps(entries, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def fingerprint_runtime() -> dict[str, object]:
    freeze_sha: str | None = None
    try:
        result = subprocess.run([sys.executable, "-m", "pip", "freeze"], capture_output=True, text=True, timeout=60)
        if result.returncode == 0:
            freeze_sha = hashlib.sha256(result.stdout.encode("utf-8")).hexdigest()
    except (OSError, subprocess.SubprocessError):
        pass
    return {
        "os": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
            "platform": platform.platform(),
        },
        "python": {
            "version": platform.python_version(),
            "implementation": platform.python_implementation(),
            "executable": sys.executable,
        },
        "package": {"candidate-qtl": {"version": __version__}},
        "dependency_lock": {"pip_freeze_sha256": freeze_sha},
    }


def git_commit_info(repo_root: Path) -> dict[str, object]:
    base = ["git", "-C", str(repo_root)]
    dot_git = repo_root / ".git"
    if dot_git.is_file():
        content = dot_git.read_text(encoding="utf-8").strip()
        if content.startswith("gitdir:"):
            git_dir = runtime_path(content.split(":", 1)[1].strip())
            base = ["git", "--git-dir", str(git_dir), "--work-tree", str(repo_root)]

    def run(*args: str) -> str:
        result = subprocess.run([*base, *args], capture_output=True, text=True)
        if result.returncode != 0:
            raise ValueError(f"git {args[0]} failed in {repo_root}: {result.stderr.strip()}")
        return result.stdout

    commit = run("rev-parse", "HEAD").strip()
    branch = run("branch", "--show-current").strip()
    status = run("status", "--porcelain")
    modified: list[str] = []
    untracked: list[str] = []
    for line in status.splitlines():
        if not line:
            continue
        code = line[:2]
        path_part = line[3:].strip()
        if code.startswith("??"):
            untracked.append(path_part)
        else:
            modified.append(path_part)
    tree = run("ls-tree", "-r", "HEAD")
    tree_fingerprint = hashlib.sha256(tree.encode("utf-8")).hexdigest()
    return {
        "commit": commit,
        "branch": branch,
        "dirty": bool(modified or untracked),
        "status": {"modified": sorted(modified), "untracked": sorted(untracked)},
        "tree_fingerprint": tree_fingerprint,
    }


def _path_host(path_wsl: str) -> str | None:
    match = re.match(r"^/mnt/([a-z])/(.*)$", path_wsl)
    if not match:
        return None
    return f"{match.group(1)}:\\{match.group(2).replace('/', chr(92))}"


def fingerprint_reference_assets(tsv_path: Path, full_hash: bool = False) -> dict[str, object]:
    if not tsv_path.is_file():
        raise ValueError(f"reference assets TSV not found: {tsv_path}")
    assets: dict[str, object] = {}
    with tsv_path.open("r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "path_wsl" not in reader.fieldnames:
            raise ValueError("reference assets TSV must have a 'path_wsl' column")
        for row in reader:
            raw = (row.get("path_wsl") or "").strip()
            if not raw:
                continue
            target = runtime_path(raw)
            key = (row.get("resource_key") or raw).strip()
            entry: dict[str, object] = {
                "path_policy": "external_absolute",
                "path_wsl": raw,
                "path_host": _path_host(raw),
                "sha256": None,
                "size": None,
                "fingerprint": None,
                "hash_status": "missing_at_build",
                "provider": row.get("provider") or None,
                "version": row.get("version") or None,
                "retrieved_at_utc": row.get("retrieved_at_utc") or None,
                "declared_in": row.get("declared_in") or str(tsv_path),
            }
            if target.is_file():
                entry["sha256"] = sha256(target)
                entry["size"] = target.stat().st_size
                entry["fingerprint"] = file_fingerprint(target, full_hash)
                entry["hash_status"] = "hashed_at_build"
            assets[key] = entry
    return assets


def _build_capsule(
    root: Path,
    out_dir: Path,
    start_utc: str,
    run_command: list[str] | None,
    reference_assets: Path | None,
    full_hash: bool,
    file_count: int,
    manifest_sha: str,
    package_tree_sha: str,
    build_run_id: str,
) -> dict[str, object]:
    resources: dict[str, object] = {}
    if reference_assets is not None:
        resources = fingerprint_reference_assets(runtime_path(reference_assets), full_hash)
    return {
        "schema_version": CAPSULE_VERSION,
        "capsule_kind": "executor_build",
        "built_at_utc": _now(),
        "candidate_qtl_version": __version__,
        "build": {
            "command": "build-executor-package",
            "argv": list(sys.argv),
            "declared_executor_command": list(run_command) if run_command else None,
            "start_utc": start_utc,
            "end_utc": _now(),
            "threads": {"os_cpu_count": os.cpu_count() or 1, "used": 1},
            "environment_allowlist": {key: value for key, value in os.environ.items() if key.startswith("CANDIDATE_QTL_")},
        },
        "git": git_commit_info(root),
        "runtime": fingerprint_runtime(),
        "resources": resources,
        "executor_package": {"manifest": "package_manifest.json", "sha256": manifest_sha, "file_count": file_count},
        "build_run_id": build_run_id,
        "package_tree_sha256": package_tree_sha,
        "truth_inputs_present": False,
        "policy": [
            "env is recorded only for allowlisted variable names; the full environment is never dumped",
            "argv is recorded verbatim; do not pass credentials as command-line arguments",
            "declared reference assets are trusted Supervisor input; they are fingerprinted in place and never copied into the package",
            "a resource missing at build time is recorded missing_at_build, not a build failure",
            "supervisor scoring opens truth only after freeze verification (existing audit-score contract)",
        ],
    }


def build_executor_package(
    repo_root: Path,
    out_dir: Path,
    include_patterns: list[str] | None = None,
    reference_assets: Path | None = None,
    full_hash: bool = False,
    run_command: list[str] | None = None,
) -> BuildResult:
    root = Path(repo_root).resolve()
    out = Path(out_dir).resolve()
    if out == root or out.is_relative_to(root):
        raise ValueError(f"output directory must be outside the repository root: {out}")
    if out.exists():
        if any(out.iterdir()):
            raise ValueError(f"output directory must not exist or be empty: {out}")
    else:
        out.mkdir(parents=True)

    patterns = list(include_patterns) if include_patterns is not None else [rule.pattern for rule in DEFAULT_ALLOWLIST]
    scan = scan_package_source(root, patterns)
    default_patterns = {rule.pattern for rule in DEFAULT_ALLOWLIST}
    user_unmatched = [pattern for pattern in scan.unmatched_includes if pattern not in default_patterns]
    if user_unmatched:
        shutil.rmtree(out, ignore_errors=True)
        raise ValueError(f"include pattern(s) matched no files: {sorted(user_unmatched)}")
    if scan.rejected:
        lines = "\n".join(_format_record(scan.files[rel]) for rel in sorted(scan.rejected))
        shutil.rmtree(out, ignore_errors=True)
        raise ValueError(f"executor package rejected {len(scan.rejected)} file(s):\n{lines}")

    start_utc = _now()
    allowed = [rel for rel, record in scan.files.items() if record.classification == "allow"]
    try:
        for rel in allowed:
            source = root / rel
            target = out / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
    except Exception:
        shutil.rmtree(out, ignore_errors=True)
        raise

    manifest_records = {rel: _manifest_record(record) for rel, record in scan.files.items() if record.classification == "allow"}
    package_tree_sha = _package_tree_digest(manifest_records)
    build_run_id = str(uuid.uuid4())
    manifest = {
        "schema_version": PACKAGE_MANIFEST_VERSION,
        "package_version": PACKAGE_VERSION,
        "package_kind": "executor",
        "built_at_utc": _now(),
        "candidate_qtl_version": __version__,
        "source_root": str(root),
        "truth_inputs_present": False,
        "files": manifest_records,
        "file_count": len(allowed),
        "package_tree_sha256": package_tree_sha,
        "build_run_id": build_run_id,
        "total_bytes": sum(record.size for record in scan.files.values() if record.classification == "allow"),
        "rules": patterns,
        "reviews": [
            {"rel": rel, "rules": list(scan.files[rel].rules), "reason": scan.files[rel].reason}
            for rel in sorted(scan.review_required)
        ],
        "run_capsule": {"file": "run_capsule.json"},
        "policy": [
            "allowlist is the primary contract; the negative scanner is the fail-closed backstop for any --include",
            "a rejected file aborts the build; undeclared files are excluded and recorded as review_required",
            "truth, source, score and freeze artifacts are never packaged (truth_inputs_present=false)",
            "relative paths use / separators; the package is relocatable",
        ],
    }
    _write_manifest(out / "package_manifest.json", manifest)

    manifest_sha = sha256(out / "package_manifest.json")
    capsule = _build_capsule(
        root=root,
        out_dir=out,
        start_utc=start_utc,
        run_command=run_command,
        reference_assets=reference_assets,
        full_hash=full_hash,
        file_count=len(allowed),
        manifest_sha=manifest_sha,
        package_tree_sha=package_tree_sha,
        build_run_id=build_run_id,
    )
    _write_manifest(out / "run_capsule.json", capsule)

    return BuildResult(
        out_dir=out,
        manifest_path=out / "package_manifest.json",
        capsule_path=out / "run_capsule.json",
        file_count=len(allowed),
        total_bytes=manifest["total_bytes"],
        rejected=scan.rejected,
        review_required=scan.review_required,
        truth_inputs_present=False,
    )


def _regular_file_under(path: Path, root: Path) -> bool:
    try:
        return path.is_file() and not path.is_symlink() and _resolved_under(path, root.resolve())
    except OSError:
        return False


def _load_json_object(path: Path, label: str) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is not a readable JSON object: {path}") from error
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a JSON object: {path}")
    return value


def _verify_packaged_inventory(package_root: Path, manifest: dict[str, object]) -> str:
    files = manifest.get("files")
    if not isinstance(files, dict) or not files:
        raise ValueError("executor package manifest has no non-empty files inventory")
    if manifest.get("file_count") != len(files):
        raise ValueError("executor package manifest file_count does not match its files inventory")
    normalized: dict[str, dict[str, object]] = {}
    for rel, raw_record in files.items():
        if not isinstance(rel, str) or not isinstance(raw_record, dict):
            raise ValueError("executor package manifest has an invalid files entry")
        candidate = package_root / rel
        if Path(rel).is_absolute() or ".." in Path(rel).parts or not _regular_file_under(candidate, package_root):
            raise ValueError(f"executor package inventory path is unsafe or missing: {rel}")
        expected_hash = raw_record.get("sha256")
        expected_size = raw_record.get("size")
        if not isinstance(expected_hash, str) or not _SHA256.fullmatch(expected_hash) or not isinstance(expected_size, int):
            raise ValueError(f"executor package inventory record is invalid: {rel}")
        if candidate.stat().st_size != expected_size or sha256(candidate) != expected_hash:
            raise ValueError(f"executor package file changed after build: {rel}")
        normalized[rel] = raw_record
    tree_sha = _package_tree_digest(normalized)
    if manifest.get("package_tree_sha256") != tree_sha:
        raise ValueError("executor package tree digest mismatch")
    return tree_sha


def _verify_frozen_predictions(freeze: dict[str, object], freeze_path: Path) -> int:
    if freeze.get("phase") != "prediction_only" or freeze.get("truth_inputs_present") is not False:
        raise ValueError("Supervisor requires a prediction_only freeze with truth_inputs_present=false")
    cases = freeze.get("cases")
    if not isinstance(cases, dict) or not cases:
        raise ValueError("prediction freeze has no non-empty case inventory")
    if freeze.get("case_count") != len(cases):
        raise ValueError("prediction freeze case_count does not match its case inventory")
    predictions = freeze_path.parent / "predictions"
    expected_paths: set[Path] = set()
    for case_id, raw_record in cases.items():
        if not isinstance(case_id, str) or not case_id or not isinstance(raw_record, dict):
            raise ValueError("prediction freeze has an invalid case record")
        declared = raw_record.get("ranking_tsv")
        digest = raw_record.get("sha256")
        if not isinstance(declared, str) or not isinstance(digest, str) or not _SHA256.fullmatch(digest):
            raise ValueError(f"prediction freeze has no valid ranking hash for case {case_id!r}")
        path = runtime_path(declared)
        if not _regular_file_under(path, freeze_path.parent):
            raise ValueError(f"prediction freeze ranking is outside the freeze or missing: {case_id!r}")
        if path.stat().st_size == 0 or sha256(path) != digest:
            raise ValueError(f"prediction freeze ranking changed after freeze: {case_id!r}")
        expected_paths.add(path.resolve())
    if not predictions.is_dir():
        raise ValueError("prediction freeze has no predictions directory")
    actual_paths = {path.resolve() for path in predictions.rglob("*") if path.is_file()}
    if actual_paths != expected_paths:
        raise ValueError("prediction freeze output inventory has missing or unrecorded ranking files")
    return len(cases)


def verify_supervisor_ready(freeze_path: Path, capsule_path: Path) -> dict[str, object]:
    """Fail closed before a Supervisor may open truth for an Executor run.

    The freeze must carry the package and capsule identities emitted at build
    time.  Merely placing a freeze next to a package is never a valid binding.
    """
    capsule_path = runtime_path(capsule_path).resolve()
    package_root = capsule_path.parent
    if not _regular_file_under(capsule_path, package_root):
        raise ValueError(f"executor run capsule is unsafe or missing: {capsule_path}")
    capsule = _load_json_object(capsule_path, "executor run capsule")
    package = capsule.get("executor_package")
    if capsule.get("capsule_kind") != "executor_build" or capsule.get("truth_inputs_present") is not False or not isinstance(package, dict):
        raise ValueError("executor run capsule does not declare a truth-free executor build")
    if package.get("manifest") != "package_manifest.json":
        raise ValueError("executor run capsule must reference the fixed package_manifest.json")
    manifest_path = package_root / "package_manifest.json"
    if not _regular_file_under(manifest_path, package_root):
        raise ValueError("executor package manifest is unsafe or missing")
    actual_manifest_sha = sha256(manifest_path)
    if package.get("sha256") != actual_manifest_sha:
        raise ValueError("executor package manifest hash mismatch")
    manifest = _load_json_object(manifest_path, "executor package manifest")
    if manifest.get("package_kind") != "executor" or manifest.get("truth_inputs_present") is not False:
        raise ValueError("executor package manifest does not declare a truth-free executor package")
    tree_sha = _verify_packaged_inventory(package_root, manifest)
    if capsule.get("package_tree_sha256") != tree_sha or capsule.get("build_run_id") != manifest.get("build_run_id"):
        raise ValueError("executor run capsule is not bound to this exact package build")

    raw_freeze = runtime_path(freeze_path)
    freeze_manifest = raw_freeze if raw_freeze.name == "freeze_manifest.json" else raw_freeze / "freeze_manifest.json"
    if not _regular_file_under(freeze_manifest, freeze_manifest.parent):
        raise ValueError(f"prediction freeze not found at {freeze_manifest}")
    freeze = _load_json_object(freeze_manifest, "prediction freeze")
    case_count = _verify_frozen_predictions(freeze, freeze_manifest)
    binding = freeze.get("executor_package")
    if not isinstance(binding, dict):
        raise ValueError("prediction freeze has no executor package binding")
    capsule_sha = sha256(capsule_path)
    if (
        binding.get("manifest_sha256") != actual_manifest_sha
        or binding.get("package_tree_sha256") != tree_sha
        or binding.get("capsule_sha256") != capsule_sha
        or binding.get("build_run_id") != capsule.get("build_run_id")
    ):
        raise ValueError("prediction freeze is not bound to this exact executor package and run capsule")
    return {
        "freeze_verified": True,
        "manifest_sha256": actual_manifest_sha,
        "package_tree_sha256": tree_sha,
        "capsule_sha256": capsule_sha,
        "case_count": case_count,
        "capsule_matches_package": True,
    }
