"""Explicit optional copy inputs, not a second evidence evaluator (#729).

Only byte identity and routing are handled here. Native family-to-copy schemas,
source/target verification and missing-data semantics remain authoritative.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from .core import runtime_path

COPY_EVIDENCE_AXES = ("orthology", "synteny", "expression", "position", "subfamily_map", "target")
CONTRACT = "canonical_copy_inputs_v1"


def inspect_copy_evidence(values: dict[str, str | None] | None) -> tuple[dict[str, str], dict[str, object]]:
    """Validate explicit paths and hash full contents; never discover local data."""
    if values is None:
        return {}, {}
    if not isinstance(values, dict) or set(values) - set(COPY_EVIDENCE_AXES):
        raise ValueError("copy evidence requires only the six named native axes")
    kwargs: dict[str, str] = {}
    inputs: dict[str, object] = {}
    for axis in COPY_EVIDENCE_AXES:
        value = values.get(axis)
        if value is None:
            continue
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"copy {axis} must be an explicit non-empty file path")
        path = runtime_path(value).resolve()
        if not path.is_file():
            raise ValueError(f"supplied copy {axis} file does not exist: {path}")
        sha = hashlib.sha256()
        size = 0
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                size += len(block)
                sha.update(block)
        kwargs[axis + "_arg"] = str(path)
        inputs[axis] = {"path": str(path), "sha256": sha.hexdigest(), "size_bytes": size}
    return kwargs, inputs


def check_copy_evidence_binding(root: Path, inputs: dict[str, object]) -> None:
    """Read-only preflight before the outer profile can modify existing outputs."""
    path = root / "b04_canonical_semantic" / "canonical_copy_inputs.json"
    if path.exists():
        record = {"contract": CONTRACT, "inputs": inputs}
        if json.loads(path.read_text(encoding="utf-8")) != record:
            raise ValueError("copy evidence changed or was omitted on resume; use a new output revision")
    elif inputs and (root / "z04_canonical_wet_handoff" / "wet_candidates.tsv").exists():
        raise ValueError("existing Wet output has no copy-input pin; use a new output revision")


def bind_copy_evidence(root: Path, values: dict[str, str | None] | None) -> tuple[dict[str, str], dict[str, object]]:
    """Pin supplied axes across continuation; changed inputs need a new output revision."""
    kwargs, inputs = inspect_copy_evidence(values)
    check_copy_evidence_binding(root, inputs)
    path = root / "b04_canonical_semantic" / "canonical_copy_inputs.json"
    record = {"contract": CONTRACT, "inputs": inputs}
    if inputs and not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        try:
            with path.open("x", encoding="utf-8") as stream:
                json.dump(record, stream, ensure_ascii=False, indent=2, sort_keys=True)
                stream.write("\n")
        except FileExistsError:
            if json.loads(path.read_text(encoding="utf-8")) != record:
                raise ValueError("concurrent copy evidence mismatch; use a new output revision") from None
    if not inputs:
        return {}, {}
    return kwargs, {"copy_evidence": {
        "contract": CONTRACT, "input_manifest": str(path), "inputs": inputs,
        "omitted_axes": [axis for axis in COPY_EVIDENCE_AXES if axis not in inputs],
        "consumer_state": "not_evaluated_before_family_to_copy",
        "fingerprint_is_biological_verification": False,
        "missing_is_negative": False,
    }}


def verify_copy_evidence(kwargs: dict[str, str], provenance: dict[str, object]) -> None:
    """Detect a supplied input changing while the downstream stage is running."""
    if not kwargs:
        return
    _, current = inspect_copy_evidence({key.removesuffix("_arg"): value for key, value in kwargs.items()})
    if current != provenance["copy_evidence"]["inputs"]:
        raise ValueError("copy evidence changed during execution; do not use the resulting artifacts")
