"""Route-C direct member ledger with projection-conflict fail-closed handling.

The v0.22 member builder remains the compatibility implementation. This facade
verifies the issue-#663 projection-evidence ledger and prevents an exact-ID row
from reaching canonical support when another projection method maps the same
source record to a different focal gene.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from .core import runtime_path, sha256
from .route_c_direct_members_legacy import (
    CONTRACT,
    LEDGER_STATUS,
    OUTPUT_COLUMNS,
    REVIEW_STATUS,
    SCHEMA_VERSION,
    SUPPORTED_STATUS,
    build_route_c_direct_members as _legacy_build_route_c_direct_members,
)
from .route_c_direct_projection import (
    PROJECTION_EVIDENCE_CONTRACT,
    validate_projection_evidence_frame,
)
from .route_c_identity_mapping import load_projection_identity_mapping


def _resolve_output_path(value: object, manifest_path: Path) -> Path:
    path = runtime_path(str(value or ""))
    if path.is_absolute():
        return path
    # Writers record outputs relative to the working directory they ran in
    # (``--out-dir c-project`` stores ``c-project/x.tsv``), so a plain
    # manifest-relative join doubles the directory. Prefer the
    # manifest-relative location, then the same-directory sibling, then the
    # recorded working-directory path; callers still verify the sealed hash.
    candidates = (manifest_path.parent / path, manifest_path.parent / path.name, path)
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def _verified_mapping_conflicts(
    projection_manifest_arg: str,
) -> tuple[set[str], dict[str, object] | None]:
    manifest_path = runtime_path(projection_manifest_arg)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    from .reference_bridge_consumers import verify_projection_reference_bridge
    verify_projection_reference_bridge(manifest, manifest_path)
    if manifest.get("inputs", {}).get("identity_mapping") is not None:
        frozen_descriptor = manifest["inputs"]["frozen"]
        frozen_path = _resolve_output_path(frozen_descriptor["path"], manifest_path)
        if not frozen_path.is_file() or sha256(frozen_path) != manifest.get("frozen_sha256"):
            raise ValueError("Route C identity mapping frozen artifact hash mismatch")
        frozen = pd.read_csv(frozen_path, sep="\t", dtype=str, keep_default_na=False).fillna("")
        load_projection_identity_mapping(
            manifest, manifest_path,
            {str(row["record_id"]): row for row in frozen.to_dict("records")},
        )
    contract = str(manifest.get("projection_evidence_contract", "")).strip()
    if not contract:
        # Backward-compatible old projection generations have no ledger. They
        # retain their historical behavior and cannot claim issue-#663 sealing.
        return set(), None
    if contract != PROJECTION_EVIDENCE_CONTRACT:
        raise ValueError(
            f"unsupported Route C projection-evidence contract: {contract}; "
            f"expected {PROJECTION_EVIDENCE_CONTRACT}"
        )
    outputs = manifest.get("outputs", {})
    if not isinstance(outputs, dict):
        raise ValueError("Route C projection manifest outputs must be an object")
    artifact = _resolve_output_path(outputs.get("projection_evidence", ""), manifest_path)
    expected_sha = (
        str(manifest.get("projection_evidence_sha256", ""))
        .strip()
        .removeprefix("sha256:")
    )
    if not artifact.exists() or not expected_sha:
        raise ValueError("Route C projection-evidence artifact or digest is missing")
    if sha256(artifact) != expected_sha:
        raise ValueError(
            "Route C projection-evidence artifact hash does not match projection manifest"
        )
    evidence = validate_projection_evidence_frame(
        pd.read_csv(artifact, sep="\t", dtype=str, keep_default_na=False)
    )
    computed = set(
        evidence.loc[
            evidence["conflict_state"].eq("mapping_conflict"), "record_id"
        ].astype(str)
    )
    declared_raw = manifest.get("mapping_conflict_record_ids", [])
    if not isinstance(declared_raw, list):
        raise ValueError("Route C mapping_conflict_record_ids must be a list")
    declared = {str(value) for value in declared_raw}
    if computed != declared:
        raise ValueError(
            "Route C mapping-conflict roster does not match projection evidence"
        )
    try:
        declared_count = int(manifest.get("mapping_conflict_count", -1))
    except (TypeError, ValueError) as error:
        raise ValueError("Route C mapping_conflict_count must be an integer") from error
    if declared_count != len(computed):
        raise ValueError(
            "Route C mapping_conflict_count does not match projection evidence"
        )
    return computed, {
        "path": str(artifact),
        "sha256": expected_sha,
        "contract": contract,
    }


def build_route_c_direct_members(
    frozen_arg: str,
    freeze_manifest_arg: str,
    projection_arg: str,
    projection_manifest_arg: str,
    atlas_arg: str,
    out_dir_arg: str,
    *,
    case_id: str = "",
) -> dict[str, object]:
    # The compatibility builder writes member artifacts. Authenticate the
    # projection ledger first so invalid evidence cannot create or overwrite
    # valid-looking output before this facade fails closed.
    conflict_ids, evidence_descriptor = _verified_mapping_conflicts(
        projection_manifest_arg
    )
    projection_manifest_path = runtime_path(projection_manifest_arg)
    projection_manifest = json.loads(projection_manifest_path.read_text(encoding="utf-8"))
    from .reference_bridge_consumers import verify_projection_reference_bridge
    bridge_pairs = verify_projection_reference_bridge(projection_manifest, projection_manifest_path)
    result = _legacy_build_route_c_direct_members(
        frozen_arg,
        freeze_manifest_arg,
        projection_arg,
        projection_manifest_arg,
        atlas_arg,
        out_dir_arg,
        case_id=case_id,
    )
    if evidence_descriptor is None:
        return result

    members_path = runtime_path(str(result["members"]))
    manifest_path = runtime_path(str(result["manifest"]))
    frame = pd.read_csv(members_path, sep="\t", dtype=str, keep_default_na=False)
    if list(frame.columns) != list(OUTPUT_COLUMNS):
        raise ValueError("Route C direct members output drifted from OUTPUT_COLUMNS")
    block = (
        (frame["record_id"].astype(str).isin(conflict_ids)
         | pd.Series([(str(row.record_id), str(row.candidate_gene_id)) in bridge_pairs
                      for row in frame.itertuples()], index=frame.index))
        & frame["support_status"].eq(SUPPORTED_STATUS)
    )
    blocked_count = int(block.sum())
    if blocked_count:
        frame.loc[block, "support_status"] = REVIEW_STATUS
        frame.loc[block, "membership_effect"] = "review_only"
        suffix = (
            "projection conflict or sequence correspondence requires review; canonical support "
            "is blocked pending source-bound locus/model and phenotype assessment"
        )
        frame.loc[block, "notes"] = frame.loc[block, "notes"].astype(str).map(
            lambda value: "; ".join(
                part for part in (value.strip(), suffix) if part
            )
        )
        frame.to_csv(
            members_path, sep="\t", index=False, lineterminator="\n"
        )

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["supported_count"] = int(
        frame["support_status"].eq(SUPPORTED_STATUS).sum()
    )
    manifest["review_count"] = int(
        frame["support_status"].eq(REVIEW_STATUS).sum()
    )
    manifest["ledger_only_count"] = int(
        frame["support_status"].eq(LEDGER_STATUS).sum()
    )
    manifest["mapping_conflict_record_ids"] = sorted(
        conflict_ids, key=str.casefold
    )
    manifest["mapping_conflict_count"] = len(conflict_ids)
    manifest["projection_conflict_support_blocked_count"] = blocked_count
    manifest["reference_bridge_review_pairs"] = [list(pair) for pair in sorted(bridge_pairs)]
    manifest["output_sha256"] = sha256(members_path)
    manifest.setdefault("inputs", {})["projection_evidence"] = evidence_descriptor
    policy = manifest.setdefault("policy", [])
    statement = (
        "a projection record whose exact/orthology/BLAST methods resolve to "
        "different focal targets is review-only even when one individual row "
        "has an exact basis"
    )
    if statement not in policy:
        policy.append(statement)
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        **result,
        "record_count": len(frame),
        "supported_count": manifest["supported_count"],
        "review_count": manifest["review_count"],
        "mapping_conflict_count": len(conflict_ids),
        "projection_conflict_support_blocked_count": blocked_count,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frozen", required=True)
    parser.add_argument("--freeze-manifest", required=True)
    parser.add_argument("--projection", required=True)
    parser.add_argument("--projection-manifest", required=True)
    parser.add_argument("--atlas", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--case-id", default="")
    args = parser.parse_args(argv)
    result = build_route_c_direct_members(
        args.frozen,
        args.freeze_manifest,
        args.projection,
        args.projection_manifest,
        args.atlas,
        args.out_dir,
        case_id=args.case_id,
    )
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))
    return 0


__all__ = [
    "CONTRACT",
    "LEDGER_STATUS",
    "OUTPUT_COLUMNS",
    "REVIEW_STATUS",
    "SCHEMA_VERSION",
    "SUPPORTED_STATUS",
    "build_route_c_direct_members",
    "main",
]


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
