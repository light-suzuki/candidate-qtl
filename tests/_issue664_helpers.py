from __future__ import annotations

import hashlib
import json
from pathlib import Path


def manifest_for(
    path: Path,
    *,
    route: str,
    case_id: str = "case1",
    qtl_id: str = "qtl1",
    interval_fingerprint: str = "interval1",
    assembly: str = "assembly1",
    manifest_dir: Path | None = None,
) -> Path:
    """Create a producer manifest whose declared bytes are real and sealed.

    Route-C fixtures include actual ancestor files rather than placeholder
    digests, so tests exercise the same lineage shape as production output.
    """
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    work_dir = manifest_dir or path.parent
    if route == "a":
        ancestor = work_dir / f"{path.stem}_route_a_input.fixture"
        ancestor.write_text("hypothesis_id\nH1\n", encoding="utf-8")
        payload = {"schema_version": "1.0", "command": "family-candidate-bridge", "contract": "family-candidate-bridge-1.0.0", "case_id": case_id,
                   "qtl_id": qtl_id, "interval_fingerprint": interval_fingerprint, "assembly": assembly,
                   "inputs": {"route_a_hypotheses": {"present": True, "path": str(ancestor),
                                                        "sha256": hashlib.sha256(ancestor.read_bytes()).hexdigest(),
                                                        "strength": "strong"}},
                   "outputs": {"family_candidate_bridge": str(path), "family_candidate_bridge_tsv_sha256": digest}}
        name = f"{path.stem}_manifest.json"
    elif route in {"b", "b_claims"}:
        claims = route == "b_claims"
        payload = {"schema_version": "1.0", "stage": "route_b_semantic_transfer_verification" if claims else "route_b_verified_candidate_freeze",
                   "contract": "route-b-semantic-transfer-verification-1.0.0" if claims else "route-b-verified-candidates-1.0.0", "case_id": case_id,
                   "qtl_id": qtl_id, "interval_fingerprint": interval_fingerprint, "assembly": assembly,
                   "artifact_sha256": digest, "output": str(path)}
        if claims:
            payload["outputs"] = {"claims": str(path), "claims_sha256": digest}
            ancestor_paths = {}
            for ancestor_name in ("source_verification", "proposals", "bundle"):
                ancestor = work_dir / f"{path.stem}_{ancestor_name}.fixture"
                ancestor.write_text(f"{ancestor_name}\n", encoding="utf-8")
                ancestor_paths[ancestor_name] = {
                    "path": str(ancestor),
                    "sha256": hashlib.sha256(ancestor.read_bytes()).hexdigest(),
                }
            payload["ancestors"] = ancestor_paths
        else:
            claims_path = work_dir / f"{path.stem}_verified_claims.fixture"
            claims_path.write_text("claim_id\nC1\n", encoding="utf-8")
            claims_manifest_path = work_dir / f"{path.stem}_verified_claims_manifest.json"
            claims_manifest_path.write_text(json.dumps({
                "schema_version": "1.0",
                "stage": "route_b_semantic_transfer_verification",
                "contract": "route-b-semantic-transfer-verification-1.0.0",
                "case_id": case_id,
                "artifact_sha256": hashlib.sha256(claims_path.read_bytes()).hexdigest(),
                "outputs": {"claims": str(claims_path), "claims_sha256": hashlib.sha256(claims_path.read_bytes()).hexdigest()},
            }, sort_keys=True) + "\n", encoding="utf-8")
            payload["inputs"] = {
                "claims": {"path": str(claims_path), "sha256": hashlib.sha256(claims_path.read_bytes()).hexdigest()},
                "claims_manifest": {"path": str(claims_manifest_path), "sha256": hashlib.sha256(claims_manifest_path.read_bytes()).hexdigest()},
            }
        name = f"{path.stem}_manifest.json"
    elif route == "c":
        frozen = path.with_name("frozen.tsv")
        projection = path.with_name("projection.tsv")
        frozen.write_text("record_id\nR1\n", encoding="utf-8")
        projection.write_text("record_id\tqtl_id\nR1\t" + qtl_id + "\n", encoding="utf-8")
        frozen_sha = hashlib.sha256(frozen.read_bytes()).hexdigest()
        projection_sha = hashlib.sha256(projection.read_bytes()).hexdigest()
        freeze_manifest = path.with_name("freeze_manifest.json")
        freeze_manifest.write_text(json.dumps({
            "schema_version": "1.0", "command": "route-c-direct-freeze",
            "case_id": case_id, "frozen_sha256": frozen_sha,
        }, sort_keys=True) + "\n", encoding="utf-8")
        projection_manifest = path.with_name("projection_manifest.json")
        projection_manifest.write_text(json.dumps({
            "schema_version": "1.0", "command": "route-c-direct-project",
            "case_id": case_id, "qtl_id": qtl_id,
            "frozen_sha256": frozen_sha, "projection_sha256": projection_sha,
        }, sort_keys=True) + "\n", encoding="utf-8")
        ancestors = {
            name: {"path": str(ancestor), "sha256": hashlib.sha256(ancestor.read_bytes()).hexdigest()}
            for name, ancestor in {
                "frozen": frozen,
                "freeze_manifest": freeze_manifest,
                "projection": projection,
                "projection_manifest": projection_manifest,
            }.items()
        }
        payload = {"schema_version": "1.1", "contract": "route-c-direct-members-1.1.0", "stage": "route_c_direct_members",
                   "case_id": case_id, "qtl_id": qtl_id, "interval_fingerprint": interval_fingerprint,
                   "assembly": assembly, "output": str(path), "output_sha256": digest, "inputs": ancestors}
        name = f"{path.stem}_manifest.json"
    else:
        raise ValueError(route)
    out = work_dir / name
    out.write_text(json.dumps(payload, sort_keys=True) + "\n", encoding="utf-8")
    return out
