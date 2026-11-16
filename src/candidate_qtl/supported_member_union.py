"""Deterministic supported-member union for production narrowing (#426).

This module bridges discovery into the Wet-candidate budget without inventing a
new score.  It takes:

* Route-A rows from the existing ``family_candidate_bridge.tsv``; rows that were
  explicitly ``tier6_not_evaluated`` are not treated as supported mappings, but
  ambiguous supported paralogs remain;
* final Route-B semantic claims; only ``verification_status=verified`` rows may
  contribute members.

The result is a union set, not a rank.  Agreement between routes is recorded as
provenance and never becomes an automatic numerical bonus.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

from .core import runtime_path

SCHEMA_VERSION = "1.0"
MANIFEST_CONTRACT = "supported-member-union-1.2.0"


def _load_sealed_input(artifact_arg: str, manifest_arg: str | None, *, label: str,
                       case_id: str = "") -> tuple[object, dict[str, object], str]:
    """Load an artifact only after its producer manifest authenticates it.

    The sibling-name lookup is a convenience for normal producer outputs; it
    never turns an unsealed TSV into a canonical input.
    """
    artifact = runtime_path(artifact_arg)
    if not artifact.exists():
        raise ValueError(f"{label} artifact does not exist: {artifact}")
    candidates = [runtime_path(manifest_arg)] if manifest_arg else [
        artifact.with_name(f"{artifact.stem}_manifest.json"),
        artifact.parent / f"{label.lower().replace(' ', '_')}_manifest.json",
    ]
    manifest_path = next((p for p in candidates if p and p.exists()), None)
    if manifest_path is None:
        raise ValueError(f"{label} requires a producer manifest (unsealed artifact rejected)")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} producer manifest is not valid JSON") from error
    if not isinstance(manifest, dict):
        raise ValueError(f"{label} producer manifest must be a JSON object")
    if not str(manifest.get("schema_version", "")).strip():
        raise ValueError(f"{label} producer manifest lacks schema_version")
    # Each producer has one authoritative output declaration. Never accept a
    # digest belonging to an unrelated input/ancestor in the same manifest.
    specs = {
        "Route A bridge": (("outputs", "family_candidate_bridge_tsv_sha256"), ("outputs", "family_candidate_bridge"), "family-candidate-bridge"),
        "Route B claims": (("outputs", "claims_sha256"), ("outputs", "claims"), None),
        "Route B candidates": (("artifact_sha256",), ("output",), "route_b_verified_candidate_freeze"),
        "Route C members": (("output_sha256",), ("output",), "route_c_direct_members"),
    }
    hash_path, output_path, required_stage = specs.get(label, (("output_sha256",), ("output",), None))
    value: object = manifest
    for key in hash_path:
        if not isinstance(value, dict): value = None; break
        value = value.get(key)
    expected = str(value or "").removeprefix("sha256:")
    if len(expected) != 64:
        raise ValueError(f"{label} producer manifest lacks its authoritative artifact SHA-256")
    output_value: object = manifest
    for key in output_path:
        if not isinstance(output_value, dict): output_value = None; break
        output_value = output_value.get(key)
    if not str(output_value or "").strip():
        raise ValueError(f"{label} producer manifest lacks its authoritative output path")
    declared_output = Path(str(output_value))
    if not declared_output.is_absolute():
        declared_output = manifest_path.parent / declared_output
    if declared_output.resolve() != artifact.resolve():
        raise ValueError(f"{label} artifact does not match producer manifest output")
    if required_stage and str(manifest.get("stage", manifest.get("command", ""))).strip() != required_stage:
        raise ValueError(f"{label} producer manifest has an unsupported contract")
    actual = hashlib.sha256(artifact.read_bytes()).hexdigest()
    if actual != expected:
        raise ValueError(f"{label} artifact hash does not match producer manifest")
    manifest_case = str(manifest.get("case_id", "")).strip()
    if not manifest_case:
        raise ValueError(f"{label} producer manifest requires non-empty case_id")
    if case_id and manifest_case and manifest_case != case_id:
        raise ValueError(f"{label} producer manifest case_id does not match requested case_id")
    return artifact, manifest, str(manifest_path)


def _verify_route_c_ancestors(manifest: dict[str, object], manifest_path: Path, case_id: str) -> None:
    inputs = manifest.get("inputs")
    if not isinstance(inputs, dict):
        raise ValueError("Route C members manifest lacks freeze/projection ancestry")
    resolved: dict[str, Path] = {}
    parsed_manifests: dict[str, dict[str, object]] = {}
    for name in ("frozen", "freeze_manifest", "projection", "projection_manifest"):
        entry = inputs.get(name)
        if not isinstance(entry, dict):
            raise ValueError(f"Route C ancestor {name} is missing")
        raw = str(entry.get("path", "")).strip()
        digest = str(entry.get("sha256", "")).removeprefix("sha256:")
        if not raw or len(digest) != 64:
            raise ValueError(f"Route C ancestor {name} lacks path or SHA-256")
        path = Path(raw)
        if not path.is_absolute(): path = manifest_path.parent / path
        if not path.exists() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f"Route C ancestor {name} hash mismatch or missing")
        resolved[name] = path
        if name.endswith("manifest"):
            try: ancestor = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as error: raise ValueError(f"Route C ancestor {name} is invalid JSON") from error
            if not str(ancestor.get("schema_version", "")).strip():
                raise ValueError(f"Route C ancestor {name} lacks schema_version")
            ancestor_case = str(ancestor.get("case_id", "")).strip()
            if ancestor_case and case_id and ancestor_case != case_id:
                raise ValueError(f"Route C ancestor {name} case_id mismatch")
            parsed_manifests[name] = ancestor

    freeze = parsed_manifests["freeze_manifest"]
    projection = parsed_manifests["projection_manifest"]
    if freeze.get("command") != "route-c-direct-freeze":
        raise ValueError("Route C freeze ancestor has an unsupported producer")
    if projection.get("command") != "route-c-direct-project":
        raise ValueError("Route C projection ancestor has an unsupported producer")
    if str(freeze.get("frozen_sha256", "")) != hashlib.sha256(resolved["frozen"].read_bytes()).hexdigest():
        raise ValueError("Route C freeze manifest does not bind the frozen artifact")
    if str(projection.get("frozen_sha256", "")) != hashlib.sha256(resolved["frozen"].read_bytes()).hexdigest():
        raise ValueError("Route C projection manifest does not bind the frozen artifact")
    if str(projection.get("projection_sha256", "")) != hashlib.sha256(resolved["projection"].read_bytes()).hexdigest():
        raise ValueError("Route C projection manifest does not bind the projection artifact")
    expected_qtl = str(manifest.get("qtl_id", "")).strip()
    if not expected_qtl or str(projection.get("qtl_id", "")).strip() != expected_qtl:
        raise ValueError("Route C qtl_id does not match the projection ancestor")


def _verify_route_b_ancestors(manifest: dict[str, object], manifest_path: Path) -> None:
    ancestors = manifest.get("ancestors")
    if not isinstance(ancestors, dict):
        raise ValueError("Route B claims producer manifest lacks source/transfer ancestor ledger")
    for name in ("source_verification", "proposals", "bundle"):
        entry = ancestors.get(name)
        if not isinstance(entry, dict):
            raise ValueError(f"Route B claims ancestor {name} is missing")
        raw = str(entry.get("path", "")).strip()
        digest = str(entry.get("sha256", "")).removeprefix("sha256:")
        if not raw or len(digest) != 64:
            raise ValueError(f"Route B claims ancestor {name} lacks path or SHA-256")
        path = Path(raw)
        if not path.is_absolute():
            path = manifest_path.parent / path
        if not path.exists() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f"Route B claims ancestor {name} hash mismatch or missing")


def _verify_route_b_candidate_ancestors(manifest: dict[str, object], manifest_path: Path) -> None:
    """Bind the frozen candidate set to the exact verified-claims generation."""
    inputs = manifest.get("inputs")
    if not isinstance(inputs, dict):
        raise ValueError("Route B candidates manifest lacks verified-claims ancestry")
    for name in ("claims", "claims_manifest"):
        entry = inputs.get(name)
        if not isinstance(entry, dict):
            raise ValueError(f"Route B candidates ancestor {name} is missing")
        raw = str(entry.get("path", "")).strip()
        digest = str(entry.get("sha256", "")).removeprefix("sha256:")
        if not raw or len(digest) != 64:
            raise ValueError(f"Route B candidates ancestor {name} lacks path or SHA-256")
        path = Path(raw)
        if not path.is_absolute():
            path = manifest_path.parent / path
        if not path.exists() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f"Route B candidates ancestor {name} hash mismatch or missing")
    claims_manifest_path = Path(str(inputs["claims_manifest"]["path"]))
    if not claims_manifest_path.is_absolute():
        claims_manifest_path = manifest_path.parent / claims_manifest_path
    try:
        claims_manifest = json.loads(claims_manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError("Route B candidates claims manifest is invalid JSON") from error
    if claims_manifest.get("stage") != "route_b_semantic_transfer_verification":
        raise ValueError("Route B candidates claims ancestor has an unsupported stage")
    if claims_manifest.get("contract") != "route-b-semantic-transfer-verification-1.0.0":
        raise ValueError("Route B candidates claims ancestor has an unsupported contract")
    if str(claims_manifest.get("case_id", "")).strip() != str(manifest.get("case_id", "")).strip():
        raise ValueError("Route B candidates claims ancestor case_id mismatch")
    claims_path = Path(str(inputs["claims"]["path"]))
    if not claims_path.is_absolute():
        claims_path = manifest_path.parent / claims_path
    claims_sha256 = hashlib.sha256(claims_path.read_bytes()).hexdigest()
    declared_claims_sha256 = str(
        claims_manifest.get("artifact_sha256")
        or (claims_manifest.get("outputs") or {}).get("claims_sha256")
        or ""
    ).removeprefix("sha256:")
    if declared_claims_sha256 != claims_sha256:
        raise ValueError("Route B candidates claims ancestor is not bound by its producer manifest")
    declared_claims_output = str((claims_manifest.get("outputs") or {}).get("claims", "")).strip()
    if not declared_claims_output:
        raise ValueError("Route B candidates claims ancestor manifest lacks its output path")
    declared_claims_path = Path(declared_claims_output)
    if not declared_claims_path.is_absolute():
        declared_claims_path = claims_manifest_path.parent / declared_claims_path
    if declared_claims_path.resolve() != claims_path.resolve():
        raise ValueError("Route B candidates claims artifact does not match its producer manifest output")


def _verify_route_a_ancestors(manifest: dict[str, object], manifest_path: Path) -> None:
    """Re-hash every declared Route-A bridge input before support is trusted."""
    inputs = manifest.get("inputs")
    if not isinstance(inputs, dict):
        raise ValueError("Route A bridge manifest lacks input ancestry")
    present = 0
    for name, entry in inputs.items():
        if not isinstance(entry, dict) or not entry.get("present"):
            continue
        present += 1
        raw = str(entry.get("path", "")).strip()
        digest = str(entry.get("sha256", "")).removeprefix("sha256:")
        if not raw or len(digest) != 64 or entry.get("strength") != "strong":
            raise ValueError(f"Route A bridge ancestor {name} lacks a strong path/SHA-256 binding")
        path = Path(raw)
        if not path.is_absolute():
            path = manifest_path.parent / path
        if not path.exists() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f"Route A bridge ancestor {name} hash mismatch or missing")
    if present == 0:
        raise ValueError("Route A bridge manifest has no sealed input ancestry")

OUTPUT_COLUMNS = [
    "candidate_gene_id",
    "route_a_supported",
    "route_b_supported",
    "route_support",
    "route_a_family_ids",
    "route_a_hypothesis_ids",
    "route_a_mapping_tiers",
    "route_a_ambiguity_statuses",
    "route_b_family_ids",
    "route_b_orthogroup_ids",
    "route_b_semantic_unit_ids",
    "route_b_claim_ids",
    "route_b_evidence_scopes",
    "route_c_supported",
    "route_c_record_ids",
    "route_c_source_ids",
    "route_c_species_distances",
    "route_c_evidence_directness",
    "route_c_projection_states",
    "route_c_projection_bases",
    "route_c_family_ids",
    "route_c_orthogroup_ids",
    "support_status",
    "rank_contribution",
    "relevance_evidence_json",
]


def _split(value: object) -> set[str]:
    return {
        token.strip()
        for token in str(value or "").replace(",", ";").split(";")
        if token.strip()
    }


def _join(values: set[str]) -> str:
    return ";".join(sorted(values))


def _evidence_json(records: list[dict[str, object]]) -> str:
    """Serialize member-bound source records deterministically.

    The union is allowed to carry evidence provenance, but it must never
    manufacture a claim by combining fields from different source rows.
    Sorting canonical objects also makes the artifact invariant to producer row
    order while retaining every distinct source record.
    """
    unique = {
        json.dumps(record, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
        for record in records
    }
    return json.dumps(
        [json.loads(value) for value in sorted(unique)],
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )


def _record_members(record: dict[str, object]) -> set[str]:
    return _split(
        record.get("member_interval_genes")
        or record.get("candidate_gene_id")
        or record.get("gene_id")
    )


def _sealed_input_rows(
    manifest: dict[str, object], manifest_path: Path, name: str,
) -> list[dict[str, object]]:
    """Read one explicit, hash-bound producer input; never search directories."""
    inputs = manifest.get("inputs")
    entry = inputs.get(name) if isinstance(inputs, dict) else None
    if not isinstance(entry, dict) or not entry.get("present", True):
        return []
    raw = str(entry.get("path", "")).strip()
    expected = str(entry.get("sha256", "")).removeprefix("sha256:")
    if not raw or len(expected) != 64:
        raise ValueError(f"sealed Route A input {name} lacks path or SHA-256")
    path = Path(raw)
    if not path.is_absolute():
        path = manifest_path.parent / path
    if not path.exists() or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError(f"sealed Route A input {name} hash mismatch or missing")
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return frame.to_dict("records")


def _member_evidence_from_route_a(
    bridge: pd.DataFrame, manifest: dict[str, object], manifest_path: Path,
) -> dict[str, list[dict[str, object]]]:
    """Recover original Route-A hypothesis/panel rows for supported members."""
    hypotheses = _sealed_input_rows(manifest, manifest_path, "route_a_hypotheses")
    panel = _sealed_input_rows(manifest, manifest_path, "route_a_bait_panel")
    panel_by_hypothesis: dict[str, list[dict[str, object]]] = {}
    for bait in panel:
        panel_by_hypothesis.setdefault(str(bait.get("hypothesis_id", "")).strip(), []).append(bait)
    result: dict[str, list[dict[str, object]]] = {}
    for row in bridge[bridge["route"].astype(str).str.casefold().eq("a")].to_dict("records"):
        gene = str(row.get("gene_id", "")).strip()
        if not gene or str(row.get("mapping_tier", "")).strip() == "tier6_not_evaluated":
            continue
        reference_accession = str(row.get("reference_accession", "")).strip()
        # The production bridge schema does not retain hypothesis_id.  Bind
        # the sealed originals to its surviving reference key instead of
        # attaching every row sharing a (possibly absent) hypothesis id.
        matching_panels = [
            bait for bait in panel
            if reference_accession and str(bait.get("bait_id", "")).strip() == reference_accession
        ]
        matching_hypotheses = []
        for hypothesis in hypotheses:
            # Match the producer's actual reference precedence. A symbol may
            # legitimately differ from its verified bait accession, while an
            # explicit proposed accession must not be overridden by that bait.
            stable = str(hypothesis.get("proposed_stable_id", "")).strip()
            baits = panel_by_hypothesis.get(str(hypothesis.get("hypothesis_id", "")).strip(), [])
            references = ({stable} if stable else
                          {str(bait.get("bait_id", "")).strip() for bait in baits} - {""})
            if not references:
                references = {str(hypothesis.get("reference_gene_name", "")).strip()
                              or str(hypothesis.get("reference_genes", "")).strip()} - {""}
            if reference_accession and reference_accession in references:
                matching_hypotheses.append(hypothesis)
        matching_hypotheses = list({
            json.dumps(item, ensure_ascii=False, sort_keys=True): item
            for item in matching_hypotheses
        }.values())
        originals = []
        # Route A is interval-blind: the original hypothesis commonly has no
        # member_interval_genes.  The sealed bridge row is the validated
        # projection binding hypothesis_id/reference to this actual member.
        for hypothesis in matching_hypotheses:
            originals.append({**hypothesis,
                              # These are lossless vocabulary aliases for the
                              # relevance annotator; the original fields stay
                              # present in the evidence object.
                              "family_id": str(hypothesis.get("family_id", "")).strip()
                              or str(hypothesis.get("gene_family", "")).strip(),
                              "reference_gene": str(hypothesis.get("reference_gene", "")).strip()
                              or str(hypothesis.get("reference_gene_name", "")).strip()
                              or str(hypothesis.get("proposed_stable_id", "")).strip(),
                              "phenotype_relation": str(hypothesis.get("phenotype_relation", "")).strip()
                              or str(hypothesis.get("expected_role", "")).strip(),
                              "candidate_gene_id": gene,
                              "evidence_route": "a", "evidence_origin": "route_a_hypothesis"})
        for bait in matching_panels:
            originals.append({**bait, "candidate_gene_id": gene,
                              "evidence_route": "a", "evidence_origin": "route_a_bait_panel"})
        if not originals:
            originals.append({"candidate_gene_id": gene, "evidence_route": "a",
                              "evidence_origin": "route_a_original_evidence_unavailable",
                              "evidence_state": "not_evaluated"})
        result.setdefault(gene, []).extend(originals)
    return result


def build_supported_member_union(
    out_dir_arg: str,
    *,
    route_a_bridge_arg: str | None = None,
    route_b_claims_arg: str | None = None,
    route_b_candidates_arg: str | None = None,
    route_c_members_arg: str | None = None,
    route_a_manifest_arg: str | None = None,
    route_b_claims_manifest_arg: str | None = None,
    route_b_candidates_manifest_arg: str | None = None,
    route_c_manifest_arg: str | None = None,
    case_id: str = "",
) -> dict[str, object]:
    """Create the production supported-member set from frozen Route A/B/C artifacts."""
    if not route_a_bridge_arg and not route_b_claims_arg and not route_b_candidates_arg and not route_c_members_arg:
        raise ValueError("supported member union requires Route A bridge and/or Route B/C verified members")
    if route_b_claims_arg and route_b_candidates_arg:
        raise ValueError("only one Route B input (claims or frozen candidates) may be supplied")

    records: dict[str, dict[str, set[str]]] = {}
    evidence_by_gene: dict[str, list[dict[str, object]]] = {}
    route_a_rows = 0
    route_a_supported_rows = 0
    route_a_not_evaluated_rows = 0
    route_b_rows = 0
    route_b_verified_rows = 0
    route_b_nonverified_rows = 0
    route_b_input_kind: str | None = None
    route_b_candidate_artifact_rows = 0
    route_b_verified_candidate_rows = 0
    route_b_nonverified_candidate_rows = 0
    route_c_rows = 0
    route_c_verified_rows = 0
    route_c_nonverified_rows = 0
    route_c_supported_candidate_rows = 0
    route_c_input_sha256: str | None = None

    def bucket(gene_id: str) -> dict[str, set[str]]:
        if gene_id not in records:
            records[gene_id] = {
                "route_a": set(), "route_b": set(), "route_c": set(),
                "a_family": set(), "a_hypothesis": set(), "a_tier": set(), "a_ambiguity": set(),
                "b_family": set(), "b_orthogroup": set(), "b_unit": set(), "b_claim": set(), "b_scope": set(),
                "c_record": set(), "c_source": set(), "c_species": set(), "c_directness": set(),
                "c_projection_state": set(), "c_projection_basis": set(), "c_family": set(), "c_orthogroup": set(),
            }
        return records[gene_id]

    input_ledger: list[dict[str, object]] = []
    bridge_path = None
    if route_a_bridge_arg:
        bridge_path, producer_manifest, producer_manifest_path = _load_sealed_input(
            route_a_bridge_arg, route_a_manifest_arg, label="Route A bridge", case_id=case_id)
        input_ledger.append({"route": "a", "artifact": str(bridge_path),
                             "artifact_sha256": hashlib.sha256(bridge_path.read_bytes()).hexdigest(),
                             "producer_manifest": producer_manifest_path,
                             "producer_manifest_sha256": hashlib.sha256(runtime_path(producer_manifest_path).read_bytes()).hexdigest(),
                             "producer_contract": producer_manifest.get("contract", producer_manifest.get("command", "")),
                             "case_id": producer_manifest.get("case_id", ""),
                             "qtl_id": producer_manifest.get("qtl_id", ""), "interval_fingerprint": producer_manifest.get("interval_fingerprint", ""), "assembly": producer_manifest.get("assembly", ""),
                             "verification": "sealed_hash_case_checked"})
        if producer_manifest.get("command") not in {"family-candidate-bridge", "family_candidate_bridge"} and producer_manifest.get("stage") not in {"family_candidate_bridge"}:
            raise ValueError("Route A bridge producer manifest has an unsupported contract")
        if producer_manifest.get("contract") != "family-candidate-bridge-1.0.0":
            raise ValueError("Route A bridge producer contract is unsupported")
        _verify_route_a_ancestors(producer_manifest, Path(producer_manifest_path))
    if bridge_path:
        if not bridge_path.exists():
            raise ValueError(f"Route A bridge does not exist: {bridge_path}")
        bridge = pd.read_csv(bridge_path, sep="\t", dtype=str).fillna("")
        required = {"gene_id", "route", "mapping_tier"}
        missing = sorted(required - set(bridge.columns))
        if missing:
            raise ValueError(f"Route A bridge is missing columns: {', '.join(missing)}")
        route_a = bridge[bridge["route"].astype(str).str.casefold().eq("a")].copy()
        for gene, evidence in _member_evidence_from_route_a(
            bridge, producer_manifest, Path(producer_manifest_path)
        ).items():
            evidence_by_gene.setdefault(gene, []).extend(evidence)
        route_a_rows = int(len(route_a))
        for row in route_a.to_dict("records"):
            gene_id = str(row.get("gene_id", "")).strip()
            if not gene_id:
                raise ValueError("Route A bridge contains empty gene_id")
            mapping_tier = str(row.get("mapping_tier", "")).strip()
            provider_status = str(row.get("provider_status", "")).strip()
            review_status = str(row.get("review_status", "")).strip()
            not_evaluated = (
                mapping_tier == "tier6_not_evaluated"
                or review_status == "not_evaluated"
                or provider_status in {"not_run", "missing", "provider_missing"}
            )
            if not_evaluated:
                route_a_not_evaluated_rows += 1
                continue
            route_a_supported_rows += 1
            item = bucket(gene_id)
            item["route_a"].add("a")
            for value, key in (
                (row.get("family_id", ""), "a_family"),
                (row.get("hypothesis_id", ""), "a_hypothesis"),
                (mapping_tier, "a_tier"),
                (row.get("ambiguity_status", ""), "a_ambiguity"),
            ):
                if str(value).strip():
                    item[key].add(str(value).strip())

    claims_path = None
    if route_b_claims_arg:
        claims_path, producer_manifest, producer_manifest_path = _load_sealed_input(
            route_b_claims_arg, route_b_claims_manifest_arg, label="Route B claims", case_id=case_id)
        input_ledger.append({"route": "b", "artifact": str(claims_path),
                             "artifact_sha256": hashlib.sha256(claims_path.read_bytes()).hexdigest(),
                             "producer_manifest": producer_manifest_path,
                             "producer_manifest_sha256": hashlib.sha256(runtime_path(producer_manifest_path).read_bytes()).hexdigest(),
                             "producer_contract": producer_manifest.get("contract", producer_manifest.get("stage", "")),
                             "case_id": producer_manifest.get("case_id", ""),
                             "qtl_id": producer_manifest.get("qtl_id", ""), "interval_fingerprint": producer_manifest.get("interval_fingerprint", ""), "assembly": producer_manifest.get("assembly", ""),
                             "verification": "sealed_hash_case_checked"})
        if producer_manifest.get("stage") != "route_b_semantic_transfer_verification" or producer_manifest.get("contract") != "route-b-semantic-transfer-verification-1.0.0":
            raise ValueError("Route B claims producer manifest has an unsupported contract")
        _verify_route_b_ancestors(producer_manifest, Path(producer_manifest_path))
    if claims_path:
        if not claims_path.exists():
            raise ValueError(f"Route B semantic claims do not exist: {claims_path}")
        claims = pd.read_csv(claims_path, sep="\t", dtype=str).fillna("")
        required = {"claim_id", "verification_status", "member_interval_genes"}
        missing = sorted(required - set(claims.columns))
        if missing:
            raise ValueError(f"Route B semantic claims are missing columns: {', '.join(missing)}")
        route_b_input_kind = "claims"
        route_b_rows = int(len(claims))
        for row in claims.to_dict("records"):
            if str(row.get("verification_status", "")).strip() != "verified":
                route_b_nonverified_rows += 1
                continue
            route_b_verified_rows += 1
            members = _split(row.get("member_interval_genes", ""))
            if not members:
                raise ValueError(f"verified Route B claim {row.get('claim_id', '')} has no member_interval_genes")
            for gene_id in members:
                item = bucket(gene_id)
                item["route_b"].add("b")
                evidence_by_gene.setdefault(gene_id, []).append({
                    **row, "candidate_gene_id": gene_id,
                    "evidence_route": "b", "evidence_origin": "route_b_verified_claim",
                })
                for value, key in (
                    (row.get("family_id", ""), "b_family"),
                    (row.get("orthogroup_id", ""), "b_orthogroup"),
                    (row.get("semantic_unit_id", ""), "b_unit"),
                    (row.get("claim_id", ""), "b_claim"),
                    (row.get("evidence_scope", ""), "b_scope"),
                ):
                    if str(value).strip():
                        item[key].add(str(value).strip())

    candidates_path = None
    if route_b_candidates_arg:
        candidates_path, producer_manifest, producer_manifest_path = _load_sealed_input(
            route_b_candidates_arg, route_b_candidates_manifest_arg, label="Route B candidates", case_id=case_id)
        input_ledger.append({"route": "b", "artifact": str(candidates_path),
                             "artifact_sha256": hashlib.sha256(candidates_path.read_bytes()).hexdigest(),
                             "producer_manifest": producer_manifest_path,
                             "producer_manifest_sha256": hashlib.sha256(runtime_path(producer_manifest_path).read_bytes()).hexdigest(),
                             "producer_contract": producer_manifest.get("contract", producer_manifest.get("stage", "")),
                             "case_id": producer_manifest.get("case_id", ""),
                             "qtl_id": producer_manifest.get("qtl_id", ""), "interval_fingerprint": producer_manifest.get("interval_fingerprint", ""), "assembly": producer_manifest.get("assembly", ""),
                             "verification": "sealed_hash_case_checked"})
    if candidates_path:
        if not candidates_path.exists():
            raise ValueError(f"Route B frozen candidates do not exist: {candidates_path}")
        if producer_manifest.get("contract") != "route-b-verified-candidates-1.0.0":
            raise ValueError("Route B candidates producer contract is unsupported")
        _verify_route_b_candidate_ancestors(producer_manifest, Path(producer_manifest_path))
        candidates = pd.read_csv(candidates_path, sep="\t", dtype=str).fillna("")
        # Frozen candidates normally carry relevance_evidence_json.  Older
        # sealed candidate artifacts are rehydrated from their exact claims
        # ancestor after the existing ancestor/hash checks above.
        candidate_claims: dict[str, list[dict[str, object]]] = {}
        claims_entry = (producer_manifest.get("inputs") or {}).get("claims")
        if isinstance(claims_entry, dict):
            claims_ancestor = Path(str(claims_entry.get("path", "")))
            if not claims_ancestor.is_absolute():
                claims_ancestor = Path(producer_manifest_path).parent / claims_ancestor
            source_claims = pd.read_csv(claims_ancestor, sep="\t", dtype=str).fillna("")
            if {"claim_id", "verification_status", "member_interval_genes"} <= set(source_claims.columns):
                for claim in source_claims.to_dict("records"):
                    if str(claim.get("verification_status", "")).strip() != "verified":
                        continue
                    for member in _record_members(claim):
                        candidate_claims.setdefault(member, []).append({
                            **claim, "candidate_gene_id": member,
                            "evidence_route": "b", "evidence_origin": "route_b_verified_claim",
                        })
        route_b_input_kind = "precomparison_frozen_candidates"
        route_b_candidate_artifact_rows = int(len(candidates))
        for row in candidates.to_dict("records"):
            gene_id = str(row.get("candidate_gene_id", "") or row.get("gene_id", "")).strip()
            if not gene_id:
                raise ValueError("Route B frozen artifact contains empty candidate_gene_id")
            route_val = str(row.get("route", "")).strip().casefold()
            if route_val != "b":
                raise ValueError(f"Route B frozen artifact contains non-Route B row (expected route=b, got {row.get('route', '')})")
            support_status = str(row.get("support_status", "")).strip()
            if support_status != "verified_semantic_support":
                route_b_nonverified_candidate_rows += 1
                continue
            route_b_verified_candidate_rows += 1
            item = bucket(gene_id)
            item["route_b"].add("b")
            saved = str(row.get("relevance_evidence_json", "")).strip()
            candidate_claim_ids = _split(row.get("claim_ids", row.get("claim_id", "")))
            if saved:
                try:
                    saved_records = json.loads(saved)
                except json.JSONDecodeError as error:
                    raise ValueError("Route B candidate relevance evidence is invalid JSON") from error
                if not isinstance(saved_records, list) or any(not isinstance(value, dict) for value in saved_records):
                    raise ValueError("Route B candidate relevance evidence must be an object list")
                for evidence in saved_records:
                    if gene_id not in _record_members(evidence):
                        raise ValueError("Route B candidate relevance evidence is bound to another member")
                    evidence_by_gene.setdefault(gene_id, []).append(evidence)
            for evidence in candidate_claims.get(gene_id, []):
                if not candidate_claim_ids or str(evidence.get("claim_id", "")).strip() in candidate_claim_ids:
                    evidence_by_gene.setdefault(gene_id, []).append(evidence)
            for value, key in (
                (row.get("family_ids", row.get("family_id", "")), "b_family"),
                (row.get("orthogroup_ids", row.get("orthogroup_id", "")), "b_orthogroup"),
                (row.get("semantic_unit_ids", row.get("semantic_unit_id", "")), "b_unit"),
                (row.get("claim_ids", row.get("claim_id", "")), "b_claim"),
                (row.get("evidence_scopes", row.get("evidence_scope", "")), "b_scope"),
            ):
                if str(value).strip():
                    for token in _split(value):
                        item[key].add(token)

    route_c_path = None
    if route_c_members_arg:
        route_c_path, producer_manifest, producer_manifest_path = _load_sealed_input(
            route_c_members_arg, route_c_manifest_arg, label="Route C members", case_id=case_id)
        if producer_manifest.get("stage") != "route_c_direct_members":
            raise ValueError("Route C members require route_c_direct_members producer manifest")
        if str(producer_manifest.get("contract", "")).strip() != "route-c-direct-members-1.1.0":
            raise ValueError("Route C members producer contract is unsupported")
        _verify_route_c_ancestors(producer_manifest, Path(producer_manifest_path), case_id)
        if not str(producer_manifest.get("qtl_id", "")).strip():
            raise ValueError("Route C members manifest lacks qtl_id")
        ancestors = producer_manifest.get("inputs")
        if not isinstance(ancestors, dict) or not all(
            isinstance(ancestors.get(name), dict) and str(ancestors[name].get("sha256", "")).strip()
            for name in ("frozen", "freeze_manifest", "projection", "projection_manifest")
        ):
            raise ValueError("Route C members manifest lacks freeze/projection ancestry")
        input_ledger.append({"route": "c", "artifact": str(route_c_path),
                             "artifact_sha256": hashlib.sha256(route_c_path.read_bytes()).hexdigest(),
                             "producer_manifest": producer_manifest_path,
                             "producer_manifest_sha256": hashlib.sha256(runtime_path(producer_manifest_path).read_bytes()).hexdigest(),
                             "producer_contract": producer_manifest.get("contract", ""),
                             "case_id": producer_manifest.get("case_id", ""),
                             "qtl_id": producer_manifest.get("qtl_id", ""),
                             "verification": "sealed_hash_case_checked"})
    if route_c_path:
        if not route_c_path.exists():
            raise ValueError(f"Route C members do not exist: {route_c_path}")
        route_c_input_sha256 = hashlib.sha256(route_c_path.read_bytes()).hexdigest()
        route_c = pd.read_csv(route_c_path, sep="\t", dtype=str).fillna("")
        required = {"candidate_gene_id", "support_status"}
        missing = sorted(required - set(route_c.columns))
        if missing:
            raise ValueError(f"Route C members are missing columns: {', '.join(missing)}")
        route_c_rows = int(len(route_c))
        # A member artifact may carry one row per source/projection.  Aggregate
        # provenance by candidate while admitting only the explicit support
        # state.  Review-only, outside, unresolved, and not-evaluated rows are
        # intentionally retained only in input accounting, never in the union.
        for row in route_c.to_dict("records"):
            support_status = str(row.get("support_status", "")).strip()
            projection_state = str(row.get("projection_states", row.get("projection_state", ""))).strip()
            if support_status != "verified_direct_support" or projection_state != "inside_current_qtl":
                route_c_nonverified_rows += 1
                continue
            gene_id = str(row.get("candidate_gene_id", "")).strip()
            if not gene_id:
                raise ValueError("verified Route C member contains empty candidate_gene_id")
            route_c_verified_rows += 1
            item = bucket(gene_id)
            item["route_c"].add("c")
            evidence_by_gene.setdefault(gene_id, []).append({
                **row, "candidate_gene_id": gene_id,
                "evidence_route": "c", "evidence_origin": "route_c_verified_member",
            })
            for value, key in (
                (row.get("record_ids", row.get("record_id", "")), "c_record"),
                (row.get("source_ids", row.get("source_id", "")), "c_source"),
                (row.get("species_distances", row.get("species_distance", "")), "c_species"),
                (row.get("evidence_directness", ""), "c_directness"),
                (row.get("projection_states", row.get("projection_state", "")), "c_projection_state"),
                (row.get("projection_bases", row.get("projection_basis", "")), "c_projection_basis"),
                (row.get("family_ids", row.get("family_id", "")), "c_family"),
                (row.get("orthogroup_ids", row.get("orthogroup_id", "")), "c_orthogroup"),
            ):
                if str(value).strip():
                    for token in _split(value):
                        item[key].add(token)
        route_c_supported_candidate_rows = sum(1 for item in records.values() if item["route_c"])

    rows: list[dict[str, object]] = []
    for gene_id in sorted(records):
        item = records[gene_id]
        a = bool(item["route_a"])
        b = bool(item["route_b"])
        c = bool(item["route_c"])
        routes = "+".join(route for route, present in (("a", a), ("b", b), ("c", c)) if present)
        rows.append({
            "candidate_gene_id": gene_id,
            "route_a_supported": a,
            "route_b_supported": b,
            "route_a_family_ids": _join(item["a_family"]),
            "route_a_hypothesis_ids": _join(item["a_hypothesis"]),
            "route_a_mapping_tiers": _join(item["a_tier"]),
            "route_a_ambiguity_statuses": _join(item["a_ambiguity"]),
            "route_b_family_ids": _join(item["b_family"]),
            "route_b_orthogroup_ids": _join(item["b_orthogroup"]),
            "route_b_semantic_unit_ids": _join(item["b_unit"]),
            "route_b_claim_ids": _join(item["b_claim"]),
            "route_b_evidence_scopes": _join(item["b_scope"]),
            "route_c_supported": c,
            "route_c_record_ids": _join(item["c_record"]),
            "route_c_source_ids": _join(item["c_source"]),
            "route_c_species_distances": _join(item["c_species"]),
            "route_c_evidence_directness": _join(item["c_directness"]),
            "route_c_projection_states": _join(item["c_projection_state"]),
            "route_c_projection_bases": _join(item["c_projection_basis"]),
            "route_c_family_ids": _join(item["c_family"]),
            "route_c_orthogroup_ids": _join(item["c_orthogroup"]),
            "support_status": "supported",
            "route_support": routes,
            "rank_contribution": "none",
            "relevance_evidence_json": _evidence_json(evidence_by_gene.get(gene_id, [])),
        })

    frame = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    # A single union cannot silently combine generations from different case,
    # interval or assembly universes. Empty metadata is retained as an
    # explicit producer limitation; conflicting metadata is always fatal.
    for key in ("case_id", "qtl_id", "interval_fingerprint", "assembly"):
        values = {str(item.get(key, "")).strip() for item in input_ledger if str(item.get(key, "")).strip()}
        if len(values) > 1:
            raise ValueError(f"supported-member union inputs disagree on {key}")
    manifest_case_ids = {str(item.get("case_id", "")).strip() for item in input_ledger}
    if "" in manifest_case_ids or len(manifest_case_ids) != 1:
        raise ValueError("supported-member union requires one non-empty producer case_id")
    resolved_case_id = next(iter(manifest_case_ids))
    if case_id and case_id != resolved_case_id:
        raise ValueError("supported-member union producer case_id does not match requested case_id")
    case_id = resolved_case_id
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    output_path = out / "supported_member_union.tsv"
    frame.to_csv(output_path, sep="\t", index=False, lineterminator="\n")
    output_sha256 = hashlib.sha256(output_path.read_bytes()).hexdigest()

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "contract": MANIFEST_CONTRACT,
        "stage": "supported_member_union",
        "case_id": case_id,
        "candidate_count": int(len(frame)),
        "route_a_candidate_count": int(frame["route_a_supported"].sum()) if len(frame) else 0,
        "route_b_candidate_count": int(frame["route_b_supported"].sum()) if len(frame) else 0,
        "route_agreement_candidate_count": int((frame["route_a_supported"] & frame["route_b_supported"]).sum()) if len(frame) else 0,
        "route_a_rows": route_a_rows,
        "route_a_supported_rows": route_a_supported_rows,
        "route_a_not_evaluated_rows": route_a_not_evaluated_rows,
        "route_c_input_status": "supplied" if route_c_path else "not_supplied",
        "route_c_input_sha256": route_c_input_sha256,
        "route_c_rows": route_c_rows,
        "route_c_verified_support_rows": route_c_verified_rows,
        "route_c_nonverified_rows": route_c_nonverified_rows,
        "route_c_supported_candidate_count": route_c_supported_candidate_rows,
        "output_sha256": output_sha256,
        "output": str(output_path),
        "inputs": input_ledger,
        "generation_digest": hashlib.sha256(json.dumps(input_ledger, sort_keys=True, ensure_ascii=False).encode()).hexdigest(),
        "rank_contribution": "none",
        "policy": [
            "Route A tier6/not_evaluated rows do not become supported members",
            "Route A ambiguous supported paralogs remain in the set rather than being collapsed",
            "only final Route B verification_status=verified claims or verified_semantic_support candidates contribute members",
            "all member_interval_genes of a verified group-level Route B claim are retained",
            "Route A/B agreement is provenance only and never a numerical score bonus",
            "the output is a union set for candidate-budget gating, not a ranked Top-N list",
            "only Route C support_status=verified_direct_support contributes members; review-only/outside/unresolved rows do not",
            "Route C is optional and never blocks an A/B-only union",
            "Route C provenance is retained without adding a numerical rank contribution",
        ],
    }
    # Every Route B counter stays present with a zero value even when that lane
    # was not supplied: a key that appears and disappears cannot be schema-checked,
    # and "not supplied" must stay readable rather than becoming an absent field.
    manifest["route_b_input_kind"] = route_b_input_kind or "none"
    manifest["route_b_claim_rows"] = route_b_rows
    manifest["route_b_verified_claim_rows"] = route_b_verified_rows
    manifest["route_b_nonverified_claim_rows"] = route_b_nonverified_rows
    manifest["route_b_candidate_artifact_rows"] = route_b_candidate_artifact_rows
    manifest["route_b_verified_candidate_rows"] = route_b_verified_candidate_rows
    manifest["route_b_nonverified_candidate_rows"] = route_b_nonverified_candidate_rows

    manifest_path = out / "supported_member_union_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "candidate_set": str(output_path),
        "manifest": str(manifest_path),
        "candidate_count": int(len(frame)),
        "route_a_candidate_count": manifest["route_a_candidate_count"],
        "route_b_candidate_count": manifest["route_b_candidate_count"],
        "route_agreement_candidate_count": manifest["route_agreement_candidate_count"],
        "inputs": input_ledger,
    }
