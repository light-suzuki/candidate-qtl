"""Resume existing A/B claim gates from immutable, independently reviewed inputs.

This module does not generate reviews or launch models. The executor receives the
frozen prompt, then supplies its real raw response and execution identity.
"""
from __future__ import annotations

import hashlib
import io
import json
import shutil
import tempfile
from pathlib import Path

import pandas as pd

from .agent_claim_review import AgentClaimReview, B_CONTRACT, B_LEGACY_CONTRACT, CONTRACT, prepare_review_prompt
from .core import runtime_path, sha256 as _sha256


def sha256(path: str | Path) -> str:
    return _sha256(Path(path))


def _read(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _write(path: Path, value: dict) -> None:
    with path.open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n")


def build_no_source_review(prompt: Path, route: str) -> dict:
    """Record missing source evidence, never an executed semantic review."""
    if route not in {"a", "b"}:
        raise ValueError("route must be a or b")
    payload = _read(prompt)
    expected = payload.get("contract")
    if expected not in ({CONTRACT} if route == "a" else {B_CONTRACT, B_LEGACY_CONTRACT}):
        raise ValueError("review prompt contract does not match route")
    sources = payload.get("source_review_records")
    if not isinstance(sources, dict) or sources.get("records") != {}:
        raise ValueError("source records are available; run an independent semantic review")
    hypotheses = payload.get("hypotheses")
    if not isinstance(hypotheses, list) or any(not isinstance(row, dict) for row in hypotheses):
        raise ValueError("review hypotheses must be a list of objects")
    key = "hypothesis_id" if route == "a" else "claim_id"
    ids = [str(row.get(key, "")).strip() for row in hypotheses]
    if any(not value for value in ids) or len(ids) != len(set(ids)):
        raise ValueError("review hypothesis IDs must be unique and non-empty")
    reason = ("No source text is present in the frozen review request; "
              "independent semantic review was not executed. This is missing evidence, "
              "not a biological negative.")
    return {
        "contract": expected,
        "execution_status": "not_executed",
        "decision_basis": "empty frozen source_review_records.records",
        "reviewer_model": "not_invoked",
        "reviews": [{key: value, "decision": "not_evaluated", "reason": reason} for value in ids],
    }


def _files(root: Path) -> dict[str, str]:
    return {str(p.relative_to(root)): sha256(p) for p in sorted(root.rglob("*")) if p.is_file()}


def _verify(root: Path, files: dict) -> None:
    for name, digest in files.items():
        path = (root / name).resolve()
        if not path.is_relative_to(root.resolve()) or not path.is_file() or sha256(path) != digest:
            raise ValueError(f"frozen review session artifact mismatch: {name}")


def _b_import_companions(paths: dict[str, Path], case_id: str) -> dict[Path, str]:
    """Bind an optional import roster to the exact frozen B response and inputs."""
    status = paths["hypotheses"].with_name("semantic_claim_proposal_status.tsv")
    if not status.is_file():
        return {}  # Legacy sessions retain explicit unknown import coverage.
    manifest_path = status.with_name("semantic_claim_proposal_manifest.json")
    source_manifest_path = paths["source_verification"].with_name("semantic_claim_source_verification_manifest.json")
    if not manifest_path.is_file() or not source_manifest_path.is_file():
        raise ValueError("Route B import status requires import and source manifests")
    manifest, source_manifest, bundle = (_read(path) for path in
                                        (manifest_path, source_manifest_path, paths["bundle"]))
    expected_case = case_id or bundle.get("case_id", "")
    for value in (manifest, source_manifest, bundle):
        if value.get("case_id") and expected_case and value["case_id"] != expected_case:
            raise ValueError("Route B import status case_id mismatch")
    outputs = manifest.get("outputs", {})
    if (manifest.get("stage") != "semantic_proposal_import"
            or manifest.get("route") != "interval_first"
            or any(not outputs.get(key) or runtime_path(outputs[key]).resolve() != path
                   for key, path in (("proposals", paths["hypotheses"]), ("status", status)))):
        raise ValueError("Route B import status is not bound to these proposals")
    for key, path in (("proposals_sha256", paths["hypotheses"]),
                      ("bundle_sha256", paths["bundle"]),
                      ("artifact_sha256", paths["source_verification"])):
        if source_manifest.get(key) != sha256(path):
            raise ValueError(f"Route B import status source binding hash mismatch: {key}")
    for key, path in (("bundle_sha256", paths["bundle"]), ("response_sha256", paths["proposal"])):
        if manifest.get(key) and manifest[key] != sha256(path):
            raise ValueError(f"Route B import status input hash mismatch: {key}")

    # Old import manifests did not hash their status output. Reuse the existing
    # deterministic importer to validate it, without retrieval or model calls.
    # Only prompt IDs/unit IDs affect import semantics; prompt text is not read.
    from .route_b_semantic_import import (
        STATUS_COLUMNS, _llm_execution_evidence, import_route_b_semantic_proposals,
    )
    companions = [status, manifest_path, source_manifest_path]
    llm_manifest = None
    if manifest.get("llm_manifest"):
        llm_path = runtime_path(manifest["llm_manifest"]).resolve()
        llm_manifest = _read(llm_path)
        for key in ("bundle_sha256", "prompt_sha256", "response_sha256"):
            digest = llm_manifest.get(key)
            if (not isinstance(digest, str) or not digest.strip()
                    or manifest.get(key) != digest.strip().casefold()
                    or (manifest.get("artifact_hashes") or {}).get(key) != digest.strip().casefold()):
                raise ValueError(f"Route B import status LLM provenance hash mismatch: {key}")
        for key in ("provider", "model"):
            if manifest.get(key, "") != str(llm_manifest.get(key, "")):
                raise ValueError(f"Route B import status LLM provenance mismatch: {key}")
        if manifest.get("llm_execution_evidence") != _llm_execution_evidence(str(llm_path)):
            raise ValueError("Route B import status LLM execution evidence mismatch")
        companions.append(llm_path)
    elif any(manifest.get(key) for key in ("provider", "model", "llm_execution_evidence", "artifact_hashes")):
        raise ValueError("Route B import status LLM provenance requires its original manifest")
    responses = _read(paths["proposal"]).get("candidates")
    if not isinstance(responses, list) or any(not isinstance(row, dict) for row in responses):
        raise ValueError("Route B import status requires the normalized proposal response")
    roster = pd.read_csv(status, sep="\t", dtype=str).fillna("")
    if not {"unit_id", "prompt_id"} <= set(roster):
        raise ValueError("Route B import status requires the frozen unit/prompt roster")
    prompt_units = dict(zip(roster["prompt_id"], roster["unit_id"]))
    if any(row.get("route_b_unit_id") and row["route_b_unit_id"] != prompt_units.get(row.get("prompt_id"))
           for row in responses):
        raise ValueError("Route B import status response unit mismatch")
    original_hashes = {path: sha256(path) for path in [*paths.values(), *companions] if path.is_file()}
    with tempfile.TemporaryDirectory(prefix="cqtl-import-status-") as temporary:
        temporary = Path(temporary)
        prompts = temporary / "prompts.jsonl"
        prompts.write_text("".join(json.dumps({"prompt_id": prompt_id, "route_b_unit_id": unit_id}) + "\n"
                                   for prompt_id, unit_id in zip(roster["prompt_id"], roster["unit_id"])
                                   if prompt_id), encoding="utf-8")
        replay_manifest = None
        if llm_manifest is not None:
            # Preserve the validated original provenance; only the synthetic
            # roster prompt differs from the frozen research prompt.
            replay_manifest = temporary / "llm_manifest.json"
            _write(replay_manifest, {**llm_manifest, "prompt_sha256": sha256(prompts)})
        replay = import_route_b_semantic_proposals(str(paths["bundle"]), str(prompts),
            str(paths["proposal"]), str(temporary / "import"), case_id=case_id,
            llm_manifest_arg=str(replay_manifest) if replay_manifest else None)
        for actual, expected, columns in (
            (status, Path(replay["status"]), STATUS_COLUMNS),
            (paths["hypotheses"], Path(replay["proposals"]), None),
        ):
            left = pd.read_csv(actual, sep="\t", dtype=str).fillna("")
            right = pd.read_csv(expected, sep="\t", dtype=str).fillna("")
            # The original prompt fingerprint remains frozen, but cannot be
            # reconstructed from a roster-only prompt used for this check.
            columns = columns or [column for column in right if column != "prompt_fingerprint"]
            if not set(columns) <= set(left) or not left[columns].equals(right[columns]):
                raise ValueError(f"Route B import status frozen response/roster mismatch: {actual.name}")
    if any(sha256(path) != digest for path, digest in original_hashes.items()):
        raise ValueError("Route B import status input changed during validation")
    return original_hashes


def prepare_session(out_arg: str, *, route: str, inputs: dict[str, str],
                    proposer_id: str, proposer_model: str, case_id: str = "",
                    focal_species: str = "", review_policy: str = "content_v1") -> dict:
    """Snapshot exact gate inputs; emit a self-contained independent-review prompt."""
    required = {"hypotheses", "sources", "proposal", "records"}
    if route == "b":
        required |= {"bundle", "source_verification"}
    elif route != "a":
        raise ValueError("route must be a or b")
    if set(inputs) != required or not proposer_id.strip() or not proposer_model.strip():
        raise ValueError("exact route input set and real proposer identity/model are required")
    paths = {key: runtime_path(value).resolve() for key, value in inputs.items()}
    if any(not path.exists() for path in paths.values()):
        raise ValueError("review input is missing")
    b_import_hashes = _b_import_companions(paths, case_id) if route == "b" else {}
    root = runtime_path(out_arg).resolve()
    root.mkdir(parents=True, exist_ok=False)
    snapshot = root / "inputs"
    snapshot.mkdir()
    copied = {}
    for key, source in paths.items():
        if source.is_dir():
            if key != "records" or any(p.is_symlink() for p in source.rglob("*")):
                raise ValueError("only a regular record-cache directory can be snapshotted")
            target = snapshot / key
            shutil.copytree(source, target)
        else:
            target = snapshot / (key + source.suffix)
            shutil.copyfile(source, target)
        copied[key] = str(target.relative_to(root))
    # Preserve the upstream Route-A policy beside the snapshotted hypotheses.
    # hypothesis_verify inherits strict benchmark mode from this artifact;
    # sessions must not silently downgrade a prepared benchmark to production.
    if route == "a":
        policy_source = paths["hypotheses"].parent / "mechanism_import_manifest.json"
        if policy_source.is_file():
            shutil.copyfile(policy_source, snapshot / "mechanism_import_manifest.json")
    else:
        source_manifest = paths["source_verification"].with_name("semantic_claim_source_verification_manifest.json")
        if source_manifest.is_file():
            shutil.copyfile(source_manifest, snapshot / source_manifest.name)
        import_companions = {
            paths["hypotheses"].with_name("semantic_claim_proposal_status.tsv"),
            paths["hypotheses"].with_name("semantic_claim_proposal_manifest.json"),
            source_manifest,
        }
        for source, digest in b_import_hashes.items():
            key = next((key for key, path in paths.items() if path == source), None)
            target = root / copied[key] if key else snapshot / source.name
            if key is None and source not in import_companions:
                # The LLM manifest's arbitrary basename must not replace inputs.
                target = snapshot / "llm_provenance" / "manifest.json"
                target.parent.mkdir(exist_ok=True)
            if key is None and source != source_manifest:
                shutil.copyfile(source, target)
            if sha256(target) != digest:
                raise ValueError("Route B import status input changed during snapshot")
    prompt = root / "prompt.json"
    from .bibliography import CONTENT_POLICY, INITIAL_POLICY
    if review_policy not in {CONTENT_POLICY, INITIAL_POLICY}:
        raise ValueError("unknown candidate review policy")
    if review_policy == INITIAL_POLICY:
        _write(prompt, {"contract": "initial-candidate-1.0", "route": route,
                       "source_review_records": _read(root / copied["sources"]),
                       "independent_review_required": False})
    else:
        prepare_review_prompt(str(root / copied["hypotheses"]), str(root / copied["sources"]), str(prompt), route=route)
    request = {"schema_version": "1.0", "route": route, "case_id": case_id,
               "focal_species": focal_species, "state": "awaiting_independent_review",
               "proposal_run": {"executor_id": proposer_id, "model": proposer_model},
               "inputs": copied, "files": _files(root), "prompt": "prompt.json",
               "next_action": "Execute prompt.json in a separate context; preserve its raw JSON response and actual executor/model."}
    request["review_policy"] = review_policy
    if review_policy == INITIAL_POLICY:
        request.update(state="ready_for_initial_candidates", next_action="Finish first-stage candidates; independent content review is optional during narrowing.")
    _write(root / "request.json", request)
    return {"request": str(root / "request.json"), "request_sha256": sha256(root / "request.json"),
            "prompt": str(prompt), "state": request["state"]}


def complete_session(session_arg: str, request_sha256: str, *, response_arg: str = "",
                     reviewer_id: str = "", reviewer_model: str = "") -> dict:
    """Validate a genuine review, then invoke existing A or B candidate consumers."""
    root = runtime_path(session_arg).resolve()
    request_path = root / "request.json"
    if sha256(request_path) != request_sha256:
        raise ValueError("review request anchor mismatch")
    request = _read(request_path)
    _verify(root, request["files"])
    from .bibliography import INITIAL_POLICY
    if request.get("review_policy") == INITIAL_POLICY:
        if response_arg or reviewer_id or reviewer_model:
            raise ValueError("initial-candidate session does not consume an independent review; use the explicit content policy for that path")
        return _complete_initial_session(root, request, request_sha256)
    if (not reviewer_id.strip() or not reviewer_model.strip()
            or reviewer_id == request["proposal_run"]["executor_id"]):
        raise ValueError("a separate actual reviewer identity and model are required")
    response = runtime_path(response_arg)
    identity = {"request_sha256": request_sha256, "response_sha256": sha256(response),
                "reviewer_id": reviewer_id, "reviewer_model": reviewer_model}
    terminal = root / "completion.json"
    if terminal.exists():
        prior = _read(terminal)
        if prior["identity"] != identity:
            raise ValueError("completed review cannot be replaced")
        _verify(root, prior["files"])
        return prior
    output = root / "completed"
    # Exclusive directory creation prevents concurrent completion or silent
    # reruns after a crash. A failed attempt and its artifacts remain intact.
    output.mkdir(exist_ok=False)
    try:
        response_payload = _read(response)
        no_source_execution = response_payload.get("execution_status") == "not_executed"
        if not no_source_execution and (reviewer_model == "not_invoked"
                                        or reviewer_id.startswith("no-source-machine-")):
            raise ValueError("machine review identity requires not_executed response")
        if no_source_execution:
            sources = _read(root / request["inputs"]["sources"])
            if (sources.get("records") or reviewer_model != "not_invoked"
                    or not reviewer_id.startswith("no-source-machine-")):
                raise ValueError("not_executed review requires empty source records and machine provenance")
            if any(row.get("decision") != "not_evaluated"
                   for row in response_payload.get("reviews", [])):
                raise ValueError("not_executed review may only record not_evaluated decisions")
        raw = output / "raw_response.json"
        shutil.copyfile(response, raw)
        if sha256(raw) != identity["response_sha256"]:
            raise ValueError("response changed while being copied")
        paths = {key: str(root / relative) for key, relative in request["inputs"].items()}
        packet = {"contract": _read(root / "prompt.json")["contract"],
                  "hypotheses_sha256": sha256(paths["hypotheses"]),
                  "source_review_records_sha256": sha256(paths["sources"]),
                  "proposal_run": {**request["proposal_run"], "raw_response": {
                      "path": paths["proposal"], "sha256": sha256(paths["proposal"])}},
                  "review_run": {"executor_id": reviewer_id, "model": reviewer_model,
                                 "raw_response": {"path": str(raw.resolve()), "sha256": sha256(raw)},
                                 "prompt": {"path": str((root / "prompt.json").resolve()),
                                            "sha256": sha256(root / "prompt.json")}}}
        packet_path = output / "packet.json"
        _write(packet_path, packet)
        checker = AgentClaimReview(str(packet_path), paths["hypotheses"], paths["sources"], route=request["route"])
        review_counts = {"requested": len(checker.hypotheses),
                         "returned": 0 if no_source_execution else len(checker.reviews),
                         "missing": (len(checker.hypotheses) if no_source_execution else
                                     len(set(checker.hypotheses) - set(checker.reviews)))}
        if request["route"] == "a":
            from .llm_pathway import bait_panel_freeze, hypothesis_verify
            verified = hypothesis_verify(paths["hypotheses"], str(output / "verify"),
                focal_species=request["focal_species"], record_cache_arg=paths["records"],
                agent_claim_review_arg=str(packet_path), source_review_records_arg=paths["sources"],
                allow_external_api=False)
            candidates = bait_panel_freeze(paths["hypotheses"], verified["validation"], verified["verified_sequences"],
                str(output / "candidates"), case_id=request["case_id"],
                reviewable_sequences_arg=verified["reviewable_sequences"])
        else:
            from .route_b_semantic_candidates import freeze_route_b_semantic_candidates
            from .route_b_semantic_transfer_verify import verify_route_b_semantic_transfer
            verified = verify_route_b_semantic_transfer(paths["hypotheses"], paths["source_verification"], paths["bundle"],
                str(output / "verify"), agent_claim_review_arg=str(packet_path),
                source_review_records_arg=paths["sources"], record_cache_arg=paths["records"], case_id=request["case_id"])
            candidates = freeze_route_b_semantic_candidates(verified["claims"], str(output / "candidates"), case_id=request["case_id"])
        _verify(root, request["files"])
        result = {"state": "completed" if not review_counts["missing"] else "completed_partial_review",
                  "identity": identity, "review_counts": review_counts,
                  "review_execution_state": ("not_executed_no_source" if no_source_execution else "executed"),
                  "verification": verified, "candidates": candidates, "files": _files(root)}
        _write(terminal, result)
        return result
    except Exception as error:
        _write(output / "failure.json", {"state": "failed", "identity": identity,
                                         "error_type": type(error).__name__, "error": str(error)})
        raise


def _complete_initial_session(root: Path, request: dict, request_sha256: str) -> dict:
    """Reuse identity/transfer consumers without inventing an executed reviewer."""
    identity = {"request_sha256": request_sha256, "review_policy": request["review_policy"]}
    terminal = root / "completion.json"
    if terminal.exists():
        prior = _read(terminal)
        if prior["identity"] != identity:
            raise ValueError("completed initial session cannot be replaced")
        _verify(root, prior["files"])
        return prior
    output = root / "completed"
    output.mkdir(exist_ok=False)
    paths = {key: str(root / relative) for key, relative in request["inputs"].items()}
    if request["route"] == "a":
        from .llm_pathway import bait_panel_freeze, hypothesis_verify
        verified = hypothesis_verify(paths["hypotheses"], str(output / "verify"),
            focal_species=request["focal_species"], record_cache_arg=paths["records"], allow_external_api=False)
        candidates = bait_panel_freeze(paths["hypotheses"], verified["validation"], verified["verified_sequences"],
            str(output / "candidates"), case_id=request["case_id"], reviewable_sequences_arg=verified["reviewable_sequences"])
    else:
        from .route_b_semantic_candidates import freeze_route_b_semantic_candidates
        from .route_b_semantic_transfer_verify import verify_route_b_semantic_transfer
        verified = verify_route_b_semantic_transfer(paths["hypotheses"], paths["source_verification"], paths["bundle"],
            str(output / "verify"), record_cache_arg=paths["records"], case_id=request["case_id"])
        candidates = freeze_route_b_semantic_candidates(verified["claims"], str(output / "candidates"), case_id=request["case_id"])
    _verify(root, request["files"])
    result = {"state": "completed", "identity": identity,
              "review_counts": {"requested": 0, "returned": 0, "missing": 0},
              "review_execution_state": "not_required_initial_policy", "verification": verified,
              "candidates": candidates, "files": _files(root)}
    _write(terminal, result)
    return result


def _handoff_review_counts(completed: dict) -> dict[str, object]:
    """Accept sealed partial results without pretending all claims were reviewed."""
    state = completed.get("state")
    if state not in {"completed", "completed_partial_review"}:
        raise ValueError("handoff requires a completed or completed_partial_review session")
    counts = completed.get("review_counts")
    keys = ("requested", "returned", "missing")
    if not isinstance(counts, dict) or any(
        type(counts.get(key)) is not int or counts[key] < 0 for key in keys
    ):
        raise ValueError("handoff review counts must be non-negative integers")
    if counts["returned"] + counts["missing"] != counts["requested"]:
        raise ValueError("handoff review counts do not cover the requested roster")
    if bool(counts["missing"]) != (state == "completed_partial_review"):
        raise ValueError("handoff review state disagrees with missing review count")
    execution = completed.get("review_execution_state", "executed")
    if execution not in {"executed", "not_executed_no_source", "not_required_initial_policy"}:
        raise ValueError("unknown review execution state")
    if execution == "not_executed_no_source" and counts["returned"]:
        raise ValueError("unexecuted review cannot have returned decisions")
    return {**{key: counts[key] for key in keys}, "review_execution_state": execution}


def _execution_artifact(value: dict, label: str) -> tuple[Path, str]:
    path_value = value.get("path") or value.get("status_path") or value.get("artifact_path") or value.get("members_path")
    digest_value = value.get("sha256") or value.get("status_sha256") or value.get("artifact_sha256") or value.get("members_sha256")
    if not path_value or not digest_value:
        raise ValueError(f"{label} execution state requires a hashed source artifact")
    path = runtime_path(str(path_value)).resolve()
    digest = str(digest_value).removeprefix("sha256:")
    if not path.is_file():
        raise ValueError(f"{label} execution source artifact is missing")
    raw = path.read_bytes()
    actual = hashlib.sha256(raw).hexdigest()
    if len(digest) != 64 or actual != digest:
        raise ValueError(f"{label} execution source artifact hash mismatch")
    if path.suffix.lower() not in {".json", ".jsonl"}:
        return path, digest
    try:
        parsed = json.loads(raw.decode("utf-8-sig"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"{label} execution source artifact is not JSON") from exc
    if isinstance(parsed, dict):
        derived = {
            "state", "path", "sha256", "status_path", "status_sha256",
            "artifact_path", "artifact_sha256", "members_path", "members_sha256",
        }
        for key, claimed in value.items():
            if key not in derived and key in parsed and parsed[key] != claimed:
                raise ValueError(f"{label} execution state is not bound to its source artifact")
    return path, digest


def _validate_route_b_execution(
    value: dict | None,
    case_id: str,
    expected_response_sha256: str = "",
    bundle_path: str | Path | None = None,
    *,
    strict_bundle_accounting: bool = True,
) -> dict:
    if value is None:
        return {"state": "not_recorded", "reason": "route_b_execution_manifest was not supplied"}
    if not isinstance(value, dict):
        raise ValueError("Route B execution context must be a JSON object")
    state = str(value.get("state") or value.get("completion_status") or "").strip()
    if state == "not_recorded":
        return {"state": state, "reason": str(value.get("reason", "not_recorded"))}
    if state not in {"complete", "partial", "not_evaluated"}:
        raise ValueError(f"invalid Route B execution state: {state}")
    if value.get("contract") != "agentic-route-b-execution-1.1.0":
        raise ValueError("Route B execution context has an unsupported contract")
    if str(value.get("case_id", "")).strip() != case_id:
        raise ValueError("Route B execution context case_id mismatch")
    prompt_count = value.get("prompt_count")
    completed = value.get("completed_count")
    not_evaluated = value.get("not_evaluated_count")
    if any(type(item) is not int or item < 0 for item in (prompt_count, completed, not_evaluated)):
        raise ValueError("Route B execution context counts must be non-negative integers")
    if prompt_count == 0:
        if state != "not_evaluated" or completed != 0 or not_evaluated != 0:
            raise ValueError("empty Route B roster must be explicitly not_evaluated")
    elif completed + not_evaluated != prompt_count:
        raise ValueError("Route B execution context counts do not cover the frozen prompt roster")
    source_prompts = str(value.get("source_prompts", "")).strip()
    prompt_digest = str(value.get("prompt_sha256", "")).removeprefix("sha256:")
    if not source_prompts or not Path(source_prompts).is_file() or len(prompt_digest) != 64 or sha256(Path(source_prompts)) != prompt_digest:
        raise ValueError("Route B execution context prompt artifact hash mismatch")
    prompt_rows = [
        json.loads(line) for line in Path(source_prompts).read_text(encoding="utf-8").splitlines() if line.strip()
    ]
    if len(prompt_rows) != prompt_count or any(row.get("case_id") != case_id for row in prompt_rows):
        raise ValueError("Route B execution context prompt roster mismatch")
    if bundle_path is None:
        raise ValueError("Route B execution context requires the frozen B bundle")
    bundle = _read(Path(str(bundle_path)).resolve())
    bundle_units = bundle.get("units")
    observability = bundle.get("compression_observability")
    has_accounting = isinstance(bundle_units, list) and type(bundle.get("semantic_unit_count")) is int \
        and type(bundle.get("execution_prompt_count")) is int
    if strict_bundle_accounting and not has_accounting:
        raise ValueError("Route B frozen bundle lacks semantic unit accounting")
    if has_accounting:
        semantic_count = bundle["semantic_unit_count"]
        execution_prompt_count = bundle["execution_prompt_count"]
        if semantic_count < 0 or execution_prompt_count < 0:
            raise ValueError("Route B frozen bundle semantic counts are invalid")
        if semantic_count != len(bundle_units) or execution_prompt_count != prompt_count:
            raise ValueError("Route B frozen bundle semantic/prompt accounting mismatch")
        if isinstance(observability, dict):
            requested_count = observability.get("requested_unit_count")
            skipped_count = observability.get("skipped_unit_count")
            if (
                type(requested_count) is not int or type(skipped_count) is not int
                or requested_count < 0 or skipped_count < 0
                or requested_count + skipped_count != semantic_count
                or requested_count != prompt_count
            ):
                raise ValueError("Route B frozen bundle semantic/prompt accounting mismatch")
    expected_fingerprints = bundle.get("execution_prompt_fingerprints")
    if not isinstance(expected_fingerprints, list) or len(expected_fingerprints) != prompt_count:
        raise ValueError("Route B frozen bundle lacks execution prompt fingerprints")
    actual_fingerprints = [
        hashlib.sha256(json.dumps(row, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()
        for row in prompt_rows
    ]
    if expected_fingerprints != actual_fingerprints:
        raise ValueError("Route B execution prompts differ from the frozen B bundle")
    response_path = Path(str(value.get("output", ""))).resolve()
    response_digest = str(value.get("response_sha256", "")).removeprefix("sha256:")
    if not response_path.is_file() or len(response_digest) != 64 or sha256(response_path) != response_digest:
        raise ValueError("Route B execution context response artifact hash mismatch")
    if expected_response_sha256 and response_digest != str(expected_response_sha256).removeprefix("sha256:"):
        raise ValueError("Route B execution response is not the B session proposal")
    response_payload = _read(response_path)
    if response_payload.get("case_id") != case_id:
        raise ValueError("Route B execution response case_id mismatch")
    actual_rows = response_payload.get("candidates")
    if not isinstance(actual_rows, list) or len(actual_rows) != prompt_count:
        raise ValueError("Route B execution response candidate roster is incomplete")
    actual_by_id: dict[str, dict] = {}
    for row in actual_rows:
        if not isinstance(row, dict) or not row.get("prompt_id") or str(row["prompt_id"]) in actual_by_id:
            raise ValueError("Route B execution response candidate IDs are invalid")
        actual_by_id[str(row["prompt_id"])] = row
    from .route_b_response_validation import route_b_response_completed, route_b_execution_breakdown

    actual_completed = sum(route_b_response_completed(row) for row in actual_rows)
    if completed != actual_completed or not_evaluated != prompt_count - actual_completed:
        raise ValueError("Route B execution context counts do not match response rows")
    roster = value.get("response_roster")
    if not isinstance(roster, list) or len(roster) != prompt_count:
        raise ValueError("Route B execution context response roster is incomplete")
    roster_by_id: dict[str, dict] = {}
    for row in roster:
        if not isinstance(row, dict) or not row.get("prompt_id") or str(row["prompt_id"]) in roster_by_id:
            raise ValueError("Route B execution context response roster IDs are invalid")
        roster_by_id[str(row["prompt_id"])] = row
    for expected in prompt_rows:
        pid = str(expected["prompt_id"])
        row = roster_by_id.get(pid)
        actual = actual_by_id.get(pid)
        if row is None or actual is None:
            raise ValueError("Route B execution context response roster is missing a prompt")
        unit = str(expected["route_b_unit_id"])
        if str(row.get("route_b_unit_id", "")) != unit or str(actual.get("route_b_unit_id", "")) != unit:
            raise ValueError("Route B execution response unit_id is not bound to the frozen prompt")
        if row.get("parse_status") != actual.get("parse_status"):
            raise ValueError("Route B execution response parse status is not bound to the frozen roster")
        raw_hash = hashlib.sha256(str(actual.get("raw_text", "")).encode()).hexdigest()
        if str(row.get("raw_text_sha256", "")) != raw_hash:
            raise ValueError("Route B execution response raw hash is not bound to the frozen roster")
        if actual.get("parse_status") != "json" and (
            actual.get("model_json") is not None or actual.get("claims")
        ):
            raise ValueError("Route B non-JSON response cannot carry claims")
        if actual.get("parse_status") == "json":
            if not isinstance(actual.get("model_json"), dict):
                raise ValueError("Route B completed response requires model_json")
            try:
                raw_json = json.loads(str(actual.get("raw_text", "")))
            except json.JSONDecodeError as exc:
                raise ValueError("Route B execution response raw text is not valid JSON") from exc
            if actual.get("model_json") != raw_json:
                raise ValueError("Route B execution response model JSON differs from raw text")
    provenance = value.get("partial_response_provenance")
    if not isinstance(provenance, list) or (prompt_count > 0 and not provenance):
        raise ValueError("Route B execution context lacks attempt provenance")
    for row in provenance:
        if (not isinstance(row, dict) or type(row.get("attempt_index")) is not int
                or row["attempt_index"] < 0 or not str(row.get("executor", "")).strip()
                or not str(row.get("model", "")).strip()):
            raise ValueError("Route B execution context attempt provenance is incomplete")
        attempt_path = Path(str(row.get("file", ""))).resolve()
        attempt_digest = str(row.get("sha256", "")).removeprefix("sha256:")
        if not attempt_path.is_file() or len(attempt_digest) != 64 or sha256(attempt_path) != attempt_digest:
            raise ValueError("Route B execution attempt provenance hash mismatch")
    retry_history = value.get("retry_history")
    if not isinstance(retry_history, list):
        raise ValueError("Route B execution context requires retry_history")
    for row in retry_history:
        if not isinstance(row, dict):
            raise ValueError("Route B execution retry history is invalid")
        for path_key, digest_key in (("replaced_response_file", "replaced_response_sha256"), ("replacement_response_file", "replacement_response_sha256")):
            retry_path = Path(str(row.get(path_key, ""))).resolve()
            retry_digest = str(row.get(digest_key, "")).removeprefix("sha256:")
            if not retry_path.is_file() or len(retry_digest) != 64 or sha256(retry_path) != retry_digest:
                raise ValueError("Route B execution retry history hash mismatch")
    expected_state = "not_evaluated" if prompt_count == 0 else ("complete" if completed == prompt_count else "partial")
    if state != expected_state:
        raise ValueError("Route B execution state does not match coverage counts")
    path, digest = _execution_artifact(value, "Route B")
    state_counts = route_b_execution_breakdown(actual_rows, manifest=value)
    return {**value, "state": state, "path": str(path), "sha256": digest, "response_state_counts": state_counts}


def _validate_empty_c_research(value: dict, case_id: str) -> None:
    """Check actual saved search artifacts; empty is an original research attestation."""
    def checked(path: object, digest: object) -> Path:
        target = runtime_path(str(path or "")).resolve()
        if (not target.is_file() or target.stat().st_size == 0
                or sha256(target) != str(digest or "").removeprefix("sha256:")):
            raise ValueError("Route C empty research artifact hash mismatch")
        if target.read_bytes().strip() in {b"", b"{}", b"[]", b"null"}:
            raise ValueError("Route C empty research cannot use empty placeholder evidence")
        return target

    if value.get("research_receipt_path"):
        receipt = _read(checked(value["research_receipt_path"], value.get("research_receipt_sha256")))
        if (receipt.get("contract") != "route-c-native-research-receipt-1.0.0"
                or receipt.get("case_id") != case_id or receipt.get("interval_blind") is not True
                or receipt.get("state") != "completed" or type(receipt.get("finding_count")) is not int
                or receipt.get("finding_count") != 0
                or any(not isinstance(receipt.get(key), str) or not receipt[key].strip()
                       for key in ("worker_id", "model", "no_findings_rationale"))):
            raise ValueError("Route C empty research requires completed interval-blind research attestation")
        queries, sources = receipt.get("queries"), receipt.get("read_sources")
        if not isinstance(queries, list) or not queries or not isinstance(sources, list):
            raise ValueError("Route C empty research requires actual successful query evidence")
        for query in queries:
            if (not isinstance(query, dict) or query.get("status") not in {"retrieved", "no_hit"}
                    or not isinstance(query.get("query"), str) or not query["query"].strip()):
                raise ValueError("Route C empty research cannot use unavailable or failed queries")
            checked(query.get("result_path"), query.get("result_sha256"))
        if any(query["status"] == "retrieved" for query in queries) and not sources:
            raise ValueError("Route C retrieved results require saved read-source evidence")
        for source in sources:
            if not isinstance(source, dict):
                raise ValueError("Route C read-source evidence is invalid")
            checked(source.get("path"), source.get("sha256"))
        return
    if not value.get("search_manifest_path"):
        raise ValueError("Route C completed_empty requires sealed search-generation lineage or native research receipt")
    manifest = _read(checked(value["search_manifest_path"], value.get("search_manifest_sha256")))
    if (manifest.get("command") != "route-c-direct-search" or manifest.get("case_id") != case_id
            or manifest.get("interval_blind") is not True):
        raise ValueError("Route C empty search manifest identity mismatch")
    digests = []
    tables = {}
    for key, hash_key in (("candidates", "candidates_sha256"), ("status", "status_sha256"),
                          ("query_provenance", "query_provenance_sha256")):
        target = checked(manifest.get("outputs", {}).get(key), manifest.get(hash_key))
        digests.append(sha256(target))
        tables[key] = pd.read_csv(target, sep="\t", dtype=str).fillna("")
    from .route_c_direct_search_cassette import _generation_sha256
    if manifest.get("search_generation_sha256") != _generation_sha256(*digests):
        raise ValueError("Route C search-generation hash mismatch")
    statuses = tables["status"]
    if (len(statuses) < 1 or len(statuses) != manifest.get("query_count")
            or not statuses["query_status"].isin({"retrieved", "no_hit"}).all()):
        raise ValueError("Route C empty research requires successful queries, not failed search")
    if not str(value.get("reason", "")).strip():
        raise ValueError("Route C empty search requires the original no-findings rationale")
    from .literature_cassette import _load_cassette
    provenance = tables["query_provenance"]
    if len(provenance) != len(statuses):
        raise ValueError("Route C query provenance roster mismatch")
    for row in statuses.to_dict("records"):
        if not str(row.get("query", "")).strip() or not row.get("cassette_path"):
            raise ValueError("Route C empty search requires actual cassette-backed query evidence")
        if Path(str(row["cassette_path"])).name != row["cassette_path"]:
            raise ValueError("Route C cassette path must be a local filename")
        cassette_path = runtime_path(str(manifest["outputs"]["cassette_directory"])) / str(row["cassette_path"])
        cassette = _read(cassette_path)
        request = {key: cassette.get(key) for key in
                   ("request_schema_version", "provider", "endpoint", "request", "parser_schema_version")}
        loaded = _load_cassette(cassette_path, request)
        matched = provenance[provenance["request_key_sha256"].eq(loaded.request_key_sha256)]
        if (loaded.query != row["query"] or loaded.retrieval_status != row["query_status"]
                or loaded.raw_response_sha256 != row["raw_response_sha256"]
                or loaded.request_key_sha256 != row["request_key_sha256"] or len(matched) != 1
                or matched.iloc[0]["raw_response_sha256"] != loaded.raw_response_sha256):
            raise ValueError("Route C query/cassette provenance mismatch")


def _validate_route_c_execution(value: dict | None, case_id: str, members_arg: str | None, producer_manifest_arg: str | None = None) -> dict:
    allowed = {
        "not_evaluated_no_input", "provider_failure", "not_evaluated", "not_attempted",
        "completed_empty", "completed_members",
    }
    if value is None:
        if not members_arg:
            return {"state": "not_evaluated_no_input"}
        member_path = runtime_path(members_arg).resolve()
        count = len(pd.read_csv(member_path, sep="\t", dtype=str))
        if count == 0:
            raise ValueError("empty Direct-C members require an explicit completed_empty status artifact")
        if not producer_manifest_arg:
            return {
                "state": "not_evaluated_no_input",
                "reason": "completed Direct-C members require a producer manifest",
            }
        producer_path = runtime_path(producer_manifest_arg).resolve()
        member_digest = sha256(member_path)
        return {
            "state": "completed_empty" if count == 0 else "completed_members",
            "case_id": case_id,
            "artifact_path": str(member_path),
            "artifact_sha256": member_digest,
            "members_path": str(member_path),
            "members_sha256": member_digest,
            "manifest_path": str(producer_path),
            "manifest_sha256": sha256(producer_path),
        }
    if not isinstance(value, dict) or str(value.get("state", "")).strip() not in allowed:
        raise ValueError("Route C execution context has an invalid state")
    state = str(value["state"]).strip()
    if state == "completed_empty":
        _validate_empty_c_research(value, case_id)
    declared_case = str(value.get("case_id", "")).strip()
    if state != "not_evaluated_no_input" and not declared_case and not (
        members_arg and state in {"completed_empty", "completed_members"} and not value.get("status_path")
    ):
        raise ValueError("Route C execution status requires case_id")
    if declared_case and declared_case != case_id:
        raise ValueError("Route C execution context case_id mismatch")
    if state in {"provider_failure", "not_evaluated", "not_attempted"} and not str(value.get("reason", "")).strip():
        raise ValueError(f"Route C {state} requires a reason")
    if members_arg and state not in {"completed_empty", "completed_members"}:
        raise ValueError("Route C member rows cannot be paired with a non-completed execution state")
    if state == "completed_empty" and not value.get("status_path"):
        raise ValueError("empty Direct-C execution requires a sealed status artifact")
    if state == "completed_empty":
        status = _read(runtime_path(str(value["status_path"])))
        if status.get("state") != state or status.get("case_id") != case_id:
            raise ValueError("Route C empty status identity/state mismatch")
        for key in ("research_receipt_path", "research_receipt_sha256", "search_manifest_path", "search_manifest_sha256"):
            if key in value and status.get(key) != value[key]:
                raise ValueError("Route C empty research lineage is not bound to status artifact")
    if members_arg and state in {"completed_empty", "completed_members"}:
        member_file = runtime_path(members_arg).resolve()
        if not member_file.is_file():
            raise ValueError("Route C member artifact is missing")
        member_bytes = member_file.read_bytes()
        member_count = len(pd.read_csv(io.BytesIO(member_bytes), sep="\t", dtype=str))
        expected_state = "completed_empty" if member_count == 0 else "completed_members"
        if state != expected_state:
            raise ValueError("Route C execution state does not match member row count")
    if state in {"completed_empty", "completed_members"} and not members_arg:
        raise ValueError("Route C completed state requires a member artifact")
    if state in {"completed_empty", "completed_members"} and not producer_manifest_arg:
        raise ValueError("Route C completed state requires a producer manifest")
    if state == "not_evaluated_no_input" and not members_arg and not (value.get("status_path") or value.get("artifact_path")):
        if set(value) == {"state"}:
            return {"state": state}
        raise ValueError("explicit Route C not_evaluated_no_input requires a hashed status artifact")
    normalized = dict(value)
    if state != "not_evaluated_no_input" or value.get("status_path") or value.get("artifact_path"):
        artifact_path, artifact_digest = _execution_artifact(value, "Route C")
        if value.get("status_path") or value.get("artifact_path"):
            normalized["status_path"] = str(artifact_path)
            normalized["status_sha256"] = artifact_digest
    if not declared_case and members_arg and state in {"completed_empty", "completed_members"}:
        normalized["case_id"] = case_id
    if members_arg and state in {"completed_empty", "completed_members"}:
        member_path = runtime_path(members_arg).resolve()
        normalized["members_path"] = str(member_path)
        normalized["members_sha256"] = sha256(member_path)
        if producer_manifest_arg:
            producer_path = runtime_path(producer_manifest_arg).resolve()
            normalized["manifest_path"] = str(producer_path)
            normalized["manifest_sha256"] = sha256(producer_path)
    return normalized


def handoff_sessions(a_session_arg: str, b_session_arg: str, out_arg: str, *,
                     genes_arg: str, proteins_arg: str, atlas_arg: str,
                     family_index_arg: str | None = None,
                     resources_arg: str | None = None,
                     route_c_members_arg: str | None = None,
                     route_c_manifest_arg: str | None = None,
                     route_b_execution: dict | None = None,
                     route_c_execution: dict | None = None,
                     allow_unrecorded_route_b: bool = False,
                     allow_provisional_wet: bool = False,
                     reference_function_arg: str | None = None, threads: int = 2,
                     max_adaptive_rounds: int = 1,
                     identity_assessment_arg: str | None = None) -> dict:
    """Consume sealed complete/partial reviews; missing claims never gain support."""
    from .bibliography import INITIAL_POLICY

    sessions = {}
    review_routes = {}
    for route, arg in (("a", a_session_arg), ("b", b_session_arg)):
        root = runtime_path(arg).resolve()
        request, completed = _read(root / "request.json"), _read(root / "completion.json")
        if request["route"] != route:
            raise ValueError("handoff requires the correct review session route")
        counts = _handoff_review_counts(completed)
        if sha256(root / "request.json") != completed["identity"]["request_sha256"]:
            raise ValueError("completed session request anchor mismatch")
        _verify(root, completed["files"])
        sessions[route] = (root, request, completed)
        review_routes[route] = {"state": completed["state"], **counts}
    claim_review = {
        "state": ("not_required_initial_policy" if all(row["review_execution_state"] == "not_required_initial_policy" for row in review_routes.values())
                  else "partial" if any(row["missing"] for row in review_routes.values()) else "complete"),
        "routes": review_routes,
    }
    aroot, areq, a = sessions["a"]
    broot, breq, b = sessions["b"]
    if not areq["case_id"] or areq["case_id"] != breq["case_id"]:
        raise ValueError("A/B review session case IDs must agree")
    case_id = areq["case_id"]
    # Resolve the same reusable reference-genome index used by guided-profile.
    # This is post-freeze orchestration; Route A's original research is untouched.
    from .family_index_resolution import _reusable_index_validation, resolve_family_index
    family_index_path, family_index_resolution = resolve_family_index(family_index_arg, resources_arg)
    family_index_manifest_path = None
    if family_index_path:
        family_index_path = str(runtime_path(family_index_path).resolve())
        validation = _reusable_index_validation(Path(family_index_path))
        if validation["status"] != "valid":
            raise ValueError("handoff requires a reference_genome family index: "
                             + str(validation.get("reason")))
        family_index_resolution = {**family_index_resolution, "path": family_index_path,
                                   "validation": validation}
        if validation.get("manifest"):
            family_index_manifest_path = str(runtime_path(validation["manifest"]).resolve())
    b_proposal = broot / str(breq["inputs"]["proposal"]) if "proposal" in breq.get("inputs", {}) else None
    b_bundle = broot / str(breq["inputs"]["bundle"]) if "bundle" in breq.get("inputs", {}) else None
    unrecorded_route_b = route_b_execution is None or route_b_execution.get("state") == "not_recorded"
    route_b_not_run = route_b_execution is not None and route_b_execution.get("state") == "not_run_optional"
    if not allow_unrecorded_route_b and unrecorded_route_b:
        raise ValueError("Route B execution manifest is required for canonical handoff")
    if route_b_not_run:
        # A/C-first (2026-11-08): optional B deliberately not run. It contributes
        # no candidates and no coverage; it is never reported as executed.
        reason = str(route_b_execution.get("reason", "")).strip()
        if not reason:
            raise ValueError("Route B not-run state requires a reason")
        if b_proposal is not None:
            from .route_b_response_validation import validate_route_b_not_run
            validate_route_b_not_run(_read(b_proposal).get("candidates"))
        validated_route_b = {"state": "not_run_optional", "reason": reason,
                             "prompt_count": route_b_execution.get("prompt_count"),
                             "completed_count": 0,
                             "not_evaluated_count": route_b_execution.get("not_evaluated_count")}
    elif unrecorded_route_b:
        validated_route_b = {"state": "not_recorded", "reason": "legacy session has no Route B execution manifest"}
    else:
        if b_proposal is None or b_bundle is None:
            raise ValueError("Route B execution context requires proposal and bundle session inputs")
        validated_route_b = _validate_route_b_execution(
            route_b_execution, case_id, sha256(b_proposal), b_bundle,
            strict_bundle_accounting=not allow_unrecorded_route_b,
        )
    provisional_screen_path = None
    provisional_screen_sha256 = None
    if allow_provisional_wet:
        if validated_route_b.get("state") in {"not_recorded", "not_run_optional"}:
            raise ValueError("provisional Wet requires validated Route B execution")
        provisional_screen_path = validated_route_b.get("screen_coverage_ledger")
        provisional_screen_sha256 = validated_route_b.get("screen_coverage_ledger_sha256")
        if not provisional_screen_path or not provisional_screen_sha256:
            raise ValueError("provisional Wet requires a hashed screen coverage ledger")
        if sha256(runtime_path(provisional_screen_path)) != str(provisional_screen_sha256).removeprefix("sha256:"):
            raise ValueError("provisional Wet screen coverage ledger hash mismatch")
    execution_context = {
        "claim_review": claim_review,
        "route_b_execution": validated_route_b,
        "route_c_execution": _validate_route_c_execution(route_c_execution, case_id, route_c_members_arg, route_c_manifest_arg),
        "family_index_resolution": family_index_resolution,
    }
    if bool(route_c_members_arg) != bool(route_c_manifest_arg):
        raise ValueError("Route C handoff requires both members TSV and manifest")
    inputs = {"genes": genes_arg, "proteins": proteins_arg, "atlas": atlas_arg}
    if family_index_path:
        inputs["family_index"] = family_index_path
    if family_index_manifest_path:
        inputs["family_index_manifest"] = family_index_manifest_path
    if reference_function_arg:
        inputs["reference_function"] = reference_function_arg
    if route_c_members_arg:
        inputs["route_c_members"] = route_c_members_arg
    if route_c_manifest_arg:
        inputs["route_c_manifest"] = route_c_manifest_arg
    if provisional_screen_path:
        inputs["screen_coverage_ledger"] = provisional_screen_path
    source_hashes = {k: {"path": str(runtime_path(v).resolve()), "sha256": sha256(runtime_path(v))}
                     for k, v in inputs.items()}
    boundary_bundle = b_bundle if b_bundle is not None and breq.get("review_policy") == INITIAL_POLICY else None
    if boundary_bundle is not None:
        source_hashes["boundary_route_b_bundle"] = {
            "path": str(boundary_bundle.resolve()), "sha256": breq["files"][breq["inputs"]["bundle"]],
        }
    output = runtime_path(out_arg).resolve()
    output.mkdir(parents=True, exist_ok=False)
    identity_descriptor = None
    if identity_assessment_arg:
        identity_path = runtime_path(identity_assessment_arg).resolve()
        identity_descriptor = {"path": str(identity_path), "sha256": sha256(identity_path)}
        source_hashes["identity_assessment"] = identity_descriptor
    if route_b_execution and route_b_execution.get("state") not in {"not_recorded", "not_run_optional"}:
        sealed_b = dict(execution_context["route_b_execution"])
        sealed_payload = _read(Path(str(sealed_b["path"])).resolve())
        prompt_target = output / "route_b_prompts.jsonl"
        response_target = output / "route_b_response.json"
        shutil.copyfile(Path(str(sealed_payload["source_prompts"])).resolve(), prompt_target)
        shutil.copyfile(Path(str(sealed_payload["output"])).resolve(), response_target)
        sealed_payload["source_prompts"] = str(prompt_target)
        sealed_payload["output"] = str(response_target)
        attempt_targets: dict[str, str] = {}

        def seal_attempt(value: str, prefix: str) -> str:
            source = Path(value).resolve()
            if str(source) not in attempt_targets:
                target = output / f"route_b_{prefix}_{len(attempt_targets)}{source.suffix}"
                shutil.copyfile(source, target)
                attempt_targets[str(source)] = str(target)
            return attempt_targets[str(source)]

        provenance = []
        for index, row in enumerate(sealed_payload.get("partial_response_provenance", [])):
            row = dict(row)
            row["file"] = seal_attempt(str(row["file"]), "attempt")
            row["sha256"] = sha256(Path(row["file"]))
            provenance.append(row)
        retry_history = []
        for row in sealed_payload.get("retry_history", []):
            row = dict(row)
            for path_key, digest_key in (("replaced_response_file", "replaced_response_sha256"), ("replacement_response_file", "replacement_response_sha256")):
                row[path_key] = seal_attempt(str(row[path_key]), "retry")
                row[digest_key] = sha256(Path(row[path_key]))
            retry_history.append(row)
        sealed_payload["partial_response_provenance"] = provenance
        sealed_payload["retry_history"] = retry_history
        target = output / "route_b_execution_manifest.json"
        _write(target, sealed_payload)
        execution_context["route_b_execution"] = {
            **sealed_b,
            "path": str(target),
            "sha256": sha256(target),
            "source_prompts": str(prompt_target),
            "output": str(response_target),
            "partial_response_provenance": provenance,
            "retry_history": retry_history,
        }
    c_context = dict(execution_context.get("route_c_execution", {}))
    if c_context.get("state") == "completed_empty":
        _validate_empty_c_research(c_context, case_id)
    if c_context.get("research_receipt_path"):
        receipt_source = runtime_path(str(c_context["research_receipt_path"])).resolve()
        receipt = _read(receipt_source)
        shutil.copyfile(receipt_source, output / "route_c_research_receipt_original.json")
        if sha256(output / "route_c_research_receipt_original.json") != str(c_context["research_receipt_sha256"]).removeprefix("sha256:"):
            raise ValueError("Route C research receipt changed while sealing")
        c_context["original_research_receipt_sha256"] = sha256(receipt_source)
        for index, query in enumerate(receipt["queries"]):
            target = output / f"route_c_query_result_{index}.txt"
            shutil.copyfile(runtime_path(query["result_path"]), target)
            if sha256(target) != str(query["result_sha256"]).removeprefix("sha256:"):
                raise ValueError("Route C query result changed while sealing")
            query["result_path"], query["result_sha256"] = str(target), sha256(target)
        for index, source in enumerate(receipt["read_sources"]):
            target = output / f"route_c_read_source_{index}.txt"
            shutil.copyfile(runtime_path(source["path"]), target)
            if sha256(target) != str(source["sha256"]).removeprefix("sha256:"):
                raise ValueError("Route C read source changed while sealing")
            source["path"], source["sha256"] = str(target), sha256(target)
        target = output / "route_c_research_receipt.json"
        _write(target, receipt)
        c_context["research_receipt_path"], c_context["research_receipt_sha256"] = str(target), sha256(target)
    if c_context.get("search_manifest_path"):
        search_source = runtime_path(c_context["search_manifest_path"])
        search_manifest = _read(search_source)
        search_root = output / "route_c_search"
        search_root.mkdir()
        shutil.copyfile(search_source, search_root / "original_manifest.json")
        c_context["original_search_manifest_sha256"] = sha256(search_root / "original_manifest.json")
        if c_context["original_search_manifest_sha256"] != str(c_context["search_manifest_sha256"]).removeprefix("sha256:"):
            raise ValueError("Route C search manifest changed while sealing")
        for key, digest_key in (("candidates", "candidates_sha256"), ("status", "status_sha256"),
                                ("query_provenance", "query_provenance_sha256")):
            target = search_root / f"{key}.tsv"
            shutil.copyfile(runtime_path(search_manifest["outputs"][key]), target)
            if sha256(target) != search_manifest[digest_key]:
                raise ValueError("Route C search evidence changed while sealing")
            search_manifest["outputs"][key] = str(target)
        cassette_root = search_root / "cassettes"
        cassette_root.mkdir()
        statuses = pd.read_csv(search_root / "status.tsv", sep="\t", dtype=str).fillna("")
        for name in set(statuses["cassette_path"]):
            original = runtime_path(search_manifest["outputs"]["cassette_directory"]) / name
            shutil.copyfile(original, cassette_root / name)
        search_manifest["outputs"]["cassette_directory"] = str(cassette_root)
        target = search_root / "manifest.json"
        search_manifest["outputs"]["manifest"] = str(target)
        _write(target, search_manifest)
        c_context["search_manifest_path"], c_context["search_manifest_sha256"] = str(target), sha256(target)
        _validate_empty_c_research(c_context, case_id)
    if c_context.get("status_path"):
        source = runtime_path(str(c_context["status_path"])).resolve()
        target = output / "route_c_status.json"
        if c_context.get("research_receipt_path") or c_context.get("search_manifest_path"):
            shutil.copyfile(source, output / "route_c_status_original.json")
            c_context["original_status_sha256"] = sha256(source)
            status_payload = _read(source)
            for key in ("research_receipt_path", "research_receipt_sha256", "search_manifest_path", "search_manifest_sha256"):
                if key in c_context:
                    status_payload[key] = c_context[key]
            _write(target, status_payload)
        else:
            shutil.copyfile(source, target)
        c_context["status_path"] = str(target)
        c_context["status_sha256"] = sha256(target)
    if c_context.get("members_path"):
        source = Path(str(c_context["members_path"])).resolve()
        target = output / "route_c_members.tsv"
        shutil.copyfile(source, target)
        c_context["members_path"] = str(target)
        c_context["members_sha256"] = sha256(target)
    if c_context.get("manifest_path"):
        source = Path(str(c_context["manifest_path"])).resolve()
        target = output / "route_c_manifest.json"
        manifest_payload = _read(source)
        copied_members = c_context.get("members_path")
        if copied_members:
            manifest_payload["output"] = str(Path(str(copied_members)).resolve())
        ancestors = manifest_payload.get("inputs")
        binding_targets: dict[str, str] = {}

        def seal_claim_binding(descriptor: object, base: Path) -> str | None:
            if not isinstance(descriptor, dict) or not descriptor.get("artifact"):
                return None
            artifact_value = str(descriptor["artifact"])
            artifact_source = Path(artifact_value)
            if not artifact_source.is_absolute():
                artifact_source = base / artifact_source
            artifact_source = artifact_source.resolve()
            expected_sha = str(descriptor.get("sha256", "")).removeprefix("sha256:")
            key = str(artifact_source)
            # Each descriptor is independently bound, including repeated
            # references to a sidecar already copied for another ancestor.
            if (
                not artifact_source.is_file()
                or len(expected_sha) != 64
                or sha256(artifact_source) != expected_sha
            ):
                raise ValueError("Route C claim-binding artifact hash mismatch before sealing")
            if key in binding_targets:
                target = Path(binding_targets[key])
            else:
                target = output / f"route_c_claim_binding_{len(binding_targets)}{artifact_source.suffix}"
                shutil.copyfile(artifact_source, target)
                binding_targets[key] = str(target.resolve())
            if sha256(target) != expected_sha:
                raise ValueError("Route C claim-binding artifact changed while sealing")
            descriptor["artifact"] = binding_targets[key]
            descriptor["sha256"] = expected_sha
            return binding_targets[key]

        top_binding = manifest_payload.get("claim_binding")
        if isinstance(top_binding, dict):
            seal_claim_binding(top_binding, source.parent)
        if isinstance(ancestors, dict):
            for index, (name, spec) in enumerate(sorted(ancestors.items())):
                if not isinstance(spec, dict) or not spec.get("path"):
                    continue
                ancestor_source = Path(str(spec["path"]))
                if not ancestor_source.is_absolute():
                    ancestor_source = source.parent / ancestor_source
                ancestor_source = ancestor_source.resolve()
                expected_ancestor_sha = str(spec.get("sha256", "")).removeprefix("sha256:")
                if (
                    not ancestor_source.is_file()
                    or len(expected_ancestor_sha) != 64
                    or sha256(ancestor_source) != expected_ancestor_sha
                ):
                    raise ValueError("Route C producer ancestor hash mismatch before sealing")
                ancestor_target = output / (
                    f"route_c_ancestor_{index}_{name}{ancestor_source.suffix}"
                )
                freeze_payload = _read(ancestor_source) if name == "freeze_manifest" else None
                rebound = bool(freeze_payload is not None and seal_claim_binding(
                    freeze_payload.get("claim_binding"), ancestor_source.parent
                ))
                if sha256(ancestor_source) != expected_ancestor_sha:
                    raise ValueError("Route C producer ancestor changed while sealing")
                if rebound:
                    # Create the relocated manifest only once. The original
                    # freeze stays immutable and the new bytes get a new hash.
                    _write(ancestor_target, freeze_payload)
                    spec["sha256"] = sha256(ancestor_target)
                else:
                    shutil.copyfile(ancestor_source, ancestor_target)
                    if sha256(ancestor_target) != expected_ancestor_sha:
                        raise ValueError("Route C producer ancestor changed while sealing")
                spec["path"] = str(ancestor_target.resolve())
        _write(target, manifest_payload)
        c_context["manifest_path"] = str(target)
        c_context["manifest_sha256"] = sha256(target)
    if c_context:
        execution_context["route_c_execution"] = c_context
    route_c_members_consumed = (
        str(c_context.get("members_path")) if c_context.get("members_path") else route_c_members_arg
    )
    route_c_manifest_consumed = (
        str(c_context.get("manifest_path")) if c_context.get("manifest_path") else route_c_manifest_arg
    )
    a_completion_target = output / "a_completion.json"
    b_completion_target = output / "b_completion.json"
    shutil.copyfile(aroot / "completion.json", a_completion_target)
    shutil.copyfile(broot / "completion.json", b_completion_target)
    _write(output / "inputs.json", {"case_id": case_id, "inputs": source_hashes,
                                    "claim_review": claim_review,
                                    "execution_context": execution_context,
                                    "post_freeze_identity_assessment": identity_descriptor,
                                    "a_completion": {"path": str(a_completion_target), "sha256": sha256(a_completion_target)},
                                    "b_completion": {"path": str(b_completion_target), "sha256": sha256(b_completion_target)}})
    try:
        from .blast_panel import freeze_blast_panel, panel_from_verified_bait_panel
        from .blast_evidence import run_panel_blastp
        from .canonical_family_recovery import build_canonical_family_recovery
        from .dual_route import phenotype_first_status
        from .family_candidate_bridge import family_candidate_bridge
        from .family_index import route_a_interval_family_match
        from .family_to_copy import family_to_copy
        from .supported_member_union import build_supported_member_union
        from .wet_candidate_handoff import wet_candidate_handoff
        av, ap = a["verification"], a["candidates"]
        panel = panel_from_verified_bait_panel(ap["projection_panel"], av["reviewable_sequences"])
        panel_path = output / "projection_panel.tsv"
        panel.to_csv(panel_path, sep="\t", index=False, lineterminator="\n")
        freeze = freeze_blast_panel(str(panel_path), str(output / "blast_freeze"), case_id=case_id)
        blast = run_panel_blastp(freeze["manifest"], proteins_arg, str(output / "blast"), num_threads=threads, runner="native")
        interval = route_a_interval_family_match(ap["panel"], str(aroot / areq["inputs"]["hypotheses"]),
            str(output / "a_interval"), family_index_arg=family_index_path,
            candidates_arg=genes_arg, case_id=case_id)
        status = phenotype_first_status(str(output / "a_status"),
            hypotheses_arg=str(aroot / areq["inputs"]["hypotheses"]), validation_arg=av["validation"],
            panel_arg=ap["panel"], case_id=case_id)
        bridge = family_candidate_bridge(str(output / "a_bridge"), route_a_hypotheses_arg=status["hypotheses"],
            route_a_family_status_arg=interval["status"],
            route_a_bait_panel_arg=ap["panel"], route_a_projection_panel_arg=ap["projection_panel"],
            route_a_blast_arg=blast["evidence"], route_a_blast_eligible_arg=blast["eligible_hits"],
            genes_arg=genes_arg, case_id=case_id)
        bc = b["candidates"]
        if route_b_not_run:
            for key in ("candidates", "review_candidates"):
                if len(pd.read_csv(bc[key], sep="\t", dtype=str)):
                    raise ValueError("Route B not-run state cannot contain candidate rows")
        b_candidates = None if route_b_not_run else bc["candidates"]
        union = build_supported_member_union(str(output / "union"), route_a_bridge_arg=bridge["bridge_tsv"],
            route_b_candidates_arg=b_candidates, route_c_members_arg=route_c_members_consumed,
            route_c_manifest_arg=route_c_manifest_consumed, case_id=case_id)
        family = build_canonical_family_recovery(union["candidate_set"], str(output / "family"),
                                                 legacy_bridge_arg=bridge["bridge_tsv"], case_id=case_id)
        copy_root = output / "copy"
        if family["family_count"]:
            family_to_copy(family["family_recovery"], str(copy_root), atlas_arg=atlas_arg,
                           bridge_arg=family["family_bridge"], case_id=case_id)
        copy_paths = {"family_copy_arg": "family_copy_discrimination.tsv",
                      "family_copy_shortlist_arg": "family_copy_shortlist.tsv"}
        wet = wet_candidate_handoff(str(output / "wet"), route_a_bridge_arg=bridge["bridge_tsv"],
            route_b_candidates_arg=b_candidates, route_c_members_arg=route_c_members_consumed,
            route_c_manifest_arg=route_c_manifest_consumed, interval_atlas_arg=atlas_arg,
            reference_function_arg=reference_function_arg, family_recovery_arg=family["family_recovery"],
            convergence_arg=family["comparison"],
            review_candidates_arg=[bridge["review_candidates"]] + ([] if route_b_not_run else [bc["review_candidates"]]),
            review_candidates_manifest_arg=[bridge["review_manifest"]] + ([] if route_b_not_run else [bc["review_manifest"]]),
            **{k: str(copy_root / v) for k, v in copy_paths.items() if (copy_root / v).exists()},
            max_adaptive_rounds=max_adaptive_rounds, case_id=case_id,
            boundary_route_b_bundle_arg=str(boundary_bundle) if boundary_bundle is not None else None,
            provisional_screen_ledger_arg=provisional_screen_path,
            provisional_screen_ledger_sha256=provisional_screen_sha256,
            identity_assessment_arg=identity_assessment_arg,
            identity_target_proteins_arg=proteins_arg if identity_assessment_arg else None)
        wm = _read(Path(wet["manifest"]))
        if wm["supported_member_union_sha256"] != sha256(union["candidate_set"]):
            raise ValueError("family/copy and Wet consumed different supported unions")
        # Add execution coverage before freezing this handoff. Do not change
        # membership or replace a stricter downstream failure/review verdict.
        wm["claim_review"] = claim_review
        wm["execution_context"] = execution_context
        if claim_review["state"] == "partial" and wm["handoff_state"] == "wet_reviewable":
            wm["handoff_state"] = "wet_reviewable_degraded_coverage"
        if wet.get("report") and wet.get("wet_candidates"):
            from .boundary_tier_review import write_wet_review_report
            current = pd.read_csv(wet["wet_candidates"], sep="\t", dtype=str).fillna("")
            pending_ids: set[str] = set()
            scope = wm.get("boundary_scope_manifest")
            if scope:
                scope_path = scope.get("path") if isinstance(scope, dict) else scope
                scope_data = _read(Path(str(scope_path)))
                pending = pd.read_csv(
                    scope_data["artifacts"]["boundary"]["path"], sep="\t", dtype=str
                ).fillna("")
                pending_ids = set(pending.get("candidate_gene_id", []))
            report = write_wet_review_report(
                wet["wet_candidates"], wet["report"]["path"],
                case_id, wm["handoff_state"], "first_stage",
                adopted_gene_ids=set(current["candidate_gene_id"]) - pending_ids,
                pending_gene_ids=pending_ids,
                execution_context=execution_context,
                route_a_rethink=wet.get("route_a_rethink"),
                next_boundary_scope_manifest=scope,
            )
            wet["report"] = wm["report"] = report
        wet["handoff_state"] = wm["handoff_state"]
        Path(wet["manifest"]).write_text(
            json.dumps(wm, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        for spec in source_hashes.values():
            if sha256(spec["path"]) != spec["sha256"]:
                raise ValueError("projection input changed during execution")
        result = {"state": "first_stage_frozen", "case_id": case_id,
                  "post_freeze_identity_assessment": identity_descriptor,
                  "claim_review": claim_review, "execution_context": execution_context,
                  "route_b_execution": execution_context["route_b_execution"],
                  "route_c_execution": execution_context["route_c_execution"],
                  "route_a_rethink": wet.get("route_a_rethink"),
                  "a": bridge, "b": bc, "union": union, "family": family, "wet": wet,
                  "second_stage": wet["second_stage_state"],
                  "files": _files(output)}
        _write(output / "completion.json", result)
        return result
    except Exception as error:
        _write(output / "failure.json", {"state": "failed", "error_type": type(error).__name__, "error": str(error)})
        raise


def finish_wet_review(handoff_arg: str, completion_sha256: str, out_arg: str, *,
                      proposals_arg: str, reviewer_id: str = "", reviewer_model: str = "",
                      verification_arg: str | None = None,
                      independent_review_arg: str | None = None,
                      executor_id: str = "", executor_model: str = "",
                      allow_provisional_wet: bool = False,
                      provisional_screen_ledger_arg: str | None = None,
                      provisional_screen_ledger_sha256: str | None = None) -> dict:
    """Import second-stage proposals and optionally consume review-inclusive Wet selection."""
    from .adaptive_synthesis import (apply_verified_refinement, consume_review_inclusive_selection,
                                     import_adaptive_synthesis_proposals)
    root = runtime_path(handoff_arg).resolve()
    anchor = root / "completion.json"
    if sha256(anchor) != completion_sha256:
        raise ValueError("first-stage completion anchor mismatch")
    first = _read(anchor)
    _verify(root, first["files"])
    # Boundary review consumes the immutable policy-bound scope through the
    # existing policy-bound selection consumer; evidence-scoped scopes do not
    # continue into unproposed lower inventory after completion.
    wet = first.get("wet", {})
    boundary_scope = wet.get("boundary_scope_manifest")
    if (wet.get("first_stage_selection_policy_version") in {"evidence_scoped_v2", "upper_evidence_scoped_v3"}
            and not boundary_scope):
        raise ValueError("no pending boundary review within the round budget; frozen Wet is retained")
    if (provisional_screen_ledger_arg or provisional_screen_ledger_sha256) and not allow_provisional_wet:
        raise ValueError("provisional screen ledger requires allow_provisional_wet")
    if allow_provisional_wet and not boundary_scope:
        raise ValueError("post-review provisional Wet requires a canonical boundary scope")
    if boundary_scope:
        from .boundary_tier_review import (
            build_boundary_review_request, consume_boundary_selection, write_wet_review_report,
        )
        if verification_arg:
            raise ValueError("boundary selection cannot bypass review through legacy verified refinement")
        scope_path = str(boundary_scope.get("path") if isinstance(boundary_scope, dict) else boundary_scope)
        scope_sha = str(wet.get("boundary_scope_sha256") or
                        (boundary_scope.get("sha256") if isinstance(boundary_scope, dict) else ""))
        if not scope_sha:
            raise ValueError("boundary review requires boundary_scope_sha256")
        adaptive_round = int(wet.get("adaptive_round", 0))
        max_adaptive_rounds = int(wet.get("max_adaptive_rounds", 1))
        request = wet.get("adaptive_synthesis_request")
        if not request or adaptive_round >= max_adaptive_rounds:
            raise ValueError("no pending boundary review within the round budget; frozen Wet is retained")
        request_data = _read(runtime_path(str(request)))
        scope_data = _read(runtime_path(scope_path))
        requested_scope = request_data.get("scope_manifest", {})
        if (request_data.get("contract") != "boundary-tier-review-request-1.0.0"
                or request_data.get("case_id") != first["case_id"]
                or requested_scope.get("sha256") != scope_sha
                or runtime_path(str(requested_scope.get("path", ""))).resolve() != runtime_path(scope_path).resolve()
                or request_data.get("adaptive_round") != adaptive_round
                or request_data.get("max_adaptive_rounds") != max_adaptive_rounds
                or request_data.get("source_wet_sha256") != scope_data["artifacts"]["boundary"]["sha256"]):
            raise ValueError("boundary request does not match its scope/case/round")
        boundary_rows = pd.read_csv(scope_data["artifacts"]["boundary"]["path"], sep="\t", dtype=str).fillna("")
        if request_data.get("candidate_gene_ids") != boundary_rows["candidate_gene_id"].tolist():
            raise ValueError("boundary request roster does not match its scope")
        for ledger in first.get("deferral_ledgers", []):
            if sha256(ledger["path"]) != ledger["sha256"]:
                raise ValueError("prior boundary deferral ledger hash mismatch")
        output = runtime_path(out_arg).resolve()
        if output.is_relative_to(root):
            raise ValueError("second-stage output must not change the frozen first-stage directory")
        output.mkdir(parents=True, exist_ok=False)
        raw = output / "raw_proposals.tsv"
        shutil.copyfile(runtime_path(proposals_arg), raw)
        independent = None
        if independent_review_arg:
            independent = output / "independent_review.tsv"
            shutil.copyfile(runtime_path(independent_review_arg), independent)
        try:
            selected = consume_boundary_selection(
                scope_path, str(raw), str(output / "boundary_selection"),
                scope_sha256=scope_sha, case_id=first["case_id"],
                reviewer_id=reviewer_id, reviewer_model=reviewer_model,
                executor_id=executor_id, executor_model=executor_model,
                independent_review_arg=str(independent) if independent else None,
                allow_provisional_wet=allow_provisional_wet,
                provisional_screen_ledger_arg=(str(runtime_path(provisional_screen_ledger_arg))
                                               if provisional_screen_ledger_arg else None),
                provisional_screen_ledger_sha256=provisional_screen_ledger_sha256,
                provisional_wet_context=first.get("provisional_wet_context"),
                prior_deferral_ledgers=first.get("deferral_ledgers", []),
                review_history_state=wet.get("review_history_state", "complete"),
            )
            _verify(root, first["files"])
            next_request = None
            next_scope = selected.get("next_scope_manifest") or {}
            if next_scope.get("path"):
                next_request = build_boundary_review_request(
                    str(next_scope["path"]), str(output / "next_boundary_request"),
                    scope_sha256=str(next_scope["sha256"]), case_id=first["case_id"],
                    adaptive_round=adaptive_round + 1,
                    max_adaptive_rounds=max_adaptive_rounds,
                )
            history_state = ("incomplete" if wet.get("review_history_state") == "incomplete"
                             or selected.get("review_history_state") == "incomplete" else "complete")
            next_wet = {
                "wet_candidates": selected["wet_candidates"]["path"],
                "wet_candidates_sha256": selected["wet_candidates"]["sha256"],
                "candidate_count": selected["output_candidate_count"],
                "boundary_scope_manifest": next_scope,
                "boundary_scope_sha256": next_scope.get("sha256"),
                "adaptive_synthesis_request": next_request.get("request") if next_request and next_request.get("request") else None,
                "adaptive_round": adaptive_round + 1,
                "max_adaptive_rounds": max_adaptive_rounds,
                "review_history_state": history_state,
                **({"assessment_mode": selected["assessment_mode"],
                    "review_history_basis": selected["review_history_basis"]}
                   if "assessment_mode" in selected else {}),
                "second_stage_state": (next_request["state"] if next_request else
                                        "review_incomplete_retry_original_scope"
                                        if history_state == "incomplete" and selected.get("resume_scope_manifest") else
                                        "not_required"),
                "first_stage_selection_policy_version": selected.get("policy_version"),
                **({"eligible_first_stage_candidate_count": len(selected["eligible_gene_ids"]),
                    "active_first_stage_candidate_count": len(selected["active_gene_ids"]),
                    "upper_priority_applied": selected["upper_priority_applied"]}
                   if "active_gene_ids" in selected else {}),
                "resume_scope_manifest": selected.get("resume_scope_manifest"),
                "resume_from_completion": ({"path": str(anchor), "sha256": completion_sha256}
                                            if selected.get("resume_scope_manifest") else None),
            }
            result_state = ("boundary_review_incomplete_retained"
                            if history_state == "incomplete" else
                            "provisional_adopted_wet" if selected["state"] == "provisional_adopted_wet" else
                            "review_required_round_budget_exhausted"
                            if next_wet["second_stage_state"] == "review_required_round_budget_exhausted" else
                            "next_boundary_review_required" if next_wet["adaptive_synthesis_request"] else
                            "completed_boundary_selection" if selected["state"] == "within_target_range" else selected["state"])
            pending_ids = set()
            if next_scope.get("path"):
                next_plan = _read(Path(next_scope["path"]))
                pending = pd.read_csv(next_plan["artifacts"]["boundary"]["path"], sep="\t", dtype=str).fillna("")
                pending_ids = set(pending.get("candidate_gene_id", []))
            decisions = pd.read_csv(selected["review_decisions"]["path"], sep="\t", dtype=str).fillna("")
            unresolved_ids = set(wet.get("unresolved_review_gene_ids", [])) | set(
                decisions.loc[decisions["applied_decision"].eq("not_evaluated"), "candidate_gene_id"])
            pending_ids |= unresolved_ids
            next_wet["unresolved_review_gene_ids"] = sorted(unresolved_ids)
            adopted = pd.read_csv(next_wet["wet_candidates"], sep="\t", dtype=str).fillna("")
            from .wet_candidate_handoff import DEFAULT_MIN_ROUTE_CANDIDATES, route_a_rethink_signal
            prior_rethink = wet.get("route_a_rethink") or {}
            rethink = {**prior_rethink, **route_a_rethink_signal(
                len(set(adopted["candidate_gene_id"]) - pending_ids),
                min_route_candidates=int(prior_rethink.get("min_route_candidates", DEFAULT_MIN_ROUTE_CANDIDATES)),
            )}
            next_wet["route_a_rethink"] = rethink
            report = write_wet_review_report(
                next_wet["wet_candidates"], str(output / "wet_review_report.md"), first["case_id"],
                result_state, "boundary_tier_review", pending_gene_ids=pending_ids,
                adopted_gene_ids=set(adopted["candidate_gene_id"]) - pending_ids,
                review_history_state=history_state,
                selection_policy_version=str(selected.get("policy_version") or ""),
                review_decisions_path=selected["review_decisions"]["path"],
                review_decisions_sha256=selected["review_decisions"]["sha256"],
                next_boundary_scope_manifest=next_scope,
                resume_scope_manifest=selected.get("resume_scope_manifest"),
                execution_context=first.get("execution_context"),
                route_a_rethink=rethink,
                provisional_wet_context=selected.get("provisional_wet_context"),
                interval_unresolved=selected.get("interval_unresolved"),
            )
            result = {
                "state": result_state,
                "first_stage_sha256": completion_sha256,
                "case_id": first["case_id"],
                "adaptive_round": adaptive_round,
                "max_adaptive_rounds": max_adaptive_rounds,
                "boundary_selection": selected,
                **({"assessment_mode": selected["assessment_mode"],
                    "review_history_basis": selected["review_history_basis"]}
                   if "assessment_mode" in selected else {}),
                "adopted_wet_candidates": selected["wet_candidates"]["path"],
                "adopted_wet_sha256": selected["wet_candidates"]["sha256"],
                "adopted_candidate_count": selected["output_candidate_count"],
                "route_a_rethink": rethink,
                "review_history_state": history_state,
                "execution_context": first.get("execution_context"),
                "deferral_ledgers": [*first.get("deferral_ledgers", []), selected["deferral_ledger"]],
                "report": report,
                "adopted_stage": "boundary_tier_review",
                "next_boundary_scope_manifest": selected.get("next_scope_manifest"),
                "resume_scope_manifest": selected.get("resume_scope_manifest"),
                "resume_from_completion": next_wet.get("resume_from_completion"),
                "resume_instruction": ("Retry from the original first-stage completion in a new output directory; this incomplete selection output is not a resume source."
                                       if next_wet.get("resume_from_completion") else None),
                "wet": next_wet,
                "files": _files(output),
            }
            if selected.get("provisional_wet_context"):
                result.update(provisional_wet_context=selected["provisional_wet_context"],
                              provisional_wet=selected["provisional_wet"],
                              interval_unresolved=selected.get("interval_unresolved"))
                next_wet.update(provisional_wet=selected["provisional_wet"],
                                interval_unresolved=selected.get("interval_unresolved"))
            _write(output / "completion.json", result)
            return result
        except Exception as error:
            _write(output / "failure.json", {"state": "failed", "error_type": type(error).__name__, "error": str(error)})
            raise
    request = first["wet"].get("adaptive_synthesis_request")
    if not request or not reviewer_id.strip() or not reviewer_model.strip():
        raise ValueError("second-stage request and actual reviewer identity/model required")
    output = runtime_path(out_arg).resolve()
    if output.is_relative_to(root):
        raise ValueError("second-stage output must not change the frozen first-stage directory")
    output.mkdir(parents=True, exist_ok=False)
    raw = output / "raw_proposals.tsv"
    shutil.copyfile(runtime_path(proposals_arg), raw)
    independent = None
    if independent_review_arg:
        independent = output / "independent_review.tsv"
        shutil.copyfile(runtime_path(independent_review_arg), independent)
    verification = None
    if verification_arg:
        verification = output / "verification.tsv"
        shutil.copyfile(runtime_path(verification_arg), verification)
    try:
        proposal_columns = set(pd.read_csv(raw, sep="\t", nrows=0).columns)
        selection_input = bool(proposal_columns.intersection({"decision", "recommendation", "selection", "disposition", "action"}))
        if selection_input and verification_arg:
            raise ValueError("verified refinement and review-inclusive selection cannot be combined in one stage")
        if {"ambiguity_kind", "evidence_kind"}.issubset(proposal_columns):
            imported = import_adaptive_synthesis_proposals(str(raw), request, str(output / "import"),
                provider=reviewer_id, model=reviewer_model, case_id=first["case_id"])
        elif selection_input:
            imported = {"proposals": str(raw), "proposal_count": len(pd.read_csv(raw, sep="\t")),
                        "provider_state": "selection_consumer_input"}
        else:
            raise ValueError("second-stage proposals lack a recognized proposal or selection schema")
        refined = None
        selection = None
        # Selection is opt-in by schema: legacy proposal-only TSVs retain the
        # historical review-only behavior and never alter Wet membership.
        if selection_input:
            selection = consume_review_inclusive_selection(
                first["wet"]["wet_candidates"], str(raw), str(output / "selection"),
                reviewer_id=reviewer_id, reviewer_model=reviewer_model,
                executor_id=executor_id, executor_model=executor_model,
                independent_review_arg=str(independent) if independent else None, case_id=first["case_id"],
                first_stage_sha256=completion_sha256)
        else:
            refined = apply_verified_refinement(first["wet"]["wet_candidates"],
                str(verification) if verification else None, str(output / "refinement"), case_id=first["case_id"])
        _verify(root, first["files"])
        result = {"state": "completed_verified_refinement" if verification else "completed_review_only",
                  "first_stage_sha256": completion_sha256, "reviewer_id": reviewer_id,
                  "reviewer_model": reviewer_model, "import": imported, "refinement": refined,
                  "files": _files(output)}
        if selection:
            result["selection"] = selection
            result["state"] = ("completed_review_inclusive_selection"
                                if selection["state"] == "completed"
                                else selection["state"])
            result["adopted_wet_candidates"] = selection["wet_candidates"]
            result["adopted_stage"] = "review_inclusive_wet_selection"
        _write(output / "completion.json", result)
        return result
    except Exception as error:
        _write(output / "failure.json", {"state": "failed", "error_type": type(error).__name__, "error": str(error)})
        raise
