"""Resume practical A/B packets plus optional Direct C through canonical Wet.

No model calls: prepare -> real B responses -> review -> finish.
External metadata retrieval is opt-in and sealed before offline verification.
Initial-candidate policy does not require independent content review; legacy
content policies retain actual reviews or explicit no-source non-execution.
Each command consumes a hash-anchored previous stage and writes
a new directory. The final projection/family/copy/Wet owner is handoff_sessions.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import shutil

from packet_to_canonical import import_a_payload, read_packet_b
from candidate_qtl.route_b_response_validation import route_b_response_completed, route_b_execution_breakdown, validate_route_b_not_run


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def write(path: Path, value: dict) -> None:
    with path.open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n")


def seal(root: Path, value: dict) -> dict:
    value["files"] = {str(p.relative_to(root)): digest(p)
                      for p in sorted(root.rglob("*")) if p.is_file()}
    anchor = root / "stage.json"
    write(anchor, value)
    return {**value, "manifest": str(anchor), "sha256": digest(anchor)}


def load_stage(anchor: Path, expected: str, state: str | tuple[str, ...]) -> dict:
    anchor = anchor.resolve()
    if digest(anchor) != expected:
        raise ValueError("stage anchor mismatch")
    value = read(anchor)
    if value["state"] not in ((state,) if isinstance(state, str) else state):
        raise ValueError(f"expected stage {state}")
    for name, sha in value["files"].items():
        path = (anchor.parent / name).resolve()
        if not path.is_relative_to(anchor.parent) or digest(path) != sha:
            raise ValueError(f"stage artifact mismatch: {name}")
    return value


def snapshot(source: Path, target: Path) -> str:
    if source.is_symlink() or not source.exists():
        raise ValueError(f"missing or symlink input: {source}")
    if source.is_dir():
        if any(p.is_symlink() for p in source.rglob("*")):
            raise ValueError("symlink within input directory")
        shutil.copytree(source, target)
    else:
        shutil.copyfile(source, target)
    return str(target)


def fresh_output(out: Path, inputs: list[Path]) -> Path:
    out = out.resolve()
    if any(out.is_relative_to(p.resolve()) or p.resolve().is_relative_to(out) for p in inputs):
        raise ValueError("output must be separate from input stages/packets")
    out.mkdir(parents=True, exist_ok=False)
    return out


def _live_bibliography(table: str | Path, out: Path, *, id_column: str, source_column: str,
                       focal_species: str = "", cache: str | None = None) -> Path:
    """Retrieve ID/title metadata for every cited paper once and return its records.

    Only bibliographic metadata is fetched (the initial_bibliography_v1 check);
    the returned records become a sealed cache that later offline verification
    consumes, so a provider failure stays explicit rather than becoming absence.
    """
    from candidate_qtl.bibliography import INITIAL_POLICY
    from candidate_qtl.source_verify import source_verify

    result = source_verify(str(table), str(out), cache_arg=cache, allow_external_api=True, id_column=id_column,
                           source_column=source_column, verification_policy=INITIAL_POLICY,
                           focal_species=focal_species)
    records = Path(str(result.get("bibliography_records", out / "bibliography_records.json")))
    if not records.is_file():
        raise ValueError("live bibliography lookup produced no metadata records")
    return records


def _merged_source_cache(target: Path, existing: str | None, live: list[Path]) -> str:
    """Copy a supplied source cache and live metadata records into one cache directory."""
    target.mkdir()
    if existing:
        source = Path(existing)
        for item in (sorted(source.glob("*.json")) if source.is_dir() else [source]):
            snapshot(item, target / ("cached-" + item.name))
    for index, records in enumerate(live):
        snapshot(records, target / f"live-{index:02d}-{records.name}")
    return str(target)


ROUTE_B_NOT_RUN_STATE = "not_run_optional"


def _route_b_not_run_context(reason: str, responses: Path) -> dict:
    """Record that optional Route B was deliberately not run (A/C-first, 2026-11-08).

    Every prompt must be answered as not evaluated; a completed B answer
    cannot be relabelled as "not run".
    """
    reason = str(reason or "").strip()
    if not reason:
        raise ValueError("--route-b-not-run requires a non-empty reason")
    payload = json.loads(responses.read_bytes().decode("utf-8-sig"))
    rows = payload.get("candidates") if isinstance(payload, dict) else None
    validate_route_b_not_run(rows)
    return {"state": ROUTE_B_NOT_RUN_STATE, "reason": reason, "prompt_count": len(rows),
            "completed_count": 0, "not_evaluated_count": len(rows)}


def _route_b_execution_context(
    manifest_arg: Path | None,
    prepared: dict,
    responses: Path,
    out: Path,
) -> dict:
    """Validate and snapshot the real Route-B shard execution manifest."""
    if manifest_arg is None:
        return {
            "state": "not_recorded",
            "reason": "route_b_execution_manifest was not supplied",
        }
    manifest_path = manifest_arg.resolve()
    manifest = read(manifest_path)
    if manifest.get("contract") != "agentic-route-b-execution-1.1.0":
        raise ValueError("unsupported Route B execution manifest contract")
    if manifest.get("case_id") != prepared["case_id"]:
        raise ValueError("Route B execution manifest case mismatch")
    prompts = Path(prepared["bundle"]["prompts"]).resolve()
    prompt_bytes = prompts.read_bytes()
    if manifest.get("prompt_sha256") != hashlib.sha256(prompt_bytes).hexdigest():
        raise ValueError("Route B execution manifest prompt hash mismatch")
    response_bytes = responses.read_bytes()
    if manifest.get("response_sha256") != hashlib.sha256(response_bytes).hexdigest():
        raise ValueError("Route B execution manifest response hash mismatch")
    prompt_text = prompt_bytes.decode("utf-8")
    prompt_count = sum(1 for line in prompt_text.splitlines() if line.strip())
    if manifest.get("prompt_count") != prompt_count:
        raise ValueError("Route B execution manifest prompt count mismatch")
    counts = (manifest.get("completed_count"), manifest.get("not_evaluated_count"))
    if any(type(value) is not int or value < 0 for value in counts) or sum(counts) != prompt_count:
        raise ValueError("Route B execution manifest coverage counts are inconsistent")
    expected_status = "not_evaluated" if prompt_count == 0 else ("complete" if counts[0] == prompt_count else "partial")
    if manifest.get("completion_status") != expected_status:
        raise ValueError("Route B execution manifest completion status is inconsistent")
    provenance = manifest.get("partial_response_provenance")
    if not isinstance(provenance, list) or (prompt_count > 0 and not provenance) or any(not isinstance(row, dict) for row in provenance):
        raise ValueError("Route B execution manifest requires attempt provenance")
    for row in provenance:
        if type(row.get("attempt_index")) is not int or row["attempt_index"] < 0:
            raise ValueError("Route B execution manifest attempt index is invalid")
        if not str(row.get("executor", "")).strip() or not str(row.get("model", "")).strip():
            raise ValueError("Route B execution manifest requires executor/model provenance")
        attempt_path = Path(str(row.get("file", ""))).resolve()
        attempt_digest = str(row.get("sha256", "")).removeprefix("sha256:")
        if not attempt_path.is_file() or len(attempt_digest) != 64 or digest(attempt_path) != attempt_digest:
            raise ValueError("Route B execution attempt provenance hash mismatch")
    retry_history = manifest.get("retry_history")
    if not isinstance(retry_history, list):
        raise ValueError("Route B execution manifest requires retry_history")
    for row in retry_history:
        if not isinstance(row, dict) or not all(str(row.get(key, "")).strip() for key in (
            "prompt_id", "route_b_unit_id", "replaced_state", "replaced_reason",
            "replaced_response_file", "replaced_response_sha256", "replacement_response_file",
            "replacement_response_sha256", "replacement_executor", "replacement_model",
        )):
            raise ValueError("Route B execution retry history is incomplete")
        for path_key, digest_key in (("replaced_response_file", "replaced_response_sha256"), ("replacement_response_file", "replacement_response_sha256")):
            retry_path = Path(str(row[path_key])).resolve()
            retry_digest = str(row[digest_key]).removeprefix("sha256:")
            if not retry_path.is_file() or len(retry_digest) != 64 or digest(retry_path) != retry_digest:
                raise ValueError("Route B execution retry history hash mismatch")
    roster = manifest.get("response_roster")
    if not isinstance(roster, list) or len(roster) != prompt_count:
        raise ValueError("Route B execution manifest response roster is incomplete")
    response_payload = json.loads(response_bytes.decode("utf-8-sig"))
    if response_payload.get("case_id") != prepared["case_id"]:
        raise ValueError("Route B response case_id does not match the frozen case")
    actual_rows = response_payload.get("candidates")
    if not isinstance(actual_rows, list) or len(actual_rows) != prompt_count:
        raise ValueError("Route B response candidate roster is incomplete")
    actual_by_id: dict[str, dict] = {}
    for row in actual_rows:
        if not isinstance(row, dict) or not row.get("prompt_id") or str(row["prompt_id"]) in actual_by_id:
            raise ValueError("Route B response candidate IDs are invalid")
        actual_by_id[str(row["prompt_id"])] = row
    actual_completed = sum(route_b_response_completed(row) for row in actual_rows)
    if manifest.get("completed_count") != actual_completed or manifest.get("not_evaluated_count") != prompt_count - actual_completed:
        raise ValueError("Route B execution manifest counts do not match response rows")
    expected_rows: list[dict] = []
    for line in prompt_text.splitlines():
        if line.strip():
            expected_rows.append(json.loads(line))
    for index, expected in enumerate(expected_rows):
        row = roster[index]
        if not isinstance(row, dict):
            raise ValueError("Route B execution manifest response roster row is invalid")
        prompt_id = str(expected["prompt_id"])
        if str(row.get("prompt_id", "")) != prompt_id or str(row.get("route_b_unit_id", "")) != str(expected["route_b_unit_id"]):
            raise ValueError("Route B execution manifest prompt/unit mapping is invalid")
        actual = actual_by_id.get(prompt_id)
        if actual is None or row.get("parse_status") != actual.get("parse_status"):
            raise ValueError("Route B execution manifest parse status is not bound to the response")
        if str(actual.get("route_b_unit_id", "")) != str(expected["route_b_unit_id"]):
            raise ValueError("Route B response unit_id is not bound to the frozen prompt")
        raw_hash = hashlib.sha256(str(actual.get("raw_text", "")).encode()).hexdigest()
        if row.get("raw_text_sha256") != raw_hash:
            raise ValueError("Route B execution manifest raw response hash is not bound to the response")
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
    target = out / "b_execution_manifest.json"
    state_counts = route_b_execution_breakdown(actual_rows, manifest=manifest)
    snapshot(manifest_path, target)
    return {
        **manifest,
        "state": manifest["completion_status"],
        "response_state_counts": state_counts,
        "path": str(target),
        "sha256": digest(target),
        "source_path": str(manifest_path),
        "source_sha256": digest(manifest_path),
    }


def finish_report(out: Path, case_id: str, state: str, handoff: dict) -> str:
    """Write a readable, provenance-linked summary from the canonical Wet handoff."""
    wet = handoff.get("wet") or {}
    handoff_state = wet.get("handoff_state", "unavailable")
    second_stage = wet.get("second_stage_state", "unavailable")
    wet_path = wet.get("wet_candidates")
    count = None
    route_counts = {}
    membership = {}
    rows = []
    if wet_path and Path(wet_path).is_file():
        wet_file = Path(wet_path)
        wet_before = digest(wet_file)
        expected_wet = str(wet.get("wet_candidates_sha256", ""))
        if expected_wet and wet_before != expected_wet.removeprefix("sha256:"):
            raise ValueError("Wet artifact changed before report generation")
        with wet_file.open(encoding="utf-8-sig", newline="") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        count = len(rows)
        for column in ("route_a_supported", "route_b_supported", "route_c_supported"):
            if rows and column in rows[0]:
                route_counts[column] = sum(str(row.get(column, "")).strip().casefold() in
                                           {"true", "1", "yes"} for row in rows)
        if rows and "membership_basis" in rows[0]:
            for row in rows:
                key = row.get("membership_basis", "") or "(empty)"
                membership[key] = membership.get(key, 0) + 1
        if digest(wet_file) != wet_before:
            raise ValueError("Wet artifact changed while reading report inputs")
    available = count is not None
    next_action = {
        "wet_reviewable": "Inspect the retained Wet genes and plan experimental confirmation.",
        "wet_reviewable_degraded_coverage": "Inspect the retained Wet genes and unresolved evidence coverage before experimental confirmation.",
        "provisional_adopted_wet": "Inspect the provisional Wet genes and the complete interval_unresolved ledger before experimental confirmation.",
        "adaptive_refinement_requested": "Complete the requested bounded refinement, then re-gate Wet.",
        "review_required_over_budget": "Complete the required second-stage review before adopting a set.",
        "review_required_no_supported_candidates": "Review the empty supported set and decide whether additional evidence is needed.",
        "boundary_review_requested": "Complete the requested boundary review using the recorded scope.",
        "review_required_round_budget_exhausted": "Review the exhausted boundary round and decide the next evidence step.",
    }.get(handoff_state, "Inspect the handoff state and resolve its stated boundary before proceeding.")
    if second_stage == "requested" and handoff_state not in {"boundary_review_requested", "adaptive_refinement_requested"}:
        next_action = "Complete the requested second-stage boundary review using the recorded scope."
    elif wet.get("adaptive_synthesis_request") and handoff_state != "adaptive_refinement_requested":
        next_action = "Complete the requested adaptive synthesis, then re-gate Wet."

    def inline(value) -> str:
        return "`" + str(value).replace("`", "\\`").replace("\n", " ") + "`"

    def cell(value: object, limit: int = 180) -> str:
        value = str(value or "").replace("\r", " ").replace("\n", " ").replace("|", "\\|")
        if not value:
            return "(not recorded)"
        return value if len(value) <= limit else value[:limit - 1] + "…"

    def artifact_path(value: object) -> str:
        if isinstance(value, dict):
            return str(value.get("path") or "")
        return str(value or "")

    handoff_manifest_path = artifact_path(wet.get("manifest"))
    membership_manifest_path = artifact_path(wet.get("wet_candidates_manifest"))
    manifest_data = {}
    membership_manifest_file = Path(membership_manifest_path) if membership_manifest_path else None
    membership_before = ""
    if membership_manifest_file and membership_manifest_file.is_file():
        membership_before = digest(membership_manifest_file)
        try:
            manifest_data = read(membership_manifest_file)
        except (OSError, ValueError):
            manifest_data = {}
        if digest(membership_manifest_file) != membership_before:
            raise ValueError("Wet membership manifest changed while reading report inputs")
    trace_path = artifact_path(manifest_data.get("candidate_trace"))
    trace_rows = []
    if trace_path and Path(trace_path).is_file():
        trace_file = Path(trace_path)
        trace_before = digest(trace_file)
        expected_trace = str(manifest_data.get("candidate_trace_sha256", ""))
        if expected_trace and trace_before != expected_trace.removeprefix("sha256:"):
            raise ValueError("Wet candidate trace changed before report generation")
        with trace_file.open(encoding="utf-8-sig", newline="") as stream:
            trace_rows = list(csv.DictReader(stream, delimiter="\t"))
        if digest(trace_file) != trace_before:
            raise ValueError("Wet candidate trace changed while reading report inputs")

    def hash_line(label: str, path: str, recorded_hash: object = "") -> str | None:
        if not path:
            return None
        artifact = Path(path)
        if not artifact.is_file():
            recorded = inline(recorded_hash or "unavailable")
            return f"- {label}: {inline(path)} (MISSING; recorded SHA-256 {recorded})"
        actual = digest(artifact)
        expected = str(recorded_hash or "")
        suffix = f"; manifest SHA-256 {inline(expected)}" if expected and expected != actual else ""
        return f"- {label}: {inline(path)} (SHA-256 {inline(actual)}{suffix})"

    lines = [f"# Practical packet run: {inline(case_id)}", "",
             f"- Execution state: {inline(state)}",
             f"- Wet handoff state: {inline(handoff_state)}",
             f"- Second-stage state: {inline(second_stage)}"]
    if any("initial_bibliography_v1" in row.get("candidate_admission_policy", "") for row in rows):
        lines += ["", "第一段の本候補として採用。初期 Tier と元の A/B の採用理由を第二段の絞り込みへ渡します。", ""]
    c_execution = handoff.get("route_c_execution") or {"state": "not_evaluated_no_input"}
    lines.append(f"- Direct C execution: {inline(c_execution['state'])}")
    if c_execution.get("reason"):
        lines.append(f"- Direct C note: {cell(c_execution['reason'])}")
    b_execution = handoff.get("route_b_execution") or {"state": "not_recorded"}
    lines.append(
        f"- Route B execution: {inline(b_execution.get('state', 'not_recorded'))}; "
        f"completed={b_execution.get('completed_count', 0)}, "
        f"not_evaluated={b_execution.get('not_evaluated_count', 0)}"
    )
    if b_execution.get("reason"):
        lines.append(f"- Route B note: {cell(b_execution['reason'])}")
    if b_execution.get("response_state_counts"):
        counts = b_execution["response_state_counts"]
        lines.append("- Route B response states: " + ", ".join(f"{cell(key)}={cell(value)}" for key, value in counts.items()))
        lines.append("- Bの内訳は回答の受領・形式状態です。読了や文献探索の完了を示すものではありません。")
    b_manifest_path = str(b_execution.get("path", "") or "")
    if b_manifest_path:
        manifest_line = hash_line("Route B execution manifest", b_manifest_path, b_execution.get("sha256", ""))
        if manifest_line:
            lines.append(manifest_line)
    c_status_path = str(c_execution.get("status_path", "") or c_execution.get("artifact_path", "") or "")
    if c_status_path:
        status_line = hash_line(
            "Route C execution status", c_status_path,
            c_execution.get("status_sha256", c_execution.get("artifact_sha256", "")),
        )
        if status_line:
            lines.append(status_line)
    c_members_path = str(c_execution.get("members_path", "") or "")
    if c_members_path:
        members_line = hash_line("Route C members", c_members_path, c_execution.get("members_sha256", ""))
        if members_line:
            lines.append(members_line)
    c_producer_path = str(c_execution.get("manifest_path", "") or "")
    if c_producer_path:
        producer_line = hash_line("Route C producer manifest", c_producer_path, c_execution.get("manifest_sha256", ""))
        if producer_line:
            lines.append(producer_line)
    if available:
        lines += [f"- Canonical Wet: present; {inline(wet_path)}",
                  f"- Wet SHA-256: {inline(digest(Path(wet_path)))}",
                  f"- Candidate rows: {count}",
                  f"- Current Wet stage: {inline(handoff.get('adopted_stage') or 'first_stage')} (this is the stage represented by the TSV)"]
        if route_counts:
            lines.append("- Route support counts: " + ", ".join(f"{key}={value}" for key, value in route_counts.items()))
        if membership:
            lines.append("- Membership basis: " + ", ".join(f"{key}={value}" for key, value in sorted(membership.items())))
        review_counts = (handoff.get("claim_review") or {}).get("routes") or {}
        review_coverage = []
        for route in ("a", "b"):
            counts = review_counts.get(route)
            if counts:
                execution = counts.get("review_execution_state", "executed")
                review_coverage.append(f"{route.upper()} {counts.get('returned', 0)}/{counts.get('requested', 0)} returned; {counts.get('missing', 0)} missing; {execution}")
        if review_coverage:
            lines.append("- Independent review coverage: " + "; ".join(review_coverage))

        # These fields are consumer displays only. Selection remains owned by
        # the canonical handoff, and absent data is never interpreted as a
        # negative result.
        lines += ["", "## Wet candidates", "",
                  "One row per gene in the recorded Wet TSV. Values are copied from that artifact; blank fields mean the handoff did not record them.", "",
                  "| Gene | Stage status | Membership | Relevance | Reason | Mapping | Copy state | Support | Source | Provenance |",
                  "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
        shown_rows = rows[:30] if count is not None and count > 25 else rows
        trace_by_id = {item.get("candidate_gene_id", ""): item for item in trace_rows}
        stage = handoff.get("adopted_stage") or "first_stage"
        for row in shown_rows:
            gene = row.get("candidate_gene_id") or row.get("gene_id") or "(ID not recorded)"
            disposition = trace_by_id.get(gene, {}).get("selection_disposition", "")
            stage_status = (f"adopted ({stage})" if stage != "first_stage" else
                            "pending boundary review" if disposition == "retained_pending_boundary_review" else
                            "first-stage Wet")
            membership_basis = row.get("membership_basis", "")
            if membership_basis == "route_review_candidate":
                reason = row.get("review_reasons") or row.get("tier_reason") or row.get("relevance_basis")
            elif membership_basis == "interval_review_coverage":
                reason = row.get("tier_reason") or row.get("relevance_basis")
            else:
                reason = row.get("review_reasons") or row.get("tier_reason") or row.get("relevance_basis")
            provenance = "; ".join(filter(None, (row.get("route_provenance", ""),
                row.get("route_support", ""), row.get("reference_gene", ""))))
            source = "; ".join(filter(None, (row.get("source_ids", ""),
                row.get("review_source_ids", ""), row.get("review_verification_statuses", ""))))
            lines.append("| " + " | ".join((cell(gene, 90), cell(stage_status),
                cell(membership_basis),
                cell("; ".join(filter(None, (row.get("relevance_tier", ""), row.get("relevance_basis", ""))))),
                cell(reason),
                cell(row.get("mapping_tier")), cell(row.get("copy_state")),
                cell(row.get("support_state")), cell(source), cell(provenance))) + " |")
        if count is not None and count > 25:
            lines += ["", f"Showing {len(shown_rows)} of {count} Wet genes. The complete membership and all per-gene fields are in the canonical Wet TSV above."]
        if stage != "first_stage":
            lines.append("- Final decision reasons: not provided in the supplied handoff artifacts.")

        if trace_rows:
            dispositions = {}
            for row in trace_rows:
                disposition = row.get("selection_disposition") or "(not recorded)"
                dispositions[disposition] = dispositions.get(disposition, 0) + 1
            lines += ["", "## Interval coverage", "",
                      f"- Trace rows: {len(trace_rows)}; unique genes: {len({row.get('candidate_gene_id', '') for row in trace_rows if row.get('candidate_gene_id')})}.",
                      "- Trace dispositions: " + "; ".join(f"{cell(key)}={value}" for key, value in sorted(dispositions.items())),
                      "- Trace is the first-stage coverage ledger. Its pending/unreached states do not override a later adopted Wet TSV stage. Neither state is a biological exclusion or a negative finding."]
            for disposition, title in (("retained_pending_boundary_review", "First-stage pending boundary genes"),
                                       ("not_reached_lower_tier", "First-stage lower-tier genes not reached")):
                genes = [row.get("candidate_gene_id", "") for row in trace_rows
                         if row.get("selection_disposition") == disposition and row.get("candidate_gene_id")]
                if genes:
                    shown = genes[:30]
                    suffix = f" … (+{len(genes) - len(shown)}; see trace)" if len(genes) > len(shown) else ""
                    lines.append(f"- {title}: {', '.join(inline(gene) for gene in shown)}{suffix}")
        else:
            lines += ["", "## Interval coverage", "", "- Candidate trace: unavailable; no interval-wide coverage is inferred."]

        for label, path, sha in (
            ("Handoff manifest", handoff_manifest_path, ""),
            ("Wet candidates manifest", membership_manifest_path, ""),
            ("Candidate trace", trace_path, manifest_data.get("candidate_trace_sha256", "")),
            ("Relevance inventory", artifact_path(manifest_data.get("relevance_inventory")), manifest_data.get("relevance_inventory_sha256", "")),
        ):
            entry = hash_line(label, path, sha)
            if entry:
                lines.append(entry)
    else:
        lines += ["- Canonical Wet: unavailable (not evaluated); no candidate count is inferred.",
                  "- Current Wet stage: none"]
    if wet.get("boundary_scope_manifest"):
        scope = wet["boundary_scope_manifest"]
        if isinstance(scope, dict):
            scope = scope.get("path", "")
        lines.append(f"- Boundary scope: {inline(scope)}")
    if wet.get("adaptive_synthesis_request"):
        lines.append(f"- Adaptive synthesis request: {inline(wet['adaptive_synthesis_request'])}")
    lines += ["", f"Next boundary: {next_action}", "",
              "This run reports workflow artifacts and coverage only. It makes no causal-accuracy claim.",
              "Unavailable or not evaluated does not mean zero candidates or a biological negative.", ""]
    path = out / "report.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    return str(path)


def prepare(config_path: Path, out: Path) -> dict:
    from Bio import SeqIO
    from candidate_qtl.route_b_prepare import RESOURCE_ENV, route_b_prepare
    from candidate_qtl.reference_gene_resolution import reference_gene_resolve
    from candidate_qtl.reference_extract import extract_candidates
    from candidate_qtl.reference_assets import (
        capture_shared_protein_db, extract_reference_interval,
        prepare_reference_asset, validate_shared_protein_db,
    )
    from candidate_qtl.route_b_semantic_bundle import build_route_b_semantic_bundle
    from candidate_qtl.source_verify import source_verify

    config_path = config_path.resolve()
    config = read(config_path)
    from candidate_qtl.bibliography import INITIAL_POLICY, CONTENT_POLICY
    admission_policy = config.get("candidate_admission_policy", INITIAL_POLICY)
    if admission_policy not in {INITIAL_POLICY, CONTENT_POLICY}:
        raise ValueError("unsupported candidate_admission_policy")
    path_keys = {"packet_a", "packet_b", "a_payload", "source_cache", "record_cache", "record_cache_manifest",
                 "family_index", "domains", "orthology", "annotations", "literature", "protein_db",
                 "reference_fasta", "annotation_gff3", "reference_corpus", "identity_review"}
    identities = {"a_proposer_id", "a_model", "b_proposer_id", "b_model"}
    ordinary_fields = {"case_id", "focal_species", "phenotype", "interval"}
    packet_mode = "packet_a" in config or "packet_b" in config
    ordinary_mode = "reference_fasta" in config or "annotation_gff3" in config
    required = identities | {"a_payload"}
    allowed = path_keys | identities | ordinary_fields | {
        "candidate_admission_policy", "a_payload_source", "external_retrieval"} | (
        {"seqid_map", "reference_storage", "reference_version", "reference_cache_dir",
         "interval_cache_dir", "protein_db_storage", "protein_db_version"} if ordinary_mode else set())
    if (not required.issubset(config) or set(config) - allowed
            or packet_mode == ordinary_mode
            or (packet_mode and not {"packet_a", "packet_b"}.issubset(config))
            or (ordinary_mode and not {"reference_fasta", "annotation_gff3", *ordinary_fields}.issubset(config))):
        raise ValueError("missing or unknown preparation config keys")
    a_payload_source = config.get("a_payload_source", "llm")
    if a_payload_source not in {"llm", "agent", "human"}:
        raise ValueError("a_payload_source must be llm, agent or human")
    external_retrieval = config.get("external_retrieval", False)
    if not isinstance(external_retrieval, bool):
        raise ValueError("external_retrieval must be true or false")
    if external_retrieval and not config.get("record_cache_manifest") and config.get("record_cache"):
        raise ValueError("external_retrieval builds the A record cache; omit record_cache or supply its manifest")
    reference_storage = config.get("reference_storage", "shared_readonly" if ordinary_mode else "snapshot")
    if not isinstance(reference_storage, str) or reference_storage not in {"snapshot", "shared_readonly", "managed_cache"}:
        raise ValueError("reference_storage must be snapshot, shared_readonly or managed_cache")
    managed_reference = ordinary_mode and reference_storage == "managed_cache"
    shared_reference = ordinary_mode and reference_storage in {"shared_readonly", "managed_cache"}
    if shared_reference and (not isinstance(config.get("reference_version"), str)
                             or not config["reference_version"].strip()):
        raise ValueError(f"{reference_storage} requires a non-empty reference_version")
    cache_paths = {}
    for key in ("reference_cache_dir", "interval_cache_dir"):
        if key in config:
            if not managed_reference or not isinstance(config[key], str) or not config[key].strip():
                raise ValueError(f"{key} requires managed_cache and a non-empty path")
            # Preserve the unresolved path for the asset helper's alias checks.
            cache_paths[key] = (config_path.parent / config[key]).absolute()
    if managed_reference and "reference_cache_dir" not in cache_paths:
        raise ValueError("managed_cache requires reference_cache_dir")
    protein_db_storage = config.get("protein_db_storage", "snapshot")
    if not isinstance(protein_db_storage, str) or protein_db_storage not in {"snapshot", "shared_readonly"}:
        raise ValueError("protein_db_storage must be snapshot or shared_readonly")
    shared_db = ordinary_mode and protein_db_storage == "shared_readonly"
    if shared_db and (not config.get("protein_db") or not isinstance(config.get("protein_db_version"), str)
                      or not config["protein_db_version"].strip()):
        raise ValueError("shared protein_db requires protein_db and a non-empty protein_db_version")
    if "protein_db_version" in config and not shared_db:
        raise ValueError("protein_db_version requires protein_db_storage=shared_readonly")
    for key in required | ({"packet_a", "packet_b"} if packet_mode else {"reference_fasta", "annotation_gff3"} | ordinary_fields):
        if not isinstance(config[key], str) or not config[key].strip():
            raise ValueError(f"non-empty {key} required")
    if shared_reference:
        for key in ("reference_fasta", "annotation_gff3"):
            raw_path = config_path.parent / config[key]
            if raw_path.is_symlink() or not raw_path.is_file():
                raise ValueError(f"shared reference must be a regular non-symlink input: {key}")
    paths = {key: (config_path.parent / value).resolve()
             for key, value in config.items() if key in path_keys and value}
    if "record_cache_manifest" in paths:
        if "record_cache" not in paths or not paths["record_cache"].is_dir():
            raise ValueError("record_cache_manifest requires a record_cache directory")
        cache_manifest = read(paths["record_cache_manifest"])
        if (cache_manifest.get("contract") != "route-a-reference-record-cache-1.0.0"
                or not cache_manifest.get("cache_usable_for_prepare")
                or Path(str(cache_manifest.get("records_dir", ""))).resolve() != paths["record_cache"]):
            raise ValueError("record_cache manifest is not usable for this directory")
        recorded = cache_manifest.get("records_sha256")
        if not isinstance(recorded, dict) or set(recorded) != {p.stem for p in paths["record_cache"].glob("*.json")}:
            raise ValueError("record_cache manifest record roster mismatch")
        for accession, sha in recorded.items():
            if digest(paths["record_cache"] / f"{accession}.json") != sha:
                raise ValueError(f"record_cache manifest SHA mismatch: {accession}")
        if cache_manifest.get("query_ledger"):
            from candidate_qtl.reference_record_cache import validate_reference_record_cache
            validate_reference_record_cache(str(paths["record_cache_manifest"]),
                                            records_dir_arg=str(paths["record_cache"]))
        elif cache_manifest.get("cache_complete") is False:
            raise ValueError("partial record_cache requires a frozen query ledger")
    seqid_map_arg = config.get("seqid_map") if ordinary_mode else None
    if seqid_map_arg is not None:
        if not isinstance(seqid_map_arg, str) or not seqid_map_arg.strip():
            raise ValueError("non-empty seqid_map required")
        seqid_map_path = (config_path.parent / seqid_map_arg).resolve()
        if seqid_map_path.exists():
            paths["seqid_map"] = seqid_map_path
    for key, variable in RESOURCE_ENV.items():
        if key != "reference_function" and os.environ.get(variable) and key not in paths:
            raise ValueError(f"declare {key} in config to snapshot {variable}")
    if packet_mode:
        a = read(paths["packet_a"] / "phenotype.json")
        b = read(paths["packet_b"] / "phenotype.json")
        for key in ("case_id", "focal_species", "phenotype"):
            if not a.get(key) or a[key] != b.get(key):
                raise ValueError(f"A/B packet mismatch: {key}")
        packet = read_packet_b(paths["packet_b"])
        case = read(paths["packet_b"] / "case.json")
        for key in ("case_id", "focal_species", "phenotype", "interval"):
            if key in case and key in b and case[key] != b[key]:
                raise ValueError(f"B case metadata mismatch: {key}")
    else:
        a = {key: config[key] for key in ("case_id", "focal_species", "phenotype")}
        b = dict(a)
    out = fresh_output(out, [*paths.values(), *cache_paths.values()])
    inputs = out / "inputs"
    inputs.mkdir()
    copied = {}
    managed_asset = None
    db_asset = None
    for key, path in paths.items():
        if ordinary_mode and key in {"packet_a", "packet_b"}:
            continue
        if shared_reference and key in {"reference_fasta", "annotation_gff3"}:
            copied[key] = str(path)
        elif shared_db and key == "protein_db":
            # Capture before B0 and revalidate after it; the complete external
            # search database stays outside the case's sealed file roster.
            db_asset = capture_shared_protein_db(str(config_path.parent / config[key]),
                                                 version=config["protein_db_version"])
            copied[key] = str(path)
        elif key == "protein_db":
            files = sorted(path.parent.glob(path.name + ".*"))
            if not files:
                raise ValueError("protein_db must name a formatted BLAST database prefix")
            target = inputs / "protein_db"
            target.mkdir()
            for file in files:
                snapshot(file, target / file.name)
            copied[key] = str(target / path.name)
        else:
            copied[key] = snapshot(path, inputs / (key if path.is_dir() else key + path.suffix))
    # A reusable family index may encode its reference_genome scope only in
    # its adjacent manifest. Snapshot that sidecar before B0 and A handoff.
    if "family_index" in paths:
        companion = paths["family_index"].parent / "family_index_manifest.json"
        if companion.is_file():
            target = Path(copied["family_index"]).parent / "family_index_manifest.json"
            snapshot(companion, target)
            copied["family_index_manifest"] = str(target)
    if "record_cache" not in copied:
        write(inputs / "records.json", {})
        copied["record_cache"] = str(inputs / "records.json")
    if ordinary_mode:
        if managed_reference:
            managed_asset = prepare_reference_asset(
                str(config_path.parent / config["reference_fasta"]),
                str(config_path.parent / config["annotation_gff3"]),
                str(cache_paths["reference_cache_dir"]), reference_version=config["reference_version"])
            copied["reference_fasta"] = managed_asset["fasta"]
            copied["annotation_gff3"] = managed_asset["gff3"]
        a_root = inputs / "packet_a"
        b_root = inputs / "packet_b"
        a_root.mkdir()
        b_root.mkdir()
        write(a_root / "phenotype.json", a)
        write(b_root / "phenotype.json", b)
        write(b_root / "case.json", {**b, "interval": config["interval"]})
        if managed_reference:
            extracted = extract_reference_interval(
                managed_asset, config["interval"], str(b_root),
                seqid_map_arg=copied.get("seqid_map", config.get("seqid_map")),
                interval_cache_dir_arg=(str(cache_paths["interval_cache_dir"])
                                        if "interval_cache_dir" in cache_paths else None))
        else:
            extracted = extract_candidates(copied["reference_fasta"], copied["annotation_gff3"],
                                           config["interval"], str(b_root),
                                           seqid_map_arg=copied.get("seqid_map", config.get("seqid_map")),
                                           **({"read_only_reference": True, "full_hash": True}
                                              if shared_reference else {}))
        if shared_reference:
            # The sealed stage owns the extracted bytes. External paths below
            # are lineage only: review/finish never reread the shared genome.
            index = Path(copied["reference_fasta"] + ".fai")
            index_snapshot = inputs / "reference.fai"
            snapshot(index, index_snapshot)
            if digest(index_snapshot) != extracted["fasta_index_fingerprint"]["sha256"]:
                raise ValueError("shared reference index changed after extraction")
            reference_manifest = inputs / "reference_asset.json"
            managed_lineage = {}
            if managed_asset is not None:
                asset_snapshot = inputs / "managed_reference_manifest.json"
                snapshot(Path(managed_asset["manifest"]), asset_snapshot)
                if digest(asset_snapshot) != managed_asset["manifest_sha256"]:
                    raise ValueError("managed reference manifest changed after extraction")
                managed_lineage = {
                    "managed_asset": {"source_manifest": managed_asset["manifest"],
                                      "snapshot": str(asset_snapshot),
                                      "sha256": managed_asset["manifest_sha256"],
                                      "cache_hit": managed_asset["cache_hit"]},
                    "original_sources": managed_asset["source_fingerprints"],
                    "interval_cache": extracted["interval_cache"],
                }
            write(reference_manifest, {
                "contract": "ordinary-shared-reference-1.0.0",
                "storage": reference_storage, "version": config["reference_version"],
                "interval": config["interval"], "gff_seqid": extracted["gff_seqid"],
                "fasta_seqid": extracted["fasta_seqid"],
                "fasta": extracted["fasta_fingerprint"], "gff3": extracted["gff3_fingerprint"],
                "index": {"source": str(index), "snapshot": str(index_snapshot),
                          "sha256": digest(index_snapshot)},
                "extraction_manifest": extracted["manifest"],
                "extraction_manifest_sha256": digest(Path(extracted["manifest"])),
                "resume_policy": "validate sealed interval artifacts; external references are lineage only",
                **managed_lineage,
            })
            copied["reference_asset_manifest"] = str(reference_manifest)
        if "reference_corpus" in copied:
            snapshot(Path(copied["reference_corpus"]), b_root / "reference_corpus.tsv")
        copied["packet_a"] = str(a_root)
        copied["packet_b"] = str(b_root)
        copied["extraction_manifest"] = extracted["manifest"]
    else:
        b_root = Path(copied["packet_b"])
    if packet_mode and not (b_root / "reference_corpus.tsv").is_file():
        raise ValueError("packet_b/reference_corpus.tsv or reference_corpus config input required")
    packet = read_packet_b(b_root)
    if str(packet.get("case_id", "")).strip() != str(a.get("case_id", "")).strip():
        raise ValueError("B packet case_id does not match the A phenotype case")
    if ordinary_mode and str(packet.get("interval", "")).strip() != str(config.get("interval", "")).strip():
        raise ValueError("B packet interval does not match the ordinary QTL config")
    proteins = b_root / "candidate_proteins.fa"
    with proteins.open(encoding="utf-8") as stream:
        protein_ids = [record.id for record in SeqIO.parse(stream, "fasta")]
    if len(protein_ids) != len(set(protein_ids)):
        raise ValueError("candidate protein FASTA requires unique gene IDs")
    if packet_mode and set(protein_ids) != set(packet["gene_ids"]):
        raise ValueError("packet protein roster must match candidate gene IDs exactly")
    if ordinary_mode and not set(protein_ids).issubset(set(packet["gene_ids"])):
        raise ValueError("extracted protein IDs must be a subset of the reference candidate roster")
    b0 = route_b_prepare(str(b_root / "candidates.tsv"), str(out / "b0"),
        proteins_arg=str(proteins) if protein_ids else None,
        reference_function_args=([str(b_root / "reference_corpus.tsv")]
                                 if (b_root / "reference_corpus.tsv").is_file() else []),
        auto_local_db=False, case_id=a["case_id"],
        **{key + "_arg": copied.get(key) for key in
           ("family_index", "domains", "orthology", "annotations", "protein_db")})
    if db_asset is not None:
        validate_shared_protein_db(db_asset)
        db_manifest = inputs / "protein_db_asset.json"
        write(db_manifest, db_asset)
        copied["protein_db_asset_manifest"] = str(db_manifest)
    bundle = build_route_b_semantic_bundle(b0["interval_gene_atlas"], str(out / "bundle"),
        phenotype=a["phenotype"], literature_arg=copied.get("literature"),
        lookup_arg=b0.get("interval_reverse_lookup"),
        identity_hits_arg=b0.get("identity_hits"), case_id=a["case_id"],
        b_review_membership_policy=INITIAL_POLICY if admission_policy == INITIAL_POLICY else "source_function_v1",
        focal_species=a["focal_species"])
    imported = import_a_payload(copied["a_payload"], out / "a_import",
                                case_id=a["case_id"], model=config["a_model"],
                                payload_source=a_payload_source)
    retrieval = {"enabled": external_retrieval}
    if external_retrieval and "record_cache_manifest" not in copied:
        # Retrieve from the exact table this stage imported, so the cache can
        # never describe a differently imported A table.
        from candidate_qtl.reference_record_cache import build_reference_record_cache
        built = build_reference_record_cache(imported["mechanisms"], str(inputs / "record_cache_live"),
                                             allow_external_api=True)
        if not built.get("cache_usable_for_prepare"):
            raise ValueError("live A reference record retrieval did not produce a usable cache")
        copied["record_cache"] = built["records_dir"]
        retrieval["record_cache"] = {"manifest": built["manifest"], "sha256": built["manifest_sha256"],
                                     "status_counts": built.get("status_counts"),
                                     "cache_complete": built.get("cache_complete")}
    if "record_cache_manifest" in copied:
        if digest(Path(imported["mechanisms"])) != cache_manifest.get("mechanisms_sha256"):
            raise ValueError("record_cache was retrieved for a different frozen A hypothesis table "
                             f"(prepare imports A with payload_source={a_payload_source}; build the cache "
                             "from that table, or set external_retrieval to let prepare build it)")
        if cache_manifest.get("query_ledger"):
            ledger_target = inputs / "record_cache_query_ledger.json"
            copied["record_cache_query_ledger"] = snapshot(Path(cache_manifest["query_ledger"]), ledger_target)
            relocated = {**cache_manifest, "records_dir": copied["record_cache"],
                         "mechanisms": imported["mechanisms"], "query_ledger": str(ledger_target),
                         "source_manifest": str(paths["record_cache_manifest"]),
                         "source_manifest_sha256": digest(paths["record_cache_manifest"])}
            Path(copied["record_cache_manifest"]).write_text(
                json.dumps(relocated, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            from candidate_qtl.reference_record_cache import validate_reference_record_cache
            validate_reference_record_cache(copied["record_cache_manifest"],
                                            mechanisms_arg=imported["mechanisms"],
                                            records_dir_arg=copied["record_cache"])
    resolved = reference_gene_resolve(imported["mechanisms"], str(out / "a_resolved"),
        record_cache_arg=copied["record_cache"], allow_external_api=False, focal_species=a["focal_species"],
        identity_review_arg=copied.get("identity_review"))
    if external_retrieval:
        # Cache-first ID/title metadata for A's cited papers, retrieved once and
        # sealed with this stage; the verification below stays offline.
        live = _live_bibliography(resolved["mechanisms_resolved"], out / "a_sources_live",
                                  id_column="hypothesis_id", source_column="source_ids",
                                  focal_species=a["focal_species"], cache=copied.get("source_cache"))
        copied["source_cache"] = _merged_source_cache(inputs / "source_cache_live",
                                                      copied.get("source_cache"), [live])
        retrieval["a_bibliography"] = {"records": str(live), "sha256": digest(live)}
    sources = source_verify(resolved["mechanisms_resolved"], str(out / "a_sources"),
                            cache_arg=copied.get("source_cache"), allow_external_api=False,
                            verification_policy=admission_policy, focal_species=a["focal_species"])
    return seal(out, {"state": "awaiting_b_response", "case_id": a["case_id"],
        "candidate_admission_policy": admission_policy, "external_retrieval": retrieval,
        "focal_species": a["focal_species"], "inputs": copied, "b0": b0, "bundle": bundle,
        "a_import": imported, "a_resolved": resolved, "a_sources": sources,
        "identities": {key: config[key] for key in identities},
        "ordinary_case": ({key: config[key] for key in ("case_id", "focal_species", "phenotype", "interval")}
                          | {"reference_storage": reference_storage}
                          | ({"reference_version": config["reference_version"]} if shared_reference else {})
                          | ({"seqid_map": config["seqid_map"]} if ordinary_mode and "seqid_map" in config else {})
                          if ordinary_mode else None),
        "prediction_calls": 0})


def review(prepared: Path, expected: str, responses: Path, out: Path,
           b_execution_manifest: Path | None = None, *,
           source_cache: Path | None = None, record_cache: Path | None = None,
           route_b_not_run: str | None = None) -> dict:
    from candidate_qtl.claim_review_session import build_no_source_review, prepare_session
    from candidate_qtl.route_b_semantic_import import import_route_b_semantic_proposals
    from candidate_qtl.route_b_semantic_source_verify import verify_route_b_semantic_sources
    from candidate_qtl.source_verify import source_verify

    prepared = prepared.resolve()
    state = load_stage(prepared, expected, "awaiting_b_response")
    from candidate_qtl.bibliography import INITIAL_POLICY
    admission_policy = state.get("candidate_admission_policy", "content_v1")
    supplied = {key: path.absolute() for key, path in
                (("source_cache", source_cache), ("record_cache", record_cache)) if path is not None}
    for key, path in supplied.items():
        if path.is_symlink() or not (path.is_file() or path.is_dir()):
            raise ValueError(f"{key} must be a regular file or directory")
        if key == "record_cache" and not path.is_dir():
            raise ValueError("record_cache must be a regular records directory")
        if path.is_dir() and any(p.is_symlink() or not (p.is_file() or p.is_dir())
                                 for p in path.rglob("*")):
            raise ValueError(f"{key} directory contains a symlink or non-regular input")
    guards = [prepared.parent, responses, *supplied.values()]
    if b_execution_manifest is not None:
        guards.append(b_execution_manifest)
    out = fresh_output(out, guards)
    caches = {key: state["inputs"].get(key) for key in ("source_cache", "record_cache")}
    materials = {}
    if supplied:
        (out / "inputs").mkdir()
    for key, path in supplied.items():
        target = out / "inputs" / (key if path.is_dir() else key + path.suffix)
        caches[key] = snapshot(path, target)
        materials[key] = {"source": str(path), "snapshot": str(target),
                          "files": {str(p.relative_to(target)): digest(p)
                                    for p in sorted(target.rglob("*")) if p.is_file()}
                          if target.is_dir() else {target.name: digest(target)}}
    raw = snapshot(responses, out / "b_response.json")
    if route_b_not_run is not None:
        if b_execution_manifest is not None:
            raise ValueError("--route-b-not-run and --b-execution-manifest are mutually exclusive")
        route_b_execution = _route_b_not_run_context(route_b_not_run, Path(raw))
    else:
        route_b_execution = _route_b_execution_context(
            b_execution_manifest, state, Path(raw), out
        )
    bundle = state["bundle"]
    imported = import_route_b_semantic_proposals(bundle["bundle"], bundle["prompts"],
        raw, str(out / "b_import"), case_id=state["case_id"])
    live_b = None
    if (state.get("external_retrieval") or {}).get("enabled"):
        # The prepared stage opted into live metadata retrieval, so B's cited
        # papers get the same sealed ID/title cache instead of depending on the
        # coordinator remembering a separate --source-cache step.
        live_b = _live_bibliography(imported["proposals"], out / "b_sources_live",
                                    id_column="claim_id", source_column="source_id",
                                    focal_species=state["focal_species"], cache=caches["source_cache"])
        (out / "inputs").mkdir(exist_ok=True)
        caches["source_cache"] = _merged_source_cache(out / "inputs" / "source_cache_live",
                                                      caches["source_cache"], [live_b])
        materials["b_bibliography"] = {"records": str(live_b), "sha256": digest(live_b),
                                       "snapshot": caches["source_cache"]}
    sources = verify_route_b_semantic_sources(imported["proposals"], str(out / "b_sources"),
        bundle_arg=bundle["bundle"], cache_arg=caches["source_cache"],
        allow_external_api=False, case_id=state["case_id"])
    records = Path(sources.get("source_review_records", out / "b_sources" / "source_verify" / "source_review_records.json"))
    if not records.exists():  # Empty/rejected proposals still get a real empty source capture.
        captured = source_verify(imported["proposals"], str(out / "b_source_capture"),
            cache_arg=caches["source_cache"], allow_external_api=False,
            id_column="claim_id", source_column="source_id")
        records = Path(captured["source_review_records"])
    a_sources = state["a_sources"]
    if source_cache is not None or live_b is not None:
        # Recheck sources of the exact resolved A rows, without rerunning its
        # proposal/import/reference resolution or changing the frozen IDs.
        a_sources = source_verify(state["a_resolved"]["mechanisms_resolved"], str(out / "a_sources"),
                                  cache_arg=caches["source_cache"], allow_external_api=False,
                                  verification_policy=admission_policy, focal_species=state["focal_species"])
    sessions = {}
    machine_reviews = {}
    independent_review_routes = []
    for route in ("a", "b"):
        session_inputs = {"records": caches["record_cache"]}
        if route == "a":
            session_inputs.update(hypotheses=a_sources["mechanisms_with_verified_sources"],
                sources=a_sources["source_review_records"], proposal=state["inputs"]["a_payload"])
        else:
            session_inputs.update(hypotheses=imported["proposals"], sources=str(records),
                proposal=raw, bundle=bundle["bundle"], source_verification=sources["source_verification"])
        sessions[route] = prepare_session(str(out / route), route=route, inputs=session_inputs,
            proposer_id=state["identities"][route + "_proposer_id"],
            proposer_model=state["identities"][route + "_model"],
            case_id=state["case_id"], focal_species=state["focal_species"], review_policy=admission_policy)
        if admission_policy == INITIAL_POLICY:
            continue
        prompt = Path(sessions[route]["prompt"])
        source_records = read(prompt).get("source_review_records")
        if not isinstance(source_records, dict) or not isinstance(source_records.get("records"), dict):
            raise ValueError("source record object is invalid; absence must not be guessed")
        if source_records["records"]:
            independent_review_routes.append(route)
        else:
            response = out / (route + "-not-executed.json")
            write(response, build_no_source_review(prompt, route))
            machine_reviews[route] = {"response": response.name,
                "reviewer_id": "no-source-machine-" + route, "reviewer_model": "not_invoked"}
    load_stage(prepared, expected, "awaiting_b_response")
    return seal(out, {"state": "ready_for_initial_candidates" if admission_policy == INITIAL_POLICY else "awaiting_independent_reviews", "case_id": state["case_id"],
        "candidate_admission_policy": admission_policy,
        "prepared": str(prepared), "prepared_sha256": expected, "sessions": sessions,
        "machine_reviews": machine_reviews, "independent_review_routes": independent_review_routes,
        "b_import": imported, "b_sources": sources, "route_b_execution": route_b_execution,
        **({"review_materials": materials, "a_sources": a_sources,
            "a_resolved_input": {"path": state["a_resolved"]["mechanisms_resolved"],
                                 "sha256": digest(Path(state["a_resolved"]["mechanisms_resolved"]))}}
           if supplied or live_b is not None else {}),
        "prediction_calls": 0})


def finish(reviewed: Path, expected: str, reviews_path: Path | None, out: Path,
           route_c_members: Path | None = None, route_c_manifest: Path | None = None,
           route_c_status: Path | None = None, max_adaptive_rounds: int = 1,
           allow_provisional_wet: bool = False,
           identity_assessment: Path | None = None) -> dict:
    from candidate_qtl.claim_review_session import build_no_source_review, complete_session, handoff_sessions

    if isinstance(max_adaptive_rounds, bool) or not isinstance(max_adaptive_rounds, int) or max_adaptive_rounds < 1:
        raise ValueError("max_adaptive_rounds must be a positive integer")
    reviewed = reviewed.resolve()
    if (route_c_members is None) != (route_c_manifest is None):
        raise ValueError("Route C members and producer manifest must be supplied together")
    c_rows = None
    c_members_digest = None
    c_manifest_digest = None
    if route_c_members:
        member_bytes = route_c_members.read_bytes()
        c_members_digest = hashlib.sha256(member_bytes).hexdigest()
        c_rows = sum(1 for _ in csv.DictReader(io.StringIO(member_bytes.decode("utf-8-sig")), delimiter="\t"))
    if route_c_manifest:
        c_manifest_digest = hashlib.sha256(route_c_manifest.read_bytes()).hexdigest()
    c_execution = {"state": "not_evaluated_no_input"}
    c_status_sha = None
    if route_c_status:
        c_status_path = route_c_status.resolve()
        status_bytes = c_status_path.read_bytes()
        c_status_sha = hashlib.sha256(status_bytes).hexdigest()
        c_execution = {
            **json.loads(status_bytes.decode("utf-8-sig")),
            "status_path": str(c_status_path),
            "status_sha256": c_status_sha,
        }
        allowed = ({"completed_members", "completed_empty"} if route_c_members else
                   {"provider_failure", "not_evaluated", "not_attempted"})
        if (c_execution.get("state") not in allowed
                or (route_c_members and c_execution["state"] != ("completed_empty" if c_rows == 0 else "completed_members"))
                or not isinstance(c_execution.get("reason", ""), str)
                or (not route_c_members and not c_execution.get("reason", "").strip())):
            raise ValueError("Route C status does not match supplied artifacts")
    elif route_c_members:
        c_execution = {"state": "completed_empty" if c_rows == 0 else "completed_members"}
    if route_c_members and route_c_manifest:
        c_execution.update({
            "members_path": str(route_c_members.resolve()),
            "members_sha256": c_members_digest,
            "manifest_path": str(route_c_manifest.resolve()),
            "manifest_sha256": c_manifest_digest,
        })
    stage = load_stage(reviewed, expected, ("ready_for_initial_candidates", "awaiting_independent_reviews"))
    if route_c_status and c_execution.get("case_id") != stage["case_id"]:
        raise ValueError("Route C status case_id does not match the frozen case")
    prepared = load_stage(Path(stage["prepared"]), stage["prepared_sha256"], "awaiting_b_response")
    reviews = read(reviews_path) if reviews_path is not None else {}
    if set(reviews) - {"a", "b"}:
        raise ValueError("unknown review route")
    # Explicit answers always win, including reviews from older frozen stages.
    # Only new stages carry sealed machine records for omitted empty-source routes.
    for route, spec in stage.get("machine_reviews", {}).items():
        if route in reviews:
            continue
        if route not in {"a", "b"}:
            raise ValueError("unknown machine review route")
        response = (reviewed.parent / spec["response"]).resolve()
        if (not response.is_relative_to(reviewed.parent)
                or str(response.relative_to(reviewed.parent)) not in stage["files"]
                or spec["reviewer_id"] != "no-source-machine-" + route
                or spec["reviewer_model"] != "not_invoked"
                or read(response) != build_no_source_review(Path(stage["sessions"][route]["prompt"]), route)):
            raise ValueError("machine review does not match the frozen source request")
        reviews[route] = {**spec, "response": str(response)}
    from candidate_qtl.bibliography import INITIAL_POLICY
    initial = stage.get("candidate_admission_policy") == INITIAL_POLICY
    if initial and reviews:
        raise ValueError("initial-candidate policy needs no independent review response")
    if not initial and set(reviews) != {"a", "b"}:
        raise ValueError("independent review responses required for routes: " + ", ".join(sorted({"a", "b"} - set(reviews))))
    # Finish in copies: neither preparation nor the frozen review requests change.
    out = fresh_output(out, [reviewed.parent, Path(stage["prepared"]).parent,
                             *([route_c_status] if route_c_status else [])])
    if route_c_status:
        status_snapshot = out / "route_c_status.json"
        snapshot(route_c_status.resolve(), status_snapshot)
        if digest(status_snapshot) != c_status_sha:
            raise ValueError("Route C status changed while being copied")
        c_execution["status_path"] = str(status_snapshot)
        c_execution["status_sha256"] = digest(status_snapshot)
    completions = {}
    for route in ("a", "b"):
        spec = reviews.get(route, {})
        response_root = reviews_path.resolve().parent if reviews_path is not None else reviewed.parent
        response = (response_root / spec["response"]).resolve() if not initial else None
        session = out / route
        snapshot(reviewed.parent / route, session)
        # Request paths are relative; snapshot hashes stay valid in the copy.
        completions[route] = complete_session(str(session), stage["sessions"][route]["request_sha256"],
            **({"response_arg": str(response), "reviewer_id": spec["reviewer_id"], "reviewer_model": spec["reviewer_model"]} if not initial else {}))
    b_execution = dict(stage.get("route_b_execution") or {"state": "not_recorded"})
    if b_execution.get("state") not in {"not_recorded", ROUTE_B_NOT_RUN_STATE}:
        b_manifest_source = Path(str(b_execution["path"])).resolve()
        b_manifest_snapshot = out / "route_b_execution_manifest.json"
        snapshot(b_manifest_source, b_manifest_snapshot)
        b_execution["path"] = str(b_manifest_snapshot)
        b_execution["sha256"] = digest(b_manifest_snapshot)
    b_root = Path(prepared["inputs"]["packet_b"])
    result = handoff_sessions(str(out / "a"), str(out / "b"), str(out / "handoff"),
        genes_arg=str(b_root / "candidates.tsv"), proteins_arg=str(b_root / "candidate_proteins.fa"),
        atlas_arg=prepared["b0"]["interval_gene_atlas"],
        family_index_arg=prepared["inputs"].get("family_index"),
        reference_function_arg=prepared["b0"].get("reference_function"),
        route_c_members_arg=str(route_c_members.resolve()) if route_c_members else None,
        route_c_manifest_arg=str(route_c_manifest.resolve()) if route_c_manifest else None,
        route_b_execution=b_execution,
        route_c_execution=c_execution,
        allow_unrecorded_route_b=not bool(prepared.get("ordinary_case")),
        allow_provisional_wet=allow_provisional_wet,
        identity_assessment_arg=str(identity_assessment.resolve()) if identity_assessment else None,
        max_adaptive_rounds=max_adaptive_rounds)
    result["route_b_execution"] = result.get("execution_context", {}).get("route_b_execution", b_execution)
    result["route_c_execution"] = result.get("execution_context", {}).get("route_c_execution", c_execution)
    load_stage(reviewed, expected, ("ready_for_initial_candidates", "awaiting_independent_reviews"))
    load_stage(Path(stage["prepared"]), stage["prepared_sha256"], "awaiting_b_response")
    state = (result["wet"]["handoff_state"] if allow_provisional_wet and
             result["wet"]["handoff_state"] == "provisional_adopted_wet" else
             result["state"] if result["wet"].get("wet_candidates") else result["wet"]["handoff_state"])
    report = finish_report(out, stage["case_id"], state, result)
    return seal(out, {"state": state, "case_id": stage["case_id"], "report": report,
        "handoff": result, "reviews": completions, "b_import": stage["b_import"],
        "review_stage": {"path": str(reviewed), "sha256": expected},
        "route_c_input": ({"members": str(route_c_members.resolve()),
                           "members_sha256": c_members_digest,
                           "manifest": str(route_c_manifest.resolve()),
                           "manifest_sha256": c_manifest_digest}
                          if route_c_members and route_c_manifest else None),
        "route_c_execution": result.get("route_c_execution", c_execution),
        "route_b_execution": result.get("route_b_execution", b_execution),
        "prediction_calls": 0})



def audit(finished: Path, expected: str) -> dict:
    """Inspect already-sealed ordinary stages without changing science artifacts.

    This is not a model run, Web/tool capability test, truth check or runtime
    measurement. Missing execution provenance stays not_recorded.
    """
    finished = finished.resolve()
    stage_name = read(finished).get("state")
    if not isinstance(stage_name, str) or not stage_name:
        raise ValueError("finished stage lacks an execution state")
    stage = load_stage(finished, expected, stage_name)
    review_ref = stage.get("review_stage")
    if not isinstance(review_ref, dict):
        raise ValueError("finished stage lacks review_stage provenance")
    reviewed_file = Path(str(review_ref["path"])).resolve()
    reviewed = load_stage(reviewed_file, str(review_ref["sha256"]), ("ready_for_initial_candidates", "awaiting_independent_reviews"))
    prepared_file = Path(str(reviewed["prepared"])).resolve()
    prepared = load_stage(prepared_file, str(reviewed["prepared_sha256"]), "awaiting_b_response")
    case_id = str(stage.get("case_id", ""))
    if not case_id or any(value.get("case_id") != case_id for value in (reviewed, prepared)):
        raise ValueError("audit stage case_id mismatch")

    source_counts = {}
    for route in ("a", "b"):
        session = reviewed.get("sessions", {}).get(route, {})
        prompt_name = session.get("prompt")
        if not prompt_name:
            source_counts[route] = None
            continue
        prompt_file = Path(str(prompt_name)).resolve()
        if not prompt_file.is_relative_to(reviewed_file.parent):
            raise ValueError("review prompt outside sealed stage")
        source_records = read(prompt_file).get("source_review_records", {}).get("records")
        if not isinstance(source_records, dict):
            raise ValueError("review source records are not a dictionary")
        source_counts[route] = len(source_records)

    b_execution = stage.get("route_b_execution") or {"state": "not_recorded"}
    b_state = str(b_execution.get("state", "not_recorded"))
    b_summary = {"state": b_state, "prompt_count": None, "completed_count": None,
                 "not_evaluated_count": None, "attempt_count": None, "retry_count": None}
    if b_state == ROUTE_B_NOT_RUN_STATE:
        b_summary.update(prompt_count=b_execution.get("prompt_count"), completed_count=0,
                         not_evaluated_count=b_execution.get("not_evaluated_count"),
                         attempt_count=0, retry_count=0, reason=b_execution.get("reason", ""))
    elif b_state != "not_recorded":
        b_manifest_file = Path(str(b_execution["path"])).resolve()
        if not b_manifest_file.is_relative_to(finished.parent):
            raise ValueError("B execution manifest outside sealed finish")
        b_expected = str(b_execution.get("sha256", "")).removeprefix("sha256:")
        if not b_expected or digest(b_manifest_file) != b_expected:
            raise ValueError("B execution manifest hash mismatch")
        b_manifest = read(b_manifest_file)
        if b_manifest.get("case_id") != case_id:
            raise ValueError("B execution manifest case_id mismatch")
        prompt_count = b_manifest.get("prompt_count")
        completed = b_manifest.get("completed_count")
        not_evaluated = b_manifest.get("not_evaluated_count")
        if any(type(value) is not int or value < 0 for value in
               (prompt_count, completed, not_evaluated)) or completed + not_evaluated != prompt_count:
            raise ValueError("B execution counts inconsistent")
        attempts = b_manifest.get("partial_response_provenance")
        retries = b_manifest.get("retry_history")
        if not isinstance(attempts, list) or not isinstance(retries, list):
            raise ValueError("B execution retry provenance missing")
        b_summary.update(prompt_count=prompt_count, completed_count=completed,
                         not_evaluated_count=not_evaluated, attempt_count=len(attempts),
                         retry_count=len(retries))

    wet = (stage.get("handoff") or {}).get("wet") or {}
    wet_file_name = wet.get("wet_candidates")
    wet_count = None
    wet_digest = None
    if wet_file_name:
        wet_file = Path(str(wet_file_name)).resolve()
        if not wet_file.is_relative_to(finished.parent):
            raise ValueError("Wet candidates outside sealed finish")
        wet_digest = digest(wet_file)
        wet_expected = str(wet.get("wet_candidates_sha256", "")).removeprefix("sha256:")
        if wet_expected and wet_digest != wet_expected:
            raise ValueError("Wet candidates hash mismatch")
        with wet_file.open(encoding="utf-8-sig", newline="") as stream:
            candidates = list(csv.DictReader(stream, delimiter="\t"))
        identifiers = [str(row.get("candidate_gene_id", "")).strip() for row in candidates]
        if any(not identifier for identifier in identifiers) or len(identifiers) != len(set(identifiers)):
            raise ValueError("Wet candidate IDs must be present and unique")
        wet_count = len(identifiers)

    return {
        "contract": "ordinary-e2e-readonly-audit-1.0",
        "case_id": case_id,
        "stage": {"prepared_sha256": str(reviewed["prepared_sha256"]),
                  "reviewed_sha256": str(review_ref["sha256"]),
                  "finished_sha256": expected, "finished_state": stage_name},
        "source_records_presented_to_review": source_counts,
        "independent_review_required_routes": reviewed.get("independent_review_routes", []),
        "no_source_machine_not_executed_routes": sorted(reviewed.get("machine_reviews", {})),
        "route_b_execution": b_summary,
        "wet": {"handoff_state": wet.get("handoff_state", "not_recorded"),
                "second_stage_state": wet.get("second_stage_state", "not_recorded"),
                "adopted_stage": (stage.get("handoff") or {}).get("adopted_stage", "not_recorded"),
                "candidate_count": wet_count, "candidate_sha256": wet_digest,
                "report": stage.get("report")},
        "not_inferred": ["actual_search_fetch_tool_capabilities", "worker_source_read_receipts",
                         "walltime", "actual_tokens", "api_cost", "truth_retention", "first_loss"],
    }

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--config", type=Path, required=True)
    for name, anchor, response in (("review", "prepared", "responses"), ("finish", "reviewed", "reviews")):
        p = sub.add_parser(name)
        p.add_argument("--" + anchor, type=Path, required=True)
        p.add_argument("--sha256", required=True)
        p.add_argument("--" + response, type=Path, required=name == "review",
                       help="Optional legacy content-review specs; initial candidates need no review responses"
                       if name == "finish" else None)
        if name == "review":
            p.add_argument("--b-execution-manifest", type=Path)
            p.add_argument("--route-b-not-run", metavar="REASON",
                           help="Optional Route B was deliberately not run (A/C-first); every response must be not evaluated")
            p.add_argument("--source-cache", type=Path,
                           help="Explicit post-proposal source cache; snapshot and recheck frozen A/B sources offline")
            p.add_argument("--record-cache", type=Path,
                           help="Explicit post-proposal reference records directory; snapshot for A/B review sessions")
    for p in sub.choices.values():
        p.add_argument("--out-dir", type=Path, required=True)
    sub.choices["finish"].add_argument("--route-c-members", type=Path)
    sub.choices["finish"].add_argument("--route-c-manifest", type=Path)
    sub.choices["finish"].add_argument("--route-c-status", type=Path)
    sub.choices["finish"].add_argument("--max-adaptive-rounds", type=int, default=1)
    sub.choices["finish"].add_argument("--allow-provisional-wet", action="store_true")
    sub.choices["finish"].add_argument("--identity-assessment", type=Path,
                                     help="Hash-bound post-freeze same-locus assessment; uses the prepared target proteins")
    p = sub.add_parser("audit", help="Read-only summary of already-sealed ordinary stages")
    p.add_argument("--finished", type=Path, required=True)
    p.add_argument("--sha256", required=True)
    args = parser.parse_args(argv)
    if args.command == "prepare":
        result = prepare(args.config, args.out_dir)
    elif args.command == "review":
        result = review(args.prepared, args.sha256, args.responses, args.out_dir,
                        args.b_execution_manifest, source_cache=args.source_cache, record_cache=args.record_cache,
                        route_b_not_run=args.route_b_not_run)
    elif args.command == "finish":
        result = finish(args.reviewed, args.sha256, args.reviews, args.out_dir,
                        args.route_c_members, args.route_c_manifest, args.route_c_status,
                        args.max_adaptive_rounds, args.allow_provisional_wet, args.identity_assessment)
    else:
        result = audit(args.finished, args.sha256)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
