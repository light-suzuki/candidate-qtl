"""Fixed-response tests of the saved Sol reader; no scientific execution."""
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from scripts.agentic_route_b import merge, prepare
from tests.test_agentic_route_b_material import semantic_prompt

REPO = Path(__file__).resolve().parents[1]
READER = REPO / "workflow/bounded_b_reader.py"


def cli(state, command, *args, ok=True):
    env = {**os.environ, "PYTHONPATH": os.pathsep.join([str(REPO), str(REPO / "src"), os.environ.get("PYTHONPATH", "")])}
    result = subprocess.run([sys.executable, "-X", "utf8", str(READER), command, "--state", str(state), *map(str, args)],
                            capture_output=True, text=True, encoding="utf-8", env=env)
    if not ok:
        assert result.returncode != 0
        return result.stderr
    assert result.returncode == 0, result.stderr
    assert len(result.stdout) <= 8000
    return json.loads(result.stdout)


@pytest.fixture
def reading(tmp_path):
    rows = [semantic_prompt(i, function='complete 漢🌱 "\\\n' * 400) for i in range(3)]
    source = tmp_path / "source.jsonl"
    source.write_bytes((''.join(json.dumps(row, ensure_ascii=False) + '\n' for row in rows)).encode())
    manifest = prepare(str(source), str(tmp_path / "view"), max_prompts=1)
    return rows, source, manifest


def init(tmp_path, manifest, name="cursor", *extra):
    state = tmp_path / name
    result = cli(state, "init", "--view-root", tmp_path / "view", "--start", 0, "--stop", 3,
                 "--expected-source-sha256", manifest["source_sha256"], "--executor", "test-worker",
                 "--model", "TEST/fixed", *extra)
    return state, result


def consume(state):
    chunks = []
    result = cli(state, "next")
    while result["state"] == "READ_ONE_CHUNK":
        chunks.append(result)
        result = cli(state, "next", "--ack", result["ack_token"])
    return result, chunks


def commit(state, row, tmp_path, raw_text='{"claims":[]}'):
    response = {"prompt_id": row["prompt_id"], "route_b_unit_id": row["route_b_unit_id"],
                "status": "completed", "raw_text": raw_text}
    decision = tmp_path / (row["prompt_id"] + "-fixed.json")
    decision.write_bytes(json.dumps(response).encode())
    cli(state, "commit", "--executor", "test-worker", "--decision-file", decision,
        "--attest-genuine-authorship", "--authorship-note", "Synthetic fixed response for transport test")
    return response


def test_shared_full_material_read_once_and_each_binding_preserved(tmp_path, reading):
    rows, source, manifest = reading
    before = source.read_bytes()
    state, _ = init(tmp_path, manifest)
    chunks = []
    for row in rows:
        status, delivered = consume(state)
        assert status["state"] == "unit_ready_requires_genuine_decision"
        assert status["prompt_id"] == row["prompt_id"]
        chunks.extend(delivered)
        commit(state, row, tmp_path)
    # Every shared raw artifact is delivered once; unit bindings remain separate.
    shared = [c for c in chunks if c["kind"] == "complete_shared_material"]
    groups = {}
    for chunk in shared:
        groups.setdefault(chunk["raw_artifact_sha256"], []).append(chunk)
    for parts in groups.values():
        assert len({p["source_index"] for p in parts}) == 1
        assert parts[0]["begin_character"] == 0
        assert parts[-1]["end_character"] == parts[-1]["total_characters"]
        assert all(a["end_character"] == b["begin_character"] for a, b in zip(parts, parts[1:]))
    wrappers = []
    for index in range(3):
        parts = [c for c in chunks if c["kind"] == "complete_unit_tables_bindings" and c["source_index"] == index]
        wrappers.append(json.loads(''.join(p["content"] for p in parts)))
    assert [w["prompts"][0]["prompt"]["prompt_id"] for w in wrappers] == [p["prompt_id"] for p in rows]
    assert len(groups) < sum(len(w["materials"]) for w in wrappers)
    assert source.read_bytes() == before
    out = tmp_path / "complete.json"
    cli(state, "export", "--out", out)
    exported = json.loads(out.read_bytes())
    assert exported["model"] == "TEST/fixed"
    assert exported["reading_provenance"]["genuine_decision_count"] == 3


def test_pending_repeat_wrong_ack_and_unread_commit_are_guarded(tmp_path, reading):
    rows, _, manifest = reading
    state, _ = init(tmp_path, manifest)
    pending = cli(state, "next")
    assert cli(state, "next") == pending
    assert "ack token" in cli(state, "next", "--ack", "wrong", ok=False)
    assert "unread" in cli(state, "commit", ok=False)
    assert cli(state, "next") == pending
    assert "all genuine" in cli(state, "export", "--out", tmp_path / "bad.json", ok=False)


def test_partial_export_and_checkpoint_retry_only_missing_ids(tmp_path, reading):
    rows, source, manifest = reading
    state, _ = init(tmp_path, manifest)
    consume(state)
    fixed = json.dumps({"claims": [{"claim_type": "mechanistic_review_seed", "biological_function": "Fixed function",
        "mechanism": "Conditional process", "phenotype_relation": "Fixed conditional bridge", "reference_gene": "ref-1",
        "reference_species": "Species one", "evidence_scope": "reference_gene", "source_ids": [],
        "alternative_hypothesis": "Other process", "uncertainty_note": "Unverified transfer"}]})
    healthy = commit(state, rows[0], tmp_path, fixed)
    cli(state, "next")  # unacknowledged partial input is not completion
    partial = tmp_path / "partial.json"
    cli(state, "export", "--allow-partial", "--out", partial)
    data = json.loads(partial.read_bytes())
    assert data["responses"][0] == healthy
    assert all(r["status"] == "not_evaluated" and r["reason"] == "reading_incomplete" for r in data["responses"][1:])
    assert all("raw_text" not in r for r in data["responses"][1:])
    checkpoint = tmp_path / "merged.json"
    first = merge(str(source), [str(partial)], str(checkpoint), case_id="case-x")
    assert first["completed_count"] == 1 and first["not_evaluated_count"] == 2
    checkpoint_before = checkpoint.read_bytes()
    retry, result = init(tmp_path, manifest, "retry", "--checkpoint", checkpoint)
    assert result["units"] == 2
    retry_data = json.loads((retry / "state.json").read_bytes())
    assert [u["source_index"] for u in retry_data["units"]] == [1, 2]
    assert retry_data["shared_fully_acknowledged"] == {}  # new context rereads sources, never another worker's ACK
    for row in rows[1:]:
        consume(retry)
        commit(retry, row, tmp_path)
    retry_export = tmp_path / "retry.json"
    cli(retry, "export", "--out", retry_export)
    result = merge(str(source), [str(checkpoint), str(retry_export)], str(tmp_path / "complete.json"), case_id="case-x")
    assert result["completed_count"] == 3 and result["not_evaluated_count"] == 0
    assert json.loads((tmp_path / "complete.json").read_bytes())["candidates"][0]["raw_text"] == healthy["raw_text"]
    assert len(json.loads((tmp_path / "complete.json").read_bytes())["candidates"][0]["claims"]) == 1
    assert checkpoint.read_bytes() == checkpoint_before


def test_checkpoint_tamper_is_not_a_retry_shortcut(tmp_path, reading):
    _, source, manifest = reading
    attempt = tmp_path / "attempt.json"
    attempt.write_bytes(json.dumps({"case_id": "case-x", "executor": "TEST", "model": "TEST", "responses": []}).encode())
    checkpoint = tmp_path / "merged.json"
    merge(str(source), [str(attempt)], str(checkpoint), case_id="case-x")
    checkpoint.write_bytes(checkpoint.read_bytes() + b' ')
    state = tmp_path / "refused"
    error = cli(state, "init", "--view-root", tmp_path / "view", "--start", 0, "--stop", 3,
                "--expected-source-sha256", manifest["source_sha256"], "--executor", "TEST", "--checkpoint", checkpoint, ok=False)
    assert "hash mismatch" in error
    assert not (state / "state.json").exists()


def test_validation_cursor_cannot_generate_or_export_science(tmp_path, reading):
    _, _, manifest = reading
    state, _ = init(tmp_path, manifest, "validation", "--validation-only")
    assert "validation cursor" in cli(state, "export", "--allow-partial", "--out", tmp_path / "bad.json", ok=False)


@pytest.fixture
def affinity_reading(tmp_path):
    rows = [semantic_prompt(i, function=("Function alpha " if i % 2 else "Function beta ") * 50,
                            source="PMID:11" if i % 2 else "PMID:22") for i in range(1, 5)]
    source = tmp_path / "source.jsonl"
    source.write_text(''.join(json.dumps(row, ensure_ascii=False) + '\n' for row in rows), encoding="utf-8")
    baseline = prepare(str(source), str(tmp_path / "baseline"), max_prompts=2)
    manifest = prepare(str(source), str(tmp_path / "view"), max_prompts=2, prefer_shared_material=True)
    return rows, source, manifest, baseline


def init_shard(tmp_path, manifest, name, shard, *extra, view="view", ok=True):
    state = tmp_path / name
    result = cli(state, "init", "--view-root", tmp_path / view, "--reading-shard", shard,
                 "--expected-source-sha256", manifest["source_sha256"], "--executor", "test-worker",
                 "--model", "TEST/fixed", *extra, ok=ok)
    return state, result


def test_affinity_assignments_reach_actual_reader_without_repeated_material(tmp_path, affinity_reading):
    rows, source, manifest, baseline = affinity_reading
    by_id = {row["prompt_id"]: row for row in rows}
    assert [sh["prompt_ids"] for sh in manifest["reading_view"]["shards"]] == [
        ["prompt-1", "prompt-3"], ["prompt-2", "prompt-4"]]
    before = source.read_bytes()
    delivered_bytes = {}
    for view, prepared in (("baseline", baseline), ("view", manifest)):
        chunks, exports = [], []
        for index, shard in enumerate(prepared["reading_view"]["shards"]):
            state, result = init_shard(tmp_path, prepared, f"{view}-worker-{index}", shard["file"], view=view)
            assert result["reading_shard"] == shard["file"] and "range" not in result
            recorded = json.loads((state / "state.json").read_bytes())
            assert [unit["prompt_id"] for unit in recorded["units"]] == shard["prompt_ids"]
            assert [unit["source_index"] for unit in recorded["units"]] == [
                rows.index(by_id[pid]) for pid in shard["prompt_ids"]]
            worker_chunks = []
            for pid in shard["prompt_ids"]:
                status, delivered = consume(state)
                assert status["prompt_id"] == pid
                worker_chunks.extend(delivered)
                commit(state, by_id[pid], tmp_path)
            # Full shared materials are ACKed only once in this worker context.
            shared = [chunk for chunk in worker_chunks if chunk["kind"] == "complete_shared_material"]
            for digest in {chunk["raw_artifact_sha256"] for chunk in shared}:
                assert len({chunk["source_index"] for chunk in shared if chunk["raw_artifact_sha256"] == digest}) == 1
            chunks.extend(worker_chunks)
            output = tmp_path / f"{view}-worker-{index}.json"
            cli(state, "export", "--out", output)
            provenance = json.loads(output.read_bytes())["reading_provenance"]
            assert provenance["assigned_reading_shard"] == {
                "file": shard["file"], "sha256": shard["sha256"], "prompt_ids": shard["prompt_ids"]}
            assert "assigned_range" not in provenance
            exports.append(str(output))
        merged = tmp_path / f"{view}-merged.json"
        joined = merge(str(source), exports, str(merged), case_id="case-x")
        assert joined["completed_count"] == 4 and joined["not_evaluated_count"] == 0
        assert [row["prompt_id"] for row in json.loads(merged.read_bytes())["candidates"]] == list(by_id)
        delivered_bytes[view] = sum(len(chunk["content"].encode("utf-8")) for chunk in chunks)
    assert delivered_bytes["view"] < delivered_bytes["baseline"]
    assert source.read_bytes() == before


def test_range_stays_canonical_and_assignment_order_keeps_original_indices(tmp_path, affinity_reading):
    rows, source, manifest, _ = affinity_reading
    old = tmp_path / "range-cursor"
    result = cli(old, "init", "--view-root", tmp_path / "view", "--start", 0, "--stop", 2,
                 "--expected-source-sha256", manifest["source_sha256"], "--executor", "test-worker", "--model", "TEST/fixed")
    assert result["range"] == [0, 2] and "reading_shard" not in result
    data = json.loads((old / "state.json").read_bytes())
    assert "reading_shard" not in data
    assert [unit["prompt_id"] for unit in data["units"]] == [rows[0]["prompt_id"], rows[1]["prompt_id"]]
    # The old schema resumes without any new selection field.
    pending = cli(old, "next")
    assert cli(old, "next") == pending
    output = tmp_path / "range-partial.json"
    cli(old, "export", "--out", output, "--allow-partial")
    assert json.loads(output.read_bytes())["reading_provenance"]["assigned_range"] == [0, 2]

    ordered_rows = [semantic_prompt(i, function=f"Function group {i % 3} " * 50, source=f"PMID:{i % 3}")
                    for i in range(1, 9)]
    ordered_source = tmp_path / "order-source.jsonl"
    ordered_source.write_text(''.join(json.dumps(row) + '\n' for row in ordered_rows), encoding="utf-8")
    ordered = prepare(str(ordered_source), str(tmp_path / "ordered"), max_prompts=4, prefer_shared_material=True)
    shard = ordered["reading_view"]["shards"][0]
    assert shard["prompt_ids"] == ["prompt-1", "prompt-4", "prompt-7", "prompt-2"]
    state, _ = init_shard(tmp_path, ordered, "noncanonical-order", shard["file"], "--validation-only", view="ordered")
    data = json.loads((state / "state.json").read_bytes())
    assert [unit["prompt_id"] for unit in data["units"]] == shard["prompt_ids"]
    assert [unit["source_index"] for unit in data["units"]] == [0, 3, 6, 1]


def test_affinity_partial_retry_keeps_only_unfinished_assigned_ids(tmp_path, affinity_reading):
    rows, source, manifest, _ = affinity_reading
    shard = manifest["reading_view"]["shards"][0]
    state, _ = init_shard(tmp_path, manifest, "first", shard["file"])
    consume(state)
    healthy = commit(state, rows[0], tmp_path)
    consume(state)  # Read prompt-3 fully, but do not invent its independent decision.
    partial = tmp_path / "partial.json"
    cli(state, "export", "--allow-partial", "--out", partial)
    data = json.loads(partial.read_bytes())
    assert data["responses"] == [healthy, {"prompt_id": "prompt-3", "route_b_unit_id": "unit-3",
                                        "status": "not_evaluated", "reason": "awaiting_independent_interpretation"}]
    checkpoint = tmp_path / "checkpoint.json"
    merge(str(source), [str(partial)], str(checkpoint), case_id="case-x")
    before = checkpoint.read_bytes()
    retry, result = init_shard(tmp_path, manifest, "retry-shard", shard["file"], "--checkpoint", checkpoint)
    assert result["units"] == 1
    recorded = json.loads((retry / "state.json").read_bytes())
    assert [(unit["source_index"], unit["prompt_id"]) for unit in recorded["units"]] == [(2, "prompt-3")]
    assert recorded["shared_fully_acknowledged"] == {}
    assert recorded["reading_shard"]["prompt_ids"] == shard["prompt_ids"]
    consume(retry)
    commit(retry, rows[2], tmp_path)
    retry_output = tmp_path / "retry.json"
    cli(retry, "export", "--out", retry_output)
    merged = tmp_path / "merged.json"
    result = merge(str(source), [str(checkpoint), str(retry_output)], str(merged), case_id="case-x")
    assert result["completed_count"] == 2 and result["not_evaluated_count"] == 2
    assert json.loads(merged.read_bytes())["candidates"][0]["raw_text"] == healthy["raw_text"]
    assert checkpoint.read_bytes() == before
    _, result = init_shard(tmp_path, manifest, "already-complete-shard", shard["file"], "--checkpoint", merged)
    assert result["units"] == 0
    # A healthy checkpoint must not waive the selected artifact hash check.
    path = tmp_path / "view" / shard["file"]
    path.write_bytes(path.read_bytes() + b' ')
    refused, error = init_shard(tmp_path, manifest, "tampered-complete-shard", shard["file"], "--checkpoint", merged, ok=False)
    assert "assignment bytes changed" in error
    assert not (refused / "state.json").exists()


@pytest.mark.parametrize("selection,extra,error", [
    ("reading/missing/assignment.json", (), "unknown or duplicate reading shard"),
    ("../assignment.json", (), "unknown or duplicate reading shard"),
    ("reading/b0000/assignment.json", ("--start", 0, "--stop", 2), "mutually exclusive"),
    ("reading/b0000/assignment.json", ("--stop", 2), "mutually exclusive"),
])
def test_reading_shard_selection_rejects_unknown_or_mixed_scope(tmp_path, affinity_reading, selection, extra, error):
    _, _, manifest, _ = affinity_reading
    state, observed = init_shard(tmp_path, manifest, "refused", selection, *extra, ok=False)
    assert error in observed
    assert not (state / "state.json").exists()


@pytest.mark.parametrize("roster", [["prompt-1", "prompt-1"], ["prompt-1", "unknown"], ["prompt-1"]])
def test_reading_shard_roster_drift_cannot_change_assigned_work(tmp_path, affinity_reading, roster):
    _, _, manifest, _ = affinity_reading
    path = tmp_path / "view" / "prepare_manifest.json"
    altered = json.loads(path.read_bytes())
    altered["reading_view"]["shards"][0]["prompt_ids"] = roster
    path.write_text(json.dumps(altered), encoding="utf-8")
    state, error = init_shard(tmp_path, manifest, "refused-roster", "reading/b0000/assignment.json", ok=False)
    assert "roster" in error
    assert not (state / "state.json").exists()
