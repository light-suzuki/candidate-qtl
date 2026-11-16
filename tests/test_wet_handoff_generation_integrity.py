"""A rejected Wet rerun must leave the original frozen generation usable."""

import hashlib
from pathlib import Path

import pytest

from candidate_qtl.boundary_tier_review import _load_scope
from candidate_qtl.cli import main
from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from tests.test_wet_relevance_selection import handoff


def _snapshot(root):
    return {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in root.rglob("*") if path.is_file()}


def _arguments(root):
    return {
        "route_b_candidates_arg": str(root / "supported.tsv"),
        "route_b_candidates_manifest_arg": str(root / "supported_manifest.json"),
        "interval_atlas_arg": str(root / "atlas.tsv"),
        "review_candidates_arg": str(root / "review.tsv"),
        "review_candidates_manifest_arg": str(root / "review_manifest.json"),
        "case_id": "case1",
    }


@pytest.mark.parametrize("changed_inputs", [False, True])
def test_rejected_rerun_preserves_all_original_artifacts_and_scope(tmp_path, changed_inputs):
    first_root = tmp_path / "first"
    original = handoff(first_root, (0, 3, 30, 0, 0))
    input_root = first_root
    if changed_inputs:
        input_root = tmp_path / "second"
        handoff(input_root, (0, 5, 30, 0, 0))
    output = first_root / "wet"
    before = _snapshot(output)

    with pytest.raises(ValueError, match="new or empty"):
        wet_candidate_handoff(str(output), **_arguments(input_root))

    assert _snapshot(output) == before
    scope, _, _ = _load_scope(original["boundary_scope_manifest"], original["boundary_scope_sha256"], "case1")
    assert scope["counts"] == {"frozen": 3, "boundary": 30, "lower": 0}


def test_cli_rejects_reusing_a_generation_before_any_write(tmp_path):
    first = handoff(tmp_path / "first", (0, 3, 0, 0, 0))
    second_root = tmp_path / "second"
    handoff(second_root, (0, 5, 0, 0, 0))
    output = tmp_path / "first" / "wet"
    before = _snapshot(output)
    command = ["wet-candidate-handoff", "--out-dir", str(output)]
    for key, value in _arguments(second_root).items():
        command.extend(["--" + key.removesuffix("_arg").replace("_", "-"), value])

    with pytest.raises(ValueError, match="preserve the previous generation"):
        main(command)

    assert _snapshot(output) == before
    _load_scope(first["boundary_scope_manifest"], first["boundary_scope_sha256"], "case1")


def test_an_existing_empty_output_directory_accepts_a_new_generation(tmp_path):
    (tmp_path / "wet").mkdir()
    result = handoff(tmp_path, (0, 3, 0, 0, 0))
    assert Path(result["wet_candidates"]).is_file()
    _load_scope(result["boundary_scope_manifest"], result["boundary_scope_sha256"], "case1")


def test_partial_generation_is_preserved_for_diagnosis(tmp_path):
    output = tmp_path / "wet"
    output.mkdir()
    partial = output / "supported_member_union"
    partial.mkdir()
    (partial / "checkpoint.tsv").write_text("candidate_gene_id\nG1\n", encoding="utf-8")
    before = _snapshot(output)
    with pytest.raises(ValueError, match="new or empty"):
        wet_candidate_handoff(str(output), case_id="case1")
    assert _snapshot(output) == before
