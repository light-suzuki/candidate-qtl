"""Ordinary Wet membership must use the supplied complete interval roster."""

import pandas as pd
import pytest

from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from tests._issue664_helpers import manifest_for


def _inputs(root, atlas_ids):
    candidates = root / "supported.tsv"
    pd.DataFrame([{"candidate_gene_id": "G_supported", "route": "b",
                   "support_status": "verified_semantic_support"}]).to_csv(
                       candidates, sep="\t", index=False)
    manifest = manifest_for(candidates, route="b", case_id="case1")
    atlas = root / "atlas.tsv"
    pd.DataFrame([{"gene_id": gene, "case_id": "case1", "family_id": "F"}
                  for gene in atlas_ids], columns=["gene_id", "case_id", "family_id"]).to_csv(
                      atlas, sep="\t", index=False)
    return {
        "route_b_candidates_arg": str(candidates),
        "route_b_candidates_manifest_arg": str(manifest),
        "interval_atlas_arg": str(atlas),
        "case_id": "case1",
    }


@pytest.mark.parametrize("atlas_ids", [["G_other_interval"], []])
def test_same_case_cannot_admit_support_outside_current_atlas(tmp_path, atlas_ids):
    arguments = _inputs(tmp_path, atlas_ids)
    output = tmp_path / "wet"
    with pytest.raises(ValueError, match="outside the supplied interval atlas: G_supported"):
        wet_candidate_handoff(str(output), **arguments)
    assert not (output / "wet_candidates.tsv").exists()
    assert not (output / "boundary_scope").exists()


@pytest.mark.parametrize("atlas_ids", [["G_supported", "G_other"], ["G_other", "G_supported"]])
def test_supported_candidate_in_current_atlas_is_retained(tmp_path, atlas_ids):
    result = wet_candidate_handoff(str(tmp_path / "wet"), **_inputs(tmp_path, atlas_ids))
    wet = pd.read_csv(result["wet_candidates"], sep="\t")
    assert wet.candidate_gene_id.tolist() == ["G_supported"]
    assert wet.support_tier.tolist() == ["tier1_verified_support"]


def test_legacy_supported_only_does_not_require_a_complete_fill_atlas(tmp_path):
    result = wet_candidate_handoff(str(tmp_path / "wet"),
                                   **_inputs(tmp_path, ["G_other_interval"]),
                                   wet_fill_mode="supported_only")
    assert pd.read_csv(result["wet_candidates"], sep="\t").candidate_gene_id.tolist() == ["G_supported"]
