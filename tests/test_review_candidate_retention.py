from pathlib import Path
import hashlib
import tempfile
import unittest
import json

import pandas as pd

from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from tests._issue664_helpers import manifest_for


def _sealed_review(path: Path, case_id: str, atlas: Path) -> Path:
    pd.DataFrame([{
        "candidate_gene_id": "G2", "route": "b",
        "support_status": "unverified_route_b_review", "rank_contribution": "none",
        "verification_status": "not_evaluated", "verification_reason": "quote incomplete",
    }]).to_csv(path, sep="\t", index=False, lineterminator="\n")
    manifest = {
        "schema_version": "1.0", "contract": "route-review-candidates-1.0.0",
        "stage": "route_b_review_candidate_freeze", "case_id": case_id,
        "output": str(path), "artifact_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "inputs": {"atlas": {"path": str(atlas), "sha256": hashlib.sha256(atlas.read_bytes()).hexdigest()}},
    }
    out = path.with_name("review_manifest.json")
    out.write_text(json.dumps(manifest), encoding="utf-8")
    return out


def _route_review_is_retained_without_generic_interval_fill(tmp_path: Path):
    supported = tmp_path / "supported.tsv"
    pd.DataFrame({"candidate_gene_id": ["G1"], "route": ["b"],
                  "support_status": ["verified_semantic_support"]}).to_csv(
        supported, sep="\t", index=False, lineterminator="\n")
    atlas = tmp_path / "atlas.tsv"
    pd.DataFrame({"gene_id": ["G1", "G2", "G3"], "family_id": ["F1", "F2", "F3"]}).to_csv(
        atlas, sep="\t", index=False, lineterminator="\n")
    review = tmp_path / "review.tsv"
    review_manifest = _sealed_review(review, "case1", atlas)
    result = wet_candidate_handoff(
        str(tmp_path / "out"), route_b_candidates_arg=str(supported),
        route_b_candidates_manifest_arg=str(manifest_for(supported, route="b", case_id="case1")),
        review_candidates_arg=str(review), review_candidates_manifest_arg=str(review_manifest),
        interval_atlas_arg=str(atlas), case_id="case1", wet_budget=3,
    )
    wet = pd.read_csv(result["wet_candidates"], sep="\t", dtype=str).fillna("")
    assert list(wet["candidate_gene_id"]) == ["G1", "G2"]
    assert wet.loc[wet.candidate_gene_id == "G2", "support_tier"].item() == "unverified_route_review"
    assert wet.loc[wet.candidate_gene_id == "G2", "membership_basis"].item() == "route_review_candidate"
    assert {"review_routes", "review_claim_ids", "review_verification_statuses", "review_reasons", "review_source_ids"} <= set(wet.columns)
    assert wet.loc[wet.candidate_gene_id == "G2", "review_verification_statuses"].item() == "not_evaluated"
    assert wet.loc[wet.candidate_gene_id == "G2", "review_reasons"].item() == "quote incomplete"
    wet_manifest = json.loads(Path(result["wet_candidates_manifest"]).read_text(encoding="utf-8"))
    trace = pd.read_csv(wet_manifest["candidate_trace"], sep="\t", dtype=str).fillna("")
    g2 = trace.set_index("candidate_gene_id").loc["G2"]
    g3 = trace.set_index("candidate_gene_id").loc["G3"]
    assert g2["explicit_review_member"] == "True"
    assert g2["adopted_first_stage"] == "True"
    assert g3["explicit_review_member"] == "False"
    assert g3["adopted_first_stage"] == "False"
    assert g3["selection_disposition"] == "not_selected_no_upstream_candidate"


def _supported_only_diagnostic_excludes_unverified_route_review(tmp_path: Path):
    supported = tmp_path / "supported.tsv"
    pd.DataFrame({"candidate_gene_id": ["G1"], "route": ["b"],
                  "support_status": ["verified_semantic_support"]}).to_csv(supported, sep="\t", index=False)
    atlas = tmp_path / "atlas.tsv"
    pd.DataFrame({"gene_id": ["G1", "G2"], "family_id": ["F1", "F2"]}).to_csv(atlas, sep="\t", index=False)
    review = tmp_path / "review.tsv"
    review_manifest = _sealed_review(review, "case1", atlas)
    result = wet_candidate_handoff(
        str(tmp_path / "out"), route_b_candidates_arg=str(supported),
        route_b_candidates_manifest_arg=str(manifest_for(supported, route="b", case_id="case1")),
        review_candidates_arg=str(review), review_candidates_manifest_arg=str(review_manifest),
        interval_atlas_arg=str(atlas), case_id="case1", wet_budget=2, wet_fill_mode="supported_only",
    )
    wet = pd.read_csv(result["wet_candidates"], sep="\t", dtype=str).fillna("")
    assert set(wet["candidate_gene_id"]) == {"G1"}


def _combined_routes_merge_disjoint_review_provenance(tmp_path: Path):
    supported = tmp_path / "supported.tsv"
    pd.DataFrame({"candidate_gene_id": ["G1"], "route": ["b"], "support_status": ["verified_semantic_support"]}).to_csv(supported, sep="\t", index=False)
    atlas = tmp_path / "atlas.tsv"
    pd.DataFrame({"gene_id": ["G1", "G2"], "family_id": ["F1", "F2"]}).to_csv(atlas, sep="\t", index=False)
    b_review = tmp_path / "b_review.tsv"
    b_manifest = _sealed_review(b_review, "case1", atlas)
    b_manifest = b_manifest.rename(tmp_path / "b_review_manifest.json")
    a_review = tmp_path / "a_review.tsv"
    a_manifest = _sealed_review(a_review, "case1", atlas)
    a_manifest = a_manifest.rename(tmp_path / "a_review_manifest.json")
    pd.DataFrame([{"candidate_gene_id": "G2", "route": "a", "support_status": "review_required",
                   "rank_contribution": "none", "review_status": "review_required", "reason": "identity pending",
                   "hypothesis_id": "H1"}]).to_csv(a_review, sep="\t", index=False)
    payload = json.loads(a_manifest.read_text(encoding="utf-8"))
    payload["artifact_sha256"] = hashlib.sha256(a_review.read_bytes()).hexdigest()
    payload["stage"] = "route_a_review_candidate_freeze"
    a_manifest.write_text(json.dumps(payload), encoding="utf-8")
    result = wet_candidate_handoff(
        str(tmp_path / "out"), route_b_candidates_arg=str(supported),
        route_b_candidates_manifest_arg=str(manifest_for(supported, route="b", case_id="case1")),
        review_candidates_arg=[str(a_review), str(b_review)],
        review_candidates_manifest_arg=[str(a_manifest), str(b_manifest)],
        interval_atlas_arg=str(atlas), case_id="case1", wet_budget=2,
    )
    wet = pd.read_csv(result["wet_candidates"], sep="\t", dtype=str).fillna("")
    row = wet.loc[wet.candidate_gene_id == "G2"].iloc[0]
    assert "review_required" in row.review_verification_statuses
    assert "not_evaluated" in row.review_verification_statuses
    assert "identity pending" in row.review_reasons
    assert "quote incomplete" in row.review_reasons


def _strict_capacity_counts_route_review_candidates(tmp_path: Path):
    supported = tmp_path / "supported.tsv"
    pd.DataFrame({"candidate_gene_id": ["G1"], "route": ["b"], "support_status": ["verified_semantic_support"]}).to_csv(supported, sep="\t", index=False)
    atlas = tmp_path / "atlas.tsv"
    pd.DataFrame({"gene_id": ["G1", "G2"], "family_id": ["F1", "F2"]}).to_csv(atlas, sep="\t", index=False)
    review = tmp_path / "review.tsv"
    review_manifest = _sealed_review(review, "case1", atlas)
    result = wet_candidate_handoff(
        str(tmp_path / "out"), route_b_candidates_arg=str(supported),
        route_b_candidates_manifest_arg=str(manifest_for(supported, route="b", case_id="case1")),
        review_candidates_arg=str(review), review_candidates_manifest_arg=str(review_manifest),
        interval_atlas_arg=str(atlas), case_id="case1", wet_budget=1, budget_policy="strict",
    )
    assert result["handoff_state"] == "review_required_over_budget"
    assert result["wet_candidates"] is None


class ReviewCandidateRetentionTests(unittest.TestCase):
    def _run(self, function):
        with tempfile.TemporaryDirectory() as directory:
            function(Path(directory))

    def test_route_review_is_retained_without_generic_interval_fill(self):
        self._run(_route_review_is_retained_without_generic_interval_fill)

    def test_supported_only_diagnostic_excludes_unverified_route_review(self):
        self._run(_supported_only_diagnostic_excludes_unverified_route_review)

    def test_combined_routes_merge_disjoint_review_provenance(self):
        self._run(_combined_routes_merge_disjoint_review_provenance)

    def test_strict_capacity_counts_route_review_candidates(self):
        self._run(_strict_capacity_counts_route_review_candidates)
