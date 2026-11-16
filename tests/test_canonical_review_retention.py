import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl.canonical_route_b_wet import run_canonical_route_b_wet_extension
from tests.test_canonical_route_b_wet import _write_b0, _write_route_a_bridge, _write_llm_manifest
from candidate_qtl.route_b_semantic_bundle import build_route_b_semantic_bundle
from tests._route_b_policy_fixtures import unversioned_bundle


def _canonical_unverified_b_proposal_reaches_wet_review_without_union_support(tmp_path: Path, policy="initial_bibliography_v1"):
    def build_policy_bundle(*args, **kwargs):
        kwargs["b_review_membership_policy"] = "source_function_v1" if policy == "legacy_grounded_review" else policy
        built = build_route_b_semantic_bundle(*args, **kwargs)
        return unversioned_bundle(built) if policy == "legacy_grounded_review" else built

    def run_extension(**kwargs):
        with patch("candidate_qtl.canonical_route_b_wet.build_route_b_semantic_bundle", side_effect=build_policy_bundle):
            return run_canonical_route_b_wet_extension(str(tmp_path), phenotype="fasciated meristem", case_id="case1", **kwargs)

    atlas_path = _write_b0(tmp_path)
    atlas = pd.read_csv(atlas_path, sep="\t", dtype=str).fillna("")
    atlas["orthology_reference_gene"] = ["", "SYNTHETIC_AT_REFERENCE"]
    atlas["reference_species"] = ["", "Arabidopsis thaliana"]
    atlas["local_annotation_description"] = ["", "synthetic meristem signaling kinase annotation"]
    atlas["broad_function_bins"] = ["", "signaling;meristem_development"]
    atlas.to_csv(atlas_path, sep="\t", index=False)
    _write_route_a_bridge(tmp_path)
    prepared = run_extension()
    prompts = [json.loads(line) for line in Path(prepared["prompts"]).read_text(encoding="utf-8").splitlines()]
    responses = Path(tmp_path / "b04_canonical_semantic" / "02_llm_run" / "llm_candidate_responses.json")
    responses.parent.mkdir(parents=True, exist_ok=True)
    candidates = []
    for prompt in prompts:
        claim = {
            "claim_type": "mechanistic_review_seed",
            "evidence_scope": "family",
            "biological_function": "meristem maintenance",
            "mechanism": "Synthetic reference material suggests meristem signaling; homologous family activity could change focal organ initiation",
            "phenotype_relation": "cross-species fasciation hypothesis, not an experiment in the focal species; trait transfer is unverified",
            "reference_species": "Arabidopsis thaliana",
            "reference_gene": "SYNTHETIC_AT_REFERENCE" if "SYNTHETIC_AT_REFERENCE" in prompt["user"] else "",
            "source_ids": [],
            "uncertainty_note": "review seed only; no direct phenotype evidence was supplied",
        }
        candidates.append({"prompt_id": prompt["prompt_id"], "parse_status": "json",
                           "model_json": {"claims": [claim]}})
    responses.write_text(json.dumps({"candidates": candidates}), encoding="utf-8")
    _write_llm_manifest(prepared["bundle"], prepared["prompts"], responses,
                        responses.parent / "llm_run_manifest.json")
    result = run_extension(wet_budget=2)
    assert result["state"] in {"wet_reviewable", "wet_reviewable_degraded_coverage"}
    wet = pd.read_csv(result["wet_candidates"], sep="\t", dtype=str).fillna("")
    review = wet[wet["membership_basis"] == "route_review_candidate"]
    extension = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
    legacy = policy == "legacy_grounded_review"
    assert set(review.candidate_gene_id) == ({"G_B"} if legacy else set())
    assert set(wet["candidate_gene_id"]) == ({"G_A", "G_B"} if legacy else {"G_A"})
    b_review = pd.read_csv(extension["route_b_review_candidates"], sep="\t", dtype=str).fillna("")
    assert set(b_review.candidate_gene_id) == ({"G_B"} if legacy else set())
    b_review_manifest = json.loads(Path(extension["route_b_review_candidates_manifest"]).read_text(encoding="utf-8"))
    assert b_review_manifest["review_membership_policy"] == policy
    assert b_review_manifest["candidate_gene_count"] == int(legacy)
    ledger = pd.read_csv(b_review_manifest["eligibility_ledger"]["path"], sep="\t", dtype=str)
    assert "G_B" in set(ledger.candidate_gene_id)
    expected = "source_absent" if policy == "source_function_v1" else "not_evaluated"
    assert set(ledger.loc[ledger.candidate_gene_id.eq("G_B"), "b_review_eligibility"]) == {expected}
    assert set(ledger.biological_exclusion) == {"not_asserted"}
    union = pd.read_csv(extension["pre_family_supported_union"], sep="\t", dtype=str).fillna("")
    assert "G_B" not in set(union["candidate_gene_id"])


class CanonicalReviewRetentionTests(unittest.TestCase):
    def test_canonical_unverified_b_proposal_reaches_wet_review_without_union_support(self):
        for policy in ("initial_bibliography_v1", "source_function_v1", "legacy_grounded_review"):
            with self.subTest(policy=policy), tempfile.TemporaryDirectory() as directory:
                _canonical_unverified_b_proposal_reaches_wet_review_without_union_support(Path(directory), policy)
