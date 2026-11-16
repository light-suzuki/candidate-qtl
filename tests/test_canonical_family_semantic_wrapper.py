"""Canonical wrapper provenance and non-blocking semantics (#505/#514)."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from candidate_qtl.canonical_family_semantic_wet import (
    run_canonical_route_b_family_semantic_wet_extension,
)


def _wet(path: Path, gene: str = "G1") -> str:
    pd.DataFrame([{"candidate_gene_id": gene, "support_status": "supported"}]).to_csv(
        path, sep="\t", index=False, lineterminator="\n"
    )
    return str(path)


class CanonicalFamilySemanticWrapperTests(unittest.TestCase):
    def test_noncanonical_family_copy_is_not_consumed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            wet = _wet(root / "wet.tsv")
            base_manifest = root / "base_manifest.json"
            base_manifest.write_text(json.dumps({
                "family_copy_producer_class": "legacy_triple_route_family_copy",
                "family_recovery_source": "triple_route_compare",
                "legacy_triple_route_family_context_consumed": True,
            }), encoding="utf-8")
            base = {
                "state": "wet_reviewable",
                "manifest": str(base_manifest),
                "wet_candidates": wet,
                "wet_handoff_manifest": str(root / "legacy_handoff.json"),
            }
            with mock.patch(
                "candidate_qtl.canonical_family_semantic_wet.run_canonical_route_b_wet_extension",
                return_value=base,
            ):
                result = run_canonical_route_b_family_semantic_wet_extension(
                    str(root), phenotype="fasciated meristem", case_id="case1"
                )

            state = result["family_semantic_resolution"]
            self.assertEqual("not_evaluated_noncanonical_family_copy_provenance", state["state"])
            self.assertEqual(wet, result["wet_candidates"])
            wrapper = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertTrue(wrapper["legacy_family_context_consumed"])
            self.assertEqual("none", wrapper["rank_contribution"])

    def test_string_false_provenance_is_noncanonical(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            wet = _wet(root / "wet.tsv")
            base_manifest = root / "base_manifest.json"
            base_manifest.write_text(json.dumps({
                "family_copy_producer_class": "canonical_supported_family_copy",
                "family_recovery_source": "supported_member_union",
                # JSON provenance must contain a boolean, not a string.
                "legacy_triple_route_family_context_consumed": "false",
            }), encoding="utf-8")
            base = {
                "state": "wet_reviewable",
                "manifest": str(base_manifest),
                "wet_candidates": wet,
                "wet_handoff_manifest": str(root / "handoff.json"),
            }
            with mock.patch(
                "candidate_qtl.canonical_family_semantic_wet.run_canonical_route_b_wet_extension",
                return_value=base,
            ):
                result = run_canonical_route_b_family_semantic_wet_extension(
                    str(root), phenotype="fasciated meristem", case_id="case1"
                )

            state = result["family_semantic_resolution"]
            self.assertEqual("not_evaluated_noncanonical_family_copy_provenance", state["state"])
            self.assertTrue(state["legacy_triple_route_family_context_consumed"])

    def test_non_string_producer_provenance_is_noncanonical(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            wet = _wet(root / "wet.tsv")
            base_manifest = root / "base_manifest.json"
            base_manifest.write_text(json.dumps({
                "family_copy_producer_class": ["canonical_supported_family_copy"],
                "family_recovery_source": "supported_member_union",
                "legacy_triple_route_family_context_consumed": False,
            }), encoding="utf-8")
            base = {
                "state": "wet_reviewable",
                "manifest": str(base_manifest),
                "wet_candidates": wet,
                "wet_handoff_manifest": str(root / "handoff.json"),
            }
            with mock.patch(
                "candidate_qtl.canonical_family_semantic_wet.run_canonical_route_b_wet_extension",
                return_value=base,
            ):
                result = run_canonical_route_b_family_semantic_wet_extension(
                    str(root), phenotype="fasciated meristem", case_id="case1"
                )

            state = result["family_semantic_resolution"]
            self.assertEqual("not_evaluated_noncanonical_family_copy_provenance", state["state"])
            self.assertEqual("", state["producer_class"])

    def test_no_wet_membership_never_blocks_base_review_state(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            base_manifest = root / "base_manifest.json"
            base_manifest.write_text(json.dumps({
                "state": "semantic_review_required",
            }), encoding="utf-8")
            base = {
                "state": "semantic_review_required",
                "manifest": str(base_manifest),
                "wet_candidates": None,
                "wet_handoff_manifest": None,
            }
            with mock.patch(
                "candidate_qtl.canonical_family_semantic_wet.run_canonical_route_b_wet_extension",
                return_value=base,
            ):
                result = run_canonical_route_b_family_semantic_wet_extension(
                    str(root), phenotype="fasciated meristem", case_id="case1"
                )

            self.assertEqual("semantic_review_required", result["state"])
            self.assertEqual(
                "not_run_without_frozen_wet_membership",
                result["family_semantic_resolution"]["state"],
            )
            self.assertIsNone(result["wet_candidates"])


if __name__ == "__main__":
    unittest.main()
