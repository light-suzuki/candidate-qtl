"""Tests for bounded family/member/copy semantic resolution (#505)."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.family_semantic_resolution import (
    build_family_semantic_request,
    import_family_semantic_responses,
    verify_family_semantic_proposals,
)


def _write_copy(root: Path, rows: list[dict[str, str]]) -> Path:
    path = root / "family_copy.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, lineterminator="\n")
    return path


def _write_atlas(root: Path, rows: list[dict[str, str]]) -> Path:
    path = root / "atlas.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, lineterminator="\n")
    return path


def _responses(path: Path, prompt_id: str, resolutions: list[dict[str, object]]) -> Path:
    path.write_text(json.dumps({
        "candidates": [{
            "prompt_id": prompt_id,
            "parse_status": "json",
            "model_json": {"resolutions": resolutions},
        }]
    }, ensure_ascii=False), encoding="utf-8")
    return path


class FamilySemanticResolutionTests(unittest.TestCase):
    def test_cheap_first_skips_deterministically_resolved_family(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            copy = _write_copy(root, [
                {"family_id": "DONE", "candidate_gene_id": "D1", "disambiguation_status": "preferred_copy", "orthology_class": "one_to_one"},
                {"family_id": "AMB", "candidate_gene_id": "A1", "disambiguation_status": "paralog_ambiguous", "orthology_class": "one_to_one"},
                {"family_id": "AMB", "candidate_gene_id": "A2", "disambiguation_status": "paralog_ambiguous", "orthology_class": "one_to_many"},
            ])
            atlas = _write_atlas(root, [
                {"gene_id": "D1", "orthology_reference_gene": "AtDONE", "orthology_relationship": "1:1"},
                {"gene_id": "A1", "orthology_reference_gene": "AtAMB", "orthology_relationship": "1:1"},
                {"gene_id": "A2", "orthology_reference_gene": "AtAMB", "orthology_relationship": "one_to_many"},
            ])
            result = build_family_semantic_request(str(copy), str(atlas), str(root / "out"), case_id="c1")
            request = json.loads(Path(result["request"]).read_text(encoding="utf-8"))
            self.assertEqual(1, result["prompt_count"])
            self.assertEqual(["AMB"], [row["family_id"] for row in request["families"]])
            skipped = {row["family_id"]: row["reason"] for row in request["skipped_families"]}
            self.assertEqual("deterministic_copy_already_resolved", skipped["DONE"])
            self.assertTrue(request["cheap_first"])

    def test_proposal_alone_does_not_change_copy_state(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            copy = _write_copy(root, [
                {"family_id": "AMB", "candidate_gene_id": "A1", "disambiguation_status": "paralog_ambiguous", "orthology_class": "one_to_one"},
                {"family_id": "AMB", "candidate_gene_id": "A2", "disambiguation_status": "paralog_ambiguous", "orthology_class": "one_to_many"},
            ])
            atlas = _write_atlas(root, [
                {"gene_id": "A1", "orthology_reference_gene": "AtAMB", "orthology_relationship": "1:1"},
                {"gene_id": "A2", "orthology_reference_gene": "AtAMB", "orthology_relationship": "one_to_many"},
            ])
            request = build_family_semantic_request(str(copy), str(atlas), str(root / "request"), case_id="c1")
            prompt_id = json.loads(Path(request["request"]).read_text(encoding="utf-8"))["families"][0]["prompt_id"]
            response = _responses(root / "responses.json", prompt_id, [{
                "family_id": "AMB",
                "semantic_relation": "likely_reference_relative",
                "candidate_gene_ids_considered": ["A1", "A2"],
                "proposed_reference_mapping": [{"candidate_gene_id": "A1", "reference_gene": "AtAMB"}],
                "reasoning_summary": "A1 is the proposed reference-relative copy",
                "evidence_to_verify": ["orthology relationship"],
                "competing_explanations_to_keep": ["A2 paralog remains possible"],
                "review_required": True,
            }])
            imported = import_family_semantic_responses(
                request["request"], request["prompts"], str(response), str(root / "import"), case_id="c1"
            )
            original = pd.read_csv(copy, sep="\t", dtype=str).fillna("")
            self.assertEqual(["paralog_ambiguous", "paralog_ambiguous"], list(original["disambiguation_status"]))
            proposals = pd.read_csv(imported["proposals"], sep="\t", dtype=str).fillna("")
            self.assertEqual("proposed", proposals.iloc[0]["proposal_status"])
            self.assertEqual("none", proposals.iloc[0]["rank_contribution"])

    def test_strong_reference_orthology_promotes_only_verified_copy(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            copy = _write_copy(root, [
                {"family_id": "AMB", "candidate_gene_id": "A1", "disambiguation_status": "paralog_ambiguous", "orthology_class": "one_to_one"},
                {"family_id": "AMB", "candidate_gene_id": "A2", "disambiguation_status": "paralog_ambiguous", "orthology_class": "one_to_many"},
            ])
            atlas = _write_atlas(root, [
                {"gene_id": "A1", "orthology_reference_gene": "AtAMB", "orthology_relationship": "1:1"},
                {"gene_id": "A2", "orthology_reference_gene": "AtAMB", "orthology_relationship": "one_to_many"},
            ])
            request = build_family_semantic_request(str(copy), str(atlas), str(root / "request"), case_id="c1")
            prompt_id = json.loads(Path(request["request"]).read_text(encoding="utf-8"))["families"][0]["prompt_id"]
            response = _responses(root / "responses.json", prompt_id, [{
                "family_id": "AMB",
                "semantic_relation": "same_subfamily",
                "candidate_gene_ids_considered": ["A1", "A2"],
                "proposed_reference_mapping": [
                    {"candidate_gene_id": "A1", "reference_gene": "AtAMB"},
                    {"candidate_gene_id": "A2", "reference_gene": "AtAMB"},
                ],
                "review_required": True,
            }])
            imported = import_family_semantic_responses(
                request["request"], request["prompts"], str(response), str(root / "import"), case_id="c1"
            )
            verified = verify_family_semantic_proposals(
                imported["proposals"], str(copy), str(atlas), str(root / "verify"), case_id="c1"
            )
            frame = pd.read_csv(verified["family_copy"], sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual("preferred_copy", frame.loc["A1", "disambiguation_status"])
            self.assertEqual("paralog_ambiguous", frame.loc["A2", "disambiguation_status"])
            self.assertEqual("verified_copy_mapping", frame.loc["A1", "semantic_verification_status"])
            self.assertEqual("", frame.loc["A2", "semantic_verification_status"])
            shortlist = pd.read_csv(verified["shortlist"], sep="\t", dtype=str).fillna("")
            self.assertEqual(["A1"], list(shortlist["candidate_gene_id"]))
            self.assertEqual("first", shortlist.iloc[0]["first_wet_priority"])
            self.assertEqual(1, verified["promoted_family_count"])

    def test_two_machine_verified_mappings_preserve_tie(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            copy = _write_copy(root, [
                {"family_id": "AMB", "candidate_gene_id": "A1", "disambiguation_status": "family_supported_copy_unresolved", "orthology_class": "one_to_one"},
                {"family_id": "AMB", "candidate_gene_id": "A2", "disambiguation_status": "family_supported_copy_unresolved", "orthology_class": "one_to_one"},
            ])
            atlas = _write_atlas(root, [
                {"gene_id": "A1", "orthology_reference_gene": "AtAMB", "orthology_relationship": "reciprocal best hit"},
                {"gene_id": "A2", "orthology_reference_gene": "AtAMB", "orthology_relationship": "1:1"},
            ])
            request = build_family_semantic_request(str(copy), str(atlas), str(root / "request"), case_id="c1")
            prompt_id = json.loads(Path(request["request"]).read_text(encoding="utf-8"))["families"][0]["prompt_id"]
            response = _responses(root / "responses.json", prompt_id, [{
                "family_id": "AMB",
                "semantic_relation": "likely_reference_relative",
                "candidate_gene_ids_considered": ["A1", "A2"],
                "proposed_reference_mapping": [
                    {"candidate_gene_id": "A1", "reference_gene": "AtAMB"},
                    {"candidate_gene_id": "A2", "reference_gene": "AtAMB"},
                ],
                "review_required": True,
            }])
            imported = import_family_semantic_responses(
                request["request"], request["prompts"], str(response), str(root / "import"), case_id="c1"
            )
            verified = verify_family_semantic_proposals(
                imported["proposals"], str(copy), str(atlas), str(root / "verify"), case_id="c1"
            )
            frame = pd.read_csv(verified["family_copy"], sep="\t", dtype=str).fillna("")
            self.assertEqual({"co_preferred_tie"}, set(frame["disambiguation_status"]))
            shortlist = pd.read_csv(verified["shortlist"], sep="\t", dtype=str).fillna("")
            self.assertEqual({"A1", "A2"}, set(shortlist["candidate_gene_id"]))
            self.assertEqual({"co_first"}, set(shortlist["first_wet_priority"]))

    def test_scoring_fields_are_rejected_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            copy = _write_copy(root, [
                {"family_id": "AMB", "candidate_gene_id": "A1", "disambiguation_status": "family_resolution_required", "orthology_class": "one_to_one"},
                {"family_id": "AMB", "candidate_gene_id": "A2", "disambiguation_status": "family_resolution_required", "orthology_class": "one_to_many"},
            ])
            atlas = _write_atlas(root, [
                {"gene_id": "A1", "orthology_reference_gene": "AtAMB", "orthology_relationship": "1:1"},
                {"gene_id": "A2", "orthology_reference_gene": "AtAMB", "orthology_relationship": "one_to_many"},
            ])
            request = build_family_semantic_request(str(copy), str(atlas), str(root / "request"), case_id="c1")
            prompt_id = json.loads(Path(request["request"]).read_text(encoding="utf-8"))["families"][0]["prompt_id"]
            response = _responses(root / "responses.json", prompt_id, [{
                "family_id": "AMB",
                "semantic_relation": "same_subfamily",
                "probability": 0.91,
                "proposed_reference_mapping": [{"candidate_gene_id": "A1", "reference_gene": "AtAMB"}],
            }])
            with self.assertRaisesRegex(ValueError, "forbidden scoring fields"):
                import_family_semantic_responses(
                    request["request"], request["prompts"], str(response), str(root / "import"), case_id="c1"
                )


if __name__ == "__main__":
    unittest.main()
