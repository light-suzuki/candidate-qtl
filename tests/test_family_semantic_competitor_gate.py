"""Final machine competitor-scan tests for family semantic resolution (#505)."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.family_semantic_competitor_gate import enforce_semantic_competitor_gate


def _write(path: Path, rows: list[dict[str, str]]) -> str:
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, lineterminator="\n")
    return str(path)


class FamilySemanticCompetitorGateTests(unittest.TestCase):
    def test_model_cannot_hide_an_equally_strong_competitor(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            original = _write(root / "original.tsv", [
                {"family_id": "AMB", "candidate_gene_id": "A1", "disambiguation_status": "paralog_ambiguous", "orthology_class": "one_to_one"},
                {"family_id": "AMB", "candidate_gene_id": "A2", "disambiguation_status": "paralog_ambiguous", "orthology_class": "one_to_one"},
            ])
            atlas = _write(root / "atlas.tsv", [
                {"gene_id": "A1", "orthology_reference_gene": "AtAMB", "orthology_relationship": "1:1"},
                {"gene_id": "A2", "orthology_reference_gene": "AtAMB", "orthology_relationship": "reciprocal best hit"},
            ])
            # Preliminary verifier saw only the model-proposed A1 mapping.
            preliminary = _write(root / "preliminary.tsv", [
                {"family_id": "AMB", "candidate_gene_id": "A1", "disambiguation_status": "preferred_copy", "orthology_class": "one_to_one", "semantic_match_basis": "llm_semantic_proposal+verified_reference_orthology", "semantic_verification_status": "verified_copy_mapping", "semantic_reference_gene": "AtAMB", "semantic_relation": "likely_reference_relative", "semantic_proposal_id": "p1"},
                {"family_id": "AMB", "candidate_gene_id": "A2", "disambiguation_status": "paralog_ambiguous", "orthology_class": "one_to_one", "semantic_match_basis": "", "semantic_verification_status": "", "semantic_reference_gene": "", "semantic_relation": "", "semantic_proposal_id": ""},
            ])

            result = enforce_semantic_competitor_gate(
                preliminary,
                original,
                atlas,
                str(root / "out"),
                case_id="case1",
            )
            final = pd.read_csv(result["family_copy"], sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual("co_preferred_tie", final.loc["A1", "disambiguation_status"])
            self.assertEqual("co_preferred_tie", final.loc["A2", "disambiguation_status"])
            self.assertEqual("verified_co_preferred_competitor", final.loc["A2", "semantic_verification_status"])
            self.assertEqual(1, result["competitor_expansion_count"])
            shortlist = pd.read_csv(result["shortlist"], sep="\t", dtype=str).fillna("")
            self.assertEqual({"A1", "A2"}, set(shortlist["candidate_gene_id"]))
            self.assertEqual({"co_first"}, set(shortlist["first_wet_priority"]))

    def test_unique_machine_relation_can_remain_preferred(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            original = _write(root / "original.tsv", [
                {"family_id": "AMB", "candidate_gene_id": "A1", "disambiguation_status": "family_resolution_required", "orthology_class": "one_to_one"},
                {"family_id": "AMB", "candidate_gene_id": "A2", "disambiguation_status": "family_resolution_required", "orthology_class": "one_to_many"},
            ])
            atlas = _write(root / "atlas.tsv", [
                {"gene_id": "A1", "orthology_reference_gene": "AtAMB", "orthology_relationship": "1:1"},
                {"gene_id": "A2", "orthology_reference_gene": "AtAMB", "orthology_relationship": "one_to_many"},
            ])
            preliminary = _write(root / "preliminary.tsv", [
                {"family_id": "AMB", "candidate_gene_id": "A1", "disambiguation_status": "preferred_copy", "orthology_class": "one_to_one", "semantic_match_basis": "llm_semantic_proposal+verified_reference_orthology", "semantic_verification_status": "verified_copy_mapping", "semantic_reference_gene": "AtAMB", "semantic_relation": "same_subfamily", "semantic_proposal_id": "p1"},
                {"family_id": "AMB", "candidate_gene_id": "A2", "disambiguation_status": "family_resolution_required", "orthology_class": "one_to_many", "semantic_match_basis": "", "semantic_verification_status": "", "semantic_reference_gene": "", "semantic_relation": "", "semantic_proposal_id": ""},
            ])

            result = enforce_semantic_competitor_gate(
                preliminary,
                original,
                atlas,
                str(root / "out"),
                case_id="case1",
            )
            final = pd.read_csv(result["family_copy"], sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual("preferred_copy", final.loc["A1", "disambiguation_status"])
            self.assertEqual("family_resolution_required", final.loc["A2", "disambiguation_status"])
            self.assertEqual(0, result["competitor_expansion_count"])


if __name__ == "__main__":
    unittest.main()
