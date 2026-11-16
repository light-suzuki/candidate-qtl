from __future__ import annotations

import unittest

from candidate_qtl import auto_triage
from candidate_qtl import homolog_panel


class ScreeningFunctionWordTests(unittest.TestCase):
    def _input_tokens(self, phrase: str) -> set[str]:
        return {
            str(row["query"]).casefold()
            for row in homolog_panel._mechanical_queries(phrase)
            if row["origin"] == "input_token"
        }

    def test_shared_helper_filters_demonstrated_structural_words(self) -> None:
        self.assertIs(auto_triage._mechanical_queries, homolog_panel._mechanical_queries)

        flowering = self._input_tokens("flowering under short days")
        self.assertNotIn("under", flowering)
        self.assertTrue({"flowering", "short", "days"}.issubset(flowering))

        measured = self._input_tokens("height measured from its base to its tip")
        self.assertNotIn("from", measured)
        self.assertNotIn("its", measured)
        self.assertIn("height", measured)

        expressed = self._input_tokens("genes that are expressed")
        self.assertNotIn("that", expressed)
        self.assertNotIn("are", expressed)
        self.assertIn("genes", expressed)

    def test_measurement_shell_words_are_not_independent_queries(self) -> None:
        tokens = self._input_tokens(
            "protein level number amount content component components type types"
        )
        self.assertFalse(
            {"level", "number", "amount", "content", "component", "components", "type", "types"}
            & tokens
        )
        self.assertIn("protein", tokens)

    def test_structural_one_token_overlap_cannot_create_support(self) -> None:
        for token in ("under", "from", "that", "are", "level", "number", "component"):
            with self.subTest(token=token):
                self.assertFalse(
                    homolog_panel._overlap(
                        token,
                        f"A source-backed annotation where {token} appears in ordinary prose.",
                    )
                )

    def test_legitimate_one_token_biological_query_still_matches(self) -> None:
        self.assertTrue(
            homolog_panel._overlap(
                "meristem",
                "CLAVATA signaling maintains the shoot apical meristem.",
            )
        )
        self.assertTrue(
            homolog_panel._overlap(
                "flowering",
                "Mutation changes flowering time under short days.",
            )
        )

    def test_structural_word_does_not_upgrade_genetic_comment(self) -> None:
        record = {
            "primaryAccession": "P00001",
            "annotationScore": 5,
            "organism": {"scientificName": "Arabidopsis thaliana", "lineage": ["Viridiplantae"]},
            "genes": [{"geneName": {"value": "EXAMPLE"}}],
            "proteinDescription": {"recommendedName": {"fullName": {"value": "example protein"}}},
            "comments": [{
                "commentType": "DISRUPTION PHENOTYPE",
                "texts": [{
                    "value": "Plants flower later under restrictive conditions.",
                    "evidences": [{"source": "PubMed", "id": "12345"}],
                }],
            }],
        }
        scored = homolog_panel._score_record(record, "under", "input_token")
        self.assertEqual("database_context", scored["evidence_level"])
        self.assertEqual([], scored["channels"])
        self.assertEqual([], scored["pmids"])


if __name__ == "__main__":
    unittest.main()
