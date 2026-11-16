"""Well-formed is not the same as real (#189).

`classify_source` was hardened in #142/#153 so a free label cannot pose as a
citable accession, but grammar says nothing about existence. A model can emit
`UniProt:Q00000` -- correct shape, no such entry -- and two of them used to
anchor two independence groups, which is exactly what `top_eligibility` reads to
promote a gene to `independent_experimental_support`.

These tests pin the end-to-end consequence, not just the import bookkeeping: the
promotion has to stay closed after the rows reach `prioritize_evidence`.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.hardening import UNVERIFIED_LLM_SOURCE_GROUP, llm_evidence_import, prioritize_evidence


CANDIDATES = "gene_id\tchrom\tstart\tend\n" "GeneA\tchr1\t100\t200\n"


def _read_optional(path: Path) -> pd.DataFrame:
    """Read a TSV that is written headerless when there is nothing to report."""
    if not path.exists() or not path.read_text(encoding="utf-8").strip():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def _claim(code: str, source_id: str, **overrides: object) -> dict[str, object]:
    claim: dict[str, object] = {
        "evidence_code": code,
        "direction": "support",
        "confidence": "medium",
        "value": "claim",
        "source_id": source_id,
        "description": "curated claim",
        "reviewed": True,
        "reviewed_by": "Reviewer A",
        "reviewed_at": "2026-08-08",
    }
    claim.update(overrides)
    return claim


class LlmEvidenceImportGateTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.candidates = self.root / "candidates.tsv"
        self.candidates.write_text(CANDIDATES, encoding="utf-8")

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def _import(self, claims: list[dict[str, object]], resolve_against: list[str] | None = None) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
        payload = self.root / "claims.json"
        payload.write_text(json.dumps([{"gene_id": "GeneA", "claims": claims}]), encoding="utf-8")
        out = self.root / "llm_evidence.tsv"
        result = llm_evidence_import(str(self.candidates), str(payload), str(out), resolve_against)
        accepted = pd.read_csv(out, sep="\t", dtype=str).fillna("")
        rejected = _read_optional(Path(str(result["rejected"])))
        return result, accepted, rejected

    def test_reviewed_true_without_attribution_is_rejected(self) -> None:
        _, accepted, rejected = self._import([_claim("trait_function_match", "PMID:123", reviewed_by="")])
        self.assertTrue(accepted.empty)
        self.assertEqual(["missing_reviewer_attribution"], list(rejected["reason"]))

    def test_attribution_is_recorded_as_a_self_declaration(self) -> None:
        result, _, _ = self._import([_claim("trait_function_match", "PMID:123")])
        self.assertEqual("self_declared_not_verified", result["human_review_attestation"])
        audit = pd.read_csv(result["review_audit"], sep="\t", dtype=str).fillna("")
        self.assertEqual(["Reviewer A"], list(audit["reviewed_by"]))
        self.assertEqual(["2026-08-08"], list(audit["reviewed_at"]))

    def test_without_a_corpus_accepted_sources_are_marked_unverified(self) -> None:
        result, accepted, _ = self._import(
            [_claim("trait_function_match", "UniProt:Q00000"), _claim("domain_supported", "UniProt:Q00001")]
        )
        self.assertEqual("not_performed", result["source_verification"])
        self.assertEqual(2, result["unverified_source_count"])
        self.assertEqual([UNVERIFIED_LLM_SOURCE_GROUP] * 2, list(accepted["independence_group"]))

    def test_with_a_corpus_an_absent_accession_is_rejected(self) -> None:
        corpus = self.root / "proteome.fa"
        corpus.write_text(">sp|Q9LZS3|REAL_ARATH real protein\nMAAA\n", encoding="utf-8")
        result, accepted, rejected = self._import(
            [_claim("trait_function_match", "UniProt:Q00000"), _claim("domain_supported", "UniProt:Q9LZS3")],
            resolve_against=[str(corpus)],
        )
        self.assertEqual("performed", result["source_verification"])
        self.assertEqual(["unresolvable_source_id"], list(rejected["reason"]))
        self.assertEqual(["UniProt:Q9LZS3"], list(accepted["source_id"]))
        # A resolved source keeps anchoring its own group, as before.
        self.assertEqual([""], list(accepted["independence_group"]))
        self.assertEqual(1, result["resolved_source_count"])

    def test_corpus_matches_a_bare_accession_inside_a_fasta_header(self) -> None:
        corpus = self.root / "proteome.fa"
        corpus.write_text(">sp|Q9LZS3|REAL_ARATH\nMAAA\n", encoding="utf-8")
        _, accepted, _ = self._import(
            [_claim("trait_function_match", "Q9LZS3")], resolve_against=[str(corpus)]
        )
        self.assertEqual(["Q9LZS3"], list(accepted["source_id"]))

    def test_missing_corpus_file_is_a_clear_error(self) -> None:
        with self.assertRaisesRegex(ValueError, "--resolve-against file does not exist"):
            self._import([_claim("trait_function_match", "PMID:123")], resolve_against=[str(self.root / "absent.fa")])


class PromotionStaysClosedTests(unittest.TestCase):
    """The bookkeeping only matters if `top_eligibility` actually stays closed."""

    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.candidates = self.root / "candidates.tsv"
        self.candidates.write_text(CANDIDATES, encoding="utf-8")

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def _eligibility(self, evidence: pd.DataFrame) -> str:
        evidence_path = self.root / "evidence.tsv"
        evidence.to_csv(evidence_path, sep="\t", index=False)
        out_dir = self.root / "rank"
        prioritize_evidence(str(self.candidates), str(evidence_path), str(out_dir))
        ranked = pd.read_csv(out_dir / "ranked_candidates.tsv", sep="\t", dtype=str).fillna("")
        return str(ranked.loc[ranked["gene_id"].eq("GeneA"), "top_eligibility"].iloc[0])

    def _imported(self, claims: list[dict[str, object]], resolve_against: list[str] | None = None) -> pd.DataFrame:
        payload = self.root / "claims.json"
        payload.write_text(json.dumps([{"gene_id": "GeneA", "claims": claims}]), encoding="utf-8")
        out = self.root / "llm_evidence.tsv"
        llm_evidence_import(str(self.candidates), str(payload), str(out), resolve_against)
        return pd.read_csv(out, sep="\t", dtype=str).fillna("")

    def test_two_hallucinated_accessions_do_not_reach_independent_support(self) -> None:
        evidence = self._imported(
            [
                _claim("trait_function_match", "UniProt:Q00000"),
                _claim("domain_supported", "TAIR:AT1G99999"),
            ]
        )
        # Both are well-formed enough for classify_source, and they are in two
        # different independence families, so before #189 this promoted.
        self.assertNotEqual(
            "independent_experimental_support",
            self._eligibility(evidence),
        )

    def test_two_resolved_accessions_still_reach_independent_support(self) -> None:
        corpus = self.root / "corpus.tsv"
        corpus.write_text("Q00000\nAT1G99999\n", encoding="utf-8")
        evidence = self._imported(
            [
                _claim("trait_function_match", "UniProt:Q00000"),
                _claim("domain_supported", "TAIR:AT1G99999"),
            ],
            resolve_against=[str(corpus)],
        )
        self.assertEqual("independent_experimental_support", self._eligibility(evidence))


if __name__ == "__main__":
    unittest.main()
