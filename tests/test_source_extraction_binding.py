"""Synthetic proposal/source binding regressions (workflow2 issue #13)."""

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.claim_evidence import EXTRACTION_COLUMNS, gate_extraction
from candidate_qtl.source_verify import source_verify


CITATION = "PMID:12345678"
QUOTES = {
    "a": "GeneA loss of function caused altered growth.",
    "b": "GeneB expression was associated with altered growth.",
}
EXTRACTIONS = [
    {"hypothesis_id": key, "citation": CITATION, "quoted_support": QUOTES[key],
     "evidence_relation": relation}
    for key, relation in (("a", "direct_genetic_test"), ("b", "expression_association"))
]


class SourceExtractionBindingTests(unittest.TestCase):
    def test_scoped_same_paper_quotes_remain_distinct_in_both_orders(self):
        outputs = []
        for extractions in (EXTRACTIONS, list(reversed(EXTRACTIONS))):
            with tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                mechanisms = root / "mechanisms.tsv"
                pd.DataFrame([
                    {"hypothesis_id": key, "reference_gene_name": f"Gene{key.upper()}",
                     "source_ids": CITATION}
                    for key in QUOTES
                ]).to_csv(mechanisms, sep="\t", index=False)
                cache = root / "cache.json"
                cache.write_text(json.dumps([{
                    "id": "12345678", "pmid": "12345678", "title": "GeneA and GeneB",
                    "abstractText": " ".join(QUOTES.values()), "pubYear": "2020",
                }]), encoding="utf-8")
                extraction = root / "extraction.json"
                extraction.write_text(json.dumps(extractions), encoding="utf-8")
                out = root / "out"
                source_verify(str(mechanisms), str(out), str(cache), claim_extraction_arg=str(extraction))
                validation = pd.read_csv(out / "source_validation.tsv", sep="\t", dtype=str).fillna("")
                validation = validation.set_index("hypothesis_id")
                self.assertEqual(QUOTES, validation["quoted_support"].to_dict())
                self.assertEqual({"a": "supported_direct_genetic_test", "b": "associative_support_only"},
                                 validation["claim_support_status"].to_dict())
                rollup = pd.read_csv(out / "hypothesis_source_status.tsv", sep="\t", dtype=str)
                self.assertEqual({"a": "mechanistic_support", "b": "associative_support_only"},
                                 rollup.set_index("hypothesis_id")["claim_support_summary"].to_dict())
                records = json.loads((out / "verified_sources.json").read_text(encoding="utf-8"))
                for key, row in validation.iterrows():
                    self.assertEqual(row["claim_support_status"], records[CITATION][
                        "candidate_qtl_hypothesis_verifications"][key]["claim_support_status"])
                outputs.append(validation[["quoted_support", "claim_support_status"]].to_dict())
        self.assertEqual(*outputs)

    def test_invalid_bindings_fail_closed_before_annotation(self):
        first, second = EXTRACTIONS
        unscoped = {key: value for key, value in first.items() if key != "hypothesis_id"}
        cases = [
            ([{**first, "hypothesis_id": "missing"}], "unknown hypothesis_id"),
            ([{**first, "citation": "PMID:999"}], "target/citation mismatch"),
            ([first, {**first, "quoted_support": second["quoted_support"]}], "conflicting"),
            ([first, {**first, "evidence_relation": "expression_association"}], "conflicting"),
            ([unscoped, {**unscoped, "quoted_support": second["quoted_support"]}], "conflicting"),
            ([first, unscoped], "scoped/unscoped"),
        ]
        for rows, reason in cases:
            for ordered in (rows, list(reversed(rows))):
                with self.subTest(rows=ordered), self.assertRaisesRegex(ValueError, reason):
                    gate_extraction(ordered, {CITATION: " ".join(QUOTES.values())},
                                    {"a": {CITATION}, "b": {CITATION}})

    def test_identical_duplicates_and_legacy_columns(self):
        scoped = gate_extraction([EXTRACTIONS[0], EXTRACTIONS[0]],
                                 {CITATION: QUOTES["a"]}, {"a": {CITATION}})
        self.assertEqual(1, len(scoped))
        self.assertEqual("a", scoped[0]["hypothesis_id"])
        unscoped = {key: value for key, value in EXTRACTIONS[0].items() if key != "hypothesis_id"}
        legacy = gate_extraction([unscoped], {CITATION: QUOTES["a"]})
        self.assertEqual([*EXTRACTION_COLUMNS, "extraction_status", "claim_support_status"], list(legacy[0]))
        self.assertEqual("accepted", legacy[0]["extraction_status"])


if __name__ == "__main__":
    unittest.main()
