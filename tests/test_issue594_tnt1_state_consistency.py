"""Issue #594: contradictory Tnt1 evidence states must fail closed."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.tnt1_evidence import normalize_tnt1_evidence


HEADER = (
    "reported_gene_id\tgenetic_background\tmutant_line\tfst_id\tinsertion_site\t"
    "insertion_observed\tpublished_mutant_phenotype\tsegregation_supported\t"
    "independent_alleles_supported\tcomplementation_or_rescue_supported\t"
    "phenotype_or_function\tgene_symbol\tontology_ids\tsource_ids\tsource_url\tnotes"
)


class Tnt1StateConsistencyTests(unittest.TestCase):
    def test_contradictory_rows_are_invalid_and_valid_strong_unknown_insertion_survives(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "tnt1.tsv"
            rows = [
                # Segregation cannot stand in for an unevaluated phenotype field.
                "medtr.R108.gnmHiC_1.ann1.C1\tR108\tNFC1\tFSTC1\tchr1:1\tsupported\t"
                "not_evaluated\tsupported\tnot_evaluated\tnot_evaluated\troot defect\t\t\t"
                "PMID:10001\thttps://example.org/c1\tsegregation without phenotype",
                # Multiple alleles cannot override an explicit phenotype contradiction.
                "medtr.R108.gnmHiC_1.ann1.C2\tR108\tNFC2\tFSTC2\tchr1:2\tsupported\t"
                "not_supported\tnot_evaluated\tsupported\tnot_evaluated\troot defect\t\t\t"
                "PMID:10002\thttps://example.org/c2\talleles contradict phenotype",
                # Rescue cannot override explicit no-insertion/no-phenotype states.
                "medtr.R108.gnmHiC_1.ann1.C3\tR108\tNFC3\tFSTC3\tchr1:3\tnot_supported\t"
                "not_supported\tnot_evaluated\tnot_evaluated\tsupported\troot defect\t\t\t"
                "PMID:10003\thttps://example.org/c3\trescue contradiction",
                # Segregation support also conflicts with explicit insertion absence.
                "medtr.R108.gnmHiC_1.ann1.C4\tR108\tNFC4\tFSTC4\tchr1:4\tnot_supported\t"
                "supported\tsupported\tnot_evaluated\tnot_evaluated\troot defect\t\t\t"
                "PMID:10004\thttps://example.org/c4\tsegregation insertion contradiction",
                # Unknown insertion remains unknown, not negative: published phenotype +
                # multiple independent alleles can still be a valid strong-causal record.
                "medtr.R108.gnmHiC_1.ann1.OK1\tR108\tNFOK1\tFSTOK1\tchr1:5\tnot_evaluated\t"
                "supported\tnot_evaluated\tsupported\tnot_evaluated\troot defect\t\t\t"
                "PMID:10005\thttps://example.org/ok\tvalid strong with unknown insertion",
            ]
            source.write_text(HEADER + "\n" + "\n".join(rows) + "\n", encoding="utf-8")

            result = normalize_tnt1_evidence(source, root / "out")
            audit = pd.read_csv(result["audit"], sep="\t", dtype=str).fillna("")
            ledger = pd.read_csv(result["ledger"], sep="\t", dtype=str).fillna("")
            functional = pd.read_csv(
                result["normalized_functional_evidence"], sep="\t", dtype=str
            ).fillna("")

            invalid = audit[audit.status.eq("invalid")]
            self.assertEqual(4, len(invalid))
            reasons = "\n".join(invalid.reason)
            self.assertIn("segregation_supported_requires_published_mutant_phenotype", reasons)
            self.assertIn("independent_alleles_supported_requires_published_mutant_phenotype", reasons)
            self.assertIn("complementation_or_rescue_supported_requires_published_mutant_phenotype", reasons)
            self.assertIn("segregation_supported_conflicts_with_insertion_not_supported", reasons)
            self.assertIn("strong_causal_flag_conflicts_with_insertion_not_supported", reasons)

            self.assertEqual(["medtr.R108.gnmHiC_1.ann1.OK1"], ledger.reported_gene_id.tolist())
            valid = ledger.iloc[0]
            self.assertEqual("not_evaluated", valid.insertion_observed)
            self.assertEqual("strong_causal", valid.evidence_class)
            self.assertEqual("true", valid.route_b_function_transfer_ready)
            self.assertEqual(["medtr.R108.gnmHiC_1.ann1.OK1"], functional.native_gene_id.tolist())
            self.assertEqual(4, result["counts"]["invalid_records"])
            self.assertEqual(1, result["counts"]["strong_causal"])
            self.assertEqual(1, result["counts"]["route_b_function_transfer_ready"])


if __name__ == "__main__":
    unittest.main()
