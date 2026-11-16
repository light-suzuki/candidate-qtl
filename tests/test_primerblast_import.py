from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main


def _caps_report(distinguishable: bool = True) -> dict:
    return {
        "target": {"name": "GeneA_snp1", "chrom": "chr1", "start": 100, "end": 600, "source": "snp"},
        "databases": ["pisum_v2"],
        "n_pairs": 2,
        "pairs": [
            {
                "name": "GeneA_snp1_P1",
                "forward": "AAAACCCCGGGGTTTTAAAA",
                "reverse": "TTTTGGGGCCCCAAAATTTT",
                "product_size": 430,
                "tm_f": 60.0, "tm_r": 60.2, "gc_f": 50.0, "gc_r": 55.0,
                "risk": "low", "risk_score": 0.0, "risk_reasons": [],
                "n_off_target": 0, "n_fr_offtarget": 0, "gel_distinguishable": True,
                "caps_enzyme": "MseI",
                "caps": {"best_enzyme": "MseI", "best_distinguishable": distinguishable, "min_gel_gap": 137, "n_candidate_enzymes": 2},
            },
            {
                "name": "GeneA_snp1_P2",
                "forward": "CCCCGGGGAAAATTTTCCCC",
                "reverse": "GGGGAAAATTTTCCCCGGGG",
                "product_size": 500,
                "tm_f": 59.0, "tm_r": 60.0, "gc_f": 55.0, "gc_r": 55.0,
                "risk": "high", "risk_score": 3.0, "risk_reasons": ["8 co-migrating off-target product(s)"],
                "n_off_target": 8, "n_fr_offtarget": 8, "gel_distinguishable": False,
                "caps_enzyme": "MseI",
                "caps": {"best_enzyme": "MseI", "best_distinguishable": distinguishable, "min_gel_gap": 5},
            },
        ],
    }


def _design_report() -> dict:
    return {
        "mode": "design",
        "template_id": "GeneB",
        "databases": ["pisum_v2", "cameor_v1"],
        "n_pairs": 1,
        "pairs": [
            {
                "index": 0,
                "forward": "ACGTACGTACGTACGTACGT",
                "reverse": "TGCATGCATGCATGCATGCA",
                "product_size": 368,
                "tm_f": 60.4, "tm_r": 60.1, "gc_f": 57.1, "gc_r": 45.5, "penalty": 3.5,
                "specificity": {"per_db": [
                    {"db": "pisum_v2", "n_off_target": 0, "n_comigrating": 0, "gel_distinguishable": True},
                    {"db": "cameor_v1", "n_off_target": 2, "n_comigrating": 1, "gel_distinguishable": True},
                ]},
                "gel_distinguishable": True,
            }
        ],
    }


class PrimerblastImportTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def _write(self, name: str, report: dict) -> Path:
        path = self.root / name
        path.write_text(json.dumps(report), encoding="utf-8")
        return path

    def test_caps_report_normalizes_and_ranks(self) -> None:
        report = self._write("caps.json", _caps_report())
        out = self.root / "out"
        self.assertEqual(0, main(["primerblast-import", "--results", str(report), "--out-dir", str(out)]))
        frame = pd.read_csv(out / "primer_candidates.tsv", sep="\t").fillna("")
        p1 = frame[frame["pair_name"] == "GeneA_snp1_P1"].iloc[0]
        p2 = frame[frame["pair_name"] == "GeneA_snp1_P2"].iloc[0]
        self.assertEqual(p1["marker_type"], "caps")
        self.assertEqual(p1["caps_enzyme"], "MseI")
        # The clean low-risk pair outranks the high-off-target one.
        self.assertEqual(int(p1["specificity_rank"]), 1)
        self.assertEqual(int(p2["specificity_rank"]), 2)
        # Predicted risk is never conflated with an experimental result.
        self.assertEqual(p1["validation_status"], "unvalidated")
        self.assertEqual(p1["experimental_status"], "pending_parent_f1_check")
        self.assertEqual(p1["parent_a_result"], "")
        self.assertIn("adopt_candidate", p1["selection_reason"])
        self.assertIn("review_or_reject", p2["selection_reason"])
        self.assertEqual(int(p1["offtarget_pisum_v2"]), 0)

    def test_dcaps_needed_when_no_distinguishing_enzyme(self) -> None:
        report = self._write("caps.json", _caps_report(distinguishable=False))
        out = self.root / "out"
        main(["primerblast-import", "--results", str(report), "--out-dir", str(out)])
        frame = pd.read_csv(out / "primer_candidates.tsv", sep="\t").fillna("")
        self.assertTrue((frame["marker_type"] == "dcaps_needed").any())
        row = frame[frame["marker_type"] == "dcaps_needed"].iloc[0]
        self.assertIn("dcaps", row["selection_reason"])

    def test_design_report_keeps_per_db_offtarget_columns(self) -> None:
        report = self._write("design.json", _design_report())
        out = self.root / "out"
        main(["primerblast-import", "--results", str(report), "--out-dir", str(out)])
        frame = pd.read_csv(out / "primer_candidates.tsv", sep="\t").fillna("")
        row = frame.iloc[0]
        self.assertEqual(row["mode"], "design")
        self.assertEqual(row["marker_type"], "sequencing_or_pcr")
        self.assertEqual(int(row["offtarget_pisum_v2"]), 0)
        self.assertEqual(int(row["offtarget_cameor_v1"]), 2)
        self.assertEqual(int(row["n_off_target_total"]), 2)

    def test_candidate_join_and_export(self) -> None:
        candidates = self.root / "candidates.tsv"
        candidates.write_text("gene_id\tmarker_id\nGeneA\tGeneA_snp1\n", encoding="utf-8")
        report = self._write("caps.json", _caps_report())
        out = self.root / "out"
        main(["primerblast-import", "--results", str(report), "--candidates", str(candidates), "--out-dir", str(out), "--export", "kasp"])
        frame = pd.read_csv(out / "primer_candidates.tsv", sep="\t").fillna("")
        self.assertTrue((frame["gene_id"] == "GeneA").all())
        self.assertTrue((out / "kasp_orders.tsv").exists())
        kasp = pd.read_csv(out / "kasp_orders.tsv", sep="\t").fillna("")
        self.assertIn("allele_specific_primer_a", kasp.columns)


if __name__ == "__main__":
    unittest.main()
