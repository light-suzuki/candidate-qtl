from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "compare_qtl_sources.py"
SPEC = importlib.util.spec_from_file_location("compare_qtl_sources", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class QtlSourceComparisonTests(unittest.TestCase):
    def test_overlap_metrics_and_peak_distance(self) -> None:
        left = {
            "qtl_source_id": "Q1", "trait": "seed size", "target_assembly": "RefA",
            "target_interval": "Chr1:100-300", "target_peak": "Chr1:200",
        }
        right = {
            "qtl_source_id": "Q2", "trait": "seed size", "target_assembly": "RefA",
            "target_interval": "Chr1:200-400", "target_peak": "Chr1:250",
        }
        result = MODULE.compare(left, right)
        self.assertIsNotNone(result)
        assert result is not None
        self.assertEqual(101, result["overlap_bp"])
        self.assertEqual(50, result["peak_distance_bp"])
        self.assertTrue(result["trait_exact_match"])
        self.assertEqual("not_assessed", result["same_locus_claim"])

    def test_different_assemblies_are_not_compared(self) -> None:
        left = {"target_assembly": "RefA", "target_interval": "Chr1:1-10"}
        right = {"target_assembly": "RefB", "target_interval": "Chr1:1-10"}
        self.assertIsNone(MODULE.compare(left, right))

    def test_cli_writes_pairwise_table(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            source = root / "qtl.tsv"
            source.write_text(
                "qtl_source_id\ttrait\ttarget_assembly\ttarget_interval\ttarget_peak\n"
                "Q1\tseed size\tRefA\tChr1:100-300\tChr1:200\n"
                "Q2\tseed size\tRefA\tChr1:200-400\tChr1:250\n",
                encoding="utf-8",
            )
            output = root / "comparison.tsv"
            self.assertEqual(0, MODULE.main(["--input", str(source), "--out", str(output)]))
            self.assertIn("Q1", output.read_text(encoding="utf-8"))

    def test_cli_writes_header_when_no_pairs_are_comparable(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            source = root / "qtl.tsv"
            source.write_text(
                "qtl_source_id\ttrait\ttarget_assembly\ttarget_interval\ttarget_peak\n"
                "Q1\tseed size\tRefA\tChr1:100-300\tChr1:200\n",
                encoding="utf-8",
            )
            output = root / "comparison.tsv"
            self.assertEqual(0, MODULE.main(["--input", str(source), "--out", str(output)]))
            self.assertTrue(output.read_text(encoding="utf-8").startswith("qtl_a\tqtl_b\t"))


if __name__ == "__main__":
    unittest.main()
