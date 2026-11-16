from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import build_parser, main


PLAN_COLUMNS = "gene_id\ttarget_interval\ttarget_length\tlocal_target_start\tlocal_target_end\toutput_json\n"


class ExonCoverageCheckTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_detects_gaps_and_plans_retries(self) -> None:
        covered_json = self.root / "w1.json"
        covered_json.write_text(json.dumps({"mode": "design", "pairs": [{"left_start": 0, "right_start": 100}]}), encoding="utf-8")
        gap_json = self.root / "w2.json"
        gap_json.write_text(json.dumps({"mode": "design", "pairs": [{"left_start": 20, "right_start": 30}]}), encoding="utf-8")

        plan = self.root / "exon_sequencing_plan.tsv"
        plan.write_text(
            PLAN_COLUMNS
            + f"GeneA\tchr1:100-140\t41\t10\t50\t{covered_json}\n"
            + f"GeneA\tchr1:200-240\t41\t10\t50\t{gap_json}\n"
            + f"GeneB\tchr1:300-340\t41\t10\t50\t{self.root / 'missing.json'}\n",
            encoding="utf-8",
        )
        out = self.root / "coverage"
        self.assertEqual(0, main(["exon-coverage-check", "--plan", str(plan), "--out-dir", str(out)]))

        manifest = json.loads((out / "coverage_check_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["covered_windows"], 1)
        self.assertEqual(manifest["coverage_gap_windows"], 1)
        self.assertEqual(manifest["result_missing_windows"], 1)
        self.assertFalse(manifest["complete"])

        by_gene = pd.read_csv(out / "coverage_by_gene.tsv", sep="\t").fillna("")
        gene_a = by_gene[by_gene["gene_id"] == "GeneA"].iloc[0]
        self.assertAlmostEqual(float(gene_a["covered_fraction"]), 0.5, places=3)
        self.assertFalse(bool(gene_a["fully_covered"]))

        retry = pd.read_csv(out / "coverage_retry_plan.tsv", sep="\t").fillna("")
        # Both the gap and the missing-result windows get a retry row; covered does not.
        self.assertEqual(len(retry), 2)
        self.assertIn("chr1:200-240", set(retry["target_interval"]))
        self.assertNotIn("chr1:100-140", set(retry["target_interval"]))


class StableCliReachabilityTests(unittest.TestCase):
    def test_extension_commands_are_in_the_stable_cli(self) -> None:
        parser = build_parser()
        # subparsers action holds the registered command names.
        choices = set()
        for action in parser._actions:
            if hasattr(action, "choices") and action.choices:
                choices.update(action.choices)
        for command in ["research-packet", "llm-evidence-import", "exon-sequencing-plan", "exon-coverage-check"]:
            self.assertIn(command, choices)

    def test_evidence_rank_exposes_rescue_n(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            (root / "cand.tsv").write_text("gene_id\nG1\n", encoding="utf-8")
            (root / "ev.tsv").write_text(
                "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\n"
                "G1\tknown_qtg\tsupport\thigh\tx\tmanual\tPMID:1\tk\n",
                encoding="utf-8",
            )
            out = root / "rank"
            self.assertEqual(
                0,
                main(["evidence-rank", "--candidates", str(root / "cand.tsv"), "--evidence", str(root / "ev.tsv"), "--out-dir", str(out), "--rescue-n", "1"]),
            )


if __name__ == "__main__":
    unittest.main()
