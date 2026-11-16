from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.orthology_staged import load_runtime_profile


def _blast_row(q, s, bitscore):
    return "\t".join([q, s, "100", "100", "100", "100", "1", "100", "1", "100", "1e-50", str(bitscore)])


class RuntimeProfileTests(unittest.TestCase):
    def test_builtin_profiles(self) -> None:
        laptop = load_runtime_profile("laptop16")
        self.assertEqual(int(laptop["threads"]), 4)
        self.assertEqual(int(laptop["max_hits_per_query"]), 25)
        accuracy = load_runtime_profile("accuracy64")
        self.assertEqual(int(accuracy["threads"]), 8)
        self.assertEqual(int(accuracy["max_hits_per_query"]), 50)

    def test_unknown_profile_raises(self) -> None:
        with self.assertRaises(ValueError):
            load_runtime_profile("phone2")


class Tier1PlanTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        (self.root / "candidates.tsv").write_text("gene_id\nG1\nG2\n", encoding="utf-8")
        (self.root / "proteins.fa").write_text(">G1\nMKMK\n>G2\nMLML\n>OUT\nMPMP\n", encoding="utf-8")
        (self.root / "focal.fa").write_text(">G1\nMKMK\n>G2\nMLML\n>PX\nMKMR\n", encoding="utf-8")

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_profile_controls_threads_not_labels(self) -> None:
        out = self.root / "tier1"
        self.assertEqual(
            0,
            main([
                "orthology-tier1-plan", "--candidates", str(self.root / "candidates.tsv"),
                "--protein-fasta", str(self.root / "proteins.fa"), "--focal-proteome", str(self.root / "focal.fa"),
                "--out-dir", str(out), "--profile", "laptop16", "--aligner", "blastp",
            ]),
        )
        manifest = json.loads((out / "orthology_tier1_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["threads"], 4)
        self.assertEqual(manifest["max_hits_per_query"], 25)
        # Query is only the two interval candidates, never the whole proteome.
        self.assertEqual(manifest["query_candidate_count"], 2)
        self.assertTrue(manifest["query_is_candidates_only"])
        self.assertFalse(manifest["runs_full_proteome_all_vs_all"])
        script = (out / "run_tier1_orthology.sh").read_text(encoding="utf-8")
        self.assertIn("focal_self", script)

        out64 = self.root / "tier1_64"
        main([
            "orthology-tier1-plan", "--candidates", str(self.root / "candidates.tsv"),
            "--protein-fasta", str(self.root / "proteins.fa"), "--focal-proteome", str(self.root / "focal.fa"),
            "--out-dir", str(out64), "--profile", "accuracy64", "--aligner", "blastp",
        ])
        manifest64 = json.loads((out64 / "orthology_tier1_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest64["threads"], 8)
        self.assertEqual(manifest64["max_hits_per_query"], 50)


class FocalParalogImportTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        (self.root / "candidates.tsv").write_text("gene_id\nG1\nG2\nG3\nG4\n", encoding="utf-8")
        rows = [
            _blast_row("G1", "G1", 200), _blast_row("G1", "PX", 160),
            _blast_row("G2", "G2", 200), _blast_row("G2", "G3", 185),
            _blast_row("G3", "G3", 200), _blast_row("G3", "G2", 185),
            _blast_row("G4", "G4", 200),
        ]
        (self.root / "focal_blast.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_out_of_interval_and_in_interval_paralogs(self) -> None:
        out = self.root / "focal_evidence.tsv"
        self.assertEqual(
            0,
            main(["focal-paralog-import", "--focal-blast", str(self.root / "focal_blast.tsv"), "--candidates", str(self.root / "candidates.tsv"), "--out", str(out)]),
        )
        summary = pd.read_csv(out.with_name("focal_evidence_paralog_summary.tsv"), sep="\t").fillna("")
        by_gene = {row["gene_id"]: row for _, row in summary.iterrows()}
        self.assertEqual(by_gene["G1"]["relationship_class"], "out_of_interval_paralog")
        self.assertEqual(int(by_gene["G1"]["out_of_interval_paralog_count"]), 1)
        self.assertEqual(by_gene["G2"]["relationship_class"], "in_interval_paralog")
        self.assertEqual(by_gene["G4"]["relationship_class"], "single_copy_in_focal")

        evidence = pd.read_csv(out, sep="\t").fillna("")
        g1 = evidence[evidence["gene_id"] == "G1"].iloc[0]
        self.assertEqual(g1["evidence_code"], "paralog_ambiguous")
        self.assertIn("out_of_interval_paralog", g1["value"])
        # A single-copy gene produces no redundancy evidence.
        self.assertTrue(evidence[evidence["gene_id"] == "G4"].empty)


if __name__ == "__main__":
    unittest.main()
