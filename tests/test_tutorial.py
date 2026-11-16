"""End-to-end tutorial and known-QTG ranking benchmark.

This test runs the offline candidate-first chain on the redistributable
docs/tutorial fixtures so CI exercises the whole workflow, and scores a small
known-QTG benchmark that records misses (a causal gene not recovered) as
failures rather than hiding them.
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main


REPO = Path(__file__).resolve().parents[1]
TUTORIAL = REPO / "docs" / "tutorial"


class TutorialEndToEndTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def _run_chain(self) -> pd.DataFrame:
        extract = self.root / "extract"
        self.assertEqual(0, main([
            "extract", "--fasta", str(TUTORIAL / "genome.fa"), "--gff3", str(TUTORIAL / "genes.gff3"),
            "--interval", "chr1:1-900", "--out-dir", str(extract),
        ]))
        candidates = extract / "candidates.tsv"
        self.assertEqual(set(pd.read_csv(candidates, sep="\t")["gene_id"]), {"GENE_CAUSAL", "DECOY1", "DECOY2"})

        evidence = self.root / "evidence"
        self.assertEqual(0, main([
            "evidence-build", "--candidates", str(candidates),
            "--eggnog", str(TUTORIAL / "candidate_eggnog.tsv"),
            "--expression", str(TUTORIAL / "expression.tsv"), "--target-tissue", "meristem",
            "--manual-evidence", str(TUTORIAL / "manual_evidence.tsv"),
            "--out-dir", str(evidence),
        ]))

        rank = self.root / "rank"
        self.assertEqual(0, main([
            "evidence-rank", "--candidates", str(candidates), "--evidence", str(evidence / "evidence.tsv"),
            "--out-dir", str(rank),
        ]))
        # Ranking metrics and weights are pinned in a reproducible file.
        self.assertTrue((rank / "weights_used.tsv").exists())
        return pd.read_csv(rank / "ranked_candidates.tsv", sep="\t").fillna("")

    def test_offline_chain_recovers_the_causal_gene(self) -> None:
        ranked = self._run_chain()
        by_gene = {row["gene_id"]: row for _, row in ranked.iterrows()}
        self.assertEqual(by_gene["GENE_CAUSAL"]["dry_group"], "top")
        self.assertEqual(int(by_gene["GENE_CAUSAL"]["dry_rank"]), 1)
        self.assertEqual(by_gene["GENE_CAUSAL"]["top_eligibility"], "strong_genetic_evidence")
        # Evidence-poor decoys never reach Top.
        self.assertNotEqual(by_gene["DECOY2"]["dry_group"], "top")

    def test_offline_side_commands_run(self) -> None:
        # Trait-term control and ontology mapping are local-only.
        terms = self.root / "terms.tsv"
        terms.write_text("term\tterm_type\tapproved\nfasciation\tdirect_phenotype\ttrue\n", encoding="utf-8")
        approved = self.root / "approved_terms.tsv"
        self.assertEqual(0, main([
            "trait-terms", "build", "--trait-terms", str(terms),
            "--ontology-map", str(REPO / "templates" / "trait_ontology_map.tsv"), "--out", str(approved),
        ]))
        extract = self.root / "extract_side"
        main(["extract", "--fasta", str(TUTORIAL / "genome.fa"), "--gff3", str(TUTORIAL / "genes.gff3"), "--interval", "chr1:1-900", "--out-dir", str(extract)])
        self.assertEqual(0, main([
            "peak-context", "--candidates", str(extract / "candidates.tsv"),
            "--interval", "chr1:1-900", "--peak", "chr1:300", "--out-dir", str(self.root / "peak"),
        ]))

    def test_two_assembly_qtl_source_comparison(self) -> None:
        out = self.root / "compare.tsv"
        completed = subprocess.run(
            [sys.executable, str(REPO / "scripts" / "compare_qtl_sources.py"), "--input", str(TUTORIAL / "qtl_sources.tsv"), "--out", str(out)],
            capture_output=True, text=True,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        compared = pd.read_csv(out, sep="\t").fillna("")
        # The same mock QTL is lifted from two assemblies onto the reference.
        self.assertGreaterEqual(len(compared), 1)
        self.assertEqual(compared.iloc[0]["qtl_a"], "TUT_asm1")

    def test_known_qtg_benchmark_records_misses(self) -> None:
        ranked = self._run_chain()
        by_gene = {row["gene_id"]: row for _, row in ranked.iterrows()}
        benchmark = pd.read_csv(TUTORIAL / "known_qtg_benchmark.tsv", sep="\t").fillna("")

        results = []
        for _, case in benchmark.iterrows():
            causal = case["causal_gene"]
            row = by_gene.get(causal)
            recovered_group = row is not None and row["dry_group"] in {"top", "reserve", "rescue"}
            top1 = row is not None and row["dry_group"] == "top" and int(row["dry_rank"]) == 1
            results.append({
                "case_id": case["case_id"],
                "causal_gene": causal,
                "variant_class": case["causal_variant_class"],
                "recovered": recovered_group,
                "top1": top1,
                "outcome": "recovered" if recovered_group else "miss",
            })

        results_frame = pd.DataFrame(results)
        coding = results_frame[results_frame["case_id"] == "tut_coding"].iloc[0]
        regulatory = results_frame[results_frame["case_id"] == "tut_regulatory"].iloc[0]

        # The coding causal gene is recovered top-1; the evidence-poor regulatory
        # case is recorded as a miss instead of being silently dropped.
        self.assertTrue(bool(coding["top1"]))
        self.assertEqual(coding["outcome"], "recovered")
        self.assertFalse(bool(regulatory["recovered"]))
        self.assertEqual(regulatory["outcome"], "miss")

        # Recall metrics are computable and a miss is present in the record.
        top1_recall = results_frame["top1"].mean()
        self.assertAlmostEqual(top1_recall, 0.5, places=3)
        self.assertIn("miss", set(results_frame["outcome"]))


if __name__ == "__main__":
    unittest.main()
