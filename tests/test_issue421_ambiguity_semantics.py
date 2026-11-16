from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from candidate_qtl.cli import main


class Issue421AmbiguitySemanticsTests(unittest.TestCase):
    """Regression contract for #421.

    A near-tied paralog is epistemic copy uncertainty.  It must remain visible
    in the evidence ledger, but it must not be converted into biological
    counterevidence or a negative numerical contribution to native dry rank.
    """

    def test_near_tied_paralog_is_context_and_never_scores_against_candidate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = root / "candidates.tsv"
            candidates.write_text(
                "gene_id\tannotation_quality\tcds_length\n"
                "GeneA\tcomplete_orf\t300\n"
                "GeneB\tcomplete_orf\t300\n",
                encoding="utf-8",
            )

            # GeneA has two near-tied similarity hits: this demonstrates that a
            # one-hit ortholog assignment is ambiguous, not that GeneA is less
            # likely to be causal for the QTL.
            blast = root / "blast.tsv"
            blast.write_text(
                "GeneA\tRelA1\t60\t80\t100\t120\t1\t80\t1\t80\t1e-30\t100\tbest hit\n"
                "GeneA\tRelA2\t58\t78\t100\t120\t1\t78\t1\t78\t1e-25\t95\tnear-tied paralog\n",
                encoding="utf-8",
            )

            # Keep one real scoring row in the fixture so score_contributions.tsv
            # is non-empty and the absence assertion below is meaningful.
            manual = root / "manual.tsv"
            manual.write_text(
                "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\n"
                "GeneB\tdirect_model_genetic\tsupport\thigh\tmutant\tliterature\tPMID:1\tdirect control evidence\n",
                encoding="utf-8",
            )

            evidence_dir = root / "evidence"
            self.assertEqual(
                0,
                main([
                    "evidence-build",
                    "--candidates", str(candidates),
                    "--blast", str(blast),
                    "--manual-evidence", str(manual),
                    "--out-dir", str(evidence_dir),
                ]),
            )
            with (evidence_dir / "evidence.tsv").open(encoding="utf-8") as handle:
                evidence_rows = list(csv.DictReader(handle, delimiter="\t"))
            ambiguity = [row for row in evidence_rows if row["evidence_code"] == "paralog_ambiguous"]
            self.assertEqual(1, len(ambiguity))
            self.assertEqual("GeneA", ambiguity[0]["gene_id"])
            self.assertEqual("context", ambiguity[0]["direction"])

            rank_dir = root / "rank"
            self.assertEqual(
                0,
                main([
                    "evidence-rank",
                    "--candidates", str(candidates),
                    "--evidence", str(evidence_dir / "evidence.tsv"),
                    "--out-dir", str(rank_dir),
                    "--top-n", "2",
                    "--reserve-n", "0",
                ]),
            )
            with (rank_dir / "score_contributions.tsv").open(encoding="utf-8") as handle:
                contributions = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(any(row["evidence_code"] == "direct_model_genetic" for row in contributions))
            self.assertFalse(any(row["evidence_code"] == "paralog_ambiguous" for row in contributions))


if __name__ == "__main__":
    unittest.main()
