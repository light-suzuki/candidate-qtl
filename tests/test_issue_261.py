"""Cross-code counterevidence conflicts must demote from Top (issue #261)."""

from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from candidate_qtl.cli import main
from candidate_qtl.hardening import prioritize_evidence


class Issue261CrossCodeConflictTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def _rank(self, evidence_text: str) -> dict[str, dict[str, str]]:
        candidates = self.root / "candidates.tsv"
        candidates.write_text("gene_id\tannotation_quality\nGeneX\tcomplete_orf\n", encoding="utf-8")
        evidence = self.root / "evidence.tsv"
        evidence.write_text(
            "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\n"
            + evidence_text,
            encoding="utf-8",
        )
        out = self.root / "rank"
        prioritize_evidence(str(candidates), str(evidence), str(out), top_n=1, reserve_n=0, rescue_n=0)
        with (out / "ranked_candidates.tsv").open(encoding="utf-8") as handle:
            return {row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}

    def test_cross_code_strong_support_with_counterevidence_is_conflicting(self) -> None:
        ranked = self._rank(
            "GeneX\tdirect_model_genetic\tsupport\thigh\tmutant\tmanual\tPMID:1\tdirect\n"
            "GeneX\tcounterevidence\tagainst\thigh\tcontradiction\tgenetics\tPMID:2\tcontradiction\n"
        )
        row = ranked["GeneX"]
        self.assertEqual("2.0", row["dry_score"])
        self.assertEqual("direct_model_genetic+counterevidence", row["evidence_conflicts"])
        self.assertEqual("conflicting_evidence", row["top_eligibility"])
        self.assertNotEqual("top", row["dry_group"])

    def test_same_code_support_and_against_is_still_conflicting(self) -> None:
        ranked = self._rank(
            "GeneX\tdirect_model_genetic\tsupport\thigh\tmutant\tmanual\tPMID:1\tdirect\n"
            "GeneX\tdirect_model_genetic\tagainst\thigh\tnormal\tmanual\tPMID:2\tcounterexample\n"
        )
        row = ranked["GeneX"]
        self.assertEqual("direct_model_genetic", row["evidence_conflicts"])
        self.assertEqual("conflicting_evidence", row["top_eligibility"])
        self.assertNotEqual("top", row["dry_group"])

    def test_counterevidence_against_without_support_is_not_conflicting(self) -> None:
        ranked = self._rank("GeneX\tcounterevidence\tagainst\thigh\tcontradiction\tgenetics\tPMID:1\tcontradiction\n")
        row = ranked["GeneX"]
        self.assertEqual("", row["evidence_conflicts"])
        self.assertNotEqual("conflicting_evidence", row["top_eligibility"])

    def test_non_counterevidence_against_does_not_conflict(self) -> None:
        ranked = self._rank(
            "GeneX\tdirect_model_genetic\tsupport\thigh\tmutant\tmanual\tPMID:1\tdirect\n"
            "GeneX\texpression_absent\tagainst\thigh\tabsent\trnaseq\tPMID:2\tno expression\n"
        )
        row = ranked["GeneX"]
        self.assertEqual("4.5", row["dry_score"])
        self.assertEqual("", row["evidence_conflicts"])
        self.assertEqual("strong_genetic_evidence", row["top_eligibility"])
        self.assertEqual("top", row["dry_group"])

    def test_case_run_demotes_cross_code_conflict_from_top(self) -> None:
        candidates = self.root / "candidates.tsv"
        candidates.write_text(
            "gene_id\tannotation_quality\nGeneClean\tcomplete_orf\nGeneConflict\tcomplete_orf\n", encoding="utf-8"
        )
        evidence = self.root / "evidence.tsv"
        evidence.write_text(
            "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\n"
            "GeneClean\tdirect_model_genetic\tsupport\thigh\tmutant\tgenetics\tPMID:1\tdirect\n"
            "GeneConflict\tdirect_model_genetic\tsupport\thigh\tmutant\tgenetics\tPMID:2\tdirect\n"
            "GeneConflict\tcounterevidence\tagainst\thigh\tcontradiction\tgenetics\tPMID:3\tcontradiction\n",
            encoding="utf-8",
        )
        case = self.root / "case.yaml"
        resources = self.root / "resources.yaml"
        out = self.root / "out"
        case.write_text(
            "case_id: cross_conflict\nsource:\n  source_id: PMID:QTL\n  target_interval: chr1:1-10\n"
            "  literature_search_status: complete\ninputs:\n  candidates_resource: candidates\n"
            "  evidence_resources: [evidence]\nassay:\n  top_n: 1\n",
            encoding="utf-8",
        )
        resources.write_text(
            f"resources:\n  candidates:\n    path: {candidates}\n  evidence:\n    path: {evidence}\n",
            encoding="utf-8",
        )
        self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources), "--out-dir", str(out)]))
        with (out / "rank" / "ranked_candidates.tsv").open(encoding="utf-8") as handle:
            ranked = {row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual("top", ranked["GeneClean"]["dry_group"])
        self.assertEqual("conflicting_evidence", ranked["GeneConflict"]["top_eligibility"])
        self.assertNotEqual("top", ranked["GeneConflict"]["dry_group"])


if __name__ == "__main__":
    unittest.main()
