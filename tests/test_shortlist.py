from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main


class ScreeningShortlistTests(unittest.TestCase):
    def test_caps_primary_twenty_and_exports_boundary_tie(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ranked = root / "ranked.tsv"
            rows = []
            for index in range(1, 26):
                rows.append({
                    "gene_id": f"Gene{index:02d}",
                    "dry_score": "1.0" if index <= 22 else "0.0",
                    "screening_rank": "1" if index <= 22 else "23",
                    "support_family_count": "1",
                    "independent_source_count": "1",
                    "dry_group": "evidence_insufficient",
                    "top_eligibility": "insufficient_independent_evidence",
                    "rescue_reason": "annotation_incomplete" if index in {24, 25} else "",
                    "rescue_priority": "3" if index in {24, 25} else "0",
                })
            pd.DataFrame(rows).to_csv(ranked, sep="\t", index=False)
            out = root / "shortlist"
            self.assertEqual(0, main(["screening-shortlist", "--ranked", str(ranked), "--out-dir", str(out)]))

            primary = pd.read_csv(out / "deterministic_shortlist.tsv", sep="\t", dtype=str).fillna("")
            boundary = pd.read_csv(out / "boundary_ties.tsv", sep="\t", dtype=str).fillna("")
            rescue = pd.read_csv(out / "rescue_appendix.tsv", sep="\t", dtype=str).fillna("")
            manifest = json.loads((out / "shortlist_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(20, len(primary))
            self.assertEqual(22, len(boundary))
            self.assertEqual({"true", "false"}, set(boundary["selected_in_primary"].str.lower()))
            self.assertEqual({"Gene24", "Gene25"}, set(rescue["gene_id"]))
            self.assertTrue(manifest["boundary_tie_truncated"])
            self.assertFalse(manifest["llm_used"])

    def test_stage_a_does_not_filter_on_top_eligibility(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ranked = root / "ranked.tsv"
            pd.DataFrame([
                {"gene_id": "InsufficientButFirst", "dry_score": 2, "screening_rank": 1,
                 "dry_group": "evidence_insufficient", "top_eligibility": "insufficient_independent_evidence"},
                {"gene_id": "EligibleButSecond", "dry_score": 1, "screening_rank": 2,
                 "dry_group": "top", "top_eligibility": "strong_genetic_evidence"},
            ]).to_csv(ranked, sep="\t", index=False)
            out = root / "out"
            self.assertEqual(0, main(["screening-shortlist", "--ranked", str(ranked), "--out-dir", str(out), "--target-size", "1"]))
            primary = pd.read_csv(out / "deterministic_shortlist.tsv", sep="\t", dtype=str)
            self.assertEqual("InsufficientButFirst", primary.iloc[0]["gene_id"])


if __name__ == "__main__":
    unittest.main()
