from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.dual_route import phenotype_first_status


class RouteACoarseIntervalMemberTests(unittest.TestCase):
    def _run(self, root: Path, rows: list[dict[str, str]]) -> pd.DataFrame:
        hypotheses = root / "hypotheses.tsv"
        falsification = root / "falsification.tsv"
        pd.DataFrame([
            {
                "hypothesis_id": "mech_004",
                "pathway_or_process": "floral transition",
                "gene_family": "PEBP/CETS family, TFL1-clade",
            }
        ]).to_csv(hypotheses, sep="\t", index=False)
        pd.DataFrame(rows).to_csv(falsification, sep="\t", index=False)

        result = phenotype_first_status(
            str(root / "out"),
            hypotheses_arg=str(hypotheses),
            falsification_arg=str(falsification),
            case_id="pea_det",
        )
        return pd.read_csv(result["hypotheses"], sep="\t", dtype=str).fillna("")

    def test_fine_empty_does_not_erase_coarse_interval_member(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            frame = self._run(Path(tmp), [
                {
                    "hypothesis_id": "mech_004",
                    "level": "coarse",
                    "verdict": "surviving_enriched",
                    "interval_gene_ids": "Psat.cameor.v2.3g257650",
                },
                {
                    "hypothesis_id": "mech_004",
                    "level": "fine",
                    "verdict": "absent_annotation_only",
                    "interval_gene_ids": "",
                },
            ])
            row = frame.iloc[0]
            self.assertEqual("Psat.cameor.v2.3g257650", row["member_interval_genes"])
            self.assertEqual("1", row["member_count"])
            self.assertEqual("interval_anchored", row["anchor_status"])

    def test_fine_resolved_member_overrides_coarse_member(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            frame = self._run(Path(tmp), [
                {
                    "hypothesis_id": "mech_004",
                    "level": "coarse",
                    "verdict": "surviving_enriched",
                    "interval_gene_ids": "coarse_gene",
                },
                {
                    "hypothesis_id": "mech_004",
                    "level": "fine",
                    "verdict": "surviving_exact",
                    "interval_gene_ids": "fine_gene",
                },
            ])
            row = frame.iloc[0]
            self.assertEqual("fine_gene", row["member_interval_genes"])
            self.assertEqual("interval_anchored", row["anchor_status"])

    def test_no_member_at_either_level_stays_unanchored(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            frame = self._run(Path(tmp), [
                {
                    "hypothesis_id": "mech_004",
                    "level": "coarse",
                    "verdict": "absent",
                    "interval_gene_ids": "",
                },
                {
                    "hypothesis_id": "mech_004",
                    "level": "fine",
                    "verdict": "absent_annotation_only",
                    "interval_gene_ids": "",
                },
            ])
            row = frame.iloc[0]
            self.assertEqual("", row["member_interval_genes"])
            self.assertEqual("0", row["member_count"])
            self.assertEqual("unanchored", row["anchor_status"])

    def test_not_computed_fine_keeps_coarse_member(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            frame = self._run(Path(tmp), [
                {
                    "hypothesis_id": "mech_004",
                    "level": "coarse",
                    "verdict": "surviving_enriched",
                    "interval_gene_ids": "coarse_gene",
                },
                {
                    "hypothesis_id": "mech_004",
                    "level": "fine",
                    "verdict": "not_computed_missing_provider",
                    "interval_gene_ids": "",
                },
            ])
            row = frame.iloc[0]
            self.assertEqual("coarse_gene", row["member_interval_genes"])
            self.assertEqual("interval_anchored", row["anchor_status"])


if __name__ == "__main__":
    unittest.main()
