"""Regression coverage for bounded annotation-aware reference selection (#647)."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.family_index import family_index
from candidate_qtl.reference_hit_selection import select_reference_hit


class ReferenceHitSelectionTests(unittest.TestCase):
    def _record(
        self,
        subject: str,
        *,
        bitscore: float,
        identity: float = 80.0,
        query_coverage: float = 0.90,
        subject_coverage: float = 0.90,
        title: str = "Protein kinase",
        usable_symbol: bool = True,
    ) -> dict[str, object]:
        return {
            "subject_id": subject,
            "bitscore": bitscore,
            "evalue": 1e-50,
            "identity": identity,
            "merged_query_coverage": query_coverage,
            "subject_coverage": subject_coverage,
            "subject_title": title,
            "usable_symbol": usable_symbol,
        }

    def test_reviewed_record_wins_only_inside_fixed_equivalence_window(self) -> None:
        records = [
            self._record("tr|TOP|TOP_SPEC", bitscore=100.0),
            self._record("sp|CURATED|CURATED_SPEC", bitscore=97.0, identity=80.5),
        ]
        result = select_reference_hit(records)
        self.assertEqual(1, result["selected_index"])
        self.assertEqual("reviewed_record_within_homology_equivalence_window", result["selection_reason"])
        self.assertEqual("tr|TOP|TOP_SPEC", result["sequence_top"]["subject_id"])

    def test_informative_title_breaks_an_exact_sequence_tie(self) -> None:
        records = [
            self._record("tr|A|A_SPEC", bitscore=100.0, title="Uncharacterized protein"),
            self._record("tr|B|B_SPEC", bitscore=100.0, title="Bowman Birk trypsin inhibitor"),
        ]
        result = select_reference_hit(records)
        self.assertEqual("tr|B|B_SPEC", result["selected"]["subject_id"])
        self.assertEqual("informative_title_within_homology_equivalence_window", result["selection_reason"])

    def test_materially_weaker_curated_record_cannot_displace_sequence_top(self) -> None:
        records = [
            self._record("tr|TOP|TOP_SPEC", bitscore=100.0),
            self._record("sp|WEAK|WEAK_SPEC", bitscore=80.0, identity=90.0),
        ]
        result = select_reference_hit(records)
        self.assertEqual(0, result["selected_index"])
        self.assertEqual("sequence_top", result["selection_reason"])

    def test_input_order_does_not_change_selected_subject(self) -> None:
        records = [
            self._record("tr|TOP|TOP_SPEC", bitscore=100.0),
            self._record("sp|CURATED|CURATED_SPEC", bitscore=97.0, identity=80.5),
        ]
        forward = select_reference_hit(records)
        reverse = select_reference_hit(list(reversed(records)))
        self.assertEqual(forward["selected"]["subject_id"], reverse["selected"]["subject_id"])
        self.assertEqual(forward["sequence_top"]["subject_id"], reverse["sequence_top"]["subject_id"])


class FamilyIndexSelectionIntegrationTests(unittest.TestCase):
    def _write_hits(self, root: Path, rows: list[list[str]]) -> Path:
        blast = root / "hits.tsv"
        blast.write_text(
            "\n".join("\t".join(row) for row in rows) + "\n",
            encoding="utf-8",
        )
        return blast

    def test_index_preserves_sequence_top_and_selected_reference(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            blast = self._write_hits(
                root,
                [
                    [
                        "q1", "tr|G7K6U2|CYCLOPS_MEDTR", "88.654", "300", "320", "320",
                        "1", "300", "1", "300", "1e-80", "929",
                        "Cyclops protein, putative GN=CYCLOPS",
                    ],
                    [
                        "q1", "sp|A7TUE1|CYCLOPS_MEDTR", "89.689", "300", "320", "320",
                        "1", "300", "1", "300", "1e-80", "914",
                        "Protein CYCLOPS GN=CYCLOPS",
                    ],
                ],
            )
            out = root / "out"
            result = family_index(str(blast), str(out))
            index = pd.read_csv(result["family_index"], sep="\t", dtype=str).fillna("")
            self.assertEqual(1, len(index))
            row = index.iloc[0]
            self.assertEqual("sp|A7TUE1|CYCLOPS_MEDTR", row["best_subject"])
            self.assertEqual("tr|G7K6U2|CYCLOPS_MEDTR", row["sequence_top_subject"])
            self.assertEqual("sp|A7TUE1|CYCLOPS_MEDTR", row["selected_reference_subject"])
            self.assertEqual("2", row["selected_reference_rank"])
            self.assertEqual("reviewed_record_within_homology_equivalence_window", row["reference_selection_reason"])

            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            policy = manifest["reference_hit_selection"]
            self.assertTrue(policy["sequence_first"])
            self.assertFalse(policy["annotation_is_causal_evidence"])
            self.assertFalse(policy["legacy_half_bitscore_symbol_promotion_active"])
            self.assertEqual(1, policy["selected_non_sequence_top_count"])

    def test_full_index_never_falls_back_to_legacy_half_bitscore_promotion(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            blast = self._write_hits(
                root,
                [
                    [
                        "q1", "tr|A0A072UD83|A0A072UD83_MEDTR", "80", "300", "320", "320",
                        "1", "300", "1", "300", "1e-80", "500",
                        "non-specific serine/threonine protein kinase OS=Medicago truncatula GN=MTR_5g034210",
                    ],
                    [
                        "q1", "sp|Q8L4H4|NORK_MEDTR", "90", "300", "320", "320",
                        "1", "300", "1", "300", "1e-80", "400",
                        "Nodulation receptor kinase OS=Medicago truncatula GN=NORK",
                    ],
                ],
            )
            result = family_index(str(blast), str(root / "out"))
            row = pd.read_csv(result["family_index"], sep="\t", dtype=str).fillna("").iloc[0]
            self.assertEqual("tr|A0A072UD83|A0A072UD83_MEDTR", row["best_subject"])
            self.assertEqual("tr|A0A072UD83|A0A072UD83_MEDTR", row["selected_reference_subject"])
            self.assertEqual("sequence_top", row["reference_selection_reason"])
            self.assertEqual("1.000000", row["selected_bitscore_ratio"])


if __name__ == "__main__":
    unittest.main()
