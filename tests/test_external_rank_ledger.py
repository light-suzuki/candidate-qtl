from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main

NATIVE = (
    "gene_id\tdry_score\tscreening_rank\n"
    "GENE_A\t9.0\t1\n"
    "GENE_B\t8.0\t2\n"
    "GENE_C\t3.0\t3\n"
    "GENE_D\t2.0\t4\n"
    "GENE_E\t1.0\t5\n"
)
EXT_HEADER = "ranker_id\tranker_family\tgene_id\trank\traw_score\tpercentile\tcoverage\tmapping_confidence\treference_species\tsource\tversion\tdependency_group\n"
EXT_BODY = (
    "R1\tfam_comp\tGENE_A\t1\t0.9\t\tcovered\thigh\tanchor\tortho\t1\tg1\n"
    "R1\tfam_comp\tGENE_C\t2\t0.7\t\tcovered\tmed\tanchor\tortho\t1\tg1\n"
    "R3\tfam_comp\tGENE_A\t1\t0.9\t\tcovered\thigh\tanchor\tortho\t1\tg1\n"
    "R3\tfam_comp\tGENE_B\t3\t0.5\t\tcovered\tmed\tanchor\tortho\t1\tg1\n"
    "R2\tfam_path\tGENE_C\t1\t20\t\tcovered\tmed\tref\tpathway\t1\tg2\n"
    "R2\tfam_path\tGENE_B\t2\t15\t\tcovered\tmed\tref\tpathway\t1\tg2\n"
    "R2\tfam_path\tGENE_D\t5\t3\t\tcovered\tlow\tref\tpathway\t1\tg2\n"
    "R2\tfam_path\tGENE_E\t\t0\t\tuncovered\tlow\tref\tpathway\t1\tg2\n"
)


class ExternalRankLedgerTests(unittest.TestCase):
    def _fixture(self, root: Path, ext_body: str = EXT_BODY) -> tuple[Path, Path]:
        native = root / "ranked_candidates.tsv"; native.write_text(NATIVE, encoding="utf-8")
        external = root / "external_ranks.tsv"; external.write_text(EXT_HEADER + ext_body, encoding="utf-8")
        return native, external

    def test_all_five_classes_and_dependency_group_dedup(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); native, external = self._fixture(root); out = root / "out"
            self.assertEqual(0, main(["external-rank", "--native-rank", str(native), "--external-ranks", str(external),
                                      "--out-dir", str(out), "--top-k", "2"]))
            rows = pd.read_csv(out / "external_rank_ledger.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            self.assertEqual("native_only", rows.loc["GENE_A", "classification"])
            self.assertEqual("consensus_top", rows.loc["GENE_B", "classification"])
            self.assertEqual("external_rescue", rows.loc["GENE_C", "classification"])
            self.assertEqual("discordant", rows.loc["GENE_D", "classification"])
            self.assertEqual("uncovered", rows.loc["GENE_E", "classification"])
            # GENE_A is ranked 1 by two rankers in the same dependency group; the
            # RRF consensus must count that group only once (not 2/(60+1)).
            self.assertAlmostEqual(round(1.0 / 61, 8), float(rows.loc["GENE_A", "consensus_rrf"]), places=8)
            self.assertEqual("2", rows.loc["GENE_A", "covered_ranker_count"])
            manifest = json.loads((out / "external_rank_ledger_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(2, manifest["dependency_group_count"])
            self.assertEqual(3, manifest["external_ranker_count"])

    def test_native_dry_score_is_never_touched(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); native, external = self._fixture(root); out = root / "out"
            before = native.read_bytes()
            self.assertEqual(0, main(["external-rank", "--native-rank", str(native), "--external-ranks", str(external), "--out-dir", str(out)]))
            self.assertEqual(before, native.read_bytes())
            ledger_cols = set(pd.read_csv(out / "external_rank_ledger.tsv", sep="\t").columns)
            self.assertNotIn("dry_score", ledger_cols)

    def test_single_external_ranker_alone_imports(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            native, external = self._fixture(root, "R1\tfam_comp\tGENE_A\t1\t0.9\t\tcovered\thigh\tanchor\tortho\t1\tg1\n")
            out = root / "out"
            self.assertEqual(0, main(["external-rank", "--native-rank", str(native), "--external-ranks", str(external), "--out-dir", str(out), "--top-k", "2"]))
            rows = pd.read_csv(out / "external_rank_ledger.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            self.assertEqual("consensus_top", rows.loc["GENE_A", "classification"])
            self.assertEqual("uncovered", rows.loc["GENE_E", "classification"])

    def test_byte_stable_regeneration(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); native, external = self._fixture(root)
            a, b = root / "a", root / "b"
            for target in (a, b):
                self.assertEqual(0, main(["external-rank", "--native-rank", str(native), "--external-ranks", str(external), "--out-dir", str(target), "--top-k", "2"]))
            for name in ("external_rank_ledger.tsv", "coverage_report.tsv", "rank_discordance.tsv", "consensus_contributions.tsv"):
                self.assertEqual((a / name).read_bytes(), (b / name).read_bytes(), name)

    def test_unknown_gene_in_external_rank_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            native, external = self._fixture(root, "R1\tfam_comp\tNOT_A_CANDIDATE\t1\t0.9\t\tcovered\thigh\tanchor\tortho\t1\tg1\n")
            with self.assertRaises(ValueError):
                main(["external-rank", "--native-rank", str(native), "--external-ranks", str(external), "--out-dir", str(root / "out")])


if __name__ == "__main__":
    unittest.main()
