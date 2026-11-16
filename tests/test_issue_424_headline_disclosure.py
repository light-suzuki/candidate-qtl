"""Issue #424: disclose how the headline homology hit was chosen.

When the closest hit (best bitscore) and the most informative hit (first
qualifying hit whose annotation states a function) differ, the rationale shows
both with their selection-rule labels (最近縁 / 機能記載) so the reader can tell
which rule chose each line. The single-line form is kept when one hit plays
both roles.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl import homology_qc
from candidate_qtl.cli import main


GENOME = ">chr1\n" + "ATG" * 300 + "\n"
GFF = (
    "##gff-version 3\n"
    "chr1\tfixture\tgene\t100\t300\t.\t+\t.\tID=GeneA;Note=candidate protein\n"
    "chr1\tfixture\tmRNA\t100\t300\t.\t+\t.\tID=GeneA.1;Parent=GeneA\n"
    "chr1\tfixture\tCDS\t100\t300\t.\t+\t0\tParent=GeneA.1\n"
)

# pea_a2-style hit set: the closest hit (Medicago, identity 94.5%) is a
# placeholder name, while the informative Swiss-Prot TTG1 entry sits lower.
CLOSEST_WD40 = {
    "subject": "tr|A0A0C3VMP6|A0A0C3VMP6_MEDTR", "identity": 94.48,
    "query_coverage": 1.003, "subject_coverage": 1.006, "evalue": 0.0,
    "bitscore": 668.0, "description": "WD40-1 protein OS=Medicago truncatula PE=2",
}
INFORMATIVE_TTG1 = {
    "subject": "sp|Q9XGN1|TTG1_ARATH", "identity": 79.19,
    "query_coverage": 1.009, "subject_coverage": 1.015, "evalue": 0.0,
    "bitscore": 560.0,
    "description": "Protein TRANSPARENT TESTA GLABRA 1 OS=Arabidopsis thaliana OX=3702 GN=TTG1",
}
INFORMATIVE_CHS = {
    "subject": "sp|P13114|CHSY_ARATH", "identity": 85.0,
    "query_coverage": 0.95, "subject_coverage": 0.95, "evalue": 1e-50,
    "bitscore": 300.0, "description": "Chalcone synthase OS=Arabidopsis thaliana",
}


class HeadlineSelectionDisclosureTests(unittest.TestCase):
    """Unit-level checks on format_grounded_rationale (issue #424)."""

    def test_divergent_hits_show_both_with_selection_labels(self) -> None:
        hits = [CLOSEST_WD40, INFORMATIVE_TTG1]
        rationale = homology_qc.format_grounded_rationale(hits=hits, include_disclaimer=False)
        self.assertIn(homology_qc.HEADLINE_LABEL_CLOSEST, rationale)
        self.assertIn(homology_qc.HEADLINE_LABEL_INFORMATIVE, rationale)
        # 最近縁 (closest) line carries species, accession, identity, annotation
        self.assertIn("Medicago truncatula (A0A0C3VMP6)", rationale)
        self.assertIn("identity 94.5%", rationale)
        self.assertIn("annotated as WD40-1 protein", rationale)
        self.assertIn(homology_qc.NO_FUNCTION_NOTE, rationale)
        # 機能記載 (informative) line carries the same fields for the other hit
        self.assertIn("Arabidopsis thaliana (Q9XGN1)", rationale)
        self.assertIn("identity 79.2%", rationale)
        self.assertIn("annotated as Protein TRANSPARENT TESTA GLABRA 1", rationale)

    def test_single_qualifying_hit_keeps_legacy_form(self) -> None:
        rationale = homology_qc.format_grounded_rationale(hits=[INFORMATIVE_CHS], include_disclaimer=False)
        self.assertIn("homologous to Arabidopsis thaliana (P13114)", rationale)
        self.assertIn("identity 85.0%", rationale)
        self.assertNotIn(homology_qc.HEADLINE_LABEL_CLOSEST, rationale)
        self.assertNotIn(homology_qc.HEADLINE_LABEL_INFORMATIVE, rationale)

    def test_closest_is_also_informative_keeps_single_line(self) -> None:
        """When the closest hit already states a function, one line is enough."""
        second = {
            "subject": "sp|Q08401|CHSY_MEDTR", "identity": 80.0,
            "query_coverage": 0.90, "subject_coverage": 0.90, "evalue": 1e-45,
            "bitscore": 250.0, "description": "Chalcone synthase 1 OS=Medicago truncatula",
        }
        rationale = homology_qc.format_grounded_rationale(hits=[INFORMATIVE_CHS, second], include_disclaimer=False)
        self.assertIn("homologous to Arabidopsis thaliana (P13114)", rationale)
        self.assertNotIn(homology_qc.HEADLINE_LABEL_CLOSEST, rationale)
        self.assertNotIn(homology_qc.HEADLINE_LABEL_INFORMATIVE, rationale)

    def test_no_informative_hit_keeps_single_line_with_closest(self) -> None:
        hits = [
            CLOSEST_WD40,
            {"subject": "tr|B4XH97|B4XH97_MEDTR", "identity": 94.48,
             "query_coverage": 1.003, "subject_coverage": 1.006, "evalue": 0.0,
             "bitscore": 668.0,
             "description": "Transparent testa glabra protein OS=Medicago truncatula"},
        ]
        rationale = homology_qc.format_grounded_rationale(hits=hits, include_disclaimer=False)
        self.assertIn("homologous to Medicago truncatula (A0A0C3VMP6)", rationale)
        self.assertNotIn(homology_qc.HEADLINE_LABEL_CLOSEST, rationale)
        self.assertNotIn(homology_qc.HEADLINE_LABEL_INFORMATIVE, rationale)

    def test_two_line_form_coexists_with_the_discordance_note(self) -> None:
        """#424: when the closest hit is uninformative and the informative hits
        conflict, the two-line disclosure and the #422 note both appear."""
        hits = [
            CLOSEST_WD40,
            {"subject": "tr|G8A360|G8A360_MEDTR", "identity": 88.74,
             "query_coverage": 1.0, "subject_coverage": 1.0, "evalue": 0.0,
             "bitscore": 515.0,
             "description": "flavonoid 3',5'-hydroxylase OS=Medicago truncatula GN=MTR_3g436540"},
            {"subject": "sp|Q9SD85|F3PH_ARATH", "identity": 48.71,
             "query_coverage": 0.977, "subject_coverage": 0.981, "evalue": 2.7e-179,
             "bitscore": 495.0,
             "description": "Flavonoid 3'-monooxygenase OS=Arabidopsis thaliana GN=CYP75B1"},
        ]
        state, _ = homology_qc.assess_annotation_concordance(hits)
        self.assertEqual(state, "discordant")
        rationale = homology_qc.format_grounded_rationale(hits=hits, include_disclaimer=False)
        self.assertIn(homology_qc.HEADLINE_LABEL_CLOSEST, rationale)
        self.assertIn(homology_qc.HEADLINE_LABEL_INFORMATIVE, rationale)
        self.assertIn("不一致", rationale)
        self.assertIn("monooxygenase", rationale)

    def test_disclaimer_and_cautious_function_phrasing_are_kept(self) -> None:
        rationale = homology_qc.format_grounded_rationale(
            hits=[CLOSEST_WD40, INFORMATIVE_TTG1], include_disclaimer=True)
        self.assertIn(homology_qc.HOMOLOGY_DISCLAIMER, rationale)
        self.assertNotIn("機能を持つ", rationale)
        self.assertNotIn("は機能がある", rationale)


class ReportAndRankingIntegrationTests(unittest.TestCase):
    """Report.md rendering and dry_score invariance through discover (issue #424)."""

    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        (self.root / "genome.fa").write_text(GENOME, encoding="utf-8")
        (self.root / "genes.gff3").write_text(GFF, encoding="utf-8")

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def _run_discover(self, hits: dict, sub_dir: str) -> tuple[Path, pd.DataFrame]:
        out = self.root / sub_dir
        args = [
            "discover",
            "--fasta", str(self.root / "genome.fa"),
            "--gff3", str(self.root / "genes.gff3"),
            "--interval", "chr1:1-900",
            "--phenotype", "xyzzy",
            "--out-dir", str(out),
            "--no-auto-local-db",
            "--allow-external-api",
        ]
        with patch("candidate_qtl.discover._best_ncbi_hits", return_value=(hits, "", True, False)):
            self.assertEqual(0, main(args))
        ranked = pd.read_csv(out / "05_rank" / "ranked_candidates.tsv", sep="\t").fillna("")
        return out, ranked

    def test_report_md_shows_both_lines_when_headline_hits_differ(self) -> None:
        out, _ = self._run_discover({"GeneA": [CLOSEST_WD40, INFORMATIVE_TTG1]}, "out_report")
        report = (out / "report.md").read_text(encoding="utf-8")
        self.assertIn("最近縁: Medicago truncatula (A0A0C3VMP6)", report)
        self.assertIn("機能記載: Arabidopsis thaliana (Q9XGN1)", report)
        self.assertIn("identity 94.5%", report)
        self.assertIn("identity 79.2%", report)
        self.assertIn("※ 全長相同性のみで判定。活性部位残基の保存は未確認。", report)
        manifest = json.loads((out / "discover_manifest.json").read_text(encoding="utf-8"))
        self.assertIn("single_informative", manifest["annotation_concordance_counts"])

    def test_ranking_and_dry_score_are_unchanged_by_the_disclosure(self) -> None:
        _, ranked_divergent = self._run_discover(
            {"GeneA": [CLOSEST_WD40, INFORMATIVE_TTG1]}, "out_divergent")
        _, ranked_single = self._run_discover({"GeneA": [INFORMATIVE_TTG1]}, "out_single")
        for column in ("dry_score", "dry_rank", "dry_group", "screening_rank"):
            self.assertEqual(
                str(ranked_divergent.loc[ranked_divergent["gene_id"] == "GeneA", column].iloc[0]),
                str(ranked_single.loc[ranked_single["gene_id"] == "GeneA", column].iloc[0]),
                f"{column} changed with the headline disclosure",
            )

    def test_evidence_tsv_does_not_carry_the_rationale_text(self) -> None:
        out, _ = self._run_discover({"GeneA": [CLOSEST_WD40, INFORMATIVE_TTG1]}, "out_evidence")
        evidence = pd.read_csv(out / "04_evidence.tsv", sep="\t").fillna("")
        text = evidence.to_csv(sep="\t", index=False)
        self.assertNotIn("最近縁", text)
        self.assertNotIn("機能記載", text)


if __name__ == "__main__":
    unittest.main()
