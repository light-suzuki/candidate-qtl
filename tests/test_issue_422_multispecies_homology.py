from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl import homology_qc
from candidate_qtl.ranking_policy import load_weight_registry


GENOME = ">chr1\n" + "ATG" * 300 + "\n"
GFF = (
    "##gff-version 3\n"
    "chr1\tfixture\tgene\t100\t300\t.\t+\t.\tID=GeneA;Note=candidate protein\n"
    "chr1\tfixture\tmRNA\t100\t300\t.\t+\t.\tID=GeneA.1;Parent=GeneA\n"
    "chr1\tfixture\tCDS\t100\t300\t.\t+\t0\tParent=GeneA.1\n"
)


class MultiSpeciesHomologyParsingTests(unittest.TestCase):
    """Unit tests for individual helper functions in homology_qc."""

    def test_extract_species(self) -> None:
        self.assertEqual(
            homology_qc.extract_species("Chalcone synthase OS=Arabidopsis thaliana OX=3702 GN=CHS PE=1 SV=1"),
            "Arabidopsis thaliana",
        )
        self.assertEqual(
            homology_qc.extract_species("Protein kinase OS=Oryza sativa subsp. japonica GN=KIN1"),
            "Oryza sativa subsp. japonica",
        )
        self.assertEqual(homology_qc.extract_species("No species annotation"), "")

    def test_is_reviewed_entry(self) -> None:
        self.assertTrue(homology_qc.is_reviewed_entry("sp|P13114|CHSY_ARATH"))
        self.assertTrue(homology_qc.is_reviewed_entry("sp|Q9LZS3|"))
        self.assertFalse(homology_qc.is_reviewed_entry("tr|A0A1I9LSP1|A0A1I9LSP1_ARATH"))
        self.assertFalse(homology_qc.is_reviewed_entry("NP_001321.1"))

    def test_extract_accession(self) -> None:
        self.assertEqual(homology_qc.extract_accession("sp|P13114|CHSY_ARATH"), "P13114")
        self.assertEqual(homology_qc.extract_accession("tr|A0A1I9LSP1|A0A1I9LSP1_ARATH"), "A0A1I9LSP1")
        self.assertEqual(homology_qc.extract_accession("UniProt:P13114"), "P13114")
        self.assertEqual(homology_qc.extract_accession("TAIR:AT1G69120"), "AT1G69120")
        self.assertEqual(homology_qc.extract_accession("P13114"), "P13114")

    def test_extract_protein_name(self) -> None:
        self.assertEqual(
            homology_qc.extract_protein_name(
                "Chalcone synthase OS=Arabidopsis thaliana OX=3702 GN=CHS PE=1 SV=1",
                "sp|P13114|CHSY_ARATH",
            ),
            "Chalcone synthase",
        )
        self.assertEqual(
            homology_qc.extract_protein_name(
                "sp|P13114|CHSY_ARATH Chalcone synthase 1 OS=Arabidopsis thaliana GN=CHS1",
            ),
            "Chalcone synthase 1",
        )

    def test_is_informative_annotation(self) -> None:
        self.assertTrue(homology_qc.is_informative_annotation("Chalcone synthase"))
        self.assertTrue(homology_qc.is_informative_annotation("Receptor-like serine/threonine kinase"))
        self.assertFalse(homology_qc.is_informative_annotation("Uncharacterized protein"))
        self.assertFalse(homology_qc.is_informative_annotation("hypothetical protein"))
        self.assertFalse(homology_qc.is_informative_annotation("putative uncharacterized protein At1g01010"))
        self.assertFalse(homology_qc.is_informative_annotation("protein of unknown function"))
        self.assertFalse(homology_qc.is_informative_annotation(""))

    def test_wd40_placeholder_is_uninformative_but_ttg1_swissprot_is_informative(self) -> None:
        """#422: 'WD40-1 protein' states no function; the TTG1 Swiss-Prot entry does."""
        self.assertFalse(homology_qc.is_informative_annotation("WD40-1 protein"))
        self.assertFalse(
            homology_qc.is_informative_annotation("WD40-1 protein", "tr|A0A0C3VMP6|A0A0C3VMP6_MEDTR")
        )
        self.assertTrue(
            homology_qc.is_informative_annotation(
                "Protein TRANSPARENT TESTA GLABRA 1", "sp|Q9XGN1|TTG1_ARATH"
            )
        )

    def test_castor_bare_name_uninformative_ion_channel_informative(self) -> None:
        """#422: bare 'CASTOR' is uninformative; 'Ion channel castor' states function."""
        self.assertFalse(homology_qc.is_informative_annotation("CASTOR"))
        self.assertFalse(homology_qc.is_informative_annotation("CASTOR", "tr|A0A072U4J9|A0A072U4J9_MEDTR"))
        self.assertTrue(homology_qc.is_informative_annotation("Ion channel castor"))
        self.assertTrue(homology_qc.is_informative_annotation("Ion channel CASTOR"))


class Dimension1MultiHitPreservationTests(unittest.TestCase):
    """Dimension 1: Multi-hit preservation in 04_evidence.tsv and homology_audit.tsv."""

    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        (self.root / "genome.fa").write_text(GENOME, encoding="utf-8")
        (self.root / "genes.gff3").write_text(GFF, encoding="utf-8")

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def _run(self, hits: dict) -> tuple[pd.DataFrame, Path]:
        out = self.root / f"out_{abs(hash(str(hits))) % 100000}"
        args = [
            "discover",
            "--fasta", str(self.root / "genome.fa"),
            "--gff3", str(self.root / "genes.gff3"),
            "--interval", "chr1:1-900",
            "--phenotype", "pigment",
            "--out-dir", str(out),
            "--no-auto-local-db",
            "--allow-external-api",
        ]
        with patch("candidate_qtl.discover._best_ncbi_hits", return_value=(hits, "", True, False)):
            self.assertEqual(0, main(args))
        return pd.read_csv(out / "04_evidence.tsv", sep="\t").fillna(""), out

    def test_multi_qualifying_hits_all_lifted_as_homology_context(self) -> None:
        hits = {
            "GeneA": [
                {"subject": "sp|P13114|CHSY_ARATH", "identity": 85.0, "query_coverage": 0.95, "subject_coverage": 0.95, "evalue": 1e-50, "bitscore": 300.0, "description": "Chalcone synthase OS=Arabidopsis thaliana GN=CHS"},
                {"subject": "sp|Q08401|CHSY_MEDTR", "identity": 80.0, "query_coverage": 0.92, "subject_coverage": 0.90, "evalue": 1e-45, "bitscore": 250.0, "description": "Chalcone synthase 1 OS=Medicago truncatula GN=CHS1"},
                {"subject": "tr|A0A1I9LSP1|A0A1I9LSP1_ORYSJ", "identity": 75.0, "query_coverage": 0.90, "subject_coverage": 0.88, "evalue": 1e-40, "bitscore": 200.0, "description": "Chalcone synthase OS=Oryza sativa GN=CHS"},
            ]
        }
        evidence, out = self._run(hits)
        homology_rows = evidence[evidence.evidence_code == "homology_context"]
        self.assertEqual(len(homology_rows), 3)
        self.assertTrue((homology_rows["direction"] == "context").all())

        # Check enriched homology_audit.tsv
        audit = pd.read_csv(out / "03_annotation" / "homology_audit.tsv", sep="\t").fillna("")
        self.assertEqual(len(audit), 3)
        self.assertIn("species", audit.columns)
        self.assertIn("accession", audit.columns)
        self.assertIn("reviewed", audit.columns)
        self.assertIn("protein_name", audit.columns)
        self.assertIn("is_informative", audit.columns)

        self.assertEqual(list(audit["species"]), ["Arabidopsis thaliana", "Medicago truncatula", "Oryza sativa"])
        self.assertEqual(list(audit["accession"]), ["P13114", "Q08401", "A0A1I9LSP1"])
        self.assertEqual(list(audit["reviewed"]), [True, True, False])


class Dimension2FourConcordanceStatesTests(unittest.TestCase):
    """Dimension 2: Correct classification of 4 concordance states."""

    def test_state_concordant(self) -> None:
        hits = [
            {"subject": "sp|P13114|CHSY_ARATH", "identity": 85.0, "query_coverage": 0.95, "evalue": 1e-50, "bitscore": 300.0, "description": "Chalcone synthase OS=Arabidopsis thaliana"},
            {"subject": "sp|Q08401|CHSY_MEDTR", "identity": 80.0, "query_coverage": 0.90, "evalue": 1e-40, "bitscore": 250.0, "description": "Chalcone synthase 1 OS=Medicago truncatula"},
        ]
        state, _ = homology_qc.assess_annotation_concordance(hits)
        self.assertEqual(state, "concordant")

    def test_state_discordant(self) -> None:
        hits = [
            {"subject": "sp|P13114|CHSY_ARATH", "identity": 85.0, "query_coverage": 0.95, "evalue": 1e-50, "bitscore": 300.0, "description": "Chalcone synthase OS=Arabidopsis thaliana"},
            {"subject": "sp|Q9LKA2|KIN1_ORYSJ", "identity": 80.0, "query_coverage": 0.90, "evalue": 1e-40, "bitscore": 250.0, "description": "Receptor-like serine/threonine kinase OS=Oryza sativa"},
        ]
        state, _ = homology_qc.assess_annotation_concordance(hits)
        self.assertEqual(state, "discordant")

    def test_state_single_informative(self) -> None:
        hits = [
            {"subject": "sp|P13114|CHSY_ARATH", "identity": 85.0, "query_coverage": 0.95, "evalue": 1e-50, "bitscore": 300.0, "description": "Chalcone synthase OS=Arabidopsis thaliana"},
            {"subject": "tr|A0A123|A0A123_ORYSJ", "identity": 80.0, "query_coverage": 0.90, "evalue": 1e-40, "bitscore": 250.0, "description": "Uncharacterized protein OS=Oryza sativa"},
        ]
        state, _ = homology_qc.assess_annotation_concordance(hits)
        self.assertEqual(state, "single_informative")

    def test_state_all_uninformative(self) -> None:
        hits = [
            {"subject": "tr|A0A111|A0A111_ARATH", "identity": 85.0, "query_coverage": 0.95, "evalue": 1e-50, "bitscore": 300.0, "description": "Uncharacterized protein OS=Arabidopsis thaliana"},
            {"subject": "tr|A0A123|A0A123_ORYSJ", "identity": 80.0, "query_coverage": 0.90, "evalue": 1e-40, "bitscore": 250.0, "description": "Hypothetical protein OS=Oryza sativa"},
        ]
        state, _ = homology_qc.assess_annotation_concordance(hits)
        self.assertEqual(state, "all_uninformative")

        state_empty, _ = homology_qc.assess_annotation_concordance([])
        self.assertEqual(state_empty, "all_uninformative")


class Dimension3ZeroScoreChangeTests(unittest.TestCase):
    """Dimension 3: Multi-species agreement is never converted into score points."""

    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        (self.root / "genome.fa").write_text(GENOME, encoding="utf-8")
        (self.root / "genes.gff3").write_text(GFF, encoding="utf-8")

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def _run_discover(self, hits: dict, sub_dir: str) -> pd.DataFrame:
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
        return pd.read_csv(out / "05_rank" / "ranked_candidates.tsv", sep="\t").fillna("")

    def test_single_hit_vs_five_concordant_hits_yields_identical_dry_score(self) -> None:
        single_hit = {
            "GeneA": [
                {"subject": "sp|P13114|CHSY_ARATH", "identity": 85.0, "query_coverage": 0.95, "subject_coverage": 0.95, "evalue": 1e-50, "bitscore": 300.0, "description": "Chalcone synthase OS=Arabidopsis thaliana"}
            ]
        }
        five_concordant_hits = {
            "GeneA": [
                {"subject": "sp|P13114|CHSY_ARATH", "identity": 85.0, "query_coverage": 0.95, "subject_coverage": 0.95, "evalue": 1e-50, "bitscore": 300.0, "description": "Chalcone synthase OS=Arabidopsis thaliana"},
                {"subject": "sp|Q08401|CHSY_MEDTR", "identity": 80.0, "query_coverage": 0.90, "subject_coverage": 0.90, "evalue": 1e-45, "bitscore": 250.0, "description": "Chalcone synthase 1 OS=Medicago truncatula"},
                {"subject": "sp|P10290|CHSY_SOLLC", "identity": 75.0, "query_coverage": 0.85, "subject_coverage": 0.85, "evalue": 1e-40, "bitscore": 200.0, "description": "Chalcone synthase 1 OS=Solanum lycopersicum"},
                {"subject": "sp|Q43210|CHSY_ORYSJ", "identity": 70.0, "query_coverage": 0.80, "subject_coverage": 0.80, "evalue": 1e-35, "bitscore": 180.0, "description": "Chalcone synthase OS=Oryza sativa"},
                {"subject": "sp|P24825|CHSY_ZEAMA", "identity": 65.0, "query_coverage": 0.75, "subject_coverage": 0.75, "evalue": 1e-30, "bitscore": 150.0, "description": "Chalcone synthase C2 OS=Zea mays"},
            ]
        }
        ranked_single = self._run_discover(single_hit, "out_single")
        ranked_multi = self._run_discover(five_concordant_hits, "out_multi")

        score_single = float(ranked_single.loc[ranked_single["gene_id"] == "GeneA", "dry_score"].iloc[0])
        score_multi = float(ranked_multi.loc[ranked_multi["gene_id"] == "GeneA", "dry_score"].iloc[0])
        self.assertEqual(score_single, score_multi)

    def test_nonzero_score_unchanged_when_hits_multiplied(self) -> None:
        """#422: with a phenotype that actually scores (trait support present),
        adding more species hits must not raise the dry_score."""
        single_hit = {
            "GeneA": [
                {"subject": "sp|P13114|CHSY_ARATH", "identity": 85.0, "query_coverage": 0.95, "subject_coverage": 0.95, "evalue": 1e-50, "bitscore": 300.0, "description": "Chalcone synthase OS=Arabidopsis thaliana"},
            ]
        }
        multi_hit = {
            "GeneA": [
                {"subject": "sp|P13114|CHSY_ARATH", "identity": 85.0, "query_coverage": 0.95, "subject_coverage": 0.95, "evalue": 1e-50, "bitscore": 300.0, "description": "Chalcone synthase OS=Arabidopsis thaliana"},
                {"subject": "sp|Q08401|CHSY_MEDTR", "identity": 80.0, "query_coverage": 0.90, "subject_coverage": 0.90, "evalue": 1e-45, "bitscore": 250.0, "description": "Chalcone synthase 1 OS=Medicago truncatula"},
                {"subject": "sp|P10290|CHSY_SOLLC", "identity": 75.0, "query_coverage": 0.85, "subject_coverage": 0.85, "evalue": 1e-40, "bitscore": 200.0, "description": "Chalcone synthase 1 OS=Solanum lycopersicum"},
            ]
        }

        def run(hits: dict, sub_dir: str) -> pd.DataFrame:
            out = self.root / sub_dir
            args = [
                "discover",
                "--fasta", str(self.root / "genome.fa"),
                "--gff3", str(self.root / "genes.gff3"),
                "--interval", "chr1:1-900",
                "--phenotype", "chalcone synthase",
                "--out-dir", str(out),
                "--no-auto-local-db",
                "--allow-external-api",
            ]
            with patch("candidate_qtl.discover._best_ncbi_hits", return_value=(hits, "", True, False)):
                self.assertEqual(0, main(args))
            return pd.read_csv(out / "05_rank" / "ranked_candidates.tsv", sep="\t").fillna("")

        ranked_single = run(single_hit, "out_scored_single")
        ranked_multi = run(multi_hit, "out_scored_multi")
        score_single = float(ranked_single.loc[ranked_single["gene_id"] == "GeneA", "dry_score"].iloc[0])
        score_multi = float(ranked_multi.loc[ranked_multi["gene_id"] == "GeneA", "dry_score"].iloc[0])
        self.assertGreater(score_single, 0.0, "fixture should produce a non-zero dry_score")
        self.assertEqual(score_single, score_multi)
        evidence = pd.read_csv(self.root / "out_scored_multi" / "04_evidence.tsv", sep="\t").fillna("")
        ctx = evidence[evidence["evidence_code"] == "homology_context"]
        self.assertEqual(len(ctx), 3)
        self.assertTrue((ctx["direction"] == "context").all())


class Dimension4NonResolutionOfDiscordanceTests(unittest.TestCase):
    """Dimension 4: Discordance is never resolved by taxonomy proximity."""

    def test_discordance_preserved_regardless_of_species_taxa(self) -> None:
        hits = [
            {"subject": "sp|P13114|CHSY_ARATH", "identity": 88.0, "query_coverage": 0.95, "evalue": 1e-50, "bitscore": 300.0, "description": "Chalcone synthase OS=Arabidopsis thaliana"},
            {"subject": "sp|Q9LKA2|KIN1_ORYSJ", "identity": 80.0, "query_coverage": 0.90, "evalue": 1e-40, "bitscore": 250.0, "description": "Receptor serine/threonine kinase OS=Oryza sativa"},
        ]
        state, details = homology_qc.assess_annotation_concordance(hits)
        self.assertEqual(state, "discordant")
        self.assertIn("discordant", details)

    def test_pea_b_case_enzyme_class_conflict_is_discordant(self) -> None:
        """#422: F3'5'H (hydroxylase) vs F3'H (monooxygenase) share 'flavonoid'
        but are different enzymes; the class-word conflict must be discordant."""
        hits = [
            {"subject": "tr|G8A360|G8A360_MEDTR", "identity": 88.74, "query_coverage": 1.0, "subject_coverage": 1.0, "evalue": 0.0, "bitscore": 948.0, "description": "flavonoid 3',5'-hydroxylase OS=Medicago truncatula GN=MTR_3g436540"},
            {"subject": "LotjaGi6g1v0187300.1", "identity": 79.09, "query_coverage": 0.984, "subject_coverage": 0.988, "evalue": 0.0, "bitscore": 862.0, "description": "LotjaGi6g1v0187300.1"},
            {"subject": "tr|A0A072UUU7|A0A072UUU7_MEDTR", "identity": 62.08, "query_coverage": 0.988, "subject_coverage": 1.006, "evalue": 0.0, "bitscore": 666.0, "description": "Flavonoid hydroxylase OS=Medicago truncatula GN=MTR_3g436390"},
            {"subject": "sp|Q9SD85|F3PH_ARATH", "identity": 48.71, "query_coverage": 0.977, "subject_coverage": 0.981, "evalue": 2.7e-179, "bitscore": 515.0, "description": "Flavonoid 3'-monooxygenase OS=Arabidopsis thaliana GN=CYP75B1"},
            {"subject": "tr|G7IWE5|G7IWE5_MEDTR", "identity": 49.58, "query_coverage": 0.917, "subject_coverage": 0.915, "evalue": 2.4e-171, "bitscore": 495.0, "description": "Flavonoid 3'-monooxygenase OS=Medicago truncatula GN=MTR_3g025230"},
        ]
        state, details = homology_qc.assess_annotation_concordance(hits)
        self.assertEqual(state, "discordant", details)

    def test_same_species_conflict_stays_discordant(self) -> None:
        """Taxonomy is never used to resolve: two Medicago hits with different
        enzyme classes are discordant just like cross-species pairs."""
        hits = [
            {"subject": "tr|A0A072X1|A0A072X1_MEDTR", "identity": 85.0, "query_coverage": 0.95, "evalue": 1e-50, "bitscore": 300.0, "description": "Beta-amyrin synthase OS=Medicago truncatula"},
            {"subject": "tr|A0A072X2|A0A072X2_MEDTR", "identity": 80.0, "query_coverage": 0.90, "evalue": 1e-40, "bitscore": 250.0, "description": "Terpene cyclase/mutase family member OS=Medicago truncatula"},
        ]
        state, details = homology_qc.assess_annotation_concordance(hits)
        self.assertEqual(state, "discordant")
        self.assertIn("discordant", details)

    def test_discordant_hit_set_has_no_taxonomy_resolution_in_rationale(self) -> None:
        """The grounded rationale for a discordant set notes the conflict but
        never picks a winner by species proximity."""
        hits = [
            {"subject": "tr|G8A360|G8A360_MEDTR", "identity": 88.74, "query_coverage": 1.0, "subject_coverage": 1.0, "evalue": 0.0, "bitscore": 948.0, "description": "flavonoid 3',5'-hydroxylase OS=Medicago truncatula GN=MTR_3g436540"},
            {"subject": "sp|Q9SD85|F3PH_ARATH", "identity": 48.71, "query_coverage": 0.977, "subject_coverage": 0.981, "evalue": 2.7e-179, "bitscore": 515.0, "description": "Flavonoid 3'-monooxygenase OS=Arabidopsis thaliana GN=CYP75B1"},
        ]
        rationale = homology_qc.format_grounded_rationale(hits=hits, include_disclaimer=False)
        self.assertIn("homologous to Medicago truncatula (G8A360)", rationale)
        self.assertIn("annotated as flavonoid 3',5'-hydroxylase", rationale)
        self.assertIn("不一致", rationale)
        self.assertIn("monooxygenase", rationale)
        self.assertNotIn("is a", rationale.casefold())
        self.assertNotIn("Medicago が正しい", rationale)
        self.assertNotIn("近い", rationale)


class Dimension5GroundingAndPhrasingTests(unittest.TestCase):
    """Dimension 5: Strict phrasing, grounding in hit data, and mandatory disclaimer."""

    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        (self.root / "genome.fa").write_text(GENOME, encoding="utf-8")
        (self.root / "genes.gff3").write_text(GFF, encoding="utf-8")

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_rationale_strict_phrasing_and_disclaimer(self) -> None:
        rationale = homology_qc.format_grounded_rationale("Arabidopsis thaliana", "Chalcone synthase", "P13114")
        self.assertIn("homologous to Arabidopsis thaliana (P13114), annotated as Chalcone synthase", rationale)
        self.assertIn("※ 全長相同性のみで判定。活性部位残基の保存は未確認。", rationale)

    def test_report_md_contains_top_candidate_homology_evidence_and_disclaimer(self) -> None:
        out = self.root / "out_report"
        hits = {
            "GeneA": [
                {"subject": "sp|P13114|CHSY_ARATH", "identity": 85.0, "query_coverage": 0.95, "subject_coverage": 0.95, "evalue": 1e-50, "bitscore": 300.0, "description": "Chalcone synthase OS=Arabidopsis thaliana"}
            ]
        }
        args = [
            "discover",
            "--fasta", str(self.root / "genome.fa"),
            "--gff3", str(self.root / "genes.gff3"),
            "--interval", "chr1:1-900",
            "--phenotype", "chalcone",
            "--out-dir", str(out),
            "--no-auto-local-db",
            "--allow-external-api",
        ]
        with patch("candidate_qtl.discover._best_ncbi_hits", return_value=(hits, "", True, False)):
            self.assertEqual(0, main(args))

        report_text = (out / "report.md").read_text(encoding="utf-8")
        self.assertIn("## Top Candidate Homology Evidence", report_text)
        self.assertIn("※ 全長相同性のみで判定。活性部位残基の保存は未確認。", report_text)
        self.assertIn(
            "homologous to Arabidopsis thaliana (P13114), identity 85.0% / "
            "query coverage 0.95 / subject coverage 0.95, "
            "annotated as Chalcone synthase",
            report_text,
        )

        # Check annotation_concordance.tsv exists and is referenced in manifest
        self.assertTrue((out / "03_annotation" / "annotation_concordance.tsv").exists())
        manifest = json.loads((out / "discover_manifest.json").read_text(encoding="utf-8"))
        self.assertIn("annotation_concordance", manifest)
        self.assertIn("annotation_concordance_counts", manifest)

    def test_rationale_includes_identity_and_coverage_from_hit_row(self) -> None:
        """#422: rationale must carry species, accession, identity, coverage and
        the annotation, all taken from the actual hit row."""
        hits = [
            {"subject": "sp|P13114|CHSY_ARATH", "identity": 85.0, "query_coverage": 0.95, "subject_coverage": 0.93, "evalue": 1e-50, "bitscore": 300.0, "description": "Chalcone synthase OS=Arabidopsis thaliana"},
        ]
        rationale = homology_qc.format_grounded_rationale(hits=hits, include_disclaimer=False)
        self.assertIn("Arabidopsis thaliana (P13114)", rationale)
        self.assertIn("identity 85.0%", rationale)
        self.assertIn("query coverage 0.95", rationale)
        self.assertIn("subject coverage 0.93", rationale)
        self.assertIn("annotated as Chalcone synthase", rationale)

    def test_rationale_uses_no_words_not_supplied_by_hit_rows(self) -> None:
        """#422: the rationale must not introduce vocabulary absent from the hit
        rows (e.g. a guessed function name)."""
        hits = [
            {"subject": "tr|A0A0C3VMP6|A0A0C3VMP6_MEDTR", "identity": 94.48, "query_coverage": 1.003, "subject_coverage": 1.006, "evalue": 0.0, "bitscore": 668.0, "description": "WD40-1 protein OS=Medicago truncatula PE=2"},
            {"subject": "sp|Q9XGN1|TTG1_ARATH", "identity": 79.19, "query_coverage": 1.009, "subject_coverage": 1.015, "evalue": 0.0, "bitscore": 560.0, "description": "Protein TRANSPARENT TESTA GLABRA 1 OS=Arabidopsis thaliana GN=TTG1"},
        ]
        rationale = homology_qc.format_grounded_rationale(hits=hits, include_disclaimer=False)
        # The informative Swiss-Prot hit is the headline; the closest hit
        # (WD40-1, higher identity) is disclosed with its selection-rule
        # label instead of being silently hidden.
        self.assertIn("Q9XGN1", rationale)
        self.assertIn("TRANSPARENT TESTA GLABRA 1", rationale)
        self.assertIn("WD40-1", rationale)
        self.assertIn("A0A0C3VMP6", rationale)
        self.assertIn(homology_qc.HEADLINE_LABEL_CLOSEST, rationale)
        self.assertIn(homology_qc.HEADLINE_LABEL_INFORMATIVE, rationale)
        self.assertNotIn("anthocyanin", rationale)
        self.assertNotIn("flavonoid", rationale)

    def test_pea_a2_ttg1_is_headline_and_closest_wd40_is_disclosed(self) -> None:
        """#422/#424: pea_a2 (TTG1) keeps the informative TTG1 annotation as the
        headline, and the closer but uninformative WD40-1 hit is disclosed."""
        hits = [
            {"subject": "tr|A0A0C3VMP6|A0A0C3VMP6_MEDTR", "identity": 94.48, "query_coverage": 1.003, "subject_coverage": 1.006, "evalue": 0.0, "bitscore": 668.0, "description": "WD40-1 protein OS=Medicago truncatula PE=2"},
            {"subject": "tr|B4XH97|B4XH97_MEDTR", "identity": 94.48, "query_coverage": 1.003, "subject_coverage": 1.006, "evalue": 0.0, "bitscore": 668.0, "description": "Transparent testa glabra protein OS=Medicago truncatula GN=MTR_3g092840"},
            {"subject": "LotjaGi1g1v0299500.1", "identity": 89.14, "query_coverage": 1.02, "subject_coverage": 1.003, "evalue": 0.0, "bitscore": 632.0, "description": "LotjaGi1g1v0299500.1"},
            {"subject": "tr|C0A1G9|C0A1G9_LOTJA", "identity": 89.14, "query_coverage": 1.02, "subject_coverage": 1.003, "evalue": 0.0, "bitscore": 632.0, "description": "WD40 repeats protein LjTTG1 OS=Lotus japonicus GN=LjTTG1"},
            {"subject": "sp|Q9XGN1|TTG1_ARATH", "identity": 79.19, "query_coverage": 1.009, "subject_coverage": 1.015, "evalue": 0.0, "bitscore": 560.0, "description": "Protein TRANSPARENT TESTA GLABRA 1 OS=Arabidopsis thaliana OX=3702 GN=TTG1"},
        ]
        state, _ = homology_qc.assess_annotation_concordance(hits)
        self.assertEqual(state, "single_informative")
        rationale = homology_qc.format_grounded_rationale(hits=hits, include_disclaimer=False)
        self.assertIn("Q9XGN1", rationale)
        self.assertIn("TRANSPARENT TESTA GLABRA 1", rationale)
        self.assertIn("WD40-1 protein", rationale)
        self.assertIn("A0A0C3VMP6", rationale)
        self.assertIn(homology_qc.HEADLINE_LABEL_CLOSEST, rationale)
        self.assertIn(homology_qc.HEADLINE_LABEL_INFORMATIVE, rationale)


class Dimension6InvarianceOfWeightAndRankingFilesTests(unittest.TestCase):
    """Dimension 6: 0 diff lines in weight/ranking files, homology_context direction is context with weight 0."""

    def test_registry_homology_context_is_weight_zero_and_context(self) -> None:
        frame, _ = load_weight_registry()
        match = frame[frame["evidence_code"] == "homology_context"]
        self.assertFalse(match.empty)
        row = match.iloc[0]
        self.assertEqual(row["base_weight_or_tier"], "")
        self.assertEqual(row["default_direction"], "context")

    def test_git_diff_invariance_against_origin_main(self) -> None:
        # Fails loudly when git is unavailable rather than silently passing.
        # Note: evidence_policy.py is excluded here because issue #421 canonically
        # normalizes paralog_ambiguous to direction=context at the evidence-build
        # boundary. The guard's goal is preventing benchmark-fitted weight/threshold
        # modifications, not blocking evidence direction semantic fixes.
        result = subprocess.run(
            ["git", "diff", "origin/main", "--",
             "src/candidate_qtl/resources/evidence_weight_registry.tsv",
             "src/candidate_qtl/ranking_policy.py"],
            capture_output=True, text=True
        )
        self.assertEqual(result.returncode, 0, f"git diff failed: {result.stderr}")
        self.assertEqual(result.stdout.strip(), "", "Invariants violated: diff found against origin/main")


if __name__ == "__main__":
    unittest.main()
