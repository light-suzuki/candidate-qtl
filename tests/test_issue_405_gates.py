"""Issue #405: confirm the existing Route A / Route B0 gates still hold.

Route A's first hypothesis generation must stay interval-blind and Route B0's
interval inventory must stay phenotype-free. The gates themselves are not
implemented here (their implementation is owned by the Route A/B lanes); this
file only confirms that the current implementation still enforces them.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import build_parser, main


BARE_CANDIDATES = "\n".join([
    "gene_id\tchrom\tstart\tend\tdescription",
    "Psat_G1\tchr5\t1000\t3000\t",
    "Psat_G2\tchr5\t9000\t11000\t",
]) + "\n"

ORTHOLOGY = "\n".join([
    "gene_id\treference_gene_id\torthogroup\trelationship\tprovider",
    "Psat_G1\tAtCLV1\tOG0002\tbest_hit_homology\tlocal_blastp:arabidopsis_uniprot",
    "Psat_G2\tAtSBE1\tOG0001\tbest_hit_homology\tlocal_blastp:arabidopsis_uniprot",
]) + "\n"

REFERENCE_CORPUS = "\n".join([
    "reference_gene_id\treference_species\tfamily_id\torthogroup_id\tdescription\tmutant_phenotype\tpathway\tsource_id",
    "AtCLV1\tArabidopsis thaliana\tCLV\tOG0002\tCLAVATA1 receptor kinase\tfasciated stem\tmeristem_maintenance\tPMID:8402901",
    "AtSBE1\tArabidopsis thaliana\tSBE\tOG0001\tstarch branching enzyme I\twrinkled seed\tstarch_biosynthesis\tPMID:1978859",
]) + "\n"


def _write(root: Path, name: str, content: str) -> Path:
    path = root / name
    path.write_text(content, encoding="utf-8")
    return path


def _parser_options(command: str) -> set[str]:
    parser = build_parser()
    action = [item for item in parser._subparsers._group_actions if hasattr(item, "choices")][0]
    options: set[str] = set()
    for sub_action in action.choices[command]._actions:
        options.update(sub_action.option_strings)
    return options


class RouteAIntervalBlindGateTests(unittest.TestCase):
    """Gate 1: Route A first hypothesis generation is interval-blind (#405)."""

    def _investigate(self, root: Path, phenotype: str = "wrinkled seed with reduced starch content",
                     candidates: Path | None = None) -> Path:
        out = root / "investigation"
        argv = [
            "phenotype-investigate", "--out-dir", str(out), "--case-id", "case_r",
            "--focal-species", "Pisum sativum", "--phenotype", phenotype,
            "--direction", "reduced", "--inheritance", "recessive",
            "--tissue", "seed", "--stage", "seed filling", "--onset", "late",
            "--penetrance", "complete", "--synonym", "wrinkled seed",
            "--ontology-term", "TO:0000123",
            "--uncertainty", "starch content was scored on a single harvest",
        ]
        if candidates is not None:
            argv += ["--candidates", str(candidates)]
        self.assertEqual(0, main(argv))
        return out

    def test_first_hypothesis_command_has_no_interval_input(self) -> None:
        for forbidden in ("--interval", "--atlas", "--gff3", "--interval-genes"):
            self.assertNotIn(forbidden, _parser_options("phenotype-investigate"))

    def test_investigation_manifest_records_the_interval_blind_gate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = self._investigate(root)
            manifest = json.loads((out / "phenotype_investigation_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual([], manifest["interval_blind"]["interval_arguments"])
            self.assertFalse(manifest["interval_blind"]["qtl_candidate_information_in_packet"])
            self.assertIn("interval_blind", manifest)

    def test_phenotype_naming_an_interval_candidate_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = _write(root, "candidates.tsv", "gene_id\nPsat_G0001\n")
            with self.assertRaises(ValueError) as context:
                self._investigate(root, phenotype="wrinkled seed caused by Psat_G0001", candidates=candidates)
            self.assertIn("Psat_G0001", str(context.exception))


class RouteB0PhenotypeFreeGateTests(unittest.TestCase):
    """Gate 2: Route B0 interval inventory is phenotype-free (#405)."""

    def test_interval_atlas_command_has_no_phenotype_option(self) -> None:
        self.assertNotIn("--phenotype", _parser_options("interval-atlas"))

    def test_atlas_manifest_records_phenotype_independence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "route_b"
            self.assertEqual(0, main([
                "route-b-prepare",
                "--candidates", str(_write(root, "candidates.tsv", BARE_CANDIDATES)),
                "--orthology", str(_write(root, "orthology.tsv", ORTHOLOGY)),
                "--reference-function", str(_write(root, "corpus.tsv", REFERENCE_CORPUS)),
                "--out-dir", str(out),
            ]))
            manifest = json.loads((out / "route_b_preparation_manifest.json").read_text(encoding="utf-8"))
            self.assertFalse(manifest["phenotype_used"])
            self.assertTrue(manifest["atlas_phenotype_independent"])
            atlas = json.loads(Path(manifest["outputs"]["interval_atlas_manifest"]).read_text(encoding="utf-8"))
            self.assertFalse(atlas["phenotype_used"])

    def test_atlas_stays_phenotype_free_even_when_b_later_sees_the_phenotype(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "route_b"
            self.assertEqual(0, main([
                "route-b-prepare",
                "--candidates", str(_write(root, "candidates.tsv", BARE_CANDIDATES)),
                "--orthology", str(_write(root, "orthology.tsv", ORTHOLOGY)),
                "--reference-function", str(_write(root, "corpus.tsv", REFERENCE_CORPUS)),
                "--phenotype", "fasciated stem with enlarged meristem",
                "--out-dir", str(out),
            ]))
            manifest = json.loads((out / "route_b_preparation_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(manifest["phenotype_used"])
            self.assertTrue(manifest["atlas_phenotype_independent"])
            atlas = json.loads(Path(manifest["outputs"]["interval_atlas_manifest"]).read_text(encoding="utf-8"))
            self.assertFalse(atlas["phenotype_used"])

    def test_family_review_queries_are_phenotype_free(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._atlas_fixture(root)
            out = root / "atlas"
            self.assertEqual(0, main([
                "interval-atlas",
                "--candidates", str(fixture["candidates.tsv"]),
                "--family-index", str(fixture["family_index.tsv"]),
                "--domains", str(fixture["domains.tsv"]),
                "--orthology", str(fixture["orthology.tsv"]),
                "--annotations", str(fixture["annotations.tsv"]),
                "--flags", str(fixture["flags.tsv"]),
                "--out-dir", str(out), "--case-id", "demo",
            ]))
            manifest = json.loads((out / "interval_atlas_manifest.json").read_text(encoding="utf-8"))
            self.assertFalse(manifest["phenotype_used"])
            queries = pd.read_csv(manifest["outputs"]["family_review_queries"], sep="\t", dtype=str).fillna("")
            self.assertTrue(len(queries))
            self.assertTrue(queries["phenotype_used"].eq("False").all())

    @staticmethod
    def _atlas_fixture(root: Path) -> dict[str, Path]:
        paths = {name: root / name for name in (
            "candidates.tsv", "family_index.tsv", "domains.tsv", "orthology.tsv", "annotations.tsv", "flags.tsv")}
        paths["candidates.tsv"].write_text("\n".join([
            "gene_id\tchrom\tstart\tend\tprotein_length",
            "G_RLK1\tchr1\t1000\t3000\t980",
            "G_RLK2\tchr1\t5000\t7000\t960",
            "G_TF1\tchr1\t20000\t22000\t420",
            "G_ORPHAN\tchr1\t40000\t40300\t85",
        ]) + "\n", encoding="utf-8")
        paths["family_index.tsv"].write_text("\n".join([
            "gene_id\tchrom\tstart\tend\tbest_subject\tbest_subject_symbol\tbest_subject_stem\tfamily_stems\tassignment_status",
            "G_RLK1\tchr1\t1000\t3000\tsp|P1|CLV1_ARATH\tCLV1\tCLV\tCLV\tassigned",
            "G_RLK2\tchr1\t5000\t7000\tsp|P2|CLV2_ARATH\tCLV2\tCLV\tCLV\tassigned",
            "G_TF1\tchr1\t20000\t22000\tsp|P4|WOX4_ARATH\tWOX4\tWOX\tWOX\tassigned",
            "G_ORPHAN\t\t\t\t\t\t\t\tdark_no_qualifying_hit",
        ]) + "\n", encoding="utf-8")
        paths["domains.tsv"].write_text("\n".join([
            "gene_id\tdomain_id\tdomain_name\tstart\tsource",
            "G_RLK1\tPF00560\tLRR\t50\tPfam",
            "G_RLK1\tPF00069\tPkinase\t600\tPfam",
            "G_RLK2\tPF00560\tLRR\t60\tPfam",
            "G_RLK2\tPF00069\tPkinase\t610\tPfam",
        ]) + "\n", encoding="utf-8")
        paths["orthology.tsv"].write_text("\n".join([
            "gene_id\torthogroup\trelationship\treference_gene_id\tprovider",
            "G_RLK1\tOG0000123\tone-to-one\tCLV1\torthofinder",
            "G_RLK2\tOG0000123\tone-to-many\tCLV2\torthofinder",
            "G_TF1\tOG0000456\tone-to-one\tWOX4\torthofinder",
        ]) + "\n", encoding="utf-8")
        paths["annotations.tsv"].write_text("\n".join([
            "gene_id\tdescription\tgo_terms",
            "G_RLK1\tleucine-rich repeat receptor kinase, meristem maintenance\tGO:0010073;GO:0004672",
            "G_RLK2\treceptor-like kinase\t",
            "G_TF1\thomeobox transcription factor\tGO:0006355",
        ]) + "\n", encoding="utf-8")
        paths["flags.tsv"].write_text("gene_id\tdark_state\nG_ORPHAN\tdark_unresolved\n", encoding="utf-8")
        return paths


if __name__ == "__main__":
    unittest.main()
