"""Issues #478/#525: accession-native mutant evidence reaches Route B reference functions.

The chain ``interval gene -> reference gene -> that reference gene's reported
mutant phenotype`` must stay traceable at accession level.  The bridge projects
``normalized_functional_evidence.tsv`` rows (LORE1/Tnt1 adapters) onto the
existing ``REFERENCE_FUNCTION_COLUMNS`` contract using an explicit
cross-accession map:

- only verified one_to_one mappings may replace the native ID;
- one_to_many keeps the native ID as ``ambiguous`` with visible targets;
- unmapped rows stay ``not_evaluated``, never absence;
- symbol-only bases are rejected fail-closed;
- native IDs are matched case-preservingly.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.mutant_reference_bridge import build_mutant_reference_functions

EVIDENCE_HEADER = (
    "native_gene_id\tnative_species\tnative_accession\tnative_assembly\t"
    "native_annotation\tnative_gene_id_namespace\tgene_symbol\t"
    "phenotype_or_function\tontology_ids\tsource_ids\tsource_confidence\t"
    "evidence_type\tmutant_or_perturbation_background\n"
)


def _evidence_row(gene_id: str, symbol: str, phenotype: str,
                  evidence_class: str, source_ids: str = "PMID:14534591") -> str:
    return "\t".join([
        gene_id, "Lotus japonicus", "Gifu B-129", "Gifu_v1.2", "Gifu_v1.3",
        "LotusBase:Gifu_v1.2", symbol, phenotype, "GO:0009877", source_ids,
        evidence_class, "curated_lore1", "Gifu B-129",
    ]) + "\n"


EVIDENCE = (
    EVIDENCE_HEADER
    + _evidence_row("LotjaGi1g1v0004000", "NODK",
                    "nodulation deficient; infection threads blocked in the epidermis",
                    "strong_causal")
    + _evidence_row("LotjaGi1g1v0005000", "NODL",
                    "delayed nodulation with pale pink nodules",
                    "published_phenotype_only")
    + _evidence_row("LotjaGi1g1v0006000", "INSX", "", "insertion_observed_only")
)

CROSSWALK = "\n".join([
    "source_full_id\ttarget_full_id\tmapping_basis\tmapping_record\tfam3_family_id\t"
    "target_species\trelationship\tmapping_state\tmapping_reason",
    "LotjaGi1g1v0004000\tlotja.MG20.gnm3.ann1.Lj2G3V3373110\texplicit_id_map\t"
    "LotusBase:TranscriptMapping:mg20_3.0-gifu_1.2:Lj2G3V3373110\tLYM\t"
    "Lotus japonicus MG20\tone_to_one\tverified\treciprocal best hit confirmed",
    "LotjaGi1g1v0005000\tlotja.MG20.gnm3.ann1.Lj2G3V1111001\tpangene_family\tfam3:LYM\tLYM\t"
    "Lotus japonicus MG20\tone_to_many\tverified\ttwo MG20 paralogs in the same fam3 family",
    "LotjaGi1g1v0005000\tlotja.MG20.gnm3.ann1.Lj2G3V2222002\tpangene_family\tfam3:LYM\tLYM\t"
    "Lotus japonicus MG20\tone_to_many\tverified\ttwo MG20 paralogs in the same fam3 family",
]) + "\n"

CANDIDATES = "\n".join([
    "gene_id\tchrom\tstart\tend\tdescription",
    "Psat_G1\tchr5\t1000\t3000\tlysM receptor kinase like protein",
    "Psat_G2\tchr5\t9000\t11000\tstorage protein",
]) + "\n"

ORTHOLOGY = "\n".join([
    "gene_id\treference_gene_id\torthogroup\trelationship\tprovider",
    "Psat_G1\tQ70KR8\tOG_LYM\tbest_hit_homology\tlocal_blastp:legume_ref",
    "Psat_G2\tP13915\tOG_CVC\tbest_hit_homology\tlocal_blastp:legume_ref",
]) + "\n"

CORPUS = "\n".join([
    "reference_gene_id\treference_species\tfamily_id\torthogroup_id\tdescription\t"
    "mutant_phenotype\tpathway\tsource_id",
    "Q70KR8\tLotus japonicus\tLYM\tOG_LYM\tLysM receptor kinase\t"
    "nodulation deficient\troot_nodulation\tPMID:14534591",
    "P13915\tPisum sativum\tCVC\tOG_CVC\tconvicilin\t\tseed_storage\tPMID:3415641",
]) + "\n"


def _write(root: Path, name: str, content: str) -> Path:
    path = root / name
    path.write_text(content, encoding="utf-8")
    return path


class BridgeProjectionTests(unittest.TestCase):
    def test_verified_one_to_one_promotes_target_and_records_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = build_mutant_reference_functions(
                [_write(root, "evidence.tsv", EVIDENCE)],
                [_write(root, "map.tsv", CROSSWALK)],
                root,
            )
            frame = pd.read_csv(root / "mutant_reference_functions.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(2, len(frame))
            promoted = frame[frame["gene_symbol"] == "NODK"].iloc[0]
            self.assertEqual("lotja.MG20.gnm3.ann1.Lj2G3V3373110", promoted["reference_gene_id"])
            self.assertEqual("verified", promoted["mapping_state"])
            self.assertEqual(
                "nodulation deficient; infection threads blocked in the epidermis",
                promoted["mutant_phenotype"])
            provenance = json.loads(promoted["provenance_detail"])
            self.assertEqual("Gifu B-129", provenance["native"]["accession"])
            self.assertEqual("explicit_id_map", provenance["mapping"]["basis"])
            self.assertEqual("strong_causal", provenance["evidence_class"])
            # assembly and annotation generations stay separate dimensions
            self.assertEqual("Gifu_v1.2", provenance["native"]["assembly"])
            self.assertEqual("Gifu_v1.3", provenance["native"]["annotation"])
            self.assertEqual("none", manifest["rank_contribution"])
            self.assertEqual(1, manifest["resolved_one_to_one_rows"])

    def test_one_to_many_keeps_the_native_id_as_ambiguous(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            build_mutant_reference_functions(
                [_write(root, "evidence.tsv", EVIDENCE)],
                [_write(root, "map.tsv", CROSSWALK)],
                root,
            )
            frame = pd.read_csv(root / "mutant_reference_functions.tsv", sep="\t", dtype=str).fillna("")
            row = frame[frame["gene_symbol"] == "NODL"].iloc[0]
            self.assertEqual("LotjaGi1g1v0005000", row["reference_gene_id"])
            self.assertEqual("ambiguous", row["mapping_state"])
            self.assertIn("one_to_many", row["mapping_relationship"])
            provenance = json.loads(row["provenance_detail"])
            self.assertEqual(2, len(provenance["mapping"]["targets"]))

    def test_unmapped_rows_stay_native_and_not_evaluated(self) -> None:
        header_only_map = (
            "source_full_id\ttarget_full_id\tmapping_basis\tmapping_record\tfam3_family_id\t"
            "target_species\trelationship\tmapping_state\tmapping_reason\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = build_mutant_reference_functions(
                [_write(root, "evidence.tsv", EVIDENCE)],
                [_write(root, "map.tsv", header_only_map)],
                root,
            )
            frame = pd.read_csv(root / "mutant_reference_functions.tsv", sep="\t", dtype=str).fillna("")
            row = frame[frame["gene_symbol"] == "NODK"].iloc[0]
            self.assertEqual("LotjaGi1g1v0004000", row["reference_gene_id"])
            self.assertEqual("not_evaluated", row["mapping_state"])
            # the phenotype text stays visible even while identity is unresolved
            self.assertTrue(row["mutant_phenotype"])
            self.assertEqual(2, manifest["unmapped_rows"])

    def test_insertion_only_rows_are_skipped_and_never_negative(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = build_mutant_reference_functions(
                [_write(root, "evidence.tsv", EVIDENCE)], [], root,
            )
            frame = pd.read_csv(root / "mutant_reference_functions.tsv", sep="\t", dtype=str).fillna("")
            self.assertNotIn("INSX", set(frame["gene_symbol"]))
            self.assertEqual(1, manifest["skipped_not_mutant_evidence"])

    def test_symbol_only_mapping_basis_is_rejected(self) -> None:
        bad_map = "\n".join([
            "source_full_id\ttarget_full_id\tmapping_basis\trelationship\tmapping_state",
            "LotjaGi1g1v0004000\tlotja.MG20.gnm3.ann1.XXX\tsymbol_match\tone_to_one\tverified",
        ]) + "\n"
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaisesRegex(ValueError, "symbol-only"):
                build_mutant_reference_functions(
                    [_write(root, "evidence.tsv", EVIDENCE)],
                    [_write(root, "bad_map.tsv", bad_map)],
                    root,
                )

    def test_native_id_matching_is_case_preserving(self) -> None:
        lowered = CROSSWALK.replace("LotjaGi1g1v0004000", "lotjagi1g1v0004000")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            build_mutant_reference_functions(
                [_write(root, "evidence.tsv", EVIDENCE)],
                [_write(root, "map.tsv", lowered)],
                root,
            )
            frame = pd.read_csv(root / "mutant_reference_functions.tsv", sep="\t", dtype=str).fillna("")
            row = frame[frame["gene_symbol"] == "NODK"].iloc[0]
            self.assertEqual("LotjaGi1g1v0004000", row["reference_gene_id"])
            self.assertEqual("not_evaluated", row["mapping_state"])


class RouteBPrepareWiringTests(unittest.TestCase):
    def test_flags_project_mutant_evidence_into_route_b_reference_functions(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "route_b"
            self.assertEqual(0, main([
                "route-b-prepare",
                "--candidates", str(_write(root, "candidates.tsv", CANDIDATES)),
                "--orthology", str(_write(root, "orthology.tsv", ORTHOLOGY)),
                "--reference-function", str(_write(root, "corpus.tsv", CORPUS)),
                "--mutant-functional-evidence", str(_write(root, "evidence.tsv", EVIDENCE)),
                "--cross-accession-map", str(_write(root, "map.tsv", CROSSWALK)),
                "--out-dir", str(out),
                "--case-id", "issue478",
            ]))
            bridge_manifest_path = out / "mutant_reference_bridge_manifest.json"
            self.assertTrue(bridge_manifest_path.exists())
            manifest = json.loads(bridge_manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(2, manifest["emitted_rows"])
            self.assertEqual(2, manifest["rows_with_mutant_phenotype"])
            table = pd.read_csv(out / "route_b_reference_function.tsv", sep="\t", dtype=str).fillna("")
            merged = table[table["source_database"] == "curated_mutant_evidence_bridge"]
            self.assertEqual(2, len(merged))
            self.assertIn("lotja.MG20.gnm3.ann1.Lj2G3V3373110", set(merged["reference_gene_id"]))
            self.assertTrue((merged["mutant_phenotype"] != "").all())
            preparation = json.loads(
                (out / "route_b_preparation_manifest.json").read_text(encoding="utf-8"))
            self.assertIn("mutant_reference_functions.tsv",
                          json.dumps(preparation["providers"]["reference_function"]["inputs"]))


if __name__ == "__main__":
    unittest.main()
