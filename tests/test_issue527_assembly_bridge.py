from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.assembly_bridge import build_assembly_bridge
from candidate_qtl.cli import main


FUNCTIONS = """source_record_id\tnative_gene_id\tnative_species\tnative_species_code\tnative_accession\tnative_assembly\tnative_annotation\tnative_gene_id_namespace\tphenotype_or_function\tsource_ids
card-1\tpissa.Cameor.gnm1.ann1.S1\tPisum sativum\tpissa\tCameor\tgnm1\tann1\tpissa.Cameor.gnm1.ann1\tseed development function\tPMID:1
card-2\tpissa.Cameor.gnm1.ann1.S2\tPisum sativum\tpissa\tCameor\tgnm1\tann1\tpissa.Cameor.gnm1.ann1\troot development function\tPMID:2
card-3\tpissa.Cameor.gnm1.ann1.S3\tPisum sativum\tpissa\tCameor\tgnm1\tann1\tpissa.Cameor.gnm1.ann1\tflower function\tPMID:3
card-4\tpissa.Cameor.gnm1.ann1.S4\tPisum sativum\tpissa\tCameor\tgnm1\tann1\tpissa.Cameor.gnm1.ann1\tunknown function\tPMID:4
card-5\tpissa.Cameor.gnm1.ann1.S5\tPisum sativum\tpissa\tCameor\tgnm1\tann1\tpissa.Cameor.gnm1.ann1\tleaf function\tPMID:5
"""
SOURCE_FAM3 = """native_gene_id\tfamily_id\tnative_gene_id_namespace
pissa.Cameor.gnm1.ann1.S1\tlegume.fam3.F1\tpissa.Cameor.gnm1.ann1
pissa.Cameor.gnm1.ann1.S2\tlegume.fam3.F2\tpissa.Cameor.gnm1.ann1
pissa.Cameor.gnm1.ann1.S3\tlegume.fam3.F3\tpissa.Cameor.gnm1.ann1
pissa.Cameor.gnm1.ann1.S4\tlegume.fam3.NO_TARGET\tpissa.Cameor.gnm1.ann1
pissa.Cameor.gnm1.ann1.S5\tlegume.fam3.F5\tpissa.Cameor.gnm1.ann1
"""
TARGET_FAM3 = """native_gene_id\tfamily_id\tnative_gene_id_namespace
pissa.ZW6.gnm1.ann1.T1\tlegume.fam3.F1\tpissa.ZW6.gnm1.ann1
pissa.ZW6.gnm1.ann1.T2\tlegume.fam3.F2\tpissa.ZW6.gnm1.ann1
pissa.ZW6.gnm1.ann1.T3\tlegume.fam3.F2\tpissa.ZW6.gnm1.ann1
pissa.ZW6.gnm1.ann1.T4\tlegume.fam3.OTHER\tpissa.ZW6.gnm1.ann1
pissa.ZW6.gnm1.ann1.T5\tlegume.fam3.F5\tpissa.ZW6.gnm1.ann1
"""
SEQUENCES = """source_full_id\ttarget_full_id\treciprocal_sequence_state\tsequence_identity\tquery_coverage\tsubject_coverage\tsynteny_state\tdomain_concordance
pissa.Cameor.gnm1.ann1.S1\tpissa.ZW6.gnm1.ann1.T1\tunique_reciprocal\t98.1\t1.0\t1.0\tconcordant\tconcordant
pissa.Cameor.gnm1.ann1.S2\tpissa.ZW6.gnm1.ann1.T2\tunique_reciprocal\t96.0\t0.99\t1.0\tconcordant\tconcordant
pissa.Cameor.gnm1.ann1.S2\tpissa.ZW6.gnm1.ann1.T3\tunique_reciprocal\t95.8\t0.98\t1.0\tconcordant\tconcordant
pissa.Cameor.gnm1.ann1.S3\tpissa.ZW6.gnm1.ann1.T4\tunique_reciprocal\t97.2\t1.0\t1.0\tconcordant\tconflict
pissa.Cameor.gnm1.ann1.S5\tpissa.ZW6.gnm1.ann1.T5\tunique_reciprocal\t97.0\t1.0\t1.0\tconflict\tconcordant
"""


class AssemblyBridgeTests(unittest.TestCase):
    def _fixture(self, root: Path) -> dict[str, Path]:
        paths = {}
        for name, text in (("functions.tsv", FUNCTIONS), ("source.tsv", SOURCE_FAM3),
                           ("target.tsv", TARGET_FAM3), ("sequence.tsv", SEQUENCES)):
            path = root / name
            path.write_text(text, encoding="utf-8")
            paths[name] = path
        return paths

    def test_mapping_preserves_native_ids_and_conservative_states(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            paths = self._fixture(root)
            result = build_assembly_bridge(paths["functions.tsv"], paths["source.tsv"], paths["target.tsv"], paths["sequence.tsv"], root / "bridge")
            mapping = pd.read_csv(result["mapping"], sep="\t", dtype=str).fillna("")
            s1 = mapping[mapping.source_full_id.eq("pissa.Cameor.gnm1.ann1.S1")].iloc[0]
            self.assertEqual("pissa.ZW6.gnm1.ann1.T1", s1.target_full_id)
            self.assertEqual("verified", s1.mapping_state)
            self.assertEqual("one_to_one", s1.relationship)
            s2 = mapping[mapping.source_full_id.eq("pissa.Cameor.gnm1.ann1.S2")]
            self.assertEqual({"one_to_many"}, set(s2.relationship))
            self.assertEqual({"ambiguous"}, set(s2.mapping_state))
            s3 = mapping[mapping.source_full_id.eq("pissa.Cameor.gnm1.ann1.S3")].iloc[0]
            self.assertEqual("provisional", s3.mapping_state)
            self.assertIn("fam3_mismatch", s3.mapping_basis)
            s4 = mapping[mapping.source_full_id.eq("pissa.Cameor.gnm1.ann1.S4")].iloc[0]
            self.assertEqual("not_evaluated", s4.mapping_state)
            self.assertIn("absence_is_not_negative", s4.mapping_reason)
            s5 = mapping[mapping.source_full_id.eq("pissa.Cameor.gnm1.ann1.S5")].iloc[0]
            self.assertEqual("ambiguous", s5.mapping_state)
            self.assertIn("synteny_conflict", s5.mapping_reason)

            orthology = pd.read_csv(result["route_b_orthology"], sep="\t", dtype=str).fillna("")
            self.assertIn("pissa.Cameor.gnm1.ann1.S1", set(orthology.reference_gene_id))
            self.assertIn("pissa.ZW6.gnm1.ann1.T1", set(orthology.gene_id))
            self.assertIn("pissa.Cameor.gnm1.ann1.S2", set(orthology.reference_gene_id))
            self.assertNotIn("pissa.Cameor.gnm1.ann1.S3", set(orthology.reference_gene_id))
            reference = pd.read_csv(result["route_b_reference_function"], sep="\t", dtype=str).fillna("")
            self.assertEqual("pissa.Cameor.gnm1.ann1.S1", reference.iloc[0].reference_gene_id)
            self.assertNotIn("pissa.Cameor.gnm1.ann1.S3", set(reference.reference_gene_id))
            one_to_many = reference[reference.reference_gene_id.eq("pissa.Cameor.gnm1.ann1.S2")]
            self.assertTrue(one_to_many.gene_id.eq("").all())
            self.assertEqual({"family_set"}, set(one_to_many.function_transfer_scope))
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual("none", manifest["rank_contribution"])
            self.assertEqual(1, manifest["counts"]["mapping_verified_one_to_one"])

    def test_cli_and_actual_route_b_prepare_reach_source_backed_reverse_lookup(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            paths = self._fixture(root)
            bridge_dir = root / "bridge"
            self.assertEqual(0, main([
                "assembly-bridge", "--functional-evidence", str(paths["functions.tsv"]),
                "--source-fam3", str(paths["source.tsv"]), "--target-fam3", str(paths["target.tsv"]),
                "--sequence-evidence", str(paths["sequence.tsv"]), "--out-dir", str(bridge_dir),
            ]))
            candidates = root / "candidates.tsv"
            pd.DataFrame({"gene_id": [f"pissa.ZW6.gnm1.ann1.T{i}" for i in range(1, 6)],
                          "chrom": ["1"] * 5, "start": list(range(1, 6)), "end": list(range(2, 7)),
                          "protein_length": ["300"] * 5}).to_csv(candidates, sep="\t", index=False)
            route = root / "route-b"
            result = main([
                "route-b-prepare", "--candidates", str(candidates), "--out-dir", str(route),
                "--orthology", str(bridge_dir / "route_b_orthology.tsv"),
                "--reference-function", str(bridge_dir / "route_b_reference_function.tsv"),
                "--no-auto-local-db", "--json",
            ])
            self.assertEqual(0, result)
            manifest = json.loads((route / "route_b_preparation_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("source_backed_reverse_lookup", manifest["route_b_readiness"])
            self.assertEqual("source_backed_reverse_lookup", manifest["reverse_lookup_status"])
            status = json.loads((route / "route_b0" / "interval_first_status.json").read_text(encoding="utf-8"))
            self.assertEqual("source_backed_reverse_lookup", status["reverse_lookup_status"])

    def test_multi_reference_union_preserves_alternate_only_and_copy_sets(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            paths = self._fixture(root)
            bridge_dir = root / "bridge"
            result = build_assembly_bridge(paths["functions.tsv"], paths["source.tsv"], paths["target.tsv"], paths["sequence.tsv"], bridge_dir)
            primary = root / "primary.tsv"
            alternate = root / "alternate.tsv"
            primary.write_text(
                "gene_id\tassembly\tchrom\tstart\tend\tstatus\n"
                "pissa.Cameor.gnm1.ann1.S1\tCameor\t1\t1\t2\tprimary\n"
                "pissa.Cameor.gnm1.ann1.S2\tCameor\t1\t3\t4\tprimary\n"
                "pissa.Cameor.gnm1.ann1.S3\tCameor\t1\t5\t6\tprimary\n"
                "pissa.Cameor.gnm1.ann1.S4\tCameor\t1\t7\t8\tannotation_rescue\n", encoding="utf-8")
            alternate.write_text(
                "gene_id\tassembly\tchrom\tstart\tend\tstatus\n"
                "pissa.ZW6.gnm1.ann1.T1\tZW6\t1\t10\t11\tnative\n"
                "pissa.ZW6.gnm1.ann1.T2\tZW6\t1\t12\t13\tnative\n"
                "pissa.ZW6.gnm1.ann1.T3\tZW6\t1\t14\t15\tnative\n"
                "pissa.ZW6.gnm1.ann1.ZONLY\tZW6\t2\t20\t21\tPAV\n", encoding="utf-8")
            union = root / "union"
            cli_result = main([
                "multi-reference-locus-union", "--primary", str(primary),
                "--alternate", f"ZW6={alternate}", "--mapping", result["mapping"],
                "--out-dir", str(union), "--json",
            ])
            self.assertEqual(0, cli_result)
            candidates = pd.read_csv(union / "candidate_universe.tsv", sep="\t", dtype=str).fillna("")
            # Source-only (S4) and alternate-only (ZONLY) are retained.
            self.assertIn("pissa.Cameor.gnm1.ann1.S4", set(candidates.native_gene_id))
            self.assertIn("pissa.ZW6.gnm1.ann1.ZONLY", set(candidates.native_gene_id))
            self.assertEqual("primary_native", candidates.loc[candidates.native_gene_id.eq("pissa.Cameor.gnm1.ann1.S1"), "candidate_status"].iloc[0])
            self.assertEqual("verified_equivalent", candidates.loc[candidates.native_gene_id.eq("pissa.ZW6.gnm1.ann1.T1"), "candidate_status"].iloc[0])
            # A verified 1:1 duplicate shares one entity; the ambiguous F2 set
            # retains S2, T2 and T3 as three separate entities/copies.
            s1_entity = candidates.loc[candidates.native_gene_id.eq("pissa.Cameor.gnm1.ann1.S1"), "locus_entity_id"].iloc[0]
            t1_entity = candidates.loc[candidates.native_gene_id.eq("pissa.ZW6.gnm1.ann1.T1"), "locus_entity_id"].iloc[0]
            self.assertEqual(s1_entity, t1_entity)
            f2 = candidates[candidates.native_gene_id.isin([
                "pissa.Cameor.gnm1.ann1.S2", "pissa.ZW6.gnm1.ann1.T2", "pissa.ZW6.gnm1.ann1.T3"])]
            self.assertEqual({"ambiguous_copy_set"}, set(f2.candidate_status))
            self.assertEqual(3, f2.locus_entity_id.nunique())
            manifest = json.loads((union / "multi_reference_locus_union_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(len(candidates), manifest["counts"]["recall_universe_count"])
            self.assertLess(manifest["counts"]["conditional_compression_input_count"], manifest["counts"]["candidate_universe_count"])
            compression = pd.read_csv(union / "compression_candidates.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("gene_id", compression.columns)
            # Verified S1/T1 representations collapse to the primary S1 row.
            self.assertIn("pissa.Cameor.gnm1.ann1.S1", set(compression.gene_id))
            self.assertNotIn("pissa.ZW6.gnm1.ann1.T1", set(compression.gene_id))
            # Ambiguous one-to-many copies and an alternate-only/PAV gene all remain.
            self.assertTrue({
                "pissa.Cameor.gnm1.ann1.S2", "pissa.ZW6.gnm1.ann1.T2",
                "pissa.ZW6.gnm1.ann1.T3", "pissa.ZW6.gnm1.ann1.ZONLY",
            }.issubset(set(compression.gene_id)))
            ledger = pd.read_csv(union / "locus_link_ledger.tsv", sep="\t", dtype=str).fillna("")
            ambiguous = ledger[ledger.relationship.eq("one_to_many")]
            self.assertEqual({"false"}, set(ambiguous.function_transfer_allowed))
            self.assertEqual({"family_set"}, set(ambiguous.function_transfer_scope))
            route = root / "route-b-union"
            self.assertEqual(0, main([
                "route-b-prepare", "--candidates", str(union / "compression_candidates.tsv"),
                "--out-dir", str(route), "--no-auto-local-db", "--json",
            ]))
            atlas = pd.read_csv(route / "interval_gene_atlas.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(set(compression.gene_id), set(atlas.gene_id))
            self.assertIn("pissa.ZW6.gnm1.ann1.ZONLY", set(atlas.gene_id))

    def test_v2_only_gene_survives_when_zw6_is_primary(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            paths = self._fixture(root)
            bridge = build_assembly_bridge(
                paths["functions.tsv"], paths["source.tsv"], paths["target.tsv"],
                paths["sequence.tsv"], root / "bridge",
            )
            zw6 = root / "zw6.tsv"
            v2 = root / "v2.tsv"
            zw6.write_text(
                "gene_id\tassembly\tchrom\tstart\tend\n"
                "pissa.ZW6.gnm1.ann1.T1\tZW6\t1\t10\t11\n",
                encoding="utf-8",
            )
            v2.write_text(
                "gene_id\tassembly\tchrom\tstart\tend\tstatus\n"
                "pissa.Cameor.gnm1.ann1.S1\tCameorV2\t1\t1\t2\tnative\n"
                "pissa.Cameor.gnm1.ann1.S4\tCameorV2\t1\t7\t8\tPAV\n",
                encoding="utf-8",
            )
            out = root / "union"
            self.assertEqual(0, main([
                "multi-reference-locus-union", "--primary", str(zw6),
                "--primary-label", "ZW6", "--alternate", f"CameorV2={v2}",
                "--mapping", bridge["mapping"], "--out-dir", str(out), "--json",
            ]))
            compression = pd.read_csv(out / "compression_candidates.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({
                "pissa.ZW6.gnm1.ann1.T1",
                "pissa.Cameor.gnm1.ann1.S4",
            }, set(compression.gene_id))
            pav = compression[compression.gene_id.eq("pissa.Cameor.gnm1.ann1.S4")].iloc[0]
            self.assertEqual("CameorV2", pav.assembly_label)
            self.assertEqual("PAV", pav.special_status)


if __name__ == "__main__":
    unittest.main()
