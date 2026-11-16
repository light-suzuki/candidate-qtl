"""Regression coverage for the #524 LIS pangene candidate-correspondence bridge."""

from __future__ import annotations

import gzip
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.assembly_bridge import build_assembly_bridge
from candidate_qtl.cli import main


FUNCTIONS = "\n".join([
    "native_gene_id\tphenotype_or_function\tsource_ids\tnative_species",
    "medtr.A17_HM341.gnm4.ann2.S1\tcurated function one\tPMID:1\tMedicago truncatula",
    "medtr.A17_HM341.gnm4.ann2.S2\tcurated function two\tPMID:2\tMedicago truncatula",
    "medtr.A17_HM341.gnm4.ann2.S3\tcurated function three\tPMID:3\tMedicago truncatula",
    "medtr.A17_HM341.gnm4.ann2.S4\tcurated function four\tPMID:4\tMedicago truncatula",
]) + "\n"

SOURCE_FAM3 = "\n".join([
    "native_gene_id\tfamily_id",
    "medtr.A17_HM341.gnm4.ann2.S1\tF1",
    "medtr.A17_HM341.gnm4.ann2.S2\tF2",
    "medtr.A17_HM341.gnm4.ann2.S3\tF3",
    "medtr.A17_HM341.gnm4.ann2.S4\tF5",
]) + "\n"

TARGET_FAM3 = "\n".join([
    "native_gene_id\tfamily_id",
    "medtr.A17.gnm5.ann1_6.T1\tF1",
    # Same family as T1 but absent from the pangene set: once pangene supplies a
    # cross-side candidate, this decoy must not be added merely by fam3.
    "medtr.A17.gnm5.ann1_6.T1_DECOY\tF1",
    "medtr.A17.gnm5.ann1_6.T2\tF2",
    "medtr.A17.gnm5.ann1_6.T3\tF2",
    "medtr.A17.gnm5.ann1_6.T4\tF4",
    # S4/T5 are deliberately absent from pangene: fam3 fallback must remain so
    # absence from the pangene snapshot never becomes negative evidence.
    "medtr.A17.gnm5.ann1_6.T5\tF5",
]) + "\n"

SEQUENCE = "\n".join([
    "source_full_id\ttarget_full_id\treciprocal_sequence_state\tsequence_identity\tquery_coverage\tsubject_coverage",
    "medtr.A17_HM341.gnm4.ann2.S1\tmedtr.A17.gnm5.ann1_6.T1\trbh\t98.2\t1.0\t1.0",
]) + "\n"

PANGENE = "\n".join([
    "gene\trepresentative",
    "medtr.A17_HM341.gnm4.ann2.S1\tmedtr.A17.gnm5.ann1_6.T1",
    "medtr.A17.gnm5.ann1_6.T1\tmedtr.A17.gnm5.ann1_6.T1",
    "medtr.A17_HM341.gnm4.ann2.S2\tPAN_SPLIT",
    "medtr.A17.gnm5.ann1_6.T2\tPAN_SPLIT",
    "medtr.A17.gnm5.ann1_6.T3\tPAN_SPLIT",
    "medtr.A17_HM341.gnm4.ann2.S3\tPAN_FAMILY_MISMATCH",
    "medtr.A17.gnm5.ann1_6.T4\tPAN_FAMILY_MISMATCH",
]) + "\n"


class MedicagoPangeneBridgeTests(unittest.TestCase):
    def _inputs(self, root: Path) -> dict[str, Path]:
        paths = {}
        for name, text in (("functions.tsv", FUNCTIONS), ("source_fam3.tsv", SOURCE_FAM3),
                           ("target_fam3.tsv", TARGET_FAM3), ("sequence.tsv", SEQUENCE)):
            path = root / name
            path.write_text(text, encoding="utf-8")
            paths[name] = path
        pangene = root / "Medicago.pan2.451K.hsh.tsv.gz"
        with gzip.open(pangene, "wt", encoding="utf-8", newline="") as handle:
            handle.write(PANGENE)
        paths["pangene"] = pangene
        return paths

    def _run(self, root: Path) -> tuple[pd.DataFrame, dict[str, object], Path]:
        paths = self._inputs(root)
        out = root / "out"
        result = build_assembly_bridge(
            paths["functions.tsv"], paths["source_fam3.tsv"], paths["target_fam3.tsv"],
            paths["sequence.tsv"], out,
            source_namespace="medtr.A17_HM341.gnm4.ann2",
            target_namespace="medtr.A17.gnm5.ann1_6",
            mapping_scope="same_accession_annotation_generation",
            pangene_hash_arg=paths["pangene"], pangene_dataset="Medicago.pan2.451K",
            case_id="medicago_generation_bridge",
        )
        mapping = pd.read_csv(result["mapping"], sep="\t", dtype=str).fillna("")
        manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
        return mapping, manifest, out

    def test_pangene_hash_records_candidates_without_weakening_verification(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            mapping, manifest, out = self._run(Path(tmp))

            exact_rows = mapping[mapping.source_full_id.eq("medtr.A17_HM341.gnm4.ann2.S1")]
            # Pangene is a candidate-generation layer: the same-fam3 decoy is
            # not considered once an explicit cross-side pangene relation exists.
            self.assertEqual({"medtr.A17.gnm5.ann1_6.T1"}, set(exact_rows.target_full_id))
            exact = exact_rows.iloc[0]
            self.assertEqual("Medicago.pan2.451K", exact.pangene_dataset)
            self.assertEqual("true", exact.pangene_relation_present)
            self.assertEqual("one_to_one_candidate", exact.pangene_relationship)
            self.assertEqual("verified", exact.mapping_state)
            self.assertEqual("one_to_one", exact.relationship)
            self.assertIn("pangene_correspondence", exact.mapping_basis)
            self.assertIn("pairwise_sequence", exact.mapping_basis)

            split = mapping[mapping.source_full_id.eq("medtr.A17_HM341.gnm4.ann2.S2")]
            self.assertEqual({"medtr.A17.gnm5.ann1_6.T2", "medtr.A17.gnm5.ann1_6.T3"}, set(split.target_full_id))
            self.assertTrue(split.mapping_state.eq("ambiguous").all())
            self.assertTrue(split.relationship.eq("one_to_many").all())
            self.assertTrue(split.pangene_relationship.eq("one_to_many_candidate").all())
            self.assertTrue(split.source_member_count_in_pangene.eq("1").all())
            self.assertTrue(split.target_member_count_in_pangene.eq("2").all())

            mismatch = mapping[
                mapping.source_full_id.eq("medtr.A17_HM341.gnm4.ann2.S3")
                & mapping.target_full_id.eq("medtr.A17.gnm5.ann1_6.T4")
            ].iloc[0]
            self.assertEqual("provisional", mismatch.mapping_state)
            self.assertIn("pangene_correspondence", mismatch.mapping_basis)
            self.assertIn("fam3_mismatch", mismatch.mapping_basis)
            self.assertIn("identity_not_verified", mismatch.mapping_reason)

            # No pangene relation for S4: old fam3 candidate discovery survives
            # as a fail-open fallback, but without reciprocal sequence it remains
            # provisional rather than being promoted to an exact identity.
            fallback = mapping[
                mapping.source_full_id.eq("medtr.A17_HM341.gnm4.ann2.S4")
                & mapping.target_full_id.eq("medtr.A17.gnm5.ann1_6.T5")
            ].iloc[0]
            self.assertEqual("false", fallback.pangene_relation_present)
            self.assertEqual("provisional", fallback.mapping_state)
            self.assertEqual("one_to_one", fallback.relationship)
            self.assertEqual("same_fam3_without_reciprocal_sequence_evidence", fallback.mapping_reason)

            # Only the reciprocal, family-concordant 1:1 mapping may create an
            # exact gene-level function transfer. Pangene-only/mismatch/fallback
            # provisional rows stay outside the copy-level Route B adapter.
            reference = pd.read_csv(out / "route_b_reference_function.tsv", sep="\t", dtype=str).fillna("")
            exact_function_rows = reference[reference.reference_gene_id.eq("medtr.A17_HM341.gnm4.ann2.S1")]
            self.assertEqual({"medtr.A17.gnm5.ann1_6.T1"}, set(exact_function_rows.gene_id))
            mismatch_function_rows = reference[reference.reference_gene_id.eq("medtr.A17_HM341.gnm4.ann2.S3")]
            fallback_function_rows = reference[reference.reference_gene_id.eq("medtr.A17_HM341.gnm4.ann2.S4")]
            self.assertTrue(mismatch_function_rows.empty)
            self.assertTrue(fallback_function_rows.empty)

            self.assertEqual("Medicago.pan2.451K", manifest["pangene_dataset"])
            coverage = manifest["pangene_namespace_coverage"]
            self.assertTrue(coverage["source_namespace_represented"])
            self.assertTrue(coverage["target_namespace_represented"])
            self.assertEqual(4, coverage["source_gene_count"])
            self.assertEqual(6, coverage["target_gene_count"])
            self.assertEqual(3, coverage["source_gene_in_pangene_count"])
            self.assertEqual(4, coverage["target_gene_in_pangene_count"])
            self.assertEqual("none", manifest["rank_contribution"])

    def test_verified_function_card_can_chain_old_a17_to_current_a17_to_r108(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _mapping, _manifest, first_out = self._run(root)

            current_fam3 = root / "current_fam3.tsv"
            current_fam3.write_text(
                "native_gene_id\tfamily_id\n"
                "medtr.A17.gnm5.ann1_6.T1\tF1\n",
                encoding="utf-8",
            )
            r108_fam3 = root / "r108_fam3.tsv"
            r108_fam3.write_text(
                "native_gene_id\tfamily_id\n"
                "medtr.R108.gnmHiC_1.ann1.R1\tF1\n",
                encoding="utf-8",
            )
            current_to_r108_sequence = root / "current_to_r108_sequence.tsv"
            current_to_r108_sequence.write_text(
                "source_full_id\ttarget_full_id\treciprocal_sequence_state\tsequence_identity\tquery_coverage\tsubject_coverage\n"
                "medtr.A17.gnm5.ann1_6.T1\tmedtr.R108.gnmHiC_1.ann1.R1\trbh\t97.5\t1.0\t1.0\n",
                encoding="utf-8",
            )
            pangene = root / "current_r108.hsh.tsv.gz"
            with gzip.open(pangene, "wt", encoding="utf-8", newline="") as handle:
                handle.write(
                    "gene\trepresentative\n"
                    "medtr.A17.gnm5.ann1_6.T1\tPAN_CHAIN\n"
                    "medtr.R108.gnmHiC_1.ann1.R1\tPAN_CHAIN\n"
                )

            second_out = root / "r108_out"
            second = build_assembly_bridge(
                first_out / "route_b_reference_function.tsv",
                current_fam3,
                r108_fam3,
                current_to_r108_sequence,
                second_out,
                source_namespace="medtr.A17.gnm5.ann1_6",
                target_namespace="medtr.R108.gnmHiC_1.ann1",
                mapping_scope="same_species_cross_assembly",
                pangene_hash_arg=pangene,
                pangene_dataset="Medicago.pan2.451K",
                case_id="old_a17_to_current_a17_to_r108",
            )

            second_mapping = pd.read_csv(second["mapping"], sep="\t", dtype=str).fillna("")
            chained = second_mapping[
                second_mapping.source_full_id.eq("medtr.A17.gnm5.ann1_6.T1")
                & second_mapping.target_full_id.eq("medtr.R108.gnmHiC_1.ann1.R1")
            ].iloc[0]
            self.assertEqual("verified", chained.mapping_state)
            self.assertEqual("one_to_one", chained.relationship)

            r108_functions = pd.read_csv(second["route_b_reference_function"], sep="\t", dtype=str).fillna("")
            self.assertEqual({"medtr.R108.gnmHiC_1.ann1.R1"}, set(r108_functions.gene_id))
            self.assertEqual({"medtr.A17.gnm5.ann1_6.T1"}, set(r108_functions.reference_gene_id))
            self.assertEqual({"curated function one"}, set(r108_functions.function_description))
            self.assertEqual({"PMID:1"}, set(r108_functions.source_id))
            self.assertEqual({"exact_gene"}, set(r108_functions.function_transfer_scope))

            # The ambiguous S2 old-A17 card has no exact current-A17 gene_id in
            # the first adapter, so it cannot silently enter the R108 bridge.
            self.assertNotIn("curated function two", set(r108_functions.function_description))

    def test_pangene_dataset_generation_is_required(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = self._inputs(root)
            with self.assertRaisesRegex(ValueError, "pangene_dataset"):
                build_assembly_bridge(
                    paths["functions.tsv"], paths["source_fam3.tsv"], paths["target_fam3.tsv"],
                    paths["sequence.tsv"], root / "out", pangene_hash_arg=paths["pangene"],
                )

    def test_headerless_lis_hsh_set_id_first_is_consumed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = self._inputs(root)
            native = root / "native.hsh.tsv.gz"
            with gzip.open(native, "wt", encoding="utf-8", newline="") as handle:
                # Native LIS hsh files are headerless: set ID first, gene second.
                handle.write(
                    "PAN_NATIVE\tmedtr.A17_HM341.gnm4.ann2.S1\n"
                    "PAN_NATIVE\tmedtr.A17.gnm5.ann1_6.T1\n"
                )
            out = root / "native_out"
            result = build_assembly_bridge(
                paths["functions.tsv"], paths["source_fam3.tsv"], paths["target_fam3.tsv"],
                paths["sequence.tsv"], out,
                source_namespace="medtr.A17_HM341.gnm4.ann2",
                target_namespace="medtr.A17.gnm5.ann1_6",
                pangene_hash_arg=native,
                pangene_dataset="Medicago.pan2.451K",
            )
            mapping = pd.read_csv(result["mapping"], sep="\t", dtype=str).fillna("")
            row = mapping[
                mapping.source_full_id.eq("medtr.A17_HM341.gnm4.ann2.S1")
                & mapping.target_full_id.eq("medtr.A17.gnm5.ann1_6.T1")
            ].iloc[0]
            self.assertEqual("PAN_NATIVE", row.pangene_set_id)
            self.assertEqual("true", row.pangene_relation_present)
            self.assertEqual("verified", row.mapping_state)

    def test_missing_target_mapping_is_not_evaluated_not_negative(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = self._inputs(root)
            functions = root / "functions_with_unresolved.tsv"
            functions.write_text(
                FUNCTIONS + "medtr.A17_HM341.gnm4.ann2.S5\tcurated function five\tPMID:5\tMedicago truncatula\n",
                encoding="utf-8",
            )
            out = root / "unresolved_out"
            result = build_assembly_bridge(
                functions, paths["source_fam3.tsv"], paths["target_fam3.tsv"], paths["sequence.tsv"], out,
                source_namespace="medtr.A17_HM341.gnm4.ann2",
                target_namespace="medtr.A17.gnm5.ann1_6",
                pangene_hash_arg=paths["pangene"],
                pangene_dataset="Medicago.pan2.451K",
            )
            mapping = pd.read_csv(result["mapping"], sep="\t", dtype=str).fillna("")
            unresolved = mapping[mapping.source_full_id.eq("medtr.A17_HM341.gnm4.ann2.S5")].iloc[0]
            self.assertEqual("", unresolved.target_full_id)
            self.assertEqual("unresolved", unresolved.relationship)
            self.assertEqual("not_evaluated", unresolved.mapping_state)
            self.assertIn("absence_is_not_negative_evidence", unresolved.mapping_reason)

    def test_cli_wires_frozen_pangene_hash_and_dataset(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = self._inputs(root)
            out = root / "cli_out"
            self.assertEqual(0, main([
                "assembly-bridge",
                "--functional-evidence", str(paths["functions.tsv"]),
                "--source-fam3", str(paths["source_fam3.tsv"]),
                "--target-fam3", str(paths["target_fam3.tsv"]),
                "--sequence-evidence", str(paths["sequence.tsv"]),
                "--pangene-hash", str(paths["pangene"]),
                "--pangene-dataset", "Medicago.pan2.451K",
                "--source-namespace", "medtr.A17_HM341.gnm4.ann2",
                "--target-namespace", "medtr.A17.gnm5.ann1_6",
                "--mapping-scope", "same_accession_annotation_generation",
                "--out-dir", str(out),
            ]))
            manifest = json.loads((out / "assembly_bridge_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("Medicago.pan2.451K", manifest["pangene_dataset"])
            self.assertIn("pangene_hash", manifest["inputs"])


if __name__ == "__main__":
    unittest.main()
