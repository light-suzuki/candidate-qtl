"""Contract tests for the Lotus Gifu/MG20 resource generations used by #525."""

from __future__ import annotations

from pathlib import Path
import unittest

from candidate_qtl.reference_resources import load_reference_resource_registry


ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "config" / "lotus_reference_registry.yml"
GIFU_ID = "lotusbase.Gifu_B-129.Gifu_v1_2.ann_v1_3"
MG20_ID = "lotja.MG20.gnm3.ann1.WF9B"


class LotusReferenceRegistryTests(unittest.TestCase):
    def _resources(self) -> dict[str, dict[str, str]]:
        return {row["resource_id"]: row for row in load_reference_resource_registry(REGISTRY)}

    def test_gifu_assembly_and_annotation_generations_remain_separate(self) -> None:
        gifu = self._resources()[GIFU_ID]
        self.assertEqual("external", gifu["namespace_scheme"])
        self.assertEqual("mutant_reference", gifu["role"])
        self.assertEqual("Gifu B-129", gifu["accession_or_ecotype"])
        self.assertEqual("Gifu_v1.2", gifu["assembly_id"])
        self.assertEqual("Gifu_v1.3", gifu["annotation_id"])
        self.assertEqual("Lotusjaponicus_Gifu_v1.2_genome.fa.gz", gifu["assembly_filename"])
        self.assertEqual("ec94ee57478d395b9dc37bbabaa96767", gifu["assembly_md5"])
        self.assertEqual("Lotusjaponicus_Gifu_v1.3_cds.fa.gz", gifu["cds_filename"])
        self.assertEqual("39fd2b5aa169e8eb18d55e4350fbc33c", gifu["cds_md5"])
        self.assertEqual("a297c83c5e213c1ca5a287316caccb36", gifu["protein_md5"])
        self.assertEqual("444d7a31a4d70b0136974a1fa3aa3204", gifu["gff3_md5"])

    def test_gifu_upstream_filenames_have_explicit_runtime_paths(self) -> None:
        gifu = self._resources()[GIFU_ID]
        payloads = ("assembly", "cds", "protein", "gff3", "go", "transcript_mapping")
        for payload in payloads:
            self.assertEqual(gifu[f"{payload}_filename"], gifu[f"{payload}_path"], payload)
        self.assertTrue(all(gifu[f"{payload}_path"] for payload in payloads))

    def test_transcript_mapping_is_candidate_correspondence_not_identity(self) -> None:
        gifu = self._resources()[GIFU_ID]
        self.assertEqual("Lotusjaponicus_TranscriptMapping.csv.gz", gifu["transcript_mapping_filename"])
        self.assertEqual("10ede41cb6356c5cbb6eed7b912d1c74", gifu["transcript_mapping_md5"])
        self.assertEqual("MG20_v3.0", gifu["transcript_mapping_source_generation"])
        self.assertEqual("Gifu_v1.2", gifu["transcript_mapping_target_generation"])
        self.assertEqual("bidirectional_blastn_candidate_correspondence", gifu["transcript_mapping_relation"])
        self.assertEqual("candidate_only_not_verified", gifu["transcript_mapping_identity_scope"])
        self.assertIn("one_to_many", gifu["transcript_mapping_multiplicity_policy"])

    def test_lore1_background_is_not_rewritten_by_mg20_coordinate_lookup(self) -> None:
        gifu = self._resources()[GIFU_ID]
        self.assertEqual("Gifu", gifu["lore1_experimental_background"])
        self.assertEqual(
            "coordinate_lookup_not_background_identity",
            gifu["lore1_mg20_coordinate_mapping_role"],
        )
        self.assertEqual("preserve_exact_upstream_case", gifu["native_id_case_policy"])
        self.assertEqual("explicit_mapping_required_no_symbol_join", gifu["cross_resource_identity_rule"])

    def test_lis_mg20_generation_and_payloads_are_explicit(self) -> None:
        mg20 = self._resources()[MG20_ID]
        self.assertEqual("lis_full_id", mg20["namespace_scheme"])
        self.assertEqual("lotja.MG20.gnm3.ann1", mg20["native_gene_id_namespace"])
        self.assertEqual("gnm3", mg20["assembly_id"])
        self.assertEqual("ann1", mg20["annotation_id"])
        self.assertEqual("WF9B", mg20["dataset_key4"])
        self.assertEqual("MG20.gnm3.QPGB", mg20["lis_genome_dataset"])
        self.assertEqual("GCA_000181115.2", mg20["lis_genome_genbank_accession"])
        for key in ("manifest_path", "checksum_path", "fam3_path", "protein_path", "gff3_path", "interpro_path"):
            self.assertIn("MG20.gnm3.ann1.WF9B", mg20[key])
        self.assertTrue(mg20["fam3_path"].endswith(".legume.fam3.VLMQ.gfa.tsv.gz"))
        self.assertTrue(mg20["protein_path"].endswith(".protein.faa.gz"))
        self.assertTrue(mg20["gff3_path"].endswith(".gene_models_main.gff3.gz"))

    def test_lotusbase_mg20_and_lis_mg20_are_not_joined_by_version_number_alone(self) -> None:
        mg20 = self._resources()[MG20_ID]
        self.assertEqual("MG20_v3.0", mg20["lotusbase_comparative_generation"])
        self.assertEqual("08bdeb9d026debc04c35672063fd6959", mg20["lotusbase_genome_md5"])
        self.assertEqual("e501f31cf850d74cb8c51ba7f684f4bf", mg20["lotusbase_protein_md5"])
        self.assertEqual("baafdee6708503c328ad9e2382e88259", mg20["lotusbase_gff3_md5"])
        self.assertEqual(
            "explicit_metadata_or_mapping_required_not_numeric_name_inference",
            mg20["generation_compatibility_rule"],
        )
        self.assertEqual("preserve_exact_upstream_case", mg20["native_id_case_policy"])
        self.assertEqual("explicit_mapping_required_no_symbol_join", mg20["cross_resource_identity_rule"])

    def test_both_resources_are_optional_until_materialized(self) -> None:
        resources = self._resources()
        self.assertEqual({GIFU_ID, MG20_ID}, set(resources))
        self.assertTrue(all(row["required"].casefold() == "false" for row in resources.values()))


if __name__ == "__main__":
    unittest.main()
