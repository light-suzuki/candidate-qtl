"""Contract tests for the Medicago LIS resource generations used by #524."""

from __future__ import annotations

from pathlib import Path
import unittest

from candidate_qtl.reference_resources import load_reference_resource_registry


ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "config" / "medicago_lis_reference_registry.yml"


class MedicagoReferenceRegistryTests(unittest.TestCase):
    def test_expected_generations_and_namespaces_are_pinned(self) -> None:
        resources = load_reference_resource_registry(REGISTRY)
        by_id = {row["resource_id"]: row for row in resources}

        self.assertEqual(
            {
                "medtr.A17_HM341.gnm4.ann2.G3ZY",
                "medtr.A17.gnm5.ann1_6.L2RX",
                "medtr.R108.gnmHiC_1.ann1.Y8NH",
            },
            set(by_id),
        )
        self.assertEqual("medtr.A17_HM341.gnm4.ann2", by_id["medtr.A17_HM341.gnm4.ann2.G3ZY"]["native_gene_id_namespace"])
        self.assertEqual("medtr.A17.gnm5.ann1_6", by_id["medtr.A17.gnm5.ann1_6.L2RX"]["native_gene_id_namespace"])
        self.assertEqual("medtr.R108.gnmHiC_1.ann1", by_id["medtr.R108.gnmHiC_1.ann1.Y8NH"]["native_gene_id_namespace"])

    def test_resources_are_optional_and_r108_is_explicit_mutant_reference(self) -> None:
        resources = load_reference_resource_registry(REGISTRY)
        by_id = {row["resource_id"]: row for row in resources}

        self.assertTrue(all(row["required"].casefold() == "false" for row in resources))
        self.assertEqual("near_reference", by_id["medtr.A17_HM341.gnm4.ann2.G3ZY"]["role"])
        self.assertEqual("near_reference", by_id["medtr.A17.gnm5.ann1_6.L2RX"]["role"])
        self.assertEqual("mutant_reference", by_id["medtr.R108.gnmHiC_1.ann1.Y8NH"]["role"])

    def test_declared_auxiliary_files_stay_generation_scoped(self) -> None:
        resources = load_reference_resource_registry(REGISTRY)
        for row in resources:
            dataset_name = row["resource_id"].split(".", 1)[1]
            key4 = row["dataset_key4"]
            self.assertTrue(row["manifest_path"].startswith(dataset_name + "/"))
            self.assertTrue(row["checksum_path"].startswith(dataset_name + "/"))
            self.assertIn(key4, row["fam3_path"])
            self.assertIn(key4, row["protein_path"])
            self.assertIn(key4, row["gff3_path"])

        r108 = next(row for row in resources if row["accession_or_ecotype"] == "R108")
        self.assertIn("Y8NH", r108["id_map_path"])
        self.assertTrue(r108["id_map_path"].endswith(".id_map.txt.gz"))
        self.assertEqual("auxiliary_published_artifact", r108["id_map_role"])
        self.assertEqual("not_evaluated", r108["id_map_identity_scope"])


if __name__ == "__main__":
    unittest.main()
