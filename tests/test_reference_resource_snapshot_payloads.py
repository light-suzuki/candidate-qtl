"""Contract tests for full declared reference-resource snapshot freezing."""

from __future__ import annotations

import hashlib
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from candidate_qtl.reference_resources import freeze_reference_resource_snapshots


DATASET = "A17.gnm5.ann1_6.L2RX"
RESOURCE_ID = f"medtr.{DATASET}"
PAYLOAD_NAMES = {
    "manifest": f"MANIFEST.{DATASET}.yml",
    "checksum": f"CHECKSUM.{DATASET}.md5",
    "fam3": f"medtr.{DATASET}.legume.fam3.VLMQ.gfa.tsv.gz",
    "protein": f"medtr.{DATASET}.protein.faa.gz",
    "gff3": f"medtr.{DATASET}.gene_models_main.gff3.gz",
}


def _registry_text() -> str:
    return f"""resources:
  - resource_id: {RESOURCE_ID}
    species: Medicago truncatula
    accession_or_ecotype: A17
    assembly_id: gnm5
    annotation_id: ann1_6
    dataset_key4: L2RX
    native_gene_id_namespace: medtr.A17.gnm5.ann1_6
    namespace_scheme: lis_full_id
    role: near_reference
    required: false
    manifest_path: {DATASET}/{PAYLOAD_NAMES['manifest']}
    checksum_path: {DATASET}/{PAYLOAD_NAMES['checksum']}
    fam3_path: {DATASET}/{PAYLOAD_NAMES['fam3']}
    protein_path: {DATASET}/{PAYLOAD_NAMES['protein']}
    gff3_path: {DATASET}/{PAYLOAD_NAMES['gff3']}
    source_snapshot_date: 2026-08-22
    source_url_or_upstream_id: https://data.legumeinfo.org/Medicago/truncatula/annotations/{DATASET}/
    license_or_redistribution_note: upstream LIS terms apply
"""


def _materialize_fixture(root: Path, *, missing: set[str] | None = None) -> dict[str, bytes]:
    missing = missing or set()
    dataset_dir = root / DATASET
    dataset_dir.mkdir(parents=True)
    contents: dict[str, bytes] = {}
    for payload, filename in PAYLOAD_NAMES.items():
        data = f"fixture:{payload}\n".encode()
        contents[payload] = data
        if payload not in missing:
            (dataset_dir / filename).write_bytes(data)
    return contents


class FullReferenceResourceSnapshotTests(unittest.TestCase):
    def test_all_declared_payloads_receive_full_sha256(self) -> None:
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            runtime_root = root / "runtime"
            expected = _materialize_fixture(runtime_root)
            registry = root / "registry.yml"
            registry.write_text(_registry_text(), encoding="utf-8")

            summary = freeze_reference_resource_snapshots(
                registry,
                root / "out",
                runtime_root_arg=runtime_root,
            )

            resource = summary["resources"][0]
            self.assertEqual("available", resource["availability"])
            self.assertEqual("complete", resource["snapshot_completeness"])
            self.assertEqual(set(PAYLOAD_NAMES), set(resource["payloads"]))
            self.assertEqual(len(PAYLOAD_NAMES), resource["declared_payload_count"])
            self.assertEqual(len(PAYLOAD_NAMES), resource["frozen_payload_count"])
            self.assertEqual(resource["payloads"]["manifest"], resource["manifest"])
            self.assertEqual(resource["payloads"]["checksum"], resource["checksum"])
            for payload, data in expected.items():
                frozen = resource["payloads"][payload]
                self.assertEqual("available", frozen["status"])
                self.assertEqual("full_sha256", frozen["method"])
                self.assertEqual(hashlib.sha256(data).hexdigest(), frozen["sha256"])

            self.assertEqual(1, summary["counts"]["full_snapshot_complete"])
            self.assertEqual(len(PAYLOAD_NAMES), summary["counts"]["declared_payloads"])
            self.assertEqual(len(PAYLOAD_NAMES), summary["counts"]["frozen_payloads"])

    def test_missing_non_manifest_payload_is_partial_without_changing_legacy_availability(self) -> None:
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            runtime_root = root / "runtime"
            _materialize_fixture(runtime_root, missing={"protein"})
            registry = root / "registry.yml"
            registry.write_text(_registry_text(), encoding="utf-8")

            summary = freeze_reference_resource_snapshots(
                registry,
                root / "out",
                runtime_root_arg=runtime_root,
            )

            resource = summary["resources"][0]
            self.assertEqual("available", resource["availability"])
            self.assertEqual("partial", resource["snapshot_completeness"])
            self.assertEqual("unavailable", resource["payloads"]["protein"]["status"])
            self.assertIsNone(resource["payloads"]["protein"]["sha256"])
            self.assertEqual(len(PAYLOAD_NAMES) - 1, resource["frozen_payload_count"])
            self.assertEqual(1, summary["counts"]["full_snapshot_partial"])
            self.assertEqual(0, summary["counts"]["required_blocked"])

    def test_relative_payloads_without_runtime_root_remain_not_evaluated(self) -> None:
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            registry = root / "registry.yml"
            registry.write_text(_registry_text(), encoding="utf-8")

            summary = freeze_reference_resource_snapshots(registry, root / "out")

            resource = summary["resources"][0]
            self.assertEqual("not_evaluated", resource["availability"])
            self.assertEqual("not_evaluated", resource["snapshot_completeness"])
            self.assertTrue(all(item["status"] == "not_evaluated" for item in resource["payloads"].values()))
            self.assertEqual(1, summary["counts"]["full_snapshot_not_evaluated"])
            self.assertEqual(0, summary["counts"]["frozen_payloads"])


if __name__ == "__main__":
    unittest.main()
