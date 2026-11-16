from __future__ import annotations

import json
import hashlib
import copy
import tempfile
import unittest
from datetime import datetime
from pathlib import Path

import pandas as pd
import yaml

from candidate_qtl import __version__
from candidate_qtl.reference_resources import (
    build_legume_reference_snapshot,
    freeze_reference_resource_snapshots,
    load_reference_resource_registry,
    parse_lis_full_gene_id,
)


REGISTRY = """resources:
  - resource_id: medtr-a17
    species: Medicago truncatula
    accession_or_ecotype: A17
    assembly_id: gnm5
    annotation_id: ann1_6
    native_gene_id_namespace: medtr.A17.gnm5.ann1_6
    namespace_scheme: lis_full_id
    role: near_reference
    source_snapshot_date: 2026-08-20
    source_url_or_upstream_id: LIS:A17.gnm5.ann1_6
    license_or_redistribution_note: upstream terms apply
  - resource_id: medtr-r108
    species: Medicago truncatula
    accession_or_ecotype: R108
    assembly_id: gnmHiC_1
    annotation_id: ann1
    native_gene_id_namespace: medtr.R108.gnmHiC_1.ann1
    namespace_scheme: lis_full_id
    role: mutant_reference
    source_snapshot_date: 2026-08-20
    source_url_or_upstream_id: LIS:R108.gnmHiC_1.ann1
    license_or_redistribution_note: upstream terms apply
  - resource_id: lotja-gifu
    species: Lotus japonicus
    accession_or_ecotype: Gifu
    assembly_id: Gifu_v1.2
    annotation_id: Gifu_v1.3
    native_gene_id_namespace: lotusbase.Gifu.v1.2-v1.3
    namespace_scheme: external
    role: mutant_reference
    source_snapshot_date: 2026-08-20
    source_url_or_upstream_id: LotusBase:Gifu-v1.2-assembly-v1.3-annotation
    license_or_redistribution_note: upstream terms apply
"""

TRAITS = """---
scientific_name: Medicago truncatula
gene_symbols:
  - MtSAME
gene_model_pub_name: MtrunA17Chr1g0000011
gene_model_full_id: medtr.A17.gnm5.ann1_6.MtrunA17Chr1g0000011
confidence: 5
phenotype_synopsis: Required for nodulation in the tested background.
traits:
  - entity_name: nodulation
    entity: GO:0009877
references:
  - citation: Example, 2026
    doi: 10.1000/a17
    pmid: 100001
---
scientific_name: Medicago truncatula
gene_symbols:
  - MtSAME
gene_model_pub_name: R108Gene0001
gene_model_full_id: medtr.R108.gnmHiC_1.ann1.R108Gene0001
confidence: 4
phenotype_synopsis: Mutant plants show an altered root phenotype.
traits:
  - entity_name: root morphology
    entity: TO:0000043
references:
  - citation: Example2, 2026
    doi: 10.1000/r108
    pmid: 100002
---
scientific_name: Pisum sativum
gene_symbols:
  - PsUNREGISTERED
gene_model_pub_name: Psat1g000001
gene_model_full_id: pissa.Cameor.gnm1.ann1.Psat1g000001
confidence: 3
phenotype_synopsis: A valid full-ID record from a namespace not in this tiny registry.
traits:
  - entity_name: shoot development
    entity: TO:0000928
references:
  - citation: Example3, 2026
    doi: 10.1000/pea
---
scientific_name: Medicago truncatula
gene_symbols:
  - MtBROKEN
confidence: 3
phenotype_synopsis: This card lacks a full LIS gene identifier and must not be guessed.
traits:
  - entity_name: nodulation
    entity: GO:0009877
references:
  - citation: Example4, 2026
    doi: 10.1000/broken
"""

FAMILY = """#gene\tfamily\tprotein\te-value\tscore\tbest_domain_score\tpct_hmm_coverage
medtr.A17.gnm5.ann1_6.MtrunA17Chr1g0000011\tlegume.fam3.VLMQ.TEST1\tmedtr.A17.gnm5.ann1_6.MtrunA17Chr1g0000011.1\t1e-80\t300.0\t299.5\t98.5
medtr.R108.gnmHiC_1.ann1.R108Gene0001\tlegume.fam3.VLMQ.TEST1\tmedtr.R108.gnmHiC_1.ann1.R108Gene0001.1\t2e-75\t280.0\t270.0\t62.0
"""


class LegumeReferenceResourceTests(unittest.TestCase):
    def test_full_lis_id_parser_requires_explicit_namespace(self) -> None:
        parsed = parse_lis_full_gene_id("medtr.A17.gnm5.ann1_6.MtrunA17Chr1g0000011")
        self.assertIsNotNone(parsed)
        assert parsed is not None
        self.assertEqual("medtr.A17.gnm5.ann1_6", parsed["namespace"])
        self.assertEqual("A17", parsed["accession"])
        self.assertIsNone(parse_lis_full_gene_id("MtSAME"))
        self.assertIsNone(parse_lis_full_gene_id("medtr.A17.MtrunA17Chr1g0000011"))

    def test_registry_rejects_lis_namespace_metadata_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            path = Path(tempdir) / "registry.yml"
            path.write_text(REGISTRY.replace("accession_or_ecotype: A17", "accession_or_ecotype: R108", 1), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "accession does not match namespace"):
                load_reference_resource_registry(path)

    def test_registry_accepts_external_native_namespace_without_lis_guessing(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            path = Path(tempdir) / "registry.yml"
            path.write_text(REGISTRY, encoding="utf-8")
            resources = load_reference_resource_registry(path)
            gifu = next(row for row in resources if row["resource_id"] == "lotja-gifu")
            self.assertEqual("external", gifu["namespace_scheme"])
            self.assertEqual("Gifu_v1.2", gifu["assembly_id"])
            self.assertEqual("Gifu_v1.3", gifu["annotation_id"])
            self.assertEqual("lotusbase.Gifu.v1.2-v1.3", gifu["native_gene_id_namespace"])
            self.assertIsNone(parse_lis_full_gene_id(gifu["native_gene_id_namespace"] + ".LjNIN"))

    def test_snapshot_preserves_accession_and_source_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            registry = root / "registry.yml"
            traits = root / "medtr.traits.yml"
            family = root / "medtr.gfa.tsv"
            out = root / "out"
            registry.write_text(REGISTRY, encoding="utf-8")
            traits.write_text(TRAITS, encoding="utf-8")
            family.write_text(FAMILY, encoding="utf-8")

            result = build_legume_reference_snapshot(str(registry), str(traits), str(family), str(out), full_hash=True)
            functions = pd.read_csv(out / "normalized_functional_evidence.tsv", sep="\t", dtype=str).fillna("")
            families = pd.read_csv(out / "family_assignment.tsv", sep="\t", dtype=str).fillna("")
            audit = pd.read_csv(out / "resource_ingest_audit.tsv", sep="\t", dtype=str).fillna("")
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))

            # Same literature symbol in two accessions must remain two native identities.
            same_symbol = functions[functions["gene_symbol"] == "MtSAME"]
            self.assertEqual(2, len(same_symbol))
            self.assertEqual({"A17", "R108"}, set(same_symbol["native_accession"]))
            self.assertEqual(2, same_symbol["native_gene_id"].nunique())

            a17 = same_symbol[same_symbol["native_accession"] == "A17"].iloc[0]
            self.assertEqual("5", a17["source_confidence"])
            self.assertIn("PMID:100001", a17["source_ids"])
            self.assertIn("DOI:10.1000/a17", a17["source_ids"])
            self.assertIn("GO:0009877", a17["ontology_ids"])
            self.assertEqual("registered_namespace", a17["resource_match_state"])

            # A valid but unregistered LIS namespace is preserved and surfaced, not dropped.
            pea = functions[functions["gene_symbol"] == "PsUNREGISTERED"].iloc[0]
            self.assertEqual("unregistered_namespace", pea["resource_match_state"])
            self.assertEqual("Cameor", pea["native_accession"])

            # External namespaces are registry metadata only for these LIS-specific adapters.
            self.assertIn("lotusbase.Gifu.v1.2-v1.3", manifest["registered_namespaces"])
            self.assertNotIn("lotusbase.Gifu.v1.2-v1.3", manifest["registered_lis_namespaces"])

            # Missing full IDs are auditable unresolved records and are never guessed from symbols.
            unresolved = audit[(audit["source_kind"] == "gene_function") & (audit["status"] == "unresolved")]
            self.assertEqual(1, len(unresolved))
            self.assertEqual("missing_or_invalid_full_lis_gene_id", unresolved.iloc[0]["reason"])
            self.assertNotIn("MtBROKEN", set(functions["gene_symbol"]))

            # Official GFA quality fields survive normalization for later fragment/fusion review.
            self.assertEqual(2, len(families))
            r108_family = families[families["native_accession"] == "R108"].iloc[0]
            self.assertEqual("legume.fam3.VLMQ.TEST1", r108_family["family_id"])
            self.assertEqual("62.0", r108_family["pct_hmm_coverage"])
            self.assertEqual("270.0", r108_family["best_domain_score"])

            self.assertFalse(manifest["cross_accession_identity_inferred"])
            self.assertEqual("none", manifest["rank_contribution"])
            self.assertEqual("2026-08-20", manifest["resources"][0]["source_snapshot_date"])
            self.assertEqual("medtr.A17.gnm5.ann1_6", manifest["resources"][0]["native_gene_id_namespace"])
            self.assertEqual("lis_full_id", manifest["resources"][0]["namespace_scheme"])
            self.assertEqual(3, manifest["resource_count"])
            self.assertEqual(3, manifest["counts"]["functional_evidence_records"])
            self.assertEqual(0, manifest["counts"]["functional_duplicate_native_gene_record_count"])
            self.assertEqual(2, manifest["counts"]["family_assignment_records"])
            self.assertEqual(2, manifest["counts"]["function_genes_with_family_assignment"])
            self.assertEqual(1, manifest["counts"]["functional_unregistered_namespace_records"])
            self.assertEqual(1, manifest["counts"]["audit_unresolved_or_invalid_records"])
            self.assertEqual("full_sha256", manifest["inputs"]["gene_function_registry"]["method"])

    def test_pisum_lis_registry_preserves_four_key4_dataset_identities(self) -> None:
        registry_path = Path(__file__).parents[1] / "config" / "pisum_lis_reference_registry.yml"
        resources = load_reference_resource_registry(registry_path)
        self.assertEqual(4, len(resources))
        self.assertEqual(
            {
                "pissa.Cameor.gnm1.ann1.7SZR",
                "pissa.Cameor.gnm2.ann1.KT6X",
                "pissa.JI2822.gnm1.ann1.D0QS",
                "pissa.ZW6.gnm1.ann1.TKZX",
            },
            {row["resource_id"] for row in resources},
        )
        self.assertEqual({"7SZR", "KT6X", "D0QS", "TKZX"}, {row["dataset_key4"] for row in resources})
        self.assertEqual(
            {
                "pissa.Cameor.gnm1.ann1",
                "pissa.Cameor.gnm2.ann1",
                "pissa.JI2822.gnm1.ann1",
                "pissa.ZW6.gnm1.ann1",
            },
            {row["native_gene_id_namespace"] for row in resources},
        )

    def test_runtime_snapshot_freezes_files_independently_and_keeps_optional_missing(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            registry_payload = yaml.safe_load(
                (Path(__file__).parents[1] / "config" / "pisum_lis_reference_registry.yml").read_text(encoding="utf-8")
            )
            registry_path = root / "registry.yml"
            registry_path.write_text(yaml.safe_dump(registry_payload, sort_keys=False), encoding="utf-8")
            runtime_root = root / "runtime"
            cameor_v1_dir = runtime_root / "Cameor.gnm1.ann1.7SZR"
            cameor_v1_dir.mkdir(parents=True)
            manifest_one = cameor_v1_dir / "MANIFEST.Cameor.gnm1.ann1.7SZR.yml"
            checksum_one = cameor_v1_dir / "CHECKSUM.Cameor.gnm1.ann1.7SZR.md5"
            manifest_one.write_text("manifest-one\n", encoding="utf-8")
            checksum_one.write_text("checksum-one\n", encoding="utf-8")
            cameor_v2_dir = runtime_root / "Cameor.gnm2.ann1.KT6X"
            cameor_v2_dir.mkdir(parents=True)
            (cameor_v2_dir / "MANIFEST.Cameor.gnm2.ann1.KT6X.yml").write_text("manifest-two\n", encoding="utf-8")

            result = freeze_reference_resource_snapshots(registry_path, root / "out", runtime_root_arg=runtime_root)
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            by_id = {row["resource_id"]: row for row in manifest["resources"]}
            cameor_v1 = by_id["pissa.Cameor.gnm1.ann1.7SZR"]
            self.assertEqual("available", cameor_v1["availability"])
            self.assertEqual("available", cameor_v1["manifest"]["status"])
            self.assertEqual("available", cameor_v1["checksum"]["status"])
            self.assertEqual(hashlib.sha256(manifest_one.read_bytes()).hexdigest(), cameor_v1["manifest"]["sha256"])
            self.assertEqual(hashlib.sha256(checksum_one.read_bytes()).hexdigest(), cameor_v1["checksum"]["sha256"])

            cameor_v2 = by_id["pissa.Cameor.gnm2.ann1.KT6X"]
            self.assertEqual("unavailable", cameor_v2["availability"])
            self.assertEqual("available", cameor_v2["manifest"]["status"])
            self.assertEqual("unavailable", cameor_v2["checksum"]["status"])
            zw6 = by_id["pissa.ZW6.gnm1.ann1.TKZX"]
            self.assertEqual("unavailable", zw6["availability"])
            self.assertEqual("unavailable", zw6["manifest"]["status"])
            self.assertEqual(1, manifest["counts"]["available"])
            self.assertEqual(3, manifest["counts"]["unavailable"])
            self.assertEqual(0, manifest["counts"]["not_evaluated"])
            self.assertEqual(3, manifest["counts"]["optional_unavailable_or_not_evaluated"])
            self.assertEqual("ready_available_resources_frozen", manifest["state"])
            self.assertFalse(manifest["cross_accession_identity_inferred"])
            self.assertEqual("full_sha256", manifest["registry"]["method"])
            self.assertEqual(hashlib.sha256(registry_path.read_bytes()).hexdigest(), manifest["registry"]["sha256"])
            self.assertEqual("1.1", manifest["schema_version"])
            self.assertEqual("1.1", manifest["registry_schema_version"])
            self.assertEqual("reference-resource-snapshot-freezer-1.1", manifest["parser_version"])
            self.assertEqual(__version__, manifest["candidate_qtl_version"])
            frozen_at = datetime.fromisoformat(manifest["frozen_at_utc"])
            self.assertIsNotNone(frozen_at.utcoffset())
            self.assertEqual("direct_target", cameor_v1["role"])
            self.assertEqual("2026-08-21", cameor_v1["source_snapshot_date"])
            self.assertIn("data.legumeinfo.org", cameor_v1["source_url_or_upstream_id"])

            checksum_one.unlink()
            blocked = freeze_reference_resource_snapshots(registry_path, root / "blocked", runtime_root_arg=runtime_root)
            self.assertEqual("blocked_required_resource_missing", blocked["state"])
            self.assertEqual(1, blocked["counts"]["required_blocked"])

    def test_pisum_dataset_key4_duplicate_or_mismatch_fails_closed(self) -> None:
        fixture = Path(__file__).parents[1] / "config" / "pisum_lis_reference_registry.yml"
        payload = yaml.safe_load(fixture.read_text(encoding="utf-8"))
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            duplicate = copy.deepcopy(payload)
            duplicate["resources"][1]["dataset_key4"] = "7SZR"
            duplicate["resources"][1]["resource_id"] = "pissa.Cameor.gnm2.ann1.7SZR"
            duplicate_path = root / "duplicate.yml"
            duplicate_path.write_text(yaml.safe_dump(duplicate, sort_keys=False), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "duplicate dataset_key4"):
                load_reference_resource_registry(duplicate_path)

            mismatch = copy.deepcopy(payload)
            mismatch["resources"][0]["dataset_key4"] = "ZZZZ"
            mismatch_path = root / "mismatch.yml"
            mismatch_path.write_text(yaml.safe_dump(mismatch, sort_keys=False), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "resource_id does not match dataset_key4"):
                load_reference_resource_registry(mismatch_path)

            wrong_file = copy.deepcopy(payload)
            wrong_file["resources"][0]["manifest_path"] = "Cameor.gnm1.ann1.7SZR/MANIFEST.ZW6.gnm1.ann1.TKZX.yml"
            wrong_file_path = root / "wrong-file.yml"
            wrong_file_path.write_text(yaml.safe_dump(wrong_file, sort_keys=False), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "manifest_path does not match dataset_key4"):
                load_reference_resource_registry(wrong_file_path)

            wrong_parent = copy.deepcopy(payload)
            expected_name = "MANIFEST.Cameor.gnm1.ann1.7SZR.yml"
            wrong_parent["resources"][0]["manifest_path"] = f"ZW6.gnm1.ann1.TKZX/{expected_name}"
            wrong_parent_path = root / "wrong-parent.yml"
            wrong_parent_path.write_text(yaml.safe_dump(wrong_parent, sort_keys=False), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "manifest_path does not match dataset_key4"):
                load_reference_resource_registry(wrong_parent_path)

            missing_required = copy.deepcopy(payload)
            del missing_required["resources"][0]["required"]
            missing_required_path = root / "missing-required.yml"
            missing_required_path.write_text(yaml.safe_dump(missing_required, sort_keys=False), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "requires explicit boolean required"):
                load_reference_resource_registry(missing_required_path)

    def test_relative_runtime_files_without_root_are_not_evaluated_not_cwd_relative(self) -> None:
        registry = Path(__file__).parents[1] / "config" / "pisum_lis_reference_registry.yml"
        with tempfile.TemporaryDirectory() as tempdir:
            result = freeze_reference_resource_snapshots(registry, Path(tempdir) / "out")
            self.assertEqual("blocked_required_resource_missing", result["state"])
            first = result["resources"][0]
            self.assertEqual("not_evaluated", first["availability"])
            self.assertEqual("runtime_root_required_for_relative_path", first["manifest"]["reason"])

    def test_optional_only_registry_with_no_available_resource_is_not_ready(self) -> None:
        fixture = Path(__file__).parents[1] / "config" / "pisum_lis_reference_registry.yml"
        payload = yaml.safe_load(fixture.read_text(encoding="utf-8"))
        payload["resources"][0]["role"] = "near_reference"
        payload["resources"][0]["required"] = False
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            registry = root / "registry.yml"
            registry.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
            result = freeze_reference_resource_snapshots(registry, root / "out", runtime_root_arg=root / "empty")
            self.assertEqual("not_ready_no_available_resources", result["state"])
            self.assertEqual(0, result["counts"]["required_blocked"])


if __name__ == "__main__":
    unittest.main()
