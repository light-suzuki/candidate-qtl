"""Acquired atlas material survives lookup overlays without reference mixing."""

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.route_b_semantic_bundle import build_route_b_semantic_bundle
from candidate_qtl.route_b_semantic_transfer_verify import verify_route_b_semantic_transfer
from candidate_qtl.semantic_claims import CLAIM_COLUMNS


class RouteBLookupMaterialPreservationTests(unittest.TestCase):
    def _bundle(self, root, atlas_rows, lookup_rows, *, policy="initial_bibliography_v1"):
        atlas, lookup = root / "atlas.tsv", root / "lookup.tsv"
        pd.DataFrame(atlas_rows).fillna("").to_csv(atlas, sep="\t", index=False)
        pd.DataFrame(lookup_rows).fillna("").to_csv(lookup, sep="\t", index=False)
        result = build_route_b_semantic_bundle(
            str(atlas), str(root / "out"), phenotype="altered development",
            lookup_arg=str(lookup), b_review_membership_policy=policy,
        )
        bundle = json.loads(Path(result["bundle"]).read_text(encoding="utf-8"))
        prompts = [json.loads(line) for line in Path(result["prompts"]).read_text(encoding="utf-8").splitlines()]
        return atlas, result, bundle, prompts

    def test_empty_and_placeholder_lookup_do_not_erase_acquired_function(self):
        for empty in ("", "unknown", "not_evaluated", "[truncated]"):
            with self.subTest(lookup_function=empty), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                atlas, result, bundle, prompts = self._bundle(root, [{
                    "gene_id": "G1", "family_id": "F001", "family_label": "F001",
                    "reference_function": "measured receptor activity",
                    "unrelated_raw_provider_text": "unrelated provider dump " * 100,
                }], [{"gene_id": "G1", "reference_gene": "REF1", "reference_function": empty}])
                self.assertEqual(1, result["prompt_count"])
                unit = bundle["units"][0]
                self.assertEqual("function_context", unit["source_material_class"])
                self.assertEqual("requested", unit["request_state"])
                self.assertEqual("not_evaluated", unit["review_status"])
                primary = unit["genes"][0]
                self.assertEqual(empty, primary["reference_function"])
                context = primary["atlas_material_context"]
                self.assertEqual("G1", context["gene_id"])
                self.assertEqual("measured receptor activity", context["reference_function"])
                self.assertNotIn("reference_gene", context)
                self.assertNotIn("unrelated_raw_provider_text", context)
                self.assertEqual(str(atlas), context["origin"]["path"])
                self.assertEqual(hashlib.sha256(atlas.read_bytes()).hexdigest(), context["origin"]["sha256"])
                self.assertEqual("full_sha256", context["origin"]["method"])
                self.assertEqual("investigation_context_not_verified_support", context["material_status"])
                self.assertEqual("none", context["rank_contribution"])
                delivered = json.loads(prompts[0]["user"].split("\n\n", 1)[1])
                self.assertEqual(context, delivered["genes"][0]["atlas_material_context"])
                self.assertNotIn("unrelated provider dump", prompts[0]["user"])
                self.assertEqual([], unit["stable_source_ids"])
                self.assertEqual(0, unit["source_backed_function_record_count"])

    def test_lookup_does_not_borrow_other_reference_species_or_citation_fields(self):
        atlas_row = {
            "gene_id": "G1", "family_id": "F001", "family_label": "F001",
            "reference_gene": "ATLAS_REF", "reference_species": "Atlas species",
            "reference_function": "measured receptor activity",
            "reference_source_id": "PMID:111", "reference_source_url": "https://example.test/111",
        }
        for lookup_row in (
            {"gene_id": "G1", "reference_gene": "OTHER_REF", "reference_source_id": "PMID:222"},
            {"gene_id": "G1", "reference_gene": "ATLAS_REF", "reference_species": "Other species", "reference_source_id": "PMID:222"},
            {"gene_id": "G1", "reference_gene": "ATLAS_REF", "reference_species": "Atlas species", "reference_source_id": "PMID:222"},
        ):
            with self.subTest(lookup=lookup_row), tempfile.TemporaryDirectory() as tmp:
                _, result, bundle, prompts = self._bundle(Path(tmp), [atlas_row], [lookup_row])
                self.assertEqual(1, result["prompt_count"])
                unit = bundle["units"][0]
                primary = unit["genes"][0]
                for field in ("reference_gene", "reference_species", "reference_function", "reference_source_id", "reference_source_url"):
                    self.assertEqual(lookup_row.get(field, ""), primary[field])
                    self.assertEqual(atlas_row[field], primary["atlas_material_context"][field])
                self.assertEqual(["PMID:222"], unit["stable_source_ids"])
                self.assertEqual(0, unit["source_backed_function_record_count"])
                self.assertEqual("OTHER_REF" if "OTHER_REF" in lookup_row.values() else "ATLAS_REF",
                                 unit["source_provenance"][0]["reference_gene"])
                delivered = json.loads(prompts[0]["user"].split("\n\n", 1)[1])
                self.assertEqual(primary, delivered["genes"][0])

    def test_other_focal_gene_lookup_cannot_supply_or_erase_material(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, result, bundle, prompts = self._bundle(Path(tmp), [
                {"gene_id": "G1", "reference_function": "measured receptor activity", "reference_gene": "ATLAS_REF"},
                {"gene_id": "G2"},
            ], [{"gene_id": "G2", "reference_gene": "OTHER_REF", "reference_function": "unknown"}])
            self.assertEqual(1, result["prompt_count"])
            units = {unit["member_interval_genes"][0]: unit for unit in bundle["units"]}
            g1 = units["G1"]["genes"][0]
            self.assertEqual("ATLAS_REF", g1["reference_gene"])
            self.assertEqual("measured receptor activity", g1["reference_function"])
            self.assertNotIn("atlas_material_context", g1)
            self.assertEqual("not_requested_no_reviewable_material", units["G2"]["request_state"])
            self.assertNotIn("measured receptor activity", json.dumps(units["G2"]))
            self.assertNotIn("OTHER_REF", prompts[0]["user"])

    def test_unchanged_reference_material_and_focal_annotations_are_not_duplicated(self):
        lookup_row = {"gene_id": "G1", "reference_gene": "REF1",
                      "reference_function": "measured receptor activity", "reference_source_id": "PMID:111"}
        for atlas_row in (
            {**lookup_row, "family_id": "F001", "family_label": "F001"},
            {"gene_id": "G1", "local_annotation_description": "annotated focal receptor"},
        ):
            with self.subTest(atlas=atlas_row), tempfile.TemporaryDirectory() as tmp:
                _, result, bundle, prompts = self._bundle(Path(tmp), [atlas_row], [lookup_row])
                self.assertEqual(1, result["prompt_count"])
                primary = bundle["units"][0]["genes"][0]
                self.assertNotIn("atlas_material_context", primary)
                self.assertEqual("measured receptor activity", primary["reference_function"])
                self.assertEqual(atlas_row.get("local_annotation_description", ""), primary["local_annotation_description"])
                self.assertEqual(1, bundle["units"][0]["source_backed_function_record_count"])
                self.assertNotIn("atlas_material_context", prompts[0]["user"])

    def test_atlas_endpoint_species_constraint_survives_independent_lookup(self):
        for atlas_reference in ("", "REF1"):
            with self.subTest(atlas_reference=atlas_reference), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                _, result, bundle, _ = self._bundle(root, [{
                    "gene_id": "G1", "orthology_reference_gene": "REF1",
                    "reference_gene": atlas_reference, "reference_species": "Atlas species",
                    "reference_accession_namespace": "UniProt", "orthology_relationship": "one_to_one",
                    "reference_function": "measured receptor activity",
                }], [{"gene_id": "G1", "reference_gene": "OTHER_REF", "reference_function": ""}],
                    policy="source_function_v1")
                rows = []
                for claim_id, species in (("same", "Atlas species"), ("wrong", "Other species")):
                    row = {key: "" for key in CLAIM_COLUMNS}
                    row.update(claim_id=claim_id, route="interval_first", claim_type="mechanistic_process",
                               biological_function="receptor activity", mechanism="signaling",
                               phenotype_relation="conditional altered development",
                               semantic_unit_id=bundle["units"][0]["unit_id"],
                               member_interval_genes="G1", reference_gene="REF1", reference_species=species,
                               evidence_scope="interval_gene", verification_status="proposed", rank_contribution="none")
                    rows.append(row)
                proposals, source = root / "proposals.tsv", root / "source.tsv"
                pd.DataFrame(rows, columns=CLAIM_COLUMNS).to_csv(proposals, sep="\t", index=False)
                pd.DataFrame([{"claim_id": row["claim_id"], "source_relation_status": "not_evaluated"} for row in rows]).to_csv(source, sep="\t", index=False)
                transfer_result = verify_route_b_semantic_transfer(
                    str(proposals), str(source), result["bundle"], str(root / "transfer"))
                transfers = pd.read_csv(transfer_result["transfer_verification"], sep="\t", dtype=str).fillna("").set_index("claim_id")
                self.assertEqual("direct_gene_orthology", transfers.loc["same", "reference_link_basis"])
                self.assertEqual("none", transfers.loc["wrong", "reference_link_basis"])
                self.assertTrue(transfers["final_verification_status"].eq("not_evaluated").all())
                primary = bundle["units"][0]["genes"][0]
                self.assertEqual("Atlas species", primary["orthology_reference_species"])
                self.assertEqual("UniProt", primary["orthology_reference_accession_namespace"])
                self.assertEqual("", primary["reference_species"])
                self.assertEqual("", primary["reference_accession_namespace"])

    def test_distinct_atlas_reference_metadata_is_not_moved_to_orthology_endpoint(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, _, bundle, _ = self._bundle(Path(tmp), [{
                "gene_id": "G1", "orthology_reference_gene": "ORTHO_REF",
                "reference_gene": "ATLAS_REF", "reference_species": "Atlas species",
                "reference_accession_namespace": "Atlas namespace",
                "reference_function": "measured receptor activity",
            }], [{"gene_id": "G1", "reference_gene": "LOOKUP_REF",
                 "reference_species": "Lookup species", "reference_accession_namespace": "Lookup namespace"}])
            primary = bundle["units"][0]["genes"][0]
            self.assertEqual("", primary["orthology_reference_species"])
            self.assertEqual("", primary["orthology_reference_accession_namespace"])
            self.assertEqual("Lookup species", primary["reference_species"])
            self.assertEqual("Lookup namespace", primary["reference_accession_namespace"])
            context = primary["atlas_material_context"]
            self.assertEqual("ATLAS_REF", context["reference_gene"])
            self.assertEqual("Atlas species", context["reference_species"])
            self.assertEqual("Atlas namespace", context["reference_accession_namespace"])


if __name__ == "__main__":
    unittest.main()
