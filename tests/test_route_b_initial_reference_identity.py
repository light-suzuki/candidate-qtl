"""Initial candidate correspondence respects each frozen reference endpoint."""

from __future__ import annotations

import copy
import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from candidate_qtl.route_b_semantic_bundle import build_route_b_semantic_bundle
from candidate_qtl.route_b_semantic_candidates import freeze_route_b_semantic_candidates
from candidate_qtl.route_b_semantic_import import import_route_b_semantic_proposals
from candidate_qtl.route_b_semantic_source_verify import verify_route_b_semantic_sources
from candidate_qtl.route_b_semantic_transfer_verify import (
    _ReferenceLink, _initial_reference_members, _reference_link, verify_route_b_semantic_transfer,
)
from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff


ARABIDOPSIS = "Arabidopsis thaliana"
MEDICAGO = "Medicago truncatula"


class RouteBInitialReferenceIdentityTests(unittest.TestCase):
    def _chain(self, root, species_values, *, literature_species=None, literature_namespace=None,
               reference_gene="REF1", gene_endpoint=True, paper_reports=None, copies=1,
               endpoint_species=ARABIDOPSIS, endpoint_taxon_id="", record_cache_arg=None):
        atlas = root / "atlas.tsv"
        gene = {
            "gene_id": "G1", "family_id": "F1",
            "local_annotation_description": "meristem signaling protein",
        }
        if gene_endpoint:
            gene.update(orthology_reference_gene=reference_gene, orthology_reference_species=endpoint_species,
                        orthology_relationship="one_to_many")
            if endpoint_taxon_id:
                gene.update(reference_gene=reference_gene, reference_species=endpoint_species,
                            reference_taxon_id=endpoint_taxon_id)
        pd.DataFrame([{**gene, "gene_id": f"G{i + 1}"} for i in range(copies)]).to_csv(atlas, sep="\t", index=False)
        literature = root / "literature.tsv"
        article = {"interval_gene_ids": ";".join(f"G{i + 1}" for i in range(copies)), "reference_gene": reference_gene, "source_id": "PMID:1234567",
                   "title": "Synthetic meristem observation", "reference_function": "regulation of meristem maintenance"}
        if literature_species is not None:
            article["reference_species"] = literature_species
        if literature_namespace is not None:
            article["reference_accession_namespace"] = literature_namespace
        pd.DataFrame([article]).to_csv(literature, sep="\t", index=False)
        bundle_result = build_route_b_semantic_bundle(
            str(atlas), str(root / "bundle"), phenotype="altered meristem development",
            literature_arg=str(literature), case_id="initial-reference-fixture", focal_species="Solanum lycopersicum",
        )
        bundle_path = Path(bundle_result["bundle"])
        bundle = json.loads(bundle_path.read_text())
        prompt = json.loads(Path(bundle_result["prompts"]).read_text())
        claims = [{
            "claim_type": "mechanistic_review_seed", "biological_function": "regulation of meristem maintenance",
            "mechanism": "altered stem cell signaling could explain altered meristem development",
            "phenotype_relation": "conditional cross-species meristem hypothesis", "reference_gene": reference_gene,
            "reference_species": species, "evidence_scope": "reference_gene", "source_ids": ["PMID:1234567"],
            "alternative_hypothesis": "another signaling protein could explain the phenotype",
            "uncertainty_note": "reference observation and focal causal hypothesis remain distinct",
            "initial_tier": "T3", "initial_tier_reason": "conditional homolog process with unverified focal causality",
        } for species in species_values]
        if paper_reports is not None:
            for claim in claims:
                claim["paper_reports_json"] = paper_reports
        response = root / "responses.json"
        response.write_text(json.dumps({"candidates": [{"prompt_id": prompt["prompt_id"], "parse_status": "json",
                                                       "model_json": {"claims": claims}}]}))
        imported = import_route_b_semantic_proposals(
            bundle_result["bundle"], bundle_result["prompts"], str(response), str(root / "import"),
        )
        cache = root / "cache.json"
        cache.write_text(json.dumps({"id": "1234567", "title": "Synthetic meristem observation"}))
        sources = verify_route_b_semantic_sources(
            imported["proposals"], str(root / "source"), bundle_arg=bundle_result["bundle"], cache_arg=str(cache),
        )
        transfer = verify_route_b_semantic_transfer(
            imported["proposals"], sources["source_verification"], bundle_result["bundle"],
            str(root / "transfer"), record_cache_arg=record_cache_arg, case_id=bundle["case_id"],
        )
        frozen = freeze_route_b_semantic_candidates(transfer["claims"], str(root / "freeze"), case_id=bundle["case_id"])
        originals = {path: path.read_bytes() for path in (bundle_path, response, Path(transfer["claims"]))}
        wet = wet_candidate_handoff(
            str(root / "wet"), route_b_candidates_arg=frozen["candidates"], route_b_candidates_manifest_arg=frozen["manifest"],
            review_candidates_arg=frozen["review_candidates"], review_candidates_manifest_arg=frozen["review_manifest"],
            interval_atlas_arg=str(atlas), case_id=bundle["case_id"],
        )
        for path, original in originals.items():
            self.assertEqual(original, path.read_bytes())
        return bundle, transfer, frozen, wet

    def test_wrong_species_is_ledger_only_and_independent_same_gene_claim_reaches_wet(self):
        for species_values, expected in (([MEDICAGO], set()), ([MEDICAGO, ARABIDOPSIS], {"G1"})):
            with self.subTest(species_values=species_values), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                _, transfer, frozen, wet = self._chain(root, species_values)
                final = pd.read_csv(transfer["claims"], sep="\t", dtype=str).fillna("")
                self.assertEqual(species_values, list(final.reference_species))
                self.assertEqual({"T3"}, set(final.initial_tier))
                self.assertEqual(0, frozen["candidate_gene_count"])
                assessed = pd.read_csv(transfer["transfer_verification"], sep="\t", dtype=str).fillna("")
                wrong_id = final.loc[final.reference_species.eq(MEDICAGO), "claim_id"].iloc[0]
                wrong = assessed[assessed.claim_id.eq(wrong_id)].iloc[0]
                self.assertEqual("literature_unit_membership", wrong.reference_link_basis)
                self.assertEqual("not_evaluated", wrong.review_eligibility)
                self.assertEqual("", wrong.review_member_genes)
                review = pd.read_csv(frozen["review_candidates"], sep="\t", dtype=str).fillna("")
                self.assertEqual(expected, set(review.candidate_gene_id))
                self.assertNotIn(wrong_id, set(review.claim_id))
                if wet["wet_candidates"] is None:
                    self.assertEqual(set(), expected)
                else:
                    candidates = pd.read_csv(wet["wet_candidates"], sep="\t", dtype=str).fillna("")
                    self.assertEqual(expected, set(candidates.candidate_gene_id))
                    inventory = pd.read_csv(root / "wet/wet_relevance_inventory.tsv", sep="\t", dtype=str).fillna("")
                    self.assertEqual({"G1"}, set(inventory.candidate_gene_id))
                ledger = pd.read_csv(root / "freeze/route_b_review_eligibility.tsv", sep="\t", dtype=str).fillna("")
                self.assertEqual({"not_asserted"}, set(ledger.biological_exclusion))

    def test_literature_only_unknown_identity_matches_missing_without_changing_strict_link(self):
        cases = [
            ("species_missing", "REF1", None, None, {"G1"}),
            ("species_unknown", "REF1", "unknown", None, {"G1"}),
            ("species_conflict", "REF1", MEDICAGO, None, set()),
            ("namespace_missing", "NS1|REF1|ENTRY", ARABIDOPSIS, None, {"G1"}),
            ("namespace_unknown", "NS1|REF1|ENTRY", ARABIDOPSIS, "unknown", {"G1"}),
            ("namespace_conflict", "NS1|REF1|ENTRY", ARABIDOPSIS, "NS2", set()),
        ]
        for name, reference, literature_species, literature_namespace, expected in cases:
            with self.subTest(case=name), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                bundle, transfer, frozen, wet = self._chain(
                    root, [ARABIDOPSIS], literature_species=literature_species,
                    literature_namespace=literature_namespace, reference_gene=reference, gene_endpoint=False,
                )
                unit = bundle["units"][0]
                article = unit["literature_link_candidates"][0]
                for key, value in (("reference_species", literature_species),
                                   ("reference_accession_namespace", literature_namespace)):
                    if value is None:
                        self.assertNotIn(key, article)
                    else:
                        self.assertEqual(value, article[key])
                self.assertEqual([], unit["reference_material_cards"])
                self.assertFalse(unit["genes"][0].get("orthology_reference_gene"))
                review = pd.read_csv(frozen["review_candidates"], sep="\t", dtype=str).fillna("")
                self.assertEqual(expected, set(review.candidate_gene_id))
                self.assertEqual(0, frozen["candidate_gene_count"])
                assessed = pd.read_csv(transfer["transfer_verification"], sep="\t", dtype=str).fillna("")
                self.assertEqual("eligible_initial_candidate" if expected else "not_evaluated", assessed.iloc[0].review_eligibility)
                final = pd.read_csv(transfer["claims"], sep="\t", dtype=str).fillna("")
                self.assertEqual("T3", final.iloc[0].initial_tier)
                self.assertEqual("not_evaluated", final.iloc[0].verification_status)
                if expected:
                    self.assertEqual({"G1"}, set(pd.read_csv(wet["wet_candidates"], sep="\t", dtype=str).candidate_gene_id))
                else:
                    self.assertTrue(pd.read_csv(wet["wet_candidates"], sep="\t", dtype=str).empty)
                    manifest = json.loads(Path(wet["manifest"]).read_text())
                    self.assertEqual(0, manifest["candidate_count_populations"]["final_wet"])
                    self.assertEqual("review_required_no_supported_candidates", manifest["handoff_state"])
                # The source/support matcher keeps the strict replay contract.
                # Only the initial membership path treats explicit unknowns as
                # missing, and its group fallback still requires the same paper.
                strict_link = _reference_link(reference, ARABIDOPSIS, "PMID:1234567", unit, unit["genes"])
                self.assertEqual("literature_unit_membership" if name.endswith("_missing") else "none", strict_link.basis)
                wrong_source = final.iloc[0].copy()
                wrong_source["source_id"] = "PMID:7654321"
                self.assertEqual(set(), _initial_reference_members(wrong_source, unit, _ReferenceLink((), "none", "none")))

    def test_unknown_literature_group_cannot_override_explicit_gene_species(self):
        genes = [{"gene_id": "G1", "orthology_reference_gene": "REF1", "orthology_reference_species": ARABIDOPSIS}]
        unit = {"genes": genes, "literature_link_candidates": [{
            "reference_gene": "REF1", "source_id": "PMID:1234567", "reference_species": "unknown",
        }]}
        claim = pd.Series({"reference_gene": "REF1", "reference_species": MEDICAGO,
                           "source_id": "PMID:1234567", "member_interval_genes": "G1"})
        link = _reference_link("REF1", MEDICAGO, "PMID:1234567", unit, genes)
        self.assertEqual("none", link.basis)
        self.assertEqual(set(), _initial_reference_members(claim, unit, link))

    def test_declared_literature_metadata_is_retained_and_legacy_transfer_stays_separate(self):
        for policy in ("initial_bibliography_v1", "source_function_v1"):
            with self.subTest(policy=policy), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                pd.DataFrame([{"gene_id": "G1", "family_id": "F1", "local_annotation_description": "meristem signaling"}]).to_csv(root / "atlas.tsv", sep="\t", index=False)
                pd.DataFrame([{"interval_gene_ids": "G1", "reference_gene": "REF1", "source_id": "PMID:1234567",
                               "reference_species": ARABIDOPSIS, "reference_accession_namespace": "fixture"}]).to_csv(root / "lit.tsv", sep="\t", index=False)
                result = build_route_b_semantic_bundle(str(root / "atlas.tsv"), str(root / "bundle"),
                    phenotype="altered development", literature_arg=str(root / "lit.tsv"), b_review_membership_policy=policy)
                unit = json.loads(Path(result["bundle"]).read_text())["units"][0]
                article = unit["literature_link_candidates"][0]
                self.assertEqual(ARABIDOPSIS, article["reference_species"])
                self.assertEqual("fixture", article["reference_accession_namespace"])
                self.assertEqual("literature_unit_membership", _reference_link("REF1", ARABIDOPSIS, "PMID:1234567", unit, unit["genes"]).basis)
                self.assertEqual("none", _reference_link("REF1", MEDICAGO, "PMID:1234567", unit, unit["genes"]).basis)
                # Old already-frozen records keep their existing group policy;
                # only new initial-member correspondence adds the conflict guard.
                old_unit = copy.deepcopy(unit)
                old_unit["literature_link_candidates"][0].pop("reference_species")
                self.assertEqual("literature_unit_membership", _reference_link("REF1", MEDICAGO, "PMID:1234567", old_unit, old_unit["genes"]).basis)

    def test_same_reference_copies_alternatives_and_unknown_metadata_are_distinct(self):
        genes = [
            {"gene_id": "wrong", "orthology_reference_gene": "REF1", "orthology_reference_species": ARABIDOPSIS},
            {"gene_id": "matching", "orthology_reference_gene": "REF1", "orthology_reference_species": MEDICAGO},
            {"gene_id": "unknown", "orthology_reference_gene": "REF1", "orthology_reference_species": "unknown"},
            {"gene_id": "other-lookup", "orthology_reference_gene": "REF1", "orthology_reference_species": ARABIDOPSIS,
             "family_index_subject_endpoints": "REF1", "reference_gene": "OTHER", "reference_species": MEDICAGO},
            {"gene_id": "matching-hit", "orthology_reference_gene": "REF1", "orthology_reference_species": ARABIDOPSIS,
             "sequence_similarity_hits": [{"hit_accession": "REF1", "hit_species": MEDICAGO}]},
            {"gene_id": "matching-card", "orthology_reference_gene": "REF1", "orthology_reference_species": ARABIDOPSIS},
        ]
        unit = {"genes": genes, "reference_material_cards": [
            {"interval_gene_id": "matching-card", "reference_gene": "REF1", "reference_species": MEDICAGO},
            {"interval_gene_id": "wrong", "reference_gene": "REF1", "reference_species": "unknown"},
        ]}
        claim = pd.Series({"reference_gene": "REF1", "reference_species": MEDICAGO,
                           "member_interval_genes": ";".join(gene["gene_id"] for gene in genes)})
        link = _ReferenceLink(tuple(genes), "literature_unit_membership", "group")
        expected = {"matching", "unknown", "matching-hit", "matching-card"}
        self.assertEqual(expected, _initial_reference_members(claim, unit, link))
        unit["genes"].reverse()
        unit["reference_material_cards"].reverse()
        self.assertEqual(expected, _initial_reference_members(claim, unit, link))

    def test_other_reference_species_does_not_borrow_a_neighboring_copy(self):
        genes = [{"gene_id": "G1", "orthology_reference_gene": "REF1", "orthology_reference_species": ARABIDOPSIS},
                 {"gene_id": "G2", "orthology_reference_gene": "REF2", "orthology_reference_species": MEDICAGO}]
        unit = {"genes": genes, "literature_link_candidates": [{"reference_gene": "REF2", "source_id": "PMID:1234567"}]}
        claim = pd.Series({"reference_gene": "REF2", "reference_species": MEDICAGO, "member_interval_genes": "G1;G2"})
        link = _reference_link("REF2", MEDICAGO, "PMID:1234567", unit, genes)
        self.assertEqual({"G2"}, _initial_reference_members(claim, unit, link))

    def test_namespace_conflict_is_not_overridden_by_untyped_family_handle(self):
        reference = "NS1|REF1|ENTRY"
        genes = [{"gene_id": "G1", "orthology_reference_gene": reference,
                  "orthology_reference_species": ARABIDOPSIS, "orthology_reference_accession_namespace": "NS2",
                  "family_index_subject_endpoints": reference}]
        unit = {"genes": genes, "reference_material_cards": []}
        claim = pd.Series({"reference_gene": reference, "reference_species": ARABIDOPSIS, "member_interval_genes": "G1"})
        link = _ReferenceLink(tuple(genes), "literature_unit_membership", "group")
        self.assertEqual(set(), _initial_reference_members(claim, unit, link))
        genes[0]["orthology_reference_accession_namespace"] = "unknown"
        self.assertEqual({"G1"}, _initial_reference_members(claim, unit, link))


if __name__ == "__main__":
    unittest.main()
