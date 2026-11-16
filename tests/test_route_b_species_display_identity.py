"""Source-bound species display checks; all genes/claims are synthetic."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from candidate_qtl.claim_review_session import prepare_session, complete_session
from candidate_qtl.llm_pathway_legacy import _load_record_cache
from candidate_qtl.reference_taxonomy import load_taxonomy_cache
from candidate_qtl.route_b_semantic_transfer_verify import _ReferenceLink, _initial_reference_members
from candidate_qtl.route_b_semantic_candidates import route_b_homolog_projection, freeze_route_b_semantic_candidates
from tests.test_route_b_initial_reference_identity import RouteBInitialReferenceIdentityTests

SCIENTIFIC = "Medicago truncatula"
DISPLAY = "Medicago truncatula (Barrel medic) (Medicago tribuloides)"
TAXONOMY = {"taxonId": 3880, "scientificName": SCIENTIFIC,
            "commonName": "Barrel medic", "synonyms": ["Medicago tribuloides"]}


class SpeciesDisplayTests(unittest.TestCase):
    def members(self, genes, *, cards=None, species=SCIENTIFIC, records=None, reference="G7KFE4"):
        claim = pd.Series({"reference_gene": reference, "reference_species": species,
                           "member_interval_genes": ";".join(gene["gene_id"] for gene in genes)})
        unit = {"genes": genes, "reference_material_cards": cards or []}
        original = copy.deepcopy((claim.to_dict(), unit, records))
        audit = []
        result = _initial_reference_members(claim, unit, _ReferenceLink((), "none", "none"),
            taxonomy_records=records, display_audit=audit)
        self.assertEqual(original, (claim.to_dict(), unit, records))
        return result, audit

    def gene(self, **extra):
        return {"gene_id": "G1", "reference_gene": "G7KFE4", "reference_species": DISPLAY,
                "reference_taxon_id": "3880", **extra}

    def test_exact_taxon_scientific_name_and_official_display_required(self):
        records = {"3880": TAXONOMY}
        self.assertEqual(set(), self.members([self.gene()])[0])
        result, audit = self.members([self.gene()], records=records)
        self.assertEqual({"G1"}, result)
        self.assertEqual("3880", audit[0]["taxon_id"])
        self.assertEqual({"G1"}, self.members([self.gene()], species="  Medicago  truncatula  ", records=records)[0])
        for species in ("Arabidopsis thaliana", "Medicago sativa", "Medicago tribuloides", "Medicago truncatula subsp. fictitious"):
            with self.subTest(species=species):
                self.assertEqual(set(), self.members([self.gene()], species=species, records=records)[0])
        for change in ({"reference_species": "Arabidopsis thaliana"},
                       {"reference_species": "Medicago truncatula (unrecorded strain)"},
                       {"reference_taxon_id": "3702"}, {"reference_taxon_id": ""}):
            self.assertEqual(set(), self.members([self.gene(**change)], records=records)[0])

    def test_endpoint_copy_and_namespace_axes_are_not_combined(self):
        genes = [self.gene(), self.gene(gene_id="G2", reference_taxon_id="3702")]
        cards = [{"interval_gene_id": "G1", "reference_gene": "G7KFE4", "reference_species": DISPLAY, "taxon_id": "3880"}]
        self.assertEqual({"G1"}, self.members(genes, cards=cards, records={"3880": TAXONOMY})[0])
        # The correct taxon for OTHER must not lend its identity to G7KFE4.
        genes = [self.gene(reference_taxon_id="3702")]
        cards = [{"interval_gene_id": "G1", "reference_gene": "OTHER", "reference_species": DISPLAY, "taxon_id": "3880"}]
        self.assertEqual(set(), self.members(genes, cards=cards, records={"3880": TAXONOMY})[0])
        typed = "uniprot|G7KFE4|ENTRY"
        genes = [self.gene(reference_gene=typed, reference_accession_namespace="WRONG")]
        cards = [{"interval_gene_id": "G1", "reference_gene": typed, "reference_species": "Arabidopsis thaliana",
                  "reference_accession_namespace": "uniprot", "taxon_id": "3702"}]
        self.assertEqual(set(), self.members(genes, cards=cards, records={"3880": TAXONOMY}, reference=typed)[0])

    def test_reserved_taxonomy_files_are_separate_from_gene_cache_and_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "taxonomy-3880.json"
            source.write_text(json.dumps(TAXONOMY))
            (root / "G7KFE4.json").write_text(json.dumps({"primaryAccession": "G7KFE4", "organism": {"scientificName": SCIENTIFIC}}))
            records, artifacts = load_taxonomy_cache(str(root))
            self.assertEqual({"3880"}, set(records))
            self.assertEqual({"G7KFE4"}, set(_load_record_cache(str(root))))
            self.assertEqual(str(source), artifacts[0]["path"])
            source.write_text(json.dumps({**TAXONOMY, "taxonId": 3702}))
            with self.assertRaisesRegex(ValueError, "mismatched reference taxonomy"):
                load_taxonomy_cache(str(root))
            source.write_text(json.dumps(TAXONOMY))
            (root / "taxonomy-3702.json").write_text(json.dumps(TAXONOMY))
            with self.assertRaisesRegex(ValueError, "mismatched reference taxonomy"):
                load_taxonomy_cache(str(root))

    def test_projection_uses_only_own_endpoint_taxon_and_preserves_raw_evidence(self):
        records = {"3880": TAXONOMY}
        gene = self.gene(orthology_reference_gene="G7KFE4", orthology_reference_species=DISPLAY,
                         orthology_relationship="one_to_many")
        claim = {"reference_gene": "G7KFE4", "reference_species": SCIENTIFIC}
        unit = {"unit_id": "synthetic", "genes": [gene]}
        original = copy.deepcopy(unit)
        self.assertEqual("", route_b_homolog_projection(claim, unit, "G1", {"G1"}))
        projected = json.loads(route_b_homolog_projection(claim, unit, "G1", {"G1"}, taxonomy_records=records))
        self.assertEqual("frozen_gene_orthology", projected["basis"])
        self.assertEqual(DISPLAY, projected["evidence"][0]["orthology_reference_species"])
        self.assertEqual(original, unit)
        hit = {"hit_accession": "G7KFE4", "hit_species": DISPLAY, "hit_taxon_id": "3880",
               "percent_identity": "90", "query_coverage": "0.9"}
        hit_unit = {"unit_id": "synthetic", "genes": [{"gene_id": "G1", "sequence_similarity_hits": [hit]}]}
        self.assertEqual("qualified_sequence_hit", json.loads(route_b_homolog_projection(
            claim, hit_unit, "G1", {"G1"}, taxonomy_records=records))["basis"])
        hit["hit_taxon_id"] = "3702"
        hit_unit["genes"][0]["reference_taxon_id"] = "3880"
        self.assertEqual("", route_b_homolog_projection(claim, hit_unit, "G1", {"G1"}, taxonomy_records=records))
        hit.update(hit_taxon_id="3880", query_coverage="0.05")
        self.assertEqual("", route_b_homolog_projection(claim, hit_unit, "G1", {"G1"}, taxonomy_records=records))
        typed = "uniprot|G7KFE4|ENTRY"
        wrong_namespace = self.gene(reference_gene=typed, orthology_reference_gene=typed,
            orthology_reference_species=DISPLAY, reference_accession_namespace="WRONG")
        self.assertEqual("", route_b_homolog_projection({**claim, "reference_gene": typed},
            {"genes": [wrong_namespace]}, "G1", {"G1"}, taxonomy_records=records))
        conflicting_hit = {"hit_accession": typed, "hit_species": DISPLAY, "hit_taxon_id": "3880",
                           "reference_accession_namespace": "WRONG", "percent_identity": "90", "query_coverage": "0.9"}
        with_card = {"genes": [{"gene_id": "G1", "sequence_similarity_hits": [conflicting_hit]}],
            "reference_material_cards": [{"interval_gene_id": "G1", "reference_gene": typed, "reference_species": DISPLAY,
                                           "taxon_id": "3880", "reference_accession_namespace": "uniprot"}]}
        self.assertEqual("", route_b_homolog_projection({**claim, "reference_gene": typed}, with_card,
                                                      "G1", {"G1"}, taxonomy_records=records))

    def test_bare_and_decorated_endpoints_keep_same_t1_through_final_wet(self):
        from tests.test_first_stage_homolog_t1 import report
        from candidate_qtl.relevance_tiers import NEW_TIERS
        for endpoint in (SCIENTIFIC, DISPLAY):
            with self.subTest(endpoint=endpoint), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                cache = root / "records"
                cache.mkdir()
                (cache / "taxonomy-3880.json").write_text(json.dumps(TAXONOMY))
                _, _, frozen, wet = RouteBInitialReferenceIdentityTests()._chain(root, [SCIENTIFIC],
                    endpoint_species=endpoint, endpoint_taxon_id="3880", record_cache_arg=str(cache),
                    paper_reports=[report()])
                for artifact in (frozen["review_candidates"], wet["wet_candidates"]):
                    frame = pd.read_csv(artifact, sep="\t", dtype=str).fillna("")
                    self.assertEqual({"G1"}, set(frame.candidate_gene_id))
                    self.assertEqual({NEW_TIERS[0]}, set(frame.relevance_tier))
                    self.assertIn("frozen_gene_orthology", frame.iloc[0]["relevance_evidence_json"])
                from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
                (cache / "taxonomy-3880.json").write_text(json.dumps({**TAXONOMY, "scientificName": "Wrong plant"}))
                with self.assertRaisesRegex(ValueError, "taxonomy cache hash mismatch"):
                    wet_candidate_handoff(str(root / "tampered-wet"),
                        route_b_candidates_arg=frozen["candidates"], route_b_candidates_manifest_arg=frozen["manifest"],
                        review_candidates_arg=frozen["review_candidates"], review_candidates_manifest_arg=frozen["review_manifest"],
                        interval_atlas_arg=str(root / "atlas.tsv"), case_id="initial-reference-fixture")

    def test_projection_namespace_missing_values_stay_unknown_with_taxonomy(self):
        reference = "uniprot|G7KFE4|ENTRY"
        claim = {"reference_gene": reference, "reference_species": SCIENTIFIC}
        for species in (SCIENTIFIC, DISPLAY):
            for namespace in ("", "unknown", "not_evaluated", "uniprot", "WRONG"):
                with self.subTest(species=species, namespace=namespace):
                    hit = {"hit_accession": reference, "hit_species": species, "hit_taxon_id": "3880",
                           "reference_accession_namespace": namespace, "percent_identity": "90", "query_coverage": "0.9"}
                    unit = {"genes": [{"gene_id": "G1", "sequence_similarity_hits": [hit]}]}
                    before = copy.deepcopy(unit)
                    projected = route_b_homolog_projection(claim, unit, "G1", {"G1"}, taxonomy_records={"3880": TAXONOMY})
                    self.assertEqual(namespace != "WRONG", bool(projected))
                    if species == SCIENTIFIC:
                        self.assertEqual(bool(projected), bool(route_b_homolog_projection(claim, unit, "G1", {"G1"})))
                    if projected:
                        self.assertEqual(namespace, json.loads(projected)["evidence"][0]["reference_accession_namespace"])
                    self.assertEqual(before, unit)

    def test_canonical_initial_session_uses_sealed_cache_and_keeps_support_unverified(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cache = root / "records"
            cache.mkdir()
            (cache / "taxonomy-3880.json").write_text(json.dumps(TAXONOMY))
            fixture = RouteBInitialReferenceIdentityTests()
            _, transfer, frozen, _ = fixture._chain(root, [SCIENTIFIC], endpoint_species=DISPLAY,
                endpoint_taxon_id="3880", record_cache_arg=str(cache))
            self.assertEqual(0, frozen["candidate_gene_count"])
            self.assertEqual(1, frozen["review_candidate_gene_count"])
            manifest = json.loads(Path(transfer["manifest"]).read_text())
            self.assertTrue(manifest["reference_species_display_audit"]["resolutions"])
            source_text = root / "source-text.json"
            source_text.write_text(json.dumps({"artifact_role": "source_text_for_review_only", "records": {}}))
            inputs = {"hypotheses": manifest["ancestors"]["proposals"]["path"],
                      "proposal": str(root / "responses.json"), "records": str(cache), "sources": str(source_text),
                      "bundle": str(root / "bundle/route_b_semantic_bundle.json"),
                      "source_verification": manifest["ancestors"]["source_verification"]["path"]}
            original = {path: path.read_bytes() for path in (Path(inputs["hypotheses"]), Path(inputs["bundle"]),
                        Path(inputs["proposal"]), cache / "taxonomy-3880.json")}
            session = root / "session"
            prepared = prepare_session(str(session), route="b", inputs=inputs, proposer_id="synthetic-actor",
                proposer_model="synthetic", case_id="initial-reference-fixture", review_policy="initial_bibliography_v1")
            result = complete_session(str(session), prepared["request_sha256"])
            self.assertEqual(1, result["candidates"]["review_candidate_gene_count"])
            self.assertEqual(0, result["candidates"]["candidate_gene_count"])
            self.assertEqual(0, result["verification"]["verified_claim_count"])
            review = pd.read_csv(result["candidates"]["review_candidates"], sep="\t", dtype=str).fillna("")
            self.assertIn("frozen_gene_orthology", review.iloc[0]["relevance_evidence_json"])
            for path, before in original.items():
                self.assertEqual(before, path.read_bytes())
            # An independently prepared session rejects changes to its frozen taxonomy bytes.
            second = root / "tampered-session"
            prepared = prepare_session(str(second), route="b", inputs=inputs, proposer_id="synthetic-actor",
                proposer_model="synthetic", case_id="initial-reference-fixture", review_policy="initial_bibliography_v1")
            (second / "inputs/records/taxonomy-3880.json").write_text(json.dumps({**TAXONOMY, "scientificName": "Wrong plant"}))
            with self.assertRaisesRegex(ValueError, "frozen review session artifact mismatch"):
                complete_session(str(second), prepared["request_sha256"])
            # The candidate consumer also verifies the taxonomy artifacts bound by transfer.
            (cache / "taxonomy-3880.json").write_text(json.dumps({**TAXONOMY, "scientificName": "Wrong plant"}))
            with self.assertRaisesRegex(ValueError, "taxonomy cache hash mismatch"):
                freeze_route_b_semantic_candidates(transfer["claims"], str(root / "tampered-freeze"),
                                                   case_id="initial-reference-fixture")


if __name__ == "__main__":
    unittest.main()
