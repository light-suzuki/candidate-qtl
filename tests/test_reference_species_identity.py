"""Retrieved reference species must agree before a sequence can be projected."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl.bibliography import INITIAL_POLICY
from candidate_qtl.llm_pathway import bait_panel_freeze, hypothesis_verify
from candidate_qtl.phenotype_investigation import mechanism_import
from candidate_qtl.reference_gene_resolution import reference_gene_resolve
from candidate_qtl.reference_name_pairs import exact_record_candidates, record_species_matches
from candidate_qtl.source_verify import source_verify


def _record(scientific_name="Arabidopsis thaliana"):
    return {
        "primaryAccession": "P12345",
        "organism": {"scientificName": scientific_name, "commonName": "thale cress"},
        "genes": [{"geneName": {"value": "REF1"}}],
        "sequence": {"value": "MAAGSVDPLA" * 12},
        "references": [{"citation": {"citationCrossReferences": [
            {"database": "PubMed", "id": "12345"},
        ]}}],
    }


def _hypothesis(hypothesis_id, species, name=""):
    return {
        "hypothesis_id": hypothesis_id,
        "pathway_or_process": "meristem growth",
        "gene_family": "reference family",
        "reference_gene_name": name,
        "reference_species": species,
        "proposed_stable_id": "P12345",
        "mechanism_statement": "Reduced meristem growth reduces branch number",
        "expected_role": "meristem proliferation",
        "source_ids": "PMID:12345",
    }


class ReferenceSpeciesIdentityTests(unittest.TestCase):
    def _run(self, root, rows, record=None):
        cache = root / "cache"
        cache.mkdir()
        (cache / "record.json").write_text(json.dumps(record or _record()), encoding="utf-8")
        payload = root / "payload.json"
        payload.write_text(json.dumps({"hypotheses": rows}), encoding="utf-8")
        imported = mechanism_import(str(payload), str(root / "import"), case_id="species")
        resolved = reference_gene_resolve(
            imported["mechanisms"], str(root / "resolve"), record_cache_arg=str(cache),
            focal_species="Pisum sativum",
        )
        publications = root / "publications.json"
        publications.write_text(json.dumps({"id": "12345", "pmid": "12345", "title": "Synthetic meristem study"}),
                                encoding="utf-8")
        sources = source_verify(
            resolved["mechanisms_resolved"], str(root / "sources"), cache_arg=str(publications),
            verification_policy=INITIAL_POLICY,
        )
        hypotheses = sources["mechanisms_with_verified_sources"]
        with patch("candidate_qtl.uniprot_error_policy.urllib.request.urlopen",
                   side_effect=AssertionError("offline fixture must not fetch")):
            verified = hypothesis_verify(hypotheses, str(root / "verified"), "Pisum sativum",
                                         record_cache_arg=str(cache))
        frozen = bait_panel_freeze(
            hypotheses, verified["validation"], verified["verified_sequences"], str(root / "frozen"),
            reviewable_sequences_arg=verified["reviewable_sequences"], case_id="species",
        )
        validation = pd.read_csv(verified["validation"], sep="\t", dtype=str).fillna("").set_index("hypothesis_id")
        projection = pd.read_csv(frozen["projection_panel"], sep="\t", dtype=str).fillna("")
        return validation, projection, resolved, verified, sources

    def test_normal_accession_path_rejects_same_genus_other_species_and_keeps_the_healthy_cached_row(self):
        for healthy_first in (True, False):
            with self.subTest(healthy_first=healthy_first), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                rows = [_hypothesis("healthy", "Arabidopsis thaliana"), _hypothesis("wrong-species", "Arabidopsis lyrata")]
                validation, projection, resolved, verified, sources = self._run(root, rows if healthy_first else rows[::-1])
                ledger = pd.read_csv(resolved["ledger"], sep="\t", dtype=str)
                self.assertEqual({"model_supplied"}, set(ledger.resolution_origin))
                bibliography = pd.read_csv(sources["mechanisms_with_verified_sources"], sep="\t", dtype=str)
                self.assertEqual({"resolved"}, set(bibliography.bibliography_state))
                self.assertEqual({INITIAL_POLICY}, set(bibliography.candidate_admission_policy))
                self.assertEqual("True", validation.loc["healthy", "species_match"])
                self.assertEqual("False", validation.loc["wrong-species", "species_match"])
                self.assertEqual("not_verified", validation.loc["wrong-species", "identity_sequence_verification_status"])
                self.assertEqual("declared_species_does_not_match_accession", validation.loc["wrong-species", "rejection_reason"])
                self.assertEqual({"healthy"}, set(projection.hypothesis_id))
                self.assertEqual({"review_only"}, set(projection.projection_eligibility))
                self.assertEqual({"healthy"}, set(json.loads(Path(verified["reviewable_sequences"]).read_text())))

    def test_name_case_whitespace_and_later_explicit_species_alternatives_still_project(self):
        scenarios = [
            ("  arabidopsis   THALIANA  ", "  ref1  "),
            ("Citrullus lanatus / Arabidopsis thaliana", ""),
            ("Citrullus lanatus; arabidopsis THALIANA | Lotus japonicus", ""),
        ]
        for species, name in scenarios:
            with self.subTest(species=species, name=name), tempfile.TemporaryDirectory() as tmp:
                validation, projection, _, _, _ = self._run(Path(tmp), [_hypothesis("matching", species, name)])
                self.assertEqual("True", validation.loc["matching", "species_match"])
                self.assertEqual("verified", validation.loc["matching", "identity_sequence_verification_status"])
                self.assertEqual({"matching"}, set(projection.hypothesis_id))

    def test_unconfirmed_names_and_subspecies_labels_do_not_gain_a_new_species_rule(self):
        record = _record()
        for species in ("Arabidopsis lyrata", "Arabidopsis", "thale cress", "Arabidopsis thaliana subsp. test"):
            with self.subTest(species=species), tempfile.TemporaryDirectory() as tmp:
                self.assertFalse(record_species_matches(record, species))
                self.assertEqual([], exact_record_candidates({"P12345": record}, "REF1", species))
                validation, projection, _, _, _ = self._run(Path(tmp), [_hypothesis("unconfirmed", species)])
                self.assertEqual("False", validation.loc["unconfirmed", "species_match"])
                self.assertEqual("not_verified", validation.loc["unconfirmed", "identity_sequence_verification_status"])
                self.assertTrue(projection.empty)

    def test_resolver_and_sequence_gate_share_the_existing_scientific_name_normalization(self):
        record = _record("Arabidopsis   thaliana")
        self.assertTrue(record_species_matches(record, "  ARABIDOPSIS thaliana  "))
        self.assertEqual(["P12345"], exact_record_candidates({"P12345": record}, "ref1", "  ARABIDOPSIS thaliana  "))
        with tempfile.TemporaryDirectory() as tmp:
            validation, projection, _, _, _ = self._run(
                Path(tmp), [_hypothesis("named", "  ARABIDOPSIS thaliana  ", "ref1")], record,
            )
            self.assertEqual("True", validation.loc["named", "species_match"])
            self.assertEqual({"named"}, set(projection.hypothesis_id))


if __name__ == "__main__":
    unittest.main()
