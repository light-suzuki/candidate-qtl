"""Graded reference identity: infraspecific records, species prefixes, provisional
declared accessions and the independent identity review ledger."""

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
from candidate_qtl.reference_name_pairs import (IDENTITY_DECLARED_ACCESSION, IDENTITY_EXACT, IDENTITY_MISMATCH,
                                                IDENTITY_REVIEWED_ALIAS, declared_record_identity,
                                                exact_record_candidates, load_identity_reviews,
                                                record_species_matches, species_prefixed_name)
from candidate_qtl.source_verify import source_verify


def _record(accession="P12345", name="CHLG", organism="Oryza sativa subsp. japonica"):
    return {
        "primaryAccession": accession,
        "organism": {"scientificName": organism},
        "genes": [{"geneName": {"value": name}}],
        "sequence": {"value": "MAAGSVDPLA" * 12},
        "references": [{"citation": {"citationCrossReferences": [{"database": "PubMed", "id": "12345"}]}}],
    }


def _hypothesis(name, species="Oryza sativa", accession="P12345"):
    return {
        "hypothesis_id": "h1", "pathway_or_process": "chlorophyll synthesis",
        "gene_family": "chlorophyll synthase", "reference_gene_name": name,
        "reference_species": species, "proposed_stable_id": f"UniProt:{accession}",
        "mechanism_statement": "Loss of chlorophyll synthesis gives yellow pods",
        "expected_role": "chlorophyll a synthesis", "source_ids": "PMID:12345",
    }


def _review_ledger(path, **overrides):
    row = {"reference_gene_name": "YGL1", "reference_species": "Oryza sativa", "accession": "UniProt:P12345",
           "record_gene_name": "CHLG", "verdict": "same_gene", "evidence_source_ids": "PMID:12345",
           "reviewer": "identity-reviewer-1", "rationale": "ygl1 is a chlorophyll synthase mutant"}
    row.update(overrides)
    pd.DataFrame([row]).to_csv(path, sep="\t", index=False)
    return path


class SpeciesAndPrefixRuleTests(unittest.TestCase):
    def test_infraspecific_record_matches_declared_binomial_but_not_the_reverse(self):
        self.assertTrue(record_species_matches(_record(), "Oryza sativa"))
        self.assertTrue(record_species_matches(_record(organism="Oryza sativa var. test"), "oryza  SATIVA"))
        self.assertFalse(record_species_matches(_record(organism="Oryza sativa"), "Oryza sativa subsp. japonica"))
        self.assertFalse(record_species_matches(_record(organism="Oryza glaberrima"), "Oryza sativa"))
        self.assertFalse(record_species_matches(_record(organism="Oryza sativaX subsp. a"), "Oryza sativa"))
        self.assertFalse(record_species_matches(_record(organism="Oryza sativa hybrid"), "Oryza sativa"))
        self.assertFalse(record_species_matches(_record(), "Oryza"))

    def test_species_prefix_is_removed_only_for_the_declared_species(self):
        self.assertEqual("GLK2", species_prefixed_name("SlGLK2", "Solanum lycopersicum"))
        self.assertEqual("", species_prefixed_name("SlGLK2", "Arabidopsis thaliana"))
        self.assertEqual("", species_prefixed_name("Slow1", "Solanum lycopersicum"))
        record = _record(name="GLK2", organism="Solanum lycopersicum")
        self.assertEqual(["P12345"], exact_record_candidates({"P12345": record}, "SlGLK2", "Solanum lycopersicum"))
        self.assertEqual([], exact_record_candidates({"P12345": record}, "AtGLK2", "Solanum lycopersicum"))


class DeclaredRecordIdentityTests(unittest.TestCase):
    def test_grades(self):
        record = _record()
        self.assertEqual(IDENTITY_EXACT, declared_record_identity(record, "P12345", ["CHLG"], "Oryza sativa"))
        self.assertEqual(IDENTITY_DECLARED_ACCESSION,
                         declared_record_identity(record, "P12345", ["YGL1"], "Oryza sativa", {"P12345"}))
        # An accession found only by search is not a declaration.
        self.assertEqual(IDENTITY_MISMATCH, declared_record_identity(record, "P12345", ["YGL1"], "Oryza sativa"))
        self.assertEqual(IDENTITY_MISMATCH,
                         declared_record_identity(record, "P12345", ["YGL1"], "Zea mays", {"P12345"}))
        self.assertEqual(IDENTITY_MISMATCH,
                         declared_record_identity(record, "Q99999", ["YGL1"], "Oryza sativa", {"Q99999"}))

    def test_review_confirms_only_with_a_gene_name_on_the_record(self):
        with tempfile.TemporaryDirectory() as tmp:
            reviews = load_identity_reviews(_review_ledger(Path(tmp) / "ok.tsv"))
            self.assertEqual(IDENTITY_REVIEWED_ALIAS, declared_record_identity(
                _record(), "P12345", ["YGL1"], "Oryza sativa", {"P12345"}, reviews))
            wrong_name = load_identity_reviews(_review_ledger(Path(tmp) / "bad.tsv", record_gene_name="OTHER"))
            self.assertEqual(IDENTITY_DECLARED_ACCESSION, declared_record_identity(
                _record(), "P12345", ["YGL1"], "Oryza sativa", {"P12345"}, wrong_name))
            for overrides in ({"verdict": "different_gene"}, {"evidence_source_ids": ""}, {"reviewer": ""}):
                with self.subTest(overrides=overrides):
                    self.assertEqual({}, load_identity_reviews(_review_ledger(Path(tmp) / "x.tsv", **overrides)))


class ProvisionalAdoptionPipelineTests(unittest.TestCase):
    def _run(self, root, hypothesis, record, identity_review=None, bound=False):
        cache = root / "cache"
        cache.mkdir()
        (cache / "record.json").write_text(json.dumps(record), encoding="utf-8")
        payload = root / "payload.json"
        payload.write_text(json.dumps({"hypotheses": [hypothesis]}), encoding="utf-8")
        imported = mechanism_import(str(payload), str(root / "import"), case_id="grading")
        resolved = reference_gene_resolve(imported["mechanisms"], str(root / "resolve"),
                                          record_cache_arg=str(cache), focal_species="Pisum sativum",
                                          identity_review_arg=str(identity_review) if identity_review else None)
        publications = root / "publications.json"
        publications.write_text(json.dumps({"id": "12345", "pmid": "12345", "title": "Synthetic study"}),
                                encoding="utf-8")
        sources = source_verify(resolved["mechanisms_resolved"], str(root / "sources"),
                                cache_arg=str(publications), verification_policy=INITIAL_POLICY)
        hypotheses = sources["mechanisms_with_verified_sources"]
        claims = None
        if bound:
            from tests.bound_claim_fixtures import bound_claim
            claims = root / "claims.tsv"
            pd.DataFrame([bound_claim(route="a", claim_id="h1", citation="PMID:12345",
                subject_gene="P12345", organism="Oryza sativa", phenotype="yellow pods")]).to_csv(
                    claims, sep="\t", index=False)
        with patch("candidate_qtl.uniprot_error_policy.urllib.request.urlopen",
                   side_effect=AssertionError("offline fixture must not fetch")):
            verified = hypothesis_verify(hypotheses, str(root / "verified"), "Pisum sativum",
                                         record_cache_arg=str(cache), bound_claims_arg=str(claims) if claims else None)
        frozen = bait_panel_freeze(hypotheses, verified["validation"], verified["verified_sequences"],
                                   str(root / "frozen"), reviewable_sequences_arg=verified["reviewable_sequences"],
                                   case_id="grading")
        validation = pd.read_csv(verified["validation"], sep="\t", dtype=str).fillna("").set_index("hypothesis_id")
        projection = pd.read_csv(frozen["projection_panel"], sep="\t", dtype=str).fillna("")
        ledger = pd.read_csv(resolved["ledger"], sep="\t", dtype=str).fillna("")
        requests = pd.read_csv(resolved["identity_review_requests"], sep="\t", dtype=str).fillna("")
        return validation, projection, ledger, requests

    def test_subspecies_record_with_matching_name_is_exact(self):
        with tempfile.TemporaryDirectory() as tmp:
            validation, projection, ledger, requests = self._run(
                Path(tmp), _hypothesis("NYC1"), _record(name="NYC1"))
            self.assertEqual({IDENTITY_EXACT}, set(ledger.reference_identity_status))
            self.assertTrue(requests.empty)
            self.assertEqual({"h1"}, set(projection.hypothesis_id))

    def test_name_mismatch_on_declared_accession_is_provisional_review_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            validation, projection, ledger, requests = self._run(
                Path(tmp), _hypothesis("YGL1"), _record(name="CHLG"))
            self.assertEqual({IDENTITY_DECLARED_ACCESSION}, set(ledger.reference_identity_status))
            self.assertEqual(["P12345"], list(requests.accession))
            self.assertEqual("CHLG", requests.loc[0, "record_gene_names"])
            self.assertNotEqual("verified", validation.loc["h1", "verification_status"])
            self.assertEqual("verified", validation.loc["h1", "identity_sequence_verification_status"])
            self.assertEqual({"h1"}, set(projection.hypothesis_id))
            self.assertEqual({"review_only"}, set(projection.projection_eligibility))

    def test_bound_claim_preserves_provisional_review_projection(self):
        with tempfile.TemporaryDirectory() as tmp:
            validation, projection, _, _ = self._run(
                Path(tmp), _hypothesis("YGL1"), _record(name="CHLG"), bound=True)
            self.assertEqual("True", validation.loc["h1", "bound_claim_authorized"])
            self.assertEqual({"h1"}, set(projection.hypothesis_id))
            self.assertEqual({"review_only"}, set(projection.projection_eligibility))

    def test_identity_review_confirms_the_provisional_record(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ledger_path = _review_ledger(root / "identity_review.tsv")
            validation, projection, ledger, requests = self._run(
                root, _hypothesis("YGL1"), _record(name="CHLG"), identity_review=ledger_path)
            self.assertEqual({IDENTITY_REVIEWED_ALIAS}, set(ledger.reference_identity_status))
            self.assertTrue(requests.empty)
            self.assertEqual({"h1"}, set(projection.hypothesis_id))

    def test_declared_accession_from_another_species_stays_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            validation, projection, ledger, requests = self._run(
                Path(tmp), _hypothesis("YGL1"), _record(name="CHLG", organism="Zea mays"))
            self.assertNotIn(IDENTITY_DECLARED_ACCESSION, set(ledger.reference_identity_status))
            self.assertTrue(requests.empty)
            self.assertTrue(projection.empty)


if __name__ == "__main__":
    unittest.main()
