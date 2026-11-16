"""Neutral provider namespaces and atomic-name parsing regressions."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from candidate_qtl import llm_pathway_legacy as pathway
from candidate_qtl.identifiers import is_usable_symbol
from candidate_qtl.reference_gene_resolution import reference_gene_resolve
from candidate_qtl.reference_name_pairs import atomic_reference_pairs, exact_record_candidates


def _record(accession: str, name: str = "NODE1", locus: str = "LOC_12345",
            species: str = "Species alpha") -> dict:
    return {"primaryAccession": accession, "organism": {"scientificName": species},
            "genes": [{"geneName": {"value": name}, "synonyms": [{"value": "ALT1"}],
                       "orderedLocusNames": [{"value": locus}]}],
            "sequence": {"value": "M" * 30},
            "references": [{"citation": {"citationCrossReferences": [{"database": "PubMed", "id": "123"}]}}]}


def _response(records: list[dict]):
    response = mock.MagicMock()
    response.read.return_value = json.dumps({"results": records}).encode("utf-8")
    response.headers = {}
    response.__enter__.return_value = response
    return response


class ReferenceNamespaceTests(unittest.TestCase):
    def test_middle_alternative_annotation_is_not_an_endpoint(self):
        self.assertEqual(atomic_reference_pairs("GENA / GENB (annotation) / GENC", "Species alpha"),
                         ([("GENA", "Species alpha"), ("GENB", "Species alpha"),
                           ("GENC", "Species alpha")], "atomic_name_species"))

    def test_annotations_do_not_create_aliases_or_species_crossproducts(self):
        self.assertEqual(atomic_reference_pairs("NODE1 (ALT9/ALT8) | ALT1 (prose)", "Species alpha")[0],
                         [("NODE1", "Species alpha"), ("ALT1", "Species alpha")])
        self.assertEqual(atomic_reference_pairs("NODE1 (prose); NODE2 (more)", "Species alpha; Species beta")[0],
                         [("NODE1", "Species alpha"), ("NODE2", "Species beta")])
        self.assertEqual(atomic_reference_pairs("NODE1/ALT1 (prose)", "Species alpha; Species beta"),
                         ([], "ambiguous_name_species_pairing"))

    def test_malformed_or_multiple_annotations_fail_closed(self):
        for name in ("NODE1 / ALT1 (open / NODE2", "NODE1 / ALT1) / NODE2",
                     "NODE1 / ALT1 ((nested)) / NODE2", "NODE1 / ALT1 (one) (two) / NODE2",
                     "NODE1 / ALT1 (prose) tail / NODE2", "NODE1 // NODE2"):
            with self.subTest(name=name):
                self.assertEqual(atomic_reference_pairs(name, "Species alpha"),
                                 ([], "ambiguous_name_species_pairing"))

    def test_all_provider_declared_names_and_loci_match_exactly(self):
        records = {"Q00001": _record("Q00001"), "Q00002": _record("Q00002")}
        for name in ("node1", "alt1", "loc_12345"):
            with self.subTest(name=name):
                self.assertEqual(exact_record_candidates(records, name, "species  alpha"), ["Q00001", "Q00002"])

    def test_partial_or_extended_provider_names_do_not_match(self):
        records = {"Q00001": _record("Q00001", "NODE1.1", "LOC_123450")}
        for name in ("NODE1", "LOC_12345", "ALT"):
            with self.subTest(name=name):
                self.assertEqual(exact_record_candidates(records, name, "Species alpha"), [])

    def test_species_and_primary_accession_guards_remain_exact(self):
        self.assertEqual(exact_record_candidates({"Q00001": _record("Q00001")}, "LOC_12345", "Species beta"), [])
        self.assertEqual(exact_record_candidates({"Q00002": _record("Q00001")}, "LOC_12345", "Species alpha"), [])

    def test_malformed_locus_metadata_is_not_a_name_source(self):
        record = _record("Q00001")
        record["genes"][0]["orderedLocusNames"] = {"value": "LOC_12345"}
        self.assertEqual(exact_record_candidates({"Q00001": record}, "LOC_12345", "Species alpha"), [])

    def test_live_full_agi_endpoint_matches_ordered_loci_and_preserves_alternatives(self):
        records = [_record("Q00002", locus="At1g12345", species="Arabidopsis thaliana"),
                   _record("Q00001", locus="AT1G12345", species="Arabidopsis thaliana"),
                   _record("Q00003", locus="AT1G12345", species="Arabidopsis lyrata"),
                   _record("Q00004", locus="AT1G123450", species="Arabidopsis thaliana")]
        with mock.patch.object(pathway.urllib.request, "urlopen", return_value=_response(records)) as provider:
            result = pathway._resolve_accession_candidates("at1g12345", "Arabidopsis thaliana")
        self.assertEqual(result, (["Q00001", "Q00002"], "ambiguous_gene_name",
                                  'gene:AT1G12345 AND organism_name:"Arabidopsis thaliana"', 2))
        provider.assert_called_once()
        self.assertFalse(is_usable_symbol("AT1G12345"))
        self.assertEqual(pathway._searchable_symbol("AT1G12345"), "")

    def test_live_symbol_does_not_accept_a_provider_prefix(self):
        records = [_record("Q00001", "NODE1 OTHER"), _record("Q00002", "NODE1.1")]
        with mock.patch.object(pathway.urllib.request, "urlopen", return_value=_response(records)):
            self.assertEqual(pathway._resolve_accession_candidates("NODE1", "Species alpha")[0], [])

    def test_live_exact_symbol_synonym_and_locus_use_the_same_record_filter(self):
        records = [_record("Q00001", locus="At1g12345", species="Arabidopsis thaliana")]
        for name in ("NODE1", "ALT1", "AT1G12345"):
            with self.subTest(name=name), mock.patch.object(pathway.urllib.request, "urlopen", return_value=_response(records)):
                self.assertEqual(pathway._resolve_accession_candidates(name, "Arabidopsis thaliana")[0], ["Q00001"])

    def test_canonical_live_resolver_preserves_each_exact_locus_record(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pd.DataFrame([{"case_id": "neutral", "hypothesis_id": "h001", "reference_gene_name": "AT1G12345",
                           "reference_species": "Arabidopsis thaliana", "source_ids": "", "proposed_stable_id": ""}]).to_csv(
                               root / "in.tsv", sep="\t", index=False)
            records = [_record(accession, locus="At1g12345", species="Arabidopsis thaliana")
                       for accession in ("Q00001", "Q00002")]
            (root / "records.json").write_text(json.dumps({"results": records}), encoding="utf-8")
            with mock.patch.object(pathway.urllib.request, "urlopen", return_value=_response(records)) as provider:
                result = reference_gene_resolve(str(root / "in.tsv"), str(root / "out"),
                                                record_cache_arg=str(root / "records.json"), allow_external_api=True)
            provider.assert_called_once()
            resolved = pd.read_csv(result["mechanisms_resolved"], sep="\t", dtype=str).fillna("")
            self.assertEqual(set(resolved["proposed_stable_id"]), {"Q00001", "Q00002"})
            self.assertTrue(all(resolved["parent_hypothesis_id"] == "h001"))
            self.assertTrue(all(resolved["source_ids"] == "PMID:123"))

    def test_canonical_offline_resolver_preserves_each_exact_locus_record(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pd.DataFrame([{"case_id": "neutral", "hypothesis_id": "h001", "reference_gene_name": "AT1G12345",
                           "reference_species": "Arabidopsis thaliana", "source_ids": "", "proposed_stable_id": ""}]).to_csv(
                               root / "in.tsv", sep="\t", index=False)
            records = [_record(accession, locus="At1g12345", species="Arabidopsis thaliana")
                       for accession in ("Q00001", "Q00002")]
            (root / "records.json").write_text(json.dumps({"results": records}), encoding="utf-8")
            with mock.patch.object(pathway.urllib.request, "urlopen", side_effect=AssertionError("offline")):
                result = reference_gene_resolve(str(root / "in.tsv"), str(root / "out"),
                                                record_cache_arg=str(root / "records.json"), allow_external_api=False)
            resolved = pd.read_csv(result["mechanisms_resolved"], sep="\t", dtype=str).fillna("")
            self.assertEqual(set(resolved["proposed_stable_id"]), {"Q00001", "Q00002"})
            self.assertTrue(all(resolved["parent_hypothesis_id"] == "h001"))
            self.assertTrue(all(resolved["source_ids"] == "PMID:123"))
