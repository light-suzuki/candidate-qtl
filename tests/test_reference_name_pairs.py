from __future__ import annotations

import unittest

from candidate_qtl.reference_name_pairs import (atomic_reference_pairs, exact_record_candidates,
                                                reference_lookup_handles, reference_lookup_names,
                                                reference_lookup_names_for_inputs)


class AtomicReferencePairsTests(unittest.TestCase):
    def test_exact_orf_names_keep_species_primary_identity_and_all_alternatives(self):
        records = {key: {"primaryAccession": key, "organism": {"scientificName": "Species alpha"},
                         "genes": [{"orfNames": [{"value": "LOC_001"}]}]}
                   for key in ("Q00002", "Q00001")}
        self.assertEqual(["Q00001", "Q00002"], exact_record_candidates(records, " loc_001 ", "Species alpha"))
        self.assertEqual([], exact_record_candidates(records, "LOC_001", "Species beta"))
        self.assertEqual([], exact_record_candidates({"Q00003": records["Q00001"]}, "LOC_001", "Species alpha"))
        self.assertEqual([], exact_record_candidates(records, "LOC001", "Species alpha"))
        self.assertEqual([], exact_record_candidates(records, "SpLOC_001", "Species alpha"))
        for malformed in (None, "LOC_001", {}, [None, "LOC_001", {}]):
            record = {**records["Q00001"], "genes": [{"orfNames": malformed}]}
            with self.subTest(malformed=malformed):
                self.assertEqual([], exact_record_candidates({"Q00001": record}, "LOC_001", "Species alpha"))

    def test_explicit_parallel_clauses_and_synonyms(self):
        self.assertEqual(atomic_reference_pairs("NODE1 (plant one); NODE2/ALT2 (plant two)",
                                               "Species alpha; Species beta"),
                         ([("NODE1", "Species alpha"), ("NODE2", "Species beta"),
                           ("ALT2", "Species beta")], "atomic_name_species"))

    def test_one_species_can_be_shared_by_explicit_name_alternatives(self):
        self.assertEqual(atomic_reference_pairs("NODE1/NODE2; NODE3", "Species alpha")[0],
                         [("NODE1", "Species alpha"), ("NODE2", "Species alpha"), ("NODE3", "Species alpha")])

    def test_unbound_species_and_empty_clauses_do_not_create_crossproduct(self):
        for name, species in [("NODE1/NODE2", "Species alpha; Species beta"),
                              ("NODE1;", "Species alpha"), ("NODE1", "Species alpha/Species beta")]:
            with self.subTest(name=name):
                self.assertEqual(atomic_reference_pairs(name, species), ([], "ambiguous_name_species_pairing"))
        self.assertEqual(atomic_reference_pairs("/".join(f"NODE{i}" for i in range(17)), "Species alpha"),
                         ([], "atomic_reference_budget_exceeded"))

    def test_authoritative_exact_names_species_and_primary_identity(self):
        record = {"primaryAccession": "Q00001", "organism": {"scientificName": "Species alpha"},
                  "genes": [{"geneName": {"value": "NODE1"}, "synonyms": [{"value": "ALT1"}]}]}
        self.assertEqual(exact_record_candidates({"Q00001": record}, "alt1", "species  alpha"), ["Q00001"])
        self.assertEqual(exact_record_candidates({"Q00001": record}, "NODE1 OTHER", "Species alpha"), [])
        self.assertEqual(exact_record_candidates({"Q00001": record}, "NODE1", "Species beta"), [])
        self.assertEqual(exact_record_candidates({"Q00002": record}, "NODE1", "Species alpha"), [])

    def test_formal_names_and_typed_lookup_handles_are_preserved_without_alias_guessing(self):
        name = "TFL1 (TERMINAL FLOWER 1)"
        self.assertEqual(reference_lookup_names(name, "TFL1", "Arabidopsis thaliana"),
                         ["TFL1", "TERMINAL FLOWER 1"])
        self.assertEqual(reference_lookup_handles([{
            "proposed_stable_id": "UniProt:P93003",
            "source_ids": "TAIR:AT5G03840; PMID:12345678",
        }]), [("uniprot", "P93003"), ("tair", "AT5G03840")])
        self.assertEqual(reference_lookup_handles([{"proposed_stable_id": "fake__TFL1",
                                                    "source_ids": "TAIR:AT5G0384"}]), [])
        self.assertEqual(reference_lookup_handles([{"proposed_stable_id": "", "source_ids": "P93003"}]), [])

    def test_composite_formal_name_queries_stay_with_the_paired_species(self):
        name = "KN1 (KNOTTED1); FT (FLOWERING LOCUS T)"
        self.assertEqual(reference_lookup_names(name, "KN1", "Zea mays"),
                         ["KN1", "KNOTTED1"])
        self.assertEqual(reference_lookup_names(name, "FT", "Arabidopsis thaliana"),
                         ["FT", "FLOWERING LOCUS T"])

    def test_short_explicit_parenthetical_symbol_is_exact_same_branch_identity(self):
        raw_name = "LONG NODE (LN1); OTHER NODE (ON2)"
        raw_species = "Species neutralis; Species secunda"
        self.assertEqual(reference_lookup_names_for_inputs(
            [{"reference_gene_name": raw_name, "reference_species": raw_species}],
            "LONG NODE", "Species neutralis"), ["LONG NODE", "LN1"])
        self.assertEqual(reference_lookup_names_for_inputs(
            [{"reference_gene_name": raw_name, "reference_species": raw_species}],
            "OTHER NODE", "Species secunda"), ["OTHER NODE", "ON2"])
        self.assertEqual(reference_lookup_names_for_inputs(
            [{"reference_gene_name": "NODE ONE (N1); NODE TWO (N2)",
              "reference_species": "Species one; Species two; Species three"}],
            "NODE ONE", "Species one"), ["NODE ONE"])
        record = {"primaryAccession": "Q00001", "organism": {"scientificName": "Species neutralis"},
                  "genes": [{"geneName": {"value": "LN1"}}]}
        self.assertEqual(exact_record_candidates({"Q00001": record}, ["LONG NODE", "LN1"],
                                                 "Species neutralis"), ["Q00001"])
        self.assertEqual(exact_record_candidates({"Q00001": record}, ["OTHER NODE", "ON2"],
                                                 "Species neutralis"), [])

    def test_parenthetical_symbols_normalize_whitespace_without_crossing_branches_or_species(self):
        for whitespace in ("  ", "\t", " \n "):
            raw_name = f"LONG{whitespace}NODE (LN1)/OTHER NODE (ON2); LONG NODE (LN3)"
            raw_species = "Species neutralis; Species secunda"
            inputs = [{"reference_gene_name": raw_name, "reference_species": raw_species}]
            pairs, status = atomic_reference_pairs(raw_name, raw_species)
            self.assertEqual("atomic_name_species", status)
            with self.subTest(whitespace=repr(whitespace)):
                self.assertEqual(reference_lookup_names_for_inputs(inputs, *pairs[0]), ["LONG NODE", "LN1"])
                self.assertEqual(reference_lookup_names_for_inputs(inputs, *pairs[1]), ["OTHER NODE", "ON2"])
                self.assertEqual(reference_lookup_names_for_inputs(inputs, *pairs[2]), ["LONG NODE", "LN3"])
                self.assertEqual(reference_lookup_names_for_inputs(inputs, "LONG NODE", "Species unknown"), ["LONG NODE"])
