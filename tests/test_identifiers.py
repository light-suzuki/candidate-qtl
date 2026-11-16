"""Unit tests for the shared identifier classifier (issue #297)."""

import unittest

from candidate_qtl import identifiers


class UniprotAccessionTests(unittest.TestCase):
    def test_6_char_accessions_are_detected(self):
        for token in ("F8T7T1", "Q9LZS3", "P0A8L9", "G7K1D0"):
            with self.subTest(token=token):
                self.assertTrue(identifiers.is_uniprot_accession(token))

    def test_10_char_accessions_are_detected(self):
        for token in ("A0ACM8R380", "A0A0C3WAY6", "A0A072UD83", "A0A1E5HLJ4"):
            with self.subTest(token=token):
                self.assertTrue(identifiers.is_uniprot_accession(token))

    def test_real_symbols_are_not_accessions(self):
        for token in ("FT", "FTa2", "CCaMK", "GA3ox", "D14", "SERKL2", "NORK",
                      "NFP", "SBE1", "CHI", "LOX", "GAI", "CIK1", "PsFTa2"):
            with self.subTest(token=token):
                self.assertFalse(identifiers.is_uniprot_accession(token))

    def test_agi_is_not_an_accession(self):
        self.assertFalse(identifiers.is_uniprot_accession("AT1G01010"))


class AgiLocusTests(unittest.TestCase):
    def test_agi_loci_are_detected(self):
        for token in ("AT1G01010", "AT5G38430", "ATMG00060", "AT4G12345.1"):
            with self.subTest(token=token):
                self.assertTrue(identifiers.is_agi_locus(token))

    def test_symbols_are_not_agi_loci(self):
        for token in ("FT", "GA3ox", "Nup85"):
            self.assertFalse(identifiers.is_agi_locus(token))


class LocusTagTests(unittest.TestCase):
    def test_legume_locus_tags_are_detected(self):
        for token in ("MTR_7g085020", "MTR.5g034210", "Medtr7g085020",
                      "LotjaGi1g1v0638100.1", "Lj1g3v0123450.1",
                      "GLYMA_01G000100", "LOC_Os01g01010", "Os01g01010",
                      "Sb01g01010", "Zm00001d027292"):
            with self.subTest(token=token):
                self.assertTrue(identifiers.is_locus_tag(token))

    def test_letter_first_locus_tag_families_are_detected(self):
        # IWGSC wheat (TraesCS...), Brassica (Brara / Bna ...D) and other
        # letter-first conventions registered in the prefix table must stay
        # detectable as locus tags (regression guard).
        for token in ("TraesCS5A02G401300", "TraesCS2A01G000100",
                      "BraraA01001", "BraraC01234",
                      "BnaA01g12340D", "BnaC02g23450D",
                      "BvuBv9_000010", "TuG000010"):
            with self.subTest(token=token):
                self.assertTrue(identifiers.is_locus_tag(token))
                self.assertFalse(identifiers.is_usable_symbol(token))

    def test_jgi_lj_observed_identifier_shapes_are_detected(self):
        # Issue #344: the pinned legume proteome (frozen DB legume_plant_pep3)
        # carries identifier-only deflines in two shapes beyond the
        # digit-leading rule, both observed in real data (45,662 rows /
        # 4,574 unique):
        # - an "_LC.<version>" low-confidence gene-model suffix, e.g.
        #   LotjaGi2g1v0260700_LC.1 (4,520 unique ids)
        # - letter-leading assembly-variant ids LotjaGiC1v... / LotjaGiM1v...
        #   ("C" = chromosome, "M" = mitochondrion), 54 unique ids
        for token in ("LotjaGi2g1v0260700_LC.1", "LotjaGi0g1v0001000_LC.2",
                      "LotjaGiC1v0000600.1", "LotjaGiM1v0000600.1"):
            with self.subTest(token=token):
                self.assertTrue(identifiers.is_locus_tag(token))
                self.assertFalse(identifiers.is_usable_symbol(token))

    def test_unobserved_identifier_shapes_stay_rejected(self):
        # Issue #344: only the shapes observed in real data are accepted; the
        # grammar must not be widened into an unbounded identifier regex.
        # Unknown _LC suffixes, unobserved assembly-variant markers and
        # suffix-less _LC tails all stay rejected.
        for token in ("LotjaGi2g1v0260700_LC.1x", "LotjaGiX9v0000000.1",
                      "LotjaGi2g1v0260700_LC", "LotjaGiC1v0000600.1x",
                      "LotjaGi0g1v0001000_LCx.2"):
            with self.subTest(token=token):
                self.assertFalse(identifiers.is_locus_tag(token))

    def test_symbols_are_not_locus_tags(self):
        for token in ("SBE1", "SBEI", "FTa2", "FT", "CCaMK", "GA3ox", "GAI",
                      "CHI", "LOX", "LUX", "WOX", "NFP", "LYK3", "NSP2",
                      "CIK1", "PsFTa2", "Nup85"):
            with self.subTest(token=token):
                self.assertFalse(identifiers.is_locus_tag(token))

    def test_accessions_are_not_locus_tags(self):
        for token in ("F8T7T1", "A0ACM8R380", "A0A0C3WAY6"):
            self.assertFalse(identifiers.is_locus_tag(token))

    def test_species_prefixed_real_symbols_are_not_locus_tags(self):
        # A real symbol with a species prefix must never be misread as a locus
        # tag just because the prefix table contains the species code.
        for token in ("OsSWEET1", "LjNFP2", "TaGW2", "SbFTL1", "ZmCCT9"):
            with self.subTest(token=token):
                self.assertFalse(identifiers.is_locus_tag(token))
                self.assertTrue(identifiers.is_usable_symbol(token))


class ClassifyTests(unittest.TestCase):
    def test_classification(self):
        cases = {
            "FT": "gene_symbol",
            "FTa2": "gene_symbol",
            "A0ACM8R380": "uniprot_accession",
            "F8T7T1": "uniprot_accession",
            "MTR_7g085020": "locus_tag",
            "LotjaGi1g1v0638100.1": "locus_tag",
            "AT1G01010": "agi_locus",
            "": "empty",
            None: "empty",
        }
        for value, expected in cases.items():
            with self.subTest(value=value):
                self.assertEqual(identifiers.classify_gene_name_field(value), expected)


class UsableSymbolTests(unittest.TestCase):
    def test_usable_symbols(self):
        for token in ("FT", "FTa2", "SBE1", "CCaMK", "GA3ox", "Flowering locus protein T"):
            with self.subTest(token=token):
                self.assertTrue(identifiers.is_usable_symbol(token))

    def test_identifiers_are_not_usable(self):
        for token in ("F8T7T1", "A0ACM8R380", "A0A0C3WAY6", "MTR_7g085020",
                      "LotjaGi1g1v0638100.1", "AT1G01010", "", None):
            with self.subTest(token=token):
                self.assertFalse(identifiers.is_usable_symbol(token))

    def test_digit_leading_fragments_are_not_usable(self):
        # Locus-tag fragments without their prefix (7g085020) must not become
        # family symbols.
        for token in ("7g085020", "01g01010", "1g1v0638100.1"):
            with self.subTest(token=token):
                self.assertFalse(identifiers.is_usable_symbol(token))


class DeflineFieldTests(unittest.TestCase):
    DEFLINE_FT = (
        "tr|F8T7T1|F8T7T1_MEDTR Flowering locus protein T "
        "OS=Medicago truncatula GN=MTR_7g085020 PE=3 SV=1"
    )
    DEFLINE_FTA2 = (
        "tr|A0ACM8R380|A0ACM8R380_MEDTR FTa2 OS=Medicago truncatula PE=2 SV=1"
    )

    def test_gene_name_field(self):
        self.assertEqual(identifiers.gene_name_field(self.DEFLINE_FT), "MTR_7g085020")
        self.assertEqual(identifiers.gene_name_field(self.DEFLINE_FTA2), None)
        self.assertEqual(identifiers.gene_name_field("no field here"), None)

    def test_gene_name_field_rejects_attribute_wrappers(self):
        # Old-style deflines use GN=Name=FT / GN=OrderedLocusNames=...; the
        # wrapper is not a symbol and must not become a high-confidence GN.
        self.assertEqual(
            identifiers.gene_name_field("... GN=Name=FT OS=x"), None)
        self.assertEqual(
            identifiers.gene_name_field("... GN=OrderedLocusNames=MTR_1g019150 OS=x"),
            None)

    def test_uniprot_name_field(self):
        self.assertEqual(
            identifiers.uniprot_name_field("tr|F8T7T1|F8T7T1_MEDTR", self.DEFLINE_FT),
            "Flowering locus protein T",
        )
        self.assertEqual(
            identifiers.uniprot_name_field("tr|A0ACM8R380|A0ACM8R380_MEDTR",
                                           self.DEFLINE_FTA2),
            "FTa2",
        )

    def test_uniprot_name_field_requires_uniprot_subject(self):
        self.assertEqual(
            identifiers.uniprot_name_field("Medicago:MtrunA17r5:NFP",
                                           "NFP-like protein OS=x"),
            None,
        )


class SymbolCandidatesTests(unittest.TestCase):
    DEFLINE_FT = (
        "tr|F8T7T1|F8T7T1_MEDTR Flowering locus protein T "
        "OS=Medicago truncatula GN=MTR_7g085020 PE=3 SV=1"
    )
    DEFLINE_FTA2 = (
        "tr|A0ACM8R380|A0ACM8R380_MEDTR FTa2 OS=Medicago truncatula PE=2 SV=1"
    )

    def test_usable_gn_wins(self):
        defline = "sp|Q9LZS3|FT_ARATH Flowering locus T OS=Arabidopsis GN=FT PE=1"
        candidates = identifiers.symbol_candidates("sp|Q9LZS3|FT_ARATH", defline)
        self.assertEqual(candidates[0]["value"], "FT")
        self.assertEqual(candidates[0]["source"], "gn")

    def test_curated_alias_resolves_protein_name(self):
        aliases = {
            "Flowering locus protein T": "ft",
            "FTa2": "ft",
        }
        candidates = identifiers.symbol_candidates(
            "tr|F8T7T1|F8T7T1_MEDTR", self.DEFLINE_FT, alias_lookup=aliases.get)
        self.assertEqual(candidates[0]["value"], "ft")
        self.assertEqual(candidates[0]["source"], "curated_alias")

    def test_symbol_like_protein_name_fallback(self):
        candidates = identifiers.symbol_candidates(
            "tr|A0ACM8R380|A0ACM8R380_MEDTR", self.DEFLINE_FTA2)
        self.assertEqual(candidates[0]["value"], "FTa2")
        self.assertEqual(candidates[0]["source"], "protein_name")

    def test_multiword_name_without_alias_yields_no_candidate(self):
        candidates = identifiers.symbol_candidates(
            "tr|F8T7T1|F8T7T1_MEDTR", self.DEFLINE_FT)
        self.assertEqual(candidates, [])

    def test_gn_locus_tag_is_skipped(self):
        candidates = identifiers.symbol_candidates(
            "tr|F8T7T1|F8T7T1_MEDTR", self.DEFLINE_FT, alias_lookup=None)
        self.assertFalse(any(c["source"] == "gn" for c in candidates))


class GrammarConsistencyTests(unittest.TestCase):
    def test_uniprot_grammar_matches_hardening(self):
        # The citable-source grammar in hardening must not drift away from the
        # shared classifier.
        from candidate_qtl.hardening import _UNIPROT
        self.assertEqual(identifiers.UNIPROT_ACCESSION_PATTERN, _UNIPROT)

    def test_grammar_valid_6_char_non_opq_accession(self):
        # F8T7T1 is a real Medicago entry and must be an accession; the shared
        # grammar covers 6-char accessions beyond the historical [OPQ] prefix.
        self.assertTrue(identifiers.is_uniprot_accession("F8T7T1"))
        self.assertTrue(identifiers.is_uniprot_accession("A0A1B2"))


if __name__ == "__main__":
    unittest.main()
