"""Shared Route A/B/C family-label normalization (issue #204)."""

from __future__ import annotations

import unittest

from candidate_qtl.family_labels import (
    ALIAS_PROVENANCE,
    alias_canonical,
    families_match,
    families_match_expanded,
    family_identity_keys,
    family_identity_keys_expanded,
    family_match_level,
    family_match_level_expanded,
    family_set_contains_expanded,
    family_set_match_level,
    family_set_match_level_expanded,
    split_respecting_brackets,
)


class FamilyLabelsTests(unittest.TestCase):
    def test_benchmark_alias_pairs_resolve_many_to_one(self) -> None:
        pairs = (
            ("Phytochrome (PHY)", "PHYA"),
            ("MLO (Mildew resistance Locus O)", "MLO"),
            ("WD40 repeat protein (TTG1-clade)", "TTG"),
            ("RWP-RK transcription factor (NIN-like)", "NLP"),
            ("LEAFY/FLORICAULA (LFY-clade)", "LFY"),
            ("PEBP (FT-clade)", "FT"),
            ("Eukaryotic translation initiation factor 4E (eIF4E-clade)", "IF4E4"),
        )
        for descriptive, symbol in pairs:
            with self.subTest(descriptive=descriptive, symbol=symbol):
                self.assertTrue(families_match(descriptive, symbol))

    def test_generic_words_do_not_form_identity_keys(self) -> None:
        self.assertFalse(families_match("receptor kinase protein", "transcription factor protein"))
        self.assertNotIn("symbol:protein", family_identity_keys("receptor kinase protein"))

    def test_descriptive_route_a_labels_meet_symbol_namespace_issue248(self) -> None:
        """Route A descriptive labels from the supervised 58-case run (issue #248).

        Each pair is (Route A `gene_family` label, truth family symbol) and
        corresponds to a case that strict string normalization scored as a
        miss even though the label names the truth family.  The expanded
        matcher (uppercase symbol tokens inside descriptive labels) must
        connect the two namespaces; strict `families_match` is deliberately
        left alone for masking and negative-control uses.
        Source artifacts:
        ~/pea_qtl_runs/supervised_58case_20260810/work/*/triple/a06_route_a_status/
        """
        pairs = (
            ("CIK LRR receptor kinases", "CIK1"),  # cik_CIK1
            ("CIK LRR receptor kinases", "CIK4"),  # cik_CIK4
            ("CIK (CLAVATA INSENSITIVE KINASE) LRR-RLK co-receptors", "CIK5"),  # cik_CIK5
            ("CIK (CLAVATA INSENSITIVE KINASE) LRR-RLK co-receptors", "CIK6"),  # cik_CIK6
            ("CLV1-like LRR-RLK / CLE ligand (CLV1, CLV3)", "CLE13"),  # cle_CLE13
            ("CLV1-like LRR-RLK / CLE ligand (CLV1, CLV3)", "CLE25"),  # cle_CLE25
            ("CLE/CLV1-like LRR-RLK (MtCLE12/MtCLE13 - SUNN/HAR1/NARK)", "CLE12"),  # cle_MtCLE12
            ("CLE/CLV1-like LRR-RLK (MtCLE12/MtCLE13 - SUNN/HAR1/NARK)", "CLE13"),  # cle_MtCLE13
            ("Leucine-rich repeat receptor-like kinases, CLV1/BAM family (BAM1/BAM2 focus)", "BAM1"),  # clv_BAM12
            ("Leucine-rich repeat receptor-like kinases, CLV1/BAM family (BAM3 focus)", "BAM3"),  # clv_BAM3
            ("CLV1-like LRR-RLK receptors of AON (HAR1/NARK/SUNN)", "HAR1"),  # clv_CLV1
            ("LRR receptor-like protein CLV2 and CORYNE (CRN) pseudokinase", "CLV2"),  # clv_CLV2
            ("CLE signaling peptide family (CLV3/CLE-class ligands)", "CLV3"),  # clv_CLV3
            ("LRR receptor-like protein CLV2 and CORYNE (CRN) pseudokinase", "CRN"),  # clv_CRN
            ("Kinase-associated protein phosphatase KAPP (PP2C-type)", "KAPP"),  # clv_KAPP
            ("RPK2 LRR receptor kinase (CLV3 signal transmitter)", "RPK2"),  # clv_RPK2
            ("WOX homeodomain transcription factors (WOX5)", "WOX5"),  # clv_WOX5
            ("WOX homeodomain transcription factors (WUSCHEL)", "WUS"),  # clv_WUS
            ("MYB-bHLH-WD40 (MBW) complex (TTG1, TT2, TT8)", "TTG1"),  # pea_A2
            ("PEBP family (TFL1/CEN/SP)", "TFL1"),  # pea_DET
            ("FT/TFL1 PEBP family", "FT"),  # pea_GIGAS
            ("FT/TFL1 PEBP family", "TFL1"),  # pea_LF
            ("BR biosynthesis enzymes (DET2/5alpha-reductase)", "DET2"),  # pea_LK
            ("phytochrome A photoreceptor", "PHYA"),  # pea_PHYA
            ("AP1/FUL MADS-box family", "AP1"),  # pea_PIM
            ("TCP transcription factor (BRC1/BRC2) family", "TCP18"),  # pea_PsBRC1
            ("CCD8 carotenoid cleavage dioxygenase", "CCD8"),  # pea_RMS1
            ("SL receptor D14 alpha/beta hydrolase and F-box protein MAX2", "D14"),  # pea_RMS3
            ("F-box protein (MAX2/D3) family", "MAX2"),  # pea_RMS4
            ("LEAFY/FLO family", "LFY"),  # pea_UNI
            ("MLO seven-transmembrane proteins", "MLO2"),  # pea_er1
            ("eIF4E family (cap-binding proteins)", "EIF4E1"),  # pea_sbm1
            ("LysM-RLK NFP and LYK3 (Medicago)", "NFP"),  # sym_SYM10
            ("LysM-RLK NFP and LYK3 (Medicago)", "LYK3"),  # sym_SYM37
            ("NSP2 GRAS transcription factor", "NSP2"),  # sym_SYM7
            ("CASTOR/POLLUX ion channels", "POLLUX"),  # sym_SYM8
        )
        for descriptive, symbol in pairs:
            with self.subTest(descriptive=descriptive, symbol=symbol):
                self.assertTrue(families_match_expanded(descriptive, symbol))

    def test_expanded_never_loosens_masking_or_decoy_separation(self) -> None:
        """Expanded keys exist only in the dedicated Route A matcher, never in
        the strict family_identity_keys used by masking and negative checks."""
        self.assertFalse(families_match("CIK LRR receptor kinases", "CIK1"))
        self.assertTrue(families_match_expanded("CIK LRR receptor kinases", "CIK1"))
        self.assertNotIn("symbol:cik", family_identity_keys("CIK LRR receptor kinases"))


class CamelCaseTokenTests(unittest.TestCase):
    """Issue #301: camelCase symbol tokens bridge Route A labels (expanded only)."""

    def test_camelcase_symbol_tokens_bridge_route_a_labels(self) -> None:
        """Lowercase-mixed symbol tokens must meet truth symbols in the expanded
        matcher while the strict matcher (masking / decoys) stays untouched."""
        pairs = (
            ("CCaMK kinase family", "CCAMK"),
            ("GA3ox hydroxylase family", "GA3OX1"),
            ("Nup85 complex", "NUP85"),
            ("Nup133 subcomplex", "NUP133"),
            ("FTa2 protein", "FT"),
        )
        for descriptive, symbol in pairs:
            with self.subTest(descriptive=descriptive, symbol=symbol):
                self.assertTrue(families_match_expanded(descriptive, symbol))
                self.assertFalse(families_match(descriptive, symbol))

    def test_generic_camel_words_never_become_symbol_keys(self) -> None:
        """Words that only look camelCase (uppercase count < 2, no digit) are
        structurally excluded: Protein, Kinase, Receptor, Factor, Family,
        Complex, Domain, Enzyme, Locus, Le, Nup, Pea."""
        for word in (
            "Protein", "Kinase", "Receptor", "Factor", "Family",
            "Complex", "Domain", "Enzyme", "Locus", "Le", "Nup", "Pea",
        ):
            with self.subTest(word=word):
                keys = family_identity_keys_expanded(f"{word} family")
                self.assertNotIn(f"symbol:{word.casefold()}", keys)
                self.assertFalse(families_match_expanded(f"{word} family", word.upper()))

    def test_no_leading_run_subkey_is_generated(self) -> None:
        """The tokenizer never emits a bare leading-run key, so "Nup85" does
        not leak "nup" and "FTa2" does not leak "ft" via the token layer."""
        self.assertNotIn("symbol:nup", family_identity_keys_expanded("Nup85 complex"))
        self.assertNotIn("symbol:ft", family_identity_keys_expanded("FTa2 protein"))


class RomanNumeralTokenTests(unittest.TestCase):
    """Issue #301: trailing Roman I/II/III converts to Arabic (expanded only)."""

    def test_roman_suffix_converts_to_arabic(self) -> None:
        pairs = (
            ("SBEI", "SBE1"),
            ("PsSBEI", "SBE1"),
            ("SBEII", "SBE2"),
            ("GBSSI", "GBSS1"),
        )
        for left, right in pairs:
            with self.subTest(left=left, right=right):
                self.assertTrue(families_match_expanded(left, right))
                self.assertFalse(families_match(left, right))

    def test_known_symbols_ending_in_i_are_never_converted(self) -> None:
        """Three-letter symbols ending in I, X-ending names and prefix-stripped
        short stems must keep their keys (spec guard list)."""
        for token in (
            "CHI", "LOX", "LUX", "WOX", "UNI", "KNOX", "POLLUX", "GA20OX",
            "LYKX", "GAI", "CLV", "PsGI", "nodX", "BBI", "CRI", "DOI", "NFP",
        ):
            with self.subTest(token=token):
                keys = family_identity_keys_expanded(token)
                self.assertNotIn(f"symbol:{token.casefold()[:-1]}1", keys)

    def test_non_trailing_roman_suffix_is_not_converted(self) -> None:
        self.assertNotIn("symbol:sbe2", family_identity_keys_expanded("SBEIIa"))
        self.assertFalse(families_match_expanded("SBEIIa", "SBE2"))

    def test_roman_iii_converts_to_3(self) -> None:
        self.assertTrue(families_match_expanded("SBEIII", "SBE3"))

    def test_real_words_with_trailing_i_letters_are_excluded(self) -> None:
        # COPI / TAXI / KARI / PSII / COII / NODI end in I as a letter, not a
        # numeral: converting them would fabricate false synonyms with real
        # genes (COPI -> COP1 pairs with the photomorphogenesis E3 ligase).
        for token in ("COPI", "TAXI", "KARI", "PSII", "COII", "NODI"):
            with self.subTest(token=token):
                keys = family_identity_keys_expanded(token)
                self.assertNotIn(f"symbol:{token.casefold()[:-1]}1", keys)
                self.assertNotIn(f"symbol:{token.casefold()[:-2]}2", keys)
        self.assertFalse(families_match_expanded("COPI", "COP1"))
        self.assertFalse(families_match_expanded("COII", "CO2"))

    def test_sbeii_never_matches_sbe1(self) -> None:
        # The SBEI -> sbe1 alias is strict; SBEII is a different enzyme and must
        # not ride the sbe1 alias nor the expanded sbe2 conversion against sbe1.
        self.assertFalse(families_match("SBEII", "SBE1"))
        self.assertFalse(families_match_expanded("SBEII", "SBE1"))
        self.assertTrue(families_match_expanded("SBEII", "SBE2"))

    def test_legumin_never_matches_leguminosae(self) -> None:
        self.assertFalse(families_match("Leguminosae family", "LEGUMIN"))
        self.assertFalse(families_match_expanded("Leguminosae family", "LEGUMIN"))


class CuratedAliasExpansionTests(unittest.TestCase):
    """Issues #297/#301: curated aliases added at the family boundary."""

    def test_storage_protein_family_aliases(self) -> None:
        pairs = (
            ("vicilin-like 7S globulin", "VICILIN"),
            ("Legumin (11S globulin)", "LEGUMIN"),
            ("legume lectin", "LECTIN"),
            ("2S seed albumin", "ALBUMIN"),
            ("convicilin 7S", "CONVICILIN"),
            ("Flowering locus protein T", "FT"),
        )
        for descriptive, symbol in pairs:
            with self.subTest(descriptive=descriptive, symbol=symbol):
                self.assertTrue(families_match(descriptive, symbol))

    def test_convicilin_stays_a_distinct_family_from_vicilin(self) -> None:
        self.assertFalse(families_match("convicilin", "VICILIN"))
        self.assertFalse(families_match_expanded("convicilin", "VICILIN"))

    def test_nsp2_descriptive_name_meets_the_symbol(self) -> None:
        self.assertTrue(families_match("Nodulation-signaling pathway 2", "NSP2"))
        self.assertTrue(families_match("NSP2 GRAS transcription factor", "NSP2"))

    def test_ft_paralog_names_meet_ft_only_in_the_expanded_matcher(self) -> None:
        self.assertTrue(families_match_expanded("FTa2", "FT"))
        self.assertFalse(families_match("FTa2", "FT"))


class AliasCanonicalTests(unittest.TestCase):
    """Issues #297/#301: alias_canonical resolves curated patterns only."""

    def test_alias_canonical_resolves_curated_patterns(self) -> None:
        cases = (
            ("Flowering locus protein T", "ft"),
            ("FTa2", "ft"),
            ("Nodulation-signaling pathway 2", "nsp2"),
            ("SBEI", "sbe1"),
        )
        for value, canonical in cases:
            with self.subTest(value=value):
                self.assertEqual(canonical, alias_canonical(value))

    def test_alias_canonical_returns_none_for_unknown_protein_names(self) -> None:
        for value in ("CCaMK kinase family", "Nup85 complex", "unknown protein", ""):
            with self.subTest(value=value):
                self.assertIsNone(alias_canonical(value))


class Brc1AliasTests(unittest.TestCase):
    """Issue #419: BRC1 / TB1 / TCP18 are the same shoot-branching TCP factor.

    Arabidopsis BRANCHED1 has locus name TCP18 (UniProtKB:A1YKT1, GN=TCP18,
    Synonyms=BRC1/TBL1), maize TEOSINTE BRANCHED 1 is the ortholog
    (UniProtKB:Q93WI2, GN=TB1), and pea PsBRC1 is the pea ortholog. The
    memorization probe answers with the common name (PsBRC1/TB1) while the
    catalog stores the Arabidopsis locus name (TCP18); string comparison alone
    therefore misreports a recalled case as a reasoning success.
    """

    def test_brc1_tb1_tcp18_resolve_to_one_family(self) -> None:
        for left, right in (
            ("BRC1", "TCP18"),
            ("TB1", "TCP18"),
            ("PsBRC1", "TCP18"),
            ("BRANCHED 1", "TCP18"),
            ("TEOSINTE BRANCHED 1", "TCP18"),
            ("TB1", "BRC1"),
        ):
            with self.subTest(left=left, right=right):
                self.assertTrue(families_match(left, right))

    def test_alias_canonical_and_provenance(self) -> None:
        self.assertEqual("brc1", alias_canonical("PsBRC1"))
        self.assertEqual("brc1", alias_canonical("TB1"))
        self.assertEqual("brc1", alias_canonical("TCP18"))
        self.assertEqual("brc1", alias_canonical("BRANCHED 1"))
        provenance = ALIAS_PROVENANCE["brc1"]
        self.assertIn("A1YKT1", provenance["accession"])
        self.assertIn("Q93WI2", provenance["accession"])

    def test_generic_tcp_family_does_not_ride_the_brc1_alias(self) -> None:
        # "TCP transcription factor" alone names the family, not the BRANCHED1
        # locus; only the specific spellings (TCP18 / BRC1 / TB1) are aliases.
        self.assertFalse(families_match("TCP transcription factor", "TCP18"))
        self.assertIsNone(alias_canonical("TCP transcription factor"))


class FamilyMatchLevelTests(unittest.TestCase):
    """Issue #268: report *how* a family match was made, not just whether."""

    def test_stem_collision_pairs_report_stem_not_exact(self) -> None:
        """The five stem collisions documented for supervised_58case_20260810."""
        pairs = (
            ("CIK1", "CIK6"), ("CIK1", "CIK4"), ("CIK4", "CIK5"),
            ("CLE12", "CLE13"), ("CLE13", "CLE25"),
            ("CLV2", "CLV3"),
            ("BAM1", "BAM3"),
            ("ELF3", "ELF4"),
        )
        for left, right in pairs:
            with self.subTest(left=left, right=right):
                self.assertEqual("stem", family_match_level(left, right))
                self.assertTrue(families_match(left, right))

    def test_exact_symbol_pairs_report_exact(self) -> None:
        self.assertEqual("exact_symbol", family_match_level("CIK1", "CIK1"))
        self.assertEqual("exact_symbol", family_match_level("Phytochrome (PHY)", "PHYA"))

    def test_no_match_reports_empty_string(self) -> None:
        self.assertEqual("", family_match_level("CIK1", "unrelated gene"))
        self.assertEqual("", family_match_level("", "CIK1"))

    def test_expanded_route_a_label_reports_level(self) -> None:
        # A bare "CIK" token in a descriptive label carries no paralog number,
        # so it can only ever establish a stem-level match against a
        # numbered truth symbol -- even though families_match_expanded
        # reports it as a hit (issue #248 behaviour, unchanged).
        self.assertEqual("stem", family_match_level_expanded("CIK LRR receptor kinases", "CIK1"))
        self.assertEqual("stem", family_match_level_expanded("CIK LRR receptor kinases", "CIK6"))
        self.assertEqual("", family_match_level_expanded("CIK LRR receptor kinases", "unrelated"))

    def test_family_set_match_level_prefers_exact_over_stem(self) -> None:
        # CIK6 is a stem-only match; CIK1 is exact. The best level across the
        # set must be reported, not just the first value's level.
        self.assertEqual("exact_symbol", family_set_match_level(["CIK6", "CIK1"], "CIK1"))
        self.assertEqual("stem", family_set_match_level(["CIK6"], "CIK1"))
        self.assertEqual("", family_set_match_level([], "CIK1"))
        self.assertEqual(
            "exact_symbol",
            family_set_match_level_expanded(["CIK LRR receptor kinases", "CIK1 co-receptor"], "CIK1"),
        )

    def test_never_changes_existing_match_functions(self) -> None:
        """families_match / family_set_contains(_expanded) keep every current hit."""
        self.assertTrue(families_match("CIK1", "CIK6"))
        self.assertTrue(family_set_contains_expanded(["CLE25 peptides"], "CLE13"))
        self.assertTrue(family_set_contains_expanded(["BAM1 receptor kinases"], "BAM3"))
        self.assertTrue(family_set_contains_expanded(["CIK1 co-receptors"], "CIK6"))
        self.assertFalse(families_match("CIK LRR receptor kinases", "CIK1"))

    def test_split_respecting_brackets(self) -> None:
        """Bracket-enclosed semicolons must not split free-text hypotheses."""
        text = "MBW complex (TT2/TT8/TTG1-type in Arabidopsis; AN2/AN1/AN11-type in Petunia); MYB activator"
        parts = split_respecting_brackets(text, ";")
        self.assertEqual([
            "MBW complex (TT2/TT8/TTG1-type in Arabidopsis; AN2/AN1/AN11-type in Petunia)",
            "MYB activator",
        ], parts)

        # Quotes and nested brackets
        quoted = 'famA [note: a; b]; "famB; with; quote"; famC (nested (x; y); z)'
        self.assertEqual([
            "famA [note: a; b]",
            '"famB; with; quote"',
            "famC (nested (x; y); z)",
        ], split_respecting_brackets(quoted, ";"))

        # Empty and simple string
        self.assertEqual([], split_respecting_brackets("", ";"))
        self.assertEqual(["single_family"], split_respecting_brackets("single_family", ";"))


if __name__ == "__main__":
    unittest.main()
