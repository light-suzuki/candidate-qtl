"""Genome-wide family index and interval hypothesis falsification."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl import family_index as family_index_module
from candidate_qtl.cli import main


def _hit(query: str, subject: str, identity: float, bitscore: float, *, qlen: int = 300) -> str:
    # outfmt 6: qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore
    return f"{query}\t{subject}\t{identity}\t{qlen}\t{qlen}\t{qlen}\t1\t{qlen}\t1\t{qlen}\t1e-60\t{bitscore}"


def _hit13(query: str, subject: str, identity: float, bitscore: float, title: str, *, qlen: int = 300) -> str:
    # outfmt 6 with a stitle (13 columns): ... evalue bitscore stitle
    return (f"{query}\t{subject}\t{identity}\t{qlen}\t{qlen}\t{qlen}\t1\t{qlen}\t1\t{qlen}\t1e-60\t{bitscore}\t{title}")


class HypergeometricTests(unittest.TestCase):
    def test_survival_matches_hand_computed_values(self) -> None:
        # 2 successes out of 4 draws from 10 items containing 5 successes.
        # P(X>=2) = 1 - P(0) - P(1) = 1 - 5/210 - 50/210 = 155/210
        self.assertAlmostEqual(155 / 210, family_index_module.hypergeom_sf(2, 10, 5, 4), places=10)
        self.assertEqual(1.0, family_index_module.hypergeom_sf(0, 10, 5, 4))
        self.assertEqual(0.0, family_index_module.hypergeom_sf(5, 10, 4, 4))
        # A singleton family is degenerate: p equals the interval fraction.
        self.assertAlmostEqual(400 / 44000, family_index_module.hypergeom_sf(1, 44000, 1, 400), places=9)

    def test_large_universe_does_not_overflow(self) -> None:
        value = family_index_module.hypergeom_sf(3, 44000, 6, 400)
        self.assertTrue(0.0 < value < 1e-3, value)

    def test_benjamini_hochberg_is_monotone(self) -> None:
        adjusted = family_index_module.benjamini_hochberg([0.001, 0.01, 0.04, 0.5])
        self.assertEqual(adjusted, sorted(adjusted))
        self.assertTrue(all(0.0 <= value <= 1.0 for value in adjusted))
        self.assertGreaterEqual(adjusted[0], 0.001)


class SymbolTests(unittest.TestCase):
    def test_symbol_and_stem_extraction(self) -> None:
        self.assertEqual("CLV1", family_index_module.subject_symbol("sp|P0CLV1|CLV1_ARATH"))
        self.assertEqual("AT1G75820", family_index_module.subject_symbol("AT1G75820"))
        self.assertEqual("CLV", family_index_module.symbol_stem("CLV1"))
        self.assertEqual("PIN", family_index_module.symbol_stem("PIN8"))
        # Locus tags have two digit runs and must not be stemmed into a shared bucket.
        self.assertEqual("AT1G75820", family_index_module.symbol_stem("AT1G75820"))
        self.assertEqual("NORK", family_index_module.symbol_stem("NORK"))

    def test_uniprot_gene_name_field_precedes_trembl_accession(self) -> None:
        title = "Phytochrome A OS=Arabidopsis thaliana GN=PHYA PE=1"
        self.assertEqual("PHYA", family_index_module.subject_symbol("tr|A0A178U8X8|A0A178U8X8_ARATH", title))


class BlastInputTests(unittest.TestCase):
    def test_valid_thirteen_column_stitle_is_retained(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            hits = Path(tmp) / "hits.tsv"
            hits.write_text(
                _hit13("G_FT", "tr|A0ACM8R380|A0ACM8R380_MEDTR", 78.0, 480,
                       "FTa2 OS=Medicago truncatula PE=2 SV=1") + "\n",
                encoding="utf-8",
            )
            frame = family_index_module._read_blast_hits(hits)
            self.assertEqual("FTa2 OS=Medicago truncatula PE=2 SV=1", frame.loc[0, "subject_title"])

    def test_non_stitle_thirteenth_column_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            hits = Path(tmp) / "hits.tsv"
            hits.write_text(
                _hit("G_BAD", "tr|A0ACM8R380|A0ACM8R380_MEDTR", 78.0, 480) + "\t99.8\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "column 13 must be .*stitle"):
                family_index_module._read_blast_hits(hits)

    def test_digit_leading_multitoken_stitle_is_retained(self) -> None:
        # Issue #403: legitimate UniProt stitles that lead with a digit
        # ("26S proteasome ...", "6-phosphogluconolactonase ...") are
        # descriptive protein text, not identifier metadata, and must pass
        # the stitle gate; only single-token identifier-shaped titles fail.
        with tempfile.TemporaryDirectory() as tmp:
            hits = Path(tmp) / "hits.tsv"
            hits.write_text(
                _hit13("G_26S", "tr|Q9T090|Q9T090_ARATH", 79.0, 628,
                       "26S proteasome regulatory particle chain RPT6-like protein OS=Arabidopsis thaliana OX=3702 GN=At4g27680 PE=1 SV=1") + "\n"
                + _hit13("G_PGL", "sp|Q84WW2|6PGL3_ARATH", 88.0, 400,
                         "6-phosphogluconolactonase 3, chloroplastic OS=Arabidopsis thaliana OX=3702 GN=PGL3 PE=1 SV=1") + "\n",
                encoding="utf-8",
            )
            frame = family_index_module._read_blast_hits(hits)
            self.assertTrue(frame.loc[0, "subject_title"].startswith("26S proteasome"))
            self.assertTrue(frame.loc[1, "subject_title"].startswith("6-phosphogluconolactonase"))

    def test_digit_leading_single_token_stitle_still_fails_closed(self) -> None:
        # A digit-leading locus-tag fragment without whitespace is identifier
        # metadata, not a description, and must keep failing the stitle gate.
        with tempfile.TemporaryDirectory() as tmp:
            hits = Path(tmp) / "hits.tsv"
            hits.write_text(
                _hit13("G_BAD", "tr|A0ACM8R380|A0ACM8R380_MEDTR", 78.0, 480, "7g085020") + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "column 13 must be .*stitle"):
                family_index_module._read_blast_hits(hits)

    def test_accession_thirteenth_column_fails_closed(self) -> None:
        # Issue #318: an accession repeated as the thirteenth column (a sacc
        # slip in -outfmt) must not be adopted as a family symbol.
        with tempfile.TemporaryDirectory() as tmp:
            hits = Path(tmp) / "hits.tsv"
            hits.write_text(
                _hit13("G_BAD", "tr|A0ACM8R380|A0ACM8R380_MEDTR", 78.0, 480, "A0ACM8R380") + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "column 13 must be .*stitle"):
                family_index_module._read_blast_hits(hits)

    def test_repeated_subject_thirteenth_column_fails_closed(self) -> None:
        # Issue #318: a duplicated sseqid in the thirteenth column carries no
        # stitle information and must fail closed instead of degrading labels.
        with tempfile.TemporaryDirectory() as tmp:
            hits = Path(tmp) / "hits.tsv"
            hits.write_text(
                _hit13("G_BAD", "tr|A0ACM8R380|A0ACM8R380_MEDTR", 78.0, 480,
                       "tr|A0ACM8R380|A0ACM8R380_MEDTR") + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "column 13 must be .*stitle"):
                family_index_module._read_blast_hits(hits)

    def test_locus_tag_only_stitle_is_valid_identifier_metadata(self) -> None:
        # The pinned legume DB contains FASTA records whose whole defline is a
        # Lotus locus tag. BLAST legitimately emits stitle == sseqid for them.
        with tempfile.TemporaryDirectory() as tmp:
            hits = Path(tmp) / "hits.tsv"
            locus = "LotjaGi6g1v0243000.1"
            hits.write_text(_hit13("G_OK", locus, 78.0, 480, locus) + "\n", encoding="utf-8")
            frame = family_index_module._read_blast_hits(hits)
            self.assertEqual(locus, frame.loc[0, "subject_title"])
            self.assertEqual(locus, family_index_module.subject_symbol(locus, locus))

    def test_locus_tag_with_observed_lc_suffix_stitle_is_valid_identifier_metadata(self) -> None:
        # Issue #344: the pinned legume proteome (frozen DB legume_plant_pep3)
        # also carries identifier-only deflines with the observed "_LC.<version>"
        # low-confidence gene-model suffix. The repeated_locus_tag gate must
        # accept them as identifier-only title metadata.
        with tempfile.TemporaryDirectory() as tmp:
            hits = Path(tmp) / "hits.tsv"
            locus = "LotjaGi2g1v0260700_LC.1"
            hits.write_text(_hit13("G_OK", locus, 78.0, 480, locus) + "\n", encoding="utf-8")
            frame = family_index_module._read_blast_hits(hits)
            self.assertEqual(locus, frame.loc[0, "subject_title"])
            self.assertEqual(locus, family_index_module.subject_symbol(locus, locus))

    def test_twelve_column_blast_input_remains_supported(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            hits = Path(tmp) / "hits.tsv"
            hits.write_text(_hit("G_OK", "sp|P0CLV1|CLV1_ARATH", 80.0, 500) + "\n", encoding="utf-8")
            frame = family_index_module._read_blast_hits(hits)
            self.assertEqual("", frame.loc[0, "subject_title"])

    def test_locus_tag_or_accession_gn_is_never_adopted(self) -> None:
        # Issue #297: a GN= that is a locus tag or accession must not become
        # the family name; the protein name resolves through the curated alias
        # registry instead.
        title = "Flowering locus protein T OS=Medicago truncatula GN=MTR_7g085020 PE=3 SV=1"
        self.assertEqual("ft", family_index_module.subject_symbol("tr|F8T7T1|F8T7T1_MEDTR", title))

    def test_protein_name_alias_without_gn_resolves_to_canonical(self) -> None:
        # Issue #297 reproduction: TrEMBL defline with no GN= and the FTa2
        # paralog protein name resolves to the ft canonical alias.
        title = "FTa2 OS=Medicago truncatula PE=2 SV=1"
        self.assertEqual("ft", family_index_module.subject_symbol("tr|A0ACM8R380|A0ACM8R380_MEDTR", title))

    def test_real_symbol_gn_is_still_adopted(self) -> None:
        title = "Flowering locus protein T OS=Medicago truncatula GN=FT PE=3 SV=1"
        self.assertEqual("FT", family_index_module.subject_symbol("tr|F8T7T1|F8T7T1_MEDTR", title))

    def test_multiword_protein_name_without_alias_falls_back_to_entry_name(self) -> None:
        # A multi-word protein name without a curated alias yields no
        # candidate: the identifier-derived fallback (UniProt entry-name
        # prefix, here the accession) is kept and the manifest records
        # source none.
        title = "non-specific serine/threonine protein kinase OS=Medicago truncatula OX=3880 GN=MTR_5g034210 PE=4 SV=1"
        self.assertEqual("A0A072UD83", family_index_module.subject_symbol("tr|A0A072UD83|A0A072UD83_MEDTR", title))

    def test_locus_tag_and_accession_are_not_usable_symbols(self) -> None:
        # Issue #297: the usable-symbol predicate lives in the shared identifier
        # grammar (identifiers.is_usable_symbol), not in a local copy. The old
        # local regexes drifted: the 10-character UniProt pattern missed F8T7T1
        # and misjudged A0ACM8R380 as a usable symbol.
        from candidate_qtl.identifiers import is_usable_symbol as usable
        self.assertTrue(usable("SERKL2"))
        self.assertTrue(usable("D14"))
        self.assertTrue(usable("NORK"))
        self.assertTrue(usable("PHYA"))
        self.assertFalse(usable("MTR_5g034210"))
        self.assertFalse(usable("MTR.5g034210"))
        self.assertFalse(usable("LotjaGi0g1v0000100.1"))
        self.assertFalse(usable("A0A072UD83"))
        self.assertFalse(usable("Q9LZS3"))
        # Drift regressions: these were usable under the old local
        # regexes.
        self.assertFalse(usable("F8T7T1"))
        self.assertFalse(usable("A0ACM8R380"))
        self.assertFalse(usable("AT1G75820"))
        # "LOC12345" (no underscore) is not in the curated prefix table, so it
        # falls through as usable; only LOC_-prefixed tags (LOC_Os01g01010)
        # are rejected. This is the identifiers grammar's documented boundary.
        self.assertTrue(usable("LOC12345"))
        self.assertFalse(usable(""))

    def test_pick_best_prefers_usable_symbol_over_locus_tag(self) -> None:
        records = [
            {"subject_id": "tr|A0A072UD83|A0A072UD83_MEDTR",
             "subject_title": "non-specific serine/threonine protein kinase OS=Medicago truncatula OX=3880 GN=MTR_5g034210 PE=4 SV=1",
             "bitscore": 500.0},
            {"subject_id": "sp|Q8L4H4|NORK_MEDTR",
             "subject_title": "Nodulation receptor kinase OS=Medicago truncatula OX=3880 GN=NORK PE=1 SV=2",
             "bitscore": 480.0},
        ]
        best = family_index_module._pick_best_record(records, {})
        self.assertEqual("sp|Q8L4H4|NORK_MEDTR", str(best["subject_id"]))

    def test_pick_best_keeps_best_when_all_symbols_are_locus_tags(self) -> None:
        records = [
            {"subject_id": "tr|A0A072UD83|A0A072UD83_MEDTR",
             "subject_title": "non-specific kinase OS=Medicago truncatula OX=3880 GN=MTR_5g034210 PE=4 SV=1",
             "bitscore": 500.0},
            {"subject_id": "tr|G7K1D0|G7K1D0_MEDTR",
             "subject_title": "non-specific kinase OS=Medicago truncatula OX=3880 GN=MTR_5g034210 PE=4 SV=1",
             "bitscore": 400.0},
        ]
        best = family_index_module._pick_best_record(records, {})
        self.assertEqual("tr|A0A072UD83|A0A072UD83_MEDTR", str(best["subject_id"]))

    def test_pick_best_promotes_an_alias_resolved_symbol_over_a_locus_tag(self) -> None:
        # Issue #297: the promotion predicate uses the *derived* symbol. The
        # FTa2 hit has no usable GN=, but its protein name resolves to the ft
        # canonical alias, so it is a usable-symbol hit and is promoted
        # (480/500 >= MIN_PROMOTED_BITSCORE_RATIO).
        records = [
            {"subject_id": "tr|A0A072UD83|A0A072UD83_MEDTR",
             "subject_title": "non-specific serine/threonine protein kinase OS=Medicago truncatula OX=3880 GN=MTR_5g034210 PE=4 SV=1",
             "bitscore": 500.0},
            {"subject_id": "tr|A0ACM8R380|A0ACM8R380_MEDTR",
             "subject_title": "FTa2 OS=Medicago truncatula PE=2 SV=1",
             "bitscore": 480.0},
        ]
        best = family_index_module._pick_best_record(records, {})
        self.assertEqual("tr|A0ACM8R380|A0ACM8R380_MEDTR", str(best["subject_id"]))
        self.assertEqual("ft", family_index_module.subject_symbol(
            str(best["subject_id"]), str(best["subject_title"])))

    def test_pick_best_does_not_promote_a_hit_whose_only_symbol_is_an_accession(self) -> None:
        # A hit whose GN= is an accession (A0ACM8R380) derives no usable
        # symbol: its protein name has no alias and the fallback entry name is
        # the accession. It must never be promoted over the locus-tag top hit.
        records = [
            {"subject_id": "tr|A0A072UD83|A0A072UD83_MEDTR",
             "subject_title": "non-specific serine/threonine protein kinase OS=Medicago truncatula OX=3880 GN=MTR_5g034210 PE=4 SV=1",
             "bitscore": 500.0},
            {"subject_id": "tr|G7K1D0|G7K1D0_MEDTR",
             "subject_title": "non-specific kinase OS=Medicago truncatula OX=3880 GN=A0ACM8R380 PE=4 SV=1",
             "bitscore": 480.0},
        ]
        best = family_index_module._pick_best_record(records, {})
        self.assertEqual("tr|A0A072UD83|A0A072UD83_MEDTR", str(best["subject_id"]))

    def test_pick_best_keeps_highest_symbol_hit(self) -> None:
        records = [
            {"subject_id": "sp|Q8L4H4|NORK_MEDTR",
             "subject_title": "Nodulation receptor kinase OS=Medicago truncatula OX=3880 GN=NORK PE=1 SV=2",
             "bitscore": 500.0},
            {"subject_id": "sp|P0CLV1|CLV1_ARATH",
             "subject_title": "Receptor kinase CLV1 OS=Arabidopsis thaliana OX=3702 GN=CLV1 PE=1 SV=2",
             "bitscore": 400.0},
        ]
        best = family_index_module._pick_best_record(records, {})
        self.assertEqual("sp|Q8L4H4|NORK_MEDTR", str(best["subject_id"]))

    def test_pick_best_keeps_true_best_when_symbol_hit_is_far_below_the_bound(self) -> None:
        # NORK's bitscore is 100/500 = 0.2 of the top hit, well under
        # MIN_PROMOTED_BITSCORE_RATIO (0.5) -- promoting it would present a
        # weak, distantly related hit as "best". The true top hit (a locus
        # tag) must be kept instead.
        records = [
            {"subject_id": "tr|A0A072UD83|A0A072UD83_MEDTR",
             "subject_title": "non-specific serine/threonine protein kinase OS=Medicago truncatula OX=3880 GN=MTR_5g034210 PE=4 SV=1",
             "bitscore": 500.0},
            {"subject_id": "sp|Q8L4H4|NORK_MEDTR",
             "subject_title": "Nodulation receptor kinase OS=Medicago truncatula OX=3880 GN=NORK PE=1 SV=2",
             "bitscore": 100.0},
        ]
        best = family_index_module._pick_best_record(records, {})
        self.assertEqual("tr|A0A072UD83|A0A072UD83_MEDTR", str(best["subject_id"]))


class FamilyIndexTests(unittest.TestCase):
    def _inputs(self, root: Path) -> tuple[Path, Path]:
        hits = root / "hits.tsv"
        hits.write_text("\n".join([
            _hit("G_CLV_A", "sp|P0CLV1|CLV1_ARATH", 62.0, 400),
            _hit("G_CLV_A", "sp|P0CLV2|CLV2_ARATH", 41.0, 210),
            _hit("G_CLV_B", "sp|P0CLV1|CLV1_ARATH", 58.0, 380),
            _hit("G_CLV_B", "sp|P0CLV2|CLV2_ARATH", 40.0, 360),
            _hit("G_CLV_FAR", "sp|P0CLV1|CLV1_ARATH", 55.0, 300),
            _hit("G_PIN_A", "sp|P0PIN1|PIN1_ARATH", 70.0, 500),
            _hit("G_KIN_A", "sp|P0KIN1|KIN1_ARATH", 45.0, 250),
            _hit("G_KIN_B", "sp|P0KIN2|KIN2_ARATH", 44.0, 240),
            _hit("G_KIN_C", "sp|P0KIN3|KIN3_ARATH", 43.0, 230),
            _hit("G_KIN_D", "sp|P0KIN4|KIN4_ARATH", 42.0, 220),
            _hit("G_LOWQ", "sp|P0XXX1|XXX1_ARATH", 12.0, 40),
        ]) + "\n", encoding="utf-8")
        genes = root / "genes.tsv"
        genes.write_text("\n".join([
            "gene_id\tchrom\tstart\tend",
            "G_CLV_A\tchr1\t1000\t2000",
            "G_CLV_B\tchr1\t1500000\t1501000",
            "G_CLV_FAR\tchr9\t100\t900",
            "G_PIN_A\tchr1\t900000\t901000",
            "G_KIN_A\tchr1\t20000\t21000",
            "G_KIN_B\tchr9\t20000\t21000",
            "G_KIN_C\tchr9\t40000\t41000",
            "G_KIN_D\tchr9\t60000\t61000",
            "G_LOWQ\tchr1\t50000\t51000",
            "G_DARK\tchr1\t60000\t61000",
        ]) + "\n", encoding="utf-8")
        return hits, genes

    def _build(self, root: Path) -> Path:
        hits, genes = self._inputs(root)
        out = root / "index"
        self.assertEqual(0, main(["family-index", "--blast-hits", str(hits), "--genes", str(genes),
                                  "--out-dir", str(out), "--reference-label", "TestRef"]))
        return out

    def test_index_assigns_fine_and_coarse_families(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._build(Path(tmp))
            index = pd.read_csv(out / "family_index.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            self.assertEqual("sp|P0CLV1|CLV1_ARATH", index.loc["G_CLV_A", "best_subject"])
            self.assertEqual("CLV1", index.loc["G_CLV_A", "best_subject_symbol"])
            self.assertEqual("CLV", index.loc["G_CLV_A", "best_subject_stem"])
            self.assertIn("CLV", index.loc["G_CLV_A", "family_stems"])
            # top2 ratio is available as the continuous axis for later tie-breaking.
            self.assertAlmostEqual(210 / 400, float(index.loc["G_CLV_A", "top2_bitscore_ratio"]), places=4)
            self.assertAlmostEqual(360 / 380, float(index.loc["G_CLV_B", "top2_bitscore_ratio"]), places=4)
            # No promotion happened for either gene (the top hit already has a
            # usable symbol), so provenance stays at its identity default.
            self.assertEqual(1, int(float(index.loc["G_CLV_A", "symbol_source_rank"])))
            self.assertAlmostEqual(1.0, float(index.loc["G_CLV_A", "symbol_source_bitscore_ratio"]), places=4)
            self.assertEqual(1, int(float(index.loc["G_CLV_B", "symbol_source_rank"])))
            self.assertAlmostEqual(1.0, float(index.loc["G_CLV_B", "symbol_source_bitscore_ratio"]), places=4)

    def test_legume_locus_tag_does_not_capture_family_name(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hits = root / "hits.tsv"
            hits.write_text("\n".join([
                _hit13("G_NORK", "tr|A0A072UD83|A0A072UD83_MEDTR", 80.0, 500,
                       "non-specific serine/threonine protein kinase OS=Medicago truncatula OX=3880 GN=MTR_5g034210 PE=4 SV=1"),
                _hit13("G_NORK", "sp|Q8L4H4|NORK_MEDTR", 78.0, 480,
                       "Nodulation receptor kinase OS=Medicago truncatula OX=3880 GN=NORK PE=1 SV=2"),
                _hit13("G_NORK", "sp|P0CLV1|CLV1_ARATH", 50.0, 300,
                       "Receptor kinase CLV1 OS=Arabidopsis thaliana OX=3702 GN=CLV1 PE=1 SV=2"),
            ]) + "\n", encoding="utf-8")
            genes = root / "genes.tsv"
            genes.write_text("gene_id\tchrom\tstart\tend\nG_NORK\tchr1\t1000\t2000\n", encoding="utf-8")
            out = root / "index"
            self.assertEqual(0, main(["family-index", "--blast-hits", str(hits), "--genes", str(genes),
                                      "--out-dir", str(out)]))
            index = pd.read_csv(out / "family_index.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            # #647 also limits the identity deficit to one percentage point.
            # NORK trails by two points here, so annotation cannot displace the
            # true sequence top even though its bitscore is near-equivalent.
            self.assertEqual("A0A072UD83", index.loc["G_NORK", "best_subject_symbol"])
            self.assertEqual("A0A072UD83", index.loc["G_NORK", "discovery_bucket_id"])
            # Issue #265, reproduction shape 1 (500 MTR locus tag / 480 NORK /
            # 300 CLV1): the promoted best (NORK, 480) is not the true top hit
            # (MTR, 500) any more, so second_bitscore must come from the
            # highest-scoring hit that is NOT the promoted record (here, the
            # excluded true top, 500) rather than a stale index into the
            # original ordering -- and the resulting ratio must land in (0, 1].
            best_bitscore = float(index.loc["G_NORK", "best_bitscore"])
            second_bitscore = float(index.loc["G_NORK", "second_bitscore"])
            ratio = float(index.loc["G_NORK", "top2_bitscore_ratio"])
            self.assertEqual(500.0, best_bitscore)
            self.assertEqual(480.0, second_bitscore)
            self.assertGreater(ratio, 0.0)
            self.assertLessEqual(ratio, 1.0)
            self.assertAlmostEqual(480 / 500, ratio, places=4)
            # Provenance: the symbol/best row came from rank 2 (NORK, the
            # second-highest-scoring original hit), at 480/500 of the true
            # top hit's bitscore.
            self.assertEqual(1, int(float(index.loc["G_NORK", "symbol_source_rank"])))
            self.assertAlmostEqual(1.0, float(index.loc["G_NORK", "symbol_source_bitscore_ratio"]), places=4)

    def test_legume_locus_tag_promotion_two_deep_keeps_ratio_bounded(self) -> None:
        # Issue #265, reproduction shape 2: two locus-tag/accession hits
        # (MTR, LotjaGi) outscore the only usable-symbol hit (NORK, rank 3).
        # Before the fix, second_bitscore stayed pinned to the original
        # records[1] (LotjaGi, 480) even though 480 > the promoted best's
        # bitscore (300), giving a ratio of 1.6 -- outside the documented
        # [0, 1] domain.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hits = root / "hits.tsv"
            hits.write_text("\n".join([
                _hit13("G_X", "tr|A0A072UD83|A0A072UD83_MEDTR", 80.0, 500,
                       "non-specific serine/threonine protein kinase OS=Medicago truncatula OX=3880 GN=MTR_5g034210 PE=4 SV=1"),
                # Issue #297: GN=LotjaGi0g1v0000100 is a locus tag and is never
                # adopted; the entry-name fallback (a realistic TrEMBL
                # accession, G7K1D0) is not a usable symbol either, so this
                # hit still cannot be promoted.
                _hit13("G_X", "tr|G7K1D0|G7K1D0_LOTJA", 79.0, 480,
                       "hypothetical protein OS=Lotus japonicus OX=34305 GN=LotjaGi0g1v0000100 PE=4 SV=1"),
                _hit13("G_X", "sp|Q8L4H4|NORK_MEDTR", 60.0, 300,
                       "Nodulation receptor kinase OS=Medicago truncatula OX=3880 GN=NORK PE=1 SV=2"),
            ]) + "\n", encoding="utf-8")
            genes = root / "genes.tsv"
            genes.write_text("gene_id\tchrom\tstart\tend\nG_X\tchr1\t1000\t2000\n", encoding="utf-8")
            out = root / "index"
            self.assertEqual(0, main(["family-index", "--blast-hits", str(hits), "--genes", str(genes),
                                      "--out-dir", str(out)]))
            index = pd.read_csv(out / "family_index.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            # #647 replaces the old half-bitscore symbol promotion with a 95%
            # homology-equivalence window. NORK retains only 60% here, so the
            # true sequence top remains selected.
            self.assertEqual("A0A072UD83", index.loc["G_X", "best_subject_symbol"])
            best_bitscore = float(index.loc["G_X", "best_bitscore"])
            second_bitscore = float(index.loc["G_X", "second_bitscore"])
            ratio = float(index.loc["G_X", "top2_bitscore_ratio"])
            self.assertEqual(500.0, best_bitscore)
            self.assertEqual(480.0, second_bitscore)
            self.assertGreater(ratio, 0.0)
            self.assertLessEqual(ratio, 1.0)
            self.assertAlmostEqual(480 / 500, ratio, places=4)
            self.assertEqual(1, int(float(index.loc["G_X", "symbol_source_rank"])))
            self.assertAlmostEqual(1.0, float(index.loc["G_X", "symbol_source_bitscore_ratio"]), places=4)

    def test_legume_locus_tag_promotion_is_bounded_by_minimum_ratio(self) -> None:
        # The only usable-symbol hit (NORK, 100) is far below
        # MIN_PROMOTED_BITSCORE_RATIO of the true top hit (MTR, 500): 0.2 < 0.5.
        # Promotion must not happen; the true best hit (the locus tag) is kept,
        # and provenance reports no promotion.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hits = root / "hits.tsv"
            hits.write_text("\n".join([
                _hit13("G_FAR", "tr|A0A072UD83|A0A072UD83_MEDTR", 80.0, 500,
                       "non-specific serine/threonine protein kinase OS=Medicago truncatula OX=3880 GN=MTR_5g034210 PE=4 SV=1"),
                _hit13("G_FAR", "sp|Q8L4H4|NORK_MEDTR", 30.0, 100,
                       "Nodulation receptor kinase OS=Medicago truncatula OX=3880 GN=NORK PE=1 SV=2"),
            ]) + "\n", encoding="utf-8")
            genes = root / "genes.tsv"
            genes.write_text("gene_id\tchrom\tstart\tend\nG_FAR\tchr1\t1000\t2000\n", encoding="utf-8")
            out = root / "index"
            self.assertEqual(0, main(["family-index", "--blast-hits", str(hits), "--genes", str(genes),
                                      "--out-dir", str(out)]))
            index = pd.read_csv(out / "family_index.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            self.assertEqual("tr|A0A072UD83|A0A072UD83_MEDTR", index.loc["G_FAR", "best_subject"])
            # Issue #297: the GN= locus tag is never adopted; without a
            # candidate the identifier-derived entry-name fallback (the
            # accession) is kept, and the manifest records source none.
            self.assertEqual("A0A072UD83", index.loc["G_FAR", "best_subject_symbol"])
            self.assertEqual(500.0, float(index.loc["G_FAR", "best_bitscore"]))
            self.assertEqual(1, int(float(index.loc["G_FAR", "symbol_source_rank"])))
            self.assertAlmostEqual(1.0, float(index.loc["G_FAR", "symbol_source_bitscore_ratio"]), places=4)
            manifest = json.loads((out / "family_index_manifest.json").read_text(encoding="utf-8"))
            provenance = manifest["best_subject_symbol_provenance"]
            self.assertEqual(1, sum(provenance["counts_by_source"].values()))
            self.assertEqual(1, provenance["counts_by_source"]["none"])
            self.assertIn("identifier_derived_bucket_rule", provenance)

    def test_accession_deflines_resolve_to_curated_alias_family(self) -> None:
        # Issue #297 reproduction (pea GIGAS family): the F8T7T1 and
        # A0ACM8R380 deflines must never put the accession into the family
        # name. Both resolve to the ft bucket via the curated alias registry
        # (GN= is a locus tag / absent); a multi-word protein name without an
        # alias falls back to the identifier-derived entry name and is
        # auditable as source none in the manifest.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hits = root / "hits.tsv"
            hits.write_text("\n".join([
                _hit13("G_FT_A", "tr|F8T7T1|F8T7T1_MEDTR", 80.0, 500,
                       "Flowering locus protein T OS=Medicago truncatula GN=MTR_7g085020 PE=3 SV=1"),
                _hit13("G_FT_B", "tr|A0ACM8R380|A0ACM8R380_MEDTR", 78.0, 480,
                       "FTa2 OS=Medicago truncatula PE=2 SV=1"),
                _hit13("G_KIN", "tr|A0ACM8R380|A0ACM8R380_MEDTR", 60.0, 300,
                       "non-specific serine/threonine protein kinase OS=Medicago truncatula OX=3880 GN=MTR_7g085020 PE=4 SV=1"),
            ]) + "\n", encoding="utf-8")
            genes = root / "genes.tsv"
            genes.write_text("gene_id\tchrom\tstart\tend\nG_FT_A\tchr1\t1000\t2000\nG_FT_B\tchr1\t3000\t4000\nG_KIN\tchr1\t5000\t6000\n",
                             encoding="utf-8")
            out = root / "index"
            self.assertEqual(0, main(["family-index", "--blast-hits", str(hits), "--genes", str(genes),
                                      "--out-dir", str(out)]))
            index = pd.read_csv(out / "family_index.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            self.assertEqual("ft", index.loc["G_FT_A", "best_subject_symbol"])
            self.assertEqual("ft", index.loc["G_FT_A", "discovery_bucket_id"])
            self.assertEqual("ft", index.loc["G_FT_B", "best_subject_symbol"])
            self.assertEqual("ft", index.loc["G_FT_B", "discovery_bucket_id"])
            self.assertEqual("A0ACM8R380", index.loc["G_KIN", "best_subject_symbol"])
            manifest = json.loads((out / "family_index_manifest.json").read_text(encoding="utf-8"))
            provenance = manifest["best_subject_symbol_provenance"]
            self.assertEqual(2, provenance["counts_by_source"]["curated_alias"])
            self.assertEqual(1, provenance["counts_by_source"]["none"])
            self.assertIn("identifier_derived_bucket_rule", provenance)
            self.assertEqual("1.4", manifest["schema_version"])
            self.assertIn("curated alias", manifest["family_definition"]["blast_symbol_policy"])

    def test_low_quality_and_missing_hits_become_an_explicit_blind_spot(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._build(Path(tmp))
            index = pd.read_csv(out / "family_index.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            # 12% identity is below the quality floor, so the gene stays unassigned...
            self.assertEqual("dark_no_qualifying_hit", index.loc["G_LOWQ", "assignment_status"])
            self.assertEqual("dark_no_qualifying_hit", index.loc["G_DARK", "assignment_status"])
            # ...but remains in the universe denominator.
            manifest = json.loads((out / "family_index_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(10, manifest["gene_count"])
            self.assertEqual(2, manifest["dark_count"])
            self.assertEqual(8, manifest["assigned_count"])
            self.assertIn("blind_spot", manifest)

    def test_curated_family_map_overrides_the_symbol_stem(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hits, genes = self._inputs(root)
            family_map = root / "families.tsv"
            family_map.write_text("subject_id\tfamily_id\tsymbol\n"
                                  "sp|P0CLV1|CLV1_ARATH\tPF00560_LRR\tCLV1\n", encoding="utf-8")
            out = root / "index"
            self.assertEqual(0, main(["family-index", "--blast-hits", str(hits), "--genes", str(genes),
                                      "--out-dir", str(out), "--family-map", str(family_map)]))
            index = pd.read_csv(out / "family_index.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            self.assertEqual("PF00560_LRR", index.loc["G_CLV_A", "best_subject_stem"])
            manifest = json.loads((out / "family_index_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(manifest["family_definition"]["override_source"])

    def test_family_index_cli_records_optional_case_id(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hits, genes = self._inputs(root)
            out = root / "index"
            self.assertEqual(0, main(["family-index", "--blast-hits", str(hits), "--genes", str(genes),
                                      "--out-dir", str(out), "--case-id", "pea_r"]))
            manifest = json.loads((out / "family_index_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("pea_r", manifest["case_id"])

    def test_issue_184_discovery_bucket_vs_evolutionary_family(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hits, genes = self._inputs(root)
            # Without a curated map the coarse value is a symbol-stem discovery
            # bucket; the evolutionary family layer stays unresolved.
            out = self._build(root)
            index = pd.read_csv(out / "family_index.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            self.assertEqual("CLV", index.loc["G_CLV_A", "discovery_bucket_id"])
            self.assertEqual("symbol_stem", index.loc["G_CLV_A", "discovery_bucket_method"])
            self.assertEqual("", index.loc["G_CLV_A", "evolutionary_family_id"])
            self.assertEqual("", index.loc["G_CLV_A", "evolutionary_family_source"])
            self.assertEqual("unresolved", index.loc["G_CLV_A", "evolutionary_family_status"])
            self.assertEqual("unresolved", index.loc["G_DARK", "evolutionary_family_status"])
            # With a curated map the same value is carried by both layers, but
            # with explicit, separated meanings: bucket method curated_keyword,
            # evolutionary family supported via source curated_family.
            family_map = root / "families.tsv"
            family_map.write_text("subject_id\tfamily_id\tsymbol\n"
                                  "sp|P0CLV1|CLV1_ARATH\tLRR-RLK_CLV\tCLV1\n"
                                  "sp|P0CLV2|CLV2_ARATH\tLRR-RLK_CLV\tCLV2\n", encoding="utf-8")
            mapped = root / "index_mapped"
            self.assertEqual(0, main(["family-index", "--blast-hits", str(hits), "--genes", str(genes),
                                      "--out-dir", str(mapped), "--family-map", str(family_map)]))
            index = pd.read_csv(mapped / "family_index.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            self.assertEqual("LRR-RLK_CLV", index.loc["G_CLV_A", "discovery_bucket_id"])
            self.assertEqual("curated_keyword", index.loc["G_CLV_A", "discovery_bucket_method"])
            self.assertEqual("LRR-RLK_CLV", index.loc["G_CLV_A", "evolutionary_family_id"])
            self.assertEqual("curated_family", index.loc["G_CLV_A", "evolutionary_family_source"])
            self.assertEqual("supported", index.loc["G_CLV_A", "evolutionary_family_status"])
            manifest = json.loads((mapped / "family_index_manifest.json").read_text(encoding="utf-8"))
            self.assertIn("semantics_184", manifest["family_definition"])

    def test_manifest_records_the_compute_budget(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._build(Path(tmp))
            manifest = json.loads((out / "family_index_manifest.json").read_text(encoding="utf-8"))
            self.assertIn("wall_seconds", manifest["compute"])
            self.assertIn("peak_python_memory_mb", manifest["compute"])
            self.assertIn("once per reference", manifest["compute"]["note"])

    def test_optional_stitle_gn_and_reviewed_tie_policy(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hits = root / "hits.tsv"
            trembl_title = "Phytochrome A OS=Arabidopsis thaliana GN=PHYA PE=1"
            swiss_title = "Phytochrome A OS=Arabidopsis thaliana GN=PHYA PE=1 SV=2"
            hits.write_text(
                _hit("G_PHY", "tr|A0A178U8X8|A0A178U8X8_ARATH", 60.0, 400) + "\t" + trembl_title + "\n"
                + _hit("G_PHY", "sp|P14712|PHYA_ARATH", 60.0, 400) + "\t" + swiss_title + "\n",
                encoding="utf-8",
            )
            genes = root / "genes.tsv"
            genes.write_text("gene_id\tchrom\tstart\tend\nG_PHY\tchr1\t1\t100\n", encoding="utf-8")
            out = root / "out"
            self.assertEqual(0, main(["family-index", "--blast-hits", str(hits), "--genes", str(genes), "--out-dir", str(out)]))
            row = pd.read_csv(out / "family_index.tsv", sep="\t", dtype=str).fillna("").iloc[0]
            self.assertEqual("sp|P14712|PHYA_ARATH", row["best_subject"])
            self.assertEqual("PHYA", row["best_subject_symbol"])
            self.assertEqual("PHYA", row["best_subject_stem"])
            self.assertTrue(family_index_module.families_match("Phytochrome (PHY)", row["best_subject_symbol"]))
            manifest = json.loads((out / "family_index_manifest.json").read_text(encoding="utf-8"))
            self.assertIn("GN=<symbol>", manifest["family_definition"]["blast_symbol_policy"])

    def test_descriptive_family_alias_matches_symbol_stem_during_falsification(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hits = root / "hits.tsv"
            hits.write_text(_hit("G_PHY", "sp|P14712|PHYA_ARATH", 60.0, 400) + "\n", encoding="utf-8")
            genes = root / "genes.tsv"
            genes.write_text("gene_id\tchrom\tstart\tend\nG_PHY\tchr1\t1\t100\n", encoding="utf-8")
            index_out = root / "index"
            self.assertEqual(0, main(["family-index", "--blast-hits", str(hits), "--genes", str(genes), "--out-dir", str(index_out)]))
            hypotheses = root / "hypotheses.tsv"
            hypotheses.write_text(
                "hypothesis_id\tpathway_or_process\tgene_family\treference_gene_name\tproposed_stable_id\n"
                "hyp_phy\tlight response\tPhytochrome (PHY)\t\t\n",
                encoding="utf-8",
            )
            out = root / "falsify"
            self.assertEqual(0, main(["pathway-interval-falsify", "--family-index", str(index_out / "family_index.tsv"),
                                      "--hypotheses", str(hypotheses), "--interval", "chr1:1-200", "--out-dir", str(out)]))
            rows = pd.read_csv(out / "hypothesis_interval_falsification.tsv", sep="\t", dtype=str).fillna("")
            coarse = rows[(rows["hypothesis_id"] == "hyp_phy") & (rows["level"] == "coarse")].iloc[0]
            self.assertEqual("1", coarse["interval_members"])
            fine = rows[(rows["hypothesis_id"] == "hyp_phy") & (rows["level"] == "fine")].iloc[0]
            self.assertEqual("not_computed_no_fine_identifier", fine["verdict"])

    def test_identifier_only_stitle_e2e_lotja_lc_and_variant_ids(self) -> None:
        # Issue #344 E2E: the pinned legume proteome (frozen DB
        # legume_plant_pep3) contains FASTA records whose whole defline is an
        # identifier -- the observed "_LC.<version>" low-confidence gene-model
        # suffix and the letter-leading C1v/M1v assembly variants. BLAST emits
        # stitle == sseqid for them, and the 13-column family-index build must
        # succeed instead of failing closed. The identifier-only title is kept
        # as metadata but never becomes a biological family symbol: the
        # fallback label is the identifier itself with provenance source none.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hits = root / "hits.tsv"
            hits.write_text("\n".join([
                _hit13("G_LC", "LotjaGi2g1v0260700_LC.1", 80.0, 500, "LotjaGi2g1v0260700_LC.1"),
                _hit13("G_LC2", "LotjaGi0g1v0001000_LC.2", 79.0, 490, "LotjaGi0g1v0001000_LC.2"),
                _hit13("G_C1V", "LotjaGiC1v0000600.1", 78.0, 480, "LotjaGiC1v0000600.1"),
                _hit13("G_M1V", "LotjaGiM1v0000600.1", 76.0, 460, "LotjaGiM1v0000600.1"),
            ]) + "\n", encoding="utf-8")
            genes = root / "genes.tsv"
            genes.write_text(
                "gene_id\tchrom\tstart\tend\n"
                "G_LC\tchr1\t1000\t2000\nG_LC2\tchr1\t3000\t4000\n"
                "G_C1V\tchr1\t5000\t6000\nG_M1V\tchr1\t7000\t8000\n",
                encoding="utf-8",
            )
            out = root / "index"
            self.assertEqual(0, main(["family-index", "--blast-hits", str(hits), "--genes", str(genes),
                                      "--out-dir", str(out)]))
            index = pd.read_csv(out / "family_index.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            self.assertEqual("assigned", index.loc["G_LC", "assignment_status"])
            self.assertEqual("assigned", index.loc["G_LC2", "assignment_status"])
            self.assertEqual("assigned", index.loc["G_C1V", "assignment_status"])
            self.assertEqual("assigned", index.loc["G_M1V", "assignment_status"])
            self.assertEqual("LotjaGi2g1v0260700_LC.1", index.loc["G_LC", "best_subject_symbol"])
            self.assertEqual("LotjaGi0g1v0001000_LC.2", index.loc["G_LC2", "best_subject_symbol"])
            self.assertEqual("LotjaGiC1v0000600.1", index.loc["G_C1V", "best_subject_symbol"])
            self.assertEqual("LotjaGiM1v0000600.1", index.loc["G_M1V", "best_subject_symbol"])
            manifest = json.loads((out / "family_index_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(4, manifest["assigned_count"])
            self.assertEqual(4, manifest["best_subject_symbol_provenance"]["counts_by_source"]["none"])


class FalsifyTests(unittest.TestCase):
    def _index_and_hypotheses(self, root: Path) -> tuple[Path, Path]:
        out = FamilyIndexTests()._build(root)
        hypotheses = root / "hypotheses.tsv"
        hypotheses.write_text("\n".join([
            "hypothesis_id\tpathway_or_process\tgene_family\treference_gene_name\tproposed_stable_id",
            # CLV: 2 of 3 genome-wide members sit in the interval -> concentrated.
            "hyp_001\tCLAVATA-WUSCHEL\tCLV\tCLV1\tsp|P0CLV1|CLV1_ARATH",
            # PIN: present in the interval but unremarkable.
            "hyp_002\tauxin transport\tPIN\tPIN1\tsp|P0PIN1|PIN1_ARATH",
            # KIN: mostly outside the interval.
            "hyp_003\tgeneric kinase signalling\tKIN\tKIN1\tsp|P0KIN1|KIN1_ARATH",
            # Absent family -> falsified by the interval at the coarse level.
            "hyp_004\tflowering time\tFT\tFT\tsp|P0FT001|FT_ARATH",
        ]) + "\n", encoding="utf-8")
        return out / "family_index.tsv", hypotheses

    def _run(self, root: Path, interval: str = "chr1:1-2000000") -> Path:
        index, hypotheses = self._index_and_hypotheses(root)
        out = root / "falsify"
        self.assertEqual(0, main(["pathway-interval-falsify", "--family-index", str(index),
                                  "--hypotheses", str(hypotheses), "--interval", interval,
                                  "--out-dir", str(out), "--case-id", "demo"]))
        return out

    def test_absent_family_is_held_not_evaluated_never_falsified(self) -> None:
        """Zero resolution of the declared identifiers in the whole index is an
        index-vocabulary gap (issue #599), so the hypothesis stays an
        evaluable hold and is never converted into a biological absence."""
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            frame = pd.read_csv(out / "hypothesis_interval_falsification.tsv", sep="\t", dtype=str).fillna("")
            coarse = frame[frame.level == "coarse"].set_index("hypothesis_id")
            self.assertEqual("not_evaluated_family_unresolved_in_index", coarse.loc["hyp_004", "verdict"])
            self.assertEqual("", coarse.loc["hyp_004", "resolved_by"])
            self.assertEqual("0", coarse.loc["hyp_004", "genome_wide_members"])
            self.assertIn("評価不能", coarse.loc["hyp_004", "notes"])
            # The cascade skips the fine computation for a held hypothesis.
            fine = frame[frame.level == "fine"].set_index("hypothesis_id")
            self.assertEqual("not_computed_coarse_family_unresolved_in_index", fine.loc["hyp_004", "verdict"])
            reviews = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            reasons = set(reviews["reason"])
            self.assertIn("hypothesis_family_not_detected_genome_wide", reasons)
            # The review row stays machine-readable-linked to the new status.
            row = reviews[reviews.reason == "hypothesis_family_not_detected_genome_wide"].iloc[0]
            self.assertEqual("not_evaluated_family_unresolved_in_index", row.status)
            manifest = json.loads((out / "pathway_interval_falsify_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(0, manifest["coarse_falsified"])
            self.assertEqual(1, manifest["coarse_family_unresolved_in_index"])

    def _absent_from_interval(self, root: Path, *extra: str) -> tuple[pd.DataFrame, dict]:
        """chr9 holds KIN_B/C/D and CLV_FAR, so an interval there excludes PIN
        entirely even though the PIN family exists in the genome."""
        index, hypotheses = self._index_and_hypotheses(root)
        out = root / f"falsify{len(extra)}"
        self.assertEqual(0, main(["pathway-interval-falsify", "--family-index", str(index),
                                  "--hypotheses", str(hypotheses), "--interval", "chr9:1-100000",
                                  "--out-dir", str(out), *extra]))
        frame = pd.read_csv(out / "hypothesis_interval_falsification.tsv", sep="\t", dtype=str).fillna("")
        manifest = json.loads((out / "pathway_interval_falsify_manifest.json").read_text(encoding="utf-8"))
        return frame, manifest

    def _rescue_table(self, root: Path, name: str, hits: int) -> Path:
        path = root / name
        path.write_text("\n".join([
            "query_protein_id\thit_count\tmodel_class\tchrom\tstart\tend",
            f"PIN1\t{hits}\tcomplete\tchr9\t5000\t7000",
        ]) + "\n", encoding="utf-8")
        return path

    def test_absence_without_a_sequence_search_is_held_not_falsified(self) -> None:
        """An annotation gap must not falsify a correct family: that is how the
        phenotype-first route drifts to ever more distant candidates."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            frame, _ = self._absent_from_interval(root)
            coarse = frame[frame.level == "coarse"].set_index("hypothesis_id")
            self.assertEqual("absent_annotation_only", coarse.loc["hyp_002", "verdict"])
            self.assertEqual("1", coarse.loc["hyp_002", "genome_wide_members"])
            self.assertEqual("0", coarse.loc["hyp_002", "interval_members"])
            # a hold is neither falsified nor surviving
            self.assertFalse(coarse.loc["hyp_002", "verdict"].startswith("falsified"))
            self.assertFalse(coarse.loc["hyp_002", "verdict"].startswith("surviving"))
            fine = frame[frame.level == "fine"].set_index("hypothesis_id")
            self.assertEqual("not_computed_coarse_absent_annotation_only", fine.loc["hyp_002", "verdict"])

    def test_absence_confirmed_by_a_sequence_search_is_falsified(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            frame, _ = self._absent_from_interval(root, "--annotation-rescue",
                                                  str(self._rescue_table(root, "rescue_miss.tsv", 0)))
            coarse = frame[frame.level == "coarse"].set_index("hypothesis_id")
            self.assertEqual("falsified_no_member_in_interval", coarse.loc["hyp_002", "verdict"])
            self.assertIn("配列レベルで確認", coarse.loc["hyp_002", "notes"])

    def test_a_sequence_hit_rescues_the_hypothesis_from_an_annotation_gap(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            frame, _ = self._absent_from_interval(root, "--annotation-rescue",
                                                  str(self._rescue_table(root, "rescue_hit.tsv", 1)))
            coarse = frame[frame.level == "coarse"].set_index("hypothesis_id")
            self.assertEqual("surviving_sequence_only_annotation_gap", coarse.loc["hyp_002", "verdict"])
            self.assertIn("chr9:5000-7000", coarse.loc["hyp_002", "notes"])
            reasons = set(pd.read_csv(root / "falsify2" / "review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("family_present_in_sequence_but_missing_from_annotation", reasons)

    def test_present_but_unenriched_family_survives(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            frame = pd.read_csv(out / "hypothesis_interval_falsification.tsv", sep="\t", dtype=str).fillna("")
            coarse = frame[frame.level == "coarse"].set_index("hypothesis_id")
            # Recall first: presence is enough to survive; enrichment is credit only.
            self.assertTrue(coarse.loc["hyp_002", "verdict"].startswith("surviving"))
            self.assertEqual("1", coarse.loc["hyp_002", "interval_members"])

    def test_concentrated_family_is_more_enriched_than_a_dispersed_one(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            frame = pd.read_csv(out / "hypothesis_interval_falsification.tsv", sep="\t", dtype=str).fillna("")
            coarse = frame[frame.level == "coarse"].set_index("hypothesis_id")
            clv = float(coarse.loc["hyp_001", "conservative_p"])
            kin = float(coarse.loc["hyp_003", "conservative_p"])
            self.assertLess(clv, kin)
            self.assertEqual("3", coarse.loc["hyp_001", "genome_wide_members"])
            self.assertEqual("2", coarse.loc["hyp_001", "interval_members"])

    def test_tandem_members_are_collapsed_before_the_conservative_p(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index, hypotheses = self._index_and_hypotheses(root)
            # A narrow interval holding two CLV genes 500 bp apart: linkage, not enrichment.
            frame = pd.read_csv(index, sep="\t", dtype=str).fillna("")
            frame.loc[frame.gene_id == "G_CLV_B", ["chrom", "start", "end"]] = ["chr1", "2500", "3400"]
            frame.to_csv(index, sep="\t", index=False, lineterminator="\n")
            out = root / "falsify"
            self.assertEqual(0, main(["pathway-interval-falsify", "--family-index", str(index),
                                      "--hypotheses", str(hypotheses), "--interval", "chr1:1-100000",
                                      "--out-dir", str(out)]))
            result = pd.read_csv(out / "hypothesis_interval_falsification.tsv", sep="\t", dtype=str).fillna("")
            clv = result[(result.level == "coarse") & (result.hypothesis_id == "hyp_001")].iloc[0]
            self.assertEqual("2", clv.interval_members)
            self.assertEqual("1", clv.interval_blocks)
            self.assertEqual("True", clv.tandem_inflated)
            # The conservative p-value uses the collapsed block count, so it is weaker.
            self.assertGreater(float(clv.conservative_p), float(clv.enrichment_p))

    def test_surviving_candidates_carry_falsifiers_and_no_score(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            survivors = pd.read_csv(out / "surviving_candidates.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(len(survivors))
            self.assertTrue(survivors["rank_contribution"].eq("none").all())
            self.assertTrue(survivors["falsifier"].ne("").all())
            self.assertNotIn("dry_score", survivors.columns)
            self.assertNotIn("score", survivors.columns)
            # G_CLV_B's second-best bitscore is 95% of its best: paralog-ambiguous.
            row = survivors.set_index("gene_id").loc["G_CLV_B"]
            self.assertEqual("paralog_ambiguous", row.paralog_ambiguity)
            self.assertIn("RBH", row.falsifier)

    def test_unassigned_interval_genes_are_reported_as_a_blind_spot(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            manifest = json.loads((out / "pathway_interval_falsify_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(2, manifest["interval_genes_without_family_assignment"])
            reasons = set(pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("interval_genes_without_family_assignment", reasons)

    def test_manifest_reports_reduction_and_compute_budget(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            manifest = json.loads((out / "pathway_interval_falsify_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(4, manifest["hypotheses_tested"])
            # hyp_004 is a zero-resolution hold, not a falsification.
            self.assertEqual(0, manifest["coarse_falsified"])
            self.assertEqual(1, manifest["coarse_family_unresolved_in_index"])
            self.assertTrue(manifest["candidate_reduction"] >= 1.0)
            self.assertTrue(manifest["human_reviewable"])
            self.assertIn("wall_seconds", manifest["compute"])
            self.assertIn("no sequence search", manifest["compute"]["note"])
            self.assertIn("hypothesis filter only", manifest["statistic"]["role"])

    def test_falsify_runs_from_a_candidate_table_without_coordinates(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index, hypotheses = self._index_and_hypotheses(root)
            candidates = root / "candidates.tsv"
            candidates.write_text("gene_id\nG_CLV_A\nG_CLV_B\nG_PIN_A\n", encoding="utf-8")
            out = root / "falsify"
            self.assertEqual(0, main(["pathway-interval-falsify", "--family-index", str(index),
                                      "--hypotheses", str(hypotheses), "--candidates", str(candidates),
                                      "--out-dir", str(out)]))
            manifest = json.loads((out / "pathway_interval_falsify_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(3, manifest["interval_genes"])
            self.assertIn("hyp_001", manifest["surviving_hypotheses"])
            self.assertNotIn("hyp_003", manifest["surviving_hypotheses"])

    def test_llm_pathway_hypotheses_feed_straight_in(self) -> None:
        """The #163 proposal table is consumed without reshaping."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index, _ = self._index_and_hypotheses(root)
            raw = root / "raw.json"
            raw.write_text(json.dumps({"hypotheses": [
                {"pathway_or_process": "CLAVATA-WUSCHEL", "gene_family": "CLV",
                 "reference_gene_name": "CLV1", "reference_species": "Arabidopsis thaliana",
                 "proposed_stable_id": "sp|P0CLV1|CLV1_ARATH", "confidence": "medium"},
                {"pathway_or_process": "flowering time", "gene_family": "FT",
                 "reference_gene_name": "FT", "reference_species": "Arabidopsis thaliana",
                 "proposed_stable_id": "sp|P0FT001|FT_ARATH", "confidence": "low"},
            ]}), encoding="utf-8")
            imported = root / "import"
            self.assertEqual(0, main(["llm-pathway-import", "--raw", str(raw), "--out-dir", str(imported), "--case-id", "demo"]))
            out = root / "falsify"
            self.assertEqual(0, main(["pathway-interval-falsify", "--family-index", str(index),
                                      "--hypotheses", str(imported / "llm_pathway_hypotheses.tsv"),
                                      "--interval", "chr1:1-2000000", "--out-dir", str(out)]))
            frame = pd.read_csv(out / "hypothesis_interval_falsification.tsv", sep="\t", dtype=str).fillna("")
            coarse = frame[frame.level == "coarse"].set_index("hypothesis_id")
            self.assertTrue(coarse.loc["hyp_001", "verdict"].startswith("surviving"))
            self.assertEqual("not_evaluated_family_unresolved_in_index", coarse.loc["hyp_002", "verdict"])


    def test_incomplete_rescue_provider_holds_absence_not_falsifies(self) -> None:
        """A failed / partial rescue provider never completes the sequence-level
        search: its queries stay held, because 'no hit from a broken search' is
        not absence (issue #355, absence/rescue contract)."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rescue = root / "rescue_failed.tsv"
            rescue.write_text("\n".join([
                "query_protein_id\thit_count\tmodel_class\tprovider_status\tchrom\tstart\tend",
                "PIN1\t0\tnot_evaluated\tfailed\tchr9\t\t",
            ]) + "\n", encoding="utf-8")
            frame, manifest = self._absent_from_interval(root, "--annotation-rescue", str(rescue))
            coarse = frame[frame.level == "coarse"].set_index("hypothesis_id")
            self.assertEqual("absent_annotation_only", coarse.loc["hyp_002", "verdict"])
            self.assertFalse(coarse.loc["hyp_002", "verdict"].startswith("falsified"))
            self.assertEqual(1, manifest["absence_verification"]["queries_not_evaluated"])
            self.assertEqual(0, manifest["absence_verification"]["queries_evaluated"])
            self.assertIn("completeness_rule", manifest["absence_verification"])

    def test_complete_rescue_provider_without_hit_still_falsifies(self) -> None:
        """A completed provider (model_class=unmapped, provider_status=complete)
        is a real search with no hit: that may falsify."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rescue = root / "rescue_complete.tsv"
            rescue.write_text("\n".join([
                "query_protein_id\thit_count\tmodel_class\tprovider_status\tchrom\tstart\tend",
                "PIN1\t0\tunmapped\tcomplete\tchr9\t\t",
            ]) + "\n", encoding="utf-8")
            frame, manifest = self._absent_from_interval(root, "--annotation-rescue", str(rescue))
            coarse = frame[frame.level == "coarse"].set_index("hypothesis_id")
            self.assertEqual("falsified_no_member_in_interval", coarse.loc["hyp_002", "verdict"])
            self.assertEqual(1, manifest["absence_verification"]["queries_evaluated"])


class RouteAIntervalMatchTests(unittest.TestCase):
    """Issue #355: verified bait -> interval family match wiring adapter."""

    PANEL_HEADER = ("case_id\tbait_id\thypothesis_id\tpathway_or_process\tgene_family\t"
                    "reference_gene_name\treference_species\texpected_role\tsequence_length\t"
                    "sequence_sha256\tsource_provider\tsource_url\tretrieved_at\tfunction_source_ids\t"
                    "rank_contribution\tpanel_status")

    def _index_and_hypotheses(self, root: Path) -> tuple[Path, Path]:
        index = FamilyIndexTests()._build(root) / "family_index.tsv"
        hypotheses = root / "hypotheses.tsv"
        hypotheses.write_text("\n".join([
            "hypothesis_id\tpathway_or_process\tgene_family\treference_gene_name\tproposed_stable_id",
            "hyp_001\tCLAVATA-WUSCHEL\tCLV\tCLV1\tsp|P0CLV1|CLV1_ARATH",
            "hyp_002\tauxin transport\tPIN\tPIN1\tsp|P0PIN1|PIN1_ARATH",
            "hyp_003\tflowering time\tFT\tFT\tsp|P0FT001|FT_ARATH",
        ]) + "\n", encoding="utf-8")
        return index, hypotheses

    def _panel(self, root: Path, verified: list[str]) -> Path:
        panel = root / "verified_bait_panel.tsv"
        rows = [self.PANEL_HEADER]
        for hypothesis_id in verified:
            bait = {"hyp_001": ("sp|P0CLV1|CLV1_ARATH", "CLV", "CLV1", "PMID:8402901"),
                    "hyp_002": ("sp|P0PIN1|PIN1_ARATH", "PIN", "PIN1", "PMID:1978859"),
                    "hyp_003": ("sp|P0FT001|FT_ARATH", "FT", "FT", "PMID:1234567")}[hypothesis_id]
            rows.append("\t".join(["case_fas", bait[0], hypothesis_id, "", bait[1], bait[2],
                                    "Arabidopsis thaliana", "", "96", "sha256" + hypothesis_id,
                                    "local_record_cache", "", "2026-01-01T00:00:00Z", bait[3], "none", "frozen"]))
        panel.write_text("\n".join(rows) + "\n", encoding="utf-8")
        return panel

    def _run(self, root: Path, *extra: str, verified: list[str] | None = None,
             interval: str = "chr1:1-2000000", rescue: Path | None = None) -> Path:
        index, hypotheses = self._index_and_hypotheses(root)
        panel = self._panel(root, verified if verified is not None else ["hyp_001"])
        out = root / f"match{len(list(root.glob('match*')))}"
        args = ["route-a-interval-match", "--panel", str(panel),
                "--hypotheses", str(hypotheses), "--out-dir", str(out),
                "--family-index", str(index), "--interval", interval, *extra]
        if rescue:
            args += ["--annotation-rescue", str(rescue)]
        self.assertEqual(0, main(args))
        return out

    def _status(self, out: Path) -> pd.DataFrame:
        return pd.read_csv(out / "route_a_interval_family_match.tsv", sep="\t", dtype=str).fillna("") \
            .set_index("hypothesis_id")

    def test_verified_bait_forward_matches_to_concrete_interval_genes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = self._run(root)
            status = self._status(out)
            self.assertEqual("concrete_interval_gene_matched", status.loc["hyp_001", "family_status"])
            self.assertEqual("full", status.loc["hyp_001", "match_level"])
            self.assertIn("G_CLV_A", status.loc["hyp_001", "interval_gene_ids"].split(";"))
            self.assertEqual("True", status.loc["hyp_001", "evaluable"])
            self.assertEqual("True", status.loc["hyp_001", "stage_completed"])
            self.assertEqual("none", status.loc["hyp_001", "rank_contribution"])
            self.assertEqual("PMID:8402901", status.loc["hyp_001", "source_ids"])
            # only the verified target was forwarded to falsification
            forwarded = pd.read_csv(out / "verified_target_hypotheses.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(["hyp_001"], forwarded["hypothesis_id"].tolist())
            self.assertTrue((out / "hypothesis_interval_falsification.tsv").exists())
            manifest = json.loads((out / "route_a_interval_family_match_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(manifest["evaluated"])
            self.assertIsNone(manifest["not_evaluated_reason"])

    def test_unverified_hypotheses_are_not_forwarded_to_interval_matching(self) -> None:
        """LLM-only / pending hypotheses never reach the interval match lane."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = self._run(root, verified=["hyp_002"])
            status = self._status(out)
            self.assertNotIn("hyp_001", status.index)
            self.assertNotIn("hyp_003", status.index)
            self.assertEqual(["hyp_002"], pd.read_csv(out / "verified_target_hypotheses.tsv",
                                                      sep="\t", dtype=str)["hypothesis_id"].tolist())

    def test_zero_verified_baits_is_not_evaluated_not_a_miss(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index, hypotheses = self._index_and_hypotheses(root)
            panel = root / "empty_panel.tsv"
            panel.write_text(self.PANEL_HEADER + "\n", encoding="utf-8")
            out = root / "match"
            self.assertEqual(0, main(["route-a-interval-match", "--panel", str(panel),
                                      "--hypotheses", str(hypotheses), "--out-dir", str(out),
                                      "--family-index", str(index), "--interval", "chr1:1-2000000"]))
            manifest = json.loads((out / "route_a_interval_family_match_manifest.json").read_text(encoding="utf-8"))
            self.assertFalse(manifest["evaluated"])
            self.assertEqual("not_evaluated_no_verified_bait", manifest["not_evaluated_reason"])
            self.assertFalse((out / "hypothesis_interval_falsification.tsv").exists())
            frame = pd.read_csv(out / "route_a_interval_family_match.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(frame.empty)

    def test_missing_reference_genome_family_index_is_not_evaluated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index, hypotheses = self._index_and_hypotheses(root)
            panel = self._panel(root, ["hyp_001"])
            out = root / "match"
            self.assertEqual(0, main(["route-a-interval-match", "--panel", str(panel),
                                      "--hypotheses", str(hypotheses), "--out-dir", str(out),
                                      "--interval", "chr1:1-2000000"]))
            status = self._status(out)
            self.assertEqual("not_evaluated_missing_reference_genome_family_index",
                             status.loc["hyp_001", "family_status"])
            self.assertEqual("False", status.loc["hyp_001", "evaluable"])
            self.assertEqual("True", status.loc["hyp_001", "stage_completed"])
            self.assertFalse((out / "hypothesis_interval_falsification.tsv").exists())

    def test_not_evaluated_rerun_removes_stale_falsification_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index, hypotheses = self._index_and_hypotheses(root)
            panel = self._panel(root, ["hyp_001"])
            out = root / "match"
            evaluated = [
                "route-a-interval-match", "--panel", str(panel),
                "--hypotheses", str(hypotheses), "--out-dir", str(out),
                "--family-index", str(index), "--interval", "chr1:1-2000000",
            ]
            self.assertEqual(0, main(evaluated))
            stale_outputs = (
                "hypothesis_interval_falsification.tsv",
                "surviving_candidates.tsv",
                "review_required.tsv",
                "pathway_interval_falsify_manifest.json",
            )
            self.assertTrue(all((out / name).exists() for name in stale_outputs))

            self.assertEqual(0, main([
                "route-a-interval-match", "--panel", str(panel),
                "--hypotheses", str(hypotheses), "--out-dir", str(out),
                "--interval", "chr1:1-2000000",
            ]))
            manifest = json.loads((out / "route_a_interval_family_match_manifest.json").read_text(encoding="utf-8"))
            self.assertFalse(manifest["evaluated"])
            self.assertTrue(all(not (out / name).exists() for name in stale_outputs))

    def test_missing_interval_is_not_evaluated_interval_unresolved(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index, hypotheses = self._index_and_hypotheses(root)
            panel = self._panel(root, ["hyp_001"])
            out = root / "match"
            self.assertEqual(0, main(["route-a-interval-match", "--panel", str(panel),
                                      "--hypotheses", str(hypotheses), "--out-dir", str(out),
                                      "--family-index", str(index)]))
            manifest = json.loads((out / "route_a_interval_family_match_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("not_evaluated_interval_unresolved", manifest["not_evaluated_reason"])

    def test_annotation_absence_without_rescue_is_held(self) -> None:
        """PIN exists in the genome but not in a chr9 interval: absent_annotation_only."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = self._run(root, verified=["hyp_002"], interval="chr9:1-100000")
            status = self._status(out)
            self.assertEqual("hold_absent_annotation_only", status.loc["hyp_002", "family_status"])
            self.assertEqual("none", status.loc["hyp_002", "match_level"])
            self.assertEqual("True", status.loc["hyp_002", "evaluable"])
            self.assertEqual("not_run", status.loc["hyp_002", "annotation_rescue_status"])
            self.assertEqual("", status.loc["hyp_002", "interval_gene_ids"])

    def test_rescue_hit_survives_as_sequence_only_annotation_gap(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rescue = root / "rescue.tsv"
            rescue.write_text("\n".join([
                "query_protein_id\thit_count\tmodel_class\tprovider_status\tchrom\tstart\tend",
                "PIN1\t1\tunmapped\tcomplete\tchr9\t5000\t7000",
            ]) + "\n", encoding="utf-8")
            out = self._run(root, verified=["hyp_002"], interval="chr9:1-100000", rescue=rescue)
            status = self._status(out)
            self.assertEqual("surviving_sequence_only_annotation_gap", status.loc["hyp_002", "family_status"])
            self.assertEqual("sequence_only", status.loc["hyp_002", "match_level"])
            self.assertEqual("complete", status.loc["hyp_002", "annotation_rescue_status"])
            self.assertNotEqual("falsified", status.loc["hyp_002", "family_status"][:9])

    def test_completed_search_without_hit_falsifies(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rescue = root / "rescue.tsv"
            rescue.write_text("\n".join([
                "query_protein_id\thit_count\tmodel_class\tprovider_status\tchrom\tstart\tend",
                "PIN1\t0\tunmapped\tcomplete\tchr9\t\t",
            ]) + "\n", encoding="utf-8")
            out = self._run(root, verified=["hyp_002"], interval="chr9:1-100000", rescue=rescue)
            status = self._status(out)
            self.assertEqual("falsified_no_member_after_complete_sequence_search",
                             status.loc["hyp_002", "family_status"])
            self.assertEqual("complete", status.loc["hyp_002", "annotation_rescue_status"])

    def test_interval_scoped_family_index_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index, hypotheses = self._index_and_hypotheses(root)
            frame = pd.read_csv(index, sep="\t", dtype=str).fillna("")
            frame["index_scope"] = "interval_candidates"
            scoped = root / "interval_scoped_index.tsv"
            frame.to_csv(scoped, sep="\t", index=False, lineterminator="\n")
            panel = self._panel(root, ["hyp_001"])
            with self.assertRaises(ValueError):
                family_index_module.route_a_interval_family_match(
                    str(panel), str(hypotheses), str(root / "match"),
                    family_index_arg=str(scoped), intervals=["chr1:1-2000000"])

    def test_paralog_tie_is_family_present_copy_unresolved(self) -> None:
        """Only G_CLV_B (top2 ratio 0.95) in the interval: family present, but
        which copy is the ortholog is unresolved."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = self._run(root, verified=["hyp_001"], interval="chr1:1500000-1501000")
            status = self._status(out)
            self.assertEqual("family_present_copy_unresolved", status.loc["hyp_001", "family_status"])
            self.assertEqual("full", status.loc["hyp_001", "match_level"])
            self.assertEqual(["G_CLV_B"], status.loc["hyp_001", "interval_gene_ids"].split(";"))
            self.assertIn("paralog ties", status.loc["hyp_001", "notes"])

    def test_coarse_family_present_but_fine_ortholog_absent_is_copy_unresolved(self) -> None:
        """G_PIN_A's best hit is CLV2 (same CLV family bucket), not the verified
        CLV1 bait: the coarse bucket survives with a real interval member, but
        the fine level finds no CLV1 ortholog. Family presence must not be
        erased by the fine-level absence: family_present_copy_unresolved."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index, hypotheses = self._index_and_hypotheses(root)
            frame = pd.read_csv(index, sep="\t", dtype=str).fillna("")
            frame.loc[frame.gene_id == "G_PIN_A", ["best_subject", "best_subject_symbol", "best_subject_stem",
                                                   "family_subjects", "family_symbols", "family_stems"]] = [
                "sp|P0CLV2|CLV2_ARATH", "CLV2", "clv", "sp|P0CLV2|CLV2_ARATH", "CLV2", "clv"]
            frame.to_csv(index, sep="\t", index=False, lineterminator="\n")
            panel = self._panel(root, ["hyp_001"])
            out = root / "match"
            # chr1:900000-901000 holds only G_PIN_A: the CLV bucket has a
            # member there, but no CLV1 ortholog.
            self.assertEqual(0, main(["route-a-interval-match", "--panel", str(panel),
                                      "--hypotheses", str(hypotheses), "--out-dir", str(out),
                                      "--family-index", str(index), "--interval", "chr1:900000-901000"]))
            status = self._status(out)
            self.assertEqual("family_present_copy_unresolved", status.loc["hyp_001", "family_status"])
            self.assertEqual("full", status.loc["hyp_001", "match_level"])
            self.assertEqual(["G_PIN_A"], status.loc["hyp_001", "interval_gene_ids"].split(";"))
            self.assertIn("fine-level verified ortholog is absent", status.loc["hyp_001", "notes"])

    def test_interval_that_resolves_to_no_genes_is_not_evaluated_not_failed(self) -> None:
        """An interval given but with no genes in the index coordinates is an
        unresolved interval, not a stage failure (issue #355)."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index, hypotheses = self._index_and_hypotheses(root)
            panel = self._panel(root, ["hyp_001"])
            out = root / "match"
            self.assertEqual(0, main(["route-a-interval-match", "--panel", str(panel),
                                      "--hypotheses", str(hypotheses), "--out-dir", str(out),
                                      "--family-index", str(index), "--interval", "chr99:1-100000"]))
            manifest = json.loads((out / "route_a_interval_family_match_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("not_evaluated_interval_unresolved", manifest["not_evaluated_reason"])
            status = self._status(out)
            self.assertEqual("not_evaluated_interval_unresolved", status.loc["hyp_001", "family_status"])
            self.assertEqual("False", status.loc["hyp_001", "evaluable"])

    def test_no_verified_bait_priority_wins_over_a_missing_index(self) -> None:
        """Zero verified baits is the first not_evaluated reason even when the
        family index is also missing (deterministic degradation order)."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            index, hypotheses = self._index_and_hypotheses(root)
            panel = root / "empty_panel.tsv"
            panel.write_text(self.PANEL_HEADER + "\n", encoding="utf-8")
            out = root / "match"
            self.assertEqual(0, main(["route-a-interval-match", "--panel", str(panel),
                                      "--hypotheses", str(hypotheses), "--out-dir", str(out),
                                      "--interval", "chr1:1-2000000"]))
            manifest = json.loads((out / "route_a_interval_family_match_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("not_evaluated_no_verified_bait", manifest["not_evaluated_reason"])

    def test_rank_contribution_is_none_and_no_rank_artifact_is_written(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = self._run(root)
            status = self._status(out)
            self.assertTrue(status["rank_contribution"].eq("none").all())
            self.assertFalse((out / "ranked_candidates.tsv").exists())
            self.assertFalse((out / "rank").exists())
            survivors = pd.read_csv(out / "surviving_candidates.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(survivors["rank_contribution"].eq("none").all())
            self.assertNotIn("dry_score", survivors.columns)


if __name__ == "__main__":
    unittest.main()
