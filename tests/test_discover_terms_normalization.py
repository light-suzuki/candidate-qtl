from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main


GENOME = ">chr1\n" + "ATG" * 300 + "\n"


def _gff(*genes: tuple[str, int, int, str]) -> str:
    lines = ["##gff-version 3"]
    for gid, start, end, note in genes:
        lines.append(f"chr1\tfixture\tgene\t{start}\t{end}\t.\t+\t.\tID={gid};Note={note}")
        lines.append(f"chr1\tfixture\tmRNA\t{start}\t{end}\t.\t+\t.\tID={gid}.1;Parent={gid}")
        lines.append(f"chr1\tfixture\tCDS\t{start}\t{end}\t.\t+\t0\tParent={gid}.1")
    return "\n".join(lines) + "\n"


class DiscoverTermsNormalizationTests(unittest.TestCase):
    """Issue #400 step 1: trait vocabulary normalization.

    Substring and shared-stem matches are additive to the exact whole-word
    boundary match, fire at low confidence, and must not create spurious hits
    from short shared prefixes or unrelated words.
    """

    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.fasta = self.root / "genome.fa"
        self.fasta.write_text(GENOME, encoding="utf-8")

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def _run(self, gff_text: str, phenotype: str) -> pd.DataFrame:
        gff = self.root / "genes.gff3"
        gff.write_text(gff_text, encoding="utf-8")
        out = self.root / f"out_{abs(hash(gff_text)) % 100000}"
        self.assertEqual(0, main(["discover", "--fasta", str(self.fasta), "--gff3", str(gff), "--interval", "chr1:1-900",
                                  "--phenotype", phenotype, "--out-dir", str(out), "--offline", "--no-auto-local-db"]))
        return pd.read_csv(out / "04_evidence.tsv", sep="\t").fillna("")

    def _match_row(self, evidence: pd.DataFrame, gene_id: str) -> pd.Series | None:
        rows = evidence[(evidence.gene_id == gene_id) & (evidence.evidence_code == "trait_function_match")]
        return rows.iloc[0] if len(rows) else None

    def test_substring_match_fires_low_confidence(self) -> None:
        # proanthocyanidin (phenotype) contains anthocyanidin (description word):
        # substring normalization must reach a match the boundary regex cannot.
        evidence = self._run(_gff(("AnrGene", 100, 300, "Anthocyanidin reductase BAN")),
                             "proanthocyanidin (condensed tannin) biosynthesis")
        row = self._match_row(evidence, "AnrGene")
        self.assertIsNotNone(row)
        self.assertEqual("low", row["confidence"])
        self.assertEqual("support", row["direction"])
        self.assertIn("proanthocyanidin", row["value"])
        self.assertIn("substring", row["description"])

    def test_stem_match_fires_low_confidence(self) -> None:
        # flavonoid / flavonone share the stem "flavon" (6 chars >= 6).
        evidence = self._run(_gff(("ChiGene", 100, 300, "Chalcone-flavonone isomerase family protein")),
                             "flavonoid and anthocyanin biosynthesis")
        row = self._match_row(evidence, "ChiGene")
        self.assertIsNotNone(row)
        self.assertEqual("low", row["confidence"])
        self.assertEqual("support", row["direction"])
        self.assertIn("flavonoid", row["value"])
        self.assertIn("stem", row["description"])

    def test_unrelated_words_do_not_fire(self) -> None:
        evidence = self._run(_gff(("SapGene", 100, 300, "saponin synthase")), "chalcone")
        self.assertIsNone(self._match_row(evidence, "SapGene"))

    def test_short_shared_prefix_does_not_fire(self) -> None:
        # colour/color share only "colo" (4 chars) < the stem minimum of 6.
        evidence = self._run(_gff(("ColGene", 100, 300, "color pigment protein")), "colour")
        self.assertIsNone(self._match_row(evidence, "ColGene"))

    def test_synth_stem_does_not_fire(self) -> None:
        # synthesis/synthase share only "synth" (5 chars) < the stem minimum of
        # 6; this regression locks out the pea_ans false-positive burst from the
        # first _MIN_STEM_LEN=5 implementation.
        evidence = self._run(_gff(("ChsGene", 100, 300, "chalcone synthase")), "anthocyanidin synthesis")
        self.assertIsNone(self._match_row(evidence, "ChsGene"))
        evidence = self._run(_gff(("SynGene", 100, 300, "flavonol synthase-like protein")),
                             "anthocyanidin synthesis")
        self.assertIsNone(self._match_row(evidence, "SynGene"))

    def test_substring_below_min_length_does_not_fire(self) -> None:
        # "colo" is contained in "color" but is below the 6-char substring guard.
        evidence = self._run(_gff(("ColGene", 100, 300, "color protein")), "colo")
        self.assertIsNone(self._match_row(evidence, "ColGene"))
        # A tiny description token ("red") must not trigger containment either.
        evidence = self._run(_gff(("RedGene", 100, 300, "red protein")), "proanthocyanidin")
        self.assertIsNone(self._match_row(evidence, "RedGene"))

    def test_exact_token_match_stays_medium(self) -> None:
        evidence = self._run(_gff(("ChsGene", 100, 300, "chalcone synthase")), "chalcone synthesis")
        row = self._match_row(evidence, "ChsGene")
        self.assertIsNotNone(row)
        self.assertEqual("medium", row["confidence"])

    def test_exact_phrase_match_stays_high(self) -> None:
        evidence = self._run(_gff(("PhraseGene", 100, 300, "meristem enlargement factor")), "meristem enlargement")
        row = self._match_row(evidence, "PhraseGene")
        self.assertIsNotNone(row)
        self.assertEqual("high", row["confidence"])

    def test_broad_anatomy_alone_stays_context(self) -> None:
        # Broad anatomy words are never normalized; only exact boundary matches
        # can produce trait_context, unchanged from before.
        evidence = self._run(_gff(("StemGene", 100, 300, "stem-associated protein")), "stem fasciation")
        rows = evidence[evidence.gene_id == "StemGene"]
        self.assertTrue((rows["evidence_code"] == "trait_context").any())
        self.assertTrue((rows["evidence_code"] == "trait_function_match").sum() == 0)

    def test_issue_negatives_bas1_and_cyp75a_do_not_fire(self) -> None:
        # Mirror the two most regex-adjacent negatives from issue #400.
        evidence = self._run(_gff(("BasGene", 100, 300, "Beta-amyrin synthase GN=bas1")),
                             "saponin (soyasaponin) content / bitterness")
        self.assertIsNone(self._match_row(evidence, "BasGene"))
        evidence = self._run(_gff(("F3hGene", 100, 300, "flavonoid 3',5'-hydroxylase")),
                             "pink vs purple/red flower colour")
        self.assertIsNone(self._match_row(evidence, "F3hGene"))

    def test_process_synthesis_does_not_fire_on_synthase_enzymes(self) -> None:
        # "synthesis" is a process word: it must not, on its own, fire substring/
        # stem normalization against unrelated synthase enzymes even though the
        # words share morphology. The discriminative token (anthocyanidin) is
        # exercised separately below.
        for description in ("HMG-CoA synthase", "ACC synthase GN=acc",
                            "chalcone synthase", "flavonol synthase-like protein"):
            evidence = self._run(_gff(("SynProcGene", 100, 300, description)),
                                 "anthocyanidin synthesis")
            self.assertIsNone(self._match_row(evidence, "SynProcGene"))

    def test_process_words_do_not_fire_substring_normalization(self) -> None:
        # Process/entity-category words carry no discriminative power on their
        # own, so a substring hit between the scored term and the description
        # (photosynthesis/biosynthesis, proteins, subfamily) must not fire.
        cases = [
            ("anthocyanidin synthesis", "Photosynthesis-related protein"),
            ("storage protein", "Seed proteins"),
            ("storage protein", "Glycoprotein-linked precursor"),
            ("R2R3 family", "MYB subfamily protein"),
        ]
        for index, (phenotype, description) in enumerate(cases):
            evidence = self._run(_gff((f"SubProcGene{index}", 100, 300, description)), phenotype)
            self.assertIsNone(self._match_row(evidence, f"SubProcGene{index}"))

    def test_process_words_do_not_fire_stem_normalization(self) -> None:
        # A shared stem between a process word and its morphological variants
        # (metabolism/metabolic, activity/activities, biosynthesis/biosynthetic,
        # synthesis/synthetic) must not fire normalization either.
        cases = [
            ("phenylpropanoid metabolism", "Metabolic enzyme"),
            ("enzyme activity", "Peroxidase activities"),
            ("chalcone biosynthesis", "Biosynthetic enzyme"),
            ("anthocyanidin synthesis", "synthetic enzyme"),
        ]
        for index, (phenotype, description) in enumerate(cases):
            evidence = self._run(_gff((f"StemProcGene{index}", 100, 300, description)), phenotype)
            self.assertIsNone(self._match_row(evidence, f"StemProcGene{index}"))

    def test_exact_match_of_process_word_still_fires(self) -> None:
        # The exclusion applies to normalization only: when the description
        # literally contains the process word, the exact token match still fires
        # at medium confidence.
        evidence = self._run(_gff(("ExactProcGene", 100, 300, "general synthesis regulator")),
                             "anthocyanidin synthesis")
        row = self._match_row(evidence, "ExactProcGene")
        self.assertIsNotNone(row)
        self.assertEqual("medium", row["confidence"])
        self.assertIn("synthesis", row["value"])

    def test_discriminative_token_with_process_word_still_fires(self) -> None:
        # A discriminative token co-occurring with a process word must keep its
        # existing matches: exact (glucosyltransferase) and stem (anthocyanin)
        # both fire despite "synthesis" being excluded.
        evidence = self._run(_gff(("DiscExactGene", 100, 300, "Anthocyanidin 3-O-glucosyltransferase")),
                             "anthocyanidin synthesis")
        row = self._match_row(evidence, "DiscExactGene")
        self.assertIsNotNone(row)
        self.assertEqual("medium", row["confidence"])
        self.assertIn("anthocyanidin", row["value"])
        evidence = self._run(_gff(("DiscStemGene", 100, 300, "anthocyanin synthase")),
                             "anthocyanidin synthesis")
        row = self._match_row(evidence, "DiscStemGene")
        self.assertIsNotNone(row)
        self.assertEqual("low", row["confidence"])
        self.assertIn("anthocyanidin", row["value"])


if __name__ == "__main__":
    unittest.main()
