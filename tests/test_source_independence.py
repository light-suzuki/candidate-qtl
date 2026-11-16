from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.hardening import classify_source, independence_group


class SourceClassificationTests(unittest.TestCase):
    def test_source_kinds(self) -> None:
        self.assertEqual(classify_source("PMID:12345"), "reference")
        self.assertEqual(classify_source("doi:10.1234/abc"), "reference")
        self.assertEqual(classify_source("GO:0006355"), "reference")
        self.assertEqual(classify_source("/data/blastp_hits.tsv"), "local_file")
        self.assertEqual(classify_source("proteins.fa"), "local_file")
        self.assertEqual(classify_source("some free text"), "label")
        self.assertEqual(classify_source(""), "none")

    def test_stable_accessions_are_citable(self) -> None:
        # Namespaced database ids.
        for source in ("UniProt:Q9LZS3", "TAIR:AT1G69120", "RefSeq:NP_001078807",
                       "EnsemblPlants:AT1G69120", "NCBIProtein:XP_003612345"):
            self.assertEqual(classify_source(source), "reference", source)
        # Bare accessions matched as a whole token.
        for source in ("Q9LZS3", "P93003", "AT1G69120", "NP_001078807"):
            self.assertEqual(classify_source(source), "reference", source)
        # Free-text labels are not citable, even if they mention a gene.
        for label in ("meristem maintenance", "CLV3 candidate", "likely ortholog"):
            self.assertEqual(classify_source(label), "label", label)

    def test_undefined_compact_tokens_are_not_citable(self) -> None:
        # #153: compact free labels that merely look like an accession must not
        # be graded citable. The old generic ``[A-Za-z]{2,4}_?[0-9]{4,}`` grammar
        # accepted these; the allowlist must reject them.
        for label in ("gene1234", "trait2026", "abc_12345", "GENE0001", "locus42_2026"):
            self.assertEqual(classify_source(label), "label", label)
        # A known namespace with a malformed local part is not citable either.
        for label in ("PMID:fixture", "UniProt:candidate", "TAIR:meristem", "GO:candidate"):
            self.assertEqual(classify_source(label), "label", label)
        # Real INSDC / PMC accessions stay citable under the tightened grammar.
        for source in ("AAR03725", "AF010190", "PMC5694491", "XP_003612345"):
            self.assertEqual(classify_source(source), "reference", source)

    def test_orthogroup_artifact_ids_are_not_citable(self) -> None:
        for artifact in ("HOG:E0991234", "OG0001234", "OG:0001234",
                         "eggnog:OG0001234", "orthodb:OG0001234", "xOG12345"):
            self.assertEqual(classify_source(artifact), "label", artifact)
        for source in ("PMID:12345", "GO:0006355", "UniProt:Q9LZS3"):
            self.assertEqual(classify_source(source), "reference", source)

    def test_independence_group_prefers_explicit(self) -> None:
        # Two different papers reusing one dataset can be tied together.
        self.assertEqual(independence_group("PMID:1", "SHARED_DATASET"), "SHARED_DATASET")
        self.assertEqual(independence_group("PMID:2", "SHARED_DATASET"), "SHARED_DATASET")
        # A bare label cannot anchor independence.
        self.assertEqual(independence_group("free text", ""), "")


class TopIndependenceTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        (self.root / "candidates.tsv").write_text("gene_id\nGeneA\n", encoding="utf-8")

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def _rank(self, evidence_text: str, header: str = "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\n") -> dict:
        evidence = self.root / "ev.tsv"; evidence.write_text(header + evidence_text, encoding="utf-8")
        out = self.root / f"rank_{abs(hash(evidence_text)) % 100000}"
        self.assertEqual(0, main(["evidence-rank", "--candidates", str(self.root / "candidates.tsv"), "--evidence", str(evidence), "--out-dir", str(out)]))
        return pd.read_csv(out / "ranked_candidates.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id").loc["GeneA"].to_dict()

    def test_same_study_position_and_expression_is_not_independent(self) -> None:
        row = self._rank(
            "GeneA\tpeak_overlap\tsupport\thigh\t0\tqtl_position\tPMID:1\tpeak\n"
            "GeneA\texpression_target_specific\tsupport\thigh\t5\texpression\tPMID:1\texpr\n"
        )
        self.assertEqual(row["support_family_count"], "2")
        self.assertEqual(row["independent_source_count"], "1")
        self.assertEqual(row["top_eligibility"], "multiple_non_redundant_channels")
        self.assertNotEqual(row["dry_group"], "top")

    def test_independent_studies_position_and_expression_reach_top(self) -> None:
        row = self._rank(
            "GeneA\tpeak_overlap\tsupport\thigh\t0\tqtl_position\tPMID:1\tpeak\n"
            "GeneA\texpression_target_specific\tsupport\thigh\t5\texpression\tPMID:2\texpr\n"
        )
        self.assertEqual(row["independent_source_count"], "2")
        self.assertEqual(row["top_eligibility"], "independent_experimental_support")
        self.assertEqual(row["dry_group"], "top")

    def test_explicit_independence_group_ties_sources_together(self) -> None:
        # Different accessions but the same declared independence group => 1 group.
        row = self._rank(
            "GeneA\tpeak_overlap\tsupport\thigh\t0\tqtl_position\tPMID:1\tpeak\tSHARED\n"
            "GeneA\texpression_target_specific\tsupport\thigh\t5\texpression\tPMID:2\texpr\tSHARED\n",
            header="gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\tindependence_group\n",
        )
        self.assertEqual(row["independent_source_count"], "1")
        self.assertEqual(row["top_eligibility"], "multiple_non_redundant_channels")

    def test_label_sources_do_not_prove_independence(self) -> None:
        row = self._rank(
            "GeneA\tpeak_overlap\tsupport\thigh\t0\tqtl_position\tinterval_note\tpeak\n"
            "GeneA\texpression_target_specific\tsupport\thigh\t5\texpression\tmatrix_note\texpr\n"
        )
        # Two channels but neither source is citable -> not independent.
        self.assertEqual(row["support_family_count"], "2")
        self.assertEqual(row["independent_source_count"], "0")
        self.assertEqual(row["top_eligibility"], "multiple_non_redundant_channels")

    def test_accession_like_labels_do_not_prove_independence(self) -> None:
        # #153 regression: two support channels whose source_ids merely look like
        # accessions (compact free tokens) must not be counted as two independent
        # citable sources, so they cannot reach independent_experimental_support
        # or Top. Before the fix these were graded ``reference``.
        row = self._rank(
            "GeneA\tpeak_overlap\tsupport\thigh\t0\tqtl_position\tgene1234\tpeak\n"
            "GeneA\texpression_target_specific\tsupport\thigh\t5\texpression\tabc_12345\texpr\n"
        )
        self.assertEqual(row["independent_source_count"], "0")
        self.assertNotEqual(row["top_eligibility"], "independent_experimental_support")
        self.assertNotEqual(row["dry_group"], "top")

    def test_orthogroup_artifact_ids_do_not_prove_independence(self) -> None:
        row = self._rank(
            "GeneA\tpeak_overlap\tsupport\thigh\t0\tqtl_position\tHOG:E0991234\tpeak\n"
            "GeneA\texpression_target_specific\tsupport\thigh\t5\texpression\tOG0001234\texpr\n"
        )
        self.assertEqual(row["independent_source_count"], "0")
        self.assertNotEqual(row["top_eligibility"], "independent_experimental_support")
        self.assertNotEqual(row["dry_group"], "top")

    def test_same_chain_homolog_accessions_are_one_independent_source(self) -> None:
        # Two UniProt hits from one BLAST run are now both citable accessions, but
        # they come from a single comparative-annotation chain, so they must count
        # as one independent source -- not inflate independent_source_count to 2.
        row = self._rank(
            "GeneA\tcurated_homolog_match\tsupport\thigh\t1\thomology\tUniProt:Q9LZS3\thit1\n"
            "GeneA\tcurated_homolog_match\tsupport\thigh\t1\thomology\tUniProt:P93003\thit2\n"
        )
        self.assertEqual(row["independent_source_count"], "1")
        self.assertNotEqual(row["top_eligibility"], "independent_experimental_support")

    def test_same_chain_homolog_and_trait_keyword_are_one_independent_source(self) -> None:
        # One blastp run + the annotation of its own best hit produces both a
        # homology row (subject accession) and a trait_function_match row (the
        # subject annotation file). Both are comparative_function evidence from a
        # single chain, so they must count as one independent source. Before the
        # fix the annotation file anchored a second group, which made a single
        # BLAST run read as ``independent_experimental_support``.
        row = self._rank(
            "GeneA\tcurated_homolog_match\tsupport\thigh\t72\tblast\tUniProt:Q9LZS3\tbest blastp hit\n"
            "GeneA\ttrait_function_match\tsupport\thigh\tstarch\tgff_annotation_keyword\t/db/subject_annotation.gff3\tkeyword in that hit's description\n"
        )
        self.assertEqual(row["independent_source_count"], "1")
        self.assertNotEqual(row["top_eligibility"], "independent_experimental_support")

    def test_same_chain_comparative_evidence_plus_position_does_not_reach_top(self) -> None:
        # Regression for the promotion path: the same single-chain pair above,
        # plus positional evidence, cleared --min-top-score and was promoted to
        # Top on what is really one comparative chain and a coordinate.
        row = self._rank(
            "GeneA\tcurated_homolog_match\tsupport\thigh\t72\tblast\tUniProt:Q9LZS3\tbest blastp hit\n"
            "GeneA\ttrait_function_match\tsupport\thigh\tstarch\tgff_annotation_keyword\t/db/subject_annotation.gff3\tkeyword in that hit's description\n"
            "GeneA\tpeak_overlap\tsupport\thigh\t0\tqtl_position\tchr1:1500\tgene overlaps the peak\n"
        )
        self.assertEqual(row["independent_source_count"], "1")
        self.assertNotEqual(row["dry_group"], "top")

    def test_citable_trait_function_still_anchors_its_own_source(self) -> None:
        # The collapse is scoped to uncited pipeline artifacts. A curated
        # trait-function claim with a real publication behind it is a separate
        # source from a domain scan, so #189's promotion path stays open.
        row = self._rank(
            "GeneA\ttrait_function_match\tsupport\thigh\tstarch\tliterature\tPMID:12345\tcurated trait function\n"
            "GeneA\tdomain_supported\tsupport\thigh\tPF02806\thmmer\tPfam:PF02806\talpha-amylase C domain\n"
        )
        self.assertEqual(row["independent_source_count"], "2")
        self.assertEqual(row["top_eligibility"], "independent_experimental_support")

    def test_uncitable_strong_code_does_not_bypass_top_gate(self) -> None:
        row = self._rank(
            "GeneA\tortholog_model_genetic\tsupport\thigh\tmodel phenotype\tmanual\tortholog_LOF_table\tuncited hypothesis\n"
        )
        self.assertEqual(row["independent_source_count"], "0")
        self.assertEqual(row["top_eligibility"], "insufficient_independent_evidence")
        self.assertNotEqual(row["dry_group"], "top")


if __name__ == "__main__":
    unittest.main()
