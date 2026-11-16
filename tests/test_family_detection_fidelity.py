from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.benchmark_metrics import aggregate_denominator_metrics
from candidate_qtl.family_index import (
    _members,
    family_index,
    pathway_interval_falsify,
)
from candidate_qtl.family_labels import ALIAS_PROVENANCE, families_match
from candidate_qtl.identifiers import canonical_accession


class AccessionCanonicalizationTests(unittest.TestCase):
    """Issue #411 / #412: Accession canonicalization handles pipe deflines,
    namespace prefixes, and raw accessions consistently."""

    def test_canonical_accession_normalization(self) -> None:
        self.assertEqual(canonical_accession("sp|P13114|CHSY_ARATH"), "P13114")
        self.assertEqual(canonical_accession("tr|A0A1I9LSP1|A0A1I9LSP1_ARATH"), "A0A1I9LSP1")
        self.assertEqual(canonical_accession("sp|P13114|"), "P13114")
        self.assertEqual(canonical_accession("UniProt:P13114"), "P13114")
        self.assertEqual(canonical_accession("uniprot:P13114"), "P13114")
        self.assertEqual(canonical_accession("TAIR:AT1G69120"), "AT1G69120")
        self.assertEqual(canonical_accession("RefSeq:NP_001323385.1"), "NP_001323385.1")
        self.assertEqual(canonical_accession("P13114"), "P13114")
        self.assertEqual(canonical_accession("AT1G69120"), "AT1G69120")
        self.assertEqual(canonical_accession(""), "")
        self.assertEqual(canonical_accession(None), "")


class SynonymAndCompoundNameTests(unittest.TestCase):
    """Issue #411: ANR (BAN) and LUP aliases match correctly with provenance."""

    def test_anr_and_lup_provenance(self) -> None:
        self.assertIn("anr", ALIAS_PROVENANCE)
        self.assertEqual(ALIAS_PROVENANCE["anr"]["accession"], "Q9SEV0")
        self.assertEqual(ALIAS_PROVENANCE["anr"]["database"], "UniProtKB")

        self.assertIn("lup", ALIAS_PROVENANCE)
        self.assertEqual(ALIAS_PROVENANCE["lup"]["accession"], "P93033")
        self.assertEqual(ALIAS_PROVENANCE["lup"]["database"], "UniProtKB")

    def test_anr_synonym_matching(self) -> None:
        self.assertTrue(families_match("BAN (ANR)", "BAN"))
        self.assertTrue(families_match("BAN (ANR)", "ANR"))
        self.assertTrue(families_match("BAN", "ANR"))
        self.assertTrue(families_match("BAN", "Anthocyanidin reductase"))
        self.assertTrue(families_match("BANYULS", "ANR"))

    def test_lup_synonym_matching(self) -> None:
        self.assertTrue(families_match("LUP1 (lupeol synthase)", "LUP1"))
        self.assertTrue(families_match("LUP1 (lupeol synthase)", "LUP"))
        self.assertTrue(families_match("LUP1", "lupeol synthase"))
        self.assertTrue(families_match("LUP2", "lupeol synthase"))


class MembersMatchingTests(unittest.TestCase):
    """Issue #411: _members matches normalized accession and symbol aliases."""

    def test_members_matches_uniprot_prefixed_and_pipe_accessions(self) -> None:
        index_records: list[dict[str, object]] = [
            {
                "gene_id": "Psat.cameor.v2.3g167300",
                "assignment_status": "assigned",
                "best_subject": "sp|P13114|CHSY_ARATH",
                "best_subject_symbol": "CHS",
                "best_subject_stem": "chs",
                "family_subjects": "sp|P13114|CHSY_ARATH;sp|P10290|CHSY_PETHY",
                "family_symbols": "CHS;CHSY",
                "family_stems": "chs",
            },
            {
                "gene_id": "Psat.cameor.v2.7g235150",
                "assignment_status": "assigned",
                "best_subject": "sp|Q9SEV0|BAN_ARATH",
                "best_subject_symbol": "BAN",
                "best_subject_stem": "anr",
                "family_subjects": "sp|Q9SEV0|BAN_ARATH",
                "family_symbols": "BAN",
                "family_stems": "anr",
            },
        ]

        # 1. Matching CHS via UniProt:P13114
        fine_chs = _members(index_records, accession="UniProt:P13114", gene_name="CHS", family_text="chalcone synthase", level="fine")
        self.assertEqual(len(fine_chs), 1)
        self.assertEqual(fine_chs[0]["gene_id"], "Psat.cameor.v2.3g167300")

        coarse_chs = _members(index_records, accession="UniProt:P13114", gene_name="CHS", family_text="chalcone synthase", level="coarse")
        self.assertEqual(len(coarse_chs), 1)
        self.assertEqual(coarse_chs[0]["gene_id"], "Psat.cameor.v2.3g167300")

        # 2. Matching ANR / BAN via compound symbol and alias
        fine_anr = _members(index_records, accession="UniProt:Q9SEV0", gene_name="BAN (ANR)", family_text="anthocyanidin reductase", level="fine")
        self.assertEqual(len(fine_anr), 1)
        self.assertEqual(fine_anr[0]["gene_id"], "Psat.cameor.v2.7g235150")

        coarse_anr = _members(index_records, accession="", gene_name="BAN (ANR)", family_text="anthocyanidin reductase", level="coarse")
        self.assertEqual(len(coarse_anr), 1)
        self.assertEqual(coarse_anr[0]["gene_id"], "Psat.cameor.v2.7g235150")


class FamilyIndexExclusionBreakdownTests(unittest.TestCase):
    """Issue #412: Family index manifest records dropped genes and exclusion breakdown."""

    def test_exclusion_breakdown_in_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            hits_path = tmp_path / "blast_hits.tsv"
            # Format: query_id, subject_id, identity, alignment_length, query_length, subject_length,
            #         query_start, query_end, subject_start, subject_end, evalue, bitscore, subject_title
            # min_identity=0.25 (25%), min_query_coverage=0.50 (50%), max_evalue=1e-5
            # Gene1: passes (identity 0.60, coverage 1.0, evalue 1e-50)
            # Gene2: fails coverage (0.336 < 0.50) -> BRC1 failure shape
            # Gene3: fails identity (0.15 < 0.25)
            hits_rows = [
                "gene1\tsp|P12345|ABC_ARATH\t60.0\t200\t200\t200\t1\t200\t1\t200\t1e-50\t250\tABC GN=ABC",
                "gene2\tsp|Q54321|TCP_ARATH\t55.0\t139\t414\t400\t1\t139\t1\t139\t1e-25\t140\tTCP GN=TCP",
                "gene3\tsp|Q99999|XYZ_ARATH\t15.0\t150\t200\t200\t1\t150\t1\t150\t1e-10\t80\tXYZ GN=XYZ",
            ]
            hits_path.write_text("\n".join(hits_rows) + "\n", encoding="utf-8")
            out_dir = tmp_path / "out"

            # Case A: without --genes (gene2 and gene3 dropped from index)
            res = family_index(blast_hits_arg=str(hits_path), out_dir_arg=str(out_dir))
            manifest = json.loads(Path(res["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual(manifest["gene_count"], 1)
            self.assertEqual(manifest["dropped_from_index_count"], 2)
            self.assertIn("exclusion_breakdown", manifest)
            breakdown = manifest["exclusion_breakdown"]
            self.assertEqual(breakdown["genes_with_hits"], 3)
            self.assertEqual(breakdown["genes_qualifying"], 1)
            self.assertEqual(breakdown["genes_no_qualifying_hits"], 2)
            self.assertEqual(breakdown["genes_below_min_query_coverage"], 1)
            self.assertEqual(breakdown["genes_below_min_identity"], 1)
            self.assertEqual(sorted(breakdown["dropped_from_index"]), ["gene2", "gene3"])

            # Case B: with --genes (gene2 and gene3 retained as dark_no_qualifying_hit)
            genes_path = tmp_path / "genes.tsv"
            genes_path.write_text("gene_id\tchrom\tstart\tend\n"
                                  "gene1\tchr1\t100\t200\n"
                                  "gene2\tchr1\t300\t400\n"
                                  "gene3\tchr1\t500\t600\n", encoding="utf-8")
            out_dir_genes = tmp_path / "out_genes"
            res_genes = family_index(blast_hits_arg=str(hits_path), genes_arg=str(genes_path), out_dir_arg=str(out_dir_genes))
            manifest_genes = json.loads(Path(res_genes["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual(manifest_genes["gene_count"], 3)
            self.assertEqual(manifest_genes["assigned_count"], 1)
            self.assertEqual(manifest_genes["dark_count"], 2)
            self.assertEqual(manifest_genes["dropped_from_index_count"], 0)


class UniverseCompletenessFalsificationTests(unittest.TestCase):
    """Issue #410: Partial universe degrades falsified_no_member_in_interval
    to absent_annotation_only."""

    def test_partial_vs_complete_universe(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)

            # Hypothesis: testing Family FOO (present on chr2 in genome, tested in interval chr1)
            hyp_path = tmp_path / "hypotheses.tsv"
            hyp_path.write_text("hypothesis_id\tproposed_stable_id\treference_gene_name\tgene_family\tpathway_or_process\n"
                                "H1\tsp|P99999|FOO_ARATH\tFOO\tfoo\tgrowth\n", encoding="utf-8")

            # Sequence rescue search evaluated with 0 hits
            rescue_path = tmp_path / "rescue.tsv"
            rescue_path.write_text("query_protein_id\thit_count\tloci\tmodel_class\tprovider_status\n"
                                   "sp|P99999|FOO_ARATH\t0\t\thit_table\tcompleted\n", encoding="utf-8")

            # Candidates table with 3 interval genes on chr1
            cand_path = tmp_path / "candidates.tsv"
            cand_path.write_text("gene_id\tchrom\tstart\tend\n"
                                 "cand1\tchr1\t100\t200\n"
                                 "cand2\tchr1\t300\t400\n"
                                 "cand3\tchr1\t500\t600\n", encoding="utf-8")

            # Partial Index: genome has genome_foo on chr2, but interval only has cand1 and cand2 in the index (cand3 missing)
            partial_index_dir = tmp_path / "partial_index"
            partial_index_dir.mkdir()
            partial_index_tsv = partial_index_dir / "family_index.tsv"
            partial_index_tsv.write_text(
                "gene_id\tprotein_id\tchrom\tstart\tend\tbest_subject\tbest_subject_symbol\tbest_subject_stem\t"
                "best_identity\tbest_query_coverage\tbest_evalue\tbest_bitscore\tsecond_bitscore\ttop2_bitscore_ratio\t"
                "symbol_source_rank\tsymbol_source_bitscore_ratio\tfamily_subjects\tfamily_symbols\tfamily_stems\t"
                "assignment_status\tdiscovery_bucket_id\tdiscovery_bucket_method\tevolutionary_family_id\t"
                "evolutionary_family_source\tevolutionary_family_status\tindex_scope\n"
                "genome_foo\tgf\tchr2\t100\t200\tsp|P99999|FOO_ARATH\tFOO\tfoo\t80.0\t1.0\t1e-50\t200\t0.0\t0.0\t1\t1.0\t"
                "sp|P99999|FOO_ARATH\tFOO\tfoo\tassigned\tfoo\tsymbol_stem\t\t\tunresolved\treference_genome\n"
                "cand1\tc1\tchr1\t100\t200\tsp|P11111|BAR_ARATH\tBAR\tbar\t80.0\t1.0\t1e-50\t200\t0.0\t0.0\t1\t1.0\t"
                "sp|P11111|BAR_ARATH\tBAR\tbar\tassigned\tbar\tsymbol_stem\t\t\tunresolved\treference_genome\n"
                "cand2\tc2\tchr1\t300\t400\tsp|P22222|BAZ_ARATH\tBAZ\tbaz\t80.0\t1.0\t1e-50\t200\t0.0\t0.0\t1\t1.0\t"
                "sp|P22222|BAZ_ARATH\tBAZ\tbaz\tassigned\tbaz\tsymbol_stem\t\t\tunresolved\treference_genome\n",
                encoding="utf-8",
            )
            partial_manifest = partial_index_dir / "family_index_manifest.json"
            partial_manifest.write_text(json.dumps({"index_scope": "reference_genome", "dropped_from_index_count": 1}), encoding="utf-8")

            # Run falsification on partial universe
            out_partial = tmp_path / "out_partial"
            res_partial = pathway_interval_falsify(
                family_index_arg=str(partial_index_tsv),
                hypotheses_arg=str(hyp_path),
                candidates_arg=str(cand_path),
                rescue_arg=str(rescue_path),
                out_dir_arg=str(out_partial),
            )
            self.assertFalse(res_partial["universe_complete"])
            fals_partial = pd.read_csv(out_partial / "hypothesis_interval_falsification.tsv", sep="\t")
            # Should degrade to absent_annotation_only instead of falsified_no_member_in_interval
            self.assertEqual(fals_partial.iloc[0]["verdict"], "absent_annotation_only")

            # Complete Index: cand1, cand2, cand3 all in index (cand3 as dark_no_qualifying_hit) + genome_foo on chr2
            complete_index_dir = tmp_path / "complete_index"
            complete_index_dir.mkdir()
            complete_index_tsv = complete_index_dir / "family_index.tsv"
            complete_index_tsv.write_text(
                "gene_id\tprotein_id\tchrom\tstart\tend\tbest_subject\tbest_subject_symbol\tbest_subject_stem\t"
                "best_identity\tbest_query_coverage\tbest_evalue\tbest_bitscore\tsecond_bitscore\ttop2_bitscore_ratio\t"
                "symbol_source_rank\tsymbol_source_bitscore_ratio\tfamily_subjects\tfamily_symbols\tfamily_stems\t"
                "assignment_status\tdiscovery_bucket_id\tdiscovery_bucket_method\tevolutionary_family_id\t"
                "evolutionary_family_source\tevolutionary_family_status\tindex_scope\n"
                "genome_foo\tgf\tchr2\t100\t200\tsp|P99999|FOO_ARATH\tFOO\tfoo\t80.0\t1.0\t1e-50\t200\t0.0\t0.0\t1\t1.0\t"
                "sp|P99999|FOO_ARATH\tFOO\tfoo\tassigned\tfoo\tsymbol_stem\t\t\tunresolved\treference_genome\n"
                "cand1\tc1\tchr1\t100\t200\tsp|P11111|BAR_ARATH\tBAR\tbar\t80.0\t1.0\t1e-50\t200\t0.0\t0.0\t1\t1.0\t"
                "sp|P11111|BAR_ARATH\tBAR\tbar\tassigned\tbar\tsymbol_stem\t\t\tunresolved\treference_genome\n"
                "cand2\tc2\tchr1\t300\t400\tsp|P22222|BAZ_ARATH\tBAZ\tbaz\t80.0\t1.0\t1e-50\t200\t0.0\t0.0\t1\t1.0\t"
                "sp|P22222|BAZ_ARATH\tBAZ\tbaz\tassigned\tbaz\tsymbol_stem\t\t\tunresolved\treference_genome\n"
                "cand3\tc3\tchr1\t500\t600\t\t\t\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t1\t1.0\t"
                "\t\t\tdark_no_qualifying_hit\t\t\t\t\tunresolved\treference_genome\n",
                encoding="utf-8",
            )
            complete_manifest = complete_index_dir / "family_index_manifest.json"
            complete_manifest.write_text(json.dumps({"index_scope": "reference_genome", "dropped_from_index_count": 0}), encoding="utf-8")

            # Run falsification on complete universe
            out_complete = tmp_path / "out_complete"
            res_complete = pathway_interval_falsify(
                family_index_arg=str(complete_index_tsv),
                hypotheses_arg=str(hyp_path),
                candidates_arg=str(cand_path),
                rescue_arg=str(rescue_path),
                out_dir_arg=str(out_complete),
            )
            self.assertTrue(res_complete["universe_complete"])
            fals_complete = pd.read_csv(out_complete / "hypothesis_interval_falsification.tsv", sep="\t")
            # Complete universe with sequence search completed confirms falsified_no_member_in_interval
            self.assertEqual(fals_complete.iloc[0]["verdict"], "falsified_no_member_in_interval")


class DenominatorEngineTests(unittest.TestCase):
    """Issue #350 / #412: not_evaluated_not_in_index is excluded from denominator."""

    def test_not_evaluated_not_in_index_excluded_from_denominator(self) -> None:
        rows = [
            {"state": "hit"},
            {"state": "hit"},
            {"state": "not_evaluated_not_in_index"},
            {"state": "not_evaluated"},
            {"state": "miss"},
        ]
        metrics = aggregate_denominator_metrics(rows)
        self.assertEqual(metrics["eligible"], 3)
        self.assertEqual(metrics["denominator"], 3)
        self.assertEqual(metrics["hit"], 2)
        self.assertEqual(metrics["miss"], 1)
        self.assertEqual(metrics["not_evaluated_not_in_index"], 1)
        self.assertEqual(metrics["not_evaluated"], 1)
        self.assertAlmostEqual(metrics["recall"], 2 / 3, places=6)


if __name__ == "__main__":
    unittest.main()
