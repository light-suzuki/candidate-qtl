"""Wiring tests for issue #230: annotation-rescue, sv-plan/multigenome-private,
promoter-review, and A/B dual-route disagreement deterministically feeding the
locus-hypothesis ledger without touching the gene rank.

These complement tests/test_hypothesis_ledger.py (unchanged behaviour) with the
four newly-wired paths, plus one end-to-end case-run check that the
annotation-rescue -> ledger wiring works through the orchestrate flow and that
--ranked stays byte-identical.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.hypothesis_ledger import HYPOTHESIS_COLUMNS, update_hypothesis_ledger


def _write(path: Path, rows: list[dict], columns: list[str] | None = None) -> Path:
    pd.DataFrame(rows, columns=columns).fillna("").to_csv(path, sep="\t", index=False, lineterminator="\n")
    return path


class HypothesisLedgerWiringTests(unittest.TestCase):
    def test_annotation_rescue_discordant_rows_seed_annotation_missing(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            rescue_rows = [
                {"query_protein_id": "Q1", "projection_id": "proj1", "model_class": "reference_rescue",
                 "chrom": "chr1", "start": "10", "end": "20", "native_gene_ids": ""},
                {"query_protein_id": "Q2", "projection_id": "proj2", "model_class": "native_discordant",
                 "chrom": "chr1", "start": "30", "end": "40", "native_gene_ids": "GeneA"},
                {"query_protein_id": "Q3", "projection_id": "proj3", "model_class": "native",
                 "chrom": "chr1", "start": "50", "end": "60", "native_gene_ids": "GeneB"},
            ]
            result = update_hypothesis_ledger(
                str(root / "out"), case_id="caseR",
                annotation_rescue_arg=str(_write(root / "rescue.tsv", rescue_rows)),
            )
            ledger = pd.read_csv(result["locus_hypotheses"], sep="\t", dtype=str).fillna("")
            self.assertEqual(2, len(ledger))
            self.assertEqual({"annotation_missing"}, set(ledger.hypothesis_type))
            by_target = ledger.set_index("target_gene_ids")
            self.assertIn("GeneA", by_target.index)
            discordant = by_target.loc["GeneA"]
            self.assertTrue(str(discordant.source_id).startswith("annotation_rescue:Q2:proj2"))
            self.assertEqual("workflow", discordant.proposal_origin)
            reference_rescue = ledger[ledger.target_gene_ids.eq("")].iloc[0]
            self.assertEqual("proj1", reference_rescue.interval_or_element_id)

    def test_sv_events_alone_seed_sv_cnv_pav_and_resolve_evidence_links(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            sv_rows = [
                {"call_id": "SV1", "gene_id": "GeneA"},
                {"call_id": "SV2", "gene_id": "GeneA"},
                {"call_id": "SV3", "gene_id": "GeneB"},
                {"call_id": "SV4", "gene_id": ""},
            ]
            result = update_hypothesis_ledger(
                str(root / "out"), case_id="caseS", sv_events_arg=str(_write(root / "sv.tsv", sv_rows)),
            )
            ledger = pd.read_csv(result["locus_hypotheses"], sep="\t", dtype=str).fillna("")
            self.assertEqual({"GeneA", "GeneB"}, set(ledger.member_gene_ids))
            self.assertTrue(ledger.hypothesis_type.eq("sv_cnv_pav").all())
            links = pd.read_csv(result["hypothesis_evidence_links"], sep="\t", dtype=str).fillna("")
            self.assertTrue(links.evidence_kind.eq("sv_event").all())
            self.assertTrue(links.link_status.eq("resolved").all())
            self.assertEqual({"SV1", "SV2"}, set(links[links.hypothesis_id.eq(ledger.set_index("member_gene_ids").loc["GeneA", "hypothesis_id"])].evidence_id))

    def test_dual_route_discordance_maps_linked_qtl_and_unknown_locus_mechanism(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            comparison_rows = [
                # Route A named a family the interval listing does not contain: the row
                # itself carries the interval-absence evidence, so linked_qtl is earned.
                {"comparison_id": "cmp_000", "gene_family": "FamX", "convergence_status": "phenotype_only_supported",
                 "shared_interval_genes": "", "route_a_anchor_status": "unanchored"},
                {"comparison_id": "cmp_001", "gene_family": "FamY", "convergence_status": "interval_only_supported",
                 "shared_interval_genes": "GeneC;GeneD", "route_a_anchor_status": ""},
                # Route A is anchored in the interval; Route B simply never reached it.
                # That is not evidence of a separate linked locus.
                {"comparison_id": "cmp_002", "gene_family": "FamW", "convergence_status": "phenotype_only_supported",
                 "shared_interval_genes": "", "route_a_anchor_status": "interval_anchored"},
                {"comparison_id": "cmp_003", "gene_family": "FamZ", "convergence_status": "mechanistic_convergence",
                 "shared_interval_genes": "GeneE", "route_a_anchor_status": "interval_anchored"},
            ]
            result = update_hypothesis_ledger(
                str(root / "out"), case_id="caseD",
                dual_route_comparison_arg=str(_write(root / "cmp.tsv", comparison_rows)),
            )
            ledger = pd.read_csv(result["locus_hypotheses"], sep="\t", dtype=str).fillna("")
            self.assertEqual(3, len(ledger))
            by_element = ledger.set_index("interval_or_element_id")
            self.assertEqual("linked_qtl", by_element.loc["FamX", "hypothesis_type"])
            self.assertEqual("", by_element.loc["FamX", "member_gene_ids"])
            self.assertEqual("unknown_locus_mechanism", by_element.loc["FamY", "hypothesis_type"])
            self.assertEqual({"GeneC", "GeneD"}, set(by_element.loc["FamY", "member_gene_ids"].split(";")))
            # An anchored phenotype-only row escalates, but must not claim a location.
            self.assertEqual("unknown_locus_mechanism", by_element.loc["FamW", "hypothesis_type"])
            self.assertIn("neither confirmed nor refuted", by_element.loc["FamW", "proposed_mechanism"])
            self.assertTrue(ledger.source_id.str.startswith("dual_route_compare:").all())

    def test_new_escalation_inputs_never_touch_ranked_gene_file(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            ranked = root / "ranked.tsv"
            ranked.write_text("gene_id\tdry_rank\nGeneA\t1\nGeneB\t2\n", encoding="utf-8")
            original = ranked.read_bytes()
            sv_rows = [{"call_id": "SV1", "gene_id": "GeneA"}]
            rescue_rows = [{"query_protein_id": "Q1", "projection_id": "p1", "model_class": "ambiguous_overlap",
                            "chrom": "chr1", "start": "1", "end": "10", "native_gene_ids": "GeneA;GeneB"}]
            comparison_rows = [{"comparison_id": "cmp_000", "gene_family": "FamX",
                               "convergence_status": "phenotype_only_supported", "shared_interval_genes": "",
                               "route_a_anchor_status": "unanchored"}]
            result = update_hypothesis_ledger(
                str(root / "out"), case_id="caseRank", ranked_arg=str(ranked),
                sv_events_arg=str(_write(root / "sv.tsv", sv_rows)),
                annotation_rescue_arg=str(_write(root / "rescue.tsv", rescue_rows)),
                dual_route_comparison_arg=str(_write(root / "cmp.tsv", comparison_rows)),
            )
            self.assertEqual(original, ranked.read_bytes())
            ledger = pd.read_csv(result["locus_hypotheses"], sep="\t", dtype=str).fillna("")
            self.assertEqual({"sv_cnv_pav", "annotation_missing", "linked_qtl"}, set(ledger.hypothesis_type))

    def test_deterministic_rerun_produces_identical_hypothesis_ids(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            sv_rows = [{"call_id": "SV1", "gene_id": "GeneA"}]
            sv_path = _write(root / "sv.tsv", sv_rows)
            first = update_hypothesis_ledger(str(root / "first"), case_id="caseX", sv_events_arg=str(sv_path))
            second = update_hypothesis_ledger(str(root / "second"), case_id="caseX", sv_events_arg=str(sv_path))
            first_ids = set(pd.read_csv(first["locus_hypotheses"], sep="\t", dtype=str)["hypothesis_id"])
            second_ids = set(pd.read_csv(second["locus_hypotheses"], sep="\t", dtype=str)["hypothesis_id"])
            self.assertEqual(first_ids, second_ids)

    def test_manual_hypotheses_still_work_without_any_new_input(self) -> None:
        # Regression guard: passing none of the new optional args must behave
        # exactly as before (annotation_rescue/dual_route_comparison default to
        # empty frames and contribute nothing).
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            row = {column: "" for column in HYPOTHESIS_COLUMNS}
            row.update({
                "hypothesis_id": "H1", "case_id": "caseM", "hypothesis_type": "single_gene_coding",
                "member_gene_ids": "GeneA", "proposed_mechanism": "test", "status": "proposed",
                "created_by": "manual:researcher", "proposal_origin": "manual", "review_status": "reviewed",
                "source_id": "local:proposal:H1",
            })
            result = update_hypothesis_ledger(
                str(root / "out"), hypotheses_arg=str(_write(root / "hyp.tsv", [row], HYPOTHESIS_COLUMNS))
            )
            ledger = pd.read_csv(result["locus_hypotheses"], sep="\t", dtype=str).fillna("")
            self.assertEqual(1, len(ledger))
            self.assertEqual("single_gene_coding", ledger.iloc[0].hypothesis_type)


class CaseRunWiringTests(unittest.TestCase):
    def _annotation_rescue_case(self, root: Path) -> tuple[Path, Path]:
        fasta = root / "genome.fa"
        fasta.write_text(">chr1\n" + "ACGT" * 100 + "\n", encoding="utf-8")
        candidates = root / "candidates.tsv"
        candidates.write_text(
            "gene_id\tchrom\tstart\tend\tstrand\tannotation_quality\n"
            "NativeGood\tchr1\t50\t130\t+\tcomplete_orf\n",
            encoding="utf-8",
        )
        queries = root / "queries.faa"
        queries.write_text(">Q_RESCUE\n" + "M" * 20 + "\n", encoding="utf-8")
        metadata = root / "query_metadata.tsv"
        metadata.write_text(
            "query_protein_id\tspecies\tassembly_id\tannotation_release\tsource_id\tquery_role\tselection_method\treview_status\n"
            "Q_RESCUE\tPisum sativum\tZW6\tr1\tUniProt:Q_RESCUE\tfunctional_homolog\tautomated_sequence_backed_search\treviewed\n",
            encoding="utf-8",
        )
        projection = root / "projection.gff3"
        projection.write_text(
            "##gff-version 3\n"
            "chr1__1_400\tminiprot\tmRNA\t145\t205\t.\t+\t.\tID=MP2;Rank=1;Identity=0.8;Positive=0.9;"
            "Frameshift=0;StopCodon=0;Target=Q_RESCUE 1 20\n"
            "chr1__1_400\tminiprot\tCDS\t145\t205\t.\t+\t0\tParent=MP2\n",
            encoding="utf-8",
        )
        case = root / "case.yaml"
        resources = root / "resources.yaml"
        case.write_text(
            "case_id: rescue_case\n"
            "trait:\n  phenotype_direction: increase\n"
            "source:\n  source_id: PMID:fixture\n  target_interval: chr1:1-400\n  literature_search_status: complete\n"
            "target:\n  fasta_resource: fasta\n"
            "inputs:\n"
            "  candidates_resource: candidates\n"
            "  annotation_rescue_query_proteins_resource: queries\n"
            "  annotation_rescue_query_metadata_resource: metadata\n"
            "  annotation_rescue_projection_gff_resource: projection\n"
            "annotation_rescue:\n  projection_complete: true\n"
            "assay:\n  top_n: 1\n",
            encoding="utf-8",
        )
        resources.write_text(
            f"resources:\n"
            f"  candidates: {{path: {candidates}}}\n"
            f"  fasta: {{path: {fasta}}}\n"
            f"  queries: {{path: {queries}}}\n"
            f"  metadata: {{path: {metadata}}}\n"
            f"  projection: {{path: {projection}}}\n",
            encoding="utf-8",
        )
        return case, resources

    def test_annotation_rescue_reaches_locus_hypotheses_through_case_run(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case, resources = self._annotation_rescue_case(root)
            out = root / "out"
            self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources), "--out-dir", str(out)]))
            ranked_before = (out / "rank" / "ranked_candidates.tsv").read_bytes()
            hypotheses = pd.read_csv(out / "hypotheses" / "locus_hypotheses.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("annotation_missing", set(hypotheses.hypothesis_type))
            rescue_rows = hypotheses[hypotheses.hypothesis_type.eq("annotation_missing")]
            self.assertTrue(rescue_rows.source_id.str.startswith("annotation_rescue:").all())
            # --ranked must stay byte-identical: rerun (resumed) and re-derive.
            self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources), "--out-dir", str(out), "--resume"]))
            self.assertEqual(ranked_before, (out / "rank" / "ranked_candidates.tsv").read_bytes())

    def test_sv_plan_with_incomplete_inputs_skips_cleanly_with_recorded_reason(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = root / "candidates.tsv"
            candidates.write_text("gene_id\tchrom\tstart\tend\tstrand\nGeneA\tchr1\t10\t100\t+\n", encoding="utf-8")
            fasta = root / "genome.fa"
            fasta.write_text(">chr1\n" + "ACGT" * 50 + "\n", encoding="utf-8")
            query_fasta = root / "query_genome.fa"
            query_fasta.write_text(">chr1\n" + "ACGT" * 50 + "\n", encoding="utf-8")
            case = root / "case.yaml"
            resources = root / "resources.yaml"
            case.write_text(
                "case_id: sv_incomplete\n"
                "source:\n  source_id: PMID:fixture\n  target_interval: chr1:1-200\n  literature_search_status: complete\n"
                "target:\n  fasta_resource: fasta\n"
                "inputs:\n"
                "  candidates_resource: candidates\n"
                "  sv_plan_query_fasta_resource: query_fasta\n",
                encoding="utf-8",
            )
            resources.write_text(
                f"resources:\n  candidates: {{path: {candidates}}}\n  fasta: {{path: {fasta}}}\n  query_fasta: {{path: {query_fasta}}}\n",
                encoding="utf-8",
            )
            out = root / "out"
            self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources), "--out-dir", str(out)]))
            review = (out / "review_required.tsv").read_text(encoding="utf-8")
            self.assertIn("sv_plan_inputs_incomplete", review)
            self.assertFalse((out / "sv_plan").exists())


if __name__ == "__main__":
    unittest.main()
