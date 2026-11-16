"""Callability-aware multigenome private-variant funnel (issue #164)."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.multigenome_private import normalized_event_key

GFF3 = (
    "##gff-version 3\n"
    "chr1\tfixture\tgene\t1000\t2000\t.\t+\t.\tID=GeneA\n"
    "chr1\tfixture\tmRNA\t1000\t2000\t.\t+\t.\tID=GeneA.1;Parent=GeneA\n"
    "chr1\tfixture\texon\t1000\t2000\t.\t+\t.\tParent=GeneA.1\n"
    "chr1\tfixture\tCDS\t1100\t1900\t.\t+\t0\tParent=GeneA.1\n"
    "chr1\tfixture\tgene\t3000\t4000\t.\t+\t.\tID=GeneB\n"
    "chr1\tfixture\tmRNA\t3000\t4000\t.\t+\t.\tID=GeneB.1;Parent=GeneB\n"
    "chr1\tfixture\texon\t3000\t4000\t.\t+\t.\tParent=GeneB.1\n"
    "chr1\tfixture\tCDS\t3100\t3900\t.\t+\t0\tParent=GeneB.1\n"
    "chr1\tfixture\tgene\t5000\t6000\t.\t+\t.\tID=GeneC\n"
    "chr1\tfixture\tmRNA\t5000\t6000\t.\t+\t.\tID=GeneC.1;Parent=GeneC\n"
    "chr1\tfixture\texon\t5000\t6000\t.\t+\t.\tParent=GeneC.1\n"
    "chr1\tfixture\tCDS\t5100\t5900\t.\t+\t0\tParent=GeneC.1\n"
    "chr1\tfixture\tgene\t7000\t8000\t.\t+\t.\tID=GeneD\n"
    "chr1\tfixture\tmRNA\t7000\t8000\t.\t+\t.\tID=GeneD.1;Parent=GeneD\n"
    "chr1\tfixture\tCDS\t7100\t7900\t.\t+\t0\tParent=GeneD.1\n"
)

VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##bcftools_normCommand=norm -f reference.fa -m -any\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
)


def _vcf(path: Path, records: list[tuple[str, int, str, str, str]], normalized: bool = True) -> None:
    header = VCF_HEADER if normalized else VCF_HEADER.replace("##bcftools_normCommand=norm -f reference.fa -m -any\n", "")
    lines = [header.rstrip("\n")]
    for chrom, pos, ref, alt, genotype in records:
        lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT\t{genotype}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _bed(path: Path, spans: list[tuple[str, int, int]]) -> None:
    # BED is 0-based half-open.
    path.write_text("".join(f"{chrom}\t{start - 1}\t{end}\n" for chrom, start, end in spans), encoding="utf-8")


class EventNormalizationTests(unittest.TestCase):
    def test_minimal_representation_matches_equivalent_spellings(self) -> None:
        self.assertEqual(("chr1", 100, "A", "T"), normalized_event_key("chr1", 100, "A", "T"))
        # Same 1 bp deletion written with and without a shared suffix.
        self.assertEqual(normalized_event_key("chr1", 100, "AGG", "AG"),
                         normalized_event_key("chr1", 100, "AG", "A"))
        # Shared prefix trimming shifts the position.
        self.assertEqual(("chr1", 101, "G", "T"), normalized_event_key("chr1", 100, "AG", "AT"))
        # A pure SNV is untouched.
        self.assertEqual(("chr1", 100, "C", "G"), normalized_event_key("chr1", 100, "c", "g"))


class MultigenomePrivateTests(unittest.TestCase):
    def _fixture(self, root: Path) -> dict[str, Path]:
        paths = {name: root / name for name in (
            "reference.gff3", "focal.vcf", "control_a.vcf", "control_b.vcf", "other.vcf",
            "control_a.bed", "control_b.bed", "projection.tsv", "reference_csq.tsv", "native_csq.tsv", "case.yaml")}
        paths["reference.gff3"].write_text(GFF3, encoding="utf-8")
        _vcf(paths["focal.vcf"], [
            ("chr1", 1200, "A", "T", "1/1"),   # GeneA CDS: private, both controls callable
            ("chr1", 3200, "C", "G", "1/1"),   # GeneB CDS: shared with control_a
            ("chr1", 5200, "G", "A", "1/1"),   # GeneC CDS: control_b uncallable -> private_possible
            ("chr1", 7050, "T", "C", "1/1"),   # GeneD non-coding: private, regulatory
        ])
        _vcf(paths["control_a.vcf"], [("chr1", 3200, "C", "G", "0/1")])
        # control_b carries the GeneB event as a homozygous reference call: a record
        # exists but the ALT allele is not present, so it is not a shared event.
        _vcf(paths["control_b.vcf"], [("chr1", 3200, "C", "G", "0/0")])
        _vcf(paths["other.vcf"], [("chr1", 1200, "A", "T", "1/1")])
        _bed(paths["control_a.bed"], [("chr1", 1000, 8000)])
        # control_b cannot call the GeneC region.
        _bed(paths["control_b.bed"], [("chr1", 1000, 5000), ("chr1", 5500, 8000)])
        paths["projection.tsv"].write_text(
            "gene_id\tprojection_state\tprojection_run_complete\n"
            "GeneA\thigh_confidence_one_to_one\ttrue\n"
            "GeneB\thigh_confidence_one_to_one\ttrue\n"
            "GeneC\thigh_confidence_one_to_one\ttrue\n"
            "GeneD\tunmapped\tfalse\n", encoding="utf-8")
        paths["reference_csq.tsv"].write_text(
            "chrom\tpos\tref\talt\tgene_id\tconsequence\n"
            "chr1\t1200\tA\tT\tGeneA\tmissense_variant\n"
            "chr1\t3200\tC\tG\tGeneB\tmissense_variant\n"
            "chr1\t5200\tG\tA\tGeneC\tstop_gained\n"
            "chr1\t7050\tT\tC\tGeneD\t5_prime_UTR_variant\n", encoding="utf-8")
        paths["native_csq.tsv"].write_text(
            "chrom\tpos\tref\talt\tgene_id\tconsequence\n"
            "chr1\t1200\tA\tT\tGeneA\tmissense_variant\n"
            "chr1\t5200\tG\tA\tGeneC\tsynonymous_variant\n", encoding="utf-8")
        paths["case.yaml"].write_text(
            "case_id: fixture_case\n"
            "source:\n"
            "  target_intervals: ['chr1:900-8100']\n"
            "private_variant_model:\n"
            "  focal_accession: focal_parent\n"
            "  phenotype_controls: [control_a, control_b]\n"
            "  comparison_accessions: [other]\n"
            "  event_normalization: left_normalized_ref_alt\n"
            "  require_all_controls_callable: true\n", encoding="utf-8")
        return paths

    def _run(self, root: Path, out: Path, extra: list[str] | None = None) -> dict[str, pd.DataFrame]:
        fixture = self._fixture(root)
        argv = ["multigenome-private", "--out-dir", str(out), "--case", str(fixture["case.yaml"]),
                "--focal-vcf", str(fixture["focal.vcf"]), "--gff3", str(fixture["reference.gff3"]),
                "--control-vcf", f"control_a={fixture['control_a.vcf']}",
                "--control-vcf", f"control_b={fixture['control_b.vcf']}",
                "--control-vcf", f"other={fixture['other.vcf']}",
                "--callable-bed", f"control_a={fixture['control_a.bed']}",
                "--callable-bed", f"control_b={fixture['control_b.bed']}",
                "--projection-ledger", str(fixture["projection.tsv"]),
                "--reference-consequence", str(fixture["reference_csq.tsv"]),
                "--native-consequence", str(fixture["native_csq.tsv"])] + (extra or [])
        self.assertEqual(0, main(argv))
        return {name: pd.read_csv(out / f"{name}.tsv", sep="\t", dtype=str).fillna("") for name in (
            "focal_event_sharing_ledger", "candidate_control_matrix", "control_callability_sensitivity",
            "event_class_summary", "gene_private_variant_ledger", "projection_and_callable_denominators",
            "working_candidate_panel", "review_required")}

    def test_event_classes_are_reproducible_per_control(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tables = self._run(root, root / "out")
            classes = dict(zip(tables["focal_event_sharing_ledger"]["event_key"],
                               tables["focal_event_sharing_ledger"]["event_class"]))
            self.assertEqual("private_confirmed", classes["chr1:1200:A>T"])
            self.assertEqual("shared_confirmed", classes["chr1:3200:C>G"])
            self.assertEqual("private_possible", classes["chr1:5200:G>A"])
            self.assertEqual("private_confirmed", classes["chr1:7050:T>C"])

    def test_missing_vcf_record_is_not_a_callable_reference_allele(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tables = self._run(root, root / "out")
            matrix = tables["candidate_control_matrix"]
            # control_b has no record for the GeneC event and cannot call the region.
            row = matrix[(matrix.event_key == "chr1:5200:G>A") & (matrix.control_accession == "control_b")].iloc[0]
            self.assertEqual("no_record", row.genotype_state)
            self.assertEqual("uncallable", row.callability_state)
            self.assertEqual("absent_uncallable", row.resolution)
            self.assertEqual("False", row.event_observed)
            # A record whose genotype is 0/0 inside a callable region is a real absence.
            reference_call = matrix[(matrix.event_key == "chr1:3200:C>G") & (matrix.control_accession == "control_b")].iloc[0]
            self.assertEqual("reference_or_other_allele", reference_call.genotype_state)
            self.assertEqual("callable", reference_call.callability_state)
            self.assertEqual("absent_callable", reference_call.resolution)

    def test_one_uncallable_required_control_blocks_private_confirmed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tables = self._run(root, root / "out")
            sensitivity = tables["control_callability_sensitivity"].set_index("control_accession")
            self.assertEqual("1", sensitivity.loc["control_b", "private_confirmed_blocked_by_callability"])
            self.assertEqual("0", sensitivity.loc["control_a", "private_confirmed_blocked_by_callability"])

    def test_comparison_accessions_do_not_change_phenotype_private_status(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tables = self._run(root, root / "out")
            ledger = tables["focal_event_sharing_ledger"].set_index("event_key")
            # `other` carries the GeneA event, but it is not a phenotype control.
            self.assertEqual("private_confirmed", ledger.loc["chr1:1200:A>T", "event_class"])
            self.assertEqual("", ledger.loc["chr1:1200:A>T", "shared_with"])
            matrix = tables["candidate_control_matrix"]
            other = matrix[(matrix.event_key == "chr1:1200:A>T") & (matrix.control_accession == "other")].iloc[0]
            self.assertEqual("comparison", other.control_role)
            self.assertEqual("shared", other.resolution)

    def test_per_locus_phenotype_controls_are_applied(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._fixture(root)
            # GeneC's locus only requires control_a, which *can* call the region.
            fixture["case.yaml"].write_text(
                "case_id: fixture_case\n"
                "private_variant_model:\n"
                "  focal_accession: focal_parent\n"
                "  phenotype_controls: [control_a, control_b]\n"
                "  loci:\n"
                "    - locus_id: locus_ab\n"
                "      interval: chr1:900-4500\n"
                "    - locus_id: locus_c\n"
                "      interval: chr1:4501-8100\n"
                "      phenotype_controls: [control_a]\n", encoding="utf-8")
            out = root / "out"
            self.assertEqual(0, main(["multigenome-private", "--out-dir", str(out), "--case", str(fixture["case.yaml"]),
                                      "--focal-vcf", str(fixture["focal.vcf"]), "--gff3", str(fixture["reference.gff3"]),
                                      "--control-vcf", f"control_a={fixture['control_a.vcf']}",
                                      "--control-vcf", f"control_b={fixture['control_b.vcf']}",
                                      "--callable-bed", f"control_a={fixture['control_a.bed']}",
                                      "--callable-bed", f"control_b={fixture['control_b.bed']}",
                                      "--projection-ledger", str(fixture["projection.tsv"])]))
            ledger = pd.read_csv(out / "focal_event_sharing_ledger.tsv", sep="\t", dtype=str).fillna("").set_index("event_key")
            self.assertEqual("locus_c", ledger.loc["chr1:5200:G>A", "locus_id"])
            self.assertEqual("control_a", ledger.loc["chr1:5200:G>A", "phenotype_controls"])
            self.assertEqual("private_confirmed", ledger.loc["chr1:5200:G>A", "event_class"])
            self.assertEqual("locus_ab", ledger.loc["chr1:1200:A>T", "locus_id"])
            self.assertEqual("control_a;control_b", ledger.loc["chr1:1200:A>T", "phenotype_controls"])
            # control_b is still reconciled for locus_c, as a comparison accession.
            matrix = pd.read_csv(out / "candidate_control_matrix.tsv", sep="\t", dtype=str).fillna("")
            locus_c = matrix[matrix.event_key == "chr1:5200:G>A"].set_index("control_accession")
            self.assertEqual("phenotype_control", locus_c.loc["control_a", "control_role"])
            self.assertEqual("comparison", locus_c.loc["control_b", "control_role"])

    def test_dual_frame_disagreement_is_reviewed_not_dropped(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tables = self._run(root, root / "out")
            ledger = tables["focal_event_sharing_ledger"].set_index("event_key")
            gene_c = ledger.loc["chr1:5200:G>A"]
            self.assertEqual("stop_gained", gene_c.reference_consequence)
            self.assertEqual("synonymous_variant", gene_c.native_consequence)
            self.assertEqual("disagree", gene_c.consequence_agreement)
            self.assertEqual("review_required", gene_c.review_status)
            reasons = set(tables["review_required"]["reason"])
            self.assertIn("dual_frame_consequence_disagreement", reasons)
            # The reference-only annotation stays visible instead of being discarded.
            self.assertEqual("reference_frame_only", ledger.loc["chr1:3200:C>G", "consequence_agreement"])

    def test_denominators_and_incomplete_projection_are_reported(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tables = self._run(root, root / "out")
            denominators = tables["projection_and_callable_denominators"].iloc[0]
            self.assertEqual("4", denominators.reference_genes_in_interval)
            self.assertEqual("3", denominators.high_confidence_one_to_one)
            # GeneD's projection run was incomplete, so its no-hit is not_evaluated.
            self.assertEqual("0", denominators.unmapped)
            self.assertEqual("1", denominators.projection_not_evaluated)
            self.assertTrue(int(denominators.intended_cds_bp) > 0)
            self.assertTrue(float(denominators.callable_cds_fraction) < 1.0)
            self.assertIn("projection", denominators.denominator_status)
            self.assertTrue(denominators.required_next_action)

    def test_private_possible_is_kept_and_multiple_co_primary_are_allowed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tables = self._run(root, root / "out")
            panel = tables["working_candidate_panel"].set_index("gene_id")
            self.assertIn("GeneC", panel.index)
            self.assertEqual("reserve", panel.loc["GeneC", "candidate_group"])
            self.assertIn("callability", panel.loc["GeneC", "required_next_action"])
            self.assertEqual("co_primary", panel.loc["GeneA", "candidate_group"])
            self.assertEqual("regulatory_sv", panel.loc["GeneD", "candidate_group"])
            self.assertTrue(panel["rank_contribution"].eq("none").all())
            self.assertTrue(set(panel["candidate_group"]) <= {
                "co_primary", "reserve", "annotation_rescue", "regulatory_sv", "evidence_insufficient"})

    def test_no_score_column_and_counts_are_not_summed_across_callsets(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tables = self._run(root, root / "out")
            for name, frame in tables.items():
                self.assertNotIn("dry_score", frame.columns, name)
                self.assertNotIn("score", frame.columns, name)
            self.assertTrue(tables["gene_private_variant_ledger"]["rank_contribution"].eq("none").all())
            # Every focal event appears once per declared accession, never merged.
            matrix = tables["candidate_control_matrix"]
            per_event = matrix.groupby("event_key")["control_accession"].nunique()
            self.assertTrue((per_event == 3).all())

    def test_missing_callability_keeps_events_at_private_possible(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._fixture(root)
            out = root / "out"
            self.assertEqual(0, main(["multigenome-private", "--out-dir", str(out), "--case", str(fixture["case.yaml"]),
                                      "--focal-vcf", str(fixture["focal.vcf"]), "--gff3", str(fixture["reference.gff3"]),
                                      "--control-vcf", f"control_a={fixture['control_a.vcf']}",
                                      "--control-vcf", f"control_b={fixture['control_b.vcf']}",
                                      "--control-vcf", f"other={fixture['other.vcf']}"]))
            ledger = pd.read_csv(out / "focal_event_sharing_ledger.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(set(), set(ledger["event_class"]) & {"private_confirmed"})
            reasons = set(pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("control_callability_missing", reasons)

    def test_missing_control_callset_yields_not_evaluated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._fixture(root)
            out = root / "out"
            self.assertEqual(0, main(["multigenome-private", "--out-dir", str(out), "--case", str(fixture["case.yaml"]),
                                      "--focal-vcf", str(fixture["focal.vcf"]), "--gff3", str(fixture["reference.gff3"]),
                                      "--control-vcf", f"control_a={fixture['control_a.vcf']}",
                                      "--callable-bed", f"control_a={fixture['control_a.bed']}"]))
            ledger = pd.read_csv(out / "focal_event_sharing_ledger.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(ledger["event_class"].eq("not_evaluated").all())
            reasons = set(pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("control_callset_missing", reasons)

    def test_unnormalized_input_is_flagged_for_review(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._fixture(root)
            _vcf(fixture["focal.vcf"], [("chr1", 1200, "A", "T", "1/1")], normalized=False)
            out = root / "out"
            self.assertEqual(0, main(["multigenome-private", "--out-dir", str(out), "--case", str(fixture["case.yaml"]),
                                      "--focal-vcf", str(fixture["focal.vcf"]), "--gff3", str(fixture["reference.gff3"]),
                                      "--control-vcf", f"control_a={fixture['control_a.vcf']}",
                                      "--control-vcf", f"control_b={fixture['control_b.vcf']}",
                                      "--callable-bed", f"control_a={fixture['control_a.bed']}",
                                      "--callable-bed", f"control_b={fixture['control_b.bed']}"]))
            manifest = json.loads((out / "private_variant_manifest.json").read_text(encoding="utf-8"))
            declared = manifest["event_normalization"]["inputs_declaring_normalization"]
            self.assertFalse(declared["focal_parent"])
            self.assertTrue(declared["control_a"])
            self.assertFalse(manifest["event_normalization"]["left_shift_applied"])
            reasons = set(pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("normalization_not_declared", reasons)

    def test_outputs_are_byte_stable_for_the_same_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first, second = root / "first", root / "second"
            self._run(root, first)
            self._run(root, second)
            for name in ("focal_event_sharing_ledger", "candidate_control_matrix", "control_callability_sensitivity",
                         "event_class_summary", "gene_private_variant_ledger",
                         "projection_and_callable_denominators", "working_candidate_panel", "review_required"):
                self.assertEqual((first / f"{name}.tsv").read_bytes(), (second / f"{name}.tsv").read_bytes(), name)

    def test_manifest_records_inputs_and_safety_boundaries(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "out"
            self._run(root, out)
            manifest = json.loads((out / "private_variant_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("focal_parent", manifest["focal_accession"])
            self.assertEqual(["control_a", "control_b"], manifest["phenotype_controls"])
            self.assertEqual(["other"], manifest["comparison_accessions"])
            self.assertEqual({"private_confirmed": 2, "private_possible": 1, "shared_confirmed": 1, "not_evaluated": 0},
                             manifest["event_class_counts"])
            self.assertIn("control_vcf:control_a", manifest["input_fingerprints"])
            self.assertIn("callable_bed:control_b", manifest["input_fingerprints"])
            self.assertTrue(any("never converted into a callable reference allele" in item for item in manifest["safety"]))
            self.assertTrue(any("never summed" in item for item in manifest["safety"]))


class CaseRunIntegrationTests(unittest.TestCase):
    """The funnel is a formal case-run stage between annotation and ranking (#164)."""

    def _case(self, root: Path, with_private: bool) -> tuple[Path, Path]:
        fixture = MultigenomePrivateTests()._fixture(root)
        candidates = root / "candidates.tsv"
        candidates.write_text("gene_id\nGeneA\nGeneB\nGeneC\nGeneD\n", encoding="utf-8")
        evidence = root / "evidence.tsv"
        evidence.write_text(
            "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\n"
            "GeneA\tknown_qtg\tsupport\thigh\tmodel\tliterature\tPMID:1\tgenetic model\n", encoding="utf-8")
        case = root / "run_case.yaml"
        resources = root / "run_resources.yaml"
        private_block = (
            "private_variant_model:\n"
            "  focal_accession: focal_parent\n"
            "  phenotype_controls: [control_a, control_b]\n"
            "  require_all_controls_callable: true\n"
            "  callability_resources:\n"
            "    control_a: callable_a\n"
            "    control_b: callable_b\n"
        ) if with_private else ""
        private_inputs = (
            "  private_focal_vcf_resource: focal_vcf\n"
            "  private_control_vcf_resources:\n"
            "    control_a: control_a_vcf\n"
            "    control_b: control_b_vcf\n"
            "  projection_ledger_resource: projection\n"
            "  reference_frame_consequence_resource: reference_csq\n"
        ) if with_private else ""
        case.write_text(
            "case_id: private_case\n"
            "source:\n"
            "  source_id: PMID:QTL\n"
            "  target_intervals: ['chr1:900-8100']\n"
            "  literature_search_status: complete\n"
            "target:\n"
            "  gff3_resource: gff3\n"
            "inputs:\n"
            "  candidates_resource: candidates\n"
            "  evidence_resources: [evidence]\n"
            f"{private_inputs}"
            "assay:\n"
            "  top_n: 2\n"
            f"{private_block}", encoding="utf-8")
        resources.write_text(
            "resources:\n"
            f"  candidates:\n    path: {candidates}\n"
            f"  evidence:\n    path: {evidence}\n"
            f"  gff3:\n    path: {fixture['reference.gff3']}\n"
            f"  focal_vcf:\n    path: {fixture['focal.vcf']}\n"
            f"  control_a_vcf:\n    path: {fixture['control_a.vcf']}\n"
            f"  control_b_vcf:\n    path: {fixture['control_b.vcf']}\n"
            f"  callable_a:\n    path: {fixture['control_a.bed']}\n"
            f"  callable_b:\n    path: {fixture['control_b.bed']}\n"
            f"  projection:\n    path: {fixture['projection.tsv']}\n"
            f"  reference_csq:\n    path: {fixture['reference_csq.tsv']}\n", encoding="utf-8")
        return case, resources

    def test_case_run_executes_the_private_stage_before_ranking(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case, resources = self._case(root, with_private=True)
            out = root / "run"
            self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources), "--out-dir", str(out)]))
            run_manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("evaluated", run_manifest["private_variant_status"])
            self.assertIn("private_variants", run_manifest["steps"])
            self.assertTrue((out / "private_variants" / "focal_event_sharing_ledger.tsv").exists())
            panel = pd.read_csv(out / "private_variants" / "working_candidate_panel.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(panel["rank_contribution"].eq("none").all())
            # The deterministic rank is unchanged: private status adds no score.
            ranked = pd.read_csv(out / "rank" / "ranked_candidates.tsv", sep="\t", dtype=str).fillna("")
            self.assertNotIn("private_confirmed", ranked.columns)
            self.assertEqual("top", ranked.set_index("gene_id").loc["GeneA", "dry_group"])
            # GeneC/GeneD carry private events but no canonical evidence, so their
            # score stays zero: private status never contributes to the rank.
            scored = ranked.set_index("gene_id")["dry_score"].astype(float).to_dict()
            self.assertEqual(0.0, scored["GeneC"])
            self.assertEqual(0.0, scored["GeneD"])
            self.assertGreater(scored["GeneA"], 0.0)
            report = (out / "report.md").read_text(encoding="utf-8")
            self.assertIn("## Private variants (multigenome)", report)
            self.assertIn("private_possible", report)

    def test_case_without_the_model_records_not_evaluated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case, resources = self._case(root, with_private=False)
            out = root / "run"
            self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources), "--out-dir", str(out)]))
            run_manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("not_requested", run_manifest["private_variant_status"])
            self.assertIsNone(run_manifest["private_variant_manifest"])
            review = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("private_variant_not_evaluated", set(review["reason"]))
            report = (out / "report.md").read_text(encoding="utf-8")
            self.assertIn("NOT biological evidence", report)

    def test_private_stage_failure_does_not_abort_the_pipeline(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case, resources = self._case(root, with_private=True)
            # Break the projection ledger schema after the resources are registered.
            (root / "projection.tsv").write_text("gene_id\nGeneA\n", encoding="utf-8")
            out = root / "run"
            self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources), "--out-dir", str(out)]))
            run_manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("not_evaluated", run_manifest["private_variant_status"])
            review = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("private_variant_stage_failed", set(review["reason"]))
            # The deterministic pipeline still produced its primary artifacts.
            self.assertTrue((out / "rank" / "ranked_candidates.tsv").exists())
            self.assertTrue((out / "report.md").exists())

    def test_resume_reuses_the_private_stage_byte_stably(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case, resources = self._case(root, with_private=True)
            out = root / "run"
            self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources), "--out-dir", str(out)]))
            first = (out / "private_variants" / "focal_event_sharing_ledger.tsv").read_bytes()
            self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources), "--out-dir", str(out), "--resume"]))
            self.assertEqual(first, (out / "private_variants" / "focal_event_sharing_ledger.tsv").read_bytes())


if __name__ == "__main__":
    unittest.main()
