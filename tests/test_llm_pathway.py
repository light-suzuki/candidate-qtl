"""Phenotype -> pathway hypothesis -> verified bait panel route (issue #163)."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.claim_binding import BOUND_CLAIM_COLUMNS
from tests.bound_claim_fixtures import bound_claim


def _uniprot_record(accession: str, species: str, sequence: str, *, reviewed: bool = True, pmid: str = "12345") -> dict:
    record: dict = {
        "primaryAccession": accession,
        "uniProtkbId": f"{accession}_TEST",
        "entryType": "UniProtKB reviewed (Swiss-Prot)" if reviewed else "UniProtKB unreviewed (TrEMBL)",
        "organism": {"scientificName": species},
        "sequence": {"value": sequence, "length": len(sequence)},
    }
    if pmid:
        record["references"] = [{"citation": {"citationCrossReferences": [{"database": "PubMed", "id": pmid}]}}]
    return record


RAW_OUTPUT = {
    "hypotheses": [
        {"pathway_or_process": "CLAVATA-WUSCHEL meristem maintenance", "gene_family": "CLE receptor kinase",
         "reference_gene_name": "CLV1", "reference_species": "Arabidopsis thaliana",
         "proposed_stable_id": "P0CLV1", "expected_role": "receptor kinase restricting stem cell pool",
         "mechanistic_link": "loss of function enlarges the meristem",
         "counterevidence": "meristem size can also change through hormone transport",
         "confidence": "medium", "uncertainty": "family member identity in the focal species"},
        {"pathway_or_process": "auxin transport", "gene_family": "PIN efflux carrier",
         "reference_gene_name": "PIN1", "reference_species": "Arabidopsis thaliana",
         "proposed_stable_id": "P0PIN1", "expected_role": "polar auxin efflux",
         "mechanistic_link": "altered phyllotaxis and fasciation",
         "counterevidence": "auxin defects usually affect more organs",
         "confidence": "low", "uncertainty": "direction of the effect"},
        {"pathway_or_process": "hallucinated pathway", "gene_family": "invented family",
         "reference_gene_name": "NOPE1", "reference_species": "Arabidopsis thaliana",
         "proposed_stable_id": "P0FAKE0", "expected_role": "unknown",
         "mechanistic_link": "none stated", "counterevidence": "", "confidence": "high", "uncertainty": ""},
        {"pathway_or_process": "species mismatch check", "gene_family": "CLE receptor kinase",
         "reference_gene_name": "HAR1", "reference_species": "Medicago truncatula",
         "proposed_stable_id": "P0LOTUS", "expected_role": "receptor kinase",
         "mechanistic_link": "nodule autoregulation", "counterevidence": "", "confidence": "medium", "uncertainty": ""},
        {"pathway_or_process": "focal species bait", "gene_family": "CLE receptor kinase",
         "reference_gene_name": "PsCLV1", "reference_species": "Pisum sativum",
         "proposed_stable_id": "P0PEA00", "expected_role": "receptor kinase",
         "mechanistic_link": "same pathway in the focal species", "counterevidence": "",
         "confidence": "high", "uncertainty": ""},
    ]
}


class PromptExportTests(unittest.TestCase):
    def test_blind_prompt_is_versioned_and_audited_against_interval_candidates(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = root / "candidates.tsv"
            candidates.write_text("gene_id\nPsat5g123456\nPsat5g123457\n", encoding="utf-8")
            out = root / "prompt"
            self.assertEqual(0, main(["llm-pathway-prompt", "--out-dir", str(out), "--case-id", "fasciation",
                                      "--focal-species", "Pisum sativum", "--phenotype", "stem fasciation",
                                      "--direction", "gain", "--inheritance", "recessive", "--organ", "shoot apex",
                                      "--stage", "vegetative", "--candidates", str(candidates)]))
            manifest = json.loads((out / "llm_input_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("blind", manifest["execution_mode"])
            self.assertFalse(manifest["qtl_candidate_information_in_model_input"])
            self.assertTrue(manifest["blind_leak_audit"]["checked"])
            self.assertEqual([], manifest["blind_leak_audit"]["leaked_identifiers"])
            self.assertEqual(2, manifest["blind_leak_audit"]["candidate_count"])
            self.assertTrue(manifest["prompt_sha256"])
            self.assertTrue(manifest["prompt_version"])
            prompt = (out / "prompt.txt").read_text(encoding="utf-8")
            self.assertNotIn("Psat5g123456", prompt)
            self.assertIn("stem fasciation", prompt)
            self.assertIn("same-species primary reference", prompt)
            self.assertIn("strict benchmark", prompt)
            record = json.loads((out / "llm_pathway_prompts.jsonl").read_text(encoding="utf-8").strip())
            self.assertEqual({"prompt_id", "gene_id", "system", "user"}, set(record))

    def test_blind_mode_rejects_researcher_context(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            context = root / "context.txt"
            context.write_text("researcher hypothesis: CLV pathway\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                main(["llm-pathway-prompt", "--out-dir", str(root / "p"), "--case-id", "c",
                      "--focal-species", "Pisum sativum", "--phenotype", "fasciation",
                      "--assisted-context", str(context)])

    def test_assisted_mode_records_the_context_separately(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            context = root / "context.txt"
            context.write_text("verified literature: PMID:1 CLV1 restricts meristem size\n", encoding="utf-8")
            out = root / "prompt"
            self.assertEqual(0, main(["llm-pathway-prompt", "--out-dir", str(out), "--case-id", "c",
                                      "--focal-species", "Pisum sativum", "--phenotype", "fasciation",
                                      "--mode", "assisted", "--assisted-context", str(context)]))
            manifest = json.loads((out / "llm_input_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("assisted", manifest["execution_mode"])
            self.assertTrue(manifest["assisted_context_supplied"])
            self.assertFalse(manifest["blind_leak_audit"]["checked"])

    def test_candidate_identifier_in_a_blind_prompt_is_a_hard_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = root / "candidates.tsv"
            candidates.write_text("gene_id\nPsat5g123456\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                main(["llm-pathway-prompt", "--out-dir", str(root / "p"), "--case-id", "c",
                      "--focal-species", "Pisum sativum",
                      "--phenotype", "fasciation near Psat5g123456", "--candidates", str(candidates)])


class ImportTests(unittest.TestCase):
    def _import(self, root: Path, payload: object) -> Path:
        raw = root / "raw_model_output.json"
        raw.write_text(json.dumps(payload), encoding="utf-8")
        out = root / "import"
        self.assertEqual(0, main(["llm-pathway-import", "--raw", str(raw), "--out-dir", str(out),
                                  "--case-id", "fasciation", "--model", "test-model-1", "--proposal-source", "unit-test"]))
        return out

    def test_raw_and_parsed_output_are_both_kept(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = self._import(root, RAW_OUTPUT)
            hypotheses = pd.read_csv(out / "llm_pathway_hypotheses.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(5, len(hypotheses))
            self.assertTrue(hypotheses["review_status"].eq("pending_verification").all())
            self.assertTrue(hypotheses["model"].eq("test-model-1").all())
            self.assertTrue((out / "raw_model_output.json").exists())
            manifest = json.loads((out / "llm_pathway_import_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(5, manifest["proposal_count"])
            self.assertTrue(manifest["raw_output_sha256"])
            self.assertGreaterEqual(manifest["distinct_pathways"], 2)

    def test_fenced_and_noisy_output_is_still_parsed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root / "raw.txt"
            raw.write_text("Here you go:\n```json\n" + json.dumps(RAW_OUTPUT) + "\n```\n", encoding="utf-8")
            out = root / "import"
            self.assertEqual(0, main(["llm-pathway-import", "--raw", str(raw), "--out-dir", str(out), "--case-id", "c"]))
            self.assertEqual(5, len(pd.read_csv(out / "llm_pathway_hypotheses.tsv", sep="\t", dtype=str)))

    def test_single_pathway_output_is_flagged_as_overconfident(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = self._import(root, {"hypotheses": [RAW_OUTPUT["hypotheses"][0]]})
            reasons = set(pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("single_pathway_overconfidence", reasons)

    def test_unparsable_output_abstains_instead_of_inventing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root / "raw.txt"
            raw.write_text("I am not sure which pathway is responsible.\n", encoding="utf-8")
            out = root / "import"
            self.assertEqual(0, main(["llm-pathway-import", "--raw", str(raw), "--out-dir", str(out), "--case-id", "c"]))
            self.assertEqual(0, len(pd.read_csv(out / "llm_pathway_hypotheses.tsv", sep="\t", dtype=str)))
            reasons = set(pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("no_parsable_proposal", reasons)


class VerificationAndFreezeTests(unittest.TestCase):
    def _pipeline(self, root: Path, *, benchmark_mode: bool = True) -> tuple[Path, Path, Path]:
        raw = root / "raw.json"
        raw.write_text(json.dumps(RAW_OUTPUT), encoding="utf-8")
        imported = root / "import"
        self.assertEqual(0, main(["llm-pathway-import", "--raw", str(raw), "--out-dir", str(imported),
                                  "--case-id", "fasciation", "--model", "test-model-1"]))
        cache = root / "records.json"
        cache.write_text(json.dumps({"results": [
            _uniprot_record("P0CLV1", "Arabidopsis thaliana", "MKLLTAVWQ" * 12),
            _uniprot_record("P0PIN1", "Arabidopsis thaliana", "MGTAAALPD" * 12),
            # Declared Medicago, actually Lotus: a species mismatch.
            _uniprot_record("P0LOTUS", "Lotus japonicus", "MSSLRRLLL" * 12),
            _uniprot_record("P0PEA00", "Pisum sativum", "MPPQQRLAA" * 12),
        ]}), encoding="utf-8")
        hypotheses = pd.read_csv(imported / "llm_pathway_hypotheses.tsv", sep="\t", dtype=str).fillna("")
        bindings = root / "bound_claims.tsv"
        pd.DataFrame([
            bound_claim(
                route="a",
                claim_id=str(row["hypothesis_id"]),
                citation="PMID:12345",
                subject_gene=str(row["proposed_stable_id"]),
                organism=str(row["reference_species"]),
                phenotype=str(row["pathway_or_process"]),
            )
            for _, row in hypotheses.iterrows()
            if str(row["proposed_stable_id"])
        ], columns=BOUND_CLAIM_COLUMNS).to_csv(bindings, sep="\t", index=False, lineterminator="\n")
        verified = root / "verify"
        verify_args = ["hypothesis-verify", "--hypotheses", str(imported / "llm_pathway_hypotheses.tsv"),
                       "--out-dir", str(verified), "--focal-species", "Pisum sativum",
                       "--record-cache", str(cache), "--bound-claims", str(bindings)]
        if benchmark_mode:
            verify_args.append("--benchmark-mode")
        self.assertEqual(0, main(verify_args))
        frozen = root / "freeze"
        self.assertEqual(0, main(["bait-panel-freeze", "--hypotheses", str(imported / "llm_pathway_hypotheses.tsv"),
                                  "--validation", str(verified / "proposal_validation.tsv"),
                                  "--sequences", str(verified / "verified_sequences.json"),
                                  "--out-dir", str(frozen), "--case-id", "fasciation"]))
        return imported, verified, frozen

    def test_hallucinated_accession_never_reaches_the_panel(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, verified, frozen = self._pipeline(root)
            validation = pd.read_csv(verified / "proposal_validation.tsv", sep="\t", dtype=str).fillna("").set_index("proposed_stable_id")
            self.assertEqual("review_required", validation.loc["P0FAKE0", "verification_status"])
            self.assertEqual("verification_not_performed_offline", validation.loc["P0FAKE0", "review_reason"])
            panel = pd.read_csv(frozen / "verified_bait_panel.tsv", sep="\t", dtype=str).fillna("")
            self.assertNotIn("P0FAKE0", set(panel["bait_id"]))

    def test_species_mismatch_and_focal_species_baits_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, verified, frozen = self._pipeline(root)
            validation = pd.read_csv(verified / "proposal_validation.tsv", sep="\t", dtype=str).fillna("").set_index("proposed_stable_id")
            self.assertEqual("rejected", validation.loc["P0LOTUS", "verification_status"])
            self.assertEqual("declared_species_does_not_match_accession", validation.loc["P0LOTUS", "rejection_reason"])
            self.assertEqual("rejected", validation.loc["P0PEA00", "verification_status"])
            self.assertIn(validation.loc["P0PEA00", "rejection_reason"],
                          {"reference_gene_is_from_the_focal_species",
                           "answer_derived_or_focal_lineage_marker_in_benchmark_mode"})
            panel = set(pd.read_csv(frozen / "verified_bait_panel.tsv", sep="\t", dtype=str).fillna("")["bait_id"])
            self.assertEqual({"P0CLV1", "P0PIN1"}, panel)
            # Rejected proposals are retained with their reason, not deleted.
            rejected = pd.read_csv(verified / "proposal_rejected.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"P0LOTUS", "P0PEA00"}, set(rejected["proposed_stable_id"]))
            self.assertTrue(rejected["rejection_reason"].ne("").all())

    def test_production_mode_allows_verified_same_species_reference(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, verified, frozen = self._pipeline(root, benchmark_mode=False)
            validation = pd.read_csv(verified / "proposal_validation.tsv", sep="\t", dtype=str).fillna("").set_index("proposed_stable_id")
            self.assertEqual("verified", validation.loc["P0PEA00", "verification_status"])
            self.assertEqual("False", validation.loc["P0PEA00", "focal_species_excluded"])
            # The fixture has no provider-side gene alias, so the independent
            # identity projection remains gated even though proposal verification
            # and sequence retrieval are accepted in production mode.
            self.assertEqual("not_verified", validation.loc["P0PEA00", "identity_sequence_verification_status"])
            panel = set(pd.read_csv(frozen / "verified_bait_panel.tsv", sep="\t", dtype=str).fillna("")["bait_id"])
            self.assertIn("P0PEA00", panel)

    def test_frozen_panel_carries_provenance_and_sequence_fingerprints(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, _, frozen = self._pipeline(root)
            panel = pd.read_csv(frozen / "verified_bait_panel.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(panel["sequence_sha256"].ne("").all())
            self.assertTrue(panel["source_provider"].ne("").all())
            self.assertTrue(panel["retrieved_at"].ne("").all())
            self.assertTrue(panel["function_source_ids"].str.contains("PMID:").all())
            self.assertTrue(panel["rank_contribution"].eq("none").all())
            # Criterion 2: the frozen panel carries a per-bait review_status
            # derived from the panel status (verified for a cleanly frozen panel).
            self.assertIn("review_status", panel.columns)
            self.assertTrue(panel["review_status"].eq("verified").all())
            manifest = json.loads((frozen / "bait_panel_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(2, manifest["bait_count"])
            self.assertEqual("frozen", manifest["panel_status"])
            self.assertEqual("verified", manifest["bait_review_status"])
            self.assertTrue(manifest["fasta_sha256"])
            self.assertEqual(set(panel["bait_id"]), set(manifest["sequence_fingerprints"]))
            fasta = (frozen / "verified_bait_panel.fasta").read_text(encoding="utf-8")
            self.assertIn(">P0CLV1 Arabidopsis thaliana", fasta)

    @unittest.skip("public release: query-eligibility is a benchmark command (not shipped)")
    def test_frozen_panel_passes_the_existing_query_eligibility_gate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, _, frozen = self._pipeline(root)
            gate = root / "gate"
            self.assertEqual(0, main(["query-eligibility", "--query-manifest", str(frozen / "verified_bait_panel_query_manifest.tsv"),
                                      "--out-dir", str(gate), "--strict"]))
            manifest = json.loads((gate / "query_eligibility_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(2, manifest["eligible_count"])
            self.assertEqual(0, manifest["excluded_count"])
            self.assertFalse(manifest["pisum_queries_executed"])

    def test_unverifiable_proposals_stay_review_required_offline(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root / "raw.json"
            raw.write_text(json.dumps(RAW_OUTPUT), encoding="utf-8")
            imported = root / "import"
            self.assertEqual(0, main(["llm-pathway-import", "--raw", str(raw), "--out-dir", str(imported), "--case-id", "c"]))
            verified = root / "verify"
            self.assertEqual(0, main(["hypothesis-verify", "--hypotheses", str(imported / "llm_pathway_hypotheses.tsv"),
                                      "--out-dir", str(verified), "--focal-species", "Pisum sativum"]))
            validation = pd.read_csv(verified / "proposal_validation.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(validation["verification_status"].eq("review_required").all())
            manifest = json.loads((verified / "hypothesis_verification_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(0, manifest["verified_count"])
            self.assertEqual("none", manifest["verification_provider"])
            self.assertFalse(manifest["external_api_allowed"])

    def test_single_bait_panel_is_review_required(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root / "raw.json"
            raw.write_text(json.dumps({"hypotheses": RAW_OUTPUT["hypotheses"][:1]}), encoding="utf-8")
            imported = root / "import"
            self.assertEqual(0, main(["llm-pathway-import", "--raw", str(raw), "--out-dir", str(imported), "--case-id", "c"]))
            cache = root / "records.json"
            cache.write_text(json.dumps({"results": [_uniprot_record("P0CLV1", "Arabidopsis thaliana", "MKLLTAVWQ" * 12)]}), encoding="utf-8")
            hypothesis = pd.read_csv(imported / "llm_pathway_hypotheses.tsv", sep="\t", dtype=str).fillna("").iloc[0]
            bindings = root / "bound_claims.tsv"
            pd.DataFrame([bound_claim(
                route="a", claim_id=str(hypothesis["hypothesis_id"]), citation="PMID:12345",
                subject_gene="P0CLV1", organism="Arabidopsis thaliana",
                phenotype=str(hypothesis["pathway_or_process"]),
            )], columns=BOUND_CLAIM_COLUMNS).to_csv(bindings, sep="\t", index=False, lineterminator="\n")
            verified = root / "verify"
            self.assertEqual(0, main(["hypothesis-verify", "--hypotheses", str(imported / "llm_pathway_hypotheses.tsv"),
                                      "--out-dir", str(verified), "--focal-species", "Pisum sativum",
                                      "--record-cache", str(cache), "--bound-claims", str(bindings)]))
            frozen = root / "freeze"
            self.assertEqual(0, main(["bait-panel-freeze", "--hypotheses", str(imported / "llm_pathway_hypotheses.tsv"),
                                      "--validation", str(verified / "proposal_validation.tsv"),
                                      "--sequences", str(verified / "verified_sequences.json"), "--out-dir", str(frozen)]))
            manifest = json.loads((frozen / "bait_panel_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("review_required_small_panel", manifest["panel_status"])
            reasons = set(pd.read_csv(frozen / "review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("panel_below_minimum_size", reasons)
            # Panel diversity is still review-required, but the bait itself
            # passed accession/source/species/sequence verification and remains
            # usable for non-scoring deterministic interval projection.
            panel = pd.read_csv(frozen / "verified_bait_panel.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("review_status", panel.columns)
            self.assertTrue(panel["review_status"].eq("verified").all())

    def test_route_never_emits_scores_or_canonical_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            imported, verified, frozen = self._pipeline(root)
            for directory in (imported, verified, frozen):
                for path in sorted(directory.glob("*.tsv")):
                    frame = pd.read_csv(path, sep="\t", dtype=str)
                    for forbidden in ("dry_score", "score", "evidence_code", "top_eligibility", "dry_group"):
                        self.assertNotIn(forbidden, frame.columns, f"{path.name}:{forbidden}")
                self.assertEqual([], sorted(p.name for p in directory.glob("evidence*.tsv")))
            # Model confidence is retained as a proposal attribute only; it never
            # appears in the frozen panel that downstream mapping consumes.
            panel = pd.read_csv(frozen / "verified_bait_panel.tsv", sep="\t", dtype=str)
            self.assertNotIn("confidence", panel.columns)


class BenchmarkTests(unittest.TestCase):
    def test_staged_metrics_are_reported_in_a_separate_lane(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pipeline = VerificationAndFreezeTests()
            imported, verified, frozen = pipeline._pipeline(root)
            ranked = root / "ranked.tsv"
            ranked.write_text("gene_id\nPsat5g000001\nPsat5g000002\nPsatCausal\n", encoding="utf-8")
            config = root / "config.tsv"
            config.write_text(
                "case_id\texpected_pathway\texpected_gene_family\tcausal_gene_id\thypotheses\tvalidation\tpanel\tranked\n"
                f"fasciation\tCLAVATA-WUSCHEL meristem maintenance\tCLE receptor kinase\tPsatCausal\t"
                f"{imported / 'llm_pathway_hypotheses.tsv'}\t{verified / 'proposal_validation.tsv'}\t"
                f"{frozen / 'verified_bait_panel.tsv'}\t{ranked}\n", encoding="utf-8")
            out = root / "benchmark"
            self.assertEqual(0, main(["phenotype-pathway-benchmark", "--config", str(config), "--out-dir", str(out)]))
            cases = pd.read_csv(out / "phenotype_pathway_cases.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("1", cases.iloc[0]["pathway_rank"])
            self.assertEqual("True", cases.iloc[0]["pathway_top1"])
            self.assertEqual("True", cases.iloc[0]["gene_family_recalled"])
            self.assertEqual("True", cases.iloc[0]["causal_in_top20"])
            self.assertEqual("", cases.iloc[0]["failure_reason"])
            manifest = json.loads((out / "phenotype_pathway_benchmark_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("phenotype_only_proposal_lane", manifest["lane"])
            self.assertIn("never be pooled", manifest["separation_rule"])
            metrics = manifest["metrics"]
            self.assertEqual(1.0, metrics["pathway_top1_recall"])
            self.assertEqual(1.0, metrics["causal_gene_family_recall"])
            self.assertEqual(0.4, metrics["verification_pass_rate"])
            self.assertEqual(0.0, metrics["abstention_rate"])

    def test_failure_reasons_are_broken_down_per_stage(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pipeline = VerificationAndFreezeTests()
            imported, verified, frozen = pipeline._pipeline(root)
            config = root / "config.tsv"
            config.write_text(
                "case_id\texpected_pathway\texpected_gene_family\tcausal_gene_id\thypotheses\tvalidation\tpanel\tranked\n"
                f"miss\tsomething else entirely\tunrelated family\t\t"
                f"{imported / 'llm_pathway_hypotheses.tsv'}\t{verified / 'proposal_validation.tsv'}\t"
                f"{frozen / 'verified_bait_panel.tsv'}\t\n", encoding="utf-8")
            out = root / "benchmark"
            self.assertEqual(0, main(["phenotype-pathway-benchmark", "--config", str(config), "--out-dir", str(out)]))
            cases = pd.read_csv(out / "phenotype_pathway_cases.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("pathway_miss", cases.iloc[0]["failure_reason"])
            manifest = json.loads((out / "phenotype_pathway_benchmark_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual({"pathway_miss": 1}, manifest["metrics"]["failure_reasons"])
            self.assertEqual(0.0, manifest["metrics"]["pathway_top5_recall"])


if __name__ == "__main__":
    unittest.main()
