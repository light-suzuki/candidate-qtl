"""T1 record warnings reach A/B Wet outputs without becoming selection rules."""
import hashlib
import json
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.observation_tiers import checklist_consistency
from candidate_qtl.paper_reports import PAPER_FIELD
from candidate_qtl.relevance_tiers import NEW_TIERS, annotate_relevance, classify_relevance
from candidate_qtl.route_b_semantic_bundle import build_route_b_semantic_bundle
from candidate_qtl.route_b_semantic_source_verify import verify_route_b_semantic_sources
from candidate_qtl.route_b_semantic_transfer_verify import verify_route_b_semantic_transfer
from candidate_qtl.route_b_semantic_candidates import freeze_route_b_semantic_candidates
from candidate_qtl.source_verify import source_verify
from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from tests import test_paper_tiers as papers
from tests import test_route_a_review_relevance as a
from tests import test_route_b_semantic_import as b
from tests import test_paper_evidence_table as table_tests


FOCAL = "Pisum sativum"
REPORTED = "Arabidopsis thaliana"


def t1_paper(**changes):
    return {**table_tests.paper(phenotype_reported="yes", equivalent_target_phenotype="yes",
                   focal_species="yes", same_gene_as_focal="yes"),
            "reported_species": REPORTED, **changes}


def read_table(path):
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


class PaperSpeciesConsistencyTests(unittest.TestCase):
    def test_t1_species_disagreement_warns_without_reclassifying_or_promoting_support(self):
        report = t1_paper()
        original = deepcopy(report)
        roster = [{"candidate_gene_id": "G1", "support_status": "unverified_route_b_review"}]
        claim = papers.projected_claim([report], focal_species=FOCAL)
        row = annotate_relevance(roster, [claim])[0]
        self.assertEqual(NEW_TIERS[0], row["relevance_tier"])
        self.assertEqual("review", row["support_state"])
        summary = json.loads(row["paper_evidence_summary_json"])[0]
        self.assertEqual("explicit_label_disagreement", summary["species_comparison_state"])
        self.assertIn("taxonomic identity, synonymy and rank were not resolved", summary["consistency_warning"])
        self.assertEqual(FOCAL, summary["focal_species_label"])
        self.assertEqual(REPORTED, summary["reported_species"])
        self.assertEqual(original, report)
        replay = annotate_relevance([row])[0]
        self.assertEqual(row["paper_evidence_summary_json"], replay["paper_evidence_summary_json"])
        self.assertEqual(row["relevance_tier"], replay["relevance_tier"])

    def test_same_species_cross_species_t2_and_focal_report_via_proxy_are_preserved(self):
        similar = t1_paper()
        similar["paper_checklist"].update(focal_species="no", same_gene_as_focal="no")
        cases = [(t1_paper(reported_species=FOCAL), "T1"), (similar, "T2"),
                 (t1_paper(reported_species=FOCAL, gene_correspondence="homolog"), "T1")]
        for report, tier in cases:
            with self.subTest(tier=tier, relation=report["gene_correspondence"]):
                claim = papers.projected_claim([report], focal_species=FOCAL, reference_species=REPORTED)
                result = annotate_relevance([{"candidate_gene_id": "G1"}], [claim])[0]
                self.assertEqual(tier, result["initial_tier"])
                self.assertFalse(json.loads(result["paper_evidence_summary_json"])[0]["consistency_warning"])

    def test_unresolved_and_formatted_species_labels_never_create_negative_evidence(self):
        for label in ("pea", "garden pea", "P. sativum", "unknown", "不明", "", "source common name"):
            with self.subTest(label=label):
                report = t1_paper(reported_species=label)
                result = checklist_consistency({**report, "focal_species": FOCAL})
                self.assertEqual("unresolved_labels", result["species_comparison_state"])
                self.assertFalse(result["consistency_warning"])
                self.assertEqual(NEW_TIERS[0], classify_relevance(papers.projected_claim([report], focal_species=FOCAL))[0])
        for label in ("  Pisum   sativum  ", "PISUM SATIVUM", "Ｐisum sativum", "Pisum sativum (pea)"):
            with self.subTest(label=label):
                result = checklist_consistency({**t1_paper(reported_species=label), "focal_species": FOCAL})
                self.assertEqual("matching_labels", result["species_comparison_state"])
                self.assertFalse(result["consistency_warning"])
        # Different scientific synonyms are only a label correction lead: this
        # audit cannot establish taxonomic identity or lower the proposed Tier.
        report = t1_paper(reported_species="Lycopersicon esculentum")
        claim = papers.projected_claim([report], focal_species="Solanum lycopersicum")
        self.assertEqual(NEW_TIERS[0], classify_relevance(claim)[0])
        self.assertIn("not resolved", checklist_consistency({**report, "focal_species": claim["focal_species"]})["consistency_warning"])

    def test_explicit_taxon_ids_compare_without_names_or_proxy_species(self):
        for focal, reported, state in (("NCBITaxon:3888", "taxid:3888", "matching_labels"),
                                        ("NCBITaxon:3888", "NCBI Taxon:3702", "explicit_label_disagreement"),
                                        ("NCBITaxon:3888", FOCAL, "unresolved_labels")):
            with self.subTest(focal=focal, reported=reported):
                row = {**t1_paper(reported_species=reported), "focal_species": focal, "reference_species": "unrelated proxy"}
                self.assertEqual(state, checklist_consistency(row)["species_comparison_state"])

    def test_warning_does_not_authorize_same_name_with_a_different_reference_id(self):
        # A species warning never repairs the existing independent anchor gate.
        report = t1_paper(reference_gene_id="GeneID:456", reported_gene="GENE_A")
        claim = papers.projected_claim([report], focal_species=FOCAL, reference_gene="GENE_A (GeneID:123)")
        self.assertEqual("U_unresolved", classify_relevance(claim)[0])
        report["reference_gene_id"] = "GeneID:123"
        self.assertEqual(NEW_TIERS[0], classify_relevance(papers.projected_claim([report], focal_species=FOCAL,
            reference_gene="GENE_A (GeneID:123)"))[0])

    def test_same_paper_display_representative_cannot_hide_another_records_warning(self):
        reports = [t1_paper(reported_species=FOCAL, reported_gene="A"), t1_paper(reported_gene="Z")]
        row = annotate_relevance([{"candidate_gene_id": "G1"}], [papers.projected_claim(reports, focal_species=FOCAL)])[0]
        summaries = json.loads(row["paper_evidence_summary_json"])
        self.assertEqual({"T1": 1}, json.loads(row["paper_tier_counts_json"]))
        self.assertEqual(1, len(summaries))
        self.assertEqual(["Z"], [item["reported_gene"] for item in summaries[0]["consistency_warnings"]])

    def test_consistency_lead_uses_only_one_existing_followup_batch(self):
        report = t1_paper()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            exporter = table_tests.PaperEvidenceTableTests()
            first, table, packet = exporter.export(root, [report], focal_species=FOCAL)
            self.assertEqual("needs_completion", first["recording_state"])
            self.assertEqual("T1", table.iloc[0].applicable_tier)
            self.assertEqual("clarify_record_consistency", table.iloc[0].followup_action)
            task = packet["requests"][0]["records"][0]
            self.assertEqual(1, packet["max_followup_batches"])
            self.assertEqual(FOCAL, task["focal_species_label"])
            self.assertIn("focal_species", task["questions"])
            report.update(followup_state="completed", unknown_reason="The supplied taxonomy labels remain unresolved after checking.")
            _, table, packet = exporter.export(root, [report], focal_species=FOCAL)
            self.assertEqual([], packet["requests"])
            self.assertEqual("T1", table.iloc[0].applicable_tier)
            self.assertTrue(table.iloc[0].consistency_warning)
            report["paper_checklist"]["focal_species"] = "no"
            _, table, packet = exporter.export(root, [report], focal_species=FOCAL)
            self.assertEqual("T2", table.iloc[0].applicable_tier)
            self.assertFalse(table.iloc[0].consistency_warning)
            self.assertEqual([], packet["requests"])

    def test_legacy_declared_and_observation_policy_do_not_gain_a_checklist_warning(self):
        for report, expected in ((papers.paper(reported_species=REPORTED), NEW_TIERS[3]),
                                 (papers.paper(reported_species=REPORTED, evidence_types=["phenotype"]), "U_unresolved")):
            claim = papers.claim([report], focal_species=FOCAL)
            before = deepcopy(claim)
            self.assertEqual({}, checklist_consistency({**report, "focal_species": FOCAL}))
            row = annotate_relevance([{"candidate_gene_id": "G1"}], [claim])[0]
            summary = json.loads(row["paper_evidence_summary_json"])[0]
            self.assertNotIn("consistency_warning", summary)
            self.assertNotIn("consistency_warnings", summary)
            self.assertEqual(before, claim)
            self.assertEqual(expected, row["relevance_tier"])

    def test_a_case_species_is_audit_only_for_legacy_and_mixed_paper_reports(self):
        legacy = papers.paper(source="PMID:456", reported_species=FOCAL, evidence_types=["phenotype"])
        for reports in ([legacy], [legacy, t1_paper()]):
            with self.subTest(report_count=len(reports)), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                investigation = root / "investigation.tsv"
                pd.DataFrame([{"case_id": "audit-only", "focal_species": FOCAL}]).to_csv(investigation, sep="\t", index=False)
                payload = {**papers.MECHANISMS[0], "reference_gene_name": "REF1", "reference_species": FOCAL,
                           "source_ids": ";".join(report["source_id"] for report in reports), PAPER_FIELD: reports}
                raw = papers._payload(root, [payload])
                original = raw.read_bytes()
                tables = []
                for name, args in (("without-case", []), ("with-case", ["--investigation", str(investigation)])):
                    self.assertEqual(0, main(["mechanism-import", "--payload", str(raw),
                        "--out-dir", str(root / name), *args]))
                    tables.append(read_table(root / name / "paper_evidence_table.tsv").set_index("source_id"))
                # Supplying the case may add a checklist warning, but cannot
                # turn the existing legacy species-unknown T2 into a focal T1.
                pd.testing.assert_series_equal(tables[0].loc["PMID:456"], tables[1].loc["PMID:456"])
                self.assertEqual("T2", tables[1].loc["PMID:456", "computed_tier"])
                self.assertEqual("T2", tables[1].loc["PMID:456", "applicable_tier"])
                if len(reports) > 1:
                    self.assertFalse(tables[0].loc["PMID:123", "consistency_warning"])
                    self.assertTrue(tables[1].loc["PMID:123", "consistency_warning"])
                    self.assertEqual("T1", tables[1].loc["PMID:123", "applicable_tier"])
                self.assertEqual(original, raw.read_bytes())

    def _assert_wet_warning(self, wet, expected_ids, expected_warning=True):
        adopted = read_table(wet["wet_candidates"])
        trace = read_table(Path(wet["wet_candidates"]).parent / "wet_candidate_trace.tsv")
        self.assertEqual(expected_ids, set(adopted.candidate_gene_id))
        self.assertEqual({NEW_TIERS[0]}, set(adopted.relevance_tier))
        self.assertEqual({"review"}, set(adopted.support_state))
        for frame in (adopted, trace):
            for summary in frame.paper_evidence_summary_json:
                report = json.loads(summary)[0]
                self.assertEqual(expected_warning, bool(report["consistency_warning"]))
        self.assertEqual({"False"}, set(trace.supported_union_member))
        self.assertEqual({"True"}, set(trace.adopted_first_stage))
        self.assertEqual({"not_established"}, set(trace.biological_exclusion))

    def test_a_cli_import_table_source_projection_wet_trace_keep_the_same_record(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            report = t1_paper(reference_gene_id="bait1")
            investigation = root / "investigation.tsv"
            pd.DataFrame([{"case_id": "route-a-review", "focal_species": FOCAL}]).to_csv(investigation, sep="\t", index=False)
            payload = {**papers.MECHANISMS[0], "hypothesis_id": "h1", "proposed_stable_id": "bait1",
                       "reference_gene_name": "REF1", "reference_species": REPORTED,
                       "source_ids": "PMID:123", PAPER_FIELD: [report]}
            raw = papers._payload(root, [payload])
            original = raw.read_bytes()
            self.assertEqual(0, main(["mechanism-import", "--payload", str(raw), "--investigation", str(investigation),
                                      "--out-dir", str(root / "import")]))
            table = read_table(root / "import/paper_evidence_table.tsv")
            self.assertEqual(FOCAL, table.iloc[0].focal_species_label)
            self.assertTrue(table.iloc[0].consistency_warning)
            self.assertEqual("T1", table.iloc[0].applicable_tier)
            frozen_table = (root / "import/paper_evidence_table.tsv").read_bytes()
            cache = root / "cache.json"
            cache.write_text(json.dumps([{"id": "123", "title": report["source_title"]}]))
            verified = source_verify(str(root / "import/mechanism_hypotheses.tsv"), str(root / "sources"),
                                     cache_arg=str(cache), focal_species=FOCAL, verification_policy="initial_bibliography_v1")
            checked = read_table(verified["mechanisms_with_verified_sources"]).iloc[0].to_dict()
            self.assertEqual("not_evaluated", checked["claim_support_summary"])
            consumer_root = root / "consumer"
            consumer_root.mkdir()
            wet, _, _, _ = a.RouteAReviewRelevanceTests()._run(consumer_root, checked, alignment_span=90)
            self._assert_wet_warning(wet, {"g1", "g2"})
            self.assertEqual(original, raw.read_bytes())
            self.assertEqual(frozen_table, (root / "import/paper_evidence_table.tsv").read_bytes())

    def test_b_cli_import_source_transfer_freeze_wet_trace_preserve_members_and_hashes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = b.RouteBSemanticImportTests()
            fixture._bundle(root)
            built = build_route_b_semantic_bundle(str(root / "atlas.tsv"), str(root / "focal-bundle"),
                phenotype="fasciated meristem", literature_arg=str(root / "literature.tsv"),
                case_id="case403", focal_species=FOCAL)
            bundle_before = Path(built["bundle"]).read_bytes()
            frozen_inputs = (built["bundle"], built["prompts"])
            runs = []
            for name, species in (("conflict", REPORTED), ("same", FOCAL)):
                report = t1_paper(reported_species=species)
                response = root / f"{name}.json"
                model_claim = {**papers.observation(), "claim_type": "mutant_phenotype", "reference_gene": "REF1",
                               "reference_species": REPORTED, "evidence_scope": "reference_gene",
                               "source_ids": ["PMID:123"], PAPER_FIELD: [report]}
                fixture._response(Path(built["prompts"]), [model_claim], response)
                hashes = {path: hashlib.sha256(Path(path).read_bytes()).hexdigest() for path in (*frozen_inputs, str(response))}
                self.assertEqual(0, main(["route-b-semantic-import", "--bundle", built["bundle"], "--prompts", built["prompts"],
                    "--responses", str(response), "--out-dir", str(root / name / "import")]))
                imported = json.loads((root / name / "import/semantic_claim_proposal_manifest.json").read_text())
                proposals = root / name / "import/semantic_claim_proposals.tsv"
                table = read_table(root / name / "import/paper_evidence_table.tsv")
                self.assertEqual(name == "conflict", bool(table.iloc[0].consistency_warning))
                self.assertEqual("T1", table.iloc[0].computed_tier)
                self.assertEqual(0, imported["verified_claim_row_count"])
                cache = root / name / "cache.json"
                cache.write_text(json.dumps([{"id": "123", "title": report["source_title"]}]))
                sources = verify_route_b_semantic_sources(str(proposals), str(root / name / "sources"),
                    bundle_arg=built["bundle"], cache_arg=str(cache))
                transfer = verify_route_b_semantic_transfer(str(proposals), sources["source_verification"], built["bundle"],
                    str(root / name / "transfer"), case_id="case403")
                frozen = freeze_route_b_semantic_candidates(transfer["claims"], str(root / name / "freeze"), case_id="case403")
                wet = wet_candidate_handoff(str(root / name / "wet"), route_b_candidates_arg=frozen["candidates"],
                    route_b_candidates_manifest_arg=frozen["manifest"], review_candidates_arg=frozen["review_candidates"],
                    review_candidates_manifest_arg=frozen["review_manifest"], interval_atlas_arg=str(root / "atlas.tsv"), case_id="case403")
                self._assert_wet_warning(wet, {"G1", "G2"}, expected_warning=name == "conflict")
                runs.append(read_table(wet["wet_candidates"])[["candidate_gene_id", "relevance_tier", "support_state", "membership_basis"]])
                self.assertEqual(hashes, {path: hashlib.sha256(Path(path).read_bytes()).hexdigest() for path in hashes})
            pd.testing.assert_frame_equal(runs[0], runs[1])
            self.assertEqual(bundle_before, Path(built["bundle"]).read_bytes())
