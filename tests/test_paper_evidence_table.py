"""Paper recording completeness is separate from biological Tier and membership."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.observation_tiers import CHECKLIST_FIELDS
from candidate_qtl.paper_evidence_table import write_paper_evidence_table
from candidate_qtl.paper_reports import PAPER_FIELD, PAPER_TABLE_POLICY, normalize_paper_reports
from candidate_qtl.relevance_tiers import classify_relevance
from candidate_qtl.route_b_semantic_import import import_route_b_semantic_proposals
from tests import test_paper_tiers as a
from tests import test_route_a_review_relevance as wet
from tests import test_route_b_semantic_import as b


def paper(**checks):
    return {"paper_table_policy": PAPER_TABLE_POLICY,
            "source_id": "PMID:123", "source_title": "Synthetic experiment", "publication_type": "original",
            "reported_gene": "GENE_A", "reported_gene_id": "NCBI Gene:12345", "reported_species": "Plant beta",
            "reference_gene_id": "REF1", "gene_correspondence": "same_gene",
            "source_excerpt": "The mutation changes meristem size.", "source_location": "Figure 3",
            "initial_tier_reason": "The paper reports a similar developmental phenotype.",
            "unknown_reason": "", "followup_state": "not_attempted",
            "paper_checklist": {**dict.fromkeys(CHECKLIST_FIELDS, "no"), "identity_resolved": "yes", **checks}}


class PaperEvidenceTableTests(unittest.TestCase):
    def export(self, root, reports, **record):
        source = root / "raw.json"
        source.write_text(json.dumps(reports))
        data = [{"hypothesis_id": "h", "reference_gene_name": "REF1", "source_ids": "PMID:123",
                 PAPER_FIELD: reports, **record}]
        result = write_paper_evidence_table(data, root / "table", route="phenotype_first", source_artifact=str(source))
        table = pd.read_csv(result["outputs"]["paper_evidence_table"], sep="\t", dtype=str).fillna("")
        packet = json.loads(Path(result["outputs"]["paper_followup_requests"]).read_text())
        return result, table, packet

    def test_blank_and_checked_unknown_survive_normalization_as_distinct_states(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            blank = paper(function_supports_hypothesis=" ", annotation_supports_hypothesis="yes")
            encoded = normalize_paper_reports([blank])
            self.assertEqual(encoded, normalize_paper_reports(encoded))
            _, table, packet = self.export(root, encoded)
            self.assertEqual("", table.iloc[0].function_supports_hypothesis)
            self.assertEqual("T5", table.iloc[0].computed_tier)
            self.assertEqual("provisional_missing_record", table.iloc[0].tier_state)
            self.assertIn("function_supports_hypothesis", packet["requests"][0]["records"][0]["missing_fields"])
            checked = paper(function_supports_hypothesis="unknown", annotation_supports_hypothesis="yes")
            checked["unknown_reason"] = "Abstract and results checked; function remains unresolved."
            result, table, packet = self.export(root, [checked])
            self.assertEqual("checked_unknown", table.iloc[0].record_state)
            self.assertEqual("assessed", table.iloc[0].tier_state)
            self.assertEqual("recorded", result["recording_state"])
            self.assertEqual([], packet["requests"])

    def test_only_concrete_upper_tier_leads_get_one_bounded_followup(self):
        p = paper(phenotype_reported="unknown", similar_cellular_or_developmental_phenotype="yes",
                  annotation_supports_hypothesis="yes")
        p["unknown_reason"] = "The available passage suggests a phenotype; experiment detail was unavailable."
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, table, packet = self.export(root, [p])
            self.assertEqual("T5", table.iloc[0].computed_tier)
            self.assertEqual("clarify_upper_tier_lead", table.iloc[0].followup_action)
            self.assertEqual(["phenotype_reported"], packet["requests"][0]["records"][0]["questions"])
            self.assertEqual(1, packet["max_followup_batches"])
            self.assertEqual("not_executed", packet["execution_state"])
            p["followup_state"] = "completed"
            _, table, packet = self.export(root, [p])
            self.assertEqual("T5", table.iloc[0].computed_tier)
            self.assertEqual([], packet["requests"])

    def test_corrected_normalized_report_can_close_missingness_without_inventing_facts(self):
        original = paper(function_supports_hypothesis="", annotation_supports_hypothesis="yes")
        normalized = json.loads(normalize_paper_reports([original]))[0]
        self.assertIn("function_supports_hypothesis", normalized["paper_checklist_missing_fields_json"])
        for answer in ("yes", "unknown"):
            corrected = json.loads(json.dumps(normalized))
            corrected["paper_checklist"] = json.loads(corrected["paper_checklist"])
            corrected["paper_checklist"]["function_supports_hypothesis"] = answer
            corrected.update(followup_state="completed", unknown_reason="Results checked; unresolved function is retained when unknown.")
            with self.subTest(answer=answer), tempfile.TemporaryDirectory() as tmp:
                _, table, packet = self.export(Path(tmp), [corrected])
                self.assertEqual("[]", table.iloc[0].missing_fields_json)
                self.assertEqual("T4" if answer == "yes" else "T5", table.iloc[0].computed_tier)
                self.assertEqual([], packet["requests"])

    def test_confirmed_t2_does_not_chase_impossible_t1_and_missing_gene_id_does_not_lower_tier(self):
        p = paper(phenotype_reported="yes", similar_cellular_or_developmental_phenotype="yes",
                  equivalent_target_phenotype="unknown")
        p.update(reported_gene_id="unknown", unknown_reason="No stable ID is given; equivalent phenotype is uncertain.")
        with tempfile.TemporaryDirectory() as tmp:
            _, table, packet = self.export(Path(tmp), [p])
            self.assertEqual("T2", table.iloc[0].applicable_tier)
            self.assertEqual([], packet["requests"])

    def test_unfilled_after_completed_attempt_remains_partial_without_automatic_retry(self):
        p = dict(paper(annotation_supports_hypothesis="yes"), source_title="", followup_state="completed")
        with tempfile.TemporaryDirectory() as tmp:
            result, table, packet = self.export(Path(tmp), [p])
            self.assertEqual("needs_completion", result["recording_state"])
            self.assertEqual("unfilled_after_followup", table.iloc[0].followup_action)
            self.assertEqual("followup_exhausted", packet["execution_state"])
            self.assertEqual([], packet["requests"])

    def test_repeated_exact_reference_checks_reuse_work_without_sharing_paper_facts(self):
        first = paper(annotation_supports_hypothesis="yes")
        second = {**paper(phenotype_reported="yes", similar_cellular_or_developmental_phenotype="yes"),
                  "source_id": "PMID:456"}
        with tempfile.TemporaryDirectory() as tmp, patch("candidate_qtl.paper_evidence_table.gene_binding", return_value="matched") as bind:
            result, table, _ = self.export(Path(tmp), [first, second], source_ids="PMID:123;PMID:456")
            self.assertEqual(1, bind.call_count)
            self.assertEqual(1, result["reference_anchor_check_count"])
            self.assertEqual(1, table.correspondence_key.nunique())
            self.assertEqual({"T2", "T5"}, set(table.computed_tier))
            self.assertEqual("T2", table.iloc[0].computed_tier)

    def test_metadata_panel_keeps_stronger_paper_ahead_of_six_annotation_papers(self):
        reports = [{**paper(annotation_supports_hypothesis="yes"), "source_id": f"PMID:{i}"} for i in range(101, 107)]
        reports.append(paper(phenotype_reported="yes", similar_cellular_or_developmental_phenotype="yes"))
        with tempfile.TemporaryDirectory() as tmp:
            _, table, _ = self.export(Path(tmp), reports, source_ids=";".join(p["source_id"] for p in reports))
            self.assertEqual(7, len(table))
            self.assertEqual(5, (table.selected_for_metadata == "yes").sum())
            self.assertEqual("yes", table[table.computed_tier == "T2"].iloc[0].selected_for_metadata)

    def test_cited_paper_without_report_is_not_annotation_only_research_completion(self):
        with tempfile.TemporaryDirectory() as tmp:
            result, table, packet = self.export(Path(tmp), [])
            self.assertEqual("needs_completion", result["recording_state"])
            self.assertEqual("missing_paper_report", table.iloc[0].record_kind)
            self.assertEqual("U", table.iloc[0].applicable_tier)
            self.assertTrue(packet["requests"])

    def test_blank_gene_correspondence_is_returned_without_authorizing_a_high_tier(self):
        for key in ("reported_gene", "reference_gene_id", "gene_correspondence", "initial_tier_reason"):
            p = paper(phenotype_reported="yes", equivalent_target_phenotype="yes", focal_species="yes", same_gene_as_focal="yes")
            p[key] = ""
            with self.subTest(key=key), tempfile.TemporaryDirectory() as tmp:
                result, table, packet = self.export(Path(tmp), [p])
                self.assertEqual("needs_completion", result["recording_state"])
                self.assertIn(key, packet["requests"][0]["records"][0]["missing_fields"])
                if key != "initial_tier_reason":
                    self.assertEqual("U", table.iloc[0].applicable_tier)
                    self.assertEqual("U_unresolved", classify_relevance(a.claim([p]))[0])

    def test_real_a_cli_exports_same_facts_and_targeted_recording_request(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            p = paper(phenotype_reported="yes", similar_cellular_or_developmental_phenotype="yes")
            p["source_title"] = ""
            payload = {**a.MECHANISMS[0], "reference_gene_name": "REF1", "source_ids": "PMID:123", PAPER_FIELD: [p]}
            self.assertEqual(0, main(["mechanism-import", "--payload", str(a._payload(root, [payload])), "--out-dir", str(root / "import")]))
            manifest = json.loads((root / "import/mechanism_import_manifest.json").read_text())
            table = pd.read_csv(manifest["outputs"]["paper_evidence_table"], sep="\t").fillna("")
            self.assertEqual("T2", table.iloc[0].computed_tier)
            self.assertEqual("not_available_before_projection", table.iloc[0].projection_state)
            self.assertEqual("needs_completion", manifest["paper_evidence"]["recording_state"])
            self.assertEqual(1, manifest["hypothesis_count"])

    def test_real_b_import_exports_frozen_members_and_retains_proposed_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = b.RouteBSemanticImportTests()
            built = fixture._bundle(root)
            claim = {"claim_type": "mutant_phenotype", "reference_gene": "REF1", "evidence_scope": "reference_gene",
                     "source_ids": ["PMID:123"], PAPER_FIELD: [paper(phenotype_reported="yes", similar_cellular_or_developmental_phenotype="yes")]}
            source = root / "response.json"
            fixture._response(Path(built["prompts"]), [claim], source)
            result = import_route_b_semantic_proposals(built["bundle"], built["prompts"], str(source), str(root / "import"))
            table = pd.read_csv(result["paper_evidence_table"], sep="\t").fillna("")
            self.assertEqual({"T2"}, set(table.computed_tier))
            self.assertTrue(table.member_interval_genes.str.len().gt(0).all())
            self.assertEqual(1, result["proposed_claim_row_count"])
            self.assertEqual(0, result["verified_claim_row_count"])

    def test_metadata_results_join_into_table_without_new_lookups(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            p = paper(phenotype_reported="yes", similar_cellular_or_developmental_phenotype="yes")
            source = root / "input.tsv"
            pd.DataFrame([{"hypothesis_id": "h", "reference_gene_name": "REF1", "source_ids": "PMID:123", PAPER_FIELD: json.dumps([p])}]).to_csv(source, sep="\t", index=False)
            with patch("candidate_qtl.bibliography.fetch_metadata", return_value={"pmid": "123", "title": p["source_title"], "doi": "10.1234/example"}) as fetch, \
                 patch("candidate_qtl.gene_metadata.fetch_gene_metadata", return_value={"state": "not_found"}):
                self.assertEqual(0, main(["source-verify", "--mechanisms", str(source), "--out-dir", str(root / "verified"), "--allow-external-api"]))
            self.assertEqual(1, fetch.call_count)
            table = pd.read_csv(root / "verified/paper_evidence_table.tsv", sep="\t").fillna("")
            self.assertEqual("id_title_matched", table.iloc[0].metadata_state)
            self.assertEqual(p["source_title"], table.iloc[0].resolved_title)
            self.assertEqual("10.1234/example", table.iloc[0].resolved_doi)

    def test_provisional_recording_never_vetoes_existing_a_wet_members(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            p = paper(phenotype_reported="yes", similar_cellular_or_developmental_phenotype="yes")
            p["source_location"] = ""
            result, _, _, _ = wet.RouteAReviewRelevanceTests()._run(root, {"source_ids": "PMID:123", PAPER_FIELD: normalize_paper_reports([p])})
            rows = pd.read_csv(result["wet_candidates"], sep="\t").fillna("")
            self.assertTrue(len(rows))
            self.assertTrue(rows.relevance_tier.str.startswith("T2_").all())
            self.assertTrue(all(json.loads(value)[0]["source_location"] == "" for value in rows.paper_evidence_summary_json))
