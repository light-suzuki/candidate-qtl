"""Initial recording requests concern scientific facts, not optional metadata."""
import json
import tempfile
import unittest
from pathlib import Path

from tests import test_paper_evidence_table as fixtures

paper = fixtures.paper


class OptionalRecordingTests(unittest.TestCase):
    def export(self, report, **record):
        with tempfile.TemporaryDirectory() as tmp:
            result, table, packet = fixtures.PaperEvidenceTableTests().export(Path(tmp), [report], **record)
            return result, table.iloc[0], packet

    def test_optional_blanks_and_unknowns_do_not_request_more_research(self):
        for value in ("", "unknown"):
            p = dict(paper(annotation_supports_hypothesis="yes"),
                     source_excerpt=value, source_location=value, publication_type=value)
            result, row, packet = self.export(p)
            self.assertEqual([], packet["requests"])
            self.assertEqual("recorded", result["recording_state"])
            self.assertEqual("T5", row.applicable_tier)
            self.assertEqual("assessed", row.tier_state)
            self.assertEqual(value, row.source_excerpt)
            fields = json.loads(row.missing_fields_json if value == "" else row.unknown_fields_json)
            self.assertIn("source_excerpt", fields)

    def test_required_facts_and_checklist_still_request_completion(self):
        for key in ("source_title", "reported_gene", "reported_gene_id", "reported_species"):
            p = dict(paper(annotation_supports_hypothesis="yes"), source_location="")
            p[key] = ""
            _, _, packet = self.export(p)
            self.assertIn(key, packet["requests"][0]["records"][0]["missing_fields"])
            self.assertNotIn("source_location", packet["requests"][0]["records"][0]["missing_fields"])
        p = dict(paper(annotation_supports_hypothesis="yes", function_supports_hypothesis=""), source_excerpt="")
        _, _, packet = self.export(p)
        self.assertIn("function_supports_hypothesis", packet["requests"][0]["records"][0]["missing_fields"])

    def test_unavailable_gene_id_is_kept_without_retry_or_tier_penalty(self):
        p = dict(paper(phenotype_reported="yes", similar_cellular_or_developmental_phenotype="yes"),
                 reported_gene_id="unknown", unknown_reason="Paper and linked record checked; no stable ID available.", source_location="")
        _, row, packet = self.export(p)
        self.assertEqual("unknown", row.reported_gene_id)
        self.assertEqual("T2", row.applicable_tier)
        self.assertEqual([], packet["requests"])

    def test_legacy_and_explicit_strict_recording_keep_optional_requirements(self):
        for record in ({"candidate_admission_policy": "source_function_v1"}, {"verification_policy": "content_v1"}):
            _, _, packet = self.export(dict(paper(annotation_supports_hypothesis="yes"), source_location=""), **record)
            self.assertIn("source_location", packet["requests"][0]["records"][0]["missing_fields"])
        p = dict(paper(annotation_supports_hypothesis="yes"), source_location="")
        p.pop("paper_table_policy")
        _, _, packet = self.export(p)
        self.assertIn("source_location", packet["requests"][0]["records"][0]["missing_fields"])
