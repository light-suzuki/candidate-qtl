"""Copy-specific unknown explanations are transported, never invented."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from scripts import agentic_route_b_light as light
from tests.test_route_b_light import CHECKLIST_FIELDS, MAIN_CACHE, PAPER, PAPER_TITLE, fixture, run_pipeline
from tests.test_route_b_light_v2 import decision, observation


class CopyUnknownTransportTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def rows(self, assessments, **observation_fields):
        payload = {"contract": light.SCIENCE_CONTRACT,
                   "paper_observations": [observation(**observation_fields)],
                   "gene_decisions": []}
        for number, fields in enumerate(assessments, 1):
            gene = decision(f"g{number}", f"REF-{chr(64 + number)}")
            gene["paper_assessments"][0].update(copy.deepcopy(fields))
            payload["gene_decisions"].append(gene)
        path = self.root / "results.jsonl"
        raw = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        path.write_bytes(raw)
        rows, issues = light._science_rows(path)
        self.assertEqual([], issues)
        self.assertEqual(raw, path.read_bytes())
        return payload, rows

    def test_copy_fields_reach_only_their_own_paper_without_marking_other_copy_attempted(self):
        reason_map = {key: f"Original exact reason: {key}." for key in CHECKLIST_FIELDS}
        payload, rows = self.rows([{"unknown_reason": "Copy one identity remains unresolved.",
                                    "followup_state": "completed", "paper_checklist_reasons": reason_map}, {}])
        first, second = rows["g1"]["papers"][0], rows["g2"]["papers"][0]
        self.assertEqual("Copy one identity remains unresolved.", first["unknown_reason"])
        self.assertEqual("completed", first["followup_state"])
        self.assertEqual(reason_map, first["paper_checklist_reasons"])
        self.assertNotIn("unknown_reason", second)
        self.assertNotIn("followup_state", second)
        self.assertEqual(payload["gene_decisions"][0]["paper_assessments"][0]["checklist"], first["checklist"])
        self.assertEqual(("T2", "T2"), (rows["g1"]["tier"], rows["g2"]["tier"]))

    def test_checked_unknown_alias_is_literal_and_nonstring_aliases_do_not_explain_unknown(self):
        literal = "  Original checked limitation.\nSecond original line.  "
        _, rows = self.rows([{}], checked_unknown_reasons=literal, unknown_reason="")
        paper = rows["g1"]["papers"][0]
        self.assertEqual(literal, paper["unknown_reason"])
        provenance = paper["research_annotations"]["copy_science_field_transport"]
        self.assertEqual(literal, provenance["input_fields"]["observation.checked_unknown_reasons"])
        for value in ({"reason": "do not invent a string"}, ["reason"], "unknown", "   "):
            with self.subTest(value=value):
                _, rows = self.rows([{"checked_unknown_reasons": value}])
                self.assertNotIn("unknown_reason", rows["g1"]["papers"][0])

    def test_nested_checked_alias_is_preserved_and_general_copy_reason_is_not_repurposed(self):
        _, rows = self.rows([{"paper_checklist_reasons": {"checked_unknown_reasons": "Pea copy expression was not established."}}])
        self.assertEqual("Pea copy expression was not established.", rows["g1"]["papers"][0]["unknown_reason"])
        _, rows = self.rows([{"reason": "A general bridge is not an unknown explanation."}])
        self.assertNotIn("unknown_reason", rows["g1"]["papers"][0])

    def test_unfilled_canonical_placeholders_do_not_hide_recorded_literal_explanations(self):
        for placeholder in ("unknown", "  ", "不明"):
            for fields in ({"checked_unknown_reasons": "Explicit original checked reason."},
                           {"unknown_reason": "Explicit original checked reason."}):
                with self.subTest(placeholder=placeholder, fields=fields):
                    _, rows = self.rows([fields], unknown_reason=placeholder)
                    paper = rows["g1"]["papers"][0]
                    self.assertEqual("Explicit original checked reason.", paper["unknown_reason"])
                    provenance = paper["research_annotations"]["copy_science_field_transport"]
                    self.assertEqual(placeholder, provenance["input_fields"]["observation.unknown_reason"])
                    self.assertEqual([], provenance["conflicts"])

    def test_empty_canonical_reason_map_does_not_mask_an_explicit_alias_map(self):
        checks = decision()["paper_assessments"][0]["checklist"]
        reasons = {key: "Original explicit reason " + key for key, value in checks.items() if value == "unknown"}
        _, rows = self.rows([{"paper_checklist_reasons": {}, "checklist_reasons": reasons}])
        paper = rows["g1"]["papers"][0]
        self.assertEqual(reasons, paper["paper_checklist_reasons"])
        self.assertEqual(reasons, json.loads(paper["unknown_reason"]))
        self.assertEqual({}, paper["research_annotations"]["copy_science_field_transport"]["input_fields"]["assessment.paper_checklist_reasons"])

    def test_full_unknown_answer_reason_map_keeps_exact_keys_and_strings(self):
        checks = decision()["paper_assessments"][0]["checklist"]
        reasons = {key: "  Exact reason for " + key + ".  " for key, answer in checks.items() if answer == "unknown"}
        _, rows = self.rows([{"paper_checklist_reasons": reasons}])
        paper = rows["g1"]["papers"][0]
        self.assertEqual(reasons, json.loads(paper["unknown_reason"]))
        self.assertEqual(checks, paper["checklist"])
        self.assertNotIn("followup_state", paper)
        for incomplete in ({key: value for key, value in reasons.items() if key != next(iter(reasons))},
                           {key: "unknown" for key in reasons}):
            _, rows = self.rows([{"paper_checklist_reasons": incomplete}])
            self.assertNotIn("unknown_reason", rows["g1"]["papers"][0])

    def test_checklist_reason_map_does_not_cover_an_unexplained_unknown_reported_identifier(self):
        reasons = {key: "Explicit original answer reason." for key in CHECKLIST_FIELDS}
        _, rows = self.rows([{"paper_checklist_reasons": reasons}], reported_gene_id="unknown")
        self.assertNotIn("unknown_reason", rows["g1"]["papers"][0])
        reasons["reported_gene_id"] = "The inspected source supplies the name but no identifier."
        _, rows = self.rows([{"paper_checklist_reasons": reasons}], reported_gene_id="unknown")
        self.assertEqual(reasons["reported_gene_id"], json.loads(rows["g1"]["papers"][0]["unknown_reason"])["reported_gene_id"])

    def test_conflicting_canonical_values_are_preserved_with_both_original_values(self):
        _, rows = self.rows([{"unknown_reason": "Copy-specific different explanation.", "followup_state": "completed",
                              "paper_checklist_reasons": {"focal_species": "Copy reason."}}],
                            unknown_reason="Original observation explanation.", followup_state="not_attempted",
                            paper_checklist_reasons={"focal_species": "Original canonical reason."})
        paper = rows["g1"]["papers"][0]
        self.assertEqual("Original observation explanation.", paper["unknown_reason"])
        self.assertEqual("not_attempted", paper["followup_state"])
        self.assertEqual({"focal_species": "Original canonical reason."}, paper["paper_checklist_reasons"])
        provenance = paper["research_annotations"]["copy_science_field_transport"]
        self.assertEqual("conflicting_original_values", provenance["state"])
        self.assertEqual({"unknown_reason", "followup_state", "paper_checklist_reasons"}, {item["field"] for item in provenance["conflicts"]})
        self.assertEqual("Copy-specific different explanation.", next(item["incoming_value"] for item in provenance["conflicts"] if item["field"] == "unknown_reason"))
        self.assertEqual("T2", rows["g1"]["tier"])

    def test_converter_import_and_freeze_keep_judgments_and_eliminate_only_recorded_unknown_blank(self):
        payload, _ = self.rows([{"unknown_reason": "The focal pea copy has not been experimentally resolved.",
                                 "followup_state": "completed",
                                 "paper_checklist_reasons": {key: "Original reason " + key for key in CHECKLIST_FIELDS}}])
        original = copy.deepcopy(payload)
        built = run_pipeline(self.root, fixture(), payload, MAIN_CACHE, "u1\n", output_contract="science-v2")
        proposals = pd.read_csv(built["imported"]["proposals"], sep="\t", dtype=str).fillna("")
        report = json.loads(proposals.iloc[0]["paper_reports_json"])[0]
        self.assertEqual("The focal pea copy has not been experimentally resolved.", report["unknown_reason"])
        self.assertEqual("completed", report["followup_state"])
        self.assertEqual(original["gene_decisions"][0]["paper_assessments"][0]["paper_checklist_reasons"], report["research_annotations"]["paper_checklist_reasons"])
        self.assertEqual((PAPER, PAPER_TITLE, "HOM1", "TAIR:AT1G12345", "Arabidopsis thaliana"),
                         tuple(report[key] for key in ("source_id", "source_title", "reported_gene", "reported_gene_id", "reported_species")))
        self.assertEqual(original["gene_decisions"][0]["paper_assessments"][0]["checklist"], json.loads(report["paper_checklist"]))
        self.assertTrue(built["review"].loc["g1", "relevance_tier"].startswith("T2"))
        self.assertEqual("unverified_route_b_review", built["review"].loc["g1", "support_status"])
        requests = json.loads((self.root / "import" / "paper_followup_requests.json").read_text(encoding="utf-8"))["requests"]
        self.assertEqual([], requests)
        self.assertEqual(original, payload)


if __name__ == "__main__":
    unittest.main()
