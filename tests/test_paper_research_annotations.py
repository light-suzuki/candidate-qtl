"""Research notes survive the actual A/B import without becoming paper facts."""
from __future__ import annotations

import json
import math
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.paper_reports import PAPER_FIELD, imported_paper_fields, normalize_paper_reports, select_bibliography_papers
from candidate_qtl.relevance_tiers import annotate_relevance
from candidate_qtl.route_b_response_validation import route_b_claim_schema_error
from candidate_qtl.route_b_semantic_import import import_route_b_semantic_proposals
from candidate_qtl.semantic_claims import narrowing_eligible_claims
from tests.test_initial_paper_prompt_examples import _filled_hypothesis
from tests.test_paper_evidence_table import paper
from tests.test_paper_tiers import claim as relevance_claim
from tests import test_route_b_semantic_import as b_fixtures


def _notes() -> dict:
    return {
        "unit_id": "synthetic-unit-1",
        "counterevidence": ["The experiment did not resolve the target copy.", {"state": "unknown", "value": None}],
        "stage": "initial_research",
        "read_scope": "abstract_and_results",
        "conditional_focal_reasoning": "Homolog transfer remains conditional.",
        "source_url": "https://example.invalid/synthetic-paper",
        "reference_species": "Plant alpha",
        "formal_claim_state": "proposed",
        "retrieved_at": "2026-01-02T03:04:05Z",
        "full_paper_read": False,
        "paper_checklist_reasons": {"phenotype_reported": "The source describes the experiment."},
    }


def _report(**extra) -> dict:
    return {**paper(phenotype_reported="yes", similar_cellular_or_developmental_phenotype="yes"), **extra}


def _claim(report: dict) -> dict:
    return {"claim_type": "mutant_phenotype", "reference_gene": "REF1",
            "reference_species": "Plant beta", "evidence_scope": "reference_gene",
            "source_ids": ["PMID:123"], PAPER_FIELD: [report]}


def _read_table(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


class PaperResearchAnnotationsTests(unittest.TestCase):
    def test_legacy_research_fields_move_losslessly_and_do_not_mutate_the_input(self):
        report = _report(**_notes())
        before = deepcopy(report)
        encoded = normalize_paper_reports([report])
        normalized = json.loads(encoded)[0]
        self.assertEqual(_notes(), normalized["research_annotations"])
        self.assertTrue(set(_notes()).isdisjoint(normalized))
        self.assertEqual(before, report)
        self.assertEqual(encoded, normalize_paper_reports(encoded))
        explicit = _report(research_annotations=_notes())
        self.assertEqual(encoded, normalize_paper_reports([explicit]))

    def test_old_normalized_bytes_have_no_new_keys_or_values(self):
        original = normalize_paper_reports([_report()])
        annotated = json.loads(normalize_paper_reports([_report(research_annotations=_notes())]))
        self.assertNotIn("research_annotations", json.loads(original)[0])
        del annotated[0]["research_annotations"]
        # The old field order and compact serialization are part of the frozen artifact contract.
        self.assertEqual(original, json.dumps(annotated, ensure_ascii=False, separators=(",", ":")))
        self.assertEqual(original, normalize_paper_reports(original))

    def test_envelope_is_deterministic_and_preserves_nested_json_values(self):
        notes = {"z": [True, False, None, 0, 1.5, {"text": " checked unknown "}], "a": {"b": "理由"}}
        first = normalize_paper_reports([_report(research_annotations=notes)])
        second = normalize_paper_reports([_report(research_annotations=dict(reversed(list(notes.items()))))])
        self.assertEqual(first, second)
        self.assertEqual(notes, json.loads(first)[0]["research_annotations"])

    def test_metadata_conflicts_reject_instead_of_overwriting_either_observation(self):
        report = _report(counterevidence=["original"], research_annotations={"counterevidence": ["different"]})
        with self.assertRaises(ValueError):
            normalize_paper_reports([report])
        report["research_annotations"]["counterevidence"] = ["original"]
        self.assertEqual({"counterevidence": ["original"]}, json.loads(normalize_paper_reports([report]))[0]["research_annotations"])

    def test_only_json_objects_and_plain_finite_nested_values_are_annotations(self):
        bad = ["text", [], None, {"nested": {"not", "json"}}, {1: "non-string key"},
               {"nested": {1: "non-string key"}}, {"nested": math.nan}, {"nested": math.inf}]
        for notes in bad:
            with self.subTest(notes=notes), self.assertRaises((ValueError, TypeError)):
                normalize_paper_reports([_report(research_annotations=notes)])

    def test_annotation_transport_does_not_relax_science_or_source_binding(self):
        invalid = [
            {"reported_gene": ["REF1"]}, {"reported_species": True},
            {"reported_genne": "REF1"}, {"paper_checklist": {"not_a_question": "yes"}},
            {"paper_checklist": {"phenotype_reported": "perhaps"}},
            {"source_id": "not-a-citation"}, {"source_identifier_aliases": ["UniProt:P12345"]},
        ]
        for changes in invalid:
            report = _report(research_annotations=_notes(), **changes)
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                normalize_paper_reports([report])
            self.assertTrue(route_b_claim_schema_error(_claim(report)))
        unbound = _claim(_report(research_annotations=_notes()))
        unbound["source_ids"] = ["PMID:999"]
        with self.assertRaisesRegex(ValueError, "not bound"):
            imported_paper_fields(unbound)
        self.assertTrue(route_b_claim_schema_error(unbound))

    def test_actual_a_cli_preserves_annotations_and_original_paper_handoff(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            frames = {}
            for label, report in (("plain", _report()), ("annotated", _report(**_notes()))):
                payload = root / f"{label}.json"
                hypothesis = _filled_hypothesis({}, "SYNTHETIC_FAMILY", [report])
                payload.write_text(json.dumps({"hypotheses": [hypothesis]}), encoding="utf-8")
                out = root / label
                self.assertEqual(0, main(["mechanism-import", "--payload", str(payload), "--out-dir", str(out)]))
                frames[label] = _read_table(out / "mechanism_hypotheses.tsv")
                table = _read_table(out / "paper_evidence_table.tsv")
                self.assertEqual(["T2"], table.computed_tier.tolist())
                self.assertEqual(["yes"], table.selected_for_handoff.tolist())
                self.assertEqual(["primary"], table.paper_handoff_state.tolist())
                if label == "annotated":
                    self.assertEqual(_notes(), json.loads(table.iloc[0].research_annotations_json))
                    self.assertEqual(_notes(), json.loads(frames[label].iloc[0][PAPER_FIELD])[0]["research_annotations"])
                else:
                    self.assertNotIn("research_annotations_json", table.columns)
            pd.testing.assert_frame_equal(frames["plain"].drop(columns=[PAPER_FIELD]), frames["annotated"].drop(columns=[PAPER_FIELD]))

    def test_actual_b_import_keeps_membership_and_support_separate_from_notes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = b_fixtures.RouteBSemanticImportTests()
            built = fixture._bundle(root)
            frames = {}
            for label, report in (("plain", _report()), ("annotated", _report(**_notes()))):
                response = root / f"{label}.json"
                fixture._response(Path(built["prompts"]), [_claim(report)], response)
                result = import_route_b_semantic_proposals(built["bundle"], built["prompts"], str(response), str(root / label))
                frames[label] = _read_table(result["proposals"])
                self.assertEqual((1, 0, 0), (result["proposed_claim_row_count"], result["verified_claim_row_count"], result["rejected_claim_row_count"]))
                self.assertEqual(["proposed"], frames[label].verification_status.tolist())
                self.assertTrue(narrowing_eligible_claims(frames[label]).empty)
                self.assertEqual({"G1", "G2"}, set(frames[label].iloc[0].member_interval_genes.split(";")))
                table = _read_table(result["paper_evidence_table"])
                self.assertEqual(["T2"], table.computed_tier.tolist())
                if label == "annotated":
                    self.assertEqual(_notes(), json.loads(frames[label].iloc[0][PAPER_FIELD])[0]["research_annotations"])
                    self.assertEqual(_notes(), json.loads(table.iloc[0].research_annotations_json))
                else:
                    self.assertNotIn("research_annotations_json", table.columns)
            ignored = [PAPER_FIELD, "output_fingerprint"]
            pd.testing.assert_frame_equal(frames["plain"].drop(columns=ignored), frames["annotated"].drop(columns=ignored))

    def test_invalid_paper_is_quarantined_without_function_only_fallback(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = b_fixtures.RouteBSemanticImportTests()
            built = fixture._bundle(root)
            bad = _claim(_report(reported_gene=["invalid"], research_annotations=_notes()))
            bad.update(biological_function="A valid-looking function must not mask an invalid paper report.")
            response = root / "response.json"
            fixture._response(Path(built["prompts"]), [_claim(_report()), bad], response)
            result = import_route_b_semantic_proposals(built["bundle"], built["prompts"], str(response), str(root / "import"))
            self.assertEqual((1, 1), (result["claim_row_count"], result["rejected_claim_row_count"]))
            self.assertEqual(0, result["verified_claim_row_count"])
            rejected = json.loads(Path(result["rejected_rows"]).read_text(encoding="utf-8"))
            self.assertEqual(bad, rejected[0]["raw_claim"])
            self.assertEqual(["partial_import"], _read_table(result["status"]).status.tolist())

    def test_summary_retains_all_notes_without_changing_tier_or_support(self):
        roster = [{"candidate_gene_id": "G1", "support_status": "unverified_route_b_review"}]
        first = _report(research_annotations={"counterevidence": ["copy unresolved"]})
        second = _report(research_annotations={"counterevidence": ["tissue unresolved"]})
        baseline = annotate_relevance(roster, [relevance_claim([_report()])])[0]
        annotated = annotate_relevance(roster, [relevance_claim([first, second])])[0]
        for field in ("relevance_tier", "initial_tier", "support_state", "tier_source_id", "paper_tier_counts_json"):
            self.assertEqual(baseline[field], annotated[field], field)
        summary = json.loads(annotated["paper_evidence_summary_json"])
        self.assertEqual(1, len(summary))
        self.assertEqual({json.dumps(p["research_annotations"], sort_keys=True) for p in (first, second)},
                         {json.dumps(notes, sort_keys=True) for notes in summary[0]["research_annotations"]})

    def test_annotation_lexical_order_cannot_change_equal_tier_five_paper_handoff(self):
        reports = [_report(source_id=f"PMID:{number}") for number in range(101, 107)]
        citations = [report["source_id"] for report in reports]
        baseline = select_bibliography_papers({PAPER_FIELD: reports}, citations)
        annotated = [{**report, "research_annotations": {"unit_id": f"note-{106 - index}"}}
                     for index, report in enumerate(reports)]
        with_notes = select_bibliography_papers({PAPER_FIELD: annotated}, citations)
        self.assertEqual(5, baseline["selected_paper_count"])
        self.assertEqual(6, baseline["candidate_paper_count"])
        self.assertEqual(baseline, with_notes)
        reversed_notes = [{**report, "research_annotations": {"unit_id": f"note-{index}"}}
                          for index, report in enumerate(reports)]
        self.assertEqual(baseline, select_bibliography_papers({PAPER_FIELD: list(reversed(reversed_notes))}, citations))

    def test_nested_scientific_names_remain_notes_and_cannot_authorize_support(self):
        notes = {"verification_status": "verified", "reference_gene_id": "OTHER",
                 "source_id": "PMID:999", "initial_tier": "T1",
                 "paper_checklist": {"phenotype_reported": "yes", "focal_species": "yes",
                                     "same_gene_as_focal": "yes", "equivalent_target_phenotype": "yes"}}
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = b_fixtures.RouteBSemanticImportTests()
            built = fixture._bundle(root)
            response = root / "response.json"
            fixture._response(Path(built["prompts"]), [_claim(_report(research_annotations=notes))], response)
            result = import_route_b_semantic_proposals(built["bundle"], built["prompts"], str(response), str(root / "import"))
            rows = _read_table(result["proposals"])
            self.assertEqual((1, 0), (result["proposed_claim_row_count"], result["verified_claim_row_count"]))
            self.assertTrue(narrowing_eligible_claims(rows).empty)
            row = rows.iloc[0]
            self.assertEqual(("T2", "proposed", "REF1", "PMID:123"),
                             (row.initial_tier, row.verification_status, row.reference_gene, row.source_id))
            self.assertEqual({"G1", "G2"}, set(row.member_interval_genes.split(";")))
            report = json.loads(row[PAPER_FIELD])[0]
            self.assertEqual(notes, report["research_annotations"])
            self.assertEqual("REF1", report["reference_gene_id"])
            self.assertEqual("no", json.loads(report["paper_checklist"])["focal_species"])
            unbound = _claim(_report(research_annotations=notes))
            unbound["source_ids"] = ["PMID:999"]
            self.assertTrue(route_b_claim_schema_error(unbound))

    def test_annotations_do_not_choose_a_different_scientific_representative(self):
        reports = [_report(source_excerpt="a"), _report(source_excerpt="z")]
        annotated = [{**report, "research_annotations": {"counterevidence": note}}
                     for report, note in zip(reports, ("z", "a"))]
        plain_record = {"source_ids": ["PMID:123"], PAPER_FIELD: reports}
        annotated_record = {"source_ids": ["PMID:123"], PAPER_FIELD: annotated}
        baseline = imported_paper_fields(plain_record, "PMID:123")
        with_notes = imported_paper_fields(annotated_record, "PMID:123")
        self.assertEqual("a", baseline["source_excerpt"])
        self.assertEqual(baseline["source_excerpt"], with_notes["source_excerpt"])
        for key in baseline.keys() - {PAPER_FIELD}:
            self.assertEqual(baseline[key], with_notes[key], key)
        roster = [{"candidate_gene_id": "G1", "support_status": "unverified_route_b_review"}]
        plain_tier = annotate_relevance(roster, [relevance_claim(reports)])[0]
        annotated_tier = annotate_relevance(roster, [relevance_claim(annotated)])[0]
        plain_summary = json.loads(plain_tier["paper_evidence_summary_json"])
        annotated_summary = json.loads(annotated_tier["paper_evidence_summary_json"])
        self.assertEqual(1, len(plain_summary))
        self.assertEqual(1, len(annotated_summary))
        for key in plain_summary[0]:
            self.assertEqual(plain_summary[0][key], annotated_summary[0][key], key)
        for key in ("relevance_tier", "initial_tier", "support_state", "tier_source_id", "tier_evidence_json"):
            self.assertEqual(plain_tier[key], annotated_tier[key], key)


if __name__ == "__main__":
    unittest.main()
