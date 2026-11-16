"""Preserve original scientific reasons and structured locations through canonical transport."""
import json
import tempfile
import unittest
from pathlib import Path
import pandas as pd
from tests.test_route_b_light import (light, CHECKLIST_FIELDS, PAPER, PAPER_TITLE, MAIN_CACHE, paper_row, gene_row, checklist, run_pipeline, fixture)


class LightChecklistReasonTransportTests(unittest.TestCase):
    def test_empty_structured_location_remains_missing_in_strict_recording(self):
        from candidate_qtl.paper_evidence_table import write_paper_evidence_table

        for location in ({}, [], ""):
            with self.subTest(location=location), tempfile.TemporaryDirectory() as temporary:
                paper = paper_row(PAPER, PAPER_TITLE, checklist("T2"), reported_gene_id="TAIR:AT1G12345",
                                  publication_type="original", source_excerpt="Original observed result.",
                                  unknown_reason="The focal copy has not been experimentally resolved.",
                                  source_location=location)
                original = json.dumps(paper, sort_keys=True)
                kept, state, _, _ = light._check_paper(paper, MAIN_CACHE.copy(), True, None, [PAPER], "Arabidopsis thaliana")
                self.assertEqual("checked", state)
                report = light._paper_report(kept, "REF-A")
                self.assertEqual("", report["source_location"])
                if not isinstance(location, str):
                    self.assertEqual(location, report["research_annotations"]["source_location_field_provenance"]["input_value"])
                out = Path(temporary) / "strict"
                records = [{"claim_id": "c1", "candidate_admission_policy": "source_function_v1",
                    "source_id": PAPER, "reference_gene_id": "REF-A", "reference_species": "Arabidopsis thaliana",
                    "paper_reports_json": [report]}]
                source = Path(temporary) / "original.json"
                source.write_text(json.dumps(records), encoding="utf-8")
                write_paper_evidence_table(records, out, route="interval_first", source_artifact=str(source))
                row = pd.read_csv(out / "paper_evidence_table.tsv", sep="\t", dtype=str).fillna("").iloc[0]
                self.assertIn("source_location", json.loads(row["missing_fields_json"]))
                self.assertEqual("complete_record", row["followup_action"])
                self.assertEqual(original, json.dumps(paper, sort_keys=True))

    def test_nonfinite_structured_location_is_a_local_invalid_paper(self):
        for location in ({"page": float("nan")}, [float("inf")]):
            paper = paper_row(PAPER, PAPER_TITLE, checklist("T2"), source_location=location)
            kept, state, detail, paper_id = light._check_paper(
                paper, MAIN_CACHE.copy(), True, None, [PAPER], "Arabidopsis thaliana")
            self.assertIsNone(kept)
            self.assertEqual("invalid_record", state)
            self.assertIn("source_location", detail)
            self.assertEqual(PAPER, paper_id)

    def test_structured_source_location_survives_converter_import_and_freeze(self):
        reasons = {key: f"Original reason for {key}." for key in CHECKLIST_FIELDS}
        locations = [{"section": "結果", "figure": {"number": 2, "panels": ["a", "b"]}, "page": 3},
                     ["Abstract", {"figure": 2, "panel": "b"}], {}, [], "  Fig. 2\n Panel B  "]
        for contract in ("legacy-v1", "science-v2"):
            for location in locations:
                with self.subTest(contract=contract, location=location), tempfile.TemporaryDirectory() as temporary:
                    root = Path(temporary)
                    paper = paper_row(PAPER, PAPER_TITLE, checklist("T2"), reported_gene_id="TAIR:AT1G12345",
                                      source_location=location, checklist_reasons=reasons,
                                      research_annotations={"source_location_field_provenance": {"note": "Original context."}})
                    rows = [gene_row("g1", "REF-A", "T2", [paper])]
                    if contract == "science-v2":
                        rows = {"contract": light.SCIENCE_CONTRACT,
                                "paper_observations": [{"observation_id": "o1", **{
                                    key: value for key, value in paper.items() if key not in {"checklist", "gene_correspondence"}}}],
                                "gene_decisions": [{"interval_gene_id": "g1", "reference_gene": "REF-A",
                                    "state": "assessed", "basis": "paper_evidence", "function": "leaf development",
                                    "reason": "Original gene bridge.", "paper_assessments": [{"observation_id": "o1",
                                        "gene_correspondence": "homolog", "checklist": paper["checklist"],
                                        "reason": "Original copy-specific application."}]}]}
                    original = json.dumps(rows, sort_keys=True)
                    expected_bytes = (json.dumps(rows) if contract == "science-v2" else
                                      "".join(json.dumps(row) + "\n" for row in rows)).encode("utf-8")
                    built = run_pipeline(root, fixture(), rows, MAIN_CACHE, "u1\n", output_contract=contract)
                    self.assertIn("g1", built["review"].index)
                    proposals = pd.read_csv(built["imported"]["proposals"], sep="\t", dtype=str).fillna("")
                    report = json.loads(proposals.iloc[0]["paper_reports_json"])[0]
                    annotations = report["research_annotations"]
                    if isinstance(location, str):
                        response = json.loads((root / "converted/agent_responses.json").read_text())["responses"][0]
                        original_report = json.loads(response["raw_text"])["claims"][0]["paper_reports_json"][0]
                        self.assertEqual(location, original_report["source_location"])
                        # Existing canonical import strips outer whitespace; the bridge does not.
                        self.assertEqual(location.strip(), report["source_location"])
                        self.assertEqual({"note": "Original context."}, annotations["source_location_field_provenance"])
                    else:
                        expected_location = json.dumps(location, ensure_ascii=False, sort_keys=True, separators=(",", ":")) if location else ""
                        self.assertEqual(expected_location, report["source_location"])
                        if location:
                            self.assertEqual(location, json.loads(report["source_location"]))
                        self.assertEqual(location, annotations["source_location_field_provenance"]["input_value"])
                        self.assertEqual({"note": "Original context."}, annotations["source_location_field_provenance"]["previous_provenance"])
                    self.assertEqual(reasons, annotations["paper_checklist_reasons"])
                    self.assertEqual(checklist("T2"), json.loads(report["paper_checklist"]))
                    self.assertEqual((PAPER, PAPER_TITLE, "HOM1", "TAIR:AT1G12345", "Arabidopsis thaliana"),
                                     tuple(report[key] for key in ("source_id", "source_title", "reported_gene", "reported_gene_id", "reported_species")))
                    self.assertTrue(built["review"].loc["g1", "relevance_tier"].startswith("T2"))
                    self.assertEqual(original, json.dumps(rows, sort_keys=True))
                    self.assertEqual(expected_bytes, (root / "light/light/b01/results.jsonl").read_bytes())

    def test_reason_alias_survives_converter_import_without_changing_raw_science(self):
        alias_reasons = {key: f"Original answer reason for {key}." for key in CHECKLIST_FIELDS}
        canonical_reasons = {key: f"Canonical answer reason for {key}." for key in CHECKLIST_FIELDS}
        for contract in ("legacy-v1", "science-v2"):
            for spelling in ("alias", "both", "annotated_canonical", "canonical"):
                with self.subTest(contract=contract, spelling=spelling), tempfile.TemporaryDirectory() as temporary:
                    root = Path(temporary)
                    paper = paper_row(PAPER, PAPER_TITLE, checklist("T2"), reported_gene_id="TAIR:AT1G12345",
                                      research_annotations={"counterevidence": "Original limitation retained.",
                                          "checklist_reason_field_provenance": {"source_note": "Original annotation retained."}})
                    if spelling != "canonical":
                        paper["checklist_reasons"] = alias_reasons
                    if spelling in {"both", "canonical"}:
                        paper["paper_checklist_reasons"] = canonical_reasons
                    if spelling == "annotated_canonical":
                        paper["research_annotations"]["paper_checklist_reasons"] = canonical_reasons
                    rows = [gene_row("g1", "REF-A", "T2", [paper])]
                    if contract == "science-v2":
                        rows = {"contract": light.SCIENCE_CONTRACT,
                                "paper_observations": [{"observation_id": "o1", **{
                                    key: value for key, value in paper.items() if key not in {"checklist", "gene_correspondence"}}}],
                                "gene_decisions": [{"interval_gene_id": "g1", "reference_gene": "REF-A",
                                    "state": "assessed", "basis": "paper_evidence", "function": "leaf development",
                                    "reason": "Original gene bridge.", "paper_assessments": [{"observation_id": "o1",
                                        "gene_correspondence": "homolog", "checklist": paper["checklist"],
                                        "reason": "Original copy-specific application."}]}]}
                    original = json.dumps(rows, sort_keys=True)
                    expected_bytes = (json.dumps(rows) if contract == "science-v2" else
                                      "".join(json.dumps(row) + "\n" for row in rows)).encode("utf-8")
                    built = run_pipeline(root, fixture(), rows, MAIN_CACHE, "u1\n", output_contract=contract)
                    proposals = pd.read_csv(built["imported"]["proposals"], sep="\t", dtype=str).fillna("")
                    proposal = proposals[proposals["member_interval_genes"] == "g1"].iloc[0]
                    report = json.loads(proposal["paper_reports_json"])[0]
                    annotations = report["research_annotations"]
                    expected = alias_reasons if spelling == "alias" else canonical_reasons
                    self.assertEqual(expected, annotations["paper_checklist_reasons"])
                    if spelling != "canonical":
                        originals = annotations["checklist_reason_field_provenance"]["input_fields"]
                        self.assertEqual(alias_reasons, originals["checklist_reasons"])
                        if spelling == "both":
                            self.assertEqual(canonical_reasons, originals["paper_checklist_reasons"])
                        self.assertEqual({"source_note": "Original annotation retained."},
                                         annotations["checklist_reason_field_provenance"]["previous_provenance"])
                    self.assertEqual("Original limitation retained.", annotations["counterevidence"])
                    self.assertEqual(checklist("T2"), json.loads(report["paper_checklist"]))
                    self.assertEqual(PAPER, report["source_id"])
                    self.assertEqual(PAPER_TITLE, report["source_title"])
                    self.assertEqual("HOM1", report["reported_gene"])
                    self.assertEqual("TAIR:AT1G12345", report["reported_gene_id"])
                    self.assertEqual("Arabidopsis thaliana", report["reported_species"])
                    self.assertTrue(built["review"].loc["g1", "relevance_tier"].startswith("T2"))
                    if contract == "science-v2":
                        self.assertEqual("Original copy-specific application.", annotations["copy_assessment"]["reason"])
                    self.assertEqual(original, json.dumps(rows, sort_keys=True))
                    self.assertEqual(expected_bytes, (root / "light/light/b01/results.jsonl").read_bytes())
