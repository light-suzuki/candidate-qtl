"""The A/B first-stage JSON examples must place paper_reports_json per record.

This extracts the concrete example from the actually generated Route A prompt
and the current Route B prompt, fills synthetic but valid report/claim fields,
and passes them through the real importers. It guards the root cause fixed by
#912: the A payload template and the B claims example omitted paper_reports_json,
so the per-record paper placement was never shown to the model even though the
shared Tier checklist prompt described it.
"""
from __future__ import annotations

import json
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path

import pandas as pd

from candidate_qtl import phenotype_investigation as pi
from candidate_qtl import route_b_semantic_bundle as rb
from candidate_qtl.cli import main
from candidate_qtl.observation_tiers import CHECKLIST_FIELDS
from candidate_qtl.paper_reports import PAPER_FIELD
from candidate_qtl.route_b_semantic_import import import_route_b_semantic_proposals
from tests import test_route_b_semantic_import as b_fixtures


def _example_after(text: str, marker: str) -> dict:
    index = text.index(marker)
    payload, _ = json.JSONDecoder().raw_decode(text[index + len(marker):].lstrip())
    return payload


def _valid_report(source_id: str = "PMID:123") -> dict:
    """A synthetic but canonical report; the checklist answers only two keys."""
    return {
        "source_id": source_id,
        "source_title": "Synthetic paper " + source_id,
        "reported_gene": "EXPERIMENTAL_GENE",
        "reported_gene_id": "NCBI Gene:111111",
        "reported_species": "Plant beta",
        "reference_gene_id": "REF1",
        "gene_correspondence": "homolog",
        "initial_tier_reason": "Synthetic source reports a similar cellular phenotype.",
        "paper_checklist": {"phenotype_reported": "yes", "identity_resolved": "yes"},
    }


def _filled_hypothesis(template: dict, family: str, report: list) -> dict:
    hypothesis = deepcopy(template)
    hypothesis.update(
        causal_layer="regulation",
        pathway_or_process=f"process {family}",
        gene_family=family,
        mechanism_statement="a lesion in this family changes the process",
        model_species_evidence="Arabidopsis mutant study",
        model_species_mutant_phenotype="altered development",
        source_ids="PMID:123",
        proposed_stable_id="P12345",
        reference_species="Arabidopsis thaliana",
        reference_gene_name="REF1",
        selection_rationale="characterised reference outside the focal species",
        counter_hypothesis="a different family could produce the same phenotype",
        discriminating_observation="measure the process directly",
        confidence="medium",
        uncertainty="copy separation is unresolved",
        paper_reports_json=report,
    )
    return hypothesis


class InitialPaperPromptExampleTests(unittest.TestCase):
    def test_reference_anchor_roles_preserve_the_original_t2_answers(self) -> None:
        from candidate_qtl.observation_tiers import classify_observation
        from candidate_qtl.paper_reports import bound_paper_rows

        report = _valid_report()
        report["paper_checklist"] = {field: "unknown" for field in CHECKLIST_FIELDS}
        report["paper_checklist"].update(phenotype_reported="yes", similar_cellular_or_developmental_phenotype="yes")
        row = {"reference_gene": "REF1", "source_ids": "PMID:123",
               "gene_id": "focal-copy-1", PAPER_FIELD: json.dumps([report])}
        bound, = bound_paper_rows(row)
        self.assertEqual("T2", classify_observation(bound)[0])
        self.assertEqual(report["paper_checklist"], json.loads(bound["paper_checklist"]))

        # A local member/profile is a different role, even with identical nine answers.
        report["reference_gene_id"] = "focal-copy-1 (profile-display)"
        row[PAPER_FIELD] = json.dumps([report])
        unresolved, = bound_paper_rows(row)
        self.assertEqual("unresolved", unresolved["tier_gene_binding"])
        self.assertEqual("U", classify_observation(unresolved)[0])
        self.assertEqual(report["paper_checklist"], json.loads(unresolved["paper_checklist"]))

    def _route_a_example(self, root: Path) -> dict:
        result = pi.phenotype_investigate(
            str(root / "inv"), "case-912", "Pisum sativum", "wrinkled seed",
        )
        text = Path(result["prompt"]).read_text(encoding="utf-8")
        return _example_after(text, "Return JSON only:")

    def _route_b_prompt(self, root: Path) -> str:
        built = b_fixtures.RouteBSemanticImportTests()._bundle(root)
        prompt = json.loads(Path(built["prompts"]).read_text(encoding="utf-8").splitlines()[0])
        return built, prompt["system"]

    def test_route_a_example_places_reports_per_hypothesis_and_imports(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            example = self._route_a_example(root)
            # Per-record placement: never at the root, always inside a hypothesis.
            self.assertNotIn(PAPER_FIELD, example)
            self.assertTrue(example["hypotheses"])
            for hypothesis in example["hypotheses"]:
                self.assertIn(PAPER_FIELD, hypothesis)
            report_example = example["hypotheses"][0][PAPER_FIELD][0]
            self.assertEqual("unknown", report_example["paper_checklist"][CHECKLIST_FIELDS[0]])
            self.assertEqual(set(CHECKLIST_FIELDS), set(report_example["paper_checklist"]))

            hypotheses = [
                _filled_hypothesis(example["hypotheses"][0], "FAMILY_A", [_valid_report()]),
                _filled_hypothesis(example["hypotheses"][0], "FAMILY_B", []),
            ]
            payload = root / "payload.json"
            payload.write_text(json.dumps({"hypotheses": hypotheses}), encoding="utf-8")
            investigation = root / "inv" / "phenotype_investigation.tsv"
            self.assertEqual(0, main([
                "mechanism-import", "--payload", str(payload), "--out-dir", str(root / "import"),
                "--investigation", str(investigation),
            ]))
            rows = pd.read_csv(root / "import" / "mechanism_hypotheses.tsv", sep="\t", dtype=str).fillna("")

        self.assertEqual(2, len(rows))
        placement = {row["gene_family"]: json.loads(row[PAPER_FIELD]) for row in rows.to_dict("records")}
        self.assertEqual(1, len(placement["FAMILY_A"]))
        self.assertEqual([], placement["FAMILY_B"])
        imported = placement["FAMILY_A"][0]
        checklist = json.loads(imported["paper_checklist"])
        # Normalization preserves every documented question; unset answers stay unknown.
        self.assertEqual(set(CHECKLIST_FIELDS), set(checklist))
        self.assertEqual(len(CHECKLIST_FIELDS) - 2, sum(value == "unknown" for value in checklist.values()))
        self.assertEqual("yes", checklist["phenotype_reported"])
        self.assertEqual("PMID:123", imported["source_id"])
        self.assertEqual("REF1", imported["reference_gene_id"])
        self.assertEqual("homolog", imported["gene_correspondence"])

    def test_route_b_current_example_places_reports_per_claim_and_imports(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            built, system = self._route_b_prompt(root)
            example = _example_after(system, "Return JSON only as")
            self.assertNotIn(PAPER_FIELD, example)
            self.assertTrue(example["claims"])
            for claim in example["claims"]:
                self.assertIn(PAPER_FIELD, claim)
            report_example = example["claims"][0][PAPER_FIELD][0]
            self.assertEqual(set(CHECKLIST_FIELDS), set(report_example["paper_checklist"]))
            self.assertTrue(all(value == "unknown" for value in report_example["paper_checklist"].values()))

            fixture = b_fixtures.RouteBSemanticImportTests()
            claim = deepcopy(example["claims"][0])
            claim.update(
                claim_type="mutant_phenotype",
                biological_function="meristem maintenance",
                mechanism="receptor signaling",
                phenotype_relation="could explain fasciation",
                reference_gene="REF1",
                reference_species="Arabidopsis thaliana",
                evidence_scope="reference_gene",
                source_ids=["PMID:123"],
                alternative_hypothesis="a related copy may act instead",
                uncertainty_note="transfer to the focal crop remains unresolved",
                paper_reports_json=[_valid_report()],
            )
            response = root / "response.json"
            fixture._response(Path(built["prompts"]), [claim], response)
            result = import_route_b_semantic_proposals(
                built["bundle"], built["prompts"], str(response), str(root / "import"),
            )
            proposals = pd.read_csv(result["proposals"], sep="\t", dtype=str).fillna("")

        self.assertEqual(1, result["claim_row_count"])
        self.assertEqual(1, result["proposed_claim_row_count"])
        row = proposals.iloc[0]
        reports = json.loads(row[PAPER_FIELD])
        self.assertEqual(1, len(reports))
        checklist = json.loads(reports[0]["paper_checklist"])
        self.assertEqual(set(CHECKLIST_FIELDS), set(checklist))
        self.assertEqual("T6", row["initial_tier"])
        self.assertEqual("PMID:123", row["source_id"])

    def test_route_b_empty_report_list_is_allowed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            built, system = self._route_b_prompt(root)
            claim = deepcopy(_example_after(system, "Return JSON only as")["claims"][0])
            claim.update(
                claim_type="mechanistic_review_seed",
                mechanism="a specific process could explain the phenotype",
                phenotype_relation="trait transfer remains unverified",
                reference_gene="REF1",
                evidence_scope="reference_gene",
                source_ids=[],
                alternative_hypothesis="another process may be causal",
                uncertainty_note="review seed only",
                paper_reports_json=[],
            )
            fixture = b_fixtures.RouteBSemanticImportTests()
            response = root / "response.json"
            fixture._response(Path(built["prompts"]), [claim], response)
            result = import_route_b_semantic_proposals(
                built["bundle"], built["prompts"], str(response), str(root / "import"),
            )
            proposals = pd.read_csv(result["proposals"], sep="\t", dtype=str).fillna("")

        self.assertEqual(1, result["claim_row_count"])
        self.assertEqual("[]", proposals.iloc[0][PAPER_FIELD])

    def test_legacy_route_b_content_example_is_unchanged(self) -> None:
        self.assertNotIn(PAPER_FIELD, rb.CONTENT_SYSTEM_PROMPT)
        example = _example_after(rb.CONTENT_SYSTEM_PROMPT, "Return JSON only as")
        self.assertEqual({"claims"}, set(example))
        for claim in example["claims"]:
            self.assertNotIn(PAPER_FIELD, claim)


if __name__ == "__main__":
    unittest.main()
