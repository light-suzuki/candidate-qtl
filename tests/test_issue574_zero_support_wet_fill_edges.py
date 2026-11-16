"""Additional edge coverage for Issue #574 zero-support Wet review filling."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from tests._issue664_helpers import manifest_for


ROUTE_B_UNVERIFIED = """candidate_gene_id\troute\tsupport_status
PsIgnored\tb\tnot_evaluated
"""
ATLAS = """gene_id\tfamily_id\tevolutionary_family_id\torthogroup\tsubfamily_id\tdomain_ids\tdomain_architecture\torthology_reference_gene
PsTier2\t\t\t\t\t\t\tREF1
PsTier3\tFAM1\t\t\t\t\t\t
PsDark\t\t\t\t\t\t\t
"""
REFERENCE_FUNCTION = """reference_gene_id\tmutant_phenotype\treference_species
REF1\treported mutant phenotype\tMedicago truncatula
"""
RANKED_EMPTY_COVERAGE = """gene_id\tdry_rank\tdry_score
OtherGene\t1\t9
"""


def _write(root: Path, name: str, text: str) -> Path:
    path = root / name
    path.write_text(text, encoding="utf-8")
    return path


class ZeroSupportWetFillEdgeTests(unittest.TestCase):
    def test_missing_interval_inventory_keeps_review_required_stop(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_b = _write(root, "route_b.tsv", ROUTE_B_UNVERIFIED)
            route_b_manifest = manifest_for(route_b, route="b", case_id="issue574-no-inventory")
            result = wet_candidate_handoff(
                str(root / "out"),
                route_b_candidates_arg=str(route_b),
                route_b_candidates_manifest_arg=str(route_b_manifest),
                wet_budget=3,
                wet_fill_mode="tiered_interval_fill",
                case_id="issue574-no-inventory",
            )
            self.assertEqual("no_supported_candidates", result["budget_status"])
            self.assertEqual("review_required_no_supported_candidates", result["handoff_state"])
            self.assertIsNone(result["wet_candidates"])
            self.assertFalse((root / "out" / "wet_candidates.tsv").exists())
            gate = json.loads(
                (root / "out" / "candidate_budget" / "candidate_budget_manifest.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(0, gate["reviewable_interval_fill_candidate_count"])
            self.assertEqual("candidate_set_empty_without_reviewable_interval_fill", gate["decision_reason"])

    def test_unproposed_inventory_does_not_fill_or_depend_on_ranked_input(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_b = _write(root, "route_b.tsv", ROUTE_B_UNVERIFIED)
            route_b_manifest = manifest_for(route_b, route="b", case_id="issue574-ordering")
            atlas = _write(root, "atlas.tsv", ATLAS)
            functions = _write(root, "functions.tsv", REFERENCE_FUNCTION)
            ranked = _write(root, "ranked.tsv", RANKED_EMPTY_COVERAGE)
            result = wet_candidate_handoff(
                str(root / "out"),
                route_b_candidates_arg=str(route_b),
                route_b_candidates_manifest_arg=str(route_b_manifest),
                ranked_arg=str(ranked),
                interval_atlas_arg=str(atlas),
                reference_function_arg=str(functions),
                wet_budget=3,
                wet_fill_mode="tiered_interval_fill",
                case_id="issue574-ordering",
            )
            self.assertEqual("no_supported_candidates", result["budget_status"])
            self.assertEqual("review_required_no_supported_candidates", result["handoff_state"])
            self.assertIsNone(result["wet_candidates"])
            self.assertFalse((root / "out" / "wet_candidates.tsv").exists())
            gate = json.loads(
                (root / "out" / "candidate_budget" / "candidate_budget_manifest.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(0, gate["reviewable_interval_fill_candidate_count"])


if __name__ == "__main__":
    unittest.main()
