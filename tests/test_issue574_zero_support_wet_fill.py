"""Issue #574: zero verified support needs an explicit eligible review seed."""

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
Ps03\t\t\t\t\t\t\tQ6UD73
Ps04\tPAL\t\t\t\t\t\t
Ps05\t\t\t\t\t\t\t
"""

REFERENCE_FUNCTION = """reference_gene_id\tmutant_phenotype\treference_species
Q6UD73\tnodulation deficient, infection threads arrested\tMedicago truncatula
"""

ALL_DARK_ATLAS = """gene_id\tfamily_id\tevolutionary_family_id\torthogroup\tsubfamily_id\tdomain_ids\tdomain_architecture\torthology_reference_gene
Ps05\t\t\t\t\t\t\t
Ps06\t\t\t\t\t\t\t
"""


def _write(root: Path, name: str, text: str) -> Path:
    path = root / name
    path.write_text(text, encoding="utf-8")
    return path


class ZeroSupportWetFillTests(unittest.TestCase):
    def _run(
        self,
        root: Path,
        *,
        fill_mode: str = "tiered_interval_fill",
        atlas: str | None = ATLAS,
        reference_function: str | None = REFERENCE_FUNCTION,
        out_name: str = "out",
    ) -> dict[str, object]:
        route_b = _write(root, f"{out_name}_route_b.tsv", ROUTE_B_UNVERIFIED)
        route_b_manifest = manifest_for(route_b, route="b", case_id="case574")
        atlas_path = _write(root, f"{out_name}_atlas.tsv", atlas) if atlas is not None else None
        function_path = (
            _write(root, f"{out_name}_reference.tsv", reference_function)
            if reference_function is not None
            else None
        )
        return wet_candidate_handoff(
            str(root / out_name),
            route_b_candidates_arg=str(route_b),
            route_b_candidates_manifest_arg=str(route_b_manifest),
            interval_atlas_arg=str(atlas_path) if atlas_path else None,
            reference_function_arg=str(function_path) if function_path else None,
            wet_budget=3,
            wet_fill_mode=fill_mode,
            case_id="case574",
        )

    def test_zero_support_does_not_promote_unproposed_interval_inventory(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = self._run(root)

            self.assertEqual("no_supported_candidates", result["budget_status"])
            self.assertEqual("review_required_no_supported_candidates", result["handoff_state"])
            self.assertEqual("route_b_declared_but_unverified", result["input_coverage_state"])
            self.assertIsNone(result["wet_candidates"])
            self.assertFalse((root / "out" / "wet_candidates.tsv").exists())

            gate_manifest = json.loads(
                (root / "out" / "candidate_budget" / "candidate_budget_manifest.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(
                "candidate_set_empty_without_reviewable_interval_fill",
                gate_manifest["decision_reason"],
            )
            self.assertEqual(0, gate_manifest["reviewable_interval_fill_candidate_count"])

    def test_supported_only_keeps_zero_support_stop(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = self._run(root, fill_mode="supported_only", out_name="supported_only")
            self.assertEqual("no_supported_candidates", result["budget_status"])
            self.assertEqual("review_required_no_supported_candidates", result["handoff_state"])
            self.assertIsNone(result["wet_candidates"])
            self.assertFalse((root / "supported_only" / "wet_candidates.tsv").exists())

    def test_all_dark_inventory_does_not_create_zero_support_candidates(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = self._run(
                root,
                atlas=ALL_DARK_ATLAS,
                reference_function=None,
                out_name="dark",
            )
            self.assertEqual("no_supported_candidates", result["budget_status"])
            self.assertEqual("review_required_no_supported_candidates", result["handoff_state"])
            self.assertIsNone(result["wet_candidates"])
            self.assertFalse((root / "dark" / "wet_candidates.tsv").exists())

            gate_manifest = json.loads(
                (root / "dark" / "candidate_budget" / "candidate_budget_manifest.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(0, gate_manifest["reviewable_interval_fill_candidate_count"])

    def test_zero_support_manifest_describes_evidence_scoped_default(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._run(root, out_name="policy")
            manifest = json.loads(
                (root / "policy" / "wet_candidate_handoff_manifest.json").read_text(encoding="utf-8")
            )
            policy = "\n".join(manifest["policy"])
            self.assertIn("ordinary target handoff retains T1 and existing eligible union/review candidates", policy)
            self.assertIn("excluding identity-only inventory", policy)


if __name__ == "__main__":
    unittest.main()
