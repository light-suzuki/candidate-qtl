"""Issue #570 regression: ordinary Wet selection does not pad short pools.

The interval inventory remains auditable, but only supported or explicitly
reviewed candidates enter ordinary Wet selection. A short pool stays short.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from tests._issue664_helpers import manifest_for

SUPPORTED = """candidate_gene_id\troute_a_supported\troute_b_supported\troute_support
Ps01\tFalse\tTrue\tb
Ps02\tFalse\tTrue\tb
"""
# Ps03/Ps04 resolve to reference genes with a reported mutant phenotype,
# Ps05/Ps06 only have a family or domain assignment, Ps07/Ps08 are dark.
ATLAS = """gene_id\tfamily_id\tevolutionary_family_id\torthogroup\tsubfamily_id\tdomain_ids\tdomain_architecture\torthology_reference_gene\tconservation_class
Ps01\tNOOT\t\t\t\t\t\tQ2HW56\trelative_enriched
Ps02\tCLE\t\t\t\t\t\tG7KXA1\trelative_enriched
Ps03\t\t\t\t\t\t\tQ6UD73\trelative_enriched
Ps04\tLYK\t\t\t\t\t\tA8WHB3\trelative_enriched
Ps05\tPAL\t\t\t\t\t\t\trelative_enriched
Ps06\t\t\t\t\tPF00069\tPkinase\t\trelative_enriched
Ps07\t\t\t\t\t\t\t\tno_detectable_homology
Ps08\t\t\t\t\t\t\t\tno_detectable_homology
"""
REFERENCE_FUNCTION = """reference_gene_id\tfunction_description\tmutant_phenotype\treference_species
Q6UD73\tLysM domain receptor-like kinase 3\tnodulation deficient, infection threads arrested\tMedicago truncatula
A8WHB3\tLysM domain receptor-like kinase\troot hair curling without infection\tMedicago truncatula
G7KXA1\tCLE peptide precursor\t\tMedicago truncatula
"""


# A ranked table built before the review slots were filled: it cannot carry the
# frozen membership, so ordering must step aside rather than reorder a subset.
RANKED_BEFORE_FILL = """gene_id\tdry_rank\tdry_score
Ps01\t1\t9
Ps02\t2\t8
"""


def _write(root: Path, name: str, text: str) -> Path:
    path = root / name
    path.write_text(text, encoding="utf-8")
    return path


def _route_b_candidates(root: Path, name: str, candidates: Path) -> Path:
    frame = pd.read_csv(candidates, sep="\t", dtype=str).fillna("")
    route_b = root / f"{name}_route_b_candidates.tsv"
    pd.DataFrame([
        {
            "candidate_gene_id": gene,
            "route": "b",
            "support_status": "verified_semantic_support",
        }
        for gene in frame["candidate_gene_id"].astype(str)
    ]).to_csv(route_b, sep="\t", index=False, lineterminator="\n")
    manifest_for(route_b, route="b", case_id="case570")
    return route_b


class TieredWetRetentionTests(unittest.TestCase):
    def _run(self, root: Path, *, wet_budget: int = 5, supported: str = SUPPORTED,
             atlas: str | None = ATLAS, reference_function: str | None = REFERENCE_FUNCTION,
             fill_mode: str = "tiered_interval_fill", out_name: str = "out"):
        supported_path = _write(root, f"{out_name}_supported.tsv", supported)
        atlas_path = _write(root, f"{out_name}_atlas.tsv", atlas) if atlas is not None else None
        function_path = (_write(root, f"{out_name}_function.tsv", reference_function)
                         if reference_function is not None else None)
        route_b = _route_b_candidates(root, out_name, supported_path)
        result = wet_candidate_handoff(
            str(root / out_name),
            route_b_candidates_arg=str(route_b),
            route_b_candidates_manifest_arg=str(route_b.with_name(route_b.stem + "_manifest.json")),
            wet_budget=wet_budget,
            interval_atlas_arg=str(atlas_path) if atlas_path else None,
            reference_function_arg=str(function_path) if function_path else None,
            wet_fill_mode=fill_mode,
            case_id="case570",
        )
        membership = pd.read_csv(root / out_name / "wet_candidates.tsv", sep="\t", dtype=str).fillna("")
        manifest = json.loads(
            (root / out_name / "wet_candidates_manifest.json").read_text(encoding="utf-8"))
        return result, membership, manifest

    def test_unproposed_interval_inventory_does_not_fill_review_slots(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, membership, manifest = self._run(root)
            self.assertEqual(["Ps01", "Ps02"],
                             membership["candidate_gene_id"].tolist())
            self.assertEqual(2, manifest["supported_candidate_count"])
            self.assertEqual(0, manifest["review_coverage_candidate_count"])
            self.assertEqual(3, manifest["unfilled_budget_slots"])
            self.assertEqual("evidence_scoped_candidates", manifest["wet_fill_state"])

    def test_unproposed_inventory_rows_are_not_membership(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, membership, _ = self._run(root)
            self.assertEqual({"Ps01", "Ps02"}, set(membership["candidate_gene_id"]))
            self.assertEqual({"verified_route_b_support"}, set(membership["membership_basis"]))
            supported = membership[membership["membership_basis"] == "verified_route_b_support"]
            self.assertEqual({"True"}, set(supported["route_b_supported"]))

    def test_unresolved_rescue_is_retained_without_inventing_support(self):
        """Identity-only and U inventory do not become candidates by exhaustion."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, membership, manifest = self._run(root, wet_budget=20)
            self.assertEqual({"Ps01", "Ps02"}, set(membership["candidate_gene_id"]))
            self.assertEqual(2, len(membership))
            self.assertEqual(18, manifest["unfilled_budget_slots"])
            self.assertNotIn("Ps07", set(membership["candidate_gene_id"]))
            self.assertEqual("U_unresolved", pd.read_csv(
                root / "out" / "wet_relevance_inventory.tsv", sep="\t", dtype=str
            ).set_index("candidate_gene_id").loc["Ps07", "relevance_tier"])

    def test_supported_only_mode_reproduces_the_previous_membership(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, membership, manifest = self._run(root, fill_mode="supported_only")
            self.assertEqual(["Ps01", "Ps02"], membership["candidate_gene_id"].tolist())
            self.assertEqual(0, manifest["review_coverage_candidate_count"])
            self.assertEqual("not_requested", manifest["wet_fill_state"])

    def test_without_an_interval_inventory_nothing_is_invented(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, membership, manifest = self._run(root, atlas=None, reference_function=None)
            self.assertEqual(["Ps01", "Ps02"], membership["candidate_gene_id"].tolist())
            self.assertEqual("no_interval_inventory", manifest["wet_fill_state"])

    def test_support_status_does_not_replace_relevance(self):
        """Verified status alone does not move an unresolved member above T4."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            supported = "candidate_gene_id\troute_a_supported\troute_b_supported\troute_support\n" + "".join(
                f"Ps{index:02d}\tFalse\tTrue\tb\n" for index in range(1, 8))
            supported_path = _write(root, "over_supported.tsv", supported)
            route_b = _route_b_candidates(root, "over", supported_path)
            route_b_manifest = route_b.with_name(route_b.stem + "_manifest.json")
            result = wet_candidate_handoff(
                str(root / "over"),
                route_b_candidates_arg=str(route_b),
                route_b_candidates_manifest_arg=str(route_b_manifest),
                wet_budget=5,
                interval_atlas_arg=str(_write(root, "over_atlas.tsv", ATLAS)),
                case_id="case570",
            )
            self.assertEqual("above_wet_target", result["budget_status"])
            wet = pd.read_csv(root / "over" / "wet_candidates.tsv", sep="\t")
            self.assertEqual(7, len(wet))
            self.assertTrue(wet["support_tier"].eq("tier1_verified_support").all())
            inventory = pd.read_csv(root / "over" / "wet_relevance_inventory.tsv", sep="\t")
            self.assertEqual(8, len(inventory))
            self.assertEqual("U_unresolved", inventory.set_index("candidate_gene_id").loc["Ps07", "relevance_tier"])

    def test_ranked_input_covers_the_unpadded_frozen_membership(self):
        """No synthetic fill rows are added beyond the ranked supported set."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ranked = _write(root, "ranked.tsv", RANKED_BEFORE_FILL)
            supported_path = _write(root, "ordering_supported.tsv", SUPPORTED)
            route_b = _route_b_candidates(root, "ordering", supported_path)
            result = wet_candidate_handoff(
                str(root / "ordering"),
                route_b_candidates_arg=str(route_b),
                route_b_candidates_manifest_arg=str(route_b.with_name(route_b.stem + "_manifest.json")),
                ranked_arg=str(ranked),
                wet_budget=5,
                interval_atlas_arg=str(_write(root, "ordering_atlas.tsv", ATLAS)),
                reference_function_arg=str(_write(root, "ordering_function.tsv", REFERENCE_FUNCTION)),
                case_id="case570",
            )
            self.assertEqual("ordering_only", result["narrowing_state"])
            self.assertIsNone(result["ordering_skip_reason"])
            membership = pd.read_csv(
                root / "ordering" / "wet_candidates.tsv", sep="	", dtype=str).fillna("")
            self.assertEqual(2, len(membership))

    def test_small_target_does_not_add_unproposed_same_or_lower_tier_members(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, membership, _ = self._run(root, wet_budget=2)
            self.assertEqual(["Ps01", "Ps02"], membership["candidate_gene_id"].tolist())
            self.assertEqual({"tier1_verified_support"}, set(membership["support_tier"]))


if __name__ == "__main__":
    unittest.main()
