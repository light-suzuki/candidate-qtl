"""Regression tests for canonical Wet membership semantics (#502, #503)."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from tests._issue664_helpers import manifest_for


def _bridge(root: Path, genes: list[str]) -> Path:
    path = root / "bridge.tsv"
    pd.DataFrame([
        {
            "route": "a",
            "gene_id": gene,
            "family_id": "FAM1",
            "hypothesis_id": "H1",
            "mapping_tier": "tier3_family_projection",
            "ambiguity_status": "unique",
            "review_status": "supported",
        }
        for gene in genes
    ]).to_csv(path, sep="\t", index=False, lineterminator="\n")
    manifest_for(path, route="a", case_id="case1")
    return path


def _ranked(root: Path, genes: list[str]) -> Path:
    path = root / "ranked.tsv"
    pd.DataFrame([
        {
            "gene_id": gene,
            "dry_rank": str(index),
            "dry_group": "top" if index <= 3 else "reserve",
            "dry_score": str(100 - index),
        }
        for index, gene in enumerate(genes, start=1)
    ]).to_csv(path, sep="\t", index=False, lineterminator="\n")
    return path


class CanonicalWetMembershipTests(unittest.TestCase):
    def test_within_budget_freezes_every_supported_candidate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            genes = ["G002", "G004", "G006"]
            result = wet_candidate_handoff(
                str(root / "out"),
                route_a_bridge_arg=str(_bridge(root, genes)),
                route_a_manifest_arg=str(root / "bridge_manifest.json"),
                wet_budget=20,
                case_id="case1",
            )

            self.assertEqual("wet_reviewable", result["handoff_state"])
            self.assertEqual("frozen", result["membership_state"])
            wet = pd.read_csv(result["wet_candidates"], sep="\t", dtype=str).fillna("")
            self.assertEqual(set(genes), set(wet["candidate_gene_id"]))
            self.assertEqual(3, len(wet))

            manifest = json.loads(Path(result["wet_candidates_manifest"]).read_text(encoding="utf-8"))
            self.assertEqual("canonical_wet_membership", manifest["producer_class"])
            self.assertEqual(3, manifest["candidate_count"])
            self.assertEqual(3, manifest["supported_candidate_count"])
            self.assertEqual(20, manifest["wet_budget"])
            # Without an interval inventory there is nothing to fill review
            # slots with, so membership is still exactly the supported set.
            self.assertEqual(0, manifest["review_coverage_candidate_count"])
            self.assertEqual("no_interval_inventory", manifest["wet_fill_state"])
            self.assertIn("never used as padding", manifest["membership_rule"])

    def test_over_budget_never_emits_a_truncated_wet_membership(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            genes = [f"G{index:03d}" for index in range(1, 26)]
            result = wet_candidate_handoff(
                str(root / "out"),
                route_a_bridge_arg=str(_bridge(root, genes)),
                route_a_manifest_arg=str(root / "bridge_manifest.json"),
                wet_budget=20,
                case_id="case1",
            )

            self.assertEqual("wet_reviewable", result["handoff_state"])
            self.assertEqual("frozen", result["membership_state"])
            wet = pd.read_csv(result["wet_candidates"], sep="\t")
            self.assertEqual(set(genes), set(wet["candidate_gene_id"]))

    def test_ordering_receives_only_frozen_wet_membership(self) -> None:
        """Genes outside canonical membership cannot consume Wet ordering slots."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            wet_genes = ["G002", "G004"]
            ranked_genes = ["G001", "G002", "G003", "G004", "G005"]
            result = wet_candidate_handoff(
                str(root / "out"),
                route_a_bridge_arg=str(_bridge(root, wet_genes)),
                route_a_manifest_arg=str(root / "bridge_manifest.json"),
                ranked_arg=str(_ranked(root, ranked_genes)),
                wet_budget=20,
                case_id="case1",
            )

            self.assertEqual("ordering_only", result["narrowing_state"])
            wet = pd.read_csv(result["wet_candidates"], sep="\t", dtype=str).fillna("")
            narrowed = pd.read_csv(result["stage_a_narrowed"], sep="\t", dtype=str).fillna("")
            self.assertEqual(set(wet["candidate_gene_id"]), set(narrowed["gene_id"]))
            self.assertEqual({"G002", "G004"}, set(narrowed["gene_id"]))
            self.assertEqual(["1", "2"], list(narrowed["narrowing_priority"]))

            handoff = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual("ordering_only", handoff["stage_a_narrowed_producer_class"])
            self.assertTrue(handoff["ordering_membership_preserved"])
            self.assertIn("cannot create, add, remove, pad, or truncate", handoff["ordering_policy"])

    def test_ordering_projects_family_and_copy_context_to_wet_set(self) -> None:
        """Global family/copy artifacts may contain genes outside final membership."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            wet_genes = ["G002", "G004"]
            ranked_genes = ["G001", "G002", "G003", "G004", "G999"]

            recovery = root / "family_recovery.tsv"
            pd.DataFrame([{
                "hypothesis_id": "H1",
                "gene_family": "FAM1",
                "member_interval_genes": "G002;G004;G999",
            }]).to_csv(recovery, sep="\t", index=False, lineterminator="\n")

            copy = root / "family_copy.tsv"
            pd.DataFrame([
                {"family_id": "FAM1", "candidate_gene_id": "G002", "disambiguation_status": "family_resolution_required"},
                {"family_id": "FAM1", "candidate_gene_id": "G004", "disambiguation_status": "family_resolution_required"},
                {"family_id": "FAM1", "candidate_gene_id": "G999", "disambiguation_status": "preferred_copy"},
            ]).to_csv(copy, sep="\t", index=False, lineterminator="\n")

            shortlist = root / "family_copy_shortlist.tsv"
            pd.DataFrame([
                {"candidate_gene_id": "G999", "first_wet_priority": "first"},
            ]).to_csv(shortlist, sep="\t", index=False, lineterminator="\n")

            result = wet_candidate_handoff(
                str(root / "out"),
                route_a_bridge_arg=str(_bridge(root, wet_genes)),
                route_a_manifest_arg=str(root / "bridge_manifest.json"),
                ranked_arg=str(_ranked(root, ranked_genes)),
                family_recovery_arg=str(recovery),
                family_copy_arg=str(copy),
                family_copy_shortlist_arg=str(shortlist),
                wet_budget=20,
                case_id="case1",
            )

            narrowed = pd.read_csv(result["stage_a_narrowed"], sep="\t", dtype=str).fillna("")
            self.assertEqual({"G002", "G004"}, set(narrowed["gene_id"]))
            self.assertNotIn("G999", set(narrowed["gene_id"]))

            handoff = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            projected_recovery = pd.read_csv(
                handoff["ordering_context_inputs"]["family_recovery"], sep="\t", dtype=str
            ).fillna("")
            self.assertEqual("G002;G004", projected_recovery.iloc[0]["member_interval_genes"])
            projected_copy = pd.read_csv(
                handoff["ordering_context_inputs"]["family_copy"], sep="\t", dtype=str
            ).fillna("")
            self.assertEqual({"G002", "G004"}, set(projected_copy["candidate_gene_id"]))
            projected_shortlist = pd.read_csv(
                handoff["ordering_context_inputs"]["family_copy_shortlist"], sep="\t", dtype=str
            ).fillna("")
            self.assertTrue(projected_shortlist.empty)


if __name__ == "__main__":
    unittest.main()
