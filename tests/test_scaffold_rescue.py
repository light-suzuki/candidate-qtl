"""Unplaced-scaffold rescue (issue #173)."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main


class ScaffoldRescueTests(unittest.TestCase):
    def _fixture(self, root: Path) -> dict[str, Path]:
        paths = {name: root / name for name in (
            "scaffolds.tsv", "markers.tsv", "synteny.tsv", "orthogroup.tsv", "placement.tsv", "genes.tsv")}
        paths["scaffolds.tsv"].write_text("\n".join([
            "scaffold_id\tlength",
            "scaf_unique\t180000",   # unique paired marker + synteny -> high
            "scaf_single\t90000",    # one unpaired unique marker -> medium
            "scaf_repeat\t40000",    # repeat marker on several scaffolds -> unresolved
            "scaf_repeat2\t35000",
            "scaf_family\t60000",    # gene-family context only -> weak
            "scaf_none\t20000",      # nothing -> not linked
        ]) + "\n", encoding="utf-8")
        paths["markers.tsv"].write_text("\n".join([
            "marker_id\tscaffold_id\tpair_id\tidentity",
            "M1_F\tscaf_unique\tM1\t100",
            "M1_R\tscaf_unique\tM1\t100",
            "M2_F\tscaf_single\t\t100",
            "REP_F\tscaf_repeat\tREP\t99",
            "REP_R\tscaf_repeat\tREP\t99",
            "REP_F\tscaf_repeat2\tREP\t99",
        ]) + "\n", encoding="utf-8")
        paths["synteny.tsv"].write_text(
            "scaffold_id\tblock_id\tanchor_species\nscaf_unique\tblk_7\tMedicago truncatula\n", encoding="utf-8")
        paths["orthogroup.tsv"].write_text("scaffold_id\torthogroup\nscaf_family\tOG0000123\n", encoding="utf-8")
        paths["placement.tsv"].write_text("scaffold_id\tsupport_type\tsource\n", encoding="utf-8")
        paths["genes.tsv"].write_text("\n".join([
            "gene_id\tscaffold_id",
            "S_G1\tscaf_unique", "S_G2\tscaf_unique", "S_G3\tscaf_single",
            "S_G4\tscaf_family", "S_G5\tscaf_none",
        ]) + "\n", encoding="utf-8")
        return paths

    def _run(self, root: Path, keys: tuple[str, ...] = ("markers", "synteny", "orthogroup", "genes")) -> Path:
        fixture = self._fixture(root)
        flags = {"markers": "--marker-hits", "synteny": "--synteny", "orthogroup": "--orthogroup-context",
                 "placement": "--placement-support", "genes": "--scaffold-genes"}
        argv = ["scaffold-rescue", "--scaffolds", str(fixture["scaffolds.tsv"]), "--out-dir", str(root / "rescue"),
                "--case-id", "demo"]
        for key in keys:
            argv += [flags[key], str(fixture[f"{key}.tsv"])]
        self.assertEqual(0, main(argv))
        return root / "rescue"

    def test_paired_marker_plus_synteny_is_high_confidence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            frame = pd.read_csv(out / "scaffold_rescue_candidates.tsv", sep="\t", dtype=str).fillna("").set_index("scaffold_id")
            self.assertEqual("scaffold_rescue_high", frame.loc["scaf_unique", "rescue_class"])
            self.assertIn("paired_marker", frame.loc["scaf_unique", "placement_support_type"])
            self.assertIn("synteny_bridge", frame.loc["scaf_unique", "placement_support_type"])
            self.assertEqual("high", frame.loc["scaf_unique", "placement_confidence"])

    def test_a_single_strong_support_is_only_medium(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            frame = pd.read_csv(out / "scaffold_rescue_candidates.tsv", sep="\t", dtype=str).fillna("").set_index("scaffold_id")
            # One unpaired unique marker is weak support, not strong.
            self.assertEqual("scaffold_rescue_weak", frame.loc["scaf_single", "rescue_class"])
            self.assertEqual("single_marker", frame.loc["scaf_single", "marker_support"])

    def test_repeat_markers_are_rejected_and_ambiguity_is_preserved(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            frame = pd.read_csv(out / "scaffold_rescue_candidates.tsv", sep="\t", dtype=str).fillna("").set_index("scaffold_id")
            # REP_F hits two scaffolds, so it cannot place either one.
            self.assertEqual("unresolved_scaffold", frame.loc["scaf_repeat", "rescue_class"])
            self.assertEqual("unresolved_scaffold", frame.loc["scaf_repeat2", "rescue_class"])
            self.assertIn("multi_mapping_markers", frame.loc["scaf_repeat", "conflicting_evidence"])
            reasons = set(pd.read_csv(out / "scaffold_rescue_review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("marker_multi_mapping", reasons)

    def test_gene_family_context_alone_is_never_high_confidence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            frame = pd.read_csv(out / "scaffold_rescue_candidates.tsv", sep="\t", dtype=str).fillna("").set_index("scaffold_id")
            self.assertEqual("scaffold_rescue_weak", frame.loc["scaf_family", "rescue_class"])
            self.assertEqual("yes", frame.loc["scaf_family", "orthogroup_context_support"])
            evidence = pd.read_csv(out / "scaffold_placement_evidence.tsv", sep="\t", dtype=str).fillna("")
            context = evidence[evidence.support_type == "orthogroup_context"].iloc[0]
            self.assertEqual("weak", context.strength)
            self.assertIn("never proves linkage", context.notes)

    def test_no_hit_is_not_negative_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            frame = pd.read_csv(out / "scaffold_rescue_candidates.tsv", sep="\t", dtype=str).fillna("").set_index("scaffold_id")
            self.assertEqual("scaffold_rejected_as_unlinked", frame.loc["scaf_none", "rescue_class"])
            self.assertIn("no placement evidence supplied", frame.loc["scaf_none", "reason"])
            manifest = json.loads((out / "scaffold_rescue_manifest.json").read_text(encoding="utf-8"))
            self.assertIn("never 'the causal gene is not on unplaced sequence'", manifest["no_hit_interpretation"])

    def test_rescued_genes_are_a_separate_class_and_carry_no_score(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            genes = pd.read_csv(out / "scaffold_gene_candidates.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(genes["candidate_class"].eq("scaffold_rescue").all())
            self.assertTrue(genes["rank_contribution"].eq("none").all())
            # Genes on an unlinked scaffold are not promoted into the candidate space.
            self.assertNotIn("S_G5", set(genes["gene_id"]))
            self.assertEqual({"S_G1", "S_G2", "S_G3", "S_G4"}, set(genes["gene_id"]))
            for path in sorted(out.glob("*.tsv")):
                columns = set(pd.read_csv(path, sep="\t", dtype=str).columns)
                self.assertFalse({"dry_score", "score", "dry_group"} & columns, path.name)

    def test_placement_confidence_is_separate_from_causal_confidence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            frame = pd.read_csv(out / "scaffold_rescue_candidates.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(frame["causal_confidence"].eq("not_assessed_here").all())
            manifest = json.loads((out / "scaffold_rescue_manifest.json").read_text(encoding="utf-8"))
            self.assertIn("not a causal probability", manifest["confidence_meaning"])

    def test_hic_support_can_reach_high_with_a_marker(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._fixture(root)
            fixture["placement.tsv"].write_text(
                "scaffold_id\tsupport_type\tsource\nscaf_single\thic_linkage\tjuicer\n", encoding="utf-8")
            out = root / "rescue"
            self.assertEqual(0, main(["scaffold-rescue", "--scaffolds", str(fixture["scaffolds.tsv"]),
                                      "--marker-hits", str(fixture["markers.tsv"]),
                                      "--placement-support", str(fixture["placement.tsv"]),
                                      "--out-dir", str(out)]))
            frame = pd.read_csv(out / "scaffold_rescue_candidates.tsv", sep="\t", dtype=str).fillna("").set_index("scaffold_id")
            # One strong support (Hi-C) plus a weak single marker: medium, not high.
            self.assertEqual("scaffold_rescue_medium", frame.loc["scaf_single", "rescue_class"])

    def test_no_placement_input_rescues_nothing_and_says_so(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = self._run(root, keys=("genes",))
            frame = pd.read_csv(out / "scaffold_rescue_candidates.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(frame["rescue_class"].eq("scaffold_rejected_as_unlinked").all())
            self.assertEqual(0, len(pd.read_csv(out / "scaffold_gene_candidates.tsv", sep="\t", dtype=str)))
            reasons = set(pd.read_csv(out / "scaffold_rescue_review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("no_placement_evidence_supplied", reasons)


if __name__ == "__main__":
    unittest.main()
