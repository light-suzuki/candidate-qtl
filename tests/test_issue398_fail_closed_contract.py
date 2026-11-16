from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main


FAMILIES = (
    "hypothesis_id\tfamily_id\tmember_interval_genes\n"
    "hyp398\tFAM398\tG_A\n"
)

ATLAS = (
    "gene_id\tchrom\tstart\tend\tfamily_label\tdomain_architecture\tdomain_ids\tdomain_source\torthology_reference_gene\n"
    "G_A\tchr1\t1000\t2000\tFAM398\t\t\t\t\n"
)

BLAST = (
    "panel_id\treview_status\tprovider_status\ttop_hit_target_id\tambiguity_status\t"
    "top_hit_identity\ttop_hit_query_coverage\ttop_hit_target_coverage\n"
    "B398\tverified\tevaluable\tG_A\tunique\t80\t0.90\t0.90\n"
)


def _bridge(rows: str) -> str:
    return "gene_id\tfamily_id\troute\tmapping_tier\n" + rows


class Issue398FailClosedContractTests(unittest.TestCase):
    """Regression contract for #398.

    A Route-A BLAST ``unique`` hit may create a single-copy promotion only after
    the interval-paralog guard was actually evaluated.  Missing bridge coverage
    or tier6-only coverage is *not* evidence that no paralog exists.
    """

    def _run(self, root: Path, bridge_rows: str | None) -> tuple[Path, pd.DataFrame, dict[str, object]]:
        families = root / "families.tsv"
        atlas = root / "atlas.tsv"
        blast = root / "blast.tsv"
        families.write_text(FAMILIES, encoding="utf-8")
        atlas.write_text(ATLAS, encoding="utf-8")
        blast.write_text(BLAST, encoding="utf-8")

        argv = [
            "family-to-copy",
            "--families", str(families),
            "--atlas", str(atlas),
            "--route-a-blast", str(blast),
            "--out-dir", str(root / "out"),
            "--case-id", "issue398",
        ]
        if bridge_rows is not None:
            bridge = root / "bridge.tsv"
            bridge.write_text(_bridge(bridge_rows), encoding="utf-8")
            argv += ["--bridge", str(bridge)]

        self.assertEqual(0, main(argv))
        out = root / "out"
        frame = pd.read_csv(
            out / "family_copy_discrimination.tsv", sep="\t", dtype=str
        ).fillna("").set_index("candidate_gene_id")
        manifest = json.loads((out / "family_copy_manifest.json").read_text(encoding="utf-8"))
        return out, frame, manifest

    def _assert_guard_not_evaluated_hold(self, out: Path, frame: pd.DataFrame, manifest: dict[str, object]) -> None:
        row = frame.loc["G_A"]
        self.assertNotIn(row["disambiguation_status"], {"preferred_copy", "co_preferred_tie"})
        self.assertEqual("", row["copy_set_evidence_source"])
        self.assertEqual("paralog_guard_not_evaluated", row["route_a_blast_guard"])
        shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
        self.assertEqual(0, len(shortlist))
        self.assertFalse(manifest["route_a_blast"]["paralog_guard_evaluated"])

    def test_no_bridge_never_turns_unique_hit_into_single_copy_promotion(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out, frame, manifest = self._run(Path(tmp), bridge_rows=None)
            self._assert_guard_not_evaluated_hold(out, frame, manifest)

    def test_tier6_only_bridge_never_counts_as_guard_evaluation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out, frame, manifest = self._run(
                Path(tmp),
                bridge_rows="G_A\tFAM398\ta\ttier6_not_evaluated\n",
            )
            self._assert_guard_not_evaluated_hold(out, frame, manifest)

    def test_evaluated_bridge_without_conflicting_paralog_keeps_positive_path(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out, frame, manifest = self._run(
                Path(tmp),
                bridge_rows="G_A\tFAM398\ta\ttier3_sequence_match_unique\n",
            )
            row = frame.loc["G_A"]
            self.assertEqual("preferred_copy", row["disambiguation_status"])
            self.assertEqual("route_a_blastp", row["copy_set_evidence_source"])
            self.assertEqual("", row["route_a_blast_guard"])
            shortlist = pd.read_csv(out / "family_copy_shortlist.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"G_A"}, set(shortlist["candidate_gene_id"]))
            self.assertTrue(manifest["route_a_blast"]["paralog_guard_evaluated"])


if __name__ == "__main__":
    unittest.main()
