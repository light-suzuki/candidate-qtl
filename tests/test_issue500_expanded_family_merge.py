"""Regression coverage for issue #500 bounded expanded A/B family matching."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.narrowing import _recovery_state_from_comparison
from candidate_qtl.retrieval_route import RETRIEVAL_COLUMNS, triple_route_compare


def _write(root: Path, name: str, content: str) -> Path:
    path = root / name
    path.write_text(content, encoding="utf-8")
    return path


def _route_a(label: str) -> str:
    return (
        "hypothesis_id\tgene_family\tmember_interval_genes\tsource_backing\n"
        f"a1\t{label}\tg1\tPMID:100\n"
    )


def _route_b(label: str) -> str:
    return (
        "hypothesis_id\tfamily_id\tgene_family\tmember_interval_genes\tsource_backing\n"
        f"b1\t{label}\t{label}\tg1\tPMID:200\n"
    )


def _empty_route_c() -> str:
    return "\t".join(RETRIEVAL_COLUMNS) + "\n"


def _run(route_a_label: str, route_b_label: str) -> tuple[pd.DataFrame, dict[str, object], pd.DataFrame]:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        a = _write(root, "a.tsv", _route_a(route_a_label))
        b = _write(root, "b.tsv", _route_b(route_b_label))
        c = _write(root, "c.tsv", _empty_route_c())
        out = root / "out"
        triple_route_compare(str(a), str(b), str(c), str(out), case_id="issue500")
        comparison = pd.read_csv(
            out / "triple_route_hypothesis_comparison.tsv", sep="\t", dtype=str
        ).fillna("")
        manifest = json.loads((out / "triple_route_convergence.json").read_text(encoding="utf-8"))
        handoff = pd.read_csv(out / "family_to_copy_input.tsv", sep="\t", dtype=str).fillna("")
        return comparison.copy(), manifest, handoff.copy()


class Issue500ExpandedFamilyMergeTests(unittest.TestCase):
    def test_descriptive_route_a_label_meets_symbol_as_non_strict_match(self) -> None:
        comparison, manifest, handoff = _run("CIK LRR receptor kinases", "CIK1")
        self.assertEqual(1, len(comparison))
        row = comparison.iloc[0]
        self.assertEqual("a1", row["route_a_hypothesis_ids"])
        self.assertEqual("b1", row["route_b_hypothesis_ids"])
        self.assertEqual("expanded_symbol_token", row["family_match_basis"])
        self.assertEqual("a_b_expanded_match", row["convergence_status"])
        self.assertEqual("a_and_b_supported", row["four_state"])
        self.assertEqual("none", row["rank_contribution"])
        self.assertEqual("source_independent", row["source_dependence_status"])

        # Expanded agreement is intentionally weaker than strict convergence.
        self.assertEqual("a_and_b_supported", _recovery_state_from_comparison(row.to_dict()))

        matching = manifest["family_matching"]
        self.assertEqual(1, matching["expanded_match_count"])
        self.assertEqual(1, matching["basis_counts"]["expanded_symbol_token"])
        self.assertEqual("g1", handoff.iloc[0]["member_interval_genes"])

    def test_exact_family_match_keeps_existing_strict_convergence(self) -> None:
        comparison, manifest, _handoff = _run("CIK1", "CIK1")
        self.assertEqual(1, len(comparison))
        row = comparison.iloc[0]
        self.assertEqual("strict_key", row["family_match_basis"])
        self.assertEqual("a_b_convergent", row["convergence_status"])
        self.assertEqual("strict_convergent", _recovery_state_from_comparison(row.to_dict()))
        self.assertEqual(0, manifest["family_matching"]["expanded_match_count"])


if __name__ == "__main__":
    unittest.main()
