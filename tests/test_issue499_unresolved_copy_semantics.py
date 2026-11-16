"""Regression coverage for issue #499 copy-state vocabulary handling."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.family_to_copy import DISAMBIGUATION_STATES
from candidate_qtl.narrowing import (
    RECOVERY_CONTEXT_COPY_UNRESOLVED_STATUSES,
    VALIDATED_COPY_STATUSES,
    integrated_narrowing,
)


RANKED = (
    "gene_id\tdry_group\n"
    "g1\tevidence_insufficient\n"
)
RECOVERY = (
    "hypothesis_id\tgene_family\tmember_interval_genes\n"
    "h1\tFAM\tg1\n"
)
CONVERGENCE = (
    "comparison_id\tfamily_id\tconvergence_status\tsource_dependence_status\tfour_state\n"
    "c1\tFAM\tphenotype_only_supported\tnot_applicable\ta_only_supported\n"
)


def _write(root: Path, name: str, content: str) -> Path:
    path = root / name
    path.write_text(content, encoding="utf-8")
    return path


def _run(copy_state: str, *, recovered: bool = True,
         bridge_mapping_tier: str | None = None,
         shortlist_supplied: bool = False) -> pd.Series:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        ranked = _write(root, "ranked.tsv", RANKED)
        copy_header = "family_id\tcandidate_gene_id\tdisambiguation_status\trank_contribution"
        copy_row = f"FAM\tg1\t{copy_state}\tnone"
        if bridge_mapping_tier is not None:
            copy_header += "\tbridge_mapping_tier"
            copy_row += f"\t{bridge_mapping_tier}"
        copy = _write(root, "copy.tsv", f"{copy_header}\n{copy_row}\n")
        kwargs: dict[str, str] = {
            "ranked_arg": str(ranked),
            "out_dir_arg": str(root / "out"),
            "family_copy_arg": str(copy),
            "case_id": "issue499",
        }
        if shortlist_supplied:
            kwargs["family_copy_shortlist_arg"] = str(
                _write(root, "shortlist.tsv", "candidate_gene_id\n"))
        if recovered:
            kwargs["family_recovery_arg"] = str(_write(root, "recovery.tsv", RECOVERY))
            kwargs["convergence_arg"] = str(_write(root, "convergence.tsv", CONVERGENCE))
        integrated_narrowing(**kwargs)
        frame = pd.read_csv(root / "out" / "integrated_narrowing.tsv", sep="\t", dtype=str).fillna("")
        return frame.iloc[0].copy()


class Issue499UnresolvedCopySemanticsTests(unittest.TestCase):
    def test_family_to_copy_vocabulary_has_explicit_narrowing_semantics(self) -> None:
        handled = (
            set(VALIDATED_COPY_STATUSES)
            | set(RECOVERY_CONTEXT_COPY_UNRESOLVED_STATUSES)
            | {"not_in_validated_copy_set"}
        )
        self.assertEqual(set(DISAMBIGUATION_STATES), handled)
        self.assertEqual(
            {"insufficient_discriminating_evidence", "family_resolution_required"},
            set(RECOVERY_CONTEXT_COPY_UNRESOLVED_STATUSES),
        )

    def test_evaluated_but_unresolved_states_preserve_recovered_family_n2(self) -> None:
        for state in sorted(RECOVERY_CONTEXT_COPY_UNRESOLVED_STATUSES):
            with self.subTest(state=state):
                row = _run(state, recovered=True)
                self.assertEqual("N2", row["narrowing_tier"])
                self.assertEqual("family_recovered_copy_unresolved", row["narrowing_reason"])
                self.assertEqual(state, row["copy_state"])
                self.assertIn("unresolved is not negative evidence", row["unresolved_reason"])

    def test_unresolved_context_state_does_not_validate_family_by_itself(self) -> None:
        for state in sorted(RECOVERY_CONTEXT_COPY_UNRESOLVED_STATUSES):
            with self.subTest(state=state):
                row = _run(state, recovered=False)
                self.assertEqual("N4", row["narrowing_tier"])
                self.assertEqual("evidence_insufficient", row["narrowing_reason"])
                self.assertEqual("not_evaluated", row["family_recovery_state"])

    def test_not_in_validated_copy_set_does_not_inherit_n2_from_coarse_family(self) -> None:
        row = _run("not_in_validated_copy_set", recovered=True)
        self.assertEqual("N4", row["narrowing_tier"])
        self.assertEqual("evidence_insufficient", row["narrowing_reason"])
        self.assertEqual("not_in_validated_copy_set", row["copy_state"])
        self.assertEqual("a_only", row["family_recovery_state"])

    def test_weak_bridge_does_not_promote_name_only_unresolved_copy(self) -> None:
        for tier in ("tier5_name_only", "tier6_not_evaluated"):
            with self.subTest(tier=tier):
                row = _run("family_resolution_required", recovered=True,
                           bridge_mapping_tier=tier, shortlist_supplied=True)
                self.assertEqual("N4", row["narrowing_tier"])
                self.assertEqual("evidence_insufficient", row["narrowing_reason"])
                self.assertEqual("a_only", row["family_recovery_state"])

    def test_source_backed_bridge_preserves_canonical_a_only_recovery(self) -> None:
        row = _run("family_resolution_required", recovered=True,
                   bridge_mapping_tier="tier1_source_backed_unique",
                   shortlist_supplied=True)
        self.assertEqual("N2", row["narrowing_tier"])
        self.assertEqual("family_recovered_copy_unresolved", row["narrowing_reason"])


if __name__ == "__main__":
    unittest.main()
