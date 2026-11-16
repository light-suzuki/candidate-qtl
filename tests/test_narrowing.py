"""Issue #303: integrated narrowing (family recovery x copy x survival -> N1..N4).

Verifies the ordinal Wet-review-order layer: ranked columns pass through
verbatim, dry_score/dry_rank/dry_group/top_eligibility/screening_rank are never
recomputed, rank_contribution stays none, no numeric pseudo-probability is
produced, unknown gene ids fail closed, missing/empty inputs degrade to
not_evaluated (never a negative), tiers only move up when family information
exists, and the output is fully deterministic.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.narrowing import integrated_narrowing

RANKED = (
    "gene_id\tchrom\tstart\tend\tdry_score\tdry_rank\tdry_group\ttop_eligibility\tscreening_rank\trescue_reason\n"
    "g_top\tchr1\t1\t100\t8.0\t1\ttop\tstrong_genetic_evidence\t1\t\n"
    "g_reserve\tchr1\t101\t200\t5.0\t2\treserve\tmultiple_non_redundant_channels\t2\t\n"
    "g_rescue\tchr1\t201\t300\t3.0\t3\trescue\t\t3\t\n"
    "g_insufficient\tchr1\t301\t400\t1.0\t4\tevidence_insufficient\t\t4\t\n"
    "g_other\tchr1\t401\t500\t0.5\t5\t\t\t5\t\n"
)

RECOVERY = (
    "hypothesis_id\tgene_family\tmember_interval_genes\n"
    "cmp_001\tSBEI\tg_top;g_insufficient\n"
    "cmp_002\tRFO\tg_reserve\n"
)

CONVERGENCE = (
    "comparison_id\tfamily_id\tconvergence_status\tsource_dependence_status\tfour_state\n"
    "cmp_001\tSBEI\ta_b_convergent\tsource_independent\ta_and_b_supported\n"
    "cmp_002\tRFO\ta_only_supported\tnot_applicable\ta_only_supported\n"
)

COPY_DISCRIMINATION = (
    "family_id\tcandidate_gene_id\tdisambiguation_status\trank_contribution\n"
    "SBEI\tg_top\tpreferred_copy\tnone\n"
    "RFO\tg_reserve\tfamily_supported_copy_unresolved\tnone\n"
)

COPY_SHORTLIST = (
    "family_id\tcandidate_gene_id\tdisambiguation_status\tfirst_wet_priority\trank_contribution\n"
    "SBEI\tg_top\tpreferred_copy\tfirst\tnone\n"
)

REMAINING = (
    "candidate_gene_id\tcurrent_gene_state\tfirst_wet_priority\trank_contribution\n"
    "g_top\tsurvives\thigh\tnone\n"
    "g_reserve\tsurvives\tmedium\tnone\n"
    "g_insufficient\tweakened\tlow\tnone\n"
)


def _write(root: Path, name: str, content: str) -> Path:
    path = root / name
    path.write_text(content, encoding="utf-8")
    return path


def _run(root: Path, out: Path, *, ranked: str = RANKED, recovery: str | None = RECOVERY,
         convergence: str | None = CONVERGENCE, copy: str | None = COPY_DISCRIMINATION,
         shortlist: str | None = COPY_SHORTLIST, remaining: str | None = REMAINING,
         case_id: str = "case_303") -> dict[str, object]:
    ranked_path = _write(root, "ranked.tsv", ranked)
    args: dict[str, object] = {"ranked_arg": str(ranked_path), "out_dir_arg": str(out), "case_id": case_id}
    for name, content, arg in (
            ("recovery.tsv", recovery, "family_recovery_arg"),
            ("convergence.tsv", convergence, "convergence_arg"),
            ("copy.tsv", copy, "family_copy_arg"),
            ("shortlist.tsv", shortlist, "family_copy_shortlist_arg"),
            ("remaining.tsv", remaining, "remaining_arg")):
        if content is None:
            args[arg] = None
        else:
            args[arg] = str(_write(root, name, content))
    return integrated_narrowing(**args)


def _frame(root: Path, name: str = "integrated_narrowing.tsv") -> pd.DataFrame:
    return pd.read_csv(root / "out" / name, sep="\t", dtype=str).fillna("")


class TierAssignmentTests(unittest.TestCase):
    """Family (recovery) x copy x base rank combination matrix."""

    def test_n1_validated_family_preferred_copy_beats_rank(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _run(root, root / "out")
            row = _frame(root)[_frame(root)["gene_id"].eq("g_top")].iloc[0]
            self.assertEqual("N1", row["narrowing_tier"])
            self.assertEqual("1", row["narrowing_priority"])
            self.assertEqual("validated_family_preferred_copy", row["narrowing_reason"])
            self.assertEqual("", row["unresolved_reason"])

    def test_n2_validated_family_with_unresolved_copy(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _run(root, root / "out")
            row = _frame(root)[_frame(root)["gene_id"].eq("g_reserve")].iloc[0]
            self.assertEqual("N2", row["narrowing_tier"])
            self.assertEqual("validated_family_copy_unresolved", row["narrowing_reason"])
            self.assertIn("copy_state=family_supported_copy_unresolved", row["unresolved_reason"])

    def test_n2_family_recovered_but_copy_not_evaluated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # SBEI is recovered via convergence; g_insufficient is its member but
            # has no copy row -> family_recovered_copy_not_evaluated (N2).
            _run(root, root / "out")
            row = _frame(root)[_frame(root)["gene_id"].eq("g_insufficient")].iloc[0]
            self.assertEqual("N2", row["narrowing_tier"])
            self.assertEqual("family_recovered_copy_not_evaluated", row["narrowing_reason"])
            self.assertEqual("not_evaluated", row["copy_state"])
            self.assertEqual("strict_convergent", row["family_recovery_state"])

    def test_n3_base_top_reserve_rescue_without_family_info(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # No family, no copy, no remaining inputs at all: g_top (base top)
            # must land on N3, not N4 -- family absence is not a demotion.
            _run(root, root / "out", recovery=None, convergence=None, copy=None,
                 shortlist=None, remaining=None)
            frame = _frame(root)
            top = frame[frame["gene_id"].eq("g_top")].iloc[0]
            self.assertEqual("N3", top["narrowing_tier"])
            self.assertEqual("base_top_reserve_rescue", top["narrowing_reason"])
            self.assertEqual("not_evaluated", top["family_recovery_state"])
            self.assertEqual("not_evaluated", top["copy_state"])

    def test_n4_evidence_insufficient_and_base_rank_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _run(root, root / "out", recovery=None, convergence=None, copy=None,
                 shortlist=None)
            frame = _frame(root)
            insufficient = frame[frame["gene_id"].eq("g_insufficient")].iloc[0]
            self.assertEqual("N4", insufficient["narrowing_tier"])
            self.assertEqual("evidence_insufficient", insufficient["narrowing_reason"])
            other = frame[frame["gene_id"].eq("g_other")].iloc[0]
            self.assertEqual("N4", other["narrowing_tier"])
            self.assertEqual("base_rank_only", other["narrowing_reason"])

    def test_not_in_remaining_set_is_reported_for_eliminated_genes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # g_other is absent from the remaining set while the remaining input
            # is present and non-empty: it was eliminated by hard evidence.
            _run(root, root / "out")
            frame = _frame(root)
            row = frame[frame["gene_id"].eq("g_other")].iloc[0]
            self.assertEqual("not_in_remaining_set", row["survival_state"])
            row_top = frame[frame["gene_id"].eq("g_top")].iloc[0]
            self.assertEqual("survives", row_top["survival_state"])

    def test_first_wet_priority_precedence_shortlist_then_remaining(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _run(root, root / "out")
            frame = _frame(root)
            by_gene = {str(row["gene_id"]): row for row in frame.to_dict("records")}
            # shortlist wins for g_top ("first" over remaining's "high")
            self.assertEqual("first", by_gene["g_top"]["first_wet_priority"])
            # g_reserve has no shortlist row -> remaining's "medium"
            self.assertEqual("medium", by_gene["g_reserve"]["first_wet_priority"])

    def test_n2_all_validated_copy_statuses_stay_n2(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            copy = (
                "family_id\tcandidate_gene_id\tdisambiguation_status\trank_contribution\n"
                "SBEI\tg_top\tpreferred_copy\tnone\n"
                "RFO\tg_reserve\tco_preferred_tie\tnone\n"
                "FT\tg_rescue\ttandem_ambiguous\tnone\n"
                "CCAMK\tg_insufficient\tparalog_ambiguous\tnone\n"
                "NSP2\tg_other\tcopy_weakened\tnone\n"
            )
            _run(root, root / "out", copy=copy, recovery=None, convergence=None,
                 shortlist=None, remaining=None)
            frame = _frame(root)
            expected = {"g_top": "N1", "g_reserve": "N2", "g_rescue": "N2",
                        "g_insufficient": "N2", "g_other": "N2"}
            by_gene = {str(row["gene_id"]): row for row in frame.to_dict("records")}
            for gene, tier in expected.items():
                self.assertEqual(tier, by_gene[gene]["narrowing_tier"], gene)

    def test_copy_only_input_validates_family_for_n1(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # No convergence input at all: a preferred_copy still means the
            # family-to-copy layer validated the copy set, so N1 holds.
            _run(root, root / "out", recovery=None, convergence=None,
                 shortlist=None, remaining=None)
            row = _frame(root)[_frame(root)["gene_id"].eq("g_top")].iloc[0]
            self.assertEqual("N1", row["narrowing_tier"])
            self.assertEqual("validated_family_preferred_copy", row["narrowing_reason"])

    def test_gate_blocked_preferred_copy_is_not_n1(self) -> None:
        """Issue #374 review: a member whose discrimination status is
        preferred_copy but that is absent from the supplied shortlist (the
        promotion gate blocked it, e.g. a bridge tier5/tier6 tier) must never
        land on N1; it falls back to the base dry-group rules while the
        copy_state pass-through stays truthful."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            copy = (
                "family_id\tcandidate_gene_id\tdisambiguation_status\trank_contribution\n"
                "SBEI\tg_top\tpreferred_copy\tnone\n"
            )
            shortlist = (
                "family_id\tcandidate_gene_id\tdisambiguation_status\tfirst_wet_priority\trank_contribution\n"
            )
            _run(root, root / "out", recovery=None, convergence=None, copy=copy,
                 shortlist=shortlist, remaining=None)
            row = _frame(root)[_frame(root)["gene_id"].eq("g_top")].iloc[0]
            self.assertNotEqual("N1", row["narrowing_tier"])
            self.assertEqual("N3", row["narrowing_tier"])
            self.assertEqual("base_top_reserve_rescue", row["narrowing_reason"])
            self.assertEqual("preferred_copy", row["copy_state"])
            self.assertEqual("", row["first_wet_priority"])

    def test_gate_blocked_co_preferred_is_not_n2(self) -> None:
        """The reconciliation covers the validated-copy tier too: a member
        blocked out of a supplied shortlist cannot claim N2 from a promotion
        status and falls back to the base dry-group rules."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            copy = (
                "family_id\tcandidate_gene_id\tdisambiguation_status\trank_contribution\n"
                "SBEI\tg_top\tco_preferred_tie\tnone\n"
            )
            shortlist = (
                "family_id\tcandidate_gene_id\tdisambiguation_status\tfirst_wet_priority\trank_contribution\n"
            )
            _run(root, root / "out", recovery=None, convergence=None, copy=copy,
                 shortlist=shortlist, remaining=None)
            row = _frame(root)[_frame(root)["gene_id"].eq("g_top")].iloc[0]
            self.assertNotEqual("N2", row["narrowing_tier"])
            self.assertEqual("N3", row["narrowing_tier"])

    def test_global_narrowing_priority_spans_tiers(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _run(root, root / "out")
            frame = _frame(root).sort_values("narrowing_priority", kind="stable")
            priorities = frame["narrowing_priority"].astype(int).tolist()
            # 1..N over the whole table, tier-major: N1 (g_top), then N2 rows in
            # ranked order (g_reserve, g_insufficient), N3 (g_rescue), N4 (g_other).
            self.assertEqual(priorities, list(range(1, len(frame) + 1)))
            self.assertEqual(frame.iloc[0]["gene_id"], "g_top")
            self.assertEqual(frame.iloc[1]["gene_id"], "g_reserve")
            self.assertEqual(frame.iloc[2]["gene_id"], "g_insufficient")
            self.assertEqual(frame.iloc[3]["gene_id"], "g_rescue")
            self.assertEqual(frame.iloc[4]["gene_id"], "g_other")

    def test_multi_family_gene_uses_strongest_recovery_state(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            recovery = RECOVERY + "cmp_003\tFT\tg_top\n"
            convergence = CONVERGENCE + (
                "cmp_003\tFT\ta_only_supported\tnot_applicable\ta_only_supported\n")
            _run(root, root / "out", recovery=recovery, convergence=convergence)
            row = _frame(root)[_frame(root)["gene_id"].eq("g_top")].iloc[0]
            # SBEI is strict_convergent, FT only a_only: the strongest state wins
            # and the family ids are joined in sorted order.
            self.assertEqual("strict_convergent", row["family_recovery_state"])
            self.assertEqual("FT;SBEI", row["family_id"])


class RecoveryStateMappingTests(unittest.TestCase):
    def test_shared_source_dependent_convergence_never_elevates(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            convergence = (
                "comparison_id\tfamily_id\tconvergence_status\tsource_dependence_status\tfour_state\n"
                "cmp_001\tSBEI\ta_b_convergent\tshared_source_dependent\ta_and_b_supported\n"
            )
            _run(root, root / "out", convergence=convergence, copy=None,
                 shortlist=None, remaining=None)
            frame = _frame(root)
            row = frame[frame["gene_id"].eq("g_top")].iloc[0]
            self.assertEqual("a_and_b_supported", row["family_recovery_state"])
            # a_and_b_supported is a validated recovery state but without a
            # preferred copy the tier stays N2 (copy not evaluated).
            self.assertEqual("N2", row["narrowing_tier"])
            self.assertEqual("family_recovered_copy_not_evaluated", row["narrowing_reason"])

    def test_source_independent_convergence_elevates_to_strict(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _run(root, root / "out")
            frame = _frame(root)
            row = frame[frame["gene_id"].eq("g_top")].iloc[0]
            self.assertEqual("strict_convergent", row["family_recovery_state"])

    def test_four_state_a_only_and_b_only_and_neither(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            convergence = (
                "comparison_id\tfamily_id\tconvergence_status\tsource_dependence_status\tfour_state\n"
                "cmp_001\tSBEI\ta_only_supported\tnot_applicable\ta_only_supported\n"
                "cmp_002\tRFO\tb_only_supported\tnot_applicable\tb_only_supported\n"
                "cmp_003\tFT\tc_only_supported\tnot_applicable\tneither_supported\n"
            )
            recovery = (
                "hypothesis_id\tgene_family\tmember_interval_genes\n"
                "cmp_001\tSBEI\tg_top\ncmp_002\tRFO\tg_reserve\ncmp_003\tFT\tg_rescue\n"
            )
            _run(root, root / "out", recovery=recovery, convergence=convergence,
                 copy=None, shortlist=None, remaining=None)
            frame = _frame(root)
            by_gene = {str(row["gene_id"]): row for row in frame.to_dict("records")}
            self.assertEqual("a_only", by_gene["g_top"]["family_recovery_state"])
            self.assertEqual("N2", by_gene["g_top"]["narrowing_tier"])
            self.assertEqual("b_only", by_gene["g_reserve"]["family_recovery_state"])
            self.assertEqual("N2", by_gene["g_reserve"]["narrowing_tier"])
            self.assertEqual("none", by_gene["g_rescue"]["family_recovery_state"])
            # "none" is family-info-absent: not a demotion, base rescue stays N3.
            self.assertEqual("N3", by_gene["g_rescue"]["narrowing_tier"])

    def test_convergence_json_degraded_input(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = {
                "schema_version": "1.1",
                "convergence_source_accounting": [
                    {"family_id": "SBEI", "routes_converged": 2, "independent_sources": 1,
                     "shared_source_ids": ["PMID:1"], "source_dependence_status": "shared_source_dependent"},
                ],
            }
            convergence_path = root / "convergence.json"
            convergence_path.write_text(json.dumps(payload), encoding="utf-8")
            ranked_path = _write(root, "ranked.tsv", RANKED)
            recovery = _write(root, "recovery.tsv", RECOVERY)
            integrated_narrowing(str(ranked_path), str(root / "out"),
                                 family_recovery_arg=str(recovery),
                                 convergence_arg=str(convergence_path),
                                 case_id="case_303")
            frame = _frame(root)
            row = frame[frame["gene_id"].eq("g_top")].iloc[0]
            # JSON input cannot split a_only/b_only; shared-source convergence
            # is capped at a_and_b_supported.
            self.assertEqual("a_and_b_supported", row["family_recovery_state"])

    def test_family_absent_from_convergence_is_not_evaluated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            convergence = (
                "comparison_id\tfamily_id\tconvergence_status\tsource_dependence_status\tfour_state\n"
                "cmp_002\tOTHER\tc_only_supported\tnot_applicable\tneither_supported\n"
            )
            _run(root, root / "out", convergence=convergence, copy=None,
                 shortlist=None, remaining=None)
            frame = _frame(root)
            row = frame[frame["gene_id"].eq("g_top")].iloc[0]
            self.assertEqual("not_evaluated", row["family_recovery_state"])
            self.assertEqual("SBEI", row["family_id"])
            self.assertEqual("N3", row["narrowing_tier"])

    def test_legacy_comparison_without_four_state_falls_back_to_status(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            convergence = (
                "comparison_id\tfamily_id\tconvergence_status\tsource_dependence_status\n"
                "cmp_001\tSBEI\tphenotype_only_supported\tnot_applicable\n"
                "cmp_002\tRFO\tb_c_convergent\tsource_independent\n"
            )
            recovery = (
                "hypothesis_id\tgene_family\tmember_interval_genes\n"
                "cmp_001\tSBEI\tg_top\ncmp_002\tRFO\tg_reserve\n"
            )
            _run(root, root / "out", recovery=recovery, convergence=convergence,
                 copy=None, shortlist=None, remaining=None)
            frame = _frame(root)
            by_gene = {str(row["gene_id"]): row for row in frame.to_dict("records")}
            self.assertEqual("a_only", by_gene["g_top"]["family_recovery_state"])
            self.assertEqual("strict_convergent", by_gene["g_reserve"]["family_recovery_state"])


class PassthroughAndSchemaGuardTests(unittest.TestCase):
    def test_ranked_columns_are_passed_through_verbatim(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ranked = _write(root, "ranked.tsv", RANKED)
            _run(root, root / "out")
            produced = _frame(root)
            expected = pd.read_csv(ranked, sep="\t", dtype=str).fillna("")
            # Values are preserved per gene (rows are re-sorted by tier).
            expected_by_gene = {str(row["gene_id"]): row for row in expected.to_dict("records")}
            produced_by_gene = {str(row["gene_id"]): row for row in produced.to_dict("records")}
            for gene_id, row in produced_by_gene.items():
                for column in expected.columns:
                    self.assertEqual(expected_by_gene[gene_id][column], row[column], f"{gene_id}.{column}")
            self.assertEqual(["gene_id", "chrom", "start", "end", "dry_score", "dry_rank",
                              "dry_group", "top_eligibility", "screening_rank", "rescue_reason",
                              "family_recovery_state", "family_id", "copy_state",
                              "first_wet_priority", "positional_constraint", "survival_state",
                              "narrowing_tier", "narrowing_priority", "narrowing_reason",
                              "unresolved_reason", "rank_contribution"],
                             list(produced.columns))

    def test_no_numeric_pseudo_probability_and_rank_contribution_none(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _run(root, root / "out")
            frame = _frame(root)
            self.assertEqual({"none"}, set(frame["rank_contribution"]))
            for column in frame.columns:
                self.assertNotIn("prob", column.casefold(), column)

    def test_deterministic_byte_identical_repeat(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ranked = _write(root, "ranked.tsv", RANKED)
            recovery = _write(root, "recovery.tsv", RECOVERY)
            convergence = _write(root, "convergence.tsv", CONVERGENCE)
            copy = _write(root, "copy.tsv", COPY_DISCRIMINATION)
            shortlist = _write(root, "shortlist.tsv", COPY_SHORTLIST)
            remaining = _write(root, "remaining.tsv", REMAINING)
            integrated_narrowing(str(ranked), str(root / "out"),
                                 family_recovery_arg=str(recovery),
                                 convergence_arg=str(convergence),
                                 family_copy_arg=str(copy),
                                 family_copy_shortlist_arg=str(shortlist),
                                 remaining_arg=str(remaining), case_id="case_303")
            first_tsv = (root / "out" / "integrated_narrowing.tsv").read_bytes()
            first_manifest = (root / "out" / "narrowing_manifest.json").read_bytes()
            # Re-run into the same out-dir and compare against the snapshot:
            # a self-comparison of the file with itself cannot detect
            # non-determinism.  (The manifest records absolute output paths, so
            # two different out-dirs are not byte-comparable.)
            integrated_narrowing(str(ranked), str(root / "out"),
                                 family_recovery_arg=str(recovery),
                                 convergence_arg=str(convergence),
                                 family_copy_arg=str(copy),
                                 family_copy_shortlist_arg=str(shortlist),
                                 remaining_arg=str(remaining), case_id="case_303")
            self.assertEqual(first_tsv, (root / "out" / "integrated_narrowing.tsv").read_bytes())
            self.assertEqual(first_manifest, (root / "out" / "narrowing_manifest.json").read_bytes())

    def test_ranked_collision_with_narrowing_columns_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ranked = (
                "gene_id\tdry_group\tfamily_id\n"
                "g_top\ttop\tFT\n"
                "g_reserve\treserve\t\n"
            )
            with self.assertRaises(ValueError) as raised:
                _run(root, root / "out", ranked=ranked, recovery=None, convergence=None,
                     copy=None, shortlist=None, remaining=None)
            self.assertIn("narrowing-owned columns", str(raised.exception))
            self.assertIn("family_id", str(raised.exception))

    def test_ranked_with_internal_sentinel_column_fails_closed(self) -> None:
        # ``_original_order`` is the internal row-order sentinel; a ranked
        # table that already carries it must fail closed instead of crashing
        # the output shaping with a raw KeyError.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ranked = (
                "gene_id\tdry_group\t_original_order\n"
                "g_top\ttop\t1\n"
                "g_reserve\treserve\t2\n"
            )
            with self.assertRaises(ValueError) as raised:
                _run(root, root / "out", ranked=ranked, recovery=None, convergence=None,
                     copy=None, shortlist=None, remaining=None)
            self.assertIn("_original_order", str(raised.exception))

    def test_positional_constraint_passthrough_when_ranked_carries_it(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ranked = (
                "gene_id\tdry_group\tpositional_constraint\n"
                "g_top\ttop\tcustom_constraint\n"
                "g_reserve\tevidence_insufficient\tcustom_too\n"
            )
            _run(root, root / "out", ranked=ranked, recovery=None, convergence=None,
                 copy=None, shortlist=None, remaining=None)
            frame = _frame(root)
            by_gene = {str(row["gene_id"]): row for row in frame.to_dict("records")}
            self.assertEqual("custom_constraint", by_gene["g_top"]["positional_constraint"])
            self.assertEqual("custom_too", by_gene["g_reserve"]["positional_constraint"])

    def test_original_row_order_is_preserved_within_tier(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # g_insufficient appears before g_other in the ranked input; both are
            # N4 here (no family inputs), so the narrowing order must keep the
            # original relative order.
            _run(root, root / "out", recovery=None, convergence=None, copy=None,
                 shortlist=None, remaining=None)
            frame = _frame(root)
            n4 = frame[frame["narrowing_tier"].eq("N4")]["gene_id"].tolist()
            self.assertEqual(["g_insufficient", "g_other"], n4)

    def test_manifest_frozen_text_and_safety(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _run(root, root / "out")
            manifest = json.loads((root / "out" / "narrowing_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("1.0", manifest["schema_version"])
            self.assertIn("never excludes candidates", manifest["safety"][0])
            self.assertIn("(elimination alone may exclude)", manifest["safety"][0])
            self.assertIn("ordinal tiers only", manifest["rank_contribution"])
            self.assertIn("N1", manifest["tier_rule"])
            self.assertIn("strict_convergent", manifest["recovery_state_mapping"])
            self.assertIn("shared-source", manifest["recovery_state_mapping"])
            self.assertIn("ranked", manifest["inputs"])
            self.assertIn("sha256", manifest["inputs"]["ranked"])


class FailClosedTests(unittest.TestCase):
    def test_unknown_gene_id_in_family_copy_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            copy = COPY_DISCRIMINATION + "OTHER\tghost\tpreferred_copy\tnone\n"
            with self.assertRaises(ValueError) as raised:
                _run(root, root / "out", copy=copy)
            self.assertIn("fail-closed", str(raised.exception))
            self.assertIn("ghost", str(raised.exception))

    def test_unknown_gene_id_in_remaining_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            remaining = REMAINING + "ghost\tsurvives\thigh\tnone\n"
            with self.assertRaises(ValueError) as raised:
                _run(root, root / "out", remaining=remaining)
            self.assertIn("fail-closed", str(raised.exception))

    def test_unknown_gene_id_in_family_recovery_members_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            recovery = RECOVERY + "cmp_003\tFT\tghost\n"
            with self.assertRaises(ValueError) as raised:
                _run(root, root / "out", recovery=recovery)
            self.assertIn("fail-closed", str(raised.exception))

    def test_remaining_without_current_gene_state_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            remaining = "candidate_gene_id\n" + "\n".join(
                line.split("\t")[0] for line in REMAINING.splitlines()[1:]) + "\n"
            with self.assertRaises(ValueError) as raised:
                _run(root, root / "out", remaining=remaining)
            self.assertIn("current_gene_state", str(raised.exception))

    def test_ranked_without_gene_id_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaises(ValueError):
                _run(root, root / "out", ranked="name\tvalue\nfoo\t1\n")

    def test_empty_optional_inputs_degrade_to_not_evaluated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            empty = "hypothesis_id\tgene_family\tmember_interval_genes\n"
            _run(root, root / "out", recovery=empty, convergence=None, copy=None,
                 shortlist=None, remaining=None)
            frame = _frame(root)
            for row in frame.to_dict("records"):
                self.assertEqual("not_evaluated", row["family_recovery_state"])
                self.assertEqual("not_evaluated", row["copy_state"])
            self.assertEqual({"N3", "N4"}, set(frame["narrowing_tier"]))


class CliCommandTests(unittest.TestCase):
    def test_integrated_narrowing_command_registered_and_runs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ranked = _write(root, "ranked.tsv", RANKED)
            recovery = _write(root, "recovery.tsv", RECOVERY)
            convergence = _write(root, "convergence.tsv", CONVERGENCE)
            copy = _write(root, "copy.tsv", COPY_DISCRIMINATION)
            out = root / "out"
            self.assertEqual(0, main(["integrated-narrowing", "--ranked", str(ranked),
                                      "--out-dir", str(out),
                                      "--family-recovery", str(recovery),
                                      "--convergence", str(convergence),
                                      "--family-copy", str(copy),
                                      "--case-id", "cli_case"]))
            frame = pd.read_csv(out / "integrated_narrowing.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(5, len(frame))
            self.assertTrue((out / "narrowing_manifest.json").exists())


if __name__ == "__main__":
    unittest.main()
