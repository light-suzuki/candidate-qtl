"""Issue #426: the A/B -> Wet-candidate handoff and its bounded adaptive round.

The contract these tests hold in place is a cost contract as much as a
scientific one: a candidate set that already fits the Wet budget must not pay
for an adaptive/LLM stage, and a set that does not fit must not be collapsed to
fit.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.adaptive_synthesis import (
    apply_verified_refinement,
    build_adaptive_synthesis_request,
    import_adaptive_synthesis_proposals,
)
from candidate_qtl.candidate_budget import candidate_budget_gate
from candidate_qtl.supported_member_union import build_supported_member_union
from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from tests._issue664_helpers import manifest_for


def _candidates(root: Path, count: int, name: str = "candidates.tsv") -> Path:
    path = root / name
    pd.DataFrame(
        [{"candidate_gene_id": f"G{index:03d}", "route_support": "route_a"} for index in range(1, count + 1)]
    ).to_csv(path, sep="\t", index=False, lineterminator="\n")
    return path


def _gate(root: Path, candidates: Path, *, wet_budget: int, adaptive_round: int = 0,
          max_adaptive_rounds: int = 1, out: str = "gate") -> dict:
    return candidate_budget_gate(
        str(candidates), str(root / out),
        wet_budget=wet_budget, adaptive_round=adaptive_round,
        budget_policy="strict",
        max_adaptive_rounds=max_adaptive_rounds, case_id="case1",
    )


def _supported_union(root: Path, count: int, *, case_id: str = "case1") -> dict[str, object]:
    candidates = root / "route_b_candidates.tsv"
    pd.DataFrame([
        {
            "candidate_gene_id": f"G{index:03d}",
            "route": "b",
            "support_status": "verified_semantic_support",
        }
        for index in range(1, count + 1)
    ]).to_csv(candidates, sep="\t", index=False, lineterminator="\n")
    manifest = manifest_for(candidates, route="b", case_id=case_id)
    return build_supported_member_union(
        str(root / "original_union"), route_b_candidates_arg=str(candidates), route_b_candidates_manifest_arg=str(manifest), case_id=case_id,
    )


class WetBudgetStopRuleTests(unittest.TestCase):
    def test_within_budget_never_reaches_the_adaptive_stage(self) -> None:
        """The cheap-first stop rule: 16 supported members need no LLM call."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = _candidates(root, 16)
            gate = _gate(root, candidates, wet_budget=20)
            self.assertEqual("within_wet_budget", gate["status"])
            self.assertFalse(gate["adaptive_refinement_required"])

            with self.assertRaises(ValueError) as ctx:
                build_adaptive_synthesis_request(
                    str(candidates), str(gate["manifest"]), str(root / "asyn"), case_id="case1",
                )
            self.assertIn("only run when the candidate budget gate requested it", str(ctx.exception))

    def test_over_budget_requests_exactly_one_bounded_round(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = _candidates(root, 30)
            first = _gate(root, candidates, wet_budget=20, adaptive_round=0)
            self.assertEqual("adaptive_refinement_required", first["status"])

            # Round budget spent: the same over-budget set now stops instead of
            # asking for another round.
            second = _gate(root, candidates, wet_budget=20, adaptive_round=1, out="gate2")
            self.assertEqual("review_required_over_budget", second["status"])
            with self.assertRaises(ValueError):
                build_adaptive_synthesis_request(
                    str(candidates), str(second["manifest"]), str(root / "asyn2"), case_id="case1",
                )


class AdaptiveSynthesisProposalTests(unittest.TestCase):
    def _request(self, root: Path, count: int = 30) -> tuple[Path, dict]:
        candidates = _candidates(root, count)
        gate = _gate(root, candidates, wet_budget=20)
        result = build_adaptive_synthesis_request(
            str(candidates), str(gate["manifest"]), str(root / "asyn"), case_id="case1",
        )
        return candidates, result

    def test_request_asks_for_ambiguity_not_a_gene_selection(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _candidates_path, result = self._request(root)
            request = json.loads(Path(result["request"]).read_text(encoding="utf-8"))
            prompt = request["system_prompt"]
            self.assertIn("Do NOT choose which genes to keep", prompt)
            self.assertIn("Do NOT rank genes", prompt)
            self.assertIn("numerical score", prompt)
            self.assertEqual(
                "candidate_count_exceeds_wet_budget", result["invocation_reason"]
            )

    def test_missing_provider_is_a_state_not_a_failure(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _candidates_path, result = self._request(root)
            imported = import_adaptive_synthesis_proposals(
                None, str(result["request"]), str(root / "imp"), case_id="case1",
            )
            self.assertEqual("not_evaluated_no_provider", imported["provider_state"])
            self.assertEqual(0, imported["proposal_count"])
            frame = pd.read_csv(imported["proposals"], sep="\t", dtype=str)
            self.assertEqual(0, len(frame))

    def test_proposal_carrying_a_score_is_rejected(self) -> None:
        """A scored proposal is a ranking attempt in proposal clothing."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _candidates_path, result = self._request(root)
            scored = root / "scored.tsv"
            pd.DataFrame([{
                "proposal_id": "p1",
                "ambiguity_kind": "paralog_copy",
                "ambiguity_target": "FAM1",
                "affected_candidate_genes": "G001;G002",
                "discriminating_evidence_request": "RBH against the reference copy",
                "evidence_kind": "sequence_check",
                "score": "0.91",
            }]).to_csv(scored, sep="\t", index=False, lineterminator="\n")
            with self.assertRaises(ValueError) as ctx:
                import_adaptive_synthesis_proposals(
                    str(scored), str(result["request"]), str(root / "imp"), case_id="case1",
                )
            self.assertIn("never scores genes", str(ctx.exception))

    def test_proposal_cannot_introduce_a_gene_outside_the_request(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _candidates_path, result = self._request(root)
            smuggled = root / "smuggled.tsv"
            pd.DataFrame([{
                "proposal_id": "p1",
                "ambiguity_kind": "family_vs_subfamily",
                "ambiguity_target": "FAM1",
                "affected_candidate_genes": "G001;NOT_IN_SET",
                "discriminating_evidence_request": "subfamily HMM",
                "evidence_kind": "sequence_check",
            }]).to_csv(smuggled, sep="\t", index=False, lineterminator="\n")
            with self.assertRaises(ValueError) as ctx:
                import_adaptive_synthesis_proposals(
                    str(smuggled), str(result["request"]), str(root / "imp"), case_id="case1",
                )
            self.assertIn("outside the requested candidate set", str(ctx.exception))

    def test_accepted_proposal_stays_non_scoring_and_review_required(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _candidates_path, result = self._request(root)
            good = root / "good.tsv"
            pd.DataFrame([{
                "proposal_id": "p1",
                "ambiguity_kind": "paralog_copy",
                "ambiguity_target": "FAM1",
                "affected_candidate_genes": "G001;G002",
                "discriminating_evidence_request": "1:1 RBH against the reference copy",
                "evidence_kind": "sequence_check",
                "expected_effect": "collapses two copies to one only if RBH is unambiguous",
                "competing_explanations_to_keep": "both copies if RBH stays ambiguous",
            }]).to_csv(good, sep="\t", index=False, lineterminator="\n")
            imported = import_adaptive_synthesis_proposals(
                str(good), str(result["request"]), str(root / "imp"),
                provider="fixture", model="fixture", case_id="case1",
            )
            frame = pd.read_csv(imported["proposals"], sep="\t", dtype=str).fillna("")
            self.assertEqual(1, len(frame))
            self.assertEqual("none", frame.iloc[0]["rank_contribution"])
            self.assertEqual("yes", frame.iloc[0]["review_required"])


class VerifiedRefinementTests(unittest.TestCase):
    def test_proposals_alone_never_shrink_the_set(self) -> None:
        """Only a machine verdict narrows; the proposal table is not evidence."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = _candidates(root, 30)
            refined = apply_verified_refinement(
                str(candidates), None, str(root / "ref"), case_id="case1",
            )
            self.assertEqual(30, refined["output_candidate_count"])
            self.assertEqual("not_evaluated_no_verification", refined["verification_state"])

    def test_not_evaluated_verdict_keeps_the_candidate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = _candidates(root, 5)
            verification = root / "verify.tsv"
            pd.DataFrame([
                {"candidate_gene_id": "G001", "verdict": "not_evaluated", "verification_source": ""},
                {"candidate_gene_id": "G002", "verdict": "support_confirmed", "verification_source": "PMID:1"},
            ]).to_csv(verification, sep="\t", index=False, lineterminator="\n")
            refined = apply_verified_refinement(
                str(candidates), str(verification), str(root / "ref"), case_id="case1",
            )
            self.assertEqual(5, refined["output_candidate_count"])
            self.assertEqual(0, refined["removed_candidate_count"])

    def test_unsourced_falsification_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = _candidates(root, 5)
            verification = root / "verify.tsv"
            pd.DataFrame([
                {"candidate_gene_id": "G001", "verdict": "support_falsified", "verification_source": ""},
            ]).to_csv(verification, sep="\t", index=False, lineterminator="\n")
            with self.assertRaises(ValueError) as ctx:
                apply_verified_refinement(
                    str(candidates), str(verification), str(root / "ref"), case_id="case1",
                )
            self.assertIn("requires a verification_source", str(ctx.exception))

    def test_sourced_falsification_removes_only_that_candidate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = _candidates(root, 5)
            verification = root / "verify.tsv"
            pd.DataFrame([
                {"candidate_gene_id": "G001", "verdict": "support_falsified",
                 "verification_source": "RBH:one_to_many_resolved"},
            ]).to_csv(verification, sep="\t", index=False, lineterminator="\n")
            refined = apply_verified_refinement(
                str(candidates), str(verification), str(root / "ref"), case_id="case1",
            )
            self.assertEqual(4, refined["output_candidate_count"])
            kept = pd.read_csv(refined["refined_candidate_set"], sep="\t", dtype=str)
            self.assertNotIn("G001", set(kept["candidate_gene_id"]))
            self.assertIn("G002", set(kept["candidate_gene_id"]))

    def test_wet_reentry_requires_verified_lineage(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            arbitrary = _candidates(root, 1)
            with self.assertRaises(ValueError) as ctx:
                wet_candidate_handoff(
                    str(root / "wet"), refined_candidates_arg=str(arbitrary),
                    wet_budget=20, adaptive_round=1, case_id="case1",
                )
            self.assertIn("requires refinement_manifest_arg", str(ctx.exception))

    def test_valid_sourced_refinement_reenters_and_preserves_lineage(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            union = _supported_union(root, 30)
            original = Path(str(union["candidate_set"]))
            verification = root / "verification.tsv"
            pd.DataFrame([
                {
                    "candidate_gene_id": f"G{index:03d}",
                    "verdict": "support_falsified",
                    "verification_source": f"RBH:resolved:{index}",
                }
                for index in range(21, 31)
            ]).to_csv(verification, sep="\t", index=False, lineterminator="\n")
            refined = apply_verified_refinement(
                str(original), str(verification), str(root / "refined"), case_id="case1",
            )
            result = wet_candidate_handoff(
                str(root / "wet"),
                refined_candidates_arg=refined["refined_candidate_set"],
                refinement_manifest_arg=refined["manifest"],
                refinement_original_candidates_arg=str(original),
                refinement_original_union_manifest_arg=str(union["manifest"]),
                wet_budget=20,
                adaptive_round=1,
                case_id="case1",
            )
            self.assertEqual("wet_reviewable", result["handoff_state"])
            self.assertEqual(20, result["candidate_count"])
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual("verified_lineage", manifest["refinement_lineage"]["verification_state"])
            self.assertEqual(10, manifest["refinement_lineage"]["removed_candidate_count"])

    def test_tampered_refined_output_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            union = _supported_union(root, 2)
            original = Path(str(union["candidate_set"]))
            verification = root / "verification.tsv"
            pd.DataFrame([{
                "candidate_gene_id": "G002",
                "verdict": "support_falsified",
                "verification_source": "RBH:resolved",
            }]).to_csv(verification, sep="\t", index=False, lineterminator="\n")
            refined = apply_verified_refinement(
                str(original), str(verification), str(root / "refined"), case_id="case1",
            )
            pd.DataFrame([{
                "candidate_gene_id": "OUTSIDE_INTERVAL",
                "route_support": "route_a",
            }]).to_csv(refined["refined_candidate_set"], sep="\t", index=False, lineterminator="\n")
            with self.assertRaises(ValueError) as ctx:
                wet_candidate_handoff(
                    str(root / "wet"),
                    refined_candidates_arg=refined["refined_candidate_set"],
                    refinement_manifest_arg=refined["manifest"],
                    refinement_original_candidates_arg=str(original),
                    refinement_original_union_manifest_arg=str(union["manifest"]),
                    wet_budget=20,
                    adaptive_round=1,
                    case_id="case1",
                )
            self.assertIn("output hash mismatch", str(ctx.exception))

    def test_internally_consistent_but_noncanonical_original_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            original = _candidates(root, 21)
            verification = root / "verification.tsv"
            pd.DataFrame([{
                "candidate_gene_id": "G021",
                "verdict": "support_falsified",
                "verification_source": "RBH:resolved",
            }]).to_csv(verification, sep="\t", index=False, lineterminator="\n")
            refined = apply_verified_refinement(
                str(original), str(verification), str(root / "refined"), case_id="case1",
            )
            fake_union_manifest = root / "fake_union_manifest.json"
            fake_union_manifest.write_text(json.dumps({
                "stage": "supported_member_union",
                "case_id": "case1",
                "output_sha256": "fabricated",
            }), encoding="utf-8")
            with self.assertRaises(ValueError) as ctx:
                wet_candidate_handoff(
                    str(root / "wet"),
                    refined_candidates_arg=refined["refined_candidate_set"],
                    refinement_manifest_arg=refined["manifest"],
                    refinement_original_candidates_arg=str(original),
                    refinement_original_union_manifest_arg=str(fake_union_manifest),
                    wet_budget=20,
                    adaptive_round=1,
                    case_id="case1",
                )
            self.assertIn("canonical supported-union contract", str(ctx.exception))


class HandoffChainTests(unittest.TestCase):
    def _bridge(self, root: Path, count: int) -> Path:
        path = root / "bridge.tsv"
        pd.DataFrame([
            {
                "route": "a",
                "gene_id": f"G{index:03d}",
                "family_id": "FAM1",
                "hypothesis_id": "H1",
                "mapping_tier": "tier3_family_projection",
                "ambiguity_status": "unique",
                "review_status": "supported",
            }
            for index in range(1, count + 1)
        ]).to_csv(path, sep="\t", index=False, lineterminator="\n")
        manifest_for(path, route="a", case_id="case1")
        return path

    def test_within_budget_handoff_is_wet_reviewable_without_heavy_stage(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = wet_candidate_handoff(
                str(root / "out"),
                route_a_bridge_arg=str(self._bridge(root, 8)),
                wet_budget=20,
                case_id="case1",
            )
            self.assertEqual("wet_reviewable", result["handoff_state"])
            self.assertFalse(result["heavy_stage_invoked"])
            self.assertEqual("not_invoked_within_wet_budget", result["heavy_stage_invocation_reason"])
            self.assertIsNone(result["adaptive_synthesis_request"])

    def test_over_budget_handoff_records_the_invocation_reason(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = wet_candidate_handoff(
                str(root / "out"),
                budget_policy="strict",
                route_a_bridge_arg=str(self._bridge(root, 30)),
                wet_budget=20,
                case_id="case1",
            )
            self.assertEqual("adaptive_refinement_requested", result["handoff_state"])
            self.assertTrue(result["heavy_stage_invoked"])
            self.assertEqual(
                "candidate_count_exceeds_wet_budget", result["heavy_stage_invocation_reason"]
            )
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual(30, manifest["candidate_count"])
            self.assertEqual(20, manifest["wet_budget"])

    def test_exhausted_rounds_stop_as_review_required_not_a_forced_top_n(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = wet_candidate_handoff(
                str(root / "out"),
                budget_policy="strict",
                route_a_bridge_arg=str(self._bridge(root, 30)),
                wet_budget=20,
                adaptive_round=1,
                max_adaptive_rounds=1,
                case_id="case1",
            )
            self.assertEqual("review_required_over_budget", result["handoff_state"])
            self.assertFalse(result["heavy_stage_invoked"])
            # The set is still 30: nothing was truncated to satisfy the budget.
            self.assertEqual(30, result["candidate_count"])

    def test_handoff_state_is_the_verdict_not_the_exit_code(self) -> None:
        """Every terminal state is explicit, including the abstentions."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = wet_candidate_handoff(
                str(root / "out"),
                route_a_bridge_arg=str(self._bridge(root, 30)),
                wet_budget=20,
                case_id="case1",
            )
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertIn(manifest["handoff_state"], {
                "wet_reviewable",
                "adaptive_refinement_requested",
                "review_required_over_budget",
                "review_required_no_supported_candidates",
            })
            self.assertEqual("none", manifest["rank_contribution"])


class DeterministicSubfamilyRefinementTests(unittest.TestCase):
    """One bounded, machine-verified subfamily pass, then re-gate (#426)."""

    def _bridge(self, root: Path, count: int) -> Path:
        path = root / "bridge.tsv"
        pd.DataFrame([
            {
                "route": "a",
                "gene_id": f"G{index:03d}",
                "family_id": "FAM1",
                "hypothesis_id": "H1",
                "mapping_tier": "tier3_family_projection",
                "ambiguity_status": "unique",
                "review_status": "supported",
            }
            for index in range(1, count + 1)
        ]).to_csv(path, sep="\t", index=False, lineterminator="\n")
        manifest_for(path, route="a", case_id="case1")
        return path

    def _family_copy(self, root: Path, count: int, mismatched: int) -> Path:
        out = root / "z02"
        out.mkdir(parents=True, exist_ok=True)
        discrimination = out / "family_copy_discrimination.tsv"
        relationships = out / "family_copy_relationships.tsv"
        disc_rows = []
        rel_rows = []
        for index in range(1, count + 1):
            gene = f"G{index:03d}"
            mismatched_row = index <= mismatched
            disc_rows.append({
                "family_id": "FAM1", "copy_set_id": "CS1", "copy_set_status": "validated",
                "candidate_gene_id": gene,
                "subfamily": "OTHER_SUBFAMILY" if mismatched_row else "TARGET_SUBFAMILY",
                "orthology_class": "", "orthology_confidence": "",
                "synteny_status": "", "domain_architecture_match": "",
                "relevant_expression_status": "", "position_discrimination": "",
                "paralog_ambiguity": "", "discriminating_axes": "", "leading_axes": "",
                "copy_specific_support": "", "copy_specific_counterevidence": "",
                "family_level_support": "", "disambiguation_status": "",
                "rank_contribution": "none", "review_reason": "",
                "target_id": "TGT1", "target_verification_status": "verified",
                "target_reference_gene": "REF1", "target_subfamily": "TARGET_SUBFAMILY",
                "target_domain_architecture": "",
                "subfamily_vs_target": "target_mismatch" if mismatched_row else "target_match",
                "domain_vs_target": "", "domain_vs_family_modal": "",
                "subfamily_family_modal": "", "target_specific_axes": "subfamily_match",
                "selection_basis": "",
            })
            rel_rows.append({
                "family_id": "FAM1", "candidate_gene_id": gene, "reference_gene": "REF1",
                "orthology_class": "", "orthology_providers": "", "evidence_chain": "",
                "subfamily": "OTHER_SUBFAMILY" if mismatched_row else "TARGET_SUBFAMILY",
                "subfamily_source": "curated_map", "paralog_count": "", "notes": "",
            })
        pd.DataFrame(disc_rows).to_csv(discrimination, sep="\t", index=False, lineterminator="\n")
        pd.DataFrame(rel_rows).to_csv(relationships, sep="\t", index=False, lineterminator="\n")
        return discrimination

    def test_curated_cross_subfamily_verdicts_narrow_below_budget_and_re_freeze(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            family_copy = self._family_copy(root, count=30, mismatched=12)
            result = wet_candidate_handoff(
                str(root / "out"),
                route_a_bridge_arg=str(self._bridge(root, 30)),
                family_copy_arg=str(family_copy),
                deterministic_subfamily_refinement=True,
                wet_budget=20,
                case_id="case1",
            )
            self.assertEqual("wet_reviewable", result["handoff_state"])
            self.assertEqual("within_wet_budget", result["budget_status"])
            # Only the sourced cross-subfamily members left the set.
            self.assertEqual(18, result["candidate_count"])
            self.assertFalse(result["heavy_stage_invoked"])
            self.assertIsNone(result["adaptive_synthesis_request"])
            wet = pd.read_csv(result["wet_candidates"], sep="\t", dtype=str).fillna("")
            self.assertEqual(18, len(wet))
            self.assertTrue((wet["support_tier"] == "tier1_verified_support").all())
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            summary = manifest["deterministic_subfamily_refinement"]
            self.assertEqual(30, summary["input_candidate_count"])
            self.assertEqual(18, summary["output_candidate_count"])
            self.assertEqual(12, summary["removed_candidate_count"])
            self.assertEqual("within_wet_budget", summary["regate_status"])
            self.assertEqual("verified_lineage", summary["lineage"]["verification_state"])

    def test_insufficient_falsifications_still_stop_as_review_required(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            family_copy = self._family_copy(root, count=30, mismatched=5)
            result = wet_candidate_handoff(
                str(root / "out"),
                budget_policy="strict",
                route_a_bridge_arg=str(self._bridge(root, 30)),
                family_copy_arg=str(family_copy),
                deterministic_subfamily_refinement=True,
                wet_budget=20,
                case_id="case1",
            )
            self.assertEqual("review_required_over_budget", result["handoff_state"])
            self.assertEqual("not_emitted_over_budget_review_required", result["membership_state"])
            self.assertIsNone(result["wet_candidates"])
            self.assertEqual(25, result["candidate_count"])
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual("review_required_over_budget", manifest["deterministic_subfamily_refinement"]["regate_status"])

    def test_flag_without_family_copy_falls_back_to_the_request(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = wet_candidate_handoff(
                str(root / "out"),
                budget_policy="strict",
                route_a_bridge_arg=str(self._bridge(root, 30)),
                deterministic_subfamily_refinement=True,
                wet_budget=20,
                case_id="case1",
            )
            self.assertEqual("adaptive_refinement_requested", result["handoff_state"])
            self.assertTrue(result["heavy_stage_invoked"])
            self.assertIn("no_family_copy_input", result["heavy_stage_invocation_reason"])

    def test_flag_off_keeps_the_existing_request_only_behavior(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            family_copy = self._family_copy(root, count=30, mismatched=12)
            result = wet_candidate_handoff(
                str(root / "out"),
                budget_policy="strict",
                route_a_bridge_arg=str(self._bridge(root, 30)),
                family_copy_arg=str(family_copy),
                wet_budget=20,
                case_id="case1",
            )
            self.assertEqual("adaptive_refinement_requested", result["handoff_state"])
            self.assertTrue(result["heavy_stage_invoked"])
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertIsNone(manifest["deterministic_subfamily_refinement"])


if __name__ == "__main__":
    unittest.main()
