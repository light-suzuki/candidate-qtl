"""Whole-roster adoption through the production handoff; not model accuracy."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from contextlib import nullcontext

import pandas as pd

from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from candidate_qtl.relevance_tiers import T1, T2, T3, T4, UNRESOLVED
from tests._issue664_helpers import manifest_for
from tests.test_relevance_tiers import known, bridge


def handoff(root, counts, *, reverse=False, max_rounds=1, provisional=False, policy="evidence_scoped_v2"):
    root.mkdir(parents=True, exist_ok=True)
    atlas_rows, reviews = [], []
    for tier, count in zip((T1, T2, T3, T4, UNRESOLVED), counts):
        for index in range(count):
            gene = f"{tier[:2]}_{index:03d}"
            atlas_rows.append({"gene_id": gene, "family_id": "F" if tier != UNRESOLVED else ""})
            if tier in (T1, T2, T3):
                record = known() if tier == T1 else bridge()
                record.update(candidate_gene_id=gene, route="b", claim_id=gene,
                              support_status="unverified_route_b_review", rank_contribution="none")
                if tier == T3:
                    record["claim_type"] = "mechanistic_review_seed"
                reviews.append(record)
    atlas = root / "atlas.tsv"
    pd.DataFrame(atlas_rows[::-1] if reverse else atlas_rows).to_csv(atlas, sep="\t", index=False)
    supported = root / "supported.tsv"
    pd.DataFrame(columns=["candidate_gene_id", "route", "support_status"]).to_csv(supported, sep="\t", index=False)
    args = {}
    if provisional:
        screen = root / "screen.tsv"
        pd.DataFrame([{"case_id": "case1", "unit_id": f"unit_{index}",
                       "member_interval_genes": row["gene_id"],
                       "screen_state": "screened", "answer_sha256": "a" * 64}
                      for index, row in enumerate(atlas_rows)]).to_csv(screen, sep="\t", index=False)
        args.update(provisional_screen_ledger_arg=str(screen),
                    provisional_screen_ledger_sha256=hashlib.sha256(screen.read_bytes()).hexdigest())
    if reviews:
        review = root / "review.tsv"
        pd.DataFrame(reviews[::-1] if reverse else reviews).fillna("").to_csv(review, sep="\t", index=False)
        meta = root / "review_manifest.json"
        meta.write_text(json.dumps({"contract": "route-review-candidates-1.0.0",
            "stage": "route_b_review_candidate_freeze", "case_id": "case1", "inputs": {},
            "output": str(review), "artifact_sha256": hashlib.sha256(review.read_bytes()).hexdigest()}))
        args.update(review_candidates_arg=str(review), review_candidates_manifest_arg=str(meta))
    # Historical fixtures explicitly exercise their frozen v2 contract.
    with patch("candidate_qtl.boundary_tier_review.DEFAULT_SELECTION_POLICY", policy) if policy else nullcontext():
        return wet_candidate_handoff(str(root / "wet"), route_b_candidates_arg=str(supported),
            route_b_candidates_manifest_arg=str(manifest_for(supported, route="b")),
            interval_atlas_arg=str(atlas), case_id="case1", max_adaptive_rounds=max_rounds, **args)


class WholeRosterHandoffTests(unittest.TestCase):
    def test_opt_in_provisional_adopts_concrete_and_preserves_full_interval(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = handoff(Path(tmp), (0, 0, 3, 40, 0), provisional=True)
            wet = pd.read_csv(result["wet_candidates"], sep="\t")
            self.assertEqual(set(wet.candidate_gene_id), {"T3_000", "T3_001", "T3_002"})
            self.assertEqual(result["handoff_state"], "provisional_adopted_wet")
            self.assertIsNone(result["boundary_scope_manifest"])
            manifest = json.loads(Path(result["wet_candidates_manifest"]).read_text())
            self.assertTrue(manifest["provisional_wet"]["eligible"])
            unresolved = pd.read_csv(manifest["interval_unresolved"], sep="\t")
            self.assertEqual(len(unresolved), 43)
            self.assertEqual((unresolved.selection_disposition == "not_evaluated_identity_only").sum(), 40)

    def test_opt_in_retains_whole_tier_when_boundary_has_comparison_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = handoff(Path(tmp), (0, 0, 40, 5, 0), provisional=True)
            wet = pd.read_csv(result["wet_candidates"], sep="\t")
            self.assertEqual(len(wet), 40)
            self.assertIsNotNone(result["boundary_scope_manifest"])
            manifest = json.loads(Path(result["wet_candidates_manifest"]).read_text())
            self.assertFalse(manifest["provisional_wet"]["eligible"])

    def test_case_identity_is_required_before_creating_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "wet"
            with self.assertRaisesRegex(ValueError, "non-empty case_id"):
                wet_candidate_handoff(str(output), interval_atlas_arg=str(Path(tmp) / "atlas.tsv"))
            self.assertFalse(output.exists())

    def test_copy_context_is_preserved_without_demoting_relevance(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            handoff(root, (0, 2, 0, 0, 0))
            copies = root / "copies.tsv"
            pd.DataFrame([
                {"candidate_gene_id": "T2_000", "family_id": "F1", "disambiguation_status": "paralog_ambiguous"},
                {"candidate_gene_id": "T2_000", "family_id": "F2", "disambiguation_status": "co_preferred_tie"},
            ]).to_csv(copies, sep="\t", index=False)
            result = wet_candidate_handoff(str(root / "with-copy"),
                route_b_candidates_arg=str(root / "supported.tsv"),
                route_b_candidates_manifest_arg=str(root / "supported_manifest.json"),
                interval_atlas_arg=str(root / "atlas.tsv"), family_copy_arg=str(copies),
                review_candidates_arg=str(root / "review.tsv"),
                review_candidates_manifest_arg=str(root / "review_manifest.json"), case_id="case1")
            wet = pd.read_csv(result["wet_candidates"], sep="\t").fillna("")
            row = wet.set_index("candidate_gene_id").loc["T2_000"]
            self.assertEqual(row.relevance_tier, T2)
            self.assertEqual(set(row.copy_state.split(";")), {"paralog_ambiguous", "co_preferred_tie"})
            self.assertEqual(len(json.loads(row.copy_evidence_json)), 2)
            self.assertEqual(len(wet), 2)

    def test_boundary_contains_entire_tier_and_keeps_higher(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = handoff(Path(tmp), (3, 8, 35, 20, 2))
            wet = pd.read_csv(result["wet_candidates"], sep="\t")
            self.assertEqual(len(wet), 46)
            scope = json.loads(Path(result["boundary_scope_manifest"]).read_text())
            self.assertEqual(scope["counts"], {"frozen": 11, "boundary": 35, "lower": 22})
            boundary = pd.read_csv(scope["artifacts"]["boundary"]["path"], sep="\t")
            self.assertEqual(set(boundary.relevance_tier), {T3})
            self.assertEqual(scope["policy_version"], "evidence_scoped_v2")
            self.assertEqual(result["second_stage_state"], "requested")
            trace = pd.read_csv(Path(tmp) / "wet/wet_candidate_trace.tsv", sep="\t")
            self.assertEqual(len(trace), 68)
            self.assertEqual((trace.selection_disposition == "not_selected_no_upstream_candidate").sum(), 22)
            self.assertTrue(trace.upstream_llm_assessment.eq("not_evaluated_by_handoff").all())

    def test_twenty_two_stays_twenty_two(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = handoff(Path(tmp), (3, 10, 9, 30, 0))
            wet = pd.read_csv(result["wet_candidates"], sep="\t")
            self.assertEqual(len(wet), 22)
            self.assertFalse(result["adaptive_synthesis_request"])

    def test_known_rescue_stays_while_later_eligible_candidates_reach_review(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = handoff(Path(tmp), (30, 2, 2, 0, 0))
            wet = pd.read_csv(result["wet_candidates"], sep="\t")
            self.assertEqual(len(wet), 34)
            self.assertEqual(int(wet.relevance_tier.eq(T1).sum()), 30)
            scope = json.loads(Path(result["boundary_scope_manifest"]).read_text())
            boundary = pd.read_csv(scope["artifacts"]["boundary"]["path"], sep="\t")
            self.assertEqual(len(boundary), 4)
            self.assertEqual(result["second_stage_state"], "requested")

    def test_unresolved_is_rescue_not_empty_or_negative(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = handoff(Path(tmp), (0, 0, 0, 0, 40))
            self.assertIsNone(result["wet_candidates"])
            self.assertEqual(result["handoff_state"], "review_required_no_supported_candidates")

    def test_evidence_scoped_overflow_carries_lower_tier_seeds_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            # T2 overflows; the already-proposed T3 seed remains comparison
            # material, while hundreds of identity-only rows stay unselected.
            result = handoff(Path(tmp), (0, 26, 1, 300, 300))
            wet = pd.read_csv(result["wet_candidates"], sep="\t")
            scope = json.loads(Path(result["boundary_scope_manifest"]).read_text())
            boundary = pd.read_csv(scope["artifacts"]["boundary"]["path"], sep="\t")
            self.assertEqual(len(wet), 27)
            self.assertEqual(len(boundary), 27)
            self.assertEqual({T2, T3}, set(boundary.relevance_tier))
            self.assertEqual(scope["policy_version"], "evidence_scoped_v2")
            self.assertEqual(scope["counts"]["lower"], 600)

    def test_evidence_scoped_19_20_21_25_26_with_same_t3_seed_and_large_identity_roster(self):
        for count in (19, 20, 21, 25, 26):
            with self.subTest(eligible_count=count), tempfile.TemporaryDirectory() as tmp:
                result = handoff(Path(tmp), (0, count, 1, 250, 250))
                scope = json.loads(Path(result["boundary_scope_manifest"]).read_text())
                if count < 25:
                    wet = pd.read_csv(result["wet_candidates"], sep="\t")
                    self.assertEqual(len(wet), count + 1)
                    self.assertEqual(scope["counts"]["boundary"], 0)
                elif count == 25:
                    wet = pd.read_csv(result["wet_candidates"], sep="\t")
                    boundary = pd.read_csv(scope["artifacts"]["boundary"]["path"], sep="\t")
                    self.assertEqual(26, len(wet))
                    self.assertEqual(25, scope["counts"]["frozen"])
                    self.assertEqual([T3], list(boundary.relevance_tier.drop_duplicates()))
                    self.assertEqual(1, len(boundary))
                else:
                    boundary = pd.read_csv(scope["artifacts"]["boundary"]["path"], sep="\t")
                    self.assertEqual(len(boundary), count + 1)
                    self.assertEqual({T2, T3}, set(boundary.relevance_tier))
                    self.assertEqual(len(pd.read_csv(result["wet_candidates"], sep="\t")), count + 1)
        with tempfile.TemporaryDirectory() as tmp:
            result = handoff(Path(tmp), (0, 8, 0, 250, 250))
            self.assertEqual(8, len(pd.read_csv(result["wet_candidates"], sep="\t")))
        with tempfile.TemporaryDirectory() as tmp:
            result = handoff(Path(tmp), (0, 8, 4, 250, 250))
            wet = pd.read_csv(result["wet_candidates"], sep="\t")
            self.assertEqual(len(wet), 12)
            scope = json.loads(Path(result["boundary_scope_manifest"]).read_text())
            self.assertEqual(scope["counts"]["boundary"], 0)

    def test_no_id_order_truncation(self):
        with tempfile.TemporaryDirectory() as tmp:
            first = handoff(Path(tmp) / "a", (0, 8, 35, 0, 0))
            second = handoff(Path(tmp) / "b", (0, 8, 35, 0, 0), reverse=True)
            self.assertEqual(Path(first["wet_candidates"]).read_bytes(), Path(second["wet_candidates"]).read_bytes())

    def test_round_budget_keeps_pending_candidates(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = handoff(Path(tmp), (0, 8, 35, 0, 0), max_rounds=0)
            self.assertEqual(len(pd.read_csv(result["wet_candidates"], sep="\t")), 43)
            self.assertFalse(result["adaptive_synthesis_request"])
            self.assertEqual(result["second_stage_state"], "review_required_round_budget_exhausted")
