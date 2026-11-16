"""Versioned scope changes reuse recorded Tier facts and the existing consumer."""
import json
import io
from contextlib import redirect_stdout
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

from candidate_qtl.boundary_tier_review import (
    UPPER_EVIDENCE_SCOPED_POLICY as V3, EVIDENCE_SCOPED_POLICY as V2,
    partition_relevance_tiers, freeze_boundary_scope, _load_scope,
    build_boundary_review_request, consume_boundary_selection,
)
from candidate_qtl.relevance_tiers import NEW_TIERS, T1, T2, T3
from tests.test_boundary_tier_review import sha, write
from tests import test_boundary_tier_review as boundary_fixtures
from tests.test_wet_relevance_selection import handoff
from candidate_qtl.claim_review_session import _files
from scripts.agentic_boundary import prepare, merge
from scripts.claim_review_session import main as session_cli


def roster(*counts):
    return [{"candidate_gene_id": f"t{tier}_{i}", "relevance_tier": NEW_TIERS[tier],
             "support_state": "review", "membership_basis": "route_review_candidate"}
            for tier, count in enumerate(counts) for i in range(count)]


class UpperScopeTests(unittest.TestCase):
    def plan(self, rows, **kwargs):
        return partition_relevance_tiers(rows, policy_version=V3,
            eligible_gene_ids=[row["candidate_gene_id"] for row in rows], **kwargs)

    def test_small_upper_pool_preserves_lower_eligibility_without_padding(self):
        rows = roster(1, 2, 3, 100)
        plan = self.plan(rows)
        self.assertEqual((6, 0, 100), tuple(len(plan[key]) for key in ("frozen", "boundary", "lower")))
        self.assertEqual(106, len(plan["eligible_gene_ids"]))
        self.assertEqual(6, len(plan["active_gene_ids"]))
        self.assertTrue(plan["upper_priority_applied"])

    def test_upper_overflow_compares_only_upper_and_zero_upper_uses_old_scope(self):
        for counts, expected in [((0, 30, 0, 100), (0, 30, 100)), ((0, 0, 0, 30), (0, 30, 0))]:
            with self.subTest(counts=counts):
                rows = roster(*counts)
                plan = self.plan(rows)
                self.assertEqual(expected, tuple(len(plan[key]) for key in ("frozen", "boundary", "lower")))
                if counts[1]:
                    self.assertEqual({NEW_TIERS[1]}, {row["relevance_tier"] for row in plan["boundary"]})

    def test_legacy_upper_meaning_t1_copies_and_protected_c_are_preserved(self):
        rows = roster(0, 0, 0, 3)
        rows += [{"candidate_gene_id": tier, "relevance_tier": tier} for tier in (T1, T2, T3)]
        plan = self.plan(rows, protected_gene_ids=["t3_0"])
        self.assertEqual({T1, T2, "t3_0"}, {row["candidate_gene_id"] for row in plan["frozen"]})
        self.assertIn(T3, {row["candidate_gene_id"] for row in plan["lower"]})
        copies = self.plan(roster(30, 0, 0, 100))
        self.assertEqual(30, len(copies["frozen"]))
        self.assertFalse(copies["boundary"])

    def test_unproposed_upper_does_not_suppress_eligible_lower(self):
        rows = roster(0, 1, 0, 30)
        plan = partition_relevance_tiers(rows, policy_version=V3,
            eligible_gene_ids=[row["candidate_gene_id"] for row in rows if row["relevance_tier"] == NEW_TIERS[3]])
        self.assertFalse(plan["upper_priority_applied"])
        self.assertEqual(30, len(plan["boundary"]))

    def test_protected_c_upper_does_not_activate_priority_or_displace_ordinary_candidates(self):
        for tier in (T1, T2, *NEW_TIERS[:3]):
            with self.subTest(tier=tier):
                rows = roster(0, 0, 0, 20) + [{"candidate_gene_id": "protected", "relevance_tier": tier}]
                plan = self.plan(rows, protected_gene_ids=["protected"])
                self.assertFalse(plan["upper_priority_applied"])
                self.assertEqual(21, len(plan["frozen"]))
                self.assertFalse(plan["boundary"])
                self.assertFalse(plan["lower"])

    def test_first_overflow_tier_excludes_later_reported_upper_from_comparison(self):
        rows = roster(1, 30, 10, 100)
        plan = self.plan(rows)
        self.assertEqual((1, 30, 110), tuple(len(plan[key]) for key in ("frozen", "boundary", "lower")))
        self.assertEqual({NEW_TIERS[1]}, {r["relevance_tier"] for r in plan["boundary"]})
        self.assertEqual((141, 31), (len(plan["eligible_gene_ids"]), len(plan["active_gene_ids"])))

    def test_frozen_v2_keeps_original_lower_scope(self):
        rows = roster(1, 2, 3, 100)
        plan = partition_relevance_tiers(rows, policy_version=V2, eligible_gene_ids=[r["candidate_gene_id"] for r in rows])
        self.assertEqual((6, 100, 0), tuple(len(plan[key]) for key in ("frozen", "boundary", "lower")))

    def test_scope_replay_rejects_changed_active_pool_even_with_new_manifest_hash(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows = roster(1, 2, 3, 100)
            source = root / "inventory.tsv"
            write(source, rows)
            scope = freeze_boundary_scope(str(source), str(root / "scope"), inventory_sha256=sha(source),
                case_id="test", policy_version=V3, eligible_gene_ids=[r["candidate_gene_id"] for r in rows])
            path = Path(scope["manifest"])
            payload = json.loads(path.read_text())
            payload["active_gene_ids"].append("t3_0")
            path.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError, "active scope differs"):
                _load_scope(str(path), sha(path), "test")

    def test_finish_or_partial_never_restores_lower_or_starts_another_round(self):
        for keep, completed in ((0, True), (30, True), (30, False)):
            with self.subTest(keep=keep, completed=completed), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                source = root / "inventory.tsv"
                rows = roster(0, 30, 0, 100)
                write(source, rows)
                scope = freeze_boundary_scope(str(source), str(root / "scope"), inventory_sha256=sha(source),
                    case_id="test", policy_version=V3, eligible_gene_ids=[r["candidate_gene_id"] for r in rows])
                _load_scope(scope["manifest"], scope["manifest_sha256"], "test")
                request = build_boundary_review_request(scope["manifest"], str(root / "request"),
                    scope_sha256=scope["manifest_sha256"], case_id="test")
                packet = json.loads(Path(request["request"]).read_text())
                self.assertEqual(30, len(packet["candidate_gene_ids"]))
                self.assertTrue(all(gene.startswith("t1_") for gene in packet["candidate_gene_ids"]))
                self.assertIn("do not add a Tier/content reviewer", " ".join(packet["policy"]))
                self.assertEqual("original_assessment_v1", packet["assessment_mode"])
                self.assertFalse(packet["reviewer_schema"]["independent"])
                helper = boundary_fixtures.ScopeTests()
                helper.root = root
                proposal, independent = helper.review(keep=keep, scope=scope)
                if not completed:
                    partial = pd.read_csv(proposal, sep="\t").fillna("")
                    partial["decision"] = "not_evaluated"
                    partial["reason"] = "Assessment unfinished"
                    partial.to_csv(proposal, sep="\t", index=False)
                with self.assertRaisesRegex(ValueError, "legacy provisional adoption"):
                    consume_boundary_selection(scope["manifest"], proposal, str(root / "invalid-provisional"),
                        scope_sha256=scope["manifest_sha256"], case_id="test", reviewer_id="reviewer", reviewer_model="model-b",
                        executor_id="executor", executor_model="model-a",
                        allow_provisional_wet=True)
                result = consume_boundary_selection(scope["manifest"], proposal, str(root / "selected"),
                    scope_sha256=scope["manifest_sha256"], case_id="test",
                    executor_id="executor", executor_model="model-a")
                self.assertEqual("original_assessment_v1", result["assessment_mode"])
                self.assertEqual("incomplete" if not completed else "complete", result["review_history_state"])
                with self.assertRaisesRegex(ValueError, "must not be represented as independent review"):
                    consume_boundary_selection(scope["manifest"], proposal, str(root / "fake-review"),
                        scope_sha256=scope["manifest_sha256"], case_id="test",
                        executor_id="executor", executor_model="model-a", reviewer_id="reviewer",
                        reviewer_model="model-b", independent_review_arg=independent)
                wet = pd.read_csv(result["wet_candidates"]["path"], sep="\t")
                self.assertEqual(keep if completed else 30, len(wet))
                self.assertFalse(result["next_scope_manifest"])
                self.assertEqual(not completed, bool(result["resume_scope_manifest"]))
                self.assertEqual(130, len(result["eligible_gene_ids"]))
                self.assertEqual(30, len(result["active_gene_ids"]))

    def test_real_handoff_defaults_v3_trace_keeps_lower_nonnegative(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = handoff(root, (1, 2, 35, 100, 0), policy=None)
            self.assertEqual(V3, result["first_stage_selection_policy_version"])
            self.assertEqual(3, len(pd.read_csv(result["wet_candidates"], sep="\t")))
            self.assertFalse(result["adaptive_synthesis_request"])
            manifest = json.loads(Path(result["manifest"]).read_text())
            self.assertEqual(3, manifest["candidate_count_populations"]["active_first_stage"])
            self.assertEqual(38, manifest["candidate_count_populations"]["eligible_first_stage"])
            trace = pd.read_csv(root / "wet/wet_candidate_trace.tsv", sep="\t")
            lower = trace[trace.relevance_tier == T3]
            self.assertEqual({"not_selected_upper_tier_priority"}, set(lower.selection_disposition))
            self.assertEqual({"not_established"}, set(lower.biological_exclusion))
            scope = json.loads(Path(result["boundary_scope_manifest"]).read_text())
            self.assertEqual((38, 3), (len(scope["eligible_gene_ids"]), len(scope["active_gene_ids"])))
            request = build_boundary_review_request(scope["manifest"] if "manifest" in scope else result["boundary_scope_manifest"],
                str(root / "no-request"), scope_sha256=result["boundary_scope_sha256"], case_id="case1")
            self.assertNotIn("request", request)

    def test_original_assessment_transport_and_finish_keep_prefix_without_lower_revival(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = root / "first"
            first.mkdir()
            source = first / "inventory.tsv"
            rows = roster(2, 0, 30, 100)
            write(source, rows)
            scope = freeze_boundary_scope(str(source), str(first / "scope"), inventory_sha256=sha(source),
                case_id="test", policy_version=V3, eligible_gene_ids=[r["candidate_gene_id"] for r in rows])
            request = build_boundary_review_request(scope["manifest"], str(first / "request"),
                scope_sha256=scope["manifest_sha256"], case_id="test")
            prepared = root / "prepared"
            plan = prepare(request["request"], str(prepared), max_candidates=40)
            self.assertEqual("original_assessment_v1", plan["assessment_mode"])
            shard = plan["shards"][0]
            response = root / "response.json"
            response.write_text(json.dumps({"case_id": "test", "source_wet_sha256": plan["source_wet_sha256"],
                "shard_id": shard["shard_id"], "executor_id": "original", "executor_model": "model",
                "responses": [{"candidate_gene_id": gene, "decision": "defer", "reason": "Existing upper evidence sufficiently prioritizes the preserved prefix"}
                              for gene in shard["candidate_gene_ids"]]}))
            proposal = root / "proposal.tsv"
            merged = merge(str(prepared), [str(response)], str(proposal), role="proposal",
                           executor_id="original", executor_model="model")
            self.assertEqual("original_assessment_v1", merged["assessment_mode"])
            with self.assertRaisesRegex(ValueError, "without independent review"):
                merge(str(prepared), [], str(root / "review.tsv"), role="review")
            completion = first / "completion.json"
            completion.write_text(json.dumps({"case_id": "test", "files": _files(first), "wet": {
                "boundary_scope_manifest": scope["manifest"], "boundary_scope_sha256": scope["manifest_sha256"],
                "adaptive_synthesis_request": request["request"], "adaptive_round": 0, "max_adaptive_rounds": 1,
                "first_stage_selection_policy_version": V3}}))
            frozen = source.read_bytes()
            printed = io.StringIO()
            with patch("sys.argv", ["claim_review_session", "finish-wet-review", "--handoff", str(first),
                       "--completion-sha256", sha(completion), "--out", str(root / "finished"),
                       "--proposals", str(proposal), "--executor-id", "original", "--executor-model", "model"]), redirect_stdout(printed):
                session_cli()
            result = json.loads(printed.getvalue())
            self.assertEqual("complete", result["review_history_state"])
            self.assertEqual("original_assessment_v1", result["review_history_basis"])
            selection = json.loads(Path(result["boundary_selection"]["boundary_selection_manifest"]["path"]).read_text())
            self.assertEqual("declared_identity_model_consistency", selection["executor_identity_basis"])
            self.assertEqual("not_required_original_assessment", selection["independent_review_state"])
            self.assertEqual("", selection["reviewer_id"])
            self.assertIsNone(selection["independent_review_artifact"])
            self.assertEqual(2, result["adopted_candidate_count"])
            self.assertFalse(result["next_boundary_scope_manifest"])
            self.assertEqual(frozen, source.read_bytes())

    def test_original_assessment_keeps_hash_actor_roster_and_reason_checks(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "inventory.tsv"
            rows = roster(0, 30)
            write(source, rows)
            scope = freeze_boundary_scope(str(source), str(root / "scope"), inventory_sha256=sha(source),
                case_id="test", policy_version=V3, eligible_gene_ids=[r["candidate_gene_id"] for r in rows])
            helper = boundary_fixtures.ScopeTests()
            helper.root = root
            proposal, _ = helper.review(keep=0, scope=scope)
            original = pd.read_csv(proposal, sep="\t").fillna("")
            for column, value, error in (("source_wet_sha256", "wrong", "different frozen Wet hash"),
                                         ("executor_id", "another", "does not match"),
                                         ("executor_model", "another", "does not match"),
                                         ("candidate_gene_id", "foreign", "higher/lower-tier candidate"),
                                         ("reason", "", "non-empty reason")):
                with self.subTest(column=column):
                    changed = original.copy()
                    changed.loc[0, column] = value
                    changed.to_csv(proposal, sep="\t", index=False)
                    with self.assertRaisesRegex(ValueError, error):
                        consume_boundary_selection(scope["manifest"], proposal, str(root / column),
                            scope_sha256=scope["manifest_sha256"], case_id="test",
                            executor_id="executor", executor_model="model-a")
