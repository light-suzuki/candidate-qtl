from __future__ import annotations

import fnmatch
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from candidate_qtl.budget_account import BudgetAccount
from candidate_qtl.llm_portfolio import (
    FileReviewBackend,
    audit_proposer_diversity,
    build_critic,
    build_review_batches,
    reduce_portfolio,
    run_cluster_review,
)
from candidate_qtl.reachability import RANK_IMPACTING_PATTERNS, _matches_rank_impacting

CASE_ID = "portfolio-case"


class RecordingBackend:
    """In-memory review backend for tests; mirrors the FakeBackend house pattern
    from tests/test_tier2_planner.py -- no network path, records every call."""

    def __init__(self, responses: dict[str, dict[str, object]] | None = None):
        self.provider = "test-provider"
        self.model = "test-model"
        self.version = "1.0"
        self.cache = None
        self.available = True
        self._responses = responses or {}
        self.calls: list[str] = []

    def review(self, batch: dict[str, object]) -> dict[str, object]:
        self.calls.append(str(batch["batch_id"]))
        entry = self._responses.get(str(batch["batch_id"]))
        if entry is None:
            return {"items": [], "input_tokens": 1, "output_tokens": 1}
        return entry


def _handoff_row(gene_id: str, cluster_id: str, allocation: str, rescue_bucket: str, review_state: str) -> dict[str, object]:
    return {
        "gene_id": gene_id, "cluster_id": cluster_id, "allocation": allocation,
        "rescue_bucket": rescue_bucket, "plan_step": "1", "query_role": "representative",
        "query_status": "executed" if review_state == "priority_review" else "not_queried",
        "expansion_condition": "representative_first", "search_result_status": "hit",
        "review_state": review_state, "handoff_note": "searched",
    }


def _plan_row(gene_id: str, cluster_id: str, ambiguity_status: str = "not_ambiguous", rescue_escalated: str = "no") -> dict[str, object]:
    return {"gene_id": gene_id, "cluster_id": cluster_id, "ambiguity_status": ambiguity_status, "rescue_escalated": rescue_escalated}


def _cluster_row(cluster_id: str, member_gene_ids: str, allocation: str, rescue_bucket: str = "", dependency_cluster_ids: str = "") -> dict[str, object]:
    return {
        "cluster_id": cluster_id, "member_gene_ids": member_gene_ids, "allocation": allocation,
        "rescue_bucket": rescue_bucket, "dependency_cluster_ids": dependency_cluster_ids,
    }


def _basic_fixture():
    """cls_A(priority: G1 rep, G2 ambiguous-homology expansion), cls_B(defer: G3, not eligible),
    cls_C(rescue/unclassified, escalated: G4), cls_D(defer/weak_annotation, not eligible: G5)."""
    handoff_rows = [
        _handoff_row("G1", "cls_A", "priority", "", "priority_review"),
        _handoff_row("G2", "cls_A", "priority", "", "deferred"),
        _handoff_row("G3", "cls_B", "defer", "", "deferred"),
        _handoff_row("G4", "cls_C", "rescue", "unclassified", "rescue_pending_review"),
        _handoff_row("G5", "cls_D", "defer", "weak_annotation", "not_evaluated"),
    ]
    plan_rows = [
        _plan_row("G1", "cls_A"),
        _plan_row("G2", "cls_A", ambiguity_status="ambiguous_homology"),
        _plan_row("G3", "cls_B"),
        _plan_row("G4", "cls_C", rescue_escalated="yes"),
        _plan_row("G5", "cls_D"),
    ]
    cluster_rows = [
        _cluster_row("cls_A", "G1;G2", "priority"),
        _cluster_row("cls_B", "G3", "defer"),
        _cluster_row("cls_C", "G4", "rescue", rescue_bucket="unclassified"),
        _cluster_row("cls_D", "G5", "defer", rescue_bucket="weak_annotation"),
    ]
    return handoff_rows, plan_rows, cluster_rows


class BuildReviewBatchesTests(unittest.TestCase):
    def test_eligibility_and_single_batch_covers_two_clusters(self) -> None:
        handoff_rows, plan_rows, cluster_rows = _basic_fixture()
        batches = build_review_batches(handoff_rows, plan_rows, cluster_rows, batch_size=5)
        self.assertEqual(1, len(batches))
        gene_ids = sorted(member["gene_id"] for member in batches[0]["members"])
        # G3 (deferred, not ambiguous, not rescue-escalated) and G5 (not_evaluated,
        # not ambiguous) are never eligible.
        self.assertEqual(["G1", "G2", "G4"], gene_ids)

    def test_never_one_call_per_gene_batching(self) -> None:
        """Case 12: one batch (one backend call) covers several genes."""
        handoff_rows, plan_rows, cluster_rows = _basic_fixture()
        batches = build_review_batches(handoff_rows, plan_rows, cluster_rows, batch_size=5)
        self.assertEqual(1, len(batches))
        self.assertGreater(len(batches[0]["members"]), 1)

    def test_batch_id_is_deterministic_content_hash(self) -> None:
        handoff_rows, plan_rows, cluster_rows = _basic_fixture()
        first = build_review_batches(handoff_rows, plan_rows, cluster_rows)
        second = build_review_batches(handoff_rows, plan_rows, cluster_rows)
        self.assertEqual([b["batch_id"] for b in first], [b["batch_id"] for b in second])
        self.assertTrue(all(b["batch_id"].startswith("batch_") for b in first))

    def test_batch_size_packs_whole_clusters_not_genes(self) -> None:
        handoff_rows, plan_rows, cluster_rows = _basic_fixture()
        batches = build_review_batches(handoff_rows, plan_rows, cluster_rows, batch_size=1)
        self.assertEqual(2, len(batches))
        cluster_sets = [set(b["cluster_ids"]) for b in batches]
        self.assertIn({"cls_A"}, cluster_sets)
        self.assertIn({"cls_C"}, cluster_sets)


class RunClusterReviewTests(unittest.TestCase):
    def test_supported_claim_round_trip(self) -> None:
        """Case 1."""
        handoff_rows, plan_rows, cluster_rows = _basic_fixture()
        batches = build_review_batches(handoff_rows, plan_rows, cluster_rows)
        batch_id = batches[0]["batch_id"]
        backend = RecordingBackend({
            batch_id: {
                "items": [{"gene_id": "G1", "supported_mechanism_ids": ["mech_001"], "rationale_source_ids": ["SRC1"]}],
                "input_tokens": 10, "output_tokens": 5,
            }
        })
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            account = BudgetAccount()
            account.declare({"llm_call_limit": 10, "token_limit": 1000, "sequence_query_limit": 0, "wall_time_limit_seconds": None, "peak_memory_limit_mb": None})
            manifest = run_cluster_review(batches, account, out, CASE_ID, backend, mechanisms={"mech_001"}, allowed_source_ids={"SRC1"})
            self.assertEqual(1, len(backend.calls))
            frame = pd.read_csv(out / "llm_review.tsv", sep="\t", dtype=str).fillna("")
            row = frame[frame["gene_id"] == "G1"].iloc[0]
            self.assertEqual("reviewed", row["review_status"])
            self.assertEqual("mech_001", row["supported_mechanism_ids"])
            self.assertEqual("", row["contradicted_mechanism_ids"])
            self.assertEqual("none", row["rank_contribution"])
            # counts are per review ROW: the batch covers three eligible genes
            # and the reviewer answered for one of them, so the other two are
            # honest empty_response rows rather than fabricated reviews.
            self.assertEqual(1, manifest["counts"]["reviewed"])
            self.assertEqual(2, manifest["counts"]["empty"])
            self.assertEqual(1, manifest["batch_counts"]["called"])
            frozen = account.freeze()
            self.assertEqual(1, frozen["usage"]["llm_call_count"])
            self.assertEqual(10, frozen["usage"]["llm_input_tokens"])
            self.assertEqual(5, frozen["usage"]["llm_output_tokens"])

    def test_contradicted_claim_never_becomes_support(self) -> None:
        """Case 2."""
        handoff_rows, plan_rows, cluster_rows = _basic_fixture()
        batches = build_review_batches(handoff_rows, plan_rows, cluster_rows)
        batch_id = batches[0]["batch_id"]
        backend = RecordingBackend({
            batch_id: {"items": [{"gene_id": "G2", "contradicted_mechanism_ids": ["mech_002"]}], "input_tokens": 4, "output_tokens": 2}
        })
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            account = BudgetAccount()
            account.declare({"llm_call_limit": 10, "token_limit": 1000, "sequence_query_limit": 0, "wall_time_limit_seconds": None, "peak_memory_limit_mb": None})
            run_cluster_review(batches, account, out, CASE_ID, backend, mechanisms={"mech_002"})
            frame = pd.read_csv(out / "llm_review.tsv", sep="\t", dtype=str).fillna("")
            row = frame[frame["gene_id"] == "G2"].iloc[0]
            self.assertEqual("reviewed", row["review_status"])
            self.assertEqual("mech_002", row["contradicted_mechanism_ids"])
            self.assertEqual("", row["supported_mechanism_ids"])

    def test_hallucinated_source_and_unknown_mechanism_rejected(self) -> None:
        """Case 3."""
        handoff_rows, plan_rows, cluster_rows = _basic_fixture()
        batches = build_review_batches(handoff_rows, plan_rows, cluster_rows)
        batch_id = batches[0]["batch_id"]
        backend = RecordingBackend({
            batch_id: {
                "items": [
                    {"gene_id": "G1", "supported_mechanism_ids": ["mech_999"], "rationale_source_ids": ["SRC1"]},
                    {"gene_id": "G4", "supported_mechanism_ids": ["mech_001"], "rationale_source_ids": ["FAKE_SRC"]},
                ],
                "input_tokens": 20, "output_tokens": 10,
            }
        })
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            account = BudgetAccount()
            account.declare({"llm_call_limit": 10, "token_limit": 1000, "sequence_query_limit": 0, "wall_time_limit_seconds": None, "peak_memory_limit_mb": None})
            manifest = run_cluster_review(
                batches, account, out, CASE_ID, backend,
                mechanisms={"mech_001"}, allowed_source_ids={"SRC1"},
            )
            frame = pd.read_csv(out / "llm_review.tsv", sep="\t", dtype=str).fillna("")
            g1 = frame[frame["gene_id"] == "G1"].iloc[0]
            g4 = frame[frame["gene_id"] == "G4"].iloc[0]
            self.assertEqual("rejected_unverified", g1["review_status"])
            self.assertEqual("", g1["supported_mechanism_ids"])
            self.assertEqual("rejected_unverified", g4["review_status"])
            self.assertEqual("", g4["supported_mechanism_ids"])
            reasons = {(entry["gene_id"], entry["reason"]) for entry in manifest["rejected_claims"]}
            self.assertIn(("G1", "unknown_mechanism_id"), reasons)
            self.assertIn(("G4", "unverified_source_id"), reasons)
            self.assertEqual(2, manifest["counts"]["rejected"])

    def test_empty_response_never_fabricated(self) -> None:
        """Case 6."""
        handoff_rows, plan_rows, cluster_rows = _basic_fixture()
        batches = build_review_batches(handoff_rows, plan_rows, cluster_rows)
        backend = RecordingBackend({})  # no canned entry for this batch_id at all
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            account = BudgetAccount()
            account.declare({"llm_call_limit": 10, "token_limit": 1000, "sequence_query_limit": 0, "wall_time_limit_seconds": None, "peak_memory_limit_mb": None})
            manifest = run_cluster_review(batches, account, out, CASE_ID, backend)
            frame = pd.read_csv(out / "llm_review.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue((frame["review_status"] == "empty_response").all())
            self.assertTrue((frame["supported_mechanism_ids"] == "").all())
            self.assertEqual(len(frame), manifest["counts"]["empty"])

    def test_budget_partial_run_defers_remaining_with_truncation_point(self) -> None:
        """Case 5."""
        handoff_rows, plan_rows, cluster_rows = _basic_fixture()
        batches = build_review_batches(handoff_rows, plan_rows, cluster_rows, batch_size=1)
        self.assertEqual(2, len(batches))
        backend = RecordingBackend({b["batch_id"]: {"items": [], "input_tokens": 1, "output_tokens": 1} for b in batches})
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            account = BudgetAccount()
            account.declare({"llm_call_limit": 1, "token_limit": 1000, "sequence_query_limit": 0, "wall_time_limit_seconds": None, "peak_memory_limit_mb": None})
            manifest = run_cluster_review(batches, account, out, CASE_ID, backend)
            self.assertEqual(1, len(backend.calls))
            frame = pd.read_csv(out / "llm_review.tsv", sep="\t", dtype=str).fillna("")
            self.assertGreater(int((frame["review_status"] == "deferred_budget").sum()), 0)
            self.assertEqual(1, len(manifest["truncation_points"]))
            self.assertEqual("llm_call", manifest["truncation_points"][0]["kind"])
            frozen = account.freeze()
            self.assertEqual(1, frozen["usage"]["llm_call_count"])
            self.assertTrue(frozen["exhausted"]["llm_call_limit"])
            self.assertGreaterEqual(frozen["remaining"]["llm_call_limit"], 0)

    def test_partial_resume_skips_reviewed_batches(self) -> None:
        """Case 8."""
        handoff_rows, plan_rows, cluster_rows = _basic_fixture()
        batches = build_review_batches(handoff_rows, plan_rows, cluster_rows, batch_size=1)
        responses = {b["batch_id"]: {"items": [], "input_tokens": 1, "output_tokens": 1} for b in batches}
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            backend1 = RecordingBackend(responses)
            account1 = BudgetAccount()
            account1.declare({"llm_call_limit": 1, "token_limit": 1000, "sequence_query_limit": 0, "wall_time_limit_seconds": None, "peak_memory_limit_mb": None})
            run_cluster_review(batches, account1, out, CASE_ID, backend1)
            self.assertEqual(1, len(backend1.calls))

            backend2 = RecordingBackend(responses)
            account2 = BudgetAccount()
            account2.declare({"llm_call_limit": 5, "token_limit": 1000, "sequence_query_limit": 0, "wall_time_limit_seconds": None, "peak_memory_limit_mb": None})
            manifest2 = run_cluster_review(batches, account2, out, CASE_ID, backend2, resume=True)
            # only the batch that was NOT reviewed on the first pass gets a call now
            self.assertEqual(1, len(backend2.calls))
            self.assertEqual(1, manifest2["resume"]["resumed_batch_count"])
            frozen2 = account2.freeze()
            self.assertEqual(1, frozen2["usage"]["llm_call_count"])
            frame = pd.read_csv(out / "llm_review.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(3, len(frame))  # every gene across both batches accounted for

    def test_warm_cache_reuse_zero_budget_rerun(self) -> None:
        """Case 9."""
        handoff_rows, plan_rows, cluster_rows = _basic_fixture()
        batches = build_review_batches(handoff_rows, plan_rows, cluster_rows)
        batch_id = batches[0]["batch_id"]
        responses = {batch_id: {"items": [{"gene_id": "G1", "supported_mechanism_ids": ["mech_001"]}], "input_tokens": 3, "output_tokens": 2}}
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            from candidate_qtl.budget_cache import ContentCache

            cache = ContentCache(root / "cache")
            out1 = root / "run1"
            backend1 = RecordingBackend(responses)
            account1 = BudgetAccount()
            account1.declare({"llm_call_limit": 10, "token_limit": 1000, "sequence_query_limit": 0, "wall_time_limit_seconds": None, "peak_memory_limit_mb": None})
            run_cluster_review(batches, account1, out1, CASE_ID, backend1, mechanisms={"mech_001"}, cache=cache)
            self.assertEqual(1, len(backend1.calls))

            out2 = root / "run2"  # fresh out-dir: no llm_review.tsv to resume from
            backend2 = RecordingBackend(responses)
            account2 = BudgetAccount()
            account2.declare({"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 0, "wall_time_limit_seconds": None, "peak_memory_limit_mb": None})
            manifest2 = run_cluster_review(batches, account2, out2, CASE_ID, backend2, mechanisms={"mech_001"}, cache=cache)
            self.assertEqual(0, len(backend2.calls))
            frozen2 = account2.freeze()
            self.assertEqual(0, frozen2["usage"]["llm_call_count"])
            frame2 = pd.read_csv(out2 / "llm_review.tsv", sep="\t", dtype=str).fillna("")
            g1_row = frame2[frame2["gene_id"] == "G1"].iloc[0]
            self.assertEqual("reviewed", g1_row["review_status"])
            self.assertEqual("mech_001", g1_row["supported_mechanism_ids"])
            self.assertEqual(1, manifest2["batch_counts"]["cache_hit"])


class ReviewAuditRemediationTests(unittest.TestCase):
    """Regressions for the PR #392 audit findings."""

    def _declare(self, **overrides) -> BudgetAccount:
        limits = {"llm_call_limit": 10, "token_limit": 1000, "sequence_query_limit": 0,
                  "wall_time_limit_seconds": None, "peak_memory_limit_mb": None}
        limits.update(overrides)
        account = BudgetAccount()
        account.declare(limits)
        return account

    def test_cache_key_includes_the_validation_allowlists(self) -> None:
        """P1: a review cached with no allowlist must not be reused as `reviewed`
        under a later, restrictive allowlist -- it would smuggle a claim past a
        validation that never ran on it."""
        handoff_rows, plan_rows, cluster_rows = _basic_fixture()
        batches = build_review_batches(handoff_rows, plan_rows, cluster_rows)
        batch_id = batches[0]["batch_id"]
        responses = {batch_id: {"items": [{"gene_id": "G1", "supported_mechanism_ids": ["mech_999"]}], "input_tokens": 3, "output_tokens": 2}}
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            from candidate_qtl.budget_cache import ContentCache

            cache = ContentCache(root / "cache")
            first = RecordingBackend(responses)
            run_cluster_review(batches, self._declare(), root / "warm", CASE_ID, first, cache=cache)
            frame = pd.read_csv(root / "warm" / "llm_review.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("reviewed", frame[frame["gene_id"] == "G1"].iloc[0]["review_status"])

            # same prompt, same provider/model -- but now mech_999 is not on the
            # allowlist, so the cache must miss and the claim must be rejected
            second = RecordingBackend(responses)
            run_cluster_review(
                batches, self._declare(), root / "strict", CASE_ID, second,
                mechanisms={"mech_001"}, cache=cache,
            )
            self.assertEqual(1, len(second.calls))
            strict = pd.read_csv(root / "strict" / "llm_review.tsv", sep="\t", dtype=str).fillna("")
            row = strict[strict["gene_id"] == "G1"].iloc[0]
            self.assertEqual("rejected_unverified", row["review_status"])
            self.assertEqual("", row["supported_mechanism_ids"])

    def test_resume_refuses_a_stale_row_when_the_request_changed(self) -> None:
        """P1: the same batch_id under a changed responses file / provider must
        be re-run, not answered from the previous run's row."""
        handoff_rows, plan_rows, cluster_rows = _basic_fixture()
        batches = build_review_batches(handoff_rows, plan_rows, cluster_rows)
        batch_id = batches[0]["batch_id"]
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            first = RecordingBackend({batch_id: {"items": [{"gene_id": "G1", "supported_mechanism_ids": ["mech_001"]}], "input_tokens": 2, "output_tokens": 1}})
            run_cluster_review(batches, self._declare(), out, CASE_ID, first)
            fingerprints = set(pd.read_csv(out / "llm_review.tsv", sep="\t", dtype=str).fillna("")["request_fingerprint"])
            self.assertEqual(1, len(fingerprints))
            self.assertNotIn("", fingerprints)

            # identical request -> genuinely resumable, no second call
            same = RecordingBackend({batch_id: {"items": [], "input_tokens": 2, "output_tokens": 1}})
            resumed = run_cluster_review(batches, self._declare(), out, CASE_ID, same, resume=True)
            self.assertEqual([], same.calls)
            self.assertEqual(1, resumed["batch_counts"]["resumed"])

            # a stricter allowlist changes the request identity -> re-run
            changed = RecordingBackend({batch_id: {"items": [], "input_tokens": 2, "output_tokens": 1}})
            rerun = run_cluster_review(
                batches, self._declare(), out, CASE_ID, changed,
                mechanisms={"mech_001"}, resume=True,
            )
            self.assertEqual(1, len(changed.calls))
            self.assertEqual(0, rerun["batch_counts"]["resumed"])

    def test_resumed_rows_are_counted_like_cache_hits(self) -> None:
        """P2: a fully resumed run must not report zero reviewed rows while
        llm_review.tsv holds them."""
        handoff_rows, plan_rows, cluster_rows = _basic_fixture()
        batches = build_review_batches(handoff_rows, plan_rows, cluster_rows)
        batch_id = batches[0]["batch_id"]
        responses = {batch_id: {"items": [{"gene_id": "G1", "supported_mechanism_ids": ["mech_001"]}], "input_tokens": 2, "output_tokens": 1}}
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            run_cluster_review(batches, self._declare(), out, CASE_ID, RecordingBackend(responses))
            second = RecordingBackend(responses)
            manifest = run_cluster_review(batches, self._declare(), out, CASE_ID, second, resume=True)
            self.assertEqual([], second.calls)
            frame = pd.read_csv(out / "llm_review.tsv", sep="\t", dtype=str).fillna("")
            counted = sum(manifest["counts"][key] for key in ("reviewed", "rejected", "empty", "deferred", "not_evaluated_no_budget", "answered_over_budget"))
            self.assertEqual(len(frame), counted)
            self.assertEqual(1, manifest["counts"]["reviewed"])

    def test_paid_response_over_token_budget_is_not_bought_twice(self) -> None:
        """P2: when the token gate refuses an answer that was already paid for,
        the claims are refused but the batch is terminal, so a resume does not
        call the provider again."""
        handoff_rows, plan_rows, cluster_rows = _basic_fixture()
        batches = build_review_batches(handoff_rows, plan_rows, cluster_rows)
        batch_id = batches[0]["batch_id"]
        responses = {batch_id: {"items": [{"gene_id": "G1", "supported_mechanism_ids": ["mech_001"]}], "input_tokens": 1, "output_tokens": 500}}
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            backend = RecordingBackend(responses)
            manifest = run_cluster_review(batches, self._declare(token_limit=10), out, CASE_ID, backend)
            self.assertEqual(1, len(backend.calls))
            frame = pd.read_csv(out / "llm_review.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue((frame["review_status"] == "answered_over_budget").all())
            # the refused claims never became support
            self.assertTrue((frame["supported_mechanism_ids"] == "").all())
            self.assertEqual("token", manifest["truncation_points"][0]["kind"])

            resumed_backend = RecordingBackend(responses)
            resumed = run_cluster_review(batches, self._declare(token_limit=10_000), out, CASE_ID, resumed_backend, resume=True)
            self.assertEqual([], resumed_backend.calls)
            self.assertEqual(1, resumed["batch_counts"]["resumed"])

    def test_wall_budget_stops_the_review_stage(self) -> None:
        """P1: the declared wall budget covers the LLM stage, not only Tier 2."""
        handoff_rows, plan_rows, cluster_rows = _basic_fixture()
        batches = build_review_batches(handoff_rows, plan_rows, cluster_rows, batch_size=1)
        self.assertEqual(2, len(batches))
        responses = {b["batch_id"]: {"items": [], "input_tokens": 1, "output_tokens": 1} for b in batches}
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            backend = RecordingBackend(responses)
            account = self._declare(wall_time_limit_seconds=0)
            manifest = run_cluster_review(batches, account, out, CASE_ID, backend)
            self.assertEqual([], backend.calls)
            self.assertEqual("wall_time", manifest["truncation_points"][0]["kind"])
            frame = pd.read_csv(out / "llm_review.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue((frame["review_status"] == "deferred_budget").all())

    def test_peak_memory_limit_stops_the_review_stage_and_reaches_the_ledger(self) -> None:
        """P1: peak memory is a hard limit for the review stage too, and the
        observed peak lands in the ledger instead of staying at 0."""
        handoff_rows, plan_rows, cluster_rows = _basic_fixture()
        batches = build_review_batches(handoff_rows, plan_rows, cluster_rows, batch_size=1)
        responses = {b["batch_id"]: {"items": [], "input_tokens": 1, "output_tokens": 1} for b in batches}
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            backend = RecordingBackend(responses)
            account = self._declare(peak_memory_limit_mb=100)
            with mock.patch("candidate_qtl.tier2_planner._current_rss_mb", return_value=(99999.0, None)):
                manifest = run_cluster_review(batches, account, out, CASE_ID, backend)
            self.assertEqual([], backend.calls)
            self.assertEqual("peak_memory", manifest["truncation_points"][0]["kind"])
            frozen = account.freeze()
            self.assertEqual(99999.0, frozen["usage"]["peak_memory_mb"])
            self.assertTrue(frozen["exhausted"]["peak_memory_limit_mb"])

    def test_bridge_rows_make_route_unique_rescue_reachable(self) -> None:
        """P2: with the canonical bridge supplied, route slots are actually
        allocated instead of always reporting no_route_information_available."""
        handoff_rows, plan_rows, cluster_rows = _basic_fixture()
        bridge_rows = [
            {"gene_id": "G1", "family_id": "cls_A", "route": "a", "mapping_tier": "tier1_source_backed_unique"},
            {"gene_id": "G4", "family_id": "cls_C", "route": "b", "mapping_tier": "tier3_sequence_match_unique"},
        ]
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            result = reduce_portfolio(
                handoff_rows, plan_rows, cluster_rows, out, CASE_ID,
                bridge_rows=bridge_rows, portfolio_size=20,
            )
            reasons = {str(entry.get("reason", "")) for entry in result["manifest"]["unfilled_slots"]}
            self.assertNotIn("no_route_information_available", reasons)
            frame = pd.read_csv(out / "deterministic_top20.tsv", sep="\t", dtype=str).fillna("")
            states = dict(zip(frame["gene_id"], frame["route_rescue_state"]))
            self.assertNotEqual("not_evaluated", states["G1"])
            self.assertNotEqual("not_evaluated", states["G4"])


class BuildCriticTests(unittest.TestCase):
    def test_fixation_finding_emitted(self) -> None:
        """Case 4."""
        handoff_rows, plan_rows, cluster_rows = _basic_fixture()
        batches = build_review_batches(handoff_rows, plan_rows, cluster_rows)
        review_rows = [
            {"batch_id": batches[0]["batch_id"], "cluster_id": "cls_A", "gene_id": "G1", "review_status": "reviewed", "supported_mechanism_ids": "mech_001", "contradicted_mechanism_ids": "", "rationale_source_ids": ""},
            {"batch_id": batches[0]["batch_id"], "cluster_id": "cls_A", "gene_id": "G2", "review_status": "reviewed", "supported_mechanism_ids": "mech_001", "contradicted_mechanism_ids": "", "rationale_source_ids": ""},
            {"batch_id": batches[0]["batch_id"], "cluster_id": "cls_C", "gene_id": "G4", "review_status": "reviewed", "supported_mechanism_ids": "mech_002", "contradicted_mechanism_ids": "", "rationale_source_ids": ""},
        ]
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            account = BudgetAccount()
            account.declare({"llm_call_limit": 10, "token_limit": 1000, "sequence_query_limit": 0, "wall_time_limit_seconds": None, "peak_memory_limit_mb": None})
            result = build_critic(review_rows, batches, handoff_rows, cluster_rows, account, out, CASE_ID)
            kinds = {finding["finding_kind"] for finding in result["findings"]}
            self.assertIn("fixation", kinds)
            fixation = next(f for f in result["findings"] if f["finding_kind"] == "fixation")
            self.assertEqual("mech_001", fixation["affected_ids"])
            self.assertEqual("none", fixation["rank_contribution"])
            frame = pd.read_csv(out / "llm_review_critic.tsv", sep="\t", dtype=str).fillna("")
            self.assertGreater(len(frame), 0)

    def test_clean_run_writes_header_only_critic(self) -> None:
        handoff_rows, plan_rows, cluster_rows = _basic_fixture()
        batches = build_review_batches(handoff_rows, plan_rows, cluster_rows)
        review_rows = [
            {"batch_id": batches[0]["batch_id"], "cluster_id": member["cluster_id"], "gene_id": member["gene_id"],
             "review_status": "reviewed", "supported_mechanism_ids": f"mech_{member['gene_id']}",
             "contradicted_mechanism_ids": "", "rationale_source_ids": ""}
            for member in batches[0]["members"]
        ]
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            account = BudgetAccount()
            account.declare({"llm_call_limit": 10, "token_limit": 1000, "sequence_query_limit": 0, "wall_time_limit_seconds": None, "peak_memory_limit_mb": None})
            result = build_critic(review_rows, batches, handoff_rows, cluster_rows, account, out, CASE_ID)
            self.assertEqual([], result["findings"])
            self.assertTrue((out / "llm_review_critic.tsv").exists())
            frame = pd.read_csv(out / "llm_review_critic.tsv", sep="\t", dtype=str)
            self.assertEqual(0, len(frame))
            self.assertEqual(list(frame.columns), ["finding_kind", "scope", "detail", "affected_ids", "severity", "rank_contribution"])


class AuditProposerDiversityTests(unittest.TestCase):
    def test_never_mutates_input_and_writes_tsv(self) -> None:
        mechanism_rows = [
            {"hypothesis_id": "mech_001", "gene_family": "FAM1", "reference_gene_name": "RefA", "source_ids": "PMID:1"},
            {"hypothesis_id": "mech_002", "gene_family": "FAM2", "reference_gene_name": "RefB", "source_ids": "PMID:2"},
            {"hypothesis_id": "mech_003", "gene_family": "FAM3", "reference_gene_name": "RefC", "source_ids": "PMID:3"},
        ]
        before = json.dumps(mechanism_rows, sort_keys=True)
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            result = audit_proposer_diversity(mechanism_rows, out, CASE_ID)
            self.assertEqual(before, json.dumps(mechanism_rows, sort_keys=True))
            self.assertEqual("pass", result["overall_verdict"])
            frame = pd.read_csv(out / "proposer_diversity_audit.tsv", sep="\t", dtype=str)
            self.assertIn("distinct_mechanism_count", set(frame["metric"]))

    def test_single_famous_gene_fixation_detected(self) -> None:
        mechanism_rows = [
            {"hypothesis_id": f"mech_{i:03d}", "gene_family": f"FAM{i}", "reference_gene_name": "SameGene", "source_ids": f"PMID:{i}"}
            for i in range(1, 5)
        ]
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            result = audit_proposer_diversity(mechanism_rows, out, CASE_ID)
            self.assertEqual("single_famous_gene_fixation", result["overall_verdict"])


class ReducePortfolioTests(unittest.TestCase):
    def test_budget_zero_deterministic_slate_and_all_not_evaluated(self) -> None:
        """Case 7 (portfolio half): zero review input -> llm_assisted identical to deterministic."""
        handoff_rows, plan_rows, cluster_rows = _basic_fixture()
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            result = reduce_portfolio(handoff_rows, plan_rows, cluster_rows, out, CASE_ID, review_rows=None, portfolio_size=20)
            # with no bridge/route information the unfilled route slot is
            # reported, never silently treated as satisfied
            self.assertIn(
                "no_route_information_available",
                {str(entry.get("reason", "")) for entry in result["manifest"]["unfilled_slots"]},
            )
            det = pd.read_csv(out / "deterministic_top20.tsv", sep="\t", dtype=str).fillna("")
            llm = pd.read_csv(out / "llm_assisted_top20.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue((llm["llm_status"] == "not_evaluated_no_budget").all())
            self.assertEqual(list(det["gene_id"]), list(llm["gene_id"]))
            self.assertEqual(sorted(det["gene_id"]), sorted(row["gene_id"] for row in handoff_rows))
            self.assertTrue((det["rank_contribution"] == "none").all())
            self.assertTrue((llm["rank_contribution"] == "none").all())
            self.assertNotIn("score", det.columns)
            self.assertNotIn("probability", det.columns)
            self.assertNotIn("weight", det.columns)

    def test_llm_reorder_never_adds_or_removes_genes(self) -> None:
        handoff_rows, plan_rows, cluster_rows = _basic_fixture()
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            review_rows = [{"gene_id": "G4", "review_status": "reviewed", "supported_mechanism_ids": "mech_001"}]
            reduce_portfolio(handoff_rows, plan_rows, cluster_rows, out, CASE_ID, review_rows=review_rows, portfolio_size=20)
            det = pd.read_csv(out / "deterministic_top20.tsv", sep="\t", dtype=str).fillna("")
            llm = pd.read_csv(out / "llm_assisted_top20.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(set(det["gene_id"]), set(llm["gene_id"]))
            g4 = llm[llm["gene_id"] == "G4"].iloc[0]
            self.assertEqual("llm_supported_promotion", g4["selection_basis"])
            self.assertEqual("reviewed", g4["llm_status"])

    def test_diversity_slot_policy_rescue_dependency_family_cap(self) -> None:
        """Case 10 (part 1): rescue slots, dependency dedup, family cap all observable."""
        handoff_rows = [
            _handoff_row("R1", "cls_r1", "rescue", "unclassified", "rescue_pending_review"),
            _handoff_row("R2", "cls_r2", "rescue", "weak_annotation", "rescue_pending_review"),
            _handoff_row("D1", "cls_d1", "priority", "", "priority_review"),
            _handoff_row("D2", "cls_d2", "priority", "", "priority_review"),
            _handoff_row("F1", "cls_fam", "priority", "", "priority_review"),
            _handoff_row("F2", "cls_fam", "priority", "", "priority_review"),
            _handoff_row("F3", "cls_fam", "priority", "", "priority_review"),
            _handoff_row("F4", "cls_fam", "priority", "", "priority_review"),
        ]
        plan_rows = [_plan_row(row["gene_id"], row["cluster_id"]) for row in handoff_rows]
        cluster_rows = [
            _cluster_row("cls_r1", "R1", "rescue", rescue_bucket="unclassified"),
            _cluster_row("cls_r2", "R2", "rescue", rescue_bucket="weak_annotation"),
            _cluster_row("cls_d1", "D1", "priority", dependency_cluster_ids="dep_x"),
            _cluster_row("cls_d2", "D2", "priority", dependency_cluster_ids="dep_x"),
            _cluster_row("cls_fam", "F1;F2;F3;F4", "priority"),
        ]
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            result = reduce_portfolio(handoff_rows, plan_rows, cluster_rows, out, CASE_ID, portfolio_size=9)
            det = pd.read_csv(out / "deterministic_top20.tsv", sep="\t", dtype=str).fillna("")
            by_gene = {row["gene_id"]: row for _, row in det.iterrows()}
            self.assertEqual("rescue_slot", by_gene["R1"]["slot_kind"])
            self.assertEqual("rescue_slot", by_gene["R2"]["slot_kind"])
            self.assertEqual("dependency_dedup", by_gene["D1"]["slot_kind"])
            self.assertEqual("dep_x", by_gene["D1"]["dependency_cluster_id"])
            self.assertEqual("dep_x", by_gene["D2"]["dependency_cluster_id"])
            fam_rows = det[det["cluster_id"] == "cls_fam"]
            self.assertEqual(4, len(fam_rows))
            non_overflow = fam_rows[fam_rows["slot_reason"].str.startswith("diversity_key:")]
            self.assertLessEqual(len(non_overflow), 3)
            overflow = fam_rows[fam_rows["slot_reason"].str.startswith("diversity_key_overflow:")]
            self.assertEqual(1, len(overflow))
            unfilled_kinds = {entry["slot_kind"] for entry in result["manifest"]["unfilled_slots"]}
            self.assertIn("route_unique_rescue", unfilled_kinds)

    def test_paralog_ambiguity_guard(self) -> None:
        """Case 10 (part 2): ambiguous paralog never SOLE at slot_index 1; co-listed case preserved."""
        handoff_rows = [
            _handoff_row("P1", "cls_p1", "priority", "", "priority_review"),
            _handoff_row("P2", "cls_p2", "priority", "", "priority_review"),
            _handoff_row("P3", "cls_p2", "priority", "", "priority_review"),
        ]
        plan_rows = [
            _plan_row("P1", "cls_p1", ambiguity_status="ambiguous_homology"),
            _plan_row("P2", "cls_p2", ambiguity_status="ambiguous_homology"),
            _plan_row("P3", "cls_p2", ambiguity_status="not_ambiguous"),
        ]
        cluster_rows = [
            _cluster_row("cls_p1", "P1", "priority"),
            _cluster_row("cls_p2", "P2;P3", "priority"),
        ]
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            reduce_portfolio(handoff_rows, plan_rows, cluster_rows, out, CASE_ID, portfolio_size=3)
            det = pd.read_csv(out / "deterministic_top20.tsv", sep="\t", dtype=str).fillna("")
            slot1 = det[det["slot_index"] == "1"].iloc[0]
            self.assertNotEqual("P1", slot1["gene_id"])
            p1 = det[det["gene_id"] == "P1"].iloc[0]
            p2 = det[det["gene_id"] == "P2"].iloc[0]
            p3 = det[det["gene_id"] == "P3"].iloc[0]
            self.assertEqual("yes", p1["paralog_ambiguity"])
            self.assertEqual("paralog_ambiguous_no_independent_support", p1["slot_reason"])
            self.assertNotEqual("1", p1["slot_index"])
            self.assertEqual("yes", p2["paralog_ambiguity"])
            self.assertEqual("paralog_ambiguous_co_listed", p2["slot_reason"])
            self.assertEqual("no", p3["paralog_ambiguity"])

    def test_rank_invariance_no_rank_impacting_files_and_rank_contribution_none(self) -> None:
        """Case 11."""
        handoff_rows, plan_rows, cluster_rows = _basic_fixture()
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            reduce_portfolio(handoff_rows, plan_rows, cluster_rows, out, CASE_ID, portfolio_size=20)
            for path in out.rglob("*"):
                if path.is_file():
                    relative = path.relative_to(out).as_posix()
                    self.assertFalse(_matches_rank_impacting(relative), f"{relative} matches a rank-impacting pattern")
            for name in ("deterministic_top20.tsv", "llm_assisted_top20.tsv"):
                # the artifact names themselves must stay off the rank-impacting
                # vocabulary the reachability gate enforces
                self.assertFalse(_matches_rank_impacting(name))
                self.assertFalse(any(fnmatch.fnmatch(name, pattern) for pattern in RANK_IMPACTING_PATTERNS))
                frame = pd.read_csv(out / name, sep="\t", dtype=str)
                self.assertTrue((frame["rank_contribution"] == "none").all())
                self.assertNotIn("score", frame.columns)


class FileReviewBackendTests(unittest.TestCase):
    def test_missing_batch_id_yields_empty_response_never_raises(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            responses_path = Path(temp) / "responses.json"
            responses_path.write_text(json.dumps({"batch_known": {"items": [{"gene_id": "G1"}], "input_tokens": 1, "output_tokens": 1}}), encoding="utf-8")
            backend = FileReviewBackend(responses_path)
            self.assertTrue(backend.available)
            result = backend.review({"batch_id": "batch_unknown"})
            self.assertEqual([], result["items"])
            self.assertEqual(0, result["input_tokens"])

    def test_known_batch_id_returns_canned_items(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            responses_path = Path(temp) / "responses.json"
            responses_path.write_text(json.dumps({"batch_known": {"items": [{"gene_id": "G1"}], "input_tokens": 4, "output_tokens": 2}}), encoding="utf-8")
            backend = FileReviewBackend(responses_path)
            result = backend.review({"batch_id": "batch_known"})
            self.assertEqual(1, len(result["items"]))
            self.assertEqual(4, result["input_tokens"])


if __name__ == "__main__":
    unittest.main()
