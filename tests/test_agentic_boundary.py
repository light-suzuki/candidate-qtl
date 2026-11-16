from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.agentic_boundary import merge, prepare
from scripts.agentic_boundary_material import restore_reading_shard, expand_reading_view


class AgenticBoundaryTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.boundary = self.root / "boundary.tsv"
        pd.DataFrame([{"candidate_gene_id": f"G{i}", "relevance_tier": "T2", "reason": "provided"} for i in range(1, 4)]).to_csv(self.boundary, sep="\t", index=False)
        self.higher = self.root / "higher.tsv"
        pd.DataFrame([{"candidate_gene_id": "H1", "relevance_tier": "T1"}]).to_csv(self.higher, sep="\t", index=False)
        self.request = self.root / "adaptive_synthesis_request.json"
        self.request.write_text(json.dumps({
            "contract": "boundary-tier-review-request-1.0.0", "stage": "boundary_review_request", "state": "awaiting_review", "case_id": "case-x",
            "source_wet_sha256": self.sha(self.boundary), "candidate_count": 3,
            "candidate_gene_ids": ["G1", "G2", "G3"], "candidate_input": {"path": str(self.boundary), "sha256": self.sha(self.boundary)},
            "higher_tier_context": {"state": "read_only", "candidate_input": {"path": str(self.higher), "sha256": self.sha(self.higher)}, "candidate_gene_ids": ["H1"]}
        }), encoding="utf-8")
        self.prepared = self.root / "prepared"
        self.manifest = prepare(str(self.request), str(self.prepared), max_candidates=2)

    def tearDown(self):
        self.tmp.cleanup()

    @staticmethod
    def sha(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def test_prepare_bounds_roster_and_keeps_higher_context_read_only(self):
        self.assertEqual(2, len(self.manifest["shards"]))
        first = json.loads((self.prepared / self.manifest["shards"][0]["file"]).read_text())
        self.assertEqual(["G1", "G2"], first["candidate_gene_ids"])
        self.assertEqual("read_only", first["higher_tier_context"]["state"])
        self.assertEqual(["H1"], first["higher_tier_context"]["candidate_gene_ids"])

    @staticmethod
    def paper_panel(scope, candidate_ids, count=5):
        return {
            "paper_handoff_policy": "tier_ranked_primary_papers_v1",
            "paper_handoff_scope": scope, "candidate_gene_ids": candidate_ids,
            "papers": [
                {"paper_handoff_rank": str(rank), "source_ids": [f"PMID:{scope}-{rank}"],
                 "candidate_gene_ids": candidate_ids, "record_state": "original_reports",
                 "reports": [
                     {"candidate_gene_ids": [cid], "paper_report": {
                         "source_id": f"PMID:{scope}-{rank}", "reported_gene": f"reference-{scope}",
                         "initial_tier": f"T{rank}", "initial_tier_reason": f"original reason {rank}",
                         "source_excerpt": f"original observation {scope} {rank}",
                         "paper_handoff_policy": "tier_ranked_primary_papers_v1",
                         "paper_handoff_scope": scope, "paper_handoff_state": "primary",
                         "paper_handoff_rank": str(rank)}} for cid in candidate_ids]}
                for rank in range(1, count + 1)]}

    def primary_request(self):
        request = json.loads(self.request.read_text())
        request["primary_paper_evidence"] = [self.paper_panel("shared", ["G1", "G3"]),
                                              self.paper_panel("other", ["G2"], 2)]
        unavailable = request["primary_paper_evidence"][1]["papers"][1]
        unavailable["record_state"] = "original_report_unavailable"
        unavailable["reports"] = []
        request["higher_tier_context"]["primary_paper_evidence"] = [self.paper_panel("higher", ["H1"])]
        self.request.write_text(json.dumps(request), encoding="utf-8")
        return request

    def test_selected_papers_reach_workers_with_original_facts_and_copy_bindings(self):
        request = self.primary_request()
        originals = {path: path.read_bytes() for path in [self.request, self.boundary, self.higher]}
        root = self.root / "with-primary-papers"
        manifest = prepare(str(self.request), str(root), max_candidates=2)
        for shard, entry in zip(manifest["shards"], manifest["reading_views"]):
            canonical = restore_reading_shard(root, manifest, shard["shard_id"])
            view = json.loads((root / entry["file"]).read_text())
            view = expand_reading_view(view, json.loads((root / "material_catalog.json").read_text())["materials"])
            assigned = set(shard["candidate_gene_ids"])
            panels = canonical["primary_paper_evidence"]
            self.assertEqual(panels, view["primary_paper_evidence"])
            self.assertEqual(request["higher_tier_context"]["primary_paper_evidence"],
                             view["higher_tier_context"]["primary_paper_evidence"])
            self.assertEqual("read_only", view["higher_tier_context"]["state"])
            expected_scopes = {panel["paper_handoff_scope"] for panel in request["primary_paper_evidence"]
                               if assigned.intersection(panel["candidate_gene_ids"])}
            self.assertEqual(expected_scopes, {panel["paper_handoff_scope"] for panel in panels})
            for panel in panels:
                original = next(p for p in request["primary_paper_evidence"]
                                if p["paper_handoff_scope"] == panel["paper_handoff_scope"])
                self.assertEqual(len(original["papers"]), len(panel["papers"]))
                self.assertEqual(assigned.intersection(original["candidate_gene_ids"]),
                                 set(panel["candidate_gene_ids"]))
                for paper, raw_paper in zip(panel["papers"], original["papers"]):
                    self.assertEqual(raw_paper["paper_handoff_rank"], paper["paper_handoff_rank"])
                    self.assertEqual(raw_paper["source_ids"], paper["source_ids"])
                    self.assertEqual(assigned.intersection(raw_paper["candidate_gene_ids"]),
                                     set(paper["candidate_gene_ids"]))
                    self.assertEqual([report for report in raw_paper["reports"]
                                      if assigned.intersection(report["candidate_gene_ids"])], paper["reports"])
        self.assertEqual(originals, {path: path.read_bytes() for path in originals})

    def test_reading_projection_rejects_lost_or_changed_primary_papers(self):
        self.primary_request()
        root = self.root / "primary-integrity"
        manifest = prepare(str(self.request), str(root), max_candidates=2)
        entry = manifest["reading_views"][0]
        path = root / entry["file"]
        original = path.read_text()
        for change in ["drop_set", "drop_paper", "higher_tier_reason"]:
            with self.subTest(change=change):
                view = json.loads(original)
                expanded = expand_reading_view(view, json.loads((root / "material_catalog.json").read_text())["materials"])
                if change == "drop_set":
                    del view["primary_paper_evidence"]
                elif change == "drop_paper":
                    expanded["primary_paper_evidence"][0]["papers"].pop()
                    view["primary_paper_evidence"] = {"$inline": expanded["primary_paper_evidence"]}
                else:
                    expanded["higher_tier_context"]["primary_paper_evidence"][0]["papers"][0]["reports"][0]["paper_report"]["initial_tier_reason"] = "rewritten"
                    view["higher_tier_context"] = {"$inline": expanded["higher_tier_context"]}
                path.write_text(json.dumps(view), encoding="utf-8")
                entry["sha256"] = self.sha(path)
                with self.assertRaisesRegex(ValueError, "does not restore canonical shard fields"):
                    restore_reading_shard(root, manifest, manifest["shards"][0]["shard_id"])

    def test_missing_original_for_one_copy_preserves_its_selected_paper_ids(self):
        request = self.primary_request()
        for paper in request["primary_paper_evidence"][0]["papers"]:
            paper["reports"] = [report for report in paper["reports"]
                                if report["candidate_gene_ids"] == ["G1"]]
        self.request.write_text(json.dumps(request), encoding="utf-8")
        root = self.root / "copy-without-original"
        manifest = prepare(str(self.request), str(root), max_candidates=2)
        canonical = restore_reading_shard(root, manifest, manifest["shards"][1]["shard_id"])
        self.assertEqual(["G3"], canonical["candidate_gene_ids"])
        papers = canonical["primary_paper_evidence"][0]["papers"]
        self.assertEqual([str(rank) for rank in range(1, 6)], [paper["paper_handoff_rank"] for paper in papers])
        for paper in papers:
            self.assertEqual(["G3"], paper["candidate_gene_ids"])
            self.assertEqual("original_report_unavailable", paper["record_state"])
            self.assertEqual([], paper["reports"])
            self.assertTrue(paper["source_ids"])

    def test_paper_material_budget_splits_work_without_dropping_papers(self):
        request = self.primary_request()
        for panel in request["primary_paper_evidence"]:
            for paper in panel["papers"]:
                for report in paper["reports"]:
                    report["paper_report"]["source_excerpt"] = "original facts " * 2000
        self.request.write_text(json.dumps(request), encoding="utf-8")
        root = self.root / "large-primary-papers"
        manifest = prepare(str(self.request), str(root), max_candidates=3, max_material_bytes=20000)
        self.assertEqual([["G1"], ["G2"], ["G3"]],
                         [shard["candidate_gene_ids"] for shard in manifest["shards"]])
        for shard, entry in zip(manifest["shards"], manifest["reading_views"]):
            canonical = restore_reading_shard(root, manifest, shard["shard_id"])
            self.assertEqual(2 if shard["candidate_gene_ids"] == ["G2"] else 5,
                             len(canonical["primary_paper_evidence"][0]["papers"]))
            self.assertTrue(entry["reading_budget_exceeded"])

    def context(self):
        path = self.root / "context.json"
        value = {"case_id": "case-x", "phenotype": "Recorded response in the developing organ",
                 "candidate_material": {cid: {"aligned_regions": [[10, 20]], "unverified": True}
                                        for cid in ["G1", "G2", "G3", "H1"]}}
        path.write_text(json.dumps(value), encoding="utf-8")
        return path, value

    def test_comparison_material_and_request_policy_reach_each_shard(self):
        context, value = self.context()
        request = json.loads(self.request.read_text())
        request["policy"] = ["independent experiment-priority comparison"]
        self.request.write_text(json.dumps(request), encoding="utf-8")
        result = prepare(str(self.request), str(self.root / "with-context"), 2,
                         comparison_context_path=str(context))
        first = json.loads((self.root / "with-context" / result["shards"][0]["file"]).read_text())
        self.assertEqual(value["phenotype"], first["phenotype"])
        self.assertEqual(request["policy"], first["request_policy"])
        self.assertEqual(self.sha(self.request), first["request_sha256"])
        self.assertEqual({cid: value["candidate_material"][cid] for cid in ["G1", "G2", "H1"]},
                         first["candidate_comparison_material"])
        self.assertEqual(self.sha(context), result["comparison_context"]["sha256"])
        self.assertEqual(self.sha(self.boundary), result["source_wet_sha256"])

    def test_comparison_context_requires_exact_case_and_scope(self):
        context, value = self.context()
        for key, wrong in [("case_id", "other"), ("phenotype", ""),
                           ("candidate_material", {"G1": {}})]:
            damaged = {**value, key: wrong}
            context.write_text(json.dumps(damaged), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "comparison context"):
                prepare(str(self.request), str(self.root / "rejected-context"),
                        comparison_context_path=str(context))
            self.assertFalse((self.root / "rejected-context").exists())

    def test_merge_binds_response_and_rejects_changed_comparison_material(self):
        context, _ = self.context()
        prepared = self.root / "bound-context"
        result = prepare(str(self.request), str(prepared), comparison_context_path=str(context))
        payload = {"case_id": "case-x", "source_wet_sha256": self.sha(self.boundary),
                   "shard_id": result["shards"][0]["shard_id"], "executor_id": "Codex",
                   "executor_model": "m", "responses": [
                       {"candidate_gene_id": cid, "decision": "defer", "reason": "comparative priority"}
                       for cid in ["G1", "G2", "G3"]]}
        response = self.root / "context-response.json"
        response.write_text(json.dumps(payload), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "comparison context SHA mismatch in response"):
            merge(str(prepared), [str(response)], str(self.root / "missing.tsv"), role="proposal",
                  executor_id="Codex", executor_model="m")
        payload["comparison_context_sha256"] = self.sha(context)
        response.write_text(json.dumps(payload), encoding="utf-8")
        sidecar = merge(str(prepared), [str(response)], str(self.root / "valid.tsv"),
                        role="proposal", executor_id="Codex", executor_model="m")
        self.assertEqual(result["comparison_context"], sidecar["comparison_context"])
        self.assertEqual(3, sidecar["state_counts"]["defer"])
        context.write_text(context.read_text() + "\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "comparison context SHA mismatch"):
            merge(str(prepared), [str(response)], str(self.root / "tampered.tsv"), role="proposal",
                  executor_id="Codex", executor_model="m")

    def test_independent_comparison_cannot_bind_to_another_context_proposal(self):
        context, value = self.context()
        prepared = self.root / "context-a"
        a = prepare(str(self.request), str(prepared), comparison_context_path=str(context))
        proposer = {"case_id": "case-x", "source_wet_sha256": self.sha(self.boundary),
                    "shard_id": a["shards"][0]["shard_id"], "executor_id": "Codex", "executor_model": "m",
                    "comparison_context_sha256": self.sha(context), "responses": [
                        {"candidate_gene_id": cid, "decision": "defer", "reason": "comparative priority"}
                        for cid in ["G1", "G2", "G3"]]}
        response = self.root / "proposer-a.json"
        response.write_text(json.dumps(proposer), encoding="utf-8")
        proposal = self.root / "proposal-a.tsv"
        merge(str(prepared), [str(response)], str(proposal), role="proposal", executor_id="Codex", executor_model="m")
        reviewer = {"case_id": "case-x", "source_wet_sha256": self.sha(self.boundary),
                    "shard_id": a["shards"][0]["shard_id"], "reviewer_id": "reviewer", "reviewer_model": "m",
                    "comparison_context_sha256": self.sha(context), "proposal_sha256": self.sha(proposal),
                    "responses": proposer["responses"]}
        review = self.root / "review.json"; review.write_text(json.dumps(reviewer), encoding="utf-8")
        merge(str(prepared), [str(review)], str(self.root / "valid-review.tsv"), role="review",
              reviewer_id="reviewer", reviewer_model="m", proposal_path=str(proposal))
        context_b = self.root / "context-b.json"
        value["candidate_material"]["G2"] = {"alternative_protein": "distinct function"}
        context_b.write_text(json.dumps(value), encoding="utf-8")
        prepared_b = self.root / "prepared-b"
        prepare(str(self.request), str(prepared_b), comparison_context_path=str(context_b))
        reviewer["comparison_context_sha256"] = self.sha(context_b)
        review.write_text(json.dumps(reviewer), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "proposal comparison context or lineage mismatch"):
            merge(str(prepared_b), [str(review)], str(self.root / "cross-context.tsv"), role="review",
                  reviewer_id="reviewer", reviewer_model="m", proposal_path=str(proposal))

    def test_proposal_merge_marks_missing_and_failed_not_evaluated(self):
        shard = self.manifest["shards"][0]
        response = self.root / "part.json"
        response.write_text(json.dumps({"case_id": "case-x", "source_wet_sha256": self.manifest["source_wet_sha256"], "shard_id": shard["shard_id"], "executor_id": "Codex", "executor_model": "m", "responses": [
            {"candidate_gene_id": "G1", "decision": "adopt", "reason": "reviewed"},
            {"candidate_gene_id": "G2", "status": "not_evaluated", "decision": "not_evaluated", "reason": "insufficient context"},
        ]}), encoding="utf-8")
        out = self.root / "proposal.tsv"
        sidecar = merge(str(self.prepared), [str(response)], str(out), role="proposal", executor_id="Codex", executor_model="m")
        rows = pd.read_csv(out, sep="\t", dtype=str).fillna("")
        self.assertEqual(["adopt", "not_evaluated", "not_evaluated"], rows["decision"].tolist())
        self.assertEqual("missing_or_failed_response", rows.iloc[2]["reason"])
        self.assertEqual("partial", sidecar["completion_status"])
        self.assertEqual("Codex", rows.iloc[0]["executor_id"])

    def test_review_mode_binds_to_proposal_and_keeps_review_distinct(self):
        proposal = self.root / "proposal.tsv"
        proposal.write_text("candidate_gene_id\tdecision\treason\tsource_wet_sha256\texecutor_id\texecutor_model\n", encoding="utf-8")
        shard = self.manifest["shards"][0]
        response = self.root / "review.json"
        response.write_text(json.dumps({"case_id": "case-x", "source_wet_sha256": self.manifest["source_wet_sha256"], "proposal_sha256": self.sha(proposal), "shard_id": shard["shard_id"], "reviewer_id": "reviewer", "reviewer_model": "rmodel", "responses": [
            {"candidate_gene_id": "G1", "decision": "defer", "reason": "independent reason"}
        ]}), encoding="utf-8")
        out = self.root / "independent.tsv"
        merge(str(self.prepared), [str(response)], str(out), role="review", reviewer_id="reviewer", reviewer_model="rmodel", proposal_path=str(proposal))
        row = pd.read_csv(out, sep="\t", dtype=str).fillna("").iloc[0]
        self.assertEqual("defer", row["decision"])
        self.assertEqual(self.sha(proposal), row["proposal_sha256"])
        self.assertEqual("reviewer", row["reviewer_id"])
        self.assertNotIn("executor_id", row.index)

    def test_unknown_duplicate_case_or_source_mismatch_is_rejected(self):
        shard = self.manifest["shards"][0]
        cases = [
            ({"case_id": "wrong", "source_wet_sha256": self.manifest["source_wet_sha256"]}, "case_id mismatch"),
            ({"case_id": "case-x", "source_wet_sha256": "bad"}, "source SHA mismatch"),
            ({"case_id": "case-x", "source_wet_sha256": self.manifest["source_wet_sha256"], "responses": [{"candidate_gene_id": "X", "decision": "adopt", "reason": "x"}]}, "unknown candidate ID"),
            ({"case_id": "case-x", "source_wet_sha256": self.manifest["source_wet_sha256"], "responses": [{"candidate_gene_id": "G1", "decision": "adopt", "reason": "x"}, {"candidate_gene_id": "G1", "decision": "defer", "reason": "y"}]}, "duplicate candidate response ID"),
        ]
        for payload, message in cases:
            payload["shard_id"] = shard["shard_id"]
            payload["executor_id"] = "Codex"
            payload["executor_model"] = "m"
            path = self.root / "bad.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, message):
                merge(str(self.prepared), [str(path)], str(self.root / "bad.tsv"), role="proposal", executor_id="Codex", executor_model="m")


if __name__ == "__main__":
    unittest.main()
