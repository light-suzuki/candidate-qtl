"""Real Wet handoff -> boundary review session integration regressions."""
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.claim_review_session import _files, finish_wet_review, sha256
from candidate_qtl.boundary_tier_review import UPPER_EVIDENCE_SCOPED_POLICY
from candidate_qtl.observation_tiers import CHECKLIST_FIELDS
from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from tests._issue664_helpers import manifest_for
from tests.test_boundary_tier_review import candidate_table_rows
from tests.test_wet_relevance_selection import handoff


def _write_completion(root: Path, wet: dict, execution_context: dict | None = None) -> Path:
    completion = root / "completion.json"
    payload = {
        "case_id": "case1", "wet": wet, "files": _files(root),
    }
    if execution_context is not None:
        payload["execution_context"] = execution_context
    completion.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return completion


def _review(scope_path: str, out: Path, *, adopt: int, disagreement: bool = False):
    scope = json.loads(Path(scope_path).read_text(encoding="utf-8"))
    boundary_path = Path(scope["artifacts"]["boundary"]["path"])
    boundary = pd.read_csv(boundary_path, sep="\t", dtype=str).fillna("")
    source_hash = hashlib.sha256(boundary_path.read_bytes()).hexdigest()
    proposal_rows = []
    for index, row in boundary.iterrows():
        proposal_rows.append({
            "candidate_gene_id": row["candidate_gene_id"],
            "decision": "adopt" if index < adopt else "defer",
            "reason": "fixture whole-tier comparison",
            "source_wet_sha256": source_hash,
            "executor_id": "fixture-executor", "executor_model": "fixture-model",
        })
    proposal = out / "proposal.tsv"
    pd.DataFrame(proposal_rows).to_csv(proposal, sep="\t", index=False)
    independent_rows = []
    for row in proposal_rows:
        independent_rows.append({
            **row,
            "decision": "adopt" if disagreement else row["decision"],
            "reason": "fixture independent comparison",
            "proposal_sha256": hashlib.sha256(proposal.read_bytes()).hexdigest(),
            "reviewer_id": "fixture-reviewer", "reviewer_model": "fixture-review-model",
        })
    independent = out / "independent.tsv"
    pd.DataFrame(independent_rows).to_csv(independent, sep="\t", index=False)
    return proposal, independent


class BoundarySessionIntegrationTests(unittest.TestCase):
    def test_six_column_report_follows_v3_handoff_and_actual_adopted_membership(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            first = base / "first"
            first.mkdir()
            genes = [f"PsC_{index:02d}" for index in range(26)]
            atlas = base / "atlas.tsv"
            pd.DataFrame([{"gene_id": gene, "gene_name": "PsC_name" if index == 0 else "",
                           "chrom": "chr05", "start": str(101 + 1000 * index), "end": str(300 + 1000 * index)}
                          for index, gene in enumerate([*genes, "PsC_lower"])]).to_csv(atlas, sep="\t", index=False)
            supported = base / "supported.tsv"
            pd.DataFrame(columns=["candidate_gene_id", "route", "support_status"]).to_csv(supported, sep="\t", index=False)
            checklist = {key: "no" for key in CHECKLIST_FIELDS}
            checklist.update(phenotype_reported="yes", equivalent_target_phenotype="yes", identity_resolved="yes")
            paper = {"source_id": "PMID:123456", "source_title": "Synthetic MtA phenotype report",
                     "reported_gene": "MtA", "reported_gene_id": "MtA.1", "reported_species": "Medicago truncatula",
                     "reference_gene_id": "AtB", "gene_correspondence": "homolog", "paper_checklist": checklist,
                     "initial_tier_reason": "Reported phenotype motivates the focal candidate hypothesis."}
            review_rows = [{"candidate_gene_id": gene, "claim_id": f"claim_{gene}", "route": "b",
                            "support_status": "unverified_route_b_review", "rank_contribution": "none",
                            "reference_gene": "AtB", "reference_species": "Arabidopsis thaliana",
                            "source_id": "PMID:123456", "paper_reports_json": json.dumps([paper])}
                           for gene in genes]
            review_rows.append({**review_rows[0], "candidate_gene_id": "PsC_lower", "claim_id": "claim_lower",
                                "paper_reports_json": json.dumps([{**paper, "paper_checklist": {
                                    **{key: "no" for key in CHECKLIST_FIELDS}, "function_supports_hypothesis": "yes"}}])})
            review = base / "review.tsv"
            pd.DataFrame(review_rows).to_csv(review, sep="\t", index=False)
            manifest = base / "review_manifest.json"
            manifest.write_text(json.dumps({"contract": "route-review-candidates-1.0.0",
                "stage": "route_b_review_candidate_freeze", "case_id": "case1", "inputs": {},
                "output": str(review), "artifact_sha256": sha256(review)}), encoding="utf-8")
            wet = wet_candidate_handoff(str(first / "wet"),
                route_b_candidates_arg=str(supported), route_b_candidates_manifest_arg=str(manifest_for(supported, route="b")),
                interval_atlas_arg=str(atlas), review_candidates_arg=str(review),
                review_candidates_manifest_arg=str(manifest), case_id="case1", max_adaptive_rounds=1)
            self.assertEqual(UPPER_EVIDENCE_SCOPED_POLICY, wet["first_stage_selection_policy_version"])
            completion = _write_completion(first, wet)
            frozen_bytes = {path: path.read_bytes() for path in first.rglob("*") if path.is_file()}
            source_bytes = {path: path.read_bytes() for path in (atlas, review)}
            first_rows = pd.read_csv(wet["wet_candidates"], sep="\t", dtype=str).fillna("")
            first_report = candidate_table_rows(Path(wet["report"]["path"]).read_text(encoding="utf-8"))
            self.assertEqual(genes, list(first_rows.candidate_gene_id))
            self.assertEqual(len(first_rows), len(first_report))
            self.assertTrue(all(row[5].startswith("保留中") for row in first_report))
            scope = json.loads(Path(wet["boundary_scope_manifest"]).read_text(encoding="utf-8"))
            boundary_path = Path(scope["artifacts"]["boundary"]["path"])
            boundary = pd.read_csv(boundary_path, sep="\t", dtype=str).fillna("")
            proposal = base / "proposal.tsv"
            pd.DataFrame([{"candidate_gene_id": gene, "decision": "adopt" if gene in genes[:2] else "defer",
                           "reason": "Recorded comparison retains these copies for experiment.",
                           "source_wet_sha256": sha256(boundary_path),
                           "executor_id": "fixture-executor", "executor_model": "fixture-model"}
                          for gene in boundary.candidate_gene_id]).to_csv(proposal, sep="\t", index=False)
            result = finish_wet_review(str(first), sha256(completion), str(base / "finished"), proposals_arg=str(proposal),
                executor_id="fixture-executor", executor_model="fixture-model")
            self.assertEqual("complete", result["review_history_state"])
            self.assertEqual("original_assessment_v1", result["review_history_basis"])
            final_rows = pd.read_csv(result["adopted_wet_candidates"], sep="\t", dtype=str).fillna("")
            final_report = candidate_table_rows(Path(result["report"]["path"]).read_text(encoding="utf-8"))
            self.assertEqual(genes[:2], list(final_rows.candidate_gene_id))
            self.assertEqual(len(final_rows), len(final_report))
            for supplied, displayed in ((first_rows, first_report), (final_rows, final_report)):
                for row, cells in zip(supplied.to_dict("records"), displayed):
                    self.assertEqual(6, len(cells))
                    self.assertIn(row["candidate_gene_id"], cells[0])
                    self.assertEqual(f"{row['chrom']}:{row['start']}-{row['end']}", cells[1])
                    self.assertEqual(row["relevance_tier"], cells[3])
                    self.assertIn("https://pubmed.ncbi.nlm.nih.gov/123456/", cells[4])
            self.assertIn("PsC_name", final_report[0][0])
            self.assertIn("MtA", final_report[0][2])
            self.assertIn("AtB", final_report[0][2])
            self.assertTrue(all("保留中" not in row[5] for row in final_report))
            self.assertEqual(sha256(result["adopted_wet_candidates"]), result["report"]["wet_candidates_sha256"])
            self.assertEqual("recommended", result["route_a_rethink"]["state"])
            self.assertEqual(2, result["route_a_rethink"]["final_wet_candidate_count"])
            self.assertEqual(result["route_a_rethink"], result["wet"]["route_a_rethink"])
            self.assertIn("Route A rethink:", Path(result["report"]["path"]).read_text(encoding="utf-8"))
            for path, content in {**frozen_bytes, **source_bytes}.items():
                self.assertEqual(content, path.read_bytes())

    def provisional_screen(self, wet, root):
        scope = json.loads(Path(wet["boundary_scope_manifest"]).read_text())
        original = pd.read_csv(scope["inventory"]["path"], sep="\t", dtype=str).fillna("")
        screen = root / "full-screen.tsv"
        # A unit crosses the future deferred/active boundary. Its original
        # members must survive verbatim; no projection can fake full coverage.
        pd.DataFrame([{"case_id": "case1", "unit_id": "all-interval",
                       "member_interval_genes": ";".join(original["candidate_gene_id"]),
                       "screen_state": "screened", "answer_sha256": "a" * 64}]).to_csv(screen, sep="\t", index=False)
        return screen, original

    def provisional_finish(self, first, completion, out, proposal, independent, screen=None, **kwargs):
        return finish_wet_review(
            str(first), sha256(completion), str(out), proposals_arg=str(proposal),
            independent_review_arg=str(independent), reviewer_id="fixture-reviewer",
            reviewer_model="fixture-review-model", executor_id="fixture-executor",
            executor_model="fixture-model", allow_provisional_wet=True,
            provisional_screen_ledger_arg=str(screen) if screen else None,
            provisional_screen_ledger_sha256=sha256(screen) if screen else None, **kwargs)

    def test_new_scope_does_not_reopen_identity_only_inventory_as_provisional_boundary(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            root = base / "first"
            wet = handoff(root, (0, 30, 4, 600, 366), max_rounds=1)
            screen, original = self.provisional_screen(wet, base)
            completion = _write_completion(root, wet)
            original_bytes = {path: path.read_bytes() for path in root.rglob("*") if path.is_file()}
            proposal, independent = _review(wet["boundary_scope_manifest"], base, adopt=2)
            result = self.provisional_finish(root, completion, base / "adopted", proposal, independent, screen)
            self.assertEqual("completed_evidence_scoped_selection", result["state"])
            self.assertEqual(2, result["adopted_candidate_count"])
            self.assertIsNone(result["wet"]["adaptive_synthesis_request"])
            self.assertIsNone(result["next_boundary_scope_manifest"])
            self.assertEqual("not_required", result["wet"]["second_stage_state"])
            adopted = pd.read_csv(result["adopted_wet_candidates"], sep="\t", dtype=str).fillna("")
            self.assertEqual(2, len(adopted))
            self.assertFalse(result["provisional_wet"]["eligible"])
            self.assertEqual("no_identity_only_boundary", result["provisional_wet"]["reason"])
            self.assertIn("evidence_scoped_v2", Path(result["report"]["path"]).read_text())
            for path, content in original_bytes.items():
                self.assertEqual(content, path.read_bytes())

    def test_new_scope_reviews_all_existing_tier_seeds_without_lower_identity_continuation(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            root = base / "first"
            wet = handoff(root, (0, 30, 40, 30, 5), max_rounds=2)
            screen, original = self.provisional_screen(wet, base)
            completion = _write_completion(root, wet)
            proposal, independent = _review(wet["boundary_scope_manifest"], base, adopt=2)
            first_root = base / "review1"
            first = self.provisional_finish(root, completion, first_root, proposal, independent, screen)
            self.assertEqual("completed_evidence_scoped_selection", first["state"])
            self.assertFalse(first["provisional_wet"]["eligible"])
            self.assertEqual("no_identity_only_boundary", first["provisional_wet"]["reason"])
            self.assertFalse(first["wet"]["boundary_scope_manifest"])
            self.assertIsNone(first["wet"]["adaptive_synthesis_request"])
            retained = pd.read_csv(first["adopted_wet_candidates"], sep="\t", dtype=str).fillna("")
            self.assertEqual(2, len(retained))

    def test_provisional_rejects_unbound_screen_and_unaccounted_prior_deferral(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            root = base / "first"
            wet = handoff(root, (0, 30, 4, 30, 5))
            screen, _ = self.provisional_screen(wet, base)
            completion = _write_completion(root, wet)
            proposal, independent = _review(wet["boundary_scope_manifest"], base, adopt=2)
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                finish_wet_review(str(root), sha256(completion), str(base / "bad-screen"),
                    proposals_arg=str(proposal), independent_review_arg=str(independent),
                    reviewer_id="fixture-reviewer", reviewer_model="fixture-review-model",
                    executor_id="fixture-executor", executor_model="fixture-model", allow_provisional_wet=True,
                    provisional_screen_ledger_arg=str(screen), provisional_screen_ledger_sha256="0" * 64)
            payload = json.loads(completion.read_text())
            payload["deferral_ledgers"] = [{"path": str(screen), "sha256": sha256(screen)}]
            completion.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError, "every prior deferral"):
                self.provisional_finish(root, completion, base / "bad-prior", proposal, independent, screen)

    def test_provisional_does_not_bypass_incomplete_review_or_identity_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            root = base / "first"
            wet = handoff(root, (0, 30, 4, 30, 5), max_rounds=2)
            screen, _ = self.provisional_screen(wet, base)
            completion = _write_completion(root, wet)
            proposal, independent = _review(wet["boundary_scope_manifest"], base, adopt=2)
            p = pd.read_csv(proposal, sep="\t")
            i = pd.read_csv(independent, sep="\t")
            p.loc[0, "decision"] = "not_evaluated"
            i.loc[0, "decision"] = "not_evaluated"
            p.to_csv(proposal, sep="\t", index=False)
            i["proposal_sha256"] = sha256(proposal)
            i.to_csv(independent, sep="\t", index=False)
            result = self.provisional_finish(root, completion, base / "incomplete", proposal, independent, screen)
            self.assertEqual("boundary_review_incomplete_retained", result["state"])
            self.assertFalse(result["provisional_wet"]["eligible"])
            self.assertIsNone(result["wet"]["adaptive_synthesis_request"])

    def test_resealed_request_cannot_change_boundary_roster(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "first"
            wet = handoff(root, (3, 8, 35, 0, 0))
            request_path = Path(wet["adaptive_synthesis_request"])
            request = json.loads(request_path.read_text())
            request["candidate_gene_ids"] = ["foreign-member"]
            request_path.write_text(json.dumps(request))
            completion = _write_completion(root, wet)
            proposal, independent = _review(wet["boundary_scope_manifest"], Path(tmp), adopt=9)
            with self.assertRaisesRegex(ValueError, "request roster"):
                finish_wet_review(str(root), sha256(completion), str(Path(tmp) / "review"),
                    proposals_arg=str(proposal), independent_review_arg=str(independent),
                    executor_id="fixture-executor", executor_model="fixture-model",
                    reviewer_id="fixture-reviewer", reviewer_model="fixture-review-model")

    def test_cli_finish_review_uses_independent_flags_and_writes_hashed_wet(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "first"
            wet = handoff(root, (3, 8, 35, 12, 5), max_rounds=2)
            completion = _write_completion(root, wet)
            proposal, independent = _review(wet["boundary_scope_manifest"], Path(tmp), adopt=8)
            first_out = Path(tmp) / "cli-review1"
            script = Path(__file__).parents[1] / "scripts" / "claim_review_session.py"
            command = [sys.executable, str(script), "finish-wet-review",
                       "--handoff", str(root), "--completion-sha256", sha256(completion),
                       "--out", str(first_out), "--proposals", str(proposal),
                       "--reviewer-id", "fixture-reviewer", "--reviewer-model", "fixture-review-model",
                       "--independent-review", str(independent), "--executor-id", "fixture-executor",
                       "--executor-model", "fixture-model"]
            completed = subprocess.run(command, cwd=Path(__file__).parents[1], check=True,
                                       capture_output=True, text=True)
            first = json.loads(completed.stdout)
            self.assertEqual("completed_evidence_scoped_selection", first["state"])
            self.assertEqual(19, first["adopted_candidate_count"])
            adopted_path = Path(first["adopted_wet_candidates"])
            self.assertEqual(first["adopted_wet_sha256"], hashlib.sha256(adopted_path.read_bytes()).hexdigest())
            report = Path(first["report"]["path"]).read_text(encoding="utf-8")
            self.assertIn("Wet SHA-256", report)
            self.assertIn("review decision artifact", report)
            self.assertIn("next action", report)

    def test_header_only_review_retains_every_boundary_member_as_incomplete(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "first"
            wet = handoff(root, (3, 8, 35, 12, 5), max_rounds=2)
            completion = _write_completion(root, wet)
            proposal = Path(tmp) / "empty.tsv"
            pd.DataFrame(columns=["candidate_gene_id", "decision", "reason",
                                  "source_wet_sha256", "executor_id", "executor_model"]).to_csv(
                proposal, sep="\t", index=False)
            result = finish_wet_review(
                str(root), sha256(completion), str(Path(tmp) / "empty-out"),
                proposals_arg=str(proposal), reviewer_id="fixture-reviewer",
                reviewer_model="fixture-review-model", executor_id="fixture-executor",
                executor_model="fixture-model")
            self.assertEqual("boundary_review_incomplete_retained", result["state"])
            self.assertEqual(46, result["adopted_candidate_count"])
            self.assertEqual("incomplete", result["review_history_state"])

    def test_explicit_not_evaluated_agreement_remains_incomplete(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "first"
            wet = handoff(root, (3, 8, 35, 12, 5), max_rounds=2)
            completion = _write_completion(root, wet)
            proposal, independent = _review(wet["boundary_scope_manifest"], Path(tmp), adopt=8)
            proposal_frame = pd.read_csv(proposal, sep="\t")
            independent_frame = pd.read_csv(independent, sep="\t")
            proposal_frame.loc[0, "decision"] = "not_evaluated"
            independent_frame.loc[0, "decision"] = "not_evaluated"
            proposal_frame.to_csv(proposal, sep="\t", index=False)
            independent_frame["proposal_sha256"] = hashlib.sha256(proposal.read_bytes()).hexdigest()
            independent_frame.to_csv(independent, sep="\t", index=False)
            result = finish_wet_review(
                str(root), sha256(completion), str(Path(tmp) / "explicit-out"),
                proposals_arg=str(proposal), reviewer_id="fixture-reviewer",
                reviewer_model="fixture-review-model", independent_review_arg=str(independent),
                executor_id="fixture-executor", executor_model="fixture-model")
            self.assertEqual("boundary_review_incomplete_retained", result["state"])
            report = Path(result["report"]["path"]).read_text(encoding="utf-8")
            self.assertIn("review_history_state: `incomplete`", report)
            self.assertIn("decision counts:", report)
            self.assertEqual("incomplete", result["review_history_state"])
    def test_finish_comparison_scope_contains_all_lower_eligible_seeds(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "first"
            wet = handoff(root, (3, 8, 35, 12, 5), max_rounds=2)
            completion = _write_completion(root, wet)
            proposal, independent = _review(Path(wet["boundary_scope_manifest"]).__str__(), Path(tmp), adopt=8)
            first = finish_wet_review(
                str(root), sha256(completion), str(Path(tmp) / "review1"),
                proposals_arg=str(proposal), reviewer_id="fixture-reviewer",
                reviewer_model="fixture-review-model", independent_review_arg=str(independent),
                executor_id="fixture-executor", executor_model="fixture-model")
            self.assertEqual("completed_evidence_scoped_selection", first["state"])
            self.assertEqual(19, first["boundary_selection"]["output_candidate_count"])
            self.assertIsNone(first["wet"]["adaptive_synthesis_request"])
            self.assertFalse(first["wet"]["boundary_scope_manifest"])
            ledger = first["boundary_selection"]["deferral_ledger"]["path"]
            self.assertEqual(27, len(pd.read_csv(ledger, sep="\t")))

    def test_incomplete_history_is_preserved_and_scope_tamper_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "first"
            wet = handoff(root, (3, 8, 35, 12, 5), max_rounds=2)
            completion = _write_completion(root, wet)
            proposal, independent = _review(wet["boundary_scope_manifest"], Path(tmp), adopt=8, disagreement=True)
            first = finish_wet_review(
                str(root), sha256(completion), str(Path(tmp) / "review1"),
                proposals_arg=str(proposal), reviewer_id="fixture-reviewer",
                reviewer_model="fixture-review-model", independent_review_arg=str(independent),
                executor_id="fixture-executor", executor_model="fixture-model")
            self.assertEqual("boundary_review_incomplete_retained", first["state"])
            self.assertEqual("incomplete", first["wet"]["review_history_state"])

            tampered = Path(wet["boundary_scope_manifest"])
            payload = json.loads(tampered.read_text(encoding="utf-8"))
            inventory = Path(payload["inventory"]["path"])
            inventory.write_text(inventory.read_text(encoding="utf-8") + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "artifact mismatch|hash mismatch"):
                finish_wet_review(
                    str(root), sha256(completion), str(Path(tmp) / "tampered"),
                    proposals_arg=str(proposal), reviewer_id="fixture-reviewer",
                    reviewer_model="fixture-review-model", independent_review_arg=str(independent),
                    executor_id="fixture-executor", executor_model="fixture-model")

    def test_agentic_boundary_merge_feeds_two_round_wet_continuation(self):
        from scripts.agentic_boundary import merge as merge_boundary
        from scripts.agentic_boundary import prepare as prepare_boundary

        def merged_round(request_path, prepared, response_dir, adopt):
            manifest = prepare_boundary(request_path, str(prepared), max_candidates=5)
            proposal_responses = []
            review_responses = []
            remaining_adopt = adopt
            for shard in manifest["shards"]:
                records = json.loads((prepared / shard["file"]).read_text(encoding="utf-8"))["candidate_input"]
                record_responses = []
                for index, row in enumerate(records):
                    decision = "adopt" if remaining_adopt > 0 else "defer"
                    if remaining_adopt > 0:
                        remaining_adopt -= 1
                    record_responses.append({
                        "candidate_gene_id": row["candidate_gene_id"],
                        "decision": decision,
                        "reason": f"fixture agentic comparison {index}",
                    })
                proposal_responses.append({
                    "case_id": manifest["case_id"],
                    "source_wet_sha256": manifest["source_wet_sha256"],
                    "shard_id": shard["shard_id"],
                    "status": "completed",
                    "executor_id": "fixture-executor",
                    "executor_model": "fixture-model",
                    "responses": record_responses,
                })
            for payload in proposal_responses:
                path = response_dir / f"proposal-{payload['shard_id']}.json"
                path.write_text(json.dumps(payload), encoding="utf-8")
            proposal = response_dir / "proposal.tsv"
            merge_boundary(
                str(prepared), [str(response_dir / f"proposal-{p['shard_id']}.json") for p in proposal_responses],
                str(proposal), role="proposal", executor_id="fixture-executor",
                executor_model="fixture-model",
            )
            proposal_sha = hashlib.sha256(proposal.read_bytes()).hexdigest()
            decisions = pd.read_csv(proposal, sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            for payload in proposal_responses:
                review_responses.append({
                    "case_id": manifest["case_id"],
                    "source_wet_sha256": manifest["source_wet_sha256"],
                    "proposal_sha256": proposal_sha,
                    "shard_id": payload["shard_id"],
                    "status": "completed",
                    "reviewer_id": "fixture-reviewer",
                    "reviewer_model": "fixture-review-model",
                    "responses": [
                        {
                            "candidate_gene_id": row["candidate_gene_id"],
                            "decision": decisions.loc[row["candidate_gene_id"], "decision"],
                            "reason": f"fixture independent comparison {row['candidate_gene_id']}",
                        }
                        for row in payload["responses"]
                    ],
                })
            for payload in review_responses:
                path = response_dir / f"review-{payload['shard_id']}.json"
                path.write_text(json.dumps(payload), encoding="utf-8")
            independent = response_dir / "independent.tsv"
            merge_boundary(
                str(prepared), [str(response_dir / f"review-{p['shard_id']}.json") for p in review_responses],
                str(independent), role="review", reviewer_id="fixture-reviewer",
                reviewer_model="fixture-review-model", proposal_path=str(proposal),
            )
            return proposal, independent

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "first"
            wet = handoff(root, (3, 8, 35, 12, 5), max_rounds=2)
            completion = _write_completion(root, wet, {
                "route_b_execution": {"state": "complete", "completed_count": 3, "not_evaluated_count": 0,
                                      "response_state_counts": {"completed_claims": 1, "completed_no_claim": 2}},
                "route_c_execution": {"state": "not_evaluated_no_input"},
            })
            first_dir = Path(tmp) / "agentic-first"
            first_dir.mkdir()
            proposal, independent = merged_round(
                wet["adaptive_synthesis_request"], Path(tmp) / "agentic-boundary-1", first_dir, 8,
            )
            first = finish_wet_review(
                str(root), sha256(completion), str(Path(tmp) / "agentic-review-1"),
                proposals_arg=str(proposal), independent_review_arg=str(independent),
                reviewer_id="fixture-reviewer", reviewer_model="fixture-review-model",
                executor_id="fixture-executor", executor_model="fixture-model",
            )
            self.assertEqual("completed_evidence_scoped_selection", first["state"])
            self.assertEqual(19, first["adopted_candidate_count"])
            report = Path(first["report"]["path"]).read_text(encoding="utf-8")
            self.assertIn("execution context", report)
            self.assertIn("not_evaluated_no_input", report)
            self.assertIn("Route B response states: completed_claims=1, completed_no_claim=2", report)

    def test_incomplete_boundary_state_wins_over_round_budget(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "first"
            wet = handoff(root, (3, 8, 35, 12, 5), max_rounds=1)
            completion = _write_completion(root, wet)
            proposal, independent = _review(wet["boundary_scope_manifest"], Path(tmp), adopt=8)
            proposal_frame = pd.read_csv(proposal, sep="\t")
            independent_frame = pd.read_csv(independent, sep="\t")
            proposal_frame.loc[0, "decision"] = "not_evaluated"
            independent_frame.loc[0, "decision"] = "not_evaluated"
            proposal_frame.to_csv(proposal, sep="\t", index=False)
            independent_frame["proposal_sha256"] = hashlib.sha256(proposal.read_bytes()).hexdigest()
            independent_frame.to_csv(independent, sep="\t", index=False)
            result = finish_wet_review(
                str(root), sha256(completion), str(Path(tmp) / "incomplete-budget-out"),
                proposals_arg=str(proposal), reviewer_id="fixture-reviewer",
                reviewer_model="fixture-review-model", independent_review_arg=str(independent),
                executor_id="fixture-executor", executor_model="fixture-model",
            )
            self.assertEqual("incomplete", result["review_history_state"])
            self.assertEqual("boundary_review_incomplete_retained", result["state"])

    def test_round_bound_retains_candidates_and_returns_exhausted(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "first"
            wet = handoff(root, (3, 8, 35, 12, 5), max_rounds=1)
            completion = _write_completion(root, wet)
            proposal, independent = _review(wet["boundary_scope_manifest"], Path(tmp), adopt=8)
            out = Path(tmp) / "review1"
            result = finish_wet_review(
                str(root), sha256(completion), str(out), proposals_arg=str(proposal),
                reviewer_id="fixture-reviewer", reviewer_model="fixture-review-model",
                independent_review_arg=str(independent), executor_id="fixture-executor",
                executor_model="fixture-model")
            self.assertIsNone(result["wet"]["adaptive_synthesis_request"])
            self.assertEqual(19, result["boundary_selection"]["output_candidate_count"])
            with self.assertRaisesRegex(ValueError, "no pending boundary review"):
                finish_wet_review(
                    str(out), sha256(out / "completion.json"), str(Path(tmp) / "exhausted"),
                    proposals_arg=str(proposal), reviewer_id="fixture-reviewer",
                    reviewer_model="fixture-review-model")

    def test_incomplete_evidence_scoped_review_retries_from_immutable_first_stage(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            root = base / "first"
            wet = handoff(root, (3, 8, 35, 12, 5), max_rounds=2)
            completion = _write_completion(root, wet)
            (base / "incomplete").mkdir()
            proposal, independent = _review(wet["boundary_scope_manifest"], base / "incomplete", adopt=9)
            proposal_frame = pd.read_csv(proposal, sep="\t")
            independent_frame = pd.read_csv(independent, sep="\t")
            proposal_frame.loc[0, "decision"] = "not_evaluated"
            independent_frame.loc[0, "decision"] = "not_evaluated"
            proposal_frame.to_csv(proposal, sep="\t", index=False)
            independent_frame["proposal_sha256"] = sha256(proposal)
            independent_frame.to_csv(independent, sep="\t", index=False)
            partial_out = base / "partial"
            partial = finish_wet_review(
                str(root), sha256(completion), str(partial_out), proposals_arg=str(proposal),
                independent_review_arg=str(independent), reviewer_id="fixture-reviewer",
                reviewer_model="fixture-review-model", executor_id="fixture-executor",
                executor_model="fixture-model")
            self.assertEqual("boundary_review_incomplete_retained", partial["state"])
            self.assertEqual("review_incomplete_retry_original_scope", partial["wet"]["second_stage_state"])
            self.assertEqual(wet["boundary_scope_manifest"], partial["resume_scope_manifest"]["path"])
            self.assertEqual(str(completion.resolve()), partial["resume_from_completion"]["path"])
            with self.assertRaisesRegex(ValueError, "no pending boundary review"):
                finish_wet_review(
                    str(partial_out), sha256(partial_out / "completion.json"), str(base / "bad-resume"),
                    proposals_arg=str(proposal), reviewer_id="fixture-reviewer",
                    reviewer_model="fixture-review-model")
            (base / "retry").mkdir()
            retry_proposal, retry_independent = _review(
                wet["boundary_scope_manifest"], base / "retry", adopt=9)
            retried = finish_wet_review(
                str(root), sha256(completion), str(base / "complete-retry"),
                proposals_arg=str(retry_proposal), independent_review_arg=str(retry_independent),
                reviewer_id="fixture-reviewer", reviewer_model="fixture-review-model",
                executor_id="fixture-executor", executor_model="fixture-model")
            self.assertEqual("completed_evidence_scoped_selection", retried["state"])
            self.assertEqual(20, retried["adopted_candidate_count"])


if __name__ == "__main__":
    unittest.main()
