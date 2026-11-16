"""Issue #782: fixed-answer consumer acceptance, never model/causal accuracy.

All production consumers are real. Only external agent answers and reference
resources are synthetic; no semantic-review receipt or verified claim is made.
The ordinary FASTA/GFF/BLAST entry regression remains owned by #758.
"""
from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.canonical_route_b_wet import run_canonical_route_b_wet_extension
from candidate_qtl.claim_review_session import finish_wet_review
from candidate_qtl.route_b_semantic_bundle import build_route_b_semantic_bundle
from candidate_qtl.route_b_semantic_import import import_route_b_semantic_proposals
from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from scripts.agentic_boundary import merge as merge_boundary, prepare as prepare_boundary
from scripts.agentic_route_b import merge
from tests._issue664_helpers import manifest_for
from tests.test_boundary_session_integration import _review, _write_completion
from tests.test_wet_relevance_selection import handoff


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: object) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return path


def write_table(path: Path, rows: list[dict]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, lineterminator="\n")
    return path


def frame(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


class FixedAnswerAcceptanceTests(unittest.TestCase):
    def test_declared_frozen_material_is_delivered_verbatim_or_rejected_before_prompt(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atlas = write_table(root / "atlas.tsv", [
                {"gene_id": "complete", "family_id": "F1", "orthology_reference_gene": "ref-full"},
                {"gene_id": "truncated", "family_id": "F2", "orthology_reference_gene": "ref-short"},
            ])
            # Distinct head and tail catch the historical short-field failure;
            # all text is generic synthetic context, not an answer-selected hint.
            long_function = "START enzyme activity " + "recorded experimental context " * 700 + "END observed substrate change"
            corpus = write_table(root / "acquired.tsv", [
                {"reference_gene_id": "ref-full", "reference_species": "Species one", "source_id": "PMID:1", "function_description": long_function},
                {"reference_gene_id": "ref-short", "reference_species": "Species two", "source_id": "PMID:2", "function_description": "[truncated]"},
                {"reference_gene_id": "unrelated", "reference_species": "Species three", "source_id": "PMID:3", "function_description": "WRONG_RESOURCE enzyme context"},
            ])
            sidecar = root / "route_b_preparation_manifest.json"
            descriptor = {"command": "route-b-prepare", "frozen_output_fingerprints": {"interval_gene_atlas": {"sha256": digest(atlas)}},
                "providers": {"reference_function": {"available": True, "path": str(corpus), "sha256": digest(corpus)}},
                "outputs": {"reference_function": str(corpus)}}
            write_json(sidecar, descriptor)
            result = build_route_b_semantic_bundle(str(atlas), str(root / "valid"), phenotype="seed color", case_id="fixture")
            bundle = json.loads(Path(result["bundle"]).read_text())
            units = {u["member_interval_genes"][0]: u for u in bundle["units"]}
            prompts = [json.loads(line) for line in Path(result["prompts"]).read_text().splitlines()]
            self.assertEqual(1, len(prompts))
            payload = json.loads(prompts[0]["user"].split("\n\n", 1)[1])
            acquired = frame(corpus).set_index("reference_gene_id", drop=False).loc["ref-full"]
            delivered = payload["reference_material_cards"][0]
            for key in ("reference_gene_id", "reference_species", "source_id", "function_description"):
                self.assertEqual(acquired[key], delivered[key], f"first divergence: acquired resource -> actual prompt field {key}")
            self.assertEqual(long_function, units["complete"]["reference_material_cards"][0]["function_description"])
            self.assertNotIn("WRONG_RESOURCE", prompts[0]["user"])
            self.assertEqual("[truncated]", units["truncated"]["reference_material_cards"][0]["function_description"])
            self.assertEqual("not_requested_no_reviewable_material", units["truncated"]["request_state"])
            manifest = json.loads(Path(result["manifest"]).read_text())
            self.assertEqual(digest(corpus), manifest["inputs"]["reference_materials"]["corpus_sha256"])
            responses = write_json(root / "answers.json", {"candidates": [
                {"prompt_id": prompts[0]["prompt_id"], "parse_status": "json", "model_json": {"claims": []}}]})
            imported = import_route_b_semantic_proposals(result["bundle"], result["prompts"], str(responses), str(root / "import"))
            status = frame(imported["status"]).set_index("unit_id")
            self.assertEqual("no_claims", status.loc[units["complete"]["unit_id"], "status"])
            self.assertEqual("not_evaluated_no_reviewable_material", status.loc[units["truncated"]["unit_id"], "status"])
            frozen_prompt = Path(result["prompts"]).read_bytes()
            mutations = [
                ("wrong-resource", {**descriptor, "outputs": {"reference_function": str(root / "other.tsv")}}, "declared reference corpus path"),
                ("missing-declared-resource", {**descriptor, "outputs": {}}, "declared reference corpus path"),
            ]
            for label, changed, error in mutations:
                write_json(sidecar, changed)
                with self.assertRaisesRegex(ValueError, error):
                    build_route_b_semantic_bundle(str(atlas), str(root / label), phenotype="seed color", case_id="fixture")
                self.assertFalse((root / label / "route_b_semantic_prompts.jsonl").exists())
            write_json(sidecar, descriptor)
            corpus.write_text(corpus.read_text().replace(long_function, "START enzyme activity [truncated]"))
            with self.assertRaisesRegex(ValueError, "corpus missing or fingerprint mismatch"):
                build_route_b_semantic_bundle(str(atlas), str(root / "truncated-after-freeze"), phenotype="seed color", case_id="fixture")
            self.assertFalse((root / "truncated-after-freeze/route_b_semantic_prompts.jsonl").exists())
            self.assertEqual(frozen_prompt, Path(result["prompts"]).read_bytes())

    def test_available_material_reaches_payload_and_missing_material_is_not_model_omission(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atlas = write_table(root / "atlas.tsv", [
                {"gene_id": "material", "family_id": "F1", "protein_name": "meristem receptor"},
                {"gene_id": "missing", "family_id": "F2", "protein_name": ""},
            ])
            # An available resource that was not explicitly supplied is not
            # silently discovered or falsely reported as delivered.
            unused = write_table(root / "available-but-unsupplied.tsv", [
                {"gene_id": "missing", "function_description": "UNDELIVERED enzyme activity"},
            ])
            bundle = build_route_b_semantic_bundle(str(atlas), str(root / "bundle"), phenotype="fasciation", case_id="fixture")
            units = json.loads(Path(bundle["bundle"]).read_text())["units"]
            by_gene = {u["member_interval_genes"][0]: u for u in units}
            prompts = [json.loads(line) for line in Path(bundle["prompts"]).read_text().splitlines()]
            self.assertEqual(1, len(prompts))
            self.assertIn("meristem receptor", prompts[0]["user"])
            self.assertNotIn("UNDELIVERED", prompts[0]["user"])
            self.assertTrue(unused.is_file())
            self.assertEqual("requested", by_gene["material"]["request_state"])
            self.assertEqual("not_requested_no_reviewable_material", by_gene["missing"]["request_state"])
            responses = write_json(root / "responses.json", {"candidates": [
                {"prompt_id": prompts[0]["prompt_id"], "parse_status": "json", "model_json": {"claims": []}},
            ]})
            imported = import_route_b_semantic_proposals(bundle["bundle"], bundle["prompts"], str(responses), str(root / "import"))
            status = frame(imported["status"]).set_index("unit_id")
            self.assertEqual("no_claims", status.loc[by_gene["material"]["unit_id"], "status"])
            self.assertEqual("not_evaluated_no_reviewable_material", status.loc[by_gene["missing"]["unit_id"], "status"])
            # Material sufficiency and biological interpretation are not proven
            # by these observables; no causal/model-failure attribution is made.

    def test_multi_claim_source_reference_and_copies_reach_actual_wet_without_support_promotion(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            b0 = root / "b01_route_b"
            atlas = write_table(b0 / "interval_gene_atlas.tsv", [
                {"gene_id": gene, "family_id": "F", "orthology_reference_gene": "refA;refB"}
                for gene in ("copy1", "copy2")
            ])
            corpus = write_table(b0 / "reference.tsv", [
                {"reference_gene_id": ref, "reference_species": "Species one", "function_description": function,
                 "source_id": source}
                for ref, function, source in (("refA", "enzyme activity", "PMID:1"), ("refB", "receptor activity", "PMID:2"))
            ])
            write_json(b0 / "route_b_preparation_manifest.json", {
                "command": "route-b-prepare", "frozen_output_fingerprints": {"interval_gene_atlas": {"sha256": digest(atlas)}},
                "providers": {"reference_function": {"available": True, "path": str(corpus), "sha256": digest(corpus)}},
                "outputs": {"reference_function": str(corpus)},
            })
            prepared = run_canonical_route_b_wet_extension(str(root), phenotype="seed color", case_id="fixture")
            prompt_rows = [json.loads(line) for line in Path(prepared["prompts"]).read_text().splitlines()]
            self.assertEqual(1, len(prompt_rows))
            for material in ("enzyme activity", "receptor activity", "PMID:1", "PMID:2", "copy1", "copy2"):
                self.assertIn(material, prompt_rows[0]["user"])
            valid = {"claim_type": "mechanistic_process_bridge", "biological_function": "pigment metabolism",
                     "mechanism": "supplies seed pigment", "phenotype_relation": "seed color changes",
                     "reference_species": "Species one", "evidence_scope": "family"}
            claims = [
                {**valid, "reference_gene": "refA", "source_ids": ["PMID:1", "PMID:2"]},
                {**valid, "reference_gene": "refB", "source_ids": ["PMID:2"]},
                {**valid, "reference_gene": "foreign", "source_ids": ["PMID:999"]},
            ]
            raw = json.dumps({"claims": claims})
            attempt = write_json(root / "attempt.json", {"case_id": "fixture", "executor": "TEST_ONLY", "model": "fixed-answer",
                "responses": [{"prompt_id": prompt_rows[0]["prompt_id"], "raw_text": raw}]})
            merged = root / "merged.json"
            execution = merge(prepared["prompts"], [str(attempt)], str(merged))
            manifest = merged.with_name("merged_manifest.json")
            execution["bundle_sha256"] = digest(Path(prepared["bundle"]))
            write_json(manifest, execution)
            frozen = {p: Path(p).read_bytes() for p in (prepared["bundle"], prepared["prompts"], str(merged))}
            result = run_canonical_route_b_wet_extension(str(root), phenotype="seed color", case_id="fixture",
                responses_arg=str(merged), llm_manifest_arg=str(manifest))
            proposals = frame(root / "b04_canonical_semantic/03_import/semantic_claim_proposals.tsv")
            self.assertEqual(4, len(proposals))
            self.assertEqual({"refA", "refB", "foreign"}, set(proposals.reference_gene))
            self.assertEqual("unsupported", proposals.loc[proposals.reference_gene.eq("foreign"), "verification_status"].iloc[0])
            review = frame(result["route_b_review_candidates"])
            self.assertTrue(review.empty)
            route_manifest = json.loads(Path(result["manifest"]).read_text())
            ledger_meta = json.loads(Path(route_manifest["route_b_review_candidates_manifest"]).read_text())
            ledger = frame(ledger_meta["eligibility_ledger"]["path"])
            for gene in ("copy1", "copy2"):
                rows = ledger[ledger.candidate_gene_id.eq(gene)]
                self.assertEqual(set(proposals.claim_id), set(rows.claim_id))
            wet = frame(result["wet_candidates"])
            self.assertTrue(wet.empty)
            self.assertTrue(frame(result["route_b_verified_candidates"]).empty)
            inventory = frame(Path(result["wet_candidates"]).parent / "wet_relevance_inventory.tsv")
            self.assertEqual({"copy1", "copy2"}, set(inventory.candidate_gene_id))
            handoff_meta = json.loads(Path(result["wet_handoff_manifest"]).read_text())
            report = Path(handoff_meta["report"]["path"]).read_text()
            self.assertIn(digest(Path(result["wet_candidates"])), report)
            self.assertIn("candidate_count: 0", report)
            self.assertEqual(execution["response_sha256"], digest(merged))
            for path, content in frozen.items():
                self.assertEqual(content, Path(path).read_bytes())

    def test_equivalent_a_b_aliases_have_same_tier_and_membership_but_weak_claim_is_distinct(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atlas = write_table(root / "atlas.tsv", [{"gene_id": "gene", "family_id": "F"}])
            supported = write_table(root / "supported.tsv", [])
            # Empty supported artifact still has the producer's declared schema.
            supported.write_text("candidate_gene_id\troute\tsupport_status\n")
            supported_meta = manifest_for(supported, route="b", case_id="fixture")
            outputs = []
            for lane, mechanism_field, claim_type in (("a", "mechanism_statement", ""), ("b", "mechanism", ""), ("b", "mechanistic_link", "mechanistic_review_seed")):
                index = len(outputs)
                review = write_table(root / f"review{index}.tsv", [{"candidate_gene_id": "gene", "family_id": "F",
                    "reference_gene": "ref", "phenotype_relation": "seed color changes", mechanism_field: "supplies seed pigment",
                    "claim_type": claim_type, "route": lane, "support_status": f"unverified_route_{lane}_review", "rank_contribution": "none"}])
                seal = write_json(root / f"review{index}_manifest.json", {"contract": "route-review-candidates-1.0.0",
                    "stage": f"route_{lane}_review_candidate_freeze", "case_id": "fixture", "inputs": {}, "output": str(review), "artifact_sha256": digest(review)})
                result = wet_candidate_handoff(str(root / f"out{index}"), route_b_candidates_arg=str(supported),
                    route_b_candidates_manifest_arg=str(supported_meta), review_candidates_arg=str(review),
                    review_candidates_manifest_arg=str(seal), interval_atlas_arg=str(atlas), case_id="fixture")
                outputs.append(frame(result["wet_candidates"]).iloc[0])
            self.assertEqual(outputs[0].relevance_tier, outputs[1].relevance_tier)
            self.assertEqual(outputs[0].membership_basis, outputs[1].membership_basis)
            self.assertTrue(outputs[0].relevance_tier.startswith("T2"))
            self.assertTrue(outputs[2].relevance_tier.startswith("T3"))
            self.assertNotEqual(outputs[1].relevance_tier, outputs[2].relevance_tier)

    def test_19_20_21_25_26_are_whole_tiers_and_only_26_requests_boundary(self):
        for count in (19, 20, 21, 25, 26):
            with self.subTest(count=count), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                result = handoff(root, (0, count, 0, 0, 0))
                wet = frame(result["wet_candidates"])
                self.assertEqual(count, wet.candidate_gene_id.nunique())
                meta = json.loads(Path(result["manifest"]).read_text())
                self.assertEqual(count, meta["final_wet_candidate_count"])
                scope = json.loads(Path(result["boundary_scope_manifest"]).read_text())
                self.assertEqual(count if count > 25 else 0, scope["counts"]["boundary"])
                report = Path(result["report"]["path"]).read_text()
                self.assertIn(f"candidate_count: {count}", report)
                self.assertIn(digest(Path(result["wet_candidates"])), report)
                if count == 26:
                    self.assertEqual("requested", result["second_stage_state"])
                    self.assertEqual("boundary_review_requested", result["handoff_state"])
                    self.assertEqual(26, scope["counts"]["boundary"])
                else:
                    self.assertEqual("skipped_at_or_below_25", result["second_stage_state"])

    def test_boundary_shards_inline_assigned_material_and_expose_sealed_whole_scope(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            wet = handoff(root / "first", (3, 26, 0, 0, 0))
            scope = json.loads(Path(wet["boundary_scope_manifest"]).read_text())
            higher = frame(scope["artifacts"]["frozen"]["path"]).candidate_gene_id.tolist()
            boundary = frame(scope["artifacts"]["boundary"]["path"]).candidate_gene_id.tolist()
            all_ids = set(higher + boundary)
            context = write_json(root / "comparison.json", {"case_id": "case1", "phenotype": "seed color",
                "candidate_material": {cid: {"function": f"recorded context for {cid}", "support_state": "unverified"} for cid in sorted(all_ids)}})
            out = root / "shards"
            prepared = prepare_boundary(wet["adaptive_synthesis_request"], str(out), max_candidates=8,
                comparison_context_path=str(context))
            self.assertGreater(len(prepared["shards"]), 1)
            assigned = []
            for descriptor in prepared["shards"]:
                shard_path = out / descriptor["file"]
                self.assertEqual(descriptor["sha256"], digest(shard_path))
                shard = json.loads(shard_path.read_text())
                own = shard["candidate_gene_ids"]
                assigned.extend(own)
                self.assertEqual(set(own + higher), set(shard["candidate_comparison_material"]))
                self.assertEqual("read_only", shard["higher_tier_context"]["state"])
                # Other shards are delivered as a descriptor, not embedded in
                # every bounded payload. Read its exact bytes as an executor
                # would; this establishes access, not that an agent compared it.
                material_path = Path(shard["comparison_context"]["path"])
                self.assertEqual(shard["comparison_context"]["sha256"], digest(material_path))
                full = json.loads(material_path.read_text())
                self.assertEqual(all_ids, set(full["candidate_material"]))
                self.assertTrue(set(boundary) - set(own))
                for cid in own + higher:
                    self.assertEqual(full["candidate_material"][cid], shard["candidate_comparison_material"][cid])
            self.assertEqual(sorted(boundary), sorted(assigned))
            response = write_json(root / "response.json", {"case_id": "case1", "source_wet_sha256": prepared["source_wet_sha256"],
                "comparison_context_sha256": digest(context), "shard_id": prepared["shards"][0]["shard_id"],
                "executor_id": "TEST_ONLY", "executor_model": "fixed-answer", "responses": [
                    {"candidate_gene_id": cid, "decision": "not_evaluated", "reason": "no semantic comparison executed"}
                    for cid in prepared["shards"][0]["candidate_gene_ids"]]})
            merged = merge_boundary(str(out), [str(response)], str(root / "partial.tsv"), role="proposal",
                executor_id="TEST_ONLY", executor_model="fixed-answer")
            self.assertEqual(26, merged["state_counts"]["not_evaluated"])
            context.write_text(context.read_text() + "\n")
            with self.assertRaisesRegex(ValueError, "comparison context SHA mismatch"):
                merge_boundary(str(out), [str(response)], str(root / "changed.tsv"), role="proposal",
                    executor_id="TEST_ONLY", executor_model="fixed-answer")

    def test_thousand_gene_inventory_reaches_adopted_25_and_reviewed_900_is_not_compression_success(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = root / "first"
            initial = handoff(first, (0, 900, 0, 100, 0), max_rounds=2)
            initial_ids = set(frame(initial["wet_candidates"]).candidate_gene_id)
            self.assertEqual(900, len(initial_ids))
            trace = frame(first / "wet/wet_candidate_trace.tsv")
            all_ids = set(trace.candidate_gene_id)
            self.assertEqual(1000, len(all_ids))
            lower_ids = set(trace.loc[trace.selection_disposition.eq("not_selected_no_upstream_candidate"), "candidate_gene_id"])
            self.assertEqual(100, len(lower_ids))
            self.assertEqual("boundary_review_requested", initial["handoff_state"])
            completion = _write_completion(first, initial)
            frozen = {path: path.read_bytes() for path in first.rglob("*") if path.is_file()}
            observed = []
            for adopt in (900, 25):
                response_root = root / f"answers-{adopt}"
                response_root.mkdir()
                # These fixed adopt/defer answers exercise the existing machine
                # consumer. They are not biological rankings or real reviews.
                proposal, independent = _review(initial["boundary_scope_manifest"], response_root, adopt=adopt)
                result = finish_wet_review(str(first), digest(completion), str(root / f"adopted-{adopt}"),
                    proposals_arg=str(proposal), independent_review_arg=str(independent),
                    executor_id="fixture-executor", executor_model="fixture-model",
                    reviewer_id="fixture-reviewer", reviewer_model="fixture-review-model")
                wet = frame(result["adopted_wet_candidates"])
                wet_ids = set(wet.candidate_gene_id)
                self.assertEqual("completed_evidence_scoped_selection", result["state"])
                self.assertEqual("complete", result["review_history_state"])
                self.assertEqual(adopt, len(wet_ids))
                self.assertEqual(adopt, result["adopted_candidate_count"])
                self.assertEqual(digest(Path(result["adopted_wet_candidates"])), result["adopted_wet_sha256"])
                deferred = frame(result["boundary_selection"]["deferral_ledger"]["path"])
                deferred_ids = set(deferred.candidate_gene_id)
                self.assertEqual(900 - adopt, len(deferred_ids))
                self.assertFalse(wet_ids & deferred_ids)
                self.assertFalse((wet_ids | deferred_ids) & lower_ids)
                self.assertEqual(all_ids, wet_ids | deferred_ids | lower_ids)
                report = Path(result["report"]["path"]).read_text()
                self.assertIn(f"candidate_count: {adopt}", report)
                self.assertIn(result["adopted_wet_sha256"], report)
                self.assertEqual(adopt, result["report"]["candidate_count"])
                # Process completion and the user's manageable Wet capacity
                # are separate. Policy never invents a Top-N fallback.
                within_requested_capacity = 20 <= len(wet_ids) <= 25
                self.assertEqual(adopt == 25, within_requested_capacity)
                observed.append({"interval_gene_count": 1000, "actual_adopted_wet_count": len(wet_ids),
                    "deferred_count": len(deferred_ids), "unreviewed_lower_count": len(lower_ids),
                    "stage": result["state"], "requested_compression_met": within_requested_capacity})
            self.assertFalse(observed[0]["requested_compression_met"])
            self.assertTrue(observed[1]["requested_compression_met"])
            for path, content in frozen.items():
                self.assertEqual(content, path.read_bytes())
            print("Issue #782 synthetic compression acceptance: " + json.dumps(observed, sort_keys=True))


if __name__ == "__main__":
    unittest.main()
