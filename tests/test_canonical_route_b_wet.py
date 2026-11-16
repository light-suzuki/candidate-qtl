"""Canonical Route-B semantic -> supported family/copy -> Wet regressions."""

from __future__ import annotations

import json
import hashlib
import tempfile
import threading
import unittest
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path

import pandas as pd

from candidate_qtl.canonical_route_b_wet import run_canonical_route_b_wet_extension
from candidate_qtl.cli import main
from tests._issue664_helpers import manifest_for


def _write_b0(root: Path) -> Path:
    out = root / "b01_route_b"
    out.mkdir(parents=True, exist_ok=True)
    path = out / "interval_gene_atlas.tsv"
    pd.DataFrame([
        {"gene_id": "G_A", "family_id": "FAM_A"},
        {"gene_id": "G_B", "family_id": "FAM_B"},
    ]).to_csv(path, sep="\t", index=False, lineterminator="\n")
    return path


def _write_route_a_bridge(root: Path) -> Path:
    out = root / "z01b_family_candidate_bridge"
    out.mkdir(parents=True, exist_ok=True)
    path = out / "family_candidate_bridge.tsv"
    pd.DataFrame([{
        "route": "a",
        "gene_id": "G_A",
        "family_id": "FAM_A",
        "hypothesis_id": "A1",
        # Must use the canonical bridge vocabulary now that #514 feeds this
        # through family_to_copy's fail-closed bridge validation.
        "mapping_tier": "tier3_sequence_match_unique",
        "ambiguity_status": "unique",
        "review_status": "supported",
    }]).to_csv(path, sep="\t", index=False, lineterminator="\n")
    manifest_for(path, route="a", case_id="case1")
    return path


def _write_legacy_b1_decoy(root: Path) -> None:
    out = root / "b02_route_b1_overlay"
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{
        "hypothesis_id": "legacy1",
        "gene_family": "DECOY",
        "member_interval_genes": "G_LEGACY_DECOY",
    }]).to_csv(out / "family_process_hypotheses.tsv", sep="\t", index=False, lineterminator="\n")


def _write_legacy_family_copy_decoy(root: Path) -> None:
    compare = root / "z01_triple_route_compare"
    compare.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{
        "hypothesis_id": "legacy-old",
        "gene_family": "F_OLD",
        "family_id": "F_OLD",
        "member_interval_genes": "G_OLD",
    }]).to_csv(compare / "family_to_copy_input.tsv", sep="\t", index=False, lineterminator="\n")
    pd.DataFrame([{
        "family_id": "F_OLD",
        "four_state": "a_and_b_supported",
        "convergence_status": "strict_convergent",
        "source_dependence_status": "source_independent",
    }]).to_csv(compare / "triple_route_hypothesis_comparison.tsv", sep="\t", index=False, lineterminator="\n")

    copy = root / "z02_family_to_copy"
    copy.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{
        "family_id": "F_OLD",
        "candidate_gene_id": "G_OLD",
        "disambiguation_status": "preferred_copy",
        "rank_contribution": "none",
    }]).to_csv(copy / "family_copy_discrimination.tsv", sep="\t", index=False, lineterminator="\n")
    pd.DataFrame([{
        "family_id": "F_OLD",
        "candidate_gene_id": "G_OLD",
        "disambiguation_status": "preferred_copy",
        "leading_axes": "orthology_class",
        "subfamily": "",
        "orthology_class": "one_to_one",
        "first_wet_priority": "first",
        "rank_contribution": "none",
        "notes": "legacy decoy",
    }]).to_csv(copy / "family_copy_shortlist.tsv", sep="\t", index=False, lineterminator="\n")


def _empty_claim_responses(prompt_path: str, destination: Path) -> Path:
    prompts = [json.loads(line) for line in Path(prompt_path).read_text(encoding="utf-8").splitlines() if line.strip()]
    payload = {
        "candidates": [
            {"prompt_id": prompt["prompt_id"], "parse_status": "json", "model_json": {"claims": []}}
            for prompt in prompts
        ]
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    return destination


def _write_llm_manifest(bundle: str, prompts: str, responses: Path, destination: Path) -> Path:
    def digest(path: str | Path) -> str:
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()

    destination.write_text(json.dumps({
        "provider": "test-provider",
        "model": "test-model",
        "bundle_sha256": digest(bundle),
        "prompt_sha256": digest(prompts),
        "response_sha256": digest(responses),
    }, indent=2) + "\n", encoding="utf-8")
    return destination


class CanonicalRouteBWetTests(unittest.TestCase):
    def test_direct_c_reaches_family_copy_and_wet_while_b_review_is_pending(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write_b0(root)
            members = root / "route_c_direct_members.tsv"
            pd.DataFrame([
                {
                    "candidate_gene_id": "G_B", "record_id": "direct-same", "source_id": "PMID:1",
                    "family_ids": "FAM_B", "orthogroup_ids": "", "species_distance": "S1_same_species",
                    "evidence_directness": "D0_causal_perturbation_or_rescue",
                    "projection_state": "inside_current_qtl", "projection_basis": "exact_focal_id",
                    "support_status": "verified_direct_support", "rank_contribution": "none",
                },
                {
                    "candidate_gene_id": "G_A", "record_id": "direct-model", "source_id": "PMID:2",
                    "family_ids": "FAM_A", "orthogroup_ids": "", "species_distance": "S5_model_or_other_plant",
                    "evidence_directness": "D0_causal_perturbation_or_rescue",
                    "projection_state": "inside_current_qtl", "projection_basis": "orthology_reference_gene",
                    "support_status": "reference_phenotype_review", "rank_contribution": "none",
                },
            ]).to_csv(members, sep="\t", index=False, lineterminator="\n")
            manifest_for(members, route="c", case_id="direct-c")

            result = run_canonical_route_b_wet_extension(
                str(root), phenotype="fasciated meristem", route_c_members_arg=str(members),
                case_id="direct-c",
            )

            self.assertIn(result["state"], {"wet_reviewable", "wet_reviewable_degraded_coverage"})
            wet = pd.read_csv(result["wet_candidates"], sep="\t", dtype=str).fillna("")
            same = wet[wet["candidate_gene_id"] == "G_B"].iloc[0]
            self.assertEqual("tier1_verified_support", same["support_tier"])
            self.assertEqual("verified_route_c_support", same["membership_basis"])
            self.assertEqual(["G_B"], wet["candidate_gene_id"].tolist())
            # This hand-written C table has no sealed projection/binding ancestry,
            # so G_A is inventory context rather than an eligible C review seed.
            handoff_manifest = json.loads(Path(result["wet_handoff_manifest"]).read_text(encoding="utf-8"))
            wet_manifest = json.loads(Path(handoff_manifest["final_wet_candidate_manifest"]).read_text(encoding="utf-8"))
            trace = pd.read_csv(wet_manifest["candidate_trace"], sep="\t", dtype=str).fillna("")
            model = trace.set_index("candidate_gene_id").loc["G_A"]
            self.assertEqual("not_selected_no_upstream_candidate", model["selection_disposition"])
            families = pd.read_csv(result["canonical_family_recovery"], sep="\t", dtype=str).fillna("")
            self.assertEqual({"G_B"}, set(families["member_interval_genes"]))
            self.assertTrue(Path(result["family_copy"]).exists())
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual("route_b_not_evaluated_no_semantic_response", manifest["provider_state"])

    def test_canonical_extension_carries_b0_source_function_into_prompt(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write_b0(root)
            identity_hits = root / "b01_route_b" / "route_b_identity_hits.tsv"
            pd.DataFrame([{
                "focal_gene_id": "G_A", "hit_accession": "REF_ALT",
                "hit_species": "Medicago truncatula", "hit_gene_symbol": "ALT1",
                "hit_protein_name": "receptor-like kinase", "hit_metadata_state": "complete",
                "hit_metadata_records_json": "[]", "raw_reference_gene_id": "REF_ALT",
                "percent_identity": "81.0", "query_coverage": "0.99", "evalue": "1e-30",
                "bitscore": "200", "provider": "local_blastp:test", "database": "test-db",
                "identity_policy": "test_v1", "rank_contribution": "none",
            }]).to_csv(identity_hits, sep="\t", index=False, lineterminator="\n")
            lookup_dir = root / "b01_route_b" / "route_b0"
            lookup_dir.mkdir(parents=True, exist_ok=True)
            pd.DataFrame([{
                "gene_id": "G_A", "reference_gene": "REF_A",
                "reference_function": "controls meristem size", "reference_pathway": "meristem maintenance",
                "reference_mutant_phenotype": "fasciation", "reference_source_id": "PMID:123",
                "lookup_status": "source_backed",
            }]).to_csv(lookup_dir / "interval_reverse_lookup.tsv", sep="\t", index=False, lineterminator="\n")

            result = run_canonical_route_b_wet_extension(
                str(root), phenotype="fasciated meristem", case_id="case1"
            )
            self.assertEqual("semantic_review_required", result["state"])
            prompts = Path(result["prompts"]).read_text(encoding="utf-8")
            self.assertIn("controls meristem size", prompts)
            self.assertIn("PMID:123", prompts)
            self.assertIn("REF_ALT", prompts)
            self.assertIn("ALT1", prompts)
            bundle_manifest = json.loads(
                (root / "b04_canonical_semantic" / "01_bundle" / "route_b_semantic_bundle_manifest.json")
                .read_text(encoding="utf-8")
            )
            self.assertEqual(1, bundle_manifest["stable_source_count"])
            self.assertEqual(
                str(identity_hits),
                bundle_manifest["inputs"]["route_b_identity_hits"]["path"],
            )

    def test_no_semantic_response_is_review_required_not_negative(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write_b0(root)
            result = run_canonical_route_b_wet_extension(
                str(root), phenotype="fasciated meristem", case_id="case1"
            )
            self.assertEqual("semantic_review_required", result["state"])
            self.assertFalse(result["legacy_b1_consumed"])
            self.assertFalse(result["legacy_triple_route_family_context_consumed"])
            self.assertTrue(Path(result["prompts"]).exists())
            self.assertIsNone(result["wet_candidates"])
            self.assertEqual(
                root / "b04_canonical_semantic" / "02_llm_run",
                Path(result["expected_semantic_run_dir"]),
            )
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual("not_evaluated_no_semantic_response", manifest["provider_state"])
            self.assertEqual("diagnostic_only", manifest["legacy_interval_reverse_lookup_role"])

    def test_conventional_response_path_preserves_a_and_ignores_all_legacy_b_context(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write_b0(root)
            _write_route_a_bridge(root)
            _write_legacy_b1_decoy(root)
            _write_legacy_family_copy_decoy(root)

            prepared = run_canonical_route_b_wet_extension(
                str(root), phenotype="fasciated meristem", case_id="case1"
            )
            responses = _empty_claim_responses(
                prepared["prompts"],
                root / "b04_canonical_semantic" / "02_llm_run" / "llm_candidate_responses.json",
            )
            _write_llm_manifest(
                prepared["bundle"], prepared["prompts"], responses,
                responses.parent / "llm_run_manifest.json",
            )
            self.assertTrue(responses.exists())
            result = run_canonical_route_b_wet_extension(
                str(root),
                phenotype="fasciated meristem",
                case_id="case1",
            )

            self.assertEqual("wet_reviewable_degraded_coverage", result["state"])
            self.assertEqual("route_b_declared_but_unverified", result["input_coverage_state"])
            self.assertFalse(result["legacy_b1_consumed"])
            self.assertFalse(result["legacy_triple_route_family_context_consumed"])
            wet = pd.read_csv(result["wet_candidates"], sep="\t", dtype=str).fillna("")
            # Supported membership is still exactly G_A. Issue #570 fills the
            # unused review slots from the frozen interval inventory, so other
            # interval genes may appear -- as labelled review coverage, never
            # as support, and never a gene from a legacy artifact.
            supported = wet[wet["membership_basis"] == "verified_route_a_support"]
            self.assertEqual({"G_A"}, set(supported["candidate_gene_id"]))
            self.assertEqual({"tier1_verified_support"}, set(supported["support_tier"]))
            filled = wet[wet["membership_basis"] == "interval_review_coverage"]
            self.assertNotIn("tier1_verified_support", set(filled["support_tier"]))
            self.assertNotIn("G_LEGACY_DECOY", set(wet["candidate_gene_id"]))
            self.assertNotIn("G_OLD", set(wet["candidate_gene_id"]))

            b_frozen = pd.read_csv(result["route_b_verified_candidates"], sep="\t", dtype=str).fillna("")
            self.assertTrue(b_frozen.empty)
            canonical_families = pd.read_csv(
                result["canonical_family_recovery"], sep="\t", dtype=str
            ).fillna("")
            self.assertEqual({"FAM_A"}, set(canonical_families["family_id"]))
            self.assertNotIn("F_OLD", set(canonical_families["family_id"]))

            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual("conventional_path", manifest["semantic_response_discovery"])
            self.assertEqual("diagnostic_only", manifest["legacy_interval_reverse_lookup_role"])
            self.assertFalse(manifest["legacy_b1_consumed"])
            self.assertFalse(manifest["legacy_triple_route_family_context_consumed"])
            self.assertTrue(manifest["supported_union_identity_verified"])
            self.assertEqual(
                manifest["pre_family_supported_union_sha256"],
                manifest["handoff_supported_union_sha256"],
            )
            self.assertEqual("supported_member_union", manifest["family_recovery_source"])

    def test_response_rerun_reuses_frozen_bundle_and_prompts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write_b0(root)
            _write_route_a_bridge(root)
            prepared = run_canonical_route_b_wet_extension(
                str(root), phenotype="fasciated meristem", case_id="case1"
            )
            responses = _empty_claim_responses(
                prepared["prompts"],
                root / "b04_canonical_semantic" / "02_llm_run" / "llm_candidate_responses.json",
            )
            _write_llm_manifest(
                prepared["bundle"], prepared["prompts"], responses,
                responses.parent / "llm_run_manifest.json",
            )
            bundle_before = Path(prepared["bundle"]).read_bytes()
            prompts_before = Path(prepared["prompts"]).read_bytes()

            result = run_canonical_route_b_wet_extension(
                str(root), phenotype="changed phenotype", case_id="case1"
            )

            self.assertEqual("wet_reviewable_degraded_coverage", result["state"])
            self.assertEqual(bundle_before, Path(prepared["bundle"]).read_bytes())
            self.assertEqual(prompts_before, Path(prepared["prompts"]).read_bytes())

    def test_standard_llm_run_manifest_allows_canonical_response_import(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write_b0(root)
            _write_route_a_bridge(root)
            prepared = run_canonical_route_b_wet_extension(
                str(root), phenotype="fasciated meristem", case_id="case1"
            )

            class Handler(BaseHTTPRequestHandler):
                def do_POST(self):  # noqa: N802
                    length = int(self.headers.get("Content-Length", "0"))
                    json.loads(self.rfile.read(length))
                    body = json.dumps({
                        "choices": [{"message": {"content": json.dumps({"claims": []})}}]
                    }).encode()
                    self.send_response(200)
                    self.send_header("Content-Type", "application/json")
                    self.send_header("Content-Length", str(len(body)))
                    self.end_headers()
                    self.wfile.write(body)

                def log_message(self, *_args):
                    return

            server = HTTPServer(("127.0.0.1", 0), Handler)
            thread = threading.Thread(target=server.serve_forever, daemon=True)
            thread.start()
            try:
                self.assertEqual(0, main([
                    "llm-run",
                    "--prompts", prepared["prompts"],
                    "--out-dir", str(root / "b04_canonical_semantic" / "02_llm_run"),
                    "--provider", "openai-compatible",
                    "--model", "fixture-model",
                    "--endpoint", f"http://127.0.0.1:{server.server_port}/v1/chat/completions",
                ]))
            finally:
                server.shutdown()
                thread.join(timeout=2)
                server.server_close()

            llm_manifest_path = root / "b04_canonical_semantic" / "02_llm_run" / "llm_run_manifest.json"
            llm_manifest = json.loads(llm_manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(
                hashlib.sha256(Path(prepared["bundle"]).read_bytes()).hexdigest(),
                llm_manifest["bundle_sha256"],
            )
            self.assertEqual("fixture-model", llm_manifest["model"])
            result = run_canonical_route_b_wet_extension(
                str(root), phenotype="changed phenotype", case_id="case1"
            )
            self.assertEqual("wet_reviewable_degraded_coverage", result["state"])
            extension_manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual(
                llm_manifest["response_sha256"],
                extension_manifest["semantic_artifact_hashes"]["response_sha256"],
            )

    def test_frozen_artifact_hash_mismatch_fails_closed(self) -> None:
        for artifact in ("bundle", "prompts", "responses"):
            with self.subTest(artifact=artifact), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                _write_b0(root)
                prepared = run_canonical_route_b_wet_extension(
                    str(root), phenotype="fasciated meristem", case_id="case1"
                )
                responses = _empty_claim_responses(
                    prepared["prompts"],
                    root / "b04_canonical_semantic" / "02_llm_run" / "llm_candidate_responses.json",
                )
                manifest_path = _write_llm_manifest(
                    prepared["bundle"], prepared["prompts"], responses,
                    responses.parent / "llm_run_manifest.json",
                )
                target = {
                    "bundle": Path(prepared["bundle"]),
                    "prompts": Path(prepared["prompts"]),
                    "responses": responses,
                }[artifact]
                target.write_bytes(target.read_bytes() + b"\n")

                with self.assertRaisesRegex(ValueError, "hash mismatch"):
                    run_canonical_route_b_wet_extension(
                        str(root), phenotype="fasciated meristem", case_id="case1"
                    )
                self.assertTrue(manifest_path.exists())

    def test_provider_model_only_legacy_manifest_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write_b0(root)
            prepared = run_canonical_route_b_wet_extension(
                str(root), phenotype="fasciated meristem", case_id="case1"
            )
            responses = _empty_claim_responses(
                prepared["prompts"],
                root / "b04_canonical_semantic" / "02_llm_run" / "llm_candidate_responses.json",
            )
            (responses.parent / "llm_run_manifest.json").write_text(
                json.dumps({"provider": "legacy", "model": "legacy-model"}) + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "missing frozen artifact hashes"):
                run_canonical_route_b_wet_extension(
                    str(root), phenotype="fasciated meristem", case_id="case1"
                )

    def test_manifest_without_response_does_not_regenerate_first_pass_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write_b0(root)
            prepared = run_canonical_route_b_wet_extension(
                str(root), phenotype="fasciated meristem", case_id="case1"
            )
            bundle_before = Path(prepared["bundle"]).read_bytes()
            manifest = root / "b04_canonical_semantic" / "02_llm_run" / "llm_run_manifest.json"
            manifest.parent.mkdir(parents=True, exist_ok=True)
            manifest.write_text(json.dumps({"provider": "legacy", "model": "legacy-model"}) + "\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "without a semantic response"):
                run_canonical_route_b_wet_extension(
                    str(root), phenotype="changed phenotype", case_id="case1"
                )
            self.assertEqual(bundle_before, Path(prepared["bundle"]).read_bytes())

    def test_missing_b0_is_not_runnable_instead_of_a_biological_miss(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            root.mkdir(exist_ok=True)
            result = run_canonical_route_b_wet_extension(
                str(root), phenotype="fasciated meristem", case_id="case1"
            )
            self.assertEqual("not_runnable_upstream_b0_missing", result["state"])
            self.assertFalse(result["legacy_triple_route_family_context_consumed"])
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual("not_evaluated_upstream_b0_missing", manifest["provider_state"])


if __name__ == "__main__":
    unittest.main()
