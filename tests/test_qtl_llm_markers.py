from __future__ import annotations

import csv
import hashlib
import json
import os
import tempfile
import threading
import unittest
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.exon_assays import collect_sanger_markers
from candidate_qtl.llm import export_llm_prompts, run_llm_prompts
from candidate_qtl.qtl_mapping import qtl2_plan


class QtlLlmMarkerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_qtl2_bundle_uses_f2_contract_and_mbp_physical_map(self) -> None:
        genotypes = self.root / "geno.tsv"
        phenotypes = self.root / "pheno.tsv"
        markers = self.root / "markers.tsv"
        marker_ids = [f"m{i}" for i in range(1, 7)]
        genotype_rows = []
        phenotype_rows = []
        codes = ["AA", "AB", "BB"]
        for i in range(1, 61):
            calls = [codes[(i // (2 + j) + j) % 3] for j in range(6)]
            genotype_rows.append(f"F2_{i}\t" + "\t".join(calls))
            dosage = {"AA": 0, "AB": 1, "BB": 2}[calls[2]]
            phenotype_rows.append(f"F2_{i}\t{10 + dosage * 5 + (i % 5) / 5}\t{1 + i % 2}")
        genotypes.write_text("individual_id\t" + "\t".join(marker_ids) + "\n" + "\n".join(genotype_rows) + "\n", encoding="utf-8")
        phenotypes.write_text(
            "individual_id\ttrait\tblock\n" + "\n".join(phenotype_rows) + "\n",
            encoding="utf-8",
        )
        markers.write_text(
            "marker_id\tchrom\tpos_bp\tpos_cm\n" + "\n".join(f"m{i}\t1\t{i * 1000000}\t{(i - 1) * 5}" for i in range(1, 7)) + "\n",
            encoding="utf-8",
        )
        execute = os.environ.get("CANDIDATE_QTL_RUN_RQTL2") == "1"
        result = qtl2_plan(
            str(genotypes), str(phenotypes), str(markers), str(self.root / "qtl"), "trait", ["block"],
            n_perm=10 if execute else 1000, execute=execute,
        )
        self.assertEqual("provided_pos_cm", result["map_source"])
        self.assertIn("covar: covar.csv", (self.root / "qtl/qtl2_control.yaml").read_text(encoding="utf-8"))
        with (self.root / "qtl/pmap.csv").open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual("1.0", rows[0]["pos"])
        if execute:
            self.assertTrue((self.root / "qtl/lod_profile.tsv").exists())
            self.assertTrue((self.root / "qtl/crossover_counts.tsv").exists())

    def test_prompt_export_and_external_api_guard(self) -> None:
        packet = self.root / "packet.json"
        packet.write_text(
            json.dumps({"species": "Pisum sativum", "trait": "flowering", "candidates": [{"gene_id": "GeneA", "required_output_schema": {}}]}),
            encoding="utf-8",
        )
        prompts = self.root / "prompts.jsonl"
        result = export_llm_prompts(str(packet), str(prompts))
        self.assertEqual(1, result["prompt_count"])
        record = json.loads(prompts.read_text(encoding="utf-8"))
        self.assertIn("Never invent a citation", record["system"])
        with self.assertRaisesRegex(ValueError, "external endpoint blocked"):
            run_llm_prompts(str(prompts), str(self.root / "llm"), "openai-compatible", "test", "https://example.com/v1/chat/completions")

        class Handler(BaseHTTPRequestHandler):
            def do_POST(self):  # noqa: N802
                length = int(self.headers.get("Content-Length", "0"))
                json.loads(self.rfile.read(length))
                body = json.dumps({"choices": [{"message": {"content": json.dumps({"gene_id": "GeneA", "claims": []})}}]}).encode()
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
            run = run_llm_prompts(
                str(prompts), str(self.root / "local_llm"), "openai-compatible", "test",
                f"http://127.0.0.1:{server.server_port}/v1/chat/completions",
            )
        finally:
            server.shutdown()
            thread.join(timeout=2)
            server.server_close()
        self.assertEqual(1, run["prompt_count"])
        normalized = json.loads(Path(run["normalized_responses"]).read_text(encoding="utf-8"))
        self.assertEqual("json", normalized["candidates"][0]["parse_status"])
        generic_manifest = json.loads(Path(run["manifest"]).read_text(encoding="utf-8"))
        self.assertNotIn("bundle_sha256", generic_manifest)

    def test_canonical_route_b_llm_run_freezes_artifact_hashes(self) -> None:
        prompts = self.root / "route_b_semantic_prompts.jsonl"
        prompts.write_text(
            json.dumps({"prompt_id": "p1", "gene_id": "GeneA", "system": "system", "user": "user"}) + "\n",
            encoding="utf-8",
        )
        bundle = self.root / "route_b_semantic_bundle.json"
        bundle.write_text("{\"schema_version\": 1}\n", encoding="utf-8")

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
            run = run_llm_prompts(
                str(prompts), str(self.root / "canonical_llm"), "openai-compatible", "test",
                f"http://127.0.0.1:{server.server_port}/v1/chat/completions",
            )
        finally:
            server.shutdown()
            thread.join(timeout=2)
            server.server_close()

        manifest = json.loads(Path(run["manifest"]).read_text(encoding="utf-8"))
        self.assertTrue(manifest["canonical_route_b_semantic"])
        self.assertEqual(
            manifest["bundle_sha256"],
            hashlib.sha256(bundle.read_bytes()).hexdigest(),
        )
        self.assertEqual(
            manifest["prompt_sha256"],
            hashlib.sha256(prompts.read_bytes()).hexdigest(),
        )
        normalized_path = Path(run["normalized_responses"])
        self.assertEqual(
            manifest["response_sha256"],
            hashlib.sha256(normalized_path.read_bytes()).hexdigest(),
        )

    def test_canonical_route_b_llm_run_requires_sibling_bundle(self) -> None:
        prompts = self.root / "route_b_semantic_prompts.jsonl"
        prompts.write_text(
            json.dumps({"prompt_id": "p1", "gene_id": "GeneA", "system": "system", "user": "user"}) + "\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(ValueError, "require sibling route_b_semantic_bundle.json"):
            run_llm_prompts(
                str(prompts), str(self.root / "canonical_llm"), "openai-compatible", "test",
                "http://127.0.0.1:9/v1/chat/completions",
            )

    def test_offline_agent_file_shapes_normalize_and_handoff(self) -> None:
        """The real CLI seam accepts agent wrappers and reaches semantic import."""
        for shape in ("array", "responses", "candidates", "jsonl-single"):
            with self.subTest(shape=shape):
                root = self.root / shape
                root.mkdir()
                prompts = root / "route_b_semantic_prompts.jsonl"
                prompt = {
                    "prompt_id": "p1", "route_b_unit_id": "u1", "gene_id": "GeneA",
                    "system": "system", "user": "user",
                }
                prompts.write_text(json.dumps(prompt) + "\n", encoding="utf-8")
                bundle = root / "route_b_semantic_bundle.json"
                bundle.write_text(json.dumps({
                    "schema_version": "1.1", "route": "interval_first",
                    "units": [{
                        "unit_id": "u1", "identity": {}, "member_interval_genes": [],
                        "stable_source_ids": [], "genes": [],
                    }],
                }) + "\n", encoding="utf-8")
                raw_text = '{"claims":[]}'
                record = {
                    "prompt_id": "p1", "route_b_unit_id": "u1", "gene_id": "GeneA",
                    "raw_text": raw_text, "model_json": {"claims": []},
                    "parse_status": "ok",
                }
                agent_file = root / "agent.json"
                if shape == "array":
                    payload = [record]
                elif shape == "jsonl-single":
                    payload = record
                else:
                    payload = {shape: [record]}
                agent_file.write_text(json.dumps(payload) + "\n", encoding="utf-8")
                output = root / "b04_canonical_semantic" / "02_llm_run"

                self.assertEqual(0, main([
                    "llm-run", "--prompts", str(prompts), "--out-dir", str(output),
                    "--provider", "codex-agent", "--model", "fixture-agent",
                    "--agent-file", str(agent_file), "--json",
                ]))
                normalized_path = output / "llm_candidate_responses.json"
                normalized = json.loads(normalized_path.read_text(encoding="utf-8"))
                self.assertEqual(["p1"], [row["prompt_id"] for row in normalized["candidates"]])
                self.assertEqual("json", normalized["candidates"][0]["parse_status"])
                self.assertEqual(raw_text, normalized["candidates"][0]["raw_text"])
                manifest = json.loads((output / "llm_run_manifest.json").read_text(encoding="utf-8"))
                self.assertEqual("codex-agent", manifest["provider"])
                self.assertEqual(hashlib.sha256(agent_file.read_bytes()).hexdigest(), manifest["raw_agent_artifact_sha256"])
                self.assertEqual(agent_file.read_bytes(), (output / "llm_agent_raw_artifact.json").read_bytes())
                self.assertEqual(hashlib.sha256(prompts.read_bytes()).hexdigest(), manifest["prompt_sha256"])
                self.assertEqual(hashlib.sha256(bundle.read_bytes()).hexdigest(), manifest["bundle_sha256"])
                self.assertEqual(hashlib.sha256(normalized_path.read_bytes()).hexdigest(), manifest["response_sha256"])

                imported = root / "import"
                self.assertEqual(0, main([
                    "route-b-semantic-import", "--bundle", str(bundle),
                    "--prompts", str(prompts), "--responses", str(normalized_path),
                    "--llm-manifest", str(output / "llm_run_manifest.json"),
                    "--out-dir", str(imported), "--json",
                ]))
                status = pd.read_csv(imported / "semantic_claim_proposal_status.tsv", sep="\t", dtype=str).fillna("")
                self.assertEqual(["no_claims"], status["status"].tolist())

    def test_offline_agent_file_fails_closed_on_mismatch_and_roster(self) -> None:
        prompts = self.root / "prompts.jsonl"
        prompts.write_text(json.dumps({
            "prompt_id": "p1", "gene_id": "GeneA", "system": "system", "user": "user",
        }) + "\n", encoding="utf-8")
        mismatch = self.root / "mismatch.json"
        mismatch.write_text(json.dumps({"responses": [{
            "prompt_id": "p1", "raw_text": '{"claims":[]}', "model_json": {"claims": [{"x": 1}]},
        }]}), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "inconsistent JSON fields"):
            run_llm_prompts(
                str(prompts), str(self.root / "mismatch-out"), "codex-agent", "fixture-agent",
                agent_file=str(mismatch),
            )

        roster = self.root / "roster.json"
        roster.write_text(json.dumps({"candidates": []}), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "roster must exactly match"):
            run_llm_prompts(
                str(prompts), str(self.root / "roster-out"), "codex-agent", "fixture-agent",
                agent_file=str(roster),
            )

        malformed_claims = self.root / "malformed-claims.json"
        malformed_claims.write_text(json.dumps({"responses": [{
            "prompt_id": "p1", "model_json": {"claims": "not-a-list"},
        }]}), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "requires a claims list"):
            run_llm_prompts(
                str(prompts), str(self.root / "malformed-claims-out"), "codex-agent", "fixture-agent",
                agent_file=str(malformed_claims),
            )

    def test_collect_sanger_markers(self) -> None:
        result_json = self.root / "result.json"
        result_json.write_text(
            json.dumps({"mode": "design", "pairs": [{"index": 1, "left_start": 5, "right_start": 150, "forward": "ACGT", "reverse": "TGCA", "product_size": 146, "tm_f": 60.1, "tm_r": 60.0, "specificity": {"rank": "A", "specific_all_db": True, "gel_distinguishable": True, "total_off_target": 0}}]}),
            encoding="utf-8",
        )
        plan = self.root / "exon_sequencing_plan.tsv"
        plan.write_text(
            "gene_id\trepresentative_transcript\ttarget_interval\tlocal_target_start\tlocal_target_end\ttemplate_fasta\toutput_json\n"
            f"GeneA\tTxA\tchr1:1-100\t10\t100\t{self.root / 'template.fa'}\t{result_json}\n",
            encoding="utf-8",
        )
        output = self.root / "candidate_sanger_markers.tsv"
        result = collect_sanger_markers(str(plan), str(output))
        self.assertEqual(1, result["marker_count"])
        self.assertIn("ACGT", output.read_text(encoding="utf-8"))

        result_json.write_text(
            json.dumps({"mode": "design", "pairs": [{"index": 1, "left_start": 5, "right_start": 80, "forward": "ACGT", "reverse": "TGCA", "specificity": {"rank": "A"}}]}),
            encoding="utf-8",
        )
        rejected = collect_sanger_markers(str(plan), str(self.root / "partial.tsv"))
        self.assertEqual(0, rejected["marker_count"])
        self.assertFalse(rejected["complete"])
        self.assertEqual(["chr1:1-100"], rejected["failed_windows"])


if __name__ == "__main__":
    unittest.main()
