"""Discovery orchestration preserves the old generation and original raw answers."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pandas as pd

from candidate_qtl.route_b_semantic_bundle import build_route_b_semantic_bundle
from candidate_qtl.route_b_semantic_import import import_route_b_semantic_proposals
from scripts.agentic_route_b import _sha, merge, preflight
from scripts.agentic_route_b_extend import extend


class RouteBExtendTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        atlas = self.root / "atlas.tsv"
        pd.DataFrame([{"gene_id": f"G{i}", "family_id": f"F{i}",
                       "family_index_subject_endpoints": f"P{i}1111",
                       "local_annotation_description": "developmental protein"}
                      for i in range(1, 4)]).to_csv(atlas, sep="\t", index=False)
        self.bundle = build_route_b_semantic_bundle(str(atlas), str(self.root / "parent"),
                                                   phenotype="altered development", case_id="case-x")
        self.prompts = [json.loads(line) for line in Path(self.bundle["prompts"]).read_text().splitlines()]
        self.parent_sha = _sha(Path(self.bundle["bundle"]))
        self.raw = self.root / "raw.json"
        self.write(self.raw, {"case_id": "case-x", "executor": "original-agent", "model": "original-model",
                              "responses": [{"prompt_id": p["prompt_id"], "route_b_unit_id": p["route_b_unit_id"],
                                             "raw_text": '{"claims":[]}'} for p in self.prompts]})
        self.checkpoint = self.root / "checkpoint.json"
        merge(self.bundle["prompts"], [str(self.raw)], str(self.checkpoint))
        hashes = self.root / "import-inputs.json"
        self.write(hashes, {"bundle_sha256": self.parent_sha, "prompt_sha256": _sha(Path(self.bundle["prompts"])),
                            "response_sha256": _sha(self.checkpoint)})
        imported = import_route_b_semantic_proposals(self.bundle["bundle"], self.bundle["prompts"],
                                                     str(self.checkpoint), str(self.root / "imported"), llm_manifest_arg=str(hashes))
        self.import_manifest = imported["manifest"]
        self.discovery = self.root / "discovery.json"
        parent = json.loads(Path(self.bundle["bundle"]).read_text())
        self.changed_unit = parent["units"][0]["unit_id"]
        self.write(self.discovery, {"case_id": "case-x", "parent_bundle_sha256": self.parent_sha,
                                    "discovered_sources": [{"unit_id": self.changed_unit,
                                        "reference_gene": parent["units"][0]["genes"][0]["family_index_subject_endpoints"],
                                        "reference_species": "Model plant", "reported_gene": "MODEL1",
                                        "reported_species": "Model plant", "source_id": "PMID:1234567",
                                        "source_title": "Synthetic development paper", "source_url": "https://example.org/paper",
                                        "retrieved_at": "2026-10-28", "source_text": "A reported developmental function.",
                                        "reference_function": "developmental signaling"}]})

    def write(self, path, data):
        path.write_text(json.dumps(data), encoding="utf-8")

    def test_alias_arrays_are_transported_without_rewriting_original_discoveries(self):
        data = json.loads(self.discovery.read_text())
        data["discovered_sources"][0]["source_identifier_aliases"] = ["PMCID:PMC1234567", "DOI:10.1234/example"]
        self.write(self.discovery, data)
        original = self.discovery.read_bytes()
        result = self.extend()
        self.assertEqual(original, self.discovery.read_bytes())
        self.assertEqual((2, 1), (result["reusable_count"], result["new_execution_required_count"]))
        normalized = json.loads(Path(result["discoveries"]).read_text())
        self.assertEqual("PMCID:PMC1234567;DOI:10.1234/example", normalized["discovered_sources"][0]["source_identifier_aliases"])
        self.assertEqual(1, len(result["transport_normalizations"]))
        self.assertEqual(_sha(self.discovery), result["transport_normalizations"][0]["sha256"])

    def test_alias_array_rejects_non_strings(self):
        data = json.loads(self.discovery.read_text())
        data["discovered_sources"][0]["source_identifier_aliases"] = [123]
        self.write(self.discovery, data)
        with self.assertRaisesRegex(ValueError, "nonempty strings"):
            self.extend()

    def extend(self, name="next", audited=True, discoveries=None):
        return extend(self.bundle["bundle"], self.parent_sha, self.bundle["prompts"],
                      discoveries or [str(self.discovery)], str(self.checkpoint), str(self.root / name),
                      import_manifest=self.import_manifest if audited else "")

    def test_changed_only_dispatch_and_reused_raw_reaches_real_merge_and_import(self):
        before = {p: _sha(p) for p in self.root.rglob("*") if p.is_file()}
        result = self.extend()
        self.assertEqual((2, 1, 0), (result["reusable_count"], result["new_execution_required_count"], result["new_model_execution_count"]))
        self.assertEqual(before, {p: _sha(p) for p in before})
        dispatch = json.loads(Path(result["redispatch"]).read_text())
        self.assertEqual(1, len(dispatch["reading_shards"]))
        shard = dispatch["reading_shards"][0]
        self.assertEqual([self.changed_unit], shard["unit_ids"])
        pid = shard["prompt_ids"][0]
        retry = self.root / "retry.json"
        self.write(retry, {"case_id": "case-x", "executor": "new-agent", "model": "new-model",
                           "responses": [{"prompt_id": pid, "raw_text": '{"claims":[]}',
                                          "read_receipts": shard["read_receipts_by_prompt"][pid]}]})
        responses = [row["file"] for row in result["reusable_responses"]] + [str(retry)]
        self.assertEqual("accepted", preflight(result["prompts"], responses, reading_manifest=result["reading_manifest"])["status"])
        output = self.root / "final.json"
        merged = merge(result["prompts"], responses, str(output), reading_manifest=result["reading_manifest"])
        self.assertEqual(2, merged["read_receipts"]["reused_prior_assessment"])
        self.assertEqual(1, merged["read_receipts"]["complete"])
        rows = json.loads(output.read_text())["candidates"]
        self.assertEqual({"original-agent", "new-agent"}, {r["attempt_executor"] for r in rows})
        imported = import_route_b_semantic_proposals(result["bundle"], result["prompts"], str(output), str(self.root / "next-import"))
        self.assertEqual(3, json.loads(Path(imported["manifest"]).read_text())["response_count"])
        resumed = merge(result["prompts"], [str(output)], str(self.root / "resumed.json"), reading_manifest=result["reading_manifest"])
        self.assertEqual(3, resumed["completed_count"])
        self.assertEqual(2, resumed["read_receipts"]["reused_prior_assessment"])

    def test_missing_audit_is_created_and_duplicates_are_losslessly_coalesced(self):
        duplicate = self.root / "same-discovery.json"
        duplicate.write_bytes(self.discovery.read_bytes())
        result = self.extend(audited=False, discoveries=[str(self.discovery), str(duplicate)])
        self.assertEqual((2, 1), (result["reusable_count"], result["new_execution_required_count"]))
        self.assertTrue((self.root / "next" / "prior_import_audit").is_dir())
        combined = json.loads(Path(result["discoveries"]).read_text())
        self.assertEqual(1, len(combined["discovered_sources"]))
        self.assertEqual(2, len(combined["discovery_provenance"]))

    def test_foreign_parent_and_raw_tampering_rejected(self):
        data = json.loads(self.discovery.read_text())
        data["parent_bundle_sha256"] = "0" * 64
        self.write(self.discovery, data)
        with self.assertRaisesRegex(ValueError, "case/parent hash"):
            self.extend()
        self.assertFalse((self.root / "next").exists())
        data["parent_bundle_sha256"] = self.parent_sha
        self.write(self.discovery, data)
        result = self.extend()
        reuse = Path(result["reusable_responses"][0]["file"])
        modified = json.loads(reuse.read_text())
        modified["model"] = "different-model"
        self.write(reuse, modified)
        with self.assertRaisesRegex(ValueError, "eligible original"):
            merge(result["prompts"], [str(reuse)], str(self.root / "bad.json"))

    def test_old_raw_changed_after_extension_is_rejected(self):
        result = self.extend()
        self.raw.write_text(self.raw.read_text() + " ", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "raw response hash mismatch"):
            merge(result["prompts"], [r["file"] for r in result["reusable_responses"]], str(self.root / "bad.json"))

    def test_cli_creates_generation_without_model_calls(self):
        repo = Path(__file__).resolve().parents[1]
        process = subprocess.run([sys.executable, str(repo / "scripts/agentic_route_b.py"), "extend",
                                  "--bundle", self.bundle["bundle"], "--bundle-sha256", self.parent_sha,
                                  "--old-prompts", self.bundle["prompts"], "--old-response", str(self.checkpoint),
                                  "--discovery", str(self.discovery), "--import-manifest", self.import_manifest,
                                  "--out", str(self.root / "cli")], cwd=repo, capture_output=True, text=True)
        self.assertEqual(0, process.returncode, process.stderr)
        self.assertEqual(2, json.loads(process.stdout)["reusable_count"])


if __name__ == "__main__":
    unittest.main()
