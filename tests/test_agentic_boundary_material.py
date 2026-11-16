from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.agentic_boundary import merge, prepare
from scripts.agentic_boundary_material import restore_reading_shard, expand_reading_view


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class BoundaryMaterialTests(unittest.TestCase):
    def test_higher_material_counts_actual_shard_placements_without_changing_raw(self):
        for max_candidates, higher_count in ((1, 1), (2, 3), (5, 0), (5, 1)):
            with self.subTest(max_candidates=max_candidates, higher_count=higher_count), tempfile.TemporaryDirectory() as td:
                root = Path(td)
                _, request_path, ids = self.fixture(root, 5)
                higher_ids = [f"H{index}" for index in range(higher_count)]
                higher = root / "higher.tsv"
                pd.DataFrame([{"candidate_gene_id": cid, "relevance_tier": "T1"} for cid in higher_ids],
                             columns=["candidate_gene_id", "relevance_tier"]).to_csv(higher, sep="\t", index=False)
                request = json.loads(request_path.read_text())
                request["higher_tier_context"]["candidate_gene_ids"] = higher_ids
                request["higher_tier_context"]["candidate_input"]["sha256"] = sha(higher)
                request_path.write_text(json.dumps(request))
                materials = {cid: {"gene_id": cid, "variant": f"variant-{cid}",
                                   "host_species": "Recorded host", "source_id": f"source-{cid}",
                                   "counterevidence": "Original unknown", "text": f"source text {cid} " * 150}
                             for cid in ids + higher_ids}
                context = root / "context.json"
                context.write_text(json.dumps({"case_id": "case-m", "phenotype": "Recorded neutral phenotype",
                                               "candidate_material": materials}))
                originals = {path: path.read_bytes() for path in (request_path, context, higher)}
                prepared = root / "prepared"
                manifest = prepare(str(request_path), str(prepared), max_candidates=max_candidates,
                                   comparison_context_path=str(context))
                catalog = json.loads((prepared / "material_catalog.json").read_text())["materials"]
                for descriptor, shard in zip(manifest["reading_views"], manifest["shards"]):
                    raw = (prepared / shard["file"]).read_bytes()
                    self.assertEqual(shard["sha256"], hashlib.sha256(raw).hexdigest())
                    canonical = json.loads(raw)
                    view = json.loads((prepared / descriptor["file"]).read_text())
                    expanded = expand_reading_view(view, catalog)
                    for cid in shard["candidate_gene_ids"] + higher_ids:
                        self.assertEqual(materials[cid], canonical["candidate_comparison_material"][cid])
                        self.assertEqual(materials[cid], expanded["candidate_material"][cid])
                    if len(manifest["shards"]) > 1:
                        for cid in higher_ids:
                            self.assertEqual({"$ref"}, set(view["candidate_material"][cid]))
                    self.assertEqual(canonical, restore_reading_shard(prepared, manifest, shard["shard_id"]))
                    self.assertEqual(raw, (prepared / shard["file"]).read_bytes())
                self.assertEqual(originals, {path: path.read_bytes() for path in originals})

    def test_structured_json_shared_delivery_preserves_bindings_and_raw_strings(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            boundary, request_path, ids = self.fixture(root, 3)
            fact = {"citation": "DOI:10.5555/neutral;x", "reported_gene": "OriginalGene",
                    "proxy_gene": "ProxyGene", "host_species": "Host species", "variant": "variant-A",
                    "counterevidence": "Known conflicting assay", "unknown": "unknown",
                    "report": "UNIQUE_SCIENTIFIC_REPORT_LONG " + "complete source text " * 100}
            shared = json.dumps([fact], ensure_ascii=False, indent=3)
            rows = [{"candidate_gene_id": cid, "shared_notes_json": shared,
                     "copy_identity_json": json.dumps({"focal_copy": cid, "source_role": "proxy"}),
                     "invalid_json": "not JSON: {", "scalar_json": '"keep original literal"',
                     "ambiguous_json": '{"same": 1, "same": 2}', "nonstandard_json": '[NaN]'} for cid in ids]
            pd.DataFrame(rows).to_csv(boundary, sep="\t", index=False)
            request = json.loads(request_path.read_text())
            request["source_wet_sha256"] = request["candidate_input"]["sha256"] = sha(boundary)
            request_path.write_text(json.dumps(request))
            prepared = root / "prepared"
            manifest = prepare(str(request_path), str(prepared), max_candidates=1)
            catalog = json.loads((prepared / "material_catalog.json").read_text())["materials"]
            self.assertEqual(1, json.dumps(catalog).count("UNIQUE_SCIENTIFIC_REPORT_LONG"))
            for descriptor in manifest["reading_views"]:
                view = expand_reading_view(json.loads((prepared / descriptor["file"]).read_text()), catalog)
                for cid, row in view["candidate_input"].items():
                    self.assertEqual([fact], row["shared_notes_json"])
                    self.assertEqual(cid, row["copy_identity_json"]["focal_copy"])
                    self.assertEqual("not JSON: {", row["invalid_json"])
                    self.assertEqual('"keep original literal"', row["scalar_json"])
                    self.assertEqual('{"same": 1, "same": 2}', row["ambiguous_json"])
                    self.assertEqual('[NaN]', row["nonstandard_json"])
                restored = restore_reading_shard(prepared, manifest, descriptor["shard_id"])
                self.assertEqual(shared, restored["candidate_input"][0]["shared_notes_json"])

    def test_paper_reading_view_keeps_all_summaries_and_restores_full_audit(self):
        for tier in ("T1_reported_focal_trait", "T2_reported_similar_phenomenon", "T3_reported_shared_process"):
            with self.subTest(tier=tier):
                self._assert_paper_reading_view(tier)

    def _assert_paper_reading_view(self, tier):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            boundary, request_path, ids = self.fixture(root, 1)
            representative = json.dumps({"source_id": "PMID:123", "initial_tier": tier[:2], "reported_gene": "REF1"})
            weak = json.dumps([{"source_id": f"PMID:{i}", "initial_tier": "T4",
                                "counterevidence": "Previously found conflict in ledger",
                                "text": "UNNECESSARY_WEAK_PAPER " * 100} for i in range(1000, 1100)])
            original = {"candidate_gene_id": ids[0], "tier_aggregation": "strongest_paper_v1",
                        "relevance_tier": tier, "tier_source_id": "PMID:123",
                        "tier_evidence_json": representative, "relevance_evidence_json": weak,
                        "counterevidence": "Known conflict must stay visible", "copy_state": "ambiguous"}
            original["paper_evidence_summary_json"] = json.dumps([
                {"source_id": "PMID:123", "tier": tier[:2], "publication_type": ["review"]},
                *[{"source_id": f"PMID:{i}", "tier": "T4", "publication_type": ["original"]} for i in range(1000, 1100)]])
            original["paper_tier_counts_json"] = json.dumps({tier[:2]: 1, "T4": 100})
            pd.DataFrame([original]).to_csv(boundary, sep="\t", index=False)
            request = json.loads(request_path.read_text())
            request["source_wet_sha256"] = request["candidate_input"]["sha256"] = sha(boundary)
            request_path.write_text(json.dumps(request))
            prepared = root / "prepared"
            manifest = prepare(str(request_path), str(prepared))
            view_path = prepared / manifest["reading_views"][0]["file"]
            view = json.loads(view_path.read_text())
            view = expand_reading_view(view, json.loads((prepared / "material_catalog.json").read_text())["materials"])
            projected = view["candidate_input"][ids[0]]
            self.assertEqual(json.loads(representative), projected["tier_evidence_json"])
            self.assertEqual(json.loads(original["paper_evidence_summary_json"]), projected["paper_evidence_summary_json"])
            self.assertEqual(json.loads(original["paper_tier_counts_json"]), projected["paper_tier_counts_json"])
            self.assertTrue(projected["relevance_evidence_json"]["audit_only"])
            self.assertEqual(original["counterevidence"], projected["counterevidence"])
            self.assertEqual(["Previously found conflict in ledger"], projected["known_counterevidence"])
            self.assertEqual("ambiguous", projected["copy_state"])
            self.assertNotIn("UNNECESSARY_WEAK_PAPER", view_path.read_text())
            self.assertLess(view_path.stat().st_size, len(weak.encode()) // 10)
            restored = restore_reading_shard(prepared, manifest, manifest["shards"][0]["shard_id"])
            self.assertEqual([original], restored["candidate_input"])

    def fixture(self, root: Path, count: int):
        boundary = root / "boundary.tsv"
        higher = root / "higher.tsv"
        ids = [f"G{i:04d}" for i in range(count)]
        pd.DataFrame([{"candidate_gene_id": cid, "relevance_tier": "T2", "reason": "provided", "family": "Kinase"}
                      for cid in ids]).to_csv(boundary, sep="\t", index=False)
        pd.DataFrame([{"candidate_gene_id": "H0", "relevance_tier": "T1", "family": "Kinase"}]).to_csv(higher, sep="\t", index=False)
        request = root / "request.json"
        request.write_text(json.dumps({"stage": "boundary_review_request", "state": "awaiting_review", "case_id": "case-m",
                                       "source_wet_sha256": sha(boundary), "candidate_count": count,
                                       "candidate_gene_ids": ids, "candidate_input": {"path": str(boundary), "sha256": sha(boundary)},
                                       "higher_tier_context": {"state": "read_only", "candidate_gene_ids": ["H0"],
                                                                "candidate_input": {"path": str(higher), "sha256": sha(higher)}}}), encoding="utf-8")
        return boundary, request, ids

    @staticmethod
    def response_files(root: Path, manifest: dict, ids: list[str]):
        paths = []
        for shard in manifest["shards"]:
            path = root / f"{shard['shard_id']}-response.json"
            path.write_text(json.dumps({"case_id": manifest["case_id"], "source_wet_sha256": manifest["source_wet_sha256"],
                                        "shard_id": shard["shard_id"], "executor_id": "test", "executor_model": "fixture",
                                        "comparison_context_sha256": manifest.get("comparison_context", {}).get("sha256"),
                                        "responses": [{"candidate_gene_id": cid, "decision": "adopt", "reason": "retained"}
                                                      for cid in shard["candidate_gene_ids"]]}), encoding="utf-8")
            paths.append(str(path))
        return paths

    def test_subtree_interning_preserves_family_source_and_copy_differences(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            boundary, request, ids = self.fixture(root, 3)
            context_path = root / "context.json"
            common_unit = {"description": "shared kinase family background " * 20, "citations": ["paper-A"]}
            materials = {
                ids[0]: {"family": "Kinase", "source": {"id": "SRC-A", "text": "distinct evidence A"},
                         "unit": common_unit, "paper_uri": "https://example.org/paper-A", "source_chunks": ["A" * 300],
                         "copy_binding": {"gene_id": ids[0], "copy": "1"}},
                ids[1]: {"family": "Kinase", "source": {"id": "SRC-B", "text": "distinct evidence B"},
                         "unit": common_unit, "paper_uri": "https://example.org/paper-B", "source_chunks": ["B" * 300],
                         "copy_binding": {"gene_id": ids[1], "copy": "2"}},
                ids[2]: {"family": "Kinase", "source": {"id": "SRC-A", "text": "distinct evidence A"},
                         "unit": common_unit, "paper_uri": "https://example.org/paper-A", "source_chunks": ["A" * 300],
                         "copy_binding": {"gene_id": ids[2], "copy": "3"}},
                "H0": {"family": "Kinase", "source": {"id": "SRC-H", "text": "read only"},
                      "unit": common_unit, "copy_binding": {"gene_id": "H0", "copy": "4"}},
            }
            context_path.write_text(json.dumps({"case_id": "case-m", "phenotype": "recorded phenotype",
                                                "candidate_material": materials}), encoding="utf-8")
            prepared = root / "prepared"
            manifest = prepare(str(request), str(prepared), max_candidates=3, comparison_context_path=str(context_path))
            catalog = json.loads((prepared / "material_catalog.json").read_text())["materials"]
            roots = [json.loads((prepared / entry["file"]).read_text())["candidate_material"] for entry in manifest["reading_views"]]
            self.assertNotEqual(roots[0][ids[0]], roots[0][ids[1]])
            self.assertNotEqual(roots[0][ids[0]], roots[0][ids[2]])
            # Exact common unit node is interned once, while distinct source and copy payloads stay bound.
            material0 = catalog[roots[0][ids[0]]["$ref"]]
            material1 = catalog[roots[0][ids[1]]["$ref"]]
            self.assertEqual(material0["value"]["unit"], material1["value"]["unit"])
            self.assertIn(material0["value"]["unit"]["$ref"], catalog)
            index = json.loads((prepared / manifest["comparison_index_artifact"]["file"]).read_text())["comparison_index"]
            self.assertEqual(ids + ["H0"], [row["candidate_gene_id"] for row in index])
            self.assertEqual({"file", "sha256"}, set(index[0]["canonical_shard"]))
            self.assertTrue(any(value.get("truncated") is True for value in index[0]["fields"].values()
                                if isinstance(value, dict)))

    def test_invalid_higher_roster_rejected_before_creating_output(self):
        for higher_ids, message in (([""], "non-empty unique"), (["H0", "H0"], "non-empty unique"),
                                    (["G0001"], "disjoint")):
            with self.subTest(higher_ids=higher_ids), tempfile.TemporaryDirectory() as td:
                root = Path(td)
                _, request, _ = self.fixture(root, 2)
                higher = root / "higher.tsv"
                pd.DataFrame([{"candidate_gene_id": cid, "relevance_tier": "T1"}
                              for cid in higher_ids]).to_csv(higher, sep="\t", index=False)
                payload = json.loads(request.read_text())
                payload["higher_tier_context"]["candidate_gene_ids"] = higher_ids
                payload["higher_tier_context"]["candidate_input"]["sha256"] = sha(higher)
                request.write_text(json.dumps(payload), encoding="utf-8")
                prepared = root / "prepared"
                with self.assertRaisesRegex(ValueError, message):
                    prepare(str(request), str(prepared), max_candidates=1)
                self.assertFalse(prepared.exists())

    def test_catalog_and_view_tampering_rejected_and_merge_keeps_full_roster(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            boundary, request, ids = self.fixture(root, 4)
            context = root / "context.json"
            context.write_text(json.dumps({"case_id": "case-m", "phenotype": "recorded phenotype",
                                          "candidate_material": {cid: {"shared": {"source": "same"}, "copy": cid}
                                                                 for cid in ids + ["H0"]}}), encoding="utf-8")
            prepared = root / "prepared"
            manifest = prepare(str(request), str(prepared), max_candidates=2, comparison_context_path=str(context))
            restored = restore_reading_shard(prepared, manifest, manifest["shards"][0]["shard_id"])
            canonical = json.loads((prepared / manifest["shards"][0]["file"]).read_text())
            self.assertEqual(canonical, restored)
            responses = self.response_files(root, manifest, ids)
            output = root / "proposal.tsv"
            merged = merge(str(prepared), responses, str(output), role="proposal", executor_id="test", executor_model="fixture")
            self.assertEqual(ids, pd.read_csv(output, sep="\t", dtype=str)["candidate_gene_id"].tolist())
            self.assertEqual(4, merged["candidate_count"])
            first_view = prepared / manifest["reading_views"][0]["file"]
            first_view.write_text(first_view.read_text() + " ", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "reading view hash mismatch"):
                merge(str(prepared), responses, str(root / "rejected.tsv"), role="proposal", executor_id="test", executor_model="fixture")

    def test_legacy_manifest_without_derived_artifacts_remains_mergeable(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            _, request, ids = self.fixture(root, 2)
            prepared = root / "prepared"
            manifest = prepare(str(request), str(prepared), max_candidates=2)
            responses = self.response_files(root, manifest, ids)
            manifest_path = prepared / "manifest.json"
            old_shape = json.loads(manifest_path.read_text())
            pages = (prepared / "reading_pages").resolve()
            self.assertTrue(pages.is_relative_to(root.resolve()))
            shutil.rmtree(pages)
            for task in old_shape["reading_pages"]["tasks"]:
                (prepared / task["file"]).unlink()
            for key in ("material_catalog", "comparison_index_artifact", "reading_views", "material_delivery_bytes",
                        "reading_pages", "reading_delivery_contract"):
                old_shape.pop(key, None)
            manifest_path.write_text(json.dumps(old_shape), encoding="utf-8")
            sidecar = merge(str(prepared), responses, str(root / "legacy.tsv"), role="proposal",
                            executor_id="test", executor_model="fixture")
            self.assertEqual(2, sidecar["candidate_count"])

    def test_direct_script_cli_prepare_and_merge(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            _, request, ids = self.fixture(root, 2)
            prepared = root / "prepared"
            script = Path(__file__).resolve().parents[1] / "scripts" / "agentic_boundary.py"
            subprocess.run([sys.executable, str(script), "prepare", "--request", str(request), "--out", str(prepared)],
                           cwd=script.parents[1], check=True, capture_output=True, text=True)
            manifest = json.loads((prepared / "manifest.json").read_text())
            responses = self.response_files(root, manifest, ids)
            for path in responses:
                payload = json.loads(Path(path).read_text())
                payload["executor_id"] = "cli"
                Path(path).write_text(json.dumps(payload), encoding="utf-8")
            output = root / "cli-proposal.tsv"
            completed = subprocess.run([sys.executable, str(script), "merge", "--prepared", str(prepared),
                                         "--role", "proposal", "--out", str(output), "--executor-id", "cli",
                                         "--executor-model", "fixture", *sum((["--response", path] for path in responses), [])],
                                        cwd=script.parents[1], check=True, capture_output=True, text=True)
            self.assertEqual(2, json.loads(completed.stdout)["candidate_count"])

    def test_material_budget_splits_without_truncation_and_marks_singleton_oversize(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            _, request, ids = self.fixture(root, 5)
            context = root / "context.json"
            payload = {cid: {"sequence_evidence": cid + "-" + ("X" * 80)} for cid in ids + ["H0"]}
            context.write_text(json.dumps({"case_id": "case-m", "phenotype": "recorded phenotype", "candidate_material": payload}), encoding="utf-8")
            prepared = root / "prepared"
            manifest = prepare(str(request), str(prepared), max_candidates=5, max_material_bytes=100, comparison_context_path=str(context))
            self.assertGreater(len(manifest["shards"]), 1)
            self.assertEqual(ids, [cid for shard in manifest["shards"] for cid in shard["candidate_gene_ids"]])
            self.assertEqual(ids, manifest["material_delivery_bytes"]["oversized_gene_ids"])
            view = json.loads((prepared / manifest["reading_views"][0]["file"]).read_text())
            self.assertTrue(view["chunk_reading_required"])
            view_text = json.dumps(view)
            catalog_text = json.dumps(json.loads((prepared / "material_catalog.json").read_text()))
            self.assertTrue(payload[view["oversized_material_gene_ids"][0]]["sequence_evidence"] in view_text
                            or payload[view["oversized_material_gene_ids"][0]]["sequence_evidence"] in catalog_text)

    def test_shared_and_singleton_payload_scaling_100_300_1000(self):
        for count in (100, 300, 1000):
            for singleton in (False, True):
                with self.subTest(count=count, singleton=singleton), tempfile.TemporaryDirectory() as td:
                    root = Path(td)
                    _, request, ids = self.fixture(root, count)
                    context = root / "context.json"
                    materials = {cid: {"family_source_unit": ({"source": "SRC-" + cid, "unit": "shared"} if singleton
                                                                          else {"source": "SRC-SHARED", "unit": "shared" + (" background " * 100)}),
                                       "copy_binding": cid} for cid in ids + ["H0"]}
                    context.write_text(json.dumps({"case_id": "case-m", "phenotype": "recorded phenotype",
                                                   "candidate_material": materials}), encoding="utf-8")
                    prepared = root / "prepared"
                    manifest = prepare(str(request), str(prepared), max_candidates=100,
                                       comparison_context_path=str(context))
                    self.assertEqual(ids, [cid for shard in manifest["shards"] for cid in shard["candidate_gene_ids"]])
                    self.assertGreater(manifest["material_delivery_bytes"]["expanded_read_bytes"], 0)
                    self.assertLess(len(manifest["shards"]), count)
                    view = json.loads((prepared / manifest["comparison_index_artifact"]["file"]).read_text())
                    index = view["comparison_index"]
                    self.assertLess(max(len(json.dumps(row, separators=(",", ":")).encode())
                                       for row in index), 2_000)
                    if not singleton:
                        self.assertLess(manifest["material_delivery_bytes"]["unique_shared_material_bytes"],
                                        len(json.dumps(materials, separators=(",", ":"))))


if __name__ == "__main__":
    unittest.main()
