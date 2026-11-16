"""Regression for candidate-local audit paths in real boundary handoffs."""
from __future__ import annotations

import copy
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.agentic_boundary import merge, prepare
from scripts.agentic_boundary_material import (
    _expand_tree, restore_reading_shard, validate_material_artifacts,
    validate_reading_projections,
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False), encoding="utf-8")


class BoundaryCandidateProjectionTests(unittest.TestCase):
    def prepare_case(self, root):
        ids = ["copy.alpha", "copy.beta", "copy.gamma"]
        higher_ids = ["higher.delta"]
        boundary, higher = root / "boundary.tsv", root / "higher.tsv"
        pd.DataFrame([{"candidate_gene_id": cid, "relevance_tier": "T2", "reason": "ambiguous copies"}
                      for cid in ids]).to_csv(boundary, sep="\t", index=False)
        pd.DataFrame([{"candidate_gene_id": higher_ids[0], "relevance_tier": "T1"}]).to_csv(higher, sep="\t", index=False)
        request = root / "request.json"
        write(request, {"stage": "boundary_review_request", "state": "awaiting_review", "case_id": "projection-case",
                        "source_wet_sha256": sha(boundary), "candidate_count": len(ids), "candidate_gene_ids": ids,
                        "candidate_input": {"path": str(boundary), "sha256": sha(boundary)},
                        "policy": ["preserve ambiguous copies and counterevidence"],
                        "higher_tier_context": {"state": "read_only", "candidate_gene_ids": higher_ids,
                                                "candidate_input": {"path": str(higher), "sha256": sha(higher)}}})
        materials = {}
        for index, cid in enumerate(ids + higher_ids):
            # Both direct and list-nested metadata exercise path-sensitive summaries.
            paper = {"tier_aggregation": "strongest_paper_v1",
                     "relevance_tier": ["T1_reported_focal_trait", "T2_reported_similar_phenomenon",
                                        "T3_reported_shared_process"][index % 3],
                     "tier_evidence_json": json.dumps({"source_id": f"PMID:{100 + index}", "reported_gene": cid}),
                     "paper_evidence_summary_json": json.dumps([
                         {"source_id": f"PMID:{100 + index}", "reason": "original perturbation"},
                         {"source_id": f"DOI:10.1234/{index}", "reason": "independent contrary observation"}]),
                     "relevance_evidence_json": json.dumps([{"source_id": f"PMID:{100 + index}",
                         "counterevidence": f"copy-specific conflict {cid}", "text": "full original observation " * 30}]),
                     "paper_reports_json": json.dumps([{"gene": cid, "source_excerpt": "unabridged source " * 40}]),
                     "source_ids": [f"PMID:{100 + index}", f"DOI:10.1234/{index}"]}
            materials[cid] = {**paper, "nested": [{"paper": copy.deepcopy(paper)}],
                              "copy_binding": {"candidate_gene_id": cid, "copy_index": index},
                              "source_text": "shared family context " * 30}
        context = root / "context.json"
        write(context, {"case_id": "projection-case", "phenotype": "recorded developmental phenotype",
                        "candidate_material": materials})
        prepared = root / "prepared"
        manifest = prepare(str(request), str(prepared), max_candidates=1, comparison_context_path=str(context))
        responses = []
        for index, shard in enumerate(manifest["shards"]):
            path = root / f"response-{index}.json"
            write(path, {"case_id": manifest["case_id"], "source_wet_sha256": manifest["source_wet_sha256"],
                         "comparison_context_sha256": manifest["comparison_context"]["sha256"],
                         "shard_id": shard["shard_id"], "executor_id": "test", "executor_model": "fixture",
                         "responses": [{"candidate_gene_id": cid, "decision": "adopt" if index != 1 else "defer",
                                        "reason": "copy comparison with retained conflicting reports"}
                                       for cid in shard["candidate_gene_ids"]]})
            responses.append(str(path))
        return prepared, manifest, materials, responses

    def test_real_prepare_validate_restore_and_proposal_merge(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            prepared, manifest, materials, responses = self.prepare_case(root)
            before = {entry["file"]: (prepared / entry["file"]).read_bytes() for entry in manifest["shards"]}
            validate_material_artifacts(prepared, manifest)
            validate_reading_projections(prepared, manifest)
            catalog = json.loads((prepared / "material_catalog.json").read_text(encoding="utf-8"))["materials"]
            for shard, descriptor in zip(manifest["shards"], manifest["reading_views"]):
                restored = restore_reading_shard(prepared, manifest, shard["shard_id"])
                expected_ids = shard["candidate_gene_ids"] + ["higher.delta"]
                self.assertEqual({cid: materials[cid] for cid in expected_ids}, restored["candidate_comparison_material"])
                view = json.loads((prepared / descriptor["file"]).read_text(encoding="utf-8"))
                self.assertEqual(set(expected_ids), set(view["candidate_material"]))
                for cid in expected_ids:
                    projection = _expand_tree(view["candidate_material"][cid], catalog)
                    self.assertEqual("paper_reports_json", projection["paper_reports_json"]["full_value_in"])
                    self.assertEqual("nested[0].paper.source_ids", projection["nested"][0]["paper"]["source_ids"]["full_value_in"])
                    self.assertEqual(json.loads(materials[cid]["paper_evidence_summary_json"]), projection["paper_evidence_summary_json"])
                    self.assertEqual(materials[cid]["copy_binding"], projection["copy_binding"])
                    self.assertEqual([f"copy-specific conflict {cid}"], projection["known_counterevidence"])
                    self.assertEqual(hashlib.sha256(json.dumps(materials[cid]["source_ids"], ensure_ascii=False,
                        sort_keys=True, separators=(",", ":")).encode()).hexdigest(), projection["source_ids"]["sha256"])
            output = root / "proposal.tsv"
            result = merge(str(prepared), responses, str(output), role="proposal", executor_id="test", executor_model="fixture")
            rows = pd.read_csv(output, sep="\t", dtype=str)
            self.assertEqual(list(materials)[:3], rows["candidate_gene_id"].tolist())
            self.assertEqual(3, result["candidate_count"])
            self.assertEqual(before, {name: (prepared / name).read_bytes() for name in before})

    def test_rehashed_view_tampering_still_fails_semantic_validation_and_merge(self):
        for tamper in ("audit_path", "source_summary", "copy_binding", "candidate_roster"):
            with self.subTest(tamper=tamper), tempfile.TemporaryDirectory() as td:
                root = Path(td)
                prepared, manifest, _, responses = self.prepare_case(root)
                catalog = json.loads((prepared / "material_catalog.json").read_text(encoding="utf-8"))["materials"]
                descriptor = manifest["reading_views"][0]
                path = prepared / descriptor["file"]
                view = json.loads(path.read_text(encoding="utf-8"))
                cid = manifest["shards"][0]["candidate_gene_ids"][0]
                material = _expand_tree(view["candidate_material"][cid], catalog)
                if tamper == "audit_path":
                    material["paper_reports_json"]["full_value_in"] = cid + ".paper_reports_json"
                elif tamper == "source_summary":
                    material["paper_evidence_summary_json"] = "[]"
                elif tamper == "copy_binding":
                    material["copy_binding"]["candidate_gene_id"] = "copy.beta"
                if tamper == "candidate_roster":
                    del view["candidate_material"][cid]
                else:
                    view["candidate_material"][cid] = {"$inline": material}
                write(path, view)
                descriptor["sha256"] = sha(path)
                write(prepared / "manifest.json", manifest)
                # Integrity checks pass deliberately: semantic projection still must reject.
                validate_material_artifacts(prepared, manifest)
                with self.assertRaisesRegex(ValueError, "does not restore canonical shard fields"):
                    validate_reading_projections(prepared, manifest)
                with self.assertRaisesRegex(ValueError, "does not restore canonical shard fields"):
                    merge(str(prepared), responses, str(root / "rejected.tsv"), role="proposal",
                          executor_id="test", executor_model="fixture")
                self.assertFalse((root / "rejected.tsv").exists())


if __name__ == "__main__":
    unittest.main()
