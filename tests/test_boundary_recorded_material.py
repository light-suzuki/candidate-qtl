"""Synthetic transport regressions; no model calls or scientific re-assessment."""
import copy
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.boundary_comparison_material import build_comparison_context
from candidate_qtl.boundary_tier_review import (
    UPPER_EVIDENCE_SCOPED_POLICY, build_boundary_review_request, freeze_boundary_scope,
)
from candidate_qtl.relevance_tiers import NEW_TIERS
from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from scripts.agentic_boundary import prepare
from scripts.agentic_boundary_material import restore_reading_shard
from tests._issue664_helpers import manifest_for
from tests.test_wet_relevance_selection import handoff


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False), encoding="utf-8")


class RecordedBoundaryMaterialTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.bundle = self.root / "bundle.json"
        self.payload = {
            "stage": "post_b0_semantic_bundle", "route": "interval_first", "case_id": "case-x",
            "phenotype": "recorded organ response",
            "units": [{"unit_id": "unit-pair", "member_interval_genes": ["G1", "G2"],
                       "genes": [{"gene_id": cid, "reference_function": "shared recorded function",
                                  "tissue_specificity": cid, "mapping_state": "unresolved",
                                  "sequence_similarity_hits": [{"hit_accession": "reference-x", "identity": 44,
                                      "query_coverage": 0.4, "qstart": 1, "qend": 40,
                                      "hit_metadata_records_json": "large provider dump" * 100}]}
                                 for cid in ["G1", "G2"]],
                       "reference_material_cards": [
                           {"interval_gene_id": "G1", "function_description": "copy-one context",
                            "source_id": "source-one", "reference_gene_id": "reference-one"},
                           {"function_description": "shared context", "source_id": "source-shared"}]},
                      {"unit_id": "unit-higher", "member_interval_genes": ["H1"],
                       "genes": [{"gene_id": "H1", "reference_function": "higher context"}]},
                      {"unit_id": "unit-lower", "member_interval_genes": ["L1"],
                       "genes": [{"gene_id": "L1", "reference_function": "lower-only material"}]}],
        }
        write(self.bundle, self.payload)
        rows = [{"candidate_gene_id": cid, "relevance_tier": tier, "support_state": "review"}
                for cid, tier in [("H1", NEW_TIERS[0]), ("G1", NEW_TIERS[1]),
                                  ("G2", NEW_TIERS[1]), ("L1", NEW_TIERS[3])]]
        self.inventory = self.root / "inventory.tsv"
        pd.DataFrame(rows).to_csv(self.inventory, sep="\t", index=False)
        self.scope = freeze_boundary_scope(str(self.inventory), str(self.root / "scope"),
            inventory_sha256=sha(self.inventory), case_id="case-x", target=1, upper=2,
            policy_version=UPPER_EVIDENCE_SCOPED_POLICY,
            eligible_gene_ids=[r["candidate_gene_id"] for r in rows])

    def request(self):
        return build_boundary_review_request(self.scope["manifest"], str(self.root / "request"),
            scope_sha256=sha(Path(self.scope["manifest"])), case_id="case-x",
            comparison_bundle_arg=str(self.bundle))

    def test_scoped_material_binding_and_source_bytes_are_preserved(self):
        before = self.bundle.read_bytes()
        context = build_comparison_context(str(self.bundle), "case-x", ["G1", "G2", "H1", "missing"])
        self.assertEqual(sha(self.bundle), context["source_bundle"]["sha256"])
        self.assertNotIn("L1", context["candidate_material"])
        g1, g2 = [context["candidate_material"][cid] for cid in ["G1", "G2"]]
        self.assertEqual("G1", g1["function_and_copy_context"]["tissue_specificity"])
        self.assertEqual("G2", g2["function_and_copy_context"]["tissue_specificity"])
        self.assertEqual(["gene", "unit_context"], [r["binding_scope"] for r in g1["reference_material_evidence"]])
        self.assertEqual(["unit_context"], [r["binding_scope"] for r in g2["reference_material_evidence"]])
        self.assertEqual(0.4, g1["alignment_evidence"][0]["query_coverage"])
        self.assertNotIn("hit_metadata_records_json", g1["alignment_evidence"][0])
        self.assertIn("sequence_similarity_hits/0", g1["alignment_evidence"][0]["source_pointer"])
        self.assertEqual({"material_state": "not_available_in_route_b_bundle"}, context["candidate_material"]["missing"])
        self.assertEqual(before, self.bundle.read_bytes())

    def test_recorded_context_automatically_reaches_real_prepare_and_restore(self):
        request = self.request()
        source = Path(request["request"])
        before = source.read_bytes()
        payload = json.loads(before)
        descriptor = payload["comparison_context"]
        self.assertEqual(sha(Path(descriptor["path"])), descriptor["sha256"])
        prepared = self.root / "prepared"
        manifest = prepare(str(source), str(prepared), max_candidates=1)
        expected = json.loads(Path(descriptor["path"]).read_text())["candidate_material"]
        for shard in manifest["shards"]:
            restored = restore_reading_shard(prepared, manifest, shard["shard_id"])
            scope = shard["candidate_gene_ids"] + ["H1"]
            self.assertEqual({cid: expected[cid] for cid in scope}, restored["candidate_comparison_material"])
            self.assertNotIn("L1", restored["candidate_comparison_material"])
        self.assertIn("not a claim of noncausality", " ".join(payload["policy"]))
        view = json.loads((prepared / manifest["reading_views"][0]["file"]).read_text())
        self.assertIn("Keep the original Tier separate", view["reader_instruction"])
        self.assertEqual(before, source.read_bytes())

    def test_additional_context_cannot_erase_recorded_material(self):
        request = self.request()
        original = json.loads(Path(request["request"]).read_text())["comparison_context"]
        expected = json.loads(Path(original["path"]).read_text())["candidate_material"]
        extra = self.root / "extra.json"
        write(extra, {"case_id": "case-x", "phenotype": self.payload["phenotype"],
                      "candidate_material": {cid: {} for cid in ["G1", "G2", "H1"]}})
        prepared = self.root / "extra-prepared"
        manifest = prepare(request["request"], str(prepared), comparison_context_path=str(extra))
        restored = restore_reading_shard(prepared, manifest, manifest["shards"][0]["shard_id"])
        for cid, material in restored["candidate_comparison_material"].items():
            self.assertEqual(expected[cid], material["recorded_source_context"])
            self.assertEqual({}, material["additional_source_context"])

    def test_additional_context_cannot_replace_recorded_phenotype(self):
        request = self.request()
        extra = self.root / "different-phenotype.json"
        write(extra, {"case_id": "case-x", "phenotype": "different observation",
                      "candidate_material": {cid: {} for cid in ["G1", "G2", "H1"]}})
        with self.assertRaisesRegex(ValueError, "recorded phenotype mismatch"):
            prepare(request["request"], str(self.root / "different"), comparison_context_path=str(extra))
        self.assertFalse((self.root / "different").exists())

    def test_recorded_context_tampering_is_rejected_before_output(self):
        request = self.request()
        descriptor = json.loads(Path(request["request"]).read_text())["comparison_context"]
        path = Path(descriptor["path"])
        path.write_text(path.read_text() + " ")
        with self.assertRaisesRegex(ValueError, "recorded comparison context SHA"):
            prepare(request["request"], str(self.root / "bad"))
        self.assertFalse((self.root / "bad").exists())

    def test_wrong_case_conflicting_card_and_duplicate_members_are_rejected(self):
        bad_cases = []
        wrong_case = copy.deepcopy(self.payload)
        wrong_case["case_id"] = "other"
        bad_cases.append(wrong_case)
        wrong_card = copy.deepcopy(self.payload)
        wrong_card["units"][0]["reference_material_cards"][0]["interval_gene_id"] = "L1"
        bad_cases.append(wrong_card)
        duplicate = copy.deepcopy(self.payload)
        duplicate["units"].append(copy.deepcopy(duplicate["units"][0]))
        bad_cases.append(duplicate)
        for payload in bad_cases:
            write(self.bundle, payload)
            with self.assertRaises(ValueError):
                build_comparison_context(str(self.bundle), "case-x", ["G1", "G2"])

    def test_real_wet_handoff_adds_material_without_changing_tiers_or_members(self):
        base = self.root / "wet-inputs"
        baseline = handoff(base, (1, 26, 2, 0, 0), policy=UPPER_EVIDENCE_SCOPED_POLICY)
        ids = pd.read_csv(base / "atlas.tsv", sep="\t")["gene_id"].tolist()
        write(self.bundle, {"stage": "post_b0_semantic_bundle", "route": "interval_first", "case_id": "case1",
            "phenotype": "recorded response", "units": [{"unit_id": "all", "member_interval_genes": ids,
                "genes": [{"gene_id": cid, "reference_function": "recorded function"} for cid in ids]}]})
        enriched = wet_candidate_handoff(str(self.root / "enriched-wet"),
            route_b_candidates_arg=str(base / "supported.tsv"),
            route_b_candidates_manifest_arg=str(manifest_for(base / "supported.tsv", route="b")),
            interval_atlas_arg=str(base / "atlas.tsv"), case_id="case1",
            review_candidates_arg=str(base / "review.tsv"), review_candidates_manifest_arg=str(base / "review_manifest.json"),
            boundary_route_b_bundle_arg=str(self.bundle))
        self.assertEqual(Path(baseline["wet_candidates"]).read_bytes(), Path(enriched["wet_candidates"]).read_bytes())
        request = json.loads(Path(enriched["adaptive_synthesis_request"]).read_text())
        self.assertEqual(26, request["candidate_count"])
        context = json.loads(Path(request["comparison_context"]["path"]).read_text())
        self.assertEqual(27, len(context["candidate_material"]))
        self.assertEqual("original_assessment_v1", request["assessment_mode"])


if __name__ == "__main__":
    unittest.main()
