"""Frozen secondary BLAST handles reach research and claims without support."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from candidate_qtl.route_b_semantic_bundle import build_route_b_semantic_bundle, refreeze_route_b_sources
from candidate_qtl.route_b_semantic_candidates import freeze_route_b_semantic_candidates
from candidate_qtl.route_b_semantic_import import import_route_b_semantic_proposals
from candidate_qtl.route_b_semantic_source_verify import verify_route_b_semantic_sources
from candidate_qtl.route_b_semantic_transfer_verify import verify_route_b_semantic_transfer


class RouteBIdentityEndpointWiringTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        atlas = self.root / "atlas.tsv"
        pd.DataFrame([
            {"gene_id": "G1", "family_id": "F1", "family_index_subject_endpoints": "P11111",
             "local_annotation_description": "meristem signaling protein"},
            {"gene_id": "G2", "family_id": "F2", "family_index_subject_endpoints": "P33333",
             "local_annotation_description": "related developmental signaling protein"},
        ]).to_csv(atlas, sep="\t", index=False)
        hits = self.root / "hits.tsv"
        pd.DataFrame([
            {"focal_gene_id": gene, "hit_accession": ref, "hit_species": "Arabidopsis thaliana",
             "hit_protein_name": "meristem signaling protein", "hit_metadata_state": "complete",
             "percent_identity": "80", "query_coverage": "0.9"}
            for gene, ref in (("G1", "P11111"), ("G1", "P22222"), ("G2", "P33333"), ("G2", "P44444"))
        ]).to_csv(hits, sep="\t", index=False)
        self.result = build_route_b_semantic_bundle(
            str(atlas), str(self.root / "bundle"), phenotype="altered meristem development",
            identity_hits_arg=str(hits), case_id="secondary-hit-fixture", focal_species="Solanum lycopersicum",
        )
        self.bundle_path = Path(self.result["bundle"])
        self.bundle_bytes = self.bundle_path.read_bytes()
        self.bundle = json.loads(self.bundle_bytes)
        self.unit = next(unit for unit in self.bundle["units"] if unit["member_interval_genes"] == ["G1"])

    def _import(self, result, reference, name, *, source_ids=()):
        prompts = [json.loads(line) for line in Path(result["prompts"]).read_text().splitlines()]
        claim = {
            "claim_type": "mechanistic_review_seed",
            "biological_function": "regulation of meristem maintenance",
            "mechanism": "altered stem cell signaling could explain altered meristem development",
            "phenotype_relation": "conditional cross-species meristem hypothesis",
            "reference_gene": reference, "reference_species": "Arabidopsis thaliana",
            "evidence_scope": "reference_gene", "source_ids": list(source_ids),
            "alternative_hypothesis": "a different signaling protein affects the phenotype",
            "uncertainty_note": "sequence identity is a retrieval handle, not verified causal support",
        }
        response = self.root / f"{name}.json"
        response.write_text(json.dumps({"candidates": [
            {"prompt_id": prompt["prompt_id"], "parse_status": "json",
             "model_json": {"claims": [claim] if prompt["route_b_unit_id"] == self.unit["unit_id"] else []}}
            for prompt in prompts
        ]}))
        imported = import_route_b_semantic_proposals(
            result["bundle"], result["prompts"], str(response), str(self.root / name),
        )
        rows = pd.read_csv(imported["proposals"], sep="\t", dtype=str).fillna("")
        row = rows[rows.reference_gene.eq(reference)].iloc[0]
        return imported, row

    def _refreeze(self, reference, name):
        parent_sha = hashlib.sha256(self.bundle_bytes).hexdigest()
        discovery = self.root / f"{name}-discovery.json"
        discovery.write_text(json.dumps({
            "case_id": self.bundle["case_id"], "parent_bundle_sha256": parent_sha,
            "discovered_sources": [{
                "unit_id": self.unit["unit_id"], "reference_gene": reference,
                "reference_species": "Arabidopsis thaliana", "reported_gene": "REF2",
                "reported_species": "Arabidopsis thaliana", "source_id": "PMID:1234567",
                "source_title": "Synthetic meristem observation", "source_url": "https://example.test/fixture",
                "retrieved_at": "2026-10-19T00:00:00Z", "source_text": "Synthetic fixture material for REF2.",
                "reference_function": "regulation of meristem maintenance",
            }],
        }))
        return refreeze_route_b_sources(
            str(self.bundle_path), str(discovery), str(self.root / name), parent_sha256=parent_sha,
        )

    def test_secondary_hit_reaches_import_refreeze_and_initial_candidate_without_support(self):
        prompt = next(json.loads(line) for line in Path(self.result["prompts"]).read_text().splitlines()
                      if json.loads(line)["route_b_unit_id"] == self.unit["unit_id"])
        self.assertIn("P22222", prompt["user"])
        self.assertNotIn("P22222", self.unit["genes"][0]["family_index_subject_endpoints"])
        imported, row = self._import(self.result, "P22222", "before-research")
        self.assertEqual("proposed", row.verification_status)
        self.assertIn("reference_gene_grounded:P22222", row.verification_reason)
        self.assertEqual("", row.source_id)

        refrozen = self._refreeze("P22222", "refrozen")
        new_bundle = json.loads(Path(refrozen["bundle"]).read_text())
        self.assertEqual(self.bundle_bytes, self.bundle_path.read_bytes())
        for old_unit, new_unit in zip(self.bundle["units"], new_bundle["units"]):
            for field in ("unit_id", "member_interval_genes", "genes", "identity"):
                self.assertEqual(old_unit[field], new_unit[field])
            if old_unit["unit_id"] != self.unit["unit_id"]:
                self.assertEqual(old_unit, new_unit)
        imported, row = self._import(refrozen, "P22222", "after-research", source_ids=["PMID:1234567"])
        self.assertEqual("proposed", row.verification_status)
        cache = self.root / "cache.json"
        cache.write_text(json.dumps({"id": "1234567", "title": "Synthetic meristem observation"}))
        sources = verify_route_b_semantic_sources(
            imported["proposals"], str(self.root / "source"), bundle_arg=refrozen["bundle"], cache_arg=str(cache),
        )
        transfer = verify_route_b_semantic_transfer(
            imported["proposals"], sources["source_verification"], refrozen["bundle"],
            str(self.root / "transfer"), case_id=self.bundle["case_id"],
        )
        frozen = freeze_route_b_semantic_candidates(
            transfer["claims"], str(self.root / "freeze"), case_id=self.bundle["case_id"],
        )
        review = pd.read_csv(frozen["review_candidates"], sep="\t").fillna("")
        self.assertEqual({"G1"}, set(review.candidate_gene_id))
        self.assertEqual({"adopted_first_stage"}, set(review.candidate_admission_state))
        self.assertEqual({"unverified_route_b_review"}, set(review.support_status))
        self.assertEqual(0, frozen["candidate_gene_count"])

    def test_unacquired_and_other_unit_hits_remain_ungrounded(self):
        for reference in ("P99999", "P44444"):
            with self.subTest(reference=reference):
                _, row = self._import(self.result, reference, f"import-{reference}")
                self.assertEqual("unsupported", row.verification_status)
                self.assertIn(f"reference_gene_not_in_frozen_bundle:{reference}", row.verification_reason)
                with self.assertRaisesRegex(ValueError, "unfrozen reference endpoint"):
                    self._refreeze(reference, f"refreeze-{reference}")
        self.assertEqual(self.bundle_bytes, self.bundle_path.read_bytes())


if __name__ == "__main__":
    unittest.main()
