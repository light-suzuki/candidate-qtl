from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.core import sha256
from candidate_qtl.route_c_claim_binding_producer import (
    prepare_claim_binding_review, seal_claim_binding_review,
)
from candidate_qtl.route_c_direct_pipeline import DIRECT_PROPOSAL_COLUMNS, direct_rescue_freeze


class CapturedSourceScopeTests(unittest.TestCase):
    def _freeze(self, root: Path, *, full_text: bool = True):
        quote = "GENE1 mutant showed altered growth in Synthetic species."
        proposal = dict.fromkeys(DIRECT_PROPOSAL_COLUMNS, "")
        proposal.update(
            record_id="r1", query_phenotype="altered growth", reported_phenotype="altered growth",
            focal_species="Synthetic species", reference_species="Synthetic species",
            reported_gene_id="GENE1", species_distance="S1_same_species",
            evidence_directness="D0_causal_perturbation_or_rescue",
            evidence_type="synthetic mutant observation", source_id="PMID:12345678",
        )
        proposals, cache, extraction = root / "proposals.tsv", root / "cache.json", root / "extraction.json"
        pd.DataFrame([proposal]).to_csv(proposals, sep="\t", index=False)
        source = {"id": "12345678", "pmid": "12345678", "title": "GENE1 growth study",
                  "abstractText": "GENE1 was studied." if full_text else quote}
        if full_text:
            source["fullText"] = quote
        cache.write_text(json.dumps(source), encoding="utf-8")
        extraction.write_text(json.dumps([{
            "citation": "PMID:12345678", "hypothesis_id": "r1",
            "evidence_relation": "genetic_association", "quoted_support": quote,
        }]), encoding="utf-8")
        result = direct_rescue_freeze(str(proposals), str(root / "freeze"), str(cache),
                                      claim_extraction_arg=str(extraction))
        manifest_path = Path(result["manifest"])
        manifest = json.loads(manifest_path.read_text())
        source_manifest = json.loads(Path(manifest["source_verification_manifest"]).read_text())
        capture = Path(source_manifest["outputs"]["source_review_records"])
        return result, manifest_path, capture

    def _prepare(self, root: Path, frozen):
        result = prepare_claim_binding_review(frozen["frozen"], frozen["manifest"], str(root / "prepared"))
        return pd.read_csv(result["review"], sep="\t", dtype=str, keep_default_na=False)

    def _seal(self, root: Path, frozen, review):
        path = root / "reviewed.tsv"
        review.to_csv(path, sep="\t", index=False)
        return seal_claim_binding_review(frozen["frozen"], frozen["manifest"], str(path), str(root / "bound"))

    def _rebind_capture(self, manifest_path: Path, capture: Path, payload):
        capture.write_text(json.dumps(payload), encoding="utf-8")
        manifest = json.loads(manifest_path.read_text())
        manifest["source_review_records_sha256"] = sha256(capture)
        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    def test_real_full_text_freeze_prepare_seal_keeps_classification_and_curation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            frozen, _, capture = self._freeze(root)
            source = json.loads(capture.read_text())["records"]["PMID:12345678"]
            self.assertEqual("full_text_cached", source["available_text_scope"])
            self.assertEqual("title_abstract", source["classification_scope"])
            validation = pd.read_csv(capture.parent / "source_validation.tsv", sep="\t")
            self.assertEqual("title_abstract", validation.iloc[0]["retrieval_scope"])
            review = self._prepare(root, frozen)
            self.assertEqual("full_text", review.iloc[0]["source_scope"])
            review.loc[0, "source_organism_state"] = "verified_match"
            review.loc[0, "reported_identity_verification_state"] = "verified_exact"
            review.loc[0, "identity_mapping_basis"] = "exact_source_text"
            review.loc[0, "experimental_target_gene_state"] = "verified_same_identity"
            review.loc[0, "phenotype_match_state"] = "target_trait"
            review.loc[0, "perturbation_to_phenotype_binding_state"] = "verified_supports"
            review.loc[0, "reviewed_by"] = "synthetic-agent"
            review.loc[0, "binding_reason"] = "source identity reviewed; no human attestation"
            sealed = self._seal(root, frozen, review)
            binding = pd.read_csv(sealed["binding"], sep="\t", dtype=str, keep_default_na=False)
            self.assertEqual("full_text", binding.iloc[0]["source_scope"])
            self.assertEqual(review.iloc[0]["quoted_support"], binding.iloc[0]["quoted_support"])
            self.assertEqual("review_required", binding.iloc[0]["binding_status"])
            self.assertFalse(sealed["complete"])

    def test_abstract_only_cannot_claim_full_text(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            frozen, _, _ = self._freeze(root, full_text=False)
            review = self._prepare(root, frozen)
            self.assertEqual("abstract", review.iloc[0]["source_scope"])
            review.loc[0, "source_scope"] = "full_text"
            with self.assertRaisesRegex(ValueError, "source_scope exceeds"):
                self._seal(root, frozen, review)

    def test_legacy_without_captured_artifact_replays_original_scope(self):
        for retrieval_scope, expected in [("title_abstract", "abstract"), ("full_text_cached", "full_text")]:
            with self.subTest(scope=retrieval_scope), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                frozen, manifest_path, _ = self._freeze(root)
                manifest = json.loads(manifest_path.read_text())
                manifest.pop("source_review_records_sha256")
                source_path = Path(manifest["source_verification_manifest"])
                source = json.loads(source_path.read_text())
                source["outputs"].pop("source_review_records")
                source_path.write_text(json.dumps(source))
                validation_path = Path(source["outputs"]["source_validation"])
                validation = pd.read_csv(validation_path, sep="\t", dtype=str, keep_default_na=False)
                validation["retrieval_scope"] = retrieval_scope
                validation.to_csv(validation_path, sep="\t", index=False)
                manifest["source_validation_sha256"] = sha256(validation_path)
                manifest_path.write_text(json.dumps(manifest))
                review = self._prepare(root, frozen)
                self.assertEqual(expected, review.iloc[0]["source_scope"])
                self.assertFalse(self._seal(root, frozen, review)["complete"])

    def test_missing_or_malformed_captured_availability_fails_closed(self):
        for value in [None, "full_text", [], "unknown"]:
            with self.subTest(value=value), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                frozen, manifest_path, capture = self._freeze(root)
                review = self._prepare(root, frozen)
                payload = json.loads(capture.read_text())
                record = payload["records"]["PMID:12345678"]
                if value is None:
                    record.pop("available_text_scope")
                else:
                    record["available_text_scope"] = value
                self._rebind_capture(manifest_path, capture, payload)
                with self.assertRaisesRegex(ValueError, "availability is missing or malformed"):
                    self._prepare(root, frozen)
                with self.assertRaisesRegex(ValueError, "availability is missing or malformed"):
                    self._seal(root, frozen, review)

    def test_availability_cannot_exceed_captured_sections(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            frozen, manifest_path, capture = self._freeze(root)
            payload = json.loads(capture.read_text())
            payload["records"]["PMID:12345678"]["sections"].pop("fullText")
            self._rebind_capture(manifest_path, capture, payload)
            with self.assertRaisesRegex(ValueError, "availability disagrees"):
                self._prepare(root, frozen)

    def test_missing_source_record_cannot_borrow_classification_scope(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            frozen, manifest_path, capture = self._freeze(root)
            self._rebind_capture(manifest_path, capture, {"records": {}})
            review = self._prepare(root, frozen)
            self.assertEqual("title", review.iloc[0]["source_scope"])
            review.loc[0, "source_scope"] = "full_text"
            with self.assertRaisesRegex(ValueError, "source_scope exceeds"):
                self._seal(root, frozen, review)

    def test_captured_artifact_hash_and_original_quote_still_bind(self):
        for tamper in ["capture", "quote", "missing_hash"]:
            with self.subTest(tamper=tamper), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                frozen, manifest_path, capture = self._freeze(root)
                review = self._prepare(root, frozen)
                if tamper == "capture":
                    capture.write_text(capture.read_text() + " ")
                    message = "captured-source artifact hash"
                elif tamper == "quote":
                    review.loc[0, "quoted_support"] = "An invented quote."
                    message = "not an accepted verbatim extraction"
                else:
                    manifest = json.loads(manifest_path.read_text())
                    manifest.pop("source_review_records_sha256")
                    manifest_path.write_text(json.dumps(manifest))
                    message = "requires artifact and hash"
                with self.assertRaisesRegex(ValueError, message):
                    self._seal(root, frozen, review)


if __name__ == "__main__":
    unittest.main()
