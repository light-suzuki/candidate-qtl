from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.semantic_claims import (
    CLAIM_COLUMNS,
    claim_content_fingerprint,
    freeze_semantic_claims,
    narrowing_eligible_claims,
    normalize_semantic_claims,
)


def _claim(**overrides: str) -> dict[str, str]:
    row = {column: "" for column in CLAIM_COLUMNS}
    row.update({
        "claim_id": "claim-1",
        "route": "interval_first",
        "semantic_unit_id": "bsem-fixture",
        "family_id": "FAM1",
        "orthogroup_id": "",
        "subfamily_or_clade": "clade-a",
        "member_interval_genes": "G_B;G_A;G_A",
        "claim_type": "mutant_phenotype",
        "biological_function": "meristem maintenance",
        "mechanism": "receptor signalling",
        "phenotype_relation": "plausibly explains fasciation",
        "reference_gene": "REF1",
        "reference_species": "Arabidopsis thaliana",
        "evidence_scope": "reference_gene",
        "source_id": "PMID:123456",
        "source_url": "",
        "alternative_hypothesis": "related receptor family",
        "uncertainty_note": "functional transfer remains family-level",
        "provider": "fixture",
        "model": "fixture-model",
        "prompt_fingerprint": "p" * 64,
        "input_fingerprint": "i" * 64,
        "output_fingerprint": "o" * 64,
        "verification_status": "verified",
        "verification_reason": "source exists and supports the stated reference-gene phenotype",
        "rank_contribution": "none",
    })
    row.update(overrides)
    return row


class SemanticClaimContractTests(unittest.TestCase):
    def test_route_b_claim_requires_frozen_semantic_unit_id(self) -> None:
        with self.assertRaisesRegex(ValueError, "semantic_unit_id"):
            normalize_semantic_claims(
                pd.DataFrame([_claim(semantic_unit_id="")], columns=CLAIM_COLUMNS)
            )

    def test_verified_claim_requires_source_and_reason(self) -> None:
        with self.assertRaisesRegex(ValueError, "source_id"):
            normalize_semantic_claims(pd.DataFrame([_claim(source_id="")], columns=CLAIM_COLUMNS))
        with self.assertRaisesRegex(ValueError, "verification_reason"):
            normalize_semantic_claims(pd.DataFrame([_claim(verification_reason="")], columns=CLAIM_COLUMNS))

    def test_proposal_and_provider_failure_never_enter_narrowing(self) -> None:
        frame = pd.DataFrame([
            _claim(claim_id="verified"),
            _claim(claim_id="proposal", verification_status="proposed", source_id="", verification_reason=""),
            _claim(claim_id="provider", verification_status="provider_failed", source_id="",
                   verification_reason="provider unavailable"),
            _claim(claim_id="not-evaluated", verification_status="not_evaluated", source_id="",
                   verification_reason="not requested"),
        ], columns=CLAIM_COLUMNS)
        eligible = narrowing_eligible_claims(frame)
        self.assertEqual(["verified"], eligible["claim_id"].tolist())

    def test_mechanistic_review_seed_is_structurally_review_only(self) -> None:
        proposed = _claim(
            claim_type="mechanistic_review_seed",
            verification_status="proposed",
            verification_reason="no_source_id_supplied; semantic_relation_verification_pending",
            source_id="",
            mechanism="altered meristem signaling could change organ initiation",
            phenotype_relation="specific fasciation hypothesis; trait support is unverified",
            uncertainty_note="review seed only; no direct phenotype evidence was supplied",
        )
        normalized = normalize_semantic_claims(pd.DataFrame([proposed], columns=CLAIM_COLUMNS))
        self.assertEqual("proposed", normalized.iloc[0]["verification_status"])
        self.assertTrue(narrowing_eligible_claims(normalized).empty)

        with self.assertRaisesRegex(ValueError, "cannot be verified support"):
            normalize_semantic_claims(pd.DataFrame([
                {**proposed, "verification_status": "verified", "source_id": "PMID:123"}
            ], columns=CLAIM_COLUMNS))

        for column in ("mechanism", "phenotype_relation", "uncertainty_note"):
            with self.subTest(column=column), self.assertRaisesRegex(ValueError, column):
                normalize_semantic_claims(pd.DataFrame([
                    {**proposed, column: ""}
                ], columns=CLAIM_COLUMNS))

    def test_semantic_claims_can_never_write_a_score(self) -> None:
        with self.assertRaisesRegex(ValueError, "rank_contribution"):
            normalize_semantic_claims(
                pd.DataFrame([_claim(rank_contribution="1.0")], columns=CLAIM_COLUMNS)
            )

    def test_unclassified_interval_gene_remains_representable_without_family(self) -> None:
        row = _claim(
            family_id="",
            reference_gene="",
            member_interval_genes="DARK_GENE_1",
            evidence_scope="interval_gene",
            claim_type="annotation_gap_hypothesis",
        )
        normalized = normalize_semantic_claims(pd.DataFrame([row], columns=CLAIM_COLUMNS))
        self.assertEqual("", normalized.iloc[0]["family_id"])
        self.assertEqual("DARK_GENE_1", normalized.iloc[0]["member_interval_genes"])

    def test_orthogroup_scope_has_a_distinct_anchor(self) -> None:
        row = _claim(
            family_id="",
            orthogroup_id="OG0001234",
            reference_gene="",
            member_interval_genes="",
            evidence_scope="orthogroup",
        )
        normalized = normalize_semantic_claims(pd.DataFrame([row], columns=CLAIM_COLUMNS))
        self.assertEqual("OG0001234", normalized.iloc[0]["orthogroup_id"])

    def test_unanchored_claim_fails_closed(self) -> None:
        row = _claim(
            family_id="",
            orthogroup_id="",
            member_interval_genes="",
            reference_gene="",
            evidence_scope="unknown",
        )
        with self.assertRaisesRegex(ValueError, "biological anchor"):
            normalize_semantic_claims(pd.DataFrame([row], columns=CLAIM_COLUMNS))

    def test_scope_specific_anchor_is_required(self) -> None:
        row = _claim(orthogroup_id="", evidence_scope="orthogroup")
        with self.assertRaisesRegex(ValueError, "orthogroup_id"):
            normalize_semantic_claims(pd.DataFrame([row], columns=CLAIM_COLUMNS))

    def test_gene_lists_and_artifact_order_are_canonical(self) -> None:
        rows = [
            _claim(claim_id="b", member_interval_genes="G_B;G_A;G_B"),
            _claim(claim_id="a", member_interval_genes="G_D,G_C"),
        ]
        normalized = normalize_semantic_claims(pd.DataFrame(rows, columns=CLAIM_COLUMNS))
        self.assertEqual(["a", "b"], normalized["claim_id"].tolist())
        self.assertEqual("G_C;G_D", normalized.iloc[0]["member_interval_genes"])
        self.assertEqual("G_A;G_B", normalized.iloc[1]["member_interval_genes"])

    def test_content_fingerprint_does_not_change_when_only_verification_state_changes(self) -> None:
        proposed = _claim(verification_status="proposed", verification_reason="")
        verified = dict(proposed)
        verified.update({
            "verification_status": "verified",
            "verification_reason": "state-only fixture",
        })
        self.assertEqual(claim_content_fingerprint(proposed), claim_content_fingerprint(verified))

    def test_freeze_writes_auditable_status_counts_and_hash(self) -> None:
        frame = pd.DataFrame([
            _claim(claim_id="verified"),
            _claim(claim_id="unsupported", verification_status="unsupported", source_id="",
                   verification_reason="source does not support relation"),
            _claim(claim_id="missing", verification_status="not_evaluated", source_id="",
                   verification_reason="provider not queried"),
        ], columns=CLAIM_COLUMNS)
        with tempfile.TemporaryDirectory() as tmp:
            result = freeze_semantic_claims(frame, tmp, case_id="case403")
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual(3, manifest["claim_count"])
            self.assertEqual(1, manifest["semantic_unit_count"])
            self.assertEqual(1, manifest["verified_claim_count"])
            self.assertEqual(1, manifest["verification_status_counts"]["not_evaluated"])
            self.assertEqual(64, len(manifest["artifact_sha256"]))
            self.assertTrue(all("score" not in key for key in manifest if key != "policy"))

    def test_cli_freeze_is_reachable(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            claims = root / "claims.tsv"
            pd.DataFrame([_claim()], columns=CLAIM_COLUMNS).to_csv(
                claims, sep="\t", index=False, lineterminator="\n"
            )
            self.assertEqual(
                0,
                main([
                    "semantic-claims-freeze",
                    "--claims", str(claims),
                    "--out-dir", str(root / "out"),
                    "--case-id", "cli403",
                ]),
            )
            manifest = json.loads((root / "out" / "semantic_claims_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("cli403", manifest["case_id"])
            self.assertEqual(1, manifest["verified_claim_count"])


if __name__ == "__main__":
    unittest.main()
