from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.semantic_claims import CLAIM_COLUMNS


class RouteBSemanticSourceVerifyCliTests(unittest.TestCase):
    def test_cli_keeps_source_support_separate_from_transfer_verification(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            row = {column: "" for column in CLAIM_COLUMNS}
            row.update({
                "claim_id": "c1", "route": "interval_first", "semantic_unit_id": "bsem-cli",
                "family_id": "EV1", "member_interval_genes": "G1", "claim_type": "mutant_phenotype",
                "biological_function": "meristem maintenance", "mechanism": "receptor signaling",
                "phenotype_relation": "fasciated shoot meristem", "reference_gene": "REF1",
                "reference_species": "Arabidopsis thaliana", "evidence_scope": "reference_gene",
                "source_id": "PMID:1234567", "provider": "fixture", "model": "fixture",
                "prompt_fingerprint": "p" * 64, "input_fingerprint": "i" * 64,
                "output_fingerprint": "o" * 64, "verification_status": "proposed",
                "verification_reason": "source relation pending", "rank_contribution": "none",
            })
            proposals = root / "proposals.tsv"
            pd.DataFrame([row], columns=CLAIM_COLUMNS).to_csv(
                proposals, sep="\t", index=False, lineterminator="\n"
            )
            cache = root / "cache"
            cache.mkdir()
            (cache / "med_1234567.json").write_text(json.dumps({
                "cache_key": "pmid_1234567", "id": "1234567", "pmid": "1234567",
                "title": "REF1 controls the shoot meristem",
                "abstractText": "Loss-of-function knockout mutants of REF1 display a fasciated shoot meristem.",
                "journalTitle": "Plant Journal", "pubYear": "2026", "doi": "10.1000/ref1",
                "authorString": "Example A",
            }), encoding="utf-8")

            self.assertEqual(0, main([
                "route-b-semantic-source-verify",
                "--proposals", str(proposals),
                "--cache", str(cache),
                "--out-dir", str(root / "out"),
                "--case-id", "cli-source",
            ]))
            frame = pd.read_csv(
                root / "out" / "semantic_claim_source_verification.tsv", sep="\t", dtype=str
            ).fillna("")
            self.assertEqual("supported_mechanistic", frame.iloc[0]["source_relation_status"])
            self.assertEqual("proposed", frame.iloc[0]["final_semantic_verification_status"])
            manifest = json.loads(
                (root / "out" / "semantic_claim_source_verification_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual("cli-source", manifest["case_id"])
            self.assertEqual(0, manifest["final_verified_claim_count"])
            self.assertEqual(1, manifest["transfer_verification_required_count"])


if __name__ == "__main__":
    unittest.main()
