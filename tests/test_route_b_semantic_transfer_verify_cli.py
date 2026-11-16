from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.claim_binding import BOUND_CLAIM_COLUMNS
from candidate_qtl.semantic_claims import CLAIM_COLUMNS
from tests.bound_claim_fixtures import bound_claim


class RouteBSemanticTransferVerifyCliTests(unittest.TestCase):
    def test_cli_transfer_verifies_claims(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            claim = {column: "" for column in CLAIM_COLUMNS}
            claim.update({
                "claim_id": "c1", "route": "interval_first", "semantic_unit_id": "bsem-cli",
                "family_id": "EV1", "orthogroup_id": "OG1", "subfamily_or_clade": "SUB1",
                "member_interval_genes": "G1;G2", "claim_type": "mutant_phenotype",
                "biological_function": "meristem maintenance", "mechanism": "receptor signaling",
                "phenotype_relation": "fasciated shoot meristem", "reference_gene": "REF1",
                "reference_species": "Arabidopsis thaliana", "evidence_scope": "family",
                "source_id": "PMID:1234567", "provider": "fixture", "model": "fixture",
                "prompt_fingerprint": "p" * 64, "input_fingerprint": "i" * 64,
                "output_fingerprint": "o" * 64, "verification_status": "proposed",
                "verification_reason": "source relation pending", "rank_contribution": "none",
            })
            proposals = root / "proposals.tsv"
            pd.DataFrame([claim], columns=CLAIM_COLUMNS).to_csv(proposals, sep="\t", index=False, lineterminator="\n")

            source = root / "source.tsv"
            pd.DataFrame([{
                "claim_id": "c1", "source_relation_status": "supported_mechanistic",
            }]).to_csv(source, sep="\t", index=False, lineterminator="\n")

            bundle = root / "bundle.json"
            bundle.write_text(json.dumps({
                "schema_version": "1.0",
                "case_id": "cli-transfer",
                "route": "interval_first",
                "units": [{
                    "unit_id": "bsem-cli", "unit_type": "family", "unit_anchor": "EV1",
                    "member_interval_genes": ["G1", "G2"],
                    "identity": {
                        "evolutionary_family_ids": ["EV1"], "family_ids": ["F1"],
                        "orthogroup_ids": ["OG1"], "subfamily_ids": ["SUB1"],
                    },
                    "genes": [
                        {"gene_id": "G1", "orthology_reference_gene": "REF1", "orthology_relationship": "1:many"},
                        {"gene_id": "G2", "orthology_reference_gene": "REF1", "orthology_relationship": "1:many"},
                    ],
                    "literature_link_candidates": [{"reference_gene": "REF1", "source_id": "PMID:1234567"}],
                    "stable_source_ids": ["PMID:1234567"], "rank_contribution": "none",
                }],
            }), encoding="utf-8")
            bindings = root / "bound.tsv"
            pd.DataFrame([bound_claim(
                route="b", claim_id="c1", citation="PMID:1234567",
                subject_gene="REF1", organism="Arabidopsis thaliana",
                phenotype="fasciated shoot meristem",
            )], columns=BOUND_CLAIM_COLUMNS).to_csv(
                bindings, sep="\t", index=False, lineterminator="\n"
            )

            self.assertEqual(0, main([
                "route-b-semantic-transfer-verify",
                "--proposals", str(proposals),
                "--source-verification", str(source),
                "--bundle", str(bundle),
                "--bound-claims", str(bindings),
                "--out-dir", str(root / "out"),
                "--case-id", "cli-transfer",
            ]))
            claims = pd.read_csv(root / "out" / "semantic_claims_verified.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("verified", claims.iloc[0]["verification_status"])
            manifest = json.loads((root / "out" / "semantic_claim_transfer_verification_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(1, manifest["verified_claim_count"])


if __name__ == "__main__":
    unittest.main()
