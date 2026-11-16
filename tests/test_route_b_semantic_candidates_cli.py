from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.semantic_claims import CLAIM_COLUMNS


class RouteBSemanticCandidatesCliTests(unittest.TestCase):
    def test_cli_freezes_route_b_verified_candidates(self) -> None:
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
                "output_fingerprint": "o" * 64, "verification_status": "verified",
                "verification_reason": "group verified", "rank_contribution": "none",
            })
            claims_file = root / "semantic_claims.tsv"
            pd.DataFrame([claim], columns=CLAIM_COLUMNS).to_csv(claims_file, sep="\t", index=False, lineterminator="\n")

            self.assertEqual(0, main([
                "route-b-semantic-candidates",
                "--claims", str(claims_file),
                "--out-dir", str(root / "out"),
                "--case-id", "cli-cand",
            ]))
            cand_path = root / "out" / "route_b_verified_candidates.tsv"
            self.assertTrue(cand_path.exists())
            cand = pd.read_csv(cand_path, sep="\t", dtype=str).fillna("")
            self.assertEqual(["G1", "G2"], cand["candidate_gene_id"].tolist())
            self.assertEqual(["b", "b"], cand["route"].tolist())
            self.assertEqual(["none", "none"], cand["rank_contribution"].tolist())
            manifest = json.loads((root / "out" / "route_b_verified_candidates_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(2, manifest["candidate_gene_count"])
            self.assertEqual(1, manifest["verified_claim_count"])


if __name__ == "__main__":
    unittest.main()
