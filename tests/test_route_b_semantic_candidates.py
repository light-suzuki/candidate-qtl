from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.route_b_semantic_candidates import freeze_route_b_semantic_candidates
from candidate_qtl.semantic_claims import CLAIM_COLUMNS


def _claim(claim_id: str, status: str, members: str, **overrides: str) -> dict[str, str]:
    row = {column: "" for column in CLAIM_COLUMNS}
    row.update({
        "claim_id": claim_id,
        "route": "interval_first",
        "semantic_unit_id": f"unit-{claim_id}",
        "family_id": "EV1",
        "orthogroup_id": "OG1",
        "subfamily_or_clade": "SUB1",
        "member_interval_genes": members,
        "claim_type": "mutant_phenotype",
        "biological_function": "meristem maintenance",
        "mechanism": "receptor signaling",
        "phenotype_relation": "fasciation",
        "reference_gene": "REF1",
        "reference_species": "Arabidopsis thaliana",
        "evidence_scope": "family",
        "source_id": "PMID:1",
        "provider": "fixture",
        "model": "fixture",
        "prompt_fingerprint": "p" * 64,
        "input_fingerprint": "i" * 64,
        "output_fingerprint": "o" * 64,
        "verification_status": status,
        "verification_reason": "fixture",
        "rank_contribution": "none",
    })
    row.update(overrides)
    return row


class RouteBSemanticCandidateFreezeTests(unittest.TestCase):
    def test_only_verified_claims_contribute_and_all_group_members_are_retained(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            claims = root / "claims.tsv"
            pd.DataFrame([
                _claim("v1", "verified", "G1;G2"),
                _claim("v2", "verified", "G2;G3", family_id="EV2", orthogroup_id="OG2"),
                _claim("n1", "not_evaluated", "G4"),
                _claim("u1", "unsupported", "G5"),
            ], columns=CLAIM_COLUMNS).to_csv(claims, sep="\t", index=False, lineterminator="\n")

            result = freeze_route_b_semantic_candidates(str(claims), str(root / "out"), case_id="case403")
            frame = pd.read_csv(result["candidates"], sep="\t", dtype=str).fillna("")
            self.assertEqual({"G1", "G2", "G3"}, set(frame["candidate_gene_id"]))
            self.assertNotIn("G4", set(frame["candidate_gene_id"]))
            self.assertNotIn("G5", set(frame["candidate_gene_id"]))
            g2 = frame.set_index("candidate_gene_id").loc["G2"]
            self.assertEqual("v1;v2", g2["claim_ids"])
            self.assertEqual("EV1;EV2", g2["family_ids"])
            self.assertEqual("b", g2["route"])
            self.assertEqual("none", g2["rank_contribution"])

            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual("interval_first", manifest["route"])
            self.assertEqual(2, manifest["verified_claim_count"])
            self.assertEqual(3, manifest["candidate_gene_count"])
            self.assertEqual(64, len(manifest["artifact_sha256"]))

    def test_no_verified_claims_produces_explicit_empty_route_b_set(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            claims = root / "claims.tsv"
            pd.DataFrame([
                _claim("n1", "not_evaluated", "G1"),
            ], columns=CLAIM_COLUMNS).to_csv(claims, sep="\t", index=False, lineterminator="\n")
            result = freeze_route_b_semantic_candidates(str(claims), str(root / "out"))
            frame = pd.read_csv(result["candidates"], sep="\t", dtype=str).fillna("")
            self.assertEqual(0, len(frame))
            self.assertEqual(0, result["candidate_gene_count"])

    def test_non_route_b_claims_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            claims = root / "claims.tsv"
            row = _claim("a1", "verified", "G1")
            row["route"] = "phenotype_first"
            row["semantic_unit_id"] = ""
            pd.DataFrame([row], columns=CLAIM_COLUMNS).to_csv(claims, sep="\t", index=False, lineterminator="\n")
            with self.assertRaisesRegex(ValueError, "interval_first"):
                freeze_route_b_semantic_candidates(str(claims), str(root / "out"))


if __name__ == "__main__":
    unittest.main()
