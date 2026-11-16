from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from tests._issue664_helpers import manifest_for


class SupportedMemberUnionCliTests(unittest.TestCase):
    def test_cli_builds_union_ready_for_candidate_budget_gate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bridge = root / "bridge.tsv"
            pd.DataFrame([{
                "gene_id": "G1", "route": "a", "mapping_tier": "tier1_source_backed_unique",
                "provider_status": "complete", "review_status": "verified", "family_id": "FA",
                "hypothesis_id": "HA", "ambiguity_status": "resolved_single_copy",
            }]).to_csv(bridge, sep="\t", index=False, lineterminator="\n")
            bridge_manifest = manifest_for(bridge, route="a", case_id="cli426")
            claims = root / "claims.tsv"
            pd.DataFrame([{
                "claim_id": "CB", "verification_status": "verified",
                "member_interval_genes": "G1;G2", "family_id": "FB", "orthogroup_id": "OGB",
                "semantic_unit_id": "UNITB", "evidence_scope": "family",
            }]).to_csv(claims, sep="\t", index=False, lineterminator="\n")
            claims_manifest = manifest_for(claims, route="b_claims", case_id="cli426")

            self.assertEqual(0, main([
                "supported-member-union",
                "--route-a-bridge", str(bridge),
                "--route-a-manifest", str(bridge_manifest),
                "--route-b-claims", str(claims),
                "--route-b-claims-manifest", str(claims_manifest),
                "--out-dir", str(root / "union"),
                "--case-id", "cli426",
            ]))
            union = root / "union" / "supported_member_union.tsv"
            frame = pd.read_csv(union, sep="\t", dtype=str).fillna("")
            self.assertEqual({"G1", "G2"}, set(frame["candidate_gene_id"]))
            self.assertEqual("a+b", frame.set_index("candidate_gene_id").loc["G1", "route_support"])

            self.assertEqual(0, main([
                "candidate-budget-gate",
                "--candidates", str(union),
                "--out-dir", str(root / "gate"),
                "--wet-budget", "20",
                "--case-id", "cli426",
            ]))
            manifest = json.loads((root / "gate" / "candidate_budget_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("within_wet_budget", manifest["status"])
            self.assertEqual(2, manifest["candidate_count"])


if __name__ == "__main__":
    unittest.main()
