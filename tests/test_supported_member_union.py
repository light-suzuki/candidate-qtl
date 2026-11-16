from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.supported_member_union import build_supported_member_union
from tests._issue664_helpers import manifest_for


class SupportedMemberUnionTests(unittest.TestCase):
    def _bridge(self, path: Path) -> None:
        pd.DataFrame([
            {
                "gene_id": "G1", "route": "a", "mapping_tier": "tier1_source_backed_unique",
                "provider_status": "complete", "review_status": "verified", "family_id": "FA",
                "hypothesis_id": "HA", "ambiguity_status": "resolved_single_copy",
            },
            {
                "gene_id": "G2", "route": "a", "mapping_tier": "tier2_source_backed_ambiguous",
                "provider_status": "complete", "review_status": "review_required", "family_id": "FA",
                "hypothesis_id": "HA", "ambiguity_status": "paralog_ambiguous",
            },
            {
                "gene_id": "G3", "route": "a", "mapping_tier": "tier6_not_evaluated",
                "provider_status": "not_run", "review_status": "not_evaluated", "family_id": "FB",
                "hypothesis_id": "HB", "ambiguity_status": "not_evaluated",
            },
            {
                "gene_id": "G9", "route": "b", "mapping_tier": "tier1_source_backed_unique",
                "provider_status": "complete", "review_status": "verified", "family_id": "OTHER",
                "hypothesis_id": "", "ambiguity_status": "resolved_single_copy",
            },
        ]).to_csv(path, sep="\t", index=False, lineterminator="\n")
        manifest_for(path, route="a", case_id="case426")

    def _claims(self, path: Path) -> None:
        pd.DataFrame([
            {
                "claim_id": "CB1", "verification_status": "verified",
                "member_interval_genes": "G2;G4", "family_id": "EVB", "orthogroup_id": "OGB",
                "semantic_unit_id": "UNITB", "evidence_scope": "family",
            },
            {
                "claim_id": "CB2", "verification_status": "not_evaluated",
                "member_interval_genes": "G5", "family_id": "EVC", "orthogroup_id": "OGC",
                "semantic_unit_id": "UNITC", "evidence_scope": "family",
            },
        ]).to_csv(path, sep="\t", index=False, lineterminator="\n")
        manifest_for(path, route="b_claims", case_id="case426")

    def test_union_keeps_ambiguous_supported_route_a_and_verified_route_b_members(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bridge = root / "bridge.tsv"
            claims = root / "claims.tsv"
            self._bridge(bridge)
            self._claims(claims)
            result = build_supported_member_union(
                str(root / "out"), route_a_bridge_arg=str(bridge), route_b_claims_arg=str(claims), route_a_manifest_arg=str(root / "bridge_manifest.json"), route_b_claims_manifest_arg=str(root / "claims_manifest.json"), case_id="case426"
            )
            frame = pd.read_csv(result["candidate_set"], sep="\t", dtype=str).fillna("")
            self.assertEqual({"G1", "G2", "G4"}, set(frame["candidate_gene_id"]))
            self.assertNotIn("G3", set(frame["candidate_gene_id"]))
            self.assertNotIn("G5", set(frame["candidate_gene_id"]))
            self.assertNotIn("G9", set(frame["candidate_gene_id"]))

            g2 = frame.set_index("candidate_gene_id").loc["G2"]
            self.assertEqual("True", g2["route_a_supported"])
            self.assertEqual("True", g2["route_b_supported"])
            self.assertEqual("a+b", g2["route_support"])
            self.assertIn("paralog_ambiguous", g2["route_a_ambiguity_statuses"])
            self.assertEqual("CB1", g2["route_b_claim_ids"])
            self.assertEqual("none", g2["rank_contribution"])

            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual(3, manifest["candidate_count"])
            self.assertEqual(2, manifest["route_a_candidate_count"])
            self.assertEqual(2, manifest["route_b_candidate_count"])
            self.assertEqual(1, manifest["route_agreement_candidate_count"])
            self.assertEqual(1, manifest["route_a_not_evaluated_rows"])
            self.assertEqual(1, manifest["route_b_nonverified_claim_rows"])

    def test_route_agreement_is_provenance_not_score(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bridge = root / "bridge.tsv"
            claims = root / "claims.tsv"
            self._bridge(bridge)
            self._claims(claims)
            result = build_supported_member_union(
                str(root / "out"), route_a_bridge_arg=str(bridge), route_b_claims_arg=str(claims), route_a_manifest_arg=str(root / "bridge_manifest.json"), route_b_claims_manifest_arg=str(root / "claims_manifest.json")
            )
            frame = pd.read_csv(result["candidate_set"], sep="\t", dtype=str).fillna("")
            self.assertFalse(any("score" in column.casefold() for column in frame.columns))
            self.assertTrue(frame["rank_contribution"].eq("none").all())

    def test_route_a_only_and_route_b_only_are_supported(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bridge = root / "bridge.tsv"
            claims = root / "claims.tsv"
            self._bridge(bridge)
            self._claims(claims)
            a = build_supported_member_union(str(root / "a"), route_a_bridge_arg=str(bridge), route_a_manifest_arg=str(root / "bridge_manifest.json"))
            b = build_supported_member_union(str(root / "b"), route_b_claims_arg=str(claims), route_b_claims_manifest_arg=str(root / "claims_manifest.json"))
            self.assertEqual({"G1", "G2"}, set(pd.read_csv(a["candidate_set"], sep="\t")["candidate_gene_id"]))
            self.assertEqual({"G2", "G4"}, set(pd.read_csv(b["candidate_set"], sep="\t")["candidate_gene_id"]))

    def test_requires_at_least_one_frozen_route_artifact(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "requires Route A bridge"):
                build_supported_member_union(str(Path(tmp) / "out"))


if __name__ == "__main__":
    unittest.main()
