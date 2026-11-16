from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.supported_member_union import build_supported_member_union
from tests._issue664_helpers import manifest_for


class SupportedMemberUnionFrozenRouteBTests(unittest.TestCase):
    def test_precomparison_route_b_artifact_is_consumed_without_reopening_claims(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_b = root / "route_b_verified_candidates.tsv"
            pd.DataFrame([
                {
                    "candidate_gene_id": "G2",
                    "family_ids": "EV1",
                    "orthogroup_ids": "OG1",
                    "subfamilies": "SUB1",
                    "semantic_unit_ids": "UNIT1",
                    "claim_ids": "C1",
                    "evidence_scopes": "family",
                    "reference_genes": "REF1",
                    "source_ids": "PMID:1",
                    "route": "b",
                    "support_status": "verified_semantic_support",
                    "rank_contribution": "none",
                },
                {
                    "candidate_gene_id": "G3",
                    "family_ids": "EV1",
                    "orthogroup_ids": "OG1",
                    "subfamilies": "SUB1",
                    "semantic_unit_ids": "UNIT1",
                    "claim_ids": "C1",
                    "evidence_scopes": "family",
                    "reference_genes": "REF1",
                    "source_ids": "PMID:1",
                    "route": "b",
                    "support_status": "verified_semantic_support",
                    "rank_contribution": "none",
                },
            ]).to_csv(route_b, sep="\t", index=False, lineterminator="\n")
            manifest = manifest_for(route_b, route="b", case_id="case426")

            result = build_supported_member_union(
                str(root / "out"), route_b_candidates_arg=str(route_b), route_b_candidates_manifest_arg=str(manifest), case_id="case426"
            )
            frame = pd.read_csv(result["candidate_set"], sep="\t", dtype=str).fillna("")
            self.assertEqual({"G2", "G3"}, set(frame["candidate_gene_id"]))
            self.assertTrue(frame["route_b_supported"].astype(str).eq("True").all())
            self.assertTrue(frame["route_support"].eq("b").all())
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual("precomparison_frozen_candidates", manifest["route_b_input_kind"])
            self.assertEqual(2, manifest["route_b_candidate_count"])

    def test_route_b_frozen_artifact_rejects_wrong_route_or_support_state(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            route_b = root / "bad.tsv"
            pd.DataFrame([{
                "candidate_gene_id": "G1",
                "route": "a",
                "support_status": "verified_semantic_support",
            }]).to_csv(route_b, sep="\t", index=False, lineterminator="\n")
            manifest = manifest_for(route_b, route="b")
            with self.assertRaisesRegex(ValueError, "route=b"):
                build_supported_member_union(str(root / "out"), route_b_candidates_arg=str(route_b), route_b_candidates_manifest_arg=str(manifest))

    def test_route_b_claims_and_frozen_candidates_cannot_both_be_supplied(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            empty = root / "empty.tsv"
            pd.DataFrame().to_csv(empty, sep="\t", index=False)
            manifest_for(empty, route="b")
            with self.assertRaisesRegex(ValueError, "only one Route B input"):
                build_supported_member_union(
                    str(root / "out"),
                    route_b_candidates_arg=str(empty),
                    route_b_claims_arg=str(empty),
                )


if __name__ == "__main__":
    unittest.main()
