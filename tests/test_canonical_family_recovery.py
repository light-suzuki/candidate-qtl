from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.canonical_family_recovery import build_canonical_family_recovery


UNION_COLUMNS = [
    "candidate_gene_id",
    "route_a_supported",
    "route_b_supported",
    "route_support",
    "route_a_family_ids",
    "route_a_hypothesis_ids",
    "route_a_mapping_tiers",
    "route_a_ambiguity_statuses",
    "route_b_family_ids",
    "route_b_orthogroup_ids",
    "route_b_semantic_unit_ids",
    "route_b_claim_ids",
    "route_b_evidence_scopes",
    "support_status",
    "rank_contribution",
]


def _union(path: Path, rows: list[dict[str, str]]) -> Path:
    frame = pd.DataFrame(rows)
    for column in UNION_COLUMNS:
        if column not in frame.columns:
            frame[column] = ""
    frame[UNION_COLUMNS].to_csv(path, sep="\t", index=False, lineterminator="\n")
    return path


def _row(
    gene: str,
    *,
    a: bool = False,
    b: bool = False,
    a_family: str = "",
    b_family: str = "",
    b_og: str = "",
    a_tier: str = "tier1_source_backed_unique",
) -> dict[str, str]:
    return {
        "candidate_gene_id": gene,
        "route_a_supported": "true" if a else "false",
        "route_b_supported": "true" if b else "false",
        "route_support": "a+b" if a and b else ("a" if a else "b"),
        "route_a_family_ids": a_family,
        "route_a_hypothesis_ids": "ha" if a else "",
        "route_a_mapping_tiers": a_tier if a else "",
        "route_a_ambiguity_statuses": "unique" if a else "",
        "route_b_family_ids": b_family,
        "route_b_orthogroup_ids": b_og,
        "route_b_semantic_unit_ids": "ub" if b else "",
        "route_b_claim_ids": "cb" if b else "",
        "route_b_evidence_scopes": "family" if b else "",
        "support_status": "supported_member",
        "rank_contribution": "none",
    }


def _bridge(path: Path, rows: list[dict[str, str]]) -> Path:
    pd.DataFrame(rows, columns=[
        "gene_id", "family_id", "route", "mapping_tier",
        "dependency_cluster_id", "ambiguity_status",
    ]).to_csv(path, sep="\t", index=False, lineterminator="\n")
    return path


class CanonicalFamilyRecoveryTests(unittest.TestCase):
    def test_exact_family_identity_merges_only_supported_members(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            union = _union(root / "union.tsv", [
                _row("G1", a=True, a_family="CLV"),
                _row("G2", b=True, b_family="CLV", b_og="OG_CLV"),
            ])
            result = build_canonical_family_recovery(str(union), str(root / "out"), case_id="exact")
            families = pd.read_csv(result["family_recovery"], sep="\t", dtype=str).fillna("")
            self.assertEqual(1, len(families))
            row = families.iloc[0]
            self.assertEqual("CLV", row["family_id"])
            self.assertEqual("G1;G2", row["member_interval_genes"])
            self.assertEqual("a+b", row["route_support"])
            self.assertEqual("exact_family_label", row["family_match_basis"])

            comparison = pd.read_csv(result["comparison"], sep="\t", dtype=str).fillna("")
            self.assertEqual("a_and_b_supported", comparison.iloc[0]["four_state"])
            self.assertEqual("not_evaluated", comparison.iloc[0]["source_dependence_status"])

    def test_unique_expanded_match_recovers_descriptive_a_vs_symbol_b(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            union = _union(root / "union.tsv", [
                _row("A1", a=True, a_family="RLK LRR receptor kinases"),
                _row("B1", b=True, b_family="RLK1", b_og="OG_RLK"),
            ])
            result = build_canonical_family_recovery(str(union), str(root / "out"), case_id="expanded")
            families = pd.read_csv(result["family_recovery"], sep="\t", dtype=str).fillna("")
            self.assertEqual(1, len(families))
            row = families.iloc[0]
            self.assertEqual("a+b", row["route_support"])
            self.assertEqual("expanded_symbol_token", row["family_match_basis"])
            self.assertEqual("A1;B1", row["member_interval_genes"])
            self.assertEqual("none", row["rank_contribution"])

    def test_ambiguous_family_key_does_not_force_merge(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            union = _union(root / "union.tsv", [
                _row("A1", a=True, a_family="RLK LRR receptor kinases"),
                _row("B1", b=True, b_family="RLK1", b_og="OG1"),
                _row("B2", b=True, b_family="RLK2", b_og="OG2"),
            ])
            result = build_canonical_family_recovery(str(union), str(root / "out"), case_id="ambiguous")
            families = pd.read_csv(result["family_recovery"], sep="\t", dtype=str).fillna("")
            self.assertEqual(3, len(families))
            self.assertNotIn("a+b", set(families["route_support"]))
            a_row = families[families["route_support"].eq("a")].iloc[0]
            self.assertEqual("identity_ambiguous_unmerged", a_row["family_match_basis"])
            self.assertEqual("yes", a_row["review_required"])
            self.assertTrue(
                (families[families["route_support"].eq("b")]["review_required"] == "yes").all()
            )

    def test_b_orthogroup_only_is_retained_without_fabricating_family_label(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            union = _union(root / "union.tsv", [
                _row("G1", b=True, b_family="", b_og="OG_DARK"),
            ])
            result = build_canonical_family_recovery(str(union), str(root / "out"), case_id="og")
            families = pd.read_csv(result["family_recovery"], sep="\t", dtype=str).fillna("")
            self.assertEqual("orthogroup:OG_DARK", families.iloc[0]["family_id"])
            self.assertEqual("route_b_orthogroup_only", families.iloc[0]["family_match_basis"])
            self.assertEqual("G1", families.iloc[0]["member_interval_genes"])

    def test_family_members_are_always_subset_of_supported_union(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            union = _union(root / "union.tsv", [
                _row("G1", a=True, a_family="FAM"),
                _row("G2", b=True, b_family="FAM"),
            ])
            result = build_canonical_family_recovery(str(union), str(root / "out"), case_id="subset")
            families = pd.read_csv(result["family_recovery"], sep="\t", dtype=str).fillna("")
            observed = {
                token
                for value in families["member_interval_genes"]
                for token in str(value).split(";")
                if token
            }
            self.assertEqual({"G1", "G2"}, observed)
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertTrue(manifest["family_members_subset_of_supported_union"])
            self.assertFalse(manifest["legacy_triple_route_family_context_consumed"])

    def test_bridge_projection_preserves_a_only_tier_but_never_invents_b_tier(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            union = _union(root / "union.tsv", [
                _row("A_ONLY", a=True, a_family="PXY", a_tier="tier5_name_only"),
                _row("AB", a=True, b=True, a_family="CLV", b_family="CLV"),
                _row("B_ONLY", b=True, b_family="WUS"),
            ])
            bridge = _bridge(root / "bridge.tsv", [
                {
                    "gene_id": "A_ONLY", "family_id": "PXY", "route": "a",
                    "mapping_tier": "tier5_name_only", "dependency_cluster_id": "dep-a",
                    "ambiguity_status": "name_only",
                },
                {
                    "gene_id": "AB", "family_id": "CLV", "route": "a",
                    "mapping_tier": "tier5_name_only", "dependency_cluster_id": "dep-ab",
                    "ambiguity_status": "name_only",
                },
                {
                    "gene_id": "AB", "family_id": "CLV", "route": "b",
                    "mapping_tier": "tier1_source_backed_unique", "dependency_cluster_id": "dep-b",
                    "ambiguity_status": "unique",
                },
            ])
            result = build_canonical_family_recovery(
                str(union), str(root / "out"), legacy_bridge_arg=str(bridge), case_id="bridge"
            )
            projected = pd.read_csv(result["family_bridge"], sep="\t", dtype=str).fillna("")
            self.assertEqual(1, len(projected))
            self.assertEqual("A_ONLY", projected.iloc[0]["gene_id"])
            self.assertEqual("PXY", projected.iloc[0]["family_id"])
            self.assertEqual("tier5_name_only", projected.iloc[0]["mapping_tier"])
            self.assertEqual({"a"}, set(projected["route"]))


if __name__ == "__main__":
    unittest.main()
