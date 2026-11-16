"""Deterministic subfamily-granularity refinement verdicts (#426 bounded pass).

The verification producer may only falsify support when the evidence chain is
fully machine-checked: curated subfamily assignment + verified target +
target-relative mismatch + Route-A-derived support. Everything else stays
``not_evaluated`` -- an unresolved question never removes a candidate.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.subfamily_refinement import build_subfamily_refinement_verification


def _candidates(root: Path, rows: list[dict[str, str]]) -> Path:
    path = root / "supported_member_union.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, lineterminator="\n")
    return path


def _discrimination(root: Path, rows: list[dict[str, str]]) -> Path:
    path = root / "family_copy_discrimination.tsv"
    base = {
        "family_id": "FAM1", "copy_set_id": "CS1", "copy_set_status": "validated",
        "orthology_class": "", "target_id": "TGT1", "target_reference_gene": "",
        "domain_vs_target": "", "disambiguation_status": "",
    }
    enriched = [{**base, **row} for row in rows]
    pd.DataFrame(enriched).to_csv(path, sep="\t", index=False, lineterminator="\n")
    return path


def _relationships(root: Path, sources: dict[str, str]) -> Path:
    path = root / "family_copy_relationships.tsv"
    rows = [
        {"family_id": "FAM1", "candidate_gene_id": gene, "reference_gene": "",
         "orthology_class": "", "orthology_providers": "", "evidence_chain": "",
         "subfamily": "", "subfamily_source": source, "paralog_count": "", "notes": ""}
        for gene, source in sorted(sources.items())
    ]
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, lineterminator="\n")
    return path


class SubfamilyRefinementVerdictTests(unittest.TestCase):
    def test_curated_cross_subfamily_route_a_only_is_falsified(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = _candidates(root, [
                {"candidate_gene_id": "G001", "route_support": "a"},
                {"candidate_gene_id": "G002", "route_support": "a"},
            ])
            discrimination = _discrimination(root, [
                {"candidate_gene_id": "G001", "subfamily": "CLV1like",
                 "subfamily_vs_target": "target_mismatch",
                 "target_verification_status": "verified"},
                {"candidate_gene_id": "G002", "subfamily": "TFL1like",
                 "subfamily_vs_target": "target_match",
                 "target_verification_status": "verified"},
            ])
            relationships = _relationships(
                root, {"G001": "curated_map", "G002": "curated_map"})
            result = build_subfamily_refinement_verification(
                str(candidates), str(discrimination), str(relationships),
                str(root / "out"), case_id="case1",
            )
            frame = pd.read_csv(result["verification_table"], sep="\t", dtype=str).fillna("")
            by_gene = dict(zip(frame["candidate_gene_id"], frame["verdict"]))
            self.assertEqual("support_falsified", by_gene["G001"])
            self.assertEqual("support_confirmed", by_gene["G002"])
            self.assertEqual(1, result["removed_candidate_count"])
            source_row = frame[frame["candidate_gene_id"] == "G001"].iloc[0]
            self.assertIn("curated_cross_subfamily_vs_verified_target", source_row["verification_source"])
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual("case1", manifest["case_id"])
            self.assertEqual(["G001"], manifest["falsified_candidate_genes"])

    def test_orthogroup_proxy_mismatch_never_falsifies(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = _candidates(root, [
                {"candidate_gene_id": "G001", "route_support": "a"}])
            discrimination = _discrimination(root, [
                {"candidate_gene_id": "G001", "subfamily": "OG_0009",
                 "subfamily_vs_target": "target_mismatch",
                 "target_verification_status": "verified"}])
            relationships = _relationships(root, {"G001": "orthogroup"})
            result = build_subfamily_refinement_verification(
                str(candidates), str(discrimination), str(relationships), str(root / "out"))
            frame = pd.read_csv(result["verification_table"], sep="\t", dtype=str).fillna("")
            self.assertEqual("not_evaluated", frame.iloc[0]["verdict"])
            self.assertEqual(0, result["removed_candidate_count"])

    def test_pending_review_target_never_falsifies(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = _candidates(root, [
                {"candidate_gene_id": "G001", "route_support": "a"}])
            discrimination = _discrimination(root, [
                {"candidate_gene_id": "G001", "subfamily": "OTHER",
                 "subfamily_vs_target": "target_mismatch",
                 "target_verification_status": "pending_review"}])
            relationships = _relationships(root, {"G001": "curated_map"})
            result = build_subfamily_refinement_verification(
                str(candidates), str(discrimination), str(relationships), str(root / "out"))
            frame = pd.read_csv(result["verification_table"], sep="\t", dtype=str).fillna("")
            self.assertEqual("not_evaluated", frame.iloc[0]["verdict"])

    def test_route_b_supported_member_is_out_of_scope(self) -> None:
        """Route B verified support does not stand on the Route A projection."""
        for support in ("b", "a+b"):
            with tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                candidates = _candidates(root, [
                    {"candidate_gene_id": "G001", "route_support": support}])
                discrimination = _discrimination(root, [
                    {"candidate_gene_id": "G001", "subfamily": "OTHER",
                     "subfamily_vs_target": "target_mismatch",
                     "target_verification_status": "verified"}])
                relationships = _relationships(root, {"G001": "curated_map"})
                result = build_subfamily_refinement_verification(
                    str(candidates), str(discrimination), str(relationships), str(root / "out"))
                frame = pd.read_csv(result["verification_table"], sep="\t", dtype=str).fillna("")
                self.assertEqual("not_evaluated", frame.iloc[0]["verdict"])

    def test_missing_discrimination_row_stays_not_evaluated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = _candidates(root, [
                {"candidate_gene_id": "G001", "route_support": "a"},
                {"candidate_gene_id": "G002", "route_support": "a"}])
            discrimination = _discrimination(root, [
                {"candidate_gene_id": "G001", "subfamily": "A",
                 "subfamily_vs_target": "target_match",
                 "target_verification_status": "verified"}])
            result = build_subfamily_refinement_verification(
                str(candidates), str(discrimination), None, str(root / "out"))
            self.assertEqual(1, result["verdict_counts"]["support_confirmed"])
            self.assertEqual(1, result["verdict_counts"]["not_evaluated"])

    def test_unresolved_relation_and_missing_files_fail_clean(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = _candidates(root, [
                {"candidate_gene_id": "G001", "route_support": "a"}])
            discrimination = _discrimination(root, [
                {"candidate_gene_id": "G001", "subfamily": "",
                 "subfamily_vs_target": "subfamily_unresolved",
                 "target_verification_status": ""}])
            result = build_subfamily_refinement_verification(
                str(candidates), str(discrimination), None, str(root / "out"))
            self.assertEqual(1, result["verdict_counts"]["not_evaluated"])

            with self.assertRaises(ValueError):
                build_subfamily_refinement_verification(
                    str(root / "missing.tsv"), str(discrimination), None, str(root / "o2"))
            with self.assertRaises(ValueError):
                build_subfamily_refinement_verification(
                    str(candidates), str(root / "missing.tsv"), None, str(root / "o3"))

    def test_duplicate_ids_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = _candidates(root, [
                {"candidate_gene_id": "G001", "route_support": "a"},
                {"candidate_gene_id": "G001", "route_support": "a"}])
            discrimination = _discrimination(root, [
                {"candidate_gene_id": "G001", "subfamily_vs_target": "subfamily_unresolved"}])
            with self.assertRaises(ValueError):
                build_subfamily_refinement_verification(
                    str(candidates), str(discrimination), None, str(root / "out"))

            candidates_ok = _candidates(root, [
                {"candidate_gene_id": "G001", "route_support": "a"}])
            discrimination_dup = _discrimination(root, [
                {"candidate_gene_id": "G001", "subfamily_vs_target": "subfamily_unresolved"},
                {"candidate_gene_id": "G001", "subfamily_vs_target": "target_match"}])
            with self.assertRaises(ValueError):
                build_subfamily_refinement_verification(
                    str(candidates_ok), str(discrimination_dup), None, str(root / "out2"))


if __name__ == "__main__":
    unittest.main()
