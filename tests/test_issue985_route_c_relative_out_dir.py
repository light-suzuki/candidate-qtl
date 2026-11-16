"""#985: a cwd-relative ``--out-dir`` for projection must not break members."""

from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.route_c_direct_members import build_route_c_direct_members
from candidate_qtl.route_c_direct_pipeline import direct_rescue_project

import test_issue_663_route_c_projection_evidence as projection_tests


class RouteCRelativeOutDirTests(unittest.TestCase):
    def test_members_reads_evidence_recorded_relative_to_cwd(self) -> None:
        helper = projection_tests.RouteCProjectionEvidenceTests()
        cwd = os.getcwd()
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            frozen, freeze = helper._frozen(root, include_binding=True)
            pd.DataFrame([{"gene_id": "UniProt:A0A000"}]).to_csv(
                root / "atlas.tsv", sep="\t", index=False, lineterminator="\n"
            )
            pd.DataFrame([helper._hit("UniProt:A0A000")]).to_csv(
                root / "hits.tsv", sep="\t", index=False, lineterminator="\n"
            )
            os.chdir(root)
            try:
                result = direct_rescue_project(
                    frozen.name, freeze.name, "atlas.tsv", "q1", "c-project", "hits.tsv"
                )
                self.assertFalse(Path(result["projection_evidence"]).is_absolute())
                members = build_route_c_direct_members(
                    frozen.name, freeze.name, result["projection"], result["manifest"],
                    "atlas.tsv", "c-members", case_id="case1",
                )
            finally:
                os.chdir(cwd)
            self.assertEqual(1, members["review_count"])

    def test_tampered_evidence_is_still_rejected_after_fallback(self) -> None:
        helper = projection_tests.RouteCProjectionEvidenceTests()
        cwd = os.getcwd()
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            frozen, freeze = helper._frozen(root, include_binding=True)
            pd.DataFrame([{"gene_id": "UniProt:A0A000"}]).to_csv(
                root / "atlas.tsv", sep="\t", index=False, lineterminator="\n"
            )
            pd.DataFrame([helper._hit("UniProt:A0A000")]).to_csv(
                root / "hits.tsv", sep="\t", index=False, lineterminator="\n"
            )
            os.chdir(root)
            try:
                result = direct_rescue_project(
                    frozen.name, freeze.name, "atlas.tsv", "q1", "c-project", "hits.tsv"
                )
                with open(result["projection_evidence"], "a", encoding="utf-8") as handle:
                    handle.write("\n")
                with self.assertRaisesRegex(ValueError, "hash does not match"):
                    build_route_c_direct_members(
                        frozen.name, freeze.name, result["projection"], result["manifest"],
                        "atlas.tsv", "c-members", case_id="case1",
                    )
            finally:
                os.chdir(cwd)


if __name__ == "__main__":
    unittest.main()
