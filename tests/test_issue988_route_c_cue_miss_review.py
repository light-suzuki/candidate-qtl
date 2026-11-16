"""#988: a verified C source without a direct-test cue stays review-only, not lost."""
from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.route_c_direct_members import build_route_c_direct_members
from candidate_qtl.route_c_direct_pipeline import (
    DIRECT_PROPOSAL_COLUMNS,
    direct_rescue_freeze,
    direct_rescue_project,
)


class RouteCCueMissReviewTests(unittest.TestCase):
    def _run(self, root: Path, abstract: str) -> tuple[dict, dict, dict]:
        proposal = {column: "" for column in DIRECT_PROPOSAL_COLUMNS}
        proposal.update({
            "record_id": "c1", "query_phenotype": "altered growth timing",
            "reported_phenotype": "altered growth timing",
            "focal_species": "Synthetic species", "reference_species": "Synthetic species",
            "reported_gene_symbol": "GENX1", "species_distance": "S1_same_species",
            "evidence_directness": "D1_mapped_gene_with_functional_validation",
            "evidence_type": "mapped locus", "source_id": "PMID:12345678",
        })
        proposals, cache = root / "proposals.tsv", root / "cache.json"
        pd.DataFrame([proposal], columns=DIRECT_PROPOSAL_COLUMNS).to_csv(proposals, sep="\t", index=False)
        cache.write_text(json.dumps([{
            "id": "12345678", "pmid": "12345678", "title": "A synthetic source",
            "abstractText": abstract, "pubYear": "2020",
        }]), encoding="utf-8")
        frozen = direct_rescue_freeze(str(proposals), str(root / "freeze"), str(cache))
        atlas = root / "atlas.tsv"
        pd.DataFrame([{"gene_id": "GENX1", "chrom": "chr1", "start": "1", "end": "100"}]).to_csv(
            atlas, sep="\t", index=False
        )
        projected = direct_rescue_project(frozen["frozen"], frozen["manifest"], str(atlas), "q1",
                                          str(root / "projection"))
        members = build_route_c_direct_members(frozen["frozen"], frozen["manifest"], projected["projection"],
                                               projected["manifest"], str(atlas), str(root / "members"),
                                               case_id="q1")
        with Path(frozen["frozen"]).open(encoding="utf-8", newline="") as stream:
            row = next(csv.DictReader(stream, delimiter="\t"))
        return row, projected, members

    def test_named_gene_without_cue_projects_as_review_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            row, projected, members = self._run(
                Path(tmp), "GENX1 contributes to regulation of growth timing in Synthetic species."
            )
            self.assertEqual("pending", row["source_verification_state"])
            self.assertEqual("verified_gene", row["gene_identity_state"])
            self.assertEqual("pending_verification", row["direct_rescue_eligibility"])
            self.assertEqual(1, projected["inside_count"])
            self.assertEqual(0, members["supported_count"])
            self.assertEqual(1, members["review_count"])

    def test_source_that_does_not_name_the_gene_is_still_not_projected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            row, projected, members = self._run(
                Path(tmp), "An unrelated locus contributes to regulation of growth timing."
            )
            self.assertEqual("unresolved", row["gene_identity_state"])
            self.assertEqual(0, projected["inside_count"])
            self.assertEqual(0, members["review_count"])


if __name__ == "__main__":
    unittest.main()
