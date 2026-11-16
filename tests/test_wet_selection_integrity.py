"""Negative-path coverage for independent selection artifacts."""
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.adaptive_synthesis import consume_review_inclusive_selection
from candidate_qtl.core import sha256


class WetSelectionIntegrityTests(unittest.TestCase):
    def test_invalid_independent_artifact_cannot_remove_a_candidate(self):
        mutations = {
            "empty_proposal_hash": ("proposal_sha256", ""),
            "wrong_proposal_hash": ("proposal_sha256", "0" * 64),
            "empty_wet_hash": ("source_wet_sha256", ""),
            "wrong_wet_hash": ("source_wet_sha256", "0" * 64),
            "missing_reviewer": ("reviewer_id", ""),
            "wrong_reviewer": ("reviewer_id", "another-agent"),
            "missing_model": ("reviewer_model", ""),
            "empty_reason": ("reason", ""),
            "unknown_gene": ("candidate_gene_id", "alien"),
            "duplicate_gene": ("candidate_gene_id", "G0"),
        }
        for label, (column, value) in mutations.items():
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                wet = root / "wet.tsv"
                pd.DataFrame({"candidate_gene_id": ["G0", "G1"], "route_a_supported": ["False", "True"]}).to_csv(wet, sep="\t", index=False)
                before = wet.read_bytes()
                proposal = root / "proposal.tsv"
                records = [{"candidate_gene_id": gene, "decision": "defer", "reason": "alternative experiment prioritized",
                            "executor_id": "executor", "executor_model": "luna",
                            "source_wet_sha256": sha256(wet)} for gene in ("G0", "G1")]
                pd.DataFrame(records).to_csv(proposal, sep="\t", index=False)
                independent = root / "independent.tsv"
                data = pd.DataFrame([{**record, "proposal_sha256": sha256(proposal),
                                      "reviewer_id": "reviewer", "reviewer_model": "luna"} for record in records])
                data.loc[1, column] = value
                data.to_csv(independent, sep="\t", index=False)
                with self.assertRaises(ValueError):
                    consume_review_inclusive_selection(str(wet), str(proposal), str(root / "out"),
                        independent_review_arg=str(independent), executor_id="executor", executor_model="luna",
                        reviewer_id="reviewer", reviewer_model="luna")
                self.assertEqual(before, wet.read_bytes())
                self.assertFalse((root / "out" / "wet_candidates.tsv").exists())

    def test_deferral_requires_a_distinct_reviewer(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            wet = root / "wet.tsv"
            pd.DataFrame({"candidate_gene_id": ["G1"]}).to_csv(wet, sep="\t", index=False)
            proposal = root / "proposal.tsv"
            pd.DataFrame([{"candidate_gene_id": "G1", "decision": "defer", "reason": "defer",
                           "source_wet_sha256": sha256(wet), "executor_id": "executor",
                           "executor_model": "luna"}]).to_csv(proposal, sep="\t", index=False)
            for executor, reviewer in (("", "reviewer"), ("executor", "reviewer"), ("same", "same")):
                with self.subTest(executor=executor, reviewer=reviewer), self.assertRaises(ValueError):
                    consume_review_inclusive_selection(str(wet), str(proposal), str(root / "out"),
                        executor_id=executor, executor_model="luna", reviewer_id=reviewer, reviewer_model="luna")
