"""Empty resolved hypotheses must remain consumable without creating support."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl.llm_pathway import bait_panel_freeze, hypothesis_verify


class EmptyHypothesisTests(unittest.TestCase):
    def test_empty_hypotheses_verify_and_freeze(self):
        for annotated in (False, True):
            with self.subTest(annotated=annotated), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                hypotheses = root / "hypotheses.tsv"
                columns = ["hypothesis_id", "case_id", "gene_family", "reference_gene_name",
                           "reference_species", "proposed_stable_id"]
                if annotated:
                    columns.append("verified_source_ids")
                pd.DataFrame(columns=columns).to_csv(hypotheses, sep="\t", index=False)
                with patch("urllib.request.urlopen", side_effect=AssertionError("unexpected network call")):
                    result = hypothesis_verify(str(hypotheses), str(root / "verify"),
                                               focal_species="Pisum sativum", allow_external_api=False)
                    panel = bait_panel_freeze(
                        str(hypotheses), result["validation"], result["verified_sequences"],
                        str(root / "panel"), case_id="empty",
                        reviewable_sequences_arg=result["reviewable_sequences"])
                manifest = json.loads(Path(result["manifest"]).read_text())
                expected_scope = ("enforced_from_source_verify_annotations" if annotated
                                  else "legacy_standalone_input_without_source_annotations")
                self.assertEqual(manifest["publication_gate_scope"], expected_scope)
                for key in ("proposal_count", "verified_count", "rejected_count", "review_required_count"):
                    self.assertEqual(manifest[key], 0)
                for key in ("validation", "rejected"):
                    self.assertTrue(pd.read_csv(result[key], sep="\t").empty)
                reviews = pd.read_csv(result["review_required"], sep="\t")
                self.assertEqual(reviews["reason"].tolist(), ["panel_below_minimum_size"])
                for key in ("verified_sequences", "reviewable_sequences"):
                    self.assertEqual(json.loads(Path(result[key]).read_text()), {})
                for key in ("panel", "projection_panel"):
                    self.assertTrue(pd.read_csv(panel[key], sep="\t").empty)
