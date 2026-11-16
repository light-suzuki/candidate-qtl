"""Count targets never remove supported copies or block the default handoff."""
from pathlib import Path
import json
import tempfile
import unittest

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from tests._issue664_helpers import manifest_for
from tests.test_review_candidate_retention import _sealed_review


class OperationalWetTargetTests(unittest.TestCase):
    def test_review_member_triggers_second_stage_and_round_limit_preserves_wet(self):
        for rounds in (0, 1):
            with self.subTest(rounds=rounds), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                source = root / "supported.tsv"
                ids = [f"G{i:03d}" for i in range(25)]
                pd.DataFrame({"candidate_gene_id": ids, "route": "b",
                              "support_status": "verified_semantic_support"}).to_csv(source, sep="\t", index=False)
                atlas = root / "atlas.tsv"
                pd.DataFrame({"gene_id": ids + ["G2"], "family_id": "F1"}).to_csv(atlas, sep="\t", index=False)
                review = root / "review.tsv"
                seal = _sealed_review(review, "mixed", atlas)
                result = wet_candidate_handoff(
                    str(root / "out"), route_b_candidates_arg=str(source),
                    route_b_candidates_manifest_arg=str(manifest_for(source, route="b", case_id="mixed")),
                    interval_atlas_arg=str(atlas), review_candidates_arg=str(review),
                    review_candidates_manifest_arg=str(seal), case_id="mixed", max_adaptive_rounds=rounds,
                )
                self.assertEqual(25, result["candidate_count"])
                self.assertEqual("requested" if rounds else "review_required_round_budget_exhausted",
                                 result["second_stage_state"])
                wet = pd.read_csv(result["wet_candidates"], sep="\t")
                self.assertEqual(set(ids + ["G2"]), set(wet.candidate_gene_id))
                if rounds:
                    request = json.loads(Path(result["adaptive_synthesis_request"]).read_text())
                    self.assertEqual(set(ids + ["G2"]), set(request["candidate_gene_ids"]))
                else:
                    self.assertIsNone(result["adaptive_synthesis_request"])

    def test_default_cli_freezes_all_supported_candidates_above_target(self):
        for count in (19, 20, 21, 25, 26, 27):
            with self.subTest(count=count), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                source = root / "route_b_verified_candidates.tsv"
                ids = [f"G{i:03d}" for i in range(count)]
                pd.DataFrame({"candidate_gene_id": ids, "route": "b",
                              "support_status": "verified_semantic_support"}).to_csv(source, sep="\t", index=False)
                manifest = manifest_for(source, route="b", case_id="soft")
                out = root / "out"
                self.assertEqual(0, main([
                    "wet-candidate-handoff", "--route-b-candidates", str(source),
                    "--route-b-candidates-manifest", str(manifest), "--wet-budget", "20",
                    "--case-id", "soft", "--out-dir", str(out), "--json",
                ]))
                result = json.loads((out / "wet_candidate_handoff_manifest.json").read_text())
                self.assertEqual("wet_reviewable", result["handoff_state"])
                self.assertEqual(count > 25, result["heavy_stage_invoked"])
                self.assertEqual(
                    "requested" if count > 25 else "skipped_at_or_below_25",
                    result["second_stage_state"],
                )
                if count > 25:
                    request = json.loads(Path(result["adaptive_synthesis_request"]).read_text())
                    self.assertEqual(set(ids), set(request["candidate_gene_ids"]))
                else:
                    self.assertIsNone(result["adaptive_synthesis_request"])
                self.assertEqual("target", result["budget_policy"])
                self.assertEqual(max(0, count - 20), result["above_target_count"])
                wet = pd.read_csv(out / "wet_candidates.tsv", sep="\t")
                self.assertEqual(set(ids), set(wet.candidate_gene_id))
                self.assertTrue(wet.support_tier.eq("tier1_verified_support").all())

    def test_explicit_strict_capacity_retains_legacy_refinement_request(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "b.tsv"
            pd.DataFrame({"candidate_gene_id": [f"G{i}" for i in range(21)], "route": "b",
                          "support_status": "verified_semantic_support"}).to_csv(source, sep="\t", index=False)
            manifest = manifest_for(source, route="b", case_id="strict")
            result = wet_candidate_handoff(str(root / "out"), route_b_candidates_arg=str(source),
                route_b_candidates_manifest_arg=str(manifest), case_id="strict", budget_policy="strict")
            self.assertEqual("adaptive_refinement_requested", result["handoff_state"])
            self.assertIsNone(result["wet_candidates"])
            self.assertEqual(21, result["candidate_count"])
