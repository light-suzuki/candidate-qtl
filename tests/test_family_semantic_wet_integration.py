"""Integration boundary: semantic copy fallback may change ordering, never Wet membership (#505)."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.canonical_family_semantic_wet import _family_semantic_fallback
from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from tests._issue664_helpers import manifest_for


def _write(path: Path, rows: list[dict[str, str]]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, lineterminator="\n")
    return str(path)


class FamilySemanticWetIntegrationTests(unittest.TestCase):
    def test_verified_semantic_copy_changes_ordering_but_not_membership(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            atlas = _write(root / "b01_route_b" / "interval_gene_atlas.tsv", [
                {"gene_id": "A1", "orthology_reference_gene": "AtAMB", "orthology_relationship": "1:1"},
                {"gene_id": "A2", "orthology_reference_gene": "AtAMB", "orthology_relationship": "one_to_many"},
            ])
            family_copy = _write(root / "canonical" / "family_copy_discrimination.tsv", [
                {"family_id": "AMB", "candidate_gene_id": "A1", "disambiguation_status": "paralog_ambiguous", "orthology_class": "one_to_one"},
                {"family_id": "AMB", "candidate_gene_id": "A2", "disambiguation_status": "paralog_ambiguous", "orthology_class": "one_to_many"},
            ])
            bridge = _write(root / "canonical" / "family_bridge.tsv", [
                {"route": "a", "gene_id": "A1", "family_id": "AMB", "hypothesis_id": "H1", "mapping_tier": "tier3_sequence_match_unique", "ambiguity_status": "paralog_ambiguous", "review_status": "supported"},
                {"route": "a", "gene_id": "A2", "family_id": "AMB", "hypothesis_id": "H1", "mapping_tier": "tier3_sequence_match_unique", "ambiguity_status": "paralog_ambiguous", "review_status": "supported"},
            ])
            bridge_manifest = manifest_for(Path(bridge), route="a", case_id="case1")
            ranked = _write(root / "ranked.tsv", [
                {"gene_id": "A2", "dry_rank": "1", "dry_group": "top", "dry_score": "5"},
                {"gene_id": "A1", "dry_rank": "2", "dry_group": "top", "dry_score": "4"},
            ])

            effective_copy, effective_shortlist, state = _family_semantic_fallback(
                root,
                atlas_arg=atlas,
                bridge_arg=bridge,
                family_copy_arg=family_copy,
                family_copy_shortlist_arg=None,
                case_id="case1",
            )
            self.assertEqual("semantic_review_available_no_response", state["state"])
            self.assertEqual(family_copy, effective_copy)
            self.assertIsNone(effective_shortlist)

            prompts = [
                json.loads(line)
                for line in Path(state["prompts"]).read_text(encoding="utf-8").splitlines()
                if line.strip()
            ]
            self.assertEqual(1, len(prompts))
            run_dir = root / "b04_canonical_semantic" / "10_family_semantic" / "02_llm_run"
            run_dir.mkdir(parents=True, exist_ok=True)
            (run_dir / "llm_candidate_responses.json").write_text(json.dumps({
                "candidates": [{
                    "prompt_id": prompts[0]["prompt_id"],
                    "gene_id": "__family_semantic_resolution__",
                    "parse_status": "json",
                    "model_json": {
                        "resolutions": [{
                            "family_id": "AMB",
                            "semantic_relation": "likely_reference_relative",
                            "candidate_gene_ids_considered": ["A1", "A2"],
                            "proposed_reference_mapping": [
                                {"candidate_gene_id": "A1", "reference_gene": "AtAMB"}
                            ],
                            "reasoning_summary": "A1 corresponds to the supplied reference biology",
                            "evidence_to_verify": ["reference orthology"],
                            "competing_explanations_to_keep": ["A2 remains a paralog candidate"],
                            "review_required": True,
                        }]
                    },
                }]
            }, ensure_ascii=False), encoding="utf-8")

            effective_copy, effective_shortlist, state = _family_semantic_fallback(
                root,
                atlas_arg=atlas,
                bridge_arg=bridge,
                family_copy_arg=family_copy,
                family_copy_shortlist_arg=None,
                case_id="case1",
            )
            self.assertEqual("verified_semantic_copy_resolution", state["state"])
            self.assertEqual(1, state["promoted_family_count"])
            copy_frame = pd.read_csv(effective_copy, sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual("preferred_copy", copy_frame.loc["A1", "disambiguation_status"])
            self.assertEqual("paralog_ambiguous", copy_frame.loc["A2", "disambiguation_status"])

            handoff = wet_candidate_handoff(
                str(root / "wet"),
                route_a_bridge_arg=bridge,
                route_a_manifest_arg=str(bridge_manifest),
                ranked_arg=ranked,
                family_copy_arg=effective_copy,
                family_copy_shortlist_arg=effective_shortlist,
                wet_budget=20,
                case_id="case1",
            )
            wet = pd.read_csv(handoff["wet_candidates"], sep="\t", dtype=str).fillna("")
            narrowed = pd.read_csv(handoff["stage_a_narrowed"], sep="\t", dtype=str).fillna("")

            self.assertEqual({"A1", "A2"}, set(wet["candidate_gene_id"]))
            self.assertEqual({"A1", "A2"}, set(narrowed["gene_id"]))
            self.assertEqual("A1", narrowed.iloc[0]["gene_id"])
            self.assertEqual("N1", narrowed.iloc[0]["narrowing_tier"])
            self.assertEqual("none", narrowed.iloc[0]["rank_contribution"])


if __name__ == "__main__":
    unittest.main()
