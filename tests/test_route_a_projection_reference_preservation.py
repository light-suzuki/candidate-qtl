"""Independent Route A references survive projection deduplication into Wet."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.family_candidate_bridge import family_candidate_bridge
from candidate_qtl.relevance_tiers import T2
from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from tests import test_issue568_sequence_keyed_projection as sequence_fixture


STRONG_MECHANISM = "Controls meristem development and branching"
COUNTER_HYPOTHESIS = "The branching effect may reflect an indirect growth response"


class RouteAProjectionReferencePreservationTests(unittest.TestCase):
    def _run(self, root, *, include_sequence=True, reverse=False, same_hypothesis=False):
        files = sequence_fixture.SequenceKeyedProjectionTests()._inputs(
            root, long_form_rows=[("Ps2", "62.5", "0.88")], ambiguity="unique", paralogs="0",
        )
        hypotheses = pd.read_csv(files["hypotheses"], sep="\t", dtype=str).fillna("")
        statuses = pd.read_csv(files["status"], sep="\t", dtype=str).fillna("")
        baits = pd.read_csv(files["bait"], sep="\t", dtype=str).fillna("")
        if not same_hypothesis:
            hypotheses = pd.concat([hypotheses, pd.DataFrame([{
                "hypothesis_id": "H_NOOT", "gene_family": "NOOT", "reference_gene_name": "NOOT1",
                "proposed_stable_id": "P12345", "member_interval_genes": "",
                "mechanism_statement": STRONG_MECHANISM, "expected_role": "meristem development",
                "source_ids": "PMID:456", "counter_hypothesis": COUNTER_HYPOTHESIS,
            }])], ignore_index=True).fillna("")
            statuses = pd.concat([statuses, pd.DataFrame([{
                "hypothesis_id": "H_NOOT", "family_status": "not_evaluated_missing_reference_genome_family_index",
                "interval_gene_ids": "", "source_ids": "PMID:456",
            }])], ignore_index=True)
            baits = pd.concat([baits, pd.DataFrame([{
                "hypothesis_id": "H_NOOT", "function_source_ids": "PMID:456",
                "source_provider": "source_verified", "bait_id": "P12345",
            }])], ignore_index=True)
        for frame, key in ((hypotheses, "hypotheses"), (statuses, "status"), (baits, "bait")):
            (frame.iloc[::-1] if reverse else frame).to_csv(files[key], sep="\t", index=False)
        index = root / "family_index.tsv"
        index.write_text(sequence_fixture.INDEX.replace("NOOT", "BOP") if same_hypothesis else sequence_fixture.INDEX,
                         encoding="utf-8")
        bridge = family_candidate_bridge(
            str(root / "bridge"), route_a_hypotheses_arg=str(files["hypotheses"]),
            route_a_family_status_arg=str(files["status"]), route_a_bait_panel_arg=str(files["bait"]),
            route_b_family_index_arg=str(index), genes_arg=str(files["genes"]), case_id="projection-references",
            route_a_blast_arg=str(files["blast"]) if include_sequence else None,
            route_a_blast_eligible_arg=str(files["eligible"]) if include_sequence else None,
        )
        wet = wet_candidate_handoff(
            str(root / "wet"), route_a_bridge_arg=bridge["bridge_tsv"], route_a_manifest_arg=bridge["manifest"],
            interval_atlas_arg=str(files["genes"]), wet_budget=20, case_id="projection-references",
        )
        return (
            pd.read_csv(bridge["bridge_tsv"], sep="\t", dtype=str).fillna(""),
            json.loads(Path(bridge["manifest"]).read_text(encoding="utf-8")),
            pd.read_csv(wet["wet_candidates"], sep="\t", dtype=str).fillna(""),
        )

    def test_an_independent_sequence_match_preserves_the_stronger_original_reference_and_wet_tier(self):
        original_strong_evidence = None
        for include_sequence in (False, True):
            for reverse in (False, True):
                with self.subTest(include_sequence=include_sequence, reverse=reverse), tempfile.TemporaryDirectory() as tmp:
                    bridge, manifest, wet = self._run(Path(tmp), include_sequence=include_sequence, reverse=reverse)
                    expected_references = {("Ps2", "P12345", "PMID:456")}
                    expected_hypotheses = {"H_NOOT"}
                    if include_sequence:
                        expected_references.add(("Ps2", "Q9M1I7", "PMID:123"))
                        expected_hypotheses.add("H_BOP")
                    self.assertEqual(expected_references, set(bridge[["gene_id", "reference_accession", "source_id"]].itertuples(index=False, name=None)))
                    self.assertEqual([], manifest["route_a_interval_family_projection"]["deduplicated_into_sequence_projection"])
                    self.assertEqual(["Ps2"], list(wet.candidate_gene_id))
                    self.assertEqual(T2, wet.iloc[0].relevance_tier)
                    evidence = json.loads(wet.iloc[0].relevance_evidence_json)
                    originals = [row for row in evidence if row.get("evidence_origin") == "route_a_hypothesis"]
                    self.assertEqual(expected_hypotheses, {row["hypothesis_id"] for row in originals})
                    self.assertEqual(len(expected_hypotheses), len(originals))
                    strong = next(row for row in originals if row["hypothesis_id"] == "H_NOOT")
                    self.assertEqual("P12345", strong["proposed_stable_id"])
                    self.assertEqual("PMID:456", strong["source_ids"])
                    self.assertEqual(STRONG_MECHANISM, strong["mechanism_statement"])
                    self.assertEqual(COUNTER_HYPOTHESIS, strong["counter_hypothesis"])
                    if original_strong_evidence is None:
                        original_strong_evidence = strong
                    self.assertEqual(original_strong_evidence, strong)

    def test_name_and_sequence_for_the_same_hypothesis_remain_one_projection_and_one_original_claim(self):
        with tempfile.TemporaryDirectory() as tmp:
            bridge, manifest, wet = self._run(Path(tmp), same_hypothesis=True)
            self.assertEqual(1, len(bridge))
            self.assertEqual("Q9M1I7", bridge.iloc[0].reference_accession)
            self.assertEqual(["Ps2"], manifest["route_a_interval_family_projection"]["deduplicated_into_sequence_projection"])
            self.assertEqual(["Ps2"], list(wet.candidate_gene_id))
            originals = [row for row in json.loads(wet.iloc[0].relevance_evidence_json)
                         if row.get("evidence_origin") == "route_a_hypothesis"]
            self.assertEqual(["H_BOP"], [row["hypothesis_id"] for row in originals])


if __name__ == "__main__":
    unittest.main()
