from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.family_candidate_bridge import family_candidate_bridge
from candidate_qtl.dual_route import phenotype_first_status
from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from tests.test_blast_evidence import FIXTURE, MockBlastTools, S, _freeze


class RouteAReviewRelevanceTests(unittest.TestCase):
    def _run(self, root: Path, hypothesis_fields=None, *, alignment_span=40, review_only=False):
        frozen = _freeze(root)
        raw = "".join(
            f"bait1\t{gene}\t{identity}\t{alignment_span}\t120\t120\t1\t{alignment_span}\t1\t{alignment_span}\t{evalue}\t90\n"
            for gene, identity, evalue in (("g1", 90, "1e-30"), ("g2", 80, "1e-20"))
        )
        with MockBlastTools(raw_hits=raw):
            search = __import__("candidate_qtl.blast_evidence", fromlist=["run_panel_blastp"]).run_panel_blastp(
                frozen["manifest"], str(FIXTURE / "interval_proteins.fa"), str(root / "search")
            )
        hypotheses = root / "hypotheses.tsv"
        pd.DataFrame([{
            "hypothesis_id": "h1", "gene_family": "trait-family",
            "reference_gene_name": "REF1", "proposed_stable_id": "bait1",
            "member_interval_genes": "", "mechanism_statement": "controls trait development",
            "expected_role": "trait1",
            **(hypothesis_fields or {}),
        }]).to_csv(hypotheses, sep="\t", index=False, lineterminator="\n")
        panel = root / "projection.tsv"
        pd.DataFrame([{
            "hypothesis_id": "h1", "bait_id": "bait1",
            "sequence_sha256": hashlib.sha256(S.encode()).hexdigest(),
            "review_status": "review_required" if review_only else "verified",
            "projection_eligibility": "review_only" if review_only else "verified_support",
        }]).to_csv(panel, sep="\t", index=False, lineterminator="\n")
        route_a = phenotype_first_status(
            str(root / "a-status"), hypotheses_arg=str(hypotheses), panel_arg=str(panel),
            case_id="route-a-review",
        )
        normalized_hypotheses = Path(route_a["hypotheses"])
        atlas = root / "genes.tsv"
        pd.DataFrame({"gene_id": ["g1", "g2"]}).to_csv(atlas, sep="\t", index=False, lineterminator="\n")
        bridge = family_candidate_bridge(
            str(root / "bridge"), route_a_hypotheses_arg=str(normalized_hypotheses),
            route_b_atlas_arg=str(atlas), genes_arg=str(atlas),
            route_a_projection_panel_arg=str(panel), route_a_blast_arg=search["evidence"],
            case_id="route-a-review",
        )
        wet = wet_candidate_handoff(
            str(root / "wet"), route_a_bridge_arg=bridge["bridge_tsv"],
            route_a_manifest_arg=bridge["manifest"],
            review_candidates_arg=bridge["review_candidates"],
            review_candidates_manifest_arg=bridge["review_manifest"],
            interval_atlas_arg=str(atlas), wet_budget=2, case_id="route-a-review",
        )
        return wet, bridge, normalized_hypotheses, root

    def test_partial_review_keeps_original_hypothesis_at_wet_consumer(self):
        with tempfile.TemporaryDirectory() as tmp:
            wet, _, normalized_hypotheses, root = self._run(Path(tmp))
            normalized = pd.read_csv(normalized_hypotheses, sep="\t", dtype=str).fillna("")
            self.assertEqual("controls trait development", normalized.loc[0, "mechanism_statement"])
            self.assertEqual("trait1", normalized.loc[0, "expected_role"])
            self.assertEqual("bait1", normalized.loc[0, "proposed_stable_id"])
            frame = pd.read_csv(root / "wet" / "wet_candidates.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual({"g1", "g2"}, set(frame["candidate_gene_id"]))
            self.assertEqual({"T2_trait_or_mechanism_linked"}, set(frame["relevance_tier"]))
            for value in frame["relevance_evidence_json"]:
                evidence = json.loads(value)
                self.assertEqual({"h1"}, {row["hypothesis_id"] for row in evidence})
                self.assertTrue(any(row.get("mechanism_statement") for row in evidence))
            self.assertEqual("wet_reviewable_degraded_coverage", wet["handoff_state"])

    def test_tampered_hypothesis_is_rejected_by_review_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, bridge, hypotheses, root = self._run(Path(tmp))
            hypotheses.write_text(hypotheses.read_text(encoding="utf-8") + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                wet_candidate_handoff(
                    str(root / "wet-tampered"), route_a_bridge_arg=bridge["bridge_tsv"],
                    route_a_manifest_arg=bridge["manifest"],
                    review_candidates_arg=bridge["review_candidates"],
                    review_candidates_manifest_arg=bridge["review_manifest"],
                    interval_atlas_arg=str(root / "genes.tsv"), wet_budget=2,
                    case_id="route-a-review",
                )

    def test_paper_reports_reach_a_projection_and_wet_gene_aggregation(self):
        from candidate_qtl.paper_reports import PAPER_FIELD, normalize_paper_reports
        from candidate_qtl.relevance_tiers import NEW_TIERS
        reports = [{"source_id": "PMID:123", "initial_tier": "T2", "initial_tier_reason": "Similar observed meristem phenotype",
                    "reported_gene": "REF1", "reported_species": "Plant alpha", "reported_phenotype": "enlarged meristem",
                    "reported_process": "meristem proliferation", "phenotype_comparison": "same_cellular_process"},
                   {"source_id": "PMID:456", "initial_tier": "T3", "initial_tier_reason": "Related root process",
                    "reported_gene": "REF1", "reported_species": "Plant alpha", "reported_phenotype": "altered root",
                    "reported_process": "meristem proliferation", "phenotype_comparison": "shared_process"}]
        with tempfile.TemporaryDirectory() as tmp:
            wet, _, _, _ = self._run(Path(tmp), {"source_ids": "PMID:123;PMID:456", PAPER_FIELD: normalize_paper_reports(reports)})
            rows = pd.read_csv(wet["wet_candidates"], sep="\t").fillna("")
            self.assertEqual({"g1", "g2"}, set(rows.candidate_gene_id))
            self.assertEqual({NEW_TIERS[1]}, set(rows.relevance_tier))
            self.assertEqual({"PMID:123"}, set(rows.tier_source_id))
            self.assertEqual({"strongest_paper_v1"}, set(rows.tier_aggregation))


if __name__ == "__main__":
    unittest.main()
