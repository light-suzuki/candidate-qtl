"""Partial homology is retained as review, without changing strict support."""
import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from candidate_qtl.blast_evidence import run_panel_blastp
from candidate_qtl.family_candidate_bridge import family_candidate_bridge
from candidate_qtl.partial_alignment_review import partial_review_hits
from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from tests.test_blast_evidence import FIXTURE, MockBlastTools, _freeze, S


class PartialAlignmentReviewTests(unittest.TestCase):
    def test_thousand_partial_copies_are_not_cut_to_wet_budget(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            frozen = _freeze(root)
            genes = [f"G{i:04d}" for i in range(1000)]
            proteins = root / "proteins.fa"
            proteins.write_text("".join(f">{gene}\n{S}\n" for gene in genes))
            raw = "".join(f"bait1\t{gene}\t80\t40\t120\t120\t1\t40\t1\t40\t1e-20\t90\n" for gene in genes)
            with MockBlastTools(raw_hits=raw):
                search = run_panel_blastp(frozen["manifest"], str(proteins), str(root / "search"))
            import hashlib
            panel = root / "projection.tsv"
            pd.DataFrame([{"hypothesis_id": "h1", "bait_id": "bait1", "sequence_sha256": hashlib.sha256(S.encode()).hexdigest(),
                           "review_status": "verified", "projection_eligibility": "verified_support"}]).to_csv(panel, sep="\t", index=False)
            atlas = root / "genes.tsv"
            pd.DataFrame({"gene_id": genes}).to_csv(atlas, sep="\t", index=False)
            result = family_candidate_bridge(str(root / "bridge"), route_b_atlas_arg=str(atlas), genes_arg=str(atlas),
                                            route_a_projection_panel_arg=str(panel), route_a_blast_arg=search["evidence"], case_id="unit")
            wet_result = wet_candidate_handoff(str(root / "wet"), route_a_bridge_arg=result["bridge_tsv"], route_a_manifest_arg=result["manifest"],
                                              review_candidates_arg=result["review_candidates"], review_candidates_manifest_arg=result["review_manifest"],
                                              interval_atlas_arg=str(atlas), case_id="unit")
            wet = pd.read_csv(wet_result["wet_candidates"], sep="\t")
            self.assertEqual(set(wet.candidate_gene_id), set(genes))
            self.assertEqual(set(wet.support_tier), {"unverified_route_review"})
            self.assertTrue(wet_result["adaptive_synthesis_request"])

    def test_partial_copies_reach_wet_without_verified_support(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            frozen = _freeze(root)
            raw = "".join(f"bait1\t{gene}\t{identity}\t{span}\t120\t120\t1\t{span}\t1\t{span}\t{evalue}\t90\n"
                          for gene, identity, span, evalue in [
                              ("g1", 90, 40, "1e-30"), ("g2", 80, 40, "1e-20"),
                              ("g4", 90, 120, "1e-40"),  # strict hit, not a partial-review row
                              ("g6", 20, 40, "1e-20"), ("g7", 90, 40, "0.1"),
                              ("g8", 90, 20, "1e-20"),
                          ])
            with MockBlastTools(raw_hits=raw):
                search = run_panel_blastp(frozen["manifest"], str(FIXTURE / "interval_proteins.fa"), str(root / "search"))
            panel = root / "projection.tsv"
            import hashlib
            records = [{"hypothesis_id": "h1", "bait_id": "bait1", "sequence_sha256": hashlib.sha256(S.encode()).hexdigest(),
                        "review_status": "verified", "projection_eligibility": "verified_support"}]
            pd.DataFrame(records).to_csv(panel, sep="\t", index=False)
            atlas = root / "genes.tsv"
            pd.DataFrame({"gene_id": ["g1", "g2", "g4", "g6", "g7", "g8"]}).to_csv(atlas, sep="\t", index=False)
            result = family_candidate_bridge(str(root / "bridge"), route_b_atlas_arg=str(atlas), genes_arg=str(atlas),
                                            route_a_projection_panel_arg=str(panel), route_a_blast_arg=search["evidence"], case_id="unit")
            review = pd.read_csv(result["review_candidates"], sep="\t")
            self.assertEqual(set(review.candidate_gene_id), {"g1", "g2"})
            self.assertEqual(set(review.reason), {"partial_alignment_requires_domain_or_copy_review"})
            wet_result = wet_candidate_handoff(str(root / "wet"), route_a_bridge_arg=result["bridge_tsv"], route_a_manifest_arg=result["manifest"], review_candidates_arg=result["review_candidates"],
                                              review_candidates_manifest_arg=result["review_manifest"], interval_atlas_arg=str(atlas), case_id="unit")
            wet = pd.read_csv(wet_result["wet_candidates"], sep="\t")
            self.assertEqual(set(wet.candidate_gene_id), {"g1", "g2"})
            self.assertEqual(set(wet.loc[wet.candidate_gene_id.isin(["g1", "g2"]), "support_tier"]), {"unverified_route_review"})
            wet_manifest = json.loads(Path(wet_result["wet_candidates_manifest"]).read_text())
            trace = pd.read_csv(wet_manifest["candidate_trace"], sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            for gene in ("g1", "g2"):
                self.assertEqual("True", trace.loc[gene, "explicit_review_member"])
                self.assertEqual("True", trace.loc[gene, "adopted_first_stage"])
            for gene in ("g4", "g6", "g7", "g8"):
                self.assertEqual("False", trace.loc[gene, "explicit_review_member"])
                self.assertEqual("False", trace.loc[gene, "adopted_first_stage"])
                self.assertEqual("not_selected_no_upstream_candidate", trace.loc[gene, "selection_disposition"])
            # A strict sequence hit is not, by itself, Route-A support or a Wet seed.
            self.assertEqual("False", trace.loc["g4", "supported_union_member"])
            strict = pd.read_csv(search["eligible_hits"], sep="\t")
            self.assertEqual(len(strict), 1)
            # Both incomplete and complete source claims retain partial matches.
            records[0].update(review_status="review_required", projection_eligibility="review_only")
            hits = partial_review_hits(search["evidence"], records, root)
            self.assertEqual({hit["target_id"] for hit in hits["bait1"]}, {"g1", "g2"})
            records[0]["sequence_sha256"] = "0" * 64
            with self.assertRaisesRegex(ValueError, "sequence differs"):
                partial_review_hits(search["evidence"], records, root)
            manifest = json.loads(Path(search["manifest"]).read_text())
            Path(manifest["outputs"]["raw_hits"]).write_text(raw + "\n")
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                partial_review_hits(search["evidence"], records, root)

    def test_legacy_or_missing_search_cannot_invent_partial_hits(self):
        with tempfile.TemporaryDirectory() as directory:
            self.assertEqual(partial_review_hits(str(Path(directory) / "evidence.tsv"), [], Path(directory)), {})
