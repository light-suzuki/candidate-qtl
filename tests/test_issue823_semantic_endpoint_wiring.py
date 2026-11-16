"""Issue #823 family-index endpoints survive semantic proposal and review wiring."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.route_b_prepare import route_b_prepare
from candidate_qtl.route_b_semantic_bundle import build_route_b_semantic_bundle
from candidate_qtl.route_b_semantic_candidates import freeze_route_b_semantic_candidates
from candidate_qtl.route_b_semantic_import import import_route_b_semantic_proposals
from candidate_qtl.route_b_semantic_source_verify import verify_route_b_semantic_sources
from candidate_qtl.route_b_semantic_transfer_verify import verify_route_b_semantic_transfer


def _write_table(path: Path, rows: list[dict[str, str]]) -> Path:
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, lineterminator="\n")
    return path


class Issue823SemanticEndpointWiringTests(unittest.TestCase):
    def test_family_only_exact_card_reaches_proposal_and_unverified_copy_review(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = _write_table(root / "candidates.tsv", [
                {"gene_id": "copyA", "chrom": "chr1", "start": "100", "end": "200"},
                {"gene_id": "copyB", "chrom": "chr1", "start": "300", "end": "400"},
            ])
            family = _write_table(root / "family_index.tsv", [
                {"gene_id": "copyA", "best_subject": "sp|P11111|KINASE_A", "best_subject_symbol": "KINASE_A",
                 "best_subject_stem": "KINASE", "family_subjects": "sp|P11111|KINASE_A;tr|P22222|KINASE_B",
                 "assignment_status": "assigned", "best_identity": "88", "best_query_coverage": "0.91",
                 "best_evalue": "1e-40", "best_bitscore": "240", "second_bitscore": "220",
                 "top2_bitscore_ratio": "0.92", "symbol_source_rank": "2",
                 "symbol_source_bitscore_ratio": "0.97"},
                {"gene_id": "copyB", "best_subject": "tr|P22222|KINASE_B", "best_subject_symbol": "KINASE_B",
                 "best_subject_stem": "KINASE", "family_subjects": "tr|P22222|KINASE_B",
                 "assignment_status": "assigned", "best_identity": "83", "best_query_coverage": "0.79",
                 "best_evalue": "1e-25", "best_bitscore": "190", "second_bitscore": "0",
                 "top2_bitscore_ratio": "0", "symbol_source_rank": "1",
                 "symbol_source_bitscore_ratio": "1"},
            ])
            corpus = _write_table(root / "reference.tsv", [
                {"reference_gene_id": "P11111", "reference_species": "Species X",
                 "function_description": "developmental kinase involved in meristem signaling",
                 "source_id": "DB:P11111", "source_url": "", "mutant_phenotype": "", "pathway": ""},
                {"reference_gene_id": "P22222", "reference_species": "Species Y",
                 "function_description": "related kinase with a distinct source card",
                 "source_id": "DB:P22222", "source_url": "", "mutant_phenotype": "", "pathway": ""},
                {"reference_gene_id": "p11111", "reference_species": "wrong-case decoy",
                 "function_description": "case mismatch must not join", "source_id": "DB:CASE",
                 "source_url": "", "mutant_phenotype": "", "pathway": ""},
            ])

            prepared = route_b_prepare(
                str(candidates), str(root / "prepared"), family_index_arg=str(family),
                reference_function_args=[str(corpus)], auto_local_db=False,
            )
            bundle_result = build_route_b_semantic_bundle(
                prepared["interval_gene_atlas"], str(root / "bundle"), phenotype="altered meristem growth",
                b_review_membership_policy="source_function_v1",
            )
            bundle = json.loads(Path(bundle_result["bundle"]).read_text(encoding="utf-8"))
            unit = next(unit for unit in bundle["units"] if set(unit["member_interval_genes"]) == {"copyA", "copyB"})
            by_gene = {gene["gene_id"]: gene for gene in unit["genes"]}
            for gene_id, expected in (("copyA", ("240", "220", "0.92")),
                                      ("copyB", ("190", "0", "0"))):
                gene = by_gene[gene_id]
                self.assertEqual("assigned", gene["family_index_assignment_status"])
                self.assertEqual(expected[0], gene["family_index_best_bitscore"])
                self.assertEqual(expected[1], gene["family_index_second_bitscore"])
                self.assertEqual(expected[2], gene["family_index_top2_bitscore_ratio"])
            self.assertEqual("sp|P11111|KINASE_A", by_gene["copyA"]["family_index_best_subject"])
            self.assertEqual("P11111;P22222", by_gene["copyA"]["family_index_subject_endpoints"])
            cards = unit["reference_material_cards"]
            self.assertEqual({"P11111", "P22222"}, {card["reference_gene"] for card in cards})
            self.assertFalse(any(card["source_id"] == "DB:CASE" for card in cards))
            self.assertTrue(all(card["material_status"] == "investigation_context_not_verified_support"
                                and card["rank_contribution"] == "none" for card in cards))
            prompt = json.loads(Path(bundle_result["prompts"]).read_text(encoding="utf-8").splitlines()[0])
            for fragment in ("sp|P11111|KINASE_A", "family_index_second_bitscore", "0.92",
                             "developmental kinase involved in meristem signaling", "DB:P11111"):
                self.assertIn(fragment, prompt["user"])

            response = root / "responses.json"
            response.write_text(json.dumps({"candidates": [{
                "prompt_id": prompt["prompt_id"], "parse_status": "json",
                "model_json": {"claims": [{
                    "claim_type": "mechanistic_review_seed",
                    "biological_function": "developmental kinase signaling",
                    "mechanism": "altered signaling could affect meristem growth",
                    "phenotype_relation": "conditional focal-phenotype hypothesis; transfer is unverified",
                    "reference_gene": "P11111", "reference_species": "Species X",
                    "evidence_scope": "family", "source_ids": [],
                    "alternative_hypothesis": "the related copy or another pathway may explain the phenotype",
                    "uncertainty_note": "family-index similarity is a retrieval handle, not orthology or support",
                }]},
            }]}), encoding="utf-8")
            imported = import_route_b_semantic_proposals(
                bundle_result["bundle"], bundle_result["prompts"], str(response), str(root / "import"),
            )
            proposal = pd.read_csv(imported["proposals"], sep="\t", dtype=str).fillna("")
            self.assertEqual(1, len(proposal))
            self.assertEqual("proposed", proposal.iloc[0]["verification_status"])
            self.assertIn("reference_gene_grounded:P11111", proposal.iloc[0]["verification_reason"])
            self.assertEqual("", proposal.iloc[0]["source_id"])
            self.assertEqual("none", proposal.iloc[0]["rank_contribution"])

            source_result = verify_route_b_semantic_sources(
                imported["proposals"], str(root / "source_verify"), bundle_arg=bundle_result["bundle"],
                allow_external_api=False,
            )
            transfer = verify_route_b_semantic_transfer(
                imported["proposals"], source_result["source_verification"], bundle_result["bundle"],
                str(root / "transfer"), case_id=str(bundle.get("case_id", "")),
            )
            final = pd.read_csv(transfer["claims"], sep="\t", dtype=str).fillna("")
            self.assertEqual("not_evaluated", final.iloc[0]["verification_status"])
            candidates_result = freeze_route_b_semantic_candidates(
                transfer["claims"], str(root / "candidate_review"), case_id=str(bundle.get("case_id", "")),
            )
            self.assertEqual(0, candidates_result["candidate_gene_count"])
            review = pd.read_csv(candidates_result["review_candidates"], sep="\t", dtype=str).fillna("")
            self.assertTrue(review.empty)
            ledger = pd.read_csv(root / "candidate_review/route_b_review_eligibility.tsv", sep="\t", dtype=str)
            self.assertEqual({"copyA", "copyB"}, set(ledger.candidate_gene_id))
            self.assertTrue(ledger.b_review_eligibility.eq("source_absent").all())
            self.assertTrue(ledger.biological_exclusion.eq("not_asserted").all())


if __name__ == "__main__":
    unittest.main()
