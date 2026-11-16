from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.claim_binding import BOUND_CLAIM_COLUMNS
from candidate_qtl.route_b_semantic_transfer_verify import verify_route_b_semantic_transfer
from candidate_qtl.route_b_semantic_candidates import freeze_route_b_semantic_candidates
from candidate_qtl.semantic_claims import CLAIM_COLUMNS, narrowing_eligible_claims


def _claim(**overrides: str) -> dict[str, str]:
    row = {column: "" for column in CLAIM_COLUMNS}
    row.update({
        "claim_id": "claim-1",
        "route": "interval_first",
        "semantic_unit_id": "unit-1",
        "family_id": "EV1",
        "orthogroup_id": "OG1",
        "subfamily_or_clade": "SUB1",
        "member_interval_genes": "G1;G2",
        "claim_type": "mutant_phenotype",
        "biological_function": "meristem maintenance",
        "mechanism": "receptor signaling",
        "phenotype_relation": "fasciated shoot meristem",
        "reference_gene": "REF1",
        "reference_species": "Arabidopsis thaliana",
        "evidence_scope": "family",
        "source_id": "PMID:1234567",
        "source_url": "",
        "alternative_hypothesis": "related receptor copy",
        "uncertainty_note": "copy identity remains unresolved",
        "provider": "fixture",
        "model": "fixture-model",
        "prompt_fingerprint": "p" * 64,
        "input_fingerprint": "i" * 64,
        "output_fingerprint": "o" * 64,
        "verification_status": "proposed",
        "verification_reason": "source relation pending",
        "rank_contribution": "none",
    })
    row.update(overrides)
    return row


def _bundle(path: Path, *, members: list[tuple[str, str]]) -> None:
    path.write_text(json.dumps({
        "schema_version": "1.0",
        "case_id": "case403",
        "route": "interval_first",
        "units": [{
            "unit_id": "unit-1",
            "unit_type": "family",
            "unit_anchor": "EV1",
            "member_interval_genes": [gene for gene, _rel in members],
            "identity": {
                "evolutionary_family_ids": ["EV1"],
                "family_ids": ["F1"],
                "orthogroup_ids": ["OG1"],
                "subfamily_ids": ["SUB1"],
            },
            "genes": [{
                "gene_id": gene,
                "orthology_reference_gene": "REF1",
                "orthology_relationship": relationship,
            } for gene, relationship in members],
            "literature_link_candidates": [{"reference_gene": "REF1", "source_id": "PMID:1234567"}],
            "stable_source_ids": ["PMID:1234567"],
            "rank_contribution": "none",
        }],
    }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def _source(path: Path, claim_id: str, status: str) -> None:
    pd.DataFrame([{
        "claim_id": claim_id,
        "source_relation_status": status,
    }]).to_csv(path, sep="\t", index=False, lineterminator="\n")


def _bound_claim(claim: dict[str, str], *, complete: bool = True) -> dict[str, str]:
    row = {column: "" for column in BOUND_CLAIM_COLUMNS}
    row.update({
        "bound_claim_id": f"bound-{claim['claim_id']}",
        "route": "b",
        "claim_id": claim["claim_id"],
        "citation": claim["source_id"],
        "source_record_id": claim["source_id"],
        "subject_gene_identity": claim["reference_gene"],
        "subject_gene_verification_state": "verified_exact",
        "experiment_perturbation": "loss-of-function allele",
        "experiment_subject_gene_identity": claim["reference_gene"],
        "experiment_verification_state": "verified_same_subject",
        "predicate_relation": "causal_perturbation_supports_target_trait",
        "relation_verification_state": "verified",
        "object_phenotype_process": claim["phenotype_relation"],
        "object_match_state": "target_trait",
        "organism": claim["reference_species"],
        "organism_verification_state": "verified_match",
        "accession": claim["reference_gene"],
        "accession_verification_state": "verified_match",
        "quoted_support": "The reference-gene mutant produced the target phenotype.",
        "locator": "abstract:sentence-3",
        "quote_binding_state": "verified_same_claim",
        "source_scope": "abstract",
        "binding_state": "verified_complete" if complete else "review_required",
        "curation_method": "curated_human",
        "reviewed_by": "fixture-reviewer",
        "binding_reason": "all relation axes verified" if complete else "direct phenotype relation not established",
        "rank_contribution": "none",
    })
    return row


class RouteBSemanticTransferVerifyTests(unittest.TestCase):
    def test_review_only_seed_blocks_group_and_copy_support_but_freezes_review(self) -> None:
        for scope, members in (("family", [("G1", "1:many"), ("G2", "1:many")]),
                               ("interval_gene", [("G1", "1:1")])):
            with self.subTest(scope=scope), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                claim = _claim(claim_type="mechanistic_review_seed", evidence_scope=scope,
                               member_interval_genes=";".join(gene for gene, _ in members))
                result = self._run(root, claim, members, "supported_mechanistic")
                final = pd.read_csv(result["claims"], sep="\t", dtype=str).fillna("")
                transfer = pd.read_csv(result["transfer_verification"], sep="\t", dtype=str).fillna("").iloc[0]
                self.assertEqual("not_evaluated", final.iloc[0]["verification_status"])
                self.assertIn("review_only_claim_type:mechanistic_review_seed", final.iloc[0]["verification_reason"])
                for column in set(CLAIM_COLUMNS) - {"verification_status", "verification_reason"}:
                    self.assertEqual(claim[column], final.iloc[0][column], column)
                self.assertEqual("verified_complete", transfer["bound_claim_state"])
                self.assertEqual("supported_mechanistic", transfer["source_relation_status"])
                self.assertEqual("source_support_insufficient", transfer["transfer_status"])
                self.assertEqual("none", transfer["transfer_scope"])
                self.assertEqual("", transfer["supported_member_genes"])
                self.assertEqual(0, len(narrowing_eligible_claims(final)))
                frozen = freeze_route_b_semantic_candidates(result["claims"], str(root / "frozen"), case_id="case403")
                self.assertEqual(0, frozen["candidate_gene_count"])
                review = pd.read_csv(frozen["review_candidates"], sep="\t", dtype=str).fillna("")
                self.assertEqual({gene for gene, _ in members}, set(review["candidate_gene_id"]))
                self.assertEqual({"unverified_route_b_review"}, set(review["support_status"]))

    def test_review_only_guard_preserves_negative_identity_and_unreviewed_verdicts(self) -> None:
        for source_status, family, expected in (("contradicted", "EV1", "contradicted"),
                                                ("supported_mechanistic", "WRONG_FAMILY", "unsupported"),
                                                ("not_evaluated", "EV1", "not_evaluated")):
            with self.subTest(source_status=source_status, family=family), tempfile.TemporaryDirectory() as tmp:
                result = self._run(Path(tmp), _claim(claim_type="mechanistic_review_seed", family_id=family),
                                   [("G1", "1:many"), ("G2", "1:many")], source_status)
                final = pd.read_csv(result["claims"], sep="\t", dtype=str).fillna("").iloc[0]
                self.assertEqual(expected, final["verification_status"])
                self.assertNotIn("semantic support promotion prohibited", final["verification_reason"])

    def _run(
        self,
        root: Path,
        claim: dict[str, str],
        members: list[tuple[str, str]],
        source_status: str,
        *,
        include_binding: bool = True,
    ):
        proposals = root / "proposals.tsv"
        pd.DataFrame([claim], columns=CLAIM_COLUMNS).to_csv(proposals, sep="\t", index=False, lineterminator="\n")
        bundle = root / "bundle.json"
        _bundle(bundle, members=members)
        source = root / "source.tsv"
        _source(source, claim["claim_id"], source_status)
        bound = root / "bound.tsv"
        pd.DataFrame(
            [_bound_claim(claim, complete=source_status != "supported_curated_function")],
            columns=BOUND_CLAIM_COLUMNS,
        ).to_csv(bound, sep="\t", index=False, lineterminator="\n")
        return verify_route_b_semantic_transfer(
            str(proposals), str(source), str(bundle), str(root / "out"),
            bound_claims_arg=str(bound) if include_binding else None, case_id="case403"
        )

    def test_mechanistic_cue_without_bound_claim_stays_not_evaluated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = self._run(
                root, _claim(), [("G1", "1:many"), ("G2", "1:many")],
                "supported_mechanistic", include_binding=False,
            )
            final = pd.read_csv(result["claims"], sep="\t", dtype=str).fillna("")
            self.assertEqual("not_evaluated", final.iloc[0]["verification_status"])
            self.assertIn("no shared bound claim", final.iloc[0]["verification_reason"])

    def test_family_level_mechanistic_support_verifies_family_but_keeps_all_copies(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = self._run(
                root, _claim(), [("G1", "1:many"), ("G2", "1:many")], "supported_mechanistic"
            )
            final = pd.read_csv(result["claims"], sep="\t", dtype=str).fillna("")
            transfer = pd.read_csv(result["transfer_verification"], sep="\t", dtype=str).fillna("")
            self.assertEqual("verified", final.iloc[0]["verification_status"])
            self.assertEqual("family_supported_copy_unresolved", transfer.iloc[0]["transfer_status"])
            self.assertEqual("G1;G2", transfer.iloc[0]["supported_member_genes"])
            self.assertEqual(1, len(narrowing_eligible_claims(final)))

    def test_curated_function_without_direct_phenotype_binding_stays_not_evaluated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = self._run(
                root, _claim(source_id="UniProt:REF1"),
                [("G1", "1:many"), ("G2", "1:many")],
                "supported_curated_function",
            )
            final = pd.read_csv(result["claims"], sep="\t", dtype=str).fillna("")
            transfer = pd.read_csv(result["transfer_verification"], sep="\t", dtype=str).fillna("")
            self.assertEqual("not_evaluated", final.iloc[0]["verification_status"])
            self.assertIn("bound claim", final.iloc[0]["verification_reason"])
            self.assertEqual("source_support_insufficient", transfer.iloc[0]["transfer_status"])
            self.assertEqual("G1;G2", transfer.iloc[0]["supported_member_genes"])

    def test_curated_function_copy_claim_still_requires_one_to_one_relation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = self._run(
                root,
                _claim(source_id="UniProt:REF1", member_interval_genes="G1", evidence_scope="interval_gene"),
                [("G1", "1:many")],
                "supported_curated_function",
            )
            final = pd.read_csv(result["claims"], sep="\t", dtype=str).fillna("")
            self.assertEqual("not_evaluated", final.iloc[0]["verification_status"])

    def test_copy_level_claim_requires_single_member_and_one_to_one_relation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = self._run(
                root,
                _claim(member_interval_genes="G1", evidence_scope="interval_gene"),
                [("G1", "1:1")],
                "supported_mechanistic",
            )
            final = pd.read_csv(result["claims"], sep="\t", dtype=str).fillna("")
            transfer = pd.read_csv(result["transfer_verification"], sep="\t", dtype=str).fillna("")
            self.assertEqual("verified", final.iloc[0]["verification_status"])
            self.assertEqual("copy_supported", transfer.iloc[0]["transfer_status"])
            self.assertEqual("G1", transfer.iloc[0]["supported_member_genes"])

    def test_copy_level_claim_does_not_collapse_one_to_many_family(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = self._run(
                root,
                _claim(evidence_scope="interval_gene"),
                [("G1", "1:many"), ("G2", "1:many")],
                "supported_mechanistic",
            )
            final = pd.read_csv(result["claims"], sep="\t", dtype=str).fillna("")
            transfer = pd.read_csv(result["transfer_verification"], sep="\t", dtype=str).fillna("")
            self.assertEqual("not_evaluated", final.iloc[0]["verification_status"])
            self.assertEqual("family_supported_copy_unresolved", transfer.iloc[0]["transfer_status"])
            self.assertEqual("G1;G2", transfer.iloc[0]["supported_member_genes"])
            self.assertEqual(0, len(narrowing_eligible_claims(final)))

    def test_associative_source_support_is_not_promoted_to_verified_mechanism(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = self._run(
                root, _claim(), [("G1", "1:1"), ("G2", "1:many")], "supported_associative"
            )
            final = pd.read_csv(result["claims"], sep="\t", dtype=str).fillna("")
            self.assertEqual("not_evaluated", final.iloc[0]["verification_status"])

    def test_identity_tampering_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = self._run(
                root,
                _claim(family_id="OTHER"),
                [("G1", "1:1"), ("G2", "1:many")],
                "supported_mechanistic",
            )
            final = pd.read_csv(result["claims"], sep="\t", dtype=str).fillna("")
            transfer = pd.read_csv(result["transfer_verification"], sep="\t", dtype=str).fillna("")
            self.assertEqual("unsupported", final.iloc[0]["verification_status"])
            # The transfer artifact must still account for the claim. An
            # unsupported claim is a recorded verdict, not a dropped row.
            self.assertEqual(list(final["claim_id"]), list(transfer["claim_id"]))

    def test_four_value_verification_statuses(self) -> None:
        """Verify that verified, unsupported, contradicted, and not_evaluated are all producible."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            claims = [
                _claim(claim_id="c_verified", member_interval_genes="G1;G2", evidence_scope="family"),
                _claim(claim_id="c_unsupported", family_id="FAKE_FAMILY"),
                _claim(claim_id="c_contradicted", member_interval_genes="G1;G2", evidence_scope="family"),
                _claim(claim_id="c_not_evaluated", member_interval_genes="G1;G2", evidence_scope="interval_gene"),
            ]
            proposals = root / "proposals.tsv"
            pd.DataFrame(claims, columns=CLAIM_COLUMNS).to_csv(proposals, sep="\t", index=False, lineterminator="\n")
            bundle = root / "bundle.json"
            _bundle(bundle, members=[("G1", "1:many"), ("G2", "1:many")])
            source = root / "source.tsv"
            pd.DataFrame([
                {"claim_id": "c_verified", "source_relation_status": "supported_mechanistic"},
                {"claim_id": "c_unsupported", "source_relation_status": "supported_mechanistic"},
                {"claim_id": "c_contradicted", "source_relation_status": "contradicted"},
                {"claim_id": "c_not_evaluated", "source_relation_status": "supported_mechanistic"},
            ]).to_csv(source, sep="\t", index=False, lineterminator="\n")
            pd.DataFrame(
                [_bound_claim(claim) for claim in claims],
                columns=BOUND_CLAIM_COLUMNS,
            ).to_csv(root / "bound.tsv", sep="\t", index=False, lineterminator="\n")

            result = verify_route_b_semantic_transfer(
                str(proposals), str(source), str(bundle), str(root / "out"),
                bound_claims_arg=str(root / "bound.tsv"), case_id="case4val"
            )
            final = pd.read_csv(result["claims"], sep="\t", dtype=str).fillna("").set_index("claim_id")
            self.assertEqual("verified", final.loc["c_verified", "verification_status"])
            self.assertEqual("unsupported", final.loc["c_unsupported", "verification_status"])
            self.assertEqual("contradicted", final.loc["c_contradicted", "verification_status"])
            self.assertEqual("not_evaluated", final.loc["c_not_evaluated", "verification_status"])


if __name__ == "__main__":
    unittest.main()
