from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.core import sha256
from candidate_qtl.route_c_claim_binding import (
    BINDING_COLUMNS,
    canonical_frozen_record_sha256,
    load_claim_binding_from_freeze_manifest,
)
from candidate_qtl.route_c_claim_binding_producer import (
    REVIEW_COLUMNS,
    prepare_claim_binding_review,
    seal_claim_binding_review,
)
from candidate_qtl.route_c_direct_rescue import DIRECT_RESCUE_COLUMNS
from candidate_qtl.route_c_direct_pipeline import (
    DIRECT_PROPOSAL_COLUMNS,
    direct_rescue_freeze,
    direct_rescue_project,
)
from candidate_qtl.route_c_direct_members import build_route_c_direct_members
from candidate_qtl.route_c_review_candidates import read_direct_c_review_candidates


class RouteCClaimBindingProducerTests(unittest.TestCase):
    def test_real_freeze_same_paper_separate_quotes_prepare_and_seal(self) -> None:
        for reverse in (False, True):
            with self.subTest(reverse=reverse), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                quotes = {
                    "r1": "GeneA loss of function caused altered growth in Synthetic species.",
                    "r2": "GeneB expression was associated with altered growth in Synthetic species.",
                }
                proposals = []
                extractions = []
                for record_id, gene, relation in (
                    ("r1", "GeneA", "direct_genetic_test"),
                    ("r2", "GeneB", "expression_association"),
                ):
                    proposal = {column: "" for column in DIRECT_PROPOSAL_COLUMNS}
                    proposal.update(
                        record_id=record_id, query_phenotype="altered growth",
                        reported_phenotype="altered growth", focal_species="Synthetic species",
                        reference_species="Synthetic species", reported_gene_id=gene,
                        species_distance="S1_same_species",
                        evidence_directness=("D0_causal_perturbation_or_rescue" if record_id == "r1"
                                             else "D2_direct_trait_gene_association"),
                        evidence_type="synthetic experiment", source_id="PMID:12345678",
                    )
                    proposals.append(proposal)
                    extractions.append({
                        "hypothesis_id": record_id, "citation": "PMID:12345678",
                        "evidence_relation": relation, "quoted_support": quotes[record_id],
                    })
                proposal_path, cache, extraction = root / "proposals.tsv", root / "cache.json", root / "ext.json"
                pd.DataFrame(proposals).to_csv(proposal_path, sep="\t", index=False)
                cache.write_text(json.dumps([{
                    "id": "12345678", "pmid": "12345678", "title": "GeneA and GeneB growth",
                    "abstractText": " ".join(quotes.values()), "pubYear": "2020",
                }]), encoding="utf-8")
                extraction.write_text(json.dumps(list(reversed(extractions)) if reverse else extractions), encoding="utf-8")
                frozen = direct_rescue_freeze(
                    str(proposal_path), str(root / "freeze"), str(cache), claim_extraction_arg=str(extraction),
                )
                prepared = prepare_claim_binding_review(frozen["frozen"], frozen["manifest"], str(root / "prepare"))
                review = pd.read_csv(prepared["review"], sep="\t", dtype=str, keep_default_na=False)
                self.assertEqual(quotes, review.set_index("record_id")["quoted_support"].to_dict())
                # A synthetic reviewer supplies only identity; no fabricated
                # human attestation or perturbation support for the association.
                review["source_organism_state"] = "verified_match"
                review["reported_identity_verification_state"] = "verified_exact"
                review["identity_mapping_basis"] = "exact_source_text"
                review["reviewed_by"] = "synthetic-agent"
                review["binding_reason"] = "identity reviewed; remaining bindings not evaluated"
                reviewed = root / "reviewed.tsv"
                review.to_csv(reviewed, sep="\t", index=False)
                sealed = seal_claim_binding_review(frozen["frozen"], frozen["manifest"], str(reviewed), str(root / "sealed"))
                bindings = pd.read_csv(sealed["binding"], sep="\t", dtype=str, keep_default_na=False)
                self.assertEqual(quotes, bindings.set_index("record_id")["quoted_support"].to_dict())
                self.assertFalse(sealed["complete"])

    def test_literal_na_identity_remains_hash_bound_through_review_consumer(self) -> None:
        """Textual NA is a source identity, including in sealed row hashes."""
        for identity in ("Neutral_G1", "NA"):
            with self.subTest(identity=identity):
                self._literal_identity_roundtrip(identity)

    def _literal_identity_roundtrip(self, identity: str) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            quote = f"The {identity} NA mutation caused altered growth in Synthetic species."
            proposal = {column: "" for column in DIRECT_PROPOSAL_COLUMNS}
            proposal.update(
                record_id="neutral-c1", query_phenotype="altered growth",
                reported_phenotype="altered growth", focal_species="Synthetic species",
                reference_species="Synthetic species", reported_gene_id=identity,
                reported_locus="NA", species_distance="S1_same_species",
                evidence_directness="D0_causal_perturbation_or_rescue",
                evidence_type="neutral experiment", source_id="PMID:12345678",
            )
            proposals, cache, extraction = (
                root / "proposals.tsv", root / "cache.json", root / "extraction.json"
            )
            pd.DataFrame([proposal], columns=DIRECT_PROPOSAL_COLUMNS).to_csv(
                proposals, sep="\t", index=False
            )
            cache.write_text(json.dumps([{
                "id": "12345678", "pmid": "12345678", "title": "Neutral source",
                "abstractText": quote, "pubYear": "2020",
            }]), encoding="utf-8")
            extraction.write_text(json.dumps([{
                "citation": "PMID:12345678", "evidence_relation": "direct_genetic_test",
                "quoted_support": quote,
            }]), encoding="utf-8")
            frozen = direct_rescue_freeze(
                str(proposals), str(root / "freeze"), str(cache),
                claim_extraction_arg=str(extraction),
            )
            prepared = prepare_claim_binding_review(
                frozen["frozen"], frozen["manifest"], str(root / "prepare")
            )
            review = pd.read_csv(prepared["review"], sep="\t", dtype=str, keep_default_na=False)
            self.assertEqual(identity, review.loc[0, "reported_identity_value"])
            review.loc[0, "source_organism_state"] = "verified_match"
            review.loc[0, "reported_identity_verification_state"] = "verified_exact"
            review.loc[0, "identity_mapping_basis"] = "exact_source_text"
            review.loc[0, "reviewed_by"] = "identity-agent"
            review.loc[0, "binding_reason"] = "source identity reviewed; trait binding not evaluated"
            self.assertEqual("", review.loc[0, "curation_method"])
            reviewed = root / "reviewed.tsv"
            review.to_csv(reviewed, sep="\t", index=False)
            sealed = seal_claim_binding_review(
                frozen["frozen"], frozen["manifest"], str(reviewed), str(root / "sealed")
            )
            bound_path = Path(sealed["freeze_manifest"])
            bound_manifest = json.loads(bound_path.read_text(encoding="utf-8"))
            record = pd.read_csv(
                frozen["frozen"], sep="\t", dtype=str, keep_default_na=False
            ).iloc[0].to_dict()
            bindings, generation = load_claim_binding_from_freeze_manifest(
                bound_manifest, bound_path, {"neutral-c1"}, {"neutral-c1": record}
            )
            binding = bindings["neutral-c1"]
            self.assertEqual(identity, binding["reported_identity_value"])
            self.assertEqual(canonical_frozen_record_sha256(record), binding["frozen_record_sha256"])
            self.assertEqual("not_evaluated", binding["perturbation_to_phenotype_binding_state"])
            self.assertEqual("not_evaluated", binding["phenotype_match_state"])
            self.assertFalse(generation["complete"])
            altered_record = {**record, "reported_locus": ""}
            with self.assertRaisesRegex(ValueError, "frozen-row hash mismatch"):
                load_claim_binding_from_freeze_manifest(
                    bound_manifest, bound_path, {"neutral-c1"}, {"neutral-c1": altered_record}
                )
            atlas = root / "atlas.tsv"
            atlas_frame = pd.DataFrame([{
                "gene_id": identity, "chrom": "chr1", "start": "1", "end": "100"
            }])
            atlas_frame.to_csv(atlas, sep="\t", index=False)
            projected = direct_rescue_project(
                frozen["frozen"], str(bound_path), str(atlas), "neutral", str(root / "projection")
            )
            members = build_route_c_direct_members(
                frozen["frozen"], str(bound_path), projected["projection"],
                projected["manifest"], str(atlas), str(root / "members"), case_id="neutral"
            )
            self.assertEqual(0, members["supported_count"])
            candidates, receipt = read_direct_c_review_candidates(
                members["members"], members["manifest"], case_id="neutral",
                atlas_frame=atlas_frame, supported_ids=set(),
            )
            self.assertEqual(1, receipt["review_record_count"])
            self.assertEqual(identity, candidates.iloc[0]["candidate_gene_id"])
            self.assertEqual("unverified_route_review", candidates.iloc[0]["support_status"])

    def _record(self, record_id: str = "r1") -> dict[str, str]:
        row = {column: "" for column in DIRECT_RESCUE_COLUMNS}
        row.update(
            {
                "record_id": record_id,
                "query_phenotype": "fasciated stem",
                "focal_species": "Pisum sativum",
                "reported_phenotype": "fasciated shoot apical meristem",
                "reference_species": "Pisum sativum",
                "reported_gene_id": f"Psat.{record_id}",
                "reported_gene_symbol": "CLV1",
                "species_distance": "S1_same_species",
                "evidence_directness": "D0_causal_perturbation_or_rescue",
                "evidence_type": "loss_of_function",
                "source_id": f"PMID:{100 + int(record_id[1:])}",
                "pmid": str(100 + int(record_id[1:])),
                "source_title": "CLV1 controls shoot apical meristem size",
                "source_year": "2020",
                "source_verification_state": "verified",
                "gene_identity_state": "verified_gene",
                "trait_relation_state": "direct_target_trait",
                "direct_rescue_eligibility": "eligible_direct_record",
                "rank_contribution": "none",
            }
        )
        return row

    def _source_row(
        self, record: dict[str, str], *, accepted: bool = True
    ) -> dict[str, str]:
        return {
            "hypothesis_id": record["record_id"],
            "citation": record["source_id"],
            "verification_status": "verified",
            "citation_verified": "True",
            "retrieval_scope": "title_abstract",
            "extraction_status": "accepted" if accepted else "",
            "quoted_support": (
                "CLV1 loss of function produced a fasciated shoot apical meristem."
                if accepted
                else ""
            ),
        }

    def _case(
        self,
        root: Path,
        records: list[dict[str, str]],
        *,
        accepted: bool = True,
    ):
        frozen = root / "route_c_direct_rescue_frozen.tsv"
        pd.DataFrame(records, columns=DIRECT_RESCUE_COLUMNS).to_csv(
            frozen, sep="\t", index=False, lineterminator="\n"
        )
        source_dir = root / "source_verification"
        source_dir.mkdir()
        validation = source_dir / "source_validation.tsv"
        pd.DataFrame(
            [self._source_row(record, accepted=accepted) for record in records]
        ).to_csv(validation, sep="\t", index=False, lineterminator="\n")
        verified_sources = source_dir / "verified_sources.json"
        verified_sources.write_text(
            json.dumps(
                {
                    record["source_id"]: {
                        "citation": record["source_id"],
                        "title": record["source_title"],
                        "abstractText": self._source_row(record)["quoted_support"],
                    }
                    for record in records
                },
                ensure_ascii=False,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        source_manifest = source_dir / "source_verify_manifest.json"
        source_manifest.write_text(
            json.dumps(
                {
                    "schema_version": "1.1",
                    "command": "source-verify",
                    "outputs": {
                        "source_validation": str(validation),
                        "verified_sources": str(verified_sources),
                    },
                },
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        freeze_manifest = root / "route_c_direct_rescue_freeze_manifest.json"
        freeze_manifest.write_text(
            json.dumps(
                {
                    "schema_version": "1.0",
                    "command": "route-c-direct-freeze",
                    "case_id": "case1",
                    "frozen_sha256": sha256(frozen),
                    "source_verification_manifest": str(source_manifest),
                    "source_validation_sha256": sha256(validation),
                    "verified_sources_sha256": sha256(verified_sources),
                    "outputs": {
                        "frozen": str(frozen),
                        "source_verification": str(source_manifest),
                    },
                },
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        return frozen, freeze_manifest, validation

    @staticmethod
    def _complete_review(review: pd.DataFrame) -> pd.DataFrame:
        review = review.copy()
        review["source_organism_state"] = "verified_match"
        review["reported_identity_verification_state"] = "verified_exact"
        review["identity_mapping_basis"] = "exact_source_text"
        review["experimental_target_gene_state"] = "verified_same_identity"
        review["phenotype_match_state"] = "target_trait"
        review["perturbation_to_phenotype_binding_state"] = "verified_supports"
        review["curation_method"] = "curated_human"
        review["reviewed_by"] = "fixture-reviewer"
        review["binding_reason"] = (
            "reviewer verified gene, organism, perturbation and phenotype"
        )
        return review

    def test_prepare_emits_exact_roster_and_never_marks_complete(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            records = [self._record("r1"), self._record("r2")]
            frozen, manifest, _ = self._case(root, records)
            result = prepare_claim_binding_review(
                str(frozen), str(manifest), str(root / "prepare")
            )
            review = pd.read_csv(
                result["review"], sep="\t", dtype=str
            ).fillna("")
            self.assertEqual(list(REVIEW_COLUMNS), list(review.columns))
            self.assertEqual({"r1", "r2"}, set(review["record_id"]))
            self.assertTrue(review["curation_method"].eq("").all())
            self.assertTrue(
                review["source_organism_state"].eq("not_evaluated").all()
            )
            self.assertTrue(review["quoted_support"].str.len().gt(0).all())

    def test_seal_creates_complete_hash_bound_sidecar_and_new_manifest(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            records = [self._record("r1")]
            frozen, manifest, _ = self._case(root, records)
            prepared = prepare_claim_binding_review(
                str(frozen), str(manifest), str(root / "prepare")
            )
            review = self._complete_review(
                pd.read_csv(
                    prepared["review"], sep="\t", dtype=str
                ).fillna("")
            )
            review_path = root / "reviewed.tsv"
            review.to_csv(
                review_path, sep="\t", index=False, lineterminator="\n"
            )
            sealed = seal_claim_binding_review(
                str(frozen),
                str(manifest),
                str(review_path),
                str(root / "sealed"),
            )
            binding = pd.read_csv(
                sealed["binding"], sep="\t", dtype=str
            ).fillna("")
            bound_manifest = json.loads(
                Path(sealed["freeze_manifest"]).read_text(encoding="utf-8")
            )
            self.assertEqual(list(BINDING_COLUMNS), list(binding.columns))
            self.assertEqual(
                "verified_complete", binding.loc[0, "binding_status"]
            )
            self.assertTrue(sealed["complete"])
            self.assertTrue(bound_manifest["claim_binding"]["complete"])
            self.assertEqual(
                sha256(Path(sealed["binding"])),
                bound_manifest["claim_binding"]["sha256"],
            )
            self.assertNotEqual(str(manifest), sealed["freeze_manifest"])

    def test_agent_label_without_curated_human_attestation_stays_review_required(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            records = [self._record("r1")]
            frozen, manifest, _ = self._case(root, records)
            prepared = prepare_claim_binding_review(
                str(frozen), str(manifest), str(root / "prepare")
            )
            review = self._complete_review(
                pd.read_csv(
                    prepared["review"], sep="\t", dtype=str
                ).fillna("")
            )
            review.loc[0, "curation_method"] = ""
            review.loc[0, "reviewed_by"] = ""
            review_path = root / "unattested.tsv"
            review.to_csv(
                review_path, sep="\t", index=False, lineterminator="\n"
            )
            sealed = seal_claim_binding_review(
                str(frozen),
                str(manifest),
                str(review_path),
                str(root / "sealed"),
            )
            binding = pd.read_csv(
                sealed["binding"], sep="\t", dtype=str
            ).fillna("")
            self.assertEqual(
                "review_required", binding.loc[0, "binding_status"]
            )
            self.assertFalse(sealed["complete"])

    def _seal_with(self, root: Path, method: str, reviewer: str, quote: str | None = None):
        records = [self._record("r1")]
        frozen, manifest, _ = self._case(root, records)
        prepared = prepare_claim_binding_review(
            str(frozen), str(manifest), str(root / "prepare")
        )
        review = self._complete_review(
            pd.read_csv(prepared["review"], sep="\t", dtype=str).fillna("")
        )
        review.loc[0, "curation_method"] = method
        review.loc[0, "reviewed_by"] = reviewer
        if quote is not None:
            review.loc[0, "quoted_support"] = quote
        review_path = root / "reviewed.tsv"
        review.to_csv(review_path, sep="\t", index=False, lineterminator="\n")
        return seal_claim_binding_review(
            str(frozen), str(manifest), str(review_path), str(root / "sealed")
        )

    def test_machine_checked_agent_review_completes_binding(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            sealed = self._seal_with(
                Path(temp), "independent_agent_review", "route-c-reviewer-agent"
            )
            binding = pd.read_csv(sealed["binding"], sep="\t", dtype=str).fillna("")
            manifest = json.loads(
                Path(sealed["freeze_manifest"]).read_text(encoding="utf-8")
            )
            self.assertEqual("verified_complete", binding.loc[0, "binding_status"])
            self.assertTrue(sealed["complete"])
            self.assertEqual(
                {"independent_agent_review": 1},
                manifest["claim_binding_producer"]["curation_method_counts"],
            )

    def test_anonymous_agent_review_stays_review_required(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            sealed = self._seal_with(Path(temp), "independent_agent_review", "")
            binding = pd.read_csv(sealed["binding"], sep="\t", dtype=str).fillna("")
            self.assertEqual("review_required", binding.loc[0, "binding_status"])
            self.assertFalse(sealed["complete"])

    def test_agent_review_cannot_bypass_verbatim_quote_check(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaisesRegex(ValueError, "accepted verbatim extraction"):
                self._seal_with(
                    Path(temp), "independent_agent_review", "route-c-reviewer-agent",
                    quote="A convincing but invented sentence.",
                )

    def test_quote_or_identity_not_bound_to_frozen_sources_is_rejected(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            records = [self._record("r1")]
            frozen, manifest, _ = self._case(root, records)
            prepared = prepare_claim_binding_review(
                str(frozen), str(manifest), str(root / "prepare")
            )
            review = self._complete_review(
                pd.read_csv(
                    prepared["review"], sep="\t", dtype=str
                ).fillna("")
            )
            review.loc[0, "quoted_support"] = (
                "A convincing but invented sentence."
            )
            review_path = root / "bad-quote.tsv"
            review.to_csv(
                review_path, sep="\t", index=False, lineterminator="\n"
            )
            with self.assertRaisesRegex(
                ValueError, "accepted verbatim extraction"
            ):
                seal_claim_binding_review(
                    str(frozen),
                    str(manifest),
                    str(review_path),
                    str(root / "bad1"),
                )

            review = self._complete_review(
                pd.read_csv(
                    prepared["review"], sep="\t", dtype=str
                ).fillna("")
            )
            review.loc[0, "reported_identity_value"] = "Psat.UNRELATED"
            review_path = root / "bad-identity.tsv"
            review.to_csv(
                review_path, sep="\t", index=False, lineterminator="\n"
            )
            with self.assertRaisesRegex(
                ValueError, "not the frozen reported_gene_id"
            ):
                seal_claim_binding_review(
                    str(frozen),
                    str(manifest),
                    str(review_path),
                    str(root / "bad2"),
                )

    def test_roster_and_source_artifact_tampering_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            records = [self._record("r1"), self._record("r2")]
            frozen, manifest, validation = self._case(root, records)
            prepared = prepare_claim_binding_review(
                str(frozen), str(manifest), str(root / "prepare")
            )
            review = self._complete_review(
                pd.read_csv(
                    prepared["review"], sep="\t", dtype=str
                ).fillna("")
            )
            review = review.iloc[:1]
            review_path = root / "missing-row.tsv"
            review.to_csv(
                review_path, sep="\t", index=False, lineterminator="\n"
            )
            with self.assertRaisesRegex(ValueError, "roster mismatch"):
                seal_claim_binding_review(
                    str(frozen),
                    str(manifest),
                    str(review_path),
                    str(root / "bad-roster"),
                )

            validation.write_text(
                validation.read_text(encoding="utf-8") + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                ValueError, "source-validation artifact hash"
            ):
                prepare_claim_binding_review(
                    str(frozen), str(manifest), str(root / "tampered")
                )


if __name__ == "__main__":
    unittest.main()
