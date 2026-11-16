import json
import contextlib
import io
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.core import sha256
from candidate_qtl.cli import main
from candidate_qtl.route_c_claim_binding import (
    BINDING_COLUMNS,
    BINDING_CONTRACT,
    canonical_frozen_record_sha256,
)
from candidate_qtl.route_c_direct_members import (
    OUTPUT_COLUMNS,
    build_route_c_direct_members,
)
from candidate_qtl.route_c_direct_rescue import DIRECT_RESCUE_COLUMNS, PROJECTION_COLUMNS


class RouteCDirectMembersTests(unittest.TestCase):
    def _binding_row(self, record):
        row = {column: "" for column in BINDING_COLUMNS}
        row.update({
            "record_id": record["record_id"],
            "frozen_record_sha256": canonical_frozen_record_sha256(record),
            "citation": record["source_id"],
            "citation_resolution_state": "verified",
            "source_organism_state": "verified_match",
            "source_accession_or_cultivar_state": "not_applicable",
            "reported_identity_type": "stable_gene_id",
            "reported_identity_value": record["reported_gene_id"],
            "reported_identity_verification_state": "verified_exact",
            "identity_mapping_basis": "exact_source_text",
            "experimental_target_gene_state": "verified_same_identity",
            "reported_variant_state": "not_applicable",
            "phenotype_match_state": "target_trait",
            "perturbation_to_phenotype_binding_state": "verified_supports",
            "source_scope": "abstract",
            "source_record_sha256": "a" * 64,
            "quoted_support": "The reported gene perturbation produced the target phenotype.",
            "locator": "abstract",
            "reviewed_by": "fixture-reviewer",
            "binding_status": "verified_complete",
            "binding_reason": "synthetic fixture binds gene, species, perturbation and phenotype",
            "rank_contribution": "none",
        })
        return row

    def _run(self, records, projections, atlas, *, include_binding=True, tamper_binding=False,
             binding_rows=None, via_cli=False):
        root = Path(tempfile.mkdtemp())
        frozen = root / "frozen.tsv"
        freeze_manifest = root / "freeze.json"
        projection = root / "projection.tsv"
        projection_manifest = root / "projection.json"
        atlas_path = root / "atlas.tsv"
        pd.DataFrame(records, columns=DIRECT_RESCUE_COLUMNS).to_csv(
            frozen, sep="\t", index=False
        )
        freeze_payload = {"frozen_sha256": sha256(frozen), "case_id": "case1"}
        if include_binding:
            binding = root / "binding.tsv"
            pd.DataFrame(
                binding_rows or [self._binding_row(record) for record in records],
                columns=BINDING_COLUMNS,
            ).to_csv(binding, sep="\t", index=False)
            binding_sha = sha256(binding)
            if tamper_binding:
                binding_sha = "0" * 64
            freeze_payload["claim_binding"] = {
                "contract": BINDING_CONTRACT,
                "artifact": str(binding),
                "sha256": binding_sha,
                "record_count": len(records),
                "complete": True,
            }
        freeze_manifest.write_text(json.dumps(freeze_payload))
        pd.DataFrame(projections, columns=PROJECTION_COLUMNS).to_csv(
            projection, sep="\t", index=False
        )
        projection_manifest.write_text(
            json.dumps({"projection_sha256": sha256(projection), "qtl_id": "q1"})
        )
        pd.DataFrame(atlas).to_csv(atlas_path, sep="\t", index=False)
        if via_cli:
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(
                    0,
                    main([
                        "route-c-direct-members",
                        "--frozen", str(frozen),
                        "--freeze-manifest", str(freeze_manifest),
                        "--projection", str(projection),
                        "--projection-manifest", str(projection_manifest),
                        "--atlas", str(atlas_path),
                        "--out-dir", str(root / "out"),
                        "--case-id", "case1",
                        "--json",
                    ]),
                )
            result = {
                "members": str(root / "out" / "route_c_direct_members.tsv"),
                "manifest": str(root / "out" / "route_c_direct_members_manifest.json"),
            }
        else:
            result = build_route_c_direct_members(
                str(frozen),
                str(freeze_manifest),
                str(projection),
                str(projection_manifest),
                str(atlas_path),
                str(root / "out"),
                case_id="case1",
            )
        return (
            pd.read_csv(result["members"], sep="\t", dtype=str).fillna(""),
            result,
        )

    def _record(
        self,
        record_id,
        distance="S1_same_species",
        directness="D0_causal_perturbation_or_rescue",
    ):
        return {
            "record_id": record_id,
            "query_phenotype": "short roots",
            "focal_species": "Pisum sativum",
            "focal_cultivar_or_accession": "",
            "reported_phenotype": "short roots",
            "reference_species": (
                "Pisum sativum" if distance.startswith("S1") else "Arabidopsis thaliana"
            ),
            "reference_cultivar_or_accession": "",
            "reported_gene_id": "Psat.X",
            "reported_gene_symbol": "GENE",
            "reported_locus": "",
            "reported_variant": "",
            "species_distance": distance,
            "evidence_directness": directness,
            "evidence_type": "mutant",
            "source_id": "PMID:1",
            "pmid": "1",
            "pmcid": "",
            "doi": "",
            "source_title": "paper",
            "source_year": "2020",
            "source_verification_state": "verified",
            "gene_identity_state": "verified_gene",
            "trait_relation_state": "direct_target_trait",
            "reference_sequence_id": "",
            "direct_rescue_eligibility": "eligible_direct_record",
            "rank_contribution": "none",
            "notes": "",
        }

    def _projection(
        self,
        record_id,
        gene="Psat.X",
        basis="exact_focal_id",
        ambiguity="unique",
        state="inside_current_qtl",
    ):
        if state == "no_interval_match":
            ambiguity = "no_threshold_hit"
            basis = "orthology_reference_gene"
        return {
            "record_id": record_id,
            "qtl_id": "q1" if state != "not_projected" else "",
            "focal_gene_id": gene if state == "inside_current_qtl" else "",
            "projection_state": state,
            "projection_basis": basis if state != "not_projected" else "not_evaluated",
            "source_reference_gene_id": "Psat.X",
            "identity": "",
            "query_coverage": "",
            "target_coverage": "",
            "ambiguity_state": ambiguity if state != "not_projected" else "not_evaluated",
            "interval_universe_fingerprint": "sha256:x",
            "rank_contribution": "none",
            "notes": "",
        }

    def test_same_species_d0_supports_and_preserves_metadata(self):
        records = [self._record("r1")]
        projections = [self._projection("r1")]
        frame, result = self._run(
            records,
            projections,
            [{"gene_id": "Psat.X", "family_id": "F1", "orthogroup_id": "OG1"}],
        )
        self.assertEqual(list(frame.columns), list(OUTPUT_COLUMNS))
        self.assertEqual(frame.loc[0, "support_status"], "verified_direct_support")
        self.assertEqual(frame.loc[0, "membership_effect"], "supported_union")
        self.assertEqual(frame.loc[0, "claim_binding_status"], "verified_complete")
        self.assertEqual(frame.loc[0, "family_ids"], "F1")
        self.assertEqual(result["supported_count"], 1)

    def test_production_cli_forwards_case_id_to_keyword_only_builder(self):
        frame, result = self._run(
            [self._record("r1")],
            [self._projection("r1")],
            [{"gene_id": "Psat.X", "family_id": "F1", "orthogroup_id": "OG1"}],
            via_cli=True,
        )
        manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
        self.assertEqual("case1", manifest["case_id"])
        self.assertEqual("verified_direct_support", frame.loc[0, "support_status"])

    def test_same_species_exact_without_structured_binding_is_review_only(self):
        records = [self._record("r1")]
        projections = [self._projection("r1")]
        frame, result = self._run(
            records,
            projections,
            [{"gene_id": "Psat.X", "family_id": "F1"}],
            include_binding=False,
        )
        self.assertEqual(frame.loc[0, "support_status"], "reference_phenotype_review")
        self.assertEqual(frame.loc[0, "membership_effect"], "review_only")
        self.assertEqual(frame.loc[0, "claim_binding_status"], "not_evaluated")
        manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
        self.assertEqual(0, manifest["supported_count"])
        self.assertEqual(1, manifest["support_blocked_incomplete_claim_binding_count"])
        self.assertEqual(
            "not_evaluated_missing_binding",
            manifest["claim_binding"]["generation_state"],
        )

    def test_corrupt_supplied_binding_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "claim-binding artifact hash"):
            self._run(
                [self._record("r1")],
                [self._projection("r1")],
                [{"gene_id": "Psat.X"}],
                tamper_binding=True,
            )

    def test_binding_rows_cannot_be_swapped_between_frozen_records(self):
        records = [self._record("r1"), self._record("r2")]
        rows = [self._binding_row(record) for record in records]
        # Keep each sidecar record_id but reuse the other frozen row's hash.
        rows[0]["frozen_record_sha256"], rows[1]["frozen_record_sha256"] = (
            rows[1]["frozen_record_sha256"], rows[0]["frozen_record_sha256"]
        )
        with self.assertRaisesRegex(ValueError, "frozen-row hash mismatch"):
            self._run(
                records,
                [self._projection("r1"), self._projection("r2")],
                [{"gene_id": "Psat.X"}],
                binding_rows=rows,
            )

    def test_frozen_source_or_trait_edit_invalidates_binding(self):
        records = [self._record("r1")]
        rows = [self._binding_row(records[0])]
        records[0]["source_id"] = "PMID:changed"
        with self.assertRaisesRegex(ValueError, "frozen-row hash mismatch"):
            self._run(
                records,
                [self._projection("r1")],
                [{"gene_id": "Psat.X"}],
                binding_rows=rows,
            )

    def test_cross_species_and_d2_are_review_only(self):
        records = [
            self._record("r1", "S5_model_or_other_plant"),
            self._record("r2", directness="D2_direct_trait_gene_association"),
        ]
        projections = [
            self._projection("r1", basis="orthology_reference_gene"),
            self._projection("r2"),
        ]
        frame, _ = self._run(records, projections, [{"gene_id": "Psat.X"}])
        self.assertTrue(frame["support_status"].eq("reference_phenotype_review").all())
        self.assertTrue(frame["membership_effect"].eq("review_only").all())

    def test_one_to_many_retains_every_copy(self):
        records = [self._record("r1")]
        projections = [
            self._projection("r1", "Psat.X1", ambiguity="one_to_many"),
            self._projection("r1", "Psat.X2", ambiguity="one_to_many"),
        ]
        frame, result = self._run(
            records,
            projections,
            [{"gene_id": "Psat.X1"}, {"gene_id": "Psat.X2"}],
        )
        self.assertEqual(list(frame["candidate_gene_id"]), ["Psat.X1", "Psat.X2"])
        self.assertTrue(frame["direct_rescue_eligibility"].eq("eligible_direct_record").all())
        self.assertTrue(frame["projection_state"].eq("inside_current_qtl").all())
        self.assertTrue(frame["projection_basis"].eq("exact_focal_id").all())
        self.assertTrue(frame["projection_ambiguity_state"].eq("one_to_many").all())
        self.assertTrue(frame["claim_binding_status"].eq("verified_complete").all())
        self.assertTrue(frame["claim_binding_identity_state"].eq("verified_exact").all())
        # The record's claim is bound to Psat.X, not either resolved copy;
        # keep both candidates visible but do not promote the ambiguous map.
        self.assertTrue(frame["claim_binding_projection_authorized"].eq("False").all())
        self.assertTrue(frame["support_status"].eq("reference_phenotype_review").all())
        self.assertTrue(frame["membership_effect"].eq("review_only").all())
        self.assertTrue(frame["rank_contribution"].eq("none").all())
        self.assertEqual(result["record_count"], 2)

    def test_no_interval_match_stays_ledger_only(self):
        frame, _ = self._run(
            [self._record("r1", "S5_model_or_other_plant")],
            [self._projection("r1", state="no_interval_match")],
            [{"gene_id": "Psat.other"}],
        )
        self.assertEqual(frame.loc[0, "support_status"], "ledger_only")
        self.assertEqual(frame.loc[0, "membership_effect"], "ledger_only")

    def test_hash_mismatch_is_rejected(self):
        root = Path(tempfile.mkdtemp())
        frozen = root / "frozen.tsv"
        frozen.write_text("record_id\nnot-used\n")
        manifest = root / "freeze.json"
        manifest.write_text(json.dumps({"frozen_sha256": "wrong"}))
        with self.assertRaisesRegex(ValueError, "hash"):
            build_route_c_direct_members(
                str(frozen),
                str(manifest),
                str(frozen),
                str(manifest),
                str(frozen),
                str(root / "out"),
            )


if __name__ == "__main__":
    unittest.main()
