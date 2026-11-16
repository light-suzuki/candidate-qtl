from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.core import sha256
from candidate_qtl.route_c_claim_binding import (
    BINDING_COLUMNS,
    BINDING_CONTRACT,
    canonical_frozen_record_sha256,
)
from candidate_qtl.route_c_direct_members import build_route_c_direct_members
from candidate_qtl.route_c_direct_pipeline import direct_rescue_project
from candidate_qtl.route_c_direct_rescue import DIRECT_RESCUE_COLUMNS
from candidate_qtl.blast_evidence import _eligible_interval_hit_rows


class RouteCProjectionEvidenceTests(unittest.TestCase):
    def _record(self) -> dict[str, str]:
        row = {column: "" for column in DIRECT_RESCUE_COLUMNS}
        row.update(
            {
                "record_id": "r1",
                "query_phenotype": "fasciated stem",
                "focal_species": "Pisum sativum",
                "reported_phenotype": "fasciated stem",
                "reference_species": "Pisum sativum",
                "reported_gene_id": "Psat.A",
                "reported_gene_symbol": "GENE_A",
                "species_distance": "S1_same_species",
                "evidence_directness": "D0_causal_perturbation_or_rescue",
                "evidence_type": "loss_of_function",
                "source_id": "PMID:1",
                "pmid": "1",
                "source_title": "Gene A controls stem morphology",
                "source_year": "2020",
                "source_verification_state": "verified",
                "gene_identity_state": "verified_gene",
                "trait_relation_state": "direct_target_trait",
                "reference_sequence_id": "UniProt:A0A000",
                "direct_rescue_eligibility": "eligible_direct_record",
                "rank_contribution": "none",
            }
        )
        return row

    def _binding_row(self, record: dict[str, str]) -> dict[str, str]:
        row = {column: "" for column in BINDING_COLUMNS}
        row.update(
            {
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
                "quoted_support": (
                    "Loss of Gene A produced the target fasciated-stem phenotype."
                ),
                "locator": "abstract",
                "reviewed_by": "issue-663-fixture",
                "binding_status": "verified_complete",
                "binding_reason": (
                    "fixture binds the source gene, perturbation and phenotype"
                ),
                "rank_contribution": "none",
            }
        )
        return row

    def _frozen(
        self, root: Path, *, include_binding: bool = False
    ) -> tuple[Path, Path]:
        record = self._record()
        frozen = root / "frozen.tsv"
        pd.DataFrame([record], columns=DIRECT_RESCUE_COLUMNS).to_csv(
            frozen, sep="\t", index=False, lineterminator="\n"
        )
        manifest = {
            "schema_version": "1.0",
            "command": "route-c-direct-freeze",
            "case_id": "case1",
            "frozen_sha256": sha256(frozen),
        }
        if include_binding:
            binding = root / "binding.tsv"
            pd.DataFrame(
                [self._binding_row(record)], columns=BINDING_COLUMNS
            ).to_csv(binding, sep="\t", index=False, lineterminator="\n")
            manifest["claim_binding"] = {
                "contract": BINDING_CONTRACT,
                "artifact": str(binding),
                "sha256": sha256(binding),
                "record_count": 1,
                "complete": True,
            }
        manifest_path = root / "freeze.json"
        manifest_path.write_text(
            json.dumps(manifest, sort_keys=True) + "\n", encoding="utf-8"
        )
        return frozen, manifest_path

    def _project(
        self,
        root: Path,
        atlas_rows: list[dict[str, str]],
        hit_rows: list[dict[str, str]] | None,
        *,
        include_binding: bool = False,
    ) -> tuple[dict[str, object], Path, Path, Path]:
        frozen, freeze_manifest = self._frozen(
            root, include_binding=include_binding
        )
        atlas = root / "atlas.tsv"
        pd.DataFrame(atlas_rows).to_csv(
            atlas, sep="\t", index=False, lineterminator="\n"
        )
        hits = None
        if hit_rows is not None:
            hits = root / "hits.tsv"
            pd.DataFrame(hit_rows).to_csv(
                hits, sep="\t", index=False, lineterminator="\n"
            )
        result = direct_rescue_project(
            str(frozen),
            str(freeze_manifest),
            str(atlas),
            "q1",
            str(root / "project"),
            str(hits) if hits else None,
        )
        return result, frozen, freeze_manifest, atlas

    @staticmethod
    def _hit(
        target: str, *, panel: str = "r1", identity: str = "60"
    ) -> dict[str, str]:
        return {
            "panel_id": panel,
            "target_id": target,
            "eligibility": "eligible_interval_hit",
            "pident": identity,
            "query_coverage": "0.80",
            "target_coverage": "0.90",
        }

    def test_exact_and_blast_same_gene_are_corroboration_and_exact_wins(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            result, *_ = self._project(
                root,
                [{"gene_id": "Psat.A", "orthology_reference_gene": ""}],
                [self._hit("Psat.A")],
            )
            projection = pd.read_csv(
                result["projection"], sep="\t", dtype=str
            ).fillna("")
            evidence = pd.read_csv(
                result["projection_evidence"], sep="\t", dtype=str
            ).fillna("")
            manifest = json.loads(
                Path(result["manifest"]).read_text(encoding="utf-8")
            )
            self.assertEqual("exact_focal_id", projection.loc[0, "projection_basis"])
            mapped = evidence[evidence["focal_gene_id"].ne("")]
            self.assertEqual(
                {"exact_focal_id", "blastp"}, set(mapped["mapping_method"])
            )
            self.assertEqual(
                {"corroborating_methods"}, set(evidence["conflict_state"])
            )
            self.assertIn(
                "all_projection_bases=exact_focal_id;blastp",
                projection.loc[0, "notes"],
            )
            self.assertEqual(0, manifest["mapping_conflict_count"])

    def test_exact_protein_id_does_not_authorize_a_different_bound_gene(self) -> None:
        """The source gene is bound, but only its protein matched the atlas ID."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            result, frozen, freeze, atlas = self._project(
                root, [{"gene_id": "UniProt:A0A000"}],
                [self._hit("UniProt:A0A000")], include_binding=True,
            )
            projection = pd.read_csv(result["projection"], sep="\t", dtype=str)
            self.assertEqual("exact_focal_id", projection.iloc[0].projection_basis)
            members = build_route_c_direct_members(
                str(frozen), str(freeze), result["projection"], result["manifest"],
                str(atlas), str(root / "members"), case_id="case1",
            )
            self.assertEqual(0, members["supported_count"])
            self.assertEqual(1, members["review_count"])

    def test_exact_vs_blast_different_targets_is_explicit_and_blocks_support(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            result, frozen, freeze_manifest, atlas = self._project(
                root,
                [
                    {"gene_id": "Psat.B", "orthology_reference_gene": ""},
                    {"gene_id": "Psat.A", "orthology_reference_gene": ""},
                ],
                [self._hit("Psat.B")],
                include_binding=True,
            )
            projection = pd.read_csv(
                result["projection"], sep="\t", dtype=str
            ).fillna("")
            manifest = json.loads(
                Path(result["manifest"]).read_text(encoding="utf-8")
            )
            self.assertEqual({"Psat.A", "Psat.B"}, set(projection["focal_gene_id"]))
            self.assertEqual({"one_to_many"}, set(projection["ambiguity_state"]))
            self.assertEqual(["r1"], manifest["mapping_conflict_record_ids"])

            members = build_route_c_direct_members(
                str(frozen),
                str(freeze_manifest),
                str(result["projection"]),
                str(result["manifest"]),
                str(atlas),
                str(root / "members"),
                case_id="case1",
            )
            frame = pd.read_csv(
                members["members"], sep="\t", dtype=str
            ).fillna("")
            self.assertTrue(
                frame["support_status"].eq("reference_phenotype_review").all()
            )
            self.assertTrue(frame["membership_effect"].eq("review_only").all())
            self.assertEqual(1, members["mapping_conflict_count"])
            self.assertEqual(
                1, members["projection_conflict_support_blocked_count"]
            )

    def test_unrelated_global_hits_file_is_not_evaluated_for_this_record(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            result, *_ = self._project(
                root,
                [{"gene_id": "Psat.other", "orthology_reference_gene": "OTHER"}],
                [self._hit("Psat.other", panel="another-record")],
            )
            evidence = pd.read_csv(
                result["projection_evidence"], sep="\t", dtype=str
            ).fillna("")
            blast = evidence[evidence["mapping_method"].eq("blastp")].iloc[0]
            projection = pd.read_csv(
                result["projection"], sep="\t", dtype=str
            ).fillna("")
            self.assertEqual(
                "not_evaluated_no_panel", blast["method_evaluation_state"]
            )
            self.assertEqual(
                "orthology_reference_gene", projection.loc[0, "projection_basis"]
            )
            self.assertEqual(
                "no_interval_match", projection.loc[0, "projection_state"]
            )

    def test_canonical_blast_producer_rows_project_and_remain_review_only(self):
        """The producer's ``eligible`` marker is consumed without support promotion."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            result, frozen, freeze_manifest, atlas = self._project(
                root,
                [{"gene_id": "Psat.B", "orthology_reference_gene": ""}],
                _eligible_interval_hit_rows(
                    "r1",
                    [{
                        "sseqid": "Psat.B",
                        "pident": 60.0,
                        "query_coverage": 0.8,
                        "target_coverage": 0.9,
                        "evalue": 1e-10,
                        "bitscore": 100.0,
                        "hsp_count": 1,
                        "qstart": 1,
                        "qend": 80,
                        "sstart": 1,
                        "send": 90,
                        "query_union_segments": "1-80",
                        "target_union_segments": "1-90",
                    }],
                    "unique",
                    review_status="review_required",
                    projection_eligibility="review_only",
                )
                + [{
                    "panel_id": "r1",
                    "target_id": "Psat.B",
                    "eligibility": "invalid",
                    "pident": "99",
                    "query_coverage": "1",
                    "target_coverage": "1",
                }],
                include_binding=True,
            )
            projection = pd.read_csv(
                result["projection"], sep="\t", dtype=str
            ).fillna("")
            evidence = pd.read_csv(
                result["projection_evidence"], sep="\t", dtype=str
            ).fillna("")
            self.assertEqual("blastp", projection.loc[0, "projection_basis"])
            blast = evidence[evidence["mapping_method"].eq("blastp")]
            self.assertEqual(1, len(blast))
            self.assertEqual("Psat.B", blast.iloc[0]["focal_gene_id"])

            members = build_route_c_direct_members(
                str(frozen),
                str(freeze_manifest),
                str(result["projection"]),
                str(result["manifest"]),
                str(atlas),
                str(root / "members"),
                case_id="case1",
            )
            frame = pd.read_csv(
                members["members"], sep="\t", dtype=str
            ).fillna("")
            self.assertEqual("verified_complete", frame.loc[0, "claim_binding_status"])
            self.assertEqual("reference_phenotype_review", frame.loc[0, "support_status"])
            self.assertEqual("review_only", frame.loc[0, "membership_effect"])

    def test_projection_and_evidence_are_byte_stable_under_input_row_permutation(
        self,
    ) -> None:
        atlas_rows = [
            {"gene_id": "Psat.B", "orthology_reference_gene": "OTHER"},
            {"gene_id": "Psat.A", "orthology_reference_gene": ""},
        ]
        hit_rows = [
            self._hit("Psat.B", identity="61"),
            self._hit("Psat.A", identity="59"),
        ]
        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            result_a, *_ = self._project(Path(first), atlas_rows, hit_rows)
            result_b, *_ = self._project(
                Path(second), list(reversed(atlas_rows)), list(reversed(hit_rows))
            )
            self.assertEqual(
                Path(result_a["projection"]).read_bytes(),
                Path(result_b["projection"]).read_bytes(),
            )
            self.assertEqual(
                Path(result_a["projection_evidence"]).read_bytes(),
                Path(result_b["projection_evidence"]).read_bytes(),
            )


if __name__ == "__main__":
    unittest.main()
