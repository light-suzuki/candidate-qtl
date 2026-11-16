"""Sealed C review reaches real Wet without a spare slot or support promotion."""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.core import sha256
from candidate_qtl.route_c_claim_binding import BINDING_COLUMNS, BINDING_CONTRACT
from candidate_qtl.route_c_direct_pipeline import direct_rescue_project
from candidate_qtl.route_c_direct_members import build_route_c_direct_members
from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from tests._issue664_helpers import manifest_for
from tests.test_issue_663_route_c_projection_evidence import RouteCProjectionEvidenceTests


class Fixture(RouteCProjectionEvidenceTests):
    __test__ = False
    # This subclass is a fixture adapter, not a second unittest suite.  The
    # inherited #663 tests otherwise run again during discovery.
    __unittest_skip__ = True
    __unittest_skip_why__ = "fixture adapter; exercised by CReviewRetentionTests"

    def _frozen(self, root, *, include_binding=False):
        frozen, path = super()._frozen(root, include_binding=include_binding)
        data = json.loads(path.read_text())
        if include_binding:
            data["claim_binding"]["complete"] = False
            path.write_text(json.dumps(data), encoding="utf-8")
        return frozen, path

    def _record(self):
        record = super()._record()
        record["evidence_directness"] = "D2_direct_trait_gene_association"
        return record

    def _binding_row(self, record):
        binding = super()._binding_row(record)
        binding.update(binding_status="review_required", phenotype_match_state="bounded_related_trait",
                       perturbation_to_phenotype_binding_state="not_evaluated")
        return binding


class CReviewRetentionTests(unittest.TestCase):
    def _run(self, root, *, count=20, binding=True, reverse=False, conflict=False, mutate=None,
             record_mutate=None, binding_mutate=None, orthology=False, blast=False,
             budget_policy="target", extra_review_routes=(), include_c=True):
        fixture = Fixture()
        if record_mutate:
            original_record = fixture._record()
            record_mutate(original_record)
            fixture._record = lambda: original_record
        if binding_mutate:
            original_binding = fixture._binding_row
            fixture._binding_row = lambda record: binding_mutate(original_binding(record))
        rows = [{"gene_id": f"G{i:04d}", "orthology_reference_gene": ""} for i in range(999)]
        rows.append({"gene_id": "Psat.A", "orthology_reference_gene": "OTHER" if orthology else ""})
        if reverse:
            rows.reverse()
        projected, frozen, freeze, atlas = fixture._project(
            root, rows, [fixture._hit("G0998")] if (conflict or blast) else None, include_binding=binding)
        members = build_route_c_direct_members(str(frozen), str(freeze), projected["projection"],
            projected["manifest"], str(atlas), str(root / "members"), case_id="case1")
        if mutate:
            mutate(root, members, atlas)
        candidates = root / "b.tsv"
        pd.DataFrame([{"candidate_gene_id": f"G{i:04d}", "route": "b",
                       "relevance_evidence_json": json.dumps([{
                           "candidate_gene_id": f"G{i:04d}", "claim_id": f"fixture{i}",
                           "family_id": "fixture-family", "claim_type": "mechanistic_process_bridge",
                           "mechanism": "fixture mechanism", "phenotype_relation": "fixture trait bridge",
                       }]),
                       "support_status": "verified_semantic_support"} for i in range(count)]).to_csv(candidates, sep="\t", index=False)
        b_manifest = manifest_for(candidates, route="b", qtl_id="q1")
        review_paths, review_manifests = [], []
        for index, route in enumerate(extra_review_routes):
            review_path = root / f"route_{index}_review.tsv"
            pd.DataFrame([{
                "candidate_gene_id": "Psat.A", "route": route,
                "support_status": "review_required" if route == "a" else "unverified_route_b_review",
                "rank_contribution": "none", "claim_id": f"fixture-{route}",
                "mechanism": "fixture stem development mechanism",
                "phenotype_relation": "fixture altered stem growth",
            }]).to_csv(review_path, sep="\t", index=False)
            review_manifest = root / f"route_{index}_review_manifest.json"
            review_manifest.write_text(json.dumps({
                "contract": "route-review-candidates-1.0.0", "stage": "review_candidate_freeze",
                "case_id": "case1", "output": str(review_path), "artifact_sha256": sha256(review_path),
                "inputs": {},
            }), encoding="utf-8")
            review_paths.append(str(review_path))
            review_manifests.append(str(review_manifest))
        result = wet_candidate_handoff(str(root / "wet"), case_id="case1",
            route_b_candidates_arg=str(candidates), route_b_candidates_manifest_arg=str(b_manifest),
            route_c_members_arg=members["members"] if include_c else None,
            route_c_manifest_arg=members["manifest"] if include_c else None,
            review_candidates_arg=review_paths or None,
            review_candidates_manifest_arg=review_manifests or None,
            interval_atlas_arg=str(atlas), budget_policy=budget_policy)
        wet = (pd.read_csv(result["wet_candidates"], sep="\t", dtype=str).fillna("")
               if result["wet_candidates"] else pd.DataFrame())
        return wet, result

    def test_same_gene_a_b_c_review_stays_protected_in_boundary_scope(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            wet, result = self._run(root, count=25, extra_review_routes=("a", "b"))
            row = wet.set_index("candidate_gene_id").loc["Psat.A"]
            self.assertEqual({"a", "b", "c"}, {value.strip() for value in row["review_routes"].split("|")})
            self.assertEqual("False", row["route_c_supported"])
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            scope = json.loads(Path(manifest["boundary_scope_manifest"]).read_text(encoding="utf-8"))
            self.assertEqual(["Psat.A"], scope["protected_gene_ids"])
            higher = pd.read_csv(scope["artifacts"]["frozen"]["path"], sep="\t", dtype=str)
            boundary = pd.read_csv(scope["artifacts"]["boundary"]["path"], sep="\t", dtype=str)
            self.assertIn("Psat.A", set(higher["candidate_gene_id"]))
            self.assertNotIn("Psat.A", set(boundary["candidate_gene_id"]))

    def test_same_gene_a_b_c_review_is_counted_once_in_strict_gate(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            wet, _ = self._run(root, count=25, extra_review_routes=("a", "b"), budget_policy="strict")
            gate = json.loads((root / "wet/candidate_budget/candidate_budget_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(1, gate.get("direct_c_review_count", 0))
            self.assertEqual(26, gate["candidate_count_including_route_review"])
            self.assertEqual("False", wet.set_index("candidate_gene_id").loc["Psat.A", "route_c_supported"])

    def test_a_b_reviews_without_accepted_c_binding_are_not_protected(self):
        for options in ({"include_c": False}, {"binding": False}):
            with self.subTest(options=options), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                _, result = self._run(root, count=25, extra_review_routes=("a", "b"), **options)
                manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
                scope = json.loads(Path(manifest["boundary_scope_manifest"]).read_text(encoding="utf-8"))
                self.assertEqual([], scope["protected_gene_ids"])
                review = pd.read_csv(scope["inventory"]["path"], sep="\t", dtype=str).fillna("")
                row = review.set_index("candidate_gene_id").loc["Psat.A"]
                self.assertEqual({"a", "b"}, {value.strip() for value in row["review_routes"].split("|")})

    def test_external_review_artifact_cannot_claim_c_route_token(self):
        for invalid in ("c", "abc", "a | c"):
            with self.subTest(route=invalid), tempfile.TemporaryDirectory() as temp:
                with self.assertRaisesRegex(ValueError, "review candidates route must be a or b"):
                    self._run(Path(temp), extra_review_routes=(invalid,), include_c=False)

    def test_locus_only_blast_source_identity_reaches_canonical_wet_review(self):
        """A BLAST-only locus record keeps its submitted locus provenance through Wet."""
        def locus_only(record):
            record.update(
                reported_gene_id="", reported_gene_symbol="", reference_sequence_id="",
                reported_locus="legacy_locus", gene_identity_state="verified_locus_only",
            )

        def locus_binding(binding):
            binding.update(
                reported_identity_type="locus", reported_identity_value="legacy_locus",
                identity_mapping_basis="exact_source_text",
                reported_identity_verification_state="verified_exact",
            )
            return binding

        with tempfile.TemporaryDirectory() as temp:
            wet, _ = self._run(
                Path(temp), record_mutate=locus_only,
                binding_mutate=locus_binding, blast=True,
            )
            row = wet.set_index("candidate_gene_id").loc["G0998"]
            self.assertEqual("route_review_candidate", row["membership_basis"])
            self.assertEqual("c", row["review_routes"])
            self.assertEqual("False", row["route_c_supported"])

    def test_locus_only_blast_without_binding_is_not_retained(self):
        def locus_only(record):
            record.update(
                reported_gene_id="", reported_gene_symbol="", reference_sequence_id="",
                reported_locus="legacy_locus", gene_identity_state="verified_locus_only",
            )

        with tempfile.TemporaryDirectory() as temp:
            wet, _ = self._run(
                Path(temp), record_mutate=locus_only, binding=False, blast=True,
            )
            self.assertNotIn("G0998", set(wet.candidate_gene_id))

    def test_supported_twenty_does_not_displace_bound_c_review(self):
        with tempfile.TemporaryDirectory() as temp:
            wet, result = self._run(Path(temp))
            self.assertEqual(len(wet), 21)
            review = wet.set_index("candidate_gene_id").loc["Psat.A"]
            self.assertEqual(review["membership_basis"], "route_review_candidate")
            self.assertEqual(review["review_routes"], "c")
            self.assertEqual(review["route_c_supported"], "False")
            self.assertEqual(int((wet["support_tier"] == "tier1_verified_support").sum()), 20)
            self.assertFalse(result.get("adaptive_synthesis_request"))

    def test_row_order_does_not_change_membership(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "a").mkdir(); (root / "b").mkdir()
            first, _ = self._run(root / "a")
            second, _ = self._run(root / "b", reverse=True)
            self.assertEqual(set(first.candidate_gene_id), set(second.candidate_gene_id))

    def test_c_addition_over_twenty_five_does_not_request_count_only_review(self):
        with tempfile.TemporaryDirectory() as temp:
            wet, result = self._run(Path(temp), count=25)
            self.assertEqual(len(wet), 26)
            self.assertFalse(result.get("adaptive_synthesis_request"))

    def test_direct_c_review_stays_in_wet_over_budget_and_boundary(self):
        with tempfile.TemporaryDirectory() as temp:
            wet, result = self._run(
                Path(temp), count=25,
            )
            self.assertEqual(26, len(wet))
            self.assertIn("Psat.A", set(wet.candidate_gene_id))
            row = wet.set_index("candidate_gene_id").loc["Psat.A"]
            self.assertEqual("route_review_candidate", row["membership_basis"])
            self.assertEqual("c", row["review_routes"])
            wet_manifest = json.loads(Path(result["wet_candidates_manifest"]).read_text(encoding="utf-8"))
            self.assertEqual("upper_evidence_scoped_v3", wet_manifest["first_stage_selection_policy_version"])
            trace = pd.read_csv(wet_manifest["candidate_trace"], sep="\t", dtype=str).fillna("")
            trace_row = trace.set_index("candidate_gene_id").loc["Psat.A"]
            self.assertEqual("True", trace_row["explicit_review_member"])
            self.assertEqual("True", trace_row["adopted_first_stage"])
            gate = json.loads((Path(temp) / "wet" / "candidate_budget" / "candidate_budget_manifest.json").read_text())
            self.assertEqual("above_wet_target", gate["status"])
            self.assertFalse(result.get("adaptive_synthesis_request"))

    def test_missing_binding_and_mapping_conflict_are_not_protected(self):
        for options in ({"binding": False}, {"conflict": True}):
            with self.subTest(options=options), tempfile.TemporaryDirectory() as temp:
                wet, _ = self._run(Path(temp), **options)
                self.assertNotIn("Psat.A", set(wet.candidate_gene_id))

    def test_cross_species_and_orthology_only_are_review_protected_not_supported(self):
        def cross_species(record):
            record["species_distance"] = "S2_same_genus"

        with tempfile.TemporaryDirectory() as temp:
            wet, _ = self._run(Path(temp), record_mutate=cross_species)
            row = wet.set_index("candidate_gene_id").loc["Psat.A"]
            self.assertEqual("route_review_candidate", row["membership_basis"])
            self.assertEqual("False", row["route_c_supported"])

        def orthology(record):
            record["reported_gene_id"] = "OTHER"
            record["reference_sequence_id"] = "UniProt:OTHER"

        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            # An orthology reference is an interval projection context, not a
            # same-ID source claim and therefore remains review-only.
            wet, _ = self._run(root, record_mutate=orthology, orthology=True)
            self.assertIn("Psat.A", set(wet.candidate_gene_id))
            row = wet.set_index("candidate_gene_id").loc["Psat.A"]
            self.assertEqual("route_review_candidate", row["membership_basis"])
            self.assertEqual("False", row["route_c_supported"])

    def test_unverified_contradicted_and_variant_only_binding_is_not_protected(self):
        mutations = (
            lambda row: (row.update(reported_variant_state="not_evaluated"), row)[1],
            lambda row: (row.update(binding_status="contradicted"), row)[1],
            lambda row: (row.update(source_organism_state="not_evaluated"), row)[1],
            lambda row: (row.update(identity_mapping_basis="none"), row)[1],
        )
        for mutation in mutations:
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as temp:
                wet, _ = self._run(Path(temp), binding_mutate=mutation)
                self.assertNotIn("Psat.A", set(wet.candidate_gene_id))

    def test_multiple_source_rows_for_one_gene_are_preserved_as_separate_review_records(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            fixture = Fixture()
            record_a = fixture._record()
            record_b = dict(record_a, record_id="r2", source_id="PMID:2", pmid="2",
                             source_title="Independent related-trait report")
            frozen = root / "frozen.tsv"
            pd.DataFrame([record_a, record_b], columns=record_a.keys()).to_csv(frozen, sep="\t", index=False)
            binding = root / "binding.tsv"
            base_fixture = RouteCProjectionEvidenceTests()
            binding_rows = [base_fixture._binding_row(record_a), base_fixture._binding_row(record_b)]
            pd.DataFrame(binding_rows, columns=BINDING_COLUMNS).to_csv(binding, sep="\t", index=False)
            freeze = root / "freeze.json"
            freeze.write_text(json.dumps({"schema_version": "1.0", "command": "route-c-direct-freeze",
                "case_id": "case1", "frozen_sha256": sha256(frozen), "claim_binding": {
                    "contract": BINDING_CONTRACT, "artifact": str(binding), "sha256": sha256(binding),
                    "record_count": 2, "complete": True}}) + "\n", encoding="utf-8")
            atlas = root / "atlas.tsv"
            pd.DataFrame([{"gene_id": "Psat.A", "orthology_reference_gene": ""}]).to_csv(atlas, sep="\t", index=False)
            projected = direct_rescue_project(str(frozen), str(freeze), str(atlas), "q1", str(root / "project"))
            members = build_route_c_direct_members(str(frozen), str(freeze), str(projected["projection"]),
                str(projected["manifest"]), str(atlas), str(root / "members"), case_id="case1")
            frame = pd.read_csv(members["members"], sep="\t", dtype=str)
            self.assertEqual({"r1", "r2"}, set(frame.record_id))
            self.assertEqual({"Psat.A"}, set(frame.candidate_gene_id))
            self.assertEqual(2, len(frame))
            result = wet_candidate_handoff(str(root / "wet"), case_id="case1",
                route_c_members_arg=members["members"], route_c_manifest_arg=members["manifest"],
                interval_atlas_arg=str(atlas))
            wet = pd.read_csv(result["wet_candidates"], sep="\t", dtype=str).fillna("")
            self.assertEqual(1, len(wet))
            self.assertEqual("PMID:1 | PMID:2", wet.iloc[0]["review_source_ids"])
            self.assertEqual("r1 | r2", wet.iloc[0]["review_claim_ids"])
            self.assertEqual("False", wet.iloc[0]["route_c_supported"])

    def test_tampered_binding_is_rejected_at_consumer(self):
        def mutate(root, members, atlas):
            with (root / "binding.tsv").open("a") as stream:
                stream.write("\n")
        with tempfile.TemporaryDirectory() as temp, self.assertRaisesRegex(ValueError, "hash"):
            self._run(Path(temp), mutate=mutate)

    def test_interval_roster_mismatch_is_rejected(self):
        def mutate(root, members, atlas):
            data = pd.read_csv(atlas, sep="\t")
            data.loc[0, "gene_id"] = "different_interval"
            data.to_csv(atlas, sep="\t", index=False)
        with tempfile.TemporaryDirectory() as temp, self.assertRaisesRegex(ValueError, "hash"):
            self._run(Path(temp), mutate=mutate)

    def test_member_atlas_cannot_be_resealed_away_from_projection_generation(self):
        def mutate(root, members, atlas):
            data = pd.read_csv(atlas, sep="\t").fillna("")
            data.loc[0, "orthology_reference_gene"] = "different_reference"
            data.to_csv(atlas, sep="\t", index=False)
            path = Path(members["manifest"])
            manifest = json.loads(path.read_text())
            manifest["inputs"]["atlas"]["sha256"] = sha256(atlas)
            path.write_text(json.dumps(manifest), encoding="utf-8")
        with tempfile.TemporaryDirectory() as temp, self.assertRaisesRegex(ValueError, "interval fingerprint"):
            self._run(Path(temp), mutate=mutate)

    def test_real_route_c_members_cli_entrypoint_writes_sealed_review_artifact(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            fixture = Fixture()
            projected, frozen, freeze, atlas = fixture._project(
                root, [{"gene_id": "Psat.A", "orthology_reference_gene": ""}], None,
                include_binding=True)
            out = root / "cli-members"
            completed = subprocess.run([
                sys.executable, "-m", "candidate_qtl.route_c_direct_members",
                "--frozen", str(frozen), "--freeze-manifest", str(freeze),
                "--projection", str(projected["projection"]),
                "--projection-manifest", str(projected["manifest"]),
                "--atlas", str(atlas), "--out-dir", str(out), "--case-id", "case1",
            ], capture_output=True, text=True, check=False)
            self.assertEqual(0, completed.returncode, completed.stderr)
            payload = json.loads(completed.stdout)
            self.assertTrue(Path(payload["members"]).is_file())
