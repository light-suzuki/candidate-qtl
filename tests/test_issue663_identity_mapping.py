"""Typed source gene -> protein -> native BLAST ID reaches review/Wet only."""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.core import sha256
from candidate_qtl.route_c_direct_members import build_route_c_direct_members
from candidate_qtl.route_c_direct_projection import direct_rescue_project
from candidate_qtl.route_c_identity_mapping import IDENTITY_MAPPING_COLUMNS
from candidate_qtl.route_c_review_candidates import read_direct_c_review_candidates
from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from candidate_qtl.wet_candidate_handoff import _read_optional, _tier_interval_candidates, TIER_DARK
from tests._issue664_helpers import manifest_for
from tests import test_issue722_c_review_retention as fixtures


class IdentityMappingTests(unittest.TestCase):
    def test_optional_reader_preserves_ids_without_turning_missing_material_into_tier_evidence(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "atlas.tsv"
            pd.DataFrame([{"gene_id": "NA", "family_id": "NA", "domain_ids": "N/A",
                           "flag": "False", "count": "0", "empty_count": ""}]).to_csv(
                path, sep="\t", index=False)
            atlas = _read_optional(str(path))
            self.assertEqual("NA", atlas.iloc[0].gene_id)
            self.assertEqual("", atlas.iloc[0].family_id)
            self.assertEqual("", atlas.iloc[0].domain_ids)
            self.assertEqual("False", atlas.iloc[0].flag)
            self.assertEqual("0", atlas.iloc[0]["count"])
            self.assertEqual("", atlas.iloc[0].empty_count)
            reference = Path(temp) / "reference.tsv"
            pd.DataFrame([{"reference_gene_id": "NA", "gene_id": "NA",
                           "mutant_phenotype": "NA"}]).to_csv(reference, sep="\t", index=False)
            functions = _read_optional(str(reference))
            self.assertEqual("NA", functions.iloc[0].reference_gene_id)
            self.assertEqual("", functions.iloc[0].mutant_phenotype)
            self.assertEqual({"NA": TIER_DARK}, _tier_interval_candidates(atlas, functions))

    def _inputs(self, root, *, mapping=True, mutate=None, conflict=False,
                multiple_copies=False, complete_binding=False):
        fixture = fixtures.Fixture()
        record = fixture._record()
        # Synthetic identity fixture; no scientific Type1/2/3 performance claim.
        record.update(reported_gene_id="AT3G27920", reported_gene_symbol="",
                      reference_sequence_id="P27900", reference_species="Arabidopsis thaliana",
                      focal_species="Arabidopsis thaliana")
        fixture._record = lambda: record
        if complete_binding:
            record["evidence_directness"] = "D0_causal_perturbation_or_rescue"
            original_binding = fixture._binding_row
            def binding_row(record):
                row = original_binding(record)
                row.update(binding_status="verified_complete", phenotype_match_state="target_trait",
                           perturbation_to_phenotype_binding_state="verified_supports")
                return row
            fixture._binding_row = binding_row
        frozen, freeze = fixture._frozen(root, include_binding=True)
        if complete_binding:
            manifest = json.loads(freeze.read_text())
            manifest["claim_binding"]["complete"] = True
            freeze.write_text(json.dumps(manifest), encoding="utf-8")
        native = "gene:AT3G27920"
        atlas = root / "atlas.tsv"
        atlas_rows = [{"gene_id": native}] + [{"gene_id": f"G{i}"} for i in range(20)]
        if conflict:
            for row in atlas_rows:
                row["orthology_reference_gene"] = "AT3G27920" if row["gene_id"] == "G19" else ""
        pd.DataFrame(atlas_rows).to_csv(
            atlas, sep="\t", index=False)
        hits = root / "hits.tsv"
        rows = [fixture._hit(native, identity="100")]
        if multiple_copies:
            rows.append(fixture._hit("G19"))
        pd.DataFrame(rows).to_csv(hits, sep="\t", index=False)
        payload = {"primaryAccession": "P27900",
                   "organism": {"scientificName": "Arabidopsis thaliana"},
                   "uniProtKBCrossReferences": [{"database": "TAIR", "id": "AT3G27920"}]}
        mapping_row = dict(zip(IDENTITY_MAPPING_COLUMNS, [
            "r1", "reported_gene_id", "AT3G27920", "reference_sequence_id", "P27900",
            "Arabidopsis thaliana", "TAIR", "curated_id_map", "uniprot.json", "pending",
            "pending", "fixture-review", "none",
        ]))
        if mutate:
            mutate(payload, mapping_row)
        source = root / "uniprot.json"
        source.write_text(json.dumps(payload, sort_keys=True), encoding="utf-8")
        mapping_row.update(source_record_sha256=sha256(source), quoted_mapping=source.read_text())
        mapping_path = root / "identity.tsv"
        pd.DataFrame([mapping_row], columns=IDENTITY_MAPPING_COLUMNS).to_csv(mapping_path, sep="\t", index=False)
        return frozen, freeze, atlas, hits, mapping_path if mapping else None

    def _project(self, root, *, cli=False, **options):
        frozen, freeze, atlas, hits, mapping = self._inputs(root, **options)
        if cli:
            command = [sys.executable, "-m", "candidate_qtl", "route-c-direct-project",
                       "--frozen", str(frozen), "--freeze-manifest", str(freeze),
                       "--atlas", str(atlas), "--qtl-id", "q1", "--out-dir", str(root / "project"),
                       "--eligible-hits", str(hits), "--identity-mapping", str(mapping), "--json"]
            output = subprocess.run(command, text=True, capture_output=True, check=True)
            projected = json.loads(output.stdout)
        else:
            projected = direct_rescue_project(str(frozen), str(freeze), str(atlas), "q1",
                                             str(root / "project"), str(hits),
                                             identity_mapping_arg=str(mapping) if mapping else None)
        members = build_route_c_direct_members(str(frozen), str(freeze), projected["projection"],
            projected["manifest"], str(atlas), str(root / "members"), case_id="case1")
        return projected, members, atlas

    def _read(self, members, atlas):
        return read_direct_c_review_candidates(members["members"], members["manifest"],
            case_id="case1", atlas_frame=pd.read_csv(atlas, sep="\t", dtype=str, keep_default_na=False),
            supported_ids=set())[0]

    def test_cli_projection_members_reader_and_canonical_wet_keep_native_review(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            projected, members, atlas = self._project(root, cli=True)
            self.assertEqual(0, members["supported_count"])
            self.assertEqual(["gene:AT3G27920"], self._read(members, atlas).candidate_gene_id.tolist())
            evidence = pd.read_csv(projected["projection_evidence"], sep="\t").fillna("")
            hit = evidence[evidence.mapping_method.eq("blastp")].iloc[0]
            self.assertEqual("P27900", hit.source_identity_value)
            self.assertEqual("reference_sequence_id", hit.source_identity_field)
            self.assertEqual("verified_structured_sequence_mapping", hit.source_identity_verification_state)
            candidates = root / "b.tsv"
            pd.DataFrame([{"candidate_gene_id": f"G{i}", "route": "b",
                           "support_status": "verified_semantic_support"} for i in range(20)]).to_csv(
                candidates, sep="\t", index=False)
            result = wet_candidate_handoff(str(root / "wet"), case_id="case1",
                route_b_candidates_arg=str(candidates),
                route_b_candidates_manifest_arg=str(manifest_for(candidates, route="b", qtl_id="q1")),
                route_c_members_arg=members["members"], route_c_manifest_arg=members["manifest"],
                interval_atlas_arg=str(atlas))
            wet = pd.read_csv(result["wet_candidates"], sep="\t", dtype=str).fillna("")
            self.assertEqual(21, len(wet))
            row = wet.set_index("candidate_gene_id").loc["gene:AT3G27920"]
            self.assertEqual("c", row.review_routes)
            self.assertEqual("False", row.route_c_supported)
            self.assertEqual("route_review_candidate", row.membership_basis)

    def test_missing_mapping_does_not_guess_prefix_or_protein_identity(self):
        with tempfile.TemporaryDirectory() as temp:
            projected, members, atlas = self._project(Path(temp), mapping=False)
            self.assertEqual(1, members["review_count"])
            self.assertTrue(self._read(members, atlas).empty)
            self.assertEqual(0, members["supported_count"])

    def test_mapping_cannot_resolve_conflicting_blast_targets(self):
        with tempfile.TemporaryDirectory() as temp:
            _, members, atlas = self._project(Path(temp), conflict=True)
            self.assertTrue(self._read(members, atlas).empty)
            self.assertEqual(0, members["supported_count"])

    def test_multiple_blast_copies_are_retained_without_support_promotion(self):
        with tempfile.TemporaryDirectory() as temp:
            _, members, atlas = self._project(Path(temp), multiple_copies=True, complete_binding=True)
            self.assertEqual({"gene:AT3G27920", "G19"}, set(self._read(members, atlas).candidate_gene_id))
            self.assertEqual(0, members["supported_count"])

    def test_exact_protein_endpoint_cannot_authorize_unmatched_source_gene_support(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            frozen, freeze, atlas, hits, mapping = self._inputs(root, complete_binding=True)
            pd.DataFrame([{"gene_id": "P27900"}]).to_csv(atlas, sep="\t", index=False)
            pd.DataFrame([fixtures.Fixture._hit("P27900")]).to_csv(hits, sep="\t", index=False)
            projected = direct_rescue_project(str(frozen), str(freeze), str(atlas), "q1",
                str(root / "project"), str(hits), identity_mapping_arg=str(mapping))
            members = build_route_c_direct_members(str(frozen), str(freeze), projected["projection"],
                projected["manifest"], str(atlas), str(root / "members"), case_id="case1")
            self.assertEqual(0, members["supported_count"])
            self.assertEqual(["P27900"], self._read(members, atlas).candidate_gene_id.tolist())

    def test_literal_na_blast_target_reaches_review_without_becoming_missing(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            frozen, freeze, atlas, hits, mapping = self._inputs(root)
            atlas_frame = pd.read_csv(atlas, sep="\t", dtype=str, keep_default_na=False)
            atlas_frame["gene_id"] = atlas_frame["gene_id"].replace("gene:AT3G27920", "NA")
            atlas_frame.to_csv(atlas, sep="\t", index=False)
            hit_frame = pd.read_csv(hits, sep="\t", dtype=str, keep_default_na=False)
            hit_frame["target_id"] = "NA"
            hit_frame.to_csv(hits, sep="\t", index=False)
            projected = direct_rescue_project(str(frozen), str(freeze), str(atlas), "q1",
                str(root / "project"), str(hits), identity_mapping_arg=str(mapping))
            self.assertEqual(1, projected["inside_count"])
            members = build_route_c_direct_members(str(frozen), str(freeze), projected["projection"],
                projected["manifest"], str(atlas), str(root / "members"), case_id="case1")
            self.assertEqual(0, members["supported_count"])
            self.assertEqual(["NA"], self._read(members, atlas).candidate_gene_id.tolist())

            candidates = root / "b.tsv"
            pd.DataFrame([{"candidate_gene_id": f"G{i}", "route": "b",
                           "support_status": "verified_semantic_support"} for i in range(20)]).to_csv(
                candidates, sep="\t", index=False)
            result = wet_candidate_handoff(str(root / "wet"), case_id="case1",
                route_b_candidates_arg=str(candidates),
                route_b_candidates_manifest_arg=str(manifest_for(candidates, route="b", qtl_id="q1")),
                route_c_members_arg=members["members"], route_c_manifest_arg=members["manifest"],
                interval_atlas_arg=str(atlas))
            wet = pd.read_csv(result["wet_candidates"], sep="\t", dtype=str, keep_default_na=False)
            self.assertEqual(21, len(wet))
            row = wet.set_index("candidate_gene_id").loc["NA"]
            self.assertEqual("c", row.review_routes)
            self.assertEqual("False", row.route_c_supported)
            self.assertEqual("route_review_candidate", row.membership_basis)

            # A genuinely empty target remains missing and does not create a hit.
            hit_frame["target_id"] = ""
            hit_frame.to_csv(hits, sep="\t", index=False)
            empty = direct_rescue_project(str(frozen), str(freeze), str(atlas), "q1",
                str(root / "empty-project"), str(hits), identity_mapping_arg=str(mapping))
            self.assertEqual(0, empty["inside_count"])

    def test_tampered_mapping_fails_before_member_artifacts_are_written(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            frozen, freeze, atlas, hits, mapping = self._inputs(root)
            projected = direct_rescue_project(str(frozen), str(freeze), str(atlas), "q1",
                str(root / "project"), str(hits), identity_mapping_arg=str(mapping))
            (root / "uniprot.json").write_text("{}", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "identity mapping source hash mismatch"):
                build_route_c_direct_members(str(frozen), str(freeze), projected["projection"],
                    projected["manifest"], str(atlas), str(root / "members"), case_id="case1")
            self.assertFalse((root / "members").exists())

    def test_wrong_provider_gene_species_protein_namespace_and_cooccurrence_rejected(self):
        mutations = [
            lambda p, r: p["uniProtKBCrossReferences"][0].update(id="AT3G99999"),
            lambda p, r: p["organism"].update(scientificName="Pisum sativum"),
            lambda p, r: p.update(primaryAccession="Q99999"),
            lambda p, r: r.update(source_namespace="PDB"),
            lambda p, r: r.update(source_identity_value="AT3G99999"),
            lambda p, r: (p.update(comment="AT3G27920 P27900 Arabidopsis thaliana"),
                           p["uniProtKBCrossReferences"][0].update(id="AT3G99999")),
        ]
        for index, mutate in enumerate(mutations):
            with self.subTest(index=index), tempfile.TemporaryDirectory() as temp:
                with self.assertRaises(ValueError):
                    self._project(Path(temp), mutate=mutate)

    def test_mapping_and_source_tampering_rejected_by_final_consumer(self):
        for artifact in ("identity.tsv", "uniprot.json"):
            with self.subTest(artifact=artifact), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                _, members, atlas = self._project(root)
                with (root / artifact).open("a", encoding="utf-8") as stream:
                    stream.write("\n")
                with self.assertRaisesRegex(ValueError, "identity mapping.*hash mismatch"):
                    self._read(members, atlas)

    def test_changed_binding_identity_still_cannot_authorize_sequence(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self._inputs(root)
            binding = pd.read_csv(root / "binding.tsv", sep="\t", dtype=str).fillna("")
            binding["reported_identity_value"] = "AT3G99999"
            binding.to_csv(root / "binding.tsv", sep="\t", index=False)
            freeze = json.loads((root / "freeze.json").read_text())
            freeze["claim_binding"]["sha256"] = sha256(root / "binding.tsv")
            (root / "freeze.json").write_text(json.dumps(freeze), encoding="utf-8")
            projected = direct_rescue_project(str(root / "frozen.tsv"), str(root / "freeze.json"),
                str(root / "atlas.tsv"), "q1", str(root / "project"), str(root / "hits.tsv"),
                identity_mapping_arg=str(root / "identity.tsv"))
            members = build_route_c_direct_members(str(root / "frozen.tsv"), str(root / "freeze.json"),
                projected["projection"], projected["manifest"], str(root / "atlas.tsv"),
                str(root / "members"), case_id="case1")
            self.assertTrue(self._read(members, root / "atlas.tsv").empty)
