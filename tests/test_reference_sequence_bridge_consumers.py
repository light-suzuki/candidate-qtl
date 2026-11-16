"""Official-record sequence correspondence reaches sealed C review and Wet.

These synthetic fixtures exercise source identity, projection and final consumers;
they are not a biological benchmark or evidence that BLAST proves a phenotype.
"""
from __future__ import annotations

import json
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from candidate_qtl.core import sha256
from candidate_qtl.observation_tiers import CHECKLIST_FIELDS
from candidate_qtl.paper_reports import PAPER_FIELD, PAPER_TABLE_POLICY, normalize_paper_reports
from candidate_qtl.relevance_tiers import classify_relevance
from candidate_qtl.route_c_direct_members import build_route_c_direct_members
from candidate_qtl.route_c_direct_projection import direct_rescue_project
from candidate_qtl.source_verify import source_verify
from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from tests._issue664_helpers import manifest_for
from tests import test_issue722_c_review_retention as fixtures


_CODONS = dict(zip("ACDEFGHIKLMNPQRSTVWY", "GCT TGT GAT GAA TTT GGT CAT ATT AAA CTG ATG AAT CCT CAA CGT TCT ACT GTT TGG TAT".split()))
_LOCAL_BLAST = bool(shutil.which("blastp") and shutil.which("makeblastdb"))


@unittest.skipUnless(_LOCAL_BLAST, "local BLAST+ required for sequence bridge consumer integration")
class ReferenceSequenceBridgeConsumerTests(unittest.TestCase):
    def _inputs(self, root: Path, *, copies: bool = False, partial_model: bool = False,
                exact_target: bool = False, complete_binding: bool = False,
                records_available: bool = True, no_hit: bool = False,
                request_gene: str = "GENE_A", outside_partial_model: bool = False,
                accession_endpoint: bool = False):
        """Build real frozen C inputs and an INSDC-shaped, local-only record."""
        from candidate_qtl.reference_sequence_bridge import build_reference_sequence_bridge

        generator = random.Random(940)
        sequence = "M" + "".join(generator.choices("ACDEFGHIKLMNPQRSTVWY", k=159))
        nucleotide = "".join(_CODONS[amino_acid] for amino_acid in sequence)
        record = SeqRecord(Seq(nucleotide), id="AA123456.1", name="AA123456",
                           description="Synthetic gene A sequence correspondence fixture")
        record.annotations.update(molecule_type="DNA", organism="Pisum sativum",
                                  accessions=["AA123456"], sequence_version=1)
        record.features = [
            SeqFeature(SimpleLocation(0, len(nucleotide)), type="source",
                       qualifiers={"organism": ["Pisum sativum"], "db_xref": ["taxon:3888"]}),
            SeqFeature(SimpleLocation(0, len(nucleotide)), type="gene", qualifiers={"gene": ["GENE_A"]}),
            SeqFeature(SimpleLocation(0, len(nucleotide)), type="CDS",
                       qualifiers={"gene": ["GENE_A"], "protein_id": ["ABC12345.1"], "translation": [sequence]}),
        ]
        cache = root / "records"
        cache.mkdir()
        if records_available:
            SeqIO.write(record, cache / "AA123456.1.gb", "genbank")
        requests = root / "requests.tsv"
        pd.DataFrame([{
            "request_id": "r1", "source_id": "PMID:1", "reported_gene": request_gene,
            "reference_accession": "AA123456.1", "reference_taxid": "3888",
            "focal_taxid": "3888", "cds_protein_id": "ABC12345.1",
        }]).to_csv(requests, sep="\t", index=False)
        proteins = root / "proteome.fa"
        model = sequence
        if no_hit:
            model = "".join(generator.choices("ACDEFGHIKLMNPQRSTVWY", k=160))
        if partial_model:
            model += "".join(generator.choices("ACDEFGHIKLMNPQRSTVWY", k=160))
        protein_rows = [("model_a", model, "GENE_A" if exact_target else "native_gene_a")]
        if copies:
            protein_rows.append(("scaffold_model_a", sequence, "unplaced_scaffold_gene_a"))
        if outside_partial_model:
            outside_sequence = sequence + "".join(generator.choices("ACDEFGHIKLMNPQRSTVWY", k=160))
            protein_rows.append(("outside_merged_model", outside_sequence, "outside_merged_gene"))
        proteins.write_text("".join(f">{protein}\n{amino_acids}*\n" for protein, amino_acids, _ in protein_rows), encoding="utf-8")
        mapping = root / "protein_gene_map.tsv"
        pd.DataFrame([{"protein_id": protein, "gene_id": gene} for protein, _, gene in protein_rows]).to_csv(mapping, sep="\t", index=False)
        bridge = build_reference_sequence_bridge(str(requests), str(proteins), str(mapping), str(root / "bridge"), records_arg=str(cache), runner="native")

        fixture = fixtures.Fixture()
        frozen_record = fixture._record()
        frozen_record.update(reported_gene_id="GENE_A", reported_gene_symbol="GENE_A", reference_sequence_id="ABC12345.1")
        if accession_endpoint:
            frozen_record.update(reported_gene_id="", reported_gene_symbol=request_gene,
                                 reference_sequence_id="AA123456.1")
            original_binding = fixture._binding_row

            def accession_binding(record):
                row = original_binding(record)
                row.update(reported_identity_type="gene_symbol", reported_identity_value=request_gene)
                return row

            fixture._binding_row = accession_binding
        if complete_binding:
            frozen_record["evidence_directness"] = "D0_causal_perturbation_or_rescue"
            original_binding = fixture._binding_row

            def complete_row(record):
                row = original_binding(record)
                row.update(binding_status="verified_complete", phenotype_match_state="target_trait",
                           perturbation_to_phenotype_binding_state="verified_supports")
                return row

            fixture._binding_row = complete_row
        fixture._record = lambda: frozen_record
        frozen, freeze = fixture._frozen(root, include_binding=True)
        if complete_binding:
            freeze_manifest = json.loads(freeze.read_text(encoding="utf-8"))
            freeze_manifest["claim_binding"]["complete"] = True
            freeze.write_text(json.dumps(freeze_manifest), encoding="utf-8")
        atlas = root / "atlas.tsv"
        pd.DataFrame([{"gene_id": gene} for _, _, gene in protein_rows if gene != "outside_merged_gene"] +
                     [{"gene_id": f"G{i}"} for i in range(20)]).to_csv(atlas, sep="\t", index=False)
        return bridge, frozen, freeze, atlas

    def _project(self, root: Path, *, cli: bool = False, **options):
        bridge, frozen, freeze, atlas = self._inputs(root, **options)
        if cli:
            command = [sys.executable, "-m", "candidate_qtl", "route-c-direct-project",
                       "--frozen", str(frozen), "--freeze-manifest", str(freeze),
                       "--atlas", str(atlas), "--qtl-id", "q1", "--out-dir", str(root / "project"),
                       "--reference-bridge", str(bridge["manifest"]), "--json"]
            output = subprocess.run(command, text=True, capture_output=True, check=True)
            projected = json.loads(output.stdout)
        else:
            projected = direct_rescue_project(str(frozen), str(freeze), str(atlas), "q1", str(root / "project"),
                                             reference_bridge_arg=str(bridge["manifest"]))
        return bridge, projected, frozen, freeze, atlas

    def _members(self, root, projected, frozen, freeze, atlas):
        return build_route_c_direct_members(str(frozen), str(freeze), projected["projection"],
                                            projected["manifest"], str(atlas), str(root / "members"), case_id="case1")

    def _wet(self, root, members, atlas):
        candidates = root / "b.tsv"
        pd.DataFrame([{"candidate_gene_id": f"G{i}", "route": "b", "support_status": "verified_semantic_support"}
                      for i in range(20)]).to_csv(candidates, sep="\t", index=False)
        result = wet_candidate_handoff(str(root / "wet"), case_id="case1",
            route_b_candidates_arg=str(candidates), route_b_candidates_manifest_arg=str(manifest_for(candidates, route="b", qtl_id="q1")),
            route_c_members_arg=members["members"], route_c_manifest_arg=members["manifest"], interval_atlas_arg=str(atlas))
        return pd.read_csv(result["wet_candidates"], sep="\t", dtype=str, keep_default_na=False)

    def test_cli_projection_real_members_and_wet_keep_source_bound_review(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, projected, frozen, freeze, atlas = self._project(root, cli=True)
            original_hashes = (sha256(frozen), sha256(freeze))
            members = self._members(root, projected, frozen, freeze, atlas)
            self.assertEqual(0, members["supported_count"])
            wet = self._wet(root, members, atlas)
            row = wet.set_index("candidate_gene_id").loc["native_gene_a"]
            self.assertEqual("c", row.review_routes)
            self.assertEqual("False", row.route_c_supported)
            self.assertEqual("route_review_candidate", row.membership_basis)
            self.assertEqual(original_hashes, (sha256(frozen), sha256(freeze)))

    def test_scaffold_competitor_and_partial_annotation_are_not_dropped(self):
        for options, expected in (({"copies": True}, {"native_gene_a", "unplaced_scaffold_gene_a"}),
                                  ({"partial_model": True}, {"native_gene_a"})):
            with self.subTest(options=options), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                _, projected, frozen, freeze, atlas = self._project(root, **options)
                members = self._members(root, projected, frozen, freeze, atlas)
                self.assertEqual(0, members["supported_count"])
                wet = self._wet(root, members, atlas)
                self.assertTrue(expected.issubset(set(wet.candidate_gene_id)))
                self.assertTrue(wet[wet.candidate_gene_id.isin(expected)].route_c_supported.eq("False").all())

    def test_unrelated_source_or_sequence_endpoint_cannot_reuse_bridge(self):
        for field, value in (("source_id", "PMID:2"), ("reference_sequence_id", "ABC99999.1"),
                             ("reported_gene_id", "UNRELATED_GENE")):
            with self.subTest(field=field), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                bridge, frozen, freeze, atlas = self._inputs(root)
                frame = pd.read_csv(frozen, sep="\t", dtype=str, keep_default_na=False)
                frame[field] = value
                frame.to_csv(frozen, sep="\t", index=False)
                manifest = json.loads(freeze.read_text(encoding="utf-8"))
                manifest["frozen_sha256"] = sha256(frozen)
                freeze.write_text(json.dumps(manifest), encoding="utf-8")
                with self.assertRaises(ValueError):
                    direct_rescue_project(str(frozen), str(freeze), str(atlas), "q1", str(root / "project"),
                                         reference_bridge_arg=str(bridge["manifest"]))
                self.assertFalse((root / "project").exists())

    def test_frozen_accession_transports_different_researcher_namespace_as_review_only(self):
        for name in ("SpeciesGeneA", "paper:gene-alpha"):
            with self.subTest(name=name), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                _, projected, frozen, freeze, atlas = self._project(
                    root, request_gene=name, accession_endpoint=True, cli=True)
                members = self._members(root, projected, frozen, freeze, atlas)
                self.assertEqual(0, members["supported_count"])
                self.assertEqual(1, members["review_count"])
                ledger = pd.read_csv(members["members"], sep="\t", dtype=str, keep_default_na=False)
                self.assertEqual({"native_gene_a"}, set(ledger.candidate_gene_id))
                evidence = pd.read_csv(projected["projection_evidence"], sep="\t", dtype=str, keep_default_na=False)
                blast = evidence[evidence.mapping_method == "blastp"]
                self.assertEqual({"reference_sequence_id"}, set(blast.source_identity_field))
                self.assertEqual({"AA123456.1"}, set(blast.source_identity_value))
                self.assertEqual({"sequence_correspondence_review_only"}, set(blast.source_identity_verification_state))

    def test_accession_endpoint_cannot_borrow_a_different_request_label(self):
        from candidate_qtl.reference_bridge_consumers import direct_bridge_hits
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bridge, frozen, _, _ = self._inputs(root, request_gene="paper:alpha", accession_endpoint=True)
            records = pd.read_csv(frozen, sep="\t", dtype=str, keep_default_na=False).to_dict("records")
            records[0]["reported_gene_symbol"] = "paper:other"
            with self.assertRaisesRegex(ValueError, "structured frozen identity"):
                direct_bridge_hits(bridge["manifest"], {r["record_id"]: r for r in records})

    def test_partial_annotation_cannot_gain_support_through_simultaneous_exact_id(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, projected, frozen, freeze, atlas = self._project(root, partial_model=True,
                                                               exact_target=True, complete_binding=True)
            members = self._members(root, projected, frozen, freeze, atlas)
            self.assertEqual(0, members["supported_count"])
            wet = self._wet(root, members, atlas)
            row = wet.set_index("candidate_gene_id").loc["GENE_A"]
            self.assertEqual("False", row.route_c_supported)
            self.assertEqual("route_review_candidate", row.membership_basis)

    def test_full_correspondence_unavailable_record_and_nohit_do_not_veto_independent_exact_support(self):
        for options, expected_blast_state in (({}, "evaluated_with_hit"),
                ({"records_available": False}, "not_evaluated_reference_bridge"),
                ({"no_hit": True}, "evaluated_no_hit")):
            with self.subTest(options=options), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                _, projected, frozen, freeze, atlas = self._project(root, exact_target=True,
                                                                  complete_binding=True, **options)
                evidence = pd.read_csv(projected["projection_evidence"], sep="\t", dtype=str, keep_default_na=False)
                self.assertIn(expected_blast_state, set(evidence[evidence.mapping_method.eq("blastp")].method_evaluation_state))
                members = self._members(root, projected, frozen, freeze, atlas)
                self.assertEqual(1, members["supported_count"])
                wet = self._wet(root, members, atlas)
                self.assertEqual("True", wet.set_index("candidate_gene_id").loc["GENE_A"].route_c_supported)

    def test_outside_partial_model_does_not_demote_inside_full_model_exact_support(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bridge, projected, frozen, freeze, atlas = self._project(root, exact_target=True,
                complete_binding=True, outside_partial_model=True)
            bridge_manifest = json.loads(Path(bridge["manifest"]).read_text(encoding="utf-8"))
            correspondence = pd.read_csv(bridge_manifest["outputs"]["correspondence"], sep="\t", dtype=str, keep_default_na=False)
            outside = correspondence[correspondence.target_id.eq("outside_merged_gene")]
            self.assertEqual({"partial_model_correspondence"}, set(outside.whole_annotation_model_identity))
            members = self._members(root, projected, frozen, freeze, atlas)
            self.assertEqual(1, members["supported_count"])
            wet = self._wet(root, members, atlas)
            self.assertEqual("True", wet.set_index("candidate_gene_id").loc["GENE_A"].route_c_supported)
            self.assertNotIn("outside_merged_gene", set(wet.candidate_gene_id))

    def test_bridge_ancestor_drift_is_rejected_by_final_wet_consumer(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, projected, frozen, freeze, atlas = self._project(root)
            members = self._members(root, projected, frozen, freeze, atlas)
            with (root / "records/AA123456.1.gb").open("a", encoding="utf-8") as handle:
                handle.write("\n")
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                self._wet(root, members, atlas)

    def _paper_table(self, root, *, id_column="hypothesis_id", source_column="source_ids",
                     source="PMID:1", gene="GENE_A", species="Pisum sativum"):
        checks = {**dict.fromkeys(CHECKLIST_FIELDS, "no"), "identity_resolved": "yes",
                  "phenotype_reported": "yes", "focal_species": "yes", "same_gene_as_focal": "unknown",
                  "equivalent_target_phenotype": "yes", "similar_cellular_or_developmental_phenotype": "yes"}
        report = {"paper_table_policy": PAPER_TABLE_POLICY, "source_id": source,
                  "source_title": "Synthetic gene phenotype experiment", "publication_type": "original",
                  "reported_gene": gene, "reported_gene_id": gene, "reported_species": species,
                  "reference_gene_id": "REF1", "gene_correspondence": "homolog",
                  "source_excerpt": "Mutation changes the reported developmental phenotype.", "source_location": "Results",
                  "initial_tier_reason": "Source phenotype is known; focal identity is unresolved.",
                  "unknown_reason": "Sequence correspondence has not been assessed.",
                  "followup_state": "completed", "paper_checklist": checks}
        row = {id_column: "h1", source_column: source, "reference_gene_name": "REF1",
               "reference_species": "Pisum sativum", PAPER_FIELD: normalize_paper_reports([report])}
        table = root / "papers.tsv"
        pd.DataFrame([row]).to_csv(table, sep="\t", index=False)
        return table, row

    def test_a_b_source_verify_attaches_postfreeze_view_without_rewriting_paper_or_tier(self):
        for id_column, source_column in (("hypothesis_id", "source_ids"), ("claim_id", "source_id")):
            for policy in ("initial_bibliography_v1", "content_v1"):
                with self.subTest(id_column=id_column, policy=policy), tempfile.TemporaryDirectory() as temporary:
                    root = Path(temporary)
                    bridge, *_ = self._inputs(root)
                    table, original = self._paper_table(root, id_column=id_column, source_column=source_column)
                    original_sha, original_tier = sha256(table), classify_relevance(original)[0]
                    result = source_verify(str(table), str(root / "verify"), id_column=id_column,
                        source_column=source_column, verification_policy=policy, reference_bridge_arg=str(bridge["manifest"]))
                    descriptor = result["paper_sequence_correspondence"]
                    view = pd.read_csv(descriptor["path"], sep="\t", dtype=str, keep_default_na=False)
                    self.assertEqual("post_freeze_parent_only", descriptor["audience"])
                    self.assertEqual("none", descriptor["membership_effect"])
                    self.assertEqual({"unknown"}, set(view.original_same_gene_as_focal))
                    self.assertEqual({"post_freeze_parent_review_required"}, set(view.assessment_state))
                    self.assertEqual(original_sha, sha256(table))
                    saved = pd.read_csv(table, sep="\t", dtype=str, keep_default_na=False).to_dict("records")[0]
                    self.assertEqual(original[PAPER_FIELD], saved[PAPER_FIELD])
                    self.assertEqual(original_tier, classify_relevance(saved)[0])

    def test_a_b_sidecar_rejects_a_different_paper_or_reported_gene(self):
        for options in ({"source": "PMID:2"}, {"gene": "UNRELATED_GENE"}):
            with self.subTest(options=options), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                bridge, *_ = self._inputs(root)
                table, _ = self._paper_table(root, **options)
                original_sha = sha256(table)
                with self.assertRaisesRegex(ValueError, "matching original paper/gene"):
                    source_verify(str(table), str(root / "verify"), verification_policy="initial_bibliography_v1",
                                  reference_bridge_arg=str(bridge["manifest"]))
                self.assertEqual(original_sha, sha256(table))
                self.assertFalse((root / "verify/paper_sequence_correspondence.tsv").exists())

    def test_requested_name_cannot_claim_a_different_official_cds_as_same_gene(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bridge, *_ = self._inputs(root, request_gene="GENE_B")
            table, _ = self._paper_table(root, gene="GENE_B")
            result = source_verify(str(table), str(root / "verify"), verification_policy="initial_bibliography_v1",
                                   reference_bridge_arg=str(bridge["manifest"]))
            view = pd.read_csv(result["paper_sequence_correspondence"]["path"], sep="\t", dtype=str, keep_default_na=False)
            self.assertEqual({"unresolved_reported_identity_to_sequence"}, set(view.paper_to_sequence_binding_state))
            self.assertTrue(all("GENE_A" in json.loads(value) for value in view.source_gene_names_json))
            self.assertEqual({"unknown"}, set(view.original_same_gene_as_focal))
            self.assertEqual({"none"}, set(view.membership_effect))

    def test_same_symbol_in_another_paper_species_is_not_a_positive_same_gene_binding(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bridge, *_ = self._inputs(root)
            table, _ = self._paper_table(root, species="Arabidopsis thaliana")
            result = source_verify(str(table), str(root / "verify"), verification_policy="initial_bibliography_v1",
                                   reference_bridge_arg=str(bridge["manifest"]))
            view = pd.read_csv(result["paper_sequence_correspondence"]["path"], sep="\t", dtype=str, keep_default_na=False)
            self.assertEqual({"structured_identity_match_taxon_requires_assessment"}, set(view.paper_to_sequence_binding_state))
            self.assertEqual({"Arabidopsis thaliana"}, set(view.paper_reported_species))
            self.assertEqual({"3888"}, set(view.reference_taxid))
            self.assertEqual({"unknown"}, set(view.original_same_gene_as_focal))


if __name__ == "__main__":
    unittest.main()
