"""Neutral frozen-paper/official-sequence edges reach review Wet, never support."""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from candidate_qtl.core import sha256
from candidate_qtl.route_c_claim_binding import canonical_frozen_record_sha256
from candidate_qtl.route_c_claim_binding_producer import prepare_claim_binding_review, seal_claim_binding_review
from candidate_qtl.route_c_direct_pipeline import DIRECT_PROPOSAL_COLUMNS, direct_rescue_freeze
from candidate_qtl.route_c_direct_projection import direct_rescue_project
from candidate_qtl.route_c_identity_mapping import load_identity_mapping
from candidate_qtl.route_c_paper_identity_mapping import prepare_paper_mapping_review
from candidate_qtl.route_c_review_candidates import read_direct_c_review_candidates
from tests import test_reference_sequence_bridge_consumers as bridge_fixtures


@unittest.skipUnless(shutil.which("blastp") and shutil.which("makeblastdb"), "local BLAST+ required")
class PaperIdentityMappingTests(unittest.TestCase):
    def replace_official_bridge(self, root, bridge, *, organism, taxid):
        from candidate_qtl.reference_sequence_bridge import build_reference_sequence_bridge
        source = root / "records" / "AA123456.1.gb"
        record = SeqIO.read(source, "genbank")
        record.annotations["organism"] = organism
        for feature in record.features:
            if feature.type == "source":
                feature.qualifiers.update(organism=[organism], db_xref=[f"taxon:{taxid}"])
        SeqIO.write(record, source, "genbank")
        requests = pd.read_csv(root / "requests.tsv", sep="\t", dtype=str)
        requests.loc[:, ["reference_taxid", "focal_taxid"]] = taxid
        request_path = root / "replacement-requests.tsv"
        requests.to_csv(request_path, sep="\t", index=False)
        bridge.update(build_reference_sequence_bridge(str(request_path), str(root / "proteome.fa"),
            str(root / "protein_gene_map.tsv"), str(root / "replacement-bridge"),
            records_arg=str(root / "records"), runner="native"))

    def add_taxonomy(self, root, edge, *, taxid="3888", scientific="Official synthetic pea", synonyms=None):
        path = root / f"taxonomy-{taxid}.json"
        path.write_text(json.dumps({"taxonId": int(taxid), "scientificName": scientific,
            "synonyms": ["Pisum sativum"] if synonyms is None else synonyms}), encoding="utf-8")
        edge.loc[:, "taxonomy_record_path"] = str(path)
        edge.loc[:, "taxonomy_record_sha256"] = sha256(path)
        return path

    def inputs(self, root, *, blank=False, copies=False, name="PaperAlpha", mutate=None,
               records_available=True, complete=False, no_hit=False):
        helper = bridge_fixtures.ReferenceSequenceBridgeConsumerTests()
        bridge, frozen, _, atlas = helper._inputs(root, request_gene=name, accession_endpoint=True,
            copies=copies, records_available=records_available, no_hit=no_hit)
        row = pd.read_csv(frozen, sep="\t", dtype=str, keep_default_na=False).to_dict("records")[0]
        if blank:
            row["reference_sequence_id"] = ""
        if complete:
            row["evidence_directness"] = "D0_causal_perturbation_or_rescue"
        proposal = root / "paper-proposal.tsv"
        pd.DataFrame([row], columns=DIRECT_PROPOSAL_COLUMNS).to_csv(proposal, sep="\t", index=False)
        quote = f"The {name} gene sequence was deposited as AA123456.1. The {name} mutant has fasciated stem in Pisum sativum."
        cache, extraction = root / "paper-cache.json", root / "extraction.json"
        cache.write_text(json.dumps([{"id": "1", "pmid": "1", "title": f"{name} stem development",
                                     "abstractText": quote, "pubYear": "2020"}]), encoding="utf-8")
        extraction.write_text(json.dumps([{"hypothesis_id": "r1", "citation": "PMID:1",
            "evidence_relation": "direct_genetic_test", "quoted_support": quote}]), encoding="utf-8")
        freeze = direct_rescue_freeze(str(proposal), str(root / "real-freeze"), str(cache), claim_extraction_arg=str(extraction))
        prepared = prepare_claim_binding_review(freeze["frozen"], freeze["manifest"], str(root / "claim-prepare"))
        review = pd.read_csv(prepared["review"], sep="\t", dtype=str, keep_default_na=False)
        review.loc[:, "source_organism_state"] = "verified_match"
        review.loc[:, "source_accession_or_cultivar_state"] = "verified_match"
        review.loc[:, "reported_identity_verification_state"] = "verified_exact"
        review.loc[:, "identity_mapping_basis"] = "exact_source_text"
        review.loc[:, "reported_variant_state"] = "not_applicable"
        review.loc[:, "phenotype_match_state"] = "bounded_related_trait"
        review.loc[:, "experimental_target_gene_state"] = "verified_same_identity"
        review.loc[:, "curation_method"] = "independent_agent_review"
        review.loc[:, "reviewed_by"] = "neutral-fixture-reviewer"
        review.loc[:, "binding_reason"] = "Synthetic source identity reviewed; no focal causal support asserted."
        if complete:
            review.loc[:, "phenotype_match_state"] = "target_trait"
            review.loc[:, "perturbation_to_phenotype_binding_state"] = "verified_supports"
        reviewed = root / "claim-reviewed.tsv"
        review.to_csv(reviewed, sep="\t", index=False)
        bound = seal_claim_binding_review(freeze["frozen"], freeze["manifest"], str(reviewed), str(root / "claim-bound"))
        template = prepare_paper_mapping_review(freeze["frozen"], bound["freeze_manifest"], str(root / "edge-prepare"))
        edge = pd.read_csv(template["review"], sep="\t", dtype=str, keep_default_na=False)
        edge.loc[:, "target_identity_value"] = "AA123456.1"
        edge.loc[:, "citation"] = "PMID:1"
        edge.loc[:, "source_scope"] = "abstract"
        edge.loc[:, "locator"] = "abstract, sentences 1-2"
        edge.loc[:, "quoted_mapping"] = quote
        edge.loc[:, "curation_method"] = "independent_agent_review"
        edge.loc[:, "reviewed_by"] = "neutral-fixture-reviewer"
        edge.loc[:, "mapping_reason"] = "Synthetic reviewer confirms this gene's deposited sequence relation, not cooccurrence."
        if mutate:
            mutate(edge, freeze, bound, bridge)
        mapping = root / "paper-edge.tsv"
        edge.to_csv(mapping, sep="\t", index=False)
        return helper, bridge, freeze["frozen"], bound["freeze_manifest"], atlas, mapping

    def project(self, root, *, cli=False, **kwargs):
        helper, bridge, frozen, freeze, atlas, mapping = self.inputs(root, **kwargs)
        if cli:
            command = [sys.executable, "-m", "candidate_qtl", "route-c-direct-project", "--frozen", str(frozen),
                "--freeze-manifest", str(freeze), "--atlas", str(atlas), "--qtl-id", "q1", "--out-dir", str(root / "project"),
                "--reference-bridge", bridge["manifest"], "--identity-mapping", str(mapping), "--json"]
            execution = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(0, execution.returncode, execution.stderr)
            project = json.loads(execution.stdout)
        else:
            project = direct_rescue_project(frozen, freeze, str(atlas), "q1", str(root / "project"),
                identity_mapping_arg=str(mapping), reference_bridge_arg=bridge["manifest"])
        members = helper._members(root, project, frozen, freeze, atlas)
        return helper, project, members, atlas, frozen, freeze, bridge, mapping

    def test_real_blast_cli_members_and_wet_keep_review_copies_and_blank_original(self):
        for blank, name in ((False, "PaperAlpha"), (True, "paper:alpha")):
            with self.subTest(blank=blank), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                helper, project, members, atlas, frozen, freeze, _, _ = self.project(
                    root, cli=True, blank=blank, copies=True, name=name)
                hashes = sha256(Path(frozen)), sha256(Path(freeze))
                self.assertEqual(0, members["supported_count"])
                self.assertEqual(2, members["review_count"])
                wet = helper._wet(root, members, atlas)
                expected = {"native_gene_a", "unplaced_scaffold_gene_a"}
                self.assertTrue(expected <= set(wet.candidate_gene_id))
                self.assertTrue(wet[wet.candidate_gene_id.isin(expected)].route_c_supported.eq("False").all())
                self.assertEqual(hashes, (sha256(Path(frozen)), sha256(Path(freeze))))
                original = pd.read_csv(frozen, sep="\t", dtype=str, keep_default_na=False).iloc[0]
                self.assertEqual("" if blank else "AA123456.1", original.reference_sequence_id)
                evidence = pd.read_csv(project["projection_evidence"], sep="\t", dtype=str)
                self.assertTrue(evidence[evidence.mapping_method.eq("blastp")].source_identity_verification_state.eq(
                    "reviewed_paper_sequence_mapping").all())

    def test_untrusted_wrong_or_incomplete_edges_fail_before_projection(self):
        changes = {
            "anonymous": ("reviewed_by", ""), "unattested": ("curation_method", ""),
            "unknown_reviewer": ("reviewed_by", "unknown"), "unknown_relation_reason": ("mapping_reason", "unknown"),
            "wrong_method": ("curation_method", "model_memory"), "unknown_relation": ("mapping_relation", "cooccurs"),
            "wrong_gene": ("source_identity_value", "OtherAlpha"), "wrong_accession": ("target_identity_value", "AA999999.1"),
            "wrong_version": ("target_identity_value", "AA123456.2"), "wrong_species": ("reference_species", "Other plant"),
            "wrong_paper": ("citation", "PMID:2"), "wrong_record": ("record_id", "other-record"),
            "title_only": ("source_scope", "title"), "unretrieved_full_text": ("source_scope", "full_text"),
            "unquoted_relation": ("quoted_mapping", "PaperAlpha AA123456.1"),
            "phenotype_only_quote": ("quoted_mapping", "The PaperAlpha mutant has fasciated stem in Pisum sativum."),
            "source_digest": ("source_record_sha256", "0" * 64), "support_rank": ("rank_contribution", "support"),
        }
        for name, (field, value) in changes.items():
            with self.subTest(name=name), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                def mutate(edge, *_):
                    edge.loc[:, field] = value
                with self.assertRaises((ValueError, subprocess.CalledProcessError)):
                    self.project(root, mutate=mutate)
                self.assertFalse((root / "project").exists())

    def test_arbitrary_source_file_with_matching_text_is_not_frozen_paper(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            def mutate(edge, *_):
                arbitrary = root / "invented.json"
                arbitrary.write_text(edge.iloc[0].quoted_mapping, encoding="utf-8")
                edge.loc[:, "source_record_path"] = str(arbitrary)
                edge.loc[:, "source_record_sha256"] = sha256(arbitrary)
            with self.assertRaisesRegex(ValueError, "frozen paper artifact"):
                self.project(root, mutate=mutate)

    def test_paper_edge_requires_official_bridge_even_with_valid_paper(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, _, frozen, freeze, _, mapping = self.inputs(root)
            records = {r["record_id"]: r for r in pd.read_csv(frozen, sep="\t", dtype=str, keep_default_na=False).to_dict("records")}
            with self.assertRaisesRegex(ValueError, "official sequence bridge"):
                load_identity_mapping(mapping, records, freeze_manifest_path=Path(freeze))

    def test_complete_phenotype_binding_still_does_not_promote_blast_support(self):
        with tempfile.TemporaryDirectory() as temporary:
            helper, _, members, atlas, *_ = self.project(Path(temporary), complete=True)
            self.assertEqual(0, members["supported_count"])
            self.assertIn("native_gene_a", set(helper._wet(Path(temporary), members, atlas).candidate_gene_id))

    def test_valid_paper_edge_without_interval_hit_stays_ledger_only(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            helper, project, members, atlas, *_ = self.project(root, no_hit=True, blank=True)
            self.assertEqual(0, project["inside_count"])
            self.assertEqual(0, members["review_count"])
            self.assertNotIn("native_gene_a", set(helper._wet(root, members, atlas).candidate_gene_id))

    def test_exact_accession_label_without_sequence_hit_cannot_admit_symbol_to_wet(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            def accession_label(edge, *_):
                atlas = pd.read_csv(root / "atlas.tsv", sep="\t", dtype=str)
                pd.concat([atlas, pd.DataFrame([{"gene_id": "AA123456.1"}])], ignore_index=True).to_csv(
                    root / "atlas.tsv", sep="\t", index=False)
            helper, project, members, atlas, *_ = self.project(root, no_hit=True, blank=True, mutate=accession_label)
            self.assertEqual(1, project["inside_count"])
            self.assertEqual(0, members["supported_count"])
            self.assertNotIn("AA123456.1", set(helper._wet(root, members, atlas).candidate_gene_id))

    def test_official_sequence_unavailable_does_not_authorize_mapping(self):
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(ValueError, "official sequence bridge"):
                self.project(Path(temporary), records_available=False)

    def test_wrong_taxid_and_conflicting_original_accession_cannot_be_enriched(self):
        from candidate_qtl.reference_sequence_bridge import build_reference_sequence_bridge
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            def wrong_taxid(edge, freeze, bound, bridge):
                request = pd.read_csv(root / "requests.tsv", sep="\t", dtype=str)
                request.loc[:, ["reference_taxid", "focal_taxid"]] = "9999"
                path = root / "wrong-taxid.tsv"
                request.to_csv(path, sep="\t", index=False)
                replacement = build_reference_sequence_bridge(str(path), str(root / "proteome.fa"),
                    str(root / "protein_gene_map.tsv"), str(root / "wrong-taxid-bridge"),
                    records_arg=str(root / "records"), runner="native")
                bridge.update(replacement)
            with self.assertRaisesRegex(ValueError, "official sequence bridge"):
                self.project(root, mutate=wrong_taxid)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, bridge, frozen, freeze, _, mapping = self.inputs(root)
            records = {r["record_id"]: r for r in pd.read_csv(frozen, sep="\t", dtype=str, keep_default_na=False).to_dict("records")}
            records["r1"]["reference_sequence_id"] = "AA999999.1"
            with self.assertRaisesRegex(ValueError, "frozen roster mismatch"):
                load_identity_mapping(mapping, records, freeze_manifest_path=Path(freeze), reference_bridge_arg=bridge["manifest"])

    def test_coherent_wrong_official_species_cannot_borrow_case_label(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            def mismatch(edge, freeze, bound, bridge):
                self.replace_official_bridge(root, bridge, organism="Arabidopsis thaliana", taxid="3702")
            with self.assertRaisesRegex(ValueError, "official organism differs"):
                self.project(root, mutate=mismatch)

    def test_bound_same_taxon_synonym_reaches_wet_without_rewriting_original(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            def synonym(edge, freeze, bound, bridge):
                self.replace_official_bridge(root, bridge, organism="Official synthetic pea", taxid="3888")
                self.add_taxonomy(root, edge)
            helper, _, members, atlas, frozen, *_ = self.project(root, blank=True, mutate=synonym)
            self.assertEqual(0, members["supported_count"])
            self.assertIn("native_gene_a", set(helper._wet(root, members, atlas).candidate_gene_id))
            raw = pd.read_csv(frozen, sep="\t", dtype=str, keep_default_na=False).iloc[0]
            self.assertEqual("Pisum sativum", raw.reference_species)
            self.assertEqual("Pisum sativum", raw.focal_species)
            self.assertEqual("", raw.reference_sequence_id)
            taxonomy = root / "taxonomy-3888.json"
            original = taxonomy.read_bytes()
            taxonomy.write_bytes(original + b"\n")
            with self.assertRaisesRegex(ValueError, "taxonomy hash mismatch"):
                read_direct_c_review_candidates(members["members"], members["manifest"], case_id="case1",
                    atlas_frame=pd.read_csv(atlas, sep="\t", dtype=str), supported_ids=set())

    def test_resealed_wrong_taxonomy_id_or_unrelated_names_do_not_authorize_synonym(self):
        scenarios = [("3702", "Official synthetic pea", ["Pisum sativum"]),
                     ("3888", "Unrelated plant", ["Different plant"]),
                     ("3888", "Official synthetic pea", ["Different plant"])]
        for taxid, scientific, synonyms in scenarios:
            with self.subTest(taxid=taxid, scientific=scientific), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                def invalid(edge, freeze, bound, bridge):
                    self.replace_official_bridge(root, bridge, organism="Official synthetic pea", taxid="3888")
                    self.add_taxonomy(root, edge, taxid=taxid, scientific=scientific, synonyms=synonyms)
                with self.assertRaisesRegex(ValueError, "taxonomy does not bind"):
                    self.project(root, mutate=invalid)

    def test_final_consumer_rechecks_every_mapping_ancestor(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, _, members, atlas, frozen, freeze, bridge, mapping = self.project(root, blank=True)
            bound = json.loads(Path(freeze).read_text())
            captured = pd.read_csv(mapping, sep="\t", dtype=str).iloc[0].source_record_path
            artifacts = [Path(mapping), Path(captured), Path(frozen), Path(freeze), Path(bridge["manifest"]),
                Path(bound["claim_binding"]["artifact"]), Path(bound["claim_binding_producer"]["source_validation"]),
                Path(bridge["manifest"]).parent / "raw" / "record_000001.txt"]
            for artifact in artifacts:
                with self.subTest(artifact=artifact.name):
                    original = artifact.read_bytes()
                    try:
                        artifact.write_bytes(original + b"\n")
                        with self.assertRaises(ValueError):
                            read_direct_c_review_candidates(members["members"], members["manifest"], case_id="case1",
                                atlas_frame=pd.read_csv(atlas, sep="\t", dtype=str), supported_ids=set())
                    finally:
                        artifact.write_bytes(original)

    def test_valid_edge_cannot_repair_unverified_source_context_variant_or_origin(self):
        scenarios = [("source_verification_state", "pending"), ("trait_relation_state", "context_only"),
                     ("evidence_directness", "D3_context_only"), ("focal_species", "Heterologous host")]
        for field, value in scenarios:
            with self.subTest(field=field), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                def mutate(edge, freeze, bound, bridge):
                    path = Path(freeze["frozen"])
                    data = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
                    data.loc[:, field] = value
                    data.loc[:, "direct_rescue_eligibility"] = "pending_verification"
                    data.to_csv(path, sep="\t", index=False)
                    manifest_path = Path(bound["freeze_manifest"])
                    manifest = json.loads(manifest_path.read_text())
                    manifest["frozen_sha256"] = sha256(path)
                    binding_path = Path(manifest["claim_binding"]["artifact"])
                    binding = pd.read_csv(binding_path, sep="\t", dtype=str, keep_default_na=False)
                    binding.loc[:, "frozen_record_sha256"] = canonical_frozen_record_sha256(data.to_dict("records")[0])
                    binding.to_csv(binding_path, sep="\t", index=False)
                    manifest["claim_binding"]["sha256"] = sha256(binding_path)
                    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
                with self.assertRaises(ValueError):
                    self.project(root, blank=True, mutate=mutate)
        with tempfile.TemporaryDirectory() as temporary:
            def bad_variant(edge, freeze, bound, bridge):
                manifest_path = Path(bound["freeze_manifest"])
                manifest = json.loads(manifest_path.read_text())
                binding_path = Path(manifest["claim_binding"]["artifact"])
                binding = pd.read_csv(binding_path, sep="\t", dtype=str, keep_default_na=False)
                binding.loc[:, "reported_variant_state"] = "not_evaluated"
                binding.to_csv(binding_path, sep="\t", index=False)
                manifest["claim_binding"]["sha256"] = sha256(binding_path)
                manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "source-bound symbol identity"):
                self.project(Path(temporary), mutate=bad_variant)
