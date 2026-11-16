"""Original-paper evidence stays separate from provider retrieval bibliography."""
import copy
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.claim_binding import route_a_authorizes
from candidate_qtl.family_candidate_bridge import _route_a_review_projection_rows
from candidate_qtl.llm_pathway import hypothesis_verify
from tests.bound_claim_fixtures import bound_claim
from tests.test_paper_handoff_consumers import paper


def fixture(root: Path, *, policy="initial_bibliography_v1", original_source="PMID:1111111"):
    provenance = {"provider_linked_source_ids": ["PMID:2222222"],
                  "original_source_ids": [original_source] if original_source else [],
                  "record_identity_state": "resolved_exact_name_species",
                  "input_proposed_stable_id": "Q96524", "candidate_count": 1,
                  "validated_reference_accession": "Q96524",
                  "validated_reference_metadata": {"ids": ["uniprot:Q96524"], "names": ["REF1", "REF1_COPY"],
                                                   "species": ["Arabidopsis thaliana"]},
                  "copy_correspondence": "explicit original copy"}
    row = {"hypothesis_id": "h1", "case_id": "synthetic", "gene_family": "meristem regulator",
           "reference_gene_name": "REF1", "reference_species": "Arabidopsis thaliana", "proposed_stable_id": "Q96524",
           "expected_role": "controls meristem cell proliferation", "mechanistic_link": "altered meristem proliferation could change organ growth",
           "member_interval_genes": "", "candidate_admission_policy": policy, "bibliography_state": "resolved",
           "source_ids": original_source, "verified_source_ids": original_source,
           "source_title": "Original meristem phenotype study",
           "reference_record_source_ids": "PMID:2222222", "retrieval_provenance": json.dumps(provenance),
           "initial_tier": "T2", "initial_tier_reason": "Original paper reports altered meristem proliferation",
           "reported_phenotype": "altered meristem proliferation", "reported_process": "meristem proliferation",
           "phenotype_comparison": "same_cellular_process", "focal_species": "Pisum sativum"}
    row["paper_reports_json"] = json.dumps([paper(1111111, "T2", source_title=row["source_title"],
        reported_gene="REF1", reported_gene_id="Q96524", reported_species="Arabidopsis thaliana",
        reference_gene_id="Q96524", gene_correspondence="homolog")] if original_source else [])
    table = root / "hypotheses.tsv"
    pd.DataFrame([row]).to_csv(table, sep="\t", index=False)
    record = {"primaryAccession": "Q96524", "entryType": "UniProtKB reviewed (Swiss-Prot)",
              "genes": [{"geneName": {"value": "REF1"}, "synonyms": [{"value": "REF1_COPY"}]}],
              "organism": {"scientificName": "Arabidopsis thaliana"}, "sequence": {"value": "M" * 100},
              "references": [{"citation": {"citationCrossReferences": [{"database": "PubMed", "id": "2222222"}]}}]}
    cache = root / "records.json"
    cache.write_text(json.dumps([record]), encoding="utf-8")
    return table, cache, row


def canonical_cli(root: Path):
    """Use declared synthetic correspondence, not a biological BLAST result."""
    table, cache, row = fixture(root)
    original = table.read_bytes(), cache.read_bytes()
    def run(*args):
        subprocess.run([sys.executable, "-m", "candidate_qtl", *map(str, args)], check=True, capture_output=True, text=True)
    literature = root / "literature.json"
    literature.write_text(json.dumps({"id": "1111111", "pmid": "1111111", "title": row["source_title"]}))
    source = root / "source"
    run("source-verify", "--mechanisms", table, "--cache", literature, "--out-dir", source,
        "--verification-policy", "initial_bibliography_v1", "--focal-species", "Pisum sativum")
    sourced = source / "mechanisms_with_verified_sources.tsv"
    assert pd.read_csv(sourced, sep="\t").iloc[0].bibliography_state == "resolved"
    verify = root / "verify"
    run("hypothesis-verify", "--hypotheses", sourced, "--record-cache", cache, "--out-dir", verify, "--focal-species", "Pisum sativum")
    validation = pd.read_csv(verify / "proposal_validation.tsv", sep="\t", dtype=str).fillna("")
    assert set(validation.iloc[0].function_source_ids.split(";")) == {"PMID:1111111", "UniProt:Q96524"}
    assert set(validation.iloc[0].reference_record_source_ids.split(";")) == {"PMID:2222222", "UniProt:Q96524"}
    panel = root / "panel"
    run("bait-panel-freeze", "--hypotheses", sourced, "--validation", verify / "proposal_validation.tsv",
        "--sequences", verify / "verified_sequences.json", "--reviewable-sequences", verify / "reviewable_sequences.json",
        "--out-dir", panel)
    blast = root / "blast.tsv"
    pd.DataFrame([{"panel_id": "Q96524", "provider_status": "evaluable", "top_hit_target_id": "G1",
                   "top_hit_identity": "90", "top_hit_query_coverage": "1", "top_hit_target_coverage": "1",
                   "top_hit_evalue": "1e-30", "paralog_count": "0", "ambiguity_status": "unique"}]).to_csv(blast, sep="\t", index=False)
    genes = root / "genes.tsv"
    pd.DataFrame([{"gene_id": "G1", "chrom": "chr1", "start": "1", "end": "300", "protein_length": "100"}]).to_csv(genes, sep="\t", index=False)
    bridge = root / "bridge"
    run("family-candidate-bridge", "--route-a-hypotheses", sourced, "--route-a-projection-panel", panel / "projection_bait_panel.tsv",
        "--route-a-blast", blast, "--genes", genes, "--out-dir", bridge, "--case-id", "synthetic")
    review = pd.read_csv(bridge / "route_a_review_candidates.tsv", sep="\t", dtype=str).fillna("")
    assert set(review.candidate_gene_id) == {"G1"}
    assert "2222222" not in review.to_json()
    evidence = json.loads(review.iloc[0].relevance_evidence_json)
    active = next(item for item in evidence if item.get("evidence_origin") == "route_a_hypothesis")
    kept_provenance = json.loads(active["retrieval_provenance"])
    original_provenance = json.loads(row["retrieval_provenance"])
    assert kept_provenance == {key: value for key, value in original_provenance.items() if key != "provider_linked_source_ids"}
    assert active["expected_role"] == row["expected_role"] and active["initial_tier"] == "T2"
    assert original == (table.read_bytes(), cache.read_bytes())
    manifest = json.loads((bridge / "route_a_review_candidates_manifest.json").read_text())
    assert manifest["inputs"]["route_a_hypotheses"]["sha256"] == hashlib.sha256(sourced.read_bytes()).hexdigest()
    wet = root / "wet"
    run("wet-candidate-handoff", "--route-a-bridge", bridge / "family_candidate_bridge.tsv",
        "--route-a-manifest", bridge / "family_candidate_bridge_manifest.json",
        "--review-candidates", bridge / "route_a_review_candidates.tsv",
        "--review-candidates-manifest", bridge / "route_a_review_candidates_manifest.json",
        "--interval-atlas", genes, "--out-dir", wet, "--case-id", "synthetic")
    adopted = pd.read_csv(wet / "wet_candidates.tsv", sep="\t", dtype=str).fillna("")
    assert set(adopted.candidate_gene_id) == {"G1"} and "2222222" not in adopted.to_json()
    assert set(adopted.tier_source_id) == {"PMID:1111111"}
    return {"review_gene_ids": sorted(review.candidate_gene_id), "active_provider_citations": 0,
            "original_inputs_unchanged": True, "identity_copy_provenance_preserved": True}


class InitialPaperProvenanceTests(unittest.TestCase):
    def test_initial_original_paper_and_identity_reach_real_cli_bridge(self):
        with tempfile.TemporaryDirectory() as tmp:
            canonical_cli(Path(tmp))

    def test_initial_missing_original_source_is_not_filled_by_provider_bibliography(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            table, cache, _ = fixture(root, original_source="")
            result = hypothesis_verify(str(table), str(root / "verify"), "Pisum sativum", record_cache_arg=str(cache))
            row = pd.read_csv(result["validation"], sep="\t", dtype=str).fillna("").iloc[0]
            self.assertEqual("UniProt:Q96524", row.function_source_ids)
            self.assertEqual("True", row.function_source_backed)
            self.assertEqual({"PMID:2222222", "UniProt:Q96524"}, set(row.reference_record_source_ids.split(";")))
            self.assertEqual("verified", row.identity_sequence_verification_status)

    def test_legacy_provider_function_sources_are_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            table, cache, _ = fixture(root, policy="content_v1")
            result = hypothesis_verify(str(table), str(root / "verify"), "Pisum sativum", record_cache_arg=str(cache))
            row = pd.read_csv(result["validation"], sep="\t", dtype=str).fillna("").iloc[0]
            self.assertEqual({"PMID:1111111", "PMID:2222222", "UniProt:Q96524"}, set(row.function_source_ids.split(";")))

    def test_missing_citation_never_authorizes_a_bound_claim(self):
        claim = bound_claim(route="a", claim_id="h1", subject_gene="Q96524", citation="PMID:1111111", organism="Arabidopsis thaliana")
        broken = copy.deepcopy(claim)
        broken.pop("citation")
        accepted, state, _ = route_a_authorizes(broken, hypothesis_id="h1", proposed_accession="Q96524", verified_source_ids={"PMID:1111111"})
        self.assertFalse(accepted)
        self.assertEqual("not_evaluated", state)

    def test_initial_provider_url_shape_stays_audit_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            table, cache, _ = fixture(root)
            records = json.loads(cache.read_text())
            records[0]["references"][0]["citation"]["citationCrossReferences"][0]["id"] = "https://pubmed.ncbi.nlm.nih.gov/2222222/"
            cache.write_text(json.dumps(records))
            result = hypothesis_verify(str(table), str(root / "verify"), "Pisum sativum", record_cache_arg=str(cache))
            row = pd.read_csv(result["validation"], sep="\t", dtype=str).fillna("").iloc[0]
            self.assertNotIn("https://", row.function_source_ids)
            self.assertIn("https://pubmed.ncbi.nlm.nih.gov/2222222/", row.reference_record_source_ids)

    def test_initial_explicit_bound_citation_is_preserved_without_unrelated_provider_paper(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            table, cache, _ = fixture(root)
            records = json.loads(cache.read_text())
            records[0]["references"][0]["citation"]["citationCrossReferences"].append({"database": "PubMed", "id": "3333333"})
            cache.write_text(json.dumps(records))
            claims = root / "bound.tsv"
            pd.DataFrame([bound_claim(route="a", claim_id="h1", subject_gene="Q96524", citation="PMID:2222222",
                                      organism="Arabidopsis thaliana")]).to_csv(claims, sep="\t", index=False)
            result = hypothesis_verify(str(table), str(root / "verify"), "Pisum sativum", record_cache_arg=str(cache), bound_claims_arg=str(claims))
            row = pd.read_csv(result["validation"], sep="\t", dtype=str).fillna("").iloc[0]
            self.assertEqual("True", row.bound_claim_authorized)
            self.assertEqual({"PMID:1111111", "PMID:2222222", "UniProt:Q96524"}, set(row.function_source_ids.split(";")))
            self.assertIn("PMID:3333333", row.reference_record_source_ids)

    def test_initial_projection_keeps_identity_and_unknown_audit_without_changing_legacy(self):
        panel = pd.DataFrame([{"hypothesis_id": "h1", "bait_id": "Q96524", "projection_eligibility": "review_only",
                               "review_status": "review_required", "function_source_ids": "PMID:1111111;UniProt:Q96524"}])
        blast = pd.DataFrame([{"panel_id": "Q96524", "provider_status": "evaluable", "top_hit_target_id": "G1",
                               "top_hit_identity": "90", "top_hit_query_coverage": "1", "top_hit_target_coverage": "1",
                               "top_hit_evalue": "1e-30", "paralog_count": "0", "ambiguity_status": "unique"}])
        for policy, provenance in [("initial_bibliography_v1", "not-json"), ("initial_bibliography_v1", "[]"),
                                   ("content_v1", "not-json")]:
            with self.subTest(policy=policy, provenance=provenance), tempfile.TemporaryDirectory() as tmp:
                _, _, original = fixture(Path(tmp), policy=policy)
                original["retrieval_provenance"] = provenance
                snapshot = copy.deepcopy(original)
                rows = _route_a_review_projection_rows(panel, blast, {}, {"G1"}, "synthetic", hypotheses=[original])
                self.assertEqual({"G1"}, {row["candidate_gene_id"] for row in rows})
                evidence = json.loads(rows[0]["relevance_evidence_json"])
                active = next(item for item in evidence if item.get("evidence_origin") == "route_a_hypothesis")
                self.assertEqual("Q96524", active["proposed_stable_id"])
                if policy == "initial_bibliography_v1":
                    self.assertEqual("audit_only_unparsed", active["retrieval_provenance_state"])
                    self.assertNotIn("reference_record_source_ids", active)
                    self.assertNotIn("retrieval_provenance", active)
                else:
                    self.assertEqual(provenance, active["retrieval_provenance"])
                    self.assertEqual("PMID:2222222", active["reference_record_source_ids"])
                self.assertEqual(snapshot, original)
