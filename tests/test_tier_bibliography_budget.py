"""Offline work-budget and source-handoff checks, not a scientific benchmark."""
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.paper_reports import PAPER_FIELD, select_bibliography_papers
from candidate_qtl.route_b_semantic_bundle import refreeze_route_b_sources
from candidate_qtl.route_b_semantic_import import import_route_b_semantic_proposals
from candidate_qtl.route_b_semantic_source_verify import verify_route_b_semantic_sources
from candidate_qtl.route_b_semantic_transfer_verify import verify_route_b_semantic_transfer
from candidate_qtl.route_b_semantic_candidates import freeze_route_b_semantic_candidates
from scripts.agentic_route_b import prepare, merge
from tests.test_paper_checklist import report
from tests import test_route_b_semantic_import as import_fixtures


class TierBibliographyBudgetTests(unittest.TestCase):
    def test_native_cli_checks_five_strongest_papers_without_content_or_observation_loss(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            papers = [dict(report(annotation_supports_hypothesis="yes"), source_id=f"PMID:{i}", source_title=f"Paper {i}") for i in range(101, 107)]
            papers.append(dict(report(phenotype_reported="yes", equivalent_target_phenotype="yes", focal_species="yes", same_gene_as_focal="yes"),
                               source_id="DOI:10.1234/strong", source_identifier_aliases="PMID:107", source_title="Strong paper"))
            raw = json.dumps(papers)
            table = root / "a.tsv"
            pd.DataFrame([{"hypothesis_id": "h", "reference_gene": "REF1", "source_ids": ";".join([f"PMID:{i}" for i in range(101, 108)]), PAPER_FIELD: raw}]).to_csv(table, sep="\t", index=False)
            def metadata(kind, identifier):
                if identifier in {"107", "10.1234/strong"}:
                    return {"pmid": "107", "doi": "10.1234/strong", "title": "Strong paper"}
                return {"pmid": identifier, "title": f"Paper {identifier}"}
            with patch("candidate_qtl.bibliography.fetch_metadata", side_effect=metadata) as fetch, \
                 patch("candidate_qtl.gene_metadata.fetch_gene_metadata", return_value={"state": "not_found"}), \
                 patch("candidate_qtl.source_verify._fetch_europe_pmc", side_effect=AssertionError("must not fetch abstracts")):
                self.assertEqual(0, main(["source-verify", "--mechanisms", str(table), "--out-dir", str(root / "out"), "--allow-external-api"]))
            self.assertEqual(5, fetch.call_count)
            selected = json.loads((root / "out/bibliography_selection.json").read_text())["panels"][0]
            self.assertEqual(5, selected["selected_paper_count"])
            self.assertEqual(7, selected["candidate_paper_count"])
            self.assertIn("DOI:10.1234/strong", selected["selected_source_ids"])
            self.assertIn("PMID:107", selected["selected_source_ids"])
            output = pd.read_csv(root / "out/mechanisms_with_verified_sources.tsv", sep="\t").fillna("")
            self.assertEqual(raw, output.iloc[0][PAPER_FIELD])
            self.assertEqual("resolved", output.iloc[0].bibliography_state)
            self.assertEqual("not_evaluated", output.iloc[0].claim_support_summary)
            from candidate_qtl.relevance_tiers import classify_relevance
            self.assertTrue(classify_relevance(output.iloc[0].to_dict())[0].startswith("T1_"))

    def test_two_papers_do_not_fill_five_and_aliases_count_once(self):
        papers = [dict(report(function_supports_hypothesis="yes"), source_id="PMID:101", source_identifier_aliases="DOI:10.1234/one"),
                  dict(report(function_supports_hypothesis="yes"), source_id="DOI:10.1234/one", source_identifier_aliases="PMID:101"),
                  dict(report(function_supports_hypothesis="yes"), source_id="PMID:102")]
        selected = select_bibliography_papers({PAPER_FIELD: papers}, ["PMID:101", "PMID:102"])
        self.assertEqual(2, selected["selected_paper_count"])
        self.assertEqual([], selected["unselected_source_ids"])


class ResearchRefreezeTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.fixture = import_fixtures.RouteBSemanticImportTests()
        self.built = self.fixture._bundle(self.root)
        self.parent_path = Path(self.built["bundle"])
        self.parent_bytes = self.parent_path.read_bytes()
        self.parent = json.loads(self.parent_bytes)
        self.parent_hash = hashlib.sha256(self.parent_bytes).hexdigest()
        self.uid = self.parent["units"][0]["unit_id"]
        self.discovery = {"unit_id": self.uid, "reference_gene": "REF1", "reference_species": "Plant beta",
            "reported_gene": "GENE_B", "reported_gene_id": "GeneID:100", "reported_species": "Plant beta",
            "source_id": "PMID:7654321", "source_title": "Receptor regulates meristem size",
            "source_url": "https://example.invalid/primary", "retrieved_at": "2026-01-01T00:00:00Z",
            "source_text": "Loss of the receptor increased meristem size in Plant beta.",
            "reference_function": "receptor signaling restricts meristem size"}

    def artifact(self, **changes):
        payload = {"case_id": "case403", "parent_bundle_sha256": self.parent_hash,
                   "discovered_sources": [{**self.discovery, **changes}], "responses": []}
        path = self.root / "discovery.json"
        path.write_text(json.dumps(payload), encoding="utf-8")
        return path

    def test_cli_refreeze_new_paper_reaches_both_copies_without_support_promotion(self):
        artifact = self.artifact()
        # A discovery must never be silently merged as an old frozen response.
        with self.assertRaisesRegex(ValueError, "re-freeze"):
            merge(self.built["prompts"], [str(artifact)], str(self.root / "bad.json"))
        self.assertEqual(0, main(["route-b-semantic-bundle", "--extend-bundle", str(self.parent_path),
            "--parent-bundle-sha256", self.parent_hash, "--research-discoveries", str(artifact), "--out-dir", str(self.root / "new")]))
        new = self.root / "new/route_b_semantic_bundle.json"
        prompts = self.root / "new/route_b_semantic_prompts.jsonl"
        self.assertEqual(self.parent_bytes, self.parent_path.read_bytes())
        unit = json.loads(new.read_text())["units"][0]
        for key in ("genes", "identity", "member_interval_genes"):
            self.assertEqual(self.parent["units"][0][key], unit[key])
        self.assertIn("PMID:7654321", unit["stable_source_ids"])
        self.assertEqual(pd.read_csv(self.built["units"], sep="\t").columns.tolist(),
                         pd.read_csv(self.root / "new/route_b_semantic_units.tsv", sep="\t").columns.tolist())
        prepared = prepare(str(prompts), str(self.root / "prepared"))
        task = next((self.root / "prepared").rglob("worker_task.txt")).read_text()
        self.assertIn(prepared["research_source_parent"]["sha256"], task)
        paper = dict(report(phenotype_reported="yes", similar_cellular_or_developmental_phenotype="yes", focal_species="no"),
                     source_id="PMID:7654321", source_title=self.discovery["source_title"])
        claim = {"claim_type": "direct_functional_bridge", "biological_function": self.discovery["reference_function"],
                 "mechanism": "receptor signaling", "phenotype_relation": "could explain fasciation", "reference_gene": "REF1",
                 "reference_species": "Plant beta", "evidence_scope": "reference_gene", "source_ids": ["PMID:7654321"], PAPER_FIELD: [paper]}
        response = self.root / "response.json"
        self.fixture._response(prompts, [claim], response)
        before = import_route_b_semantic_proposals(self.built["bundle"], self.built["prompts"], str(response), str(self.root / "before"))
        self.assertEqual(1, before["unsupported_claim_row_count"])
        imported = import_route_b_semantic_proposals(str(new), str(prompts), str(response), str(self.root / "import"))
        self.assertEqual(0, imported["unsupported_claim_row_count"])
        cache = self.root / "cache.json"
        cache.write_text(json.dumps({"pmid": "7654321", "id": "7654321", "title": self.discovery["source_title"]}))
        sources = verify_route_b_semantic_sources(imported["proposals"], str(self.root / "sources"), bundle_arg=str(new), cache_arg=str(cache))
        transfer = verify_route_b_semantic_transfer(imported["proposals"], sources["source_verification"], str(new), str(self.root / "transfer"), case_id="case403")
        frozen = freeze_route_b_semantic_candidates(transfer["claims"], str(self.root / "freeze"), case_id="case403")
        self.assertEqual({"G1", "G2"}, set(pd.read_csv(frozen["review_candidates"], sep="\t").candidate_gene_id))
        self.assertEqual(0, frozen["candidate_gene_count"])

    def test_rejects_unknown_endpoint_unit_and_parent_mismatch_without_writing(self):
        for change in ({"unit_id": "unknown"}, {"reference_gene": "invented"}, {"source_text": ""}):
            with self.subTest(change=change):
                with self.assertRaises(ValueError):
                    refreeze_route_b_sources(str(self.parent_path), str(self.artifact(**change)), str(self.root / "rejected"), parent_sha256=self.parent_hash)
                self.assertFalse((self.root / "rejected").exists())
        with self.assertRaisesRegex(ValueError, "SHA-256"):
            refreeze_route_b_sources(str(self.parent_path), str(self.artifact()), str(self.root / "rejected"), parent_sha256="bad")

    def test_b_per_source_expansion_still_checks_only_five_and_retains_all_claims(self):
        discoveries, claims = [], []
        for number in range(701, 708):
            source, title = f"PMID:{number}", f"Synthetic paper {number}"
            discoveries.append({**self.discovery, "source_id": source, "source_title": title})
            answers = ({"phenotype_reported": "yes", "equivalent_target_phenotype": "yes", "focal_species": "yes", "same_gene_as_focal": "yes"}
                       if number == 707 else {"annotation_supports_hypothesis": "yes"})
            claims.append({"claim_type": "direct_functional_bridge", "biological_function": self.discovery["reference_function"],
                "mechanism": "receptor signaling", "phenotype_relation": "could explain fasciation", "reference_gene": "REF1",
                "reference_species": "Plant beta", "evidence_scope": "reference_gene", "source_ids": [source],
                PAPER_FIELD: [dict(report(**answers), source_id=source, source_title=title)]})
        artifact = self.artifact()
        payload = json.loads(artifact.read_text())
        payload["discovered_sources"] = discoveries
        artifact.write_text(json.dumps(payload))
        new = refreeze_route_b_sources(str(self.parent_path), str(artifact), str(self.root / "new"), parent_sha256=self.parent_hash)
        response = self.root / "response.json"
        self.fixture._response(Path(new["prompts"]), claims, response)
        imported = import_route_b_semantic_proposals(new["bundle"], new["prompts"], str(response), str(self.root / "import"))
        original = Path(imported["proposals"]).read_bytes()
        with patch("candidate_qtl.bibliography.fetch_metadata", side_effect=lambda kind, identifier: {"pmid": identifier, "title": f"Synthetic paper {identifier}"}) as fetch, \
             patch("candidate_qtl.gene_metadata.fetch_gene_metadata", return_value={"state": "not_found"}), \
             patch("candidate_qtl.source_verify._fetch_europe_pmc", side_effect=AssertionError("must not fetch abstracts")):
            sources = verify_route_b_semantic_sources(imported["proposals"], str(self.root / "sources"), bundle_arg=new["bundle"], allow_external_api=True)
        self.assertEqual(5, fetch.call_count)
        self.assertIn("707", [call.args[1] for call in fetch.call_args_list])
        rows = pd.read_csv(sources["source_verification"], sep="\t").fillna("")
        self.assertEqual(7, len(rows))
        self.assertEqual(2, rows.bibliography_state.eq("paper_handoff_audit_only").sum())
        self.assertEqual({"not_evaluated"}, set(rows.source_relation_status))
        self.assertEqual(original, Path(imported["proposals"]).read_bytes())
        transfer = verify_route_b_semantic_transfer(imported["proposals"], sources["source_verification"], new["bundle"], str(self.root / "transfer"), case_id="case403")
        frozen = freeze_route_b_semantic_candidates(transfer["claims"], str(self.root / "freeze"), case_id="case403")
        self.assertEqual({"G1", "G2"}, set(pd.read_csv(frozen["review_candidates"], sep="\t").candidate_gene_id))
        self.assertEqual(0, frozen["candidate_gene_count"])
