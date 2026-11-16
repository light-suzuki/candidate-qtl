"""Meaning-preserving formats and negative identity controls at the A/B boundary."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl.bibliography import normalize_title, verify_bibliography
from candidate_qtl.cli import main
from candidate_qtl.observation_tiers import gene_binding
from candidate_qtl.paper_evidence_table import write_paper_evidence_table
from candidate_qtl.paper_reports import PAPER_FIELD, imported_paper_fields, normalize_paper_reports
from candidate_qtl.relevance_tiers import classify_relevance
from candidate_qtl.route_b_semantic_import import import_route_b_semantic_proposals
from candidate_qtl.route_b_semantic_source_verify import verify_route_b_semantic_sources
from candidate_qtl.source_verify import source_verify
from tests.test_paper_evidence_table import paper
from tests import test_paper_tiers as a
from tests import test_route_b_semantic_import as b


class PaperMatchingEdgeTests(unittest.TestCase):
    def test_same_citation_formats_bind_without_a_declared_cross_paper_alias(self):
        pairs = [("DOI:10.1234/AbC", "https://doi.org/10.1234/abc"),
                 ("https://doi.org/10.1234/abc", "10.1234/ABC"),
                 ("PMID:12345", "12345"), ("PMID:12345", "https://pubmed.ncbi.nlm.nih.gov/12345/")]
        for reported, declared in pairs:
            with self.subTest(reported=reported, declared=declared):
                p = {**paper(phenotype_reported="yes", similar_cellular_or_developmental_phenotype="yes"), "source_id": reported}
                record = {"source_ids": [declared], PAPER_FIELD: [p]}
                result = imported_paper_fields(record, declared)
                self.assertEqual("T2", result["initial_tier"])
                claim = {**a.claim([p]), "source_id": declared, "source_ids": declared}
                self.assertEqual("T2_reported_similar_phenomenon", classify_relevance(claim)[0])

    def test_unknown_gene_never_authorizes_t1_in_any_display_case(self):
        for missing in ("UNKNOWN", " Unknown ", "unresolved", "N/A", "不明"):
            with self.subTest(missing=missing):
                p = paper(phenotype_reported="yes", equivalent_target_phenotype="yes", focal_species="yes", same_gene_as_focal="yes")
                p.update(reported_gene=missing, unknown_reason="The reported gene cannot be identified from the inspected passage.")
                self.assertEqual("U_unresolved", classify_relevance(a.claim([p]))[0])
                self.assertNotEqual("matched", gene_binding({"reference_gene_name": missing}, missing))

    def test_typed_gene_id_namespaces_and_versions_are_not_discarded(self):
        for reference, reported in (("GeneID:12345", "NCBI Gene:12345"),
                                    ("UniProt:P12345", "sp|P12345|SYNTHETIC"),
                                    ("TAIR:AT1G12345", "AT1G12345"),
                                    ("GenBank:AB123456.1", "EMBL:AB123456.1")):
            with self.subTest(reference=reference, reported=reported):
                self.assertEqual("matched", gene_binding({"reference_gene": reference}, reported))
        for reference, reported in (("UniProt:P12345", "RefSeq:P12345"),
                                    ("GeneID:12345", "GeneID:12346"),
                                    ("RefSeq:NP_123456.1", "RefSeq:NP_123456.2")):
            with self.subTest(reference=reference, reported=reported):
                self.assertEqual("unresolved", gene_binding({"reference_gene": reference}, reported))

    def test_shared_symbol_cannot_hide_a_conflicting_stable_accession(self):
        row = {"reference_gene_name": "GENEA", "proposed_stable_id": "P12345"}
        for reported in ("GENEA|Q12345", "GENEA (UniProt:Q12345)", "GENEA (RefSeq:P12345)"):
            with self.subTest(reported=reported):
                self.assertEqual("unresolved", gene_binding(row, reported))
        self.assertEqual("matched", gene_binding(row, "GENEA"))

    def test_omitted_null_and_explicit_unknown_answers_have_different_recording_states(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "raw.json"
            source.write_text("{}")
            for answer, expected in ((None, "complete_record"), (" UNKNOWN ", "none"), ("不明", "none")):
                with self.subTest(answer=answer):
                    p = paper(function_supports_hypothesis=answer, annotation_supports_hypothesis="yes")
                    p["unknown_reason"] = "Results checked; function remains unknown."
                    result = write_paper_evidence_table([{"hypothesis_id": "h", "reference_gene_name": "REF1",
                        "source_ids": "PMID:123", PAPER_FIELD: [p]}], root, route="phenotype_first", source_artifact=str(source))
                    row = pd.read_csv(result["outputs"]["paper_evidence_table"], sep="\t").fillna("").iloc[0]
                    self.assertEqual(expected, row.followup_action)
                    self.assertEqual("T5", row.computed_tier)

    def test_formatted_enums_and_optional_nulls_reach_same_t2(self):
        p = paper(phenotype_reported=" TRUE ", similar_cellular_or_developmental_phenotype=True)
        p.update(gene_correspondence=" SAME-GENE ", followup_state=" Completed ",
                 publication_type=" Review ", unknown_reason=None, followup_question=None)
        result = json.loads(normalize_paper_reports([p]))[0]
        self.assertEqual("T2", result["initial_tier"])
        self.assertEqual("same_gene", result["gene_correspondence"])
        self.assertEqual("completed", result["followup_state"])

    def test_title_markup_is_display_format_but_different_title_stays_different(self):
        self.assertEqual(normalize_title("Regulation of <i>GENEA</i> &amp; growth"),
                         normalize_title("Regulation of GENEA & growth"))
        self.assertEqual(normalize_title("ＡＢＣ — gene regulation"), normalize_title("ABC: gene regulation"))
        self.assertNotEqual(normalize_title("GENEA promotes growth"), normalize_title("GENEA inhibits growth"))
        self.assertEqual(normalize_title("<i>GENEA</i>1 function"), normalize_title("GENEA1 function"))

    def test_different_citation_cannot_borrow_another_papers_tier(self):
        p = dict(paper(phenotype_reported="yes", similar_cellular_or_developmental_phenotype="yes"), source_id="DOI:10.1234/one")
        for other in ("DOI:10.1234/two", "PMID:12345"):
            with self.subTest(other=other), self.assertRaisesRegex(ValueError, "not bound"):
                imported_paper_fields({"source_ids": [other], PAPER_FIELD: [p]})

    def test_live_bibliography_path_reuses_cached_doi_url_and_title_markup(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            p = dict(paper(phenotype_reported="yes", similar_cellular_or_developmental_phenotype="yes"),
                     source_id="https://doi.org/10.1234/ABC", source_title="<i>GENEA</i>1 &amp; development")
            raw = a._payload(root, [{**a.MECHANISMS[0], "reference_gene_name": "REF1",
                                    "source_ids": "DOI:10.1234/abc", PAPER_FIELD: [p]}])
            self.assertEqual(0, main(["mechanism-import", "--payload", str(raw), "--out-dir", str(root / "import")]))
            source = root / "import/mechanism_hypotheses.tsv"
            cache = root / "cache.json"
            cache.write_text(json.dumps([{"doi": "https://doi.org/10.1234/abc", "title": "GENEA1 & development"}]))
            result = verify_bibliography(str(source), str(root / "checked"), cache_arg=str(cache), allow_external_api=False,
                                         id_column="hypothesis_id", source_column="source_ids")
            row = pd.read_csv(result["paper_evidence_table"], sep="\t").iloc[0]
            self.assertEqual("id_title_matched", row.metadata_state)
            self.assertEqual("T2", row.applicable_tier)

    def test_live_metadata_aliases_do_not_repeat_same_paper_lookup(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            p = dict(paper(phenotype_reported="yes", similar_cellular_or_developmental_phenotype="yes"),
                     source_id="PMID:12345", source_identifier_aliases=["https://doi.org/10.1234/ABC"])
            source = root / "input.tsv"
            pd.DataFrame([{"hypothesis_id": "h", "reference_gene_name": "REF1", "source_ids": "PMID:12345",
                           PAPER_FIELD: json.dumps([p])}]).to_csv(source, sep="\t", index=False)
            with patch("candidate_qtl.bibliography.fetch_metadata", return_value={"pmid": "12345", "doi": "https://doi.org/10.1234/abc", "title": p["source_title"]}) as fetch, \
                 patch("candidate_qtl.gene_metadata.fetch_gene_metadata", return_value={"status": "unconfirmed"}):
                result = verify_bibliography(str(source), str(root / "checked"), cache_arg=None, allow_external_api=True,
                                             id_column="hypothesis_id", source_column="source_ids")
            self.assertEqual(1, fetch.call_count)
            rows = pd.read_csv(result["source_validation"], sep="\t")
            self.assertEqual({"id_title_matched"}, set(rows.bibliography_state))

    def test_cached_review_material_keeps_doi_url_without_promoting_claim_support(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "input.tsv"
            pd.DataFrame([{"hypothesis_id": "h", "source_ids": "", "review_only_source_ids": "DOI:10.1234/one"}]).to_csv(source, sep="\t", index=False)
            cache = root / "cache.json"
            cache.write_text(json.dumps([{"doi": "https://doi.org/10.1234/one", "title": "Synthetic review material",
                                          "abstractText": "An experimentally observed developmental change."}]))
            result = source_verify(str(source), str(root / "checked"), cache_arg=str(cache), verification_policy="content_v1")
            records = json.loads(Path(result["source_review_records"]).read_text())
            self.assertIn("DOI:10.1234/one", records["records"])
            self.assertEqual("none", records["membership_effect"])
            self.assertEqual(0, result["verified_citation_count"])

    def test_same_checklist_meaning_preserves_every_tier_across_input_formats(self):
        cases = [(dict(phenotype_reported="yes", equivalent_target_phenotype="yes", focal_species="yes", same_gene_as_focal="yes"), "T1"),
                 (dict(phenotype_reported="yes", similar_cellular_or_developmental_phenotype="yes"), "T2"),
                 (dict(phenotype_reported="yes", shared_specific_process="yes"), "T3"),
                 (dict(function_supports_hypothesis="yes"), "T4"), (dict(annotation_supports_hypothesis="yes"), "T5"),
                 ({}, "T6")]
        for answers, expected in cases:
            original = paper(**answers)
            for form in ("booleans", "strings", "japanese"):
                p = json.loads(json.dumps(original))
                conversions = {"booleans": {"yes": True, "no": False},
                               "strings": {"yes": " TRUE ", "no": " False "},
                               "japanese": {"yes": "はい", "no": "いいえ"}}[form]
                p["paper_checklist"] = {key: conversions.get(value, value) for key, value in p["paper_checklist"].items()}
                p["gene_correspondence"] = " SAME GENE "
                with self.subTest(tier=expected, form=form):
                    result = json.loads(normalize_paper_reports([p]))[0]
                    self.assertEqual(expected, result["initial_tier"])
                    self.assertEqual(expected, classify_relevance(a.projected_claim([p]))[0].split("_", 1)[0])

    def test_b_reference_formats_map_only_to_the_same_frozen_typed_id(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = b.RouteBSemanticImportTests()
            built = fixture._bundle(root)
            bundle = json.loads(Path(built["bundle"]).read_text())
            for unit in bundle["units"]:
                for gene in unit["genes"]:
                    gene["orthology_reference_gene"] = "sp|P12345|SYNTHETIC"
            Path(built["bundle"]).write_text(json.dumps(bundle))
            for name, reference, expected in (("same", "UniProt:P12345", 1), ("different", "RefSeq:P12345", 0)):
                claim = {"claim_type": "mutant_phenotype", "reference_gene": reference,
                         "evidence_scope": "reference_gene", "source_ids": ["PMID:123"], PAPER_FIELD: []}
                response = root / f"{name}.json"
                fixture._response(Path(built["prompts"]), [claim], response)
                result = import_route_b_semantic_proposals(built["bundle"], built["prompts"], str(response), str(root / name))
                self.assertEqual(expected, result["proposed_claim_row_count"])

    def test_b_frozen_source_formats_are_reused_without_admitting_an_unfrozen_paper(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = b.RouteBSemanticImportTests()
            built = fixture._bundle(root)
            bundle = json.loads(Path(built["bundle"]).read_text())
            bundle["b_review_membership_policy"] = "initial_bibliography_v1"
            Path(built["bundle"]).write_text(json.dumps(bundle))
            p = dict(paper(phenotype_reported="yes", similar_cellular_or_developmental_phenotype="yes"),
                     source_id="https://pubmed.ncbi.nlm.nih.gov/123/")
            claim = {"claim_type": "mutant_phenotype", "reference_gene": "REF1", "evidence_scope": "reference_gene",
                     "source_ids": ["PMID:123", "https://pubmed.ncbi.nlm.nih.gov/123/"], PAPER_FIELD: [p]}
            source = root / "response.json"
            fixture._response(Path(built["prompts"]), [claim], source)
            result = import_route_b_semantic_proposals(built["bundle"], built["prompts"], str(source), str(root / "import"))
            self.assertEqual(1, result["proposed_claim_row_count"])
            proposals = pd.read_csv(result["proposals"], sep="\t")
            self.assertEqual({"PMID:123"}, set(proposals.source_id))
            cache = root / "cache.json"
            cache.write_text(json.dumps([{"id": "123", "title": p["source_title"]}]))
            checked = verify_route_b_semantic_sources(result["proposals"], str(root / "checked"),
                bundle_arg=built["bundle"], cache_arg=str(cache))
            rows = pd.read_csv(root / "checked/source_verify/source_validation.tsv", sep="\t")
            self.assertEqual(1, len(rows))
            self.assertEqual({"id_title_matched"}, set(rows.bibliography_state))
            self.assertEqual(0, checked["final_verified_claim_count"])
            claim.update(source_ids=["PMID:999"], paper_reports_json=[dict(p, source_id="PMID:999")])
            fixture._response(Path(built["prompts"]), [claim], root / "unfrozen.json")
            unsupported = import_route_b_semantic_proposals(built["bundle"], built["prompts"], str(root / "unfrozen.json"), str(root / "unfrozen"))
            self.assertEqual(0, unsupported["proposed_claim_row_count"])
            self.assertEqual(1, unsupported["unsupported_claim_row_count"])

    def test_real_a_targeted_reimport_closes_request_without_touching_original(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            p = paper(phenotype_reported="unknown", similar_cellular_or_developmental_phenotype="yes", annotation_supports_hypothesis="yes")
            p["unknown_reason"] = "Experiment details not yet available."
            payload = {**a.MECHANISMS[0], "reference_gene_name": "REF1", "source_ids": "PMID:123", PAPER_FIELD: [p]}
            raw = a._payload(root, [payload])
            original = raw.read_bytes()
            self.assertEqual(0, main(["mechanism-import", "--payload", str(raw), "--out-dir", str(root / "first")]))
            first_table = (root / "first/paper_evidence_table.tsv").read_bytes()
            packet = json.loads((root / "first/paper_followup_requests.json").read_text())
            self.assertEqual(["phenotype_reported"], packet["requests"][0]["records"][0]["questions"])
            p["paper_checklist"]["phenotype_reported"] = "yes"
            p.update(followup_state="completed", unknown_reason="")
            corrected = root / "corrected.json"
            corrected.write_text(json.dumps({"hypotheses": [payload]}))
            self.assertEqual(0, main(["mechanism-import", "--payload", str(corrected), "--out-dir", str(root / "corrected")]))
            final = pd.read_csv(root / "corrected/paper_evidence_table.tsv", sep="\t").iloc[0]
            self.assertEqual("T2", final.applicable_tier)
            self.assertEqual("none", final.followup_action)
            self.assertEqual(original, raw.read_bytes())
            self.assertEqual(first_table, (root / "first/paper_evidence_table.tsv").read_bytes())

    def test_real_b_checked_unknown_reimport_stops_without_faking_t2(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = b.RouteBSemanticImportTests()
            built = fixture._bundle(root)
            p = paper(phenotype_reported="unknown", similar_cellular_or_developmental_phenotype="yes", annotation_supports_hypothesis="yes")
            p["unknown_reason"] = "Experiment details unavailable."
            claim = {"claim_type": "mutant_phenotype", "reference_gene": "REF1", "evidence_scope": "reference_gene",
                     "source_ids": ["PMID:123"], PAPER_FIELD: [p]}
            source = root / "response.json"
            fixture._response(Path(built["prompts"]), [claim], source)
            initial = import_route_b_semantic_proposals(built["bundle"], built["prompts"], str(source), str(root / "first"))
            original = source.read_bytes()
            packet = json.loads(Path(initial["paper_followup_requests"]).read_text())
            self.assertTrue(packet["requests"])
            p.update(followup_state="completed", unknown_reason="Abstract and results inspected; the attribution is still uncertain.")
            corrected = root / "corrected.json"
            fixture._response(Path(built["prompts"]), [claim], corrected)
            result = import_route_b_semantic_proposals(built["bundle"], built["prompts"], str(corrected), str(root / "corrected"))
            row = pd.read_csv(result["paper_evidence_table"], sep="\t").iloc[0]
            self.assertEqual("T5", row.applicable_tier)
            self.assertEqual("none", row.followup_action)
            self.assertEqual(0, result["verified_claim_row_count"])
            self.assertEqual(original, source.read_bytes())
