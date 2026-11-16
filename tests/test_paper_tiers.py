"""Paper-owned Tier, strongest-paper reduction, and actual A/B consumers."""
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.bibliography import INITIAL_POLICY
from candidate_qtl.boundary_tier_review import EVIDENCE_SCOPED_POLICY, freeze_boundary_scope, build_boundary_review_request
from candidate_qtl.cli import main
from candidate_qtl.paper_reports import PAPER_FIELD, imported_paper_fields, normalize_paper_reports
from candidate_qtl.relevance_tiers import NEW_TIERS, annotate_relevance, classify_relevance
from candidate_qtl.route_b_semantic_import import import_route_b_semantic_proposals
from candidate_qtl.route_b_semantic_source_verify import verify_route_b_semantic_sources
from candidate_qtl.route_b_semantic_transfer_verify import verify_route_b_semantic_transfer
from candidate_qtl.route_b_semantic_candidates import freeze_route_b_semantic_candidates
from candidate_qtl.source_verify import source_verify
from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from tests.test_initial_tiers import observation
from tests import test_route_b_semantic_import as b_fixtures
from tests.test_phenotype_investigation import MECHANISMS, _payload


def paper(source="PMID:123", tier="T1", **overrides):
    return {"source_id": source, "initial_tier": tier, "initial_tier_reason": "Observation in this paper",
            "reported_gene": "REF1", "reported_species": "Plant alpha", "reported_phenotype": "enlarged meristem",
            "reported_process": "meristem proliferation", "phenotype_comparison": "same_equivalent" if tier == "T1" else "shared_process",
            "source_title": "Paper " + source, **overrides}


def claim(reports, **overrides):
    return observation("T2", source_id="", source_ids=";".join(r["source_id"] for r in reports),
                       **{PAPER_FIELD: normalize_paper_reports(reports)}, **overrides)


def projected_claim(reports, gene="G1", **overrides):
    """A member-bound claim carrying the qualified projection #7 requires for T1.

    Paper-side checklist T1 is conditional on a member projection; this mirrors
    the output of the existing qualified-sequence-hit producer for one copy.
    """
    reference = next((r["reference_gene_id"] for r in reports if r.get("reference_gene_id")), "REF1")
    projection = {"basis": "qualified_sequence_hit", "route": "a", "candidate_gene_id": gene,
                  "member_genes": [gene], "copy_count": 1, "reference_gene_id": reference,
                  "evidence": {"pident": 90.0, "query_coverage": 0.9, "target_coverage": 0.9}}
    return claim(reports, candidate_gene_id=gene, homolog_projection_json=json.dumps(projection), **overrides)


class PaperTierTests(unittest.TestCase):
    def test_t1_review_only_and_original_with_three_reviews_have_no_count_condition(self):
        review = paper(publication_type="review")
        for reports in ([review], [paper(publication_type="original"),
                                  *[paper(f"PMID:{i}", publication_type="review", underlying_source_ids="PMID:123")
                                    for i in range(501, 504)]]):
            with self.subTest(count=len(reports)):
                row = annotate_relevance([{"candidate_gene_id": "G1"}], [claim(reports)])[0]
                self.assertEqual(NEW_TIERS[0], row["relevance_tier"])
                self.assertEqual({"T1": len(reports)}, json.loads(row["paper_tier_counts_json"]))

    def test_multiple_originals_and_reviews_keep_tier_and_show_distinct_counts(self):
        reports = [paper(publication_type="original"),
                   paper("PMID:456", publication_type="review", underlying_source_ids="PMID:123"),
                   *[paper(f"PMID:{i}", "T2", phenotype_comparison="same_cellular_process") for i in range(701, 704)]]
        first = annotate_relevance([{"candidate_gene_id": "G1"}], [claim(reports), claim(reports, route="a")])[0]
        self.assertEqual(NEW_TIERS[0], first["relevance_tier"])
        self.assertEqual({"T1": 2, "T2": 3}, json.loads(first["paper_tier_counts_json"]))
        summaries = json.loads(first["paper_evidence_summary_json"])
        review = next(row for row in summaries if row["source_id"] == "PMID:456")
        self.assertEqual(["review"], review["publication_type"])
        self.assertEqual(["PMID:123"], review["underlying_source_ids"])
        self.assertEqual("not_assessed", review["independence_state"])
        replay = annotate_relevance([first])[0]
        self.assertEqual(first["paper_tier_counts_json"], replay["paper_tier_counts_json"])
        self.assertEqual(first["paper_evidence_summary_json"], replay["paper_evidence_summary_json"])

    def test_pmid_doi_alias_and_repeated_route_are_counted_as_one_paper(self):
        reports = [paper(source_identifier_aliases="DOI:10.1234/first"),
                   paper("doi:10.1234/FIRST", publication_type="review")]
        row = annotate_relevance([{"candidate_gene_id": "G1"}], [claim(reports), claim(reports, route="b")])[0]
        self.assertEqual({"T1": 1}, json.loads(row["paper_tier_counts_json"]))
        self.assertEqual(1, len(json.loads(row["paper_evidence_summary_json"])))

    def test_strongest_is_invariant_to_weak_paper_count_route_duplicates_and_order(self):
        strong = paper()
        weak = [paper(f"PMID:{i}", "T3") for i in range(1000, 1100)]
        roster = [{"candidate_gene_id": "G1", "support_status": "unverified_route_b_review"}]
        alone = annotate_relevance(roster, [claim([strong])])[0]
        mixed = annotate_relevance(roster, [claim([*weak, strong]), claim([strong], route="a")])[0]
        reversed_row = annotate_relevance(roster, [claim([strong], route="a"), claim([strong, *reversed(weak)])])[0]
        for field in ("relevance_tier", "initial_tier", "tier_reason", "tier_source_id", "tier_evidence_json", "support_state"):
            self.assertEqual(alone[field], mixed[field], field)
            self.assertEqual(mixed[field], reversed_row[field], field)
        self.assertEqual(NEW_TIERS[0], mixed["relevance_tier"])
        self.assertEqual("review", mixed["support_state"])
        self.assertEqual({"T1": 1, "T3": 100}, json.loads(mixed["paper_tier_counts_json"]))
        self.assertEqual(101, len(json.loads(mixed["paper_evidence_summary_json"])))
        self.assertEqual(101, max(len(json.loads(row[PAPER_FIELD])) for row in json.loads(mixed["relevance_evidence_json"])))

    def test_rehydration_preserves_tier_and_does_not_expand_to_other_copies(self):
        row = annotate_relevance([{"candidate_gene_id": "G1"}], [claim([paper(), paper("PMID:456", "T3")])])[0]
        replay = annotate_relevance([{"candidate_gene_id": "G1"}, {"candidate_gene_id": "G2", "family_id": "F"}], [row])
        self.assertEqual(row["tier_evidence_json"], replay[0]["tier_evidence_json"])
        self.assertEqual(NEW_TIERS[0], replay[0]["relevance_tier"])
        self.assertEqual(NEW_TIERS[5], replay[1]["relevance_tier"])

    def test_per_source_binding_prevents_tier_and_title_broadcast(self):
        record = {"source_ids": ["PMID:123", "PMID:456", "UniProt:P12345"],
                  PAPER_FIELD: [paper(), paper("PMID:456", "T3")]}
        for source, tier in [("PMID:123", "T1"), ("PMID:456", "T3"), ("UniProt:P12345", "")]:
            bound = imported_paper_fields(record, source)
            self.assertEqual(tier, bound["initial_tier"])
            if tier:
                self.assertEqual("Paper " + source, bound["source_title"])
            else:
                self.assertEqual("[]", bound[PAPER_FIELD])
                self.assertEqual("", bound["source_title"])

    def test_empty_or_unbound_paper_cannot_inherit_a_hypothesis_t1(self):
        for source, expected in [("PMID:123", NEW_TIERS[5]), ("UniProt:P12345", NEW_TIERS[4])]:
            row = observation("T1", source_id=source, phenotype_comparison="same_equivalent", **{PAPER_FIELD: "[]"})
            result = annotate_relevance([{"candidate_gene_id": "G1"}], [row])[0]
            self.assertEqual(expected, result["relevance_tier"])
            self.assertEqual(expected.split("_", 1)[0], result["initial_tier"])
        row = claim([paper()])
        row["source_ids"] = "PMID:999"
        self.assertEqual(NEW_TIERS[5], classify_relevance(row)[0])

    def test_reported_species_is_not_borrowed_from_the_projection_reference(self):
        row = claim([paper(reported_species="Plant beta")])
        self.assertEqual(NEW_TIERS[3], classify_relevance(row)[0])
        row = claim([paper(tier="T6")])
        self.assertEqual(NEW_TIERS[5], classify_relevance(row)[0])

    def test_invalid_or_undeclared_paper_is_rejected(self):
        for reports in [[paper(source="UniProt:P12345")], [paper(initial_tier_reason="")],
                        [paper(reported_gene="")], [paper(source_identifier_aliases="UniProt:P12345")], "not json"]:
            with self.subTest(reports=reports), self.assertRaises(ValueError):
                normalize_paper_reports(reports)
        with self.assertRaisesRegex(ValueError, "not bound"):
            imported_paper_fields({"source_ids": ["PMID:999"], PAPER_FIELD: [paper()]})

    def test_a_cli_preserves_separate_reports_and_checks_each_papers_title(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reports = [paper(), paper("PMID:456", "T3", source_title="Wrong second title")]
            record = {**MECHANISMS[0], "source_ids": "PMID:123;PMID:456", PAPER_FIELD: reports}
            self.assertEqual(0, main(["mechanism-import", "--payload", str(_payload(root, [record])), "--out-dir", str(root / "import")]))
            table = root / "import/mechanism_hypotheses.tsv"
            rows = pd.read_csv(table, sep="\t").fillna("")
            self.assertEqual(2, len(json.loads(rows.iloc[0][PAPER_FIELD])))
            cache = root / "cache.json"
            cache.write_text(json.dumps([{"id": "123", "title": "Paper PMID:123"}, {"id": "456", "title": "Paper PMID:456"}]))
            checked = source_verify(str(table), str(root / "sources"), cache_arg=str(cache), verification_policy=INITIAL_POLICY)
            validation = pd.read_csv(checked["source_validation"], sep="\t").set_index("citation")
            self.assertEqual("id_title_matched", validation.loc["PMID:123", "bibliography_state"])
            self.assertEqual("title_mismatch", validation.loc["PMID:456", "bibliography_state"])

    def test_paper_aliases_share_identity_checks_without_combining_different_papers(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reports = [paper(source_identifier_aliases="DOI:10.1234/first"), paper("PMID:456", "T3")]
            table = root / "reports.tsv"
            pd.DataFrame([{"hypothesis_id": "h1", "source_ids": "PMID:123;PMID:456",
                           PAPER_FIELD: normalize_paper_reports(reports)}]).to_csv(table, sep="\t", index=False)
            cache = root / "cache.json"
            cache.write_text(json.dumps([{"id": "123", "doi": "10.1234/first", "title": "Paper PMID:123"},
                                         {"id": "456", "title": "Paper PMID:456"}]))
            checked = source_verify(str(table), str(root / "good"), cache_arg=str(cache), verification_policy=INITIAL_POLICY)
            validation = pd.read_csv(checked["source_validation"], sep="\t")
            self.assertEqual(3, len(validation))
            self.assertEqual({"id_title_matched"}, set(validation.bibliography_state))
            reports[0]["source_identifier_aliases"] = "PMID:456"
            reports[1]["source_title"] = reports[0]["source_title"]
            # No titles makes this specifically an alias identity test.
            for report in reports:
                report["source_title"] = ""
            pd.DataFrame([{"hypothesis_id": "h1", "source_ids": "PMID:123;PMID:456",
                           PAPER_FIELD: normalize_paper_reports(reports)}]).to_csv(table, sep="\t", index=False)
            checked = source_verify(str(table), str(root / "bad"), cache_arg=str(cache), verification_policy=INITIAL_POLICY)
            validation = pd.read_csv(checked["source_validation"], sep="\t")
            self.assertEqual({"multiple_identifiers_mismatch"}, set(validation.bibliography_state))

    def test_b_model_response_reaches_source_transfer_and_freeze_with_paper_tiers(self):
        self._assert_b_model_response_reaches_wet()

    def _assert_b_model_response_reaches_wet(self, machine=False, semantic=False, checklist=False):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = b_fixtures.RouteBSemanticImportTests()
            built = fixture._bundle(root)
            bundle = json.loads(Path(built["bundle"]).read_text())
            for unit in bundle["units"]:
                unit["stable_source_ids"].append("PMID:456")
            Path(built["bundle"]).write_text(json.dumps(bundle))
            response = root / "response.json"
            reports = [paper(tier="T2", phenotype_comparison="same_cellular_process"), paper("PMID:456", "T3")]
            if semantic:
                self.assertIn("reference_gene_id", Path(built["prompts"]).read_text())
                self.assertIn("gene_correspondence", Path(built["prompts"]).read_text())
                for report in reports:
                    report.update(reference_gene_id="REF1", gene_correspondence="same_gene", reported_gene="NCBI Gene:111111",
                                  reported_species="source common name", evidence_types=["phenotype"])
            if machine:
                for report in reports:
                    del report["initial_tier"]
                    del report["initial_tier_reason"]
                    report["evidence_types"] = ["phenotype"]
            if checklist:
                for report in reports:
                    label = report.pop("initial_tier")
                    report.update(reference_gene_id="REF1", gene_correspondence="homolog",
                                  reported_gene="Brassica locus X", reported_gene_id="GenBank:AB123456",
                                  reported_species="Brassica rapa",
                                  initial_tier_reason="Brassica phenotype report is connected through supplied REF1 as a sequence proxy.",
                                  paper_checklist={"phenotype_reported": "yes", "identity_resolved": "yes",
                                                   "similar_cellular_or_developmental_phenotype" if label == "T2" else "shared_specific_process": "yes"})
            model_claim = {**observation(), "claim_type": "mutant_phenotype", "reference_gene": "REF1",
                           "evidence_scope": "reference_gene", "source_ids": ["PMID:123", "PMID:456"], PAPER_FIELD: reports}
            fixture._response(Path(built["prompts"]), [model_claim], response)
            imported = import_route_b_semantic_proposals(built["bundle"], built["prompts"], str(response), str(root / "import"))
            rows = pd.read_csv(imported["proposals"], sep="\t").fillna("")
            self.assertEqual({""} if machine else {"T2", "T3"}, set(rows.initial_tier))
            self.assertTrue(all(len(json.loads(value)) == 1 for value in rows[PAPER_FIELD]))
            cache = root / "cache.json"
            cache.write_text(json.dumps([{"id": "123", "title": "Paper PMID:123"}, {"id": "456", "title": "Paper PMID:456"}]))
            if checklist:
                cache.write_text(json.dumps({"records": json.loads(cache.read_text()), "gene_records": {
                    "embl:AB123456": {"ids": ["GenBank:AB123456"], "names": ["Brassica locus X"], "species": ["Brassica rapa"]}}}))
            sources = verify_route_b_semantic_sources(imported["proposals"], str(root / "sources"), bundle_arg=built["bundle"], cache_arg=str(cache))
            if checklist:
                source_rows = pd.read_csv(sources["source_verification"], sep="\t").fillna("")
                self.assertTrue(all(json.loads(value)[0]["state"] == "matched" for value in source_rows.gene_identity_checks_json))
            transfer = verify_route_b_semantic_transfer(imported["proposals"], sources["source_verification"], built["bundle"], str(root / "transfer"), case_id="case403")
            frozen = freeze_route_b_semantic_candidates(transfer["claims"], str(root / "freeze"), case_id="case403")
            review = pd.read_csv(frozen["review_candidates"], sep="\t").fillna("")
            self.assertEqual({"G1", "G2"}, set(review.candidate_gene_id))
            self.assertEqual({NEW_TIERS[1], NEW_TIERS[2]}, set(review.relevance_tier))
            self.assertEqual({"PMID:123", "PMID:456"}, set(review.tier_source_id))
            self.assertEqual({"review"}, set(review.support_state))
            self.assertEqual(0, frozen["candidate_gene_count"])
            atlas = root / "interval1000.tsv"
            pd.DataFrame([{"gene_id": f"G{i}", "family_id": "F"} for i in range(1, 1001)]).to_csv(atlas, sep="\t", index=False)
            wet = wet_candidate_handoff(str(root / "wet"), route_b_candidates_arg=frozen["candidates"],
                route_b_candidates_manifest_arg=frozen["manifest"], review_candidates_arg=frozen["review_candidates"],
                review_candidates_manifest_arg=frozen["review_manifest"], interval_atlas_arg=str(atlas), case_id="case403")
            adopted = pd.read_csv(wet["wet_candidates"], sep="\t").fillna("")
            self.assertEqual({"G1", "G2"}, set(adopted.candidate_gene_id))
            self.assertEqual({NEW_TIERS[1]}, set(adopted.relevance_tier))
            self.assertEqual({"PMID:123"}, set(adopted.tier_source_id))
            self.assertTrue(all(json.loads(value) == {"T2": 1, "T3": 1} for value in adopted.paper_tier_counts_json))

    def test_1000_gene_inventory_reaches_hash_bound_boundary_with_paper_policy(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            roster = [{"candidate_gene_id": f"G{i}", "family_id": "F"} for i in range(1000)]
            claims = [claim([paper(tier="T3")], candidate_gene_id=f"G{i}") for i in range(30)]
            rows = annotate_relevance(roster, claims)
            self.assertEqual(30, sum(row["relevance_tier"] == NEW_TIERS[2] for row in rows))
            table = root / "inventory.tsv"
            pd.DataFrame(rows).fillna("").to_csv(table, sep="\t", index=False)
            scoped = freeze_boundary_scope(str(table), str(root / "scope"),
                inventory_sha256=hashlib.sha256(table.read_bytes()).hexdigest(), case_id="synthetic-paper-policy",
                policy_version=EVIDENCE_SCOPED_POLICY, eligible_gene_ids={f"G{i}" for i in range(30)})
            request = build_boundary_review_request(scoped["manifest"], str(root / "request"),
                scope_sha256=scoped["manifest_sha256"], case_id="synthetic-paper-policy")
            value = json.loads(Path(request["request"]).read_text())
            self.assertEqual(30, value["candidate_count"])
            self.assertTrue(any("Paper counts" in rule for rule in value["policy"]))
            bound = pd.read_csv(value["candidate_input"]["path"], sep="\t")
            self.assertEqual({"PMID:123"}, set(bound.tier_source_id))


if __name__ == "__main__":
    unittest.main()
