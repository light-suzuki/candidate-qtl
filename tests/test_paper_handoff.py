"""Import owns the five-paper handoff; all observations retain their old meaning."""
import json
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.paper_reports import (
    PAPER_FIELD, PAPER_HANDOFF_FIELDS, PAPER_HANDOFF_POLICY,
    imported_paper_fields, normalize_paper_reports, select_bibliography_papers,
    source_tokens, stamp_paper_handoff,
)
from candidate_qtl.relevance_tiers import annotate_relevance
from candidate_qtl.route_b_semantic_bundle import build_route_b_semantic_bundle
from candidate_qtl.route_b_semantic_import import import_route_b_semantic_proposals
from tests import test_paper_tiers as a
from tests import test_route_b_semantic_import as b
from tests.test_paper_evidence_table import paper


def tier_paper(tier, number):
    answers = {
        "T1": dict(phenotype_reported="yes", equivalent_target_phenotype="yes", focal_species="yes", same_gene_as_focal="yes"),
        "T2": dict(phenotype_reported="yes", similar_cellular_or_developmental_phenotype="yes"),
        "T3": dict(phenotype_reported="yes", shared_specific_process="yes"),
        "T4": dict(function_supports_hypothesis="yes"),
        "T5": dict(annotation_supports_hypothesis="yes"),
        "T6": {},
    }
    return {**paper(**answers[tier]), "source_id": f"PMID:{number}",
            "source_title": f"Synthetic paper {number}",
            "initial_tier_reason": f"Original recorded observations satisfy {tier}."}


def twenty_papers():
    return [tier_paper(tier, 10001 + index) for index, tier in enumerate(["T6"] * 15 + ["T5", "T4", "T3", "T2", "T1"])]


def a_row(reports, **changes):
    return {"hypothesis_id": "h1", "case_id": "synthetic", "proposed_stable_id": "REF1",
            "reference_gene_name": "REF1", "source_ids": ";".join(report["source_id"] for report in reports),
            "counter_hypothesis": "A distinct regulatory mechanism remains possible.",
            PAPER_FIELD: reports, **changes}


def reports_in(rows):
    return [report for row in rows for report in json.loads(row[PAPER_FIELD])]


def primary_ids(reports):
    return {report["source_id"] for report in reports if report.get("paper_handoff_state") == "primary"}


def read_table(path):
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def b_bundle(root, papers, policy="initial_bibliography_v1"):
    atlas, literature = root / "atlas.tsv", root / "literature.tsv"
    pd.DataFrame([{"gene_id": gene, "family_id": "FAM1", "orthogroup": "OG1",
                   "subfamily_id": "SUB1", "orthology_reference_gene": "REF1",
                   "domain_architecture": "kinase", "broad_function_bins": "signaling"}
                  for gene in ("G1", "G2")]).to_csv(atlas, sep="\t", index=False)
    pd.DataFrame([{"interval_gene_ids": "G1;G2", "reference_gene": "REF1",
                   "reference_function": "meristem receptor", "query_status": "retrieved",
                   "article_rank": str(index), "source_id": report["source_id"], "title": report["source_title"],
                   "claim_type": "literature_link_candidate"} for index, report in enumerate(papers, 1)]).to_csv(literature, sep="\t", index=False)
    return build_route_b_semantic_bundle(str(atlas), str(root / "bundle"),
        phenotype="fasciated meristem", literature_arg=str(literature), case_id="synthetic",
        b_review_membership_policy=policy)


class PaperHandoffTests(unittest.TestCase):
    def test_initial_tiers_choose_five_without_losing_twenty_or_changing_gene_tier(self):
        papers = twenty_papers()
        original = a_row(papers)
        untouched = deepcopy(original)
        stamped = stamp_paper_handoff([original], route="phenotype_first")[0]
        reports = json.loads(stamped[PAPER_FIELD])
        self.assertEqual(untouched, original)
        self.assertEqual(20, len(reports))
        self.assertEqual({f"PMID:{number}" for number in range(10016, 10021)}, primary_ids(reports))
        primary = sorted((report for report in reports if report["paper_handoff_state"] == "primary"),
                         key=lambda report: int(report["paper_handoff_rank"]))
        self.assertEqual(["T1", "T2", "T3", "T4", "T5"], [report["initial_tier"] for report in primary])
        self.assertEqual(["1", "2", "3", "4", "5"], [report["paper_handoff_rank"] for report in primary])
        self.assertEqual(normalize_paper_reports(papers), normalize_paper_reports(reports, discard_handoff=True))
        self.assertEqual(original["source_ids"], stamped["source_ids"])
        self.assertEqual(original["counter_hypothesis"], stamped["counter_hypothesis"])
        before = annotate_relevance([{"candidate_gene_id": "G1"}], [a.claim(papers)])[0]
        after = annotate_relevance([{"candidate_gene_id": "G1"}], [a.claim(reports)])[0]
        for key in ("relevance_tier", "support_state", "paper_tier_counts_json", "tier_source_id"):
            self.assertEqual(before[key], after[key])
        self.assertEqual(20, len(json.loads(after["paper_evidence_summary_json"])))

    def test_aliases_and_multiple_reports_share_one_paper_rank(self):
        papers = twenty_papers()
        papers[-1]["source_identifier_aliases"] = "DOI:10.1234/primary"
        papers += [dict(papers[-1], source_id="DOI:10.1234/primary", source_identifier_aliases="PMCID:PMC10020", reported_gene="SECOND_GENE"),
                   dict(papers[-1], source_id="PMCID:PMC10020", source_identifier_aliases="", reported_gene="THIRD_GENE")]
        row = stamp_paper_handoff([a_row(papers)], route="phenotype_first")[0]
        reports = json.loads(row[PAPER_FIELD])
        top = [report for report in reports if report["initial_tier"] == "T1"]
        self.assertEqual(3, len(top))
        self.assertEqual({"1"}, {report["paper_handoff_rank"] for report in top})
        self.assertEqual({"primary"}, {report["paper_handoff_state"] for report in top})
        selected = select_bibliography_papers(row, row["source_ids"].split(";"))
        self.assertEqual(20, selected["candidate_paper_count"])
        self.assertEqual(5, selected["selected_paper_count"])
        self.assertTrue(source_tokens([report["source_id"] for report in top]) <= source_tokens(selected["selected_source_ids"]))

    def test_five_top_tier_papers_keep_other_evidence_if_one_report_is_misgraded(self):
        papers = [tier_paper("T1" if number < 10007 else "T5", number) for number in range(10001, 10021)]
        row = stamp_paper_handoff([a_row(papers)], route="phenotype_first")[0]
        reports = json.loads(row[PAPER_FIELD])
        primary = [report for report in reports if report["paper_handoff_state"] == "primary"]
        self.assertEqual(["T1"] * 5, [report["initial_tier"] for report in primary])
        selected_ids = primary_ids(reports)
        # Simulate a later erroneous report copy, without rewriting a real freeze.
        # The other four T1 reports remain in the original five-paper handoff.
        primary[0]["paper_checklist"] = tier_paper("T5", 1)["paper_checklist"]
        primary[0]["initial_tier"] = "T5"
        row[PAPER_FIELD] = normalize_paper_reports(reports)
        selection = select_bibliography_papers(row, row["source_ids"].split(";"))
        self.assertEqual(selected_ids, set(selection["selected_source_ids"]))
        result = annotate_relevance([{"candidate_gene_id": "G1"}], [a.projected_claim(json.loads(row[PAPER_FIELD]))])[0]
        self.assertTrue(result["relevance_tier"].startswith("T1"))

    def test_reader_consumes_recorded_selection_without_reranking(self):
        row = stamp_paper_handoff([a_row(twenty_papers())], route="phenotype_first")[0]
        selected = select_bibliography_papers(row, row["source_ids"].split(";"))
        reports = json.loads(row[PAPER_FIELD])
        audit = next(report for report in reports if report["source_id"] == "PMID:10001")
        audit["paper_checklist"] = tier_paper("T1", 10001)["paper_checklist"]
        # This local mutation is only a reader test: a new import is required to
        # revise a real frozen generation. Metadata must not silently reselect.
        changed = select_bibliography_papers({**row, PAPER_FIELD: reports}, row["source_ids"].split(";"))
        self.assertEqual(PAPER_HANDOFF_POLICY, changed["policy"])
        self.assertEqual(selected["selected_source_ids"], changed["selected_source_ids"])
        self.assertIn("no_metadata_reselection", changed["selection_basis"])
        self.assertNotIn("PMID:10001", changed["selected_source_ids"])

    def test_import_discards_model_markers_and_each_reimport_recomputes(self):
        papers = twenty_papers()
        papers[0].update(paper_handoff_policy={"untrusted": True}, paper_handoff_scope="model-selected",
                         paper_handoff_state="primary", paper_handoff_rank=1)
        imported = imported_paper_fields(a_row(papers))
        self.assertTrue(all(not set(PAPER_HANDOFF_FIELDS) & report.keys() for report in json.loads(imported[PAPER_FIELD])))
        row = stamp_paper_handoff([{**a_row(papers), **imported}], route="phenotype_first")[0]
        self.assertNotIn("PMID:10001", primary_ids(json.loads(row[PAPER_FIELD])))
        self.assertEqual([row], stamp_paper_handoff([row], route="phenotype_first"))
        revised = json.loads(row[PAPER_FIELD])
        next(report for report in revised if report["source_id"] == "PMID:10001")["paper_checklist"] = tier_paper("T1", 10001)["paper_checklist"]
        again = stamp_paper_handoff([{**row, PAPER_FIELD: revised}], route="phenotype_first")[0]
        self.assertIn("PMID:10001", primary_ids(json.loads(again[PAPER_FIELD])))

    def test_each_original_panel_keeps_its_own_five_and_summary_retains_all_scopes(self):
        first = a_row(twenty_papers())
        second = a_row([tier_paper("T6", 10001)], hypothesis_id="h2")
        rows = stamp_paper_handoff([first, second], route="phenotype_first")
        self.assertEqual([5, 1], [len(primary_ids(json.loads(row[PAPER_FIELD]))) for row in rows])
        self.assertEqual(2, len({report["paper_handoff_scope"] for report in reports_in(rows)}))
        candidate = annotate_relevance([{"candidate_gene_id": "G1"}],
            [a.claim(json.loads(row[PAPER_FIELD]), hypothesis_id=row["hypothesis_id"]) for row in rows])[0]
        summary = next(report for report in json.loads(candidate["paper_evidence_summary_json"]) if report["source_id"] == "PMID:10001")
        self.assertEqual({"primary", "audit_only"}, {panel["paper_handoff_state"] for panel in summary["paper_handoff_panels"]})
        self.assertEqual(20, len(json.loads(candidate["paper_evidence_summary_json"])))
        replay = annotate_relevance([candidate])[0]
        self.assertEqual(candidate["paper_evidence_summary_json"], replay["paper_evidence_summary_json"])

    def test_old_observation_and_unmarked_selection_keep_their_representation(self):
        papers = [a.paper(f"PMID:{number}", tier="T2") for number in range(101, 108)]
        encoded = normalize_paper_reports(papers)
        self.assertTrue(all(not set(PAPER_HANDOFF_FIELDS) & report.keys() for report in json.loads(encoded)))
        original = {**a_row(papers), PAPER_FIELD: encoded}
        self.assertEqual([original], stamp_paper_handoff([original], route="phenotype_first"))
        self.assertEqual({"policy": "tier_relevant_papers_max5_v1", "max_papers": 5,
            "candidate_paper_count": 7, "selected_paper_count": 5,
            "selected_source_ids": [f"PMID:{number}" for number in range(101, 106)],
            "unselected_source_ids": ["PMID:106", "PMID:107"],
            "selection_basis": "original_paper_tier_then_stable_report_order; aliases_share_one_slot",
            "unselected_reason": "metadata_budget_exceeded; original observations retained",
            "membership_effect": "none; original paper reports retained"},
            select_bibliography_papers(original, original["source_ids"].split(";")))
        candidate = annotate_relevance([{"candidate_gene_id": "G1"}], [a.claim(papers)])[0]
        self.assertTrue(all("paper_handoff_panels" not in report for report in json.loads(candidate["paper_evidence_summary_json"])))

    def test_invalid_recorded_selection_is_rejected_without_inventing_another_panel(self):
        stamped = stamp_paper_handoff([a_row([tier_paper("T1", 101)])], route="phenotype_first")[0]
        first = json.loads(stamped[PAPER_FIELD])[0]
        malformed = []
        for field in PAPER_HANDOFF_FIELDS:
            missing = {**first, "paper_handoff_state": "audit_only", "paper_handoff_rank": "6"}
            del missing[field]
            malformed.append((f"missing {field}", [missing]))
        malformed.extend([
            ("same paper, two primary ranks", [first, {**first, "reported_gene": "SECOND_GENE", "paper_handoff_rank": "2"}]),
            ("same alias paper, conflicting state", [first, {**first, "source_id": "DOI:10.1234/alias", "source_identifier_aliases": "PMID:101",
                "paper_handoff_state": "audit_only", "paper_handoff_rank": "6"}]),
            ("audit rank inside primary slots", [{**first, "paper_handoff_state": "audit_only"}]),
            ("two papers, one primary rank", [first, {**first, "source_id": "PMID:102"}]),
        ])
        for reason, reports in malformed:
            with self.subTest(reason=reason), self.assertRaises(ValueError):
                select_bibliography_papers({PAPER_FIELD: reports}, [report["source_id"] for report in reports])

    def test_same_paper_rank_conflict_across_bound_claims_cannot_reach_summary(self):
        row = stamp_paper_handoff([a_row([tier_paper("T1", 101)])], route="phenotype_first")[0]
        first = json.loads(row[PAPER_FIELD])[0]
        second = {**first, "reported_gene": "SECOND_GENE", "paper_handoff_rank": "2"}
        # Each claim alone has a well-shaped stamp; the combined owner must
        # still prevent the same paper occupying two slots in one scope.
        claims = [a.claim([first], claim_id="c1"), a.claim([second], claim_id="c2")]
        with self.assertRaisesRegex(ValueError, "one paper has inconsistent"):
            annotate_relevance([{"candidate_gene_id": "G1"}], claims)

    def test_rejected_empty_audit_rank_can_coexist_with_the_same_selected_paper(self):
        row = stamp_paper_handoff([a_row([tier_paper("T1", 101)])], route="phenotype_first")[0]
        first = json.loads(row[PAPER_FIELD])[0]
        rejected = {**first, "reported_gene": "REJECTED_RECORD", "paper_handoff_state": "audit_only", "paper_handoff_rank": ""}
        papers = [first, rejected]
        selected = select_bibliography_papers({PAPER_FIELD: papers}, ["PMID:101"])
        self.assertEqual(1, selected["selected_paper_count"])
        candidate = annotate_relevance([{"candidate_gene_id": "G1"}], [a.claim(papers)])[0]
        summary = json.loads(candidate["paper_evidence_summary_json"])
        self.assertEqual(1, len(summary))
        self.assertEqual({"primary", "audit_only"}, {panel["paper_handoff_state"] for panel in summary[0]["paper_handoff_panels"]})

    def test_real_a_import_persists_twenty_reports_primary_five_and_readable_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            papers = twenty_papers()
            papers[0].update(paper_handoff_policy="model", paper_handoff_scope="model-scope", paper_handoff_state="primary", paper_handoff_rank="1")
            payload = {**a.MECHANISMS[0], **a_row(papers)}
            raw = a._payload(root, [payload])
            before = raw.read_bytes()
            self.assertEqual(0, main(["mechanism-import", "--payload", str(raw), "--out-dir", str(root / "import")]))
            rows = read_table(root / "import/mechanism_hypotheses.tsv")
            reports = reports_in(rows.to_dict("records"))
            self.assertEqual(20, len(reports))
            self.assertEqual(5, len(primary_ids(reports)))
            self.assertNotIn("PMID:10001", primary_ids(reports))
            table = read_table(root / "import/paper_evidence_table.tsv")
            self.assertEqual(20, len(table))
            self.assertEqual(5, table.selected_for_handoff.eq("yes").sum())
            self.assertEqual(table.selected_for_handoff.tolist(), table.selected_for_metadata.tolist())
            self.assertEqual({PAPER_HANDOFF_POLICY}, set(table.paper_handoff_policy))
            self.assertEqual(before, raw.read_bytes())
            self.assertEqual(before, (root / "import/raw_payload.json").read_bytes())

    def test_real_b_import_selects_once_across_source_expansion_and_keeps_rejected_paper(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            papers = twenty_papers()
            built = b_bundle(root, papers)
            claims = [{"claim_type": "direct_functional_bridge", "reference_gene": "REF1", "reference_species": "Plant beta",
                       "biological_function": "meristem receptor", "mechanism": "receptor signaling",
                       "phenotype_relation": "could explain fasciation", "evidence_scope": "reference_gene",
                       "source_ids": [report["source_id"] for report in group], PAPER_FIELD: group}
                      for group in (papers[:10], papers[10:])]
            rejected = tier_paper("T1", 99999)
            claims.append({**claims[0], "source_ids": [rejected["source_id"]], PAPER_FIELD: [rejected]})
            response = root / "response.json"
            b.RouteBSemanticImportTests()._response(Path(built["prompts"]), claims, response)
            frozen = {str(path): path.read_bytes() for path in (response, Path(built["bundle"]), Path(built["prompts"]))}
            result = import_route_b_semantic_proposals(built["bundle"], built["prompts"], str(response), str(root / "import"))
            rows = read_table(result["proposals"])
            reports = reports_in(rows.to_dict("records"))
            self.assertEqual(21, len(rows))
            self.assertEqual(20, result["proposed_claim_row_count"])
            self.assertEqual(1, result["unsupported_claim_row_count"])
            self.assertEqual({f"PMID:{number}" for number in range(10016, 10021)}, primary_ids(reports))
            self.assertEqual(1, len({report["paper_handoff_scope"] for report in reports}))
            retained = next(report for report in reports if report["source_id"] == "PMID:99999")
            self.assertEqual(("audit_only", "", "T1"), (retained["paper_handoff_state"], retained["paper_handoff_rank"], retained["initial_tier"]))
            self.assertEqual({"G1;G2"}, set(rows.member_interval_genes))
            table = read_table(result["paper_evidence_table"])
            self.assertEqual(5, table.selected_for_handoff.eq("yes").sum())
            for path, data in frozen.items():
                self.assertEqual(data, Path(path).read_bytes())

    def test_old_b_source_function_import_does_not_gain_handoff_fields(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            papers = [tier_paper("T2", 123)]
            built = b_bundle(root, papers, policy="source_function_v1")
            claim = {"claim_type": "mutant_phenotype", "reference_gene": "REF1", "evidence_scope": "reference_gene",
                     "source_ids": ["PMID:123"], PAPER_FIELD: papers}
            response = root / "response.json"
            b.RouteBSemanticImportTests()._response(Path(built["prompts"]), [claim], response)
            result = import_route_b_semantic_proposals(built["bundle"], built["prompts"], str(response), str(root / "import"))
            rows = read_table(result["proposals"])
            self.assertEqual(normalize_paper_reports(papers), rows.iloc[0][PAPER_FIELD])
            self.assertTrue(all(not set(PAPER_HANDOFF_FIELDS) & report.keys() for report in reports_in(rows.to_dict("records"))))
            self.assertNotIn("selected_for_handoff", read_table(result["paper_evidence_table"]).columns)


if __name__ == "__main__":
    unittest.main()
