"""Imported paper panels reach source audits and the real boundary request."""
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.boundary_tier_review import (
    EVIDENCE_SCOPED_POLICY, build_boundary_review_request, freeze_boundary_scope,
)
from candidate_qtl.cli import main
from candidate_qtl.paper_reports import (
    PAPER_FIELD, PAPER_HANDOFF_POLICY,
    normalize_paper_reports, stamp_paper_handoff,
)
from candidate_qtl.relevance_tiers import NEW_TIERS, annotate_relevance
from candidate_qtl.route_b_semantic_bundle import build_route_b_semantic_bundle
from candidate_qtl.route_b_semantic_candidates import freeze_route_b_semantic_candidates
from candidate_qtl.route_b_semantic_source_verify import _initial_bibliography_budget, verify_route_b_semantic_sources
from candidate_qtl.route_b_semantic_transfer_verify import verify_route_b_semantic_transfer
from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff
from scripts.agentic_boundary import prepare
from scripts.agentic_boundary_material import (
    expand_reading_view, paper_reading_projection, restore_reading_shard,
)
from tests.test_paper_checklist import report
from tests import test_paper_tiers as fixtures
from tests import test_route_b_semantic_import as b_fixtures


def paper(number, tier="T5", **changes):
    answers = {
        "T1": dict(phenotype_reported="yes", equivalent_target_phenotype="yes",
                   focal_species="yes", same_gene_as_focal="yes"),
        "T2": dict(phenotype_reported="yes", equivalent_target_phenotype="yes", focal_species="no"),
        "T4": dict(function_supports_hypothesis="yes"),
        "T5": dict(annotation_supports_hypothesis="yes"),
    }[tier]
    return {**report(**answers), "source_id": f"PMID:{number}", "source_title": f"Paper {number}",
            "reported_gene": "REF1", "reported_species": "Plant alpha",
            "source_excerpt": f"Original observation in paper {number}", **changes}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class PaperHandoffConsumerTests(unittest.TestCase):
    def test_a_and_b_metadata_consume_stamps_without_refilling_or_regrading(self):
        # A consumer obeys a frozen selection, even if re-sorting the displayed
        # Tiers would choose another paper. Only a new import selects a panel.
        reports = [paper(number, "T1" if number == 106 else "T5") for number in range(101, 107)]
        stamped = [{**value, "paper_handoff_policy": PAPER_HANDOFF_POLICY,
                    "paper_handoff_scope": "frozen-panel", "paper_handoff_rank": str(index),
                    "paper_handoff_state": "primary" if index <= 5 else "audit_only"}
                   for index, value in enumerate(reports, 1)]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cache = root / "cache.json"
            cache.write_text(json.dumps([{"id": str(i), "title": f"Paper {i}"} for i in range(101, 107)]))
            for name, values in (("marked", stamped), ("legacy", reports)):
                with self.subTest(policy=name):
                    raw = normalize_paper_reports(values)
                    source = root / f"{name}.tsv"
                    pd.DataFrame([{"hypothesis_id": "h", "reference_gene": "REF1",
                        "source_ids": ";".join(value["source_id"] for value in values), PAPER_FIELD: raw}]).to_csv(source, sep="\t", index=False)
                    original = source.read_bytes()
                    output = root / name
                    self.assertEqual(0, main(["source-verify", "--mechanisms", str(source),
                        "--out-dir", str(output), "--cache", str(cache)]))
                    selected = json.loads((output / "bibliography_selection.json").read_text())["panels"][0]
                    self.assertEqual(name == "legacy", "PMID:106" in selected["selected_source_ids"])
                    self.assertEqual(5, selected["selected_paper_count"])
                    self.assertEqual(PAPER_HANDOFF_POLICY if name == "marked" else "tier_relevant_papers_max5_v1", selected["policy"])
                    annotated = pd.read_csv(output / "mechanisms_with_verified_sources.tsv", sep="\t", dtype=str).fillna("").iloc[0]
                    self.assertEqual(raw, annotated[PAPER_FIELD])
                    self.assertEqual("not_evaluated", annotated.claim_support_summary)
                    self.assertEqual(original, source.read_bytes())
                    manifest = json.loads((output / "source_verification_manifest.json").read_text())
                    self.assertEqual(name == "marked", "handoff_consumer" in manifest["paper_selection"])

                    # The B adapter sees per-source rows, but consumes the same
                    # shared reference panel before any copy expansion.
                    rows = pd.DataFrame([{"claim_id": f"c{index}", "semantic_unit_id": "unit", "reference_gene": "REF1",
                        "reference_species": "Plant alpha", "source_id": value["source_id"],
                        PAPER_FIELD: normalize_paper_reports([value])} for index, value in enumerate(values)])
                    keep, skipped, panels = _initial_bibliography_budget(rows)
                    self.assertEqual(set(selected["selected_source_ids"]), set(keep.source_id))
                    self.assertEqual(1, len(skipped))
                    self.assertEqual(name == "marked", "audit_only_claim_ids" in panels[0])

    def test_boundary_payload_keeps_scope_specific_reports_and_shares_copies(self):
        reports = [paper(i, "T2" if i < 104 else "T4" if i < 106 else "T5") for i in range(101, 121)]
        proposals = [
            {"hypothesis_id": "panel-a", "reference_gene": "REF1", "source_ids": ";".join(p["source_id"] for p in reports), PAPER_FIELD: reports},
            # The same paper has a different reported gene in another scope.
            {"hypothesis_id": "panel-b", "reference_gene": "REF1", "source_ids": "PMID:101",
             PAPER_FIELD: [paper(101, "T4", reported_gene="OTHER_REPORTED_GENE", source_excerpt="A different experiment in the same paper")]},
        ]
        imported = stamp_paper_handoff(proposals, route="phenotype_first")
        scopes = [json.loads(item[PAPER_FIELD])[0]["paper_handoff_scope"] for item in imported]
        claims = [fixtures.claim(json.loads(item[PAPER_FIELD]), candidate_gene_id=gene, hypothesis_id=item["hypothesis_id"])
                  for gene in ("G1", "G2") for item in imported]
        rows = annotate_relevance([{"candidate_gene_id": gene, "support_status": "unverified_route_b_review"}
                                   for gene in ("G1", "G2")], claims)
        self.assertEqual({NEW_TIERS[1]}, {row["relevance_tier"] for row in rows})
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inventory = root / "inventory.tsv"
            pd.DataFrame(rows).to_csv(inventory, sep="\t", index=False)
            original = inventory.read_bytes()
            scope = freeze_boundary_scope(str(inventory), str(root / "scope"), inventory_sha256=sha(inventory),
                case_id="panel-case", target=1, upper=1, policy_version=EVIDENCE_SCOPED_POLICY,
                eligible_gene_ids=["G1", "G2"])
            frozen_bytes = {role: Path(value["path"]).read_bytes() for role, value in scope["artifacts"].items()}
            built = build_boundary_review_request(scope["manifest"], str(root / "request"),
                scope_sha256=scope["manifest_sha256"], case_id="panel-case")
            request = json.loads(Path(built["request"]).read_text())
            panels = {panel["paper_handoff_scope"]: panel for panel in request["primary_paper_evidence"]}
            self.assertEqual({"G1", "G2"}, set(request["candidate_gene_ids"]))
            self.assertEqual([5, 1], [len(panels[key]["papers"]) for key in scopes])
            self.assertEqual(["1", "2", "3", "4", "5"], [p["paper_handoff_rank"] for p in panels[scopes[0]]["papers"]])
            for panel in panels.values():
                self.assertEqual(["G1", "G2"], panel["candidate_gene_ids"])
                for value in panel["papers"]:
                    self.assertEqual("original_reports", value["record_state"])
                    self.assertEqual(1, len(value["reports"]))
                    self.assertEqual(["G1", "G2"], value["reports"][0]["candidate_gene_ids"])
            other_report = panels[scopes[1]]["papers"][0]["reports"][0]["paper_report"]
            self.assertEqual("OTHER_REPORTED_GENE", other_report["reported_gene"])
            self.assertEqual("A different experiment in the same paper", other_report["source_excerpt"])
            self.assertEqual("T4", other_report["initial_tier"])
            self.assertIn("Read primary_paper_evidence first", " ".join(request["policy"]))
            self.assertEqual(original, inventory.read_bytes())
            for role, value in scope["artifacts"].items():
                self.assertEqual(frozen_bytes[role], Path(value["path"]).read_bytes())
            audit = pd.read_csv(request["candidate_input"]["path"], sep="\t", dtype=str).fillna("")
            for value in audit.relevance_evidence_json:
                self.assertEqual(21, sum(len(json.loads(claim[PAPER_FIELD])) for claim in json.loads(value)))

    def test_unmarked_boundary_keeps_legacy_payload_and_missing_reports_are_explicit(self):
        from candidate_qtl.boundary_tier_review import _primary_paper_evidence
        row = annotate_relevance([{"candidate_gene_id": "G1"}], [fixtures.claim([paper(101, "T2")])])[0]
        self.assertEqual([], _primary_paper_evidence([row]))
        original = stamp_paper_handoff([{"hypothesis_id": "h", "reference_gene": "REF1",
            "source_ids": "PMID:101", PAPER_FIELD: [paper(101, "T2")]}], route="phenotype_first")[0]
        marked = annotate_relevance([{"candidate_gene_id": "G1"}], [fixtures.claim(json.loads(original[PAPER_FIELD]))])[0]
        marked["relevance_evidence_json"] = ""
        result = _primary_paper_evidence([marked])[0]["papers"][0]
        self.assertEqual("original_report_unavailable", result["record_state"])
        self.assertEqual([], result["reports"])
        self.assertEqual(["G1"], result["candidate_gene_ids"])

    def test_b_twenty_papers_reach_real_wet_boundary_workers_as_the_same_five(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = b_fixtures.RouteBSemanticImportTests()
            fixture._bundle(root)
            template = pd.read_csv(root / "atlas.tsv", sep="\t", dtype=str).fillna("").iloc[0].to_dict()
            genes = [f"G{index:02d}" for index in range(26)]
            atlas = root / "shared-copy-atlas.tsv"
            pd.DataFrame([{**template, "gene_id": gene} for gene in genes]).to_csv(atlas, sep="\t", index=False)
            reports = [paper(number, "T2" if number < 704 else "T4" if number < 706 else "T5",
                reported_gene=f"REPORTED_GENE_{number}", reported_gene_id=f"GeneID:{number}",
                reported_species="Plant beta", initial_tier_reason=f"Original reason for source {number}")
                for number in range(701, 721)]
            literature = root / "twenty-papers.tsv"
            pd.DataFrame([{"interval_gene_ids": ";".join(genes), "reference_gene": "REF1",
                "reference_function": "meristem receptor controls developmental signaling", "query_status": "retrieved",
                "article_rank": str(index), "source_id": value["source_id"], "title": value["source_title"],
                "claim_type": "literature_link_candidate"} for index, value in enumerate(reports, 1)]).to_csv(literature, sep="\t", index=False)
            bundle = build_route_b_semantic_bundle(str(atlas), str(root / "twenty-paper-bundle"),
                phenotype="fasciated meristem", literature_arg=str(literature), case_id="case403", focal_species="Plant alpha")
            response = root / "response.json"
            fixture._response(Path(bundle["prompts"]), [{"claim_type": "direct_functional_bridge",
                "biological_function": "meristem receptor controls developmental signaling", "mechanism": "receptor signaling",
                "phenotype_relation": "changed receptor signaling could explain fasciated meristem",
                "reference_gene": "REF1", "reference_species": "Plant beta", "evidence_scope": "reference_gene",
                "source_ids": [value["source_id"]], PAPER_FIELD: [value]} for value in reports], response)
            frozen_hashes = {path: sha(path) for path in (str(atlas), str(literature), bundle["bundle"], bundle["prompts"], str(response))}
            imported = root / "import"
            self.assertEqual(0, main(["route-b-semantic-import", "--bundle", bundle["bundle"], "--prompts", bundle["prompts"],
                "--responses", str(response), "--out-dir", str(imported)]))
            proposals = imported / "semantic_claim_proposals.tsv"
            frozen_hashes[str(proposals)] = sha(proposals)
            original = pd.read_csv(proposals, sep="\t", dtype=str).fillna("")
            stored = [value for encoded in original[PAPER_FIELD] for value in json.loads(encoded)]
            expected = {value["source_id"]: value for value in stored if value["paper_handoff_state"] == "primary"}
            self.assertEqual(20, len(stored))
            self.assertEqual({f"PMID:{number}" for number in range(701, 706)}, set(expected))
            scope_ids = {value["paper_handoff_scope"] for value in stored}
            self.assertEqual(1, len(scope_ids))
            cache = root / "cache.json"
            cache.write_text(json.dumps([{"id": str(number), "title": f"Paper {number}"} for number in range(701, 721)]))
            sources = verify_route_b_semantic_sources(str(proposals), str(root / "source"), bundle_arg=bundle["bundle"], cache_arg=str(cache))
            source_rows = pd.read_csv(sources["source_verification"], sep="\t", dtype=str).fillna("")
            self.assertEqual(15, source_rows.bibliography_state.eq("paper_handoff_audit_only").sum())
            selected = json.loads((root / "source/bibliography_selection.json").read_text())["panels"][0]
            self.assertEqual(set(expected), set(selected["selected_source_ids"]))
            transfer = verify_route_b_semantic_transfer(str(proposals), sources["source_verification"], bundle["bundle"],
                str(root / "transfer"), case_id="case403")
            frozen = freeze_route_b_semantic_candidates(transfer["claims"], str(root / "freeze"), case_id="case403")
            self.assertEqual(0, frozen["candidate_gene_count"])
            self.assertEqual(set(genes), set(pd.read_csv(frozen["review_candidates"], sep="\t").candidate_gene_id))
            frozen_hashes.update({path: sha(path) for path in (sources["source_verification"], transfer["claims"],
                frozen["candidates"], frozen["manifest"], frozen["review_candidates"], frozen["review_manifest"])})
            wet = wet_candidate_handoff(str(root / "wet"), route_b_candidates_arg=frozen["candidates"],
                route_b_candidates_manifest_arg=frozen["manifest"], review_candidates_arg=frozen["review_candidates"],
                review_candidates_manifest_arg=frozen["review_manifest"], interval_atlas_arg=str(atlas), case_id="case403")
            # All twenty observations remain in the existing hash-bound audit
            # ancestors. Unselected claims are not injected into a candidate's
            # applicable evidence, which would change its classification scope.
            wet_manifest = json.loads(Path(wet["wet_candidates_manifest"]).read_text())
            self.assertEqual(frozen["review_manifest"], wet_manifest["review_candidates_manifest"])
            producer = json.loads(Path(wet_manifest["review_candidates_manifest"]).read_text())
            audit_claims = producer["inputs"]["claims"]
            self.assertEqual(audit_claims["sha256"], sha(audit_claims["path"]))
            audit_rows = pd.read_csv(audit_claims["path"], sep="\t", dtype=str).fillna("")
            audit_reports = [value for encoded in audit_rows[PAPER_FIELD] for value in json.loads(encoded)]
            self.assertEqual({value["source_id"]: value for value in stored},
                             {value["source_id"]: value for value in audit_reports})
            transfer_descriptor = producer["inputs"]["claims_manifest"]
            self.assertEqual(transfer_descriptor["sha256"], sha(transfer_descriptor["path"]))
            transfer_manifest = json.loads(Path(transfer_descriptor["path"]).read_text())
            proposals_descriptor = transfer_manifest["ancestors"]["proposals"]
            self.assertEqual(str(proposals), proposals_descriptor["path"])
            self.assertEqual(frozen_hashes[str(proposals)], proposals_descriptor["sha256"])
            request_path = Path(wet["adaptive_synthesis_request"])
            request = json.loads(request_path.read_text())
            self.assertEqual(set(genes), set(request["candidate_gene_ids"]))
            self.assertEqual(1, len(request["primary_paper_evidence"]))
            frozen_hashes[str(request_path)] = sha(request_path)
            for descriptor in (request["candidate_input"], request["higher_tier_context"]["candidate_input"]):
                frozen_hashes[descriptor["path"]] = sha(descriptor["path"])
            prepared = root / "prepared"
            manifest = prepare(str(request_path), str(prepared), max_candidates=7, max_material_bytes=10_000_000)
            seen = set()
            reading_views = {value["shard_id"]: value for value in manifest["reading_views"]}
            catalog = json.loads((prepared / manifest["material_catalog"]["file"]).read_text())["materials"]
            for shard in manifest["shards"]:
                assigned = set(shard["candidate_gene_ids"])
                seen.update(assigned)
                canonical = json.loads((prepared / shard["file"]).read_text())
                view = json.loads((prepared / reading_views[shard["shard_id"]]["file"]).read_text())
                # Validate hashes, page delivery and the complete canonical
                # projection before resolving the worker's shared catalog refs.
                restored = restore_reading_shard(prepared, manifest, shard["shard_id"])
                self.assertEqual(canonical, restored)
                view = expand_reading_view(view, catalog)
                self.assertEqual(paper_reading_projection(canonical["primary_paper_evidence"], structured_json=True),
                                 view["primary_paper_evidence"])
                raw_received = {bound["paper_report"]["source_id"]: bound["paper_report"]
                                for value in canonical["primary_paper_evidence"][0]["papers"]
                                for bound in value["reports"]}
                self.assertEqual(expected, raw_received)
                panel = view["primary_paper_evidence"][0]
                self.assertEqual(scope_ids, {panel["paper_handoff_scope"]})
                self.assertEqual(assigned, set(panel["candidate_gene_ids"]))
                self.assertEqual(5, len(panel["papers"]))
                self.assertEqual(["1", "2", "3", "4", "5"], [value["paper_handoff_rank"] for value in panel["papers"]])
                received = {}
                for value in panel["papers"]:
                    self.assertEqual(assigned, set(value["candidate_gene_ids"]))
                    self.assertEqual(1, len(value["reports"]))
                    bound = value["reports"][0]
                    self.assertEqual(assigned, set(bound["candidate_gene_ids"]))
                    received[bound["paper_report"]["source_id"]] = bound["paper_report"]
                # Original Tier, reason, reported gene, reference and answers,
                # including losslessly structured *_json values in the view.
                self.assertEqual(paper_reading_projection(expected, structured_json=True), received)
                self.assertEqual(canonical["candidate_input"], restored["candidate_input"])
                for row in restored["candidate_input"]:
                    all_reports = [value for claim in json.loads(row["relevance_evidence_json"])
                                   for value in json.loads(claim.get(PAPER_FIELD) or "[]")]
                    self.assertEqual(set(expected), {value["source_id"] for value in all_reports})
                    self.assertEqual(NEW_TIERS[1], row["relevance_tier"])
                    self.assertEqual("review", row["support_state"])
            self.assertEqual(set(genes), seen)
            self.assertEqual(4, len(manifest["shards"]))
            self.assertEqual(frozen_hashes, {path: sha(path) for path in frozen_hashes})
