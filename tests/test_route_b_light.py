"""Light Route B research path (#1004): prepare, machine checks, conversion."""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
import urllib.error
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from unittest import mock

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.route_b_semantic_bundle import refreeze_route_b_sources
from candidate_qtl.route_b_semantic_candidates import freeze_route_b_semantic_candidates
from candidate_qtl.route_b_semantic_import import import_route_b_semantic_proposals
from candidate_qtl.route_b_semantic_source_verify import verify_route_b_semantic_sources
from candidate_qtl.route_b_semantic_transfer_verify import verify_route_b_semantic_transfer
from scripts import agentic_route_b_light as light
from scripts.agentic_route_b_light import (
    CHECKLIST_FIELDS,
    MAX_PAPERS,
    light_convert,
    light_prepare,
)

ROOT = Path(__file__).resolve().parents[1]
PAPER = "PMID:11111111"
PAPER_TITLE = "A homolog controls leaflet development"
PAPER_ABSTRACT = "Loss of the homolog produced fused leaflets in Arabidopsis thaliana."


def digest(data: bytes) -> str:
    import hashlib

    return hashlib.sha256(data).hexdigest()


def unit(uid: str, genes: list[tuple[str, str]], *, requested: bool = True,
         gene_fields: dict[str, dict] | None = None,
         literature: list[dict] | None = None) -> dict:
    records, cards = [], []
    for gid, reference in genes:
        gene = {"gene_id": gid, "reference_gene": reference, "reference_species": "Arabidopsis thaliana",
                "reference_function": f"{reference} leaf development process",
                "protein_name": f"{reference} protein", "go_terms": "leaf development"}
        gene.update((gene_fields or {}).get(gid, {}))
        records.append(gene)
        cards.append({"interval_gene_id": gid, "reference_gene": reference,
                      "reference_species": "Arabidopsis thaliana",
                      "function_description": f"{reference} controls leaflet development",
                      "mutant_phenotype": f"{reference} mutant has fused leaflets"})
    return {
        "unit_id": uid, "unit_type": "family", "unit_anchor": genes[0][1],
        "member_interval_genes": [gid for gid, _ in genes],
        "identity": {"evolutionary_family_ids": [f"EV-{uid}"], "family_ids": [f"F-{uid}"],
                     "orthogroup_ids": [f"OG-{uid}"], "subfamily_ids": [f"SUB-{uid}"]},
        "genes": records, "reference_material_cards": cards,
        "literature_link_candidates": literature or [], "source_provenance": [],
        "stable_source_ids": [], "source_backed_function_record_count": 1,
        "request_state": "requested" if requested else "not_requested_no_reviewable_material",
        "review_status": "not_evaluated", "source_material_class": "function_context",
        "estimated_review_cost": {"estimated_input_tokens": 10, "estimated_model_call_count": 1},
    }


def fixture() -> dict:
    return {
        "case_id": "fixture-light", "route": "interval_first", "stage": "post_b0_semantic_bundle",
        "phenotype": "fused leaflets", "focal_species": "Pisum sativum",
        "b_review_membership_policy": "initial_bibliography_v1",
        "interval_gene_count": 7, "semantic_unit_count": 7,
        "units": [
            unit("u1", [("g1", "REF-A")], gene_fields={"g1": {
                "orthology_reference_gene": "orth-1",
                "family_index_subject_endpoints": "fam-1",
                "sequence_similarity_hits": [{"hit_accession": "hit-1", "hit_species": "Medicago truncatula"}],
            }}, literature=[{"reference_gene": "lit-1", "reference_species": "Medicago truncatula",
                             "source_id": "PMID:42"}]),
            unit("u2", [("g2", "REF-B")]),
            unit("u3", [("g3", "REF-C")]),
            unit("u4", [("g4", "REF-D")]),
            unit("u5", [("g5", "REF-E")]),
            unit("u6", [("g6", "REF-F")]),
            unit("u7", [("g7", "REF-G")]),
        ],
    }


def checklist(tier: str) -> dict[str, str]:
    checks = {key: "unknown" for key in CHECKLIST_FIELDS}
    if tier == "T2":
        checks.update(phenotype_reported="yes", equivalent_target_phenotype="yes",
                      focal_species="no", identity_resolved="yes")
    elif tier == "T3":
        checks.update(phenotype_reported="yes", shared_specific_process="yes", identity_resolved="yes")
    return checks


def cache_entry(paper_id: str, title: str) -> dict:
    kind, identifier = paper_id.split(":", 1)
    record = {"pmid": "", "pmcid": "", "doi": "", "title": title, "abstract": PAPER_ABSTRACT,
              "url": "https://doi.org/10.1000/leaflet", "retrieved_at": "2026-01-01T00:00:00+00:00",
              "journal": "Plant Journal", "year": "2024", "authors": "Example A"}
    record[kind.casefold()] = identifier
    record["aliases"] = [f"{kind.casefold()}:{identifier}"]
    return record


class LightPrepareTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.bundle = self.root / "route_b_semantic_bundle.json"
        self.bundle.write_text(json.dumps(fixture()), encoding="utf-8")
        self.sha = digest(self.bundle.read_bytes())

    def prepare(self, **kwargs):
        return light_prepare(str(self.bundle), self.sha, str(self.root / "light"), **kwargs)

    def index(self) -> dict:
        return json.loads((self.root / "light/light/index.json").read_text(encoding="utf-8"))

    def test_requested_units_are_packed_within_the_limit(self):
        result = self.prepare(batch_units=2)
        self.assertEqual(7, result["requested_unit_count"])
        batches = self.index()["batches"]
        self.assertEqual([2, 2, 2, 1], [len(batch["unit_ids"]) for batch in batches])
        self.assertEqual([["u1", "u2"], ["u3", "u4"], ["u5", "u6"], ["u7"]],
                         [batch["unit_ids"] for batch in batches])
        self.assertEqual(["g1", "g2"], batches[0]["gene_ids"])

    def test_units_file_restricts_to_the_listed_subset_in_canonical_order(self):
        listed = self.root / "units.txt"
        listed.write_text("u4\nu2\n", encoding="utf-8")
        self.prepare(batch_units=20, units_file=str(listed))
        batches = self.index()["batches"]
        self.assertEqual([["u2", "u4"]], [batch["unit_ids"] for batch in batches])
        self.assertEqual(2, self.index()["requested_unit_count"])

    def test_units_file_rejects_unknown_or_unrequested_units(self):
        listed = self.root / "units.txt"
        listed.write_text("u2\nnope\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "not requested"):
            self.prepare(units_file=str(listed))
        self.bundle.write_text(json.dumps(fixture()), encoding="utf-8")
        value = json.loads(self.bundle.read_text())
        value["units"][6]["request_state"] = "not_requested_no_reviewable_material"
        self.bundle.write_text(json.dumps(value), encoding="utf-8")
        self.sha = digest(self.bundle.read_bytes())
        listed.write_text("u7\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "not requested"):
            self.prepare(units_file=str(listed))

    def test_material_contains_each_requested_gene_and_its_reference_cards(self):
        self.prepare(batch_units=20)
        material = (self.root / "light/light/b01/material.md").read_text(encoding="utf-8")
        self.assertIn("## g1", material)
        self.assertIn("- unit_id: u1", material)
        self.assertIn("ref REF-A [Arabidopsis thaliana]", material)
        self.assertIn("function_description: REF-A controls leaflet development", material)
        self.assertIn("mutant_phenotype: REF-A mutant has fused leaflets", material)
        self.assertIn("- gene_reference_gene: REF-A", material)
        self.assertIn("- gene_orthology_reference_gene: orth-1", material)
        for gene in ("g2", "g6"):
            self.assertIn(f"## {gene}", material)
        task = (self.root / "light/light/b01/task.md").read_text(encoding="utf-8")
        self.assertIn("fused leaflets", task)
        self.assertIn("`not_assessed` means the gene was not judged", task)
        self.assertIn("`reference_gene`: the exact `ref` accession", task)
        self.assertIn("`results.jsonl` in the batch directory", task)
        self.assertIn("up to about 20", task)
        self.assertIn("up to 5 useful papers", task)
        self.assertIn("homologs in other plants", task)

    def test_roster_endpoints_match_the_importer_grounding_sources(self):
        self.prepare(batch_units=20)
        roster = json.loads((self.root / "light/light/b01/roster.json").read_text(encoding="utf-8"))
        self.assertEqual({"REF-A", "orth-1", "fam-1", "hit-1", "lit-1"}, set(roster["g1"]["reference_endpoints"]))
        self.assertEqual("u1", roster["g1"]["unit_id"])
        self.assertIn("REF-A", roster["g1"]["names"])
        self.assertEqual({"REF-B"}, set(roster["g2"]["reference_endpoints"]))

    def test_material_preserves_original_family_and_conflict_context_per_gene(self):
        value = fixture()
        context = {
            "g1": {"family_id": "F-A", "family_label": "Family A", "reference_family_id": "RF-A",
                   "evidence_completeness": "conflicting", "notes": "Function agrees; family assignment conflicts."},
            "g2": {"family_id": "F-B", "family_label": "Family B", "reference_family_id": "RF-B",
                   "evidence_completeness": "family_assignment", "notes": "Original copy B annotation."},
        }
        value["units"] = [unit("u-context", [("g1", "REF-A"), ("g2", "REF-A")], gene_fields=context)]
        value["interval_gene_count"], value["semantic_unit_count"] = 2, 1
        self.bundle.write_text(json.dumps(value), encoding="utf-8")
        self.sha = digest(self.bundle.read_bytes())
        self.prepare(batch_units=20)
        directory = self.root / "light/light/b01"
        material = (directory / "material.md").read_text(encoding="utf-8")
        self.assertIn("# Original unit identity: u-context\n" + json.dumps(
            value["units"][0]["identity"], sort_keys=True, ensure_ascii=False), material)
        sections = material.split("## g1\n", 1)[1].split("## g2\n", 1)
        for section, gene_id, other_id in ((sections[0], "g1", "g2"), (sections[1], "g2", "g1")):
            for field, original in context[gene_id].items():
                self.assertIn(f"- gene {field}: {original}\n", section)
            self.assertNotIn(context[other_id]["notes"], section)
        full = json.loads((directory / "material.full.json").read_text(encoding="utf-8"))
        self.assertEqual(value["units"], full["units"])
        self.assertEqual(digest((directory / "material.md").read_bytes()), self.index()["batches"][0]["material_sha256"])

    def test_gene_context_does_not_change_shared_cards_or_copy_bindings(self):
        item = unit("u-copies", [("g1", "REF-A"), ("g2", "REF-A"), ("g3", "REF-B")],
                    gene_fields={"g1": {"family_id": "F-A", "notes": "conflict for copy A"},
                                 "g2": {"family_id": "F-B", "notes": "copy B only"}})
        material = light._render_material("b01", [item])
        self.assertEqual(1, material.count("## card-1\n"))
        self.assertEqual(1, material.count("## card-2\n"))
        sections = material.split("# Shared biological reference cards", 1)[0]
        for gene_id, next_heading, expected, excluded in (
            ("g1", "## g2", "card-1", "card-2"),
            ("g2", "## g3", "card-1", "card-2"),
            ("g3", "# Shared", "card-2", "card-1"),
        ):
            section = sections.split(f"## {gene_id}\n", 1)[1].split(next_heading, 1)[0]
            self.assertIn(f"- complete biological card: {expected}\n", section)
            self.assertNotIn(f"- complete biological card: {excluded}\n", section)

    def test_legacy_missing_context_does_not_fabricate_family_or_identity(self):
        item = unit("u-old", [("g1", "REF-A")])
        item.pop("identity")
        material = light._render_material("b01", [item])
        self.assertIn("## g1\n- unit_id: u-old\n", material)
        self.assertIn("function_description: REF-A controls leaflet development", material)
        for absent in ("Original unit identity", "gene family_id:", "gene family_label:",
                       "gene reference_family_id:", "gene evidence_completeness:", "gene notes:"):
            self.assertNotIn(absent, material)

    def test_index_hashes_describe_the_written_batch_files(self):
        self.prepare(batch_units=3)
        batch = self.index()["batches"][0]
        directory = self.root / "light/light" / batch["bid"]
        self.assertEqual(digest((directory / "material.md").read_bytes()), batch["material_sha256"])
        self.assertEqual(digest((directory / "task.md").read_bytes()), batch["task_sha256"])
        self.assertEqual(digest((directory / "roster.json").read_bytes()), batch["roster_sha256"])
        self.assertEqual("agentic-route-b-light-1.0.0", self.index()["contract"])
        self.assertEqual(self.sha, self.index()["bundle"]["sha256"])

    def test_bundle_hash_and_stage_are_verified(self):
        with self.assertRaisesRegex(ValueError, "SHA-256"):
            light_prepare(str(self.bundle), "0" * 64, str(self.root / "bad"))
        self.assertFalse((self.root / "bad").exists())
        value = json.loads(self.bundle.read_text())
        value["stage"] = "other"
        self.bundle.write_text(json.dumps(value), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "post-B0"):
            light_prepare(str(self.bundle), digest(self.bundle.read_bytes()), str(self.root / "bad-stage"))

    def test_cli_prepare_writes_the_same_index(self):
        result = subprocess.run(
            [sys.executable, str(ROOT / "scripts/agentic_route_b.py"), "light-prepare",
             "--bundle", str(self.bundle), "--bundle-sha256", self.sha,
             "--out", str(self.root / "cli"), "--batch-units", "20"],
            capture_output=True, text=True, encoding="utf-8")
        self.assertEqual(0, result.returncode, result.stderr)
        summary = json.loads(result.stdout)
        self.assertEqual(7, summary["requested_unit_count"])
        self.assertTrue((self.root / "cli/light/b01/material.md").is_file())


class LightConvertTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.bundle = self.root / "route_b_semantic_bundle.json"
        self.bundle.write_text(json.dumps(fixture()), encoding="utf-8")
        self.sha = digest(self.bundle.read_bytes())
        self.units = self.root / "units.txt"
        self.units.write_text("u1\nu2\nu3\nu4\nu5\nu6\n", encoding="utf-8")
        light_prepare(str(self.bundle), self.sha, str(self.root / "light"),
                      batch_units=20, units_file=str(self.units))
        self.index = self.root / "light/light/index.json"
        self.index_sha = digest(self.index.read_bytes())
        self.cache = self.root / "cache.json"
        entry = cache_entry(PAPER, PAPER_TITLE)
        other = cache_entry("DOI:10.1000/other", "A different study of root architecture")
        self.cache.write_text(json.dumps({PAPER.casefold(): entry, "doi:10.1000/other": other}), encoding="utf-8")
        self.write_results()

    def write_results(self, rows: list[dict] | None = None) -> None:
        rows = rows if rows is not None else [
            {"interval_gene_id": "g1", "reference_gene": "ref-a", "tier": "T2",
             "function": "leaflet boundary specification", "reason": "Homolog loss fuses leaflets, matching the target.",
             "papers": [{"paper_id": PAPER, "title": PAPER_TITLE, "reported_gene": "HOM1",
                         "reported_species": "Arabidopsis thaliana", "gene_correspondence": "homolog",
                         "finding": "HOM1 loss fused leaflets in Arabidopsis.",
                         "checklist": checklist("T2")}], "literature": "screened 3"},
            {"interval_gene_id": "g2", "reference_gene": "WRONG-REF", "tier": "T2",
             "function": "cell wall process", "reason": "A bridge that relies on an unfrozen reference.", "papers": []},
            {"interval_gene_id": "g3", "reference_gene": "REF-C", "tier": "T3",
             "function": "vascular patterning", "reason": "Shared process claim.",
             "papers": [{"paper_id": "PMID:99999999", "title": "Unknown paper", "reported_gene": "X",
                         "reported_species": "Arabidopsis thaliana", "gene_correspondence": "homolog",
                         "finding": "unverified", "checklist": checklist("T3")},
                        {"paper_id": "DOI:10.1000/other", "title": "A study of leaflet development",
                         "reported_gene": "Y", "reported_species": "Arabidopsis thaliana",
                         "gene_correspondence": "homolog", "finding": "mismatched title",
                         "checklist": checklist("T3")}]},
            {"interval_gene_id": "g4", "reference_gene": "REF-D", "tier": "T6",
             "function": "identity only", "reason": "Identity resolved but no functional bridge.", "papers": []},
            {"interval_gene_id": "g6", "reference_gene": "", "tier": "not_assessed",
             "function": "", "reason": "", "papers": []},
            {"interval_gene_id": "gX", "reference_gene": "REF-A", "tier": "T1",
             "function": "invented", "reason": "A gene outside the batch.", "papers": []},
        ]
        target = self.root / "light/light/b01/results.jsonl"
        target.write_text("".join(json.dumps(row) + "\n" for row in rows), encoding="utf-8")

    def convert(self, **kwargs):
        return light_convert(str(self.bundle), self.sha, str(self.index), self.index_sha,
                             str(self.root / "converted"), cache_arg=str(self.cache),
                             offline=True, **kwargs)

    def responses(self) -> dict[str, dict]:
        payload = json.loads((self.root / "converted/agent_responses.json").read_text(encoding="utf-8"))
        return {row["prompt_id"]: row for row in payload["responses"]}

    def test_cli_convert_writes_the_canonical_artifacts(self):
        result = subprocess.run(
            [sys.executable, str(ROOT / "scripts/agentic_route_b.py"), "light-convert",
             "--bundle", str(self.bundle), "--bundle-sha256", self.sha,
             "--light-index", str(self.index), "--light-index-sha256", self.index_sha,
             "--out", str(self.root / "cli-converted"), "--cache", str(self.cache), "--offline"],
            capture_output=True, text=True, encoding="utf-8")
        self.assertEqual(0, result.returncode, result.stderr)
        self.assertTrue((self.root / "cli-converted/agent_responses.json").is_file())
        self.assertTrue((self.root / "cli-converted/source_verify_cache.json").is_file())
        self.assertTrue((self.root / "cli-converted/light_convert_receipt.json").is_file())

    def test_machine_checks_reject_bad_gene_reference_and_papers(self):
        self.convert()
        check = json.loads((self.root / "converted/light_check.json").read_text(encoding="utf-8"))
        gene_issues = {(row["gene_id"], row["issue"]) for row in check["gene_issues"]}
        self.assertIn(("gX", "gene_not_in_batch"), gene_issues)
        self.assertIn(("g2", "reference_not_in_frozen_endpoints"), gene_issues)
        paper_issues = {(row["paper_id"], row["state"]) for row in check["paper_issues"]}
        self.assertIn(("PMID:99999999", "not_attempted"), paper_issues)
        self.assertIn(("DOI:10.1000/other", "title_mismatch"), paper_issues)
        self.assertEqual(1, check["totals"]["invalid_gene_id"])
        self.assertEqual(1, check["totals"]["invalid_reference"])
        self.assertEqual(1, check["totals"]["papers_checked"])
        self.assertEqual("light-route-b-convert", json.loads(
            (self.root / "converted/agent_responses.json").read_text(encoding="utf-8"))["executor"])

    def test_claim_keeps_the_frozen_reference_spelling_and_drops_failed_papers(self):
        self.convert()
        responses = self.responses()
        claims = json.loads(responses["route-b-semantic-u1"]["raw_text"])["claims"]
        self.assertEqual(1, len(claims))
        self.assertEqual("REF-A", claims[0]["reference_gene"])
        self.assertEqual("literature_bridge", claims[0]["claim_type"])
        self.assertEqual(["HOM1"], [report["reported_gene"] for report in claims[0]["paper_reports_json"]])
        self.assertEqual([PAPER], claims[0]["source_ids"])
        self.assertEqual(9, len(claims[0]["paper_reports_json"][0]["paper_checklist"]))
        self.assertEqual("T2", claims[0]["initial_tier"])

    def test_units_stay_not_evaluated_instead_of_no_claims(self):
        self.convert()
        responses = self.responses()
        self.assertEqual("not_evaluated", responses["route-b-semantic-u2"]["status"])
        self.assertEqual("not_evaluated", responses["route-b-semantic-u3"]["status"])
        self.assertEqual([], json.loads(responses["route-b-semantic-u4"]["raw_text"])["claims"])
        self.assertEqual("not_evaluated", responses["route-b-semantic-u5"]["status"])
        self.assertEqual("not_evaluated", responses["route-b-semantic-u6"]["status"])
        self.assertEqual("outside light run", responses["route-b-semantic-u7"]["reason"])

    def test_discoveries_only_contain_new_ids_and_pass_refreeze(self):
        result = self.convert()
        discoveries = json.loads((self.root / "converted/discoveries.json").read_text(encoding="utf-8"))
        self.assertEqual(self.sha, discoveries["parent_bundle_sha256"])
        self.assertEqual(1, len(discoveries["discovered_sources"]))
        entry = discoveries["discovered_sources"][0]
        self.assertEqual(("u1", "REF-A", PAPER), (entry["unit_id"], entry["reference_gene"], entry["source_id"]))
        self.assertEqual("Arabidopsis thaliana", entry["reference_species"])
        self.assertEqual(PAPER_ABSTRACT, entry["source_text"])
        extended = refreeze_route_b_sources(str(self.bundle), str(self.root / "converted/discoveries.json"),
                                            str(self.root / "extended"), parent_sha256=self.sha)
        self.assertEqual(["u1"], extended["changed_unit_ids"])
        self.assertEqual(str(self.root / "converted/discoveries.json"), result["discoveries"])

    def test_both_cache_formats_are_written_from_the_same_metadata(self):
        self.convert()
        metadata = json.loads((self.root / "converted/metadata_cache.json").read_text(encoding="utf-8"))
        self.assertEqual(PAPER_TITLE, metadata[PAPER.casefold()]["title"])
        self.assertEqual(PAPER_ABSTRACT, metadata[PAPER.casefold()]["abstract"])
        source = json.loads((self.root / "converted/source_verify_cache.json").read_text(encoding="utf-8"))
        record = next(row for row in source["records"] if row["title"] == PAPER_TITLE)
        self.assertEqual(PAPER_TITLE, record["title"])
        self.assertEqual(PAPER_ABSTRACT, record["abstractText"])
        self.assertTrue(record["cache_key"].startswith("pmid_"))

    def test_receipt_lists_the_canonical_next_commands(self):
        self.convert()
        receipt = json.loads((self.root / "converted/light_convert_receipt.json").read_text(encoding="utf-8"))
        steps = {row["step"]: row["command"] for row in receipt["next_commands"]}
        self.assertIn("route-b-semantic-bundle", steps["extend"])
        self.assertIn("--research-discoveries", steps["extend"])
        self.assertIn("llm-run --provider codex-agent --model light", steps["llm-run"])
        self.assertIn("agent_responses.json", steps["llm-run"])
        self.assertIn("route-b-semantic-source-verify", steps["source-verify"])
        self.assertIn("source_verify_cache.json", steps["source-verify"])
        self.assertIn("route-b-semantic-candidates", steps["freeze"])
        self.assertEqual(1, receipt["counts"]["claims"])

    def test_online_mode_fetches_and_extends_the_cache_on_a_resolve(self):
        self.write_results([{"interval_gene_id": "g1", "reference_gene": "REF-A", "tier": "T2",
                             "function": "leaflet boundary specification",
                             "reason": "A fetched homolog paper reports fused leaflets.",
                             "papers": [{"paper_id": "PMID:22222222", "title": "A fetched homolog study",
                                         "reported_gene": "HOM2",
                                         "reported_species": "Arabidopsis thaliana",
                                         "gene_correspondence": "homolog", "finding": "Fused leaflets.",
                                         "checklist": checklist("T2")}]}])
        fetch = lambda kind, identifier, timeout=30: cache_entry(f"{kind}:{identifier}", "A fetched homolog study")  # noqa: E731
        result = light_convert(str(self.bundle), self.sha, str(self.index), self.index_sha,
                               str(self.root / "online"), offline=False, fetch=fetch)
        self.assertEqual(1, result["counts"]["claims"])
        metadata = json.loads((self.root / "online/metadata_cache.json").read_text(encoding="utf-8"))
        self.assertEqual("A fetched homolog study", metadata["pmid:22222222"]["title"])

    def test_exact_duplicates_do_not_consume_distinct_paper_slots(self):
        papers = [{"paper_id": PAPER, "title": PAPER_TITLE, "reported_gene": "HOM1",
                   "reported_species": "Arabidopsis thaliana", "gene_correspondence": "homolog",
                   "finding": "Fused leaflets.", "checklist": checklist("T2")}] * (MAX_PAPERS + 1)
        self.write_results([{"interval_gene_id": "g1", "reference_gene": "REF-A", "tier": "T2",
                             "function": "leaflet boundary specification", "reason": "Homolog loss fuses leaflets.",
                             "papers": papers}])
        self.convert()
        check = json.loads((self.root / "converted/light_check.json").read_text(encoding="utf-8"))
        self.assertNotIn("more_than_five_papers_truncated",
                         {row["issue"] for row in check["gene_issues"]})
        self.assertEqual(MAX_PAPERS, len(check["paper_issues"]))

    def test_minimal_four_field_cache_entry_is_still_resolvable(self):
        paper = "PMID:33333333"
        cache = self.root / "minimal-cache.json"
        cache.write_text(json.dumps({"pmid:33333333": {
            "title": "Minimal cached paper", "abstract": "A cached abstract.",
            "url": "https://pubmed.ncbi.nlm.nih.gov/33333333/",
            "retrieved_at": "2026-01-01T00:00:00+00:00"}}), encoding="utf-8")
        self.write_results([{"interval_gene_id": "g1", "reference_gene": "REF-A", "tier": "T2",
                             "function": "leaflet boundary specification",
                             "reason": "Cached homolog loss fuses leaflets.",
                             "papers": [{"paper_id": paper, "title": "Minimal cached paper",
                                         "reported_gene": "HOM3",
                                         "reported_species": "Arabidopsis thaliana",
                                         "gene_correspondence": "homolog", "finding": "Fused leaflets.",
                                         "checklist": checklist("T2")}]}])
        light_convert(str(self.bundle), self.sha, str(self.index), self.index_sha,
                      str(self.root / "minimal"), cache_arg=str(cache), offline=True)
        source = json.loads((self.root / "minimal/source_verify_cache.json").read_text(encoding="utf-8"))
        record = next(row for row in source["records"] if row["pmid"] == "33333333")
        self.assertEqual("pmid_33333333", record["cache_key"])
        self.assertEqual("Minimal cached paper", record["title"])

    def test_no_new_paper_writes_no_discoveries_and_says_so(self):
        value = fixture()
        value["units"][0]["stable_source_ids"] = [PAPER]
        bundle = self.root / "already-frozen.json"
        bundle.write_text(json.dumps(value), encoding="utf-8")
        sha = digest(bundle.read_bytes())
        light_prepare(str(bundle), sha, str(self.root / "frozen-light"), batch_units=20,
                      units_file=str(self.units))
        index = self.root / "frozen-light/light/index.json"
        (self.root / "frozen-light/light/b01/results.jsonl").write_text(
            (self.root / "light/light/b01/results.jsonl").read_text(encoding="utf-8"), encoding="utf-8")
        light_convert(str(bundle), sha, str(index), digest(index.read_bytes()),
                      str(self.root / "frozen-converted"), cache_arg=str(self.cache), offline=True)
        self.assertFalse((self.root / "frozen-converted/discoveries.json").exists())
        receipt = json.loads((self.root / "frozen-converted/light_convert_receipt.json").read_text(encoding="utf-8"))
        self.assertTrue(receipt["discovery_state"].startswith("none_required"))
        self.assertIsNone(receipt["outputs"]["discoveries"])
        self.assertNotIn("extend", {row["step"] for row in receipt["next_commands"]})


class LightAgentFileTests(unittest.TestCase):
    def test_not_evaluated_declaration_survives_llm_run(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        prompts = root / "prompts.jsonl"
        prompts.write_text(json.dumps({"prompt_id": "p1", "gene_id": "", "system": "s", "user": "u"}) + "\n",
                           encoding="utf-8")
        agent = root / "agent.json"
        agent.write_text(json.dumps([{"prompt_id": "p1", "status": "not_evaluated",
                                      "reason": "not judged in this run"}]), encoding="utf-8")
        with redirect_stdout(StringIO()):
            self.assertEqual(0, main(["llm-run", "--prompts", str(prompts), "--out-dir", str(root / "run"),
                                      "--provider", "codex-agent", "--model", "light",
                                      "--agent-file", str(agent), "--json"]))
        normalized = json.loads((root / "run/llm_candidate_responses.json").read_text(encoding="utf-8"))
        candidate = normalized["candidates"][0]
        self.assertEqual("not_evaluated", candidate["parse_status"])
        self.assertEqual("not judged in this run", candidate["not_evaluated_reason"])
        self.assertEqual([], candidate["claims"])

    def test_contradictory_declaration_is_rejected_not_repaired(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        prompts = root / "prompts.jsonl"
        prompts.write_text(json.dumps({"prompt_id": "p1", "gene_id": "", "system": "s", "user": "u"}) + "\n",
                           encoding="utf-8")
        agent = root / "agent.json"
        agent.write_text(json.dumps([{"prompt_id": "p1", "status": "not_evaluated", "reason": "x",
                                      "raw_text": "{\"claims\":[]}"}]), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "also supplies raw_text"):
            with redirect_stdout(StringIO()):
                main(["llm-run", "--prompts", str(prompts), "--out-dir", str(root / "run"),
                      "--provider", "codex-agent", "--model", "light", "--agent-file", str(agent), "--json"])


class LightEndToEndTests(unittest.TestCase):
    """light results -> convert -> extend -> llm-run -> import -> verify -> freeze."""

    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.bundle = self.root / "route_b_semantic_bundle.json"
        self.bundle.write_text(json.dumps(fixture()), encoding="utf-8")
        self.sha = digest(self.bundle.read_bytes())

    def build(self) -> dict:
        units = self.root / "units.txt"
        units.write_text("u1\nu2\nu3\nu4\nu5\nu6\n", encoding="utf-8")
        light_prepare(str(self.bundle), self.sha, str(self.root / "light"),
                      batch_units=20, units_file=str(units))
        index = self.root / "light/light/index.json"
        cache = self.root / "cache.json"
        cache.write_text(json.dumps({PAPER.casefold(): cache_entry(PAPER, PAPER_TITLE)}), encoding="utf-8")
        rows = [
            {"interval_gene_id": "g1", "reference_gene": "REF-A", "tier": "T2",
             "function": "leaflet boundary specification",
             "reason": "Homolog loss fuses leaflets, matching the observed phenotype.",
             "papers": [{"paper_id": PAPER, "title": PAPER_TITLE, "reported_gene": "HOM1",
                         "reported_species": "Arabidopsis thaliana", "gene_correspondence": "homolog",
                         "finding": "HOM1 loss fused leaflets in Arabidopsis.",
                         "checklist": checklist("T2")}], "literature": "screened 3"},
            {"interval_gene_id": "g4", "reference_gene": "REF-D", "tier": "T6",
             "function": "identity only", "reason": "Identity resolved but no functional bridge.", "papers": []},
            {"interval_gene_id": "g5", "reference_gene": "", "tier": "not_assessed",
             "function": "", "reason": "", "papers": []},
            {"interval_gene_id": "g6", "reference_gene": "", "tier": "not_assessed",
             "function": "", "reason": "", "papers": []},
        ]
        (self.root / "light/light/b01/results.jsonl").write_text(
            "".join(json.dumps(row) + "\n" for row in rows), encoding="utf-8")
        converted = light_convert(str(self.bundle), self.sha, str(index), digest(index.read_bytes()),
                                  str(self.root / "converted"), cache_arg=str(cache), offline=True)
        refreeze_route_b_sources(str(self.bundle), str(self.root / "converted/discoveries.json"),
                                 str(self.root / "extended"), parent_sha256=self.sha)
        extended = self.root / "extended"
        run = self.root / "run"
        with redirect_stdout(StringIO()):
            self.assertEqual(0, main([
                "llm-run", "--prompts", str(extended / "route_b_semantic_prompts.jsonl"),
                "--out-dir", str(run), "--provider", "codex-agent", "--model", "light",
                "--agent-file", str(self.root / "converted/agent_responses.json"), "--json",
            ]))
        imported = import_route_b_semantic_proposals(
            str(extended / "route_b_semantic_bundle.json"),
            str(extended / "route_b_semantic_prompts.jsonl"),
            str(run / "llm_candidate_responses.json"), str(self.root / "import"),
            llm_manifest_arg=str(run / "llm_run_manifest.json"), case_id="fixture-light")
        sources = verify_route_b_semantic_sources(
            imported["proposals"], str(self.root / "sources"),
            bundle_arg=str(extended / "route_b_semantic_bundle.json"),
            cache_arg=str(self.root / "converted/source_verify_cache.json"))
        transfer = verify_route_b_semantic_transfer(
            imported["proposals"], sources["source_verification"],
            str(extended / "route_b_semantic_bundle.json"), str(self.root / "transfer"),
            case_id="fixture-light")
        frozen = freeze_route_b_semantic_candidates(transfer["claims"], str(self.root / "freeze"),
                                                    case_id="fixture-light")
        return {"converted": converted, "imported": imported, "sources": sources,
                "transfer": transfer, "frozen": frozen, "run": run}

    def test_t2_paper_claim_reaches_adopted_first_stage(self):
        built = self.build()
        review = pd.read_csv(built["frozen"]["review_candidates"], sep="\t", dtype=str).fillna("")
        row = review[review["candidate_gene_id"] == "g1"].iloc[0]
        self.assertEqual("adopted_first_stage", row["candidate_admission_state"])
        self.assertEqual(0, built["frozen"]["candidate_gene_count"])
        transfer = pd.read_csv(built["transfer"]["transfer_verification"], sep="\t", dtype=str).fillna("")
        self.assertEqual({"eligible_initial_candidate"}, set(transfer["review_eligibility"]))
        self.assertIn("fused leaflets", row["initial_tier_reason"].lower())

    def test_not_assessed_unit_imports_as_not_evaluated_not_no_claims(self):
        built = self.build()
        status = pd.read_csv(built["imported"]["status"], sep="\t", dtype=str).fillna("").set_index("unit_id")
        self.assertEqual("not_evaluated", status.loc["u5", "parse_status"])
        self.assertNotEqual("no_claims", status.loc["u5", "status"])
        self.assertEqual("unparsed", status.loc["u5", "status"])
        normalized = json.loads((built["run"] / "llm_candidate_responses.json").read_text(encoding="utf-8"))
        response = next(row for row in normalized["candidates"]
                        if row["prompt_id"] == "route-b-semantic-u5")
        self.assertEqual("not_evaluated", response["parse_status"])
        self.assertEqual("not_assessed_or_missing_gene:g5", response["not_evaluated_reason"])

    def test_t6_genes_produce_no_claim(self):
        built = self.build()
        proposals = pd.read_csv(built["imported"]["proposals"], sep="\t", dtype=str).fillna("")
        self.assertEqual({"g1"}, set(proposals["member_interval_genes"]))
        status = pd.read_csv(built["imported"]["status"], sep="\t", dtype=str).fillna("").set_index("unit_id")
        self.assertEqual("no_claims", status.loc["u4", "status"])
        self.assertEqual("not_evaluated", status.loc["u6", "parse_status"])
        self.assertNotIn("u6", set(proposals["semantic_unit_id"]))



def with_database_source(value: dict, unit_id: str, gene_id: str, reference: str,
                         accession: str = "UniProtKB:Q9XXXX") -> dict:
    target = next(item for item in value["units"] if item["unit_id"] == unit_id)
    target["stable_source_ids"] = [accession]
    target["source_provenance"] = [{"source_ids": [accession], "reference_gene": reference,
                                    "source_database": "UniProt"}]
    gene = next(item for item in target["genes"] if item["gene_id"] == gene_id)
    gene.update(reference_source_id=accession, lookup_status="source_backed")
    return value


def paper_row(paper_id: str, title: str, tier_checks: dict, **extra) -> dict:
    row = {"paper_id": paper_id, "title": title, "reported_gene": "HOM1",
           "reported_species": "Arabidopsis thaliana", "gene_correspondence": "homolog",
           "finding": "HOM1 loss fused leaflets in Arabidopsis.", "checklist": tier_checks}
    row.update(extra)
    return row


def gene_row(gene_id: str, reference: str, tier: str, papers: list | None = None, **extra) -> dict:
    row = {"interval_gene_id": gene_id, "reference_gene": reference, "tier": tier,
           "function": "leaflet boundary specification",
           "reason": "The reference function bridges to the fused-leaflet phenotype.",
           "papers": papers or [], "literature": "screened 3"}
    row.update(extra)
    return row


def checks_for(**answers: str) -> dict[str, str]:
    checks = {key: "unknown" for key in CHECKLIST_FIELDS}
    checks.update(answers)
    return checks


class ConvertHarness(unittest.TestCase):
    """Prepare a light run on the fixture, write results, convert offline."""

    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def convert(self, rows: list[dict], *, mutate=None, cache: dict | None = None, offline: bool = True,
                fetch=None, unit_ids: str = "u1\nu2\nu3\nu4\nu5\nu6\nu7\n", name: str = "converted") -> dict:
        value = fixture()
        if mutate:
            mutate(value)
        bundle = self.root / "route_b_semantic_bundle.json"
        bundle.write_text(json.dumps(value), encoding="utf-8")
        self.sha = digest(bundle.read_bytes())
        units = self.root / "units.txt"
        units.write_text(unit_ids, encoding="utf-8")
        light_prepare(str(bundle), self.sha, str(self.root / "light"), batch_units=20, units_file=str(units))
        index = self.root / "light/light/index.json"
        (self.root / "light/light/b01/results.jsonl").write_bytes(
            "".join(json.dumps(row) + "\n" for row in rows).encode("utf-8"))
        cache_path = ""
        if cache is not None:
            cache_path = str(self.root / "cache.json")
            Path(cache_path).write_text(json.dumps(cache), encoding="utf-8")
        self.bundle, self.index = bundle, index
        return light_convert(str(bundle), self.sha, str(index), digest(index.read_bytes()),
                             str(self.root / name), cache_arg=cache_path, offline=offline, fetch=fetch)

    def read(self, name: str, folder: str = "converted"):
        return json.loads((self.root / folder / name).read_text(encoding="utf-8"))

    def responses(self) -> dict[str, dict]:
        return {row["prompt_id"]: row for row in self.read("agent_responses.json")["responses"]}

    def claims(self, unit_id: str) -> list[dict]:
        return json.loads(self.responses()[f"route-b-semantic-{unit_id}"]["raw_text"])["claims"]

    def paper_states(self) -> set[tuple[str, str]]:
        return {(row["paper_id"], row["state"]) for row in self.read("light_check.json")["paper_issues"]}


MAIN_CACHE = {PAPER.casefold(): cache_entry(PAPER, PAPER_TITLE)}


class LightFileBytesTests(ConvertHarness):
    def test_hashed_files_are_lf_even_when_text_mode_translates_newlines(self):
        original = Path.write_text

        def crlf_write_text(path, data, encoding=None, errors=None, newline=None):
            return original(path, data.replace("\n", "\r\n"), encoding=encoding, errors=errors, newline="")

        with mock.patch.object(Path, "write_text", crlf_write_text):
            self.convert([gene_row("g1", "REF-A", "T2", [paper_row(PAPER, PAPER_TITLE, checklist("T2"))])],
                         cache=MAIN_CACHE)
        directory = self.root / "light/light/b01"
        index = json.loads((self.root / "light/light/index.json").read_bytes())
        batch = index["batches"][0]
        for name, key in (("material.md", "material_sha256"), ("task.md", "task_sha256"),
                          ("roster.json", "roster_sha256")):
            data = (directory / name).read_bytes()
            self.assertNotIn(b"\r", data, name)
            self.assertEqual(digest(data), batch[key])
        self.assertNotIn(b"\r", (self.root / "light/light/index.json").read_bytes())
        self.assertNotIn(b"\r", (self.root / "converted/agent_responses.json").read_bytes())
        receipt = self.read("light_convert_receipt.json")
        for key in ("agent_responses", "metadata_cache", "source_verify_cache", "light_check"):
            entry = receipt["outputs"][key]
            self.assertEqual(digest((self.root / "converted" / entry["file"]).read_bytes()), entry["sha256"])


class LightCheckFixTests(ConvertHarness):
    def test_strong_paper_dropped_keeps_weaker_claim_and_records_the_tier_difference(self):
        t1 = checks_for(phenotype_reported="yes", focal_species="yes", equivalent_target_phenotype="yes")
        rows = [gene_row("g1", "REF-A", "T1", [
            paper_row("PMID:99999999", "A missing strong paper", t1),
            paper_row(PAPER, PAPER_TITLE, checklist("T3"))])]
        self.convert(rows, cache=MAIN_CACHE)
        (claim,) = self.claims("u1")
        self.assertEqual([PAPER], claim["source_ids"])
        self.assertNotIn("initial_tier", claim)
        self.assertNotIn("original researcher", claim["uncertainty_note"])
        (change,) = self.read("light_check.json")["tier_comparison"]
        self.assertEqual(("g1", "T1", "T3"), (change["gene_id"], change["declared_tier"], change["computed_tier"]))
        self.assertEqual([{"paper_id": "PMID:99999999", "state": "not_attempted"}], change["dropped_papers"])
        self.assertEqual(1, self.read("light_check.json")["totals"]["tier_differs_from_declared"])

    def test_equal_tier_records_no_difference(self):
        self.convert([gene_row("g1", "REF-A", "T2", [paper_row(PAPER, PAPER_TITLE, checklist("T2"))])],
                     cache=MAIN_CACHE)
        self.assertEqual([], self.read("light_check.json")["tier_comparison"])
        self.assertEqual("T2", self.claims("u1")[0]["initial_tier"])

    def test_new_paper_without_abstract_or_url_or_species_is_dropped_and_recorded(self):
        no_abstract = {**cache_entry("PMID:51", "A paper without an abstract"), "abstract": ""}
        http_url = {**cache_entry("PMID:52", "A paper with a bad url"), "url": "http://example.org/x"}
        no_retrieved = {**cache_entry("PMID:53", "A paper without retrieval time"), "retrieved_at": ""}
        cache = {"pmid:51": no_abstract, "pmid:52": http_url, "pmid:53": no_retrieved,
                 PAPER.casefold(): cache_entry(PAPER, PAPER_TITLE)}
        rows = [
            gene_row("g1", "REF-A", "T2", [
                paper_row("PMID:51", "A paper without an abstract", checklist("T2")),
                paper_row("PMID:52", "A paper with a bad url", checklist("T2")),
                paper_row("PMID:53", "A paper without retrieval time", checklist("T2")),
                paper_row(PAPER, PAPER_TITLE, checklist("T2"))]),
            gene_row("g2", "REF-B", "T2", [paper_row(PAPER, PAPER_TITLE, checklist("T2"),
                                                     reported_species="")]),
        ]
        self.convert(rows, cache=cache)
        self.assertEqual({("PMID:51", "missing_abstract"), ("PMID:52", "cache_entry_invalid"),
                          ("PMID:53", "cache_entry_invalid"), (PAPER, "discovery_not_formable")},
                         self.paper_states())
        self.assertEqual([PAPER], self.claims("u1")[0]["source_ids"])
        discoveries = self.read("discoveries.json")["discovered_sources"]
        self.assertEqual([("u1", PAPER)], [(entry["unit_id"], entry["source_id"]) for entry in discoveries])
        refreeze_route_b_sources(str(self.bundle), str(self.root / "converted/discoveries.json"),
                                 str(self.root / "extended"), parent_sha256=self.sha)

    def test_endpoint_without_a_frozen_species_cannot_form_a_discovery(self):
        self.convert([gene_row("g1", "fam-1", "T2", [paper_row(PAPER, PAPER_TITLE, checklist("T2"))])],
                     cache=MAIN_CACHE)
        self.assertIn((PAPER, "discovery_not_formable"), self.paper_states())
        self.assertEqual("not_evaluated", self.responses()["route-b-semantic-u1"]["status"])

    def test_species_is_normalised_and_comes_from_the_chosen_endpoint(self):
        def mutate(value):
            unit_value = value["units"][0]
            unit_value["reference_material_cards"][0]["reference_species"] = "Arabidopsis thaliana (thale cress)"
            unit_value["reference_material_cards"].append({
                "interval_gene_id": "g1", "reference_gene": "OTHER", "reference_species": "Oryza sativa",
                "function_description": "unrelated"})
            unit_value["genes"][0]["reference_species"] = "Arabidopsis thaliana (thale cress)"

        self.convert([gene_row("g1", "REF-A", "T2", [paper_row(PAPER, PAPER_TITLE, checklist("T2"))])],
                     mutate=mutate, cache=MAIN_CACHE)
        self.assertEqual("Arabidopsis thaliana", self.claims("u1")[0]["reference_species"])
        (entry,) = self.read("discoveries.json")["discovered_sources"]
        self.assertEqual("Arabidopsis thaliana", entry["reference_species"])

    def test_frozen_paper_needs_no_abstract_and_is_cited_by_its_frozen_spelling(self):
        entry = {**cache_entry("PMID:61", "A frozen paper"), "doi": "10.1000/frozen", "abstract": ""}
        cache = {"pmid:61": entry}
        self.convert([gene_row("g1", "REF-A", "T2", [paper_row("PMID:61", "A frozen paper", checklist("T2"))])],
                     mutate=lambda value: value["units"][0].update(stable_source_ids=["DOI:10.1000/frozen"]),
                     cache=cache)
        self.assertEqual(["DOI:10.1000/frozen"], self.claims("u1")[0]["source_ids"])
        self.assertEqual("DOI:10.1000/frozen", self.claims("u1")[0]["paper_reports_json"][0]["source_id"])
        self.assertFalse((self.root / "converted/discoveries.json").exists())

    def test_paper_ids_are_emitted_in_canonical_form(self):
        self.convert([gene_row("g1", "REF-A", "T2", [
            paper_row("https://pubmed.ncbi.nlm.nih.gov/11111111/", PAPER_TITLE, checklist("T2"))])],
            cache=MAIN_CACHE)
        self.assertEqual([PAPER], self.claims("u1")[0]["source_ids"])
        (entry,) = self.read("discoveries.json")["discovered_sources"]
        self.assertEqual(PAPER, entry["source_id"])

    def test_aliases_with_separators_are_not_emitted(self):
        record = {**cache_entry("PMID:71", "A paper with a semicolon DOI"), "doi": "10.1002/(SICI)1097;2-B"}
        self.convert([gene_row("g1", "REF-A", "T2", [
            paper_row("PMID:71", "A paper with a semicolon DOI", checklist("T2"))])],
            cache={"pmid:71": record})
        (entry,) = self.read("discoveries.json")["discovered_sources"]
        self.assertEqual("", entry["source_identifier_aliases"])
        refreeze_route_b_sources(str(self.bundle), str(self.root / "converted/discoveries.json"),
                                 str(self.root / "extended"), parent_sha256=self.sha)

    def test_a_primary_doi_with_a_separator_is_dropped_not_split(self):
        record = {**cache_entry("DOI:10.1002/(SICI)1097;2-B", "A paper with a semicolon DOI")}
        self.convert([gene_row("g1", "REF-A", "T2", [
            paper_row("DOI:10.1002/(SICI)1097;2-B", "A paper with a semicolon DOI", checklist("T2"))])],
            cache={"doi:10.1002/(sici)1097;2-b": record})
        self.assertIn(("DOI:10.1002/(SICI)1097;2-B", "identifier_not_representable"), self.paper_states())

    def test_no_paper_claim_without_a_frozen_database_source_is_dropped_and_recorded(self):
        self.convert([gene_row("g2", "REF-B", "T5")])
        check = self.read("light_check.json")
        self.assertEqual("no_paper_and_no_frozen_database_source_for_this_reference",
                         check["dropped_claims"][0]["reason"])
        self.assertEqual("not_evaluated", self.responses()["route-b-semantic-u2"]["status"])
        self.assertIn("machine_check_dropped_claim:g2", self.responses()["route-b-semantic-u2"]["reason"])

    def test_incomplete_checklist_is_dropped_not_filled_with_unknown(self):
        partial = {key: "yes" for key in CHECKLIST_FIELDS if key != "identity_resolved"}
        blank = checklist("T2") | {"identity_resolved": ""}
        shouting = checklist("T2") | {"identity_resolved": "YES"}
        cache = {**MAIN_CACHE, "pmid:81": cache_entry("PMID:81", "Paper A"),
                 "pmid:82": cache_entry("PMID:82", "Paper B"), "pmid:83": cache_entry("PMID:83", "Paper C")}
        self.convert([gene_row("g1", "REF-A", "T2", [
            paper_row("PMID:81", "Paper A", partial), paper_row("PMID:82", "Paper B", blank),
            paper_row("PMID:83", "Paper C", shouting),
            paper_row(PAPER, PAPER_TITLE, checklist("T2"))])], cache=cache)
        self.assertEqual({("PMID:81", "paper_report_incomplete"), ("PMID:82", "paper_report_incomplete"),
                          ("PMID:83", "paper_report_incomplete")}, self.paper_states())
        (claim,) = self.claims("u1")
        self.assertEqual([PAPER], claim["source_ids"])
        self.assertEqual(checklist("T2"), claim["paper_reports_json"][0]["paper_checklist"])

    def test_title_match_strips_markup_and_entities_and_stores_the_clean_title(self):
        cache = {PAPER.casefold(): cache_entry(PAPER, "A <i>homolog</i> controls leaflet &amp; stem development")}
        self.convert([gene_row("g1", "REF-A", "T2", [paper_row(
            PAPER, "A homolog controls leaflet & stem development", checklist("T2"))])], cache=cache)
        (claim,) = self.claims("u1")
        self.assertEqual("A homolog controls leaflet & stem development",
                         claim["paper_reports_json"][0]["source_title"])
        (entry,) = self.read("discoveries.json")["discovered_sources"]
        self.assertEqual("A homolog controls leaflet & stem development", entry["source_title"])

    def test_title_similarity_does_not_use_autojunk(self):
        long_title = " ".join(f"word{index}" for index in range(300))
        self.assertGreaterEqual(light._title_similarity(long_title, long_title + " extra"), 0.99)
        self.assertLess(light._title_similarity(long_title, long_title.replace("word", "term")), 0.85)

    def test_conflicting_duplicate_rows_leave_the_gene_not_assessed_identical_ones_collapse(self):
        first = gene_row("g1", "REF-A", "T2", [paper_row(PAPER, PAPER_TITLE, checklist("T2"))])
        conflicting = gene_row("g1", "REF-A", "T3", [paper_row(PAPER, PAPER_TITLE, checklist("T3"))])
        same = gene_row("g4", "REF-D", "T6")
        self.convert([first, conflicting, same, dict(same)], cache=MAIN_CACHE)
        issues = {(row["gene_id"], row["issue"]) for row in self.read("light_check.json")["gene_issues"]}
        self.assertIn(("g1", "duplicate_conflicting_rows_not_assessed"), issues)
        self.assertIn(("g4", "duplicate_identical_row_collapsed"), issues)
        self.assertEqual("not_evaluated", self.responses()["route-b-semantic-u1"]["status"])
        self.assertEqual([], self.claims("u4"))
        self.assertEqual(1, self.read("light_check.json")["totals"]["genes_not_assessed"])

    def test_results_file_hashes_are_recorded_per_batch(self):
        self.convert([gene_row("g4", "REF-D", "T6")])
        expected = digest((self.root / "light/light/b01/results.jsonl").read_bytes())
        receipt = self.read("light_convert_receipt.json")
        self.assertEqual(expected, receipt["input"]["results"]["b01"]["sha256"])
        self.assertEqual(expected, receipt["checks"]["batches"]["b01"]["results_sha256"])

    def test_provider_failure_is_counted_apart_from_not_found(self):
        def fetch(kind, identifier, timeout=30):
            if identifier == "1":
                raise OSError("network down")
            return None

        rows = [gene_row("g1", "REF-A", "T2", [paper_row("PMID:1", "Down", checklist("T2")),
                                               paper_row("PMID:2", "Absent", checklist("T2"))])]
        self.convert(rows, offline=False, fetch=fetch)
        totals = self.read("light_check.json")["totals"]
        self.assertEqual((1, 1), (totals["papers_provider_failure"], totals["papers_not_resolved"]))
        self.assertEqual({("PMID:1", "provider_failure"), ("PMID:2", "not_resolved")}, self.paper_states())

    def test_fetched_record_must_carry_the_requested_identifier(self):
        wrong = cache_entry("PMID:999", PAPER_TITLE)
        self.convert([gene_row("g1", "REF-A", "T2", [paper_row(PAPER, PAPER_TITLE, checklist("T2"))])],
                     offline=False, fetch=lambda kind, identifier, timeout=30: wrong)
        self.assertEqual({(PAPER, "not_resolved")}, self.paper_states())


class LightNetworkTests(unittest.TestCase):
    class Response:
        def __init__(self, payload):
            self.payload = payload

        def read(self):
            return json.dumps(self.payload).encode("utf-8")

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    def http_error(self, code):
        return urllib.error.HTTPError("https://example.org", code, "error", {}, None)

    def test_transient_errors_are_retried_with_backoff(self):
        calls = [self.http_error(429), self.http_error(503), self.Response({"ok": True})]
        with mock.patch.object(light.urllib.request, "urlopen", side_effect=calls) as opened, \
                mock.patch.object(light, "_sleep") as sleep:
            self.assertEqual({"ok": True}, light._http_json("https://example.org"))
        self.assertEqual(3, opened.call_count)
        self.assertEqual([1.0, 2.0], [call.args[0] for call in sleep.call_args_list])

    def test_persistent_failure_raises_after_three_tries_and_client_errors_do_not_retry(self):
        with mock.patch.object(light.urllib.request, "urlopen", side_effect=self.http_error(503)) as opened, \
                mock.patch.object(light, "_sleep"):
            with self.assertRaises(urllib.error.HTTPError):
                light._http_json("https://example.org")
        self.assertEqual(3, opened.call_count)
        with mock.patch.object(light.urllib.request, "urlopen", side_effect=self.http_error(404)) as opened, \
                mock.patch.object(light, "_sleep") as sleep:
            with self.assertRaises(urllib.error.HTTPError):
                light._http_json("https://example.org")
        self.assertEqual(1, opened.call_count)
        sleep.assert_not_called()

    def test_europe_pmc_result_for_another_identifier_is_not_accepted(self):
        payload = {"resultList": {"result": [
            {"pmid": "999", "title": "Another paper", "abstractText": "x"},
            {"pmid": "111", "title": "<i>Wanted</i> paper", "abstractText": "y"}]}}
        with mock.patch.object(light, "_http_json", return_value=payload):
            record = light._fetch_paper("pmid", "111")
            self.assertEqual(("111", "Wanted paper"), (record["pmid"], record["title"]))
            self.assertIsNone(light._fetch_paper("pmid", "222"))


class DeclaredStateTests(unittest.TestCase):
    def test_declared_not_evaluated_is_counted_by_the_route_b_state_reader(self):
        from candidate_qtl.route_b_response_validation import route_b_response_state

        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        prompts = root / "prompts.jsonl"
        prompts.write_text(json.dumps({"prompt_id": "p1", "gene_id": "", "system": "s", "user": "u"}) + "\n",
                           encoding="utf-8")
        agent = root / "agent.json"
        agent.write_text(json.dumps([{"prompt_id": "p1", "status": "not_evaluated", "reason": "not judged"}]),
                         encoding="utf-8")
        with redirect_stdout(StringIO()):
            self.assertEqual(0, main(["llm-run", "--prompts", str(prompts), "--out-dir", str(root / "run"),
                                      "--provider", "codex-agent", "--model", "light",
                                      "--agent-file", str(agent), "--json"]))
        candidate = json.loads((root / "run/llm_candidate_responses.json").read_text(encoding="utf-8"))["candidates"][0]
        self.assertEqual("declared_not_evaluated", route_b_response_state(candidate))


def run_pipeline(root: Path, value: dict, rows: list[dict] | dict, cache: dict, unit_ids: str,
                 *, output_contract: str = "legacy-v1") -> dict:
    bundle = root / "route_b_semantic_bundle.json"
    bundle.write_text(json.dumps(value), encoding="utf-8")
    sha = digest(bundle.read_bytes())
    units = root / "units.txt"
    units.write_text(unit_ids, encoding="utf-8")
    light_prepare(str(bundle), sha, str(root / "light"), batch_units=20, units_file=str(units), output_contract=output_contract)
    index = root / "light/light/index.json"
    (root / "cache.json").write_text(json.dumps(cache), encoding="utf-8")
    (root / "light/light/b01/results.jsonl").write_bytes(
        (json.dumps(rows) if output_contract == "science-v2" else "".join(json.dumps(row) + "\n" for row in rows)).encode("utf-8"))
    converted = light_convert(str(bundle), sha, str(index), digest(index.read_bytes()),
                              str(root / "converted"), cache_arg=str(root / "cache.json"), offline=True)
    refreeze_route_b_sources(str(bundle), str(root / "converted/discoveries.json"),
                             str(root / "extended"), parent_sha256=sha)
    extended, run = root / "extended", root / "run"
    with redirect_stdout(StringIO()):
        code = main(["llm-run", "--prompts", str(extended / "route_b_semantic_prompts.jsonl"),
                     "--out-dir", str(run), "--provider", "codex-agent", "--model", "light",
                     "--agent-file", str(root / "converted/agent_responses.json"), "--json"])
    assert code == 0
    imported = import_route_b_semantic_proposals(
        str(extended / "route_b_semantic_bundle.json"), str(extended / "route_b_semantic_prompts.jsonl"),
        str(run / "llm_candidate_responses.json"), str(root / "import"),
        llm_manifest_arg=str(run / "llm_run_manifest.json"), case_id="fixture-light")
    sources = verify_route_b_semantic_sources(
        imported["proposals"], str(root / "sources"), bundle_arg=str(extended / "route_b_semantic_bundle.json"),
        cache_arg=str(root / "converted/source_verify_cache.json"))
    transfer = verify_route_b_semantic_transfer(
        imported["proposals"], sources["source_verification"], str(extended / "route_b_semantic_bundle.json"),
        str(root / "transfer"), case_id="fixture-light")
    frozen = freeze_route_b_semantic_candidates(transfer["claims"], str(root / "freeze"), case_id="fixture-light")
    review = pd.read_csv(frozen["review_candidates"], sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
    return {"converted": converted, "imported": imported, "transfer": transfer, "frozen": frozen,
            "review": review, "run": run}


class LightNoPaperEndToEndTests(unittest.TestCase):
    """Tier survival after freeze for claims that carry no paper."""

    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def build(self) -> dict:
        value = fixture()
        with_database_source(value, "u2", "g2", "REF-B")
        with_database_source(value, "u3", "g3", "REF-C", accession="UniProtKB:Q9YYYY")
        t4_paper = "PMID:44444444"
        cache = {**MAIN_CACHE, t4_paper.casefold(): cache_entry(t4_paper, "A function paper")}
        rows = [
            gene_row("g1", "REF-A", "T2", [paper_row(PAPER, PAPER_TITLE, checklist("T2"))]),
            gene_row("g2", "REF-B", "T5"),
            gene_row("g3", "REF-C", "T4"),
            gene_row("g7", "REF-G", "T4", [paper_row(t4_paper, "A function paper", checks_for(
                function_supports_hypothesis="yes", identity_resolved="yes"))]),
        ]
        return run_pipeline(self.root, value, rows, cache, "u1\nu2\nu3\nu7\n")

    def test_annotation_only_t5_keeps_t5_through_freeze(self):
        built = self.build()
        review = built["review"]
        self.assertTrue(review.loc["g2", "relevance_tier"].startswith("T5"), review.loc["g2", "relevance_tier"])
        self.assertEqual("adopted_first_stage", review.loc["g2", "candidate_admission_state"])
        proposals = pd.read_csv(built["imported"]["proposals"], sep="\t", dtype=str).fillna("")
        row = proposals[proposals["member_interval_genes"] == "g2"].iloc[0]
        self.assertEqual("", row["paper_reports_json"])
        self.assertEqual('["annotation"]', row["evidence_types"])
        self.assertEqual("UniProtKB:Q9XXXX", row["source_id"])

    def test_database_reported_function_keeps_t4_through_freeze_without_support_promotion(self):
        built = self.build()
        self.assertTrue(built["review"].loc["g3", "relevance_tier"].startswith("T4"))
        self.assertEqual("adopted_first_stage", built["review"].loc["g3", "candidate_admission_state"])
        self.assertEqual("unverified_route_b_review", built["review"].loc["g3", "support_status"])
        check = json.loads((self.root / "converted/light_check.json").read_text(encoding="utf-8"))
        self.assertEqual([], check["tier_comparison"])
        claim = json.loads(json.loads((self.root / "converted/agent_responses.json").read_text(
            encoding="utf-8"))["responses"][2]["raw_text"])["claims"][0]
        self.assertEqual("T4", claim["initial_tier"])
        self.assertEqual(["function"], claim["evidence_types"])

    def test_t4_with_a_paper_and_t2_paper_claims_keep_their_tiers(self):
        review = self.build()["review"]
        self.assertTrue(review.loc["g7", "relevance_tier"].startswith("T4"), review.loc["g7", "relevance_tier"])
        self.assertTrue(review.loc["g1", "relevance_tier"].startswith("T2"), review.loc["g1", "relevance_tier"])


if __name__ == "__main__":
    unittest.main()
