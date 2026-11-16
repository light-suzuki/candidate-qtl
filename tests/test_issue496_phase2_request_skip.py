"""Issue #496 Phase 2: do not spend a request on a unit that cannot answer.

IDs and bare family/domain labels are retrieval handles, not reviewable
function. A unit with no descriptive annotation, reference function or
literature therefore has exactly one legal answer: an empty claim list.
Skipping that request costs no information; the unit stays in the bundle, stays
in the candidate universe, and is imported as an explicit not_evaluated row.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.route_b_semantic_bundle import build_route_b_semantic_bundle
from candidate_qtl.route_b_semantic_import import import_route_b_semantic_proposals

ATLAS_COLUMNS = [
    "gene_id", "evolutionary_family_id", "evolutionary_family_status", "orthogroup",
    "family_id", "domain_architecture", "domain_ids", "orthology_reference_gene",
    "subfamily_id", "broad_function_bins", "conservation_class",
    "evidence_completeness", "notes", "local_annotation_description",
]
ROWS = [
    # An ID alone is not functional material.
    {"gene_id": "G_FAMILY", "family_id": "NOOT", "conservation_class": "relative_enriched"},
    # A named domain architecture is reviewable annotation.
    {"gene_id": "G_DOMAIN", "domain_architecture": "Pkinase", "conservation_class": "relative_enriched"},
    # A reference accession alone does not describe function.
    {"gene_id": "G_REFERENCE", "orthology_reference_gene": "Q2HW56", "conservation_class": "relative_enriched"},
    # Descriptive annotation supplies reviewable material.
    {"gene_id": "G_ANNOTATED", "local_annotation_description": "receptor kinase regulates meristem development",
     "conservation_class": "relative_enriched"},
    # Not reviewable: nothing at all.
    {"gene_id": "G_DARK1", "conservation_class": "no_detectable_homology"},
    {"gene_id": "G_DARK2", "conservation_class": "no_detectable_homology"},
]


def _write_atlas(path: Path, rows: list[dict[str, str]]) -> None:
    path.write_text(
        "\t".join(ATLAS_COLUMNS) + "\n"
        + "".join("\t".join(str(row.get(column, "")) for column in ATLAS_COLUMNS) + "\n"
                  for row in rows),
        encoding="utf-8")


class RequestSkipTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        atlas = self.root / "atlas.tsv"
        _write_atlas(atlas, ROWS)
        self.result = build_route_b_semantic_bundle(
            str(atlas), str(self.root / "out"), phenotype="determinate inflorescence",
            case_id="issue496p2")
        self.bundle = json.loads(Path(self.result["bundle"]).read_text(encoding="utf-8"))
        self.manifest = json.loads(Path(self.result["manifest"]).read_text(encoding="utf-8"))

    def test_units_without_material_are_not_requested(self):
        states = {unit["unit_id"]: unit["request_state"] for unit in self.bundle["units"]}
        by_gene = {unit["member_interval_genes"][0]: unit["request_state"]
                   for unit in self.bundle["units"] if len(unit["member_interval_genes"]) == 1}
        self.assertEqual("not_requested_no_reviewable_material", by_gene["G_FAMILY"])
        self.assertEqual("requested", by_gene["G_DOMAIN"])
        self.assertEqual("not_requested_no_reviewable_material", by_gene["G_REFERENCE"])
        self.assertEqual("requested", by_gene["G_ANNOTATED"])
        self.assertEqual("not_requested_no_reviewable_material", by_gene["G_DARK1"])
        self.assertEqual("not_requested_no_reviewable_material", by_gene["G_DARK2"])
        self.assertEqual(6, len(states))

    def test_prompts_cover_exactly_the_requested_units(self):
        prompts = [json.loads(line) for line in
                   Path(self.result["prompts"]).read_text(encoding="utf-8").splitlines() if line.strip()]
        requested = {unit["unit_id"] for unit in self.bundle["units"]
                     if unit["request_state"] == "requested"}
        self.assertEqual(requested, {prompt["route_b_unit_id"] for prompt in prompts})
        self.assertEqual(2, len(prompts))

    def test_skipped_units_are_ledgered_and_counted(self):
        ledger = pd.read_csv(self.result["units_not_requested"], sep="\t", dtype=str).fillna("")
        self.assertEqual(4, len(ledger))
        self.assertEqual({"not_evaluated"}, set(ledger["review_status"]))
        self.assertEqual({"none"}, set(ledger["rank_contribution"]))
        observability = self.manifest["compression_observability"]
        self.assertEqual(6, observability["semantic_unit_count"])
        self.assertEqual(2, observability["requested_unit_count"])
        self.assertEqual(4, observability["skipped_unit_count"])
        self.assertEqual(4, observability["skipped_gene_count"])
        self.assertEqual(2, observability["expected_llm_request_count"])
        self.assertEqual(33.333333, observability["expected_requests_per_100_genes"])

    def test_no_candidate_leaves_the_universe(self):
        members = [gene for unit in self.bundle["units"] for gene in unit["member_interval_genes"]]
        self.assertEqual(sorted(row["gene_id"] for row in ROWS), sorted(members))

    def test_placeholder_values_do_not_count_as_material(self):
        """`unknown` fills the column but says nothing.

        B0 writes `primary_broad_function=unknown` for every unclassified gene.
        Counting that as material made the skip a no-op on real data: the
        200-gene bench reported a 0% reduction until placeholders were excluded.
        """
        atlas = self.root / "placeholder_atlas.tsv"
        _write_atlas(atlas, [
            {"gene_id": "G_PLACEHOLDER", "conservation_class": "no_detectable_homology",
             "evidence_completeness": "unknown", "notes": "unknown"},
            {"gene_id": "G_REAL", "family_id": "LYK", "conservation_class": "relative_enriched",
             "local_annotation_description": "receptor kinase regulates root hair development"},
        ])
        result = build_route_b_semantic_bundle(
            str(atlas), str(self.root / "placeholder_out"), phenotype="trait", case_id="issue496p2")
        bundle = json.loads(Path(result["bundle"]).read_text(encoding="utf-8"))
        states = {unit["member_interval_genes"][0]: unit["request_state"] for unit in bundle["units"]}
        self.assertEqual("not_requested_no_reviewable_material", states["G_PLACEHOLDER"])
        self.assertEqual("requested", states["G_REAL"])
        self.assertEqual(1, result["prompt_count"])

    def test_import_records_skipped_units_as_not_evaluated(self):
        prompts = [json.loads(line) for line in
                   Path(self.result["prompts"]).read_text(encoding="utf-8").splitlines() if line.strip()]
        run_dir = self.root / "llm_run"
        run_dir.mkdir()
        responses = {"candidates": [
            {"prompt_id": prompt["prompt_id"], "gene_id": prompt["gene_id"], "claims": [],
             "model_json": {"claims": []}, "raw_text": '{"claims":[]}', "parse_status": "json",
             "review_status": "not_reviewed"}
            for prompt in prompts]}
        response_path = run_dir / "llm_candidate_responses.json"
        response_path.write_text(json.dumps(responses, ensure_ascii=False), encoding="utf-8")
        imported = import_route_b_semantic_proposals(
            self.result["bundle"], self.result["prompts"], str(response_path),
            str(self.root / "import"), case_id="issue496p2")
        status = pd.read_csv(imported["status"], sep="\t", dtype=str).fillna("")
        skipped = status[status["status"] == "not_evaluated_no_reviewable_material"]
        self.assertEqual(4, len(skipped))
        self.assertEqual({"not_requested"}, set(skipped["parse_status"]))
        # Every unit is still accounted for, requested or not.
        self.assertEqual(6, len(status))


if __name__ == "__main__":
    unittest.main()
