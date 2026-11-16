"""Offline source-handle regressions; no retrieval or biological scoring."""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from scripts.agentic_route_b_research import plan_research_view


def fixture():
    """Mirror exporter fields using synthetic IDs, not a known QTL answer."""
    return {
        "case_id": "source-fixture", "stage": "post_b0_semantic_bundle",
        "route": "interval_first", "phenotype": "synthetic observation",
        "interval_gene_count": 2, "semantic_unit_count": 1,
        "execution_prompt_count": 0, "execution_prompt_fingerprints": [],
        "units": [{
            "unit_id": "u1", "member_interval_genes": ["g1", "g2"],
            "request_state": "not_requested_no_reviewable_material",
            "source_material_class": "retrieval_handles_only", "identity": {},
            "genes": [{"gene_id": "g1"}, {"gene_id": "g2", "reference_gene": "ref-other"}],
            "reference_material_cards": [],
            "literature_link_candidates": [{
                "reference_gene": "ref-paper", "reference_function": "",
                "source_id": "", "source_url": "https://example.invalid/paper",
                "doi": "TEST:doi", "pmid": "TEST:pmid", "pmcid": "TEST:pmcid",
                "title": "Synthetic transport observation", "query_status": "metadata_only",
            }],
            "source_provenance": [{
                "gene_id": "g2", "reference_gene": "ref-other",
                "source_ids": ["TEST:record", "TEST:paper"],
                "source_url": "https://example.invalid/record",
                "evidence_type": "fixture-only", "reference_species": "test species",
            }],
            "literature_query_status": [],
        }],
    }


class ResearchSourceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.bundle = fixture()

    def plan(self, max_bytes=262144):
        raw = json.dumps(self.bundle, ensure_ascii=False).encode()
        source = self.root / "bundle.json"
        source.write_bytes(raw)
        files, metadata = plan_research_view(
            str(source), hashlib.sha256(raw).hexdigest(), [], "source-fixture", max_bytes,
        )
        self.assertEqual(raw, source.read_bytes())
        rows = [json.loads(line) for name, data in sorted(files.items())
                if name.endswith(".jsonl") for line in data.splitlines()]
        return files, metadata, rows

    def resolve(self, pointer):
        value = self.bundle
        for token in pointer.strip("/").split("/"):
            value = value[int(token)] if isinstance(value, list) else value[token]
        return value

    def test_title_only_paper_is_visible_without_formal_prompt(self):
        _, metadata, rows = self.plan()
        paper = rows[0]["literature_link_previews"][0]
        original = self.bundle["units"][0]["literature_link_candidates"][0]
        for key in ("doi", "pmid", "pmcid", "source_url", "reference_gene"):
            self.assertEqual(original[key], paper[key])
        self.assertEqual(original["title"], paper["cues"]["title"])
        self.assertEqual(original, self.resolve(paper["source_pointer"]))
        self.assertEqual("navigation_only_not_evaluated", rows[0]["view_status"])
        self.assertEqual("not_started", metadata["execution_state"])
        self.assertNotIn("claims", rows[0])

    def test_all_papers_and_duplicate_occurrences_keep_their_own_pointer(self):
        articles = self.bundle["units"][0]["literature_link_candidates"]
        articles.extend([dict(articles[0], reference_gene="ref-second"), copy.deepcopy(articles[0])])
        _, _, rows = self.plan()
        previews = rows[0]["literature_link_previews"]
        self.assertEqual(3, len(previews))
        for index, preview in enumerate(previews):
            self.assertEqual(f"/units/0/literature_link_candidates/{index}", preview["source_pointer"])
            self.assertEqual(articles[index], self.resolve(preview["source_pointer"]))

    def test_long_title_is_marked_but_source_identifiers_are_not_cut(self):
        article = self.bundle["units"][0]["literature_link_candidates"][0]
        article.update(title="題" * 700, doi="TEST:" + "x" * 500)
        _, _, rows = self.plan()
        preview = rows[0]["literature_link_previews"][0]
        self.assertEqual("題" * 240, preview["cues"]["title"])
        self.assertIn("title", preview["truncated_fields"])
        self.assertEqual(article["doi"], preview["doi"])
        self.assertEqual(article["title"], self.resolve(preview["source_pointer"])["title"])

    def test_provenance_binding_does_not_spread_to_other_gene(self):
        _, _, rows = self.plan()
        preview = rows[0]["source_provenance_previews"][0]
        original = self.bundle["units"][0]["source_provenance"][0]
        self.assertEqual(original["source_ids"], preview["source_ids"])
        self.assertEqual("g2", preview["gene_id"])
        self.assertEqual("ref-other", preview["reference_gene"])
        self.assertEqual(original, self.resolve(preview["source_pointer"]))
        self.assertNotIn("reference_gene", rows[0]["genes"][0])
        self.assertNotIn("support", rows[0])

    def test_query_failures_no_hit_and_missing_history_remain_different(self):
        histories = [[], [{"reference_gene": "ref-paper", "query_status": "no_hit"}],
                     [{"reference_gene": "ref-paper", "query_status": "provider_failure",
                       "reason": "synthetic timeout", "unrecognized_provider_field": "preserve me"}]]
        for history in histories:
            with self.subTest(history=history):
                self.bundle["units"][0]["literature_query_status"] = history
                _, _, rows = self.plan()
                records = rows[0]["literature_query_status_records"]
                self.assertEqual(history, [item["record"] for item in records])
                for index, item in enumerate(records):
                    self.assertEqual(history[index], self.resolve(item["source_pointer"]))
                self.assertEqual("navigation_only_not_evaluated", rows[0]["view_status"])

    def test_reference_card_ids_and_species_are_preserved(self):
        card = {"interval_gene_id": "g1", "reference_gene_id": "namespace:ref-card",
                "reference_species": "card species", "reference_accession_namespace": "namespace",
                "source_ids": ["TEST:card"], "title": "card title"}
        self.bundle["units"][0]["reference_material_cards"] = [card]
        _, _, rows = self.plan()
        preview = rows[0]["reference_card_previews"][0]
        for key in ("reference_gene_id", "reference_species", "reference_accession_namespace", "source_ids"):
            self.assertEqual(card[key], preview[key])
        self.assertEqual(card, self.resolve(preview["source_pointer"]))

    def test_canonical_annotation_leads_and_full_hit_pointer_are_visible(self):
        gene = self.bundle["units"][0]["genes"][0]
        fields = ("local_annotation_description", "local_annotation_go_terms",
                  "reference_domain_description", "protein_family", "family_label",
                  "keywords", "primary_broad_function")
        for field in fields:
            gene[field] = "specific functional lead " + field
        gene["domain_ids"] = "DOMAIN:" + "x" * 500
        gene["reference_domain_ids"] = "DOMAIN:ref"
        gene["sequence_similarity_hits"] = [{"hit_accession": "exact-accession",
                                              "hit_description": "long functional clue " * 100}]
        _, _, rows = self.plan()
        preview = rows[0]["genes"][0]
        for field in fields:
            self.assertEqual(gene[field], preview["cues"][field])
        self.assertEqual(gene["domain_ids"], preview["domain_ids"])
        self.assertEqual(gene["reference_domain_ids"], preview["reference_domain_ids"])
        hit = preview["sequence_handles"][0]
        self.assertEqual(240, len(hit["hit_description"]))
        self.assertEqual(["hit_description"], hit["truncated_fields"])
        self.assertEqual(gene["sequence_similarity_hits"][0], self.resolve(hit["source_pointer"]))
        self.assertEqual("exact-accession", hit["hit_accession"])

    def test_invalid_source_tables_are_not_silently_hidden(self):
        for field in ("literature_link_candidates", "source_provenance", "literature_query_status"):
            for invalid in (None, {}, ["invalid"]):
                with self.subTest(field=field, invalid=invalid):
                    self.bundle = fixture()
                    self.bundle["units"][0][field] = invalid
                    with self.assertRaises(ValueError):
                        self.plan()

    def test_legacy_absent_tables_are_empty_not_negative(self):
        for field in ("literature_link_candidates", "source_provenance", "literature_query_status"):
            self.bundle["units"][0].pop(field)
        _, _, rows = self.plan()
        for field in ("literature_link_previews", "source_provenance_previews", "literature_query_status_records"):
            self.assertEqual([], rows[0][field])
        self.assertEqual("navigation_only_not_evaluated", rows[0]["view_status"])

    def test_prepare_cli_with_sources_preserves_formal_artifacts(self):
        from scripts.agentic_route_b_material import restore_reading_shard

        for requested in (False, True):
            for affinity in (False, True):
                with self.subTest(requested=requested, affinity=affinity):
                    bundle = copy.deepcopy(self.bundle)
                    prompts = []
                    if requested:
                        unit = bundle["units"][0]
                        unit["request_state"] = "requested"
                        unit["source_material_class"] = "annotation_context"
                        unit["genes"][0]["reference_function"] = "synthetic cue"
                        prompts = [{"case_id": "source-fixture", "prompt_id": "route-b-semantic-u1",
                                    "route_b_unit_id": "u1", "system": "TEST ONLY",
                                    "user": "Phenotype: synthetic observation\n\n" + json.dumps(unit, ensure_ascii=False, sort_keys=True)}]
                    bundle["execution_prompt_count"] = len(prompts)
                    bundle["execution_prompt_fingerprints"] = [
                        hashlib.sha256(json.dumps(p, ensure_ascii=False, sort_keys=True,
                                                  separators=(",", ":")).encode()).hexdigest()
                        for p in prompts]
                    source = self.root / "bundle.json"
                    source.write_text(json.dumps(bundle), encoding="utf-8")
                    sha = hashlib.sha256(source.read_bytes()).hexdigest()
                    prompt_file = self.root / "prompts.jsonl"
                    prompt_file.write_text("".join(json.dumps(p) + "\n" for p in prompts), encoding="utf-8")
                    before = source.read_bytes(), prompt_file.read_bytes()
                    outputs, manifests = [], []
                    for research in (False, True):
                        out = self.root / f"cli-{requested}-{affinity}-{research}"
                        command = [sys.executable, str(Path(__file__).resolve().parents[1] / "scripts/agentic_route_b.py"),
                                   "prepare", "--prompts", str(prompt_file), "--out", str(out),
                                   "--case-id", "source-fixture"]
                        if affinity:
                            command.append("--prefer-shared-material")
                        if research:
                            command.extend(["--research-bundle", str(source), "--research-bundle-sha256", sha])
                        result = subprocess.run(command, capture_output=True, text=True, check=True)
                        outputs.append(out)
                        manifests.append(json.loads((out / "prepare_manifest.json").read_text()))
                        self.assertEqual(requested, bool(json.loads(result.stdout)["prompt_count"]))
                    self.assertEqual(manifests[0], {k: v for k, v in manifests[1].items() if k != "research_view"})
                    for file in outputs[0].rglob("*"):
                        if file.is_file() and file.name != "prepare_manifest.json":
                            self.assertEqual(file.read_bytes(), (outputs[1] / file.relative_to(outputs[0])).read_bytes())
                    restored = [p for shard in manifests[1]["reading_view"]["shards"]
                                for p in restore_reading_shard(outputs[1], shard["file"])]
                    self.assertEqual(prompts, restored)
                    view = json.loads((outputs[1] / "research/manifest.json").read_text())
                    row = json.loads((outputs[1] / view["pages"][0]["file"]).read_text())
                    self.assertEqual("TEST:doi", row["literature_link_previews"][0]["doi"])
                    self.assertEqual(before, (source.read_bytes(), prompt_file.read_bytes()))

    def test_added_source_material_obeys_existing_page_accounting(self):
        files, metadata, rows = self.plan(max_bytes=1)
        manifest = json.loads(files["research/manifest.json"])
        self.assertEqual(2, metadata["gene_count"])
        self.assertEqual(["g1", "g2"], rows[0]["member_interval_genes"])
        page = manifest["pages"][0]
        self.assertEqual(len(files[page["file"]]), page["bytes"])
        self.assertEqual(page["bytes"] - 1, page["over_bound_bytes"])
        self.assertEqual(hashlib.sha256(files[page["file"]]).hexdigest(), page["sha256"])
        self.assertEqual(1, len(rows[0]["literature_link_previews"]))


if __name__ == "__main__":
    unittest.main()
