"""Hit-only canonical material delivery; no functional/support promotion."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from candidate_qtl.route_b_semantic_bundle import INITIAL_MATERIAL_PROJECTION, build_route_b_semantic_bundle, refreeze_route_b_sources
from scripts.agentic_route_b_initial import collect_material, plan_initial_material, summarize_initial
from scripts.agentic_route_b_research import merge_initial_batches, plan_research_view


class HitInitialMaterialTests(unittest.TestCase):
    def test_batched_hit_material_reaches_merged_initial_ledger(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            atlas, hits = root / "atlas.tsv", root / "hits.tsv"
            pd.DataFrame({"gene_id": ["G1", "G2"]}).to_csv(atlas, sep="\t", index=False)
            pd.DataFrame({"focal_gene_id": ["G1", "G2"], "hit_accession": ["P11111", "P22222"],
                          "hit_protein_name": "DNA repair helicase"}).to_csv(hits, sep="\t", index=False)
            result = build_route_b_semantic_bundle(str(atlas), str(root / "semantic"), phenotype="neutral",
                                                   identity_hits_arg=str(hits), case_id="batched-hit")
            built = json.loads(Path(result["bundle"]).read_bytes())
            prompts = [json.loads(line) for line in Path(result["prompts"]).read_text().splitlines()]
            for legacy in (False, True):
                with self.subTest(legacy=legacy):
                    bundle = json.loads(json.dumps(built))
                    if legacy:
                        del bundle["initial_material_projection"]
                    source = root / f"bundle-{legacy}.json"
                    source.write_text(json.dumps(bundle), encoding="utf-8")
                    original = source.read_bytes()
                    sha = hashlib.sha256(original).hexdigest()
                    files, _ = plan_research_view(str(source), sha, prompts, bundle["case_id"], 2000, batch_units=1)
                    prepared = root / f"prepared-{legacy}"
                    for name, data in files.items():
                        target = prepared / name
                        target.parent.mkdir(parents=True, exist_ok=True)
                        target.write_bytes(data)
                    index_path = prepared / "research/batches/index.json"
                    index = json.loads(index_path.read_bytes())
                    self.assertEqual(2, len(index["batches"]))
                    for batch in index["batches"]:
                        manifest = json.loads(files[batch["manifest"]["file"]])
                        material = manifest["initial_material"]
                        self.assertEqual(int(not legacy), material["delivered_gene_count"])
                        groups = [g for page in material["pages"] for g in json.loads(files[page["file"]])]
                        self.assertEqual(int(not legacy), len(groups))
                        if legacy:
                            continue
                        binding, = groups[0]["bindings"]
                        self.assertEqual(batch["gene_ids"], binding["member_interval_genes"])
                        self.assertTrue(binding["source_pointer"].endswith("/sequence_similarity_hits/0"))
                        ledger = {"bundle_sha256": sha, "assessments": [{
                            "unit_id": batch["unit_ids"][0], "gene_id": batch["gene_ids"][0],
                            "state": "assessed_no_specific_bridge", "reason": "Original supplied annotation judgment.",
                            "used_material_pointers": [binding["source_pointer"]], "assessed_fields": ["hit_protein_name"],
                        }]}
                        (prepared / manifest["ledger"]["file"]).write_text(json.dumps(ledger), encoding="utf-8")
                    if not legacy:
                        merged = root / "merged"
                        merge_initial_batches(str(source), sha, str(index_path), hashlib.sha256(index_path.read_bytes()).hexdigest(), str(merged))
                        receipt = json.loads((merged / "merge_receipt.json").read_bytes())
                        self.assertEqual(2, receipt["totals"]["assessments"])
                        self.assertEqual(0, receipt["totals"]["dropped_record_count"])
                        summary = summarize_initial(str(source), sha, str(merged / "initial_assessments.json"), str(root / "batch-summary"))
                        self.assertEqual({"present": 2}, summary["material_counts"])
                        self.assertEqual({"interpreted": 2}, summary["interpretation_counts"])
                        self.assertEqual(0, summary["genes_with_fields_without_read_receipt"])
                    self.assertEqual(original, source.read_bytes())

    def test_refreeze_inherits_marker_and_annotation_accounting(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            atlas, hits = root / "atlas.tsv", root / "hits.tsv"
            pd.DataFrame({"gene_id": ["G1"]}).to_csv(atlas, sep="\t", index=False)
            pd.DataFrame({"focal_gene_id": ["G1"], "hit_accession": ["P11111"],
                          "hit_protein_name": ["DNA repair helicase"]}).to_csv(hits, sep="\t", index=False)
            result = build_route_b_semantic_bundle(str(atlas), str(root / "semantic"), phenotype="neutral",
                                                   identity_hits_arg=str(hits), case_id="refreeze-marker")
            built = json.loads(Path(result["bundle"]).read_bytes())
            for legacy in (False, True):
                with self.subTest(legacy=legacy):
                    bundle = json.loads(json.dumps(built))
                    if legacy:
                        del bundle["initial_material_projection"]
                        bundle["units"][0]["investigation_material"]["annotation_record_count"] = 0
                    parent = root / f"parent-{legacy}.json"
                    parent.write_text(json.dumps(bundle), encoding="utf-8")
                    original = parent.read_bytes()
                    sha = hashlib.sha256(original).hexdigest()
                    discovery = root / f"discovery-{legacy}.json"
                    discovery.write_text(json.dumps({"case_id": bundle["case_id"], "parent_bundle_sha256": sha,
                        "discovered_sources": [{
                            "unit_id": bundle["units"][0]["unit_id"], "reference_gene": "P11111",
                            "reference_species": "Species one", "reported_gene": "REF1", "reported_species": "Species one",
                            "source_id": "PMID:1234567", "source_title": "Synthetic observation",
                            "source_url": "https://example.test/observation", "retrieved_at": "2026-11-01T00:00:00Z",
                            "source_text": "Synthetic source context.", "reference_function": "DNA repair",
                        }]}), encoding="utf-8")
                    result = refreeze_route_b_sources(str(parent), str(discovery), str(root / f"refrozen-{legacy}"), parent_sha256=sha)
                    child = json.loads(Path(result["bundle"]).read_bytes())
                    self.assertEqual(not legacy, "initial_material_projection" in child)
                    if not legacy:
                        self.assertEqual(INITIAL_MATERIAL_PROJECTION, child["initial_material_projection"])
                    self.assertEqual(int(not legacy), child["units"][0]["investigation_material"]["annotation_record_count"])
                    self.assertEqual(1, child["units"][0]["investigation_material"]["functional_record_count"])
                    hit_groups = [g for g in collect_material(child)[0] if "hit_protein_name" in g["context"]]
                    self.assertEqual(int(not legacy), len(hit_groups))
                    self.assertEqual(original, parent.read_bytes())

    def test_builder_delivery_and_read_receipts_preserve_each_endpoint(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            atlas, hits_path = root / "atlas.tsv", root / "hits.tsv"
            genes = ["G1", "G2", "G3", "G4", "G5", "G6", "G7"]
            pd.DataFrame({"gene_id": genes, "family_id": "F1"}).to_csv(atlas, sep="\t", index=False)
            hits = [
                {"focal_gene_id": gene, "hit_accession": ref, "raw_reference_gene_id": raw,
                 "hit_protein_name": name, "hit_description": description,
                 "hit_species": species, "provider": "blastp", "database": database,
                 "percent_identity": "55", "query_coverage": "0.8", "evalue": "1e-20"}
                for gene, ref, raw, name, description, species, database in [
                    ("G1", "P11111", "sp|P11111|REF1", "DNA repair helicase", "sp|P11111|REF1 DNA repair helicase OS=Species one", "Species one", "db1"),
                    ("G1", "P22222", "P22222", "DNA repair helicase", "P22222 DNA repair helicase [Species two]", "Species two", "db2"),
                    ("G2", "P11111", "sp|P11111|REF1", "DNA repair helicase", "sp|P11111|REF1 DNA repair helicase OS=Species one", "Species one", "db1"),
                    ("G2", "P33333", "P33333", "unknown", "P33333 RNA helicase [Species three]", "Species three", "db3"),
                    ("G3", "P44444", "P44444", "", "", "Species four", "db4"),
                    ("G4", "P55555", "P55555", "hypothetical protein", "P55555 hypothetical protein [Species five]", "Species five", "db5"),
                    ("G5", "P66666", "P66666", "Uncharacterized protein", "P66666 Uncharacterized protein [Species six]", "Species six", "db6"),
                    ("G6", "P77777", "P77777", "P77777", "P77777", "Species seven", "db7"),
                    ("G7", "P88888", "P88888", "unknown", "P88888 OS=Species eight OX=8 GN=ABC", "Species eight", "db8"),
                ]
            ]
            pd.DataFrame(hits).to_csv(hits_path, sep="\t", index=False)
            result = build_route_b_semantic_bundle(str(atlas), str(root / "semantic"),
                                                   phenotype="altered development", identity_hits_arg=str(hits_path), case_id="hit-only")
            source = Path(result["bundle"])
            raw_bundle = source.read_bytes()
            sha = hashlib.sha256(raw_bundle).hexdigest()
            bundle = json.loads(raw_bundle)
            self.assertEqual(INITIAL_MATERIAL_PROJECTION, bundle["initial_material_projection"])
            unit = bundle["units"][0]
            self.assertEqual(4, unit["investigation_material"]["annotation_record_count"])
            self.assertEqual(0, unit["investigation_material"]["functional_record_count"])
            self.assertEqual(0, unit["source_backed_function_record_count"])
            self.assertEqual("annotation_context", unit["source_material_class"])
            self.assertEqual("requested", unit["request_state"])
            self.assertEqual([], unit["stable_source_ids"])
            prompt = json.loads(Path(result["prompts"]).read_text().splitlines()[0])
            self.assertIn("P33333 RNA helicase [Species three]", prompt["user"])
            groups, roster = collect_material(bundle)
            self.assertEqual(set(genes), set(roster))
            self.assertEqual([1, 1, 2], sorted(len(g["bindings"]) for g in groups))
            by_gene = {}
            for group in groups:
                self.assertEqual("UNVERIFIED_sequence_similarity_annotation_context", group["material_status"])
                self.assertNotIn("reference_function", group["context"])
                for binding in group["bindings"]:
                    gid, = binding["member_interval_genes"]
                    by_gene.setdefault(gid, []).append((group, binding))
                    # Every pointer resolves to the exact original frozen hit.
                    record = bundle
                    for part in binding["source_pointer"].strip("/").split("/"):
                        record = record[int(part)] if isinstance(record, list) else record[part]
                    self.assertEqual(record, binding["reference"])
                    self.assertEqual("sequence_similarity_handle_only", record["interpretation"])
                    self.assertEqual("none", record["rank_contribution"])
            self.assertEqual({"G1", "G2"}, set(by_gene))
            self.assertEqual({"P11111", "P22222"}, {b["reference"]["hit_accession"] for _, b in by_gene["G1"]})
            self.assertEqual({"P11111", "P33333"}, {b["reference"]["hit_accession"] for _, b in by_gene["G2"]})
            files, manifest = plan_initial_material(bundle, 1000)
            self.assertTrue(files)
            self.assertEqual(2, manifest["genes_with_material"])
            self.assertEqual(genes[2:], manifest["genes_without_material"])
            prompts = [json.loads(line) for line in Path(result["prompts"]).read_text().splitlines()]
            views, _ = plan_research_view(str(source), sha, prompts, bundle["case_id"], 1000)
            task = views["research/task.md"].decode("utf-8")
            self.assertIn("hit_protein_name/hit_description", task)
            self.assertIn("/units/{i}/genes/{j}/sequence_similarity_hits/{k}", task)
            self.assertIn("UNVERIFIED similarity context", task)
            pointer = by_gene["G1"][0][1]["source_pointer"]
            rows = [{"gene_id": "G1", "unit_id": unit["unit_id"], "state": "assessed_no_specific_bridge",
                     "reason": "Original annotation-context judgment only.", "used_material_pointers": [pointer],
                     "assessed_fields": ["hit_protein_name"]}]
            answers = root / "answers.json"
            answers.write_text(json.dumps({"bundle_sha256": sha, "assessments": rows}), encoding="utf-8")
            summary = summarize_initial(str(source), sha, str(answers), str(root / "summary"))
            self.assertEqual({"present": 2, "absent": 5}, summary["material_counts"])
            output = [json.loads(line) for line in (root / "summary/initial_assessment_genes.jsonl").read_text().splitlines()]
            gene = next(g for g in output if g["gene_id"] == "G1")
            self.assertEqual([pointer], gene["used_material_pointers"])
            self.assertEqual(["hit_description"], gene["fields_without_read_receipt"][pointer])
            self.assertFalse(gene["has_original_claim"])
            self.assertEqual(raw_bundle, source.read_bytes())
            files, followup = plan_initial_material(bundle, 1000, required_pointers={pointer})
            self.assertEqual(1, followup["delivered_gene_count"])
            self.assertTrue(all(b["source_pointer"] == pointer for data in files.values()
                                for group in json.loads(data) for b in group["bindings"]))
            rows[0]["gene_id"] = "G2"
            answers.write_text(json.dumps({"bundle_sha256": sha, "assessments": rows}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "pointer/field mismatch"):
                summarize_initial(str(source), sha, str(answers), str(root / "wrong-member"))
            self.assertFalse((root / "wrong-member").exists())

    def test_placeholder_only_units_remain_unscheduled(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            atlas, hits = root / "atlas.tsv", root / "hits.tsv"
            descriptions = ["P11111 hypothetical protein [Species one]", "P11111 [Arabidopsis thaliana]",
                            "sp|P11111|X [Source:UniProtKB/TrEMBL;Acc:P11111]",
                            "sp|P11111|X", "P11111 OS=Arabidopsis thaliana OX=3702", "P11111 hypothetical"]
            genes = [f"G{i}" for i in range(len(descriptions))]
            pd.DataFrame({"gene_id": genes}).to_csv(atlas, sep="\t", index=False)
            pd.DataFrame({"focal_gene_id": genes, "hit_accession": "P11111",
                          "hit_description": descriptions}).to_csv(hits, sep="\t", index=False)
            result = build_route_b_semantic_bundle(str(atlas), str(root / "semantic"), phenotype="neutral",
                                                   identity_hits_arg=str(hits))
            bundle = json.loads(Path(result["bundle"]).read_bytes())
            self.assertEqual(0, result["prompt_count"])
            for unit in bundle["units"]:
                self.assertEqual("retrieval_handles_only", unit["source_material_class"])
                self.assertEqual(0, unit["investigation_material"]["annotation_record_count"])
            self.assertEqual([], collect_material(bundle)[0])


if __name__ == "__main__":
    unittest.main()
