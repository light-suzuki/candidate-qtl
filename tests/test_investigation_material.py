"""Generic provider-material regression; no biological answer fixtures."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from candidate_qtl.route_b_prepare import route_b_prepare
from candidate_qtl.route_b_semantic_bundle import build_route_b_semantic_bundle
from candidate_qtl.uniprot_reference_material import read_reference_directory


class InvestigationMaterialTests(unittest.TestCase):
    @unittest.skip("public release: uses the benchmark harness route_b_shards.py (not shipped)")
    def test_thousand_gene_local_records_reach_bounded_packets_without_support(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            records = root / "records"
            records.mkdir()
            record = {
                "primaryAccession": "P12345",
                "proteinDescription": {"recommendedName": {"fullName": {"value": "example protein"}}},
                "organism": {"scientificName": "Synthetic species", "taxonId": 123},
                "comments": [{"commentType": "SIMILARITY", "texts": [{"value": "Example family."}]}],
                "uniProtKBCrossReferences": [
                    {"database": "Pfam", "id": "PF12345", "properties": [{"key": "EntryName", "value": "Example domain"}]},
                    {"database": "GO", "id": "GO:1234567", "properties": [{"key": "GoTerm", "value": "F:example activity"}]},
                ],
                "references": [{"citation": {"citationCrossReferences": [{"database": "PubMed", "id": "123"}]}}],
            }
            (records / "P12345.json").write_text(json.dumps(record))
            (records / "failures.json").write_text('{"failed": ["R_missing"]}')
            genes = [f"G{i:04d}" for i in range(1000)]
            candidates, orthology = root / "genes.tsv", root / "orthology.tsv"
            pd.DataFrame({"gene_id": genes, "chrom": "chr1", "start": range(1000), "end": range(1, 1001)}).to_csv(candidates, sep="\t", index=False)
            pd.DataFrame({"gene_id": genes, "reference_gene_id": "P12345", "relationship": "best_hit_homology",
                          "orthogroup": [f"OG{i // 5:03d}" for i in range(1000)]}).to_csv(orthology, sep="\t", index=False)
            prepared = route_b_prepare(str(candidates), str(root / "b0"), orthology_arg=str(orthology),
                                       reference_function_args=[str(records)], auto_local_db=False)
            lookup = pd.read_csv(prepared["interval_reverse_lookup"], sep="\t", dtype=str).fillna("")
            self.assertEqual(set(lookup.gene_id), set(genes))
            for column, expected in {"protein_name": "example protein", "reference_species": "Synthetic species",
                                     "reference_domain_ids": "Pfam:PF12345", "protein_family": "Example family.",
                                     "verified_source_ids": "", "phenotype_source_ids": ""}.items():
                self.assertEqual(set(lookup[column]), {expected}, column)
            self.assertTrue(lookup.go_terms.str.contains("example activity").all())
            self.assertTrue(lookup.reference_function.eq("").all())
            semantic = build_route_b_semantic_bundle(prepared["interval_gene_atlas"], str(root / "semantic"),
                                                     phenotype="altered development", lookup_arg=prepared["interval_reverse_lookup"], case_id="scale")
            bundle = json.loads(Path(semantic["bundle"]).read_text())
            self.assertEqual(set(g for unit in bundle["units"] for g in unit["member_interval_genes"]), set(genes))
            self.assertTrue(all(unit["investigation_material"]["state"] == "annotation_only" for unit in bundle["units"]))
            self.assertTrue(all(unit["request_state"] == "requested" for unit in bundle["units"]))
            self.assertTrue(all(unit["source_material_class"] == "annotation_context" for unit in bundle["units"]))
            self.assertEqual(len(bundle["units"]), bundle["review_cost_estimate"]["estimated_model_call_count"])
            self.assertIsNone(bundle["review_cost_estimate"]["actual_input_tokens"])
            self.assertTrue(all(row["protein_name"] == "example protein" for unit in bundle["units"] for row in unit["genes"]))
            # Exercise the existing executor boundary without calling any model.
            path = Path(__file__).parents[1] / "workflow/route_b_shards.py"
            spec = importlib.util.spec_from_file_location("material_shards", path)
            sharder = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(sharder)
            max_units, max_chars = 8, 60000
            sharder.prepare("scale", root / "semantic", root / "shards", max_units=max_units, max_chars=max_chars)
            manifest = json.loads((root / "shards/prepare_manifest.json").read_text())
            self.assertEqual(manifest["candidate_count"], 1000)
            self.assertEqual(manifest["requested_unit_count"], len(bundle["units"]))
            self.assertGreater(len(manifest["shards"]), 1)
            canonical_prompts = {prompt["prompt_id"]: prompt for prompt in
                                 (json.loads(line) for line in Path(semantic["prompts"]).read_text().splitlines())}
            assigned_units, assigned_prompts, assigned_genes = [], [], []
            for shard in manifest["shards"]:
                packet_path = root / "shards" / shard["packet"]
                self.assertEqual(shard["packet_sha256"], hashlib.sha256(packet_path.read_bytes()).hexdigest())
                packet = json.loads(packet_path.read_text())
                self.assertEqual(shard["unit_count"], len(packet["assignments"]))
                self.assertLessEqual(shard["unit_count"], max_units)
                actual_chars = 0
                for assignment in packet["assignments"]:
                    path = packet_path.parent / assignment["file"]
                    text = path.read_text(encoding="utf-8")
                    prompt = json.loads(text)
                    self.assertEqual(packet["files"][assignment["file"]], hashlib.sha256(path.read_bytes()).hexdigest())
                    self.assertEqual(assignment["serialized_chars"], len(text))
                    self.assertEqual(canonical_prompts[assignment["prompt_id"]], prompt)
                    self.assertEqual(assignment["unit_id"], prompt["route_b_unit_id"])
                    actual_chars += len(text)
                    assigned_units.append(assignment["unit_id"])
                    assigned_prompts.append(assignment["prompt_id"])
                    assigned_genes.extend(json.loads(prompt["user"][prompt["user"].index("{"):])["member_interval_genes"])
                self.assertEqual(packet["serialized_chars"], actual_chars)
                oversized = actual_chars > max_chars
                self.assertEqual(packet["oversized_single_unit"], oversized)
                self.assertEqual(shard["oversized_single_unit"], oversized)
                if oversized:
                    # Preserve a full oversized unit; never truncate its material.
                    self.assertEqual(shard["unit_count"], 1)
            self.assertCountEqual(assigned_units, [unit["unit_id"] for unit in bundle["units"]])
            self.assertCountEqual(assigned_prompts, canonical_prompts)
            self.assertCountEqual(assigned_genes, genes)

    def test_duplicate_records_fail_and_citation_roles_stay_separate(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            record = {"primaryAccession": "P12345", "comments": [{"commentType": "DISRUPTION PHENOTYPE",
                      "texts": [{"value": "altered development", "evidences": [{"source": "PubMed", "id": "456"}]}]}],
                      "references": [{"citation": {"citationCrossReferences": [{"database": "PubMed", "id": "123"}]}}]}
            (root / "a.json").write_text(json.dumps(record))
            rows, _ = read_reference_directory(root)
            self.assertEqual(rows[0]["phenotype_source_ids"], "PMID:456")
            self.assertEqual(rows[0]["literature_source_ids"], "PMID:123")
            self.assertEqual(rows[0]["verified_source_ids"], "")
            (root / "b.json").write_text(json.dumps(record))
            with self.assertRaisesRegex(ValueError, "duplicate UniProt"):
                read_reference_directory(root)

    def test_uncharacterized_local_record_and_citation_are_retrieval_handles_only(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            records = root / "records"
            records.mkdir()
            record = {
                "primaryAccession": "P12345",
                "proteinDescription": {"recommendedName": {"fullName": {"value": "Uncharacterized protein"}}},
                "organism": {"scientificName": "Synthetic species", "taxonId": 123},
                "references": [{"citation": {"citationCrossReferences": [{"database": "PubMed", "id": "123"}]}}],
            }
            (records / "P12345.json").write_text(json.dumps(record), encoding="utf-8")
            candidates, orthology = root / "genes.tsv", root / "orthology.tsv"
            genes = ["G1", "G2", "G3"]
            pd.DataFrame({"gene_id": genes, "chrom": "chr1", "start": [0, 10, 20], "end": [9, 19, 29]}).to_csv(candidates, sep="\t", index=False)
            pd.DataFrame({"gene_id": genes, "reference_gene_id": "P12345", "relationship": "best_hit_homology", "orthogroup": "OG001"}).to_csv(orthology, sep="\t", index=False)
            prepared = route_b_prepare(str(candidates), str(root / "b0"), orthology_arg=str(orthology),
                                       reference_function_args=[str(records)], auto_local_db=False)
            result = build_route_b_semantic_bundle(prepared["interval_gene_atlas"], str(root / "semantic"),
                                                  phenotype="neutral development", lookup_arg=prepared["interval_reverse_lookup"], case_id="handles")
            bundle = json.loads(Path(result["bundle"]).read_text(encoding="utf-8"))
            self.assertEqual(set(genes), {gene for unit in bundle["units"] for gene in unit["member_interval_genes"]})
            self.assertEqual(0, result["prompt_count"])
            self.assertEqual("", Path(result["prompts"]).read_text(encoding="utf-8"))
            self.assertEqual(0, bundle["review_cost_estimate"]["estimated_input_tokens"])
            for unit in bundle["units"]:
                self.assertEqual("needs_material_retrieval", unit["investigation_material"]["state"])
                self.assertEqual("retrieval_handles_only", unit["source_material_class"])
                self.assertEqual("not_requested_no_reviewable_material", unit["request_state"])
                self.assertEqual("not_evaluated", unit["review_status"])
                self.assertEqual(0, unit["source_backed_function_record_count"])
                self.assertIn("PMID:123", unit["stable_source_ids"])
                self.assertTrue(unit["reference_material_cards"])
                for gene in unit["genes"]:
                    self.assertEqual("", gene["verified_source_ids"])
                    self.assertEqual("", gene["phenotype_source_ids"])
