"""Compact display must preserve copy-specific fields, duplicates and sources."""

import json
from pathlib import Path
import tempfile
import unittest
from scripts.agentic_route_b import prepare as prepare_reading
from scripts.agentic_route_b_compact import prepare, restore


class TestCompactRouteB(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.original = self.root / "original"
        self.compact = self.root / "compact"
        prompts = []
        for i in range(2):
            gene = f"gene-copy-{i}"
            card = {"interval_gene_id": gene, "reference_gene": "REF", "reference_species": "Species", "function_description": "Full nonempty evidence 日本語", "source_id": "PMID:1"}
            unit = {
                "unit_id": f"unit-{i}",
                "member_interval_genes": [gene],
                "genes": [
                    {"gene_id": gene, "local_annotation_description": f"copy effect {i}", "reference_gene": "REF", "reference_species": "Species", "empty": "", "null": None, "flag": False, "zero": 0}
                ],
                "reference_material_cards": [card, dict(card)],
                "source_provenance": [],
                "literature_link_candidates": [],
            }
            prompts.append(
                {
                    "case_id": "case-x",
                    "prompt_id": f"prompt-{i}",
                    "route_b_unit_id": f"unit-{i}",
                    "gene_id": gene,
                    "system": "Read every copy",
                    "user": "Phenotype test\n\n" + json.dumps(unit, ensure_ascii=False, sort_keys=True),
                }
            )
        source = self.root / "prompts.jsonl"
        source.write_text("".join(json.dumps(x, ensure_ascii=False) + "\n" for x in prompts), encoding="utf-8")
        prepare_reading(str(source), str(self.original), max_prompts=1)

    def tearDown(self):
        self.temp.cleanup()

    def test_full_assignments_restore_with_duplicates_copies_and_empty_types(self):
        m = prepare(self.original, self.compact)
        self.assertTrue(m["exact_roundtrip_verified"])
        self.assertEqual(m["prompt_count"], 2)
        genes = []
        for s in m["shards"]:
            a = restore(self.compact, s["file"])
            original = json.loads((self.original / s["original_file"]).read_text(encoding="utf-8"))
            self.assertEqual(a, original)
            for p in a["prompts"]:
                genes.extend(p["unit"]["member_interval_genes"])
                binding = p["unit_tables"]["genes"][0]["binding"]
                self.assertEqual(binding["empty"], "")
                self.assertIsNone(binding["null"])
                self.assertIs(binding["flag"], False)
                self.assertEqual(binding["zero"], 0)
                self.assertEqual(len(p["unit_tables"]["reference_material_cards"]), 2)
        self.assertEqual(genes, ["gene-copy-0", "gene-copy-1"])

    def test_tampered_display_rejected(self):
        m = prepare(self.original, self.compact)
        s = m["shards"][0]
        p = self.compact / s["file"]
        p.write_bytes(p.read_bytes() + b" ")
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            restore(self.compact, s["file"])

    def test_tampered_material_rejected(self):
        m = prepare(self.original, self.compact)
        material = m["materials"][0]
        p = self.compact / material["file"]
        p.write_bytes(b"{}")
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            restore(self.compact, m["shards"][0]["file"])

    def test_unverified_source_rejected_before_output(self):
        p = self.original / "prepare_manifest.json"
        m = json.loads(p.read_text())
        m["reading_view"]["exact_roundtrip_verified"] = False
        p.write_text(json.dumps(m))
        with self.assertRaisesRegex(ValueError, "not verified"):
            prepare(self.original, self.compact)
        self.assertFalse(self.compact.exists())

    def test_existing_output_preserved(self):
        prepare(self.original, self.compact)
        with self.assertRaisesRegex(ValueError, "already exists"):
            prepare(self.original, self.compact)

    def test_unknown_shard_rejected(self):
        prepare(self.original, self.compact)
        with self.assertRaisesRegex(ValueError, "unknown compact shard"):
            restore(self.compact, "reading/unknown.json")
