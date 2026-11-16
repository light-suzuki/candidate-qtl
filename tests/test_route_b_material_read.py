"""Bounded local reads preserve exact facts and report unselected material unread."""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from scripts import agentic_route_b_light as light
from scripts.agentic_route_b_material_read import material_read, source_read
from tests.test_route_b_light import digest, fixture


class MaterialReadTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.bundle = fixture()
        self.long = "Original phenotype evidence " * 40
        self.bundle["units"][0]["genes"][0]["reference_function"] = self.long
        self.bundle_path = self.root / "bundle.json"
        self.bundle_path.write_text(json.dumps(self.bundle), encoding="utf-8")
        self.bundle_sha = digest(self.bundle_path.read_bytes())
        light.light_prepare(str(self.bundle_path), self.bundle_sha, str(self.root / "prepared"), batch_units=2)
        self.index_path = self.root / "prepared/light/index.json"
        self.index_sha = digest(self.index_path.read_bytes())

    def read(self, **kwargs):
        return material_read(str(self.bundle_path), self.bundle_sha, str(self.index_path), self.index_sha,
                             batch_id="b01", **kwargs)

    def test_default_is_compact_index_without_material_content(self):
        packet = self.read()
        self.assertEqual("index_only", packet["state"])
        self.assertEqual("not_read", packet["read_state"])
        self.assertNotIn("selections", packet)
        self.assertNotIn(self.long, json.dumps(packet))

    def test_exact_gene_and_trait_fields_keep_over_400_characters(self):
        packet = self.read(gene_ids=["g1"], fields=["reference_function"])
        self.assertEqual(self.long, packet["selections"][0]["value"])
        self.assertEqual("/units/0/genes/0/reference_function", packet["selected_pointers"][0])
        self.assertIn("not_read", packet["omitted"])
        trait = self.read(gene_ids=["g2"], fields=["mutant_phenotype"])
        self.assertEqual(["/units/1/reference_material_cards/0/mutant_phenotype"], trait["selected_pointers"])
        self.assertNotIn("g1", json.dumps(trait["selections"]))

    def test_unit_pointer_and_bare_field_select_without_clipping(self):
        packet = self.read(unit_ids=["u1"], fields=["/genes/0/reference_function"])
        self.assertEqual(self.long, packet["selections"][0]["value"])
        self.assertEqual("u1", self.read(unit_ids=["u1"], fields=["unit_id"])["selections"][0]["value"])

    def test_oversize_returns_no_partial_content(self):
        packet = self.read(gene_ids=["g1"], fields=["reference_function"], max_chars=400)
        self.assertEqual("narrow_selection_required", packet["state"])
        self.assertEqual("not_read", packet["read_state"])
        self.assertNotIn("selections", packet)
        self.assertNotIn(self.long[:400], json.dumps(packet))
        self.assertEqual(self.long, self.read(gene_ids=["g1"], fields=["reference_function"], max_chars=2000)["selections"][0]["value"])

    def test_unknown_and_out_of_batch_selectors_fail(self):
        for selection in ({"gene_ids": ["g3"]}, {"gene_ids": ["unknown"]}, {"unit_ids": ["u3"]},
                          {"fields": ["unknown"]}, {"gene_ids": ["g1"], "fields": ["/genes/0/unknown"]},
                          {"unit_ids": ["u2"], "gene_ids": ["g1"]}):
            with self.subTest(selection=selection), self.assertRaises(ValueError):
                self.read(**selection)

    def test_tampered_bundle_index_material_and_roster_fail(self):
        for path in (self.bundle_path, self.index_path,
                     self.index_path.parent / "b01/material.full.json", self.index_path.parent / "b01/roster.json"):
            raw = path.read_bytes()
            path.write_bytes(raw + b" ")
            with self.subTest(path=path), self.assertRaises(ValueError):
                self.read()
            path.write_bytes(raw)

    def test_rehashed_wrong_case_roster_or_sidecar_fail(self):
        original = self.index_path.read_bytes()
        index = json.loads(original)
        for key, value in (("case_id", "wrong"),):
            index[key] = value
            self.index_path.write_text(json.dumps(index))
            self.index_sha = digest(self.index_path.read_bytes())
            with self.assertRaises(ValueError):
                self.read()
        self.index_path.write_bytes(original)
        index = json.loads(original)
        index["batches"][0]["gene_ids"] = ["g9"]
        self.index_path.write_text(json.dumps(index))
        self.index_sha = digest(self.index_path.read_bytes())
        with self.assertRaises(ValueError):
            self.read()
        self.index_path.write_bytes(original)
        sidecar = self.index_path.parent / "b01/material.full.json"
        data = json.loads(sidecar.read_text())
        data["units"][0]["genes"][0]["reference_function"] = "changed"
        sidecar.write_text(json.dumps(data))
        index = json.loads(original)
        index["batches"][0]["material_full_sha256"] = digest(sidecar.read_bytes())
        self.index_path.write_text(json.dumps(index))
        self.index_sha = digest(self.index_path.read_bytes())
        with self.assertRaises(ValueError):
            self.read()

    def test_no_other_batch_files_are_opened(self):
        other = self.index_path.parent / "b02/material.full.json"
        other.unlink()
        self.assertEqual("selected", self.read(gene_ids=["g1"], fields=["reference_function"])["state"])

    def test_science_v2_and_read_only_inputs(self):
        light.light_prepare(str(self.bundle_path), self.bundle_sha, str(self.root / "v2"),
                            batch_units=2, output_contract="science-v2")
        index = self.root / "v2/light/index.json"
        tracked = [self.bundle_path, index, index.parent / "b01/material.full.json", index.parent / "b01/roster.json"]
        before = {path: path.read_bytes() for path in tracked}
        packet = material_read(str(self.bundle_path), self.bundle_sha, str(index), digest(index.read_bytes()),
                               batch_id="b01", gene_ids=["g1"], fields=["reference_function"])
        self.assertEqual(self.long, packet["selections"][0]["value"])
        self.assertEqual(before, {path: path.read_bytes() for path in tracked})

    def test_rehashed_wrong_roster_is_rejected(self):
        index = json.loads(self.index_path.read_text())
        roster = self.index_path.parent / "b01/roster.json"
        roster.write_text(json.dumps({"g1": {"wrong": "value"}}))
        index["batches"][0]["roster_sha256"] = digest(roster.read_bytes())
        self.index_path.write_text(json.dumps(index))
        self.index_sha = digest(self.index_path.read_bytes())
        with self.assertRaises(ValueError):
            self.read()


class SavedSourceReadTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / "source.json"

    def write(self, text):
        self.path.write_bytes(text.encode("utf-8"))
        return digest(self.path.read_bytes())

    def test_json_metadata_abstract_and_distinct_observations(self):
        facts = {"title": "Original title", "source_id": "PMID:42", "abstract": "x" * 900,
                 "observations": [{"gene": "A", "trait": "leaf"}, {"gene": "B", "trait": "root"}]}
        sha = self.write(json.dumps(facts))
        packet = source_read(str(self.path), sha, pointers=["/title", "/source_id", "/abstract", "/observations/1"])
        self.assertEqual([facts[key] for key in ("title", "source_id", "abstract")] + [facts["observations"][1]],
                         [item["value"] for item in packet["selections"]])
        self.assertNotIn(facts["observations"][0], [item["value"] for item in packet["selections"]])
        self.assertEqual("index_only", source_read(str(self.path), sha)["state"])

    def test_inclusive_text_range_keeps_exact_content(self):
        sha = self.write("title\n" + "original fact " * 60 + "\nnot requested\n")
        packet = source_read(str(self.path), sha, line_start=2, line_end=2)
        self.assertEqual("original fact " * 60 + "\n", packet["selections"][0]["value"])
        self.assertEqual(["lines:2-2"], packet["selected_pointers"])
        self.assertEqual("narrow_selection_required", source_read(str(self.path), sha, line_start=2, line_end=2, max_chars=400)["state"])

    def test_bad_hash_unknown_pointer_and_invalid_range_fail(self):
        sha = self.write('{"abstract":"fact"}\n')
        for kwargs in ({"pointers": ["/missing"]}, {"pointers": ["/abstract~2"]},
                       {"line_start": 0, "line_end": 1}, {"line_start": 1, "line_end": 2},
                       {"pointers": ["/abstract"], "line_start": 1, "line_end": 1}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                source_read(str(self.path), sha, **kwargs)
        with self.assertRaises(ValueError):
            source_read(str(self.path), "0" * 64, pointers=["/abstract"])

    def test_absolute_cli_is_reachable_from_unrelated_working_directory(self):
        sha = self.write('{"title":"Original title","abstract":"Exact original fact"}')
        script = Path(__file__).resolve().parents[1] / "scripts/agentic_route_b.py"
        result = subprocess.run(
            [sys.executable, str(script), "light-read-source", "--source", str(self.path),
             "--source-sha256", sha, "--pointer", "/abstract"],
            cwd=self.path.parent, capture_output=True, text=True, encoding="utf-8")
        self.assertEqual(0, result.returncode, result.stderr)
        packet = json.loads(result.stdout)
        self.assertEqual("Exact original fact", packet["selections"][0]["value"])
        self.assertEqual(["/abstract"], packet["selected_pointers"])
