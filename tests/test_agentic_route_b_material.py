"""Readable transport preserves the entire closed-world semantic assignment."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from scripts.agentic_route_b import merge, prepare
from scripts.agentic_route_b_material import (
    DISPLAY_MAX_LINE_CHARS, _display_sidecar, _encoded, plan_reading_views,
    restore_display, restore_reading_shard,
)


def semantic_prompt(index, function="Complete function text", source="PMID:1", **extra):
    gene = f"gene-{index}"
    card = {"interval_gene_id": gene, "gene_id": "reference-1", "reference_gene": "ref-1",
            "reference_species": "Species one", "family_id": "same-family",
            "function_description": function, "source_id": source, "source_url": "https://example.org/record",
            "optional_null": None, "optional_empty": "", **extra}
    unit = {"unit_id": f"unit-{index}", "unit_type": "family", "unit_anchor": "same-family",
            "member_interval_genes": [gene], "request_state": "requested", "review_status": "not_evaluated",
            "stable_source_ids": [source], "genes": [{"gene_id": gene, "start": str(index),
                "reference_gene": "ref-1", "reference_function": function, "reference_source_id": source,
                "local_annotation_description": "copy-specific observation", "sequence_similarity_hits": [{"percent_identity": str(index)}]}],
            "reference_material_cards": [card, dict(card)],
            "source_provenance": [{"gene_id": gene, "reference_gene": "ref-1", "source_ids": [source]}],
            "literature_link_candidates": []}
    return {"case_id": "case-x", "prompt_id": f"prompt-{index}", "route_b_unit_id": f"unit-{index}",
            "gene_id": "__route_b_semantic_unit__", "system": "Complete instructions for every unit",
            "user": "Phenotype: test\nRoute: interval_first\n\n" + json.dumps(unit, ensure_ascii=False, sort_keys=True)}


class TestSharedRouteBMaterial(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.source = self.root / "source.jsonl"
        self.out = self.root / "prepared"

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, rows):
        self.source.write_text("".join(json.dumps(r, ensure_ascii=False) + "\n" for r in rows), encoding="utf-8")

    def restore(self, manifest):
        return [r for shard in manifest["reading_view"]["shards"] for r in restore_reading_shard(self.out, shard["file"])]

    def test_full_fields_duplicates_and_distinct_same_family_material(self):
        rows = [semantic_prompt(1), semantic_prompt(2), semantic_prompt(3, "Contradictory function", "PMID:2", additional_field=False)]
        self.write(rows)
        manifest = prepare(str(self.source), str(self.out), max_prompts=2)
        self.assertEqual(rows, self.restore(manifest))
        self.assertTrue(manifest["reading_view"]["exact_roundtrip_verified"])
        material = [json.loads(p.read_text(encoding="utf-8")) for p in (self.out / "shared-material").glob("*.json")]
        cards = [m["value"] for m in material if m["kind"] == "reference-material-cards"]
        self.assertEqual(2, len(cards))  # Identical cross-gene copies shared; different evidence retained.
        self.assertEqual({"PMID:1", "PMID:2"}, {c["source_id"] for c in cards})
        self.assertEqual(1, sum(m["kind"] == "system" for m in material))
        self.assertEqual(1, sum(m["kind"] == "user-prefix" for m in material))
        assignment = json.loads((self.out / manifest["reading_view"]["shards"][0]["file"]).read_text(encoding="utf-8"))
        bindings = assignment["prompts"][0]["unit_tables"]["reference_material_cards"]
        self.assertEqual(2, len(bindings))
        self.assertEqual("gene-1", bindings[0]["binding"]["interval_gene_id"])
        self.assertEqual("ref-1", bindings[0]["reference_handle"]["reference_gene"])
        self.assertEqual(bindings[0], bindings[1])
        self.assertEqual("copy-specific observation", assignment["prompts"][0]["unit_tables"]["genes"][0]["binding"]["local_annotation_description"])
        # Preparation preserves the source bytes, canonical rows and every assigned ID.
        self.assertEqual(hashlib.sha256(self.source.read_bytes()).hexdigest(), manifest["source_sha256"])
        canonical = [json.loads(line) for s in manifest["shards"] for line in (self.out / s["file"]).read_text(encoding="utf-8").splitlines()]
        self.assertEqual(rows, canonical)

    def test_display_roundtrip_handles_unicode_long_keys_and_literal_markers(self):
        marker = "__cqtl_display_string_chunks__"
        long_key = "長い/key~" * 500
        value = {"__display_key_0__": "original key", marker: ["original", "array"],
                 "nested": [{marker: ["not", "encoded"]}, None, True, 1, -0.5],
                 long_key: {"source": "生物🌱\n\\\"\u2028\u0085" * 900,
                            "controls": "\x00" * 1200}, "empty": ""}
        canonical = _encoded(value)
        data, descriptor = _display_sidecar(canonical, "shared-material/fixture.json")
        self.assertEqual(value, restore_display(data, canonical))
        self.assertTrue(all(len(line) <= DISPLAY_MAX_LINE_CHARS for line in data.decode().splitlines()))
        self.assertEqual(len(data.decode().splitlines()), descriptor["line_count"])
        self.assertEqual(hashlib.sha256(canonical).hexdigest(), descriptor["source_sha256"])
        displayed = json.loads(data)
        self.assertEqual("__display_key_1__", displayed["renamed_keys"][0]["display_key"])
        self.assertEqual(["original", "array"], displayed["value"][marker])
        self.assertEqual({marker: ["not", "encoded"]}, displayed["value"]["nested"][0])

    def test_canonical_assignment_and_material_hashes_do_not_change_for_display(self):
        # Captured from the same fixture at a71746b, before display sidecars.
        files, shards = plan_reading_views([semantic_prompt(1), semantic_prompt(2)], "case-x", "0" * 64, 2, 262144)
        self.assertEqual("f09bc8102d296e9cad9c4f055d8b33f2c225104cb16032dfd4c84ffcf5758b9e",
                         hashlib.sha256(files["reading/b0000/assignment.json"]).hexdigest())
        assignment = json.loads(files[shards[0]["file"]])
        for material in assignment["materials"]:
            self.assertEqual(material["sha256"], hashlib.sha256(files[material["file"]]).hexdigest())
        self.assertEqual(len(files[shards[0]["file"]]) + sum(m["bytes"] for m in assignment["materials"]),
                         shards[0]["total_reading_bytes"])

    def test_display_bounds_count_unicode_separator_escapes(self):
        for separator in ("\u0085", "\u2028", "\u2029"):
            with self.subTest(separator=ord(separator)):
                value = {separator * 200: [separator * 200, separator * 1200]}
                canonical = _encoded(value)
                data, _ = _display_sidecar(canonical, "shared-material/separators.json")
                self.assertEqual(value, restore_display(data, canonical))
                self.assertLessEqual(max(map(len, data.decode().splitlines())), DISPLAY_MAX_LINE_CHARS)
                self.assertNotIn(separator, data.decode())

    def test_direct_cli_long_line_display_reaches_every_material_and_preserves_source(self):
        rows = [semantic_prompt(1, "complete source 日本語🌱 " * 1000), semantic_prompt(2)]
        rows[0]["system"] = "Original system instructions. " * 1000
        rows[1]["user"] = "Unknown layout retains every character. 日本語 " * 1000
        self.write(rows)
        before = self.source.read_bytes()
        script = Path(__file__).resolve().parents[1] / "scripts" / "agentic_route_b.py"
        subprocess.run([sys.executable, str(script), "prepare", "--prompts", str(self.source),
                        "--out", str(self.out), "--max-prompts", "1"],
                       check=True, capture_output=True, text=True, encoding="utf-8")
        manifest = json.loads((self.out / "prepare_manifest.json").read_bytes())
        self.assertEqual(before, self.source.read_bytes())
        self.assertEqual(rows, self.restore(manifest))
        all_displays = set()
        for shard in manifest["reading_view"]["shards"]:
            view = shard["display_view"]
            assignment = json.loads((self.out / shard["file"]).read_bytes())
            self.assertEqual({m["file"] for m in assignment["materials"]},
                             {d["source_file"] for d in view["materials"]})
            for descriptor in [view["assignment"], *view["materials"]]:
                original = (self.out / descriptor["source_file"]).read_bytes()
                data = (self.out / descriptor["file"]).read_bytes()
                self.assertEqual(json.loads(original), restore_display(data, original))
                self.assertLessEqual(max(map(len, data.decode().splitlines())), DISPLAY_MAX_LINE_CHARS)
                all_displays.add(descriptor["file"])
            task = (self.out / shard["worker_task"]["file"]).read_text(encoding="utf-8")
            self.assertIn(view["assignment"]["file"], task)
            self.assertIn("display/<materials[].file>", task)
            self.assertLessEqual(max(map(len, task.splitlines())), DISPLAY_MAX_LINE_CHARS)
        self.assertEqual(all_displays, {str(p.relative_to(self.out)).replace("\\", "/")
                                       for p in (self.out / "display").rglob("*.json")})

    def test_display_tampering_missing_links_and_legacy_generations(self):
        self.write([semantic_prompt(1)])
        manifest = prepare(str(self.source), str(self.out))
        shard = manifest["reading_view"]["shards"][0]
        display = shard["display_view"]["assignment"]
        path = self.out / display["file"]
        original = path.read_bytes()
        path.write_bytes(original + b" ")
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            self.restore(manifest)
        path.write_bytes(original)
        material = self.out / shard["display_view"]["materials"][0]["file"]
        saved = material.read_bytes()
        material.unlink()
        with self.assertRaises(FileNotFoundError):
            self.restore(manifest)
        material.write_bytes(saved)
        changed = json.loads(original)
        changed["value"]["case_id"] = "changed-science"
        data = (json.dumps(changed, ensure_ascii=False, indent=2) + "\n").encode()
        path.write_bytes(data)
        display.update(sha256=hashlib.sha256(data).hexdigest(), bytes=len(data),
                       line_count=len(data.decode().splitlines()),
                       max_line_characters=max(map(len, data.decode().splitlines())))
        (self.out / "prepare_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "display differs from canonical"):
            self.restore(manifest)
        shard.pop("display_view")
        shard.pop("worker_task")
        (self.out / "prepare_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        self.assertEqual([semantic_prompt(1)], self.restore(manifest))

    def test_byte_bound_counts_all_linked_material_and_splits_without_loss(self):
        rows = [semantic_prompt(i, "full " + str(i) + " " + "x" * 4000) for i in range(4)]
        self.write(rows)
        _, one = plan_reading_views(rows, "case-x", hashlib.sha256(self.source.read_bytes()).hexdigest(), 1, 100000)
        bound = max(s["total_reading_bytes"] for s in one)
        manifest = prepare(str(self.source), str(self.out), max_prompts=4, max_material_bytes=bound)
        self.assertEqual(4, len(manifest["reading_view"]["shards"]))
        self.assertTrue(all(s["total_reading_bytes"] <= bound for s in manifest["reading_view"]["shards"]))
        self.assertEqual(rows, self.restore(manifest))
        self.assertEqual(1, len(manifest["shards"]))  # Original max-prompts contract unchanged.

    def test_refrozen_paper_shared_across_units_preserves_local_bindings(self):
        rows = [semantic_prompt(i) for i in range(1, 4)]
        for index, row in enumerate(rows, 1):
            prefix, encoded = row["user"].split("\n\n", 1)
            unit = json.loads(encoded)
            paper = {"unit_id": unit["unit_id"], "reference_gene": "ref-1",
                     "reference_species": "Species one", "source_id": "PMID:1",
                     "source_text": "Complete shared paper text. " * 300}
            if index == 3:
                paper["source_id"] = "PMID:2"
            unit["literature_link_candidates"] = [paper, dict(paper)]
            unit["source_provenance"] = [{key: value for key, value in paper.items() if key != "source_text"}]
            row["user"] = prefix + "\n\n" + json.dumps(unit, ensure_ascii=False, sort_keys=True)
        self.write(rows)
        manifest = prepare(str(self.source), str(self.out))
        self.assertEqual(rows, self.restore(manifest))
        assignment = json.loads((self.out / manifest["reading_view"]["shards"][0]["file"]).read_bytes())
        for field in ("literature_link_candidates", "source_provenance"):
            links = [prompt["unit_tables"][field][0] for prompt in assignment["prompts"]]
            self.assertEqual(links[0]["material_id"], links[1]["material_id"])
            self.assertNotEqual(links[0]["material_id"], links[2]["material_id"])
            self.assertEqual(["unit-1", "unit-2", "unit-3"], [link["binding"]["unit_id"] for link in links])
        for prompt in assignment["prompts"]:
            papers = prompt["unit_tables"]["literature_link_candidates"]
            self.assertEqual(2, len(papers))
            self.assertEqual(papers[0], papers[1])
        # Even a self-consistent assignment hash cannot authorize a changed binding.
        path = self.out / manifest["reading_view"]["shards"][0]["file"]
        assignment["prompts"][0]["unit_tables"]["literature_link_candidates"][0]["binding"]["unit_id"] = "unit-2"
        changed = json.dumps(assignment, ensure_ascii=False).encode("utf-8")
        path.write_bytes(changed)
        entry = manifest["reading_view"]["shards"][0]
        entry["sha256"] = hashlib.sha256(changed).hexdigest()
        entry["total_reading_bytes"] = len(changed) + entry["material_bytes"]
        (self.out / "prepare_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "restored prompt hash mismatch"):
            self.restore(manifest)

    def test_invalid_bound_fails_before_output_creation(self):
        self.write([semantic_prompt(1)])
        for bound in (0, -1):
            with self.subTest(bound=bound), self.assertRaisesRegex(ValueError, "positive"):
                prepare(str(self.source), str(self.out), max_material_bytes=bound)
            self.assertFalse(self.out.exists())

    def test_launch_task_reaches_linked_schema_and_full_material_from_root(self):
        rows = [semantic_prompt(1), semantic_prompt(2)]
        self.write(rows)
        manifest = prepare(str(self.source), str(self.out), max_prompts=1)
        for entry in manifest["reading_view"]["shards"]:
            task = entry["worker_task"]
            data = (self.out / task["file"]).read_bytes()
            self.assertEqual(task["bytes"], len(data))
            self.assertEqual(task["sha256"], hashlib.sha256(data).hexdigest())
            launch = data.decode("utf-8")
            self.assertIn(entry["file"], launch)
            self.assertIn("relative to the prepare output root", launch)
            self.assertIn("including system and user_prefix", launch)
            self.assertIn('"case_id":"case-x","responses":[...]', launch)
            self.assertIn("raw_text", launch)
            self.assertIn('status="not_evaluated"', launch)
            assignment = json.loads((self.out / entry["file"]).read_bytes())
            for link in assignment["materials"]:
                self.assertTrue((self.out / link["file"]).is_file())
                self.assertFalse((self.out / Path(entry["file"]).parent / link["file"]).exists())
        self.assertEqual(rows, self.restore(manifest))

    def test_launch_task_tampering_rejected_and_legacy_generation_replays(self):
        self.write([semantic_prompt(1)])
        manifest = prepare(str(self.source), str(self.out))
        entry = manifest["reading_view"]["shards"][0]
        task = entry["worker_task"]
        path = self.out / task["file"]
        path.write_bytes(b"Read no other local paths.\n")
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            self.restore(manifest)
        task.update(sha256=hashlib.sha256(path.read_bytes()).hexdigest(), bytes=path.stat().st_size + 1)
        (self.out / "prepare_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "worker task"):
            self.restore(manifest)
        entry.pop("worker_task")
        (self.out / "prepare_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        self.assertEqual([semantic_prompt(1)], self.restore(manifest))

    def test_oversized_unit_keeps_full_material_without_blocking_healthy_units(self):
        rows = [semantic_prompt(1), semantic_prompt(2, "x" * 100000), semantic_prompt(3)]
        self.write(rows)
        manifest = prepare(str(self.source), str(self.out), max_material_bytes=20000)
        shards = manifest["reading_view"]["shards"]
        self.assertEqual(["bounded", "oversized_full_material", "bounded"], [s["reading_state"] for s in shards])
        self.assertEqual(["prompt-2"], shards[1]["prompt_ids"])
        self.assertGreater(shards[1]["over_bound_bytes"], 0)
        self.assertEqual(1, manifest["reading_view"]["oversized_unit_count"])
        self.assertEqual(rows, self.restore(manifest))

    def test_unknown_user_layout_and_missing_fields_remain_exact(self):
        rows = [{"case_id": "case-x", "prompt_id": "p1", "route_b_unit_id": "u1", "user": "unformatted 日本語 {literal}"},
                {"case_id": "case-x", "prompt_id": "p2", "route_b_unit_id": "u2", "system": None, "user": None},
                {"case_id": "case-x", "prompt_id": "p3", "route_b_unit_id": "u3"}]
        self.write(rows)
        self.assertEqual(rows, self.restore(prepare(str(self.source), str(self.out))))

    def test_material_and_canonical_hash_tampering_are_rejected(self):
        self.write([semantic_prompt(1)])
        manifest = prepare(str(self.source), str(self.out))
        shard = manifest["reading_view"]["shards"][0]["file"]
        material = next((self.out / "shared-material").glob("*.json"))
        original = material.read_bytes()
        material.write_bytes(original + b" ")
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            restore_reading_shard(self.out, shard)
        material.write_bytes(original)
        canonical = self.out / manifest["shards"][0]["file"]
        canonical.write_bytes(canonical.read_bytes() + b" ")
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            restore_reading_shard(self.out, shard)

    def test_assignment_case_id_source_id_and_binding_tampering_are_rejected(self):
        self.write([semantic_prompt(1)])
        manifest = prepare(str(self.source), str(self.out))
        shard = manifest["reading_view"]["shards"][0]["file"]
        path = self.out / shard
        original = path.read_bytes()
        cases = (
            (lambda a: a.update(case_id="another-case"), "source/case mismatch"),
            (lambda a: a.update(source_sha256="0" * 64), "source/case mismatch"),
            (lambda a: a["prompts"][0]["prompt"].update(prompt_id="unknown-prompt"), "restored prompt hash mismatch"),
            (lambda a: a["prompts"][0]["unit_tables"]["reference_material_cards"][0]["binding"].update(interval_gene_id="another-gene"), "restored prompt hash mismatch"),
            (lambda a: a["prompts"][0]["unit_tables"]["reference_material_cards"][0]["reference_handle"].update(source_id="PMID:999"), "reference/source handle mismatch"),
        )
        for mutate, error in cases:
            with self.subTest(error=error):
                assignment = json.loads(original)
                mutate(assignment)
                changed = json.dumps(assignment, ensure_ascii=False).encode("utf-8")
                path.write_bytes(changed)
                entry = manifest["reading_view"]["shards"][0]
                entry["sha256"] = hashlib.sha256(changed).hexdigest()
                entry["total_reading_bytes"] = len(changed) + entry["material_bytes"]
                (self.out / "prepare_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
                with self.assertRaisesRegex(ValueError, error):
                    restore_reading_shard(self.out, shard)

    def test_empty_roster_preserves_case_and_no_assignments(self):
        self.write([])
        manifest = prepare(str(self.source), str(self.out), case_id="empty-case")
        self.assertEqual("empty-case", manifest["case_id"])
        self.assertEqual([], manifest["reading_view"]["shards"])
        self.assertEqual([], self.restore(manifest))

    def test_existing_direct_cli_reaches_reading_view(self):
        rows = [semantic_prompt(1), semantic_prompt(2)]
        self.write(rows)
        script = Path(__file__).resolve().parents[1] / "scripts" / "agentic_route_b.py"
        run = subprocess.run([sys.executable, str(script), "prepare", "--prompts", str(self.source),
                              "--out", str(self.out), "--max-prompts", "1", "--max-material-bytes", "100000"],
                             check=True, capture_output=True, text=True, encoding="utf-8")
        summary = json.loads(run.stdout)
        self.assertEqual(2, summary["reading_view"]["shard_count"])
        self.assertTrue(summary["reading_view"]["exact_roundtrip_verified"])
        manifest = json.loads((self.out / "prepare_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(rows, self.restore(manifest))


    def test_opt_in_affinity_colocates_identical_reference_cards(self):
        # Interleaved same-function evidence: transport may co-locate, not merge,
        # original distinct gene/copy/source choices.
        rows = [
            semantic_prompt(1, "Function alpha", "PMID:11"),
            semantic_prompt(2, "Function beta", "PMID:22"),
            semantic_prompt(3, "Function alpha", "PMID:11"),
            semantic_prompt(4, "Function beta", "PMID:22"),
        ]
        self.write(rows)
        baseline_root = self.root / "baseline"
        baseline = prepare(str(self.source), str(baseline_root), max_prompts=2)
        optimized = prepare(
            str(self.source), str(self.out), max_prompts=2,
            prefer_shared_material=True,
        )

        assignments = optimized["reading_view"]["shards"]
        self.assertEqual(
            [["prompt-1", "prompt-3"], ["prompt-2", "prompt-4"]],
            [shard["prompt_ids"] for shard in assignments],
        )
        self.assertTrue(optimized["reading_view"]["exact_roundtrip_verified"])
        self.assertEqual(
            "exact_shared_reference_material_affinity_opt_in",
            optimized["reading_view"]["packing_strategy"],
        )
        restored = self.restore(optimized)
        self.assertEqual(
            {row["prompt_id"]: row for row in rows},
            {row["prompt_id"]: row for row in restored},
        )
        self.assertEqual(len(rows), len(restored))
        # The frozen canonical prompt shards are unaffected byte for byte.
        for base, packed in zip(baseline["shards"], optimized["shards"]):
            self.assertEqual(
                (baseline_root / base["file"]).read_bytes(),
                (self.out / packed["file"]).read_bytes(),
            )
        baseline_bytes = sum(s["total_reading_bytes"] for s in baseline["reading_view"]["shards"])
        packed_bytes = sum(s["total_reading_bytes"] for s in assignments)
        self.assertLess(packed_bytes, baseline_bytes)
        self.assertEqual(packed_bytes, optimized["reading_view"]["total_reading_bytes"])
        self.assertEqual(
            {row["route_b_unit_id"] for row in rows},
            {row["route_b_unit_id"] for row in restored},
        )
        # Responses may arrive in reader order; canonical merge is still
        # prompt-ID keyed and reports the complete original roster.
        responses = self.root / "affinity-responses.json"
        responses.write_text(json.dumps({
            "case_id": "case-x", "executor": "test", "model": "stub",
            "responses": [
                {"prompt_id": row["prompt_id"], "route_b_unit_id": row["route_b_unit_id"],
                 "raw_text": '{"claims":[]}'}
                for row in reversed(restored)
            ],
        }), encoding="utf-8")
        joined = self.root / "affinity-joined.json"
        result = merge(str(self.source), [str(responses)], str(joined))
        self.assertEqual(4, result["completed_count"])
        self.assertEqual("complete", result["completion_status"])
        self.assertEqual(
            [row["prompt_id"] for row in rows],
            [row["prompt_id"] for row in json.loads(joined.read_text(encoding="utf-8"))["candidates"]],
        )

    def test_opt_in_does_not_fold_same_family_with_distinct_sources(self):
        rows = [
            semantic_prompt(1, "Function alpha", "PMID:11"),
            semantic_prompt(2, "Contradictory function", "PMID:22"),
            {"case_id": "case-x", "prompt_id": "p-inline", "route_b_unit_id": "u-inline",
             "user": "unknown unformatted layout 日本語"},
        ]
        self.write(rows)
        optimized = prepare(
            str(self.source), str(self.out), max_prompts=2,
            prefer_shared_material=True,
        )
        self.assertEqual(rows, self.restore(optimized))
        self.assertEqual(3, optimized["prompt_count"])
        self.assertEqual(2, len(optimized["reading_view"]["shards"]))

    def test_opt_in_does_not_drop_oversized_copy_or_material(self):
        rows = [
            semantic_prompt(1, "Function alpha", "PMID:11"),
            semantic_prompt(2, "Large copy-specific " + "x" * 100000, "PMID:22"),
            semantic_prompt(3, "Function alpha", "PMID:11"),
        ]
        self.write(rows)
        optimized = prepare(
            str(self.source), str(self.out), max_prompts=3,
            max_material_bytes=20000, prefer_shared_material=True,
        )
        self.assertEqual(
            {row["prompt_id"]: row for row in rows},
            {row["prompt_id"]: row for row in self.restore(optimized)},
        )
        self.assertEqual(1, optimized["reading_view"]["oversized_unit_count"])
        self.assertEqual(
            ["prompt-2"],
            next(s["prompt_ids"] for s in optimized["reading_view"]["shards"]
                 if s["reading_state"] == "oversized_full_material"),
        )

    def test_opt_in_falls_back_when_small_cards_break_large_reference_sharing(self):
        rows = []
        for index in range(4):
            prompt = semantic_prompt(index + 1, "short card", f"PMID:{index % 2}")
            prefix, payload = prompt["user"].split("\n\n", 1)
            unit = json.loads(payload)
            unit["genes"][0]["reference_function"] = ("X" if index < 2 else "Y") * 20000
            unit["genes"][0]["reference_source_id"] = f"REF:{index // 2}"
            prompt["user"] = prefix + "\n\n" + json.dumps(unit, ensure_ascii=False, sort_keys=True)
            rows.append(prompt)

        # Canonical pairs share the large records. Card affinity alone would
        # split those pairs, and at the smaller bound double the worker count.
        for bound in (262144, 30000):
            with self.subTest(bound=bound):
                baseline = plan_reading_views(rows, "case-x", "0" * 64, 2, bound)
                optimized = plan_reading_views(
                    rows, "case-x", "0" * 64, 2, bound,
                    prefer_shared_material=True,
                )
                self.assertEqual(baseline, optimized)
                self.assertEqual(2, len(optimized[1]))

        self.write(rows)
        subprocess.run([
            sys.executable, "scripts/agentic_route_b.py", "prepare",
            "--prompts", str(self.source), "--out", str(self.out),
            "--max-prompts", "2", "--max-material-bytes", "30000",
            "--prefer-shared-material",
        ], check=True, capture_output=True, text=True)
        manifest = json.loads((self.out / "prepare_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual("canonical_order_no_safe_affinity_gain", manifest["reading_view"]["packing_strategy"])
        self.assertEqual(rows, self.restore(manifest))

    def test_opt_in_falls_back_for_overlapping_reference_cards(self):
        rows = []
        for index, names in enumerate((("A", "X"), ("X",), ("A", "Y"), ("A", "Y")), 1):
            prompt = semantic_prompt(index, "shared reference", "source")
            prefix, payload = prompt["user"].split("\n\n", 1)
            unit = json.loads(payload)
            unit["genes"] = []
            unit["source_provenance"] = []
            unit["reference_material_cards"] = [
                {"interval_gene_id": f"gene-{index}", "reference_gene": name,
                 "source_id": f"source-{name}", "family_id": "same-family",
                 "function_description": name * (2500 if name == "A" else 3000)}
                for name in names
            ]
            unit["stable_source_ids"] = [f"source-{name}" for name in names]
            prompt["user"] = prefix + "\n\n" + json.dumps(unit, ensure_ascii=False, sort_keys=True)
            rows.append(prompt)
        baseline = plan_reading_views(rows, "case-x", "0" * 64, 2, 262144)
        optimized = plan_reading_views(rows, "case-x", "0" * 64, 2, 262144, prefer_shared_material=True)
        self.assertEqual(baseline, optimized)


if __name__ == "__main__":
    unittest.main()
