"""Batched fresh-context Route B research (#1004): splitting and merge accounting."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from scripts import agentic_route_b as route_b
from scripts.agentic_route_b_initial import summarize_initial
from scripts.agentic_route_b_research import (
    CONTRACT,
    FOLLOWUP_TASK,
    TASK,
    merge_initial_batches,
    plan_research_view,
)
from tests.test_agentic_route_b_research import seal_prompts

ROOT = Path(__file__).resolve().parents[1]


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def unit(uid: str, gid: str, reference_gene: str = "", requested: bool = True) -> dict:
    gene = {"gene_id": gid}
    record = {"unit_id": uid, "member_interval_genes": [gid], "identity": {},
              "request_state": "requested" if requested else "not_requested_no_reviewable_material",
              "source_material_class": "annotation_context"}
    if reference_gene:
        gene["reference_gene"] = reference_gene
        gene["reference_species"] = "Arabidopsis thaliana"
        gene["reference_function"] = f"{reference_gene} documented process"
        record["reference_material_cards"] = [
            {"interval_gene_id": gid, "reference_gene": reference_gene,
             "function_description": f"{reference_gene} documented process"}]
    record["genes"] = [gene]
    return record


def fixture() -> tuple[dict, list[dict]]:
    units = [
        unit("u1", "g1", "REF-A"), unit("u2", "g2", "REF-A"), unit("u3", "g3", "REF-A"),
        unit("u4", "g4", "REF-B"), unit("u5", "g5", "REF-B"), unit("u6", "g6"),
        unit("u7", "g7", "", requested=False),
    ]
    prompts = [{"case_id": "fixture", "prompt_id": f"route-b-semantic-{u['unit_id']}",
                "route_b_unit_id": u["unit_id"], "system": "TEST ONLY",
                "user": "Phenotype: test\n\n" + json.dumps(u, ensure_ascii=False, sort_keys=True)}
               for u in units if u["request_state"] == "requested"]
    bundle = {"case_id": "fixture", "route": "interval_first", "stage": "post_b0_semantic_bundle",
              "phenotype": "test", "units": units, "interval_gene_count": 7, "semantic_unit_count": 7}
    seal_prompts(bundle, prompts)
    return bundle, prompts


class BatchedResearchViewTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.bundle, self.prompts = fixture()

    def save(self) -> str:
        self.source = self.root / "bundle.json"
        self.source.write_text(json.dumps(self.bundle, ensure_ascii=False), encoding="utf-8")
        self.prompts_path = self.root / "prompts.jsonl"
        self.prompts_path.write_text("".join(json.dumps(p, ensure_ascii=False) + "\n" for p in self.prompts),
                                     encoding="utf-8")
        return digest(self.source.read_bytes())

    def plan(self, **kwargs):
        sha = self.save()
        return plan_research_view(str(self.source), sha, self.prompts, "fixture", 262144, **kwargs)

    @staticmethod
    def records(files: dict[str, bytes], prefix: str) -> list[dict]:
        return [json.loads(line) for name, data in sorted(files.items())
                if name.startswith(prefix) and name.endswith(".jsonl") for line in data.splitlines()]

    def test_default_output_is_byte_identical_to_no_batching(self):
        default, default_meta = self.plan()
        explicit, explicit_meta = self.plan(batch_units=0)
        self.assertEqual(default, explicit)
        self.assertEqual(default_meta, explicit_meta)
        self.assertFalse(any("batches/" in name for name in default))
        self.assertNotIn("research_batches", json.loads(default["research/manifest.json"]))
        self.assertEqual(TASK.encode("utf-8"), default["research/task.md"])

    def test_every_requested_unit_appears_once_and_batches_stay_within_the_limit(self):
        files, _ = self.plan(batch_units=2)
        batches = json.loads(files["research/batches/index.json"])["batches"]
        seen = [uid for batch in batches for uid in batch["unit_ids"]]
        self.assertEqual(["u1", "u2", "u3", "u4", "u5", "u6"], seen)
        self.assertTrue(all(batch["requested_unit_count"] <= 2 for batch in batches))
        self.assertNotIn("u7", seen)
        rows = self.records(files, "research/batches/")
        self.assertEqual(["u1", "u2", "u3", "u4", "u5", "u6"], [row["unit_id"] for row in rows])

    def test_units_sharing_a_reference_endpoint_share_a_batch(self):
        files, _ = self.plan(batch_units=3)
        batches = json.loads(files["research/batches/index.json"])["batches"]
        by_unit = {uid: batch["batch_id"] for batch in batches for uid in batch["unit_ids"]}
        self.assertEqual(by_unit["u1"], by_unit["u2"])
        self.assertEqual(by_unit["u1"], by_unit["u3"])
        self.assertEqual(by_unit["u4"], by_unit["u5"])
        self.assertEqual(["b001", "b002"], [batch["batch_id"] for batch in batches])
        self.assertEqual([], json.loads(files["research/batches/index.json"])["split_groups"])

    def test_oversized_group_is_split_and_the_split_is_recorded(self):
        files, _ = self.plan(batch_units=2)
        index = json.loads(files["research/batches/index.json"])
        batches = index["batches"]
        self.assertEqual([["u1", "u2"], ["u3"], ["u4", "u5"], ["u6"]],
                         [batch["unit_ids"] for batch in batches])
        self.assertEqual(1, len(index["split_groups"]))
        split = index["split_groups"][0]
        self.assertEqual(["u1", "u2", "u3"], split["group_unit_ids"])
        self.assertEqual(2, split["chunk_count"])
        self.assertEqual(["b001", "b002"], split["batch_ids"])
        self.assertEqual({"chunk_index": 0, "chunk_count": 2, "group_unit_ids": ["u1", "u2", "u3"]},
                         batches[0]["split"])
        self.assertEqual(1, batches[1]["split"]["chunk_index"])

    def test_batch_task_manifest_and_material_bind_only_that_batch(self):
        files, metadata = self.plan(batch_units=3)
        index = json.loads(files["research/batches/index.json"])
        batch = index["batches"][0]
        bid = batch["batch_id"]
        self.assertEqual("b001", bid)
        task = files["research/batches/b001/task.md"].decode("utf-8")
        self.assertTrue(task.startswith(TASK))
        for text in (bid, "新しい文脈", "research/batches/b001/manifest.json",
                     "research/batches/b001/initial_assessments.json", "b001:", "merge-initial",
                     "他のbatch", "このbatchの外のgeneについて記録を書かない"):
            self.assertIn(text, task)
        self.assertEqual(digest(task.encode("utf-8")), batch["task"]["sha256"])
        manifest = json.loads(files["research/batches/b001/manifest.json"])
        self.assertEqual(["u1", "u2", "u3"], manifest["unit_ids"])
        self.assertEqual(["g1", "g2", "g3"], manifest["gene_ids"])
        self.assertEqual("b001:", manifest["ledger"]["id_prefix"])
        self.assertEqual("batch_only", manifest["delivery_scope"])
        self.assertNotIn("g7", manifest["gene_ids"])
        self.assertNotIn("g4", json.dumps(manifest))
        page_files = [page["file"] for page in manifest["pages"]]
        self.assertTrue(all(name.startswith("research/batches/b001/") for name in page_files))
        self.assertFalse(any("research/initial-material" in name for name in page_files))
        material_files = [name for name in files if name.startswith("research/batches/b001/initial-material")]
        self.assertTrue(material_files)
        for name in material_files:
            for group in json.loads(files[name]):
                for binding in group["bindings"]:
                    self.assertIn(binding["unit_id"], manifest["unit_ids"])
                    self.assertTrue(binding["source_pointer"].startswith(("/units/0/", "/units/1/", "/units/2/")))
        self.assertEqual(3, manifest["initial_material"]["gene_count"])
        self.assertIn("research_batches", json.loads(files["research/manifest.json"]))
        self.assertEqual(batch["manifest"]["sha256"], digest(files["research/batches/b001/manifest.json"]))
        self.assertEqual(2, metadata["research_batches"]["batch_count"])
        # Global pages and task stay as they are without batching.
        default, _ = self.plan()
        self.assertEqual(default["research/task.md"], files["research/task.md"])
        self.assertEqual(default["research/index-0000.jsonl"], files["research/index-0000.jsonl"])

    def test_followup_plan_cannot_be_combined_with_batching(self):
        with self.assertRaisesRegex(ValueError, "initial phase"):
            self.plan(batch_units=2, followup_plan="plan.json", followup_plan_sha256="0" * 64)


class MergeInitialTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.bundle, self.prompts = fixture()
        self.source = self.root / "bundle.json"
        self.source.write_text(json.dumps(self.bundle, ensure_ascii=False), encoding="utf-8")
        self.sha = digest(self.source.read_bytes())
        self.prompts_path = self.root / "prompts.jsonl"
        self.prompts_path.write_text("".join(json.dumps(p, ensure_ascii=False) + "\n" for p in self.prompts),
                                     encoding="utf-8")
        self.view, _ = plan_research_view(str(self.source), self.sha, self.prompts, "fixture", 262144, batch_units=3)
        self.out = self.root / "view"
        for name, data in self.view.items():
            target = self.out / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
        self.index = self.out / "research/batches/index.json"
        self.index_sha = digest(self.index.read_bytes())

    def ledger(self, bid: str, assessments: list[dict], shared: list[dict] | None = None,
               screens: list[dict] | None = None, gene_rows: list[dict] | None = None,
               bundle_sha256: str | None = None) -> dict:
        return {"bundle_sha256": bundle_sha256 or self.sha, "assessments": assessments,
                "shared_assessments": shared or [],
                "research_screening": {"reference_screens": screens or [],
                                       "gene_assessments": gene_rows or []}}

    def write(self, bid: str, payload: dict) -> None:
        target = self.out / f"research/batches/{bid}/initial_assessments.json"
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")

    def happy(self) -> None:
        self.write("b001", self.ledger(
            "b001",
            [{"assessment_id": "b001:A1", "gene_id": "g1", "unit_id": "u1",
              "state": "assessed_no_specific_bridge", "reason": "Original judgment for this batch.",
              "used_material_pointers": [], "assessed_fields": []},
             {"assessment_id": "b001:A2", "gene_id": "g2", "unit_id": "u2",
              "state": "not_assessed", "reason": "Batch scope left this copy unread.",
              "assessed_fields": []}],
            [],
            [{"screen_id": "b001:S1", "reference_gene": "REF-A", "reference_species": "Arabidopsis thaliana",
              "member_gene_ids": ["g1", "g2"], "state": "no_concrete_upper_lead",
              "queries": ["REF-A mutant phenotype"], "retrieval_state": "completed",
              "checked_species_scopes": ["focal_species", "other_plants"], "screened_source_ids": [],
              "reason": "No concrete upper lead in this actual inspected scope."}],
            [{"gene_id": "g1", "unit_id": "u1", "provisional_tier": "T4",
              "tier_reason": "A provisional planning entry.", "screen_ids": ["b001:S1"]}]))
        self.write("b002", self.ledger(
            "b002",
            [{"assessment_id": "b002:A1", "gene_id": "g4", "unit_id": "u4",
              "state": "assessed_concrete_lead", "reason": "Original copy-specific lead.",
              "used_material_pointers": [], "assessed_fields": []}],
            [{"assessment_id": "b002:J1", "material_id": "", "member_gene_ids": ["g6"],
              "state": "not_assessed", "reason": "No functional context in the supplied input for this copy.",
              "assessed_fields": []}],
            [{"screen_id": "b002:S1", "reference_gene": "REF-B", "reference_species": "Arabidopsis thaliana",
              "member_gene_ids": ["g4"], "state": "no_concrete_upper_lead",
              "queries": ["REF-B mutant phenotype"], "retrieval_state": "completed",
              "checked_species_scopes": ["focal_species", "other_plants"], "screened_source_ids": [],
              "reason": "No concrete upper lead in this actual inspected scope."}]))

    def invalid_error(self, out: Path, bid: str = "b001") -> str:
        """The whole-batch error, or the reasons of the records that were left out."""
        receipt = json.loads((out / "merge_receipt.json").read_text(encoding="utf-8"))
        batch = next(b for b in receipt["batches"] if b["batch_id"] == bid)
        if batch["state"] == "invalid":
            return batch["error"]
        self.assertTrue(batch["dropped_records"])
        return " | ".join(record["error"] for record in batch["dropped_records"])

    def merge(self, out: Path, **kwargs):
        return merge_initial_batches(str(self.source), kwargs.pop("bundle_sha256", self.sha),
                                     str(kwargs.pop("index", self.index)),
                                     kwargs.pop("index_sha256", self.index_sha), str(out))

    def test_happy_path_output_is_accepted_by_summarize_initial(self):
        self.happy()
        result = self.merge(self.root / "merged")
        self.assertEqual(2, result["merged_batch_count"])
        self.assertEqual(0, result["missing_batch_count"])
        payload = json.loads((self.root / "merged/initial_assessments.json").read_text(encoding="utf-8"))
        self.assertEqual(self.sha, payload["bundle_sha256"])
        self.assertEqual(["b001:A1", "b001:A2", "b002:A1"],
                         [row["assessment_id"] for row in payload["assessments"]])
        self.assertEqual(["b002:J1"], [row["assessment_id"] for row in payload["shared_assessments"]])
        self.assertEqual(["b001:S1", "b002:S1"],
                         [row["screen_id"] for row in payload["research_screening"]["reference_screens"]])
        self.assertTrue(json.loads((self.root / "merged/merge_receipt.json").read_text(encoding="utf-8"))["limitations"]
                        .startswith("Merge only concatenates original batch ledgers; no judgment is created, changed or inferred."))
        summary = summarize_initial(str(self.source), self.sha,
                                    str(self.root / "merged/initial_assessments.json"), str(self.root / "summary"))
        self.assertEqual(7, summary["gene_count"])
        self.assertEqual(4, summary["assessment_encoding"]["expanded_count"])

    def test_foreign_gene_is_rejected(self):
        self.write("b001", self.ledger("b001", [
            {"assessment_id": "b001:A1", "gene_id": "g4", "unit_id": "u4",
             "state": "not_assessed", "reason": "Outside this batch.", "assessed_fields": []}]))
        result = self.merge(self.root / "merged")
        self.assertEqual((0, 1, 0), (result["invalid_batch_count"], result["dropped_record_count"],
                                     result["assessments"]))
        self.assertIn("outside the batch", self.invalid_error(self.root / "merged"))

    def test_foreign_shared_and_screen_members_are_rejected(self):
        cases = {
            "shared": self.ledger("b001", [], [{"assessment_id": "b001:J1", "material_id": "",
                                                "member_gene_ids": ["g4"], "state": "not_assessed",
                                                "reason": "Outside this batch.", "assessed_fields": []}]),
            "screen": self.ledger("b001", [], None, [{"screen_id": "b001:S1", "member_gene_ids": ["g4"],
                                                      "state": "no_concrete_upper_lead", "reason": "x"}]),
        }
        for name, payload in cases.items():
            with self.subTest(name=name):
                self.write("b001", payload)
                self.merge(self.root / f"merged-{name}")
                self.assertIn("outside the batch", self.invalid_error(self.root / f"merged-{name}"))

    def test_record_ids_must_carry_the_batch_prefix(self):
        self.write("b001", self.ledger("b001", [
            {"assessment_id": "A1", "gene_id": "g1", "unit_id": "u1",
             "state": "not_assessed", "reason": "Prefix is missing.", "assessed_fields": []}]))
        self.merge(self.root / "merged")
        self.assertIn("must start with b001:", self.invalid_error(self.root / "merged"))

    def test_per_gene_rows_without_an_id_are_merged(self):
        # The existing ledger schema does not require IDs on per-gene assessment rows.
        self.write("b001", self.ledger("b001", [
            {"gene_id": "g1", "unit_id": "u1", "state": "assessed_no_specific_bridge",
             "reason": "Original judgment without a row ID.", "used_material_pointers": [], "assessed_fields": []}]))
        result = self.merge(self.root / "merged")
        self.assertEqual(1, result["assessments"])
        summarize_initial(str(self.source), self.sha,
                          str(self.root / "merged/initial_assessments.json"), str(self.root / "summary"))

    def test_one_row_id_may_label_a_judgment_applied_to_several_genes(self):
        # Seen in a real a101 batch: one assessment_id repeated across the genes it covers.
        self.write("b001", self.ledger("b001", [
            {"assessment_id": "b001:A1", "gene_id": "g1", "unit_id": "u1",
             "state": "not_assessed", "reason": "First record.", "assessed_fields": []},
            {"assessment_id": "b001:A1", "gene_id": "g2", "unit_id": "u2",
             "state": "not_assessed", "reason": "Second record.", "assessed_fields": []}]))
        self.assertEqual(2, self.merge(self.root / "merged")["assessments"])

    def test_duplicate_shared_and_screen_ids_are_rejected(self):
        self.write("b001", self.ledger("b001", [], [
            {"assessment_id": "b001:J1", "material_id": "", "member_gene_ids": ["g1"],
             "state": "not_assessed", "reason": "First.", "assessed_fields": []},
            {"assessment_id": "b001:J1", "material_id": "", "member_gene_ids": ["g2"],
             "state": "not_assessed", "reason": "Second.", "assessed_fields": []}]))
        self.merge(self.root / "merged")
        self.assertIn("duplicate research record ID", self.invalid_error(self.root / "merged"))

    def test_bundle_and_index_hashes_are_verified(self):
        self.happy()
        with self.assertRaisesRegex(ValueError, "bundle SHA-256 mismatch"):
            merge_initial_batches(str(self.source), "0" * 64, str(self.index), self.index_sha,
                                  str(self.root / "bad-bundle"))
        with self.assertRaisesRegex(ValueError, "batch index SHA-256 mismatch"):
            self.merge(self.root / "bad-index", index_sha256="0" * 64)
        self.assertFalse((self.root / "bad-bundle").exists())
        self.assertFalse((self.root / "bad-index").exists())

    def test_batch_ledger_must_bind_the_bundle(self):
        self.write("b001", self.ledger("b001", [
            {"assessment_id": "b001:A1", "gene_id": "g1", "unit_id": "u1",
             "state": "not_assessed", "reason": "Wrong bundle.", "assessed_fields": []}],
            bundle_sha256="1" * 64))
        self.merge(self.root / "merged")
        self.assertIn("does not bind the research bundle SHA-256", self.invalid_error(self.root / "merged"))

    def test_missing_batch_ledger_is_recorded_not_invented(self):
        self.happy()
        (self.out / "research/batches/b002/initial_assessments.json").unlink()
        result = self.merge(self.root / "merged")
        self.assertEqual(1, result["merged_batch_count"])
        self.assertEqual(1, result["missing_batch_count"])
        receipt = json.loads((self.root / "merged/merge_receipt.json").read_text(encoding="utf-8"))
        self.assertEqual(["merged", "missing"], [batch["state"] for batch in receipt["batches"]])
        payload = json.loads((self.root / "merged/initial_assessments.json").read_text(encoding="utf-8"))
        self.assertEqual(["b001:A1", "b001:A2"],
                         [row["assessment_id"] for row in payload["assessments"]])

    def test_records_downstream_would_reject_are_left_out(self):
        screen = {"screen_id": "b001:S1", "reference_gene": "REF-A", "reference_species": "Arabidopsis thaliana",
                  "member_gene_ids": ["g1"], "state": "no_concrete_upper_lead", "queries": ["REF-A mutant"],
                  "retrieval_state": "completed", "checked_species_scopes": ["focal_species"],
                  "screened_source_ids": [], "reason": "Inspected scope."}
        row = {"gene_id": "g1", "unit_id": "u1", "state": "assessed_no_specific_bridge",
               "reason": "Original judgment.", "used_material_pointers": [], "assessed_fields": []}
        cases = {
            "invented Tier label": ([{**row, "research_planning": {"provisional_tier": "tier_b_conditional"}}], [screen],
                                    "invalid provisional Tier"),
            "free-form screen": ([row], [{"screen_id": "b001:S1", "question": "q", "linked_genes": ["g1"]}],
                                 "invalid members"),
            "unknown screen state": ([row], [{**screen, "state": "done"}], "invalid state"),
            "borrowed reference": ([row], [{**screen, "reference_gene": "REF-B"}], "borrows another gene"),
            "unknown material pointer": ([{**row, "used_material_pointers": ["/units/9/genes/0"]}], [screen],
                                         "material pointer/field mismatch"),
        }
        for name, (rows, screens, error) in cases.items():
            with self.subTest(name=name):
                self.happy()
                self.write("b001", self.ledger("b001", rows, None, screens))
                out = self.root / f"merged-{name.replace(' ', '-')}"
                result = self.merge(out)
                self.assertEqual((0, 2), (result["invalid_batch_count"], result["merged_batch_count"]))
                self.assertIn(error, self.invalid_error(out))
                summarize_initial(str(self.source), self.sha, str(out / "initial_assessments.json"),
                                  str(out / "summary"))

    def test_one_bad_record_does_not_discard_the_rest_of_the_batch(self):
        # Seen in d001: one gene's pointer slip set aside a batch that held the true gene.
        good = {"gene_id": "g1", "unit_id": "u1", "state": "assessed_concrete_lead",
                "reason": "Concrete lead kept.", "used_material_pointers": [], "assessed_fields": [],
                "research_planning": {"provisional_tier": "T2", "tier_reason": "Lead.", "screen_ids": []}}
        bad = {"gene_id": "g2", "unit_id": "u2", "state": "assessed_no_specific_bridge",
               "reason": "Pointer slip.", "used_material_pointers": ["/units/9/genes/0"], "assessed_fields": []}
        depends = {"gene_id": "g2", "unit_id": "u2", "state": "assessed_no_specific_bridge",
                   "reason": "Names a screen that is dropped.", "used_material_pointers": [], "assessed_fields": [],
                   "research_planning": {"provisional_tier": "T5", "tier_reason": "x", "screen_ids": ["b001:S9"]}}
        broken_screen = {"screen_id": "b001:S9", "member_gene_ids": ["g2"], "state": "done"}
        self.write("b001", self.ledger("b001", [good, bad, depends], None, [broken_screen]))
        result = self.merge(self.root / "merged")
        self.assertEqual((0, 1, 3), (result["invalid_batch_count"], result["assessments"],
                                     result["dropped_record_count"]))
        payload = json.loads((self.root / "merged/initial_assessments.json").read_text(encoding="utf-8"))
        self.assertEqual([good], payload["assessments"])
        errors = self.invalid_error(self.root / "merged")
        for text in ("material pointer/field mismatch", "invalid state", "missing, dropped or unbound screen"):
            self.assertIn(text, errors)
        summarize_initial(str(self.source), self.sha, str(self.root / "merged/initial_assessments.json"),
                          str(self.root / "summary"))

    def test_fields_read_from_any_listed_material_count(self):
        # Seen in d001: the true gene's record listed fields read across several cards.
        # reference_function is only on the gene record, function_description only on the card.
        row = {"gene_id": "g1", "unit_id": "u1", "state": "assessed_concrete_lead", "reason": "Read gene and card.",
               "used_material_pointers": ["/units/0/genes/0", "/units/0/reference_material_cards/0"],
               "assessed_fields": ["reference_function", "function_description"]}
        self.write("b001", self.ledger("b001", [row]))
        result = self.merge(self.root / "merged")
        self.assertEqual((1, 0), (result["assessments"], result["dropped_record_count"]))
        summarize_initial(str(self.source), self.sha, str(self.root / "merged/initial_assessments.json"),
                          str(self.root / "summary"))

    def test_task_states_the_vocabularies_downstream_checks(self):
        from scripts.agentic_route_b_initial import CONTEXT_FIELDS, STATES
        from scripts.agentic_route_b_triage import SCREEN_STATES, SPECIES_SCOPES, TIERS
        for value in sorted(STATES | SCREEN_STATES | SPECIES_SCOPES | TIERS | set(CONTEXT_FIELDS)):
            self.assertIn(value, TASK)

    def test_tampered_gene_roster_is_rejected(self):
        index = json.loads(self.index.read_text(encoding="utf-8"))
        index["batches"][0]["gene_ids"] = ["g1", "g2", "g4"]
        self.index.write_text(json.dumps(index, ensure_ascii=False), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "gene roster does not match the bundle"):
            self.merge(self.root / "merged", index_sha256=digest(self.index.read_bytes()))

    def test_existing_output_directory_is_refused(self):
        self.happy()
        (self.root / "merged").mkdir()
        with self.assertRaisesRegex(ValueError, "already exists"):
            self.merge(self.root / "merged")


class BatchedResearchCliTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.bundle, self.prompts = fixture()
        self.source = self.root / "route_b_semantic_bundle.json"
        self.source.write_text(json.dumps(self.bundle, ensure_ascii=False), encoding="utf-8")
        self.sha = digest(self.source.read_bytes())
        self.prompts_path = self.root / "route_b_semantic_prompts.jsonl"
        self.prompts_path.write_text("".join(json.dumps(p, ensure_ascii=False) + "\n" for p in self.prompts),
                                     encoding="utf-8")
        self.script = ROOT / "scripts" / "agentic_route_b.py"

    def test_prepare_rejects_batching_without_a_research_bundle(self):
        out = self.root / "invalid"
        result = subprocess.run([sys.executable, str(self.script), "prepare", "--prompts", str(self.prompts_path),
                                 "--out", str(out), "--research-batch-units", "2"],
                                capture_output=True, text=True, encoding="utf-8")
        self.assertNotEqual(0, result.returncode)
        self.assertIn("--research-batch-units requires --research-bundle", result.stderr)
        self.assertFalse(out.exists())

    def test_prepare_rejects_batching_with_a_followup_plan(self):
        out = self.root / "invalid-plan"
        with self.assertRaisesRegex(ValueError, "cannot be combined with a follow-up plan"):
            route_b.prepare(str(self.prompts_path), str(out), research_bundle=str(self.source),
                            research_bundle_sha256=self.sha, followup_plan="plan.json",
                            followup_plan_sha256="0" * 64, research_batch_units=2)
        self.assertFalse(out.exists())

    def test_prepare_writes_batches_and_merge_initial_accepts_them(self):
        out = self.root / "view"
        result = subprocess.run([sys.executable, str(self.script), "prepare", "--prompts", str(self.prompts_path),
                                 "--out", str(out), "--case-id", "fixture", "--research-bundle", str(self.source),
                                 "--research-bundle-sha256", self.sha, "--research-batch-units", "3"],
                                capture_output=True, text=True, encoding="utf-8", check=True)
        manifest = json.loads(result.stdout)
        self.assertEqual(2, manifest["research_view"]["research_batches"]["batch_count"])
        self.assertTrue((out / "research/batches/b001/task.md").is_file())
        index_sha = digest((out / "research/batches/index.json").read_bytes())
        (out / "research/batches/b001/initial_assessments.json").write_text(json.dumps({
            "bundle_sha256": self.sha,
            "assessments": [{"assessment_id": "b001:A1", "gene_id": "g1", "unit_id": "u1",
                             "state": "not_assessed", "reason": "Original record.", "assessed_fields": []}],
        }, ensure_ascii=False), encoding="utf-8")
        merged = self.root / "merged"
        result = subprocess.run([sys.executable, str(self.script), "merge-initial", "--bundle", str(self.source),
                                 "--bundle-sha256", self.sha, "--batches", str(out / "research/batches/index.json"),
                                 "--batches-sha256", index_sha, "--out", str(merged)],
                                capture_output=True, text=True, encoding="utf-8", check=True)
        summary = json.loads(result.stdout)
        self.assertEqual(1, summary["merged_batch_count"])
        self.assertEqual(1, summary["missing_batch_count"])
        self.assertTrue((merged / "initial_assessments.json").is_file())


class FollowupTaskLiteratureSentenceTests(unittest.TestCase):
    def test_followup_task_collects_up_to_twenty_and_selects_five_by_tier(self):
        self.assertIn("上限約20報までできるだけ論文を集める（原著に限らずレビューも可）", FOLLOWUP_TASK)
        self.assertIn("集めた中から高Tier順に有用な最大5報を選ぶ", FOLLOWUP_TASK)
        self.assertIn("尽きたら「文献が尽きた」と記録するだけで詳しい理由は不要", FOLLOWUP_TASK)
        self.assertNotIn("約20報/有用5報は目標で", FOLLOWUP_TASK)
        self.assertNotIn("最強の適用Tier", FOLLOWUP_TASK)
        self.assertNotIn("候補件数の目標", FOLLOWUP_TASK)
        self.assertEqual("agentic-route-b-research-view-1.7.0", CONTRACT)


if __name__ == "__main__":
    unittest.main()
