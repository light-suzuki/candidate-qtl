"""Preflight uses the original A transport and never creates import artifacts."""
from __future__ import annotations

import hashlib
import io
import json
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from copy import deepcopy
from pathlib import Path

from candidate_qtl.cli import run


def _hypothesis() -> dict:
    return {
        "case_id": "transport-case",
        "pathway_or_process": "cellular process",
        "gene_family": "FAMILY",
        "mechanism_statement": "The family could alter the cellular process.",
        "source_ids": "PMID:123",
        "paper_reports_json": [{
            "source_id": "PMID:123",
            "source_title": "Synthetic primary paper",
            "reported_gene": "EXAMPLE_GENE",
            "reported_gene_id": "NCBI Gene:111111",
            "reported_species": "Plant beta",
            "reference_gene_id": "REFERENCE_1",
            "gene_correspondence": "homolog",
            "initial_tier_reason": "The reported observation may inform the process.",
            "paper_checklist": {"phenotype_reported": "yes", "identity_resolved": "yes"},
        }],
    }


class MechanismImportPreflightTests(unittest.TestCase):
    def _run(self, args: list[str]) -> tuple[int, dict]:
        stdout = io.StringIO()
        with redirect_stdout(stdout), redirect_stderr(io.StringIO()):
            code = run(["mechanism-import", *args, "--preflight", "--json"])
        return code, json.loads(stdout.getvalue())

    def test_valid_payload_is_read_only_and_reports_the_raw_digest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = root / "payload.json"
            raw = json.dumps({"hypotheses": [_hypothesis()]}, indent=2).encode()
            payload.write_bytes(raw)
            out = root / "must-not-exist"
            code, result = self._run(["--payload", str(payload), "--out-dir", str(out)])
            self.assertEqual(0, code)
            self.assertEqual("valid", result["status"])
            self.assertEqual(1, result["rows_read"])
            self.assertEqual(hashlib.sha256(raw).hexdigest(), result["raw_payload_sha256"])
            self.assertEqual(raw, payload.read_bytes())
            self.assertFalse(out.exists())
            self.assertEqual([payload], list(root.iterdir()))

    def test_name_only_reference_is_advised_without_invalidating_payload(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            named = {**_hypothesis(), "reference_gene_name": "NODE1", "reference_species": "Plant alpha",
                     "proposed_stable_id": ""}
            with_id = {**_hypothesis(), "reference_gene_name": "NODE2", "proposed_stable_id": "UniProt:Q00002"}
            payload = root / "payload.json"
            payload.write_text(json.dumps({"hypotheses": [named, with_id, _hypothesis()]}), encoding="utf-8")
            code, result = self._run(["--payload", str(payload)])
            self.assertEqual(0, code)
            self.assertEqual("valid", result["status"])
            self.assertEqual([("rows[0]", "reference_gene_without_stable_id", "NODE1")],
                             [(a["path"], a["advisory"], a["reference_gene_name"]) for a in result["advisories"]])

    def test_invalid_second_record_reports_its_path_without_partial_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = root / "payload.json"
            bad = _hypothesis()
            bad["paper_reports_json"][0]["reported_gene"] = ["unsupported scientific array"]
            raw = json.dumps({"hypotheses": [_hypothesis(), bad]}).encode()
            payload.write_bytes(raw)
            code, result = self._run(["--payload", str(payload)])
            self.assertEqual(1, code)
            self.assertEqual("invalid", result["status"])
            self.assertEqual("$.hypotheses[1].paper_reports_json", result["errors"][0]["path"])
            self.assertIn("must be strings", result["errors"][0]["error"])
            self.assertEqual(raw, payload.read_bytes())
            self.assertEqual([payload], list(root.iterdir()))

    def test_scientific_and_source_binding_errors_stay_rejected(self) -> None:
        for changed in (
            {"source_id": "not-a-citation"},
            {"source_id": "PMID:456"},
            {"paper_checklist": {"phenotype_reported": "invented_answer"}},
            {"invented_scientific_fact": "unsupported"},
        ):
            with self.subTest(changed=changed), tempfile.TemporaryDirectory() as tmp:
                payload = Path(tmp) / "payload.json"
                record = _hypothesis()
                record["paper_reports_json"][0].update(changed)
                payload.write_text(json.dumps({"hypotheses": [record]}), encoding="utf-8")
                code, result = self._run(["--payload", str(payload)])
                self.assertEqual(1, code)
                self.assertEqual("$.hypotheses[0].paper_reports_json", result["errors"][0]["path"])

    def test_tsv_cells_reach_the_actual_paper_normalizer(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            payload = Path(tmp) / "payload.tsv"
            record = _hypothesis()
            record["paper_reports_json"][0]["reported_gene"] = ["invalid"]
            payload.write_text("source_ids\tpaper_reports_json\nPMID:123\t"
                               + json.dumps(record["paper_reports_json"]) + "\n", encoding="utf-8")
            code, result = self._run(["--payload", str(payload)])
            self.assertEqual(1, code)
            self.assertEqual("rows[0].paper_reports_json", result["errors"][0]["path"])

    def test_nested_annotations_survive_preflight_and_real_import(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = root / "payload.json"
            annotations = {"counterevidence": ["Contrary observation remains unresolved."],
                           "full_paper_read": False,
                           "paper_checklist_reasons": {"identity_resolved": "Reported accession."}}
            record = _hypothesis()
            record["paper_reports_json"][0]["research_annotations"] = deepcopy(annotations)
            raw = json.dumps({"hypotheses": [record]}).encode()
            payload.write_bytes(raw)
            code, result = self._run(["--payload", str(payload)])
            self.assertEqual(0, code)
            self.assertEqual("valid", result["status"])
            with redirect_stdout(io.StringIO()):
                self.assertEqual(0, run(["mechanism-import", "--payload", str(payload),
                                         "--out-dir", str(root / "import")]))
            import csv
            with (root / "import" / "mechanism_hypotheses.tsv").open(encoding="utf-8", newline="") as stream:
                imported = next(csv.DictReader(stream, delimiter="\t"))
            reports = json.loads(imported["paper_reports_json"])
            self.assertEqual(annotations, reports[0]["research_annotations"])
            self.assertEqual(raw, payload.read_bytes())

    def test_freeze_hash_is_checked_without_writing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = root / "payload.json"
            payload.write_text(json.dumps({"hypotheses": [_hypothesis()]}), encoding="utf-8")
            freeze = root / "payload_freeze.tsv"
            freeze.write_text("case_id\tmechanism_payload_sha256\ntransport-case\twrong-hash\n", encoding="utf-8")
            code, result = self._run(["--payload", str(payload), "--mechanism-freeze-manifest", str(freeze),
                                      "--case-id", "transport-case"])
            self.assertEqual(1, code)
            self.assertIn("freeze mismatch", result["errors"][0]["error"])
            self.assertEqual(2, len(list(root.iterdir())))

    def test_normal_import_still_requires_an_output_directory(self) -> None:
        with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            self.assertEqual(1, run(["mechanism-import", "--payload", "unused.json"]))


if __name__ == "__main__":
    unittest.main()
