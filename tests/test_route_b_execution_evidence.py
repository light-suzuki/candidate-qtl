from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from candidate_qtl.cli import main
from candidate_qtl.llm import (
    EXECUTION_EVIDENCE_SMALL_ROSTER,
    EXECUTION_EVIDENCE_UNIFORM,
    EXECUTION_EVIDENCE_VARIED,
    UNIFORM_RESPONSE_MIN_PROMPTS,
    response_execution_evidence,
)

HARNESS = Path(__file__).resolve().parents[1] / "benchmarks" / "harnesses" / "type12-harness"


def _responses(bodies: list[str]) -> list[dict[str, object]]:
    return [
        {"raw_text": body, "claims": json.loads(body).get("claims", [])}
        for body in bodies
    ]


class ResponseExecutionEvidenceTests(unittest.TestCase):
    """A uniform roster is the shape a runner produces, and must be named."""

    def test_uniform_full_roster_is_not_countable_as_per_prompt_evaluation(self) -> None:
        evidence = response_execution_evidence(_responses(['{"claims":[]}'] * 12))
        self.assertEqual(EXECUTION_EVIDENCE_UNIFORM, evidence["execution_evidence"])
        self.assertEqual(1, evidence["distinct_response_count"])
        self.assertEqual(12, evidence["empty_claim_response_count"])
        self.assertFalse(evidence["counts_as_per_prompt_evaluation"])

    def test_varied_roster_counts_as_per_prompt_evaluation(self) -> None:
        bodies = ['{"claims":[]}'] * 11 + ['{"claims":[{"claim_type":"direct_functional_bridge"}]}']
        evidence = response_execution_evidence(_responses(bodies))
        self.assertEqual(EXECUTION_EVIDENCE_VARIED, evidence["execution_evidence"])
        self.assertEqual(2, evidence["distinct_response_count"])
        self.assertEqual(11, evidence["empty_claim_response_count"])
        self.assertTrue(evidence["counts_as_per_prompt_evaluation"])

    def test_small_roster_carries_no_signal_either_way(self) -> None:
        """Three thin units genuinely can all be abstentions; do not accuse."""
        bodies = ['{"claims":[]}'] * (UNIFORM_RESPONSE_MIN_PROMPTS - 1)
        evidence = response_execution_evidence(_responses(bodies))
        self.assertEqual(EXECUTION_EVIDENCE_SMALL_ROSTER, evidence["execution_evidence"])
        self.assertFalse(evidence["counts_as_per_prompt_evaluation"])


class ExecutionEvidenceReachesTheImportManifestTests(unittest.TestCase):
    """The verdict has to travel to the first artifact readable as "Route B ran"."""

    def _case(self, root: Path, bodies: list[str]) -> tuple[Path, Path, Path]:
        prompts = root / "route_b_semantic_prompts.jsonl"
        units = []
        prompt_lines = []
        records = []
        for index, body in enumerate(bodies, 1):
            unit_id = f"u{index}"
            prompt_id = f"p{index}"
            units.append({
                "unit_id": unit_id, "identity": {}, "member_interval_genes": [],
                "stable_source_ids": [], "genes": [],
            })
            prompt_lines.append(json.dumps({
                "prompt_id": prompt_id, "route_b_unit_id": unit_id, "gene_id": "",
                "system": "system", "user": "user",
            }))
            records.append({
                "prompt_id": prompt_id, "route_b_unit_id": unit_id, "gene_id": "",
                "raw_text": body, "model_json": json.loads(body), "parse_status": "ok",
            })
        prompts.write_text("\n".join(prompt_lines) + "\n", encoding="utf-8")
        bundle = root / "route_b_semantic_bundle.json"
        bundle.write_text(json.dumps({
            "schema_version": "1.4", "route": "interval_first", "units": units,
        }) + "\n", encoding="utf-8")
        agent_file = root / "agent.json"
        agent_file.write_text(json.dumps(records) + "\n", encoding="utf-8")
        return prompts, bundle, agent_file

    def test_uniform_responses_are_recorded_through_llm_run_and_import(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            prompts, bundle, agent_file = self._case(root, ['{"claims":[]}'] * 8)
            out = root / "02_llm_run"
            self.assertEqual(0, main([
                "llm-run", "--prompts", str(prompts), "--out-dir", str(out),
                "--provider", "codex-agent", "--model", "fixture-runner",
                "--agent-file", str(agent_file), "--json",
            ]))
            manifest = json.loads((out / "llm_run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(EXECUTION_EVIDENCE_UNIFORM, manifest["execution_evidence"])
            self.assertFalse(manifest["counts_as_per_prompt_evaluation"])

            imported = root / "import"
            self.assertEqual(0, main([
                "route-b-semantic-import", "--bundle", str(bundle), "--prompts", str(prompts),
                "--responses", str(out / "llm_candidate_responses.json"),
                "--llm-manifest", str(out / "llm_run_manifest.json"),
                "--out-dir", str(imported), "--json",
            ]))
            import_manifest = json.loads(
                (imported / "semantic_claim_proposal_manifest.json").read_text(encoding="utf-8")
            )
            evidence = import_manifest["llm_execution_evidence"]
            self.assertEqual(EXECUTION_EVIDENCE_UNIFORM, evidence["execution_evidence"])
            self.assertEqual(1, evidence["distinct_response_count"])

    def test_import_without_a_manifest_records_no_verdict_rather_than_a_false_one(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            prompts, bundle, agent_file = self._case(root, ['{"claims":[]}'] * 8)
            out = root / "02_llm_run"
            self.assertEqual(0, main([
                "llm-run", "--prompts", str(prompts), "--out-dir", str(out),
                "--provider", "codex-agent", "--model", "fixture-runner",
                "--agent-file", str(agent_file), "--json",
            ]))
            imported = root / "import"
            self.assertEqual(0, main([
                "route-b-semantic-import", "--bundle", str(bundle), "--prompts", str(prompts),
                "--responses", str(out / "llm_candidate_responses.json"),
                "--out-dir", str(imported), "--json",
            ]))
            import_manifest = json.loads(
                (imported / "semantic_claim_proposal_manifest.json").read_text(encoding="utf-8")
            )
            self.assertIsNone(import_manifest["llm_execution_evidence"])


class FreezeRefusesAnUnexecutedRouteBTests(unittest.TestCase):
    """The freeze is the last honest place to record a lane that never ran."""

    def _generation(self, root: Path, evidence: str, executed: bool) -> None:
        cell = "type1_n200"
        (root / "executor" / cell).mkdir(parents=True)
        (root / "executor" / cell / "cases.tsv").write_text(
            "type\tcase_id\tphenotype\n" + "".join(f"type1\tq{n:03d}\tp\n" for n in (1, 2)),
            encoding="utf-8",
        )
        for index in (1, 2):
            case_id = f"q{index:03d}"
            run = root / "runs" / cell / case_id
            (run / "b04_canonical_semantic" / "07_supported_union_for_family").mkdir(parents=True)
            (run / "b04_canonical_semantic" / "02_llm_run").mkdir(parents=True)
            (run / "z04_canonical_wet_handoff").mkdir(parents=True)
            (run / "b04_canonical_semantic" / "canonical_route_b_extension_manifest.json").write_text(
                json.dumps({"state": "wet_reviewable"}), encoding="utf-8"
            )
            (run / "z04_canonical_wet_handoff" / "wet_candidates.tsv").write_text(
                "candidate_gene_id\nG1\n", encoding="utf-8"
            )
            (run / "b04_canonical_semantic" / "07_supported_union_for_family" / "supported_member_union.tsv").write_text(
                "candidate_gene_id\nG1\n", encoding="utf-8"
            )
            (run / "b04_canonical_semantic" / "02_llm_run" / "llm_run_manifest.json").write_text(
                json.dumps({
                    "execution_evidence": evidence,
                    "counts_as_per_prompt_evaluation": executed,
                    "distinct_response_count": 1 if not executed else 4,
                }),
                encoding="utf-8",
            )

    def _freeze(self, root: Path, *args: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(HARNESS / "freeze_predictions.py"), "type1_n200", *args],
            capture_output=True, text=True, env={"CQTL_BENCH_ROOT": str(root), "PATH": "/usr/bin:/bin"},
            check=False,
        )

    @unittest.skip("public release: uses the benchmark harness freeze_predictions.py (not shipped)")
    def test_uniform_route_b_is_recorded_and_can_be_made_fatal(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._generation(root, EXECUTION_EVIDENCE_UNIFORM, executed=False)

            recorded = self._freeze(root)
            self.assertEqual(0, recorded.returncode, recorded.stderr)
            manifest = json.loads(
                (root / "freeze" / "prediction_freeze_manifest.type1_n200.json").read_text(encoding="utf-8")
            )
            self.assertEqual(["q001", "q002"], manifest["route_b_not_countable_as_executed"])
            self.assertEqual(
                [EXECUTION_EVIDENCE_UNIFORM] * 2,
                [row["route_b_execution_evidence"] for row in manifest["cases"]],
            )

            refused = self._freeze(root, "--require-route-b-executed")
            self.assertNotEqual(0, refused.returncode)
            self.assertIn("not countable as executed", refused.stderr)

    @unittest.skip("public release: uses the benchmark harness freeze_predictions.py (not shipped)")
    def test_varied_route_b_freezes_under_the_strict_flag(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._generation(root, EXECUTION_EVIDENCE_VARIED, executed=True)
            result = self._freeze(root, "--require-route-b-executed")
            self.assertEqual(0, result.returncode, result.stderr)
            manifest = json.loads(
                (root / "freeze" / "prediction_freeze_manifest.type1_n200.json").read_text(encoding="utf-8")
            )
            self.assertEqual([], manifest["route_b_not_countable_as_executed"])


if __name__ == "__main__":
    unittest.main()
