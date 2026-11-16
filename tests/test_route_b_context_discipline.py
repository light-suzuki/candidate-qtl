"""Context handoff and output discipline for the ordinary Route B flow (#1004).

These tests pin already delivered launch bytes. A failure means an existing
run's worker instruction artifact would change silently; make that a new
generation instead of updating these constants.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from scripts import agentic_route_b as route_b
from scripts.agentic_route_b_material import _worker_task
from scripts.agentic_route_b_research import CONTRACT, FOLLOWUP_TASK, TASK
from tests.test_agentic_route_b_material import semantic_prompt

ROOT = Path(__file__).resolve().parents[1]

# origin/main fe943d3 outputs for fixture_rows() below, computed with
# `git show origin/main:scripts/...` copies outside the worktree.
WORKER_TASK_SHA256 = "e47ff25deae9b9fe1da16fddfa67ea97af11dc77b576196c82a5a3954e839ebe"
WORKER_TASK_WITH_PARENT_SHA256 = "926dcb5fc5a7b5668812df60208b1427bca1501d83a04c1e930e3aff8670e29b"
PREPARE_TREE_SHA256 = "2073cc58b5735b2ad6b5d73779c223f47bf52210870344c57ce81345b4dd0953"


def fixture_rows() -> list[dict]:
    return [semantic_prompt(1), semantic_prompt(2, "Second complete function " * 40, "PMID:2")]


def _newline_stable(data: bytes) -> bytes:
    """Text-mode writers spell "\\n" as "\\r\\n" on Windows; compare content."""
    return data.replace(b"\r\n", b"\n")


def _normalize_manifest(data: bytes, root: Path) -> bytes:
    """Drop the temporary directory prefix and the platform's newline spelling."""
    root = root.resolve()
    bases = (root, root.parent)

    def walk(value):
        if isinstance(value, str) and os.path.isabs(value):
            path = Path(value).resolve()
            for base in bases:
                try:
                    return "<ROOT>/" + path.relative_to(base).as_posix()
                except ValueError:
                    continue
            return "<ABS>/" + path.name
        if isinstance(value, str):
            # Windows writes relative shard paths with backslashes.
            return value.replace("\\", "/")
        if isinstance(value, list):
            return [walk(item) for item in value]
        if isinstance(value, dict):
            return {key: walk(item) for key, item in value.items()}
        return value

    manifest = json.loads(data)
    # These two fields record the newline spelling of the files they hash.
    for shard in manifest.get("shards", []):
        shard["sha256"] = hashlib.sha256(
            _newline_stable((root / str(shard["file"]).replace("\\", "/")).read_bytes())).hexdigest()
    source = Path(str(manifest.get("source_prompts", "")))
    if source.is_file():
        manifest["source_sha256"] = hashlib.sha256(_newline_stable(source.read_bytes())).hexdigest()
    return json.dumps(walk(manifest), ensure_ascii=False, sort_keys=True).encode("utf-8")


def tree_digest(root: Path) -> str:
    """Digest a prepare output directory independently of location and newlines."""
    digest = hashlib.sha256()
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        data = _newline_stable(path.read_bytes())
        if path.name == "prepare_manifest.json":
            data = _normalize_manifest(path.read_bytes(), root)
        digest.update(path.relative_to(root).as_posix().encode("utf-8") + b"\0")
        digest.update(hashlib.sha256(data).digest())
    return digest.hexdigest()


def write_rows(path: Path, rows: list[dict]) -> None:
    text = "".join(json.dumps(row, ensure_ascii=False) + "\n" for row in rows)
    path.write_text(text, encoding="utf-8", newline="\n")


class ContextDisciplineTaskTests(unittest.TestCase):
    def test_research_task_states_the_output_and_handoff_discipline(self):
        for text in ("4,000字", "8,000字", "truncated", "SHA-256", "artifact", "新しい文脈",
                     "lane", "provenance", "予算消化", "連続した行範囲"):
            self.assertIn(text, TASK)
        self.assertIn("実行環境（ローカル/クラウド、provider）に依存しない", TASK)
        # The replaced abstract line must not survive as the only output rule.
        self.assertNotIn("長い既読原文を繰り返し会話へ出力しない。詳細はartifactへ保存し", TASK)

    def test_followup_task_resumes_from_artifacts_not_conversation_history(self):
        self.assertNotIn("同じ原担当の追加調査・正式回答への継続", FOLLOWUP_TASK)
        for text in ("research/manifest.json", "path/SHA-256", "新しい文脈", "原担当lane",
                     "4,000字", "truncated", "実行環境", "研究台帳を指定した別の新規文脈"):
            self.assertIn(text, FOLLOWUP_TASK)
        self.assertIn("既読bundle全体", FOLLOWUP_TASK)

    def test_contract_version_reflects_the_new_instruction_generation(self):
        self.assertEqual("agentic-route-b-research-view-1.7.0", CONTRACT)


class FrozenLaunchBytesTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.prompts = self.root / "prompts.jsonl"
        write_rows(self.prompts, fixture_rows())

    def test_worker_task_without_ledger_matches_origin_main_bytes(self):
        self.assertEqual(
            WORKER_TASK_SHA256,
            hashlib.sha256(_worker_task("reading/b0000/assignment.json", "case-x")).hexdigest(),
        )
        self.assertEqual(
            WORKER_TASK_WITH_PARENT_SHA256,
            hashlib.sha256(_worker_task("reading/b0000/assignment.json", "case-x", "a" * 64)).hexdigest(),
        )

    def test_prepare_without_ledger_matches_origin_main_output_tree(self):
        out = self.root / "prepared"
        route_b.prepare(str(self.prompts), str(out))
        self.assertEqual(PREPARE_TREE_SHA256, tree_digest(out))


class ResearchLedgerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.prompts = self.root / "prompts.jsonl"
        write_rows(self.prompts, fixture_rows())
        self.ledger = self.root / "initial_assessments.json"
        self.ledger.write_text(json.dumps({"bundle_sha256": "0" * 64, "assessments": [
            {"unit_id": "unit-1", "gene_id": "gene-1", "state": "assessed_claim",
             "reason": "original recorded judgment", "used_material_pointers": ["/units/0/genes/0"]},
        ]}, ensure_ascii=False), encoding="utf-8")
        self.ledger_sha = hashlib.sha256(self.ledger.read_bytes()).hexdigest()

    def launch_texts(self, out: Path) -> list[str]:
        manifest = json.loads((out / "prepare_manifest.json").read_text(encoding="utf-8"))
        return [(out / shard["worker_task"]["file"]).read_text(encoding="utf-8")
                for shard in manifest["reading_view"]["shards"]]

    def test_ledger_rules_are_absent_unless_requested(self):
        out = self.root / "plain"
        manifest = route_b.prepare(str(self.prompts), str(out))
        self.assertNotIn("research_ledgers", manifest)
        for text in self.launch_texts(out):
            self.assertNotIn("台帳に判断なし", text)
            self.assertNotIn("Research ledger transcription mode", text)

    def test_ledger_is_recorded_and_adds_transcription_rules(self):
        out = self.root / "ledger"
        manifest = route_b.prepare(str(self.prompts), str(out), research_ledger=[str(self.ledger)],
                                   research_ledger_sha256=[self.ledger_sha])
        name = f"research-ledger/00-{self.ledger_sha[:16]}.json"
        self.assertEqual([{"file": name, "sha256": self.ledger_sha, "source_path": str(self.ledger.resolve())}],
                         manifest["research_ledgers"])
        self.assertEqual(self.ledger.read_bytes(), (out / name).read_bytes())
        for text in self.launch_texts(out):
            self.assertIn(name, text)
            self.assertNotIn(str(self.root), text)
            self.assertIn(self.ledger_sha, text)
            self.assertIn("台帳に判断なし", text)
            self.assertIn("no new research", text)
            self.assertIn("read_receipt", text)
            self.assertIn("4,000 characters", text)
            self.assertIn("truncated", text)
            self.assertIn("do not depend on the execution environment", text)
            self.assertLessEqual(max(map(len, text.splitlines())), 1000)

    def test_mismatched_or_unpaired_ledger_is_rejected_before_output(self):
        cases = {
            "sha": ({"research_ledger": [str(self.ledger)], "research_ledger_sha256": ["0" * 64]}, "SHA-256"),
            "missing_pair": ({"research_ledger": [str(self.ledger)]}, "together"),
            "missing_file": ({"research_ledger": [str(self.root / "absent.json")],
                              "research_ledger_sha256": ["0" * 64]}, "readable file"),
        }
        for name, (kwargs, message) in cases.items():
            with self.subTest(name=name):
                out = self.root / f"rejected-{name}"
                with self.assertRaisesRegex(ValueError, message):
                    route_b.prepare(str(self.prompts), str(out), **kwargs)
                self.assertFalse(out.exists())

    def test_cli_rejects_mismatch_and_records_a_match(self):
        script = ROOT / "scripts" / "agentic_route_b.py"
        rejected = self.root / "cli-rejected"
        result = subprocess.run([sys.executable, str(script), "prepare", "--prompts", str(self.prompts),
                                 "--out", str(rejected), "--research-ledger", str(self.ledger),
                                 "--research-ledger-sha256", "0" * 64],
                                capture_output=True, text=True, encoding="utf-8")
        self.assertNotEqual(0, result.returncode)
        self.assertIn("research ledger SHA-256 mismatch", result.stderr)
        self.assertFalse(rejected.exists())
        accepted = self.root / "cli-accepted"
        result = subprocess.run([sys.executable, str(script), "prepare", "--prompts", str(self.prompts),
                                 "--out", str(accepted), "--research-ledger", str(self.ledger),
                                 "--research-ledger-sha256", self.ledger_sha],
                                capture_output=True, text=True, encoding="utf-8", check=True)
        self.assertEqual(self.ledger_sha, json.loads(result.stdout)["research_ledgers"][0]["sha256"])


if __name__ == "__main__":
    unittest.main()
