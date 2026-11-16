"""Rank invariance around a blast-panel freeze (issue #374 / #355).

Three guarantees are pinned here:

1. ``candidate_qtl.blast_panel`` is a deterministic island: importing it never
   pulls in ``hardening`` / ``shortlist`` / ``family_index`` / ``llm_pathway`` /
   ``ranking_policy`` (sys.modules snapshot, plus an AST check on the source),
   so a freeze cannot mutate ranking state through an import side effect.
2. ``evidence-rank`` on the committed pea_tl fixture produces byte-identical
   ``ranked_candidates.tsv`` before and after ``freeze_blast_panel`` runs in
   the same process.
3. The same holds when ``blast-panel-search`` (mocked blastp) runs in between:
   a search emits ``blast_evidence.tsv`` whose rows all carry
   ``rank_contribution=none`` and ranking output stays byte-identical.
"""

from __future__ import annotations

import ast
import hashlib
import importlib
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.blast_panel import freeze_blast_panel
from candidate_qtl.blast_evidence import run_panel_blastp
from candidate_qtl.cli import main
from tests.test_blast_evidence import MockBlastTools

FIXTURE = Path(__file__).resolve().parent / "fixtures" / "route_a"
PEA_TL = Path(__file__).resolve().parent / "fixtures" / "pea_tl"

BANNED_MODULES = (
    "candidate_qtl.hardening",
    "candidate_qtl.shortlist",
    "candidate_qtl.family_index",
    "candidate_qtl.llm_pathway",
    "candidate_qtl.ranking_policy",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class BlastPanelImportArchitectureTests(unittest.TestCase):
    def test_importing_blast_panel_loads_no_ranking_or_llm_modules(self) -> None:
        before = set(sys.modules)
        importlib.import_module("candidate_qtl.blast_panel")
        loaded_new = set(sys.modules) - before
        for banned in BANNED_MODULES:
            self.assertNotIn(banned, loaded_new, f"blast_panel import pulled in {banned}")

    def test_blast_panel_source_never_imports_ranking_or_llm_modules(self) -> None:
        source_path = Path(importlib.import_module("candidate_qtl.blast_panel").__file__)
        tree = ast.parse(source_path.read_text(encoding="utf-8"), filename=str(source_path))
        imported: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                base = node.module or ""
                if node.level:
                    base = ("candidate_qtl." + base) if base else "candidate_qtl"
                imported.extend(f"{base}.{alias.name}" for alias in node.names)
        for banned in BANNED_MODULES:
            self.assertFalse(
                any(name == banned or name.startswith(banned + ".") for name in imported),
                f"blast_panel source imports {banned}",
            )


class BlastPanelRankInvarianceTests(unittest.TestCase):
    def _evidence_rank(self, root: Path, candidates: Path, evidence: Path, out_dir: Path) -> str:
        self.assertEqual(
            0,
            main([
                "evidence-rank",
                "--candidates", str(candidates),
                "--evidence", str(evidence),
                "--out-dir", str(out_dir),
            ]),
        )
        return _sha256(out_dir / "ranked_candidates.tsv")

    def test_evidence_rank_output_unchanged_around_a_panel_freeze(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            candidates = root / "candidates.tsv"
            candidates.write_bytes((PEA_TL / "candidates.tsv").read_bytes())
            evidence = root / "evidence.tsv"
            evidence.write_bytes((PEA_TL / "evidence.tsv").read_bytes())

            before = self._evidence_rank(root, candidates, evidence, root / "rank_before")

            frozen = freeze_blast_panel(
                str(FIXTURE / "panel_ok.tsv"),
                str(root / "freeze"),
                case_id="rank-invariance",
            )
            self.assertEqual("verified", frozen["panel_status"])

            after = self._evidence_rank(root, candidates, evidence, root / "rank_after")

            self.assertEqual(before, after, "evidence-rank output must not depend on a panel freeze")

    def test_evidence_rank_output_unchanged_around_a_panel_search(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            candidates = root / "candidates.tsv"
            candidates.write_bytes((PEA_TL / "candidates.tsv").read_bytes())
            evidence = root / "evidence.tsv"
            evidence.write_bytes((PEA_TL / "evidence.tsv").read_bytes())

            before = self._evidence_rank(root, candidates, evidence, root / "rank_before")

            frozen = freeze_blast_panel(
                str(FIXTURE / "panel_ok.tsv"),
                str(root / "freeze"),
                case_id="rank-invariance-search",
            )
            self.assertEqual("verified", frozen["panel_status"])

            with MockBlastTools():
                search = run_panel_blastp(
                    str(frozen["manifest"]),
                    str(FIXTURE / "interval_proteins.fa"),
                    str(root / "search"),
                )
            self.assertEqual("evaluable", search["provider_status"])

            blast_evidence = pd.read_csv(search["evidence"], sep="\t", dtype=str).fillna("")
            self.assertEqual(3, len(blast_evidence))
            self.assertTrue((blast_evidence["rank_contribution"] == "none").all())

            after = self._evidence_rank(root, candidates, evidence, root / "rank_after")

            self.assertEqual(before, after, "evidence-rank output must not depend on a panel blast search")


if __name__ == "__main__":
    unittest.main()
