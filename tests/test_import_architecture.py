from __future__ import annotations

import ast
import importlib
import pkgutil
import unittest
from pathlib import Path


class ImportArchitectureTests(unittest.TestCase):
    def test_package_import_is_side_effect_free(self) -> None:
        # Importing the package must never patch functions or shared maps in other
        # modules. Pin the surviving production surfaces after the ranking merge
        #: the single ranking implementation lives in ``hardening`` and the
        # shared building blocks live in ``evidence_io``.
        evidence = importlib.import_module("candidate_qtl.evidence")
        evidence_io = importlib.import_module("candidate_qtl.evidence_io")
        hardening = importlib.import_module("candidate_qtl.hardening")
        core = importlib.import_module("candidate_qtl.core")
        original = {
            "evidence_build": evidence.evidence_build,
            "orthology_plan": evidence.orthology_plan,
            "prioritize_evidence": hardening.prioritize_evidence,
            "candidate_pipeline": hardening.candidate_pipeline,
            "candidate_frame": evidence_io.candidate_frame,
            "evidence_row": evidence_io.evidence_row,
            "extract_candidates": core.extract_candidates,
            "confidence": dict(evidence_io.CONFIDENCE_MULTIPLIERS),
        }

        package = importlib.import_module("candidate_qtl")
        importlib.reload(package)

        self.assertIs(original["evidence_build"], evidence.evidence_build)
        self.assertIs(original["orthology_plan"], evidence.orthology_plan)
        self.assertIs(original["prioritize_evidence"], hardening.prioritize_evidence)
        self.assertIs(original["candidate_pipeline"], hardening.candidate_pipeline)
        self.assertIs(original["candidate_frame"], evidence_io.candidate_frame)
        self.assertIs(original["evidence_row"], evidence_io.evidence_row)
        self.assertIs(original["extract_candidates"], core.extract_candidates)
        self.assertEqual(original["confidence"], evidence_io.CONFIDENCE_MULTIPLIERS)

    def test_single_production_ranking_implementation(self) -> None:
        # The old duplicate ranking routes on ``evidence`` are gone; the one
        # production implementation lives on ``hardening``.
        evidence = importlib.import_module("candidate_qtl.evidence")
        hardening = importlib.import_module("candidate_qtl.hardening")
        for removed in ("prioritize_evidence", "candidate_pipeline", "DEFAULT_WEIGHTS", "EVIDENCE_FAMILIES", "_load_weights"):
            self.assertFalse(hasattr(evidence, removed), f"evidence.{removed} should have been removed")
        self.assertTrue(callable(hardening.prioritize_evidence))
        self.assertTrue(callable(hardening.candidate_pipeline))

    def test_no_module_imports_private_evidence_helpers(self) -> None:
        # The shared I/O / provenance helpers moved to ``evidence_io``; no other
        # module may reach into ``evidence`` for a private ``_``-prefixed name
        #. Public functions (evidence_build, orthology_plan, ...) stay
        # importable from ``evidence``.
        package = importlib.import_module("candidate_qtl")
        package_dir = Path(package.__file__).parent
        offenders: list[str] = []
        for module_info in pkgutil.walk_packages([str(package_dir)], prefix="candidate_qtl."):
            source_path = package_dir / (module_info.name.split(".", 1)[1].replace(".", "/") + ".py")
            if not source_path.exists() or source_path.name == "evidence.py":
                continue
            tree = ast.parse(source_path.read_text(encoding="utf-8"), filename=str(source_path))
            for node in ast.walk(tree):
                if not isinstance(node, ast.ImportFrom):
                    continue
                is_evidence = (node.level == 1 and node.module == "evidence") or node.module == "candidate_qtl.evidence"
                if not is_evidence:
                    continue
                for alias in node.names:
                    if alias.name.startswith("_"):
                        offenders.append(f"{source_path.name}: from .evidence import {alias.name}")
        self.assertEqual(offenders, [], f"private evidence helpers imported across module boundary: {offenders}")

    def test_shared_helpers_resolve_to_evidence_io(self) -> None:
        # Modules that reuse the shared helpers get the same objects that
        # ``evidence_io`` defines, proving the boundary is real and not a copy.
        evidence_io = importlib.import_module("candidate_qtl.evidence_io")
        hardening = importlib.import_module("candidate_qtl.hardening")
        evidence = importlib.import_module("candidate_qtl.evidence")
        ranking_policy = importlib.import_module("candidate_qtl.ranking_policy")
        self.assertIs(hardening._candidate_frame, evidence_io.candidate_frame)
        self.assertIs(hardening._row, evidence_io.evidence_row)
        self.assertIs(evidence._row, evidence_io.evidence_row)
        self.assertIs(ranking_policy._row, evidence_io.evidence_row)

    def test_command_modules_import_hardened_routes_explicitly(self) -> None:
        extraction_commands = importlib.import_module("candidate_qtl.commands.extraction")
        ranking_commands = importlib.import_module("candidate_qtl.commands.ranking")
        evidence_commands = importlib.import_module("candidate_qtl.commands.evidence")
        hardening = importlib.import_module("candidate_qtl.hardening")
        evidence_policy = importlib.import_module("candidate_qtl.evidence_policy")

        self.assertIs(extraction_commands.annotation_keyword_evidence, hardening.annotation_keyword_evidence)
        self.assertIs(ranking_commands.prioritize_evidence, hardening.prioritize_evidence)
        self.assertIs(ranking_commands.candidate_pipeline, hardening.candidate_pipeline)
        self.assertIs(evidence_commands.evidence_build, evidence_policy.evidence_build)


if __name__ == "__main__":
    unittest.main()
