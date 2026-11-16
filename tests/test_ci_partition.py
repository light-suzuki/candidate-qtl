from __future__ import annotations

import importlib.util
import importlib
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_test_partition.py"


def _module():
    spec = importlib.util.spec_from_file_location("run_test_partition", SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class CiPartitionTests(unittest.TestCase):
    @unittest.skip("public release: benchmark partition tests are not shipped")
    def test_partitions_are_disjoint_and_complete(self) -> None:
        module = _module()
        all_tests = set(module.modules_for("all"))
        core = set(module.modules_for("core"))
        benchmark = set(module.modules_for("benchmark"))
        self.assertFalse(core & benchmark)
        self.assertEqual(all_tests, core | benchmark)
        self.assertIn("tests.test_discovery_chain_benchmark", benchmark)
        self.assertIn("tests.test_supported_member_union", core)

    def test_fixture_and_truth_harnesses_never_enter_core(self) -> None:
        module = _module()
        core = set(module.modules_for("core"))
        for name in module.BENCHMARK_EXACT:
            self.assertNotIn(f"tests.{name}", core)

    def test_partition_runner_preserves_discovery_harness_imports(self) -> None:
        _module()
        imported = importlib.import_module("tests.test_issue_611_curated_function_gate")
        self.assertTrue(hasattr(imported, "CuratedFunctionGateTests"))

    @unittest.skip("public release: benchmark partition/job is not shipped")
    def test_local_ci_uses_the_same_partitions_as_hosted_ci(self) -> None:
        script = (ROOT / "scripts" / "local_ci.sh").read_text(encoding="utf-8")
        self.assertIn("--partition core", script)
        self.assertIn("--partition benchmark", script)
        self.assertNotIn("coverage run -m unittest discover", script)


if __name__ == "__main__":
    unittest.main()
