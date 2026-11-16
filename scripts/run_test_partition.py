#!/usr/bin/env python3
"""Run deterministic core/benchmark unittest partitions for issue #660."""

from __future__ import annotations

import argparse
import sys
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
TESTS = REPO / "tests"
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
if str(TESTS) not in sys.path:
    # ``unittest discover -s tests`` historically exposes this directory for
    # a few reusable test harness modules; preserve that import contract.
    sys.path.insert(0, str(TESTS))

# Filename classification is intentionally conservative.  Benchmark truth,
# frozen generations, scorers, audit harnesses and fixture-governance tests run
# outside ordinary core CI.  Production regressions remain in the core lane
# even when they were originally discovered by a benchmark.
BENCHMARK_EXACT = frozenset({
    "test_audit_freeze",
    "test_blind_e2e_contract",
    "test_command_log",
    "test_command_log_validation",
    "test_corpus_masking",
    "test_lane_coverage_declaration",
    "test_query_eligibility",
    "test_opencode_campaign",
    "test_opencode_manager_science_terminal",
    "test_opencode_science_adapter",
    "test_opencode_science_deadlines",
    "test_opencode_science_pack",
    "test_opencode_science_packets",
    "test_opencode_science_pipeline",
    "test_opencode_science_result",
    "test_opencode_type1_pilot",
    "test_type12_harness_freeze",
    "test_type2_legume10_fixture",
    "test_type2_lmp_fixture_and_scorer",
    "test_type2_lmp_route_a_only",
    "test_type2_observation_distance",
    "test_restored_type2_fixtures",
})
BENCHMARK_MARKERS = ("benchmark",)


def test_modules() -> list[str]:
    return [f"tests.{path.stem}" for path in sorted(TESTS.glob("test_*.py"))]


def is_benchmark_module(module: str) -> bool:
    stem = module.rsplit(".", 1)[-1]
    return stem in BENCHMARK_EXACT or any(marker in stem for marker in BENCHMARK_MARKERS)


def modules_for(partition: str) -> list[str]:
    modules = test_modules()
    if partition == "all":
        return modules
    return [module for module in modules if is_benchmark_module(module) == (partition == "benchmark")]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--partition", choices=("core", "benchmark", "all"), required=True)
    parser.add_argument("--list", action="store_true", help="Print selected modules without running them.")
    parser.add_argument("--verbosity", type=int, default=1)
    parser.add_argument("--durations", type=int, default=10, help="Report the N slowest tests (default: 10).")
    parser.add_argument(
        "--no-buffer",
        action="store_true",
        help="Stream stdout/stderr from successful tests (buffered by default to keep CI logs readable).",
    )
    args = parser.parse_args(argv)
    modules = modules_for(args.partition)
    if args.list:
        print("\n".join(modules))
        return 0
    suite = unittest.defaultTestLoader.loadTestsFromNames(modules)
    result = unittest.TextTestRunner(
        verbosity=args.verbosity,
        buffer=not args.no_buffer,
        durations=args.durations,
    ).run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
