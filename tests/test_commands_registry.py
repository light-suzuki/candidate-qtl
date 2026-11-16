"""The feature-command registry is the complete CLI ownership contract."""

from __future__ import annotations

import argparse
from types import SimpleNamespace
import unittest

from candidate_qtl.cli import build_parser
from candidate_qtl.commands import COMMAND_NAMES, _collect_handlers, command_names, handler_for, register_all


EXPECTED_COMMANDS = {
    "assembly-bridge",
    "adaptive-refinement-apply",
    "adaptive-synthesis-import",
    "adaptive-synthesis-request",
    "allele-compare",
    "annotate-plan",
    "annotation-keywords",
    "annotation-rescue",
    "assay-plan",
    "auto-triage",
    "bait-panel-freeze",
    "benchmark-reference",
    "blast-panel-freeze",
    "blast-panel-search",
    "blast-panel-verify",
    "budget-llm-review",
    "budget-run",
    "build-executor-package",
    "candidate-budget-gate",
    "candidate-elimination",
    "candidate-flags",
    "candidate-pipeline",
    "case-run",
    "claim-verification-benchmark",
    "collect-sanger-markers",
    "comparative-locus-benchmark",
    "comparative-locus-reconstruct",
    "corpus-mask",
    "cosegregation",
    "dark-candidates",
    "direction-evidence",
    "discover",
    "discover-promote",
    "doctor",
    "dry-prioritize",
    "dual-route-compare",
    "evidence-build",
    "evidence-rank",
    "evidence-term-promote",
    "evidence-tier-benchmark",
    "evidence-tiers",
    "exon-coverage-check",
    "exon-sequencing-plan",
    "external-rank",
    "extract",
    "family-candidate-bridge",
    "family-index",
    "family-to-copy",
    "focal-paralog-import",
    "gff-preflight",
    "guided-profile",
    "homolog-linking-counterfactual",
    "homolog-linking-export",
    "homolog-linking-import",
    "homolog-linking-probe-export",
    "homolog-linking-score",
    "homolog-panel-auto",
    "homolog-panel-benchmark",
    "homolog-panel-match",
    "homology-concordance",
    "hypothesis-ledger",
    "hypothesis-verify",
    "integrated-narrowing",
    "interval-atlas",
    "interval-reverse-lookup",
    "investigation-attribution",
    "liftover-markers",
    "literature",
    "llm-evidence-import",
    "llm-pathway-import",
    "llm-pathway-prompt",
    "llm-prompt-export",
    "llm-review-benchmark",
    "llm-run",
    "llm-shortlist-import",
    "llm-shortlist-prompt",
    "llm-usage",
    "mapping-gate",
    "mechanism-coverage-audit",
    "mechanism-import",
    "multigenome-private",
    "multi-reference-locus-union",
    "next-experiment",
    "omamer-import",
    "orchestrate",
    "orthofinder-import",
    "orthology-import",
    "orthology-plan",
    "orthology-tier1-plan",
    "pathway-interval-falsify",
    "peak-context",
    "phenotype-first-status",
    "phenotype-localization-freeze",
    "phenotype-localization-packet",
    "phenotype-investigate",
    "phenotype-pathway-benchmark",
    "phenotype-variant-benchmark",
    "primer-plan",
    "primerblast-import",
    "promoter-review",
    "propagate-terms",
    "provenance-ledger",
    "qtl-source-init",
    "qtl2-plan",
    "rank",
    "ranking-qa",
    "ranking-stability-benchmark",
    "reference-candidate-table",
    "reference-corpus-build",
    "reference-corpus-compare",
    "reference-corpus-fetch",
    "reference-gene-cache",
    "reference-gene-resolve",
    "reference-panel",
    "research",
    "research-packet",
    "retrieval-hybrid-merge",
    "retrieval-import",
    "retrieval-index-build",
    "retrieval-interval-match",
    "retrieval-merge",
    "retrieval-run",
    "reverse-literature",
    "route-a-interval-match",
    "route-b-prepare",
    "route-b-semantic-bundle",
    "route-b-semantic-candidates",
    "route-b-semantic-import",
    "route-b-semantic-source-verify",
    "route-b-semantic-transfer-verify",
    "route-c-benchmark",
    "route-c-direct-benchmark",
    "route-c-direct-freeze",
    "route-c-direct-members",
    "route-c-direct-project",
    "route-c-direct-search",
    "route-c-provider-bait",
    "route-c-provider-comparison",
    "scaffold-rescue",
    "scenario-create",
    "scenario-report",
    "screening-shortlist",
    "semantic-claims-freeze",
    "source-verify",
    "supported-member-union",
    "sv-plan",
    "synteny-import",
    "trait-ontology-map",
    "trait-terms",
    "tnt1-evidence-normalize",
    "triple-route-compare",
    "variant-evidence",
    "weight-policy-benchmark",
    "wet-candidate-handoff",
    "wet-state",
}


class CommandsRegistryTests(unittest.TestCase):
    def test_registry_covers_every_stable_command(self) -> None:
        self.assertEqual(EXPECTED_COMMANDS, set(COMMAND_NAMES))
        for command in EXPECTED_COMMANDS:
            self.assertIsNotNone(handler_for(command), command)

    def test_register_all_matches_the_handler_registry(self) -> None:
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers(dest="command")
        register_all(subparsers)
        self.assertEqual(EXPECTED_COMMANDS, set(subparsers.choices))

    def test_benchmark_surface_is_a_superset_of_the_core_surface(self) -> None:
        benchmark_commands = set(command_names("benchmark"))
        core_commands = set(command_names("core"))
        self.assertTrue(core_commands)
        self.assertTrue(core_commands < benchmark_commands)
        self.assertEqual(EXPECTED_COMMANDS, benchmark_commands)
        parser = build_parser(surface="core")
        choices = [action for action in parser._subparsers._group_actions if hasattr(action, "choices")][0].choices
        self.assertEqual(core_commands, set(choices))

    def test_duplicate_handler_registration_is_rejected(self) -> None:
        first = SimpleNamespace(HANDLERS={"doctor": lambda args: {}})
        second = SimpleNamespace(HANDLERS={"doctor": lambda args: {}})
        with self.assertRaisesRegex(RuntimeError, "duplicate CLI command handler"):
            _collect_handlers([first, second])

    def test_full_parser_has_exactly_the_stable_commands(self) -> None:
        parser = build_parser(surface="all")
        choices = [action for action in parser._subparsers._group_actions if hasattr(action, "choices")][0].choices
        self.assertEqual(EXPECTED_COMMANDS, set(choices))


if __name__ == "__main__":
    unittest.main()
