"""Canonical feature-reachability gate (issue #373).

The registry maps features claiming pipeline integration to their canonical
CLI entrypoint + fixtures + expected artifacts + manifest keys + rank policy;
these tests pin the fail-closed schema, the runner layer (subprocess
execution + artifact/key assertions), the negative orphan rejection, the
determinism / no-rank-artifacts / fixture-immutability guarantees, the
clean-install fallback guard and the committed registry itself.
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl import reachability
from candidate_qtl.core import sha256
from candidate_qtl.reachability import (
    REQUIRED_INTEGRATED_FIELDS,
    STATUS_CORE_ONLY,
    STATUS_INTEGRATED,
    ReachabilityError,
    load_registry,
    render_argv,
    run_entry,
    run_gate,
    validate_registry,
)

REPO = Path(__file__).resolve().parents[1]
REGISTRY = REPO / "tests" / "fixtures" / "reachability" / "registry.json"
ORPHAN_REGISTRY = REPO / "tests" / "fixtures" / "reachability" / "negative" / "orphan_registry.json"
PYTHON = sys.executable

FAMILY_FIXTURES = [
    "tests/fixtures/reachability/features/family_index/blast_hits.tsv",
    "tests/fixtures/reachability/features/family_index/genes.tsv",
]


def _registry_entry(feature_id: str) -> dict:
    """Look up one committed registry entry by feature_id (registry order is
    alphabetical, not stable API)."""
    return next(item for item in load_registry(REGISTRY)["entries"] if item["feature_id"] == feature_id)


# Declared core_only_non_closing_slice entries: no entrypoint, artifacts or
# fixtures by contract, so execution-based gate tests apply to integrated
# entries only. Wiring this feature belongs to its owning issues, not here.
CORE_ONLY_FEATURES = ["route-b-light-local-research"]


def _integrated_registry_entries() -> list[dict]:
    return [entry for entry in load_registry(REGISTRY)["entries"] if entry["integration_status"] == STATUS_INTEGRATED]


def _integrated_entry(**overrides: object) -> dict:
    """A schema-valid integrated entry (family-index shaped) for negatives."""
    entry: dict = {
        "feature_id": "family-index",
        "integration_status": STATUS_INTEGRATED,
        "entrypoint": {
            "command": "family-index",
            "argv": [
                "--blast-hits",
                "{fixture:tests/fixtures/reachability/features/family_index/blast_hits.tsv}",
                "--genes",
                "{fixture:tests/fixtures/reachability/features/family_index/genes.tsv}",
                "--out-dir",
                "{out}",
            ],
        },
        "fixtures": list(FAMILY_FIXTURES),
        "expected_artifacts": ["family_index.tsv", "family_index_manifest.json"],
        "expected_manifest_keys": [{"artifact": "family_index_manifest.json", "keys": ["schema_version", "gene_count"]}],
        "rank_contribution": "none",
        "rank_assertion": {"type": "deterministic_and_no_rank_artifacts", "reason": "test entry"},
        "clean_install_smoke": True,
        "cache_fingerprint_assertion": False,
        "owning_issue": "#373",
        "doc_reference": "docs/REACHABILITY-GATE.ja.md",
        "assert_values": {"artifact": "family_index_manifest.json", "values": {"gene_count": 6}},
    }
    entry.update(overrides)
    return entry


def _core_only_entry(**overrides: object) -> dict:
    entry = {
        "feature_id": "synthetic-core-only-feature",
        "integration_status": STATUS_CORE_ONLY,
        "remaining_wiring": "wire the canonical CLI entrypoint and fixtures",
        "doc_reference": "docs/REACHABILITY-GATE.ja.md#core-only",
    }
    entry.update(overrides)
    return entry


def _validate(entry: dict) -> list[str]:
    return validate_registry({"schema_version": 1, "entries": [entry]}, repo_root=REPO)


class SchemaValidationTests(unittest.TestCase):
    """Fail-closed registry schema: every negative must be an error."""

    def test_duplicate_feature_id_is_rejected(self) -> None:
        errors = validate_registry(
            {"schema_version": 1, "entries": [_integrated_entry(), _integrated_entry()]}, repo_root=REPO
        )
        self.assertTrue(any("duplicate feature_id" in error for error in errors))

    def test_unknown_integration_status_is_rejected(self) -> None:
        errors = _validate(_integrated_entry(integration_status="integrated_maybe"))
        self.assertTrue(any("unknown integration_status" in error for error in errors))

    def test_unknown_top_level_key_is_rejected(self) -> None:
        errors = _validate(_integrated_entry(mystery_wiring=True))
        self.assertTrue(any("unknown key" in error for error in errors))

    def test_missing_required_integrated_fields_are_rejected(self) -> None:
        for field in sorted(REQUIRED_INTEGRATED_FIELDS):
            with self.subTest(field=field):
                entry = _integrated_entry()
                del entry[field]
                errors = _validate(entry)
                self.assertTrue(
                    any(f"missing required field {field!r}" in error for error in errors),
                    f"expected a missing-field error for {field}",
                )

    def test_absolute_fixture_paths_are_rejected(self) -> None:
        for bad in ("C:\\x\\y.tsv", "/abs.tsv", "\\\\srv\\share"):
            with self.subTest(path=bad):
                errors = _validate(_integrated_entry(fixtures=[bad]))
                self.assertTrue(any("repo-relative" in error for error in errors), bad)

    def test_backslash_fixture_path_is_rejected(self) -> None:
        errors = _validate(_integrated_entry(fixtures=["tests\\fixtures\\x.tsv"]))
        self.assertTrue(any("forward slashes" in error for error in errors))

    def test_dotdot_fixture_path_is_rejected(self) -> None:
        errors = _validate(_integrated_entry(fixtures=["features/../secret.tsv"]))
        self.assertTrue(any(".. component" in error for error in errors))

    def test_fixture_escaping_the_repo_root_is_rejected(self) -> None:
        errors = _validate(_integrated_entry(fixtures=["tests/fixtures/../../secret"]))
        self.assertTrue(any(".. component" in error for error in errors))

    def test_nonexistent_fixture_is_rejected_when_repo_root_given(self) -> None:
        errors = _validate(_integrated_entry(fixtures=["tests/fixtures/reachability/features/family_index/nope.tsv"]))
        self.assertTrue(any("fixture does not exist" in error for error in errors))

    def test_unknown_command_name_is_rejected(self) -> None:
        errors = _validate(_integrated_entry(entrypoint={"command": "no-such-command", "argv": ["--out-dir", "{out}"]}))
        self.assertTrue(any("unknown command" in error for error in errors))

    def test_unsorted_entries_are_rejected(self) -> None:
        first = _integrated_entry(feature_id="zeta-feature")
        second = _integrated_entry(feature_id="alpha-feature")
        errors = validate_registry({"schema_version": 1, "entries": [first, second]}, repo_root=REPO)
        self.assertTrue(any("sorted by feature_id" in error for error in errors))

    def test_invalid_rank_assertion_type_is_rejected(self) -> None:
        errors = _validate(
            _integrated_entry(rank_assertion={"type": "maybe_deterministic", "reason": "test"})
        )
        self.assertTrue(any("rank_assertion.type" in error for error in errors))

    def test_empty_rank_assertion_reason_is_rejected(self) -> None:
        errors = _validate(_integrated_entry(rank_assertion={"type": "deterministic_and_no_rank_artifacts", "reason": "  "}))
        self.assertTrue(any("rank_assertion.reason" in error for error in errors))

    def test_invalid_rank_contribution_is_rejected(self) -> None:
        errors = _validate(_integrated_entry(rank_contribution="sometimes"))
        self.assertTrue(any("rank_contribution" in error for error in errors))
        self.assertEqual([], _validate(_integrated_entry(rank_contribution="policy:dry-v1")))

    def test_core_only_entry_with_entrypoint_is_rejected(self) -> None:
        entry = _core_only_entry(entrypoint={"command": "family-index", "argv": ["--out-dir", "{out}"]})
        errors = _validate(entry)
        self.assertTrue(any("must not carry entrypoint" in error for error in errors))

    def test_core_only_entry_with_expected_artifacts_is_rejected(self) -> None:
        errors = _validate(_core_only_entry(expected_artifacts=["family_index.tsv"]))
        self.assertTrue(any("must not carry expected_artifacts" in error for error in errors))

    def test_core_only_missing_remaining_wiring_is_rejected(self) -> None:
        entry = _core_only_entry()
        del entry["remaining_wiring"]
        errors = _validate(entry)
        self.assertTrue(any("remaining_wiring" in error for error in errors))

    def test_core_only_missing_doc_reference_is_rejected(self) -> None:
        entry = _core_only_entry()
        del entry["doc_reference"]
        errors = _validate(entry)
        self.assertTrue(any("doc_reference" in error for error in errors))

    def test_registry_root_not_an_object_is_rejected(self) -> None:
        for bad in (["not", "a", "dict"], "not-a-dict"):
            with self.subTest(bad=bad):
                errors = validate_registry(bad)  # type: ignore[arg-type]
                self.assertEqual(["registry must be a JSON object"], errors)

    def test_integrated_entry_with_empty_argv_is_rejected(self) -> None:
        errors = _validate(_integrated_entry(entrypoint={"command": "family-index", "argv": []}))
        self.assertTrue(any("argv must be a non-empty list" in error for error in errors))

    def test_undeclared_fixture_placeholder_is_rejected(self) -> None:
        entry = _integrated_entry()
        entry["entrypoint"]["argv"] = [
            "--blast-hits",
            "{fixture:tests/fixtures/reachability/features/route_a_interval_match/hypotheses.tsv}",
        ]
        errors = _validate(entry)
        self.assertTrue(any("undeclared fixture" in error for error in errors))

    def test_unknown_placeholder_is_rejected(self) -> None:
        entry = _integrated_entry()
        entry["entrypoint"]["argv"] = ["--blast-hits", "{blah}"]
        errors = _validate(entry)
        self.assertTrue(any("unknown argv placeholder" in error for error in errors))

    def test_inline_placeholder_form_is_rejected(self) -> None:
        entry = _integrated_entry()
        entry["entrypoint"]["argv"] = ["--blast-hits", "x.tsv", "--out-dir={out}"]
        errors = _validate(entry)
        self.assertTrue(any("full-match token" in error for error in errors))

    def test_argv_with_non_string_element_is_rejected(self) -> None:
        entry = _integrated_entry()
        entry["entrypoint"]["argv"] = ["--blast-hits", 42]
        errors = _validate(entry)
        self.assertTrue(any("non-empty strings" in error for error in errors))

    def test_duplicate_expected_artifact_is_rejected(self) -> None:
        errors = _validate(_integrated_entry(expected_artifacts=["family_index.tsv", "family_index.tsv"]))
        self.assertTrue(any("duplicate expected artifact" in error for error in errors))

    def test_expected_artifact_must_be_safe_out_dir_relative_path(self) -> None:
        for bad in ("/absolute.tsv", "C:/absolute.tsv", "../escape.tsv", "rank/../escape.tsv", "rank\\output.tsv"):
            with self.subTest(path=bad):
                errors = _validate(_integrated_entry(expected_artifacts=[bad]))
                self.assertTrue(any("expected artifact" in error for error in errors), errors)

    def test_manifest_key_dict_shape_errors(self) -> None:
        for item, expected_fragment in (
            ({"keys": ["status"]}, "artifact name"),
            ({"artifact": "", "keys": ["status"]}, "artifact name"),
            ({"artifact": "x.json"}, "key paths"),
            ({"artifact": "x.json", "keys": []}, "key paths"),
            ({"artifact": "x.json", "keys": [""]}, "key paths"),
        ):
            with self.subTest(item=item):
                errors = _validate(_integrated_entry(expected_manifest_keys=[item]))
                self.assertTrue(any(expected_fragment in error for error in errors), item)

    def test_manifest_key_dict_artifact_not_declared_is_rejected(self) -> None:
        entry = _integrated_entry()
        entry["expected_manifest_keys"] = [{"artifact": "typo_manifest.json", "keys": ["status"]}]
        errors = _validate(entry)
        self.assertTrue(
            any("manifest-key artifact 'typo_manifest.json' must be an exact declared JSON artifact" in error for error in errors),
        )

    def test_manifest_keys_must_bind_an_exact_declared_json_artifact(self) -> None:
        for item in (
            "schema_version",
            {"artifact": "family_index.tsv", "keys": ["schema_version"]},
            {"artifact": "*.json", "keys": ["schema_version"]},
            {"artifact": "../family_index_manifest.json", "keys": ["schema_version"]},
        ):
            with self.subTest(item=item):
                errors = _validate(_integrated_entry(expected_manifest_keys=[item]))
                self.assertTrue(any("manifest" in error and "artifact" in error for error in errors), errors)

    def test_assert_values_must_bind_the_declared_json_manifest(self) -> None:
        errors = _validate(_integrated_entry(assert_values={"artifact": "other.json", "values": {"gene_count": 6}}))
        self.assertTrue(any("assert_values.artifact" in error for error in errors), errors)

    def test_expected_exit_codes_shape_errors(self) -> None:
        for bad in ("0", [], [0, "x"]):
            with self.subTest(bad=bad):
                errors = _validate(_integrated_entry(expected_exit_codes=bad))
                self.assertTrue(any("expected_exit_codes" in error for error in errors), bad)

    def test_boolean_declarations_must_be_bool(self) -> None:
        errors = _validate(_integrated_entry(clean_install_smoke="yes"))
        self.assertTrue(any("clean_install_smoke must be a boolean" in error for error in errors))
        errors = _validate(_integrated_entry(cache_fingerprint_assertion=1))
        self.assertTrue(any("cache_fingerprint_assertion must be a boolean" in error for error in errors))

    def test_empty_owning_issue_is_rejected(self) -> None:
        errors = _validate(_integrated_entry(owning_issue="   "))
        self.assertTrue(any("owning_issue must be a non-empty string" in error for error in errors))

    def test_empty_doc_reference_is_rejected(self) -> None:
        errors = _validate(_integrated_entry(doc_reference=""))
        self.assertTrue(any("doc_reference must be a non-empty string" in error for error in errors))

    def test_assert_values_wrong_shape_is_rejected(self) -> None:
        errors = _validate(_integrated_entry(assert_values=["schema_version"]))
        self.assertTrue(any("assert_values must be an object" in error for error in errors))

    def test_duplicate_fixture_is_rejected(self) -> None:
        errors = _validate(_integrated_entry(fixtures=[FAMILY_FIXTURES[0], FAMILY_FIXTURES[0]]))
        self.assertTrue(any("duplicate fixture path" in error for error in errors))

    def test_bad_schema_version_is_rejected(self) -> None:
        errors = validate_registry({"schema_version": 99, "entries": [_integrated_entry()]}, repo_root=REPO)
        self.assertTrue(any("schema_version" in error for error in errors))

    def test_entries_missing_is_rejected(self) -> None:
        self.assertTrue(any("entries must be a list" in error for error in validate_registry({"schema_version": 1})))

    def test_error_messages_carry_the_feature_id(self) -> None:
        errors = _validate(_integrated_entry(feature_id="bad-entry", integration_status="nope"))
        self.assertTrue(all("bad-entry" in error for error in errors))


class CoreOnlySemanticsTests(unittest.TestCase):
    """core_only declarations record remaining wiring; they never count as coverage."""

    def test_valid_core_only_entry_returns_no_errors(self) -> None:
        self.assertEqual([], _validate(_core_only_entry()))

    def test_run_gate_counts_but_does_not_execute_core_only_entries(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            integrated = _registry_entry("family-index")
            registry = root / "registry.json"
            registry.write_text(
                json.dumps({"schema_version": 1, "entries": [integrated, _core_only_entry()]}),
                encoding="utf-8",
            )
            out_root = root / "gate"
            result = run_gate(registry, REPO, "source", PYTHON, out_root=out_root)
            self.assertTrue(result["valid"])
            self.assertEqual(1, result["integrated_count"])
            self.assertEqual(1, result["core_only_count"])
            self.assertTrue((out_root / "family-index").exists())
            self.assertFalse((out_root / "synthetic-core-only-feature").exists())

    def test_registry_with_zero_integrated_entries_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            registry = Path(tmp) / "registry.json"
            registry.write_text(
                json.dumps({"schema_version": 1, "entries": [_core_only_entry()]}), encoding="utf-8"
            )
            with self.assertRaises(ReachabilityError) as ctx:
                run_gate(registry, REPO, "source", PYTHON)
            self.assertIn("no integrated entries", str(ctx.exception))


class OrphanAndRunnerLayerTests(unittest.TestCase):
    """A schema-shaped entry with no canonical wiring must be rejected, and the
    runner layer must fail on missing artifacts and missing manifest keys."""

    def test_synthetic_orphan_is_rejected_with_feature_id_and_command(self) -> None:
        data = load_registry(ORPHAN_REGISTRY)
        errors = validate_registry(data, repo_root=REPO)
        self.assertTrue(errors)
        combined = "; ".join(errors)
        self.assertIn("synthetic-orphan-feature", combined)
        self.assertIn("synthetic-orphan-cmd", combined)

    def test_real_cli_command_without_canonical_orchestration_is_rejected(self) -> None:
        entry = _integrated_entry(feature_id="synthetic-real-command")
        errors = _validate(entry)
        combined = "; ".join(errors)
        self.assertIn("synthetic-real-command", combined)
        self.assertIn("not wired into canonical orchestration", combined)

    def test_missing_expected_artifact_raises_naming_the_artifact(self) -> None:
        entry = _registry_entry("family-index")
        entry = json.loads(json.dumps(entry))
        entry["expected_artifacts"] = ["never_written.tsv"]
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ReachabilityError) as ctx:
                run_entry(entry, REPO, Path(tmp) / "out", PYTHON, cwd=REPO)
            self.assertEqual("family-index", ctx.exception.feature_id)
            self.assertIn("never_written.tsv", str(ctx.exception))

    def test_missing_manifest_key_raises_naming_the_key(self) -> None:
        entry = _registry_entry("family-index")
        entry = json.loads(json.dumps(entry))
        entry["expected_manifest_keys"] = [{"artifact": "family_index_manifest.json", "keys": ["bogus.path.0.key"]}]
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ReachabilityError) as ctx:
                run_entry(entry, REPO, Path(tmp) / "out", PYTHON, cwd=REPO)
            self.assertEqual("family-index", ctx.exception.feature_id)
            self.assertIn("bogus.path.0.key", str(ctx.exception))

    def test_assert_values_mismatch_raises_naming_the_path(self) -> None:
        entry = _registry_entry("family-index")
        entry = json.loads(json.dumps(entry))
        entry["assert_values"] = {"artifact": "family_index_manifest.json", "values": {"gene_count": 999}}
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ReachabilityError) as ctx:
                run_entry(entry, REPO, Path(tmp) / "out", PYTHON, cwd=REPO)
            self.assertEqual("family-index", ctx.exception.feature_id)
            self.assertIn("gene_count", str(ctx.exception))

    def test_manifest_key_dict_artifact_typo_fails_validation_and_run(self) -> None:
        entry = _registry_entry("family-index")
        entry = json.loads(json.dumps(entry))
        entry["expected_artifacts"] = ["family_index.tsv", "family_index_manifest.json"]
        entry["expected_manifest_keys"] = [{"artifact": "typo_manifest.json", "keys": ["status"]}]
        errors = validate_registry({"schema_version": 1, "entries": [entry]}, repo_root=REPO)
        self.assertTrue(
            any("manifest-key artifact 'typo_manifest.json' must be an exact declared JSON artifact" in error for error in errors)
        )
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ReachabilityError) as ctx:
                run_entry(entry, REPO, Path(tmp) / "out", PYTHON, cwd=REPO)
            self.assertEqual("family-index", ctx.exception.feature_id)
            self.assertIn("typo_manifest.json", str(ctx.exception))

    def test_render_argv_resolves_placeholders_to_absolute_paths(self) -> None:
        entry = _registry_entry("family-index")
        out_dir = Path(tempfile.mkdtemp(prefix="reachability-argv-")) / "out"
        rendered = render_argv(entry, REPO, out_dir)
        flags = {"--blast-hits", "--genes", "--out-dir"}
        positional = [item for item in rendered if item not in flags]
        self.assertEqual([str(out_dir)], positional[-1:])
        self.assertTrue(all(Path(item).is_absolute() for item in positional))


class PositiveE2EGateTests(unittest.TestCase):
    """The committed registry must pass the gate from a source checkout."""

    def test_committed_registry_passes_source_mode_gate(self) -> None:
        result = _cached_gate()
        self.assertTrue(result["valid"])
        self.assertEqual(10, result["integrated_count"])
        # Core-only declarations are explicit non-closing slices, never counted
        # as integrated coverage. Pin the declared roster so a new one cannot
        # appear silently.
        self.assertEqual(len(CORE_ONLY_FEATURES), result["core_only_count"])
        self.assertEqual(CORE_ONLY_FEATURES, sorted(entry["feature_id"] for entry in load_registry(REGISTRY)["entries"]
                                                    if entry["integration_status"] == STATUS_CORE_ONLY))
        self.assertEqual([], [entry for entry in result["entries"] if entry.get("error")])
        by_id = {entry["feature_id"]: entry for entry in result["entries"]}
        self.assertEqual(
            ["budget-run", "family-candidate-bridge",
             "family-index", "mechanism-import", "mechanism-import-coverage", "route-a-interval-match",
             "route-b-research-refreeze", "tier-bibliography-budget-a",
             "tier-bibliography-budget-b", "triple-route-e2e"],
            sorted(by_id),
        )
        for entry in result["entries"]:
            self.assertEqual(0, entry["exit_code"])
            self.assertTrue(entry["artifacts"])
            self.assertTrue(entry["manifest_keys"])

    def test_tier_bibliography_and_refreeze_reach_the_real_cli(self) -> None:
        by_id = {entry["feature_id"]: entry for entry in _cached_gate()["entries"]}
        for feature in ("tier-bibliography-budget-a", "tier-bibliography-budget-b"):
            out = Path(by_id[feature]["out_dir"])
            selection = json.loads((out / "bibliography_selection.json").read_text())["panels"][0]
            self.assertEqual((7, 5), (selection["candidate_paper_count"], selection["selected_paper_count"]))
            self.assertIn("PMID:107", selection["selected_source_ids"])
        b = pd.read_csv(Path(by_id["tier-bibliography-budget-b"]["out_dir"]) / "semantic_claim_source_verification.tsv", sep="\t")
        self.assertEqual(7, len(b))
        self.assertEqual(2, b.bibliography_state.eq("metadata_budget_not_selected").sum())
        parent = json.loads((REPO / "tests/fixtures/reachability/features/tier_bibliography_budget/parent.json").read_text())
        new = json.loads((Path(by_id["route-b-research-refreeze"]["out_dir"]) / "route_b_semantic_bundle.json").read_text())
        for old_unit, new_unit in zip(parent["units"], new["units"], strict=True):
            for key in ("genes", "identity", "member_interval_genes"):
                self.assertEqual(old_unit[key], new_unit[key])
        self.assertIn("PMID:107", new["units"][0]["stable_source_ids"])

    def test_mechanism_import_output_preserves_legacy_values(self) -> None:
        by_id = {entry["feature_id"]: entry for entry in _cached_gate()["entries"]}
        out_dir = Path(by_id["mechanism-import"]["out_dir"])
        hypotheses = (out_dir / "mechanism_hypotheses.tsv").read_text(encoding="utf-8")
        self.assertIn("Chalcone synthase (CHS)", hypotheses)
        self.assertIn("mechanism_statement", hypotheses)
        self.assertIn("selection_rationale", hypotheses)
        adapter_map = (out_dir / "mechanism_payload_adapter_map.tsv").read_text(encoding="utf-8")
        self.assertIn("pea_a_M1", adapter_map)
        self.assertIn("corpus_ortholog_model_plant", adapter_map)
        manifest = json.loads((out_dir / "mechanism_import_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual("legacy_v1_7col", manifest["accepted_payload_schema"])
        self.assertEqual(3, manifest["rows_imported"])

    def test_tier_bibliography_outputs_repeat_without_membership_changes(self) -> None:
        for feature in ("tier-bibliography-budget-a", "tier-bibliography-budget-b", "route-b-research-refreeze"):
            with self.subTest(feature=feature):
                first, second = _fresh_runs(feature)
                snapshots = [{path.relative_to(root).as_posix(): path.read_bytes() for path in root.rglob("*.tsv")} for root in (first, second)]
                self.assertTrue(snapshots[0])
                self.assertEqual(snapshots[0], snapshots[1])

    def test_budget_run_manifest_reports_completed_zero_profile(self) -> None:
        by_id = {entry["feature_id"]: entry for entry in _cached_gate()["entries"]}
        out_dir = Path(by_id["budget-run"]["out_dir"])
        manifest = json.loads((out_dir / "budget_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual("1.0", manifest["schema_version"])
        self.assertEqual("completed", manifest["status"])
        self.assertEqual("zero", manifest["budget_profile"])
        self.assertTrue(manifest["budgets"]["sequence_query_limit"]["exhausted"])
        self.assertEqual("tier0_inventory.tsv", manifest["outputs"]["tier0_inventory_tsv"])
        self.assertEqual("tier1_clusters.tsv", manifest["outputs"]["tier1_clusters_tsv"])
        handoff = (out_dir / "budget_handoff.tsv").read_text(encoding="utf-8")
        self.assertIn("ShSBEI", (out_dir / "tier0_inventory.tsv").read_text(encoding="utf-8"))
        self.assertIn("rescue_pending_review", handoff)
        self.assertIn("not_evaluated_provider_missing", handoff)

    def test_family_index_output_is_substantive(self) -> None:
        by_id = {entry["feature_id"]: entry for entry in _cached_gate()["entries"]}
        out_dir = Path(by_id["family-index"]["out_dir"])
        index = (out_dir / "family_index.tsv").read_text(encoding="utf-8")
        self.assertGreater(len(index), 0)
        self.assertIn("CLV1", index)
        self.assertIn("reference_genome", index)

    def test_family_candidate_bridge_output_carries_the_mapping_tiers(self) -> None:
        by_id = {entry["feature_id"]: entry for entry in _cached_gate()["entries"]}
        out_dir = Path(by_id["family-candidate-bridge"]["out_dir"])
        bridge = (out_dir / "family_candidate_bridge.tsv").read_text(encoding="utf-8")
        self.assertIn("tier1_source_backed_unique", bridge)
        self.assertIn("tier5_name_only", bridge)
        self.assertIn("rank_contribution", bridge)
        manifest = json.loads((out_dir / "family_candidate_bridge_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual("none", manifest["rank_contribution"])
        self.assertEqual(2, manifest["row_counts"]["total"])

    def test_route_a_output_carries_the_concrete_match_status(self) -> None:
        by_id = {entry["feature_id"]: entry for entry in _cached_gate()["entries"]}
        out_dir = Path(by_id["route-a-interval-match"]["out_dir"])
        status = (out_dir / "route_a_interval_family_match.tsv").read_text(encoding="utf-8")
        self.assertIn("concrete_interval_gene_matched", status)
        self.assertIn("rank_contribution", status)
        self.assertIn("none", status)

    @unittest.skip("public release: discovery-chain-benchmark is a benchmark command (not shipped)")
    def test_discovery_chain_summary_reports_the_expected_recall(self) -> None:
        by_id = {entry["feature_id"]: entry for entry in _cached_gate()["entries"]}
        out_dir = Path(by_id["discovery-chain-benchmark"]["out_dir"])
        summary = (out_dir / "discovery_chain_summary.tsv").read_text(encoding="utf-8")
        self.assertIn("a_only_recall\t1\t1\t1.0", summary)
        self.assertIn("convergence_recall\t1\t1\t1.0", summary)

    def test_gate_fails_closed_on_all_validation_errors(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bad = root / "bad.json"
            bad.write_text(
                json.dumps({"schema_version": 99, "entries": [_core_only_entry()]}), encoding="utf-8"
            )
            with self.assertRaises(ReachabilityError) as ctx:
                run_gate(bad, REPO, "source", PYTHON)
            self.assertIn("schema_version", str(ctx.exception))


class DeterminismTests(unittest.TestCase):
    """Every declared TSV artifact must be byte-identical across runs.

    Three independent executions per entry (the module-cached gate run plus
    two fresh runs into separate out-dirs) are compared pairwise, so a
    same-second timestamp or counter cannot slip through a two-run check.
    """

    def test_all_tsv_artifacts_are_byte_identical_across_three_runs(self) -> None:
        cached = {entry["feature_id"]: entry for entry in _cached_gate()["entries"]}
        for entry in _integrated_registry_entries():
            with self.subTest(feature_id=entry["feature_id"]):
                tsv_names = [name for name in entry["expected_artifacts"] if name.endswith(".tsv")]
                self.assertTrue(tsv_names, entry["feature_id"])
                first, second = _fresh_runs(entry["feature_id"])
                third = Path(cached[entry["feature_id"]]["out_dir"])
                for name in tsv_names:
                    for run_a, run_b in ((first, second), (first, third), (second, third)):
                        self.assertEqual(
                            (run_a / name).read_bytes(),
                            (run_b / name).read_bytes(),
                            f"{entry['feature_id']}: {name} differs between {run_a} and {run_b}",
                        )

    def test_committed_route_a_family_index_matches_canonical_family_index_output(self) -> None:
        """The committed static family_index.tsv must be byte-identical to a
        fresh family-index run over the committed blast_hits/genes fixtures."""
        first, _ = _fresh_runs("family-index")
        committed = REPO / "tests" / "fixtures" / "reachability" / "features" / "route_a_interval_match" / "family_index.tsv"
        self.assertEqual(
            (first / "family_index.tsv").read_bytes(),
            committed.read_bytes(),
        )


class NoRankArtifactsAndFixtureImmutabilityTests(unittest.TestCase):
    """rank_contribution=none entries must write no rank artifact and may not
    modify their declared fixtures."""

    def test_rank_impacting_pattern_matcher(self) -> None:
        for path, expected in {
            "rank/ranked_candidates.tsv": True,
            "score_contributions.tsv": True,
            "candidates.tsv": True,
            "flagged_candidates.tsv": True,
            "stage_a/anything.tsv": True,
            "dry_rank.tsv": True,
            "ranked_candidates_top10.tsv": True,
            "family_index.tsv": False,
            "surviving_candidates.tsv": False,
            "hypothesis_interval_falsification.tsv": False,
            "discovery_chain_cases.tsv": False,
        }.items():
            with self.subTest(path=path):
                self.assertEqual(expected, reachability._matches_rank_impacting(path))

    def test_gate_outputs_contain_no_rank_artifacts(self) -> None:
        registry = {entry["feature_id"]: entry for entry in load_registry(REGISTRY)["entries"]}
        for entry in _cached_gate()["entries"]:
            # rank_contribution=policy entries run the canonical rank stage and
            # are asserted by the companion test below; the no-rank-artifact
            # contract binds entries that declare rank_contribution=none.
            if registry[entry["feature_id"]]["rank_contribution"] != "none":
                continue
            out_dir = Path(entry["out_dir"])
            offending = [
                path.relative_to(out_dir).as_posix()
                for path in out_dir.rglob("*")
                if path.is_file() and reachability._matches_rank_impacting(path.relative_to(out_dir).as_posix())
            ]
            self.assertEqual([], offending, entry["feature_id"])

    def test_triple_route_e2e_rank_outputs_are_confined_to_the_case_run_subtree(self) -> None:
        """The triple-route-e2e policy entry emits the canonical deterministic
        rank only under its z03_case_run subtree (the case-run stage), never
        scattered elsewhere in the out-dir: the family layers between the
        bridge and case-run keep writing rank_contribution=none artifacts."""
        by_id = {entry["feature_id"]: entry for entry in _cached_gate()["entries"]}
        out_dir = Path(by_id["triple-route-e2e"]["out_dir"])
        offending = [
            path.relative_to(out_dir).as_posix()
            for path in out_dir.rglob("*")
            if path.is_file() and reachability._matches_rank_impacting(path.relative_to(out_dir).as_posix())
        ]
        self.assertTrue(offending, "the E2E entry must emit the canonical rank table")
        for relative in offending:
            self.assertTrue(relative.startswith("z03_case_run/"), relative)
        for canonical in ("z03_case_run/rank/ranked_candidates.tsv",
                          "z03_case_run/rank/score_contributions.tsv",
                          "z03_case_run/stage_a/deterministic_shortlist.tsv"):
            self.assertTrue((out_dir / canonical).is_file(), canonical)

    def test_triple_route_e2e_bridge_reaches_family_copy(self) -> None:
        """Issue #374 acceptance: one guided-profile triple-route run produces
        the full bridge -> family-to-copy -> integrated narrowing -> stage A
        narrowed -> next-experiment chain, and family-to-copy consumed the
        bridge mapping tier (family_copy_manifest.bridge.supplied=true)."""
        by_id = {entry["feature_id"]: entry for entry in _cached_gate()["entries"]}
        out_dir = Path(by_id["triple-route-e2e"]["out_dir"])
        guided = json.loads((out_dir / "guided_profile_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual("partial", guided["overall_status"])
        self.assertFalse(guided["final"])
        self.assertEqual("completed", guided["legacy_profile"]["overall_status"])
        self.assertTrue(guided["legacy_profile"]["final"])
        canonical = guided["canonical_execution"]
        self.assertEqual("semantic_review_required", canonical["state"])
        self.assertTrue(canonical["required_work_remaining"])
        self.assertIsNone(canonical["wet_candidates"])
        self.assertTrue(Path(canonical["expected_semantic_run_dir"]).parent.exists())
        self.assertEqual(sha256(Path(canonical["manifest"]["path"])), canonical["manifest"]["sha256"])
        self.assertEqual("explicit", guided["family_index_resolution"]["status"])
        route_a_status = pd.read_csv(
            out_dir / "a05b_interval_family_match" / "route_a_interval_family_match.tsv", sep="\t", dtype=str
        ).fillna("")
        self.assertTrue(route_a_status["evaluable"].eq("True").all())
        self.assertIn("Psat_CLV1", set(route_a_status["interval_gene_ids"]))
        by_command = {record["command"]: record for record in guided["execution_results"]}
        for command in ("family-candidate-bridge", "family-to-copy", "case-run"):
            self.assertEqual("completed", by_command[command]["status"], command)
        copy_manifest = json.loads((out_dir / "z02_family_to_copy" / "family_copy_manifest.json")
                                   .read_text(encoding="utf-8"))
        self.assertTrue(copy_manifest["bridge"]["supplied"])
        self.assertEqual("case_fas", copy_manifest["case_id"])
        bridge = json.loads((out_dir / "z01b_family_candidate_bridge" / "family_candidate_bridge_manifest.json")
                            .read_text(encoding="utf-8"))
        self.assertEqual("none", bridge["rank_contribution"])
        self.assertGreaterEqual(bridge["row_counts"]["total"], 1)
        for artifact in (
            "z03_case_run/integrated_narrowing/integrated_narrowing.tsv",
            "z03_case_run/stage_a_narrowed/deterministic_shortlist.tsv",
            "z03_case_run/next_experiment/next_experiment_options.tsv",
        ):
            self.assertTrue((out_dir / artifact).is_file(), artifact)

    def test_declared_fixtures_are_unchanged_after_the_gate(self) -> None:
        gate = _cached_gate()
        by_id = {entry["feature_id"]: entry for entry in gate["entries"]}
        for entry in _integrated_registry_entries():
            for fixture in entry["fixtures"]:
                recorded = by_id[entry["feature_id"]]["fixture_hashes"][fixture]
                self.assertEqual(
                    recorded,
                    reachability._fixture_sha256(REPO, fixture),
                    f"{entry['feature_id']}: fixture changed {fixture}",
                )

    def test_canonical_case_run_baseline_preserves_rank_score_candidate_and_stage_a(self) -> None:
        baseline = _cached_gate()["canonical_baseline"]
        self.assertIsNone(baseline["error"])
        self.assertEqual(
            list(reachability.CANONICAL_BASELINE_ARTIFACTS),
            baseline["artifacts"],
        )
        self.assertEqual(baseline["before"], baseline["after"])
        self.assertEqual(
            set(reachability.CANONICAL_BASELINE_ARTIFACTS),
            set(baseline["before"]),
        )


class FeatureEnvTests(unittest.TestCase):
    """Feature subprocesses must never inherit a caller PYTHONPATH."""

    def test_feature_env_removes_pythonpath_from_the_process_environment(self) -> None:
        with patch.dict(os.environ, {"PYTHONPATH": "/some/other/dir"}, clear=False):
            env = reachability._feature_env(None)
            self.assertNotIn("PYTHONPATH", env)

    def test_feature_env_merges_callers_and_still_removes_pythonpath(self) -> None:
        env = reachability._feature_env({"PYTHONPATH": "/x", "CANDIDATE_QTL_CUSTOM": "1"})
        self.assertNotIn("PYTHONPATH", env)
        self.assertEqual("1", env["CANDIDATE_QTL_CUSTOM"])


class CleanInstallFallbackGuardTests(unittest.TestCase):
    """assert_no_source_fallback must reject an interpreter that resolves
    candidate_qtl into the source checkout (editable install / PYTHONPATH)."""

    @unittest.skipUnless(importlib.util.find_spec("candidate_qtl"), "candidate_qtl not installed")
    def test_dev_python_cannot_resolve_into_the_source_checkout(self) -> None:
        try:
            resolved = reachability.assert_no_source_fallback(PYTHON, REPO)
        except ReachabilityError as error:
            # Editable / PYTHONPATH-based dev environments: the guard fires.
            self.assertEqual("<clean-install>", error.feature_id)
            message = str(error)
            self.assertTrue("source checkout" in message or "site-packages" in message, message)
            return
        # A genuinely non-editable dev environment: the guard's contract is
        # still satisfied -- the resolved path is in site-packages and outside
        # the source checkout.
        self.assertIn("site-packages", resolved)
        self.assertFalse(Path(resolved).resolve().is_relative_to(REPO.resolve()))


class CommittedRegistrySelfCheckTests(unittest.TestCase):
    """The committed registry must stay valid; this test fails CI on regressions."""

    def test_committed_registry_is_valid(self) -> None:
        self.assertEqual([], validate_registry(load_registry(REGISTRY), repo_root=REPO))


_cache: dict = {}


def _cached_gate() -> dict:
    """Run the committed registry gate once per test module and reuse it."""
    if "result" not in _cache:
        out_root = Path(tempfile.mkdtemp(prefix="reachability-e2e-"))
        _cache["result"] = run_gate(REGISTRY, REPO, "source", PYTHON, out_root=out_root)
    return _cache["result"]


def _fresh_runs(feature_id: str) -> tuple[Path, Path]:
    """Two fresh independent runs of one entry, cached per test module."""
    if feature_id not in _cache:
        entry = next(item for item in load_registry(REGISTRY)["entries"] if item["feature_id"] == feature_id)
        first = Path(tempfile.mkdtemp(prefix=f"reachability-fresh-{feature_id}-")) / "first"
        second = Path(tempfile.mkdtemp(prefix=f"reachability-fresh-{feature_id}-")) / "second"
        run_entry(entry, REPO, first, PYTHON, cwd=REPO)
        run_entry(entry, REPO, second, PYTHON, cwd=REPO)
        _cache[feature_id] = (first, second)
    return _cache[feature_id]


if __name__ == "__main__":
    unittest.main()
