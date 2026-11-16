"""Issue #363 E2E: case-run wires A/B discordance into the causal scenario ledger.

The ``causal_scenarios`` step appends review-required discordance rows whenever
``inputs.dual_route_comparison_resource`` is registered. The comparison input
never changes rank / stage-A / score outputs, the ledger is rebuilt from scratch
every run (resume cannot duplicate scenarios), and a route-B manifest declaring
``not_evaluated_*`` semantic status collapses every comparison row into one
``provider_missing`` / ``not_evaluated`` row instead of mixing into ``neither``.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import tempfile
import unittest
from pathlib import Path

import pandas as pd
import yaml

from candidate_qtl.causal_scenarios import SCENARIO_COLUMNS
from candidate_qtl.cli import main
from candidate_qtl.discordance_scenarios import (
    PROVIDER_MISSING_SCENARIO_TYPE,
    comparison_fingerprint,
)

COMPARISON_HEADER = (
    "comparison_id\troute_a_hypothesis_id\troute_b_hypothesis_id\tconvergence_status\t"
    "route_a_source_ids\troute_b_source_ids\tshared_source_ids"
)

FOUR_STATE_ROWS = (
    COMPARISON_HEADER + "\n"
    "cmp_000\tha1\thb1\tconvergent\tSRC1\tSRC2\tSRC1\n"
    "cmp_001\tha2\t\tphenotype_only_supported\tSRC3\t\t\n"
    "cmp_002\t\thb3\thypothesis_only_supported\t\tSRC4\t\n"
    "cmp_003\t\t\tno_support\t\t\t\n"
)

SCENARIO_ROWS = (
    "case_id\tscenario_type\tsupporting_evidence\tcontradicting_evidence\tassessment_state\t"
    "recommended_next_experiment\treason\n"
    "bundled_fas_like\tcoding_gene\tcandidate coding variant\t\tsupported\t"
    "sequence candidate CDS\tcoding hypothesis remains plausible\n"
    "bundled_fas_like\tregulatory_element\t\t\tunassessed\tpromoter review\t"
    "noncoding space not evaluated\n"
)


class CaseRunDiscordanceScenarioWiringTests(unittest.TestCase):
    """End-to-end case-run wiring of the issue #363 adapter, per PR #372 harness."""

    def _work(self, root: Path) -> Path:
        fixture = Path(__file__).parent / "fixtures" / "case_run"
        work = root / "case"
        shutil.copytree(fixture, work)
        candidates = pd.read_csv(work / "candidates.tsv", sep="\t", dtype=str).fillna("")
        candidates["strand"] = "+"
        extra = pd.DataFrame([{
            "gene_id": "GeneExcl", "chrom": "chr1", "start": "800", "end": "980",
            "annotation_quality": "complete_orf", "protein_length": "50", "strand": "+",
        }])
        pd.concat([candidates, extra], ignore_index=True).to_csv(
            work / "candidates.tsv", sep="\t", index=False, lineterminator="\n")
        return work

    def _write(self, work: Path, name: str, content: str) -> Path:
        path = work / name
        path.write_text(content, encoding="utf-8")
        return path

    def _run_case(self, work: Path, out: Path, *, inputs: dict | None = None,
                  resources: dict | None = None, resume: bool = False) -> Path:
        case = yaml.safe_load((work / "case.yaml").read_text(encoding="utf-8"))
        if inputs:
            case["inputs"].update(inputs)
        case_path = work / "case_run.yaml"
        case_path.write_text(yaml.safe_dump(case, sort_keys=False), encoding="utf-8")
        registry = yaml.safe_load((work / "resources.yaml").read_text(encoding="utf-8"))
        if resources:
            registry["resources"].update(resources)
        (work / "resources.yaml").write_text(yaml.safe_dump(registry, sort_keys=False), encoding="utf-8")
        args = ["case-run", "--case", str(case_path), "--resources", str(work / "resources.yaml"),
                "--out-dir", str(out)]
        if resume:
            args.append("--resume")
        self.assertEqual(0, main(args))
        return out

    def _comparison_resources(self, comparison_name: str = "comparison.tsv",
                              manifest_name: str | None = None,
                              route_a_manifest_name: str | None = None) -> tuple[dict, dict]:
        inputs = {"dual_route_comparison_resource": "comparison"}
        resources = {"comparison": {"type": "dual_route_comparison", "path": comparison_name}}
        if manifest_name:
            inputs["dual_route_route_b_manifest_resource"] = "route_b_manifest"
            resources["route_b_manifest"] = {"type": "dual_route_route_b_manifest", "path": manifest_name}
        if route_a_manifest_name:
            inputs["dual_route_route_a_manifest_resource"] = "route_a_manifest"
            resources["route_a_manifest"] = {"type": "dual_route_route_a_manifest", "path": route_a_manifest_name}
        return inputs, resources

    def test_discordance_rows_are_appended_to_ledger(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            comparison_path = self._write(work, "comparison.tsv", FOUR_STATE_ROWS)
            inputs, resources = self._comparison_resources()
            out = self._run_case(work, root / "out", inputs=inputs, resources=resources)
            ledger = pd.read_csv(out / "locus_causal_scenarios.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(SCENARIO_COLUMNS, list(ledger.columns))
            self.assertEqual(11, len(ledger))
            self.assertTrue(ledger.assessment_state.eq("review_required").all())
            self.assertTrue(ledger.rank_contribution.eq("none").all())
            self.assertEqual({
                "coding_gene", "regulatory_element", "mapping_or_liftover_error",
                "reference_parent_mismatch", "sv_te_cnv_pav", "novel_or_unannotated_mechanism",
                "multi_gene_or_epistasis", "phenotype_context_review",
                "mapping_adequacy_review", "environment_specific",
            }, set(ledger.scenario_type))
            fingerprint = comparison_fingerprint(comparison_path)
            self.assertTrue(all(
                "comparison=cmp_" in value and f"fingerprint={fingerprint}" in value
                for value in ledger.source_ids
            ))
            regulatory = ledger[ledger.scenario_type.eq("regulatory_element")].iloc[0]
            self.assertIn("route_a_source_ids=SRC1", regulatory.source_ids)
            report = (out / "report.md").read_text(encoding="utf-8")
            self.assertIn("## Locus causal scenarios", report)
            self.assertIn("`review_required`", report)
            self.assertIn("discordance-derived: 11 (present)", report)
            manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("complete", manifest["steps"]["causal_scenarios"]["status"])
            self.assertEqual(11, manifest["steps"]["causal_scenarios"]["discordance_rows"])
            self.assertEqual("present", manifest["steps"]["causal_scenarios"]["discordance_input"])
            self.assertEqual("absent", manifest["steps"]["causal_scenarios"]["route_b_status"])
            self.assertEqual("present", manifest["scenario_ledger_status"]["discordance_input"])
            self.assertEqual(11, manifest["scenario_ledger_status"]["discordance_rows"])
            self.assertEqual("absent", manifest["scenario_ledger_status"]["route_b_status"])
            self.assertEqual(11, manifest["scenario_ledger_status"]["scenario_count"])
            state = json.loads((out / "case_state.json").read_text(encoding="utf-8"))
            self.assertEqual(11, state["scenario_count"])
            self.assertEqual(str(out / "locus_causal_scenarios.tsv"), state["scenario_ledger"])

    def test_discordance_input_never_changes_rank_or_stage_a(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root / "with")
            self._write(work, "comparison.tsv", FOUR_STATE_ROWS)
            inputs, resources = self._comparison_resources()
            out = self._run_case(work, root / "out", inputs=inputs, resources=resources)
            baseline = self._run_case(self._work(root / "baseline"), root / "baseline_out")
            for relative in ("rank/ranked_candidates.tsv", "rank/score_contributions.tsv",
                             "stage_a/deterministic_shortlist.tsv", "candidates.tsv"):
                baseline_path = baseline / relative
                if baseline_path.exists():
                    self.assertEqual(baseline_path.read_bytes(), (out / relative).read_bytes())

    def test_ledger_append_is_idempotent_across_resume(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            self._write(work, "comparison.tsv", FOUR_STATE_ROWS)
            inputs, resources = self._comparison_resources()
            out = self._run_case(work, root / "out", inputs=inputs, resources=resources)
            ledger_bytes = (out / "locus_causal_scenarios.tsv").read_bytes()
            ledger = pd.read_csv(out / "locus_causal_scenarios.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(11, len(ledger))
            self.assertFalse(ledger["scenario_id"].duplicated().any())
            state = json.loads((out / "case_state.json").read_text(encoding="utf-8"))
            first_hash = state["steps"]["causal_scenarios"]["input_hash"]
            out = self._run_case(work, root / "out", inputs=inputs, resources=resources, resume=True)
            self.assertEqual(ledger_bytes, (out / "locus_causal_scenarios.tsv").read_bytes())
            ledger2 = pd.read_csv(out / "locus_causal_scenarios.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(11, len(ledger2))
            self.assertFalse(ledger2["scenario_id"].duplicated().any())
            state = json.loads((out / "case_state.json").read_text(encoding="utf-8"))
            self.assertEqual(first_hash, state["steps"]["causal_scenarios"]["input_hash"])
            self.assertEqual("complete", state["steps"]["causal_scenarios"]["status"])
            self.assertEqual(11, state["steps"]["causal_scenarios"]["discordance_rows"])
            self.assertEqual(11, state["scenario_count"])

    def test_route_b_provider_missing_yields_not_evaluated_row(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            self._write(work, "comparison.tsv",
                        COMPARISON_HEADER + "\n"
                        "cmp_000\tha1\thb1\tconvergent\tSRC1\tSRC2\tSRC1\n")
            self._write(work, "interval_first_status.json", json.dumps({
                "schema_version": "1.4", "route": "interval_first", "case_id": "bundled_fas_like",
                "command_completed": True, "lane_variant": "b0_phenotype_free", "lane": "interval_informed",
                "phenotype_used": False, "blind_discovery_claim": False,
                "semantically_evaluable": False, "semantic_status": "not_evaluated_provider_missing",
                "semantic_evaluability_reason": "the command completed, but no interval gene has a usable "
                                                "Route B semantic provider annotation",
                "provider_coverage": {
                    "family_label_gene_count": 0, "orthogroup_gene_count": 0,
                    "domain_architecture_gene_count": 0, "broad_function_gene_count": 0,
                    "source_backed_function_gene_count": 0,
                },
            }, indent=2, sort_keys=True) + "\n")
            inputs, resources = self._comparison_resources(manifest_name="interval_first_status.json")
            out = self._run_case(work, root / "out", inputs=inputs, resources=resources)
            ledger = pd.read_csv(out / "locus_causal_scenarios.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(1, len(ledger))
            row = ledger.iloc[0]
            self.assertEqual(PROVIDER_MISSING_SCENARIO_TYPE, row.scenario_type)
            self.assertEqual("not_evaluated", row.assessment_state)
            self.assertEqual("", row.supporting_evidence)
            self.assertIn("route_b provider missing", row.reason)
            self.assertEqual(
                "provide route B provider or record route B as intentionally absent",
                row.recommended_next_experiment,
            )
            self.assertIn("comparison=cmp_000", row.source_ids)
            self.assertNotIn("phenotype_context_review", row.scenario_type)
            manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("not_evaluated", manifest["scenario_ledger_status"]["route_b_status"])
            self.assertEqual("not_evaluated", manifest["steps"]["causal_scenarios"]["route_b_status"])
            self.assertEqual(1, manifest["scenario_ledger_status"]["discordance_rows"])

    def test_route_a_provider_missing_never_becomes_b_only_or_neither(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            self._write(work, "comparison.tsv",
                        COMPARISON_HEADER + "\n"
                        "cmp_a_only\tha1\t\tphenotype_only_supported\tSRC1\t\t\n"
                        "cmp_neither\t\t\tno_support\t\t\t\n")
            self._write(work, "phenotype_first_status.json", json.dumps({
                "schema_version": "1.4", "route": "phenotype_first", "case_id": "bundled_fas_like",
                "hypothesis_count": 2, "verified_bait_count": 0,
            }, indent=2, sort_keys=True) + "\n")
            inputs, resources = self._comparison_resources(
                route_a_manifest_name="phenotype_first_status.json",
            )
            out = self._run_case(work, root / "out", inputs=inputs, resources=resources)
            ledger = pd.read_csv(out / "locus_causal_scenarios.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(2, len(ledger))
            self.assertTrue(ledger.scenario_type.eq(PROVIDER_MISSING_SCENARIO_TYPE).all())
            self.assertTrue(ledger.assessment_state.eq("not_evaluated").all())
            self.assertTrue(ledger.reason.str.contains("route_a provider missing").all())
            self.assertFalse(ledger.scenario_type.isin({"coding_gene", "phenotype_context_review"}).any())
            manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("not_evaluated", manifest["scenario_ledger_status"]["route_a_status"])
            self.assertEqual("not_evaluated", manifest["steps"]["causal_scenarios"]["route_a_status"])

    def test_provider_gap_rows_keep_distinct_ids_and_scenario_rows_survive(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            self._write(work, "comparison.tsv",
                        COMPARISON_HEADER + "\n"
                        "cmp_000\tha1\thb1\tconvergent\tSRC1\tSRC2\tSRC1\n"
                        "cmp_001\tha2\thb2\tconvergent\tSRC5\tSRC6\t\n")
            self._write(work, "interval_first_status.json", json.dumps({
                "schema_version": "1.4", "route": "interval_first", "case_id": "bundled_fas_like",
                "semantically_evaluable": False, "semantic_status": "not_evaluated_provider_missing",
                "provider_coverage": {"family_label_gene_count": 0, "orthogroup_gene_count": 0},
            }, indent=2, sort_keys=True) + "\n")
            self._write(work, "scenario_rows.tsv", SCENARIO_ROWS)
            inputs, resources = self._comparison_resources(manifest_name="interval_first_status.json")
            inputs["scenario_rows_resource"] = "scenario_rows"
            resources["scenario_rows"] = {"type": "causal_scenario_rows", "path": "scenario_rows.tsv"}
            out = self._run_case(work, root / "out", inputs=inputs, resources=resources)
            ledger = pd.read_csv(out / "locus_causal_scenarios.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(4, len(ledger))
            self.assertFalse(ledger["scenario_id"].duplicated().any())
            not_evaluated = ledger[ledger.assessment_state.eq("not_evaluated")]
            self.assertEqual(2, len(not_evaluated))
            self.assertEqual(2, not_evaluated["scenario_id"].nunique())
            self.assertTrue(not_evaluated.scenario_type.eq(PROVIDER_MISSING_SCENARIO_TYPE).all())
            self.assertEqual(
                {"comparison=cmp_000", "comparison=cmp_001"},
                {value.split(";", 1)[0] for value in not_evaluated.source_ids},
            )
            supported = ledger[ledger.assessment_state.eq("supported")]
            self.assertEqual(1, len(supported))
            self.assertEqual("coding_gene", supported.iloc[0].scenario_type)
            self.assertIn("candidate coding variant", supported.iloc[0].supporting_evidence)
            manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("not_evaluated", manifest["steps"]["causal_scenarios"]["route_b_status"])
            self.assertEqual(2, manifest["steps"]["causal_scenarios"]["discordance_rows"])
            self.assertEqual("present", manifest["steps"]["causal_scenarios"]["rows_input"])
            self.assertEqual("present", manifest["steps"]["causal_scenarios"]["discordance_input"])

    def test_invalid_comparison_degrades_without_losing_scenario_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            self._write(work, "comparison.tsv",
                        "family_id\tconvergence_status\nfam1\tphenotype_only_supported\n")
            self._write(work, "scenario_rows.tsv", SCENARIO_ROWS)
            inputs, resources = self._comparison_resources("comparison.tsv")
            inputs["scenario_rows_resource"] = "scenario_rows"
            resources["scenario_rows"] = {"type": "causal_scenario_rows", "path": "scenario_rows.tsv"}
            out = self._run_case(work, root / "out", inputs=inputs, resources=resources)
            ledger = pd.read_csv(out / "locus_causal_scenarios.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(2, len(ledger))
            self.assertEqual({"coding_gene", "regulatory_element"}, set(ledger.scenario_type))
            self.assertTrue(ledger.assessment_state.isin({"supported", "unassessed"}).all())
            manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            step = manifest["steps"]["causal_scenarios"]
            self.assertEqual("complete", step["status"])
            self.assertEqual("invalid", step["discordance_input"])
            self.assertEqual("present", step["rows_input"])
            self.assertIn("not a canonical comparison artifact", step["error"])
            self.assertEqual("invalid", manifest["scenario_ledger_status"]["discordance_input"])
            report = (out / "report.md").read_text(encoding="utf-8")
            self.assertIn("discordance input invalid", report)
            self.assertIn("not a canonical comparison artifact", report)
            self.assertIn("discordance rows not added to the ledger", report)

    def test_comparison_row_order_does_not_change_ledger(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            lines = FOUR_STATE_ROWS.rstrip("\n").splitlines()
            reversed_tsv = lines[0] + "\n" + "\n".join(lines[1:][::-1]) + "\n"
            self._write(work, "comparison_forward.tsv", FOUR_STATE_ROWS)
            self._write(work, "comparison_reversed.tsv", reversed_tsv)
            inputs_fwd, resources_fwd = self._comparison_resources("comparison_forward.tsv")
            inputs_rev, resources_rev = self._comparison_resources("comparison_reversed.tsv")
            out_fwd = self._run_case(work, root / "out_fwd", inputs=inputs_fwd, resources=resources_fwd)
            out_rev = self._run_case(work, root / "out_rev", inputs=inputs_rev, resources=resources_rev)
            self.assertEqual((out_fwd / "locus_causal_scenarios.tsv").read_bytes(),
                             (out_rev / "locus_causal_scenarios.tsv").read_bytes())

    def test_no_comparison_resource_leaves_ledger_untouched(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            self._write(work, "scenario_rows.tsv", SCENARIO_ROWS)
            out = self._run_case(
                work, root / "out",
                inputs={"scenario_rows_resource": "scenario_rows"},
                resources={"scenario_rows": {"type": "causal_scenario_rows", "path": "scenario_rows.tsv"}})
            ledger = pd.read_csv(out / "locus_causal_scenarios.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(2, len(ledger))
            self.assertEqual({"coding_gene", "regulatory_element"}, set(ledger.scenario_type))
            self.assertTrue(ledger.assessment_state.isin({"supported", "unassessed"}).all())
            self.assertTrue(ledger.source_ids.eq("unknown").all())
            self.assertTrue(ledger.reason.ne("").all())
            manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("absent", manifest["scenario_ledger_status"]["discordance_input"])
            self.assertEqual(0, manifest["scenario_ledger_status"]["discordance_rows"])
            self.assertEqual("absent", manifest["scenario_ledger_status"]["route_b_status"])
            self.assertEqual("complete", manifest["steps"]["causal_scenarios"]["status"])
            report = (out / "report.md").read_text(encoding="utf-8")
            self.assertIn("discordance-derived: 0 (absent)", report)


class ComparisonFingerprintTests(unittest.TestCase):
    """Fixed contract of the comparison input fingerprint helper."""

    def test_comparison_fingerprint_is_full_sha256_of_canonical_content(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "comparison.tsv"
            path.write_text(
                "comparison_id\troute_a_hypothesis_id\ncmp_002\tha2\ncmp_001\tha1\n",
                encoding="utf-8",
            )
            canonical = "comparison_id\troute_a_hypothesis_id\ncmp_001\tha1\ncmp_002\tha2"
            self.assertEqual(
                hashlib.sha256(canonical.encode("utf-8")).hexdigest(),
                comparison_fingerprint(path),
            )
            self.assertEqual(64, len(comparison_fingerprint(path)))

    def test_comparison_fingerprint_is_order_insensitive(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            forward = root / "forward.tsv"
            forward.write_text(
                "comparison_id\troute_a_hypothesis_id\ncmp_001\tha1\ncmp_002\tha2\ncmp_003\tha3\n",
                encoding="utf-8",
            )
            reversed_path = root / "reversed.tsv"
            reversed_path.write_text(
                "comparison_id\troute_a_hypothesis_id\ncmp_003\tha3\ncmp_002\tha2\ncmp_001\tha1\n",
                encoding="utf-8",
            )
            self.assertEqual(comparison_fingerprint(forward), comparison_fingerprint(reversed_path))

    def test_comparison_fingerprint_missing_path_raises(self) -> None:
        with self.assertRaises(FileNotFoundError):
            comparison_fingerprint(Path("does_not_exist_363.tsv"))


if __name__ == "__main__":
    unittest.main()
