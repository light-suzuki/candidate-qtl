"""PR #372 (issue #362): case-run wires the locus causal scenario ledger.

The new ``causal_scenarios`` step runs after next-experiment and before the
report. It never touches the rank pipeline (every row keeps
``rank_contribution=none``); old cases without scenario inputs produce an
honest header-only unassessed ledger, and the optional #360
qtl_effect_contract input degrades to a recorded status instead of crashing.
"""

from __future__ import annotations

import json
import shutil
import tempfile
import unittest
from pathlib import Path

import pandas as pd
import yaml

from candidate_qtl.causal_scenarios import (
    SCENARIO_COLUMNS,
    build_scenario_ledger,
    load_qtl_effect_contract,
    load_scenario_rows,
    sign_concordance_evidence,
)
from candidate_qtl.cli import main


class CaseRunScenarioLedgerWiringTests(unittest.TestCase):
    """End-to-end case-run wiring, following the case_run_elimination harness."""

    SCENARIO_ROWS = (
        "case_id\tscenario_type\tsupporting_evidence\tcontradicting_evidence\tassessment_state\t"
        "recommended_next_experiment\treason\n"
        "bundled_fas_like\tcoding_gene\tcandidate coding variant\t\tsupported\t"
        "sequence candidate CDS\tcoding hypothesis remains plausible\n"
        "bundled_fas_like\tregulatory_element\t\t\tunassessed\tpromoter review\t"
        "noncoding space not evaluated\n"
    )
    CONTRACT_INCOMPATIBLE = (
        "case_id\tsign_concordance\tfavorable_parent\teffect_direction\n"
        "bundled_fas_like\tincompatible\tReferenceParent\tincrease\n"
    )
    CONTRACT_CONCORDANT = (
        "case_id\tsign_concordance\tfavorable_parent\teffect_direction\n"
        "bundled_fas_like\tconcordant\tReferenceParent\tincrease\n"
    )

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

    def _scenario_resources(self) -> tuple[dict, dict]:
        return (
            {"scenario_rows_resource": "scenario_rows"},
            {"scenario_rows": {"type": "causal_scenario_rows", "path": "scenario_rows.tsv"}},
        )

    def test_scenario_rows_resource_produces_ledger_and_report_section(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            self._write(work, "scenario_rows.tsv", self.SCENARIO_ROWS)
            inputs, resources = self._scenario_resources()
            out = self._run_case(work, root / "out", inputs=inputs, resources=resources)
            ledger = pd.read_csv(out / "locus_causal_scenarios.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(SCENARIO_COLUMNS, list(ledger.columns))
            self.assertEqual(2, len(ledger))
            self.assertEqual({"coding_gene", "regulatory_element"}, set(ledger.scenario_type))
            self.assertTrue(ledger.rank_contribution.eq("none").all())
            report = (out / "report.md").read_text(encoding="utf-8")
            self.assertIn("## Locus causal scenarios", report)
            self.assertIn("bundled_fas_like/", report)
            manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            step = manifest["steps"]["causal_scenarios"]
            self.assertEqual("complete", step["status"])
            self.assertEqual("present", step["rows_input"])
            self.assertEqual(2, step["scenario_count"])
            self.assertEqual("absent", step["sign_contract_status"])
            self.assertEqual("present", manifest["scenario_ledger_status"]["rows_input"])
            self.assertEqual("absent", manifest["scenario_ledger_status"]["sign_contract_status"])
            self.assertEqual(2, manifest["scenario_ledger_status"]["scenario_count"])
            state = json.loads((out / "case_state.json").read_text(encoding="utf-8"))
            self.assertIn("causal_scenarios", state["steps"])
            self.assertEqual(str(out / "locus_causal_scenarios.tsv"), state["scenario_ledger"])
            self.assertEqual(2, state["scenario_count"])

    def test_scenario_inputs_never_change_the_deterministic_rank(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            scen_work = self._work(root / "scen")
            self._write(scen_work, "scenario_rows.tsv", self.SCENARIO_ROWS)
            inputs, resources = self._scenario_resources()
            out = self._run_case(scen_work, root / "out", inputs=inputs, resources=resources)
            baseline = self._run_case(self._work(root / "baseline"), root / "baseline_out")
            self.assertEqual((baseline / "rank" / "ranked_candidates.tsv").read_bytes(),
                             (out / "rank" / "ranked_candidates.tsv").read_bytes())
            self.assertEqual((baseline / "rank" / "score_contributions.tsv").read_bytes(),
                             (out / "rank" / "score_contributions.tsv").read_bytes())
            if (baseline / "stage_a" / "deterministic_shortlist.tsv").exists():
                self.assertEqual((baseline / "stage_a" / "deterministic_shortlist.tsv").read_bytes(),
                                 (out / "stage_a" / "deterministic_shortlist.tsv").read_bytes())
            self.assertEqual((baseline / "candidates.tsv").read_bytes(),
                             (out / "candidates.tsv").read_bytes())

    def test_invalid_scenario_rows_degrades_without_crashing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            self._write(work, "scenario_rows.tsv", "scenario_type\tassessment_state\nunknown_type\tunassessed\n")
            inputs, resources = self._scenario_resources()
            out = self._run_case(work, root / "out", inputs=inputs, resources=resources)
            ledger = pd.read_csv(out / "locus_causal_scenarios.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(0, len(ledger))
            self.assertEqual(SCENARIO_COLUMNS, list(ledger.columns))
            manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("invalid", manifest["scenario_ledger_status"]["rows_input"])
            self.assertEqual("invalid", manifest["steps"]["causal_scenarios"]["rows_input"])
            self.assertEqual(0, manifest["scenario_ledger_status"]["scenario_count"])
            report = (out / "report.md").read_text(encoding="utf-8")
            self.assertIn("scenario_rows input invalid", report)
            self.assertIn("unsupported scenario_type", report)
            self.assertIn("ledger left unassessed", report)

    def test_old_case_without_scenario_inputs_stays_honestly_unassessed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            out = self._run_case(work, root / "out")
            ledger = pd.read_csv(out / "locus_causal_scenarios.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(0, len(ledger))
            self.assertEqual(SCENARIO_COLUMNS, list(ledger.columns))
            report = (out / "report.md").read_text(encoding="utf-8")
            self.assertIn("No causal scenarios have been registered", report)
            self.assertIn("unassessed", report)
            state = json.loads((out / "case_state.json").read_text(encoding="utf-8"))
            self.assertEqual(0, state["scenario_count"])
            self.assertEqual(str(out / "locus_causal_scenarios.tsv"), state["scenario_ledger"])
            manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("complete", manifest["steps"]["causal_scenarios"]["status"])
            self.assertEqual("absent", manifest["steps"]["causal_scenarios"]["rows_input"])
            self.assertEqual("absent", manifest["scenario_ledger_status"]["sign_contract_status"])
            self.assertEqual(0, manifest["scenario_ledger_status"]["scenario_count"])

    def test_resume_rerun_keeps_ledger_and_step_hash_stable(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            self._write(work, "scenario_rows.tsv", self.SCENARIO_ROWS)
            inputs, resources = self._scenario_resources()
            out = self._run_case(work, root / "out", inputs=inputs, resources=resources)
            ledger_bytes = (out / "locus_causal_scenarios.tsv").read_bytes()
            state = json.loads((out / "case_state.json").read_text(encoding="utf-8"))
            first_hash = state["steps"]["causal_scenarios"]["input_hash"]
            out = self._run_case(work, root / "out", inputs=inputs, resources=resources, resume=True)
            self.assertEqual(ledger_bytes, (out / "locus_causal_scenarios.tsv").read_bytes())
            ledger = pd.read_csv(out / "locus_causal_scenarios.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(2, len(ledger))
            self.assertFalse(ledger["scenario_id"].duplicated().any())
            state = json.loads((out / "case_state.json").read_text(encoding="utf-8"))
            self.assertEqual(first_hash, state["steps"]["causal_scenarios"]["input_hash"])
            self.assertEqual("complete", state["steps"]["causal_scenarios"]["status"])
            self.assertEqual(2, state["scenario_count"])

    def test_written_ledger_roundtrips_through_default_pandas_read(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            self._write(work, "scenario_rows.tsv", self.SCENARIO_ROWS)
            inputs, resources = self._scenario_resources()
            out = self._run_case(work, root / "out", inputs=inputs, resources=resources)
            written = pd.read_csv(out / "locus_causal_scenarios.tsv", sep="\t", dtype=str).fillna("")
            reloaded = pd.read_csv(out / "locus_causal_scenarios.tsv", sep="\t")
            rebuilt = build_scenario_ledger(reloaded, case_id="bundled_fas_like")
            self.assertEqual(written["scenario_id"].tolist(), rebuilt["scenario_id"].tolist())
            regulatory = rebuilt[rebuilt.scenario_type.eq("regulatory_element")].iloc[0]
            self.assertEqual("unassessed", regulatory.assessment_state)
            self.assertTrue(pd.isna(
                reloaded[reloaded.scenario_type.eq("regulatory_element")].iloc[0].supporting_evidence))

    def test_sign_contract_incompatible_appends_contradicting_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            self._write(work, "scenario_rows.tsv", self.SCENARIO_ROWS)
            self._write(work, "contract.tsv", self.CONTRACT_INCOMPATIBLE)
            inputs, resources = self._scenario_resources()
            inputs["qtl_effect_contract_resource"] = "sign_contract"
            resources["sign_contract"] = {"type": "qtl_effect_contract", "path": "contract.tsv"}
            out = self._run_case(work, root / "out", inputs=inputs, resources=resources)
            ledger = pd.read_csv(out / "locus_causal_scenarios.tsv", sep="\t", dtype=str).fillna("")
            coding = ledger[ledger.scenario_type.eq("coding_gene")].iloc[0]
            self.assertEqual("sign discordance: ReferenceParent/increase",
                             coding.contradicting_evidence)
            self.assertEqual("qtl_effect_contract:bundled_fas_like", coding.source_ids)
            regulatory = ledger[ledger.scenario_type.eq("regulatory_element")].iloc[0]
            self.assertEqual("sign discordance: ReferenceParent/increase",
                             regulatory.contradicting_evidence)
            manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("present", manifest["scenario_ledger_status"]["sign_contract_status"])
            self.assertEqual("present", manifest["steps"]["causal_scenarios"]["sign_contract_status"])

    def test_sign_contract_concordant_appends_supporting_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            self._write(work, "scenario_rows.tsv", self.SCENARIO_ROWS)
            self._write(work, "contract.tsv", self.CONTRACT_CONCORDANT)
            inputs, resources = self._scenario_resources()
            inputs["qtl_effect_contract_resource"] = "sign_contract"
            resources["sign_contract"] = {"type": "qtl_effect_contract", "path": "contract.tsv"}
            out = self._run_case(work, root / "out", inputs=inputs, resources=resources)
            ledger = pd.read_csv(out / "locus_causal_scenarios.tsv", sep="\t", dtype=str).fillna("")
            coding = ledger[ledger.scenario_type.eq("coding_gene")].iloc[0]
            self.assertIn("sign concordance: ReferenceParent/increase", coding.supporting_evidence)
            self.assertEqual("candidate coding variant; sign concordance: ReferenceParent/increase",
                             coding.supporting_evidence)
            self.assertIn("qtl_effect_contract:", coding.source_ids)

    def test_missing_or_broken_declared_contract_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # (a) declared resource whose file does not exist
            work = self._work(root / "missing")
            self._write(work, "scenario_rows.tsv", self.SCENARIO_ROWS)
            inputs, resources = self._scenario_resources()
            inputs["qtl_effect_contract_resource"] = "sign_contract"
            resources["sign_contract"] = {"type": "qtl_effect_contract", "path": "missing_contract.tsv"}
            with self.assertRaises(ValueError):
                self._run_case(work, root / "out", inputs=inputs, resources=resources)
            # (b) unparsable YAML
            work2 = self._work(root / "broken")
            self._write(work2, "scenario_rows.tsv", self.SCENARIO_ROWS)
            self._write(work2, "contract.yaml", "{{{{ not valid yaml")
            inputs2, resources2 = self._scenario_resources()
            inputs2["qtl_effect_contract_resource"] = "sign_contract"
            resources2["sign_contract"] = {"type": "qtl_effect_contract", "path": "contract.yaml"}
            with self.assertRaises((ValueError, yaml.YAMLError)):
                self._run_case(work2, root / "out2", inputs=inputs2, resources=resources2)

    def test_without_contract_input_nothing_changes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            self._write(work, "scenario_rows.tsv", self.SCENARIO_ROWS)
            inputs, resources = self._scenario_resources()
            out = self._run_case(work, root / "out", inputs=inputs, resources=resources)
            ledger = pd.read_csv(out / "locus_causal_scenarios.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(ledger.contradicting_evidence.eq("").all())
            self.assertTrue(ledger.source_ids.eq("unknown").all())
            report = (out / "report.md").read_text(encoding="utf-8")
            self.assertNotIn("sign discordance", report)
            self.assertNotIn("qtl_effect_contract:", report)
            manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("absent", manifest["scenario_ledger_status"]["sign_contract_status"])

    def test_yaml_scenario_rows_and_keyed_yaml_contract_flow_through_case_run(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            self._write(work, "rows.yaml",
                        "- scenario_type: coding_gene\n"
                        "  supporting_evidence: candidate coding variant\n"
                        "  assessment_state: supported\n"
                        "  recommended_next_experiment: sequence candidate CDS\n"
                        "  reason: coding hypothesis remains plausible\n"
                        "- scenario_type: environment_specific\n")
            self._write(work, "contract.yaml",
                        "qtl_effect_contract:\n"
                        "  case_id: bundled_fas_like\n"
                        "  sign_concordance: concordant\n"
                        "  favorable_parent: ReferenceParent\n"
                        "  effect_direction: increase\n")
            out = self._run_case(
                work, root / "out",
                inputs={"scenario_rows_resource": "rows", "qtl_effect_contract_resource": "contract"},
                resources={"rows": {"type": "causal_scenario_rows", "path": "rows.yaml"},
                           "contract": {"type": "qtl_effect_contract", "path": "contract.yaml"}})
            ledger = pd.read_csv(out / "locus_causal_scenarios.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(2, len(ledger))
            coding = ledger[ledger.scenario_type.eq("coding_gene")].iloc[0]
            self.assertIn("sign concordance: ReferenceParent/increase", coding.supporting_evidence)
            self.assertIn("qtl_effect_contract:bundled_fas_like", coding.source_ids)
            env = ledger[ledger.scenario_type.eq("environment_specific")].iloc[0]
            self.assertEqual("unassessed", env.assessment_state)
            manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("present", manifest["scenario_ledger_status"]["rows_input"])
            self.assertEqual("present", manifest["scenario_ledger_status"]["sign_contract_status"])


class ScenarioWiringHelperTests(unittest.TestCase):
    """Unit contracts for the case-run wiring helpers."""

    def test_load_scenario_rows_tsv_preserves_empty_cells(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "rows.tsv"
            path.write_text("scenario_type\tsupporting_evidence\ncoding_gene\tE1\nregulatory_element\t\n",
                            encoding="utf-8")
            rows = load_scenario_rows(path)
        self.assertEqual(2, len(rows))
        self.assertEqual("coding_gene", rows[0]["scenario_type"])
        self.assertEqual("E1", rows[0]["supporting_evidence"])
        self.assertEqual("regulatory_element", rows[1]["scenario_type"])
        self.assertEqual("", rows[1]["supporting_evidence"])

    def test_load_scenario_rows_yaml_list_and_rejections(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "rows.yaml"
            path.write_text("- scenario_type: coding_gene\n  supporting_evidence: E1\n"
                            "- scenario_type: environment_specific\n", encoding="utf-8")
            rows = load_scenario_rows(path)
            self.assertEqual(2, len(rows))
            self.assertEqual("E1", rows[0]["supporting_evidence"])
            empty = root / "empty.yaml"
            empty.write_text("", encoding="utf-8")
            self.assertEqual([], load_scenario_rows(empty))
            bad = root / "bad.yaml"
            bad.write_text("{{{{", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "could not be parsed"):
                load_scenario_rows(bad)
            not_list = root / "not_list.yaml"
            not_list.write_text("scenario_type: coding_gene\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "list"):
                load_scenario_rows(not_list)
            self.assertRaises(FileNotFoundError, load_scenario_rows, root / "absent.tsv")

    def test_load_qtl_effect_contract_yaml_forms_and_tsv(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            keyed = root / "contract.yaml"
            keyed.write_text("qtl_effect_contract:\n- case_id: c1\n  sign_concordance: incompatible\n",
                             encoding="utf-8")
            self.assertEqual("incompatible", load_qtl_effect_contract(keyed)[0]["sign_concordance"])
            top = root / "top.yaml"
            top.write_text("- case_id: c1\n  sign_concordance: concordant\n", encoding="utf-8")
            self.assertEqual("concordant", load_qtl_effect_contract(top)[0]["sign_concordance"])
            scalar = root / "scalar.yaml"
            scalar.write_text("hello\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "list"):
                load_qtl_effect_contract(scalar)
            mapping = root / "mapping.yaml"
            mapping.write_text("case_id: c1\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "qtl_effect_contract"):
                load_qtl_effect_contract(mapping)
            garbage = root / "garbage.yaml"
            garbage.write_text("{{{{", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "could not be parsed"):
                load_qtl_effect_contract(garbage)
            tsv = root / "contract.tsv"
            tsv.write_text("case_id\tsign_concordance\nc1\tunknown\n", encoding="utf-8")
            self.assertEqual("unknown", load_qtl_effect_contract(tsv)[0]["sign_concordance"])

    def test_sign_concordance_evidence_binds_only_incompatible_or_concordant(self) -> None:
        rows = [
            {"case_id": "c1", "sign_concordance": "incompatible",
             "favorable_parent": "ParentA", "effect_direction": "increase"},
            {"case_id": "c2", "sign_concordance": "concordant",
             "favorable_parent": "ParentB", "effect_direction": "decrease"},
            {"case_id": "c3", "sign_concordance": "unknown"},
            {"case_id": "c4", "sign_concordance": "incompatible"},
        ]
        self.assertEqual(("incompatible", "sign discordance: ParentA/increase"),
                         sign_concordance_evidence(rows, "c1"))
        self.assertEqual(("concordant", "sign concordance: ParentB/decrease"),
                         sign_concordance_evidence(rows, "c2"))
        self.assertEqual((None, None), sign_concordance_evidence(rows, "c3"))
        self.assertEqual(("incompatible", "sign discordance"), sign_concordance_evidence(rows, "c4"))
        self.assertEqual((None, None), sign_concordance_evidence(rows, "missing"))
        self.assertEqual((None, None), sign_concordance_evidence([], "c1"))
        self.assertEqual(("incompatible", "sign discordance: ParentA/increase"),
                         sign_concordance_evidence(pd.DataFrame(rows), "c1"))


if __name__ == "__main__":
    unittest.main()
