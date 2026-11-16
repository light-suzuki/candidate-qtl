"""CLI wiring tests for the locus causal scenario ledger (issue #362)."""

from __future__ import annotations

import csv
import io
import json
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from candidate_qtl.causal_scenarios import SCENARIO_COLUMNS
from candidate_qtl.cli import EXIT_EXPECTED_FAILURE, EXIT_OK, main, run


class ScenarioCliTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def _run(self, arguments: list[str]) -> tuple[int, str, str]:
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            code = run(arguments)
        return code, out.getvalue(), err.getvalue()

    def _write_rows_tsv(self) -> Path:
        path = self.root / "rows.tsv"
        path.write_text(
            "scenario_type\tsupporting_evidence\tcontradicting_evidence\tassessment_state\t"
            "recommended_next_experiment\treason\n"
            "coding_gene\tcandidate coding variant\t\tsupported\tsequence candidate CDS\t"
            "coding hypothesis remains plausible\n"
            "regulatory_element\t\t\tunassessed\tpromoter review\tnoncoding space not evaluated\n",
            encoding="utf-8",
        )
        return path

    def _read_ledger(self, out_dir: Path) -> list[dict[str, str]]:
        with (out_dir / "locus_causal_scenarios.tsv").open(encoding="utf-8") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_scenario_create_writes_canonical_ledger(self) -> None:
        rows = self._write_rows_tsv()
        out_dir = self.root / "ledger"
        self.assertEqual(
            EXIT_OK,
            main(["scenario-create", "--input", str(rows), "--out-dir", str(out_dir), "--case", "case1"]),
        )
        ledger_path = out_dir / "locus_causal_scenarios.tsv"
        self.assertTrue(ledger_path.is_file())
        with ledger_path.open(encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            header = reader.fieldnames
            content = list(reader)
        self.assertEqual(SCENARIO_COLUMNS, header)
        self.assertEqual(2, len(content))
        self.assertTrue(all(row["rank_contribution"] == "none" for row in content))

    def test_case_flag_overrides_row_case_ids(self) -> None:
        rows = self.root / "rows.tsv"
        rows.write_text(
            "case_id\tscenario_type\tsupporting_evidence\tassessment_state\n"
            "unknown\tcoding_gene\tE1\tsupported\n",
            encoding="utf-8",
        )
        out_dir = self.root / "ledger"
        self.assertEqual(
            EXIT_OK,
            main(["scenario-create", "--input", str(rows), "--out-dir", str(out_dir), "--case", "forced_case"]),
        )
        rows_out = self._read_ledger(out_dir)
        self.assertEqual(1, len(rows_out))
        self.assertEqual("forced_case", rows_out[0]["case_id"])

    def test_scenario_report_prints_case_and_scenario_ids(self) -> None:
        rows = self._write_rows_tsv()
        out_dir = self.root / "ledger"
        self.assertEqual(
            EXIT_OK,
            main(["scenario-create", "--input", str(rows), "--out-dir", str(out_dir), "--case", "case1"]),
        )
        buffer = io.StringIO()
        with redirect_stdout(buffer):
            code = main(["scenario-report", "--ledger", str(out_dir / "locus_causal_scenarios.tsv")])
        self.assertEqual(EXIT_OK, code)
        output = buffer.getvalue()
        self.assertIn("## Locus causal scenarios", output)
        self.assertIn("case1", output)
        self.assertIn("SCN-", output)
        self.assertIn("`coding_gene` / `supported`", output)
        self.assertIn("`regulatory_element` / `unassessed`", output)

    def test_scenario_create_json_reports_machine_readable_breakdown(self) -> None:
        rows = self._write_rows_tsv()
        buffer = io.StringIO()
        with redirect_stdout(buffer):
            code = main(
                ["scenario-create", "--input", str(rows), "--out-dir", str(self.root / "ledger"),
                 "--case", "case1", "--json"]
            )
        self.assertEqual(EXIT_OK, code)
        result = json.loads(buffer.getvalue())
        self.assertEqual(["case1"], result["case_ids"])
        self.assertEqual(2, result["scenario_count"])
        self.assertEqual(1, result["assessment_breakdown"]["supported"])
        self.assertEqual(1, result["assessment_breakdown"]["unassessed"])
        self.assertEqual(0, result["assessment_breakdown"]["contradicted"])
        self.assertEqual(0, result["assessment_breakdown"]["mixed"])
        self.assertEqual(0, result["assessment_breakdown"]["review_required"])
        self.assertTrue(result["out_path"].endswith("locus_causal_scenarios.tsv"))

    def test_scenario_report_json_returns_summary_records(self) -> None:
        rows = self._write_rows_tsv()
        out_dir = self.root / "ledger"
        main(["scenario-create", "--input", str(rows), "--out-dir", str(out_dir), "--case", "case1"])
        buffer = io.StringIO()
        with redirect_stdout(buffer):
            code = main(
                ["scenario-report", "--ledger", str(out_dir / "locus_causal_scenarios.tsv"), "--json"]
            )
        self.assertEqual(EXIT_OK, code)
        result = json.loads(buffer.getvalue())
        self.assertEqual(2, result["scenario_count"])
        by_type = {record["scenario_type"]: record for record in result["summary"]}
        self.assertEqual(1, by_type["coding_gene"]["supported_rows"])
        self.assertEqual(1, by_type["regulatory_element"]["unassessed_rows"])

    def test_yaml_list_input_builds_ledger(self) -> None:
        rows = self.root / "rows.yaml"
        rows.write_text(
            "- scenario_type: sv_te_cnv_pav\n"
            "  assessment_state: unassessed\n"
            "  recommended_next_experiment: inspect long-read/SV evidence\n"
            "  reason: SV/PAV space has not been evaluated\n",
            encoding="utf-8",
        )
        out_dir = self.root / "ledger"
        self.assertEqual(
            EXIT_OK,
            main(["scenario-create", "--input", str(rows), "--out-dir", str(out_dir), "--case", "case1"]),
        )
        rows_out = self._read_ledger(out_dir)
        self.assertEqual(1, len(rows_out))
        self.assertEqual("sv_te_cnv_pav", rows_out[0]["scenario_type"])
        self.assertEqual("unassessed", rows_out[0]["assessment_state"])

    def test_invalid_scenario_type_exits_one_without_traceback(self) -> None:
        rows = self.root / "bad.tsv"
        rows.write_text("scenario_type\tassessment_state\nunknown_type\tunassessed\n", encoding="utf-8")
        code, _, err = self._run(
            ["scenario-create", "--input", str(rows), "--out-dir", str(self.root / "ledger"), "--case", "case1"]
        )
        self.assertEqual(EXIT_EXPECTED_FAILURE, code)
        self.assertIn("unsupported scenario_type", err)
        self.assertNotIn("Traceback", err)

    def test_invalid_yaml_exits_one_without_traceback(self) -> None:
        rows = self.root / "bad.yaml"
        rows.write_text("{{{{", encoding="utf-8")
        code, _, err = self._run(
            ["scenario-create", "--input", str(rows), "--out-dir", str(self.root / "ledger"), "--case", "case1"]
        )
        self.assertEqual(EXIT_EXPECTED_FAILURE, code)
        self.assertIn("could not be parsed", err)
        self.assertNotIn("Traceback", err)

    def test_missing_input_is_an_expected_error(self) -> None:
        code, _, err = self._run(
            ["scenario-create", "--input", str(self.root / "absent.tsv"),
             "--out-dir", str(self.root / "ledger"), "--case", "case1"]
        )
        self.assertEqual(EXIT_EXPECTED_FAILURE, code)
        self.assertIn("scenario input does not exist", err)
        self.assertNotIn("Traceback", err)

    def test_empty_ledger_report_is_unassessed(self) -> None:
        ledger = self.root / "locus_causal_scenarios.tsv"
        ledger.write_text("\t".join(SCENARIO_COLUMNS) + "\n", encoding="utf-8")
        buffer = io.StringIO()
        with redirect_stdout(buffer):
            code = main(["scenario-report", "--ledger", str(ledger)])
        self.assertEqual(EXIT_OK, code)
        self.assertIn("No causal scenarios have been registered; status is unassessed.", buffer.getvalue())


if __name__ == "__main__":
    unittest.main()
