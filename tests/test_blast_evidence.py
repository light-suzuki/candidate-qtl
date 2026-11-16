"""Unit tests for the BLAST+ execution layer (issue #374 / #355).

blastp / makeblastdb are mocked (unittest.mock); no real BLAST+ is required.
The freeze itself is produced by the real ``freeze_blast_panel`` so the
verify-then-search ordering is exercised against genuine manifests. Covers:
the not_evaluated_* abstention ladder (no verified bait / BLAST+ missing /
makeblastdb failure / blastp failure), fail-closed freeze verification,
threshold filtering and ambiguity classification, short-query mode arguments,
evidence row shape, and search determinism.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from candidate_qtl.blast_evidence import (
    AMBIGUITY_VOCAB,
    BLAST_EVIDENCE_COLUMNS,
    ELIGIBLE_INTERVAL_HIT_COLUMNS,
    PROVIDER_STATUS_VOCAB,
    _classify_ambiguity,
    _evaluate_bait,
    _parse_blast_hits,
    run_panel_blastp,
)
from candidate_qtl.blast_panel import DEFAULT_THRESHOLDS, freeze_blast_panel
from candidate_qtl.core import sha256

FIXTURE = Path(__file__).resolve().parent / "fixtures" / "route_a"

S = "CDLYQWLICSGRWNGGYHCLKFRVGRYVDSYIMRQYFNMFLEVAWCGNRTDNRMNHMNNQDKTTEQPNHVNKQKNLGSRNSEYNGGSQYKGNLRRCKLTIYEEPFSWSDAPQCWSCIFET"

# Fixed, deterministic blastp output: bait1 hits g1/g2/g4 (g4 within threshold),
# bait2 ties g6/g7, bait3 has a single low-coverage hit that must be filtered.
RAW_HITS = (
    "bait1\tg1\t100.000\t120\t120\t120\t1\t120\t1\t120\t0.0\t250.0\n"
    "bait1\tg2\t80.000\t120\t120\t120\t1\t120\t1\t120\t1e-40\t200.0\n"
    "bait1\tg4\t70.000\t120\t120\t120\t1\t120\t1\t120\t1e-30\t170.0\n"
    "bait2\tg6\t100.000\t116\t116\t116\t1\t116\t1\t116\t0.0\t240.0\n"
    "bait2\tg7\t100.000\t116\t116\t116\t1\t116\t1\t116\t0.0\t240.0\n"
    "bait3\tg1\t100.000\t48\t48\t120\t1\t48\t1\t48\t1e-10\t90.0\n"
)


def _which_side_effect(name: str) -> str | None:
    return {"blastp": "/usr/bin/blastp", "makeblastdb": "/usr/bin/makeblastdb"}.get(name)


class MockBlastTools:
    """Deterministic mock blastp/makeblastdb: no real BLAST+ is ever invoked.

    makeblastdb materializes the ``.pdb`` file the real tool would create;
    blastp writes its ``-out`` file with the fixed hit table.
    """

    def __init__(
        self,
        raw_hits: str = RAW_HITS,
        makeblastdb_returncode: int = 0,
        blastp_returncode: int = 0,
    ) -> None:
        self.raw_hits = raw_hits
        self.makeblastdb_returncode = makeblastdb_returncode
        self.blastp_returncode = blastp_returncode
        self.blastp_args: list[str] = []

    def __enter__(self) -> "MockBlastTools":
        self._which = mock.patch("candidate_qtl.blast_evidence.shutil.which", side_effect=_which_side_effect)
        self._run_tool = mock.patch("candidate_qtl.blast_evidence.run_tool", side_effect=self._fake_run_tool)
        self._which.start()
        self._run_tool.start()
        return self

    def __exit__(self, *exc_info: object) -> None:
        self._run_tool.stop()
        self._which.stop()

    def _fake_run_tool(self, command: list[str], **kwargs: object) -> subprocess.CompletedProcess:
        args = list(command)
        if args and args[0] == "wsl":
            args = args[1:]
        tool = args[0]
        if tool == "makeblastdb":
            db_out = Path(args[args.index("-out") + 1])
            db_out.parent.mkdir(parents=True, exist_ok=True)
            Path(str(db_out) + ".pdb").write_bytes(b"mock-db")
            return subprocess.CompletedProcess(command, self.makeblastdb_returncode, "", "makeblastdb boom")
        if tool == "blastp":
            self.blastp_args = args
            out = Path(args[args.index("-out") + 1])
            out.parent.mkdir(parents=True, exist_ok=True)
            query = Path(args[args.index("-query") + 1])
            query_ids = {line[1:].strip() for line in query.read_text(encoding="utf-8").splitlines() if line.startswith(">")}
            raw = "".join(line for line in self.raw_hits.splitlines(keepends=True) if line.split("\t", 1)[0] in query_ids)
            out.write_text(raw, encoding="utf-8")
            return subprocess.CompletedProcess(command, self.blastp_returncode, "", "blastp boom")
        raise AssertionError(f"unexpected external tool: {command}")


class ModeAwareBlastTools(MockBlastTools):
    """Synthetic provider whose output depends only on the query FASTA."""

    def __init__(self, fail_mode: str | None = None) -> None:
        super().__init__(raw_hits="")
        self.fail_mode = fail_mode
        self.calls: list[tuple[str, list[str]]] = []

    def _fake_run_tool(self, command: list[str], **kwargs: object) -> subprocess.CompletedProcess:
        args = list(command)
        if args and args[0] == "wsl":
            args = args[1:]
        if args[0] == "makeblastdb":
            db_out = Path(args[args.index("-out") + 1])
            db_out.parent.mkdir(parents=True, exist_ok=True)
            Path(str(db_out) + ".pdb").write_bytes(b"mock-db")
            return subprocess.CompletedProcess(command, 0, "", "")
        query_path = Path(args[args.index("-query") + 1])
        mode = "short" if "-matrix" in args else "standard"
        self.calls.append((mode, args))
        out = Path(args[args.index("-out") + 1])
        if mode == self.fail_mode:
            return subprocess.CompletedProcess(command, 2, "", f"{mode} failed")
        rows = []
        for query_id, sequence in zip(query_path.read_text(encoding="utf-8").splitlines()[::2], query_path.read_text(encoding="utf-8").splitlines()[1::2]):
            qid = query_id.removeprefix(">").strip()
            qlen = len(sequence.strip())
            rows.append(f"{qid}\tg1\t100.000\t{qlen}\t{qlen}\t120\t1\t{qlen}\t1\t{qlen}\t0.0\t250.0")
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("\n".join(rows) + ("\n" if rows else ""), encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")


class MissingRawBlastTools(MockBlastTools):
    def _fake_run_tool(self, command: list[str], **kwargs: object) -> subprocess.CompletedProcess:
        args = list(command)
        if args and args[0] == "wsl":
            args = args[1:]
        if args[0] == "makeblastdb":
            db_out = Path(args[args.index("-out") + 1])
            db_out.parent.mkdir(parents=True, exist_ok=True)
            Path(str(db_out) + ".pdb").write_bytes(b"mock-db")
            return subprocess.CompletedProcess(command, 0, "", "")
        return subprocess.CompletedProcess(command, 0, "", "")


def _hit(**overrides: object) -> dict[str, object]:
    hit: dict[str, object] = {
        "qseqid": "bait1", "sseqid": "g1", "pident": 100.0, "length": 120,
        "qlen": 120, "slen": 120, "qstart": 1, "qend": 120, "sstart": 1,
        "send": 120, "evalue": 0.0, "bitscore": 250.0,
    }
    hit.update(overrides)
    return hit


def _freeze(root: Path, panel: Path | None = None) -> dict[str, object]:
    return freeze_blast_panel(str(panel or (FIXTURE / "panel_ok.tsv")), str(root / "freeze"), case_id="unit")


def _all_review_required_panel(root: Path) -> Path:
    """The fixture panel with every accession made invalid so all rows fail V-ACC."""
    text = (FIXTURE / "panel_ok.tsv").read_text(encoding="utf-8")
    for accession in ("Q9TES1", "P9TES2", "Q9TES3", "Q9TES4"):
        text = text.replace(accession, "P0FAKE")
    path = root / "all_review_required.tsv"
    path.write_text(text, encoding="utf-8")
    return path


def _panel_with_ids(root: Path, panel_ids: set[str]) -> Path:
    lines = (FIXTURE / "panel_ok.tsv").read_text(encoding="utf-8").splitlines()
    path = root / ("panel_" + "_".join(sorted(panel_ids)) + ".tsv")
    path.write_text("\n".join([lines[0]] + [line for line in lines[1:] if line.split("\t", 1)[0] in panel_ids]) + "\n", encoding="utf-8")
    return path


class NoVerifiedBaitTests(unittest.TestCase):
    def test_no_verified_bait_is_not_evaluated_with_header_only_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            freeze = _freeze(root, _all_review_required_panel(root))
            self.assertEqual(0, freeze["verified_count"])

            result = run_panel_blastp(str(freeze["manifest"]), str(FIXTURE / "interval_proteins.fa"), str(root / "search"))

            self.assertEqual("not_evaluated_no_verified_bait", result["provider_status"])
            self.assertEqual(0, result["rows_evaluated"])
            self.assertEqual(0, result["exit_code"])
            self.assertEqual("post_freeze_similarity_only", result["interval_usage"])

            evidence = pd.read_csv(result["evidence"], sep="\t", dtype=str)
            self.assertEqual(0, len(evidence))
            self.assertEqual(BLAST_EVIDENCE_COLUMNS, list(evidence.columns))
            eligible = pd.read_csv(result["eligible_interval_hits"], sep="\t", dtype=str)
            self.assertEqual(ELIGIBLE_INTERVAL_HIT_COLUMNS, list(eligible.columns))
            self.assertEqual(0, len(eligible))

            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual("not_evaluated_no_verified_bait", manifest["provider_status"])
            self.assertEqual(["bait1", "bait2", "bait3", "bait4"], manifest["skipped_review_required"])
            self.assertTrue(manifest["reason"])
            self.assertEqual("post_freeze_similarity_only", manifest["interval_usage"])
            self.assertNotIn("interval_accessed", manifest)
            self.assertEqual(0, manifest["counts"]["rows_evaluated"])


class NotEvaluatedStickyTests(unittest.TestCase):
    def test_not_evaluated_row_never_becomes_a_query(self) -> None:
        # F1: an input not_evaluated row passes every mechanical check, stays
        # not_evaluated in the frozen panel, and run_panel_blastp must not
        # turn it into a query or an evidence row (verified filter path).
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            lines = (FIXTURE / "panel_ok.tsv").read_text(encoding="utf-8").splitlines()
            bait2 = lines[2].split("\t")
            bait2[0] = "bait5n"
            bait2[21] = "not_evaluated"
            bait2[22] = ""
            panel = root / "not_evaluated_panel.tsv"
            panel.write_text("\n".join([lines[0], lines[1], "\t".join(bait2), lines[3], lines[4]]) + "\n", encoding="utf-8")

            freeze = _freeze(root, panel)
            self.assertEqual(2, freeze["verified_count"])
            frozen = pd.read_csv(freeze["panel"], sep="\t", dtype=str, keep_default_na=False).fillna("")
            bait5n = frozen.loc[frozen["panel_id"] == "bait5n", "review_status"].iloc[0]
            self.assertEqual("not_evaluated", bait5n)

            with MockBlastTools():
                result = run_panel_blastp(str(freeze["manifest"]), str(FIXTURE / "interval_proteins.fa"), str(root / "search"))

            self.assertEqual("evaluable", result["provider_status"])
            evidence = pd.read_csv(result["evidence"], sep="\t", dtype=str).fillna("")
            self.assertEqual(["bait1", "bait3"], sorted(evidence["panel_id"]))
            self.assertNotIn("bait5n", list(evidence["panel_id"]))
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertIn("bait5n", manifest["skipped_review_required"])
            queries = (root / "search" / "queries.fa").read_text(encoding="utf-8")
            self.assertNotIn("bait5n", queries)


class AvailabilityTests(unittest.TestCase):
    def test_missing_blastp_is_blast_unavailable(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            freeze = _freeze(root)
            with mock.patch("candidate_qtl.blast_evidence.shutil.which", return_value=None):
                result = run_panel_blastp(str(freeze["manifest"]), str(FIXTURE / "interval_proteins.fa"), str(root / "search"))

            self.assertEqual("not_evaluated_blast_unavailable", result["provider_status"])
            self.assertEqual(0, result["exit_code"])
            self.assertTrue(Path(result["eligible_hits"]).exists())
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual("not_evaluated_blast_unavailable", manifest["provider_status"])
            self.assertTrue(manifest["reason"])

    def test_makeblastdb_failure_is_provider_failure(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            freeze = _freeze(root)
            with MockBlastTools(makeblastdb_returncode=1):
                result = run_panel_blastp(str(freeze["manifest"]), str(FIXTURE / "interval_proteins.fa"), str(root / "search"))

            self.assertEqual("not_evaluated_provider_failure", result["provider_status"])
            self.assertEqual(0, result["exit_code"])
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual("not_evaluated_provider_failure", manifest["provider_status"])
            self.assertIn("boom", manifest["reason"])

    def test_blastp_failure_is_provider_failure_not_a_miss(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            freeze = _freeze(root)
            with MockBlastTools(blastp_returncode=1):
                result = run_panel_blastp(str(freeze["manifest"]), str(FIXTURE / "interval_proteins.fa"), str(root / "search"))

            self.assertEqual("not_evaluated_provider_failure", result["provider_status"])
            self.assertEqual(0, result["exit_code"])
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertIn("boom", manifest["reason"])

    def test_unparseable_raw_hits_is_provider_failure_not_a_crash(self) -> None:
        # F6: a provider that writes malformed raw hits (wrong column count)
        # must degrade to not_evaluated_provider_failure with a written
        # manifest and header-only evidence, never raise.
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            freeze = _freeze(root)
            with MockBlastTools(raw_hits="bait1\tg1\t100.000\n"):
                result = run_panel_blastp(str(freeze["manifest"]), str(FIXTURE / "interval_proteins.fa"), str(root / "search"))

            self.assertEqual("not_evaluated_provider_failure", result["provider_status"])
            self.assertEqual(0, result["exit_code"])
            self.assertEqual(0, result["rows_evaluated"])
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual("not_evaluated_provider_failure", manifest["provider_status"])
            self.assertIn("unparseable blastp raw hits", manifest["reason"])
            self.assertTrue(manifest["reason"])
            evidence = pd.read_csv(result["evidence"], sep="\t", dtype=str)
            self.assertEqual(0, len(evidence))
            self.assertEqual(BLAST_EVIDENCE_COLUMNS, list(evidence.columns))

    def test_verify_failure_propagates_before_any_tool_runs(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            freeze = _freeze(root)
            panel_path = Path(freeze["panel"])
            panel_path.write_text(panel_path.read_text(encoding="utf-8").replace("CDLYQ", "CDLYR", 1), encoding="utf-8")

            with mock.patch("candidate_qtl.blast_evidence.subprocess.run") as sub_run, mock.patch(
                "candidate_qtl.blast_evidence.run_tool"
            ) as tool_run:
                with self.assertRaises(ValueError):
                    run_panel_blastp(str(freeze["manifest"]), str(FIXTURE / "interval_proteins.fa"), str(root / "search"))
            sub_run.assert_not_called()
            tool_run.assert_not_called()

    def test_provider_status_vocab(self) -> None:
        self.assertEqual(
            ("evaluable", "not_evaluated_no_verified_bait", "not_evaluated_blast_unavailable", "not_evaluated_provider_failure"),
            PROVIDER_STATUS_VOCAB,
        )


class HappyPathTests(unittest.TestCase):
    @unittest.skipUnless(shutil.which("blastp") and shutil.which("makeblastdb"), "BLAST+ not installed")
    def test_native_mixed_panel_preserves_long_bait_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            long_freeze = _freeze(root / "long", _panel_with_ids(root, {"bait1"}))
            mixed_freeze = _freeze(root / "mixed", _panel_with_ids(root, {"bait1", "bait3"}))
            proteins = FIXTURE / "interval_proteins.fa"
            alone = run_panel_blastp(str(long_freeze["manifest"]), str(proteins), str(root / "alone"), runner="native")
            mixed = run_panel_blastp(str(mixed_freeze["manifest"]), str(proteins), str(root / "mixed"), runner="native")
            self.assertEqual("evaluable", alone["provider_status"])
            self.assertEqual("evaluable", mixed["provider_status"])
            alone_row = pd.read_csv(alone["evidence"], sep="\t", dtype=str).fillna("").query("panel_id == 'bait1'").to_dict("records")[0]
            mixed_row = pd.read_csv(mixed["evidence"], sep="\t", dtype=str).fillna("").query("panel_id == 'bait1'").to_dict("records")[0]
            self.assertEqual(alone_row, mixed_row)

    def test_mixed_panel_keeps_each_query_in_its_own_blast_mode(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            long_freeze = _freeze(root / "long", _panel_with_ids(root, {"bait1"}))
            mixed_freeze = _freeze(root / "mixed", _panel_with_ids(root, {"bait1", "bait3"}))
            with ModeAwareBlastTools() as alone_tools:
                alone = run_panel_blastp(str(long_freeze["manifest"]), str(FIXTURE / "interval_proteins.fa"), str(root / "alone"))
            with ModeAwareBlastTools() as mixed_tools:
                mixed = run_panel_blastp(str(mixed_freeze["manifest"]), str(FIXTURE / "interval_proteins.fa"), str(root / "mixed"))
            alone_rows = pd.read_csv(alone["evidence"], sep="\t", dtype=str).fillna("")
            mixed_rows = pd.read_csv(mixed["evidence"], sep="\t", dtype=str).fillna("")
            self.assertEqual(alone_rows.to_dict("records")[0], mixed_rows.loc[mixed_rows["panel_id"] == "bait1"].to_dict("records")[0])
            self.assertEqual(["standard"], [mode for mode, _ in alone_tools.calls])
            self.assertEqual(["standard", "short"], [mode for mode, _ in mixed_tools.calls])
            self.assertNotIn("-matrix", alone_tools.calls[0][1])
            self.assertIn("-matrix", mixed_tools.calls[1][1])
            manifest = json.loads(Path(mixed["manifest"]).read_text(encoding="utf-8"))
            self.assertIn("sha256", manifest["query_fasta_fingerprint"])
            self.assertIn("sha256", manifest["groups"]["standard"]["query_fasta_fingerprint"])
            self.assertRegex(manifest["groups"]["standard"]["query_fasta_sha256"], r"^[0-9a-f]{64}$")
            self.assertEqual("standard", manifest["query_modes"]["bait1"])
            self.assertEqual("short", manifest["query_modes"]["bait3"])
            self.assertEqual(2, len(manifest["groups"]))
            self.assertEqual(sha256(Path(mixed["manifest"]).parent / "panel_vs_interval.tsv"), manifest["raw_hits_sha256"])
            self.assertEqual(manifest["raw_hits_sha256"], manifest["groups"]["standard"]["merged_raw_hits_sha256"])

    def test_homogeneous_panel_makes_no_empty_group_call(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            freeze = _freeze(root, _panel_with_ids(root, {"bait1"}))
            with ModeAwareBlastTools() as tools:
                result = run_panel_blastp(str(freeze["manifest"]), str(FIXTURE / "interval_proteins.fa"), str(root / "search"))
            self.assertEqual("evaluable", result["provider_status"])
            self.assertEqual(["standard"], [mode for mode, _ in tools.calls])
            self.assertFalse((root / "search" / "queries.short.fa").exists())

    def test_mixed_group_failure_is_provider_failure_with_group_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            freeze = _freeze(root, _panel_with_ids(root, {"bait1", "bait3"}))
            with ModeAwareBlastTools(fail_mode="short"):
                result = run_panel_blastp(str(freeze["manifest"]), str(FIXTURE / "interval_proteins.fa"), str(root / "search"))
            self.assertEqual("not_evaluated_provider_failure", result["provider_status"])
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual("evaluable", manifest["groups"]["standard"]["provider_status"])
            self.assertEqual("not_evaluated_provider_failure", manifest["groups"]["short"]["provider_status"])
            self.assertEqual(0, result["rows_evaluated"])

    def test_zero_exit_without_raw_output_is_provider_failure(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            freeze = _freeze(root, _panel_with_ids(root, {"bait1"}))
            with MissingRawBlastTools():
                result = run_panel_blastp(str(freeze["manifest"]), str(FIXTURE / "interval_proteins.fa"), str(root / "search"))
            self.assertEqual("not_evaluated_provider_failure", result["provider_status"])
            self.assertIn("without creating raw output", result["reason"])

    def test_evaluable_search_with_mock_blastp(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            freeze = _freeze(root)
            with MockBlastTools():
                result = run_panel_blastp(str(freeze["manifest"]), str(FIXTURE / "interval_proteins.fa"), str(root / "search"))

            self.assertEqual("evaluable", result["provider_status"])
            self.assertEqual(3, result["rows_evaluated"])
            self.assertEqual(1, result["unique_count"])
            self.assertEqual(1, result["ambiguous_count"])
            self.assertEqual(1, result["no_threshold_count"])
            self.assertEqual(0, result["exit_code"])

            evidence = pd.read_csv(result["evidence"], sep="\t", dtype=str).fillna("")
            self.assertEqual(BLAST_EVIDENCE_COLUMNS, list(evidence.columns))
            by_id = {str(row["panel_id"]): row for row in evidence.to_dict("records")}

            bait1 = by_id["bait1"]
            self.assertEqual("verified", bait1["review_status"])
            self.assertEqual("evaluable", bait1["provider_status"])
            self.assertEqual("unique", bait1["ambiguity_status"])
            self.assertEqual("g1", bait1["top_hit_target_id"])
            self.assertEqual("3", bait1["hit_count"])
            self.assertEqual("2", bait1["paralog_count"])
            self.assertEqual("none", bait1["rank_contribution"])
            # F8: in the mixed panel (short bait3 forces PAM30 mode), the
            # >=100 aa rows still record a non-empty raw evalue.
            self.assertNotEqual("", bait1["top_hit_evalue"])

            bait2 = by_id["bait2"]
            self.assertEqual("ambiguous_paralog", bait2["ambiguity_status"])
            self.assertEqual("2", bait2["hit_count"])
            self.assertEqual("1", bait2["paralog_count"])

            bait3 = by_id["bait3"]
            self.assertEqual("no_threshold_hit", bait3["ambiguity_status"])
            self.assertEqual("0", bait3["hit_count"])
            self.assertEqual("0", bait3["paralog_count"])

            self.assertNotIn("bait4", by_id)
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual(sha256(Path(result["eligible_interval_hits"])), result["eligible_interval_hits_sha256"])
            self.assertEqual(result["eligible_interval_hits_sha256"], manifest["eligible_interval_hits_sha256"])
            self.assertEqual(["bait4"], manifest["skipped_review_required"])
            self.assertEqual(3, manifest["counts"]["rows_evaluated"])
            self.assertEqual(1, manifest["counts"]["ambiguous_paralog"])
            self.assertEqual(1, manifest["counts"]["no_threshold_hit"])
            self.assertEqual("post_freeze_similarity_only", manifest["interval_usage"])
            self.assertNotIn("interval_accessed", manifest)
            self.assertEqual("unit", manifest["case_id"])

    def test_short_query_mode_arguments(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            freeze = _freeze(root)
            with MockBlastTools() as tools:
                run_panel_blastp(str(freeze["manifest"]), str(FIXTURE / "interval_proteins.fa"), str(root / "search"))

            args = tools.blastp_args
            self.assertEqual("2", args[args.index("-word_size") + 1])
            self.assertEqual("PAM30", args[args.index("-matrix") + 1])
            self.assertEqual("9", args[args.index("-gapopen") + 1])
            self.assertEqual("1", args[args.index("-gapextend") + 1])
            self.assertEqual("no", args[args.index("-seg") + 1])
            self.assertEqual("1000", args[args.index("-evalue") + 1])
            self.assertEqual("5000", args[args.index("-max_target_seqs") + 1])
            # Multi-HSP rows are retained and union-aggregated by the parser;
            # limiting the provider to one HSP would lose target coverage.
            self.assertNotIn("-max_hsps", args)

    def test_search_is_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            freeze = _freeze(root)
            with MockBlastTools():
                first = run_panel_blastp(str(freeze["manifest"]), str(FIXTURE / "interval_proteins.fa"), str(root / "s1"))
            with MockBlastTools():
                second = run_panel_blastp(str(freeze["manifest"]), str(FIXTURE / "interval_proteins.fa"), str(root / "s2"))

            self.assertEqual(sha256(Path(first["evidence"])), sha256(Path(second["evidence"])))
            self.assertEqual(first["evidence"], second["evidence"].replace("s2", "s1"))

    def test_duplicate_sseqid_hsp_rows_are_union_aggregated(self) -> None:
        # Two HSP rows for the same target remain one target hit.  Identity is
        # conservative (minimum), while E-value/bitscore retain the strongest
        # value and the long-form row records hsp_count and union segments.
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            freeze = _freeze(root)
            raw = (
                "bait1\tg1\t100.000\t120\t120\t120\t1\t120\t1\t120\t0.0\t120.0\n"
                "bait1\tg1\t90.000\t90\t120\t120\t1\t90\t1\t90\t1e-30\t90.0\n"
            )
            with MockBlastTools(raw_hits=raw):
                result = run_panel_blastp(str(freeze["manifest"]), str(FIXTURE / "interval_proteins.fa"), str(root / "search"))

            self.assertEqual("evaluable", result["provider_status"])
            evidence = pd.read_csv(result["evidence"], sep="\t", dtype=str).fillna("")
            by_id = {str(row["panel_id"]): row for row in evidence.to_dict("records")}
            bait1 = by_id["bait1"]
            self.assertEqual("1", bait1["hit_count"])
            self.assertEqual("0", bait1["paralog_count"])
            self.assertEqual("120.0", bait1["top_hit_bitscore"])
            self.assertEqual("90.0", bait1["top_hit_identity"])
            eligible = pd.read_csv(result["eligible_interval_hits"], sep="\t", dtype=str).fillna("")
            g1 = eligible.loc[eligible["target_id"] == "g1"].iloc[0]
            self.assertEqual("2", g1["hsp_count"])
            self.assertEqual("1.0", g1["query_coverage"])
            self.assertEqual("1-120", g1["query_union_segments"])

    def test_disjoint_hsps_use_coordinate_union_for_threshold(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            raw_path = Path(tempdir) / "raw.tsv"
            raw_path.write_text(
                "bait1\ttarget\t99.0\t584\t1200\t1192\t1\t584\t1\t584\t1e-30\t400\n"
                "bait1\ttarget\t98.0\t584\t1200\t1192\t586\t1169\t586\t1169\t1e-20\t390\n",
                encoding="utf-8",
            )
            hits = _parse_blast_hits(raw_path)
            self.assertEqual(1, len(hits))
            hit = hits[0]
            self.assertEqual(2, hit["hsp_count"])
            self.assertAlmostEqual(0.979865, float(hit["target_coverage"]), places=5)
            self.assertEqual("1-584;586-1169", hit["target_union_segments"])
            self.assertEqual(98.0, hit["pident"])
            self.assertEqual(1e-30, hit["evalue"])
            self.assertEqual(400.0, hit["bitscore"])


class ClassifyAmbiguityTests(unittest.TestCase):
    def test_ambiguity_vocab(self) -> None:
        self.assertEqual(("unique", "ambiguous_paralog", "no_threshold_hit", "not_evaluated"), AMBIGUITY_VOCAB)

    def test_no_hits_is_no_threshold_hit(self) -> None:
        status, top = _classify_ambiguity([], DEFAULT_THRESHOLDS)
        self.assertEqual("no_threshold_hit", status)
        self.assertIsNone(top)

    def test_single_hit_is_unique(self) -> None:
        hit = _hit()
        status, top = _classify_ambiguity([hit], DEFAULT_THRESHOLDS)
        self.assertEqual("unique", status)
        self.assertIs(top, hit)

    def test_second_within_tie_ratio_is_ambiguous(self) -> None:
        hits = [_hit(sseqid="g1", bitscore=200.0), _hit(sseqid="g2", bitscore=185.0)]
        status, top = _classify_ambiguity(hits, DEFAULT_THRESHOLDS)
        self.assertEqual("ambiguous_paralog", status)
        self.assertEqual("g1", top["sseqid"])

    def test_clear_winner_is_unique(self) -> None:
        hits = [_hit(sseqid="g1", bitscore=200.0), _hit(sseqid="g2", bitscore=120.0)]
        status, top = _classify_ambiguity(hits, DEFAULT_THRESHOLDS)
        self.assertEqual("unique", status)
        self.assertEqual("g1", top["sseqid"])


class EvaluateBaitTests(unittest.TestCase):
    ROW = {"panel_id": "bait1", "query_sequence": S, "query_sequence_sha256": "a" * 64}
    TOOL_META = {"blastp_version": "blastp: 2.12.0+", "db_fingerprint": "db" * 32, "interval_proteins_sha256": "iv" * 32}

    def _evaluate(self, hits: list[dict[str, object]], short_query_mode: bool = False) -> dict[str, object]:
        return _evaluate_bait(self.ROW, hits, DEFAULT_THRESHOLDS, self.TOOL_META, short_query_mode)

    def test_cross_species_identity_at_threshold_is_retained(self) -> None:
        row = self._evaluate([_hit(pident=30.0)])
        self.assertEqual("unique", row["ambiguity_status"])
        self.assertEqual("g1", row["top_hit_target_id"])

    def test_identity_below_threshold_is_dropped(self) -> None:
        row = self._evaluate([_hit(pident=29.999)])
        self.assertEqual("no_threshold_hit", row["ambiguity_status"])
        self.assertEqual("", row["top_hit_target_id"])
        self.assertEqual("", row["top_hit_identity"])

    def test_query_coverage_below_threshold_is_dropped(self) -> None:
        row = self._evaluate([_hit(qend=71)])  # 71/120 = 0.592 < 0.6
        self.assertEqual("no_threshold_hit", row["ambiguity_status"])

    def test_target_coverage_below_threshold_is_dropped(self) -> None:
        row = self._evaluate([_hit(send=72, slen=122)])  # 72/122 = 0.590 < 0.6
        self.assertEqual("no_threshold_hit", row["ambiguity_status"])

    def test_evalue_above_threshold_is_dropped_for_full_length_query(self) -> None:
        row = self._evaluate([_hit(evalue=1e-3)], short_query_mode=False)  # qlen 120 >= 100
        self.assertEqual("no_threshold_hit", row["ambiguity_status"])

    def test_evalue_above_threshold_is_dropped_for_short_query(self) -> None:
        # PAM30 improves sensitivity but must not turn a high-E-value random
        # short hit into a unique Route A positive.
        row = self._evaluate([_hit(evalue=1e-3, qlen=48, slen=48, length=48)], short_query_mode=True)
        self.assertEqual("no_threshold_hit", row["ambiguity_status"])

    def test_gappy_alignment_uses_query_and_target_residue_spans(self) -> None:
        row = self._evaluate([_hit(length=120, qend=90, sstart=10, send=99)])
        self.assertEqual("unique", row["ambiguity_status"])
        self.assertEqual(0.75, row["top_hit_query_coverage"])
        self.assertEqual(0.75, row["top_hit_target_coverage"])

    def test_hit_metrics_and_paralog_count(self) -> None:
        row = self._evaluate(
            [
                _hit(sseqid="g1", bitscore=250.0),
                _hit(sseqid="g2", bitscore=235.0),
                _hit(sseqid="g3", bitscore=230.0),
            ]
        )
        self.assertEqual(3, row["hit_count"])
        self.assertEqual(2, row["paralog_count"])
        self.assertEqual("ambiguous_paralog", row["ambiguity_status"])
        self.assertEqual(100.0, row["top_hit_identity"])
        self.assertEqual(1.0, row["top_hit_query_coverage"])
        self.assertEqual(1.0, row["top_hit_target_coverage"])
        self.assertEqual(250.0, row["top_hit_bitscore"])
        self.assertEqual("none", row["rank_contribution"])

    def test_evidence_row_has_19_columns_in_declared_order(self) -> None:
        self.assertEqual(19, len(BLAST_EVIDENCE_COLUMNS))
        row = self._evaluate([_hit()])
        self.assertEqual(BLAST_EVIDENCE_COLUMNS, list(row.keys()))
        self.assertEqual(120, row["query_length"])
        self.assertEqual("verified", row["review_status"])
        self.assertEqual("evaluable", row["provider_status"])
        self.assertEqual(self.TOOL_META["blastp_version"], row["blastp_version"])
        self.assertEqual(self.TOOL_META["db_fingerprint"], row["db_fingerprint"])
        self.assertEqual(self.TOOL_META["interval_proteins_sha256"], row["interval_proteins_sha256"])
        self.assertEqual("none", row["rank_contribution"])
        self.assertEqual("", row["note"])


if __name__ == "__main__":
    unittest.main()
