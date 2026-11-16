from __future__ import annotations

import io
import json
import os
import tempfile
import time
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock

import pandas as pd

from candidate_qtl.budget_account import (
    BUDGET_PROFILES,
    LEDGER_INTEGER_FIELDS,
    LEDGER_LIMIT_FIELDS,
    LEDGER_STATUSES,
    LEDGER_USAGE_FIELDS,
)
from candidate_qtl.budget_engine import GATING_HEADER_COMMENT, BlastpBackend, _write_json, run_budget_engine
from candidate_qtl.cli import main
from candidate_qtl.core import ExternalDependencyError
from candidate_qtl.tier0_inventory import load_provider_table

CASE_ID = "engine-case"


class FakeBackend:
    def __init__(self, hits=None, cpu=0.01, hit_gene_ids=None):
        self.available = True
        self.cache = None
        self.search_calls = 0
        self.seen_timeouts = []
        self.cpu = cpu
        self.hits = hits
        # None = hit for every queried gene; a set restricts hits to those ids.
        self.hit_gene_ids = set(hit_gene_ids) if hit_gene_ids is not None else None

    def search(self, query_fasta, db, bin, timeout=None):
        self.search_calls += 1
        self.seen_timeouts.append(timeout)
        if self.hits is not None:
            return self.hits, self.cpu
        qseqid = Path(query_fasta).stem
        if self.hit_gene_ids is None or qseqid in self.hit_gene_ids:
            return [{"qseqid": qseqid, "sseqid": f"hit_{qseqid}", "pident": "95", "evalue": "1e-50", "bitscore": "200"}], self.cpu
        return [], self.cpu


def _engine_genes():
    return [
        {
            "id": "G1", "start": 11, "end": 22, "cds_start": 11, "cds_end": 22, "cds": "ATGGAAAACTAA", "tx": "Tx1",
            "attrs": "ID=G1;Name=GeneA;Note=annotated%20gene;Dbxref=GeneID:1;Ontology_term=GO:0000001",
        },
        {"id": "G2", "start": 31, "end": 45, "cds_start": 31, "cds_end": 45, "cds": "ATGGAAAACGCCTAA", "tx": "Tx2", "attrs": "ID=G2"},
        {"id": "G3", "start": 51, "end": 60, "cds": None, "tx": None, "attrs": "ID=G3"},
        {"id": "G4", "start": 71, "end": 80, "cds_start": 71, "cds_end": 80, "cds": "ATGGAAAACA", "tx": "Tx4", "attrs": "ID=G4"},
        {
            "id": "G5", "start": 91, "end": 105, "cds_start": 91, "cds_end": 105, "cds": "ATGGAAAACGATTAA", "tx": "Tx5",
            "attrs": "ID=G5;Name=GeneC;Note=desc5",
        },
        {"id": "G6", "start": 111, "end": 125, "cds_start": 111, "cds_end": 125, "cds": "ATGGAAAACGATTAA", "tx": "Tx6", "attrs": "ID=G6;Name=GeneA"},
        {"id": "G7", "start": 131, "end": 145, "cds_start": 131, "cds_end": 145, "cds": "ATGGAAAACGATTAA", "tx": "Tx7", "attrs": "ID=G7;Name=GeneA"},
        {"id": "G8", "start": 151, "end": 165, "cds_start": 151, "cds_end": 165, "cds": "ATGGAAAACGATTAA", "tx": "Tx8", "attrs": "ID=G8;Name=GeneA"},
        {
            "id": "G9", "start": 171, "end": 185, "cds_start": 171, "cds_end": 185, "cds": "ATGGAAAACGATTAA", "tx": "Tx9",
            "attrs": "ID=G9;Name=GeneB;Note=desc9",
        },
    ]


def _write_extract_fixture(root: Path, genes, length: int = 220) -> tuple[Path, Path]:
    seq = ["G"] * length
    for gene in genes:
        if gene["cds"]:
            start = gene["cds_start"] - 1
            for offset, base in enumerate(gene["cds"]):
                seq[start + offset] = base
    fasta = root / "genome.fa"
    fasta.write_text(">chr1\n" + "".join(seq) + "\n", encoding="utf-8")
    gff = root / "genes.gff3"
    lines = ["##gff-version 3"]
    for gene in genes:
        lines.append(f"chr1\tfixture\tgene\t{gene['start']}\t{gene['end']}\t.\t+\t.\t{gene['attrs']}")
        if gene["cds"]:
            lines.append(f"chr1\tfixture\tmRNA\t{gene['start']}\t{gene['end']}\t.\t+\t.\tID={gene['tx']};Parent={gene['id']}")
            lines.append(f"chr1\tfixture\tCDS\t{gene['cds_start']}\t{gene['cds_end']}\t.\t+\t0\tParent={gene['tx']}")
    gff.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return fasta, gff


def _engine_fixture(root: Path, testcase: unittest.TestCase) -> tuple[Path, dict[str, dict[str, object]], Path]:
    fasta, gff = _write_extract_fixture(root, _engine_genes())
    extract_out = root / "extract"
    testcase.assertEqual(
        0,
        main(["extract", "--fasta", str(fasta), "--gff3", str(gff), "--interval", "chr1:1-1000", "--out-dir", str(extract_out)]),
    )
    (root / "expression.tsv").write_text("gene_id\tvalue\tsource_id\nG5\t10\tE1\nG5\t20\tE1\nG2\t5\tE2\n", encoding="utf-8")
    (root / "domains.tsv").write_text("gene_id\tdomain_id\nG5\tPF001\nG5\tPF002\n", encoding="utf-8")
    (root / "orthogroup.tsv").write_text("gene_id\torthogroup\nG5\tOG0001\n", encoding="utf-8")
    (root / "variant.tsv").write_text("gene_id\teffect\tsource_id\nG5\tmissense\tV1\n", encoding="utf-8")
    (root / "family.tsv").write_text("gene_id\tfamily_id\nG5\tFAM001\n", encoding="utf-8")
    providers = {
        "expression": load_provider_table("expression", str(root / "expression.tsv"), "expr-db", "v2"),
        "domains": load_provider_table("domains", str(root / "domains.tsv"), "dom-db", "v1"),
        "orthogroup": load_provider_table("orthogroup", str(root / "orthogroup.tsv"), "egg-db", "v1"),
        "variant": load_provider_table("variant", str(root / "variant.tsv"), "var-db", "v1"),
        "family": load_provider_table("family", str(root / "family.tsv"), "fam-db", "v1"),
    }
    protein_db = root / "proteome.fa"
    protein_db.write_text(">ref1\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ\n", encoding="utf-8")
    return extract_out, providers, protein_db


class BudgetEngineTests(unittest.TestCase):
    def test_determinism_two_runs_byte_identical(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out, providers, protein_db = _engine_fixture(root, self)
            out_a = root / "run_a"
            out_b = root / "run_b"
            run_budget_engine(
                extract_out, out_a, CASE_ID, BUDGET_PROFILES["student_cpu"], providers,
                protein_db=protein_db, search_backend=FakeBackend(),
            )
            run_budget_engine(
                extract_out, out_b, CASE_ID, BUDGET_PROFILES["student_cpu"], providers,
                protein_db=protein_db, search_backend=FakeBackend(),
            )
            for name in (
                "tier0_inventory.tsv", "tier1_clusters.tsv", "tier2_plan.tsv", "budget_handoff.tsv",
                "tier0_manifest.json", "tier1_manifest.json", "tier2_manifest.json",
                "deterministic_top20.tsv", "llm_assisted_top20.tsv",
            ):
                self.assertEqual((out_a / name).read_bytes(), (out_b / name).read_bytes(), name)
            # budget_manifest.json now carries runtime wall/memory consumption
            # (charged per query by the Tier 2 loop) and runtime
            # truncation-point remaining_before floats: compare structurally
            # with the wall/peak consumed/remaining and every
            # truncation_points[].remaining_before zeroed, like the benchmark
            # carve-out.
            left_manifest = json.loads((out_a / "budget_manifest.json").read_text(encoding="utf-8"))
            right_manifest = json.loads((out_b / "budget_manifest.json").read_text(encoding="utf-8"))
            for entry in (left_manifest, right_manifest):
                for limit in ("wall_time_limit_seconds", "peak_memory_limit_mb"):
                    entry["budgets"][limit]["consumed"] = 0
                    entry["budgets"][limit]["remaining"] = 0
                for point in entry["unprocessed"]["truncation_points"]:
                    point["remaining_before"] = 0
            self.assertEqual(left_manifest, right_manifest)
            left = json.loads((out_a / "benchmark_resource_usage.json").read_text(encoding="utf-8"))
            right = json.loads((out_b / "benchmark_resource_usage.json").read_text(encoding="utf-8"))
            for entry in (left, right):
                for field in ("wall_time_seconds", "peak_memory_mb", "sequence_search_cpu_seconds"):
                    entry["usage"][field] = 0
            self.assertEqual(left, right)

    def test_zero_budget_full_inventory_and_completed_status(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out, providers, protein_db = _engine_fixture(root, self)
            out = root / "run"
            manifest = run_budget_engine(
                extract_out, out, CASE_ID, BUDGET_PROFILES["zero"], providers,
                protein_db=protein_db, search_backend=FakeBackend(),
            )
            inventory = pd.read_csv(out / "tier0_inventory.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(9, len(inventory))
            handoff = pd.read_csv(out / "budget_handoff.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(9, len(handoff))
            self.assertEqual("completed", manifest["status"])
            self.assertEqual("zero", manifest["budget_profile"])
            usage = json.loads((out / "benchmark_resource_usage.json").read_text(encoding="utf-8"))["usage"]
            self.assertEqual(0, usage["sequence_query_count"])
            self.assertEqual(0, usage["searched_protein_count"])
            self.assertEqual(9, usage["unreviewed_candidate_count"])
            self.assertEqual(0, usage["reviewed_candidate_count"])
            # cache-first execution: exactly one lookup happens before the
            # zero-budget gate stops the loop
            self.assertEqual({"hit_count": 0, "miss_count": 1, "corrupt_reject_count": 0, "stale_reject_count": 0, "cache_status": "enabled"}, manifest["cache"])
            self.assertEqual(9, manifest["unprocessed"]["unreviewed_candidate_count"])
            benchmark = json.loads((out / "benchmark_resource_usage.json").read_text(encoding="utf-8"))
            self.assertEqual("deterministic", benchmark["benchmark_arm"])
            self.assertEqual("completed", benchmark["status"])
            self.assertEqual({name: 0 for name in LEDGER_LIMIT_FIELDS}, benchmark["limits"])

    def test_mid_run_exhaustion_reports_unprocessed(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out, providers, protein_db = _engine_fixture(root, self)
            limits = {"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 1, "wall_time_limit_seconds": 3600, "peak_memory_limit_mb": 8192}
            out = root / "run"
            manifest = run_budget_engine(
                extract_out, out, CASE_ID, limits, providers,
                protein_db=protein_db, flag_genes=["G1"], search_backend=FakeBackend(),
            )
            self.assertEqual("budget_exhausted", manifest["status"])
            plan = pd.read_csv(out / "tier2_plan.tsv", sep="\t", dtype=str).fillna("")
            deferred = plan[plan["query_status"] == "deferred_budget"]
            self.assertGreater(len(deferred), 0)
            self.assertTrue((deferred["deferred_reason"] == "budget_exhausted").all())
            self.assertEqual("executed", plan[plan["query_status"] == "executed"]["query_status"].iloc[0])
            self.assertGreater(manifest["unprocessed"]["deferred_query_count"], 0)
            self.assertGreater(len(manifest["unprocessed"]["truncation_points"]), 0)
            self.assertEqual("sequence_query", manifest["unprocessed"]["truncation_points"][0]["kind"])
            self.assertGreater(manifest["unprocessed"]["unreviewed_candidate_count"], 0)
            benchmark = json.loads((out / "benchmark_resource_usage.json").read_text(encoding="utf-8"))
            self.assertEqual("budget_exhausted", benchmark["status"])
            self.assertEqual("custom", manifest["budget_profile"])
            self.assertEqual(1, benchmark["usage"]["sequence_query_count"])

    def test_corrupt_cache_entry_recomputed_identically(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out, providers, protein_db = _engine_fixture(root, self)
            out = root / "run"
            run_budget_engine(
                extract_out, out, CASE_ID, BUDGET_PROFILES["student_cpu"], providers,
                protein_db=protein_db, search_backend=FakeBackend(),
            )
            before = {name: (out / name).read_bytes() for name in ("tier0_inventory.tsv", "tier1_clusters.tsv", "tier2_plan.tsv", "budget_handoff.tsv", "tier2_manifest.json")}
            payloads = sorted((out / "budget_cache").glob("*.payload"))
            self.assertGreater(len(payloads), 0)
            corrupted = payloads[0]
            data = bytearray(corrupted.read_bytes())
            data[-1] ^= 0xFF
            corrupted.write_bytes(bytes(data))
            run_budget_engine(
                extract_out, out, CASE_ID, BUDGET_PROFILES["student_cpu"], providers,
                protein_db=protein_db, search_backend=FakeBackend(),
            )
            manifest = json.loads((out / "budget_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(1, manifest["cache"]["corrupt_reject_count"])
            for name, payload in before.items():
                self.assertEqual(payload, (out / name).read_bytes(), name)

    def test_cache_reuse_skips_searches(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out, providers, protein_db = _engine_fixture(root, self)
            fake = FakeBackend()
            out = root / "run"
            run_budget_engine(
                extract_out, out, CASE_ID, BUDGET_PROFILES["student_cpu"], providers,
                protein_db=protein_db, search_backend=fake,
            )
            self.assertEqual(3, fake.search_calls)
            cold_usage = json.loads((out / "benchmark_resource_usage.json").read_text(encoding="utf-8"))["usage"]
            self.assertEqual(3 * fake.cpu, cold_usage["sequence_search_cpu_seconds"])
            run_budget_engine(
                extract_out, out, CASE_ID, BUDGET_PROFILES["student_cpu"], providers,
                protein_db=protein_db, search_backend=fake,
            )
            self.assertEqual(3, fake.search_calls)
            benchmark = json.loads((out / "benchmark_resource_usage.json").read_text(encoding="utf-8"))
            self.assertEqual(3, benchmark["usage"]["cache_hit_count"])
            self.assertEqual(0, benchmark["usage"]["cache_miss_count"])
            self.assertEqual(0, benchmark["usage"]["sequence_query_count"])
            # warm rerun records the same measured cpu cost as the cold run
            self.assertEqual(cold_usage["sequence_search_cpu_seconds"], benchmark["usage"]["sequence_search_cpu_seconds"])
            self.assertIn("cache_hit_count", benchmark["usage"])
            self.assertIn("cache_miss_count", benchmark["usage"])

    def test_zero_provider_flags_statuses_not_evaluated_handoff_complete(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out, _providers, protein_db = _engine_fixture(root, self)
            out = root / "run"
            manifest = run_budget_engine(
                extract_out, out, CASE_ID, BUDGET_PROFILES["student_cpu"], {},
                protein_db=protein_db, search_backend=FakeBackend(),
            )
            tier0 = json.loads((out / "tier0_manifest.json").read_text(encoding="utf-8"))
            for column in ("domain_status", "orthogroup_status", "expression_status", "variant_status", "sv_status", "pav_status", "promoter_status"):
                self.assertEqual({"not_evaluated": 9}, tier0["counts"][column])
            handoff = pd.read_csv(out / "budget_handoff.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(9, len(handoff))
            self.assertTrue((handoff["review_state"] != "").all())
            benchmark = json.loads((out / "benchmark_resource_usage.json").read_text(encoding="utf-8"))
            self.assertEqual(4, benchmark["usage"]["reviewed_candidate_count"])
            self.assertEqual("completed", manifest["status"])

    def test_exhaustive_arm_without_provider_all_not_run(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out, providers, _protein_db = _engine_fixture(root, self)
            out = root / "run"
            manifest = run_budget_engine(
                extract_out, out, CASE_ID, BUDGET_PROFILES["student_cpu"], providers,
                exhaustive=True, search_backend=FakeBackend(),
            )
            comparison = (out / "budget_gating_comparison.tsv").read_text(encoding="utf-8")
            self.assertTrue(comparison.startswith(GATING_HEADER_COMMENT + "\n"))
            frame = pd.read_csv(io.StringIO(comparison), sep="\t", dtype=str, skiprows=1).fillna("")
            self.assertEqual(9, len(frame))
            self.assertTrue((frame["hierarchical_hit_presence"] == "not_run").all())
            self.assertTrue((frame["exhaustive_hit_presence"] == "not_run").all())
            self.assertTrue((frame["gating_gap"] == "no").all())
            self.assertEqual(0, manifest["gating_added_false_negative_count"])
            self.assertEqual("exhaustive_arm_provider_missing", manifest["gating_note"])
            exhaustive_benchmark = json.loads((out / "exhaustive" / "benchmark_resource_usage.json").read_text(encoding="utf-8"))
            self.assertEqual("exhaustive", exhaustive_benchmark["benchmark_arm"])
            self.assertEqual("completed", exhaustive_benchmark["status"])

    def test_exhaustive_arm_gating_gap_detected(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out, providers, protein_db = _engine_fixture(root, self)
            limits = {"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 1, "wall_time_limit_seconds": 3600, "peak_memory_limit_mb": 8192}
            out = root / "run"
            manifest = run_budget_engine(
                extract_out, out, CASE_ID, limits, providers,
                protein_db=protein_db, flag_genes=["G1"], exhaustive=True, search_backend=FakeBackend(),
            )
            frame = pd.read_csv(out / "budget_gating_comparison.tsv", sep="\t", dtype=str, skiprows=1).fillna("")
            gaps = frame[frame["gating_gap"] == "yes"]
            self.assertEqual(6, len(gaps))
            self.assertEqual(len(gaps), manifest["gating_added_false_negative_count"])
            self.assertTrue((gaps["hierarchical_query_status"] != "executed").all())
            self.assertTrue((gaps["exhaustive_query_status"] == "executed").all())
            self.assertTrue((gaps["exhaustive_hit_presence"] == "hit").all())
            self.assertEqual("budget_exhausted", manifest["status"])

    def test_gating_gap_yes_when_exhaustive_searched_and_hit(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out, providers, protein_db = _engine_fixture(root, self)
            limits = {"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 1, "wall_time_limit_seconds": 3600, "peak_memory_limit_mb": 8192}
            out = root / "run"
            manifest = run_budget_engine(
                extract_out, out, CASE_ID, limits, providers,
                protein_db=protein_db, flag_genes=["G1"], exhaustive=True, search_backend=FakeBackend(),
            )
            frame = pd.read_csv(out / "budget_gating_comparison.tsv", sep="\t", dtype=str, skiprows=1).fillna("")
            # G2 is deferred by the hierarchical arm (sequence limit 1) but the
            # exhaustive arm searched it and found a hit -> a real gap.
            g2 = frame[frame["gene_id"] == "G2"].iloc[0]
            self.assertEqual("deferred_budget", g2["hierarchical_query_status"])
            self.assertEqual("executed", g2["exhaustive_query_status"])
            self.assertEqual("hit", g2["exhaustive_hit_presence"])
            self.assertEqual("yes", g2["gating_gap"])
            self.assertEqual(len(frame[frame["gating_gap"] == "yes"]), manifest["gating_added_false_negative_count"])

    def test_gating_gap_no_without_exhaustive_hit(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out, providers, protein_db = _engine_fixture(root, self)
            limits = {"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 1, "wall_time_limit_seconds": 3600, "peak_memory_limit_mb": 8192}
            out = root / "run"
            # backend hits only G1: every other exhaustive search returns no
            # hit, so those deferred genes must NOT be counted as gaps
            manifest = run_budget_engine(
                extract_out, out, CASE_ID, limits, providers,
                protein_db=protein_db, flag_genes=["G1"], exhaustive=True,
                search_backend=FakeBackend(hit_gene_ids={"G1"}),
            )
            frame = pd.read_csv(out / "budget_gating_comparison.tsv", sep="\t", dtype=str, skiprows=1).fillna("")
            gaps = frame[frame["gating_gap"] == "yes"]
            self.assertEqual(1, len(gaps))
            g1 = frame[frame["gene_id"] == "G1"].iloc[0]
            self.assertEqual("deferred_budget", g1["hierarchical_query_status"])
            self.assertEqual("executed", g1["exhaustive_query_status"])
            self.assertEqual("hit", g1["exhaustive_hit_presence"])
            self.assertEqual("yes", g1["gating_gap"])
            self.assertEqual(1, manifest["gating_added_false_negative_count"])
            searched_no_hit = frame[(frame["exhaustive_query_status"] == "executed") & (frame["exhaustive_hit_presence"] == "no_hit")]
            self.assertGreater(len(searched_no_hit), 0)
            self.assertTrue((searched_no_hit["gating_gap"] == "no").all())

    def test_rank_invariance_no_extra_files(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out, providers, protein_db = _engine_fixture(root, self)
            before = {str(path.relative_to(extract_out)): path.read_bytes() for path in extract_out.rglob("*") if path.is_file()}
            out = root / "engine_out"
            run_budget_engine(
                extract_out, out, CASE_ID, BUDGET_PROFILES["student_cpu"], providers,
                protein_db=protein_db, search_backend=FakeBackend(),
            )
            after = {str(path.relative_to(extract_out)): path.read_bytes() for path in extract_out.rglob("*") if path.is_file()}
            self.assertEqual(before, after)
            artifacts = [path.name for path in out.rglob("*") if path.is_file()]
            self.assertFalse([name for name in artifacts if name.startswith("weights_used") or name.startswith("ranked")])

    def test_cli_smoke_budget_run_zero_profile(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out, _providers, _protein_db = _engine_fixture(root, self)
            stdout = io.StringIO()
            with redirect_stdout(stdout):
                exit_code = main(["budget-run", "--extract-dir", str(extract_out), "--budget-profile", "zero", "--json"])
            self.assertEqual(0, exit_code)
            emitted = json.loads(stdout.getvalue())
            self.assertEqual("1.0", emitted["schema_version"])
            self.assertEqual("completed", emitted["status"])
            self.assertEqual(0, emitted["exit_code"])
            self.assertTrue(emitted["outputs"])
            self.assertTrue((extract_out / "budget_engine" / "budget_manifest.json").exists())

    def test_ledger_validation(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out, providers, protein_db = _engine_fixture(root, self)
            out = root / "run"
            manifest = run_budget_engine(
                extract_out, out, CASE_ID, BUDGET_PROFILES["student_cpu"], providers,
                protein_db=protein_db, search_backend=FakeBackend(),
            )
            benchmark = json.loads((out / "benchmark_resource_usage.json").read_text(encoding="utf-8"))
            self.assertEqual("1.0", benchmark["schema_version"])
            self.assertEqual(CASE_ID, benchmark["case_id"])
            self.assertEqual(set(LEDGER_LIMIT_FIELDS), set(benchmark["limits"]))
            self.assertEqual(set(LEDGER_USAGE_FIELDS), set(benchmark["usage"]))
            for field in LEDGER_INTEGER_FIELDS:
                self.assertIsInstance(benchmark["usage"][field], int)
            for name, value in benchmark["usage"].items():
                self.assertGreaterEqual(value, 0)
                self.assertTrue(pd.notna(value))
            self.assertIn(benchmark["status"], LEDGER_STATUSES)
            self.assertIn(manifest["status"], LEDGER_STATUSES)
            self.assertEqual(BUDGET_PROFILES["student_cpu"]["sequence_query_limit"], benchmark["limits"]["sequence_query_limit"])
            self.assertEqual(3, benchmark["usage"]["searched_protein_count"])
            self.assertEqual("budgeted_hierarchical", benchmark["benchmark_arm"])
            self.assertNotIn("timestamp", (out / "budget_manifest.json").read_text(encoding="utf-8"))

    def test_handoff_review_states(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out, providers, protein_db = _engine_fixture(root, self)
            out = root / "run"
            manifest = run_budget_engine(
                extract_out, out, CASE_ID, BUDGET_PROFILES["student_cpu"], providers,
                protein_db=protein_db, search_backend=FakeBackend(),
            )
            handoff = pd.read_csv(out / "budget_handoff.tsv", sep="\t", dtype=str).fillna("")
            counts = dict(handoff["review_state"].value_counts())
            self.assertEqual({"priority_review": 3, "deferred": 3, "rescue_pending_review": 3}, counts)
            self.assertEqual(6, manifest["unprocessed"]["unreviewed_candidate_count"])

    def test_zero_budget_all_deferred_still_completed(self) -> None:
        """Fix 2: a zero-budget run with every planned query deferred stays
        completed -- a genuinely zero declared budget is by-design
        work-free, not an exhaustion (see _status_for)."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out, providers, protein_db = _engine_fixture(root, self)
            out = root / "run"
            manifest = run_budget_engine(
                extract_out, out, CASE_ID, BUDGET_PROFILES["zero"], providers,
                protein_db=protein_db, search_backend=FakeBackend(),
            )
            plan = pd.read_csv(out / "tier2_plan.tsv", sep="\t", dtype=str).fillna("")
            planned = plan[plan["queries_planned"] == "1"]
            self.assertGreater(len(planned), 0)
            self.assertTrue((planned["query_status"] == "deferred_budget").all())
            self.assertEqual(0, int((plan["query_status"] == "executed").sum()))
            self.assertEqual("completed", manifest["status"])

    def test_positive_budget_start_time_overrun_is_budget_exhausted(self) -> None:
        """Fix 2: a run that starts over a positive declared peak-memory
        limit -- so every planned query is deferred before any executes --
        must be budget_exhausted, not completed. Previously the status rule
        only looked at queries_executed > 0 and mislabelled this
        completed."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out, providers, protein_db = _engine_fixture(root, self)
            limits = {"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 50, "wall_time_limit_seconds": 3600, "peak_memory_limit_mb": 100}
            out = root / "run"
            with mock.patch("candidate_qtl.tier2_planner._current_rss_mb", return_value=(99999.0, None)):
                manifest = run_budget_engine(
                    extract_out, out, CASE_ID, limits, providers,
                    protein_db=protein_db, search_backend=FakeBackend(),
                )
            plan = pd.read_csv(out / "tier2_plan.tsv", sep="\t", dtype=str).fillna("")
            planned = plan[plan["queries_planned"] == "1"]
            self.assertGreater(len(planned), 0)
            self.assertEqual(0, int((plan["query_status"] == "executed").sum()))
            self.assertGreater(int((plan["query_status"] == "deferred_budget").sum()), 0)
            self.assertEqual("budget_exhausted", manifest["status"])
            self.assertTrue(manifest["budgets"]["peak_memory_limit_mb"]["exhausted"])

    def test_peak_memory_ledger_and_resource_report_agree(self) -> None:
        """Fix 1 (completion): the peak in budget_manifest.json and the peak in
        benchmark_resource_usage.json must be the same quantity. Before this,
        the report was overwritten with ru_maxrss while the ledger held the
        observed RSS, so a run stopped for memory could report a peak
        comfortably *under* its own limit."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out, providers, protein_db = _engine_fixture(root, self)
            limits = {"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 50, "wall_time_limit_seconds": 3600, "peak_memory_limit_mb": 100}
            out = root / "run"
            with mock.patch("candidate_qtl.tier2_planner._current_rss_mb", return_value=(9999.0, None)):
                manifest = run_budget_engine(
                    extract_out, out, CASE_ID, limits, providers,
                    protein_db=protein_db, search_backend=FakeBackend(),
                )
            report = json.loads((out / "benchmark_resource_usage.json").read_text(encoding="utf-8"))
            reported_peak = report["usage"]["peak_memory_mb"]
            budget = manifest["budgets"]["peak_memory_limit_mb"]
            # the report shows the observed maximum, not a smaller process
            # high-water mark that would contradict the stop
            self.assertEqual(9999.0, reported_peak)
            self.assertEqual(budget["declared"] - budget["remaining"], min(reported_peak, budget["declared"]))
            self.assertTrue(budget["exhausted"])
            self.assertGreaterEqual(budget["remaining"], 0)
            # and the report's own status agrees with the manifest's
            self.assertEqual(manifest["status"], report["status"])

    def test_pre_run_process_high_water_mark_never_exhausts_this_run(self) -> None:
        """A lifetime ru_maxrss left by an earlier operation in a long-lived
        process is not this run's peak: with every in-run RSS sample under the
        limit, a huge pre-existing high-water mark must not mark the peak
        budget exhausted."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out, providers, protein_db = _engine_fixture(root, self)
            limits = {"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 50, "wall_time_limit_seconds": 3600, "peak_memory_limit_mb": 4096}
            out = root / "run"
            # ru_maxrss reports 50 GB and never grows during the run (the mark
            # predates it); every in-run RSS sample is a modest 120 MB.
            with mock.patch("candidate_qtl.budget_engine._peak_memory_mb", return_value=(51200.0, None)), \
                 mock.patch("candidate_qtl.tier2_planner._current_rss_mb", return_value=(120.0, None)):
                manifest = run_budget_engine(
                    extract_out, out, CASE_ID, limits, providers,
                    protein_db=protein_db, search_backend=FakeBackend(),
                )
            budget = manifest["budgets"]["peak_memory_limit_mb"]
            self.assertFalse(budget["exhausted"])
            self.assertEqual(120.0, budget["declared"] - budget["remaining"])
            self.assertEqual("completed", manifest["status"])
            report = json.loads((out / "benchmark_resource_usage.json").read_text(encoding="utf-8"))
            self.assertEqual(120.0, report["usage"]["peak_memory_mb"])

    def test_high_water_mark_raised_by_this_run_does_reach_the_ledger(self) -> None:
        """The other side of the same rule: when the run itself raises the
        process high-water mark, that growth IS this run's peak and is folded
        in."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out, providers, protein_db = _engine_fixture(root, self)
            limits = {"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 50, "wall_time_limit_seconds": 3600, "peak_memory_limit_mb": 4096}
            out = root / "run"
            marks = iter([100.0, 900.0])  # baseline at entry, then after the run
            with mock.patch("candidate_qtl.budget_engine._peak_memory_mb", side_effect=lambda: (next(marks), None)), \
                 mock.patch("candidate_qtl.tier2_planner._current_rss_mb", return_value=(120.0, None)):
                manifest = run_budget_engine(
                    extract_out, out, CASE_ID, limits, providers,
                    protein_db=protein_db, search_backend=FakeBackend(),
                )
            budget = manifest["budgets"]["peak_memory_limit_mb"]
            self.assertEqual(900.0, budget["declared"] - budget["remaining"])
            report = json.loads((out / "benchmark_resource_usage.json").read_text(encoding="utf-8"))
            self.assertEqual(900.0, report["usage"]["peak_memory_mb"])

    def test_peak_memory_report_keeps_process_high_water_mark_without_a_limit(self) -> None:
        """Folding ru_maxrss into the ledger must not lose it: with no declared
        peak limit (nothing samples RSS during the run) the report still
        carries a real measurement rather than 0."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out, providers, protein_db = _engine_fixture(root, self)
            limits = {"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 50, "wall_time_limit_seconds": 3600, "peak_memory_limit_mb": None}
            out = root / "run"
            run_budget_engine(
                extract_out, out, CASE_ID, limits, providers,
                protein_db=protein_db, search_backend=FakeBackend(),
            )
            report = json.loads((out / "benchmark_resource_usage.json").read_text(encoding="utf-8"))
            self.assertGreater(report["usage"]["peak_memory_mb"], 0)

    def test_provider_failed_handoff_not_evaluated_even_under_rescue_allocation(self) -> None:
        """Fix 4: a provider-failed row is always review_state=not_evaluated
        / handoff_note=provider_timeout, checked before the allocation
        branches -- even a rescue-escalated row (which would otherwise be
        searched and reviewed) must never be mislabelled
        rescue_pending_review, which would hide the provider failure behind
        ordinary rescue triage."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out, providers, protein_db = _engine_fixture(root, self)
            # G9 lands in a rescue cluster in this fixture (a provider-gap
            # singleton) and has a valid translated protein; sv evidence
            # escalates it (rescue_escalated="yes"), so it is actually
            # planned and searched instead of staying not_queried.
            (root / "sv.tsv").write_text("gene_id\tsv_type\nG9\tDEL\n", encoding="utf-8")
            providers = dict(providers)
            providers["sv"] = load_provider_table("sv", str(root / "sv.tsv"), "sv-db", "v1")
            out = root / "run"
            backend = FakeBackend(hit_gene_ids=set())
            original_search = backend.search

            def failing_search(query_fasta, db, bin, timeout=None):
                if Path(query_fasta).stem == "G9":
                    raise ExternalDependencyError("timeout", "fake provider timeout")
                return original_search(query_fasta, db, bin, timeout=timeout)

            backend.search = failing_search
            manifest = run_budget_engine(
                extract_out, out, CASE_ID, BUDGET_PROFILES["student_cpu"], providers,
                protein_db=protein_db, search_backend=backend,
            )
            plan = pd.read_csv(out / "tier2_plan.tsv", sep="\t", dtype=str).fillna("")
            g9 = plan[plan["gene_id"] == "G9"].iloc[0]
            self.assertEqual("rescue_escalation", g9["expansion_condition"])
            self.assertEqual("provider_failed", g9["query_status"])
            handoff = pd.read_csv(out / "budget_handoff.tsv", sep="\t", dtype=str).fillna("")
            g9_handoff = handoff[handoff["gene_id"] == "G9"].iloc[0]
            self.assertEqual("rescue", g9_handoff["allocation"])
            self.assertEqual("not_evaluated", g9_handoff["review_state"])
            self.assertEqual("provider_timeout", g9_handoff["handoff_note"])
            self.assertNotEqual("rescue_pending_review", g9_handoff["review_state"])
            # the provider failure stays unreviewed: it is never counted as a
            # reviewed candidate and never becomes a biological no_hit.
            reviewed = int((handoff["review_state"] == "priority_review").sum())
            self.assertEqual(len(handoff) - reviewed, manifest["unprocessed"]["unreviewed_candidate_count"])
            self.assertNotIn(g9["search_result_status"], {"hit", "no_hit"})


def _write_script(root: Path, name: str, body: str) -> Path:
    script = root / name
    script.write_text("#!/bin/sh\n" + body, encoding="utf-8")
    script.chmod(0o755)
    return script


class BudgetEngineRemediationTests(unittest.TestCase):
    def test_blastp_backend_probes_actual_bin(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            script = _write_script(root, "fake_blastp", 'if [ "$1" = "-version" ]; then echo "blastp 2.14.0+"; exit 0; fi\nsleep 60\n')
            backend = BlastpBackend(cache=None, bin=str(script))
            self.assertTrue(backend.available)
            self.assertEqual(os.path.realpath(str(script)), backend.resolved_bin)

    def test_blastp_backend_custom_bin_works_without_blastp_in_path(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            script = _write_script(root, "fake_blastp", 'if [ "$1" = "-version" ]; then echo "blastp 2.14.0+"; exit 0; fi\nsleep 60\n')
            with mock.patch("candidate_qtl.budget_engine.shutil.which", return_value=None):
                backend = BlastpBackend(cache=None, bin=str(script))
            self.assertTrue(backend.available)
            self.assertEqual(os.path.realpath(str(script)), backend.resolved_bin)

    def test_blastp_backend_missing_tool_unavailable(self) -> None:
        with tempfile.TemporaryDirectory():
            backend = BlastpBackend(cache=None, bin="definitely-not-a-real-blastp-binary-zzz")
            self.assertFalse(backend.available)
            self.assertIsNone(backend.resolved_bin)
            self.assertEqual("unknown", backend.version)
            with self.assertRaises(RuntimeError):
                backend.search(Path("q.fa"), "db", "blastp")

    def test_blastp_backend_captures_version(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            script = _write_script(root, "ver_blastp", 'echo "blastp 2.14.0+"; exit 0\n')
            backend = BlastpBackend(cache=None, bin=str(script))
            self.assertTrue(backend.available)
            self.assertEqual("blastp 2.14.0+", backend.version)

    def test_blastp_backend_nonzero_version_exit_unknown(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            script = _write_script(root, "bad_blastp", 'echo "blastp 2.14.0+"; exit 1\n')
            backend = BlastpBackend(cache=None, bin=str(script))
            self.assertEqual("unknown", backend.version)

    def test_blastp_backend_search_timeout_contract(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            hang = _write_script(root, "hang_blastp", 'if [ "$1" = "-version" ]; then echo "fake 1.0"; exit 0; fi\nsleep 60\n')
            backend = BlastpBackend(cache=None, bin=str(hang))
            query = root / "q.fa"
            query.write_text(">q\nMK\n", encoding="utf-8")
            start = time.monotonic()
            with self.assertRaises(ExternalDependencyError) as raised:
                backend.search(query, str(root / "db"), str(hang), timeout=0.5)
            self.assertEqual("timeout", raised.exception.kind)
            self.assertLess(time.monotonic() - start, 5.0)

    def test_write_json_rejects_nan(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "nan.json"
            with self.assertRaises(ValueError):
                _write_json(path, {"value": float("nan")})
            self.assertFalse(path.exists())


if __name__ == "__main__":
    unittest.main()
