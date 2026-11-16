from __future__ import annotations

import json
import os
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from candidate_qtl.budget_account import BUDGET_PROFILES, BudgetAccount
from candidate_qtl.budget_cache import ContentCache
from candidate_qtl.budget_engine import BlastpBackend
from candidate_qtl.cli import main
from candidate_qtl.core import ExternalDependencyError
from candidate_qtl.tier0_inventory import TIER0_COLUMNS, build_tier0_inventory, load_provider_table
from candidate_qtl.tier1_cluster import build_tier1_clusters
from candidate_qtl.tier2_planner import (
    BULK_MAX_GENES,
    PER_QUERY_COST,
    PLAN_MODES,
    QUERY_ROLES,
    QUERY_STATUSES,
    TIER2_COLUMNS,
    db_fingerprint,
    plan_sequence_search,
)

CASE_ID = "tier2-case"


class FakeBackend:
    def __init__(self, hits=None, cpu=0.01, available=True, resolved_bin="blastp", version="1.0", sleep_seconds=0.0):
        self.available = available
        self.cache = None
        self.cpu = cpu
        self.hits = hits if hits is not None else [{"qseqid": "G1", "sseqid": "hit1", "pident": "95", "evalue": "1e-50", "bitscore": "200"}]
        self.resolved_bin = resolved_bin
        self.version = version
        self.sleep_seconds = sleep_seconds
        self.search_calls = 0
        self.seen_timeouts = []

    def search(self, query_fasta, db, bin, timeout=None):
        self.search_calls += 1
        self.seen_timeouts.append(timeout)
        if self.sleep_seconds:
            time.sleep(self.sleep_seconds)
        return self.hits, self.cpu


class TimeoutBackend:
    """Backend whose search always hangs (sleep) then raises a provider timeout."""

    def __init__(self, sleep_seconds: float = 0.2, kind: str = "timeout"):
        self.available = True
        self.cache = None
        self.kind = kind
        self.search_calls = 0
        self.seen_timeouts = []
        self.sleep_seconds = sleep_seconds

    def search(self, query_fasta, db, bin, timeout=None):
        self.search_calls += 1
        self.seen_timeouts.append(timeout)
        if self.sleep_seconds:
            time.sleep(self.sleep_seconds)
        raise ExternalDependencyError(self.kind, "fake provider timeout")


def _tier2_genes():
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


def _run_extract(root: Path, testcase: unittest.TestCase) -> Path:
    fasta, gff = _write_extract_fixture(root, _tier2_genes())
    extract_out = root / "extract"
    testcase.assertEqual(
        0,
        main(["extract", "--fasta", str(fasta), "--gff3", str(gff), "--interval", "chr1:1-1000", "--out-dir", str(extract_out)]),
    )
    return extract_out


def _write_provider_tables(root: Path) -> dict[str, dict[str, object]]:
    (root / "expression.tsv").write_text("gene_id\tvalue\tsource_id\nG5\t10\tE1\nG5\t20\tE1\nG2\t5\tE2\n", encoding="utf-8")
    (root / "domains.tsv").write_text("gene_id\tdomain_id\nG5\tPF001\nG5\tPF002\n", encoding="utf-8")
    (root / "orthogroup.tsv").write_text("gene_id\torthogroup\nG5\tOG0001\n", encoding="utf-8")
    (root / "variant.tsv").write_text("gene_id\teffect\tsource_id\nG5\tmissense\tV1\n", encoding="utf-8")
    (root / "family.tsv").write_text("gene_id\tfamily_id\nG5\tFAM001\n", encoding="utf-8")
    return {
        "expression": load_provider_table("expression", str(root / "expression.tsv"), "expr-db", "v2"),
        "domains": load_provider_table("domains", str(root / "domains.tsv"), "dom-db", "v1"),
        "orthogroup": load_provider_table("orthogroup", str(root / "orthogroup.tsv"), "egg-db", "v1"),
        "variant": load_provider_table("variant", str(root / "variant.tsv"), "var-db", "v1"),
        "family": load_provider_table("family", str(root / "family.tsv"), "fam-db", "v1"),
    }


def _tier0_tier1(root: Path, testcase: unittest.TestCase, providers: dict[str, dict[str, object]] | None = None) -> tuple[list[dict[str, str]], list[dict[str, str]], Path]:
    extract_out = _run_extract(root, testcase)
    if providers is None:
        providers = _write_provider_tables(root)
    tier0_out = root / "tier0"
    build_tier0_inventory(extract_out, tier0_out, providers, CASE_ID)
    inventory_rows = pd.read_csv(tier0_out / "tier0_inventory.tsv", sep="\t", dtype=str).fillna("").to_dict("records")
    tier1_out = root / "tier1"
    build_tier1_clusters(inventory_rows, tier1_out, CASE_ID, providers=providers)
    cluster_rows = pd.read_csv(tier1_out / "tier1_clusters.tsv", sep="\t", dtype=str).fillna("").to_dict("records")
    return inventory_rows, cluster_rows, extract_out


def _hand_row(gene_id: str, **overrides) -> dict[str, str]:
    row = {column: "" for column in TIER0_COLUMNS}
    row.update(
        {
            "gene_id": gene_id,
            "model_status": "complete_orf",
            "annotation_status": "present",
            "symbol_status": "present",
            "provider_applicability": "not_evaluated",
        }
    )
    for cell in ("domain_status", "orthogroup_status", "expression_status", "variant_status", "sv_status", "pav_status", "promoter_status"):
        row[cell] = "not_evaluated"
    row.update(overrides)
    return row


def _hand_cluster(cluster_id: str, allocation: str, members: list[str], representative: str, rescue_bucket: str = "not_applicable") -> dict[str, str]:    return {
        "cluster_id": cluster_id,
        "seed_reason": "hand_built",
        "representative_gene_id": representative,
        "representative_reason": "hand_built",
        "member_gene_ids": ";".join(members),
        "member_count": str(len(members)),
        "rule_key": "annotation",
        "cluster_key": cluster_id,
        "family_ids": "",
        "domain_ids": "",
        "orthogroup_ids": "",
        "allocation": allocation,
        "allocation_reason": "hand_built",
        "rescue_bucket": rescue_bucket,
        "dependency_cluster_ids": "",
    }


def _plan_rows(out: Path) -> dict[str, dict[str, str]]:
    frame = pd.read_csv(out / "tier2_plan.tsv", sep="\t", dtype=str).fillna("")
    return {row["gene_id"]: dict(row) for row in frame.to_dict("records")}


def _all_proteins() -> dict[str, str]:
    return {f"G{index}": "MKA" for index in range(1, 10)}


class Tier2PlannerTests(unittest.TestCase):
    def test_priority_representatives_first_and_expansion_triggers(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inventory = [
                _hand_row("G1", gene_name="GeneA", conflict_flags="expression:G1:value"),
                _hand_row("G6", gene_name="GeneA"),
                _hand_row("G8", gene_name="GeneA"),
                _hand_row("G9", gene_name="GeneB"),
            ]
            clusters = [
                _hand_cluster("cls_a", "priority", ["G1", "G6", "G8"], "G6"),
                _hand_cluster("cls_b", "priority", ["G9"], "G9"),
            ]
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            account = BudgetAccount(BUDGET_PROFILES["student_cpu"])
            out = root / "plan"
            plan_sequence_search(inventory, clusters, account, out, CASE_ID, protein_db=root / "db.fa", search_backend=FakeBackend())
            rows = _plan_rows(out)
            rep = rows["G9"]
            self.assertEqual("1", rep["plan_step"])
            self.assertEqual("representative", rep["query_role"])
            self.assertEqual("representative_first", rep["expansion_condition"])
            self.assertEqual("1", rep["queries_planned"])
            self.assertEqual("planned", rep["query_status"])
            ambiguous = rows["G1"]
            self.assertEqual("2", ambiguous["plan_step"])
            self.assertEqual("expansion_member", ambiguous["query_role"])
            self.assertEqual("ambiguity_status", ambiguous["expansion_condition"])
            self.assertEqual("1", ambiguous["queries_planned"])
            self.assertEqual("ambiguous_annotation", ambiguous["ambiguity_status"])
            untriggered = rows["G8"]
            self.assertEqual("2", untriggered["plan_step"])
            self.assertEqual("not_queried", untriggered["query_role"])
            self.assertEqual("0", untriggered["queries_planned"])
            self.assertEqual("not_planned", untriggered["query_status"])
            self.assertEqual("expansion_not_triggered", untriggered["deferred_reason"])
            self.assertEqual("not_ambiguous", untriggered["ambiguity_status"])
            self.assertEqual("no", untriggered["rescue_escalated"])
            self.assertEqual("no", untriggered["user_flagged"])

    def test_user_flagged_triggers_expansion(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inventory = [_hand_row("G1", gene_name="GeneA"), _hand_row("G6", gene_name="GeneA"), _hand_row("G8", gene_name="GeneA")]
            clusters = [_hand_cluster("cls_a", "priority", ["G1", "G6", "G8"], "G6")]
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            account = BudgetAccount(BUDGET_PROFILES["student_cpu"])
            out = root / "plan"
            plan_sequence_search(inventory, clusters, account, out, CASE_ID, user_flagged={"G1"}, protein_db=root / "db.fa", search_backend=FakeBackend())
            rows = _plan_rows(out)
            flagged = rows["G1"]
            self.assertEqual("yes", flagged["user_flagged"])
            self.assertEqual("user_flagged", flagged["expansion_condition"])
            self.assertEqual("1", flagged["queries_planned"])
            self.assertEqual("expansion_member", flagged["query_role"])

    def test_rescue_escalation_plans_only_escalated_genes(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inventory = [
                _hand_row("R1", conflict_flags="sv:R1:sv_type"),
                _hand_row("R2"),
                _hand_row("R3", sv_evidence="DEL", sv_status="present"),
            ]
            clusters = [
                _hand_cluster("cls_r1", "rescue", ["R1"], "R1", rescue_bucket="unclassified"),
                _hand_cluster("cls_r2", "rescue", ["R2"], "R2", rescue_bucket="unclassified"),
                _hand_cluster("cls_r3", "rescue", ["R3"], "R3", rescue_bucket="sv_pav_regulatory"),
            ]
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            account = BudgetAccount(BUDGET_PROFILES["student_cpu"])
            out = root / "plan"
            plan_sequence_search(inventory, clusters, account, out, CASE_ID, protein_db=root / "db.fa", search_backend=FakeBackend())
            rows = _plan_rows(out)
            escalated = rows["R1"]
            self.assertEqual("yes", escalated["rescue_escalated"])
            self.assertEqual("4", escalated["plan_step"])
            self.assertEqual("representative", escalated["query_role"])
            self.assertEqual("rescue_escalation", escalated["expansion_condition"])
            self.assertEqual("1", escalated["queries_planned"])
            sv_escalated = rows["R3"]
            self.assertEqual("yes", sv_escalated["rescue_escalated"])
            self.assertEqual("1", sv_escalated["queries_planned"])
            pending = rows["R2"]
            self.assertEqual("no", pending["rescue_escalated"])
            self.assertEqual("not_planned", pending["query_status"])
            self.assertEqual("0", pending["queries_planned"])
            self.assertEqual("rescue_pending_review", pending["deferred_reason"])
            self.assertEqual("not_queried", pending["query_role"])

    def test_ambiguous_homology_escalates_rescue(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inventory = [_hand_row("H1", orthogroup="OG1;OG2", orthogroup_status="present"), _hand_row("H2")]
            clusters = [_hand_cluster("cls_h", "rescue", ["H1", "H2"], "H1", rescue_bucket="unclassified")]
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            account = BudgetAccount(BUDGET_PROFILES["student_cpu"])
            out = root / "plan"
            plan_sequence_search(inventory, clusters, account, out, CASE_ID, protein_db=root / "db.fa", search_backend=FakeBackend())
            rows = _plan_rows(out)
            self.assertEqual("ambiguous_homology", rows["H1"]["ambiguity_status"])
            self.assertEqual("yes", rows["H1"]["rescue_escalated"])
            self.assertEqual("1", rows["H1"]["queries_planned"])

    def test_defer_representatives_planned_after_priority(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inventory = [_hand_row("D1"), _hand_row("D2")]
            clusters = [_hand_cluster("cls_d", "defer", ["D1", "D2"], "D1")]
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            account = BudgetAccount(BUDGET_PROFILES["student_cpu"])
            out = root / "plan"
            plan_sequence_search(inventory, clusters, account, out, CASE_ID, protein_db=root / "db.fa", search_backend=FakeBackend())
            rows = _plan_rows(out)
            self.assertEqual("3", rows["D1"]["plan_step"])
            self.assertEqual("representative", rows["D1"]["query_role"])
            self.assertEqual("1", rows["D1"]["queries_planned"])
            self.assertEqual("5", rows["D2"]["plan_step"])
            self.assertEqual("not_selected", rows["D2"]["deferred_reason"])
            self.assertEqual("0", rows["D2"]["queries_planned"])

    def test_a_b_discordance_always_not_evaluated(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inventory_rows, cluster_rows, _extract = _tier0_tier1(root, self)
            account = BudgetAccount(BUDGET_PROFILES["student_cpu"])
            out = root / "plan"
            plan_sequence_search(inventory_rows, cluster_rows, account, out, CASE_ID, search_backend=FakeBackend())
            rows = _plan_rows(out)
            self.assertEqual(9, len(rows))
            for row in rows.values():
                self.assertEqual("not_evaluated", row["a_b_discordance"])

    def test_nine_gene_fixture_plan_shapes(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inventory_rows, cluster_rows, _extract = _tier0_tier1(root, self)
            account = BudgetAccount(BUDGET_PROFILES["student_cpu"])
            out = root / "plan"
            plan_sequence_search(inventory_rows, cluster_rows, account, out, CASE_ID, search_backend=FakeBackend())
            rows = _plan_rows(out)
            step1 = {gene for gene, row in rows.items() if row["plan_step"] == "1"}
            # tier1's unique-name rule (audit finding G1) drops G9's singleton
            # GeneB annotation to a provider-gap rescue, so step 1 holds G5/G6
            self.assertEqual({"G5", "G6"}, step1)
            self.assertTrue(all(rows[gene]["query_role"] == "representative" for gene in step1))
            for gene in ("G1", "G7", "G8"):
                self.assertEqual("not_queried", rows[gene]["query_role"])
                self.assertEqual("expansion_not_triggered", rows[gene]["deferred_reason"])
            self.assertEqual("3", rows["G2"]["plan_step"])
            self.assertEqual("representative", rows["G2"]["query_role"])
            self.assertEqual("1", rows["G2"]["queries_planned"])
            for gene in ("G3", "G4", "G9"):
                self.assertEqual("rescue_pending_review", rows[gene]["deferred_reason"])
                self.assertEqual("not_planned", rows[gene]["query_status"])
            self.assertEqual("hierarchical", json.loads((out / "tier2_manifest.json").read_text(encoding="utf-8"))["plan_mode"])

    def test_bulk_selected_when_hierarchical_plans_every_gene(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inventory = [
                _hand_row("G1", gene_name="GeneA", conflict_flags="expression:G1:value"),
                _hand_row("G6", gene_name="GeneA"),
                _hand_row("G9", gene_name="GeneB"),
            ]
            clusters = [
                _hand_cluster("cls_a", "priority", ["G1", "G6"], "G6"),
                _hand_cluster("cls_b", "priority", ["G9"], "G9"),
            ]
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            account = BudgetAccount(BUDGET_PROFILES["student_cpu"])
            out = root / "plan"
            manifest = plan_sequence_search(inventory, clusters, account, out, CASE_ID, protein_db=root / "db.fa", search_backend=FakeBackend())
            rows = _plan_rows(out)
            self.assertEqual(3, len(rows))
            self.assertTrue(manifest["bulk_selected"])
            self.assertEqual("bulk", manifest["plan_mode"])
            self.assertIn("bulk", PLAN_MODES)
            for row in rows.values():
                self.assertEqual("representative", row["query_role"])
                self.assertEqual("representative_first", row["expansion_condition"])
                self.assertEqual("1", row["queries_planned"])
                self.assertEqual("planned", row["query_status"])
            self.assertLessEqual(3, BULK_MAX_GENES)
            self.assertEqual(3 * PER_QUERY_COST, manifest["cost_estimate_hierarchical"])
            self.assertEqual(3 * PER_QUERY_COST, manifest["cost_estimate_bulk"])

    def test_nine_gene_fixture_not_bulk(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inventory_rows, cluster_rows, _extract = _tier0_tier1(root, self)
            account = BudgetAccount(BUDGET_PROFILES["student_cpu"])
            out = root / "plan"
            manifest = plan_sequence_search(inventory_rows, cluster_rows, account, out, CASE_ID, search_backend=FakeBackend())
            self.assertFalse(manifest["bulk_selected"])
            self.assertEqual("hierarchical", manifest["plan_mode"])
            self.assertEqual(3, manifest["queries_planned"])

    def test_ambiguous_homology_precedence_under_annotation(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inventory = [_hand_row("P1", orthogroup="OG1;OG2", conflict_flags="")]
            clusters = [_hand_cluster("cls_p", "rescue", ["P1"], "P1", rescue_bucket="unclassified")]
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            account = BudgetAccount(BUDGET_PROFILES["student_cpu"])
            out = root / "plan"
            plan_sequence_search(inventory, clusters, account, out, CASE_ID, protein_db=root / "db.fa", search_backend=FakeBackend())
            rows = _plan_rows(out)
            self.assertEqual("ambiguous_homology", rows["P1"]["ambiguity_status"])
            self.assertEqual("yes", rows["P1"]["rescue_escalated"])

    def test_zero_budget_plans_everything_then_defers_all(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inventory_rows, cluster_rows, _extract = _tier0_tier1(root, self)
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            account = BudgetAccount(BUDGET_PROFILES["zero"])
            out = root / "plan"
            manifest = plan_sequence_search(inventory_rows, cluster_rows, account, out, CASE_ID, protein_db=root / "db.fa", search_backend=FakeBackend(), proteins=_all_proteins())
            rows = _plan_rows(out)
            self.assertGreater(manifest["queries_planned"], 0)
            self.assertEqual(0, manifest["queries_executed"])
            deferred = [row for row in rows.values() if row["query_status"] == "deferred_budget"]
            self.assertEqual(3, len(deferred))
            for row in deferred:
                self.assertEqual("budget_exhausted", row["deferred_reason"])
                self.assertEqual("not_run", row["search_result_status"])
            self.assertEqual(0, account.freeze()["usage"]["sequence_query_count"])
            self.assertEqual(1, len(manifest["truncation_points"]))
            self.assertEqual("sequence_query", manifest["truncation_points"][0]["kind"])
            self.assertEqual(0, account.freeze()["usage"]["searched_protein_count"])

    def test_provider_missing_never_consumes(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inventory_rows, cluster_rows, _extract = _tier0_tier1(root, self)
            account = BudgetAccount(BUDGET_PROFILES["student_cpu"])
            out = root / "plan"
            manifest = plan_sequence_search(inventory_rows, cluster_rows, account, out, CASE_ID, search_backend=FakeBackend(available=False))
            rows = _plan_rows(out)
            self.assertEqual("not_evaluated_provider_missing", manifest["search_provider"])
            planned = [row for row in rows.values() if row["queries_planned"] == "1"]
            self.assertEqual(3, len(planned))
            for row in planned:
                self.assertEqual("not_evaluated", row["search_result_status"])
            self.assertEqual(0, manifest["queries_executed"])
            self.assertEqual(0, account.freeze()["usage"]["sequence_query_count"])
            self.assertEqual(0, account.freeze()["usage"]["cache_miss_count"])

    def test_manifest_written_without_truncation_points_and_no_timestamps(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inventory_rows, cluster_rows, _extract = _tier0_tier1(root, self)
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            account = BudgetAccount(BUDGET_PROFILES["student_cpu"])
            out = root / "plan"
            returned = plan_sequence_search(inventory_rows, cluster_rows, account, out, CASE_ID, protein_db=root / "db.fa", search_backend=FakeBackend(), proteins=_all_proteins())
            written = json.loads((out / "tier2_manifest.json").read_text(encoding="utf-8"))
            # returned carries provider_errors, written only the deterministic
            # provider_error_count (which the returned dict does not repeat)
            returned_without_errors = {key: value for key, value in returned.items() if key != "provider_errors"}
            written_without_count = {key: value for key, value in written.items() if key != "provider_error_count"}
            self.assertEqual(returned_without_errors, {**written_without_count, "truncation_points": returned["truncation_points"]})
            self.assertEqual([], returned["provider_errors"])
            self.assertEqual(0, written["provider_error_count"])
            self.assertIn("truncation_points", returned)
            self.assertNotIn("truncation_points", written)
            self.assertNotIn("provider_errors", written)
            self.assertIn("peak_memory_enforcement", written)
            self.assertEqual("1.0", written["schema_version"])
            self.assertEqual(CASE_ID, written["case_id"])
            self.assertEqual(9, written["gene_count"])
            self.assertEqual(6, written["cluster_count"])
            self.assertEqual(3, written["queries_planned"])
            self.assertEqual(3, written["queries_executed"])
            self.assertEqual(0, written["queries_deferred_budget"])
            self.assertEqual(0, written["queries_deferred_plan"])
            self.assertEqual(6, written["queries_not_planned"])
            self.assertEqual(0, written["queries_provider_failed"])
            self.assertEqual("blastp", written["search_provider"])
            self.assertEqual(3 * PER_QUERY_COST, written["cost_estimate_hierarchical"])
            self.assertEqual(9 * PER_QUERY_COST, written["cost_estimate_bulk"])
            self.assertFalse(written["bulk_selected"])
            self.assertEqual({"a_b_discordance": 0, "ambiguity_status": {"ambiguous_annotation": 1, "not_ambiguous": 8}, "rescue_escalation": 0, "user_flagged": 0}, written["expansion_counts"])
            self.assertNotIn("timestamp", str(written))
            for row in _plan_rows(out).values():
                self.assertIn(row["query_role"], QUERY_ROLES)
                self.assertIn(row["query_status"], QUERY_STATUSES)
            self.assertEqual(TIER2_COLUMNS, tuple(pd.read_csv(out / "tier2_plan.tsv", sep="\t", dtype=str).columns))

    def test_cache_hit_recharges_stored_cpu_seconds(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inventory_rows, cluster_rows, _extract = _tier0_tier1(root, self)
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            backend = FakeBackend(cpu=2.5)
            backend.cache = ContentCache(root / "cache")
            cold = BudgetAccount(BUDGET_PROFILES["student_cpu"])
            plan_sequence_search(inventory_rows, cluster_rows, cold, root / "cold", CASE_ID, protein_db=root / "db.fa", search_backend=backend, proteins=_all_proteins())
            warm = BudgetAccount(BUDGET_PROFILES["student_cpu"])
            plan_sequence_search(inventory_rows, cluster_rows, warm, root / "warm", CASE_ID, protein_db=root / "db.fa", search_backend=backend, proteins=_all_proteins())
            cold_usage = cold.freeze()["usage"]
            warm_usage = warm.freeze()["usage"]
            # warm rerun recharges the stored measured cpu cost from the cache
            self.assertEqual(3 * 2.5, cold_usage["sequence_search_cpu_seconds"])
            self.assertEqual(cold_usage["sequence_search_cpu_seconds"], warm_usage["sequence_search_cpu_seconds"])
            self.assertEqual(3, cold_usage["sequence_query_count"])
            self.assertEqual(0, warm_usage["sequence_query_count"])
            self.assertEqual(3, cold_usage["searched_protein_count"])
            self.assertEqual(0, warm_usage["searched_protein_count"])
            self.assertEqual(3, warm_usage["cache_hit_count"])
            self.assertEqual(0, warm_usage["cache_miss_count"])

    def test_exhaustive_override_plans_every_gene(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inventory_rows, cluster_rows, _extract = _tier0_tier1(root, self)
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            account = BudgetAccount(BUDGET_PROFILES["exhaustive"])
            out = root / "plan"
            manifest = plan_sequence_search(
                inventory_rows, cluster_rows, account, out, CASE_ID,
                protein_db=root / "db.fa", search_backend=FakeBackend(), exhaustive_override=True, proteins=_all_proteins(),
            )
            rows = _plan_rows(out)
            self.assertEqual("exhaustive", manifest["plan_mode"])
            self.assertEqual(9, manifest["queries_planned"])
            self.assertEqual(9, manifest["queries_executed"])
            self.assertFalse(manifest["bulk_selected"])
            for row in rows.values():
                self.assertEqual("representative", row["query_role"])
                self.assertEqual("representative_first", row["expansion_condition"])
                self.assertEqual("1", row["queries_planned"])
                self.assertEqual("executed", row["query_status"])
            self.assertEqual(9, account.freeze()["usage"]["sequence_query_count"])


def _simple_fixture(root: Path, gene_ids: list[str]) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, str]]:
    inventory = [_hand_row(gene_id, gene_name=f"name_{gene_id}") for gene_id in gene_ids]
    clusters = [_hand_cluster(f"cls_{index}", "priority", [gene_id], gene_id) for index, gene_id in enumerate(gene_ids)]
    return inventory, clusters, {gene_id: "MKA" for gene_id in gene_ids}


def _write_script(root: Path, name: str, body: str) -> Path:
    script = root / name
    script.write_text("#!/bin/sh\n" + body, encoding="utf-8")
    script.chmod(0o755)
    return script


class BudgetRemediationPlannerTests(unittest.TestCase):
    def test_db_fingerprint_full_content_detects_tail_differences(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            size = 2 * 1024 * 1024
            first_dir = root / "db_a"
            second_dir = root / "db_b"
            first_dir.mkdir()
            second_dir.mkdir()
            first = first_dir / "proteome.fa"
            second = second_dir / "proteome.fa"
            first.write_bytes(b"A" * (size - 8) + b"TAILAAAA")
            second.write_bytes(b"A" * (size - 8) + b"TAILBBBB")
            fp_a = db_fingerprint(first)
            fp_b = db_fingerprint(second)
            self.assertIsNotNone(fp_a)
            self.assertNotEqual(fp_a, fp_b)

    def test_db_fingerprint_blast_prefix_without_prefix_file(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            for name in ("db.pin", "db.phr", "db.psq"):
                (root / name).write_bytes(b"data-" + name.encode("utf-8"))
            self.assertIsNotNone(db_fingerprint(root / "db"))
            self.assertEqual(db_fingerprint(root / "db"), db_fingerprint(root / "db"))

    def test_db_fingerprint_missing_path_and_empty_dir_return_none(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.assertIsNone(db_fingerprint(root / "missing"))
            self.assertIsNone(db_fingerprint(None))
            empty = root / "empty"
            empty.mkdir()
            self.assertIsNone(db_fingerprint(empty))

    def test_db_fingerprint_order_invariant_but_name_sensitive(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            dir_a = root / "a"
            dir_b = root / "b"
            dir_a.mkdir()
            dir_b.mkdir()
            for name in ("one.fa", "two.fa"):
                (dir_a / name).write_bytes(b"content-" + name.encode("utf-8"))
            (dir_b / "two.fa").write_bytes(b"content-two.fa")
            (dir_b / "one.fa").write_bytes(b"content-one.fa")
            self.assertEqual(db_fingerprint(dir_a), db_fingerprint(dir_b))
            renamed = root / "c"
            renamed.mkdir()
            (renamed / "one.fa").write_bytes(b"content-one.fa")
            (renamed / "two2.fa").write_bytes(b"content-two.fa")
            self.assertNotEqual(db_fingerprint(dir_a), db_fingerprint(renamed))

    def test_db_missing_path_marks_rows_not_evaluated_no_backend_calls(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inventory, clusters, proteins = _simple_fixture(root, ["G1", "G2"])
            backend = FakeBackend()
            account = BudgetAccount(BUDGET_PROFILES["student_cpu"])
            out = root / "plan"
            manifest = plan_sequence_search(inventory, clusters, account, out, CASE_ID, protein_db=root / "missing-db", search_backend=backend, proteins=proteins)
            rows = _plan_rows(out)
            self.assertEqual(0, backend.search_calls)
            self.assertEqual(0, manifest["queries_executed"])
            self.assertEqual(0, manifest["queries_deferred_budget"])
            for row in rows.values():
                self.assertEqual("not_evaluated", row["search_result_status"])
            self.assertEqual(0, account.freeze()["usage"]["sequence_query_count"])
            self.assertEqual(0, account.freeze()["usage"]["cache_miss_count"])

    def test_provider_version_change_invalidates_warm_cache(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            inventory, clusters, proteins = _simple_fixture(root, ["G1"])
            cache_dir = root / "cache"
            cold_backend = FakeBackend(version="1.0")
            cold_backend.cache = ContentCache(cache_dir)
            plan_sequence_search(inventory, clusters, BudgetAccount(BUDGET_PROFILES["student_cpu"]), root / "cold", CASE_ID, protein_db=root / "db.fa", search_backend=cold_backend, proteins=proteins)
            self.assertEqual(1, cold_backend.search_calls)
            warm_backend = FakeBackend(version="2.0")
            warm_backend.cache = ContentCache(cache_dir)
            warm = BudgetAccount(BUDGET_PROFILES["student_cpu"])
            plan_sequence_search(inventory, clusters, warm, root / "warm", CASE_ID, protein_db=root / "db.fa", search_backend=warm_backend, proteins=proteins)
            # the version participates in the content key, so the "1.0" entry
            # is a plain miss (different key), not a stale entry, and the
            # search re-executes instead of reusing the warm hits
            self.assertEqual(1, warm_backend.search_calls)
            self.assertEqual(0, warm.freeze()["usage"]["cache_hit_count"])
            self.assertEqual(1, warm.freeze()["usage"]["cache_miss_count"])
            self.assertEqual(0, warm_backend.cache.stats()["stale_reject_count"])

    def test_stale_provider_version_metadata_rejected_and_recomputed(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            inventory, clusters, proteins = _simple_fixture(root, ["G1"])
            cache_dir = root / "cache"
            first = FakeBackend(version="1.0")
            first.cache = ContentCache(cache_dir)
            plan_sequence_search(inventory, clusters, BudgetAccount(BUDGET_PROFILES["student_cpu"]), root / "cold", CASE_ID, protein_db=root / "db.fa", search_backend=first, proteins=proteins)
            # simulate an entry whose recorded provider_version disagrees with
            # the key (e.g. written by a cache layout without the version):
            # same key file, bogus recorded version
            for meta in cache_dir.glob("*.json"):
                data = json.loads(meta.read_text(encoding="utf-8"))
                data["provider_version"] = "0.9"
                meta.write_text(json.dumps(data, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            second = FakeBackend(version="1.0")
            second.cache = ContentCache(cache_dir)
            warm = BudgetAccount(BUDGET_PROFILES["student_cpu"])
            plan_sequence_search(inventory, clusters, warm, root / "warm", CASE_ID, protein_db=root / "db.fa", search_backend=second, proteins=proteins)
            self.assertEqual(1, second.cache.stats()["stale_reject_count"])
            self.assertEqual(1, second.search_calls)
            self.assertEqual(1, warm.freeze()["usage"]["cache_miss_count"])
            self.assertEqual(1, warm.freeze()["usage"]["searched_protein_count"])

    def test_warm_cache_hits_bypass_zero_budget(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inventory_rows, cluster_rows, _extract = _tier0_tier1(root, self)
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            cache_dir = root / "cache"
            cold_backend = FakeBackend()
            cold_backend.cache = ContentCache(cache_dir)
            cold = BudgetAccount(BUDGET_PROFILES["student_cpu"])
            plan_sequence_search(inventory_rows, cluster_rows, cold, root / "cold", CASE_ID, protein_db=root / "db.fa", search_backend=cold_backend, proteins=_all_proteins())
            self.assertEqual(3, cold_backend.search_calls)
            cold_plan = _plan_rows(root / "cold")
            warm_backend = FakeBackend()
            warm_backend.cache = ContentCache(cache_dir)
            warm = BudgetAccount(BUDGET_PROFILES["zero"])
            manifest = plan_sequence_search(inventory_rows, cluster_rows, warm, root / "warm", CASE_ID, protein_db=root / "db.fa", search_backend=warm_backend, proteins=_all_proteins())
            rows = _plan_rows(root / "warm")
            executed = [gene for gene, row in rows.items() if row["query_status"] == "executed"]
            self.assertEqual(manifest["queries_planned"], len(executed))
            self.assertEqual(0, warm_backend.search_calls)
            self.assertEqual(0, warm.freeze()["usage"]["sequence_query_count"])
            self.assertEqual(3, warm.freeze()["usage"]["cache_hit_count"])
            for gene, row in rows.items():
                self.assertEqual(cold_plan[gene]["search_result_status"], row["search_result_status"])

    def test_cold_zero_budget_defers_all_without_backend_calls(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inventory_rows, cluster_rows, _extract = _tier0_tier1(root, self)
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            backend = FakeBackend()
            backend.cache = ContentCache(root / "cache")
            manifest = plan_sequence_search(inventory_rows, cluster_rows, BudgetAccount(BUDGET_PROFILES["zero"]), root / "plan", CASE_ID, protein_db=root / "db.fa", search_backend=backend, proteins=_all_proteins())
            rows = _plan_rows(root / "plan")
            deferred = [row for row in rows.values() if row["query_status"] == "deferred_budget"]
            self.assertEqual(3, len(deferred))
            self.assertEqual(0, backend.search_calls)
            self.assertEqual("sequence_query", manifest["truncation_points"][0]["kind"])

    def test_corrupt_payload_with_zero_budget_deferred_no_recompute(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inventory_rows, cluster_rows, _extract = _tier0_tier1(root, self)
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            cache_dir = root / "cache"
            cold_backend = FakeBackend()
            cold_backend.cache = ContentCache(cache_dir)
            plan_sequence_search(inventory_rows, cluster_rows, BudgetAccount(BUDGET_PROFILES["student_cpu"]), root / "cold", CASE_ID, protein_db=root / "db.fa", search_backend=cold_backend, proteins=_all_proteins())
            payloads = sorted(cache_dir.glob("*.payload"))
            self.assertGreater(len(payloads), 0)
            # corrupt every payload so the first planned row's lookup is the
            # only one before the zero-budget gate stops the loop
            for payload in payloads:
                data = bytearray(payload.read_bytes())
                data[-1] ^= 0xFF
                payload.write_bytes(bytes(data))
            warm_backend = FakeBackend()
            warm_backend.cache = ContentCache(cache_dir)
            plan_sequence_search(inventory_rows, cluster_rows, BudgetAccount(BUDGET_PROFILES["zero"]), root / "warm", CASE_ID, protein_db=root / "db.fa", search_backend=warm_backend, proteins=_all_proteins())
            self.assertEqual(1, warm_backend.cache.stats()["corrupt_reject_count"])
            rows = _plan_rows(root / "warm")
            planned = [row for row in rows.values() if row["queries_planned"] == "1"]
            self.assertTrue(all(row["query_status"] == "deferred_budget" for row in planned))
            self.assertEqual(0, warm_backend.search_calls)

    def test_stale_provider_version_with_zero_budget_deferred(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            inventory, clusters, proteins = _simple_fixture(root, ["G1"])
            cache_dir = root / "cache"
            first = FakeBackend(version="1.0")
            first.cache = ContentCache(cache_dir)
            plan_sequence_search(inventory, clusters, BudgetAccount(BUDGET_PROFILES["student_cpu"]), root / "cold", CASE_ID, protein_db=root / "db.fa", search_backend=first, proteins=proteins)
            for meta in cache_dir.glob("*.json"):
                data = json.loads(meta.read_text(encoding="utf-8"))
                data["provider_version"] = "0.9"
                meta.write_text(json.dumps(data, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            second = FakeBackend(version="1.0")
            second.cache = ContentCache(cache_dir)
            manifest = plan_sequence_search(inventory, clusters, BudgetAccount(BUDGET_PROFILES["zero"]), root / "warm", CASE_ID, protein_db=root / "db.fa", search_backend=second, proteins=proteins)
            self.assertEqual(1, second.cache.stats()["stale_reject_count"])
            rows = _plan_rows(root / "warm")
            self.assertEqual("deferred_budget", rows["G1"]["query_status"])
            self.assertEqual(0, second.search_calls)
            self.assertEqual("sequence_query", manifest["truncation_points"][0]["kind"])

    def test_wall_time_budget_stops_after_first_expensive_query(self) -> None:
        """Fix 3: the wall-time charge uses the measured elapsed time, not the
        backend's self-reported cpu_seconds -- a backend that (mis)reports a
        huge cpu_seconds (5.0s, versus a 1s wall limit) must not blow the
        budget instantly; the real ~0.4s sleeps are what get charged."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            inventory, clusters, proteins = _simple_fixture(root, ["G1", "G2", "G3", "G4"])
            limits = {"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 50, "wall_time_limit_seconds": 1, "peak_memory_limit_mb": None}
            account = BudgetAccount(limits)
            backend = FakeBackend(cpu=5.0, sleep_seconds=0.4)
            manifest = plan_sequence_search(inventory, clusters, account, root / "plan", CASE_ID, protein_db=root / "db.fa", search_backend=backend, proteins=proteins)
            rows = _plan_rows(root / "plan")
            executed = {gene for gene, row in rows.items() if row["query_status"] == "executed"}
            deferred = {gene for gene, row in rows.items() if row["query_status"] == "deferred_budget"}
            # G1 and G2 fit inside the 1s budget at ~0.4s of measured elapsed
            # each; G3's post-charge overdraws the ~0.2s left and is refused
            # (a refused post-charge never reverts a finished query, so G3
            # stays "executed"); G4 is never even attempted because the loop
            # already stopped -- if cpu_seconds (5.0) had been charged
            # instead, G1 alone would have exhausted the budget.
            self.assertEqual({"G1", "G2", "G3"}, executed)
            self.assertEqual({"G4"}, deferred)
            self.assertEqual(3, backend.search_calls)
            self.assertEqual("wall_time", manifest["truncation_points"][0]["kind"])
            usage = account.freeze()["usage"]
            self.assertGreater(usage["wall_time_seconds"], 0.6)
            self.assertLess(usage["wall_time_seconds"], 2.0)

    def test_wall_time_remaining_never_negative_and_ledger_consistent(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            inventory, clusters, proteins = _simple_fixture(root, ["G1", "G2", "G3"])
            limits = {"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 50, "wall_time_limit_seconds": 1, "peak_memory_limit_mb": None}
            account = BudgetAccount(limits)
            # cpu is far larger than the real sleep and the wall limit; if it
            # were charged instead of the measured elapsed, the very first
            # query would refuse and usage would stay at 0.
            backend = FakeBackend(cpu=5.0, sleep_seconds=0.4)
            plan_sequence_search(inventory, clusters, account, root / "plan", CASE_ID, protein_db=root / "db.fa", search_backend=backend, proteins=proteins)
            frozen = account.freeze()
            remaining = frozen["remaining"]["wall_time_limit_seconds"]
            usage = frozen["usage"]
            self.assertGreaterEqual(remaining, 0)
            self.assertTrue(frozen["exhausted"]["wall_time_limit_seconds"])
            self.assertTrue(account.exhausted("wall_time_budget"))
            # two successful ~0.4s charges (~0.8s total), well under the 5.0s
            # cpu_seconds the backend claims -- proves elapsed, not
            # cpu_seconds, was charged.
            self.assertGreater(usage["wall_time_seconds"], 0.6)
            self.assertLess(usage["wall_time_seconds"], 2.0)
            self.assertEqual(frozen["limits"]["wall_time_limit_seconds"] - remaining, usage["wall_time_seconds"])

    def test_wall_time_budget_charges_measured_elapsed_even_when_backend_reports_zero_cpu(self) -> None:
        """Fix 3 acceptance test: a backend that reports cpu_seconds=0.0 but
        really sleeps must still burn the wall budget -- a custom/test
        backend cannot bypass the wall hard budget by misreporting
        cpu_seconds."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            gene_ids = [f"G{i}" for i in range(1, 11)]
            inventory, clusters, proteins = _simple_fixture(root, gene_ids)
            limits = {"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 50, "wall_time_limit_seconds": 1, "peak_memory_limit_mb": None}
            account = BudgetAccount(limits)
            backend = FakeBackend(cpu=0.0, sleep_seconds=0.3)
            plan_sequence_search(inventory, clusters, account, root / "plan", CASE_ID, protein_db=root / "db.fa", search_backend=backend, proteins=proteins)
            usage = account.freeze()["usage"]
            self.assertGreater(usage["wall_time_seconds"], 0)
            # not every planned query ran: the wall budget stopped the loop
            # even though the backend claimed zero cpu cost for every query.
            self.assertLess(backend.search_calls, len(gene_ids))
            rows = _plan_rows(root / "plan")
            self.assertGreater(sum(1 for row in rows.values() if row["query_status"] == "deferred_budget"), 0)

    def test_effective_query_timeout_capped_by_wall_budget(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            inventory, clusters, proteins = _simple_fixture(root, ["G1"])
            limits = {"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 50, "wall_time_limit_seconds": 3600, "peak_memory_limit_mb": None}
            backend = FakeBackend()
            # limits are integer-only; the sub-second timeout comes from the
            # external-tool timeout env override (capped by the wall budget)
            with mock.patch.dict(os.environ, {"CANDIDATE_QTL_TOOL_TIMEOUT": "0.5"}):
                plan_sequence_search(inventory, clusters, BudgetAccount(limits), root / "plan", CASE_ID, protein_db=root / "db.fa", search_backend=backend, proteins=proteins)
            self.assertEqual([0.5], backend.seen_timeouts)

    def test_hanging_provider_timeout_records_provider_errors_deterministically(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            hang = _write_script(root, "hang_blastp", 'if [ "$1" = "-version" ]; then echo "fake 1.0"; exit 0; fi\nsleep 60\n')
            inventory, clusters, proteins = _simple_fixture(root, ["G1"])
            with mock.patch.dict(os.environ, {"CANDIDATE_QTL_TOOL_TIMEOUT": "0.5"}):
                backend = BlastpBackend(cache=None, bin=str(hang))
                out = root / "plan"
                returned = plan_sequence_search(inventory, clusters, BudgetAccount(BUDGET_PROFILES["student_cpu"]), out, CASE_ID, protein_db=root / "db.fa", search_backend=backend, proteins=proteins)
            self.assertTrue(backend.available)
            rows = _plan_rows(out)
            self.assertEqual("provider_timeout", rows["G1"]["search_result_status"])
            self.assertEqual(1, len(returned["provider_errors"]))
            self.assertEqual("G1", returned["provider_errors"][0]["gene_id"])
            self.assertEqual("timeout", returned["provider_errors"][0]["kind"])
            written = json.loads((out / "tier2_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(1, written["provider_error_count"])
            self.assertNotIn("provider_errors", written)
            out2 = root / "plan2"
            with mock.patch.dict(os.environ, {"CANDIDATE_QTL_TOOL_TIMEOUT": "0.5"}):
                plan_sequence_search(inventory, clusters, BudgetAccount(BUDGET_PROFILES["student_cpu"]), out2, CASE_ID, protein_db=root / "db.fa", search_backend=BlastpBackend(cache=None, bin=str(hang)), proteins=proteins)
            self.assertEqual((out / "tier2_manifest.json").read_bytes(), (out2 / "tier2_manifest.json").read_bytes())

    def test_peak_memory_budget_blocks_first_planned_row(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            inventory, clusters, proteins = _simple_fixture(root, ["G1", "G2"])
            limits = {"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 50, "wall_time_limit_seconds": 3600, "peak_memory_limit_mb": 1024}
            backend = FakeBackend()
            account = BudgetAccount(limits)
            with mock.patch("candidate_qtl.tier2_planner._current_rss_mb", return_value=(99999.0, None)):
                manifest = plan_sequence_search(inventory, clusters, account, root / "plan", CASE_ID, protein_db=root / "db.fa", search_backend=backend, proteins=proteins)
            rows = _plan_rows(root / "plan")
            self.assertEqual(0, backend.search_calls)
            deferred = [row for row in rows.values() if row["query_status"] == "deferred_budget"]
            self.assertEqual(2, len(deferred))
            self.assertEqual("peak_memory", manifest["truncation_points"][0]["kind"])
            self.assertEqual("enforced", manifest["peak_memory_enforcement"])
            # Fix 1: the ledger now records the true observed peak and marks
            # the limit exhausted, instead of staying at 0/false while the
            # run was actually stopped for memory.
            frozen = account.freeze()
            self.assertEqual(99999.0, frozen["usage"]["peak_memory_mb"])
            self.assertTrue(frozen["exhausted"]["peak_memory_limit_mb"])
            self.assertEqual(0, frozen["remaining"]["peak_memory_limit_mb"])

    def test_peak_memory_not_obtainable_continues_execution(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            inventory, clusters, proteins = _simple_fixture(root, ["G1", "G2"])
            limits = {"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 50, "wall_time_limit_seconds": 3600, "peak_memory_limit_mb": 1024}
            backend = FakeBackend()
            with mock.patch("candidate_qtl.tier2_planner._current_rss_mb", return_value=(None, "not_obtainable")):
                manifest = plan_sequence_search(inventory, clusters, BudgetAccount(limits), root / "plan", CASE_ID, protein_db=root / "db.fa", search_backend=backend, proteins=proteins)
            rows = _plan_rows(root / "plan")
            self.assertEqual(2, backend.search_calls)
            executed = [row for row in rows.values() if row["query_status"] == "executed"]
            self.assertEqual(2, len(executed))
            self.assertEqual("not_obtainable", manifest["peak_memory_enforcement"])

    def test_peak_memory_post_check_stops_after_overrun(self) -> None:
        """Fix 4: peak memory is a hard limit. A run that starts under the
        declared limit but ends over it stops after that query: the post-check
        records the truncation point and the remaining rows are deferred."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            inventory, clusters, proteins = _simple_fixture(root, ["G1", "G2"])
            limits = {"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 50, "wall_time_limit_seconds": None, "peak_memory_limit_mb": 1024}
            backend = FakeBackend()
            account = BudgetAccount(limits)
            rss_values = iter([100.0, 3000.0])

            def fake_rss():
                return next(rss_values), None

            with mock.patch("candidate_qtl.tier2_planner._current_rss_mb", side_effect=fake_rss):
                manifest = plan_sequence_search(inventory, clusters, account, root / "plan", CASE_ID, protein_db=root / "db.fa", search_backend=backend, proteins=proteins)
            rows = _plan_rows(root / "plan")
            self.assertEqual(1, backend.search_calls)
            self.assertEqual("executed", rows["G1"]["query_status"])
            self.assertEqual("deferred_budget", rows["G2"]["query_status"])
            self.assertEqual("budget_exhausted", rows["G2"]["deferred_reason"])
            self.assertEqual(1, len(manifest["truncation_points"]))
            point = manifest["truncation_points"][0]
            self.assertEqual("peak_memory", point["kind"])
            self.assertEqual(1, point["plan_step"])
            self.assertEqual(3000.0, point["remaining_before"])
            self.assertEqual("enforced", manifest["peak_memory_enforcement"])
            self.assertEqual(1, manifest["queries_executed"])
            self.assertEqual(1, manifest["queries_deferred_budget"])
            frozen = account.freeze()
            # Fix 1: the ledger records the observed maximum (3000.0, the
            # post-run overrun measurement), not 0, and marks the limit
            # exhausted; remaining is clamped at 0, never negative.
            self.assertEqual(3000.0, frozen["usage"]["peak_memory_mb"])
            self.assertTrue(frozen["exhausted"]["peak_memory_limit_mb"])
            self.assertEqual(0, frozen["remaining"]["peak_memory_limit_mb"])

    def test_provider_timeout_charges_wall_budget_when_declared(self) -> None:
        """Fix 1: a provider timeout charges the actual elapsed wall time, so a
        run with a declared wall limit cannot burn unbounded wall seconds."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            inventory, clusters, proteins = _simple_fixture(root, ["G1", "G2", "G3"])
            limits = {"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 50, "wall_time_limit_seconds": 1, "peak_memory_limit_mb": None}
            backend = TimeoutBackend(sleep_seconds=0.15)
            account = BudgetAccount(limits)
            manifest = plan_sequence_search(inventory, clusters, account, root / "plan", CASE_ID, protein_db=root / "db.fa", search_backend=backend, proteins=proteins)
            rows = _plan_rows(root / "plan")
            self.assertEqual(3, backend.search_calls)
            for row in rows.values():
                self.assertEqual("provider_timeout", row["search_result_status"])
                # Fix 4: provider_timeout rows are now the terminal
                # "provider_failed" query_status, so they land in exactly one
                # accounted-for bucket instead of dangling as "planned".
                self.assertEqual("provider_failed", row["query_status"])
                self.assertEqual("provider_timeout", row["deferred_reason"])
            self.assertEqual(3, len(manifest["provider_errors"]))
            self.assertTrue(all(error["kind"] == "timeout" for error in manifest["provider_errors"]))
            self.assertEqual([], manifest["truncation_points"])
            self.assertEqual(3, manifest["queries_provider_failed"])
            frozen = account.freeze()
            self.assertGreater(frozen["usage"]["wall_time_seconds"], 0.3)
            self.assertLess(frozen["remaining"]["wall_time_limit_seconds"], 1.0)
            self.assertFalse(frozen["exhausted"]["wall_time_limit_seconds"])
            written = json.loads((root / "plan" / "tier2_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(3, written["provider_error_count"])
            self.assertEqual(3, written["queries_provider_failed"])

    def test_provider_timeout_charge_refusal_stops_loop(self) -> None:
        """Fix 1 (b): when a provider hang's elapsed time overdraws the wall
        budget, the charge is refused, the loop stops, and the remaining rows
        become deferred_budget with a wall_time truncation point."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            inventory, clusters, proteins = _simple_fixture(root, ["G1", "G2", "G3"])
            limits = {"llm_call_limit": 0, "token_limit": 0, "sequence_query_limit": 50, "wall_time_limit_seconds": 1, "peak_memory_limit_mb": None}
            backend = TimeoutBackend(sleep_seconds=1.2)
            account = BudgetAccount(limits)
            manifest = plan_sequence_search(inventory, clusters, account, root / "plan", CASE_ID, protein_db=root / "db.fa", search_backend=backend, proteins=proteins)
            rows = _plan_rows(root / "plan")
            self.assertEqual(1, backend.search_calls)
            self.assertEqual("provider_failed", rows["G1"]["query_status"])
            self.assertEqual("provider_timeout", rows["G1"]["search_result_status"])
            self.assertEqual("provider_timeout", rows["G1"]["deferred_reason"])
            for gene in ("G2", "G3"):
                self.assertEqual("deferred_budget", rows[gene]["query_status"])
                self.assertEqual("budget_exhausted", rows[gene]["deferred_reason"])
                self.assertEqual("not_run", rows[gene]["search_result_status"])
            self.assertEqual(1, len(manifest["truncation_points"]))
            point = manifest["truncation_points"][0]
            self.assertEqual("wall_time", point["kind"])
            self.assertEqual(1, point["plan_step"])
            self.assertGreater(point["remaining_before"], 0.9)
            self.assertEqual(0, manifest["queries_executed"])
            self.assertEqual(2, manifest["queries_deferred_budget"])
            self.assertEqual(1, manifest["queries_provider_failed"])
            frozen = account.freeze()
            # the refused charge charged nothing; the budget is exhausted
            self.assertEqual(0, frozen["usage"]["wall_time_seconds"])
            self.assertTrue(frozen["exhausted"]["wall_time_limit_seconds"])

    def test_win32_reserved_device_names_fallback_to_gene_stem(self) -> None:
        """Fix 3: query stems that are win32 reserved device names fall back to
        "gene" (plus the collision suffix) so no OSError can occur on Windows."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            gene_ids = ["con", "CON", "nul", "aux", "com1", "lpt9", "G1"]
            inventory, clusters, proteins = _simple_fixture(root, gene_ids)
            out = root / "plan"
            manifest = plan_sequence_search(inventory, clusters, BudgetAccount(BUDGET_PROFILES["student_cpu"]), out, CASE_ID, protein_db=root / "db.fa", search_backend=FakeBackend(), proteins=proteins)
            queries_dir = (out / "queries").resolve()
            frame = pd.read_csv(out / "queries" / "query_file_map.tsv", sep="\t", dtype=str)
            self.assertEqual(len(gene_ids), len(frame))
            self.assertEqual(set(gene_ids), set(frame["gene_id"]))
            reserved = {"con", "nul", "aux", "prn"} | {f"com{i}" for i in range(1, 10)} | {f"lpt{i}" for i in range(1, 10)}
            for _, row in frame.iterrows():
                stem = row["stem"]
                self.assertNotIn(stem.rstrip(" .").lower(), reserved)
                query_file = (queries_dir / f"{stem}.fa").resolve()
                self.assertTrue(query_file.is_file())
                self.assertTrue(query_file.is_relative_to(queries_dir))
            self.assertEqual(7, manifest["queries_executed"])
            self.assertEqual(7, len(sorted((out / "queries").glob("*.fa"))))
            headers = {path.read_text(encoding="utf-8").splitlines()[0][1:] for path in (out / "queries").glob("*.fa")}
            self.assertEqual(set(gene_ids), headers)

    def test_safe_query_paths_keep_original_identity(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            gene_ids = ["../evil", "C:/abs", "a/b", "..", ".", "G1"]
            inventory, clusters, proteins = _simple_fixture(root, gene_ids)
            out = root / "plan"
            plan_sequence_search(inventory, clusters, BudgetAccount(BUDGET_PROFILES["student_cpu"]), out, CASE_ID, protein_db=root / "db.fa", search_backend=FakeBackend(), proteins=proteins)
            queries_dir = (out / "queries").resolve()
            fastas = sorted((out / "queries").glob("*.fa"))
            self.assertEqual(len(gene_ids), len(fastas))
            for path in fastas:
                self.assertTrue(path.resolve().is_relative_to(queries_dir))
            headers = {path.read_text(encoding="utf-8").splitlines()[0][1:] for path in fastas}
            self.assertEqual(set(gene_ids), headers)
            rows = _plan_rows(out)
            for gene in gene_ids:
                self.assertEqual(gene, rows[gene]["gene_id"])
                self.assertEqual("hit", rows[gene]["search_result_status"])

    def test_colliding_stems_get_deterministic_suffixes(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            inventory, clusters, proteins = _simple_fixture(root, ["a/b", "a_b"])
            for label in ("plan", "plan2"):
                plan_sequence_search(inventory, clusters, BudgetAccount(BUDGET_PROFILES["student_cpu"]), root / label, CASE_ID, protein_db=root / "db.fa", search_backend=FakeBackend(), proteins=proteins)
            names = sorted(path.name for path in (root / "plan" / "queries").glob("*.fa"))
            names2 = sorted(path.name for path in (root / "plan2" / "queries").glob("*.fa"))
            self.assertEqual(names, names2)
            self.assertEqual(2, len(names))
            self.assertIn("a_b.fa", names)
            self.assertTrue(any(name.startswith("a_b-") and len(name) == len("a_b-") + 8 + len(".fa") for name in names))

    def test_query_file_map_tsv_columns_and_coverage(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inventory_rows, cluster_rows, _extract = _tier0_tier1(root, self)
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            plan_sequence_search(inventory_rows, cluster_rows, BudgetAccount(BUDGET_PROFILES["student_cpu"]), root / "plan", CASE_ID, protein_db=root / "db.fa", search_backend=FakeBackend(), proteins=_all_proteins())
            frame = pd.read_csv(root / "plan" / "queries" / "query_file_map.tsv", sep="\t", dtype=str)
            self.assertEqual(["gene_id", "stem", "query_file", "payload_file"], list(frame.columns))
            plan = _plan_rows(root / "plan")
            planned = sorted(gene for gene, row in plan.items() if row["queries_planned"] == "1")
            self.assertEqual(planned, sorted(frame["gene_id"]))
            for _, row in frame.iterrows():
                self.assertTrue((root / "plan" / row["query_file"]).is_file())
                self.assertTrue(row["payload_file"].startswith("queries/"))
                self.assertTrue(row["payload_file"].endswith(".search_payload.json"))

    def test_plan_sequence_search_lazy_backend_uses_custom_blastp_bin(self) -> None:
        """Fix 5: when plan_sequence_search constructs its own BlastpBackend
        (search_backend=None), a custom blastp_bin must actually be used --
        not silently dropped in favor of the default "blastp" on PATH."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "db.fa").write_text(">db\nMK\n", encoding="utf-8")
            script = _write_script(
                root, "fake_blastp",
                'if [ "$1" = "-version" ]; then echo "fake-blastp 9.9"; exit 0; fi\n'
                'printf "G1\\thit1\\t95\\t1e-50\\t200\\n"\n',
            )
            inventory, clusters, proteins = _simple_fixture(root, ["G1"])
            out = root / "plan"
            manifest = plan_sequence_search(
                inventory, clusters, BudgetAccount(BUDGET_PROFILES["student_cpu"]), out, CASE_ID,
                protein_db=root / "db.fa", blastp_bin=str(script), search_backend=None, proteins=proteins,
            )
            rows = _plan_rows(out)
            self.assertEqual("executed", rows["G1"]["query_status"])
            self.assertEqual("hit", rows["G1"]["search_result_status"])
            self.assertEqual(1, manifest["queries_executed"])
            self.assertEqual("blastp", manifest["search_provider"])


if __name__ == "__main__":
    unittest.main()
