from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.tier0_inventory import TIER0_COLUMNS, build_tier0_inventory, load_provider_table
from candidate_qtl.tier1_cluster import (
    ALLOCATIONS,
    RESCUE_BUCKETS,
    RULE_KEYS,
    TIER1_COLUMNS,
    build_tier1_clusters,
    dep_id,
)

CASE_ID = "tier1-case"


def _tier1_genes():
    return [
        {"id": "G1", "start": 11, "end": 22, "cds_start": 11, "cds_end": 22, "cds": "ATGGAAAACTAA", "tx": "Tx1", "attrs": "ID=G1;Name=GeneA;Note=desc1"},
        {"id": "G2", "start": 31, "end": 45, "cds_start": 31, "cds_end": 45, "cds": "ATGGAAAACGCCTAA", "tx": "Tx2", "attrs": "ID=G2;Note=desc2"},
        {"id": "G3", "start": 51, "end": 60, "cds": None, "tx": None, "attrs": "ID=G3"},
        {"id": "G4", "start": 71, "end": 80, "cds_start": 71, "cds_end": 80, "cds": "ATGGAAAACA", "tx": "Tx4", "attrs": "ID=G4;Note=desc4"},
        {"id": "G5", "start": 91, "end": 105, "cds_start": 91, "cds_end": 105, "cds": "ATGGAAAACGATTAA", "tx": "Tx5", "attrs": "ID=G5;Name=GeneC;Note=desc5"},
        {"id": "G6", "start": 111, "end": 125, "cds_start": 111, "cds_end": 125, "cds": "ATGGAAAACGATTAA", "tx": "Tx6", "attrs": "ID=G6;Name=GeneA;Note=desc6"},
        {"id": "G7", "start": 131, "end": 145, "cds_start": 131, "cds_end": 145, "cds": "ATGGAAAACGATTAA", "tx": "Tx7", "attrs": "ID=G7;Name=GeneA;Note=desc7"},
        {"id": "G8", "start": 151, "end": 165, "cds_start": 151, "cds_end": 165, "cds": "ATGGAAAACGATTAA", "tx": "Tx8", "attrs": "ID=G8;Name=GeneA;Note=desc8"},
        {"id": "G9", "start": 171, "end": 185, "cds_start": 171, "cds_end": 185, "cds": "ATGGAAAACGATTAA", "tx": "Tx9", "attrs": "ID=G9;Name=GeneB;Note=desc9"},
        {"id": "G10", "start": 191, "end": 205, "cds_start": 191, "cds_end": 205, "cds": "ATGGAAAACGATTAA", "tx": "Tx10", "attrs": "ID=G10;Note=desc10"},
        {"id": "G11", "start": 211, "end": 225, "cds_start": 211, "cds_end": 225, "cds": "ATGGAAAACGATTAA", "tx": "Tx11", "attrs": "ID=G11;Note=desc11"},
        {"id": "G12", "start": 231, "end": 245, "cds_start": 231, "cds_end": 245, "cds": "ATGGAAAACGATTAA", "tx": "Tx12", "attrs": "ID=G12;Note=desc12"},
        {"id": "G13", "start": 251, "end": 265, "cds_start": 251, "cds_end": 265, "cds": "ATGGAAAACGATTAA", "tx": "Tx13", "attrs": "ID=G13;Note=desc13"},
    ]


def _write_extract_fixture(root: Path, genes, length: int = 300) -> tuple[Path, Path]:
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


def _write_provider_tables(root: Path) -> dict[str, dict[str, object]]:
    (root / "domains.tsv").write_text("gene_id\tdomain_id\nG2\tPF001\nG4\tPF001\n", encoding="utf-8")
    (root / "orthogroup.tsv").write_text("gene_id\torthogroup\nG10\tOG001\nG11\tOG001\n", encoding="utf-8")
    (root / "variant.tsv").write_text("gene_id\teffect\tsource_id\nG2\tmissense\tX1\nG4\tmissense\tX1\n", encoding="utf-8")
    (root / "sv.tsv").write_text("gene_id\tsv_type\tsource_id\nG12\tDEL\tSVDB1\nG12\tINS\tSVDB2\nG12\tDUP\tSVDB1\n", encoding="utf-8")
    (root / "family.tsv").write_text("gene_id\tfamily_id\nG13\tFAM0001\nG5\tFAM0002\n", encoding="utf-8")
    return {
        "domains": load_provider_table("domains", str(root / "domains.tsv"), "pfam", "v1"),
        "orthogroup": load_provider_table("orthogroup", str(root / "orthogroup.tsv"), "eggnog", "v1"),
        "variant": load_provider_table("variant", str(root / "variant.tsv"), "var-db", "v1"),
        "sv": load_provider_table("sv", str(root / "sv.tsv"), "sv-db", "v1"),
        "family": load_provider_table("family", str(root / "family.tsv"), "fam-db", "v1"),
    }


def _read_clusters(path: Path) -> dict[frozenset[str], dict[str, str]]:
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return {frozenset(row["member_gene_ids"].split(";")): dict(row) for row in frame.to_dict("records")}


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


class Tier1ClusterTests(unittest.TestCase):
    def _build(self, root: Path, out_name: str) -> dict[str, object]:
        fasta, gff = _write_extract_fixture(root, _tier1_genes())
        extract_out = root / "extract"
        self.assertEqual(
            0,
            main(["extract", "--fasta", str(fasta), "--gff3", str(gff), "--interval", "chr1:1-1000", "--out-dir", str(extract_out)]),
        )
        providers = _write_provider_tables(root)
        tier0_out = root / "tier0"
        build_tier0_inventory(extract_out, tier0_out, providers, CASE_ID)
        frame = pd.read_csv(tier0_out / "tier0_inventory.tsv", sep="\t", dtype=str).fillna("")
        out = root / out_name
        return build_tier1_clusters(frame.to_dict("records"), out, CASE_ID, providers=providers), out

    def test_rule_order_exhaustive_and_exclusive(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest, out = self._build(root, "tier1")
            clusters = _read_clusters(out / "tier1_clusters.tsv")
            self.assertEqual(8, len(clusters))

            annotation_a = clusters[frozenset({"G1", "G6", "G7", "G8"})]
            self.assertEqual("annotation", annotation_a["rule_key"])
            self.assertEqual("GeneA", annotation_a["cluster_key"])
            self.assertEqual("priority", annotation_a["allocation"])
            self.assertEqual("shared_annotation", annotation_a["allocation_reason"])
            self.assertEqual("4", annotation_a["member_count"])
            self.assertEqual("annotation:GeneA", annotation_a["seed_reason"])
            self.assertEqual("G6", annotation_a["representative_gene_id"])
            self.assertEqual("model_status=complete_orf;cds_length=15;protein_length=4", annotation_a["representative_reason"])
            self.assertRegex(annotation_a["cluster_id"], r"^cls_[0-9a-f]{20}$")

            annotation_c = clusters[frozenset({"G5"})]
            self.assertEqual("family", annotation_c["rule_key"])
            self.assertEqual("FAM0002", annotation_c["cluster_key"])
            self.assertEqual("family:FAM0002", annotation_c["seed_reason"])
            self.assertEqual("priority", annotation_c["allocation"])
            self.assertEqual("shared_family", annotation_c["allocation_reason"])
            self.assertEqual("FAM0002", annotation_c["family_ids"])

            provider_gap = clusters[frozenset({"G9"})]
            self.assertEqual("model_state_rescue", provider_gap["rule_key"])
            self.assertEqual("rescue", provider_gap["allocation"])
            self.assertEqual("provider_gap", provider_gap["rescue_bucket"])
            self.assertEqual("model_state_rescue:provider_gap", provider_gap["allocation_reason"])
            self.assertEqual("model_state_rescue:G9", provider_gap["seed_reason"])

            domain = clusters[frozenset({"G2", "G4"})]
            self.assertEqual("domain", domain["rule_key"])
            self.assertEqual("PF001", domain["cluster_key"])
            self.assertEqual("priority", domain["allocation"])
            self.assertEqual("shared_domain", domain["allocation_reason"])
            self.assertEqual("PF001", domain["domain_ids"])
            self.assertEqual(dep_id("X1"), domain["dependency_cluster_ids"])
            self.assertEqual("G2", domain["representative_gene_id"])
            self.assertEqual("model_status=complete_orf;cds_length=15;protein_length=4", domain["representative_reason"])

            orthogroup = clusters[frozenset({"G10", "G11"})]
            self.assertEqual("orthogroup", orthogroup["rule_key"])
            self.assertEqual("OG001", orthogroup["cluster_key"])
            self.assertEqual("defer", orthogroup["allocation"])
            self.assertEqual("shared_orthogroup", orthogroup["allocation_reason"])
            self.assertEqual("OG001", orthogroup["orthogroup_ids"])
            self.assertEqual("G10", orthogroup["representative_gene_id"])

            rescue_missing_model = clusters[frozenset({"G3"})]
            self.assertEqual("model_state_rescue", rescue_missing_model["rule_key"])
            self.assertEqual("rescue", rescue_missing_model["allocation"])
            self.assertEqual("missing_model", rescue_missing_model["rescue_bucket"])
            self.assertEqual("model_state_rescue:missing_model", rescue_missing_model["allocation_reason"])
            self.assertEqual("model_state_rescue:G3", rescue_missing_model["seed_reason"])

            rescue_sv = clusters[frozenset({"G12"})]
            self.assertEqual("model_state_rescue", rescue_sv["rule_key"])
            self.assertEqual("sv_pav_regulatory", rescue_sv["rescue_bucket"])
            self.assertEqual("model_state_rescue:sv_pav_regulatory", rescue_sv["allocation_reason"])
            self.assertEqual(";".join(sorted([dep_id("SVDB1"), dep_id("SVDB2")])), rescue_sv["dependency_cluster_ids"])

            singleton = clusters[frozenset({"G13"})]
            self.assertEqual("family", singleton["rule_key"])
            self.assertEqual("FAM0001", singleton["cluster_key"])
            self.assertEqual("family:FAM0001", singleton["seed_reason"])
            self.assertEqual("priority", singleton["allocation"])
            self.assertEqual("shared_family", singleton["allocation_reason"])
            self.assertEqual("FAM0001", singleton["family_ids"])

            members = {gene for cluster in clusters.values() for gene in cluster["member_gene_ids"].split(";")}
            self.assertEqual(13, len(members))
            self.assertEqual(13, manifest["genes_in_cluster_count"])
            self.assertEqual(13, manifest["gene_count"])

    def test_manifest_counts_and_deterministic_ids(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest, out = self._build(root, "tier1")
            self.assertEqual(8, manifest["cluster_count"])
            # G9 mixes "missing" (domains/orthogroup/variant/sv provided) and
            # "not_evaluated" (expression/pav/promoter not provided) cells ->
            # provider_gap rescue under the new semantics; nothing is left for
            # the singleton rule.
            self.assertEqual({"defer": 1, "priority": 4, "rescue": 3}, manifest["allocation_counts"])
            self.assertEqual(
                {"missing_model": 1, "not_applicable": 5, "provider_gap": 1, "sv_pav_regulatory": 1},
                manifest["rescue_bucket_counts"],
            )
            self.assertEqual(
                {"annotation": 1, "domain": 1, "family": 2, "model_state_rescue": 3, "orthogroup": 1},
                manifest["rule_key_counts"],
            )
            self.assertEqual(3, manifest["dependency_cluster_count"])
            self.assertEqual(6, manifest["empty_dependency_source_count"])
            self.assertEqual("1.0", manifest["schema_version"])
            self.assertEqual(CASE_ID, manifest["case_id"])
            clusters = _read_clusters(out / "tier1_clusters.tsv")
            for cluster in clusters.values():
                self.assertRegex(cluster["cluster_id"], r"^cls_[0-9a-f]{20}$")
                self.assertIn(cluster["allocation"], ALLOCATIONS)
                self.assertIn(cluster["rescue_bucket"], RESCUE_BUCKETS | {"not_applicable"})
                self.assertIn(cluster["rule_key"], RULE_KEYS)
            self.assertEqual(TIER1_COLUMNS, tuple(pd.read_csv(out / "tier1_clusters.tsv", sep="\t", dtype=str).columns))

    def test_byte_identical_double_run(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            first, out_first = self._build(root, "first")
            second, out_second = self._build(root, "second")
            self.assertEqual(first, second)
            self.assertEqual((out_first / "tier1_clusters.tsv").read_bytes(), (out_second / "tier1_clusters.tsv").read_bytes())
            self.assertEqual((out_first / "tier1_manifest.json").read_bytes(), (out_second / "tier1_manifest.json").read_bytes())

    def test_no_file_writes_outside_out_dir(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest, out = self._build(root, "tier1")
            self.assertEqual({"tier1_clusters.tsv", "tier1_manifest.json"}, {path.name for path in out.iterdir()})
            self.assertTrue(manifest)

    def test_hand_built_rows_rescue_buckets_and_singletons(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            rows = [
                _hand_row("H1", annotation_status="missing"),
                _hand_row("H2"),
                _hand_row("H3", sv_status="present", sv_evidence="DEL"),
            ]
            out = root / "tier1"
            manifest = build_tier1_clusters(rows, out, CASE_ID)
            clusters = _read_clusters(out / "tier1_clusters.tsv")
            weak = clusters[frozenset({"H1"})]
            self.assertEqual("model_state_rescue", weak["rule_key"])
            self.assertEqual("weak_annotation", weak["rescue_bucket"])
            self.assertEqual("model_state_rescue:weak_annotation", weak["allocation_reason"])
            rescue_sv = clusters[frozenset({"H3"})]
            self.assertEqual("sv_pav_regulatory", rescue_sv["rescue_bucket"])
            singleton = clusters[frozenset({"H2"})]
            self.assertEqual("singleton", singleton["rule_key"])
            self.assertEqual("not_applicable", singleton["rescue_bucket"])
            self.assertEqual(3, manifest["genes_in_cluster_count"])
            self.assertEqual(3, manifest["gene_count"])
            self.assertEqual({"rescue": 2, "defer": 1}, manifest["allocation_counts"])

    def test_hand_built_representative_defaults_without_cds_length(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            rows = [
                _hand_row("X1", gene_name="GeneZ"),
                _hand_row("X2", gene_name="GeneZ"),
            ]
            out = root / "tier1"
            manifest = build_tier1_clusters(rows, out, CASE_ID)
            clusters = _read_clusters(out / "tier1_clusters.tsv")
            cluster = clusters[frozenset({"X1", "X2"})]
            self.assertEqual("annotation", cluster["rule_key"])
            self.assertEqual("GeneZ", cluster["cluster_key"])
            self.assertEqual("2", cluster["member_count"])
            self.assertEqual("X1", cluster["representative_gene_id"])
            self.assertEqual("model_status=complete_orf;cds_length=0;protein_length=0", cluster["representative_reason"])
            self.assertEqual(1, manifest["cluster_count"])
            self.assertEqual({"priority": 1}, manifest["allocation_counts"])

    def test_builder_provider_gap_bucket_reachable(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            genes = [
                {"id": "P1", "start": 11, "end": 22, "cds_start": 11, "cds_end": 22, "cds": "ATGGAAAACTAA", "tx": "Tx1", "attrs": "ID=P1;Note=desc1"},
                {"id": "P2", "start": 31, "end": 45, "cds_start": 31, "cds_end": 45, "cds": "ATGGAAAACGCCTAA", "tx": "Tx2", "attrs": "ID=P2;Note=desc2"},
                {"id": "P3", "start": 51, "end": 65, "cds_start": 51, "cds_end": 65, "cds": "ATGGAAAACGATTAA", "tx": "Tx3", "attrs": "ID=P3;Note=desc3"},
            ]
            fasta, gff = _write_extract_fixture(root, genes)
            extract_out = root / "extract"
            self.assertEqual(
                0,
                main(["extract", "--fasta", str(fasta), "--gff3", str(gff), "--interval", "chr1:1-1000", "--out-dir", str(extract_out)]),
            )
            # every provider kind is provided, but rows exist only for P1:
            # P2/P3 have provider_applicability "applicable" with every
            # resource status cell "missing" -> provider_gap bucket.
            for kind, value_column in (
                ("domains", "domain_id"),
                ("orthogroup", "orthogroup"),
                ("expression", "value"),
                ("variant", "effect"),
                ("sv", "sv_type"),
                ("pav", "pav_type"),
                ("promoter", "motif"),
            ):
                (root / f"{kind}.tsv").write_text(f"gene_id\t{value_column}\nP1\thit1\n", encoding="utf-8")
            providers = {
                kind: load_provider_table(kind, str(root / f"{kind}.tsv"), "prov-db", "v1")
                for kind in ("domains", "orthogroup", "expression", "variant", "sv", "pav", "promoter")
            }
            tier0_out = root / "tier0"
            build_tier0_inventory(extract_out, tier0_out, providers, CASE_ID)
            frame = pd.read_csv(tier0_out / "tier0_inventory.tsv", sep="\t", dtype=str).fillna("")
            rows = {row["gene_id"]: row for row in frame.to_dict("records")}
            p3 = rows["P3"]
            self.assertEqual("complete_orf", p3["model_status"])
            self.assertEqual("present", p3["annotation_status"])
            self.assertEqual("applicable", p3["provider_applicability"])
            for cell in ("domain_status", "orthogroup_status", "expression_status", "variant_status", "sv_status", "pav_status", "promoter_status"):
                self.assertEqual("missing", p3[cell])
            out = root / "tier1"
            manifest = build_tier1_clusters(frame.to_dict("records"), out, CASE_ID, providers=providers)
            clusters = _read_clusters(out / "tier1_clusters.tsv")
            provider_gap = clusters[frozenset({"P3"})]
            self.assertEqual("model_state_rescue", provider_gap["rule_key"])
            self.assertEqual("rescue", provider_gap["allocation"])
            self.assertEqual("provider_gap", provider_gap["rescue_bucket"])
            self.assertEqual("model_state_rescue:provider_gap", provider_gap["allocation_reason"])
            self.assertEqual({"not_applicable": 1, "provider_gap": 2}, manifest["rescue_bucket_counts"])
            self.assertEqual({"priority": 1, "rescue": 2}, manifest["allocation_counts"])

    def test_builder_unclassified_bucket_reachable_for_unevaluated_model(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out = root / "extract"
            extract_out.mkdir()
            # hand-written candidates.tsv: U2 has an empty/absent
            # annotation_quality -> model_status "not_evaluated"
            (extract_out / "candidates.tsv").write_text(
                "gene_id\tchrom\tstart\tend\tstrand\tgene_name\tdescription\tcds_length\tannotation_quality\n"
                "U1\tchr1\t1\t10\t+\t\tName1\t12\tcomplete_orf\n"
                "U2\tchr1\t21\t30\t+\t\tName2\t12\n",
                encoding="utf-8",
            )
            tier0_out = root / "tier0"
            build_tier0_inventory(extract_out, tier0_out, {}, CASE_ID)
            frame = pd.read_csv(tier0_out / "tier0_inventory.tsv", sep="\t", dtype=str).fillna("")
            rows = {row["gene_id"]: row for row in frame.to_dict("records")}
            self.assertEqual("not_evaluated", rows["U2"]["model_status"])
            self.assertEqual("present", rows["U2"]["annotation_status"])
            self.assertEqual("not_evaluated", rows["U2"]["provider_applicability"])
            out = root / "tier1"
            manifest = build_tier1_clusters(frame.to_dict("records"), out, CASE_ID)
            clusters = _read_clusters(out / "tier1_clusters.tsv")
            unclassified = clusters[frozenset({"U2"})]
            self.assertEqual("model_state_rescue", unclassified["rule_key"])
            self.assertEqual("rescue", unclassified["allocation"])
            self.assertEqual("unclassified", unclassified["rescue_bucket"])
            self.assertEqual("model_state_rescue:unclassified", unclassified["allocation_reason"])
            self.assertEqual("singleton", clusters[frozenset({"U1"})]["rule_key"])
            self.assertEqual({"not_applicable": 1, "unclassified": 1}, manifest["rescue_bucket_counts"])
            self.assertEqual({"rescue": 1, "defer": 1}, manifest["allocation_counts"])

    def test_dep_id_formula_and_dedup(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.assertEqual("dep_" + hashlib.sha256("X1".encode("utf-8")).hexdigest()[:20], dep_id("X1"))
            self.assertEqual(dep_id("X1"), dep_id("X1"))
            self.assertNotEqual(dep_id("X1"), dep_id("X2"))
            rows = [
                _hand_row("D1", gene_name="GeneD", variant_source_id="B1;A1"),
                _hand_row("D2", gene_name="GeneD", variant_source_id="A1"),
            ]
            out = root / "tier1"
            manifest = build_tier1_clusters(rows, out, CASE_ID)
            clusters = _read_clusters(out / "tier1_clusters.tsv")
            cluster = clusters[frozenset({"D1", "D2"})]
            self.assertEqual(";".join(sorted([dep_id("A1"), dep_id("B1")])), cluster["dependency_cluster_ids"])
            self.assertEqual(2, manifest["dependency_cluster_count"])
            self.assertEqual(0, manifest["empty_dependency_source_count"])

    def test_dependency_manifest_json(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest, out = self._build(root, "tier1")
            written = json.loads((out / "tier1_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest, written)
            self.assertNotIn("timestamp", str(written))


class Tier1RemediationTests(unittest.TestCase):
    """Audit findings G1/G2/G4/G5: unique gene names never auto-advantage,
    family is a real cluster key, mixed missing/not_evaluated cells still
    rescue as provider_gap, and representative reasons preserve conflicts."""

    def _family_table(self, root: Path, content: str, name: str = "family.tsv") -> dict[str, object]:
        path = root / name
        path.write_text(content, encoding="utf-8")
        return load_provider_table("family", str(path), "fam-db", "v1")

    def test_unique_name_singleton_not_priority(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            rows = [_hand_row("S1", gene_name="UniqueName")]
            out = root / "tier1"
            manifest = build_tier1_clusters(rows, out, CASE_ID)
            clusters = _read_clusters(out / "tier1_clusters.tsv")
            cluster = clusters[frozenset({"S1"})]
            self.assertEqual("singleton", cluster["rule_key"])
            self.assertEqual("defer", cluster["allocation"])
            self.assertEqual("singleton:S1", cluster["seed_reason"])
            self.assertEqual(1, manifest["cluster_count"])
            self.assertEqual({"defer": 1}, manifest["allocation_counts"])

    def test_annotated_singleton_no_priority_over_unannotated(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            rows = [_hand_row("A1", gene_name="OnlyName"), _hand_row("A2")]
            out = root / "tier1"
            manifest = build_tier1_clusters(rows, out, CASE_ID)
            clusters = _read_clusters(out / "tier1_clusters.tsv")
            named = clusters[frozenset({"A1"})]
            unnamed = clusters[frozenset({"A2"})]
            for cluster in (named, unnamed):
                self.assertEqual("singleton", cluster["rule_key"])
                self.assertEqual("defer", cluster["allocation"])
                self.assertEqual("not_applicable", cluster["rescue_bucket"])
            self.assertEqual({"defer": 2}, manifest["allocation_counts"])

    def test_family_cluster_actual_key_two_members_one_cluster(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            family = self._family_table(root, "gene_id\tfamily_id\nG1\tFAM001\nG2\tFAM001\n")
            rows = [_hand_row("G1"), _hand_row("G2")]
            out = root / "tier1"
            manifest = build_tier1_clusters(rows, out, CASE_ID, providers={"family": family})
            clusters = _read_clusters(out / "tier1_clusters.tsv")
            self.assertEqual(1, len(clusters))
            cluster = clusters[frozenset({"G1", "G2"})]
            self.assertEqual("family", cluster["rule_key"])
            self.assertEqual("FAM001", cluster["cluster_key"])
            self.assertEqual("family:FAM001", cluster["seed_reason"])
            self.assertEqual("priority", cluster["allocation"])
            self.assertEqual("shared_family", cluster["allocation_reason"])
            self.assertEqual("2", cluster["member_count"])
            self.assertEqual("FAM001", cluster["family_ids"])
            self.assertIn(cluster["representative_gene_id"], {"G1", "G2"})
            self.assertEqual(1, manifest["cluster_count"])

    def test_family_precedence_over_annotation_and_domain(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            family = self._family_table(root, "gene_id\tfamily_id\nG1\tFAM100\nG2\tFAM100\n")
            rows = [
                _hand_row("G1", gene_name="GeneX", domain_hits="PF001"),
                _hand_row("G2", gene_name="GeneX", domain_hits="PF002"),
            ]
            out = root / "tier1"
            manifest = build_tier1_clusters(rows, out, CASE_ID, providers={"family": family})
            clusters = _read_clusters(out / "tier1_clusters.tsv")
            self.assertEqual(1, len(clusters))
            cluster = clusters[frozenset({"G1", "G2"})]
            self.assertEqual("family", cluster["rule_key"])
            self.assertEqual("FAM100", cluster["cluster_key"])
            self.assertEqual("priority", cluster["allocation"])
            self.assertEqual("shared_family", cluster["allocation_reason"])
            self.assertEqual("2", cluster["member_count"])
            self.assertEqual({"family": 1}, manifest["rule_key_counts"])

    def test_rule_order_family_annotation_domain_orthogroup_rescue_singleton(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            family = self._family_table(root, "gene_id\tfamily_id\nF1\tFAMF\nF2\tFAMF\n")
            rows = [
                _hand_row("F1"),
                _hand_row("F2"),
                _hand_row("A1", gene_name="GeneA"),
                _hand_row("A2", gene_name="GeneA"),
                _hand_row("D1", domain_hits="PF001"),
                _hand_row("D2", domain_hits="PF001"),
                _hand_row("O1", orthogroup="OG001"),
                _hand_row("O2", orthogroup="OG001"),
                _hand_row("R1", model_status="no_cds"),
                _hand_row("S1"),
            ]
            out = root / "tier1"
            manifest = build_tier1_clusters(rows, out, CASE_ID, providers={"family": family})
            clusters = _read_clusters(out / "tier1_clusters.tsv")
            self.assertEqual("family", clusters[frozenset({"F1", "F2"})]["rule_key"])
            self.assertEqual("annotation", clusters[frozenset({"A1", "A2"})]["rule_key"])
            self.assertEqual("domain", clusters[frozenset({"D1", "D2"})]["rule_key"])
            self.assertEqual("orthogroup", clusters[frozenset({"O1", "O2"})]["rule_key"])
            rescue = clusters[frozenset({"R1"})]
            self.assertEqual("model_state_rescue", rescue["rule_key"])
            self.assertEqual("missing_model", rescue["rescue_bucket"])
            singleton = clusters[frozenset({"S1"})]
            self.assertEqual("singleton", singleton["rule_key"])
            self.assertEqual("defer", singleton["allocation"])
            self.assertEqual(
                {"annotation": 1, "domain": 1, "family": 1, "model_state_rescue": 1, "orthogroup": 1, "singleton": 1},
                manifest["rule_key_counts"],
            )
            self.assertEqual(10, manifest["genes_in_cluster_count"])
            self.assertEqual(10, manifest["gene_count"])

    def test_every_gene_exactly_one_cluster(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            fasta, gff = _write_extract_fixture(root, _tier1_genes())
            extract_out = root / "extract"
            self.assertEqual(
                0,
                main(["extract", "--fasta", str(fasta), "--gff3", str(gff), "--interval", "chr1:1-1000", "--out-dir", str(extract_out)]),
            )
            providers = _write_provider_tables(root)
            tier0_out = root / "tier0"
            build_tier0_inventory(extract_out, tier0_out, providers, CASE_ID)
            frame = pd.read_csv(tier0_out / "tier0_inventory.tsv", sep="\t", dtype=str).fillna("")
            out = root / "tier1"
            manifest = build_tier1_clusters(frame.to_dict("records"), out, CASE_ID, providers=providers)
            clusters = _read_clusters(out / "tier1_clusters.tsv")
            self.assertEqual(13, manifest["gene_count"])
            self.assertEqual(13, manifest["genes_in_cluster_count"])
            members = [gene for cluster in clusters.values() for gene in cluster["member_gene_ids"].split(";")]
            self.assertEqual(13, len(members))
            self.assertEqual(13, len(set(members)))

    def test_mixed_missing_and_not_evaluated_keeps_provider_gap(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out = root / "extract"
            extract_out.mkdir()
            (extract_out / "candidates.tsv").write_text(
                "gene_id\tchrom\tstart\tend\tstrand\tgene_name\tdescription\tcds_length\n"
                "M1\tchr1\t1\t10\t+\t\tName1\t12\n"
                "M2\tchr1\t21\t30\t+\t\tName2\t12\n",
                encoding="utf-8",
            )
            (root / "domains.tsv").write_text("gene_id\tdomain_id\nM2\tPF001\n", encoding="utf-8")
            providers = {"domains": load_provider_table("domains", str(root / "domains.tsv"), "dom-db", "v1")}
            tier0_out = root / "tier0"
            build_tier0_inventory(extract_out, tier0_out, providers, CASE_ID)
            frame = pd.read_csv(tier0_out / "tier0_inventory.tsv", sep="\t", dtype=str).fillna("")
            out = root / "tier1"
            manifest = build_tier1_clusters(frame.to_dict("records"), out, CASE_ID, providers=providers)
            clusters = _read_clusters(out / "tier1_clusters.tsv")
            m1 = clusters[frozenset({"M1"})]
            self.assertEqual("model_state_rescue", m1["rule_key"])
            self.assertEqual("rescue", m1["allocation"])
            self.assertEqual("provider_gap", m1["rescue_bucket"])
            self.assertEqual("model_state_rescue:provider_gap", m1["allocation_reason"])
            m2 = clusters[frozenset({"M2"})]
            self.assertEqual("domain", m2["rule_key"])
            self.assertEqual({"priority": 1, "rescue": 1}, manifest["allocation_counts"])

    def test_not_evaluated_never_converted_to_negative(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out = root / "extract"
            extract_out.mkdir()
            (extract_out / "candidates.tsv").write_text(
                "gene_id\tchrom\tstart\tend\tstrand\tgene_name\tdescription\tcds_length\n"
                "M1\tchr1\t1\t10\t+\t\tName1\t12\n",
                encoding="utf-8",
            )
            (root / "domains.tsv").write_text("gene_id\tdomain_id\nM2\tPF001\n", encoding="utf-8")
            providers = {"domains": load_provider_table("domains", str(root / "domains.tsv"), "dom-db", "v1")}
            tier0_out = root / "tier0"
            build_tier0_inventory(extract_out, tier0_out, providers, CASE_ID)
            frame = pd.read_csv(tier0_out / "tier0_inventory.tsv", sep="\t", dtype=str).fillna("")
            rows = {row["gene_id"]: row for row in frame.to_dict("records")}
            self.assertEqual("missing", rows["M1"]["domain_status"])
            for cell in ("orthogroup_status", "expression_status", "variant_status", "sv_status", "pav_status", "promoter_status"):
                with self.subTest(cell=cell):
                    self.assertEqual("not_evaluated", rows["M1"][cell])
            out = root / "tier1"
            build_tier1_clusters(frame.to_dict("records"), out, CASE_ID, providers=providers)
            clusters = _read_clusters(out / "tier1_clusters.tsv")
            self.assertEqual("provider_gap", clusters[frozenset({"M1"})]["rescue_bucket"])
            # a gene with no providers at all (fully not_evaluated) is never
            # demoted into provider_gap: unchanged behavior.
            hand_rows = [_hand_row("N1")]
            out2 = root / "tier1b"
            build_tier1_clusters(hand_rows, out2, CASE_ID)
            clusters2 = _read_clusters(out2 / "tier1_clusters.tsv")
            self.assertEqual("singleton", clusters2[frozenset({"N1"})]["rule_key"])
            self.assertEqual("not_applicable", clusters2[frozenset({"N1"})]["rescue_bucket"])

    def test_ambiguity_preserved_in_representative_reason_and_members(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            rows = [
                _hand_row("C1", gene_name="GeneC", conflict_flags="domains:C1:domain_id;expression:C1:value"),
                _hand_row("C2", gene_name="GeneC"),
            ]
            out = root / "tier1"
            manifest = build_tier1_clusters(rows, out, CASE_ID)
            clusters = _read_clusters(out / "tier1_clusters.tsv")
            cluster = clusters[frozenset({"C1", "C2"})]
            self.assertEqual("C1", cluster["representative_gene_id"])
            self.assertIn("conflict_flags=domains:C1:domain_id;expression:C1:value", cluster["representative_reason"])
            self.assertIn(";member_conflicts=C1", cluster["representative_reason"])
            self.assertEqual("C1;C2", cluster["member_gene_ids"])
            self.assertEqual(2, manifest["genes_in_cluster_count"])

    def test_multi_family_gene_joins_smallest_cluster_with_union_family_ids(self) -> None:
        """Fix 7 (a): a gene with several family ids joins the cluster of its
        lexicographically smallest family id; family_ids is the union of every
        member's family values."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            family = self._family_table(root, "gene_id\tfamily_id\nG1\tFAM_B\nG1\tFAM_A\nG2\tFAM_A\n")
            rows = [_hand_row("G1"), _hand_row("G2")]
            out = root / "tier1"
            manifest = build_tier1_clusters(rows, out, CASE_ID, providers={"family": family})
            clusters = _read_clusters(out / "tier1_clusters.tsv")
            self.assertEqual(1, len(clusters))
            cluster = clusters[frozenset({"G1", "G2"})]
            self.assertEqual("family", cluster["rule_key"])
            self.assertEqual("FAM_A", cluster["cluster_key"])
            self.assertEqual("shared_family", cluster["allocation_reason"])
            self.assertEqual("2", cluster["member_count"])
            self.assertEqual("FAM_A;FAM_B", cluster["family_ids"])
            self.assertEqual(1, manifest["cluster_count"])
    def test_name_sharing_gene_of_family_consumed_gene_stays_singleton(self) -> None:
        """Fix 7 (b): a gene sharing a gene_name with a family-consumed gene is
        never pulled into an annotation cluster; it remains a singleton."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            family = self._family_table(root, "gene_id\tfamily_id\nG1\tFAM001\n")
            rows = [_hand_row("G1", gene_name="GeneX"), _hand_row("G2", gene_name="GeneX")]
            out = root / "tier1"
            build_tier1_clusters(rows, out, CASE_ID, providers={"family": family})
            clusters = _read_clusters(out / "tier1_clusters.tsv")
            self.assertEqual("family", clusters[frozenset({"G1"})]["rule_key"])
            g2 = clusters[frozenset({"G2"})]
            self.assertEqual("singleton", g2["rule_key"])
            self.assertEqual("defer", g2["allocation"])
            self.assertEqual("singleton:G2", g2["seed_reason"])

    def test_sv_pav_bucket_wins_over_missing_model_overlap(self) -> None:
        """Fix 7 (c): first-match rescue bucket priority -- a gene with both
        missing model and sv/pav evidence lands in sv_pav_regulatory."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            rows = [_hand_row("R1", model_status="no_cds", sv_status="present", sv_evidence="DEL")]
            out = root / "tier1"
            manifest = build_tier1_clusters(rows, out, CASE_ID)
            clusters = _read_clusters(out / "tier1_clusters.tsv")
            cluster = clusters[frozenset({"R1"})]
            self.assertEqual("model_state_rescue", cluster["rule_key"])
            self.assertEqual("sv_pav_regulatory", cluster["rescue_bucket"])
            self.assertEqual("model_state_rescue:sv_pav_regulatory", cluster["allocation_reason"])
            self.assertEqual({"sv_pav_regulatory": 1}, manifest["rescue_bucket_counts"])

    def test_header_only_provider_table_no_index_error_all_missing(self) -> None:
        """Fix 7 (d): a 0-row (header-only) provider table never raises
        IndexError anywhere; every gene stays missing."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out = root / "extract"
            extract_out.mkdir()
            (extract_out / "candidates.tsv").write_text(
                "gene_id\tchrom\tstart\tend\tstrand\tgene_name\tdescription\tcds_length\n"
                "H1\tchr1\t1\t10\t+\t\tName1\t12\n",
                encoding="utf-8",
            )
            (root / "expression.tsv").write_text("gene_id\tvalue\n", encoding="utf-8")
            table = load_provider_table("expression", str(root / "expression.tsv"), "expr-db", "v1")
            tier0_out = root / "tier0"
            manifest = build_tier0_inventory(extract_out, tier0_out, {"expression": table}, CASE_ID)
            frame = pd.read_csv(tier0_out / "tier0_inventory.tsv", sep="\t", dtype=str).fillna("")
            row = {r["gene_id"]: r for r in frame.to_dict("records")}["H1"]
            self.assertEqual("missing", row["expression_status"])
            self.assertEqual("", row["expression_evidence"])
            self.assertEqual("", manifest["provider_status"]["expression"]["value_column"])
            out = root / "tier1"
            manifest1 = build_tier1_clusters(frame.to_dict("records"), out, CASE_ID, providers={"expression": table})
            self.assertEqual(1, manifest1["gene_count"])

    def test_padded_evidence_values_no_phantom_conflict(self) -> None:
        """Fix 2 regression (Fix 7 e): whitespace-padded evidence values dedupe
        into one value and produce no conflict entry, end to end."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out = root / "extract"
            extract_out.mkdir()
            (extract_out / "candidates.tsv").write_text(
                "gene_id\tchrom\tstart\tend\tstrand\tgene_name\tdescription\tcds_length\n"
                "P1\tchr1\t1\t10\t+\t\tName1\t12\n",
                encoding="utf-8",
            )
            (root / "domains.tsv").write_text("gene_id\tdomain_id\nP1\t PF001 \nP1\tPF001\n", encoding="utf-8")
            providers = {"domains": load_provider_table("domains", str(root / "domains.tsv"), "dom-db", "v1")}
            tier0_out = root / "tier0"
            build_tier0_inventory(extract_out, tier0_out, providers, CASE_ID)
            frame = pd.read_csv(tier0_out / "tier0_inventory.tsv", sep="\t", dtype=str).fillna("")
            row = {r["gene_id"]: r for r in frame.to_dict("records")}["P1"]
            self.assertEqual("PF001", row["domain_hits"])
            self.assertEqual("", row["conflict_flags"])
            self.assertEqual("present", row["domain_status"])
            out = root / "tier1"
            build_tier1_clusters(frame.to_dict("records"), out, CASE_ID, providers=providers)
            clusters = _read_clusters(out / "tier1_clusters.tsv")
            cluster = clusters[frozenset({"P1"})]
            self.assertEqual("domain", cluster["rule_key"])
            self.assertEqual("PF001", cluster["cluster_key"])

    def test_family_ids_dedupe_padded_values(self) -> None:
        """Fix 2: whitespace-padded family values dedupe into one family id."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            family = self._family_table(root, "gene_id\tfamily_id\nG1\tFAM001 \nG1\tFAM001\nG2\tFAM001\n")
            rows = [_hand_row("G1"), _hand_row("G2")]
            out = root / "tier1"
            build_tier1_clusters(rows, out, CASE_ID, providers={"family": family})
            clusters = _read_clusters(out / "tier1_clusters.tsv")
            cluster = clusters[frozenset({"G1", "G2"})]
            self.assertEqual("FAM001", cluster["family_ids"])
            self.assertEqual("FAM001", cluster["cluster_key"])


if __name__ == "__main__":
    unittest.main()
