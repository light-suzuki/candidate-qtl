from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.tier0_inventory import (
    PROVIDER_KINDS,
    TIER0_ANNOTATION_STATUSES,
    TIER0_APPLICABILITY,
    TIER0_COLUMNS,
    TIER0_EVIDENCE_STATUSES,
    TIER0_MODEL_STATUSES,
    build_tier0_inventory,
    load_provider_table,
)

CASE_ID = "tier0-case"
EXTRACT_COLUMNS = [*TIER0_COLUMNS, "cds_length"]
KIND_STATUS_CELLS = {
    "domains": "domain_status",
    "orthogroup": "orthogroup_status",
    "expression": "expression_status",
    "variant": "variant_status",
    "sv": "sv_status",
    "pav": "pav_status",
    "promoter": "promoter_status",
}
KIND_EVIDENCE_CELLS = {
    "domains": "domain_hits",
    "orthogroup": "orthogroup",
    "expression": "expression_evidence",
    "variant": "variant_evidence",
    "sv": "sv_evidence",
    "pav": "pav_evidence",
    "promoter": "promoter_evidence",
}
STATUS_COLUMNS = ["model_status", "annotation_status", "symbol_status", *KIND_STATUS_CELLS.values()]


def _tier0_genes():
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
    fasta, gff = _write_extract_fixture(root, _tier0_genes())
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
    (root / "variant.tsv").write_text("gene_id\teffect\tsource_id\nG5\tmissense\tV1\nG5\tsplice\tV2\n", encoding="utf-8")
    (root / "family.tsv").write_text("gene_id\tfamily_id\nG5\tFAM001\n", encoding="utf-8")
    return {
        "expression": load_provider_table("expression", str(root / "expression.tsv"), "expr-db", "v2"),
        "domains": load_provider_table("domains", str(root / "domains.tsv"), "dom-db", "v1"),
        "orthogroup": load_provider_table("orthogroup", str(root / "orthogroup.tsv"), "egg-db", "v1"),
        "variant": load_provider_table("variant", str(root / "variant.tsv"), "var-db", "v1"),
        "family": load_provider_table("family", str(root / "family.tsv"), "fam-db", "v1"),
    }


class Tier0InventoryTests(unittest.TestCase):
    def test_no_providers_all_statuses_not_evaluated(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out = _run_extract(root, self)
            out = root / "tier0"
            manifest = build_tier0_inventory(extract_out, out, {}, CASE_ID)
            frame = pd.read_csv(out / "tier0_inventory.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(EXTRACT_COLUMNS, list(frame.columns))
            self.assertEqual(["G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9"], list(frame["gene_id"]))
            for row in frame.to_dict("records"):
                for kind in ("domains", "orthogroup", "expression", "variant", "sv", "pav", "promoter"):
                    self.assertEqual("not_evaluated", row[KIND_STATUS_CELLS[kind]])
                    self.assertEqual("", row[KIND_EVIDENCE_CELLS[kind]])
                self.assertEqual("not_evaluated", row["provider_applicability"])
                self.assertEqual("", row["conflict_flags"])
                self.assertEqual("", row["variant_source_id"])
            self.assertEqual(9, manifest["gene_count"])
            self.assertEqual("1.0", manifest["schema_version"])
            self.assertEqual(CASE_ID, manifest["case_id"])
            self.assertEqual({}, manifest["provider_status"])
            self.assertEqual([], manifest["weak_input_fingerprints"])
            for column in KIND_STATUS_CELLS.values():
                self.assertEqual({"not_evaluated": 9}, manifest["counts"][column])
            self.assertEqual({"not_evaluated": 9}, manifest["counts"]["provider_applicability"])
            self.assertNotIn("timestamp", (out / "tier0_manifest.json").read_text(encoding="utf-8"))

    def test_with_providers_statuses_present_and_missing(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out = _run_extract(root, self)
            providers = _write_provider_tables(root)
            out = root / "tier0"
            manifest = build_tier0_inventory(extract_out, out, providers, CASE_ID)
            frame = pd.read_csv(out / "tier0_inventory.tsv", sep="\t", dtype=str).fillna("")
            rows = {row["gene_id"]: row for row in frame.to_dict("records")}
            g1, g2, g5 = rows["G1"], rows["G2"], rows["G5"]
            self.assertEqual("applicable", g1["provider_applicability"])
            for kind in ("expression", "domains", "orthogroup", "variant"):
                self.assertEqual("missing", g1[KIND_STATUS_CELLS[kind]])
                self.assertEqual("", g1[KIND_EVIDENCE_CELLS[kind]])
            self.assertEqual("present", g5["expression_status"])
            self.assertEqual("10;20", g5["expression_evidence"])
            self.assertEqual("expr-db@v2", g5["expression_provider"])
            self.assertEqual("present", g5["domain_status"])
            self.assertEqual("PF001;PF002", g5["domain_hits"])
            self.assertEqual("dom-db@v1", g5["domain_provider"])
            self.assertEqual("present", g5["orthogroup_status"])
            self.assertEqual("OG0001", g5["orthogroup"])
            self.assertEqual("present", g5["variant_status"])
            self.assertEqual("missense;splice", g5["variant_evidence"])
            self.assertEqual("var-db@v1", g5["variant_provider"])
            self.assertEqual("V1;V2", g5["variant_source_id"])
            for kind in ("sv", "pav", "promoter"):
                self.assertEqual("not_evaluated", g5[f"{kind}_status"])
                self.assertEqual("", g5[f"{kind}_evidence"])
            self.assertEqual("present", g2["expression_status"])
            self.assertEqual("5", g2["expression_evidence"])
            self.assertEqual("missing", g1["expression_status"])
            for row in frame.to_dict("records"):
                for column in STATUS_COLUMNS:
                    self.assertNotIn(row[column], ("", "0"))
                self.assertNotIn(row["provider_applicability"], ("", "0"))
            self.assertEqual(9, manifest["gene_count"])

    def test_conflicts_and_sorted_distinct_values(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out = _run_extract(root, self)
            providers = _write_provider_tables(root)
            out = root / "tier0"
            build_tier0_inventory(extract_out, out, providers, CASE_ID)
            frame = pd.read_csv(out / "tier0_inventory.tsv", sep="\t", dtype=str).fillna("")
            rows = {row["gene_id"]: row for row in frame.to_dict("records")}
            self.assertEqual("domains:G5:domain_id;expression:G5:value;variant:G5:effect;variant:G5:source_id", rows["G5"]["conflict_flags"])
            self.assertEqual("", rows["G1"]["conflict_flags"])
            self.assertEqual("", rows["G2"]["conflict_flags"])
            self.assertEqual("10;20", rows["G5"]["expression_evidence"])
            self.assertEqual("PF001;PF002", rows["G5"]["domain_hits"])
            self.assertEqual("missense;splice", rows["G5"]["variant_evidence"])
            self.assertEqual("V1;V2", rows["G5"]["variant_source_id"])
            self.assertEqual(1, providers["expression"]["conflicts"])
            self.assertEqual(1, providers["domains"]["conflicts"])
            self.assertEqual(0, providers["orthogroup"]["conflicts"])

    def test_source_id_only_populated_when_column_present(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out = _run_extract(root, self)
            (root / "variant.tsv").write_text("gene_id\teffect\tsource_id\nG2\tmissense\tV9\n", encoding="utf-8")
            (root / "sv.tsv").write_text("gene_id\tsv_type\nG2\tDEL\n", encoding="utf-8")
            (root / "pav.tsv").write_text("gene_id\tpav_type\tsource_id\nG2\tABS\tP9\n", encoding="utf-8")
            (root / "promoter.tsv").write_text("gene_id\tmotif\nG2\tTATA\n", encoding="utf-8")
            providers = {
                "variant": load_provider_table("variant", str(root / "variant.tsv"), "var-db", "v1"),
                "sv": load_provider_table("sv", str(root / "sv.tsv"), "sv-db", "v1"),
                "pav": load_provider_table("pav", str(root / "pav.tsv"), "pav-db", "v1"),
                "promoter": load_provider_table("promoter", str(root / "promoter.tsv"), "prom-db", "v1"),
            }
            out = root / "tier0"
            build_tier0_inventory(extract_out, out, providers, CASE_ID)
            frame = pd.read_csv(out / "tier0_inventory.tsv", sep="\t", dtype=str).fillna("")
            g2 = {row["gene_id"]: row for row in frame.to_dict("records")}["G2"]
            self.assertEqual("V9", g2["variant_source_id"])
            self.assertEqual("", g2["sv_source_id"])
            self.assertEqual("P9", g2["pav_source_id"])
            self.assertEqual("", g2["promoter_source_id"])
            self.assertEqual("DEL", g2["sv_evidence"])
            self.assertEqual("TATA", g2["promoter_evidence"])
            self.assertEqual("sv-db@v1", g2["sv_provider"])

    def test_annotation_semantics_no_advantage(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out = _run_extract(root, self)
            out = root / "tier0"
            build_tier0_inventory(extract_out, out, {}, CASE_ID)
            frame = pd.read_csv(out / "tier0_inventory.tsv", sep="\t", dtype=str).fillna("")
            rows = {row["gene_id"]: row for row in frame.to_dict("records")}
            self.assertEqual("present", rows["G1"]["symbol_status"])
            self.assertEqual("present", rows["G6"]["symbol_status"])
            self.assertEqual("present", rows["G1"]["annotation_status"])
            self.assertEqual("partial", rows["G6"]["annotation_status"])
            self.assertEqual("missing", rows["G3"]["symbol_status"])
            self.assertEqual("no_cds", rows["G3"]["model_status"])
            self.assertEqual("partial_or_unverified_orf", rows["G4"]["model_status"])
            self.assertEqual("complete_orf", rows["G1"]["model_status"])
            self.assertEqual("", rows["G3"]["variant_source_id"])

    def test_manifest_keys_fingerprints_provider_status_counts(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out = _run_extract(root, self)
            providers = _write_provider_tables(root)
            out = root / "tier0"
            manifest = build_tier0_inventory(extract_out, out, providers, CASE_ID)
            self.assertEqual(
                ["candidates_columns_observed", "case_id", "counts", "gene_count", "input_fingerprints", "provider_status", "schema_version", "weak_input_fingerprints"],
                sorted(manifest),
            )
            self.assertEqual(
                ["candidate_proteins", "candidates", "domains", "expression", "extract_manifest", "family", "orthogroup", "transcripts", "variant"],
                sorted(manifest["input_fingerprints"]),
            )
            for key in ("candidates", "transcripts", "candidate_proteins", "extract_manifest"):
                fingerprint = manifest["input_fingerprints"][key]
                self.assertIn("path", fingerprint)
                self.assertEqual("strong", fingerprint["strength"])
                self.assertEqual("full_sha256", fingerprint["method"])
            self.assertEqual([], manifest["weak_input_fingerprints"])
            self.assertEqual(
                {"provided": True, "provider": "expr-db", "version": "v2", "rows": 3, "genes_covered": 2, "conflicts": 1, "value_column": "value"},
                manifest["provider_status"]["expression"],
            )
            self.assertEqual(
                {"provided": True, "provider": "dom-db", "version": "v1", "rows": 2, "genes_covered": 1, "conflicts": 1, "value_column": "domain_id"},
                manifest["provider_status"]["domains"],
            )
            self.assertEqual(
                {"provided": True, "provider": "fam-db", "version": "v1", "rows": 1, "genes_covered": 1, "conflicts": 0, "value_column": "family_id"},
                manifest["provider_status"]["family"],
            )
            self.assertEqual({"complete_orf": 7, "no_cds": 1, "partial_or_unverified_orf": 1}, manifest["counts"]["model_status"])
            self.assertEqual({"partial": 6, "present": 3}, manifest["counts"]["annotation_status"])
            self.assertEqual({"missing": 3, "present": 6}, manifest["counts"]["symbol_status"])
            self.assertEqual({"missing": 7, "present": 2}, manifest["counts"]["expression_status"])
            self.assertEqual({"missing": 8, "present": 1}, manifest["counts"]["domain_status"])
            self.assertEqual({"missing": 8, "present": 1}, manifest["counts"]["orthogroup_status"])
            self.assertEqual({"missing": 8, "present": 1}, manifest["counts"]["variant_status"])
            self.assertEqual({"not_evaluated": 9}, manifest["counts"]["sv_status"])
            self.assertEqual({"not_evaluated": 9}, manifest["counts"]["pav_status"])
            self.assertEqual({"not_evaluated": 9}, manifest["counts"]["promoter_status"])
            self.assertEqual({"applicable": 9}, manifest["counts"]["provider_applicability"])
            self.assertEqual(25, len(manifest["candidates_columns_observed"]))
            self.assertEqual(manifest["candidates_columns_observed"], sorted(manifest["candidates_columns_observed"]))

    def test_byte_identical_double_run(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out = _run_extract(root, self)
            providers = _write_provider_tables(root)
            first = root / "first"
            second = root / "second"
            build_tier0_inventory(extract_out, first, providers, CASE_ID)
            build_tier0_inventory(extract_out, second, providers, CASE_ID)
            self.assertEqual((first / "tier0_inventory.tsv").read_bytes(), (second / "tier0_inventory.tsv").read_bytes())
            self.assertEqual((first / "tier0_manifest.json").read_bytes(), (second / "tier0_manifest.json").read_bytes())

    def test_no_file_writes_outside_out_dir(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out = _run_extract(root, self)
            out = root / "tier0"
            build_tier0_inventory(extract_out, out, {}, CASE_ID)
            self.assertEqual({"tier0_inventory.tsv", "tier0_manifest.json"}, {path.name for path in out.iterdir()})

    def test_duplicate_gene_id_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out = root / "extract"
            extract_out.mkdir()
            (extract_out / "candidates.tsv").write_text(
                "gene_id\tchrom\tstart\tend\tstrand\tgene_name\tcds_length\n"
                "DUP1\tchr1\t1\t10\t+\tGeneD\t12\n"
                "DUP2\tchr1\t21\t30\t+\t\t12\n"
                "DUP1\tchr1\t41\t50\t+\t\t12\n",
                encoding="utf-8",
            )
            with self.assertRaises(ValueError) as caught:
                build_tier0_inventory(extract_out, root / "tier0", {}, CASE_ID)
            self.assertIn("DUP1", str(caught.exception))
            self.assertIn("duplicate gene_id", str(caught.exception))
            self.assertFalse((root / "tier0" / "tier0_inventory.tsv").exists())

    def test_load_provider_table_shape_and_errors(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            table = root / "expression.tsv"
            table.write_text("gene_id\tvalue\tsource_id\nG1\t10\tE1\nG1\t10\tE2\n", encoding="utf-8")
            loaded = load_provider_table("expression", str(table), "db", "v3")
            self.assertEqual("expression", loaded["kind"])
            self.assertEqual("db", loaded["provider"])
            self.assertEqual("v3", loaded["version"])
            self.assertEqual(str(table), loaded["path"])
            self.assertEqual(["G1"], loaded["gene_id_set"])
            self.assertEqual(1, loaded["conflicts"])
            self.assertEqual("10", loaded["rows"][0]["value"])
            self.assertIsInstance(loaded["rows"][0]["value"], str)
            with self.assertRaises(ValueError):
                load_provider_table("bogus", str(table), "db", "v1")
            with self.assertRaises(ValueError):
                load_provider_table("expression", str(root / "missing.tsv"), "db", "v1")
            no_gene_id = root / "no_gene_id.tsv"
            no_gene_id.write_text("value\tsource_id\n10\tE1\n", encoding="utf-8")
            with self.assertRaises(ValueError) as caught:
                load_provider_table("expression", str(no_gene_id), "db", "v1")
            self.assertIn("gene_id", str(caught.exception))
            no_value_column = root / "no_value.tsv"
            no_value_column.write_text("gene_id\tfoo\nG1\tx\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                load_provider_table("expression", str(no_value_column), "db", "v1")
            self.assertEqual(PROVIDER_KINDS, {"expression", "domains", "orthogroup", "family", "variant", "sv", "pav", "promoter"})
            for kind in ("expression", "domains", "orthogroup", "family", "variant", "sv", "pav", "promoter"):
                for status in ("not_evaluated",):
                    self.assertIn(status, TIER0_EVIDENCE_STATUSES)
            self.assertEqual({"complete_orf", "partial_or_unverified_orf", "no_cds", "not_evaluated"}, TIER0_MODEL_STATUSES)
            self.assertEqual({"present", "partial", "missing", "not_evaluated"}, TIER0_ANNOTATION_STATUSES)
            self.assertEqual({"applicable", "not_evaluated"}, TIER0_APPLICABILITY)


class Tier0EmptyValueTests(unittest.TestCase):
    """Audit finding F: empty/whitespace-only evidence values are never
    positive; a row present in a provided table only counts as evidence when
    its kind's value cell is non-empty after stripping."""

    def test_empty_evidence_value_is_missing_not_present(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out = _run_extract(root, self)
            (root / "expression.tsv").write_text("gene_id\tvalue\nG1\t\nG5\t10\n", encoding="utf-8")
            (root / "domains.tsv").write_text("gene_id\tdomain_id\nG1\t  \nG5\tPF001\n", encoding="utf-8")
            providers = {
                "expression": load_provider_table("expression", str(root / "expression.tsv"), "expr-db", "v1"),
                "domains": load_provider_table("domains", str(root / "domains.tsv"), "dom-db", "v1"),
            }
            out = root / "tier0"
            manifest = build_tier0_inventory(extract_out, out, providers, CASE_ID)
            frame = pd.read_csv(out / "tier0_inventory.tsv", sep="\t", dtype=str).fillna("")
            rows = {row["gene_id"]: row for row in frame.to_dict("records")}
            self.assertEqual("present", rows["G5"]["expression_status"])
            self.assertEqual("present", rows["G5"]["domain_status"])
            for kind in ("expression", "domains"):
                with self.subTest(kind=kind):
                    status_cell = KIND_STATUS_CELLS[kind]
                    evidence_cell = KIND_EVIDENCE_CELLS[kind]
                    self.assertEqual("missing", rows["G1"][status_cell])
                    self.assertEqual("", rows["G1"][evidence_cell])
                    self.assertEqual({"missing": 8, "present": 1}, manifest["counts"][status_cell])
                    self.assertEqual(1, manifest["provider_status"][kind]["genes_covered"])

    def test_source_id_alone_never_positive(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out = _run_extract(root, self)
            (root / "variant.tsv").write_text("gene_id\teffect\tsource_id\nG1\t\tV9\n", encoding="utf-8")
            providers = {"variant": load_provider_table("variant", str(root / "variant.tsv"), "var-db", "v1")}
            out = root / "tier0"
            manifest = build_tier0_inventory(extract_out, out, providers, CASE_ID)
            frame = pd.read_csv(out / "tier0_inventory.tsv", sep="\t", dtype=str).fillna("")
            g1 = {row["gene_id"]: row for row in frame.to_dict("records")}["G1"]
            self.assertEqual("missing", g1["variant_status"])
            self.assertEqual("", g1["variant_evidence"])
            self.assertEqual("V9", g1["variant_source_id"])
            self.assertEqual(0, manifest["provider_status"]["variant"]["genes_covered"])
            self.assertEqual({"missing": 9}, manifest["counts"]["variant_status"])

    def test_kind_not_provided_stays_not_evaluated(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out = _run_extract(root, self)
            (root / "expression.tsv").write_text("gene_id\tvalue\nG1\t\nG5\t10\n", encoding="utf-8")
            providers = {"expression": load_provider_table("expression", str(root / "expression.tsv"), "expr-db", "v1")}
            out = root / "tier0"
            build_tier0_inventory(extract_out, out, providers, CASE_ID)
            frame = pd.read_csv(out / "tier0_inventory.tsv", sep="\t", dtype=str).fillna("")
            rows = {row["gene_id"]: row for row in frame.to_dict("records")}
            for kind in ("sv", "pav", "promoter", "orthogroup"):
                with self.subTest(kind=kind):
                    self.assertEqual("not_evaluated", rows["G1"][KIND_STATUS_CELLS[kind]])
                    self.assertEqual("", rows["G1"][KIND_EVIDENCE_CELLS[kind]])

    def test_genes_covered_excludes_empty_value_only_genes(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out = _run_extract(root, self)
            (root / "expression.tsv").write_text("gene_id\tvalue\nG1\t\nG2\t   \nG5\t10\n", encoding="utf-8")
            providers = {"expression": load_provider_table("expression", str(root / "expression.tsv"), "expr-db", "v1")}
            out = root / "tier0"
            manifest = build_tier0_inventory(extract_out, out, providers, CASE_ID)
            self.assertEqual(3, manifest["provider_status"]["expression"]["rows"])
            self.assertEqual(1, manifest["provider_status"]["expression"]["genes_covered"])
            self.assertEqual({"missing": 8, "present": 1}, manifest["counts"]["expression_status"])

    def test_padded_evidence_values_stripped_deduped_no_phantom_conflict(self) -> None:
        """Fix 2: " 10 " and "10" are the same evidence value; whitespace
        padding must never create a distinct set entry or a phantom conflict."""
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            extract_out = _run_extract(root, self)
            (root / "expression.tsv").write_text("gene_id\tvalue\nG5\t 10 \nG5\t10\nG2\t  5 \n", encoding="utf-8")
            providers = {"expression": load_provider_table("expression", str(root / "expression.tsv"), "expr-db", "v1")}
            out = root / "tier0"
            manifest = build_tier0_inventory(extract_out, out, providers, CASE_ID)
            frame = pd.read_csv(out / "tier0_inventory.tsv", sep="\t", dtype=str).fillna("")
            rows = {row["gene_id"]: row for row in frame.to_dict("records")}
            self.assertEqual("10", rows["G5"]["expression_evidence"])
            self.assertEqual("5", rows["G2"]["expression_evidence"])
            self.assertEqual("", rows["G5"]["conflict_flags"])
            self.assertEqual("", rows["G2"]["conflict_flags"])
            self.assertEqual(0, providers["expression"]["conflicts"])
            self.assertEqual("present", rows["G5"]["expression_status"])
            self.assertEqual("present", rows["G2"]["expression_status"])
            self.assertEqual({"missing": 7, "present": 2}, manifest["counts"]["expression_status"])
            self.assertEqual("value", manifest["provider_status"]["expression"]["value_column"])

    def test_header_only_provider_table_succeeds_all_missing(self) -> None:
        """Fix 7 (d): a 0-row (header-only) provider table loads, never raises
        IndexError, leaves every gene missing, and records value_column ""."""
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
            self.assertEqual([], table["rows"])
            out = root / "tier0"
            manifest = build_tier0_inventory(extract_out, out, {"expression": table}, CASE_ID)
            frame = pd.read_csv(out / "tier0_inventory.tsv", sep="\t", dtype=str).fillna("")
            row = frame.iloc[0]
            self.assertEqual("missing", row["expression_status"])
            self.assertEqual("", row["expression_evidence"])
            self.assertEqual("", manifest["provider_status"]["expression"]["value_column"])
            self.assertEqual(0, manifest["provider_status"]["expression"]["genes_covered"])


if __name__ == "__main__":
    unittest.main()
