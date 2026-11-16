"""Unit tests for the legacy 7-column mechanism payload adapter (issue #380)."""

from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from candidate_qtl.mechanism_payload_adapter import (
    ADAPTER_VERSION,
    FORBIDDEN_TOKENS,
    LEGACY_V1_7COL_COLUMNS,
    OUTPUT_SCHEMA_VERSION,
    TEMPLATE_24_COLUMNS,
    convert_legacy_v1_7col,
    detect_payload_schema,
    forbidden_columns,
    normalize_source_ids,
    validate_payload_file,
    write_adapter_map,
)
from candidate_qtl.phenotype_investigation import MECHANISM_COLUMNS

LEGACY_HEADER = "\t".join(LEGACY_V1_7COL_COLUMNS)
PEA_A_ROWS = [
    ("pea_a", "pea_a_M1", "Chalcone synthase (CHS) first committed step of anthocyanin biosynthesis",
     "corpus_ortholog_model_plant", "source_log.tsv#pea_a_M1", "verified",
     "Flower colour and seed-coat anthocyanin phenotype text maps to anthocyanin biosynthesis; CHS is the committed first step (corpus: At CHS, Mt CHS-class)."),
    ("pea_a", "pea_a_M2", "Late anthocyanin pathway enzymes (DFR / LDOX/ANS)",
     "corpus_ortholog_model_plant", "source_log.tsv#pea_a_M2", "verified",
     "Anthocyanin pigmentation loss in flower/seed coat follows late-pathway enzyme defects; DFR and LDOX corpus-verified (At)."),
    ("pea_a", "pea_a_M3", "MYB-bHLH-WD40 (MBW) transcriptional regulation of anthocyanin genes",
     "corpus_ortholog_model_plant", "source_log.tsv#pea_a_M3", "verified",
     "Flower colour and seed coat anthocyanin are classically controlled by MBW complexes; TT2/TT8/TTG1 entries corpus-verified in Arabidopsis and Lotus (LjTT2a-c/LjTT8/LjTTG1)."),
]
TWELVE_COL_SUBSET = ["pathway_or_process", "gene_family", "mechanism_statement", "source_ids", "counter_hypothesis",
                     "discriminating_observation", "model_species_evidence", "model_species_mutant_phenotype",
                     "proposed_stable_id", "reference_species", "reference_gene_name", "selection_rationale"]


def _legacy_payload(root: Path, name: str = "mechanism_payload.tsv", rows: list[tuple[str, ...]] | None = None) -> Path:
    path = root / name
    body = rows if rows is not None else PEA_A_ROWS
    path.write_text("\n".join([LEGACY_HEADER, *["\t".join(row) for row in body]]) + "\n", encoding="utf-8")
    return path


def _records() -> list[dict[str, str]]:
    return [dict(zip(LEGACY_V1_7COL_COLUMNS, row)) for row in PEA_A_ROWS]


class SchemaDetectionTests(unittest.TestCase):
    def test_adapter_constants_are_stable(self) -> None:
        self.assertEqual("mechanism-payload-adapter-1.0.0", ADAPTER_VERSION)
        self.assertEqual("mechanism_hypotheses_29col_v2", OUTPUT_SCHEMA_VERSION)
        self.assertEqual(9, len(FORBIDDEN_TOKENS))
        self.assertEqual(24, len(TEMPLATE_24_COLUMNS))
        self.assertTrue(TEMPLATE_24_COLUMNS <= set(MECHANISM_COLUMNS))

    def test_detect_legacy_7col_exact_header(self) -> None:
        self.assertEqual("legacy_v1_7col", detect_payload_schema(list(LEGACY_V1_7COL_COLUMNS)))
        self.assertEqual("legacy_v1_7col", detect_payload_schema(list(reversed(LEGACY_V1_7COL_COLUMNS))))

    def test_detect_canonical_29_24_subset(self) -> None:
        self.assertEqual("canonical_29_col", detect_payload_schema(list(MECHANISM_COLUMNS)))
        self.assertEqual("canonical_29_col", detect_payload_schema(sorted(MECHANISM_COLUMNS)))
        self.assertEqual("canonical_24_col_template", detect_payload_schema(sorted(TEMPLATE_24_COLUMNS)))
        self.assertEqual("canonical_subset", detect_payload_schema(TWELVE_COL_SUBSET))
        self.assertEqual("canonical_subset", detect_payload_schema(["case_id"]))
        self.assertEqual("unknown", detect_payload_schema([]))
        with tempfile.TemporaryDirectory() as tmp:
            payload = Path(tmp) / "canonical.tsv"
            payload.write_text("\t".join(TWELVE_COL_SUBSET)
                               + "\nstarch\tSBE\tstatement\tPMID:1\tcounter\tdiscriminating\tevidence\tphenotype\tP30924\tPisum\tsbei\trationale\n",
                               encoding="utf-8")
            result = validate_payload_file(payload, None, "pea_a")
            self.assertEqual("canonical_subset", result["schema"])
            self.assertEqual(1, result["rows_read"])
            self.assertEqual(64, len(result["raw_payload_sha256"]))

    def test_detect_unknown_reports_missing_and_unknown_columns(self) -> None:
        header = ["case_id", "pathway_or_process", "mystery_column"]
        self.assertEqual("unknown", detect_payload_schema(header))
        with tempfile.TemporaryDirectory() as tmp:
            payload = Path(tmp) / "payload.tsv"
            payload.write_text("\t".join(header) + "\npea_a\tstarch\tx\n", encoding="utf-8")
            with self.assertRaises(ValueError) as context:
                validate_payload_file(payload, None, "pea_a")
            message = str(context.exception)
            self.assertIn("schema=unknown", message)
            self.assertIn("unknown column(s): [mystery_column]", message)
            self.assertIn("missing canonical column(s): [", message)
            self.assertIn("hypothesis_id", message)
            self.assertNotIn("forbidden column", message)

    def test_detect_duplicate_header_columns(self) -> None:
        header = ["case_id", "case_id", "pathway_or_process"]
        self.assertEqual("unknown", detect_payload_schema(header))
        with tempfile.TemporaryDirectory() as tmp:
            payload = Path(tmp) / "payload.tsv"
            payload.write_text("\t".join(header) + "\npea_a\tpea_a\tstarch\n", encoding="utf-8")
            with self.assertRaises(ValueError) as context:
                validate_payload_file(payload, None, "pea_a")
            self.assertIn("duplicate column(s): [case_id, case_id]", str(context.exception))

    def test_forbidden_columns_never_match_canonical_names(self) -> None:
        self.assertEqual(["truth_gene", "candidate_rank"],
                         forbidden_columns(["rank_contribution", "interval_informed", "truth_gene", "candidate_rank"]))
        self.assertEqual([], forbidden_columns(list(MECHANISM_COLUMNS)))
        with tempfile.TemporaryDirectory() as tmp:
            payload = Path(tmp) / "payload.tsv"
            payload.write_text("case_id\ttruth_gene\npea_a\tx\n", encoding="utf-8")
            with self.assertRaises(ValueError) as context:
                validate_payload_file(payload, None, "pea_a")
            message = str(context.exception)
            self.assertIn("schema=unknown", message)
            self.assertIn("forbidden column 'truth_gene' (leakage class: truth|candidate|rank|score|interval|accession|bait); "
                          "forbidden columns are never accepted", message)

    def test_source_ids_normalize_losslessly_and_keep_string_compatibility(self) -> None:
        self.assertEqual(
            "PMID:1;DOI:10.1000/example;PMID:1",
            normalize_source_ids(["PMID:1", "DOI:10.1000/example", "PMID:1"]),
        )
        self.assertEqual("PMID:1; DOI:10.1000/example",
                         normalize_source_ids("PMID:1; DOI:10.1000/example"))

    def test_source_ids_reject_nested_objects_and_non_string_members(self) -> None:
        for value in ([["PMID:1"]], [{"id": "PMID:1"}], ["PMID:1", 2], ["PMID:1", None]):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, r"must be a string or list\[str\]"):
                    normalize_source_ids(value)

    def test_json_source_ids_shape_is_validated_before_import(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = root / "invalid_source_ids.json"
            payload.write_text(json.dumps({
                "hypotheses": [{"source_ids": {"citation": "PMID:1"}}],
            }), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, r"source_ids must be a string or list\[str\]"):
                validate_payload_file(payload, None, "")


class LegacyConversionTests(unittest.TestCase):
    def test_legacy_row_validation_missing_ids_and_duplicates(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            missing_id = _legacy_payload(root, "missing_id.tsv", [
                PEA_A_ROWS[0],
                ("pea_a", "", "no id", "corpus_ortholog_model_plant", "source_log.tsv#x", "verified", "r"),
            ])
            with self.assertRaises(ValueError) as context:
                validate_payload_file(missing_id, None, "pea_a")
            self.assertIn("row 2: mechanism_id is empty", str(context.exception))
            missing_case = _legacy_payload(root, "missing_case.tsv", [
                PEA_A_ROWS[0],
                ("", "pea_a_M9", "no case", "corpus_ortholog_model_plant", "source_log.tsv#x", "verified", "r"),
            ])
            with self.assertRaises(ValueError) as context:
                validate_payload_file(missing_case, None, "pea_a")
            self.assertIn("row 2: case_id is empty", str(context.exception))
            duplicate = _legacy_payload(root, "duplicate.tsv", [PEA_A_ROWS[0], PEA_A_ROWS[1], PEA_A_ROWS[1]])
            with self.assertRaises(ValueError) as context:
                validate_payload_file(duplicate, None, "pea_a")
            self.assertIn("row 3: duplicate mechanism_id 'pea_a_M2' (first at row 2)", str(context.exception))

    def test_legacy_conversion_maps_only_known_fields(self) -> None:
        records = _records()
        canonical, mapping = convert_legacy_v1_7col(records, case_id="pea_a", payload_source="llm")
        self.assertEqual(3, len(canonical))
        self.assertEqual(3, len(mapping))
        self.assertEqual(set(MECHANISM_COLUMNS), set(canonical[0]))
        first = canonical[0]
        self.assertEqual("pea_a", first["case_id"])
        self.assertEqual(records[0]["mechanism_name"], first["mechanism_statement"])
        self.assertEqual(records[0]["source_log_entry"], first["source_ids"])
        self.assertEqual(records[0]["selection_reason"], first["selection_rationale"])
        for column in ("pathway_or_process", "gene_family", "expected_role", "direction_consistency",
                       "model_species_evidence", "model_species_mutant_phenotype", "proposed_stable_id",
                       "reference_species", "reference_gene_name", "counter_hypothesis",
                       "discriminating_observation", "uncertainty", "model", "retrieval_provenance", "notes"):
            self.assertEqual("", first[column], column)
        self.assertEqual("unknown", first["confidence"])
        self.assertEqual("", first["hypothesis_id"])
        self.assertEqual("pending_verification", first["review_status"])
        self.assertEqual("none", first["rank_contribution"])
        self.assertEqual(1, first["proposal_round"])
        self.assertFalse(first["interval_informed"])
        self.assertEqual("pea_a", canonical[2]["case_id"])
        self.assertEqual(records[2]["mechanism_name"], canonical[2]["mechanism_statement"])

    def test_legacy_conversion_rejects_requested_case_mismatch(self) -> None:
        with self.assertRaises(ValueError) as context:
            convert_legacy_v1_7col(_records(), case_id="pea_b")
        self.assertIn("requested case_id 'pea_b'", str(context.exception))
        self.assertIn("payload case_id 'pea_a'", str(context.exception))

    def test_adapter_map_written_losslessly(self) -> None:
        _, mapping = convert_legacy_v1_7col(_records(), case_id="pea_a")
        with tempfile.TemporaryDirectory() as tmp:
            path = write_adapter_map(mapping, Path(tmp))
            self.assertEqual(Path(tmp) / "mechanism_payload_adapter_map.tsv", path)
            lines = [line.split("\t") for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
            self.assertEqual(["hypothesis_id", "source_row", "mechanism_id", "mechanism_name",
                              "evidence_class", "source_log_entry", "status"], lines[0])
            self.assertEqual(4, len(lines))
            for line, row in zip(lines[1:], mapping):
                self.assertEqual([row["hypothesis_id"], str(row["source_row"]), row["mechanism_id"],
                                  row["mechanism_name"], row["evidence_class"], row["source_log_entry"],
                                  row["status"]], line)
            self.assertEqual("pea_a_M1", lines[1][2])
            self.assertEqual("source_log.tsv#pea_a_M1", lines[1][5])
            self.assertEqual("verified", lines[1][6])
            self.assertEqual("mech_003", lines[3][0])
            self.assertEqual("3", lines[3][1])


class ValidationAndFreezeTests(unittest.TestCase):
    def test_valid_legacy_payload_returns_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            payload = _legacy_payload(Path(tmp))
            result = validate_payload_file(payload, None, "pea_a")
            self.assertEqual("legacy_v1_7col", result["schema"])
            self.assertEqual(3, result["rows_read"])
            self.assertEqual(hashlib.sha256(payload.read_bytes()).hexdigest(), result["raw_payload_sha256"])
            self.assertIsNone(result["freeze_verification"])
            self.assertEqual([], result["errors"])

    def test_freeze_manifest_match_records_verification(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = _legacy_payload(root)
            digest = hashlib.sha256(payload.read_bytes()).hexdigest()
            freeze = root / "freeze.tsv"
            freeze.write_text(f"case_id\tmechanism_payload_sha256\tnotes\npea_a\t{digest}\textra\n", encoding="utf-8")
            result = validate_payload_file(payload, freeze, "pea_a")
            self.assertEqual({"source": str(freeze), "expected_sha256": digest, "actual_sha256": digest, "matched": True},
                             result["freeze_verification"])

    def test_freeze_manifest_uses_unique_payload_case_when_unspecified(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = _legacy_payload(root)
            digest = hashlib.sha256(payload.read_bytes()).hexdigest()
            freeze = root / "freeze.tsv"
            freeze.write_text(f"case_id\tmechanism_payload_sha256\npea_a\t{digest}\n", encoding="utf-8")
            result = validate_payload_file(payload, freeze, "")
            self.assertEqual("pea_a", result["payload_case_id"])
            self.assertTrue(result["freeze_verification"]["matched"])

    def test_payload_case_mismatch_and_multiple_cases_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = _legacy_payload(root)
            with self.assertRaises(ValueError) as context:
                validate_payload_file(payload, None, "pea_b")
            self.assertIn("requested case_id 'pea_b'", str(context.exception))
            mixed = _legacy_payload(root, "mixed.tsv", [PEA_A_ROWS[0], ("pea_b", *PEA_A_ROWS[1][1:])])
            with self.assertRaises(ValueError) as context:
                validate_payload_file(mixed, None, "")
            self.assertIn("multiple case_id values", str(context.exception))

    def test_json_forbidden_and_unknown_keys_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            forbidden = root / "forbidden.json"
            forbidden.write_text(json.dumps({"hypotheses": [{"case_id": "pea_a"}], "truth_gene": "Psat1"}), encoding="utf-8")
            with self.assertRaises(ValueError) as context:
                validate_payload_file(forbidden, None, "")
            self.assertIn("forbidden JSON key 'truth_gene'", str(context.exception))
            unknown = root / "unknown.json"
            unknown.write_text(json.dumps({"hypotheses": [{"case_id": "pea_a", "mystery": "x"}]}), encoding="utf-8")
            with self.assertRaises(ValueError) as context:
                validate_payload_file(unknown, None, "")
            self.assertIn("unknown JSON key(s)", str(context.exception))

    def test_freeze_manifest_mismatch_and_missing_case_fail_closed(self) -> None:
        zeros = "0" * 64
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = _legacy_payload(root)
            freeze = root / "freeze.tsv"
            freeze.write_text(f"case_id\tmechanism_payload_sha256\npea_a\t{zeros}\n", encoding="utf-8")
            with self.assertRaises(ValueError) as context:
                validate_payload_file(payload, freeze, "pea_a")
            message = str(context.exception)
            self.assertIn("mechanism payload freeze mismatch: ", message)
            self.assertIn(f"expected sha256 {zeros} from ", message)
            self.assertIn("; actual ", message)
            missing = root / "freeze_missing_case.tsv"
            missing.write_text(f"case_id\tmechanism_payload_sha256\npea_b\t{zeros}\n", encoding="utf-8")
            with self.assertRaises(ValueError) as context:
                validate_payload_file(payload, missing, "pea_a")
            self.assertIn("mechanism payload freeze mismatch: ", str(context.exception))
            self.assertIn("has no row for case_id 'pea_a'", str(context.exception))

    def test_empty_and_missing_payload_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaises(ValueError) as context:
                validate_payload_file(root / "missing.tsv", None, "pea_a")
            self.assertIn("file not found or unreadable", str(context.exception))
            empty = _legacy_payload(root, "empty.tsv", [])
            empty.write_bytes(b"")
            with self.assertRaises(ValueError) as context:
                validate_payload_file(empty, None, "pea_a")
            self.assertIn("file is empty", str(context.exception))
            headers_only = _legacy_payload(root, "headers_only.tsv", [])
            result = validate_payload_file(headers_only, None, "pea_a")
            self.assertEqual("legacy_v1_7col", result["schema"])
            self.assertEqual(0, result["rows_read"])


if __name__ == "__main__":
    unittest.main()
