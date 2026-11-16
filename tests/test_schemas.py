from __future__ import annotations

import unittest

import pandas as pd

from candidate_qtl import schemas


def _evidence(rows: list[dict[str, str]]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=schemas.EVIDENCE_COLUMNS)


def _good_evidence_row(**overrides: str) -> dict[str, str]:
    row = {
        "gene_id": "GeneA",
        "evidence_code": "trait_function_match",
        "direction": "support",
        "confidence": "high",
        "value": "1",
        "source_type": "publication",
        "source_id": "PMID:1",
        "description": "ok",
    }
    row.update(overrides)
    return row


def _candidates(rows: list[dict[str, object]]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=schemas.CANDIDATE_COLUMNS)


def _good_candidate_row(**overrides: object) -> dict[str, object]:
    row = {col: "" for col in schemas.CANDIDATE_COLUMNS}
    row.update({"gene_id": "GeneA", "chrom": "chr1", "start": 100, "end": 200, "strand": "+"})
    row.update(overrides)
    return row


class EvidenceSchemaTests(unittest.TestCase):
    def test_valid_frame_passes(self) -> None:
        frame = _evidence([_good_evidence_row(), _good_evidence_row(gene_id="GeneB")])
        self.assertEqual(schemas.collect_errors(frame, schemas.EVIDENCE_SCHEMA), [])

    def test_invalid_direction_is_rejected(self) -> None:
        frame = _evidence([_good_evidence_row(direction="maybe")])
        with self.assertRaises(schemas.SchemaValidationError) as ctx:
            schemas.validate(frame, schemas.EVIDENCE_SCHEMA)
        codes = {err.code for err in ctx.exception.errors}
        self.assertIn("invalid_enum", codes)

    def test_invalid_confidence_is_rejected(self) -> None:
        frame = _evidence([_good_evidence_row(confidence="0.9")])
        errors = schemas.collect_errors(frame, schemas.EVIDENCE_SCHEMA)
        self.assertTrue(any(err.column == "confidence" and err.code == "invalid_enum" for err in errors))

    def test_empty_gene_id_is_rejected(self) -> None:
        frame = _evidence([_good_evidence_row(gene_id="  ")])
        errors = schemas.collect_errors(frame, schemas.EVIDENCE_SCHEMA)
        self.assertTrue(any(err.column == "gene_id" and err.code == "empty_value" for err in errors))

    def test_missing_column_reports_table_level(self) -> None:
        frame = _evidence([_good_evidence_row()]).drop(columns=["confidence"])
        errors = schemas.collect_errors(frame, schemas.EVIDENCE_SCHEMA)
        self.assertEqual([err.code for err in errors], ["missing_column"])

    def test_non_raising_returns_errors(self) -> None:
        frame = _evidence([_good_evidence_row(direction="bad")])
        errors = schemas.validate(frame, schemas.EVIDENCE_SCHEMA, raise_on_error=False)
        self.assertTrue(errors)


class CandidateSchemaTests(unittest.TestCase):
    def test_valid_candidates_pass(self) -> None:
        frame = _candidates([_good_candidate_row(), _good_candidate_row(gene_id="GeneB")])
        self.assertEqual(schemas.collect_errors(frame, schemas.CANDIDATE_SCHEMA), [])

    def test_duplicate_gene_id_is_rejected(self) -> None:
        frame = _candidates([_good_candidate_row(), _good_candidate_row(start=300, end=400)])
        errors = schemas.collect_errors(frame, schemas.CANDIDATE_SCHEMA)
        self.assertTrue(any(err.code == "duplicate_key" for err in errors))

    def test_start_after_end_is_rejected(self) -> None:
        frame = _candidates([_good_candidate_row(start=500, end=200)])
        errors = schemas.collect_errors(frame, schemas.CANDIDATE_SCHEMA)
        self.assertTrue(any(err.code == "invalid_coordinate" and "must be <=" in err.message for err in errors))

    def test_non_positive_coordinate_is_rejected(self) -> None:
        frame = _candidates([_good_candidate_row(start=0, end=10)])
        errors = schemas.collect_errors(frame, schemas.CANDIDATE_SCHEMA)
        self.assertTrue(any(err.code == "invalid_coordinate" for err in errors))

    def test_empty_seqid_is_rejected(self) -> None:
        frame = _candidates([_good_candidate_row(chrom="")])
        errors = schemas.collect_errors(frame, schemas.CANDIDATE_SCHEMA)
        self.assertTrue(any(err.column == "chrom" and err.code == "empty_value" for err in errors))

    def test_non_numeric_coordinate_is_rejected(self) -> None:
        frame = _candidates([_good_candidate_row(start="abc")])
        errors = schemas.collect_errors(frame, schemas.CANDIDATE_SCHEMA)
        self.assertTrue(any(err.code == "invalid_coordinate" for err in errors))


class WetSchemaTests(unittest.TestCase):
    def test_assay_and_candidate_states_are_separate_and_never_causal(self) -> None:
        self.assertIn("completed", schemas.ASSAY_STATES)
        self.assertIn("requires_independent_validation", schemas.CANDIDATE_WET_STATES)
        self.assertNotIn("causal", schemas.CANDIDATE_WET_STATES)

    def test_wet_assay_schema_rejects_unknown_state(self) -> None:
        row = {column: "" for column in schemas.WET_ASSAY_COLUMNS}
        row.update({"assay_id": "A1", "gene_id": "GeneA", "assay_type": "sanger", "assay_state": "causal"})
        errors = schemas.collect_errors(pd.DataFrame([row]), schemas.WET_ASSAY_SCHEMA)
        self.assertTrue(any(error.code == "invalid_enum" for error in errors))

class CopySetEvidenceSourceVocabularyTests(unittest.TestCase):
    def test_route_a_blast_columns_and_vocabulary_are_registered_and_shared(self) -> None:
        from candidate_qtl import family_to_copy
        self.assertEqual(["copy_set_evidence_source", "route_a_blast_guard"],
                         schemas.FAMILY_COPY_ROUTE_A_BLAST_COLUMNS)
        self.assertEqual(
            frozenset({"atlas_orthogroup", "atlas_curated_family", "atlas_domain", "route_a_blastp"}),
            schemas.COPY_SET_EVIDENCE_SOURCES)
        self.assertEqual(
            frozenset({"interval_paralog_present", "paralog_guard_not_evaluated"}),
            schemas.ROUTE_A_BLAST_GUARD_REASONS,
        )
        # family_to_copy must import the same constants, not its own copies.
        self.assertIs(family_to_copy.COPY_SET_EVIDENCE_SOURCES, schemas.COPY_SET_EVIDENCE_SOURCES)
        self.assertIs(family_to_copy.ROUTE_A_BLAST_GUARD_REASONS, schemas.ROUTE_A_BLAST_GUARD_REASONS)
        self.assertEqual(family_to_copy.ROUTE_A_BLAST_DISCRIMINATION_COLUMNS,
                         schemas.FAMILY_COPY_ROUTE_A_BLAST_COLUMNS)
        self.assertEqual(schemas.FAMILY_COPY_EVIDENCE_SOURCE_SCHEMA.name, "family_copy_evidence_source")
        self.assertEqual(schemas.FAMILY_COPY_EVIDENCE_SOURCE_SCHEMA.columns,
                         schemas.FAMILY_COPY_ROUTE_A_BLAST_COLUMNS)
        self.assertEqual(schemas.FAMILY_COPY_EVIDENCE_SOURCE_SCHEMA.enums["route_a_blast_guard"],
                         schemas.ROUTE_A_BLAST_GUARD_REASONS)

    def test_route_a_blast_guard_holds_on_canonical_vocabulary_only(self) -> None:
        from candidate_qtl import family_to_copy
        # A valid hold passes; an off-vocabulary reason is rejected fail-closed.
        family_to_copy._validate_route_a_blast_guard([{"route_a_blast_guard": "interval_paralog_present"}])
        family_to_copy._validate_route_a_blast_guard([{"route_a_blast_guard": "paralog_guard_not_evaluated"}])
        family_to_copy._validate_route_a_blast_guard([{"route_a_blast_guard": ""}])
        with self.assertRaises(schemas.SchemaValidationError):
            family_to_copy._validate_route_a_blast_guard([{"route_a_blast_guard": "bogus_hold"}])


class UseCaseTypeSchemaTests(unittest.TestCase):
    def test_use_case_types_vocabulary_registered(self) -> None:
        expected = {
            "type1_reference_or_id_mismatch",
            "type2_semantic_or_mechanism_mismatch",
            "type3_cross_species_transferable",
            "type4_species_specific_novel",
            "unclassified",
        }
        self.assertEqual(schemas.USE_CASE_TYPES, frozenset(expected))
        self.assertEqual(schemas.RETENTION_SCOPE_USE_CASE_TYPES, frozenset({
            "type1_reference_or_id_mismatch",
            "type2_semantic_or_mechanism_mismatch",
            "type3_cross_species_transferable",
        }))

    def test_require_use_case_type_validation(self) -> None:
        for uct in schemas.USE_CASE_TYPES:
            self.assertEqual(schemas.require_use_case_type(uct), uct)
        with self.assertRaises(ValueError):
            schemas.require_use_case_type("type5_unknown")
        with self.assertRaises(ValueError):
            schemas.require_use_case_type("")


class PanderaBackendTests(unittest.TestCase):
    @unittest.skipUnless(schemas.pandera_available(), "pandera extra not installed")
    def test_pandera_schema_matches_builtin(self) -> None:
        import pandera as pa

        pandera_schema = schemas.pandera_schema(schemas.EVIDENCE_SCHEMA)
        good = _evidence([_good_evidence_row()])
        pandera_schema.validate(good)  # should not raise
        bad = _evidence([_good_evidence_row(direction="nope")])
        with self.assertRaises(pa.errors.SchemaError):
            pandera_schema.validate(bad)


if __name__ == "__main__":
    unittest.main()
