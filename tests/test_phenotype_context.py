from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from candidate_qtl.phenotype_context import (
    PHENOTYPE_CONTEXT_COLUMNS,
    context_frame,
    load_phenotype_context,
    normalize_phenotype_context,
    structured_context_block,
    structured_context_status,
)


class PhenotypeContextTests(unittest.TestCase):
    def _write_context(self, root: Path, content: str) -> Path:
        path = root / "phenotype_context.yaml"
        path.write_text(content, encoding="utf-8")
        return path
    def test_legacy_phenotype_text_is_backward_compatible(self) -> None:
        context = load_phenotype_context(None, case_id="legacy", phenotype_text="seed dormancy")
        self.assertEqual("legacy", context["case_id"])
        self.assertEqual("seed dormancy", context["phenotype_text"])
        self.assertEqual("raw_text", context["term_usage"])
        self.assertEqual("unknown", context["developmental_stage"])
        self.assertEqual("unknown", context["environment_or_treatment"])
        status = structured_context_status(context)
        self.assertEqual("structured_context_missing", status["status"])
        self.assertEqual("none", status["present_fields"])
        self.assertEqual("none", status["rank_contribution"])

    def test_nested_yaml_tracks_normalized_ontology_term_and_context(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "phenotype_context.yaml"
            path.write_text(
                "case_id: flowering_case\n"
                "trait:\n"
                "  phenotype_text: delayed flowering\n"
                "phenotype_context:\n"
                "  trait_ontology_id: TO:0002616\n"
                "  trait_ontology_term: flowering time trait\n"
                "  developmental_stage: reproductive transition\n"
                "  tissue: shoot apical meristem\n"
                "  environment_or_treatment: long day\n"
                "  measurement_method: days to first flower\n"
                "  comparison_or_control: wild type\n"
                "  effect_direction: increase\n"
                "  replicate_info: 3 biological replicates\n"
                "  source_id: PMID:fixture\n"
                "  source_detail: Figure 2\n",
                encoding="utf-8",
            )
            context = load_phenotype_context(path)
        self.assertEqual("flowering_case", context["case_id"])
        self.assertEqual("delayed flowering", context["phenotype_text"])
        self.assertEqual("TO:0002616", context["trait_ontology_id"])
        self.assertEqual("normalized_term", context["term_usage"])
        status = structured_context_status(context)
        self.assertEqual("structured_context_used", status["status"])
        self.assertIn("developmental_stage", status["present_fields"])
        self.assertIn("environment_or_treatment", status["present_fields"])

    def test_unmapped_term_is_preserved_without_forced_normalization(self) -> None:
        context = normalize_phenotype_context(
            {
                "case_id": "c1",
                "phenotype_text": "abnormal top cluster",
                "unmapped_term": "頂端花束",
                "developmental_stage": "adult reproductive",
            }
        )
        self.assertEqual("頂端花束", context["unmapped_term"])
        self.assertEqual("unknown", context["trait_ontology_id"])
        self.assertEqual("unknown", context["trait_ontology_term"])
        self.assertEqual("unmapped_term", context["term_usage"])

    def test_same_free_text_different_stage_or_environment_has_distinct_fingerprint(self) -> None:
        base = {"case_id": "case_a", "phenotype_text": "reduced branching"}
        seedling = normalize_phenotype_context(
            {**base, "developmental_stage": "seedling", "environment_or_treatment": "control"}
        )
        adult = normalize_phenotype_context(
            {**base, "developmental_stage": "adult", "environment_or_treatment": "drought"}
        )
        self.assertEqual(seedling["phenotype_text"], adult["phenotype_text"])
        self.assertNotEqual(seedling["context_fingerprint"], adult["context_fingerprint"])

    def test_fingerprint_is_deterministic_and_ignores_source_provenance(self) -> None:
        first = normalize_phenotype_context(
            {"phenotype_text": "short stem", "tissue": "stem", "source_id": "PMID:1"}
        )
        second = normalize_phenotype_context(
            {"phenotype_text": "short stem", "tissue": "stem", "source_id": "PMID:2"}
        )
        self.assertEqual(first["context_fingerprint"], second["context_fingerprint"])

    def test_case_mismatch_fails_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "does not match"):
            normalize_phenotype_context({"case_id": "case_a"}, case_id="case_b")

    def test_context_frame_has_canonical_columns(self) -> None:
        frame = context_frame({"case_id": "c1", "phenotype_text": "leaf angle", "tissue": "leaf"})
        self.assertEqual(PHENOTYPE_CONTEXT_COLUMNS, list(frame.columns))
        self.assertEqual("raw_text", frame.iloc[0]["term_usage"])

    def test_non_mapping_fields_never_imply_a_normalized_term(self) -> None:
        context = normalize_phenotype_context(
            {
                "phenotype_text": "pod wall thickness",
                "measurement_method": "caliper",
                "comparison_or_control": "thin-wall parent",
                "replicate_info": "n=20",
            }
        )
        self.assertEqual("raw_text", context["term_usage"])
        self.assertEqual("unknown", context["trait_ontology_id"])
        self.assertEqual("unknown", context["trait_ontology_term"])

    def test_top_level_case_id_mismatch_fails_closed_before_normalization(self) -> None:
        """P1: a YAML root declaring case_id: case_a must not silently attach
        case_a's context to a caller passing case_id="case_b"."""
        with tempfile.TemporaryDirectory() as tmp:
            path = self._write_context(
                Path(tmp),
                "case_id: case_a\n"
                "phenotype_context:\n"
                "  tissue: leaf\n"
                "  developmental_stage: seedling\n",
            )
            with self.assertRaisesRegex(ValueError, "does not match"):
                load_phenotype_context(path, case_id="case_b")

    def test_top_level_case_id_matching_the_caller_is_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = self._write_context(
                Path(tmp),
                "case_id: case_a\n"
                "phenotype_context:\n"
                "  tissue: leaf\n",
            )
            context = load_phenotype_context(path, case_id="case_a")
        self.assertEqual("case_a", context["case_id"])
        self.assertEqual("leaf", context["tissue"])

    def test_nested_case_id_mismatch_is_still_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = self._write_context(
                Path(tmp),
                "case_id: case_a\n"
                "phenotype_context:\n"
                "  case_id: case_a\n"
                "  tissue: leaf\n",
            )
            with self.assertRaisesRegex(ValueError, "does not match"):
                load_phenotype_context(path, case_id="case_b")

    def test_falsy_yaml_roots_are_rejected_not_silently_emptied(self) -> None:
        """P2: [], false, 0 and \"\" must raise, never become {} via `or {}`."""
        for content in ("[]\n", "false\n", "0\n", '""\n'):
            with self.subTest(content=content.strip()):
                with tempfile.TemporaryDirectory() as tmp:
                    path = self._write_context(Path(tmp), content)
                    with self.assertRaisesRegex(ValueError, "must contain an object"):
                        load_phenotype_context(path, case_id="c1")

    def test_empty_document_and_empty_mapping_are_the_intended_empty_context(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            empty = self._write_context(Path(tmp), "")
            context = load_phenotype_context(empty, case_id="c1", phenotype_text="seed dormancy")
        self.assertEqual("seed dormancy", context["phenotype_text"])
        self.assertEqual("raw_text", context["term_usage"])
        self.assertEqual("structured_context_missing", structured_context_status(context)["status"])
        context = normalize_phenotype_context({}, case_id="c1")
        self.assertEqual("unknown", context["phenotype_text"])

    def test_non_mapping_normalize_input_is_rejected(self) -> None:
        for payload in ([], False, 0, ""):
            with self.subTest(payload=payload):
                with self.assertRaisesRegex(ValueError, "must be a mapping"):
                    normalize_phenotype_context(payload)  # type: ignore[arg-type]

    def test_missing_context_file_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "not found"):
                load_phenotype_context(Path(tmp) / "missing.yaml", case_id="c1")

    def test_pandas_and_numpy_missing_scalars_normalize_to_unknown(self) -> None:
        context = normalize_phenotype_context({
            "phenotype_text": "dwarf plant",
            "developmental_stage": pd.NA,
            "tissue": pd.NaT,
            "environment_or_treatment": np.nan,
            "measurement_method": np.float64("nan"),
        })
        for field in ("developmental_stage", "tissue", "environment_or_treatment", "measurement_method"):
            self.assertEqual("unknown", context[field], field)
        self.assertEqual("dwarf plant", context["phenotype_text"])
        self.assertNotIn("<NA>", str(context))
        self.assertNotIn("NaT", str(context))

    def test_array_like_values_are_never_run_through_pd_isna(self) -> None:
        context = normalize_phenotype_context({
            "phenotype_text": "seed size",
            "measurement_method": ["caliper", "scale"],
            "replicate_info": np.array(["rep1", "rep2"]),
            "comparison_or_control": {"parent": "thin-wall"},
        })
        self.assertEqual("['caliper', 'scale']", context["measurement_method"])
        self.assertEqual("['rep1' 'rep2']", context["replicate_info"])
        self.assertEqual("{'parent': 'thin-wall'}", context["comparison_or_control"])
        self.assertEqual("raw_text", context["term_usage"])

    def test_structured_context_block_is_manifest_ready(self) -> None:
        block = structured_context_block(None, case_id="c1", phenotype_text="dwarf plant")
        self.assertEqual("structured_context_missing", block["status"])
        self.assertEqual("raw_text", block["term_usage"])
        self.assertEqual("dwarf plant", block["phenotype_text"])
        self.assertEqual("unknown", block["unmapped_term"])
        self.assertEqual("none", block["present_fields"])
        self.assertEqual("none", block["rank_contribution"])
        self.assertEqual("stored_only", block["query_usage"]["usage"])
        self.assertEqual(64, len(block["context_fingerprint"]))

    def test_structured_context_block_records_used_context_and_honest_query_usage(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = self._write_context(
                Path(tmp),
                "phenotype_context:\n"
                "  trait_ontology_id: TO:0002616\n"
                "  developmental_stage: reproductive transition\n",
            )
            block = structured_context_block(
                path, case_id="c1", phenotype_text="delayed flowering",
                query_usage={"usage": "stored_only",
                             "reason": "the investigation prompt is built from the free text exactly as before"},
            )
        self.assertEqual("structured_context_used", block["status"])
        self.assertEqual("normalized_term", block["term_usage"])
        self.assertIn("developmental_stage", block["present_fields"])
        self.assertEqual("delayed flowering", block["phenotype_text"])
        self.assertEqual("stored_only", block["query_usage"]["usage"])
        self.assertEqual("none", block["rank_contribution"])


if __name__ == "__main__":
    unittest.main()
