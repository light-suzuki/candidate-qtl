"""Issue #590: connect frozen ontology resolution to existing term propagation."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.term_propagation import propagate_terms


CANDIDATES = "gene_id\nPs1\n"
ORTHOLOGY = (
    "candidate_gene_id\treference_gene_id\trelationship\tconfidence\t"
    "reference_role\tsource_genome\treference_species\n"
    "Ps1\tREF1\tone_to_many\t0.9\tfunctional_reference\tmedtr.A17\tMedicago truncatula\n"
)


def _write(root: Path, name: str, text: str) -> Path:
    path = root / name
    path.write_text(text, encoding="utf-8")
    return path


def _resolution(rows: str) -> str:
    return (
        "trait_query_id\traw_text\trequested_ontology_id\trequested_ontology_namespace\t"
        "resolution_status\tontology_id\tontology_name\tontology_namespace\tmatch_basis\t"
        "candidate_ontology_ids\tparent_ids\tobsolete_source_id\t"
        "replacement_or_consider_ids\trank_contribution\n"
        + rows
    )


class OntologyTermPropagationTests(unittest.TestCase):
    def _run(self, root: Path, annotation: str, trait_terms: str, out_name: str) -> tuple[pd.DataFrame, dict]:
        candidates = _write(root, f"{out_name}_candidates.tsv", CANDIDATES)
        orthology = _write(root, f"{out_name}_orthology.tsv", ORTHOLOGY)
        annotations = _write(root, f"{out_name}_annotation.tsv", annotation)
        terms = _write(root, f"{out_name}_terms.tsv", trait_terms)
        result = propagate_terms(
            str(candidates), str(orthology), str(annotations), str(root / out_name),
            trait_terms_arg=str(terms), full_hash=True,
        )
        frame = pd.read_csv(result["propagated_terms"], sep="\t", dtype=str).fillna("")
        manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
        return frame, manifest

    def test_exact_ontology_id_matches_without_changing_orthology_confidence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            annotation = (
                "reference_gene_id\tterm\tontology_id\tontology_name\tsource\n"
                "REF1\tdifferent source wording\tTO:0001\tRoot architecture trait\tGFR\n"
            )
            resolution = _resolution(
                "q1\troot architecture\t\tTO\tresolved\tTO:0001\tRoot architecture trait\tTO\t"
                "exact_label\tTO:0001\tTO:PARENT\t\t\tnone\n"
            )
            frame, manifest = self._run(root, annotation, resolution, "id_match")
            self.assertEqual("trait_matched", frame.iloc[0].trait_relevance)
            self.assertEqual("0.6", frame.iloc[0].term_confidence)
            self.assertEqual("ontology_resolution", manifest["trait_vocabulary_mode"])
            self.assertEqual(1, manifest["resolved_trait_concept_count"])
            self.assertEqual(1, manifest["trait_ontology_id_count"])
            self.assertEqual(0, manifest["ambiguous_trait_canonical_name_count"])

    def test_canonical_name_fallback_matches_when_annotation_id_is_missing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            annotation = (
                "reference_gene_id\tterm\tontology_id\tontology_name\n"
                "REF1\tdifferent source wording\t\tRoot architecture trait\n"
            )
            resolution = _resolution(
                "q1\tuser wording\t\tTO\tresolved\tTO:0001\tRoot architecture trait\tTO\t"
                "exact_label\tTO:0001\t\t\t\tnone\n"
            )
            frame, _manifest = self._run(root, annotation, resolution, "name_match")
            self.assertEqual("trait_matched", frame.iloc[0].trait_relevance)

    def test_same_canonical_name_across_namespaces_disables_name_only_fallback(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            annotation = (
                "reference_gene_id\tterm\tontology_id\tontology_name\n"
                "REF1\tShared canonical label\t\tShared canonical label\n"
            )
            resolution = _resolution(
                "q1\tTO concept\tTO:0001\tTO\tresolved\tTO:0001\tShared canonical label\tTO\t"
                "exact_id\tTO:0001\t\t\t\tnone\n"
                "q2\tPO concept\tPO:0002\tPO\tresolved\tPO:0002\tShared canonical label\tPO\t"
                "exact_id\tPO:0002\t\t\t\tnone\n"
            )
            frame, manifest = self._run(root, annotation, resolution, "namespace_name_ambiguity")
            self.assertEqual("unrelated", frame.iloc[0].trait_relevance)
            self.assertEqual(1, manifest["ambiguous_trait_canonical_name_count"])
            self.assertEqual(2, manifest["trait_ontology_id_count"])

    def test_conflicting_annotation_id_cannot_be_overridden_by_matching_name(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            annotation = (
                "reference_gene_id\tterm\tontology_id\tontology_name\n"
                "REF1\tRoot architecture trait\tTO:OTHER\tRoot architecture trait\n"
            )
            resolution = _resolution(
                "q1\troot architecture\t\tTO\tresolved\tTO:0001\tRoot architecture trait\tTO\t"
                "exact_label\tTO:0001\t\t\t\tnone\n"
            )
            frame, _manifest = self._run(root, annotation, resolution, "id_conflict")
            self.assertEqual("unrelated", frame.iloc[0].trait_relevance)

    def test_same_term_with_distinct_ontology_ids_preserves_both_provenance_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            annotation = (
                "reference_gene_id\tterm\tontology_id\tontology_name\n"
                "REF1\tshared display term\tGO:OTHER\tShared display term\n"
                "REF1\tshared display term\tTO:0001\tShared display term\n"
            )
            resolution = _resolution(
                "q1\tshared display term\t\tTO\tresolved\tTO:0001\tShared display term\tTO\t"
                "exact_label\tTO:0001\t\t\t\tnone\n"
            )
            frame, manifest = self._run(root, annotation, resolution, "distinct_ids")
            self.assertEqual(2, len(frame))
            observed = set(zip(frame.ontology_id, frame.trait_relevance, strict=True))
            self.assertEqual(
                {("GO:OTHER", "unrelated"), ("TO:0001", "trait_matched")},
                observed,
            )
            self.assertEqual(1, manifest["trait_matched_terms"])

    def test_parent_and_ambiguous_concepts_are_not_promoted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            annotation = (
                "reference_gene_id\tterm\tontology_id\tontology_name\n"
                "REF1\tBroad parent trait\tTO:PARENT\tBroad parent trait\n"
                "REF1\tAmbiguous label\tTO:AMBIG\tAmbiguous label\n"
            )
            resolution = _resolution(
                "q1\tchild\t\tTO\tresolved\tTO:CHILD\tSpecific child trait\tTO\texact_label\tTO:CHILD\tTO:PARENT\t\t\tnone\n"
                "q2\tambiguous label\t\t\tambiguous\t\t\t\t\tTO:AMBIG;PATO:AMBIG\t\t\t\tnone\n"
            )
            frame, manifest = self._run(root, annotation, resolution, "context_only")
            self.assertEqual({"unrelated"}, set(frame.trait_relevance))
            self.assertEqual(1, manifest["resolved_trait_concept_count"])

    def test_zero_resolved_concepts_leave_relevance_not_evaluated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            annotation = "reference_gene_id\tterm\tontology_id\nREF1\tAnything\tTO:X\n"
            resolution = _resolution(
                "q1\tunknown\t\tTO\tunresolved\t\t\t\t\t\t\t\t\tnone\n"
            )
            frame, manifest = self._run(root, annotation, resolution, "none_resolved")
            self.assertEqual("not_evaluated", frame.iloc[0].trait_relevance)
            self.assertFalse(manifest["trait_terms_provided"])
            self.assertEqual("ontology_resolution", manifest["trait_vocabulary_mode"])

    def test_scoring_shaped_resolution_artifact_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = _write(root, "candidates.tsv", CANDIDATES)
            orthology = _write(root, "orthology.tsv", ORTHOLOGY)
            annotation = _write(root, "annotation.tsv", "reference_gene_id\tterm\nREF1\tRoot trait\n")
            resolution = _write(
                root,
                "resolution.tsv",
                _resolution(
                    "q1\troot\t\tTO\tresolved\tTO:1\tRoot trait\tTO\texact_label\tTO:1\t\t\t\tscore\n"
                ),
            )
            with self.assertRaisesRegex(ValueError, "rank_contribution=none"):
                propagate_terms(
                    str(candidates), str(orthology), str(annotation), str(root / "out"),
                    trait_terms_arg=str(resolution),
                )

    def test_legacy_plain_term_list_keeps_existing_substring_behavior(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            annotation = "reference_gene_id\tterm\nREF1\troot growth\n"
            frame, manifest = self._run(root, annotation, "term\nroot\n", "legacy")
            self.assertEqual("trait_matched", frame.iloc[0].trait_relevance)
            self.assertEqual("legacy_terms", manifest["trait_vocabulary_mode"])


if __name__ == "__main__":
    unittest.main()
