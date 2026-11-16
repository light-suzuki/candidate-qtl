from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8"); return path


class TermPropagationTests(unittest.TestCase):
    def _fixture(self, root: Path) -> tuple[Path, Path, Path]:
        candidates = _write(root / "candidates.tsv", "gene_id\nGENE_A\nGENE_B\nGENE_C\n")
        orthology = _write(root / "orthology.tsv",
            "candidate_gene_id\treference_gene_id\trelationship\tsource_genome\treference_role\treference_species\tconfidence\n"
            "GENE_A\trefX\tone_to_one\tAT_genome\tfunctional_reference\tArabidopsis\t\n"
            "GENE_B\trefY\tone_to_many\tanchor_genome\tclose_anchor\tMedicago\t\n"
            "GENE_A\trefZ\tortholog\tanchor_genome\tclose_anchor\tMedicago\t\n"
            "NOT_A_CANDIDATE\trefX\tone_to_one\tAT_genome\tfunctional_reference\tArabidopsis\t\n")
        annotation = _write(root / "ref_annotation.tsv",
            "reference_gene_id\tterm\tontology_id\tontology_name\tsource\n"
            "refX\tmeristem maintenance\tGO:0010073\tmeristem maintenance\tTAIR\n"
            "refX\tCLAVATA signaling\tGO:0009934\tregulation of meristem\tTAIR\n"
            "refY\tstem elongation\tGO:0009826\tunidimensional cell growth\tMtGDB\n")
        return candidates, orthology, annotation

    def test_propagates_reference_ortholog_terms_with_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); candidates, orthology, annotation = self._fixture(root); out = root / "out"
            self.assertEqual(0, main(["propagate-terms", "--candidates", str(candidates), "--orthology", str(orthology),
                                      "--reference-annotation", str(annotation), "--out-dir", str(out)]))
            terms = pd.read_csv(out / "propagated_terms.tsv", sep="\t", dtype=str).fillna("")
            a_terms = set(terms[terms["candidate_gene_id"].eq("GENE_A")]["term_normalized"])
            self.assertEqual({"meristem maintenance", "clavata signaling"}, a_terms)
            self.assertEqual({"stem elongation"}, set(terms[terms["candidate_gene_id"].eq("GENE_B")]["term_normalized"]))
            self.assertEqual({"ortholog_propagated"}, set(terms["match_type"]))
            self.assertEqual({"pending"}, set(terms["review_status"]))
            self.assertIn("GO:0010073", set(terms["ontology_id"]))
            # interval-scoped: the non-candidate ortholog row is ignored.
            self.assertNotIn("NOT_A_CANDIDATE", set(terms["candidate_gene_id"]))
            manifest = json.loads((out / "term_propagation_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(3, manifest["candidate_count"])
            self.assertEqual(2, manifest["candidates_with_propagated_terms"])
            self.assertIn("GENE_C", manifest["candidates_without_ortholog"])
            self.assertTrue(manifest["interval_scoped"])

    def test_trait_terms_tag_relevance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); candidates, orthology, annotation = self._fixture(root)
            traits = _write(root / "traits.tsv", "term\nmeristem maintenance\n")
            out = root / "out"
            self.assertEqual(0, main(["propagate-terms", "--candidates", str(candidates), "--orthology", str(orthology),
                                      "--reference-annotation", str(annotation), "--trait-terms", str(traits), "--out-dir", str(out)]))
            terms = pd.read_csv(out / "propagated_terms.tsv", sep="\t", dtype=str).fillna("").set_index("term_normalized")
            self.assertEqual("trait_matched", terms.loc["meristem maintenance", "trait_relevance"])
            self.assertEqual("unrelated", terms.loc["stem elongation", "trait_relevance"])

    def test_native_only_input_degrades_gracefully(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); candidates, _, _ = self._fixture(root); out = root / "out"
            self.assertEqual(0, main(["propagate-terms", "--candidates", str(candidates), "--out-dir", str(out)]))
            terms = pd.read_csv(out / "propagated_terms.tsv", sep="\t", dtype=str)
            self.assertEqual(0, len(terms))
            manifest = json.loads((out / "term_propagation_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(manifest["native_only"])
            self.assertIn("no_orthology_mapping_native_only", manifest["skip_reasons"])

    def test_ortholog_without_reference_function_is_separated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = _write(root / "candidates.tsv", "gene_id\nGENE_A\n")
            orthology = _write(root / "orthology.tsv",
                "candidate_gene_id\treference_gene_id\trelationship\n" "GENE_A\trefNoFunc\tone_to_one\n")
            annotation = _write(root / "ann.tsv", "reference_gene_id\tterm\nrefOther\tsomething\n")
            out = root / "out"
            self.assertEqual(0, main(["propagate-terms", "--candidates", str(candidates), "--orthology", str(orthology),
                                      "--reference-annotation", str(annotation), "--out-dir", str(out)]))
            self.assertEqual(0, len(pd.read_csv(out / "propagated_terms.tsv", sep="\t")))
            manifest = json.loads((out / "term_propagation_manifest.json").read_text(encoding="utf-8"))
            self.assertIn("GENE_A", manifest["orthologs_without_reference_function"])

    def test_byte_stable_regeneration(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); candidates, orthology, annotation = self._fixture(root)
            a, b = root / "a", root / "b"
            for target in (a, b):
                self.assertEqual(0, main(["propagate-terms", "--candidates", str(candidates), "--orthology", str(orthology),
                                          "--reference-annotation", str(annotation), "--out-dir", str(target)]))
            for name in ("propagated_terms.tsv", "search_term_expansion.tsv"):
                self.assertEqual((a / name).read_bytes(), (b / name).read_bytes(), name)


if __name__ == "__main__":
    unittest.main()
