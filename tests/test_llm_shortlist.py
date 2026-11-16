import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.hardening import _candidate_literature_query
from candidate_qtl.llm import export_shortlist_review_prompt, import_shortlist_review


class LlmShortlistTests(unittest.TestCase):
    def test_query_uses_homolog_not_local_gene_id(self) -> None:
        query = _candidate_literature_query([{
            "direction": "support",
            "value": "Mitogen kinase OS=Arabidopsis thaliana GN=MEKK1 PE=1",
        }], "er2 recessive powdery mildew resistance")
        self.assertIn("MEKK1", query)
        self.assertIn('"powdery mildew"', query)

    def test_frozen_set_round_trip(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            genes = ["G1", "G2"]
            packet = {"species": "pea", "trait": "trait", "candidates": [{"gene_id": gene} for gene in genes]}
            (root / "packet.json").write_text(json.dumps(packet), encoding="utf-8")
            export_shortlist_review_prompt(str(root / "packet.json"), str(root / "prompt.jsonl"))
            pd.DataFrame({"gene_id": genes, "deterministic_review_order": [1, 2]}).to_csv(
                root / "shortlist.tsv", sep="\t", index=False
            )
            ranking = [
                {"gene_id": "G2", "rank": 1, "rationale": "x", "source_ids": [], "wet_test": "sequence"},
                {"gene_id": "G1", "rank": 2, "rationale": "y", "source_ids": ["PMID:1"], "wet_test": "sequence"},
            ]
            response = {"candidates": [{"model_json": {"ranking": ranking}}]}
            (root / "responses.json").write_text(json.dumps(response), encoding="utf-8")
            result = import_shortlist_review(str(root / "shortlist.tsv"), str(root / "responses.json"), str(root / "out"))
            output = pd.read_csv(result["ranking"], sep="\t")
            self.assertEqual(["G2", "G1"], list(output["gene_id"]))
            self.assertFalse(output["reviewed"].any())

    def test_rejects_changed_candidate_set(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pd.DataFrame({"gene_id": ["G1", "G2"]}).to_csv(root / "shortlist.tsv", sep="\t", index=False)
            bad = {"candidates": [{"model_json": {"ranking": [
                {"gene_id": "G1", "rank": 1, "rationale": "x", "source_ids": [], "wet_test": "x"},
                {"gene_id": "G3", "rank": 2, "rationale": "x", "source_ids": [], "wet_test": "x"},
            ]}}]}
            (root / "bad.json").write_text(json.dumps(bad), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "exactly the frozen"):
                import_shortlist_review(str(root / "shortlist.tsv"), str(root / "bad.json"), str(root / "out"))
