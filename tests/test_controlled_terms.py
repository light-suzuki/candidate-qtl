import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.controlled_terms import promote_controlled_terms


class ControlledTermTests(unittest.TestCase):
    def test_promotes_exact_approved_term_without_inventing_source(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pd.DataFrame([{
                "gene_id": "G1", "evidence_code": "homology_supported", "direction": "support",
                "confidence": "low", "value": "Transcription factor TCP4 GN=TCP4",
                "source_type": "blastp", "source_id": "hits.tsv", "description": "best hit",
            }]).to_csv(root / "evidence.tsv", sep="\t", index=False)
            result = promote_controlled_terms(str(root / "evidence.tsv"), str(root / "out"), ["TCP4"])
            output = pd.read_csv(result["evidence"], sep="\t")
            promoted = output[output["evidence_code"].eq("trait_function_match")].iloc[0]
            self.assertEqual("low", promoted["confidence"])
            self.assertEqual("hits.tsv", promoted["source_id"])
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertFalse(manifest["llm_used"])

    def test_rejects_too_short_broad_term(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "evidence.tsv"
            path.write_text("gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "at least three"):
                promote_controlled_terms(str(path), str(Path(tmp) / "out"), ["R"])

    def test_demotes_false_generic_resistance_match(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pd.DataFrame([{
                "gene_id": "G1", "evidence_code": "trait_function_match", "direction": "support",
                "confidence": "medium", "value": "Camphor resistance CrcB family protein GN=FEX",
                "source_type": "blastp", "source_id": "hits.tsv", "description": "generic resistance token",
            }]).to_csv(root / "evidence.tsv", sep="\t", index=False)
            result = promote_controlled_terms(
                str(root / "evidence.tsv"), str(root / "out"), ["MLO", "MEKK1"], True
            )
            output = pd.read_csv(result["evidence"], sep="\t")
            self.assertEqual("trait_context", output.iloc[0]["evidence_code"])
            self.assertEqual("context", output.iloc[0]["direction"])
