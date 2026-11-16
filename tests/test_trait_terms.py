from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main


class TraitTermsBuildTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.terms = self.root / "terms.tsv"
        self.terms.write_text(
            "term\tterm_type\tsource_id\tapproved\n"
            "fasciation\tdirect_phenotype\tpaper1\ttrue\n"
            "meristem\torgan_stage\tpaper1\ttrue\n"
            "stem cell maintenance\tbiological_process\tpaper2\tfalse\n",
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def _build(self, *extra: str) -> pd.DataFrame:
        out = self.root / "approved.tsv"
        self.assertEqual(
            0,
            main(["trait-terms", "build", "--trait-terms", str(self.terms), "--out", str(out), *extra]),
        )
        return pd.read_csv(out, sep="\t").fillna("")

    def test_classification_and_approval_gate(self) -> None:
        approved = self._build()
        by_term = {row["term"]: row for _, row in approved.iterrows()}
        # Unapproved term is excluded from the approved export.
        self.assertNotIn("stem cell maintenance", by_term)
        self.assertEqual(by_term["fasciation"]["weight_class"], "direct")
        self.assertEqual(by_term["meristem"]["weight_class"], "context")

    def test_manual_synonyms_are_approved_and_llm_pending(self) -> None:
        synonyms = self.root / "syn.tsv"
        synonyms.write_text("term\tsynonym\tsource_id\nfasciation\tstem fasciation\tcurator\n", encoding="utf-8")
        proposals = self.root / "llm.json"
        proposals.write_text(
            json.dumps([{"term": "fasciation", "synonym": "band widening", "source_id": "gpt"}]),
            encoding="utf-8",
        )
        out = self.root / "approved.tsv"
        self.assertEqual(
            0,
            main(
                [
                    "trait-terms", "build", "--trait-terms", str(self.terms), "--out", str(out),
                    "--manual-synonyms", str(synonyms), "--llm-proposals", str(proposals),
                ]
            ),
        )
        approved = pd.read_csv(out, sep="\t").fillna("")
        self.assertIn("stem fasciation", set(approved["term"]))
        # A model-proposed synonym never reaches the approved export.
        self.assertNotIn("band widening", set(approved["term"]))

        audit = pd.read_csv(out.with_name("approved_audit.tsv"), sep="\t").fillna("")
        band = audit[audit["term"] == "band widening"].iloc[0]
        self.assertEqual(band["match_type"], "llm_suggested")
        self.assertEqual(band["approval_status"], "pending_review")

    def test_benchmark_mode_rejects_answer_leakage(self) -> None:
        leaky = self.root / "leaky.tsv"
        leaky.write_text(
            "term\tterm_type\tsource_id\tapproved\n"
            "fasciation\tdirect_phenotype\tpaper1\ttrue\n"
            "PsCLV3\tmechanism_hypothesis\tpaper1\ttrue\n",
            encoding="utf-8",
        )
        out = self.root / "bench.tsv"
        self.assertEqual(
            0,
            main(
                [
                    "trait-terms", "build", "--trait-terms", str(leaky), "--out", str(out),
                    "--benchmark", "--blocked-term", "PsCLV3",
                ]
            ),
        )
        approved = pd.read_csv(out, sep="\t").fillna("")
        self.assertIn("fasciation", set(approved["term"]))
        self.assertNotIn("PsCLV3", set(approved["term"]))
        audit = pd.read_csv(out.with_name("bench_audit.tsv"), sep="\t").fillna("")
        rejected = audit[audit["term"] == "PsCLV3"].iloc[0]
        self.assertTrue(rejected["exclusion_reason"].startswith("answer_leakage"))

    def test_ontology_ids_attached_from_template(self) -> None:
        repo = Path(__file__).resolve().parents[1]
        approved = self._build("--ontology-map", str(repo / "templates" / "trait_ontology_map.tsv"))
        by_term = {row["term"]: row for _, row in approved.iterrows()}
        self.assertEqual(by_term["fasciation"]["ontology_id"], "TO:0000208")
        self.assertEqual(by_term["meristem"]["ontology_id"], "PO:0009013")


if __name__ == "__main__":
    unittest.main()
