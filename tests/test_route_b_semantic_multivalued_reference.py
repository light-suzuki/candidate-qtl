"""A one-to-many family names several reference genes, and must still ground (issue #483).

The bundle hands the model ``orthology_reference_gene`` as a semicolon-joined
field, so a paralog-rich family arrives as
``sp|P93003|TFL1_ARATH;sp|Q9ZNV5|CEN_ARATH``. A model obeying "use only the
supplied unit" copies it back whole. Matching that joined string against the
frozen roster as a single token found nothing, so every one-to-many family --
exactly the ones this lane exists to disambiguate -- was rejected as
``reference_gene_not_in_frozen_bundle`` before any citation was looked at.

Grounding is still a closed-world exact-membership test. These tests pin that:
tokens the bundle never supplied ground nothing, and a claim whose tokens are
all foreign stays unsupported.
"""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main

JOINED = "sp|P93003|TFL1_ARATH;sp|Q9ZNV5|CEN_ARATH"
UNIT_ID = "bsem-multi"
PROMPT_ID = "route-b-semantic-bsem-multi"


def _bundle() -> dict:
    return {
        "schema_version": "1.0", "case_id": "multi", "route": "interval_first",
        "stage": "route_b_semantic_bundle", "phenotype": "determinate inflorescence",
        "interval_gene_count": 1, "semantic_unit_count": 1,
        "candidate_identity_preserved": True, "policy": "non_scoring", "b0_coverage": {},
        "units": [{
            "unit_id": UNIT_ID, "unit_type": "family", "unit_anchor": "TFL",
            "identity": {"family_ids": ["TFL"], "orthogroup_ids": [], "subfamily_ids": [],
                         "evolutionary_family_ids": []},
            "member_interval_genes": ["G_TFL"],
            "genes": [{"gene_id": "G_TFL", "orthology_reference_gene": JOINED,
                       "orthology_relationship": "one_to_many"}],
            "literature_link_candidates": [],
            "literature_query_status": [],
            "stable_source_ids": [],
            "rank_contribution": "none",
        }],
    }


def _prompt() -> dict:
    return {"case_id": "multi", "gene_id": "__route_b_semantic_unit__", "prompt_id": PROMPT_ID,
            "route_b_unit_id": UNIT_ID, "system": "system", "user": "user"}


def _claim(reference_gene: str) -> dict:
    return {"claim_type": "family_function", "biological_function": "meristem identity",
            "mechanism": "represses floral transition", "phenotype_relation": "plausible",
            "reference_gene": reference_gene, "reference_species": "Arabidopsis thaliana",
            "evidence_scope": "family", "source_ids": [],
            "alternative_hypothesis": "another family member carries the effect",
            "uncertainty_note": "copy-level assignment unresolved"}


class MultiValuedReferenceGeneTests(unittest.TestCase):
    def _import(self, reference_gene: str) -> pd.DataFrame:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bundle = root / "bundle.json"
            bundle.write_text(json.dumps(_bundle()), encoding="utf-8")
            prompts = root / "prompts.jsonl"
            prompts.write_text(json.dumps(_prompt()) + "\n", encoding="utf-8")
            responses = root / "responses.json"
            responses.write_text(json.dumps({"candidates": [{
                "prompt_id": PROMPT_ID, "parse_status": "json",
                "model_json": {"claims": [_claim(reference_gene)]},
            }]}), encoding="utf-8")
            out = root / "import"
            self.assertEqual(0, main(["route-b-semantic-import", "--bundle", str(bundle),
                                      "--prompts", str(prompts), "--responses", str(responses),
                                      "--out-dir", str(out), "--case-id", "multi"]))
            return pd.read_csv(out / "semantic_claim_proposals.tsv", sep="\t", dtype=str).fillna("")

    def test_joined_field_echoed_verbatim_is_grounded(self) -> None:
        proposals = self._import(JOINED)
        row = proposals.iloc[0]
        self.assertEqual("proposed", row["verification_status"])
        self.assertIn("reference_gene_grounded:", row["verification_reason"])

    def test_one_token_of_the_pair_is_grounded(self) -> None:
        proposals = self._import("sp|Q9ZNV5|CEN_ARATH")
        row = proposals.iloc[0]
        self.assertEqual("proposed", row["verification_status"])
        self.assertIn("sp|Q9ZNV5|CEN_ARATH", row["verification_reason"])

    def test_a_foreign_token_does_not_ground_the_claim(self) -> None:
        # Closed world: a gene the bundle never supplied cannot ground anything.
        proposals = self._import("sp|Q9SXZ2|FT_ARATH")
        row = proposals.iloc[0]
        self.assertEqual("unsupported", row["verification_status"])
        self.assertIn("reference_gene_not_in_frozen_bundle", row["verification_reason"])

    def test_mixing_a_foreign_token_in_reports_only_the_supplied_one(self) -> None:
        proposals = self._import("sp|P93003|TFL1_ARATH;sp|Q9SXZ2|FT_ARATH")
        row = proposals.iloc[0]
        self.assertEqual("proposed", row["verification_status"])
        self.assertIn("sp|P93003|TFL1_ARATH", row["verification_reason"])
        # The unsupplied one must not be presented as grounded.
        grounded = row["verification_reason"].split("reference_gene_grounded:", 1)[1].split(";")[0]
        self.assertNotIn("FT_ARATH", grounded)

    def test_empty_reference_gene_is_not_a_grounding_failure(self) -> None:
        proposals = self._import("")
        self.assertEqual("proposed", proposals.iloc[0]["verification_status"])


if __name__ == "__main__":
    unittest.main()
