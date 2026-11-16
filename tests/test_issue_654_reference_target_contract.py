"""Synthetic contract tests for non-focal reference targets (#654)."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.phenotype_investigation import (
    mechanism_import,
    phenotype_investigate,
)


def _portfolio() -> list[dict[str, object]]:
    layers = ["perception", "signaling", "regulation", "execution", "alternative", "unclear"]
    rows: list[dict[str, object]] = []
    for index in range(8):
        rows.append(
            {
                "causal_layer": layers[index % len(layers)],
                "pathway_or_process": f"synthetic_process_{index}",
                "gene_family": f"SYNTHETIC_FAMILY_{index}",
                "mechanism_statement": "a lesion in this family could alter the observed developmental transition",
                "expected_role": "regulatory input",
                "direction_consistency": "consistent",
                "model_species_evidence": "a named model-plant mutant study",
                "model_species_mutant_phenotype": "altered developmental timing",
                "source_ids": f"PMID:{10000000 + index}",
                "proposed_stable_id": "P12345",
                "reference_species": "Arabidopsis thaliana",
                "reference_gene_name": f"REF{index}",
                "selection_rationale": "characterised non-focal reference target",
                "counter_hypothesis": "the phenotype may arise through an independent tissue pathway",
                "discriminating_observation": "measure the transition in a staged tissue series",
                "confidence": "medium",
                "uncertainty": "copy-level transfer remains unresolved",
            }
        )
    return rows


class ReferenceTargetPromptTests(unittest.TestCase):
    def test_prompt_requires_real_reference_target_and_forbids_fabrication(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            result = phenotype_investigate(
                tmp,
                "case-654",
                "Pisum sativum",
                "late transition from vegetative to reproductive development",
            )
            prompt = Path(result["prompt"]).read_text(encoding="utf-8")
            self.assertIn("one real, named reference gene or locus", prompt)
            self.assertIn("focal-species primary reference is allowed", prompt)
            self.assertIn("strict benchmark route", prompt)
            self.assertIn("primary approach is to bridge", prompt)
            self.assertIn("Focal-species direct reports and other-species mechanistic comparators", prompt)
            self.assertIn("Do not invent genes, loci, accessions, species, or source ids", prompt)
            self.assertIn("missing target is review-required", prompt)
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            contract = manifest["payload_contract"]["reference_target_contract"]
            self.assertFalse(contract["target_must_be_non_focal"])
            self.assertTrue(contract["provider_failure_is_separate"])
            self.assertIn("verified_reference_species_relation", contract["portfolio_ready_row_requires"])

    def test_strict_prompt_contains_policy_and_matches_manifest_hash(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            result = phenotype_investigate(
                tmp, "case-654", "Pisum sativum", "late transition", benchmark_mode=True
            )
            prompt = Path(result["prompt"]).read_text(encoding="utf-8")
            self.assertIn("STRICT BENCHMARK POLICY", prompt)
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            import hashlib
            self.assertEqual(manifest["prompt_sha256"], hashlib.sha256(prompt.encode()).hexdigest())


class ReferenceTargetImportTests(unittest.TestCase):
    def _run(
        self,
        rows: list[dict[str, object]],
        root: Path,
        *,
        include_investigation: bool = True,
        focal_species: str = "Pisum sativum",
        benchmark_mode: bool = False,
    ) -> tuple[dict[str, object], dict[str, object], pd.DataFrame]:
        payload = root / "payload.json"
        payload.write_text(json.dumps({"hypotheses": rows}), encoding="utf-8")
        investigation_arg: str | None = None
        if include_investigation:
            investigation = phenotype_investigate(
                str(root / "investigation"),
                "case-654",
                focal_species,
                "late transition from vegetative to reproductive development",
            )
            investigation_arg = str(investigation["investigation"])
        result = mechanism_import(
            str(payload),
            str(root / "out"),
            investigation_arg=investigation_arg,
            case_id="case-654",
            benchmark_mode=benchmark_mode,
        )
        manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
        reviews = pd.read_csv(result["review_required"], sep="\t", dtype=str).fillna("")
        return result, manifest, reviews

    def test_complete_targeted_portfolio_remains_ready(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            _result, manifest, reviews = self._run(_portfolio(), Path(tmp))
            self.assertTrue(manifest["reference_target_contract_met"])
            self.assertTrue(manifest["family_node_portfolio_met"])
            self.assertEqual([], manifest["hypotheses_without_complete_reference_target"])
            self.assertNotIn("mechanism_without_reference_target", set(reviews["reason"]))
            self.assertEqual(
                "evaluated_from_investigation",
                manifest["reference_target_contract"]["focal_species_check_state"],
            )

    def test_p1_like_payload_without_any_target_is_explicit_review_missingness(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            rows = _portfolio()
            rows[0]["reference_gene_name"] = ""
            rows[0]["reference_species"] = ""
            rows[0]["proposed_stable_id"] = ""
            rows[0]["uncertainty"] = "no defensible non-focal reference target is known"
            result, manifest, reviews = self._run(rows, Path(tmp))
            self.assertFalse(manifest["reference_target_contract_met"])
            self.assertFalse(manifest["family_node_portfolio_met"])
            self.assertEqual(["mech_001"], manifest["hypotheses_without_complete_reference_target"])
            target_reviews = reviews[reviews["target"] == "mech_001"]
            self.assertIn("mechanism_without_reference_target", set(target_reviews["reason"]))
            self.assertFalse(target_reviews["reason"].str.contains("provider_unavailable").any())
            mechanisms = pd.read_csv(result["mechanisms"], sep="\t", dtype=str).fillna("")
            row = mechanisms[mechanisms["hypothesis_id"] == "mech_001"].iloc[0]
            self.assertEqual("", row["proposed_stable_id"])
            self.assertEqual("", row["reference_gene_name"])
            self.assertEqual("", row["reference_species"])

    def test_identity_without_species_has_its_own_review_reason(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            rows = _portfolio()
            rows[1]["reference_species"] = ""
            _result, manifest, reviews = self._run(rows, Path(tmp), benchmark_mode=True)
            self.assertEqual(1, manifest["reference_target_status_counts"]["missing_reference_species"])
            target_reviews = reviews[reviews["target"] == "mech_002"]
            self.assertIn("mechanism_reference_target_without_species", set(target_reviews["reason"]))

    def test_species_without_identity_has_its_own_review_reason(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            rows = _portfolio()
            rows[2]["reference_gene_name"] = ""
            rows[2]["proposed_stable_id"] = ""
            _result, manifest, reviews = self._run(rows, Path(tmp), benchmark_mode=True)
            self.assertEqual(1, manifest["reference_target_status_counts"]["missing_reference_gene_or_locus"])
            target_reviews = reviews[reviews["target"] == "mech_003"]
            self.assertIn("mechanism_reference_species_without_gene_or_locus", set(target_reviews["reason"]))

    def test_focal_species_target_never_satisfies_non_focal_contract(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            rows = _portfolio()
            rows[0]["reference_species"] = "Pisum sativum"
            _result, manifest, reviews = self._run(rows, Path(tmp), benchmark_mode=True)
            self.assertFalse(manifest["reference_target_contract_met"])
            self.assertEqual(
                1,
                manifest["reference_target_status_counts"]["focal_species_used_as_reference_target"],
            )
            target_reviews = reviews[reviews["target"] == "mech_001"]
            self.assertIn("mechanism_reference_target_is_focal_species", set(target_reviews["reason"]))

    def test_production_focal_reference_is_portfolio_ready(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            rows = _portfolio()
            rows[0]["reference_species"] = "Pisum sativum"
            _result, manifest, reviews = self._run(rows, Path(tmp))
            self.assertTrue(manifest["reference_target_contract_met"])
            self.assertEqual(1, manifest["reference_target_status_counts"]["declared_reference_target"])
            self.assertNotIn("mechanism_reference_target_is_focal_species", set(reviews["reason"]))

    def test_strict_prepare_policy_is_inherited_by_import(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows = _portfolio()
            rows[0]["reference_species"] = "Pisum sativum"
            payload = root / "payload.json"
            payload.write_text(json.dumps({"hypotheses": rows}), encoding="utf-8")
            investigation = phenotype_investigate(
                str(root / "investigation"), "case-654", "Pisum sativum",
                "late transition from vegetative to reproductive development",
                benchmark_mode=True,
            )
            result = mechanism_import(
                str(payload), str(root / "out"),
                investigation_arg=str(investigation["investigation"]), case_id="case-654",
            )
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertTrue(manifest["benchmark_mode"])
            self.assertEqual(1, manifest["reference_target_status_counts"]["focal_species_used_as_reference_target"])
            self.assertFalse(manifest["reference_target_contract_met"])

    def test_production_ambiguous_species_is_not_promoted_by_focal_match(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            rows = _portfolio()
            rows[0]["reference_species"] = "Pisum sativum / Arabidopsis thaliana"
            _result, manifest, reviews = self._run(rows, Path(tmp))
            self.assertFalse(manifest["reference_target_contract_met"])
            self.assertEqual(1, manifest["reference_target_status_counts"]["ambiguous_reference_species"])
            self.assertIn("mechanism_reference_species_is_ambiguous", set(reviews["reason"]))

    def test_malformed_strict_policy_is_not_silently_downgraded(self) -> None:
        from candidate_qtl.llm_pathway import hypothesis_verify
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            policy = root / "mechanism_import_manifest.json"
            policy.write_text('{"benchmark_mode": "false"}', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "benchmark_mode must be boolean"):
                hypothesis_verify(str(root / "hypotheses.tsv"), str(root / "out"), focal_species="Pisum sativum")

    def test_missing_investigation_keeps_non_focal_relation_not_evaluated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            _result, manifest, reviews = self._run(
                _portfolio(), Path(tmp), include_investigation=False
            )
            self.assertFalse(manifest["reference_target_contract_met"])
            self.assertFalse(manifest["family_node_portfolio_met"])
            self.assertEqual(
                8,
                manifest["reference_target_status_counts"][
                    "focal_species_not_available_for_reference_check"
                ],
            )
            self.assertEqual(
                "not_evaluated_investigation_not_supplied",
                manifest["reference_target_contract"]["focal_species_check_state"],
            )
            self.assertEqual(
                8,
                int(
                    reviews["reason"]
                    .eq("mechanism_non_focal_target_not_evaluated")
                    .sum()
                ),
            )

    def test_multiple_reference_species_are_review_required_not_silently_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            rows = _portfolio()
            rows[3]["reference_species"] = "Arabidopsis thaliana / Medicago truncatula"
            _result, manifest, reviews = self._run(rows, Path(tmp))
            self.assertFalse(manifest["reference_target_contract_met"])
            self.assertEqual(
                1,
                manifest["reference_target_status_counts"]["ambiguous_reference_species"],
            )
            target_reviews = reviews[reviews["target"] == "mech_004"]
            self.assertIn("mechanism_reference_species_is_ambiguous", set(target_reviews["reason"]))


if __name__ == "__main__":
    unittest.main()
