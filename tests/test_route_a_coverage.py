"""Truth-free tests for the Route A causal-layer coverage contract."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl.core import sha256
from candidate_qtl.phenotype_investigation import MECHANISM_COLUMNS, mechanism_import
from candidate_qtl.route_a_coverage import (
    mechanism_coverage_audit,
    phenotype_localization_freeze,
    phenotype_localization_packet,
)


def _investigation(root: Path, phenotype: str = "short roots with delayed gravitropic bending") -> Path:
    path = root / "phenotype_investigation.tsv"
    pd.DataFrame([{
        "case_id": "opaque_case", "focal_species": "Test plant", "phenotype": phenotype,
        "direction": "reduced", "inheritance": "recessive", "affected_tissue": "root",
        "developmental_stage": "seedling", "onset": "early", "penetrance": "complete",
        "synonyms": "", "ontology_terms": "", "known_exclusions": "",
        "researcher_constraints": "", "uncertainty": "assay timing varies",
    }]).to_csv(path, sep="\t", index=False, lineterminator="\n")
    return path


def _freeze(root: Path, investigation: Path, location: str, confidence: str) -> Path:
    payload = root / f"{location}_{confidence}.json"
    payload.write_text(json.dumps({
        "causal_localization": location,
        "localization_confidence": confidence,
        "supporting_phrase": "delayed gravitropic bending",
    }), encoding="utf-8")
    result = phenotype_localization_freeze(str(payload), str(investigation), str(root / f"freeze_{location}_{confidence}"))
    return Path(str(result["manifest"]))


def _mechanisms(root: Path, name: str, layers: list[str], family_prefix: str = "F", proposal_round: int = 1) -> Path:
    rows = []
    for index, layer in enumerate(layers, 1):
        row = {column: "" for column in MECHANISM_COLUMNS}
        row.update({
            "hypothesis_id": f"{name}_{index:02d}", "case_id": "opaque_case", "causal_layer": layer,
            "pathway_or_process": f"process {index}", "gene_family": f"{family_prefix}{index}",
            "mechanism_statement": f"mechanism {index}", "payload_source": "agent",
            "proposal_round": str(proposal_round), "interval_informed": "False", "review_status": "pending_verification",
            "rank_contribution": "none",
        })
        rows.append(row)
    path = root / f"{name}.tsv"
    pd.DataFrame(rows, columns=MECHANISM_COLUMNS).to_csv(path, sep="\t", index=False, lineterminator="\n")
    return path


def _import_manifest(root: Path, mechanisms: Path, localization: Path, proposal_round: int = 1) -> Path:
    frame = pd.read_csv(mechanisms, sep="\t", dtype=str).fillna("")
    distinct_families = frame["gene_family"].str.casefold().nunique()
    path = root / f"{mechanisms.stem}.import.json"
    path.write_text(json.dumps({
        "case_id": "opaque_case", "proposal_round": proposal_round,
        "interval_informed": False, "blind_discovery_claim": True,
        "family_node_portfolio_met": 8 <= len(frame) <= 15 and distinct_families == len(frame),
        "localization_freeze": {"path": str(localization), "sha256": sha256(localization)},
        "output_sha256": sha256(mechanisms), "outputs": {"mechanisms": str(mechanisms)},
    }), encoding="utf-8")
    return path


class RouteACoverageTests(unittest.TestCase):
    def test_degraded_portfolios_reach_the_guided_consumer_unchanged(self) -> None:
        from candidate_qtl.guided_profiles import _resolve_chained
        for layers in (["signaling"] * 12, ["regulation"]):
            with self.subTest(layers=layers), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                investigation = _investigation(root)
                freeze = _freeze(root, investigation, "broad", "low")
                initial = _mechanisms(root, "initial", layers)
                result = mechanism_coverage_audit(
                    [str(initial)], str(freeze), str(investigation), str(root / "coverage"),
                    mechanism_import_manifests_args=[str(_import_manifest(root, initial, freeze))])
                self.assertEqual("degraded_coverage", result["coverage_state"])
                manifest = json.loads(Path(result["manifest"]).read_text())
                self.assertTrue(manifest["freeze_ready"])
                self.assertEqual(result["mechanisms"], result["candidate_mechanisms"])
                self.assertEqual(result["hypotheses"], result["candidate_hypotheses"])
                pd.testing.assert_frame_equal(pd.read_csv(initial, sep="\t"),
                                              pd.read_csv(result["mechanisms"], sep="\t"))
                stage = {"arguments": {"mechanisms": "<from mechanism-coverage-audit.mechanisms>"}, "inputs_from": {"mechanisms": {
                    "stage": "mechanism-coverage-audit", "result_key": "mechanisms"}}}
                self.assertFalse(_resolve_chained(stage, [{"command": "mechanism-coverage-audit",
                                                          "status": "completed", "result": result}]))
                self.assertEqual(result["mechanisms"], stage["arguments"]["mechanisms"])

    def test_gap_prompt_is_self_contained_and_imports_with_round_specific_ids(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            investigation = _investigation(root)
            freeze = _freeze(root, investigation, "signaling", "high")
            initial = _mechanisms(root, "initial", ["signaling"] * 12)
            request = mechanism_coverage_audit(
                [str(initial)], str(freeze), str(investigation), str(root / "coverage"),
                mechanism_import_manifests_args=[str(_import_manifest(root, initial, freeze))])
            prompt = Path(request["gap_prompt"]).read_text(encoding="utf-8")
            template = json.loads(prompt.split("Return JSON only:\n", 1)[1])
            self.assertNotIn("same schema as the initial", prompt)
            audit = json.loads(Path(request["audit"]).read_text(encoding="utf-8"))
            rows = []
            for index, layer in enumerate(audit["requested_causal_layers"]):
                row = dict(template["hypotheses"][0])
                row.update(causal_layer=layer, confidence="low", gene_family=f"GAP{index}",
                           pathway_or_process=f"alternative process {index}",
                           mechanism_statement="A testable alternative for root development",
                           reference_species="Other plant", reference_gene_name=f"REF{index}",
                           counter_hypothesis="The defect could be downstream",
                           discriminating_observation="Measure root bending after stimulation")
                rows.append(row)
            payload = root / "gap.json"
            payload.write_text(json.dumps({"hypotheses": rows}), encoding="utf-8")
            result = mechanism_import(str(payload), str(root / "import"),
                                      investigation_arg=str(investigation), payload_source="agent",
                                      model="synthetic-test", proposal_round=2,
                                      localization_freeze_manifest=str(freeze))
            imported = pd.read_csv(result["mechanisms"], sep="\t").fillna("")
            self.assertEqual(len(rows), len(imported))
            self.assertTrue(imported.hypothesis_id.str.startswith("mech_r2_").all())
            self.assertEqual(list(audit["requested_causal_layers"]), imported.causal_layer.tolist())
            self.assertTrue(imported.source_ids.eq("").all())
            self.assertTrue(imported.proposed_stable_id.eq("").all())

    def test_initial_portfolio_accepts_flexible_counts_with_target_twelve_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            investigation = _investigation(root)
            freeze = _freeze(root, investigation, "signaling", "high")
            for count in (7, 8, 12, 13, 15, 16):
                with self.subTest(count=count):
                    layers = ["signaling"] * (count - 3) + ["regulation", "alternative", "alternative"]
                    initial = _mechanisms(root, f"initial_{count}", layers)
                    result = mechanism_coverage_audit(
                        [str(initial)], str(freeze), str(investigation), str(root / f"audit_{count}"),
                        mechanism_import_manifests_args=[str(_import_manifest(root, initial, freeze))])
                    self.assertEqual("pass", result["status"])
                    audit = json.loads(Path(str(result["audit"])).read_text(encoding="utf-8"))
                    self.assertEqual(12, audit["initial_hypothesis_target"])
                    self.assertEqual(count, audit["merged_row_count"])

    def test_initial_portfolio_rejects_empty_portfolio(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            investigation = _investigation(root)
            freeze = _freeze(root, investigation, "signaling", "high")
            initial = _mechanisms(root, "initial_empty", [])
            with self.assertRaisesRegex(ValueError, "at least one"):
                mechanism_coverage_audit(
                    [str(initial)], str(freeze), str(investigation), str(root / "audit_empty"),
                    mechanism_import_manifests_args=[str(_import_manifest(root, initial, freeze))])

    def test_localization_freeze_resume_preserves_bytes_and_rejects_changes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            investigation = _investigation(root)
            freeze = _freeze(root, investigation, "signaling", "high")
            before = freeze.read_bytes()
            payload = root / "signaling_high.json"
            with patch("candidate_qtl.route_a_coverage._now", return_value="different-time"):
                phenotype_localization_freeze(str(payload), str(investigation), str(freeze.parent))
            self.assertEqual(before, freeze.read_bytes())
            value = json.loads(payload.read_text(encoding="utf-8"))
            value["localization_confidence"] = "medium"
            payload.write_text(json.dumps(value), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "cannot be replaced"):
                phenotype_localization_freeze(str(payload), str(investigation), str(freeze.parent))
            self.assertEqual(before, freeze.read_bytes())
            value["localization_confidence"] = "high"
            payload.write_text(json.dumps(value), encoding="utf-8")
            changed = pd.read_csv(investigation, sep="\t", dtype=str).fillna("")
            changed.loc[0, "direction"] = "increased"
            changed.to_csv(investigation, sep="\t", index=False, lineterminator="\n")
            with self.assertRaisesRegex(ValueError, "cannot be replaced"):
                phenotype_localization_freeze(str(payload), str(investigation), str(freeze.parent))
            _investigation(root)
            frozen = freeze.parent / "phenotype_localization.frozen.json"
            frozen.write_text("{}", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "cannot be replaced"):
                phenotype_localization_freeze(str(payload), str(investigation), str(freeze.parent))
            self.assertEqual("{}", frozen.read_text(encoding="utf-8"))
            self.assertEqual(before, freeze.read_bytes())

    def test_gap_request_hash_is_reconstructed_and_same_directory_resume_preserves_history(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            investigation = _investigation(root)
            freeze = _freeze(root, investigation, "signaling", "high")
            initial = _mechanisms(root, "initial", ["signaling"] * 12)
            initial_manifest = _import_manifest(root, initial, freeze)
            out = root / "coverage"
            request = mechanism_coverage_audit(
                [str(initial)], str(freeze), str(investigation), str(out),
                mechanism_import_manifests_args=[str(initial_manifest)])
            audit_path = Path(str(request["audit"]))
            original_audit = audit_path.read_bytes()
            original_hash = sha256(audit_path)
            prompt_path = Path(str(request["gap_prompt"]))
            original_prompt = prompt_path.read_bytes()
            gap = _mechanisms(root, "gap", ["regulation", "alternative", "alternative"], "G", 2)
            gap_manifest = _import_manifest(root, gap, freeze, 2)
            args = ([str(initial), str(gap)], str(freeze), str(investigation), str(out))
            kwargs = {"gap_request_audit_arg": str(audit_path),
                      "mechanism_import_manifests_args": [str(initial_manifest), str(gap_manifest)]}
            altered = json.loads(original_audit)
            altered["gap_prompt_sha256"] = "a" * 64
            audit_path.write_text(json.dumps(altered), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "deterministic coverage request"):
                mechanism_coverage_audit(*args, **kwargs)
            audit_path.write_bytes(original_audit)
            prompt_path.write_text("unrelated content", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "existing gap-fill prompt"):
                mechanism_coverage_audit(*args, **kwargs)
            self.assertEqual(original_audit, audit_path.read_bytes())
            self.assertEqual("unrelated content", prompt_path.read_text(encoding="utf-8"))
            prompt_path.write_bytes(original_prompt)
            result = mechanism_coverage_audit(*args, **kwargs)
            self.assertEqual("pass", result["status"])
            audit = json.loads(audit_path.read_text(encoding="utf-8"))
            preserved = Path(audit["gap_request_audit"])
            self.assertNotEqual(audit_path, preserved)
            self.assertEqual(original_audit, preserved.read_bytes())
            self.assertEqual(original_hash, audit["gap_request_audit_sha256"])
            self.assertEqual(original_hash, sha256(preserved))
            self.assertEqual(original_prompt, Path(audit["gap_request_prompt"]).read_bytes())
            self.assertFalse(prompt_path.exists())
            self.assertIsNone(result["gap_prompt"])
            kwargs["gap_request_audit_arg"] = str(preserved)
            self.assertEqual("pass", mechanism_coverage_audit(*args, **kwargs)["status"])

    def test_initial_portfolio_counts_do_not_trust_import_self_report(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            investigation = _investigation(root)
            freeze = _freeze(root, investigation, "signaling", "high")
            for column in ("gene_family", "pathway_or_process"):
                with self.subTest(column=column):
                    initial = _mechanisms(root, column, ["signaling"] * 12)
                    frame = pd.read_csv(initial, sep="\t", dtype=str).fillna("")
                    frame[column] = "same"
                    frame.to_csv(initial, sep="\t", index=False, lineterminator="\n")
                    manifest = _import_manifest(root, initial, freeze)
                    value = json.loads(manifest.read_text(encoding="utf-8"))
                    value["family_node_portfolio_met"] = True
                    manifest.write_text(json.dumps(value), encoding="utf-8")
                    if column == "gene_family":
                        with self.assertRaisesRegex(ValueError, "distinct family"):
                            mechanism_coverage_audit(
                                [str(initial)], str(freeze), str(investigation), str(root / "audit"),
                                mechanism_import_manifests_args=[str(manifest)])
                    else:
                        result = mechanism_coverage_audit(
                            [str(initial)], str(freeze), str(investigation), str(root / "audit"),
                            mechanism_import_manifests_args=[str(manifest)])
                        audit = json.loads(Path(str(result["audit"])).read_text(encoding="utf-8"))
                        self.assertIn("distinct_mechanism_breadth_below_advisory_minimum", audit["initial_reasons"])

    def test_localization_packet_is_phenotype_only_and_rejects_candidate_leak(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            investigation = _investigation(root)
            result = phenotype_localization_packet(str(investigation), str(root / "packet"))
            prompt = Path(str(result["prompt"])).read_text(encoding="utf-8")
            self.assertIn("Do not name genes or gene families", prompt)
            candidates = root / "candidates.tsv"
            candidates.write_text("gene_id\nsecret_gene_77\n", encoding="utf-8")
            leaky = _investigation(root, "short roots involving secret_gene_77")
            with self.assertRaisesRegex(ValueError, "candidate identifier"):
                phenotype_localization_packet(str(leaky), str(root / "leaky"), str(candidates))

    def test_localization_freeze_requires_verbatim_phenotype_phrase(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            investigation = _investigation(root)
            bad = root / "bad.json"
            bad.write_text(json.dumps({"causal_localization": "signaling", "localization_confidence": "high",
                                       "supporting_phrase": "invented phrase"}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "verbatim"):
                phenotype_localization_freeze(str(bad), str(investigation), str(root / "freeze"))
            bad.write_text(json.dumps({"causal_localization": "signaling", "localization_confidence": "high",
                                       "supporting_phrase": "Delayed gravitropic bending"}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "verbatim"):
                phenotype_localization_freeze(str(bad), str(investigation), str(root / "freeze_case_changed"))

    def test_localization_manifest_must_remain_interval_blind_and_case_bound(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            investigation = _investigation(root)
            freeze = _freeze(root, investigation, "signaling", "high")
            manifest = json.loads(freeze.read_text(encoding="utf-8"))
            manifest["interval_informed"] = True
            freeze.write_text(json.dumps(manifest), encoding="utf-8")
            initial = _mechanisms(root, "initial", ["signaling"] * 9 + ["regulation", "alternative", "alternative"])
            with self.assertRaisesRegex(ValueError, "interval_informed=false"):
                mechanism_coverage_audit([str(initial)], str(freeze), str(investigation), str(root / "audit"),
                                         mechanism_import_manifests_args=[str(_import_manifest(root, initial, freeze))])

    def test_broad_low_confidence_requests_missing_layers_without_equalizing_counts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            investigation = _investigation(root)
            freeze = _freeze(root, investigation, "broad", "low")
            initial = _mechanisms(root, "initial", ["regulation"] * 8 + ["execution"] * 4)
            result = mechanism_coverage_audit([str(initial)], str(freeze), str(investigation), str(root / "audit"),
                                              mechanism_import_manifests_args=[str(_import_manifest(root, initial, freeze))])
            audit = json.loads(Path(str(result["audit"])).read_text(encoding="utf-8"))
            self.assertEqual("gap_fill_requested", audit["status"])
            self.assertEqual(["perception", "signaling", "execution"], audit["requested_causal_layers"])
            self.assertIn("single_layer_dominance_above_60_percent", audit["reasons"])

    def test_high_confidence_localization_allows_bias_with_adjacent_and_two_alternatives(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            investigation = _investigation(root)
            freeze = _freeze(root, investigation, "signaling", "high")
            initial = _mechanisms(root, "initial", ["signaling"] * 9 + ["regulation", "alternative", "alternative"])
            result = mechanism_coverage_audit([str(initial)], str(freeze), str(investigation), str(root / "audit"),
                                              mechanism_import_manifests_args=[str(_import_manifest(root, initial, freeze))])
            audit = json.loads(Path(str(result["audit"])).read_text(encoding="utf-8"))
            self.assertEqual("pass", audit["status"])
            self.assertEqual([], audit["requested_causal_layers"])

    def test_medium_confidence_requests_other_layers(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            investigation = _investigation(root)
            freeze = _freeze(root, investigation, "regulation", "medium")
            initial = _mechanisms(root, "initial", ["regulation"] * 11 + ["unclear"])
            result = mechanism_coverage_audit([str(initial)], str(freeze), str(investigation), str(root / "audit"),
                                              mechanism_import_manifests_args=[str(_import_manifest(root, initial, freeze))])
            audit = json.loads(Path(str(result["audit"])).read_text(encoding="utf-8"))
            self.assertEqual(["signaling", "alternative"], audit["requested_causal_layers"])

    def test_gap_fill_merge_deduplicates_family_and_stays_within_fifteen(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            investigation = _investigation(root)
            freeze = _freeze(root, investigation, "broad", "low")
            initial = _mechanisms(root, "initial", ["regulation"] * 8 + ["execution"] * 4)
            initial_manifest = _import_manifest(root, initial, freeze)
            request = mechanism_coverage_audit(
                [str(initial)], str(freeze), str(investigation), str(root / "request"),
                mechanism_import_manifests_args=[str(initial_manifest)])
            gap = _mechanisms(root, "gap", ["perception", "signaling", "execution"], family_prefix="G", proposal_round=2)
            gap_frame = pd.read_csv(gap, sep="\t", dtype=str).fillna("")
            gap_frame.loc[0, "gene_family"] = "F1"
            gap_frame.to_csv(gap, sep="\t", index=False, lineterminator="\n")
            gap_manifest = _import_manifest(root, gap, freeze, 2)
            result = mechanism_coverage_audit(
                [str(initial), str(gap)], str(freeze), str(investigation), str(root / "final"),
                gap_request_audit_arg=str(request["audit"]),
                mechanism_import_manifests_args=[str(initial_manifest), str(gap_manifest)])
            audit = json.loads(Path(str(result["audit"])).read_text(encoding="utf-8"))
            self.assertEqual(14, audit["merged_row_count"])
            self.assertEqual(["gap_01"], audit["dropped_duplicate_hypothesis_ids"])
            self.assertEqual("review_required", audit["status"])
            self.assertEqual("degraded_coverage", result["coverage_state"])
            self.assertEqual(14, len(pd.read_csv(result["mechanisms"], sep="\t")))
            self.assertTrue(json.loads(Path(result["manifest"]).read_text())["freeze_ready"])

    def test_interval_informed_rows_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            investigation = _investigation(root)
            freeze = _freeze(root, investigation, "signaling", "high")
            initial = _mechanisms(root, "initial", ["signaling"] * 9 + ["regulation", "alternative", "alternative"])
            frame = pd.read_csv(initial, sep="\t", dtype=str).fillna("")
            frame["interval_informed"] = "True"
            frame.to_csv(initial, sep="\t", index=False, lineterminator="\n")
            with self.assertRaisesRegex(ValueError, "interval-informed input is forbidden"):
                mechanism_coverage_audit([str(initial)], str(freeze), str(investigation), str(root / "audit"),
                                         mechanism_import_manifests_args=[str(_import_manifest(root, initial, freeze))])

    def test_candidate_leak_and_duplicate_hypothesis_ids_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            investigation = _investigation(root)
            freeze = _freeze(root, investigation, "signaling", "high")
            initial = _mechanisms(root, "initial", ["signaling"] * 9 + ["regulation", "alternative", "alternative"])
            frame = pd.read_csv(initial, sep="\t", dtype=str).fillna("")
            frame.loc[1, "hypothesis_id"] = frame.loc[0, "hypothesis_id"]
            frame.to_csv(initial, sep="\t", index=False, lineterminator="\n")
            with self.assertRaisesRegex(ValueError, "unique non-empty hypothesis_id"):
                mechanism_coverage_audit(
                    [str(initial)], str(freeze), str(investigation), str(root / "duplicate"),
                    mechanism_import_manifests_args=[str(_import_manifest(root, initial, freeze))])

            frame.loc[1, "hypothesis_id"] = "initial_02"
            frame.loc[1, "mechanism_statement"] = "secret_gene_77 controls the phenotype"
            frame.to_csv(initial, sep="\t", index=False, lineterminator="\n")
            candidates = root / "candidates.tsv"
            candidates.write_text("gene_id\nsecret_gene_77\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "candidate identifier"):
                mechanism_coverage_audit(
                    [str(initial)], str(freeze), str(investigation), str(root / "leak"), str(candidates),
                    mechanism_import_manifests_args=[str(_import_manifest(root, initial, freeze))])

    def test_initial_duplicate_families_and_swapped_localization_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            investigation = _investigation(root)
            freeze = _freeze(root, investigation, "signaling", "high")
            initial = _mechanisms(root, "initial", ["signaling"] * 9 + ["regulation", "alternative", "alternative"])
            frame = pd.read_csv(initial, sep="\t", dtype=str).fillna("")
            frame.loc[:8, "gene_family"] = "DUPLICATE"
            frame.to_csv(initial, sep="\t", index=False, lineterminator="\n")
            with self.assertRaisesRegex(ValueError, "family-node portfolio"):
                mechanism_coverage_audit(
                    [str(initial)], str(freeze), str(investigation), str(root / "duplicates"),
                    mechanism_import_manifests_args=[str(_import_manifest(root, initial, freeze))])

            valid = _mechanisms(root, "valid", ["signaling"] * 9 + ["regulation", "alternative", "alternative"])
            valid_manifest = _import_manifest(root, valid, freeze)
            swapped = _freeze(root, investigation, "regulation", "high")
            with self.assertRaisesRegex(ValueError, "not bound to the supplied localization freeze"):
                mechanism_coverage_audit(
                    [str(valid)], str(swapped), str(investigation), str(root / "swapped"),
                    mechanism_import_manifests_args=[str(valid_manifest)])
            changed = pd.read_csv(investigation, sep="\t", dtype=str).fillna("")
            changed.loc[0, "direction"] = "increased"
            changed.to_csv(investigation, sep="\t", index=False, lineterminator="\n")
            with self.assertRaisesRegex(ValueError, "supplied phenotype packet"):
                mechanism_coverage_audit(
                    [str(valid)], str(freeze), str(investigation), str(root / "changed_context"),
                    mechanism_import_manifests_args=[str(valid_manifest)])


if __name__ == "__main__":
    unittest.main()
