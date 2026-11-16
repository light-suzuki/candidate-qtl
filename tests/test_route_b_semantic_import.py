from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl import route_b_semantic_import as semantic_import
from candidate_qtl.route_b_semantic_bundle import build_route_b_semantic_bundle
from candidate_qtl.route_b_semantic_import import import_route_b_semantic_proposals
from candidate_qtl.semantic_claims import CLAIM_COLUMNS, narrowing_eligible_claims, normalize_semantic_claims
from scripts.agentic_route_b import merge


def _single_row_errors(rows: list[dict[str, str]]) -> dict[int, str]:
    """The importer's previous one-row validation, used as a parity oracle."""
    errors = {}
    for index, row in enumerate(rows):
        try:
            normalize_semantic_claims(pd.DataFrame([row], columns=CLAIM_COLUMNS))
        except ValueError as error:
            errors[index] = str(error)
    return errors


class RouteBSemanticImportTests(unittest.TestCase):
    def test_invalid_claim_type_is_quarantined_and_mixed_claims_survive(self) -> None:
        good = {"claim_type": "mutant_phenotype", "reference_gene": "REF1",
                "evidence_scope": "reference_gene", "source_ids": ["PMID:123"]}
        for index, bad in enumerate(({}, {"claim_type": ""}, {"claim_type": "  "},
                                     {"claim_type": None}, {"claim_type": 12}, "bad claim")):
            with self.subTest(bad=bad), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                built = self._bundle(root)
                frames = {}
                for name, claims in (("valid", [good]), ("mixed", [good, bad]), ("invalid", [bad])):
                    response = root / f"{name}.json"
                    self._response(Path(built["prompts"]), claims, response)
                    result = import_route_b_semantic_proposals(built["bundle"], built["prompts"], str(response), str(root / name))
                    frames[name] = pd.read_csv(result["proposals"], sep="\t").drop(columns=["output_fingerprint"])
                    status = pd.read_csv(result["status"], sep="\t").iloc[0]
                    self.assertEqual(len(claims), status["model_claim_count"])
                    self.assertEqual(int(name != "valid"), result["rejected_claim_row_count"])
                    self.assertEqual({"valid": "proposals_imported", "mixed": "partial_import",
                                      "invalid": "not_evaluated_all_claim_rows_rejected"}[name], status["status"])
                    if name != "valid":
                        rejected = json.loads(Path(result["rejected_rows"]).read_text())
                        self.assertEqual(bad, rejected[0]["raw_claim"])
                        self.assertTrue(rejected[0]["reason"])
                    self.assertEqual(0, result["verified_claim_row_count"])
                pd.testing.assert_frame_equal(frames["valid"], frames["mixed"])
                self.assertTrue(frames["invalid"].empty)

    def test_merge_import_checkpoint_retries_invalid_claim_type_and_keeps_valid_proposal(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            built = self._bundle(root)
            prompt = json.loads(Path(built["prompts"]).read_text().splitlines()[0])
            invalid = {"mechanism": "receptor process", "reference_gene": "REF1", "evidence_scope": "reference_gene"}
            good = {**invalid, "claim_type": "pathway_relation"}
            for name, claims in (("invalid", [invalid]), ("mixed", [good, invalid]), ("corrected", [good])):
                path = root / f"{name}.json"
                raw = json.dumps({"claims": claims})
                path.write_text(json.dumps({"case_id": "case403", "executor": "fixture", "model": "fixture",
                                            "responses": [{"prompt_id": prompt["prompt_id"], "raw_text": raw}]}))
            for name in ("invalid", "mixed"):
                checkpoint = root / f"{name}-checkpoint.json"
                initial = merge(built["prompts"], [str(root / f"{name}.json")], str(checkpoint))
                self.assertEqual(("partial", 0, 1), (initial["completion_status"], initial["completed_count"], initial["not_evaluated_count"]))
                imported = import_route_b_semantic_proposals(built["bundle"], built["prompts"], str(checkpoint), str(root / f"import-{name}"))
                self.assertEqual(1, imported["rejected_claim_row_count"])
                self.assertEqual(int(name == "mixed"), imported["claim_row_count"])
                rejected = json.loads(Path(imported["rejected_rows"]).read_text())
                self.assertEqual(invalid, rejected[0]["raw_claim"])
                original_bytes = checkpoint.read_bytes()
                resumed_path = root / f"{name}-resumed.json"
                resumed = merge(built["prompts"], [str(checkpoint), str(root / "corrected.json")], str(resumed_path))
                self.assertEqual(("complete", 1, 0), (resumed["completion_status"], resumed["completed_count"], resumed["not_evaluated_count"]))
                self.assertEqual(1, len(resumed["retry_history"]))
                self.assertEqual("model_response_invalid_claim_schema", resumed["retry_history"][0]["replaced_reason"])
                self.assertIn("claim_type", resumed["retry_history"][0]["replaced_response_format_error"])
                self.assertEqual(original_bytes, checkpoint.read_bytes())
                final = import_route_b_semantic_proposals(built["bundle"], built["prompts"], str(resumed_path), str(root / f"final-{name}"))
                self.assertEqual(1, final["proposed_claim_row_count"])
                self.assertEqual(0, final["rejected_claim_row_count"])
                self.assertEqual(0, final["verified_claim_row_count"])

    def test_bad_scope_anchor_is_quarantined_without_losing_valid_claim(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            built = self._bundle(root)
            bundle = json.loads(Path(built["bundle"]).read_text())
            for unit in bundle["units"]:
                unit["identity"].pop("evolutionary_family_ids", None)
                unit["identity"].pop("family_ids", None)
            Path(built["bundle"]).write_text(json.dumps(bundle))
            good = {"claim_type": "mutant_phenotype", "evidence_scope": "reference_gene",
                    "source_ids": ["PMID:123"], "reference_gene": "REF1"}
            bad = {"claim_type": "mutant_phenotype", "evidence_scope": "family",
                   "source_ids": ["PMID:123"], "reference_gene": "REF1"}
            frames = {}
            for name, claims in (("alone", [good]), ("mixed", [good, bad]), ("invalid", [bad])):
                response = root / f"{name}.json"
                self._response(Path(built["prompts"]), claims, response)
                result = import_route_b_semantic_proposals(built["bundle"], built["prompts"], str(response), str(root / name))
                frames[name] = pd.read_csv(result["proposals"], sep="\t").drop(columns=["output_fingerprint"])
                self.assertEqual(int(name != "alone"), result["rejected_claim_row_count"])
                if name == "invalid":
                    status = pd.read_csv(result["status"], sep="\t")
                    self.assertEqual("not_evaluated_all_claim_rows_rejected", status.iloc[0]["status"])
                    rejected = json.loads(Path(result["rejected_rows"]).read_text())
                    self.assertEqual(bad, rejected[0]["raw_claim"])
            pd.testing.assert_frame_equal(frames["alone"], frames["mixed"])
            self.assertTrue(frames["invalid"].empty)

    def _bundle(self, root: Path) -> dict[str, object]:
        atlas = root / "atlas.tsv"
        pd.DataFrame([
            {
                "gene_id": "G1", "family_id": "FAM1", "evolutionary_family_id": "EV1",
                "orthogroup": "OG1", "subfamily_id": "SUB1", "orthology_reference_gene": "REF1",
                "domain_architecture": "kinase", "broad_function_bins": "signaling",
            },
            {
                "gene_id": "G2", "family_id": "FAM1", "evolutionary_family_id": "EV1",
                "orthogroup": "OG1", "subfamily_id": "SUB1", "orthology_reference_gene": "REF1",
                "domain_architecture": "kinase", "broad_function_bins": "signaling",
            },
        ]).to_csv(atlas, sep="\t", index=False, lineterminator="\n")
        literature = root / "literature.tsv"
        pd.DataFrame([
            {
                "interval_gene_ids": "G1;G2", "reference_gene": "REF1", "reference_function": "meristem receptor",
                "query_status": "retrieved", "article_rank": "1", "source_id": "PMID:123",
                "title": "Reference receptor controls meristem development", "doi": "10.1/example",
                "journal": "Plant", "year": "2020", "matched_phenotype_terms": "fasciation",
                "claim_type": "literature_link_candidate",
            }
        ]).to_csv(literature, sep="\t", index=False, lineterminator="\n")
        return build_route_b_semantic_bundle(
            str(atlas), str(root / "bundle"), phenotype="fasciated meristem",
            literature_arg=str(literature), case_id="case403",
        )

    def _response(self, prompts: Path, claims: list[dict[str, object]], out: Path) -> None:
        prompt = json.loads(prompts.read_text(encoding="utf-8").splitlines()[0])
        out.write_text(json.dumps({
            "candidates": [{
                "prompt_id": prompt["prompt_id"],
                "gene_id": prompt["gene_id"],
                "claims": claims,
                "model_json": {"claims": claims},
                "raw_text": json.dumps({"claims": claims}),
                "parse_status": "json",
                "review_status": "not_reviewed",
            }]
        }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    def test_known_source_is_only_proposed_until_relation_verification(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            built = self._bundle(root)
            prompts = Path(built["prompts"])
            responses = root / "responses.json"
            self._response(prompts, [{
                "claim_type": "mutant_phenotype",
                "biological_function": "meristem maintenance",
                "mechanism": "receptor signaling",
                "phenotype_relation": "could explain fasciation",
                "reference_gene": "REF1",
                "reference_species": "Arabidopsis thaliana",
                "evidence_scope": "reference_gene",
                "source_ids": ["PMID:123"],
                "alternative_hypothesis": "related receptor copy",
                "uncertainty_note": "transfer to focal crop remains unresolved",
            }], responses)
            result = import_route_b_semantic_proposals(
                built["bundle"], built["prompts"], str(responses), str(root / "import")
            )
            frame = pd.read_csv(result["proposals"], sep="\t", dtype=str).fillna("")
            self.assertEqual(["proposed"], frame["verification_status"].tolist())
            self.assertIn("verification_pending", frame.iloc[0]["verification_reason"])
            self.assertEqual(0, len(narrowing_eligible_claims(frame)))
            self.assertEqual(0, result["verified_claim_row_count"])

    def test_source_free_mechanistic_review_seed_is_proposed_but_never_supported(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            built = self._bundle(root)
            responses = root / "responses.json"
            self._response(Path(built["prompts"]), [{
                "claim_type": "mechanistic_review_seed",
                "biological_function": "developmental receptor signaling",
                "mechanism": "altered meristem boundary signaling could change organ initiation",
                "phenotype_relation": "specific hypothesis for fasciation; trait support is unverified",
                "reference_gene": "REF1",
                "reference_species": "Arabidopsis thaliana",
                "evidence_scope": "family",
                "source_ids": [],
                "alternative_hypothesis": "an unrelated meristem pathway may be causal",
                "uncertainty_note": "review seed only; no direct phenotype evidence was supplied",
            }], responses)
            result = import_route_b_semantic_proposals(
                built["bundle"], built["prompts"], str(responses), str(root / "import")
            )
            frame = pd.read_csv(result["proposals"], sep="\t", dtype=str).fillna("")
            self.assertEqual(["mechanistic_review_seed"], frame["claim_type"].tolist())
            self.assertEqual(["proposed"], frame["verification_status"].tolist())
            self.assertEqual([""], frame["source_id"].tolist())
            self.assertIn("no_source_id_supplied", frame.iloc[0]["verification_reason"])
            self.assertEqual(0, len(narrowing_eligible_claims(frame)))
            self.assertEqual(0, result["verified_claim_row_count"])

    def test_invented_source_or_reference_gene_is_unsupported(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            built = self._bundle(root)
            prompts = Path(built["prompts"])
            responses = root / "responses.json"
            self._response(prompts, [{
                "claim_type": "pathway_relation",
                "biological_function": "invented function",
                "mechanism": "invented mechanism",
                "phenotype_relation": "invented relation",
                "reference_gene": "INVENTED_REF",
                "reference_species": "Arabidopsis thaliana",
                "evidence_scope": "reference_gene",
                "source_ids": ["PMID:999999999"],
                "alternative_hypothesis": "",
                "uncertainty_note": "",
            }], responses)
            result = import_route_b_semantic_proposals(
                built["bundle"], built["prompts"], str(responses), str(root / "import")
            )
            frame = pd.read_csv(result["proposals"], sep="\t", dtype=str).fillna("")
            self.assertEqual(["unsupported"], frame["verification_status"].tolist())
            self.assertIn("source_id_not_in_frozen_bundle", frame.iloc[0]["verification_reason"])
            self.assertIn("reference_gene_not_in_frozen_bundle", frame.iloc[0]["verification_reason"])
            self.assertEqual(1, result["unsupported_claim_row_count"])

    def test_response_roster_must_exactly_match_frozen_prompts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            built = self._bundle(root)
            responses = root / "responses.json"
            responses.write_text(json.dumps({"candidates": []}) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "roster mismatch"):
                import_route_b_semantic_proposals(
                    built["bundle"], built["prompts"], str(responses), str(root / "import")
                )

    def test_empty_claims_is_valid_abstention(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            built = self._bundle(root)
            responses = root / "responses.json"
            self._response(Path(built["prompts"]), [], responses)
            result = import_route_b_semantic_proposals(
                built["bundle"], built["prompts"], str(responses), str(root / "import")
            )
            status = pd.read_csv(result["status"], sep="\t", dtype=str).fillna("")
            self.assertEqual("no_claims", status.iloc[0]["status"])
            self.assertEqual(0, result["claim_row_count"])

    def test_manifest_hash_mismatch_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            built = self._bundle(root)
            responses = root / "responses.json"
            self._response(Path(built["prompts"]), [], responses)
            manifest = root / "llm_run_manifest.json"

            def digest(path: str | Path) -> str:
                return hashlib.sha256(Path(path).read_bytes()).hexdigest()

            manifest.write_text(json.dumps({
                "provider": "test-provider",
                "model": "test-model",
                "bundle_sha256": digest(built["bundle"]),
                "prompt_sha256": digest(built["prompts"]),
                "response_sha256": digest(responses),
            }) + "\n", encoding="utf-8")
            responses.write_bytes(responses.read_bytes() + b"\n")

            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                import_route_b_semantic_proposals(
                    built["bundle"], built["prompts"], str(responses), str(root / "import"),
                    llm_manifest_arg=str(manifest),
                )

    def test_provider_model_only_manifest_is_not_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            built = self._bundle(root)
            responses = root / "responses.json"
            self._response(Path(built["prompts"]), [], responses)
            manifest = root / "llm_run_manifest.json"
            manifest.write_text(json.dumps({"provider": "legacy", "model": "legacy-model"}) + "\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "missing frozen artifact hashes"):
                import_route_b_semantic_proposals(
                    built["bundle"], built["prompts"], str(responses), str(root / "import"),
                    llm_manifest_arg=str(manifest),
                )

    def test_batched_validation_preserves_single_row_error_text(self) -> None:
        valid = dict.fromkeys(CLAIM_COLUMNS, "")
        valid.update({
            "claim_id": "valid", "route": "interval_first", "semantic_unit_id": "unit",
            "family_id": " FAM1 ", "orthogroup_id": "OG1", "member_interval_genes": " G2, G1;G2 ",
            "reference_gene": "REF1", "claim_type": "mutant_phenotype", "evidence_scope": "unknown",
            "verification_status": "proposed", "rank_contribution": "none", "source_id": "PMID:123",
            "mechanism": "mechanism", "phenotype_relation": "relation", "uncertainty_note": "uncertain",
        })
        invalid_overrides = [
            {"claim_id": " "}, {"route": " "}, {"semantic_unit_id": " "}, {"claim_type": " "},
            {"family_id": "", "orthogroup_id": "", "member_interval_genes": None, "reference_gene": ""},
            {"verification_status": "invalid"}, {"evidence_scope": "invalid"},
            {"claim_type": "mechanistic_review_seed", "verification_status": "verified"},
            *({"claim_type": "mechanistic_review_seed", column: " "}
              for column in ("mechanism", "phenotype_relation", "uncertainty_note")),
            {"evidence_scope": "family", "family_id": ""},
            {"evidence_scope": "orthogroup", "orthogroup_id": ""},
            {"evidence_scope": "interval_gene", "member_interval_genes": " ,; "},
            {"evidence_scope": "reference_gene", "reference_gene": ""},
            {"rank_contribution": "1"},
            {"verification_status": "verified", "source_id": ""},
            {"verification_status": "verified", "verification_reason": None},
            {"verification_status": "provider_failed", "verification_reason": " "},
            {"claim_type": " ", "evidence_scope": "invalid", "rank_contribution": "1"},
        ]
        invalid = [{**valid, "claim_id": f"bad-{index}", **changes}
                   for index, changes in enumerate(invalid_overrides)]
        mixed = [item for index, row in enumerate(invalid)
                 for item in ({**valid, "claim_id": f"good-{index}"}, row)]
        for rows in ([], [valid], [invalid[0]], invalid, mixed, [valid, valid], [valid, {}, valid]):
            with self.subTest(row_count=len(rows)):
                self.assertEqual(_single_row_errors(rows), semantic_import._semantic_claim_row_errors(rows))

    def test_batch_import_artifacts_match_single_row_validation_and_keep_quarantine_order(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            built = self._bundle(root)
            bundle = json.loads(Path(built["bundle"]).read_text())
            for unit in bundle["units"]:
                unit["identity"].pop("evolutionary_family_ids", None)
                unit["identity"].pop("family_ids", None)
            Path(built["bundle"]).write_text(json.dumps(bundle))
            good = {"claim_type": " mutant_phenotype ", "reference_gene": " REF1;OTHER ",
                    "evidence_scope": "reference_gene", "mechanism": None, "biological_function": " 生育 ",
                    "source_ids": ["PMID:987", " PMID:123 ", "PMID:123", " "]}
            bad_scope = {**good, "evidence_scope": "family"}
            bad_review = {"claim_type": "mechanistic_review_seed", "source_ids": [], "mechanism": " "}
            batches = [
                ("valid", [good, {**good, "source_ids": [], "evidence_scope": " "}]),
                ("mixed", [bad_scope, {}, good, bad_review, "invalid claim"]),
                ("invalid", [bad_scope, bad_review]),
                ("empty", []),
            ]
            for name, claims in batches:
                with self.subTest(name=name):
                    responses = root / f"{name}.json"
                    self._response(Path(built["prompts"]), claims, responses)
                    out = root / name
                    with patch.object(semantic_import, "_semantic_claim_row_errors", side_effect=_single_row_errors):
                        expected = import_route_b_semantic_proposals(
                            built["bundle"], built["prompts"], str(responses), str(out))
                    expected_bytes = {path.name: path.read_bytes() for path in out.iterdir()}
                    actual = import_route_b_semantic_proposals(
                        built["bundle"], built["prompts"], str(responses), str(out))
                    self.assertEqual(expected, actual)
                    self.assertEqual(expected_bytes, {path.name: path.read_bytes() for path in out.iterdir()})
                    if name == "mixed":
                        rejected = json.loads(Path(actual["rejected_rows"]).read_text())
                        self.assertEqual([(1, 1), (1, 2), (2, 0), (4, 1), (5, 0)],
                                         [(row["claim_index"], row["source_index"]) for row in rejected])

            prompt = json.loads(Path(built["prompts"]).read_text().splitlines()[0])
            pre_rejected = [{"claim_index": index, "reason": "invalid shape", "raw_claim": {}}
                            for index in (9, 2)]
            responses = root / "pre-rejected.json"
            responses.write_text(json.dumps({"candidates": [{"prompt_id": prompt["prompt_id"],
                                  "parse_status": "not_evaluated", "rejected_claims": pre_rejected}]}))
            result = import_route_b_semantic_proposals(
                built["bundle"], built["prompts"], str(responses), str(root / "pre-rejected"))
            rejected = json.loads(Path(result["rejected_rows"]).read_text())
            self.assertEqual([9, 2], [row["claim_index"] for row in rejected])

    def test_unit_provenance_and_validation_work_are_reused_without_changing_frozen_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            built = self._bundle(root)
            bundle = json.loads(Path(built["bundle"]).read_text())
            unit = bundle["units"][0]
            unit["stable_source_ids"] = ["PMID:123", "PMID:456"]
            other = deepcopy(unit)
            other["unit_id"] += "-second"
            bundle["units"].append(other)
            Path(built["bundle"]).write_text(json.dumps(bundle))
            prompt = json.loads(Path(built["prompts"]).read_text().splitlines()[0])
            prompts = [prompt, {**prompt, "prompt_id": prompt["prompt_id"] + "-second",
                                "route_b_unit_id": other["unit_id"]}]
            Path(built["prompts"]).write_text("".join(json.dumps(record) + "\n" for record in prompts))
            records = [{"prompt_id": record["prompt_id"], "parse_status": "json", "model_json": {"claims": [
                {"claim_type": "mutant_phenotype", "reference_gene": "REF1", "evidence_scope": "reference_gene",
                 "source_ids": ["PMID:123", "PMID:456"], "mechanism": f"{record['prompt_id']} mechanism {index}"}
                for index in range(3)]}} for record in prompts]
            responses = root / "responses.json"
            responses.write_text(json.dumps({"candidates": records}))
            frozen = {Path(path): Path(path).read_bytes() for path in (built["bundle"], built["prompts"], responses)}
            with patch.object(semantic_import, "_sha256_json", wraps=semantic_import._sha256_json) as digest, \
                    patch.object(semantic_import, "normalize_semantic_claims", wraps=normalize_semantic_claims) as normalize:
                result = import_route_b_semantic_proposals(
                    built["bundle"], built["prompts"], str(responses), str(root / "import"))
            self.assertEqual(12, result["claim_row_count"])
            self.assertEqual(6, digest.call_count)  # Three frozen JSON objects for each of two units.
            self.assertEqual(3, normalize.call_count)  # Two units plus the unchanged final table gate.
            self.assertEqual(frozen, {path: path.read_bytes() for path in frozen})
            frame = pd.read_csv(result["proposals"], sep="\t", dtype=str).fillna("")
            for unit, prompt, response in zip(bundle["units"], prompts, records):
                rows = frame[frame["semantic_unit_id"].eq(unit["unit_id"])]
                for column, value in (("prompt_fingerprint", prompt), ("input_fingerprint", unit),
                                      ("output_fingerprint", response["model_json"])):
                    encoded = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
                    self.assertEqual({hashlib.sha256(encoded.encode("utf-8")).hexdigest()}, set(rows[column]))

    def test_repeated_unit_responses_still_fail_final_claim_id_uniqueness_gate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            built = self._bundle(root)
            prompt = json.loads(Path(built["prompts"]).read_text().splitlines()[0])
            prompts = [prompt, {**prompt, "prompt_id": prompt["prompt_id"] + "-second"}]
            Path(built["prompts"]).write_text("".join(json.dumps(record) + "\n" for record in prompts))
            responses = root / "responses.json"
            responses.write_text(json.dumps({"candidates": [
                {"prompt_id": record["prompt_id"], "parse_status": "json", "model_json": {"claims": [
                    {"claim_type": "mutant_phenotype", "reference_gene": "REF1", "evidence_scope": "reference_gene"}]}}
                for record in prompts]}))
            with self.assertRaisesRegex(ValueError, "semantic claim_id values must be unique"):
                import_route_b_semantic_proposals(
                    built["bundle"], built["prompts"], str(responses), str(root / "import"))


if __name__ == "__main__":
    unittest.main()
