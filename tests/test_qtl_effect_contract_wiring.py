"""Issue #360 E2E: case-run wiring of the optional QTL-effect contract.

The contract is an audit-only input: it annotates a copy of the locus-hypothesis
ledger (hypotheses/locus_hypothesis_sign_audit.tsv), records a
``qtl_effect_contract`` block in run_manifest.json / case_state.json, and never
touches rank/score/weight outputs (rank_contribution=none). ``incompatible`` is
review-visible only: no exclusion, no penalty, no reviewed marking.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import tempfile
import unittest
from pathlib import Path

import pandas as pd
import yaml

from candidate_qtl.cli import main
from candidate_qtl.hypothesis_ledger import HYPOTHESIS_COLUMNS

POSITIVE_CONTRACT = (
    "case_id: bundled_fas_like\n"
    "qtl_effect_contract:\n"
    "  favorable_parent: ParentA\n"
    "  effect_direction: increase\n"
    "  additive_effect: '3.2'\n"
    "  dominance_class: additive\n"
    "  environment_id: field_2026\n"
    "  stable_across_environments: true\n"
    "  source_id: PMID:360\n"
    "  source_detail: Table 2\n"
)

AUDIT_COLUMNS = [
    "expected_mechanism_sign",
    "sign_concordance",
    "sign_concordance_reason",
    "qtl_effect_source_id",
    "qtl_effect_source_detail",
    "qtl_effect_environment_id",
]


def _hypotheses(work: Path, sign: str = "increase") -> Path:
    row = {column: "" for column in HYPOTHESIS_COLUMNS}
    row.update({
        "hypothesis_id": "h_fas_like", "case_id": "bundled_fas_like", "hypothesis_type": "single_gene_coding",
        "member_gene_ids": "GeneStrong", "proposed_mechanism": "FAS-like gain-of-function",
        "status": "proposed", "created_by": "manual:researcher", "proposal_origin": "manual",
        "review_status": "reviewed", "source_id": "local:proposal:h_fas_like",
        "gene_disposition": "retained", "rejection_scope": "none",
        "expected_mechanism_sign": sign,
    })
    path = work / "hypotheses.tsv"
    pd.DataFrame([row]).to_csv(path, sep="\t", index=False, lineterminator="\n")
    return path


class QtlEffectContractWiringTests(unittest.TestCase):
    def _work(self, root: Path) -> Path:
        fixture = Path(__file__).parent / "fixtures" / "case_run"
        work = root / "case"
        shutil.copytree(fixture, work)
        return work

    def _run_case(self, work: Path, out: Path, *, inputs: dict | None = None,
                  resources: dict | None = None, resume: bool = False,
                  expect_fail: bool = False) -> Path:
        case = yaml.safe_load((work / "case.yaml").read_text(encoding="utf-8"))
        if inputs:
            case["inputs"].update(inputs)
        case_path = work / "case_run.yaml"
        case_path.write_text(yaml.safe_dump(case, sort_keys=False), encoding="utf-8")
        registry = yaml.safe_load((work / "resources.yaml").read_text(encoding="utf-8"))
        if resources:
            registry["resources"].update(resources)
        (work / "resources.yaml").write_text(yaml.safe_dump(registry, sort_keys=False), encoding="utf-8")
        argv = ["case-run", "--case", str(case_path), "--resources", str(work / "resources.yaml"),
                "--out-dir", str(out)]
        if resume:
            argv.append("--resume")
        if expect_fail:
            with self.assertRaises(ValueError):
                main(argv)
        else:
            self.assertEqual(0, main(argv))
        return out

    def _contract_run(self, work: Path, out: Path, *, contract: str = POSITIVE_CONTRACT,
                      sign: str = "increase", resume: bool = False,
                      expect_fail: bool = False) -> tuple[Path, Path]:
        contract_path = work / "qtl_effect_contract.yaml"
        contract_path.write_text(contract, encoding="utf-8")
        _hypotheses(work, sign=sign)
        inputs = {"hypotheses_resource": "hypotheses", "qtl_effect_contract_resource": "qtl_contract"}
        resources = {
            "hypotheses": {"type": "hypotheses", "path": "hypotheses.tsv"},
            "qtl_contract": {"type": "qtl_effect_contract", "path": "qtl_effect_contract.yaml"},
        }
        self._run_case(work, out, inputs=inputs, resources=resources, resume=resume, expect_fail=expect_fail)
        return contract_path, work / "hypotheses.tsv"

    def _audit(self, out: Path) -> pd.DataFrame:
        return pd.read_csv(out / "hypotheses" / "locus_hypothesis_sign_audit.tsv", sep="\t", dtype=str).fillna("")

    def _qtl_block(self, path: Path) -> dict:
        return json.loads(path.read_text(encoding="utf-8"))["qtl_effect_contract"]

    def test_positive_contract_annotates_concordant_and_never_changes_rank(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            baseline_work = self._work(root / "baseline")
            _hypotheses(baseline_work)
            baseline = self._run_case(baseline_work, root / "baseline_out",
                                      inputs={"hypotheses_resource": "hypotheses"},
                                      resources={"hypotheses": {"type": "hypotheses", "path": "hypotheses.tsv"}})
            contract_path, _ = self._contract_run(self._work(root / "contract"), root / "out")
            audit = self._audit(root / "out")
            row = audit[audit.hypothesis_id.eq("h_fas_like")].iloc[0]
            self.assertEqual("increase", row["expected_mechanism_sign"])
            self.assertEqual("concordant", row["sign_concordance"])
            self.assertIn("expected_mechanism_sign=positive; qtl_effect_sign=positive", row["sign_concordance_reason"])
            self.assertEqual("PMID:360", row["qtl_effect_source_id"])
            self.assertEqual("Table 2", row["qtl_effect_source_detail"])
            self.assertEqual("field_2026", row["qtl_effect_environment_id"])
            ledger = pd.read_csv(root / "out" / "hypotheses" / "locus_hypotheses.tsv", sep="\t", dtype=str).fillna("")
            self.assertNotIn("sign_concordance", ledger.columns)
            self.assertEqual(list(ledger.columns) + AUDIT_COLUMNS, list(audit.columns))
            self.assertEqual(["h_fas_like"], list(audit["hypothesis_id"]))
            # Rank/score outputs stay byte-identical to the no-contract baseline.
            for relative in ("rank/ranked_candidates.tsv", "rank/score_contributions.tsv"):
                self.assertEqual((baseline / relative).read_bytes(), (root / "out" / relative).read_bytes())
            artifact = pd.read_csv(root / "out" / "audit" / "qtl_effect_contract.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("bundled_fas_like", artifact.iloc[0]["case_id"])
            self.assertEqual("increase", artifact.iloc[0]["effect_direction"])
            expected_fp = hashlib.sha256(contract_path.read_bytes()).hexdigest()
            for manifest_path in (root / "out" / "run_manifest.json", root / "out" / "case_state.json"):
                block = self._qtl_block(manifest_path)
                self.assertEqual("used", block["status"])
                self.assertEqual(expected_fp, block["fingerprint"])
                self.assertEqual({"concordant": 1, "incompatible": 0, "unknown": 0}, block["sign_concordance_counts"])
                self.assertFalse(block["sign_concordance_review_required"])
                self.assertEqual("none", block["rank_contribution"])
            step = json.loads((root / "out" / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("complete", step["steps"]["hypothesis_ledger"]["status"])
            self.assertEqual(expected_fp, step["steps"]["hypothesis_ledger"]["qtl_effect_contract"]["fingerprint"])

    def test_incompatible_contract_is_review_visible_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "out"
            self._contract_run(self._work(root), out,
                               contract=POSITIVE_CONTRACT.replace("effect_direction: increase",
                                                                  "effect_direction: decrease"))
            audit = self._audit(out)
            self.assertEqual("incompatible", audit.iloc[0]["sign_concordance"])
            for manifest_path in (out / "run_manifest.json", out / "case_state.json"):
                block = self._qtl_block(manifest_path)
                self.assertEqual({"concordant": 0, "incompatible": 1, "unknown": 0}, block["sign_concordance_counts"])
                self.assertTrue(block["sign_concordance_review_required"])
            # No exclusion: the candidate stays fully ranked.
            ranked = pd.read_csv(out / "rank" / "ranked_candidates.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("GeneStrong", set(ranked["gene_id"]))
            self.assertIn("top", set(ranked["dry_group"]))

    def test_unknown_contract_fixture_is_valid(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "out"
            self._contract_run(self._work(root), out,
                               contract=POSITIVE_CONTRACT.replace("  effect_direction: increase\n", ""))
            self.assertEqual("unknown", self._audit(out).iloc[0]["sign_concordance"])
            block = self._qtl_block(out / "run_manifest.json")
            self.assertEqual({"concordant": 0, "incompatible": 0, "unknown": 1}, block["sign_concordance_counts"])
            self.assertFalse(block["sign_concordance_review_required"])

    def test_resume_rerun_is_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "out"
            self._contract_run(self._work(root), out)
            ranked_bytes = (out / "rank" / "ranked_candidates.tsv").read_bytes()
            audit_bytes = (out / "hypotheses" / "locus_hypothesis_sign_audit.tsv").read_bytes()
            manifest_block = self._qtl_block(out / "run_manifest.json")
            state_block = self._qtl_block(out / "case_state.json")
            self._contract_run(self._work(root / "resume_work"), out, resume=True)
            self.assertEqual(ranked_bytes, (out / "rank" / "ranked_candidates.tsv").read_bytes())
            self.assertEqual(audit_bytes, (out / "hypotheses" / "locus_hypothesis_sign_audit.tsv").read_bytes())
            self.assertEqual(manifest_block, self._qtl_block(out / "run_manifest.json"))
            self.assertEqual(state_block, self._qtl_block(out / "case_state.json"))

    def test_contract_edit_invalidates_the_hypothesis_step_cache(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "out"
            work = self._work(root)
            self._contract_run(work, out)
            manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            first_hash = manifest["steps"]["hypothesis_ledger"]["input_hash"]
            first_audit = (out / "hypotheses" / "locus_hypothesis_sign_audit.tsv").read_bytes()
            incompatible = POSITIVE_CONTRACT.replace("effect_direction: increase", "effect_direction: decrease")
            self._contract_run(work, out, resume=True, contract=incompatible)
            manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            second_hash = manifest["steps"]["hypothesis_ledger"]["input_hash"]
            self.assertNotEqual(first_hash, second_hash)
            second_audit = (out / "hypotheses" / "locus_hypothesis_sign_audit.tsv").read_bytes()
            self.assertNotEqual(first_audit, second_audit)
            self.assertTrue(self._qtl_block(out / "run_manifest.json")["sign_concordance_review_required"])

    def test_case_id_mismatch_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._contract_run(self._work(root), root / "out",
                               contract=POSITIVE_CONTRACT.replace("case_id: bundled_fas_like", "case_id: other_case"),
                               expect_fail=True)

    def test_legacy_case_without_contract_is_all_unknown(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._work(root)
            _hypotheses(work)
            out = self._run_case(work, root / "out",
                                 inputs={"hypotheses_resource": "hypotheses"},
                                 resources={"hypotheses": {"type": "hypotheses", "path": "hypotheses.tsv"}})
            audit = self._audit(out)
            row = audit[audit.hypothesis_id.eq("h_fas_like")].iloc[0]
            self.assertEqual("increase", row["expected_mechanism_sign"])
            self.assertEqual("unknown", row["sign_concordance"])
            self.assertEqual("unknown", row["qtl_effect_source_id"])
            self.assertFalse((out / "audit").exists())
            block = self._qtl_block(out / "run_manifest.json")
            self.assertEqual("missing", block["status"])
            self.assertIsNone(block["artifact"])
            self.assertIsNone(block["fingerprint"])
            self.assertEqual({"concordant": 0, "incompatible": 0, "unknown": 1}, block["sign_concordance_counts"])
            self.assertFalse(block["sign_concordance_review_required"])
            self.assertEqual(block, self._qtl_block(out / "case_state.json"))


if __name__ == "__main__":
    unittest.main()
