"""Candidate elimination / hypothesis survival layer (issue #171)."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl import elimination as elimination_module
from candidate_qtl.cli import main


class RuleRegistryTests(unittest.TestCase):
    def test_default_registry_loads_and_ships_in_the_wheel(self) -> None:
        path = elimination_module.default_rules_path()
        self.assertTrue(path.exists())
        rules, meta = elimination_module.load_rules()
        self.assertGreater(len(rules), 10)
        self.assertIn("sha256", meta)
        repo_root = Path(__file__).resolve().parents[1]
        vendored = repo_root / "src" / "candidate_qtl" / "resources" / "elimination_rules.tsv"
        self.assertEqual((repo_root / "templates" / "elimination_rules.tsv").read_bytes(), vendored.read_bytes())

    def _registry(self, root: Path, row: str) -> Path:
        path = root / "rules.tsv"
        path.write_text(
            "elimination_rule_id\tscope\tevidence_code\tminimum_quality\tminimum_coverage\tstrength\t"
            "default_transition\trequired_review\trationale\tregistry_version\n" + row, encoding="utf-8")
        return path

    def test_a_non_hard_rule_may_not_exclude_a_gene(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = self._registry(Path(tmp), "r1\tgene\tx\tany\t0.0\tmoderate\texcluded_by_genetics\tfalse\twishful\t1\n")
            with self.assertRaises(ValueError) as context:
                elimination_module.load_rules(str(path))
            self.assertIn("only a hard rule", str(context.exception))

    def test_a_non_evidence_rule_may_not_transition_anywhere_meaningful(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = self._registry(Path(tmp), "r1\tcoding\tx\tany\t0.0\tnone\tcoding_rejected\tfalse\twishful\t1\n")
            with self.assertRaises(ValueError) as context:
                elimination_module.load_rules(str(path))
            self.assertIn("non-evidence rule", str(context.exception))

    def test_a_mechanism_transition_requires_a_mechanism_scope(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = self._registry(Path(tmp), "r1\tgene\tx\tany\t0.0\tstrong\tcoding_rejected\tfalse\twishful\t1\n")
            with self.assertRaises(ValueError):
                elimination_module.load_rules(str(path))


class EliminationTests(unittest.TestCase):
    def _fixture(self, root: Path) -> dict[str, Path]:
        paths = {name: root / name for name in ("candidates.tsv", "evidence.tsv", "hypotheses.tsv")}
        paths["candidates.tsv"].write_text(
            "gene_id\tdry_group\nG_OUT\ttop\nG_CODING\ttop\nG_WEAK\treserve\nG_FAILED\treserve\nG_KEEP\ttop\n", encoding="utf-8")
        paths["hypotheses.tsv"].write_text(
            "hypothesis_id\thypothesis_type\tmember_gene_ids\n"
            "h_coding\tcoding\tG_CODING\n"
            "h_reg\tregulatory\tG_CODING\n", encoding="utf-8")
        paths["evidence.tsv"].write_text("\n".join([
            "gene_id\thypothesis_id\tevidence_code\tevidence_id\tquality\tcoverage",
            # Hard positional exclusion.
            "G_OUT\t\tinformative_recombinant_excludes\tREC:1\thigh\t1.0",
            # Full CDS read at full coverage: rejects the coding hypothesis only.
            "G_CODING\th_coding\tfull_cds_no_parental_difference\tSANGER:1\thigh\t1.0",
            # Moderate counterevidence: weakens, never excludes.
            "G_WEAK\t\talternate_paralog_strongly_supported\tORTH:1\tmedium\t0.0",
            # Non-evidence: an assay that did not run.
            "G_FAILED\t\tassay_failed\tPCR:1\tany\t0.0",
        ]) + "\n", encoding="utf-8")
        return paths

    def _run(self, root: Path, extra: list[str] | None = None) -> Path:
        fixture = self._fixture(root)
        out = root / "elim"
        self.assertEqual(0, main(["candidate-elimination", "--candidates", str(fixture["candidates.tsv"]),
                                  "--evidence", str(fixture["evidence.tsv"]),
                                  "--hypotheses", str(fixture["hypotheses.tsv"]),
                                  "--out-dir", str(out), "--case-id", "demo"] + (extra or [])))
        return out

    def test_hard_positional_evidence_excludes_the_gene(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            survival = pd.read_csv(out / "candidate_survival.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual("excluded_by_recombinant", survival.loc["G_OUT", "current_gene_state"])
            self.assertEqual("REC:1", survival.loc["G_OUT", "strongest_exclusion_evidence_id"])
            remaining = pd.read_csv(out / "remaining_candidate_set.tsv", sep="\t", dtype=str).fillna("")
            self.assertNotIn("G_OUT", set(remaining["candidate_gene_id"]))

    def test_coding_negative_rejects_the_mechanism_not_the_gene(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            mechanisms = pd.read_csv(out / "hypothesis_survival.tsv", sep="\t", dtype=str).fillna("")
            coding = mechanisms[(mechanisms.candidate_gene_id == "G_CODING") & (mechanisms.hypothesis_id == "h_coding")].iloc[0]
            self.assertEqual("coding_rejected", coding.mechanism_state)
            regulatory = mechanisms[(mechanisms.candidate_gene_id == "G_CODING") & (mechanisms.hypothesis_id == "h_reg")].iloc[0]
            self.assertEqual("regulatory_survives", regulatory.mechanism_state)
            survival = pd.read_csv(out / "candidate_survival.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            # The gene keeps its other branches.
            self.assertEqual("survives", survival.loc["G_CODING", "current_gene_state"])
            self.assertEqual("0", survival.loc["G_CODING", "hard_exclusion_count"])
            remaining = pd.read_csv(out / "remaining_candidate_set.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertIn("regulatory", remaining.loc["G_CODING", "surviving_mechanisms"])

    def test_moderate_counterevidence_only_weakens(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            survival = pd.read_csv(out / "candidate_survival.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual("weakened", survival.loc["G_WEAK", "current_gene_state"])
            self.assertIn("alternate_paralog_strongly_supported", survival.loc["G_WEAK", "weakened_reason"])
            remaining = pd.read_csv(out / "remaining_candidate_set.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("G_WEAK", set(remaining["candidate_gene_id"]))

    def test_assay_failure_never_removes_a_candidate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            survival = pd.read_csv(out / "candidate_survival.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual("unresolved", survival.loc["G_FAILED", "current_gene_state"])
            events = pd.read_csv(out / "elimination_events.tsv", sep="\t", dtype=str).fillna("")
            failure = events[events.candidate_gene_id == "G_FAILED"].iloc[0]
            self.assertEqual("none", failure.strength)
            self.assertEqual("technical_failure", failure.new_state)
            remaining = pd.read_csv(out / "remaining_candidate_set.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("G_FAILED", set(remaining["candidate_gene_id"]))

    def test_partial_coverage_cannot_reject_a_coding_hypothesis(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._fixture(root)
            fixture["evidence.tsv"].write_text(
                "gene_id\thypothesis_id\tevidence_code\tevidence_id\tquality\tcoverage\n"
                "G_CODING\th_coding\tfull_cds_no_parental_difference\tSANGER:1\thigh\t0.6\n", encoding="utf-8")
            out = root / "elim"
            self.assertEqual(0, main(["candidate-elimination", "--candidates", str(fixture["candidates.tsv"]),
                                      "--evidence", str(fixture["evidence.tsv"]),
                                      "--hypotheses", str(fixture["hypotheses.tsv"]), "--out-dir", str(out)]))
            mechanisms = pd.read_csv(out / "hypothesis_survival.tsv", sep="\t", dtype=str).fillna("")
            coding = mechanisms[(mechanisms.candidate_gene_id == "G_CODING") & (mechanisms.hypothesis_id == "h_coding")].iloc[0]
            self.assertEqual("coding_survives", coding.mechanism_state)
            events = pd.read_csv(out / "elimination_events.tsv", sep="\t", dtype=str).fillna("")
            # The strict rule's coverage requirement was not met, so the registry's
            # partial-coverage counterpart fires instead and eliminates nothing.
            self.assertEqual("partial_cds_no_parental_difference", events.iloc[0].rule_id)
            self.assertEqual("none", events.iloc[0].strength)
            self.assertEqual("coverage_incomplete", events.iloc[0].new_state)
            self.assertEqual("review_required", events.iloc[0].review_status)

    def test_low_quality_cosegregation_does_not_exclude(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._fixture(root)
            fixture["evidence.tsv"].write_text(
                "gene_id\thypothesis_id\tevidence_code\tevidence_id\tquality\tcoverage\n"
                "G_KEEP\t\tcosegregation_discordant\tF2:1\tlow\t0.0\n", encoding="utf-8")
            out = root / "elim"
            self.assertEqual(0, main(["candidate-elimination", "--candidates", str(fixture["candidates.tsv"]),
                                      "--evidence", str(fixture["evidence.tsv"]), "--out-dir", str(out)]))
            survival = pd.read_csv(out / "candidate_survival.tsv", sep="\t", dtype=str).fillna("").set_index("candidate_gene_id")
            self.assertEqual("unresolved", survival.loc["G_KEEP", "current_gene_state"])
            remaining = pd.read_csv(out / "remaining_candidate_set.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("G_KEEP", set(remaining["candidate_gene_id"]))

    def test_history_is_append_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = self._run(root)
            events_before = len(pd.read_csv(first / "elimination_events.tsv", sep="\t", dtype=str))
            second = self._run(root, extra=["--previous-history", str(first / "elimination_events.tsv")])
            events_after = len(pd.read_csv(second / "elimination_events.tsv", sep="\t", dtype=str))
            self.assertEqual(events_before * 2, events_after)

    def test_manifest_reports_compression_and_false_exclusion_risk(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            manifest = json.loads((out / "elimination_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(5, manifest["initial_candidates"])
            self.assertEqual(4, manifest["remaining_candidates"])
            self.assertEqual(1, manifest["hard_excluded_count"])
            self.assertEqual(1, manifest["mechanism_only_rejected_count"])
            self.assertTrue(manifest["false_exclusion_risk_flags"])
            self.assertIn("hard", manifest["strength_contract"])
            self.assertIn("only hard positional or genetic evidence may exclude a gene",
                          " ".join(manifest["safety"]))
            self.assertIn("discriminate the most remaining candidates", manifest["next_experiment_hint"])

    def _candidates(self, root: Path, count: int) -> Path:
        path = root / f"candidates_{count}.tsv"
        path.write_text("gene_id\tdry_group\n" + "".join(f"G{i:03d}\ttop\n" for i in range(count)), encoding="utf-8")
        return path

    def _exclusion_evidence(self, root: Path, gene_ids: list[str]) -> Path:
        path = root / "evidence_exclude.tsv"
        path.write_text("gene_id\thypothesis_id\tevidence_code\tevidence_id\tquality\tcoverage\n" + "".join(
            f"{gene}\t\tinformative_recombinant_excludes\tREC:{gene}\thigh\t1.0\n" for gene in gene_ids), encoding="utf-8")
        return path

    def test_manifest_potentially_eliminated_is_eliminated_not_remaining(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # 20 candidates, 16 removed by hard positional evidence -> 4 remaining.
            candidates = self._candidates(root, 20)
            evidence = self._exclusion_evidence(root, [f"G{i:03d}" for i in range(16)])
            out = root / "elim20"
            self.assertEqual(0, main(["candidate-elimination", "--candidates", str(candidates),
                                      "--evidence", str(evidence), "--out-dir", str(out)]))
            manifest = json.loads((out / "elimination_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(20, manifest["initial_candidates"])
            self.assertEqual(4, manifest["remaining_candidates"])
            self.assertEqual(16, manifest["candidates_potentially_eliminated"])
            self.assertEqual(4, len(pd.read_csv(out / "remaining_candidate_set.tsv", sep="\t", dtype=str)))

    def test_manifest_potentially_eliminated_is_zero_when_nothing_is_eliminated(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = self._candidates(root, 5)
            out = root / "elim5"
            self.assertEqual(0, main(["candidate-elimination", "--candidates", str(candidates), "--out-dir", str(out)]))
            manifest = json.loads((out / "elimination_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(5, manifest["initial_candidates"])
            self.assertEqual(5, manifest["remaining_candidates"])
            self.assertEqual(0, manifest["candidates_potentially_eliminated"])

    def test_manifest_potentially_eliminated_covers_all_when_nothing_remains(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = self._candidates(root, 3)
            evidence = self._exclusion_evidence(root, ["G000", "G001", "G002"])
            out = root / "elim0"
            self.assertEqual(0, main(["candidate-elimination", "--candidates", str(candidates),
                                      "--evidence", str(evidence), "--out-dir", str(out)]))
            manifest = json.loads((out / "elimination_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(3, manifest["initial_candidates"])
            self.assertEqual(0, manifest["remaining_candidates"])
            self.assertEqual(3, manifest["candidates_potentially_eliminated"])
            self.assertEqual(0, len(pd.read_csv(out / "remaining_candidate_set.tsv", sep="\t", dtype=str)))

    def test_no_evidence_leaves_everything_untested_and_flagged(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = self._fixture(root)
            out = root / "elim"
            self.assertEqual(0, main(["candidate-elimination", "--candidates", str(fixture["candidates.tsv"]),
                                      "--out-dir", str(out)]))
            survival = pd.read_csv(out / "candidate_survival.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(survival["current_gene_state"].eq("not_tested").all())
            self.assertEqual(5, len(pd.read_csv(out / "remaining_candidate_set.tsv", sep="\t", dtype=str)))
            reasons = set(pd.read_csv(out / "elimination_review_required.tsv", sep="\t", dtype=str).fillna("")["reason"])
            self.assertIn("no_elimination_evidence_supplied", reasons)

    def test_layer_writes_no_score(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run(Path(tmp))
            for path in sorted(out.glob("*.tsv")):
                columns = set(pd.read_csv(path, sep="\t", dtype=str).columns)
                self.assertFalse({"dry_score", "score", "top_eligibility"} & columns, path.name)
            survival = pd.read_csv(out / "candidate_survival.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(survival["rank_contribution"].eq("none").all())


if __name__ == "__main__":
    unittest.main()
