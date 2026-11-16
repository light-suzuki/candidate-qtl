"""Issue #364: canonical evidence schema + case-run evidence dependency wiring.

The canonical evidence schema gains two OPTIONAL provenance columns
(evidence_kind / root_source_id). case-run keeps them through the evidence
step, reports dependency clusters in run_manifest.json / report.md, and the
reporting never touches rank/score/weight/keyword/candidate-set bytes.
"""

from __future__ import annotations

import json
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import yaml

from candidate_qtl.cli import main
from candidate_qtl.evidence_dependency import (
    annotate_evidence_dependencies,
    dependency_reporting_counts,
)
from candidate_qtl.schemas import EVIDENCE_COLUMNS

FIXTURES = Path(__file__).parent / "fixtures" / "case_run"

# 8-core + phenotype_direction header of the bundled evidence resource; the
# dependency columns are appended only in the with-columns variant.
EVIDENCE_HEADER = "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\tphenotype_direction"
ORIGINAL_ROWS = [
    ["GeneStrong", "known_qtg", "support", "high", "validated model phenotype", "literature", "PMID:12345", "Source-backed genetic model evidence", "increase"],
    ["GeneDark", "annotation_incomplete", "rescue", "medium", "no homology", "annotation", "fixture:gff3", "Short complete ORF without homology", "unknown"],
]
ORIGINAL_DEP = [("primary_observation", "ROOT1"), ("curated_assertion", "ROOT2")]

SECOND_HEADER = "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription"
SECOND_ROWS = [
    ["GeneStrong", "literature_curated", "support", "high", "curated model statement", "publication", "PMID:1001", "ROOT1 primary observation"],
    ["GeneDark", "orthology_transfer", "support", "medium", "ortholog function", "orthology", "OMAC:7", "ROOT1 orthology transfer"],
    ["GeneStrong", "sequence_domain_call", "support", "medium", "conserved domain", "interpro", "IPR:1", "ROOT2 domain call"],
    ["GeneDark", "annotation_incomplete", "rescue", "medium", "no homology", "annotation", "fixture:gff3", "unrooted row"],
]
SECOND_DEP = [
    ("primary_observation", "ROOT1"),
    ("orthology_transfer", "ROOT1"),
    ("sequence_domain_call", "ROOT2"),
    ("unknown", ""),
]


class EvidenceDependencyNullableDtypeTests(unittest.TestCase):
    """Direct-API contract: nullable-dtype frames (pd.NA cells) stay clean."""

    def test_nullable_dtype_frame_annotates_without_na_strings(self) -> None:
        evidence = pd.DataFrame({
            "source_id": ["S1", "S2", "S3"],
            "evidence_kind": pd.array(["primary_observation", "curated_assertion", "prediction"], dtype="string"),
            "root_source_id": pd.array(["ROOT1", "ROOT1", pd.NA], dtype="string"),
        })
        annotated = annotate_evidence_dependencies(evidence)
        self.assertFalse(annotated.astype(str).apply(lambda column: column.str.contains("<NA>")).any().any())
        self.assertEqual("", annotated.iloc[2]["root_source_id"])
        self.assertEqual("dependency_unknown", annotated.iloc[2]["dependency_status"])
        counts = dependency_reporting_counts(evidence)
        self.assertEqual(3, counts["evidence_count"])
        self.assertEqual(1, counts["independent_evidence_cluster_count"])
        self.assertEqual(1, counts["dependency_unknown_evidence_count"])
        self.assertTrue(counts["dependency_unknown_present"])


class CaseRunEvidenceDependencyWiringTests(unittest.TestCase):
    def _make_work(self, root: Path, with_dependency_columns: bool = True,
                   original_dep: list[tuple[str, str]] | None = None) -> Path:
        work = root / "case"
        shutil.copytree(FIXTURES, work)
        original_dep = ORIGINAL_DEP if original_dep is None else original_dep
        original_rows = []
        for row, (kind, root_source) in zip(ORIGINAL_ROWS, original_dep):
            original_rows.append(row + ([kind, root_source] if with_dependency_columns else []))
        header = EVIDENCE_HEADER + ("\tevidence_kind\troot_source_id" if with_dependency_columns else "")
        (work / "evidence.tsv").write_text(
            header + "\n" + "\n".join("\t".join(row) for row in original_rows) + "\n",
            encoding="utf-8")
        second_rows = []
        for row, (kind, root_source) in zip(SECOND_ROWS, SECOND_DEP):
            second_rows.append(row + ([kind, root_source] if with_dependency_columns else []))
        header2 = SECOND_HEADER + ("\tevidence_kind\troot_source_id" if with_dependency_columns else "")
        (work / "dependency_evidence.tsv").write_text(
            header2 + "\n" + "\n".join("\t".join(row) for row in second_rows) + "\n",
            encoding="utf-8")
        return work

    def _wire_case(self, work: Path) -> None:
        case = yaml.safe_load((work / "case.yaml").read_text(encoding="utf-8"))
        case["inputs"]["evidence_resources"] = ["evidence", "dependency_evidence"]
        case_path = work / "case_run.yaml"
        case_path.write_text(yaml.safe_dump(case, sort_keys=False), encoding="utf-8")
        registry = yaml.safe_load((work / "resources.yaml").read_text(encoding="utf-8"))
        registry["resources"]["dependency_evidence"] = {"type": "canonical_evidence", "path": "dependency_evidence.tsv"}
        (work / "resources.yaml").write_text(yaml.safe_dump(registry, sort_keys=False), encoding="utf-8")

    def _run_case(self, work: Path, out: Path, *, resume: bool = False) -> Path:
        args = ["case-run", "--case", str(work / "case_run.yaml"),
                "--resources", str(work / "resources.yaml"), "--out-dir", str(out)]
        if resume:
            args.append("--resume")
        self.assertEqual(0, main(args))
        return out

    def _run_with_dep(self, root: Path, out_name: str = "out", *, resume: bool = False) -> Path:
        work = self._make_work(root)
        self._wire_case(work)
        return self._run_case(work, root / out_name, resume=resume)

    def test_case_run_reports_dependency_clusters_and_keeps_columns(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run_with_dep(Path(tmp))
            evidence_out = pd.read_csv(out / "evidence.tsv", sep="\t", dtype=str).fillna("")
            # 8 core columns first, then the optional provenance columns.
            self.assertEqual(EVIDENCE_COLUMNS, list(evidence_out.columns[:8]))
            self.assertEqual(["evidence_kind", "root_source_id"], list(evidence_out.columns[8:]))
            self.assertEqual(6, len(evidence_out))
            # Run manifest: the evidence step record carries the reporting
            # counts beside status/input_hash/outputs (complete() flattens
            # extra into the step dict).
            manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            evidence_step = manifest["steps"]["evidence"]
            self.assertEqual(6, evidence_step["evidence_count"])
            self.assertEqual(2, evidence_step["independent_evidence_cluster_count"])
            self.assertEqual(1, evidence_step["dependency_unknown_evidence_count"])
            self.assertTrue(evidence_step["dependency_unknown_present"])
            self.assertTrue(evidence_step["evidence_dependency_manifest"].endswith("evidence_dependency_manifest.json"))
            self.assertTrue(evidence_step["evidence_dependency_clusters"].endswith("evidence_dependency_clusters.tsv"))
            # Top-level block: reporting-only contract.
            block = manifest["evidence_dependency"]
            self.assertEqual("none", block["rank_contribution"])
            self.assertEqual("evaluated", block["dependency_status"])
            self.assertEqual("evidence_dependency_clusters.tsv", block["clusters"])
            self.assertEqual("evidence_dependency_manifest.json", block["manifest"])
            self.assertEqual(2, len(block["input_fingerprint"]))
            self.assertIn("never contribute to the candidate rank", block["note"])
            # Cluster table: ROOT1(3) / ROOT2(2) / unknown(1) -> exactly 3 rows.
            clusters = pd.read_csv(out / "evidence_dependency_clusters.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(3, len(clusters))
            self.assertEqual(3, int(clusters.set_index("root_source_id").loc["ROOT1", "evidence_count"]))
            self.assertEqual(2, int(clusters.set_index("root_source_id").loc["ROOT2", "evidence_count"]))
            unknown = clusters[clusters.dependency_status.eq("dependency_unknown")].iloc[0]
            self.assertEqual(1, int(unknown["evidence_count"]))
            # Dependency manifest mirrors the counts.
            dep_manifest = json.loads((out / "evidence_dependency_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("evaluated", dep_manifest["dependency_status"])
            self.assertEqual(6, dep_manifest["evidence_count"])
            self.assertEqual("none", dep_manifest["rank_contribution"])
            self.assertEqual(2, len(dep_manifest["input_fingerprint"]))
            # report.md carries the section and the cluster table.
            report = (out / "report.md").read_text(encoding="utf-8")
            self.assertIn("## Evidence dependency", report)
            self.assertIn("- Evidence rows: 6", report)
            self.assertIn("- Known independent dependency clusters: 2", report)
            self.assertIn("- Dependency-unknown rows: 1 (present)", report)
            self.assertIn("| dependency_cluster | dependency_status | root_source_id |", report)
            self.assertIn("dependency_unknown", report)

    def test_dependency_columns_never_change_rank_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with_dep = self._run_with_dep(root, "out_with_dep")
            without_work = self._make_work(root / "plain", with_dependency_columns=False)
            self._wire_case(without_work)
            without_dep = self._run_case(without_work, root / "out_without_dep")
            for name in ("ranked_candidates.tsv", "score_contributions.tsv"):
                self.assertEqual(
                    (without_dep / "rank" / name).read_bytes(),
                    (with_dep / "rank" / name).read_bytes(),
                )

    def test_identical_inputs_produce_byte_identical_dependency_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._make_work(root)
            self._wire_case(work)
            out_a = self._run_case(work, root / "out_a")
            out_b = self._run_case(work, root / "out_b")
            # Pure computations over the same inputs: byte-identical across
            # separate out-dirs (they embed no out-dir paths).
            for name in ("evidence_dependency_clusters.tsv", "evidence_dependency_manifest.json"):
                self.assertEqual((out_a / name).read_bytes(), (out_b / name).read_bytes())
            # run_manifest embeds absolute artifact paths and mtime-based step
            # hashes, so byte-identity is only defined for the same out-dir: a
            # resumed recomposition must reproduce it byte-for-byte.
            first_manifest = (out_a / "run_manifest.json").read_bytes()
            self._run_case(work, out_a, resume=True)
            self.assertEqual(first_manifest, (out_a / "run_manifest.json").read_bytes())

    def test_resume_reuses_evidence_step_and_column_change_invalidates(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._make_work(root)
            self._wire_case(work)
            out = self._run_case(work, root / "out")
            first_hash = json.loads((out / "case_state.json").read_text(encoding="utf-8"))["steps"]["evidence"]["input_hash"]
            first_clusters = (out / "evidence_dependency_clusters.tsv").read_bytes()
            # Unchanged out-dir + --resume: the evidence step is reused, so the
            # dependency computation never runs again.
            with patch("candidate_qtl.workflow_extensions.dependency_reporting_counts",
                       side_effect=AssertionError("evidence step rerun on resume")) as counted:
                self._run_case(work, out, resume=True)
                counted.assert_not_called()
            state = json.loads((out / "case_state.json").read_text(encoding="utf-8"))
            self.assertEqual(first_hash, state["steps"]["evidence"]["input_hash"])
            self.assertEqual(first_clusters, (out / "evidence_dependency_clusters.tsv").read_bytes())
            # Changing the columns of an input resource changes the evidence
            # content, so the digest invalidates and the step recomputes.
            (work / "dependency_evidence.tsv").write_text(
                SECOND_HEADER + "\n" + "\n".join("\t".join(row) for row in SECOND_ROWS) + "\n",
                encoding="utf-8")
            self._run_case(work, out, resume=True)
            state = json.loads((out / "case_state.json").read_text(encoding="utf-8"))
            self.assertNotEqual(first_hash, state["steps"]["evidence"]["input_hash"])
            self.assertNotEqual(first_clusters, (out / "evidence_dependency_clusters.tsv").read_bytes())
            # All four rows now lack provenance -> unknown bucket grows to 4.
            clusters = pd.read_csv(out / "evidence_dependency_clusters.tsv", sep="\t", dtype=str).fillna("")
            unknown = clusters[clusters.dependency_status.eq("dependency_unknown")].iloc[0]
            self.assertEqual(4, int(unknown["evidence_count"]))

    def test_evidence_row_backward_compatibility(self) -> None:
        from candidate_qtl.evidence_io import evidence_row
        core = evidence_row("G1", "known_qtg", "support", "high", 1, "literature", "PMID:1", "d")
        self.assertEqual(EVIDENCE_COLUMNS, list(core))
        extended = evidence_row("G1", "known_qtg", "support", "high", 1, "literature", "PMID:1", "d",
                                evidence_kind="primary_observation", root_source_id="ROOT1")
        self.assertEqual(EVIDENCE_COLUMNS + ["evidence_kind", "root_source_id"], list(extended))
        partial = evidence_row("G1", "known_qtg", "support", "high", 1, "literature", "PMID:1", "d",
                               evidence_kind="primary_observation")
        self.assertEqual(EVIDENCE_COLUMNS + ["evidence_kind"], list(partial))

    def test_space_form_evidence_kind_is_normalized_before_write(self) -> None:
        # The documented-acceptable space-separated form must be normalized in
        # the evidence step, BEFORE evidence.tsv is written: the rank step
        # validates the strict EVIDENCE_DEPENDENCY_KINDS enum and would crash
        # on the raw value (adversarial review fix, F1).
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._make_work(root, original_dep=[("primary observation", "ROOT1"), ("curated_assertion", "ROOT2")])
            self._wire_case(work)
            out = self._run_case(work, root / "out")
            evidence_out = pd.read_csv(out / "evidence.tsv", sep="\t", dtype=str).fillna("")
            self.assertNotIn("primary observation", set(evidence_out["evidence_kind"]))
            strong = evidence_out[evidence_out.gene_id.eq("GeneStrong")]
            self.assertIn("primary_observation", set(strong["evidence_kind"]))
            # The run completed with the usual counts.
            manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            block = manifest["evidence_dependency"]
            self.assertEqual("evaluated", block["dependency_status"])
            self.assertEqual(6, block["evidence_count"])
            self.assertEqual(2, block["independent_evidence_cluster_count"])
            self.assertEqual(1, block["dependency_unknown_evidence_count"])

    def test_invalid_evidence_kind_provider_failed_path(self) -> None:
        # An out-of-vocabulary evidence_kind must not crash the case: it
        # degrades to the evidence_kind_invalid review item, the on-disk
        # dependency manifest says not_evaluated with zero counts, and the
        # run_manifest block (rebuilt from that manifest) never shows stale
        # counts (adversarial review fixes F1/F2/F4/F5).
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            work = self._make_work(root, original_dep=[("invalid_kind", "ROOT1"), ("curated_assertion", "ROOT2")])
            self._wire_case(work)
            out = self._run_case(work, root / "out")
            # The placeholder evidence.tsv keeps the attempted column shape
            # (8 core columns + the optional columns present in the inputs).
            placeholder = pd.read_csv(out / "evidence.tsv", sep="\t", nrows=0)
            self.assertEqual(EVIDENCE_COLUMNS + ["evidence_kind", "root_source_id"], list(placeholder.columns))
            dep_manifest = json.loads((out / "evidence_dependency_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("not_evaluated", dep_manifest["dependency_status"])
            self.assertEqual(0, dep_manifest["evidence_count"])
            manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            block = manifest["evidence_dependency"]
            self.assertEqual("not_evaluated", block["dependency_status"])
            self.assertEqual(0, block["evidence_count"])
            self.assertEqual(0, block["independent_evidence_cluster_count"])
            self.assertEqual(0, block["dependency_unknown_evidence_count"])
            self.assertFalse(block["dependency_unknown_present"])
            self.assertEqual("none", block["rank_contribution"])
            report = (out / "report.md").read_text(encoding="utf-8")
            self.assertIn("- Status: not_evaluated (evidence provider failed)", report)
            review = pd.read_csv(out / "review_required.tsv", sep="\t", dtype=str).fillna("")
            invalid = review[review.reason.eq("evidence_kind_invalid")]
            self.assertEqual(1, len(invalid))
            row = invalid.iloc[0]
            self.assertIn("fix evidence_kind values", row["required_next_action"])
            self.assertIn("primary_observation", row["required_next_action"])
            self.assertIn("keep the rank lane positional-fallback", row["allowed_choices"])

    def test_manifest_clusters_and_report_counts_are_consistent(self) -> None:
        # Cross-surface consistency: the run_manifest evidence_dependency
        # block (rebuilt from the on-disk manifest) must equal the cluster
        # table totals and the report.md numbers.
        with tempfile.TemporaryDirectory() as tmp:
            out = self._run_with_dep(Path(tmp))
            manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
            block = manifest["evidence_dependency"]
            clusters = pd.read_csv(out / "evidence_dependency_clusters.tsv", sep="\t", dtype=str).fillna("")
            total_rows = int(clusters["evidence_count"].astype(int).sum())
            independent = int(clusters.dependency_status.eq("root_source_known").sum())
            unknown_rows = int(clusters.loc[clusters.dependency_status.eq("dependency_unknown"), "evidence_count"].astype(int).sum())
            self.assertEqual(block["evidence_count"], total_rows)
            self.assertEqual(block["independent_evidence_cluster_count"], independent)
            self.assertEqual(block["dependency_unknown_evidence_count"], unknown_rows)
            self.assertEqual(block["dependency_unknown_present"], unknown_rows > 0)
            report = (out / "report.md").read_text(encoding="utf-8")
            self.assertIn(f"- Evidence rows: {block['evidence_count']}", report)
            self.assertIn(f"- Known independent dependency clusters: {block['independent_evidence_cluster_count']}", report)
            present = "present" if block["dependency_unknown_present"] else "absent"
            self.assertIn(f"- Dependency-unknown rows: {block['dependency_unknown_evidence_count']} ({present})", report)


if __name__ == "__main__":
    unittest.main()
