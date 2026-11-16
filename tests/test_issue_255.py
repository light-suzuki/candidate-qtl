from __future__ import annotations

import hashlib
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl import workflow_extensions
from candidate_qtl.cli import main


class Issue255ContentDigestTests(unittest.TestCase):
    def _evidence_case(self, root: Path) -> tuple[Path, Path, Path, Path]:
        candidates = root / "candidates.tsv"
        candidates.write_text("gene_id\nGeneA\nGeneB\n", encoding="utf-8")
        evidence = root / "evidence.tsv"
        evidence.write_text(
            "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\n"
            "GeneA\tknown_qtg\tsupport\thigh\tmodel\tliterature\tPMID:1\tgenetic model\n"
            "GeneB\tknown_qtg\tsupport\tlow\tmodel\tliterature\tPMID:2\tweak model\n",
            encoding="utf-8",
        )
        case = root / "case.yaml"
        resources = root / "resources.yaml"
        case.write_text(
            "case_id: content_digest\n"
            "trait:\n  phenotype_direction: increase\n"
            "source:\n  source_id: PMID:QTL\n  target_interval: chr1:1-10\n  literature_search_status: complete\n"
            "inputs:\n  candidates_resource: candidates\n  evidence_resources: [evidence]\n"
            "assay:\n  top_n: 2\n",
            encoding="utf-8",
        )
        resources.write_text(
            f"resources:\n  candidates: {{path: {candidates}}}\n  evidence: {{path: {evidence}}}\n",
            encoding="utf-8",
        )
        return case, resources, candidates, evidence

    def test_resume_rebuilds_damaged_outputs_and_reuses_healthy_steps(self) -> None:
        for strict in (False, True):
            for damage in ("same_size_same_mtime", "truncated", "missing_companion"):
                with self.subTest(strict=strict, damage=damage), tempfile.TemporaryDirectory() as tmp:
                    root = Path(tmp)
                    case, resources, _, _ = self._evidence_case(root)
                    out = root / "out"
                    args = ["case-run", "--case", str(case), "--resources", str(resources), "--out-dir", str(out)]
                    if strict:
                        args.append("--strict-fingerprints")
                    self.assertEqual(0, main(args))
                    before = json.loads((out / "case_state.json").read_text(encoding="utf-8"))
                    rank_outputs = {
                        Path(path): Path(path).read_bytes()
                        for path in before["steps"]["rank"]["outputs"]
                    }
                    healthy_outputs = {
                        Path(path): (Path(path).read_bytes(), Path(path).stat().st_mtime_ns)
                        for name in ("coordinate", "extract", "candidate_flags", "evidence")
                        for path in before["steps"][name]["outputs"]
                    }
                    ranked = out / "rank" / "ranked_candidates.tsv"
                    if damage == "same_size_same_mtime":
                        stat = ranked.stat()
                        damaged = rank_outputs[ranked].replace(b"\t4.0\t", b"\t9.0\t", 1)
                        self.assertNotEqual(rank_outputs[ranked], damaged)
                        ranked.write_bytes(damaged)
                        os.utime(ranked, ns=(stat.st_atime_ns, stat.st_mtime_ns))
                        self.assertEqual(stat.st_size, ranked.stat().st_size)
                        self.assertEqual(stat.st_mtime_ns, ranked.stat().st_mtime_ns)
                    elif damage == "truncated":
                        ranked.write_bytes(b"")
                    else:
                        (out / "rank" / "score_contributions.tsv").unlink()

                    with patch.object(workflow_extensions, "prioritize_evidence", wraps=workflow_extensions.prioritize_evidence) as ranker:
                        self.assertEqual(0, main(args + ["--resume"]))
                        self.assertEqual(1, ranker.call_count)
                    after = json.loads((out / "case_state.json").read_text(encoding="utf-8"))
                    self.assertEqual("REPORT_READY", after["state"])
                    self.assertEqual(before["steps"]["rank"]["input_hash"], after["steps"]["rank"]["input_hash"])
                    for path, expected in rank_outputs.items():
                        self.assertEqual(expected, path.read_bytes())
                    for path, expected in healthy_outputs.items():
                        self.assertEqual(expected, (path.read_bytes(), path.stat().st_mtime_ns))

                    rank_mtimes = {path: path.stat().st_mtime_ns for path in rank_outputs}
                    with patch.object(workflow_extensions, "prioritize_evidence", wraps=workflow_extensions.prioritize_evidence) as ranker:
                        self.assertEqual(0, main(args + ["--resume"]))
                        ranker.assert_not_called()
                    self.assertEqual(rank_mtimes, {path: path.stat().st_mtime_ns for path in rank_outputs})

    def test_legacy_manifest_rebuilds_once_before_reusing_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case, resources, _, _ = self._evidence_case(root)
            out = root / "out"
            args = ["case-run", "--case", str(case), "--resources", str(resources), "--out-dir", str(out), "--strict-fingerprints"]
            self.assertEqual(0, main(args))
            state_path = out / "case_state.json"
            legacy = json.loads(state_path.read_text(encoding="utf-8"))
            for step in legacy["steps"].values():
                step.pop("output_fingerprints", None)
            state_path.write_text(json.dumps(legacy), encoding="utf-8")
            ranked = out / "rank" / "ranked_candidates.tsv"
            original = ranked.read_bytes()
            stat = ranked.stat()
            ranked.write_bytes(original.replace(b"\t4.0\t", b"\t9.0\t", 1))
            self.assertNotEqual(original, ranked.read_bytes())
            os.utime(ranked, ns=(stat.st_atime_ns, stat.st_mtime_ns))

            with patch.object(workflow_extensions, "prioritize_evidence", wraps=workflow_extensions.prioritize_evidence) as ranker:
                self.assertEqual(0, main(args + ["--resume"]))
                self.assertEqual(1, ranker.call_count)
            self.assertEqual(original, ranked.read_bytes())
            rebuilt = json.loads(state_path.read_text(encoding="utf-8"))
            self.assertTrue(all("output_fingerprints" in step for step in rebuilt["steps"].values()))
            rebuilt_manifest = (out / "run_manifest.json").read_bytes()
            with patch.object(workflow_extensions, "prioritize_evidence", wraps=workflow_extensions.prioritize_evidence) as ranker:
                self.assertEqual(0, main(args + ["--resume"]))
                ranker.assert_not_called()
            self.assertEqual(rebuilt_manifest, (out / "run_manifest.json").read_bytes())

    def test_resume_detects_same_size_same_mtime_content_change(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case, resources, _, evidence = self._evidence_case(root)
            out = root / "out"
            self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources), "--out-dir", str(out)]))
            state = json.loads((out / "case_state.json").read_text(encoding="utf-8"))
            evidence_hash_before = state["steps"]["evidence"]["input_hash"]
            evidence_before = (out / "evidence.tsv").read_bytes()
            ranked_before = (out / "rank" / "ranked_candidates.tsv").read_bytes()
            stat_before = evidence.stat()
            evidence.write_text(
                "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\n"
                "GeneA\tknown_qtg\tsupport\tlow\tmodel\tliterature\tPMID:3\tweak model\n"
                "GeneB\tknown_qtg\tsupport\thigh\tmodel\tliterature\tPMID:2\tgenetic model\n",
                encoding="utf-8",
            )
            os.utime(evidence, ns=(stat_before.st_atime_ns, stat_before.st_mtime_ns))
            self.assertEqual(stat_before.st_size, evidence.stat().st_size)
            self.assertEqual(stat_before.st_mtime_ns, evidence.stat().st_mtime_ns)
            self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources), "--out-dir", str(out), "--resume"]))
            state = json.loads((out / "case_state.json").read_text(encoding="utf-8"))
            self.assertNotEqual(evidence_hash_before, state["steps"]["evidence"]["input_hash"])
            self.assertNotEqual(evidence_before, (out / "evidence.tsv").read_bytes())
            self.assertIn("PMID:3", (out / "evidence.tsv").read_text())
            self.assertNotEqual(ranked_before, (out / "rank" / "ranked_candidates.tsv").read_bytes())
            rerun_evidence = (out / "evidence.tsv").read_bytes()
            self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources), "--out-dir", str(out), "--resume"]))
            self.assertEqual(rerun_evidence, (out / "evidence.tsv").read_bytes())

    def test_rank_hash_covers_min_top_score_and_rank_parameters(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case, resources, _, _ = self._evidence_case(root)
            out = root / "out"
            self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources), "--out-dir", str(out)]))
            state = json.loads((out / "case_state.json").read_text(encoding="utf-8"))
            rank_hash_before = state["steps"]["rank"]["input_hash"]
            ranked_before = (out / "rank" / "ranked_candidates.tsv").read_bytes()
            case.write_text(case.read_text(encoding="utf-8") + "min_top_score: 99.0\n", encoding="utf-8")
            self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources), "--out-dir", str(out), "--resume"]))
            state = json.loads((out / "case_state.json").read_text(encoding="utf-8"))
            self.assertNotEqual(rank_hash_before, state["steps"]["rank"]["input_hash"])
            ranked = pd.read_csv(out / "rank" / "ranked_candidates.tsv", sep="\t", dtype=str).fillna("")
            self.assertNotEqual("top", ranked.set_index("gene_id").loc["GeneA", "dry_group"])
            self.assertNotEqual(ranked_before, (out / "rank" / "ranked_candidates.tsv").read_bytes())
            case.write_text(case.read_text(encoding="utf-8").replace("min_top_score: 99.0\n", "").replace("  top_n: 2\n", "  top_n: 1\n"), encoding="utf-8")
            self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources), "--out-dir", str(out), "--resume"]))
            state = json.loads((out / "case_state.json").read_text(encoding="utf-8"))
            self.assertNotEqual(rank_hash_before, state["steps"]["rank"]["input_hash"])
            ranked = pd.read_csv(out / "rank" / "ranked_candidates.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(1, int((ranked["dry_group"] == "top").sum()))

    def test_strict_fingerprints_force_full_content_digest_on_resume(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates = root / "candidates.tsv"
            candidates.write_text("gene_id\n" + "\n".join(["GeneA"] + [f"Gene{i:05d}" for i in range(1, 2000)]) + "\n", encoding="utf-8")
            evidence = root / "evidence.tsv"
            evidence.write_text(
                "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\n"
                "GeneA\tknown_qtg\tsupport\thigh\tmodel\tliterature\tPMID:1\tgenetic model\n",
                encoding="utf-8",
            )
            case = root / "case.yaml"
            resources = root / "resources.yaml"
            out = root / "out"
            case.write_text(
                "case_id: strict_digest\n"
                "trait:\n  phenotype_direction: increase\n"
                "source:\n  source_id: PMID:QTL\n  target_interval: chr1:1-10\n  literature_search_status: complete\n"
                "inputs:\n  candidates_resource: candidates\n  evidence_resources: [evidence]\n",
                encoding="utf-8",
            )
            resources.write_text(f"resources:\n  candidates: {{path: {candidates}}}\n  evidence: {{path: {evidence}}}\n", encoding="utf-8")
            self.assertGreater(candidates.stat().st_size, 10 * 1024)
            calls: list[tuple[str, bool]] = []

            def tiered_fingerprint(path: Path, full_hash: bool = False) -> dict[str, object]:
                calls.append((str(path), full_hash))
                stat = path.stat()
                result: dict[str, object] = {"path": str(path), "size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}
                if full_hash or stat.st_size <= 10 * 1024:
                    result.update({"sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "method": "full_sha256", "strength": "strong"})
                else:
                    result.update({"sha256": None, "method": "metadata_only", "strength": "weak"})
                return result

            with patch.object(workflow_extensions, "file_fingerprint", side_effect=tiered_fingerprint):
                self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources), "--out-dir", str(out)]))
                evidence_output = (out / "evidence.tsv").read_bytes()
                stat_before = candidates.stat()
                candidates.write_text(candidates.read_text(encoding="utf-8").replace("GeneA", "GenAA"), encoding="utf-8")
                os.utime(candidates, ns=(stat_before.st_atime_ns, stat_before.st_mtime_ns))
                self.assertEqual(stat_before.st_size, candidates.stat().st_size)
                self.assertEqual(stat_before.st_mtime_ns, candidates.stat().st_mtime_ns)
                self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources), "--out-dir", str(out), "--resume"]))
                self.assertEqual(evidence_output, (out / "evidence.tsv").read_bytes())
                self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources), "--out-dir", str(out), "--resume", "--strict-fingerprints"]))
                self.assertNotEqual(evidence_output, (out / "evidence.tsv").read_bytes())
                self.assertNotIn("GeneA", (out / "evidence.tsv").read_text())
            self.assertTrue(any(path == str(candidates) and full_hash for path, full_hash in calls))


if __name__ == "__main__":
    unittest.main()
