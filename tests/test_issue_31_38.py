from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl.cli import main


class Issue3138Tests(unittest.TestCase):
    def test_bundled_case_runs_to_ranked_primer_report(self) -> None:
        fixture = Path(__file__).parent / "fixtures" / "case_run"
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "out"
            self.assertEqual(0, main(["case-run", "--case", str(fixture / "case.yaml"), "--resources", str(fixture / "resources.yaml"), "--out-dir", str(out)]))
            ranked = pd.read_csv(out / "rank" / "ranked_candidates.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("top", ranked.set_index("gene_id").loc["GeneStrong", "dry_group"])
            primers = pd.read_csv(out / "primer" / "primer_candidates.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual("in_silico_unique", primers.iloc[0]["pair_specificity"])
            self.assertEqual("target_covered", primers.iloc[0]["coverage_status"])
            self.assertEqual("not_tested", primers.iloc[0]["wet_validation_status"])
            self.assertIn("Full source-level trace", (out / "report.md").read_text())
            history_path = out / "hypotheses" / "hypothesis_history.tsv"
            first_history = history_path.read_bytes()
            copied_case = Path(tmp) / "case.yaml"
            copied_case.write_text((fixture / "case.yaml").read_text(encoding="utf-8"), encoding="utf-8")
            self.assertEqual(0, main(["case-run", "--case", str(copied_case), "--resources", str(fixture / "resources.yaml"), "--out-dir", str(out)]))
            self.assertEqual(first_history, history_path.read_bytes())

    def test_case_run_uses_source_backed_evidence_not_gene_order(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); candidates = root / "candidates.tsv"; evidence = root / "evidence.tsv"
            candidates.write_text("gene_id\nA_first\nZ_supported\n", encoding="utf-8")
            evidence.write_text(
                "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\tphenotype_direction\n"
                "Z_supported\tknown_qtg\tsupport\thigh\tmodel\tliterature\tPMID:1\tgenetic model\tincrease\n"
                "A_first\tknown_qtg\tsupport\thigh\tmodel\tliterature\t\tmissing source\tincrease\n",
                encoding="utf-8",
            )
            llm = root / "llm.tsv"
            llm.write_text("gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\treviewed\nA_first\tknown_qtg\tsupport\thigh\tcurated\tliterature\tPMID:2\treviewed claim\ttrue\n", encoding="utf-8")
            case = root / "case.yaml"; resources = root / "resources.yaml"; out = root / "out"
            case.write_text("case_id: ranked\ntrait:\n  phenotype_direction: increase\nsource:\n  source_id: PMID:QTL\n  target_interval: chr1:1-10\n  literature_search_status: complete\ninputs:\n  candidates_resource: candidates\n  evidence_resources: [evidence]\n  llm_evidence_resources: [llm]\nassay:\n  top_n: 1\n", encoding="utf-8")
            resources.write_text(f"resources:\n  candidates:\n    path: {candidates}\n  evidence:\n    path: {evidence}\n  llm:\n    path: {llm}\n", encoding="utf-8")
            self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources), "--out-dir", str(out)]))
            ranked = list(csv.DictReader((out / "rank" / "ranked_candidates.tsv").open(), delimiter="\t"))
            self.assertEqual("Z_supported", ranked[0]["gene_id"]); self.assertEqual("top", ranked[0]["dry_group"])
            self.assertIn("A_first", (out / "evidence_excluded.tsv").read_text())
            hashes = [(path.name, path.read_bytes()) for path in [out / "evidence.tsv", out / "rank" / "ranked_candidates.tsv", out / "report.md"]]
            self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources), "--out-dir", str(out), "--resume"])); self.assertEqual(hashes, [(path.name, path.read_bytes()) for path in [out / "evidence.tsv", out / "rank" / "ranked_candidates.tsv", out / "report.md"]])
            deterministic_rank = (out / "rank" / "ranked_candidates.tsv").read_bytes()
            self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources), "--out-dir", str(out), "--resume", "--mode", "llm-assisted"])); self.assertEqual(deterministic_rank, (out / "rank" / "ranked_candidates.tsv").read_bytes()); self.assertTrue((out / "llm_assisted" / "rank" / "ranked_candidates.tsv").exists())

    def test_marker_pairs_split_direction_dark_and_case_resume(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            hits = root / "hits.tsv"
            hits.write_text(
                "marker\tprimer\tchrom\tstart\tend\tstrand\tidentity\n"
                "M1\tF\tchr1\t100\t120\t+\t100\nM1\tR\tchr1\t180\t200\t-\t100\n"
                "M2\tF\tchr1\t100000\t100020\t+\t100\nM2\tR\tchr1\t100080\t100100\t-\t100\n",
                encoding="utf-8",
            )
            out = root / "markers"
            self.assertEqual(0, main(["liftover-markers", "--hits", str(hits), "--out-dir", str(out), "--max-block-gap", "10000"]))
            result = json.loads((out / "marker_liftover.json").read_text())
            self.assertEqual("split", result["synteny_state"]); self.assertFalse(result["single_span_allowed"])

            evidence = root / "evidence.tsv"
            evidence.write_text("gene_id\tdirection\tphenotype_direction\nA\tsupport\tdecrease\nB\tsupport\tincrease\n", encoding="utf-8")
            directed = root / "directed.tsv"
            self.assertEqual(0, main(["direction-evidence", "--evidence", str(evidence), "--expected", "increase", "--out", str(directed)]))
            rows = list(csv.DictReader(directed.open(), delimiter="\t")); self.assertEqual("against", rows[0]["direction"])

            candidates = root / "candidates.tsv"
            candidates.write_text("gene_id\tannotation_quality\tprotein_length\thomology_hit\nD\tcomplete_orf\t90\t\nN\tcomplete_orf\t500\t\n", encoding="utf-8")
            dark = root / "dark.tsv"; self.assertEqual(0, main(["dark-candidates", "--candidates", str(candidates), "--out", str(dark)]))
            self.assertIn("D", dark.read_text()); self.assertNotIn("\nN\t", dark.read_text())

            case = root / "case.yaml"; resources = root / "resources.yaml"; run = root / "run"
            case.write_text("case_id: unresolved\nsource:\n  literature_search_status: complete\n", encoding="utf-8")
            resources.write_text("resources: {}\n", encoding="utf-8")
            self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources), "--out-dir", str(run)]))
            state1 = (run / "case_state.json").read_bytes()
            self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources), "--out-dir", str(run), "--resume"])); self.assertEqual(state1, (run / "case_state.json").read_bytes())
            self.assertIn("rescue", (run / "case_state.json").read_text())

    def test_qtlift_local_library_token_and_clear_missing_ref(self) -> None:
        import importlib.util
        script = Path(__file__).parents[1] / "scripts" / "qtlift_client.py"
        spec = importlib.util.spec_from_file_location("qtlift_client", script); module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
        calls = []
        def fake(base, method, path, payload=None):
            calls.append((path, payload)); return {"genomes": [{"name": "v1"}, {"name": "v2"}]} if path == "/api/genomes/scan" else {"job_id": "x"}
        with patch.object(module, "request_json", side_effect=fake), patch.object(module, "write_json"):
            self.assertEqual(0, module.main(["submit", "--local-library", "--source-ref", "v1", "--target-ref", "v2", "--contig", "c", "--start", "1", "--end", "2"]))
        self.assertEqual("@this-pc", calls[-1][1]["genome_root"])


if __name__ == "__main__": unittest.main()
