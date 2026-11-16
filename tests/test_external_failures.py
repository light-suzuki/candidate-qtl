from __future__ import annotations

import io
import json
import os
import subprocess
import tempfile
import time
import unittest
import urllib.error
from contextlib import redirect_stderr
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import candidate_qtl.discover as discover
from candidate_qtl.cli import EXIT_EXPECTED_FAILURE, run
from candidate_qtl.core import ExternalDependencyError, europe_pmc_search, run_tool


def _fixture(root: Path) -> tuple[Path, Path]:
    fasta = root / "genome.fa"
    gff = root / "genes.gff3"
    fasta.write_text(">chr1\n" + "ATG" * 200 + "\n", encoding="utf-8")
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tfixture\tgene\t100\t300\t.\t+\t.\tID=MeristemGene;Note=meristem regulator\n"
        "chr1\tfixture\tmRNA\t100\t300\t.\t+\t.\tID=MeristemGene.1;Parent=MeristemGene\n"
        "chr1\tfixture\tCDS\t100\t300\t.\t+\t0\tParent=MeristemGene.1\n",
        encoding="utf-8",
    )
    return fasta, gff


def _args(fasta: Path, gff: Path, out: Path) -> list[str]:
    return ["discover", "--fasta", str(fasta), "--gff3", str(gff), "--interval", "chr1:91-350",
            "--phenotype", "meristem enlargement", "--out-dir", str(out), "--no-auto-local-db",
            "--allow-external-api"]


class ExternalFailureTests(unittest.TestCase):
    def test_qblast_network_error_exits_one_with_kind_network_and_failure_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta, gff = _fixture(root)
            out = root / "out"
            with patch("Bio.Blast.NCBIWWW.qblast", side_effect=urllib.error.URLError("no route to host")) as qblast_mock, \
                    patch.object(discover, "_NCBI_RETRY_BACKOFF", 0.0), redirect_stderr(io.StringIO()) as err:
                code = run(_args(fasta, gff, out))
            self.assertEqual(EXIT_EXPECTED_FAILURE, code)
            self.assertIn("external dependency failure (kind=network)", err.getvalue())
            self.assertNotIn("defect in candidate-qtl", err.getvalue())
            self.assertEqual(discover._NCBI_ATTEMPTS, qblast_mock.call_count)
            manifest = json.loads((out / "discover_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("ncbi_blast_unavailable", manifest["annotation_provider"])
            self.assertEqual("network", manifest["external_api_failure"]["kind"])
            self.assertEqual("annotation_failed", manifest["workflow_stage"])
            self.assertFalse(manifest["external_results_used"])

    def test_qblast_hang_is_classified_as_timeout(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta, gff = _fixture(root)
            out = root / "out"

            def hanging(*_args, **_kwargs) -> None:
                time.sleep(2.0)
                raise AssertionError("hanging qblast must not finish")

            with patch.dict(os.environ, {"CANDIDATE_QTL_NCBI_TIMEOUT": "0.05"}), \
                    patch.object(discover, "_NCBI_RETRY_BACKOFF", 0.0), \
                    patch("Bio.Blast.NCBIWWW.qblast", side_effect=hanging), redirect_stderr(io.StringIO()) as err:
                code = run(_args(fasta, gff, out))
            self.assertEqual(EXIT_EXPECTED_FAILURE, code)
            self.assertIn("external dependency failure (kind=timeout)", err.getvalue())
            self.assertNotIn("defect in candidate-qtl", err.getvalue())

    def test_external_tool_timeout_exits_one_with_kind_timeout(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta, gff = _fixture(root)
            out = root / "out"
            args = _args(fasta, gff, out) + ["--protein-db", str(root / "absent_db")]
            with patch("candidate_qtl.core.subprocess.run", side_effect=subprocess.TimeoutExpired(["blastp"], 300)), \
                    redirect_stderr(io.StringIO()) as err:
                code = run(args)
            self.assertEqual(EXIT_EXPECTED_FAILURE, code)
            self.assertIn("external dependency failure (kind=timeout)", err.getvalue())
            self.assertNotIn("defect in candidate-qtl", err.getvalue())

    def test_missing_external_tool_keeps_file_not_found_classification(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fasta, gff = _fixture(root)
            out = root / "out"
            args = _args(fasta, gff, out) + ["--protein-db", str(root / "absent_db")]
            with patch("candidate_qtl.core.subprocess.run", side_effect=FileNotFoundError(2, "simulated missing blastp")), \
                    redirect_stderr(io.StringIO()) as err:
                code = run(args)
            self.assertEqual(EXIT_EXPECTED_FAILURE, code)
            self.assertIn("simulated missing blastp", err.getvalue())
            self.assertNotIn("external dependency failure", err.getvalue())
            self.assertNotIn("defect in candidate-qtl", err.getvalue())

    def test_run_tool_defaults_to_300s_and_honors_env_override(self) -> None:
        with patch("candidate_qtl.core.subprocess.run", return_value=SimpleNamespace(returncode=0)) as mock_run:
            run_tool(["blastp"])
            self.assertEqual(300, mock_run.call_args.kwargs["timeout"])
        with patch.dict(os.environ, {"CANDIDATE_QTL_TOOL_TIMEOUT": "0.25"}), \
                patch("candidate_qtl.core.subprocess.run", return_value=SimpleNamespace(returncode=0)) as mock_run:
            run_tool(["blastp"])
            self.assertEqual(0.25, mock_run.call_args.kwargs["timeout"])

    def test_run_tool_wraps_timeout_expired_with_kind_timeout(self) -> None:
        with patch("candidate_qtl.core.subprocess.run", side_effect=subprocess.TimeoutExpired(["sleep"], 7)) as mock_run:
            with self.assertRaises(ExternalDependencyError) as raised:
                run_tool(["sleep"], timeout=7)
        self.assertEqual("timeout", raised.exception.kind)
        self.assertIn("external dependency failure (kind=timeout)", str(raised.exception))
        self.assertIn("sleep", str(raised.exception))
        mock_run.assert_called_once_with(["sleep"], timeout=7)

    def test_europe_pmc_network_failure_is_classified_network(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "literature.json"
            with patch("candidate_qtl.core.urlopen", side_effect=urllib.error.URLError("dns failure")):
                with self.assertRaises(ExternalDependencyError) as raised:
                    europe_pmc_search("query terms", str(out), 10)
            self.assertEqual("network", raised.exception.kind)
            self.assertIn("external dependency failure (kind=network)", str(raised.exception))


if __name__ == "__main__":
    unittest.main()
