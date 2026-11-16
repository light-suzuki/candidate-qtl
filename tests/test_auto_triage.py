from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl.auto_triage import _tool_version, auto_triage


class AutoTriageTests(unittest.TestCase):
    def test_missing_version_command_is_recorded_instead_of_failing(self) -> None:
        with patch("candidate_qtl.auto_triage.shutil.which", return_value=None), patch(
            "candidate_qtl.core.subprocess.run", side_effect=FileNotFoundError
        ):
            self.assertEqual("version_unavailable", _tool_version("auto"))

    def test_all_candidates_receive_mapping_or_failure_without_curated_gene_panel(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            candidates = root / "candidates.tsv"
            candidates.write_text("gene_id\tchrom\tstart\tend\nGeneA\tchr1\t1\t300\nGeneB\tchr1\t500\t800\n", encoding="utf-8")
            candidate_proteins = root / "candidates.faa"
            candidate_proteins.write_text(">GeneA\n" + "ACDEFGHIKLMNPQRSTVWY" * 6 + "\n>GeneB\n" + "YWVTSRQPNMLKIHGFEDCA" * 6 + "\n", encoding="utf-8")
            native = root / "native.faa"
            native.write_text(candidate_proteins.read_text(encoding="utf-8"), encoding="utf-8")
            reference = root / "reference.faa"
            reference.write_text(">RefA dwarf protein\n" + "ACDEFGHIKLMNPQRSTVWY" * 6 + "\n>RefB other protein\n" + "YWVTSRQPNMLKIHGFEDCA" * 6 + "\n", encoding="utf-8")
            annotations = root / "annotations.tsv"
            annotations.write_text(
                "reference_gene_id\tgene_description\tmutant_phenotype\tsource_id\n"
                "RefA\tgibberellin enzyme\tdwarf stem\tPMID:1\n"
                "RefB\tmembrane protein\t\tPMID:2\n",
                encoding="utf-8",
            )
            panel = root / "panel.yaml"
            panel.write_text(
                "panel_id: test\n"
                "native:\n  species: Pisum sativum\n  proteome: native.faa\n"
                "references:\n"
                "  - label: model\n    role: functional_reference\n    species: Arabidopsis thaliana\n    proteome: reference.faa\n    annotation: annotations.tsv\n",
                encoding="utf-8",
            )

            calls = 0

            def fake_blast(query_path, _subject_path, output, _runner, _threads):
                nonlocal calls
                calls += 1
                output.parent.mkdir(parents=True, exist_ok=True)
                if calls == 1:
                    output.write_text(
                        "GeneA\tRefA\t100\t120\t120\t120\t1\t120\t1\t120\t1e-80\t300\n"
                        "GeneB\tRefB\t100\t120\t120\t120\t1\t120\t1\t120\t1e-80\t300\n",
                        encoding="utf-8",
                    )
                else:
                    output.write_text(
                        "RefA\tGeneA\t100\t120\t120\t120\t1\t120\t1\t120\t1e-80\t300\n"
                        "RefB\tGeneB\t100\t120\t120\t120\t1\t120\t1\t120\t1e-80\t300\n",
                        encoding="utf-8",
                    )
                return {"output": str(output), "runner": "mock", "command": ["blastp"]}

            with patch("candidate_qtl.auto_triage._run_interval_blastp", side_effect=fake_blast):
                result = auto_triage(
                    "dwarf stem", str(panel), str(root / "out"),
                    str(candidates), str(candidate_proteins), target_size=2,
                )
            orthology = pd.read_csv(result["orthologs"], sep="\t", dtype=str).fillna("")
            self.assertEqual({"GeneA", "GeneB"}, set(orthology["candidate_gene_id"]))
            self.assertTrue(orthology["mapping_status"].eq("mapped").all())
            self.assertTrue(orthology["relationship"].eq("one_to_one_candidate").all())
            functions = pd.read_csv(result["functions"], sep="\t", dtype=str).fillna("")
            self.assertEqual({"source_retrieved"}, set(functions["retrieval_status"]))
            ranked = pd.read_csv(result["ranked"], sep="\t", dtype=str).fillna("")
            self.assertEqual("GeneA", ranked.iloc[0]["gene_id"])
            self.assertEqual("unreviewed", ranked.iloc[0]["review_status"])
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual(2, manifest["candidate_with_any_mapping_count"])
            self.assertFalse(manifest["causal_probability"])
            self.assertEqual("not_supplied", manifest["parental_variant_status"])
            for key in (
                "entity_aliases", "provider_runs", "candidate_provider_coverage",
                "normalized_observations", "evidence_lineage", "resource_manifest", "workflow_state",
            ):
                self.assertTrue(Path(result[key]).exists(), key)
            providers = pd.read_csv(result["provider_runs"], sep="\t", dtype=str).fillna("")
            self.assertEqual({"blastp:model", "annotation:model"}, set(providers["provider_name"]))

            with patch("candidate_qtl.auto_triage._run_interval_blastp", side_effect=fake_blast):
                auto_triage(
                    "dwarf stem", str(panel), str(root / "out"),
                    str(candidates), str(candidate_proteins), target_size=2, resume=True,
                )
            self.assertEqual(2, calls, "same input must reuse both forward and reciprocal provider outputs")

            reference.write_text(reference.read_text(encoding="utf-8").replace("dwarf", "DWARF"), encoding="utf-8")
            calls = 0
            with patch("candidate_qtl.auto_triage._run_interval_blastp", side_effect=fake_blast):
                auto_triage(
                    "dwarf stem", str(panel), str(root / "out"),
                    str(candidates), str(candidate_proteins), target_size=2, resume=True,
                )
            self.assertEqual(2, calls, "reference content change must invalidate forward and reciprocal cache")

            vcf = root / "parents.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "##contig=<ID=chr1,length=1000>\n"
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tparent1\tparent2\n"
                "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t1/1\n",
                encoding="utf-8",
            )
            calls = 0
            with patch("candidate_qtl.auto_triage._run_interval_blastp", side_effect=fake_blast):
                with_vcf = auto_triage(
                    "dwarf stem", str(panel), str(root / "out_with_vcf"),
                    str(candidates), str(candidate_proteins), parent_vcf=str(vcf), target_size=2,
                )
            with_vcf_manifest = json.loads(Path(with_vcf["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual("completed", with_vcf_manifest["parental_variant_status"])
            variant_rows = pd.read_csv(root / "out_with_vcf" / "05_parental_variants.tsv", sep="\t")
            self.assertEqual("GeneA", variant_rows.iloc[0]["gene_id"])
            self.assertTrue(bool(variant_rows.iloc[0]["genotype_differs"]))


if __name__ == "__main__":
    unittest.main()
