from __future__ import annotations

import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import yaml

from candidate_qtl.cli import main


class DiscoverTests(unittest.TestCase):
    def test_source_backed_direction_and_inheritance_affect_ranking_and_promote(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fasta = root / "genome.fa"; gff = root / "genes.gff3"; supplied = root / "evidence.tsv"
            fasta.write_text(">chr1\n" + "ATG" * 300 + "\n", encoding="utf-8")
            gff.write_text(
                "##gff-version 3\n"
                "chr1\tx\tgene\t100\t250\t.\t+\t.\tID=MatchGene\nchr1\tx\tmRNA\t100\t250\t.\t+\t.\tID=MatchGene.1;Parent=MatchGene\nchr1\tx\tCDS\t100\t249\t.\t+\t0\tParent=MatchGene.1\n"
                "chr1\tx\tgene\t300\t450\t.\t+\t.\tID=OppositeGene\nchr1\tx\tmRNA\t300\t450\t.\t+\t.\tID=OppositeGene.1;Parent=OppositeGene\nchr1\tx\tCDS\t300\t449\t.\t+\t0\tParent=OppositeGene.1\n",
                encoding="utf-8",
            )
            supplied.write_text(
                "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\tphenotype_direction\tinheritance_model\tmolecular_effect\n"
                "MatchGene\tdirect_model_genetic\tsupport\thigh\tmutant\tliterature\tPMID:1\tmatching mutant\tincrease\trecessive\tloss_of_function\n"
                "OppositeGene\tdirect_model_genetic\tsupport\thigh\tmutant\tliterature\tPMID:2\topposite mutant\tdecrease\tdominant\tgain_of_function\n"
                "MatchGene\ttrait_function_match\tsupport\thigh\tclaim\tmanual\t\tmissing provenance\tincrease\trecessive\tloss_of_function\n",
                encoding="utf-8",
            )
            out = root / "discover"
            self.assertEqual(0, main(["discover", "--fasta", str(fasta), "--gff3", str(gff), "--interval", "chr1:1-500",
                                      "--phenotype", "architecture", "--direction", "increase", "--inheritance", "recessive",
                                      "--evidence", str(supplied), "--out-dir", str(out), "--offline", "--no-auto-local-db"]))
            evidence = pd.read_csv(out / "04_evidence.tsv", sep="\t", dtype=str).fillna("")
            by_gene = evidence[evidence["evidence_code"].eq("direct_model_genetic")].set_index("gene_id")
            self.assertEqual("support", by_gene.loc["MatchGene", "direction"])
            self.assertEqual("against", by_gene.loc["OppositeGene", "direction"])
            self.assertEqual("match", by_gene.loc["MatchGene", "inheritance_match"])
            self.assertEqual("mismatch", by_gene.loc["OppositeGene", "inheritance_match"])
            self.assertIn("missing provenance", (out / "evidence_excluded.tsv").read_text(encoding="utf-8"))
            manifest = json.loads((out / "discover_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(manifest["evidence_constraints"]["direction"]["used_for_ranking"])
            self.assertTrue(manifest["evidence_constraints"]["inheritance"]["used_for_ranking"])
            self.assertEqual("preliminary_screening", manifest["workflow_stage"])
            self.assertIn("case_run_source_gate", manifest["unrun_stages"])
            promoted = root / "promoted"
            self.assertEqual(0, main(["discover-promote", "--discover-dir", str(out), "--out-dir", str(promoted), "--case-id", "case1"]))
            self.assertIn("REVIEW_REQUIRED", (promoted / "case.yaml").read_text(encoding="utf-8"))
            self.assertTrue((promoted / "resources.yaml").exists())
            case = yaml.safe_load((promoted / "case.yaml").read_text(encoding="utf-8"))
            case["source"]["source_id"] = "PMID:QTL"; case["source"]["literature_search_status"] = "complete"
            case["target"]["assembly"] = "fixture_v1"
            (promoted / "case.yaml").write_text(yaml.safe_dump(case, sort_keys=False), encoding="utf-8")
            case_run = root / "case_run"
            self.assertEqual(0, main(["case-run", "--case", str(promoted / "case.yaml"), "--resources", str(promoted / "resources.yaml"), "--out-dir", str(case_run)]))
            self.assertTrue((case_run / "report.md").exists())

    def test_unknown_direction_and_unbacked_inheritance_are_explicitly_unused(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fasta = root / "genome.fa"; gff = root / "genes.gff3"; out = root / "out"
            fasta.write_text(">chr1\n" + "ATG" * 100 + "\n", encoding="utf-8")
            gff.write_text("##gff-version 3\nchr1\tx\tgene\t1\t90\t.\t+\t.\tID=G1\n", encoding="utf-8")
            self.assertEqual(0, main(["discover", "--fasta", str(fasta), "--gff3", str(gff), "--interval", "chr1:1-100",
                                      "--phenotype", "trait", "--inheritance", "recessive", "--out-dir", str(out), "--offline", "--no-auto-local-db"]))
            manifest = json.loads((out / "discover_manifest.json").read_text(encoding="utf-8"))
            self.assertFalse(manifest["evidence_constraints"]["direction"]["used_for_ranking"])
            self.assertFalse(manifest["evidence_constraints"]["inheritance"]["used_for_ranking"])
            self.assertIn("molecular_effect", manifest["evidence_constraints"]["inheritance"]["unused_reason"])

    def test_discover_uses_extract_safety_options_and_matches_extract(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fasta = root / "genome.fa"; gff = root / "genes.gff3"
            fasta.write_text(">NC_test\n" + "ATG" * 300 + "\n", encoding="utf-8")
            gff.write_text(
                "##gff-version 3\n"
                "chrAlias\tfixture\tgene\t100\t300\t.\t+\t.\tID=Gene1;Note=meristem regulator\n"
                "chrAlias\tfixture\tmRNA\t100\t300\t.\t+\t.\tID=Gene1.1;Parent=Gene1\n"
                "chrAlias\tfixture\tCDS\t100\t300\t.\t+\t0\tParent=Gene1.1\n",
                encoding="utf-8",
            )
            qc = root / "gff_qc.txt"; qc.write_text("validation passed\n", encoding="utf-8")
            provenance = root / "qtl_source_preflight.json"
            provenance.write_text(json.dumps({"positional_ranking_allowed": True, "position_status": "provisional_interval"}), encoding="utf-8")
            extract_out, discover_out = root / "extract", root / "discover"
            common = ["--fasta", str(fasta), "--gff3", str(gff), "--interval", "chrAlias:91-350",
                      "--seqid-map", "chrAlias=NC_test", "--gff-qc-report", str(qc),
                      "--qtl-source-manifest", str(provenance), "--full-hash"]
            self.assertEqual(0, main(["extract", *common, "--flank", "500", "--out-dir", str(extract_out)]))
            self.assertEqual(0, main(["discover", *common, "--phenotype", "meristem enlargement",
                                      "--out-dir", str(discover_out), "--offline", "--no-auto-local-db"] ))
            expected = pd.read_csv(extract_out / "candidates.tsv", sep="\t", dtype=str).fillna("")
            actual = pd.read_csv(discover_out / "01_extract" / "candidates.tsv", sep="\t", dtype=str).fillna("")
            pd.testing.assert_frame_equal(expected, actual)
            manifest = json.loads((discover_out / "discover_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("chrAlias=NC_test", manifest["extraction"]["seqid_map"])
            self.assertTrue(manifest["extraction"]["full_hash"])
            extract_manifest = json.loads((discover_out / "01_extract" / "block_01" / "extract_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(extract_manifest["seqid_mapping_applied"])
            self.assertEqual(["validation passed"], extract_manifest["external_gff_qc"]["excerpt"])

    def test_discover_rejects_unresolved_and_requires_split_intervals(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fasta = root / "genome.fa"; gff = root / "genes.gff3"
            fasta.write_text(">chr1\n" + "ATG" * 300 + "\n", encoding="utf-8")
            gff.write_text("##gff-version 3\nchr1\tx\tgene\t100\t200\t.\t+\t.\tID=Gene1\n", encoding="utf-8")
            base = ["discover", "--fasta", str(fasta), "--gff3", str(gff), "--interval", "chr1:91-250",
                    "--phenotype", "trait", "--offline", "--no-auto-local-db"]
            unresolved = root / "unresolved.json"
            unresolved.write_text(json.dumps({"positional_ranking_allowed": False, "position_status": "unresolved"}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "coordinate provenance gate rejected"):
                main([*base, "--qtl-source-manifest", str(unresolved), "--out-dir", str(root / "rejected")])
            split = root / "split.json"
            split.write_text(json.dumps({"positional_ranking_allowed": True, "position_status": "provisional_interval", "synteny_state": "split"}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "requires explicit multiple"):
                main([*base, "--qtl-source-manifest", str(split), "--out-dir", str(root / "split_rejected")])

    def test_minimal_inputs_produce_source_traced_candidates(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fasta = root / "genome.fa"; gff = root / "genes.gff3"; out = root / "out"
            fasta.write_text(">chr1\n" + "ATG" * 200 + "\n", encoding="utf-8")
            gff.write_text(
                "##gff-version 3\n"
                "chr1\tfixture\tgene\t1\t ninety\t.\t+\t.\tID=bad\n".replace(" ninety", "90")
                + "chr1\tfixture\tgene\t100\t300\t.\t+\t.\tID=MeristemGene;Note=meristem regulator\n"
                + "chr1\tfixture\tmRNA\t100\t300\t.\t+\t.\tID=MeristemGene.1;Parent=MeristemGene\n"
                + "chr1\tfixture\tCDS\t100\t300\t.\t+\t0\tParent=MeristemGene.1\n",
                encoding="utf-8",
            )
            self.assertEqual(0, main(["discover", "--fasta", str(fasta), "--gff3", str(gff), "--interval", "chr1:91-350",
                                      "--phenotype", "meristem enlargement", "--direction", "increase", "--inheritance", "recessive",
                                      "--out-dir", str(out), "--top-n", "3", "--offline", "--no-auto-local-db"]))
            ranked = pd.read_csv(out / "05_rank" / "ranked_candidates.tsv", sep="\t", dtype=str).fillna("")
            # Interval membership + a single keyword match is not enough for Top.
            self.assertEqual("evidence_insufficient", ranked.set_index("gene_id").loc["MeristemGene", "dry_group"])
            evidence = (out / "04_evidence.tsv").read_text(encoding="utf-8")
            self.assertIn("trait_function_match", evidence); self.assertIn(str(gff), evidence)
            # Interval membership is recorded as context, never scored support_interval.
            self.assertIn("interval_membership\tcontext", evidence)
            manifest = json.loads((out / "discover_manifest.json").read_text())
            self.assertFalse(manifest["external_sequence_submission_this_run"]); self.assertEqual("gff3_only", manifest["annotation_provider"])
            self.assertIn("protein_annotation_not_configured", (out / "review_required.tsv").read_text())

            external_out = root / "external"
            fake_hits = {"MeristemGene": [{"subject":"sp|TEST","identity":80.0,"query_coverage":1.0,"subject_coverage":0.9,"evalue":1e-20,"bitscore":100.0,"description":"meristem regulator"}]}
            with patch("candidate_qtl.discover._best_ncbi_hits", return_value=(fake_hits, "", True, False)):
                self.assertEqual(0, main(["discover", "--fasta", str(fasta), "--gff3", str(gff), "--interval", "chr1:91-350",
                                          "--phenotype", "meristem enlargement", "--out-dir", str(external_out), "--no-auto-local-db", "--allow-external-api"]))
            external_manifest = json.loads((external_out / "discover_manifest.json").read_text())
            self.assertTrue(external_manifest["external_sequence_submission_this_run"])
            self.assertTrue(external_manifest["external_results_used"])


    def test_external_transmission_is_fail_closed_with_warning_and_dry_run(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fasta = root / "genome.fa"; gff = root / "genes.gff3"
            fasta.write_text(">chr1\n" + "ATG" * 200 + "\n", encoding="utf-8")
            gff.write_text(
                "##gff-version 3\n"
                "chr1\tfixture\tgene\t100\t300\t.\t+\t.\tID=MeristemGene;Note=meristem regulator\n"
                "chr1\tfixture\tmRNA\t100\t300\t.\t+\t.\tID=MeristemGene.1;Parent=MeristemGene\n"
                "chr1\tfixture\tCDS\t100\t300\t.\t+\t0\tParent=MeristemGene.1\n",
                encoding="utf-8",
            )

            def discover_args(out_dir: Path, *extra: str) -> list[str]:
                return ["discover", "--fasta", str(fasta), "--gff3", str(gff), "--interval", "chr1:91-350",
                        "--phenotype", "meristem enlargement", "--out-dir", str(out_dir), "--no-auto-local-db", *extra]

            # Fail-closed: the NCBI fallback never runs without --allow-external-api.
            with patch("candidate_qtl.discover._best_ncbi_hits", return_value=({}, "", True, False)) as ncbi:
                default_out = root / "default"
                self.assertEqual(0, main(discover_args(default_out)))
                ncbi.assert_not_called()
            default_manifest = json.loads((default_out / "discover_manifest.json").read_text(encoding="utf-8"))
            self.assertFalse(default_manifest["external_api_allowed"])
            self.assertFalse(default_manifest["external_transmission_warning"])
            self.assertFalse(default_manifest["external_sequence_submission_this_run"])
            self.assertEqual("gff3_only", default_manifest["annotation_provider"])

            # Explicit opt-in transmits and warns on stderr and in the manifest.
            with patch("candidate_qtl.discover._best_ncbi_hits", return_value=({}, "", True, False)) as ncbi, \
                    patch("sys.stderr", new_callable=io.StringIO) as stderr:
                external_out = root / "external"
                self.assertEqual(0, main(discover_args(external_out, "--allow-external-api")))
                ncbi.assert_called_once()
            external_manifest = json.loads((external_out / "discover_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(external_manifest["external_transmission_warning"])
            self.assertTrue(external_manifest["external_sequence_submission_this_run"])
            self.assertIn("candidate protein sequences may be transmitted to the external NCBI BLAST API", stderr.getvalue())

            # Dry-run reports the plan without any transmission or search.
            with patch("candidate_qtl.discover._best_ncbi_hits", return_value=({}, "", True, False)) as ncbi, \
                    patch("sys.stderr", new_callable=io.StringIO) as stderr:
                dry_out = root / "dry"
                self.assertEqual(0, main(discover_args(dry_out, "--allow-external-api", "--dry-run")))
                ncbi.assert_not_called()
            dry_manifest = json.loads((dry_out / "discover_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(dry_manifest["dry_run"])
            self.assertTrue(dry_manifest["dry_run_plan"]["external_transmission_would_occur"])
            self.assertEqual("ncbi_blast_swissprot", dry_manifest["dry_run_plan"]["planned_provider"])
            self.assertFalse(dry_manifest["external_transmission_warning"])
            self.assertFalse(dry_manifest["external_sequence_submission_this_run"])
            self.assertIn("DRY RUN", stderr.getvalue())

    def test_top_requires_more_than_interval_and_keyword(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fasta = root / "genome.fa"; gff = root / "genes.gff3"
            fasta.write_text(">chr1\n" + "ATG" * 200 + "\n", encoding="utf-8")
            gff.write_text(
                "##gff-version 3\n"
                "chr1\tfixture\tgene\t100\t300\t.\t+\t.\tID=MeristemGene;Note=meristem regulator\n"
                "chr1\tfixture\tmRNA\t100\t300\t.\t+\t.\tID=MeristemGene.1;Parent=MeristemGene\n"
                "chr1\tfixture\tCDS\t100\t300\t.\t+\t0\tParent=MeristemGene.1\n"
                "chr1\tfixture\tgene\t100\t180\t.\t+\t.\tID=PlainGene;Note=conserved hypothetical protein\n"
                "chr1\tfixture\tmRNA\t100\t180\t.\t+\t.\tID=PlainGene.1;Parent=PlainGene\n"
                "chr1\tfixture\tCDS\t100\t180\t.\t+\t0\tParent=PlainGene.1\n",
                encoding="utf-8",
            )
            def args(out_dir: Path) -> list[str]:
                return ["discover", "--fasta", str(fasta), "--gff3", str(gff), "--interval", "chr1:91-350",
                        "--phenotype", "meristem enlargement", "--out-dir", str(out_dir), "--top-n", "5"]

            # Cases 1 & 2: interval-only and keyword+interval, both offline.
            offline = root / "offline"
            self.assertEqual(0, main(args(offline) + ["--offline", "--no-auto-local-db"]))
            ranked = pd.read_csv(offline / "05_rank" / "ranked_candidates.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            # Interval membership never creates a support family; keyword+interval is not Top.
            self.assertEqual(int(ranked.loc["PlainGene", "support_family_count"]), 0)   # interval only
            self.assertEqual(int(ranked.loc["MeristemGene", "support_family_count"]), 1)  # keyword + interval
            self.assertNotEqual(ranked.loc["MeristemGene", "dry_group"], "top")
            self.assertNotEqual(ranked.loc["PlainGene", "dry_group"], "top")

            # Case 3: keyword + homology ranks the candidate above an interval-only
            # gene, but homology and the keyword share the comparative_function
            # channel, so they do not fabricate a second independent family.
            external = root / "external"
            fake_hits = {"MeristemGene": [{"subject": "sp|REG", "identity": 80.0, "query_coverage": 1.0, "subject_coverage": 0.9, "evalue": 1e-25, "bitscore": 120.0, "description": "meristem regulator kinase"}]}
            with patch("candidate_qtl.discover._best_ncbi_hits", return_value=(fake_hits, "", True, False)):
                self.assertEqual(0, main(args(external) + ["--no-auto-local-db", "--allow-external-api"]))
            ranked_ext = pd.read_csv(external / "05_rank" / "ranked_candidates.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            # Generic homology is unified as context (homology_context), like evidence-build.
            self.assertIn("homology_context", (external / "04_evidence.tsv").read_text(encoding="utf-8"))
            self.assertEqual(int(ranked_ext.loc["MeristemGene", "screening_rank"]), 1)
            self.assertGreater(int(ranked_ext.loc["PlainGene", "screening_rank"]), 1)
            self.assertEqual(int(ranked_ext.loc["MeristemGene", "support_family_count"]), 1)


    def test_promote_persists_ranking_defaults_and_preserves_dry_groups(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fasta = root / "genome.fa"; gff = root / "genes.gff3"; supplied = root / "evidence.tsv"
            fasta.write_text(">chr1\n" + "ATG" * 300 + "\n", encoding="utf-8")
            gff.write_text(
                "##gff-version 3\n"
                "chr1\tx\tgene\t100\t250\t.\t+\t.\tID=StrongGene\nchr1\tx\tmRNA\t100\t250\t.\t+\t.\tID=StrongGene.1;Parent=StrongGene\nchr1\tx\tCDS\t100\t249\t.\t+\t0\tParent=StrongGene.1\n"
                "chr1\tx\tgene\t300\t450\t.\t+\t.\tID=WeakGene\nchr1\tx\tmRNA\t300\t450\t.\t+\t.\tID=WeakGene.1;Parent=WeakGene\nchr1\tx\tCDS\t300\t449\t.\t+\t0\tParent=WeakGene.1\n",
                encoding="utf-8",
            )
            # WeakGene scores 1.25 (trait_function_match 1.0 + expression_detected 0.25):
            # above the discover lane's lenient gate (1.0) but below the case-run
            # fallback gate (2.0), so forgetting to persist the defaults drops it.
            supplied.write_text(
                "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\n"
                "StrongGene\tdirect_model_genetic\tsupport\thigh\tmutant\tliterature\tPMID:1\tmatching mutant\n"
                "WeakGene\ttrait_function_match\tsupport\thigh\tclaim\tmanual\tPMID:2\tmeristem regulator\n"
                "WeakGene\texpression_detected\tsupport\tmedium\tdetected\texpression\tPMID:3\texpressed in meristem\n",
                encoding="utf-8",
            )
            out = root / "discover"
            self.assertEqual(0, main(["discover", "--fasta", str(fasta), "--gff3", str(gff), "--interval", "chr1:1-500",
                                      "--phenotype", "architecture", "--evidence", str(supplied), "--top-n", "1",
                                      "--out-dir", str(out), "--offline", "--no-auto-local-db"]))
            discover_groups = pd.read_csv(out / "05_rank" / "ranked_candidates.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")["dry_group"].to_dict()
            manifest = json.loads((out / "discover_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(1.0, manifest["ranking_defaults"]["min_top_score"])
            self.assertEqual({"top_n": 1, "reserve_n": 1, "rescue_n": 1, "min_top_score": 1.0},
                             {key: manifest["ranking_defaults"][key] for key in ("top_n", "reserve_n", "rescue_n", "min_top_score")})
            promoted = root / "promoted"
            self.assertEqual(0, main(["discover-promote", "--discover-dir", str(out), "--out-dir", str(promoted), "--case-id", "case1"]))
            case = yaml.safe_load((promoted / "case.yaml").read_text(encoding="utf-8"))
            self.assertEqual(1.0, case["min_top_score"])
            self.assertEqual({"top_n": 1, "reserve_n": 1, "rescue_n": 1}, case["assay"])
            case["source"]["source_id"] = "PMID:QTL"; case["source"]["literature_search_status"] = "complete"
            case["target"]["assembly"] = "fixture_v1"
            (promoted / "case.yaml").write_text(yaml.safe_dump(case, sort_keys=False), encoding="utf-8")
            case_run = root / "case_run"
            self.assertEqual(0, main(["case-run", "--case", str(promoted / "case.yaml"), "--resources", str(promoted / "resources.yaml"), "--out-dir", str(case_run)]))
            case_groups = pd.read_csv(case_run / "rank" / "ranked_candidates.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")["dry_group"].to_dict()
            self.assertEqual(discover_groups, case_groups)
            self.assertEqual("top", case_groups["StrongGene"])
            self.assertEqual("reserve", case_groups["WeakGene"])


if __name__ == "__main__": unittest.main()
