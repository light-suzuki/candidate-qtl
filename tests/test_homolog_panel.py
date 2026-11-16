from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from candidate_qtl.cli import build_parser
from candidate_qtl.homolog_panel import _known_gene_queries, _llm_hypotheses, _mechanical_queries, _overlap, _repair_hypothesis_json, _score_record, automated_homolog_match, automated_homolog_panel, benchmark_automated_homolog_panels


def _record(accession: str, species: str, comment_text: str, pmid: str | None, sequence: str = "M" * 120) -> dict[str, object]:
    evidences = [{"source": "PubMed", "id": pmid}] if pmid else []
    return {
        "primaryAccession": accession,
        "uniProtkbId": f"{accession}_ENTRY",
        "annotationScore": 5,
        "organism": {"scientificName": species, "lineage": ["Viridiplantae", "Fabaceae"]},
        "genes": [{"geneName": {"value": f"GENE_{accession}"}}],
        "proteinDescription": {"recommendedName": {"fullName": {"value": f"protein {accession}"}}},
        "comments": [{"commentType": "DISRUPTION PHENOTYPE", "texts": [{"value": comment_text, "evidences": evidences}]}],
        "entryAudit": {"entryVersion": 3, "lastAnnotationUpdateDate": "2026-01-01"},
        "sequence": {"value": sequence, "length": len(sequence)},
    }


class AutomatedHomologPanelTests(unittest.TestCase):
    def test_no_llm_query_expansion_is_deterministic(self) -> None:
        rows = _mechanical_queries("axial versus terminal or fasciated flowering")
        queries = [row["query"] for row in rows]
        self.assertIn("fasciated flowering", queries)
        self.assertIn("fasciated", queries)
        self.assertIn("fasciation", queries)
        self.assertEqual(["SBEII", "GBSSI"], _known_gene_queries("SBEII, GBSSI"))
        self.assertEqual([], _known_gene_queries("None"))
        self.assertFalse(_overlap("auxin transport regulation determining inflorescence architecture polarity", "general regulation of membrane transport"))
        self.assertTrue(_overlap("WUSCHEL CLAVATA signaling shoot apical meristem", "CLAVATA signaling maintains the shoot apical meristem"))
        repaired = _repair_hypothesis_json(
            '{\n"search_hypotheses": [{\n"query": "meristem identity",\n"mechanism": "CLV signaling",\n"known_gene": "WUS",\n"rationale":\nMissing opening quote but ends here."\n}]}\n'
        )
        self.assertEqual("Missing opening quote but ends here.", repaired["search_hypotheses"][0]["rationale"])
        repaired_closing = _repair_hypothesis_json(
            '{\n"search_hypotheses": [{\n"query": "meristem identity",\n'
            '"mechanism": "CLV signaling",\n"known_gene": "WUS",\n'
            '"rationale": "Missing closing quote.\n}]}\n'
        )
        self.assertEqual("Missing closing quote.", repaired_closing["search_hypotheses"][0]["rationale"])

    def test_free_phenotype_requires_explicit_uniprot_transmission(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            with self.assertRaisesRegex(ValueError, "allow-external-api"):
                automated_homolog_panel("wrinkled seed", "Pisum sativum", tempdir)

    def test_panel_is_sequence_backed_pending_review_and_auditable(self) -> None:
        payload = {
            "results": [
                _record("P00001", "Medicago truncatula", "Mutants have wrinkled seed morphology.", "12345"),
                _record("P00002", "Arabidopsis thaliana", "Mutants flower late.", None),
                _record("PISUM1", "Pisum sativum", "Mutants have wrinkled seed morphology.", "999"),
            ]
        }
        with tempfile.TemporaryDirectory() as tempdir, patch("candidate_qtl.homolog_panel._http_json", return_value=payload):
            result = automated_homolog_panel(
                "wrinkled seed", "Pisum sativum", tempdir,
                allow_external_api=True, panel_size=2, per_query=3,
            )
            panel = pd.read_csv(result["panel"], sep="\t", dtype=str).fillna("")
            self.assertEqual(["P00001", "P00002"], list(panel["reference_gene_id"]))
            self.assertEqual("direct_model_genetic", panel.iloc[0]["evidence_level"])
            self.assertEqual("database_context", panel.iloc[1]["evidence_level"])
            self.assertTrue(panel["review_status"].eq("pending_review").all())
            self.assertIn("PMID:12345", panel.iloc[0]["source_id"])
            self.assertIn(">P00001", Path(result["fasta"]).read_text(encoding="utf-8"))
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertFalse(manifest["external_transmission"]["qtl_coordinates"])
            self.assertFalse(manifest["external_transmission"]["candidate_ids"])
            self.assertFalse(manifest["llm"]["used"])

    def test_unrelated_free_text_does_not_inherit_an_entrys_genetic_comment(self) -> None:
        scored = _score_record(
            _record("P00003", "Arabidopsis thaliana", "Mutants flower late.", "55"),
            "wrinkled seed", "input",
        )
        self.assertEqual("database_context", scored["evidence_level"])
        self.assertNotIn("DISRUPTION PHENOTYPE", scored["channels"])

    def test_llm_receives_trait_metadata_and_only_proposes_queries(self) -> None:
        payload = {"results": [_record("P00004", "Arabidopsis thaliana", "Mutants have short stems.", "77")]}
        seen_prompt: dict[str, str] = {}

        def fake_llm(prompts_arg, out_dir_arg, *_args, **kwargs):
            self.assertFalse(kwargs["allow_external_api"])
            self.assertEqual(1800, kwargs["max_output_tokens"])
            prompt_path = Path(prompts_arg)
            seen_prompt.update(json.loads(prompt_path.read_text(encoding="utf-8")))
            out_dir = Path(out_dir_arg)
            out_dir.mkdir(parents=True, exist_ok=True)
            normalized = out_dir / "responses.json"
            normalized.write_text(json.dumps({"candidates": [{"model_json": {"search_hypotheses": [{"query": "short stem", "mechanism": "growth", "known_gene": "", "rationale": "alternative"}]}}]}), encoding="utf-8")
            manifest = out_dir / "manifest.json"
            manifest.write_text("{}", encoding="utf-8")
            return {"normalized_responses": str(normalized), "manifest": str(manifest)}

        with tempfile.TemporaryDirectory() as tempdir, \
                patch("candidate_qtl.homolog_panel.run_llm_prompts", side_effect=fake_llm), \
                patch("candidate_qtl.homolog_panel._http_json", return_value=payload):
            result = automated_homolog_panel(
                "dwarf plant", "Pisum sativum", tempdir,
                llm_provider="openai-compatible", model="local-test",
                endpoint="http://127.0.0.1:1234/v1/chat/completions",
                allow_external_api=True, panel_size=1,
            )
            self.assertIn("dwarf plant", seen_prompt["user"])
            self.assertNotIn("candidate", seen_prompt["user"].casefold())
            hypotheses = pd.read_csv(Path(tempdir) / "search_hypotheses.tsv", sep="\t")
            self.assertIn("short stem", set(hypotheses["query"]))
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertTrue(manifest["llm"]["used"])

    def test_llm_resume_requires_matching_prompt_provider_and_model(self) -> None:
        calls: list[str] = []

        def fake_llm(_prompts, out_dir_arg, _provider, model, *_args, **_kwargs):
            calls.append(model)
            out_dir = Path(out_dir_arg)
            out_dir.mkdir(parents=True, exist_ok=True)
            normalized = out_dir / "llm_candidate_responses.json"
            normalized.write_text(
                json.dumps({"candidates": [{"model_json": {"search_hypotheses": [{
                    "query": "short stem", "mechanism": "growth", "known_gene": "",
                    "rationale": "alternative",
                }]}}]}),
                encoding="utf-8",
            )
            manifest = out_dir / "llm_run_manifest.json"
            manifest.write_text(json.dumps({"external_api": False}), encoding="utf-8")
            return {"normalized_responses": str(normalized), "manifest": str(manifest)}

        with tempfile.TemporaryDirectory() as tempdir, patch(
            "candidate_qtl.homolog_panel.run_llm_prompts", side_effect=fake_llm
        ):
            root = Path(tempdir)
            first = _llm_hypotheses(
                "dwarf plant", "Pisum sativum", "unknown", "unknown", root,
                "openai-compatible", "model-a", None, False, True,
            )
            repeated = _llm_hypotheses(
                "dwarf plant", "Pisum sativum", "unknown", "unknown", root,
                "openai-compatible", "model-a", None, False, True,
            )
            changed = _llm_hypotheses(
                "dwarf plant", "Pisum sativum", "unknown", "unknown", root,
                "openai-compatible", "model-b", None, False, True,
            )

        self.assertEqual(["model-a", "model-b"], calls)
        self.assertFalse(first[1]["reused"])
        self.assertTrue(repeated[1]["reused"])
        self.assertFalse(changed[1]["reused"])
        self.assertNotEqual(first[1]["request_fingerprint"], changed[1]["request_fingerprint"])

    def test_interval_match_caps_database_context_at_low_confidence(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            candidates = root / "candidates.tsv"
            candidates.write_text("gene_id\nGeneA\nGeneB\n", encoding="utf-8")
            proteins = root / "proteins.faa"
            proteins.write_text(">GeneA\n" + "M" * 120 + "\n>GeneB\n" + "A" * 120 + "\n", encoding="utf-8")
            panel = root / "panel.tsv"
            panel.write_text(
                "reference_gene_id\tsource_id\tevidence_level\n"
                "P00001\tUniProt:P00001\tdatabase_context\n",
                encoding="utf-8",
            )
            panel_fasta = root / "panel.faa"
            panel_fasta.write_text(">P00001\n" + "M" * 120 + "\n", encoding="utf-8")
            blast = root / "blast.tsv"
            blast.write_text("P00001\tGeneA\t80\t120\t120\t120\t1\t120\t1\t120\t1e-50\t250\n", encoding="utf-8")
            with patch("candidate_qtl.homolog_panel._run_interval_blastp", return_value={"output": str(blast), "runner": "mock"}):
                result = automated_homolog_match(
                    str(candidates), str(proteins), str(panel), str(panel_fasta), str(root / "out")
                )
            evidence = pd.read_csv(result["evidence"], sep="\t", dtype=str).fillna("")
            self.assertEqual("automated_homolog_match", evidence.iloc[0]["evidence_code"])
            self.assertEqual("low", evidence.iloc[0]["confidence"])
            self.assertEqual("UniProt:P00001", evidence.iloc[0]["source_id"])

    def test_cli_exposes_provider_neutral_panel_route(self) -> None:
        args = build_parser().parse_args([
            "homolog-panel-auto", "--phenotype", "dwarf", "--species", "Pisum sativum",
            "--out-dir", "out", "--allow-external-api", "--llm-provider", "openai-responses", "--model", "gpt-test",
        ])
        self.assertEqual("openai-responses", args.llm_provider)

    def test_benchmark_uses_truth_only_after_panel_and_interval_match(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            config = root / "config.tsv"
            config.write_text(
                "case_id\ttrait\tanswer_gene_id\tphenotype_terms\n"
                "case1\tdwarf plant\tGeneTruth\tshort stem|plant height\n",
                encoding="utf-8",
            )
            reference = root / "reference.json"
            reference.write_text(json.dumps({
                "config": str(config),
                "cases": [{
                    "locus_id": "case1", "term_profile": "phenotype_only",
                    "extract": {"outputs": {"candidates_tsv": "candidate-input.tsv"}},
                    "official_protein_subset": {"output": "interval-proteins.faa"},
                    "evidence": {"evidence": "base-evidence.tsv"},
                }],
            }), encoding="utf-8")

            def fake_panel(phenotype, species, out_dir, **kwargs):
                self.assertEqual("dwarf plant", phenotype)
                self.assertNotIn("GeneTruth", json.dumps(kwargs))
                self.assertEqual([], kwargs["terms"])
                return {"panel": "panel.tsv", "fasta": "panel.faa", "manifest": "panel.json", "panel_count": 3}

            def fake_match(*args, **_kwargs):
                self.assertNotIn("GeneTruth", json.dumps(args))
                out = root / "fake_match"
                out.mkdir(exist_ok=True)
                ranked = out / "ranked.tsv"
                ranked.write_text(
                    "gene_id\tdry_score\tdry_group\tscore_breakdown\ttop_eligibility\n"
                    "GeneTruth\t1\tevidence_insufficient\tautomated_homolog_match/support:+0.75\tinsufficient_independent_evidence\n"
                    "GeneOther\t0\tevidence_insufficient\tno_scored_evidence:0\tinsufficient_independent_evidence\n",
                    encoding="utf-8",
                )
                shortlist = out / "shortlist.tsv"
                shortlist.write_text("gene_id\nGeneTruth\n", encoding="utf-8")
                evidence = out / "evidence.tsv"
                evidence.write_text("gene_id\nGeneTruth\n", encoding="utf-8")
                return {"ranked": str(ranked), "shortlist": str(shortlist), "evidence": str(evidence), "manifest": "match.json"}

            with patch("candidate_qtl.homolog_panel.automated_homolog_panel", side_effect=fake_panel), \
                    patch("candidate_qtl.homolog_panel.automated_homolog_match", side_effect=fake_match):
                result = benchmark_automated_homolog_panels(
                    str(reference), str(root / "out"), allow_external_api=True,
                )
            self.assertEqual(1.0, result["fixed_stage_a_20_recall"])
            self.assertFalse(result["answer_gene_used_before_evaluation"])


if __name__ == "__main__":
    unittest.main()
