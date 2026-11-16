from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from candidate_qtl.hardening import (
    annotation_keyword_evidence,
    llm_evidence_import,
    prioritize_evidence,
    research_packet,
)
from candidate_qtl.cli import main


class EvidenceHardeningTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_keyword_evidence_requires_quality_homology(self) -> None:
        candidates = self.root / "candidates.tsv"
        candidates.write_text("gene_id\trepresentative_transcript\nGeneA\tTxA\n", encoding="utf-8")
        blast = self.root / "blast.tsv"
        blast.write_text(
            "GeneA\tTargetA\t20\t20\t100\t100\t1\t20\t1\t20\t1e-2\t50\tweak hit\n",
            encoding="utf-8",
        )
        gff = self.root / "target.gff3"
        gff.write_text("chr1\tx\tgene\t1\t100\t.\t+\t.\tID=TargetA;Note=starch%20branching%20enzyme\n", encoding="utf-8")
        output = self.root / "keyword.tsv"
        result = annotation_keyword_evidence(str(candidates), str(blast), str(gff), ["starch"], str(output))
        self.assertEqual(0, result["matched_candidates"])
        self.assertEqual(0, result["quality_pass_count"])
        audit = (self.root / "keyword_audit.tsv").read_text(encoding="utf-8")
        self.assertIn("fail", audit)

    def test_conflicts_block_top_and_rescue_batch_is_bounded(self) -> None:
        candidates = self.root / "candidates.tsv"
        candidates.write_text(
            "gene_id\tannotation_quality\tcds_length\n"
            "GeneConflict\tcomplete_orf\t300\n"
            "GeneWeak\tcomplete_orf\t300\n"
            "GenePeak\tcomplete_orf\t300\n"
            "GeneZero\tcomplete_orf\t300\n",
            encoding="utf-8",
        )
        evidence = self.root / "evidence.tsv"
        evidence.write_text(
            "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\n"
            "GeneConflict\tdirect_model_genetic\tsupport\thigh\tmutant\tmanual\tPMID:1\tsupport\n"
            "GeneConflict\tdirect_model_genetic\tagainst\thigh\tnormal\tmanual\tPMID:2\tcounterexample\n"
            "GeneWeak\thomology_supported\tsupport\thigh\tHit1\tblast\tblast.tsv\tconserved unknown\n"
            "GenePeak\thomology_supported\tsupport\thigh\tHit2\tblast\tblast.tsv\tconserved unknown\n"
            "GenePeak\tpeak_proximal\tsupport\tmedium\t100\tqtl_position\tchr1:500\tnear peak\n",
            encoding="utf-8",
        )
        out = self.root / "rank"
        result = prioritize_evidence(str(candidates), str(evidence), str(out), top_n=1, reserve_n=0, rescue_n=1)
        self.assertEqual(1, result["groups"].get("rescue", 0))
        with (out / "ranked_candidates.tsv").open(encoding="utf-8") as handle:
            ranked = {row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual("conflicting_evidence", ranked["GeneConflict"]["top_eligibility"])
        self.assertNotEqual("top", ranked["GeneConflict"]["dry_group"])
        self.assertEqual("rescue", ranked["GenePeak"]["dry_group"])
        self.assertEqual("evidence_insufficient", ranked["GeneWeak"]["dry_group"])
        self.assertEqual("1", ranked["GenePeak"]["screening_rank"])
        self.assertIn("screening_rank_max", ranked["GeneWeak"])
        self.assertIn("screening_tie_size", ranked["GeneWeak"])
        contributions = (out / "score_contributions.tsv").read_text(encoding="utf-8")
        self.assertIn("support", contributions)
        self.assertIn("against", contributions)

    def test_evidence_rank_cli_records_optional_case_id(self) -> None:
        candidates = self.root / "candidates.tsv"
        candidates.write_text("gene_id\tannotation_quality\nGeneA\tcomplete_orf\n", encoding="utf-8")
        evidence = self.root / "evidence.tsv"
        evidence.write_text(
            "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\n"
            "GeneA\tknown_qtg\tsupport\thigh\tknown\tliterature\tPMID:1\tdirect evidence\n",
            encoding="utf-8",
        )
        out = self.root / "rank"
        self.assertEqual(0, main(["evidence-rank", "--candidates", str(candidates), "--evidence", str(evidence),
                                  "--out-dir", str(out), "--case-id", "pea_r"]))
        decision = json.loads((out / "candidate_path_decision.json").read_text(encoding="utf-8"))
        self.assertEqual("pea_r", decision["case_id"])

    def test_llm_import_requires_review_and_source(self) -> None:
        candidates = self.root / "candidates.tsv"
        candidates.write_text("gene_id\nGeneA\n", encoding="utf-8")
        payload = self.root / "llm.json"
        payload.write_text(
            json.dumps(
                [
                    {
                        "gene_id": "GeneA",
                        "claims": [
                            {
                                "evidence_code": "ortholog_model_genetic",
                                "direction": "support",
                                "confidence": "high",
                                "value": "mutant phenotype",
                                "source_id": "PMID:123",
                                "description": "ortholog mutant matches trait",
                                "reviewed": True,
                                "reviewed_by": "Reviewer A",
                            },
                            {
                                "evidence_code": "known_qtg",
                                "direction": "support",
                                "confidence": "high",
                                "value": "claim",
                                "source_id": "",
                                "description": "no source",
                                "reviewed": True,
                                "reviewed_by": "Reviewer A",
                            },
                            {
                                "evidence_code": "trait_function_match",
                                "direction": "support",
                                "confidence": "medium",
                                "value": "claim",
                                "source_id": "PMID:999",
                                "description": "not reviewed",
                                "reviewed": False,
                            },
                            {
                                # reviewed=true is a self-declaration; without a
                                # reviewer it is not attributable to anyone.
                                "evidence_code": "trait_function_match",
                                "direction": "support",
                                "confidence": "medium",
                                "value": "claim",
                                "source_id": "PMID:456",
                                "description": "reviewed by nobody in particular",
                                "reviewed": True,
                            },
                        ],
                    }
                ]
            ),
            encoding="utf-8",
        )
        output = self.root / "reviewed.tsv"
        result = llm_evidence_import(str(candidates), str(payload), str(output))
        self.assertEqual(1, result["accepted_count"])
        self.assertEqual(3, result["rejected_count"])
        self.assertIn("PMID:123", output.read_text(encoding="utf-8"))
        self.assertNotIn("PMID:456", output.read_text(encoding="utf-8"))

    def test_research_packet_embeds_canonical_evidence(self) -> None:
        candidates = self.root / "ranked.tsv"
        candidates.write_text("gene_id\tdry_group\tdry_rank\tdry_score\nGeneA\ttop\t1\t3\n", encoding="utf-8")
        evidence = self.root / "evidence.tsv"
        evidence.write_text(
            "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\n"
            "GeneA\trbh_ortholog\tsupport\thigh\tAt1g01010\treciprocal_blast\trbh.tsv\tRBH\n",
            encoding="utf-8",
        )
        contributions = self.root / "contributions.tsv"
        contributions.write_text(
            "gene_id\tevidence_code\tdirection\tcontribution\tconfidence\tsource_id\n"
            "GeneA\trbh_ortholog\tsupport\t2\thigh\trbh.tsv\n",
            encoding="utf-8",
        )
        result = research_packet(
            str(candidates), str(self.root / "packet"), "Pisum sativum", "seed shape", 1, False, 10,
            str(evidence), str(contributions),
        )
        payload = json.loads(Path(result["packet"]).read_text(encoding="utf-8"))
        self.assertEqual("rbh_ortholog", payload["candidates"][0]["evidence"][0]["evidence_code"])
        self.assertEqual("2", payload["candidates"][0]["score_contributions"][0]["contribution"])


if __name__ == "__main__":
    unittest.main()
