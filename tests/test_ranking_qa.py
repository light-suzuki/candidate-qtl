from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.hardening import prioritize_evidence

CANDIDATES = (
    "gene_id\tannotation_quality\tcds_length\n"
    + "\n".join(f"Gene{name}\tcomplete_orf\t300" for name in
                ("Stable", "Source", "Provider", "Mid", "TieA", "TieB", "Edge", "Solo", "Generic"))
    + "\n"
)

EVIDENCE = (
    "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\n"
    "GeneStable\tdirect_model_genetic\tsupport\thigh\tmutant\tgenetics\tPMID:10\tdirect mutant\n"
    "GeneStable\tpeak_overlap\tsupport\thigh\t0\tmapping\tPMID:11\tpeak\n"
    "GeneStable\texpression_target_specific\tsupport\thigh\t4x\trnaseq\tPMID:12\tRNA\n"
    "GeneStable\texpression_target_specific\tsupport\thigh\t4x\teggnog\tPMID:13\tRNA\n"
    "GeneSource\tpeak_overlap\tsupport\thigh\t0\tmapping\tPMID:20\tpeak\n"
    "GeneSource\texpression_target_specific\tsupport\thigh\t4x\trnaseq\tPMID:21\tRNA\n"
    "GeneProvider\texpression_target_specific\tsupport\thigh\t4x\teggnog\tPMID:30\tRNA\n"
    "GeneProvider\tpeak_overlap\tsupport\thigh\t0\teggnog\tPMID:31\tpeak\n"
    "GeneProvider\tdomain_supported\tsupport\thigh\tPF1\teggnog\tPMID:32\tdomain\n"
    "GeneMid\tpeak_proximal\tsupport\thigh\t500\tmapping\tPMID:40\tproximal\n"
    "GeneTieA\texpression_target_specific\tsupport\thigh\t4x\trnaseq\tPMID:50\tRNA\n"
    "GeneTieB\texpression_target_specific\tsupport\thigh\t4x\trnaseq\tPMID:51\tRNA\n"
    "GeneEdge\texpression_target_specific\tsupport\thigh\t4x\trnaseq\tPMID:52\tRNA\n"
    "GeneSolo\texpression_target_specific\tsupport\thigh\t4x\trnaseq\tPMID:53\tRNA\n"
    "GeneGeneric\thomology_context\tcontext\thigh\thit1\tblast\tDB:1\tgeneric\n"
    "GeneGeneric\tdomain_context\tcontext\thigh\tPF2\thmmer\tDB:2\tgeneric\n"
)


def _run_qa(root: Path, evidence_text: str = EVIDENCE, candidates_text: str = CANDIDATES, top_k: int = 3,
            extra_args: list[str] | None = None) -> tuple[Path, Path, Path]:
    candidates = root / "candidates.tsv"
    evidence = root / "evidence.tsv"
    out = root / "qa"
    candidates.write_text(candidates_text, encoding="utf-8")
    evidence.write_text(evidence_text, encoding="utf-8")
    args = ["ranking-qa", "--candidates", str(candidates), "--evidence", str(evidence), "--out-dir", str(out), "--top-k", str(top_k)]
    if extra_args:
        args.extend(extra_args)
    main(args)
    return candidates, evidence, out


class RankingQaStabilityClassTests(unittest.TestCase):
    def test_adversarial_ranking_paths_are_reported_without_overwriting_default(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); candidates = root / "candidates.tsv"; evidence = root / "evidence.tsv"; out = root / "qa"
            candidates.write_text(
                "gene_id\tannotation_quality\tcds_length\n"
                "GeneGenetic\tcomplete_orf\t300\nGenePosition\tcomplete_orf\t300\n"
                "GeneExpression\tcomplete_orf\t300\nGeneGeneric\tcomplete_orf\t300\nGeneConflict\tcomplete_orf\t300\n",
                encoding="utf-8",
            )
            evidence.write_text(
                "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\n"
                "GeneGenetic\tdirect_model_genetic\tsupport\thigh\tmutant\tgenetics\tPMID:1\tdirect mutant\n"
                "GenePosition\tpeak_overlap\tsupport\thigh\t0\tmapping\tPMID:2\tpeak only\n"
                "GeneExpression\texpression_target_specific\tsupport\thigh\t4x\trnaseq\tPMID:3\texpression only\n"
                "GeneGeneric\thomology_context\tcontext\thigh\thit1\tblast\tDB:1\tgeneric hit\n"
                "GeneGeneric\tdomain_context\tcontext\thigh\tPF1\thmmer\tDB:2\tgeneric domain\n"
                "GeneGeneric\tgo_annotation_context\tcontext\thigh\tGO:1\teggnog\tDB:3\tgeneric GO\n"
                "GeneConflict\ttrait_function_match\tsupport\thigh\ttrait\tcuration\tPMID:4\ttrait support\n"
                "GeneConflict\tcounterevidence\tagainst\thigh\tcontradiction\tgenetics\tPMID:5\tcontradiction\n",
                encoding="utf-8",
            )
            self.assertEqual(0, main(["ranking-qa", "--candidates", str(candidates), "--evidence", str(evidence),
                                      "--out-dir", str(out), "--top-k", "2", "--benchmark-mode"]))
            for name in ("rank_concordance.tsv", "ranking_sensitivity.tsv", "candidate_rank_distribution.tsv",
                         "consensus_shortlist.tsv", "ranking_discordance_review.tsv", "ranking_qa_manifest.json",
                         "rank_stability.tsv", "source_dropout_stability.tsv", "weight_perturbation_summary.tsv"):
                self.assertTrue((out / name).exists(), name)
            distribution = pd.read_csv(out / "candidate_rank_distribution.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            self.assertEqual("1", distribution.loc["GeneGenetic", "default_rank"])
            self.assertIn("single_family_dependent", distribution.loc["GenePosition", "instability_reason"])
            self.assertGreater(int(distribution.loc["GenePosition", "family_dropout_max_shift"]), 0)
            manifest = json.loads((out / "ranking_qa_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(manifest["benchmark_mode"])
            self.assertIn("never overwrites", manifest["interpretation"])
            self.assertIn("probabilit", manifest["interpretation"])
            self.assertIn("recommended weights", manifest["interpretation"])
            sensitivity = pd.read_csv(out / "ranking_sensitivity.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(sensitivity["scenario"].str.startswith("perturb:").any())
            self.assertTrue(sensitivity["scenario"].str.startswith("provider:").any())

    def test_stability_classes_are_distinguished(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, _, out = _run_qa(root)
            distribution = pd.read_csv(out / "candidate_rank_distribution.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            expected = {"GeneStable": "stable", "GeneSource": "source_sensitive",
                        "GeneProvider": "provider_sensitive", "GeneTieA": "underdetermined_ties"}
            for gene_id, stability_class in expected.items():
                self.assertEqual(stability_class, distribution.loc[gene_id, "stability_class"], gene_id)
            self.assertEqual("4", distribution.loc["GeneTieA", "screening_tie_size"])
            self.assertEqual("1", distribution.loc["GeneStable", "default_rank"])
            self.assertEqual("", distribution.loc["GeneStable", "instability_reason"].replace("provider_discordant", ""))
            self.assertIn("source_group_dependent", distribution.loc["GeneSource", "instability_reason"])
            self.assertIn("provider_dependent", distribution.loc["GeneProvider", "instability_reason"])
            self.assertNotIn("source_group_dependent", distribution.loc["GeneProvider", "instability_reason"])

    def test_weight_sensitive_class_appears_when_perturbation_crosses_top_k(self) -> None:
        # GeneW's two peak_overlap rows come from two sources and two providers,
        # so no single source/provider dropout can displace it -- only a weight
        # change (or its position-family dropout) pushes it out of top_k=1.
        candidates = "gene_id\tannotation_quality\tcds_length\nGeneW\tcomplete_orf\t300\nGeneB\tcomplete_orf\t300\nGeneC\tcomplete_orf\t300\n"
        evidence = (
            "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\n"
            "GeneW\tpeak_overlap\tsupport\thigh\t0\tmapping\tPMID:1\tpeak\n"
            "GeneW\tpeak_overlap\tsupport\thigh\t0\tqtl_position\tPMID:2\tpeak\n"
            "GeneB\tpeak_proximal\tsupport\thigh\t500\tmapping\tPMID:3\tproximal\n"
            "GeneC\texpression_target_specific\tsupport\thigh\t4x\trnaseq\tPMID:4\tRNA\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, _, out = _run_qa(root, candidates_text=candidates, evidence_text=evidence, top_k=1)
            distribution = pd.read_csv(out / "candidate_rank_distribution.tsv", sep="\t", dtype=str).fillna("").set_index("gene_id")
            self.assertEqual("1", distribution.loc["GeneW", "default_rank"])
            self.assertEqual("weight_sensitive", distribution.loc["GeneW", "stability_class"])
            self.assertIn("weight_sensitive", distribution.loc["GeneW", "instability_reason"])
            self.assertNotIn("source_group_dependent", distribution.loc["GeneW", "instability_reason"])
            self.assertNotIn("provider_dependent", distribution.loc["GeneW", "instability_reason"])
            dropout = pd.read_csv(out / "source_dropout_stability.tsv", sep="\t", dtype=str).fillna("")
            w_dropout = dropout[dropout["gene_id"].eq("GeneW")]
            self.assertTrue((w_dropout["crossed_top_k"] == "False").all())

    def test_source_dropout_uses_resolved_groups_and_skips_empty_group(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, _, out = _run_qa(root)
            dropout = pd.read_csv(out / "source_dropout_stability.tsv", sep="\t", dtype=str).fillna("")
            labels = set(dropout["dropped_label"])
            # The two uncited context rows from one BLAST/domain pipeline collapse
            # to one shared group, never to their individual accession-like ids.
            self.assertIn("comparative_annotation_chain", labels)
            self.assertNotIn("DB:1", labels)
            self.assertNotIn("DB:2", labels)
            # Every citable source anchors its own group.
            self.assertIn("PMID:10", labels)
            # The empty group ("") is never a dropout target.
            self.assertNotIn("", labels)
            scenario_names = set(dropout["scenario"])
            self.assertFalse(any(name.endswith(":") for name in scenario_names))
            sensitivity = pd.read_csv(out / "ranking_sensitivity.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(sensitivity["scenario"].str.startswith("source_dropout:").any())
            self.assertTrue(sensitivity["scenario"].str.startswith("provider_dropout:").any())

    def test_explicit_independence_group_dropout_removes_whole_group(self) -> None:
        candidates = "gene_id\tannotation_quality\tcds_length\nGeneA\tcomplete_orf\t300\nGeneB\tcomplete_orf\t300\nGeneC\tcomplete_orf\t300\nGeneD\tcomplete_orf\t300\n"
        evidence = (
            "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tindependence_group\tdescription\n"
            "GeneA\tdirect_model_genetic\tsupport\thigh\tmutant\tgenetics\tPMID:60\tstudyA\tdirect mutant\n"
            "GeneB\tpeak_overlap\tsupport\thigh\t0\tmapping\tPMID:61\tstudyA\tpeak\n"
            "GeneB\texpression_detected\tsupport\thigh\t4x\trnaseq\tPMID:63\tstudyB\tRNA\n"
            "GeneC\texpression_target_specific\tsupport\thigh\t4x\trnaseq\tPMID:62\tstudyB\tRNA\n"
            "GeneD\tdomain_supported\tsupport\thigh\tPF1\thmmer\tPMID:64\tstudyB\tdomain\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, _, out = _run_qa(root, candidates_text=candidates, evidence_text=evidence, top_k=2)
            dropout = pd.read_csv(out / "source_dropout_stability.tsv", sep="\t", dtype=str).fillna("")
            group_rows = dropout[(dropout["scenario"].eq("source_dropout:studyA")) & dropout["gene_id"].isin({"GeneA", "GeneB"})]
            self.assertEqual(2, len(group_rows), group_rows.to_string())
            # One group removal drops both genes' studyA rows below top_k=2.
            for gene_id in ("GeneA", "GeneB"):
                row = group_rows[group_rows["gene_id"].eq(gene_id)].iloc[0]
                self.assertEqual("studyA", row["dropped_label"])
                self.assertEqual("True", row["crossed_top_k"])
                self.assertGreater(int(row["scenario_rank"]), 2)
            # The other explicit group is untouched by the studyA scenario, and
            # vice versa: studyB dropout leaves GeneA/GeneB exactly at default.
            other = dropout[dropout["scenario"].eq("source_dropout:studyB")].set_index("gene_id")
            for gene_id, expected_rank in (("GeneA", "1"), ("GeneB", "2")):
                self.assertEqual(expected_rank, other.loc[gene_id, "default_rank"], gene_id)
                self.assertEqual(expected_rank, other.loc[gene_id, "scenario_rank"], gene_id)
                self.assertEqual("False", other.loc[gene_id, "crossed_top_k"], gene_id)

    def test_label_sources_never_anchor_a_dropout_group(self) -> None:
        candidates = "gene_id\tannotation_quality\tcds_length\nGeneLabel\tcomplete_orf\t300\nGenePeer\tcomplete_orf\t300\n"
        evidence = (
            "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\n"
            "GeneLabel\ttrait_function_match\tsupport\thigh\ttrait\tcuration\tlit1\tlabel only\n"
            "GenePeer\texpression_target_specific\tsupport\thigh\t4x\trnaseq\tPMID:70\tRNA\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, _, out = _run_qa(root, candidates_text=candidates, evidence_text=evidence, top_k=1)
            dropout = pd.read_csv(out / "source_dropout_stability.tsv", sep="\t", dtype=str).fillna("")
            labels = set(dropout["dropped_label"])
            self.assertNotIn("lit1", labels)
            self.assertNotIn("", labels)
            self.assertIn("PMID:70", labels)

    def test_dropout_filters_restrict_scenarios(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, _, out = _run_qa(root, extra_args=["--source-dropout", "comparative_annotation_chain",
                                                  "--provider-dropout", "blast"])
            dropout = pd.read_csv(out / "source_dropout_stability.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(["provider_dropout:blast", "source_dropout:comparative_annotation_chain"],
                             sorted(set(dropout["scenario"])))
            manifest = json.loads((out / "ranking_qa_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(["comparative_annotation_chain"], manifest["source_dropout_groups"])
            self.assertEqual(["blast"], manifest["provider_dropout_types"])

    def test_new_output_schemas_and_perturbation_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, _, out = _run_qa(root)
            rank_stability = pd.read_csv(out / "rank_stability.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(["gene_id", "default_rank", "rank_min", "rank_max", "rank_range", "tie_size",
                              "stability_class", "instability_reason", "top_k_rate", "scenario_shift_max"],
                             list(rank_stability.columns))
            by_gene = rank_stability.set_index("gene_id")
            self.assertEqual("1", by_gene.loc["GeneStable", "default_rank"])
            self.assertEqual("stable", by_gene.loc["GeneStable", "stability_class"])
            self.assertEqual(int(by_gene.loc["GeneStable", "rank_max"]) - int(by_gene.loc["GeneStable", "rank_min"]),
                             int(by_gene.loc["GeneStable", "rank_range"]))
            self.assertGreater(int(by_gene.loc["GeneProvider", "rank_range"]), 0)
            dropout = pd.read_csv(out / "source_dropout_stability.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(["gene_id", "scenario", "scenario_type", "dropped_label", "default_rank",
                              "scenario_rank", "rank_shift", "crossed_top_k"], list(dropout.columns))
            self.assertGreater(len(dropout), 0)
            summary = pd.read_csv(out / "weight_perturbation_summary.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(["evidence_code", "perturbation_fraction", "direction", "max_rank_shift",
                              "median_rank_shift", "candidates_affected", "worst_affected_gene"],
                             list(summary.columns))
            # 5 weighted codes x 3 fractions x 2 directions.
            self.assertEqual(30, len(summary))
            midpoint = summary[(summary["evidence_code"].eq("peak_proximal")) & (summary["direction"].eq("minus"))
                               & (summary["perturbation_fraction"].eq("0.5"))].iloc[0]
            self.assertEqual("GeneMid", midpoint["worst_affected_gene"])
            self.assertGreater(int(midpoint["max_rank_shift"]), 0)
            self.assertGreater(int(midpoint["candidates_affected"]), 0)
            # The summary reports rank movement only -- never weight deltas or scores.
            self.assertNotIn("weight", "".join(summary.columns))

    def test_canonical_default_rank_is_unchanged(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            candidates, evidence, out = _run_qa(root)
            before = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in (candidates, evidence)}
            sensitivity = pd.read_csv(out / "ranking_sensitivity.tsv", sep="\t", dtype=str).fillna("")
            qa_default = sensitivity[sensitivity["scenario"].eq("default_weighted")].set_index("gene_id")["rank"].astype(int)
            rank_dir = root / "fresh_rank"
            prioritize_evidence(str(candidates), str(evidence), str(rank_dir), None, 3, 3, 0.0, 3)
            fresh = pd.read_csv(rank_dir / "ranked_candidates.tsv", sep="\t", dtype=str).fillna("")
            fresh_rank = fresh.set_index("gene_id")["screening_rank"].astype(int)
            self.assertEqual(set(fresh_rank.index), set(qa_default.index))
            self.assertEqual({gene: int(value) for gene, value in fresh_rank.items()},
                             {gene: int(value) for gene, value in qa_default.items()})
            after = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in (candidates, evidence)}
            self.assertEqual(before, after)


if __name__ == "__main__":
    unittest.main()
