from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from candidate_qtl.cli import main
from candidate_qtl.ranking_policy import LITERATURE_INFORMED_WEIGHTS, load_weight_registry, load_weights


class RankingPolicyTests(unittest.TestCase):
    def _rank(self, evidence_rows: str, candidates_rows: str | None = None) -> dict[str, dict[str, str]]:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            candidates = root / "candidates.tsv"
            candidates.write_text(
                candidates_rows
                or "gene_id\tannotation_quality\tcds_length\nGeneA\tcomplete_orf\t300\n",
                encoding="utf-8",
            )
            evidence = root / "evidence.tsv"
            evidence.write_text(
                "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\n"
                + evidence_rows,
                encoding="utf-8",
            )
            out = root / "rank"
            self.assertEqual(
                0,
                main(
                    [
                        "evidence-rank",
                        "--candidates",
                        str(candidates),
                        "--evidence",
                        str(evidence),
                        "--out-dir",
                        str(out),
                        "--top-n",
                        "3",
                        "--reserve-n",
                        "0",
                    ]
                ),
            )
            with (out / "ranked_candidates.tsv").open(encoding="utf-8") as handle:
                return {
                    row["gene_id"]: row
                    for row in csv.DictReader(handle, delimiter="\t")
                }

    def test_correlated_comparative_annotations_are_one_evidence_family(self) -> None:
        rows = self._rank(
            "GeneA\thomology_supported\tsupport\thigh\thit\tblast\tDB:1\thomology\n"
            "GeneA\torthogroup_supported\tsupport\thigh\tOG1\teggnog\tDB:2\torthogroup\n"
            "GeneA\tdomain_supported\tsupport\thigh\tPF1\thmmer\tDB:3\tdomain\n"
            "GeneA\ttrait_function_match\tsupport\thigh\tmeristem\tannotation\tDB:4\tfunction\n"
            "GeneA\tsynteny_supported\tsupport\thigh\tblock1\tjcvi\tDB:5\tsynteny\n"
        )
        self.assertEqual("1", rows["GeneA"]["support_family_count"])
        self.assertEqual("1.0", rows["GeneA"]["dry_score"])
        self.assertEqual("insufficient_independent_evidence", rows["GeneA"]["top_eligibility"])
        self.assertNotEqual("top", rows["GeneA"]["dry_group"])

    def test_comparative_function_plus_peak_from_independent_sources_can_enter_top(self) -> None:
        rows = self._rank(
            "GeneA\tdomain_supported\tsupport\thigh\tPF1\thmmer\tPMID:111\tdomain\n"
            "GeneA\tpeak_overlap\tsupport\thigh\t0\tqtl_position\tdoi:10.0/qtl\tpeak\n"
        )
        self.assertEqual("2", rows["GeneA"]["support_family_count"])
        self.assertEqual("independent_experimental_support", rows["GeneA"]["top_eligibility"])
        self.assertEqual("top", rows["GeneA"]["dry_group"])

    def test_two_channels_from_one_source_are_not_independent(self) -> None:
        # Same study/source for both channels: non-redundant, but not independent.
        rows = self._rank(
            "GeneA\tdomain_supported\tsupport\thigh\tPF1\thmmer\tPMID:111\tdomain\n"
            "GeneA\tpeak_overlap\tsupport\thigh\t0\tqtl_position\tPMID:111\tpeak\n"
        )
        self.assertEqual("2", rows["GeneA"]["support_family_count"])
        self.assertEqual("1", rows["GeneA"]["independent_source_count"])
        self.assertEqual("multiple_non_redundant_channels", rows["GeneA"]["top_eligibility"])
        self.assertNotEqual("top", rows["GeneA"]["dry_group"])

    def test_weight_profile_keeps_family_homology_weak(self) -> None:
        self.assertGreater(LITERATURE_INFORMED_WEIGHTS["direct_model_genetic"], LITERATURE_INFORMED_WEIGHTS["peak_overlap"])
        self.assertGreater(LITERATURE_INFORMED_WEIGHTS["peak_overlap"], LITERATURE_INFORMED_WEIGHTS["homology_supported"])
        self.assertEqual(0.0, LITERATURE_INFORMED_WEIGHTS["boundary_proximal"])

    def test_versioned_registry_is_default_and_contributions_are_traceable(self) -> None:
        registry, metadata = load_weight_registry()
        self.assertEqual("1.1.0", metadata["version"])
        self.assertEqual(64, len(metadata["sha256"]))
        self.assertEqual(LITERATURE_INFORMED_WEIGHTS, load_weights(None))
        self.assertFalse(registry["rationale"].str.strip().eq("").any())
        self.assertFalse(registry["reference_id"].str.strip().eq("").any())
        generic = registry.set_index("evidence_code")
        self.assertEqual("", generic.loc["homology_context", "base_weight_or_tier"])
        self.assertEqual("", generic.loc["domain_context", "base_weight_or_tier"])
        self.assertEqual("", generic.loc["go_annotation_context", "base_weight_or_tier"])

        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            candidates = root / "candidates.tsv"
            candidates.write_text("gene_id\tannotation_quality\tcds_length\nGeneA\tcomplete_orf\t300\n", encoding="utf-8")
            evidence = root / "evidence.tsv"
            evidence.write_text(
                "gene_id\tevidence_code\tdirection\tconfidence\tvalue\tsource_type\tsource_id\tdescription\n"
                "GeneA\tdirect_model_genetic\tsupport\thigh\tmutant\tliterature\tPMID:1\tcomplemented mutant\n",
                encoding="utf-8",
            )
            out = root / "rank"
            self.assertEqual(0, main(["evidence-rank", "--candidates", str(candidates), "--evidence", str(evidence), "--out-dir", str(out)]))
            contribution = next(csv.DictReader((out / "score_contributions.tsv").open(encoding="utf-8"), delimiter="\t"))
            self.assertEqual("causal_test", contribution["directness"])
            self.assertIn("ECO:0000016", contribution["reference_id"])
            decision = json.loads((out / "candidate_path_decision.json").read_text(encoding="utf-8"))
            self.assertEqual("1.1.0", decision["weight_registry"]["version"])
            self.assertIn("not a causal probability", decision["weight_interpretation"])

    def test_explicit_peak_window_has_higher_confidence_than_auto_window(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            candidates = root / "candidates.tsv"
            candidates.write_text(
                "gene_id\tchrom\tstart\tend\tinterval\tannotation_quality\tcds_length\n"
                "GeneA\tchr1\t490\t510\tchr1:1-1000\tcomplete_orf\t300\n"
                "GeneB\tchr1\t540\t560\tchr1:1-1000\tcomplete_orf\t300\n",
                encoding="utf-8",
            )
            explicit = root / "explicit"
            self.assertEqual(
                0,
                main(
                    [
                        "evidence-build",
                        "--candidates",
                        str(candidates),
                        "--out-dir",
                        str(explicit),
                        "--peak",
                        "chr1:500",
                        "--interval",
                        "chr1:1-1000",
                        "--peak-window",
                        "100",
                    ]
                ),
            )
            with (explicit / "evidence.tsv").open(encoding="utf-8") as handle:
                explicit_rows = list(csv.DictReader(handle, delimiter="\t"))
            overlap = next(row for row in explicit_rows if row["gene_id"] == "GeneA")
            proximal = next(row for row in explicit_rows if row["gene_id"] == "GeneB")
            self.assertEqual("peak_overlap", overlap["evidence_code"])
            self.assertEqual("high", overlap["confidence"])
            self.assertEqual("medium", proximal["confidence"])

            automatic = root / "automatic"
            self.assertEqual(
                0,
                main(
                    [
                        "evidence-build",
                        "--candidates",
                        str(candidates),
                        "--out-dir",
                        str(automatic),
                        "--peak",
                        "chr1:500",
                        "--interval",
                        "chr1:1-1000",
                    ]
                ),
            )
            with (automatic / "evidence.tsv").open(encoding="utf-8") as handle:
                automatic_rows = list(csv.DictReader(handle, delimiter="\t"))
            automatic_overlap = next(row for row in automatic_rows if row["gene_id"] == "GeneA")
            self.assertEqual("peak_overlap", automatic_overlap["evidence_code"])
            self.assertEqual("medium", automatic_overlap["confidence"])
            # Without an explicit window, the nearby gene gets an unscored context
            # row, not scored peak_proximal support.
            automatic_geneb = [row for row in automatic_rows if row["gene_id"] == "GeneB"]
            self.assertFalse(any(row["evidence_code"] == "peak_proximal" for row in automatic_geneb))
            heuristic = next(row for row in automatic_geneb if row["evidence_code"] == "peak_heuristic_proximal")
            self.assertEqual("context", heuristic["direction"])

    def test_peak_only_proximity_does_not_scale_with_interval_width(self) -> None:
        # The same peak precision must not give more scored support just because
        # the reported interval is wider.
        def proximal_support_codes(interval: str) -> set[str]:
            with tempfile.TemporaryDirectory() as tempdir:
                root = Path(tempdir)
                candidates = root / "candidates.tsv"
                candidates.write_text(
                    "gene_id\tchrom\tstart\tend\tinterval\tannotation_quality\tcds_length\n"
                    f"GeneNear\tchr1\t9000\t9100\t{interval}\tcomplete_orf\t300\n",
                    encoding="utf-8",
                )
                out = root / "out"
                self.assertEqual(0, main(["evidence-build", "--candidates", str(candidates), "--out-dir", str(out), "--peak", "chr1:5000", "--interval", interval]))
                with (out / "evidence.tsv").open(encoding="utf-8") as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                return {row["evidence_code"] for row in rows if row["direction"] == "support" and row["evidence_code"].startswith("peak")}

        # Narrow vs wide interval, same peak and gene: neither yields scored peak support.
        self.assertEqual(proximal_support_codes("chr1:4900-5100"), set())
        self.assertEqual(proximal_support_codes("chr1:1-100000"), set())


if __name__ == "__main__":
    unittest.main()
