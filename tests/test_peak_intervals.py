from __future__ import annotations

import io
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.peak_intervals import resolve_peak_context, support_interval_from_lod


CANDIDATES = (
    "gene_id\tchrom\tstart\tend\trepresentative_transcript\n"
    "GenePeak\tchr1\t2900\t3100\tGenePeak.1\n"
    "GeneSupport\tchr1\t2000\t2100\tGeneSupport.1\n"
    "GeneCore\tchr1\t2450\t2550\tGeneCore.1\n"
    "GeneOut\tchr1\t8000\t8100\tGeneOut.1\n"
)

LOD_PROFILE = (
    "chrom\tposition\tlod\n"
    "chr1\t1000\t1.0\n"
    "chr1\t2000\t3.0\n"
    "chr1\t3000\t5.0\n"
    "chr1\t4000\t3.0\n"
    "chr1\t5000\t1.0\n"
)


class SupportIntervalTests(unittest.TestCase):
    def test_lod_drop_support_interval(self) -> None:
        profile = pd.read_csv(io.StringIO(LOD_PROFILE), sep="\t")
        result = support_interval_from_lod(profile, drop=2.0)
        self.assertEqual(result["peak_position"], 3000)
        self.assertEqual(result["support_interval_start"], 2000)
        self.assertEqual(result["support_interval_end"], 4000)
        self.assertEqual(result["warnings"], [])

    def test_second_chromosome_peak_warns_linked_or_ghost(self) -> None:
        text = LOD_PROFILE + "chr2\t1000\t4.0\n"
        profile = pd.read_csv(io.StringIO(text), sep="\t")
        result = support_interval_from_lod(profile, drop=2.0)
        self.assertIn("linked_or_ghost_qtl_on_other_chromosome", result["warnings"])

    def test_interval_only_gives_low_confidence(self) -> None:
        context = resolve_peak_context(interval="chr1:1-100000")
        self.assertEqual(context.peak_confidence, "low")
        self.assertEqual(context.position_evidence_reason, "broad_interval_only")


class PeakContextCliTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.candidates = self.root / "candidates.tsv"
        self.candidates.write_text(CANDIDATES, encoding="utf-8")
        self.lod = self.root / "lod.tsv"
        self.lod.write_text(LOD_PROFILE, encoding="utf-8")

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_lod_profile_grades_position_evidence(self) -> None:
        out = self.root / "peak"
        self.assertEqual(
            0,
            main(
                [
                    "peak-context",
                    "--candidates", str(self.candidates),
                    "--lod-profile", str(self.lod),
                    "--lod-drop", "2.0",
                    "--out-dir", str(out),
                ]
            ),
        )
        context = (out / "peak_context.tsv").read_text(encoding="utf-8")
        self.assertIn("2_lod", context)
        self.assertIn("3000", context)

        evidence = pd.read_csv(out / "position_evidence.tsv", sep="\t").fillna("")
        codes = dict(zip(evidence["gene_id"], evidence["evidence_code"]))
        self.assertEqual(codes["GenePeak"], "peak_overlap")
        self.assertEqual(codes["GeneSupport"], "support_interval")
        # A gene only inside the broad interval receives no positive position row.
        self.assertNotIn("GeneOut", codes)

    def test_recombinant_core_is_high_confidence_and_audits_in_rank(self) -> None:
        metadata = self.root / "peak_meta.tsv"
        metadata.write_text(
            "chrom\tpeak_position\tpeak_method\trecombinant_core_start\trecombinant_core_end\tsource_id\n"
            "chr1\t3000\tinterval_mapping\t2400\t2600\tdoi:10.0/example\n",
            encoding="utf-8",
        )
        out = self.root / "peak_meta_out"
        self.assertEqual(
            0,
            main(
                [
                    "peak-context",
                    "--candidates", str(self.candidates),
                    "--peak-metadata", str(metadata),
                    "--out-dir", str(out),
                ]
            ),
        )
        evidence = pd.read_csv(out / "position_evidence.tsv", sep="\t").fillna("")
        core_rows = evidence[evidence["gene_id"] == "GeneCore"]
        self.assertEqual(list(core_rows["evidence_code"]), ["recombinant_core"])
        self.assertEqual(list(core_rows["confidence"]), ["high"])
        self.assertTrue((core_rows["source_id"] == "doi:10.0/example").all())

        # The recombinant-core contribution must be auditable in evidence-rank.
        rank_out = self.root / "rank"
        self.assertEqual(
            0,
            main(
                [
                    "evidence-rank",
                    "--candidates", str(self.candidates),
                    "--evidence", str(out / "position_evidence.tsv"),
                    "--out-dir", str(rank_out),
                ]
            ),
        )
        contributions = pd.read_csv(rank_out / "score_contributions.tsv", sep="\t").fillna("")
        core_contrib = contributions[
            (contributions["gene_id"] == "GeneCore") & (contributions["evidence_code"] == "recombinant_core")
        ]
        self.assertEqual(len(core_contrib), 1)
        self.assertGreater(float(core_contrib.iloc[0]["contribution"]), 0.0)

    def test_manifest_records_method_and_sources(self) -> None:
        out = self.root / "peak"
        main(
            [
                "peak-context",
                "--candidates", str(self.candidates),
                "--lod-profile", str(self.lod),
                "--lod-drop", "2.0",
                "--out-dir", str(out),
            ]
        )
        manifest = json.loads((out / "peak_context_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["peak_method"], "lod_profile")
        self.assertTrue(any(source.startswith("lod_profile:") for source in manifest["source_ids"]))
        self.assertIn("peak_overlap", manifest["position_evidence_counts"])


if __name__ == "__main__":
    unittest.main()
