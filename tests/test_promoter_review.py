from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from candidate_qtl import promoter_review as promoter_review_module
from candidate_qtl.cli import main


def _brute_force_intersect(regions: pd.DataFrame, features: pd.DataFrame, feature_type: str) -> list[dict[str, object]]:
    """The pre-#161 all-pairs scan, kept as the parity reference."""
    observations: list[dict[str, object]] = []
    if regions.empty or features.empty:
        return observations
    for region in regions.to_dict("records"):
        for feature in features.to_dict("records"):
            if str(feature["chrom"]) != str(region["chrom"]):
                continue
            overlap = promoter_review_module._overlap(
                int(region["start"]), int(region["end"]), int(feature["start"]), int(feature["end"]))
            if overlap <= 0:
                continue
            observations.append({"region_id": region["region_id"], "feature_id": str(feature.get("feature_id", "")),
                                 "start": int(feature["start"]), "end": int(feature["end"]), "overlap_bp": overlap})
    return observations


class PromoterReviewTests(unittest.TestCase):
    def _fixture(self, root: Path) -> dict[str, Path]:
        paths = {name: root / name for name in ("genes.tsv", "genes.gff3", "coverage.tsv", "variants.tsv", "sv.tsv", "motifs.tsv", "hypotheses.tsv")}
        paths["genes.tsv"].write_text("gene_id\trepresentative_transcript\nPlusGene\tPlus.1\nMinusGene\tMinus.1\n", encoding="utf-8")
        paths["genes.gff3"].write_text(
            "##gff-version 3\n"
            "chr1\tx\tgene\t100\t500\t.\t+\t.\tID=PlusGene\n"
            "chr1\tx\tmRNA\t100\t500\t.\t+\t.\tID=Plus.1;Parent=PlusGene\n"
            "chr1\tx\tfive_prime_UTR\t100\t129\t.\t+\t.\tParent=Plus.1\n"
            "chr1\tx\tCDS\t130\t450\t.\t+\t0\tParent=Plus.1\n"
            "chr1\tx\tthree_prime_UTR\t451\t500\t.\t+\t.\tParent=Plus.1\n"
            "chr1\tx\tmRNA\t120\t500\t.\t+\t.\tID=Plus.2;Parent=PlusGene\n"
            "chr1\tx\tfive_prime_UTR\t120\t129\t.\t+\t.\tParent=Plus.2\n"
            "chr1\tx\tgene\t1000\t1500\t.\t-\t.\tID=MinusGene\n"
            "chr1\tx\tmRNA\t1000\t1500\t.\t-\t.\tID=Minus.1;Parent=MinusGene\n"
            "chr1\tx\tfive_prime_UTR\t1470\t1500\t.\t-\t.\tParent=Minus.1\n"
            "chr1\tx\tCDS\t1050\t1469\t.\t-\t0\tParent=Minus.1\n", encoding="utf-8",
        )
        paths["coverage.tsv"].write_text("gene_id\tregion_class\tcoverage_status\nPlusGene\tcoding_exon\tpartial\n", encoding="utf-8")
        paths["variants.tsv"].write_text("feature_id\tchrom\tstart\tend\tflags\tsource_id\nSNV1\tchr1\t90\t90\t\tlocal:vcf\n", encoding="utf-8")
        paths["sv.tsv"].write_text("feature_id\tchrom\tstart\tend\tobservation_class\tflags\tsource_id\nSV1\tchr1\t80\t110\tsv\tlow_MAPQ\tlocal:sv\nTE1\tchr1\t1510\t1520\tte\trepeat\tlocal:te\n", encoding="utf-8")
        paths["motifs.tsv"].write_text("feature_id\tchrom\tstart\tend\tsource_id\nM1\tchr1\t95\t98\tlocal:motif\n", encoding="utf-8")
        paths["hypotheses.tsv"].write_text(
            "hypothesis_id\thypothesis_type\tmember_gene_ids\ttarget_gene_ids\n"
            "H_PROM\tpromoter_or_utr\tPlusGene\t\nH_SV\tsv_cnv_pav\tPlusGene\t\nH_TE\tte_regulatory\tMinusGene\t\n", encoding="utf-8",
        )
        return paths

    def test_strand_isoforms_common_intersection_gate_and_nonpromotion(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fixture = self._fixture(root); out = root / "out"
            argv = ["promoter-review", "--genes", str(fixture["genes.tsv"]), "--gff3", str(fixture["genes.gff3"]),
                    "--coverage", str(fixture["coverage.tsv"]), "--parental-variants", str(fixture["variants.tsv"]),
                    "--sv-calls", str(fixture["sv.tsv"]), "--motif-hits", str(fixture["motifs.tsv"]),
                    "--hypotheses", str(fixture["hypotheses.tsv"]), "--promoter-window", "50", "--out-dir", str(out)]
            self.assertEqual(0, main(argv))
            regions = pd.read_csv(out / "promoter_utr_regions.tsv", sep="\t", dtype=str).fillna("")
            plus = regions[(regions.transcript_id == "Plus.1") & (regions.region_class == "promoter")].iloc[0]
            minus = regions[(regions.transcript_id == "Minus.1") & (regions.region_class == "promoter")].iloc[0]
            self.assertEqual(("50", "99"), (plus.start, plus.end))
            self.assertEqual(("1501", "1550"), (minus.start, minus.end))
            self.assertEqual("20", regions[regions.gene_id == "PlusGene"].iloc[0].tss_uncertainty_bp)
            self.assertEqual({"Plus.1", "Plus.2", "Minus.1"}, set(regions.transcript_id))
            observations = pd.read_csv(out / "review_observations.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(observations.causal_promotion.eq("none").all())
            self.assertEqual("context", observations.loc[observations.feature_type == "reference_motif", "review_status"].iloc[0])
            self.assertIn("sv_cnv_pav", set(observations.links_to_hypothesis_type))
            self.assertIn("te_regulatory", set(observations.links_to_hypothesis_type))
            links = pd.read_csv(out / "hypothesis_observation_links.tsv", sep="\t", dtype=str).fillna("")
            self.assertTrue(links.hypothesis_id.ne("").all())
            self.assertTrue(links.direction.eq("context").all())
            manifest = json.loads((out / "promoter_review_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(manifest["route_gate"]["route_started"])
            self.assertEqual(["MinusGene"], manifest["route_gate"]["unexamined_genes"])
            self.assertNotIn("dry_score", observations.columns)

    def test_missing_coding_coverage_stops_without_silent_force(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fixture = self._fixture(root); out = root / "out"
            self.assertEqual(0, main(["promoter-review", "--genes", str(fixture["genes.tsv"]), "--gff3", str(fixture["genes.gff3"]), "--out-dir", str(out)]))
            manifest = json.loads((out / "promoter_review_manifest.json").read_text(encoding="utf-8"))
            self.assertFalse(manifest["started"])
            self.assertFalse((out / "promoter_utr_regions.tsv").exists())


class FeatureIndexTests(unittest.TestCase):
    """The interval index must reproduce the all-pairs scan exactly (#161)."""

    def _synthetic(self, region_count: int, feature_count: int) -> tuple[pd.DataFrame, pd.DataFrame]:
        chroms = ["chr1", "chr2", "chr3"]
        regions = pd.DataFrame([{
            "region_id": f"reg_{i:04d}", "chrom": chroms[i % len(chroms)],
            "start": 1 + i * 977, "end": 1 + i * 977 + (50 if i % 3 else 4000),
        } for i in range(region_count)])
        features = pd.DataFrame([{
            "feature_id": f"f_{i:06d}", "chrom": chroms[i % len(chroms)],
            # A mix of point calls, short indels, and rare long SVs, plus a
            # deliberately reversed pair that the scan silently skipped.
            "start": 1 + (i * 37) % 100_000,
            "end": 1 + (i * 37) % 100_000 + (0 if i % 4 else (3 if i % 7 else 9000)) - (5 if i % 5000 == 4999 else 0),
        } for i in range(feature_count)])
        return regions, features

    def test_index_matches_all_pairs_scan(self) -> None:
        regions, features = self._synthetic(120, 20_000)
        expected = _brute_force_intersect(regions, features, "parental_variant")
        index = promoter_review_module._FeatureIndex(features)
        actual: list[dict[str, object]] = []
        for region in regions.to_dict("records"):
            start, end = int(region["start"]), int(region["end"])
            for feature in index.overlapping(str(region["chrom"]), start, end):
                overlap = promoter_review_module._overlap(start, end, int(feature["start"]), int(feature["end"]))
                if overlap <= 0:
                    continue
                actual.append({"region_id": region["region_id"], "feature_id": str(feature["feature_id"]),
                               "start": int(feature["start"]), "end": int(feature["end"]), "overlap_bp": overlap})
        self.assertEqual(expected, actual)
        self.assertTrue(expected, "the synthetic fixture must produce overlaps")

    def test_intersection_does_not_test_every_region_feature_pair(self) -> None:
        regions, features = self._synthetic(100, 60_000)
        regions = regions.assign(gene_id="G", transcript_id="T", region_class="promoter")
        real_overlap = promoter_review_module._overlap
        calls = 0

        def counting_overlap(*args: int) -> int:
            nonlocal calls
            calls += 1
            return real_overlap(*args)

        with mock.patch.object(promoter_review_module, "_overlap", side_effect=counting_overlap):
            observations = promoter_review_module._intersect(regions, features, "parental_variant")
        all_pairs = len(regions) * len(features)
        self.assertTrue(observations)
        # The all-pairs scan would evaluate 6,000,000 pairs; the index must stay
        # orders of magnitude below that.
        self.assertLess(calls, all_pairs // 100)


class FeatureStreamingTests(unittest.TestCase):
    """Genome-wide feature tables must not be materialised in full (#161)."""

    def _table(self, path: Path, rows: int) -> None:
        lines = ["feature_id\tchrom\tstart\tend\tsource_id"]
        for i in range(rows):
            chrom = "chr1" if i % 2 == 0 else "chrOther"
            position = 10 + i * 13
            lines.append(f"F{i}\t{chrom}\t{position}\t{position}\tlocal:vcf")
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    def test_rows_outside_the_reviewed_span_are_dropped_while_streaming(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "variants.tsv"
            self._table(table, 5_000)
            regions = pd.DataFrame([{"chrom": "chr1", "start": 100, "end": 400}])
            bounds = promoter_review_module._region_bounds(regions)
            with mock.patch.object(promoter_review_module, "FEATURE_CHUNK_ROWS", 128):
                frame = promoter_review_module._load_features(str(table), "parental_variant", bounds)
            self.assertEqual(5_000, frame.attrs["rows_read"])
            self.assertLess(len(frame), 30)
            self.assertEqual({"chr1"}, set(frame["chrom"]))
            self.assertTrue((pd.to_numeric(frame["start"]) <= 400).all())
            self.assertTrue((pd.to_numeric(frame["end"]) >= 100).all())
            # Unfiltered loading is unchanged for callers without region context.
            full = promoter_review_module._load_features(str(table), "parental_variant")
            self.assertEqual(5_000, len(full))

    def test_malformed_coordinates_still_raise(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "variants.tsv"
            table.write_text("feature_id\tchrom\tstart\tend\nF1\tchr1\tnot_a_number\t120\n", encoding="utf-8")
            bounds = {"chr1": (1, 1000)}
            with self.assertRaises(ValueError):
                promoter_review_module._load_features(str(table), "parental_variant", bounds)

    def test_manifest_records_the_intersection_strategy(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixture = PromoterReviewTests()._fixture(root)
            out = root / "out"
            self.assertEqual(0, main(["promoter-review", "--genes", str(fixture["genes.tsv"]), "--gff3", str(fixture["genes.gff3"]),
                                      "--coverage", str(fixture["coverage.tsv"]), "--parental-variants", str(fixture["variants.tsv"]),
                                      "--sv-calls", str(fixture["sv.tsv"]), "--promoter-window", "50", "--out-dir", str(out)]))
            manifest = json.loads((out / "promoter_review_manifest.json").read_text(encoding="utf-8"))
            intersection = manifest["feature_intersection"]
            self.assertIn("chr1", intersection["region_bounds"])
            self.assertEqual(1, intersection["row_scan"]["parental_variants"]["rows_read"])
            self.assertEqual(2, intersection["row_scan"]["sv_calls"]["rows_read"])
            self.assertIn("not an all-pairs scan", intersection["method"])
            self.assertFalse(intersection["tabix_index_used"])


if __name__ == "__main__":
    unittest.main()
