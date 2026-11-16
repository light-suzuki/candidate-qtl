"""Real BLAST+ end-to-end tests for the blast-panel flow (issue #374 / #355).

Runs ``blast-panel-freeze`` -> ``blast-panel-search`` on a COPY of the
committed route_a fixtures (the fixture dir is never mutated). Skipped unless
BLAST+ (blastp + makeblastdb) is on PATH; the CI `integration-blast` job
installs ncbi-blast+ and runs it for real. Assertions are limited to stable
mechanics: provider status, ambiguity classes, hit/paralog counts and the
rank_contribution=none contract -- never evalue or bitscore absolute values.
"""

from __future__ import annotations

import json
import shutil
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main
from candidate_qtl.blast_panel import PANEL_COLUMNS, panel_from_verified_bait_panel, validate_panel_row

FIXTURE = Path(__file__).resolve().parent / "fixtures" / "route_a"

INVALID_ACCESSIONS = ("Q9TES1", "P9TES2", "Q9TES3", "Q9TES4")


class PanelConversionTests(unittest.TestCase):
    """The bait-panel -> blast-panel converter (guided triple-route #374).

    Pure conversion, no BLAST+ required: the llm_pathway verified bait panel is
    re-normalized into the blast-panel freeze input schema, and every converted
    row passes the real row validation (so a freeze yields verified baits).
    """

    def _panel_tsv(self, root: Path) -> Path:
        path = root / "verified_bait_panel.tsv"
        path.write_text("\t".join([
            "hypothesis_id", "bait_id", "pathway_or_process", "gene_family", "reference_gene_name",
            "reference_species", "expected_role", "sequence_length", "sequence_sha256",
            "source_provider", "source_url", "retrieved_at", "function_source_ids",
            "rank_contribution", "panel_status", "review_status",
        ]) + "\n" + "\t".join([
            "hyp_001", "P30924", "meristem maintenance", "CLV", "CLV1", "Arabidopsis thaliana",
            "receptor kinase", "120", "a" * 64, "local_record_cache", "file://fixture",
            "2026-08-12T00:00:00+00:00", "PMID:8402901", "none", "frozen", "verified",
        ]) + "\n", encoding="utf-8")
        return path

    def test_converts_to_validated_blast_panel_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            panel = self._panel_tsv(root)
            sequence = ("CDLYQWLICSGRWNGGYHCLKFRVGRYVDSYIMRQYFNMFLEVAWCGNRTDNRMNHMNNQDKTTEQPNHVN"
                        "KQKNLGSRNSEYNGGSQYKGNLRRCKLTIYEEPFSWSDAPQCWSCIFET")
            sequences = root / "verified_sequences.json"
            sequences.write_text(json.dumps({"hyp_001": sequence}), encoding="utf-8")
            frame = panel_from_verified_bait_panel(str(panel), str(sequences))
            self.assertEqual(["P30924"], list(frame["panel_id"]))
            self.assertEqual(PANEL_COLUMNS, list(frame.columns))
            self.assertEqual("none", frame.iloc[0]["rank_contribution"])
            for record in frame.to_dict("records"):
                verdict = validate_panel_row({column: str(record.get(column, "")) for column in PANEL_COLUMNS})
                self.assertEqual("verified", verdict["review_status"], verdict["review_reason"])

    def test_empty_panel_or_sequences_yields_no_verified_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            panel = self._panel_tsv(root)
            empty_sequences = root / "empty_sequences.json"
            empty_sequences.write_text("{}", encoding="utf-8")
            frame = panel_from_verified_bait_panel(str(panel), str(empty_sequences))
            self.assertEqual(0, len(frame))
            # A panel with no rows still converts (header only): freeze reports
            # review_required_small_panel, never a failure.
            empty_panel = root / "empty_panel.tsv"
            empty_panel.write_text("\t".join([
                "hypothesis_id", "bait_id", "sequence_sha256", "reference_gene_name",
                "reference_species", "pathway_or_process", "gene_family", "source_provider",
                "source_url", "retrieved_at",
            ]) + "\n", encoding="utf-8")
            frame = panel_from_verified_bait_panel(str(empty_panel), str(empty_sequences))
            self.assertEqual(0, len(frame))

    def test_missing_required_column_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            bad = root / "bad_panel.tsv"
            bad.write_text("hypothesis_id\tbait_id\nh1\tACC1\n", encoding="utf-8")
            with self.assertRaises(ValueError) as caught:
                panel_from_verified_bait_panel(str(bad), str(root / "x.json"))
            self.assertIn("sequence_sha256", str(caught.exception))

    def test_upstream_review_required_row_is_never_promoted(self) -> None:
        """Issue #374 review: the converter inherits the upstream per-row
        review_status. A bait the upstream freeze marked review_required for a
        bait-level reason stays review_required -- never VERIFIED -- so the
        downstream freeze reports review_required_small_panel and the search
        skips the row."""
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            panel = root / "small_panel.tsv"
            panel.write_text("\t".join([
                "hypothesis_id", "bait_id", "pathway_or_process", "gene_family", "reference_gene_name",
                "reference_species", "expected_role", "sequence_length", "sequence_sha256",
                "source_provider", "source_url", "retrieved_at", "function_source_ids",
                "rank_contribution", "panel_status", "review_status",
            ]) + "\n" + "\t".join([
                "hyp_001", "P30924", "meristem maintenance", "CLV", "CLV1", "Arabidopsis thaliana",
                "receptor kinase", "120", "a" * 64, "local_record_cache", "file://fixture",
                "2026-08-12T00:00:00+00:00", "PMID:8402901", "none", "review_required_small_panel",
                "review_required",
            ]) + "\n", encoding="utf-8")
            sequence = ("CDLYQWLICSGRWNGGYHCLKFRVGRYVDSYIMRQYFNMFLEVAWCGNRTDNRMNHMNNQDKTTEQPNHVN"
                        "KQKNLGSRNSEYNGGSQYKGNLRRCKLTIYEEPFSWSDAPQCWSCIFET")
            sequences = root / "verified_sequences.json"
            sequences.write_text(json.dumps({"hyp_001": sequence}), encoding="utf-8")
            frame = panel_from_verified_bait_panel(str(panel), str(sequences))
            self.assertEqual(["P30924"], list(frame["panel_id"]))
            self.assertEqual("review_required", frame.iloc[0]["review_status"])
            verdict = validate_panel_row({column: str(frame.iloc[0].get(column, "")) for column in PANEL_COLUMNS})
            self.assertEqual("review_required", verdict["review_status"], verdict["review_reason"])

    def test_duplicate_accession_with_identical_sequence_is_frozen_once(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            panel = self._panel_tsv(root)
            lines = panel.read_text(encoding="utf-8").splitlines()
            lines.append(lines[1].replace("hyp_001", "hyp_002", 1))
            panel.write_text("\n".join(lines) + "\n", encoding="utf-8")
            sequence = ("CDLYQWLICSGRWNGGYHCLKFRVGRYVDSYIMRQYFNMFLEVAWCGNRTDNRMNHMNNQDKTTEQPNHVN"
                        "KQKNLGSRNSEYNGGSQYKGNLRRCKLTIYEEPFSWSDAPQCWSCIFET")
            sequences = root / "verified_sequences.json"
            sequences.write_text(
                json.dumps({"hyp_001": sequence, "hyp_002": sequence}), encoding="utf-8"
            )

            frame = panel_from_verified_bait_panel(str(panel), str(sequences))
            self.assertEqual(["P30924"], list(frame["panel_id"]))
            converted = root / "blast_panel.tsv"
            frame.to_csv(converted, sep="\t", index=False, lineterminator="\n", encoding="utf-8")
            self.assertEqual(
                0,
                main([
                    "blast-panel-freeze", "--panel", str(converted),
                    "--out-dir", str(root / "freeze"), "--case-id", "duplicate-accession",
                ]),
            )
            manifest = json.loads(
                (root / "freeze" / "blast_panel_freeze_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(1, manifest["verified_count"])

    def test_duplicate_accession_with_different_sequences_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            panel = self._panel_tsv(root)
            lines = panel.read_text(encoding="utf-8").splitlines()
            lines.append(lines[1].replace("hyp_001", "hyp_002", 1))
            panel.write_text("\n".join(lines) + "\n", encoding="utf-8")
            sequences = root / "verified_sequences.json"
            sequences.write_text(
                json.dumps({"hyp_001": "ACDEFG", "hyp_002": "ACDEFA"}), encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "one bait_id to different sequences"):
                panel_from_verified_bait_panel(str(panel), str(sequences))


@unittest.skipUnless(shutil.which("blastp") and shutil.which("makeblastdb"), "BLAST+ (blastp/makeblastdb) not installed")
class BlastPanelIntegrationTests(unittest.TestCase):
    def _freeze(self, root: Path, panel: Path) -> dict[str, object]:
        self.assertEqual(
            0,
            main(["blast-panel-freeze", "--panel", str(panel), "--out-dir", str(root / "freeze"), "--case-id", "integration"]),
        )
        return json.loads((root / "freeze" / "blast_panel_freeze_manifest.json").read_text(encoding="utf-8"))

    def _search(self, root: Path) -> dict[str, object]:
        self.assertEqual(
            0,
            main(
                [
                    "blast-panel-search",
                    "--freeze", str(root / "freeze" / "blast_panel_freeze_manifest.json"),
                    "--proteins", str(root / "fixtures" / "interval_proteins.fa"),
                    "--out-dir", str(root / "search"),
                    "--threads", "2",
                ]
            ),
        )
        return json.loads((root / "search" / "blast_panel_search_manifest.json").read_text(encoding="utf-8"))

    def _evidence(self, root: Path) -> pd.DataFrame:
        return pd.read_csv(root / "search" / "blast_evidence.tsv", sep="\t", dtype=str).fillna("")

    def test_panel_search_e2e(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            shutil.copytree(FIXTURE, root / "fixtures")

            freeze = self._freeze(root, root / "fixtures" / "panel_ok.tsv")
            self.assertEqual("verified", freeze["panel_status"])
            self.assertEqual(3, freeze["verified_count"])

            search = self._search(root)
            self.assertEqual("evaluable", search["provider_status"])
            self.assertEqual("post_freeze_similarity_only", search["interval_usage"])
            self.assertEqual(3, search["counts"]["rows_evaluated"])
            self.assertEqual(1, search["counts"]["unique"])
            self.assertEqual(1, search["counts"]["ambiguous_paralog"])
            self.assertEqual(1, search["counts"]["no_threshold_hit"])
            self.assertEqual(str(root / "search" / "eligible_interval_hits.tsv"), search["outputs"]["eligible_interval_hits"])
            self.assertTrue((root / "search" / "eligible_interval_hits.tsv").exists())
            eligible = pd.read_csv(root / "search" / "eligible_interval_hits.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(search["counts"]["eligible_interval_hits"], len(eligible))
            self.assertTrue(all(eligible["eligibility"] == "eligible"))
            self.assertEqual(["bait4"], search["skipped_review_required"])

            evidence = self._evidence(root)
            self.assertEqual(["bait1", "bait2", "bait3"], list(evidence["panel_id"]))
            by_id = {str(row["panel_id"]): row for row in evidence.to_dict("records")}
            self.assertTrue(all(str(row["rank_contribution"]) == "none" for row in by_id.values()))

            bait1 = by_id["bait1"]
            self.assertEqual("verified", bait1["review_status"])
            self.assertEqual("evaluable", bait1["provider_status"])
            self.assertEqual("unique", bait1["ambiguity_status"])
            self.assertEqual("g1", bait1["top_hit_target_id"])
            self.assertGreaterEqual(float(bait1["top_hit_identity"]), 99.0)
            self.assertEqual("3", bait1["hit_count"])
            self.assertEqual("2", bait1["paralog_count"])

            bait2 = by_id["bait2"]
            self.assertEqual("ambiguous_paralog", bait2["ambiguity_status"])
            self.assertEqual("2", bait2["hit_count"])
            self.assertEqual("1", bait2["paralog_count"])

            bait3 = by_id["bait3"]
            self.assertEqual("no_threshold_hit", bait3["ambiguity_status"])
            self.assertEqual("0", bait3["hit_count"])
            self.assertEqual("", bait3["top_hit_target_id"])

            self.assertNotIn("bait4", by_id)
            self.assertEqual(["bait4"], search["skipped_review_required"])
            self.assertTrue((root / "search" / "queries.fa").exists())
            self.assertTrue((root / "search" / "panel_vs_interval.tsv").exists())

    def test_verified_single_bait_is_evaluated_with_panel_review_warning(self) -> None:
        """Issue #538: panel diversity and per-bait verification are separate."""
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            shutil.copytree(FIXTURE, root / "fixtures")
            source = pd.read_csv(root / "fixtures" / "panel_ok.tsv", sep="\t", dtype=str).fillna("")
            single = source[source["review_status"] == "verified"].head(1)
            panel = root / "fixtures" / "panel_single_verified.tsv"
            single.to_csv(panel, sep="\t", index=False, lineterminator="\n", encoding="utf-8")

            freeze = self._freeze(root, panel)
            self.assertEqual(1, freeze["verified_count"])
            self.assertEqual("review_required_small_panel", freeze["panel_status"])

            search = self._search(root)
            self.assertEqual("evaluable", search["provider_status"])
            self.assertEqual("review_required_small_panel", search["freeze_panel_status"])
            self.assertEqual(1, search["counts"]["rows_evaluated"])
            self.assertEqual([], search["skipped_review_required"])
            evidence = self._evidence(root)
            self.assertEqual(1, len(evidence))
            self.assertEqual("verified", evidence.iloc[0]["review_status"])
            self.assertEqual("none", evidence.iloc[0]["rank_contribution"])

    def test_safe_abstention_all_review_required(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            shutil.copytree(FIXTURE, root / "fixtures")
            panel = root / "fixtures" / "panel_all_review_required.tsv"
            text = (root / "fixtures" / "panel_ok.tsv").read_text(encoding="utf-8")
            for accession in INVALID_ACCESSIONS:
                text = text.replace(accession, "P0FAKE")
            panel.write_text(text, encoding="utf-8")

            freeze = self._freeze(root, panel)
            self.assertEqual(0, freeze["verified_count"])
            self.assertEqual("review_required_small_panel", freeze["panel_status"])

            search = self._search(root)
            self.assertEqual("not_evaluated_no_verified_bait", search["provider_status"])
            self.assertEqual(0, search["counts"]["rows_evaluated"])

            evidence = self._evidence(root)
            self.assertEqual(0, len(evidence))
            eligible = pd.read_csv(root / "search" / "eligible_interval_hits.tsv", sep="\t", dtype=str)
            self.assertEqual(0, len(eligible))
            self.assertTrue((root / "search" / "blast_panel_search_manifest.json").exists())

    def test_post_freeze_mutation_is_refused_by_verify_and_search(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            shutil.copytree(FIXTURE, root / "fixtures")
            self._freeze(root, root / "fixtures" / "panel_ok.tsv")

            panel_path = root / "freeze" / "frozen_bait_panel.tsv"
            panel_path.write_text(panel_path.read_text(encoding="utf-8").replace("CDLYQ", "CDLYR", 1), encoding="utf-8")
            manifest = root / "freeze" / "blast_panel_freeze_manifest.json"

            with self.assertRaises(ValueError):
                main(["blast-panel-verify", "--freeze", str(manifest)])
            with self.assertRaises(ValueError):
                main(
                    [
                        "blast-panel-search",
                        "--freeze", str(manifest),
                        "--proteins", str(root / "fixtures" / "interval_proteins.fa"),
                        "--out-dir", str(root / "search_mutated"),
                    ]
                )
            self.assertFalse((root / "search_mutated" / "queries.fa").exists())


if __name__ == "__main__":
    unittest.main()
