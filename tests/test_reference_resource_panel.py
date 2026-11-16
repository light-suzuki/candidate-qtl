from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import main

HEADER = "resource_id\trole\tresource_type\ttaxon\tsource\turi_or_path\texternal_api\n"


class ReferenceResourcePanelTests(unittest.TestCase):
    def _write_panel(self, root: Path, body: str, *, with_native_files: bool = True) -> Path:
        if with_native_files:
            (root / "native.fa").write_text(">chr1\nACGT\n", encoding="utf-8")
            (root / "native.gff3").write_text("##gff-version 3\n", encoding="utf-8")
        panel = root / "panel.tsv"
        panel.write_text(HEADER + body, encoding="utf-8")
        return panel

    def test_native_only_panel_is_local_only_and_completes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            panel = self._write_panel(
                root,
                "native_genome\tnative\tgenome_fasta\ttxid1\thttps://example.org/a\tnative.fa\tfalse\n"
                "native_annotation\tnative\tannotation_gff3\ttxid1\thttps://example.org/a\tnative.gff3\tfalse\n",
            )
            out = root / "out"
            self.assertEqual(0, main(["reference-panel", "--panel", str(panel), "--out-dir", str(out)]))
            manifest = json.loads((out / "reference_resource_panel_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(manifest["native_available"])
            self.assertTrue(manifest["local_only"])
            self.assertFalse(manifest["external_api_used"])
            self.assertEqual(sorted(["close_anchor", "functional_reference", "outgroup"]), manifest["optional_roles_absent"])

    def test_missing_optional_resource_records_machine_readable_skip_reason(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            panel = self._write_panel(
                root,
                "native_genome\tnative\tgenome_fasta\ttxid1\thttps://example.org/a\tnative.fa\tfalse\n"
                "anchor\tclose_anchor\tgenome_fasta\ttxid2\thttps://example.org/b\tmissing_anchor.fa\tfalse\n",
            )
            out = root / "out"
            self.assertEqual(0, main(["reference-panel", "--panel", str(panel), "--out-dir", str(out)]))
            rows = pd.read_csv(out / "reference_resource_panel.tsv", sep="\t", dtype=str).fillna("").set_index("resource_id")
            self.assertEqual("skipped", rows.loc["anchor", "availability"])
            self.assertEqual("local_resource_missing", rows.loc["anchor", "skip_reason"])
            manifest = json.loads((out / "reference_resource_panel_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("local_resource_missing", manifest["skip_reasons"]["anchor"])
            self.assertTrue(manifest["native_available"])

    def test_external_resource_is_declared_not_fetched_and_flips_local_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            panel = self._write_panel(
                root,
                "native_genome\tnative\tgenome_fasta\ttxid1\thttps://example.org/a\tnative.fa\tfalse\n"
                "ref_api\tfunctional_reference\trest_api\ttxid3\thttps://rest.uniprot.org\thttps://rest.uniprot.org/uniprotkb\ttrue\n",
            )
            out = root / "out"
            self.assertEqual(0, main(["reference-panel", "--panel", str(panel), "--out-dir", str(out)]))
            rows = pd.read_csv(out / "reference_resource_panel.tsv", sep="\t", dtype=str).fillna("").set_index("resource_id")
            self.assertEqual("declared_external", rows.loc["ref_api", "availability"])
            self.assertEqual("external_resource_not_fetched_local_first", rows.loc["ref_api", "skip_reason"])
            manifest = json.loads((out / "reference_resource_panel_manifest.json").read_text(encoding="utf-8"))
            self.assertTrue(manifest["external_api_used"])
            self.assertFalse(manifest["local_only"])
            self.assertIn("https://rest.uniprot.org/uniprotkb", manifest["external_destinations"])

    def test_missing_native_resource_is_a_hard_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            panel = self._write_panel(
                root,
                "native_genome\tnative\tgenome_fasta\ttxid1\thttps://example.org/a\tabsent.fa\tfalse\n",
                with_native_files=False,
            )
            with self.assertRaises(ValueError):
                main(["reference-panel", "--panel", str(panel), "--out-dir", str(root / "out")])

    def test_unknown_role_and_duplicate_id_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bad_role = self._write_panel(
                root, "native_genome\tsister_taxon\tgenome_fasta\ttxid1\thttps://example.org/a\tnative.fa\tfalse\n"
            )
            with self.assertRaises(ValueError):
                main(["reference-panel", "--panel", str(bad_role), "--out-dir", str(root / "o1")])
            dup = self._write_panel(
                root,
                "native_genome\tnative\tgenome_fasta\ttxid1\thttps://example.org/a\tnative.fa\tfalse\n"
                "native_genome\tnative\tannotation_gff3\ttxid1\thttps://example.org/a\tnative.gff3\tfalse\n",
            )
            with self.assertRaises(ValueError):
                main(["reference-panel", "--panel", str(dup), "--out-dir", str(root / "o2")])


if __name__ == "__main__":
    unittest.main()
