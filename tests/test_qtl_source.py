from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from candidate_qtl.cli import main


class QtlSourceTests(unittest.TestCase):
    def test_routes_physical_interval_without_markers_to_qtlift(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "case.yaml"
            manifest.write_text(
                "case_id: fas\nsource:\n  source_id: PMID:1\n  source_reference_name: v1\n"
                "  source_contig: chr5\n  source_start: 10\n  source_end: 20\n"
                "  coordinate_basis: physical\n  literature_search_status: complete\n",
                encoding="utf-8",
            )
            self.assertEqual(0, main(["qtl-source-init", "--manifest", str(manifest), "--out-dir", str(root / "out")]))
            result = json.loads((root / "out" / "qtl_source_preflight.json").read_text(encoding="utf-8"))
            self.assertEqual("qtlift_interval", result["route"])
            self.assertTrue(result["positional_ranking_allowed"])

    def test_unresolved_is_rescue_only_and_emits_review(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "case.yaml"
            manifest.write_text("case_id: sin\nsource:\n  literature_search_status: partial\n", encoding="utf-8")
            self.assertEqual(0, main(["qtl-source-init", "--manifest", str(manifest), "--out-dir", str(root / "out")]))
            result = json.loads((root / "out" / "qtl_source_preflight.json").read_text(encoding="utf-8"))
            self.assertEqual("unresolved_position", result["route"])
            self.assertEqual("rescue", result["candidate_scope"])
            self.assertIn("complete source interval", (root / "out" / "review_required.tsv").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
