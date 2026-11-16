"""C rescue survives two real boundary-consumer generations without crowding AB."""
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.boundary_tier_review import freeze_boundary_scope, _load_scope
from tests import test_boundary_tier_review as support
from tests.test_boundary_tier_review import inventory, sha, write


class ProtectedContinuationTests(unittest.TestCase):
    review = support.ScopeTests.review
    consume = support.ScopeTests.consume

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.source = self.root / "inventory.tsv"
        rows = inventory(3, 8, 35, 30, 30)
        self.protected = {row["candidate_gene_id"] for row in rows
                          if row["candidate_gene_id"].startswith("t4_")}
        write(self.source, rows)
        self.source_hash = sha(self.source)
        self.scope = freeze_boundary_scope(
            str(self.source), str(self.root / "scope"),
            inventory_sha256=self.source_hash, case_id="test",
            protected_gene_ids=self.protected,
        )

    def tearDown(self):
        self.temp.cleanup()

    def test_c_remains_through_two_real_selections_and_next_scope(self):
        proposal, review = self.review(keep=0)
        first = self.consume(proposal, review)
        next_scope = json.loads(Path(first["next_scope_manifest"]["path"]).read_text())
        next_scope.update(manifest=first["next_scope_manifest"]["path"],
                          manifest_sha256=first["next_scope_manifest"]["sha256"])
        self.assertEqual(self.protected, set(next_scope["protected_gene_ids"]))
        _, _, plan = _load_scope(next_scope["manifest"], next_scope["manifest_sha256"], "test")
        self.assertTrue(self.protected <= set(plan["frozen_gene_ids"]))
        self.assertFalse(self.protected & {row["candidate_gene_id"] for row in plan["boundary"]})
        proposal, review = self.review(keep=9, scope=next_scope, name="second")
        second = self.consume(proposal, review, scope=next_scope, name="second-selected")
        wet = pd.read_csv(second["wet_candidates"]["path"], sep="\t", dtype=str)
        ids = set(wet.candidate_gene_id)
        self.assertTrue(self.protected <= ids)
        self.assertEqual(50, len(ids))  # 20 ordinary selections plus all 30 C copies.
        self.assertEqual(50, len(wet))
        self.assertEqual(self.source_hash, sha(self.source))

    def test_boundary_cannot_defer_protected_c(self):
        proposal, review = self.review(keep=0)
        rows = pd.read_csv(proposal, sep="\t").to_dict("records")
        rows[0]["candidate_gene_id"] = sorted(self.protected)[0]
        write(Path(proposal), rows)
        with self.assertRaisesRegex(ValueError, "higher/lower"):
            self.consume(proposal, review)


if __name__ == "__main__":
    unittest.main()
