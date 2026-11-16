"""The default weight registry must ship inside the wheel and stay in sync.

A clean ``pip install`` has no source tree, so ranking would fail if the default
registry were only reachable at the repo-root ``templates/`` path. A vendored
copy under ``candidate_qtl/resources/`` is packaged; these tests keep it
byte-identical to the authoritative repo-root copy and confirm it resolves.
"""

from __future__ import annotations

import unittest
from pathlib import Path

from candidate_qtl.ranking_policy import default_registry_path, load_weight_registry


class WeightRegistryPackagingTests(unittest.TestCase):
    def test_default_registry_path_exists_and_loads(self) -> None:
        path = default_registry_path()
        self.assertTrue(path.exists(), f"default registry not found at {path}")
        registry, meta = load_weight_registry()
        self.assertGreater(len(registry), 0)
        self.assertIn("sha256", meta)

    def test_vendored_copy_is_byte_identical_to_source_tree(self) -> None:
        repo_root = Path(__file__).resolve().parents[1]
        source_tree = repo_root / "templates" / "evidence_weight_registry.tsv"
        vendored = repo_root / "src" / "candidate_qtl" / "resources" / "evidence_weight_registry.tsv"
        self.assertTrue(vendored.exists(), "vendored registry copy is missing from the package")
        self.assertTrue(source_tree.exists(), "authoritative repo-root registry is missing")
        self.assertEqual(
            source_tree.read_bytes(), vendored.read_bytes(),
            "templates/evidence_weight_registry.tsv and candidate_qtl/data/ copy have drifted; re-sync them",
        )


if __name__ == "__main__":
    unittest.main()
