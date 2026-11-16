"""Regression coverage for triple-route reusable family-index auto-wiring (#551)."""

from __future__ import annotations

import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import yaml

from candidate_qtl.cli import main
from candidate_qtl.family_index_resolution import FAMILY_INDEX_ENV
from candidate_qtl.family_index_resolution import _default_candidates
from candidate_qtl.family_index_resolution import resolve_family_index


class FamilyIndexAutowiringTests(unittest.TestCase):
    @staticmethod
    def _write_family_index(path: Path, *, scope: str = "reference_genome") -> None:
        path.write_text(
            "gene_id\tassignment_status\tbest_subject\tfamily_stems\tindex_scope\n"
            f"g1\tassigned\tREF1\tFAM1\t{scope}\n",
            encoding="utf-8",
        )
        (path.parent / "family_index_manifest.json").write_text(json.dumps({
            "schema_version": "1.4",
            "index_scope": scope,
            "reference_label": "test_reference",
            "reference_species": "Test species",
        }) + "\n", encoding="utf-8")

    def _dry_run(self, root: Path, *extra: str) -> dict[str, object]:
        out = root / "plan"
        self.assertEqual(0, main([
            "guided-profile", "--profile", "triple-route", "--out-dir", str(out), *extra,
        ]))
        return json.loads((out / "guided_profile_manifest.json").read_text(encoding="utf-8"))

    @staticmethod
    def _family_index_arguments(manifest: dict[str, object]) -> dict[str, object]:
        stages = manifest["stages"]
        return {
            str(stage["command"]): stage["arguments"].get("family_index")
            for stage in stages
            if stage["command"] in {"route-a-interval-match", "route-b-prepare"}
        }

    def test_resources_registry_relative_family_index_is_wired_to_both_routes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            family_index = root / "family_index.tsv"
            self._write_family_index(family_index)
            resources = root / "resources.yaml"
            resources.write_text(yaml.safe_dump({
                "family_index": {"type": "family_index_tsv", "path": "family_index.tsv"},
            }, sort_keys=False), encoding="utf-8")

            with patch.dict(os.environ, {FAMILY_INDEX_ENV: ""}, clear=False), \
                    patch("candidate_qtl.family_index_resolution._default_candidates", return_value=[]):
                manifest = self._dry_run(root, "--resources", str(resources))

            self.assertTrue(manifest["family_index_auto_resolve_attempted"])
            resolution = manifest["family_index_resolution"]
            self.assertEqual("resolved", resolution["status"])
            self.assertEqual("resources", resolution["source"])
            self.assertEqual(str(family_index), resolution["path"])
            self.assertEqual({
                "route-a-interval-match": str(family_index),
                "route-b-prepare": str(family_index),
            }, self._family_index_arguments(manifest))

    def test_nested_resources_registry_is_supported(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            family_index = root / "index" / "family_index.tsv"
            family_index.parent.mkdir()
            self._write_family_index(family_index)
            resources = root / "resources.yaml"
            resources.write_text(yaml.safe_dump({
                "resources": {"family_index": {"path": "index/family_index.tsv"}},
            }, sort_keys=False), encoding="utf-8")

            with patch.dict(os.environ, {FAMILY_INDEX_ENV: ""}, clear=False), \
                    patch("candidate_qtl.family_index_resolution._default_candidates", return_value=[]):
                manifest = self._dry_run(root, "--resources", str(resources))

            self.assertEqual("resources", manifest["family_index_resolution"]["source"])
            self.assertEqual(str(family_index), manifest["family_index_resolution"]["path"])

    def test_environment_family_index_is_wired_when_registry_has_none(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            family_index = root / "env_family_index.tsv"
            self._write_family_index(family_index)

            with patch.dict(os.environ, {FAMILY_INDEX_ENV: str(family_index)}, clear=False), \
                    patch("candidate_qtl.family_index_resolution._default_candidates", return_value=[]):
                manifest = self._dry_run(root)

            self.assertTrue(manifest["family_index_auto_resolve_attempted"])
            self.assertEqual("environment", manifest["family_index_resolution"]["source"])
            self.assertEqual({
                "route-a-interval-match": str(family_index),
                "route-b-prepare": str(family_index),
            }, self._family_index_arguments(manifest))

    def test_missing_auto_sources_are_recorded_without_fabricating_an_index(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            missing = root / "missing_family_index.tsv"

            with patch.dict(os.environ, {FAMILY_INDEX_ENV: ""}, clear=False), \
                    patch("candidate_qtl.family_index_resolution._default_candidates", return_value=[missing]):
                manifest = self._dry_run(root)

            self.assertTrue(manifest["family_index_auto_resolve_attempted"])
            resolution = manifest["family_index_resolution"]
            self.assertEqual("not_found", resolution["status"])
            self.assertIsNone(resolution["path"])
            self.assertEqual({
                "route-a-interval-match": None,
                "route-b-prepare": None,
            }, self._family_index_arguments(manifest))
            self.assertTrue(any(
                attempt["source"] == "default_cache" and attempt["status"] == "file_missing"
                for attempt in resolution["attempts"]
            ))

    def test_explicit_cli_family_index_keeps_precedence_and_skips_auto_resolution(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            explicit = root / "explicit.tsv"
            explicit.write_text("gene_id\tfamily_id\n", encoding="utf-8")
            env_index = root / "env.tsv"
            env_index.write_text("gene_id\tfamily_id\n", encoding="utf-8")

            with patch.dict(os.environ, {FAMILY_INDEX_ENV: str(env_index)}, clear=False):
                manifest = self._dry_run(root, "--family-index", str(explicit))

            self.assertFalse(manifest["family_index_auto_resolve_attempted"])
            resolution = manifest["family_index_resolution"]
            self.assertEqual("explicit", resolution["status"])
            self.assertEqual("explicit_cli", resolution["source"])
            self.assertEqual({
                "route-a-interval-match": str(explicit),
                "route-b-prepare": str(explicit),
            }, self._family_index_arguments(manifest))

    def test_interval_scoped_auto_candidate_is_rejected_for_both_routes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            interval_index = root / "family_index.tsv"
            self._write_family_index(interval_index, scope="interval_candidates")

            with patch.dict(os.environ, {FAMILY_INDEX_ENV: str(interval_index)}, clear=False), \
                    patch("candidate_qtl.family_index_resolution._default_candidates", return_value=[]):
                manifest = self._dry_run(root)

            resolution = manifest["family_index_resolution"]
            self.assertEqual("not_found", resolution["status"])
            env_attempt = next(attempt for attempt in resolution["attempts"]
                               if attempt["source"] == "environment")
            self.assertEqual("invalid_reusable_index", env_attempt["status"])
            self.assertEqual("not_reusable_reference_genome_index", env_attempt["validation"]["reason"])
            self.assertEqual({
                "route-a-interval-match": None,
                "route-b-prepare": None,
            }, self._family_index_arguments(manifest))

    def test_default_candidates_are_only_conventional_cache_locations(self) -> None:
        """Default discovery must be portable, not tied to a maintainer home layout."""
        candidates = _default_candidates()
        self.assertEqual(2, len(candidates))
        self.assertTrue(str(candidates[0]).replace("\\", "/").endswith("candidate-qtl/family_index.tsv"))
        self.assertTrue(str(candidates[1]).replace("\\", "/").endswith(".candidate-qtl/family_index.tsv"))
        self.assertFalse(any("famidx_v2" in str(candidate) for candidate in candidates))

    def test_patched_default_candidate_resolves_when_valid(self) -> None:
        """A valid nonstandard location still works when supplied by the resolver seam."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            family_index = root / "custom-cache" / "family_index.tsv"
            family_index.parent.mkdir(parents=True)
            self._write_family_index(family_index)

            with patch.dict(os.environ, {FAMILY_INDEX_ENV: ""}, clear=False), \
                    patch("candidate_qtl.family_index_resolution._default_candidates",
                          return_value=[family_index]):
                path, resolution = resolve_family_index(None, None)

            self.assertEqual("resolved", resolution["status"])
            self.assertEqual("default_cache", resolution["source"])
            self.assertEqual(str(family_index), path)


if __name__ == "__main__":
    unittest.main()
