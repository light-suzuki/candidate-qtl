from __future__ import annotations

import json
import os
import tempfile
import unittest
from pathlib import Path

from candidate_qtl.cli import main
from candidate_qtl.core import file_fingerprint, fingerprint_is_weak


REPO = Path(__file__).resolve().parents[1]
TUTORIAL = REPO / "docs" / "tutorial"


class FingerprintStrengthTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def test_small_file_is_strong_full_sha256(self) -> None:
        path = self.root / "small.txt"
        path.write_text("hello", encoding="utf-8")
        fp = file_fingerprint(path)
        self.assertEqual(fp["method"], "full_sha256")
        self.assertEqual(fp["strength"], "strong")
        self.assertIsNotNone(fp["sha256"])
        self.assertFalse(fingerprint_is_weak(fp))

    def test_metadata_only_is_weak(self) -> None:
        path = self.root / "big.fa"
        path.write_text(">c\nACGT\n", encoding="utf-8")
        fp = file_fingerprint(path, full_hash_max_bytes=0)
        self.assertEqual(fp["method"], "metadata_only")
        self.assertEqual(fp["strength"], "weak")
        self.assertIsNone(fp["sha256"])
        self.assertTrue(fingerprint_is_weak(fp))

    def test_strict_full_hash_detects_same_size_different_content(self) -> None:
        a = self.root / "a.fa"
        b = self.root / "b.fa"
        a.write_text(">x\nAAAA\n", encoding="utf-8")
        b.write_text(">x\nGGGG\n", encoding="utf-8")  # same byte length, different content
        # Pin identical mtimes so a metadata-only fingerprint cannot tell them apart.
        os.utime(a, ns=(1_000, 1_000))
        os.utime(b, ns=(1_000, 1_000))
        weak_a = file_fingerprint(a, full_hash_max_bytes=0)
        weak_b = file_fingerprint(b, full_hash_max_bytes=0)
        self.assertEqual(weak_a["size_bytes"], weak_b["size_bytes"])
        self.assertEqual(weak_a["mtime_ns"], weak_b["mtime_ns"])
        # Full content hashing distinguishes them; the weak fingerprints cannot.
        self.assertNotEqual(file_fingerprint(a, full_hash=True)["sha256"], file_fingerprint(b, full_hash=True)["sha256"])

    def test_extract_manifest_reports_fingerprint_strength(self) -> None:
        out = self.root / "extract"
        self.assertEqual(
            0,
            main(["extract", "--fasta", str(TUTORIAL / "genome.fa"), "--gff3", str(TUTORIAL / "genes.gff3"), "--interval", "chr1:1-900", "--out-dir", str(out)]),
        )
        manifest = json.loads((out / "extract_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["fasta_fingerprint"]["strength"], "strong")
        self.assertIn("weak_input_fingerprints", manifest)
        self.assertEqual(manifest["weak_input_fingerprints"], [])


if __name__ == "__main__":
    unittest.main()
