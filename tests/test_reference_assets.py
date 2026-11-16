from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import json
import os
from pathlib import Path
import tempfile
import threading
import unittest
from unittest.mock import patch

from candidate_qtl import reference_assets as assets
from candidate_qtl.reference_extract import extract_candidates


class ReferenceAssetTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.source = self.root / "original"
        self.source.mkdir()
        self.fasta = self.source / "genome.fa"
        self.gff = self.source / "genes.gff3"
        sequence = list("N" * 90)
        sequence[0:9] = "ATGAAATAA"
        sequence[19:28] = "TTATTTCAT"
        sequence[59:68] = "ATGAAATAA"
        self.fasta.write_text(">1\n" + "".join(sequence) + "\n", encoding="utf-8")
        lines = ["##gff-version 3"]
        for gene, start, end, strand in (("Plus", 1, 9, "+"), ("Minus", 20, 28, "-"), ("Outside", 60, 68, "+")):
            lines.extend((
                f"chr1\ttest\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gene}",
                f"chr1\ttest\tmRNA\t{start}\t{end}\t.\t{strand}\t.\tID={gene}.t1;Parent={gene}",
                f"chr1\ttest\tCDS\t{start}\t{end}\t.\t{strand}\t0\tParent={gene}.t1",
            ))
        lines.extend(("chr1\ttest\tgene\t35\t39\t.\t+\t.\tID=Dark",
                      "chr1\ttest\tpseudogene\t40\t48\t.\t+\t.\tID=Pseudo"))
        self.gff.write_text("\n".join(lines) + "\n", encoding="utf-8")
        self.cache = self.root / "cache"

    def prepare(self, cache=None):
        return assets.prepare_reference_asset(str(self.fasta), str(self.gff), str(cache or self.cache),
                                              reference_version="synthetic-v1")

    def extract(self, asset, name, **kwargs):
        return assets.extract_reference_interval(asset, "chr1:5-48", str(self.root / name),
                                                 seqid_map_arg="chr1=1", **kwargs)

    def test_registration_copies_once_builds_only_managed_index_and_hashes_every_time(self):
        import pysam
        original_index = Path(str(self.fasta) + ".fai")
        original_index.write_text("old invalid index retained untouched\n", encoding="utf-8")
        originals = (self.fasta, self.gff, original_index)
        before = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in originals}
        with patch.object(assets, "_copy", wraps=assets._copy) as copied, \
                patch.object(pysam, "faidx", wraps=pysam.faidx) as indexed, \
                patch.object(assets, "sha256", wraps=assets.sha256) as hashed:
            first = self.prepare()
            second = self.prepare()
        self.assertFalse(first["cache_hit"])
        self.assertTrue(second["cache_hit"])
        self.assertEqual(first["manifest_sha256"], second["manifest_sha256"])
        self.assertEqual(2, copied.call_count)
        self.assertEqual(1, indexed.call_count)
        for path in (self.fasta, self.gff):
            self.assertEqual(2, sum(call.args[0] == path for call in hashed.call_args_list))
        self.assertEqual(before, {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in originals})
        self.assertNotEqual(self.fasta.stat().st_ino, Path(first["fasta"]).stat().st_ino)
        self.assertEqual(0, Path(first["fasta"]).stat().st_mode & 0o222)
        self.assertEqual("strong", first["source_fingerprints"]["fasta"]["strength"])

    def test_missing_original_index_is_never_created(self):
        asset = self.prepare()
        self.assertFalse(Path(str(self.fasta) + ".fai").exists())
        self.assertTrue(Path(asset["fasta"] + ".fai").is_file())

    def test_gzip_annotation_keeps_decoding_suffix_and_original_bytes(self):
        compressed = self.source / "genes.gff3.gz"
        with gzip.open(compressed, "wb") as stream:
            stream.write(self.gff.read_bytes())
        before = compressed.read_bytes(), compressed.stat().st_mtime_ns
        asset = assets.prepare_reference_asset(str(self.fasta), str(compressed), str(self.cache), reference_version="v1")
        self.assertTrue(asset["gff3"].endswith(".gff3.gz"))
        extracted = self.extract(asset, "compressed")
        self.assertEqual(4, extracted["candidate_count"])
        self.assertEqual(before, (compressed.read_bytes(), compressed.stat().st_mtime_ns))

    def test_compressed_fasta_is_rejected_without_creating_an_index(self):
        compressed = self.source / "genome.fa.gz"
        with gzip.open(compressed, "wb") as stream:
            stream.write(self.fasta.read_bytes())
        with self.assertRaisesRegex(ValueError, "uncompressed FASTA"):
            assets.prepare_reference_asset(str(compressed), str(self.gff), str(self.cache), reference_version="v1")
        self.assertFalse(Path(str(compressed) + ".fai").exists())
        self.assertFalse(self.cache.exists())

    def test_interval_hit_matches_fresh_science_and_relocates_case_manifest(self):
        asset = self.prepare()
        for name in ("first", "second"):
            out = self.root / name
            out.mkdir()
            (out / "phenotype.json").write_text(json.dumps({"case_id": name, "phenotype": name}), encoding="utf-8")
        with patch.object(assets, "extract_candidates", wraps=assets.extract_candidates) as called:
            first = self.extract(asset, "first")
            second = self.extract(asset, "second")
        self.assertEqual(1, called.call_count)
        self.assertFalse(first["interval_cache"]["hit"])
        self.assertTrue(second["interval_cache"]["hit"])
        fresh = extract_candidates(asset["fasta"], asset["gff3"], "chr1:5-48", str(self.root / "fresh"),
                                   seqid_map_arg="chr1=1", full_hash=True, read_only_reference=True)
        for name in assets.INTERVAL_OUTPUTS.values():
            paths = [self.root / label / name for label in ("first", "second", "fresh")]
            self.assertEqual(paths[0].read_bytes(), paths[1].read_bytes())
            self.assertEqual(paths[1].read_bytes(), paths[2].read_bytes())
            self.assertNotEqual(paths[0].stat().st_ino, paths[1].stat().st_ino)
        self.assertEqual(fresh["gff_dialect_qc"], second["gff_dialect_qc"])
        with (self.root / "second/candidates.tsv").open() as stream:
            rows = {row["gene_id"]: row for row in csv.DictReader(stream, delimiter="\t")}
        self.assertEqual({"Plus", "Minus", "Dark", "Pseudo"}, set(rows))
        self.assertEqual("1", rows["Plus"]["start"])
        self.assertEqual("9", rows["Plus"]["end"])
        self.assertEqual("0", rows["Dark"]["protein_length"])
        self.assertEqual("pseudogene", rows["Pseudo"]["gene_biotype"])
        self.assertTrue(all(Path(value).parent == self.root / "second" for value in second["outputs"].values()))
        cache_dir = Path(second["interval_cache"]["manifest"]).parent
        self.assertFalse((cache_dir / "phenotype.json").exists())
        self.assertEqual(assets.INTERVAL_OUTPUTS, json.loads((cache_dir / assets.EXTRACT_MANIFEST).read_text())["outputs"])
        cached = (cache_dir / "candidate_proteins.fa").read_bytes()
        (self.root / "first/candidate_proteins.fa").write_text(">changed\nM\n", encoding="utf-8")
        self.assertEqual(cached, (cache_dir / "candidate_proteins.fa").read_bytes())

    def test_same_size_same_mtime_source_changes_never_hit(self):
        previous = self.prepare()
        for path, old, new in ((self.fasta, "ATGAAATAA", "ATGAAGTAA"), (self.gff, "ID=Dark", "ID=DArk")):
            with self.subTest(path=path.name):
                value = path.stat()
                text = path.read_text()
                path.write_text(text.replace(old, new, 1), encoding="utf-8")
                os.utime(path, ns=(value.st_atime_ns, value.st_mtime_ns))
                self.assertEqual(value.st_size, path.stat().st_size)
                current = self.prepare()
                self.assertFalse(current["cache_hit"])
                self.assertNotEqual(previous["manifest"], current["manifest"])
                previous = current

    def test_interval_mapping_and_implementation_changes_invalidate(self):
        asset = self.prepare()
        first = self.extract(asset, "first")
        changed_interval = assets.extract_reference_interval(asset, "chr1:6-47", str(self.root / "interval"), seqid_map_arg="chr1=1")
        changed_mapping = assets.extract_reference_interval(asset, "chr1:5-48", str(self.root / "mapping"),
                                                            seqid_map_arg="unused=1,chr1=1")
        implementation = assets._extraction_implementation()
        changed = json.loads(json.dumps(implementation))
        changed["code_sha256"]["core.py"] = "different-core"
        with patch.object(assets, "_extraction_implementation", return_value=changed):
            changed_code = self.extract(asset, "code")
        with patch.object(assets, "EXTRACTION_POLICY_VERSION", "new-extractor-contract"):
            changed_policy = self.extract(asset, "policy")
        changed = json.loads(json.dumps(implementation))
        changed["dependencies"]["pandas"] = "different-pandas"
        with patch.object(assets, "_extraction_implementation", return_value=changed):
            changed_dependency = self.extract(asset, "dependency")
        all_runs = (first, changed_interval, changed_mapping, changed_code, changed_policy, changed_dependency)
        self.assertEqual(len(all_runs), len({run["interval_cache"]["key"] for run in all_runs}))
        self.assertTrue(all(not run["interval_cache"]["hit"] for run in all_runs))

    def test_mapping_file_reuses_cache_but_raw_interval_format_is_preserved(self):
        asset = self.prepare()
        first = self.extract(asset, "first")
        mapping = self.root / "mapping.tsv"
        mapping.write_text("chr1\t1\n", encoding="utf-8")
        second = assets.extract_reference_interval(asset, "chr1:5-48", str(self.root / "second"), seqid_map_arg=str(mapping))
        self.assertTrue(second["interval_cache"]["hit"])
        self.assertEqual(first["interval_cache"]["key"], second["interval_cache"]["key"])
        raw = " chr1:005-048 "
        third = assets.extract_reference_interval(asset, raw, str(self.root / "raw"), seqid_map_arg=str(mapping))
        self.assertFalse(third["interval_cache"]["hit"])
        self.assertNotEqual(first["interval_cache"]["key"], third["interval_cache"]["key"])
        extract_candidates(asset["fasta"], asset["gff3"], raw, str(self.root / "raw-fresh"),
                           seqid_map_arg=str(mapping), full_hash=True, read_only_reference=True)
        for name in assets.INTERVAL_OUTPUTS.values():
            self.assertEqual((self.root / "raw" / name).read_bytes(), (self.root / "raw-fresh" / name).read_bytes())

    def test_managed_files_and_manifest_tampering_are_rejected(self):
        for name in ("reference.fa", "reference.fa.fai", "annotation.gff3", assets.ASSET_MANIFEST):
            with self.subTest(name=name):
                asset = self.prepare(self.root / ("cache-" + name))
                path = Path(asset["manifest"]).parent / name
                path.chmod(0o644)
                path.write_bytes(path.read_bytes() + b"changed")
                with self.assertRaises(ValueError):
                    self.prepare(self.root / ("cache-" + name))

    def test_incomplete_and_extra_reference_entries_fail_without_repair(self):
        for alteration in ("missing", "extra"):
            with self.subTest(alteration=alteration):
                cache = self.root / alteration
                asset = self.prepare(cache)
                entry = Path(asset["manifest"]).parent
                if alteration == "missing":
                    (entry / "reference.fa.fai").unlink()
                else:
                    (entry / "unexpected.txt").write_text("unexpected", encoding="utf-8")
                with patch.object(assets, "_copy", side_effect=AssertionError("must not repair a corrupted entry")):
                    with self.assertRaisesRegex(ValueError, "incomplete or unexpected"):
                        self.prepare(cache)

    def test_interval_corruption_missing_and_extra_files_fail_closed(self):
        asset = self.prepare()
        for alteration in ("changed", "missing", "extra"):
            with self.subTest(alteration=alteration):
                cache = self.root / ("interval-" + alteration)
                first = self.extract(asset, "first-" + alteration, interval_cache_dir_arg=str(cache))
                entry = Path(first["interval_cache"]["manifest"]).parent
                target = entry / "candidate_cds.fa"
                if alteration == "changed":
                    target.chmod(0o644)
                    target.write_bytes(target.read_bytes() + b"changed")
                elif alteration == "missing":
                    target.unlink()
                else:
                    (entry / "unexpected.txt").write_text("unexpected", encoding="utf-8")
                with patch.object(assets, "extract_candidates", side_effect=AssertionError("must not silently recompute corrupt cache")):
                    with self.assertRaises(ValueError):
                        self.extract(asset, "second-" + alteration, interval_cache_dir_arg=str(cache))

    def test_source_aliases_nested_symlinks_and_cache_overlap_are_rejected(self):
        original_bytes = self.fasta.read_bytes()
        alias = self.root / "alias.fa"
        alias.symlink_to(self.fasta)
        with self.assertRaisesRegex(ValueError, "symlink"):
            assets.prepare_reference_asset(str(alias), str(self.gff), str(self.cache), reference_version="v1")
        alias.unlink()
        alias.hardlink_to(self.fasta)
        with self.assertRaisesRegex(ValueError, "hardlinked"):
            self.prepare()
        alias.unlink()
        linked_directory = self.root / "linked-source"
        linked_directory.symlink_to(self.source, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, "symlink"):
            assets.prepare_reference_asset(str(linked_directory / self.fasta.name), str(self.gff), str(self.cache), reference_version="v1")
        with self.assertRaisesRegex(ValueError, "separate"):
            self.prepare(self.source)
        self.assertEqual(original_bytes, self.fasta.read_bytes())
        self.assertFalse(self.cache.exists())

    def test_output_aliases_and_mapping_symlink_are_rejected_before_writes(self):
        asset = self.prepare()
        for kind in ("symlink", "hardlink"):
            with self.subTest(kind=kind):
                out = self.root / kind
                out.mkdir()
                alias = out / "candidate_genomic.fa"
                if kind == "symlink":
                    alias.symlink_to(self.fasta)
                else:
                    alias.hardlink_to(self.fasta)
                with self.assertRaisesRegex(ValueError, "fresh output"):
                    self.extract(asset, kind)
                self.assertFalse((out / "candidates.tsv").exists())
                alias.unlink()
        mapping = self.root / "mapping.tsv"
        mapping.write_text("chr1\t1\n", encoding="utf-8")
        link = self.root / "mapping-link.tsv"
        link.symlink_to(mapping)
        with self.assertRaisesRegex(ValueError, "symlink"):
            assets.extract_reference_interval(asset, "chr1:5-48", str(self.root / "linked-map"), seqid_map_arg=str(link))
        with self.assertRaisesRegex(ValueError, "separate"):
            assets.extract_reference_interval(asset, "chr1:5-48", str(self.cache), seqid_map_arg="chr1=1")

    def test_failed_creation_never_publishes_or_reuses_partial_entry(self):
        with patch.object(assets, "_publish", side_effect=RuntimeError("simulated interruption")):
            with self.assertRaisesRegex(RuntimeError, "interruption"):
                self.prepare()
        self.assertEqual([], list((self.cache / "references").iterdir()))
        asset = self.prepare()
        with patch.object(assets, "_publish", side_effect=RuntimeError("simulated interruption")):
            with self.assertRaisesRegex(RuntimeError, "interruption"):
                self.extract(asset, "interrupted")
        self.assertEqual([], list((self.cache / "intervals").iterdir()))
        self.assertFalse(self.extract(asset, "retried")["interval_cache"]["hit"])

    def test_original_mutation_during_copy_never_publishes_an_asset(self):
        original_copy = assets._copy
        def mutate_after_copy(source, target):
            original_copy(source, target)
            if source == self.fasta:
                self.fasta.write_text(self.fasta.read_text().replace("ATGAAATAA", "ATGAAGTAA", 1), encoding="utf-8")
        with patch.object(assets, "_copy", side_effect=mutate_after_copy):
            with self.assertRaisesRegex(ValueError, "changed during use"):
                self.prepare()
        self.assertEqual([], list((self.cache / "references").iterdir()))

    def test_managed_reference_mutation_during_extract_never_publishes_cache(self):
        asset = self.prepare()
        original_extract = assets.extract_candidates
        def mutate_after_extract(*args, **kwargs):
            result = original_extract(*args, **kwargs)
            fasta = Path(asset["fasta"])
            fasta.chmod(0o644)
            fasta.write_text(fasta.read_text().replace("ATGAAATAA", "ATGAAGTAA", 1), encoding="utf-8")
            return result
        with patch.object(assets, "extract_candidates", side_effect=mutate_after_extract):
            with self.assertRaisesRegex(ValueError, "changed during use"):
                self.extract(asset, "concurrently-changed")
        self.assertEqual([], list((self.cache / "intervals").iterdir()))
        self.assertFalse((self.root / "concurrently-changed" / assets.EXTRACT_MANIFEST).exists())

    def test_simultaneous_creators_leave_one_complete_asset(self):
        entered, release = threading.Event(), threading.Event()
        original_copy = assets._copy
        def delayed_copy(source, target):
            if target.name == "reference.fa":
                entered.set()
                if not release.wait(5):
                    raise AssertionError("test creation barrier timed out")
            return original_copy(source, target)
        with ThreadPoolExecutor(max_workers=1) as executor, patch.object(assets, "_copy", side_effect=delayed_copy):
            future = executor.submit(self.prepare)
            try:
                self.assertTrue(entered.wait(5))
                with self.assertRaisesRegex(ValueError, "already in progress"):
                    self.prepare()
            finally:
                release.set()
            first = future.result(timeout=5)
        self.assertTrue(self.prepare()["cache_hit"])
        self.assertEqual([Path(first["manifest"]).parent], list((self.cache / "references").iterdir()))

    def test_pinned_asset_manifest_is_required_for_extraction(self):
        asset = self.prepare()
        altered = dict(asset, manifest_sha256="0" * 64)
        with self.assertRaisesRegex(ValueError, "identity/hash"):
            self.extract(altered, "bad-anchor")
        with self.assertRaisesRegex(ValueError, "paths do not match"):
            self.extract(dict(asset, fasta=str(self.fasta)), "bad-source")


class SharedProteinDatabaseTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.prefix = self.root / "proteins"
        for suffix in (".pin", ".phr", ".psq"):
            Path(str(self.prefix) + suffix).write_bytes(b"fake database bytes\n")

    def capture(self):
        return assets.capture_shared_protein_db(str(self.prefix), version="synthetic-db-v1")

    def test_direct_database_snapshot_is_json_serializable_and_read_only(self):
        files = list(self.root.iterdir())
        before = {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in files}
        snapshot = self.capture()
        restored = json.loads(json.dumps(snapshot))
        self.assertEqual(restored, assets.validate_shared_protein_db(restored))
        self.assertEqual(sorted(before), snapshot["file_roster"])
        self.assertTrue(all(value["strength"] == "strong" for value in snapshot["files"].values()))
        self.assertEqual(before, {p.name: (p.read_bytes(), p.stat().st_mtime_ns) for p in self.root.iterdir()})

    def test_multivolume_roster_includes_every_prefix_file(self):
        for path in list(self.root.iterdir()):
            path.unlink()
        for volume in ("00", "01"):
            for suffix in ("pin", "phr", "psq"):
                (self.root / f"proteins.{volume}.{suffix}").write_bytes(b"volume")
        (self.root / "proteins.pdb").write_bytes(b"v5 metadata")
        (self.root / "unrelated.psq").write_bytes(b"other DB")
        snapshot = self.capture()
        self.assertEqual(7, len(snapshot["file_roster"]))
        assets.validate_shared_protein_db(snapshot)

    def test_same_size_same_mtime_database_change_is_rejected(self):
        snapshot = self.capture()
        path = Path(str(self.prefix) + ".psq")
        before = path.stat()
        path.write_bytes(b"Fake database bytes\n")
        os.utime(path, ns=(before.st_atime_ns, before.st_mtime_ns))
        self.assertEqual(before.st_size, path.stat().st_size)
        with self.assertRaisesRegex(ValueError, "changed during use"):
            assets.validate_shared_protein_db(snapshot)

    def test_roster_addition_removal_and_metadata_changes_are_rejected(self):
        snapshot = self.capture()
        extra = Path(str(self.prefix) + ".pdb")
        extra.write_bytes(b"new file")
        with self.assertRaisesRegex(ValueError, "changed during use"):
            assets.validate_shared_protein_db(snapshot)
        extra.unlink()
        path = Path(str(self.prefix) + ".pin")
        value = path.stat()
        os.utime(path, ns=(value.st_atime_ns, value.st_mtime_ns + 1000))
        with self.assertRaisesRegex(ValueError, "changed during use"):
            assets.validate_shared_protein_db(snapshot)
        path.unlink()
        with self.assertRaisesRegex(ValueError, "complete"):
            assets.validate_shared_protein_db(snapshot)

    def test_alias_database_and_file_aliases_are_rejected(self):
        alias = Path(str(self.prefix) + ".pal")
        alias.write_text("DBLIST another_database\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "aliases"):
            self.capture()
        alias.unlink()
        target = Path(str(self.prefix) + ".psq")
        target.unlink()
        target.symlink_to(Path(str(self.prefix) + ".pin"))
        with self.assertRaisesRegex(ValueError, "symlink"):
            self.capture()


if __name__ == "__main__":
    unittest.main()
