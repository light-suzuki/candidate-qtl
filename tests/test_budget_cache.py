from __future__ import annotations

import json
import tempfile
import threading
import time
import unittest
from pathlib import Path
from unittest import mock

from candidate_qtl import budget_cache as budget_cache_mod
from candidate_qtl.budget_cache import CACHE_SCHEMA_VERSION, ContentCache
from candidate_qtl.core import sha256


class ContentCacheTests(unittest.TestCase):
    def _inputs(self) -> tuple[dict[str, str], str, str, str, dict[str, object]]:
        return {"fa": "a" * 64}, "diamond", "local", "1.2.3", {"evalue": 1e-5}

    def test_put_get_round_trip(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cache_dir = Path(tmp) / "cache"
            payload = Path(tmp) / "payload.bin"
            payload.write_bytes(b"result payload bytes")
            cache = ContentCache(cache_dir)
            input_hashes, tool, provider, provider_version, params = self._inputs()
            key = cache.put(input_hashes, tool, provider, provider_version, params, payload)
            self.assertEqual(64, len(key))
            hit, result_path = cache.get(input_hashes, tool, provider, provider_version, params)
            self.assertTrue(hit)
            self.assertIsInstance(result_path, Path)
            self.assertTrue(result_path.is_file())
            self.assertEqual(b"result payload bytes", result_path.read_bytes())
            stats = cache.stats()
            self.assertEqual(1, stats["hit_count"])
            self.assertEqual(0, stats["miss_count"])
            self.assertEqual(0, stats["corrupt_reject_count"])
            self.assertEqual(0, stats["stale_reject_count"])
            self.assertEqual("enabled", stats["cache_status"])

    def test_miss_on_absent_entry(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cache = ContentCache(Path(tmp) / "cache")
            input_hashes, tool, provider, provider_version, params = self._inputs()
            hit, result_path = cache.get(input_hashes, tool, provider, provider_version, params)
            self.assertFalse(hit)
            self.assertIsNone(result_path)
            self.assertEqual(1, cache.stats()["miss_count"])
            self.assertEqual(0, cache.stats()["hit_count"])

    def test_corrupt_payload_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cache_dir = Path(tmp) / "cache"
            payload = Path(tmp) / "payload.bin"
            payload.write_bytes(b"original payload")
            cache = ContentCache(cache_dir)
            input_hashes, tool, provider, provider_version, params = self._inputs()
            key = cache.put(input_hashes, tool, provider, provider_version, params, payload)
            payload_path = cache_dir / f"{key}.payload"
            raw = bytearray(payload_path.read_bytes())
            raw[0] ^= 0xFF
            payload_path.write_bytes(bytes(raw))
            hit, result_path = cache.get(input_hashes, tool, provider, provider_version, params)
            self.assertFalse(hit)
            self.assertIsNone(result_path)
            stats = cache.stats()
            self.assertEqual(1, stats["corrupt_reject_count"])
            self.assertEqual(1, stats["miss_count"])
            self.assertEqual(0, stats["hit_count"])

    def test_stale_provider_version_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cache_dir = Path(tmp) / "cache"
            payload = Path(tmp) / "payload.bin"
            payload.write_bytes(b"payload")
            cache = ContentCache(cache_dir)
            input_hashes, tool, provider, provider_version, params = self._inputs()
            key = cache.put(input_hashes, tool, provider, provider_version, params, payload)
            meta_path = cache_dir / f"{key}.json"
            metadata = json.loads(meta_path.read_text(encoding="utf-8"))
            metadata["provider_version"] = "9.9.9"
            meta_path.write_text(json.dumps(metadata, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            hit, result_path = cache.get(input_hashes, tool, provider, provider_version, params)
            self.assertFalse(hit)
            self.assertIsNone(result_path)
            stats = cache.stats()
            self.assertEqual(1, stats["stale_reject_count"])
            self.assertEqual(0, stats["corrupt_reject_count"])
            self.assertEqual(0, stats["miss_count"])

    def test_schema_version_mismatch_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cache_dir = Path(tmp) / "cache"
            payload = Path(tmp) / "payload.bin"
            payload.write_bytes(b"payload")
            cache = ContentCache(cache_dir)
            input_hashes, tool, provider, provider_version, params = self._inputs()
            key = cache.put(input_hashes, tool, provider, provider_version, params, payload)
            meta_path = cache_dir / f"{key}.json"
            metadata = json.loads(meta_path.read_text(encoding="utf-8"))
            metadata["schema_version"] = "0.9"
            meta_path.write_text(json.dumps(metadata, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            self.assertEqual(CACHE_SCHEMA_VERSION, "1.0")
            hit, result_path = cache.get(input_hashes, tool, provider, provider_version, params)
            self.assertFalse(hit)
            self.assertIsNone(result_path)
            stats = cache.stats()
            self.assertEqual(1, stats["stale_reject_count"])
            self.assertEqual(1, stats["miss_count"])

    def test_input_hashes_change_is_miss_not_stale(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cache_dir = Path(tmp) / "cache"
            payload = Path(tmp) / "payload.bin"
            payload.write_bytes(b"payload")
            cache = ContentCache(cache_dir)
            input_hashes, tool, provider, provider_version, params = self._inputs()
            cache.put(input_hashes, tool, provider, provider_version, params, payload)
            changed = dict(input_hashes)
            changed["fa"] = "b" * 64
            hit, result_path = cache.get(changed, tool, provider, provider_version, params)
            self.assertFalse(hit)
            self.assertIsNone(result_path)
            stats = cache.stats()
            self.assertEqual(1, stats["miss_count"])
            self.assertEqual(0, stats["stale_reject_count"])
            self.assertEqual(0, stats["hit_count"])

    def test_tool_provider_mismatch_rejected_as_stale(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cache_dir = Path(tmp) / "cache"
            payload = Path(tmp) / "payload.bin"
            payload.write_bytes(b"payload")
            cache = ContentCache(cache_dir)
            input_hashes, tool, provider, provider_version, params = self._inputs()
            key = cache.put(input_hashes, tool, provider, provider_version, params, payload)
            for field, wrong in (("tool", "other-tool"), ("provider", "other-provider")):
                meta_path = cache_dir / f"{key}.json"
                metadata = json.loads(meta_path.read_text(encoding="utf-8"))
                metadata[field] = wrong
                meta_path.write_text(json.dumps(metadata, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
                hit, result_path = cache.get(input_hashes, tool, provider, provider_version, params)
                self.assertFalse(hit)
                self.assertIsNone(result_path)
                self.assertEqual(0, cache.stats()["hit_count"])
            stats = cache.stats()
            self.assertEqual(2, stats["stale_reject_count"])
            self.assertEqual(0, stats["miss_count"])
            self.assertEqual(0, stats["corrupt_reject_count"])

    def test_key_determinism(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cache = ContentCache(Path(tmp) / "cache")
            input_hashes, tool, provider, provider_version, params = self._inputs()
            first = cache.key(input_hashes, tool, provider, provider_version, params)
            second = cache.key(input_hashes, tool, provider, provider_version, params)
            self.assertEqual(first, second)
            different = cache.key(input_hashes, tool, provider, provider_version, {"evalue": 1e-4})
            self.assertNotEqual(first, different)
            self.assertEqual(64, len(first))
            self.assertTrue(set(first) <= set("0123456789abcdef"))

    def test_no_leftover_tmp_files_after_put(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cache_dir = Path(tmp) / "cache"
            payload = Path(tmp) / "payload.bin"
            payload.write_bytes(b"payload")
            cache = ContentCache(cache_dir)
            input_hashes, tool, provider, provider_version, params = self._inputs()
            key = cache.put(input_hashes, tool, provider, provider_version, params, payload)
            self.assertEqual([], [p.name for p in cache_dir.glob("*.tmp")])
            self.assertTrue((cache_dir / f"{key}.json").is_file())
            self.assertTrue((cache_dir / f"{key}.payload").is_file())

    def test_cache_dir_created_lazily_on_first_put(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cache_dir = Path(tmp) / "deep" / "nested" / "cache"
            cache = ContentCache(cache_dir)
            self.assertFalse(cache_dir.exists())
            payload = Path(tmp) / "payload.bin"
            payload.write_bytes(b"x")
            input_hashes, tool, provider, provider_version, params = self._inputs()
            cache.put(input_hashes, tool, provider, provider_version, params, payload)
            self.assertTrue(cache_dir.is_dir())

    def test_disabled_cache_is_noop(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cache_dir = Path(tmp) / "cache"
            payload = Path(tmp) / "payload.bin"
            payload.write_bytes(b"payload")
            cache = ContentCache(cache_dir, enabled=False)
            input_hashes, tool, provider, provider_version, params = self._inputs()
            key = cache.put(input_hashes, tool, provider, provider_version, params, payload)
            self.assertEqual(64, len(key))
            self.assertFalse(cache_dir.exists())
            hit, result_path = cache.get(input_hashes, tool, provider, provider_version, params)
            self.assertFalse(hit)
            self.assertIsNone(result_path)
            stats = cache.stats()
            self.assertEqual("disabled", stats["cache_status"])
            self.assertEqual(0, stats["hit_count"])
            self.assertEqual(0, stats["miss_count"])
            self.assertEqual(0, stats["corrupt_reject_count"])
            self.assertEqual(0, stats["stale_reject_count"])


class ContentCacheConcurrencyTests(unittest.TestCase):
    def _inputs(self) -> tuple[dict[str, str], str, str, str, dict[str, object]]:
        return {"fa": "a" * 64}, "diamond", "local", "1.2.3", {"evalue": 1e-5}

    def test_concurrent_put_same_key_16_threads(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cache_dir = Path(tmp) / "cache"
            payload = Path(tmp) / "payload.bin"
            expected = bytes(range(256)) * 64
            payload.write_bytes(expected)
            cache = ContentCache(cache_dir)
            input_hashes, tool, provider, provider_version, params = self._inputs()
            barrier = threading.Barrier(16)
            errors: list[BaseException] = []

            def worker() -> None:
                try:
                    barrier.wait(timeout=30)
                    cache.put(input_hashes, tool, provider, provider_version, params, payload)
                except BaseException as exc:
                    errors.append(exc)

            threads = [threading.Thread(target=worker) for _ in range(16)]
            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join(timeout=60)
            self.assertFalse(any(thread.is_alive() for thread in threads))
            self.assertEqual([], errors)
            stats = cache.stats()
            self.assertEqual(0, stats["hit_count"])
            self.assertEqual(0, stats["miss_count"])
            hit, result_path = cache.get(input_hashes, tool, provider, provider_version, params)
            self.assertTrue(hit)
            self.assertIsInstance(result_path, Path)
            self.assertEqual(expected, result_path.read_bytes())
            self.assertEqual([], [p.name for p in cache_dir.glob("*.tmp")])
            self.assertTrue((cache_dir / f"{cache.key(input_hashes, tool, provider, provider_version, params)}.json").is_file())

    def test_reader_never_sees_partial_payload_as_hit(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cache_dir = Path(tmp) / "cache"
            payload = Path(tmp) / "payload.bin"
            expected = b"complete payload bytes"
            payload.write_bytes(expected)
            cache = ContentCache(cache_dir)
            input_hashes, tool, provider, provider_version, params = self._inputs()
            real_replace = budget_cache_mod.os.replace

            def slow_replace(src, dst):
                if Path(dst).name.endswith(".payload"):
                    time.sleep(0.2)
                return real_replace(src, dst)

            with mock.patch("candidate_qtl.budget_cache.os.replace", side_effect=slow_replace):
                put_thread = threading.Thread(target=cache.put, args=(input_hashes, tool, provider, provider_version, params, payload))
                put_thread.start()
                deadline = time.monotonic() + 10
                while time.monotonic() < deadline:
                    hit, result_path = cache.get(input_hashes, tool, provider, provider_version, params)
                    if hit:
                        self.assertEqual(expected, Path(result_path).read_bytes())
                    else:
                        self.assertIsNone(result_path)
                    if not put_thread.is_alive():
                        break
                    time.sleep(0.005)
                put_thread.join(timeout=30)
                self.assertFalse(put_thread.is_alive())
            hit, result_path = cache.get(input_hashes, tool, provider, provider_version, params)
            self.assertTrue(hit)
            self.assertIsInstance(result_path, Path)
            self.assertEqual(expected, result_path.read_bytes())

    def test_put_writer_unique_temp_no_fixed_names(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cache_dir = Path(tmp) / "cache"
            payload = Path(tmp) / "payload.bin"
            payload.write_bytes(b"payload")
            cache = ContentCache(cache_dir)
            input_hashes, tool, provider, provider_version, params = self._inputs()
            key = cache.key(input_hashes, tool, provider, provider_version, params)
            real_mkstemp = tempfile.mkstemp
            created: list[str] = []

            def record_mkstemp(*args, **kwargs):
                fd, name = real_mkstemp(*args, **kwargs)
                created.append(name)
                return fd, name

            with mock.patch("candidate_qtl.budget_cache.tempfile.mkstemp", side_effect=record_mkstemp):
                cache.put(input_hashes, tool, provider, provider_version, params, payload)
            self.assertGreaterEqual(2, len(created))
            self.assertEqual(len(created), len(set(created)))
            for name in created:
                self.assertTrue(name.startswith(str(cache_dir)), name)
                self.assertNotEqual(f"{key}.json.tmp", Path(name).name)
                self.assertTrue(Path(name).name.endswith(".tmp"), name)
            self.assertEqual([], [p.name for p in cache_dir.glob("*.tmp")])
            self.assertTrue((cache_dir / f"{key}.json").is_file())
            self.assertTrue((cache_dir / f"{key}.payload").is_file())

    def test_put_cleans_temps_on_forced_failure(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cache_dir = Path(tmp) / "cache"
            payload = Path(tmp) / "payload.bin"
            payload.write_bytes(b"payload")
            cache = ContentCache(cache_dir)
            input_hashes, tool, provider, provider_version, params = self._inputs()
            key = cache.key(input_hashes, tool, provider, provider_version, params)
            with mock.patch.object(Path, "read_bytes", side_effect=OSError("boom")):
                with self.assertRaises(OSError):
                    cache.put(input_hashes, tool, provider, provider_version, params, payload)
            self.assertEqual([], [p.name for p in cache_dir.glob("*.tmp")])
            self.assertFalse((cache_dir / f"{key}.payload").exists())
            self.assertFalse((cache_dir / f"{key}.json").exists())

    def test_concurrent_put_then_digest_verified_hit(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cache_dir = Path(tmp) / "cache"
            payload = Path(tmp) / "payload.bin"
            expected = b"digest check payload"
            payload.write_bytes(expected)
            cache = ContentCache(cache_dir)
            input_hashes, tool, provider, provider_version, params = self._inputs()
            barrier = threading.Barrier(16)
            errors: list[BaseException] = []

            def worker() -> None:
                try:
                    barrier.wait(timeout=30)
                    cache.put(input_hashes, tool, provider, provider_version, params, payload)
                except BaseException as exc:
                    errors.append(exc)

            threads = [threading.Thread(target=worker) for _ in range(16)]
            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join(timeout=60)
            self.assertEqual([], errors)
            self.assertEqual(0, cache.stats()["hit_count"])
            self.assertEqual(0, cache.stats()["miss_count"])
            key = cache.key(input_hashes, tool, provider, provider_version, params)
            hit, result_path = cache.get(input_hashes, tool, provider, provider_version, params)
            self.assertTrue(hit)
            self.assertEqual(expected, result_path.read_bytes())
            metadata = json.loads((cache_dir / f"{key}.json").read_text(encoding="utf-8"))
            self.assertEqual(sha256(result_path), metadata["payload_sha256"])
            stats = cache.stats()
            self.assertEqual(1, stats["hit_count"])
            self.assertEqual(0, stats["miss_count"])
            self.assertEqual(0, stats["corrupt_reject_count"])
            self.assertEqual(0, stats["stale_reject_count"])


if __name__ == "__main__":
    unittest.main()
