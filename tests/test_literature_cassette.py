from __future__ import annotations

import io
import json
import tempfile
import unittest
import urllib.error
from pathlib import Path

import pandas as pd

from candidate_qtl.literature_cassette import (
    CASSETTE_CONTRACT,
    request_key_sha256,
    retrieve_europe_pmc,
    write_query_provenance,
)


class _Response:
    def __init__(
        self, payload: object, status: int = 200, raw: bytes | None = None
    ) -> None:
        self.status = status
        self.raw = raw if raw is not None else json.dumps(payload).encode("utf-8")

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return False

    def read(self) -> bytes:
        return self.raw


class LiteratureCassetteTests(unittest.TestCase):
    @staticmethod
    def _payload(*articles: dict[str, object]) -> dict[str, object]:
        return {"resultList": {"result": list(articles)}}

    def test_live_then_offline_replay_preserves_provider_timestamp_and_bytes(
        self,
    ) -> None:
        calls: list[tuple[str, int]] = []

        def opener(request, timeout=0):
            calls.append((request.full_url, timeout))
            return _Response(self._payload({"pmid": "1", "title": "Paper"}))

        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            first = retrieve_europe_pmc(
                '"CLV1"',
                5,
                write_dir_arg=root,
                allow_external_api=True,
                opener=opener,
                now=lambda: "2026-09-03T00:00:00+00:00",
            )
            replay = retrieve_europe_pmc(
                '"CLV1"',
                5,
                cache_arg=str(root),
                write_dir_arg=root / "unused",
                allow_external_api=False,
                opener=lambda *_args, **_kwargs: (_ for _ in ()).throw(
                    AssertionError("network used")
                ),
                now=lambda: "2099-01-01T00:00:00+00:00",
            )
            self.assertEqual(1, len(calls))
            self.assertEqual("retrieved", first.retrieval_status)
            self.assertEqual("live_europe_pmc", first.retrieval_source)
            self.assertEqual("live_europe_pmc", replay.retrieval_source)
            self.assertEqual("live_then_persisted", first.execution_source)
            self.assertEqual("immutable_cassette", replay.execution_source)
            self.assertEqual(
                first.provider_retrieved_at, replay.provider_retrieved_at
            )
            self.assertEqual(
                first.raw_response_sha256, replay.raw_response_sha256
            )
            self.assertEqual(first.articles, replay.articles)
            cassette = json.loads(
                Path(first.cassette_path).read_text(encoding="utf-8")
            )
            self.assertEqual(CASSETTE_CONTRACT, cassette["cassette_contract"])

    def test_no_hit_is_persisted_and_replayed(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            first = retrieve_europe_pmc(
                "nothing",
                3,
                write_dir_arg=root,
                allow_external_api=True,
                opener=lambda *_args, **_kwargs: _Response(self._payload()),
                now=lambda: "2026-09-03T00:00:01+00:00",
            )
            replay = retrieve_europe_pmc(
                "nothing",
                3,
                cache_arg=str(root),
                write_dir_arg=root,
                allow_external_api=False,
            )
            self.assertEqual("no_hit", first.retrieval_status)
            self.assertEqual("no_hit", replay.retrieval_status)
            self.assertEqual((), replay.articles)

    def test_request_identity_includes_page_size(self) -> None:
        self.assertNotEqual(
            request_key_sha256("query", 5),
            request_key_sha256("query", 10),
        )

    def test_rate_limit_provider_failure_and_parse_failure_remain_distinct(
        self,
    ) -> None:
        def http_error(code: int):
            def opener(request, timeout=0):
                raise urllib.error.HTTPError(
                    request.full_url,
                    code,
                    "error",
                    hdrs=None,
                    fp=io.BytesIO(b"provider body"),
                )

            return opener

        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            rate = retrieve_europe_pmc(
                "rate",
                1,
                write_dir_arg=root / "rate",
                allow_external_api=True,
                opener=http_error(429),
                max_attempts=2,
                sleep=lambda _seconds: None,
            )
            server = retrieve_europe_pmc(
                "server",
                1,
                write_dir_arg=root / "server",
                allow_external_api=True,
                opener=http_error(503),
                max_attempts=1,
            )
            parse = retrieve_europe_pmc(
                "parse",
                1,
                write_dir_arg=root / "parse",
                allow_external_api=True,
                opener=lambda *_args, **_kwargs: _Response(
                    {}, raw=b"not-json"
                ),
            )
            self.assertEqual("rate_limited", rate.retrieval_status)
            self.assertEqual(2, rate.attempt_count)
            self.assertEqual("provider_unavailable", server.retrieval_status)
            self.assertEqual("response_parse_error", parse.retrieval_status)
            parse_payload = json.loads(
                Path(parse.cassette_path).read_text(encoding="utf-8")
            )
            self.assertEqual(
                "not-json",
                parse_payload["raw_response"]["unparsed_response_text"],
            )

    def test_legacy_cache_is_read_without_rewriting_current_time(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            legacy = root / "legacy.json"
            legacy.write_text(
                json.dumps(
                    {
                        "query": "legacy query",
                        "retrieved_at": "2020-01-01T00:00:00+00:00",
                        "results": [{"pmid": "2"}],
                    }
                ),
                encoding="utf-8",
            )
            result = retrieve_europe_pmc(
                "legacy query",
                5,
                cache_arg=str(legacy),
                write_dir_arg=root / "new",
                allow_external_api=False,
                now=lambda: "2099-01-01T00:00:00+00:00",
            )
            self.assertTrue(result.legacy_cache)
            self.assertEqual("legacy_query_cache", result.retrieval_source)
            self.assertEqual(
                "2020-01-01T00:00:00+00:00",
                result.provider_retrieved_at,
            )

    def test_query_provenance_is_canonically_sorted(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            first_retrieval = retrieve_europe_pmc(
                "b",
                1,
                write_dir_arg=root / "b",
                allow_external_api=False,
            )
            second_retrieval = retrieve_europe_pmc(
                "a",
                1,
                write_dir_arg=root / "a",
                allow_external_api=False,
            )
            first = Path(
                write_query_provenance(
                    [first_retrieval, second_retrieval], root / "first.tsv"
                )
            )
            second = Path(
                write_query_provenance(
                    [second_retrieval, first_retrieval], root / "second.tsv"
                )
            )
            self.assertEqual(first.read_bytes(), second.read_bytes())
            frame = pd.read_csv(first, sep="\t", dtype=str).fillna("")
            self.assertEqual(2, len(frame))


if __name__ == "__main__":
    unittest.main()
