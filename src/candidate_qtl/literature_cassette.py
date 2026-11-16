"""Immutable request/response cassettes for Europe PMC literature retrieval.

The cassette is the shared reproducibility boundary for Route B literature
lookup and Direct Route C search (issue #662). It separates provider retrieval
time from current execution time, preserves raw response content and errors,
and makes per-query cache/live provenance explicit.
"""

from __future__ import annotations

import hashlib
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

from .core import runtime_path, safe_id

CASSETTE_CONTRACT = "literature-cassette-1.0.0"
REQUEST_SCHEMA_VERSION = "europe-pmc-search-request-1.0.0"
PARSER_SCHEMA_VERSION = "europe-pmc-search-parser-1.0.0"
EUROPE_PMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
EUROPE_PMC_RESULT_TYPE = "core"
DEFAULT_USER_AGENT = "candidate-qtl literature-cassette (read-only literature lookup)"


class ResponseParseError(ValueError):
    """Provider returned bytes that were not valid JSON."""

    def __init__(self, message: str, raw_text: str) -> None:
        super().__init__(message)
        self.raw_text = raw_text


RETRIEVAL_STATES = frozenset(
    {
        "retrieved",
        "no_hit",
        "not_evaluated_offline",
        "invalid_request",
        "rate_limited",
        "provider_unavailable",
        "response_parse_error",
    }
)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _canonical_json(value: object) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )


def _sha256_json(value: object) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def normalise_request(query: str, limit: int) -> dict[str, object]:
    return {
        "request_schema_version": REQUEST_SCHEMA_VERSION,
        "provider": "europe_pmc_rest",
        "endpoint": EUROPE_PMC_SEARCH,
        "request": {
            "query": str(query),
            "format": "json",
            "resultType": EUROPE_PMC_RESULT_TYPE,
            "pageSize": int(limit),
            "sort": "",
        },
        "parser_schema_version": PARSER_SCHEMA_VERSION,
    }


def request_key_sha256(query: str, limit: int) -> str:
    return _sha256_json(normalise_request(query, limit))


def _articles_from_raw(raw_response: object) -> list[dict[str, object]]:
    if not isinstance(raw_response, dict):
        raise ValueError("Europe PMC raw response must be a JSON object")
    result_list = raw_response.get("resultList")
    if not isinstance(result_list, dict):
        raise ValueError("Europe PMC response lacks resultList")
    results = result_list.get("result")
    if not isinstance(results, list):
        raise ValueError("Europe PMC response resultList.result must be a list")
    if any(not isinstance(record, dict) for record in results):
        raise ValueError("Europe PMC response contains a non-object result")
    return [record for record in results if isinstance(record, dict)]


@dataclass(frozen=True)
class LiteratureRetrieval:
    query: str
    request_key_sha256: str
    retrieval_status: str
    articles: tuple[dict[str, object], ...]
    provider: str
    retrieval_source: str
    provider_retrieved_at: str
    raw_response_sha256: str
    cassette_path: str
    execution_source: str
    http_status: int | None
    attempt_count: int
    error: str
    legacy_cache: bool = False

    def provenance_row(self) -> dict[str, object]:
        return {
            "query": self.query,
            "request_key_sha256": self.request_key_sha256,
            "retrieval_status": self.retrieval_status,
            "article_count": len(self.articles),
            "provider": self.provider,
            "retrieval_source": self.retrieval_source,
            "provider_retrieved_at": self.provider_retrieved_at,
            "raw_response_sha256": self.raw_response_sha256,
            "cassette_path": Path(self.cassette_path).name if self.cassette_path else "",
            "http_status": "" if self.http_status is None else self.http_status,
            "attempt_count": self.attempt_count,
            "error": self.error,
            "legacy_cache": self.legacy_cache,
        }


def _cassette_payload(
    request: dict[str, object],
    *,
    status: str,
    retrieved_at: str,
    http_status: int | None,
    raw_response: object,
    attempt_count: int,
    error: str,
) -> dict[str, object]:
    key = _sha256_json(request)
    return {
        "cassette_contract": CASSETTE_CONTRACT,
        **request,
        "request_key_sha256": key,
        "retrieval_status": status,
        "provider_retrieved_at": retrieved_at,
        "http_status": http_status,
        "attempt_count": int(attempt_count),
        "raw_response": raw_response,
        "raw_response_sha256": _sha256_json(raw_response),
        "error": str(error),
        "provider_origin": "live_europe_pmc",
    }


def _atomic_write(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rendered = json.dumps(
        payload, ensure_ascii=False, indent=2, sort_keys=True
    ) + "\n"
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(rendered, encoding="utf-8")
    if path.exists():
        existing = path.read_text(encoding="utf-8")
        if existing != rendered:
            temporary.unlink(missing_ok=True)
            raise ValueError(
                f"immutable literature cassette already exists with different bytes: {path}"
            )
        temporary.unlink(missing_ok=True)
        return
    temporary.replace(path)


def _load_cassette(path: Path, expected_request: dict[str, object]) -> LiteratureRetrieval:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"literature cassette is not valid JSON: {path}") from error
    if not isinstance(payload, dict):
        raise ValueError("literature cassette must contain a JSON object")
    if payload.get("cassette_contract") != CASSETTE_CONTRACT:
        raise ValueError("unsupported literature cassette contract")
    actual_request = {
        "request_schema_version": payload.get("request_schema_version"),
        "provider": payload.get("provider"),
        "endpoint": payload.get("endpoint"),
        "request": payload.get("request"),
        "parser_schema_version": payload.get("parser_schema_version"),
    }
    if actual_request != expected_request:
        raise ValueError("literature cassette request does not match current request")
    key = _sha256_json(actual_request)
    if str(payload.get("request_key_sha256", "")) != key:
        raise ValueError("literature cassette request-key hash mismatch")
    raw = payload.get("raw_response")
    if str(payload.get("raw_response_sha256", "")) != _sha256_json(raw):
        raise ValueError("literature cassette raw-response hash mismatch")
    status = str(payload.get("retrieval_status", ""))
    if status not in RETRIEVAL_STATES:
        raise ValueError(f"invalid literature cassette retrieval_status: {status}")
    articles: list[dict[str, object]] = []
    if status in {"retrieved", "no_hit"}:
        articles = _articles_from_raw(raw)
        if (status == "retrieved") != bool(articles):
            raise ValueError("literature cassette status/article-count mismatch")
    http_status_value = payload.get("http_status")
    try:
        http_status = None if http_status_value in (None, "") else int(http_status_value)
    except (TypeError, ValueError) as error:
        raise ValueError("literature cassette http_status must be an integer or null") from error
    return LiteratureRetrieval(
        query=str(expected_request["request"]["query"]),  # type: ignore[index]
        request_key_sha256=key,
        retrieval_status=status,
        articles=tuple(articles),
        provider="europe_pmc_rest",
        retrieval_source=str(payload.get("provider_origin", "live_europe_pmc")),
        provider_retrieved_at=str(payload.get("provider_retrieved_at", "")),
        raw_response_sha256=str(payload.get("raw_response_sha256", "")),
        cassette_path=str(path),
        execution_source="immutable_cassette",
        http_status=http_status,
        attempt_count=int(payload.get("attempt_count", 0) or 0),
        error=str(payload.get("error", "")),
    )


def _new_cassette_candidates(
    cache_arg: str | None,
    write_dir: Path,
    key: str,
) -> list[Path]:
    candidates: list[Path] = []
    if cache_arg:
        cache = runtime_path(cache_arg)
        if cache.is_dir():
            candidates.append(cache / f"{key}.json")
        elif cache.is_file():
            candidates.append(cache)
    candidates.append(write_dir / f"{key}.json")
    unique: list[Path] = []
    seen: set[str] = set()
    for path in candidates:
        token = str(path.resolve()) if path.exists() else str(path.absolute())
        if token not in seen:
            seen.add(token)
            unique.append(path)
    return unique


def _legacy_query_cache(
    cache_arg: str | None,
    query: str,
    key: str,
) -> LiteratureRetrieval | None:
    if not cache_arg:
        return None
    path = runtime_path(cache_arg)
    files = sorted(path.glob("*.json")) if path.is_dir() else [path]
    for item in files:
        try:
            payload = json.loads(item.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if isinstance(payload, dict) and payload.get("cassette_contract"):
            continue
        articles: list[dict[str, object]] | None = None
        retrieved_at = ""
        if isinstance(payload, dict) and str(payload.get("query", "")) == query:
            values = payload.get("results", [])
            if isinstance(values, list):
                articles = [value for value in values if isinstance(value, dict)]
                retrieved_at = str(
                    payload.get("provider_retrieved_at", payload.get("retrieved_at", ""))
                )
        elif isinstance(payload, dict):
            for legacy_query, values in payload.items():
                if (
                    str(legacy_query) == query
                    or safe_id(str(legacy_query)) == safe_id(query)
                ) and isinstance(values, list):
                    articles = [value for value in values if isinstance(value, dict)]
                    break
        if articles is None:
            continue
        return LiteratureRetrieval(
            query=query,
            request_key_sha256=key,
            retrieval_status="retrieved" if articles else "no_hit",
            articles=tuple(articles),
            provider="europe_pmc_rest",
            retrieval_source="legacy_query_cache",
            provider_retrieved_at=retrieved_at,
            raw_response_sha256=_sha256_json(
                {"query": query, "results": articles}
            ),
            cassette_path=str(item),
            execution_source="legacy_query_cache",
            http_status=None,
            attempt_count=0,
            error="",
            legacy_cache=True,
        )
    return None


def _http_fetch(
    request_payload: dict[str, object],
    *,
    timeout: int,
    user_agent: str,
    opener: Callable[..., object],
) -> tuple[int, object]:
    params = request_payload["request"]
    if not isinstance(params, dict):
        raise ValueError("normalised literature request lacks request parameters")
    query_params = {key: value for key, value in params.items() if value != ""}
    url = f"{EUROPE_PMC_SEARCH}?" + urllib.parse.urlencode(query_params)
    request = urllib.request.Request(url, headers={"User-Agent": user_agent})
    with opener(request, timeout=timeout) as response:  # type: ignore[attr-defined]
        status = int(getattr(response, "status", 200) or 200)
        raw_bytes = response.read()
    try:
        raw_text = raw_bytes.decode("utf-8")
    except UnicodeDecodeError:
        raw_text = raw_bytes.decode("utf-8", errors="replace")
        raise ResponseParseError("Europe PMC returned non-UTF-8 JSON", raw_text)
    try:
        return status, json.loads(raw_text)
    except json.JSONDecodeError as error:
        raise ResponseParseError("Europe PMC returned invalid JSON", raw_text) from error


def retrieve_europe_pmc(
    query: str,
    limit: int,
    *,
    cache_arg: str | None = None,
    write_dir_arg: str | Path,
    allow_external_api: bool = False,
    timeout: int = 30,
    max_attempts: int = 3,
    user_agent: str = DEFAULT_USER_AGENT,
    opener: Callable[..., object] | None = None,
    now: Callable[[], str] = _now,
    sleep: Callable[[float], None] = time.sleep,
) -> LiteratureRetrieval:
    """Load an exact request cassette or retrieve and atomically persist one."""
    query_text = str(query).strip()
    opener = opener or urllib.request.urlopen
    write_dir = runtime_path(str(write_dir_arg))
    request_payload = normalise_request(query_text, limit)
    key = _sha256_json(request_payload)
    if not query_text or limit < 1:
        return LiteratureRetrieval(
            query=query_text,
            request_key_sha256=key,
            retrieval_status="invalid_request",
            articles=(),
            provider="europe_pmc_rest",
            retrieval_source="validation",
            provider_retrieved_at="",
            raw_response_sha256="",
            cassette_path="",
            execution_source="validation",
            http_status=None,
            attempt_count=0,
            error="query must be non-empty and limit must be at least 1",
        )

    for candidate in _new_cassette_candidates(cache_arg, write_dir, key):
        if not candidate.is_file():
            continue
        try:
            payload = json.loads(candidate.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if (
            isinstance(payload, dict)
            and payload.get("cassette_contract") == CASSETTE_CONTRACT
        ):
            return _load_cassette(candidate, request_payload)

    legacy = _legacy_query_cache(cache_arg, query_text, key)
    if legacy is not None:
        return legacy
    if not allow_external_api:
        return LiteratureRetrieval(
            query=query_text,
            request_key_sha256=key,
            retrieval_status="not_evaluated_offline",
            articles=(),
            provider="europe_pmc_rest",
            retrieval_source="offline",
            provider_retrieved_at="",
            raw_response_sha256="",
            cassette_path="",
            execution_source="offline",
            http_status=None,
            attempt_count=0,
            error="",
        )

    status = "provider_unavailable"
    http_status: int | None = None
    raw_response: object = {}
    error_text = ""
    attempts = 0
    retrieved_at = now()
    for attempt in range(1, max(1, max_attempts) + 1):
        attempts = attempt
        try:
            http_status, raw_response = _http_fetch(
                request_payload,
                timeout=timeout,
                user_agent=user_agent,
                opener=opener,
            )
            articles = _articles_from_raw(raw_response)
            status = "retrieved" if articles else "no_hit"
            error_text = ""
            break
        except urllib.error.HTTPError as error:
            http_status = int(error.code)
            try:
                body = error.read().decode("utf-8", errors="replace")
            except (AttributeError, OSError):
                body = ""
            raw_response = {"http_error_body": body}
            error_text = f"Europe PMC HTTP {error.code} for query {query_text!r}"
            if error.code == 429:
                status = "rate_limited"
            elif 400 <= error.code < 500:
                status = "invalid_request"
            else:
                status = "provider_unavailable"
        except (urllib.error.URLError, TimeoutError, OSError) as error:
            status = "provider_unavailable"
            error_text = f"Europe PMC connection failed for query {query_text!r}: {error}"
        except ResponseParseError as error:
            status = "response_parse_error"
            raw_response = {"unparsed_response_text": error.raw_text}
            error_text = str(error)
        except ValueError as error:
            status = "response_parse_error"
            raw_response = {"parse_error": str(error)}
            error_text = str(error)
        if (
            status not in {"rate_limited", "provider_unavailable"}
            or attempt >= max_attempts
        ):
            break
        sleep(min(2 ** (attempt - 1), 4))

    cassette = _cassette_payload(
        request_payload,
        status=status,
        retrieved_at=retrieved_at,
        http_status=http_status,
        raw_response=raw_response,
        attempt_count=attempts,
        error=error_text,
    )
    cassette_path = write_dir / f"{key}.json"
    _atomic_write(cassette_path, cassette)
    loaded = _load_cassette(cassette_path, request_payload)
    return LiteratureRetrieval(
        **{
            **loaded.__dict__,
            "execution_source": "live_then_persisted",
        }
    )


def write_query_provenance(
    retrievals: list[LiteratureRetrieval],
    output_arg: str | Path,
) -> str:
    import pandas as pd

    columns = (
        "query",
        "request_key_sha256",
        "retrieval_status",
        "article_count",
        "provider",
        "retrieval_source",
        "provider_retrieved_at",
        "raw_response_sha256",
        "cassette_path",
        "http_status",
        "attempt_count",
        "error",
        "legacy_cache",
    )
    frame = pd.DataFrame(
        [retrieval.provenance_row() for retrieval in retrievals],
        columns=columns,
    )
    if not frame.empty:
        frame = frame.sort_values(
            ["request_key_sha256", "query"], kind="stable"
        ).reset_index(drop=True)
    output = runtime_path(str(output_arg))
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output, sep="\t", index=False, lineterminator="\n")
    return str(output)
