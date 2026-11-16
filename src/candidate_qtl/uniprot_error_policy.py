"""Typed UniProt entry-fetch outcomes for Route A.

A submitted identifier that is syntactically invalid, or that UniProt rejects
with HTTP 400, is an input problem. It is not provider unavailability. HTTP 404
is a definitive not-found result. Rate limits, 5xx responses, timeouts and
connection failures remain provider-unavailable/not-evaluated (issue #655).
"""

from __future__ import annotations

import contextlib
import contextvars
import json
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from typing import Iterator

from .identifiers import canonical_accession, is_uniprot_accession

UNIPROT_ENTRY_ENDPOINT = "https://rest.uniprot.org/uniprotkb/"

INVALID_ACCESSION_STATES = frozenset({"invalid_syntax", "invalid_http_400"})
PROVIDER_FAILURE_STATES = frozenset(
    {"rate_limited", "provider_unavailable", "response_parse_error"}
)


@dataclass(frozen=True)
class UniProtEntryEvent:
    accession: str
    state: str
    http_status: int | None
    network_attempted: bool
    detail: str


_CURRENT_EVENTS: contextvars.ContextVar[list[UniProtEntryEvent] | None] = (
    contextvars.ContextVar("candidate_qtl_uniprot_entry_events", default=None)
)


def _event_key(accession: object) -> str:
    value = str(accession or "").strip()
    canonical = canonical_accession(value)
    return (canonical or value).casefold()


def _record(
    accession: str,
    state: str,
    *,
    http_status: int | None,
    network_attempted: bool,
    detail: str,
) -> None:
    events = _CURRENT_EVENTS.get()
    if events is not None:
        events.append(
            UniProtEntryEvent(
                accession=str(accession),
                state=state,
                http_status=http_status,
                network_attempted=network_attempted,
                detail=detail,
            )
        )


@contextlib.contextmanager
def capture_uniprot_entry_events() -> Iterator[list[UniProtEntryEvent]]:
    events: list[UniProtEntryEvent] = []
    token = _CURRENT_EVENTS.set(events)
    try:
        yield events
    finally:
        _CURRENT_EVENTS.reset(token)


def invalid_events_by_accession(
    events: list[UniProtEntryEvent],
) -> dict[str, UniProtEntryEvent]:
    result: dict[str, UniProtEntryEvent] = {}
    for event in events:
        if event.state in INVALID_ACCESSION_STATES:
            result[_event_key(event.accession)] = event
    return result


def fetch_uniprot_entry(
    accession: str,
    timeout: int = 60,
) -> dict[str, object] | None:
    """Fetch one UniProt entry with fail-closed error taxonomy.

    Invalid syntax is rejected before any request. HTTP 400 is returned as an
    invalid submitted identifier event, HTTP 404 as a normal not-found result,
    and only transient/provider failures raise ``RuntimeError``.
    """
    submitted = str(accession or "").strip()
    normalised = canonical_accession(submitted)
    if not normalised or not is_uniprot_accession(normalised):
        _record(
            submitted,
            "invalid_syntax",
            http_status=None,
            network_attempted=False,
            detail="submitted value is not a canonical UniProt accession",
        )
        return None

    url = UNIPROT_ENTRY_ENDPOINT + urllib.parse.quote(normalised) + ".json"
    request = urllib.request.Request(
        url, headers={"User-Agent": "candidate-qtl hypothesis-verify"}
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as error:
        if error.code == 400:
            _record(
                normalised,
                "invalid_http_400",
                http_status=400,
                network_attempted=True,
                detail="UniProt rejected the submitted accession request",
            )
            return None
        if error.code == 404:
            _record(
                normalised,
                "not_found_http_404",
                http_status=404,
                network_attempted=True,
                detail="UniProt has no entry for the submitted accession",
            )
            return None
        state = "rate_limited" if error.code == 429 else "provider_unavailable"
        _record(
            normalised,
            state,
            http_status=int(error.code),
            network_attempted=True,
            detail=f"UniProt API HTTP {error.code}",
        )
        raise RuntimeError(
            f"UniProt API HTTP {error.code} for {normalised}"
        ) from error
    except (urllib.error.URLError, TimeoutError, OSError) as error:
        reason = getattr(error, "reason", error)
        _record(
            normalised,
            "provider_unavailable",
            http_status=None,
            network_attempted=True,
            detail=str(reason),
        )
        raise RuntimeError(
            f"UniProt API connection failed for {normalised}: {reason}"
        ) from error
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        _record(
            normalised,
            "response_parse_error",
            http_status=None,
            network_attempted=True,
            detail=str(error),
        )
        raise RuntimeError(
            f"UniProt API returned invalid JSON for {normalised}"
        ) from error

    if (
        not isinstance(payload, dict)
        or not isinstance(payload.get("primaryAccession"), str)
        or not isinstance(payload.get("organism"), dict)
        or not isinstance(payload.get("sequence"), dict)
    ):
        _record(
            normalised,
            "response_parse_error",
            http_status=None,
            network_attempted=True,
            detail="record payload lacks required UniProt fields",
        )
        raise RuntimeError(
            f"UniProt API returned an invalid record payload for {normalised}"
        )
    returned = canonical_accession(str(payload["primaryAccession"]))
    if returned != normalised:
        _record(
            normalised,
            "response_parse_error",
            http_status=None,
            network_attempted=True,
            detail=f"requested {normalised}, received {returned}",
        )
        raise RuntimeError(
            f"UniProt API returned accession {returned} for request {normalised}"
        )
    _record(
        normalised,
        "retrieved",
        http_status=200,
        network_attempted=True,
        detail="",
    )
    return payload
