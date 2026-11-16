"""Public Route-B literature facade with immutable cassette retrieval.

The legacy helper functions remain exported for compatibility and request-level
tests. The production ``reverse_literature`` entry point uses the shared
cassette implementation introduced for issue #662.
"""

from __future__ import annotations

import json
# Existing tests patch ``candidate_qtl.reverse_literature.urllib.request`` while
# calling the legacy low-level fetch helper. Re-exporting the shared urllib
# module preserves that public test seam without keeping the old pipeline live.
import urllib
from pathlib import Path

from .core import sha256
from .reverse_literature_cassette import (
    QUERY_STATES,
    SCHEMA_VERSION,
    reverse_literature as _reverse_literature,
)
from .reverse_literature_legacy import (
    EUROPE_PMC_RESULT_TYPE,
    EUROPE_PMC_SEARCH,
    LITERATURE_COLUMNS,
    STATUS_COLUMNS,
    USER_AGENT,
    _fetch,
    _load_cache,
    _phenotype_terms,
    _query_for,
    _searchable_reference,
)


def _provider_summary(manifest: dict[str, object]) -> str:
    """Keep the historical run-level provider label when it is unambiguous.

    Exact mixed-source provenance remains in ``retrieval_source_counts`` and the
    query-provenance TSV.  Existing callers, however, use the top-level provider
    field to distinguish offline, local-cache, and Europe-PMC runs.
    """
    raw_counts = manifest.get("retrieval_source_counts", {})
    sources = {
        str(source)
        for source, count in raw_counts.items()
        if isinstance(raw_counts, dict) and int(count) > 0
    }
    if not sources or sources <= {"offline", "validation"}:
        return "none"
    if sources <= {"legacy_query_cache"}:
        return "local_query_cache"
    if sources <= {"live_europe_pmc"}:
        return "europe_pmc_rest"
    return "per_query"


def reverse_literature(
    lookup_arg: str,
    out_dir_arg: str,
    phenotype: str = "",
    cache_arg: str | None = None,
    allow_external_api: bool = False,
    limit: int = 5,
    case_id: str = "",
    full_hash: bool = False,
) -> dict[str, object]:
    """Run cassette retrieval while preserving the public manifest contract."""
    result = _reverse_literature(
        lookup_arg,
        out_dir_arg,
        phenotype=phenotype,
        cache_arg=cache_arg,
        allow_external_api=allow_external_api,
        limit=limit,
        case_id=case_id,
        full_hash=full_hash,
    )
    manifest_path = Path(str(result["manifest"]))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict):
        raise ValueError("reverse-literature manifest must be a JSON object")
    manifest["provider"] = _provider_summary(manifest)
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    run_envelope_value = result.get("run_envelope")
    if run_envelope_value:
        run_envelope_path = Path(str(run_envelope_value))
        run_envelope = json.loads(run_envelope_path.read_text(encoding="utf-8"))
        if not isinstance(run_envelope, dict):
            raise ValueError("reverse-literature run envelope must be a JSON object")
        run_envelope["scientific_manifest_sha256"] = sha256(manifest_path)
        run_envelope_path.write_text(
            json.dumps(run_envelope, ensure_ascii=False, indent=2, sort_keys=True)
            + "\n",
            encoding="utf-8",
        )
    return result


__all__ = [
    "EUROPE_PMC_RESULT_TYPE",
    "EUROPE_PMC_SEARCH",
    "LITERATURE_COLUMNS",
    "QUERY_STATES",
    "SCHEMA_VERSION",
    "STATUS_COLUMNS",
    "USER_AGENT",
    "_fetch",
    "_load_cache",
    "_phenotype_terms",
    "_query_for",
    "_searchable_reference",
    "reverse_literature",
    "urllib",
]
