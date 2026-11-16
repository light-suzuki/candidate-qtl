"""Lossless source snapshots for review, never a claim-support producer."""

from __future__ import annotations

import hashlib
import json


def content_sha256(value: object) -> str:
    return hashlib.sha256(json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":"),
    ).encode("utf-8")).hexdigest()


def capture_source(
    record: dict[str, object], *, citation: str, citation_type: str,
    resolved_id: str, classification_scope: str,
) -> dict[str, object]:
    # Keep original values, whitespace and markup; normalised bibliographic
    # fields and concatenated classifier text cannot reconstruct quote scope.
    sections = {key: record[key] for key in (
        "title", "abstractText", "fullText", "journalTitle", "keywordList",
    ) if key in record}
    def has_text(key: str) -> bool:
        value = sections.get(key)
        return isinstance(value, str) and bool(value.strip())
    available_scope = (
        "full_text_cached" if has_text("fullText") else
        "title_abstract" if has_text("abstractText") else
        "title_only" if has_text("title") else "none"
    )
    return {
        "citation": citation,
        "citation_type": citation_type,
        "source_record_id": resolved_id,
        "source_identifiers": {key: record[key] for key in (
            "id", "pmid", "pmcid", "doi", "source",
        ) if key in record},
        "sections": sections,
        "sections_sha256": content_sha256(sections),
        "provider_record_sha256": content_sha256(record),
        "available_text_scope": available_scope,
        "classification_scope": classification_scope,
        "hypothesis_ids": [],
        "semantic_binding_state": "not_evaluated",
    }
