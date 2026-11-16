"""Split citation lists without changing punctuation inside a DOI suffix."""
from __future__ import annotations

import re

_DOI_START = re.compile(r"^(?:(?:doi\s*:\s*)|https?://(?:dx\.)?doi\.org/)?10\.\d{4,9}/", re.I)
_EXPLICIT_START = re.compile(
    r"^(?:(?:pmid|pubmed|medline|med|pmcid|pmc|doi)\s*[: ]|PMC\d{4,9}\b|10\.\d{4,9}/"
    r"|https?://(?:(?:dx\.)?doi\.org/|pubmed\.ncbi\.nlm\.nih\.gov/|pmc\.ncbi\.nlm\.nih\.gov/))", re.I)


def split_citations(value: object) -> list[str]:
    """Keep DOI semicolons/commas unless another explicit citation starts.

    Pipe remains an explicit list separator. A DOI followed by an ambiguous
    bare number is kept intact; use a PMID prefix to express that boundary.
    Other citation types retain the legacy semicolon/comma/pipe separators.
    """
    if isinstance(value, (list, tuple)):
        return [token for item in value for token in split_citations(item)]
    text = str(value or "")
    result, start = [], 0
    for match in re.finditer(r"[;,|]", text):
        token = text[start:match.start()].strip()
        following = text[match.end():].lstrip()
        if _DOI_START.match(token) and match.group() != "|" and not _EXPLICIT_START.match(following):
            continue
        if token:
            result.append(token)
        start = match.end()
    if text[start:].strip():
        result.append(text[start:].strip())
    return result
