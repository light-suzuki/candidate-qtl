"""Lossless bibliography references; scientific observations are never shared."""
from __future__ import annotations

import copy
import json


def expand_paper_observations(value: object, bibliography: dict | None = None, *,
                              id_field: str = "source_id", title_field: str = "source_title") -> object:
    """Expand optional bibliography/observations envelopes before existing checks.

    Catalogs contain only the LLM's paper ID and title. Gene, phenotype, context,
    correspondence and checklist stay on every observation, including reviews.
    No alias inference, cross-route cache, scientific merge or default answers.
    """
    if value is None or value == "":
        return value
    decoded = json.loads(value) if isinstance(value, str) else value
    catalog = bibliography or {}
    if isinstance(decoded, dict):
        if set(decoded) != {"bibliography", "observations"}:
            raise ValueError("paper envelope requires only bibliography and observations")
        if bibliography:
            raise ValueError("paper envelope cannot override an external bibliography")
        catalog, decoded = decoded["bibliography"], decoded["observations"]
    if not isinstance(catalog, dict):
        raise ValueError("bibliography must be an object")
    for key, entry in catalog.items():
        if (not isinstance(key, str) or not key.strip() or not isinstance(entry, dict)
                or set(entry) != {"source_id", "source_title"}
                or any(not isinstance(v, str) or not v.strip() for v in entry.values())):
            raise ValueError("bibliography entries require only nonempty source_id and source_title")
    if not isinstance(decoded, list):
        raise ValueError("paper observations must be an array")
    result = []
    for observation in decoded:
        if not isinstance(observation, dict):
            raise ValueError("each paper observation must be an object")
        row = copy.deepcopy(observation)
        if "paper_ref" in row:
            key = row.pop("paper_ref")
            if not isinstance(key, str) or key not in catalog:
                raise ValueError("unknown paper_ref in this payload bibliography")
            entry = catalog[key]
            for target, source in ((id_field, "source_id"), (title_field, "source_title")):
                if target in row and row[target] != entry[source]:
                    raise ValueError("paper_ref conflicts with the observation's bibliography")
                row[target] = entry[source]
        result.append(row)
    return result
