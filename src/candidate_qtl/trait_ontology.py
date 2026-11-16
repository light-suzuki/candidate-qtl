"""Map free trait terms onto ontology IDs using only local tables.

This importer never calls an external ontology service. It joins the analyst's
trait terms against a local synonym/ontology table (Plant Trait Ontology, Plant
Ontology, Crop Ontology, or a hand-curated sheet) and preserves the original
trait expression alongside the mapped ID. Every mapping keeps its match type
(exact / synonym / manual / llm_suggested), an auto-match confidence, and a
human approval status so that unreviewed or model-suggested synonyms cannot
silently drive scoring.
"""

from __future__ import annotations


import pandas as pd

from .core import nonempty, runtime_path, sha256


TRAIT_TERM_COLUMNS = {"term"}
ONTOLOGY_MAP_COLUMNS = {"term", "ontology_id"}

VALID_MATCH_TYPES = {"exact", "synonym", "manual", "llm_suggested"}
DEFAULT_MATCH_CONFIDENCE = {
    "exact": "high",
    "manual": "high",
    "synonym": "medium",
    "llm_suggested": "low",
}

OUTPUT_COLUMNS = [
    "original_term",
    "term_type",
    "ontology",
    "ontology_id",
    "ontology_name",
    "match_type",
    "match_confidence",
    "approval_status",
    "source_id",
    "notes",
]


def _bool(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "approved", "y"}


def build_ontology_index(path_arg: str) -> dict[str, list[dict[str, str]]]:
    """Load a local ontology/synonym table keyed by casefolded term."""
    mapping = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    missing = sorted(ONTOLOGY_MAP_COLUMNS.difference(mapping.columns))
    if missing:
        raise ValueError(f"ontology map is missing required columns: {', '.join(missing)}")
    index: dict[str, list[dict[str, str]]] = {}
    for _, row in mapping.iterrows():
        key = str(row["term"]).strip().casefold()
        index.setdefault(key, []).append({str(k): str(v) for k, v in row.items()})
    return index


def map_trait_terms(
    trait_terms_arg: str,
    ontology_map_arg: str,
    out_arg: str,
) -> dict[str, object]:
    """Join trait terms to a local ontology/synonym table.

    Unmapped terms are retained with ``match_type=unmapped`` so a crop without
    ontology coverage never halts the workflow.
    """

    terms_path = runtime_path(trait_terms_arg)
    terms = pd.read_csv(terms_path, sep="\t", dtype=str).fillna("")
    missing = sorted(TRAIT_TERM_COLUMNS.difference(terms.columns))
    if missing:
        raise ValueError(f"trait terms table is missing required columns: {', '.join(missing)}")
    terms = terms[terms["term"].map(nonempty)]

    map_path = runtime_path(ontology_map_arg)
    mapping = pd.read_csv(map_path, sep="\t", dtype=str).fillna("")
    missing = sorted(ONTOLOGY_MAP_COLUMNS.difference(mapping.columns))
    if missing:
        raise ValueError(f"ontology map is missing required columns: {', '.join(missing)}")
    bad_types = set(mapping.get("match_type", pd.Series(dtype=str)).str.strip()).difference(VALID_MATCH_TYPES | {""})
    if bad_types:
        raise ValueError(f"ontology map has unknown match_type values: {', '.join(sorted(bad_types))}")

    index: dict[str, list[dict[str, str]]] = {}
    for _, row in mapping.iterrows():
        key = str(row["term"]).strip().casefold()
        index.setdefault(key, []).append({str(k): str(v) for k, v in row.items()})

    records: list[dict[str, str]] = []
    for _, term_row in terms.iterrows():
        original = str(term_row["term"]).strip()
        term_type = str(term_row.get("term_type", ""))
        term_approved = _bool(term_row.get("approved", ""))
        term_source = str(term_row.get("source_id", "")) or str(terms_path)
        notes = str(term_row.get("notes", ""))
        matches = index.get(original.casefold(), [])
        if not matches:
            records.append(
                {
                    "original_term": original,
                    "term_type": term_type,
                    "ontology": "",
                    "ontology_id": "",
                    "ontology_name": "",
                    "match_type": "unmapped",
                    "match_confidence": "",
                    "approval_status": "approved" if term_approved else "pending_review",
                    "source_id": term_source,
                    "notes": notes or "No local ontology entry; original term retained.",
                }
            )
            continue
        for match in matches:
            match_type = (match.get("match_type", "") or "exact").strip()
            confidence = match.get("confidence", "").strip() or DEFAULT_MATCH_CONFIDENCE.get(match_type, "low")
            # A model-suggested synonym stays pending until a human approves it,
            # regardless of the trait term's own approved flag.
            map_approved = _bool(match.get("approved", ""))
            if match_type == "llm_suggested":
                approval = "approved" if map_approved else "pending_review"
            else:
                approval = "approved" if (term_approved or map_approved) else "pending_review"
            records.append(
                {
                    "original_term": original,
                    "term_type": term_type,
                    "ontology": match.get("ontology", ""),
                    "ontology_id": match.get("ontology_id", ""),
                    "ontology_name": match.get("ontology_name", ""),
                    "match_type": match_type,
                    "match_confidence": confidence,
                    "approval_status": approval,
                    "source_id": match.get("source_id", "") or term_source,
                    "notes": notes,
                }
            )

    result = pd.DataFrame(records, columns=OUTPUT_COLUMNS)
    out_path = runtime_path(out_arg)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(out_path, sep="\t", index=False)

    manifest = {
        "trait_ontology_map": str(out_path),
        "term_count": int(terms.shape[0]),
        "mapped_term_count": int((result["match_type"] != "unmapped").groupby(result["original_term"]).any().sum()),
        "unmapped_term_count": int(result[result["match_type"] == "unmapped"]["original_term"].nunique()),
        "approved_count": int((result["approval_status"] == "approved").sum()),
        "pending_review_count": int((result["approval_status"] == "pending_review").sum()),
        "match_type_counts": {str(k): int(v) for k, v in result["match_type"].value_counts().items()},
        "input_fingerprints": [
            {"path": str(terms_path), "sha256": sha256(terms_path)},
            {"path": str(map_path), "sha256": sha256(map_path)},
        ],
        "principle": (
            "Ontology mapping is local-only; original trait terms are preserved and "
            "llm_suggested or unapproved synonyms stay pending_review until a human approves them."
        ),
    }
    return manifest
