"""Build a controlled, approval-gated trait-term set for keyword evidence.

Free-text keywords are easy to leak answers with (once the causal gene is
known) and easy to over-match with generic words. This builder turns raw trait
terms plus a local ontology table, optional human synonyms, and optional
model-proposed synonyms into a single audited set. Only human-approved,
non-leaking terms reach the approved output that ``annotation-keywords``
consumes; the free ``--keyword`` path stays available as a manual escape hatch.
"""

from __future__ import annotations

import json

import pandas as pd

from .core import nonempty, runtime_path, sha256
from .trait_ontology import build_ontology_index


TERM_CLASSES = {"direct_phenotype", "organ_stage", "biological_process", "mechanism_hypothesis"}
DIRECT_CLASSES = {"direct_phenotype"}
MATCH_METHODS = {"exact_phrase", "token", "stem"}

APPROVED_COLUMNS = [
    "term",
    "term_type",
    "weight_class",
    "match_method",
    "match_type",
    "ontology_id",
    "approval_status",
    "source_id",
]
AUDIT_COLUMNS = APPROVED_COLUMNS + ["included", "exclusion_reason"]


def _bool(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "approved", "y"}


def _weight_class(term_type: str, explicit: str) -> str:
    if explicit in {"direct", "context"}:
        return explicit
    return "direct" if term_type in DIRECT_CLASSES else "context"


def _normalize(text: str) -> str:
    return "".join(ch for ch in text.strip().casefold() if ch.isalnum() or ch.isspace()).strip()


def _leaks_answer(term: str, blocked: list[str]) -> str | None:
    term_norm = _normalize(term)
    if not term_norm:
        return None
    for name in blocked:
        name_norm = _normalize(name)
        if not name_norm:
            continue
        if name_norm == term_norm or name_norm in term_norm.split() or name_norm in term_norm:
            return f"answer_leakage:{name}"
    return None


def build_trait_terms(
    trait_terms_arg: str,
    out_arg: str,
    ontology_map_arg: str | None = None,
    manual_synonyms_arg: str | None = None,
    llm_proposals_arg: str | None = None,
    blocked_terms: list[str] | None = None,
    benchmark: bool = False,
) -> dict[str, object]:
    terms_path = runtime_path(trait_terms_arg)
    terms = pd.read_csv(terms_path, sep="\t", dtype=str).fillna("")
    if "term" not in terms.columns:
        raise ValueError("trait terms table requires a term column")
    terms = terms[terms["term"].map(nonempty)]

    bad_classes = set(terms.get("term_type", pd.Series(dtype=str)).str.strip()).difference(TERM_CLASSES | {""})
    if bad_classes:
        raise ValueError(f"trait terms have unknown term_type values: {', '.join(sorted(bad_classes))}")

    ontology_index = build_ontology_index(ontology_map_arg) if ontology_map_arg else {}

    def ontology_id_for(term: str) -> str:
        matches = ontology_index.get(term.strip().casefold(), [])
        return matches[0].get("ontology_id", "") if matches else ""

    collected: list[dict[str, str]] = []

    for _, row in terms.iterrows():
        term = str(row["term"]).strip()
        term_type = str(row.get("term_type", "")).strip()
        match_method = str(row.get("match_method", "")).strip() or "exact_phrase"
        if match_method not in MATCH_METHODS:
            raise ValueError(f"unknown match_method {match_method!r}; use one of {sorted(MATCH_METHODS)}")
        collected.append(
            {
                "term": term,
                "term_type": term_type,
                "weight_class": _weight_class(term_type, str(row.get("weight_class", "")).strip()),
                "match_method": match_method,
                "match_type": "curated",
                "ontology_id": ontology_id_for(term),
                "approval_status": "approved" if _bool(row.get("approved", "")) else "pending_review",
                "source_id": str(row.get("source_id", "")) or str(terms_path),
            }
        )

    parent_type = {str(row["term"]).strip().casefold(): str(row.get("term_type", "")).strip() for _, row in terms.iterrows()}

    if manual_synonyms_arg:
        synonyms = pd.read_csv(runtime_path(manual_synonyms_arg), sep="\t", dtype=str).fillna("")
        if not {"term", "synonym"}.issubset(synonyms.columns):
            raise ValueError("manual synonyms table requires term and synonym columns")
        for _, row in synonyms.iterrows():
            synonym = str(row["synonym"]).strip()
            if not synonym:
                continue
            term_type = str(row.get("term_type", "")).strip() or parent_type.get(str(row["term"]).strip().casefold(), "")
            collected.append(
                {
                    "term": synonym,
                    "term_type": term_type,
                    "weight_class": _weight_class(term_type, str(row.get("weight_class", "")).strip()),
                    "match_method": str(row.get("match_method", "")).strip() or "exact_phrase",
                    "match_type": "manual",
                    "ontology_id": ontology_id_for(synonym),
                    # A human curated this synonym table, so it is approved.
                    "approval_status": "approved",
                    "source_id": str(row.get("source_id", "")) or str(runtime_path(manual_synonyms_arg)),
                }
            )

    if llm_proposals_arg:
        proposals = json.loads(runtime_path(llm_proposals_arg).read_text(encoding="utf-8"))
        for proposal in proposals:
            synonym = str(proposal.get("synonym") or proposal.get("term") or "").strip()
            if not synonym:
                continue
            term_type = str(proposal.get("term_type", "")).strip() or parent_type.get(str(proposal.get("term", "")).strip().casefold(), "")
            collected.append(
                {
                    "term": synonym,
                    "term_type": term_type,
                    "weight_class": _weight_class(term_type, str(proposal.get("weight_class", "")).strip()),
                    "match_method": str(proposal.get("match_method", "")).strip() or "exact_phrase",
                    "match_type": "llm_suggested",
                    "ontology_id": str(proposal.get("ontology_id", "")) or ontology_id_for(synonym),
                    # Model proposals are never auto-approved.
                    "approval_status": "pending_review",
                    "source_id": str(proposal.get("source_id", "")) or str(runtime_path(llm_proposals_arg)),
                }
            )

    blocked = [value for value in (blocked_terms or []) if nonempty(value)]
    audit: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for entry in collected:
        key = (entry["term"].casefold(), entry["match_type"])
        if key in seen:
            continue
        seen.add(key)
        if entry["approval_status"] != "approved":
            reason = "not_approved"
        elif benchmark and (leak := _leaks_answer(entry["term"], blocked)):
            reason = leak
        else:
            reason = ""
        included = reason == ""
        audit.append({**entry, "included": "yes" if included else "no", "exclusion_reason": reason})

    audit_frame = pd.DataFrame(audit, columns=AUDIT_COLUMNS)
    approved_frame = audit_frame[audit_frame["included"] == "yes"][APPROVED_COLUMNS].reset_index(drop=True)

    out_path = runtime_path(out_arg)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    approved_frame.to_csv(out_path, sep="\t", index=False)
    audit_path = out_path.with_name(out_path.stem + "_audit.tsv")
    audit_frame.to_csv(audit_path, sep="\t", index=False)

    direct_keywords = sorted(set(approved_frame[approved_frame["weight_class"] == "direct"]["term"]))
    context_keywords = sorted(set(approved_frame[approved_frame["weight_class"] == "context"]["term"]))
    rejected_leakage = sorted(
        set(audit_frame[audit_frame["exclusion_reason"].str.startswith("answer_leakage")]["term"])
    )

    manifest = {
        "approved_terms": str(out_path),
        "audit": str(audit_path),
        "approved_count": int(len(approved_frame)),
        "pending_review_count": int((audit_frame["approval_status"] == "pending_review").sum()),
        "rejected_answer_leakage": rejected_leakage,
        "direct_keywords": direct_keywords,
        "context_keywords": context_keywords,
        "benchmark_mode": benchmark,
        "input_fingerprints": [
            {"path": str(path), "sha256": sha256(path)}
            for path in [terms_path]
        ],
        "principle": (
            "Only human-approved, non-leaking terms are exported; llm_suggested and unapproved "
            "terms stay in the audit table. Keyword matches remain a low auxiliary functional signal."
        ),
    }
    return manifest
