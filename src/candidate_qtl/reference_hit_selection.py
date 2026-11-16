"""Truth-independent reference-hit selection within a narrow homology window.

Sequence evidence remains primary. Annotation quality may only break a record
choice when the alternative stays inside predeclared identity, coverage and
bitscore bounds relative to the strongest sequence hit. The selected record and
the true sequence-top record are returned separately for provenance.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

MIN_EQUIVALENT_BITSCORE_RATIO = 0.95
MAX_IDENTITY_DEFICIT_PERCENTAGE_POINTS = 1.0
MAX_QUERY_COVERAGE_DEFICIT = 0.02
MAX_SUBJECT_COVERAGE_DEFICIT = 0.02

_GENERIC_TITLE_MARKERS = (
    "uncharacterized protein",
    "hypothetical protein",
    "unknown protein",
)
_EPSILON = 1e-12


def is_reviewed_subject(subject: object) -> bool:
    """Return True only for an explicit UniProtKB/Swiss-Prot ``sp|`` header."""
    return str(subject or "").strip().casefold().startswith("sp|")


def is_informative_title(title: object) -> bool:
    """Conservatively reject only conventional no-function placeholders.

    ``putative`` by itself is intentionally retained: a title such as
    ``Cyclops protein, putative`` still communicates an identity and must not be
    treated like an uncharacterized record.
    """
    text = " ".join(str(title or "").split()).casefold()
    if not text:
        return False
    return not any(marker in text for marker in _GENERIC_TITLE_MARKERS)


def _number(record: Mapping[str, Any], key: str) -> float:
    try:
        return float(record.get(key, 0.0) or 0.0)
    except (TypeError, ValueError):
        return 0.0


def _query_coverage(record: Mapping[str, Any]) -> float:
    if "merged_query_coverage" in record:
        return _number(record, "merged_query_coverage")
    return _number(record, "query_coverage")


def _subject_coverage(record: Mapping[str, Any]) -> float:
    if "subject_coverage" in record:
        return _number(record, "subject_coverage")
    length = _number(record, "subject_length")
    if length <= 0:
        return 0.0
    start = _number(record, "subject_start")
    end = _number(record, "subject_end")
    return min(1.0, (abs(end - start) + 1.0) / length)


def _subject(record: Mapping[str, Any]) -> str:
    return str(record.get("subject_id", record.get("subject", "")) or "")


def _sequence_sort_key(record: Mapping[str, Any]) -> tuple[float, float, str]:
    return (-_number(record, "bitscore"), _number(record, "evalue"), _subject(record))


def within_homology_equivalence_window(
    candidate: Mapping[str, Any],
    sequence_top: Mapping[str, Any],
) -> bool:
    """Return whether annotation may legitimately break this sequence choice."""
    top_bitscore = _number(sequence_top, "bitscore")
    candidate_bitscore = _number(candidate, "bitscore")
    if top_bitscore > 0 and candidate_bitscore / top_bitscore + _EPSILON < MIN_EQUIVALENT_BITSCORE_RATIO:
        return False
    if _number(candidate, "identity") + MAX_IDENTITY_DEFICIT_PERCENTAGE_POINTS + _EPSILON < _number(sequence_top, "identity"):
        return False
    if _query_coverage(candidate) + MAX_QUERY_COVERAGE_DEFICIT + _EPSILON < _query_coverage(sequence_top):
        return False
    top_subject_coverage = _subject_coverage(sequence_top)
    candidate_subject_coverage = _subject_coverage(candidate)
    if top_subject_coverage > 0 and candidate_subject_coverage + MAX_SUBJECT_COVERAGE_DEFICIT + _EPSILON < top_subject_coverage:
        return False
    return True


def _annotation_sort_key(record: Mapping[str, Any]) -> tuple[int, int, int, float, float, str]:
    usable_symbol = bool(record.get("usable_symbol", True))
    reviewed = bool(record.get("reviewed", is_reviewed_subject(_subject(record))))
    informative = is_informative_title(record.get("subject_title", record.get("title", "")))
    return (
        -int(usable_symbol),
        -int(reviewed),
        -int(informative),
        -_number(record, "bitscore"),
        _number(record, "evalue"),
        _subject(record),
    )


def select_reference_hit(records: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Select one record and preserve the exact sequence-top alternative."""
    if not records:
        raise ValueError("reference-hit selection requires at least one record")

    sequence_top_index = min(range(len(records)), key=lambda index: _sequence_sort_key(records[index]))
    sequence_top = records[sequence_top_index]
    eligible_indices = [
        index
        for index, record in enumerate(records)
        if within_homology_equivalence_window(record, sequence_top)
    ]
    selected_index = min(eligible_indices, key=lambda index: _annotation_sort_key(records[index]))
    selected = records[selected_index]

    if selected_index == sequence_top_index:
        reason = "sequence_top"
    else:
        top_reviewed = bool(sequence_top.get("reviewed", is_reviewed_subject(_subject(sequence_top))))
        selected_reviewed = bool(selected.get("reviewed", is_reviewed_subject(_subject(selected))))
        top_informative = is_informative_title(sequence_top.get("subject_title", sequence_top.get("title", "")))
        selected_informative = is_informative_title(selected.get("subject_title", selected.get("title", "")))
        if selected_reviewed and not top_reviewed:
            reason = "reviewed_record_within_homology_equivalence_window"
        elif selected_informative and not top_informative:
            reason = "informative_title_within_homology_equivalence_window"
        else:
            reason = "annotation_tiebreak_within_homology_equivalence_window"

    return {
        "sequence_top_index": sequence_top_index,
        "selected_index": selected_index,
        "sequence_top": dict(sequence_top),
        "selected": dict(selected),
        "selection_reason": reason,
        "equivalent_hit_count": len(eligible_indices),
        "policy": "sequence first; annotation quality only within fixed bitscore/identity/coverage bounds",
        "bounds": {
            "min_bitscore_ratio": MIN_EQUIVALENT_BITSCORE_RATIO,
            "max_identity_deficit_percentage_points": MAX_IDENTITY_DEFICIT_PERCENTAGE_POINTS,
            "max_query_coverage_deficit": MAX_QUERY_COVERAGE_DEFICIT,
            "max_subject_coverage_deficit": MAX_SUBJECT_COVERAGE_DEFICIT,
        },
    }
