from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd

from .core import nonempty, parse_interval, runtime_path
from .evidence_io import evidence_row as _row


# These are ordinal assay-prioritization priors, not causal probabilities.
# The profile deliberately keeps family-level homology/function weak, gives
# position meaningful weight only when a peak is supplied, and reserves the
# largest values for genetic evidence. It must still be calibrated on known
# QTGs before being interpreted across species.
LITERATURE_INFORMED_WEIGHTS = {
    "direct_model_genetic": 5.0,
    "known_qtg": 4.0,
    "ortholog_model_genetic": 3.0,
    "recombinant_core": 2.0,
    "peak_overlap": 2.0,
    "peak_proximal": 1.25,
    "support_interval": 0.5,
    "rbh_ortholog": 1.25,
    "curated_homolog_match": 1.25,
    "automated_homolog_match": 0.75,
    "automated_candidate_function_match": 1.25,
    "synteny_supported": 1.0,
    "expression_target_specific": 1.0,
    "orthogroup_supported": 0.75,
    "domain_supported": 0.75,
    "trait_function_match": 1.0,
    "paralog_ambiguous": 0.75,
    "counterevidence": 3.0,
    "expression_absent": 0.5,
    "homology_supported": 0.25,
    "expression_detected": 0.25,
    "boundary_proximal": 0.0,
}


# Evidence produced from the same comparative annotation chain is intentionally
# one family. BLAST + domain + annotation keyword must not masquerade as three
# independent lines of evidence.
INDEPENDENCE_FAMILIES = {
    "direct_model_genetic": "genetic",
    "ortholog_model_genetic": "genetic",
    "known_qtg": "genetic",
    "rbh_ortholog": "comparative_function",
    "curated_homolog_match": "comparative_function",
    "automated_homolog_match": "comparative_function",
    "automated_candidate_function_match": "comparative_function",
    "orthogroup_supported": "comparative_function",
    "homology_supported": "comparative_function",
    "paralog_ambiguous": "comparative_function",
    "synteny_supported": "comparative_function",
    "domain_supported": "comparative_function",
    "homology_context": "comparative_function",
    "domain_context": "comparative_function",
    "go_annotation_context": "comparative_function",
    "trait_function_match": "comparative_function",
    "trait_context": "comparative_function",
    "expression_target_specific": "expression",
    "expression_detected": "expression",
    "expression_absent": "expression",
    "expression_comparison_unavailable": "expression",
    "recombinant_core": "position",
    "peak_overlap": "position",
    "peak_proximal": "position",
    "peak_heuristic_proximal": "position",
    "support_interval": "position",
    "interval_membership": "position",
    "boundary_proximal": "position",
    "counterevidence": "counterevidence",
}

REGISTRY_VERSION = "1.1.0"
REGISTRY_COLUMNS = {
    "evidence_code", "weight_profile", "base_weight_or_tier", "default_direction",
    "allowed_confidence", "evidence_channel", "directness", "trait_specificity",
    "independence_rule", "rationale", "reference_type", "reference_id",
    "review_status", "registry_version",
}


def default_registry_path() -> Path:
    """Return the registry used by the default ranker.

    The authoritative, human-editable copy lives at the repo-root
    ``templates/evidence_weight_registry.tsv`` and is used from a source tree or
    editable install. A byte-identical copy is vendored into the package at
    ``candidate_qtl/resources/`` so a clean ``pip install`` of the wheel can
    still rank without the source tree present (a drift test keeps them in sync).
    """
    here = Path(__file__).resolve()
    source_tree = here.parents[2] / "templates" / "evidence_weight_registry.tsv"
    if source_tree.exists():
        return source_tree
    return here.parent / "resources" / "evidence_weight_registry.tsv"


def load_weight_registry(path_arg: str | None = None) -> tuple[pd.DataFrame, dict[str, str]]:
    path = runtime_path(path_arg) if path_arg else default_registry_path()
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    missing = sorted(REGISTRY_COLUMNS.difference(frame.columns))
    if missing:
        raise ValueError(f"weight registry is missing columns: {', '.join(missing)}")
    if frame["evidence_code"].duplicated().any():
        duplicates = sorted(frame.loc[frame["evidence_code"].duplicated(False), "evidence_code"].unique())
        raise ValueError(f"weight registry contains duplicate evidence_code values: {', '.join(duplicates)}")
    versions = sorted(set(frame["registry_version"]) - {""})
    if len(versions) != 1:
        raise ValueError("weight registry must contain exactly one non-empty registry_version")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return frame, {"path": str(path), "version": versions[0], "sha256": digest}


def load_weights(path_arg: str | None) -> dict[str, float]:
    if not path_arg:
        frame, _ = load_weight_registry()
        weighted = frame[frame["base_weight_or_tier"].str.strip().ne("")]
        return {str(row["evidence_code"]): float(row["base_weight_or_tier"]) for _, row in weighted.iterrows()}
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    if "base_weight_or_tier" in frame.columns:
        weighted = frame[frame["base_weight_or_tier"].str.strip().ne("")]
        return {str(row["evidence_code"]): float(row["base_weight_or_tier"]) for _, row in weighted.iterrows()}
    if not {"evidence_code", "weight"}.issubset(frame.columns):
        raise ValueError("weight table must contain evidence_code and weight, or be a full evidence weight registry")
    return {
        str(row["evidence_code"]): float(row["weight"])
        for _, row in frame.iterrows()
    }


def position_evidence(
    candidates: pd.DataFrame,
    peak: str | None,
    interval: str | None,
    peak_window: int,
    boundary_window: int,
) -> tuple[list[dict[str, str]], pd.DataFrame]:
    """Create positional evidence while retaining peak precision.

    An explicit --peak-window is treated as an analyst-supplied uncertainty or
    support half-width and is the only way proximity becomes scored support. A
    bare peak (no window) only scores exact overlap; the interval-length 5%
    heuristic is reported as context (peak_heuristic_proximal) but never scored,
    so mapping precision -- not interval width -- decides positional support.
    Merely being close to a broad QTL boundary is context, not positive evidence.
    """
    if not {"chrom", "start", "end"}.issubset(candidates.columns):
        return [], pd.DataFrame()

    inferred_interval = interval
    if not inferred_interval and "interval" in candidates.columns:
        values = [value for value in candidates["interval"].unique() if nonempty(value)]
        inferred_interval = values[0] if len(values) == 1 else None
    interval_data = parse_interval(inferred_interval) if inferred_interval else None

    peak_chrom = None
    peak_pos = None
    if peak:
        if ":" in peak:
            peak_chrom, peak_text = peak.split(":", 1)
            peak_pos = int(peak_text)
        else:
            peak_pos = int(peak)

    explicit_peak_window = peak_window > 0
    if interval_data and peak_window <= 0:
        peak_window = max(
            1000,
            int((interval_data[2] - interval_data[1] + 1) * 0.05),
        )
    peak_window_source = (
        "explicit_uncertainty_window"
        if explicit_peak_window
        else "automatic_5_percent_interval_heuristic"
    )

    rows: list[dict[str, str]] = []
    summaries: list[dict[str, object]] = []
    for _, candidate in candidates.iterrows():
        gene_id = str(candidate["gene_id"])
        chrom = str(candidate["chrom"])
        start = int(candidate["start"])
        end = int(candidate["end"])
        distance_peak: int | str = ""
        distance_boundary: int | str = ""
        positional_class = "interval_only"

        if peak_pos is not None and (not peak_chrom or peak_chrom == chrom):
            distance_peak = (
                0
                if start <= peak_pos <= end
                else min(abs(start - peak_pos), abs(end - peak_pos))
            )
            if int(distance_peak) == 0:
                positional_class = "peak_overlap"
                confidence = "high" if explicit_peak_window else "medium"
                rows.append(
                    _row(
                        gene_id,
                        "peak_overlap",
                        "support",
                        confidence,
                        distance_peak,
                        "qtl_position",
                        peak or "",
                        (
                            "Gene overlaps the reported QTL peak. Confidence reflects "
                            f"the supplied peak-window policy: {peak_window_source}."
                        ),
                    )
                )
            elif explicit_peak_window and peak_window > 0 and int(distance_peak) <= peak_window:
                positional_class = "peak_support_window"
                rows.append(
                    _row(
                        gene_id,
                        "peak_proximal",
                        "support",
                        "medium",
                        distance_peak,
                        "qtl_position",
                        peak or "",
                        (
                            f"Gene is within the analyst-supplied {peak_window} bp peak uncertainty window."
                        ),
                    )
                )
            elif not explicit_peak_window and peak_window > 0 and int(distance_peak) <= peak_window:
                # A peak with no explicit uncertainty only gives an exact-overlap
                # signal. Proximity under the interval-length heuristic is reported
                # as context (distance is auditable) but is never scored, so wide
                # QTLs do not hand out position support to more genes.
                positional_class = "peak_heuristic_context"
                rows.append(
                    _row(
                        gene_id,
                        "peak_heuristic_proximal",
                        "context",
                        "low",
                        distance_peak,
                        "qtl_position",
                        peak or "",
                        (
                            f"Gene is within {peak_window} bp of the peak under the {peak_window_source}; "
                            "context only. Supply --peak-window, a LOD support interval, or a recombinant core for scored proximity."
                        ),
                    )
                )

        if interval_data and interval_data[0] == chrom:
            distance_boundary = min(
                abs(start - interval_data[1]),
                abs(end - interval_data[2]),
            )
            if boundary_window > 0 and int(distance_boundary) <= boundary_window:
                rows.append(
                    _row(
                        gene_id,
                        "boundary_proximal",
                        "context",
                        "low",
                        distance_boundary,
                        "qtl_position",
                        inferred_interval or "",
                        "Near a mapped interval boundary; retained as context, not positive evidence.",
                    )
                )

        summaries.append(
            {
                "gene_id": gene_id,
                "distance_to_peak_bp": distance_peak,
                "distance_to_boundary_bp": distance_boundary,
                "position_class": positional_class,
                "peak_window_bp": peak_window if peak_pos is not None else "",
                "peak_window_source": peak_window_source if peak_pos is not None else "",
            }
        )
    return rows, pd.DataFrame(summaries)
