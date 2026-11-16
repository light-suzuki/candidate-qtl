from __future__ import annotations

import hashlib
import json

import pandas as pd

from .core import runtime_path


MAPPING_GATE_VERSION = "1.0.0"
DEFAULT_THRESHOLDS = {
    "min_markers_per_mb": 1.0,
    "min_informative_recombinants": 3,
    "low_recombination_bp_per_cm": 5_000_000,
    "broad_interval_bp": 20_000_000,
    "broad_interval_gene_count": 200,
}
REQUIRED_METRICS = {
    "interval_size_bp",
    "marker_count",
    "informative_recombinant_count",
    "lod_peak_count",
    "peak_stability",
    "phenotype_direction_consistent",
}


def _bool(value: object) -> bool | None:
    text = str(value).strip().casefold()
    if text in {"true", "1", "yes", "y"}:
        return True
    if text in {"false", "0", "no", "n"}:
        return False
    return None


def mapping_adequacy_gate(
    metrics_arg: str,
    out_dir_arg: str,
    threshold_arg: str | None = None,
) -> dict[str, object]:
    """Gate candidate triage on mapping adequacy without changing gene scores."""
    metrics_path = runtime_path(metrics_arg)
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(metrics_path, sep="\t", dtype=str).fillna("")
    if not {"metric", "value", "source"}.issubset(frame.columns):
        raise ValueError("mapping metrics require metric, value, and source columns")
    if frame["metric"].duplicated().any() or frame["metric"].eq("").any():
        raise ValueError("mapping metric names must be unique and non-empty")
    values = dict(zip(frame["metric"], frame["value"]))
    sources = dict(zip(frame["metric"], frame["source"]))
    thresholds = dict(DEFAULT_THRESHOLDS)
    threshold_fingerprint = None
    if threshold_arg:
        threshold_path = runtime_path(threshold_arg)
        supplied = pd.read_csv(threshold_path, sep="\t", dtype=str).fillna("")
        if not {"threshold", "value"}.issubset(supplied.columns):
            raise ValueError("threshold table requires threshold and value columns")
        for row in supplied.itertuples():
            if row.threshold not in thresholds:
                raise ValueError(f"unknown mapping threshold: {row.threshold}")
            thresholds[row.threshold] = float(row.value)
        threshold_fingerprint = {
            "path": str(threshold_path),
            "sha256": hashlib.sha256(threshold_path.read_bytes()).hexdigest(),
        }

    missing = sorted(metric for metric in REQUIRED_METRICS if not str(values.get(metric, "")).strip())
    statuses: list[str] = []
    warnings: list[str] = []
    review_rows: list[dict[str, str]] = []

    def number(name: str) -> float | None:
        if not str(values.get(name, "")).strip():
            return None
        try:
            return float(values[name])
        except ValueError as error:
            raise ValueError(f"mapping metric {name} must be numeric") from error

    interval_bp = number("interval_size_bp")
    marker_count = number("marker_count")
    recombinants = number("informative_recombinant_count")
    peak_count = number("lod_peak_count")
    candidate_count = number("candidate_gene_count")
    interval_cm = number("interval_size_cm")
    peak_stable = _bool(values.get("peak_stability", ""))
    direction_consistent = _bool(values.get("phenotype_direction_consistent", ""))
    if values.get("peak_stability", "") and peak_stable is None:
        raise ValueError("peak_stability must be true or false")
    if values.get("phenotype_direction_consistent", "") and direction_consistent is None:
        raise ValueError("phenotype_direction_consistent must be true or false")

    if missing:
        statuses.append("mapping_model_review_required")
        for metric in missing:
            review_rows.append({"item": metric, "reason": "missing_required_input", "source": ""})
    if direction_consistent is False:
        statuses.append("phenotype_remeasurement_required")
    if peak_count is not None and peak_count > 1:
        statuses.append("linked_qtl_suspected")
    if peak_stable is False:
        statuses.append("mapping_model_review_required")
        warnings.append("lod_peak_not_stable_across_resampling_or_environments")
    if interval_bp and interval_cm and interval_cm > 0:
        bp_per_cm = interval_bp / interval_cm
        if bp_per_cm >= thresholds["low_recombination_bp_per_cm"]:
            statuses.append("low_recombination_region")
    else:
        bp_per_cm = None
    if recombinants is not None and recombinants < thresholds["min_informative_recombinants"]:
        statuses.append("insufficient_informative_recombinants")
    if interval_bp and marker_count is not None:
        markers_per_mb = marker_count / (interval_bp / 1_000_000)
        if markers_per_mb < thresholds["min_markers_per_mb"]:
            statuses.append("marker_density_insufficient")
    else:
        markers_per_mb = None
    if (
        (interval_bp is not None and interval_bp > thresholds["broad_interval_bp"])
        or (candidate_count is not None and candidate_count > thresholds["broad_interval_gene_count"])
    ):
        statuses.append("fine_mapping_required")
    statuses = list(dict.fromkeys(statuses))
    if not statuses:
        statuses = ["ready_for_candidate_triage"]

    priority = [
        "phenotype_remeasurement_required", "linked_qtl_suspected",
        "mapping_model_review_required", "low_recombination_region",
        "insufficient_informative_recombinants", "marker_density_insufficient",
        "fine_mapping_required", "ready_for_candidate_triage",
    ]
    primary = next(status for status in priority if status in statuses)
    allowed = statuses == ["ready_for_candidate_triage"]
    action = {
        "phenotype_remeasurement_required": "remeasure phenotype and verify effect direction before remapping",
        "linked_qtl_suspected": "fit linked/multiple-QTL models and collect discriminating recombinants",
        "mapping_model_review_required": "complete missing mapping inputs and review model/peak stability",
        "low_recombination_region": "increase informative population size or use association/haplotype evidence",
        "insufficient_informative_recombinants": "genotype additional recombinants before candidate exclusion",
        "marker_density_insufficient": "add markers across the support interval",
        "fine_mapping_required": "fine-map before treating a short candidate list as exhaustive",
        "ready_for_candidate_triage": "freeze the deterministic Stage-A candidate shortlist",
    }[primary]
    adequacy_rows = [{
        "primary_status": primary,
        "status": status,
        "candidate_triage_allowed": allowed,
        "recommended_next_action": action if status == primary else "",
    } for status in statuses]
    adequacy_path = out_dir / "mapping_adequacy.tsv"
    pd.DataFrame(adequacy_rows).to_csv(adequacy_path, sep="\t", index=False)
    review_path = out_dir / "mapping_review_required.tsv"
    pd.DataFrame(review_rows, columns=["item", "reason", "source"]).to_csv(review_path, sep="\t", index=False)
    recombinant_path = out_dir / "informative_recombinants.tsv"
    pd.DataFrame(columns=["individual_id", "left_marker", "right_marker", "breakpoint_interval", "phenotype"] ).to_csv(
        recombinant_path, sep="\t", index=False
    )
    gate = {
        "schema_version": "1.0",
        "mapping_gate_version": MAPPING_GATE_VERSION,
        "candidate_triage_allowed": allowed,
        "primary_status": primary,
        "statuses": statuses,
        "warnings": warnings,
        "missing_required_inputs": missing,
        "recommended_next_action": action,
        "derived_metrics": {"markers_per_mb": markers_per_mb, "bp_per_cm": bp_per_cm},
        "thresholds": thresholds,
        "threshold_profile": threshold_fingerprint or {"profile": "default", "version": MAPPING_GATE_VERSION},
        "evidence_used": [
            {"metric": metric, "value": value, "source": sources.get(metric, "")}
            for metric, value in values.items()
        ],
        "input": {"path": str(metrics_path), "sha256": hashlib.sha256(metrics_path.read_bytes()).hexdigest()},
        "ranking_score_modified": False,
        "causal_probability": False,
        "interpretation": "A failed mapping gate does not mean that no causal gene exists.",
    }
    gate_path = out_dir / "mapping_gate.json"
    gate_path.write_text(json.dumps(gate, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return {
        "adequacy": str(adequacy_path), "gate": str(gate_path),
        "informative_recombinants": str(recombinant_path), "review": str(review_path),
        "candidate_triage_allowed": allowed, "primary_status": primary,
    }
