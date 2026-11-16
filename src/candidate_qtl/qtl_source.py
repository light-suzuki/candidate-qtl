from __future__ import annotations

import csv
import json

import yaml

from .core import runtime_path


def qtl_source_init(manifest_arg: str, out_dir_arg: str) -> dict[str, object]:
    """Validate coordinate provenance and choose the strongest safe entry route."""
    manifest_path = runtime_path(manifest_arg)
    data = yaml.safe_load(manifest_path.read_text(encoding="utf-8")) or {}
    if not isinstance(data, dict):
        raise ValueError("QTL source manifest must be a YAML/JSON object")

    source = data.get("source", data)
    if not isinstance(source, dict):
        raise ValueError("source must be an object")

    physical = str(source.get("coordinate_basis", "")).lower() == "physical"
    interval_ready = all(source.get(key) not in (None, "") for key in ("source_reference_name", "source_contig", "source_start", "source_end"))
    marker_sequences = bool(source.get("marker_sequences_available"))
    marker_names = bool(source.get("left_marker") or source.get("right_marker") or source.get("flanking_markers"))
    search_status = str(source.get("literature_search_status", "unknown")).lower()

    if physical and interval_ready:
        route = "qtlift_interval"
        allowed = True
        next_action = "run source-to-target interval liftover"
    elif marker_sequences:
        route = "marker_liftover"
        allowed = True
        next_action = "run pair-aware marker liftover"
    elif marker_names:
        route = "marker_recovery"
        allowed = False
        next_action = "recover flanking marker or primer sequences"
    else:
        route = "unresolved_position"
        allowed = False
        next_action = "complete source interval and marker literature search"

    warnings: list[str] = []
    if search_status != "complete":
        warnings.append("literature_search_status is not complete")
    if str(source.get("coordinate_basis", "")).lower() in {"genetic", "inferred", "inferred_linear_cm_to_bp"}:
        warnings.append("genetic or inferred coordinates cannot contribute formal peak-proximity evidence")

    result = {
        "case_id": str(data.get("case_id", source.get("source_id", manifest_path.stem))),
        "source_id": str(source.get("source_id", "")),
        "route": route,
        "positional_ranking_allowed": allowed,
        "candidate_scope": "positional" if allowed else "rescue",
        "position_status": "provisional_interval" if allowed else "unresolved",
        "source_search_status": search_status,
        "required_next_action": next_action,
        "warnings": warnings,
        "input_manifest": str(manifest_path),
    }
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "qtl_source_preflight.json"
    json_path.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    review_path = out_dir / "review_required.tsv"
    with review_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["case_id", "reason", "required_next_action", "allowed_choices"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        if not allowed or warnings:
            writer.writerow({
                "case_id": result["case_id"],
                "reason": "; ".join(warnings) or route,
                "required_next_action": next_action,
                "allowed_choices": "provide physical interval | provide marker sequences | keep rescue-only",
            })
    result["output"] = str(json_path)
    result["review_required"] = str(review_path)
    return result
