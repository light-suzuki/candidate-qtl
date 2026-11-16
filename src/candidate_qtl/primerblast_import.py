"""Import PrimerBlast OSS JSON results into a candidate/marker ledger.

`primer-plan` writes PrimerBlast OSS commands; after they run, this importer
normalizes the JSON reports (design reports and CAPS/InDel/assay marker reports)
into one orderable primer table. It separates *predicted* specificity/risk from
*experimentally confirmed* results, distinguishes CAPS / dCAPS / InDel / plain
PCR marker types, keeps per-reference off-target counts, and preserves the
parent/F1 validation gate that must pass before a marker is expanded across F2.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from .core import runtime_path, sha256


BASE_COLUMNS = [
    "marker_id",
    "gene_id",
    "source_file",
    "mode",
    "marker_type",
    "pair_name",
    "forward",
    "reverse",
    "product_size",
    "tm_f",
    "tm_r",
    "gc_f",
    "gc_r",
    "specificity_rank",
    "pair_product_count",
    "pair_specificity",
    "coverage_status",
    "n_off_target_total",
    "n_comigrating_total",
    "gel_distinguishable",
    "predicted_risk",
    "risk_score",
    "risk_reasons",
    "caps_enzyme",
    "caps_min_gel_gap",
    "caps_distinguishable",
    "validation_status",
    "wet_validation_status",
    "parent_a_result",
    "parent_b_result",
    "f1_result",
    "experimental_status",
    "selection_reason",
    "notes",
]

EXPORT_FORMATS = {"kasp", "asp", "hrm"}


def _num(value: object) -> object:
    return "" if value is None else value


def _report_mode(report: dict) -> str:
    if "target" in report:
        return str(report.get("target", {}).get("source") or "assay")
    return str(report.get("mode") or "design")


def _marker_type(report: dict, pair: dict) -> str:
    if "target" not in report:
        return "sequencing_or_pcr"
    source = str(report.get("target", {}).get("source") or "").lower()
    caps = pair.get("caps") or {}
    if source == "snp":
        if caps.get("best_distinguishable"):
            return "caps"
        if pair.get("dcaps_candidates") or caps.get("dcaps_candidates"):
            return "dcaps_candidate"
        return "dcaps_needed"
    if source == "indel":
        return "indel_size"
    return "assay_pcr"


def _offtarget_by_db(report: dict, pair: dict) -> dict[str, int]:
    per_db = (pair.get("specificity") or {}).get("per_db")
    if per_db:
        return {str(entry.get("db")): int(entry.get("n_off_target") or 0) for entry in per_db}
    databases = report.get("databases") or []
    total = int(pair.get("n_off_target") or 0)
    if len(databases) == 1:
        return {str(databases[0]): total}
    # Total is not resolvable per database; record it under a synthetic key.
    return {"_total": total}


def _predicted_risk(pair: dict, off_total: int, gel_ok: bool) -> str:
    if pair.get("risk"):
        return str(pair["risk"])
    if off_total == 0 and gel_ok:
        return "low"
    if not gel_ok:
        return "high"
    return "medium"


def _selection_reason(marker_type: str, risk: str, gel_ok: bool, caps_ok: bool) -> str:
    if marker_type == "dcaps_needed":
        return "no_distinguishing_caps_design_dcaps_or_alternative_marker"
    if marker_type == "caps" and not caps_ok:
        return "caps_enzyme_not_gel_distinguishable_review"
    if risk == "high" or not gel_ok:
        return "high_offtarget_or_comigration_review_or_reject"
    if risk == "low":
        return "adopt_candidate_pending_parent_validation"
    return "review_before_ordering"


def _normalize_report(report: dict, source_file: str, candidate_map: dict[str, str]) -> tuple[list[dict[str, object]], set[str]]:
    mode = _report_mode(report)
    target = report.get("target") or {}
    marker_base = str(target.get("name") or report.get("template_id") or Path(source_file).stem)
    gene_id = candidate_map.get(marker_base, "")

    pairs = report.get("pairs") or []
    scored: list[tuple[tuple, dict, dict[str, int]]] = []
    for pair in pairs:
        off_map = _offtarget_by_db(report, pair)
        off_total = sum(off_map.values())
        gel_ok = bool(pair.get("gel_distinguishable", True))
        # Rank key: fewer off-targets, gel-distinguishable, lower risk, lower penalty.
        risk_order = {"low": 0, "medium": 1, "high": 2}.get(str(pair.get("risk", "")), 1)
        penalty = float(pair.get("penalty", pair.get("risk_score", 0)) or 0)
        scored.append(((off_total, 0 if gel_ok else 1, risk_order, penalty), pair, off_map))
    scored.sort(key=lambda item: item[0])

    rows: list[dict[str, object]] = []
    db_keys: set[str] = set()
    for rank, (_, pair, off_map) in enumerate(scored, 1):
        db_keys.update(off_map)
        off_total = sum(off_map.values())
        gel_ok = bool(pair.get("gel_distinguishable", True))
        marker_type = _marker_type(report, pair)
        caps = pair.get("caps") or {}
        caps_ok = bool(caps.get("best_distinguishable"))
        risk = _predicted_risk(pair, off_total, gel_ok)
        pair_name = str(pair.get("name") or f"pair_{pair.get('index', rank)}")
        row = {
            "marker_id": f"{marker_base}:{pair_name}",
            "gene_id": gene_id,
            "source_file": source_file,
            "mode": mode,
            "marker_type": marker_type,
            "pair_name": pair_name,
            "forward": pair.get("forward", ""),
            "reverse": pair.get("reverse", ""),
            "product_size": _num(pair.get("product_size")),
            "tm_f": _num(pair.get("tm_f")),
            "tm_r": _num(pair.get("tm_r")),
            "gc_f": _num(pair.get("gc_f")),
            "gc_r": _num(pair.get("gc_r")),
            "specificity_rank": rank,
            "pair_product_count": int((pair.get("specificity") or {}).get("pair_product_count", pair.get("n_pair_products", 1 + off_total)) or 0),
            "pair_specificity": "in_silico_unique" if int((pair.get("specificity") or {}).get("pair_product_count", pair.get("n_pair_products", 1 + off_total)) or 0) == 1 else "multi_product",
            "coverage_status": str(pair.get("coverage_status") or ("target_covered" if pair.get("product_size") else "unknown")),
            "n_off_target_total": off_total,
            "n_comigrating_total": _num(pair.get("n_comigrating", pair.get("n_fr_offtarget"))),
            "gel_distinguishable": "yes" if gel_ok else "no",
            "predicted_risk": risk,
            "risk_score": _num(pair.get("risk_score")),
            "risk_reasons": ";".join(pair.get("risk_reasons") or []),
            "caps_enzyme": pair.get("caps_enzyme") or caps.get("best_enzyme") or "",
            "caps_min_gel_gap": _num(caps.get("min_gel_gap")),
            "caps_distinguishable": "yes" if caps_ok else ("no" if caps else ""),
            # Predicted specificity is never confused with a wet result.
            "validation_status": "unvalidated",
            "wet_validation_status": "not_tested",
            "parent_a_result": "",
            "parent_b_result": "",
            "f1_result": "",
            "experimental_status": "pending_parent_f1_check",
            "selection_reason": _selection_reason(marker_type, risk, gel_ok, caps_ok),
            "notes": "Predicted from assembly/reference; confirm in parents and F1 before F2 expansion.",
            **{f"offtarget_{db}": count for db, count in off_map.items()},
        }
        rows.append(row)
    return rows, {f"offtarget_{db}" for db in db_keys}


def _load_candidate_map(candidates_arg: str | None) -> dict[str, str]:
    if not candidates_arg:
        return {}
    frame = pd.read_csv(runtime_path(candidates_arg), sep="\t", dtype=str).fillna("")
    mapping: dict[str, str] = {}
    for _, row in frame.iterrows():
        gene_id = str(row.get("gene_id", ""))
        if not gene_id:
            continue
        mapping[gene_id] = gene_id
        for key in ("representative_transcript", "marker_id", "gene_name"):
            value = str(row.get(key, ""))
            if value:
                mapping.setdefault(value, gene_id)
    return mapping


def _write_export(out_dir: Path, fmt: str, rows: pd.DataFrame) -> str:
    """Emit an assay-specific ordering scaffold; allele-specific fields stay blank."""
    if fmt == "kasp":
        columns = ["marker_id", "gene_id", "allele_specific_primer_a", "allele_specific_primer_b", "common_primer", "caps_enzyme", "notes"]
        data = [
            {
                "marker_id": row["marker_id"],
                "gene_id": row["gene_id"],
                "allele_specific_primer_a": "",
                "allele_specific_primer_b": "",
                "common_primer": row["reverse"],
                "caps_enzyme": row["caps_enzyme"],
                "notes": "Enter allele-specific 3' bases from confirmed parental SNP.",
            }
            for _, row in rows.iterrows()
        ]
    elif fmt == "asp":
        columns = ["marker_id", "gene_id", "forward", "reverse", "allele", "notes"]
        data = [
            {"marker_id": row["marker_id"], "gene_id": row["gene_id"], "forward": row["forward"], "reverse": row["reverse"], "allele": "", "notes": "Allele-specific PCR; set the discriminating 3' base per parent."}
            for _, row in rows.iterrows()
        ]
    else:  # hrm
        columns = ["marker_id", "gene_id", "forward", "reverse", "product_size", "notes"]
        data = [
            {"marker_id": row["marker_id"], "gene_id": row["gene_id"], "forward": row["forward"], "reverse": row["reverse"], "product_size": row["product_size"], "notes": "Keep amplicon short for HRM; confirm Tm shift between parents."}
            for _, row in rows.iterrows()
        ]
    path = out_dir / f"{fmt}_orders.tsv"
    pd.DataFrame(data, columns=columns).to_csv(path, sep="\t", index=False)
    return str(path)


def primerblast_import(
    results_args: list[str],
    out_dir_arg: str,
    candidates_arg: str | None = None,
    export_formats: list[str] | None = None,
) -> dict[str, object]:
    if not results_args:
        raise ValueError("at least one --results JSON file is required")
    exports = [fmt for fmt in (export_formats or []) if fmt]
    bad = sorted(set(exports).difference(EXPORT_FORMATS))
    if bad:
        raise ValueError(f"unknown export format(s): {', '.join(bad)}; choose from {sorted(EXPORT_FORMATS)}")

    candidate_map = _load_candidate_map(candidates_arg)
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict[str, object]] = []
    extra_columns: set[str] = set()
    fingerprints = []
    for arg in results_args:
        path = runtime_path(arg)
        report = json.loads(path.read_text(encoding="utf-8"))
        rows, db_columns = _normalize_report(report, str(path), candidate_map)
        all_rows.extend(rows)
        extra_columns.update(db_columns)
        fingerprints.append({"path": str(path), "sha256": sha256(path)})

    columns = BASE_COLUMNS + sorted(extra_columns)
    frame = pd.DataFrame(all_rows).reindex(columns=columns).fillna("")
    table_path = out_dir / "primer_candidates.tsv"
    frame.to_csv(table_path, sep="\t", index=False)

    export_paths = {fmt: _write_export(out_dir, fmt, frame) for fmt in exports}

    manifest = {
        "primer_candidates": str(table_path),
        "pair_count": int(len(frame)),
        "marker_type_counts": {str(k): int(v) for k, v in frame["marker_type"].value_counts().items()} if not frame.empty else {},
        "reference_databases": sorted(column[len("offtarget_"):] for column in extra_columns),
        "exports": export_paths,
        "input_fingerprints": fingerprints,
        "validation_gate": (
            "Every row is unvalidated with predicted risk only; parent_a/parent_b/F1 results and "
            "experimental_status must be filled before a marker is expanded across F2."
        ),
    }
    manifest_path = out_dir / "primerblast_import_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest | {"manifest": str(manifest_path)}
