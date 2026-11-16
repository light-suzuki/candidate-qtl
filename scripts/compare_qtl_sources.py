#!/usr/bin/env python3
"""Compare QTL records after QTLift harmonization to a common assembly."""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

INTERVAL = re.compile(r"^([^:]+):(\d+)-(\d+)$")
PEAK = re.compile(r"^(?:([^:]+):)?(\d+)$")
MISSING = {"", "na", "none", "pending", "unknown", "required"}
REQUIRED_COLUMNS = {"qtl_source_id", "trait", "target_assembly", "target_interval", "target_peak"}
OUTPUT_COLUMNS = [
    "qtl_a", "qtl_b", "target_assembly", "chrom", "interval_a", "interval_b",
    "overlap_bp", "jaccard", "reciprocal_overlap_a", "reciprocal_overlap_b",
    "overlap_class", "peak_distance_bp", "trait_exact_match", "same_locus_claim", "note",
]


def missing(value: object) -> bool:
    return str(value).strip().lower() in MISSING


def parse_interval(value: str) -> tuple[str, int, int]:
    match = INTERVAL.fullmatch(value.strip())
    if not match:
        raise ValueError(f"invalid interval: {value}")
    chrom, start, end = match.groups()
    start_i, end_i = int(start), int(end)
    if start_i < 1 or end_i < start_i:
        raise ValueError(f"invalid interval: {value}")
    return chrom, start_i, end_i


def parse_peak(value: str, chrom: str) -> int | None:
    if missing(value):
        return None
    match = PEAK.fullmatch(value.strip())
    if not match or (match.group(1) and match.group(1) != chrom):
        raise ValueError(f"invalid target peak: {value}")
    return int(match.group(2))


def compare(left: dict[str, str], right: dict[str, str]) -> dict[str, object] | None:
    if missing(left.get("target_interval")) or missing(right.get("target_interval")):
        return None
    if missing(left.get("target_assembly")) or missing(right.get("target_assembly")):
        return None
    if left.get("target_assembly") != right.get("target_assembly"):
        return None
    chrom_a, start_a, end_a = parse_interval(left["target_interval"])
    chrom_b, start_b, end_b = parse_interval(right["target_interval"])
    if chrom_a != chrom_b:
        return None
    overlap = max(0, min(end_a, end_b) - max(start_a, start_b) + 1)
    len_a, len_b = end_a - start_a + 1, end_b - start_b + 1
    union = len_a + len_b - overlap
    reciprocal_a, reciprocal_b = overlap / len_a, overlap / len_b
    if not overlap:
        overlap_class = "none"
    elif reciprocal_a >= 0.8 and reciprocal_b >= 0.8:
        overlap_class = "high_reciprocal"
    elif reciprocal_a >= 0.8 or reciprocal_b >= 0.8:
        overlap_class = "nested"
    elif reciprocal_a >= 0.5 and reciprocal_b >= 0.5:
        overlap_class = "moderate_reciprocal"
    else:
        overlap_class = "partial"
    peak_a = parse_peak(left.get("target_peak", ""), chrom_a)
    peak_b = parse_peak(right.get("target_peak", ""), chrom_b)
    trait_a = left.get("trait", "").strip().casefold()
    trait_b = right.get("trait", "").strip().casefold()
    return {
        "qtl_a": left.get("qtl_source_id", ""), "qtl_b": right.get("qtl_source_id", ""),
        "target_assembly": left.get("target_assembly", ""), "chrom": chrom_a,
        "interval_a": left["target_interval"], "interval_b": right["target_interval"],
        "overlap_bp": overlap, "jaccard": round(overlap / union, 6),
        "reciprocal_overlap_a": round(reciprocal_a, 6),
        "reciprocal_overlap_b": round(reciprocal_b, 6), "overlap_class": overlap_class,
        "peak_distance_bp": abs(peak_a - peak_b) if peak_a is not None and peak_b is not None else "",
        "trait_exact_match": bool(trait_a and trait_a == trait_b),
        "same_locus_claim": "not_assessed",
        "note": "Overlap is not proof of the same causal locus.",
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--same-trait-only", action="store_true")
    args = parser.parse_args(argv)
    try:
        with Path(args.input).open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing_columns = sorted(REQUIRED_COLUMNS.difference(reader.fieldnames or []))
            if missing_columns:
                raise ValueError(f"input is missing required columns: {', '.join(missing_columns)}")
            rows = list(reader)
        ids = [row["qtl_source_id"].strip() for row in rows]
        if any(not value for value in ids):
            raise ValueError("qtl_source_id must not be empty")
        if len(ids) != len(set(ids)):
            raise ValueError("qtl_source_id values must be unique")
        results = []
        for index, left in enumerate(rows):
            for right in rows[index + 1:]:
                result = compare(left, right)
                if result and (not args.same_trait_only or result["trait_exact_match"]):
                    results.append(result)
        out = Path(args.out); out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
            writer.writeheader(); writer.writerows(results)
        print(json.dumps({"records": len(rows), "comparisons": len(results), "output": str(out)}, indent=2))
        return 0
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
