"""Import QTL peak uncertainty, support intervals, and recombinant cores.

This module is an *importer*, not a QTL detector. It normalizes the position
information a mapping study already produced -- a LOD profile, an explicit
support interval, or a recombinant core -- into a single ``PeakContext`` and
emits canonical evidence rows for the genes that fall inside those features.

The position evidence stays deliberately conservative:

* A gene overlapping a recombinant core is the strongest positional signal.
* A gene overlapping the exact peak point is strong.
* A gene inside a statistical (LOD/Bayesian) support interval is weak support.
* A gene that only sits inside the broad QTL interval receives no positive
  positional evidence -- broad-interval membership is context, not evidence.

Every emitted row keeps the mapping method and the input file as its source id
so the position contribution remains auditable in ``score_contributions.tsv``.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from .core import nonempty, parse_interval, runtime_path, sha256
from .evidence_io import EVIDENCE_COLUMNS, evidence_row as _row


# Precedence order: only the single strongest positional code is emitted per
# gene so the position family is never double counted in the ranked score.
POSITION_PRECEDENCE = ("recombinant_core", "peak_overlap", "peak_proximal", "support_interval")


@dataclass
class PeakContext:
    """Normalized positional context for one QTL."""

    chrom: str = ""
    peak_position: int | None = None
    peak_method: str = "interval_only"
    support_interval_start: int | None = None
    support_interval_end: int | None = None
    support_interval_type: str = ""
    peak_confidence: str = "low"
    marker_density_near_peak: str = ""
    recombinant_core_start: int | None = None
    recombinant_core_end: int | None = None
    position_evidence_reason: str = "broad_interval_only"
    peak_window: int = 0
    interval: str | None = None
    source_ids: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def as_record(self) -> dict[str, str]:
        return {
            "chrom": self.chrom,
            "peak_method": self.peak_method,
            "peak_position": _int_or_blank(self.peak_position),
            "support_interval_start": _int_or_blank(self.support_interval_start),
            "support_interval_end": _int_or_blank(self.support_interval_end),
            "support_interval_type": self.support_interval_type,
            "peak_confidence": self.peak_confidence,
            "marker_density_near_peak": str(self.marker_density_near_peak),
            "recombinant_core_start": _int_or_blank(self.recombinant_core_start),
            "recombinant_core_end": _int_or_blank(self.recombinant_core_end),
            "position_evidence_reason": self.position_evidence_reason,
            "peak_window": str(self.peak_window) if self.peak_window else "",
            "interval": self.interval or "",
            "source_ids": ";".join(self.source_ids),
            "warnings": ";".join(self.warnings),
        }


def _int_or_blank(value: int | None) -> str:
    return "" if value is None else str(int(value))


def _split_chrom_pos(text: str) -> tuple[str | None, int]:
    if ":" in text:
        chrom, pos = text.split(":", 1)
        return chrom, int(pos)
    return None, int(text)


def support_interval_from_lod(
    profile: pd.DataFrame,
    drop: float,
) -> dict[str, object]:
    """Compute the peak and a LOD-drop support interval from a LOD profile.

    ``profile`` must contain ``chrom``, ``position``, and ``lod`` columns. The
    support interval is the contiguous run of positions around the peak where
    ``lod >= peak_lod - drop``. Additional above-threshold runs on the peak
    chromosome are reported as a multiple-peak warning; above-threshold signal
    on other chromosomes is reported as a linked/ghost QTL warning.
    """

    required = {"chrom", "position", "lod"}
    if not required.issubset(profile.columns):
        raise ValueError("LOD profile must contain chrom, position, and lod columns")
    frame = profile.copy()
    frame["position"] = pd.to_numeric(frame["position"], errors="coerce")
    frame["lod"] = pd.to_numeric(frame["lod"], errors="coerce")
    frame = frame.dropna(subset=["position", "lod"])
    if frame.empty:
        raise ValueError("LOD profile has no numeric position/lod rows")

    peak_row = frame.loc[frame["lod"].idxmax()]
    peak_chrom = str(peak_row["chrom"])
    peak_pos = int(peak_row["position"])
    peak_lod = float(peak_row["lod"])
    threshold = peak_lod - drop

    chrom_frame = frame[frame["chrom"].astype(str) == peak_chrom].sort_values("position")
    positions = chrom_frame["position"].astype(int).tolist()
    lods = chrom_frame["lod"].astype(float).tolist()

    # Contiguous above-threshold runs on the peak chromosome.
    runs: list[tuple[int, int]] = []
    run_start: int | None = None
    for pos, lod in zip(positions, lods):
        if lod >= threshold:
            if run_start is None:
                run_start = pos
            run_end = pos
        else:
            if run_start is not None:
                runs.append((run_start, run_end))
                run_start = None
    if run_start is not None:
        runs.append((run_start, run_end))

    support_start, support_end = peak_pos, peak_pos
    for start, end in runs:
        if start <= peak_pos <= end:
            support_start, support_end = start, end
            break

    warnings: list[str] = []
    extra_runs = [run for run in runs if not (run[0] <= peak_pos <= run[1])]
    if extra_runs:
        warnings.append("multiple_peaks_on_chromosome")
    other_chrom = frame[frame["chrom"].astype(str) != peak_chrom]
    if not other_chrom.empty and float(other_chrom["lod"].max()) >= threshold:
        warnings.append("linked_or_ghost_qtl_on_other_chromosome")

    markers_in_interval = int(
        ((chrom_frame["position"] >= support_start) & (chrom_frame["position"] <= support_end)).sum()
    )
    span_mb = max((support_end - support_start + 1) / 1_000_000.0, 1e-9)
    density = round(markers_in_interval / span_mb, 3)

    return {
        "chrom": peak_chrom,
        "peak_position": peak_pos,
        "peak_lod": peak_lod,
        "support_interval_start": support_start,
        "support_interval_end": support_end,
        "marker_density_near_peak": density,
        "warnings": warnings,
    }


def _load_peak_metadata(path: Path) -> dict[str, str]:
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if frame.empty:
        raise ValueError("peak metadata table is empty")
    if len(frame) > 1:
        # Multiple QTL rows are allowed upstream, but this importer resolves one
        # peak at a time; require the caller to pass a single-QTL slice.
        raise ValueError("peak metadata must describe one QTL; supply a single-row table")
    return {key: str(value).strip() for key, value in frame.iloc[0].items()}


def resolve_peak_context(
    lod_profile_arg: str | None = None,
    peak_metadata_arg: str | None = None,
    peak: str | None = None,
    interval: str | None = None,
    lod_drop: float = 1.5,
    peak_window: int = 0,
) -> PeakContext:
    """Merge available position inputs into one PeakContext.

    Precedence for each field: explicit recombinant core > explicit metadata >
    LOD-derived support interval > bare peak point > interval only. Confidence
    reflects how the position was constrained, not how likely the gene is causal.
    """

    context = PeakContext(interval=interval, peak_window=peak_window)
    reasons: list[str] = []

    if interval:
        chrom, _start, _end = parse_interval(interval)
        context.chrom = chrom

    if lod_profile_arg:
        path = runtime_path(lod_profile_arg)
        profile = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        derived = support_interval_from_lod(profile, lod_drop)
        context.chrom = str(derived["chrom"])
        context.peak_position = int(derived["peak_position"])
        context.support_interval_start = int(derived["support_interval_start"])
        context.support_interval_end = int(derived["support_interval_end"])
        context.support_interval_type = f"{lod_drop:g}_lod"
        context.marker_density_near_peak = derived["marker_density_near_peak"]
        context.peak_method = "lod_profile"
        context.warnings.extend(str(w) for w in derived["warnings"])
        context.source_ids.append(f"lod_profile:{path}")
        reasons.append(f"lod_{lod_drop:g}_support_interval")

    if peak_metadata_arg:
        path = runtime_path(peak_metadata_arg)
        meta = _load_peak_metadata(path)
        if nonempty(meta.get("chrom", "")):
            context.chrom = meta["chrom"]
        if nonempty(meta.get("peak_position", "")):
            context.peak_position = int(float(meta["peak_position"]))
        if nonempty(meta.get("peak_method", "")):
            context.peak_method = meta["peak_method"]
        if nonempty(meta.get("support_interval_start", "")) and nonempty(meta.get("support_interval_end", "")):
            context.support_interval_start = int(float(meta["support_interval_start"]))
            context.support_interval_end = int(float(meta["support_interval_end"]))
            context.support_interval_type = meta.get("support_interval_type", "") or "provided"
            reasons.append(f"{context.support_interval_type}_support_interval")
        if nonempty(meta.get("recombinant_core_start", "")) and nonempty(meta.get("recombinant_core_end", "")):
            context.recombinant_core_start = int(float(meta["recombinant_core_start"]))
            context.recombinant_core_end = int(float(meta["recombinant_core_end"]))
            reasons.append("recombinant_core_defined")
        if nonempty(meta.get("marker_density_near_peak", "")):
            context.marker_density_near_peak = meta["marker_density_near_peak"]
        source = meta.get("source_id", "") or f"peak_metadata:{path}"
        context.source_ids.append(source)

    if context.peak_position is None and peak:
        chrom, pos = _split_chrom_pos(peak)
        if chrom:
            context.chrom = chrom
        context.peak_position = pos
        if context.peak_method == "interval_only":
            context.peak_method = "point_peak"
        context.source_ids.append(f"peak:{peak}")
        reasons.append("exact_peak_point")

    # Confidence policy.
    if context.recombinant_core_start is not None:
        context.peak_confidence = "high"
    elif context.support_interval_start is not None:
        context.peak_confidence = "medium"
    elif context.peak_position is not None:
        context.peak_confidence = "medium"
    else:
        context.peak_confidence = "low"

    if not reasons:
        reasons.append("broad_interval_only" if interval else "no_position_input")
    context.position_evidence_reason = ";".join(dict.fromkeys(reasons))
    if not context.source_ids and interval:
        context.source_ids.append(f"interval:{interval}")
    return context


def _overlaps(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start <= b_end and b_start <= a_end


def position_evidence_rows(
    candidates: pd.DataFrame,
    context: PeakContext,
) -> list[dict[str, str]]:
    """Emit at most one positional evidence row per candidate gene."""

    if not {"gene_id", "chrom", "start", "end"}.issubset(candidates.columns):
        return []
    rows: list[dict[str, str]] = []
    source = context.source_ids[0] if context.source_ids else "qtl_position"
    for _, candidate in candidates.iterrows():
        gene_id = str(candidate["gene_id"])
        chrom = str(candidate["chrom"])
        start = int(candidate["start"])
        end = int(candidate["end"])
        if context.chrom and chrom != context.chrom:
            continue

        code: str | None = None
        confidence = "low"
        description = ""

        if (
            context.recombinant_core_start is not None
            and _overlaps(start, end, context.recombinant_core_start, context.recombinant_core_end)
        ):
            code = "recombinant_core"
            confidence = "high"
            description = (
                f"Gene overlaps the recombinant core "
                f"{context.chrom}:{context.recombinant_core_start}-{context.recombinant_core_end}."
            )
        elif context.peak_position is not None and start <= context.peak_position <= end:
            code = "peak_overlap"
            confidence = "high" if context.peak_method in {"lod_profile", "point_peak"} else "medium"
            description = f"Gene overlaps the reported QTL peak at {context.chrom}:{context.peak_position} ({context.peak_method})."
        elif (
            context.peak_position is not None
            and context.peak_window > 0
            and min(abs(start - context.peak_position), abs(end - context.peak_position)) <= context.peak_window
        ):
            code = "peak_proximal"
            confidence = "medium"
            description = f"Gene is within {context.peak_window} bp of the QTL peak ({context.peak_method})."
        elif (
            context.support_interval_start is not None
            and _overlaps(start, end, context.support_interval_start, context.support_interval_end)
        ):
            code = "support_interval"
            confidence = "low"
            description = (
                f"Gene lies inside the {context.support_interval_type or 'reported'} support interval "
                f"{context.chrom}:{context.support_interval_start}-{context.support_interval_end}; weak positional support only."
            )

        if code is None:
            continue
        rows.append(
            _row(
                gene_id,
                code,
                "support",
                confidence,
                context.position_evidence_reason,
                "qtl_position",
                source,
                description,
            )
        )
    return rows


def peak_context_import(
    candidates_arg: str,
    out_dir_arg: str,
    lod_profile_arg: str | None = None,
    peak_metadata_arg: str | None = None,
    peak: str | None = None,
    interval: str | None = None,
    lod_drop: float = 1.5,
    peak_window: int = 0,
) -> dict[str, object]:
    """CLI entry point: resolve a PeakContext and emit auditable position evidence."""

    context = resolve_peak_context(
        lod_profile_arg=lod_profile_arg,
        peak_metadata_arg=peak_metadata_arg,
        peak=peak,
        interval=interval,
        lod_drop=lod_drop,
        peak_window=peak_window,
    )
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)

    context_path = out_dir / "peak_context.tsv"
    pd.DataFrame([context.as_record()]).to_csv(context_path, sep="\t", index=False)

    evidence_path = out_dir / "position_evidence.tsv"
    evidence_counts: dict[str, int] = {}
    if candidates_arg:
        candidates = pd.read_csv(runtime_path(candidates_arg), sep="\t", dtype=str).fillna("")
        rows = position_evidence_rows(candidates, context)
        evidence = pd.DataFrame(rows, columns=EVIDENCE_COLUMNS)
        evidence.to_csv(evidence_path, sep="\t", index=False)
        evidence_counts = {
            str(code): int(count) for code, count in evidence["evidence_code"].value_counts().items()
        }
    else:
        pd.DataFrame(columns=EVIDENCE_COLUMNS).to_csv(evidence_path, sep="\t", index=False)

    fingerprints = []
    for value in [candidates_arg, lod_profile_arg, peak_metadata_arg]:
        if not value:
            continue
        path = runtime_path(value)
        if path.exists():
            fingerprints.append({"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256(path)})

    manifest = {
        "peak_context": str(context_path),
        "position_evidence": str(evidence_path),
        "peak_method": context.peak_method,
        "support_interval_type": context.support_interval_type,
        "peak_confidence": context.peak_confidence,
        "position_evidence_reason": context.position_evidence_reason,
        "position_evidence_counts": evidence_counts,
        "warnings": context.warnings,
        "source_ids": context.source_ids,
        "input_fingerprints": fingerprints,
        "principle": (
            "Positional membership is graded (recombinant core > peak overlap > support interval); "
            "broad-interval-only genes receive no positive positional evidence and this contribution "
            "is auditable in score_contributions.tsv."
        ),
    }
    manifest_path = out_dir / "peak_context_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest | {"manifest": str(manifest_path)}
