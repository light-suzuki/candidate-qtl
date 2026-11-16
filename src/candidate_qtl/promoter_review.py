"""Promoter / 5'UTR / 3'UTR / SV-PAV review route (issue #157).

A limited, **non-promoting** review route that runs *after* the coding hypothesis
has been examined. It extracts strand-aware promoter / 5'UTR / 3'UTR regions per
transcript, intersects them with parental SNV/indel and SV/PAV/TE calls under one
common schema, and separates repeat / TE / reference-gap / low-quality calls as
review-required.

Deliberate safety boundaries:

* Observations are ``context`` or ``pending_review`` by default; they NEVER add to
  a gene's score and NEVER promote a gene to Top. ``causal_promotion`` is always
  ``none``.
* A reference-only motif match is context, not causal evidence.
* The route reads coding-exon coverage / Wet state and only starts once the
  coding hypothesis is examined; it records which coding scope is still
  unexamined instead of silently rejecting the coding hypothesis.
* Observation ids link into the locus-hypothesis ledger as typed context
  (``promoter_or_utr`` / ``te_regulatory`` / ``sv_cnv_pav``); the number of
  observations or hypotheses never changes a gene score.
"""

from __future__ import annotations

import bisect
import json

import pandas as pd

from .core import attributes, file_fingerprint, open_text, runtime_path

# Feature tables are genome-wide (hundreds of thousands to millions of parental
# variant rows) while the reviewed regions are a handful of QTL windows, so the
# tables are streamed in chunks and pre-filtered to the reviewed chrom/ranges
# instead of being materialised in full.
FEATURE_CHUNK_ROWS = 200_000

REGION_COLUMNS = [
    "region_id", "gene_id", "transcript_id", "region_class", "chrom", "start", "end",
    "strand", "tss", "tss_selection", "tss_uncertainty_bp", "window_label", "review_status", "notes",
]
OBSERVATION_COLUMNS = [
    "observation_id", "region_id", "gene_id", "transcript_id", "region_class", "feature_type",
    "feature_id", "chrom", "start", "end", "overlap_bp", "call_flags", "observation_class",
    "review_status", "links_to_hypothesis_type", "causal_promotion", "source_id", "notes",
]
HYPOTHESIS_LINK_COLUMNS = [
    "link_id", "hypothesis_id", "hypothesis_type", "observation_id", "gene_id",
    "direction", "source_table", "link_status", "notes",
]
# Region class -> locus-hypothesis type it can inform (never a gene score).
_REGION_HYPOTHESIS_TYPE = {"promoter": "promoter_or_utr", "five_prime_utr": "promoter_or_utr", "three_prime_utr": "promoter_or_utr"}
_REVIEW_FLAG = ("repeat", "te", "transposon", "gap", "unplaced", "low", "lowmapq", "low_mapq", "boundary", "masked")


def _parse_gff3(path, gene_ids: set[str]) -> dict[str, dict[str, object]]:
    """Parse gene / mRNA / exon / UTR / CDS features for the requested genes."""
    genes: dict[str, dict[str, object]] = {}
    transcripts: dict[str, dict[str, object]] = {}
    tx_to_gene: dict[str, str] = {}
    with open_text(path) as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                continue
            chrom, feature, start, end, strand = fields[0], fields[2], fields[3], fields[4], fields[6]
            attrs = attributes(fields[8])
            if feature == "gene":
                gene_id = attrs.get("ID") or attrs.get("Name")
                if gene_id in gene_ids:
                    genes[gene_id] = {"chrom": chrom, "start": int(start), "end": int(end), "strand": strand}
            elif feature in {"mRNA", "transcript"}:
                parent = attrs.get("Parent", "")
                tx_id = attrs.get("ID", "")
                if parent in gene_ids and tx_id:
                    transcripts[tx_id] = {"gene_id": parent, "chrom": chrom, "start": int(start), "end": int(end), "strand": strand,
                                          "exons": [], "five_prime_utr": [], "three_prime_utr": [], "cds": []}
                    tx_to_gene[tx_id] = parent
            elif feature in {"exon", "five_prime_UTR", "three_prime_UTR", "CDS"}:
                parent = attrs.get("Parent", "")
                key = {"exon": "exons", "five_prime_UTR": "five_prime_utr", "three_prime_UTR": "three_prime_utr", "CDS": "cds"}[feature]
                for tx_id in parent.split(","):
                    if tx_id in transcripts:
                        transcripts[tx_id][key].append((int(start), int(end)))
    for tx_id, tx in transcripts.items():
        genes.setdefault(tx["gene_id"], {"chrom": tx["chrom"], "start": tx["start"], "end": tx["end"], "strand": tx["strand"]})
        genes[tx["gene_id"]].setdefault("transcripts", {})[tx_id] = tx
    return genes


def _promoter_and_utr_regions(
    gene_id: str, gene: dict[str, object], windows: list[int], representative: str,
) -> list[dict[str, object]]:
    transcripts: dict[str, dict[str, object]] = gene.get("transcripts", {})  # type: ignore[assignment]
    if not transcripts:
        return []
    strand = str(gene["strand"])
    chrom = str(gene["chrom"])
    # TSS per transcript is strand-aware: transcript start on +, transcript end on -.
    tss_by_tx = {tx_id: (tx["start"] if strand != "-" else tx["end"]) for tx_id, tx in transcripts.items()}
    tss_values = sorted(set(tss_by_tx.values()))
    tss_uncertainty = tss_values[-1] - tss_values[0]
    rows: list[dict[str, object]] = []
    for tx_id, tx in transcripts.items():
        tss = int(tss_by_tx[tx_id])
        selection = "representative_transcript" if tx_id == representative else "isoform"
        note = f"tss_from_{'transcript_start' if strand != '-' else 'transcript_end'}; isoform_tss_spread={tss_uncertainty}bp"
        for window in windows:
            # promoter is upstream of the TSS, strand-aware. Multiple windows are
            # separate scopes, not one absolute promoter length.
            if strand != "-":
                start, end = max(1, tss - window), max(1, tss - 1)
            else:
                start, end = tss + 1, tss + window
            rows.append({
                "gene_id": gene_id, "transcript_id": tx_id, "region_class": "promoter", "chrom": chrom,
                "start": start, "end": end, "strand": strand, "tss": tss, "tss_selection": selection,
                "tss_uncertainty_bp": tss_uncertainty, "window_label": f"promoter_{window}bp", "notes": note,
            })
        for region_class, key in (("five_prime_utr", "five_prime_utr"), ("three_prime_utr", "three_prime_utr")):
            for (start, end) in tx.get(key, []):  # type: ignore[union-attr]
                rows.append({
                    "gene_id": gene_id, "transcript_id": tx_id, "region_class": region_class, "chrom": chrom,
                    "start": int(start), "end": int(end), "strand": strand, "tss": tss, "tss_selection": selection,
                    "tss_uncertainty_bp": tss_uncertainty, "window_label": region_class, "notes": "from_gff3_utr_feature",
                })
    return rows


def _overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def _flag_review(flags: str) -> bool:
    lowered = str(flags).casefold()
    return any(marker in lowered for marker in _REVIEW_FLAG)


class _FeatureIndex:
    """Per-chromosome, start-sorted interval index over a feature table (#161).

    Replaces the O(regions x features) Python double loop. A lookup binary-searches
    the window of features whose start can still reach the region, then applies the
    same exact overlap test as before. Matches are returned in original input order,
    so observation rows and ids are identical to the full-scan implementation.
    """

    __slots__ = ("_chroms",)

    def __init__(self, features: pd.DataFrame) -> None:
        grouped: dict[str, list[tuple[int, int, int, dict[str, object]]]] = {}
        for order, feature in enumerate(features.to_dict("records")):
            # Coordinates are parsed here exactly as the scanning implementation did,
            # so malformed values still raise instead of being silently dropped.
            start, end = int(feature["start"]), int(feature["end"])
            grouped.setdefault(str(feature["chrom"]), []).append((start, end, order, feature))
        # chrom -> (sorted starts, rows sorted by start, longest feature span)
        self._chroms: dict[str, tuple[list[int], list[tuple[int, int, int, dict[str, object]]], int]] = {}
        for chrom, rows in grouped.items():
            rows.sort(key=lambda row: row[0])
            # max(0, ...) keeps the window bound valid for reversed-coordinate rows,
            # which the scan silently skipped and which stay skipped here.
            span = max(0, max((row[1] - row[0] for row in rows), default=0)) + 1
            self._chroms[chrom] = ([row[0] for row in rows], rows, span)

    def overlapping(self, chrom: str, start: int, end: int) -> list[dict[str, object]]:
        bucket = self._chroms.get(str(chrom))
        if not bucket:
            return []
        starts, rows, span = bucket
        # A feature reaching `start` from the left must begin at or after
        # `start - span`, where span bounds the longest feature on this chromosome.
        low = bisect.bisect_left(starts, start - span)
        high = bisect.bisect_right(starts, end)
        window = [row for row in rows[low:high] if row[1] >= start]
        window.sort(key=lambda row: row[2])
        return [row[3] for row in window]

    def __len__(self) -> int:
        return sum(len(rows) for _, rows, _ in self._chroms.values())


def _intersect(regions: pd.DataFrame, features: pd.DataFrame, feature_type: str) -> list[dict[str, object]]:
    observations: list[dict[str, object]] = []
    if regions.empty or features.empty:
        return observations
    index = _FeatureIndex(features)
    for region in regions.to_dict("records"):
        region_start, region_end = int(region["start"]), int(region["end"])
        for feature in index.overlapping(str(region["chrom"]), region_start, region_end):
            overlap = _overlap(region_start, region_end, int(feature["start"]), int(feature["end"]))
            if overlap <= 0:
                continue
            flags = str(feature.get("flags", ""))
            observation_class = str(feature.get("observation_class", feature_type))
            review = "review_required" if _flag_review(flags) or observation_class in {"te", "repeat_context"} else "pending_review"
            if observation_class in {"te", "transposon", "transposable_element"}:
                hypothesis_type = "te_regulatory"
            elif feature_type == "sv_pav" or observation_class in {"sv", "pav", "cnv", "sv_pav"}:
                hypothesis_type = "sv_cnv_pav"
            else:
                hypothesis_type = _REGION_HYPOTHESIS_TYPE.get(str(region["region_class"]), "promoter_or_utr")
            observations.append({
                "region_id": region["region_id"], "gene_id": region["gene_id"], "transcript_id": region["transcript_id"],
                "region_class": region["region_class"], "feature_type": feature_type, "feature_id": str(feature.get("feature_id", "")),
                "chrom": region["chrom"], "start": int(feature["start"]), "end": int(feature["end"]), "overlap_bp": overlap,
                "call_flags": flags, "observation_class": observation_class, "review_status": review,
                "links_to_hypothesis_type": hypothesis_type,
                "causal_promotion": "none", "source_id": str(feature.get("source_id", "")), "notes": str(feature.get("notes", "")),
            })
    return observations


def _coverage_gate(coverage: pd.DataFrame, gene_ids: set[str]) -> dict[str, object]:
    """Read coding-exon coverage / Wet state; the route starts only once coding is examined."""
    if coverage.empty:
        return {"route_started": False, "reason": "no_coding_coverage_supplied", "examined_genes": [], "unexamined_genes": sorted(gene_ids)}
    coding = coverage[coverage.get("region_class", pd.Series("", index=coverage.index)).astype(str).str.contains("coding|cds|exon", case=False, regex=True)]
    examined = set(coding["gene_id"].astype(str)) & gene_ids if "gene_id" in coding.columns else set()
    unexamined = sorted(gene_ids - examined)
    return {
        "route_started": bool(examined),
        "reason": "coding_examined_for_at_least_one_gene" if examined else "coding_hypothesis_not_yet_examined",
        "examined_genes": sorted(examined),
        "unexamined_genes": unexamined,
        "note": "partial or one-parent coverage does NOT reject the coding hypothesis; it only records unexamined scope",
    }


def promoter_review(
    genes_arg: str,
    gff_arg: str,
    out_dir_arg: str,
    parental_variants_arg: str | None = None,
    sv_calls_arg: str | None = None,
    coverage_arg: str | None = None,
    motif_hits_arg: str | None = None,
    hypotheses_arg: str | None = None,
    promoter_windows: list[int] | None = None,
    case_id: str = "",
    force: bool = False,
) -> dict[str, object]:
    windows = sorted({int(w) for w in (promoter_windows or [2000]) if int(w) > 0}) or [2000]
    genes_frame = pd.read_csv(runtime_path(genes_arg), sep="\t", dtype=str).fillna("")
    if "gene_id" not in genes_frame.columns:
        raise ValueError("target genes table must contain gene_id")
    gene_ids = set(genes_frame["gene_id"].astype(str))
    representative = {str(row["gene_id"]): str(row.get("representative_transcript", "")) for _, row in genes_frame.iterrows()}
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)

    coverage = pd.read_csv(runtime_path(coverage_arg), sep="\t", dtype=str).fillna("") if coverage_arg else pd.DataFrame()
    gate = _coverage_gate(coverage, gene_ids)
    if not gate["route_started"] and not force:
        gate["note"] = "route not started: examine the coding hypothesis first, or pass force=True to run the review scope explicitly"
        (out_dir / "promoter_review_manifest.json").write_text(
            json.dumps({"schema_version": "1.0", "case_id": case_id, "route_gate": gate, "started": False}, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        return {"started": False, "route_gate": gate, "manifest": str(out_dir / "promoter_review_manifest.json")}

    parsed = _parse_gff3(runtime_path(gff_arg), gene_ids)
    region_rows: list[dict[str, object]] = []
    for gene_id in sorted(gene_ids):
        gene = parsed.get(gene_id)
        if gene:
            region_rows.extend(_promoter_and_utr_regions(gene_id, gene, windows, representative.get(gene_id, "")))
    regions = pd.DataFrame(region_rows)
    if not regions.empty:
        regions.insert(0, "region_id", [f"reg_{i:04d}" for i in range(len(regions))])
        regions["review_status"] = "context"
        regions = regions[REGION_COLUMNS]
    else:
        regions = pd.DataFrame(columns=REGION_COLUMNS)
    regions.to_csv(out_dir / "promoter_utr_regions.tsv", sep="\t", index=False, lineterminator="\n")
    _write_bed(regions, out_dir / "promoter_utr_regions.bed")

    bounds = _region_bounds(regions) if not regions.empty else {}
    feature_scan: dict[str, dict[str, int]] = {}

    def _load_scoped(path_arg: str, default_class: str, label: str) -> pd.DataFrame:
        frame = _load_features(path_arg, default_class=default_class, bounds=bounds)
        feature_scan[label] = {"rows_read": int(frame.attrs.get("rows_read", len(frame))),
                               "rows_retained": int(frame.attrs.get("rows_retained", len(frame)))}
        return frame

    observations: list[dict[str, object]] = []
    if parental_variants_arg:
        variants = _load_scoped(parental_variants_arg, "parental_variant", "parental_variants")
        observations.extend(_intersect(regions, variants, "parental_variant"))
    if sv_calls_arg:
        sv = _load_scoped(sv_calls_arg, "sv_pav", "sv_calls")
        observations.extend(_intersect(regions, sv, "sv_pav"))
    if motif_hits_arg:
        # Reference-only motif hits are context, never causal, never promoting.
        motifs = _load_scoped(motif_hits_arg, "motif_context", "motif_hits")
        for observation in _intersect(regions, motifs, "reference_motif"):
            observation["review_status"] = "context"
            observation["observation_class"] = "reference_only_motif"
            observation["notes"] = "reference-only motif match: context, not causal evidence"
            observations.append(observation)

    obs_frame = pd.DataFrame(observations)
    if not obs_frame.empty:
        obs_frame.insert(0, "observation_id", [f"obs_{i:04d}" for i in range(len(obs_frame))])
        obs_frame = obs_frame[OBSERVATION_COLUMNS]
    else:
        obs_frame = pd.DataFrame(columns=OBSERVATION_COLUMNS)
    obs_frame.to_csv(out_dir / "review_observations.tsv", sep="\t", index=False, lineterminator="\n")
    review_required = obs_frame[obs_frame["review_status"] == "review_required"]
    review_required.to_csv(out_dir / "review_required.tsv", sep="\t", index=False, lineterminator="\n")
    hypothesis_links = _link_hypotheses(obs_frame, hypotheses_arg)
    hypothesis_links.to_csv(out_dir / "hypothesis_observation_links.tsv", sep="\t", index=False, lineterminator="\n")

    manifest = {
        "schema_version": "1.0",
        "case_id": case_id,
        "route_type": "limited non-promoting promoter/UTR/SV review (post coding-hypothesis)",
        "started": True,
        "route_gate": gate,
        "promoter_windows_bp": windows,
        "tss_policy": "TSS is strand-aware (transcript start on +, transcript end on -); every used transcript and its TSS selection are recorded per region, isoform TSS spread is kept as tss_uncertainty_bp",
        "gene_count": len(gene_ids),
        "region_count": len(regions),
        "transcripts_used": sorted(set(regions["transcript_id"])) if not regions.empty else [],
        "observation_count": len(obs_frame),
        "review_required_count": len(review_required),
        "observation_class_counts": {str(k): int(v) for k, v in obs_frame["observation_class"].value_counts().items()} if not obs_frame.empty else {},
        "ledger_link": "observation ids link to locus-hypothesis types (promoter_or_utr / te_regulatory / sv_cnv_pav) as typed context",
        "feature_intersection": {
            "input_format": "tab-separated feature table (chrom/start/end), 1-based inclusive; emitted BED is 0-based half-open",
            "region_bounds": {chrom: {"min_start": low, "max_end": high} for chrom, (low, high) in sorted(bounds.items())},
            "row_scan": feature_scan,
            "streaming": f"read in {FEATURE_CHUNK_ROWS}-row chunks and pre-filtered to the reviewed chrom/ranges before intersection",
            "method": "per-chromosome start-sorted interval index (binary-searched window), not an all-pairs scan",
            "tabix_index_used": False,
        },
        "regions": str(out_dir / "promoter_utr_regions.tsv"),
        "bed": str(out_dir / "promoter_utr_regions.bed"),
        "observations": str(out_dir / "review_observations.tsv"),
        "review_required": str(out_dir / "review_required.tsv"),
        "hypothesis_observation_links": str(out_dir / "hypothesis_observation_links.tsv"),
        "hypothesis_link_count": int(hypothesis_links["hypothesis_id"].ne("").sum()) if len(hypothesis_links) else 0,
        "gff3": file_fingerprint(runtime_path(gff_arg)),
        "safety": [
            "observations are context or pending_review; none adds to a gene score or promotes a gene to Top (causal_promotion is always none)",
            "a reference-only motif match is context, not causal evidence",
            "partial coverage, PCR failure, or one-parent non-amplification does not reject the coding hypothesis",
            "gene rank and locus mechanism hypothesis are kept separate; this route never writes a dry_score",
        ],
    }
    manifest_path = out_dir / "promoter_review_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return {
        "started": True,
        "manifest": str(manifest_path),
        "regions": str(out_dir / "promoter_utr_regions.tsv"),
        "bed": str(out_dir / "promoter_utr_regions.bed"),
        "observations": str(out_dir / "review_observations.tsv"),
        "review_required": str(out_dir / "review_required.tsv"),
        "hypothesis_observation_links": str(out_dir / "hypothesis_observation_links.tsv"),
        "region_count": len(regions),
        "observation_count": len(obs_frame),
        "review_required_count": len(review_required),
        "route_gate": gate,
    }


def _region_bounds(regions: pd.DataFrame) -> dict[str, tuple[int, int]]:
    """Per-chromosome min start / max end of the reviewed regions.

    A feature outside every chromosome's reviewed span cannot overlap any region,
    so it can be discarded while streaming instead of being kept in memory (#161).
    """
    bounds: dict[str, tuple[int, int]] = {}
    for region in regions.to_dict("records"):
        chrom = str(region["chrom"])
        start, end = int(region["start"]), int(region["end"])
        low, high = bounds.get(chrom, (start, end))
        bounds[chrom] = (min(low, start), max(high, end))
    return bounds


def _keep_in_bounds(chunk: pd.DataFrame, bounds: dict[str, tuple[int, int]]) -> pd.DataFrame:
    chrom = chunk["chrom"].astype(str)
    # Not coerced: a malformed coordinate must still fail loudly, as it did when
    # every row was parsed with int() during intersection.
    start = pd.to_numeric(chunk["start"])
    end = pd.to_numeric(chunk["end"])
    lows = chrom.map({name: low for name, (low, _) in bounds.items()})
    highs = chrom.map({name: high for name, (_, high) in bounds.items()})
    return chunk[lows.notna() & (end >= lows) & (start <= highs)]


def _load_features(path_arg: str, default_class: str, bounds: dict[str, tuple[int, int]] | None = None) -> pd.DataFrame:
    """Stream a feature table, keeping only rows that can reach a reviewed region.

    ``bounds`` is the per-chromosome reviewed span. When it is omitted the whole
    table is loaded, which keeps the function usable for callers that have no
    region context yet.
    """
    path = runtime_path(path_arg)
    header = pd.read_csv(path, sep="\t", dtype=str, nrows=0)
    required = {"chrom", "start", "end"}
    missing = sorted(required - set(header.columns))
    if missing:
        raise ValueError(f"feature table {path_arg} is missing columns: {', '.join(missing)}")
    kept: list[pd.DataFrame] = []
    total_rows = 0
    with pd.read_csv(path, sep="\t", dtype=str, chunksize=FEATURE_CHUNK_ROWS) as reader:
        for chunk in reader:
            chunk = chunk.fillna("")
            total_rows += len(chunk)
            if bounds is not None:
                if not bounds:
                    continue
                chunk = _keep_in_bounds(chunk, bounds)
            if not chunk.empty:
                kept.append(chunk)
    frame = pd.concat(kept, ignore_index=True) if kept else header.fillna("")
    if "observation_class" not in frame.columns:
        frame["observation_class"] = default_class
    frame.attrs["rows_read"] = total_rows
    frame.attrs["rows_retained"] = len(frame)
    return frame


def _tokens(value: object) -> set[str]:
    return {token.strip() for token in str(value).replace("|", ";").split(";") if token.strip()}


def _link_hypotheses(observations: pd.DataFrame, hypotheses_arg: str | None) -> pd.DataFrame:
    """Build typed observation links consumable beside hypothesis-ledger output.

    Missing hypotheses remain review-required; the route never invents or
    promotes a hypothesis merely because an overlap was observed.
    """
    if observations.empty:
        return pd.DataFrame(columns=HYPOTHESIS_LINK_COLUMNS)
    hypotheses = pd.read_csv(runtime_path(hypotheses_arg), sep="\t", dtype=str).fillna("") if hypotheses_arg else pd.DataFrame()
    rows: list[dict[str, str]] = []
    for observation in observations.to_dict("records"):
        gene_id = str(observation["gene_id"])
        requested_type = str(observation["links_to_hypothesis_type"])
        matched: list[dict[str, object]] = []
        if not hypotheses.empty and {"hypothesis_id", "hypothesis_type"}.issubset(hypotheses.columns):
            for hypothesis in hypotheses.to_dict("records"):
                genes = _tokens(hypothesis.get("member_gene_ids", "")) | _tokens(hypothesis.get("target_gene_ids", ""))
                if str(hypothesis.get("hypothesis_type", "")) == requested_type and gene_id in genes:
                    matched.append(hypothesis)
        if not matched:
            matched = [{"hypothesis_id": "", "hypothesis_type": requested_type}]
        for hypothesis in matched:
            hypothesis_id = str(hypothesis.get("hypothesis_id", ""))
            rows.append({
                "link_id": f"prlink_{len(rows):04d}", "hypothesis_id": hypothesis_id,
                "hypothesis_type": requested_type, "observation_id": str(observation["observation_id"]),
                "gene_id": gene_id, "direction": "context", "source_table": "review_observations.tsv",
                "link_status": "resolved" if hypothesis_id else "review_required",
                "notes": "typed context only; observation count does not affect gene score" if hypothesis_id else "add or review a matching locus hypothesis before linking",
            })
    return pd.DataFrame(rows, columns=HYPOTHESIS_LINK_COLUMNS)


def _write_bed(regions: pd.DataFrame, path) -> None:
    lines: list[str] = []
    for row in regions.to_dict("records"):
        # BED is 0-based half-open; the TSV keeps 1-based inclusive coordinates.
        start = max(0, int(row["start"]) - 1)
        name = f"{row['gene_id']}|{row['transcript_id']}|{row['region_class']}|{row['window_label']}"
        lines.append(f"{row['chrom']}\t{start}\t{int(row['end'])}\t{name}\t0\t{row['strand']}")
    path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
