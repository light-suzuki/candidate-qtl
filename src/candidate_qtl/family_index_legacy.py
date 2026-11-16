"""Genome-wide family index and interval-based hypothesis falsification.

The core configuration of this workflow is *reference genome + well-studied
model/relative genomes + QTL interval*. Everything else (accession assemblies,
parental resequencing, multi-line phenotypes) is optional enrichment. In that
core configuration the remaining bottleneck is not gene ranking - given a correct
bait, homology already recovers the causal gene - but **how many competing
pathway hypotheses have to be carried forward**.

This module attacks that, and is deliberately shaped around a laptop compute
budget:

* ``family_index`` is the only expensive step. It runs **once per reference
  genome**, assigns every gene to its best-hit reference gene plus a broader
  family bucket, and is cached with an input fingerprint. Per-case work after
  that is table lookup only.
* ``pathway_interval_falsify`` is a cascade. The coarse level (broad family) is a
  table scan that kills hypotheses with no member in the interval at all; the
  fine level (specific ortholog) is computed only for the survivors.

Deliberate boundaries:

* Enrichment is a **hypothesis filter, not a gene ranker**. The hypergeometric
  statistic is degenerate for single-member families (m=1, k=1 gives p = n/N for
  every gene in the interval), so it is never written as a per-gene score. Gene
  ordering stays with orthology confidence.
* Only ``interval_members == 0`` falsifies a hypothesis. A present-but-unenriched
  family survives; enrichment is credit, not a gate. Recall comes first.
* Tandem arrays are collapsed into blocks before the conservative p-value,
  because gene families cluster physically and linkage would otherwise look like
  enrichment.
* Genes with no assignment stay in the universe denominator and are reported as
  an explicit blind spot; a family's true size can only be larger than measured.
"""

from __future__ import annotations

import json
import math
import re
import time
import tracemalloc
from collections import Counter
from pathlib import Path

import pandas as pd

from .core import file_fingerprint, parse_interval, runtime_path
from .family_labels import alias_canonical, families_match
from .identifiers import UNIPROT_ID, canonical_accession, is_locus_tag, is_usable_symbol, symbol_candidates

SCHEMA_VERSION = "1.4"
DEFAULT_TANDEM_DISTANCE = 100_000
DEFAULT_MAX_Q = 0.05
DEFAULT_FAMILY_HITS = 5

# Issue #184: the coarse "family" value is a *discovery bucket* (naming /
# broad functional search), which must stay separate from an evolutionary
# homologous family. best_subject_stem / family_stems remain as legacy aliases
# during migration; the discovery_bucket_* / evolutionary_family_* columns
# carry the separated semantics.
FAMILY_INDEX_COLUMNS = [
    "gene_id", "protein_id", "chrom", "start", "end",
    "best_subject", "best_subject_symbol", "best_subject_stem",
    "best_identity", "best_query_coverage", "best_evalue", "best_bitscore",
    "second_bitscore", "top2_bitscore_ratio",
    # Issue #265 provenance: which original bitscore-ranked hit the symbol/best
    # row actually came from, and how far it was demoted from the true top hit.
    "symbol_source_rank", "symbol_source_bitscore_ratio",
    "family_subjects", "family_symbols", "family_stems", "assignment_status",
    "discovery_bucket_id", "discovery_bucket_method",
    "evolutionary_family_id", "evolutionary_family_source", "evolutionary_family_status",
    "index_scope",
]
FALSIFICATION_COLUMNS = [
    "hypothesis_id", "level", "pathway_or_process", "gene_family", "reference_gene_name",
    "proposed_stable_id", "resolved_by", "genome_wide_members", "genome_wide_blocks",
    "interval_members", "interval_blocks", "expected_members", "enrichment_p",
    "conservative_p", "enrichment_q", "tandem_inflated", "verdict", "interval_gene_ids", "notes",
]
SURVIVING_COLUMNS = [
    "gene_id", "chrom", "start", "end", "hypothesis_ids", "pathways", "gene_families",
    "resolved_level", "best_subject", "best_subject_symbol", "best_identity",
    "top2_bitscore_ratio", "paralog_ambiguity", "tandem_block", "orthology_hint",
    "rank_contribution", "falsifier", "review_status",
]
REVIEW_COLUMNS = ["reason", "target", "required_next_action", "allowed_choices", "status"]

_BLAST_COLUMNS = ["query_id", "subject_id", "identity", "alignment_length", "query_length",
                  "subject_length", "query_start", "query_end", "subject_start", "subject_end",
                  "evalue", "bitscore"]


# --------------------------------------------------------------------------- #
# symbols and family stems
# --------------------------------------------------------------------------- #

def _best_symbol(subject_id: str, subject_title: str, override_symbol: str = "") -> tuple[str, str, str | None]:
    """Chosen family symbol for a hit plus its provenance ``(symbol, source, confidence)``.

    Priority: a curated family-map override symbol, then the first candidate
    from the shared identifier grammar (:func:`symbol_candidates` -- a usable
    ``GN=`` symbol, a curated alias of the protein name, a symbol-like
    single-token protein name), then the legacy fallback (UniProt entry-name
    prefix, else the raw subject id).  The fallback is identifier-derived and
    recorded with source ``none`` so identifier-derived discovery buckets stay
    auditable in the manifest.
    """
    if override_symbol:
        return override_symbol, "family_map_override", "high"
    candidates = symbol_candidates(subject_id, subject_title, alias_lookup=alias_canonical)
    if candidates:
        return str(candidates[0]["value"]), str(candidates[0]["source"]), str(candidates[0]["confidence"])
    match = UNIPROT_ID.match(subject_id)
    if match:
        name = match.group(2)
        return (name.split("_")[0] if "_" in name else name), "none", None
    return subject_id, "none", None


def subject_symbol(subject_id: str, subject_title: str = "") -> str:
    """Gene symbol for a hit: first candidate from the shared identifier
    grammar (usable ``GN=``, curated alias of the protein name, symbol-like
    protein name), else the legacy fallback (UniProt entry-name prefix, else
    the subject id).
    """
    return _best_symbol(subject_id, subject_title)[0]


# Issue #265: floor on how far _pick_best_record is allowed to demote the true
# top hit in order to reach a hit with a usable gene symbol. Promotion trades a
# stronger locus-tag/accession hit for a weaker symbol hit purely so the family
# name is readable; left unbounded, a symbol hit could sit arbitrarily far
# below the true best match and still get presented as "best", overstating
# homology confidence to a reviewer. Bitscore is roughly log-linear in
# alignment quality (identity x length), so requiring the promoted hit to keep
# at least half of the top hit's bitscore is a deliberately round, conservative
# number: it screens out the pathological far-downstream promotions this issue
# was filed over while still allowing the common case (a locus tag narrowly
# outscoring a curated symbol by a few bitscore points). This is a review
# heuristic to bound an editorial substitution, not a statistically fitted
# homology cutoff -- do not tune it against benchmark pass rates.
MIN_PROMOTED_BITSCORE_RATIO = 0.5


def _pick_best_record(records: list[dict[str, object]], overrides: dict[str, dict[str, str]]) -> dict[str, object]:
    """Best hit, preferring a hit whose symbol is a usable gene symbol.

    ``records`` is already sorted by bitscore (descending).  A TrEMBL locus tag
    (MTR_5g034210) can outrank a curated SwissProt symbol (NORK) in legume DBs;
    when any top hit carries a usable symbol we promote the best such hit so
    the family name stays in the symbol namespace.

    The promotion is bounded by :data:`MIN_PROMOTED_BITSCORE_RATIO`: the first
    (highest-bitscore) usable-symbol hit is only promoted if its bitscore is
    still within that fraction of the true top hit's bitscore. Because
    ``records`` is bitscore-sorted, that first usable-symbol hit is always the
    strongest candidate for promotion, so if it fails the floor no later
    (lower-scoring) usable-symbol hit could pass it either -- the true top hit
    (locus tag/accession symbol and all) is kept instead. Callers that need to
    know whether/how far a promotion happened should not re-derive it from
    this return value alone; the row builder in :func:`family_index` records
    that separately as ``symbol_source_rank`` / ``symbol_source_bitscore_ratio``.
    """
    top_bitscore = float(records[0]["bitscore"])
    for record in records:
        override = overrides.get(str(record["subject_id"]), {})
        symbol = override.get("symbol") or subject_symbol(
            str(record["subject_id"]), str(record.get("subject_title", ""))
        )
        if is_usable_symbol(symbol):
            if record is records[0]:
                return record
            ratio = float(record["bitscore"]) / top_bitscore if top_bitscore else 0.0
            return record if ratio >= MIN_PROMOTED_BITSCORE_RATIO else records[0]
    return records[0]


def symbol_stem(symbol: str) -> str:
    """Coarse family bucket: CLV1 -> CLV, PIN8 -> PIN.

    Locus tags (two or more digit runs, e.g. AT1G75820) are left alone; stemming
    them would merge unrelated genes into one bucket.
    """
    if len(re.findall(r"\d+", symbol)) >= 2:
        return symbol
    stem = re.sub(r"[-_.]?\d+$", "", symbol)
    return stem if len(stem) >= 3 else symbol


# --------------------------------------------------------------------------- #
# hypergeometric (dependency-free, log-space)
# --------------------------------------------------------------------------- #

def _log_choose(total: int, chosen: int) -> float:
    if chosen < 0 or chosen > total:
        return -math.inf
    return math.lgamma(total + 1) - math.lgamma(chosen + 1) - math.lgamma(total - chosen + 1)


def hypergeom_sf(observed: int, universe: int, successes: int, draws: int) -> float:
    """P(X >= observed) for a hypergeometric draw. Returns 1.0 when undefined."""
    if observed <= 0 or universe <= 0 or successes <= 0 or draws <= 0:
        return 1.0
    upper = min(draws, successes)
    if observed > upper:
        return 0.0
    denominator = _log_choose(universe, draws)
    if denominator == -math.inf:
        return 1.0
    total = 0.0
    for count in range(observed, upper + 1):
        term = _log_choose(successes, count) + _log_choose(universe - successes, draws - count) - denominator
        if term > -math.inf:
            total += math.exp(term)
    return min(1.0, max(0.0, total))


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    if not p_values:
        return []
    indexed = sorted(enumerate(p_values), key=lambda item: item[1])
    total = len(p_values)
    adjusted = [1.0] * total
    running = 1.0
    for rank, (position, value) in enumerate(reversed(indexed), start=1):
        candidate = value * total / (total - rank + 1)
        running = min(running, candidate)
        adjusted[position] = min(1.0, running)
    return adjusted


# --------------------------------------------------------------------------- #
# family index build
# --------------------------------------------------------------------------- #

def _read_blast_hits(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", header=None, dtype=str, comment="#").fillna("")
    if frame.shape[1] < 12:
        raise ValueError(f"{path}: expected BLAST/DIAMOND outfmt 6 with at least 12 columns")
    keep = min(frame.shape[1], 13)
    frame = frame.iloc[:, :keep]
    frame.columns = _BLAST_COLUMNS + (["subject_title"] if keep == 13 else [])
    if "subject_title" not in frame.columns:
        frame["subject_title"] = ""
    else:
        # A thirteenth outfmt field is only useful here when it is ``stitle``.
        # Without this check a numeric field such as qcovs, or a repeated
        # accession, is fed to the identifier resolver as if it were a protein
        # description and silently degrades the family label.  Keep this
        # deliberately narrow: valid one-token protein names (for example
        # FTa2) remain accepted.
        title = frame["subject_title"].astype(str).str.strip()
        numeric = title.str.fullmatch(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?")
        repeated_subject = title.eq(frame["subject_id"].astype(str).str.strip())
        # Some pinned proteomes contain legitimate FASTA records whose entire
        # defline is a locus tag (for example LotjaGi6g1v0243000.1). BLAST then
        # emits stitle == sseqid.  Keep rejecting repeated UniProt ids and bare
        # accessions (the #318 wrong-column cases), but accept a repeated,
        # curated-grammar locus tag as identifier-only title metadata. Any
        # fallback label retains provenance source ``none`` and remains unusable
        # as a biological gene symbol.
        repeated_locus_tag = repeated_subject & title.map(is_locus_tag)
        # Identifier-shaped titles are single tokens (digit-leading locus-tag
        # fragments such as 7g085020, bare accessions).  Multi-word stitle text
        # is a legitimate protein description even when it leads with a digit
        # (e.g. "26S proteasome ..."), so whitespace separates descriptive
        # stitles from identifier metadata instead of rejecting the former.
        single_token = ~title.str.contains(r"\s")
        identifier = single_token & ~title.map(is_usable_symbol)
        invalid = title.eq("") | numeric | (repeated_subject & ~repeated_locus_tag) | (identifier & ~repeated_locus_tag)
        if invalid.any():
            rows = ", ".join(str(index + 1) for index in frame.index[invalid][:5])
            raise ValueError(
                f"{path}: column 13 must be BLAST/DIAMOND stitle; invalid value(s) at data row(s) {rows}"
            )
    for column in _BLAST_COLUMNS[2:]:
        frame[column] = pd.to_numeric(frame[column], errors="raise")
    frame["query_coverage"] = ((frame["query_end"] - frame["query_start"]).abs() + 1) / frame["query_length"].clip(lower=1)
    frame["merged_query_coverage"] = _merged_query_coverage(frame)
    frame["reviewed_subject"] = frame["subject_id"].astype(str).str.casefold().str.startswith("sp|")
    return frame


def _merged_query_coverage(frame: pd.DataFrame) -> pd.Series:
    """Union of every HSP of one (query, subject) pair, over the query length.

    A real homolog is often reported as several HSPs rather than one alignment
    (issue #567): a duplicated domain arrangement, a gap in the gene model, or
    simply a segmented alignment. Judging each HSP on its own then hides a
    98%-covered, 84%-identity match behind a 50% single-HSP threshold and the
    gene leaves the index as ``no_detectable_homology``.

    Overlapping HSPs are merged, so repeating the same region never inflates
    coverage; the per-HSP ``query_coverage`` stays on the frame for provenance.
    """
    if frame.empty:
        return pd.Series(dtype=float, index=frame.index)
    work = pd.DataFrame({
        "key": frame["query_id"].astype(str) + "\x1f" + frame["subject_id"].astype(str),
        "low": frame[["query_start", "query_end"]].min(axis=1),
        "high": frame[["query_start", "query_end"]].max(axis=1),
        "query_length": frame["query_length"].clip(lower=1),
    }).sort_values(["key", "low", "high"], kind="stable")
    # Highest HSP end seen earlier in the same pair; NaN on a pair's first row
    # so a new pair always opens a new block.
    previous_high = work.groupby("key", sort=False)["high"].cummax().shift()
    previous_high[work["key"].ne(work["key"].shift())] = float("nan")
    opens_block = previous_high.isna() | (work["low"] > previous_high + 1)
    work["block"] = opens_block.cumsum()
    blocks = work.groupby("block", sort=False).agg(
        key=("key", "first"), low=("low", "min"), high=("high", "max"))
    covered = (blocks["high"] - blocks["low"] + 1).groupby(blocks["key"]).sum()
    lengths = work.groupby("key", sort=False)["query_length"].first()
    ratio = (covered / lengths).clip(upper=1.0)
    return work["key"].map(ratio).reindex(frame.index).astype(float)


def _gene_coordinates(path_arg: str | None) -> dict[str, dict[str, object]]:
    """gene_id -> chrom/start/end, from an extract-style candidates/gene table."""
    if not path_arg:
        return {}
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    if "gene_id" not in frame.columns:
        raise ValueError("gene coordinate table requires gene_id")
    result: dict[str, dict[str, object]] = {}
    for row in frame.to_dict("records"):
        result[str(row["gene_id"])] = {
            "chrom": str(row.get("chrom", "")),
            "start": int(row["start"]) if str(row.get("start", "")).strip() else 0,
            "end": int(row["end"]) if str(row.get("end", "")).strip() else 0,
            "protein_id": str(row.get("protein_id", "")),
        }
    return result


def _protein_to_gene(path_arg: str | None) -> dict[str, str]:
    if not path_arg:
        return {}
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    columns = set(frame.columns)
    protein_column = next((name for name in ("protein_id", "fasta_record_id", "transcript_id") if name in columns), "")
    if not protein_column or "gene_id" not in columns:
        raise ValueError("id map requires gene_id and one of protein_id / fasta_record_id / transcript_id")
    return {str(row[protein_column]): str(row["gene_id"]) for row in frame.to_dict("records") if row[protein_column]}


def _family_overrides(path_arg: str | None) -> dict[str, dict[str, str]]:
    """Optional curated family labels per subject id (Pfam, orthogroup, ...)."""
    if not path_arg:
        return {}
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    if "subject_id" not in frame.columns or "family_id" not in frame.columns:
        raise ValueError("family map requires subject_id and family_id")
    return {str(row["subject_id"]): {"family_id": str(row["family_id"]),
                                     "symbol": str(row.get("symbol", "")) or subject_symbol(str(row["subject_id"]))}
            for row in frame.to_dict("records")}


def family_index(
    blast_hits_arg: str,
    out_dir_arg: str,
    genes_arg: str | None = None,
    id_map_arg: str | None = None,
    family_map_arg: str | None = None,
    reference_label: str = "",
    reference_species: str = "",
    family_hits: int = DEFAULT_FAMILY_HITS,
    min_identity: float = 25.0,
    min_query_coverage: float = 0.5,
    max_evalue: float = 1e-5,
    full_hash: bool = False,
    case_id: str = "",
    index_scope: str = "reference_genome",
) -> dict[str, object]:
    """Build the cached genome-wide family index for one reference genome.

    This is the only heavy step and it is amortised across every case run on the
    same reference. It consumes a precomputed BLAST/DIAMOND outfmt-6 table, so the
    search tool stays the operator's choice.
    """
    if family_hits < 1:
        raise ValueError("family_hits must be at least 1")
    if index_scope not in {"reference_genome", "interval_candidates"}:
        raise ValueError("index_scope must be reference_genome or interval_candidates")
    started = time.perf_counter()
    tracemalloc.start()
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    hits_path = runtime_path(blast_hits_arg)
    hits = _read_blast_hits(hits_path)
    overrides = _family_overrides(family_map_arg)
    coordinates = _gene_coordinates(genes_arg)
    protein_gene = _protein_to_gene(id_map_arg)

    hits["gene_id"] = [protein_gene.get(value, value) for value in hits["query_id"].astype(str)]
    raw_hit_genes = set(hits["gene_id"].astype(str))
    # Coverage is judged on the union of a pair's HSPs. Single-HSP
    # coverage stays on the frame, so a hit that only clears the threshold once
    # its HSPs are merged is visible as such rather than silently promoted.
    single_hsp_quality = (hits["identity"] >= min_identity) & (hits["query_coverage"] >= min_query_coverage) \
        & (hits["evalue"] <= max_evalue)
    quality = (hits["identity"] >= min_identity) & (hits["merged_query_coverage"] >= min_query_coverage) \
        & (hits["evalue"] <= max_evalue)
    qualifying_hits = hits[quality].copy()
    qualifying_genes = set(qualifying_hits["gene_id"].astype(str))
    non_qualifying_hit_genes = raw_hit_genes - qualifying_genes
    merged_hsp_recovered_genes = sorted(
        qualifying_genes - set(hits[single_hsp_quality]["gene_id"].astype(str)))

    below_identity_genes = set(hits[hits["identity"] < min_identity]["gene_id"].astype(str)) & non_qualifying_hit_genes
    below_coverage_genes = set(hits[hits["merged_query_coverage"] < min_query_coverage]["gene_id"].astype(str)) & non_qualifying_hit_genes
    above_evalue_genes = set(hits[hits["evalue"] > max_evalue]["gene_id"].astype(str)) & non_qualifying_hit_genes

    coord_genes = set(coordinates.keys())
    input_genes_total = len(coord_genes) if coord_genes else len(raw_hit_genes)
    no_hit_genes = coord_genes - raw_hit_genes if coord_genes else set()
    dropped_from_index = sorted(non_qualifying_hit_genes if not coordinates else set())

    hits = qualifying_hits
    hits = hits.sort_values(["gene_id", "bitscore", "evalue", "reviewed_subject", "subject_id"],
                            ascending=[True, False, True, False, True], kind="stable")

    rows: list[dict[str, object]] = []
    symbol_source_counts = Counter(
        {"gn": 0, "curated_alias": 0, "protein_name": 0, "family_map_override": 0, "none": 0}
    )
    for gene_id, group in hits.groupby("gene_id", sort=True):
        records = group.head(family_hits).to_dict("records")
        best = _pick_best_record(records, overrides)
        subject = str(best["subject_id"])
        override = overrides.get(subject, {})
        symbol, symbol_source, _ = _best_symbol(
            subject, str(best.get("subject_title", "")), override.get("symbol", "")
        )
        symbol_source_counts[symbol_source] += 1
        # Legacy alias: the coarse family value (curated map wins over the
        # symbol stem). Its meaning is now made explicit by the separated
        # columns below instead of silently switching with input presence.
        stem = override.get("family_id") or symbol_stem(symbol)
        curated_family = override.get("family_id", "")
        best_bitscore = float(best["bitscore"])
        top_bitscore = float(records[0]["bitscore"])
        # Issue #265: `best` may be a promoted (lower-scoring) record, not
        # records[0]. Identify it by position, not by value, so ties in
        # bitscore never misattribute provenance, and derive second_bitscore
        # from everything *except* that exact record -- never from a fixed
        # index -- so it can no longer equal or exceed a demoted best via a
        # stale index into the original ordering.
        symbol_source_rank = next(idx for idx, record in enumerate(records) if record is best) + 1
        symbol_source_bitscore_ratio = round(best_bitscore / top_bitscore, 4) if top_bitscore else 1.0
        second = max((float(record["bitscore"]) for record in records if record is not best), default=0.0)
        # min/max keeps the ratio inside (0, 1] even when the excluded true
        # top hit outscores a demoted `best`, instead of assuming best is
        # always the larger of the two (only true pre-promotion).
        top2_bitscore_ratio = (
            round(min(second, best_bitscore) / max(second, best_bitscore), 4)
            if second and best_bitscore else 0.0
        )
        subjects = [str(record["subject_id"]) for record in records]
        symbols = [
            overrides.get(str(record["subject_id"]), {}).get("symbol")
            or subject_symbol(str(record["subject_id"]), str(record.get("subject_title", "")))
            for record in records
        ]
        stems = [overrides.get(value, {}).get("family_id") or symbol_stem(name) for value, name in zip(subjects, symbols)]
        coordinate = coordinates.get(str(gene_id), {})
        rows.append({
            "gene_id": str(gene_id), "protein_id": str(best["query_id"]),
            "chrom": str(coordinate.get("chrom", "")), "start": coordinate.get("start", 0), "end": coordinate.get("end", 0),
            "best_subject": subject, "best_subject_symbol": symbol, "best_subject_stem": stem,
            "best_identity": round(float(best["identity"]), 3),
            # Coverage of the whole (query, subject) pair, not of one HSP, so a
            # segmented alignment is not reported as half a match downstream.
            "best_query_coverage": round(float(best["merged_query_coverage"]), 4),
            "best_evalue": float(best["evalue"]), "best_bitscore": best_bitscore,
            "second_bitscore": second,
            "top2_bitscore_ratio": top2_bitscore_ratio,
            "symbol_source_rank": symbol_source_rank,
            "symbol_source_bitscore_ratio": symbol_source_bitscore_ratio,
            "family_subjects": ";".join(subjects), "family_symbols": ";".join(sorted(set(symbols))),
            "family_stems": ";".join(sorted(set(stems))), "assignment_status": "assigned",
            "discovery_bucket_id": stem,
            "discovery_bucket_method": "curated_keyword" if curated_family else "symbol_stem",
            "evolutionary_family_id": curated_family,
            "evolutionary_family_source": "curated_family" if curated_family else "",
            "evolutionary_family_status": "supported" if curated_family else "unresolved",
            "index_scope": index_scope,
        })

    assigned = {str(row["gene_id"]) for row in rows}
    # Genes with coordinates but no qualifying hit stay in the universe and are
    # reported as an explicit blind spot: a family can only be larger than measured.
    for gene_id, coordinate in sorted(coordinates.items()):
        if gene_id in assigned:
            continue
        rows.append({
            "gene_id": gene_id, "protein_id": str(coordinate.get("protein_id", "")),
            "chrom": str(coordinate.get("chrom", "")), "start": coordinate.get("start", 0), "end": coordinate.get("end", 0),
            "best_subject": "", "best_subject_symbol": "", "best_subject_stem": "",
            "best_identity": 0.0, "best_query_coverage": 0.0, "best_evalue": 0.0, "best_bitscore": 0.0,
            "second_bitscore": 0.0, "top2_bitscore_ratio": 0.0,
            "symbol_source_rank": 1, "symbol_source_bitscore_ratio": 1.0,
            "family_subjects": "", "family_symbols": "", "family_stems": "",
            "assignment_status": "dark_no_qualifying_hit",
            "discovery_bucket_id": "", "discovery_bucket_method": "",
            "evolutionary_family_id": "", "evolutionary_family_source": "",
            "evolutionary_family_status": "unresolved",
            "index_scope": index_scope,
        })

    index = pd.DataFrame(rows, columns=FAMILY_INDEX_COLUMNS).sort_values("gene_id", kind="stable")
    index_path = out / "family_index.tsv"
    index.to_csv(index_path, sep="\t", index=False, lineterminator="\n")
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    dark = int(index["assignment_status"].ne("assigned").sum())
    manifest = {
        "schema_version": SCHEMA_VERSION, "case_id": case_id, "index_scope": index_scope,
        "reference_label": reference_label, "reference_species": reference_species,
        "gene_count": int(len(index)), "assigned_count": int(len(index)) - dark, "dark_count": dark,
        "input_genes_count": input_genes_total,
        "dropped_from_index_count": len(dropped_from_index),
        "exclusion_breakdown": {
            "genes_with_hits": len(raw_hit_genes),
            "genes_qualifying": len(qualifying_genes),
            "genes_no_qualifying_hits": len(non_qualifying_hit_genes),
            "genes_below_min_identity": len(below_identity_genes),
            "genes_below_min_query_coverage": len(below_coverage_genes),
            "genes_above_max_evalue": len(above_evalue_genes),
            "genes_no_blast_hits": len(no_hit_genes),
            "dropped_from_index": dropped_from_index,
            # Genes that qualify only once the HSPs of a pair are merged
            #. Listed, not just counted: the threshold did not
            # move, so a reader must be able to check each recovered gene.
            "genes_recovered_by_merged_hsp_coverage": merged_hsp_recovered_genes,
        },
        "distinct_best_subjects": int(index[index.best_subject.ne("")]["best_subject"].nunique()),
        "distinct_family_stems": int(index[index.best_subject_stem.ne("")]["best_subject_stem"].nunique()),
        "family_hits_kept": family_hits,
        "quality_thresholds": {
            "min_identity": min_identity, "min_query_coverage": min_query_coverage,
            "max_evalue": max_evalue,
            "query_coverage_basis": "merged_hsps_per_query_subject_pair",
            "query_coverage_basis_note": (
                "coverage is the union of a pair's HSPs over the query length; overlapping HSPs are "
                "merged so a repeated region never inflates coverage, and best_query_coverage reports "
                "the same merged value"
            ),
        },
        "family_definition": {
            "fine": "best-hit reference gene (co-ortholog group proxy)",
            "coarse": "gene-symbol stem of any top hit, or a curated family map when supplied",
            "override_source": str(runtime_path(family_map_arg)) if family_map_arg else None,
            "blast_symbol_policy": (
                "symbols come from the shared identifier grammar (issue #297): a usable GN=<symbol> from the "
                "optional outfmt stitle, else a curated alias of the protein name, else a symbol-like single-token "
                "protein name; GN= values that are locus tags / accessions are never adopted; without any candidate "
                "the UniProt entry-name prefix or the raw sseqid is kept as an identifier-derived fallback "
                "(provenance source none); for exact bitscore/evalue ties prefer reviewed sp| records over tr| "
                "records; promotion of a lower-scoring symbol hit is bounded by MIN_PROMOTED_BITSCORE_RATIO"
            ),
            "semantics_184": {
                "discovery_bucket": "discovery_bucket_id is a naming / broad functional search bucket (symbol_stem or curated_keyword). It is a recall-first exploration index and alone never defines a homologous family",
                "evolutionary_family": "evolutionary_family_id is the homologous-family layer this step can supply: a curated family map assignment per reference subject (source=curated_family). orthogroup/HOG and domain layers are separate inputs consumed downstream (interval-atlas)",
                "legacy_alias": "best_subject_stem / family_stems remain as legacy aliases; their value may be a symbol-stem bucket or a curated family id depending on the family map. Use the discovery_bucket_* / evolutionary_family_* columns for the separated semantics",
                "family_to_copy_gate": "family-to-copy requires orthogroup/HOG, curated family, or domain+provenance support before copy discrimination; a discovery bucket alone yields copy_set_status provisional / unresolved (issue #184)",
            },
        },
        "best_subject_symbol_provenance": {
            "vocabulary": {
                "gn": "usable GN=<symbol> from the optional outfmt stitle",
                "curated_alias": "protein name resolved through the curated alias registry (family_labels.alias_canonical)",
                "protein_name": "symbol-like single-token protein name from the UniProt defline",
                "family_map_override": "explicit symbol from the curated family map",
                "none": "no candidate: identifier-derived fallback (UniProt entry-name prefix / raw subject id) or an unassigned gene",
            },
            "confidence_by_source": {
                "gn": "high", "curated_alias": "medium", "protein_name": "medium",
                "family_map_override": "high", "none": None,
            },
            "counts_by_source": dict(symbol_source_counts),
            "unassigned_genes_without_symbol": dark,
            "identifier_derived_bucket_rule": (
                "assigned rows with source none carry an identifier-derived discovery bucket (the symbol is the "
                "UniProt entry-name prefix / raw subject id); audit identifier-derived buckets via counts_by_source.none"
            ),
        },
        "inputs": {
            "blast_hits": file_fingerprint(hits_path, full_hash=full_hash),
            "genes": file_fingerprint(runtime_path(genes_arg), full_hash=full_hash) if genes_arg else None,
            "id_map": file_fingerprint(runtime_path(id_map_arg), full_hash=full_hash) if id_map_arg else None,
            "family_map": file_fingerprint(runtime_path(family_map_arg), full_hash=full_hash) if family_map_arg else None,
        },
        "compute": {"wall_seconds": round(time.perf_counter() - started, 3), "peak_python_memory_mb": round(peak / 1_048_576, 2),
                    "note": ("built once per reference genome; per-case falsification is table lookup only"
                             if index_scope == "reference_genome" else
                             "interval-scoped fallback generated from the supplied candidate proteins; not a reusable genome-wide index")},
        "blind_spot": "genes without a qualifying hit stay in the universe denominator but cannot be assigned to a family",
        "outputs": {"family_index": str(index_path)},
    }
    manifest_path = out / "family_index_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), "family_index": str(index_path), "gene_count": int(len(index)),
            "assigned_count": int(len(index)) - dark, "dark_count": dark,
            "wall_seconds": manifest["compute"]["wall_seconds"]}


# --------------------------------------------------------------------------- #
# interval falsification
# --------------------------------------------------------------------------- #

def _tandem_blocks(members: list[dict[str, object]], distance: int) -> int:
    """Collapse physically adjacent members of the same family into one block."""
    if not members:
        return 0
    blocks = 0
    by_chrom: dict[str, list[int]] = {}
    for member in members:
        by_chrom.setdefault(str(member["chrom"]), []).append(int(member["start"] or 0))
    for starts in by_chrom.values():
        starts.sort()
        blocks += 1
        for previous, current in zip(starts, starts[1:]):
            if current - previous > distance:
                blocks += 1
    return blocks


def _interval_gene_ids(index: pd.DataFrame, intervals: list[str], candidates_arg: str | None) -> tuple[set[str], list[str]]:
    ids: set[str] = set()
    resolved: list[str] = []
    for value in intervals or []:
        chrom, start, end = parse_interval(value)
        resolved.append(value)
        scope = index[(index["chrom"].astype(str) == chrom)
                      & (pd.to_numeric(index["end"], errors="coerce").fillna(0) >= start)
                      & (pd.to_numeric(index["start"], errors="coerce").fillna(0) <= end)]
        ids.update(scope["gene_id"].astype(str))
    if candidates_arg:
        frame = pd.read_csv(runtime_path(candidates_arg), sep="\t", dtype=str).fillna("")
        if "gene_id" not in frame.columns:
            raise ValueError("candidate table requires gene_id")
        ids.update(frame["gene_id"].astype(str))
        resolved.append(f"candidates:{runtime_path(candidates_arg).name}")
    if not ids:
        raise ValueError("no interval genes resolved; supply --interval matching the index coordinates, or --candidates")
    return ids, resolved


def _members(index_records: list[dict[str, object]], accession: str, gene_name: str, family_text: str, level: str) -> list[dict[str, object]]:
    accession_key = accession.strip().casefold()
    canonical_acc_key = canonical_accession(accession).strip().casefold()
    name_key = gene_name.strip().casefold()
    stem_key = symbol_stem(gene_name).strip().casefold() if gene_name else ""
    family_key = family_text.strip().casefold()
    matched: list[dict[str, object]] = []
    for record in index_records:
        if str(record["assignment_status"]) != "assigned":
            continue
        best_subject = str(record["best_subject"]).casefold()
        best_subject_acc = canonical_accession(str(record["best_subject"])).casefold()
        if level == "fine":
            if (accession_key and (best_subject == accession_key or (canonical_acc_key and best_subject_acc == canonical_acc_key))) or \
               (name_key and str(record["best_subject_symbol"]).casefold() == name_key):
                matched.append(record)
            continue
        subjects = {value.casefold() for value in str(record["family_subjects"]).split(";") if value}
        subject_accs = {canonical_accession(value).casefold() for value in str(record["family_subjects"]).split(";") if value}
        symbols = {value.casefold() for value in str(record["family_symbols"]).split(";") if value}
        stems = {value.casefold() for value in str(record["family_stems"]).split(";") if value}
        labels = stems | symbols
        if (accession_key and (accession_key in subjects or (canonical_acc_key and canonical_acc_key in subject_accs))) \
                or (name_key and name_key in symbols) \
                or (stem_key and stem_key in stems) or (family_key and family_key in labels) \
                or any(families_match(query, label) for query in (gene_name, family_text) if query for label in labels):
            matched.append(record)
    return matched


def _falsifier(record: dict[str, object], level: str) -> str:
    reasons = ["区間が狭まってこの遺伝子が外れる"]
    if float(record.get("top2_bitscore_ratio") or 0) >= 0.9:
        reasons.append("RBH確認で1:1関係が否定される（second-best と拮抗）")
    if level == "coarse":
        reasons.append("fine level のortholog判定で別サブファミリーと判明する")
    if float(record.get("best_identity") or 0) < 40:
        reasons.append("遠縁ヒットのみで機能保存が支持されない")
    return " / ".join(reasons)


def _rescue_lookup(rescue_arg: str | None) -> dict[str, dict[str, object]]:
    """Sequence-level projection results keyed by the query protein.

    ``annotation-rescue`` projects a reference protein onto the interval's
    genomic sequence with miniprot, so it sees a family member the gene
    annotation missed entirely -- a broken model, a fragment, or a gene that is
    simply not in the GFF.
    """
    if not rescue_arg:
        return {}
    frame = pd.read_csv(runtime_path(rescue_arg), sep="\t", dtype=str).fillna("")
    missing = sorted({"query_protein_id", "hit_count"} - set(frame.columns))
    if missing:
        raise ValueError(f"annotation-rescue table is missing columns: {', '.join(missing)}")
    lookup: dict[str, dict[str, object]] = {}
    for row in frame.to_dict("records"):
        key = str(row["query_protein_id"]).strip().casefold()
        if not key:
            continue
        hits = int(pd.to_numeric(row.get("hit_count", 0), errors="coerce") or 0)
        record = lookup.setdefault(key, {"searched": True, "hits": 0, "loci": set(),
                                         "model_classes": set(), "evaluated": True})
        record["hits"] = int(record["hits"]) + hits
        # A rescue provider that failed, was unavailable, ran partially or never
        # ran is not a completed sequence-level search; neither is a query
        # explicitly marked not_evaluated. "No hit" from such a table must not
        # falsify: absence stays unproven (absent_annotation_only hold). A
        # legacy table without provider_status is treated as complete.
        provider = str(row.get("provider_status", "")).casefold()
        model_class = str(row.get("model_class", ""))
        if model_class == "not_evaluated" or provider in {"failed", "unavailable", "imported_partial", "not_run"}:
            record["evaluated"] = False
        if hits:
            locus = f"{row.get('chrom', '')}:{row.get('start', '')}-{row.get('end', '')}"
            record["loci"].add(locus)  # type: ignore[union-attr]
            record["model_classes"].add(str(row.get("model_class", "")))  # type: ignore[union-attr]
    return lookup


def pathway_interval_falsify(
    family_index_arg: str,
    hypotheses_arg: str,
    out_dir_arg: str,
    intervals: list[str] | None = None,
    candidates_arg: str | None = None,
    tandem_distance: int = DEFAULT_TANDEM_DISTANCE,
    max_q: float = DEFAULT_MAX_Q,
    case_id: str = "",
    rescue_arg: str | None = None,
) -> dict[str, object]:
    """Falsify and rank pathway hypotheses against the QTL interval.

    Coarse level first (a table scan that removes hypotheses with no member in the
    interval at all), fine level only for the survivors.

    **"Not in the annotated gene list" is not "absent".** One missed family
    member -- a fragmented model, a gene the annotation never called -- would
    otherwise falsify the correct hypothesis and push the phenotype-first route
    toward ever more distant candidates, which is exactly how that route fails.
    So an interval miss is three-valued, not two:

    * ``surviving_sequence_only_annotation_gap`` -- the annotation has no member,
      but a sequence-level projection of the reference protein hits the interval.
      The hypothesis lives, and the annotation gap is flagged.
    * ``falsified_no_member_in_interval`` -- the annotation has no member *and* a
      sequence-level search was performed and found nothing. Only now may the
      hypothesis be called falsified.
    * ``absent_annotation_only`` -- the annotation has no member and no sequence
      search was run. This is a **hold**, counted as neither falsified nor
      surviving, so the route cannot move on from an unproven absence.
    """
    started = time.perf_counter()
    tracemalloc.start()
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    index = pd.read_csv(runtime_path(family_index_arg), sep="\t", dtype=str).fillna("")
    scopes = {value for value in index.get("index_scope", pd.Series(dtype=str)).astype(str) if value}
    if not scopes:
        manifest_path = runtime_path(family_index_arg).parent / "family_index_manifest.json"
        try:
            manifest_scope = str(json.loads(manifest_path.read_text(encoding="utf-8")).get("index_scope") or "")
        except (OSError, ValueError):
            manifest_scope = ""
        if manifest_scope:
            scopes = {manifest_scope}
    if scopes and scopes != {"reference_genome"}:
        raise ValueError(
            "pathway interval falsification requires a reference_genome family index; "
            f"received index_scope={','.join(sorted(scopes))}"
        )
    missing = sorted({"gene_id", "assignment_status", "best_subject", "family_stems"} - set(index.columns))
    if missing:
        raise ValueError(f"family index is missing columns: {', '.join(missing)}")
    hypotheses = pd.read_csv(runtime_path(hypotheses_arg), sep="\t", dtype=str).fillna("")
    if "hypothesis_id" not in hypotheses.columns:
        raise ValueError("hypothesis table requires hypothesis_id")

    rescue = _rescue_lookup(rescue_arg)
    interval_ids, resolved_intervals = _interval_gene_ids(index, intervals or [], candidates_arg)
    index_records = index.to_dict("records")
    universe = len(index_records)
    draws = len(interval_ids)
    interval_records = [record for record in index_records if str(record["gene_id"]) in interval_ids]
    dark_in_interval = sum(1 for record in interval_records if str(record["assignment_status"]) != "assigned")
    indexed_gene_ids = {str(record["gene_id"]) for record in index_records}

    # Universe completeness determination:
    candidate_total = draws
    missing_from_index: list[str] = []
    if candidates_arg:
        cand_frame = pd.read_csv(runtime_path(candidates_arg), sep="\t", dtype=str).fillna("")
        if "gene_id" in cand_frame.columns:
            cand_ids = set(cand_frame["gene_id"].astype(str))
            candidate_total = len(cand_ids)
            missing_from_index = sorted(cand_ids - indexed_gene_ids)
            universe_complete = len(missing_from_index) == 0 and candidate_total > 0
        else:
            universe_complete = True
    else:
        manifest_path = runtime_path(family_index_arg).parent / "family_index_manifest.json"
        try:
            fam_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            dropped_count = int(fam_manifest.get("dropped_from_index_count", 0))
            genes_input = fam_manifest.get("inputs", {}).get("genes")
            universe_complete = bool(dropped_count == 0 and genes_input is not None)
        except (OSError, ValueError):
            has_coords = any(bool(record.get("chrom")) for record in interval_records)
            universe_complete = has_coords

    review_rows: list[dict[str, str]] = []
    if dark_in_interval:
        review_rows.append({"reason": "interval_genes_without_family_assignment", "target": f"{dark_in_interval}/{draws}",
                            "required_next_action": "これらの遺伝子はどの家系にも割り当てられていないため反証も支持もできない。annotation rescue かより近縁のリファレンスプロテオームを足す",
                            "allowed_choices": "annotation rescue | プロテオーム追加 | 盲点として明示のまま進む"})

    rows: list[dict[str, object]] = []
    surviving: dict[str, dict[str, object]] = {}
    for hypothesis in hypotheses.to_dict("records"):
        hypothesis_id = str(hypothesis["hypothesis_id"])
        accession = str(hypothesis.get("proposed_stable_id", ""))
        gene_name = str(hypothesis.get("reference_gene_name", ""))
        family_text = str(hypothesis.get("gene_family", ""))
        pathway = str(hypothesis.get("pathway_or_process", ""))
        coarse_alive = True
        coarse_verdict = ""
        for level in ("coarse", "fine"):
            if not coarse_alive:
                # Cascade: the fine level is never computed for a hypothesis the
                # coarse table scan already falsified or put on hold. The cascade
                # reason must carry the coarse verdict: holds and falsifications
                # look identical here otherwise.
                cascade_hold = coarse_verdict in ("absent_annotation_only",
                                                  "not_evaluated_family_unresolved_in_index")
                rows.append({
                    "hypothesis_id": hypothesis_id, "level": level, "pathway_or_process": pathway,
                    "gene_family": family_text, "reference_gene_name": gene_name, "proposed_stable_id": accession,
                    "resolved_by": "", "genome_wide_members": 0, "genome_wide_blocks": 0,
                    "interval_members": 0, "interval_blocks": 0, "expected_members": 0.0,
                    "enrichment_p": 1.0, "conservative_p": 1.0, "enrichment_q": 1.0,
                    "tandem_inflated": False,
                    "verdict": ("not_computed_coarse_absent_annotation_only"
                                if coarse_verdict == "absent_annotation_only"
                                else ("not_computed_coarse_family_unresolved_in_index"
                                      if cascade_hold else "not_computed_coarse_falsified")),
                    "interval_gene_ids": "",
                    "notes": ("粗い階層で区間内メンバー0だが配列検索未実施のため保留。精密計算は省略"
                              if coarse_verdict == "absent_annotation_only"
                              else ("粗い階層で宣言identifierがindex全体で未解決のため評価不能。精密計算は省略"
                                    if cascade_hold else "粗い階層で区間内メンバー0のため精密計算を省略")),
                })
                continue
            if level == "fine" and not (accession or gene_name):
                rows.append({
                    "hypothesis_id": hypothesis_id, "level": level, "pathway_or_process": pathway,
                    "gene_family": family_text, "reference_gene_name": gene_name, "proposed_stable_id": accession,
                    "resolved_by": "", "genome_wide_members": 0, "genome_wide_blocks": 0,
                    "interval_members": 0, "interval_blocks": 0, "expected_members": 0.0,
                    "enrichment_p": 1.0, "conservative_p": 1.0, "enrichment_q": 1.0,
                    "tandem_inflated": False, "verdict": "not_computed_no_fine_identifier",
                    "interval_gene_ids": "", "notes": "family label only; fine ortholog evaluation requires a reference gene or stable accession",
                })
                continue
            genome_members = _members(index_records, accession, gene_name, family_text, level)
            interval_members = [record for record in genome_members if str(record["gene_id"]) in interval_ids]
            successes = len(genome_members)
            observed = len(interval_members)
            interval_blocks = _tandem_blocks(interval_members, tandem_distance)
            genome_blocks = _tandem_blocks(genome_members, tandem_distance)
            expected = (successes * draws / universe) if universe else 0.0
            p_value = hypergeom_sf(observed, universe, successes, draws)
            # Conservative variant: tandem-collapsed observations against the
            # uncollapsed family size, so linkage cannot masquerade as enrichment.
            conservative = hypergeom_sf(interval_blocks, universe, successes, draws)
            resolved_by = ""
            if accession and any(str(record["best_subject"]).casefold() == accession.casefold() for record in genome_members):
                resolved_by = "accession"
            elif gene_name:
                resolved_by = "gene_symbol"
            elif family_text:
                resolved_by = "family_text"
            note = "enrichment は仮説フィルタであり遺伝子スコアではない"
            if successes == 0:
                # Zero resolution: none of the declared identifiers matched
                # anything anywhere in the index. This is most plausibly an
                # index-vocabulary gap (literature-standard names are often
                # stored under different BLAST-symbol stems), so the hypothesis
                # stays an evaluable hold and is never converted into a
                # biological absence claim.
                resolved_by = ""
                verdict = "not_evaluated_family_unresolved_in_index"
                note = ("宣言されたidentifierがindex全体のどこにも解決しないため評価不能。"
                        "index語彙の不足(文献標準名≠index登録語彙)の可能性があり、不在の反証にはしない")
                if level == "coarse":
                    coarse_alive = False
                    review_rows.append({"reason": "hypothesis_family_not_detected_genome_wide", "target": hypothesis_id,
                                        "required_next_action": f"{gene_name or family_text} の家系がリファレンス全体で検出されていない。baitの妥当性か、より近縁のリファレンスプロテオームを確認する",
                                        "allowed_choices": "baitを見直す | プロテオームを追加 | 真の不在として受け入れる",
                                        "status": "not_evaluated_family_unresolved_in_index"})
            elif observed == 0:
                # The annotation has no member. Whether that means "absent" depends
                # on whether anything looked at the sequence (see the docstring).
                projection = rescue.get(accession.casefold()) or rescue.get(gene_name.casefold())
                if projection and int(projection["hits"]):
                    verdict = "surviving_sequence_only_annotation_gap"
                    loci = ";".join(sorted(str(value) for value in projection["loci"]))  # type: ignore[union-attr]
                    note = (f"注釈済み遺伝子には無いが、参照タンパクの配列投影が区間に当たっている({loci})。"
                            "注釈の取りこぼしであり不在ではない")
                    review_rows.append({"reason": "family_present_in_sequence_but_missing_from_annotation",
                                        "target": hypothesis_id,
                                        "required_next_action": "投影された座位を注釈のギャップとして扱い、候補として拾い直す",
                                        "allowed_choices": "rescue候補として追加 | 遺伝子モデルを修正 | 手動確認"})
                elif projection and bool(projection.get("evaluated", True)):
                    if not universe_complete:
                        # Issue #410: Partial universe cannot conclude falsified_no_member_in_interval.
                        # Missing members might be located in unindexed candidate genes.
                        verdict = "absent_annotation_only"
                        note = ("部分的な universe（未収録の区間候補遺伝子が存在）のため、配列投影で不在であっても "
                                "falsified_no_member_in_interval とせず保留とする (absent_annotation_only)")
                        if level == "coarse":
                            coarse_alive = False
                        review_rows.append({"reason": "partial_universe_absence_unverified",
                                            "target": hypothesis_id,
                                            "required_next_action": "--genes を渡して全候補遺伝子を index に収録するか、完全な universe で再評価する",
                                            "allowed_choices": "全遺伝子で index を再構築 | 保留のまま進む"})
                    else:
                        verdict = "falsified_no_member_in_interval"
                        note = "注釈にも配列投影にもメンバーが無い。配列レベルで確認した上での不在"
                        if level == "coarse":
                            coarse_alive = False
                else:
                    # No sequence-level search was run, or the rescue provider did
                    # not complete: absence is unproven, so the hypothesis is held
                    # rather than falsified. Letting an annotation gap falsify a
                    # correct family is how the phenotype-first route drifts to
                    # ever more distant candidates.
                    verdict = "absent_annotation_only"
                    note = ("注釈にメンバーが無いが配列検索は未実施または完了していない。不在は未証明なので反証扱いにしない"
                            "(annotation-rescue の完了した出力を通すこと)")
                    if level == "coarse":
                        coarse_alive = False
                    review_rows.append({"reason": "interval_absence_not_verified_at_sequence_level",
                                        "target": hypothesis_id,
                                        "required_next_action": "annotation-rescue で参照タンパクを区間配列に投影し、真の不在か注釈の取りこぼしかを確定する",
                                        "allowed_choices": "annotation-rescue を実行 | 保留のまま進む | 真の不在として明示的に受け入れる"})
            else:
                verdict = "surviving"
            if level == "coarse":
                coarse_verdict = verdict
            rows.append({
                "hypothesis_id": hypothesis_id, "level": level, "pathway_or_process": pathway,
                "gene_family": family_text, "reference_gene_name": gene_name, "proposed_stable_id": accession,
                "resolved_by": resolved_by, "genome_wide_members": successes, "genome_wide_blocks": genome_blocks,
                "interval_members": observed, "interval_blocks": interval_blocks,
                "expected_members": round(expected, 4), "enrichment_p": p_value, "conservative_p": conservative,
                "enrichment_q": 1.0, "tandem_inflated": bool(observed > interval_blocks),
                "verdict": verdict,
                "interval_gene_ids": ";".join(sorted(str(record["gene_id"]) for record in interval_members)),
                "notes": note,
            })
            for record in interval_members:
                gene_id = str(record["gene_id"])
                bucket = surviving.setdefault(gene_id, {
                    "gene_id": gene_id, "chrom": str(record["chrom"]), "start": record["start"], "end": record["end"],
                    "hypothesis_ids": set(), "pathways": set(), "gene_families": set(), "resolved_level": level,
                    "best_subject": str(record["best_subject"]), "best_subject_symbol": str(record["best_subject_symbol"]),
                    "best_identity": record["best_identity"], "top2_bitscore_ratio": record["top2_bitscore_ratio"],
                })
                bucket["hypothesis_ids"].add(hypothesis_id)
                if pathway:
                    bucket["pathways"].add(pathway)
                if family_text:
                    bucket["gene_families"].add(family_text)
                if level == "fine":
                    bucket["resolved_level"] = "fine"

    # Multiple testing across every hypothesis actually evaluated at this level.
    for level in ("coarse", "fine"):
        scope = [row for row in rows if row["level"] == level and row["verdict"] == "surviving"]
        adjusted = benjamini_hochberg([float(row["conservative_p"]) for row in scope])
        for row, value in zip(scope, adjusted):
            row["enrichment_q"] = value
            if value <= max_q:
                row["verdict"] = "surviving_enriched"
            else:
                row["verdict"] = "surviving_uninformative"

    falsification = pd.DataFrame(rows, columns=FALSIFICATION_COLUMNS) if rows else pd.DataFrame(columns=FALSIFICATION_COLUMNS)

    tandem_flags: dict[str, bool] = {}
    for gene_id in surviving:
        record = next(item for item in interval_records if str(item["gene_id"]) == gene_id)
        same_family = [item for item in interval_records
                       if str(item["best_subject_stem"]) and str(item["best_subject_stem"]) == str(record["best_subject_stem"])]
        tandem_flags[gene_id] = _tandem_blocks(same_family, tandem_distance) < len(same_family)

    surviving_rows: list[dict[str, object]] = []
    for gene_id in sorted(surviving):
        bucket = surviving[gene_id]
        ratio = float(bucket["top2_bitscore_ratio"] or 0)
        ambiguity = "paralog_ambiguous" if ratio >= 0.9 else ("one_to_one_candidate" if ratio <= 0.6 else "intermediate")
        surviving_rows.append({
            "gene_id": gene_id, "chrom": bucket["chrom"], "start": bucket["start"], "end": bucket["end"],
            "hypothesis_ids": ";".join(sorted(bucket["hypothesis_ids"])),
            "pathways": ";".join(sorted(bucket["pathways"])), "gene_families": ";".join(sorted(bucket["gene_families"])),
            "resolved_level": bucket["resolved_level"], "best_subject": bucket["best_subject"],
            "best_subject_symbol": bucket["best_subject_symbol"], "best_identity": bucket["best_identity"],
            "top2_bitscore_ratio": bucket["top2_bitscore_ratio"], "paralog_ambiguity": ambiguity,
            "tandem_block": tandem_flags.get(gene_id, False),
            "orthology_hint": "gene ordering belongs to orthology confidence, not to family enrichment",
            "rank_contribution": "none",
            "falsifier": _falsifier(bucket, str(bucket["resolved_level"])),
            "review_status": "human_reviewable",
        })
    survivors = pd.DataFrame(surviving_rows, columns=SURVIVING_COLUMNS) if surviving_rows else pd.DataFrame(columns=SURVIVING_COLUMNS)

    if len(survivors) > 30:
        review_rows.append({"reason": "candidate_set_above_human_review_size", "target": str(len(survivors)),
                            "required_next_action": "候補が人手で全件見られる規模(~30)を超えている。仮説を絞るか、発現・区間狭小化など別クラスの証拠を足す",
                            "allowed_choices": "仮説を絞る | 発現データを足す | 区間を狭める | このまま全件レビュー"})
    reviews = pd.DataFrame(review_rows, columns=REVIEW_COLUMNS).drop_duplicates().sort_values(["reason", "target"], kind="stable") \
        if review_rows else pd.DataFrame(columns=REVIEW_COLUMNS)

    paths = {"hypothesis_interval_falsification": out / "hypothesis_interval_falsification.tsv",
             "surviving_candidates": out / "surviving_candidates.tsv",
             "review_required": out / "review_required.tsv"}
    falsification.to_csv(paths["hypothesis_interval_falsification"], sep="\t", index=False, lineterminator="\n")
    survivors.to_csv(paths["surviving_candidates"], sep="\t", index=False, lineterminator="\n")
    reviews.to_csv(paths["review_required"], sep="\t", index=False, lineterminator="\n")
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    coarse_falsified = int(falsification[(falsification.level == "coarse")
                                         & falsification.verdict.astype(str).str.startswith("falsified")].shape[0]) if not falsification.empty else 0
    coarse_family_unresolved = int(falsification[(falsification.level == "coarse")
                                                 & falsification.verdict.astype(str).eq("not_evaluated_family_unresolved_in_index")].shape[0]) \
        if not falsification.empty else 0
    held = int(falsification[falsification.verdict.astype(str).str.startswith("absent_annotation_only")].shape[0]) \
        if not falsification.empty else 0
    rescued = int(falsification[falsification.verdict.astype(str)
                                .eq("surviving_sequence_only_annotation_gap")].shape[0]) if not falsification.empty else 0
    manifest = {
        "schema_version": SCHEMA_VERSION, "case_id": case_id,
        "intervals": resolved_intervals, "universe_genes": universe, "interval_genes": draws,
        "interval_genes_without_family_assignment": dark_in_interval,
        "hypotheses_tested": int(len(hypotheses)),
        "coarse_falsified": coarse_falsified,
        "coarse_family_unresolved_in_index": coarse_family_unresolved,
        "absence_verification": {
            "sequence_search_supplied": bool(rescue_arg),
            "queries_with_sequence_search": len(rescue),
            "queries_evaluated": int(sum(1 for item in rescue.values() if bool(item.get("evaluated", True)))),
            "queries_not_evaluated": int(sum(1 for item in rescue.values() if not bool(item.get("evaluated", True)))),
            "held_absent_annotation_only": held,
            "surviving_by_sequence_only_annotation_gap": rescued,
            "universe_complete": universe_complete,
            "rule": (
                "an interval miss is three-valued: surviving_sequence_only_annotation_gap (the annotation missed a "
                "member the sequence projection found), falsified_no_member_in_interval (annotation and sequence "
                "search both found nothing on a complete universe), absent_annotation_only (no sequence search ran, "
                "the rescue provider did not complete, or run was on a partial universe -- a hold, counted as neither "
                "falsified nor surviving)"
            ),
            "completeness_rule": (
                "a rescue row with model_class=not_evaluated or a provider_status of failed/unavailable/"
                "imported_partial/not_run never counts as a completed search: its queries stay held. "
                "A partial universe similarly keeps unobserved queries held"
            ),
            "why": (
                "one missed family member would otherwise falsify the correct hypothesis and push the "
                "phenotype-first route toward ever more distant candidates"
            ),
        },
        "universe_completeness": {
            "universe_complete": universe_complete,
            "candidate_total": candidate_total,
            "universe_genes_in_interval": draws,
            "missing_from_index_count": len(missing_from_index),
            "missing_from_index": missing_from_index,
            "rule": (
                "falsified_no_member_in_interval requires a complete universe (all candidate genes accounted for "
                "in the index); partial universes degrade unobserved families to absent_annotation_only"
            ),
        },
        "surviving_hypotheses": sorted({str(row["hypothesis_id"]) for row in rows if str(row["verdict"]).startswith("surviving")}),
        "candidate_count": int(len(survivors)),
        "candidate_reduction": round(draws / len(survivors), 2) if len(survivors) else None,
        "human_reviewable": bool(len(survivors) and len(survivors) <= 30),
        "tandem_distance": tandem_distance, "max_q": max_q,
        "statistic": {
            "test": "hypergeometric survival, Benjamini-Hochberg adjusted per level",
            "conservative_rule": "tandem-collapsed interval blocks are used as the observed count against the uncollapsed family size",
            "role": "hypothesis filter only; degenerate for single-member families (m=1,k=1 gives p=n/N for every interval gene) and never written as a per-gene score",
        },
        "cascade": {"coarse": "broad family bucket, table scan", "fine": "best-hit reference gene, computed only for coarse survivors"},
        "compute": {"wall_seconds": round(time.perf_counter() - started, 3),
                    "peak_python_memory_mb": round(peak / 1_048_576, 2),
                    "note": "no sequence search in this stage; the family index is reused"},
        "inputs": {"family_index": file_fingerprint(runtime_path(family_index_arg)),
                   "hypotheses": file_fingerprint(runtime_path(hypotheses_arg))},
        "safety": [
            "only an interval member count of zero falsifies a hypothesis; present-but-unenriched families survive",
            "family enrichment never contributes to a gene score (rank_contribution is always none)",
            "tandem arrays are collapsed before the conservative p-value so linkage cannot look like enrichment",
            "unassigned interval genes stay in the denominator and are reported as a blind spot",
        ],
        "outputs": {name: str(path) for name, path in sorted(paths.items())},
    }
    manifest_path = out / "pathway_interval_falsify_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), "candidate_count": int(len(survivors)),
            "interval_genes": draws, "coarse_falsified": coarse_falsified,
            "candidate_reduction": manifest["candidate_reduction"],
            "human_reviewable": manifest["human_reviewable"],
            "universe_complete": universe_complete,
            "wall_seconds": manifest["compute"]["wall_seconds"],
            **{name: str(path) for name, path in paths.items()}}


# --------------------------------------------------------------------------- #
# Route A: verified bait -> interval family match
# --------------------------------------------------------------------------- #

# Route A per-hypothesis status vocabulary. ``completed`` (the stage ran) and
# ``evaluable`` (it had a verified bait, a reference-genome index and a
# resolvable interval) are separate: a stage that ran with zero verified baits
# is completed but never evaluated, and never a miss.
ROUTE_A_FAMILY_STATUS_COLUMNS = [
    "hypothesis_id", "gene_family", "reference_gene_name", "proposed_stable_id",
    "reference_species", "source_ids", "verdict", "family_status", "match_level",
    "interval_gene_ids", "annotation_rescue_status", "rank_contribution",
    "stage_completed", "evaluable", "notes",
]

# Deterministic mapping from the falsification verdict (fine level preferred,
# matching phenotype-first-status) to the Route A family status. Additional
# statuses beyond the minimum contract distinguish honest holds from
# confirmed interval absence: zero resolution of the declared identifiers in
# the index (``not_evaluated_family_unresolved_in_index``) is an index-
# vocabulary gap, never a biological absence claim.
# ``not_computed_*`` rows are skipped by ``_effective_verdicts`` and never
# reach this mapping.
_VERDICT_TO_FAMILY_STATUS = {
    "surviving": "surviving_family_present_in_interval",
    "surviving_enriched": "surviving_family_present_in_interval",
    "surviving_uninformative": "surviving_family_present_in_interval",
    "surviving_sequence_only_annotation_gap": "surviving_sequence_only_annotation_gap",
    "falsified_no_member_in_interval": "falsified_no_member_after_complete_sequence_search",
    "not_evaluated_family_unresolved_in_index": "not_evaluated_family_unresolved_in_index",
    "absent_annotation_only": "hold_absent_annotation_only",
}

NOT_EVALUATED_NO_VERIFIED_BAIT = "not_evaluated_no_verified_bait"
NOT_EVALUATED_MISSING_FAMILY_INDEX = "not_evaluated_missing_reference_genome_family_index"
NOT_EVALUATED_INTERVAL_UNRESOLVED = "not_evaluated_interval_unresolved"
NOT_EVALUATED_FAMILY_UNRESOLVED_IN_INDEX = "not_evaluated_family_unresolved_in_index"


def _effective_verdicts(falsification: pd.DataFrame) -> dict[str, dict[str, dict[str, str]]]:
    """Per-hypothesis falsification verdicts, keeping coarse and fine separate.

    Mirrors phenotype-first-status: iterate coarse then fine, skip
    ``not_computed_*`` rows, and prefer the fine level. The coarse level is
    kept separately because it answers a different question: whether the
    *family bucket* has any member in the interval at all, regardless of the
    verified target's specific ortholog.
    """
    effective: dict[str, dict[str, dict[str, str]]] = {}
    if falsification.empty:
        return effective
    levels = ("coarse", "fine") if "level" in falsification.columns else ("coarse",)
    for level in levels:
        scope = falsification[falsification["level"].eq(level)]
        for row in scope.to_dict("records"):
            verdict = str(row["verdict"])
            if verdict.startswith("not_computed_"):
                continue
            effective.setdefault(str(row["hypothesis_id"]), {})[level] = {
                "verdict": verdict,
                "interval_gene_ids": str(row.get("interval_gene_ids", "")),
            }
    return effective


def _family_status(verdict: str, interval_gene_ids: str, paralog_ambiguous: bool) -> str:
    """Map an effective falsification verdict to the Route A status vocabulary."""
    status = _VERDICT_TO_FAMILY_STATUS.get(verdict, "not_evaluated")
    if interval_gene_ids and status == "surviving_family_present_in_interval":
        # A concrete interval gene was matched, unless every matched gene is a
        # paralog tie: then the family is present but the copy is unresolved.
        status = "family_present_copy_unresolved" if paralog_ambiguous else "concrete_interval_gene_matched"
    return status


def route_a_interval_family_match(
    panel_arg: str,
    hypotheses_arg: str,
    out_dir_arg: str,
    family_index_arg: str | None = None,
    intervals: list[str] | None = None,
    candidates_arg: str | None = None,
    rescue_arg: str | None = None,
    tandem_distance: int = DEFAULT_TANDEM_DISTANCE,
    max_q: float = DEFAULT_MAX_Q,
    case_id: str = "",
) -> dict[str, object]:
    """Route A stage: forward-match verified baits/families into the interval.

    This is the wiring adapter between ``bait-panel-freeze`` and
    ``phenotype-first-status`` (issue #355). It deliberately does no new family
    inference: only the verified/source-backed targets frozen by the bait panel
    are forwarded to ``pathway_interval_falsify``, and only a
    ``reference_genome``-scoped family index (fail-closed) is accepted.

    When the stage cannot evaluate -- no verified bait, no family index, or no
    resolvable interval -- it still completes and records the reason as a
    ``not_evaluated_*`` family status; a verified-bait count of zero is never a
    miss and never an interval match.
    """
    started = time.perf_counter()
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    falsification_path = out / "hypothesis_interval_falsification.tsv"
    status_path = out / "route_a_interval_family_match.tsv"

    panel = pd.read_csv(runtime_path(panel_arg), sep="\t", dtype=str).fillna("") if panel_arg else pd.DataFrame()
    hypotheses = pd.read_csv(runtime_path(hypotheses_arg), sep="\t", dtype=str).fillna("") if hypotheses_arg else pd.DataFrame()
    for required, name in ((panel, "verified bait panel"), (hypotheses, "hypothesis table")):
        if not required.empty and "hypothesis_id" not in required.columns:
            raise ValueError(f"{name} requires hypothesis_id")

    verified_ids = set(panel["hypothesis_id"].astype(str)) if not panel.empty else set()
    verified_hypotheses = hypotheses[hypotheses["hypothesis_id"].astype(str).isin(verified_ids)] \
        if not hypotheses.empty else pd.DataFrame()
    if not verified_hypotheses.empty:
        verified_hypotheses.to_csv(out / "verified_target_hypotheses.tsv", sep="\t", index=False, lineterminator="\n")

    source_ids_by_hypothesis: dict[str, str] = {}
    panel_info: dict[str, dict[str, str]] = {}
    if not panel.empty:
        for row in panel.to_dict("records"):
            hypothesis_id = str(row["hypothesis_id"])
            source_ids_by_hypothesis.setdefault(hypothesis_id, str(row.get("function_source_ids", "")))
            panel_info.setdefault(hypothesis_id, {
                "gene_family": str(row.get("gene_family", "")),
                "reference_gene_name": str(row.get("reference_gene_name", "")),
                "proposed_stable_id": str(row.get("bait_id", "") or row.get("reference_gene_name", "")),
                "reference_species": str(row.get("reference_species", "")),
            })

    # Failure-free degradation: record why the stage is not evaluated and stop.
    not_evaluated_reason = ""
    if not verified_ids or verified_hypotheses.empty:
        not_evaluated_reason = NOT_EVALUATED_NO_VERIFIED_BAIT
    elif not family_index_arg:
        not_evaluated_reason = NOT_EVALUATED_MISSING_FAMILY_INDEX
    elif not intervals and not candidates_arg:
        not_evaluated_reason = NOT_EVALUATED_INTERVAL_UNRESOLVED
    else:
        # An interval/candidate table that resolves to no genes is an
        # unresolved interval, not a falsification: degrade to
        # not_evaluated_interval_unresolved instead of failing the lane.
        # Other index errors (e.g. an interval-scoped index) still fail
        # closed inside pathway_interval_falsify.
        try:
            _interval_gene_ids(pd.read_csv(runtime_path(family_index_arg), sep="\t", dtype=str).fillna(""),
                               list(intervals or []), candidates_arg)
        except ValueError as error:
            if str(error).startswith("no interval genes resolved"):
                not_evaluated_reason = NOT_EVALUATED_INTERVAL_UNRESOLVED
            else:
                raise

    rescue_lookup = _rescue_lookup(rescue_arg)
    rows: list[dict[str, object]] = []
    falsify_result: dict[str, object] = {}
    if not_evaluated_reason:
        # A rerun into the same output directory must not expose artifacts from
        # an earlier evaluated run as if they belonged to this not-evaluated
        # result.  These files are all owned by pathway_interval_falsify().
        for stale_path in (
            falsification_path,
            out / "surviving_candidates.tsv",
            out / "review_required.tsv",
            out / "pathway_interval_falsify_manifest.json",
        ):
            stale_path.unlink(missing_ok=True)
        for hypothesis_id in sorted(verified_ids):
            info = panel_info.get(hypothesis_id, {})
            rows.append({
                "hypothesis_id": hypothesis_id, "gene_family": info.get("gene_family", ""),
                "reference_gene_name": info.get("reference_gene_name", ""),
                "proposed_stable_id": info.get("proposed_stable_id", ""),
                "reference_species": info.get("reference_species", ""),
                "source_ids": source_ids_by_hypothesis.get(hypothesis_id, ""),
                "verdict": "", "family_status": not_evaluated_reason, "match_level": "none",
                "interval_gene_ids": "", "annotation_rescue_status": _rescue_status(rescue_lookup, hypothesis_id, verified_hypotheses),
                "rank_contribution": "none", "stage_completed": True, "evaluable": False,
                "notes": f"Route A interval family match not evaluated: {not_evaluated_reason}",
            })
    else:
        verified_hypotheses_path = out / "verified_target_hypotheses.tsv"
        falsify_result = pathway_interval_falsify(
            family_index_arg=str(family_index_arg), hypotheses_arg=str(verified_hypotheses_path),
            out_dir_arg=str(out_dir_arg), intervals=list(intervals or []), candidates_arg=candidates_arg,
            tandem_distance=tandem_distance, max_q=max_q, case_id=case_id, rescue_arg=rescue_arg,
        )
        falsification = pd.read_csv(falsification_path, sep="\t", dtype=str).fillna("")
        effective = _effective_verdicts(falsification)
        paralog_by_gene: dict[str, bool] = {}
        survivors_path = out / "surviving_candidates.tsv"
        if survivors_path.exists():
            survivors = pd.read_csv(survivors_path, sep="\t", dtype=str).fillna("")
            if not survivors.empty and "paralog_ambiguity" in survivors.columns:
                for record in survivors.to_dict("records"):
                    paralog_by_gene[str(record["gene_id"])] = str(record["paralog_ambiguity"]) == "paralog_ambiguous"
        for hypothesis_id in sorted(verified_ids):
            info = panel_info.get(hypothesis_id, {})
            levels = effective.get(hypothesis_id, {})
            coarse = levels.get("coarse", {})
            fine = levels.get("fine", {})
            verdict = str(fine.get("verdict") or coarse.get("verdict") or "")
            gene_ids = str(fine.get("interval_gene_ids", "")) if fine.get("verdict") \
                else str(coarse.get("interval_gene_ids", ""))
            matched_genes = [gene for gene in gene_ids.split(";") if gene]
            paralog_ambiguous = bool(matched_genes) and all(paralog_by_gene.get(gene, False) for gene in matched_genes)
            family_status = _family_status(verdict, gene_ids, paralog_ambiguous)
            note = (f"interval family falsification completed; effective verdict: {verdict}"
                    if verdict else "interval family falsification did not evaluate this hypothesis")
            # The family bucket is present in the interval (coarse surviving)
            # but the verified target's fine-level ortholog is absent: the
            # family is there, the specific copy is unresolved. The coarse
            # members stay visible in interval_gene_ids; the fine absence is
            # recorded in the note instead of erasing the family presence.
            fine_verdict = str(fine.get("verdict") or "")
            if (fine_verdict and str(coarse.get("verdict", "")).startswith("surviving")
                    and not fine_verdict.startswith("surviving")):
                family_status = "family_present_copy_unresolved"
                gene_ids = str(coarse.get("interval_gene_ids", ""))
                note += "; coarse family present in the interval but the fine-level verified ortholog is absent (copy unresolved)"
            elif paralog_ambiguous:
                note += "; all matched interval genes are paralog ties (copy unresolved)"
            if family_status in ("concrete_interval_gene_matched", "family_present_copy_unresolved",
                                 "surviving_family_present_in_interval"):
                match_level = "full"
            elif family_status == "surviving_sequence_only_annotation_gap":
                match_level = "sequence_only"
            else:
                match_level = "none"
            rows.append({
                "hypothesis_id": hypothesis_id, "gene_family": info.get("gene_family", ""),
                "reference_gene_name": info.get("reference_gene_name", ""),
                "proposed_stable_id": info.get("proposed_stable_id", ""),
                "reference_species": info.get("reference_species", ""),
                "source_ids": source_ids_by_hypothesis.get(hypothesis_id, ""),
                "verdict": verdict, "family_status": family_status, "match_level": match_level,
                "interval_gene_ids": gene_ids,
                "annotation_rescue_status": _rescue_status(rescue_lookup, hypothesis_id, verified_hypotheses),
                "rank_contribution": "none", "stage_completed": True, "evaluable": True,
                "notes": note,
            })

    status_frame = pd.DataFrame(rows, columns=ROUTE_A_FAMILY_STATUS_COLUMNS) \
        if rows else pd.DataFrame(columns=ROUTE_A_FAMILY_STATUS_COLUMNS)
    status_frame.to_csv(status_path, sep="\t", index=False, lineterminator="\n")
    counts = status_frame["family_status"].value_counts().to_dict() if not status_frame.empty else {}
    manifest = {
        "schema_version": SCHEMA_VERSION, "case_id": case_id,
        "route": "a", "stage": "interval_family_match",
        "verified_bait_hypotheses": sorted(verified_ids),
        "verified_target_hypotheses_tsv": str(out / "verified_target_hypotheses.tsv") if not verified_hypotheses.empty else None,
        "family_status_counts": {str(key): int(value) for key, value in sorted(counts.items())},
        "evaluated": not bool(not_evaluated_reason),
        "not_evaluated_reason": not_evaluated_reason or None,
        "universe_complete": falsify_result.get("universe_complete") if falsify_result else None,
        "falsification_manifest": falsify_result.get("manifest") if falsify_result else None,
        "inputs": {
            "panel": file_fingerprint(runtime_path(panel_arg)) if panel_arg else None,
            "hypotheses": file_fingerprint(runtime_path(hypotheses_arg)) if hypotheses_arg else None,
            "family_index": file_fingerprint(runtime_path(family_index_arg)) if family_index_arg else None,
            "candidates": file_fingerprint(runtime_path(candidates_arg)) if candidates_arg else None,
            "annotation_rescue": file_fingerprint(runtime_path(rescue_arg)) if rescue_arg else None,
        },
        "safety": [
            "only verified/source-backed targets frozen by the bait panel are forwarded to falsification",
            "only a reference_genome-scoped family index is accepted; interval-scoped fallbacks fail closed",
            "verified bait count of zero is not_evaluated_no_verified_bait, never a miss and never an interval match",
            "annotation member absence alone is held (absent_annotation_only) until a completed sequence search",
            "rank_contribution is always none; this stage writes no rank artifact",
        ],
        "wall_seconds": round(time.perf_counter() - started, 3),
    }
    manifest_path = out / "route_a_interval_family_match_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    result: dict[str, object] = dict(falsify_result)
    result.update({
        "manifest": str(manifest_path), "status": str(status_path),
        "hypothesis_interval_falsification": str(falsification_path) if falsify_result else None,
        "family_status_counts": manifest["family_status_counts"],
        "evaluated": manifest["evaluated"], "not_evaluated_reason": manifest["not_evaluated_reason"],
        "wall_seconds": manifest["wall_seconds"],
    })
    return result


def _rescue_status(rescue_lookup: dict[str, dict[str, object]], hypothesis_id: str,
                   hypotheses: pd.DataFrame) -> str:
    """annotation_rescue_status for one verified hypothesis.

    ``complete`` only when the query has rescue rows and every one of them came
    from a completed provider; ``incomplete`` when any row was not evaluated;
    ``not_run`` when no sequence-level search was supplied for the query.
    """
    if not rescue_lookup:
        return "not_run"
    accession = ""
    gene_name = ""
    if not hypotheses.empty:
        scope = hypotheses[hypotheses["hypothesis_id"].astype(str).eq(hypothesis_id)]
        if not scope.empty:
            row = scope.iloc[0]
            accession = str(row.get("proposed_stable_id", ""))
            gene_name = str(row.get("reference_gene_name", ""))
    projection = rescue_lookup.get(accession.casefold()) or rescue_lookup.get(gene_name.casefold())
    if not projection:
        return "not_run"
    return "complete" if bool(projection.get("evaluated", True)) else "incomplete"
