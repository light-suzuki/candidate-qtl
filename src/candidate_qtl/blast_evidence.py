"""BLAST+ execution layer for the deterministic blast-panel flow (issue #374 / #355).

Issue #374 ("recall-first Top-20 benchmark objective") asks for the blast
evidence path to be a closed, auditable layer. ``blast_panel`` freezes the
bait panel byte-for-byte; this module is the only place that may open the
interval proteins and turn the frozen baits into BLAST queries:

* ``run_panel_blastp`` first verifies the freeze manifest (a one-character
  drift in ``frozen_bait_panel.tsv`` is a ValueError *before* any BLAST tool
  runs -- fail closed), then opens the interval proteins.
* canonical ``verified`` rows and identity-confirmed ``review_only`` rows
  become queries. Review-only projections remain explicitly unverified and
  are emitted only for the separate human-review handoff; other
  ``review_required`` / ``rejected`` rows are skipped and recorded in the
  search manifest.
* a missing BLAST+ binary, a failed ``makeblastdb``, or a failed ``blastp``
  degrades to ``not_evaluated_*`` (never a miss), with evidence output
  restricted to the header and exit code 0.
* every evidence row records ``rank_contribution=none`` and the search
  manifest records ``interval_usage=post_freeze_similarity_only``: blast
  similarity is evidence of similarity only, never causal proof and never a
  gene score. The interval proteins are opened by this stage, so the search
  manifest reports their usage instead of the freeze-time ``interval_accessed``
  flag.

Import policy (pinned by ``tests/test_blast_panel_rank_invariance.py``):
``hardening`` / ``shortlist`` / ``family_index`` / ``llm_pathway`` /
``ranking_policy`` must never be imported here; ``blast_panel`` (the
deterministic freeze core) is allowed.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

import pandas as pd

from . import __version__
from .blast_panel import DEFAULT_THRESHOLDS, VERIFIED, verify_blast_panel_freeze
from .core import ExternalDependencyError, file_fingerprint, runtime_path, run_tool, safe_id, sha256, write_fasta

COMMAND = "blast-panel-search"
SCHEMA_VERSION = "1.0"

EVIDENCE_NAME = "blast_evidence.tsv"
ELIGIBLE_INTERVAL_HITS_NAME = "eligible_interval_hits.tsv"
SEARCH_MANIFEST_NAME = "blast_panel_search_manifest.json"
QUERIES_NAME = "queries.fa"
RAW_HITS_NAME = "panel_vs_interval.tsv"
DB_NAME = "interval_protein_db"

BLAST_EVIDENCE_COLUMNS = [
    "panel_id", "review_status", "provider_status", "query_sequence_sha256",
    "query_length", "top_hit_target_id", "top_hit_identity",
    "top_hit_query_coverage", "top_hit_target_coverage", "top_hit_evalue",
    "top_hit_bitscore", "hit_count", "paralog_count", "ambiguity_status",
    "blastp_version", "db_fingerprint", "interval_proteins_sha256",
    "rank_contribution", "note",
]

# One row per threshold-passing (query, interval-target) pair.  This is a
# deliberately separate artifact from ``blast_evidence.tsv``: the latter is
# the backwards-compatible one-row-per-bait summary, while this table retains
# every eligible competitor for copy/family review.
ELIGIBLE_INTERVAL_HIT_COLUMNS = [
    "panel_id", "target_id", "review_status", "provider_status",
    "ambiguity_status", "pident", "query_coverage", "target_coverage", "evalue", "bitscore",
    "hsp_count", "query_union_start", "query_union_end",
    "target_union_start", "target_union_end", "query_union_segments",
    "target_union_segments", "eligibility", "rank_contribution", "note",
]

PROVIDER_STATUS_VOCAB = (
    "evaluable",
    "not_evaluated_no_verified_bait",
    "not_evaluated_blast_unavailable",
    "not_evaluated_provider_failure",
)

AMBIGUITY_VOCAB = ("unique", "ambiguous_paralog", "no_threshold_hit", "not_evaluated")

# Same outfmt as benchmarks.py (12 columns); the frozen thresholds are the only
# caller-defined filter applied on top.
BLAST_OUTFMT = "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"

POLICY_APPENDIX = "blast hits are evidence of similarity only"


# --------------------------------------------------------------------------- #
# runner selection (mirrors benchmarks._run_interval_blastp)
# --------------------------------------------------------------------------- #

def _wsl_path(path: Path) -> str:
    resolved = path.resolve()
    drive = resolved.drive.rstrip(":").casefold()
    if drive:
        suffix = resolved.as_posix().split(":", 1)[1]
        return f"/mnt/{drive}{suffix}"
    return resolved.as_posix()


def _resolve_runner(runner: str, tools: tuple[str, ...]) -> tuple[str, str | None]:
    """Pick the BLAST runner; ``auto`` = native when the tool is on PATH, else WSL.

    Returns ``(selected, error)`` where error is None when the selected runner
    can actually find every required tool.
    """
    if runner not in {"auto", "native", "wsl"}:
        return "", "blastp runner must be auto, native, or wsl"
    selected = runner
    if runner == "auto":
        selected = "native" if shutil.which(tools[0]) else "wsl"
    if selected == "native":
        missing = [tool for tool in tools if not shutil.which(tool)]
        if missing:
            return selected, f"{', '.join(missing)} not found on PATH"
    if selected == "wsl" and not shutil.which("wsl"):
        return selected, "WSL was not found; install BLAST+ natively or run with WSL"
    return selected, None


def _tool_version(tool: str, runner: str = "native") -> str:
    """Best-effort ``<tool> -version`` first output line; "unknown" on failure."""
    command = (["wsl"] if runner == "wsl" else []) + [tool, "-version"]
    try:
        completed = subprocess.run(command, capture_output=True, text=True, timeout=10, check=False)
    except (OSError, subprocess.TimeoutExpired):
        return "unknown"
    for line in (completed.stdout + completed.stderr).splitlines():
        if line.strip():
            return line.strip()
    return "unknown"


def _blastp_version(runner: str = "native") -> str:
    return _tool_version("blastp", runner)


# --------------------------------------------------------------------------- #
# protein DB
# --------------------------------------------------------------------------- #

def build_protein_db(
    proteins_arg: str,
    db_out_arg: str,
    num_threads: int = 4,
    runner: str = "auto",
) -> dict[str, object]:
    """Build a makeblastdb protein database; failures are reported, never raised.

    A missing makeblastdb binary returns ``{"available": False, "reason":
    "makeblastdb not found on PATH"}``; a non-zero exit reports the stderr.
    """
    proteins = runtime_path(proteins_arg)
    db_out = runtime_path(db_out_arg)
    db_out.parent.mkdir(parents=True, exist_ok=True)
    selected, error = _resolve_runner(runner, ("makeblastdb",))
    if error:
        return {"available": False, "reason": "makeblastdb not found on PATH"}
    command = (["wsl"] if selected == "wsl" else []) + [
        "makeblastdb",
        "-in", _wsl_path(proteins) if selected == "wsl" else str(proteins),
        "-dbtype", "prot",
        "-parse_seqids",
        "-out", _wsl_path(db_out) if selected == "wsl" else str(db_out),
    ]
    try:
        completed = run_tool(command, capture_output=True, text=True, check=False)
    except ExternalDependencyError as error:
        return {"available": False, "reason": f"{error.kind}: {error}"}
    if completed.returncode != 0:
        stderr = (completed.stderr or completed.stdout or "").strip()
        return {"available": False, "reason": stderr or f"makeblastdb exited {completed.returncode}"}
    return {
        "available": True,
        "db": str(db_out),
        "db_fingerprint": sha256(Path(str(db_out) + ".pdb")),
        "source_fingerprint": file_fingerprint(proteins, full_hash=True),
        "tool_version": _tool_version("makeblastdb", selected),
    }


# --------------------------------------------------------------------------- #
# hit parsing and evaluation
# --------------------------------------------------------------------------- #

def _parse_blast_hits(path: Path) -> list[dict[str, object]]:
    """Parse a 12-column BLAST outfmt 6 TSV into a dict list.

    Multiple HSP rows for the same (qseqid, sseqid) are aggregated into one
    target hit.  Query and subject coverage are calculated from coordinate
    unions (overlapping HSPs are counted once); identity is the conservative
    minimum HSP identity, while E-value and bitscore retain the strongest
    evidence (minimum E-value and maximum bitscore).  The HSP count and union
    coordinates remain available to the long-form audit artifact.
    """
    hits: list[dict[str, object]] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != 12:
                raise ValueError(f"{path}:{line_number}: expected 12 BLAST outfmt columns")
            hits.append(
                {
                    "qseqid": fields[0],
                    "sseqid": fields[1],
                    "pident": float(fields[2]),
                    "length": int(fields[3]),
                    "qlen": int(fields[4]),
                    "slen": int(fields[5]),
                    "qstart": int(fields[6]),
                    "qend": int(fields[7]),
                    "sstart": int(fields[8]),
                    "send": int(fields[9]),
                    "evalue": float(fields[10]),
                    "bitscore": float(fields[11]),
                }
            )
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for hit in hits:
        grouped[(str(hit["qseqid"]), str(hit["sseqid"]))].append(hit)

    aggregated: list[dict[str, object]] = []
    for (qseqid, sseqid), group in grouped.items():
        query_intervals = [(int(hit["qstart"]), int(hit["qend"])) for hit in group]
        target_intervals = [(int(hit["sstart"]), int(hit["send"])) for hit in group]
        query_union = _merge_intervals(query_intervals)
        target_union = _merge_intervals(target_intervals)
        # BLAST normally repeats identical qlen/slen on every HSP.  Taking
        # the maximum is conservative for malformed/mixed provider output and
        # prevents coverage inflation when lengths disagree.
        qlen = max(int(hit["qlen"]) for hit in group)
        slen = max(int(hit["slen"]) for hit in group)
        aggregated.append({
            "qseqid": qseqid,
            "sseqid": sseqid,
            "pident": min(float(hit["pident"]) for hit in group),
            "length": sum(int(hit["length"]) for hit in group),
            "qlen": qlen,
            "slen": slen,
            "qstart": min(interval[0] for interval in query_union),
            "qend": max(interval[1] for interval in query_union),
            "sstart": min(interval[0] for interval in target_union),
            "send": max(interval[1] for interval in target_union),
            "evalue": min(float(hit["evalue"]) for hit in group),
            "bitscore": max(float(hit["bitscore"]) for hit in group),
            "hsp_count": len(group),
            "query_union_length": _interval_length(query_union),
            "target_union_length": _interval_length(target_union),
            "query_union_segments": _format_intervals(query_union),
            "target_union_segments": _format_intervals(target_union),
            "query_coverage": _interval_length(query_union) / max(1, qlen),
            "target_coverage": _interval_length(target_union) / max(1, slen),
        })
    return sorted(aggregated, key=lambda hit: (str(hit["qseqid"]), -float(hit["bitscore"]), str(hit["sseqid"])))


def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Return sorted, inclusive coordinate unions without double-counting overlap."""
    normalized = sorted((min(start, end), max(start, end)) for start, end in intervals)
    merged: list[list[int]] = []
    for start, end in normalized:
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def _interval_length(intervals: list[tuple[int, int]]) -> int:
    return sum(end - start + 1 for start, end in intervals)


def _format_intervals(intervals: list[tuple[int, int]]) -> str:
    return ";".join(f"{start}-{end}" for start, end in intervals)


def _threshold_hits(
    hits: list[dict[str, object]],
    thresholds: dict[str, float | int],
    short_query_mode: bool,
) -> list[dict[str, object]]:
    """Apply the frozen thresholds; returns bitscore-descending, tie-broken by sseqid.

    Identity / query coverage / target coverage and the E-value cap always
    apply. Coverage uses the query/subject coordinate union, rather than BLAST
    alignment length, because gapped or multi-HSP alignments can otherwise
    overstate residue coverage. ``short_query_mode`` is retained for call-signature symmetry
    with ``_evaluate_bait``; PAM30 changes search sensitivity, not evidence
    significance requirements.
    """
    passed: list[dict[str, object]] = []
    for hit in hits:
        if float(hit["pident"]) < float(thresholds["min_identity"]):
            continue
        query_coverage = float(hit.get("query_coverage", (abs(int(hit["qend"]) - int(hit["qstart"])) + 1) / max(1, int(hit["qlen"]))))
        target_coverage = float(hit.get("target_coverage", (abs(int(hit["send"]) - int(hit["sstart"])) + 1) / max(1, int(hit["slen"]))))
        if query_coverage < float(thresholds["min_query_coverage"]):
            continue
        if target_coverage < float(thresholds["min_target_coverage"]):
            continue
        if float(hit["evalue"]) > float(thresholds["max_evalue"]):
            continue
        passed.append(hit)
    passed.sort(key=lambda hit: (float(hit["bitscore"]), str(hit["sseqid"])), reverse=True)
    return passed


def _classify_ambiguity(
    hits: list[dict[str, object]],
    thresholds: dict[str, float | int],
) -> tuple[str, dict[str, object] | None]:
    """Classify threshold-passed, bitscore-descending hits.

    0 hits -> (``no_threshold_hit``, None); 1 hit -> (``unique``, hits[0]);
    2+ hits: the second-best bitscore at least ``paralog_tie_bitscore_ratio``
    of the best means ``ambiguous_paralog``, otherwise ``unique``.
    """
    if not hits:
        return "no_threshold_hit", None
    if len(hits) == 1:
        return "unique", hits[0]
    ratio = float(thresholds["paralog_tie_bitscore_ratio"])
    if float(hits[1]["bitscore"]) >= ratio * float(hits[0]["bitscore"]):
        return "ambiguous_paralog", hits[0]
    return "unique", hits[0]


def _evaluate_bait(
    row: dict[str, object],
    hits: list[dict[str, object]],
    thresholds: dict[str, float | int],
    tool_meta: dict[str, object],
    short_query_mode: bool,
) -> dict[str, object]:
    """One evidence row for one verified bait; hits are the raw BLAST hits.

    Threshold filtering runs inside this function, so callers may pass the raw
    per-query hit list. ``top_hit_evalue`` always records the raw evalue.
    Short queries (30-99 aa) run in PAM30 mode but retain the same frozen
    E-value threshold as other queries; their raw E-value is recorded for
    audit.
    """
    passed = _threshold_hits(hits, thresholds, short_query_mode)
    ambiguity, top = _classify_ambiguity(passed, thresholds)
    if top is None:
        target_id, identity, query_coverage, target_coverage, evalue, bitscore = "", "", "", "", "", ""
    else:
        target_id = str(top["sseqid"])
        identity = round(float(top["pident"]), 3)
        query_coverage = round(float(top.get("query_coverage", (abs(int(top["qend"]) - int(top["qstart"])) + 1) / max(1, int(top["qlen"])))), 6)
        target_coverage = round(float(top.get("target_coverage", (abs(int(top["send"]) - int(top["sstart"])) + 1) / max(1, int(top["slen"])))), 6)
        evalue = float(top["evalue"])
        bitscore = float(top["bitscore"])
    return {
        "panel_id": str(row["panel_id"]),
        # This is the upstream projection eligibility, not a statement that a
        # source claim is verified.  Review-only hits are kept for the separate
        # human-review handoff and cannot enter supported_member_union.
        "review_status": str(row.get("review_status", VERIFIED)),
        "provider_status": "evaluable",
        "query_sequence_sha256": str(row["query_sequence_sha256"]),
        "query_length": len(str(row["query_sequence"])),
        "top_hit_target_id": target_id,
        "top_hit_identity": identity,
        "top_hit_query_coverage": query_coverage,
        "top_hit_target_coverage": target_coverage,
        "top_hit_evalue": evalue,
        "top_hit_bitscore": bitscore,
        "hit_count": len(passed),
        "paralog_count": max(0, len(passed) - 1),
        "ambiguity_status": ambiguity,
        "blastp_version": str(tool_meta["blastp_version"]),
        "db_fingerprint": str(tool_meta["db_fingerprint"]),
        "interval_proteins_sha256": str(tool_meta["interval_proteins_sha256"]),
        "rank_contribution": "none",
        "note": "",
    }


def _eligible_interval_hit_rows(
    panel_id: str,
    hits: list[dict[str, object]],
    ambiguity_status: str,
    review_status: str = VERIFIED,
    projection_eligibility: str = "verified_support",
) -> list[dict[str, object]]:
    """Serialize every threshold-passing aggregated hit for copy review."""
    rows: list[dict[str, object]] = []
    for hit in hits:
        rows.append({
            "panel_id": panel_id,
            "target_id": str(hit["sseqid"]),
            "review_status": review_status,
            "provider_status": "evaluable",
            "ambiguity_status": ambiguity_status,
            "pident": round(float(hit["pident"]), 3),
            "query_coverage": round(float(hit["query_coverage"]), 6),
            "target_coverage": round(float(hit["target_coverage"]), 6),
            "evalue": float(hit["evalue"]),
            "bitscore": float(hit["bitscore"]),
            "hsp_count": int(hit["hsp_count"]),
            "query_union_start": int(hit["qstart"]),
            "query_union_end": int(hit["qend"]),
            "target_union_start": int(hit["sstart"]),
            "target_union_end": int(hit["send"]),
            "query_union_segments": str(hit["query_union_segments"]),
            "target_union_segments": str(hit["target_union_segments"]),
            "eligibility": "eligible",
            "rank_contribution": "none",
            "note": f"ambiguity_status={ambiguity_status}",
        })
    return rows


# --------------------------------------------------------------------------- #
# search
# --------------------------------------------------------------------------- #

def run_panel_blastp(
    freeze_arg: str,
    proteins_arg: str,
    out_dir_arg: str,
    num_threads: int = 4,
    runner: str = "auto",
) -> dict[str, object]:
    """BLAST frozen canonical and review-only baits against interval proteins.

    Fail-closed ordering: (1) verify the freeze manifest, (2) read the frozen
    panel and query only canonical verified or identity-confirmed review-only
    rows, (3) never a miss -- a missing
    BLAST+ binary or a provider failure is ``not_evaluated_*`` with exit code
    0 and header-only evidence. The search manifest records
    ``interval_usage=post_freeze_similarity_only``: this stage does open the
    interval proteins, and a blast hit is evidence of similarity only.
    """
    freeze_path = runtime_path(freeze_arg)
    proteins_path = runtime_path(proteins_arg)
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1) verify before anything else: a drifted freeze is a hard error and no
    #    BLAST tool ever runs.
    verified = verify_blast_panel_freeze(freeze_arg)
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    thresholds = dict(freeze.get("thresholds") or DEFAULT_THRESHOLDS)
    case_id = str(freeze.get("case_id") or "")
    if not proteins_path.exists():
        raise ValueError(f"interval protein FASTA does not exist: {proteins_path}")
    interval_fingerprint = file_fingerprint(proteins_path, full_hash=True)
    interval_sha256 = sha256(proteins_path)

    # 2) panel read: canonical verified rows and identity-verified review-only
    # rows become queries.  The latter are explicitly carried as review-only
    # evidence and never become canonical support.
    panel = pd.read_csv(runtime_path(str(verified["panel"])), sep="\t", dtype=str, keep_default_na=False).fillna("")
    records = panel.to_dict("records")
    verified_rows = [record for record in records if (
        str(record["review_status"]) == VERIFIED
        or (str(record.get("projection_eligibility", "")).casefold() == "review_only"
            and str(record.get("verification_status", "")) == VERIFIED)
    )]
    skipped = [str(record["panel_id"]) for record in records if record not in verified_rows]

    evidence_path = out_dir / EVIDENCE_NAME
    eligible_hits_path = out_dir / ELIGIBLE_INTERVAL_HITS_NAME
    manifest_path = out_dir / SEARCH_MANIFEST_NAME
    queries_path = out_dir / QUERIES_NAME
    raw_hits_path = out_dir / RAW_HITS_NAME
    db_path = out_dir / DB_NAME

    def emit(
        provider_status: str,
        rows: list[dict[str, object]],
        counts: dict[str, int],
        eligible_rows: list[dict[str, object]] | None = None,
        reason: str | None = None,
        selected: str = "",
        short_query: bool = False,
        db_fingerprint: str | None = None,
        query_fingerprint: dict[str, object] | None = None,
        blastp_version: str = "",
        groups: dict[str, object] | None = None,
        query_modes: dict[str, str] | None = None,
        raw_hits_fingerprint: str | None = None,
    ) -> dict[str, object]:
        pd.DataFrame(rows, columns=BLAST_EVIDENCE_COLUMNS).to_csv(
            evidence_path, sep="\t", index=False, lineterminator="\n", encoding="utf-8"
        )
        pd.DataFrame(eligible_rows or [], columns=ELIGIBLE_INTERVAL_HIT_COLUMNS).to_csv(
            eligible_hits_path, sep="\t", index=False, lineterminator="\n", encoding="utf-8"
        )
        counts = {**counts, "eligible_interval_hits": len(eligible_rows or [])}
        manifest = {
            "schema_version": SCHEMA_VERSION,
            "command": COMMAND,
            "candidate_qtl_version": __version__,
            "case_id": case_id,
            "runner": selected,
            "blastp_version": blastp_version,
            "short_query_mode": short_query,
            "thresholds": thresholds,
            "provider_status": provider_status,
            "reason": reason,
            "freeze_verified": True,
            "freeze_panel_status": str(freeze.get("panel_status") or ""),
            "freeze": str(freeze_path),
            "db_fingerprint": db_fingerprint,
            "interval_proteins_fingerprint": interval_fingerprint,
            "query_fasta_fingerprint": query_fingerprint,
            "query_modes": query_modes or {},
            "groups": groups or {},
            "raw_hits_sha256": raw_hits_fingerprint,
            "eligible_interval_hits_sha256": sha256(eligible_hits_path),
            "skipped_review_required": skipped,
            "counts": counts,
            "policy": list(freeze.get("policy") or []) + [POLICY_APPENDIX],
            "interval_usage": "post_freeze_similarity_only",
            "outputs": {
                "evidence": str(evidence_path),
                "eligible_interval_hits": str(eligible_hits_path),
                "eligible_hits": str(eligible_hits_path),
                "manifest": str(manifest_path),
                "raw_hits": str(raw_hits_path),
                "queries": str(queries_path),
                "db": str(db_path),
            },
        }
        manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return {
            "evidence": str(evidence_path),
            "eligible_interval_hits": str(eligible_hits_path),
            "eligible_hits": str(eligible_hits_path),
            "manifest": str(manifest_path),
            "provider_status": provider_status,
            "reason": reason,
            "freeze_panel_status": str(freeze.get("panel_status") or ""),
            "rows_evaluated": counts["rows_evaluated"],
            "unique_count": counts["unique"],
            "ambiguous_count": counts["ambiguous_paralog"],
            "no_threshold_count": counts["no_threshold_hit"],
            "skipped_review_required": skipped,
            "interval_usage": "post_freeze_similarity_only",
            "runner": selected,
            "short_query_mode": short_query,
            "groups": groups or {},
            "query_modes": query_modes or {},
            "raw_hits_sha256": raw_hits_fingerprint,
            "eligible_hit_count": len(eligible_rows or []),
            "eligible_interval_hits_sha256": sha256(eligible_hits_path),
            "exit_code": 0,
        }

    empty_counts = {"rows_evaluated": 0, "unique": 0, "ambiguous_paralog": 0, "no_threshold_hit": 0}

    # 3) verified 0 rows: nothing may become a query; never a miss.
    if not verified_rows:
        return emit(
            "not_evaluated_no_verified_bait",
            [],
            empty_counts,
            reason="no verified bait rows in the frozen panel; nothing was searched",
        )

    # 4) BLAST+ availability under the requested runner.
    selected, error = _resolve_runner(runner, ("blastp", "makeblastdb"))
    if error:
        return emit("not_evaluated_blast_unavailable", [], empty_counts, reason=error, selected=selected)

    # 5) Partition by the existing short-query policy.  A mixed panel must not
    # let PAM30 parameters selected for one query affect its long queries.
    groups_by_mode: dict[str, list[dict[str, object]]] = {"standard": [], "short": []}
    for record in verified_rows:
        mode = "short" if len(str(record["query_sequence"])) < 100 else "standard"
        groups_by_mode[mode].append(record)
    active_modes = [mode for mode in ("standard", "short") if groups_by_mode[mode]]
    short_query = "short" in active_modes
    query_modes = {safe_id(str(record["panel_id"])): ("short" if len(str(record["query_sequence"])) < 100 else "standard") for record in verified_rows}
    # Keep the historical single-group paths byte-compatible.  Mixed panels
    # use explicit per-group files and are merged into these canonical outputs.
    query_paths: dict[str, Path] = {}
    raw_paths: dict[str, Path] = {}
    query_fingerprints: dict[str, object] = {}
    for mode in active_modes:
        query_path = queries_path if len(active_modes) == 1 else out_dir / f"queries.{mode}.fa"
        raw_path = raw_hits_path if len(active_modes) == 1 else out_dir / f"panel_vs_interval.{mode}.tsv"
        write_fasta(query_path, [(safe_id(str(record["panel_id"])), str(record["query_sequence"])) for record in groups_by_mode[mode]])
        query_paths[mode] = query_path
        raw_paths[mode] = raw_path
        query_fingerprints[mode] = file_fingerprint(query_path, full_hash=True)
    # The public query artifact remains the complete frozen panel, including
    # mixed panels where BLAST itself consumed the per-mode files.
    if len(active_modes) > 1:
        write_fasta(queries_path, [(safe_id(str(record["panel_id"])), str(record["query_sequence"])) for record in verified_rows])
    query_fingerprint = file_fingerprint(queries_path, full_hash=True)

    # 6) protein DB (provider failure is not_evaluated, not a miss).
    db = build_protein_db(str(proteins_path), str(db_path), num_threads, runner=selected)
    if not db["available"]:
        return emit(
            "not_evaluated_provider_failure",
            [],
            empty_counts,
            reason=db["reason"],
            selected=selected,
            short_query=short_query,
            query_fingerprint=query_fingerprint,
        )

    blastp_version = _blastp_version(selected)

    # 7) Run each active group against the same frozen DB.  Any group failure
    # fails closed for the panel; partial hits must never look like success.
    group_meta: dict[str, object] = {}
    raw_parts: list[str] = []
    for mode in active_modes:
        blast_args = ["blastp", "-query", _wsl_path(query_paths[mode]) if selected == "wsl" else str(query_paths[mode]),
                      "-db", _wsl_path(db_path) if selected == "wsl" else str(db_path), "-outfmt", BLAST_OUTFMT,
                      "-max_target_seqs", "5000", "-num_threads", str(num_threads),
                      "-out", _wsl_path(raw_paths[mode]) if selected == "wsl" else str(raw_paths[mode])]
        if mode == "short":
            blast_args.extend(["-word_size", "2", "-matrix", "PAM30", "-gapopen", "9", "-gapextend", "1", "-comp_based_stats", "0", "-seg", "no", "-evalue", "1000"])
        command = (["wsl"] + blast_args) if selected == "wsl" else blast_args
        meta: dict[str, object] = {
            "mode": mode,
            "query_fasta": str(query_paths[mode]),
            "query_fasta_fingerprint": query_fingerprints[mode],
            "query_fasta_sha256": sha256(query_paths[mode]),
            "command": command,
        }
        try:
            completed = run_tool(command, capture_output=True, text=True, check=False)
        except ExternalDependencyError as error:
            meta.update({"provider_status": "not_evaluated_provider_failure", "reason": f"{error.kind}: {error}"})
            group_meta[mode] = meta
            return emit("not_evaluated_provider_failure", [], empty_counts, reason=f"{mode} group: {error.kind}: {error}", selected=selected, short_query=short_query, db_fingerprint=db["db_fingerprint"], query_fingerprint=query_fingerprint, blastp_version=blastp_version, groups=group_meta, query_modes=query_modes)
        if completed.returncode != 0:
            stderr = (completed.stderr or completed.stdout or "").strip()
            meta.update({"provider_status": "not_evaluated_provider_failure", "reason": stderr or f"blastp exited {completed.returncode}"})
            group_meta[mode] = meta
            return emit("not_evaluated_provider_failure", [], empty_counts, reason=f"{mode} group: {meta['reason']}", selected=selected, short_query=short_query, db_fingerprint=db["db_fingerprint"], query_fingerprint=query_fingerprint, blastp_version=blastp_version, groups=group_meta, query_modes=query_modes)
        if not raw_paths[mode].exists():
            meta.update({"provider_status": "not_evaluated_provider_failure", "reason": "blastp exited 0 without creating raw output"})
            group_meta[mode] = meta
            return emit("not_evaluated_provider_failure", [], empty_counts, reason=f"{mode} group: {meta['reason']}", selected=selected, short_query=short_query, db_fingerprint=db["db_fingerprint"], query_fingerprint=query_fingerprint, blastp_version=blastp_version, groups=group_meta, query_modes=query_modes)
        raw_text = raw_paths[mode].read_text(encoding="utf-8")
        meta.update({"provider_status": "evaluable", "raw_hits": str(raw_paths[mode]), "raw_hits_sha256": sha256(raw_paths[mode])})
        group_meta[mode] = meta
        raw_parts.append(raw_text)
    # Canonical merged raw output: stable group order, retaining every HSP row.
    if len(active_modes) > 1:
        merged_lines = [line for text in raw_parts for line in text.splitlines() if line.strip()]
        raw_hits_path.write_text(("\n".join(merged_lines) + "\n") if merged_lines else "", encoding="utf-8")
    raw_hits_fingerprint = sha256(raw_hits_path)
    for meta in group_meta.values():
        if isinstance(meta, dict):
            meta["merged_raw_hits_sha256"] = raw_hits_fingerprint

    # 8) evaluate: threshold-filter, classify, emit one row per verified bait.
    #    A provider that writes unparseable raw hits is not a miss: degrade to
    #    not_evaluated_provider_failure with header-only evidence and exit 0.
    try:
        parsed_hits = _parse_blast_hits(raw_hits_path)
    except (ValueError, RuntimeError, OSError) as error:
        return emit(
            "not_evaluated_provider_failure",
            [],
            empty_counts,
            reason=f"unparseable blastp raw hits: {error}",
            selected=selected,
            short_query=short_query,
            db_fingerprint=db["db_fingerprint"],
            query_fingerprint=query_fingerprint,
            blastp_version=blastp_version,
            groups=group_meta, query_modes=query_modes, raw_hits_fingerprint=raw_hits_fingerprint,
        )
    hits_by_query: dict[str, list[dict[str, object]]] = defaultdict(list)
    for hit in parsed_hits:
        hits_by_query[str(hit["qseqid"])].append(hit)

    tool_meta = {
        "blastp_version": blastp_version,
        "db_fingerprint": db["db_fingerprint"],
        "interval_proteins_sha256": interval_sha256,
    }
    rows: list[dict[str, object]] = []
    eligible_interval_rows: list[dict[str, object]] = []
    for record in verified_rows:
        panel_id = safe_id(str(record["panel_id"]))
        hits = hits_by_query.get(panel_id, [])
        mode = query_modes[panel_id]
        row = _evaluate_bait(record, hits, thresholds, tool_meta, mode == "short")
        rows.append(row)
        passed = _threshold_hits(hits, thresholds, mode == "short")
        eligible_interval_rows.extend(
            _eligible_interval_hit_rows(
                str(record["panel_id"]), passed, str(row["ambiguity_status"]),
                str(record.get("review_status", VERIFIED)),
                str(record.get("projection_eligibility", "verified_support")),
            )
        )
    counts = {
        "rows_evaluated": len(rows),
        "unique": sum(1 for row in rows if row["ambiguity_status"] == "unique"),
        "ambiguous_paralog": sum(1 for row in rows if row["ambiguity_status"] == "ambiguous_paralog"),
        "no_threshold_hit": sum(1 for row in rows if row["ambiguity_status"] == "no_threshold_hit"),
    }
    return emit(
        "evaluable",
        rows,
        counts,
        eligible_rows=eligible_interval_rows,
        selected=selected,
        short_query=short_query,
        db_fingerprint=db["db_fingerprint"],
        query_fingerprint=query_fingerprint,
        blastp_version=blastp_version,
        groups=group_meta, query_modes=query_modes, raw_hits_fingerprint=raw_hits_fingerprint,
    )
