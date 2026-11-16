"""Residue-aligned BLAST evidence for Route B, not functional/orthology support.

Policy query_residue_consensus_v1 uses the same unique query residues for
identity and coverage. A subject gap is a mismatch; query-gap columns do not
represent query residues. Disagreeing mappings at an overlapping residue count
as nonidentical. This deliberately is not BLAST's alignment-column pident.
"""
from __future__ import annotations

import csv
import hashlib
import math
import subprocess
from collections.abc import Callable
from pathlib import Path
from typing import Any

from .route_b_identity import METADATA_COLUMNS, merge_header_metadata

IDENTITY_POLICY = "query_residue_consensus_v1"
RETRIEVAL_LIMIT = 5
OUTFMT = ("6 qseqid sseqid pident length qlen slen qstart qend sstart send "
          "evalue bitscore qseq sseq stitle staxids sscinames")
EVIDENCE_COLUMNS = [
    "gene_id", "raw_reference_gene_id", "reference_gene_id", "percent_identity",
    "max_hsp_identity", "query_coverage", "covered_query_residues",
    "identical_query_residues", "conflicting_query_residues", "query_length",
    "evalue", "bitscore", "identity_policy", "threshold_pass", *METADATA_COLUMNS,
]


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def aggregate_hsps(path: Path, max_evalue: float) -> tuple[list[dict[str, Any]], list[str]]:
    """Read the extended outfmt, preserving a legacy 13-column projection.

    The residue observation set is independent of row order, duplicate HSPs,
    and segmentation of the same alignment into multiple HSPs. Raw input is
    never overwritten. Malformed output is a technical failure, not a no-hit.
    """
    pairs: dict[tuple[str, str], dict[str, Any]] = {}
    legacy: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for lineno, line in enumerate(handle, 1):
            if not line.strip():
                continue
            f = line.rstrip("\r\n").split("\t", 16)
            try:
                if len(f) not in {15, 17} or not f[0] or not f[1]:
                    raise ValueError("expected 15 fields including qseq/sseq, or 17 with taxonomy")
                pident, length, qlen, slen = float(f[2]), int(f[3]), int(f[4]), int(f[5])
                qs, qe, ss, se = map(int, f[6:10])
                evalue, bitscore = float(f[10]), float(f[11])
                qseq, sseq = f[12].upper(), f[13].upper()
                if not all(math.isfinite(x) for x in (pident, evalue, bitscore)):
                    raise ValueError("nonfinite measurement")
                if not (0 <= pident <= 100 and evalue >= 0 and bitscore >= 0):
                    raise ValueError("invalid measurement range")
                if not (1 <= qs <= qe <= qlen and 1 <= ss <= se <= slen):
                    raise ValueError("invalid blastp coordinates")
                if length != len(qseq) or length != len(sseq):
                    raise ValueError("alignment length mismatch")
                if (len(qseq.replace("-", "")) != qe - qs + 1
                        or len(sseq.replace("-", "")) != se - ss + 1):
                    raise ValueError("alignment/coordinate mismatch")
                if any(q == s == "-" for q, s in zip(qseq, sseq)):
                    raise ValueError("double-gap alignment column")
                legacy.add("\t".join(f[:12] + [f[14]]))
                if evalue > max_evalue:
                    continue
                pair = pairs.setdefault((f[0], f[1]), {
                    "gene_id": f[0], "raw_reference_gene_id": f[1],
                    "query_length": qlen, "subject_length": slen,
                    "max_hsp_identity": pident, "evalue": evalue,
                    "bitscore": bitscore, "observations": {}, "headers": set(),
                })
                if pair["query_length"] != qlen or pair["subject_length"] != slen:
                    raise ValueError("inconsistent sequence length within pair")
                pair["max_hsp_identity"] = max(pair["max_hsp_identity"], pident)
                pair["bitscore"] = max(pair["bitscore"], bitscore)
                pair["evalue"] = min(pair["evalue"], evalue)
                pair["headers"].add((f[14], f[15] if len(f) == 17 else "", f[16] if len(f) == 17 else ""))
                qp, sp = qs, ss
                for q, s in zip(qseq, sseq):
                    if q != "-":
                        pair["observations"].setdefault(qp, set()).add(
                            (sp if s != "-" else None, q, s))
                        qp += 1
                    if s != "-":
                        sp += 1
            except (ValueError, OverflowError) as exc:
                raise ValueError(f"invalid BLAST evidence at {path}:{lineno}: {exc}") from exc
    rows = []
    for key in sorted(pairs):
        pair = pairs[key]
        observations = pair.pop("observations")
        pair.update(merge_header_metadata(pair.pop("headers")))
        identical = sum(len(values) == 1 and next(iter(values))[1] == next(iter(values))[2]
                        for values in observations.values())
        covered = len(observations)
        pair.update({
            "covered_query_residues": covered,
            "identical_query_residues": identical,
            "conflicting_query_residues": sum(len(values) > 1 for values in observations.values()),
            "percent_identity": 100.0 * identical / covered,
            "query_coverage": covered / pair["query_length"],
            "identity_policy": IDENTITY_POLICY,
        })
        rows.append(pair)
    return rows, sorted(legacy)


def run_local_blastp(
    proteins: Path, database: Path, out_dir: Path, threads: int,
    min_identity: float, min_coverage: float, max_evalue: float, *,
    run_tool: Callable[..., Any], find_binary: Callable[[str], str | None],
    normalise_subject: Callable[[str], str],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Run BLAST and retain all eligible retrieved references (#712/#713)."""
    if not (0 <= min_identity <= 100 and 0 <= min_coverage <= 1
            and math.isfinite(max_evalue) and max_evalue >= 0):
        raise ValueError("invalid homology thresholds")
    binary = find_binary("blastp")
    if not binary:
        return [], {"available": False, "reason": "blastp is not on PATH",
                    "evaluation_state": "not_evaluated_provider_unavailable",
                    "missing_provider_is_not_negative_evidence": True}
    raw = out_dir / "route_b_local_blastp_alignment.tsv"
    output = out_dir / "route_b_local_blastp.tsv"
    command = [binary, "-query", str(proteins), "-db", str(database),
               "-evalue", str(max_evalue), "-max_target_seqs", str(RETRIEVAL_LIMIT),
               "-num_threads", str(threads), "-outfmt", OUTFMT, "-out", str(raw)]
    try:
        completed = run_tool(command, text=True, capture_output=True)
    except (OSError, subprocess.SubprocessError) as error:
        return [], {"available": False, "evaluation_state": "provider_failure",
                    "reason": str(error), "database": str(database),
                    "missing_provider_is_not_negative_evidence": True}
    if completed.returncode:
        return [], {"available": False, "evaluation_state": "provider_failure",
                    "reason": f"local blastp failed: {completed.stderr.strip()}",
                    "returncode": completed.returncode, "database": str(database),
                    "missing_provider_is_not_negative_evidence": True}
    try:
        pairs, legacy = aggregate_hsps(raw, max_evalue)
    except (OSError, ValueError) as error:
        return [], {"available": False, "evaluation_state": "failed_output_parse",
                    "reason": str(error), "database": str(database),
                    "missing_provider_is_not_negative_evidence": True}
    # Keep the existing family-index input layout; extended alignments have
    # their own immutable-in-the-run source path and manifest fingerprint.
    output.write_text("".join(line + "\n" for line in legacy), encoding="utf-8")
    eligible: list[dict[str, Any]] = []
    for pair in pairs:
        pair["reference_gene_id"] = normalise_subject(pair["raw_reference_gene_id"])
        pair["hit_accession"] = pair["reference_gene_id"]
        pair["threshold_pass"] = (pair["percent_identity"] >= min_identity
                                  and pair["query_coverage"] >= min_coverage)
        if not pair["threshold_pass"]:
            continue
        eligible.append(pair)
    rows = []
    for pair in sorted(eligible, key=lambda row: (row["gene_id"], row["reference_gene_id"], row["raw_reference_gene_id"])):
        rows.append({
            **{key: pair[key] for key in EVIDENCE_COLUMNS if key in pair},
            "orthogroup": "", "relationship": "best_hit_homology",
            "provider": f"local_blastp:{database.name}",
        })
    evidence = out_dir / "route_b_local_blastp_pair_evidence.tsv"
    with evidence.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=EVIDENCE_COLUMNS,
                                delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(pairs)
    meta = {
        "available": True, "evaluation_state": "completed", "tool": binary, "database": str(database),
        "hits_table": str(output), "command": " ".join(command),
        "thresholds": {"min_percent_identity": min_identity,
                       "min_query_coverage": min_coverage, "max_evalue": max_evalue},
        "genes_with_a_usable_hit": len({row["gene_id"] for row in rows}),
        "retained_reference_pair_count": len(rows),
        "relationship_claim": "best_hit_homology; not evidence of one-to-one orthology or trait support",
        "subject_id_normalization": "UniProt-style sp|ACC|NAME headers are reduced to ACC",
        "query_coverage_basis": "merged_hsps_per_query_subject_pair",
        "identity_policy": IDENTITY_POLICY,
        "identity_denominator": "unique covered query residues, including subject-gap mismatches",
        "identity_overlap_policy": "conflicting residue mappings count as nonidentical; duplicate mappings count once",
        "query_gap_policy": "query-gap columns have no query residue and are excluded from identity denominator",
        "max_hsp_identity_is_gate_input": False,
        "retrieval_max_target_seqs": RETRIEVAL_LIMIT,
        "hit_metadata_columns": METADATA_COLUMNS,
        "metadata_policy": "explicit database/header values only; conflicting records retained",
        "species_quota": None,
        "retention_policy": "all_threshold_passing_retrieved_subjects_v1",
        "retrieval_limit_is_not_retention_limit": True,
        "legacy_relationship_note": "best_hit_homology names local homology, not a sole best reference",
        "raw_alignment_table": str(raw), "raw_alignment_sha256": _digest(raw),
        "pair_evidence_table": str(evidence), "pair_evidence_sha256": _digest(evidence),
        "hits_table_sha256": _digest(output),
    }
    return rows, meta
