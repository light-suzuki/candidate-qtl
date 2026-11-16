"""Recover significant partial alignments as review context, never support."""
from __future__ import annotations

import json
import math
from pathlib import Path

from Bio import SeqIO

from .blast_evidence import _parse_blast_hits, _threshold_hits
from .blast_panel import verify_blast_panel_freeze
from .core import runtime_path, safe_id, sha256


def partial_review_hits(evidence_arg: str | None, panel_records: list[dict], out: Path) -> dict[str, list[dict]]:
    if not evidence_arg:
        return {}
    source = runtime_path(evidence_arg)
    manifest_path = source.parent / "blast_panel_search_manifest.json"
    if not manifest_path.exists():
        return {}  # legacy compact-only artifacts cannot reconstruct partial hits
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("provider_status") != "evaluable" or not manifest.get("raw_hits_sha256"):
        return {}
    raw = runtime_path(manifest["outputs"]["raw_hits"])
    if sha256(raw) != manifest["raw_hits_sha256"]:
        raise ValueError("partial review raw BLAST hash mismatch")
    freeze_path = runtime_path(manifest["freeze"])
    verified = verify_blast_panel_freeze(str(freeze_path))
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    thresholds = freeze["thresholds"]
    if thresholds != manifest["thresholds"]:
        raise ValueError("partial review thresholds differ from frozen search")
    proteins = manifest["interval_proteins_fingerprint"]
    if not proteins.get("sha256") or sha256(runtime_path(proteins["path"])) != proteins["sha256"]:
        raise ValueError("partial review interval proteins hash mismatch")
    target_lengths = {}
    with runtime_path(proteins["path"]).open(encoding="utf-8") as handle:
        for sequence in SeqIO.parse(handle, "fasta"):
            if sequence.id in target_lengths:
                raise ValueError("partial review duplicate interval protein ID")
            target_lengths[sequence.id] = len(sequence.seq)
    import pandas as pd
    frozen = pd.read_csv(runtime_path(verified["panel"]), sep="\t", dtype=str).fillna("")
    valid = {safe_id(str(r["panel_id"])): r for r in frozen.to_dict("records")
             if r.get("verification_status") == "verified" and r.get("review_status") in {"verified", "review_required"}}
    expected = {safe_id(str(r.get("bait_id", ""))): r for r in panel_records}
    result: dict[str, list[dict]] = {}
    for hit in _parse_blast_hits(raw):
        query = str(hit["qseqid"])
        record, panel = valid.get(query), expected.get(query)
        if record is None or panel is None:
            continue
        if not panel.get("sequence_sha256") or str(record.get("query_sequence_sha256", "")) != str(panel["sequence_sha256"]):
            raise ValueError("partial review bait sequence differs from frozen search")
        if not all(math.isfinite(float(hit[key])) for key in ("pident", "evalue", "bitscore", "query_coverage", "target_coverage")):
            raise ValueError("partial review non-finite BLAST measurement")
        if _threshold_hits([hit], thresholds, False):
            continue  # leave strict-support eligibility unchanged
        if (float(hit["pident"]) < float(thresholds["min_identity"])
                or float(hit["evalue"]) > float(thresholds["max_evalue"])):
            continue
        # Reuse the already-frozen minimum query length as the minimum span on
        # EACH sequence. No case/family-specific thresholds or candidate cap.
        qspan = int(hit["query_union_length"])
        tspan = int(hit["target_union_length"])
        if min(qspan, tspan) < int(thresholds["min_query_length"]):
            continue
        if int(hit["qlen"]) != len(record["query_sequence"]):
            raise ValueError("partial review query length mismatch")
        if target_lengths.get(str(hit["sseqid"])) != int(hit["slen"]):
            raise ValueError("partial review target length or identity mismatch")
        if not (0 < float(hit["query_coverage"]) <= 1 and 0 < float(hit["target_coverage"]) <= 1):
            raise ValueError("partial review invalid alignment coverage")
        result.setdefault(str(panel["bait_id"]), []).append({
            "target_id": str(hit["sseqid"]), "identity": str(hit["pident"]),
            "query_coverage": str(hit["query_coverage"]), "target_coverage": str(hit["target_coverage"]),
            "evalue": str(hit["evalue"]), "bitscore": str(hit["bitscore"]),
            "projection_basis": "partial_alignment_review_only",
        })
    payload = {"state": "review_only", "rank_contribution": "none", "search_manifest_sha256": sha256(manifest_path),
               "raw_hits_sha256": sha256(raw), "freeze_sha256": sha256(freeze_path), "hits": result,
               "interval_proteins_sha256": proteins["sha256"],
               "policy": "significant partial sequence match is not full-length homology or trait support; preserve competing copies for review"}
    payload["inputs"] = {name: {"path": str(path), "sha256": sha256(path)} for name, path in (
        ("search_manifest", manifest_path), ("raw_hits", raw), ("freeze", freeze_path),
        ("frozen_panel", runtime_path(verified["panel"])), ("interval_proteins", runtime_path(proteins["path"])),
    )}
    (out / "partial_alignment_review.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result
