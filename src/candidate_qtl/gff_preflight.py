from __future__ import annotations

import json
import shlex
import shutil
from collections import Counter
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from .core import file_fingerprint, runtime_path, run_tool
from .workflow_extensions import extract_multi


def _feature_counts(path: Path) -> Counter[str]:
    counts: Counter[str] = Counter()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"): continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) == 9: counts[fields[2].casefold()] += 1
    return counts


def _sequences(path: Path | None) -> dict[str, str]:
    if not path or not path.exists(): return {}
    with path.open(encoding="utf-8") as handle:
        return {record.id: str(record.seq).upper() for record in SeqIO.parse(handle, "fasta")}


def gff_preflight(fasta_arg: str, gff_arg: str, intervals: list[str], out_dir_arg: str,
                  seqid_map: str | None = None, normalized_gff_arg: str | None = None,
                  agat_cds_arg: str | None = None, gffread_cds_arg: str | None = None,
                  execute: bool = False) -> dict[str, object]:
    fasta, raw, out = runtime_path(fasta_arg), runtime_path(gff_arg), runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    normalized = runtime_path(normalized_gff_arg) if normalized_gff_arg else out / "agat_normalized.gff3"
    agat_cds = runtime_path(agat_cds_arg) if agat_cds_arg else out / "agat_cds.fa"
    gffread_cds = runtime_path(gffread_cds_arg) if gffread_cds_arg else out / "gffread_cds.fa"
    commands = {
        "agat_normalize": f"agat_convert_sp_gxf2gxf.pl -g {shlex.quote(str(raw))} -o {shlex.quote(str(normalized))}",
        "agat_extract_cds": f"agat_sp_extract_sequences.pl -g {shlex.quote(str(normalized))} -f {shlex.quote(str(fasta))} -t cds -o {shlex.quote(str(agat_cds))}",
        "gffread_extract_cds": f"gffread {shlex.quote(str(normalized if normalized.exists() else raw))} -g {shlex.quote(str(fasta))} -x {shlex.quote(str(gffread_cds))}",
    }
    versions: dict[str, str] = {}
    if execute:
        for tool in ("agat_convert_sp_gxf2gxf.pl", "agat_sp_extract_sequences.pl", "gffread"):
            if not shutil.which(tool): raise ValueError(f"required optional GFF cross-check tool not found: {tool}")
        for name in ("agat_normalize", "agat_extract_cds", "gffread_extract_cds"):
            completed = run_tool(commands[name], shell=True, text=True, capture_output=True)
            if completed.returncode: raise RuntimeError(f"{name} failed: {completed.stderr.strip()}")
        for tool in ("agat_convert_sp_gxf2gxf.pl", "gffread"):
            completed = run_tool([tool, "--version"], text=True, capture_output=True)
            versions[tool] = (completed.stdout or completed.stderr).strip().splitlines()[0] if (completed.stdout or completed.stderr).strip() else "unknown"

    internal = out / "candidate_qtl_extract"
    result = extract_multi(str(fasta), str(raw), intervals, str(internal), seqid_map, False, None, 0)
    transcripts = pd.read_csv(internal / "transcripts.tsv", sep="\t", dtype=str).fillna("") if (internal / "transcripts.tsv").exists() else pd.DataFrame()
    internal_seq = _sequences(internal / "candidate_transcript_cds.fa")
    provider_sequences = {"agat": _sequences(agat_cds if agat_cds.exists() else None), "gffread": _sequences(gffread_cds if gffread_cds.exists() else None)}
    # A missing per-transcript provider sequence only signals a problem when a
    # cross-check was actually requested (execute, or provider CDS supplied).
    # Without one, every transcript would be "missing" simply because no
    # comparison ran, which is plan-mode noise, not a review item.
    cross_check_expected = bool(execute or agat_cds_arg or gffread_cds_arg)
    rows = []
    for transcript_id in sorted(set(internal_seq) | set(provider_sequences["agat"]) | set(provider_sequences["gffread"])):
        internal_value = internal_seq.get(transcript_id, "")
        agat_value = provider_sequences["agat"].get(transcript_id, "")
        gffread_value = provider_sequences["gffread"].get(transcript_id, "")
        available = [value for value in (internal_value, agat_value, gffread_value) if value]
        sequence_match = len(set(available)) <= 1 if len(available) >= 2 else None
        provider_missing = cross_check_expected and len(available) < 2
        tx = transcripts[transcripts["transcript_id"].eq(transcript_id)] if not transcripts.empty else pd.DataFrame()
        rows.append({"transcript_id":transcript_id,"gene_id":str(tx.iloc[0].get("gene_id", "")) if not tx.empty else "",
                     "candidate_qtl_length":len(internal_value) if internal_value else "","agat_length":len(agat_value) if agat_value else "",
                     "gffread_length":len(gffread_value) if gffread_value else "","sequence_match":sequence_match if sequence_match is not None else "not_compared",
                     "phase_status":str(tx.iloc[0].get("cds_phase_status", "")) if not tx.empty else "",
                     "review_required":bool((len(available) >= 2 and not sequence_match) or provider_missing),
                     "reason":"sequence_or_length_mismatch" if len(available) >= 2 and not sequence_match else ("provider_missing_transcript" if provider_missing else "")})
    concordance = pd.DataFrame(rows)
    concordance.to_csv(out / "gff_concordance.tsv", sep="\t", index=False)
    reviews = concordance[concordance["review_required"].eq(True)] if not concordance.empty else concordance
    reviews.to_csv(out / "review_required.tsv", sep="\t", index=False)
    mismatch_count = int((concordance["reason"] == "sequence_or_length_mismatch").sum()) if not concordance.empty else 0
    provider_missing_count = int((concordance["reason"] == "provider_missing_transcript").sum()) if not concordance.empty else 0
    raw_counts, normalized_counts = _feature_counts(raw), (_feature_counts(normalized) if normalized.exists() else Counter())
    feature_types = sorted(set(raw_counts) | set(normalized_counts))
    pd.DataFrame([{"feature_type":key,"raw_count":raw_counts[key],"normalized_count":normalized_counts[key] if normalized.exists() else "",
                   "delta":normalized_counts[key]-raw_counts[key] if normalized.exists() else ""} for key in feature_types]).to_csv(out / "gff_feature_diff.tsv", sep="\t", index=False)
    manifest = {"schema_version":"1.0","raw_gff":file_fingerprint(raw, False),
                "normalized_gff":file_fingerprint(normalized, False) if normalized.exists() else None,
                "reference_fasta":file_fingerprint(fasta, False),"intervals":intervals,"seqid_map":seqid_map,
                "candidate_extraction_gff":str(raw),"cross_check_only":True,"execute":execute,"commands":commands,"tool_versions":versions,
                "candidate_count":result["candidate_count"],"transcript_count":len(transcripts),"compared_transcripts":int((concordance["sequence_match"] != "not_compared").sum()) if not concordance.empty else 0,
                "discordant_transcripts":mismatch_count,"provider_missing_transcripts":provider_missing_count,"review_required_transcripts":len(reviews),"cross_check_requested":cross_check_expected,
                "interpretation":"AGAT/gffread agreement is parser robustness, not evidence for causality. Normalized GFF never replaces raw input automatically."}
    manifest_path = out / "gff_preflight_manifest.json"; manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    (out / "run_gff_crosscheck.sh").write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + "\n".join(commands.values()) + "\n", encoding="utf-8")
    return {"manifest":str(manifest_path),"concordance":str(out / "gff_concordance.tsv"),"review_required":str(out / "review_required.tsv"),"discordant_transcripts":mismatch_count,"review_required_transcripts":len(reviews)}
