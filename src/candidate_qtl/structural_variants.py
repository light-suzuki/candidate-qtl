from __future__ import annotations

import json
import shlex
from pathlib import Path

import pandas as pd

from .core import file_fingerprint, parse_interval, runtime_path


STRUCTURAL_TYPES = {"INV","TRANS","INVTR","DUP","INVDP","NOTAL","CPG","CPL","HDR","TDM","INS","DEL"}


def _syri_calls(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str).fillna("")
    if raw.shape[1] < 11: raise ValueError("SyRI TSV requires at least 11 columns")
    raw = raw.iloc[:, :12]; raw.columns = ["chrom","start","end","ref_sequence","query_sequence","query_chrom","query_start","query_end","call_id","parent_id","svtype","copy_status"]
    raw = raw[raw["svtype"].isin(STRUCTURAL_TYPES)].copy()
    raw["provider"] = "syri"; raw["support"] = "assembly_alignment"; raw["mapq"] = ""; raw["sample"] = "assembly_pair"
    raw["analysis_chain"] = "assembly:minimap2:syri"; raw["flags"] = ""; raw["source_file"] = str(path)
    return raw


def _provider_calls(specs: list[str] | None) -> list[pd.DataFrame]:
    frames = []
    for spec in specs or []:
        if "=" not in spec: raise ValueError("--call-table must be PROVIDER=PATH")
        provider, value = spec.split("=", 1); path = runtime_path(value)
        frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        required = {"chrom","start","end","svtype"}
        if not required.issubset(frame.columns): raise ValueError(f"{path}: call table missing {', '.join(sorted(required-set(frame.columns)))}")
        frame["provider"] = provider
        for column, default in (("call_id", ""),("support", ""),("mapq", ""),("sample", ""),("analysis_chain", f"read:{provider}"),("flags", ""),("source_file", str(path))):
            if column not in frame.columns: frame[column] = default
        frames.append(frame)
    return frames


def _overlap_class(row: pd.Series, candidates: pd.DataFrame, qtl: tuple[str,int,int]) -> tuple[str, str]:
    chrom, start, end = str(row["chrom"]), int(row["start"]), int(row["end"])
    if chrom != qtl[0] or end < qtl[1] or start > qtl[2]: return "outside_qtl", ""
    for candidate in candidates.itertuples():
        if str(getattr(candidate, "chrom", "")) != chrom: continue
        gene_start, gene_end = int(candidate.start), int(candidate.end)
        if not (end < gene_start or start > gene_end):
            for label, left_name, right_name in (("five_prime_UTR","five_utr_start","five_utr_end"),("three_prime_UTR","three_utr_start","three_utr_end")):
                left, right = getattr(candidate, left_name, ""), getattr(candidate, right_name, "")
                if str(left) and str(right) and not (end < int(left) or start > int(right)): return label, str(candidate.gene_id)
            return "gene_body", str(candidate.gene_id)
        strand = str(getattr(candidate, "strand", "+")); promoter_start, promoter_end = ((max(1, gene_start-2000), gene_start-1) if strand != "-" else (gene_end+1, gene_end+2000))
        if not (end < promoter_start or start > promoter_end): return "promoter", str(candidate.gene_id)
    return "intergenic_qtl", ""


def sv_plan(reference_arg: str, query_arg: str, qtl_interval: str, out_dir_arg: str,
            candidates_arg: str | None = None, syri_arg: str | None = None,
            call_specs: list[str] | None = None, breakpoint_tolerance: int = 100,
            decision_arg: str | None = None) -> dict[str, object]:
    if breakpoint_tolerance < 0: raise ValueError("breakpoint_tolerance must be non-negative")
    reference, query, out = runtime_path(reference_arg), runtime_path(query_arg), runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True); qtl = parse_interval(qtl_interval)
    trigger = "manual"
    if decision_arg:
        decision = json.loads(runtime_path(decision_arg).read_text(encoding="utf-8")); trigger = str(decision.get("mode", ""))
        if trigger != "broaden_to_regulatory_sv_te_cnv":
            raise ValueError("--decision must request broaden_to_regulatory_sv_te_cnv")
    candidates = pd.read_csv(runtime_path(candidates_arg), sep="\t", dtype=str).fillna("") if candidates_arg else pd.DataFrame(columns=["gene_id","chrom","start","end","strand"])
    frames = []
    if syri_arg: frames.append(_syri_calls(runtime_path(syri_arg)))
    frames.extend(_provider_calls(call_specs))
    columns = ["chrom","start","end","svtype","provider","call_id","support","mapq","sample","analysis_chain","flags","source_file"]
    calls = pd.concat(frames, ignore_index=True, sort=False).fillna("") if frames else pd.DataFrame(columns=columns)
    for column in columns:
        if column not in calls.columns: calls[column] = ""
    if not calls.empty:
        calls[["start","end"]] = calls[["start","end"]].astype(int)
        overlaps = [_overlap_class(row, candidates, qtl) for _, row in calls.iterrows()]
        calls["overlap_class"] = [item[0] for item in overlaps]; calls["gene_id"] = [item[1] for item in overlaps]
        numeric_mapq = pd.to_numeric(calls["mapq"], errors="coerce")
        risk = calls["flags"].str.contains("gap|unplaced|low|repeat", case=False, regex=True) | numeric_mapq.lt(20).fillna(False)
        calls["review_required"] = risk
    else:
        calls["overlap_class"] = []; calls["gene_id"] = []; calls["review_required"] = []
    calls.to_csv(out / "canonical_sv.tsv", sep="\t", index=False)

    concordance_rows = []
    used: set[int] = set()
    for index, call in calls.iterrows():
        if index in used: continue
        same_chrom = calls[calls["chrom"].eq(call["chrom"])]
        close = same_chrom[(same_chrom["start"].sub(call["start"]).abs() <= breakpoint_tolerance) & (same_chrom["end"].sub(call["end"]).abs() <= breakpoint_tolerance)]
        if len(close) == 1:
            broad = same_chrom[(same_chrom["svtype"].eq(call["svtype"])) &
                               (same_chrom["start"].sub(call["start"]).abs() <= max(1000, breakpoint_tolerance * 10)) &
                               (same_chrom["end"].sub(call["end"]).abs() <= max(1000, breakpoint_tolerance * 10)) &
                               (~same_chrom.index.isin(used))]
            if any(provider == "syri" for provider in broad["provider"]) and any(provider != "syri" for provider in broad["provider"]): close = broad
        used.update(close.index)
        providers = sorted(set(close["provider"])); types = sorted(set(close["svtype"])); chains = sorted(set(close["analysis_chain"]))
        assembly = any(provider == "syri" for provider in providers); read = any(provider != "syri" for provider in providers)
        max_delta = max(int(close["start"].max()-close["start"].min()),int(close["end"].max()-close["end"].min()))
        if assembly and read and len(types) == 1 and max_delta <= breakpoint_tolerance: status = "assembly_read_concordant"
        elif assembly and read and len(types) == 1: status = "breakpoint_discordant"
        elif assembly and read: status = "type_discordant"
        elif assembly: status = "assembly_only"
        else: status = "read_only"
        concordance_rows.append({"cluster_id":f"SVCLUSTER{len(concordance_rows)+1}","chrom":call["chrom"],"start_min":int(close["start"].min()),"start_max":int(close["start"].max()),
                                 "end_min":int(close["end"].min()),"end_max":int(close["end"].max()),"svtypes":";".join(types),"providers":";".join(providers),
                                 "analysis_chains":";".join(chains),"status":status,"type_match":len(types)==1,"breakpoint_max_delta":max_delta,
                                 "independent_biological_evidence":False,"review_required":status in {"type_discordant","breakpoint_discordant"} or bool(close["review_required"].any())})
    concordance = pd.DataFrame(concordance_rows, columns=["cluster_id","chrom","start_min","start_max","end_min","end_max","svtypes","providers","analysis_chains","status","type_match","breakpoint_max_delta","independent_biological_evidence","review_required"])
    concordance.to_csv(out / "sv_concordance.tsv", sep="\t", index=False)
    concordance[concordance["review_required"].eq(True)].to_csv(out / "review_required.tsv", sep="\t", index=False)
    evidence = calls[calls["gene_id"].ne("")].groupby("gene_id")["svtype"].apply(lambda values:";".join(sorted(set(values)))).reset_index() if not calls.empty else pd.DataFrame(columns=["gene_id","svtype"])
    if not evidence.empty:
        evidence = evidence.rename(columns={"svtype":"value"}); evidence["evidence_code"]="structural_variant_context"; evidence["direction"]="context"; evidence["confidence"]="low"
        evidence["source_type"]="sv_concordance"; evidence["source_id"]=str(out / "sv_concordance.tsv"); evidence["description"]="SV overlaps a candidate feature; parental validation and functional/segregation evidence are still required."
    evidence.to_csv(out / "sv_evidence.tsv", sep="\t", index=False)
    commands = {"minimap2":f"minimap2 -ax asm5 {shlex.quote(str(reference))} {shlex.quote(str(query))} | samtools sort -o assembly_alignment.bam",
                "index":"samtools index assembly_alignment.bam",
                "syri":f"syri -c assembly_alignment.bam -r {shlex.quote(str(reference))} -q {shlex.quote(str(query))} -F B --prefix assembly_"}
    (out / "run_assembly_sv.sh").write_text("#!/usr/bin/env bash\nset -euo pipefail\n"+"\n".join(commands.values())+"\n", encoding="utf-8")
    manifest = {"schema_version":"1.0","reference":file_fingerprint(reference, False),"query":file_fingerprint(query, False),"qtl_interval":qtl_interval,
                "commands":commands,"breakpoint_tolerance":breakpoint_tolerance,"provider_interface":{"required":["chrom","start","end","svtype"],"optional":["call_id","support","mapq","sample","analysis_chain","flags"]},
                "trigger":trigger,"call_count":len(calls),"concordance_count":len(concordance),"independence_rule":"Multiple callers from the same sample/alignment chain are robustness checks, not independent biological evidence.",
                "ranking_rule":"QTL overlap or caller agreement alone emits context only and cannot promote a candidate to Top.",
                "small_variant_boundary":"SNV and short indel consequences use variant-evidence; rearrangements/CNV/PAV and large or split events use sv-plan."}
    manifest_path=out/"sv_manifest.json"; manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    return {"manifest":str(manifest_path),"canonical_sv":str(out/"canonical_sv.tsv"),"concordance":str(out/"sv_concordance.tsv"),"evidence":str(out/"sv_evidence.tsv"),"call_count":len(calls)}
