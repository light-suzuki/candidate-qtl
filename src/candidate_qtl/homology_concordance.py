from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from . import homology_qc
from .core import file_fingerprint, runtime_path


def _specs(values: list[str] | None) -> list[tuple[str, Path]]:
    result=[]
    for value in values or []:
        if "=" not in value: raise ValueError("provider input must use PROVIDER=PATH")
        provider,path=value.split("=",1); result.append((provider,runtime_path(path)))
    return result


def _hits(provider: str, path: Path, top_n: int) -> pd.DataFrame:
    raw=pd.read_csv(path,sep="\t",header=None,comment="#",dtype=str).fillna("")
    if raw.shape[1] < 12: raise ValueError(f"{path}: expected 12+ BLAST/DIAMOND/MMseqs columns")
    raw=raw.iloc[:,:13]; raw.columns=["query_id","subject_id","identity","alignment_length","query_length","subject_length","query_start","query_end","subject_start","subject_end","evalue","bitscore"] + (["description"] if raw.shape[1]==13 else [])
    if "description" not in raw.columns: raw["description"]=""
    for column in ("identity","alignment_length","query_length","subject_length","evalue","bitscore"): raw[column]=pd.to_numeric(raw[column],errors="raise")
    raw["query_coverage"]=raw["alignment_length"]/raw["query_length"].clip(lower=1); raw["subject_coverage"]=raw["alignment_length"]/raw["subject_length"].clip(lower=1)
    raw=raw.sort_values(["query_id","bitscore","evalue","subject_id"],ascending=[True,False,True,True],kind="stable")
    raw["rank"]=raw.groupby("query_id").cumcount()+1; raw=raw[raw["rank"]<=top_n].copy(); raw["provider"]=provider
    raw["quality_class"]=[homology_qc.quality_class(r.identity,r.query_coverage,r.subject_coverage) for r in raw.itertuples()]
    raw["quality_pass"]=[homology_qc.passes_quality(r.identity,r.query_coverage,r.evalue) for r in raw.itertuples()]
    return raw


def _family(subject: str) -> str:
    return subject.split("|")[0].split(".")[0].casefold()


def homology_concordance(hit_specs: list[str], out_dir_arg: str, top_n: int=5,
                         domain_specs: list[str] | None=None, orthology_specs: list[str] | None=None,
                         candidates_arg: str | None=None, metadata_specs: list[str] | None=None) -> dict[str,object]:
    if top_n<1: raise ValueError("top_n must be at least 1")
    out=runtime_path(out_dir_arg); out.mkdir(parents=True,exist_ok=True)
    hit_inputs=_specs(hit_specs); frames=[_hits(provider,path,top_n) for provider,path in hit_inputs]
    hits=pd.concat(frames,ignore_index=True) if frames else pd.DataFrame()
    hits.to_csv(out/"canonical_homology_hits.tsv",sep="\t",index=False)
    providers=sorted(hits["provider"].unique()) if not hits.empty else []
    concordance_rows=[]
    for query,group in hits.groupby("query_id",sort=True) if not hits.empty else []:
        top=group[group["rank"].eq(1)].set_index("provider")
        subjects={provider:str(top.loc[provider,"subject_id"]) for provider in top.index}
        qualities={provider:str(top.loc[provider,"quality_class"]) for provider in top.index}
        if len(subjects)<len(providers): status="provider_missing"
        elif len(set(subjects.values()))==1: status="exact_hit"
        elif len({_family(value) for value in subjects.values()})==1: status="same_family"
        elif all(value=="conserved_domain_only" for value in qualities.values()): status="domain_only"
        else: status="discordant"
        paralog=[]
        for provider,pgroup in group.groupby("provider"):
            ordered=pgroup.sort_values("rank")
            if len(ordered)>1 and homology_qc.is_paralog_ambiguous(ordered.iloc[0]["bitscore"],ordered.iloc[1]["bitscore"]): paralog.append(provider)
        concordance_rows.append({"gene_id":query,"providers_present":";".join(sorted(subjects)),"top_subjects":";".join(f"{p}:{subjects[p]}" for p in sorted(subjects)),
                                 "top_quality_classes":";".join(f"{p}:{qualities[p]}" for p in sorted(qualities)),"concordance":status,
                                 "paralog_ambiguous_providers":";".join(sorted(paralog)),"review_required":status in {"provider_missing","discordant"} or bool(paralog),
                                 "possible_reason":"database_or_isoform_difference" if status=="discordant" else ("partial_or_domain_only_sequence" if status=="domain_only" else "")})
    concordance=pd.DataFrame(concordance_rows,columns=["gene_id","providers_present","top_subjects","top_quality_classes","concordance","paralog_ambiguous_providers","review_required","possible_reason"])
    concordance.to_csv(out/"homology_concordance.tsv",sep="\t",index=False)

    domain_rows=[]
    for provider,path in _specs(domain_specs):
        frame=pd.read_csv(path,sep="\t",dtype=str).fillna(""); required={"gene_id","domain_id","start","end"}
        if not required.issubset(frame.columns): raise ValueError(f"{path}: domain table missing {', '.join(sorted(required-set(frame.columns)))}")
        frame["provider"]=provider; domain_rows.append(frame)
    domains=pd.concat(domain_rows,ignore_index=True,sort=False) if domain_rows else pd.DataFrame(columns=["gene_id","domain_id","start","end","provider"])
    domain_concordance=[]
    for gene,group in domains.groupby("gene_id",sort=True):
        architectures={provider:";".join(f"{r.domain_id}:{r.start}-{r.end}" for r in p.sort_values(["start","end"]).itertuples()) for provider,p in group.groupby("provider")}
        sets={provider:set(p["domain_id"]) for provider,p in group.groupby("provider")}
        status="exact_architecture" if len(set(architectures.values()))==1 else ("same_domain_set" if len({tuple(sorted(v)) for v in sets.values()})==1 else "discordant")
        domain_concordance.append({"gene_id":gene,"architectures":";".join(f"{p}:{architectures[p]}" for p in sorted(architectures)),"concordance":status,"trait_linked":False,"review_required":status=="discordant"})
    pd.DataFrame(domain_concordance,columns=["gene_id","architectures","concordance","trait_linked","review_required"]).to_csv(out/"domain_concordance.tsv",sep="\t",index=False)

    orthology_rows=[]
    for provider,path in _specs(orthology_specs):
        frame=pd.read_csv(path,sep="\t",dtype=str).fillna(""); required={"gene_id","orthogroup","relationship"}
        if not required.issubset(frame.columns): raise ValueError(f"{path}: orthology table missing {', '.join(sorted(required-set(frame.columns)))}")
        frame["provider"]=provider; frame["source_file"]=str(path); orthology_rows.append(frame)
    orthology=pd.concat(orthology_rows,ignore_index=True,sort=False) if orthology_rows else pd.DataFrame(columns=["gene_id","orthogroup","relationship","provider","source_file"])
    orth_concordance=[]
    for gene,group in orthology.groupby("gene_id",sort=True):
        groups={provider:set(p["orthogroup"])-{""} for provider,p in group.groupby("provider")}; relationships=sorted(set(group["relationship"])-{""})
        status="exact_orthogroup" if len({tuple(sorted(v)) for v in groups.values()})==1 else "orthogroup_discordant"
        ambiguous=any(value in {"one_to_many","many_to_many"} for value in relationships)
        orth_concordance.append({"gene_id":gene,"provider_orthogroups":";".join(f"{p}:{','.join(sorted(groups[p]))}" for p in sorted(groups)),"relationships":";".join(relationships),
                                 "concordance":status,"paralog_ambiguous":ambiguous,"review_required":status!="exact_orthogroup" or ambiguous})
    orth_frame=pd.DataFrame(orth_concordance,columns=["gene_id","provider_orthogroups","relationships","concordance","paralog_ambiguous","review_required"])
    orth_frame.to_csv(out/"orthology_concordance.tsv",sep="\t",index=False)
    review=pd.concat([concordance[concordance["review_required"].eq(True)].assign(review_area="homology"),orth_frame[orth_frame["review_required"].eq(True)].assign(review_area="orthology")],ignore_index=True,sort=False)
    review.to_csv(out/"review_required.tsv",sep="\t",index=False)
    evidence=pd.DataFrame([{"gene_id":row.gene_id,"evidence_code":"homology_context","direction":"context","confidence":"low","value":row.concordance,
                            "source_type":"homology_concordance","source_id":str(out/"homology_concordance.tsv"),"description":"Multi-provider homology concordance is computational robustness, not independent trait evidence.",
                            "independence_group":"comparative_annotation_chain"} for row in concordance.itertuples()])
    evidence.to_csv(out/"homology_evidence.tsv",sep="\t",index=False)
    escalation=[]
    if candidates_arg:
        candidates=pd.read_csv(runtime_path(candidates_arg),sep="\t",dtype=str).fillna("")
        selected=candidates[candidates.get("dry_group",pd.Series("",index=candidates.index)).isin({"top","reserve","rescue"})]
        escalation=[{"gene_id":str(row.gene_id),"next_tools":"InterProScan;RBH;OMAmer;OrthoFinder;MAFFT+IQ-TREE_if_ambiguous"} for row in selected.itertuples()]
    pd.DataFrame(escalation,columns=["gene_id","next_tools"]).to_csv(out/"homology_escalation_plan.tsv",sep="\t",index=False)
    metadata={}
    for spec in metadata_specs or []:
        if "=" not in spec: raise ValueError("--provider-metadata must use PROVIDER=JSON")
        provider,value=spec.split("=",1); metadata[provider]=json.loads(runtime_path(value).read_text(encoding="utf-8"))
    manifest={"schema_version":"1.0","top_n":top_n,"hit_inputs":{provider:file_fingerprint(path,False) for provider,path in hit_inputs},"provider_metadata":metadata,
              "candidate_count":int(hits["query_id"].nunique()) if not hits.empty else 0,"review_count":len(review),
              "independence_rule":"BLAST, DIAMOND, MMseqs2, domain, and orthology outputs derived from the same sequences are one comparative-analysis chain, not independent biological evidence.",
              "escalation_rule":"Heavy domain/phylogeny tools are planned only for Top/Reserve/Rescue candidates; generic domains remain context unless a source-backed trait link is supplied."}
    manifest_path=out/"homology_concordance_manifest.json"; manifest_path.write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    return {"manifest":str(manifest_path),"homology_concordance":str(out/"homology_concordance.tsv"),"orthology_concordance":str(out/"orthology_concordance.tsv"),"review_count":len(review)}
