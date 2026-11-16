"""Production-owned Direct-C registrations with lazy implementation imports."""
from __future__ import annotations
from argparse import _SubParsersAction
from typing import Any
from ._shared import add_json_flag

def register(p: _SubParsersAction) -> None:
    specs = {
      "route-c-direct-search": [("--phenotype",{"required":True}),("--focal-species",{"required":True}),("--species-scope",{"required":True}),("--out-dir",{"required":True}),("--cache",{}),("--allow-external-api",{"action":"store_true"}),("--limit",{"type":int,"default":10}),("--case-id",{"default":""}),("--full-hash",{"action":"store_true"})],
      "route-c-direct-freeze": [("--proposals",{"required":True}),("--out-dir",{"required":True}),("--citation-cache",{}),("--allow-external-api",{"action":"store_true"}),("--claim-extraction",{}),("--full-hash",{"action":"store_true"})],
      "route-c-provider-bait": [("--frozen",{"required":True}),("--freeze-manifest",{"required":True}),("--out-dir",{"required":True}),("--cache",{"help":"Directory of cached UniProt literature-search JSON (uniprot_lit_<PMID>.json)."}),("--allow-external-api",{"action":"store_true","help":"Query UniProt by public PubMed ID when the cache lacks it."})],
      "route-c-direct-project": [("--frozen",{"required":True}),("--freeze-manifest",{"required":True}),("--atlas",{"required":True}),("--qtl-id",{"required":True}),("--eligible-hits",{}),("--identity-mapping",{}),("--reference-bridge",{"help":"Hash-verified public sequence bridge; keeps all matched copies review-only."}),("--out-dir",{"required":True}),("--full-hash",{"action":"store_true"})],
      "route-c-direct-members": [("--frozen",{"required":True}),("--freeze-manifest",{"required":True}),("--projection",{"required":True}),("--projection-manifest",{"required":True}),("--atlas",{"required":True}),("--out-dir",{"required":True}),("--case-id",{"default":""})],
      "triple-route-compare": [("--route-a",{"required":True}),("--route-b",{"required":True}),("--route-c",{"required":True}),("--out-dir",{"required":True}),("--case-id",{"default":""}),("--route-c-interval-match",{})],
    }
    for name,args in specs.items():
        sub=p.add_parser(name)
        for flag,kw in args: sub.add_argument(flag,**kw)
        add_json_flag(sub)

def _run(name: str, args: Any) -> dict:
    from .. import route_c_direct_pipeline as pipeline
    if name == "route-c-direct-search": return pipeline.direct_literature_search(args.phenotype,args.focal_species,args.species_scope,args.out_dir,args.cache,args.allow_external_api,args.limit,args.case_id,args.full_hash)
    if name == "route-c-direct-freeze": return pipeline.direct_rescue_freeze(args.proposals,args.out_dir,args.citation_cache,args.allow_external_api,args.claim_extraction,args.full_hash)
    if name == "route-c-direct-project":
        extra = {"reference_bridge_arg": args.reference_bridge} if getattr(args, "reference_bridge", None) else {}
        return pipeline.direct_rescue_project(args.frozen,args.freeze_manifest,args.atlas,args.qtl_id,args.out_dir,args.eligible_hits,args.full_hash,args.identity_mapping,**extra)
    if name == "route-c-provider-bait":
        from ..route_c_provider_bait import build_provider_bait
        return build_provider_bait(args.frozen,args.freeze_manifest,args.out_dir,cache_arg=args.cache,allow_external_api=args.allow_external_api)
    if name == "triple-route-compare":
        from ..retrieval_route import triple_route_compare
        return triple_route_compare(args.route_a, args.route_b, args.route_c, args.out_dir, args.case_id, args.route_c_interval_match)
    from ..route_c_direct_members import build_route_c_direct_members
    return build_route_c_direct_members(
        args.frozen,
        args.freeze_manifest,
        args.projection,
        args.projection_manifest,
        args.atlas,
        args.out_dir,
        case_id=args.case_id,
    )

HANDLERS={name:(lambda args,n=name:_run(n,args)) for name in ("route-c-direct-search","route-c-direct-freeze","route-c-provider-bait","route-c-direct-project","route-c-direct-members","triple-route-compare")}
