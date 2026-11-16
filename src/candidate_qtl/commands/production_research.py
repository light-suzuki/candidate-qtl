"""Production-owned Route A command registrations (#669).

Only argparse and shared constants are imported at parser construction time;
the research implementations are loaded when a command is actually run.
"""
from __future__ import annotations
from argparse import _SubParsersAction
from typing import Any
from ._shared import add_json_flag

def register(p: _SubParsersAction) -> None:
    specs = {
      "phenotype-investigate": [("--out-dir", {"required": True}), ("--case-id", {"required": True}), ("--focal-species", {"required": True}), ("--phenotype", {"required": True}), ("--direction", {"default":"unknown"}), ("--inheritance", {"default":"unknown"}), ("--tissue", {"default":""}), ("--stage", {"default":""}), ("--onset", {"default":""}), ("--penetrance", {"default":""}), ("--synonym", {"action":"append"}), ("--ontology-term", {"action":"append"}), ("--exclusion", {"action":"append"}), ("--constraint", {"action":"append"}), ("--uncertainty", {"default":""}), ("--candidates", {}), ("--benchmark-mode", {"action":"store_true"})],
      "phenotype-localization-packet": [("--investigation", {"required": True}), ("--out-dir", {"required": True}), ("--candidates", {})],
      "phenotype-localization-freeze": [("--payload", {"required": True}), ("--investigation", {"required": True}), ("--out-dir", {"required": True})],
      "mechanism-import": [("--payload", {"required": True}), ("--out-dir", {}), ("--preflight", {"action":"store_true", "help":"Validate the original payload and paper reports without creating import artifacts."}), ("--investigation", {}), ("--payload-source", {"choices":("llm", "agent", "human"), "default":"llm"}), ("--model", {"default":""}), ("--case-id", {"default":""}), ("--retrieval-provenance", {"default":""}), ("--candidates", {}), ("--min-mechanisms", {"type":int,"default":5}), ("--proposal-round", {"type":int,"default":1}), ("--interval-informed", {"action":"store_true"}), ("--mechanism-freeze-manifest", {}), ("--localization-freeze-manifest", {}), ("--benchmark-mode", {"action":"store_true"})],
      "mechanism-coverage-audit": [("--mechanisms", {"action":"append", "required":True}), ("--mechanism-import-manifest", {"action":"append", "required":True}), ("--localization-freeze-manifest", {"required":True}), ("--investigation", {"required":True}), ("--out-dir", {"required":True}), ("--candidates", {}), ("--gap-request-audit", {})],
      "investigation-attribution": [("--out-dir", {"required":True}), ("--mechanisms", {}), ("--validation", {}), ("--panel", {}), ("--interval-match", {}), ("--ranked", {}), ("--case-id", {"default":""}), ("--min-mechanisms", {"type":int, "default":5})],
      "reference-gene-resolve": [("--mechanisms", {}), ("--out", {"required": True}), ("--allow-external-api", {"action":"store_true"}), ("--records", {}), ("--timeout", {"type":int,"default":30}), ("--focal-species", {"default":""}), ("--identity-review", {"help":"Identity review ledger TSV (reference_gene_name, reference_species, accession, record_gene_name, verdict, evidence_source_ids, reviewer, rationale). A same_gene row with a citation confirms a provisional declared accession whose record carries record_gene_name."}), ("--sequence-requests", {"help":"Post-freeze public INSDC accession requests; separate from the blind mechanism resolver."}), ("--post-freeze", {"action":"store_true"}), ("--proteome", {}), ("--protein-gene-map", {}), ("--runner", {"choices":("auto", "native", "wsl"), "default":"auto"}), ("--threads", {"type":int, "default":4})],
      "reference-gene-cache": [("--mechanisms", {"required": True}), ("--out", {"required": True}), ("--allow-external-api", {"action":"store_true"}), ("--timeout", {"type":int,"default":60}), ("--max-queries", {"type":int,"default":32}), ("--max-records", {"type":int,"default":128}), ("--max-elapsed-seconds", {"type":int,"default":1800}), ("--resume-manifest", {})],
      "source-verify": [("--mechanisms", {"required": True}), ("--out-dir", {"required": True}), ("--cache", {}), ("--allow-external-api", {"action":"store_true"}), ("--id-column", {"default":"hypothesis_id"}), ("--source-column", {"default":"source_ids"}), ("--verification-policy", {"choices":("initial_bibliography_v1", "content_v1"), "default":"initial_bibliography_v1", "help":"Default: ID/title metadata checks for the retained evidence papers. Explicit content_v1 preserves strict content verification."}), ("--focal-species", {"default":""}), ("--claim-extraction", {}), ("--full-hash", {"action":"store_true"}), ("--reference-bridge", {"help":"Hash-verified sequence correspondence for post-freeze parent review; never rewrites paper facts or promotes support."})],
    }
    for name, args in specs.items():
        sub = p.add_parser(name)
        for flag, kw in args: sub.add_argument(flag, **kw)
        add_json_flag(sub)

def _run(name: str, args: Any) -> dict:
    from .. import phenotype_investigation as pi
    if name == "phenotype-investigate":
        return pi.phenotype_investigate(args.out_dir,args.case_id,args.focal_species,args.phenotype,args.direction,args.inheritance,args.tissue,args.stage,args.onset,args.penetrance,args.synonym,args.ontology_term,args.exclusion,args.constraint,args.uncertainty,args.candidates, benchmark_mode=args.benchmark_mode)
    if name == "phenotype-localization-packet":
        from ..route_a_coverage import phenotype_localization_packet; return phenotype_localization_packet(args.investigation,args.out_dir,args.candidates)
    if name == "phenotype-localization-freeze":
        from ..route_a_coverage import phenotype_localization_freeze; return phenotype_localization_freeze(args.payload,args.investigation,args.out_dir)
    if name == "mechanism-import":
        if args.preflight:
            from ..mechanism_payload_adapter import preflight_payload_file
            return preflight_payload_file(args.payload, args.mechanism_freeze_manifest, args.case_id)
        if not args.out_dir:
            raise ValueError("mechanism-import requires --out-dir unless --preflight is supplied")
        return pi.mechanism_import(payload_arg=args.payload,out_dir_arg=args.out_dir,investigation_arg=args.investigation,payload_source=args.payload_source,model=args.model,case_id=args.case_id,retrieval_provenance=args.retrieval_provenance,candidates_arg=args.candidates,min_mechanisms=args.min_mechanisms,proposal_round=args.proposal_round,interval_informed=args.interval_informed,mechanism_freeze_manifest=args.mechanism_freeze_manifest,localization_freeze_manifest=args.localization_freeze_manifest, benchmark_mode=args.benchmark_mode)
    if name == "mechanism-coverage-audit":
        from ..route_a_coverage import mechanism_coverage_audit
        return mechanism_coverage_audit(args.mechanisms,args.localization_freeze_manifest,args.investigation,args.out_dir,args.candidates,args.gap_request_audit,args.mechanism_import_manifest)
    if name == "investigation-attribution":
        return pi.investigation_attribution(args.out_dir,args.mechanisms,args.validation,args.panel,args.interval_match,args.ranked,args.case_id,args.min_mechanisms)
    if name == "reference-gene-resolve":
        return run_reference_gene_resolve(args)
    if name == "reference-gene-cache":
        from ..reference_record_cache import build_reference_record_cache
        return build_reference_record_cache(args.mechanisms,args.out,allow_external_api=args.allow_external_api,timeout=args.timeout,max_queries=args.max_queries,max_records=args.max_records,max_elapsed_seconds=args.max_elapsed_seconds,resume_manifest_arg=args.resume_manifest)
    from ..source_verify import source_verify; return source_verify(args.mechanisms,args.out_dir,args.cache,args.allow_external_api,args.id_column,args.source_column,args.full_hash,args.claim_extraction,verification_policy=args.verification_policy,focal_species=args.focal_species,reference_bridge_arg=args.reference_bridge)


def run_reference_gene_resolve(args: Any) -> dict:
    if getattr(args, "sequence_requests", None):
        if args.mechanisms or not args.post_freeze or not args.proteome or not args.protein_gene_map:
            raise ValueError("sequence mode requires --post-freeze, --proteome and --protein-gene-map; do not supply --mechanisms")
        from ..reference_sequence_bridge import build_reference_sequence_bridge
        return build_reference_sequence_bridge(args.sequence_requests, args.proteome, args.protein_gene_map, args.out,
            records_arg=args.records, allow_external_api=args.allow_external_api, timeout=args.timeout,
            threads=args.threads, runner=args.runner)
    if not args.mechanisms or any(getattr(args, key, None) for key in ("post_freeze", "proteome", "protein_gene_map")):
        raise ValueError("blind reference resolution requires --mechanisms; sequence inputs require --sequence-requests --post-freeze")
    from ..reference_gene_resolution import reference_gene_resolve
    return reference_gene_resolve(args.mechanisms, args.out, args.allow_external_api, args.records, args.timeout, args.focal_species,
                                  identity_review_arg=getattr(args, "identity_review", None))

HANDLERS = {name: (lambda args, n=name: _run(n,args)) for name in ("phenotype-investigate","phenotype-localization-packet","phenotype-localization-freeze","mechanism-import","mechanism-coverage-audit","investigation-attribution","reference-gene-resolve","reference-gene-cache","source-verify")}
