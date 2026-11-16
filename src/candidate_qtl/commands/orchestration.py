"""CLI commands for the orchestration feature group (issue #144)."""

from __future__ import annotations

import json
from pathlib import Path

from ..canonical_family_semantic_wet import run_canonical_route_b_family_semantic_wet_extension
from ..family_index_resolution import resolve_family_index
from ..copy_evidence_inputs import COPY_EVIDENCE_AXES, inspect_copy_evidence, check_copy_evidence_binding
from ..orchestration import orchestrate_plan
from ..guided_profiles import PROFILE_NAMES, guided_profile
from ..phenotype_investigation import PAYLOAD_SOURCES
from ..route_c_direct_pipeline import direct_rescue_project
from ..route_c_direct_members import build_route_c_direct_members
from ..workflow_extensions import case_run
from ..second_stage_recheck import attach_guided_recheck
from ..core import sha256
from argparse import _SubParsersAction
from typing import Any
from ._shared import add_json_flag


def register(subparsers: _SubParsersAction) -> None:
    sub = subparsers.add_parser("guided-profile", help="Start from a screening, full-case, answerless audit, or triple-route profile without hiding the stable commands or review gates.")
    sub.add_argument("--profile", required=True, choices=PROFILE_NAMES)
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--execute", action="store_true", help="Execute runnable stages. Default is a dry-run manifest.")
    sub.add_argument("--fasta"); sub.add_argument("--gff3"); sub.add_argument("--interval", action="append"); sub.add_argument("--phenotype")
    sub.add_argument("--phenotype-context", help="Optional phenotype_context.yaml (case_id, developmental_stage, tissue, environment_or_treatment, ontology term ...). Recorded for audit in the Route A/C manifests; the free text is never replaced or inferred, and no structured field is inserted into any query.")
    sub.add_argument("--allow-external-api", action="store_true", help="Permit profile-specific external providers: NCBI BLAST for screening, and UniProt/Europe PMC verification or retrieval for audit-predict/triple-route. Review the profile manifest's external_transmission fields before execution.")
    sub.add_argument(
        "--allow-reverse-literature-api",
        action="store_true",
        help=("Triple-route only: separately opt in to external Europe PMC lookups for the legacy "
              "diagnostic Route B1 reverse-literature stage. Canonical Route B does not require this."),
    )
    sub.add_argument("--protein-db", help="Screening only: preferred local protein BLAST database (same meaning as `discover --protein-db`).")
    sub.add_argument("--no-auto-local-db", dest="auto_local_db", action="store_false", help="Screening only: disable discovery of the configured local protein DB.")
    sub.set_defaults(auto_local_db=True)
    sub.add_argument("--screening-dir"); sub.add_argument("--case"); sub.add_argument("--resources"); sub.add_argument("--case-id")
    sub.add_argument("--query-manifest"); sub.add_argument("--query-fasta"); sub.add_argument("--blind-manifest"); sub.add_argument("--blind-route")
    sub.add_argument("--predict-manifest", help="audit-benchmark / audit-predict: answerless prediction manifest (case_id/phenotype_text/interval/interval_genes). An answer/truth column is rejected.")
    sub.add_argument("--focal-species", help="focal species for the phenotype-investigation and verification stages (audit-predict / triple-route).")
    sub.add_argument("--direction", default="unknown"); sub.add_argument("--inheritance", default="unknown")
    sub.add_argument("--mechanism-payload", help="audit-predict / triple-route: LLM or hand-authored mechanism payload for mechanism-import (canonical 29/24-column TSV, canonical-subset TSV, JSON/JSONL, or the frozen legacy 7-column TSV, which is converted with an explicit adapter).")
    sub.add_argument("--phenotype-localization-payload", help="audit-predict / triple-route: phenotype-only localization JSON, created before the mechanism payload and frozen by the profile.")
    sub.add_argument("--gap-fill-mechanisms", help="Optional already-imported blind round-2 mechanism_hypotheses.tsv responding to a prior coverage request.")
    sub.add_argument("--gap-fill-mechanism-import-manifest", help="Required with --gap-fill-mechanisms; matching blind round-2 mechanism_import_manifest.json.")
    sub.add_argument("--gap-request-audit", help="Required with --gap-fill-mechanisms; prior mechanism_coverage_audit.json that requested the exact labels.")
    sub.add_argument("--mechanism-freeze-manifest", help="TSV with case_id/mechanism_payload_sha256 (coordination/payload_freeze.tsv format); the payload bytes must match the frozen hash or the import fails closed.")
    sub.add_argument("--payload-source", choices=list(PAYLOAD_SOURCES), default="llm",
                     help="Who produced the mechanism payload: an API completion (llm), this agent driving the flow with its own literature tools (agent), or a person (human). All three pass the same gate.")
    sub.add_argument("--model", default="", help="Model identifier recorded with an LLM payload.")
    sub.add_argument("--record-cache", "--records", dest="record_cache", help="audit-predict / triple-route: local UniProt-shaped record cache for offline verification.")
    sub.add_argument("--source-cache", help="triple-route only: cached Europe PMC records for offline source-verify.")
    sub.add_argument("--bound-claims", help="audit-predict / triple-route: shared source-bound claim TSV required before Route-A bait promotion.")
    sub.add_argument("--agent-claim-review", help="triple-route: pre-reviewed independent-agent packet for the exact frozen source/hypothesis inputs; alternative to --bound-claims.")
    sub.add_argument("--route-b-agent-claim-review", help="Route-B-only independent review packet; never reused by Route A.")
    sub.add_argument("--route-b-source-review-records", help="Route-B-only source snapshot matching the semantic claims.")
    sub.add_argument("--reverse-lit-cache", help="triple-route only: query-keyed cache directory for reverse-literature (offline use).")
    sub.add_argument("--shared-route", help="triple-route only: pre-fixed phenotype-to-pathway route used by the non-scoring B1 overlay.")
    sub.add_argument("--retrieval-index", help="triple-route only: index_manifest.json from retrieval-index-build (Route C Phase 1).")
    sub.add_argument("--proteins", help="triple-route only: candidate protein FASTA for Route B homology and, when --family-index is absent, an interval-scoped fallback family index.")
    sub.add_argument("--reference-function", help="triple-route only: locally cached source-backed reference records for Route B.")
    sub.add_argument("--orthology"); sub.add_argument("--domains"); sub.add_argument("--family-index", help="Preferred reusable genome-wide family_index.tsv. For triple-route, omission triggers lookup in resources.yaml, CANDIDATE_QTL_FAMILY_INDEX, then the per-user cache; Route B may still generate an interval-scoped fallback from --proteins + --protein-db if none is found."); sub.add_argument("--annotations")
    for axis in COPY_EVIDENCE_AXES:
        sub.add_argument("--copy-" + axis.replace("_", "-"),
                         help="Triple-route canonical copy stage only: explicit native family-to-copy "
                              + axis + " evidence file. Omission is not negative evidence.")
    sub.add_argument("--annotation-rescue", help="triple-route only: annotation_rescue.tsv from annotation-rescue, passed to the Route A interval family match stage; a rescue provider that did not complete holds absence rather than falsifying it.")
    sub.add_argument("--benchmark-mode", action="store_true", help="Mark verification as benchmark mode (rejects focal-lineage / answer-backmapped markers).")
    sub.add_argument("--deterministic-subfamily-refinement", action="store_true",
                     help="triple-route --execute only: when the supported member union exceeds the Wet budget, apply one machine-verified curated cross-subfamily refinement pass and re-gate instead of stopping at an adaptive synthesis request.")
    sub.add_argument("--route-c-direct-frozen", help="Triple-route only: pre-interval frozen Direct-C rescue TSV. Search/freeze must have completed before this profile opens the interval.")
    sub.add_argument("--route-c-direct-freeze-manifest", help="Required with --route-c-direct-frozen; binds the interval-blind frozen bytes.")
    sub.add_argument("--route-c-direct-qtl-id", help="Required with --route-c-direct-frozen; QTL identifier recorded during post-freeze projection.")
    sub.add_argument("--route-c-direct-eligible-hits", help="Optional precomputed eligible BLAST hits used only after Direct-C freeze verification.")
    sub.add_argument("--candidates"); sub.add_argument("--evidence"); sub.add_argument("--causal-gene")
    add_json_flag(sub)

    sub = subparsers.add_parser("case-run", help="Run/resume a deterministic manifest-driven case and emit a provenance report.")
    sub.add_argument("--case", required=True, help="Case manifest (YAML). Ranking falls back to top_n=3, reserve_n=3, rescue_n=3, min_top_score=2.0 when the manifest declares none.")
    sub.add_argument("--resources", required=True)
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--mode", choices=["deterministic", "llm-assisted"], default="deterministic")
    sub.add_argument("--resume", action="store_true")
    sub.add_argument("--strict-fingerprints", action="store_true", help="Hash every input by full content and record fingerprint strength.")
    add_json_flag(sub)

    sub = subparsers.add_parser(
        "orchestrate",
        help="Generate a laptop-safe OSS workflow for candidate annotation and primer/marker design.",
    )
    sub.add_argument("--candidates", required=True, help="Dry candidate table or ranked TSV.")
    sub.add_argument("--out-dir", required=True)
    sub.add_argument("--protein-fasta")
    sub.add_argument("--blastp-db", action="append")
    sub.add_argument("--diamond-db", action="append")
    sub.add_argument("--hmmer-db")
    sub.add_argument("--eggnog-data-dir")
    sub.add_argument("--primerblast-root")
    sub.add_argument("--genome")
    sub.add_argument("--gff3")
    sub.add_argument("--nucleotide-db", action="append", help="BLAST nucleotide DB for primerblast-oss; repeatable.")
    sub.add_argument("--candidate-genomic-fasta")
    sub.add_argument("--interval", help="Optional QTL interval for coarse grid-marker design.")
    sub.add_argument("--variants", help="Optional confirmed parental variant TSV.")
    sub.add_argument("--vcf")
    sub.add_argument("--top-n", type=int, default=10)
    sub.add_argument("--num-threads", type=int, default=4)
    add_json_flag(sub)


def _run_case_run(args: Any) -> dict:
    result = case_run(args.case, args.resources, args.out_dir, args.resume, args.mode, args.strict_fingerprints)
    return result


def _run_orchestrate(args: Any) -> dict:
    result = orchestrate_plan(
        args.candidates,
        args.out_dir,
        args.protein_fasta,
        args.blastp_db,
        args.diamond_db,
        args.hmmer_db,
        args.eggnog_data_dir,
        args.primerblast_root,
        args.genome,
        args.gff3,
        args.nucleotide_db,
        args.candidate_genomic_fasta,
        args.interval,
        args.variants,
        args.vcf,
        args.top_n,
        args.num_threads,
    )
    return result


def _record_canonical_execution(result: dict) -> None:
    """Persist canonical progress after the legacy profile, without editing Wet.

    The canonical artifacts own the production state; legacy finality remains
    visible separately and is never a substitute for missing required review.
    """
    canonical = result.get("canonical_route_b_wet")
    if not canonical or not result.get("manifest"):
        return

    def reference(value):
        if not value or not Path(value).is_file():
            return None
        path = Path(value).resolve()
        return {"path": str(path), "sha256": sha256(path)}

    manifest = reference(canonical.get("manifest"))
    canonical_manifest = json.loads(Path(manifest["path"]).read_text(encoding="utf-8")) if manifest else {}
    canonical_state = canonical_manifest.get("state", canonical.get("state"))
    handoff_ref = reference(canonical.get("wet_handoff_manifest"))
    handoff = json.loads(Path(handoff_ref["path"]).read_text(encoding="utf-8")) if handoff_ref else {}
    wet = reference(canonical.get("wet_candidates"))
    second_stage = handoff.get("second_stage_state")
    boundary = result.get("boundary_review")
    remaining = (not (manifest and handoff_ref and wet)
                 or canonical_state == "semantic_review_required"
                 or bool(canonical.get("route_c_members") and not canonical.get("route_b_verified_candidates"))
                 or second_stage in {"requested", "awaiting_review", "review_required_round_budget_exhausted"})
    record = {"manifest": manifest, "state": canonical_state,
              "wet_handoff_manifest": handoff_ref, "wet_candidates": wet,
              "membership_state": handoff.get("membership_state"), "input_coverage_state": handoff.get("input_coverage_state"),
              "second_stage_state": second_stage, "required_work_remaining": bool(remaining),
              "completed_stage": None if remaining else "first_stage",
              "next_action": (boundary["next_action"] if boundary else
                              "continue the recorded canonical requests" if remaining else None),
              "expected_semantic_run_dir": canonical.get("expected_semantic_run_dir"),
              "prompts": reference(canonical.get("prompts")), "bundle": reference(canonical.get("bundle"))}
    path = Path(result["manifest"])
    payload = json.loads(path.read_text(encoding="utf-8"))
    legacy = {key: payload.get(key) for key in ("overall_status", "final", "incomplete_stages")}
    status = "completed" if not remaining else "failed" if legacy["overall_status"] == "failed" else "partial"
    incomplete = ["canonical_route_b_wet"] if remaining else []
    payload.update(canonical_execution=record, legacy_profile=legacy, overall_status=status, final=not remaining,
                   incomplete_stages=incomplete)
    result.update(canonical_execution=record, legacy_profile=legacy, overall_status=status, final=not remaining,
                  incomplete_stages=incomplete)
    temporary = path.with_name(path.name + ".canonical.tmp")
    temporary.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def _record_family_index_resolution(result: dict, resolution: dict[str, object]) -> None:
    """Attach resolver provenance to the current guided-profile manifest."""
    manifest_value = result.get("manifest")
    if not manifest_value:
        return
    manifest_path = Path(str(manifest_value))
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError):
        return
    if not isinstance(manifest, dict):
        return
    manifest["family_index_auto_resolve_attempted"] = bool(resolution.get("attempted"))
    manifest["family_index_resolution"] = resolution
    tmp = manifest_path.with_name(manifest_path.name + ".tmp")
    tmp.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    tmp.replace(manifest_path)


def _run_guided_profile(args: Any) -> dict:
    copy_inputs = {axis: getattr(args, "copy_" + axis, None) for axis in COPY_EVIDENCE_AXES
                   if getattr(args, "copy_" + axis, None) is not None}
    if copy_inputs and args.profile != "triple-route":
        raise ValueError("--copy-* inputs require --profile triple-route; they are not screening evidence")
    if args.profile == "triple-route" and args.execute:
        _, copy_fingerprints = inspect_copy_evidence(copy_inputs)
        check_copy_evidence_binding(Path(args.out_dir), copy_fingerprints)
    family_index = args.family_index
    family_index_resolution = None
    if args.profile == "triple-route":
        family_index, family_index_resolution = resolve_family_index(args.family_index, args.resources)
    direct_values = (
        bool(args.route_c_direct_frozen), bool(args.route_c_direct_freeze_manifest),
        bool(args.route_c_direct_qtl_id),
    )
    if any(direct_values) and not all(direct_values):
        raise ValueError(
            "Direct C guided wiring requires --route-c-direct-frozen, "
            "--route-c-direct-freeze-manifest and --route-c-direct-qtl-id together"
        )
    if any(direct_values) and args.profile != "triple-route":
        raise ValueError("Direct C guided wiring is available only for --profile triple-route")

    result = guided_profile(
        args.profile, args.out_dir, execute=args.execute, fasta=args.fasta, gff3=args.gff3,
        intervals=args.interval, phenotype=args.phenotype, allow_external_api=args.allow_external_api,
        allow_reverse_literature_api=args.allow_reverse_literature_api,
        protein_db=args.protein_db, auto_local_db=args.auto_local_db,
        screening_dir=args.screening_dir, case=args.case, resources=args.resources, case_id=args.case_id,
        query_manifest=args.query_manifest, query_fasta=args.query_fasta, blind_manifest=args.blind_manifest,
        blind_route=args.blind_route, predict_manifest=args.predict_manifest,
        focal_species=args.focal_species, direction=args.direction, inheritance=args.inheritance,
        mechanism_payload=args.mechanism_payload, payload_source=args.payload_source, model=args.model,
        phenotype_localization_payload=args.phenotype_localization_payload,
        gap_fill_mechanisms=args.gap_fill_mechanisms,
        gap_fill_mechanism_import_manifest=args.gap_fill_mechanism_import_manifest,
        gap_request_audit=args.gap_request_audit,
        mechanism_freeze_manifest=args.mechanism_freeze_manifest,
        phenotype_context=args.phenotype_context,
        record_cache=args.record_cache, source_cache=args.source_cache,
        bound_claims=args.bound_claims,
        agent_claim_review=args.agent_claim_review,
        route_b_agent_claim_review=args.route_b_agent_claim_review,
        route_b_source_review_records=args.route_b_source_review_records,
        reverse_lit_cache=args.reverse_lit_cache,
        shared_route=args.shared_route,
        retrieval_index=args.retrieval_index, proteins=args.proteins,
        reference_function=args.reference_function, orthology=args.orthology, domains=args.domains,
        family_index=family_index, annotations=args.annotations, benchmark_mode=args.benchmark_mode,
        candidates=args.candidates, evidence=args.evidence, causal_gene=args.causal_gene,
        annotation_rescue=args.annotation_rescue,
    )
    if family_index_resolution is not None:
        _record_family_index_resolution(result, family_index_resolution)
        result["family_index_auto_resolve_attempted"] = bool(family_index_resolution.get("attempted"))
        result["family_index_resolution"] = family_index_resolution
    if args.profile == "triple-route" and args.execute:
        route_c_members = None
        if all(direct_values):
            atlas = Path(args.out_dir) / "b01_route_b" / "interval_gene_atlas.tsv"
            projection = direct_rescue_project(
                args.route_c_direct_frozen, args.route_c_direct_freeze_manifest,
                str(atlas), args.route_c_direct_qtl_id,
                str(Path(args.out_dir) / "c04_direct_projection"),
                args.route_c_direct_eligible_hits,
            )
            members = build_route_c_direct_members(
                args.route_c_direct_frozen, args.route_c_direct_freeze_manifest,
                str(projection["projection"]), str(projection["manifest"]), str(atlas),
                str(Path(args.out_dir) / "c05_direct_members"), case_id=args.case_id or "",
            )
            route_c_members = str(members["members"])
            result["route_c_direct_projection"] = projection
            result["route_c_direct_members"] = members
        result["canonical_route_b_wet"] = run_canonical_route_b_family_semantic_wet_extension(
            args.out_dir,
            phenotype=args.phenotype or "",
            source_cache_arg=args.source_cache,
            bound_claims_arg=args.bound_claims,
            agent_claim_review_arg=args.route_b_agent_claim_review,
            source_review_records_arg=args.route_b_source_review_records,
            record_cache_arg=args.record_cache,
            allow_external_api=args.allow_external_api,
            deterministic_subfamily_refinement=bool(getattr(args, "deterministic_subfamily_refinement", False)),
            case_id=args.case_id or "",
            route_c_members_arg=route_c_members,
            **({"copy_evidence_args": copy_inputs} if copy_inputs else {}),
        )
    if copy_inputs:
        copy_state = {"state": "planned_not_consumed", "explicit_inputs": copy_inputs}
        if args.execute:
            canonical = result["canonical_route_b_wet"]
            payload = json.loads(Path(canonical["manifest"]).read_text(encoding="utf-8"))
            copy_state = payload["copy_evidence"]
        result["copy_evidence"] = copy_state
        manifest_path = Path(result["manifest"])
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        payload["copy_evidence"] = copy_state
        temporary = manifest_path.with_name(manifest_path.name + ".copy.tmp")
        temporary.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        temporary.replace(manifest_path)
    attach_guided_recheck(result, args)
    _record_canonical_execution(result)
    return result


HANDLERS = {
    "guided-profile": _run_guided_profile,
    "case-run": _run_case_run,
    "orchestrate": _run_orchestrate,
}
