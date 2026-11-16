"""Guided entry profiles over the stable command surface (issue #156).

Profiles are deliberately thin orchestration.  They expose the commands and
arguments that will run, retain review gates, and never implement ranking rules.
"""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

from . import __version__
from .blast_evidence import run_panel_blastp
from .blast_panel import freeze_blast_panel, panel_from_verified_bait_panel
from .core import runtime_path
from .discover import discover_candidates, promote_discover
from .dual_route import interval_reverse_lookup, phenotype_first_status
from .family_candidate_bridge import family_candidate_bridge
from .family_to_copy import family_to_copy
from .family_index import route_a_interval_family_match
from .llm_pathway import bait_panel_freeze, hypothesis_verify
from .mechanism_payload_adapter import validate_payload_file
from .phenotype_investigation import investigation_attribution, mechanism_import, phenotype_investigate
from .phenotype_context import structured_context_block
from .query_eligibility import query_eligibility
from .ranking_qa import ranking_qa
from .reference_gene_resolution import reference_gene_resolve
from .route_a_coverage import (
    mechanism_coverage_audit,
    phenotype_localization_freeze,
    phenotype_localization_packet,
)
from .reverse_literature import reverse_literature
from .route_b_prepare import route_b_prepare
from .source_verify import source_verify
from .workflow_extensions import case_run


PROFILE_NAMES = ("screening", "full-case", "audit-benchmark", "audit-predict", "triple-route")
TRUTH_DEPENDENT_LANES = ("robustness_benchmark", "answer_recall", "answer_derived_eligibility_gate")
# audit lanes are independent by contract, so one failed lane must not cancel the
# rest; the case-oriented profiles stay fail-fast but still persist a manifest.
_INDEPENDENT_LANE_PROFILES = ("audit-benchmark", "audit-predict", "triple-route")
# Profiles that must run without any truth input.
_PREDICTION_ONLY_PROFILES = ("audit-benchmark", "audit-predict")
# Profiles that run A/B/C as separate lanes and only meet at the comparison.
_TRIPLE_ROUTE_PROFILES = ("triple-route",)


_CHAIN_PREFIX = "<from "


def _stage(command: str, arguments: dict[str, object], required: list[str], *,
           external_transmission: str = "none", optional_tools: list[str] | None = None,
           stop_condition: str = "command_failure",
           provider_policy: str | None = None,
           inputs_from: dict[str, tuple[str, str]] | None = None,
           optional_chain: bool = False,
           optional_chain_inputs: list[str] | None = None) -> dict[str, object]:
    """Describe one stage.

    ``inputs_from`` names arguments produced by an earlier stage in the same
    profile ({argument: (stage command, result key)}). They are not "missing" at
    plan time -- they are pending -- and they are filled from the earlier stage's
    result at execution, so a chained Route A run does not require the operator
    to copy paths between commands (#215).

    ``optional_chain_inputs`` marks individual chained arguments whose upstream
    may legitimately produce nothing (e.g. a provider that was not resolved);
    those fill as ``None`` instead of blocking the stage, while the other chains
    stay required.
    """
    inputs_from = inputs_from or {}
    for name, (source, key) in inputs_from.items():
        if not arguments.get(name):
            arguments[name] = f"{_CHAIN_PREFIX}{source}.{key}>"
    missing = [name for name in required if not arguments.get(name)]
    stage: dict[str, object] = {
        "command": command,
        "arguments": arguments,
        "candidate_qtl_version": __version__,
        "required_inputs": required,
        "missing_inputs": missing,
        "status": "ready" if not missing else "not_runnable_missing_input",
        "external_transmission": external_transmission,
        "optional_tools": optional_tools or [],
        "stop_condition": stop_condition,
    }
    if inputs_from:
        stage["inputs_from"] = {name: {"stage": source, "result_key": key}
                                for name, (source, key) in inputs_from.items()}
        # An optional chain reports on whatever upstream stages did produce. That
        # is the point of the attribution stage: "nothing ran" is the finding.
        stage["optional_chain"] = optional_chain
        if optional_chain_inputs:
            stage["optional_chain_inputs"] = optional_chain_inputs
    if provider_policy:
        stage["provider_policy"] = provider_policy
    return stage


# Which lane a stage belongs to. Lanes may not read each other's intermediate
# output before the comparison: A and B traverse the same literature, so the only
# independence that can be protected is procedural -- two investigators who do
# not confer until they compare notes.
_LANE_OF = {
    "phenotype-investigate": "a", "phenotype-localization-packet": "a",
    "phenotype-localization-freeze": "a", "mechanism-import": "a",
    "mechanism-coverage-audit": "a", "reference-gene-resolve": "a",
    "source-verify": "a",
    "hypothesis-verify": "a", "bait-panel-freeze": "a", "route-a-interval-match": "a",
    "blast-panel-freeze": "a", "blast-panel-search": "a",
    "phenotype-first-status": "a",
    "route-b-prepare": "b", "interval-reverse-lookup": "b", "reverse-literature": "b",
    "interval-atlas": "b",
    "retrieval-run": "c", "retrieval-import": "c",
}
# Stages that legitimately read more than one lane (or run after the meeting
# point and therefore may read lane outputs).
_MEETING_POINT_COMMANDS = frozenset({"triple-route-compare", "retrieval-interval-match", "investigation-attribution",
                                     "family-to-copy", "family-candidate-bridge"})


def _assert_lane_isolation(stage: dict[str, object], out_dir: Path) -> None:
    """Refuse a stage whose input path lives under another lane's output.

    This is the structural half of the independence rule: the profile lays the
    lanes out in sibling directories, so an argument that reaches across is
    detectable without knowing what any command means.
    """
    command = str(stage["command"])
    lane = _LANE_OF.get(command)
    if lane is None or command in _MEETING_POINT_COMMANDS:
        return
    arguments = stage["arguments"]
    assert isinstance(arguments, dict)
    root = out_dir.resolve()
    for name, value in arguments.items():
        if name == "out_dir" or not isinstance(value, str) or not value:
            continue
        try:
            candidate = Path(value).resolve()
        except (OSError, ValueError):
            continue
        if root not in candidate.parents:
            continue
        relative = candidate.relative_to(root).parts
        if not relative:
            continue
        other = _LANE_OF.get("", "")
        prefix = str(relative[0])[:1]
        if prefix in {"a", "b", "c"} and prefix != lane:
            other = prefix
            raise ValueError(
                f"lane isolation violated: stage {command} (lane {lane}) was given {name}={value}, "
                f"which is lane {other}'s output. The lanes must not read each other before the comparison"
            )


def _resolve_chained(stage: dict[str, object], results: list[dict[str, object]]) -> str:
    """Fill pending arguments from completed upstream stages.

    Returns an empty string when the stage is runnable, or the status to record
    when an upstream stage did not produce what this one needs. An unavailable
    upstream is an unavailable stage, not a failed one.
    """
    chained = stage.get("inputs_from") or {}
    assert isinstance(chained, dict)
    optional = bool(stage.get("optional_chain"))
    optional_names = set(stage.get("optional_chain_inputs") or [])
    arguments = stage["arguments"]
    assert isinstance(arguments, dict)
    for name, spec in chained.items():
        if not str(arguments.get(name, "")).startswith(_CHAIN_PREFIX):
            continue  # explicitly supplied by the operator
        source, key = str(spec["stage"]), str(spec["result_key"])
        record = next((item for item in results
                       if item.get("command") == source and item.get("status") == "completed"), None)
        value = (record.get("result") or {}).get(key) if record else None
        if not value:
            if not (optional or name in optional_names):
                return (f"not_runnable_upstream_unavailable:{source}" if record is None
                        else f"not_runnable_upstream_output_missing:{source}.{key}")
            arguments[name] = None
            continue
        arguments[name] = value
    return ""


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _plan_time_rank_lane_preflight(values: dict[str, object]) -> dict[str, object] | None:
    """Preflight the rank lane of a planned case-run stage (issue #346).

    At plan time only the case manifest is known, not the evidence rows. When
    ``inputs.evidence_resources`` is absent or empty the rank lane is already
    known to be non-informative: every candidate scores 0.0 and screening_rank
    is a full tie that degenerates to coordinate order. When evidence is
    declared, row-level confirmation happens at execution, and the executed
    case-run result replaces this preflight's state fields with the
    authoritative values (``plan_time_detection: false``). Returns None when
    the profile has no case manifest to read at plan time.
    """
    case_arg = values.get("case")
    if not case_arg:
        return None
    try:
        case = yaml.safe_load(runtime_path(str(case_arg)).read_text(encoding="utf-8")) or {}
    except (OSError, ValueError):
        return None
    if not isinstance(case, dict):
        return None
    declared = list((case.get("inputs") or {}).get("evidence_resources") or [])
    if declared:
        return {
            "plan_time_detection": True,
            "evidence_resources_declared": True,
            "rank_semantic_status": "evaluable",
            "rank_is_positional_fallback": False,
            "rank_metric_eligible": True,
            "stage_a_is_budget_only": False,
            "warning": "",
        }
    return {
        "plan_time_detection": True,
        "evidence_resources_declared": False,
        "rank_semantic_status": "not_evaluated_no_scoring_evidence",
        "rank_is_positional_fallback": True,
        "rank_metric_eligible": False,
        "stage_a_is_budget_only": True,
        "warning": "rank lane is non-informative: no scoring evidence (positional fallback)",
    }


def _plan_time_payload_schema_preflight(values: dict[str, object]) -> dict[str, object] | None:
    """Fail-closed schema preflight of the mechanism payload at plan time (#380).

    Uses the same validate_payload_file as mechanism-import (single source of
    truth) so a payload that will be rejected at execution is rejected before
    any stage runs. Returns None when the profile has no mechanism payload;
    otherwise the detection record: status accepted carries the accepted schema
    version, data-row count and raw SHA; status rejected carries the fail-closed
    error message(s) with the payload still never executed. freeze_verification
    is present only when a freeze manifest was supplied and matched.
    """
    payload = values.get("mechanism_payload")
    if not payload:
        return None
    try:
        validated = validate_payload_file(
            payload, values.get("mechanism_freeze_manifest"), str(values.get("case_id") or ""))
    except ValueError as error:
        return {
            "plan_time_detection": True,
            "status": "rejected",
            "payload_file": str(payload),
            "accepted_schema_version": None,
            "row_count": None,
            "raw_payload_sha256": None,
            "errors": [str(error)],
        }
    preflight: dict[str, object] = {
        "plan_time_detection": True,
        "status": "accepted",
        "payload_file": str(payload),
        "accepted_schema_version": str(validated["schema"]),
        "row_count": int(validated["rows_read"]),
        "raw_payload_sha256": str(validated["raw_payload_sha256"]),
        "errors": list(validated["errors"]),
    }
    if int(validated["rows_read"]) == 0:
        preflight["row_count_note"] = "zero_rows_schema_accepted"
    freeze = validated.get("freeze_verification")
    if freeze is not None:
        preflight["freeze_verification"] = freeze
    return preflight


def _write_manifest(path: Path, manifest: dict[str, object]) -> None:
    """Persist the manifest with an atomic replace so a mid-run failure or an
    interrupt still leaves a complete, readable audit record (#160)."""
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(tmp, path)


_GUIDED_PROFILE_MANIFEST = "guided_profile_manifest.json"


def _new_run_id() -> str:
    """Unique, sortable identifier for one guided-profile run (issue #302)."""
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")


def _archive_superseded_manifest(out_dir: Path, manifest_path: Path) -> str | None:
    """Preserve a manifest that a new run is about to overwrite, marked superseded.

    Every guided-profile run replaces ``guided_profile_manifest.json`` in place,
    so an unhandled rerun would silently destroy the previous run's record
    (issue #302). When a manifest already exists it is copied to
    ``guided_profile_manifest.<run_id>.json`` -- named after its own run_id (an
    mtime-derived id for a legacy manifest without one) -- and gains
    ``final: false`` and ``superseded_by`` pointing back at the canonical path,
    so the old version is machine-detectable as not final. Returns the archive
    path, or None when no previous manifest existed. An unreadable previous
    manifest is skipped: it is overwritten exactly as before, never an error.
    """
    if not manifest_path.exists():
        return None
    try:
        previous = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if not isinstance(previous, dict):
        return None
    previous_id = str(previous.get("run_id") or f"legacy_{manifest_path.stat().st_mtime_ns}")
    archive_path = out_dir / f"guided_profile_manifest.{previous_id}.json"
    previous["final"] = False
    previous["superseded_by"] = str(manifest_path)
    _write_manifest(archive_path, previous)
    return str(archive_path)


def resolve_final_manifest(out_dir_arg: str) -> dict[str, object]:
    """Machine-determine the final guided-profile manifest of an out-dir (#302).

    Rule: a manifest with ``final: true`` and no ``superseded_by`` is the final
    version; a manifest carrying ``superseded_by`` was replaced; a manifest with
    an explicit ``final: false`` is not final; a legacy manifest without either
    field cannot be judged and is reported as ``unknown`` -- never an error.
    """
    out_dir = runtime_path(out_dir_arg)
    final_candidates: list[str] = []
    superseded: list[str] = []
    not_final: list[str] = []
    unknown: list[str] = []
    for path in sorted(out_dir.glob("guided_profile_manifest*.json")):
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            unknown.append(str(path))
            continue
        if not isinstance(data, dict) or "final" not in data:
            unknown.append(str(path))
            continue
        if data.get("final") is True and not data.get("superseded_by"):
            final_candidates.append(str(path))
        elif data.get("superseded_by"):
            superseded.append(str(path))
        else:
            not_final.append(str(path))
    return {
        "final": final_candidates[0] if len(final_candidates) == 1 else None,
        "final_candidates": final_candidates,
        "superseded": superseded,
        "not_final": not_final,
        "unknown": unknown,
    }


def _plan(profile: str, values: dict[str, object], out_dir: Path) -> tuple[list[dict[str, object]], str]:
    if profile == "screening":
        # Provider selection must match `discover` exactly: a configured local
        # protein DB is auto-discovered unless the operator disables it.
        auto_local_db = bool(values.get("auto_local_db", True))
        stages = [_stage(
            "discover",
            {
                "fasta": values.get("fasta"), "gff3": values.get("gff3"),
                "interval": values.get("intervals") or [], "phenotype": values.get("phenotype"),
                "out_dir": str(out_dir / "screening"), "offline": not bool(values.get("allow_external_api")),
                "allow_external_api": bool(values.get("allow_external_api")),
                "protein_db": values.get("protein_db"), "auto_local_db": auto_local_db,
            },
            ["fasta", "gff3", "interval", "phenotype"],
            external_transmission="candidate protein sequences may be sent to NCBI BLAST" if values.get("allow_external_api") else "none (offline)",
            optional_tools=["local BLASTP database"],
            stop_condition="coordinate gate rejected, extraction failure, or preliminary screening complete",
            provider_policy=(
                "--protein-db explicit local BLASTP database" if values.get("protein_db")
                else ("auto-discover the configured local protein DB (same default as `discover`)" if auto_local_db
                      else "local protein DB discovery disabled via --no-auto-local-db")
            ),
        )]
        promotion = "Run full-case with --screening-dir <out>/screening. It will create reviewable case/resources manifests and stop; REVIEW_REQUIRED fields are never auto-filled."
        return stages, promotion

    if profile == "full-case":
        screening_dir = values.get("screening_dir")
        if screening_dir:
            stages = [_stage(
                "discover-promote",
                {"discover_dir": screening_dir, "out_dir": str(out_dir / "promotion"), "case_id": values.get("case_id") or ""},
                ["discover_dir"],
                stop_condition="reviewable manifests emitted; stop for explicit review of REVIEW_REQUIRED fields",
            )]
            return stages, "Review promotion/case.yaml and promotion/resources.yaml, then rerun full-case with --case and --resources."
        stages = [_stage(
            "case-run",
            {"case": values.get("case"), "resources": values.get("resources"), "out_dir": str(out_dir / "case_run"), "mode": "deterministic"},
            ["case", "resources"],
            optional_tools=["R/qtl2", "BLAST/DIAMOND", "HMMER", "primerblast-oss"],
            stop_condition="source/mapping/review gate, missing required input, or report completion",
        )]
        return stages, "The standard promotion path is screening -> full-case --screening-dir -> human review -> full-case --case/--resources."

    if profile == "triple-route":
        allow_external = bool(values.get("allow_external_api"))
        allow_reverse_literature = bool(values.get("allow_reverse_literature_api"))
        stages = [
            # --- Route A: phenotype-first, interval-blind while generating -----
            _stage("phenotype-investigate",
                   {"case_id": values.get("case_id"), "focal_species": values.get("focal_species"),
                    "phenotype": values.get("phenotype"), "direction": values.get("direction") or "unknown",
                    "inheritance": values.get("inheritance") or "unknown", "candidates": values.get("candidates"),
                    "phenotype_context": values.get("phenotype_context"),
                    "out_dir": str(out_dir / "a01_phenotype_investigation")},
                   ["case_id", "focal_species", "phenotype"]),
            _stage("phenotype-localization-packet",
                   {"candidates": values.get("candidates"),
                    "out_dir": str(out_dir / "a01b_phenotype_localization_packet")},
                   [],
                   inputs_from={"investigation": ("phenotype-investigate", "investigation")},
                   stop_condition="phenotype-only localization prompt emitted before mechanism generation"),
            _stage("phenotype-localization-freeze",
                   {"payload": values.get("phenotype_localization_payload"),
                    "out_dir": str(out_dir / "a01c_phenotype_localization_freeze")},
                   ["payload"],
                   inputs_from={"investigation": ("phenotype-investigate", "investigation")},
                   stop_condition="localization response validated and hash-frozen before mechanism generation"),
            _stage("mechanism-import",
                   {"payload": values.get("mechanism_payload"), "payload_source": values.get("payload_source") or "llm",
                    "model": values.get("model") or "", "case_id": values.get("case_id") or "",
                    "mechanism_freeze_manifest": values.get("mechanism_freeze_manifest"),
                    "candidates": values.get("candidates"), "out_dir": str(out_dir / "a02_mechanism_import")},
                   ["payload"],
                   inputs_from={"investigation": ("phenotype-investigate", "investigation"),
                                "localization_freeze_manifest": ("phenotype-localization-freeze", "manifest")},
                   stop_condition="well-formed initial labeled hypotheses imported after localization freeze (target 12; count remains advisory)"),
            _stage("mechanism-coverage-audit",
                   {"gap_fill_mechanisms": values.get("gap_fill_mechanisms"),
                    "gap_fill_import_manifest": values.get("gap_fill_mechanism_import_manifest"),
                    "gap_request_audit": values.get("gap_request_audit"),
                    "candidates": values.get("candidates"),
                    "out_dir": str(out_dir / "a02c_mechanism_coverage")},
                   (["gap_fill_mechanisms", "gap_fill_import_manifest", "gap_request_audit"]
                    if values.get("gap_fill_mechanisms") else []),
                   inputs_from={"initial_mechanisms": ("mechanism-import", "mechanisms"),
                                "initial_import_manifest": ("mechanism-import", "manifest"),
                                "localization_freeze_manifest": ("phenotype-localization-freeze", "manifest"),
                                "investigation": ("phenotype-investigate", "investigation")},
                   stop_condition="valid hypotheses continue with explicit degraded coverage; a blind gap prompt is advisory"),
            _stage("reference-gene-resolve",
                   {"allow_external_api": allow_external, "record_cache": values.get("record_cache"),
                    "focal_species": values.get("focal_species") or "",
                    "out_dir": str(out_dir / "a02b_reference_gene_resolution")},
                   [],
                   inputs_from={"mechanisms": ("mechanism-coverage-audit", "mechanisms")},
                   external_transmission=("UniProt REST gene/accession resolution"
                                          if allow_external else "none (offline record cache only)"),
                   stop_condition=("every missing reference gene citation is provider-resolved or left "
                                   "review-required; malformed accessions are never queried")),
            _stage("source-verify",
                   {"cache": values.get("source_cache"), "allow_external_api": allow_external,
                    "out_dir": str(out_dir / "a03_source_verify")},
                   [],
                   inputs_from={"mechanisms": ("reference-gene-resolve", "mechanisms_resolved"),
                                "claim_extraction": ("mechanism-import", "claim_extractions")},
                   optional_chain_inputs=["claim_extraction"],
                   external_transmission="citation identifiers only (Europe PMC)" if allow_external else "none (offline cache only)",
                   stop_condition="every citation resolved, rejected, or left review-required"),
            _stage("hypothesis-verify",
                   {"focal_species": values.get("focal_species"), "record_cache": values.get("record_cache"),
                    "bound_claims": values.get("bound_claims"),
                    "agent_claim_review": values.get("agent_claim_review"),
                    "allow_external_api": allow_external, "benchmark_mode": bool(values.get("benchmark_mode")),
                    "out_dir": str(out_dir / "a04_hypothesis_verify")},
                   ["focal_species"],
                   inputs_from={"hypotheses": ("source-verify", "mechanisms_with_verified_sources"),
                                "source_records": ("source-verify", "verified_sources"),
                                "source_review_records": ("source-verify", "source_review_records")},
                   external_transmission="UniProt REST accession/species/sequence verification" if allow_external else "none (offline record cache only)"),
            _stage("bait-panel-freeze",
                   {"case_id": values.get("case_id") or "", "out_dir": str(out_dir / "a05_bait_panel_freeze")},
                   [],
                   inputs_from={"hypotheses": ("source-verify", "mechanisms_with_verified_sources"),
                                "validation": ("hypothesis-verify", "validation"),
                                "sequences": ("hypothesis-verify", "verified_sequences"),
                                "reviewable_sequences": ("hypothesis-verify", "reviewable_sequences")}),
            # Route A interval family match: only after the verified bait panel
            # is frozen. Forward-matches verified families to concrete genes
            # inside the interval via the reference-genome family index; the
            # interval is never opened for hypothesis generation.
            _stage("route-a-interval-match",
                   {"case_id": values.get("case_id") or "", "out_dir": str(out_dir / "a05b_interval_family_match"),
                    "family_index": values.get("family_index"), "interval": values.get("intervals") or [],
                    "candidates": values.get("candidates"), "annotation_rescue": values.get("annotation_rescue")},
                   [],
                   inputs_from={"panel": ("bait-panel-freeze", "panel"),
                                "hypotheses": ("mechanism-coverage-audit", "hypotheses")},
                   optional_chain=True,
                   optional_chain_inputs=["panel", "hypotheses"],
                   stop_condition=("verified targets only are falsified against the reference-genome family index; "
                                   "zero verified baits or a missing index degrade to not_evaluated, never to a miss")),
            _stage("phenotype-first-status",
                   {"case_id": values.get("case_id") or "", "out_dir": str(out_dir / "a06_route_a_status"),
                    "phenotype_context": values.get("phenotype_context"),
                    "phenotype": values.get("phenotype") or ""},
                   [],
                   inputs_from={"hypotheses": ("mechanism-coverage-audit", "hypotheses"),
                                "validation": ("hypothesis-verify", "validation"),
                                "panel": ("bait-panel-freeze", "panel"),
                                "falsification": ("route-a-interval-match", "hypothesis_interval_falsification")},
                   optional_chain=True,
                   optional_chain_inputs=["falsification"],
                   stop_condition="Route A hypothesis table emitted for the comparison"),
            # --- Route C Phase 1: phenotype-only, before any interval artifact --
            _stage("retrieval-run",
                   {"index": values.get("retrieval_index"), "phenotype": values.get("phenotype"),
                    "case_id": values.get("case_id") or "", "out_dir": str(out_dir / "c01_retrieval_run"),
                    "phenotype_context": values.get("phenotype_context")},
                   ["index", "phenotype"],
                   stop_condition="phenotype-only retrieval written; this stage has no interval argument"),
            _stage("retrieval-import",
                   {"phenotype": values.get("phenotype") or "", "case_id": values.get("case_id") or "",
                    "out_dir": str(out_dir / "c02_retrieval_freeze"),
                    "phenotype_context": values.get("phenotype_context")},
                   [],
                   inputs_from={"retrieval": ("retrieval-run", "retrieval")},
                   stop_condition="Phase 1 frozen; the interval may be opened only after this"),
            # --- Route B: interval-first, phenotype-blind ------------------------
            _stage("route-b-prepare",
                   {"candidates": values.get("candidates"), "proteins": values.get("proteins"),
                    "protein_db": values.get("protein_db"), "reference_function": values.get("reference_function"),
                    "orthology": values.get("orthology"), "domains": values.get("domains"),
                    "family_index": values.get("family_index"), "annotations": values.get("annotations"),
                    "case_id": values.get("case_id") or "", "out_dir": str(out_dir / "b01_route_b")},
                   ["candidates"],
                   optional_tools=["local BLASTP database"],
                   stop_condition="B0 complete; readiness reported as source_backed / partial / not_evaluated"),
            # B1 overlay: the B0 atlas is frozen before the phenotype is applied
            _stage("interval-reverse-lookup",
                   {"phenotype": values.get("phenotype"), "case_id": values.get("case_id") or "",
                    "shared_route": values.get("shared_route"),
                    "out_dir": str(out_dir / "b02_route_b1_overlay")},
                   ["phenotype"],
                   inputs_from={"atlas": ("route-b-prepare", "interval_gene_atlas"),
                                "reference_function": ("route-b-prepare", "reference_function")},
                   optional_chain_inputs=["reference_function"],
                   stop_condition="B1 phenotype overlay against the frozen B0 atlas; the atlas is not modified"),
            # B1 -> literature: what happens when an interval gene is broken
            _stage("reverse-literature",
                   {"phenotype": values.get("phenotype"), "cache": values.get("reverse_lit_cache"),
                    "allow_external_api": allow_reverse_literature, "case_id": values.get("case_id") or "",
                    "out_dir": str(out_dir / "b03_reverse_literature")},
                   [],
                   inputs_from={"lookup": ("route-b-prepare", "interval_reverse_lookup")},
                   external_transmission=("reference gene names and generic phenotype words (Europe PMC)"
                                          if allow_reverse_literature else "none (offline cache only)"),
                   stop_condition="query status retrieved/no_hit/not_evaluated_offline recorded; link candidates wait for judgment"),
            # --- Route C Phase 2: interval matching, after the freeze ------------
            _stage("retrieval-interval-match",
                   {"case_id": values.get("case_id") or "", "out_dir": str(out_dir / "c03_interval_match")},
                   [],
                   inputs_from={"retrieval": ("retrieval-import", "frozen"),
                                "atlas": ("route-b-prepare", "interval_gene_atlas")},
                   stop_condition="frozen retrieval matched against the interval atlas; the retrieval is not modified"),
            # --- the only place the three lanes meet ----------------------------
            _stage("triple-route-compare",
                   {"case_id": values.get("case_id") or "", "out_dir": str(out_dir / "z01_triple_route_compare")},
                   [],
                   inputs_from={"route_a": ("phenotype-first-status", "hypotheses"),
                                "route_b": ("route-b-prepare", "family_process_hypotheses"),
                                "route_c": ("retrieval-import", "frozen"),
                                "route_c_interval_match": ("retrieval-interval-match", "interval_match"),
                                "literature": ("reverse-literature", "literature"),
                                "literature_status": ("reverse-literature", "status")},
                   optional_chain_inputs=["literature", "literature_status"],
                   stop_condition="per-family failure-mode accounting emitted; convergence adds no score"),
            # --- after the meeting: canonical family candidate bridge ---
            # Normalizes the three lanes' family evidence (Route A verified-bait
            # family match, Route B orthology/atlas, Route C frozen retrieval
            # match) into one reporting table with a strict mapping tier. The
            # bridge never scores: rank_contribution is always none and it
            # writes no rank artifact. It runs after the comparison, so reading
            # lane outputs is allowed; every route input is optional because the
            # bridge runs as soon as at least one route produced evidence
            # (route-a-interval-match absent means empty family status ->
            # tier6 rows, never a miss). The bridge's route A hypotheses input
            # is phenotype_first_hypotheses.tsv (the only A output that carries
            # member_interval_genes).
            _stage("family-candidate-bridge",
                   {"case_id": values.get("case_id") or "", "out_dir": str(out_dir / "z01b_family_candidate_bridge")},
                   [],
                   inputs_from={
                       "route_a_hypotheses": ("phenotype-first-status", "hypotheses"),
                       "route_a_family_status": ("route-a-interval-match", "status"),
                        "route_a_bait_panel": ("bait-panel-freeze", "panel"),
                        "route_a_projection_panel": ("bait-panel-freeze", "projection_panel"),
                        "route_a_blast": ("blast-panel-search", "evidence"),
                        "route_a_blast_hits": ("blast-panel-search", "eligible_hits"),
                       "route_b_orthology": ("route-b-prepare", "orthology"),
                       "route_b_atlas": ("route-b-prepare", "interval_gene_atlas"),
                       # Positive-only family projection after the A/B freezes:
                       # B0 may generate this interval-scoped index locally.
                       # It can retain matching interval members for Route A,
                       # but cannot prove that an unmatched family is absent.
                       "route_b_family_index": ("route-b-prepare", "family_index"),
                       "route_c_retrieval": ("retrieval-import", "frozen"),
                       "route_c_match": ("retrieval-interval-match", "interval_match"),
                       "genes": ("route-b-prepare", "interval_gene_atlas"),
                       # Route B evidence-dependency clusters: one
                       # cluster per orthogroup/HOG observed by route-b-prepare.
                       # Optional: when no orthogroup/HOG data exists the
                       # upstream result is None and the bridge rows keep
                       # dependency_unknown -- nothing is fabricated. Clusters
                       # label dependency only and never add score.
                       "dependency_clusters": ("route-b-prepare", "dependency_clusters"),
                   },
                   optional_chain=True,
                   optional_chain_inputs=["route_a_hypotheses", "route_a_family_status", "route_a_bait_panel",
                                           "route_a_projection_panel",
                                           "route_a_blast",
                                           "route_a_blast_hits",
                                          "route_b_orthology", "route_b_atlas", "route_b_family_index",
                                          "route_c_retrieval",
                                          "route_c_match", "genes", "dependency_clusters"],
                   stop_condition=("three independently produced route artifacts normalized into one reporting table; "
                                   "rank_contribution stays none and no rank artifact is written")),
            # --- after the meeting: narrowing via the canonical handoff ---
            # Runs after the comparison, so reading the meeting output and the
            # frozen B0 atlas is allowed; orthology / synteny / expression /
            # target are never fabricated. If the comparison produced no handoff
            # the stage is upstream-unavailable, not a biological failure. The
            # bridge is an optional promotion gate: when the bridge stage did
            # not produce its table, the legacy behavior (no bridge columns)
            # applies unchanged.
            _stage("family-to-copy",
                   {"case_id": values.get("case_id") or "", "out_dir": str(out_dir / "z02_family_to_copy")},
                   [],
                   inputs_from={"families": ("triple-route-compare", "family_to_copy_input"),
                                "atlas": ("route-b-prepare", "interval_gene_atlas"),
                                "bridge": ("family-candidate-bridge", "bridge_tsv"),
                                "route_a_blast": ("blast-panel-search", "evidence"),
                                "route_a_blast_hits": ("blast-panel-search", "eligible_hits")},
                   optional_chain_inputs=["bridge", "route_a_blast", "route_a_blast_hits"],
                   stop_condition="canonical handoff narrowed by the existing family-to-copy logic; score unchanged"),
        ]
        # Route A blast-panel layer (#374 Criterion 3). Only planned when
        # interval proteins are supplied: the re-freeze and the search are both
        # meaningless without them, and leaving the search as a not_runnable
        # stage would mask an otherwise complete run. Inserted immediately
        # after the bait freeze so the frozen panel is re-frozen and searched
        # before any interval family matching reads the same panel.
        if values.get("proteins"):
            blast_stages = [
                # blast-panel-freeze: re-freezes the verified bait panel into
                # the auditable blast-panel schema. Still freezing -- the
                # interval proteins are only opened by blast-panel-search,
                # which verifies this freeze first.
                _stage(
                    "blast-panel-freeze",
                    {"case_id": values.get("case_id") or "", "out_dir": str(out_dir / "a05c_blast_panel_freeze")},
                    ["panel", "sequences"],
                    inputs_from={"panel": ("bait-panel-freeze", "projection_panel"),
                                 "sequences": ("hypothesis-verify", "reviewable_sequences")},
                    stop_condition="verified bait panel re-frozen into the blast-panel schema; the interval proteins are not opened here",
                ),
                # blast-panel-search: BLASTs only the verified frozen baits
                # against the interval proteins, after verifying the blast-panel
                # freeze. The interval proteins are opened by this stage; a
                # blast hit is evidence of similarity only (rank_contribution
                # stays none). Missing BLAST+ or a provider failure degrades to
                # not_evaluated_* with exit code 0, never a miss.
                _stage(
                    "blast-panel-search",
                    {"proteins": values.get("proteins"), "runner": "auto", "threads": 4,
                     "out_dir": str(out_dir / "a05d_blast_panel_search")},
                    ["freeze", "proteins"],
                    inputs_from={"freeze": ("blast-panel-freeze", "manifest")},
                    optional_tools=["BLAST+ (blastp/makeblastdb)"],
                    external_transmission="the frozen bait sequences are sent to the local BLAST+ runner only",
                    stop_condition="verified frozen baits BLASTed against the interval proteins post-freeze; BLAST+ absence or provider failure degrades to not_evaluated_* with exit code 0",
                ),
            ]
            anchor = next(index for index, stage in enumerate(stages)
                          if stage["command"] == "bait-panel-freeze")
            stages[anchor + 1:anchor + 1] = blast_stages
        if values.get("case") and values.get("resources"):
            stages.append(_stage(
                "case-run",
                {"case": values.get("case"), "resources": values.get("resources"),
                 "out_dir": str(out_dir / "z03_case_run"), "mode": "deterministic",
                 "bridge_candidates": values.get("candidates")},
                ["case", "resources"],
                inputs_from={
                    "narrowing_family_recovery": ("triple-route-compare", "family_to_copy_input"),
                    "narrowing_convergence": ("triple-route-compare", "triple_route_hypothesis_comparison"),
                    "narrowing_convergence_manifest": ("triple-route-compare", "convergence"),
                    "narrowing_family_copy": ("family-to-copy", "family_copy_discrimination"),
                    "narrowing_family_copy_shortlist": ("family-to-copy", "family_copy_shortlist"),
                    "narrowing_family_copy_manifest": ("family-to-copy", "manifest"),
                    "bridge_atlas": ("route-b-prepare", "interval_gene_atlas"),
                    # The same canonical comparison that feeds narrowing also
                    # feeds the causal-scenario ledger.  Route status artifacts
                    # are forwarded alongside it so a provider gap stays
                    # not_evaluated rather than being classified as B-only or
                    # neither inside case-run.
                    "discordance_comparison": ("triple-route-compare", "triple_route_hypothesis_comparison"),
                    "route_a_manifest": ("phenotype-first-status", "status"),
                    "route_b_manifest": ("route-b-prepare", "route_b0_status"),
                },
                stop_condition="case-run consumes triple-route family artifacts; biological rank remains unchanged",
            ))
        return stages, ("The three lanes were generated separately and meet only in "
                        "z01_triple_route_compare. Read triple_route_hypothesis_comparison.tsv: a family with "
                        "failure_mode_coverage=all_applicable_modes_neutralised has each supporting route's "
                        "characteristic error checked off. That lowers risk; it is not proof and adds no score.")

    if profile in _PREDICTION_ONLY_PROFILES:
        # Prediction-only: nothing in this lane may receive a truth/answer input.
        # `audit-predict` rejects an answer-shaped manifest column itself, and the
        # truth-dependent lanes are declared not_evaluated rather than attempted.
        allow_external = bool(values.get("allow_external_api"))
        stages = [
            _stage("query-eligibility",
                   {"query_manifest": values.get("query_manifest"), "query_fasta": values.get("query_fasta"),
                    "out_dir": str(out_dir / "01_query_eligibility")},
                   ["query_manifest"]),
            # --- Route A phenotype research, before any interval access ---
            _stage("phenotype-investigate",
                   {"case_id": values.get("case_id"), "focal_species": values.get("focal_species"),
                    "phenotype": values.get("phenotype"), "direction": values.get("direction") or "unknown",
                    "inheritance": values.get("inheritance") or "unknown",
                    "candidates": values.get("candidates"),
                    "phenotype_context": values.get("phenotype_context"),
                    "out_dir": str(out_dir / "02_phenotype_investigation")},
                   ["case_id", "focal_species", "phenotype"],
                   stop_condition="research packet emitted; run the prompt or author the payload, then import it"),
            _stage("phenotype-localization-packet",
                   {"candidates": values.get("candidates"),
                    "out_dir": str(out_dir / "02b_phenotype_localization_packet")},
                   [],
                   inputs_from={"investigation": ("phenotype-investigate", "investigation")},
                   stop_condition="phenotype-only localization prompt emitted before mechanism generation"),
            _stage("phenotype-localization-freeze",
                   {"payload": values.get("phenotype_localization_payload"),
                    "out_dir": str(out_dir / "02c_phenotype_localization_freeze")},
                   ["payload"],
                   inputs_from={"investigation": ("phenotype-investigate", "investigation")},
                   stop_condition="localization response validated and hash-frozen before mechanism generation"),
            _stage("mechanism-import",
                   {"payload": values.get("mechanism_payload"), "payload_source": values.get("payload_source") or "llm",
                    "model": values.get("model") or "", "case_id": values.get("case_id") or "",
                    "mechanism_freeze_manifest": values.get("mechanism_freeze_manifest"),
                    "candidates": values.get("candidates"), "out_dir": str(out_dir / "03_mechanism_import")},
                   ["payload"],
                   inputs_from={"investigation": ("phenotype-investigate", "investigation"),
                                "localization_freeze_manifest": ("phenotype-localization-freeze", "manifest")},
                   stop_condition="well-formed initial labeled hypotheses imported after localization freeze (target 12; count remains advisory)"),
            _stage("mechanism-coverage-audit",
                   {"gap_fill_mechanisms": values.get("gap_fill_mechanisms"),
                    "gap_fill_import_manifest": values.get("gap_fill_mechanism_import_manifest"),
                    "gap_request_audit": values.get("gap_request_audit"),
                    "candidates": values.get("candidates"),
                    "out_dir": str(out_dir / "03c_mechanism_coverage")},
                   (["gap_fill_mechanisms", "gap_fill_import_manifest", "gap_request_audit"]
                    if values.get("gap_fill_mechanisms") else []),
                   inputs_from={"initial_mechanisms": ("mechanism-import", "mechanisms"),
                                "initial_import_manifest": ("mechanism-import", "manifest"),
                                "localization_freeze_manifest": ("phenotype-localization-freeze", "manifest"),
                                "investigation": ("phenotype-investigate", "investigation")},
                   stop_condition="valid hypotheses continue with explicit degraded coverage; a blind gap prompt is advisory"),
            _stage("reference-gene-resolve",
                   {"allow_external_api": allow_external, "record_cache": values.get("record_cache"),
                    "focal_species": values.get("focal_species") or "",
                    "out_dir": str(out_dir / "03b_reference_gene_resolution")},
                   [],
                   inputs_from={"mechanisms": ("mechanism-coverage-audit", "mechanisms")},
                   external_transmission=("UniProt REST gene/accession resolution"
                                          if allow_external else "none (offline record cache only)"),
                   stop_condition=("every missing reference gene citation is provider-resolved or left "
                                   "review-required; malformed accessions are never queried")),
            _stage("hypothesis-verify",
                   {"focal_species": values.get("focal_species"), "record_cache": values.get("record_cache"),
                    "bound_claims": values.get("bound_claims"),
                    "allow_external_api": allow_external, "benchmark_mode": True,
                    "out_dir": str(out_dir / "04_hypothesis_verify")},
                   ["focal_species"],
                   inputs_from={"hypotheses": ("reference-gene-resolve", "mechanisms_resolved")},
                   external_transmission="UniProt REST lookups for accession/species/sequence verification" if allow_external else "none (offline record cache only)",
                   stop_condition="every proposal verified, rejected, or left review-required"),
            _stage("bait-panel-freeze",
                   {"case_id": values.get("case_id") or "", "out_dir": str(out_dir / "05_bait_panel_freeze")},
                   [],
                   inputs_from={"hypotheses": ("reference-gene-resolve", "mechanisms_resolved"),
                                "validation": ("hypothesis-verify", "validation"),
                                "sequences": ("hypothesis-verify", "verified_sequences"),
                                "reviewable_sequences": ("hypothesis-verify", "reviewable_sequences")},
                   stop_condition="panel frozen with SHA-256 and UTC time; a small panel stays review-required"),
            # --- deterministic prediction, still answerless ---
            _stage("audit-predict",
                   {"manifest": values.get("predict_manifest") or values.get("blind_manifest"),
                    "route": values.get("blind_route"), "case_id": values.get("case_id") or "",
                    "out_dir": str(out_dir / "06_prediction_freeze")},
                   ["manifest", "route"],
                   stop_condition="every case ranked and frozen, or a rejected (answer-carrying) manifest"),
            _stage("investigation-attribution",
                   {"case_id": values.get("case_id") or "", "ranked": values.get("candidates"),
                    "out_dir": str(out_dir / "07_attribution")},
                   [],
                   inputs_from={"mechanisms": ("reference-gene-resolve", "mechanisms_resolved"),
                                "validation": ("hypothesis-verify", "validation"),
                                "panel": ("bait-panel-freeze", "panel")},
                   optional_chain=True,
                   stop_condition="stage-by-stage attribution emitted; truth is never consulted"),
            _stage("ranking-qa",
                   {"candidates": values.get("candidates"), "evidence": values.get("evidence"),
                    "out_dir": str(out_dir / "08_ranking_qa"), "benchmark_mode": True},
                   ["candidates", "evidence"]),
        ]
        return stages, ("Predictions are frozen. Open the truth only in a separate process: "
                        "`candidate-qtl audit-score --freeze <out>/06_prediction_freeze/freeze_manifest.json "
                        "--truth <truth.tsv> --out-dir <scoring dir outside the freeze>`.")



def _execute_stage(stage: dict[str, object]) -> dict[str, object]:
    from .retrieval_local import retrieval_run
    from .retrieval_route import retrieval_import, retrieval_interval_match, triple_route_compare
    """Call the stable command with keyword arguments only.

    Positional forwarding silently changed meaning whenever a command grew a
    parameter, which is how guided screening ended up disabling local-BLAST
    discovery (#159). Keywords keep guided profiles pinned to the CLI defaults.
    """
    command = str(stage["command"])
    args = stage["arguments"]
    assert isinstance(args, dict)
    if command == "discover":
        return discover_candidates(
            fasta_arg=str(args["fasta"]),
            gff3_arg=str(args["gff3"]),
            intervals=list(args["interval"]),
            phenotype=str(args["phenotype"]),
            out_dir_arg=str(args["out_dir"]),
            protein_db_arg=str(args.get("protein_db") or "") or None,
            allow_external_api=bool(args["allow_external_api"]),
            auto_local_db=bool(args.get("auto_local_db", True)),
        )
    if command == "discover-promote":
        return promote_discover(
            discover_dir_arg=str(args["discover_dir"]), out_dir_arg=str(args["out_dir"]),
            case_id=str(args.get("case_id") or "") or None,
        )
    if command == "case-run":
        overrides = {
            "family_recovery": str(args.get("narrowing_family_recovery") or ""),
            "convergence": str(args.get("narrowing_convergence") or ""),
            "family_copy": str(args.get("narrowing_family_copy") or ""),
            "family_copy_shortlist": str(args.get("narrowing_family_copy_shortlist") or ""),
            "convergence_manifest": str(args.get("narrowing_convergence_manifest") or ""),
            "family_copy_manifest": str(args.get("narrowing_family_copy_manifest") or ""),
            "bridge_atlas": str(args.get("bridge_atlas") or ""),
            "bridge_candidates": str(args.get("bridge_candidates") or ""),
            "discordance_comparison": str(args.get("discordance_comparison") or ""),
            "route_a_manifest": str(args.get("route_a_manifest") or ""),
            "route_b_manifest": str(args.get("route_b_manifest") or ""),
        }
        return case_run(
            case_arg=str(args["case"]), resources_arg=str(args["resources"]),
            out_dir_arg=str(args["out_dir"]), mode="deterministic",
            narrowing_overrides={key: value for key, value in overrides.items() if value},
        )
    if command == "query-eligibility":
        return query_eligibility(
            manifest_arg=str(args["query_manifest"]), out_arg=str(args["out_dir"]),
            query_fasta_arg=str(args.get("query_fasta") or "") or None,
        )
    if command == "audit-predict":
        from .audit_freeze import audit_predict
        return audit_predict(
            manifest_arg=str(args["manifest"]), route_arg=str(args["route"]),
            out_arg=str(args["out_dir"]), case_id=str(args.get("case_id") or ""),
        )
    if command == "phenotype-investigate":
        return phenotype_investigate(
            out_dir_arg=str(args["out_dir"]), case_id=str(args["case_id"]),
            focal_species=str(args["focal_species"]), phenotype=str(args["phenotype"]),
            direction=str(args.get("direction") or "unknown"),
            inheritance=str(args.get("inheritance") or "unknown"),
            candidates_arg=str(args.get("candidates") or "") or None,
            phenotype_context=str(args.get("phenotype_context") or "") or None,
        )
    if command == "phenotype-localization-packet":
        return phenotype_localization_packet(
            str(args["investigation"]), str(args["out_dir"]),
            str(args.get("candidates") or "") or None,
        )
    if command == "phenotype-localization-freeze":
        return phenotype_localization_freeze(
            str(args["payload"]), str(args["investigation"]), str(args["out_dir"]),
        )
    if command == "mechanism-import":
        return mechanism_import(
            payload_arg=str(args["payload"]), out_dir_arg=str(args["out_dir"]),
            investigation_arg=str(args.get("investigation") or "") or None,
            payload_source=str(args.get("payload_source") or "llm"), model=str(args.get("model") or ""),
            case_id=str(args.get("case_id") or ""),
            mechanism_freeze_manifest=str(args.get("mechanism_freeze_manifest") or "") or None,
            localization_freeze_manifest=str(args.get("localization_freeze_manifest") or "") or None,
            candidates_arg=str(args.get("candidates") or "") or None,
        )
    if command == "mechanism-coverage-audit":
        mechanisms = [str(args["initial_mechanisms"])]
        manifests = [str(args["initial_import_manifest"])]
        if args.get("gap_fill_mechanisms"):
            mechanisms.append(str(args["gap_fill_mechanisms"]))
            manifests.append(str(args["gap_fill_import_manifest"]))
        return mechanism_coverage_audit(
            mechanisms, str(args["localization_freeze_manifest"]),
            str(args["investigation"]), str(args["out_dir"]),
            str(args.get("candidates") or "") or None,
            str(args.get("gap_request_audit") or "") or None,
            manifests,
        )
    if command == "reference-gene-resolve":
        return reference_gene_resolve(
            mechanisms_arg=str(args["mechanisms"]), out_dir_arg=str(args["out_dir"]),
            allow_external_api=bool(args.get("allow_external_api")),
            record_cache_arg=str(args.get("record_cache") or "") or None,
            focal_species=str(args.get("focal_species") or ""),
            identity_review_arg=str(args.get("identity_review") or "") or None,
        )
    if command == "hypothesis-verify":
        return hypothesis_verify(
            hypotheses_arg=str(args["hypotheses"]), out_dir_arg=str(args["out_dir"]),
            focal_species=str(args["focal_species"]),
            record_cache_arg=str(args.get("record_cache") or "") or None,
            allow_external_api=bool(args.get("allow_external_api")),
            benchmark_mode=bool(args.get("benchmark_mode")),
            source_records_arg=str(args.get("source_records") or "") or None,
            bound_claims_arg=str(args.get("bound_claims") or "") or None,
            agent_claim_review_arg=str(args.get("agent_claim_review") or "") or None,
            source_review_records_arg=str(args.get("source_review_records") or "") or None,
        )
    if command == "bait-panel-freeze":
        return bait_panel_freeze(
            hypotheses_arg=str(args["hypotheses"]), validation_arg=str(args["validation"]),
            sequences_arg=str(args["sequences"]), out_dir_arg=str(args["out_dir"]),
            case_id=str(args.get("case_id") or ""),
            reviewable_sequences_arg=str(args.get("reviewable_sequences") or "") or None,
        )
    if command == "blast-panel-freeze":
        # Guided Route A wiring: re-freeze the verified bait panel (a05) into
        # the blast-panel schema, then freeze it exactly like the CLI command.
        # The converted panel is an intermediate inside this stage's out-dir;
        # blast-panel-search consumes the freeze manifest produced here.
        converted = panel_from_verified_bait_panel(str(args["panel"]), str(args["sequences"]))
        out_dir = runtime_path(str(args["out_dir"]))
        out_dir.mkdir(parents=True, exist_ok=True)
        converted_panel = out_dir / "converted_bait_panel.tsv"
        converted.to_csv(converted_panel, sep="\t", index=False, lineterminator="\n", encoding="utf-8")
        return freeze_blast_panel(
            panel_arg=str(converted_panel), out_dir_arg=str(args["out_dir"]),
            case_id=str(args.get("case_id") or ""),
        )
    if command == "blast-panel-search":
        return run_panel_blastp(
            freeze_arg=str(args["freeze"]), proteins_arg=str(args["proteins"]),
            out_dir_arg=str(args["out_dir"]),
            num_threads=int(args.get("threads") or 4),
            runner=str(args.get("runner") or "auto"),
        )
    if command == "route-a-interval-match":
        return route_a_interval_family_match(
            panel_arg=str(args["panel"]) if args.get("panel") else "",
            hypotheses_arg=str(args["hypotheses"]) if args.get("hypotheses") else "",
            out_dir_arg=str(args["out_dir"]),
            family_index_arg=str(args.get("family_index") or "") or None,
            intervals=list(args.get("interval") or []),
            candidates_arg=str(args.get("candidates") or "") or None,
            rescue_arg=str(args.get("annotation_rescue") or "") or None,
            case_id=str(args.get("case_id") or ""),
        )
    if command == "source-verify":
        return source_verify(
            mechanisms_arg=str(args["mechanisms"]), out_dir_arg=str(args["out_dir"]),
            cache_arg=str(args.get("cache") or "") or None,
            allow_external_api=bool(args.get("allow_external_api")),
            claim_extraction_arg=str(args.get("claim_extraction") or "") or None,
        )
    if command == "phenotype-first-status":
        return phenotype_first_status(
            out_dir_arg=str(args["out_dir"]),
            hypotheses_arg=str(args.get("hypotheses") or "") or None,
            validation_arg=str(args.get("validation") or "") or None,
            panel_arg=str(args.get("panel") or "") or None,
            falsification_arg=str(args.get("falsification") or "") or None,
            case_id=str(args.get("case_id") or ""),
            phenotype_context=str(args.get("phenotype_context") or "") or None,
            phenotype_text=str(args.get("phenotype") or ""),
        )
    if command == "route-b-prepare":
        return route_b_prepare(
            candidates_arg=str(args["candidates"]), out_dir_arg=str(args["out_dir"]),
            proteins_arg=str(args.get("proteins") or "") or None,
            protein_db_arg=str(args.get("protein_db") or "") or None,
            reference_function_args=[str(args["reference_function"])] if args.get("reference_function") else None,
            orthology_arg=str(args.get("orthology") or "") or None,
            domains_arg=str(args.get("domains") or "") or None,
            family_index_arg=str(args.get("family_index") or "") or None,
            annotations_arg=str(args.get("annotations") or "") or None,
            case_id=str(args.get("case_id") or ""),
        )
    if command == "interval-reverse-lookup":
        return interval_reverse_lookup(
            atlas_arg=str(args["atlas"]), out_dir_arg=str(args["out_dir"]),
            reference_function_arg=str(args.get("reference_function") or "") or None,
            phenotype=str(args.get("phenotype") or ""),
            case_id=str(args.get("case_id") or ""),
            shared_route_arg=str(args.get("shared_route") or "") or None,
        )
    if command == "reverse-literature":
        return reverse_literature(
            lookup_arg=str(args["lookup"]), out_dir_arg=str(args["out_dir"]),
            phenotype=str(args.get("phenotype") or ""),
            cache_arg=str(args.get("cache") or "") or None,
            allow_external_api=bool(args.get("allow_external_api")),
            case_id=str(args.get("case_id") or ""),
        )
    if command == "retrieval-run":
        return retrieval_run(
            index_arg=str(args["index"]), phenotype=str(args["phenotype"]),
            out_dir_arg=str(args["out_dir"]), case_id=str(args.get("case_id") or ""),
            phenotype_context=str(args.get("phenotype_context") or "") or None,
        )
    if command == "retrieval-import":
        return retrieval_import(
            retrieval_arg=str(args["retrieval"]), out_dir_arg=str(args["out_dir"]),
            phenotype=str(args.get("phenotype") or ""), case_id=str(args.get("case_id") or ""),
            phenotype_context=str(args.get("phenotype_context") or "") or None,
        )
    if command == "retrieval-interval-match":
        return retrieval_interval_match(
            retrieval_arg=str(args["retrieval"]), atlas_arg=str(args["atlas"]),
            out_dir_arg=str(args["out_dir"]), case_id=str(args.get("case_id") or ""),
        )
    if command == "triple-route-compare":
        return triple_route_compare(
            route_a_arg=str(args["route_a"]), route_b_arg=str(args["route_b"]),
            route_c_arg=str(args["route_c"]), out_dir_arg=str(args["out_dir"]),
            case_id=str(args.get("case_id") or ""),
            route_c_match_arg=str(args.get("route_c_interval_match") or "") or None,
            literature_arg=str(args.get("literature") or "") or None,
            literature_status_arg=str(args.get("literature_status") or "") or None,
        )
    if command == "family-candidate-bridge":
        # Issue #374 review: every route input is optional in this profile, so
        # a run where no lane produced any evidence would hand the bridge eight
        # empty paths and make it fail with "requires at least one route input".
        # That is an unavailable stage, not a defect: record it as not_runnable
        # (same vocabulary as the plan-time / upstream-missing states) without
        # calling the command. With at least one input the command runs as-is.
        route_inputs = (
            str(args.get("route_a_hypotheses") or ""), str(args.get("route_a_family_status") or ""),
            str(args.get("route_a_bait_panel") or ""), str(args.get("route_a_blast") or ""),
            str(args.get("route_a_projection_panel") or ""),
            str(args.get("route_a_blast_hits") or ""),
            str(args.get("route_b_orthology") or ""), str(args.get("route_b_atlas") or ""),
            str(args.get("route_b_family_index") or ""),
            str(args.get("route_c_retrieval") or ""), str(args.get("route_c_match") or ""),
            str(args.get("genes") or ""),
        )
        if not any(route_inputs):
            return {
                "guided_stage_status": "not_runnable_no_route_input",
                "reason": "family-candidate-bridge requires at least one route input; "
                          "no route lane produced evidence in this run",
            }
        return family_candidate_bridge(
            out_dir_arg=str(args["out_dir"]),
            route_a_hypotheses_arg=str(args.get("route_a_hypotheses") or "") or None,
            route_a_family_status_arg=str(args.get("route_a_family_status") or "") or None,
            route_a_bait_panel_arg=str(args.get("route_a_bait_panel") or "") or None,
            route_a_projection_panel_arg=str(args.get("route_a_projection_panel") or "") or None,
            route_a_blast_arg=str(args.get("route_a_blast") or "") or None,
            route_a_blast_eligible_arg=str(args.get("route_a_blast_hits") or "") or None,
            route_b_orthology_arg=str(args.get("route_b_orthology") or "") or None,
            route_b_atlas_arg=str(args.get("route_b_atlas") or "") or None,
            route_b_family_index_arg=str(args.get("route_b_family_index") or "") or None,
            route_c_retrieval_arg=str(args.get("route_c_retrieval") or "") or None,
            route_c_match_arg=str(args.get("route_c_match") or "") or None,
            dependency_clusters_arg=str(args.get("dependency_clusters") or "") or None,
            genes_arg=str(args.get("genes") or "") or None,
            case_id=str(args.get("case_id") or ""),
        )
    if command == "family-to-copy":
        # Issue #284 minimal chaining: canonical handoff + the frozen B0 atlas
        # only. Missing orthology / synteny / expression / target stay absent and
        # the existing unresolved / review semantics of family-to-copy apply.
        # Issue #374: the family-candidate-bridge table (when produced) gates
        # shortlist promotion; without it the legacy behavior is byte-identical.
        return family_to_copy(
            families_arg=str(args["families"]), out_dir_arg=str(args["out_dir"]),
            atlas_arg=str(args.get("atlas") or "") or None,
            case_id=str(args.get("case_id") or ""),
            bridge_arg=str(args.get("bridge") or "") or None,
            route_a_blast_arg=str(args.get("route_a_blast") or "") or None,
            route_a_blast_hits_arg=str(args.get("route_a_blast_hits") or "") or None,
        )
    if command == "investigation-attribution":
        return investigation_attribution(
            out_dir_arg=str(args["out_dir"]),
            mechanisms_arg=str(args.get("mechanisms") or "") or None,
            validation_arg=str(args.get("validation") or "") or None,
            panel_arg=str(args.get("panel") or "") or None,
            interval_match_arg=str(args.get("interval_match") or "") or None,
            ranked_arg=str(args.get("ranked") or "") or None,
            case_id=str(args.get("case_id") or ""),
        )
    if command == "ranking-qa":
        return ranking_qa(
            candidates_arg=str(args["candidates"]), evidence_arg=str(args["evidence"]),
            out_dir_arg=str(args["out_dir"]), benchmark_mode=True,
        )
    raise ValueError(f"unsupported guided stage: {command}")


def _resolved_provider(command: str, result: dict[str, object]) -> dict[str, object]:
    """Record which annotation provider a stage actually used, so the manifest
    shows the resolved provider next to the declared policy (#159)."""
    if command != "discover":
        return {}
    return {
        "annotation_provider": result.get("annotation_provider"),
        "local_protein_db": result.get("local_protein_db"),
        "homology_searched": result.get("homology_searched"),
    }


def _overall_status(results: list[dict[str, object]]) -> str:
    statuses = {str(record.get("status")) for record in results}
    if not statuses:
        return "no_stages"
    if "interrupted" in statuses:
        return "interrupted"
    completed = "completed" in statuses
    incomplete = bool(statuses - {"completed"})
    if not incomplete:
        return "completed"
    return "partial" if completed else "failed"


def guided_profile(profile: str, out_dir_arg: str, *, execute: bool = False, **values: Any) -> dict[str, object]:
    if profile not in PROFILE_NAMES:
        raise ValueError(f"profile must be one of: {', '.join(PROFILE_NAMES)}")
    if values.get("agent_claim_review") and profile != "triple-route":
        raise ValueError("guided agent claim review requires the triple-route source-verify chain")
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    run_id = _new_run_id()
    manifest_path = out_dir / _GUIDED_PROFILE_MANIFEST
    superseded_manifest = _archive_superseded_manifest(out_dir, manifest_path)
    stages, promotion_path = _plan(profile, values, out_dir)
    independent_lanes = profile in _INDEPENDENT_LANE_PROFILES
    manifest: dict[str, object] = {
        "schema_version": "1.2", "profile": profile, "candidate_qtl_version": __version__,
        "run_id": run_id, "final": False,
        "execution_requested": execute, "dry_run": not execute, "stages": stages,
        "review_policy": "review-required inputs are never auto-filled or implicitly approved",
        "ranking_policy": "no profile-specific weights or hidden ranking logic; stages call stable commands",
        "failure_policy": (
            "each lane is independent: a failed lane is recorded and the remaining lanes still run"
            if independent_lanes else
            "fail-fast: the first stage failure stops the profile, and the manifest is written before stopping"
        ),
        "screening_to_full_case": promotion_path,
        "overall_status": "dry_run",
    }
    if superseded_manifest:
        # The previous manifest survives under the archive path and points back
        # here via superseded_by; this run records the old version it replaced
        # so the pair is machine-determinable.
        manifest["supersedes"] = superseded_manifest
    if profile in _PREDICTION_ONLY_PROFILES:
        # A lane that needs truth is unavailable here by design. Recording it as
        # not_evaluated keeps it distinct from a lane that ran and failed.
        manifest["prediction_only"] = True
        manifest["truth_inputs_accepted"] = False
        manifest["not_evaluated_lanes"] = {
            lane: "truth-dependent; unavailable before the freeze is scored by `audit-score`"
            for lane in TRUTH_DEPENDENT_LANES
        }
        manifest["post_freeze_scoring_command"] = "audit-score"
    if profile in _TRIPLE_ROUTE_PROFILES:
        manifest["route_independence"] = {
            "kind": "procedural",
            "rule": "no lane reads another lane's output before the comparison; enforced by path check per stage",
            "route_a": "phenotype-first; generated without any interval, atlas or candidate input (the packet and payload are leak-audited)",
            "route_b": "interval-first; the atlas is built and frozen phenotype-free (B0), and only then do the B1 phenotype overlay and the reverse-literature lookups run against those frozen artifacts",
            "route_c": "phenotype-only retrieval frozen in Phase 1 before any interval artifact is opened; matched only in Phase 2",
            "meeting_point": ("triple-route-compare is the only stage that reads more than one lane; the "
                              "family-candidate-bridge stage after it normalizes the lanes' family evidence into one "
                              "reporting table, and the family-to-copy stage reads the meeting output, the frozen B0 "
                              "atlas and the bridge, still adding no score"),
            "not_claimed": (
                "A and B traverse the same literature, so their agreement is corroboration over a shared evidence "
                "base and is never reported as independent confirmation"
            ),
        }
        manifest["route_c_status"] = "experimental"
        manifest["route_c_rule"] = "Route C is experimental: it may inform review but must not affect candidate ordering"
    rank_preflight = _plan_time_rank_lane_preflight(values)
    if rank_preflight is not None:
        manifest["rank_lane_preflight"] = rank_preflight
    payload_preflight = _plan_time_payload_schema_preflight(values)
    if payload_preflight is not None:
        manifest["payload_schema_preflight"] = payload_preflight
        if payload_preflight["status"] == "rejected":
            manifest["overall_status"] = "failed"
            _write_manifest(manifest_path, manifest)
            return {"profile": profile, "dry_run": not execute, "manifest": str(manifest_path),
                    "stage_count": len(stages), "results": [], "overall_status": "failed",
                    "failed_stages": [], "exit_code": 1}
    manifest["structured_context"] = structured_context_block(
        values.get("phenotype_context"),
        case_id=str(values.get("case_id") or ""),
        phenotype_text=str(values.get("phenotype") or ""),
        query_usage={
            "usage": "stored_only",
            "reason": "profile-level audit record; per-stage manifests carry each consumer stage's own query_usage",
        },
    )
    results: list[dict[str, object]] = []
    if not execute:
        _write_manifest(manifest_path, manifest)
        return {"profile": profile, "dry_run": True, "manifest": str(manifest_path),
                "stage_count": len(stages), "results": results, "overall_status": "dry_run"}

    manifest["execution_results"] = results
    manifest["overall_status"] = "running"
    _write_manifest(manifest_path, manifest)
    stop_after_failure = False
    for stage in stages:
        command = str(stage["command"])
        if stop_after_failure:
            results.append({"command": command, "status": "skipped_after_failure",
                            "note": "a preceding stage failed in a fail-fast profile"})
            _write_manifest(manifest_path, manifest)
            continue
        if stage["status"] != "ready":
            results.append({"command": command, "status": stage["status"], "missing_inputs": stage["missing_inputs"]})
            _write_manifest(manifest_path, manifest)
            continue
        if profile in _TRIPLE_ROUTE_PROFILES:
            _assert_lane_isolation(stage, out_dir)
        blocked = _resolve_chained(stage, results)
        if blocked:
            stage["status"] = blocked
            results.append({"command": command, "status": blocked,
                            "note": "an earlier stage in this profile did not produce this stage's input; the lane is unavailable, not failed"})
            _write_manifest(manifest_path, manifest)
            continue
        record: dict[str, object] = {"command": command, "status": "running", "started_at": _now()}
        results.append(record)
        _write_manifest(manifest_path, manifest)
        try:
            result = _execute_stage(stage)
        except (KeyboardInterrupt, SystemExit):
            # Operator-driven termination is not a stage defect; record it and re-raise.
            record.update({"status": "interrupted", "finished_at": _now()})
            manifest["overall_status"] = "interrupted"
            _write_manifest(manifest_path, manifest)
            raise
        except Exception as error:  # noqa: BLE001 - one lane's defect must not erase the audit record
            record.update({"status": "failed", "finished_at": _now(),
                           "error_type": type(error).__name__, "error": str(error)[:2000]})
            stop_after_failure = not independent_lanes
        else:
            # A stage may decline to run after all (e.g. every optional route
            # input of the family-candidate-bridge stayed empty). The decline is
            # recorded as the stage's status -- not_runnable vocabulary, never a
            # failure -- the same way the plan-time and upstream-missing states
            # are recorded above.
            guided_status = (str(result.get("guided_stage_status"))
                             if isinstance(result, dict) and result.get("guided_stage_status") else "")
            if guided_status:
                stage["status"] = guided_status
                record.update({"status": guided_status, "finished_at": _now(),
                               "note": str(result.get("reason", "")), "result": result})
            else:
                record.update({"status": "completed", "finished_at": _now(), "result": result})
            record.update(_resolved_provider(command, result))
        _write_manifest(manifest_path, manifest)
    # Issue #346: an executed case-run stage owns the authoritative rank-lane
    # state (row-level evidence counts); its result replaces the plan-time
    # preflight so the final manifest never claims a plan-time guess as fact.
    for record in results:
        if record.get("command") == "case-run" and record.get("status") == "completed":
            executed_rank = record.get("result") or {}
            if "rank_semantic_status" in executed_rank:
                preflight = dict(manifest.get("rank_lane_preflight") or {"plan_time_detection": False})
                preflight.update({
                    "plan_time_detection": False,
                    "rank_semantic_status": executed_rank["rank_semantic_status"],
                    "rank_is_positional_fallback": executed_rank["rank_is_positional_fallback"],
                    "rank_metric_eligible": executed_rank["rank_metric_eligible"],
                    "stage_a_is_budget_only": executed_rank["stage_a_is_budget_only"],
                    "warning": executed_rank.get("rank_lane_warning", ""),
                })
                # Issue #397: the executed case-run declares rank_is_penalty_only
                # only when true; the preflight never guesses it at plan time, so
                # it is propagated here only when the executed run set it.
                if executed_rank.get("rank_is_penalty_only"):
                    preflight["rank_is_penalty_only"] = True
                manifest["rank_lane_preflight"] = preflight
    manifest["overall_status"] = _overall_status(results)
    # A run is the final version only when every stage completed; failed,
    # partial, interrupted and dry-run records stay final=false.
    manifest["final"] = manifest["overall_status"] == "completed"
    _write_manifest(manifest_path, manifest)
    failed = [str(record["command"]) for record in results if record.get("status") == "failed"]
    return {"profile": profile, "dry_run": False, "manifest": str(manifest_path), "stage_count": len(stages),
            "results": results, "overall_status": manifest["overall_status"],
            "failed_stages": failed, "exit_code": 1 if failed else 0}
