"""Phenotype-only causal-layer coverage contract for blind Route A (issue #658)."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from . import __version__
from .core import runtime_path, sha256
from .relevance_tiers import INITIAL_FIELDS
from .llm_pathway import HYPOTHESIS_COLUMNS, _blind_leak_audit
from .phenotype_investigation import (
    CAUSAL_LAYERS,
    MECHANISM_COLUMNS,
    MECHANISM_PAYLOAD_TEMPLATE,
    MIN_COMPETING_MECHANISMS,
    MIN_FAMILY_NODE_PORTFOLIO,
    MAX_FAMILY_NODE_PORTFOLIO,
    TARGET_FAMILY_NODE_PORTFOLIO,
)

LOCALIZATIONS = ("broad", "upstream", "signaling", "regulation", "execution", "unknown")
LOCALIZATION_CONFIDENCES = ("high", "medium", "low")
CORE_LAYERS = ("perception", "signaling", "regulation", "execution")
LOCALIZED_LAYER = {
    "upstream": "perception",
    "signaling": "signaling",
    "regulation": "regulation",
    "execution": "execution",
}
ADJACENT_LAYER = {
    "upstream": "signaling",
    "signaling": "regulation",
    "regulation": "signaling",
    "execution": "regulation",
}
LOCALIZATION_PROMPT_VERSION = "phenotype-localization-1.0.0"
COVERAGE_CONTRACT_VERSION = "route-a-causal-coverage-1.1.0"


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _normal(value: object) -> str:
    return " ".join(str(value or "").casefold().split())


def _read_investigation(path_arg: str) -> dict[str, str]:
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    if len(frame) != 1:
        raise ValueError("phenotype investigation must contain exactly one row")
    return {str(key): str(value) for key, value in frame.iloc[0].to_dict().items()}


def localization_prompt(investigation: dict[str, str]) -> str:
    """Build a phenotype-only localization prompt with no hypothesis context."""
    return f"""You are localizing a plant phenotype before any mechanism hypotheses are generated.
You may use only the phenotype packet below. Do not name genes or gene families. Do not use a QTL
interval, candidate list, Route B output, benchmark answer, score, or previous generation.

Choose causal_localization from: {', '.join(LOCALIZATIONS)}.
Choose localization_confidence from: {', '.join(LOCALIZATION_CONFIDENCES)}.
supporting_phrase must be a verbatim non-empty phrase copied from the Phenotype field.

Case: {investigation.get('case_id', '')}
Focal species: {investigation.get('focal_species', '')}
Phenotype: {investigation.get('phenotype', '')}
Direction: {investigation.get('direction', '')}
Inheritance: {investigation.get('inheritance', '')}
Affected tissue/organ: {investigation.get('affected_tissue', '')}
Developmental stage: {investigation.get('developmental_stage', '')}
Onset: {investigation.get('onset', '')}
Penetrance/expressivity: {investigation.get('penetrance', '')}
Stated uncertainty: {investigation.get('uncertainty', '')}

Return JSON only:
{{"causal_localization":"", "localization_confidence":"", "supporting_phrase":""}}
"""


def phenotype_localization_packet(investigation_arg: str, out_dir_arg: str,
                                  candidates_arg: str | None = None) -> dict[str, object]:
    investigation = _read_investigation(investigation_arg)
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    prompt = localization_prompt(investigation)
    leak = _blind_leak_audit(prompt, candidates_arg)
    if leak["leaked_identifiers"]:
        raise ValueError("phenotype localization packet contains QTL interval candidate identifier(s): "
                         + ", ".join(str(value) for value in leak["leaked_identifiers"][:5]))
    prompt_path = out / "phenotype_localization_prompt.txt"
    prompt_path.write_text(prompt, encoding="utf-8")
    manifest = {
        "schema_version": "1.0", "candidate_qtl_version": __version__, "frozen_at_utc": _now(),
        "case_id": investigation.get("case_id", ""), "prompt_version": LOCALIZATION_PROMPT_VERSION,
        "prompt_sha256": hashlib.sha256(prompt.encode("utf-8")).hexdigest(),
        "interval_blind": True, "leak_audit": leak,
        "rule": "run this fresh phenotype-only localization before mechanism generation, then freeze its response",
        "outputs": {"prompt": str(prompt_path)},
    }
    manifest_path = out / "phenotype_localization_packet_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"prompt": str(prompt_path), "manifest": str(manifest_path), "interval_blind": True}


def phenotype_localization_freeze(payload_arg: str, investigation_arg: str, out_dir_arg: str) -> dict[str, object]:
    """Validate and hash-freeze the fresh localization response before generation."""
    payload_path = runtime_path(payload_arg)
    investigation = _read_investigation(investigation_arg)
    try:
        value = json.loads(payload_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise ValueError(f"localization payload is not valid JSON: {payload_path}") from error
    required = {"causal_localization", "localization_confidence", "supporting_phrase"}
    if not isinstance(value, dict) or set(value) != required:
        raise ValueError("localization payload must be one JSON object with exactly: " + ", ".join(sorted(required)))
    location = _normal(value["causal_localization"])
    confidence = _normal(value["localization_confidence"])
    phrase = str(value["supporting_phrase"]).strip()
    if location not in LOCALIZATIONS:
        raise ValueError("causal_localization must be one of: " + ", ".join(LOCALIZATIONS))
    if confidence not in LOCALIZATION_CONFIDENCES:
        raise ValueError("localization_confidence must be one of: " + ", ".join(LOCALIZATION_CONFIDENCES))
    if not phrase or phrase not in str(investigation.get("phenotype", "")):
        raise ValueError("supporting_phrase must be a verbatim non-empty phrase from the phenotype field")
    frozen = {"causal_localization": location, "localization_confidence": confidence,
              "supporting_phrase": phrase}
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    frozen_path = out / "phenotype_localization.frozen.json"
    manifest_path = out / "phenotype_localization_freeze_manifest.json"
    prompt_hash = hashlib.sha256(localization_prompt(investigation).encode("utf-8")).hexdigest()
    if frozen_path.exists() or manifest_path.exists():
        try:
            existing, existing_manifest = _load_frozen_localization(str(manifest_path))
            if (existing != frozen
                    or existing_manifest.get("case_id") != investigation.get("case_id", "")
                    or existing_manifest.get("source_payload_sha256") != sha256(payload_path)
                    or existing_manifest.get("prompt_sha256") != prompt_hash
                    or existing_manifest.get("prompt_version") != LOCALIZATION_PROMPT_VERSION
                    or runtime_path(str(existing_manifest.get("frozen_localization", ""))).resolve()
                    != frozen_path.resolve()):
                raise ValueError("existing localization freeze does not match this payload and phenotype packet")
        except (OSError, ValueError, TypeError) as error:
            raise ValueError(f"existing localization freeze cannot be replaced: {error}") from error
        return {"localization": str(frozen_path), "manifest": str(manifest_path), **existing}
    frozen_path.write_text(json.dumps(frozen, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = {
        "schema_version": "1.0", "candidate_qtl_version": __version__, "frozen_at_utc": _now(),
        "case_id": investigation.get("case_id", ""), "prompt_version": LOCALIZATION_PROMPT_VERSION,
        "prompt_sha256": prompt_hash,
        "source_payload": str(payload_path), "source_payload_sha256": sha256(payload_path),
        "frozen_localization": str(frozen_path), "frozen_localization_sha256": sha256(frozen_path),
        "interval_informed": False,
        "rule": "this localization is frozen before mechanism generation and never changes category counts or scores",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"localization": str(frozen_path), "manifest": str(manifest_path), **frozen}


def _load_frozen_localization(manifest_arg: str) -> tuple[dict[str, str], dict[str, object]]:
    manifest_path = runtime_path(manifest_arg)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict) or manifest.get("interval_informed") is not False:
        raise ValueError("localization freeze manifest must explicitly declare interval_informed=false")
    frozen_value = Path(str(manifest.get("frozen_localization", "")))
    frozen_path = runtime_path(frozen_value if frozen_value.is_absolute() else manifest_path.parent / frozen_value)
    expected = str(manifest.get("frozen_localization_sha256", ""))
    if not frozen_path.is_file() or not expected or sha256(frozen_path) != expected:
        raise ValueError("frozen localization hash mismatch")
    value = json.loads(frozen_path.read_text(encoding="utf-8"))
    required = {"causal_localization", "localization_confidence", "supporting_phrase"}
    if not isinstance(value, dict) or set(value) != required:
        raise ValueError("frozen localization must contain exactly the localization contract fields")
    if _normal(value["causal_localization"]) not in LOCALIZATIONS:
        raise ValueError("frozen causal_localization is invalid")
    if _normal(value["localization_confidence"]) not in LOCALIZATION_CONFIDENCES:
        raise ValueError("frozen localization_confidence is invalid")
    return {key: str(value[key]) for key in ("causal_localization", "localization_confidence", "supporting_phrase")}, manifest


def _requested_layers(counts: Counter[str], localization: dict[str, str], total: int) -> tuple[list[str], list[str]]:
    location = localization["causal_localization"]
    confidence = localization["localization_confidence"]
    reasons: list[str] = []
    requested: list[str] = []
    if confidence == "low" or location in {"broad", "unknown"}:
        missing = [layer for layer in CORE_LAYERS if counts[layer] == 0]
        requested.extend(missing)
        if missing:
            reasons.append("broad_or_uncertain_core_layer_gap")
        dominant = max((counts[layer] for layer in CAUSAL_LAYERS), default=0)
        if total and dominant / total > 0.60:
            reasons.append("single_layer_dominance_above_60_percent")
            for layer in sorted(CORE_LAYERS, key=lambda item: (counts[item], CORE_LAYERS.index(item))):
                if layer not in requested and counts[layer] < dominant:
                    requested.append(layer)
    elif confidence == "high" and location in LOCALIZED_LAYER:
        adjacent = ADJACENT_LAYER[location]
        if counts[adjacent] == 0:
            requested.append(adjacent)
            reasons.append("high_confidence_adjacent_layer_gap")
        alternative_gap = max(0, 2 - counts["alternative"])
        requested.extend(["alternative"] * alternative_gap)
        if alternative_gap:
            reasons.append("high_confidence_alternative_gap")
    elif confidence == "medium" and location in LOCALIZED_LAYER:
        focal = LOCALIZED_LAYER[location]
        other_count = sum(count for layer, count in counts.items() if layer not in {focal, "unclear"})
        gap = max(0, 2 - other_count)
        if gap:
            choices = [ADJACENT_LAYER[location], "alternative", *CORE_LAYERS]
            requested.extend([layer for layer in choices if layer != focal][:gap])
            reasons.append("medium_confidence_other_layer_gap")
    return requested, reasons


def _pathway_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    return [{
        "hypothesis_id": row["hypothesis_id"], "case_id": row["case_id"],
        "pathway_or_process": row["pathway_or_process"], "gene_family": row["gene_family"],
        "reference_gene_name": row["reference_gene_name"], "reference_species": row["reference_species"],
        "proposed_stable_id": row["proposed_stable_id"], "expected_role": row["expected_role"],
        "mechanistic_link": row["mechanism_statement"], "counterevidence": row["counter_hypothesis"],
        "confidence": row["confidence"], "uncertainty": row["uncertainty"],
        "proposal_source": f"phenotype_investigation:{row['payload_source']}", "model": row["model"],
        "prompt_version": row["prompt_version"], "parser_version": row["parser_version"],
        "execution_mode": "blind", "review_status": row["review_status"], "notes": row["notes"],
    } for row in rows]


def _gap_prompt(investigation: dict[str, str], rows: list[dict[str, object]], requested: list[str]) -> str:
    family_names = [str(row["gene_family"]) for row in rows]
    return f"""Generate exactly {len(requested)} complementary plant mechanism hypotheses.
Use only the phenotype packet, requested causal layers, and already-used family names below.
Do not use an interval, candidate list, Route B, truth, scores, benchmark history, or prior prose.
Each new row must use a different family not in the existing list. Requested labels are positional:
{json.dumps(requested)}

Case: {investigation.get('case_id', '')}
Focal species: {investigation.get('focal_species', '')}
Phenotype: {investigation.get('phenotype', '')}
Direction: {investigation.get('direction', '')}
Inheritance: {investigation.get('inheritance', '')}
Affected tissue/organ: {investigation.get('affected_tissue', '')}
Developmental stage: {investigation.get('developmental_stage', '')}
Stated uncertainty: {investigation.get('uncertainty', '')}
Already-used family names: {json.dumps(family_names, ensure_ascii=False)}

Use the field names in the JSON template below; no earlier payload or conversation is needed.
For causal_layer choose the requested positional label; for confidence choose low, medium, or high.
Supply a real literature-supported reference gene/locus and species when known. Prioritize bridging
functional and mechanistic evidence from other species. A focal-species primary reference is allowed
as primary evidence in production; its absence must not block a hypothesis. Do not optimize
for a particular species, example or benchmark answer. An accession is optional. A strict benchmark route
may require a non-focal reference explicitly.
Never invent genes, accessions or citations. Leave unknown fields empty and explain the gap in uncertainty.
source_ids must be a string of PMID:/DOI:/UniProt:/TAIR: identifiers separated by semicolons, or empty.
Include a counter_hypothesis and discriminating_observation. Do not emit hypothesis_id: the importer
assigns round-specific IDs. Labels and confidence never contribute to a score.

Return JSON only:
{MECHANISM_PAYLOAD_TEMPLATE}
"""


def _preserve_artifact(path: Path, content: bytes) -> None:
    if path.exists():
        if path.read_bytes() != content:
            raise ValueError(f"preserved artifact content mismatch: {path}")
    else:
        path.write_bytes(content)


def mechanism_coverage_audit(mechanisms_args: list[str], localization_manifest_arg: str,
                             investigation_arg: str, out_dir_arg: str,
                             candidates_arg: str | None = None,
                             gap_request_audit_arg: str | None = None,
                             mechanism_import_manifests_args: list[str] | None = None) -> dict[str, object]:
    """Merge the initial blind portfolio plus at most three blind gap rows."""
    if not mechanisms_args or len(mechanisms_args) > 2:
        raise ValueError("supply one initial mechanism table and optionally one gap-fill table")
    if not mechanism_import_manifests_args or len(mechanism_import_manifests_args) != len(mechanisms_args):
        raise ValueError("supply one mechanism-import manifest for every mechanism table")
    investigation = _read_investigation(investigation_arg)
    localization, localization_manifest = _load_frozen_localization(localization_manifest_arg)
    if str(localization_manifest.get("case_id", "")).strip() != str(investigation.get("case_id", "")).strip():
        raise ValueError("localization freeze case_id does not match the phenotype investigation")
    expected_prompt_hash = hashlib.sha256(localization_prompt(investigation).encode("utf-8")).hexdigest()
    if localization_manifest.get("prompt_sha256") != expected_prompt_hash:
        raise ValueError("localization freeze does not bind the supplied phenotype packet")
    if localization["supporting_phrase"] not in str(investigation.get("phenotype", "")):
        raise ValueError("frozen supporting_phrase is not verbatim in the phenotype investigation")
    frames = [pd.read_csv(runtime_path(path), sep="\t", dtype=str).fillna("") for path in mechanisms_args]
    if not len(frames[0]):
        raise ValueError("initial mechanism portfolio must contain at least one well-formed row")
    if len(frames) == 2 and len(frames[1]) > 3:
        raise ValueError(f"gap-fill portfolio may contain at most 3 rows; got {len(frames[1])}")
    for frame in frames:
        optional_initial = {"source_title", "source_identifier_aliases", *INITIAL_FIELDS}
        missing = [column for column in MECHANISM_COLUMNS
                   if column not in optional_initial and column not in frame.columns]
        if missing:
            raise ValueError("mechanism table is missing canonical columns: " + ", ".join(missing))
    for index, frame in enumerate(frames, 1):
        expected_round = str(index)
        rounds = {str(value).strip() for value in frame["proposal_round"]}
        if rounds != {expected_round}:
            raise ValueError(
                f"mechanism table {index} must contain proposal_round={expected_round}; got {sorted(rounds)}"
            )
        interval_states = {_normal(value) for value in frame["interval_informed"]}
        if not interval_states <= {"", "0", "false"}:
            raise ValueError("mechanism coverage audit accepts phenotype-only blind rows; interval-informed input is forbidden")
        import_manifest_path = runtime_path(mechanism_import_manifests_args[index - 1])
        import_manifest = json.loads(import_manifest_path.read_text(encoding="utf-8"))
        if not isinstance(import_manifest, dict) or import_manifest.get("interval_informed") is not False:
            raise ValueError("mechanism-import manifest must explicitly declare interval_informed=false")
        if import_manifest.get("blind_discovery_claim") is not True:
            raise ValueError("mechanism-import manifest must retain the blind discovery claim")
        if str(import_manifest.get("case_id", "")).strip() != str(investigation.get("case_id", "")).strip():
            raise ValueError("mechanism-import manifest case_id does not match the phenotype investigation")
        if str(import_manifest.get("proposal_round")) != expected_round:
            raise ValueError("mechanism-import manifest proposal_round does not match its table position")
        if import_manifest.get("output_sha256") != sha256(runtime_path(mechanisms_args[index - 1])):
            raise ValueError("mechanism-import manifest does not bind the supplied mechanism table")
        localization_binding = import_manifest.get("localization_freeze")
        if not isinstance(localization_binding, dict) or localization_binding.get("sha256") != sha256(
                runtime_path(localization_manifest_arg)):
            raise ValueError("mechanism-import manifest is not bound to the supplied localization freeze")
        # The historical 8-15 portfolio flag is diagnostic metadata.  Coverage
        # evaluates the actual rows and does not make that advisory count a gate.
    source_rows = pd.concat(frames, ignore_index=True).to_dict("records")
    mechanism_leak = _blind_leak_audit(json.dumps(source_rows, ensure_ascii=False, sort_keys=True), candidates_arg)
    if mechanism_leak["leaked_identifiers"]:
        raise ValueError("mechanism portfolio contains QTL interval candidate identifier(s): "
                         + ", ".join(str(value) for value in mechanism_leak["leaked_identifiers"][:5]))
    hypothesis_ids = [str(row["hypothesis_id"]).strip() for row in source_rows]
    duplicate_ids = sorted(value for value, count in Counter(hypothesis_ids).items() if value and count > 1)
    if any(not value for value in hypothesis_ids) or duplicate_ids:
        detail = ", ".join(duplicate_ids) if duplicate_ids else "empty hypothesis_id"
        raise ValueError("mechanism portfolio requires unique non-empty hypothesis_id values: " + detail)
    case_ids = {str(row["case_id"]).strip() for row in source_rows}
    if case_ids != {str(investigation.get("case_id", "")).strip()}:
        raise ValueError("mechanism case_id does not match the phenotype investigation")
    initial_rows: list[dict[str, object]] = []
    initial_families: set[str] = set()
    for row in frames[0].to_dict("records"):
        family = _normal(row["gene_family"])
        if family and family not in initial_families:
            initial_families.add(family)
            initial_rows.append(row)
    initial_pathways = {_normal(row["pathway_or_process"]) for row in initial_rows} - {""}
    if len(initial_families) != len(frames[0]):
        raise ValueError("initial family-node portfolio requires one distinct family per row ("
                         f"target {TARGET_FAMILY_NODE_PORTFOLIO}); duplicate or empty families are not allowed")
    initial_counts = Counter(_normal(row["causal_layer"]) for row in initial_rows)
    initial_requested, initial_reasons = _requested_layers(initial_counts, localization, len(initial_rows))
    initial_bounded_request = initial_requested[:3]
    if len(initial_pathways) < MIN_COMPETING_MECHANISMS:
        initial_reasons.append("distinct_mechanism_breadth_below_advisory_minimum")
    prior_audit: dict[str, object] | None = None
    prior_audit_bytes: bytes | None = None
    prior_audit_hash: str | None = None
    prior_prompt: str | None = None
    if len(frames) == 2:
        if not gap_request_audit_arg:
            raise ValueError("a gap-fill table requires the prior gap-request audit")
        prior_audit_path = runtime_path(gap_request_audit_arg)
        prior_audit_bytes = prior_audit_path.read_bytes()
        prior_audit_hash = hashlib.sha256(prior_audit_bytes).hexdigest()
        prior_audit = json.loads(prior_audit_bytes)
        if not isinstance(prior_audit, dict) or prior_audit.get("status") != "gap_fill_requested":
            raise ValueError("gap-request audit must have status=gap_fill_requested")
        if prior_audit.get("contract_version") != COVERAGE_CONTRACT_VERSION:
            raise ValueError("gap-request audit contract version mismatch")
        if prior_audit.get("initial_mechanisms_sha256") != sha256(runtime_path(mechanisms_args[0])):
            raise ValueError("gap-request audit does not bind the supplied initial mechanisms")
        if prior_audit.get("localization_freeze_sha256") != sha256(runtime_path(localization_manifest_arg)):
            raise ValueError("gap-request audit does not bind the supplied localization freeze")
        prior_prompt = _gap_prompt(investigation, initial_rows, initial_bounded_request)
        if (prior_audit.get("requested_causal_layers") != initial_bounded_request
                or prior_audit.get("gap_prompt_sha256") != hashlib.sha256(prior_prompt.encode("utf-8")).hexdigest()):
            raise ValueError("gap-request audit does not match the deterministic coverage request")
        gap_layers = [_normal(value) for value in frames[1]["causal_layer"]]
        if gap_layers != initial_bounded_request:
            raise ValueError("gap-fill causal_layer rows must exactly match the prior positional request")
    elif gap_request_audit_arg:
        raise ValueError("gap-request audit is only valid when a gap-fill table is supplied")

    rows: list[dict[str, object]] = []
    seen_families: set[str] = set()
    dropped_duplicates: list[str] = []
    for row in source_rows:
        layer = _normal(row["causal_layer"])
        if layer not in CAUSAL_LAYERS:
            raise ValueError(f"invalid or empty causal_layer on {row['hypothesis_id']}: {row['causal_layer']}")
        if str(row["rank_contribution"]).strip() != "none":
            raise ValueError("causal-layer coverage labels must never contribute to ranking")
        family = _normal(row["gene_family"])
        if not family:
            raise ValueError(f"empty gene_family on {row['hypothesis_id']}")
        if family in seen_families:
            dropped_duplicates.append(str(row["hypothesis_id"]))
            continue
        seen_families.add(family)
        row["causal_layer"] = layer
        rows.append(row)
    counts = Counter(str(row["causal_layer"]) for row in rows)
    requested, reasons = _requested_layers(counts, localization, len(rows))
    bounded_request = requested[:3]
    unresolved = requested[len(bounded_request):]
    portfolio_met = bool(rows)
    status = ("pass" if not requested and portfolio_met else
              ("gap_fill_requested" if len(frames) == 1 and bounded_request and portfolio_met else
               "review_required"))
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    preserved_audit_path: Path | None = None
    preserved_prompt_path: Path | None = None
    active_prompt_path = out / "mechanism_gap_fill_prompt.txt"
    if prior_audit_bytes is not None and prior_prompt is not None:
        if active_prompt_path.exists() and active_prompt_path.read_text(encoding="utf-8") != prior_prompt:
            raise ValueError("existing gap-fill prompt does not match the validated prior request")
        preserved_audit_path = out / f"mechanism_gap_request_audit.{prior_audit_hash}.json"
        _preserve_artifact(preserved_audit_path, prior_audit_bytes)
        prompt_bytes = prior_prompt.encode("utf-8")
        preserved_prompt_path = out / f"mechanism_gap_request_prompt.{hashlib.sha256(prompt_bytes).hexdigest()}.txt"
        _preserve_artifact(preserved_prompt_path, prompt_bytes)
        if active_prompt_path.exists():
            active_prompt_path.unlink()
    merged_path = out / "mechanism_hypotheses.tsv"
    pd.DataFrame(rows, columns=MECHANISM_COLUMNS).to_csv(merged_path, sep="\t", index=False, lineterminator="\n")
    pathway_path = out / "llm_pathway_hypotheses.tsv"
    pd.DataFrame(_pathway_rows(rows), columns=HYPOTHESIS_COLUMNS).to_csv(
        pathway_path, sep="\t", index=False, lineterminator="\n")
    audit = {
        "contract_version": COVERAGE_CONTRACT_VERSION, "status": status,
        "initial_mechanisms_sha256": sha256(runtime_path(mechanisms_args[0])),
        "localization_freeze_sha256": sha256(runtime_path(localization_manifest_arg)),
        "localization": localization, "causal_layer_counts": {layer: counts[layer] for layer in CAUSAL_LAYERS},
        "input_row_count": len(source_rows), "merged_row_count": len(rows),
        "dropped_duplicate_hypothesis_ids": dropped_duplicates,
        "reasons": reasons, "requested_causal_layers": bounded_request,
        "unresolved_requested_causal_layers": unresolved,
        "initial_reasons": initial_reasons,
        "gap_request_audit": str(preserved_audit_path) if preserved_audit_path else None,
        "gap_request_audit_sha256": prior_audit_hash,
        "gap_request_prompt": str(preserved_prompt_path) if preserved_prompt_path else None,
        "rank_contribution": "none",
        "initial_hypothesis_target": TARGET_FAMILY_NODE_PORTFOLIO,
        "initial_count_guidance": {"historical_min": MIN_FAMILY_NODE_PORTFOLIO,
                                    "historical_max": MAX_FAMILY_NODE_PORTFOLIO,
                                    "enforced": False},
    }
    audit_path = out / "mechanism_coverage_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    gap_prompt_path: Path | None = None
    gap_prompt_leak: dict[str, object] | None = None
    if status == "gap_fill_requested":
        prompt = _gap_prompt(investigation, rows, bounded_request)
        gap_prompt_leak = _blind_leak_audit(prompt, candidates_arg)
        if gap_prompt_leak["leaked_identifiers"]:
            raise ValueError("gap-fill prompt contains QTL interval candidate identifier(s): "
                             + ", ".join(str(value) for value in gap_prompt_leak["leaked_identifiers"][:5]))
        gap_prompt_path = out / "mechanism_gap_fill_prompt.txt"
        gap_prompt_path.write_text(prompt, encoding="utf-8")
        audit["gap_prompt_sha256"] = hashlib.sha256(prompt.encode("utf-8")).hexdigest()
        audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = {
        "schema_version": "1.0", "candidate_qtl_version": __version__, "created_at_utc": _now(),
        "case_id": investigation.get("case_id", ""), "contract_version": COVERAGE_CONTRACT_VERSION,
        "source_mechanisms": [{"path": str(runtime_path(path)), "sha256": sha256(runtime_path(path))}
                              for path in mechanisms_args],
        "mechanism_import_manifests": [
            {"path": str(runtime_path(path)), "sha256": sha256(runtime_path(path))}
            for path in mechanism_import_manifests_args
        ],
        "localization_freeze_manifest": str(runtime_path(localization_manifest_arg)),
        "localization_freeze_sha256": sha256(runtime_path(localization_manifest_arg)),
        "localization_source_payload_sha256": localization_manifest.get("source_payload_sha256"),
        "interval_blind": True, "mechanism_leak_audit": mechanism_leak,
        "gap_prompt_leak_audit": gap_prompt_leak,
        "rank_contribution": "none", "status": status,
        "coverage_state": "pass" if status == "pass" else "degraded_coverage",
        "freeze_ready": bool(rows),
        "outputs": {"mechanisms": str(merged_path), "pathway_hypotheses": str(pathway_path),
                    "mechanisms_sha256": sha256(merged_path),
                    "pathway_hypotheses_sha256": sha256(pathway_path),
                    "audit": str(audit_path), "audit_sha256": sha256(audit_path),
                    "gap_prompt": str(gap_prompt_path) if gap_prompt_path else None},
    }
    manifest_path = out / "mechanism_coverage_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), "audit": str(audit_path), "status": status,
            "gap_prompt": str(gap_prompt_path) if gap_prompt_path else None,
            "coverage_state": "pass" if status == "pass" else "degraded_coverage",
            "mechanisms": str(merged_path) if rows else None,
            "hypotheses": str(pathway_path) if rows else None,
            "candidate_mechanisms": str(merged_path), "candidate_hypotheses": str(pathway_path)}
