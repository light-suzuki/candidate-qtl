"""Phenotype investigation in front of blind Route A (issue #215).

The monitored blind audit recovered 1/6 known cases, and the reason was upstream
of every ranking rule: Route A's entire phenotype->pathway step was a five-token
static table (``gibberellin``, ``starch``, ``branching``, ``chlorophyll``,
``cotyledon``). Two phenotypes matched nothing at all, three matched a plausible
pathway that no interval annotation carried, and only the mechanism-rich starch
case produced a positive deterministic hit. Improving the rank could not have
fixed any of those: the hypothesis never existed.

The repository already had the pieces -- ``llm-pathway-prompt`` /
``llm-pathway-import`` / ``hypothesis-verify`` / ``bait-panel-freeze`` /
``query-eligibility`` -- but nothing orchestrated them, and the proposal step
asked only for "pathways worth testing". This module adds the missing research
layer and makes it auditable:

``phenotype-investigate``
    Build a phenotype research packet from a structured description (direction,
    inheritance, tissue, stage, onset, penetrance, synonyms/ontology terms,
    stated uncertainty). Emits a prompt package **and** a blank payload template,
    because the point is a schema an LLM and a human fill identically, not an
    LLM feature.

``mechanism-import``
    Normalise either payload into one strict mechanism-hypothesis table:
    a bounded family-node portfolio covering competing mechanisms, with
    model-species mutant/functional evidence, source ids, a proposed non-focal
    accession, a selection rationale, a counter-hypothesis and a discriminating
    observation for every node. Fewer than five
    competing mechanisms, a missing counter-hypothesis, a missing discriminating
    observation, or an unsourced claim is review-required, never silently
    accepted. It also emits the existing ``llm_pathway_hypotheses.tsv`` shape so
    ``hypothesis-verify`` and ``bait-panel-freeze`` consume it unchanged.

``investigation-attribution``
    Say **where** the chain stopped -- phenotype research, verification, bait
    freeze, interval match, or candidate recovery -- so a failure is attributed
    to a stage instead of disappearing into a recall number.

Boundaries, enforced here and covered by tests:

* The generation phase is phenotype-only and interval-blind: no interval, atlas
  or candidate argument reaches the prompt, and when a candidate table is passed
  for auditing, a leaked candidate id aborts the packet.
* An LLM payload and a human payload go through exactly one schema and one gate.
  ``payload_source`` records which, and neither is trusted more.
* Model prose and model confidence never reach a gene score: every row carries
  ``rank_contribution=none``, and promotion still requires the existing
  accession/species/sequence/source verification and the frozen bait panel.
"""

from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from . import __version__
from .core import file_fingerprint, runtime_path, safe_id, sha256
from .claim_evidence import EXTRACTION_COLUMNS
from .llm_pathway import HYPOTHESIS_COLUMNS, _blind_leak_audit, _extract_json
from .phenotype_context import structured_context_block
from .source_verify import classify_citation

SCHEMA_VERSION = "1.0"
INVESTIGATION_PROMPT_VERSION = "phenotype-investigation-1.3.0"
PARSER_VERSION = "mechanism-hypothesis-parser-1.3.0"
MIN_FAMILY_NODE_PORTFOLIO = 8
TARGET_FAMILY_NODE_PORTFOLIO = 12
MAX_FAMILY_NODE_PORTFOLIO = 15
# How the payload was produced. An `agent` (Claude Code / Codex driving the flow
# with its own literature tools) is recorded separately from a raw `llm` API
# completion: they have different provenance and different failure modes, and
# neither is trusted more than a `human` payload. All three pass one gate.
PAYLOAD_SOURCES = ("llm", "agent", "human")
MIN_COMPETING_MECHANISMS = 5
CAUSAL_LAYERS = ("perception", "signaling", "regulation", "execution", "alternative", "unclear")

INVESTIGATION_COLUMNS = [
    "case_id", "focal_species", "phenotype", "direction", "inheritance", "affected_tissue",
    "developmental_stage", "onset", "penetrance", "synonyms", "ontology_terms",
    "known_exclusions", "researcher_constraints", "uncertainty",
]
MECHANISM_COLUMNS = [
    "hypothesis_id", "case_id", "causal_layer", "pathway_or_process", "gene_family", "mechanism_statement",
    "expected_role", "direction_consistency", "model_species_evidence", "model_species_mutant_phenotype",
    "source_ids", "proposed_stable_id", "reference_species", "reference_gene_name", "selection_rationale",
    "counter_hypothesis", "discriminating_observation", "confidence", "uncertainty",
    "payload_source", "model", "prompt_version", "parser_version", "retrieval_provenance",
    "proposal_round", "interval_informed", "review_status", "rank_contribution", "notes",
]
REVIEW_COLUMNS = ["reason", "target", "required_next_action", "allowed_choices"]
from .relevance_tiers import INITIAL_FIELDS
MECHANISM_COLUMNS += ["source_title", "source_identifier_aliases", *INITIAL_FIELDS]
ATTRIBUTION_STAGES = ("phenotype_research", "verification", "bait_freeze", "interval_match", "candidate_recovery")

SYSTEM_PROMPT = """You are researching a plant QTL phenotype BEFORE anyone looks at the interval.

You are NOT ranking candidate genes, NOT naming a causal gene, and you will not be
shown the interval. Your job is to produce several COMPETING mechanism hypotheses
that could explain this phenotype, each one testable and sourced.

For every hypothesis give:
- one causal-layer label chosen only from perception, signaling, regulation,
  execution, alternative, or unclear;
- the pathway or process, and the gene family that carries it;
- a mechanism statement: how a lesion in that family produces THIS phenotype, in
  THIS direction, with THIS inheritance;
- model-species evidence: a named mutant or functional study in a species other
  than the focal species, with a source id (PMID/DOI/UniProt/TAIR accession);
- one stable database accession (UniProt preferred) for a reference gene in a
  non-focal species, plus why you selected that gene rather than a paralogue;
- a counter-hypothesis: the strongest reason this mechanism is wrong;
- a discriminating observation: what you would measure to tell this hypothesis
  apart from the others you listed.

Rules:
- Return a target of 12 initial family-level hypotheses while remaining
  interval-blind. Use the target when useful, but do not pad or reject a
  well-formed portfolio solely because its row count differs. It should contain
  at least five competing mechanisms spanning
  causal-level diversity. Cover the applicable levels among sensing/timing,
  synthesis/catabolism, signal integration/regulation, downstream effector/transport,
  and a developmental/tissue alternative. If one level is
  not biologically applicable, say so explicitly and replace it with an
  independent causal alternative at another level.
- For each major plausible mechanism, include 2-3 functionally non-redundant
  gene-family nodes when the biology supports them (for example upstream input,
  core regulator, and downstream output). Emit each node as its own row with its
  own reference gene and evidence. Do not collapse a multi-protein pathway to one
  convenient representative.
- Do not satisfy this requirement by swapping gene names, paralogues, or
  reference accessions within the same family. Such swaps add no family coverage
  and do not count as distinct mechanisms.
- Do not invent accessions or source ids. Leave the field empty instead. An empty
  field is recoverable; a fabricated one is not.
- Everything you return is a proposal. Every accession, species, sequence and
  function is re-verified against an authoritative database before use, and your
  confidence never contributes to any gene score."""

MECHANISM_PAYLOAD_TEMPLATE = """{"hypotheses": [{"causal_layer": "perception|signaling|regulation|execution|alternative|unclear",
  "pathway_or_process": "", "gene_family": "", "mechanism_statement": "",
  "expected_role": "", "direction_consistency": "", "model_species_evidence": "",
  "model_species_mutant_phenotype": "", "source_ids": "", "proposed_stable_id": "",
  "reference_species": "", "reference_gene_name": "", "selection_rationale": "",
  "counter_hypothesis": "", "discriminating_observation": "",
  "confidence": "low|medium|high", "uncertainty": ""}]}"""

USER_TEMPLATE = """Case: {case_id}
Focal species: {focal_species}
Phenotype: {phenotype}
Direction: {direction}
Inheritance: {inheritance}
Affected tissue/organ: {tissue}
Developmental stage: {stage}
Onset: {onset}
Penetrance/expressivity: {penetrance}
Phenotype synonyms: {synonyms}
Ontology terms: {ontology_terms}
Known exclusions: {exclusions}
Researcher constraints: {constraints}
Stated uncertainty: {uncertainty}

Return JSON only:
""" + MECHANISM_PAYLOAD_TEMPLATE.replace("{", "{{").replace("}", "}}")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _frame(rows: list[dict[str, object]], columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame(rows)[columns] if rows else pd.DataFrame(columns=columns)


def _write(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


# --------------------------------------------------------------------------- #
# 1. investigation packet
# --------------------------------------------------------------------------- #

def phenotype_investigate(
    out_dir_arg: str,
    case_id: str,
    focal_species: str,
    phenotype: str,
    direction: str = "unknown",
    inheritance: str = "unknown",
    tissue: str = "",
    stage: str = "",
    onset: str = "",
    penetrance: str = "",
    synonyms: list[str] | None = None,
    ontology_terms: list[str] | None = None,
    exclusions: list[str] | None = None,
    constraints: list[str] | None = None,
    uncertainty: str = "",
    candidates_arg: str | None = None,
    phenotype_context: str | None = None,
) -> dict[str, object]:
    """Build the interval-blind phenotype research packet.

    ``candidates_arg`` is used only as a leak check: the interval candidate ids
    must not appear anywhere in the packet, and are never inserted into it.

    ``phenotype_context`` (optional ``phenotype_context.yaml``) is recorded in
    the manifest for audit only.  It never replaces the free text and no
    structured field is inserted into the investigation prompt.
    """
    if not case_id or not focal_species or not phenotype:
        raise ValueError("case_id, focal_species and phenotype are required")
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)

    record = {
        "case_id": case_id, "focal_species": focal_species, "phenotype": phenotype,
        "direction": direction or "unknown", "inheritance": inheritance or "unknown",
        "affected_tissue": tissue, "developmental_stage": stage, "onset": onset, "penetrance": penetrance,
        "synonyms": ";".join(synonyms or []), "ontology_terms": ";".join(ontology_terms or []),
        "known_exclusions": ";".join(exclusions or []), "researcher_constraints": ";".join(constraints or []),
        "uncertainty": uncertainty,
    }
    investigation_path = out / "phenotype_investigation.tsv"
    _write(_frame([record], INVESTIGATION_COLUMNS), investigation_path)

    user = USER_TEMPLATE.format(
        case_id=case_id, focal_species=focal_species, phenotype=phenotype,
        direction=record["direction"], inheritance=record["inheritance"],
        tissue=tissue or "unspecified", stage=stage or "unspecified", onset=onset or "unspecified",
        penetrance=penetrance or "unspecified",
        synonyms="; ".join(synonyms or []) or "none",
        ontology_terms="; ".join(ontology_terms or []) or "none",
        exclusions="; ".join(exclusions or []) or "none",
        constraints="; ".join(constraints or []) or "none",
        uncertainty=uncertainty or "not stated",
    )
    prompt_text = f"{SYSTEM_PROMPT}\n\n{user}\n"
    prompt_path = out / "investigation_prompt.txt"
    prompt_path.write_text(prompt_text, encoding="utf-8")
    prompts_path = out / "investigation_prompts.jsonl"
    prompts_path.write_text(json.dumps(
        {"prompt_id": f"investigation::{safe_id(case_id)}", "gene_id": f"phenotype_investigation:{case_id}",
         "system": SYSTEM_PROMPT, "user": user}, ensure_ascii=False) + "\n", encoding="utf-8")

    # The same schema a model must produce, empty, so a human can author the
    # payload by hand and import it through the identical gate.
    template_path = out / "mechanism_payload_template.tsv"
    _write(pd.DataFrame(columns=[column for column in MECHANISM_COLUMNS
                                 if column not in {"hypothesis_id", "parser_version", "review_status",
                                                   "rank_contribution", "notes"}]), template_path)

    leak = _blind_leak_audit(prompt_text, candidates_arg)
    if leak["leaked_identifiers"]:
        raise ValueError("phenotype investigation packet contains QTL interval candidate identifier(s): "
                         + ", ".join(str(value) for value in leak["leaked_identifiers"][:5]))

    manifest = {
        "schema_version": SCHEMA_VERSION, "case_id": case_id, "candidate_qtl_version": __version__,
        "created_at_utc": _now(),
        "stage": "phenotype_research",
        "prompt_version": INVESTIGATION_PROMPT_VERSION,
        "prompt_sha256": hashlib.sha256(prompt_text.encode("utf-8")).hexdigest(),
        "investigation": record,
        "structured_context": structured_context_block(
            phenotype_context, case_id=case_id, phenotype_text=phenotype,
            query_usage={
                "usage": "stored_only",
                "reason": ("the investigation prompt is built from the free-text phenotype exactly as before; "
                           "no structured field is inserted into the prompt"),
            },
        ),
        "interval_blind": {
            "interval_arguments": [],
            "qtl_candidate_information_in_packet": False,
            "leak_audit": leak,
            "rule": "the research packet is built from the phenotype description only; no interval, atlas or candidate content enters it",
        },
        "payload_contract": {
            "min_competing_mechanisms": MIN_COMPETING_MECHANISMS,
            "initial_hypothesis_target": TARGET_FAMILY_NODE_PORTFOLIO,
            "initial_count_guidance": {"historical_min": MIN_FAMILY_NODE_PORTFOLIO,
                                        "historical_max": MAX_FAMILY_NODE_PORTFOLIO,
                                        "enforced": False},
            "initial_hypothesis_count": TARGET_FAMILY_NODE_PORTFOLIO,
            "causal_layers": list(CAUSAL_LAYERS),
            "bounded_gap_fill_max": 3,
            "accepted_payload_sources": list(PAYLOAD_SOURCES),
            "template": str(template_path),
            "rule": "an LLM payload and a human payload are imported through one schema and one gate; neither is trusted more",
        },
        "external_transmission": "none: this command writes a prompt package; executing it is a separate, recorded step",
        "outputs": {"investigation": str(investigation_path), "prompt": str(prompt_path),
                    "prompts_jsonl": str(prompts_path), "payload_template": str(template_path)},
        "next_step": "run the prompt (or author the template by hand), then `mechanism-import`",
    }
    manifest_path = out / "phenotype_investigation_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), "investigation": str(investigation_path), "prompt": str(prompt_path),
            "prompts_jsonl": str(prompts_path), "payload_template": str(template_path),
            "prompt_version": INVESTIGATION_PROMPT_VERSION, "interval_blind": True}


# --------------------------------------------------------------------------- #
# 2. mechanism payload import (LLM or human, one schema)
# --------------------------------------------------------------------------- #

def _payload_records(path: Path, detected_schema: str | None = None) -> list[dict[str, object]]:
    """Accept a TSV (human), a JSON object/array, or JSONL (model output)."""
    if path.suffix.casefold() == ".tsv":
        return pd.read_csv(path, sep="\t", dtype=str).fillna("").to_dict("records")
    text = path.read_text(encoding="utf-8")
    payloads: list[object] = []
    if detected_schema == "jsonl" or path.suffix.casefold() == ".jsonl":
        payloads = [json.loads(line) for line in text.splitlines() if line.strip()]
    else:
        parsed = _extract_json(text)
        if parsed is not None:
            payloads.append(parsed)
    records: list[dict[str, object]] = []
    for payload in payloads:
        if isinstance(payload, dict):
            for key in ("hypotheses", "mechanisms", "proposals", "results"):
                value = payload.get(key)
                if isinstance(value, list):
                    records.extend(item for item in value if isinstance(item, dict))
                    break
            else:
                nested_found = False
                for key in ("response", "parsed", "output", "text", "content", "output_text"):
                    value = payload.get(key)
                    if isinstance(value, str):
                        nested = _extract_json(value)
                        if isinstance(nested, dict):
                            records.extend(item for item in nested.get("hypotheses", []) if isinstance(item, dict))
                            nested_found = True
                    elif isinstance(value, (dict, list)):
                        records.extend(item for item in (value if isinstance(value, list) else [value]) if isinstance(item, dict))
                        nested_found = True
                if not nested_found:
                    records.append(payload)
        elif isinstance(payload, list):
            records.extend(item for item in payload if isinstance(item, dict))
    return records


def _unrecognised_citation_tokens(source_ids: str) -> list[str]:
    """source_ids tokens that source-verify would reject as unrecognised.

    Splits with the same separator expression source-verify uses, so the early
    gate and the downstream rejection can never disagree about a token.
    """
    return [token for token in (value.strip() for value in re.split(r"[;,|]+", str(source_ids or "")) if value.strip())
            if classify_citation(token)[0] == "unrecognised"]


def mechanism_import(
    payload_arg: str,
    out_dir_arg: str,
    investigation_arg: str | None = None,
    payload_source: str = "llm",
    model: str = "",
    case_id: str = "",
    retrieval_provenance: str = "",
    candidates_arg: str | None = None,
    min_mechanisms: int = MIN_COMPETING_MECHANISMS,
    proposal_round: int = 1,
    interval_informed: bool = False,
    mechanism_freeze_manifest: str | None = None,
    localization_freeze_manifest: str | None = None,
) -> dict[str, object]:
    """Normalise an LLM or human payload into the mechanism-hypothesis schema.

    ``proposal_round`` and ``interval_informed`` are independent. A bounded
    phenotype-only complement may be round 2 and remain blind. Any payload that
    saw an interval verdict is explicitly marked and must never be reported as
    blind discovery.

    Payloads are schema-validated fail-closed before anything is written:
    an unknown schema, duplicate or forbidden columns, a legacy row without
    case_id/mechanism_id, a duplicate mechanism_id, a payload whose bytes do
    not match ``mechanism_freeze_manifest`` (coordination/payload_freeze.tsv
    format) raises ValueError. Canonical rows with unresolvable citation formats
    are quarantined with their original records; other rows remain proposals. A
    frozen legacy 7-column payload is converted
    through the explicit adapter, whose out-of-band fields (mechanism_id,
    evidence_class, status) land only in the adapter map, never in the
    canonical rows.
    """
    if payload_source not in PAYLOAD_SOURCES:
        raise ValueError(f"payload_source must be one of: {', '.join(PAYLOAD_SOURCES)}")
    if proposal_round < 1:
        raise ValueError("proposal round must be 1 or greater")
    # Lazy import: mechanism_payload_adapter imports this module at module
    # level, so a top-level import here would be a cycle.
    from .mechanism_payload_adapter import (ADAPTER_VERSION, OUTPUT_SCHEMA_VERSION, convert_legacy_v1_7col,
                                            normalize_source_ids, validate_payload_file, write_adapter_map)
    payload_path = runtime_path(payload_arg)
    out = runtime_path(out_dir_arg)

    investigation: dict[str, object] = {}
    if investigation_arg:
        frame = pd.read_csv(runtime_path(investigation_arg), sep="\t", dtype=str).fillna("")
        if frame.empty:
            raise ValueError("phenotype investigation table is empty")
        investigation = frame.to_dict("records")[0]
    investigation_case_id = str(investigation.get("case_id", ""))

    localization_provenance: dict[str, object] | None = None
    if localization_freeze_manifest:
        localization_manifest_path = runtime_path(localization_freeze_manifest)
        localization_manifest = json.loads(localization_manifest_path.read_text(encoding="utf-8"))
        if not isinstance(localization_manifest, dict) or localization_manifest.get("interval_informed") is not False:
            raise ValueError("localization freeze manifest must explicitly declare interval_informed=false")
        localization_case_id = str(localization_manifest.get("case_id", "")).strip()
        expected_case_id = str(case_id or investigation_case_id).strip()
        if expected_case_id and localization_case_id != expected_case_id:
            raise ValueError("localization freeze manifest case_id does not match mechanism import")
        localization_provenance = {
            "path": str(localization_manifest_path),
            "sha256": sha256(localization_manifest_path),
            "frozen_at_utc": localization_manifest.get("frozen_at_utc"),
            "external_generation_order_verified": False,
            "verification_boundary": "import binds to an existing freeze; external generation order requires executor provenance",
        }

    # Fail-closed before anything is written: a rejected schema, legacy row or
    # freeze mismatch raises ValueError here, never a silent 0-row import.
    validation = validate_payload_file(payload_path, mechanism_freeze_manifest,
                                       case_id or investigation_case_id)
    resolved_case_id = case_id or investigation_case_id or str(validation.get("payload_case_id", ""))
    out.mkdir(parents=True, exist_ok=True)
    records = _payload_records(payload_path, str(validation["schema"]))
    adapter_version = "none"
    adapter_map_rows: list[dict[str, object]] = []
    if validation["schema"] == "legacy_v1_7col":
        canonical_rows, adapter_map_rows = convert_legacy_v1_7col(
            records, case_id=resolved_case_id,
            payload_source=payload_source, proposal_round=proposal_round,
            interval_informed=interval_informed)
        records = canonical_rows
        adapter_version = ADAPTER_VERSION

    raw_copy = out / f"raw_payload{payload_path.suffix or '.txt'}"
    raw_copy.write_bytes(payload_path.read_bytes())
    copied_digest = sha256(raw_copy)
    if copied_digest != validation["raw_payload_sha256"]:
        raise ValueError(
            f"mechanism payload freeze mismatch: {payload_path}; payload changed between validation and copy "
            f"(validated sha256 {validation['raw_payload_sha256']}; copied sha256 {copied_digest})")
    adapter_map_path: Path | None = None
    if validation["schema"] == "legacy_v1_7col":
        adapter_map_path = write_adapter_map(adapter_map_rows, out)

    # Agent-extracted claim spans ride in the payload as a top-level
    # "claim_extractions" list. They are written out verbatim for source-verify's
    # deterministic gate; the gate, not this import, decides what they mean.
    claim_extraction_rows: list[dict[str, object]] = []
    if validation["schema"] in ("json", "jsonl"):
        parsed_top = _extract_json(payload_path.read_text(encoding="utf-8"))
        if isinstance(parsed_top, list) and parsed_top and isinstance(parsed_top[0], dict):
            parsed_top = parsed_top[0]
        if isinstance(parsed_top, dict) and isinstance(parsed_top.get("claim_extractions"), list):
            claim_extraction_rows = [item for item in parsed_top["claim_extractions"] if isinstance(item, dict)]

    case_id = resolved_case_id
    # A legacy frozen payload keeps its historical `file#pointer` tokens, so the
    # machine-verifiable-citation gate applies only to canonical schemas.
    citation_gate_exempt = validation["schema"] == "legacy_v1_7col"
    rejected_rows: list[dict[str, object]] = []
    review_rows: list[dict[str, str]] = []
    rows: list[dict[str, object]] = []
    for index, record in enumerate(records, 1):
        generated_id = (f"mech_{index:03d}" if proposal_round == 1
                        else f"mech_r{proposal_round}_{index:03d}")
        hypothesis_id = str(record.get("hypothesis_id", "") or generated_id)
        values: dict[str, str] = {}
        offending_tokens: list[str] = []
        for column in MECHANISM_COLUMNS:
            raw_value = record.get(column, "")
            if column == "source_ids":
                normalized_source_ids = normalize_source_ids(
                    raw_value, location=f"hypotheses[{index}].source_ids"
                )
                values[column] = (
                    normalized_source_ids
                    if isinstance(raw_value, list)
                    else normalized_source_ids.strip()
                )
                if not citation_gate_exempt and values[column]:
                    offending_tokens.extend(
                        f"{hypothesis_id}:'{token}'"
                        for token in _unrecognised_citation_tokens(values[column]))
            else:
                values[column] = "" if raw_value is None else str(raw_value).strip()
        from .paper_reports import imported_paper_fields
        values.update(imported_paper_fields(record))
        pathway = values["pathway_or_process"]
        family = values["gene_family"]
        causal_layer = values["causal_layer"].casefold()
        if causal_layer and causal_layer not in CAUSAL_LAYERS:
            raise ValueError(
                f"mechanism payload rejected: {hypothesis_id} causal_layer '{values['causal_layer']}' is invalid; "
                f"choose one of: {', '.join(CAUSAL_LAYERS)}"
            )
        values["causal_layer"] = causal_layer
        if offending_tokens:
            rejected_rows.append({"row_index": index, "hypothesis_id": hypothesis_id,
                                  "reason": "citation_is_not_a_recognisable_pmid_pmcid_or_doi",
                                  "tokens": offending_tokens, "raw_record": record})
            review_rows.append({"reason": "unrecognised_source_id_quarantined", "target": hypothesis_id,
                                "required_next_action": "verify the original source identifiers; the row remains not_evaluated",
                                "allowed_choices": "supply independently verified identifiers | leave not_evaluated"})
            continue
        missing = [name for name in ("pathway_or_process", "gene_family", "mechanism_statement") if not values[name]]
        if missing:
            review_rows.append({"reason": "mechanism_missing_required_field", "target": hypothesis_id,
                                "required_next_action": f"supply {', '.join(missing)} or discard the hypothesis",
                                "allowed_choices": "re-run the prompt | complete the payload manually | discard"})
        for name, reason in (("counter_hypothesis", "mechanism_without_counter_hypothesis"),
                             ("discriminating_observation", "mechanism_without_discriminating_observation"),
                             ("source_ids", "mechanism_without_source_id"),
                             ("model_species_evidence", "mechanism_without_model_species_evidence")):
            if not values[name]:
                review_rows.append({"reason": reason, "target": hypothesis_id,
                                    "required_next_action": f"a hypothesis without {name.replace('_', ' ')} cannot be treated as researched",
                                    "allowed_choices": "supply it | keep as an unverifiable proposal | discard"})
        if not values["proposed_stable_id"]:
            review_rows.append({"reason": "mechanism_without_stable_accession", "target": hypothesis_id,
                                "required_next_action": "supply a stable accession for the reference gene, or keep the hypothesis unverifiable",
                                "allowed_choices": "add accession | keep as unverifiable proposal"})
        rows.append({
            **values,
            "hypothesis_id": hypothesis_id, "case_id": case_id or values["case_id"],
            "pathway_or_process": pathway, "gene_family": family,
            "confidence": values["confidence"] or "unknown",
            "payload_source": payload_source, "model": model,
            "prompt_version": str(investigation.get("prompt_version", "") or INVESTIGATION_PROMPT_VERSION),
            "parser_version": PARSER_VERSION,
            "retrieval_provenance": values["retrieval_provenance"] or retrieval_provenance,
            "proposal_round": proposal_round, "interval_informed": bool(interval_informed),
            "review_status": "pending_verification", "rank_contribution": "none",
            "notes": ("unknown: " + ", ".join(missing)) if missing else "",
        })

    distinct_mechanisms = len({str(row["pathway_or_process"]).casefold() for row in rows if row["pathway_or_process"]})
    distinct_families = len({str(row["gene_family"]).casefold() for row in rows if row["gene_family"]})
    family_counts: dict[str, int] = {}
    for row in rows:
        family = str(row["gene_family"]).strip().casefold()
        if family:
            family_counts[family] = family_counts.get(family, 0) + 1
    duplicate_families = sorted(family for family, count in family_counts.items() if count > 1)
    if not rows:
        review_rows.append({"reason": "no_parsable_mechanism", "target": str(raw_copy),
                            "required_next_action": "the payload produced no schema-conformant hypothesis; re-run or author it manually",
                            "allowed_choices": "re-run the prompt | author manually | abstain"})
    if interval_informed:
        review_rows.append({"reason": "interval_informed_proposal_round", "target": case_id or "case",
                            "required_next_action": (f"round {proposal_round} was authored after an interval verdict; "
                                                     "keep it in its own lane and never report it as blind discovery"),
                            "allowed_choices": "別laneとして保持 | round 1 のみで評価 | 破棄"})
    for family in duplicate_families:
        review_rows.append({"reason": "duplicate_family_node", "target": family,
                            "required_next_action": ("collapse paralogue/accession swaps within this family or replace them "
                                                     "with a functionally non-redundant family node"),
                            "allowed_choices": "collapse duplicate family rows | replace with another family node | abstain"})
    family_node_portfolio_met = bool(rows and distinct_families == len(rows))
    leak = _blind_leak_audit(json.dumps([*rows, *rejected_rows], ensure_ascii=False), candidates_arg)
    if leak["leaked_identifiers"]:
        raise ValueError("mechanism payload names QTL interval candidate identifier(s), so it was not interval-blind: "
                         + ", ".join(str(value) for value in leak["leaked_identifiers"][:5]))

    from .paper_reports import stamp_paper_handoff
    rows = stamp_paper_handoff(rows, route="phenotype_first")
    optional_initial = {"source_title", "source_identifier_aliases", *INITIAL_FIELDS}
    mechanisms = _frame(rows, [column for column in MECHANISM_COLUMNS
        if column not in optional_initial or any(row.get(column, "") for row in rows)])
    quarantine_path = out / "mechanism_rejected_rows.json"
    quarantine_path.write_text(json.dumps(rejected_rows, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    mechanisms_path = out / "mechanism_hypotheses.tsv"
    _write(mechanisms, mechanisms_path)

    # The existing verification chain (hypothesis-verify -> bait-panel-freeze)
    # consumes the llm_pathway hypothesis shape; emit it rather than forking the
    # downstream commands.
    pathway_rows = [{
        **{key: row.get(key, "") for key in ("source_title", "source_identifier_aliases", *INITIAL_FIELDS)},
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
    pathway_path = out / "llm_pathway_hypotheses.tsv"
    _write(_frame(pathway_rows, HYPOTHESIS_COLUMNS), pathway_path)

    reviews = _frame(review_rows, REVIEW_COLUMNS).drop_duplicates().sort_values(["reason", "target"], kind="stable")
    reviews_path = out / "review_required.tsv"
    _write(reviews, reviews_path)

    claim_extractions_path: Path | None = None
    if claim_extraction_rows:
        claim_extractions_path = out / "claim_extractions.tsv"
        extraction_columns = list(EXTRACTION_COLUMNS)
        if any(str(row.get("hypothesis_id", "")).strip() for row in claim_extraction_rows):
            extraction_columns.append("hypothesis_id")
        pd.DataFrame(claim_extraction_rows, columns=extraction_columns).fillna("").to_csv(
            claim_extractions_path, sep="\t", index=False, lineterminator="\n")

    from .paper_evidence_table import write_paper_evidence_table
    paper_evidence = write_paper_evidence_table(
        rows, out, route="phenotype_first", source_artifact=str(raw_copy),
        audit_focal_species=str(investigation.get("focal_species", "")))
    manifest = {
        "paper_evidence": paper_evidence,
        "schema_version": SCHEMA_VERSION, "case_id": case_id, "candidate_qtl_version": __version__,
        "stage": "phenotype_research", "imported_at_utc": _now(),
        "payload_source": payload_source, "model": model,
        "proposal_round": proposal_round,
        "interval_informed": bool(interval_informed),
        "localization_freeze": localization_provenance,
        "blind_discovery_claim": not bool(interval_informed),
        "round_rule": (
            "proposal round and interval exposure are independent: a bounded phenotype-only complement may be a "
            "blind later round; only an explicitly interval-informed payload loses the blind discovery claim"
        ),
        "parser_version": PARSER_VERSION, "prompt_version": str(investigation.get("prompt_version", INVESTIGATION_PROMPT_VERSION)),
        "raw_payload": str(raw_copy), "raw_payload_sha256": sha256(raw_copy),
        "investigation": str(runtime_path(investigation_arg)) if investigation_arg else None,
        "accepted_payload_schema": validation["schema"],
        "adapter_version": adapter_version,
        "input_schema_version": validation["schema"],
        "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "rows_read": validation["rows_read"],
        "rows_imported": len(rows),
        "rows_rejected": len(rejected_rows),
        "import_state": "not_evaluated_all_rows_rejected" if rejected_rows and not rows else (
            "partial_import" if rejected_rows else "imported"),
        "rejected_rows": str(quarantine_path),
        "rejected_rows_sha256": sha256(quarantine_path),
        "converted_at_utc": _now(),
        "output_sha256": sha256(mechanisms_path),
        "freeze_verification": validation["freeze_verification"],
        "adapter_map": ({"path": str(adapter_map_path), "sha256": sha256(adapter_map_path)}
                        if adapter_map_path else None),
        "hypothesis_count": len(rows),
        "distinct_mechanisms": distinct_mechanisms,
        "distinct_gene_families": distinct_families,
        "duplicate_gene_families": duplicate_families,
        "causal_layer_counts": {
            layer: sum(1 for row in rows if row["causal_layer"] == layer) for layer in CAUSAL_LAYERS
        },
        "target_family_node_portfolio": TARGET_FAMILY_NODE_PORTFOLIO,
        "initial_count_guidance": {"historical_min": MIN_FAMILY_NODE_PORTFOLIO,
                                    "historical_max": MAX_FAMILY_NODE_PORTFOLIO,
                                    "enforced": False},
        "family_node_portfolio_met": family_node_portfolio_met,
        "min_competing_mechanisms": min_mechanisms,
        "competing_mechanisms_met": bool(distinct_mechanisms >= min_mechanisms),
        "sourced_hypotheses": sum(1 for row in rows if row["source_ids"]),
        "hypotheses_with_counter_hypothesis": sum(1 for row in rows if row["counter_hypothesis"]),
        "hypotheses_with_discriminating_observation": sum(1 for row in rows if row["discriminating_observation"]),
        "hypotheses_with_accession": sum(1 for row in rows if row["proposed_stable_id"]),
        "review_required_count": int(len(reviews)),
        "interval_blind": {"leak_audit": leak, "interval_arguments": [],
                           "rule": "a payload naming an interval candidate id is rejected: the research phase must be blind"},
        "scoring_boundary": (
            "mechanism prose and model confidence never enter a gene score; rank_contribution is none on every row, "
            "and promotion still requires hypothesis-verify plus a frozen bait panel"
        ),
        "outputs": {**paper_evidence["outputs"], "mechanisms": str(mechanisms_path), "pathway_hypotheses": str(pathway_path),
                    "review_required": str(reviews_path), "raw_payload": str(raw_copy)},
        "claim_extractions": {"supplied": bool(claim_extractions_path), "row_count": len(claim_extraction_rows),
                              "path": str(claim_extractions_path) if claim_extractions_path else None},
        "next_step": (paper_evidence["next_action"] if paper_evidence["recording_state"] == "needs_completion" else
                      "hypothesis-verify --hypotheses " + str(pathway_path)),
    }
    manifest_path = out / "mechanism_import_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"paper_evidence": paper_evidence, **paper_evidence["outputs"],
            "manifest": str(manifest_path), "mechanisms": str(mechanisms_path),
            "hypotheses": str(pathway_path), "review_required": str(reviews_path),
            "hypothesis_count": len(rows), "distinct_mechanisms": distinct_mechanisms,
            "rows_rejected": len(rejected_rows), "import_state": manifest["import_state"],
            "rejected_rows": str(quarantine_path),
            "competing_mechanisms_met": bool(distinct_mechanisms >= min_mechanisms),
            "claim_extractions": str(claim_extractions_path) if claim_extractions_path else None}


# --------------------------------------------------------------------------- #
# 3. stage-by-stage failure attribution
# --------------------------------------------------------------------------- #

def _count(path_arg: str | None, predicate=None) -> int:
    if not path_arg:
        return 0
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    return int(len(frame if predicate is None else frame[frame.apply(predicate, axis=1)]))


def investigation_attribution(
    out_dir_arg: str,
    mechanisms_arg: str | None = None,
    validation_arg: str | None = None,
    panel_arg: str | None = None,
    interval_match_arg: str | None = None,
    ranked_arg: str | None = None,
    case_id: str = "",
    min_mechanisms: int = MIN_COMPETING_MECHANISMS,
) -> dict[str, object]:
    """Attribute a Route A outcome to the stage where it actually stopped.

    Truth-free by construction: it reports how far the chain got, not whether the
    answer was right. That keeps it usable inside a prediction-only run.
    """
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)

    mechanisms = pd.read_csv(runtime_path(mechanisms_arg), sep="\t", dtype=str).fillna("") if mechanisms_arg else pd.DataFrame()
    distinct = len({str(value).casefold() for value in mechanisms.get("pathway_or_process", pd.Series(dtype=str)) if str(value).strip()})
    sourced = int(mechanisms.get("source_ids", pd.Series(dtype=str)).astype(str).str.strip().ne("").sum()) if len(mechanisms) else 0
    verified = _count(validation_arg, lambda row: str(row.get("verification_status", "")) == "verified")
    rejected = _count(validation_arg, lambda row: str(row.get("verification_status", "")) == "rejected")
    baits = _count(panel_arg)
    matched_genes: set[str] = set()
    if interval_match_arg:
        match = pd.read_csv(runtime_path(interval_match_arg), sep="\t", dtype=str).fillna("")
        column = "matched_interval_gene" if "matched_interval_gene" in match.columns else "member_interval_genes"
        matched_genes = {gene for value in match.get(column, pd.Series(dtype=str))
                         for gene in str(value).split(";") if gene.strip()}
    recovered = 0
    if ranked_arg and matched_genes:
        ranked = pd.read_csv(runtime_path(ranked_arg), sep="\t", dtype=str).fillna("")
        recovered = int(sum(1 for value in ranked.get("gene_id", pd.Series(dtype=str)) if str(value) in matched_genes))

    stages = {
        "phenotype_research": {
            "status": "complete" if distinct >= min_mechanisms and sourced else
                      ("partial" if len(mechanisms) else "not_reached"),
            "distinct_mechanisms": distinct, "required": min_mechanisms,
            "sourced_hypotheses": sourced, "hypothesis_count": int(len(mechanisms)),
            "reason": (
                "fewer than the required competing mechanisms" if distinct < min_mechanisms else
                ("no hypothesis carries a source id" if not sourced else "")
            ),
        },
        "verification": {
            "status": "complete" if verified else ("attempted" if validation_arg else "not_reached"),
            "verified": verified, "rejected": rejected,
            "reason": "" if verified else ("every proposal was rejected or unverifiable" if validation_arg else "verification was not run"),
        },
        "bait_freeze": {
            "status": "complete" if baits else ("attempted" if panel_arg else "not_reached"),
            "frozen_baits": baits,
            "reason": "" if baits else ("no verified proposal survived to the panel" if panel_arg else "the panel was not frozen"),
        },
        "interval_match": {
            "status": "complete" if matched_genes else ("attempted" if interval_match_arg else "not_reached"),
            "matched_interval_genes": len(matched_genes),
            "reason": "" if matched_genes else ("no bait reached an interval gene" if interval_match_arg else "interval matching was not run"),
        },
        "candidate_recovery": {
            "status": "complete" if recovered else ("attempted" if ranked_arg else "not_reached"),
            "recovered_ranked_genes": recovered,
            "reason": "" if recovered else ("matched interval genes did not appear in the ranked table" if ranked_arg else "no ranked table supplied"),
        },
    }
    first_incomplete = next((stage for stage in ATTRIBUTION_STAGES if stages[stage]["status"] != "complete"), "")
    manifest = {
        "schema_version": SCHEMA_VERSION, "case_id": case_id, "candidate_qtl_version": __version__,
        "attributed_at_utc": _now(), "stage_order": list(ATTRIBUTION_STAGES), "stages": stages,
        "first_incomplete_stage": first_incomplete or "none",
        "outcome": "chain_complete" if not first_incomplete else f"stopped_at:{first_incomplete}",
        "truth_used": False,
        "interpretation": [
            "this attribution says where the chain stopped, not whether the answer was right",
            "a miss at phenotype_research cannot be repaired by any ranking change downstream",
            "a not_reached stage is an unavailable stage, not a failed one",
        ],
        "inputs": {name: (file_fingerprint(runtime_path(value)) if value else None) for name, value in (
            ("mechanisms", mechanisms_arg), ("validation", validation_arg), ("panel", panel_arg),
            ("interval_match", interval_match_arg), ("ranked", ranked_arg))},
    }
    manifest_path = out / "investigation_attribution.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    table_path = out / "investigation_attribution.tsv"
    _write(pd.DataFrame([{"stage": stage, "status": stages[stage]["status"], "reason": stages[stage]["reason"]}
                         for stage in ATTRIBUTION_STAGES], columns=["stage", "status", "reason"]), table_path)
    return {"manifest": str(manifest_path), "attribution": str(table_path),
            "first_incomplete_stage": manifest["first_incomplete_stage"], "outcome": manifest["outcome"]}
