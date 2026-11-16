"""Route-A phenotype investigation with an explicit reference-target contract.

The established implementation remains in ``phenotype_investigation_legacy``.
This facade strengthens the prompt and import manifest for issue #654 without
inventing identifiers or weakening downstream source/provider verification.
"""

from __future__ import annotations

import json as _json
import hashlib as _hashlib
import re as _re
from pathlib import Path as _Path

import pandas as _pd

from . import phenotype_investigation_legacy as _legacy
from .core import runtime_path as _runtime_path
from .core import sha256 as _sha256
from .paper_reports import LITERATURE_RESEARCH_POLICY, LITERATURE_SEARCH_GUIDANCE, render_initial_paper_report_example

for _name in dir(_legacy):
    if not _name.startswith("__"):
        globals()[_name] = getattr(_legacy, _name)

_LEGACY_INVESTIGATION_PROMPT_VERSION = _legacy.INVESTIGATION_PROMPT_VERSION
_LEGACY_PARSER_VERSION = _legacy.PARSER_VERSION
MIN_FAMILY_NODE_PORTFOLIO = _legacy.MIN_FAMILY_NODE_PORTFOLIO
TARGET_FAMILY_NODE_PORTFOLIO = _legacy.TARGET_FAMILY_NODE_PORTFOLIO
MAX_FAMILY_NODE_PORTFOLIO = _legacy.MAX_FAMILY_NODE_PORTFOLIO
MIN_COMPETING_MECHANISMS = _legacy.MIN_COMPETING_MECHANISMS
INVESTIGATION_PROMPT_VERSION = "phenotype-investigation-1.17.0"
PARSER_VERSION = "mechanism-hypothesis-parser-1.11.0"

SYSTEM_PROMPT = _legacy.SYSTEM_PROMPT.replace(
    "- one stable database accession (UniProt preferred) for a reference gene in a\n"
    "  non-focal species, plus why you selected that gene rather than a paralogue;",
    "- one real, named reference gene or locus from a literature-supported species whenever the literature provides one, "
    "plus its reference species and, when known, a stable database accession (UniProt preferred); explain why "
    "you selected that target rather than a paralogue;\n"
    "- when no defensible literature-supported reference target is known, leave reference_gene_name, reference_species and "
    "proposed_stable_id empty and state that gap in uncertainty; the importer will preserve it as explicit "
    "review-required missingness, not provider failure;",
).replace(
    "- Do not invent accessions or source ids. Leave the field empty instead. An empty\n"
    "  field is recoverable; a fabricated one is not.",
    "- Do not invent genes, loci, accessions, species, or source ids. A gene-family label alone is not a reference "
    "target. Leave unknown reference-target fields empty instead; a missing target is review-required and "
    "recoverable, while a fabricated pseudo-accession is not.\n"
    "- In a production/assisted investigation, a focal-species primary reference is allowed when explicitly supported "
    "by literature; other species remain valid mechanistic comparators. A strict benchmark route may require a "
    "non-focal reference, but that policy is supplied by the benchmark route.",
)
SYSTEM_PROMPT = SYSTEM_PROMPT.replace(
    "- model-species evidence: a named mutant or functional study in a species other\n"
    "  than the focal species, with a source id (PMID/DOI/UniProt/TAIR accession);",
    "- reported evidence: a named mutant or functional study with its experimental species "
    "and source id (PMID/DOI/UniProt/TAIR accession) when available; the legacy "
    "model_species_evidence field may hold focal or non-focal evidence in production;",
)
SYSTEM_PROMPT += """

The primary approach is to bridge functional and mechanistic evidence from other
species to the focal phenotype. Investigate primary literature with Web search
when available, using the phenotype, process, family and reference gene. Explain
the reported observation, the homolog/family correspondence and the proposed
focal mechanism separately; a plausible transfer is a hypothesis, not a reported
experimental fact. Identical phenotypes and focal-species experimental proof are
not prerequisites for a meaningful candidate hypothesis.
Focal-species direct reports and other-species mechanistic comparators are both
valid research inputs in production; neither is required for every hypothesis.
Use the actual case context, without species-specific or example-specific
shortlists, a required model species, or rules tailored to a benchmark answer.
"""

SYSTEM_PROMPT += LITERATURE_SEARCH_GUIDANCE

# The legacy USER_TEMPLATE embeds the old MECHANISM_PAYLOAD_TEMPLATE verbatim
# (braces doubled for the later ``str.format``). Rebuild it with the same case
# fields and hypothesis schema, but with paper_reports_json shown inside EACH
# hypothesis. This keeps the inherited case fields and the existing parser
# contract while making the per-record paper placement explicit.
_MECHANISM_PAYLOAD_TEMPLATE = (
    '{"hypotheses": [{"causal_layer": "perception|signaling|regulation|execution|alternative|unclear",\n'
    '  "pathway_or_process": "", "gene_family": "", "mechanism_statement": "",\n'
    '  "expected_role": "", "direction_consistency": "", "model_species_evidence": "",\n'
    '  "model_species_mutant_phenotype": "", "source_ids": "", "proposed_stable_id": "",\n'
    '  "reference_species": "", "reference_gene_name": "", "selection_rationale": "",\n'
    '  "counter_hypothesis": "", "discriminating_observation": "",\n'
    '  "confidence": "low|medium|high", "uncertainty": "",\n'
    '  "paper_reports_json": ' + render_initial_paper_report_example() + "}]}"
)
_LEGACY_PAYLOAD_DOUBLED = _legacy.MECHANISM_PAYLOAD_TEMPLATE.replace("{", "{{").replace("}", "}}")
_NEW_PAYLOAD_DOUBLED = _MECHANISM_PAYLOAD_TEMPLATE.replace("{", "{{").replace("}", "}}")
USER_TEMPLATE = _legacy.USER_TEMPLATE.replace(_LEGACY_PAYLOAD_DOUBLED, _NEW_PAYLOAD_DOUBLED)
if USER_TEMPLATE == _legacy.USER_TEMPLATE:
    raise RuntimeError("failed to place paper_reports_json in the Route A payload example")
from .relevance_tiers import INITIAL_TIER_PROMPT
SYSTEM_PROMPT += INITIAL_TIER_PROMPT
SYSTEM_PROMPT += (
    "\nThe JSON example below places paper_reports_json inside each hypothesis, never at "
    "the root: one report object per paper you actually read. Optional research_annotations belongs inside each paper report, never at hypothesis or payload root; keep paper-specific read scope, query notes and checklist reasons there. Replace every placeholder "
    "value with that paper's real observations; never return a placeholder as a fact. Use "
    "paper_reports_json: [] when no paper defensibly supports a hypothesis, and never invent "
    "a citation. Leave an unverifiable checklist answer as unknown rather than no. A paper "
    "report is not required for a hypothesis to be proposed. reference_gene_id selects an "
    "existing reference endpoint, while reported_gene/reported_gene_id remain the "
    "experimentally reported gene.\n"
)

_legacy.SYSTEM_PROMPT = SYSTEM_PROMPT
_legacy.USER_TEMPLATE = USER_TEMPLATE
_legacy.INVESTIGATION_PROMPT_VERSION = INVESTIGATION_PROMPT_VERSION
_legacy.PARSER_VERSION = PARSER_VERSION

_ORIGINAL_PHENOTYPE_INVESTIGATE = _legacy.phenotype_investigate
_ORIGINAL_MECHANISM_IMPORT = _legacy.mechanism_import


def _normalise_species(value: object) -> str:
    return " ".join(str(value or "").strip().casefold().split())


def _target_state(
    row: _pd.Series | dict[str, object],
    focal_species: str,
    benchmark_mode: bool = False,
) -> str:
    gene_or_locus = bool(
        str(row.get("reference_gene_name", "")).strip()
        or str(row.get("proposed_stable_id", "")).strip()
    )
    reference_value = str(row.get("reference_species", "")).strip()
    species = bool(reference_value)
    if not gene_or_locus and not species:
        return "missing_reference_target"
    if gene_or_locus and not species:
        return "missing_reference_species"
    if species and not gene_or_locus:
        return "missing_reference_gene_or_locus"

    focal = _normalise_species(focal_species)
    if not focal:
        return "focal_species_not_available_for_reference_check"
    references = {
        _normalise_species(value)
        for value in _re.split(r"\s*[/;|]\s*", reference_value)
        if _normalise_species(value)
    }
    if len(references) != 1:
        return "ambiguous_reference_species"
    if focal in references and benchmark_mode:
        return "focal_species_used_as_reference_target"
    if focal in references:
        return "declared_reference_target"
    return "declared_non_focal_target"


def _positional_or_keyword(
    args: tuple[object, ...],
    kwargs: dict[str, object],
    name: str,
    position: int,
) -> object | None:
    if name in kwargs:
        return kwargs[name]
    return args[position] if len(args) > position else None


def _focal_species_for_import(
    args: tuple[object, ...],
    kwargs: dict[str, object],
) -> tuple[str, str]:
    investigation_arg = _positional_or_keyword(
        args, kwargs, "investigation_arg", 2
    )
    if not investigation_arg:
        return "", "not_evaluated_investigation_not_supplied"
    investigation_path = _runtime_path(str(investigation_arg))
    investigation = _pd.read_csv(
        investigation_path, sep="\t", dtype=str
    ).fillna("")
    if investigation.empty or "focal_species" not in investigation.columns:
        return "", "not_evaluated_focal_species_missing"
    values = {
        str(value).strip()
        for value in investigation["focal_species"].astype(str)
        if str(value).strip()
    }
    if len(values) != 1:
        return "", "not_evaluated_focal_species_ambiguous"
    return next(iter(values)), "evaluated_from_investigation"


def _benchmark_policy_from_investigation(
    args: tuple[object, ...], kwargs: dict[str, object]
) -> bool:
    investigation_arg = _positional_or_keyword(args, kwargs, "investigation_arg", 2)
    if not investigation_arg:
        return False
    path = _runtime_path(str(investigation_arg))
    manifest_path = path.parent / "phenotype_investigation_manifest.json"
    try:
        manifest = _json.loads(manifest_path.read_text(encoding="utf-8"))
        if not isinstance(manifest, dict):
            raise ValueError("investigation policy manifest must be an object")
        contract = manifest.get("payload_contract", {})
        policy = contract.get("reference_target_contract", {}) if isinstance(contract, dict) else {}
        mode = policy.get("benchmark_mode", False) if isinstance(policy, dict) else False
        if not isinstance(mode, bool):
            raise ValueError("investigation benchmark_mode must be boolean")
        return mode
    except FileNotFoundError:
        return False


def phenotype_investigate(*args: object, **kwargs: object) -> dict[str, object]:
    """Write the strengthened prompt and record its target requirements."""
    benchmark_mode = bool(kwargs.pop("benchmark_mode", False))
    result = _ORIGINAL_PHENOTYPE_INVESTIGATE(*args, **kwargs)
    manifest_path = _Path(str(result["manifest"]))
    manifest = _json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict):
        raise ValueError("phenotype-investigation manifest must be a JSON object")
    manifest["prompt_version"] = INVESTIGATION_PROMPT_VERSION
    if benchmark_mode:
        strict_policy = (
            "\n\nSTRICT BENCHMARK POLICY:\n"
            "Use only a reference gene or locus from a species other than the focal species. "
            "Same-species references are excluded by this benchmark route. Do not infer or name "
            "the benchmark answer, interval candidates, or a causal gene.\n"
        )
        prompt_path = _Path(str(result["prompt"]))
        prompt_text = prompt_path.read_text(encoding="utf-8")
        if strict_policy.strip() not in prompt_text:
            prompt_text += strict_policy
            prompt_path.write_text(prompt_text, encoding="utf-8")
            prompts_path = _Path(str(result["prompts_jsonl"]))
            prompt_record = _json.loads(prompts_path.read_text(encoding="utf-8").strip())
            prompt_record["system"] = str(prompt_record.get("system", "")) + strict_policy
            prompts_path.write_text(_json.dumps(prompt_record, ensure_ascii=False) + "\n", encoding="utf-8")
        manifest["prompt_sha256"] = _hashlib.sha256(prompt_text.encode("utf-8")).hexdigest()
    contract = manifest.setdefault("payload_contract", {})
    if isinstance(contract, dict):
        contract["literature_research"] = dict(LITERATURE_RESEARCH_POLICY)
        contract["initial_hypothesis_target"] = TARGET_FAMILY_NODE_PORTFOLIO
        contract["initial_count_guidance"] = {
            "historical_min": MIN_FAMILY_NODE_PORTFOLIO,
            "historical_max": MAX_FAMILY_NODE_PORTFOLIO,
            "enforced": False,
        }
        contract["reference_target_contract"] = {
            "portfolio_ready_row_requires": [
                "reference_species",
                "reference_gene_name_or_proposed_stable_id",
                "verified_reference_species_relation",
            ],
            "target_must_be_non_focal": benchmark_mode,
            "benchmark_mode": benchmark_mode,
            "unknown_target_policy": "leave fields empty and emit review-required missingness; never fabricate",
            "provider_failure_is_separate": True,
        }
    manifest_path.write_text(
        _json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    result["prompt_version"] = INVESTIGATION_PROMPT_VERSION
    return result


def _append_review(
    reviews: _pd.DataFrame,
    hypothesis_id: str,
    state: str,
) -> _pd.DataFrame:
    if state == "missing_reference_species":
        reason = "mechanism_reference_target_without_species"
        action = "supply the non-focal reference species, or clear the unsupported target identity and keep the row review-required"
    elif state == "missing_reference_gene_or_locus":
        reason = "mechanism_reference_species_without_gene_or_locus"
        action = "supply a real reference gene/locus or stable accession for this species, or clear the species and keep the row review-required"
    elif state == "focal_species_used_as_reference_target":
        reason = "mechanism_reference_target_is_focal_species"
        action = "replace the target with a real gene or locus from a non-focal reference species, or keep the hypothesis review-required"
    elif state == "ambiguous_reference_species":
        reason = "mechanism_reference_species_is_ambiguous"
        action = "choose one explicit non-focal reference species for this target; do not combine alternatives in a portfolio-ready row"
    elif state == "focal_species_not_available_for_reference_check":
        reason = "mechanism_non_focal_target_not_evaluated"
        action = "supply the phenotype investigation table so the declared reference species can be checked against the focal species"
    else:
        reason = "mechanism_without_reference_target"
        action = "supply a real non-focal reference gene/locus plus species when known; otherwise retain this explicit unresolved target"
    addition = _pd.DataFrame(
        [
            {
                "reason": reason,
                "target": hypothesis_id,
                "required_next_action": action,
                "allowed_choices": "supply verified target | keep as unresolved reference-target gap | discard hypothesis",
            }
        ]
    )
    return _pd.concat([reviews, addition], ignore_index=True)


def mechanism_import(*args: object, **kwargs: object) -> dict[str, object]:
    """Require a declared, checked literature reference target for portfolio readiness."""
    benchmark_mode = bool(kwargs.pop("benchmark_mode", False))
    benchmark_mode = benchmark_mode or _benchmark_policy_from_investigation(args, kwargs)
    result = _ORIGINAL_MECHANISM_IMPORT(*args, **kwargs)
    mechanisms_path = _Path(str(result["mechanisms"]))
    mechanisms = _pd.read_csv(mechanisms_path, sep="\t", dtype=str).fillna("")
    focal_species, focal_check_state = _focal_species_for_import(args, kwargs)
    initial_tier_present = "initial_tier" in mechanisms and mechanisms["initial_tier"].astype(str).str.strip().ne("").any()
    standalone_legacy = focal_check_state == "not_evaluated_investigation_not_supplied" and not initial_tier_present
    if standalone_legacy:
        # Preserve byte-compatible standalone imports. The strengthened parser
        # contract is attached to an investigation generation; without that
        # artifact the non-focal check is explicitly not evaluated rather than
        # silently rewriting historical payload provenance.
        mechanisms["prompt_version"] = _LEGACY_INVESTIGATION_PROMPT_VERSION
        mechanisms["parser_version"] = _LEGACY_PARSER_VERSION
        mechanisms.to_csv(
            mechanisms_path, sep="\t", index=False, lineterminator="\n"
        )
    complete_states = {"declared_non_focal_target"}
    if not benchmark_mode:
        complete_states.add("declared_reference_target")
    states = {
        str(row.get("hypothesis_id", "")): _target_state(row, focal_species, benchmark_mode)
        for row in mechanisms.to_dict("records")
    }
    missing_ids = sorted(
        hypothesis_id
        for hypothesis_id, state in states.items()
        if state not in complete_states
    )

    reviews_path = _Path(str(result["review_required"]))
    reviews = _pd.read_csv(reviews_path, sep="\t", dtype=str).fillna("")
    for hypothesis_id in missing_ids:
        reviews = _append_review(reviews, hypothesis_id, states[hypothesis_id])
    if not reviews.empty:
        reviews = reviews.drop_duplicates().sort_values(["reason", "target"], kind="stable")
    reviews.to_csv(reviews_path, sep="\t", index=False, lineterminator="\n")

    manifest_path = _Path(str(result["manifest"]))
    manifest = _json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict):
        raise ValueError("mechanism-import manifest must be a JSON object")
    counts: dict[str, int] = {}
    for state in states.values():
        counts[state] = counts.get(state, 0) + 1
    contract_met = not missing_ids
    manifest["parser_version"] = (
        _LEGACY_PARSER_VERSION
        if standalone_legacy
        else PARSER_VERSION
    )
    manifest["prompt_version"] = (
        _LEGACY_INVESTIGATION_PROMPT_VERSION
        if standalone_legacy
        else INVESTIGATION_PROMPT_VERSION
    )
    manifest["reference_target_status_counts"] = dict(sorted(counts.items()))
    manifest["hypotheses_without_complete_reference_target"] = missing_ids
    manifest["reference_target_contract_met"] = contract_met
    manifest["benchmark_mode"] = benchmark_mode
    manifest["reference_target_contract"] = {
        "required_for_portfolio_ready": True,
        "complete_state": "declared_non_focal_target" if benchmark_mode else "declared_reference_target",
        "identity_requirement": "reference_gene_name or proposed_stable_id",
        "species_requirement": ("exactly one reference_species distinct from focal_species"
                                if benchmark_mode else "exactly one literature-supported reference_species"),
        "benchmark_mode": benchmark_mode,
        "focal_species": focal_species,
        "focal_species_check_state": focal_check_state,
        "unknown_target_policy": "review_required_missing_target",
        "provider_failure_state": "not assessed at mechanism import; recorded only by downstream resolver/provider",
        "fabricated_identifier_policy": "never generated or repaired by this stage",
    }
    prior_portfolio = bool(manifest.get("family_node_portfolio_met"))
    manifest["family_node_portfolio_met"] = bool(prior_portfolio and contract_met)
    manifest["competing_mechanisms_met"] = bool(
        int(manifest.get("distinct_mechanisms", 0) or 0)
        >= int(manifest.get("min_competing_mechanisms", MIN_COMPETING_MECHANISMS) or MIN_COMPETING_MECHANISMS)
    )
    manifest["review_required_count"] = int(len(reviews))
    manifest["output_sha256"] = _sha256(mechanisms_path)
    manifest_path.write_text(
        _json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    result["competing_mechanisms_met"] = manifest["competing_mechanisms_met"]
    result["reference_target_contract_met"] = contract_met
    return result


_legacy.phenotype_investigate = phenotype_investigate
_legacy.mechanism_import = mechanism_import
