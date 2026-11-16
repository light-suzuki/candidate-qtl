"""Cheap Tier derivation from original observations, never a content reviewer."""
from __future__ import annotations

import json
import re
import unicodedata
from collections.abc import Mapping

from .identifiers import canonical_accession, classify_gene_name_field
from .reference_name_pairs import reference_handle

POLICY = "observation_rules_v1"
JUDGMENT_POLICY = "reported_paper_judgment_v1"
CHECKLIST_POLICY = "paper_checklist_v1"
CHECKLIST_FIELDS = (
    "focal_species", "same_gene_as_focal", "phenotype_reported",
    "equivalent_target_phenotype", "similar_cellular_or_developmental_phenotype",
    "shared_specific_process", "function_supports_hypothesis",
    "annotation_supports_hypothesis", "identity_resolved",
)
# A post-freeze identity assessment is a stricter same-species member
# correspondence than a qualified sequence hit (unique exact whole-proteome
# model plus explicit review). Its projection basis is accepted only on rows
# that carry the identity contract, so no other producer can claim it.
POST_FREEZE_IDENTITY_CONTRACT = "post-freeze-identity-assessment-1.0.0"
POST_FREEZE_IDENTITY_BASIS = "post_freeze_exact_identity"


def display_text(value: object) -> str:
    return unicodedata.normalize("NFKC", str(value or "")).strip()


def is_unknown(value: object) -> bool:
    return display_text(value).casefold() in {"unknown", "none", "unresolved", "n/a", "-", "不明"}


def unfilled_answer(value: object) -> bool:
    return value is None or (isinstance(value, str) and not display_text(value))


def _answer(value: object) -> object:
    if value is None:
        return "unknown"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, str):
        text = display_text(value).casefold()
        return {"true": "yes", "false": "no", "はい": "yes", "いいえ": "no", "不明": "unknown", "": "unknown"}.get(text, text)
    return value


def normalize_checklist(value: object) -> str:
    if value is None or value == "":
        return ""
    checks = json.loads(value) if isinstance(value, str) else value
    if not isinstance(checks, dict) or set(checks) - set(CHECKLIST_FIELDS):
        raise ValueError("paper_checklist must contain only the documented questions")
    # JSON booleans and harmless display case carry the same factual answers.
    checks = {key: _answer(answer) for key, answer in checks.items()}
    if any(answer not in ("yes", "no", "unknown") for answer in checks.values()):
        raise ValueError("paper_checklist answers must be yes, no or unknown")
    return json.dumps({key: checks.get(key, "unknown") for key in CHECKLIST_FIELDS}, sort_keys=True, separators=(",", ":"))


def checklist_tier(value: object) -> str:
    """One short evidence decision table; no semantic second pass.

    ``same_gene_as_focal`` is intentionally not a Tier-table predicate: the
    downstream member-bound ``classify_observation`` checks qualified member
    projection separately. Preserve the original ``same_gene_as_focal`` answer.
    """
    checks = json.loads(normalize_checklist(value) or "{}")
    yes = lambda key: checks.get(key) == "yes"
    if yes("phenotype_reported"):
        if yes("focal_species") and (yes("equivalent_target_phenotype") or yes("similar_cellular_or_developmental_phenotype")):
            return "T1"
        if yes("equivalent_target_phenotype") or yes("similar_cellular_or_developmental_phenotype"):
            return "T2"
        if yes("shared_specific_process"):
            return "T3"
    if yes("function_supports_hypothesis"):
        return "T4"
    if yes("annotation_supports_hypothesis"):
        return "T5"
    return "T6" if yes("identity_resolved") else "U"
GENE_RELATIONS = {"same_gene", "homolog", "different_gene", "unknown"}
OBSERVATION_FIELDS = ("evidence_types", "reported_function", "reported_annotation")
COMPARISONS = {"same_equivalent", "same_cellular_process", "shared_process", "different", "unknown"}


def normalize_evidence_types(value: object) -> str:
    if value is None or value == "":
        return ""
    types = json.loads(value) if isinstance(value, str) else value
    if not isinstance(types, list) or any(not isinstance(item, str) for item in types):
        raise ValueError("evidence_types must be an array of strings")
    types = [display_text(item).casefold() for item in types]
    if set(types) - {"phenotype", "function", "annotation", "unknown"}:
        raise ValueError("unknown evidence_types value")
    return json.dumps(sorted(set(types)), separators=(",", ":"))


def has_observations(row: Mapping) -> bool:
    return bool(normalize_checklist(row.get("paper_checklist")) or normalize_evidence_types(row.get("evidence_types")))


def _text(value: object) -> str:
    text = display_text(value)
    return "" if is_unknown(text) else text


def gene_endpoint_key(value: object) -> str:
    """Normalize identifier syntax, retaining namespaces and accession versions."""
    text = _text(value)
    if not text:
        return ""
    handle = reference_handle(text, allow_untyped=True)
    if handle is None and ":" not in text:
        handle = reference_handle(canonical_accession(text), allow_untyped=True)
    if handle:
        kind, identifier = handle
        return f"{kind}:{identifier}"  # the existing parser owns namespace case rules
    return "name:" + " ".join(text.split()).casefold()


def _species_key(value: object) -> str:
    """Normalize display formatting only; never guess taxonomic synonyms."""
    return re.sub(r"\s+\([^()]+\)$", "", " ".join(_text(value).split())).casefold()


CHECKLIST_CONSISTENCY_FIELDS = ("focal_species_label", "species_comparison_state", "consistency_warning")


def checklist_consistency(row: Mapping) -> dict[str, str]:
    """Audit a T1 answer against supplied labels, without resolving taxonomy.

    A paper's reported species and the case species are the only endpoints.
    Reference/proxy species, gene names and homology never decide this audit.
    Different complete binomial-shaped labels or explicit NCBI taxon IDs are
    a concrete correction lead, not proof of different biological species:
    synonyms and taxonomic ranks have not been resolved. Labels outside those
    shapes stay unknown; shape recognition does not verify a scientific name.
    No lookup, Tier change or candidate/support decision is performed here.
    """
    if not normalize_checklist(row.get("paper_checklist")):
        return {}
    focal = display_text(row.get("focal_species"))
    reported = display_text(row.get("reported_species"))
    result = {"focal_species_label": focal, "species_comparison_state": "not_required_for_non_t1",
              "consistency_warning": ""}
    if checklist_tier(row.get("paper_checklist")) != "T1":
        return result
    result["species_comparison_state"] = "unresolved_labels"
    if not _text(focal) or not _text(reported):
        return result
    if _species_key(focal) == _species_key(reported):
        return {**result, "species_comparison_state": "matching_labels"}

    def explicit_label(value: str) -> tuple[str, str] | None:
        value = " ".join(value.split())
        taxon = re.fullmatch(r"(?:NCBI\s*)?(?:Taxon|Taxid)\s*:\s*([1-9][0-9]*)", value, re.I)
        if taxon:
            return "ncbi_taxon_id", taxon[1]
        # Restrict disagreement warnings to fully written binomial-shaped
        # labels. Never expand abbreviated genera or guess a common-name alias.
        if re.fullmatch(r"[A-Z][a-z]+ [a-z][a-z-]+", value):
            return "binomial_label", value.casefold()
        return None

    left, right = explicit_label(focal), explicit_label(reported)
    if not left or not right or left[0] != right[0]:
        return result
    if left == right:
        return {**result, "species_comparison_state": "matching_labels"}
    return {**result, "species_comparison_state": "explicit_label_disagreement",
            "consistency_warning": "T1 focal_species=yes disagrees with the supplied reported/focal species labels; "
                "taxonomic identity, synonymy and rank were not resolved. Recheck the original record; "
                "Tier, candidate membership and support are unchanged."}


def gene_binding(row: Mapping, reported_gene: str, reported_species: str = "") -> str:
    """Use existing named endpoints; no guessed aliases or online lookup."""
    def endpoints(value):
        result, explicit = set(), set()
        for text in display_text(value).split(";"):
            text = text.strip()
            parts = [text] if re.match(r"^(sp|tr)\|", text, re.I) else text.split("|")
            for part in parts:
                decorated = re.fullmatch(r"(.+?)\s*\(([^()]+)\)", part.strip())
                if decorated:
                    # An explicit accession wins over a display symbol. Do not
                    # match a shared symbol when its supplied accession differs.
                    for candidate in (decorated[2].strip(), decorated[1].strip()):
                        accession = canonical_accession(candidate)
                        if reference_handle(candidate, allow_untyped=True) or classify_gene_name_field(accession) in {"uniprot_accession", "agi_locus", "locus_tag"}:
                            part = candidate if reference_handle(candidate, allow_untyped=True) else accession
                            break
                key = gene_endpoint_key(part)
                if key:
                    result.add(key)
                    if not key.startswith("name:") or classify_gene_name_field(canonical_accession(part)) in {"locus_tag"}:
                        explicit.add(key)
        return result, explicit
    references = set()
    for key in ("reference_gene", "reference_gene_id", "reference_gene_name", "source_reference_gene_id", "proposed_stable_id"):
        references.update(endpoints(row.get(key))[0])
    references.discard("")
    all_reported, explicit = endpoints(reported_gene)
    reported = explicit or all_reported
    reported.discard("")
    if not references or not reported:
        return "unknown"
    reference_species = {_species_key(species) for species in str(row.get("reference_species", "")).split(";") if _text(species)}
    if _text(reported_species) and reference_species and _species_key(reported_species) not in reference_species:
        return "unresolved"
    return "matched" if references & reported else "unresolved"


def classify_observation(row: Mapping) -> tuple[str, str, str, list[str]]:
    """Return short Tier, basis, reused observation text, unresolved fields.

    A lower supported interpretation can coexist with missing stronger facts.
    Missing data never triggers another research task or a membership veto.
    """
    from .bibliography import concrete_bridge
    from .source_verify import classify_citation
    types = set(json.loads(normalize_evidence_types(row.get("evidence_types")) or "[]"))
    missing = []
    source = _text(row.get("source_id") or row.get("source_ids"))
    kind = classify_citation(source)[0]
    paper = kind in {"pmid", "pmcid", "doi"}
    identity = any(_text(row.get(key)) for key in
                   ("reference_gene", "reference_gene_name", "proposed_stable_id", "family_id", "domain_ids"))
    if row.get("verification_status") == "contradicted" or row.get("source_verification_state") == "contradicted":
        return "T6" if identity else "U", "observation_contradicted", "Recorded contradiction retained.", ["contradiction"]
    if row.get("projection_state") and row["projection_state"] != "inside_current_qtl":
        return "T6" if identity else "U", "observation_projection_unresolved", "Projection is not inside the current QTL.", ["projection"]
    if row.get("tier_classification_policy") == JUDGMENT_POLICY:
        # The original research LLM owns semantic correspondence and Tier.
        # Machines check the supplied reference anchor, not names or synonyms.
        if not paper or row.get("tier_gene_binding") != "matched":
            return "U", "reported_judgment_unresolved", "Original paper/reference correspondence is unresolved.", ["gene_correspondence"]
        label = _text(row.get("initial_tier")).split("_", 1)[0]
        reason = _text(row.get("initial_tier_reason"))
        if label not in {"T1", "T2", "T3", "T4", "T5", "T6", "U"} or not reason:
            return "U", "reported_judgment_incomplete", "Original paper Tier or reason is missing.", ["paper_judgment"]
        if label == "T1" and row.get("gene_correspondence") == "homolog":
            return "U", "reported_judgment_inconsistent", "Homology alone cannot authorize a focal-gene T1.", ["gene_correspondence"]
        return label, "original_paper_judgment", reason, (["paper_judgment"] if label == "U" else [])
    if row.get("tier_classification_policy") == CHECKLIST_POLICY:
        if not paper or row.get("tier_gene_binding") != "matched":
            return "U", "checklist_anchor_unresolved", "Paper/reference correspondence is unresolved.", ["gene_correspondence"]
        label = checklist_tier(row.get("paper_checklist"))
        if label == "T1":
            members = _text(row.get("candidate_gene_id") or row.get("member_interval_genes"))
            projection = row.get("homolog_projection_json") or "{}"
            projection = json.loads(projection) if isinstance(projection, str) else projection
            bases = {"qualified_sequence_hit", "frozen_gene_orthology"}
            if row.get("post_freeze_identity_contract") == POST_FREEZE_IDENTITY_CONTRACT:
                bases.add(POST_FREEZE_IDENTITY_BASIS)
            linked = (projection.get("basis") in bases
                      and projection.get("candidate_gene_id") in members.split(";")
                      and gene_endpoint_key(projection.get("reference_gene_id")) == gene_endpoint_key(row.get("reference_gene_id"))
                      and projection.get("candidate_gene_id") in projection.get("member_genes", [])
                      and projection.get("evidence"))
            # Paper-side T1 is conditional on projection. At the member boundary,
            # only the existing verified producer may supply that correspondence.
            # An original same-gene answer cannot replace the member projection.
            if not _text(row.get("initial_tier_reason")) or (members and not linked):
                return "T2", "focal_report_projection_unresolved", _text(row.get("initial_tier_reason")), ["homolog_projection"]
            return "T1", "same_species_reported_homolog_candidate", _text(row.get("initial_tier_reason")), []
        # gene_correspondence joins the reported gene to a reference/bait;
        # same_gene_as_focal answers a different question about the focal locus.
        # A proxy bait must not erase a grounded same-locus correspondence.
        return label, "original_paper_checklist", _text(row.get("initial_tier_reason")), (["paper_checklist"] if label == "U" else [])
    if paper and row.get("tier_gene_binding") != "matched":
        return "U", "observation_gene_binding_unresolved", "Reported gene is not resolved to this reference endpoint.", ["gene_correspondence"]
    if not source and types:
        missing.append("source_id")
    if "unknown" in types:
        missing.append("evidence_types")
    if "phenotype" in types:
        reported = _text(row.get("reported_phenotype"))
        comparison = _text(row.get("phenotype_comparison"))
        if not paper:
            missing.append("paper_source")
        if not reported:
            missing.append("reported_phenotype")
        if not comparison:
            missing.append("phenotype_comparison")
        if paper and reported:
            focal = _text(row.get("focal_species"))
            species = _text(row.get("reported_species"))
            if comparison == "same_equivalent":
                if focal and species and _species_key(focal) == _species_key(species):
                    return "T1", "observed_focal_trait", reported, sorted(set(missing))
                if not focal or not species:
                    missing.append("species_comparison")
                return "T2", "observed_similar_phenomenon", reported, sorted(set(missing))
            if comparison == "same_cellular_process":
                return "T2", "observed_similar_phenomenon", reported, sorted(set(missing))
            if comparison == "shared_process":
                process = _text(row.get("reported_process"))
                if process:
                    return "T3", "observed_shared_process", f"{reported}; shared process: {process}", sorted(set(missing))
                missing.append("reported_process")
    for evidence_type, field, label in (("function", "reported_function", "T4"), ("annotation", "reported_annotation", "T5")):
        if evidence_type not in types:
            continue
        observation = _text(row.get(field))
        if not observation:
            missing.append(field)
        function_source = paper or kind == "database_accession"
        if evidence_type == "function" and not function_source:
            missing.append("paper_source")
        if observation and source and (function_source or evidence_type == "annotation"):
            if concrete_bridge({**row, "biological_function": observation}):
                return label, f"observed_{evidence_type}_hypothesis", observation, sorted(set(missing))
            missing.append("specific_phenotype_bridge")
    if missing:
        return "U", "observation_unresolved", "Original observation is incomplete; no additional review requested.", sorted(set(missing))
    return ("T6", "observed_identity_only", "Identity retained without a specific phenotype bridge.", []) if identity else (
        "U", "observation_identity_unresolved", "No usable identity supplied.", ["identity"])
