"""Trait relevance annotations, independent of support, mapping and copy state.

This module consumes already identity/lineage-checked records; it is not a source
verifier or a membership selector. Structured bridge text remains a biological
proposal when its source verification is incomplete. No score is calculated.
"""
from __future__ import annotations

import json
import re
from collections.abc import Iterable, Mapping
from typing import Any
from .paper_reports import PAPER_FIELD, PAPER_POLICY, REPORT_FIELDS, RESEARCH_ANNOTATIONS, bound_paper_rows, paper_summary
from .observation_tiers import CHECKLIST_POLICY, JUDGMENT_POLICY, OBSERVATION_FIELDS, POLICY as OBSERVATION_POLICY, classify_observation, has_observations

CONTRACT = "candidate-relevance-2.5.0"
INITIAL_TIER_POLICY = "reported_function_v2"
INITIAL_FIELDS = ("initial_tier", "initial_tier_reason", "reported_phenotype", "reported_process",
                  "phenotype_comparison", "reported_organ", "reported_stage", PAPER_FIELD, *OBSERVATION_FIELDS)
INITIAL_TIER_PROMPT = """
During the original investigation, fill the factual checklist for each paper.
The common A/B output is the paper title/ID, this assessment and short reason,
and the reported gene's name/ID and species. Use source_title, source_id,
paper_checklist, initial_tier_reason, reported_gene, reported_gene_id and
reported_species. Preserve gene IDs with their original database namespace;
leave unavailable names/IDs empty or unknown instead of inventing them or
starting an extra search just to fill fields.
Do not decide Tier from scratch: a small decision table derives it from your
yes/no/unknown answers. Supply paper_checklist as a JSON object with these keys:
focal_species: is the reported gene from the focal species?
same_gene_as_focal: is it the same focal gene, including a grounded assembly/
reference-ID correspondence, rather than homology alone?
phenotype_reported: does this paper actually report a phenotype for the gene?
equivalent_target_phenotype: is that reported phenotype equivalent to the target?
Judge equivalence at the specificity actually supplied in the target description.
Do not require matching direction, developmental stage or allele when the input
does not specify them, and do not invent those details. This is correspondence
of a reported phenotype, not proof of a causal allele in this particular case.
An organ-name overlap or generic function alone is not sufficient: explain how
the reported observation matches the supplied phenotype at that stated scale.
Benchmark Type labels do not determine Tier; any Type may yield T1 when the
actual paper and gene correspondence satisfy the same checklist.
similar_cellular_or_developmental_phenotype: is the reported phenotype very close
to the target at a specific cellular/developmental scale? Record the observed
phenomenon and the comparison reason. Organ overlap or generic function is not a yes.
shared_specific_process: does another reported phenotype share a specific process?
function_supports_hypothesis: does reported function support the proposed mechanism?
annotation_supports_hypothesis: does specific annotation support that mechanism?
identity_resolved: is the reported gene/reference identity available?
Every answer is exactly yes, no or unknown; missing answers mean unknown, not no.
Record all nine answers explicitly. A blank is unfilled recording work, not a
checked unknown. The importer preserves this distinction while calculating the
same factual Tier; missing fields never become positive evidence.
T1 = phenotype_reported + focal_species + (equivalent_target_phenotype or
similar_cellular_or_developmental_phenotype), conditional on a grounded homolog
projection into the QTL. This is a high-priority first-stage Wet candidate,
not proven same-locus identity or verified causal support. Keep every eligible
copy; do not require one-to-one orthology, reciprocal best hit, the same
annotation model or a full-length 100-percent match. Keep same_gene_as_focal
unknown/no when appropriate; it is not a T1 prerequisite.
Otherwise T2 = phenotype_reported + (equivalent_target_phenotype or
similar_cellular_or_developmental_phenotype); T3 = phenotype_reported +
shared_specific_process; T4 = function_supports_hypothesis;
T5 = annotation_supports_hypothesis; T6 = identity_resolved; otherwise U.
Choose the first satisfied rule. Do not supply initial_tier for checklist records;
the machine fills it. Keep one short initial_tier_reason explaining the source
facts and correspondence, with reported_phenotype/function when available.
For each paper/reported gene, put a record in paper_reports_json (JSON array;
JSON text in TSV). Required: source_id (PMID/PMCID/DOI), reported_gene, and
reference_gene_id, gene_correspondence, paper_checklist and initial_tier_reason.
To avoid repeating a review's title, paper_reports_json also accepts
{"bibliography":{"p1":{"source_id":"PMID:...","source_title":"exact LLM-copied title"}},
 "observations":[{"paper_ref":"p1", ...each observation's own fields...}]}.
Only ID/title are shared within that payload. Keep reported_gene/reported_gene_id,
species, phenotype, experiment/context, reason and all nine answers on each
observation. Different findings in the same review remain separate; no semantic
deduplication, extra matching research or Tier inheritance is needed.
reported_gene may be the gene ID/name in the paper (NCBI Gene, TAIR, UniProt or
another namespace); reference_gene_id selects an existing reference endpoint in
this claim, not an invented ID. Judge gene_correspondence from the source and
supplied reference context: same_gene, homolog, different_gene or unknown.
Names, accession namespaces and reference releases may differ: no alias
dictionary or exact display-name match is required. Explain correspondence and
the source facts in initial_tier_reason. Do not confuse the paper/reference relation with
the separate reference-to-focal homology/projection or claim verified support.
Keep the reported gene distinct from a proxy reference used to obtain sequence.
For example, a Brassica phenotype paper can be connected through an Arabidopsis
homolog when the Brassica sequence is unavailable. reported_gene/id/species
remain Brassica, reference_gene_id selects the existing Arabidopsis reference,
gene_correspondence is homolog, and the short reason records the chain and its
basis. The phenotype evidence remains the Brassica paper (T2 for a similar
phenomenon in a pea case); do not replace it with an Arabidopsis function paper
or require a phenotype report for the proxy itself. same_gene_as_focal concerns
the reported gene versus the focal locus, not the proxy. Grounded same-locus
mapping across reference versions can answer yes; merely falling somewhere in
the same broad QTL cannot. Distinct focal paralogs remain distinct.
Use evidence_types as an array with any of phenotype,
function, annotation; multiple types are allowed. Use ["unknown"] for unknown
basis, [] for identity only. Never infer evidence_types from a desired Tier.
Reuse reported_species, reported_phenotype, reported_process and
phenotype_comparison: same_equivalent, same_cellular_process, shared_process,
different, unknown. For function/annotation evidence supply a short
reported_function/reported_annotation from the source; keep the proposed target
mechanism in the existing hypothesis fields. Do not turn an inferred phenotype
into an observed one. Unknown is not false: omit unavailable fields or use unknown.
T1: a focal-species gene/trait report with the same/equivalent or very close
specific cellular/developmental phenotype and a grounded interval homolog.
If the reported gene sequence is unavailable, explain its correspondence to an
available sequence proxy and project that proxy. Preserve reported gene/species,
actual proxy ID/species/sequence, projection and ambiguous copies separately.
No extra proxy phenotype paper or separate Tier reviewer is required. Generic
domains, short local similarity and an unexplained proxy chain do not suffice.
T2: reported same/similar cellular or
developmental phenotype in this or another plant with applicable correspondence.
T3: another phenotype linked by a specific shared process. T4: reported function
supporting a new mechanism hypothesis. T5: specific annotation supporting a
hypothesis. T6: identity only. U: unresolved. Homology alone supplies no phenotype.
An original function judgment can cite an identified database source instead of
a paper: retain evidence_types=["function"], reported_function and the concrete
phenotype bridge. A domain/annotation alone remains ["annotation"] (T5); an ID
alone is not reported function. Source/identity checks and review status remain.
Investigate alternate reference/assembly IDs and model-plant homologs during the
original research, using available authoritative cross-references and sequence/
homology context. Keep GenBank/EMBL, RefSeq, NCBI Gene, UniProt, TAIR,
Ensembl Plants, Gramene or Phytozome identifiers in their original namespace.
Different identifier spelling is not evidence of a different gene; a shared name
alone does not establish the same copy. Unknown mappings stay explicit.
The machine derives each paper's initial Tier from the decision table before
selecting the strongest-Tier evidence papers to carry forward. Bibliography
ID/title checks consume that same selected set. Existing reference-anchor and
member checks still apply, and the gene keeps its strongest applicable paper
Tier. The machine does not re-interpret names, phenotype meaning or the source's
biology.

The same paper reports also produce paper_evidence_table.tsv, one row per
paper/reported gene within a proposal. Set paper_table_policy=paper_evidence_table_v1.
Include publication_type (original or
review, or unknown), reported_gene_id in its namespace (unknown if unavailable),
source_excerpt (a short factual observation) and source_location (abstract,
section, figure or page actually inspected). Keep proposed mechanisms separate
from source observations. Explain explicit unknowns in unknown_reason, including
what source/context was checked. A missing gene ID alone is not weaker biological
evidence; keep the known gene name/species and explain the unavailable ID.
Use followup_question only for a concrete source-backed lead whose clarification
could change the applicable Tier, not an imagined positive answer to every unknown.
Use followup_state=not_attempted initially. The existing research owner receives
at most one targeted follow-up batch for unfilled records and concrete upper-Tier
leads; answer unknown with checked sources/reasons and followup_state=completed if
still unresolved. Do not repeat resolved work or chase annotation-only T5 merely
to raise it. Prefer already applicable T1/T2 papers over incomplete weaker papers
when selecting the useful five; never pad the panel or assign a speculative Tier.
Reuse already established exact gene/ID correspondence within this route, keeping
reported gene, reference endpoint and projected focal copies separate. Repeated
gene names alone do not establish the same ID. Every paper retains its own facts
and checklist. DOI/title verification stays in the existing cache-first metadata
path; a generated correction request is not evidence of model execution.

Keep uncertain or conflicting correspondence explicit; never guess a sole copy.
Each source must be declared among source_ids.
Supply source_title and source_identifier_aliases when known; aliases must be
IDs of the SAME paper. Annotation-only claims can use top-level evidence_types
and reported_annotation; use paper_reports_json=[] if there is no paper report.
Originals and reviews are both acceptable, with no count/type score bonus.
Ordinary literature investigation: screen about 20 distinct candidate papers
per original reference gene or shared mechanism, using titles/abstracts first;
this is not 20 full-text reads per interval gene. Reuse shared research across
copies while preserving each paper/gene binding. Record each examined paper's
factual checklist and short reason first. The importer derives its initial Tier,
then selects up to five useful distinct papers in strongest-Tier order for the
next stage: T1, T2, T3, T4, T5, T6, U. Keep all five selected papers' observations,
gene correspondence, Tier and reason available, so a mistake or inapplicable
finding in one report does not erase the other retained evidence. This selection
cannot repair an incorrect initial checklist by itself.
One found paper is not a completed literature investigation. If the strongest
paper is inaccessible or inapplicable, try another access path or the next
paper. If searches yield only T5 annotation, change queries toward homologs,
other species, mutant/complementation phenotypes, specific processes and reviews;
100 annotation-only citations add no evidence strength. Do not fabricate a Tier
or pad the five reports. If fewer useful papers are obtainable, record the
screened IDs, queries, access failures and shortfall in the research handoff;
keep usable evidence and otherwise eligible candidates.
Search screening is provisional; record unknowns and retrieval failures rather
than inventing a paper's Tier. Five is the evidence handoff target per original
reference gene or shared mechanism. DOI/PMID aliases share one paper slot.
The importer marks the selected reports primary and the remainder audit_only;
bibliography verification consumes the selected set with the existing
cache-first ID/title checks, without choosing another set of five. Preserve all
original observations, unselected reports and counterevidence in the audit
record. The next stage reads the primary evidence for every original panel;
copies reuse their panel rather than opening five new paper slots each.
Known publication_type, underlying_source_ids and corroboration_note are optional;
unknown independence does not trigger extra checking. Counts are traceability
metadata, not independent replication or causal scores. Repeated routes/aliases
add no bonus. The gene uses its strongest applicable paper interpretation.
No Tier-only search, rereading, separate agent or extra review stage. Requesting
a DOI/reason makes a judgment traceable; it does not prove its biological truth.
Use the existing A/B literature checks at their existing strength; Tier adds no
dedicated lookup, reviewer or content-verification gate. Incomplete
fields do not veto an otherwise eligible first-stage candidate.
A grounded proposal is adopted as a first-stage candidate after bibliography and existing
identity/member correspondence checks; this does not promote verified support.
Reuse the original observations in stage two and preserve known counterevidence.
"""
T1 = "T1_same_species_known"
T2 = "T2_trait_or_mechanism_linked"
T3 = "T3_functionally_plausible"
T4 = "T4_identity_only"
UNRESOLVED = "U_unresolved"
LEGACY_RELEVANCE_TIERS = (T1, T2, T3, T4, UNRESOLVED)
NEW_TIERS = ("T1_reported_focal_trait", "T2_reported_similar_phenomenon",
             "T3_reported_shared_process", "T4_reported_function_hypothesis",
             "T5_annotation_hypothesis", "T6_identity_only", UNRESOLVED)
# Distinct labels preserve old scopes and old meanings. Relative legacy order
# is unchanged; no frozen generation is relabelled as new scientific evidence.
RELEVANCE_TIERS = (T1, NEW_TIERS[0], T2, NEW_TIERS[1], NEW_TIERS[2], T3,
                   NEW_TIERS[3], NEW_TIERS[4], T4, NEW_TIERS[5], UNRESOLVED)
RELEVANCE_COLUMNS = (
    "relevance_tier", "relevance_basis", "tier_reason", "support_state",
    "mapping_tier", "identity_level", "subfamily_or_orthogroup", "copy_state",
    "route_provenance", "family_id", "reference_gene", "reference_species", "source_ids", "legacy_support_tier",
    "relevance_migration_state", "relevance_evidence_json",
    "relevance_policy", "initial_tier", "initial_tier_reason", "candidate_admission_policy",
    "candidate_admission_state",
    "tier_aggregation", "tier_source_id", "tier_evidence_json",
    "paper_tier_counts_json", "paper_evidence_summary_json",
    "tier_classification_policy", "tier_classification_state", "tier_missing_fields_json", "tier_gene_binding",
)
_MISSING = {"", "unknown", "unresolved", "not_evaluated", "none", "nan", "n/a", "-"}
_IDENTITY_FIELDS = (
    "family_id", "family_ids", "evolutionary_family_id", "orthogroup", "orthogroup_id",
    "orthogroup_ids", "subfamily_id", "subfamily_or_clade", "subfamilies", "domain_ids",
    "domain_architecture", "reference_gene", "reference_genes", "reference_gene_id",
    "orthology_reference_gene", "hit_accession", "identity_accessions",
    "source_reference_gene_id", "biological_function", "function_description",
)
_STRONG_BRIDGES = {"direct_functional_bridge", "mutant_phenotype_bridge", "mechanistic_process_bridge"}
_VERIFIED_SUPPORT = {"verified", "verified_semantic_support", "verified_direct_support"}


def text(value: object) -> str:
    value = str(value if value is not None else "").strip()
    return "" if value.casefold() in _MISSING else value


def _first(row: Mapping[str, Any], *keys: str) -> str:
    return next((text(row.get(key)) for key in keys if text(row.get(key))), "")


def _tokens(value: object) -> set[str]:
    return {part for token in re.split(r";|,|\s\|\s", text(value)) if (part := text(token))}


def _join(values: Iterable[object]) -> str:
    return ";".join(sorted({token for value in values for token in _tokens(value)}))


def _json(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _known_same_trait(row: Mapping[str, Any]) -> bool:
    """Require one bound record, not a synthetic join of unrelated assertions.

    direct_related_trait is deliberately NOT same/equivalent trait. The existing
    Direct-C vocabulary does not establish equivalence for that broader label.
    Support status is never used as a shortcut for these checks.
    """
    focal, reference = text(row.get("focal_species")), text(row.get("reference_species"))
    return bool(
        focal and reference and focal.casefold() == reference.casefold()
        and row.get("species_distance") in {"S0_same_cultivar_or_accession", "S1_same_species"}
        and row.get("trait_relation_state") == "direct_target_trait"
        and text(row.get("query_phenotype")) and text(row.get("reported_phenotype"))
        and text(row.get("source_id")) and text(row.get("source_reference_gene_id"))
        and row.get("source_verification_state") == "verified"
        and row.get("gene_identity_state") in {"verified_gene", "verified_locus_only"}
        and row.get("direct_rescue_eligibility") == "eligible_direct_record"
        and row.get("projection_state") == "inside_current_qtl"
        and row.get("projection_basis") in {"exact_focal_id", "assembly_or_alias_mapping"}
        and str(row.get("claim_binding_projection_authorized", "")).casefold() in {"true", "1"}
        and row.get("claim_binding_status") == "verified_complete"
    )


def classify_relevance(row: Mapping[str, Any]) -> tuple[str, str, str]:
    """Classify ONE member-bound claim without changing its verification state."""
    if text(row.get(PAPER_FIELD)):
        papers = bound_paper_rows(row)
        if papers:
            return min((_classify_initial(paper) for paper in papers),
                       key=lambda result: (RELEVANCE_TIERS.index(result[0]), result[1:]))
        if has_observations(row):
            return _classify_initial(row)
        # A DB annotation can remain T5, but cannot borrow a paper's T1-T4.
        from .source_verify import classify_citation
        source = _first(row, "source_id", "source_ids")
        annotation = classify_citation(source)[0] == "database_accession"
        if not annotation:
            identity = any(text(row.get(key)) for key in _IDENTITY_FIELDS)
            return ((NEW_TIERS[5], "paper_report_unresolved", "No applicable paper report is bound to this evidence.")
                    if identity else (UNRESOLVED, "identity_unresolved", "No bound paper or usable identity."))
        fallback = {**row, "initial_tier": "T5",
                    "initial_tier_reason": "Source-backed annotation; no bound paper report.",
                    "reported_phenotype": "", "reported_process": "", "phenotype_comparison": "unknown"}
        return _classify_initial(fallback)
    if has_observations(row) or row.get("candidate_admission_policy") == "initial_bibliography_v1" or row.get("relevance_policy") == INITIAL_TIER_POLICY:
        return _classify_initial(row)
    identity = any(text(row.get(key)) for key in _IDENTITY_FIELDS)
    fallback = (T4, "identity_without_trait_bridge", "Identity/function handle only; trait relevance is not established.") if identity else (
        UNRESOLVED, "identity_unresolved", "No usable identity/function handle or trait bridge was supplied.")
    if _first(row, "verification_status", "source_verification_state") == "contradicted":
        return fallback
    if row.get("projection_state") and row.get("projection_state") != "inside_current_qtl":
        return fallback
    if _known_same_trait(row):
        return T1, "same_species_bound_target_trait_report", "One source/identity/trait-bound same-species record maps to this interval member."
    relation = _first(row, "phenotype_relation", "trait_relevance_reason", "expected_role")
    mechanism = _first(row, "mechanism", "mechanism_statement", "mechanistic_link")
    function = _first(row, "biological_function", "function_description", "model_species_evidence")
    source = _first(row, "source_id", "source_ids")
    direct = row.get("trait_relation_state") in {"direct_target_trait", "direct_related_trait"}
    if (direct and text(row.get("reported_phenotype")) and text(row.get("query_phenotype"))
            and source and identity and row.get("evidence_directness") in {
                "D0_causal_perturbation_or_rescue", "D1_mapped_gene_with_functional_validation",
                "D2_direct_trait_gene_association"}):
        return T2, "reported_trait_bridge", "Reported trait relation retained as a candidate bridge; not transferred focal-gene support."
    raw_claim_type = str(row.get("claim_type") or "").strip()
    claim_type = text(raw_claim_type)
    if relation and mechanism and identity and (claim_type in _STRONG_BRIDGES or not raw_claim_type):
        return T2, "explicit_mechanism_bridge", f"Mechanism: {mechanism}; trait relation: {relation}"
    if relation and identity and (mechanism or function):
        return T3, "explicit_functional_plausibility", f"Functional rationale: {mechanism or function}; trait relation: {relation}"
    return fallback


def _classify_initial(row: Mapping[str, Any]) -> tuple[str, str, str]:
    if has_observations(row) or row.get("tier_classification_policy") == JUDGMENT_POLICY:
        label, basis, reason, _ = classify_observation(row)
        return (UNRESOLVED if label == "U" else NEW_TIERS[int(label[1:]) - 1]), basis, reason
    from .bibliography import concrete_bridge
    identity = any(text(row.get(key)) for key in _IDENTITY_FIELDS)
    fallback = (NEW_TIERS[5], "identity_without_trait_bridge", "Identity only; no specific phenotype bridge.") if identity else (
        UNRESOLVED, "identity_unresolved", "Identity or functional interpretation unresolved.")
    if _first(row, "verification_status", "source_verification_state") == "contradicted":
        return fallback
    if row.get("projection_state") and row.get("projection_state") != "inside_current_qtl":
        return fallback
    if row.get("tier_aggregation") != PAPER_POLICY and _known_same_trait(row):
        return NEW_TIERS[0], "bound_focal_trait_report", "Strict bound same-species target-trait record."
    if row.get("tier_aggregation") == PAPER_POLICY:
        proposed = text(row.get("initial_tier")).split("_", 1)[0]
        if proposed == "U":
            return UNRESOLVED, "paper_report_unresolved", text(row.get("initial_tier_reason"))
        if proposed == "T6":
            return fallback
    if not concrete_bridge(dict(row)):
        return fallback
    proposed = text(row.get("initial_tier"))
    if proposed in NEW_TIERS:
        proposed = proposed.split("_", 1)[0]
    reason = text(row.get("initial_tier_reason"))
    reported = text(row.get("reported_phenotype"))
    process = text(row.get("reported_process"))
    comparison = text(row.get("phenotype_comparison"))
    source = _first(row, "source_id", "source_ids", "function_source_ids")
    if proposed == "T1" and reason and reported and source and comparison == "same_equivalent":
        focal, reference = text(row.get("focal_species")), text(row.get("reference_species"))
        if focal and focal.casefold() == reference.casefold():
            return NEW_TIERS[0], "original_reported_focal_trait", reason
    if proposed == "T2" and reason and reported and source and comparison in {"same_equivalent", "same_cellular_process"}:
        return NEW_TIERS[1], "original_reported_similar_phenomenon", reason
    if proposed == "T3" and reason and reported and process and source and comparison == "shared_process":
        return NEW_TIERS[2], "original_reported_shared_process", reason
    if proposed == "T5" and reason and source:
        return NEW_TIERS[4], "original_annotation_hypothesis", reason
    # An underspecified/missing initial Tier never blocks an otherwise eligible
    # candidate. Functional plausibility alone does not become a phenotype report.
    if source:
        return NEW_TIERS[3], "original_function_hypothesis", reason or "Specific functional hypothesis; reported-phenotype Tier unresolved."
    return fallback


def _support_state(row: Mapping[str, Any]) -> str:
    status = text(row.get("support_status"))
    legacy = text(row.get("support_tier"))
    if status in _VERIFIED_SUPPORT or legacy == "tier1_verified_support":
        return "verified"
    if status.startswith("unverified_") or "review" in status or legacy == "unverified_route_review":
        return "review"
    if text(row.get("support_state")) in {"verified", "review", "context", "not_evaluated"}:
        return str(row["support_state"])
    return "context" if any(text(row.get(key)) for key in _IDENTITY_FIELDS) else "not_evaluated"


def _identity_level(row: Mapping[str, Any]) -> str:
    if text(row.get("identity_level")):
        return text(row["identity_level"])
    for level, fields in (
        ("subfamily_or_orthogroup", ("subfamily_id", "subfamily_or_clade", "subfamilies", "orthogroup", "orthogroup_id", "orthogroup_ids")),
        ("family", ("family_id", "family_ids", "evolutionary_family_id")),
        ("reference_homolog", ("reference_gene", "reference_genes", "orthology_reference_gene", "hit_accession", "identity_accessions", "source_reference_gene_id")),
        ("domain", ("domain_ids", "domain_architecture")),
    ):
        if any(text(row.get(field)) for field in fields):
            return level
    return "unresolved"


def annotate_relevance(
    candidates: Iterable[Mapping[str, Any]], claims: Iterable[Mapping[str, Any]] = (),
    *, policy: str = "",
) -> list[dict[str, Any]]:
    """Annotate a roster, preserving each copy and every member-bound claim.

    Claims use candidate_gene_id/gene_id or the frozen member_interval_genes
    roster. A row lacking a member is an error, never global evidence. Unknown
    claim members are ignored (for a supported/review subset), never introduced.
    Existing relevance labels are not trusted as proof. Legacy labels survive
    in their own column and do not become relevance tiers by name substitution.
    """
    roster = [dict(row) for row in candidates]
    ids = [_first(row, "candidate_gene_id", "gene_id") for row in roster]
    if any(not gene for gene in ids) or len(set(ids)) != len(ids):
        raise ValueError("relevance roster requires non-empty unique gene IDs")
    records: dict[str, dict[str, dict[str, Any]]] = {gene: {} for gene in ids}
    accumulated = list(claims)
    for gene, candidate in zip(ids, roster):
        saved = candidate.get("relevance_evidence_json", "")
        if not text(saved):
            continue
        saved = json.loads(str(saved))
        if not isinstance(saved, list) or any(not isinstance(row, dict) for row in saved):
            raise ValueError("relevance evidence must be an object list")
        for row in saved:
            members = _tokens(_first(row, "candidate_gene_id", "gene_id", "member_interval_genes"))
            if gene not in members:
                raise ValueError("saved relevance evidence does not bind its candidate")
        accumulated.extend(saved)
    expanded = []
    for raw in accumulated:
        saved = text(raw.get("relevance_evidence_json"))
        if saved:
            nested = json.loads(saved)
            if not isinstance(nested, list) or any(not isinstance(row, dict) for row in nested):
                raise ValueError("relevance evidence must be an object list")
            owner = _tokens(_first(raw, "candidate_gene_id", "gene_id", "member_interval_genes"))
            for row in nested:
                if not owner & _tokens(_first(row, "candidate_gene_id", "gene_id", "member_interval_genes")):
                    raise ValueError("saved relevance evidence does not bind its candidate")
                if raw.get("candidate_admission_policy") == "initial_bibliography_v1":
                    # A unit claim may mention other references/copies. The
                    # candidate freeze has already limited this application;
                    # rehydration must not re-expand it into a T1 rescue.
                    application = owner & _tokens(_first(row, "candidate_gene_id", "gene_id", "member_interval_genes"))
                    row = {**row, "candidate_gene_id": ";".join(sorted(application)),
                           **{key: raw[key] for key in ("candidate_admission_policy", "candidate_admission_state") if raw.get(key)}}
                expanded.append(row)
        else:
            expanded.append(raw)
    for raw in expanded:
        row = {str(key): value for key, value in dict(raw).items()}
        members = _tokens(_first(row, "candidate_gene_id", "gene_id", "member_interval_genes"))
        if not members:
            raise ValueError("relevance claim requires a bound interval member roster")
        encoded = _json(row)
        for gene in members & set(records):
            records[gene][encoded] = row
    result = []
    initial_policy = policy == INITIAL_TIER_POLICY or any(has_observations(row) or text(row.get(PAPER_FIELD)) or row.get("candidate_admission_policy") == "initial_bibliography_v1"
                         or row.get("relevance_policy") == INITIAL_TIER_POLICY for row in expanded)
    for gene, candidate in sorted(zip(ids, roster)):
        evidence = [records[gene][key] for key in sorted(records[gene])]
        # Candidate metadata is identity/context, never an implicit claim. Only
        # evidence records are allowed to supply a target-trait bridge.
        fallback = {key: candidate.get(key, "") for key in _IDENTITY_FIELDS}
        if initial_policy or candidate.get("candidate_admission_policy") == "initial_bibliography_v1":
            fallback["relevance_policy"] = INITIAL_TIER_POLICY
        paper_policy = any(text(row.get(PAPER_FIELD)) for row in evidence)
        classified = []
        observations = []
        for row in evidence:
            for part in (bound_paper_rows(row) or [row]):
                if has_observations(part) or part.get("tier_classification_policy") == JUDGMENT_POLICY:
                    _, _, _, missing = classify_observation(part)
                    part = {**part, "tier_classification_policy": part.get("tier_classification_policy") or OBSERVATION_POLICY,
                            "tier_classification_state": "partial" if missing else "classified",
                            "tier_missing_fields_json": _json(missing)}
                    observations.append(part)
                classified.append((classify_relevance(part), part))
        classified.append((classify_relevance(fallback), fallback))
        best = min(RELEVANCE_TIERS.index(item[0][0]) for item in classified)
        chosen = [item[0] for item in classified if RELEVANCE_TIERS.index(item[0][0]) == best]
        provenance = [candidate, *evidence]
        old = _first(candidate, "legacy_support_tier", "support_tier")
        candidate.update({
            **{key: _join(row.get(key, "") for row in evidence) for key in INITIAL_FIELDS if key in RELEVANCE_COLUMNS},
            **{key: _join(row.get(key, "") for row in provenance) for key in
               ("candidate_admission_policy", "candidate_admission_state")},
            "relevance_policy": INITIAL_TIER_POLICY if (initial_policy or fallback.get("relevance_policy")) else "legacy_relevance_v1",
            "relevance_tier": RELEVANCE_TIERS[best],
            "relevance_basis": _join(item[1] for item in chosen),
            "tier_reason": " | ".join(sorted({item[2] for item in chosen})),
            "support_state": _support_state(candidate),
            "mapping_tier": _join(row.get("mapping_tier", "") for row in provenance),
            "identity_level": _identity_level(candidate),
            "subfamily_or_orthogroup": _join(_first(row, "subfamily_or_orthogroup", "subfamily_or_clade", "subfamily_id", "subfamilies", "orthogroup_id", "orthogroup_ids", "orthogroup") for row in provenance),
            "copy_state": _join(_first(row, "copy_state", "copy_discrimination_state", "projection_ambiguity_state") for row in provenance),
            "route_provenance": _join(_first(row, "route_provenance", "route", "review_routes", "route_support") for row in provenance),
            "family_id": _join(_first(row, "family_id", "family_ids", "evolutionary_family_id") for row in provenance),
            "reference_gene": _join(_first(row, "reference_gene", "reference_genes", "reference_gene_id", "source_reference_gene_id", "orthology_reference_gene", "identity_accessions") for row in provenance),
            "reference_species": _join(row.get("reference_species", "") for row in provenance),
            "source_ids": _join(_first(row, "source_ids", "source_id") for row in provenance),
            "legacy_support_tier": old,
            "relevance_migration_state": "from_bound_claims" if evidence else "legacy_identity_only_no_trait_proof",
            "relevance_evidence_json": _json(evidence),
        })
        if paper_policy or observations:
            # One deterministic representative for display. All other reports
            # stay in the evidence ledger; multiplicity never affects Tier.
            winners = [row for tier, row in classified if RELEVANCE_TIERS.index(tier[0]) == best]
            representative = min(winners, key=lambda row: (
                not bool(row.get("tier_aggregation") == PAPER_POLICY),
                row.get("tier_classification_state") == "partial",
                _first(row, "source_id", "source_ids"), _json({k: v for k, v in row.items() if k != RESEARCH_ANNOTATIONS})))
            representative_tier = classify_relevance(representative)
            paper_counts, paper_summaries = paper_summary(classified)
            display_evidence = {key: representative.get(key, "") for key in ("source_id", *REPORT_FIELDS,
                                "homolog_projection_json", "proxy_reference_species")}
            if has_observations(representative) or representative.get("tier_classification_policy") == JUDGMENT_POLICY:
                display_evidence.update(initial_tier=representative_tier[0].split("_", 1)[0], initial_tier_reason=representative_tier[2],
                                        **{key: representative[key] for key in
                                           ("tier_classification_policy", "tier_classification_state", "tier_missing_fields_json", "tier_gene_binding")
                                           if key in representative})
            candidate.update(tier_aggregation=PAPER_POLICY,
                             paper_tier_counts_json=_json(paper_counts),
                             paper_evidence_summary_json=_json(paper_summaries),
                             tier_source_id=_first(representative, "source_id", "source_ids"),
                             tier_evidence_json=_json(display_evidence)
                             if representative.get("tier_aggregation") == PAPER_POLICY else "{}",
                             initial_tier=representative_tier[0].split("_", 1)[0],
                             initial_tier_reason=representative_tier[2],
                             tier_reason=representative_tier[2])
            if observations:
                selected = [row for tier, row in classified if RELEVANCE_TIERS.index(tier[0]) == best
                            and row.get("tier_classification_policy") in {OBSERVATION_POLICY, JUDGMENT_POLICY, CHECKLIST_POLICY}]
                status_row = min(selected, key=lambda row: (row.get("tier_classification_state") == "partial",
                                      _json({k: v for k, v in row.items() if k != RESEARCH_ANNOTATIONS}))) if selected else None
                unresolved = (json.loads(status_row["tier_missing_fields_json"]) if status_row else
                              sorted({field for row in observations for field in json.loads(row["tier_missing_fields_json"])}))
                candidate.update(tier_classification_policy=(status_row or {}).get("tier_classification_policy", observations[0]["tier_classification_policy"]),
                                 tier_classification_state="partial" if unresolved else "classified",
                                 tier_missing_fields_json=_json(unresolved),
                                 tier_gene_binding=(status_row or {}).get("tier_gene_binding", ""))
        result.append(candidate)
    return result
