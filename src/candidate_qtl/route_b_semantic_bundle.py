"""Deterministic Route-B semantic-investigation bundles for issue #403.

The interval atlas is built and frozen before phenotype-aware semantic review.
This module takes that frozen B0 inventory plus optional reverse-literature link
candidates and compresses the interval into auditable semantic units.  One unit
becomes one provider-neutral prompt; no LLM is called here.

Grouping is deliberately identity-preserving:

* evolutionary_family_id when available;
* otherwise orthogroup;
* otherwise legacy/coarse family_id;
* otherwise the individual interval gene.

The original family / orthogroup / gene identities remain separate fields in the
bundle.  Unclassified genes are never dropped merely because no family could be
assigned.
"""

from __future__ import annotations

import hashlib
import json
import re
from copy import deepcopy

import pandas as pd

from .core import file_fingerprint, runtime_path
from .route_b_reference_materials import load_reference_materials
from .paper_reports import LITERATURE_RESEARCH_POLICY

SCHEMA_VERSION = "1.18"
INITIAL_MATERIAL_PROJECTION = "similarity_annotation_context_v1"
B0_COVERAGE_STATUSES = ("assigned", "unassigned", "unresolved", "failed")

# IDs and bare similarity/source links are retrieval handles, not reviewable
# function. Classify supplied material before scheduling; keep every unit and
# record skipped review as not evaluated, never as a biological negative.
REQUEST_STATE_REQUESTED = "requested"
REQUEST_STATE_SKIPPED = "not_requested_no_reviewable_material"
_FUNCTION_MATERIAL_FIELDS = (
    "reference_function", "reference_mutant_phenotype", "reference_pathway",
    "function_description", "mutant_phenotype", "pathway",
)
_GENE_MATERIAL_FIELDS = (
    # Functional annotation, including named family/domain context. Arbitrary
    # family/orthogroup/reference IDs and citation rosters are excluded.
    "protein_family", "family_label", "keywords", "go_terms",
    "protein_name", "reference_domain_ids", "reference_domain_description",
    "domain_ids", "domain_architecture", "broad_function_bins", "primary_broad_function",
    "local_annotation_description", "local_annotation_go_terms",
)
_LOOKUP_CONTEXT_FIELDS = (
    "reference_gene", "annotation_origin", "reference_function",
    "reference_pathway", "reference_pathway_id",
    "reference_mutant_phenotype", "reference_source_id", "reference_source_url", "source_url",
    "literature_source_ids", "phenotype_source_ids", "verified_source_ids",
    "evidence_codes", "evidence_type", "provenance_tier", "provenance_detail",
    "source_database", "source_release", "retrieved_at", "corpus_version",
    "tissue_specificity", "developmental_stage", "protein_family", "go_terms", "keywords",
    "reference_taxon_id", "reference_gene_symbol", "reference_family_id",
    "reference_orthogroup_id", "reference_species", "reference_accession_namespace",
    "protein_name", "reference_domain_ids", "reference_domain_description",
    "mapping_state", "mapping_relationship", "function_transfer_scope", "lookup_status",
)
_ATLAS_CONTEXT_FIELDS = _LOOKUP_CONTEXT_FIELDS + (
    "gene_id", "family_id", "evolutionary_family_id", "subfamily_id", "orthogroup",
    "orthology_reference_gene", "orthology_reference_species", "orthology_reference_accession_namespace",
)
# Placeholders that occupy a field without saying anything. `unknown` is what
# B0 writes when it could not classify a gene, so treating it as material would
# make every unit look reviewable and the skip a no-op (caught by the 200-gene
# bench, which showed a 0% reduction before this was excluded).
_EMPTY_TOKENS = frozenset({
    "", "unknown", "none", "null", "na", "n/a", "unclassified", "not_applicable",
    "not_evaluated", "unresolved", "false", "0", "nan", "-", "...", "…", "[truncated]",
})


def _has_material(value: object) -> bool:
    from .route_b_identity import annotation_material

    return any(text.strip().casefold() not in _EMPTY_TOKENS and annotation_material(text)
               for text in str(value or "").replace(";", "\n").splitlines())


def similarity_annotation_context(hit: dict[str, object]) -> dict[str, object]:
    """Return original informative hit annotations, never verified function.

    Header identifiers and attribution are stripped only for classification;
    delivered values and their endpoint/source binding remain untouched.
    """
    from .route_b_identity import header_metadata

    context = {}
    identifier = r"(?:[A-Za-z]{1,8}[_:.|-]*)?\d[\w.:|-]*|\S*\|\S*"
    endpoints = {str(hit.get(key, "")).strip() for key in ("hit_accession", "raw_reference_gene_id")}
    for field in ("hit_protein_name", "hit_description"):
        raw = hit.get(field, "")
        text = str(raw or "").strip()
        if field == "hit_description":
            parsed = header_metadata(text)
            if parsed["hit_protein_name"] != "unknown":
                text = parsed["hit_protein_name"]
            else:
                text = re.split(r"\s(?:OS|OX|GN|PE|SV)=", text, maxsplit=1)[0]
            # Species and database-attribution brackets alone are not annotation.
            text = re.sub(r"\s*\[[^\[\]]*\]\s*$", "", text).strip()
            parts = text.split(None, 1)
            if parts and (parts[0] in endpoints or re.fullmatch(identifier, parts[0])):
                text = parts[1] if len(parts) == 2 else ""
        if (text not in endpoints and text.casefold() != "hypothetical" and _has_material(text)
                and not re.fullmatch(identifier, text)):
            context[field] = raw
    return context


def _investigation_material(unit: dict[str, object], *, include_similarity_annotations: bool = False) -> dict[str, object]:
    """Classify available context without inferring trait support from IDs."""
    records = list(unit["genes"]) + list(unit["reference_material_cards"]) + list(unit["literature_link_candidates"])
    records.extend(row["atlas_material_context"] for row in unit["genes"] if "atlas_material_context" in row)
    functional_count = sum(any(_has_material(row.get(field, "")) for field in _FUNCTION_MATERIAL_FIELDS) for row in records)
    annotation_count = sum(any(
        _has_material(row.get(field, ""))
        and not (field == "family_label" and row.get(field) in {row.get(key) for key in ("family_id", "evolutionary_family_id", "subfamily_id")})
        for field in _GENE_MATERIAL_FIELDS
    ) for row in records)
    if include_similarity_annotations:
        annotation_count += sum(bool(similarity_annotation_context(hit))
                                for gene in unit["genes"] for hit in gene.get("sequence_similarity_hits", []))
    handles = bool(unit["stable_source_ids"] or unit["literature_link_candidates"]
                   or any(unit["identity"].values()) or unit["reference_material_cards"]
                   or any(row.get("reference_gene") or row.get("orthology_reference_gene") or row.get("sequence_similarity_hits")
                          for row in unit["genes"]))
    reviewable = bool(functional_count or annotation_count)
    return {
        "annotation_record_count": annotation_count,
        "functional_record_count": functional_count,
        "source_backed_function_record_count": unit["source_backed_function_record_count"],
        "literature_link_count": len(unit["literature_link_candidates"]),
        "state": ("functional_or_literature_material" if functional_count else "annotation_only" if annotation_count else "needs_material_retrieval"),
        "source_material_class": ("function_context" if functional_count else "annotation_context" if annotation_count
                                  else "retrieval_handles_only" if handles else "no_reviewable_material"),
        "semantic_review_recommended": reviewable,
        "policy": "function and annotation are investigation context, not causal support; retrieve material before spending model calls on bare identity/hit/source IDs",
    }


def _unit_request_state(unit: dict[str, object]) -> str:
    """Use the material ledger, not a second identity-based fallback."""
    return REQUEST_STATE_REQUESTED if unit["investigation_material"]["semantic_review_recommended"] else REQUEST_STATE_SKIPPED


def _classify_gene_b0_coverage(row: dict[str, object] | pd.Series) -> str:
    """Classify an interval gene's B0 coverage into assigned/unassigned/unresolved/failed."""
    evo_status = str(row.get("evolutionary_family_status", "")).strip().casefold()
    if evo_status == "failed":
        return "failed"
    if evo_status == "unresolved" or str(row.get("conservation_class", "")).strip().casefold() == "unresolved":
        return "unresolved"
    notes = str(row.get("notes", "")).strip().casefold()
    if "failed" in notes:
        return "failed"
    if "conflict" in notes or str(row.get("evidence_completeness", "")).strip().casefold() == "conflicting_annotation":
        return "unresolved"
    if (
        evo_status == "assigned"
        or str(row.get("evolutionary_family_id", "")).strip()
        or str(row.get("orthogroup", "")).strip()
        or str(row.get("family_id", "")).strip()
        or str(row.get("domain_architecture", "")).strip()
        or str(row.get("orthology_reference_gene", "")).strip()
    ):
        return "assigned"
    return "unassigned"


IDENTITY_SIGNAL_COLUMNS = (
    "evolutionary_family_id",
    "orthogroup",
    "family_id",
    "domain_architecture",
    "orthology_reference_gene",
)


def _identity_signal_counts(atlas: pd.DataFrame) -> dict[str, int]:
    """Count interval genes carrying each B0 identity signal.

    ``_classify_gene_b0_coverage`` short-circuits on an unresolved
    ``evolutionary_family_status`` or ``conservation_class`` before it looks at
    the identity columns, so a run where every gene is evolutionarily
    unresolved reports ``assigned: 0`` even when most genes do carry a
    ``family_id``.  That count is correct for what it measures and misleading
    as a summary, so publish the per-signal counts next to it and let the
    reader see which layer is actually missing.
    """
    counts: dict[str, int] = {}
    for column in IDENTITY_SIGNAL_COLUMNS:
        if column not in atlas.columns:
            counts[column] = 0
            continue
        counts[column] = int(atlas[column].astype(str).str.strip().ne("").sum())
    return counts

SYSTEM_PROMPT = """You are the semantic investigation stage of an interval-first plant-QTL workflow.
The interval inventory was frozen before this phenotype-aware review. This is the interpretation step of a Web-enabled investigation: the research agent may search primary literature by reference gene, homolog, family, process and phenotype, then freeze the retrieved source-link candidates in a new bundle before this step. Use only the supplied unit and source-link candidates for this frozen response; this source boundary is not a ban on upstream Web research.
In ordinary investigations, the Route B research agent may use available general Web search, paper search and source retrieval from the initial assessment and throughout reasoning. Move between annotation, homolog function, model knowledge, searches and source reading as useful; do not impose a search-free preliminary gate, one-paper limit or uniform per-gene search quota. Supplied-material-only benchmark assignments are a separate explicit execution mode, not the ordinary research policy.
When research adds a source or changes the material used for a claim, return the retrieved material and its gene/species/source binding to the coordinator for a new frozen bundle before submitting that claim. Preserve the earlier bundle and responses. Model knowledge and search snippets are leads, not verified paper evidence. Record unavailable tools, failed retrieval and inaccessible sources separately from a completed search with no relevant result; none proves biological absence. The offline verifier consumes the supplied source cache and does not perform the research agent's searches.
Start with the supplied interval genes/copies and their reference endpoints: determine the actually annotated or reported homolog/family functions, perturbation phenotypes, biological processes and evidence limitations before comparing them with the observed focal phenotype. Do not presume that every unit must produce a candidate claim.
Your primary task is a cross-species biological bridge, assessed neutrally: does the actual reference-gene/family biology, together with the supported homolog/copy correspondence, provide a specific plausible explanation of the observed focal phenotype? Consider evidence that the proposed link does not fit, not only ways it might fit; do not start by advocating a candidate or inventing a mechanism for every gene. Same-species literature is supplementary evidence that may strengthen a candidate's rationale; it must not replace this main approach with retrieval of known focal causal genes. A defensible cross-species hypothesis may become a review candidate without an experiment in the focal species or an identical reference phenotype. Apply this reasoning to the actual case, without species-specific or example-specific selection rules or a required model species.
State the chain explicitly: reported reference observation -> homolog/family correspondence -> proposed focal mechanism, with organ, developmental stage, effect direction and competing explanations where known. Homology grounds the correspondence; the conditional mechanistic inference is your reasoning, not a reported experimental fact. Do not reject a specific bridge merely because focal trait transfer remains unverified.
Assess biological specificity, not mere logical possibility: generic kinase or transcription-factor involvement is not a specific bridge, and a chain of unobserved expression changes, regulatory partners or perturbation directions must not be invented merely to connect an unrelated function to the phenotype. A multi-step but specifically grounded cross-species mechanism remains eligible for review when its uncertain steps are stated as assumptions; neither a one-step pathway nor a published focal causal gene is required.
A retrieved paper is only a link candidate, not verified evidence. Never invent a citation or source identifier.
When a supplied curated reference-material card supports the stated reference-function observation, preserve its source_id in source_ids (for example a UniProt record); a focal-species causal paper is not required to cite that reference function. Read the actual function content and its provenance before citing it. A list of literature_source_ids is a set of retrieval leads, not proof that those papers were read or support the observation. Do not silently drop the curated function's source merely because no phenotype paper was retrieved; if its content or identity cannot be checked, record that limitation and keep the claim unresolved.
reference_gene must exactly match a supplied unit reference endpoint; do not substitute a paper's gene symbol for that identifier. Keep unresolved identity mappings explicit rather than guessing an accession.
atlas_material_context, when present, is a separate original atlas record with its own origin fingerprint. Do not attach its function, species or citations to the containing gene's lookup reference; any correspondence or trait transfer remains unverified.
Sequence-similarity hits, when supplied, are complete retrieval-handle alternatives for each interval gene. They are not confirmed orthologs, function annotations, or trait evidence. Do not turn their alignment metrics into a causal score; use their accession/species/symbol as search handles and keep any biological bridge conditional on retrieved source evidence.
In your analysis, do not rank genes numerically, do not force one copy, and do not convert missing data into negative evidence.
Return JSON only as {\"claims\":[{\"claim_type\":\"direct_functional_bridge|mechanistic_process_bridge|mutant_phenotype_bridge|mechanistic_review_seed\",\"biological_function\":\"...\",\"mechanism\":\"...\",\"phenotype_relation\":\"...\",\"reference_gene\":\"exact supplied reference endpoint\",\"reference_species\":\"...\",\"evidence_scope\":\"reference_gene\",\"source_ids\":[\"only supplied source ids\"],\"alternative_hypothesis\":\"...\",\"uncertainty_note\":\"...\"}]}.
The JSON example uses reference_gene scope, which requires a non-empty reference_gene matching a supplied unit reference endpoint. Choose one evidence_scope from family, subfamily, orthogroup, reference_gene, interval_gene, process, unknown according to what the stated observation or hypothesis concerns; the example is not a default scope for every claim.
Scope anchors are mandatory: family requires a non-empty family_id derived from unit.identity.evolutionary_family_ids, or unit.identity.family_ids when evolutionary_family_ids is empty; orthogroup requires orthogroup_id derived from unit.identity.orthogroup_ids; interval_gene requires member_interval_genes from the frozen unit. These fields are filled by the importer from the unit, not created by adding IDs to your response.
subfamily transfer verification requires subfamily_or_clade from unit.identity.subfamily_ids. process and unknown have no additional scope-specific ID requirement, but every claim still needs a supplied biological anchor (family_id, orthogroup_id, member_interval_genes or reference_gene). A scope label does not itself establish source support or trait transfer.
If the unit has no family IDs, do not emit family-scoped claims or invent a family_id from a protein name, domain, reference label or your own inference. Reconsider the actual scope and supplied anchors: a reference observation may use reference_gene only with its exact supplied endpoint; a claim about a focal interval gene may use interval_gene only when its content concerns that gene. Do not relabel evidence merely to bypass a missing anchor; preserve limitations and return an empty claims list if no meaningful anchored bridge can be stated.
Use direct_functional_bridge only when the supplied curated function directly supports the trait relation, mechanistic_process_bridge only when the supplied process/pathway supports a concrete mechanism, and mutant_phenotype_bridge only when a supplied mutant/perturbation phenotype supports the relation.
Use mechanistic_review_seed when the supplied family, domain, annotation or reference function suggests a specific phenotype-directed mechanism worth preserving for later review, but the supplied material does not directly support the trait relation. State the proposed mechanism and why it could connect to the focal phenotype, keep source_ids empty unless a supplied source supports the stated observation, and explicitly say that trait transfer/support is unverified. Do not use this type for a generic protein class, housekeeping function, or mere statement that relevance is unknown.
mechanistic_review_seed is a review candidate only: it must not be described as verified support, must not select one copy, and must not receive a numerical score. If no specific bridge can be stated from the supplied unit, return an empty claims list.
For newly frozen source_function_v1 bundles, a raw seed with missing sources remains in the investigation ledger. B-only Wet review admission requires checked external reference-function content, frozen correspondence and a specific conditional phenotype bridge; it never requires a focal causal paper or an identical reference phenotype. Do not fabricate sources to gain admission; A/C-independent eligibility is unaffected.
A reviewed unit may legitimately return an empty claims list when the examined reference functions do not currently offer a concrete focal-phenotype bridge. This is a decision not to propose a claim from the current material, not biological evidence against the gene. Missing annotation, inaccessible literature, unsearched alternative homologs and provider failure remain unresolved/not_evaluated rather than negative. Do not turn missing evidence or an ambiguous copy relationship into an exclusion.
Except for mechanistic_review_seed, keep each claim to one independently source-supported observation or process. Split observations involving different experimental genes, organisms, or phenotype components into separate claims; do not combine them into a stronger conjunction.
Each source_ids entry must support that claim's stated observation, not merely another part of a multi-paper story. A genome-resource citation or homology alone does not establish a functional experiment.
Use biological_function and mechanism for the reported function/process. In phenotype_relation distinguish what was observed from the conditional relevance to the focal phenotype; put untested extensions and transfer limitations in uncertainty_note or alternative_hypothesis.
A relevant source-backed component need not explain every component of the focal phenotype. Do not describe a plausible extension as an experimentally demonstrated phenotype, or claim that a focal interval copy is causal.
Do not emit a claim merely to restate a generic annotation or to say that trait relevance is unknown/not established; return {\"claims\":[]} for that unit.
If the supplied material does not support a meaningful bridge, return {\"claims\":[]} rather than manufacturing one.
All returned claims remain proposals until source/identity checks and semantic review; mechanical checks do not prove biological meaning. Review membership and verified trait support are separate outcomes."""

from .relevance_tiers import INITIAL_TIER_PROMPT
CONTENT_SYSTEM_PROMPT = SYSTEM_PROMPT
SYSTEM_PROMPT = "\n".join(line for line in SYSTEM_PROMPT.splitlines() if not line.startswith((
    "mechanistic_review_seed is a review candidate only:",
    "For newly frozen source_function_v1 bundles,",
    "All returned claims remain proposals until source/identity checks and semantic review;",
))) + "\n" + (
    "Your original biological judgment supplies first-stage candidate admission. "
    "A concrete reference-function bridge with resolved bibliography and frozen correspondence "
    "is adopted as a first-stage candidate; a raw source-free seed stays in the investigation ledger. "
    "Do not fabricate sources or transfer reported observations to a different copy. "
    "A/C-independent eligibility is unaffected. Mechanical checks do not establish causal biology. "
    "First-stage admission and verified trait support are separate outcomes. "
    "If no specific bridge can be stated from the supplied unit, return an empty claims list.\n"
) + (
    "Put paper_reports_json inside each claim: one report per paper read, or []. Never invent citations or return placeholders. "
    "reported_gene/reported_gene_id name the paper's gene; reference_gene_id is the exact selected claim.reference_gene, "
    "without focal-gene/profile labels. member_interval_genes lists focal copies. Keep unchecked answers and unresolved "
    "correspondence unknown with reasons; do not change answers to match. Missing reports do not veto claims.\n"
) + INITIAL_TIER_PROMPT

from .paper_reports import LITERATURE_SEARCH_GUIDANCE, render_initial_paper_report_example
# Only the current initial_bibliography_v1 prompt gains the concrete per-claim
# paper_reports_json example; CONTENT_SYSTEM_PROMPT (legacy source_function_v1)
# is captured above and intentionally left unchanged.
SYSTEM_PROMPT = SYSTEM_PROMPT.replace(
    "These fields are filled by the importer from the unit, not created by adding IDs to your response.",
    "Identity anchors are filled by the importer from the unit. An optional member_interval_genes declaration may narrow applicability using exact IDs already in that unit.",
)
SYSTEM_PROMPT += LITERATURE_SEARCH_GUIDANCE + "\nOptional claim.member_interval_genes is a non-empty list of exact frozen unit gene IDs (or a semicolon-separated string). Use it to declare the precise claim applicability subset; omit it only when the claim applies to the whole frozen unit. Retain all plausible copies within that scope; do not force a sole preferred copy without copy-specific evidence. Explicit empty or outside-unit membership is invalid.\n"
_B_PAPER_REPORTS_EXAMPLE = render_initial_paper_report_example()
SYSTEM_PROMPT = SYSTEM_PROMPT.replace(
    '"uncertainty_note":"..."}]}',
    '"uncertainty_note":"...","paper_reports_json":' + _B_PAPER_REPORTS_EXAMPLE + "}]}",
)
if '"paper_reports_json"' not in SYSTEM_PROMPT:
    raise RuntimeError("failed to place paper_reports_json in the Route B claim example")

UNIT_COLUMNS = [
    "unit_id",
    "unit_type",
    "unit_anchor",
    "member_count",
    "member_interval_genes",
    "evolutionary_family_ids",
    "orthogroup_ids",
    "family_ids",
    "subfamily_ids",
    "reference_genes",
    "domain_architectures",
    "broad_function_bins",
    "unclassified_member_count",
    # Two independent bridge materials, kept side by side on purpose.  The
    # literature column is filled only by the legacy B1 reverse-literature
    # lane, which is diagnostic-only and offline by default, so it reads zero
    # in a canonical run.  The curated-function column is what canonical Route
    # B actually supplies and what the frozen system prompt keys on.  Showing
    # only the first invites a reader to conclude the unit carries nothing to
    # bridge from when it does.
    "literature_link_count",
    "source_backed_function_record_count",
    "stable_source_count",
    "source_material_class",
    "request_state",
    "review_status",
    "estimated_input_tokens",
    "estimated_model_call_count",
    "rank_contribution",
]


def _tokens(values: list[object] | pd.Series) -> list[str]:
    found: set[str] = set()
    for value in values:
        for token in str(value).replace(",", ";").split(";"):
            token = token.strip()
            if token:
                found.add(token)
    return sorted(found)


def _unit_key(row: pd.Series) -> tuple[str, str]:
    evolutionary = str(row.get("evolutionary_family_id", "")).strip()
    if evolutionary:
        return "family", evolutionary
    orthogroup = str(row.get("orthogroup", "")).strip()
    if orthogroup:
        return "orthogroup", orthogroup
    family = str(row.get("family_id", "")).strip()
    if family:
        return "family", family
    return "interval_gene", str(row["gene_id"]).strip()


def _unit_id(unit_type: str, anchor: str, members: list[str]) -> str:
    payload = json.dumps(
        {"unit_type": unit_type, "anchor": anchor, "members": members},
        sort_keys=True,
        separators=(",", ":"),
    )
    return f"bsem-{hashlib.sha256(payload.encode('utf-8')).hexdigest()[:16]}"


def _stable_source_ids(rows: list[dict[str, object]]) -> list[str]:
    """Build one de-duplicated source roster from literature/corpus fields.

    A PMID repeated as a phenotype source, literature source, and verified source
    is one source, not three votes.  Database accessions remain in the roster
    alongside PMID/PMCID/DOI tokens because the semantic provider needs the
    complete closed world, while the downstream verifier decides what each
    source actually supports.
    """
    found: dict[str, str] = {}
    fields = (
        "source_id", "reference_source_id", "literature_source_ids", "phenotype_source_ids", "verified_source_ids",
        "doi", "pmid", "pmcid",
    )
    for row in rows:
        for field in fields:
            for token in str(row.get(field, "") or "").replace(",", ";").split(";"):
                token = token.strip()
                if token:
                    found.setdefault(token.casefold(), token)
    return sorted(found.values(), key=lambda item: (item.casefold(), item))


def _prompt_record(unit: dict[str, object], phenotype: str, case_id: str,
                   policy: str = "initial_bibliography_v1") -> dict[str, object]:
    """Build the same bounded prompt for estimates and later export."""
    prompt_unit = dict(unit)
    prompt_genes = []
    for gene in unit["genes"]:
        prompt_gene = dict(gene)
        prompt_gene["sequence_similarity_hits"] = [
            {key: hit.get(key, "") for key in (
                "hit_accession", "hit_species", "hit_gene_symbol", "hit_protein_name", "hit_description", "hit_metadata_state",
                "percent_identity", "query_coverage", "interpretation", "rank_contribution",
            )}
            for hit in gene.get("sequence_similarity_hits", [])
        ]
        prompt_genes.append(prompt_gene)
    prompt_unit["genes"] = prompt_genes
    return {
        "prompt_id": f"route-b-semantic-{unit['unit_id']}",
        **({"literature_research": dict(LITERATURE_RESEARCH_POLICY)} if policy == "initial_bibliography_v1" else {}),
        "gene_id": "__route_b_semantic_unit__",
        "system": SYSTEM_PROMPT if policy == "initial_bibliography_v1" else CONTENT_SYSTEM_PROMPT,
        "user": (
            f"Phenotype: {str(phenotype).strip()}\n"
            "Route: interval_first; B0 inventory already frozen.\n"
            "Investigate this semantic unit and return proposal claims only.\n\n"
            + json.dumps(prompt_unit, ensure_ascii=False, sort_keys=True)
        ),
        "route_b_unit_id": unit["unit_id"],
        "case_id": case_id,
    }


def _unit_summary(unit: dict) -> dict:
    """Keep the same navigation-table contract for initial and extended bundles."""
    row = {key: unit[key] for key in ("unit_id", "unit_type", "unit_anchor", "source_backed_function_record_count",
                                     "source_material_class", "request_state", "review_status")}
    row.update(member_count=len(unit["member_interval_genes"]), member_interval_genes=";".join(unit["member_interval_genes"]),
               unclassified_member_count=len(unit["member_interval_genes"]) if unit["unit_type"] == "interval_gene" else 0,
               literature_link_count=len(unit["literature_link_candidates"]), stable_source_count=len(unit["stable_source_ids"]),
               estimated_input_tokens=unit["estimated_review_cost"]["estimated_input_tokens"],
               estimated_model_call_count=unit["estimated_review_cost"]["estimated_model_call_count"], rank_contribution="none")
    row.update({key: ";".join(unit["identity"][key]) for key in ("evolutionary_family_ids", "family_ids", "orthogroup_ids", "subfamily_ids")})
    for column, field in (("reference_genes", "orthology_reference_gene"), ("domain_architectures", "domain_architecture"), ("broad_function_bins", "broad_function_bins")):
        row[column] = ";".join(_tokens([gene.get(field, "") for gene in unit["genes"]]))
    return row


def refreeze_route_b_sources(parent_arg: str, discoveries_arg: str, out_dir_arg: str, *, parent_sha256: str) -> dict:
    """Append explicitly returned Web material to a NEW frozen generation.

    No retrieval, Tier judgment, candidate selection or support promotion runs
    here. Reference endpoints and member rosters remain those of the parent.
    """
    from .source_verify import classify_citation
    parent_path, discovery_path, out = map(runtime_path, (parent_arg, discoveries_arg, out_dir_arg))
    actual_hash = hashlib.sha256(parent_path.read_bytes()).hexdigest()
    if not parent_sha256 or actual_hash != parent_sha256:
        raise ValueError("research parent bundle SHA-256 mismatch")
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError("research re-freeze requires a new or empty output directory")
    parent = json.loads(parent_path.read_text(encoding="utf-8"))
    payload = json.loads(discovery_path.read_text(encoding="utf-8"))
    if (parent.get("stage") != "post_b0_semantic_bundle" or parent.get("route") != "interval_first"
            or not parent.get("case_id") or not parent.get("phenotype")):
        raise ValueError("research requires a frozen post-B0 bundle")
    if payload.get("case_id") != parent["case_id"] or payload.get("parent_bundle_sha256") != actual_hash:
        raise ValueError("research discovery case/parent hash mismatch")
    discoveries = payload.get("discovered_sources")
    if not isinstance(discoveries, list) or not discoveries:
        raise ValueError("research requires non-empty discovered_sources")
    bundle = deepcopy(parent)
    bundle["schema_version"] = SCHEMA_VERSION
    units = {unit["unit_id"]: unit for unit in bundle["units"]}
    changed = set()
    required = ("unit_id", "reference_gene", "reference_species", "reported_gene", "reported_species",
                "source_id", "source_title", "source_url", "retrieved_at", "source_text", "reference_function")
    for discovery in discoveries:
        if not isinstance(discovery, dict) or any(not isinstance(discovery.get(key), str) or not discovery[key].strip() for key in required):
            raise ValueError("research discovery requires source text, title, retrieval and gene/species bindings")
        unit = units.get(discovery["unit_id"])
        if unit is None:
            raise ValueError("research discovery names an unknown unit")
        references = set()
        for record in [*unit["genes"], *unit.get("reference_material_cards", []), *unit.get("literature_link_candidates", [])]:
            for key in ("reference_gene", "orthology_reference_gene", "family_index_subject_endpoints"):
                references.update(_tokens([record.get(key, "")]))
            # A secondary frozen BLAST hit is also an exact retrieval handle.
            # New material may describe it without turning similarity into
            # orthology, changing the unit roster, or authorizing support.
            for hit in record.get("sequence_similarity_hits", []):
                if isinstance(hit, dict) and str(hit.get("hit_accession", "")).strip():
                    references.add(str(hit["hit_accession"]).strip())
        from .reference_name_pairs import reference_handle
        handle = reference_handle(discovery["reference_gene"], allow_untyped=True)
        if not any(discovery["reference_gene"].casefold() == ref.casefold()
                   or (handle is not None and handle == reference_handle(ref, allow_untyped=True)) for ref in references):
            raise ValueError("research discovery names an unfrozen reference endpoint")
        aliases = discovery.get("source_identifier_aliases", "")
        if not isinstance(aliases, str):
            raise ValueError("research source aliases must be a string")
        ids = list(dict.fromkeys([discovery["source_id"], *_tokens([aliases])]))
        if any(classify_citation(identifier)[0] not in {"pmid", "pmcid", "doi"} for identifier in ids):
            raise ValueError("research discovery requires literature identifiers")
        if not discovery["source_url"].startswith(("https://", "http://")):
            raise ValueError("research discovery requires an HTTP source URL")
        material = {key: discovery.get(key, "") for key in (*required, "reported_gene_id", "source_identifier_aliases")}
        material.update(title=discovery["source_title"], query_status="agent_retrieved_unverified",
                        source_ids=ids, rank_contribution="none",
                        source_text_sha256=hashlib.sha256(discovery["source_text"].encode("utf-8")).hexdigest())
        unit["literature_link_candidates"].append(material)
        unit["stable_source_ids"] = sorted(set(unit["stable_source_ids"]) | set(ids))
        unit["source_provenance"].append({**{key: value for key, value in material.items() if key != "source_text"},
                                          "semantic_binding_state": "original_research_claim_not_verified"})
        unit["request_state"] = REQUEST_STATE_REQUESTED
        unit["review_status"] = "not_evaluated"
        unit["not_evaluated_reason"] = "awaiting_semantic_review_after_source_refreeze"
        unit["investigation_material"] = _investigation_material(
            unit, include_similarity_annotations=bundle.get("initial_material_projection") == INITIAL_MATERIAL_PROJECTION)
        unit["source_material_class"] = unit["investigation_material"]["source_material_class"]
        unit.pop("estimated_review_cost", None)
        estimate = _prompt_record(unit, bundle["phenotype"], bundle["case_id"], bundle.get("b_review_membership_policy", "source_function_v1"))
        input_bytes = len(estimate["system"].encode()) + len(estimate["user"].encode())
        unit["estimated_review_cost"] = {"status": "estimate_not_actual_usage", "input_utf8_bytes": input_bytes,
            "estimated_input_tokens": (input_bytes + 3) // 4, "estimated_model_call_count": 1, "actual_input_tokens": None,
            "method": "ceil(system_and_user_utf8_bytes/4); not measured token usage"}
        changed.add(unit["unit_id"])
    policy = bundle.get("b_review_membership_policy", "source_function_v1")
    prompts = [_prompt_record(unit, bundle["phenotype"], bundle["case_id"], policy)
               for unit in bundle["units"] if unit["request_state"] == REQUEST_STATE_REQUESTED]
    bundle["execution_prompt_count"] = len(prompts)
    bundle["review_cost_estimate"] = {"status": "estimate_not_actual_usage", "actual_input_tokens": None,
        "estimated_model_call_count": len(prompts),
        "estimated_input_tokens": sum(unit["estimated_review_cost"]["estimated_input_tokens"] for unit in bundle["units"]
                                      if unit["request_state"] == REQUEST_STATE_REQUESTED)}
    bundle["execution_prompt_fingerprints"] = [hashlib.sha256(json.dumps(prompt, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
                                                for prompt in prompts]
    bundle["research_source_lineage"] = {"parent_bundle_sha256": actual_hash,
        "discovery_sha256": hashlib.sha256(discovery_path.read_bytes()).hexdigest(), "changed_unit_ids": sorted(changed),
        "membership_effect": "none", "source_binding_verified": False}
    out.mkdir(parents=True, exist_ok=True)
    bundle_path, prompts_path = out / "route_b_semantic_bundle.json", out / "route_b_semantic_prompts.jsonl"
    bundle_path.write_text(json.dumps(bundle, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    prompts_path.write_text("".join(json.dumps(prompt, ensure_ascii=False, sort_keys=True) + "\n" for prompt in prompts), encoding="utf-8")
    # Tables are navigation only; the complete frozen unit stays in the bundle.
    units_path, skipped_path = out / "route_b_semantic_units.tsv", out / "route_b_semantic_units_not_requested.tsv"
    table = pd.DataFrame([_unit_summary(unit) for unit in bundle["units"]], columns=UNIT_COLUMNS)
    table.to_csv(units_path, sep="\t", index=False)
    table[table.request_state.ne(REQUEST_STATE_REQUESTED)].to_csv(skipped_path, sep="\t", index=False)
    manifest_path = out / "route_b_semantic_bundle_manifest.json"
    outputs = {"bundle": str(bundle_path), "prompts": str(prompts_path), "units": str(units_path),
               "units_not_requested": str(skipped_path), "manifest": str(manifest_path)}
    manifest = {"schema_version": SCHEMA_VERSION, "case_id": bundle["case_id"], "route": "interval_first",
        "stage": "post_b0_semantic_bundle", "interval_gene_count": bundle["interval_gene_count"],
        "semantic_unit_count": len(units), "prompt_count": len(prompts), "candidate_identity_preserved": True,
        "research_source_lineage": bundle["research_source_lineage"], "outputs": outputs,
        "bundle_sha256": hashlib.sha256(bundle_path.read_bytes()).hexdigest(),
        "prompt_sha256": hashlib.sha256(prompts_path.read_bytes()).hexdigest(),
        "next_action": "prepare new prompts; regenerate changed units; reuse only byte-identical old prompts"}
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {**outputs, "prompt_count": len(prompts), "interval_gene_count": bundle["interval_gene_count"],
            "semantic_unit_count": len(units), "changed_unit_ids": sorted(changed)}


def build_route_b_semantic_bundle(
    atlas_arg: str,
    out_dir_arg: str,
    *,
    phenotype: str,
    literature_arg: str | None = None,
    literature_status_arg: str | None = None,
    lookup_arg: str | None = None,
    identity_hits_arg: str | None = None,
    case_id: str = "",
    full_hash: bool = False,
    b_review_membership_policy: str = "initial_bibliography_v1",
    focal_species: str = "",
) -> dict[str, object]:
    """Build frozen semantic units and provider-neutral prompts from Route B0.

    ``phenotype`` is intentionally introduced only here, after the caller has a
    frozen atlas.  This function does not consume Route-A artifacts or benchmark
    truth and never changes candidate membership.
    """
    if not str(phenotype).strip():
        raise ValueError("phenotype is required for the post-B0 semantic stage")
    if b_review_membership_policy not in {"source_function_v1", "initial_bibliography_v1"}:
        raise ValueError("unsupported Route B membership policy")
    atlas_path = runtime_path(atlas_arg)
    if not atlas_path.exists():
        raise ValueError(f"Route B atlas does not exist: {atlas_path}")
    atlas = pd.read_csv(atlas_path, sep="\t", dtype=str).fillna("")
    if "gene_id" not in atlas.columns:
        raise ValueError("Route B atlas requires gene_id")
    if atlas["gene_id"].astype(str).str.strip().eq("").any() or atlas["gene_id"].duplicated().any():
        raise ValueError("Route B atlas gene_id values must be non-empty and unique")

    literature_path = runtime_path(literature_arg) if literature_arg else None
    if literature_path and not literature_path.exists():
        raise ValueError(f"reverse literature table does not exist: {literature_path}")
    literature = (
        pd.read_csv(literature_path, sep="\t", dtype=str).fillna("")
        if literature_path
        else pd.DataFrame()
    )
    if not literature.empty and "interval_gene_ids" not in literature.columns:
        raise ValueError("reverse literature table requires interval_gene_ids")

    literature_status_path = runtime_path(literature_status_arg) if literature_status_arg else None
    if literature_status_path and not literature_status_path.exists():
        raise ValueError(f"reverse literature status table does not exist: {literature_status_path}")
    literature_status = (
        pd.read_csv(literature_status_path, sep="\t", dtype=str).fillna("")
        if literature_status_path
        else pd.DataFrame()
    )

    lookup_path = runtime_path(lookup_arg) if lookup_arg else None
    if lookup_path and not lookup_path.exists():
        raise ValueError(f"interval reverse lookup table does not exist: {lookup_path}")
    lookup = (
        pd.read_csv(lookup_path, sep="\t", dtype=str).fillna("")
        if lookup_path
        else pd.DataFrame()
    )
    if not lookup.empty:
        if "gene_id" not in lookup.columns:
            raise ValueError("interval reverse lookup table requires gene_id")
        if lookup["gene_id"].astype(str).str.strip().eq("").any() or lookup["gene_id"].duplicated().any():
            raise ValueError("interval reverse lookup gene_id values must be non-empty and unique")
    lookup_by_gene = {
        str(row["gene_id"]): row for row in lookup.to_dict("records")
    } if not lookup.empty else {}

    # The identity layer deliberately retains every threshold-passing BLAST
    # handle. Carry those alternatives into the semantic unit so an upstream
    # best-hit projection cannot hide a plausible family member. These are
    # retrieval handles only: sequence similarity does not establish
    # orthology, function, or trait support.
    identity_hits_path = runtime_path(identity_hits_arg) if identity_hits_arg else None
    if identity_hits_path and not identity_hits_path.exists():
        raise ValueError(f"Route B identity hits table does not exist: {identity_hits_path}")
    identity_hits = (
        pd.read_csv(identity_hits_path, sep="\t", dtype=str).fillna("")
        if identity_hits_path else pd.DataFrame()
    )
    identity_hits_by_gene: dict[str, list[dict[str, object]]] = {}
    if not identity_hits.empty:
        if "focal_gene_id" not in identity_hits.columns:
            raise ValueError("Route B identity hits table requires focal_gene_id")
        unknown_genes = sorted(set(identity_hits["focal_gene_id"].astype(str)) - set(atlas["gene_id"].astype(str)))
        if unknown_genes:
            raise ValueError(f"Route B identity hits contain genes outside the frozen atlas: {unknown_genes[:10]}")
        for gene, rows in identity_hits.groupby("focal_gene_id", sort=True):
            candidates = []
            for row in rows.to_dict("records"):
                candidates.append({
                    key: str(row.get(key, "")) for key in (
                        "hit_accession", "hit_species", "hit_taxon_id", "hit_gene_symbol",
                        "hit_protein_name", "hit_description", "hit_metadata_state",
                        "hit_metadata_records_json",
                        "raw_reference_gene_id", "percent_identity", "query_coverage",
                        "evalue", "bitscore", "provider", "database", "identity_policy",
                    )
                } | {"interpretation": "sequence_similarity_handle_only",
                    "rank_contribution": "none"})
            # Stable input order is not guaranteed by arbitrary TSV producers.
            candidates.sort(key=lambda row: (
                str(row.get("hit_species", "")).casefold(),
                str(row.get("hit_accession", "")).casefold(),
                str(row.get("raw_reference_gene_id", "")).casefold(),
            ))
            identity_hits_by_gene[str(gene)] = candidates

    reference_materials, material_provenance = load_reference_materials(
        atlas_path, atlas.to_dict("records"))

    gene_set = set(atlas["gene_id"].astype(str))
    literature_records: list[dict[str, object]] = []
    for row in literature.to_dict("records"):
        linked_genes = [gene for gene in _tokens([row.get("interval_gene_ids", "")]) if gene in gene_set]
        if not linked_genes:
            continue
        literature_records.append({**row, "_linked_genes": linked_genes})

    work = atlas.copy()
    keys = work.apply(_unit_key, axis=1)
    work["_unit_type"] = [item[0] for item in keys]
    work["_unit_anchor"] = [item[1] for item in keys]
    atlas_origin = None

    units: list[dict[str, object]] = []
    unit_rows: list[dict[str, object]] = []
    for (unit_type, anchor), group in work.groupby(["_unit_type", "_unit_anchor"], sort=True):
        members = sorted(group["gene_id"].astype(str).tolist())
        unit_id = _unit_id(str(unit_type), str(anchor), members)
        member_set = set(members)
        linked_articles = [
            row for row in literature_records
            if member_set.intersection(set(row["_linked_genes"]))
        ]
        linked_articles.sort(key=lambda row: (
            str(row.get("reference_gene", "")),
            int(str(row.get("article_rank", "0") or "0")),
            str(row.get("source_id", "")),
            str(row.get("title", "")),
        ))
        lookup_records = [lookup_by_gene[gene] for gene in members if gene in lookup_by_gene]
        material_cards = [card for gene in members for card in reference_materials.get(gene, [])]
        sources = _stable_source_ids(linked_articles + lookup_records + material_cards)

        gene_records = []
        for record in group.sort_values("gene_id", kind="stable").to_dict("records"):
            gene_record = {
                key: str(record.get(key, ""))
                for key in (
                    "gene_id", "chrom", "start", "end", "protein_length",
                    "local_annotation_description", "local_annotation_go_terms",
                    "family_id", "family_label", "evolutionary_family_id",
                    "evolutionary_family_status", "orthogroup",
                    "orthology_relationship", "orthology_reference_gene",
                    # Family-index BLAST endpoints are retrieval handles, not
                    # orthology. Keep their exact IDs and ambiguity metrics in
                    # the frozen unit/prompt so semantic review can see the
                    # original assignment context without upgrading support.
                    "family_index_best_subject", "family_index_subjects",
                    "family_index_subject_endpoints", "family_index_assignment_status",
                    "family_index_best_identity", "family_index_best_query_coverage",
                    "family_index_best_evalue", "family_index_best_bitscore",
                    "family_index_second_bitscore", "family_index_top2_bitscore_ratio",
                    "family_index_symbol_source_rank", "family_index_symbol_source_bitscore_ratio",
                    # Keep endpoint namespace/species when the atlas supplies
                    # it; transfer verification uses these only when present.
                    "orthology_reference_species", "reference_species",
                    "orthology_reference_accession_namespace",
                    "reference_accession_namespace", "accession_namespace",
                    "subfamily_id", "domain_ids", "domain_architecture",
                    "broad_function_bins", "primary_broad_function",
                    "conservation_class", "evidence_completeness",
                    "annotation_rescue_needed", "gene_model_completeness",
                    "classification_basis", "notes",
                )
            }
            source_record = lookup_by_gene.get(str(record.get("gene_id", "")), {})
            # Keep an atlas endpoint's explicit species/namespace constraints
            # on that endpoint when an independent lookup replaces the generic
            # reference fields. Never move another atlas reference's metadata.
            atlas_endpoint = str(record.get("orthology_reference_gene", "")).strip()
            atlas_reference = str(record.get("reference_gene", "")).strip()
            if source_record and atlas_endpoint and atlas_reference in {"", atlas_endpoint}:
                for endpoint_field, reference_field in (
                    ("orthology_reference_species", "reference_species"),
                    ("orthology_reference_accession_namespace", "reference_accession_namespace"),
                ):
                    if not gene_record[endpoint_field]:
                        gene_record[endpoint_field] = str(record.get(reference_field, ""))
            # A lookup is an independent reference record. Filling its absent
            # columns from the atlas would attach another endpoint's function,
            # species or citation to it. Preserve acquired atlas material as a
            # separate investigation context rather than overwrite or join it.
            reference_record = source_record if source_record else record
            atlas_has_reference_material = any(
                _has_material(record.get(key, ""))
                for key in _FUNCTION_MATERIAL_FIELDS + _GENE_MATERIAL_FIELDS if key in _LOOKUP_CONTEXT_FIELDS)
            if source_record and atlas_has_reference_material and any(
                str(record.get(key, "")) != str(source_record.get(key, "")) for key in _LOOKUP_CONTEXT_FIELDS
            ):
                if atlas_origin is None:
                    atlas_origin = file_fingerprint(atlas_path, full_hash=True)
                gene_record["atlas_material_context"] = {
                    **{key: str(record[key]) for key in _ATLAS_CONTEXT_FIELDS if key in record},
                    "origin": atlas_origin,
                    "material_status": "investigation_context_not_verified_support",
                    "rank_contribution": "none",
                }
            gene_record.update({
                key: str(reference_record.get(key, ""))
                for key in _LOOKUP_CONTEXT_FIELDS
            })
            gene_record["sequence_similarity_hits"] = identity_hits_by_gene.get(
                str(record.get("gene_id", "")), [])
            gene_records.append(gene_record)

        article_payload = [
            {
                "reference_gene": str(row.get("reference_gene", "")),
                "reference_function": str(row.get("reference_function", "")),
                "source_id": str(row.get("source_id", "")),
                "source_url": str(row.get("source_url", "")),
                "doi": str(row.get("doi", "")),
                "pmid": str(row.get("pmid", "")),
                "pmcid": str(row.get("pmcid", "")),
                "title": str(row.get("title", "")),
                "journal": str(row.get("journal", "")),
                "year": str(row.get("year", "")),
                "matched_phenotype_terms": str(row.get("matched_phenotype_terms", "")),
                "query_status": str(row.get("query_status", "")),
                "claim_type": str(row.get("claim_type", "")),
                # Keep declared endpoint constraints with their own source;
                # older tables without these fields retain their old payload.
                **{key: str(row[key]) for key in ("reference_species", "reference_accession_namespace") if key in row},
            }
            for row in linked_articles
        ]

        # Keep source URLs and evidence metadata attached to the identity they
        # describe.  The roster above is only a closed-world identifier list;
        # this payload is the auditable explanation of where each identifier
        # came from and what kind of corpus record supplied it.
        source_provenance = []
        for row in lookup_records:
            source_provenance.append({
                "gene_id": str(row.get("gene_id", "")),
                "reference_gene": str(row.get("reference_gene", "")),
                "source_ids": _stable_source_ids([row]),
                "source_url": str(row.get("reference_source_url", "") or row.get("source_url", "")),
                "evidence_codes": str(row.get("evidence_codes", "")),
                "evidence_type": str(row.get("evidence_type", "")),
                "provenance_tier": str(row.get("provenance_tier", "")),
                "provenance_detail": str(row.get("provenance_detail", "")),
                "source_database": str(row.get("source_database", "")),
                "source_release": str(row.get("source_release", "")),
                "retrieved_at": str(row.get("retrieved_at", "")),
                "corpus_version": str(row.get("corpus_version", "")),
                "mapping_state": str(row.get("mapping_state", "")),
                "mapping_relationship": str(row.get("mapping_relationship", "")),
                "function_transfer_scope": str(row.get("function_transfer_scope", "")),
            })

        for card in material_cards:
            source_provenance.append({
                "gene_id": card["interval_gene_id"],
                "reference_gene": card["reference_gene"],
                "reference_species": card.get("reference_species", ""),
                "source_ids": _stable_source_ids([card]),
                "source_url": card.get("source_url", ""),
                "material_status": card["material_status"],
            })

        status_payload: list[dict[str, str]] = []
        if not literature_status.empty and "reference_gene" in literature_status.columns:
            references = set(_tokens(group.get("orthology_reference_gene", pd.Series(dtype=str)).tolist()))
            for status_row in literature_status.to_dict("records"):
                if str(status_row.get("reference_gene", "")) in references:
                    status_payload.append({key: str(status_row.get(key, "")) for key in literature_status.columns})

        unit = {
            "unit_id": unit_id,
            "unit_type": str(unit_type),
            "unit_anchor": str(anchor),
            "member_interval_genes": members,
            "identity": {
                "evolutionary_family_ids": _tokens(group.get("evolutionary_family_id", pd.Series(dtype=str)).tolist()),
                "family_ids": _tokens(group.get("family_id", pd.Series(dtype=str)).tolist()),
                "orthogroup_ids": _tokens(group.get("orthogroup", pd.Series(dtype=str)).tolist()),
                "subfamily_ids": _tokens(group.get("subfamily_id", pd.Series(dtype=str)).tolist()),
            },
            "genes": gene_records,
            "reference_material_cards": material_cards,
            "literature_link_candidates": article_payload,
            "literature_query_status": status_payload,
            "stable_source_ids": sources,
            "source_provenance": source_provenance,
            "source_backed_function_record_count": sum(
                1 for row in lookup_records
                if str(row.get("reference_source_id", "")).strip()
                and any(_has_material(row.get(key, "")) for key in (
                    "reference_function", "reference_pathway", "reference_mutant_phenotype"
                ))
            ),
            "rank_contribution": "none",
        }
        unit["additional_reference_material_count"] = len(material_cards)
        unit["reference_material_function_record_count"] = sum(
            bool(_stable_source_ids([card])) and any(_has_material(card.get(key, ""))
                for key in ("function_description", "pathway", "mutant_phenotype"))
            for card in material_cards)
        unit["investigation_material"] = _investigation_material(unit, include_similarity_annotations=True)
        unit["source_material_class"] = unit["investigation_material"]["source_material_class"]
        unit["request_state"] = _unit_request_state(unit)
        unit["review_status"] = "not_evaluated"
        unit["not_evaluated_reason"] = "awaiting_semantic_review" if unit["request_state"] == REQUEST_STATE_REQUESTED else "no_reviewable_material"
        estimated_prompt = _prompt_record(unit, phenotype, case_id, b_review_membership_policy)
        input_bytes = len(estimated_prompt["system"].encode("utf-8")) + len(estimated_prompt["user"].encode("utf-8"))
        unit["estimated_review_cost"] = {
            "status": "estimate_not_actual_usage",
            "method": "ceil(system_and_user_utf8_bytes/4); tokenizer-independent input heuristic; excludes estimate metadata",
            "input_utf8_bytes": input_bytes,
            "estimated_input_tokens": (input_bytes + 3) // 4,
            "estimated_model_call_count": int(unit["request_state"] == REQUEST_STATE_REQUESTED),
            "actual_input_tokens": None,
        }
        units.append(unit)
        unit_rows.append(_unit_summary(unit))

    covered_genes = {gene for unit in units for gene in unit["member_interval_genes"]}
    if covered_genes != gene_set:
        missing = sorted(gene_set - covered_genes)
        extra = sorted(covered_genes - gene_set)
        raise RuntimeError(f"Route B semantic bundle lost candidate identity: missing={missing}, extra={extra}")

    coverage_by_gene: dict[str, str] = {}
    for record in atlas.to_dict("records"):
        coverage_by_gene[str(record["gene_id"])] = _classify_gene_b0_coverage(record)
    b0_coverage_counts = {
        status: sum(1 for gene_status in coverage_by_gene.values() if gene_status == status)
        for status in B0_COVERAGE_STATUSES
    }
    b0_coverage_total = sum(b0_coverage_counts.values())
    if b0_coverage_total != len(gene_set):
        raise RuntimeError(
            f"B0 coverage sum ({b0_coverage_total}) does not match interval gene count ({len(gene_set)})"
        )
    identity_signal_counts = _identity_signal_counts(atlas)

    interval_gene_count = len(gene_set)
    semantic_unit_count = len(units)
    single_member_unit_count = sum(1 for unit in units if len(unit["member_interval_genes"]) == 1)
    multi_member_unit_count = semantic_unit_count - single_member_unit_count
    compression_observability = {
        "status": "measured" if interval_gene_count and semantic_unit_count else "not_applicable_empty_interval",
        "compression_gene_per_unit": (
            round(interval_gene_count / semantic_unit_count, 6) if semantic_unit_count else 0.0
        ),
        "unit_to_gene_ratio": (
            round(semantic_unit_count / interval_gene_count, 6) if interval_gene_count else 0.0
        ),
        "single_member_unit_count": single_member_unit_count,
        "single_member_unit_fraction": (
            round(single_member_unit_count / semantic_unit_count, 6) if semantic_unit_count else 0.0
        ),
        "multi_member_unit_count": multi_member_unit_count,
        "multi_member_unit_fraction": (
            round(multi_member_unit_count / semantic_unit_count, 6) if semantic_unit_count else 0.0
        ),
        "coverage_gene_counts": dict(b0_coverage_counts),
        "identity_signal_counts": dict(identity_signal_counts),
        "rank_contribution": "none",
    }

    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    unit_table_path = out / "route_b_semantic_units.tsv"
    pd.DataFrame(unit_rows, columns=UNIT_COLUMNS).to_csv(
        unit_table_path, sep="\t", index=False, lineterminator="\n"
    )

    bundle = {
        "schema_version": SCHEMA_VERSION,
        "initial_material_projection": INITIAL_MATERIAL_PROJECTION,
        "case_id": case_id,
        "route": "interval_first",
        "stage": "post_b0_semantic_bundle",
        "b_review_membership_policy": b_review_membership_policy,
        "focal_species": focal_species,
        "phenotype": str(phenotype).strip(),
        "interval_gene_count": len(gene_set),
        "semantic_unit_count": len(units),
        "candidate_identity_preserved": True,
        "b0_coverage": {
            "counts": b0_coverage_counts,
            "total": b0_coverage_total,
            "sum_matches_total": True,
            "by_gene": coverage_by_gene,
            "identity_signal_counts": dict(identity_signal_counts),
            "rule": "every interval gene is accounted for as assigned, unassigned, unresolved, or failed; sum must equal total interval genes",
            "reading_note": "an unresolved evolutionary_family_status or conservation_class classifies a gene as unresolved before its identity columns are read, so assigned=0 does not mean no family was assigned; identity_signal_counts shows which layer is actually missing",
        },
        "units": units,
        "policy": [
            "the B0 interval inventory must be frozen before phenotype is introduced here",
            "Route A artifacts and benchmark truth are not inputs",
            "every interval gene appears in exactly one semantic unit",
            "unclassified genes remain explicit gene-level units",
            "literature rows are link candidates, not verified evidence",
            "source-backed B0 function rows are frozen inputs, not model-authored evidence",
            "no unit or gene receives a numerical score",
        ],
    }
    bundle_path = out / "route_b_semantic_bundle.json"
    prompt_path = out / "route_b_semantic_prompts.jsonl"
    prompt_fingerprints: list[str] = []
    requested_units = [unit for unit in units if unit["request_state"] == REQUEST_STATE_REQUESTED]
    skipped_units = [unit for unit in units if unit["request_state"] != REQUEST_STATE_REQUESTED]
    requested_unit_count = len(requested_units)
    skipped_gene_count = sum(len(unit["member_interval_genes"]) for unit in skipped_units)
    compression_observability.update({
        # Phase 1 defined the expected request count as the unit count. Phase 2
        # skips units that carry no reviewable material, so the two numbers are
        # reported separately rather than one silently redefining the other.
        "semantic_unit_count": semantic_unit_count,
        "requested_unit_count": requested_unit_count,
        "skipped_unit_count": len(skipped_units),
        "skipped_gene_count": skipped_gene_count,
        "skipped_unit_reason": REQUEST_STATE_SKIPPED if skipped_units else None,
        "expected_llm_request_count": requested_unit_count,
        "expected_requests_per_100_genes": (
            round(requested_unit_count / interval_gene_count * 100.0, 6) if interval_gene_count else 0.0
        ),
        "request_skip_policy": (
            "schedule from investigation_material.semantic_review_recommended; function, domain and "
            "annotation context permit review, while identity IDs, bare sequence hits and source links "
            "remain retrieval handles. All units and genes remain in the bundle with not_evaluated "
            "review status until actual interpretation; a skipped request is never negative evidence"
        ),
    })
    review_cost_estimate = {
        "status": "estimate_not_actual_usage",
        "estimated_model_call_count": requested_unit_count,
        "estimated_input_tokens": sum(unit["estimated_review_cost"]["estimated_input_tokens"] for unit in requested_units),
        "actual_input_tokens": None,
        "method": "sum of requested per-unit input estimates; excludes output tokens and pricing",
    }
    bundle["review_cost_estimate"] = review_cost_estimate
    skipped_path = out / "route_b_semantic_units_not_requested.tsv"
    pd.DataFrame(
        [{"unit_id": unit["unit_id"], "unit_type": unit["unit_type"],
          "member_count": len(unit["member_interval_genes"]),
          "member_interval_genes": ";".join(unit["member_interval_genes"]),
          "request_state": unit["request_state"],
          "source_material_class": unit["source_material_class"],
          "not_evaluated_reason": unit["not_evaluated_reason"],
          "estimated_input_tokens": unit["estimated_review_cost"]["estimated_input_tokens"],
          "review_status": "not_evaluated", "rank_contribution": "none"}
         for unit in skipped_units],
        columns=["unit_id", "unit_type", "member_count", "member_interval_genes",
                 "request_state", "source_material_class", "not_evaluated_reason", "estimated_input_tokens", "review_status", "rank_contribution"],
    ).to_csv(skipped_path, sep="\t", index=False, lineterminator="\n")
    with prompt_path.open("w", encoding="utf-8") as handle:
        for unit in requested_units:
            # The JSON bundle retains the full hit rows for audit. Prompts use
            # a compact view to avoid repeating header JSON and BLAST provider
            # metadata for every member of a large semantic unit.
            record = _prompt_record(unit, phenotype, case_id, b_review_membership_policy)
            encoded = json.dumps(record, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
            prompt_fingerprints.append(hashlib.sha256(encoded.encode("utf-8")).hexdigest())
            handle.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")

    bundle["execution_prompt_fingerprints"] = prompt_fingerprints
    bundle["execution_prompt_count"] = len(prompt_fingerprints)
    bundle_path.write_text(
        json.dumps(bundle, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id,
        "route": "interval_first",
        "stage": "post_b0_semantic_bundle",
        "interval_gene_count": len(gene_set),
        "semantic_unit_count": len(units),
        "gene_level_unclassified_unit_count": int(sum(1 for unit in units if unit["unit_type"] == "interval_gene")),
        "literature_link_count": int(sum(len(unit["literature_link_candidates"]) for unit in units)),
        "stable_source_count": len({source for unit in units for source in unit["stable_source_ids"]}),
        "candidate_identity_preserved": True,
        "b0_coverage": {
            "counts": b0_coverage_counts,
            "total": b0_coverage_total,
            "sum_matches_total": True,
            "identity_signal_counts": dict(identity_signal_counts),
            "rule": "every interval gene is accounted for as assigned, unassigned, unresolved, or failed; sum must equal total interval genes",
            "reading_note": "an unresolved evolutionary_family_status or conservation_class classifies a gene as unresolved before its identity columns are read, so assigned=0 does not mean no family was assigned; identity_signal_counts shows which layer is actually missing",
        },
        "compression_observability": compression_observability,
        "prompt_count": len(requested_units),
        "prompt_fingerprints": prompt_fingerprints,
        "review_cost_estimate": review_cost_estimate,
        "rank_contribution": "none",
        "inputs": {
            "atlas": file_fingerprint(atlas_path, full_hash),
            "interval_reverse_lookup": file_fingerprint(lookup_path, full_hash) if lookup_path else None,
            "route_b_identity_hits": file_fingerprint(identity_hits_path, full_hash) if identity_hits_path else None,
            "reference_materials": material_provenance,
            "literature": file_fingerprint(literature_path, full_hash) if literature_path else None,
            "literature_status": file_fingerprint(literature_status_path, full_hash) if literature_status_path else None,
        },
        "outputs": {
            "bundle": str(bundle_path),
            "units": str(unit_table_path),
            "prompts": str(prompt_path),
            "units_not_requested": str(skipped_path),
            "manifest": str(out / "route_b_semantic_bundle_manifest.json"),
        },
        "policy": bundle["policy"],
    }
    manifest_path = out / "route_b_semantic_bundle_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        "bundle": str(bundle_path),
        "units": str(unit_table_path),
        "prompts": str(prompt_path),
        "manifest": str(manifest_path),
        "interval_gene_count": len(gene_set),
        "semantic_unit_count": len(units),
        "prompt_count": len(requested_units),
        "units_not_requested": str(skipped_path),
        "skipped_unit_count": len(skipped_units),
        "b0_coverage_counts": b0_coverage_counts,
    }
