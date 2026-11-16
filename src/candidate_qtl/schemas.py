"""Canonical table schemas, shared enums, and validators.

This module is the single source of truth for the column contracts, enum
vocabularies, and cross-column checks that the rest of the package relies on.
Previously each producer (``discover``, ``evidence-build``, ``case-run``, the
external importers) validated its own frames ad hoc, so the same column name
could drift in nullability or meaning between commands. Centralising the
definitions here keeps those contracts consistent and gives one place to bump
:data:`SCHEMA_VERSION` when a contract changes.

The validators are dependency-free by default. ``pandera`` is an *optional*
extra (``pip install candidate-qtl[schema]``): when it is installed,
:func:`pandera_schema` derives a ``pandera.DataFrameSchema`` from the very same
:class:`TableSchema` descriptor, so the pandera view never diverges from the
built-in checks. The local-first runtime never requires pandera to be present.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import NamedTuple, Optional

import pandas as pd

# Bump when any canonical column contract below changes in a
# backwards-incompatible way. Recorded into run manifests so an artifact can be
# traced back to the contract it was produced under.
SCHEMA_VERSION = "1.8"

# ---------------------------------------------------------------------------
# Shared enum vocabularies
# ---------------------------------------------------------------------------
# Evidence direction. ``context`` rows are never scored; they carry provenance
# and caveats only.
DIRECTIONS = frozenset({"support", "against", "rescue", "context"})
# Qualitative confidence attached to each evidence row.
CONFIDENCE_LEVELS = frozenset({"high", "medium", "low"})
# Evidence provenance/dependence vocabulary. ``unknown`` is the
# bucket for rows whose provenance kind is not stated; it is never counted as
# an independent evidence cluster.
EVIDENCE_DEPENDENCY_KINDS = frozenset({
    "primary_observation", "curated_assertion", "orthology_transfer",
    "sequence_domain_call", "prediction", "unknown",
})
# Dry-prioritisation bucket assigned during ranking.
DRY_GROUPS = frozenset({"top", "reserve", "rescue", "not_selected_yet"})
ASSAY_STATES = frozenset({"not_planned", "planned", "primer_designed", "ordered", "wet_pending", "failed", "partial_result", "completed", "needs_redesign"})
CANDIDATE_WET_STATES = frozenset({"not_tested", "partially_covered", "coding_region_cleared", "no_parental_difference", "candidate_variant_found_unconfirmed", "candidate_variant_confirmed", "fails_cosegregation", "excluded_by_recombinant", "supports_qtl", "requires_independent_validation", "inconclusive"})
HYPOTHESIS_TYPES = frozenset({
    "single_gene_coding", "single_gene_splice", "promoter_or_utr", "te_regulatory",
    "sv_cnv_pav", "annotation_missing", "tandem_or_paralog_redundancy",
    "multi_gene_combination", "linked_qtl", "intergenic_regulatory", "unknown_locus_mechanism",
})
HYPOTHESIS_STATUSES = frozenset({
    "proposed", "pending_review", "ready_to_test", "testing", "partially_tested", "supported",
    "weakened", "rejected", "inconclusive", "superseded", "requires_independent_validation",
})
PROPOSAL_ORIGINS = frozenset({"manual", "llm", "workflow"})
ANNOTATION_RESCUE_CLASSES = frozenset({"native", "native_discordant", "reference_rescue", "ambiguous_overlap", "unmapped", "not_evaluated"})
ANNOTATION_RESCUE_REVIEW_STATES = frozenset({"not_required", "pending"})
# Annotation-rescue model classes that deterministically seed an annotation_missing
# locus hypothesis: the native annotation disagrees with, or is
# entirely absent for, a rescue projection. ``native`` and ``unmapped``/``not_evaluated``
# are deliberately excluded: ``native`` means the projection agrees with the existing
# annotation, and ``unmapped``/``not_evaluated`` mean the query was never resolved in
# this interval, neither of which is evidence of an annotation gap.
ANNOTATION_RESCUE_HYPOTHESIS_ESCALATION_CLASSES = frozenset({"native_discordant", "reference_rescue", "ambiguous_overlap"})
# Comparative-genomics resource roles. Deliberately generic: no crop name or
# gene-ID system is hard-coded, and this panel is distinct from the LLM-proposed
# known-gene panel (reference_gene_panel).
REFERENCE_RESOURCE_ROLES = frozenset({"native", "close_anchor", "functional_reference", "outgroup"})
REFERENCE_RESOURCE_TYPES = frozenset({"genome_fasta", "annotation_gff3", "proteome_fasta", "sequence_db", "annotation_db", "rest_api", "other"})
REFERENCE_RESOURCE_AVAILABILITY = frozenset({"available", "skipped", "declared_external"})
# External comparative-ranking ledger. External ranks are kept in a separate
# ledger and never summed into the native ``dry_score``; the consensus is a
# reciprocal-rank-fusion screening view, not a causal probability.
EXTERNAL_RANK_COVERAGE = frozenset({"covered", "uncovered", "na"})
EXTERNAL_RANK_CLASSES = frozenset({"consensus_top", "native_only", "external_rescue", "discordant", "uncovered"})
# Cross-crop / cross-resource benchmark conditions used to check that the
# meta-ranking is not overfit to one crop or one resource configuration.
CROSS_CROP_RESOURCE_CONDITIONS = frozenset({"full_native", "anchor_functional", "local_only", "external_uncovered"})
# Query eligibility gate: the known-gene benchmark on Cameor V2 must use only
# non-Pisum functional homologs (ordinary map-based cloning). Any Pisum-derived
# query (another Pisum assembly, a projection, a renamed record, or a sequence
# backmapped from the answer) would be searching the answer itself and is barred.
QUERY_SELECTION_BASIS = frozenset({
    "phenotype_literature_search", "mechanism_literature_search", "mutant_phenotype_similarity",
    "pathway_member", "domain_hypothesis", "reference_ortholog_propagation",
    "answer_backmapped", "pisum_sequence_based",
})
QUERY_SEQUENCE_ORIGIN = frozenset({
    "native_non_pisum_database", "non_pisum_reference_proteome",
    "pisum_reference_proteome", "pisum_derived_projection", "unknown",
})
BENCHMARK_ELIGIBILITY = frozenset({"eligible", "excluded_pisum", "excluded_answer_derived", "excluded_unknown_provenance"})
# Deterministic trait-search-term propagation. Function vocabulary is derived
# mechanically from reference-genome annotations of interval-candidate
# orthologs; every term keeps provenance and stays pending human review.
TERM_MATCH_TYPES = frozenset({"ortholog_propagated"})
TERM_REVIEW_STATES = frozenset({"pending", "approved", "rejected"})
TERM_TRAIT_RELEVANCE = frozenset({"trait_matched", "unrelated", "not_evaluated"})
TERM_REFERENCE_ROLES = frozenset({"native", "close_anchor", "functional_reference", "outgroup", "unspecified"})
# Rank-lane semantic status. ``evaluable`` means the rank was
# scored against supporting evidence; the other three values are distinct
# not_evaluated buckets recorded in the case-run / guided-profile manifests:
# ``not_evaluated_no_scoring_evidence`` covers a lane with zero evidence rows,
# ``not_evaluated_evidence_provider_failed`` a provider that could not be read,
# and ``not_evaluated_no_supporting_evidence`` a lane whose evidence
# rows exist but carry no ``support`` direction row -- the run is penalty-only
# and is excluded from the rank-metric denominator exactly like the
# no-scoring-evidence bucket.
RANK_SEMANTIC_STATUSES = frozenset({
    "evaluable",
    "not_evaluated_no_scoring_evidence",
    "not_evaluated_evidence_provider_failed",
    "not_evaluated_no_supporting_evidence",
})
# Memorization-probe stratification. The probe answer is matched
# against the catalog family through the shared family_labels identity-key
# registry (curated aliases with provenance included). "match" / "mismatch"
# are decided only via a shared known key or an explicit probe abstention;
# every other non-matching relation is "not_evaluated_probe_ambiguous" so a
# probe answered with an un-catalogued synonym is never misreported as
# recall_negative (memory must not be misclassified as reasoning).
PROBE_STRATIFICATION_STATUSES = frozenset({
    "match",
    "mismatch",
    "not_evaluated_probe_ambiguous",
})
# QTL candidate workflow 4 use-case strata.
# Categorizes benchmark/real-world QTL candidate triage situations into 4 distinct types
# plus 'unclassified'. Type 1–3 form the primary evaluable universe for the conditional
# Top-20 retention contract. Type 4 represents novel or species-specific cases with no
# functional or homology prior, which are out of Top-20 guarantee scope and must never
# be mixed into the Type 1–3 recall denominator.
USE_CASE_TYPES = frozenset({
    "type1_reference_or_id_mismatch",
    "type2_semantic_or_mechanism_mismatch",
    "type3_cross_species_transferable",
    "type4_species_specific_novel",
    "unclassified",
})
RETENTION_SCOPE_USE_CASE_TYPES = frozenset({
    "type1_reference_or_id_mismatch",
    "type2_semantic_or_mechanism_mismatch",
    "type3_cross_species_transferable",
})
USE_CASE_TYPE_DESCRIPTIONS = {
    "type1_reference_or_id_mismatch": "Causal gene reported; missed due to retrieval gaps, legacy IDs, old references, or assembly mismatch. Evaluates ID integration, coordinate provenance, liftover, and homology recovery.",
    "type2_semantic_or_mechanism_mismatch": "Causal gene reported; missed due to divergent phenotype or mechanism vocabularies. Evaluates phenotype reframing, multi-mechanism hypotheses, and semantic retrieval.",
    "type3_cross_species_transferable": "Unknown in focal species; homologous genes or functional family known in other plants. Evaluates multi-bait, family/orthogroup/domain, and cross-species function transfer.",
    "type4_species_specific_novel": "Novel or species-specific locus with no homologous gene or functional prior. Out of Top-20 guarantee scope; evaluates honest abstention and escalation to mapping/Wet/rescue branches.",
    "unclassified": "Unclassified use case; retained explicitly without guessing.",
}

# Non-gene causal architecture (issue #444 "Non-gene causal architecture"): an
# axis ORTHOGONAL to use_case_type. It is not a fifth use-case type, and it is
# not a conversion of a non-gene cause into a gene truth. Deliberately reuses
# the non-gene members of docs/LOCUS-CAUSAL-SCENARIOS.ja.md's ``scenario_type``
# vocabulary instead of introducing a parallel one, so a locus is described
# once. A case with gene-level truth records "not_applicable" here; it is never
# forced to pick a non-gene value, and recording a non-gene architecture never
# changes which use_case_type the case belongs to.
NON_GENE_CAUSAL_ARCHITECTURES = frozenset({
    "sv_te_cnv_pav",
    "regulatory_element",
    "expression",
    "multi_gene_or_epistasis",
    "not_applicable",
})


def require_rank_semantic_status(status: str) -> str:
    """Fail-closed gate for the rank-lane semantic-status vocabulary (issue #397).

    Mirrors how the table enums are checked (``DIRECTIONS`` etc. through
    :func:`collect_errors`): a status outside :data:`RANK_SEMANTIC_STATUSES`
    is a programming or vocabulary error and raises instead of silently
    propagating a novel literal into a manifest.
    """
    if status not in RANK_SEMANTIC_STATUSES:
        raise ValueError(
            f"rank_semantic_status {status!r} is not in the registered vocabulary "
            f"{sorted(RANK_SEMANTIC_STATUSES)}"
        )
    return status


def require_probe_stratification(status: str) -> str:
    """Fail-closed gate for the memorization-probe stratification vocabulary.

    Mirrors :func:`require_rank_semantic_status`: a status outside
    :data:`PROBE_STRATIFICATION_STATUSES` is a programming or vocabulary error
    and raises instead of silently propagating a novel literal into a manifest.
    """
    if status not in PROBE_STRATIFICATION_STATUSES:
        raise ValueError(
            f"probe_match_status {status!r} is not in the registered vocabulary "
            f"{sorted(PROBE_STRATIFICATION_STATUSES)}"
        )
    return status


def require_use_case_type(value: str) -> str:
    """Fail-closed gate for the 4 QTL use-case types vocabulary (issue #444).

    Mirrors :func:`require_rank_semantic_status`: a value outside
    :data:`USE_CASE_TYPES` is a programming or schema error and raises
    instead of silently propagating a novel literal into a manifest.
    """
    if value not in USE_CASE_TYPES:
        raise ValueError(
            f"use_case_type {value!r} is not in the registered vocabulary "
            f"{sorted(USE_CASE_TYPES)}"
        )
    return value


def normalize_non_gene_causal_architecture(value: object) -> str:
    """Normalize the non-gene causal-architecture axis (issue #444).

    Empty/missing means "no non-gene architecture recorded for this locus",
    which is ``not_applicable`` -- a gene-truth case is never forced to pick a
    non-gene value. An explicit value outside
    :data:`NON_GENE_CAUSAL_ARCHITECTURES` raises, so a novel literal cannot
    reach a manifest and later be mistaken for a fifth use-case type.
    """
    text = str(value or "").strip().lower()
    if not text:
        return "not_applicable"
    if text not in NON_GENE_CAUSAL_ARCHITECTURES:
        raise ValueError(
            f"non_gene_causal_architecture {value!r} is not in the registered vocabulary "
            f"{sorted(NON_GENE_CAUSAL_ARCHITECTURES)}"
        )
    return text

# Provenance DAG ledger. Nodes record evidence sources and route
# outputs; edges record how each node depends on another. A source shared by
# several routes stays one node, so dependence is explicit rather than
# implicit in row counts.
PROVENANCE_NODE_KINDS = frozenset({"primary_experiment", "database_annotation", "review", "route_output", "literature"})
PROVENANCE_EDGE_KINDS = frozenset({"derived_from", "annotates", "reviewed_in", "produced"})
PROVENANCE_NODE_COLUMNS = ["node_id", "node_kind", "source_id", "label", "source_ref"]
PROVENANCE_EDGE_COLUMNS = ["child_node", "parent_node", "edge_kind"]
PROVENANCE_GRAPH_COLUMNS = ["node_id", "node_kind", "source_id", "label", "source_ref", "parent_node", "edge_kind"]

WET_ASSAY_COLUMNS = [
    "assay_id", "gene_id", "transcript_id", "hypothesis_id", "assay_type", "planned_region",
    "region_class", "assay_state", "failure_reason", "source_ref", "review_status", "reviewer",
]
WET_COVERAGE_COLUMNS = [
    "assay_id", "gene_id", "transcript_id", "planned_region", "actual_covered_region",
    "coverage_fraction", "region_class", "required_region", "quality_status", "read_direction",
    "off_target_status", "multi_product_status", "source_ref",
]
WET_RESULT_COLUMNS = [
    "result_id", "assay_id", "gene_id", "sample_id", "sample_role", "result_type", "result_value",
    "variant_id", "variant_confirmation_status", "validation_result", "quality_status", "source_ref",
    "review_status", "reviewer", "requested_candidate_action", "generation", "individual_id",
    "amplification_status", "genotype_call", "phenotype_value", "cosegregation_result",
    "recombinant_result", "sequence_quality", "read_direction",
]
WET_HISTORY_COLUMNS = [
    "history_id", "event_index", "entity_type", "entity_id", "assay_id", "gene_id", "from_state",
    "to_state", "transition_status", "transition_reason", "source_ref", "review_status", "reviewer",
]
LOCUS_HYPOTHESIS_COLUMNS = [
    "hypothesis_id", "case_id", "hypothesis_type", "interval_or_element_id", "member_gene_ids",
    "target_gene_ids", "proposed_mechanism", "predicted_molecular_effect",
    "predicted_phenotype_direction", "supporting_evidence_ids", "counterevidence_ids", "context_evidence_ids",
    "required_experiment_ids", "status", "requested_status", "status_guardrail", "status_reason", "created_by", "proposal_origin",
    "review_status", "source_id", "parent_hypothesis_id", "supersedes_hypothesis_id",
    "gene_disposition", "rejection_scope",
]
ANNOTATION_RESCUE_COLUMNS = [
    "query_protein_id", "query_sequence_sha256", "query_source_id", "query_species", "query_assembly_id", "query_annotation_release", "query_role", "query_selection_method", "query_review_status",
    "projection_id", "projection_rank", "hit_count", "near_tie", "model_class", "quality_tier", "provider_status", "chrom", "start", "end", "strand",
    "native_gene_ids", "native_annotation_quality", "query_coverage", "identity", "positive_fraction",
    "projected_cds_count", "projected_cds_length", "frameshift_count", "inframe_stop_count", "orf_audit",
    "splice_audit", "noncanonical_splice_count", "terminal_stop_feature", "interval_boundary_truncation",
    "native_overlap_fraction", "projection_overlap_fraction", "review_status", "promotion_status",
    "source_type", "source_id", "notes",
]
REFERENCE_RESOURCE_PANEL_COLUMNS = [
    "resource_id", "role", "resource_type", "taxon", "assembly_id", "annotation_version",
    "proteome_version", "db_version", "source", "uri_or_path", "external_api", "review_status",
    "availability", "skip_reason", "resource_sha256", "resource_size_bytes", "notes",
]
EXTERNAL_RANK_INPUT_COLUMNS = [
    "ranker_id", "ranker_family", "gene_id", "rank", "raw_score", "percentile", "coverage",
    "mapping_confidence", "reference_species", "source", "version", "dependency_group",
]
EXTERNAL_RANK_LEDGER_COLUMNS = [
    "gene_id", "native_screening_rank", "native_top", "covered_ranker_count", "covering_families",
    "consensus_rrf", "consensus_rank", "consensus_top", "classification", "external_rank_min",
    "external_rank_max", "external_rank_spread", "notes",
]
PROPAGATED_TERM_COLUMNS = [
    "candidate_gene_id", "term", "term_normalized", "ontology_id", "ontology_name",
    "reference_gene_id", "source_genome", "reference_role", "reference_species", "relationship",
    "term_source", "match_type", "term_confidence", "trait_relevance", "review_status", "notes",
]
QUERY_ELIGIBILITY_COLUMNS = [
    "case_id", "homolog_query_id", "homolog_source_id", "query_species", "query_genus", "query_taxon_id",
    "source_assembly", "selection_basis", "sequence_origin", "audited_species", "audited_genus",
    "query_sha256", "provenance_method", "benchmark_eligibility", "exclusion_reason",
]

# Confidence -> numeric multiplier used by the scoring layer. Kept beside the
# enum it belongs to so the two cannot drift apart.
CONFIDENCE_MULTIPLIERS = {"high": 1.0, "medium": 0.75, "low": 0.5}

# ---------------------------------------------------------------------------
# Canonical column orders
# ---------------------------------------------------------------------------
EVIDENCE_COLUMNS = [
    "gene_id",
    "evidence_code",
    "direction",
    "confidence",
    "value",
    "source_type",
    "source_id",
    "description",
]

# Optional provenance columns: which root source an evidence row
# derives from and what kind of assertion it carries. Absent columns and blank
# cells are valid; rows without a root are reported as dependency_unknown and
# never counted as independent evidence.
EVIDENCE_DEPENDENCY_COLUMNS = ["evidence_kind", "root_source_id"]


# ---------------------------------------------------------------------------
# Descriptor types
# ---------------------------------------------------------------------------
class CoordinateSpec(NamedTuple):
    """Genomic-coordinate columns that must satisfy ``start <= end`` etc."""

    seqid: str
    start: str
    end: str
    one_based: bool = True


@dataclass(frozen=True)
class TableSchema:
    """A declarative contract for one canonical table.

    Attributes
    ----------
    name:
        Human-readable table name, surfaced in error messages.
    columns:
        Canonical column order. Producers should emit exactly these columns.
    required:
        Columns that must be present for the frame to be usable.
    enums:
        Mapping of column name -> allowed value set. Empty strings are ignored
        here (use ``non_empty`` to forbid them) so a blank optional cell does
        not trip an enum check.
    non_empty:
        Columns whose every cell must be a non-blank string.
    unique_key:
        Column whose values must be unique across rows, if any.
    coordinates:
        Optional genomic-coordinate contract.
    """

    name: str
    columns: list[str]
    required: frozenset[str] = frozenset()
    enums: dict[str, frozenset[str]] = field(default_factory=dict)
    non_empty: frozenset[str] = frozenset()
    unique_key: Optional[str] = None
    coordinates: Optional[CoordinateSpec] = None


class SchemaError(NamedTuple):
    """A single validation failure. ``row`` is ``None`` for table-level issues."""

    table: str
    code: str
    message: str
    column: Optional[str] = None
    row: Optional[int] = None


class SchemaValidationError(ValueError):
    """Raised when a frame fails its :class:`TableSchema`.

    Carries the full list of :class:`SchemaError` entries on ``.errors`` so a
    caller (or a test) can inspect each row/column failure, mirroring pandera's
    aggregated ``lazy`` reporting.
    """

    def __init__(self, schema_name: str, errors: list[SchemaError]) -> None:
        self.schema_name = schema_name
        self.errors = errors
        preview = "; ".join(_format_error(err) for err in errors[:10])
        if len(errors) > 10:
            preview += f"; ... (+{len(errors) - 10} more)"
        super().__init__(f"{schema_name}: {len(errors)} schema error(s): {preview}")


def _format_error(err: SchemaError) -> str:
    location = err.column or "<table>"
    if err.row is not None:
        location = f"{location}[row {err.row}]"
    return f"{location}: {err.message}"


# ---------------------------------------------------------------------------
# Canonical table descriptors
# ---------------------------------------------------------------------------
EVIDENCE_SCHEMA = TableSchema(
    name="evidence",
    columns=[*EVIDENCE_COLUMNS, *EVIDENCE_DEPENDENCY_COLUMNS],
    required=frozenset(EVIDENCE_COLUMNS),
    enums={"direction": DIRECTIONS, "confidence": CONFIDENCE_LEVELS, "evidence_kind": EVIDENCE_DEPENDENCY_KINDS},
    non_empty=frozenset({"gene_id", "evidence_code", "direction", "confidence"}),
)

# Columns produced by ``extract`` (reference_extract.write_candidates). Only the
# structural core is enforced; descriptive columns stay optional.
CANDIDATE_COLUMNS = [
    "gene_id",
    "gene_name",
    "gene_biotype",
    "is_placeholder_gene",
    "description",
    "dbxref",
    "ontology_terms",
    "interval",
    "coordinate_basis",
    "gff_seqid",
    "fasta_seqid",
    "chrom",
    "start",
    "end",
    "strand",
    "gene_length",
    "transcript_count",
    "representative_transcript",
    "representative_reason",
    "cds_length",
    "protein_length",
    "cds_phase_status",
    "translation_status",
    "annotation_quality",
    "overlaps_qtl_boundary",
]

CANDIDATE_SCHEMA = TableSchema(
    name="candidates",
    columns=CANDIDATE_COLUMNS,
    required=frozenset({"gene_id", "chrom", "start", "end", "strand"}),
    non_empty=frozenset({"gene_id"}),
    unique_key="gene_id",
    coordinates=CoordinateSpec(seqid="chrom", start="start", end="end"),
)

WET_ASSAY_SCHEMA = TableSchema(
    name="wet_assay_registry",
    columns=WET_ASSAY_COLUMNS,
    required=frozenset(WET_ASSAY_COLUMNS),
    enums={"assay_state": ASSAY_STATES},
    non_empty=frozenset({"assay_id", "gene_id", "assay_type", "assay_state"}),
    unique_key="assay_id",
)

PROVENANCE_NODE_SCHEMA = TableSchema(
    name="provenance_nodes",
    columns=PROVENANCE_NODE_COLUMNS,
    required=frozenset({"node_id", "node_kind", "source_id"}),
    enums={"node_kind": PROVENANCE_NODE_KINDS},
    non_empty=frozenset({"node_id", "node_kind", "source_id"}),
    unique_key="node_id",
)

PROVENANCE_EDGE_SCHEMA = TableSchema(
    name="provenance_edges",
    columns=PROVENANCE_EDGE_COLUMNS,
    required=frozenset(PROVENANCE_EDGE_COLUMNS),
    enums={"edge_kind": PROVENANCE_EDGE_KINDS},
    non_empty=frozenset(PROVENANCE_EDGE_COLUMNS),
)

# Join view of the ledger: parent_node/edge_kind are legitimately blank for
# root nodes, so they are enum-checked but not required to be non-empty.
PROVENANCE_GRAPH_SCHEMA = TableSchema(
    name="provenance_graph",
    columns=PROVENANCE_GRAPH_COLUMNS,
    required=frozenset(PROVENANCE_GRAPH_COLUMNS),
    enums={"node_kind": PROVENANCE_NODE_KINDS, "edge_kind": PROVENANCE_EDGE_KINDS},
    non_empty=frozenset({"node_id", "node_kind", "source_id"}),
)

WET_HISTORY_SCHEMA = TableSchema(
    name="wet_state_history",
    columns=WET_HISTORY_COLUMNS,
    required=frozenset(WET_HISTORY_COLUMNS),
    non_empty=frozenset({"history_id", "event_index", "entity_type", "entity_id", "from_state", "to_state", "transition_status"}),
    unique_key="history_id",
)

LOCUS_HYPOTHESIS_SCHEMA = TableSchema(
    name="locus_hypotheses",
    columns=LOCUS_HYPOTHESIS_COLUMNS,
    required=frozenset(LOCUS_HYPOTHESIS_COLUMNS),
    enums={
        "hypothesis_type": HYPOTHESIS_TYPES,
        "status": HYPOTHESIS_STATUSES,
        "proposal_origin": PROPOSAL_ORIGINS,
    },
    non_empty=frozenset({"hypothesis_id", "case_id", "hypothesis_type", "status", "created_by", "proposal_origin"}),
    unique_key="hypothesis_id",
)

ANNOTATION_RESCUE_SCHEMA = TableSchema(
    name="annotation_rescue",
    columns=ANNOTATION_RESCUE_COLUMNS,
    required=frozenset(ANNOTATION_RESCUE_COLUMNS),
    enums={"model_class": ANNOTATION_RESCUE_CLASSES, "review_status": ANNOTATION_RESCUE_REVIEW_STATES},
    non_empty=frozenset({"query_protein_id", "query_sequence_sha256", "query_source_id", "query_role", "query_selection_method", "query_review_status", "model_class", "quality_tier", "provider_status", "chrom", "review_status", "promotion_status", "source_type"}),
)

QUERY_ELIGIBILITY_SCHEMA = TableSchema(
    name="query_eligibility",
    columns=QUERY_ELIGIBILITY_COLUMNS,
    required=frozenset(QUERY_ELIGIBILITY_COLUMNS),
    enums={"selection_basis": QUERY_SELECTION_BASIS, "sequence_origin": QUERY_SEQUENCE_ORIGIN, "benchmark_eligibility": BENCHMARK_ELIGIBILITY},
    non_empty=frozenset({"case_id", "homolog_query_id", "benchmark_eligibility"}),
)

REFERENCE_RESOURCE_PANEL_SCHEMA = TableSchema(
    name="reference_resource_panel",
    columns=REFERENCE_RESOURCE_PANEL_COLUMNS,
    required=frozenset(REFERENCE_RESOURCE_PANEL_COLUMNS),
    enums={"role": REFERENCE_RESOURCE_ROLES, "resource_type": REFERENCE_RESOURCE_TYPES, "availability": REFERENCE_RESOURCE_AVAILABILITY},
    non_empty=frozenset({"resource_id", "role", "resource_type", "taxon", "source", "uri_or_path", "availability"}),
    unique_key="resource_id",
)

EXTERNAL_RANK_INPUT_SCHEMA = TableSchema(
    name="external_rank_input",
    columns=EXTERNAL_RANK_INPUT_COLUMNS,
    required=frozenset({"ranker_id", "ranker_family", "gene_id", "coverage"}),
    enums={"coverage": EXTERNAL_RANK_COVERAGE},
    non_empty=frozenset({"ranker_id", "ranker_family", "gene_id", "coverage"}),
)

EXTERNAL_RANK_LEDGER_SCHEMA = TableSchema(
    name="external_rank_ledger",
    columns=EXTERNAL_RANK_LEDGER_COLUMNS,
    required=frozenset(EXTERNAL_RANK_LEDGER_COLUMNS),
    enums={"classification": EXTERNAL_RANK_CLASSES},
    non_empty=frozenset({"gene_id", "classification"}),
    unique_key="gene_id",
)

PROPAGATED_TERM_SCHEMA = TableSchema(
    name="propagated_terms",
    columns=PROPAGATED_TERM_COLUMNS,
    required=frozenset(PROPAGATED_TERM_COLUMNS),
    enums={"match_type": TERM_MATCH_TYPES, "review_status": TERM_REVIEW_STATES,
           "trait_relevance": TERM_TRAIT_RELEVANCE, "reference_role": TERM_REFERENCE_ROLES},
    non_empty=frozenset({"candidate_gene_id", "term", "term_normalized", "reference_gene_id", "match_type", "review_status", "reference_role", "trait_relevance"}),
)

# ---------------------------------------------------------------------------
# Family candidate bridge (route A/B/C -> one reporting table). The bridge is
# a parallel reporting artifact only: it never feeds ranking, and
# family_to_copy / narrowing wiring is untouched. All 22 columns are emitted
# by family_candidate_bridge; mapping_tier is a tuple (strongest first) so
# enum validation via ``Series.isin`` and pandera's ``sorted()`` both accept it.
# ---------------------------------------------------------------------------
BRIDGE_COLUMNS = [
    "case_id", "gene_id", "family_id", "route", "reference_accession", "source_id",
    "match_type", "identity", "query_coverage", "target_coverage", "evalue",
    "paralog_count", "ambiguity_status", "provider_status", "provider_version",
    "input_fingerprint", "output_fingerprint", "review_status",
    "dependency_cluster_id", "mapping_tier", "rank_contribution", "reason",
]
BRIDGE_MATCH_TYPES = frozenset({"curated", "rbh", "orthogroup", "domain", "blastp", "name_only"})
# Strongest -> weakest; ``tier6_not_evaluated`` rows may carry an empty
# match_type, which is why match_type is not in ``non_empty``.
BRIDGE_MAPPING_TIERS = (
    "tier1_source_backed_unique", "tier2_source_backed_ambiguous",
    "tier3_sequence_match_unique", "tier4_sequence_match_ambiguous",
    "tier5_name_only", "tier6_not_evaluated",
)
BRIDGE_AMBIGUITY_STATUSES = frozenset({
    "resolved_single_copy", "paralog_ambiguous", "tandem_ambiguous",
    "coverage_insufficient", "unclassified_gene", "not_evaluated",
})
BRIDGE_PROVIDER_STATUSES = frozenset({"complete", "failed", "unavailable",
                                      "imported_partial", "imported_complete", "not_run"})
BRIDGE_REVIEW_STATUSES = frozenset({"not_required", "review_required", "review_complete", "not_evaluated"})
BRIDGE_NOT_EVALUATED_REASONS = frozenset({
    "not_evaluated_provider_missing", "not_evaluated_provider_failed",
    "not_evaluated_provider_unavailable", "not_evaluated_no_verified_bait",
    "not_evaluated_no_reference_protein", "not_evaluated_coverage_below_threshold",
    "not_evaluated_interval_unresolved", "not_evaluated_unclassified_gene",
    "not_evaluated_not_in_index", "not_evaluated_family_unresolved_in_index",
})

# Family-to-copy route A evidence-source provenance (issue #374 route A layer).
# The ``copy_set_evidence_source`` cell of ``family_copy_discrimination.tsv`` is
# a sorted ';'-joined subset of this vocabulary, or "" when no layer supplied
# the row. family_to_copy enforces the vocabulary fail-closed (each ';'-split
# token must be a member or "") only when the route A blast input is active;
# the column itself exists only then, mirroring the bridge columns.
COPY_SET_EVIDENCE_SOURCES = frozenset({
    "atlas_orthogroup", "atlas_curated_family", "atlas_domain", "route_a_blastp",
})
# Route A blast paralog guard: ``ambiguity_status=unique`` only
# means "the single hit", not "the right hit". When the hit member gene has an
# in-interval relative recorded under the same bridge family (same route, same
# family_id -- a ``tier6_not_evaluated`` bridge row is a provider gap and is
# never an assignment), the true copy may be missing from the annotation and
# only a paralog remains, so the ``route_a_blast_guard`` cell records the hold
# reason ("" when the row injected). family_to_copy enforces the vocabulary
# fail-closed while the route A blast input is active, mirroring the
# ``copy_set_evidence_source`` gate.
ROUTE_A_BLAST_GUARD_REASONS = frozenset({"interval_paralog_present", "paralog_guard_not_evaluated"})
FAMILY_COPY_ROUTE_A_BLAST_COLUMNS = ["copy_set_evidence_source", "route_a_blast_guard"]
FAMILY_COPY_EVIDENCE_SOURCE_SCHEMA = TableSchema(
    name="family_copy_evidence_source",
    columns=FAMILY_COPY_ROUTE_A_BLAST_COLUMNS,
    required=frozenset(),
    enums={"route_a_blast_guard": ROUTE_A_BLAST_GUARD_REASONS},
    non_empty=frozenset(),
    unique_key=None,
)

FAMILY_CANDIDATE_BRIDGE_SCHEMA = TableSchema(
    name="family_candidate_bridge",
    columns=BRIDGE_COLUMNS,
    required=frozenset({"case_id", "gene_id", "family_id", "route", "ambiguity_status",
                        "provider_status", "mapping_tier", "rank_contribution"}),
    enums={
        "route": frozenset({"a", "b", "c"}),
        "match_type": BRIDGE_MATCH_TYPES,
        "ambiguity_status": BRIDGE_AMBIGUITY_STATUSES,
        "provider_status": BRIDGE_PROVIDER_STATUSES,
        "review_status": BRIDGE_REVIEW_STATUSES,
        "mapping_tier": BRIDGE_MAPPING_TIERS,
        "rank_contribution": frozenset({"none"}),
        "reason": BRIDGE_NOT_EVALUATED_REASONS,
    },
    non_empty=frozenset({"case_id", "gene_id", "family_id", "route", "ambiguity_status",
                         "provider_status", "mapping_tier", "rank_contribution"}),
    unique_key=None,
    coordinates=None,
)


# ---------------------------------------------------------------------------
# Pure-python validation (canonical, dependency-free)
# ---------------------------------------------------------------------------
def _series_str(frame: pd.DataFrame, column: str) -> pd.Series:
    return frame[column].astype(str)


def collect_errors(frame: pd.DataFrame, schema: TableSchema) -> list[SchemaError]:
    """Return every schema violation in ``frame`` (lazy, non-raising)."""

    errors: list[SchemaError] = []

    missing = [col for col in sorted(schema.required) if col not in frame.columns]
    for col in missing:
        errors.append(
            SchemaError(schema.name, "missing_column", f"required column '{col}' is absent", column=col)
        )
    if missing:
        # Row-level checks below assume the required columns exist.
        return errors

    for column in schema.non_empty:
        if column not in frame.columns:
            continue
        blank = _series_str(frame, column).str.strip().eq("") | frame[column].isna()
        for idx in frame.index[blank]:
            errors.append(
                SchemaError(schema.name, "empty_value", "value must be non-empty", column=column, row=int(idx))
            )

    for column, allowed in schema.enums.items():
        if column not in frame.columns:
            continue
        values = _series_str(frame, column)
        stripped = values.str.strip()
        bad = ~stripped.isin(allowed) & stripped.ne("")
        for idx in frame.index[bad]:
            errors.append(
                SchemaError(
                    schema.name,
                    "invalid_enum",
                    f"'{values[idx]}' is not one of {sorted(allowed)}",
                    column=column,
                    row=int(idx),
                )
            )

    if schema.unique_key and schema.unique_key in frame.columns:
        key = _series_str(frame, schema.unique_key)
        duplicated = key.duplicated(keep=False) & key.str.strip().ne("")
        for idx in frame.index[duplicated]:
            errors.append(
                SchemaError(
                    schema.name,
                    "duplicate_key",
                    f"duplicate {schema.unique_key} '{key[idx]}'",
                    column=schema.unique_key,
                    row=int(idx),
                )
            )

    if schema.coordinates:
        errors.extend(_coordinate_errors(frame, schema))

    return errors


def _coordinate_errors(frame: pd.DataFrame, schema: TableSchema) -> list[SchemaError]:
    spec = schema.coordinates
    assert spec is not None
    errors: list[SchemaError] = []
    for column in (spec.start, spec.end, spec.seqid):
        if column not in frame.columns:
            errors.append(
                SchemaError(schema.name, "missing_column", f"coordinate column '{column}' is absent", column=column)
            )
    if any(err.code == "missing_column" for err in errors):
        return errors

    starts = pd.to_numeric(frame[spec.start], errors="coerce")
    ends = pd.to_numeric(frame[spec.end], errors="coerce")
    seqids = _series_str(frame, spec.seqid).str.strip()
    floor = 1 if spec.one_based else 0

    for idx in frame.index:
        start = starts[idx]
        end = ends[idx]
        if pd.isna(start) or pd.isna(end):
            errors.append(
                SchemaError(schema.name, "invalid_coordinate", "start/end must be integers", column=spec.start, row=int(idx))
            )
            continue
        if start < floor or end < floor:
            errors.append(
                SchemaError(
                    schema.name,
                    "invalid_coordinate",
                    f"coordinates must be >= {floor}",
                    column=spec.start,
                    row=int(idx),
                )
            )
        if start > end:
            errors.append(
                SchemaError(
                    schema.name,
                    "invalid_coordinate",
                    f"start ({int(start)}) must be <= end ({int(end)})",
                    column=spec.start,
                    row=int(idx),
                )
            )
        if seqids[idx] == "":
            errors.append(
                SchemaError(schema.name, "empty_value", "seqid must be non-empty", column=spec.seqid, row=int(idx))
            )
    return errors


def validate(frame: pd.DataFrame, schema: TableSchema, *, raise_on_error: bool = True) -> list[SchemaError]:
    """Validate ``frame`` against ``schema``.

    Returns the list of :class:`SchemaError` entries. When ``raise_on_error`` is
    true (the default) and any error is found, raises
    :class:`SchemaValidationError` carrying the same list.
    """

    errors = collect_errors(frame, schema)
    if errors and raise_on_error:
        raise SchemaValidationError(schema.name, errors)
    return errors


# ---------------------------------------------------------------------------
# Optional pandera backend (derived from the same descriptor)
# ---------------------------------------------------------------------------
def pandera_available() -> bool:
    try:
        import pandera  # noqa: F401
    except Exception:
        return False
    return True


def pandera_schema(schema: TableSchema):
    """Derive a ``pandera.DataFrameSchema`` from ``schema``.

    Requires the optional ``pandera`` extra. The result is generated from the
    same :class:`TableSchema` descriptor the built-in validator uses, so the two
    views cannot drift.
    """

    import pandera as pa

    columns: dict[str, pa.Column] = {}
    for column in schema.columns:
        checks = []
        if column in schema.enums:
            checks.append(pa.Check.isin(sorted(schema.enums[column])))
        if column in schema.non_empty:
            checks.append(pa.Check.str_length(min_value=1))
        columns[column] = pa.Column(
            str,
            checks=checks or None,
            required=column in schema.required,
            unique=(column == schema.unique_key),
            nullable=column not in schema.non_empty,
            coerce=True,
        )

    dataframe_checks = []
    if schema.coordinates:
        spec = schema.coordinates
        floor = 1 if spec.one_based else 0

        def _coords_ok(df: pd.DataFrame, spec=spec, floor=floor) -> pd.Series:
            starts = pd.to_numeric(df[spec.start], errors="coerce")
            ends = pd.to_numeric(df[spec.end], errors="coerce")
            return starts.notna() & ends.notna() & (starts >= floor) & (starts <= ends)

        dataframe_checks.append(pa.Check(_coords_ok, element_wise=False, error="start<=end and >= floor"))

    return pa.DataFrameSchema(columns, checks=dataframe_checks or None, strict=False, name=schema.name)
