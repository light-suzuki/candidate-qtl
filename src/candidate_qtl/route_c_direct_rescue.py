"""Route C direct trait-gene rescue contract (issue #640).

Route C is a source-backed *known-answer rescue / reference-evidence lane*.
It retrieves direct reported trait--gene/locus/variant relationships from the
focal species, close relatives, or model plants. It does not infer a new
mechanism or transfer a phenotype to a focal homolog by prose alone.

The central safety boundary is structural:

1. ``DIRECT_RESCUE_COLUMNS`` contain no QTL interval or focal-candidate fields.
   Direct reports are verified and frozen while interval-blind.
2. Only the separate ``PROJECTION_COLUMNS`` table may say whether a frozen
   reference record maps inside/outside the current QTL after exact-ID,
   assembly, orthology, BLAST, or family/copy projection.
3. Projection evidence is not itself trait-specific transfer evidence. A
   cross-species record remains a reference fact until existing machine
   verification gates justify a focal member.
4. ``no_interval_match`` and ``outside_current_qtl`` are different claims. A
   completed interval-only BLAST with no eligible hit does not prove that the
   focal homolog exists elsewhere in the genome. ``outside_current_qtl`` is
   reserved for a resolved focal gene/coordinate that is actually outside.

This module intentionally remains the schema/validator layer. The executable
Phase-2a search, verification, freeze, projection, and benchmark stages live in
``route_c_direct_pipeline``. Neither module adds supported-union membership, a
causal score, or Wet promotion.
"""

from __future__ import annotations

import pandas as pd

SCHEMA_VERSION = "0.1"
ROUTE_C_KIND = "direct_trait_gene_rescue"
RANK_CONTRIBUTION = "none"

# Search distance and evidence directness are deliberately separate axes. A D0
# Arabidopsis knockout can be highly informative while still requiring a
# cross-species projection before it says anything trait-specific about pea.
SPECIES_DISTANCE = frozenset({
    "S0_same_cultivar_or_accession",
    "S1_same_species",
    "S2_same_genus",
    "S3_close_clade",
    "S4_broader_legume",
    "S5_model_or_other_plant",
})

EVIDENCE_DIRECTNESS = frozenset({
    "D0_causal_perturbation_or_rescue",
    "D1_mapped_gene_with_functional_validation",
    "D2_direct_trait_gene_association",
    "D3_context_only",
})

SOURCE_VERIFICATION_STATES = frozenset({
    "verified",
    "pending",
    "failed",
    "not_evaluated",
})

GENE_IDENTITY_STATES = frozenset({
    "verified_gene",
    "verified_locus_only",
    "unresolved",
    "not_evaluated",
})

TRAIT_RELATION_STATES = frozenset({
    "direct_target_trait",
    "direct_related_trait",
    "context_only",
    "unresolved",
})

DIRECT_RESCUE_ELIGIBILITY = frozenset({
    "eligible_direct_record",
    "pending_verification",
    "context_only",
    "rejected",
})

# This exact-column contract is intentionally interval-blind. Focal species and
# cultivar are permitted because they define the biological query and make the
# S0..S5 distance auditable; QTL/candidate fields are not permitted here.
DIRECT_RESCUE_COLUMNS = (
    "record_id",
    "query_phenotype",
    "focal_species",
    "focal_cultivar_or_accession",
    "reported_phenotype",
    "reference_species",
    "reference_cultivar_or_accession",
    "reported_gene_id",
    "reported_gene_symbol",
    "reported_locus",
    "reported_variant",
    "species_distance",
    "evidence_directness",
    "evidence_type",
    "source_id",
    "pmid",
    "pmcid",
    "doi",
    "source_title",
    "source_year",
    "source_verification_state",
    "gene_identity_state",
    "trait_relation_state",
    "reference_sequence_id",
    "direct_rescue_eligibility",
    "rank_contribution",
    "notes",
)

PROJECTION_STATES = frozenset({
    "not_projected",
    "inside_current_qtl",
    "outside_current_qtl",
    "no_interval_match",
    "mapping_unresolved",
})

PROJECTION_BASES = frozenset({
    "not_evaluated",
    "exact_focal_id",
    "reported_coordinate",
    "assembly_or_alias_mapping",
    "orthology_reference_gene",
    "blastp",
    "family_or_copy_mapping",
    "unresolved",
})

PROJECTION_AMBIGUITY = frozenset({
    "not_evaluated",
    "unique",
    "one_to_many",
    "paralog_ambiguous",
    "no_threshold_hit",
    "unresolved",
})

PROJECTION_COLUMNS = (
    "record_id",
    "qtl_id",
    "focal_gene_id",
    "projection_state",
    "projection_basis",
    "source_reference_gene_id",
    "identity",
    "query_coverage",
    "target_coverage",
    "ambiguity_state",
    "interval_universe_fingerprint",
    "rank_contribution",
    "notes",
)


def _require_exact_columns(frame: pd.DataFrame, expected: tuple[str, ...], label: str) -> None:
    missing = sorted(set(expected) - set(frame.columns))
    extra = sorted(set(frame.columns) - set(expected))
    if missing:
        raise ValueError(f"{label} is missing columns: {', '.join(missing)}")
    if extra:
        raise ValueError(f"{label} has columns outside the contract: {', '.join(extra)}")


def _require_nonempty(frame: pd.DataFrame, columns: tuple[str, ...], label: str) -> None:
    for column in columns:
        if frame[column].astype(str).str.strip().eq("").any():
            raise ValueError(f"{label} requires non-empty {column}")


def _require_enum(frame: pd.DataFrame, column: str, allowed: frozenset[str], label: str) -> None:
    values = {str(value).strip() for value in frame[column]}
    invalid = sorted(values - allowed)
    if invalid:
        raise ValueError(
            f"{label} has invalid {column}: {', '.join(invalid)}; "
            f"allowed={', '.join(sorted(allowed))}"
        )


def validate_direct_rescue_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate an interval-blind frozen Route-C direct-report table.

    The function is strict by design: unknown columns are rejected so QTL or
    candidate-universe information cannot leak into the pre-freeze record set.
    It validates epistemic consistency but never turns a report into focal-gene
    support or a numerical score.
    """

    _require_exact_columns(frame, DIRECT_RESCUE_COLUMNS, "Route C direct-rescue table")
    if frame.empty:
        raise ValueError("Route C direct-rescue table is empty")

    normalized = frame.copy().fillna("")
    _require_nonempty(
        normalized,
        (
            "record_id",
            "query_phenotype",
            "focal_species",
            "reported_phenotype",
            "reference_species",
            "species_distance",
            "evidence_directness",
            "evidence_type",
            "source_id",
            "source_verification_state",
            "gene_identity_state",
            "trait_relation_state",
            "direct_rescue_eligibility",
            "rank_contribution",
        ),
        "Route C direct-rescue table",
    )

    normalized["record_id"] = normalized["record_id"].astype(str).str.strip()
    if normalized["record_id"].duplicated().any():
        raise ValueError("Route C direct-rescue record_id must be unique")

    identity_present = normalized[
        ["reported_gene_id", "reported_gene_symbol", "reported_locus", "reported_variant"]
    ].astype(str).apply(lambda row: row.str.strip().ne("").any(), axis=1)
    if not identity_present.all():
        raise ValueError(
            "Route C direct-rescue requires at least one reported gene/locus/variant identifier per row"
        )

    _require_enum(normalized, "species_distance", SPECIES_DISTANCE, "Route C direct-rescue table")
    _require_enum(normalized, "evidence_directness", EVIDENCE_DIRECTNESS, "Route C direct-rescue table")
    _require_enum(
        normalized,
        "source_verification_state",
        SOURCE_VERIFICATION_STATES,
        "Route C direct-rescue table",
    )
    _require_enum(normalized, "gene_identity_state", GENE_IDENTITY_STATES, "Route C direct-rescue table")
    _require_enum(normalized, "trait_relation_state", TRAIT_RELATION_STATES, "Route C direct-rescue table")
    _require_enum(
        normalized,
        "direct_rescue_eligibility",
        DIRECT_RESCUE_ELIGIBILITY,
        "Route C direct-rescue table",
    )

    if normalized["rank_contribution"].astype(str).str.strip().ne(RANK_CONTRIBUTION).any():
        raise ValueError("Route C direct-rescue rank_contribution must always be 'none'")

    for row in normalized.to_dict("records"):
        directness = str(row["evidence_directness"]).strip()
        source_state = str(row["source_verification_state"]).strip()
        identity_state = str(row["gene_identity_state"]).strip()
        trait_state = str(row["trait_relation_state"]).strip()
        eligibility = str(row["direct_rescue_eligibility"]).strip()

        if directness == "D3_context_only" and eligibility == "eligible_direct_record":
            raise ValueError("D3_context_only cannot be eligible_direct_record")
        if trait_state == "context_only" and eligibility == "eligible_direct_record":
            raise ValueError("context_only trait relation cannot be eligible_direct_record")
        if eligibility == "eligible_direct_record":
            if source_state != "verified":
                raise ValueError("eligible_direct_record requires source_verification_state=verified")
            if identity_state not in {"verified_gene", "verified_locus_only"}:
                raise ValueError("eligible_direct_record requires verified gene/locus identity")
            if trait_state not in {"direct_target_trait", "direct_related_trait"}:
                raise ValueError("eligible_direct_record requires a direct trait relation")
            if directness not in {
                "D0_causal_perturbation_or_rescue",
                "D1_mapped_gene_with_functional_validation",
                "D2_direct_trait_gene_association",
            }:
                raise ValueError("eligible_direct_record requires D0-D2 directness")

    return normalized


def validate_projection_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate the post-freeze Route-C reference -> focal/QTL projection table.

    ``outside_current_qtl`` is an evaluated mapped-gene result, not a biological
    negative and not a reason to delete the reference record. ``no_interval_match``
    is weaker: an interval-only search completed but found no eligible target;
    it never claims that a focal homolog exists elsewhere in the genome.
    """

    _require_exact_columns(frame, PROJECTION_COLUMNS, "Route C projection table")
    if frame.empty:
        raise ValueError("Route C projection table is empty")

    normalized = frame.copy().fillna("")
    _require_nonempty(
        normalized,
        ("record_id", "projection_state", "projection_basis", "ambiguity_state", "rank_contribution"),
        "Route C projection table",
    )

    _require_enum(normalized, "projection_state", PROJECTION_STATES, "Route C projection table")
    _require_enum(normalized, "projection_basis", PROJECTION_BASES, "Route C projection table")
    _require_enum(normalized, "ambiguity_state", PROJECTION_AMBIGUITY, "Route C projection table")

    if normalized["rank_contribution"].astype(str).str.strip().ne(RANK_CONTRIBUTION).any():
        raise ValueError("Route C projection rank_contribution must always be 'none'")

    for row in normalized.to_dict("records"):
        state = str(row["projection_state"]).strip()
        basis = str(row["projection_basis"]).strip()
        qtl_id = str(row["qtl_id"]).strip()
        focal_gene = str(row["focal_gene_id"]).strip()
        ambiguity = str(row["ambiguity_state"]).strip()

        if state in {"inside_current_qtl", "outside_current_qtl"}:
            if not qtl_id:
                raise ValueError(f"{state} requires qtl_id")
            if not focal_gene:
                raise ValueError(f"{state} requires focal_gene_id")
            if basis in {"not_evaluated", "unresolved"}:
                raise ValueError(f"{state} requires an evaluated projection_basis")
        elif state == "no_interval_match":
            if not qtl_id:
                raise ValueError("no_interval_match requires qtl_id")
            if focal_gene:
                raise ValueError("no_interval_match must not assert focal_gene_id")
            if basis in {"not_evaluated", "unresolved"}:
                raise ValueError("no_interval_match requires an evaluated interval-search basis")
            if ambiguity != "no_threshold_hit":
                raise ValueError("no_interval_match requires ambiguity_state=no_threshold_hit")
        elif state == "not_projected" and basis != "not_evaluated":
            raise ValueError("not_projected requires projection_basis=not_evaluated")
        elif state == "mapping_unresolved" and basis not in {"unresolved", "not_evaluated"}:
            raise ValueError("mapping_unresolved requires unresolved/not_evaluated projection_basis")

    return normalized


def direct_record_can_seed_projection(row: dict[str, object]) -> bool:
    """Return whether a verified direct record may enter post-freeze projection.

    This is intentionally *not* a Wet/support predicate. It only says that the
    literature/reference fact is strong enough to justify resolving its gene or
    locus against the focal genome using the existing projection machinery.

    A ``pending_verification`` record whose cited source resolved and named the
    reported gene, but whose claim cue was not detected, may also be projected
    (#988). Members keep it review-only because support still requires
    ``eligible_direct_record``.
    """

    identity = str(row.get("gene_identity_state", "")).strip()
    directness = str(row.get("evidence_directness", "")).strip()
    if identity not in {"verified_gene", "verified_locus_only"} or directness not in {
        "D0_causal_perturbation_or_rescue",
        "D1_mapped_gene_with_functional_validation",
        "D2_direct_trait_gene_association",
    }:
        return False
    eligibility = str(row.get("direct_rescue_eligibility", "")).strip()
    source_state = str(row.get("source_verification_state", "")).strip()
    if eligibility == "pending_verification":
        return source_state == "pending"
    return (
        eligibility == "eligible_direct_record"
        and source_state == "verified"
        and str(row.get("trait_relation_state", "")).strip()
        in {"direct_target_trait", "direct_related_trait"}
    )
