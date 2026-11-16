"""Family-to-copy disambiguation (issue #170).

Once the forward and reverse routes have recovered a plausible gene family, the
next question is plant-QTL specific: **which copy?** Paralogs, tandem duplications
and subfamily expansions mean family-level homology cannot, on its own, pick the
member to test first.

This stage compares the members of one family against each other and records what
does and does not discriminate between them. It is a tie-break / narrowing layer,
not a causal-probability estimator: ``preferred_copy`` means *first Wet / first
sequencing priority*, never *this is the causal gene*.

Deliberate boundaries:

* A discovery bucket alone never forms a copy set (issue #184). Copy
  discrimination runs only when the members of a family share a homology layer
  (same orthogroup/HOG, same curated homologous family, or the same domain
  architecture with provenance); with only a symbol-stem / naming bucket the
  family stays provisional (``family_resolution_required``) and no single
  ``preferred_copy`` is emitted. Members whose homology differs from the rest of
  the bucket are recorded as ``not_in_validated_copy_set`` and never compared as
  paralogs.
* Family-level support alone never promotes a single copy. With no discriminating
  axis the whole family stays ``family_supported_copy_unresolved``.
* Copy selection is **target-relative** (issue #185): the family-internal modal
  subfamily / modal domain architecture is descriptive statistics and outlier
  QC only -- it never picks the preferred copy. When a ``--family-target``
  artifact (verified / source-backed / curated) supplies an expected subfamily
  or domain architecture, copies are scored by how closely they match the
  target (``subfamily_vs_target`` / ``domain_vs_target``). Without a target the
  subfamily and domain axes are unresolved and no copy is promoted on them.
* A generic BLAST best hit is not a one-to-one ortholog. Orthology class comes
  only from a supplied orthology table; the family index alone yields
  ``family_level_only``.
* Multiple orthology providers agreeing is robustness, not extra evidence. An axis
  derived from the same evidence chain as another is counted once.
* Expression discriminates only when it *differs* between copies of the family,
  and ``unmeasured`` is never converted into ``not expressed``.
* Broad interval membership does not discriminate between copies, so it earns
  nothing.
* When tandem or paralog ambiguity cannot be resolved, no single copy is forced
  into first place.
* An unverified (``pending_review``) target never becomes the sole basis of a
  preferred copy, and a reference target duplicated in the focal genome
  (several copies matching it) never yields a sole preferred copy from target
  match alone -- locus-specific evidence is required.
* Unknown / orphan / dark genes are not removed at this stage.
* A supplied ``family_candidate_bridge`` mapping tier (issue #374) is recorded
  and gates shortlist promotion: a member whose bridge tier is
  ``tier5_name_only`` / ``tier6_not_evaluated`` is never promoted into the
  shortlist. The block is recorded in ``bridge_promotion_gate``
  (``blocked_tier5`` / ``blocked_tier6``) without changing any existing
  disambiguation status. The bridge is looked up by ``(gene_id, family_id)``,
  never by gene alone: a gene that also has a bridge row under a different
  family never inherits that other family's tier, route or dependency
  cluster (issue #374 T3 F1); with no row under this gene's own family the
  lookup reads as ``bridge_missingness="no_bridge_row"``. Route support
  (``rescue_a_only`` .. ``multi_route_supported``, or ``not_evaluated`` when
  every bridge row a family has is ``tier6_not_evaluated``) and bridge-backed
  copy-selection loss are reported in separate columns; a ``tier6_not_evaluated``
  row is a provider gap and never contributes a route (issue #374 T3 F2).
  The bridge columns are validated against the canonical vocabulary
  (``schemas.BRIDGE_MAPPING_TIERS``; route must be one of ``a``/``b``/``c``)
  fail-closed. Without a bridge the TSV outputs (members / relationships /
  discrimination / shortlist) are byte-identical; the manifest gains a bridge
  section whether or not a bridge was supplied.
* A supplied ``blast_evidence.tsv`` from the route A blast-panel search adds a
  ``route_a_blastp`` homology layer for qualifying rows. When the separate
  ``eligible_interval_hits.tsv`` is supplied, every threshold-passing target
  is adapted into that layer; it is never reduced to the compact top-hit
  summary (issue #532). Without either route A blast input the TSV outputs are
  byte-identical; the discrimination table gains the
  ``copy_set_evidence_source`` provenance column only when an input is
  supplied, and the manifest records the route_a_blast section either way.
* Issue #398 paralog guard: ``unique`` only means the single hit, not the right
  hit. When a qualifying hit member has an in-interval relative recorded under
  the same bridge family (same route, same ``family_id``; ``tier6_not_evaluated``
  bridge rows are a provider gap and never an assignment) and that relative is
  not one of the family's own members, the layer is not injected and the
  machine-readable reason is recorded in the ``route_a_blast_guard``
  discrimination column (present only when the input is supplied). The guard
  sits inside the layer qualification, below the bridge promotion gate.
"""

from __future__ import annotations

import json
import math

import pandas as pd

from .core import file_fingerprint, runtime_path
from .family_labels import canonical_family_key
from .schemas import (BRIDGE_MAPPING_TIERS, COPY_SET_EVIDENCE_SOURCES, FAMILY_COPY_ROUTE_A_BLAST_COLUMNS,
                      ROUTE_A_BLAST_GUARD_REASONS, SchemaError, SchemaValidationError)

SCHEMA_VERSION = "1.3"
DEFAULT_TANDEM_DISTANCE = 100_000

ORTHOLOGY_CLASSES = ("one_to_one", "one_to_many", "many_to_many", "close_paralog",
                     "family_level_only", "unmapped_ambiguous")
DISAMBIGUATION_STATES = ("preferred_copy", "co_preferred_tie", "family_supported_copy_unresolved",
                         "paralog_ambiguous", "tandem_ambiguous", "copy_weakened",
                         "insufficient_discriminating_evidence", "family_resolution_required",
                         "not_in_validated_copy_set")
# Target artifact columns. ``verification_status`` outside the
# trusted set (verified / source_backed / curated) is treated as
# ``pending_review``: an LLM-only or otherwise unverified target may inform the
# axes but never becomes the sole basis of a preferred copy.
TARGET_COLUMNS = ["target_id", "hypothesis_id", "reference_gene_id", "reference_species",
                  "expected_family_id", "expected_orthogroup", "expected_subfamily",
                  "expected_domain_architecture", "expected_synteny_anchor",
                  "expected_expression_context", "source_ids", "verification_status"]
VERIFIED_TARGET_STATUSES = frozenset({"verified", "source_backed", "curated"})
# Axis -> ordered tiers, best first. A copy is "leading" on an axis when it sits
# at the best tier *observed within its own family*, and the axis only counts when
# the family's members do not all share that tier.
AXIS_TIERS: dict[str, tuple[str, ...]] = {
    "orthology_class": ORTHOLOGY_CLASSES,
    "subfamily_match": ("target_match", "target_mismatch", "subfamily_unresolved"),
    "synteny_status": ("synteny_supported", "synteny_partial", "synteny_absent", "synteny_not_evaluated"),
    "domain_architecture_match": ("architecture_match", "architecture_partial", "architecture_mismatch", "architecture_not_evaluated"),
    "domain_vs_target": ("target_match", "target_mismatch", "target_not_evaluated"),
    "relevant_expression_status": ("expressed_in_relevant_context", "expressed_elsewhere_only",
                                   "not_detected_in_measured_context", "unmeasured"),
    "position_discrimination": ("recombinant_core", "peak_overlap", "support_interval",
                                "broad_interval_only", "position_not_evaluated"),
}
# Axes that decide copy selection. ``domain_architecture_match`` (the family
# modal comparison) is deliberately excluded: it is outlier / truncation QC
# only, never a selection axis.
SELECTION_AXES = ("orthology_class", "subfamily_match", "synteny_status", "domain_vs_target",
                  "relevant_expression_status", "position_discrimination")
# Axes whose values are defined only against the target artifact.
TARGET_DERIVED_AXES = frozenset({"subfamily_match", "domain_vs_target"})
# Tiers that carry no information and therefore never make an axis discriminating.
_UNINFORMATIVE = {"subfamily_unresolved", "synteny_not_evaluated", "architecture_not_evaluated",
                  "unmeasured", "broad_interval_only", "position_not_evaluated", "family_level_only",
                  "target_not_evaluated", "unmapped_ambiguous"}

MEMBER_COLUMNS = ["family_id", "candidate_gene_id", "chrom", "start", "end", "family_label",
                  "reference_gene", "tandem_block", "member_index", "notes"]
RELATIONSHIP_COLUMNS = ["family_id", "candidate_gene_id", "reference_gene", "orthology_class",
                        "orthology_providers", "evidence_chain", "subfamily", "subfamily_source",
                        "paralog_count", "notes"]
DISCRIMINATION_COLUMNS = [
    "family_id", "copy_set_id", "copy_set_status", "candidate_gene_id", "subfamily",
    "orthology_class", "orthology_confidence",
    "synteny_status", "domain_architecture_match", "relevant_expression_status", "position_discrimination",
    "paralog_ambiguity", "discriminating_axes", "leading_axes", "copy_specific_support",
    "copy_specific_counterevidence", "family_level_support", "disambiguation_status",
    "rank_contribution", "review_reason",
    # Issue #185 target-relative columns. Kept after the legacy columns so the
    # original schema stays positionally backward compatible.
    "target_id", "target_verification_status", "target_reference_gene", "target_subfamily",
    "target_domain_architecture", "subfamily_vs_target", "domain_vs_target",
    "domain_vs_family_modal", "subfamily_family_modal", "target_specific_axes",
    "selection_basis",
]
SHORTLIST_COLUMNS = ["family_id", "candidate_gene_id", "disambiguation_status", "leading_axes",
                     "subfamily", "orthology_class", "first_wet_priority", "rank_contribution", "notes"]
REVIEW_COLUMNS = ["reason", "target", "required_next_action", "allowed_choices"]

# Issue #374 bridge wiring: the canonical family_candidate_bridge.tsv supplies
# mapping tier / dependency cluster / ambiguity / route support per member gene.
# All bridge columns are appended to the row dicts only when a bridge is
# supplied, and the frame column lists below are extended accordingly, so the
# legacy schema stays byte-identical without one.
BRIDGE_TIER_ORDER = ("tier1_source_backed_unique", "tier2_source_backed_ambiguous",
                     "tier3_sequence_match_unique", "tier4_sequence_match_ambiguous",
                     "tier5_name_only", "tier6_not_evaluated")
BRIDGE_SOURCE_TIERS = frozenset(BRIDGE_TIER_ORDER[:4])
BRIDGE_BLOCKED_TIERS = frozenset({"tier5_name_only", "tier6_not_evaluated"})
BRIDGE_REQUIRED_COLUMNS = frozenset({"gene_id", "family_id", "route", "mapping_tier"})
BRIDGE_ROUTE_STATES = ("rescue_a_only", "rescue_b_only", "rescue_c_only", "multi_route_supported", "not_evaluated")
BRIDGE_DISCRIMINATION_COLUMNS = [
    "bridge_mapping_tier", "bridge_dependency_cluster_id", "bridge_ambiguity_status",
    "bridge_missingness", "bridge_promotion_gate", "route_rescue_state", "copy_selection_loss",
]
BRIDGE_SHORTLIST_COLUMNS = ["bridge_mapping_tier"]

# Issue #374 route A verified BLASTP layer: qualifying blast_evidence.tsv rows
# are supplied as a fourth homology layer and their provenance is recorded in
# the discrimination column ``copy_set_evidence_source`` (present only when the
# input is supplied, mirroring the bridge columns). The thresholds are
# independent defensive floors with the same values as the bridge's MIN_*
# fractions; the real pipeline's frozen panel thresholds are stricter and were
# already applied when blast_evidence.tsv was written.
ROUTE_A_BLAST_REQUIRED_COLUMNS = frozenset({
    "panel_id", "review_status", "provider_status", "top_hit_target_id", "ambiguity_status",
    "top_hit_identity", "top_hit_query_coverage", "top_hit_target_coverage",
})
ROUTE_A_BLAST_HITS_REQUIRED_COLUMNS = frozenset({
    "panel_id", "target_id", "pident", "query_coverage", "target_coverage", "evalue",
    "review_status", "provider_status", "ambiguity_status", "eligibility", "rank_contribution",
})
ROUTE_A_BLAST_DISCRIMINATION_COLUMNS = list(FAMILY_COPY_ROUTE_A_BLAST_COLUMNS)
ROUTE_A_BLAST_LAYER_PREFIX = "route_a_blast:"
ROUTE_A_BLAST_MIN_IDENTITY = 0.3
ROUTE_A_BLAST_MIN_COVERAGE = 0.5


def _frame(rows: list[dict[str, object]], columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame(rows)[columns] if rows else pd.DataFrame(columns=columns)


def _read(path_arg: str | None, required: set[str], label: str) -> pd.DataFrame:
    if not path_arg:
        return pd.DataFrame(columns=sorted(required))
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{label} table is missing columns: {', '.join(missing)}")
    return frame


def _normalize_route_a_blast_hits(frame: pd.DataFrame) -> pd.DataFrame:
    """Adapt the long-form eligible-hit artifact to the legacy layer shape.

    ``blast_evidence.tsv`` is intentionally a one-row-per-bait summary.  The
    sibling ``eligible_interval_hits.tsv`` is one row per threshold-passing
    interval target and must remain complete here: every target becomes a
    separate homology-layer row.  The compact summary remains the default when
    no long-form artifact is supplied.
    """
    if frame.empty:
        return pd.DataFrame(columns=sorted(ROUTE_A_BLAST_REQUIRED_COLUMNS))
    rows: list[dict[str, object]] = []
    for record in frame.to_dict("records"):
        if str(record.get("review_status", "")).strip().casefold() != "verified":
            continue
        if str(record.get("provider_status", "")).strip().casefold() != "evaluable":
            continue
        if str(record.get("rank_contribution", "")).strip().casefold() != "none":
            continue
        eligibility = str(record.get("eligibility", "eligible")).strip().casefold()
        if eligibility and eligibility not in {"eligible", "threshold_pass", "threshold_passing", "passed", "true", "1"}:
            continue
        target = str(record.get("target_id", "")).strip()
        panel = str(record.get("panel_id", "")).strip()
        if not panel or not target:
            continue
        rows.append({
            "panel_id": panel,
            "review_status": "verified",
            "provider_status": "evaluable",
            "top_hit_target_id": target,
            "ambiguity_status": str(record.get("ambiguity_status", "")).strip(),
            "top_hit_identity": str(record.get("pident", "")).strip(),
            "top_hit_query_coverage": str(record.get("query_coverage", "")).strip(),
            "top_hit_target_coverage": str(record.get("target_coverage", "")).strip(),
            "top_hit_evalue": str(record.get("evalue", "")).strip(),
            "top_hit_bitscore": str(record.get("bitscore", "")).strip(),
            "hit_count": str(record.get("hit_count", "")).strip() or "1",
            "paralog_count": str(record.get("paralog_count", "")).strip() or "0",
            "rank_contribution": str(record.get("rank_contribution", "none")).strip() or "none",
        })
    return pd.DataFrame(rows, columns=sorted(ROUTE_A_BLAST_REQUIRED_COLUMNS | {
        "top_hit_bitscore", "hit_count", "paralog_count", "rank_contribution"}))


def _tier_rank(axis: str, value: str) -> int:
    tiers = AXIS_TIERS[axis]
    return tiers.index(value) if value in tiers else len(tiers)


def _verification_status(value: object) -> str:
    """Normalize a target verification status; anything unverified is pending_review."""
    text = str(value or "").strip().casefold().replace("-", "_").replace(" ", "_")
    return text if text in VERIFIED_TARGET_STATUSES else "pending_review"


def _split_values(value: object) -> set[str]:
    return {part.strip() for part in str(value or "").split(";") if part.strip()}


def _read_targets(path_arg: str | None) -> list[dict[str, str]]:
    frame = _read(path_arg, {"target_id"}, "family target")
    return [{"target_id": str(row["target_id"]), **{column: str(row.get(column, "")) for column in TARGET_COLUMNS if column != "target_id"}}
            for row in frame.to_dict("records")]


def _target_order(target: dict[str, str]) -> tuple[int, str]:
    return (0 if _verification_status(target.get("verification_status", "")) != "pending_review" else 1,
            str(target.get("target_id", "")))


def _match_target(targets: list[dict[str, str]], hypothesis_id: str, family_id: str) -> dict[str, str] | None:
    """Best target for a family: hypothesis_id match, then expected_family_id.

    Within one specificity the trusted (verified / source-backed / curated)
    targets win over ``pending_review`` ones, with target_id as the
    deterministic tie-break.
    """
    if not targets:
        return None
    for specificity, column in ((str(hypothesis_id), "hypothesis_id"), (str(family_id), "expected_family_id")):
        if not specificity:
            continue
        candidates = [target for target in targets if str(target.get(column, "")).strip() == specificity]
        if candidates:
            return min(candidates, key=_target_order)
    return None


def _target_for_component(family_target: dict[str, str] | None, targets: list[dict[str, str]],
                          shared_layer_ids: str) -> dict[str, str] | None:
    """Family-level target, else a target whose expected_orthogroup is this copy set's layer."""
    if family_target is not None:
        return family_target
    if not targets:
        return None
    layer_ids = {value.strip() for value in shared_layer_ids.split(";") if value.strip()}
    candidates = [target for target in targets if str(target.get("expected_orthogroup", "")).strip() in layer_ids]
    return min(candidates, key=_target_order) if candidates else None


def _target_specific_axes(target: dict[str, str] | None) -> str:
    """The selection axes whose evaluation involved the target artifact."""
    if not target:
        return ""
    axes = []
    if _split_values(target.get("expected_subfamily")):
        axes.append("subfamily_match")
    if _split_values(target.get("expected_domain_architecture")):
        axes.append("domain_vs_target")
    if target.get("expected_synteny_anchor"):
        axes.append("synteny_status")
    if target.get("reference_gene_id"):
        axes.append("orthology_class")
    if target.get("expected_expression_context"):
        axes.append("relevant_expression_status")
    return ";".join(axes)


def _subfamily_vs_target(subfamily: str, target: dict[str, str] | None) -> str:
    if not target:
        return "subfamily_unresolved"
    expected = _split_values(target.get("expected_subfamily"))
    if not expected or not subfamily:
        return "subfamily_unresolved"
    return "target_match" if subfamily in expected else "target_mismatch"


def _domain_vs_target(architecture: str, target: dict[str, str] | None) -> str:
    if not target:
        return "target_not_evaluated"
    expected = _split_values(target.get("expected_domain_architecture"))
    if not expected or not architecture:
        return "target_not_evaluated"
    return "target_match" if architecture in expected else "target_mismatch"


def _target_fields(target: dict[str, str] | None) -> dict[str, str]:
    return {
        "target_id": str(target.get("target_id", "")) if target else "",
        "target_verification_status": _verification_status(target.get("verification_status", "")) if target else "",
        "target_reference_gene": str(target.get("reference_gene_id", "")) if target else "",
        "target_subfamily": str(target.get("expected_subfamily", "")) if target else "",
        "target_domain_architecture": str(target.get("expected_domain_architecture", "")) if target else "",
    }


def _orthology_records(frame: pd.DataFrame) -> dict[str, dict[str, object]]:
    """Per gene orthology class, providers and evidence chain.

    Providers are recorded but never counted as separate votes: the class is the
    strongest single relationship claimed, and the chain identifies the shared
    lineage of that claim so a derived axis cannot double count it.
    """
    records: dict[str, dict[str, object]] = {}
    if frame.empty:
        return records
    for gene_id, group in frame.groupby("gene_id", sort=True):
        classes = [str(value).strip().casefold().replace("-", "_") for value in group.get("relationship", pd.Series(dtype=str))]
        normalized = [value for value in classes if value in ORTHOLOGY_CLASSES]
        best = min(normalized, key=lambda value: ORTHOLOGY_CLASSES.index(value)) if normalized else "unmapped_ambiguous"
        providers = sorted({str(value) for value in group.get("provider", pd.Series(dtype=str)) if str(value)})
        chains = sorted({str(value) for value in group.get("evidence_chain", pd.Series(dtype=str)) if str(value)})
        records[str(gene_id)] = {
            "orthology_class": best, "orthology_providers": ";".join(providers),
            "evidence_chain": ";".join(chains) or (";".join(providers) or "orthology"),
            "reference_gene": ";".join(sorted({str(value) for value in group.get("reference_gene_id", pd.Series(dtype=str)) if str(value)})),
            "orthogroup": ";".join(sorted({str(value) for value in group.get("orthogroup", pd.Series(dtype=str)) if str(value)})),
            "paralog_count": str(group.get("paralog_count", pd.Series(dtype=str)).iloc[0]) if "paralog_count" in group.columns else "",
        }
    return records


def _keyed(frame: pd.DataFrame, key: str = "gene_id") -> dict[str, dict[str, object]]:
    return {str(row[key]): row for row in frame.to_dict("records")} if not frame.empty else {}


def _bridge_tier_rank(tier: str) -> int:
    return BRIDGE_TIER_ORDER.index(tier) if tier in BRIDGE_TIER_ORDER else len(BRIDGE_TIER_ORDER)


def _bridge_index(frame: pd.DataFrame) -> tuple[dict[tuple[str, str], dict[str, str]], dict[tuple[str, str], set[str]]]:
    """Best bridge row per ``(gene_id, family_id)`` (lowest mapping tier wins, file order breaks ties).

    Keyed by the pair, never by gene alone (issue #374 T3 F1): the bridge TSV
    carries one ``family_id`` column per row and a gene can legitimately have
    rows under several families (e.g. a route B/C axis keyed on domain
    architecture or orthogroup rather than the mechanism hypothesis's
    ``gene_family``). Keying by gene alone let one family's strongest tier leak
    into another family's promotion gate for the same gene; keying by the pair
    means a family with no bridge row of its own for a gene reads as
    ``bridge_missingness="no_bridge_row"`` and never inherits another family's
    tier, dependency cluster or ambiguity status.

    Returns ``(best_row_per_key, routes_per_key)``; ``mapping_tier`` is the
    strongest tier claimed by any route of the ``(gene_id, family_id)`` pair.
    The canonical vocabulary is enforced fail-closed (issue #374 review): a
    ``mapping_tier`` value outside ``schemas.BRIDGE_MAPPING_TIERS`` or a
    ``route`` value outside ``a/b/c`` is rejected instead of being silently
    carried, so a typo or an old vocabulary can never read as ``passed`` and
    get promoted.

    A ``tier6_not_evaluated`` row never contributes a route (issue #374 T3
    F2): it is a provider gap, not evidence that a route reached this gene, so
    it must never manufacture ``rescue_*_only`` / ``multi_route_supported``
    route support downstream. Tier1-tier5 rows (tier5_name_only included --
    a weak but real mapping) still contribute their route.
    """
    canonical_tiers = frozenset(BRIDGE_MAPPING_TIERS)
    invalid_tiers = sorted({str(value).strip() for value in frame["mapping_tier"]}
                           - canonical_tiers)
    if invalid_tiers:
        raise ValueError(
            "family candidate bridge mapping_tier values must be on the canonical "
            f"vocabulary {sorted(canonical_tiers)}; got {invalid_tiers}")
    invalid_routes = sorted({str(value).strip() for value in frame["route"]}
                            - {"a", "b", "c"})
    if invalid_routes:
        raise ValueError(
            f"family candidate bridge route values must be one of a/b/c; got {invalid_routes}")
    by_key: dict[tuple[str, str], dict[str, str]] = {}
    routes: dict[tuple[str, str], set[str]] = {}
    for row in frame.to_dict("records"):
        gene_id = str(row["gene_id"])
        fam = str(row.get("family_id", ""))
        tier = str(row.get("mapping_tier", ""))
        record_dict = {
            "mapping_tier": tier,
            "dependency_cluster_id": str(row.get("dependency_cluster_id", "")),
            "ambiguity_status": str(row.get("ambiguity_status", "")),
        }
        route = str(row.get("route", "")).strip()
        for k in {(gene_id, fam), (gene_id, canonical_family_key(fam))}:
            current = by_key.get(k)
            if current is None or _bridge_tier_rank(tier) < _bridge_tier_rank(str(current.get("mapping_tier", ""))):
                by_key[k] = record_dict
            if route and tier != "tier6_not_evaluated":
                routes.setdefault(k, set()).add(route)
    return by_key, routes


def _bridge_related_by_gene(bridge_frame: pd.DataFrame) -> dict[str, set[str]]:
    """In-interval genes recorded under the same bridge family as each gene.

    Two genes count as same-family when they share a non-``tier6`` bridge row
    under the same ``(route, family_id)`` pair -- the route-scoped namespace
    rule of ``_bridge_index`` (issue #374 T3 F1). ``tier6_not_evaluated`` rows
    are a provider gap, never a family assignment, so they are excluded. The
    map feeds the issue #398 paralog guard: a ``unique`` route A blast hit on a
    member that has an in-interval relative (outside the current family) is
    held back, because ``unique`` only means the hit was alone, not that the
    hit is the right copy.
    """
    by_key: dict[tuple[str, str], set[str]] = {}
    for row in bridge_frame.to_dict("records"):
        if str(row.get("mapping_tier", "")) == "tier6_not_evaluated":
            continue
        route = str(row.get("route", "")).strip()
        family_id = str(row.get("family_id", "")).strip()
        gene = str(row.get("gene_id", "")).strip()
        if not (route and family_id and gene):
            continue
        by_key.setdefault((route, family_id), set()).add(gene)
    related: dict[str, set[str]] = {}
    for genes in by_key.values():
        for gene in genes:
            related.setdefault(gene, set()).update(genes - {gene})
    return related


def _route_rescue_state(routes: set[str], has_bridge_coverage: bool) -> str:
    """Family-level route aggregation from the route-CONTRIBUTING bridge rows only.

    ``routes`` already excludes any route whose only bridge rows for this
    family are ``tier6_not_evaluated`` (issue #374 T3 F2): a provider gap
    must never be converted into route support. One route -> rescue_<route>_only;
    more than one -> multi_route_supported. When the family has no bridge
    coverage at all the state stays ``""`` (``bridge_missingness`` already
    reports that absence). When the family DOES have bridge rows but none of
    them contribute a route -- every row is tier6_not_evaluated -- the gap is
    kept visible as ``"not_evaluated"`` instead of collapsing to the same
    ``""`` as true absence.
    """
    if routes:
        return "multi_route_supported" if len(routes) > 1 else f"rescue_{next(iter(routes))}_only"
    return "not_evaluated" if has_bridge_coverage else ""


def _promotion_gate(bridge_active: bool, tier: str) -> str:
    """Gate value at a promotion point (preferred / co-preferred status).

    tier5/tier6 block shortlist entry; anything else passes -- a member with no
    bridge row at all is reported via ``bridge_missingness``, never blocked.
    """
    if not bridge_active:
        return ""
    if tier == "tier5_name_only":
        return "blocked_tier5"
    if tier == "tier6_not_evaluated":
        return "blocked_tier6"
    return "passed"


def _bridge_member_fields(bridge_active: bool, gene_id: str, family_id: str,
                          bridge_by_gene: dict[tuple[str, str], dict[str, str]],
                          promotion_gate: str, family_fields: dict[str, str]) -> dict[str, str]:
    """Bridge-recorded columns for one discrimination row; empty without a bridge.

    Looked up by ``(gene_id, family_id)`` (issue #374 T3 F1): a bridge row this
    gene carries under a different family is never read here, so it can never
    stand in for this family's own (missing) coverage.
    """
    if not bridge_active:
        return {}
    record = bridge_by_gene.get((gene_id, family_id)) or bridge_by_gene.get((gene_id, canonical_family_key(family_id)))
    return {
        "bridge_mapping_tier": str(record.get("mapping_tier", "")) if record else "",
        "bridge_dependency_cluster_id": str(record.get("dependency_cluster_id", "")) if record else "",
        "bridge_ambiguity_status": str(record.get("ambiguity_status", "")) if record else "",
        "bridge_missingness": "covered" if record else "no_bridge_row",
        "bridge_promotion_gate": promotion_gate,
        "route_rescue_state": str(family_fields.get("route_rescue_state", "")),
        "copy_selection_loss": str(family_fields.get("copy_selection_loss", "false")),
    }


def _subfamily(gene_id: str, orthology: dict[str, object], curated: dict[str, str], reference: str) -> tuple[str, str]:
    """Subfamily is derived, never invented.

    A curated map wins; otherwise the orthogroup is the finest available proxy.
    With neither, the answer is an honest ``subfamily_unresolved`` rather than a
    guess built from a gene-symbol stem, which cannot separate e.g. CLV1/BAM-like
    from PXY/TDR-like.
    """
    for key in (gene_id, reference):
        if key and key in curated:
            return curated[key], "curated_map"
    orthogroup = str(orthology.get("orthogroup", ""))
    if orthogroup:
        return orthogroup, "orthogroup"
    return "", "unresolved"


def _homology_layers(record: dict[str, object], ortholog: dict[str, object],
                     route_a_blast: tuple[str, ...] = ()) -> dict[str, set[str]]:
    """Homology layers carried by one member, from the atlas, orthology and route A tables.

    A layer counts only when it is an explicit evolutionary or sequence-domain
    assignment with provenance: an orthogroup/HOG id, a curated evolutionary
    family id (``evolutionary_family_id``), a domain architecture backed by a
    non-empty domain source, or a route A verified BLASTP hit. A symbol-stem
    naming bucket is not a layer, which is the whole point of the #184 gate.
    The ``route_a_blast`` key is always present; when nothing was supplied for
    the member it is an empty set, so downstream behaviour is unchanged.
    """
    orthogroups: set[str] = set()
    for value in (str(ortholog.get("orthogroup", "")) + ";" + str(record.get("orthogroup", ""))).split(";"):
        if value.strip():
            orthogroups.add(value.strip())
    curated = {str(value).strip() for value in str(record.get("evolutionary_family_id", "")).split(";") if str(value).strip()}
    architecture = str(record.get("domain_architecture", "")).strip()
    domain_arch = {architecture} if (
        architecture and str(record.get("domain_ids", "")).strip() and str(record.get("domain_source", "")).strip()
    ) else set()
    return {"orthogroups": orthogroups, "curated": curated, "domain_arch": domain_arch,
            "route_a_blast": set(route_a_blast)}


def _has_homology(layers: dict[str, set[str]]) -> bool:
    return bool(layers["orthogroups"] or layers["curated"] or layers["domain_arch"] or layers["route_a_blast"])


def _copy_set_components(gene_ids: list[str], layers: dict[str, dict[str, set[str]]]) -> list[set[str]]:
    """Connected components of members that share a homology layer.

    Two members are connected when they share an orthogroup, a curated
    evolutionary family id, or a domain architecture with provenance. Members
    that share nothing are their own singleton components, so a family bucket
    mixing several homologous families splits into one component per family.
    """
    parent = {gene_id: gene_id for gene_id in gene_ids}

    def find(gene_id: str) -> str:
        while parent[gene_id] != gene_id:
            parent[gene_id] = parent[parent[gene_id]]
            gene_id = parent[gene_id]
        return gene_id

    def union(left: str, right: str) -> None:
        root_left, root_right = find(left), find(right)
        if root_left != root_right:
            parent[root_right] = root_left

    by_value: dict[str, list[str]] = {}
    for gene_id in gene_ids:
        for values in layers[gene_id].values():
            for value in values:
                by_value.setdefault(value, []).append(gene_id)
    for members in by_value.values():
        if len(members) > 1:
            for other in members[1:]:
                union(members[0], other)
    components: dict[str, set[str]] = {}
    for gene_id in gene_ids:
        components.setdefault(find(gene_id), set()).add(gene_id)
    return list(components.values())


def _shared_layer_ids(component: set[str], layers: dict[str, dict[str, set[str]]]) -> str:
    """The layer ids present in every member of the component: the copy-set id."""
    members = list(component)
    shared: set[str] = set()
    for key in ("orthogroups", "curated", "domain_arch", "route_a_blast"):
        common = set.intersection(*[layers[member][key] for member in members]) if members else set()
        shared |= common
    return ";".join(sorted(shared))


def _own_layer_ids(gene_id: str, layers: dict[str, dict[str, set[str]]]) -> str:
    """The member's own homology layer ids (its own homologous family)."""
    return ";".join(sorted(value for values in layers[gene_id].values() for value in values))


_LAYER_SOURCE_LABELS = {
    "orthogroups": "atlas_orthogroup",
    "curated": "atlas_curated_family",
    "domain_arch": "atlas_domain",
    "route_a_blast": "route_a_blastp",
}


def _identity_fraction(value: str) -> float | None:
    """Parse a ``top_hit_identity`` cell (BLAST outfmt6 pident, a percent).

    The fraction form is the percent divided by 100, and the defensive floor
    is checked against that fraction. Non-positive, blank, unparseable or
    non-finite values fail closed (return None).
    """
    try:
        number = float(str(value).strip())
    except ValueError:
        return None
    if not math.isfinite(number):
        return None
    if number <= 0.0:
        return None
    return number / 100.0


def _coverage_fraction(value: str) -> float | None:
    """Parse a ``top_hit_*_coverage`` cell (span / query-or-target length, a 0-1 fraction).

    Values in [0, 1] are used as-is; values above 1 are read as legacy
    percentage cells and divided by 100 (never produced by the blast_evidence
    writer). Negative, blank, unparseable or non-finite values fail closed
    (return None).
    """
    try:
        number = float(str(value).strip())
    except ValueError:
        return None
    if not math.isfinite(number):
        return None
    if number < 0.0:
        return None
    if number > 1.0:
        number /= 100.0
    return number


def _route_a_blast_layers(frame: pd.DataFrame, member_gene_ids: list[str],
                          interval_related: dict[str, set[str]] | None = None,
                          ) -> tuple[dict[str, set[str]], int, dict[str, str]]:
    """Route A verified BLASTP rows that qualify as a homology layer, per hit member gene.

    Qualification is the full AND of the issue #374 route A rules: provider
    must be ``evaluable``, review ``verified`` (a skipped_review_required or
    synthetic row never injects), ambiguity ``unique``, the top hit must be one
    of this family's members, and identity / query coverage / target coverage
    must parse as numbers passing the defensive floors (identity is a BLAST
    pident percent read as >= 0.3 after dividing by 100; coverages are
    span/query-length fractions read as >= 0.5, with legacy values above 1
    read as percentages). Everything else --
    ambiguous_paralog, no_threshold_hit, not_evaluated, unparseable cells, a
    hit outside the family -- fails closed and injects nothing. The real
    pipeline's frozen panel thresholds are stricter and were already applied
    when blast_evidence.tsv was written.

    Issue #398 paralog guard: ``unique`` only means the single hit, not the
    right hit. When ``interval_related`` (the bridge same-family map) shows the
    hit member has an in-interval relative that is NOT one of this family's own
    members, the row is held back -- no layer is injected and the reason is
    recorded in the returned ``held`` map -- because the true copy may be
    missing from the annotation and only a paralog was hit (same treatment as
    ``ambiguous_paralog``). Relatives that ARE family members are covered by
    the component logic and can never promote a singleton, so they do not
    trigger the guard. Returns (layers, eligible_row_count, held_by_gene).
    """
    members = set(member_gene_ids)
    layers: dict[str, set[str]] = {}
    eligible = 0
    held: dict[str, str] = {}
    if frame.empty:
        return layers, eligible, held
    for row in frame.to_dict("records"):
        if str(row.get("provider_status", "")) != "evaluable":
            continue
        if str(row.get("review_status", "")) != "verified":
            continue
        if str(row.get("ambiguity_status", "")) != "unique":
            continue
        target = str(row.get("top_hit_target_id", "")).strip()
        if not target or target not in members:
            continue
        identity = _identity_fraction(str(row.get("top_hit_identity", "")))
        query_coverage = _coverage_fraction(str(row.get("top_hit_query_coverage", "")))
        target_coverage = _coverage_fraction(str(row.get("top_hit_target_coverage", "")))
        if identity is None or query_coverage is None or target_coverage is None:
            continue
        if identity < ROUTE_A_BLAST_MIN_IDENTITY:
            continue
        if query_coverage < ROUTE_A_BLAST_MIN_COVERAGE or target_coverage < ROUTE_A_BLAST_MIN_COVERAGE:
            continue
        if interval_related is None or target not in interval_related:
            held[target] = "paralog_guard_not_evaluated"
            continue
        paralogs = interval_related.get(target, set()) - members - {target}
        if paralogs:
            held[target] = "interval_paralog_present"
            continue
        eligible += 1
        panel_id = str(row.get("panel_id", "")).strip()
        if panel_id:
            layers.setdefault(target, set()).add(ROUTE_A_BLAST_LAYER_PREFIX + panel_id)
    return layers, eligible, held


def _copy_set_evidence_sources(gene_ids: list[str] | set[str],
                               layers: dict[str, dict[str, set[str]]]) -> str:
    """Provenance of one copy-set row: which homology layers supplied it.

    A source counts when every listed gene carries that layer (for a validated
    component row the component members, for a lone ``not_in_validated_copy_set``
    member the gene itself). The cell is the sorted ';'-joined subset of
    ``schemas.COPY_SET_EVIDENCE_SOURCES`` labels, or "" when no layer is present.
    """
    members = list(gene_ids)
    labels: set[str] = set()
    for key, label in _LAYER_SOURCE_LABELS.items():
        if members and all(bool(layers[member][key]) for member in members):
            labels.add(label)
    return ";".join(sorted(labels))


def _validate_copy_set_evidence_sources(rows: list[dict[str, object]]) -> None:
    """Fail-closed vocabulary gate for the ``copy_set_evidence_source`` cells.

    Every ';'-split token must be a member of
    ``schemas.COPY_SET_EVIDENCE_SOURCES`` (an empty cell is fine); anything
    else raises ``schemas.SchemaValidationError`` so an unknown provenance
    label can never reach the artifact. Enforced only while the route A blast
    input is active, i.e. exactly when the column exists.
    """
    allowed = COPY_SET_EVIDENCE_SOURCES
    errors: list[SchemaError] = []
    for index, row in enumerate(rows):
        for token in str(row.get("copy_set_evidence_source", "")).split(";"):
            if token and token not in allowed:
                errors.append(SchemaError(
                    "family_copy_discrimination", "invalid_enum",
                    f"copy_set_evidence_source token not on the canonical vocabulary {sorted(allowed)}: {token!r}",
                    column="copy_set_evidence_source", row=index))
    if errors:
        raise SchemaValidationError("family_copy_discrimination", errors)


def _validate_route_a_blast_guard(rows: list[dict[str, object]]) -> None:
    """Fail-closed vocabulary gate for the ``route_a_blast_guard`` cells.

    Every non-empty cell must be a member of
    ``schemas.ROUTE_A_BLAST_GUARD_REASONS``; anything else raises
    ``schemas.SchemaValidationError`` so an unknown hold reason can never reach
    the artifact. Enforced only while the route A blast input is active,
    mirroring ``_validate_copy_set_evidence_sources``.
    """
    allowed = ROUTE_A_BLAST_GUARD_REASONS
    errors: list[SchemaError] = []
    for index, row in enumerate(rows):
        value = str(row.get("route_a_blast_guard", ""))
        if value and value not in allowed:
            errors.append(SchemaError(
                "family_copy_discrimination", "invalid_enum",
                f"route_a_blast_guard value not on the canonical vocabulary {sorted(allowed)}: {value!r}",
                column="route_a_blast_guard", row=index))
    if errors:
        raise SchemaValidationError("family_copy_discrimination", errors)


def _tandem_blocks(members: list[dict[str, object]], distance: int) -> dict[str, int]:
    blocks: dict[str, int] = {}
    by_chrom: dict[str, list[tuple[int, str]]] = {}
    for member in members:
        start = str(member.get("start", ""))
        by_chrom.setdefault(str(member.get("chrom", "")), []).append((int(start) if start.isdigit() else 0, str(member["candidate_gene_id"])))
    index = 0
    for _, entries in sorted(by_chrom.items()):
        entries.sort()
        index += 1
        blocks[entries[0][1]] = index
        for (previous, _), (current, gene_id) in zip(entries, entries[1:]):
            if current - previous > distance:
                index += 1
            blocks[gene_id] = index
    return blocks


def family_to_copy(
    families_arg: str,
    out_dir_arg: str,
    atlas_arg: str | None = None,
    orthology_arg: str | None = None,
    synteny_arg: str | None = None,
    expression_arg: str | None = None,
    position_arg: str | None = None,
    subfamily_map_arg: str | None = None,
    target_arg: str | None = None,
    tandem_distance: int = DEFAULT_TANDEM_DISTANCE,
    case_id: str = "",
    bridge_arg: str | None = None,
    route_a_blast_arg: str | None = None,
    route_a_blast_hits_arg: str | None = None,
) -> dict[str, object]:
    """Compare the members of each recovered family and record what separates them."""
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    families = _read(families_arg, {"hypothesis_id", "member_interval_genes"}, "family hypothesis")
    atlas = _keyed(_read(atlas_arg, {"gene_id"}, "interval atlas"))
    orthology = _orthology_records(_read(orthology_arg, {"gene_id"}, "orthology"))
    synteny = _keyed(_read(synteny_arg, {"gene_id"}, "synteny"))
    expression = _keyed(_read(expression_arg, {"gene_id"}, "expression"))
    position = _keyed(_read(position_arg, {"gene_id"}, "position"))
    targets = _read_targets(target_arg)
    bridge_active = bool(bridge_arg)
    bridge_by_gene: dict[tuple[str, str], dict[str, str]] = {}
    bridge_routes: dict[tuple[str, str], set[str]] = {}
    bridge_related: dict[str, set[str]] = {}
    if bridge_active:
        bridge_frame = _read(bridge_arg, BRIDGE_REQUIRED_COLUMNS, "family candidate bridge")
        bridge_by_gene, bridge_routes = _bridge_index(bridge_frame)
        bridge_related = _bridge_related_by_gene(bridge_frame)
    route_a_blast_hits_active = bool(route_a_blast_hits_arg)
    route_a_blast_hits_frame = _read(route_a_blast_hits_arg, set(ROUTE_A_BLAST_HITS_REQUIRED_COLUMNS),
                                     "route A eligible blast hits")
    route_a_blast_active = bool(route_a_blast_arg or route_a_blast_hits_arg)
    route_a_blast_frame = _read(route_a_blast_arg, ROUTE_A_BLAST_REQUIRED_COLUMNS,
                                "route A blast evidence")
    if route_a_blast_hits_active:
        # The long-form artifact is authoritative for layer membership.  Do
        # not reduce it to the compact top-hit summary; every eligible target
        # remains available to the copy/paralog review.
        route_a_blast_frame = _normalize_route_a_blast_hits(route_a_blast_hits_frame)
    route_a_blast_eligible_rows = 0
    route_a_blast_injected_genes: set[str] = set()
    route_a_blast_held_paralog_genes: set[str] = set()
    route_a_blast_held_not_evaluated_genes: set[str] = set()
    curated_subfamily: dict[str, str] = {}
    if subfamily_map_arg:
        frame = _read(subfamily_map_arg, {"subfamily"}, "subfamily map")
        for row in frame.to_dict("records"):
            for key in ("gene_id", "reference_gene_id", "family_id"):
                if str(row.get(key, "")):
                    curated_subfamily[str(row[key])] = str(row["subfamily"])

    axes_supplied = {
        "orthology_class": bool(orthology_arg), "subfamily_match": bool(subfamily_map_arg or orthology_arg),
        "synteny_status": bool(synteny_arg), "domain_architecture_match": bool(atlas_arg),
        "domain_vs_target": bool(atlas_arg), "relevant_expression_status": bool(expression_arg),
        "position_discrimination": bool(position_arg), "target": bool(target_arg),
    }
    review_rows: list[dict[str, str]] = []
    member_rows: list[dict[str, object]] = []
    relationship_rows: list[dict[str, object]] = []
    discrimination_rows: list[dict[str, object]] = []
    shortlist_rows: list[dict[str, object]] = []
    resolution_required: set[str] = set()
    mixed_family: set[str] = set()

    for hypothesis in families.to_dict("records"):
        family_id = (
            str(hypothesis.get("family_id", "")).strip()
            or str(hypothesis.get("gene_family", "")).strip()
            or str(hypothesis["hypothesis_id"])
        )
        family_target = _match_target(targets, str(hypothesis.get("hypothesis_id", "")), family_id)
        family_target_axes = _target_specific_axes(family_target)
        target_applied = family_target is not None
        made_preferred = False
        gene_ids = [value.strip() for value in str(hypothesis["member_interval_genes"]).replace(",", ";").split(";") if value.strip()]
        if not gene_ids:
            continue
        family_routes = {
            route for gene_id in gene_ids
            for route in (bridge_routes.get((gene_id, family_id), set())
                          | bridge_routes.get((gene_id, canonical_family_key(family_id)), set()))
        }
        family_has_bridge_coverage = any(
            (gene_id, family_id) in bridge_by_gene or (gene_id, canonical_family_key(family_id)) in bridge_by_gene
            for gene_id in gene_ids
        )
        family_rescue_state = _route_rescue_state(family_routes, family_has_bridge_coverage)
        members: list[dict[str, object]] = []
        for order, gene_id in enumerate(sorted(gene_ids), 1):
            record = atlas.get(gene_id, {})
            members.append({
                "family_id": family_id, "candidate_gene_id": gene_id,
                "chrom": str(record.get("chrom", "")), "start": str(record.get("start", "")),
                "end": str(record.get("end", "")), "family_label": str(record.get("family_label", "")),
                "reference_gene": str(orthology.get(gene_id, {}).get("reference_gene", ""))
                or str(record.get("orthology_reference_gene", "")).split(";")[0],
                "tandem_block": 0, "member_index": order,
                "notes": "" if record else "gene not present in the supplied atlas; only family membership is known",
            })
        blocks = _tandem_blocks(members, tandem_distance)
        for member in members:
            member["tandem_block"] = blocks.get(str(member["candidate_gene_id"]), 0)
        member_rows.extend(members)

        no_axes = not any(axes_supplied.values())
        route_a_blast_layers, route_a_blast_eligible, route_a_blast_held = _route_a_blast_layers(
            route_a_blast_frame, gene_ids, bridge_related)
        route_a_blast_eligible_rows += route_a_blast_eligible
        for target_id, reason in route_a_blast_held.items():
            if reason == "interval_paralog_present":
                route_a_blast_held_paralog_genes.add(target_id)
            elif reason == "paralog_guard_not_evaluated":
                route_a_blast_held_not_evaluated_genes.add(target_id)
        route_a_blast_injected_genes.update(
            gene_id for gene_id, values in route_a_blast_layers.items() if values)
        layers = {
            str(member["candidate_gene_id"]): _homology_layers(atlas.get(str(member["candidate_gene_id"]), {}),
                                                               orthology.get(str(member["candidate_gene_id"]), {}),
                                                               tuple(route_a_blast_layers.get(str(member["candidate_gene_id"]), ())))
            for member in members
        }
        components = _copy_set_components([str(member["candidate_gene_id"]) for member in members], layers)
        component_of = {gene_id: component for component in components for gene_id in component}
        has_homology = {str(member["candidate_gene_id"]): _has_homology(layers[str(member["candidate_gene_id"])])
                        for member in members}
        single_component = len(components) == 1

        axis_values: dict[str, dict[str, str]] = {axis: {} for axis in AXIS_TIERS}
        chains: dict[str, dict[str, str]] = {}
        subfamily_by_gene: dict[str, tuple[str, str]] = {}
        for member in members:
            gene_id = str(member["candidate_gene_id"])
            record = atlas.get(gene_id, {})
            ortholog = orthology.get(gene_id, {})
            subfamily, subfamily_source = _subfamily(gene_id, ortholog, curated_subfamily, str(member["reference_gene"]))
            subfamily_by_gene[gene_id] = (subfamily, subfamily_source)
            relationship_rows.append({
                "family_id": family_id, "candidate_gene_id": gene_id,
                "reference_gene": str(member["reference_gene"]),
                "orthology_class": str(ortholog.get("orthology_class", "family_level_only")),
                "orthology_providers": str(ortholog.get("orthology_providers", "")),
                "evidence_chain": str(ortholog.get("evidence_chain", "")),
                "subfamily": subfamily, "subfamily_source": subfamily_source,
                "paralog_count": str(ortholog.get("paralog_count", "")),
                "notes": "providers are recorded for robustness; they are not independent votes",
            })
            architecture = str(record.get("domain_architecture", ""))
            axis_values["orthology_class"][gene_id] = str(ortholog.get("orthology_class", "family_level_only"))
            # Target-relative by design: without a target the
            # subfamily axis stays unresolved and can never drive selection.
            axis_values["subfamily_match"][gene_id] = _subfamily_vs_target(subfamily, family_target)
            axis_values["synteny_status"][gene_id] = str(synteny.get(gene_id, {}).get("synteny_status", "") or "synteny_not_evaluated")
            axis_values["domain_architecture_match"][gene_id] = (
                "architecture_not_evaluated" if not architecture else "architecture_match")
            axis_values["domain_vs_target"][gene_id] = _domain_vs_target(architecture, family_target)
            axis_values["relevant_expression_status"][gene_id] = str(
                expression.get(gene_id, {}).get("relevant_expression_status", "") or "unmeasured")
            axis_values["position_discrimination"][gene_id] = str(
                position.get(gene_id, {}).get("position_class", "") or "position_not_evaluated")
            chains[gene_id] = {
                "orthology_class": str(ortholog.get("evidence_chain", "")) or "orthology",
                "subfamily_match": "curated" if subfamily_source == "curated_map" else (str(ortholog.get("evidence_chain", "")) or "orthology"),
                # Synteny anchored through the same provider as orthology is the
                # same evidence chain, so it must not be counted twice.
                "synteny_status": str(synteny.get(gene_id, {}).get("provider", "")) or "synteny",
                "domain_architecture_match": "domain", "domain_vs_target": "domain",
                "relevant_expression_status": "expression",
                "position_discrimination": "position",
            }

        # Descriptive family-level modals: recorded as statistics,
        # never as a selection signal for the copies.
        family_architectures = [str(atlas.get(str(member["candidate_gene_id"]), {}).get("domain_architecture", "")) for member in members]
        family_modal_architecture = max((value for value in family_architectures if value), key=family_architectures.count, default="")
        family_subfamilies = [subfamily_by_gene[str(member["candidate_gene_id"])][0] for member in members]
        family_modal_subfamily = max((value for value in family_subfamilies if value), key=family_subfamilies.count, default="")

        # Issue #184 gate: only members that share a homology layer (same
        # orthogroup/HOG, same curated homologous family, or the same domain
        # architecture with provenance) belong to one validated copy set and get
        # compared as copies. Members connected only by a symbol-stem naming
        # bucket never are, and no single preferred_copy is made from a
        # bucket-only family.
        status_of: dict[str, str] = {}
        copy_set_of: dict[str, tuple[str, str]] = {}
        for member in members:
            gene_id = str(member["candidate_gene_id"])
            if no_axes:
                status_of[gene_id] = "insufficient_discriminating_evidence"
            elif has_homology[gene_id] and len(component_of[gene_id]) == 1 and single_component:
                status_of[gene_id] = "preferred_copy"
            elif has_homology[gene_id]:
                status_of[gene_id] = "not_in_validated_copy_set"
            else:
                status_of[gene_id] = "family_resolution_required"
            if status_of[gene_id] == "not_in_validated_copy_set":
                # This member has its own homologous family (its own layer ids)
                # but it is not part of this family's validated copy set, so the
                # copy-set label here is not_applicable.
                copy_set_of[gene_id] = (_own_layer_ids(gene_id, layers), "not_applicable")
            elif has_homology[gene_id]:
                copy_set_of[gene_id] = (_own_layer_ids(gene_id, layers), "validated")
            elif not atlas.get(gene_id):
                copy_set_of[gene_id] = ("", "unresolved")
            else:
                bucket = str(atlas.get(gene_id, {}).get("discovery_bucket_id", "")) or family_id
                copy_set_of[gene_id] = (bucket, "provisional")

        # Issue #374: copy-selection loss is a family-level state kept in a
        # column separate from route rescue: the bridge has a tier1-4 hit for at
        # least one member, yet family_to_copy could not narrow the copies
        # (existing family_resolution_required / insufficient_discriminating_evidence).
        family_bridge_fields: dict[str, str] = {}
        if bridge_active:
            family_bridge_fields = {
                "route_rescue_state": family_rescue_state,
                "copy_selection_loss": (
                    "true" if (
                        any(str(bridge_by_gene.get((gene_id, family_id), {}).get("mapping_tier", "")) in BRIDGE_SOURCE_TIERS
                            for gene_id in gene_ids)
                        and any(status_of[gene_id] in {"family_resolution_required", "insufficient_discriminating_evidence"}
                                for gene_id in gene_ids)
                    ) else "false"),
            }

        # Discrimination runs within each validated component (>= 2 members that
        # share a homology layer). Members outside the component are never
        # compared as paralogs.
        handled: set[str] = set()
        family_row_start = len(discrimination_rows)
        for component in components:
            comp = sorted(component)
            if len(comp) < 2 or no_axes:
                continue
            handled.update(comp)
            architectures = [str(atlas.get(gene_id, {}).get("domain_architecture", "")) for gene_id in comp]
            modal = max((value for value in architectures if value), key=architectures.count, default="")
            subfamilies = [subfamily_by_gene[gene_id][0] for gene_id in comp]
            modal_subfamily = max((value for value in subfamilies if value), key=subfamilies.count, default="")
            component_target = _target_for_component(family_target, targets, _shared_layer_ids(component, layers))
            if component_target is not None:
                target_applied = True
            target_axes = _target_specific_axes(component_target)
            subfamily_vs_target_by_gene: dict[str, str] = {}
            domain_vs_target_by_gene: dict[str, str] = {}
            target_match_count = 0
            for gene_id in comp:
                subfamily, _ = subfamily_by_gene[gene_id]
                # Target-relative subfamily: the family modal is
                # descriptive only and never a selection signal.
                subfamily_vs_target = _subfamily_vs_target(subfamily, component_target)
                axis_values["subfamily_match"][gene_id] = subfamily_vs_target
                subfamily_vs_target_by_gene[gene_id] = subfamily_vs_target
                architecture = str(atlas.get(gene_id, {}).get("domain_architecture", ""))
                domain_vs_target = _domain_vs_target(architecture, component_target)
                axis_values["domain_vs_target"][gene_id] = domain_vs_target
                domain_vs_target_by_gene[gene_id] = domain_vs_target
                # The modal comparison stays as outlier / truncation QC only.
                axis_values["domain_architecture_match"][gene_id] = (
                    "architecture_not_evaluated" if not modal or not architecture
                    else ("architecture_match" if architecture == modal else "architecture_mismatch"))
                if subfamily_vs_target == "target_match" or domain_vs_target == "target_match":
                    target_match_count += 1
            # The reference target is duplicated in the focal genome when more
            # than one copy matches it: target match alone then cannot single
            # out one copy (issue #185, requirement D).
            duplicated_target = target_match_count >= 2

            # Issue #282 hard positional gate: when at least one member of the
            # validated copy set lies inside the recombinant-defined core, only
            # core members may become preferred. Comparative axes may still rank
            # core members against each other, but no amount of comparative
            # support lets a non-core copy overtake the core.
            core_members = [gene_id for gene_id in comp
                            if axis_values["position_discrimination"][gene_id] == "recombinant_core"]
            # The selection pool is the promotion universe: top / best-count /
            # preferred / tie are decided inside it, never outside. With at
            # least two core members the whole comparison (which axes
            # discriminate, what the best tier is, who leads) also runs inside
            # the pool, so a non-core copy's values can neither steal an axis
            # from the core nor turn an axis discriminating against it. With a
            # single core member there is nothing to compare inside the pool,
            # so the descriptive family-wide comparison still runs -- the gate
            # alone decides promotion for the sole core.
            selection_pool = core_members if core_members else comp
            pool_set = set(selection_pool)
            comparison = selection_pool if len(selection_pool) >= 2 else comp

            discriminating: list[str] = []
            for axis in SELECTION_AXES:
                # Issue #281: an axis only discriminates between copies when at
                # least two members of the comparison are EVALUABLE on it and
                # their evaluable values differ. A lone evaluable member
                # compared only against unmeasured / *_not_evaluated / unresolved
                # values never leads, and missingness is never a comparative loss
                # vote for the other members.
                informative = [axis_values[axis][gene_id] for gene_id in comparison
                               if axis_values[axis][gene_id] not in _UNINFORMATIVE]
                if len(informative) >= 2 and len(set(informative)) > 1:
                    discriminating.append(axis)

            leaders: dict[str, list[str]] = {}
            for gene_id in comparison:
                leading: list[str] = []
                counted_chains: set[str] = set()
                for axis in discriminating:
                    best = min(_tier_rank(axis, axis_values[axis][other]) for other in comparison)
                    if _tier_rank(axis, axis_values[axis][gene_id]) != best:
                        continue
                    chain = chains[gene_id][axis]
                    if chain in counted_chains:
                        continue
                    counted_chains.add(chain)
                    leading.append(axis)
                leaders[gene_id] = leading
            best_count = max((len(value) for gene_id, value in leaders.items() if gene_id in pool_set), default=0)
            top = sorted(gene_id for gene_id, value in leaders.items()
                         if gene_id in pool_set and len(value) == best_count and best_count > 0)
            # Issue #282 rules 2/4: when the comparison cannot rank the
            # recombinant-core members, the core itself keeps the running -- a
            # sole core member stays the first-Wet priority candidate and
            # indistinguishable cores keep their tie, even with no
            # discriminating evidence at all. Comparative support may narrow
            # inside the core; it can never demote the core members out of the
            # running, and the core gate is a positional hold, not evidence for
            # or against any copy (issue #281 semantics are untouched).
            gate_fallback = False
            if not top and core_members:
                top = sorted(selection_pool)
                gate_fallback = True

            for gene_id in comp:
                leading = leaders.get(gene_id, [])
                ortholog_class = axis_values["orthology_class"][gene_id]
                member = next(item for item in members if str(item["candidate_gene_id"]) == gene_id)
                tandem_shared = sum(1 for other in comp if other != gene_id
                                    and next(item for item in members if str(item["candidate_gene_id"]) == other)["tandem_block"]
                                    == member["tandem_block"]) > 0
                paralog_ambiguous = ortholog_class in {"one_to_many", "many_to_many", "close_paralog", "unmapped_ambiguous"}
                # The gate-promoted top (sole core, or indistinguishable cores
                # held by the gate fallback) may exist with no discriminating
                # axis at all, so being in top is checked first: the core gate
                # keeps the core in the running even when nothing compares.
                if gene_id in top:
                    status = "preferred_copy" if len(top) == 1 else "co_preferred_tie"
                elif not discriminating:
                    status = ("tandem_ambiguous" if tandem_shared else
                              ("paralog_ambiguous" if paralog_ambiguous else "family_supported_copy_unresolved"))
                elif leading:
                    status = "family_supported_copy_unresolved"
                elif all(axis_values[axis][gene_id] in _UNINFORMATIVE for axis in discriminating):
                    # Issue #281: a member whose own values are all unmeasured /
                    # not-evaluated on the discriminating axes is unresolved, not
                    # weakened. Only an evaluable value can be counterevidence.
                    status = "family_supported_copy_unresolved"
                else:
                    status = "copy_weakened"
                # Issue #185 promotion guards: when the only axes that single
                # this copy out are target-derived, the target itself must be
                # verified and must not be duplicated in the focal genome,
                # otherwise no sole preferred copy is made.
                review_reason = ""
                if status == "preferred_copy":
                    target_only_leading = bool(leading) and all(axis in TARGET_DERIVED_AXES for axis in leading)
                    if target_only_leading and component_target is not None:
                        if _verification_status(component_target.get("verification_status", "")) == "pending_review":
                            status = "family_supported_copy_unresolved"
                            review_reason = ("target が pending_review（未検証 / LLM proposal）のため、target 由来の軸だけでは単独 preferred を作らない")
                            review_rows.append({
                                "reason": "copy_target_pending_review", "target": f"{family_id}:{gene_id}",
                                "required_next_action": "target を verified / source_backed / curated に検証してから再実行する（pending_review は単独 preferred の根拠にしない）",
                                "allowed_choices": "targetを検証する | 証拠を足す | 全copyを並行検証"})
                        elif duplicated_target:
                            status = "paralog_ambiguous"
                            review_reason = ("reference target が focal genome で duplicated（複数 copy が target に一致）のため、target match だけで単独 preferred を作らない")
                            review_rows.append({
                                "reason": "copy_target_duplicated", "target": f"{family_id}:{gene_id}",
                                "required_next_action": "synteny / position / expression / allele の locus-specific 証拠で複数 copy を分離する。不足なら co-preferred で進める",
                                "allowed_choices": "locus-specific証拠を足す | co-preferredで進める | family単位で保留"})
                if not review_reason and status in {"tandem_ambiguous", "paralog_ambiguous"}:
                    review_reason = "区別できる軸が無いため単独1位を作らない。追加の識別証拠が要る"
                    review_rows.append({"reason": f"copy_{status}", "target": f"{family_id}:{gene_id}",
                                        "required_next_action": "target解決（expected subfamily / domain / reference gene）か synteny / 対象組織のexpression / 組換え体位置のいずれかを足す",
                                        "allowed_choices": "証拠を足す | 全copyを並行検証 | family単位で保留"})
                # Issue #282 audit: members outside the core pool can never
                # become first-Wet candidates, and the gate is recorded in the
                # review_reason so the hold is auditable without inventing
                # evidence against them. When the gate itself promoted a core
                # (no discriminating axis), that too is recorded -- the core
                # gate is a positional hold, not a comparative vote.
                if core_members and gene_id not in pool_set:
                    if len(core_members) == 1:
                        gate_note = ("recombinant_core gate: この copy は組換え体コア（recombinant_core）の外にあり、"
                                     "唯一のコア copy が存在するため first-Wet 候補にならない。"
                                     "コアが単独の場合は比較は family 全体で行われる（コア内比較は不可）")
                    else:
                        gate_note = ("recombinant_core gate: この copy は組換え体コア（recombinant_core）の外にあり、"
                                     "コア内 copy が存在するため first-Wet 候補にならない。比較軸はコア内でしか数えない")
                    review_reason = f"{review_reason} {gate_note}".strip() if review_reason else gate_note
                elif gate_fallback and status in {"preferred_copy", "co_preferred_tie"} and not review_reason:
                    review_reason = ("recombinant_core gate: コア内の比較軸が識別できないため、"
                                     "組換え体コア copy を first-Wet 候補として保持する（軸投票でコアを覆さない）")
                status_of[gene_id] = status
                copy_set_of[gene_id] = (_shared_layer_ids(component, layers), "validated")
                counter = [f"{axis}={axis_values[axis][gene_id]}" for axis in discriminating
                           if axis_values[axis][gene_id] not in _UNINFORMATIVE
                           and _tier_rank(axis, axis_values[axis][gene_id]) > min(
                               _tier_rank(axis, axis_values[axis][other]) for other in comparison)]
                # Issue #374 promotion gate: recorded only at promotion points;
                # the block never rewrites the disambiguation status itself.
                promotion_gate = (_promotion_gate(bridge_active, str(bridge_by_gene.get((gene_id, family_id), {}).get("mapping_tier", "")))
                                  if status in {"preferred_copy", "co_preferred_tie"} else "")
                if status in {"preferred_copy", "co_preferred_tie"} and not promotion_gate.startswith("blocked"):
                    made_preferred = True
                discrimination_rows.append({
                    "family_id": family_id, "copy_set_id": copy_set_of[gene_id][0],
                    "copy_set_status": copy_set_of[gene_id][1], "candidate_gene_id": gene_id,
                    "subfamily": subfamily_by_gene[gene_id][0],
                    "orthology_class": ortholog_class,
                    "orthology_confidence": "high" if ortholog_class == "one_to_one" else ("low" if paralog_ambiguous else "medium"),
                    "synteny_status": axis_values["synteny_status"][gene_id],
                    "domain_architecture_match": axis_values["domain_architecture_match"][gene_id],
                    "relevant_expression_status": axis_values["relevant_expression_status"][gene_id],
                    "position_discrimination": axis_values["position_discrimination"][gene_id],
                    "paralog_ambiguity": "tandem" if tandem_shared else ("paralog" if paralog_ambiguous else "none"),
                    "discriminating_axes": ";".join(discriminating),
                    "leading_axes": ";".join(leading),
                    "copy_specific_support": ";".join(f"{axis}={axis_values[axis][gene_id]}" for axis in leading),
                    "copy_specific_counterevidence": ";".join(counter),
                    # Family membership is discovery support, kept apart from anything
                    # that separates this copy from its siblings.
                    "family_level_support": f"member of {family_id} ({len(members)} copies in interval)",
                    "disambiguation_status": status, "rank_contribution": "none", "review_reason": review_reason,
                    **_target_fields(component_target),
                    "subfamily_vs_target": subfamily_vs_target_by_gene[gene_id],
                    "domain_vs_target": domain_vs_target_by_gene[gene_id],
                    "domain_vs_family_modal": axis_values["domain_architecture_match"][gene_id],
                    "subfamily_family_modal": modal_subfamily,
                    "target_specific_axes": target_axes,
                    "selection_basis": "",
                    **({"copy_set_evidence_source": _copy_set_evidence_sources(comp, layers),
                        "route_a_blast_guard": route_a_blast_held.get(gene_id, "")}
                       if route_a_blast_active else {}),
                    **_bridge_member_fields(bridge_active, gene_id, family_id, bridge_by_gene, promotion_gate, family_bridge_fields),
                })
                if status in {"preferred_copy", "co_preferred_tie"} and not promotion_gate.startswith("blocked"):
                    shortlist_rows.append({
                        "family_id": family_id, "candidate_gene_id": gene_id, "disambiguation_status": status,
                        "leading_axes": ";".join(leading),
                        "subfamily": subfamily_by_gene[gene_id][0], "orthology_class": ortholog_class,
                        "first_wet_priority": "first" if status == "preferred_copy" else "co_first",
                        "rank_contribution": "none",
                        "notes": "first Wet / first sequencing priority, not a causal claim",
                        **({"bridge_mapping_tier": str(bridge_by_gene.get((gene_id, family_id), {}).get("mapping_tier", ""))}
                           if bridge_active else {}),
                    })

        for member in members:
            gene_id = str(member["candidate_gene_id"])
            if gene_id in handled:
                continue
            status = status_of[gene_id]
            if status == "family_resolution_required":
                resolution_required.add(family_id)
            elif status == "not_in_validated_copy_set":
                mixed_family.add(family_id)
            promotion_gate = (_promotion_gate(bridge_active, str(bridge_by_gene.get((gene_id, family_id), {}).get("mapping_tier", "")))
                              if status == "preferred_copy" else "")
            ortholog_class = axis_values["orthology_class"][gene_id]
            paralog_ambiguous = ortholog_class in {"one_to_many", "many_to_many", "close_paralog", "unmapped_ambiguous"}
            tandem_shared = sum(1 for other in members if other["tandem_block"] == member["tandem_block"]) > 1
            if status == "family_resolution_required":
                if not atlas.get(gene_id):
                    review_reason = "gene not present in the supplied atlas; only family membership is known"
                else:
                    review_reason = "homologous family が未確定（naming bucket のみ）のため copy 比較を実行しない"
            elif status == "not_in_validated_copy_set":
                review_reason = "この member は別の homologous family（orthogroup 等）に属し、この family bucket の copy set に含まれない"
            elif status == "insufficient_discriminating_evidence":
                review_reason = "識別用の入力が1つも供給されていない"
            else:
                review_reason = ""
            discrimination_rows.append({
                "family_id": family_id, "copy_set_id": copy_set_of[gene_id][0],
                "copy_set_status": copy_set_of[gene_id][1], "candidate_gene_id": gene_id,
                "subfamily": subfamily_by_gene[gene_id][0],
                "orthology_class": ortholog_class,
                "orthology_confidence": "high" if ortholog_class == "one_to_one" else ("low" if paralog_ambiguous else "medium"),
                "synteny_status": axis_values["synteny_status"][gene_id],
                "domain_architecture_match": axis_values["domain_architecture_match"][gene_id],
                "relevant_expression_status": axis_values["relevant_expression_status"][gene_id],
                "position_discrimination": axis_values["position_discrimination"][gene_id],
                "paralog_ambiguity": "tandem" if tandem_shared else ("paralog" if paralog_ambiguous else "none"),
                "discriminating_axes": "", "leading_axes": "",
                "copy_specific_support": "", "copy_specific_counterevidence": "",
                "family_level_support": f"member of {family_id} ({len(members)} copies in interval)",
                "disambiguation_status": status, "rank_contribution": "none", "review_reason": review_reason,
                **_target_fields(family_target),
                "subfamily_vs_target": _subfamily_vs_target(subfamily_by_gene[gene_id][0], family_target),
                "domain_vs_target": _domain_vs_target(
                    str(atlas.get(gene_id, {}).get("domain_architecture", "")), family_target),
                "domain_vs_family_modal": (
                    "architecture_not_evaluated"
                    if not family_modal_architecture or not str(atlas.get(gene_id, {}).get("domain_architecture", ""))
                    else ("architecture_match"
                          if str(atlas.get(gene_id, {}).get("domain_architecture", "")) == family_modal_architecture
                          else "architecture_mismatch")),
                "subfamily_family_modal": family_modal_subfamily,
                "target_specific_axes": family_target_axes,
                "selection_basis": "",
                **({"copy_set_evidence_source": (
                        _copy_set_evidence_sources([gene_id], layers)
                        if status in {"preferred_copy", "not_in_validated_copy_set"} else ""),
                    "route_a_blast_guard": route_a_blast_held.get(gene_id, "")}
                   if route_a_blast_active else {}),
                **_bridge_member_fields(bridge_active, gene_id, family_id, bridge_by_gene, promotion_gate, family_bridge_fields),
            })
            if status == "preferred_copy" and not promotion_gate.startswith("blocked"):
                made_preferred = True
                shortlist_rows.append({
                    "family_id": family_id, "candidate_gene_id": gene_id, "disambiguation_status": status,
                    "leading_axes": "", "subfamily": subfamily_by_gene[gene_id][0], "orthology_class": ortholog_class,
                    "first_wet_priority": "first", "rank_contribution": "none",
                    "notes": "first Wet / first sequencing priority, not a causal claim",
                    **({"bridge_mapping_tier": str(bridge_by_gene.get((gene_id, family_id), {}).get("mapping_tier", ""))}
                       if bridge_active else {}),
                })

        # Per-family selection basis: target-relative when a target
        # artifact was applied, locus_specific_only when copies were picked on
        # independent evidence without a target, unresolved_no_target otherwise.
        family_basis = ("target_relative" if target_applied
                        else ("locus_specific_only" if made_preferred else "unresolved_no_target"))
        for row in discrimination_rows[family_row_start:]:
            row["selection_basis"] = family_basis
        if len(members) >= 2 and not target_applied and not made_preferred:
            review_rows.append({
                "reason": "copy_target_missing", "target": family_id,
                "required_next_action": "verified target specification（expected_subfamily / expected_domain_architecture / reference_gene_id）か locus-specific 証拠（synteny / position / expression）を供給しないと copy は絞れない（family 内 modal は単独 preferred の根拠にしない）",
                "allowed_choices": "target表を足す | locus-specific証拠を足す | 全copyを並行検証",
            })

        if resolution_required:
            review_rows.append({
                "reason": "family_resolution_required", "target": family_id,
                "required_next_action": "orthogroup/HOG、curated family map、または domain+provenance で homologous family を確定してから copy 比較する（symbol stem / naming bucket のみでは copy 集合を定義しない）",
                "allowed_choices": "orthology表を足す | curated family mapを足す | domain表を足す | family単位で保留",
            })
        if family_id in mixed_family:
            review_rows.append({
                "reason": "family_bucket_mixes_homologous_families", "target": family_id,
                "required_next_action": "family仮説を homologous family（orthogroup 等）単位に分解し、同じ homologous family 内の copy だけを比較する",
                "allowed_choices": "family仮説を分解 | 全memberを並行検証 | family単位で保留",
            })

    if not any(axes_supplied.values()):
        review_rows.append({"reason": "no_discriminating_input_supplied", "target": case_id or "case",
                            "required_next_action": "orthology / synteny / subfamily map / expression / position のいずれかを供給しないとcopyは絞れない",
                            "allowed_choices": "証拠を足す | family単位のまま進む"})

    if route_a_blast_active:
        _validate_copy_set_evidence_sources(discrimination_rows)
        _validate_route_a_blast_guard(discrimination_rows)

    members_frame = _frame(member_rows, MEMBER_COLUMNS).sort_values(["family_id", "candidate_gene_id"], kind="stable")
    relationships = _frame(relationship_rows, RELATIONSHIP_COLUMNS).sort_values(["family_id", "candidate_gene_id"], kind="stable")
    discrimination = _frame(discrimination_rows,
                            DISCRIMINATION_COLUMNS
                            + (BRIDGE_DISCRIMINATION_COLUMNS if bridge_active else [])
                            + (ROUTE_A_BLAST_DISCRIMINATION_COLUMNS if route_a_blast_active else [])
                            ).sort_values(["family_id", "candidate_gene_id"], kind="stable")
    shortlist = _frame(shortlist_rows, SHORTLIST_COLUMNS + (BRIDGE_SHORTLIST_COLUMNS if bridge_active else [])
                       ).sort_values(["family_id", "candidate_gene_id"], kind="stable")
    reviews = _frame(review_rows, REVIEW_COLUMNS).drop_duplicates().sort_values(["reason", "target"], kind="stable")

    paths = {
        "family_copy_members": out / "family_copy_members.tsv",
        "family_copy_relationships": out / "family_copy_relationships.tsv",
        "family_copy_discrimination": out / "family_copy_discrimination.tsv",
        "family_copy_shortlist": out / "family_copy_shortlist.tsv",
        "family_copy_review_required": out / "family_copy_review_required.tsv",
    }
    for frame, path in ((members_frame, paths["family_copy_members"]), (relationships, paths["family_copy_relationships"]),
                        (discrimination, paths["family_copy_discrimination"]), (shortlist, paths["family_copy_shortlist"]),
                        (reviews, paths["family_copy_review_required"])):
        frame.to_csv(path, sep="\t", index=False, lineterminator="\n")

    manifest = {
        "schema_version": SCHEMA_VERSION, "case_id": case_id,
        "family_count": int(members_frame["family_id"].nunique()) if not members_frame.empty else 0,
        "family_member_count": int(len(members_frame)),
        "first_wet_candidate_count": int(len(shortlist)),
        "family_to_first_wet_reduction": round(len(members_frame) / len(shortlist), 2) if len(shortlist) else None,
        "status_counts": {str(k): int(v) for k, v in discrimination["disambiguation_status"].value_counts().sort_index().items()} if not discrimination.empty else {},
        "copy_set_status_counts": {str(k): int(v) for k, v in discrimination["copy_set_status"].value_counts().sort_index().items()} if not discrimination.empty else {},
        "homology_gate": {
            "basis": "a family's copy set is validated only when its members share a homology layer: the same orthogroup/HOG, the same curated homologous family (evolutionary_family_id), or the same domain architecture with provenance. A symbol-stem / naming discovery bucket alone never forms a copy set and never yields a single preferred_copy; such families require family resolution first (issue #184)",
            "copy_set_status_meaning": "validated: members share a homology layer and were compared as copies; provisional: discovery-bucket only, no copy set defined, family resolution required; not_applicable: member has its own homologous family but is outside this family's validated copy set; unresolved: no atlas record / dark member, kept without a claim",
            "families_resolution_required": sorted(resolution_required),
            "families_mixing_homologous_families": sorted(mixed_family),
        },
        "axes_supplied": axes_supplied,
        "selection_basis": ({str(row["family_id"]): str(row["selection_basis"]) for row in discrimination.to_dict("records")}
                            if not discrimination.empty else {}),
        "selection_basis_counts": {str(k): int(v) for k, v in discrimination["selection_basis"].value_counts().sort_index().items()} if not discrimination.empty else {},
        "target_relative": {
            "basis": "copy selection is scored against a verified / source-backed / curated target specification (expected subfamily, expected domain architecture, reference gene) instead of the family-internal modal majority; the modal subfamily / modal domain architecture remain descriptive statistics and outlier QC only. A pending_review target never becomes the sole basis of a preferred copy, and a target duplicated in the focal genome (several copies matching it) never yields a sole preferred copy from target match alone (issue #185)",
            "targets_applied": sorted({str(row["family_id"]) for row in discrimination.to_dict("records") if str(row.get("target_id", ""))}),
            "pending_review_targets_held_back": sorted({str(row["family_id"]) for row in discrimination.to_dict("records")
                                                        if str(row.get("target_verification_status", "")) == "pending_review"}),
        },
        "subfamily_resolution": {str(k): int(v) for k, v in relationships["subfamily_source"].value_counts().sort_index().items()} if not relationships.empty else {},
        "bridge": {
            "supplied": bridge_active,
            "bridge_tsv_sha256": (file_fingerprint(runtime_path(bridge_arg)).get("sha256") if bridge_active else None),
            "basis": "a supplied family_candidate_bridge mapping tier gates shortlist promotion: members with tier5_name_only / tier6_not_evaluated never enter the shortlist (bridge_promotion_gate=blocked_tier5 / blocked_tier6) while their disambiguation status is left unchanged; the bridge is looked up by (gene_id, family_id), never by gene alone, so a gene's bridge row under a different family can never promote or block it here; route support (rescue_a_only .. multi_route_supported, or not_evaluated when every bridge row a family has is tier6_not_evaluated) and copy_selection_loss are reported in separate columns",
            "blocked_count": int(sum(1 for row in discrimination_rows
                                     if str(row.get("bridge_promotion_gate", "")).startswith("blocked_"))),
            "bridge_rescue_counts": {state: int(len({str(row["family_id"]) for row in discrimination_rows
                                                     if row.get("route_rescue_state") == state}))
                                     for state in BRIDGE_ROUTE_STATES},
            "copy_selection_loss_count": int(len({str(row["family_id"]) for row in discrimination_rows
                                                  if row.get("copy_selection_loss") == "true"})),
        },
        "route_a_blast": {
            "supplied": route_a_blast_active,
            "blast_evidence_tsv_sha256": (file_fingerprint(runtime_path(route_a_blast_arg)).get("sha256")
                                          if route_a_blast_arg else None),
            "eligible_hits_supplied": route_a_blast_hits_active,
            "eligible_interval_hits_tsv_sha256": (file_fingerprint(runtime_path(route_a_blast_hits_arg)).get("sha256")
                                                   if route_a_blast_hits_active else None),
            "eligible_interval_hit_rows": int(len(route_a_blast_hits_frame)) if route_a_blast_hits_active else 0,
            "eligible_interval_hit_rows_consumed": int(len(route_a_blast_frame)) if route_a_blast_hits_active else 0,
            "evidence_rows": int(len(route_a_blast_frame)),
            "eligible_unique_rows": route_a_blast_eligible_rows,
            "injected_member_count": len(route_a_blast_injected_genes),
            "held_interval_paralog_gene_count": len(route_a_blast_held_paralog_genes),
            "held_guard_not_evaluated_gene_count": len(route_a_blast_held_not_evaluated_genes),
            # Issue #398: the guard reads the in-interval relative map out of the
            # bridge, so without evaluated non-tier6 bridge coverage it cannot run.
            # Record whether the check actually ran so a reader can tell "no paralog
            # found" apart from "no paralog looked for".
            "paralog_guard_evaluated": bool(bridge_active and bridge_related),
            "basis": "a supplied blast_evidence.tsv adds a route_a_blastp homology layer per qualifying row: provider_status=evaluable, review_status=verified, ambiguity_status=unique, the top hit is a member of the family in question, and identity / query coverage / target coverage pass the defensive floors (identity >= 0.3 with values >= 1 read as percentages; coverages >= 0.5). A qualifying row lets a single-member family reach preferred_copy and records its provenance in the copy_set_evidence_source discrimination column (present only when this input is supplied); ambiguous_paralog / no_threshold_hit / not_evaluated / skipped_review_required / below-threshold rows inject nothing, and the bridge promotion gate still wins over the layer. Issue #398 paralog guard: a unique row whose hit member has an in-interval relative recorded under the same bridge family (same route, same family_id; tier6_not_evaluated bridge rows are never an assignment) is held back -- no layer is injected and the reason is recorded in the route_a_blast_guard discrimination column -- because unique only means the single hit, not the right hit, when the true copy may be missing from the annotation. The guard needs the bridge relative map, so it fails closed when that map is unavailable for the hit gene: with no bridge input at all, or with only tier6_not_evaluated coverage for that gene, no layer is injected, the row is held with route_a_blast_guard=paralog_guard_not_evaluated, and no promotion is recorded. paralog_guard_evaluated reports whether the run had any evaluated (non-tier6) bridge coverage at all; per-gene evaluability is carried by the route_a_blast_guard column and the held_guard_not_evaluated_gene_count counter",
        },
        "review_required_count": int(len(reviews)),
        "meaning_of_preferred_copy": "first Wet / first sequencing priority; never a causal claim",
        "inputs": {name: (file_fingerprint(runtime_path(value)) if value else None) for name, value in (
            ("families", families_arg), ("atlas", atlas_arg), ("orthology", orthology_arg),
            ("synteny", synteny_arg), ("expression", expression_arg), ("position", position_arg),
            ("subfamily_map", subfamily_map_arg), ("target", target_arg), ("bridge", bridge_arg),
            ("route_a_blast", route_a_blast_arg),
            ("route_a_blast_hits", route_a_blast_hits_arg))},
        "safety": [
            "a discovery bucket alone (symbol stem / naming) never forms a validated copy set; bucket-only families stay provisional and produce no single preferred_copy (issue #184)",
            "members whose homology differs from the rest of the family bucket are never compared as paralogs and never promoted",
            "family-level support alone never promotes a single copy",
            "a generic best hit is not a one-to-one ortholog; orthology class comes only from a supplied table",
            "multiple orthology providers are robustness, not independent votes; one evidence chain is counted once",
            "expression discriminates only when it differs between copies, and unmeasured is never treated as not expressed",
            "broad interval membership does not discriminate between copies and earns nothing",
            "unresolved tandem or paralog ambiguity is preserved instead of forcing a single first place",
            "subfamily is derived from a curated map or an orthogroup, never invented from a symbol stem",
            "the family-internal modal subfamily / modal domain architecture are descriptive statistics and outlier QC only; they never pick the preferred copy (issue #185)",
            "copy selection is target-relative: subfamily and domain axes are scored against the target specification when one is supplied (issue #185)",
            "a pending_review (unverified / LLM-only) target never becomes the sole basis of a preferred copy (issue #185)",
            "a reference target duplicated in the focal genome (several copies matching it) never yields a sole preferred copy from target match alone; locus-specific evidence is required (issue #185)",
            "a route A unique hit whose member has an in-interval relative under the same bridge family is never promoted on that hit alone; the layer is held back because unique only means the single hit, not the right hit (issue #398)",
        ],
        "outputs": {name: str(path) for name, path in sorted(paths.items())},
    }
    manifest_path = out / "family_copy_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), "family_member_count": int(len(members_frame)),
            "first_wet_candidate_count": int(len(shortlist)),
            "family_to_first_wet_reduction": manifest["family_to_first_wet_reduction"],
            **{name: str(path) for name, path in paths.items()}}
