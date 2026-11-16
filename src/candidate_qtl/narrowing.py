"""Integrated candidate narrowing (issue #303).

Purpose: compress ``interval + phenotype`` into a Wet-reviewable candidate set.
The narrowing layer re-orders the already-ranked candidate list for Wet work by
combining *family recovery* (which families the phenotype/interval routes
recovered, and how strongly) with *copy discrimination* (which copy of a
validated family family-to-copy selected for first Wet work) and *survival
state* (which candidates candidate-elimination left standing).

The layer is a review-order / Wet-priority layer only:

- it never changes ``dry_score`` / ``dry_rank`` / ``dry_group`` /
  ``top_eligibility`` / ``screening_rank`` (every ranked column is passed
  through verbatim and nothing is recomputed),
- it produces ordinal tiers (N1..N4), never a numeric pseudo-probability,
- ``rank_contribution`` stays ``none`` on every row,
- it never excludes a candidate: only candidate-elimination (hard positional /
  genetic evidence) may exclude, and elimination alone may do that.

Family recovery is used as candidate-set narrowing / experimental routing only;
it is never causal evidence. ``not_evaluated`` / ``none`` recovery states mean
"no family information", never a negative, so they never demote a candidate:
tiers only move *up* when information exists.
"""

from __future__ import annotations

import json

import pandas as pd

from .core import file_fingerprint, runtime_path

SCHEMA_VERSION = "1.0"

# recovery_state vocabulary (ordinal). ``not_evaluated`` = no convergence input
# for this family; ``none`` = evaluated but no knowledge-route recovery.
RECOVERY_STATES = ("strict_convergent", "a_and_b_supported", "a_only", "b_only",
                   "none", "not_evaluated")
_RECOVERY_ORDINAL = {"strict_convergent": 5, "a_and_b_supported": 4, "a_only": 3,
                     "b_only": 2, "none": 1, "not_evaluated": 0}

# convergence_status values that record structural agreement between routes
# (retrieval_route.py / dual_route.py vocabularies).
CONVERGENT_STATUSES = frozenset({
    "a_b_c_convergent", "a_c_convergent", "b_c_convergent", "a_b_convergent",
    "convergent",
})

# disambiguation_status values that only arise when family-to-copy validated
# the family's copy set on a homology layer: these statuses imply
# the family is validated even when no convergence input was supplied.  The
# full vocabulary family_to_copy emits for validated sets is included, so a
# validated family is never silently demoted by the status vocabulary
# (tandem_ambiguous / paralog_ambiguous / copy_weakened all mean "validated
# family, copy not a single preferred choice").
VALIDATED_COPY_STATUSES = frozenset({
    "preferred_copy", "co_preferred_tie", "family_supported_copy_unresolved",
    "tandem_ambiguous", "paralog_ambiguous", "copy_weakened",
})
# family-to-copy also emits states that mean "we evaluated the bridge/copy
# question but could not resolve it".  These two states are deliberately NOT
# copy-only promotion states.  When an independent Route A/B recovery already
# supports the family, running family-to-copy must not make that recovered
# candidate worse than leaving copy discrimination unmeasured.
#
# ``not_in_validated_copy_set`` is intentionally excluded: issue #184 defines
# it as a member that does not belong to this validated homologous copy set.
# That rejected coarse family identity cannot, by itself, inherit N2 from the
# recovered family; the candidate stays in membership and falls back to other
# independent/base evidence.
RECOVERY_CONTEXT_COPY_UNRESOLVED_STATUSES = frozenset({
    "insufficient_discriminating_evidence",
    "family_resolution_required",
})
# Name-only and provider-gap bridge rows are not independent family support.
# They may retain a recovered candidate, but must not promote an unresolved
# copy state to N2.  Keep this local so legacy copy artifacts without the
# optional bridge column retain their existing semantics.
BRIDGE_NARROWING_BLOCKED_TIERS = frozenset({"tier5_name_only", "tier6_not_evaluated"})
# recovery states that validate the family on the recovery side.
VALIDATED_RECOVERY_STATES = frozenset({"strict_convergent", "a_and_b_supported"})
# recovery states that count as "family recovered" for tier assignment.
RECOVERED_RECOVERY_STATES = frozenset(
    {"strict_convergent", "a_and_b_supported", "a_only", "b_only"})

TIER_ORDER = ("N1", "N2", "N3", "N4")

# Structural narrowing_reason codes (frozen vocabulary).
NARROWING_REASONS = (
    "validated_family_preferred_copy",
    "validated_family_copy_unresolved",
    "family_recovered_copy_not_evaluated",
    "family_recovered_copy_unresolved",
    "base_top_reserve_rescue",
    "base_rank_only",
    "evidence_insufficient",
)

OUTPUT_SUFFIX_COLUMNS = (
    "family_recovery_state", "family_id", "copy_state", "first_wet_priority",
    "positional_constraint", "survival_state", "narrowing_tier",
    "narrowing_priority", "narrowing_reason", "unresolved_reason",
    "rank_contribution",
)

TIER_RULE_TEXT = (
    "tier is assigned left-to-right from the row's recovery x copy x base-rank state; "
    "only the first matching rule applies (ordinal tiers, never a numeric score). "
    "N1 = validated family (recovery_state in {strict_convergent, a_and_b_supported} "
    "or the family-copy input validated the copy set on a homology layer) AND "
    "copy_state == preferred_copy. "
    "N2 = validated family AND copy_state in the validated vocabulary other than "
    "preferred_copy (co_preferred_tie / family_supported_copy_unresolved / "
    "tandem_ambiguous / paralog_ambiguous / copy_weakened), or family recovered "
    "(strict_convergent / a_and_b_supported / a_only / b_only) AND copy_state is "
    "not_evaluated or an evaluated-but-unresolved context state "
    "(insufficient_discriminating_evidence / family_resolution_required). These "
    "states never validate or promote a family by themselves; they only preserve "
    "an independently recovered family. not_in_validated_copy_set is a copy-set "
    "identity exclusion, not an unresolved-state synonym, so it falls back to "
    "independent/base evidence rather than inheriting N2 from that family. An "
    "explicit bridge_mapping_tier of tier5_name_only or tier6_not_evaluated is "
    "a weak/provider-gap bridge and also blocks this unresolved-family N2 "
    "preservation; source-backed tiers retain the canonical recovered-family "
    "behavior. "
    "N3 = not N1/N2 and the base dry_group is in {top, reserve, rescue}. "
    "N4 = everything else (including base evidence_insufficient). "
    "When the family-copy shortlist is supplied, a member whose discrimination "
    "status is a promotion state (preferred_copy / co_preferred_tie) but that is "
    "absent from the shortlist was gate-blocked (e.g. bridge tier5/tier6) and is "
    "never N1/N2; it falls back to the base dry-group rules. "
    "not_evaluated / none family states carry no family information and never demote: "
    "they are not negatives, and a tier only moves up when information exists."
)

RECOVERY_MAPPING_TEXT = (
    "gene -> family comes from the family-recovery handoff (family_to_copy_input.tsv: "
    "member_interval_genes lists the family's interval genes). family -> recovery "
    "state comes from the convergence input: triple_route_hypothesis_comparison.tsv "
    "(four_state: a_and_b_supported / a_only_supported / b_only_supported / "
    "neither_supported) or, degraded, triple_route_convergence.json (only convergent "
    "families are resolvable from the JSON; everything else is not_evaluated). "
    "a_and_b_supported elevates to strict_convergent only when the convergence "
    "record shows source-independent agreement (convergence_status a_b_convergent / "
    "a_b_c_convergent AND source_dependence_status == source_independent). "
    "A shared-source-dependent convergence never elevates the recovery_state TO "
    "strict_convergent (it stays a_and_b_supported, which is a recovery "
    "(operational) state, not an independence claim: N routes citing the same "
    "reference provide one independent source, not N); a tier N1 placement may "
    "still follow from a validated preferred copy, which is an operational "
    "placement, not an independence claim. a_only_supported -> a_only; "
    "b_only_supported -> b_only; neither_supported / c_only_supported / "
    "gene_overlap_only / text_only_match -> none."
)

SAFETY_TEXT = (
    "narrowing is an experimental review-order / Wet-priority layer; it is not a "
    "causal claim, adds nothing to dry_score, and never excludes candidates "
    "(elimination alone may exclude)"
)

RECOVERY_DEGRADED_TEXT = (
    "a missing or empty optional input degrades that feature to not_evaluated; "
    "unmeasured is never treated as a negative"
)

FIRST_WET_PRIORITY_PRECEDENCE_TEXT = (
    "first_wet_priority is passed through from family_copy_shortlist.tsv "
    "(first / co_first) when present, else from remaining_candidate_set.tsv "
    "(high / medium / low), else empty"
)

POSITIONAL_CONSTRAINT_TEXT = (
    "positional_constraint is a review-order descriptor, never a score: it is "
    "passed through from the ranked table when the ranked table already carries "
    "the column, otherwise it is derived from the base dry_group "
    "(positionally_constrained for top/reserve/rescue, "
    "not_positionally_constrained for evidence_insufficient, empty when the "
    "ranked table has no dry_group)"
)


def _read_optional(path_arg: str | None, required: set[str], label: str) -> pd.DataFrame | None:
    """Read an optional narrowing input; an absent or empty file degrades to None."""
    if not path_arg:
        return None
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    if frame.empty:
        return None
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"{label} table is missing columns: {', '.join(missing)}")
    return frame


def _genes(text: str) -> set[str]:
    return {token.strip() for token in text.split(";") if token.strip()}


def _recovery_state_from_comparison(row: dict[str, object]) -> str:
    """Map one triple_route_hypothesis_comparison.tsv row to a recovery_state."""
    four_state = str(row.get("four_state", "")).strip()
    convergence_status = str(row.get("convergence_status", "")).strip()
    source_dependence = str(row.get("source_dependence_status", "")).strip()
    if four_state == "a_and_b_supported":
        if convergence_status in CONVERGENT_STATUSES and source_dependence == "source_independent":
            return "strict_convergent"
        return "a_and_b_supported"
    if four_state == "a_only_supported":
        return "a_only"
    if four_state == "b_only_supported":
        return "b_only"
    if four_state == "neither_supported":
        return "none"
    # Legacy artifacts without four_state fall back to the convergence_status.
    if convergence_status in CONVERGENT_STATUSES:
        return "strict_convergent" if source_dependence == "source_independent" else "a_and_b_supported"
    if convergence_status == "phenotype_only_supported":
        return "a_only"
    if convergence_status == "interval_only_supported":
        return "b_only"
    if convergence_status in {"c_only_supported", "gene_overlap_only", "text_only_match"}:
        return "none"
    return "not_evaluated"


def _read_convergence(path_arg: str | None) -> tuple[dict[str, str] | None, str]:
    """Resolve family -> recovery_state from the convergence input.

    Accepts ``triple_route_hypothesis_comparison.tsv`` (full per-family states)
    or ``triple_route_convergence.json`` (degraded: only convergent families are
    resolvable, and their a_only / b_only split cannot be recovered from the
    JSON alone). Returns (family_state_map, source_kind) where the map is None
    when the input is absent.
    """
    if not path_arg:
        return None, "absent"
    path = runtime_path(path_arg)
    if not path.exists():
        raise ValueError(f"convergence input does not exist: {path}")
    if path.suffix.casefold() == ".json":
        payload = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            raise ValueError("triple_route_convergence.json must contain an object")
        accounting = payload.get("convergence_source_accounting", [])
        if not isinstance(accounting, list) or not accounting:
            return {}, "json"
        states: dict[str, str] = {}
        for account in accounting:
            if not isinstance(account, dict):
                continue
            family = str(account.get("family_id", "")).strip()
            if not family:
                continue
            dependence = str(account.get("source_dependence_status", "")).strip()
            states[family] = ("strict_convergent" if dependence == "source_independent"
                              else "a_and_b_supported")
        return states, "json"
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if frame.empty:
        return {}, "tsv"
    if "family_id" not in frame.columns:
        raise ValueError("convergence table must contain a family_id column")
    states = {}
    for row in frame.to_dict("records"):
        family = str(row.get("family_id", "")).strip()
        if family:
            states.setdefault(family, _recovery_state_from_comparison(row))
    return states, "tsv"


def integrated_narrowing(
    ranked_arg: str,
    out_dir_arg: str,
    family_recovery_arg: str | None = None,
    convergence_arg: str | None = None,
    family_copy_arg: str | None = None,
    family_copy_shortlist_arg: str | None = None,
    remaining_arg: str | None = None,
    case_id: str = "",
    family_copy_manifest_arg: str | None = None,
) -> dict[str, object]:
    """Compress a ranked candidate set into an ordinal Wet-review order.

    ``ranked_arg`` (ranked_candidates.tsv) is required and every one of its
    columns is passed through verbatim, in the original column order, without
    recomputing dry_score / dry_rank / dry_group / top_eligibility /
    screening_rank. All other inputs are optional and fail-closed: a gene id
    that appears in any supplied input but not in the ranked table raises
    ValueError; a missing or empty optional input degrades its feature to
    ``not_evaluated`` (never a negative).
    """
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)

    ranked_path = runtime_path(ranked_arg)
    if not ranked_path.exists():
        raise ValueError(f"ranked candidates do not exist: {ranked_path}")
    ranked = pd.read_csv(ranked_path, sep="\t", dtype=str).fillna("")
    if "gene_id" not in ranked.columns:
        raise ValueError("ranked candidates require gene_id")
    if ranked["gene_id"].str.strip().eq("").any() or ranked["gene_id"].duplicated().any():
        raise ValueError("ranked candidate gene_id values must be non-empty and unique")
    if ranked.empty:
        raise ValueError("ranked candidates is empty")
    ranked_gene_ids = set(ranked["gene_id"].astype(str))

    family_recovery = _read_optional(
        family_recovery_arg, {"hypothesis_id", "gene_family", "member_interval_genes"},
        "family recovery handoff")
    family_copy = _read_optional(
        family_copy_arg, {"family_id", "candidate_gene_id", "disambiguation_status"},
        "family copy discrimination")
    family_copy_shortlist = _read_optional(
        family_copy_shortlist_arg, {"candidate_gene_id"}, "family copy shortlist")
    remaining = _read_optional(
        remaining_arg, {"candidate_gene_id", "current_gene_state"},
        "remaining candidate set")
    family_states, _convergence_kind = _read_convergence(convergence_arg)

    # Fail-closed: every gene id in any supplied input must exist in the ranked
    # table (mirrors next_experiment's unknown-remaining-gene guard).
    unknown: dict[str, list[str]] = {}
    if family_recovery is not None:
        for row in family_recovery.to_dict("records"):
            for gene in _genes(str(row.get("member_interval_genes", ""))):
                if gene not in ranked_gene_ids:
                    unknown.setdefault("family recovery handoff", []).append(gene)
    if family_copy is not None:
        for gene in family_copy["candidate_gene_id"].astype(str):
            if gene and gene not in ranked_gene_ids:
                unknown.setdefault("family copy discrimination", []).append(gene)
    if family_copy_shortlist is not None:
        for gene in family_copy_shortlist["candidate_gene_id"].astype(str):
            if gene and gene not in ranked_gene_ids:
                unknown.setdefault("family copy shortlist", []).append(gene)
    if remaining is not None:
        for gene in remaining["candidate_gene_id"].astype(str):
            if gene and gene not in ranked_gene_ids:
                unknown.setdefault("remaining candidate set", []).append(gene)
    if unknown:
        details = "; ".join(f"{label}: {', '.join(sorted(set(genes)))}"
                            for label, genes in sorted(unknown.items()))
        raise ValueError(f"narrowing inputs contain gene ids absent from the ranked table "
                         f"(fail-closed): {details}")

    # gene -> family resolution (family-recovery handoff only; the convergence
    # input never re-sources gene membership).
    gene_families: dict[str, set[str]] = {}
    if family_recovery is not None:
        for row in family_recovery.to_dict("records"):
            family = str(row.get("gene_family", "")).strip()
            if not family:
                continue
            for gene in _genes(str(row.get("member_interval_genes", ""))):
                gene_families.setdefault(gene, set()).add(family)

    # copy_state pass-through (family_copy_discrimination.disambiguation_status).
    copy_state: dict[str, str] = {}
    bridge_mapping_tier: dict[str, str] = {}
    if family_copy is not None:
        for row in family_copy.to_dict("records"):
            gene = str(row["candidate_gene_id"])
            copy_state[gene] = str(row["disambiguation_status"]).strip()
            # This column is optional on legacy discrimination artifacts.
            # Empty/absent preserves the pre-bridge unresolved-copy behavior.
            bridge_mapping_tier[gene] = str(row.get("bridge_mapping_tier", "")).strip()

    # Issue #374 promotion-gate reconciliation: family_copy_shortlist.tsv is the
    # promotion universe. family-to-copy keeps a gate-blocked member's
    # disambiguation_status (e.g. preferred_copy with bridge tier5/tier6) in the
    # discrimination artifact, but such a member is absent from the shortlist.
    # A promotion-status member that is missing from a supplied shortlist was
    # never promoted, so it must not claim the N1 preferred-copy tier (or the
    # N2 validated-copy tier); it falls back to the base dry-group rules. With
    # no shortlist input the legacy behavior is unchanged (copy-only inputs
    # still validate the family). A supplied shortlist file that is empty is
    # still a supplied, empty promotion universe: ``_read_optional`` would
    # degrade it to None, which would silently re-open the false-promotion
    # path for gate-blocked members, so the "was the shortlist supplied" fact
    # is tracked independently.
    shortlist_supplied = bool(family_copy_shortlist_arg)
    shortlist_gene_ids: set[str] = set()
    if shortlist_supplied and family_copy_shortlist is not None:
        shortlist_gene_ids = set(family_copy_shortlist["candidate_gene_id"].astype(str))
    promotion_blocked = {
        gene for gene, state in copy_state.items()
        if state in {"preferred_copy", "co_preferred_tie"} and gene not in shortlist_gene_ids
    } if shortlist_supplied else set()
    # When the promotion universe is supplied, an unresolved member outside
    # that universe cannot inherit N2 from a weak/provider-gap bridge.  With
    # no shortlist input, preserve legacy copy-artifact behavior.
    weak_bridge_unresolved = {
        gene for gene, state in copy_state.items()
        if state in RECOVERY_CONTEXT_COPY_UNRESOLVED_STATUSES
        and bridge_mapping_tier.get(gene, "") in BRIDGE_NARROWING_BLOCKED_TIERS
        and shortlist_supplied and gene not in shortlist_gene_ids
    }

    # first_wet_priority: shortlist first, remaining second (frozen precedence).
    first_wet_priority: dict[str, str] = {}
    if family_copy_shortlist is not None:
        for row in family_copy_shortlist.to_dict("records"):
            gene = str(row["candidate_gene_id"])
            if "first_wet_priority" in row and str(row.get("first_wet_priority", "")).strip():
                first_wet_priority.setdefault(gene, str(row["first_wet_priority"]).strip())
    if remaining is not None:
        for row in remaining.to_dict("records"):
            gene = str(row["candidate_gene_id"])
            if "first_wet_priority" in row and str(row.get("first_wet_priority", "")).strip():
                first_wet_priority.setdefault(gene, str(row["first_wet_priority"]).strip())

    # survival_state from remaining_candidate_set.current_gene_state.
    survival_state: dict[str, str] = {}
    remaining_gene_ids: set[str] = set()
    if remaining is not None:
        remaining_gene_ids = set(remaining["candidate_gene_id"].astype(str))
        survival_state = {str(row["candidate_gene_id"]): str(row["current_gene_state"]).strip()
                          for row in remaining.to_dict("records")}

    # Fail-closed: a ranked table that already carries a narrowing-owned column
    # would silently lose that base value (it is replaced by the narrowing
    # calculation), so refuse instead of overwriting (issue #303 review).
    # ``positional_constraint`` is the documented exception: the ranked table
    # may legitimately carry its own constraint and it is passed through
    # verbatim, never recomputed.  The internal ``_original_order`` sentinel is
    # protected too: a ranked column of that name would otherwise crash the
    # output shaping with a raw KeyError.
    collisions = sorted(
        set(ranked.columns) & (set(OUTPUT_SUFFIX_COLUMNS) | {"_original_order"})
        - {"positional_constraint"})
    if collisions:
        raise ValueError(
            f"ranked table already carries narrowing-owned columns (fail-closed): "
            f"{', '.join(collisions)}; integrated-narrowing never overwrites a base column")
    base_columns = list(ranked.columns)
    rows: list[dict[str, object]] = []
    for position, record in enumerate(ranked.to_dict("records")):
        gene_id = str(record["gene_id"])
        families = sorted(gene_families.get(gene_id, set()))
        family_id = ";".join(families)
        if families:
            if family_states is None:
                family_recovery_state = "not_evaluated"
            else:
                # A gene may belong to several families; take the strongest
                # recovery state by the explicit ordinal, never by string order.
                family_recovery_state = max(
                    (family_states.get(family, "not_evaluated") for family in families),
                    key=_RECOVERY_ORDINAL.get,
                )
        else:
            family_recovery_state = "not_evaluated"
        copy_state_value = copy_state.get(gene_id, "not_evaluated")
        dry_group = str(record.get("dry_group", "")).strip()
        if "positional_constraint" in record:
            # A ranked table that already carries its own positional constraint
            # column is passed through verbatim, never recomputed.
            positional_constraint = str(record.get("positional_constraint", ""))
        elif dry_group in {"top", "reserve", "rescue"}:
            positional_constraint = "positionally_constrained"
        elif dry_group == "evidence_insufficient":
            positional_constraint = "not_positionally_constrained"
        else:
            positional_constraint = ""
        if remaining is not None and remaining_gene_ids and gene_id not in remaining_gene_ids:
            survival_value = "not_in_remaining_set"
        else:
            survival_value = survival_state.get(gene_id, "not_evaluated")

        # Tier rules (ordinal; see TIER_RULE_TEXT). A family is validated by the
        # recovery state (strict_convergent / a_and_b_supported) or by the copy
        # input's validated statuses; preferred_copy itself implies a validated
        # copy set, so N1 is exactly copy_state == preferred_copy. A member the
        # promotion gate blocked (absent from the supplied shortlist despite a
        # promotion status) is never N1/N2 (promotion_blocked). Evaluated but
        # unresolved context states preserve N2 only when the family was already
        # recovered independently; they never promote a copy-only row. An
        # explicit tier5/tier6 bridge is too weak to preserve that N2 state.
        if copy_state_value == "preferred_copy" and gene_id not in promotion_blocked:
            narrowing_tier = "N1"
            narrowing_reason = "validated_family_preferred_copy"
        elif copy_state_value in VALIDATED_COPY_STATUSES and gene_id not in promotion_blocked:
            narrowing_tier = "N2"
            narrowing_reason = "validated_family_copy_unresolved"
        elif family_recovery_state in RECOVERED_RECOVERY_STATES and copy_state_value == "not_evaluated":
            narrowing_tier = "N2"
            narrowing_reason = "family_recovered_copy_not_evaluated"
        elif (family_recovery_state in RECOVERED_RECOVERY_STATES
              and copy_state_value in RECOVERY_CONTEXT_COPY_UNRESOLVED_STATUSES
              and gene_id not in weak_bridge_unresolved):
            narrowing_tier = "N2"
            narrowing_reason = "family_recovered_copy_unresolved"
        elif dry_group in {"top", "reserve", "rescue"}:
            narrowing_tier = "N3"
            narrowing_reason = "base_top_reserve_rescue"
        elif dry_group == "evidence_insufficient":
            narrowing_tier = "N4"
            narrowing_reason = "evidence_insufficient"
        else:
            narrowing_tier = "N4"
            narrowing_reason = "base_rank_only"

        if narrowing_reason == "validated_family_copy_unresolved":
            unresolved_reason = (
                f"family validated but copy discrimination resolved to a non-single state "
                f"(copy_state={copy_state_value}); no single preferred copy")
        elif narrowing_reason == "family_recovered_copy_not_evaluated":
            unresolved_reason = (
                f"family recovered (family_recovery_state={family_recovery_state}) but copy "
                f"discrimination was not evaluated")
        elif narrowing_reason == "family_recovered_copy_unresolved":
            unresolved_reason = (
                f"family recovered (family_recovery_state={family_recovery_state}) and copy "
                f"discrimination was evaluated but remains unresolved "
                f"(copy_state={copy_state_value}); unresolved is not negative evidence")
        elif narrowing_reason == "base_top_reserve_rescue":
            unresolved_reason = (
                "no validated family with a preferred copy; the base dry group is "
                "top/reserve/rescue only")
        elif narrowing_reason == "evidence_insufficient":
            unresolved_reason = (
                "no validated family with a preferred copy and the base dry evidence "
                "is insufficient")
        elif narrowing_reason == "base_rank_only":
            unresolved_reason = (
                "no validated family with a preferred copy; the base dry group is "
                "outside top/reserve/rescue")
        else:
            unresolved_reason = ""

        narrowed_row: dict[str, object] = {column: record.get(column, "") for column in base_columns}
        narrowed_row.update({
            "family_recovery_state": family_recovery_state,
            "family_id": family_id,
            "copy_state": copy_state_value,
            "first_wet_priority": first_wet_priority.get(gene_id, ""),
            "positional_constraint": positional_constraint,
            "survival_state": survival_value,
            "narrowing_tier": narrowing_tier,
            "narrowing_priority": 0,
            "narrowing_reason": narrowing_reason,
            "unresolved_reason": unresolved_reason,
            "rank_contribution": "none",
        })
        narrowed_row["_original_order"] = position
        rows.append(narrowed_row)

    narrowed = pd.DataFrame(rows)
    narrowed = narrowed.sort_values(
        ["narrowing_tier", "_original_order", "gene_id"], kind="stable"
    ).drop(columns="_original_order")
    narrowed["narrowing_priority"] = range(1, len(narrowed) + 1)

    remaining_base = [column for column in base_columns if column != "gene_id"]
    suffix_columns = [column for column in OUTPUT_SUFFIX_COLUMNS if column not in base_columns]
    output_columns = ["gene_id", *remaining_base, *suffix_columns]
    narrowed = narrowed[output_columns].reset_index(drop=True)

    narrowing_path = out / "integrated_narrowing.tsv"
    narrowed.to_csv(narrowing_path, sep="\t", index=False, lineterminator="\n")

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id,
        "candidate_count": int(len(narrowed)),
        "tier_counts": {tier: int((narrowed["narrowing_tier"] == tier).sum())
                        for tier in TIER_ORDER},
        "tier_rule": TIER_RULE_TEXT,
        "recovery_state_mapping": RECOVERY_MAPPING_TEXT,
        "meaning": (
            "integrated narrowing orders the ranked candidate set for Wet review by "
            "family recovery and copy discrimination; it is a review-order / "
            "Wet-priority layer, never a causal claim and never a score change"
        ),
        "narrowing_priority_meaning": (
            "narrowing_priority is the global 1..N experimental review order "
            "(tier-major, then the original ranked row order). It is an ordinal "
            "operational sequence for Wet handling, not a probability and not a "
            "replacement for dry_rank."
        ),
        "rank_contribution": (
            "none on every row; no numeric pseudo-probability is produced (ordinal "
            "tiers only), and dry_score / dry_rank / dry_group / top_eligibility / "
            "screening_rank are never generated, recomputed, or overwritten"
        ),
        "first_wet_priority_precedence": FIRST_WET_PRIORITY_PRECEDENCE_TEXT,
        "positional_constraint_meaning": POSITIONAL_CONSTRAINT_TEXT,
        "safety": [SAFETY_TEXT, RECOVERY_DEGRADED_TEXT],
        "inputs": {
            name: (file_fingerprint(runtime_path(value)) if value else None)
            for name, value in (
                ("ranked", ranked_arg), ("family_recovery", family_recovery_arg),
                ("convergence", convergence_arg), ("family_copy", family_copy_arg),
                ("family_copy_shortlist", family_copy_shortlist_arg),
                ("remaining", remaining_arg),
                ("family_copy_manifest", family_copy_manifest_arg),
            )
        },
        "outputs": {"integrated_narrowing": str(narrowing_path),
                    "narrowing_manifest": str(out / "narrowing_manifest.json")},
    }
    manifest_path = out / "narrowing_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
                             encoding="utf-8")
    return {
        "integrated_narrowing": str(narrowing_path),
        "narrowing_manifest": str(manifest_path),
        "candidate_count": int(len(narrowed)),
        "tier_counts": manifest["tier_counts"],
        "rank_contribution": "none",
    }
