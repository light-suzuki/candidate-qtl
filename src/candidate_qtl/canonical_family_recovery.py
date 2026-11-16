"""Canonical supported-family recovery for issue #514.

This module projects the already-frozen ``supported_member_union.tsv`` into the
family-level handoff consumed by ``family_to_copy`` / ``integrated_narrowing``.
It deliberately does **not** re-read Route-A/Route-B candidate universes to make
membership decisions: every emitted member must already belong to the supported
union.

Family identity is cheap-first and bounded:

1. exact label equality;
2. mutually-unique deterministic canonical-family-key equality;
3. mutually-unique expanded Route-A-label -> Route-B-symbol matching;
4. ambiguous/nonmatching identities remain separate route-specific families.

No agreement bonus or numerical score is created.  The comparison artifact is
only a recovery-state adapter for the existing narrowing contract.

The optional bridge projection is intentionally conservative.  It preserves the
original Route-A mapping tier only for A-only families.  A final semantic Route-B
member has already passed source+transfer verification, which is a different
proof axis from the historical bridge tier; this module therefore never invents
a Route-B bridge tier and never lets a weak Route-A tier block an independently
verified Route-B-supported family.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable

import pandas as pd

from .core import file_fingerprint, runtime_path
from .family_labels import canonical_family_key, families_match_expanded

SCHEMA_VERSION = "1.0"

FAMILY_COLUMNS = [
    "hypothesis_id",
    "gene_family",
    "family_id",
    "member_interval_genes",
    "recovery_state",
    "family_match_basis",
    "route_support",
    "route_a_family_ids",
    "route_b_family_ids",
    "route_b_orthogroup_ids",
    "route_c_family_ids",
    "route_c_orthogroup_ids",
    "identity_ambiguity_status",
    "review_required",
    "rank_contribution",
]

COMPARISON_COLUMNS = [
    "comparison_id",
    "family_id",
    "convergence_status",
    "source_dependence_status",
    "four_state",
    "family_match_basis",
    "rank_contribution",
]

CANONICAL_BRIDGE_COLUMNS = [
    "gene_id",
    "family_id",
    "route",
    "mapping_tier",
    "dependency_cluster_id",
    "ambiguity_status",
    "source_family_id",
    "projection_basis",
    "rank_contribution",
]

UNION_REQUIRED = frozenset({
    "candidate_gene_id",
    "route_a_supported",
    "route_b_supported",
    "route_a_family_ids",
    "route_b_family_ids",
    "route_b_orthogroup_ids",
})

BRIDGE_REQUIRED = frozenset({"gene_id", "family_id", "route", "mapping_tier"})


def _tokens(value: object) -> list[str]:
    return sorted({token.strip() for token in str(value or "").replace(",", ";").split(";") if token.strip()})


def _truthy(value: object) -> bool:
    return str(value or "").strip().casefold() in {"1", "true", "yes", "y"}


def _norm(value: str) -> str:
    return " ".join(str(value).strip().casefold().split())


def _identity_maps(union: pd.DataFrame) -> tuple[dict[str, set[str]], dict[str, set[str]], dict[str, set[str]]]:
    route_a: dict[str, set[str]] = {}
    route_b_family: dict[str, set[str]] = {}
    route_b_orthogroup_only: dict[str, set[str]] = {}
    for row in union.to_dict("records"):
        gene = str(row["candidate_gene_id"]).strip()
        if _truthy(row.get("route_a_supported", "")):
            for family in _tokens(row.get("route_a_family_ids", "")):
                route_a.setdefault(family, set()).add(gene)
        if _truthy(row.get("route_b_supported", "")):
            b_families = _tokens(row.get("route_b_family_ids", ""))
            if b_families:
                for family in b_families:
                    route_b_family.setdefault(family, set()).add(gene)
            else:
                for orthogroup in _tokens(row.get("route_b_orthogroup_ids", "")):
                    route_b_orthogroup_only.setdefault(orthogroup, set()).add(gene)
    return route_a, route_b_family, route_b_orthogroup_only


def _route_c_identity_maps(union: pd.DataFrame) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    families: dict[str, set[str]] = {}
    orthogroups: dict[str, set[str]] = {}
    for row in union.to_dict("records"):
        if not _truthy(row.get("route_c_supported", "")):
            continue
        gene = str(row["candidate_gene_id"]).strip()
        family_ids = _tokens(row.get("route_c_family_ids", ""))
        for family in family_ids:
            families.setdefault(family, set()).add(gene)
        if not family_ids:
            for orthogroup in _tokens(row.get("route_c_orthogroup_ids", "")):
                orthogroups.setdefault(orthogroup, set()).add(gene)
    return families, orthogroups


def _unique_pairs(
    a_values: set[str],
    b_values: set[str],
    predicate: Callable[[str, str], bool],
) -> tuple[list[tuple[str, str]], set[str], set[str]]:
    """Return mutually-unique pairs plus nodes participating in ambiguity."""
    a_candidates = {a: sorted(b for b in b_values if predicate(a, b)) for a in sorted(a_values)}
    b_candidates = {b: sorted(a for a in a_values if predicate(a, b)) for b in sorted(b_values)}
    pairs: list[tuple[str, str]] = []
    ambiguous_a: set[str] = set()
    ambiguous_b: set[str] = set()
    for a, candidates in a_candidates.items():
        if not candidates:
            continue
        if len(candidates) == 1 and len(b_candidates[candidates[0]]) == 1:
            pairs.append((a, candidates[0]))
        else:
            ambiguous_a.add(a)
            ambiguous_b.update(candidates)
    for b, candidates in b_candidates.items():
        if len(candidates) > 1:
            ambiguous_b.add(b)
            ambiguous_a.update(candidates)
    return pairs, ambiguous_a, ambiguous_b


def _family_row(
    *,
    case_id: str,
    family_id: str,
    members: set[str],
    route_support: str,
    match_basis: str,
    route_a_ids: list[str],
    route_b_ids: list[str],
    route_b_orthogroups: list[str],
    ambiguous: bool,
) -> dict[str, str]:
    payload = "|".join([
        case_id,
        family_id,
        route_support,
        ";".join(sorted(route_a_ids)),
        ";".join(sorted(route_b_ids)),
        ";".join(sorted(route_b_orthogroups)),
        ";".join(sorted(members)),
    ])
    hypothesis_id = "canonical-family-" + hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]
    recovery_state = {
        "a+b": "a_and_b_supported",
        "a": "a_only",
        "b": "b_only",
    }.get(route_support, f"{route_support.replace('+', '_')}_supported")
    return {
        "hypothesis_id": hypothesis_id,
        "gene_family": family_id,
        "family_id": family_id,
        "member_interval_genes": ";".join(sorted(members)),
        "recovery_state": recovery_state,
        "family_match_basis": match_basis,
        "route_support": route_support,
        "route_a_family_ids": ";".join(sorted(route_a_ids)),
        "route_b_family_ids": ";".join(sorted(route_b_ids)),
        "route_b_orthogroup_ids": ";".join(sorted(route_b_orthogroups)),
        "route_c_family_ids": "",
        "route_c_orthogroup_ids": "",
        "identity_ambiguity_status": "ambiguous_unmerged" if ambiguous else "resolved_or_route_specific",
        "review_required": "yes" if ambiguous else "no",
        "rank_contribution": "none",
    }


def _comparison_row(family: dict[str, str]) -> dict[str, str]:
    route_support = family["route_support"]
    four_state = {
        "a+b": "a_and_b_supported",
        "a": "a_only_supported",
        "b": "b_only_supported",
    }.get(route_support, f"routes_{route_support.replace('+', '_')}_supported")
    convergence_status = {
        "a+b": "a_b_supported",
        "a": "phenotype_only_supported",
        "b": "interval_only_supported",
    }.get(route_support, f"{route_support.replace('+', '_')}_supported")
    return {
        "comparison_id": "canonical-comparison-" + hashlib.sha256(
            family["hypothesis_id"].encode("utf-8")
        ).hexdigest()[:16],
        "family_id": family["family_id"],
        "convergence_status": convergence_status,
        # Same-family recovery across A/B is provenance here, not proof that
        # the underlying biological sources are independent.
        "source_dependence_status": "not_evaluated" if route_support == "a+b" else "not_applicable",
        "four_state": four_state,
        "family_match_basis": family["family_match_basis"],
        "rank_contribution": "none",
    }


def _project_a_only_bridge(
    families: list[dict[str, str]],
    bridge_arg: str | None,
) -> pd.DataFrame:
    """Project only original A-only bridge rows onto the canonical family id.

    Route-B semantic verification is not translated into the old bridge-tier
    vocabulary.  For A+B/B families no historical mapping-tier gate is emitted;
    their support is represented by the canonical recovery/verified-B artifacts.
    """
    if not bridge_arg:
        return pd.DataFrame(columns=CANONICAL_BRIDGE_COLUMNS)
    path = runtime_path(bridge_arg)
    if not path.exists():
        raise ValueError(f"family bridge does not exist: {path}")
    bridge = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    missing = sorted(BRIDGE_REQUIRED - set(bridge.columns))
    if missing:
        raise ValueError(f"family bridge is missing columns: {', '.join(missing)}")

    projected: list[dict[str, str]] = []
    for family in families:
        if family["route_support"] != "a":
            continue
        source_ids = set(_tokens(family["route_a_family_ids"]))
        members = set(_tokens(family["member_interval_genes"]))
        if not source_ids or not members:
            continue
        for row in bridge.to_dict("records"):
            if str(row.get("route", "")) != "a":
                continue
            source_family = str(row.get("family_id", "")).strip()
            gene = str(row.get("gene_id", "")).strip()
            if source_family not in source_ids or gene not in members:
                continue
            projected.append({
                "gene_id": gene,
                "family_id": family["family_id"],
                "route": "a",
                "mapping_tier": str(row.get("mapping_tier", "")),
                "dependency_cluster_id": str(row.get("dependency_cluster_id", "")),
                "ambiguity_status": str(row.get("ambiguity_status", "")),
                "source_family_id": source_family,
                "projection_basis": family["family_match_basis"],
                "rank_contribution": "none",
            })
    frame = pd.DataFrame(projected, columns=CANONICAL_BRIDGE_COLUMNS)
    if len(frame):
        frame = frame.drop_duplicates().sort_values(
            ["family_id", "gene_id", "mapping_tier", "source_family_id"], kind="stable"
        ).reset_index(drop=True)
    return frame


def build_canonical_family_recovery(
    supported_union_arg: str,
    out_dir_arg: str,
    *,
    legacy_bridge_arg: str | None = None,
    case_id: str = "",
) -> dict[str, object]:
    """Project frozen supported A/B members into canonical family handoff rows."""
    union_path = runtime_path(supported_union_arg)
    if not union_path.exists():
        raise ValueError(f"supported member union does not exist: {union_path}")
    union = pd.read_csv(union_path, sep="\t", dtype=str).fillna("")
    missing = sorted(UNION_REQUIRED - set(union.columns))
    if missing:
        raise ValueError(f"supported member union is missing columns: {', '.join(missing)}")
    for column in ("route_c_supported", "route_c_family_ids", "route_c_orthogroup_ids"):
        if column not in union.columns:
            union[column] = ""
    if union["candidate_gene_id"].astype(str).str.strip().eq("").any():
        raise ValueError("supported member union contains an empty candidate_gene_id")
    if union["candidate_gene_id"].duplicated().any():
        raise ValueError("supported member union candidate_gene_id values must be unique")

    union_ids = set(union["candidate_gene_id"].astype(str).str.strip())
    a_map, b_map, b_og_only = _identity_maps(union)
    c_map, c_og_only = _route_c_identity_maps(union)
    unmatched_a = set(a_map)
    unmatched_b = set(b_map)
    locked_ambiguous_a: set[str] = set()
    locked_ambiguous_b: set[str] = set()
    families: list[dict[str, str]] = []

    # 1. Exact raw family labels.
    exact_pairs, _, _ = _unique_pairs(
        unmatched_a, unmatched_b, lambda a, b: _norm(a) == _norm(b)
    )
    for a, b in exact_pairs:
        families.append(_family_row(
            case_id=case_id,
            family_id=a,
            members=a_map[a] | b_map[b],
            route_support="a+b",
            match_basis="exact_family_label",
            route_a_ids=[a],
            route_b_ids=[b],
            route_b_orthogroups=[],
            ambiguous=False,
        ))
        unmatched_a.discard(a)
        unmatched_b.discard(b)

    # 2. Existing deterministic canonical family-key relation, but only when
    # mutually unique.  A stem shared by several paralog labels is not merged.
    deterministic_pairs, ambiguous_a, ambiguous_b = _unique_pairs(
        unmatched_a,
        unmatched_b,
        lambda a, b: bool(canonical_family_key(a))
        and canonical_family_key(a) == canonical_family_key(b),
    )
    for a, b in deterministic_pairs:
        families.append(_family_row(
            case_id=case_id,
            family_id=a,
            members=a_map[a] | b_map[b],
            route_support="a+b",
            match_basis="deterministic_family_key",
            route_a_ids=[a],
            route_b_ids=[b],
            route_b_orthogroups=[],
            ambiguous=False,
        ))
        unmatched_a.discard(a)
        unmatched_b.discard(b)
    locked_ambiguous_a |= (ambiguous_a & unmatched_a)
    locked_ambiguous_b |= (ambiguous_b & unmatched_b)

    # 3. Expanded Route-A descriptive label -> Route-B symbol relation.  Do not
    # use a looser tier to resolve an identity already ambiguous at tier 2.
    expanded_a = unmatched_a - locked_ambiguous_a
    expanded_b = unmatched_b - locked_ambiguous_b
    expanded_pairs, expanded_ambiguous_a, expanded_ambiguous_b = _unique_pairs(
        expanded_a, expanded_b, families_match_expanded
    )
    for a, b in expanded_pairs:
        families.append(_family_row(
            case_id=case_id,
            family_id=a,
            members=a_map[a] | b_map[b],
            route_support="a+b",
            match_basis="expanded_symbol_token",
            route_a_ids=[a],
            route_b_ids=[b],
            route_b_orthogroups=[],
            ambiguous=False,
        ))
        unmatched_a.discard(a)
        unmatched_b.discard(b)
    locked_ambiguous_a |= (expanded_ambiguous_a & unmatched_a)
    locked_ambiguous_b |= (expanded_ambiguous_b & unmatched_b)

    # 4. Route-specific family identities stay explicit; ambiguous cross-route
    # relations are never force-merged.
    for a in sorted(unmatched_a):
        families.append(_family_row(
            case_id=case_id,
            family_id=a,
            members=a_map[a],
            route_support="a",
            match_basis="identity_ambiguous_unmerged" if a in locked_ambiguous_a else "route_a_only",
            route_a_ids=[a],
            route_b_ids=[],
            route_b_orthogroups=[],
            ambiguous=a in locked_ambiguous_a,
        ))
    for b in sorted(unmatched_b):
        families.append(_family_row(
            case_id=case_id,
            family_id=b,
            members=b_map[b],
            route_support="b",
            match_basis="identity_ambiguous_unmerged" if b in locked_ambiguous_b else "route_b_only",
            route_a_ids=[],
            route_b_ids=[b],
            route_b_orthogroups=[],
            ambiguous=b in locked_ambiguous_b,
        ))
    for orthogroup in sorted(b_og_only):
        families.append(_family_row(
            case_id=case_id,
            family_id=f"orthogroup:{orthogroup}",
            members=b_og_only[orthogroup],
            route_support="b",
            match_basis="route_b_orthogroup_only",
            route_a_ids=[],
            route_b_ids=[],
            route_b_orthogroups=[orthogroup],
            ambiguous=False,
        ))

    # Attach Direct-C evidence only when family identity is mutually unique.
    # Ambiguous identities remain C-specific; no family agreement score exists.
    for c_family in sorted(c_map):
        candidates: list[int] = []
        for index, family in enumerate(families):
            labels = [family["family_id"], *_tokens(family["route_a_family_ids"]),
                      *_tokens(family["route_b_family_ids"])]
            if any(
                _norm(label) == _norm(c_family)
                or (
                    bool(canonical_family_key(label))
                    and canonical_family_key(label) == canonical_family_key(c_family)
                )
                or families_match_expanded(label, c_family)
                for label in labels if label
            ):
                candidates.append(index)
        if len(candidates) == 1:
            family = families[candidates[0]]
            family["member_interval_genes"] = ";".join(sorted(
                set(_tokens(family["member_interval_genes"])) | c_map[c_family]
            ))
            routes = set(_tokens(family["route_support"].replace("+", ";"))) | {"c"}
            family["route_support"] = "+".join(route for route in ("a", "b", "c") if route in routes)
            family["route_c_family_ids"] = ";".join(sorted(
                set(_tokens(family["route_c_family_ids"])) | {c_family}
            ))
            family["recovery_state"] = f"{family['route_support'].replace('+', '_')}_supported"
            family["family_match_basis"] += "+route_c_deterministic_identity"
        else:
            row = _family_row(
                case_id=case_id, family_id=c_family, members=c_map[c_family], route_support="c",
                match_basis="route_c_only" if not candidates else "route_c_identity_ambiguous_unmerged",
                route_a_ids=[], route_b_ids=[], route_b_orthogroups=[], ambiguous=bool(candidates),
            )
            row["route_c_family_ids"] = c_family
            families.append(row)
    for orthogroup in sorted(c_og_only):
        row = _family_row(
            case_id=case_id, family_id=f"orthogroup:{orthogroup}", members=c_og_only[orthogroup],
            route_support="c", match_basis="route_c_orthogroup_only", route_a_ids=[], route_b_ids=[],
            route_b_orthogroups=[], ambiguous=False,
        )
        row["route_c_orthogroup_ids"] = orthogroup
        families.append(row)

    for family in families:
        if "c" not in family["route_support"].split("+"):
            continue
        payload = "|".join([
            case_id, family["family_id"], family["route_support"],
            family["route_a_family_ids"], family["route_b_family_ids"],
            family["route_b_orthogroup_ids"], family["route_c_family_ids"],
            family["route_c_orthogroup_ids"], family["member_interval_genes"],
        ])
        family["hypothesis_id"] = "canonical-family-" + hashlib.sha256(
            payload.encode("utf-8")
        ).hexdigest()[:16]

    # Deterministic row order and a load-bearing membership-subset invariant.
    families.sort(key=lambda row: (row["family_id"], row["route_support"], row["hypothesis_id"]))
    family_member_ids = {
        gene
        for family in families
        for gene in _tokens(family["member_interval_genes"])
    }
    outside = sorted(family_member_ids - union_ids)
    if outside:
        raise RuntimeError(
            "canonical family recovery introduced genes outside supported membership: "
            + ", ".join(outside)
        )

    comparisons = [_comparison_row(family) for family in families]
    bridge_frame = _project_a_only_bridge(families, legacy_bridge_arg)

    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    family_path = out / "canonical_family_recovery.tsv"
    comparison_path = out / "canonical_family_comparison.tsv"
    bridge_path = out / "canonical_family_bridge.tsv"
    pd.DataFrame(families, columns=FAMILY_COLUMNS).to_csv(
        family_path, sep="\t", index=False, lineterminator="\n"
    )
    pd.DataFrame(comparisons, columns=COMPARISON_COLUMNS).to_csv(
        comparison_path, sep="\t", index=False, lineterminator="\n"
    )
    bridge_frame.to_csv(bridge_path, sep="\t", index=False, lineterminator="\n")

    basis_counts: dict[str, int] = {}
    for family in families:
        basis = family["family_match_basis"]
        basis_counts[basis] = basis_counts.get(basis, 0) + 1
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id,
        "producer_class": "canonical_supported_family_recovery",
        "membership_source": "supported_member_union",
        "supported_member_count": int(len(union)),
        "family_row_count": len(families),
        "family_member_union_count": len(family_member_ids),
        "family_members_subset_of_supported_union": True,
        "family_match_basis_counts": basis_counts,
        "canonical_bridge_policy": (
            "only original Route-A mapping tiers for A-only families are projected; "
            "Route-B semantic verification is never translated into historical bridge tiers, "
            "and weak Route-A tiers never gate an independently verified Route-B family"
        ),
        "legacy_triple_route_family_context_consumed": False,
        "rank_contribution": "none",
        "inputs": {
            "supported_member_union": file_fingerprint(union_path, full_hash=True),
            "legacy_family_bridge": (
                file_fingerprint(runtime_path(legacy_bridge_arg), full_hash=True)
                if legacy_bridge_arg else None
            ),
        },
        "outputs": {
            "family_recovery": str(family_path),
            "comparison": str(comparison_path),
            "family_bridge": str(bridge_path),
            "manifest": str(out / "canonical_family_recovery_manifest.json"),
        },
        "policy": [
            "supported_member_union.tsv is the sole candidate-membership authority",
            "family grouping may duplicate a gene across unresolved identity hypotheses but never adds a new gene",
            "exact/deterministic identity is attempted before bounded expanded matching",
            "deterministic or expanded ambiguity is preserved as separate route-specific families",
            "A/B family agreement is provenance only and never a numerical score bonus",
            "Route-B semantic support is not rewritten into the historical bridge-tier vocabulary",
        ],
    }
    manifest_path = out / "canonical_family_recovery_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        "family_recovery": str(family_path),
        "comparison": str(comparison_path),
        "family_bridge": str(bridge_path),
        "manifest": str(manifest_path),
        "family_count": len(families),
        "family_member_union_count": len(family_member_ids),
    }
