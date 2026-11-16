"""Route A/B discordance into review-required causal scenario rows (issue #363).

A/B discordance is an entrance to alternative causal scenarios, not a failure.
This adapter maps each row of a canonical comparison artifact
(``dual_route_hypothesis_comparison.tsv`` or
``triple_route_hypothesis_comparison.tsv``) deterministically onto review
candidates of the #362 causal-scenario ledger:

- A∩B     -> coding_gene, regulatory_element
- A-only  -> mapping_or_liftover_error, reference_parent_mismatch, sv_te_cnv_pav
- B-only  -> coding_gene, novel_or_unannotated_mechanism, multi_gene_or_epistasis
- neither -> phenotype_context_review, mapping_adequacy_review, environment_specific

Every four-state row stops at ``review_required`` with a judgment reason; the
adapter never excludes or promotes a candidate. When a route's provider was
missing (``route_a_evaluated=False`` / ``route_b_evaluated=False``) the state
is indeterminate, so the row collapses to a single ``not_evaluated`` row with
``scenario_type=provider_missing`` instead of being mixed into the ``neither``
bucket. Output rows can be handed to ``build_scenario_ledger`` unchanged (the
ledger fills ``rank_contribution=none`` and sorts/dedups).
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Mapping

import pandas as pd

from .causal_scenarios import _text

FOUR_STATE_LABELS = {
    "a_and_b_supported": "A∩B",
    "a_only_supported": "A-only",
    "b_only_supported": "B-only",
    "neither_supported": "neither",
}

# Fixed candidate order per state (issue #363 rules, verbatim).
STATE_SCENARIO_TYPES = {
    "A∩B": ("coding_gene", "regulatory_element"),
    "A-only": ("mapping_or_liftover_error", "reference_parent_mismatch", "sv_te_cnv_pav"),
    "B-only": ("coding_gene", "novel_or_unannotated_mechanism", "multi_gene_or_epistasis"),
    "neither": ("phenotype_context_review", "mapping_adequacy_review", "environment_specific"),
}

# The single not_evaluated row emitted for a provider gap names no causal form
# (the state is indeterminate), so it carries its own scenario type rather than
# being rounded into any review category.
PROVIDER_MISSING_SCENARIO_TYPE = "provider_missing"

RECOMMENDED_NEXT_EXPERIMENTS = {
    "coding_gene": "sequence candidate CDS and test coding variants",
    "regulatory_element": "promoter review and regulatory element assay",
    "sv_te_cnv_pav": "inspect long-read/SV/PAV evidence in the interval",
    "reference_parent_mismatch": "recheck marker pairs against both reference parents",
    "mapping_or_liftover_error": "recheck marker pairs and liftover chain",
    "multi_gene_or_epistasis": "test multiple interval genes and epistatic combinations",
    "environment_specific": "repeat in additional environments",
    "novel_or_unannotated_mechanism": "search for unannotated or novel mechanism candidates in the interval",
    "phenotype_context_review": "review phenotype measurement context and trait definition",
    "mapping_adequacy_review": "review mapping and marker adequacy before further causal analysis",
    PROVIDER_MISSING_SCENARIO_TYPE: "provide the missing route provider or record the route as intentionally absent",
}

PROVENANCE_COLUMNS = (
    "route_a_source_ids", "route_b_source_ids", "route_c_source_ids",
    "shared_source_ids", "source_dependence_status", "source_independence_status",
)


def classify_comparison_state(row: Mapping[str, object]) -> str:
    """Classify one canonical comparison row into A∩B / A-only / B-only / neither.

    Triple artifacts carry ``four_state`` (values from
    ``discovery_chain_benchmark.JOINT_STATES``; unknown values fail closed).
    Dual artifacts carry ``route_a_hypothesis_id`` / ``route_b_hypothesis_id``:
    both non-empty -> A∩B, only A -> A-only, only B -> B-only, both empty ->
    neither. A row with neither schema is not a canonical comparison artifact
    and raises ``ValueError``.
    """
    if "four_state" in row:
        value = _text(row.get("four_state")).casefold()
        if value not in FOUR_STATE_LABELS:
            raise ValueError(f"unsupported four_state: {value!r}")
        return FOUR_STATE_LABELS[value]
    if "route_a_hypothesis_id" in row or "route_b_hypothesis_id" in row:
        route_a = _text(row.get("route_a_hypothesis_id"), "")
        route_b = _text(row.get("route_b_hypothesis_id"), "")
        if route_a and route_b:
            return "A∩B"
        if route_a:
            return "A-only"
        if route_b:
            return "B-only"
        return "neither"
    raise ValueError(
        "comparison row is not a canonical comparison artifact "
        "(needs four_state or route_a/route_b_hypothesis_id columns)"
    )


def comparison_fingerprint(path: Path) -> str:
    """Fixed SHA-256 hex digest of a comparison file's canonical content.

    The digest covers the file's full content bytes with the data rows sorted
    lexicographically (the header row stays first, line endings normalized), so
    reordering the same comparison rows cannot churn the recorded fingerprint or
    the resulting ledger bytes. The case-run wiring computes this over the exact
    comparison artifact that produced the discordance rows and records it in
    each row's ``source_ids`` (``fingerprint=<hex>``, 64 hex chars). A missing
    file raises ``FileNotFoundError`` (no sentinel: callers wrap the wiring in
    their own guard).
    """
    lines = Path(path).read_bytes().splitlines()
    if len(lines) > 1:
        lines = lines[:1] + sorted(lines[1:])
    return hashlib.sha256(b"\n".join(lines)).hexdigest()


def _read_comparison(comparison: Path | str | pd.DataFrame) -> list[dict[str, object]]:
    if isinstance(comparison, pd.DataFrame):
        frame = comparison
    else:
        frame = pd.read_csv(comparison, sep="\t", dtype=str, keep_default_na=False)
    return frame.to_dict(orient="records")


def _basis(label: str, row: Mapping[str, object]) -> str:
    """Human-readable evidence basis for one state, filled from comparison columns."""
    pathway = _text(row.get("pathway_or_process"), "")
    family = _text(row.get("gene_family"), "")
    if label == "A∩B":
        base = f"route A and route B both support {pathway}" if pathway else "route A and route B both support the same hypothesis"
    elif label == "A-only":
        base = f"route A supports {pathway}; route B did not reach it" if pathway else "route A supports a hypothesis; route B did not reach it"
    elif label == "B-only":
        base = f"route B supports {pathway}; route A did not reach it" if pathway else "route B supports a hypothesis; route A did not reach it"
    else:
        base = f"no route supports {pathway}" if pathway else "no route supports the hypothesis"
    if family:
        return f"{base}; gene_family={family}"
    return base


def _dependence_caveats(row: Mapping[str, object]) -> list[str]:
    """Same-root / source-dependence caveats for a convergent (A∩B) row."""
    caveats: list[str] = []
    shared = _text(row.get("shared_source_ids"), "")
    if shared:
        caveats.append(
            f"same-root convergence: shared_source_ids={shared}; "
            "shared source ids are one independent source, not N"
        )
    for column in ("source_dependence_status", "source_independence_status"):
        status = _text(row.get(column), "")
        if status.casefold().startswith("shared_source"):
            caveats.append(f"{column}={status}")
    return caveats


def _source_ids(row: Mapping[str, object], comparison_id: str, fingerprint: str) -> str:
    parts = [f"comparison={comparison_id}", f"fingerprint={fingerprint}"]
    for column in PROVENANCE_COLUMNS:
        value = _text(row.get(column), "")
        if value:
            parts.append(f"{column}={value}")
    return ";".join(parts)


def _provider_gap_row(
    row: Mapping[str, object],
    *,
    case_id: str,
    fingerprint: str,
    comparison_id: str,
    missing_routes: tuple[str, ...],
) -> dict[str, str]:
    if len(missing_routes) == 2:
        reason = "route_a and route_b providers missing; discordance cannot be classified"
        next_experiment = "provide route A and route B providers or record either as intentionally absent"
    else:
        route = missing_routes[0]
        reason = f"route_{route} provider missing; discordance cannot be classified"
        next_experiment = (
            f"provide route {route.upper()} provider or record route {route.upper()} as intentionally absent"
        )
    return {
        "case_id": case_id,
        "scenario_type": PROVIDER_MISSING_SCENARIO_TYPE,
        "supporting_evidence": "",
        "contradicting_evidence": "",
        "assessment_state": "not_evaluated",
        "recommended_next_experiment": next_experiment,
        "reason": reason,
        "source_ids": _source_ids(row, comparison_id, fingerprint),
    }


def discordance_review_rows(
    comparison: Path | str | pd.DataFrame,
    *,
    case_id: str,
    input_fingerprint: str,
    route_a_evaluated: bool = True,
    route_b_evaluated: bool = True,
) -> list[dict[str, str]]:
    """Map comparison rows into review-required (or not_evaluated) scenario rows.

    ``comparison`` is a Path/str to a comparison TSV (read with string dtype
    and no NaN rewriting) or a DataFrame. ``input_fingerprint`` is computed by
    the caller (SHA-256 of the input, say); this function never computes it and
    an empty string works. When a route's provider was missing, every row of
    that comparison collapses to one ``not_evaluated`` row with reason
    ``route_X provider missing; discordance cannot be classified`` (an empty
    comparison yields ``[]``). Output is deterministic per input order; the
    ledger builder sorts and dedups.
    """
    rows = _read_comparison(comparison)
    if not rows:
        return []
    missing_routes = tuple(
        route for route, evaluated in (("a", route_a_evaluated), ("b", route_b_evaluated))
        if not evaluated
    )
    output: list[dict[str, str]] = []
    for index, row in enumerate(rows):
        label = classify_comparison_state(row)
        comparison_id = _text(row.get("comparison_id"), "") or f"row{index}"
        if missing_routes:
            output.append(_provider_gap_row(
                row, case_id=case_id, fingerprint=input_fingerprint,
                comparison_id=comparison_id, missing_routes=missing_routes,
            ))
            continue
        basis = _basis(label, row)
        reason = f"four-state={label}; judgment: {basis}"
        caveats = _dependence_caveats(row)
        if caveats:
            reason += "; " + "; ".join(caveats)
        for scenario_type in STATE_SCENARIO_TYPES[label]:
            output.append({
                "case_id": case_id,
                "scenario_type": scenario_type,
                "supporting_evidence": f"{label}: {basis}",
                "contradicting_evidence": "",
                "assessment_state": "review_required",
                "recommended_next_experiment": RECOMMENDED_NEXT_EXPERIMENTS[scenario_type],
                "reason": reason,
                "source_ids": _source_ids(row, comparison_id, input_fingerprint),
            })
    return output
