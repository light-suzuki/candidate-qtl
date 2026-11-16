"""Bounded adaptive synthesis for an over-budget Wet-candidate set (issue #426).

This stage exists only for the case where machine-first family -> member
expansion left more supported candidates than a Wet reviewer can carry.  It is
deliberately *not* a ranker:

    the LLM is asked which ambiguity to resolve, never which genes to keep.

The three questions it may answer are the ones a supervisor would ask a
colleague who knows the biology:

* which family / subfamily ambiguity, if resolved, would shrink the candidate
  set most safely;
* which additional material (a specific homolog, a paper, an expression set, an
  annotation rescue) would actually discriminate;
* which competing explanations must be kept rather than collapsed.

Everything the stage emits is a *proposal*.  Proposals never carry a numerical
gene score, never name a "final" candidate set, and never remove a candidate on
their own.  Narrowing happens only after a separate machine verification step
confirms the targeted evidence, which is why
:func:`apply_verified_refinement` accepts a verification table rather than the
proposal table.

The stage is gated on the budget decision itself: building a request against a
candidate set that is already within budget raises, so a caller cannot spend an
LLM call on a set that did not need one (issue #434's cheap-first contract).
"""

from __future__ import annotations

import hashlib
import json

import pandas as pd

from .core import runtime_path

SCHEMA_VERSION = "1.0"

# The budget statuses that legitimately reach this stage. Anything else means
# the caller skipped the gate or ignored its verdict.
REFINEMENT_ELIGIBLE_STATUSES = frozenset({"adaptive_refinement_required"})

PROPOSAL_COLUMNS = [
    "proposal_id",
    "case_id",
    "ambiguity_kind",
    "ambiguity_target",
    "affected_candidate_genes",
    "discriminating_evidence_request",
    "evidence_kind",
    "expected_effect",
    "competing_explanations_to_keep",
    "provider",
    "model",
    "prompt_fingerprint",
    "input_fingerprint",
    "output_fingerprint",
    "rank_contribution",
    "review_required",
]

AMBIGUITY_KINDS = frozenset({
    "family_vs_subfamily",
    "paralog_copy",
    "orthogroup_membership",
    "annotation_gap",
    "mechanism_alternative",
})

EVIDENCE_KINDS = frozenset({
    "specific_homolog",
    "literature",
    "expression",
    "annotation_rescue",
    "sequence_check",
    "synteny",
})

# A verification verdict may only shrink the candidate set when the machine
# actually falsified a candidate's support. Everything else keeps the candidate.
REFINEMENT_VERDICTS = frozenset({
    "support_falsified",
    "support_confirmed",
    "not_evaluated",
})
REFINEMENT_CONTRACT = "adaptive-refinement-1.1.0"
SELECTION_CONTRACT = "wet-review-inclusive-selection-1.0.0"

SYSTEM_PROMPT = """You are reviewing an over-budget QTL candidate set.

Do NOT choose which genes to keep. Do NOT rank genes. Do NOT assign any
numerical score, probability or confidence to a gene.

Answer only these questions:
1. Which family/subfamily/orthogroup ambiguity, if resolved, would most safely
   reduce the number of candidates?
2. What specific additional material would actually discriminate between the
   ambiguous members?
3. Which competing explanations must be preserved rather than collapsed?

If several paralogs are genuinely supported by the same evidence, say so and
propose keeping all of them. A smaller candidate set is not the goal; a
defensible one is.
"""


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _sha256_path(path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _candidate_id_column(frame: pd.DataFrame) -> str:
    for column in ("candidate_gene_id", "gene_id"):
        if column in frame.columns:
            return column
    raise ValueError("adaptive synthesis input requires candidate_gene_id or gene_id")


def _selection_column(frame: pd.DataFrame) -> str:
    """Return the explicit review decision column, never infer one from prose."""
    for column in ("decision", "recommendation", "selection", "disposition", "action"):
        if column in frame.columns:
            return column
    raise ValueError("review-inclusive selection requires an explicit decision column")


_SELECTION_DECISION_ALIASES = {
    "adopt": "adopt", "adopted": "adopt", "retain": "adopt", "keep": "adopt",
    "defer": "defer", "deferred": "defer",
    "review_required": "not_evaluated", "not_evaluated": "not_evaluated", "failed": "not_evaluated",
}


def _normalize_selection_decision(value: object, gene: str, source: str) -> str:
    decision = str(value).strip().lower()
    try:
        return _SELECTION_DECISION_ALIASES[decision]
    except KeyError as error:
        raise ValueError(f"unsupported {source} selection decision {decision!r} for {gene!r}") from error


def consume_review_inclusive_selection(
    wet_candidates_arg: str, review_arg: str, out_dir_arg: str, *,
    reviewer_id: str = "", reviewer_model: str = "", executor_id: str = "", executor_model: str = "",
    independent_review_arg: str | None = None, case_id: str = "", first_stage_sha256: str | None = None,
    original_assessment_scope_arg: str | None = None, original_assessment_scope_sha256: str | None = None,
) -> dict[str, object]:
    """Consume explicit second-stage adoption/deferral proposals mechanically.

    Rows are copied from the frozen Wet set. Only an explicit ``defer`` decision
    removes a row into the ledger; support flags and biological status are never
    rewritten, and missing/failed review decisions retain the original row.
    """
    source, review_path = runtime_path(wet_candidates_arg).resolve(), runtime_path(review_arg).resolve()
    if not source.is_file() or not review_path.is_file():
        raise ValueError("selection inputs must be existing files")
    original_assessment = bool(original_assessment_scope_arg)
    if original_assessment != bool(original_assessment_scope_sha256):
        raise ValueError("original assessment requires a sealed scope path and SHA")
    if original_assessment:
        from .boundary_tier_review import _load_scope, UPPER_EVIDENCE_SCOPED_POLICY
        scope, _, scope_plan = _load_scope(original_assessment_scope_arg, original_assessment_scope_sha256, case_id)
        if (scope.get("policy_version") != UPPER_EVIDENCE_SCOPED_POLICY
                or scope["inventory"]["sha256"] != first_stage_sha256
                or runtime_path(scope["artifacts"]["boundary"]["path"]).resolve() != source
                or scope["artifacts"]["boundary"]["sha256"] != _sha256_path(source)
                or not scope_plan["boundary"]):
            raise ValueError("original assessment requires the exact sealed v3 boundary")
        if reviewer_id or reviewer_model or independent_review_arg:
            raise ValueError("original assessment must not be represented as independent review")
    original = pd.read_csv(source, sep="\t", dtype=str).fillna("")
    id_column = _candidate_id_column(original)
    ids = original[id_column].astype(str).str.strip().tolist()
    if any(not value for value in ids) or len(ids) != len(set(ids)):
        raise ValueError("frozen Wet candidate IDs must be non-empty and unique")
    review = pd.read_csv(review_path, sep="\t", dtype=str).fillna("")
    review_id, decision_column = _candidate_id_column(review), _selection_column(review)
    review_ids = review[review_id].astype(str).str.strip().tolist()
    if any(not value for value in review_ids) or len(review_ids) != len(set(review_ids)):
        raise ValueError("selection review IDs must be non-empty and unique")
    unknown = sorted(set(review_ids) - set(ids))
    if unknown:
        raise ValueError("selection review references unknown Wet candidates: " + ", ".join(unknown))
    for required in ("executor_id", "executor_model", "reason"):
        if required not in review.columns or any(not str(value).strip() for value in review[required]):
            raise ValueError(f"selection proposal requires non-empty {required} on every row")
    if any(str(value).strip() != executor_id.strip() for value in review["executor_id"]):
        raise ValueError("selection proposal executor_id does not match the supplied executor")
    if any(str(value).strip() != executor_model.strip() for value in review["executor_model"]):
        raise ValueError("selection proposal executor_model does not match the supplied executor")
    if not original_assessment and (not reviewer_id.strip() or not reviewer_model.strip()):
        raise ValueError("an actual independent reviewer identity and model are required")
    if not executor_id.strip() or not executor_model.strip():
        raise ValueError("selection proposal executor identity and model are required")
    if not original_assessment and reviewer_id.strip() == executor_id.strip():
        raise ValueError("selection reviewer must differ from the selection proposal executor")
    source_hash_values = [str(value).strip() for value in review.get("source_wet_sha256", pd.Series(dtype=str))]
    if len(review) and (not source_hash_values or any(value != _sha256_path(source) for value in source_hash_values)):
        raise ValueError("selection review is bound to a different frozen Wet hash")
    if "source_wet_sha256" not in review.columns:
        raise ValueError("selection review must bind the frozen Wet hash")
    independent_agreement = True
    review_coverage = "complete" if set(review_ids) == set(ids) else "partial"
    proposal_records = {str(row[review_id]).strip(): row for row in review.to_dict("records")}
    proposal_decisions = {
        gene: _normalize_selection_decision(row[decision_column], gene, "selection proposal")
        for gene, row in proposal_records.items()
    }
    if not original_assessment and not independent_review_arg and "defer" in proposal_decisions.values():
        raise ValueError("deferral requires an independent review artifact and non-empty executor_id")
    independent = None
    independent_records = {}
    independent_decisions = {}
    if independent_review_arg:
        independent = pd.read_csv(runtime_path(independent_review_arg), sep="\t", dtype=str).fillna("")
        independent_id, independent_decision = _candidate_id_column(independent), _selection_column(independent)
        independent_ids = independent[independent_id].astype(str).str.strip().tolist()
        if any(not value for value in independent_ids) or len(independent_ids) != len(set(independent_ids)):
            raise ValueError("independent review IDs must be non-empty and unique")
        if set(independent_ids) != set(review_ids):
            raise ValueError("independent review roster must exactly match the proposal roster")
        independent_records = {str(row[independent_id]).strip(): row for row in independent.to_dict("records")}
        independent_decisions = {
            gene: _normalize_selection_decision(row[independent_decision], gene, "independent review")
            for gene, row in independent_records.items()
        }
        independent_wet_hashes = [str(value).strip() for value in independent.get("source_wet_sha256", pd.Series(dtype=str))]
        if not independent_wet_hashes or any(value != _sha256_path(source) for value in independent_wet_hashes):
            raise ValueError("independent review is bound to a different frozen Wet hash")
        proposal_hashes = [str(value).strip() for value in independent.get("proposal_sha256", pd.Series(dtype=str))]
        if not proposal_hashes or any(value != _sha256_path(review_path) for value in proposal_hashes):
            raise ValueError("independent review is not bound to the proposal artifact hash")
        artifact_reviewers = [str(value).strip() for value in independent.get("reviewer_id", pd.Series(dtype=str))]
        if not artifact_reviewers or any(value != reviewer_id.strip() for value in artifact_reviewers):
            raise ValueError("independent review artifact identity does not match reviewer_id")
        artifact_models = [str(value).strip() for value in independent.get("reviewer_model", pd.Series(dtype=str))]
        if not artifact_models or any(value != reviewer_model.strip() for value in artifact_models):
            raise ValueError("independent review artifact model does not match reviewer_model")
        independent_agreement = all(proposal_decisions[gene] == independent_decisions[gene] for gene in review_ids)
    decisions = {}
    review_decision_rows = []
    for gene in ids:
        proposal_row = proposal_records.get(gene, {})
        independent_row = independent_records.get(gene, {})
        proposal_decision = proposal_decisions.get(gene, "")
        independent_decision = independent_decisions.get(gene, "")
        if original_assessment:
            resolution_state = "original_assessment" if gene in proposal_records else "unassessed"
            applied_decision = proposal_decision if gene in proposal_records else "not_evaluated"
        elif not independent_review_arg:
            resolution_state = "not_evaluated_no_independent_review" if gene in proposal_records else "unreviewed"
            applied_decision = "not_evaluated"
        elif gene not in proposal_records:
            resolution_state = "unreviewed"
            applied_decision = "not_evaluated"
        elif proposal_decision == independent_decision:
            resolution_state = "agreement"
            applied_decision = proposal_decision
        else:
            resolution_state = "unresolved_disagreement"
            applied_decision = "not_evaluated"
            independent_agreement = False
        if gene in proposal_records:
            decisions[gene] = applied_decision
        review_decision_rows.append({
            "candidate_gene_id": gene,
            "proposal_decision": proposal_decision,
            "independent_decision": independent_decision,
            "proposal_reason": str(proposal_row.get("reason", proposal_row.get("rationale", ""))).strip(),
            "independent_reason": str(independent_row.get("reason", independent_row.get("rationale", ""))).strip(),
            "resolution_state": resolution_state,
            "applied_decision": applied_decision,
        })
    deferred = {gene for gene, decision in decisions.items() if decision in {"defer", "deferred"}}
    if deferred and not original_assessment and (not independent_review_arg or not executor_id.strip()):
        raise ValueError("deferral requires an independent review artifact and non-empty executor_id")
    if independent_review_arg:
        independent_reason = "reason" if "reason" in independent.columns else ("rationale" if "rationale" in independent.columns else None)
        if independent_reason is None or any(not str(row.get(independent_reason, "")).strip()
                                             for row in independent.to_dict("records")):
            raise ValueError("each independent review row requires a non-empty reason")
    kept = original[~original[id_column].astype(str).str.strip().isin(deferred)].copy()
    output = runtime_path(out_dir_arg).resolve()
    output.mkdir(parents=True, exist_ok=False)
    wet_path = output / "wet_candidates.tsv"
    kept.to_csv(wet_path, sep="\t", index=False, lineterminator="\n")
    review_decisions_path = output / "review_decisions.tsv"
    pd.DataFrame(review_decision_rows, columns=[
        "candidate_gene_id", "proposal_decision", "independent_decision",
        "proposal_reason", "independent_reason", "resolution_state", "applied_decision",
    ]).to_csv(review_decisions_path, sep="\t", index=False, lineterminator="\n")
    ledger_rows = []
    for row in review.to_dict("records"):
        gene = str(row[review_id]).strip()
        if gene in deferred:
            ledger_rows.append({"candidate_gene_id": gene, "disposition": "deferred",
                                "reason": str(row.get("reason", row.get("rationale", ""))).strip(),
                                "independent_reason": str(independent_records.get(gene, {}).get("reason", independent_records.get(gene, {}).get("rationale", ""))).strip(),
                                "reviewer_id": reviewer_id, "reviewer_model": reviewer_model,
                                "case_id": case_id, "stage": "second_stage_review",
                                "support_state_preserved": "yes",
                                **({"assessment_mode": "original_assessment_v1", "executor_id": executor_id, "executor_model": executor_model}
                                   if original_assessment else {})})
    ledger_path = output / "deferral_ledger.tsv"
    columns = ["candidate_gene_id", "disposition", "reason", "independent_reason", "reviewer_id", "reviewer_model", "case_id", "stage", "support_state_preserved"]
    if original_assessment:
        columns += ["assessment_mode", "executor_id", "executor_model"]
    pd.DataFrame(ledger_rows, columns=columns).to_csv(ledger_path, sep="\t", index=False, lineterminator="\n")
    disagreement_count = sum(row["resolution_state"] == "unresolved_disagreement" for row in review_decision_rows)
    manifest = {"contract": SELECTION_CONTRACT, "case_id": case_id, "stage": "review_inclusive_wet_selection",
                "source_wet_candidates": str(source), "source_wet_sha256": _sha256_path(source),
                "proposal_artifact": str(review_path),
                "independent_review_artifact": str(runtime_path(independent_review_arg).resolve()) if independent_review_arg else None,
                "independent_review_state": ("not_required_original_assessment" if original_assessment else
                                             ("agreement" if independent_agreement else "disagreement_retained") if independent_review_arg else "not_supplied"),
                **({"assessment_mode": "original_assessment_v1",
                    "executor_identity_basis": "declared_identity_model_consistency",
                    "original_assessment_scope": {
                    "path": str(runtime_path(original_assessment_scope_arg).resolve()), "sha256": original_assessment_scope_sha256}}
                   if original_assessment else {}),
                "input_review_sha256": _sha256_path(review_path),
                "first_stage_completion_sha256": first_stage_sha256,
                "input_candidate_count": len(ids), "reviewed_candidate_count": len(review_ids),
                "unreviewed_candidate_count": len(set(ids) - set(review_ids)),
                "review_coverage": review_coverage,
                "output_candidate_count": len(kept), "deferred_candidate_count": len(deferred),
                "disagreement_candidate_count": disagreement_count,
                "deferred_candidate_genes": sorted(deferred), "reviewer_id": reviewer_id, "reviewer_model": reviewer_model,
                "executor_id": executor_id, "executor_model": executor_model,
                "independent_review_sha256": _sha256_path(runtime_path(independent_review_arg).resolve()) if independent_review_arg else None,
                "review_decisions_artifact": str(review_decisions_path),
                "review_decisions_sha256": _sha256_path(review_decisions_path),
                "wet_candidates_sha256": _sha256_path(wet_path),
                "deferral_ledger_sha256": _sha256_path(ledger_path), "support_state": "preserved", "rank_contribution": "none"}
    manifest_path = output / "review_inclusive_selection_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"wet_candidates": str(wet_path), "deferral_ledger": str(ledger_path), "manifest": str(manifest_path),
            **({"assessment_mode": "original_assessment_v1"} if original_assessment else {}),
            "review_decisions": str(review_decisions_path),
            "output_candidate_count": len(kept), "deferred_candidate_count": len(deferred),
            "disagreement_candidate_count": disagreement_count,
            "review_coverage": review_coverage,
            "state": (("completed" if review_coverage == "complete" and all(value != "not_evaluated" for value in decisions.values())
                       else "completed_partial_assessment") if original_assessment else
                      "not_evaluated_no_independent_review" if not independent_review_arg else
                      "review_unresolved_disagreement" if not independent_agreement else
                      "completed_partial_review" if review_coverage == "partial" else "completed")}


apply_review_inclusive_selection = consume_review_inclusive_selection


def _read_budget_manifest(path_arg: str):
    path = runtime_path(path_arg)
    if not path.exists():
        raise ValueError(f"candidate budget manifest does not exist: {path}")
    return path, json.loads(path.read_text(encoding="utf-8"))


def build_adaptive_synthesis_request(
    candidates_arg: str,
    budget_manifest_arg: str,
    out_dir_arg: str,
    *,
    ambiguity_arg: str | None = None,
    case_id: str = "",
) -> dict[str, object]:
    """Emit the proposal request for an over-budget candidate set.

    Fail-closed on the budget verdict: a set that is within budget, empty, or
    already out of adaptive rounds must not reach an LLM. That check is the
    cheap-first stop rule (#434) expressed as code rather than as guidance.
    """
    manifest_path, budget = _read_budget_manifest(budget_manifest_arg)
    status = str(budget.get("status", "")).strip()
    if status not in REFINEMENT_ELIGIBLE_STATUSES:
        raise ValueError(
            "adaptive synthesis may only run when the candidate budget gate requested it; "
            f"gate status is {status!r} (expected one of {sorted(REFINEMENT_ELIGIBLE_STATUSES)})"
        )
    if not budget.get("adaptive_refinement_required", False):
        raise ValueError(
            "candidate budget manifest does not require adaptive refinement; refusing to "
            "spend an adaptive/LLM stage on a set that already fits the Wet budget"
        )

    source = runtime_path(candidates_arg)
    if not source.exists():
        raise ValueError(f"adaptive synthesis candidate input does not exist: {source}")
    candidates = pd.read_csv(source, sep="\t", dtype=str).fillna("")
    id_column = _candidate_id_column(candidates)
    candidate_ids = [str(value).strip() for value in candidates[id_column] if str(value).strip()]
    if not candidate_ids:
        raise ValueError("adaptive synthesis candidate input has no candidate ids")

    ambiguity_rows: list[dict[str, str]] = []
    ambiguity_path = runtime_path(ambiguity_arg) if ambiguity_arg else None
    ambiguity_state = "not_evaluated_no_ambiguity_input"
    if ambiguity_path is not None:
        if not ambiguity_path.exists():
            raise ValueError(f"ambiguity context does not exist: {ambiguity_path}")
        ambiguity = pd.read_csv(ambiguity_path, sep="\t", dtype=str).fillna("")
        unknown = sorted(
            {
                str(row.get("candidate_gene_id", "")).strip()
                for row in ambiguity.to_dict("records")
                if str(row.get("candidate_gene_id", "")).strip()
            }
            - set(candidate_ids)
        )
        if unknown:
            raise ValueError(
                "ambiguity context references genes outside the candidate set: " + ", ".join(unknown)
            )
        ambiguity_rows = ambiguity.to_dict("records")
        ambiguity_state = "supplied"

    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)

    request = {
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id,
        "stage": "adaptive_synthesis_request",
        "candidate_count": len(candidate_ids),
        "wet_budget": budget.get("wet_budget"),
        "adaptive_round": budget.get("adaptive_round"),
        "max_adaptive_rounds": budget.get("max_adaptive_rounds"),
        "candidate_gene_ids": candidate_ids,
        "ambiguity_context_state": ambiguity_state,
        "ambiguity_context": ambiguity_rows,
        "system_prompt": SYSTEM_PROMPT,
        "allowed_ambiguity_kinds": sorted(AMBIGUITY_KINDS),
        "allowed_evidence_kinds": sorted(EVIDENCE_KINDS),
        "policy": [
            "the reviewer is asked which ambiguity to resolve, never which genes to keep",
            "no numerical gene score, probability or confidence may be produced by this stage",
            "several supported paralogs may all be kept; a smaller set is not the objective",
            "proposals are not evidence until a separate machine verification confirms them",
        ],
    }
    request_path = out / "adaptive_synthesis_request.json"
    request_path.write_text(
        json.dumps(request, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id,
        "stage": "adaptive_synthesis_request",
        "candidate_count": len(candidate_ids),
        "candidate_input_sha256": _sha256_path(source),
        "budget_manifest_sha256": _sha256_path(manifest_path),
        "prompt_fingerprint": _sha256_text(SYSTEM_PROMPT),
        "output_sha256": _sha256_path(request_path),
        "ambiguity_context_state": ambiguity_state,
        "ambiguity_row_count": len(ambiguity_rows),
        "invocation_reason": "candidate_count_exceeds_wet_budget",
        "rank_contribution": "none",
    }
    manifest_path_out = out / "adaptive_synthesis_request_manifest.json"
    manifest_path_out.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        "request": str(request_path),
        "manifest": str(manifest_path_out),
        "candidate_count": len(candidate_ids),
        "invocation_reason": "candidate_count_exceeds_wet_budget",
    }


def import_adaptive_synthesis_proposals(
    proposals_arg: str | None,
    request_arg: str,
    out_dir_arg: str,
    *,
    provider: str = "",
    model: str = "",
    case_id: str = "",
) -> dict[str, object]:
    """Normalize reviewer proposals into the non-scoring proposal table.

    A missing provider is a state, not a failure: the table is written empty
    with ``provider_state=not_evaluated_no_provider`` so the surrounding
    workflow still completes and the over-budget set stops as review_required.
    """
    request_path = runtime_path(request_arg)
    if not request_path.exists():
        raise ValueError(f"adaptive synthesis request does not exist: {request_path}")
    request = json.loads(request_path.read_text(encoding="utf-8"))
    known_genes = {str(value).strip() for value in request.get("candidate_gene_ids", [])}

    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str]] = []
    provider_state = "not_evaluated_no_provider"
    rejected: list[str] = []

    if proposals_arg:
        source = runtime_path(proposals_arg)
        if not source.exists():
            raise ValueError(f"adaptive synthesis proposals do not exist: {source}")
        supplied = pd.read_csv(source, sep="\t", dtype=str).fillna("")
        provider_state = "imported"
        for index, record in enumerate(supplied.to_dict("records"), 1):
            kind = str(record.get("ambiguity_kind", "")).strip()
            evidence_kind = str(record.get("evidence_kind", "")).strip()
            if kind not in AMBIGUITY_KINDS:
                raise ValueError(
                    f"ambiguity_kind {kind!r} is not in the registered vocabulary {sorted(AMBIGUITY_KINDS)}"
                )
            if evidence_kind not in EVIDENCE_KINDS:
                raise ValueError(
                    f"evidence_kind {evidence_kind!r} is not in the registered vocabulary {sorted(EVIDENCE_KINDS)}"
                )
            # A proposal that carries a score is a ranking attempt wearing a
            # proposal's clothes. Reject rather than silently dropping the
            # column, so the contract breach is visible.
            for banned in ("score", "dry_score", "rank", "priority_score", "confidence_score"):
                if str(record.get(banned, "")).strip():
                    raise ValueError(
                        f"adaptive synthesis proposals must not carry {banned!r}; this stage never scores genes"
                    )
            affected = {
                token.strip()
                for token in str(record.get("affected_candidate_genes", "")).replace(",", ";").split(";")
                if token.strip()
            }
            unknown = sorted(affected - known_genes)
            if unknown:
                # Silently keeping a gene the request never contained would let
                # a proposal smuggle a candidate into the set.
                rejected.extend(unknown)
                raise ValueError(
                    "proposal references genes outside the requested candidate set: " + ", ".join(unknown)
                )
            rows.append({
                "proposal_id": str(record.get("proposal_id", "")).strip() or f"asyn-{index:03d}",
                "case_id": case_id,
                "ambiguity_kind": kind,
                "ambiguity_target": str(record.get("ambiguity_target", "")).strip(),
                "affected_candidate_genes": ";".join(sorted(affected)),
                "discriminating_evidence_request": str(record.get("discriminating_evidence_request", "")).strip(),
                "evidence_kind": evidence_kind,
                "expected_effect": str(record.get("expected_effect", "")).strip(),
                "competing_explanations_to_keep": str(record.get("competing_explanations_to_keep", "")).strip(),
                "provider": provider,
                "model": model,
                "prompt_fingerprint": str(request.get("system_prompt", "") and _sha256_text(request["system_prompt"])),
                "input_fingerprint": _sha256_path(request_path),
                "output_fingerprint": _sha256_path(source),
                "rank_contribution": "none",
                "review_required": "yes",
            })

    frame = pd.DataFrame(rows, columns=PROPOSAL_COLUMNS)
    table_path = out / "adaptive_synthesis_proposals.tsv"
    frame.to_csv(table_path, sep="\t", index=False, lineterminator="\n")

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id,
        "stage": "adaptive_synthesis_proposals",
        "provider_state": provider_state,
        "provider": provider,
        "model": model,
        "proposal_count": len(rows),
        "rejected_gene_count": len(rejected),
        "request_sha256": _sha256_path(request_path),
        "output_sha256": _sha256_path(table_path),
        "rank_contribution": "none",
        "policy": [
            "a proposal is a request for evidence, not evidence",
            "no proposal may narrow the candidate set before machine verification",
            "a missing provider is not_evaluated_no_provider, never an empty-but-clean result",
        ],
    }
    manifest_path = out / "adaptive_synthesis_proposals_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        "proposals": str(table_path),
        "manifest": str(manifest_path),
        "proposal_count": len(rows),
        "provider_state": provider_state,
    }


def apply_verified_refinement(
    candidates_arg: str,
    verification_arg: str | None,
    out_dir_arg: str,
    *,
    case_id: str = "",
) -> dict[str, object]:
    """Narrow the candidate set using machine-verified verdicts only.

    Only ``support_falsified`` removes a candidate, and only when the row also
    carries a non-empty ``verification_source``. ``not_evaluated`` and a missing
    verification table both leave the set untouched -- an unanswered question is
    not a reason to drop a gene.
    """
    source = runtime_path(candidates_arg).resolve()
    if not source.exists():
        raise ValueError(f"refinement candidate input does not exist: {source}")
    candidates = pd.read_csv(source, sep="\t", dtype=str).fillna("")
    id_column = _candidate_id_column(candidates)

    out = runtime_path(out_dir_arg).resolve()
    out.mkdir(parents=True, exist_ok=True)

    removed: list[str] = []
    verdict_counts = {verdict: 0 for verdict in sorted(REFINEMENT_VERDICTS)}
    verification_state = "not_evaluated_no_verification"
    verification_path = None

    if verification_arg:
        verification_path = runtime_path(verification_arg).resolve()
        if not verification_path.exists():
            raise ValueError(f"refinement verification does not exist: {verification_path}")
        verification = pd.read_csv(verification_path, sep="\t", dtype=str).fillna("")
        for required in ("candidate_gene_id", "verdict"):
            if required not in verification.columns:
                raise ValueError(f"refinement verification requires {required}")
        verification_state = "applied"
        if verification["candidate_gene_id"].astype(str).str.strip().duplicated().any():
            raise ValueError("refinement verification contains duplicate candidate_gene_id values")
        known = set(candidates[id_column].astype(str).str.strip())
        for record in verification.to_dict("records"):
            gene = str(record.get("candidate_gene_id", "")).strip()
            verdict = str(record.get("verdict", "")).strip()
            if verdict not in REFINEMENT_VERDICTS:
                raise ValueError(
                    f"verdict {verdict!r} is not in the registered vocabulary {sorted(REFINEMENT_VERDICTS)}"
                )
            if gene not in known:
                raise ValueError(f"refinement verification references unknown candidate {gene!r}")
            verdict_counts[verdict] += 1
            if verdict == "support_falsified":
                if not str(record.get("verification_source", "")).strip():
                    raise ValueError(
                        f"support_falsified for {gene!r} requires a verification_source; "
                        "an unsourced falsification may not remove a candidate"
                    )
                removed.append(gene)

    kept = candidates[~candidates[id_column].astype(str).str.strip().isin(set(removed))]
    refined_path = out / "refined_candidate_set.tsv"
    kept.to_csv(refined_path, sep="\t", index=False, lineterminator="\n")

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "contract": REFINEMENT_CONTRACT,
        "case_id": case_id,
        "stage": "adaptive_refinement",
        "verification_state": verification_state,
        "input_candidate_count": int(len(candidates)),
        "output_candidate_count": int(len(kept)),
        "removed_candidate_count": len(removed),
        "removed_candidate_genes": sorted(set(removed)),
        "verdict_counts": verdict_counts,
        "input_sha256": _sha256_path(source),
        "input_candidate_set": str(source),
        "verification_table": str(verification_path) if verification_path else None,
        "verification_sha256": _sha256_path(verification_path) if verification_path else None,
        "output_sha256": _sha256_path(refined_path),
        "refined_candidate_set": str(refined_path),
        "rank_contribution": "none",
        "policy": [
            "only a machine-verified, sourced support_falsified verdict removes a candidate",
            "not_evaluated leaves the candidate in the set; an unanswered question is not a negative",
            "supported paralogs are never collapsed to satisfy the budget",
        ],
    }
    manifest_path = out / "adaptive_refinement_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        "refined_candidate_set": str(refined_path),
        "manifest": str(manifest_path),
        "input_candidate_count": int(len(candidates)),
        "output_candidate_count": int(len(kept)),
        "removed_candidate_count": len(removed),
        "verification_state": verification_state,
    }


def verify_adaptive_refinement(
    refined_candidates_arg: str,
    refinement_manifest_arg: str,
    original_candidates_arg: str,
    *,
    case_id: str = "",
) -> dict[str, object]:
    """Verify that a refined set is an exact, sourced subset of its A/B union."""
    refined_path = runtime_path(refined_candidates_arg)
    manifest_path = runtime_path(refinement_manifest_arg)
    original_path = runtime_path(original_candidates_arg)
    for label, path in (
        ("refined candidate set", refined_path),
        ("adaptive refinement manifest", manifest_path),
        ("original supported candidate set", original_path),
    ):
        if not path.exists():
            raise ValueError(f"{label} does not exist: {path}")

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("contract") != REFINEMENT_CONTRACT:
        raise ValueError("adaptive refinement manifest contract is missing or stale")
    if manifest.get("stage") != "adaptive_refinement":
        raise ValueError("adaptive refinement manifest has the wrong stage")
    if manifest.get("verification_state") != "applied":
        raise ValueError("adaptive refinement was not backed by an applied verification table")
    if case_id and str(manifest.get("case_id", "")) != case_id:
        raise ValueError("adaptive refinement case_id does not match the Wet handoff case")
    if manifest.get("input_sha256") != _sha256_path(original_path):
        raise ValueError("adaptive refinement original candidate hash mismatch")
    if manifest.get("output_sha256") != _sha256_path(refined_path):
        raise ValueError("adaptive refinement output hash mismatch")

    verification_value = str(manifest.get("verification_table") or "").strip()
    if not verification_value:
        raise ValueError("adaptive refinement manifest does not identify its verification table")
    verification_path = runtime_path(verification_value)
    if not verification_path.exists():
        raise ValueError(f"adaptive refinement verification table does not exist: {verification_path}")
    if manifest.get("verification_sha256") != _sha256_path(verification_path):
        raise ValueError("adaptive refinement verification table hash mismatch")

    original = pd.read_csv(original_path, sep="\t", dtype=str).fillna("")
    refined = pd.read_csv(refined_path, sep="\t", dtype=str).fillna("")
    original_id_column = _candidate_id_column(original)
    refined_id_column = _candidate_id_column(refined)
    if original_id_column != refined_id_column or list(original.columns) != list(refined.columns):
        raise ValueError("adaptive refinement changed the candidate-set schema")
    for label, frame, id_column in (
        ("original", original, original_id_column),
        ("refined", refined, refined_id_column),
    ):
        ids = frame[id_column].astype(str).str.strip()
        if ids.eq("").any() or ids.duplicated().any():
            raise ValueError(f"adaptive refinement {label} candidate ids must be non-empty and unique")

    original_by_id = {
        str(record[original_id_column]).strip(): record
        for record in original.to_dict("records")
    }
    refined_by_id = {
        str(record[refined_id_column]).strip(): record
        for record in refined.to_dict("records")
    }
    injected = sorted(set(refined_by_id) - set(original_by_id))
    if injected:
        raise ValueError("adaptive refinement injected candidates outside the supported union: " + ", ".join(injected))
    mutated = sorted(
        gene for gene, record in refined_by_id.items()
        if record != original_by_id[gene]
    )
    if mutated:
        raise ValueError("adaptive refinement mutated supported candidate rows: " + ", ".join(mutated))

    verification = pd.read_csv(verification_path, sep="\t", dtype=str).fillna("")
    for required in ("candidate_gene_id", "verdict", "verification_source"):
        if required not in verification.columns:
            raise ValueError(f"adaptive refinement verification requires {required}")
    verification_ids = verification["candidate_gene_id"].astype(str).str.strip()
    if verification_ids.eq("").any() or verification_ids.duplicated().any():
        raise ValueError("adaptive refinement verification candidate ids must be non-empty and unique")
    allowed_removed: set[str] = set()
    for record in verification.to_dict("records"):
        gene = str(record.get("candidate_gene_id", "")).strip()
        verdict = str(record.get("verdict", "")).strip()
        if gene not in original_by_id:
            raise ValueError(f"adaptive refinement verification references unknown candidate {gene!r}")
        if verdict not in REFINEMENT_VERDICTS:
            raise ValueError(f"adaptive refinement verification has invalid verdict {verdict!r}")
        if verdict == "support_falsified":
            if not str(record.get("verification_source", "")).strip():
                raise ValueError(f"support_falsified for {gene!r} is not source-backed")
            allowed_removed.add(gene)

    actual_removed = set(original_by_id) - set(refined_by_id)
    if actual_removed != allowed_removed:
        raise ValueError("adaptive refinement removals do not match sourced support_falsified verdicts")
    if int(manifest.get("input_candidate_count", -1)) != len(original):
        raise ValueError("adaptive refinement input candidate count mismatch")
    if int(manifest.get("output_candidate_count", -1)) != len(refined):
        raise ValueError("adaptive refinement output candidate count mismatch")
    if set(manifest.get("removed_candidate_genes", [])) != actual_removed:
        raise ValueError("adaptive refinement manifest removal roster mismatch")

    return {
        "manifest": str(manifest_path),
        "original_candidate_set": str(original_path),
        "verification_table": str(verification_path),
        "input_candidate_count": len(original),
        "output_candidate_count": len(refined),
        "removed_candidate_count": len(actual_removed),
        "verification_state": "verified_lineage",
    }
