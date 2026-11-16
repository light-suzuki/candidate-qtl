from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from .core import runtime_path
from .schemas import HYPOTHESIS_STATUSES, PROPOSAL_ORIGINS


POLICY_VERSION = "next-experiment-v2"
EXPERIMENTS = (
    "additional_grid_markers",
    "top_n_parental_sanger",
    "multiplex_amplicon_sequencing",
    "target_capture",
    "sv_te_cnv_pav_branch",
    "expression_assay",
    "annotation_rescue",
    "non_gene_follow_up",
)
ASSAY_STATES = {"not_started", "partial", "failed", "complete"}
ASSAY_OUTCOMES = {"technical_failure", "informative_negative", "positive", "inconclusive", "unknown"}
INFORMATION_ORDER = {"high": 0, "medium": 1, "low": 2, "unknown": 3}
COST_ORDER = {"low": 0, "medium": 1, "high": 2, "unknown": 3}
SAMPLE_REQUIREMENTS = {
    "additional_grid_markers": ("informative_recombinants", "segregating_population"),
    "top_n_parental_sanger": ("parents",),
    "multiplex_amplicon_sequencing": ("parents",),
    "target_capture": ("parents",),
    "sv_te_cnv_pav_branch": ("parents",),
    "expression_assay": ("target_tissue",),
    "annotation_rescue": ("reference_proteome", "comparative_annotation_resources"),
    "non_gene_follow_up": ("association_cohort", "population_phenotypes"),
}
HYPOTHESES = {
    "additional_grid_markers": "QTL interval and recombinant breakpoint refinement",
    "top_n_parental_sanger": "coding and targeted candidate-region parental variants",
    "multiplex_amplicon_sequencing": "parallel targeted variation across multiple candidates",
    "target_capture": "broad coding and selected regulatory variation across the locus",
    "sv_te_cnv_pav_branch": "structural variation, TE, CNV, PAV, and reference-missing alleles",
    "expression_assay": "tissue- and stage-specific expression differences",
    "annotation_rescue": "comparative annotation projection for reference-missing loci",
    "non_gene_follow_up": "linked-QTL, association/haplotype, and non-gene locus mechanisms",
}


def _load_structured(path: Path | None) -> Any:
    if path is None:
        return None
    if path.suffix.casefold() in {".yaml", ".yml"}:
        return yaml.safe_load(path.read_text(encoding="utf-8"))
    if path.suffix.casefold() == ".json":
        return json.loads(path.read_text(encoding="utf-8"))
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _text(value: Any) -> str:
    if isinstance(value, pd.DataFrame):
        return value.to_csv(sep="\t", index=False).casefold()
    return json.dumps(value, ensure_ascii=False, sort_keys=True, default=str).casefold() if value is not None else ""


def _known_bool(value: Any) -> str:
    normalized = str(value).strip().casefold()
    if normalized in {"true", "yes", "1", "available", "complete"}:
        return "yes"
    if normalized in {"false", "no", "0", "unavailable", "missing"}:
        return "no"
    return "unknown"


def _cost_status(cost: str, budget: str) -> str:
    if cost not in COST_ORDER or budget not in COST_ORDER or cost == "unknown" or budget == "unknown":
        return "unknown"
    return "allowed" if COST_ORDER[cost] <= COST_ORDER[budget] else "blocked"


def _assay_history(frame: pd.DataFrame | None) -> dict[str, tuple[str, str, str]]:
    if frame is None or frame.empty:
        return {}
    experiment_column = next((name for name in ("experiment_id", "experiment_class", "experiment") if name in frame.columns), None)
    if not experiment_column or "state" not in frame.columns or "outcome" not in frame.columns:
        raise ValueError("assay outcomes require experiment_id/experiment_class, state, and outcome")
    history: dict[str, tuple[str, str, str]] = {}
    for row in frame.itertuples(index=False):
        experiment = str(getattr(row, experiment_column)).strip()
        state = str(getattr(row, "state")).strip().casefold()
        outcome = str(getattr(row, "outcome")).strip().casefold()
        if state not in ASSAY_STATES:
            raise ValueError(f"unsupported assay state: {state}")
        if outcome not in ASSAY_OUTCOMES:
            raise ValueError(f"unsupported assay outcome: {outcome}")
        reason = str(getattr(row, "failure_reason", "")).strip()
        if experiment in EXPERIMENTS:
            history[experiment] = (state, outcome, reason)
    return history


def _elimination_from_recombinants(frame: pd.DataFrame | None) -> tuple[str, str]:
    if frame is None or frame.empty:
        return "unknown", "no recombinant candidate-count calculation supplied"
    for before_name, after_name in (("candidates_before", "candidates_after"), ("candidate_count_before", "candidate_count_after")):
        if before_name in frame.columns and after_name in frame.columns:
            values = []
            for before, after in zip(frame[before_name], frame[after_name], strict=False):
                try:
                    values.append(max(0, int(float(before)) - int(float(after))))
                except (TypeError, ValueError):
                    continue
            if values:
                return str(max(values)), f"observed {before_name}-{after_name} from informative recombinant table"
    return "unknown", "recombinant table lacks observed before/after candidate counts"


def _mapping_gate_text(mapping_gate: Any) -> str:
    if isinstance(mapping_gate, pd.DataFrame):
        return " ".join(mapping_gate.astype(str).to_numpy().ravel()).casefold()
    return _text(mapping_gate)


def _eligible_hypotheses(artifacts: dict[str, Any]) -> pd.DataFrame:
    """Reviewed non-LLM locus hypotheses; LLM proposals stay pending_review."""
    hypotheses = artifacts.get("locus_hypotheses")
    if not isinstance(hypotheses, pd.DataFrame) or "hypothesis_type" not in hypotheses.columns:
        return pd.DataFrame()
    eligible = hypotheses.copy()
    if "proposal_origin" in eligible.columns:
        eligible = eligible[~eligible["proposal_origin"].str.casefold().eq("llm")]
    if "review_status" in eligible.columns:
        eligible = eligible[eligible["review_status"].str.casefold().eq("reviewed")]
    if "status" in eligible.columns:
        eligible = eligible[~eligible["status"].str.casefold().isin({"rejected", "superseded"})]
    return eligible


def _hypotheses_requesting(eligible: pd.DataFrame, experiment: str) -> int:
    if eligible.empty or "required_experiment_ids" not in eligible.columns:
        return 0
    return int(eligible["required_experiment_ids"].map(
        lambda value: experiment in {token.strip() for token in str(value).replace("|", ";").split(";")}
    ).sum())


NON_GENE_HYPOTHESIS_TYPES = frozenset({"unknown_locus_mechanism", "linked_qtl", "multi_gene_combination"})
NON_GENE_MAPPING_TOKENS = ("linked_qtl_suspected", "low_recombination_region")


def _has_structural_signal(artifacts: dict[str, Any]) -> bool:
    decision_text = _text(artifacts.get("candidate_path_decision"))
    if any(token in decision_text for token in ("broaden_to_regulatory_sv", "structural_variant", "sv_te_cnv_pav")):
        return True
    flags = artifacts.get("candidate_flags")
    if isinstance(flags, pd.DataFrame):
        for column in flags.columns:
            if not any(token in column.casefold() for token in ("tandem", "repeat", "paralog", "structural")):
                continue
            values = flags[column].astype(str).str.strip().str.casefold()
            if values.ne("").any() and (~values.isin({"0", "false", "no", "none", "na", "unknown"})).any():
                return True
        note_columns = [column for column in flags.columns if column.casefold() in {"notes", "reason", "failure_reason"}]
        if note_columns and any(token in " ".join(flags[note_columns].astype(str).to_numpy().ravel()).casefold() for token in ("tandem", "repeat", "paralog", "non-ampl")):
            return True
    coverage = artifacts.get("exon_coverage")
    if isinstance(coverage, pd.DataFrame):
        values = " ".join(coverage.astype(str).to_numpy().ravel()).casefold()
        if any(token in values for token in ("coverage_gap", "partial", "failed", "non-ampl")):
            return True
    outcomes = artifacts.get("assay_outcomes")
    if isinstance(outcomes, pd.DataFrame) and "failure_reason" in outcomes.columns:
        reasons = " ".join(outcomes["failure_reason"].astype(str)).casefold()
        if any(token in reasons for token in ("non-ampl", "repeat", "paralog", "structural")):
            return True
    wet_state = artifacts.get("wet_candidate_state")
    if isinstance(wet_state, pd.DataFrame):
        values = " ".join(wet_state.astype(str).to_numpy().ravel()).casefold()
        if any(token in values for token in ("sv_pav", "regulatory_sv", "primer_mismatch_warning")):
            return True
    eligible = _eligible_hypotheses(artifacts)
    if not eligible.empty and eligible["hypothesis_type"].str.casefold().isin({"sv_cnv_pav", "te_regulatory"}).any():
        return True
    return False


def _has_annotation_gap(artifacts: dict[str, Any]) -> bool:
    """A reviewed annotation gap favors comparative rescue before coding sequencing."""
    eligible = _eligible_hypotheses(artifacts)
    if not eligible.empty and eligible["hypothesis_type"].str.casefold().eq("annotation_missing").any():
        return True
    return _hypotheses_requesting(eligible, "annotation_rescue") > 0


def _has_non_gene_signal(artifacts: dict[str, Any]) -> bool:
    """Linked-QTL / low-recombination / non-gene locus evidence favors follow-up."""
    if any(token in _mapping_gate_text(artifacts.get("mapping_gate")) for token in NON_GENE_MAPPING_TOKENS):
        return True
    eligible = _eligible_hypotheses(artifacts)
    if not eligible.empty and eligible["hypothesis_type"].str.casefold().isin(NON_GENE_HYPOTHESIS_TYPES).any():
        return True
    return False


def _needs_mapping_refinement(mapping_gate: Any, recombinants: pd.DataFrame | None) -> bool:
    if recombinants is not None and len(recombinants):
        return True
    mapping_text = _mapping_gate_text(mapping_gate)
    # localization-insufficiency tokens route to interval refinement. Tokens are
    # matched deliberately: "fine_mapping_required" (broad interval) and
    # "linked_qtl_suspected" (multiple peaks need discriminating recombinants)
    # both imply the interval is not yet localized. "low_recombination_region"
    # is NOT matched: adding markers there cannot resolve the interval, so it
    # stays a review-only signal for the non_gene_follow_up branch.
    return any(token in mapping_text for token in (
        "insufficient", "additional_grid", "low_marker", "marker_density_low",
        "marker_density_insufficient", "recombinant_required",
        "fine_mapping_required", "linked_qtl_suspected",
    ))


def select_next_experiment(
    ranked_arg: str,
    case_manifest_arg: str,
    out_dir_arg: str,
    *,
    mapping_gate_arg: str | None = None,
    candidate_path_decision_arg: str | None = None,
    candidate_flags_arg: str | None = None,
    assay_outcomes_arg: str | None = None,
    exon_coverage_arg: str | None = None,
    informative_recombinants_arg: str | None = None,
    wet_candidate_state_arg: str | None = None,
    hypotheses_arg: str | None = None,
    remaining_candidates_arg: str | None = None,
    narrowing_arg: str | None = None,
) -> dict[str, object]:
    """Rank experimental actions without modifying biological candidate rank."""
    ranked_path = runtime_path(ranked_arg)
    manifest_path = runtime_path(case_manifest_arg)
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    ranked_bytes = ranked_path.read_bytes()
    ranked = pd.read_csv(ranked_path, sep="\t", dtype=str).fillna("")
    if "gene_id" not in ranked.columns:
        raise ValueError("ranked candidates require gene_id")
    manifest = _load_structured(manifest_path) or {}
    if not isinstance(manifest, dict):
        raise ValueError("case manifest must be an object")
    support = manifest.get("decision_support") or {}
    if not isinstance(support, dict):
        raise ValueError("decision_support must be an object")

    optional_paths = {
        "mapping_gate": runtime_path(mapping_gate_arg) if mapping_gate_arg else None,
        "candidate_path_decision": runtime_path(candidate_path_decision_arg) if candidate_path_decision_arg else None,
        "candidate_flags": runtime_path(candidate_flags_arg) if candidate_flags_arg else None,
        "assay_outcomes": runtime_path(assay_outcomes_arg) if assay_outcomes_arg else None,
        "exon_coverage": runtime_path(exon_coverage_arg) if exon_coverage_arg else None,
        "informative_recombinants": runtime_path(informative_recombinants_arg) if informative_recombinants_arg else None,
        "wet_candidate_state": runtime_path(wet_candidate_state_arg) if wet_candidate_state_arg else None,
        "locus_hypotheses": runtime_path(hypotheses_arg) if hypotheses_arg else None,
        "remaining_candidates": runtime_path(remaining_candidates_arg) if remaining_candidates_arg else None,
        "integrated_narrowing": runtime_path(narrowing_arg) if narrowing_arg else None,
    }
    for name, path in optional_paths.items():
        if path is not None and not path.exists():
            raise ValueError(f"{name} does not exist: {path}")
    artifacts = {name: _load_structured(path) for name, path in optional_paths.items()}
    hypotheses_input = artifacts["locus_hypotheses"]
    if isinstance(hypotheses_input, pd.DataFrame) and len(hypotheses_input):
        required_hypothesis_columns = {"hypothesis_type", "proposal_origin", "review_status", "status"}
        missing_hypothesis_columns = sorted(required_hypothesis_columns - set(hypotheses_input.columns))
        if missing_hypothesis_columns:
            raise ValueError(
                "locus hypotheses require fail-closed origin/review/status columns: "
                + ", ".join(missing_hypothesis_columns)
            )
        unknown_origins = sorted(set(hypotheses_input["proposal_origin"].str.casefold()) - PROPOSAL_ORIGINS)
        unknown_statuses = sorted(set(hypotheses_input["status"].str.casefold()) - HYPOTHESIS_STATUSES)
        if unknown_origins or unknown_statuses:
            raise ValueError(f"invalid hypothesis origin/status: origins={unknown_origins}; statuses={unknown_statuses}")
        llm_rows = hypotheses_input["proposal_origin"].str.casefold().eq("llm")
        if (
            ~hypotheses_input.loc[llm_rows, "review_status"].str.casefold().isin({"pending_review", "unreviewed"})
        ).any() or (~hypotheses_input.loc[llm_rows, "status"].str.casefold().eq("pending_review")).any():
            raise ValueError("LLM hypothesis action context must remain pending_review")
    history = _assay_history(artifacts["assay_outcomes"] if isinstance(artifacts["assay_outcomes"], pd.DataFrame) else None)
    wet_frame = artifacts["wet_candidate_state"] if isinstance(artifacts["wet_candidate_state"], pd.DataFrame) else None
    wet_counts = wet_frame["candidate_state"].value_counts().to_dict() if wet_frame is not None and "candidate_state" in wet_frame.columns else {}
    hypothesis_frame = artifacts["locus_hypotheses"] if isinstance(artifacts["locus_hypotheses"], pd.DataFrame) else None

    allowed_raw = support.get("allowed_technologies")
    allowed = {str(value) for value in allowed_raw} if isinstance(allowed_raw, list) else None
    samples = support.get("available_samples") if isinstance(support.get("available_samples"), dict) else {}
    budget = str(support.get("budget_class", "unknown")).casefold()
    if budget not in COST_ORDER:
        raise ValueError("budget_class must be low, medium, high, or unknown")
    local_costs = support.get("local_cost_classes") if isinstance(support.get("local_cost_classes"), dict) else {}
    amplicon_counts = support.get("estimated_amplicon_counts") if isinstance(support.get("estimated_amplicon_counts"), dict) else {}

    recombinants = artifacts["informative_recombinants"] if isinstance(artifacts["informative_recombinants"], pd.DataFrame) else None
    mapping_first = _needs_mapping_refinement(artifacts["mapping_gate"], recombinants)
    sv_first = _has_structural_signal(artifacts)
    eligible_hypotheses = _eligible_hypotheses(artifacts)
    annotation_count = int(eligible_hypotheses["hypothesis_type"].str.casefold().eq("annotation_missing").sum()) if not eligible_hypotheses.empty else 0
    non_gene_count = int(eligible_hypotheses["hypothesis_type"].str.casefold().isin(NON_GENE_HYPOTHESIS_TYPES).sum()) if not eligible_hypotheses.empty else 0
    annotation_first = _has_annotation_gap(artifacts)
    non_gene_first = _has_non_gene_signal(artifacts)
    elimination_count, elimination_basis = _elimination_from_recombinants(recombinants)
    # Issue #285: an optional post-elimination remaining set switches the
    # candidate counts to the survivors. It is consumed as-is: surviving genes
    # are never declared biologically true, absent ones never declared false.
    remaining_frame = artifacts["remaining_candidates"] if isinstance(artifacts["remaining_candidates"], pd.DataFrame) else None
    remaining_ids: set[str] | None = None
    if remaining_frame is not None:
        if "candidate_gene_id" not in remaining_frame.columns:
            raise ValueError("remaining candidates require candidate_gene_id")
        remaining_ids = {str(value).strip() for value in remaining_frame["candidate_gene_id"]}
        ranked_ids = {str(value).strip() for value in ranked["gene_id"]}
        unknown_ids = sorted(remaining_ids - ranked_ids)
        if unknown_ids:
            raise ValueError(
                "remaining candidates contain gene ids absent from the ranked table (fail-closed): "
                + ", ".join(unknown_ids[:5]) + ("..." if len(unknown_ids) > 5 else "")
            )
    if remaining_ids is not None:
        ranked_subset = ranked[ranked["gene_id"].astype(str).str.strip().isin(remaining_ids)]
        top_count = int(ranked_subset["dry_group"].isin(["top", "reserve"]).sum()) if "dry_group" in ranked_subset.columns else len(ranked_subset)
    else:
        top_count = int(ranked["dry_group"].isin(["top", "reserve"]).sum()) if "dry_group" in ranked.columns else len(ranked)

    narrowing_frame = artifacts["integrated_narrowing"] if isinstance(artifacts["integrated_narrowing"], pd.DataFrame) else None
    narrowing_counts: dict[str, int] | None = None
    if narrowing_frame is not None:
        required_narrowing = {"gene_id", "narrowing_tier", "narrowing_priority"}
        missing_narrowing = sorted(required_narrowing - set(narrowing_frame.columns))
        if missing_narrowing:
            raise ValueError("integrated narrowing requires: " + ", ".join(missing_narrowing))
        narrowed_ids = {str(value).strip() for value in narrowing_frame["gene_id"]}
        ranked_ids = {str(value).strip() for value in ranked["gene_id"]}
        unknown_ids = sorted(narrowed_ids - ranked_ids)
        if unknown_ids:
            raise ValueError("integrated narrowing contains genes absent from ranked candidates: " + ", ".join(unknown_ids[:5]))
        eligible_narrowing = narrowing_frame
        if remaining_ids is not None:
            eligible_narrowing = narrowing_frame[narrowing_frame["gene_id"].astype(str).str.strip().isin(remaining_ids)]
        narrowing_counts = {
            str(key): int(value)
            for key, value in eligible_narrowing["narrowing_tier"].value_counts().to_dict().items()
        }
        # N1/N2 are the existing Wet-priority tiers. Consuming their count here
        # changes no biological rank; it only sizes the next Wet action.
        top_count = int(eligible_narrowing["narrowing_tier"].isin(["N1", "N2"]).sum())

    option_rows: list[dict[str, object]] = []
    dependency_rows: list[dict[str, str]] = []
    review_rows: list[dict[str, str]] = []
    for experiment in EXPERIMENTS:
        if remaining_ids is not None:
            facts: list[str] = [f"candidate_count={len(remaining_ids)}",
                                f"candidate_count_before_elimination={len(ranked)}",
                                f"candidate_count_after_elimination={len(remaining_ids)}",
                                f"top_or_reserve_count={top_count}"]
        else:
            facts: list[str] = [f"candidate_count={len(ranked)}", f"top_or_reserve_count={top_count}"]
        if narrowing_counts is not None:
            facts.extend([
                "integrated_narrowing_status=evaluated",
                f"narrowing_wet_priority_count={top_count}",
                "narrowing_tier_counts=" + ",".join(
                    f"{key}:{narrowing_counts.get(key, 0)}" for key in ("N1", "N2", "N3", "N4")
                ),
            ])
        rules: list[str] = []
        blocking: list[str] = []
        review: list[str] = []
        technology_status = "unknown" if allowed is None else ("allowed" if experiment in allowed else "blocked")
        dependency_rows.append({"experiment_id": experiment, "dependency": "technology_allowed", "status": technology_status, "source": "decision_support.allowed_technologies"})
        if technology_status == "blocked":
            blocking.append("technology_not_allowed")
        elif technology_status == "unknown":
            review.append("allowed_technologies_missing")

        sample_states = []
        for requirement in SAMPLE_REQUIREMENTS[experiment]:
            state = _known_bool(samples.get(requirement, "unknown"))
            sample_states.append(f"{requirement}={state}")
            dependency_rows.append({"experiment_id": experiment, "dependency": f"sample:{requirement}", "status": state, "source": "decision_support.available_samples"})
            if state == "no":
                blocking.append(f"sample_unavailable:{requirement}")
            elif state == "unknown":
                review.append(f"sample_unknown:{requirement}")

        cost = str(local_costs.get(experiment, "unknown")).casefold()
        if cost not in COST_ORDER:
            raise ValueError(f"local cost class for {experiment} must be low, medium, high, or unknown")
        cost_source = "case_manifest:decision_support.local_cost_classes" if experiment in local_costs else "not_supplied"
        cost_status = _cost_status(cost, budget)
        dependency_rows.append({"experiment_id": experiment, "dependency": "budget", "status": cost_status, "source": f"budget={budget};cost={cost}"})
        if cost_status == "blocked":
            blocking.append("cost_exceeds_budget_class")
        elif cost_status == "unknown":
            review.append("cost_or_budget_unknown")

        state, outcome, failure_reason = history.get(experiment, ("not_started", "unknown", ""))
        facts.extend((f"assay_state={state}", f"assay_outcome={outcome}"))
        if hypothesis_frame is not None and "required_experiment_ids" in hypothesis_frame.columns:
            linked_count = hypothesis_frame["required_experiment_ids"].map(
                lambda value: experiment in {token.strip() for token in str(value).replace("|", ";").split(";")}
            ).sum()
            if linked_count:
                facts.append(f"linked_hypotheses={int(linked_count)}")
        if state == "complete":
            rules.append("R_HISTORY_COMPLETE")
        if outcome == "technical_failure":
            rules.append("R_TECHNICAL_FAILURE_IS_NOT_NEGATIVE")
            blocking.append("same_technology_technical_failure_review")
        elif outcome == "informative_negative":
            rules.append("R_INFORMATIVE_NEGATIVE_ONLY_ADDRESSES_MEASURED_SCOPE")
        elif state == "partial":
            rules.append("R_PARTIAL_COVERAGE_REQUIRES_COMPLETION_OR_ALTERNATIVE")
            review.append("partial_assay_scope")
        if state == "failed" and outcome != "technical_failure":
            rules.append("R_FAILED_ASSAY_REQUIRES_OUTCOME_CLASSIFICATION")
            review.append("failed_assay_outcome_requires_classification")

        info = "medium"
        priority = 30
        risk = "medium"
        runtime = "medium"
        reuse = "medium"
        reason = "compare this action after prerequisites and local constraints are reviewed"
        if experiment == "additional_grid_markers":
            info = "high" if mapping_first else "medium"
            priority = 1 if mapping_first else 20
            risk, runtime, reuse = "low", "low", "high"
            if mapping_first:
                rules.append("R_MAPPING_INSUFFICIENT_OR_RECOMBINANTS_FAVOR_INTERVAL_REFINEMENT")
                reason = "mapping evidence indicates that breakpoint refinement may reduce the candidate interval before sequencing"
        elif experiment == "sv_te_cnv_pav_branch":
            info = "high" if sv_first else "low"
            priority = 0 if sv_first else 40
            risk, runtime, reuse = "medium", "high", "high"
            if sv_first:
                rules.append("R_STRUCTURAL_OR_AMPLIFICATION_SIGNAL_FAVORS_SV_BRANCH")
                reason = "repeat/paralog/coverage or amplification signals make a structural-variant route more informative than repeated Sanger"
        elif experiment == "top_n_parental_sanger":
            info = "medium" if top_count <= 20 else "low"
            priority = 10
            risk, runtime, reuse = "medium", "medium", "medium"
            reason = "targeted parental sequencing tests nominated regions but does not exclude unmeasured regulatory or structural causes"
            incomplete_wet = sum(int(wet_counts.get(state_name, 0)) for state_name in ("not_tested", "partially_covered", "no_parental_difference", "inconclusive"))
            coding_cleared = int(wet_counts.get("coding_region_cleared", 0))
            if incomplete_wet:
                rules.append("R_WET_STATE_INCOMPLETE_SCOPE_REQUIRES_COMPLETION")
                review.append("partial_wet_candidate_scope")
                facts.append(f"wet_candidates_incomplete={incomplete_wet}")
            if coding_cleared:
                rules.append("R_WET_CODING_NEGATIVE_IS_SCOPE_LIMITED")
                facts.append(f"wet_candidates_coding_cleared={coding_cleared}")
        elif experiment == "multiplex_amplicon_sequencing":
            info, priority, risk, runtime, reuse = ("high" if top_count > 5 else "medium"), 12, "medium", "medium", "high"
            reason = "parallel targeted sequencing may cover several candidates when single-amplicon Sanger scales poorly"
        elif experiment == "target_capture":
            info, priority, risk, runtime, reuse = "high", 14, "medium", "high", "high"
            reason = "capture can broaden locus coverage when many coding or selected regulatory targets remain"
        elif experiment == "expression_assay":
            info, priority, risk, runtime, reuse = "medium", 25, "medium", "medium", "medium"
            reason = "expression can prioritize mechanisms only when the relevant tissue and stage are available"
        elif experiment == "annotation_rescue":
            info = "high" if annotation_first else "low"
            priority = 2 if annotation_first else 35
            risk, runtime, reuse = "low", "medium", "high"
            if annotation_first:
                rules.append("R_ANNOTATION_GAP_FAVORS_RESCUE")
                facts.append(f"annotation_gap_hypotheses={annotation_count}")
                reason = "a reviewed annotation gap (reference-missing locus) favors comparative annotation rescue before further coding sequencing"
        elif experiment == "non_gene_follow_up":
            info = "high" if non_gene_first else "low"
            priority = 3 if non_gene_first else 45
            risk, runtime, reuse = "medium", "medium", "high"
            if non_gene_first:
                rules.append("R_NON_GENE_LOCUS_FAVORS_FOLLOW_UP")
                facts.append(f"non_gene_locus_hypotheses={non_gene_count}")
                reason = "linked-QTL, low-recombination, or non-gene locus evidence favors association/haplotype or multi-locus follow-up over coding-gene continuation"

        if state == "complete":
            priority += 50
            reason = "this action is already complete; interpret its outcome before repeating it"
        feasibility = "blocked" if blocking else ("review_required" if review else "feasible")
        if feasibility == "blocked":
            priority += 100
        elif feasibility == "review_required":
            priority += 30
        option_rows.append({
            "experiment_id": experiment,
            "action_rank": 0,
            "feasibility_status": feasibility,
            "hypotheses_addressed": HYPOTHESES[experiment],
            "candidates_potentially_eliminated": elimination_count if experiment == "additional_grid_markers" else "unknown",
            "elimination_basis": elimination_basis if experiment == "additional_grid_markers" else "not calculated from observed breakpoints",
            "expected_information_class": info,
            "estimated_cost_class": cost,
            "cost_source": cost_source,
            "estimated_amplicon_count": str(amplicon_counts.get(experiment, "unknown")),
            "sample_requirement": ";".join(SAMPLE_REQUIREMENTS[experiment]),
            "sample_available": ";".join(sample_states),
            "runtime_class": runtime,
            "technical_failure_risk": risk,
            "downstream_reusability": reuse,
            "assay_state": state,
            "assay_outcome": outcome,
            "blocking_dependencies": ";".join(sorted(blocking)),
            "input_facts": ";".join(facts + ([f"failure_reason={failure_reason}"] if failure_reason else [])),
            "rules_fired": ";".join(sorted(rules)) or "R_DEFAULT_COMPARISON",
            "recommendation_reason": reason,
            "human_approval": "pending",
            "_sort": (priority, 0 if feasibility == "feasible" else 1 if feasibility == "review_required" else 2, INFORMATION_ORDER[info], COST_ORDER[cost], experiment),
            "_reviews": review,
        })

    option_rows.sort(key=lambda row: row["_sort"])
    for index, row in enumerate(option_rows, 1):
        row["action_rank"] = index
        for reason in row.pop("_reviews"):
            review_rows.append({"experiment_id": str(row["experiment_id"]), "review_reason": reason, "required_input": reason.split(":", 1)[0], "decision_effect": "recommendation remains provisional or blocked until reviewed"})
        row.pop("_sort")
    if remaining_ids is not None and not remaining_ids:
        # Issue #285: zero survivors is an explicit state, not a claim that no
        # causal gene exists. Nothing is resurrected.
        review_rows.append({
            "experiment_id": "no_current_gene_candidate",
            "review_reason": "remaining_candidates_empty",
            "required_input": "remaining_candidates",
            "decision_effect": "no current gene candidate survives elimination; do not resurrect any candidate without new evidence",
        })
    options = pd.DataFrame(option_rows)
    dependencies = pd.DataFrame(dependency_rows).sort_values(["experiment_id", "dependency"], kind="stable")
    reviews = pd.DataFrame(review_rows, columns=["experiment_id", "review_reason", "required_input", "decision_effect"]).drop_duplicates().sort_values(["experiment_id", "review_reason"], kind="stable")

    options_path = out_dir / "next_experiment_options.tsv"
    dependencies_path = out_dir / "experiment_dependencies.tsv"
    reviews_path = out_dir / "experiment_review_required.tsv"
    decision_path = out_dir / "next_experiment_decision.json"
    options.to_csv(options_path, sep="\t", index=False, lineterminator="\n")
    dependencies.to_csv(dependencies_path, sep="\t", index=False, lineterminator="\n")
    reviews.to_csv(reviews_path, sep="\t", index=False, lineterminator="\n")
    primary = option_rows[0]
    inputs = {"ranked_candidates": _sha256(ranked_path), "case_manifest": _sha256(manifest_path)}
    inputs.update({name: _sha256(path) for name, path in optional_paths.items() if path is not None})
    decision = {
        "schema_version": "1.0",
        "policy_version": POLICY_VERSION,
        "decision_type": "rule_based_ordinal_experiment_support",
        "primary_recommendation": {key: primary[key] for key in ("experiment_id", "action_rank", "feasibility_status", "recommendation_reason", "rules_fired")},
        "recommendation_status": "ready_for_human_approval" if primary["feasibility_status"] == "feasible" else "provisional_review_required",
        "alternatives": [{"experiment_id": row["experiment_id"], "action_rank": row["action_rank"], "feasibility_status": row["feasibility_status"]} for row in option_rows[1:]],
        "human_approval": "pending",
        "causal_probability": False,
        "bayesian_information_gain": False,
        "automatic_wet_execution": False,
        "unknown_values_imputed": False,
        "input_fingerprints_sha256": inputs,
    }
    decision_path.write_text(json.dumps(decision, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if ranked_path.read_bytes() != ranked_bytes:
        raise RuntimeError("biological ranked-candidate input was modified")
    return {
        "options": str(options_path),
        "decision": str(decision_path),
        "dependencies": str(dependencies_path),
        "review_required": str(reviews_path),
        "primary_recommendation": str(primary["experiment_id"]),
        "human_approval": "pending",
    }
