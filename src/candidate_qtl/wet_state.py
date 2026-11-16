from __future__ import annotations

import hashlib
import json

import pandas as pd

from .core import runtime_path
from .schemas import (
    ASSAY_STATES,
    CANDIDATE_WET_STATES,
    WET_ASSAY_COLUMNS as ASSAY_COLUMNS,
    WET_COVERAGE_COLUMNS as COVERAGE_COLUMNS,
    WET_HISTORY_COLUMNS as HISTORY_COLUMNS,
    WET_RESULT_COLUMNS as RESULT_COLUMNS,
)


CANDIDATE_TRANSITIONS = {
    "not_tested": CANDIDATE_WET_STATES - {"not_tested"},
    "partially_covered": {"coding_region_cleared", "no_parental_difference", "candidate_variant_found_unconfirmed", "candidate_variant_confirmed", "inconclusive"},
    "no_parental_difference": {"partially_covered", "coding_region_cleared", "candidate_variant_found_unconfirmed", "candidate_variant_confirmed", "inconclusive"},
    "candidate_variant_found_unconfirmed": {"candidate_variant_confirmed", "inconclusive"},
    "candidate_variant_confirmed": {"fails_cosegregation", "excluded_by_recombinant", "supports_qtl", "requires_independent_validation", "inconclusive"},
    "supports_qtl": {"requires_independent_validation", "inconclusive"},
    "requires_independent_validation": {"fails_cosegregation", "excluded_by_recombinant", "inconclusive"},
    "coding_region_cleared": {"candidate_variant_found_unconfirmed", "candidate_variant_confirmed", "inconclusive"},
    "inconclusive": CANDIDATE_WET_STATES - {"not_tested", "inconclusive"},
    "fails_cosegregation": set(),
    "excluded_by_recombinant": set(),
}
ALLOWED_TRANSITIONS = {
    "not_planned": {"planned", "primer_designed", "ordered", "wet_pending", "failed", "partial_result", "completed", "needs_redesign"},
    "planned": {"primer_designed", "ordered", "wet_pending", "failed", "needs_redesign"},
    "primer_designed": {"ordered", "wet_pending", "failed", "needs_redesign"},
    "ordered": {"wet_pending", "failed", "needs_redesign"},
    "wet_pending": {"failed", "partial_result", "completed", "needs_redesign"},
    "failed": {"needs_redesign"},
    "partial_result": {"wet_pending", "completed", "needs_redesign"},
    "needs_redesign": {"planned", "primer_designed"},
    "completed": set(),
}
REVIEW_COLUMNS = ["gene_id", "assay_id", "review_reason", "required_action", "source_ref"]


def _read(path_arg: str | None, columns: list[str]) -> pd.DataFrame:
    if not path_arg:
        return pd.DataFrame(columns=columns)
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    for column in columns:
        if column not in frame.columns:
            frame[column] = ""
    return frame[columns]


def _require_nonempty(frame: pd.DataFrame, columns: tuple[str, ...], name: str) -> None:
    for column in columns:
        if column not in frame.columns or frame[column].astype(str).str.strip().eq("").any():
            raise ValueError(f"{name} requires non-empty {column}")


def _stable_id(*values: object) -> str:
    payload = "\x1f".join(str(value) for value in values)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:20]


def _yes(value: object) -> bool:
    return str(value).strip().casefold() in {"yes", "true", "1", "required"}


def _pass(value: object) -> bool:
    return str(value).strip().casefold() in {"pass", "passed", "high", "good", "reviewed_pass"}


def _primer_adapter(path_arg: str | None) -> pd.DataFrame:
    if not path_arg:
        return pd.DataFrame(columns=ASSAY_COLUMNS)
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    rows = []
    for index, row in frame.iterrows():
        assay_id = str(row.get("marker_id", "") or f"primer:{row.get('gene_id', '')}:{index + 1}")
        design_text = " ".join(str(row.get(column, "")) for column in ("pair_specificity", "coverage_status", "validation_status")).casefold()
        redesign = any(token in design_text for token in ("non-specific", "nonspecific", "coverage_gap", "no_pair"))
        state = "needs_redesign" if redesign else "primer_designed"
        rows.append({
            "assay_id": assay_id, "gene_id": row.get("gene_id", ""), "transcript_id": "", "hypothesis_id": "",
            "assay_type": "primerblast_marker", "planned_region": row.get("target_interval", ""), "region_class": "marker",
            "assay_state": state, "failure_reason": design_text if redesign else "", "source_ref": row.get("source_file", path_arg),
            "review_status": "unreviewed", "reviewer": "",
        })
    return pd.DataFrame(rows, columns=ASSAY_COLUMNS)


def _exon_design_adapter(path_arg: str | None) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not path_arg:
        return pd.DataFrame(columns=ASSAY_COLUMNS), pd.DataFrame(columns=COVERAGE_COLUMNS)
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    required = {"gene_id", "target_interval", "status"}
    if missing := sorted(required.difference(frame.columns)):
        raise ValueError(f"exon coverage adapter is missing: {', '.join(missing)}")
    assays, coverage = [], []
    for index, row in frame.iterrows():
        assay_id = f"exon-design:{row['gene_id']}:{index + 1}"
        status = str(row["status"])
        state = "primer_designed" if status == "covered" else "needs_redesign" if status == "coverage_gap" else "planned"
        assays.append({
            "assay_id": assay_id, "gene_id": row["gene_id"], "transcript_id": row.get("transcript_id", ""),
            "hypothesis_id": "", "assay_type": "exon_sanger", "planned_region": row["target_interval"],
            "region_class": "cds", "assay_state": state, "failure_reason": "" if status == "covered" else status,
            "source_ref": row.get("output_json", path_arg), "review_status": "unreviewed", "reviewer": "",
        })
        coverage.append({
            "assay_id": assay_id, "gene_id": row["gene_id"], "transcript_id": row.get("transcript_id", ""),
            "planned_region": row["target_interval"], "actual_covered_region": "",
            "coverage_fraction": "0", "region_class": "cds", "required_region": "yes",
            "quality_status": "design_only", "read_direction": "", "off_target_status": "unknown",
            "multi_product_status": "unknown", "source_ref": row.get("output_json", path_arg),
        })
    return pd.DataFrame(assays, columns=ASSAY_COLUMNS), pd.DataFrame(coverage, columns=COVERAGE_COLUMNS)


def _transition_history(
    assays: pd.DataFrame,
    previous: pd.DataFrame,
    reviews: list[dict[str, str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    history = previous.copy()
    for column in HISTORY_COLUMNS:
        if column not in history.columns:
            history[column] = ""
    history = history[HISTORY_COLUMNS]
    next_index = int(pd.to_numeric(history["event_index"], errors="coerce").max()) + 1 if len(history) else 1
    latest: dict[str, str] = {}
    if len(history):
        ordered = history.assign(_order=pd.to_numeric(history["event_index"], errors="coerce").fillna(0)).sort_values("_order", kind="stable")
        assay_events = ordered[(ordered["entity_type"].isin({"", "assay"})) & ordered["assay_id"].ne("")]
        latest = dict(zip(assay_events["assay_id"], assay_events["to_state"], strict=False))
    new_rows = []
    for index, row in assays.iterrows():
        assay_id = str(row["assay_id"])
        requested = str(row["assay_state"]).strip()
        if requested not in ASSAY_STATES:
            raise ValueError(f"unsupported assay state: {requested}")
        if assay_id not in latest:
            status, reason, effective = "accepted", "initial_import", requested
            from_state = "not_planned"
        else:
            from_state = latest[assay_id]
            if requested == from_state:
                continue
            allowed = requested in ALLOWED_TRANSITIONS.get(from_state, set())
            status = "accepted" if allowed else "rejected_forbidden_transition"
            reason = "allowed_transition" if allowed else f"forbidden:{from_state}->{requested}"
            effective = requested if allowed else from_state
            if not allowed:
                reviews.append({"gene_id": str(row["gene_id"]), "assay_id": assay_id, "review_reason": "forbidden_assay_state_transition", "required_action": reason, "source_ref": str(row["source_ref"])})
        assays.at[index, "assay_state"] = effective
        new_rows.append({
            "history_id": "hist_" + _stable_id(assay_id, next_index, from_state, requested, row["source_ref"]),
            "event_index": next_index, "entity_type": "assay", "entity_id": assay_id, "assay_id": assay_id, "gene_id": row["gene_id"], "from_state": from_state,
            "to_state": effective, "transition_status": status, "transition_reason": reason, "source_ref": row["source_ref"],
            "review_status": row["review_status"], "reviewer": row["reviewer"],
        })
        latest[assay_id] = effective
        next_index += 1
    combined = pd.concat([history, pd.DataFrame(new_rows, columns=HISTORY_COLUMNS)], ignore_index=True)
    return assays, combined


def _candidate_states(
    assays: pd.DataFrame,
    coverage: pd.DataFrame,
    results: pd.DataFrame,
    min_coverage: float,
    required_read_direction: str,
    reviews: list[dict[str, str]],
) -> pd.DataFrame:
    genes = sorted(set(assays["gene_id"]) | set(coverage["gene_id"]) | set(results["gene_id"]))
    rows = []
    for gene in genes:
        gene_assays = assays[assays["gene_id"].eq(gene)]
        gene_coverage = coverage[coverage["gene_id"].eq(gene)].copy()
        gene_results = results[results["gene_id"].eq(gene)]
        assay_ids = ";".join(sorted(set(gene_assays["assay_id"])))
        hypothesis_ids = ";".join(sorted(value for value in set(gene_assays["hypothesis_id"]) if value))
        state, coding_state, regulatory_state = "not_tested", "not_tested", "not_tested"
        variant_state, validation_state, negative_scope = "not_tested", "not_tested", ""
        reason = "no Wet result supplied"
        fractions = pd.to_numeric(gene_coverage["coverage_fraction"], errors="coerce") if len(gene_coverage) else pd.Series(dtype=float)
        min_fraction = float(fractions.min()) if fractions.notna().any() else None
        required_assays = set(gene_assays[gene_assays["region_class"].isin({"cds", "splice_boundary"})]["assay_id"])
        completed_required = set(gene_assays[gene_assays["assay_state"].eq("completed")]["assay_id"])
        wet_coverage = gene_coverage[gene_coverage["quality_status"].map(_pass)] if len(gene_coverage) else gene_coverage
        if required_read_direction != "any" and len(wet_coverage):
            wet_coverage = wet_coverage[wet_coverage["read_direction"].str.casefold().eq(required_read_direction)]
        covered_required = set(wet_coverage[pd.to_numeric(wet_coverage["coverage_fraction"], errors="coerce").ge(min_coverage)]["assay_id"])
        coverage_gate = bool(required_assays) and required_assays.issubset(covered_required)
        if len(gene_coverage):
            state, coding_state, reason = "partially_covered", "partially_covered", "actual coverage does not clear every required coding/splice assay"
        if gene_assays["assay_state"].isin({"failed", "needs_redesign"}).any():
            state, reason = "inconclusive", "assay failure or redesign is not negative evidence"
        result_text = " ".join(gene_results["result_value"].astype(str)).casefold()
        failure_text = " ".join(gene_assays["failure_reason"].astype(str)).casefold()
        non_amplification = any(token in result_text + " " + failure_text for token in ("non-ampl", "non_ampl", "no amplification", "failed amplification"))
        if non_amplification:
            regulatory_state = "sv_pav_or_primer_mismatch_warning"
            state, reason = "inconclusive", "parent non-amplification may indicate SV/PAV or primer mismatch"
            reviews.append({"gene_id": gene, "assay_id": assay_ids, "review_reason": "parent_non_amplification_not_negative", "required_action": "redesign primers and evaluate SV/PAV before candidate exclusion", "source_ref": ";".join(gene_results["source_ref"])})
        reviewed_quality = gene_results["review_status"].str.casefold().eq("reviewed") & gene_results["quality_status"].map(_pass) & gene_results["source_ref"].str.strip().ne("")
        validated_results = gene_results[reviewed_quality]
        confirmed = validated_results["variant_confirmation_status"].str.casefold().eq("confirmed").any()
        unconfirmed = gene_results["variant_confirmation_status"].str.casefold().isin({"unconfirmed", "predicted", "candidate"}).any()
        if gene_results["variant_confirmation_status"].str.casefold().eq("confirmed").any() and not confirmed:
            reviews.append({"gene_id": gene, "assay_id": assay_ids, "review_reason": "variant_confirmation_without_reviewed_quality_source", "required_action": "supply reviewed orthogonal confirmation with quality and source", "source_ref": ";".join(gene_results["source_ref"])})
        if confirmed:
            state, variant_state, reason = "candidate_variant_confirmed", "confirmed", "a parental variant is confirmed but causality remains unproven"
        elif unconfirmed:
            state, variant_state, reason = "candidate_variant_found_unconfirmed", "unconfirmed", "a candidate variant requires orthogonal confirmation"
        validations = set(validated_results["validation_result"].str.casefold())
        unvalidated_claims = set(gene_results["validation_result"].str.casefold()) - {""} - validations
        if unvalidated_claims:
            reviews.append({"gene_id": gene, "assay_id": assay_ids, "review_reason": "genetic_validation_without_reviewed_quality_source", "required_action": "retain as unconfirmed until individual-level reviewed evidence is supplied", "source_ref": ";".join(gene_results["source_ref"])})
        if "excluded_by_recombinant" in validations:
            state, validation_state, reason = "excluded_by_recombinant", "excluded_by_recombinant", "reviewed recombinant breakpoint excludes the candidate locus"
        elif "fails_cosegregation" in validations:
            state, validation_state, reason = "fails_cosegregation", "fails_cosegregation", "candidate marker fails cosegregation; distinguish this from recombinant exclusion"
        elif "supports_qtl" in validations:
            state, validation_state, reason = "requires_independent_validation", "supports_qtl", "cosegregation supports the QTL but is not causal proof"
        negative_mask = gene_results["result_type"].str.casefold().eq("parental_difference") & gene_results["result_value"].str.casefold().isin({"none", "absent", "no_difference", "no parental difference"})
        no_difference = negative_mask.any()
        reviewed_negative = bool((negative_mask & reviewed_quality).any())
        reviewed_negative_rows = gene_results[negative_mask & reviewed_quality]
        negative_assays = {
            assay_id
            for assay_id, group in reviewed_negative_rows.groupby("assay_id", sort=False)
            if {"parent_a", "parent_b"}.issubset(set(group["sample_role"].str.casefold()))
        }
        negative_gate = bool(required_assays) and required_assays.issubset(negative_assays) and required_assays.issubset(completed_required)
        protected_state = state in {"candidate_variant_confirmed", "candidate_variant_found_unconfirmed", "excluded_by_recombinant", "fails_cosegregation", "requires_independent_validation"}
        if reviewed_negative and coverage_gate and negative_gate and not protected_state:
            state, coding_state = "coding_region_cleared", "coding_region_cleared"
            negative_scope = ";".join(sorted(set(wet_coverage["actual_covered_region"])))
            reason = "all declared required coding/splice assays pass coverage and quality with no parental difference"
            regulatory_state = regulatory_state if regulatory_state != "not_tested" else "not_tested_retained"
        elif no_difference and not protected_state:
            state = "no_parental_difference"
            negative_scope = ";".join(sorted(value for value in set(gene_coverage["actual_covered_region"]) if value))
            reason = "no parental difference is limited to the measured region; every required completed assay needs reviewed parent A and B negatives"
            reviews.append({"gene_id": gene, "assay_id": assay_ids, "review_reason": "negative_result_without_complete_coverage_quality_gate", "required_action": "retain Reserve/regulatory/SV hypotheses and complete required coverage", "source_ref": ";".join(gene_results["source_ref"])})
        requested_actions = set(gene_results["requested_candidate_action"].str.casefold())
        genetic_exclusion = state in {"excluded_by_recombinant", "fails_cosegregation"}
        reserve_allowed = "yes_genetic_exclusion" if genetic_exclusion else "scope_limited_coding_only" if state == "coding_region_cleared" else "no"
        reject_allowed = "yes_genetic_exclusion" if genetic_exclusion else "no"
        if requested_actions & {"reserve", "reject", "exclude"} and reserve_allowed == "no":
            reviews.append({"gene_id": gene, "assay_id": assay_ids, "review_reason": "candidate_transition_requested_without_coverage_quality_gate", "required_action": "complete coverage/quality or provide reviewed genetic exclusion", "source_ref": ";".join(gene_results["source_ref"])})
        rows.append({
            "gene_id": gene, "hypothesis_ids": hypothesis_ids, "candidate_state": state, "coding_region_state": coding_state,
            "regulatory_sv_state": regulatory_state, "minimum_actual_coverage_fraction": "" if min_fraction is None else round(min_fraction, 6),
            "coverage_quality_gate": "pass" if coverage_gate else "fail", "variant_state": variant_state,
            "validation_state": validation_state, "negative_result_scope": negative_scope,
            "reserve_transition_allowed": reserve_allowed, "hypothesis_reject_allowed": reject_allowed,
            "state_reason": reason, "source_assay_ids": assay_ids,
        })
    return pd.DataFrame(rows)


def _append_candidate_history(
    states: pd.DataFrame,
    history: pd.DataFrame,
    reviews: list[dict[str, str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    next_index = int(pd.to_numeric(history["event_index"], errors="coerce").max()) + 1 if len(history) else 1
    candidate_events = history[history["entity_type"].eq("candidate")].copy()
    latest: dict[str, str] = {}
    if len(candidate_events):
        candidate_events["_order"] = pd.to_numeric(candidate_events["event_index"], errors="coerce").fillna(0)
        candidate_events = candidate_events.sort_values("_order", kind="stable")
        latest = dict(zip(candidate_events["entity_id"], candidate_events["to_state"], strict=False))
    new_rows = []
    for index, row in states.iterrows():
        gene = str(row["gene_id"])
        requested = str(row["candidate_state"])
        previous = latest.get(gene, "not_tested")
        if requested == previous:
            continue
        allowed = gene not in latest or requested in CANDIDATE_TRANSITIONS.get(previous, set())
        effective = requested if allowed else previous
        status = "accepted_derived" if allowed else "rejected_forbidden_transition"
        reason = str(row["state_reason"]) if allowed else f"forbidden_candidate_transition:{previous}->{requested}"
        if not allowed:
            states.at[index, "candidate_state"] = effective
            reviews.append({"gene_id": gene, "assay_id": str(row["source_assay_ids"]), "review_reason": "forbidden_candidate_state_transition", "required_action": reason, "source_ref": "wet_state_history"})
        new_rows.append({
            "history_id": "hist_" + _stable_id("candidate", gene, next_index, previous, requested),
            "event_index": next_index, "entity_type": "candidate", "entity_id": gene, "assay_id": "", "gene_id": gene,
            "from_state": previous, "to_state": effective, "transition_status": status,
            "transition_reason": reason, "source_ref": str(row["source_assay_ids"]), "review_status": "derived_unreviewed", "reviewer": "",
        })
        latest[gene] = effective
        next_index += 1
    return states, pd.concat([history, pd.DataFrame(new_rows, columns=HISTORY_COLUMNS)], ignore_index=True)


def update_wet_state(
    out_dir_arg: str,
    *,
    assay_registry_arg: str | None = None,
    coverage_arg: str | None = None,
    results_arg: str | None = None,
    previous_history_arg: str | None = None,
    primer_results_arg: str | None = None,
    exon_coverage_arg: str | None = None,
    min_coverage: float = 1.0,
    required_read_direction: str = "bidirectional",
) -> dict[str, object]:
    if not 0 < min_coverage <= 1:
        raise ValueError("min_coverage must be within (0, 1]")
    if required_read_direction not in {"bidirectional", "forward", "reverse", "any"}:
        raise ValueError("required_read_direction must be bidirectional, forward, reverse, or any")
    assays = _read(assay_registry_arg, ASSAY_COLUMNS)
    coverage = _read(coverage_arg, COVERAGE_COLUMNS)
    results = _read(results_arg, RESULT_COLUMNS)
    primer_assays = _primer_adapter(primer_results_arg)
    exon_assays, exon_design_coverage = _exon_design_adapter(exon_coverage_arg)
    assays = pd.concat([assays, primer_assays, exon_assays], ignore_index=True).drop_duplicates(subset=["assay_id"], keep="first")
    coverage = pd.concat([coverage, exon_design_coverage], ignore_index=True).drop_duplicates()
    if assays.empty:
        raise ValueError("at least one assay registry, PrimerBlast result, or exon coverage input is required")
    _require_nonempty(assays, ("assay_id", "gene_id", "assay_type", "assay_state"), "assay registry")
    if assays["assay_id"].duplicated().any():
        raise ValueError("assay_id must be unique")
    unknown_states = sorted(set(assays["assay_state"]) - ASSAY_STATES)
    if unknown_states:
        raise ValueError(f"unsupported assay states: {', '.join(unknown_states)}")
    if len(coverage):
        _require_nonempty(coverage, ("assay_id", "gene_id", "planned_region", "region_class", "quality_status"), "assay coverage")
        numeric = pd.to_numeric(coverage["coverage_fraction"], errors="coerce")
        if numeric.isna().any() or numeric.lt(0).any() or numeric.gt(1).any():
            raise ValueError("coverage_fraction must be numeric within [0, 1]")
        unknown_assays = sorted(set(coverage["assay_id"]) - set(assays["assay_id"]))
        if unknown_assays:
            raise ValueError(f"coverage references unknown assay_id: {', '.join(unknown_assays)}")
        assay_genes = assays.set_index("assay_id")["gene_id"].to_dict()
        mismatched = coverage[coverage.apply(lambda row: assay_genes.get(row["assay_id"]) != row["gene_id"], axis=1)]
        if len(mismatched):
            raise ValueError("coverage gene_id must match its assay registry gene_id")
        actual_pass = coverage[coverage["quality_status"].map(_pass)]
        if actual_pass["actual_covered_region"].str.strip().eq("").any():
            raise ValueError("quality-passing coverage requires actual_covered_region")
    if len(results):
        _require_nonempty(results, ("result_id", "assay_id", "gene_id", "sample_role", "result_type", "source_ref"), "wet results")
        for index, row in results.iterrows():
            if not str(row["validation_result"]).strip():
                recombinant = str(row["recombinant_result"]).strip().casefold()
                cosegregation = str(row["cosegregation_result"]).strip().casefold()
                if recombinant == "excluded_by_recombinant":
                    results.at[index, "validation_result"] = recombinant
                elif cosegregation in {"fails_cosegregation", "supports_qtl"}:
                    results.at[index, "validation_result"] = cosegregation
        allowed_roles = {"parent_a", "parent_b", "f1", "f2", "f2_3", "recombinant", "other"}
        unknown_roles = sorted(set(results["sample_role"].str.casefold()) - allowed_roles)
        if unknown_roles:
            raise ValueError(f"unsupported sample roles: {', '.join(unknown_roles)}")
        individual_roles = results["sample_role"].str.casefold().isin({"f2", "f2_3", "recombinant"})
        if results.loc[individual_roles, "individual_id"].str.strip().eq("").any():
            raise ValueError("F2/F2:3/recombinant results require individual_id")
        validation_values = results["validation_result"].str.casefold()
        unknown_validations = sorted(set(validation_values) - {"", "supports_qtl", "fails_cosegregation", "excluded_by_recombinant"})
        if unknown_validations:
            raise ValueError(f"unsupported validation_result: {', '.join(unknown_validations)}")
        genetic_rows = validation_values.ne("")
        if (~results.loc[genetic_rows, "sample_role"].str.casefold().isin({"f2", "f2_3", "recombinant"})).any():
            raise ValueError("genetic validation requires F2/F2:3/recombinant individual-level result")
        recombinant_claims = validation_values.eq("excluded_by_recombinant")
        if (~results.loc[recombinant_claims, "sample_role"].str.casefold().eq("recombinant")).any():
            raise ValueError("excluded_by_recombinant requires sample_role=recombinant")
        unknown_assays = sorted(set(results["assay_id"]) - set(assays["assay_id"]))
        if unknown_assays:
            raise ValueError(f"results reference unknown assay_id: {', '.join(unknown_assays)}")
        assay_genes = assays.set_index("assay_id")["gene_id"].to_dict()
        mismatched = results[results.apply(lambda row: assay_genes.get(row["assay_id"]) != row["gene_id"], axis=1)]
        if len(mismatched):
            raise ValueError("result gene_id must match its assay registry gene_id")
    reviews: list[dict[str, str]] = []
    previous = _read(previous_history_arg, HISTORY_COLUMNS)
    history_ids = previous["history_id"]
    if history_ids[history_ids.ne("")].duplicated().any():
        raise ValueError("previous history contains duplicate history_id")
    assays, history = _transition_history(assays, previous, reviews)
    candidate_states = _candidate_states(assays, coverage, results, min_coverage, required_read_direction, reviews)
    candidate_states, history = _append_candidate_history(candidate_states, history, reviews)
    review_frame = pd.DataFrame(reviews, columns=REVIEW_COLUMNS).drop_duplicates().sort_values(["gene_id", "assay_id", "review_reason"], kind="stable")
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "assay_registry": out_dir / "assay_registry.tsv",
        "assay_target_coverage": out_dir / "assay_target_coverage.tsv",
        "wet_results": out_dir / "wet_results.tsv",
        "candidate_wet_state": out_dir / "candidate_wet_state.tsv",
        "wet_state_history": out_dir / "wet_state_history.tsv",
        "wet_review_required": out_dir / "wet_review_required.tsv",
    }
    assays.sort_values(["gene_id", "assay_id"], kind="stable").to_csv(paths["assay_registry"], sep="\t", index=False, lineterminator="\n")
    coverage.sort_values(["gene_id", "assay_id", "planned_region"], kind="stable").to_csv(paths["assay_target_coverage"], sep="\t", index=False, lineterminator="\n")
    results.sort_values(["gene_id", "assay_id", "result_id"], kind="stable").to_csv(paths["wet_results"], sep="\t", index=False, lineterminator="\n")
    candidate_states.sort_values(["gene_id"], kind="stable").to_csv(paths["candidate_wet_state"], sep="\t", index=False, lineterminator="\n")
    history.assign(_order=pd.to_numeric(history["event_index"], errors="coerce").fillna(0)).sort_values("_order", kind="stable").drop(columns=["_order"]).to_csv(paths["wet_state_history"], sep="\t", index=False, lineterminator="\n")
    review_frame.to_csv(paths["wet_review_required"], sep="\t", index=False, lineterminator="\n")
    manifest = {
        "schema_version": "1.1", "minimum_coverage_fraction": min_coverage, "required_read_direction": required_read_direction,
        "principles": [
            "assay state and candidate/mechanism state are separate",
            "assay failure and parent non-amplification are not negative evidence",
            "coding clearance is scope-limited and retains regulatory/SV hypotheses",
            "cosegregation support is not causal proof",
            "history is append-only; forbidden transitions are rejected and reviewed",
        ],
        "outputs": {key: str(value) for key, value in paths.items()},
    }
    manifest_path = out_dir / "wet_state_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {key: str(value) for key, value in paths.items()} | {"manifest": str(manifest_path), "candidate_count": len(candidate_states), "review_count": len(review_frame)}
