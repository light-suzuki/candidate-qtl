from __future__ import annotations

import hashlib
import json

import pandas as pd

from .core import runtime_path
from .schemas import (
    ANNOTATION_RESCUE_HYPOTHESIS_ESCALATION_CLASSES,
    HYPOTHESIS_STATUSES,
    HYPOTHESIS_TYPES,
    LOCUS_HYPOTHESIS_COLUMNS as HYPOTHESIS_COLUMNS,
    PROPOSAL_ORIGINS,
)

# dual_route.py's dual_route_compare() convergence_status values that represent a
# family-level A/B disagreement: one route reached the family, the other did not
# (dual_route.py, _convergence_status/dual_route_compare, ~L890-1102). A row that
# both routes agree on (mechanistic/family/gene-overlap/process/text convergence)
# is not disagreement and is left alone here.
_DUAL_ROUTE_DISCORDANT_STATUSES = frozenset({"phenotype_only_supported", "interval_only_supported"})


TRANSITIONS = {
    "proposed": HYPOTHESIS_STATUSES - {"proposed", "rejected"},
    "pending_review": {"proposed", "ready_to_test", "inconclusive", "superseded"},
    "ready_to_test": {"testing", "partially_tested", "weakened", "inconclusive", "superseded"},
    "testing": {"partially_tested", "supported", "weakened", "rejected", "inconclusive", "requires_independent_validation"},
    "partially_tested": {"testing", "supported", "weakened", "rejected", "inconclusive", "requires_independent_validation"},
    "supported": {"weakened", "rejected", "requires_independent_validation", "superseded"},
    "weakened": {"testing", "partially_tested", "rejected", "inconclusive", "superseded"},
    "inconclusive": {"ready_to_test", "testing", "partially_tested", "weakened", "superseded"},
    "requires_independent_validation": {"supported", "weakened", "rejected", "inconclusive"},
    "rejected": set(),
    "superseded": set(),
}
EVIDENCE_LINK_COLUMNS = ["link_id", "hypothesis_id", "evidence_kind", "evidence_id", "direction", "source_table", "link_status"]
EXPERIMENT_LINK_COLUMNS = ["link_id", "hypothesis_id", "experiment_id", "relationship", "link_status"]
HISTORY_COLUMNS = [
    "history_id", "event_index", "hypothesis_id", "from_status", "requested_status", "to_status", "transition_status",
    "status_reason", "evidence_ids", "wet_result_ids", "source_id", "review_status", "reviewer",
]
REVIEW_COLUMNS = ["hypothesis_id", "review_reason", "required_action", "source_id"]


def _read(path_arg: str | None, columns: list[str] | None = None) -> pd.DataFrame:
    if not path_arg:
        return pd.DataFrame(columns=columns or [])
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    if columns:
        for column in columns:
            if column not in frame.columns:
                frame[column] = ""
        frame = frame[columns]
    return frame


def _stable_id(prefix: str, *values: object) -> str:
    payload = "\x1f".join(str(value) for value in values)
    return f"{prefix}_{hashlib.sha256(payload.encode('utf-8')).hexdigest()[:20]}"


def _tokens(value: object) -> list[str]:
    return sorted(set(token.strip() for token in str(value).replace("|", ";").split(";") if token.strip()))


def _gene_members(row: pd.Series) -> set[str]:
    return set(_tokens(row["member_gene_ids"])) | set(_tokens(row["target_gene_ids"]))


def _review(reviews: list[dict[str, str]], hypothesis_id: str, reason: str, action: str, source: str) -> None:
    reviews.append({"hypothesis_id": hypothesis_id, "review_reason": reason, "required_action": action, "source_id": source})


def _normalize_hypotheses(frame: pd.DataFrame, case_id: str, reviews: list[dict[str, str]]) -> pd.DataFrame:
    for column in HYPOTHESIS_COLUMNS:
        if column not in frame.columns:
            frame[column] = ""
    frame = frame[HYPOTHESIS_COLUMNS].copy()
    for index, row in frame.iterrows():
        row_case = str(row["case_id"]).strip() or case_id
        hypothesis_type = str(row["hypothesis_type"]).strip()
        if hypothesis_type not in HYPOTHESIS_TYPES:
            raise ValueError(f"unsupported hypothesis_type: {hypothesis_type}")
        origin = str(row["proposal_origin"]).strip().casefold()
        if not origin:
            raise ValueError("proposal_origin is required and must explicitly be manual, workflow, or llm")
        if origin not in PROPOSAL_ORIGINS:
            raise ValueError(f"unsupported proposal_origin: {origin}")
        hypothesis_id = str(row["hypothesis_id"]).strip() or _stable_id(
            "hyp", row_case, hypothesis_type, row["interval_or_element_id"], row["member_gene_ids"],
            row["target_gene_ids"], row["proposed_mechanism"], origin,
        )
        status = str(row["status"]).strip().casefold() or ("pending_review" if origin == "llm" else "proposed")
        requested_status = status
        status_guardrail = ""
        if status not in HYPOTHESIS_STATUSES:
            raise ValueError(f"unsupported hypothesis status: {status}")
        review_status = str(row["review_status"]).strip().casefold() or "pending_review"
        if origin == "llm" and (status != "pending_review" or review_status not in {"pending_review", "unreviewed"}):
            _review(reviews, hypothesis_id, "llm_proposal_cannot_self_promote", "human must verify source and mechanism before status promotion", str(row["source_id"]))
            status, review_status = "pending_review", "pending_review"
            status_guardrail = "llm_proposal_cannot_self_promote"
        genes = _gene_members(row)
        if hypothesis_type == "multi_gene_combination" and len(genes) < 2:
            raise ValueError("multi_gene_combination requires at least two member/target genes")
        non_gene = hypothesis_type in {"annotation_missing", "linked_qtl", "intergenic_regulatory", "unknown_locus_mechanism"}
        if not genes and not str(row["interval_or_element_id"]).strip() and not non_gene:
            raise ValueError(f"{hypothesis_type} requires a gene member/target or interval_or_element_id")
        if non_gene and not str(row["interval_or_element_id"]).strip():
            _review(reviews, hypothesis_id, "non_gene_hypothesis_without_interval", "supply a locus interval or element ID", str(row["source_id"]))
        if not str(row["source_id"]).strip():
            _review(reviews, hypothesis_id, "hypothesis_source_missing", "add a stable publication/database/local artifact reference", "")
        frame.at[index, "hypothesis_id"] = hypothesis_id
        frame.at[index, "case_id"] = row_case
        frame.at[index, "proposal_origin"] = origin
        frame.at[index, "status"] = status
        frame.at[index, "requested_status"] = requested_status
        frame.at[index, "status_guardrail"] = status_guardrail
        frame.at[index, "review_status"] = review_status
        frame.at[index, "created_by"] = str(row["created_by"]).strip() or origin
        frame.at[index, "gene_disposition"] = str(row["gene_disposition"]).strip() or "retained"
        frame.at[index, "rejection_scope"] = str(row["rejection_scope"]).strip() or "none"
    if frame["hypothesis_id"].duplicated().any():
        raise ValueError("hypothesis_id must be unique")
    known = set(frame["hypothesis_id"])
    for row in frame.itertuples(index=False):
        for relation in (row.parent_hypothesis_id, row.supersedes_hypothesis_id):
            if relation and relation not in known:
                _review(reviews, row.hypothesis_id, "related_hypothesis_missing", f"add related hypothesis {relation} to the same ledger", row.source_id)
    return frame


def _workflow_branches(
    hypotheses: pd.DataFrame,
    wet_states: pd.DataFrame,
    wet_results: pd.DataFrame,
    sv_events: pd.DataFrame,
    annotation_rescue: pd.DataFrame,
    dual_route_comparison: pd.DataFrame,
    case_id: str,
    reviews: list[dict[str, str]],
) -> pd.DataFrame:
    rows = hypotheses.to_dict(orient="records")
    known_keys = {(row["hypothesis_type"], tuple(sorted(_gene_members(pd.Series(row))))) for row in rows}

    # --- Structural-variant / PAV events (sv-plan, multigenome-private) -> sv_cnv_pav (issue #230.2/.3) ---
    # ``sv_events`` reuses the exact table shape ``_build_links`` below already reads
    # (one of sv_event_id/call_id/cluster_id, plus gene_id): sv-plan's canonical_sv.tsv
    # and the multigenome-private adapter both produce that shape directly, so no
    # parallel schema is introduced here. A gene with at least one such event
    # deterministically gets one sv_cnv_pav hypothesis, and the underlying event ids
    # become context_evidence_ids so _build_links resolves them from the same table
    # instead of copying payloads.
    if not sv_events.empty and "gene_id" in sv_events.columns:
        sv_id_column = next((column for column in ("sv_event_id", "call_id", "cluster_id") if column in sv_events.columns), None)
        if sv_id_column:
            events_by_gene: dict[str, list[str]] = {}
            for record in sv_events.to_dict(orient="records"):
                gene = str(record.get("gene_id", "")).strip()
                event_id = str(record.get(sv_id_column, "")).strip()
                if gene and event_id:
                    events_by_gene.setdefault(gene, []).append(event_id)
            for gene in sorted(events_by_gene):
                key = ("sv_cnv_pav", (gene,))
                if key in known_keys:
                    continue
                child = {column: "" for column in HYPOTHESIS_COLUMNS}
                child.update({
                    "hypothesis_id": _stable_id("hyp", case_id, gene, "sv_cnv_pav", "sv_events"), "case_id": case_id,
                    "hypothesis_type": "sv_cnv_pav", "member_gene_ids": gene,
                    "proposed_mechanism": f"structural variant/PAV event overlapping {gene} or its local regulatory context",
                    "status": "proposed", "status_reason": "workflow seed from sv-plan/multigenome-private structural event evidence",
                    "created_by": "workflow:sv_events", "proposal_origin": "workflow", "review_status": "pending_review",
                    "source_id": f"sv_events:{gene}", "gene_disposition": "retained", "rejection_scope": "none",
                    "context_evidence_ids": ";".join(_tokens(";".join(events_by_gene[gene]))), "requested_status": "proposed",
                })
                rows.append(child)
                known_keys.add(key)

    # --- Annotation rescue -> annotation_missing (issue #230.1) ---
    # native_discordant / reference_rescue / ambiguous_overlap rows mean the native
    # annotation disagrees with, or is entirely missing for, a rescue projection:
    # evidence of a possible annotation gap, never evidence that no causal locus
    # exists there.
    if not annotation_rescue.empty and "model_class" in annotation_rescue.columns:
        for record in annotation_rescue.to_dict(orient="records"):
            model_class = str(record.get("model_class", "")).strip()
            query_id = str(record.get("query_protein_id", "")).strip()
            if model_class not in ANNOTATION_RESCUE_HYPOTHESIS_ESCALATION_CLASSES or not query_id:
                continue
            projection_id = str(record.get("projection_id", "")).strip()
            native_genes = _tokens(record.get("native_gene_ids", ""))
            chrom, start, end = str(record.get("chrom", "")), str(record.get("start", "")), str(record.get("end", ""))
            element_id = projection_id or (f"{chrom}:{start}-{end}" if chrom and start and end else query_id)
            key = ("annotation_missing", tuple(sorted(native_genes)) if native_genes else (f"rescue:{query_id}:{projection_id}",))
            if key in known_keys:
                continue
            child = {column: "" for column in HYPOTHESIS_COLUMNS}
            child.update({
                "hypothesis_id": _stable_id("hyp", case_id, "annotation_missing", query_id, projection_id, model_class),
                "case_id": case_id, "hypothesis_type": "annotation_missing",
                "interval_or_element_id": element_id, "target_gene_ids": ";".join(sorted(native_genes)),
                "proposed_mechanism": f"annotation-rescue projection of query {query_id} classified {model_class}; native annotation may be missing, incomplete, or discordant",
                "status": "proposed", "status_reason": f"workflow seed from annotation-rescue model_class={model_class}",
                "created_by": "workflow:annotation_rescue", "proposal_origin": "workflow", "review_status": "pending_review",
                "source_id": f"annotation_rescue:{query_id}:{projection_id or 'unmapped'}", "gene_disposition": "retained",
                "rejection_scope": "none", "requested_status": "proposed",
            })
            rows.append(child)
            known_keys.add(key)

    # --- Route A/B disagreement -> linked_qtl / unknown_locus_mechanism (issue #230.4) ---
    # dual_route_compare's per-row convergence_status of phenotype_only_supported /
    # interval_only_supported IS the family-level A/B disagreement: one route
    # reached the family, the other did not. Route B (interval-first) only ever
    # proposes genes that are actually present in the interval, so
    # interval_only_supported guarantees interval membership and becomes
    # unknown_locus_mechanism.
    #
    # phenotype_only_supported only says "the interval-first route did not reach
    # this hypothesis" (dual_route.py's own note on that row). That is NOT the same
    # claim as "the family is absent from the interval": a thin or incomplete
    # interval annotation produces the identical status. Calling every such row a
    # linked_qtl would assert a separate causal locus that the comparison never
    # established, so linked_qtl is used only where the row carries the evidence
    # for it -- Route A's own anchoring guard reporting ``unanchored``, i.e. the
    # family is missing from the interval listing itself. Everything else stays
    # unknown_locus_mechanism: a real escalation, without a fabricated location.
    if not dual_route_comparison.empty and "convergence_status" in dual_route_comparison.columns:
        for index, record in enumerate(dual_route_comparison.to_dict(orient="records")):
            status = str(record.get("convergence_status", "")).strip()
            if status not in _DUAL_ROUTE_DISCORDANT_STATUSES:
                continue
            family = str(record.get("gene_family", "")).strip()
            comparison_id = str(record.get("comparison_id", "")).strip() or f"row{index}"
            if status == "interval_only_supported":
                hypothesis_type = "unknown_locus_mechanism"
                genes = _tokens(record.get("shared_interval_genes", ""))
                mechanism = (
                    f"interval-first (Route B) reverse lookup located {family or 'a candidate family'} in the QTL "
                    "interval with no independently corroborated phenotype-first (Route A) mechanism"
                )
                reason = "interval_only_supported_no_phenotype_link"
            elif str(record.get("route_a_anchor_status", "")).strip() == "unanchored":
                hypothesis_type = "linked_qtl"
                genes = set()
                mechanism = (
                    f"phenotype-first (Route A) proposed {family or 'a candidate family'}, which the interval "
                    "listing does not contain (anchor_status=unanchored); if it is causal it sits outside this "
                    "interval"
                )
                reason = "phenotype_only_supported_family_absent_from_interval"
            else:
                hypothesis_type = "unknown_locus_mechanism"
                genes = set()
                mechanism = (
                    f"phenotype-first (Route A) proposed {family or 'a candidate family'}; the interval-first "
                    "(Route B) reverse lookup did not independently reach it. Interval membership is neither "
                    "confirmed nor refuted here -- a thin interval annotation produces the same signal"
                )
                reason = "phenotype_only_supported_not_corroborated_by_interval_route"
            element_id = family or f"dual_route:{comparison_id}"
            key = (hypothesis_type, tuple(sorted(genes)) if genes else (f"family:{element_id}",))
            if key in known_keys:
                continue
            child = {column: "" for column in HYPOTHESIS_COLUMNS}
            child.update({
                "hypothesis_id": _stable_id("hyp", case_id, hypothesis_type, family, comparison_id, "dual_route"),
                "case_id": case_id, "hypothesis_type": hypothesis_type,
                "interval_or_element_id": element_id, "member_gene_ids": ";".join(sorted(genes)),
                "proposed_mechanism": mechanism, "status": "proposed",
                "status_reason": f"workflow escalation from A/B dual-route disagreement ({reason})",
                "created_by": "workflow:dual_route_compare", "proposal_origin": "workflow", "review_status": "pending_review",
                "source_id": f"dual_route_compare:{comparison_id}", "gene_disposition": "retained", "rejection_scope": "none",
                "requested_status": "proposed",
            })
            rows.append(child)
            known_keys.add(key)

    if wet_states.empty or "gene_id" not in wet_states.columns:
        result = pd.DataFrame(rows, columns=HYPOTHESIS_COLUMNS)
        for row in result.itertuples(index=False):
            if row.status == "weakened" and row.rejection_scope != "coding_mechanism_only":
                _review(reviews, row.hypothesis_id, "coding_negative_scope_not_preserved", "retain gene and record coding_mechanism_only scope", row.source_id)
        return result
    for wet in wet_states.to_dict(orient="records"):
        gene = str(wet.get("gene_id", "")).strip()
        if not gene:
            continue
        candidate_state = str(wet.get("candidate_state", ""))
        regulatory_state = str(wet.get("regulatory_sv_state", ""))
        source = f"wet_candidate_state:{gene}"
        gene_results = wet_results[wet_results["gene_id"].eq(gene)] if "gene_id" in wet_results.columns else wet_results.iloc[0:0]
        negative_ids: list[str] = []
        non_amplification_ids: list[str] = []
        if len(gene_results) and "result_id" in gene_results.columns:
            result_types = gene_results["result_type"].str.casefold() if "result_type" in gene_results.columns else pd.Series("", index=gene_results.index)
            result_values = gene_results["result_value"].str.casefold() if "result_value" in gene_results.columns else pd.Series("", index=gene_results.index)
            negative_mask = result_types.eq("parental_difference") & result_values.isin({"none", "absent", "no_difference", "no parental difference"})
            non_amp_mask = result_values.str.contains("non-ampl|non_ampl|no amplification|failed amplification", regex=True)
            negative_ids = _tokens(";".join(gene_results.loc[negative_mask, "result_id"]))
            non_amplification_ids = _tokens(";".join(gene_results.loc[non_amp_mask, "result_id"]))
        coding_rows = [row for row in rows if row["hypothesis_type"] in {"single_gene_coding", "single_gene_splice"} and gene in _gene_members(pd.Series(row))]
        if not coding_rows:
            coding_rows.append({column: "" for column in HYPOTHESIS_COLUMNS})
            coding_rows[-1].update({
                "hypothesis_id": _stable_id("hyp", case_id, gene, "single_gene_coding"), "case_id": case_id,
                "hypothesis_type": "single_gene_coding", "member_gene_ids": gene,
                "proposed_mechanism": f"coding change in {gene}", "status": "proposed",
                "status_reason": "workflow seed from Wet-tested candidate", "created_by": "workflow:wet_state",
                "proposal_origin": "workflow", "review_status": "pending_review", "source_id": source,
                "gene_disposition": "retained", "rejection_scope": "none", "requested_status": "proposed",
            })
            rows.extend(coding_rows)
            known_keys.add(("single_gene_coding", (gene,)))
        if candidate_state == "coding_region_cleared":
            for row in coding_rows:
                if row["status"] != "rejected":
                    row["status"] = "weakened"
                    row["requested_status"] = "weakened"
                row["status_reason"] = "required coding/splice scope was quality-cleared; this does not reject the gene"
                row["gene_disposition"] = "retained"
                row["rejection_scope"] = "coding_mechanism_only"
                row["counterevidence_ids"] = ";".join(_tokens(f"{row['counterevidence_ids']};{';'.join(negative_ids)}"))
            for hypothesis_type, mechanism in (
                ("promoter_or_utr", f"regulatory change affecting {gene}"),
                ("sv_cnv_pav", f"SV/CNV/PAV affecting {gene} or its local regulatory context"),
            ):
                key = (hypothesis_type, (gene,))
                if key not in known_keys:
                    parent = coding_rows[0]["hypothesis_id"]
                    child = {column: "" for column in HYPOTHESIS_COLUMNS}
                    child.update({
                        "hypothesis_id": _stable_id("hyp", case_id, gene, hypothesis_type), "case_id": case_id,
                        "hypothesis_type": hypothesis_type, "member_gene_ids": gene,
                        "proposed_mechanism": mechanism, "status": "proposed",
                        "status_reason": "workflow branch retained after scoped coding-negative result",
                        "created_by": "workflow:wet_state", "proposal_origin": "workflow",
                        "review_status": "pending_review", "source_id": source,
                        "parent_hypothesis_id": parent, "gene_disposition": "retained", "rejection_scope": "none",
                        "context_evidence_ids": ";".join(negative_ids), "requested_status": "proposed",
                    })
                    rows.append(child)
                    known_keys.add(key)
        if "sv_pav" in regulatory_state or "primer_mismatch" in regulatory_state:
            key = ("sv_cnv_pav", (gene,))
            if key not in known_keys:
                child = {column: "" for column in HYPOTHESIS_COLUMNS}
                child.update({
                    "hypothesis_id": _stable_id("hyp", case_id, gene, "sv_cnv_pav"), "case_id": case_id,
                    "hypothesis_type": "sv_cnv_pav", "member_gene_ids": gene,
                    "proposed_mechanism": f"SV/PAV or primer-site divergence near {gene}", "status": "proposed",
                    "status_reason": "parent non-amplification is a branch warning, not negative evidence",
                    "created_by": "workflow:wet_state", "proposal_origin": "workflow", "review_status": "pending_review",
                    "source_id": source, "gene_disposition": "retained", "rejection_scope": "none",
                    "context_evidence_ids": ";".join(non_amplification_ids), "requested_status": "proposed",
                })
                rows.append(child)
                known_keys.add(key)
    result = pd.DataFrame(rows, columns=HYPOTHESIS_COLUMNS)
    for row in result.itertuples(index=False):
        if row.status == "weakened" and row.rejection_scope != "coding_mechanism_only":
            _review(reviews, row.hypothesis_id, "coding_negative_scope_not_preserved", "retain gene and record coding_mechanism_only scope", row.source_id)
    return result


def _validate_relationship_graph(hypotheses: pd.DataFrame, reviews: list[dict[str, str]]) -> None:
    known = set(hypotheses["hypothesis_id"])
    graph: dict[str, list[str]] = {}
    for row in hypotheses.itertuples(index=False):
        relations = [relation for relation in (row.parent_hypothesis_id, row.supersedes_hypothesis_id) if relation]
        graph[row.hypothesis_id] = relations
        for relation in relations:
            if relation not in known:
                _review(reviews, row.hypothesis_id, "related_hypothesis_missing", f"add related hypothesis {relation} to the same ledger", row.source_id)

    active: set[str] = set()
    complete: set[str] = set()

    def visit(hypothesis_id: str) -> None:
        if hypothesis_id in active:
            raise ValueError(f"hypothesis relationship cycle detected at {hypothesis_id}")
        if hypothesis_id in complete:
            return
        active.add(hypothesis_id)
        for relation in graph.get(hypothesis_id, []):
            if relation in known:
                visit(relation)
        active.remove(hypothesis_id)
        complete.add(hypothesis_id)

    for hypothesis_id in sorted(known):
        visit(hypothesis_id)


def _enforce_rejection_gate(
    hypotheses: pd.DataFrame,
    wet_states: pd.DataFrame,
    reviews: list[dict[str, str]],
) -> pd.DataFrame:
    state_by_gene = {
        str(row.get("gene_id", "")): str(row.get("candidate_state", ""))
        for row in wet_states.to_dict(orient="records")
    }
    for index, row in hypotheses.iterrows():
        if row["status"] != "rejected":
            continue
        genes = _gene_members(row)
        genetically_excluded = bool(genes) and all(state_by_gene.get(gene) == "excluded_by_recombinant" for gene in genes)
        if genetically_excluded:
            hypotheses.at[index, "gene_disposition"] = "rejected"
            hypotheses.at[index, "rejection_scope"] = "locus_excluded_by_recombinant"
            continue
        coding_only = row["hypothesis_type"] in {"single_gene_coding", "single_gene_splice"}
        coding_scope_cleared = coding_only and bool(genes) and all(
            state_by_gene.get(gene) == "coding_region_cleared" for gene in genes
        )
        if coding_scope_cleared:
            hypotheses.at[index, "gene_disposition"] = "retained"
            hypotheses.at[index, "rejection_scope"] = "coding_mechanism_only"
            continue
        hypotheses.at[index, "status"] = "weakened" if coding_only else "inconclusive"
        hypotheses.at[index, "status_guardrail"] = "hypothesis_rejection_without_genetic_exclusion_gate"
        hypotheses.at[index, "gene_disposition"] = "retained"
        hypotheses.at[index, "rejection_scope"] = "coding_mechanism_only" if coding_only else "none"
        _review(
            reviews, str(row["hypothesis_id"]), "hypothesis_rejection_without_genetic_exclusion_gate",
            "retain the gene/hypothesis until reviewed recombinant exclusion or equivalent scoped evidence is supplied",
            str(row["source_id"]),
        )
    return hypotheses


def _build_links(
    hypotheses: pd.DataFrame,
    evidence: pd.DataFrame,
    wet_results: pd.DataFrame,
    sv_events: pd.DataFrame,
    reviews: list[dict[str, str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    evidence_column = next((column for column in ("observation_id", "evidence_id") if column in evidence.columns), None)
    evidence_ids = set(evidence[evidence_column]) if evidence_column else set()
    wet_ids = set(wet_results["result_id"]) if "result_id" in wet_results.columns else set()
    sv_id_column = next((column for column in ("sv_event_id", "call_id", "cluster_id") if column in sv_events.columns), None)
    sv_ids = set(sv_events[sv_id_column]) if sv_id_column else set()
    collisions = (evidence_ids & wet_ids) | (evidence_ids & sv_ids) | (wet_ids & sv_ids)
    if collisions:
        raise ValueError(f"evidence ID namespace collision: {', '.join(sorted(collisions))}")
    evidence_genes = dict(zip(evidence[evidence_column], evidence["gene_id"], strict=False)) if evidence_column and "gene_id" in evidence.columns else {}
    wet_genes = dict(zip(wet_results["result_id"], wet_results["gene_id"], strict=False)) if "result_id" in wet_results.columns and "gene_id" in wet_results.columns else {}
    sv_genes = dict(zip(sv_events[sv_id_column], sv_events["gene_id"], strict=False)) if sv_id_column and "gene_id" in sv_events.columns else {}
    evidence_rows, experiment_rows = [], []
    for row in hypotheses.itertuples(index=False):
        hypothesis_genes = _gene_members(pd.Series(row._asdict()))
        for direction, raw in (
            ("support", row.supporting_evidence_ids),
            ("against", row.counterevidence_ids),
            ("context", row.context_evidence_ids),
        ):
            for evidence_id in _tokens(raw):
                if evidence_id in evidence_ids:
                    kind, source, status = "canonical_evidence", "evidence", "resolved"
                    linked_gene = str(evidence_genes.get(evidence_id, ""))
                elif evidence_id in wet_ids:
                    kind, source, status = "wet_result", "wet_results", "resolved"
                    linked_gene = str(wet_genes.get(evidence_id, ""))
                elif evidence_id in sv_ids:
                    kind, source, status = "sv_event", "structural_variants", "resolved"
                    linked_gene = str(sv_genes.get(evidence_id, ""))
                else:
                    kind, source, status = "unresolved", "unknown", "review_required"
                    linked_gene = ""
                    _review(reviews, row.hypothesis_id, "evidence_reference_unresolved", f"provide canonical evidence or Wet result ID {evidence_id}", row.source_id)
                if linked_gene and hypothesis_genes and linked_gene not in hypothesis_genes:
                    status = "review_required"
                    _review(
                        reviews, row.hypothesis_id, "evidence_gene_mismatch",
                        f"evidence {evidence_id} belongs to {linked_gene}, outside hypothesis members",
                        row.source_id,
                    )
                evidence_rows.append({
                    "link_id": _stable_id("helink", row.hypothesis_id, direction, evidence_id),
                    "hypothesis_id": row.hypothesis_id, "evidence_kind": kind, "evidence_id": evidence_id,
                    "direction": direction, "source_table": source, "link_status": status,
                })
        for experiment_id in _tokens(row.required_experiment_ids):
            experiment_rows.append({
                "link_id": _stable_id("hxlink", row.hypothesis_id, experiment_id),
                "hypothesis_id": row.hypothesis_id, "experiment_id": experiment_id,
                "relationship": "required_to_test", "link_status": "declared",
            })
    return (
        pd.DataFrame(evidence_rows, columns=EVIDENCE_LINK_COLUMNS),
        pd.DataFrame(experiment_rows, columns=EXPERIMENT_LINK_COLUMNS),
    )


def _append_history(
    hypotheses: pd.DataFrame,
    previous: pd.DataFrame,
    evidence_links: pd.DataFrame,
    reviews: list[dict[str, str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    for column in HISTORY_COLUMNS:
        if column not in previous.columns:
            previous[column] = ""
    previous = previous[HISTORY_COLUMNS]
    history_ids = previous["history_id"]
    if history_ids[history_ids.ne("")].duplicated().any():
        raise ValueError("previous hypothesis history contains duplicate history_id")
    ordered = previous.assign(_order=pd.to_numeric(previous["event_index"], errors="coerce").fillna(0)).sort_values("_order", kind="stable")
    latest = dict(zip(ordered["hypothesis_id"], ordered["to_status"], strict=False))
    next_index = int(ordered["_order"].max()) + 1 if len(ordered) else 1
    new_rows = []
    for index, row in hypotheses.iterrows():
        hypothesis_id = str(row["hypothesis_id"])
        requested = str(row["status"])
        raw_requested = str(row["requested_status"] or requested)
        guardrail = str(row["status_guardrail"])
        old = latest.get(hypothesis_id, "")
        if requested == old and raw_requested == requested:
            continue
        allowed = not old or requested in TRANSITIONS.get(old, set())
        effective = requested if allowed else old
        status = (
            "rejected_guardrail_transition"
            if guardrail or raw_requested != requested
            else "accepted_initial" if not old else "accepted" if allowed else "rejected_forbidden_transition"
        )
        reason = guardrail or (str(row["status_reason"]) if allowed else f"forbidden:{old}->{requested}")
        if not allowed:
            hypotheses.at[index, "status"] = effective
            _review(reviews, hypothesis_id, "forbidden_hypothesis_status_transition", reason, str(row["source_id"]))
        links = evidence_links[evidence_links["hypothesis_id"].eq(hypothesis_id)] if len(evidence_links) else evidence_links
        canonical = ";".join(sorted(links[links["evidence_kind"].eq("canonical_evidence")]["evidence_id"]))
        wet = ";".join(sorted(links[links["evidence_kind"].eq("wet_result")]["evidence_id"]))
        new_rows.append({
            "history_id": _stable_id("hhist", hypothesis_id, next_index, old, requested),
            "event_index": next_index, "hypothesis_id": hypothesis_id, "from_status": old,
            "requested_status": raw_requested, "to_status": effective, "transition_status": status, "status_reason": reason,
            "evidence_ids": canonical, "wet_result_ids": wet, "source_id": row["source_id"],
            "review_status": row["review_status"], "reviewer": "",
        })
        latest[hypothesis_id] = effective
        next_index += 1
    return hypotheses, pd.concat([previous, pd.DataFrame(new_rows, columns=HISTORY_COLUMNS)], ignore_index=True)


def update_hypothesis_ledger(
    out_dir_arg: str,
    *,
    hypotheses_arg: str | None = None,
    case_id: str = "",
    evidence_arg: str | None = None,
    wet_results_arg: str | None = None,
    sv_events_arg: str | None = None,
    wet_candidate_state_arg: str | None = None,
    annotation_rescue_arg: str | None = None,
    dual_route_comparison_arg: str | None = None,
    previous_history_arg: str | None = None,
    ranked_arg: str | None = None,
) -> dict[str, object]:
    reviews: list[dict[str, str]] = []
    hypotheses = _read(hypotheses_arg)
    wet_states = _read(wet_candidate_state_arg)
    evidence = _read(evidence_arg)
    wet_results = _read(wet_results_arg)
    sv_events = _read(sv_events_arg)
    # Deterministic non-gene escalation inputs: annotation-rescue's
    # projection audit and dual_route_compare's A/B comparison are read the same
    # way wet_candidate_state is -- optional, and absence degrades to simply not
    # generating that branch's hypotheses rather than failing the run.
    annotation_rescue = _read(annotation_rescue_arg)
    dual_route_comparison = _read(dual_route_comparison_arg)
    previous = _read(previous_history_arg, HISTORY_COLUMNS)
    ranked_path = runtime_path(ranked_arg) if ranked_arg else None
    ranked_bytes = ranked_path.read_bytes() if ranked_path else None
    if hypotheses.empty and wet_states.empty and sv_events.empty and annotation_rescue.empty and dual_route_comparison.empty:
        raise ValueError("provide --hypotheses and/or --wet-candidate-state")
    effective_case = case_id or (str(hypotheses.iloc[0].get("case_id", "")) if len(hypotheses) else "")
    if not effective_case:
        raise ValueError("case_id is required when it is absent from the hypothesis table")
    hypotheses = _normalize_hypotheses(hypotheses, effective_case, reviews)
    hypotheses = _workflow_branches(hypotheses, wet_states, wet_results, sv_events, annotation_rescue, dual_route_comparison, effective_case, reviews)
    _validate_relationship_graph(hypotheses, reviews)
    hypotheses = _enforce_rejection_gate(hypotheses, wet_states, reviews)
    evidence_links, experiment_links = _build_links(hypotheses, evidence, wet_results, sv_events, reviews)
    hypotheses, history = _append_history(hypotheses, previous, evidence_links, reviews)
    review_frame = pd.DataFrame(reviews, columns=REVIEW_COLUMNS).drop_duplicates().sort_values(
        ["hypothesis_id", "review_reason"], kind="stable"
    )
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "locus_hypotheses": out_dir / "locus_hypotheses.tsv",
        "hypothesis_evidence_links": out_dir / "hypothesis_evidence_links.tsv",
        "hypothesis_experiment_links": out_dir / "hypothesis_experiment_links.tsv",
        "hypothesis_history": out_dir / "hypothesis_history.tsv",
        "hypothesis_review_required": out_dir / "hypothesis_review_required.tsv",
    }
    hypotheses.sort_values(["case_id", "hypothesis_id"], kind="stable").to_csv(paths["locus_hypotheses"], sep="\t", index=False, lineterminator="\n")
    evidence_links.sort_values(["hypothesis_id", "direction", "evidence_id"], kind="stable").to_csv(paths["hypothesis_evidence_links"], sep="\t", index=False, lineterminator="\n")
    experiment_links.sort_values(["hypothesis_id", "experiment_id"], kind="stable").to_csv(paths["hypothesis_experiment_links"], sep="\t", index=False, lineterminator="\n")
    history.assign(_order=pd.to_numeric(history["event_index"], errors="coerce").fillna(0)).sort_values("_order", kind="stable").drop(columns="_order").to_csv(paths["hypothesis_history"], sep="\t", index=False, lineterminator="\n")
    review_frame.to_csv(paths["hypothesis_review_required"], sep="\t", index=False, lineterminator="\n")
    if ranked_path and ranked_path.read_bytes() != ranked_bytes:
        raise RuntimeError("deterministic ranked-candidate input was modified")
    manifest = {
        "schema_version": "1.0", "case_id": effective_case,
        "rank_is_read_only": bool(ranked_path), "rank_overwritten": False,
        "hypothesis_count": len(hypotheses), "review_count": len(review_frame),
        "principles": [
            "gene rank and locus-mechanism hypotheses are separate",
            "links reference canonical evidence/Wet IDs without copying evidence payloads",
            "LLM proposals remain pending review",
            "coding-negative status weakens only the coding mechanism and retains the gene",
            "annotation-rescue discordant/rescue/ambiguous rows, sv-plan/multigenome-private structural events, and "
            "A/B dual-route disagreement deterministically seed non-gene locus hypotheses without touching the gene score",
        ],
        "outputs": {key: str(value) for key, value in paths.items()},
    }
    manifest_path = out_dir / "hypothesis_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {key: str(value) for key, value in paths.items()} | {
        "manifest": str(manifest_path), "hypothesis_count": len(hypotheses), "review_count": len(review_frame),
    }
