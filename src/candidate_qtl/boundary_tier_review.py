"""Scope existing Wet review to a versioned relevance pool; never rank within a tier.

Pure planning plus an adapter around adaptive_synthesis's existing independent
review consumer. This module does not run a model or change support states.
Callers must supply an already validated, relevance-annotated interval roster.
"""
from __future__ import annotations

import csv
import hashlib
import html
import json
import shutil
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any
from urllib.parse import quote

from .relevance_tiers import RELEVANCE_TIERS, T1, T2, NEW_TIERS

CONTRACT = "boundary-tier-review-1.0.0"
EVIDENCE_SCOPED_POLICY = "evidence_scoped_v2"
UPPER_EVIDENCE_SCOPED_POLICY = "upper_evidence_scoped_v3"
EVIDENCE_SCOPED_POLICIES = (EVIDENCE_SCOPED_POLICY, UPPER_EVIDENCE_SCOPED_POLICY)
DEFAULT_SELECTION_POLICY = UPPER_EVIDENCE_SCOPED_POLICY
REPORTED_UPPER_TIERS = frozenset((T1, T2, *NEW_TIERS[:3]))
LEGACY_POLICY = "legacy_whole_tier_v1"
REQUEST_CONTRACT = "boundary-tier-review-request-1.0.0"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _id(row: Mapping[str, Any]) -> str:
    return str(row.get("candidate_gene_id") or row.get("gene_id") or "").strip()


def _validated(rows: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows = [dict(row) for row in rows]
    ids = [_id(row) for row in rows]
    if any(not gene for gene in ids) or len(ids) != len(set(ids)):
        raise ValueError("boundary roster requires non-empty unique candidate IDs")
    for row in rows:
        if row.get("relevance_tier") not in RELEVANCE_TIERS:
            raise ValueError("boundary review requires canonical relevance_tier; legacy support_tier is not a substitute")
    # This controls serialization only. No within-tier ordering admits/excludes
    # a member: the planner always moves an entire tier at once.
    return sorted(rows, key=lambda row: (RELEVANCE_TIERS.index(row["relevance_tier"]), _id(row)))


def partition_relevance_tiers(
    rows: Iterable[Mapping[str, Any]], *, target: int = 20, upper: int = 25,
    frozen_gene_ids: Iterable[str] = (), protected_gene_ids: Iterable[str] = (),
    policy_version: str = LEGACY_POLICY, eligible_gene_ids: Iterable[str] = (),
    active_gene_ids: Iterable[str] | None = None,
) -> dict[str, Any]:
    """Build a legacy whole-tier or evidence-scoped comparison boundary.

    Legacy mode preserves the historic continuation rule. Evidence-scoped v2
    retains all T1 plus the already eligible candidate pool and sends the first
    overflowing tier and every lower eligible seed into one comparison scope.
    T1 is always frozen in full, including when T1 alone exceeds upper.
    V3 prefers reported upper evidence and compares only the first overflowing
    tier, preserving original eligibility for later tiers in the ledger.
    A supplied active scope locks post-review planning to the original pool.
    """
    if isinstance(target, bool) or isinstance(upper, bool) or not isinstance(target, int) or not isinstance(upper, int) or not 0 < target <= upper:
        raise ValueError("require integer 0 < target <= upper")
    if policy_version not in {LEGACY_POLICY, *EVIDENCE_SCOPED_POLICIES}:
        raise ValueError("unknown boundary selection policy")
    records = _validated(rows)
    ids = {_id(row) for row in records}
    prefix = set(frozen_gene_ids)
    if not prefix <= ids:
        raise ValueError("frozen prefix contains unknown candidates")
    protected = set(protected_gene_ids)
    if not protected <= ids:
        raise ValueError("protected candidates contain unknown candidates")
    eligible = set(eligible_gene_ids)
    if policy_version in EVIDENCE_SCOPED_POLICIES and not eligible <= ids:
        raise ValueError("eligible candidate IDs contain unknown candidates")
    if policy_version == LEGACY_POLICY and eligible:
        raise ValueError("legacy boundary selection does not accept eligible candidate IDs")
    if policy_version in EVIDENCE_SCOPED_POLICIES:
        # Evidence-scoped review compares already-proposed candidates only.
        # Known-trait T1 rows remain unconditional rescue; lower identity-only
        # inventory rows stay auditable in lower, never fill a short set.
        eligible = (eligible | protected) & ids
    known = {_id(row) for row in records if row["relevance_tier"] in {T1, NEW_TIERS[0]}}
    active = eligible
    upper_priority_applied = False
    if policy_version == UPPER_EVIDENCE_SCOPED_POLICY:
        reported = ({_id(row) for row in records if row["relevance_tier"] in REPORTED_UPPER_TIERS}
                    & (eligible | known)) - protected
        upper_priority_applied = bool(reported)
        active = (reported | protected) if reported else eligible
        if active_gene_ids is not None:
            active = set(active_gene_ids)
            if not active <= eligible | known:
                raise ValueError("active candidate IDs exceed original eligibility")
    elif active_gene_ids is not None:
        raise ValueError("active candidate scope requires upper evidence policy")
    # Run the ordinary whole-tier selector first. Protected C rows are then
    # appended from whichever normal partition they landed in; this prevents
    # rescue rows from creating extra ordinary fill or displacing A/B tiers.
    prefix -= protected
    prefix |= known
    frozen = [row for row in records if _id(row) in prefix]
    remaining = [row for row in records if _id(row) not in prefix]
    validation_remaining = [row for row in remaining if _id(row) not in protected
                            and (policy_version == LEGACY_POLICY or _id(row) in active)]
    if frozen and validation_remaining:
        last = max(RELEVANCE_TIERS.index(row["relevance_tier"]) for row in frozen)
        if any(RELEVANCE_TIERS.index(row["relevance_tier"]) <= last for row in validation_remaining):
            raise ValueError("frozen prefix would skip unreviewed members at a higher/equal tier")
    boundary: list[dict[str, Any]] = []
    for tier in RELEVANCE_TIERS[1:]:
        if policy_version == LEGACY_POLICY and len(frozen) >= target:
            break
        group = [row for row in remaining if row["relevance_tier"] == tier]
        if policy_version in EVIDENCE_SCOPED_POLICIES:
            group = [row for row in group if _id(row) in active]
            if not group:
                continue
        # C rescue may exceed capacity, but must not create an A/B review
        # boundary solely by being added to an otherwise admissible tier.
        ordinary_count = sum(_id(row) not in protected for row in frozen + group)
        if ordinary_count > upper:
            boundary = (group if policy_version == UPPER_EVIDENCE_SCOPED_POLICY else
                        [row for row in remaining
                         if RELEVANCE_TIERS.index(row["relevance_tier"]) >= RELEVANCE_TIERS.index(tier)
                         and _id(row) in active]
                        if policy_version in EVIDENCE_SCOPED_POLICIES else group)
            break
        frozen.extend(group)
    admitted = {_id(row) for row in frozen + boundary}
    lower = [row for row in records if _id(row) not in admitted]
    frozen_ids = {_id(row) for row in frozen}
    protected_rows = [row for row in records if _id(row) in protected - frozen_ids]
    frozen.extend(protected_rows)
    protected_ids = {_id(row) for row in protected_rows}
    boundary = [row for row in boundary if _id(row) not in protected_ids]
    lower = [row for row in lower if _id(row) not in protected_ids]
    if policy_version == UPPER_EVIDENCE_SCOPED_POLICY:
        active = {_id(row) for row in frozen + boundary}
    state = ("boundary_review_required" if boundary else
             "known_rescue_over_upper_retained" if len(known) > upper else
             "reviewed_over_upper_retained" if len(frozen) > upper else
             "within_target_range" if len(frozen) >= target else "under_target_pool_exhausted")
    return {
        "frozen": frozen, "boundary": boundary, "lower": lower,
        "boundary_tier": boundary[0]["relevance_tier"] if boundary else None,
        "target": target, "upper": upper, "state": state,
        "frozen_gene_ids": sorted(_id(row) for row in frozen),
        "rank_contribution": "none",
        "policy_version": policy_version,
        "eligible_gene_ids": sorted(eligible) if policy_version in EVIDENCE_SCOPED_POLICIES else [],
        **({"active_gene_ids": sorted(active), "upper_priority_applied": upper_priority_applied}
           if policy_version == UPPER_EVIDENCE_SCOPED_POLICY else {}),
    }


def _read(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    # Frozen provenance JSON can exceed csv's small default field limit. A
    # field cannot contain more characters than the entire UTF-8 file has
    # bytes, so raise the process-local parser limit to the actual input size.
    # Never lower it: another reader in this process may need the larger limit.
    csv.field_size_limit(max(csv.field_size_limit(), path.stat().st_size))
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = list(reader.fieldnames or [])
        if not columns or len(columns) != len(set(columns)):
            raise ValueError("candidate table requires unique column names")
        rows = list(reader)
        if any(None in row or any(value is None for value in row.values()) for row in rows):
            raise ValueError("malformed candidate table row")
        return columns, rows


def _write(path: Path, columns: list[str], rows: Iterable[Mapping[str, Any]]) -> dict[str, str]:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return {"path": str(path.resolve()), "sha256": _sha(path)}


def _descriptor(path: Path) -> dict[str, str]:
    return {"path": str(path.resolve()), "sha256": _sha(path)}


def _checked(entry: Mapping[str, Any]) -> Path:
    path = Path(str(entry.get("path", "")))
    if not path.is_absolute() or not path.is_file() or _sha(path) != entry.get("sha256"):
        raise ValueError("boundary artifact hash mismatch or missing")
    return path


def _wet_summary_cells(
    row: Mapping[str, str], disposition: str, decision: Mapping[str, str],
) -> list[str]:
    """Render saved facts only; never join independently aggregated endpoints."""
    from .observation_tiers import gene_binding
    from .paper_reports import source_tokens
    from .source_verify import classify_citation

    def text(value: Any) -> str:
        return " ".join(str(value if value is not None else "").split())

    def first(record: Mapping, *keys: str) -> str:
        return next((text(record.get(key)) for key in keys if text(record.get(key))), "")

    def cell(value: Any) -> str:
        value = html.escape(text(value), quote=False).replace("\\", "\\\\")
        for char in "|[]`*_":
            value = value.replace(char, "\\" + char)
        return value

    def records(value: Any) -> list[dict]:
        try:
            value = json.loads(value) if isinstance(value, str) and value else value
        except (ValueError, TypeError):
            return []
        return [item for item in value if isinstance(item, dict)] if isinstance(value, list) else []

    def endpoint(record: Mapping, role: str) -> str:
        keys = (("reported_gene", "reported_gene_id") if role == "reported" else
                ("reference_gene_name", "reference_gene", "reference_gene_id", "source_reference_gene_id"))
        names = list(dict.fromkeys(text(record.get(key)) for key in keys if text(record.get(key))))
        if not names:
            return ""
        species = text(record.get(role + "_species"))
        # A multi-species aggregate cannot assign a species to any one ID.
        species = (species if species and ";" not in species else
                   "種との個別対応は未記録" if species else "種は未記録")
        return " / ".join(names) + f"（{species}）"

    gene = _id(row)
    name = first(row, "gene_name", "candidate_gene_name", "focal_gene_name")
    gene_label = f"{name}（{gene}）" if name and name != gene else gene
    chrom, start, end = (text(row.get(key)) for key in ("chrom", "start", "end"))
    position = (f"{chrom or '未記録'}:{start or '未記録'}-{end or '未記録'}"
                if chrom or start or end else "未記録")
    build = first(row, "reference_build", "genome_assembly", "assembly")
    if build:
        position += f"（{build}）"

    claims = records(row.get("relevance_evidence_json"))
    homologs: list[str] = []
    for claim in claims:
        reference = endpoint(claim, "reference")
        bound_sources = source_tokens(claim.get("source_id") or claim.get("source_ids"))
        paired = False
        for paper in records(claim.get("paper_reports_json")):
            paper_sources = source_tokens(paper.get("source_id")) | source_tokens(paper.get("source_identifier_aliases"))
            if not bound_sources & paper_sources:
                continue
            reported = endpoint(paper, "reported")
            reference_id = text(paper.get("reference_gene_id"))
            if not reported:
                continue
            if not reference_id:
                homologs.append(f"報告: {reported}（参照との対応は未記録）")
                continue
            matched = gene_binding(claim, reference_id) == "matched"
            if not matched:
                homologs.append(f"報告: {reported}; 原報の対応先ID: {reference_id}（参照との対応は未記録）")
                continue
            # A claim can contain several references. The paper selects only
            # its own endpoint; displaying all claim IDs would invent edges.
            proxy = endpoint({"reference_gene_id": reference_id,
                              "reference_species": claim.get("reference_species", "")},
                             "reference")
            relation = text(paper.get("gene_correspondence"))
            if relation in {"homolog", "same_gene"}:
                homologs.append(f"報告: {reported} → 参照: {proxy}")
            else:
                homologs.append(f"報告: {reported}; 参照: {proxy}（対応: {relation or '未記録'}）")
            paired = True
        if reference and not paired:
            homologs.append("参照: " + reference)
    if not homologs:
        # The top-level fields are independent set joins in annotate_relevance.
        # Keep any saved IDs visible, without inventing a gene/species pairing.
        reference = first(row, "reference_gene", "reference_gene_id", "orthology_reference_gene")
        if reference:
            homologs.append(f"参照: {reference}（種との個別対応は未記録）")

    summaries = records(row.get("paper_evidence_summary_json"))
    # Structured single-paper IDs may themselves contain semicolons (DOIs).
    # Use the legacy set-valued column only when no paper summary is present.
    citations = [text(row.get("tier_source_id"))]
    citations += [text(paper.get("source_id")) for paper in summaries]
    if not summaries and text(row.get("source_ids")) != text(row.get("tier_source_id")):
        citations += sorted(source_tokens(row.get("source_ids")))
    links, seen = [], set()
    for citation in citations:
        kind, identifier = classify_citation(citation)
        if kind not in {"doi", "pmid", "pmcid"} or (kind, identifier) in seen:
            continue
        seen.add((kind, identifier))
        prefix = {"doi": "https://doi.org/", "pmid": "https://pubmed.ncbi.nlm.nih.gov/",
                  "pmcid": "https://pmc.ncbi.nlm.nih.gov/articles/"}[kind]
        url = prefix + quote(identifier, safe="/._-:") + ("/" if kind != "doi" else "")
        links.append(f"[{cell(citation)}]({url})")

    reason = first(row, "tier_reason", "initial_tier_reason", "reason", "relevance_basis")
    reason = reason or first(decision, "proposal_reason", "independent_reason", "reason") or "未記録"
    memo = reason if len(reason) <= 180 else reason[:179] + "…"
    copy_state = text(row.get("copy_state"))
    if any(token in copy_state for token in ("ambiguous", "co_preferred", "unresolved")):
        memo += f"; copy: {copy_state}"
    if disposition == "pending":
        memo = "保留中。" + memo
    elif disposition != "adopted":
        memo = "保持（採用状態未指定）。" + memo
    return [cell(gene_label), cell(position), cell("; ".join(dict.fromkeys(homologs)) or "未記録"),
            cell(row.get("relevance_tier", "") or "未記録"), "; ".join(links) or "未記録", cell(memo)]


def write_wet_review_report(
    wet_path: str,
    out_path: str,
    case_id: str,
    state: str,
    stage: str,
    *,
    adopted_gene_ids: Iterable[str] = (),
    pending_gene_ids: Iterable[str] = (),
    execution_context: Mapping[str, Any] | None = None,
    route_a_rethink: Mapping[str, Any] | None = None,
    review_history_state: str = "",
    selection_policy_version: str = "",
    review_decisions_path: str | Path | None = None,
    review_decisions_sha256: str | None = None,
    next_boundary_scope_manifest: Mapping[str, Any] | None = None,
    resume_scope_manifest: Mapping[str, Any] | None = None,
    provisional_wet_context: Mapping[str, Any] | None = None,
    interval_unresolved: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Write a diagnosis-oriented Markdown view of an adopted Wet artifact.

    The report is descriptive only: it records rows and review state and makes
    no accuracy or causal attribution claim.  Values are escaped as Markdown
    text so IDs/reasons cannot alter the table structure.
    """
    source = Path(wet_path).resolve()
    source_sha = _sha(source)
    columns, rows = _read(source)
    if _sha(source) != source_sha:
        raise ValueError("Wet artifact changed while building the review report")
    adopted, pending = set(adopted_gene_ids), set(pending_gene_ids)
    known = {_id(row) for row in rows}
    if not adopted <= known or not pending <= known:
        raise ValueError("review report IDs must belong to the Wet artifact")
    decision_rows: list[dict[str, str]] = []
    decision_path = None
    decision_sha = ""
    if review_decisions_path:
        if not review_decisions_sha256:
            raise ValueError("review decision artifact requires a sealed SHA-256")
        decision_path = Path(review_decisions_path).resolve()
        if not decision_path.is_file():
            raise ValueError("review decision artifact is missing")
        decision_sha = _sha(decision_path)
        if review_decisions_sha256 and decision_sha != str(review_decisions_sha256).removeprefix("sha256:"):
            raise ValueError("review decision artifact does not match the sealed selection hash")
        _, decision_rows = _read(decision_path)
        if _sha(decision_path) != decision_sha:
            raise ValueError("review decision artifact changed while reading")
    decision_by_id: dict[str, dict[str, str]] = {}
    for decision in decision_rows:
        gene = _id(decision)
        if gene and gene not in decision_by_id:
            decision_by_id[gene] = decision
    decision_counts: dict[str, int] = {}
    for decision in decision_rows:
        key = str(decision.get("applied_decision", "not_recorded")).strip() or "not_recorded"
        decision_counts[key] = decision_counts.get(key, 0) + 1
    report = Path(out_path).resolve()
    report.parent.mkdir(parents=True, exist_ok=True)

    def cell(value: object) -> str:
        return str(value if value is not None else "").replace("\\", "\\\\").replace("|", "\\|").replace("\n", " ")

    def joined(row: dict[str, str], *keys: str) -> str:
        values: list[str] = []
        for key in keys:
            value = str(row.get(key, "") or "").strip()
            if value and value not in values:
                values.append(value)
        return "; ".join(values)

    next_scope_path = ""
    if isinstance(next_boundary_scope_manifest, Mapping):
        raw_scope_path = str(next_boundary_scope_manifest.get("path", "") or "")
        if raw_scope_path:
            scope_path = Path(raw_scope_path).resolve()
            expected_scope_sha = str(next_boundary_scope_manifest.get("sha256", "") or "")
            if not expected_scope_sha:
                raise ValueError("next boundary scope manifest requires a sealed SHA-256")
            if not scope_path.is_file():
                raise ValueError("next boundary scope manifest is missing at report time")
            actual_scope_sha = _sha(scope_path)
            if expected_scope_sha and actual_scope_sha != expected_scope_sha.removeprefix("sha256:"):
                raise ValueError("next boundary scope manifest changed after selection")
            next_scope_path = str(scope_path)
    resume_scope_path = ""
    if isinstance(resume_scope_manifest, Mapping):
        raw_resume_path = str(resume_scope_manifest.get("path", "") or "")
        if raw_resume_path:
            resume_path = Path(raw_resume_path).resolve()
            expected_resume_sha = str(resume_scope_manifest.get("sha256", "") or "")
            if not expected_resume_sha or not resume_path.is_file() or _sha(resume_path) != expected_resume_sha.removeprefix("sha256:"):
                raise ValueError("resume scope manifest is missing or changed")
            resume_scope_path = str(resume_path)
    next_action = {
        "boundary_review_requested": "Complete the recorded boundary review using the recorded scope; request creation is not review execution.",
        "boundary_review_incomplete_retained": "Complete or explicitly review the unresolved boundary candidates before treating the set as final.",
        "next_boundary_review_required": "Complete the next recorded boundary review; higher-tier membership is already frozen.",
        "review_required_round_budget_exhausted": "The bounded review budget is exhausted; decide the next evidence or experiment step explicitly.",
        "completed_boundary_selection": "Inspect the adopted Wet set and plan experimental confirmation.",
        "completed_evidence_scoped_selection": "Inspect the adopted Wet set and plan experimental confirmation.",
        "provisional_adopted_wet": "Use the concrete Wet set for provisional experiments while investigating every unresolved interval candidate in the sealed ledger.",
    }.get(state, "Inspect the recorded review state before Wet adoption.")
    context_text = ""
    if execution_context:
        context_text = json.dumps(execution_context, ensure_ascii=False, sort_keys=True)
        if len(context_text) > 1800:
            context_text = context_text[:1799] + "…"
    lines = [
        "# Wet review report", "", f"- case_id: `{cell(case_id)}`",
        f"- state: `{cell(state)}`", f"- stage: `{cell(stage)}`",
        f"- candidate_count: {len(rows)}",
        f"- adopted_candidate_count: {len(adopted)}",
        f"- pending_candidate_count: {len(pending)}",
        f"- Wet SHA-256: `{cell(source_sha)}`",
    ]
    if route_a_rethink:
        lines.extend(["", "Route A rethink: " + json.dumps(dict(route_a_rethink), ensure_ascii=False, sort_keys=True), ""])
    if selection_policy_version:
        lines.insert(4, f"- selection_policy_version: `{cell(selection_policy_version)}`")
    if review_history_state:
        lines.append(f"- review_history_state: `{cell(review_history_state)}`")
    if decision_counts:
        lines.append("- decision counts: " + ", ".join(f"{key}={value}" for key, value in sorted(decision_counts.items())))
    if decision_path:
        lines.append(f"- review decision artifact: `{cell(str(decision_path))}` (SHA-256 `{cell(decision_sha)}`)")
    if next_scope_path:
        lines.append(f"- next boundary scope: `{cell(next_scope_path)}`")
    if resume_scope_path:
        lines.append(f"- retry from original boundary scope: `{cell(resume_scope_path)}`")
    if provisional_wet_context:
        for role in ("original_inventory", "screen_ledger"):
            entry = provisional_wet_context[role]
            _checked(entry)
            lines.append(f"- provisional {role}: `{cell(entry['path'])}` (SHA-256 `{cell(entry['sha256'])}`)")
        lines.append(f"- provisional deferral generations: {len(provisional_wet_context['deferral_provenance'])}")
    if interval_unresolved:
        _checked(interval_unresolved)
        lines.append(f"- interval unresolved: `{cell(interval_unresolved['path'])}` (SHA-256 `{cell(interval_unresolved['sha256'])}`); not_evaluated and deferred are not biological negatives")
    if context_text:
        lines.append(f"- execution context: `{cell(context_text)}`")
    if isinstance(execution_context, Mapping):
        route_b = execution_context.get("route_b_execution")
        if isinstance(route_b, Mapping):
            lines.append(
                "- Route B execution state: "
                f"`{cell(route_b.get('state', 'not_recorded'))}`; "
                f"completed=`{cell(route_b.get('completed_count', 0))}`, "
                f"not_evaluated=`{cell(route_b.get('not_evaluated_count', 0))}`"
            )
            if isinstance(route_b.get("response_state_counts"), Mapping):
                counts = route_b["response_state_counts"]
                lines.append("- Route B response states: " + ", ".join(f"{cell(key)}={cell(value)}" for key, value in counts.items()))
                lines.append("- Bの内訳は回答の受領・形式状態です。読了や文献探索の完了を示すものではありません。")
            b_artifact = str(route_b.get("path", "") or "")
            if b_artifact:
                b_path = Path(b_artifact)
                if not b_path.is_file():
                    raise ValueError("Route B execution manifest is missing at report time")
                b_digest = _sha(b_path)
                if route_b.get("sha256") and b_digest != str(route_b["sha256"]).removeprefix("sha256:"):
                    raise ValueError("Route B execution manifest changed after the first stage")
                lines.append(f"- Route B execution artifact: {cell(b_artifact)} (SHA-256 `{cell(b_digest)}`)")
        route_c = execution_context.get("route_c_execution")
        if isinstance(route_c, Mapping):
            lines.append(f"- Route C execution state: `{cell(route_c.get('state', 'not_evaluated_no_input'))}`")
            if route_c.get("reason"):
                lines.append(f"- Route C execution note: {cell(route_c['reason'])}")
            for label, path_key, digest_key in (
                ("status", "status_path", "status_sha256"),
                ("members", "members_path", "members_sha256"),
                ("producer manifest", "manifest_path", "manifest_sha256"),
            ):
                c_artifact = str(route_c.get(path_key, "") or "")
                if c_artifact:
                    c_path = Path(c_artifact)
                    if not c_path.is_file():
                        raise ValueError(f"Route C {label} artifact is missing at report time")
                    c_digest = _sha(c_path)
                    if route_c.get(digest_key) and c_digest != str(route_c[digest_key]).removeprefix("sha256:"):
                        raise ValueError(f"Route C {label} artifact changed after the first stage")
                    lines.append(f"- Route C {label} artifact: {cell(c_artifact)} (SHA-256 `{cell(c_digest)}`)")
    lines += [
        f"- next action: {cell(next_action)}", "",
        "This report describes the supplied artifact and review state; it does not assert causal accuracy or gene attribution.", "",
        "## 候補一覧", "",
        "入力Wetの行順・ID・Tierを保持した表示です。位置は保存値のまま、未記録は不存在を意味しません。"
        "論文リンクは保存済みIDから作成し、一言メモは元の理由の冒頭を表示します。全文と採用状態は下の詳細表を参照してください。", "",
        "| 遺伝子名・ID | 位置 | 相同遺伝子（種・ID） | Tier | 論文 | 一言メモ |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in rows:
        gene = _id(row)
        disposition = "adopted" if gene in adopted else "pending" if gene in pending else "retained"
        lines.append("| " + " | ".join(_wet_summary_cells(row, disposition, decision_by_id.get(gene, {}))) + " |")
    lines += [
        "", "## 判断・来歴の詳細", "",
        "| candidate_gene_id | relevance_tier | support_state | membership_basis | mapping_tier | copy_state | route_provenance | source | decision | reason | disposition |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in rows:
        gene = _id(row)
        decision = decision_by_id.get(gene, {})
        disposition = "adopted" if gene in adopted else "pending" if gene in pending else "retained"
        reason = joined(decision, "proposal_reason", "independent_reason", "reason")
        reason = reason or joined(row, "reason", "tier_reason", "relevance_basis")
        applied_decision = decision.get("applied_decision", "")
        source_text = joined(row, "source_ids", "review_source_ids", "review_verification_statuses")
        route = joined(row, "route_provenance", "route_support", "reference_gene")
        lines.append("| " + " | ".join(cell(value) for value in (
            gene,
            row.get("relevance_tier", ""),
            row.get("support_state", ""),
            row.get("membership_basis", ""),
            row.get("mapping_tier", ""),
            row.get("copy_state", ""),
            route,
            source_text,
            applied_decision or "not_recorded",
            reason,
            disposition,
        )) + " |")
    if _sha(source) != source_sha:
        raise ValueError("Wet artifact changed while writing the review report")
    report.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    result = {
        "path": str(report), "sha256": _sha(report), "candidate_count": len(rows),
        "adopted_candidate_count": len(adopted), "pending_candidate_count": len(pending),
        "wet_candidates_sha256": source_sha, "state": state, "stage": stage,
        "review_history_state": review_history_state, "decision_counts": decision_counts,
        "review_decisions": {"path": str(decision_path), "sha256": decision_sha} if decision_path else None,
        "next_boundary_scope_manifest": next_boundary_scope_manifest,
        "execution_context": execution_context, "next_action": next_action,
        "accuracy_asserted": False,
    }
    if provisional_wet_context:
        result.update(provisional_wet_context=provisional_wet_context, interval_unresolved=interval_unresolved)
    return result


def _primary_paper_evidence(rows: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Materialize the imported handoff, sharing each paper panel across copies.

    Summary markers identify the selected papers. Their actual observations
    come from the original bound reports, because a summary's strongest display
    representative can belong to another scope/reported gene in the same paper.
    This reads stamped ranks; it never grades, selects again, or fills slots.
    """
    from .paper_reports import PAPER_FIELD, PAPER_HANDOFF_POLICY

    def records(value: Any) -> list[dict]:
        try:
            value = json.loads(value) if isinstance(value, str) and value else value
        except (ValueError, TypeError):
            return []
        return [item for item in value if isinstance(item, dict)] if isinstance(value, list) else []

    def primary_key(marker: Mapping) -> tuple[str, int] | None:
        rank = str(marker.get("paper_handoff_rank", ""))
        scope = str(marker.get("paper_handoff_scope", "")).strip()
        if (marker.get("paper_handoff_policy") == PAPER_HANDOFF_POLICY
                and marker.get("paper_handoff_state") == "primary" and scope
                and rank.isdecimal() and 1 <= int(rank) <= 5):
            return scope, int(rank)
        return None

    panels: dict[str, dict] = {}
    for row in rows:
        gene = _id(row)
        selected = set()
        for summary in records(row.get("paper_evidence_summary_json")):
            for marker in records(summary.get("paper_handoff_panels")):
                key = primary_key(marker)
                if key is None:
                    continue
                scope, rank = key
                selected.add(key)
                panel = panels.setdefault(scope, {"candidate_gene_ids": set(), "papers": {}})
                panel["candidate_gene_ids"].add(gene)
                paper = panel["papers"].setdefault(rank, {"candidate_gene_ids": set(), "source_ids": set(), "reports": {}})
                paper["candidate_gene_ids"].add(gene)
                if summary.get("source_id"):
                    paper["source_ids"].add(str(summary["source_id"]))
        for claim in records(row.get("relevance_evidence_json")):
            for report in records(claim.get(PAPER_FIELD)):
                key = primary_key(report)
                if key not in selected:
                    continue
                scope, rank = key
                paper = panels[scope]["papers"][rank]
                if report.get("source_id"):
                    paper["source_ids"].add(str(report["source_id"]))
                encoded = json.dumps(report, ensure_ascii=False, sort_keys=True)
                bound = paper["reports"].setdefault(encoded, {"paper_report": report, "candidate_gene_ids": set()})
                bound["candidate_gene_ids"].add(gene)
    result = []
    for scope, panel in sorted(panels.items()):
        papers = []
        for rank, paper in sorted(panel["papers"].items()):
            reports = [{"paper_report": bound["paper_report"], "candidate_gene_ids": sorted(bound["candidate_gene_ids"])}
                       for _, bound in sorted(paper["reports"].items())]
            papers.append({"paper_handoff_rank": str(rank), "source_ids": sorted(paper["source_ids"]),
                           "candidate_gene_ids": sorted(paper["candidate_gene_ids"]), "reports": reports,
                           "record_state": "original_reports" if reports else "original_report_unavailable"})
        result.append({"paper_handoff_policy": PAPER_HANDOFF_POLICY, "paper_handoff_scope": scope,
                       "candidate_gene_ids": sorted(panel["candidate_gene_ids"]), "papers": papers})
    return result


def build_boundary_review_request(
    scope_manifest_arg: str,
    out_dir_arg: str,
    *,
    scope_sha256: str,
    case_id: str,
    adaptive_round: int = 0,
    max_adaptive_rounds: int = 1,
    comparison_bundle_arg: str | None = None,
) -> dict[str, Any]:
    """Create an immutable, boundary-only review request.

    The old adaptive request asks for ambiguity resolution and intentionally has
    no candidate selection schema.  A boundary review is a separate seam: the
    reviewer may propose adopt/defer/not_evaluated for members of the current
    tier, while higher tiers are supplied as read-only context.  This function
    only writes the request; it never executes a model.
    """
    if (isinstance(adaptive_round, bool) or isinstance(max_adaptive_rounds, bool)
            or not isinstance(adaptive_round, int) or not isinstance(max_adaptive_rounds, int)
            or adaptive_round < 0 or max_adaptive_rounds < 0):
        raise ValueError("adaptive rounds must be non-negative integers")
    manifest, columns, plan = _load_scope(scope_manifest_arg, scope_sha256, case_id)
    if not plan["boundary"]:
        terminal_state = (plan["state"] if len(plan["frozen"]) > manifest["upper"]
                          else "skipped_at_or_below_25")
        return {
            "contract": REQUEST_CONTRACT, "stage": "boundary_review_exhausted",
            "case_id": case_id, "adaptive_round": adaptive_round,
            "max_adaptive_rounds": max_adaptive_rounds, "state": terminal_state,
            "scope_manifest": _descriptor(Path(scope_manifest_arg)),
        }
    if adaptive_round >= max_adaptive_rounds:
        return {
            "contract": REQUEST_CONTRACT, "stage": "boundary_review_exhausted",
            "case_id": case_id, "adaptive_round": adaptive_round,
            "max_adaptive_rounds": max_adaptive_rounds, "state": "review_required_round_budget_exhausted",
            "scope_manifest": _descriptor(Path(scope_manifest_arg)),
            "boundary_tier": plan["boundary_tier"],
            "boundary_candidate_count": len(plan["boundary"]),
        }
    output = Path(out_dir_arg).resolve()
    output.mkdir(parents=True, exist_ok=False)
    boundary_path = Path(manifest["artifacts"]["boundary"]["path"])
    frozen_path = Path(manifest["artifacts"]["frozen"]["path"])
    policy_version = manifest.get("policy_version", LEGACY_POLICY)
    boundary_tiers = sorted({row["relevance_tier"] for row in plan["boundary"]},
                            key=RELEVANCE_TIERS.index)
    request = {
        "contract": REQUEST_CONTRACT, "schema_version": "1.0",
        "stage": "boundary_review_request", "state": "awaiting_review",
        "case_id": case_id, "adaptive_round": adaptive_round,
        "max_adaptive_rounds": max_adaptive_rounds,
        "scope_manifest": _descriptor(Path(scope_manifest_arg)),
        "source_wet_sha256": _sha(boundary_path),
        "policy_version": policy_version,
        "boundary_tier": ("evidence_scoped_mixed_tiers"
                          if policy_version in EVIDENCE_SCOPED_POLICIES else plan["boundary_tier"]),
        "boundary_tiers": boundary_tiers,
        "candidate_count": len(plan["boundary"]),
        "candidate_gene_ids": [_id(row) for row in plan["boundary"]],
        "candidate_input": _descriptor(boundary_path),
        "higher_tier_context": {
            "state": "read_only",
            "candidate_input": _descriptor(frozen_path),
            "candidate_gene_ids": [_id(row) for row in plan["frozen"]],
            "rows": plan["frozen"],
        },
        "proposal_schema": {
            "required": ["candidate_gene_id", "decision", "reason",
                          "source_wet_sha256", "executor_id", "executor_model"],
            "decisions": ["adopt", "defer", "not_evaluated"],
            "note": ("Each row must bind this boundary generation; only IDs in the sealed comparison scope are allowed."
                     if policy_version in EVIDENCE_SCOPED_POLICIES else
                     "Each row must bind this boundary generation; no higher/lower tier IDs."),
        },
        "reviewer_schema": {
            "required": ["candidate_gene_id", "decision", "reason",
                          "source_wet_sha256", "proposal_sha256", "reviewer_id", "reviewer_model"],
            "independent": True,
        },
        "policy": [
            ("compare only the first overflowing relevance tier; preserve higher tiers and keep later tiers in the ledger without another comparison or within-tier ranking"
             if policy_version == UPPER_EVIDENCE_SCOPED_POLICY else
             "compare every eligible candidate already present from the first overflowing tier downward; non-proposed identity-only inventory is excluded"
             if policy_version in EVIDENCE_SCOPED_POLICIES else
             "compare the whole boundary tier without within-tier ranking"),
            "preserve missing or disagreeing review as not_evaluated",
            "never modify higher-tier context or support state",
            "do not infer gene attribution from missing evidence",
            "adopt means retain the review candidate in Wet membership; it does not promote support, verify attribution, or choose a single copy",
            "incomplete support or unresolved copy alone is not a reason for not_evaluated or defer; retain all plausible copies when the review packet supports the relevance tier",
            ("use not_evaluated for unfinished original assessment; use defer only for the original actor's explicit reasoned experiment-priority decision"
             if policy_version == UPPER_EVIDENCE_SCOPED_POLICY else
             "use not_evaluated for missing/insufficient review material and unresolved disagreement; use defer only for an explicit agreed experiment-priority decision"),
            "compare explanations against the whole recorded phenotype; separate supplied observations from additional assumptions about tissue, timing, regulation or perturbation direction",
            "inspect reference-specific function, aligned regions and available domain annotations alongside alternative protein identities; shared family similarity alone does not verify transfer of a particular function",
            "a plausible mechanism is a review lead, not by itself an experiment-priority comparison; explain adopt/defer relative to the other mechanisms in the whole tier while retaining indistinguishable plausible copies and regulatory or structural routes",
        ],
        "rank_contribution": "none",
    }
    request_path = output / "adaptive_synthesis_request.json"
    if comparison_bundle_arg:
        from .boundary_comparison_material import build_comparison_context
        context = build_comparison_context(
            comparison_bundle_arg, case_id,
            request["candidate_gene_ids"] + request["higher_tier_context"]["candidate_gene_ids"],
        )
        context_path = output / "comparison_context.json"
        context_path.write_text(json.dumps(context, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        request["comparison_context"] = _descriptor(context_path)
    if policy_version == UPPER_EVIDENCE_SCOPED_POLICY:
        request["assessment_mode"] = "original_assessment_v1"
        request["reviewer_schema"] = {"independent": False, "required": [], "state": "not_required_original_assessment"}
        request["policy"].append(
            "Compare experiment priority only within the sealed active scope using existing recorded evidence. "
            "Reuse the original Tier and phenotype observations; do not add a Tier/content reviewer or a separate "
            "phenotype reassessment. Lower original eligibility remains in the ledger, never padding or a next round. "
            "Twenty-five is not a truncation cap and does not require additional searches or forced reduction.")
        request["policy"].append(
            "Keep the original Tier unchanged and compare source-to-copy correspondence separately. "
            "Use the supplied function, aligned regions, available domains, alternative identities and known conflicts "
            "to explain differences between candidates; repeated annotations or A/B agreement are not independent support. "
            "Defer is a reasoned decision to postpone this experiment, not a claim of noncausality: it does not require "
            "proof that the gene cannot be causal. State the comparative reason and what unresolved question would "
            "change that decision. Do not automatically defer partial alignments, ambiguous copies or missing material, "
            "and do not force the count to 20 or 25. If assessment is unfinished, use not_evaluated rather than a "
            "biological negative. Use the same original actor, without an added Tier reviewer or new search quota.")
    primary_papers = _primary_paper_evidence(plan["boundary"])
    higher_primary_papers = _primary_paper_evidence(plan["frozen"])
    if primary_papers:
        request["primary_paper_evidence"] = primary_papers
    if higher_primary_papers:
        request["higher_tier_context"]["primary_paper_evidence"] = higher_primary_papers
    if primary_papers or higher_primary_papers:
        request["policy"].append(
            "Read primary_paper_evidence first: it contains up to five already selected strongest distinct papers per original "
            "reference/shared-mechanism scope, chosen after the original per-paper Tier assessment. Reuse each scope across "
            "its listed copies; reports retain their own candidate bindings, reported genes, reference endpoints and checklist. "
            "This is the evidence handoff, not a new metadata budget or another paper selection. Do not refill, regrade or "
            "require five papers. Candidate-bound reports remain in paper_evidence_summary_json and "
            "relevance_evidence_json in the sealed candidate inputs. Unselected source records remain in the upstream "
            "import/freeze audit artifacts; do not expand candidate applicability from those records. Read audit-only records for a concrete question and "
            "retain known counterevidence. original_report_unavailable is missing material, never negative evidence. "
            "Several correctly assessed strong papers reduce reliance on one mistaken report; sorting cannot correct "
            "the original mistaken Tier or establish independent replication. Candidate membership, strongest applicable "
            "Tier and support policies are unchanged; higher-tier context remains read-only.")
    if any(row.get("tier_aggregation") == "strongest_paper_v1" for row in plan["boundary"] + plan["frozen"]):
        request["policy"].append(
            "For strongest_paper_v1, Tier belongs to each paper's report and the gene uses its strongest applicable paper. "
            "Paper counts, repeated A/B citations and weaker papers add no Tier or priority bonus. "
            "One T1 plus many T3 papers has the same Tier as one T1 alone. "
            "For paper_checklist_v1, reuse the Tier derived from the original A/B yes/no/unknown checklist; "
            "reported genes, proxy references and focal loci are separate roles. Metadata checks are audit-only; "
            "name/display differences or missing aliases never authorize a candidate veto. "
            "For observation_rules_v1, Tier is derived from original evidence_types, reported observations and comparison; "
            "For reported_paper_judgment_v1, reuse the original LLM's paper Tier/reason and semantic gene correspondence; "
            "the supplied reference_gene_id is an existing anchor, not a name-dictionary match or proof of verified support. "
            "legacy records retain their declared-Tier policy. Do not ask an agent to regrade Tier or fill missing fields. "
            "Read tier_missing_fields_json with the established interpretation; unknown is not negative evidence. "
            "Original studies and reviews are both allowed, including one T1 original plus three T1 reviews. "
            "No original-paper requirement, fixed paper quota, mandatory corroboration search, or forced T2 search after T1. "
            "Keep the existing T1 retention policy without extra count/type/independence conditions. " +
            ("Read the selected primary_paper_evidence alongside paper_tier_counts_json; all summaries remain audit context. "
             if primary_papers or higher_primary_papers else
             "Read paper_evidence_summary_json across Tiers and paper_tier_counts_json alongside the strongest paper. ") +
            "Counts are traceability metadata only: no score, priority bonus or independent-replication claim. "
            "Known review-to-study links may be recorded; unknown links do not trigger independence checks. "
            "Keep T3/lower literature when available; absent literature is explicit, not fabricated or an automatic veto. "
            "Open additional sources only for a concrete unresolved question; stop merely repetitive searches. "
            "Full records stay in the audit ledger for reading when needed, without an automatic extra review stage. "
            "Retain already known counterevidence and evaluate other genes/copies independently.")
    request_path.write_text(json.dumps(request, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    request_manifest = {
        "contract": REQUEST_CONTRACT, "stage": "boundary_review_request",
        "case_id": case_id, "scope_manifest_sha256": _sha(Path(scope_manifest_arg)),
        "source_wet_sha256": _sha(boundary_path), "candidate_count": len(plan["boundary"]),
        "policy_version": policy_version,
        "boundary_tier": request["boundary_tier"], "boundary_tiers": boundary_tiers,
        "adaptive_round": adaptive_round,
        "max_adaptive_rounds": max_adaptive_rounds, "output_sha256": _sha(request_path),
        "rank_contribution": "none",
    }
    manifest_path = output / "boundary_review_request_manifest.json"
    manifest_path.write_text(json.dumps(request_manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"request": str(request_path), "manifest": str(manifest_path), "state": "requested",
            "boundary_tier": plan["boundary_tier"], "candidate_count": len(plan["boundary"])}


def freeze_boundary_scope(
    inventory_arg: str, out_dir_arg: str, *, inventory_sha256: str, case_id: str,
    target: int = 20, upper: int = 25, frozen_gene_ids: Iterable[str] = (),
    protected_gene_ids: Iterable[str] = (),
    policy_version: str = LEGACY_POLICY, eligible_gene_ids: Iterable[str] = (),
) -> dict[str, Any]:
    """Snapshot a validated inventory into immutable, hash-bound review scopes.

    This is not a replacement canonical Wet producer. The inventory is not
    rewritten; the sealed policy determines which rows are frozen, compared, or
    left outside the active scope. Existing paths are never overwritten.
    """
    source = Path(inventory_arg).resolve()
    if not source.is_file() or _sha(source) != inventory_sha256:
        raise ValueError("inventory hash mismatch or missing")
    if not case_id.strip():
        raise ValueError("boundary scope requires case_id")
    columns, rows = _read(source)
    prefix = sorted(set(frozen_gene_ids))
    protected = sorted(set(protected_gene_ids))
    eligible = sorted(set(eligible_gene_ids))
    plan = partition_relevance_tiers(rows, target=target, upper=upper,
                                     frozen_gene_ids=prefix, protected_gene_ids=protected,
                                     policy_version=policy_version, eligible_gene_ids=eligible)
    output = Path(out_dir_arg).resolve()
    output.mkdir(parents=True, exist_ok=False)
    artifacts = {
        role: _write(output / f"{role}_tier_candidates.tsv", columns, plan[role])
        for role in ("frozen", "boundary", "lower")
    }
    manifest = {
        "contract": CONTRACT, "case_id": case_id, "stage": "boundary_scope_freeze",
        "inventory": _descriptor(source), "target": target, "upper": upper,
        "input_frozen_gene_ids": prefix, "protected_gene_ids": protected,
        "policy_version": policy_version,
        "eligible_gene_ids": plan["eligible_gene_ids"],
        **({"active_gene_ids": plan["active_gene_ids"], "upper_priority_applied": plan["upper_priority_applied"]}
           if policy_version == UPPER_EVIDENCE_SCOPED_POLICY else {}),
        "boundary_tier": plan["boundary_tier"],
        "state": plan["state"], "artifacts": artifacts,
        "counts": {role: len(plan[role]) for role in artifacts},
        "rank_contribution": "none",
    }
    path = output / "boundary_scope_manifest.json"
    path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {**manifest, "manifest": str(path), "manifest_sha256": _sha(path)}


def _load_scope(path_arg: str, expected_sha: str, case_id: str) -> tuple[dict[str, Any], list[str], dict[str, Any]]:
    path = Path(path_arg).resolve()
    if not path.is_file() or _sha(path) != expected_sha:
        raise ValueError("boundary scope manifest hash mismatch")
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("contract") != CONTRACT or manifest.get("stage") != "boundary_scope_freeze" or manifest.get("case_id") != case_id:
        raise ValueError("boundary scope contract/stage/case mismatch")
    source = _checked(manifest["inventory"])
    columns, rows = _read(source)
    policy_version = manifest.get("policy_version", LEGACY_POLICY)
    plan = partition_relevance_tiers(rows, target=manifest["target"], upper=manifest["upper"],
                                     frozen_gene_ids=manifest["input_frozen_gene_ids"],
                                     protected_gene_ids=manifest.get("protected_gene_ids", []),
                                     policy_version=policy_version,
                                     eligible_gene_ids=manifest.get("eligible_gene_ids", []))
    if plan["boundary_tier"] != manifest["boundary_tier"] or plan["state"] != manifest["state"]:
        raise ValueError("boundary scope no longer matches its inventory")
    if policy_version == UPPER_EVIDENCE_SCOPED_POLICY and any(
        manifest.get(key) != plan[key] for key in ("active_gene_ids", "upper_priority_applied")
    ):
        raise ValueError("active scope differs from original eligibility and inventory")
    for role in ("frozen", "boundary", "lower"):
        table = _checked(manifest["artifacts"][role])
        actual_columns, actual_rows = _read(table)
        if actual_columns != columns or actual_rows != plan[role] or manifest["counts"][role] != len(actual_rows):
            raise ValueError("boundary scope rows differ from the original inventory")
    return manifest, columns, plan


def consume_boundary_selection(
    scope_manifest_arg: str, review_arg: str, out_dir_arg: str, *,
    scope_sha256: str, case_id: str, reviewer_id: str = "", reviewer_model: str = "",
    executor_id: str, executor_model: str, independent_review_arg: str | None = None,
    allow_provisional_wet: bool = False,
    provisional_screen_ledger_arg: str | None = None,
    provisional_screen_ledger_sha256: str | None = None,
    provisional_wet_context: Mapping[str, Any] | None = None,
    prior_deferral_ledgers: Iterable[Mapping[str, Any]] = (),
    review_history_state: str = "complete",
) -> dict[str, Any]:
    """Delegate scoped decisions to the existing policy-bound selection consumer.

    Rows outside the sealed comparison scope cannot be mentioned in the review.
    Missing/disagreeing review keeps the scoped member under the existing
    contract. Legacy scopes may continue through lower tiers; evidence-scoped
    scopes complete with the already-reviewed eligible pool and never add an
    unproposed lower roster merely because the adopted set is under target.
    """
    manifest, columns, plan = _load_scope(scope_manifest_arg, scope_sha256, case_id)
    if manifest.get("policy_version") == UPPER_EVIDENCE_SCOPED_POLICY and (allow_provisional_wet or provisional_wet_context):
        raise ValueError("legacy provisional adoption cannot override the sealed upper evidence scope")
    if (provisional_screen_ledger_arg or provisional_screen_ledger_sha256) and not allow_provisional_wet:
        raise ValueError("provisional screen ledger requires allow_provisional_wet")
    if not plan["boundary"]:
        raise ValueError("no boundary tier requires review")
    boundary_ids = {_id(row) for row in plan["boundary"]}
    for path_arg in (review_arg, independent_review_arg):
        if path_arg:
            _, rows = _read(Path(path_arg))
            if any(_id(row) not in boundary_ids for row in rows):
                raise ValueError("review attempted to change a higher/lower-tier candidate")
    # Local import avoids a module cycle and keeps pure scope tests provider-free.
    from .adaptive_synthesis import consume_review_inclusive_selection

    output = Path(out_dir_arg).resolve()
    output.mkdir(parents=True, exist_ok=False)
    selected = consume_review_inclusive_selection(
        manifest["artifacts"]["boundary"]["path"], review_arg, str(output / "boundary_selection"),
        reviewer_id=reviewer_id, reviewer_model=reviewer_model,
        executor_id=executor_id, executor_model=executor_model,
        independent_review_arg=independent_review_arg, case_id=case_id,
        first_stage_sha256=manifest["inventory"]["sha256"],
        **({"original_assessment_scope_arg": scope_manifest_arg, "original_assessment_scope_sha256": scope_sha256}
           if manifest.get("policy_version") == UPPER_EVIDENCE_SCOPED_POLICY else {}),
    )
    selection_manifest = Path(str(selected["manifest"]))
    decision = json.loads(selection_manifest.read_text(encoding="utf-8"))
    chosen_path = Path(str(selected["wet_candidates"]))
    if decision.get("source_wet_sha256") != manifest["artifacts"]["boundary"]["sha256"] or decision.get("wet_candidates_sha256") != _sha(chosen_path):
        raise ValueError("boundary selection does not bind the reviewed generation")
    chosen_columns, chosen = _read(chosen_path)
    by_id = {_id(row): row for row in plan["boundary"]}
    if chosen_columns != columns or len({_id(row) for row in chosen}) != len(chosen) or any(by_id.get(_id(row)) != row for row in chosen):
        raise ValueError("boundary selection changed a candidate row or introduced a member")
    # The existing consumer reports ``completed`` when proposal and
    # independent review agree, including agreement on an explicit
    # ``not_evaluated`` decision.  That is still an incomplete biological
    # review state and must remain visible across generations.
    _, decision_rows = _read(Path(str(selected["review_decisions"])))
    explicit_unresolved = any(str(row.get("applied_decision", "")).strip() == "not_evaluated"
                              for row in decision_rows)
    # Recheck all snapshots after the consumer, including the original source.
    _load_scope(scope_manifest_arg, scope_sha256, case_id)
    prefix = plan["frozen"] + chosen
    remaining = prefix + plan["lower"]
    policy_version = manifest.get("policy_version", LEGACY_POLICY)
    continued = partition_relevance_tiers(
        remaining, target=plan["target"], upper=plan["upper"],
        frozen_gene_ids=(_id(row) for row in prefix),
        protected_gene_ids=manifest.get("protected_gene_ids", []),
        policy_version=policy_version,
        eligible_gene_ids=(manifest.get("eligible_gene_ids", []) if policy_version == LEGACY_POLICY else
                           set(manifest.get("eligible_gene_ids", [])) & {_id(row) for row in remaining}),
        active_gene_ids=(set(manifest["active_gene_ids"]) & {_id(row) for row in remaining}
                         if policy_version == UPPER_EVIDENCE_SCOPED_POLICY else None),
    )
    provisional = None
    context = None
    unresolved = None
    if allow_provisional_wet or provisional_wet_context:
        from .provisional_wet import assess_provisional_wet

        prior_ledgers = list(prior_deferral_ledgers)
        previous = dict(provisional_wet_context or {})
        provenance = list(previous.get("deferral_provenance", []))
        if [entry["deferral_ledger"] for entry in provenance] != prior_ledgers:
            raise ValueError("provisional provenance must bind every prior deferral ledger")
        original_entry = previous.get("original_inventory", manifest["inventory"])
        original_path = _checked(original_entry)
        original_columns, original_rows = _read(original_path)
        if original_columns != columns:
            raise ValueError("provisional original inventory columns changed")
        if previous and previous.get("case_id") != case_id:
            raise ValueError("provisional context case mismatch")
        screen_entry = previous.get("screen_ledger")
        if provisional_screen_ledger_arg or provisional_screen_ledger_sha256:
            if not provisional_screen_ledger_arg or not provisional_screen_ledger_sha256:
                raise ValueError("provisional screen ledger requires path and SHA-256")
            supplied = {"path": str(Path(provisional_screen_ledger_arg).resolve()),
                        "sha256": provisional_screen_ledger_sha256}
            _checked(supplied)
            if screen_entry and supplied["sha256"] != screen_entry["sha256"]:
                raise ValueError("provisional screen ledger changed between generations")
            screen_entry = supplied
        if not screen_entry:
            raise ValueError("provisional Wet requires a hash-bound full screen ledger")
        screen_path = _checked(screen_entry)
        _, screen_rows = _read(screen_path)
        provenance.append({"deferral_ledger": _descriptor(Path(str(selected["deferral_ledger"]))),
                           "review_decisions": _descriptor(Path(str(selected["review_decisions"]))),
                           "selection_manifest": _descriptor(selection_manifest)})
        deferred = []
        for entry in provenance:
            ledger_path = _checked(entry["deferral_ledger"])
            decisions_path = _checked(entry["review_decisions"])
            selection_path = _checked(entry["selection_manifest"])
            bound = json.loads(selection_path.read_text(encoding="utf-8"))
            if (bound.get("case_id") != case_id
                    or bound.get("deferral_ledger_sha256") != entry["deferral_ledger"]["sha256"]
                    or bound.get("review_decisions_sha256") != entry["review_decisions"]["sha256"]):
                raise ValueError("provisional deferrals are not bound to their selection")
            for path_key, sha_key in (("proposal_artifact", "input_review_sha256"),
                                      ("independent_review_artifact", "independent_review_sha256")):
                if bound.get(path_key):
                    _checked({"path": bound[path_key], "sha256": bound.get(sha_key)})
            _, ledger_rows = _read(ledger_path)
            _, applied_rows = _read(decisions_path)
            agreed = {_id(row) for row in applied_rows if row.get("resolution_state") == "agreement"
                      and row.get("proposal_decision") in {"defer", "deferred"}
                      and row.get("independent_decision") == row.get("proposal_decision")
                      and row.get("applied_decision") == row.get("proposal_decision")
                      and row.get("proposal_reason", "").strip() and row.get("independent_reason", "").strip()}
            ledger_ids = [_id(row) for row in ledger_rows]
            decisions_by_id = {_id(row): row for row in applied_rows}
            if (set(ledger_ids) != agreed or len(ledger_ids) != len(set(ledger_ids))
                    or set(ledger_ids) != set(bound.get("deferred_candidate_genes", []))
                    or any(row.get("case_id") != case_id or row.get("disposition") != "deferred"
                           or not row.get("reason", "").strip() or not row.get("independent_reason", "").strip()
                           or row.get("reason") != decisions_by_id.get(_id(row), {}).get("proposal_reason")
                           or row.get("independent_reason") != decisions_by_id.get(_id(row), {}).get("independent_reason")
                           for row in ledger_rows)):
                raise ValueError("provisional deferral ledger must contain exactly independently agreed deferrals")
            if ledger_ids and not bound.get("independent_review_artifact"):
                raise ValueError("provisional deferrals require their sealed independent review")
            deferred.extend(ledger_ids)
        preserved = set(previous.get("original_frozen_gene_ids", plan["frozen_gene_ids"]))
        protected = set(previous.get("protected_gene_ids", manifest.get("protected_gene_ids", [])))
        if protected != set(manifest.get("protected_gene_ids", [])) or (preserved | protected) & set(deferred):
            raise ValueError("provisional deferrals cannot remove frozen or protected candidates")
        # Keep the complete source bytes and original unit roster; never project
        # the screen ledger down to the active candidates.
        original_copy = output / "original_interval_inventory.tsv"
        screen_copy = output / "provisional_screen_ledger.tsv"
        shutil.copyfile(original_path, original_copy)
        shutil.copyfile(screen_path, screen_copy)
        if _sha(original_copy) != original_entry["sha256"] or _sha(screen_copy) != screen_entry["sha256"]:
            raise ValueError("provisional source changed while snapshotting")
        context = {"case_id": case_id, "original_inventory": _descriptor(original_copy),
                   "screen_ledger": _descriptor(screen_copy), "deferral_provenance": provenance,
                   "original_frozen_gene_ids": sorted(preserved), "protected_gene_ids": sorted(protected)}
        # Account the original-minus-active roster even if the current gate is
        # ineligible. The assessor shares all first-stage eligibility rules.
        provisional = assess_provisional_wet(
            remaining, continued, screen_rows, case_id=case_id,
            protected_gene_ids=protected | {_id(row) for row in prefix}, upper=plan["upper"],
            coverage_inventory=original_rows, agreed_deferred_gene_ids=deferred,
        )
        if (not allow_provisional_wet or review_history_state == "incomplete"
                or selected["state"] != "completed" or explicit_unresolved):
            provisional = {"eligible": False, "reason": "review_history_incomplete" if
                           review_history_state == "incomplete" or selected["state"] != "completed" or explicit_unresolved
                           else "provisional_not_requested"}
        if provisional["eligible"]:
            unresolved_ids = set(provisional["unresolved_gene_ids"])
            unresolved_columns = [*columns, "provisional_disposition", "material_observation_state", "biological_exclusion"]
            unresolved = _write(output / "interval_unresolved.tsv", unresolved_columns, [
                {**row, "provisional_disposition": "deferred_not_negative" if _id(row) in set(deferred)
                 else "not_evaluated_identity_only", "material_observation_state": provisional["material_by_gene"][_id(row)],
                 "biological_exclusion": "not_established"}
                for row in original_rows if _id(row) in unresolved_ids
            ])
    # A pending next boundary is still retained as first-stage material. It is
    # never advertised as a completed ~20-gene final answer.
    wet_rows = (plan["frozen"] + chosen if policy_version in EVIDENCE_SCOPED_POLICIES
                else continued["frozen"] + continued["boundary"])
    if provisional and provisional["eligible"]:
        wet_rows = [row for row in _validated(remaining) if _id(row) in set(provisional["adopted_gene_ids"])]
    wet = _write(output / "wet_candidates.tsv", columns, wet_rows)
    continuation_inventory = _write(output / "continuation_inventory.tsv", columns, _validated(remaining))
    continuation = None if provisional and provisional["eligible"] or policy_version in EVIDENCE_SCOPED_POLICIES else freeze_boundary_scope(
        continuation_inventory["path"], str(output / "next_scope"),
        inventory_sha256=continuation_inventory["sha256"], case_id=case_id,
        target=plan["target"], upper=plan["upper"], frozen_gene_ids=(_id(row) for row in prefix),
        protected_gene_ids=manifest.get("protected_gene_ids", []),
    )
    state = ("boundary_review_incomplete_retained" if selected["state"] != "completed" or explicit_unresolved else
             "provisional_adopted_wet" if provisional and provisional["eligible"] else
             "completed_evidence_scoped_selection" if policy_version in EVIDENCE_SCOPED_POLICIES else
             "next_boundary_review_required" if continued["boundary"] else continued["state"])
    result = {
        "contract": CONTRACT, "stage": "boundary_selection_merge", "case_id": case_id,
        "policy_version": policy_version,
        "state": state, "inventory": manifest["inventory"],
        **({"assessment_mode": selected["assessment_mode"], "review_history_basis": selected["assessment_mode"]}
           if selected.get("assessment_mode") else {}),
        **({"eligible_gene_ids": manifest["eligible_gene_ids"], "active_gene_ids": manifest["active_gene_ids"],
            "upper_priority_applied": manifest["upper_priority_applied"]}
           if policy_version == UPPER_EVIDENCE_SCOPED_POLICY else {}),
        "selection_state": selected["state"],
        "review_history_state": ("incomplete" if explicit_unresolved or selected["state"] in {
            "not_evaluated_no_independent_review", "review_unresolved_disagreement",
            "completed_partial_review", "completed_partial_assessment"} else "complete"),
        "scope_manifest": _descriptor(Path(scope_manifest_arg)),
        "resume_scope_manifest": (_descriptor(Path(scope_manifest_arg))
                                  if policy_version in EVIDENCE_SCOPED_POLICIES
                                  and (selected["state"] != "completed" or explicit_unresolved) else None),
        "boundary_selection_manifest": _descriptor(selection_manifest),
        "higher_tier_rows_preserved": True, "frozen_higher_gene_ids": plan["frozen_gene_ids"],
        "output_candidate_count": len(wet_rows), "wet_candidates": wet,
        "next_scope_manifest": {"path": continuation["manifest"], "sha256": continuation["manifest_sha256"]} if continuation else None,
        "next_boundary_tier": continued["boundary_tier"],
        "deferral_ledger": _descriptor(Path(str(selected["deferral_ledger"]))),
        "review_decisions": _descriptor(Path(str(selected["review_decisions"]))),
        "support_state": "preserved", "rank_contribution": "none",
    }
    if context:
        result.update(provisional_wet_context=context, provisional_wet=provisional, interval_unresolved=unresolved)
    result_path = output / "boundary_selection_manifest.json"
    result_path.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {**result, "manifest": str(result_path)}
