"""Paper-bound observations; no evidence-count score or support verification."""
from __future__ import annotations

import hashlib
import json
import math
import re
from collections.abc import Mapping
from .observation_tiers import CHECKLIST_FIELDS, COMPARISONS, GENE_RELATIONS, JUDGMENT_POLICY, OBSERVATION_FIELDS, POLICY, gene_binding, has_observations, normalize_evidence_types
from .observation_tiers import CHECKLIST_POLICY, checklist_tier, normalize_checklist
from .observation_tiers import display_text, is_unknown, unfilled_answer
from .observation_tiers import CHECKLIST_CONSISTENCY_FIELDS, checklist_consistency

PAPER_FIELD = "paper_reports_json"
PAPER_POLICY = "strongest_paper_v1"
MAX_BIBLIOGRAPHY_PAPERS = 5
PAPER_TABLE_POLICY = "paper_evidence_table_v1"
PAPER_HANDOFF_POLICY = "tier_ranked_primary_papers_v1"
PAPER_HANDOFF_FIELDS = ("paper_handoff_policy", "paper_handoff_scope", "paper_handoff_state", "paper_handoff_rank")
RESEARCH_ANNOTATIONS = "research_annotations"
# Compatibility for research ledgers copied into a paper report. These names
# describe the investigation, never canonical paper facts or verification.
LEGACY_RESEARCH_FIELDS = frozenset({
    "unit_id", "counterevidence", "stage", "read_scope", "conditional_focal_reasoning",
    "source_url", "reference_species", "formal_claim_state", "retrieved_at",
    "full_paper_read", "paper_checklist_reasons",
})
PAPER_TABLE_FIELDS = ("paper_table_policy", "source_excerpt", "source_location", "unknown_reason", "followup_question", "followup_state",
                      "paper_checklist_missing_fields_json")
LITERATURE_RESEARCH_POLICY = {
    "candidate_paper_target": 20,
    "evidence_paper_target": MAX_BIBLIOGRAPHY_PAPERS,
    "unit": "original_reference_gene_or_shared_mechanism",
    "priority": "T1 T2 T3 T4 T5 T6 U",
    "publication_types": ["original", "review"],
    "shortfall_policy": "record screened sources, queries and retrieval failures; never fabricate or veto candidates",
    "membership_effect": "none",
    "paper_table_policy": PAPER_TABLE_POLICY,
    "followup_budget": "one targeted batch to the original researcher; never a membership gate",
}
LITERATURE_SEARCH_GUIDANCE = """
During the existing mechanism/reference investigation, check direct focal-species
phenotype evidence as well as functional evidence in other plants. Use the
phenotype and cellular/developmental synonyms, the focal species name, and any
source-grounded gene names/synonyms encountered; record queries, sources checked,
coverage and access failures. Broaden the mechanistic search when direct reports
are absent or inaccessible. This is not a Tier-only search or a T1 quota, and a
focal-species paper is not required to propose a mechanism. Do not use benchmark
answer hints, case-specific rescue rules, truth or other routes' results before freeze.
Respect the assigned route's interval visibility; A remains interval-blind, while B
uses its own frozen interval units.
Keep the experimentally reported gene/species, the sequence/reference proxy and
the focal locus distinct. Record observed facts separately from proposed transfer;
homology alone does not establish that two records are the same focal locus.
"""

REPORT_FIELDS = (
    "initial_tier", "initial_tier_reason", "reported_gene", "reported_gene_id", "reported_species",
    "reported_phenotype", "reported_process", "phenotype_comparison",
    "reported_organ", "reported_stage", "source_title", "source_identifier_aliases",
    "publication_type", "underlying_source_ids", "corroboration_note",
    "reference_gene_id", "gene_correspondence", "paper_checklist",
    *PAPER_TABLE_FIELDS,
    *OBSERVATION_FIELDS,
)


def initial_paper_report_example() -> dict:
    """One concrete per-record ``paper_reports_json`` element for the A/B prompts.

    This is a shared prompt example, never an admission requirement. It uses the
    canonical report keys and places every documented ``paper_checklist`` question
    explicitly with the ``unknown`` default, so a model cannot turn a missing
    answer into ``no``. ``reference_gene_id`` selects an existing reference
    endpoint; ``reported_gene``/``reported_gene_id`` stay the experimentally
    reported gene. The values are recognizable placeholders: the model must
    replace them with the real paper's observations, or drop the report.
    """
    return {
        "paper_table_policy": PAPER_TABLE_POLICY,
        "source_id": "<PMID/PMCID/DOI of the paper you actually read, or drop this report>",
        "source_title": "<the paper title>",
        "publication_type": "<original|review|unknown, when available>",
        "reported_gene": "<experimental gene name/ID as reported in that paper>",
        "reported_gene_id": "<the paper's own gene ID in its namespace, or unknown with a checked reason>",
        "reported_species": "<species of the reported experiment>",
        "reference_gene_id": "<selected reference/bait ID; B: exact claim.reference_gene>",
        "gene_correspondence": "<same_gene|homolog|different_gene|unknown>",
        "initial_tier_reason": "<one short reason tying the source facts to the reported gene>",
        "source_excerpt": "<short factual observation from the source when available, or empty>",
        "source_location": "<abstract, section, figure or page actually inspected when available, or empty>",
        "unknown_reason": "<what was checked and why any explicit unknown remains, or empty when none>",
        "followup_question": "<specific unresolved upper-Tier lead from this source, or empty>",
        "followup_state": "not_attempted",
        "paper_checklist": {field: "unknown" for field in CHECKLIST_FIELDS},
        RESEARCH_ANNOTATIONS: {
            "read_scope": "<source sections actually inspected>",
            "paper_checklist_reasons": {"phenotype_reported": "<short source-based reason, or checked unknown>"},
            "counterevidence": "<limitations or contrary observations, or empty>",
        },
    }


def render_initial_paper_report_example() -> str:
    """Compact JSON text of the one-element ``paper_reports_json`` example."""
    return json.dumps(
        [initial_paper_report_example()], ensure_ascii=False, separators=(",", ":")
    )


def source_tokens(value: object) -> set[str]:
    from .citation_tokens import split_citations
    from .source_verify import classify_citation
    values = split_citations(value)
    result = set()
    for value in values:
        if value is None or not str(value).strip():
            continue
        kind, identifier = classify_citation(str(value))
        result.add(f"{kind}:{identifier}".casefold() if kind in {"pmid", "pmcid", "doi"} else str(value).strip().casefold())
    return result


def normalize_source_aliases(value: object) -> str:
    """Accept JSON arrays and TSV text without changing identifier contents."""
    if value is None:
        return ""
    if isinstance(value, list):
        if any(not isinstance(item, str) for item in value):
            raise ValueError("source identifier aliases must be strings")
        return ";".join(item.strip() for item in value if item.strip())
    if not isinstance(value, str):
        raise ValueError("source identifier aliases must be text or an array")
    return value.strip()


def _ranked_bibliography_groups(row: Mapping, citations: list[str], *, recorded: bool = False) -> list[dict]:
    """Group aliases using either original Tiers or the recorded handoff order."""
    from .citation_tokens import split_citations
    from .source_verify import classify_citation
    reports = json.loads(normalize_paper_reports(row.get(PAPER_FIELD)) or "[]")
    declared = source_tokens(citations)
    groups = []
    for position, report in enumerate(reports):
        ids = [report["source_id"], *split_citations(report["source_identifier_aliases"])]
        if not source_tokens(ids) & declared:
            continue
        if recorded:
            # A downstream reader must not choose again from current report
            # contents. The importer already assigned each paper its position.
            rank = int(report.get("paper_handoff_rank") or 1000000)
        else:
            tier = report.get("initial_tier", "U").split("_", 1)[0]
            if not tier and has_observations(report):
                from .observation_tiers import classify_observation
                tier = classify_observation({**row, **report, "source_id": report["source_id"],
                                             "tier_gene_binding": "matched"})[0]
            rank = "T1 T2 T3 T4 T5 T6 U".split().index(tier) if tier in "T1 T2 T3 T4 T5 T6 U".split() else 6
        keys = {classify_citation(identifier) for identifier in ids}
        overlaps = [group for group in groups if group["keys"] & keys]
        group = {"keys": keys, "ids": ids, "rank": rank, "position": position, "reports": [report]}
        for other in overlaps:
            group["keys"].update(other["keys"])
            group["ids"].extend(other["ids"])
            group["rank"] = min(group["rank"], other["rank"])
            group["position"] = min(group["position"], other["position"])
            group["reports"].extend(other["reports"])
            groups.remove(other)
        groups.append(group)
    # Legacy proposals have no paper checklist. Retain their source order, but
    # never expand a provider's bibliography merely to fill five slots.
    if not reports:
        aliases = source_tokens(row.get("source_identifier_aliases"))
        for position, citation in enumerate(citations):
            kind, identifier = classify_citation(citation)
            if kind not in {"doi", "pmid", "pmcid"}:
                continue
            ids = [citation]
            if position == 0:
                ids.extend(sorted(aliases))
            keys = {classify_citation(value) for value in ids}
            overlaps = [group for group in groups if group["keys"] & keys]
            if overlaps:
                group = {"keys": keys, "ids": ids, "rank": 6, "position": position}
                for other in overlaps:
                    group["ids"].extend(other["ids"])
                    group["keys"].update(other["keys"])
                    group["position"] = min(group["position"], other["position"])
                    groups.remove(other)
                groups.append(group)
            else:
                groups.append({"keys": keys, "ids": ids, "rank": 6, "position": position})
    return sorted(groups, key=lambda group: (group["rank"], group["position"]))


def select_bibliography_papers(row: Mapping, citations: list[str]) -> dict:
    """Read the imported primary papers; retain the old budget for old records.

    The same five-paper handoff is consumed by metadata checks and stage two.
    No downstream Tier sort, new model judgment or candidate gate is added.
    """
    reports = json.loads(normalize_paper_reports(row.get(PAPER_FIELD)) or "[]")
    stamped = any(any(field in report for field in PAPER_HANDOFF_FIELDS) for report in reports)
    if stamped:
        if any(report.get("paper_handoff_policy") != PAPER_HANDOFF_POLICY for report in reports):
            raise ValueError("paper handoff cannot mix imported and unselected reports")
        ordered = _ranked_bibliography_groups(row, citations, recorded=True)
        selected = [group for group in ordered if any(report["paper_handoff_state"] == "primary" for report in group["reports"])]
        unselected = [group for group in ordered if group not in selected]
        return {"policy": PAPER_HANDOFF_POLICY, "max_papers": MAX_BIBLIOGRAPHY_PAPERS,
                "candidate_paper_count": len(ordered), "selected_paper_count": len(selected),
                "selected_source_ids": list(dict.fromkeys(identifier for group in selected for identifier in group["ids"])),
                "unselected_source_ids": list(dict.fromkeys(identifier for group in unselected for identifier in group["ids"])),
                "selection_basis": "imported_primary_handoff; no_metadata_reselection; aliases_share_one_slot",
                "unselected_reason": "outside_primary_handoff; original observations retained",
                "membership_effect": "none; original paper reports retained"}
    ordered = _ranked_bibliography_groups(row, citations)
    selected = ordered[:MAX_BIBLIOGRAPHY_PAPERS]
    ids = list(dict.fromkeys(identifier for group in selected for identifier in group["ids"] if identifier.strip()))
    return {"policy": "tier_relevant_papers_max5_v1", "max_papers": MAX_BIBLIOGRAPHY_PAPERS,
            "candidate_paper_count": len(ordered), "selected_paper_count": len(selected),
            "selected_source_ids": ids,
            "unselected_source_ids": list(dict.fromkeys(identifier for group in ordered[MAX_BIBLIOGRAPHY_PAPERS:] for identifier in group["ids"])),
            "selection_basis": "original_paper_tier_then_stable_report_order; aliases_share_one_slot",
            "unselected_reason": "metadata_budget_exceeded; original observations retained",
            "membership_effect": "none; original paper reports retained"}


def stamp_paper_handoff(rows: list[dict], *, route: str) -> list[dict]:
    """Select once after import Tiers, retaining all original observations.

    A uses the original hypothesis, B its existing unit/reference/species panel.
    A checklist/table panels and the caller's current B policy are stamped;
    old A observation policies keep their unmarked representation. Rejected B claims remain
    audit records and cannot displace the papers eligible for bibliography.
    """
    from .citation_tokens import split_citations
    if route not in {"phenotype_first", "interval_first"}:
        raise ValueError("paper handoff is an A/B import operation")
    result = [dict(row) for row in rows]
    panels: dict[tuple, list[int]] = {}
    reports_by_row = {}
    for index, row in enumerate(result):
        encoded = normalize_paper_reports(row.get(PAPER_FIELD), discard_handoff=True)
        if not encoded:
            continue
        reports_by_row[index] = json.loads(encoded)
        row[PAPER_FIELD] = encoded
        key = ((str(row.get("hypothesis_id", "")) or str(index),)
               if route == "phenotype_first" else
               tuple(str(row.get(field, "")) for field in ("semantic_unit_id", "reference_gene", "reference_species")))
        panels.setdefault(key, []).append(index)
    for key, indices in panels.items():
        if route == "phenotype_first" and not any(report.get("paper_checklist") or report.get("paper_table_policy") == PAPER_TABLE_POLICY
                                                  for index in indices for report in reports_by_row[index]):
            continue
        eligible = [index for index in indices if route != "interval_first" or result[index].get("verification_status") == "proposed"]
        reports = [report for index in eligible for report in reports_by_row[index]]
        citations = [source for index in eligible for source in source_tokens(result[index].get("source_id") or result[index].get("source_ids"))]
        context = result[indices[0]] if route == "phenotype_first" else {}
        groups = _ranked_bibliography_groups({**context, PAPER_FIELD: reports}, citations)
        positions = {token: str(rank) for rank, group in enumerate(groups, 1) for token in source_tokens(group["ids"])}
        scope = "paper_panel_" + hashlib.sha256(json.dumps([route, *key], ensure_ascii=False).encode()).hexdigest()[:20]
        for index in indices:
            stamped_reports = []
            for report in reports_by_row[index]:
                ids = source_tokens([report["source_id"], *split_citations(report["source_identifier_aliases"])])
                rank = min((positions[token] for token in ids if token in positions), key=int, default="") if index in eligible else ""
                stamped_reports.append({**report, "paper_handoff_policy": PAPER_HANDOFF_POLICY,
                    "paper_handoff_scope": scope, "paper_handoff_state": "primary" if rank and int(rank) <= MAX_BIBLIOGRAPHY_PAPERS else "audit_only",
                    "paper_handoff_rank": rank})
            result[index][PAPER_FIELD] = normalize_paper_reports(stamped_reports)
    return result


def _validate_paper_handoff(reports: list[Mapping]) -> None:
    """Validate a recorded panel, never infer a replacement selection.

    Multiple reports/aliases for one paper share its state and rank within each
    scope. Rejected B records may retain an empty audit rank beside that paper.
    Different scopes retain their own selections, including different ranks.
    """
    marked = [report for report in reports if any(field in report for field in PAPER_HANDOFF_FIELDS)]
    if not marked:
        return
    groups = []
    for report in marked:
        if (any(field not in report or not isinstance(report[field], str) for field in PAPER_HANDOFF_FIELDS)
                or report["paper_handoff_policy"] != PAPER_HANDOFF_POLICY
                or not report["paper_handoff_scope"].strip()
                or report["paper_handoff_state"] not in {"primary", "audit_only"}):
            raise ValueError("paper handoff requires all four imported selection fields")
        rank, state = report["paper_handoff_rank"], report["paper_handoff_state"]
        if (rank and (not rank.isdecimal() or int(rank) < 1)) or (
                state == "primary" and (not rank or int(rank) > MAX_BIBLIOGRAPHY_PAPERS)) or (
                state == "audit_only" and rank and int(rank) <= MAX_BIBLIOGRAPHY_PAPERS):
            raise ValueError("paper handoff rank is inconsistent with its selection state")
        ids = source_tokens(report.get("source_id")) | source_tokens(report.get("source_identifier_aliases"))
        overlapping = [group for group in groups if group[0] & ids]
        records = [report]
        for group in overlapping:
            ids.update(group[0])
            records.extend(group[1])
            groups.remove(group)
        groups.append((ids, records))
    slots = {}
    for index, (_, records) in enumerate(groups):
        scopes = {}
        for report in records:
            if not report["paper_handoff_rank"]:
                continue
            scope = report["paper_handoff_scope"]
            selection = (int(report["paper_handoff_rank"]), report["paper_handoff_state"])
            if scope in scopes and scopes[scope] != selection:
                raise ValueError("one paper has inconsistent handoff ranks or states in the same scope")
            scopes[scope] = selection
            slot = (scope, selection[0])
            if slot in slots and slots[slot] != index:
                raise ValueError("different papers share one handoff rank in the same scope")
            slots[slot] = index


def _validate_annotation_json(value: object, path: str) -> None:
    """Annotations retain JSON types; they never supply canonical facts."""
    if value is None or isinstance(value, (str, bool, int)):
        return
    if isinstance(value, float) and math.isfinite(value):
        return
    if isinstance(value, list):
        for index, item in enumerate(value):
            _validate_annotation_json(item, f"{path}[{index}]")
        return
    if isinstance(value, dict) and all(isinstance(key, str) for key in value):
        for key, item in value.items():
            _validate_annotation_json(item, f"{path}.{key}")
        return
    raise ValueError(f"{path} must contain only finite JSON data")


def _separate_research_annotations(report: dict, index: int) -> dict:
    """Copy recognized misplaced metadata without relaxing scientific fields."""
    path = f"paper_reports_json[{index}].{RESEARCH_ANNOTATIONS}"
    present = RESEARCH_ANNOTATIONS in report
    annotations = report.get(RESEARCH_ANNOTATIONS, {})
    if not isinstance(annotations, dict):
        raise ValueError(f"{path} must be an object")
    _validate_annotation_json(annotations, path)
    annotations = dict(annotations)
    result = dict(report)
    for key in sorted(LEGACY_RESEARCH_FIELDS & report.keys()):
        value = result.pop(key)
        _validate_annotation_json(value, f"{path}.{key}")
        if key in annotations and json.dumps(annotations[key], sort_keys=True) != json.dumps(value, sort_keys=True):
            raise ValueError(f"{path}.{key} conflicts with misplaced research field")
        annotations[key] = value
        present = True
    if present:
        # Canonicalize nested key order too, without sharing mutable input data.
        result[RESEARCH_ANNOTATIONS] = json.loads(json.dumps(annotations, ensure_ascii=False, sort_keys=True))
    return result


def normalize_paper_reports(value: object, *, discard_handoff: bool = False) -> str:
    """Accept a JSON array in a model response or JSON text in a TSV cell.

    Empty text is legacy/absent. An explicit [] declares the paper policy with
    no reported paper evidence, rather than inheriting a hypothesis-level Tier.
    """
    if value is None or value == "":
        return ""
    from .paper_transport import expand_paper_observations
    reports = expand_paper_observations(value)
    if not isinstance(reports, list):
        raise ValueError("paper_reports_json must be an array of paper observations")
    from .source_verify import classify_citation
    normalized = {}
    for index, report in enumerate(reports):
        if not isinstance(report, dict):
            raise ValueError("each paper report must be an object")
        report = _separate_research_annotations(report, index)
        if discard_handoff:
            report = {key: value for key, value in report.items() if key not in PAPER_HANDOFF_FIELDS}
        report = {key: "" if value is None and key in REPORT_FIELDS and key not in {"paper_checklist", "evidence_types"}
                  else value for key, value in report.items()}
        report = {**report, "source_identifier_aliases": normalize_source_aliases(report.get("source_identifier_aliases"))}
        if isinstance(report.get("paper_checklist_missing_fields_json"), list):
            report["paper_checklist_missing_fields_json"] = json.dumps(report["paper_checklist_missing_fields_json"])
        invalid = sorted(k for k, v in report.items()
                         if k not in {"evidence_types", "paper_checklist", RESEARCH_ANNOTATIONS} and not isinstance(v, str))
        if invalid:
            raise ValueError(f"paper report fields must be strings: paper_reports_json[{index}].{invalid[0]}")
        unknown = set(report) - {"source_id", *REPORT_FIELDS, *PAPER_HANDOFF_FIELDS, RESEARCH_ANNOTATIONS}
        if unknown:
            raise ValueError(f"unknown paper report field: paper_reports_json[{index}].{sorted(unknown)[0]}")
        for key in ("paper_table_policy", "gene_correspondence", "followup_state", "publication_type", "phenotype_comparison"):
            if key in report:
                report[key] = "unknown" if is_unknown(report[key]) else re.sub(r"[-\s]+", "_", display_text(report[key]).casefold())
        for key in ("reported_gene", "reported_gene_id", "reported_species", "reference_gene_id"):
            if is_unknown(report.get(key)):
                report[key] = "unknown"
        source = report.get("source_id", "").strip()
        if classify_citation(source)[0] not in {"pmid", "pmcid", "doi"}:
            raise ValueError("paper report requires one PMID, PMCID or DOI")
        aliases = source_tokens(report.get("source_identifier_aliases", ""))
        if any(classify_citation(alias)[0] not in {"pmid", "pmcid", "doi"} for alias in aliases):
            raise ValueError("paper aliases must be PMID, PMCID or DOI, never a gene accession")
        table_policy = report.get("paper_table_policy", "")
        if table_policy not in {"", PAPER_TABLE_POLICY}:
            raise ValueError("unknown paper table policy")
        if not report.get("reported_gene", "").strip() and not has_observations(report) and not table_policy:
            raise ValueError("paper report requires its reported gene")
        observations = normalize_evidence_types(report.get("evidence_types"))
        if observations and report.get("phenotype_comparison", "unknown") not in COMPARISONS | {""}:
            raise ValueError("unknown phenotype_comparison")
        tier = report.get("initial_tier", "").split("_", 1)[0]
        if table_policy and not report.get("paper_checklist"):
            report = {**report, "paper_checklist": {}}
        checklist = normalize_checklist(report.get("paper_checklist"))
        inherited_missing = json.loads(report.get("paper_checklist_missing_fields_json") or "[]")
        if not isinstance(inherited_missing, list) or any(not isinstance(key, str) for key in inherited_missing):
            raise ValueError("missing checklist fields must be an array of question names")
        missing_checks = set(inherited_missing)
        if missing_checks - set(CHECKLIST_FIELDS):
            raise ValueError("unknown missing checklist field")
        if checklist:
            raw_checks = report.get("paper_checklist")
            raw_checks = json.loads(raw_checks) if isinstance(raw_checks, str) else raw_checks
            answers = json.loads(checklist)
            checked_unknown = report.get("followup_state") == "completed" and bool(report.get("unknown_reason", "").strip()) and not is_unknown(report.get("unknown_reason"))
            missing_checks = {key for key in missing_checks if answers[key] == "unknown" and not checked_unknown}
            missing_checks.update(key for key in CHECKLIST_FIELDS
                                  if key not in raw_checks or unfilled_answer(raw_checks[key]))
        if report.get("followup_state", "") not in {"", "not_attempted", "completed"}:
            raise ValueError("followup_state must be not_attempted or completed")
        if checklist:
            tier = checklist_tier(checklist)
        relation = report.get("gene_correspondence", "").strip()
        if checklist and not relation and not table_policy:
            raise ValueError("paper checklist requires gene_correspondence")
        if relation and relation not in GENE_RELATIONS:
            raise ValueError("unknown gene_correspondence")
        if relation in {"same_gene", "homolog"} and not report.get("reference_gene_id", "").strip() and not table_policy:
            raise ValueError("paper judgment requires its existing reference_gene_id")
        if (relation or not observations) and tier not in {"T1", "T2", "T3", "T4", "T5", "T6", "U"}:
            raise ValueError("paper report requires T1-T6/U initial_tier")
        if (relation or not observations) and not report.get("initial_tier_reason", "").strip() and not table_policy:
            raise ValueError("paper report requires its own initial_tier_reason")
        item = {"source_id": source, **{k: report.get(k, "").strip() for k in REPORT_FIELDS if k not in {"evidence_types", "paper_checklist"}}}
        item["paper_checklist"] = checklist
        item["paper_checklist_missing_fields_json"] = json.dumps(sorted(missing_checks), separators=(",", ":"))
        if checklist:
            item["initial_tier"] = tier
        item["evidence_types"] = observations
        # Optional machine annotations must not add empty keys to old records.
        item.update({key: report[key] for key in PAPER_HANDOFF_FIELDS if key in report})
        if RESEARCH_ANNOTATIONS in report:
            item[RESEARCH_ANNOTATIONS] = report[RESEARCH_ANNOTATIONS]
        encoded = json.dumps(item, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
        normalized[encoded] = item
    # Metadata must not reorder equally ranked scientific reports at the
    # five-paper handoff. Keep distinct notes, but order by the original facts.
    ordered = sorted(normalized.values(), key=lambda item: (
        json.dumps({k: v for k, v in item.items() if k != RESEARCH_ANNOTATIONS}, ensure_ascii=False, sort_keys=True, separators=(",", ":")),
        json.dumps(item, ensure_ascii=False, sort_keys=True, separators=(",", ":"))))
    _validate_paper_handoff(ordered)
    return json.dumps(ordered, ensure_ascii=False, separators=(",", ":"))


def imported_paper_fields(record: Mapping, source_id: str | None = None) -> dict:
    """Keep each B citation bound to its own report, and A reports together."""
    top_fields = {"source_identifier_aliases": normalize_source_aliases(record.get("source_identifier_aliases"))}
    if has_observations(record):
        top_fields["evidence_types"] = normalize_evidence_types(record.get("evidence_types"))
    encoded = normalize_paper_reports(record.get(PAPER_FIELD), discard_handoff=True)
    if not encoded:
        return top_fields
    reports = json.loads(encoded)
    declared = source_tokens(record.get("source_ids", record.get("source_id", "")))
    for report in reports:
        ids = source_tokens([report["source_id"]]) | source_tokens(report["source_identifier_aliases"])
        if not ids & declared:
            raise ValueError("paper report is not bound to a declared source ID")
    if source_id is None:
        return {**top_fields, PAPER_FIELD: encoded}
    matching = [r for r in reports if source_tokens([source_id]) &
                (source_tokens([r["source_id"]]) | source_tokens(r["source_identifier_aliases"]))]
    fields = {PAPER_FIELD: normalize_paper_reports(matching)}
    # Do not copy an unrelated paper's title, phenotype, or Tier to this source.
    fields.update({k: "" for k in REPORT_FIELDS})
    if matching:
        best = min(matching, key=lambda r: ("T1 T2 T3 T4 T5 T6 U".split().index(r["initial_tier"].split("_", 1)[0]) if not has_observations(r) else 6,
                                            r["source_id"], json.dumps({k: v for k, v in r.items() if k != RESEARCH_ANNOTATIONS}, sort_keys=True)))
        fields.update({k: best[k] for k in REPORT_FIELDS})
    elif source_tokens([source_id]) & declared:
        from .source_verify import classify_citation
        if classify_citation(source_id)[0] == "database_accession" and not reports and has_observations(record):
            # An explicit empty paper list must not erase an original DB-backed
            # function/annotation judgment. Never borrow facts from a paper.
            fields.update(evidence_types=top_fields.get("evidence_types", ""),
                          reported_function=str(record.get("reported_function") or "").strip(),
                          reported_annotation=str(record.get("reported_annotation") or "").strip())
        elif classify_citation(source_id)[0] == "database_accession" and "annotation" in json.loads(top_fields.get("evidence_types", "[]")):
            # Keep only this declared DB source's annotation. Paper functions,
            # phenotypes and Tier declarations remain attached to their paper.
            fields.update(evidence_types=normalize_evidence_types(["annotation"]),
                          reported_annotation=str(record.get("reported_annotation") or "").strip())
    return fields


def bound_paper_rows(row: Mapping) -> list[dict]:
    """Expand only reports belonging to this already bound claim/source row."""
    from .reference_shared_gene import validate_shared_reference_transports
    validate_shared_reference_transports(row)
    encoded = normalize_paper_reports(row.get(PAPER_FIELD))
    if not encoded:
        return []
    bound = source_tokens(row.get("source_id") or row.get("source_ids", ""))
    result = []
    for report in json.loads(encoded):
        ids = source_tokens([report["source_id"]]) | source_tokens(report["source_identifier_aliases"])
        if not bound & ids:
            continue
        result.append({**row, **report, PAPER_FIELD: "",
                       "source_id": report["source_id"], "source_ids": report["source_id"],
                       "reference_species": report["reported_species"],
                       "proxy_reference_species": row.get("reference_species", ""),
                       "relevance_policy": "reported_function_v2",
                       "tier_aggregation": PAPER_POLICY})
        if report.get("gene_correspondence") or (report.get("paper_table_policy") and report.get("paper_checklist")):
            relation = report.get("gene_correspondence")
            anchor = gene_binding(row, report["reference_gene_id"])
            result[-1].update(tier_classification_policy=CHECKLIST_POLICY if report.get("paper_checklist") else JUDGMENT_POLICY,
                             tier_gene_binding="matched" if anchor == "matched" and relation in {"same_gene", "homolog"}
                             and report.get("reported_gene", "").strip() and not is_unknown(report.get("reported_gene")) else "unresolved")
            result[-1].update(checklist_consistency(result[-1]))
        elif has_observations(report):
            result[-1].update(tier_classification_policy=POLICY,
                             tier_gene_binding=gene_binding(row, report["reported_gene"], report["reported_species"]))
    return result


def paper_summary(classified: list[tuple[tuple, Mapping]]) -> tuple[dict, list[dict]]:
    """Count distinct declared papers, never claim replication or independence.

    Alias groups are transitive; a review remains a separate publication even
    when it cites an original study. Underlying-study links stay explicit.
    """
    from .source_verify import classify_citation
    _validate_paper_handoff([row for _, row in classified if row.get("tier_aggregation") == PAPER_POLICY])
    groups = []
    for tier, row in classified:
        if row.get("tier_aggregation") != PAPER_POLICY:
            continue
        ids = {classify_citation(source) for source in
               source_tokens(row.get("source_id")) | source_tokens(row.get("source_identifier_aliases"))}
        overlapping = [group for group in groups if group[0] & ids]
        records = [(tier, row)]
        for group in overlapping:
            ids.update(group[0])
            records.extend(group[1])
            groups.remove(group)
        groups.append((ids, records))
    counts, summaries = {}, []
    order = "T1 T2 T3 T4 T5 T6 U".split()
    for _, records in groups:
        tier, best = min(records, key=lambda item: (
            order.index(item[0][0].split("_", 1)[0]), item[1].get("tier_classification_state") == "partial", item[1]["source_id"],
            json.dumps({k: v for k, v in item[1].items() if k != RESEARCH_ANNOTATIONS}, sort_keys=True)))
        label = tier[0].split("_", 1)[0]
        counts[label] = counts.get(label, 0) + 1
        summary = {key: best.get(key, "") for key in
                   ("source_id", "source_title", "reported_gene", "reported_gene_id", "reported_species", "reference_gene_id", "gene_correspondence", "paper_checklist", "evidence_types",
                    "reported_phenotype", "reported_process", "reported_function", "reported_annotation", "phenotype_comparison",
                    "tier_classification_policy", "tier_classification_state", "tier_missing_fields_json", "tier_gene_binding",
                    "homolog_projection_json", "proxy_reference_species", *PAPER_TABLE_FIELDS)}
        summary.update(tier=label, tier_reason=tier[2], independence_state="not_assessed")
        if best.get("tier_classification_policy") == CHECKLIST_POLICY:
            summary.update({key: best.get(key, "") for key in CHECKLIST_CONSISTENCY_FIELDS})
        handoff_panels = {json.dumps({key: row[key] for key in PAPER_HANDOFF_FIELDS}, ensure_ascii=False, sort_keys=True)
                         for _, row in records if row.get("paper_handoff_policy") == PAPER_HANDOFF_POLICY}
        if handoff_panels:
            # A display representative can belong to another original panel.
            # Preserve every panel's selection without changing this paper's Tier.
            summary["paper_handoff_panels"] = sorted((json.loads(value) for value in handoff_panels), key=lambda panel: (
                panel["paper_handoff_scope"], int(panel["paper_handoff_rank"] or 1000000), panel["paper_handoff_state"]))
        warnings = {json.dumps({key: row.get(key, "") for key in
                    ("source_id", "reported_gene", "reported_gene_id", "reported_species", *CHECKLIST_CONSISTENCY_FIELDS)},
                    ensure_ascii=False, sort_keys=True) for _, row in records if row.get("consistency_warning")}
        if warnings:
            # One publication may contain several gene reports. A consistent
            # display representative must not hide another record's warning.
            summary["consistency_warnings"] = [json.loads(value) for value in sorted(warnings)]
        annotations = {json.dumps(row[RESEARCH_ANNOTATIONS], ensure_ascii=False, sort_keys=True, separators=(",", ":"))
                       for _, row in records if RESEARCH_ANNOTATIONS in row}
        if annotations:
            summary[RESEARCH_ANNOTATIONS] = [json.loads(value) for value in sorted(annotations)]
        # Preserve differing source descriptions instead of inferring which is true.
        for field in ("publication_type", "underlying_source_ids", "corroboration_note"):
            summary[field] = sorted({str(row[field]) for _, row in records if row.get(field)})
        summaries.append(summary)
    return counts, sorted(summaries, key=lambda item: (order.index(item["tier"]), item["source_id"]))
