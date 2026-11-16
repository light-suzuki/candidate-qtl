"""Readable view of the original A/B paper reports and bounded correction requests.

This does not review source meaning, choose members, or execute another model.
The report JSON remains the input to the existing Tier decision table.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

from .observation_tiers import CHECKLIST_FIELDS, CHECKLIST_POLICY, JUDGMENT_POLICY, POLICY, checklist_tier, classify_observation, gene_binding, is_unknown
from .observation_tiers import CHECKLIST_CONSISTENCY_FIELDS, checklist_consistency
from .paper_reports import PAPER_FIELD, PAPER_TABLE_POLICY, PAPER_HANDOFF_FIELDS, PAPER_HANDOFF_POLICY, RESEARCH_ANNOTATIONS, normalize_paper_reports, select_bibliography_papers, source_tokens
from .source_verify import classify_citation

_ORDER = "T1 T2 T3 T4 T5 T6 U".split()
_ANCHOR_FIELDS = ("reference_gene", "reference_gene_id", "reference_gene_name", "source_reference_gene_id", "proposed_stable_id")
_FACT_FIELDS = ("source_title", "publication_type", "reported_gene", "reported_gene_id", "reported_species", "reference_gene_id", "gene_correspondence",
                "initial_tier_reason", "source_excerpt", "source_location")
_OPTIONAL_INITIAL_FIELDS = frozenset({"publication_type", "source_excerpt", "source_location"})
TABLE_COLUMNS = (
    "record_id", "proposal_id", "record_kind", "paper_table_policy", "source_id", "source_identifier_aliases", "source_title", "publication_type",
    "reported_gene", "reported_gene_id", "reported_species", "reference_gene_id", "gene_correspondence",
    "correspondence_key", "reference_anchor_state", "member_interval_genes", "projection_state",
    "source_excerpt", "source_location", "reported_phenotype", "reported_function", "reported_annotation",
    *CHECKLIST_FIELDS, "unknown_reason", "missing_fields_json", "unknown_fields_json", "record_state",
    "computed_tier", "applicable_tier", "tier_state", "initial_tier_reason", "selected_for_metadata",
    *CHECKLIST_CONSISTENCY_FIELDS,
    "metadata_state", "resolved_title", "resolved_doi", "bibliography_checks_json",
    "followup_question", "followup_state", "followup_action",
)


def _text(value: object) -> str:
    return str(value or "").strip()


def _digest(value: object) -> str:
    return hashlib.sha256(json.dumps(value, ensure_ascii=False, sort_keys=True).encode()).hexdigest()


def _upper_tier_questions(checks: dict, tier: str) -> list[str]:
    """Only explicit phenotype facts can trigger an upper-Tier follow-up.

    Merely imagining that every unknown could be yes is not a research lead.
    """
    if tier == "T1":
        return []
    if tier == "T2":
        keys = ("focal_species",)
        return [] if any(checks.get(key) == "no" for key in keys) else [key for key in keys if checks.get(key) == "unknown"]
    if checks.get("phenotype_reported") == "yes":
        return [key for key in ("equivalent_target_phenotype", "similar_cellular_or_developmental_phenotype")
                if checks.get(key) == "unknown"]
    if checks.get("phenotype_reported") == "unknown" and any(checks.get(key) == "yes" for key in
            ("equivalent_target_phenotype", "similar_cellular_or_developmental_phenotype")):
        return ["phenotype_reported"]
    return []


def write_paper_evidence_table(records: list[dict], out: Path, *, route: str, source_artifact: str,
                             metadata_rows: list[dict] | None = None, audit_focal_species: str = "") -> dict:
    """Export one row per proposal/paper/reported gene, without changing records.

    Blanks remain distinct from explicit unknowns, including after normalization.
    Exact frozen reference-anchor checks are memoized; this is not biological
    identity verification or a cross-route cache of semantic judgments.
    Case context supplied only for the checklist audit must not change the
    inputs or Tier of legacy observations, including mixed-policy proposals.
    """
    out.mkdir(parents=True, exist_ok=True)
    table, tasks, panels, seen, anchors, selection_reports = [], {}, {}, set(), {}, {}
    for index, record in enumerate(records, 1):
        from .reference_shared_gene import validate_shared_reference_transports
        validate_shared_reference_transports(record)
        owner = _text(record.get("hypothesis_id") or record.get("claim_id")) or f"proposal_{index:03d}"
        if route == "interval_first":
            owner = owner.rsplit("-s", 1)[0]
        reports = json.loads(normalize_paper_reports(record.get(PAPER_FIELD)) or "[]")
        declared = source_tokens(record.get("source_id") or record.get("source_ids"))
        covered = set().union(*(source_tokens([r["source_id"], *r["source_identifier_aliases"].split(";")]) for r in reports))
        # A cited paper with an absent report is recording work, not a T5 result.
        for source in sorted(declared - covered):
            if classify_citation(source)[0] in {"doi", "pmid", "pmcid"}:
                reports.append({"source_id": source, "source_title": _text(record.get("source_title")),
                                "paper_checklist_missing_fields_json": json.dumps(list(CHECKLIST_FIELDS)),
                                "missing_report": True})
        for report in reports:
            fingerprint = _digest([owner, report])
            if fingerprint in seen:
                continue
            seen.add(fingerprint)
            record_id = "paper_" + fingerprint[:20]
            checks = json.loads(report.get("paper_checklist") or "{}")
            missing = list(json.loads(report.get("paper_checklist_missing_fields_json") or "[]"))
            missing += [key for key in _FACT_FIELDS if not _text(report.get(key))]
            initial_record = (report.get("paper_table_policy") == PAPER_TABLE_POLICY
                              or record.get("candidate_admission_policy") == "initial_bibliography_v1")
            if (record.get("candidate_admission_policy") == "source_function_v1"
                    or record.get("verification_policy") == "content_v1"):
                initial_record = False
            optional = _OPTIONAL_INITIAL_FIELDS if initial_record else frozenset()
            unknown = [key for key in CHECKLIST_FIELDS if checks.get(key) == "unknown" and key not in missing]
            unknown += [key for key in _FACT_FIELDS if is_unknown(report.get(key))]
            if any(key not in optional for key in unknown) and (not _text(report.get("unknown_reason")) or is_unknown(report.get("unknown_reason"))):
                missing.append("unknown_reason")
            required_missing = [key for key in missing if key not in optional]
            relation = report.get("gene_correspondence")
            target = report.get("reference_gene_id") if relation else report.get("reported_gene")
            species = "" if relation else _text(report.get("reported_species"))
            anchor_key = tuple(_text(record.get(key)) for key in (*_ANCHOR_FIELDS, "reference_species")) + (_text(target), species)
            if anchor_key not in anchors:
                anchors[anchor_key] = gene_binding(record, target, species)
            anchor = anchors[anchor_key]
            applicable = anchor == "matched" and (not relation or relation in {"same_gene", "homolog"})
            if report.get("paper_table_policy") and checks:
                applicable = applicable and relation in {"same_gene", "homolog"} and bool(_text(report.get("reported_gene"))) and not is_unknown(report.get("reported_gene"))
            tier = checklist_tier(checks) if checks else _text(report.get("initial_tier")).split("_", 1)[0] or "U"
            applied_tier = classify_observation({**record, **report,
                "source_id": report["source_id"], "tier_gene_binding": "matched" if applicable else "unresolved",
                "tier_classification_policy": CHECKLIST_POLICY if checks else JUDGMENT_POLICY if relation else POLICY})[0]
            if not checks and not report.get("missing_report"):
                tier = applied_tier
            consistency = checklist_consistency({**record, **report,
                "focal_species": audit_focal_species or record.get("focal_species", "")})
            warning = consistency.get("consistency_warning", "")
            questions = _upper_tier_questions(checks, tier)
            if warning:
                questions.extend(("reported_species", "focal_species"))
            if tier in {"T1", "T2"} and not applicable:
                questions.extend(key for key in ("reported_gene", "reference_gene_id", "gene_correspondence")
                                 if not _text(report.get(key)) or is_unknown(report.get(key)))
                if anchor != "matched":
                    questions.append("reference_gene_id")
                questions = sorted(set(questions))
            explicit_question = _text(report.get("followup_question"))
            attempted = report.get("followup_state") == "completed"
            action = ("unfilled_after_followup" if attempted else "complete_record") if required_missing else (
                "clarify_record_consistency" if not attempted and warning else
                "clarify_upper_tier_lead" if not attempted and (questions or explicit_question) else "none")
            state = "unfilled" if missing else "checked_unknown" if unknown else "recorded"
            row = {key: _text(report.get(key)) for key in TABLE_COLUMNS}
            if RESEARCH_ANNOTATIONS in report:
                row["research_annotations_json"] = json.dumps(report[RESEARCH_ANNOTATIONS], ensure_ascii=False, sort_keys=True, separators=(",", ":"))
            row.update({key: "" if key in missing else checks.get(key, "") for key in CHECKLIST_FIELDS})
            # Names/IDs, species and relation stay separate; repeated names alone
            # do not authorize reuse of a biological correspondence decision.
            identity = [_text(report.get(key)) for key in
                        ("reported_gene_id", "reported_gene", "reported_species", "reference_gene_id", "gene_correspondence")]
            row.update(record_id=record_id, proposal_id=owner,
                       **consistency,
                       record_kind="missing_paper_report" if report.get("missing_report") else "paper_report",
                       correspondence_key=_digest(identity), reference_anchor_state=anchor,
                       member_interval_genes=_text(record.get("member_interval_genes")),
                       projection_state="frozen_unit_members" if route == "interval_first" else "not_available_before_projection",
                       missing_fields_json=json.dumps(sorted(set(missing))), unknown_fields_json=json.dumps(unknown),
                       record_state=state, computed_tier=tier, applicable_tier=applied_tier,
                       tier_state="provisional_missing_record" if required_missing else "correspondence_unresolved" if not applicable else "assessed",
                       metadata_state="not_checked_at_proposal_import; use existing bibliography verification",
                       followup_action=action)
            table.append(row)
            selection_reports[record_id] = report
            panels.setdefault(owner, []).append(row)
            if action in {"complete_record", "clarify_upper_tier_lead", "clarify_record_consistency"}:
                task = tasks.setdefault(owner, [])
                task.append({"record_id": record_id, "source_id": report["source_id"], "action": action,
                             "missing_fields": sorted(set(required_missing)), "questions": questions,
                             "followup_question": explicit_question, "current_report": report,
                             **consistency})
                if tier != "T1" and checks.get("focal_species") == "yes" and checks.get("same_gene_as_focal") == "unknown":
                    task[-1]["sequence_identity_followup"] = {
                        "state": "public_sequence_lookup_required",
                        "before_projection": "Retrieve the reported gene's official accession/version and CDS protein_id; "
                            "save raw record, taxid, cultivar and source location. Missing sequence remains not_evaluated.",
                        "after_route_freeze_parent_only": "Run reference-gene-resolve --sequence-requests --post-freeze "
                            "against the supplied whole proteome and exact protein-gene map. Compare both coverages, "
                            "all competing copies including unplaced scaffolds, and annotation model boundaries. "
                            "Use source-verify --reference-bridge for the paper-bound view; never send interval results to blind A.",
                        "assessment": "Keep locus identity separate from whole annotation model identity. "
                            "A sequence match alone does not answer phenotype questions, change Tier, or verify support."}
    for panel in panels.values():
        # Read the imported primary selection; old records retain their budget.
        reports = [selection_reports[row["record_id"]]
                   for row in panel if row["record_kind"] == "paper_report"]
        if reports:
            selected = source_tokens(select_bibliography_papers({PAPER_FIELD: reports},
                                      [r["source_id"] for r in reports])["selected_source_ids"])
        else:
            selected = set()
        for row in panel:
            row["selected_for_metadata"] = "yes" if source_tokens(row["source_id"]) & selected else "no"
            report = selection_reports[row["record_id"]]
            if report.get("paper_handoff_policy") == PAPER_HANDOFF_POLICY:
                row.update({field: report[field] for field in PAPER_HANDOFF_FIELDS})
                row["selected_for_handoff"] = "yes" if report["paper_handoff_state"] == "primary" else "no"
    if metadata_rows is not None:
        audits = {}
        for audit in metadata_rows:
            owner = _text(audit.get("hypothesis_id"))
            if route == "interval_first":
                owner = owner.rsplit("-s", 1)[0]
            audits.setdefault(owner, []).append(audit)
        for row in table:
            ids = source_tokens([row["source_id"], *row["source_identifier_aliases"].split(";")])
            matches = [audit for audit in audits.get(row["proposal_id"], []) if source_tokens(audit["citation"]) & ids]
            row["bibliography_checks_json"] = json.dumps(matches, ensure_ascii=False, sort_keys=True)
            row["metadata_state"] = ";".join(sorted({audit["bibliography_state"] for audit in matches})) or "not_checked"
            row["resolved_title"] = ";".join(sorted({audit["title"] for audit in matches if audit.get("title")}))
            row["resolved_doi"] = ";".join(sorted({audit["doi"] for audit in matches if audit.get("doi")}))
    table.sort(key=lambda row: (row["proposal_id"], _ORDER.index(row["applicable_tier"]), row["record_id"]))
    path = out / "paper_evidence_table.tsv"
    columns = (*TABLE_COLUMNS, *PAPER_HANDOFF_FIELDS, "selected_for_handoff") if any(row.get("paper_handoff_policy") for row in table) else TABLE_COLUMNS
    if any("research_annotations_json" in row for row in table):
        columns = (*columns, "research_annotations_json")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(table)
    request_path = out / "paper_followup_requests.json"
    source = Path(source_artifact)
    exhausted = any(row["followup_action"] == "unfilled_after_followup" for row in table)
    packet = {"policy": PAPER_TABLE_POLICY, "route": route, "max_followup_batches": 1,
              "source_artifact": source_artifact, "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
              "execution_state": "not_executed" if tasks else "followup_exhausted" if exhausted else "not_required",
              "instructions": "Return these questions to the original researcher with the original route inputs. "
                  "Correct only the listed paper records in a new copy of the original payload/response; keep other facts, IDs and scope unchanged. "
                  "Use unknown with checked sources/reasons if unresolved; mark followup_state=completed after the bounded attempt. "
                  "Re-import into a new output generation. New B sources must be extended/frozen before import. "
                  "Do not fabricate answers, force five reports, promote support, expose B/interval inputs to blind A, or rewrite a freeze.",
              "requests": [{"proposal_id": owner, "records": rows} for owner, rows in sorted(tasks.items())]}
    request_path.write_text(json.dumps(packet, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    state = "needs_completion" if tasks or exhausted else "no_paper_reports" if not table else "recorded"
    return {"policy": PAPER_TABLE_POLICY, "recording_state": state, "paper_row_count": len(table),
            "followup_proposal_count": len(tasks), "reference_anchor_check_count": len(anchors),
            "next_action": "return targeted requests to the original researcher before claiming recording completion"
                if tasks else "bounded followup exhausted; preserve incomplete recording without another automatic retry"
                if exhausted else "retain checked unknowns; continue existing source/identity/projection workflow",
            "membership_effect": "none", "research_execution": "not measured by record completeness",
            "outputs": {"paper_evidence_table": str(path), "paper_followup_requests": str(request_path)},
            "sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest(),
                       request_path.name: hashlib.sha256(request_path.read_bytes()).hexdigest()}}
