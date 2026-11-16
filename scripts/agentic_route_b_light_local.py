"""Offline submission checks and explicit, hash-bound original-researcher corrections."""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
from typing import Any

if __package__:
    from . import agentic_route_b_light as light
    from . import agentic_route_b_material_read as material
else:
    import agentic_route_b_light as light
    import agentic_route_b_material_read as material

PREFLIGHT_CONTRACT = "agentic-route-b-light-preflight-1.0.0"
PATCH_CONTRACT = "agentic-route-b-light-patch-1.0.0"
REVIEW_REQUEST_CONTRACT = "agentic-route-b-light-review-request-1.0.0"
REVIEW_CONTRACT = "agentic-route-b-light-batch-review-1.0.0"
DECISIONS_CONTRACT = "agentic-route-b-light-batch-decisions-1.0.0"
BATCH_APPLY_CONTRACT = "agentic-route-b-light-batch-apply-1.0.0"
POSTPROCESS_CONTRACT = "agentic-route-b-light-postprocess-1.0.0"
MISSING = {"$missing": True}
PAPER_FIELDS = {"paper_id", "title", "reported_gene", "reported_gene_id", "reported_species", "finding",
                "gene_correspondence", "checklist", "reported_phenotype", "reported_process", "reported_organ", "reported_stage",
                "source_excerpt", "source_location", "unknown_reason", "publication_type", "counterevidence", "research_annotations",
                "phenotype_comparison", "underlying_source_ids", "corroboration_note", "followup_question", "followup_state",
                "reported_function", "reported_annotation", "paper_checklist_reasons", "read_scope", "conditional_focal_reasoning"}
GENE_FIELDS = {"reference_gene", "state", "basis", "function", "reason", "tier", "literature"}
ASSESSMENT_FIELDS = {"observation_id", "gene_correspondence", "checklist", "reason", "counterevidence",
                     "paper_checklist_reasons", "research_annotations", "conditional_focal_reasoning"}


def _hash(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _assessments(record):
    value = record.get("paper_assessments", [])
    return value if isinstance(value, list) else []


def _write_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def _binding(bundle_arg, bundle_sha256, light_index_arg, light_index_sha256, batch_id):
    _, bundle = light._load_bundle(bundle_arg, bundle_sha256)
    index, batches = light._validated_index(bundle, bundle_sha256, Path(light_index_arg).resolve(), light_index_sha256)
    batch = next((row for row in batches if row["bid"] == batch_id), None)
    if batch is None:
        raise ValueError("unknown light batch_id")
    units = {unit["unit_id"]: unit for unit in bundle["units"]}
    roster = {gene: {**row, "unit_id": uid} for uid in batch["unit_ids"] for gene, row in light._roster(units[uid]).items()}
    return bundle, index, batch, roster


def _records(path: Path, science: bool):
    if science:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(payload, dict) or payload.get("contract") != light.SCIENCE_CONTRACT:
            raise ValueError("science results contract mismatch")
        if not isinstance(payload.get("gene_decisions"), list) or not isinstance(payload.get("paper_observations"), list):
            raise ValueError("science results require record arrays")
        records = [("gene_decision", row.get("interval_gene_id"), row) for row in payload["gene_decisions"] if isinstance(row, dict)]
        records += [("paper_observation", row.get("observation_id"), row) for row in payload["paper_observations"] if isinstance(row, dict)]
        if len(records) != len(payload["gene_decisions"]) + len(payload["paper_observations"]):
            raise ValueError("results record must be an object")
    else:
        payload = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
        if any(not isinstance(row, dict) for row in payload):
            raise ValueError("legacy results record must be an object")
        records = [("legacy_gene", row.get("interval_gene_id"), row) for row in payload if set(row) != {"bibliography"}]
        records += [("bibliography_header", str(number), row) for number, row in enumerate(payload) if set(row) == {"bibliography"}]
    keyed = {}
    for kind, key, record in records:
        if not isinstance(key, str) or not key.strip() or (kind, key) in keyed:
            raise ValueError("missing or duplicate record ID")
        keyed[kind, key] = record
    return payload, keyed


def _packet(bundle_arg, bundle_sha256, light_index_arg, light_index_sha256, batch_id, results_arg):
    bundle, index, batch, roster = _binding(bundle_arg, bundle_sha256, light_index_arg, light_index_sha256, batch_id)
    path = Path(results_arg).resolve()
    science = index["contract"] == light.V2_CONTRACT
    findings, states = [], {}

    def finding(gene, field, issue, observation=None):
        item = {"gene_id": gene, "field": field, "issue": issue}
        if observation is not None:
            item["observation_id"] = observation
        findings.append(item)

    try:
        _, records = _records(path, science)
    except (ValueError, TypeError, KeyError) as error:
        records = {}
        finding("", "results", str(error))
    try:
        rows, issues = (light._science_rows if science else light._results_rows)(path)
    except (ValueError, TypeError, KeyError) as error:
        rows, issues = {}, []
        finding("", "results", str(error))
    for issue in issues:
        finding(issue.get("gene_id", ""), "results", issue.get("detail", issue["issue"]))
    for gene in rows:
        if gene not in roster:
            finding(gene, "interval_gene_id", "gene_not_in_batch")
    if science:
        from candidate_qtl.source_verify import classify_citation
        # Inspect the original copy answers before the transport parser rejects
        # an incomplete decision and removes its expanded paper list.
        for (kind, gene), decision in records.items():
            if kind != "gene_decision" or decision.get("state") != "assessed":
                continue
            for number, assessment in enumerate(_assessments(decision)):
                if not isinstance(assessment, dict):
                    continue
                answers = assessment.get("checklist")
                observation_id = assessment.get("observation_id")
                for key in light.CHECKLIST_FIELDS:
                    if not isinstance(answers, dict) or answers.get(key) not in ("yes", "no", "unknown"):
                        finding(gene, f"paper_assessments[{number}].checklist.{key}",
                                "missing_or_invalid_original_answer", observation_id)
                if isinstance(answers, dict) and set(answers) - set(light.CHECKLIST_FIELDS):
                    finding(gene, f"paper_assessments[{number}].checklist", "unknown_checklist_key", observation_id)
        for (kind, key), observation in records.items():
            if kind != "paper_observation":
                continue
            users = [gene for (record_kind, gene), record in records.items() if record_kind == "gene_decision"
                     and any(isinstance(item, dict) and item.get("observation_id") == key for item in _assessments(record))]
            for field in ("paper_id", "title", "reported_gene", "reported_gene_id", "reported_species", "finding"):
                if not isinstance(observation.get(field), str) or not observation[field].strip():
                    for gene in users or [""]:
                        finding(gene, field, "required_original_source_field", key)
            if classify_citation(str(observation.get("paper_id", "")))[0] not in {"pmid", "pmcid", "doi"}:
                for gene in users or [""]:
                    finding(gene, "paper_id", "invalid_paper_identifier", key)
    for gene, member in roster.items():
        row = rows.get(gene)
        states[gene] = "not_assessed"
        if row is None:
            finding(gene, "interval_gene_id", "missing_member_response")
            continue
        if row.get("_science_state") in {"not_assessed", "retrieval_failed"} or row.get("tier") == "not_assessed" and not row.get("_conflict"):
            states[gene] = row.get("_science_state", "not_assessed")
            continue
        tier = row.get("tier")
        if tier not in light.STATE_ORDER:
            finding(gene, "tier", "invalid_or_missing_tier_or_science_decision")
        if str(row.get("reference_gene", "")).casefold() not in {ref.casefold() for ref in member["reference_endpoints"]}:
            # v1 intentionally permits no-reference T6/U bookkeeping.
            if science or tier not in {"T6", "U"}:
                finding(gene, "reference_gene", "reference_not_in_frozen_endpoints")
        if tier not in {"T6", "U"} or science and row.get("papers"):
            for field in ("function", "reason"):
                if not isinstance(row.get(field), str) or not row[field].strip():
                    finding(gene, field, "required_nonempty_original_judgment")
        # Preserve old v1 T6/U's deliberate non-use of ancillary papers.
        papers = row.get("papers", []) if science or tier not in {"T6", "U"} else []
        if not isinstance(papers, list):
            finding(gene, "papers", "papers_must_be_array")
            papers = []
        if tier in light.PAPER_TIERS and not papers:
            finding(gene, "papers", "paper_basis_requires_original_observation")
        if tier in {"T4", "T5"} and not papers:
            target_unit = next(unit for unit in bundle["units"] if unit["unit_id"] == member["unit_id"])
            if not light._frozen_database_source(target_unit, row.get("reference_gene", "")):
                finding(gene, "reference_gene", "no_frozen_database_source_for_paperless_basis")
        for number, paper in enumerate(papers):
            prefix = f"papers[{number}]"
            if not isinstance(paper, dict):
                finding(gene, prefix, "paper_must_be_object")
                continue
            for field in ("paper_id", "title", "reported_gene", "reported_gene_id", "reported_species", "finding"):
                if not isinstance(paper.get(field), str) or not paper[field].strip():
                    finding(gene, f"{prefix}.{field}", "required_original_source_field")
            from candidate_qtl.source_verify import classify_citation
            if classify_citation(str(paper.get("paper_id", "")))[0] not in {"pmid", "pmcid", "doi"}:
                finding(gene, f"{prefix}.paper_id", "invalid_paper_identifier")
            if paper.get("gene_correspondence") not in {"same_gene", "homolog", "different_gene", "unknown"}:
                finding(gene, f"{prefix}.gene_correspondence", "invalid_correspondence")
            answers = paper.get("checklist")
            if not isinstance(answers, dict):
                finding(gene, f"{prefix}.checklist", "all_nine_original_answers_required")
            else:
                for key in light.CHECKLIST_FIELDS:
                    if answers.get(key) not in ("yes", "no", "unknown"):
                        finding(gene, f"{prefix}.checklist.{key}", "missing_or_invalid_original_answer")
                if set(answers) - set(light.CHECKLIST_FIELDS):
                    finding(gene, f"{prefix}.checklist", "unknown_checklist_key")
        states[gene] = "invalid" if any(item["gene_id"] == gene for item in findings) else "assessed"
    packet = {"contract": PREFLIGHT_CONTRACT, "case_id": bundle["case_id"],
              "base": {"bundle_sha256": bundle_sha256, "light_index_sha256": light_index_sha256,
                       "batch_id": batch_id, "results_sha256": light._sha(path)},
              "roster_sha256": batch["roster_sha256"], "findings": findings,
              "member_states": states, "counts": {state: list(states.values()).count(state)
                                                       for state in ("assessed", "not_assessed", "invalid", "retrieval_failed")},
              "scientific_completion": "complete" if all(value == "assessed" for value in states.values()) and not findings else "partial",
              "scope": "local schema/identity/coverage only; no metadata retrieval or biological verification"}
    return packet, records, science


def light_preflight(bundle_arg: str, bundle_sha256: str, light_index_arg: str, light_index_sha256: str,
                    batch_id: str, results_arg: str, out_arg: str) -> dict[str, Any]:
    packet, _, _ = _packet(bundle_arg, bundle_sha256, light_index_arg, light_index_sha256, batch_id, results_arg)
    out = Path(out_arg)
    if out.exists():
        raise ValueError("preflight output must be new")
    _write_new(out, packet)
    return {"packet": str(out.resolve()), "packet_sha256": light._sha(out),
            "counts": packet["counts"], "scientific_completion": packet["scientific_completion"], "finding_count": len(packet["findings"]),
            "issue_count": len(packet["findings"]), "state": "needs_correction" if packet["findings"] else "ready"}


def _allowed(kind, path):
    if kind not in {"legacy_gene", "gene_decision", "paper_observation"}:
        return False
    if not isinstance(path, list) or not path or any(not isinstance(key, (str, int)) or isinstance(key, bool) for key in path):
        return False
    if kind == "paper_observation":
        root, tail = path[0], path[1:]
        allowed = PAPER_FIELDS - {"gene_correspondence", "checklist"}
    elif path[0] in {"papers", "paper_assessments"}:
        if len(path) < 3 or not isinstance(path[1], int) or path[1] < 0:
            return False
        if (kind == "legacy_gene") != (path[0] == "papers"):
            return False
        root, tail = path[2], path[3:]
        allowed = PAPER_FIELDS if path[0] == "papers" else ASSESSMENT_FIELDS
    else:
        root, tail, allowed = path[0], path[1:], GENE_FIELDS
    if root not in allowed:
        return False
    if root == "checklist":
        return len(tail) == 1 and tail[0] in light.CHECKLIST_FIELDS
    if root in {"research_annotations", "paper_checklist_reasons"}:
        return bool(tail) and all(isinstance(key, str) and key for key in tail)
    return not tail


def _prepare_apply(bundle_arg, bundle_sha256, light_index_arg, light_index_sha256, batch_id, results_arg, patch):
    before, records, science = _packet(bundle_arg, bundle_sha256, light_index_arg, light_index_sha256, batch_id, results_arg)
    if not isinstance(patch, dict) or patch.get("contract") != PATCH_CONTRACT or patch.get("base") != before["base"]:
        raise ValueError("patch contract/base hash binding mismatch")
    payload, records = _records(Path(results_arg).resolve(), science)
    operations = patch.get("operations")
    if not isinstance(operations, list) or not operations:
        raise ValueError("patch requires explicit operations")
    before_hashes = {key: _hash(record) for key, record in records.items()}
    touched, paths = set(), []
    roster = _binding(bundle_arg, bundle_sha256, light_index_arg, light_index_sha256, batch_id)[3]
    for operation in operations:
        if not isinstance(operation, dict) or set(operation) - {"op", "record", "path", "expected", "value", "affected_gene_ids"}:
            raise ValueError("unknown patch operation field")
        selector, path = operation.get("record"), operation.get("path")
        if (not isinstance(selector, dict) or set(selector) != {"kind", "id"}
                or any(not isinstance(value, str) for value in selector.values())):
            raise ValueError("operation requires exact record kind/id")
        key = (selector["kind"], selector["id"])
        if key not in records or operation.get("op") != "set" or not _allowed(key[0], path):
            raise ValueError("unknown record, deletion or forbidden key change")
        if "expected" not in operation or "value" not in operation or operation["value"] == MISSING:
            raise ValueError("set operation requires expected and value; deletion is forbidden")
        for previous_key, previous_path in paths:
            if key == previous_key and (path[:len(previous_path)] == previous_path or previous_path[:len(path)] == path):
                raise ValueError("duplicate/conflicting patch operations")
        paths.append((key, path))
        if key[0] == "paper_observation":
            users = sorted({record["interval_gene_id"] for (kind, _), record in records.items() if kind == "gene_decision"
                            and any(assessment.get("observation_id") == key[1] for assessment in _assessments(record) if isinstance(assessment, dict))})
            acknowledgement = operation.get("affected_gene_ids")
            if len(users) > 1 and acknowledgement != users:
                raise ValueError("shared observation change requires exact affected_gene_ids acknowledgement")
        elif key[1] not in roster:
            raise ValueError("cannot patch an unknown roster gene")
        record = records[key]
        container = record
        try:
            for part in path[:-1]:
                container = container[part]
            leaf = path[-1]
            actual = container[leaf] if isinstance(container, list) or leaf in container else MISSING
        except (KeyError, IndexError, TypeError) as error:
            raise ValueError("patch path does not resolve") from error
        if _hash(actual) != _hash(operation["expected"]):
            raise ValueError("patch expected value conflict")
        if not isinstance(container, dict):
            raise ValueError("only explicit object keys may change")
        if leaf == "reference_gene" and str(operation["value"]).casefold() not in {ref.casefold() for ref in roster[key[1]]["reference_endpoints"]}:
            raise ValueError("patch reference is outside frozen roster")
        if leaf == "observation_id" and ("paper_observation", operation["value"]) not in records:
            raise ValueError("patch observation reference is unknown")
        container[leaf] = copy.deepcopy(operation["value"])
        touched.add(key)
    unchanged = {f"{kind}:{key}": value for (kind, key), value in before_hashes.items() if (kind, key) not in touched}
    if any(_hash(records[key]) != value for key, value in before_hashes.items() if key not in touched):
        raise ValueError("non-target record changed")
    return before, payload, science, touched, unchanged


def light_apply(bundle_arg: str, bundle_sha256: str, light_index_arg: str, light_index_sha256: str,
                batch_id: str, results_arg: str, patch_arg: str, patch_sha256: str,
                out_results_arg: str, receipt_arg: str) -> dict[str, Any]:
    source, patch_path = Path(results_arg).resolve(), Path(patch_arg).resolve()
    out, receipt_path = Path(out_results_arg).resolve(), Path(receipt_arg).resolve()
    if out.exists() or receipt_path.exists() or out == receipt_path or out in {source, patch_path} or receipt_path in {source, patch_path}:
        raise ValueError("patch outputs must be distinct new paths")
    if not patch_sha256 or light._sha(patch_path) != patch_sha256:
        raise ValueError("patch SHA-256 mismatch")
    patch = json.loads(patch_path.read_text(encoding="utf-8"))
    before, payload, science, touched, unchanged = _prepare_apply(
        bundle_arg, bundle_sha256, light_index_arg, light_index_sha256, batch_id, results_arg, patch)
    # Exclusive creation preserves the old generation even on accidental reuse.
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("x", encoding="utf-8", newline="\n") as stream:
        if science:
            stream.write(json.dumps(payload, ensure_ascii=False, separators=(",", ":")) + "\n")
        else:
            stream.write("".join(json.dumps(row, ensure_ascii=False, separators=(",", ":")) + "\n" for row in payload))
    after, _, _ = _packet(bundle_arg, bundle_sha256, light_index_arg, light_index_sha256, batch_id, str(out))
    if light._sha(source) != before["base"]["results_sha256"]:
        raise ValueError("original results changed during patch apply")
    receipt = {"contract": PATCH_CONTRACT, "base": before["base"], "patch_sha256": patch_sha256,
               "input_paths": {"bundle": str(Path(bundle_arg).resolve()), "light_index": str(Path(light_index_arg).resolve()),
                               "results": str(source), "patch": str(patch_path)},
               "output": {"path": str(out), "sha256": light._sha(out)}, "original_results_unchanged": True,
               "changed_records": [f"{kind}:{key}" for kind, key in sorted(touched)],
               "unchanged_record_hashes": unchanged, "postpreflight": after,
               "scientific_judgment_owner": "original researcher payload; apply performs no scientific inference"}
    _write_new(receipt_path, receipt)
    return {"results": str(out), "results_sha256": light._sha(out), "receipt": str(receipt_path),
            "receipt_sha256": light._sha(receipt_path),
            "postpreflight": {"counts": after["counts"], "scientific_completion": after["scientific_completion"],
                              "state": "needs_correction" if after["findings"] else "ready", "issue_count": len(after["findings"])}}


def _verified_json(path_arg, expected):
    path, raw = material._verified(path_arg, expected)
    return path, json.loads(raw)


def _exact_keys(value, required, optional=()):
    if not isinstance(value, dict) or set(value) - set(required) - set(optional) or set(required) - set(value):
        raise ValueError("missing or unknown manifest fields")


def _relative_path(parent, value):
    if not isinstance(value, str) or not value.strip():
        raise ValueError("manifest paths must be nonempty strings")
    return str((parent / value).resolve())


def _review_request(path_arg, expected):
    path, request = _verified_json(path_arg, expected)
    _exact_keys(request, ("contract", "binding", "batches"), ("max_chars",))
    if request["contract"] != REVIEW_REQUEST_CONTRACT:
        raise ValueError("review request contract mismatch")
    binding = request["binding"].copy() if isinstance(request["binding"], dict) else request["binding"]
    _exact_keys(binding, ("bundle", "bundle_sha256", "light_index", "light_index_sha256"))
    for key in ("bundle", "light_index"):
        binding[key] = _relative_path(path.parent, binding[key])
    if not isinstance(request["batches"], list) or not request["batches"]:
        raise ValueError("review requires explicit batches")
    max_chars = request.get("max_chars", 12000)
    if not isinstance(max_chars, int) or isinstance(max_chars, bool) or max_chars < 1:
        raise ValueError("max_chars must be a positive integer")
    batches, seen = [], set()
    for row in request["batches"]:
        _exact_keys(row, ("batch_id", "results"), ("sources",))
        bid = row["batch_id"]
        if not isinstance(bid, str) or not bid or bid in seen:
            raise ValueError("batch IDs must be unique nonempty strings")
        seen.add(bid)
        results = _relative_path(path.parent, row["results"])
        packet, _, _ = _packet(*_binding_args(binding), bid, results)
        batches.append({"batch_id": bid, "results": results, "packet": packet, "sources": []})
        sources = row.get("sources", [])
        if not isinstance(sources, list):
            raise ValueError("sources must be an array")
        for selection in sources:
            _exact_keys(selection, ("source", "source_sha256"), ("pointers", "line_start", "line_end", "max_chars"))
            source = _relative_path(path.parent, selection["source"])
            options = {key: value for key, value in selection.items() if key not in {"source", "source_sha256"}}
            source_packet = material.source_read(source, selection["source_sha256"], **options)
            batches[-1]["sources"].append(source_packet)
    return path, binding, batches, max_chars


def _binding_args(binding):
    return tuple(binding[key] for key in ("bundle", "bundle_sha256", "light_index", "light_index_sha256"))


def _batch_counts(batches):
    counts = {state: 0 for state in ("assessed", "not_assessed", "invalid", "retrieval_failed")}
    issues = 0
    for row in batches:
        packet = row["packet"]
        for state in counts:
            counts[state] += packet["counts"][state]
        issues += len(packet["findings"])
    return {"batch_count": len(batches), "counts": counts, "issue_count": issues,
            "scientific_completion": "complete" if all(row["packet"]["scientific_completion"] == "complete" for row in batches) else "partial"}


def light_review_batch(request_arg: str, request_sha256: str, out_arg: str) -> dict[str, Any]:
    """One offline call checks batches and returns only corrections/explicit source selections."""
    out = Path(out_arg).resolve()
    if out.exists():
        raise ValueError("review output directory must be new")
    request_path, binding, batches, max_chars = _review_request(request_arg, request_sha256)
    review_batches, corrections, selected = [], [], []
    for number, row in enumerate(batches, 1):
        packet = row["packet"]
        packet_path = out / f"{number:04d}" / "preflight.json"
        _write_new(packet_path, packet)
        review_batches.append({"batch_id": row["batch_id"], "results": row["results"],
                               "results_sha256": packet["base"]["results_sha256"],
                               "preflight": {"path": str(packet_path), "sha256": light._sha(packet_path)},
                               "sources": row["sources"]})
        if packet["findings"]:
            corrections.append({"batch_id": row["batch_id"], "findings": packet["findings"]})
        for source in row["sources"]:
            item = {"batch_id": row["batch_id"], "source": source["source"], "state": source["state"]}
            if source.get("selections") is not None:
                item["selections"] = source["selections"]
            else:
                item["next_action"] = source.get("next_action", "specify exact source selections")
            selected.append(item)
    summary = _batch_counts(batches)
    review = {"contract": REVIEW_CONTRACT, "request": {"path": str(request_path), "sha256": request_sha256},
              "binding": binding, "batches": review_batches, "summary": summary,
              "scientific_judgment_owner": "original researcher; local review performs no scientific inference"}
    path = out / "review.json"
    _write_new(path, review)
    response = {"review": str(path), "review_sha256": light._sha(path), **summary}
    if corrections:
        corrections_path = out / "corrections.json"
        _write_new(corrections_path, {"batches": corrections})
        response["corrections_artifact"] = {"path": str(corrections_path), "sha256": light._sha(corrections_path),
                                            "batch_count": len(corrections), "finding_count": summary["issue_count"]}
        size = len(json.dumps(corrections, ensure_ascii=False))
        response["corrections"] = corrections if size <= max_chars else {
            "state": "narrow_selection_required", "selected_characters": size,
            "finding_count": summary["issue_count"],
            "next_action": "read exact batch findings from corrections_artifact with light-read-source; none were discarded"}
    if selected:
        size = len(json.dumps(selected, ensure_ascii=False))
        if size > max_chars:
            response["sources"] = {"state": "narrow_selection_required", "selected_characters": size,
                                   "next_action": "request narrower selections or explicitly raise request max_chars"}
        else:
            response["sources"] = selected
    return response


def light_apply_batch(review_arg: str, review_sha256: str, decisions_arg: str,
                      decisions_sha256: str, out_arg: str) -> dict[str, Any]:
    """Apply exact original-researcher operations; keep detailed receipts off context."""
    out = Path(out_arg).resolve()
    if out.exists():
        raise ValueError("apply output directory must be new")
    review_path, review = _verified_json(review_arg, review_sha256)
    if not isinstance(review, dict) or review.get("contract") != REVIEW_CONTRACT:
        raise ValueError("batch review contract mismatch")
    _, binding, batches, _ = _review_request(review["request"]["path"], review["request"]["sha256"])
    if review["binding"] != binding or len(review["batches"]) != len(batches):
        raise ValueError("review request binding mismatch")
    for saved, current in zip(review["batches"], batches):
        _, packet = _verified_json(saved["preflight"]["path"], saved["preflight"]["sha256"])
        if (saved["batch_id"] != current["batch_id"] or saved["results"] != current["results"]
                or saved["results_sha256"] != current["packet"]["base"]["results_sha256"]
                or packet != current["packet"] or saved["sources"] != current["sources"]):
            raise ValueError("review base/source hash binding mismatch")
    decisions_path, decisions = _verified_json(decisions_arg, decisions_sha256)
    _exact_keys(decisions, ("contract", "batches"))
    if decisions["contract"] != DECISIONS_CONTRACT or not isinstance(decisions["batches"], list):
        raise ValueError("batch decisions contract mismatch")
    by_id = {}
    for row in decisions["batches"]:
        _exact_keys(row, ("batch_id",), ("operations", "unchanged"))
        bid = row["batch_id"]
        if not isinstance(bid, str) or bid in by_id:
            raise ValueError("duplicate or invalid decision batch")
        if (set(row) == {"batch_id", "unchanged"} and row["unchanged"] is True):
            by_id[bid] = None
        elif set(row) == {"batch_id", "operations"}:
            if not isinstance(row["operations"], list) or not row["operations"]:
                raise ValueError("operations must be a nonempty array; use explicit unchanged:true")
            by_id[bid] = row["operations"]
        else:
            raise ValueError("each batch requires explicit operations or unchanged:true")
    if set(by_id) != {row["batch_id"] for row in batches}:
        raise ValueError("decisions must identify every reviewed batch exactly; missing or foreign batch")
    patches = {}
    # Check every operation before writing any batch: a later invalid decision
    # must not leave an earlier batch looking successfully applied.
    for row in batches:
        operations = by_id[row["batch_id"]]
        if operations is not None:
            patch = {"contract": PATCH_CONTRACT, "base": row["packet"]["base"], "operations": operations}
            _prepare_apply(*_binding_args(binding), row["batch_id"], row["results"], patch)
            patches[row["batch_id"]] = patch
    applied, changed = [], 0
    for number, row in enumerate(batches, 1):
        bid = row["batch_id"]
        if bid not in patches:
            applied.append({"batch_id": bid, "unchanged": True, "results": row["results"],
                            "results_sha256": row["packet"]["base"]["results_sha256"]})
            continue
        directory = out / f"{number:04d}"
        patch_path = directory / "patch.json"
        _write_new(patch_path, patches[bid])
        result = light_apply(*_binding_args(binding), bid, row["results"], str(patch_path), light._sha(patch_path),
                             str(directory / "results.jsonl"), str(directory / "receipt.json"))
        row["packet"] = _packet(*_binding_args(binding), bid, result["results"])[0]
        applied.append({"batch_id": bid, **result})
        changed += 1
    summary = {**_batch_counts(batches), "changed_batch_count": changed, "unchanged_batch_count": len(batches) - changed}
    receipt = {"contract": BATCH_APPLY_CONTRACT, "review": {"path": str(review_path), "sha256": review_sha256},
               "decisions": {"path": str(decisions_path), "sha256": decisions_sha256}, "batches": applied,
               "summary": summary, "scientific_judgment_owner": "original researcher decisions; no inference or automatic sharing"}
    receipt_path = out / "apply.json"
    _write_new(receipt_path, receipt)
    return {"receipt": str(receipt_path), "receipt_sha256": light._sha(receipt_path), **summary}


def light_postprocess(request_arg: str, request_sha256: str, out_arg: str, *,
                      decisions_arg: str = "", decisions_sha256: str = "",
                      cache_arg: str = "", offline: bool = False) -> dict[str, Any]:
    """Transport saved science through existing checks/conversion without another model turn."""
    out = Path(out_arg).resolve()
    if out.exists():
        raise ValueError("postprocess output directory must be new")
    if bool(decisions_arg) != bool(decisions_sha256):
        raise ValueError("decisions require their exact SHA-256")
    _, binding, batches, _ = _review_request(request_arg, request_sha256)
    _, bundle = light._load_bundle(binding["bundle"], binding["bundle_sha256"])
    original_index = Path(binding["light_index"])
    _, declared = light._validated_index(bundle, binding["bundle_sha256"], original_index, binding["light_index_sha256"])
    if {row["batch_id"] for row in batches} != {row["bid"] for row in declared}:
        raise ValueError("postprocess requires every batch in the bound light index; use a separately frozen scoped index")
    review = light_review_batch(request_arg, request_sha256, str(out / "review"))
    results = {row["batch_id"]: {"path": row["results"], "sha256": row["packet"]["base"]["results_sha256"]}
               for row in batches}
    applied = None
    if decisions_arg:
        applied = light_apply_batch(review["review"], review["review_sha256"], decisions_arg,
                                    decisions_sha256, str(out / "apply"))
        _, apply_receipt = _verified_json(applied["receipt"], applied["receipt_sha256"])
        results = {row["batch_id"]: {"path": row["results"], "sha256": row["results_sha256"]}
                   for row in apply_receipt["batches"]}
    staged_index = out / "input" / "light" / "index.json"
    staged_index.parent.mkdir(parents=True, exist_ok=True)
    staged_index.write_bytes(original_index.read_bytes())
    snapshots = []
    for batch in declared:
        bid = batch["bid"]
        target = staged_index.parent / bid
        target.mkdir()
        for filename in ("task.md", "material.md", "roster.json", "material.full.json"):
            source = original_index.parent / bid / filename
            if filename == "material.full.json" and "material_full_sha256" not in batch:
                continue
            raw = source.read_bytes()
            (target / filename).write_bytes(raw)
        source, raw = material._verified(results[bid]["path"], results[bid]["sha256"])
        staged_results = target / "results.jsonl"
        staged_results.write_bytes(raw)
        snapshots.append({"batch_id": bid, "source": {"path": str(source), "sha256": results[bid]["sha256"]},
                          "snapshot": {"path": str(staged_results), "sha256": light._sha(staged_results)}})
    converted = light.light_convert(binding["bundle"], binding["bundle_sha256"], str(staged_index),
                                    binding["light_index_sha256"], str(out / "converted"),
                                    cache_arg=cache_arg, offline=offline)
    gate = light.verify_light_completion(converted["receipt"], converted["receipt_sha256"])
    # Recheck the original inputs after conversion; snapshots never authorize changed originals.
    _review_request(request_arg, request_sha256)
    for row in batches:
        material._verified(row["results"], row["packet"]["base"]["results_sha256"])
    if decisions_arg:
        material._verified(decisions_arg, decisions_sha256)
    output_paths = {"review": review["review"], "conversion_receipt": converted["receipt"],
                    "agent_responses": converted["agent_responses"],
                    "member_coverage": converted["member_coverage"]["path"]}
    if converted["discoveries"]:
        output_paths["discoveries"] = converted["discoveries"]
    if applied:
        output_paths["apply_receipt"] = applied["receipt"]
    artifacts = {key: {"path": path, "sha256": light._sha(Path(path))} for key, path in output_paths.items()}
    response = {"contract": POSTPROCESS_CONTRACT, "case_id": bundle["case_id"],
                "state": "converted" if gate["state"] == "complete" else "converted_partial",
                "scientific_completion": gate["state"], "member_counts": gate["member_counts"],
                "may_claim_all_genes_assessed": gate["may_claim_all_genes_assessed"],
                "counts": converted["counts"], "artifacts": artifacts,
                "next_action": ("extend/refreeze discoveries before canonical import" if converted["discoveries"]
                                else "pass saved agent_responses through the existing llm-run/import consumer"),
                "additional_model_calls": 0}
    if review.get("corrections_artifact"):
        response["original_findings"] = review["corrections_artifact"]
        if not applied:
            response["correction_action"] = "return only original_findings to the original researcher; preserve healthy records"
    receipt = {**response, "request": {"path": str(Path(request_arg).resolve()), "sha256": request_sha256},
               "binding": binding, "result_snapshots": snapshots,
               "original_results_unchanged": True,
               "scientific_judgment_owner": "original researcher; no inference, new source reading, or model execution"}
    path = out / "postprocess.json"
    _write_new(path, receipt)
    return {**response, "receipt": str(path), "receipt_sha256": light._sha(path)}
