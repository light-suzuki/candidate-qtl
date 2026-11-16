"""Bounded family/member/copy semantic resolution with a mechanical gate (#505).

This module fills the deliberate seam between deterministic family-to-copy
resolution and production Wet ordering.  It never calls a model itself.  It:

1. selects only families whose deterministic copy state is still unresolved;
2. writes one bounded provider-neutral prompt per eligible family;
3. imports ``llm-run``/agent responses as structured proposals;
4. promotes a proposal only when the proposed candidate/reference mapping is
   independently supported by interval membership plus a strong orthology
   relation already present in the frozen machine artifacts.

A semantic proposal cannot delete a candidate, cannot create a numerical score,
and cannot turn missing evidence into counterevidence.  If verification fails,
the original deterministic ``disambiguation_status`` is preserved byte-for-
meaning and every unresolved copy remains in the set.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import pandas as pd

from .core import file_fingerprint, runtime_path

SCHEMA_VERSION = "1.0"

ELIGIBLE_UNRESOLVED_STATES = frozenset({
    "family_supported_copy_unresolved",
    "paralog_ambiguous",
    "tandem_ambiguous",
    "insufficient_discriminating_evidence",
    "family_resolution_required",
})
DETERMINISTIC_RESOLVED_STATES = frozenset({"preferred_copy", "co_preferred_tie"})
SEMANTIC_RELATIONS = frozenset({
    "same_family",
    "same_subfamily",
    "likely_reference_relative",
    "distinct_subfamily",
    "unresolved",
})
STRONG_ORTHOLOGY = frozenset({
    "1_1",
    "one_to_one",
    "one_one",
    "reciprocal_best_hit",
    "rbh",
})
FORBIDDEN_MODEL_KEYS = frozenset({
    "score", "causal_score", "probability", "causal_probability", "confidence_score",
})

PROPOSAL_COLUMNS = [
    "proposal_id", "prompt_id", "family_id", "semantic_relation",
    "candidate_gene_ids_considered", "proposed_reference_mapping",
    "reasoning_summary", "evidence_to_verify", "competing_explanations_to_keep",
    "review_required", "proposal_status", "provider", "model",
    "input_sha256", "output_sha256", "rank_contribution",
]
VERIFICATION_EXTRA_COLUMNS = [
    "semantic_match_basis", "semantic_verification_status",
    "semantic_reference_gene", "semantic_relation", "semantic_proposal_id",
]
SHORTLIST_COLUMNS = [
    "family_id", "candidate_gene_id", "disambiguation_status", "leading_axes",
    "subfamily", "orthology_class", "first_wet_priority", "rank_contribution", "notes",
]

SYSTEM_PROMPT = """You are resolving a bounded family/member/copy identity ambiguity in a plant-QTL workflow.
The candidate family members and machine facts were frozen before this request.
Do not rank genes numerically and do not assign probabilities/confidence scores.
Do not invent genes, accessions, references, sources, domains, or orthology relations.
Your output is a semantic proposal only; a deterministic verifier will decide whether any proposed mapping is usable.
If the supplied evidence cannot distinguish the copies, return semantic_relation=unresolved and keep competing explanations.
Return JSON only as {\"resolutions\":[{\"family_id\":\"...\",\"semantic_relation\":\"same_family|same_subfamily|likely_reference_relative|distinct_subfamily|unresolved\",\"candidate_gene_ids_considered\":[\"...\"],\"proposed_reference_mapping\":[{\"candidate_gene_id\":\"...\",\"reference_gene\":\"...\"}],\"reasoning_summary\":\"...\",\"evidence_to_verify\":[\"...\"],\"competing_explanations_to_keep\":[\"...\"],\"review_required\":true}]}.
A proposal never deletes a candidate and never itself establishes preferred_copy."""


def _split(value: object) -> list[str]:
    return [part.strip() for part in str(value or "").replace(",", ";").split(";") if part.strip()]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _norm_relation(value: object) -> str:
    text = str(value or "").strip().casefold()
    return re.sub(r"[^a-z0-9]+", "_", text).strip("_")


def _strong_relation(value: object) -> bool:
    return _norm_relation(value) in STRONG_ORTHOLOGY


def _atlas_reference_genes(row: dict[str, object]) -> list[str]:
    value = row.get("orthology_reference_gene", row.get("reference_gene", ""))
    return sorted(set(_split(value)))


def _atlas_relation(row: dict[str, object]) -> str:
    return str(
        row.get("orthology_relationship", "")
        or row.get("relationship", "")
        or row.get("orthology_class", "")
    ).strip()


def _forbidden_keys(value: object) -> set[str]:
    found: set[str] = set()
    if isinstance(value, dict):
        for key, item in value.items():
            if str(key).casefold() in FORBIDDEN_MODEL_KEYS:
                found.add(str(key))
            found |= _forbidden_keys(item)
    elif isinstance(value, list):
        for item in value:
            found |= _forbidden_keys(item)
    return found


def _prompt_family_record(
    family_id: str,
    rows: pd.DataFrame,
    atlas_by_gene: dict[str, dict[str, object]],
    bridge: pd.DataFrame,
) -> dict[str, object]:
    members: list[dict[str, object]] = []
    allowed_refs: set[str] = set()
    bridge_rows = bridge[bridge["family_id"].astype(str).eq(family_id)] if not bridge.empty and "family_id" in bridge.columns else pd.DataFrame()
    for record in rows.to_dict("records"):
        gene_id = str(record["candidate_gene_id"])
        atlas = atlas_by_gene.get(gene_id, {})
        references = _atlas_reference_genes(atlas)
        allowed_refs.update(references)
        gene_bridge = (
            bridge_rows[bridge_rows["gene_id"].astype(str).eq(gene_id)].to_dict("records")
            if not bridge_rows.empty and "gene_id" in bridge_rows.columns else []
        )
        members.append({
            "candidate_gene_id": gene_id,
            "deterministic_disambiguation_status": str(record.get("disambiguation_status", "")),
            "copy_set_status": str(record.get("copy_set_status", "")),
            "orthology_class": str(record.get("orthology_class", "")),
            "subfamily": str(record.get("subfamily", "")),
            "subfamily_vs_target": str(record.get("subfamily_vs_target", "")),
            "domain_vs_target": str(record.get("domain_vs_target", "")),
            "domain_architecture": str(atlas.get("domain_architecture", "")),
            "orthogroup": str(atlas.get("orthogroup", "")),
            "orthology_reference_genes": references,
            "orthology_relationship": _atlas_relation(atlas),
            "description": str(atlas.get("description", atlas.get("annotation", ""))),
            "bridge": [{
                "route": str(item.get("route", "")),
                "mapping_tier": str(item.get("mapping_tier", "")),
                "ambiguity_status": str(item.get("ambiguity_status", "")),
                "provider_status": str(item.get("provider_status", "")),
            } for item in gene_bridge],
        })
    return {
        "family_id": family_id,
        "member_interval_genes": [str(row["candidate_gene_id"]) for row in rows.to_dict("records")],
        "allowed_reference_genes": sorted(allowed_refs),
        "members": members,
    }


def build_family_semantic_request(
    family_copy_arg: str,
    atlas_arg: str,
    out_dir_arg: str,
    *,
    bridge_arg: str | None = None,
    case_id: str = "",
    max_members_per_family: int = 20,
) -> dict[str, object]:
    """Write bounded semantic prompts only for deterministically unresolved families."""
    if max_members_per_family < 2:
        raise ValueError("max_members_per_family must be >= 2")
    copy_path = runtime_path(family_copy_arg)
    atlas_path = runtime_path(atlas_arg)
    if not copy_path.exists() or not atlas_path.exists():
        raise ValueError("family semantic request requires existing family-copy and atlas artifacts")
    copy = pd.read_csv(copy_path, sep="\t", dtype=str).fillna("")
    required_copy = {"family_id", "candidate_gene_id", "disambiguation_status"}
    missing = sorted(required_copy - set(copy.columns))
    if missing:
        raise ValueError(f"family-copy table is missing columns: {', '.join(missing)}")
    atlas = pd.read_csv(atlas_path, sep="\t", dtype=str).fillna("")
    if "gene_id" not in atlas.columns:
        raise ValueError("family semantic request atlas requires gene_id")
    atlas_by_gene = {str(row["gene_id"]): row for row in atlas.to_dict("records")}
    bridge = (
        pd.read_csv(runtime_path(bridge_arg), sep="\t", dtype=str).fillna("")
        if bridge_arg else pd.DataFrame()
    )

    families: list[dict[str, object]] = []
    skipped: list[dict[str, object]] = []
    for family_id, group in copy.groupby("family_id", sort=True):
        family_id = str(family_id)
        statuses = set(group["disambiguation_status"].astype(str))
        if statuses & DETERMINISTIC_RESOLVED_STATES:
            skipped.append({"family_id": family_id, "reason": "deterministic_copy_already_resolved"})
            continue
        if not statuses & ELIGIBLE_UNRESOLVED_STATES:
            skipped.append({"family_id": family_id, "reason": "no_semantic_fallback_eligible_state"})
            continue
        if len(group) > max_members_per_family:
            skipped.append({"family_id": family_id, "reason": "family_exceeds_bounded_member_limit", "member_count": int(len(group))})
            continue
        missing_genes = sorted(set(group["candidate_gene_id"].astype(str)) - set(atlas_by_gene))
        if missing_genes:
            skipped.append({"family_id": family_id, "reason": "atlas_membership_incomplete", "missing_genes": missing_genes})
            continue
        record = _prompt_family_record(family_id, group, atlas_by_gene, bridge)
        record["prompt_id"] = f"family-semantic-{hashlib.sha256((case_id + '|' + family_id).encode('utf-8')).hexdigest()[:16]}"
        families.append(record)

    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    request = {
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id,
        "stage": "family_semantic_resolution_request",
        "cheap_first": True,
        "eligible_family_count": len(families),
        "skipped_family_count": len(skipped),
        "families": families,
        "skipped_families": skipped,
        "policy": [
            "families with deterministic preferred/co-preferred copies never reach the LLM fallback",
            "only unresolved family/member/copy states are prompt-eligible",
            "every prompt is bounded to one family and at most max_members_per_family copies",
            "model output is proposal-only and cannot alter membership or score",
        ],
    }
    request_path = out / "family_semantic_request.json"
    request_path.write_text(json.dumps(request, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    prompts_path = out / "family_semantic_prompts.jsonl"
    prompt_fingerprints: list[str] = []
    with prompts_path.open("w", encoding="utf-8") as handle:
        for family in families:
            prompt = {
                "prompt_id": family["prompt_id"],
                "gene_id": "__family_semantic_resolution__",
                "system": SYSTEM_PROMPT,
                "user": "Resolve only this frozen family/copy ambiguity.\n\n" + json.dumps(family, ensure_ascii=False, sort_keys=True),
                "family_id": family["family_id"],
                "case_id": case_id,
            }
            encoded = json.dumps(prompt, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
            prompt_fingerprints.append(hashlib.sha256(encoded.encode("utf-8")).hexdigest())
            handle.write(json.dumps(prompt, ensure_ascii=False, sort_keys=True) + "\n")

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id,
        "status": "requested" if families else "not_required_or_not_evaluable",
        "cheap_first": True,
        "max_members_per_family": max_members_per_family,
        "eligible_family_count": len(families),
        "skipped_family_count": len(skipped),
        "prompt_count": len(families),
        "prompt_fingerprints": prompt_fingerprints,
        "rank_contribution": "none",
        "inputs": {
            "family_copy": file_fingerprint(copy_path, full_hash=True),
            "atlas": file_fingerprint(atlas_path, full_hash=True),
            "bridge": file_fingerprint(runtime_path(bridge_arg), full_hash=True) if bridge_arg else None,
        },
        "outputs": {
            "request": str(request_path),
            "prompts": str(prompts_path),
            "manifest": str(out / "family_semantic_request_manifest.json"),
        },
        "policy": request["policy"],
    }
    manifest_path = out / "family_semantic_request_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "request": str(request_path),
        "prompts": str(prompts_path),
        "manifest": str(manifest_path),
        "status": manifest["status"],
        "prompt_count": len(families),
    }


def _load_llm_identity(llm_manifest_arg: str | None) -> tuple[str, str]:
    if not llm_manifest_arg:
        return "", ""
    path = runtime_path(llm_manifest_arg)
    if not path.exists():
        raise ValueError(f"LLM manifest does not exist: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("LLM manifest must contain a JSON object")
    return str(payload.get("provider", "")), str(payload.get("model", ""))


def _mapping_pairs(value: object) -> list[tuple[str, str]]:
    if not isinstance(value, list):
        return []
    pairs: list[tuple[str, str]] = []
    for item in value:
        if not isinstance(item, dict):
            continue
        gene = str(item.get("candidate_gene_id", "")).strip()
        reference = str(item.get("reference_gene", "")).strip()
        if gene or reference:
            pairs.append((gene, reference))
    return pairs


def import_family_semantic_responses(
    request_arg: str,
    prompts_arg: str,
    responses_arg: str,
    out_dir_arg: str,
    *,
    llm_manifest_arg: str | None = None,
    case_id: str = "",
) -> dict[str, object]:
    """Validate provider responses into proposal-only rows."""
    request_path = runtime_path(request_arg)
    prompts_path = runtime_path(prompts_arg)
    responses_path = runtime_path(responses_arg)
    for label, path in (("request", request_path), ("prompts", prompts_path), ("responses", responses_path)):
        if not path.exists():
            raise ValueError(f"family semantic {label} does not exist: {path}")
    request = json.loads(request_path.read_text(encoding="utf-8"))
    if not isinstance(request, dict):
        raise ValueError("family semantic request must be a JSON object")
    family_by_prompt = {
        str(record["prompt_id"]): record
        for record in request.get("families", [])
        if isinstance(record, dict) and record.get("prompt_id")
    }
    prompt_ids = {
        str(json.loads(line)["prompt_id"])
        for line in prompts_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    if prompt_ids != set(family_by_prompt):
        raise ValueError("family semantic prompt roster differs from frozen request")

    payload = json.loads(responses_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not isinstance(payload.get("candidates", []), list):
        raise ValueError("family semantic responses require a candidates list")
    provider, model = _load_llm_identity(llm_manifest_arg)
    input_sha = _sha256(prompts_path)
    output_sha = _sha256(responses_path)
    rows: list[dict[str, object]] = []
    seen_prompt_ids: set[str] = set()

    for response in payload.get("candidates", []):
        if not isinstance(response, dict):
            continue
        prompt_id = str(response.get("prompt_id", "")).strip()
        if prompt_id not in family_by_prompt:
            raise ValueError(f"unknown family semantic prompt_id in response: {prompt_id}")
        if prompt_id in seen_prompt_ids:
            raise ValueError(f"duplicate family semantic response for prompt_id: {prompt_id}")
        seen_prompt_ids.add(prompt_id)
        family = family_by_prompt[prompt_id]
        family_id = str(family["family_id"])
        allowed_genes = set(str(value) for value in family.get("member_interval_genes", []))
        allowed_refs = set(str(value) for value in family.get("allowed_reference_genes", []))
        parse_status = str(response.get("parse_status", "")).strip()
        model_json = response.get("model_json", {})
        if parse_status != "json" or not isinstance(model_json, dict):
            rows.append({
                "proposal_id": f"{prompt_id}:provider-gap", "prompt_id": prompt_id, "family_id": family_id,
                "semantic_relation": "unresolved", "candidate_gene_ids_considered": "",
                "proposed_reference_mapping": "", "reasoning_summary": "",
                "evidence_to_verify": "[]", "competing_explanations_to_keep": "[]",
                "review_required": "true", "proposal_status": "not_evaluated_provider_response",
                "provider": provider, "model": model, "input_sha256": input_sha,
                "output_sha256": output_sha, "rank_contribution": "none",
            })
            continue
        forbidden = _forbidden_keys(model_json)
        if forbidden:
            raise ValueError("family semantic response contains forbidden scoring fields: " + ", ".join(sorted(forbidden)))
        resolutions = model_json.get("resolutions", [])
        if not isinstance(resolutions, list):
            raise ValueError("family semantic model_json.resolutions must be a list")
        if not resolutions:
            resolutions = [{"family_id": family_id, "semantic_relation": "unresolved", "review_required": True}]
        for index, resolution in enumerate(resolutions):
            if not isinstance(resolution, dict):
                continue
            if str(resolution.get("family_id", "")).strip() != family_id:
                raise ValueError(f"family semantic response changed frozen family_id for prompt {prompt_id}")
            relation = str(resolution.get("semantic_relation", "unresolved")).strip()
            if relation not in SEMANTIC_RELATIONS:
                raise ValueError(f"unsupported semantic_relation: {relation}")
            considered_raw = resolution.get("candidate_gene_ids_considered", [])
            considered = [str(value).strip() for value in considered_raw] if isinstance(considered_raw, list) else []
            unknown_genes = sorted(set(considered) - allowed_genes)
            pairs = _mapping_pairs(resolution.get("proposed_reference_mapping", []))
            unknown_genes += sorted({gene for gene, _ in pairs if gene and gene not in allowed_genes})
            unknown_refs = sorted({reference for _, reference in pairs if reference and reference not in allowed_refs})
            proposal_status = "proposed"
            if unknown_genes:
                proposal_status = "unsupported_unknown_candidate"
            elif unknown_refs:
                proposal_status = "unsupported_reference_not_supplied"
            mapping_text = ";".join(f"{gene}={reference}" for gene, reference in pairs if gene or reference)
            rows.append({
                "proposal_id": f"{prompt_id}:{index}",
                "prompt_id": prompt_id,
                "family_id": family_id,
                "semantic_relation": relation,
                "candidate_gene_ids_considered": ";".join(sorted(set(considered))),
                "proposed_reference_mapping": mapping_text,
                "reasoning_summary": str(resolution.get("reasoning_summary", "")),
                "evidence_to_verify": json.dumps(resolution.get("evidence_to_verify", []), ensure_ascii=False, sort_keys=True),
                "competing_explanations_to_keep": json.dumps(resolution.get("competing_explanations_to_keep", []), ensure_ascii=False, sort_keys=True),
                "review_required": "true" if bool(resolution.get("review_required", True)) else "false",
                "proposal_status": proposal_status,
                "provider": provider,
                "model": model,
                "input_sha256": input_sha,
                "output_sha256": output_sha,
                "rank_contribution": "none",
            })

    missing_responses = sorted(set(family_by_prompt) - seen_prompt_ids)
    for prompt_id in missing_responses:
        family_id = str(family_by_prompt[prompt_id]["family_id"])
        rows.append({
            "proposal_id": f"{prompt_id}:missing", "prompt_id": prompt_id, "family_id": family_id,
            "semantic_relation": "unresolved", "candidate_gene_ids_considered": "",
            "proposed_reference_mapping": "", "reasoning_summary": "",
            "evidence_to_verify": "[]", "competing_explanations_to_keep": "[]",
            "review_required": "true", "proposal_status": "not_evaluated_missing_response",
            "provider": provider, "model": model, "input_sha256": input_sha,
            "output_sha256": output_sha, "rank_contribution": "none",
        })

    frame = pd.DataFrame(rows, columns=PROPOSAL_COLUMNS)
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    proposals_path = out / "family_semantic_proposals.tsv"
    frame.to_csv(proposals_path, sep="\t", index=False, lineterminator="\n")
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id or str(request.get("case_id", "")),
        "proposal_count": int(len(frame)),
        "provider": provider,
        "model": model,
        "input_sha256": input_sha,
        "output_sha256": output_sha,
        "proposal_status_counts": frame["proposal_status"].value_counts().to_dict() if len(frame) else {},
        "rank_contribution": "none",
        "inputs": {
            "request": file_fingerprint(request_path, full_hash=True),
            "prompts": file_fingerprint(prompts_path, full_hash=True),
            "responses": file_fingerprint(responses_path, full_hash=True),
            "llm_manifest": file_fingerprint(runtime_path(llm_manifest_arg), full_hash=True) if llm_manifest_arg else None,
        },
        "outputs": {"proposals": str(proposals_path), "manifest": str(out / "family_semantic_proposals_manifest.json")},
        "policy": [
            "responses are proposals only and never alter copy status or candidate membership",
            "unknown candidate/reference identifiers are unsupported rather than guessed",
            "numerical causal score/probability fields are rejected fail-closed",
        ],
    }
    manifest_path = out / "family_semantic_proposals_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"proposals": str(proposals_path), "manifest": str(manifest_path), "proposal_count": int(len(frame))}


def _parse_mapping(text: object) -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    for token in str(text or "").split(";"):
        if "=" not in token:
            continue
        gene, reference = token.split("=", 1)
        gene, reference = gene.strip(), reference.strip()
        if gene and reference:
            pairs.append((gene, reference))
    return pairs


def verify_family_semantic_proposals(
    proposals_arg: str,
    family_copy_arg: str,
    atlas_arg: str,
    out_dir_arg: str,
    *,
    base_shortlist_arg: str | None = None,
    case_id: str = "",
) -> dict[str, object]:
    """Promote copy status only when frozen machine facts verify the proposal."""
    proposals_path = runtime_path(proposals_arg)
    copy_path = runtime_path(family_copy_arg)
    atlas_path = runtime_path(atlas_arg)
    proposals = pd.read_csv(proposals_path, sep="\t", dtype=str).fillna("")
    copy = pd.read_csv(copy_path, sep="\t", dtype=str).fillna("")
    atlas = pd.read_csv(atlas_path, sep="\t", dtype=str).fillna("")
    required = {"family_id", "candidate_gene_id", "disambiguation_status"}
    missing = sorted(required - set(copy.columns))
    if missing:
        raise ValueError(f"family-copy table is missing columns: {', '.join(missing)}")
    if "gene_id" not in atlas.columns:
        raise ValueError("family semantic verifier atlas requires gene_id")
    atlas_by_gene = {str(row["gene_id"]): row for row in atlas.to_dict("records")}

    verification: dict[tuple[str, str], dict[str, str]] = {}
    family_verified_genes: dict[str, set[str]] = {}
    for proposal in proposals.to_dict("records"):
        if str(proposal.get("proposal_status", "")) != "proposed":
            continue
        family_id = str(proposal.get("family_id", ""))
        relation = str(proposal.get("semantic_relation", "unresolved"))
        if relation not in {"same_subfamily", "likely_reference_relative"}:
            continue
        family_rows = copy[copy["family_id"].astype(str).eq(family_id)]
        eligible_genes = {
            str(row["candidate_gene_id"])
            for row in family_rows.to_dict("records")
            if str(row.get("disambiguation_status", "")) in ELIGIBLE_UNRESOLVED_STATES
        }
        for gene_id, reference in _parse_mapping(proposal.get("proposed_reference_mapping", "")):
            if gene_id not in eligible_genes:
                continue
            atlas_row = atlas_by_gene.get(gene_id, {})
            refs = set(_atlas_reference_genes(atlas_row))
            relation_value = _atlas_relation(atlas_row)
            # family-to-copy has already normalized orthology_class; use it as
            # a fallback when the B0 atlas stores no relationship text.
            copy_row = family_rows[family_rows["candidate_gene_id"].astype(str).eq(gene_id)]
            copy_orthology = str(copy_row.iloc[0].get("orthology_class", "")) if len(copy_row) else ""
            relation_ok = _strong_relation(relation_value) or _strong_relation(copy_orthology)
            if reference in refs and relation_ok:
                verification[(family_id, gene_id)] = {
                    "semantic_match_basis": "llm_semantic_proposal+verified_reference_orthology",
                    "semantic_verification_status": "verified_copy_mapping",
                    "semantic_reference_gene": reference,
                    "semantic_relation": relation,
                    "semantic_proposal_id": str(proposal.get("proposal_id", "")),
                }
                family_verified_genes.setdefault(family_id, set()).add(gene_id)

    out_rows: list[dict[str, object]] = []
    for row in copy.to_dict("records"):
        family_id = str(row["family_id"])
        gene_id = str(row["candidate_gene_id"])
        updated = dict(row)
        verified_genes = sorted(family_verified_genes.get(family_id, set()))
        info = verification.get((family_id, gene_id))
        if info and len(verified_genes) == 1:
            updated["disambiguation_status"] = "preferred_copy"
        elif info and len(verified_genes) > 1:
            updated["disambiguation_status"] = "co_preferred_tie"
        for column in VERIFICATION_EXTRA_COLUMNS:
            updated[column] = info.get(column, "") if info else ""
        out_rows.append(updated)

    output_columns = list(copy.columns) + [column for column in VERIFICATION_EXTRA_COLUMNS if column not in copy.columns]
    verified_frame = pd.DataFrame(out_rows, columns=output_columns)

    if base_shortlist_arg:
        base_shortlist_path = runtime_path(base_shortlist_arg)
        base_shortlist = pd.read_csv(base_shortlist_path, sep="\t", dtype=str).fillna("")
    else:
        base_shortlist = pd.DataFrame(columns=SHORTLIST_COLUMNS)
    for column in SHORTLIST_COLUMNS:
        if column not in base_shortlist.columns:
            base_shortlist[column] = ""
    shortlist = base_shortlist[SHORTLIST_COLUMNS].copy()
    promoted_families = set(family_verified_genes)
    if promoted_families and len(shortlist):
        shortlist = shortlist[~shortlist["family_id"].astype(str).isin(promoted_families)].copy()
    additions: list[dict[str, object]] = []
    copy_index = {(str(row["family_id"]), str(row["candidate_gene_id"])): row for row in copy.to_dict("records")}
    for family_id, genes in sorted(family_verified_genes.items()):
        ordered = sorted(genes)
        for gene_id in ordered:
            base = copy_index[(family_id, gene_id)]
            additions.append({
                "family_id": family_id,
                "candidate_gene_id": gene_id,
                "disambiguation_status": "preferred_copy" if len(ordered) == 1 else "co_preferred_tie",
                "leading_axes": "semantic_reference_mapping_verified",
                "subfamily": str(base.get("subfamily", "")),
                "orthology_class": str(base.get("orthology_class", "")),
                "first_wet_priority": "first" if len(ordered) == 1 else "co_first",
                "rank_contribution": "none",
                "notes": "LLM/agent proposed the semantic mapping; interval/reference/orthology facts were machine-verified before promotion",
            })
    if additions:
        shortlist = pd.concat([shortlist, pd.DataFrame(additions, columns=SHORTLIST_COLUMNS)], ignore_index=True)
    shortlist = shortlist.drop_duplicates(subset=["family_id", "candidate_gene_id"], keep="last")

    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    verified_path = out / "family_semantic_copy_verification.tsv"
    shortlist_path = out / "family_semantic_copy_shortlist.tsv"
    verified_frame.to_csv(verified_path, sep="\t", index=False, lineterminator="\n")
    shortlist.to_csv(shortlist_path, sep="\t", index=False, lineterminator="\n")
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id,
        "verified_copy_mapping_count": len(verification),
        "promoted_family_count": len(family_verified_genes),
        "preferred_copy_count": int(verified_frame["disambiguation_status"].eq("preferred_copy").sum()),
        "co_preferred_tie_count": int(verified_frame["disambiguation_status"].eq("co_preferred_tie").sum()),
        "rank_contribution": "none",
        "inputs": {
            "proposals": file_fingerprint(proposals_path, full_hash=True),
            "family_copy": file_fingerprint(copy_path, full_hash=True),
            "atlas": file_fingerprint(atlas_path, full_hash=True),
            "base_shortlist": file_fingerprint(runtime_path(base_shortlist_arg), full_hash=True) if base_shortlist_arg else None,
        },
        "outputs": {
            "family_copy": str(verified_path),
            "shortlist": str(shortlist_path),
            "manifest": str(out / "family_semantic_verification_manifest.json"),
        },
        "policy": [
            "semantic proposal alone never changes disambiguation_status",
            "copy promotion requires a supplied reference gene already linked to that interval member plus one-to-one/RBH-like machine orthology",
            "one verified mapping yields preferred_copy; several yield co_preferred_tie",
            "all non-verified rows preserve their original deterministic status",
            "no candidate is deleted and rank_contribution remains none",
        ],
    }
    manifest_path = out / "family_semantic_verification_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "family_copy": str(verified_path),
        "shortlist": str(shortlist_path),
        "manifest": str(manifest_path),
        "verified_copy_mapping_count": len(verification),
        "promoted_family_count": len(family_verified_genes),
    }
