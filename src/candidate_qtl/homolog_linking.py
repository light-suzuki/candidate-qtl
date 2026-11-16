"""Agent-in-the-loop homolog linking module (issue #400 step 2).

Provides provider-neutral prompt export, mechanical validation and import of
agent linking proposals, counterfactual Stage-A rank evaluation, and
memorization probe stratification.

Invariant contracts:
* rank_contribution is ALWAYS "none" in the production artifacts.
* Mechanical gates enforce closed-world candidates, substring grounding of quotes,
  strict upper budgets, and first-class abstention (not_evaluated_no_link).
* Deterministic weights, score policies, and canonical ranking tables remain untouched.
* Memorization-probe stratification is three-valued (issue #419): an answer is
  "match" only through a shared family_labels identity key (curated alias with
  provenance, symbol, normalized label, or stem), "mismatch" only on an explicit
  probe abstention, and any other non-matching relation is
  "not_evaluated_probe_ambiguous" -- memory is never reported as reasoning.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd

from .core import runtime_path
from .family_labels import (
    ALIAS_PROVENANCE,
    family_identity_keys,
    family_identity_keys_expanded,
    family_match_level,
    family_match_level_expanded,
)
from .schemas import require_probe_stratification

SCHEMA_VERSION = "1.0"
RANK_CONTRIBUTION = "none"
REVIEW_REQUIRED = "yes"

PROBE_MATCH = "match"
PROBE_MISMATCH = "mismatch"
PROBE_AMBIGUOUS = "not_evaluated_probe_ambiguous"

# Columns of the per-case memorization_probe_score.tsv artifact. The evidence
# cell is a JSON list of {"term", "matched_key", "match_level", "provenance"}
# entries so the alias used for a match and its source stay machine-readable.
PROBE_SCORE_COLUMNS = [
    "case_id",
    "probe_match_status",
    "stratification",
    "verdict",
    "linking_hit",
    "probe_abstained",
    "probe_match_evidence",
]

VALID_CONFIDENCE_TIERS = frozenset({
    "direct_catalytic_enzyme",
    "pathway_regulator_or_tf",
    "broad_family_member",
    "hypothetical_link",
    "no_link",
})

PROPOSAL_TSV_COLUMNS = [
    "case_id",
    "proposal_id",
    "gene_id",
    "confidence_tier",
    "evidence_quote",
    "mechanism_rationale",
    "review_required",
    "rank_contribution",
    "status",
]


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def extract_homology_candidates(
    evidence_path: str | Path,
    discover_manifest_path: str | Path | None = None,
) -> tuple[str, list[dict[str, str]]]:
    """Extract candidate gene IDs and their homology descriptions.

    Returns (phenotype, list_of_candidates).
    """
    ev_path = runtime_path(evidence_path)
    phenotype = ""
    if discover_manifest_path:
        mf_path = runtime_path(discover_manifest_path)
        if mf_path.is_file():
            mf = json.loads(mf_path.read_text(encoding="utf-8"))
            phenotype = mf.get("inputs", {}).get("phenotype", "")

    edf = pd.read_csv(ev_path, sep="\t", dtype=str).fillna("")
    homology_rows = edf[edf["evidence_code"] == "homology_context"]
    
    candidates: list[dict[str, str]] = []
    seen = set()
    for _, row in homology_rows.iterrows():
        gid = row["gene_id"]
        desc = row["description"]
        if gid not in seen:
            seen.add(gid)
            candidates.append({"gene_id": gid, "homology_description": desc})

    all_genes = edf["gene_id"].unique()
    for gid in all_genes:
        if gid not in seen:
            seen.add(gid)
            candidates.append({"gene_id": gid, "homology_description": "no homology hit found"})

    return phenotype, candidates


def export_homolog_linking_prompts(
    out_dir: str | Path,
    case_id: str,
    phenotype: str,
    evidence_path: str | Path,
    discover_manifest_path: str | Path | None = None,
) -> dict[str, Any]:
    """Export self-contained prompt package for homolog linking."""
    out = runtime_path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    derived_pheno, candidates = extract_homology_candidates(evidence_path, discover_manifest_path)
    final_phenotype = phenotype or derived_pheno
    if not final_phenotype:
        raise ValueError("phenotype text must be provided explicitly or via discover_manifest.json")

    prompt_payload = {
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id,
        "phenotype": final_phenotype,
        "instructions": (
            "You are a plant genetics researcher evaluating candidate genes in a QTL interval for a target phenotype.\n"
            "Evaluate whether any candidate gene's homology description indicates a direct or indirect causal mechanism for the phenotype.\n"
            "Requirements:\n"
            "1. Closed-world: You may ONLY propose gene_ids listed in the candidates universe below.\n"
            "2. Grounding: You MUST provide an exact substring quote from the candidate's homology description as `evidence_quote`.\n"
            "3. Confidence tier: Must be one of ['direct_catalytic_enzyme', 'pathway_regulator_or_tf', 'broad_family_member', 'hypothetical_link', 'no_link'].\n"
            "4. Abstention: If no candidate mechanistically links to the phenotype, set abstain=true and provide empty proposals.\n"
            "5. Budget: Propose at most 3 candidates.\n"
            "6. Alternatives: Explain why other close candidates (if any) were rejected in `rejected_alternatives`."
        ),
        "allowed_confidence_tiers": sorted(VALID_CONFIDENCE_TIERS),
        "candidates": candidates,
        "expected_response_format": {
            "case_id": case_id,
            "abstain": False,
            "proposals": [
                {
                    "gene_id": "<exact_gene_id_from_candidates>",
                    "evidence_quote": "<exact_substring_from_homology_description>",
                    "confidence_tier": "direct_catalytic_enzyme",
                    "mechanism_rationale": "<explanation_of_biochemical_or_biological_link>",
                }
            ],
            "rejected_alternatives": [
                {
                    "gene_id": "<nearby_or_related_gene_id>",
                    "reason": "<why_this_was_not_promoted>",
                }
            ],
        },
    }

    prompt_json_str = json.dumps(prompt_payload, ensure_ascii=False, indent=2) + "\n"
    prompt_json_path = out / "homolog_linking_prompt.json"
    prompt_json_path.write_text(prompt_json_str, encoding="utf-8")

    md_lines = [
        f"# Homolog Linking Task: Case `{case_id}`",
        "",
        f"**Phenotype**: {final_phenotype}",
        "",
        "## Instructions",
        prompt_payload["instructions"],
        "",
        "## Candidate Genes & Homology Descriptions",
        f"Total candidates in interval universe: {len(candidates)}",
        "",
        "| Gene ID | Homology Description |",
        "|---|---|",
    ]
    for c in candidates:
        sanitized_desc = c["homology_description"].replace("|", "\\|")
        md_lines.append(f"| {c['gene_id']} | {sanitized_desc} |")

    prompt_md_path = out / "homolog_linking_prompt.md"
    prompt_md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    prompt_sha256 = _sha256_text(prompt_json_str)

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id,
        "phenotype": final_phenotype,
        "candidate_count": len(candidates),
        "prompt_json": str(prompt_json_path),
        "prompt_md": str(prompt_md_path),
        "prompt_sha256": prompt_sha256,
        "inputs": {
            "evidence": str(runtime_path(evidence_path)),
            "discover_manifest": str(runtime_path(discover_manifest_path)) if discover_manifest_path else None,
        },
    }
    manifest_path = out / "homolog_linking_export_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    return {
        "manifest": str(manifest_path),
        "prompt_json": str(prompt_json_path),
        "prompt_md": str(prompt_md_path),
        "prompt_sha256": prompt_sha256,
        "candidate_count": len(candidates),
    }


def validate_and_import_linking_proposals(
    prompt_export_manifest_or_json: str | Path,
    response_payload_path: str | Path,
    out_dir: str | Path,
    reviewed_by: str,
    payload_source: str = "agent",
    provider: str = "local_subagent",
    model: str = "claude-or-gemini",
    max_proposals: int = 3,
) -> dict[str, Any]:
    """Mechanically validate and import agent linking proposals with strict gates."""
    if not reviewed_by:
        raise ValueError("reviewed_by is required for provenance and independent verification")

    out = runtime_path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    prompt_input_path = runtime_path(prompt_export_manifest_or_json)
    if prompt_input_path.name == "homolog_linking_export_manifest.json":
        exp_mf = json.loads(prompt_input_path.read_text(encoding="utf-8"))
        prompt_json_path = Path(exp_mf["prompt_json"])
    else:
        prompt_json_path = prompt_input_path

    prompt_content = prompt_json_path.read_text(encoding="utf-8")
    prompt_sha256 = _sha256_text(prompt_content)
    prompt_data = json.loads(prompt_content)
    case_id = prompt_data.get("case_id", "")

    candidate_universe = {
        c["gene_id"]: c.get("homology_description", "")
        for c in prompt_data.get("candidates", [])
    }
    if not candidate_universe:
        raise ValueError("Candidate universe in prompt package is empty")

    resp_path = runtime_path(response_payload_path)
    raw_response_text = resp_path.read_text(encoding="utf-8")
    response_sha256 = _sha256_text(raw_response_text)

    cleaned_json = raw_response_text.strip()
    cleaned_json = re.sub(r"^```(?:json)?\s*", "", cleaned_json, flags=re.IGNORECASE)
    cleaned_json = re.sub(r"\s*```$", "", cleaned_json)

    try:
        resp_data = json.loads(cleaned_json)
    except json.JSONDecodeError as err:
        raise ValueError(f"Invalid JSON in response payload: {err}") from err

    resp_case_id = resp_data.get("case_id", "")
    if resp_case_id and resp_case_id != case_id:
        raise ValueError(f"case_id mismatch: prompt had '{case_id}', response had '{resp_case_id}'")

    is_abstain = bool(resp_data.get("abstain", False))
    raw_proposals = resp_data.get("proposals", [])

    if len(raw_proposals) > max_proposals:
        raise ValueError(
            f"Budget exceeded: proposed {len(raw_proposals)} candidates, maximum allowed is {max_proposals}"
        )

    rows: list[dict[str, str]] = []
    validated_proposals: list[dict[str, Any]] = []

    if is_abstain or len(raw_proposals) == 0:
        overall_status = "not_evaluated_no_link"
    else:
        overall_status = "completed"
        for idx, prop in enumerate(raw_proposals, 1):
            gid = str(prop.get("gene_id", "")).strip()
            quote = str(prop.get("evidence_quote", "")).strip()
            tier = str(prop.get("confidence_tier", "")).strip()
            rationale = str(prop.get("mechanism_rationale", "")).strip()

            if gid not in candidate_universe:
                raise ValueError(
                    f"Closed-world violation: proposed gene_id '{gid}' does not exist in candidate universe"
                )

            homology_desc = candidate_universe[gid]
            if not quote:
                raise ValueError(f"Grounding violation: evidence_quote is empty for candidate '{gid}'")
            if quote not in homology_desc:
                raise ValueError(
                    f"Grounding violation: evidence_quote '{quote}' is not a substring of homology description for '{gid}': '{homology_desc}'"
                )

            if tier not in VALID_CONFIDENCE_TIERS:
                raise ValueError(
                    f"Invalid confidence_tier '{tier}'. Must be one of {sorted(VALID_CONFIDENCE_TIERS)}"
                )

            proposal_id = f"hl_{case_id}_{idx:02d}"
            rows.append({
                "case_id": case_id,
                "proposal_id": proposal_id,
                "gene_id": gid,
                "confidence_tier": tier,
                "evidence_quote": quote,
                "mechanism_rationale": rationale,
                "review_required": REVIEW_REQUIRED,
                "rank_contribution": RANK_CONTRIBUTION,
                "status": "validated_proposal",
            })
            validated_proposals.append({
                "proposal_id": proposal_id,
                "gene_id": gid,
                "confidence_tier": tier,
                "evidence_quote": quote,
                "mechanism_rationale": rationale,
            })

    proposals_df = pd.DataFrame(rows, columns=PROPOSAL_TSV_COLUMNS)
    proposals_tsv_path = out / "homolog_linking_proposals.tsv"
    proposals_df.to_csv(proposals_tsv_path, sep="\t", index=False, lineterminator="\n")

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id,
        "status": overall_status,
        "payload_source": payload_source,
        "provider": provider,
        "model": model,
        "reviewed_by": reviewed_by,
        "prompt_sha256": prompt_sha256,
        "response_sha256": response_sha256,
        "abstain": is_abstain,
        "proposal_count": len(rows),
        "max_proposals": max_proposals,
        "rank_contribution": RANK_CONTRIBUTION,
        "proposals_tsv": str(proposals_tsv_path),
        "rejected_alternatives": resp_data.get("rejected_alternatives", []),
        "policy": [
            "rank_contribution is always none: proposals inform review and integrated narrowing only",
            "deterministic ranking weights and scores remain untouched",
            "closed-world and exact substring grounding enforced mechanically",
        ],
    }
    manifest_path = out / "homolog_linking_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    return {
        "manifest": str(manifest_path),
        "proposals_tsv": str(proposals_tsv_path),
        "status": overall_status,
        "proposal_count": len(rows),
        "proposals": validated_proposals,
    }


def calculate_counterfactual_ranks(
    ranked_tsv_path: str | Path,
    linking_proposals_tsv_path: str | Path,
    truth_gene_ids: list[str] | set[str] | str,
    out_tsv_path: str | Path | None = None,
    counterfactual_bonus: float = 0.5,
) -> dict[str, Any]:
    """Calculate what the truth gene rank WOULD BE if linking proposals acted as support (+bonus).

    Leaves the original ranked_candidates.tsv file completely untouched.
    """
    ranked_p = runtime_path(ranked_tsv_path)
    proposals_p = runtime_path(linking_proposals_tsv_path)

    if isinstance(truth_gene_ids, str):
        truth_set = set(g.strip() for g in truth_gene_ids.split(";") if g.strip())
    else:
        truth_set = set(truth_gene_ids)

    rdf = pd.read_csv(ranked_p, sep="\t", dtype=str).fillna("")
    if proposals_p.is_file():
        pdf = pd.read_csv(proposals_p, sep="\t", dtype=str).fillna("")
        proposed_genes = set(pdf["gene_id"].unique()) if not pdf.empty else set()
    else:
        proposed_genes = set()

    rdf["actual_dry_score"] = pd.to_numeric(rdf["dry_score"], errors="coerce").fillna(0.0)
    rdf["actual_screening_rank"] = pd.to_numeric(rdf["screening_rank"], errors="coerce").fillna(999999).astype(int)

    rdf["counterfactual_bonus"] = rdf["gene_id"].apply(lambda g: counterfactual_bonus if g in proposed_genes else 0.0)
    rdf["counterfactual_score"] = rdf["actual_dry_score"] + rdf["counterfactual_bonus"]

    rdf = rdf.sort_values(by=["counterfactual_score", "actual_screening_rank"], ascending=[False, True]).reset_index(drop=True)
    
    ranks = []
    cur_rank = 1
    for i in range(len(rdf)):
        if i > 0 and rdf.loc[i, "counterfactual_score"] < rdf.loc[i - 1, "counterfactual_score"]:
            cur_rank = i + 1
        ranks.append(cur_rank)
    rdf["counterfactual_screening_rank"] = ranks

    truth_rows = rdf[rdf["gene_id"].isin(truth_set)]
    
    if not truth_rows.empty:
        best_truth = truth_rows.iloc[0]
        actual_rank = int(best_truth["actual_screening_rank"])
        cf_rank = int(best_truth["counterfactual_screening_rank"])
        actual_score = float(best_truth["actual_dry_score"])
        cf_score = float(best_truth["counterfactual_score"])
        truth_hit = best_truth["gene_id"] in proposed_genes
    else:
        actual_rank = None
        cf_rank = None
        actual_score = None
        cf_score = None
        truth_hit = False

    result_summary = {
        "truth_genes": sorted(truth_set),
        "proposed_genes": sorted(proposed_genes),
        "truth_proposed": truth_hit,
        "actual_rank": actual_rank,
        "counterfactual_rank": cf_rank,
        "actual_score": actual_score,
        "counterfactual_score": cf_score,
        "stage_a_top20_actual": (actual_rank is not None and actual_rank <= 20),
        "stage_a_top20_counterfactual": (cf_rank is not None and cf_rank <= 20),
    }

    if out_tsv_path:
        out_p = runtime_path(out_tsv_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        export_cols = [
            "gene_id", "actual_dry_score", "counterfactual_score",
            "actual_screening_rank", "counterfactual_screening_rank",
            "counterfactual_bonus",
        ]
        rdf[export_cols].to_csv(out_p, sep="\t", index=False, lineterminator="\n")

    return result_summary


def _key_provenance(key: str, match_level: str) -> str:
    """Source string for a matched family_labels identity key.

    Curated alias keys resolve to their ALIAS_PROVENANCE entry; symbol and
    label keys name the equality tier that produced them. The provenance is
    recorded verbatim in the probe-match evidence so the artifact column can
    be audited without re-deriving the relation.
    """
    if key.startswith("alias:"):
        prov = ALIAS_PROVENANCE.get(key[len("alias:"):], {})
        db = prov.get("database", "family_labels")
        accession = prov.get("accession", "")
        field = prov.get("field", "")
        base = f"{db}:{accession}" if accession else db
        return f"{base} ({field})" if field else base
    if key.startswith("symbol:"):
        if match_level == "stem":
            return "stem symbol equality (family_labels stem key; may conflate paralogs)"
        return "direct symbol equality (family_labels)"
    return "normalized label equality (family_labels)"


def probe_family_match_status(
    probe_data: dict[str, Any],
    truth_family_or_symbol: str,
) -> tuple[str, list[dict[str, str]]]:
    """Classify a memorization-probe answer against the catalog family.

    Returns ``(probe_match_status, evidence)`` where ``probe_match_status`` is
    one of :data:`PROBE_STRATIFICATION_STATUSES`:

    * ``match`` -- the answer shares a known family_labels identity key
      (curated alias with provenance, symbol, normalized label, or stem) with
      the truth family. ``evidence`` records every matched term, the key, the
      match level, and the provenance of the key.
    * ``mismatch`` -- the probe explicitly abstained (``unknown`` / empty
      answer). Only an explicit abstention is a genuine non-recall: without
      alias-table coverage we cannot tell an un-catalogued synonym from a
      real miss, so any other non-match must not fall to recall_negative.
    * ``not_evaluated_probe_ambiguous`` -- the answer is concrete but shares
      no known identity key with the truth family.

    The truth side stays strict (``family_identity_keys``); only the probe's
    free-text family label is matched with the expanded (Route A) matcher, so
    descriptive answers can meet symbol truths without loosening strict keys
    elsewhere.
    """
    raw_symbols = [str(s).strip() for s in probe_data.get("predicted_gene_symbols", []) if str(s).strip()]
    _unknown_tokens = {"unknown", "n/a", "na", "none"}
    symbols = [s for s in raw_symbols if s.casefold() not in _unknown_tokens]
    family_text = str(probe_data.get("predicted_enzyme_or_family", "")).strip()

    truth_keys = family_identity_keys(truth_family_or_symbol)
    if not truth_keys:
        return PROBE_AMBIGUOUS, []

    for sym in symbols:
        level = family_match_level(sym, truth_family_or_symbol)
        if not level:
            continue
        shared = family_identity_keys(sym) & truth_keys
        key = _preferred_shared_key(shared)
        return PROBE_MATCH, [{
            "term": sym,
            "matched_key": key,
            "match_level": level,
            "provenance": _key_provenance(key, level),
        }]

    if family_text and family_text.casefold() not in _unknown_tokens:
        level = family_match_level_expanded(family_text, truth_family_or_symbol)
        if level:
            shared = family_identity_keys_expanded(family_text) & truth_keys
            key = _preferred_shared_key(shared)
            return PROBE_MATCH, [{
                "term": family_text,
                "matched_key": key,
                "match_level": level,
                "provenance": _key_provenance(key, level),
            }]

    family_empty = not family_text or family_text.casefold() in _unknown_tokens
    if not symbols and family_empty:
        return PROBE_MISMATCH, []
    return PROBE_AMBIGUOUS, []


def _preferred_shared_key(shared_keys: set[str]) -> str:
    """Most informative shared key: curated alias, then the most specific key.

    Within the symbol/label prefixes the longest key wins so an exact symbol
    ("symbol:bas1") is reported instead of its stem ("symbol:bas") when both
    are shared; a stem key only appears in the evidence when the match itself
    is stem-level.
    """
    for prefix in ("alias:", "symbol:", "label:"):
        candidates = sorted(
            (key for key in shared_keys if key.startswith(prefix)),
            key=lambda key: (-len(key), key),
        )
        if candidates:
            return candidates[0]
    return sorted(shared_keys)[0]


def export_memorization_probe_prompt(
    out_dir: str | Path,
    case_id: str,
    phenotype: str,
) -> dict[str, Any]:
    """Export control probe prompt with phenotype ONLY (no interval, no homology)."""
    out = runtime_path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    payload = {
        "schema_version": SCHEMA_VERSION,
        "probe_type": "memorization_control",
        "case_id": case_id,
        "phenotype": phenotype,
        "instructions": (
            "Based ONLY on your general pretraining knowledge of plant genetics, name the known gene symbol(s) "
            "or gene family / enzyme responsible for the given phenotype in pea (Pisum sativum) or related legumes.\n"
            "Do NOT access any candidate lists. If you do not know with certainty, state 'unknown'."
        ),
        "expected_response_format": {
            "case_id": case_id,
            "predicted_gene_symbols": ["<e.g. BAS1, CYP75A>"],
            "predicted_enzyme_or_family": "<e.g. Beta-amyrin synthase>",
            "confidence": "high / medium / low / unknown",
        },
    }

    content = json.dumps(payload, ensure_ascii=False, indent=2) + "\n"
    prompt_path = out / "memorization_probe_prompt.json"
    prompt_path.write_text(content, encoding="utf-8")
    prompt_sha256 = _sha256_text(content)

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id,
        "probe_type": "memorization_control",
        "prompt_path": str(prompt_path),
        "prompt_sha256": prompt_sha256,
    }
    mf_path = out / "memorization_probe_manifest.json"
    mf_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    return {"manifest": str(mf_path), "prompt_path": str(prompt_path), "prompt_sha256": prompt_sha256}


def score_linking_and_probe_results(
    linking_manifest_path: str | Path,
    memorization_probe_response_path: str | Path | None,
    truth_gene_ids: list[str] | set[str] | str,
    truth_family_or_symbol: str,
) -> dict[str, Any]:
    """Score linking performance stratified by memorization probe.

    ``probe_match_status`` is three-valued (issue #419): "match" maps to
    ``recall_positive``, "mismatch" to ``recall_negative``, and
    "not_evaluated_probe_ambiguous" to neither -- an answer whose relation to
    the truth family is not covered by the curated alias registry must never
    be reported as reasoning (the error direction would overestimate
    performance). The aliases used for a match and their provenance are
    recorded in ``probe_match_evidence`` and in the
    ``memorization_probe_score.tsv`` artifact written next to the linking
    manifest, so headline numbers can be traced to per-case strata
    mechanically.
    """
    l_mf = json.loads(runtime_path(linking_manifest_path).read_text(encoding="utf-8"))

    if isinstance(truth_gene_ids, str):
        truth_set = set(g.strip() for g in truth_gene_ids.split(";") if g.strip())
    else:
        truth_set = set(truth_gene_ids)

    proposals_tsv = Path(l_mf.get("proposals_tsv", ""))
    linking_hit = False
    proposed_genes = []
    if proposals_tsv.is_file():
        pdf = pd.read_csv(proposals_tsv, sep="\t", dtype=str).fillna("")
        if not pdf.empty:
            proposed_genes = list(pdf["gene_id"].unique())
            linking_hit = bool(set(proposed_genes) & truth_set)

    probe_match_status = PROBE_AMBIGUOUS
    probe_evidence: list[dict[str, str]] = []
    probe_abstained = False
    if memorization_probe_response_path:
        p_path = runtime_path(memorization_probe_response_path)
        if p_path.is_file():
            p_text = p_path.read_text(encoding="utf-8").strip()
            p_text = re.sub(r"^```(?:json)?\s*", "", p_text, flags=re.IGNORECASE)
            p_text = re.sub(r"\s*```$", "", p_text)
            try:
                p_data = json.loads(p_text)
            except json.JSONDecodeError:
                p_data = None
            if isinstance(p_data, dict):
                probe_match_status, probe_evidence = probe_family_match_status(p_data, truth_family_or_symbol)
                probe_abstained = probe_match_status == PROBE_MISMATCH and not probe_evidence
            # A response that is missing or not parseable is NOT an explicit
            # abstention: it stays not_evaluated_probe_ambiguous so an
            # unreadable answer can never become reasoned_link_success.
    require_probe_stratification(probe_match_status)

    if probe_match_status == PROBE_MATCH:
        stratification = "recall_positive"
    elif probe_match_status == PROBE_MISMATCH:
        stratification = "recall_negative"
    else:
        stratification = "not_evaluated_probe_ambiguous"

    if l_mf.get("status") == "not_evaluated_no_link":
        verdict = "not_evaluated_no_link"
    elif probe_match_status == PROBE_MATCH:
        verdict = "not_evaluated_memorization" if linking_hit else "miss_or_unlinked"
    elif probe_match_status == PROBE_MISMATCH:
        verdict = "reasoned_link_success" if linking_hit else "miss_or_unlinked"
    else:
        verdict = "not_evaluated_probe_ambiguous"

    score_row = {
        "case_id": l_mf.get("case_id", ""),
        "probe_match_status": probe_match_status,
        "stratification": stratification,
        "verdict": verdict,
        "linking_hit": str(linking_hit),
        "probe_abstained": str(probe_abstained),
        "probe_match_evidence": json.dumps(probe_evidence, ensure_ascii=False),
    }
    score_path = runtime_path(linking_manifest_path).parent / "memorization_probe_score.tsv"
    pd.DataFrame([score_row], columns=PROBE_SCORE_COLUMNS).to_csv(
        score_path, sep="\t", index=False, lineterminator="\n"
    )

    return {
        "case_id": l_mf.get("case_id", ""),
        "truth_genes": sorted(truth_set),
        "truth_family": truth_family_or_symbol,
        "proposed_genes": proposed_genes,
        "linking_hit": linking_hit,
        "probe_match_status": probe_match_status,
        "probe_match_evidence": probe_evidence,
        "probe_abstained": probe_abstained,
        "stratification": stratification,
        "verdict": verdict,
        "linking_status": l_mf.get("status", ""),
        "probe_score_tsv": str(score_path),
    }
