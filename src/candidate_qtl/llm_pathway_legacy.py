"""Phenotype -> pathway hypothesis -> verified bait panel route (issue #163).

The workflow is strong at map-based prioritization *once* a legitimate non-focal
pathway homolog is supplied as bait, and weak at deciding **which** pathway to
look at from phenotype alone. This route adds a verifiable proposal layer in
front of the existing machinery, without replacing any of it:

    phenotype + direction + inheritance + organ/stage + focal species
      -> llm-pathway-prompt   (versioned, provider-neutral input package)
      -> llm-pathway-import   (raw output -> strict schema, raw kept verbatim)
      -> hypothesis-verify    (accession / species / source / sequence gate)
      -> bait-panel-freeze    (immutable panel + FASTA + query manifest)
      -> query-eligibility -> interval-only mapping -> deterministic rank

Safety boundaries, enforced here and covered by tests:

* An LLM proposal is never canonical evidence. Nothing in this module writes an
  ``evidence.tsv`` row, a ``dry_score``, or a Top-eligibility change.
* Model confidence is recorded but never added to a gene score.
* An unverified or model-invented accession never reaches the bait panel.
* Model memory is not a citable source: a proposal is promoted only after the
  accession, its species, and its sequence are re-retrieved from an authoritative
  database and fingerprinted.
* Blind mode must not see QTL-interval candidates. The prompt package is audited
  against the interval candidate ids, and the manifest records the verdict.
* A single-bait panel is review-required by default, while a bait that passed
  identity/source/species/sequence verification remains usable for non-scoring
  deterministic interval projection.
"""

from __future__ import annotations

import hashlib
import json
import re
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from .claim_binding import load_bound_claims, route_a_authorizes
from .core import file_fingerprint, runtime_path, safe_id
from .identifiers import AGI_PATTERN, canonical_accession, is_agi_locus, is_uniprot_accession, is_usable_symbol
from .query_eligibility import PISUM_MARKERS
from .reference_name_pairs import (IDENTITY_DECLARED_ACCESSION, PROJECTABLE_IDENTITIES,
                                   exact_record_candidates, record_species_matches)
from .source_verify import classify_citation
from .uniprot_error_policy import UniProtEntryEvent, capture_uniprot_entry_events

PROMPT_VERSION = "pathway-hypothesis-1.1.0"
SCHEMA_VERSION = "1.1"
PARSER_VERSION = "pathway-proposal-parser-1.0.0"
UNIPROT_ENTRY_ENDPOINT = "https://rest.uniprot.org/uniprotkb/"
UNIPROT_SEARCH_ENDPOINT = "https://rest.uniprot.org/uniprotkb/search"
# How a proposal without a stable id was turned into one, if at all.
ACCESSION_RESOLUTION_STATES = (
    "supplied_by_proposal", "resolved_from_gene_name", "ambiguous_gene_name",
    "unresolved_gene_name", "not_attempted_offline", "no_gene_name_to_resolve",
    "resolved_from_tair_locus", "resolved_from_tair_locus_cache",
    "ambiguous_tair_locus", "unresolved_tair_locus",
    "resolved_from_declared_source_accession", "ambiguous_declared_source_accessions",
    "provider_unavailable",
    "resolved_from_reference_gene_agi",
    "ambiguous_reference_gene_agi",
    "conflicting_reference_gene_agi",
    "unresolved_reference_gene_correspondence",
    "resolved_from_agi_correspondence",
)
EXECUTION_MODES = ("blind", "assisted")

HYPOTHESIS_COLUMNS = [
    "hypothesis_id", "case_id", "pathway_or_process", "gene_family", "reference_gene_name",
    "reference_species", "proposed_stable_id", "expected_role", "mechanistic_link",
    "counterevidence", "confidence", "uncertainty", "proposal_source", "model",
    "prompt_version", "parser_version", "execution_mode", "review_status", "notes",
]
from .relevance_tiers import INITIAL_FIELDS
HYPOTHESIS_COLUMNS += ["source_title", "source_identifier_aliases", *INITIAL_FIELDS]
VALIDATION_COLUMNS = [
    "hypothesis_id", "case_id", "proposed_stable_id", "declared_species", "verification_status",
    "accession_resolution_status", "accession_resolution_query", "accession_resolution_hit_count",
    "accession_exists", "retrieved_species", "species_match", "sequence_retrieved", "sequence_length",
    "sequence_sha256", "identity_sequence_verification_status", "function_source_backed", "function_source_ids", "focal_species_excluded",
    "proposal_literature_verified", "proposal_verified_source_ids", "reference_record_source_ids",
    "bound_claim_authorized", "bound_claim_state", "bound_claim_reason",
    "bound_claim_curation_method",
    "publication_source_provenance", "publication_records_sha256",
    "answer_derived_check", "paralog_alternatives", "source_provider", "source_url", "retrieved_at",
    "source_record_version", "accession_provider_error",
    "rejection_reason", "review_reason",
]
PANEL_COLUMNS = [
    "case_id", "bait_id", "hypothesis_id", "pathway_or_process", "gene_family", "reference_gene_name",
    "reference_species", "expected_role", "sequence_length", "sequence_sha256", "source_provider",
    "source_url", "retrieved_at", "function_source_ids", "rank_contribution", "panel_status",
    "review_status", "accession_resolution_status", "accession_resolution_query",
    "accession_resolution_hit_count", "proposal_verified_source_ids", "source_record_version",
    "publication_source_provenance", "publication_records_sha256",
    "bound_claim_curation_method", "projection_eligibility", "projection_reason",
]
QUERY_MANIFEST_COLUMNS = [
    "case_id", "homolog_query_id", "homolog_source_id", "query_species", "query_genus",
    "sequence_origin", "selection_basis",
]
REVIEW_COLUMNS = ["reason", "target", "required_next_action", "allowed_choices"]

from .paper_reports import LITERATURE_SEARCH_GUIDANCE

SYSTEM_PROMPT = """You propose *competing* pathway hypotheses for a plant QTL phenotype.

You are NOT ranking candidate genes and NOT naming a causal gene. Your only job is
to propose which biological pathways / processes and which gene families are worth
testing, plus literature-supported reference genes that could serve as sequence
baits. The primary approach is to bridge mechanistic and functional evidence
from other species to the focal species. Focal-species direct reports are also
valid primary research inputs in production; neither focal reports nor one
particular reference species are required for candidate generation.
A strict benchmark may impose a non-focal-reference
restriction, but that restriction must be stated by the benchmark route rather
than inferred for every production case.

Investigate primary literature with Web search when available. Your central role
is to bridge species and form hypotheses: a reported phenotype or process in a
reference plant suggests a mechanism/family whose homologs could explain the
focal phenotype. Explain the reported observation, the proposed mechanistic
analogy, and the uncertainty separately. Identical phenotypes and experimental
validation in the focal species are not prerequisites for proposing a candidate
mechanism. Reference identity/sequence and interval correspondence are checked
downstream; homology itself is not proof of transferred biological function.
Use the case's actual trait and biological context, not a species-specific or
example-specific shortlist. No single model species or trait category is required.

Rules:
- Return several competing hypotheses, not one. A single-pathway answer is treated
  as overconfident and is rejected downstream.
- Every proposal must name a gene family and, where you can, one stable database
  accession (UniProt preferred) for a literature-supported reference gene. Do not
  exclude a same-species primary reference merely because it is from the focal
  species; use the explicit strict-benchmark rule when that route is active.
- Do not invent accessions. If you are not confident an accession exists, leave
  proposed_stable_id empty; a missing id is far better than a fabricated one.
- State the mechanistic link to the phenotype, at least one counterevidence or
  alternative explanation, your confidence, and what you are uncertain about.
- Your output is a proposal only. Every accession, species and sequence is
  re-verified against an authoritative database before anything is used."""
SYSTEM_PROMPT += LITERATURE_SEARCH_GUIDANCE

USER_TEMPLATE = """Case: {case_id}
Focal species: {focal_species}
Phenotype: {phenotype}
Direction: {direction}
Inheritance: {inheritance}
Affected organ/tissue: {organ}
Developmental stage: {stage}
Researcher constraints: {constraints}
Known exclusions: {exclusions}
Execution mode: {mode}
{assisted_context}
Return JSON only:
{{"hypotheses": [{{"pathway_or_process": "", "gene_family": "", "reference_gene_name": "",
  "reference_species": "", "proposed_stable_id": "", "expected_role": "",
  "mechanistic_link": "", "counterevidence": "", "confidence": "low|medium|high",
  "uncertainty": ""}}]}}"""


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _frame(rows: list[dict[str, object]], columns: list[str]) -> pd.DataFrame:
    if not rows:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(rows)[columns]


def _write(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


# --------------------------------------------------------------------------- #
# 1. prompt export
# --------------------------------------------------------------------------- #

def llm_pathway_prompt(
    out_dir_arg: str,
    case_id: str,
    focal_species: str,
    phenotype: str,
    direction: str = "unknown",
    inheritance: str = "unknown",
    organ: str = "",
    stage: str = "",
    mode: str = "blind",
    constraints: list[str] | None = None,
    exclusions: list[str] | None = None,
    assisted_context_arg: str | None = None,
    candidates_arg: str | None = None,
) -> dict[str, object]:
    """Export a versioned, provider-neutral prompt package.

    ``candidates_arg`` is used only as a *leak check*: the interval candidate ids
    must not appear anywhere in the prompt. It is never inserted into the prompt.
    """
    if mode not in EXECUTION_MODES:
        raise ValueError(f"mode must be one of: {', '.join(EXECUTION_MODES)}")
    if not case_id or not focal_species or not phenotype:
        raise ValueError("case_id, focal_species and phenotype are required")
    assisted_context = ""
    if assisted_context_arg:
        if mode == "blind":
            raise ValueError("blind mode must not receive researcher context; use --mode assisted")
        assisted_context = runtime_path(assisted_context_arg).read_text(encoding="utf-8").strip()
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)

    user = USER_TEMPLATE.format(
        case_id=case_id, focal_species=focal_species, phenotype=phenotype, direction=direction,
        inheritance=inheritance, organ=organ or "unspecified", stage=stage or "unspecified",
        constraints="; ".join(constraints or []) or "none",
        exclusions="; ".join(exclusions or []) or "none",
        mode=mode,
        assisted_context=(f"Researcher context (assisted mode):\n{assisted_context}\n" if assisted_context else ""),
    )
    prompt_text = f"{SYSTEM_PROMPT}\n\n{user}\n"
    prompt_path = out / "prompt.txt"
    prompt_path.write_text(prompt_text, encoding="utf-8")
    # The JSONL contract is shared with `llm-run`, which requires a gene_id key.
    # This route is not gene-scoped, so the field carries the case scope instead.
    record = {"prompt_id": f"pathway::{safe_id(case_id)}", "gene_id": f"pathway_hypothesis:{case_id}",
              "system": SYSTEM_PROMPT, "user": user}
    prompts_path = out / "llm_pathway_prompts.jsonl"
    prompts_path.write_text(json.dumps(record, ensure_ascii=False) + "\n", encoding="utf-8")

    leak = _blind_leak_audit(prompt_text, candidates_arg) if mode == "blind" else {
        "checked": False, "reason": "assisted mode may include researcher context", "leaked_identifiers": []}
    if mode == "blind" and leak["leaked_identifiers"]:
        raise ValueError("blind prompt contains QTL interval candidate identifier(s): "
                         + ", ".join(str(value) for value in leak["leaked_identifiers"][:5]))

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id,
        "execution_mode": mode,
        "focal_species": focal_species,
        "input": {"phenotype": phenotype, "direction": direction, "inheritance": inheritance,
                  "organ": organ, "developmental_stage": stage,
                  "researcher_constraints": constraints or [], "known_exclusions": exclusions or []},
        "prompt_version": PROMPT_VERSION,
        "prompt_sha256": hashlib.sha256(prompt_text.encode("utf-8")).hexdigest(),
        "prompt": str(prompt_path),
        "prompts_jsonl": str(prompts_path),
        "assisted_context_supplied": bool(assisted_context),
        "qtl_candidate_information_in_model_input": False,
        "blind_leak_audit": leak,
        "external_transmission": "none until the prompt is executed; `llm-run` records the endpoint it used",
        "responsibility": "the model proposes pathways, families and reference baits only; it never ranks interval candidates",
        "created_at": _now(),
    }
    manifest_path = out / "llm_input_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), "prompt": str(prompt_path), "prompts_jsonl": str(prompts_path),
            "execution_mode": mode, "prompt_version": PROMPT_VERSION,
            "blind_leak_checked": bool(leak["checked"])}


def _blind_leak_audit(prompt_text: str, candidates_arg: str | None) -> dict[str, object]:
    if not candidates_arg:
        return {"checked": False, "reason": "no candidate table supplied for the leak check",
                "leaked_identifiers": [], "candidate_count": 0}
    frame = pd.read_csv(runtime_path(candidates_arg), sep="\t", dtype=str).fillna("")
    if "gene_id" not in frame.columns:
        raise ValueError("candidate table for the blind leak check requires gene_id")
    identifiers = sorted({value for value in frame["gene_id"].astype(str) if value})
    lowered = prompt_text.casefold()
    leaked = [value for value in identifiers if value.casefold() in lowered]
    return {"checked": True, "reason": "interval candidate ids must not appear in a blind prompt",
            "leaked_identifiers": leaked, "candidate_count": len(identifiers)}


# --------------------------------------------------------------------------- #
# 2. raw output import
# --------------------------------------------------------------------------- #

def _extract_json(text: str) -> object | None:
    cleaned = re.sub(r"\x1b\[[0-?]*[ -/]*[@-~]", "", text).strip()
    cleaned = re.sub(r"^```(?:json)?\s*", "", cleaned, flags=re.IGNORECASE)
    cleaned = re.sub(r"\s*```$", "", cleaned)
    try:
        return json.loads(cleaned)
    except json.JSONDecodeError:
        pass
    start, end = cleaned.find("{"), cleaned.rfind("}")
    if start >= 0 and end > start:
        try:
            return json.loads(cleaned[start:end + 1])
        except json.JSONDecodeError:
            return None
    return None


def _proposal_records(payload: object) -> list[dict[str, object]]:
    if isinstance(payload, dict):
        for key in ("hypotheses", "proposals", "results"):
            value = payload.get(key)
            if isinstance(value, list):
                return [item for item in value if isinstance(item, dict)]
        # `llm-run` wraps each prompt response; unwrap one level when present.
        for key in ("response", "parsed", "output"):
            value = payload.get(key)
            if isinstance(value, (dict, list)):
                nested = _proposal_records(value)
                if nested:
                    return nested
        return []
    if isinstance(payload, list):
        records: list[dict[str, object]] = []
        for item in payload:
            records.extend(_proposal_records(item))
        return records
    return []


def llm_pathway_import(
    raw_arg: str,
    out_dir_arg: str,
    input_manifest_arg: str | None = None,
    model: str = "",
    proposal_source: str = "",
    case_id: str = "",
    execution_mode: str = "",
) -> dict[str, object]:
    """Normalize raw model output into the strict proposal schema.

    The raw payload is copied verbatim next to the parsed table so both survive;
    a manually written proposal file can be imported through the same path.
    """
    raw_path = runtime_path(raw_arg)
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    raw_text = raw_path.read_text(encoding="utf-8")
    raw_copy = out / f"raw_model_output{raw_path.suffix or '.txt'}"
    raw_copy.write_text(raw_text, encoding="utf-8")

    manifest_input: dict[str, object] = {}
    if input_manifest_arg:
        manifest_input = json.loads(runtime_path(input_manifest_arg).read_text(encoding="utf-8"))
    case_id = case_id or str(manifest_input.get("case_id", ""))
    execution_mode = execution_mode or str(manifest_input.get("execution_mode", "")) or "blind"
    if execution_mode not in EXECUTION_MODES:
        raise ValueError(f"execution_mode must be one of: {', '.join(EXECUTION_MODES)}")

    payloads: list[object] = []
    if raw_path.suffix == ".jsonl":
        for line in raw_text.splitlines():
            if line.strip():
                payloads.append(json.loads(line))
    else:
        parsed = _extract_json(raw_text)
        if parsed is not None:
            payloads.append(parsed)
    proposals: list[dict[str, object]] = []
    for payload in payloads:
        proposals.extend(_proposal_records(payload))
        if isinstance(payload, dict):
            for key in ("text", "content", "output_text"):
                value = payload.get(key)
                if isinstance(value, str):
                    nested = _extract_json(value)
                    if nested is not None:
                        proposals.extend(_proposal_records(nested))

    review_rows: list[dict[str, str]] = []
    rows: list[dict[str, object]] = []
    for index, proposal in enumerate(proposals, 1):
        hypothesis_id = f"hyp_{index:03d}"
        pathway = str(proposal.get("pathway_or_process", "")).strip()
        family = str(proposal.get("gene_family", "")).strip()
        accession = str(proposal.get("proposed_stable_id", "")).strip()
        species = str(proposal.get("reference_species", "")).strip()
        missing = [name for name, value in (("pathway_or_process", pathway), ("gene_family", family)) if not value]
        notes = ""
        if missing:
            notes = "unknown: " + ", ".join(missing)
            review_rows.append({"reason": "proposal_missing_required_field", "target": hypothesis_id,
                                "required_next_action": f"supply {', '.join(missing)} or discard the proposal",
                                "allowed_choices": "re-run the prompt | complete the proposal manually | discard"})
        if not accession:
            review_rows.append({"reason": "proposal_without_stable_accession", "target": hypothesis_id,
                                "required_next_action": "supply a stable accession for this reference gene, or keep the proposal unverifiable",
                                "allowed_choices": "add accession | keep as unverifiable proposal"})
        from .paper_reports import imported_paper_fields, normalize_source_aliases
        rows.append({
            "hypothesis_id": hypothesis_id, "case_id": case_id, "pathway_or_process": pathway,
            **{key: str(proposal.get(key, "")).strip() for key in INITIAL_FIELDS},
            **imported_paper_fields(proposal),
            "source_title": str(proposal.get("source_title", "")).strip(),
            "source_identifier_aliases": normalize_source_aliases(proposal.get("source_identifier_aliases")),
            "gene_family": family, "reference_gene_name": str(proposal.get("reference_gene_name", "")).strip(),
            "reference_species": species, "proposed_stable_id": accession,
            "expected_role": str(proposal.get("expected_role", "")).strip(),
            "mechanistic_link": str(proposal.get("mechanistic_link", proposal.get("rationale", ""))).strip(),
            "counterevidence": str(proposal.get("counterevidence", "")).strip(),
            "confidence": str(proposal.get("confidence", "unknown")).strip() or "unknown",
            "uncertainty": str(proposal.get("uncertainty", "")).strip(),
            "proposal_source": proposal_source or str(manifest_input.get("proposal_source", "")) or "llm",
            "model": model or str(manifest_input.get("model", "")),
            "prompt_version": str(manifest_input.get("prompt_version", PROMPT_VERSION)),
            "parser_version": PARSER_VERSION, "execution_mode": execution_mode,
            "review_status": "pending_verification", "notes": notes,
        })

    if not rows:
        review_rows.append({"reason": "no_parsable_proposal", "target": str(raw_copy),
                            "required_next_action": "the raw output produced no schema-conformant proposal; re-run or write proposals manually",
                            "allowed_choices": "re-run the prompt | write proposals manually | abstain"})
    distinct_pathways = len({str(row["pathway_or_process"]).casefold() for row in rows if row["pathway_or_process"]})
    if rows and distinct_pathways < 2:
        review_rows.append({"reason": "single_pathway_overconfidence", "target": case_id or "case",
                            "required_next_action": "a single-pathway proposal set is treated as overconfident; obtain competing hypotheses",
                            "allowed_choices": "re-run for competing hypotheses | accept explicitly as a single-hypothesis run"})

    optional_initial = {"source_title", "source_identifier_aliases", *INITIAL_FIELDS}
    hypotheses = _frame(rows, [column for column in HYPOTHESIS_COLUMNS
        if column not in optional_initial or any(row.get(column, "") for row in rows)])
    reviews = _frame(review_rows, REVIEW_COLUMNS).drop_duplicates().sort_values(["reason", "target"], kind="stable")
    hypotheses_path = out / "llm_pathway_hypotheses.tsv"
    _write(hypotheses, hypotheses_path)
    _write(reviews, out / "review_required.tsv")
    manifest = {
        "schema_version": SCHEMA_VERSION, "case_id": case_id, "execution_mode": execution_mode,
        "parser_version": PARSER_VERSION, "prompt_version": str(manifest_input.get("prompt_version", PROMPT_VERSION)),
        "model": model or str(manifest_input.get("model", "")),
        "proposal_source": proposal_source or "llm",
        "raw_output": str(raw_copy), "raw_output_sha256": hashlib.sha256(raw_text.encode("utf-8")).hexdigest(),
        "input_manifest": str(runtime_path(input_manifest_arg)) if input_manifest_arg else None,
        "proposal_count": len(rows), "distinct_pathways": distinct_pathways,
        "review_required_count": int(len(reviews)),
        "outputs": {"hypotheses": str(hypotheses_path), "review_required": str(out / "review_required.tsv")},
        "safety": [
            "every proposal is review_status=pending_verification until the accession, species and sequence are verified",
            "raw and parsed output are both retained; the parser never repairs a proposal silently",
            "no proposal is written to canonical evidence and no proposal changes a gene score",
        ],
    }
    manifest_path = out / "llm_pathway_import_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), "hypotheses": str(hypotheses_path),
            "raw_output": str(raw_copy), "proposal_count": len(rows),
            "distinct_pathways": distinct_pathways, "review_required": str(out / "review_required.tsv")}


# --------------------------------------------------------------------------- #
# 3. verification gate
# --------------------------------------------------------------------------- #

def _load_record_cache(path_arg: str | None) -> dict[str, dict[str, object]]:
    """Local UniProt-shaped records, keyed by accession (offline verification)."""
    if not path_arg:
        return {}
    path = runtime_path(path_arg)
    records: dict[str, dict[str, object]] = {}
    files = sorted(path.glob("*.json")) if path.is_dir() else [path]
    for item in files:
        payload = json.loads(item.read_text(encoding="utf-8"))
        entries = payload.get("results", payload) if isinstance(payload, dict) else payload
        for entry in entries if isinstance(entries, list) else [entries]:
            if isinstance(entry, dict) and entry.get("primaryAccession"):
                records[str(entry["primaryAccession"])] = entry
    return records


def _cached_tair_index(records: dict[str, dict[str, object]]) -> dict[str, list[str]]:
    """Index exact TAIR cross-references from frozen UniProt records.

    A cache is an authoritative record snapshot, not an alias guess.  Only
    explicit ``TAIR`` cross-references with canonical AGI loci are accepted,
    and all matching accessions are retained so callers can fail closed on an
    ambiguous cache rather than selecting one.
    """
    index: dict[str, set[str]] = {}
    for accession, record in records.items():
        cross_references = record.get("uniProtKBCrossReferences", [])
        for cross_reference in cross_references if isinstance(cross_references, list) else []:
            if not isinstance(cross_reference, dict):
                continue
            if str(cross_reference.get("database", "")).casefold() != "tair":
                continue
            locus = canonical_accession(str(cross_reference.get("id", "")))
            if is_agi_locus(locus):
                index.setdefault(locus, set()).add(accession)
    return {locus: sorted(accessions) for locus, accessions in index.items()}


def _searchable_symbol(gene_name: str) -> str:
    """Pull the gene symbol out of a free-text reference-gene field.

    A responder writes what a person writes: "TFL1 (TERMINAL FLOWER 1)",
    "IPT (isopentenyltransferase) family". The symbol is the leading token; the
    parenthetical is a gloss. Anything that is not usable as a symbol -- an
    accession, a locus tag, a digit-led fragment -- yields nothing, because
    searching on it would either fail or match the wrong record.
    """
    for token in re.findall(r"[A-Za-z][A-Za-z0-9_.-]*", str(gene_name or "")):
        if is_usable_symbol(token):
            return token
    return ""


_EXPLICIT_AGI_TOKEN = re.compile(
    rf"(?<![A-Za-z0-9_.-])(?P<locus>{AGI_PATTERN})(?![A-Za-z0-9_.-])", re.IGNORECASE
)


def _explicit_agi_loci(gene_name: str) -> list[str]:
    """Return distinct, explicitly delimited Arabidopsis loci in a name field.

    This is intentionally narrower than symbol extraction: an AGI must be a
    complete token, and only its case is normalized.  Thus a decorated name
    can provide a deterministic TAIR lookup without treating an embedded or
    malformed identifier as evidence.
    """
    return sorted({match.group("locus").upper() for match in _EXPLICIT_AGI_TOKEN.finditer(
        str(gene_name or ""))})


def _resolve_accession(gene_name: str, species: str, timeout: int = 60) -> tuple[str, str, str, int]:
    """Resolve a gene symbol + species to one reviewed UniProt accession.

    The investigation prompt tells the responder not to invent accessions and to
    leave the field empty instead, promising that an empty field is recoverable.
    Nothing recovered it: a proposal naming a real, correct gene was rejected
    ``no_stable_accession_proposed`` while the accession sat one query away on the
    same API this stage is already permitted to call (issue #492).

    Resolution is not an exemption. The accession returned here goes on to the
    same existence / species / sequence / source checks as one the proposal
    supplied itself. A symbol matching several reviewed entries is NOT picked
    between -- that is a question for a reviewer, not a tie-break.

    Returns ``(accession, resolution_status, query, hit_count)``.
    """
    symbol = _searchable_symbol(gene_name)
    if not symbol or not str(species or "").strip():
        return "", "no_gene_name_to_resolve", "", 0
    query = f'gene:{symbol} AND organism_name:"{species.strip()}" AND reviewed:true'
    url = UNIPROT_SEARCH_ENDPOINT + "?" + urllib.parse.urlencode(
        {"query": query, "format": "json", "size": 5, "fields": "accession"})
    request = urllib.request.Request(url, headers={"User-Agent": "candidate-qtl hypothesis-verify"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as error:
        raise RuntimeError(f"UniProt API HTTP {error.code} for gene-name query {query!r}") from error
    except urllib.error.URLError as error:
        raise RuntimeError(f"UniProt API connection failed for {query!r}: {error.reason}") from error
    if (not isinstance(payload, dict) or "results" not in payload
            or not isinstance(payload.get("results"), list)):
        raise ValueError(f"UniProt API returned an invalid search payload for {query!r}")
    results = payload.get("results", [])
    accessions = [str(item.get("primaryAccession", "")).strip()
                  for item in results if isinstance(item, dict) and item.get("primaryAccession")]
    if len(accessions) == 1:
        return accessions[0], "resolved_from_gene_name", query, 1
    if len(accessions) > 1:
        return "", "ambiguous_gene_name", query, len(accessions)
    return "", "unresolved_gene_name", query, 0


def _resolve_accession_candidates(gene_name: str, species: str, timeout: int = 60, *,
                                  query_trace: dict | None = None) -> tuple[list[str], str, str, int]:
    """Return every accession from the exact gene/species provider query.

    This is the non-selecting companion to ``_resolve_accession``.  It keeps
    the provider query and ordering deterministic while allowing Route A to
    carry all exact alternatives downstream for independent verification.
    """
    # Explicit whole AGI endpoints are identifiers, not general usable symbols.
    locus = str(gene_name or "").strip().upper()
    symbol = locus if is_agi_locus(locus) else _searchable_symbol(gene_name)
    if not symbol or not str(species or "").strip():
        return [], "no_gene_name_to_resolve", "", 0
    query = f'gene:{symbol} AND organism_name:"{species.strip()}"'
    if query_trace is not None:
        query_trace["query"] = query
        query_trace.setdefault("query_urls", [])
        query_trace.setdefault("page_accessions", [])
        query_trace.setdefault("accessions", [])
    accessions: set[str] = set()
    cursor = "*"
    seen_cursors: set[str] = set()
    while cursor:
        if cursor in seen_cursors:
            raise RuntimeError(f"UniProt API pagination repeated cursor for gene-name query {query!r}")
        seen_cursors.add(cursor)
        params = {"query": query, "format": "json", "size": 500}
        if cursor != "*":
            params["cursor"] = cursor
        url = UNIPROT_SEARCH_ENDPOINT + "?" + urllib.parse.urlencode(params)
        if query_trace is not None:
            query_trace["query_urls"].append(url)
        request = urllib.request.Request(url, headers={"User-Agent": "candidate-qtl hypothesis-verify"})
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                payload = json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as error:
            raise RuntimeError(f"UniProt API HTTP {error.code} for gene-name query {query!r}") from error
        except urllib.error.URLError as error:
            raise RuntimeError(f"UniProt API connection failed for {query!r}") from error
        if not isinstance(payload, dict) or not isinstance(payload.get("results"), list):
            raise ValueError(f"UniProt API returned an invalid search payload for {query!r}")
        page_accessions: set[str] = set()
        for item in payload.get("results", []):
            if not isinstance(item, dict):
                raise ValueError(f"UniProt API returned a malformed result for {query!r}")
            organism = item.get("organism")
            if not isinstance(organism, dict):
                raise ValueError(f"UniProt API returned a malformed organism for {query!r}")
            accession = item.get("primaryAccession")
            if isinstance(accession, str) and accession:
                page_accessions.update(exact_record_candidates({accession: item}, symbol, species))
        accessions.update(page_accessions)
        if query_trace is not None:
            query_trace["page_accessions"].append(sorted(page_accessions))
            query_trace["accessions"] = sorted(accessions)
        next_link = ""
        link_header = response.headers.get("Link", "")
        for link in link_header.split(","):
            if 'rel="next"' in link:
                next_link = link.split("<", 1)[1].split(">", 1)[0] if "<" in link and ">" in link else ""
                if not next_link:
                    raise ValueError("UniProt API returned a malformed next-page link")
                break
        cursor = urllib.parse.parse_qs(urllib.parse.urlparse(next_link).query).get("cursor", [""])[0] if next_link else ""
        if next_link and not cursor:
            raise ValueError("UniProt API next-page link has no cursor; search is incomplete")
    accessions = sorted(accessions)
    return accessions, ("resolved_from_gene_name" if len(accessions) == 1 else
                        "ambiguous_gene_name" if accessions else "unresolved_gene_name"), query, len(accessions)


def _resolve_tair_locus(locus: str, species: str, timeout: int = 60) -> tuple[str, str, str, int]:
    """Resolve one supplied TAIR/AGI locus to a reviewed UniProt accession.

    ``phenotype-investigate`` explicitly permits TAIR stable identifiers, while
    the verification provider is UniProt.  A TAIR locus must therefore be
    translated through UniProt's TAIR cross-reference index before the normal
    record, species, sequence and source gates run.  Multiple hits are never
    tie-broken.
    """
    normalised = canonical_accession(locus)
    if not is_agi_locus(normalised) or not str(species or "").strip():
        return "", "unresolved_tair_locus", "", 0
    query = f'xref:tair-{normalised} AND organism_name:"{species.strip()}" AND reviewed:true'
    url = UNIPROT_SEARCH_ENDPOINT + "?" + urllib.parse.urlencode(
        {"query": query, "format": "json", "size": 5, "fields": "accession"})
    request = urllib.request.Request(url, headers={"User-Agent": "candidate-qtl hypothesis-verify"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as error:
        raise RuntimeError(f"UniProt API HTTP {error.code} for TAIR query {query!r}") from error
    except urllib.error.URLError as error:
        raise RuntimeError(f"UniProt API connection failed for {query!r}: {error.reason}") from error
    if (not isinstance(payload, dict) or "results" not in payload
            or not isinstance(payload.get("results"), list)):
        raise ValueError(f"UniProt API returned an invalid TAIR search payload for {query!r}")
    results = payload.get("results", [])
    accessions = [str(item.get("primaryAccession", "")).strip()
                  for item in results if isinstance(item, dict) and item.get("primaryAccession")]
    if len(accessions) == 1:
        return accessions[0], "resolved_from_tair_locus", query, 1
    if len(accessions) > 1:
        return "", "ambiguous_tair_locus", query, len(accessions)
    return "", "unresolved_tair_locus", query, 0


def _agi_corresponds_to_accession(
    locus: str, accession: str, species: str, *, allow_external_api: bool,
    cached_tair: dict[str, list[str]],
) -> tuple[str, str, int, str]:
    """Prove an AGI/accession correspondence without selecting either side."""
    normalised = canonical_accession(locus).upper()
    cached_accessions = cached_tair.get(normalised, [])
    if len(cached_accessions) == 1 and accession.casefold() == cached_accessions[0].casefold():
        return "resolved_from_agi_correspondence", f"record_cache:xref:tair-{normalised}", 1, ""
    if len(cached_accessions) > 1:
        return "unresolved_reference_gene_correspondence", f"record_cache:xref:tair-{normalised}", len(cached_accessions), ""
    if not allow_external_api:
        return "unresolved_reference_gene_correspondence", f"record_cache:xref:tair-{normalised}", 0, ""
    try:
        resolved, status, query, hits = _resolve_tair_locus(normalised, species)
    except (RuntimeError, TimeoutError, OSError, ValueError) as error:
        return "provider_unavailable", "", 0, str(error)
    if status == "resolved_from_tair_locus" and resolved.casefold() == accession.casefold():
        return "resolved_from_agi_correspondence", query, hits, ""
    return "unresolved_reference_gene_correspondence", query, hits, ""


def _fetch_uniprot(accession: str, timeout: int = 60) -> dict[str, object] | None:
    url = UNIPROT_ENTRY_ENDPOINT + urllib.parse.quote(accession) + ".json"
    request = urllib.request.Request(url, headers={"User-Agent": "candidate-qtl hypothesis-verify"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = json.loads(response.read().decode("utf-8"))
        if (not isinstance(payload, dict)
                or not isinstance(payload.get("primaryAccession"), str)
                or not isinstance(payload.get("organism"), dict)
                or not isinstance(payload.get("sequence"), dict)):
            raise ValueError(f"UniProt API returned an invalid record payload for {accession}")
        return payload
    except urllib.error.HTTPError as error:
        if error.code == 404:
            return None
        raise RuntimeError(f"UniProt API HTTP {error.code} for {accession}") from error
    except urllib.error.URLError as error:
        raise RuntimeError(f"UniProt API connection failed for {accession}: {error.reason}") from error


def _record_species(record: dict[str, object]) -> str:
    organism = record.get("organism", {})
    return str(organism.get("scientificName", "")) if isinstance(organism, dict) else ""


def _record_sequence(record: dict[str, object]) -> str:
    sequence = record.get("sequence", {})
    return str(sequence.get("value", "")) if isinstance(sequence, dict) else ""


def _record_version(record: dict[str, object]) -> str:
    audit = record.get("entryAudit", {})
    if isinstance(audit, dict) and audit.get("entryVersion") is not None:
        return str(audit["entryVersion"])
    return str(record.get("candidate_qtl_record_version", ""))


def _record_source_ids(record: dict[str, object]) -> list[str]:
    ids: set[str] = set()
    for reference in record.get("references", []) if isinstance(record.get("references"), list) else []:
        citation = reference.get("citation", {}) if isinstance(reference, dict) else {}
        for cross in citation.get("citationCrossReferences", []) if isinstance(citation.get("citationCrossReferences"), list) else []:
            if isinstance(cross, dict) and str(cross.get("database", "")).casefold() == "pubmed":
                ids.add(f"PMID:{cross.get('id')}")
    if str(record.get("entryType", "")).casefold().startswith("uniprotkb reviewed"):
        ids.add(f"UniProt:{record.get('primaryAccession')}")
    return sorted(ids)


def _record_gene_aliases(record: dict[str, object]) -> set[str]:
    """Provider-declared gene symbols for direct (non-resolver) review use."""
    aliases: set[str] = set()
    genes = record.get("genes", [])
    for gene in genes if isinstance(genes, list) else []:
        if not isinstance(gene, dict):
            continue
        for key in ("geneName", "orderedLocusNames", "orfNames", "synonyms"):
            value = gene.get(key)
            values = value if isinstance(value, list) else [value]
            for item in values:
                text = str(item.get("value", "")) if isinstance(item, dict) else str(item or "")
                symbol = _searchable_symbol(text)
                if symbol:
                    aliases.add(symbol.casefold())
    return aliases


def _declared_source_accessions(row: dict[str, object]) -> list[str]:
    """Return unique accession-like source ids that this verifier can resolve.

    The source verifier owns literature verification, while this stage owns
    database identity. Keeping that split lets an LLM cite a paper and a stable
    accession in the same ``source_ids`` field without pretending that the
    accession itself proves the biological claim.
    """
    accessions: set[str] = set()
    for token in re.split(r"[;,|]+", str(row.get("source_ids", ""))):
        citation_type, _ = classify_citation(token.strip())
        if citation_type != "database_accession":
            continue
        accession = canonical_accession(token)
        if is_uniprot_accession(accession) or is_agi_locus(accession):
            accessions.add(accession)
    return sorted(accessions)


def _declared_species_alternatives(declared: str) -> list[str]:
    """Retain each explicitly separated reference species for exact matching."""
    return [
        alternative.strip()
        for alternative in re.split(r"\s*[/;|]\s*", declared.strip())
        if alternative.strip()
    ]


def _inherited_benchmark_mode(hypotheses_arg: str) -> bool:
    """Carry a strict policy from the mechanism-import artifact when present."""
    manifest_path = runtime_path(hypotheses_arg).parent / "mechanism_import_manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if not isinstance(manifest, dict):
            raise ValueError("mechanism import policy manifest must be an object")
        mode = manifest.get("benchmark_mode", False)
        if not isinstance(mode, bool):
            raise ValueError("mechanism import benchmark_mode must be boolean")
        return mode
    except FileNotFoundError:
        return False


def hypothesis_verify(
    hypotheses_arg: str,
    out_dir_arg: str,
    focal_species: str,
    record_cache_arg: str | None = None,
    allow_external_api: bool = False,
    benchmark_mode: bool = False,
    min_panel_size: int = 2,
    source_records_arg: str | None = None,
    bound_claims_arg: str | None = None,
    agent_claim_review_arg: str | None = None,
    source_review_records_arg: str | None = None,
    *,
    _entry_events_by_hypothesis: dict[tuple[str, str], list[UniProtEntryEvent]] | None = None,
) -> dict[str, object]:
    """Verify every proposal against an authoritative database before promotion."""
    benchmark_mode = bool(benchmark_mode) or _inherited_benchmark_mode(hypotheses_arg)
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    hypotheses = pd.read_csv(runtime_path(hypotheses_arg), sep="\t", dtype=str).fillna("")
    missing = sorted({"hypothesis_id", "proposed_stable_id", "reference_species"} - set(hypotheses.columns))
    if missing:
        raise ValueError(f"hypothesis table is missing columns: {', '.join(missing)}")
    cache = _load_record_cache(record_cache_arg)
    source_records: dict[str, object] = {}
    if source_records_arg:
        loaded_source_records = json.loads(runtime_path(source_records_arg).read_text(encoding="utf-8"))
        if not isinstance(loaded_source_records, dict):
            raise ValueError("source records must be a JSON object keyed by citation id")
        source_records = loaded_source_records
    bound_claims = load_bound_claims(bound_claims_arg, "a")
    agent_review = None
    if agent_claim_review_arg:
        if bound_claims_arg or not source_review_records_arg:
            raise ValueError("agent claim review requires source review records and cannot mix with --bound-claims")
        from .agent_claim_review import AgentClaimReview
        agent_review = AgentClaimReview(agent_claim_review_arg, hypotheses_arg, source_review_records_arg)
    cached_tair = _cached_tair_index(cache)
    focal_genus = focal_species.split(" ")[0].casefold() if focal_species else ""
    source_annotation_columns = {
        "verified_source_ids", "claim_support_summary", "mechanistic_support_citation_count",
        "publication_source_provenance", "publication_records_sha256",
    }
    source_annotations_present = bool(source_annotation_columns & set(hypotheses.columns))

    family_members: dict[str, list[str]] = {}
    for row in hypotheses.to_dict("records"):
        family = str(row.get("gene_family", "")).strip().casefold()
        if family:
            family_members.setdefault(family, []).append(str(row["hypothesis_id"]))

    rows: list[dict[str, object]] = []
    review_rows: list[dict[str, str]] = []
    # ``verified_sequences`` remains the strict canonical-support payload.
    # A second payload may retain a sequence whose provider identity, species
    # and fingerprint are verified but whose literature/claim binding still
    # needs review.  It is only usable by the explicitly review-only Route-A
    # projection path; it never changes verification_status to ``verified``.
    sequences: dict[str, str] = {}
    reviewable_sequences: dict[str, str] = {}
    for row in hypotheses.to_dict("records"):
        hypothesis_id = str(row["hypothesis_id"])
        accession = str(row["proposed_stable_id"]).strip()
        declared = str(row["reference_species"]).strip()
        record: dict[str, object] | None = None
        provider, source_url, retrieved_at = "", "", ""
        provider_failure = ""
        # A proposal that obeyed "do not invent an accession" arrives with a real
        # gene name and an empty id. Turn the name into an id before verifying,
        # then verify exactly as if the proposal had supplied it.
        resolution_query, resolution_hits = "", 0
        resolution_status = ""
        supplied_accession = canonical_accession(accession)
        explicit_agi_loci = _explicit_agi_loci(str(row.get("reference_gene_name", "")))
        declared_source_accessions = _declared_source_accessions(row) if not accession else []
        if len(explicit_agi_loci) > 1:
            resolution_status = "ambiguous_reference_gene_agi"
            resolution_query = "reference_gene_name:explicit_agi=" + ";".join(explicit_agi_loci)
            resolution_hits = len(explicit_agi_loci)
        elif explicit_agi_loci and accession and supplied_accession.casefold() == explicit_agi_loci[0].casefold():
            accession = supplied_accession = explicit_agi_loci[0]
            resolution_status = "resolved_from_reference_gene_agi"
            resolution_query = "reference_gene_name:explicit_agi=" + supplied_accession
            resolution_hits = 1
        elif explicit_agi_loci and accession and supplied_accession.casefold() != explicit_agi_loci[0].casefold():
            resolution_query = "reference_gene_name:explicit_agi=" + explicit_agi_loci[0]
            resolution_hits = 1
            if is_agi_locus(supplied_accession):
                resolution_status = "conflicting_reference_gene_agi"
            else:
                resolution_status, resolution_query, resolution_hits, provider_failure = (
                    _agi_corresponds_to_accession(
                        explicit_agi_loci[0], supplied_accession, declared,
                        allow_external_api=allow_external_api, cached_tair=cached_tair))
        elif explicit_agi_loci and not accession:
            if len(declared_source_accessions) > 1:
                resolution_status = "ambiguous_declared_source_accessions"
                resolution_query = (
                    "reference_gene_name:explicit_agi=" + explicit_agi_loci[0]
                    + ";source_ids:database_accession"
                )
                resolution_hits = len(declared_source_accessions)
            elif len(declared_source_accessions) == 1:
                declared_source = declared_source_accessions[0]
                if declared_source.casefold() == explicit_agi_loci[0].casefold():
                    supplied_accession = explicit_agi_loci[0]
                    accession = supplied_accession
                    resolution_status = "resolved_from_reference_gene_agi"
                    resolution_query = "reference_gene_name:explicit_agi=" + supplied_accession
                    resolution_hits = 1
                else:
                    correspondence_status, correspondence_query, correspondence_hits, correspondence_error = (
                        _agi_corresponds_to_accession(
                            explicit_agi_loci[0], declared_source, declared,
                            allow_external_api=allow_external_api, cached_tair=cached_tair))
                    if correspondence_status == "resolved_from_agi_correspondence":
                        accession = supplied_accession = declared_source
                        resolution_status = correspondence_status
                        resolution_query = correspondence_query
                        resolution_hits = correspondence_hits
                    elif correspondence_status == "provider_unavailable":
                        resolution_status = correspondence_status
                        provider_failure = correspondence_error
                        resolution_query = correspondence_query
                        resolution_hits = correspondence_hits
                    else:
                        resolution_status = correspondence_status
                        resolution_query = correspondence_query
                        resolution_hits = correspondence_hits
            else:
                supplied_accession = explicit_agi_loci[0]
                accession = supplied_accession
                resolution_status = "resolved_from_reference_gene_agi"
                resolution_query = "reference_gene_name:explicit_agi=" + supplied_accession
                resolution_hits = 1
        blocked_name_resolution = {
            "ambiguous_reference_gene_agi", "conflicting_reference_gene_agi",
            "unresolved_reference_gene_correspondence", "provider_unavailable",
        }
        if resolution_status in blocked_name_resolution:
            accession = ""
            supplied_accession = ""
        elif not accession and len(declared_source_accessions) == 1:
            accession = declared_source_accessions[0]
            supplied_accession = accession
            resolution_status = "resolved_from_declared_source_accession"
            resolution_query = "source_ids:database_accession"
            resolution_hits = 1
        elif not accession and len(declared_source_accessions) > 1:
            resolution_status = "ambiguous_declared_source_accessions"
            resolution_query = "source_ids:database_accession"
            resolution_hits = len(declared_source_accessions)
        if accession and is_agi_locus(supplied_accession) and allow_external_api:
            try:
                accession, resolution_status, resolution_query, resolution_hits = _resolve_tair_locus(
                    supplied_accession, declared)
            except (RuntimeError, TimeoutError, OSError, ValueError) as error:
                accession, resolution_status, provider_failure = "", "provider_unavailable", str(error)
        elif accession and is_agi_locus(supplied_accession):
            cached_accessions = cached_tair.get(supplied_accession, [])
            resolution_query = f"record_cache:xref:tair-{supplied_accession}"
            resolution_hits = len(cached_accessions)
            if len(cached_accessions) == 1:
                accession = cached_accessions[0]
                resolution_status = "resolved_from_tair_locus_cache"
            elif len(cached_accessions) > 1:
                accession = ""
                resolution_status = "ambiguous_tair_locus"
            else:
                accession = supplied_accession
                resolution_status = "not_attempted_offline"
        elif accession and resolution_status in {"resolved_from_declared_source_accession", "resolved_from_agi_correspondence"}:
            pass
        elif resolution_status in blocked_name_resolution:
            pass
        elif accession:
            resolution_status = "supplied_by_proposal"
        elif resolution_status == "ambiguous_declared_source_accessions":
            pass
        elif not allow_external_api:
            resolution_status = "not_attempted_offline"
        else:
            try:
                accession, resolution_status, resolution_query, resolution_hits = _resolve_accession(
                    str(row.get("reference_gene_name", "")), declared)
            except (RuntimeError, TimeoutError, OSError, ValueError) as error:
                accession, resolution_status, provider_failure = "", "provider_unavailable", str(error)
        if accession and resolution_status not in blocked_name_resolution:
            if accession in cache:
                record, provider = cache[accession], "local_record_cache"
                source_url = str(runtime_path(record_cache_arg)) if record_cache_arg else ""
                retrieved_at = str(record.get("candidate_qtl_retrieved_at", "")) or _now()
            elif allow_external_api:
                # Keep the concrete provider outcome attached to this fetch.
                # Several hypotheses may legitimately use the same accession;
                # a later failed fetch must not rewrite an earlier healthy row.
                with capture_uniprot_entry_events() as entry_events:
                    try:
                        record, provider = _fetch_uniprot(accession), "uniprot_rest"
                        source_url = UNIPROT_ENTRY_ENDPOINT + accession + ".json"
                        retrieved_at = _now()
                    except (RuntimeError, TimeoutError, OSError, ValueError) as error:
                        provider_failure = str(error)
                if _entry_events_by_hypothesis is not None:
                    _entry_events_by_hypothesis[(str(row.get("case_id", "")), hypothesis_id)] = entry_events
        retrieved_species = _record_species(record) if record else ""
        sequence = _record_sequence(record) if record else ""
        record_version = _record_version(record) if record else ""
        record_source_ids = _record_source_ids(record) if record else []
        proposal_verified_source_ids = sorted({
            value.strip() for value in re.split(r"[;,|]+", str(row.get("verified_source_ids", "")))
            if value.strip() and classify_citation(value.strip())[0] in {"pmid", "pmcid", "doi"}
        })
        try:
            mechanistic_citation_count = int(str(row.get("mechanistic_support_citation_count", "0")) or "0")
        except ValueError:
            mechanistic_citation_count = 0
        publication_source_provenance = str(row.get("publication_source_provenance", "")).strip()
        publication_records_sha256 = str(row.get("publication_records_sha256", "")).strip()
        publication_record_ids_match = bool(proposal_verified_source_ids) and all(
            source_id in source_records for source_id in proposal_verified_source_ids
        )
        publication_records_actual_sha256 = ""
        if publication_record_ids_match:
            publication_records_actual_sha256 = hashlib.sha256(json.dumps(
                [source_records[source_id] for source_id in proposal_verified_source_ids],
                ensure_ascii=False, sort_keys=True, separators=(",", ":"),
            ).encode("utf-8")).hexdigest()
        provenance_ids_match = False
        provenance_matches_source_records = False
        source_record_mechanistic_count = 0
        try:
            provenance_rows = json.loads(publication_source_provenance) if publication_source_provenance else []
            provenance_rows_are_exact = (
                isinstance(provenance_rows, list)
                and all(isinstance(item, dict) for item in provenance_rows)
                and len(provenance_rows) == len(proposal_verified_source_ids)
            )
            provenance_ids_match = (
                provenance_rows_are_exact
                and {str(item.get("citation", "")) for item in provenance_rows if isinstance(item, dict)}
                == set(proposal_verified_source_ids)
            )
            proposed_provenance = {
                str(item.get("citation", "")): item
                for item in provenance_rows if isinstance(item, dict)
            } if isinstance(provenance_rows, list) else {}
            frozen_provenance: dict[str, object] = {}
            for source_id in proposal_verified_source_ids:
                source_record = source_records.get(source_id, {})
                verifications = (source_record.get("candidate_qtl_hypothesis_verifications", {})
                                 if isinstance(source_record, dict) else {})
                verification = verifications.get(hypothesis_id) if isinstance(verifications, dict) else None
                if isinstance(verification, dict):
                    frozen_provenance[source_id] = verification
                    if str(verification.get("claim_support_status", "")) in {
                        "supported_direct_genetic_test", "supported_perturbation_phenotype",
                    }:
                        source_record_mechanistic_count += 1
            provenance_matches_source_records = (
                provenance_ids_match
                and set(frozen_provenance) == set(proposal_verified_source_ids)
                and proposed_provenance == frozen_provenance
            )
        except json.JSONDecodeError:
            provenance_ids_match = False
            provenance_matches_source_records = False
            frozen_provenance = {}
        # Curated-function alternative path (sibling of issue #611): an
        # authoritative database record whose curated function is itself
        # source-backed may substitute for mechanistic literature support,
        # provided the payload still supplies at least one provider-resolved
        # citation and every anti-tamper hash/provenance check holds.
        proposal_literature_verified = bool(
            proposal_verified_source_ids
            and str(row.get("claim_support_summary", "")) == "mechanistic_support"
            and mechanistic_citation_count > 0
            and provenance_ids_match
            and provenance_matches_source_records
            and mechanistic_citation_count == source_record_mechanistic_count
            and publication_record_ids_match
            and publication_records_actual_sha256 == publication_records_sha256
        )
        # Curated-function alternative path (sibling of issue #611): an
        # authoritative database record whose curated function is itself
        # source-backed may substitute for mechanistic literature support,
        # provided the payload still supplies at least one provider-resolved
        # citation and every anti-tamper hash/provenance check holds.
        curated_function_with_resolved_citation = bool(
            record_source_ids
            and proposal_verified_source_ids
            and provenance_ids_match
            and provenance_matches_source_records
            and any(
                str((frozen_provenance.get(source_id, {}) or {}).get("resolved_id", "")).strip()
                for source_id in proposal_verified_source_ids
            )
            and publication_record_ids_match
            and publication_records_actual_sha256 == publication_records_sha256
        )
        proposal_literature_verified = bool(
            proposal_literature_verified
            or curated_function_with_resolved_citation
        )
        source_ids = sorted(set(record_source_ids) | set(proposal_verified_source_ids))
        bound_claim = bound_claims.get(hypothesis_id)
        if agent_review is not None:
            bound_claim = agent_review.bind(row, accession, record, record_source_ids)
        binding_source_ids = set(source_ids)
        if agent_review is not None and bound_claim.get("_agent_review_receipt") is not None:
            binding_source_ids.add(str(bound_claim["citation"]))
        bound_ok, bound_state, bound_reason = route_a_authorizes(
            bound_claim,
            hypothesis_id=hypothesis_id,
            proposed_accession=accession,
            verified_source_ids=binding_source_ids,
        )
        if agent_review is not None and bound_ok:
            # A separately reviewed relation can resolve lexical triage failure.
            # Do not rewrite the original source classifier or its provenance.
            proposal_literature_verified = True
            source_ids = sorted(binding_source_ids)
        function_source_ids = source_ids
        from .bibliography import INITIAL_POLICY
        if row.get("candidate_admission_policy") == INITIAL_POLICY:
            # Provider bibliographies are retrieval context, not new evidence
            # for the original hypothesis. Keep explicit reviewed bindings.
            # _record_source_ids emits only PubMed citations and this exact
            # reviewed UniProt record handle. Unknown citation shapes stay audit-only.
            record_handles = {f"UniProt:{accession}"} & set(record_source_ids)
            function_source_ids = sorted(set(proposal_verified_source_ids) | record_handles | (
                {str(bound_claim["citation"])} if bound_ok and bound_claim else set()))
        retrieved_genus = retrieved_species.split(" ")[0].casefold() if retrieved_species else ""
        species_match = bool(retrieved_species) and any(
            record_species_matches(record or {}, alternative)
            for alternative in _declared_species_alternatives(declared)
        )
        focal_excluded = not (retrieved_genus and focal_genus and retrieved_genus == focal_genus)
        blob = " ".join([accession, declared, retrieved_species, str(row.get("reference_gene_name", ""))]).casefold()
        answer_derived = benchmark_mode and any(marker in blob for marker in PISUM_MARKERS)

        status, rejection, review = "verified", "", ""
        if provider_failure:
            status, review = "review_required", "accession_provider_unavailable"
        elif resolution_status in blocked_name_resolution:
            status, review = "review_required", (
                "multiple_explicit_agi_in_reference_gene_name"
                if resolution_status == "ambiguous_reference_gene_agi"
                else ("explicit_agi_conflicts_with_proposed_stable_id"
                      if resolution_status == "conflicting_reference_gene_agi"
                      else "explicit_agi_requires_identifier_correspondence_review")
            )
        elif not accession and resolution_status in {
            "ambiguous_gene_name", "ambiguous_tair_locus", "ambiguous_declared_source_accessions",
        }:
            # Several reviewed entries carry this symbol for this species. Choosing
            # one would be a guess dressed as verification.
            status, review = "review_required", (
                "tair_locus_resolves_to_multiple_accessions"
                if resolution_status == "ambiguous_tair_locus"
                else ("declared_source_ids_contain_multiple_accessions"
                      if resolution_status == "ambiguous_declared_source_accessions"
                      else "gene_name_resolves_to_multiple_accessions"))
        elif not accession:
            status, rejection = "rejected", "no_stable_accession_proposed"
        elif record is None and not allow_external_api and accession not in cache:
            status, review = "review_required", "verification_not_performed_offline"
        elif record is None:
            status, rejection = "rejected", "accession_not_found_in_authoritative_database"
        elif not sequence:
            status, rejection = "rejected", "sequence_not_retrievable"
        elif not retrieved_species:
            status, review = "review_required", "retrieved_record_has_no_species"
        elif not species_match:
            status, rejection = "rejected", "declared_species_does_not_match_accession"
        elif not focal_excluded and benchmark_mode:
            status, rejection = "rejected", "reference_gene_is_from_the_focal_species"
        elif answer_derived:
            status, rejection = "rejected", "answer_derived_or_focal_lineage_marker_in_benchmark_mode"
        elif source_annotations_present and not proposal_literature_verified:
            status, review = "review_required", "proposal_literature_is_not_verified"
        elif not source_ids:
            status, review = "review_required", "reference_gene_function_is_not_source_backed"
        elif agent_review is not None and not bound_ok:
            status, review = "review_required", "independent_agent_claim_binding_not_verified"

        named_reference = bool(str(row.get("reference_gene_name", "")).strip()
                               and str(row.get("reference_species", "")).strip())
        resolved_identity = str(row.get("reference_identity_status", "")).strip()
        direct_name = _searchable_symbol(str(row.get("reference_gene_name", ""))).casefold()
        # Resolver-backed named rows carry an explicit exact name/species
        # receipt.  Direct calls have no receipt, so require that the provider
        # record itself declares the requested symbol or alias.  An accession
        # alone never authorizes a named review projection.
        if named_reference:
            name_identity_verified = (
                resolved_identity in PROJECTABLE_IDENTITIES
                if resolved_identity else bool(direct_name and direct_name in _record_gene_aliases(record or {}))
            )
        else:
            # An unnamed provider-resolved accession has no symbol claim to
            # cross-check.  Keep the normal accession/species/sequence gates,
            # while rejecting a resolver state that explicitly says the
            # identity is unresolved.
            name_identity_verified = resolved_identity in {"", "accession_only"}
        identity_sequence_verified = bool(
            accession and record is not None and sequence and retrieved_species
            and species_match and (focal_excluded or not benchmark_mode) and not provider_failure
            and resolution_status not in blocked_name_resolution and not answer_derived
            and name_identity_verified
        )
        if status == "verified" and resolved_identity == IDENTITY_DECLARED_ACCESSION:
            # The declared accession and species agree but the gene name is not
            # on the record and no identity review confirmed it: review-only.
            status, review = "review_required", "reference_name_unconfirmed_on_declared_accession"
        if status == "verified":
            sequences[hypothesis_id] = sequence
        if identity_sequence_verified:
            reviewable_sequences[hypothesis_id] = sequence
        if status == "rejected":
            review_rows.append({"reason": f"proposal_rejected:{rejection}", "target": hypothesis_id,
                                "required_next_action": "the proposal is kept in the rejected ledger and never enters the bait panel",
                                "allowed_choices": "propose another reference gene | keep rejected"})
        if status == "review_required":
            review_rows.append({"reason": f"proposal_review_required:{review}", "target": hypothesis_id,
                                "required_next_action": "resolve the verification gap before this proposal can become a bait",
                                "allowed_choices": "supply a record cache | allow the database API | replace the proposal"})
        family_key = str(row.get("gene_family", "")).strip().casefold()
        rows.append({
            "hypothesis_id": hypothesis_id, "case_id": str(row.get("case_id", "")),
            "proposed_stable_id": accession, "declared_species": declared, "verification_status": status,
            "accession_resolution_status": resolution_status,
            "accession_resolution_query": resolution_query,
            "accession_resolution_hit_count": resolution_hits,
            "accession_exists": bool(record), "retrieved_species": retrieved_species, "species_match": species_match,
            "sequence_retrieved": bool(sequence), "sequence_length": len(sequence),
            "sequence_sha256": hashlib.sha256(sequence.encode("ascii", "ignore")).hexdigest() if sequence else "",
            "identity_sequence_verification_status": "verified" if identity_sequence_verified else "not_verified",
            "function_source_backed": bool(function_source_ids), "function_source_ids": ";".join(function_source_ids),
            "reference_record_source_ids": ";".join(record_source_ids),
            "focal_species_excluded": focal_excluded,
            "proposal_literature_verified": proposal_literature_verified,
            "proposal_verified_source_ids": ";".join(proposal_verified_source_ids),
            "bound_claim_authorized": bound_ok,
            "bound_claim_state": bound_state,
            "bound_claim_reason": bound_reason,
            "bound_claim_curation_method": str((bound_claim or {}).get("curation_method", "")),
            "publication_source_provenance": publication_source_provenance,
            "publication_records_sha256": publication_records_sha256,
            "answer_derived_check": "flagged" if answer_derived else ("checked" if benchmark_mode else "not_applicable"),
            "paralog_alternatives": ";".join(value for value in family_members.get(family_key, []) if value != hypothesis_id),
            "source_provider": provider, "source_url": source_url, "retrieved_at": retrieved_at,
            "source_record_version": record_version,
            "accession_provider_error": provider_failure,
            "rejection_reason": rejection, "review_reason": review,
        })

    validation = _frame(rows, VALIDATION_COLUMNS)
    verified = validation[validation["verification_status"] == "verified"] if not validation.empty else validation
    rejected = validation[validation["verification_status"] == "rejected"] if not validation.empty else validation
    needs_review = validation[validation["verification_status"] == "review_required"] if not validation.empty else validation
    if len(verified) < max(1, min_panel_size):
        review_rows.append({"reason": "panel_below_minimum_size", "target": "verified bait panel",
                            "required_next_action": f"a panel of fewer than {min_panel_size} verified baits is review-required, not automatically accepted",
                            "allowed_choices": "verify more proposals | accept explicitly as a single-bait panel | abstain"})

    reviews = _frame(review_rows, REVIEW_COLUMNS).drop_duplicates().sort_values(["reason", "target"], kind="stable")
    paths = {"proposal_validation": out / "proposal_validation.tsv",
             "proposal_rejected": out / "proposal_rejected.tsv",
             "proposal_review_required": out / "review_required.tsv",
             "verified_sequences": out / "verified_sequences.json",
             "reviewable_sequences": out / "reviewable_sequences.json"}
    _write(validation, paths["proposal_validation"])
    _write(rejected, paths["proposal_rejected"])
    _write(reviews, paths["proposal_review_required"])
    paths["verified_sequences"].write_text(
        json.dumps({key: sequences[key] for key in sorted(sequences)}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    paths["reviewable_sequences"].write_text(
        json.dumps({key: reviewable_sequences[key] for key in sorted(reviewable_sequences)}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if agent_review is not None:
        paths["agent_claim_review_audit"] = agent_review.write_audit(out)

    manifest = {
        "schema_version": SCHEMA_VERSION, "focal_species": focal_species, "benchmark_mode": benchmark_mode,
        "verification_provider": "uniprot_rest" if allow_external_api else ("local_record_cache" if cache else "none"),
        "external_api_allowed": allow_external_api,
        "record_cache": str(runtime_path(record_cache_arg)) if record_cache_arg else None,
        "source_records": (file_fingerprint(runtime_path(source_records_arg), True)
                           if source_records_arg else None),
        "bound_claims": (file_fingerprint(runtime_path(bound_claims_arg), True)
                         if bound_claims_arg else None),
        "agent_claim_review": (file_fingerprint(runtime_path(agent_claim_review_arg), True)
                               if agent_claim_review_arg else None),
        "source_review_records": (file_fingerprint(runtime_path(source_review_records_arg), True)
                                  if source_review_records_arg else None),
        "agent_claim_review_audit": (file_fingerprint(paths["agent_claim_review_audit"], True)
                                     if agent_review is not None else None),
        "agent_review_validation_sha256": (hashlib.sha256(paths["proposal_validation"].read_bytes()).hexdigest()
                                           if agent_review is not None else None),
        "agent_review_hypotheses_sha256": (hashlib.sha256(runtime_path(hypotheses_arg).read_bytes()).hexdigest()
                                           if agent_review is not None else None),
        "hypotheses_input": file_fingerprint(runtime_path(hypotheses_arg), False),
        "publication_gate_scope": (
            "enforced_from_source_verify_annotations"
            if source_annotations_present else "legacy_standalone_input_without_source_annotations"
        ),
        "proposal_count": int(len(validation)), "verified_count": int(len(verified)),
        "rejected_count": int(len(rejected)), "review_required_count": int(len(needs_review)),
        "rejection_reasons": {str(k): int(v) for k, v in rejected["rejection_reason"].value_counts().sort_index().items()} if not rejected.empty else {},
        "gates": [
            "the stable accession must exist in an authoritative database",
            "the retrieved species must match the declared species",
            "the sequence must be re-retrievable and is fingerprinted",
            ("a reference gene from the focal species is rejected in benchmark mode"
             if benchmark_mode else "same-species references remain eligible subject to identity, sequence, species and source gates"),
            "in benchmark mode, focal-lineage / answer-backmapped markers are rejected",
            "a reference gene without a source-backed function stays review_required",
            "when source verification annotations are present, either mechanistic support from a verified proposal paper or an authoritative curated function plus one provider-resolved proposal paper is required independently of accession identity",
            "canonical support also requires a shared bound claim joining citation, gene, experiment, phenotype, organism and scope",
            "independent_agent_review requires distinct recorded executor contexts plus replayed text/ID/species checks; semantic attestation is not deterministic proof or human curation",
            "a single declared source accession may be resolved; multiple accessions are never tie-broken",
            "rejected proposals are retained with their reason and never enter the panel",
        ],
        "outputs": {name: str(path) for name, path in sorted(paths.items())},
    }
    manifest_path = out / "hypothesis_verification_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), "validation": str(paths["proposal_validation"]),
            "rejected": str(paths["proposal_rejected"]), "review_required": str(paths["proposal_review_required"]),
            "verified_sequences": str(paths["verified_sequences"]),
            "reviewable_sequences": str(paths["reviewable_sequences"]),
            **({"agent_claim_review_audit": str(paths["agent_claim_review_audit"])} if agent_review is not None else {}),
            "verified_count": int(len(verified)), "rejected_count": int(len(rejected)),
            "review_required_count": int(len(needs_review))}


# --------------------------------------------------------------------------- #
# 4. bait panel freeze
# --------------------------------------------------------------------------- #

def bait_panel_freeze(
    hypotheses_arg: str,
    validation_arg: str,
    sequences_arg: str,
    out_dir_arg: str,
    case_id: str = "",
    min_panel_size: int = 2,
    reviewable_sequences_arg: str | None = None,
) -> dict[str, object]:
    """Freeze canonical baits plus identity-verified review-only baits.

    ``verified_bait_panel.tsv`` remains the canonical-support input.  The
    companion ``projection_bait_panel.tsv`` additionally carries baits whose
    provider identity/species/sequence checks passed but whose source claim is
    incomplete.  Those rows are labelled review-only and are never support.
    """
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    hypotheses_frame = pd.read_csv(runtime_path(hypotheses_arg), sep="\t", dtype=str).fillna("")
    if "hypothesis_id" not in hypotheses_frame.columns:
        raise ValueError("hypothesis table requires hypothesis_id")
    if hypotheses_frame["hypothesis_id"].duplicated().any():
        raise ValueError("hypothesis table requires unique hypothesis_id")
    hypotheses = {str(row["hypothesis_id"]): row for row in hypotheses_frame.to_dict("records")}
    validation = pd.read_csv(runtime_path(validation_arg), sep="\t", dtype=str).fillna("")
    sequences = json.loads(runtime_path(sequences_arg).read_text(encoding="utf-8"))
    reviewable_sequences = (json.loads(runtime_path(reviewable_sequences_arg).read_text(encoding="utf-8"))
                            if reviewable_sequences_arg else sequences)
    agent_audit_fingerprint = None
    if "bound_claim_curation_method" in validation and validation["bound_claim_curation_method"].eq("independent_agent_review").any():
        verification_manifest = runtime_path(validation_arg).parent / "hypothesis_verification_manifest.json"
        if not verification_manifest.is_file():
            raise ValueError("agent-reviewed bait requires the hypothesis verification manifest and audit")
        sealed = json.loads(verification_manifest.read_text(encoding="utf-8"))
        agent_audit_fingerprint = sealed.get("agent_claim_review_audit") or {}
        audit_path = runtime_path(str(agent_audit_fingerprint.get("path", "")))
        if (sealed.get("agent_review_validation_sha256") != hashlib.sha256(runtime_path(validation_arg).read_bytes()).hexdigest()
                or sealed.get("agent_review_hypotheses_sha256") != hashlib.sha256(runtime_path(hypotheses_arg).read_bytes()).hexdigest()
                or not audit_path.is_file()
                or agent_audit_fingerprint.get("sha256") != hashlib.sha256(audit_path.read_bytes()).hexdigest()):
            raise ValueError("agent-reviewed bait inputs/audit differ from the verified generation")

    review_rows: list[dict[str, str]] = []
    rows: list[dict[str, object]] = []
    projection_rows: list[dict[str, object]] = []
    query_rows: list[dict[str, object]] = []
    fasta_lines: list[str] = []
    for record in validation.sort_values("hypothesis_id", kind="stable").to_dict("records"):
        if str(record["verification_status"]) not in {"verified", "review_required"}:
            continue
        hypothesis_id = str(record["hypothesis_id"])
        bound_authorized = str(record.get("bound_claim_authorized", "")).strip().casefold() in {"true", "1"}
        # Claim binding cannot upgrade a provisional reference identity.
        bound_authorized = bound_authorized and str(record["verification_status"]) == "verified"
        if (not bound_authorized
                and str(record.get("identity_sequence_verification_status", "")).strip().casefold() != "verified"):
            review_rows.append({
                "reason": "canonical_support_blocked:source_claim_binding_not_verified",
                "target": hypothesis_id,
                "required_next_action": "supply a complete shared bound claim joining citation, gene, experiment and phenotype",
                "allowed_choices": "curate the bound claim | keep the proposal ledger-only",
            })
            review_rows.append({
                "reason": "review_projection_blocked:identity_species_sequence_not_verified",
                "target": hypothesis_id,
                "required_next_action": "re-run hypothesis-verify with provider-confirmed identity, species and sequence",
                "allowed_choices": "re-verify | keep the proposal ledger-only",
            })
            continue
        sequence_source = sequences if bound_authorized else reviewable_sequences
        sequence = str(sequence_source.get(hypothesis_id, ""))
        if not sequence:
            review_rows.append({"reason": "verified_proposal_without_sequence", "target": hypothesis_id,
                                "required_next_action": "re-run hypothesis-verify; a verified bait must carry its retrieved sequence",
                                "allowed_choices": "re-verify | drop the bait"})
            continue
        digest = hashlib.sha256(sequence.encode("ascii", "ignore")).hexdigest()
        if digest != str(record["sequence_sha256"]):
            review_rows.append({"reason": "sequence_fingerprint_mismatch", "target": hypothesis_id,
                                "required_next_action": "the stored sequence does not match its recorded fingerprint; re-verify before freezing",
                                "allowed_choices": "re-verify | drop the bait"})
            continue
        source = hypotheses.get(hypothesis_id, {})
        from .bibliography import INITIAL_POLICY, concrete_bridge
        if source.get("candidate_admission_policy") == INITIAL_POLICY and not (
            source.get("bibliography_state") == "resolved" and concrete_bridge(source)
        ):
            # Missing metadata/bridge is a ledger state, not biological rejection.
            # Record it: a sequence-verified bait must never leave projection silently.
            state = str(source.get("bibliography_state", "")).strip() or "not_evaluated"
            review_rows.append({
                "reason": ("projection_blocked:bibliography_" + state
                           if state != "resolved" else "projection_blocked:no_concrete_bridge"),
                "target": hypothesis_id,
                "required_next_action": (
                    "supply a stable reference ID for the proposal (or a paper reference_gene_id the provider "
                    "record cross-references) so its papers bind to this endpoint, then re-run"),
                "allowed_choices": "supply the reference ID | keep the proposal ledger-only",
            })
            continue
        bait_id = str(record["proposed_stable_id"])
        case = case_id or str(record.get("case_id", "")) or str(source.get("case_id", ""))
        entry = {
            "case_id": case, "bait_id": bait_id, "hypothesis_id": hypothesis_id,
            "pathway_or_process": str(source.get("pathway_or_process", "")),
            "gene_family": str(source.get("gene_family", "")),
            "reference_gene_name": str(source.get("reference_gene_name", "")),
            "reference_species": str(record.get("retrieved_species", "")),
            "expected_role": str(source.get("expected_role", "")),
            "sequence_length": len(sequence), "sequence_sha256": digest,
            "source_provider": str(record.get("source_provider", "")), "source_url": str(record.get("source_url", "")),
            "retrieved_at": str(record.get("retrieved_at", "")), "function_source_ids": str(record.get("function_source_ids", "")),
            "rank_contribution": "none", "panel_status": "frozen",
            "review_status": "verified" if bound_authorized else "review_required",
            "accession_resolution_status": str(record.get("accession_resolution_status", "")),
            "accession_resolution_query": str(record.get("accession_resolution_query", "")),
            "accession_resolution_hit_count": str(record.get("accession_resolution_hit_count", "")),
            "proposal_verified_source_ids": str(record.get("proposal_verified_source_ids", "")),
            "source_record_version": str(record.get("source_record_version", "")),
            "publication_source_provenance": str(record.get("publication_source_provenance", "")),
            "publication_records_sha256": str(record.get("publication_records_sha256", "")),
            "bound_claim_curation_method": str(record.get("bound_claim_curation_method", "")),
            "projection_eligibility": "verified_support" if bound_authorized else "review_only",
            "projection_reason": "" if bound_authorized else "source_claim_binding_not_verified",
        }
        projection_rows.append(entry)
        if not bound_authorized:
            review_rows.append({
                "reason": "canonical_support_blocked:source_claim_binding_not_verified",
                "target": hypothesis_id,
                "required_next_action": (
                    "supply a complete shared bound claim joining the citation, subject gene, "
                    "experiment and target phenotype"
                ),
                "allowed_choices": "curate the bound claim | keep the proposal review-only",
            })
            continue
        rows.append(entry)
        species = str(record["retrieved_species"])
        query_rows.append({
            "case_id": case, "homolog_query_id": bait_id, "homolog_source_id": f"UniProt:{bait_id}",
            "query_species": species, "query_genus": species.split(" ")[0] if species else "",
            # Controlled vocabulary owned by `schemas`; the freeze reuses it rather
            # than inventing a parallel one, so `query-eligibility` stays the gate.
            "sequence_origin": "non_pisum_reference_proteome",
            "selection_basis": "pathway_member",
        })
        fasta_lines.append(f">{bait_id} {str(record['retrieved_species'])}")
        fasta_lines.extend(sequence[index:index + 60] for index in range(0, len(sequence), 60))

    panel = _frame(rows, PANEL_COLUMNS)
    query_manifest = _frame(query_rows, QUERY_MANIFEST_COLUMNS)
    panel_status = "frozen"
    if len(panel) < max(1, min_panel_size):
        panel_status = "review_required_small_panel"
        review_rows.append({"reason": "panel_below_minimum_size", "target": case_id or "case",
                            "required_next_action": f"the frozen panel has {len(panel)} bait(s); fewer than {min_panel_size} is review-required",
                            "allowed_choices": "verify more proposals | accept explicitly | abstain"})
    # Per-bait verification and panel-level mechanism coverage are distinct.
    # Every row here already passed the accession/source/species/sequence gates;
    # a small panel therefore remains review-required as a *panel* without
    # invalidating those rows for non-scoring deterministic interval projection.
    # Copy ambiguity remains downstream and the panel warning is preserved.
    row_review_status = "verified" if len(panel) else "review_required"
    for row in rows:
        row["review_status"] = row_review_status
    # Rebuild the frame so the derived per-row review_status reaches the TSV.
    panel = _frame(rows, PANEL_COLUMNS)
    reviews = _frame(review_rows, REVIEW_COLUMNS).drop_duplicates().sort_values(["reason", "target"], kind="stable")

    panel_path = out / "verified_bait_panel.tsv"
    projection_panel_path = out / "projection_bait_panel.tsv"
    fasta_path = out / "verified_bait_panel.fasta"
    query_path = out / "verified_bait_panel_query_manifest.tsv"
    _write(panel, panel_path)
    _write(_frame(projection_rows, PANEL_COLUMNS), projection_panel_path)
    _write(query_manifest, query_path)
    _write(reviews, out / "review_required.tsv")
    fasta_path.write_text("\n".join(fasta_lines) + ("\n" if fasta_lines else ""), encoding="utf-8")

    bait_id_counts: dict[str, int] = {}
    for row in rows:
        bait_id_counts[str(row["bait_id"])] = bait_id_counts.get(str(row["bait_id"]), 0) + 1

    def manifest_bait_key(row: dict[str, object]) -> str:
        bait_id = str(row["bait_id"])
        return bait_id if bait_id_counts[bait_id] == 1 else f"{bait_id}::{row['hypothesis_id']}"

    manifest = {
        "schema_version": SCHEMA_VERSION, "case_id": case_id, "panel_status": panel_status,
        "bait_count": int(len(panel)),
        "bait_review_status": row_review_status,
        "pathways": sorted({str(row["pathway_or_process"]) for row in rows if row["pathway_or_process"]}),
        "gene_families": sorted({str(row["gene_family"]) for row in rows if row["gene_family"]}),
        "panel_sha256": hashlib.sha256(panel_path.read_bytes()).hexdigest(),
        "projection_panel": projection_panel_path.name,
        "projection_panel_sha256": hashlib.sha256(projection_panel_path.read_bytes()).hexdigest(),
        "fasta_sha256": hashlib.sha256(fasta_path.read_bytes()).hexdigest(),
        "sequence_fingerprints": {manifest_bait_key(row): str(row["sequence_sha256"]) for row in rows},
        "provenance": {manifest_bait_key(row): {"provider": row["source_provider"], "url": row["source_url"],
                                             "retrieved_at": row["retrieved_at"],
                                             "record_version": row["source_record_version"],
                                             "function_source_ids": row["function_source_ids"],
                                             "proposal_verified_source_ids": row["proposal_verified_source_ids"],
                                             "publication_source_provenance": row["publication_source_provenance"],
                                             "publication_records_sha256": row["publication_records_sha256"],
                                             "accession_resolution_status": row["accession_resolution_status"],
                                             "accession_resolution_query": row["accession_resolution_query"],
                                             "accession_resolution_hit_count": row["accession_resolution_hit_count"]}
                       for row in rows},
        "inputs": {"hypotheses": file_fingerprint(runtime_path(hypotheses_arg)),
                   "validation": file_fingerprint(runtime_path(validation_arg)),
                   "reviewable_sequences": (file_fingerprint(runtime_path(reviewable_sequences_arg))
                                            if reviewable_sequences_arg else None)},
        "agent_claim_review_audit": agent_audit_fingerprint,
        "downstream_contract": "verified_bait_panel_query_manifest.tsv is consumable by `query-eligibility`; interval mapping and ranking are unchanged",
        "safety": [
            "only bound-claim-authorized rows enter verified_bait_panel.tsv and canonical support",
            "identity/species/sequence-verified rows with incomplete claim binding may enter projection_bait_panel.tsv as review_only only",
            "each bait carries its retrieval provenance and sequence fingerprint",
            "panel diversity review does not overwrite per-bait machine verification",
            "the panel contributes no score: rank_contribution is none for every bait",
        ],
    }
    manifest_path = out / "bait_panel_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), "panel": str(panel_path), "projection_panel": str(projection_panel_path), "fasta": str(fasta_path),
            "query_manifest": str(query_path), "review_required": str(out / "review_required.tsv"),
            "bait_count": int(len(panel)), "panel_status": panel_status,
            "review_status": row_review_status}


# --------------------------------------------------------------------------- #
# 5. staged benchmark
# --------------------------------------------------------------------------- #

BENCHMARK_CASE_COLUMNS = [
    "case_id", "expected_pathway", "expected_gene_family", "causal_gene_id",
    "hypotheses", "validation", "panel", "ranked",
]
BENCHMARK_RESULT_COLUMNS = [
    "case_id", "expected_pathway", "pathway_rank", "pathway_top1", "pathway_top3", "pathway_top5",
    "gene_family_recalled", "valid_nonfocal_bait", "proposal_count", "accession_proposed",
    "accession_exists", "accession_species_match", "verified_count", "panel_size",
    "causal_gene_rank", "causal_in_top20", "single_pathway_overconfident", "abstained", "failure_reason",
]


def _rank_of(values: list[str], expected: str) -> int:
    wanted = expected.strip().casefold()
    for index, value in enumerate(values, 1):
        if value.strip().casefold() == wanted:
            return index
    return 0


def phenotype_pathway_benchmark(config_arg: str, out_dir_arg: str, top_n: int = 20) -> dict[str, object]:
    """Score the proposal layer stage by stage, separately from the map-based lane."""
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    config = pd.read_csv(runtime_path(config_arg), sep="\t", dtype=str).fillna("")
    missing = sorted({"case_id", "expected_pathway", "expected_gene_family"} - set(config.columns))
    if missing:
        raise ValueError(f"benchmark config is missing columns: {', '.join(missing)}")

    rows: list[dict[str, object]] = []
    for case in config.to_dict("records"):
        case_id = str(case["case_id"])
        hypotheses = pd.read_csv(runtime_path(case["hypotheses"]), sep="\t", dtype=str).fillna("") if case.get("hypotheses") else pd.DataFrame()
        validation = pd.read_csv(runtime_path(case["validation"]), sep="\t", dtype=str).fillna("") if case.get("validation") else pd.DataFrame()
        panel = pd.read_csv(runtime_path(case["panel"]), sep="\t", dtype=str).fillna("") if case.get("panel") else pd.DataFrame()
        ranked = pd.read_csv(runtime_path(case["ranked"]), sep="\t", dtype=str).fillna("") if case.get("ranked") else pd.DataFrame()

        pathways = list(hypotheses["pathway_or_process"]) if "pathway_or_process" in hypotheses.columns else []
        families = {value.strip().casefold() for value in hypotheses.get("gene_family", pd.Series(dtype=str))}
        pathway_rank = _rank_of(pathways, str(case["expected_pathway"]))
        family_recalled = str(case["expected_gene_family"]).strip().casefold() in families
        accession_proposed = int(hypotheses["proposed_stable_id"].astype(str).str.strip().ne("").sum()) if "proposed_stable_id" in hypotheses.columns else 0
        accession_exists = int(validation["accession_exists"].astype(str).eq("True").sum()) if "accession_exists" in validation.columns else 0
        species_match = int(validation["species_match"].astype(str).eq("True").sum()) if "species_match" in validation.columns else 0
        verified = int(validation["verification_status"].eq("verified").sum()) if "verification_status" in validation.columns else 0
        causal_rank = 0
        if not ranked.empty and "gene_id" in ranked.columns and case.get("causal_gene_id"):
            gene_ids = list(ranked["gene_id"].astype(str))
            causal_rank = _rank_of(gene_ids, str(case["causal_gene_id"]))
        distinct = len({value.strip().casefold() for value in pathways if value.strip()})
        abstained = not len(hypotheses)
        failure = ""
        if abstained:
            failure = "abstention_or_unparsable_output"
        elif not pathway_rank:
            failure = "pathway_miss"
        elif not family_recalled:
            failure = "family_miss"
        elif not accession_proposed:
            failure = "retrieval_miss"
        elif accession_exists == 0:
            failure = "accession_hallucination"
        elif verified == 0:
            failure = "verification_miss"
        elif len(panel) == 0:
            failure = "panel_miss"
        elif case.get("causal_gene_id") and not (0 < causal_rank <= top_n):
            failure = "mapping_or_ranking_miss"
        rows.append({
            "case_id": case_id, "expected_pathway": str(case["expected_pathway"]), "pathway_rank": pathway_rank,
            "pathway_top1": bool(pathway_rank == 1), "pathway_top3": bool(0 < pathway_rank <= 3),
            "pathway_top5": bool(0 < pathway_rank <= 5), "gene_family_recalled": family_recalled,
            "valid_nonfocal_bait": bool(len(panel)), "proposal_count": int(len(hypotheses)),
            "accession_proposed": accession_proposed, "accession_exists": accession_exists,
            "accession_species_match": species_match, "verified_count": verified, "panel_size": int(len(panel)),
            "causal_gene_rank": causal_rank, "causal_in_top20": bool(0 < causal_rank <= top_n),
            "single_pathway_overconfident": bool(len(hypotheses) and distinct < 2),
            "abstained": abstained, "failure_reason": failure,
        })

    detail = _frame(rows, BENCHMARK_RESULT_COLUMNS)
    detail_path = out / "phenotype_pathway_cases.tsv"
    _write(detail, detail_path)

    total = max(1, len(rows))
    proposed_total = sum(int(row["accession_proposed"]) for row in rows)
    summary = {
        "case_count": len(rows),
        "pathway_top1_recall": round(sum(1 for row in rows if row["pathway_top1"]) / total, 4),
        "pathway_top3_recall": round(sum(1 for row in rows if row["pathway_top3"]) / total, 4),
        "pathway_top5_recall": round(sum(1 for row in rows if row["pathway_top5"]) / total, 4),
        "causal_gene_family_recall": round(sum(1 for row in rows if row["gene_family_recalled"]) / total, 4),
        "cases_with_valid_nonfocal_bait": round(sum(1 for row in rows if row["valid_nonfocal_bait"]) / total, 4),
        "stable_accession_existence_rate": round(sum(int(row["accession_exists"]) for row in rows) / max(1, proposed_total), 4),
        "accession_species_match_rate": round(sum(int(row["accession_species_match"]) for row in rows) / max(1, proposed_total), 4),
        "verification_pass_rate": round(sum(int(row["verified_count"]) for row in rows) / max(1, proposed_total), 4),
        "causal_recovered_in_top20_from_verified_baits": round(sum(1 for row in rows if row["causal_in_top20"]) / total, 4),
        "wrong_single_pathway_overconfidence_rate": round(
            sum(1 for row in rows if row["single_pathway_overconfident"] and not row["pathway_top1"]) / total, 4),
        "abstention_rate": round(sum(1 for row in rows if row["abstained"]) / total, 4),
        "failure_reasons": {reason: sum(1 for row in rows if row["failure_reason"] == reason)
                            for reason in sorted({str(row["failure_reason"]) for row in rows if row["failure_reason"]})},
    }
    summary_path = out / "phenotype_pathway_summary.tsv"
    _write(pd.DataFrame([{"metric": key, "value": json.dumps(value, ensure_ascii=False) if isinstance(value, dict) else value}
                         for key, value in summary.items()]), summary_path)

    manifest = {
        "schema_version": SCHEMA_VERSION, "lane": "phenotype_only_proposal_lane", "top_n": top_n,
        "separation_rule": "this lane is scored separately from the map-based homolog-informed lane and must never be pooled with it",
        "metrics": summary, "config": file_fingerprint(runtime_path(config_arg)),
        "outputs": {"cases": str(detail_path), "summary": str(summary_path)},
        "interpretation": [
            "pathway and family recall measure the proposal layer, not the deterministic rank",
            "a verified bait recovering the causal homolog is still map-based prioritization, not autonomous discovery",
            "abstention is reported, not penalised into a forced answer",
        ],
    }
    manifest_path = out / "phenotype_pathway_benchmark_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), "cases": str(detail_path), "summary": str(summary_path),
            "case_count": len(rows), "metrics": summary}
