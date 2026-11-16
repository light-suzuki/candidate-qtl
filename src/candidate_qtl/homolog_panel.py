"""Automate phenotype-to-reference-homolog discovery before Stage-A review.

The semantic proposal step may use no model, a local model, or an explicitly
allowed external model.  UniProt records, stable source IDs, protein sequences,
and interval-only BLASTP are then validated deterministically.  The resulting
rows remain pending review and contribute only weak comparative-function
screening evidence; they are never causal calls.
"""

from __future__ import annotations

import hashlib
import json
import re
import time
import urllib.error
import urllib.request
from datetime import UTC, datetime
from pathlib import Path
from urllib.parse import urlencode

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from .core import file_fingerprint, runtime_path, safe_id
from .homolog_shared import _run_interval_blastp, homolog_blast_evidence
from .hardening import prioritize_evidence
from .llm import run_llm_prompts
from .shortlist import screening_shortlist


UNIPROT_ENDPOINT = "https://rest.uniprot.org/uniprotkb/search"
PANEL_VERSION = "1.0.0"
PANEL_COLUMNS = [
    "reference_gene_id", "reference_gene_name", "reference_species",
    "trait_relation", "evidence_level", "source_id", "discovery_method",
    "review_status", "notes", "auto_score", "search_hypotheses",
    "evidence_channels", "uniprot_entry_version", "sequence_length",
]
# Trait-independent grammatical and measurement-shell words cannot carry
# phenotype-specific support on their own.  Keep this shared because both
# mechanical query generation and lexical overlap consume STOPWORDS.
STOPWORDS = {
    "and", "or", "the", "a", "an", "of", "in", "to", "for", "with",
    "trait", "phenotype", "plant", "pea", "pisum", "sativum", "versus",
    "from", "under", "over", "during", "after", "before", "between", "through", "into", "onto",
    "that", "which", "when", "where", "while", "because",
    "is", "are", "was", "were", "be", "been", "being", "has", "have", "had",
    "its", "their", "this", "these", "those",
    "level", "levels", "number", "numbers", "amount", "content", "component", "components", "type", "types",
}
MORPHOLOGY_VARIANTS = {
    "fasciated": "fasciation",
    "wrinkled": "wrinkle",
    "dwarfed": "dwarfism",
    "dwarf": "dwarfism",
    "elongated": "elongation",
    "branched": "branching",
    "pigmented": "pigmentation",
    "resistant": "resistance",
    "susceptible": "susceptibility",
}
PANEL_SYSTEM_PROMPT = """You propose search hypotheses for an automated plant-QTL homolog panel.
Use only the supplied phenotype, inheritance, direction, and species; no QTL genes or sequences are available.
Return JSON only as {"search_hypotheses":[{"query":"short phenotype or mechanism phrase","mechanism":"short mechanism","known_gene":"optional public model-gene symbol","rationale":"one sentence"}]}.
Return 6-10 diverse, database-ready hypotheses. Keep query to 2-5 words. Put exactly one canonical model-plant gene symbol in known_gene, or use an empty string.
Prioritize phenotype-tested model genes and alternative mechanisms, not generic stress explanations.
Do not invent accessions or citations. The next stage will validate every proposal against UniProt and PubMed-linked annotations."""


def _http_json(url: str, timeout: int = 60, retries: int = 3) -> dict[str, object]:
    request = urllib.request.Request(url, headers={"User-Agent": "candidate-qtl/0.6 homolog-panel"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as error:
            if error.code not in {429, 500, 502, 503, 504} or attempt + 1 == retries:
                detail = error.read().decode("utf-8", errors="replace")
                raise RuntimeError(f"UniProt API HTTP {error.code}: {detail[:300]}") from error
        except urllib.error.URLError as error:
            if attempt + 1 == retries:
                raise RuntimeError(f"UniProt API connection failed: {error.reason}") from error
        time.sleep(1.0 * (attempt + 1))
    raise RuntimeError("UniProt API retry loop ended unexpectedly")


def _terms(text: str) -> set[str]:
    return {
        token for token in re.findall(r"[A-Za-z][A-Za-z0-9_-]{2,}", text.casefold())
        if token not in STOPWORDS
    }


def _overlap(term: str, text: str) -> bool:
    wanted = _terms(term)
    observed = _terms(text)
    if not wanted:
        return False
    if wanted.issubset(observed):
        return True
    required = min(2, len(wanted)) if len(wanted) <= 5 else max(3, (2 * len(wanted) + 4) // 5)
    return len(wanted & observed) >= required


def _mechanical_queries(phrase: str) -> list[dict[str, str]]:
    """Expand phenotype prose reproducibly without semantic model judgment."""
    cleaned = re.sub(r"\s+", " ", phrase).strip()
    if not cleaned:
        return []
    queries = [(cleaned, "input_full")]
    segments = [
        value.strip(" ,/()")
        for value in re.split(r"(?i)\b(?:versus|vs\.?|and|or)\b|[|;,/]", cleaned)
        if value.strip(" ,/()")
    ]
    queries.extend((value, "input_segment") for value in segments if value.casefold() != cleaned.casefold())
    tokens = [
        token for token in re.findall(r"[A-Za-z][A-Za-z0-9_-]{2,}", cleaned)
        if token.casefold() not in STOPWORDS
    ]
    queries.extend((token, "input_token") for token in tokens)
    for token in tokens:
        variant = MORPHOLOGY_VARIANTS.get(token.casefold())
        if variant:
            queries.append((variant, "input_variant"))
    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    for query, origin in queries:
        key = query.casefold()
        if key in seen:
            continue
        seen.add(key)
        rows.append({"query": query, "origin": origin, "mechanism": "", "known_gene": "", "rationale": "deterministic phenotype-text expansion"})
    return rows[:12]


def _known_gene_queries(value: str) -> list[str]:
    if value.strip().casefold() in {"", "none", "null", "n/a", "na", "unknown"}:
        return []
    tokens = re.findall(r"[A-Za-z][A-Za-z0-9_.-]{1,30}", value)
    blocked = {"and", "or", "gene", "protein", "none", "unknown"}
    return list(dict.fromkeys(token for token in tokens if token.casefold() not in blocked))[:3]


def _repair_hypothesis_json(raw_text: str) -> dict[str, object] | None:
    """Repair narrow, line-local quote omissions emitted by small LLMs."""
    cleaned = re.sub(r"\x1b\[[0-?]*[ -/]*[@-~]", "", raw_text).strip()
    start, end = cleaned.find("{"), cleaned.rfind("}")
    if start < 0 or end <= start:
        return None
    lines = cleaned[start:end + 1].splitlines()
    repaired: list[str] = []
    index = 0
    changed = False
    value_keys = {"query", "mechanism", "known_gene", "rationale"}
    while index < len(lines):
        line = lines[index]
        scalar_match = re.match(r'^(\s*)"([A-Za-z_]+)"\s*:\s*"(.*?)(,?)\s*$', line)
        if (
            scalar_match
            and scalar_match.group(2) in value_keys
            and len(re.findall(r'(?<!\\)"', line)) == 3
        ):
            line = line[:-1] + '",' if scalar_match.group(4) else line + '"'
            repaired.append(line)
            index += 1
            changed = True
            continue
        match = re.search(r'"([A-Za-z_]+)"\s*:\s*$', line)
        if match and match.group(1) in value_keys and index + 1 < len(lines):
            value_line = lines[index + 1].strip()
            if value_line and not value_line.startswith(('"', "{", "[")):
                trailing_comma = value_line.endswith(",")
                value = value_line[:-1].rstrip() if trailing_comma else value_line
                if value.endswith('"'):
                    value = value[:-1]
                repaired.append(line + " " + json.dumps(value) + ("," if trailing_comma else ""))
                index += 2
                changed = True
                continue
        repaired.append(line)
        index += 1
    if not changed:
        return None
    try:
        payload = json.loads("\n".join(repaired))
    except json.JSONDecodeError:
        return None
    return payload if isinstance(payload, dict) else None


def _comment_rows(record: dict[str, object]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for comment in record.get("comments", []):
        if not isinstance(comment, dict):
            continue
        channel = str(comment.get("commentType", ""))
        for item in comment.get("texts", []):
            if not isinstance(item, dict):
                continue
            pmids = [
                str(ev.get("id")) for ev in item.get("evidences", [])
                if isinstance(ev, dict) and str(ev.get("source", "")).casefold() == "pubmed" and ev.get("id")
            ]
            rows.append({"channel": channel, "text": str(item.get("value", "")), "pmids": pmids})
    for feature in record.get("features", []):
        if not isinstance(feature, dict) or str(feature.get("type", "")) != "Mutagenesis":
            continue
        pmids = [
            str(ev.get("id")) for ev in feature.get("evidences", [])
            if isinstance(ev, dict) and str(ev.get("source", "")).casefold() == "pubmed" and ev.get("id")
        ]
        rows.append({"channel": "MUTAGENESIS", "text": str(feature.get("description", "")), "pmids": pmids})
    return rows


def _protein_name(record: dict[str, object]) -> str:
    description = record.get("proteinDescription", {})
    if not isinstance(description, dict):
        return ""
    recommended = description.get("recommendedName", {})
    if isinstance(recommended, dict):
        full = recommended.get("fullName", {})
        if isinstance(full, dict) and full.get("value"):
            return str(full["value"])
    return ""


def _gene_name(record: dict[str, object]) -> str:
    genes = record.get("genes", [])
    if genes and isinstance(genes[0], dict):
        name = genes[0].get("geneName", {})
        if isinstance(name, dict):
            return str(name.get("value", ""))
    return ""


def _score_record(record: dict[str, object], hypothesis: str, query_origin: str) -> dict[str, object]:
    comments = _comment_rows(record)
    matched = [row for row in comments if _overlap(hypothesis, str(row["text"]))]
    gene_name = _gene_name(record)
    protein_name = _protein_name(record)
    # A symbol proposed by the model may match the entry name rather than be
    # repeated in every experimental comment.  Free phenotype queries do not
    # receive that exception: unrelated comments on a generally well studied
    # protein must not be upgraded to trait evidence.
    name_match = _overlap(hypothesis, f"{gene_name} {protein_name}")
    if not matched and query_origin == "llm_known_gene" and name_match:
        matched = comments
    channels = {str(row["channel"]) for row in matched}
    pmids = sorted({pmid for row in matched for pmid in row["pmids"]})
    score = min(float(record.get("annotationScore", 0.0) or 0.0), 5.0) / 5.0
    if "DISRUPTION PHENOTYPE" in channels:
        score += 6.0
    if "MUTAGENESIS" in channels:
        score += 5.0
    if "FUNCTION" in channels:
        score += 3.0
    if "PATHWAY" in channels:
        score += 2.0
    if pmids:
        score += 1.0
    if name_match:
        score += 2.0
    lineage = set(record.get("organism", {}).get("lineage", [])) if isinstance(record.get("organism"), dict) else set()
    if "Fabaceae" in lineage:
        score += 0.75
    direct = bool(pmids) and bool(channels & {"DISRUPTION PHENOTYPE", "MUTAGENESIS"})
    functional = bool(pmids) and "FUNCTION" in channels
    evidence_level = "direct_model_genetic" if direct else "model_functional_evidence" if functional else "database_context"
    excerpts = [str(row["text"]) for row in matched if row["text"]][:2]
    return {
        "score": round(score, 4), "channels": sorted(channels), "pmids": pmids,
        "gene_name": gene_name, "protein_name": protein_name,
        "evidence_level": evidence_level, "excerpts": excerpts,
    }


def _llm_hypotheses(
    phenotype: str,
    species: str,
    direction: str,
    inheritance: str,
    out_dir: Path,
    provider: str | None,
    model: str | None,
    endpoint: str | None,
    allow_external_llm: bool,
    reuse_existing: bool = True,
) -> tuple[list[dict[str, str]], dict[str, object]]:
    if bool(provider) != bool(model):
        raise ValueError("--llm-provider and --model must be supplied together")
    if not provider:
        return [], {"used": False, "status": "not_requested"}
    out_dir.mkdir(parents=True, exist_ok=True)
    prompt_path = out_dir / "homolog_panel_prompt.jsonl"
    prompt = {
        "prompt_id": "phenotype-to-homolog-search-hypotheses",
        "gene_id": "__trait__",
        "system": PANEL_SYSTEM_PROMPT,
        "user": f"Species: {species}\nPhenotype: {phenotype}\nDirection: {direction}\nInheritance: {inheritance}",
    }
    prompt_path.write_text(json.dumps(prompt, ensure_ascii=False) + "\n", encoding="utf-8")
    llm_dir = out_dir / "llm"
    llm_dir.mkdir(parents=True, exist_ok=True)
    normalized_existing = llm_dir / "llm_candidate_responses.json"
    manifest_existing = llm_dir / "llm_run_manifest.json"
    request_meta_path = llm_dir / "llm_request_fingerprint.json"
    request_payload = {
        "prompt": prompt, "provider": provider, "model": model, "endpoint": endpoint or "",
        "allow_external_llm": bool(allow_external_llm), "max_output_tokens": 1800,
    }
    request_encoded = json.dumps(request_payload, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    request_fingerprint = hashlib.sha256(request_encoded.encode("utf-8")).hexdigest()
    previous_request = json.loads(request_meta_path.read_text(encoding="utf-8")) if request_meta_path.exists() else {}
    if (
        reuse_existing and normalized_existing.exists() and manifest_existing.exists()
        and previous_request.get("request_fingerprint") == request_fingerprint
    ):
        run = {"normalized_responses": str(normalized_existing), "manifest": str(manifest_existing)}
        reused = True
    else:
        run = run_llm_prompts(
            str(prompt_path), str(llm_dir), provider, str(model), endpoint,
            allow_external_api=allow_external_llm, max_output_tokens=1800,
        )
        reused = False
        request_meta_path.write_text(
            json.dumps({**request_payload, "request_fingerprint": request_fingerprint}, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    payload = json.loads(runtime_path(str(run["normalized_responses"])).read_text(encoding="utf-8"))
    candidate_response = payload.get("candidates", [{}])[0]
    model_json = candidate_response.get("model_json")
    parse_status = "native_json"
    if not isinstance(model_json, dict):
        model_json = _repair_hypothesis_json(str(candidate_response.get("raw_text", "")))
        parse_status = "conservative_repair" if isinstance(model_json, dict) else "unparsed"
    raw = model_json.get("search_hypotheses", []) if isinstance(model_json, dict) else []
    hypotheses = []
    for item in raw[:12] if isinstance(raw, list) else []:
        if not isinstance(item, dict):
            continue
        query = re.sub(r"\s+", " ", str(item.get("query", "")).strip())[:120]
        if not query:
            continue
        known_gene = str(item.get("known_gene", "")).strip()[:60]
        if not _known_gene_queries(known_gene):
            known_gene = ""
        hypotheses.append({
            "query": query,
            "mechanism": str(item.get("mechanism", "")).strip()[:120],
            "known_gene": known_gene,
            "rationale": str(item.get("rationale", "")).strip()[:300],
        })
    llm_run_meta = json.loads(runtime_path(str(run["manifest"])).read_text(encoding="utf-8"))
    return hypotheses, {
        "used": True, "provider": provider, "model": model,
        "external": bool(llm_run_meta.get("external_api", False)),
        "status": "parsed" if hypotheses else "no_valid_hypotheses",
        "parse_status": parse_status,
        "reused": reused,
        "request_fingerprint": request_fingerprint,
        "manifest": str(run["manifest"]),
    }


def automated_homolog_panel(
    phenotype: str,
    species: str,
    out_dir_arg: str,
    direction: str = "unknown",
    inheritance: str = "unknown",
    terms: list[str] | None = None,
    llm_provider: str | None = None,
    model: str | None = None,
    endpoint: str | None = None,
    allow_external_api: bool = False,
    allow_external_llm: bool = False,
    panel_size: int = 20,
    per_query: int = 10,
) -> dict[str, object]:
    """Create an unreviewed, sequence-backed reference panel without human triage."""
    if not allow_external_api:
        raise ValueError("UniProt search is external; pass --allow-external-api after reviewing phenotype transmission")
    if not phenotype.strip() or not species.strip():
        raise ValueError("phenotype and species are required")
    if panel_size < 1 or per_query < 1:
        raise ValueError("panel-size and per-query must be positive")
    out_dir = runtime_path(out_dir_arg)
    raw_dir = out_dir / "uniprot_raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    hypotheses, llm_meta = _llm_hypotheses(
        phenotype, species, direction, inheritance, out_dir,
        llm_provider, model, endpoint, allow_external_llm,
    )
    search_rows: list[dict[str, str]] = []
    for phrase in [phenotype, *(terms or [])]:
        search_rows.extend(_mechanical_queries(phrase))
    # Canonical model-gene symbols are the most precise semantic output, so
    # allocate the fixed query budget to them before longer mechanism phrases.
    for item in hypotheses:
        for gene_query in _known_gene_queries(item["known_gene"]):
            search_rows.append({"query": gene_query, "origin": "llm_known_gene", **{key: item[key] for key in ("mechanism", "known_gene", "rationale")}})
    for item in hypotheses:
        search_rows.append({"query": item["query"], "origin": "llm", **{key: item[key] for key in ("mechanism", "known_gene", "rationale")}})
    deduplicated: list[dict[str, str]] = []
    seen_terms: set[str] = set()
    for row in search_rows:
        key = row["query"].casefold()
        if key not in seen_terms:
            seen_terms.add(key)
            deduplicated.append(row)
    search_rows = deduplicated[:24]
    pd.DataFrame(search_rows).to_csv(out_dir / "search_hypotheses.tsv", sep="\t", index=False)

    audit_rows: list[dict[str, object]] = []
    raw_files: list[str] = []
    for index, hypothesis in enumerate(search_rows, 1):
        query = f'(reviewed:true) AND (taxonomy_id:33090) AND ("{hypothesis["query"].replace(chr(34), " ")}")'
        url = UNIPROT_ENDPOINT + "?" + urlencode({"query": query, "format": "json", "size": per_query})
        payload = _http_json(url)
        raw_path = raw_dir / f"{index:02d}_{safe_id(hypothesis['query'])}.json"
        raw_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        raw_files.append(str(raw_path))
        for record in payload.get("results", []):
            if not isinstance(record, dict):
                continue
            organism = record.get("organism", {}) if isinstance(record.get("organism"), dict) else {}
            if str(organism.get("scientificName", "")).casefold() == species.casefold():
                continue
            scored = _score_record(record, hypothesis["query"], hypothesis["origin"])
            sequence = record.get("sequence", {}) if isinstance(record.get("sequence"), dict) else {}
            accession = str(record.get("primaryAccession", ""))
            if not accession or not sequence.get("value"):
                continue
            audit_rows.append({
                "reference_gene_id": accession,
                "reference_gene_name": scored["gene_name"] or str(record.get("uniProtkbId", "")),
                "reference_species": str(organism.get("scientificName", "")),
                "protein_name": scored["protein_name"],
                "query": hypothesis["query"], "query_origin": hypothesis["origin"],
                "mechanism": hypothesis["mechanism"], "rationale": hypothesis["rationale"],
                "auto_score": scored["score"], "evidence_level": scored["evidence_level"],
                "evidence_channels": ";".join(scored["channels"]),
                "source_id": ";".join([*(f"PMID:{value}" for value in scored["pmids"]), f"UniProt:{accession}"]),
                "evidence_excerpt": " | ".join(scored["excerpts"]),
                "uniprot_entry_version": str(record.get("entryAudit", {}).get("entryVersion", "")),
                "annotation_update": str(record.get("entryAudit", {}).get("lastAnnotationUpdateDate", "")),
                "sequence_length": int(sequence.get("length", len(str(sequence.get("value", ""))))),
                "sequence": str(sequence.get("value", "")),
            })
    audit = pd.DataFrame(audit_rows)
    if audit.empty:
        raise RuntimeError("UniProt search returned no sequence-backed plant homolog proposals")
    audit = audit.sort_values(["auto_score", "reference_gene_id"], ascending=[False, True], kind="stable")
    audit.to_csv(out_dir / "homolog_proposal_audit.tsv", sep="\t", index=False)

    selected_rows: list[pd.Series] = []
    selected_ids: set[str] = set()
    for _, group in audit.groupby("query", sort=False):
        eligible = group[group["evidence_level"].ne("database_context")]
        if not eligible.empty:
            row = eligible.iloc[0]
            if row["reference_gene_id"] not in selected_ids:
                selected_rows.append(row)
                selected_ids.add(str(row["reference_gene_id"]))
                if len(selected_rows) >= panel_size:
                    break
    if len(selected_rows) < panel_size:
        for _, row in audit.iterrows():
            accession = str(row["reference_gene_id"])
            if accession in selected_ids:
                continue
            selected_rows.append(row)
            selected_ids.add(accession)
            if len(selected_rows) >= panel_size:
                break
    selected = pd.DataFrame(selected_rows)
    discovery = "uniprot_api" if not llm_meta["used"] else f"{llm_provider}+uniprot_api"
    panel_rows = []
    records = []
    for row in selected.itertuples():
        relation = row.evidence_excerpt or f"UniProt search match for {row.query}"
        panel_rows.append({
            "reference_gene_id": row.reference_gene_id,
            "reference_gene_name": row.reference_gene_name,
            "reference_species": row.reference_species,
            "trait_relation": relation,
            "evidence_level": row.evidence_level,
            "source_id": row.source_id,
            "discovery_method": discovery,
            "review_status": "pending_review",
            "notes": "Automatically selected before interval inspection; screening hypothesis only.",
            "auto_score": row.auto_score,
            "search_hypotheses": row.query,
            "evidence_channels": row.evidence_channels,
            "uniprot_entry_version": row.uniprot_entry_version,
            "sequence_length": row.sequence_length,
        })
        records.append(SeqRecord(Seq(str(row.sequence)), id=str(row.reference_gene_id), description=f"{row.reference_gene_name} {row.reference_species}"))
    panel_path = out_dir / "reference_gene_panel.tsv"
    pd.DataFrame(panel_rows, columns=PANEL_COLUMNS).to_csv(panel_path, sep="\t", index=False)
    fasta_path = out_dir / "reference_gene_panel.faa"
    SeqIO.write(records, str(fasta_path), "fasta")
    review_path = out_dir / "reference_gene_panel_review.tsv"
    pd.DataFrame(panel_rows, columns=PANEL_COLUMNS).to_csv(review_path, sep="\t", index=False)
    manifest = {
        "schema_version": "1.0", "panel_version": PANEL_VERSION,
        "phenotype": phenotype, "species": species, "direction": direction, "inheritance": inheritance,
        "input_terms": list(terms or []), "panel_size_requested": panel_size, "per_query": per_query,
        "search_hypothesis_count": len(search_rows), "proposal_count": len(audit), "panel_count": len(panel_rows),
        "panel": str(panel_path), "panel_fasta": str(fasta_path), "review": str(review_path),
        "llm": llm_meta, "uniprot_endpoint": UNIPROT_ENDPOINT,
        "external_transmission": {"phenotype": True, "qtl_coordinates": False, "candidate_ids": False, "sequences": False},
        "raw_api_files": raw_files,
        "review_status": "pending_review", "causal_probability": False,
        "selection_policy": "one trait-linked experimentally annotated hit per search hypothesis, then score fill; no interval inspection",
        "warning": "Automated panel membership is weak screening evidence until a person reviews the final Stage-A set.",
        "created_at": datetime.now(UTC).isoformat(),
    }
    manifest_path = out_dir / "homolog_panel_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return {"panel": str(panel_path), "fasta": str(fasta_path), "audit": str(out_dir / "homolog_proposal_audit.tsv"), "review": str(review_path), "manifest": str(manifest_path), "panel_count": len(panel_rows)}


def automated_homolog_match(
    candidates_arg: str,
    candidate_proteins_arg: str,
    panel_arg: str,
    panel_fasta_arg: str,
    out_dir_arg: str,
    base_evidence_arg: str | None = None,
    target_size: int = 20,
    runner: str = "auto",
    threads: int = 4,
) -> dict[str, object]:
    """Match an automated panel to interval proteins and optionally freeze Stage A."""

    candidates = runtime_path(candidates_arg)
    candidate_proteins = runtime_path(candidate_proteins_arg)
    panel_path = runtime_path(panel_arg)
    panel_fasta = runtime_path(panel_fasta_arg)
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    panel = pd.read_csv(panel_path, sep="\t", dtype=str).fillna("")
    if not {"reference_gene_id", "source_id"}.issubset(panel.columns):
        raise ValueError("panel requires reference_gene_id and source_id")
    source_map = dict(zip(panel["reference_gene_id"], panel["source_id"]))
    confidence_map = {
        str(row.reference_gene_id): "low" if str(row.evidence_level) == "database_context" else "medium"
        for row in panel.itertuples()
    }
    blast = _run_interval_blastp(panel_fasta, candidate_proteins, out_dir / "panel_vs_interval.tsv", runner, threads)
    evidence_result = homolog_blast_evidence(
        str(candidates), str(blast["output"]), str(out_dir / "automated_homolog_evidence.tsv"),
        "automated_homolog_panel", evidence_code="automated_homolog_match",
        query_source_map=source_map, query_confidence_map=confidence_map, max_confidence="medium",
    )
    result: dict[str, object] = {"evidence": evidence_result["output"], "audit": evidence_result["audit"], "blast": blast["output"]}
    if base_evidence_arg:
        base_path = runtime_path(base_evidence_arg)
        base = pd.read_csv(base_path, sep="\t", dtype=str).fillna("")
        automatic = pd.read_csv(evidence_result["output"], sep="\t", dtype=str).fillna("")
        combined_path = out_dir / "combined_evidence.tsv"
        pd.concat([base, automatic], ignore_index=True).drop_duplicates().to_csv(combined_path, sep="\t", index=False)
        rank = prioritize_evidence(str(candidates), str(combined_path), str(out_dir / "rank"))
        shortlist = screening_shortlist(str(rank["ranked_candidates"]), str(out_dir / "stage_a"), target_size=target_size)
        result.update({"combined_evidence": str(combined_path), "ranked": rank["ranked_candidates"], "shortlist": shortlist["shortlist"]})
    manifest = {
        "schema_version": "1.0", "panel": file_fingerprint(panel_path, False),
        "panel_fasta": file_fingerprint(panel_fasta, False), "candidates": file_fingerprint(candidates, False),
        "candidate_proteins": file_fingerprint(candidate_proteins, False),
        "blast": blast, "evidence_code": "automated_homolog_match",
        "llm_writes_gene_score": False, "review_status": "pending_review", "causal_probability": False,
        "warning": "Automated homolog similarity is one comparative-function family and cannot establish orthology or causality.",
    }
    manifest_path = out_dir / "automated_homolog_match_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    result["manifest"] = str(manifest_path)
    return result


def benchmark_automated_homolog_panels(
    reference_manifest_arg: str,
    out_dir_arg: str,
    species: str = "Pisum sativum",
    case_ids: list[str] | None = None,
    llm_provider: str | None = None,
    model: str | None = None,
    endpoint: str | None = None,
    allow_external_api: bool = False,
    allow_external_llm: bool = False,
    panel_size: int = 20,
    per_query: int = 10,
    target_size: int = 20,
    runner: str = "auto",
    threads: int = 4,
    resume: bool = False,
    use_curated_terms: bool = False,
) -> dict[str, object]:
    """Blindly benchmark phenotype-to-panel retrieval on frozen reference cases."""
    from .homolog_shared import evaluate_rank as _evaluate_rank

    manifest_path = runtime_path(reference_manifest_arg)
    reference = json.loads(manifest_path.read_text(encoding="utf-8"))
    config_path = runtime_path(str(reference.get("config", "")))
    if not config_path.exists() and reference.get("config"):
        config_path = Path(__file__).resolve().parents[2] / str(reference["config"])
    config = pd.read_csv(config_path, sep="\t", dtype=str).fillna("")
    if not {"case_id", "trait", "answer_gene_id"}.issubset(config.columns):
        raise ValueError("reference benchmark config requires case_id, trait, and answer_gene_id")
    config_by_case = config.set_index("case_id", drop=False)
    selected_cases = [
        case for case in reference.get("cases", [])
        if isinstance(case, dict)
        and case.get("term_profile") == "phenotype_only"
        and (not case_ids or str(case.get("locus_id")) in set(case_ids))
    ]
    if not selected_cases:
        raise ValueError("no phenotype_only benchmark cases matched --case")
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    failures: list[dict[str, str]] = []
    for case in selected_cases:
        locus_id = str(case["locus_id"])
        if locus_id not in config_by_case.index:
            failures.append({"case_id": locus_id, "error": "case absent from benchmark config"})
            continue
        case_config = config_by_case.loc[locus_id]
        if isinstance(case_config, pd.DataFrame):
            case_config = case_config.iloc[0]
        phenotype = str(case_config["trait"])
        phenotype_terms = [value for value in str(case_config.get("phenotype_terms", "")).split("|") if value] if use_curated_terms else []
        case_dir = out_dir / safe_id(locus_id)
        try:
            panel_dir = case_dir / "panel"
            panel_manifest_path = panel_dir / "homolog_panel_manifest.json"
            panel_reused = False
            if resume and panel_manifest_path.exists() and (panel_dir / "reference_gene_panel.tsv").exists() and (panel_dir / "reference_gene_panel.faa").exists():
                previous_panel = json.loads(panel_manifest_path.read_text(encoding="utf-8"))
                expected_llm = bool(llm_provider)
                previous_llm = previous_panel.get("llm", {})
                if (
                    previous_panel.get("phenotype") == phenotype
                    and previous_panel.get("species") == species
                    and previous_panel.get("input_terms", []) == phenotype_terms
                    and int(previous_panel.get("panel_size_requested", -1)) == panel_size
                    and int(previous_panel.get("per_query", -1)) == per_query
                    and bool(previous_llm.get("used")) == expected_llm
                    and (not expected_llm or (
                        previous_llm.get("provider") == llm_provider
                        and previous_llm.get("model") == model
                    ))
                ):
                    panel = {
                        "panel": str(panel_dir / "reference_gene_panel.tsv"),
                        "fasta": str(panel_dir / "reference_gene_panel.faa"),
                        "manifest": str(panel_manifest_path),
                        "panel_count": int(previous_panel.get("panel_count", 0)),
                    }
                    panel_reused = True
            if not panel_reused:
                panel = automated_homolog_panel(
                    phenotype, species, str(panel_dir), terms=phenotype_terms,
                    llm_provider=llm_provider, model=model, endpoint=endpoint,
                    allow_external_api=allow_external_api, allow_external_llm=allow_external_llm,
                    panel_size=panel_size, per_query=per_query,
                )
            match_dir = case_dir / "match"
            match_manifest_path = match_dir / "automated_homolog_match_manifest.json"
            ranked_path = match_dir / "rank" / "ranked_candidates.tsv"
            shortlist_path = match_dir / "stage_a" / "deterministic_shortlist.tsv"
            evidence_path = match_dir / "automated_homolog_evidence.tsv"
            match_reused = bool(resume and match_manifest_path.exists() and ranked_path.exists() and shortlist_path.exists() and evidence_path.exists())
            if match_reused:
                match = {
                    "ranked": str(ranked_path), "shortlist": str(shortlist_path),
                    "evidence": str(evidence_path), "manifest": str(match_manifest_path),
                }
            else:
                match = automated_homolog_match(
                    str(case["extract"]["outputs"]["candidates_tsv"]),
                    str(case["official_protein_subset"]["output"]),
                    str(panel["panel"]), str(panel["fasta"]), str(match_dir),
                    str(case["evidence"]["evidence"]), target_size, runner, threads,
                )
            answer = str(case_config["answer_gene_id"])
            evaluation = _evaluate_rank(runtime_path(str(match["ranked"])), answer)
            shortlist = pd.read_csv(runtime_path(str(match["shortlist"])), sep="\t", dtype=str).fillna("")
            accepted = {value.strip() for value in answer.split("|") if value.strip()}
            selected_truth = sorted(accepted.intersection(set(shortlist["gene_id"])))
            automatic_evidence = pd.read_csv(runtime_path(str(match["evidence"])), sep="\t", dtype=str).fillna("")
            truth_automatic = sorted(accepted.intersection(set(automatic_evidence["gene_id"])))
            rows.append({
                "case_id": locus_id,
                "trait": phenotype,
                "answer_gene_id": answer,
                "panel_count": panel["panel_count"],
                "homolog_supported_gene_count": len(automatic_evidence),
                "truth_automated_homolog_supported": bool(truth_automatic),
                "truth_automated_homolog_gene_id": "|".join(truth_automatic),
                "rank_min": evaluation["rank_min"],
                "rank_max": evaluation["rank_max"],
                "top20_conservative": evaluation["top20"],
                "in_fixed_stage_a_20": bool(selected_truth),
                "matched_answer_gene_id": "|".join(selected_truth),
                "score_breakdown": evaluation.get("score_breakdown", ""),
                "panel_manifest": panel["manifest"],
                "match_manifest": match["manifest"],
                "shortlist": match["shortlist"],
                "panel_reused": panel_reused,
                "match_reused": match_reused,
            })
        except Exception as error:
            failures.append({"case_id": locus_id, "error": str(error)})
    summary = pd.DataFrame(rows)
    summary_path = out_dir / "automated_homolog_benchmark.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    failure_path = out_dir / "benchmark_failures.tsv"
    pd.DataFrame(failures, columns=["case_id", "error"]).to_csv(failure_path, sep="\t", index=False)
    completed = len(rows)
    fixed_recall = (sum(bool(row["in_fixed_stage_a_20"]) for row in rows) / completed) if completed else 0.0
    conservative_recall = (sum(bool(row["top20_conservative"]) for row in rows) / completed) if completed else 0.0
    truth_homolog_recall = (sum(bool(row["truth_automated_homolog_supported"]) for row in rows) / completed) if completed else 0.0
    run_manifest = {
        "schema_version": "1.0", "benchmark_type": "blind_automated_homolog_selection",
        "reference_manifest": file_fingerprint(manifest_path, False),
        "species": species, "requested_case_ids": case_ids or [],
        "curated_phenotype_terms_used": use_curated_terms,
        "resume_requested": resume,
        "provider_lane": "no_llm" if not llm_provider else llm_provider,
        "llm_provider": llm_provider, "model": model, "endpoint": endpoint,
        "llm_role": "search_hypothesis_generation_only",
        "answer_gene_used_before_evaluation": False,
        "candidate_ids_sent_to_llm_or_uniprot": False,
        "case_count": len(selected_cases), "completed_case_count": completed,
        "failed_case_count": len(failures), "fixed_stage_a_20_recall": round(fixed_recall, 6),
        "conservative_top20_recall": round(conservative_recall, 6),
        "truth_automated_homolog_recall": round(truth_homolog_recall, 6),
        "summary": str(summary_path), "failures": str(failure_path),
        "warning": "This measures retrieval on known loci; it does not estimate causal probability for unknown traits.",
        "created_at": datetime.now(UTC).isoformat(),
    }
    run_manifest_path = out_dir / "automated_homolog_benchmark_manifest.json"
    run_manifest_path.write_text(json.dumps(run_manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return run_manifest | {"manifest": str(run_manifest_path)}
