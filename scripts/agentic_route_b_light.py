"""Light Route B research handoff (#1004).

An opt-in replacement for the research ledger plus formal transcription stages.
The original researcher judges every interval gene against the frozen reference
material and records one Tier, one reason and the paper checklist. The machine
checks only hallucination: the gene id belongs to the batch, the reference
endpoint is frozen for that gene, the Tier is in the vocabulary, and every
paper resolves with a matching title and a retrievable abstract. The output is
converted back into the existing artifacts (discoveries, llm-run agent file) so
that semantic import, source verification, transfer verification, candidate
freeze and Wet remain unchanged.
"""
from __future__ import annotations

import copy
import difflib
import hashlib
import html
import json
import re
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

CONTRACT = "agentic-route-b-light-1.0.0"
V2_CONTRACT = "agentic-route-b-light-2.0.0"
SCIENCE_CONTRACT = "agentic-route-b-science-2.0.0"


def _paper_failure_state(dropped: list[dict[str, str]]) -> str:
    retrieval_states = {"not_attempted", "provider_failure", "not_resolved", "missing_abstract", "cache_entry_invalid"}
    return "retrieval_failed" if dropped and all(item["state"] in retrieval_states for item in dropped) else "invalid"
BATCH_ID_FORMAT = "b{:02d}"
TIERS = ("T1", "T2", "T3", "T4", "T5", "T6", "U", "not_assessed")
PAPER_TIERS = ("T1", "T2", "T3")
MAX_PAPERS = 5
HTTP_TRIES = 3
HTTP_BACKOFF_SECONDS = 1.0
STATE_ORDER = ("T1", "T2", "T3", "T4", "T5", "T6", "U")
# A recorded function and an annotation remain different evidence types.
# Both need an actual frozen source and the original concrete phenotype bridge.
NO_PAPER_DOWNSTREAM_TIER = {"T4": "T4", "T5": "T5"}
TITLE_THRESHOLD = 0.85
EUROPE_PMC_REST = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
CROSSREF_WORKS = "https://api.crossref.org/works/"
NCBI_IDCONV = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
USER_AGENT = "cqtl-light-check/0.1 (local research)"
CHECKLIST_FIELDS = (
    "focal_species", "same_gene_as_focal", "phenotype_reported",
    "equivalent_target_phenotype", "similar_cellular_or_developmental_phenotype",
    "shared_specific_process", "function_supports_hypothesis",
    "annotation_supports_hypothesis", "identity_resolved",
)
CARD_FIELDS = (
    "gene_symbol", "protein_name", "function_description", "mutant_phenotype", "pathway",
    "go_terms", "keywords", "tissue_specificity", "developmental_stage", "protein_family",
    "literature_source_ids", "domain_description", "domain_ids", "evidence_codes", "evidence_type", "source_id",
)

TASK = """# Route B — per-gene Tier judgment (light format)

You are an independent Route B researcher. The focal species is `{FOCAL_SPECIES}`. The target phenotype is `{PHENOTYPE}`.

Judge every gene in `material.md` in this directory. That file is the frozen reference material, one section per interval gene with its unit id. Only `material.md` and the assigned `material.full.json` are research inputs in this directory tree. Do not use subagents.

Use these absolute paths regardless of the tool's working directory:
- assigned material: `{MATERIAL_PATH}`
- complete raw material (only open the needed unit/field if provenance or details are needed): `{FULL_MATERIAL_PATH}`
- output: `{RESULTS_PATH}`
- assigned interval genes: {ASSIGNED_GENES}
Do not glob, search for another material.md, or open sibling batches. If the assigned
path is unavailable, report that exact failure and stop; do not substitute another batch.
Shared biological cards in this material are untruncated and appear once; follow
each gene's card IDs. The raw JSON preserves all source/projection fields; do not
dump the full JSON into the conversation or re-read the full batch to find one field.

## What to do for each gene
1. Read that gene's reference material. Decide whether its known function bridges to the phenotype.
   - Judge each gene (copy) individually; shared material may be reused.
   - Do not transfer a judgment by family or name similarity alone.
   - No keyword-only exclusion.
   - Unknown is not "no". Not investigated is not negative.
2. Where a bridge is possible (a candidate for T1–T5), search the literature: PubMed / Europe PMC / PMC via web fetch.
   - Collect papers up to about 20 per original reference gene/shared mechanism as far as possible, not per copy. Originals and reviews are both fine.
   - From them, keep up to 5 useful papers, strongest Tier first.
   - If the literature runs out, write "literature exhausted". No further reason is needed.
   - A bridge is often only visible through other species. When a gene's function could relate to the phenotype, look for papers on its homologs in other plants that report a similar phenotype (T2) or a shared specific process (T3), not only papers on the pea gene.
   - Genes without a possible bridge need no literature.
3. Answer the 9-question checklist for each kept paper with yes / no / unknown:
   - `focal_species`
   - `same_gene_as_focal`
   - `phenotype_reported`
   - `equivalent_target_phenotype`
   - `similar_cellular_or_developmental_phenotype`
   - `shared_specific_process`
   - `function_supports_hypothesis`
   - `annotation_supports_hypothesis`
   - `identity_resolved`
4. The paper Tier follows this table:
   - When `phenotype_reported` = yes:
     - T1 = `focal_species` is yes and (`equivalent_target_phenotype` or `similar_cellular_or_developmental_phenotype`) is yes
     - otherwise T2 = `equivalent_target_phenotype` or `similar_cellular_or_developmental_phenotype` is yes
     - otherwise T3 = `shared_specific_process` is yes
   - Otherwise:
     - T4 = `function_supports_hypothesis` is yes
     - T5 = `annotation_supports_hypothesis` is yes
     - T6 = `identity_resolved` is yes
     - U = none of the above
5. The gene Tier is the Tier of the strongest kept paper. Without a paper, use T4/T5/T6 from the material itself, or U.

## Output (the only file you write): `results.jsonl` in the batch directory
Paths are relative to the batch directory that contains this task and `material.md`.
Write one JSON object per gene in `material.md`. Keys:
- `interval_gene_id`
- `reference_gene`: the exact `ref` accession from `material.md` that the judgment relies on
- `tier`: `T1` `T2` `T3` `T4` `T5` `T6` `U`, or `not_assessed`
- `function`: one concrete phrase for the reference function relied on
- `reason`: one or two sentences linking function to phenotype
- `papers`: observations from up to 5 distinct papers, strongest first; `[]` when none. Keep multiple observations from the same paper separate when gene, phenotype, experiment, tissue/stage, finding or checklist differs. Each entry:
  - `paper_id` (`PMID:…`, `PMCID:…` or `DOI:…`)
  - `title`: copied exactly from the source
  - `reported_gene`, `reported_gene_id` (original namespace, or explicit `unknown` if unavailable), and `reported_species`; never substitute the reference/proxy ID for the paper's gene ID
  - `gene_correspondence`: `same_gene` / `homolog` / `different_gene` / `unknown`
  - `finding`: one sentence on what the paper reports
  - `checklist`: the 9 keys → yes / no / unknown
  - when available: `reported_phenotype`, `reported_process`, `reported_organ`, `reported_stage`, `source_excerpt`, `source_location`, `unknown_reason`, `publication_type`
- `literature`: one of `"none needed"`, `"exhausted"` or `"screened N"`

To avoid repeating titles, optionally write one first JSONL line:
{"bibliography":{"p1":{"source_id":"PMID:...","source_title":"exact title copied by you"}}}
Then an observation may use `"paper_ref":"p1"` instead of paper_id/title.
Only paper ID/title are shared. Write the reported gene name/ID, species, finding,
correspondence and checklist separately on EVERY observation. A review's different
genes/phenotypes must never inherit one another's facts or Tier. No extra matching
research is needed: when in doubt, keep separate full observations.

`not_assessed` means the gene was not judged. `U` means the gene was judged and no bridge was found.
Unknown is not "no". Do not invent IDs or titles.

A machine will check afterwards that each `paper_id` exists and its title matches, and that `interval_gene_id` and `reference_gene` match the material. A paper that fails is dropped; an unverified gene keeps no claim.

## Discipline
Keep tool output short and do not print long JSON. Your final reply gives only the path and the number of genes written.
"""


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_text(path: Path, text: str) -> None:
    """Write UTF-8 with LF on every platform so the recorded SHA-256 is the file's SHA-256."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(text.encode("utf-8"))


def _write_json(path: Path, value: Any) -> None:
    _write_text(path, json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha(path: Path) -> str:
    return _sha256(path.read_bytes())


def _tokens(value: object) -> list[str]:
    found: list[str] = []
    seen: set[str] = set()
    for token in str(value or "").replace(",", ";").split(";"):
        token = token.strip()
        if token and token.casefold() not in seen:
            found.append(token)
            seen.add(token.casefold())
    return found


def _plain(value: object) -> str:
    """Strip markup tags, unescape entities and collapse whitespace (provider titles carry <i>, &amp;)."""
    text = html.unescape(re.sub(r"<[^>]+>", "", str(value or "")))
    return " ".join(text.split())


def _title_key(value: object) -> str:
    text = unicodedata.normalize("NFKC", _plain(value)).casefold()
    return " ".join("".join(char if char.isalnum() else " " for char in text).split())


def _title_similarity(applied: object, resolved: object) -> float:
    return difflib.SequenceMatcher(None, _title_key(applied), _title_key(resolved), autojunk=False).ratio()


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _citation(paper_id: object) -> tuple[str, str]:
    from candidate_qtl.source_verify import classify_citation

    kind, identifier = classify_citation(str(paper_id or ""))
    return (kind, identifier) if kind in {"pmid", "pmcid", "doi"} else (kind, str(identifier or "").strip())


def _cache_key(paper_id: object) -> str:
    kind, identifier = _citation(paper_id)
    if kind == "doi":
        identifier = identifier.casefold()
    elif kind == "pmcid":
        identifier = identifier.upper()
    return f"{kind}:{identifier}"


def _cache_aliases(paper_id: object, record: dict[str, Any]) -> list[str]:
    aliases = []
    for kind, value in (("pmid", record.get("pmid")), ("pmcid", record.get("pmcid")), ("doi", record.get("doi"))):
        if str(value or "").strip():
            aliases.append(_cache_key(f"{kind}:{value}"))
    if not aliases:
        aliases.append(_cache_key(paper_id))
    return list(dict.fromkeys(aliases))


_LABELS = {"pmid": "PMID", "pmcid": "PMCID", "doi": "DOI"}


def _canonical_id(kind: str, identifier: str, record: dict[str, Any] | None = None) -> str:
    """PMID:n / PMCID:PMCn / DOI:x, with the DOI spelled as the provider record spells it."""
    if kind == "doi" and record and str(record.get("doi") or "").strip().casefold() == identifier.casefold():
        identifier = str(record["doi"]).strip()
    return f"{_LABELS[kind]}:{identifier.upper() if kind == 'pmcid' else identifier}"


def _record_tokens(record: dict[str, Any]) -> list[str]:
    tokens = []
    for kind in ("pmid", "pmcid", "doi"):
        value = str(record.get(kind) or "").strip()
        if value:
            tokens.append(_canonical_id(kind, value))
    return tokens


# ---------------------------------------------------------------------------
# bundle validation and material rendering
# ---------------------------------------------------------------------------

def _load_bundle(path_arg: str, expected_sha256: str) -> tuple[Path, dict[str, Any]]:
    path = Path(path_arg).resolve()
    raw = path.read_bytes()
    if not expected_sha256 or _sha256(raw) != expected_sha256:
        raise ValueError("light bundle SHA-256 mismatch")
    bundle = json.loads(raw)
    if (not isinstance(bundle, dict) or bundle.get("stage") != "post_b0_semantic_bundle"
            or bundle.get("route") != "interval_first" or not str(bundle.get("case_id", "")).strip()):
        raise ValueError("light requires a post-B0 semantic bundle")
    if not isinstance(bundle.get("phenotype"), str) or not bundle["phenotype"].strip():
        raise ValueError("light bundle phenotype is required")
    units = bundle.get("units")
    if not isinstance(units, list):
        raise ValueError("light bundle requires units")
    seen_units: set[str] = set()
    seen_genes: set[str] = set()
    for unit in units:
        if not isinstance(unit, dict):
            raise ValueError("light bundle contains an invalid unit")
        unit_id, members, genes = unit.get("unit_id"), unit.get("member_interval_genes"), unit.get("genes")
        if not isinstance(unit_id, str) or not unit_id.strip() or unit_id in seen_units:
            raise ValueError("light bundle requires unique non-empty unit IDs")
        seen_units.add(unit_id)
        if (not isinstance(members, list) or not members
                or any(not isinstance(gene, str) or not gene.strip() for gene in members)
                or len(set(members)) != len(members) or seen_genes.intersection(members)):
            raise ValueError("light bundle requires unique non-empty member gene IDs")
        if (not isinstance(genes, list)
                or not all(isinstance(gene, dict) and isinstance(gene.get("gene_id"), str) for gene in genes)
                or sorted(gene.get("gene_id", "") for gene in genes) != sorted(members)):
            raise ValueError("light bundle gene records differ from the member roster")
        seen_genes.update(members)
    return path, bundle


def _requested_units(bundle: dict[str, Any]) -> list[dict[str, Any]]:
    return [unit for unit in bundle["units"] if str(unit.get("request_state", "requested")) == "requested"]


def _unit_gene(unit: dict[str, Any], gene_id: str) -> dict[str, Any]:
    return next(gene for gene in unit["genes"] if gene["gene_id"] == gene_id)


def _cards_for_gene(unit: dict[str, Any], gene_id: str) -> list[dict[str, Any]]:
    cards = unit.get("reference_material_cards", [])
    if not isinstance(cards, list):
        raise ValueError("light bundle reference_material_cards must be a list")
    bound = []
    for card in cards:
        if not isinstance(card, dict):
            raise ValueError("light bundle reference material card is invalid")
        owner = str(card.get("interval_gene_id") or "").strip()
        if owner in {"", gene_id}:
            bound.append(card)
    return bound


def _literature_links(unit: dict[str, Any]) -> list[dict[str, Any]]:
    rows = unit.get("literature_link_candidates", [])
    if not isinstance(rows, list) or not all(isinstance(row, dict) for row in rows):
        raise ValueError("light bundle literature_link_candidates must be a list of records")
    return rows


def _gene_endpoints(unit: dict[str, Any], gene: dict[str, Any]) -> list[str]:
    """The exact endpoints the semantic importer grounds against."""
    values: list[object] = [
        gene.get("reference_gene"), gene.get("orthology_reference_gene"),
        gene.get("family_index_subject_endpoints"),
    ]
    hits = gene.get("sequence_similarity_hits", [])
    if not isinstance(hits, list) or not all(isinstance(hit, dict) for hit in hits):
        raise ValueError("light bundle sequence_similarity_hits must be a list of records")
    values.extend(hit.get("hit_accession") for hit in hits)
    values.extend(card.get("reference_gene") for card in _cards_for_gene(unit, gene["gene_id"]))
    values.extend(row.get("reference_gene") for row in _literature_links(unit))
    return _tokens(";".join(str(value) for value in values if value))


def _gene_names(unit: dict[str, Any], gene: dict[str, Any]) -> list[str]:
    values = [gene.get("gene_symbol"), gene.get("reference_gene")]
    for card in _cards_for_gene(unit, gene["gene_id"]):
        values.extend((card.get("gene_symbol"), card.get("reference_gene")))
    return _tokens(";".join(str(value) for value in values if value))


def _species(value: object) -> str:
    """One spelling for species everywhere: drop a trailing common-name parenthesis."""
    return " ".join(str(value or "").split(" (")[0].split())


def _endpoint_species(unit: dict[str, Any], gene: dict[str, Any]) -> dict[str, str]:
    mapping: dict[str, str] = {}

    def add(references: object, species: object) -> None:
        for reference in _tokens(references):
            if reference not in mapping or (species and not mapping[reference]):
                mapping[reference] = _species(species)

    orthology = str(gene.get("orthology_reference_gene", "")).strip()
    same_reference = str(gene.get("reference_gene", "")).strip() in {"", orthology}
    add(orthology, gene.get("orthology_reference_species") or gene.get("orthology_species")
        or (gene.get("reference_species", "") if same_reference else ""))
    add(gene.get("reference_gene"), gene.get("reference_species"))
    for hit in gene.get("sequence_similarity_hits", []) or []:
        if isinstance(hit, dict):
            add(hit.get("hit_accession"), hit.get("hit_species"))
    for card in _cards_for_gene(unit, gene["gene_id"]):
        add(card.get("reference_gene"), card.get("reference_species"))
    for row in _literature_links(unit):
        add(row.get("reference_gene"), row.get("reference_species"))
    return {reference: species for reference, species in mapping.items() if species}


def _roster(unit: dict[str, Any]) -> dict[str, dict[str, Any]]:
    roster = {}
    for gene_id in unit["member_interval_genes"]:
        gene = _unit_gene(unit, gene_id)
        roster[gene_id] = {
            "unit_id": unit["unit_id"],
            "reference_endpoints": _gene_endpoints(unit, gene),
            "names": _gene_names(unit, gene),
        }
    return roster


def _render_material(batch_id: str, units: list[dict[str, Any]]) -> str:
    lines = [f"# Material for batch {batch_id}", ""]
    cards_by_key: dict[str, dict[str, Any]] = {}
    card_ids: dict[str, str] = {}
    for unit in units:
        if unit.get("identity"):
            lines.extend([
                f"# Original unit identity: {unit['unit_id']}",
                json.dumps(unit["identity"], sort_keys=True, ensure_ascii=False),
                "",
            ])
        for gene_id in unit["member_interval_genes"]:
            gene = _unit_gene(unit, gene_id)
            lines.append(f"## {gene_id}")
            lines.append(f"- unit_id: {unit['unit_id']}")
            if str(gene.get("reference_gene", "")).strip():
                lines.append(f"- gene_reference_gene: {gene['reference_gene']}")
            if str(gene.get("orthology_reference_gene", "")).strip():
                lines.append(f"- gene_orthology_reference_gene: {gene['orthology_reference_gene']}")
            cards = _cards_for_gene(unit, gene_id)
            if not cards:
                lines.append("- (no reference card)")
            for card in cards:
                # Share only exactly equal displayed biology/reference/source
                # fields. Per-gene provenance stays intact in material.full.json.
                content = {key: card[key] for key in ("reference_gene", "reference_species", *CARD_FIELDS) if key in card}
                key = _sha256(json.dumps(content, sort_keys=True, ensure_ascii=False).encode("utf-8"))
                cards_by_key[key] = content
                card_ids.setdefault(key, f"card-{len(card_ids) + 1}")
                lines.append(f"- complete biological card: {card_ids[key]}")
            # Gene-specific context must not disappear behind a shared card.
            for field in (
                "reference_function", "reference_mutant_phenotype", "reference_pathway",
                "family_id", "family_label", "reference_family_id", "evidence_completeness", "notes", *CARD_FIELDS,
            ):
                if gene.get(field):
                    lines.append(f"- gene {field}: {gene[field]}")
            lines.append("")
    lines.append("# Shared biological reference cards (no field truncation)")
    for key, card in cards_by_key.items():
        lines.extend([f"## {card_ids[key]}", f"- ref {card.get('reference_gene', '')} [{_species(card.get('reference_species', ''))}]"])
        for field in CARD_FIELDS:
            value = card.get(field)
            if value:
                lines.append(f"- {field}: " + (value if isinstance(value, str) else json.dumps(value, ensure_ascii=False)))
        lines.append("")
    return "\n".join(lines)


SCIENCE_TASK = """# Route B science-only output v2
Focal species: {FOCAL_SPECIES}. Target phenotype: {PHENOTYPE}.
Read only assigned {MATERIAL_PATH} and needed fields in {FULL_MATERIAL_PATH}.
Assigned genes: {ASSIGNED_GENES}. Never substitute another batch. Do not use subagents.
Investigate every assigned gene and plausible reference mechanism, retaining copies,
unknowns and counterevidence. Screen about 20 distinct papers per original reference
gene/shared mechanism when possible and carry five useful distinct papers, strongest
applicable evidence first. Keep all distinct observations from those papers.
Save ONE JSON object to {RESULTS_PATH}:
{"contract":"agentic-route-b-science-2.0.0","paper_observations":[],"gene_decisions":[]}
Each paper_observations record needs a unique observation_id, paper_id (PMID/PMCID/DOI),
title copied exactly from the source, reported_gene, reported_gene_id (original
namespace or explicit unknown), reported_species, and finding. Retain source_excerpt,
source_location, reported_process/reported_organ/reported_stage, counterevidence and checked unknown
reasons when available. One paper may contain different genes and observations:
give each a distinct observation_id. Share ONLY by explicit observation_id reference,
never by inferred similarity, paper ID alone, gene symbol or family.
Each gene_decisions record needs interval_gene_id, exact frozen reference_gene,
state (assessed/not_assessed/retrieval_failed), basis
(paper_evidence/source_bound_function/annotation/identity_only/no_bridge/not_assessed), function,
reason, and paper_assessments (empty if none). Function/reason must describe the
original reference mechanism for source_bound_function/annotation; do not copy a
different paper's observations to a database source. A database function can coexist
with papers. For assessed no_bridge explain the judgment; unknown is never no.
Each paper_assessments entry needs observation_id, gene_correspondence
(same_gene/homolog/different_gene/unknown), checklist with ALL nine yes/no/unknown
answers, and a short reason for this target gene/copy. Record copy-specific
paper_checklist_reasons/counterevidence separately on this assessment. Checklist keys:
{CHECKLIST_KEYS}.
Do not fill unknown by guessing. Missing assessment differs from an explicit unknown.
Do not write Tier numbers or per-gene literature administration fields. Machines
derive Tier with the existing checklist/binding rules and check bibliography later.
Paper identity checks do not establish biological support. Never invent sources.
"""


def _affinity_key(unit: dict[str, Any]) -> str:
    """Exact frozen endpoint/species/material snapshot; ordering only, never adoption."""
    snapshots = []
    for gene_id in unit["member_interval_genes"]:
        gene = _unit_gene(unit, gene_id)
        endpoints = _endpoint_species(unit, gene)
        if not endpoints or any(not str(species).strip() for species in endpoints.values()):
            return "unshared:" + unit["unit_id"]
        snapshots.append({"endpoints": endpoints,
                          "context": {k: v for k, v in gene.items() if k in CARD_FIELDS or k.startswith("reference_") or k == "atlas_material_context"}})
    cards = [{k: v for k, v in card.items() if k != "interval_gene_id"}
             for card in unit.get("reference_material_cards", [])]
    return json.dumps({"genes": snapshots, "cards": cards,
                       "sources": unit.get("source_provenance", []), "stable": unit.get("stable_source_ids", []),
                       "context": unit.get("atlas_material_context", {})}, sort_keys=True, ensure_ascii=False)


def light_prepare(bundle_arg: str, bundle_sha256: str, out_dir_arg: str,
                  batch_units: int = 20, units_file: str = "", *,
                  output_contract: str = "legacy-v1", packing: str = "original") -> dict[str, Any]:
    """Pack requested units into bounded light batches with material and roster."""
    if output_contract not in {"legacy-v1", "science-v2"} or packing not in {"original", "reference-affinity"}:
        raise ValueError("unknown light output contract or packing")
    if type(batch_units) is not int or batch_units < 1:
        raise ValueError("light batch size must be a positive integer")
    path, bundle = _load_bundle(bundle_arg, bundle_sha256)
    requested = _requested_units(bundle)
    if units_file:
        listed: list[str] = []
        for line in Path(units_file).read_text(encoding="utf-8").splitlines():
            value = line.strip()
            if value and not value.startswith("#") and value not in listed:
                listed.append(value)
        available = {unit["unit_id"] for unit in requested}
        unknown = [unit_id for unit_id in listed if unit_id not in available]
        if unknown:
            raise ValueError(f"light units file names units that are not requested: {unknown}")
        requested = [unit for unit in requested if unit["unit_id"] in set(listed)]
    if not requested:
        raise ValueError("light prepare requires at least one requested unit")
    if packing == "reference-affinity":
        groups: dict[str, list[dict[str, Any]]] = {}
        for unit in requested:
            groups.setdefault(_affinity_key(unit), []).append(unit)
        requested = [unit for group in groups.values() for unit in group]
    out = Path(out_dir_arg)
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError("light prepare requires a new or empty output directory")
    out.mkdir(parents=True, exist_ok=True)
    batches = []
    for number, start in enumerate(range(0, len(requested), batch_units), start=1):
        bid = BATCH_ID_FORMAT.format(number)
        chunk = requested[start:start + batch_units]
        target = out / "light" / bid
        target.mkdir(parents=True, exist_ok=True)
        material = _render_material(bid, chunk)
        task = (SCIENCE_TASK if output_contract == "science-v2" else TASK).replace("{PHENOTYPE}", str(bundle["phenotype"]).strip())
        task += """
Save original judgments once. The coordinator batches local checks with
light-review-batch; inspect only its compact findings and selected original facts.
Write explicit corrections only when original research supports them; never guess
missing facts. The coordinator runs light-apply-batch on those saved decisions.
Do not spend model turns wrapping hashes, rebuilding full results, reading full
receipts/member rosters, or repeating checks already recorded by apply. Existing
recorded decisions can be transported by Python without another model invocation.
Shared observation edits require all affected gene IDs; one paper is not one
biological observation. Complete sources remain available through selective reads.
Local readiness does not mean scientific completion or verified support.
"""
        task = task.replace("{CHECKLIST_KEYS}", ", ".join(CHECKLIST_FIELDS))
        task = task.replace("{FOCAL_SPECIES}", str(bundle.get("focal_species") or "unknown; do not assume a species"))
        task = task.replace("{MATERIAL_PATH}", str((target / "material.md").resolve()))
        task = task.replace("{FULL_MATERIAL_PATH}", str((target / "material.full.json").resolve()))
        task = task.replace("{RESULTS_PATH}", str((target / "results.jsonl").resolve()))
        task = task.replace("{ASSIGNED_GENES}", ", ".join(g for unit in chunk for g in unit["member_interval_genes"]))
        roster = {gene_id: row for unit in chunk for gene_id, row in _roster(unit).items()}
        _write_text(target / "material.md", material)
        _write_json(target / "material.full.json", {"units": chunk})
        _write_text(target / "task.md", task)
        _write_json(target / "roster.json", roster)
        batches.append({
            "bid": bid,
            "unit_ids": [unit["unit_id"] for unit in chunk],
            "gene_ids": [gene for unit in chunk for gene in unit["member_interval_genes"]],
            "task_sha256": _sha(target / "task.md"),
            "material_sha256": _sha(target / "material.md"),
            "material_full_sha256": _sha(target / "material.full.json"),
            "roster_sha256": _sha(target / "roster.json"),
        })
    index = {
        "contract": V2_CONTRACT if output_contract == "science-v2" else CONTRACT,
        "output_contract": output_contract,
        "packing": packing,
        "case_id": bundle["case_id"],
        "phenotype": bundle["phenotype"],
        "bundle": {"path": str(path), "sha256": bundle_sha256},
        "batch_units": batch_units,
        "requested_unit_count": len(requested),
        "units_file": str(Path(units_file).resolve()) if units_file else "",
        "batches": batches,
        "limitations": "Light batches are delivery only; no gene is judged, reviewed or selected here.",
    }
    _write_json(out / "light" / "index.json", index)
    return {"contract": index["contract"], "case_id": bundle["case_id"], "out": str(out.resolve()),
            "batch_count": len(batches), "requested_unit_count": len(requested),
            "gene_count": sum(len(batch["gene_ids"]) for batch in batches),
            "index": {"file": "light/index.json", "sha256": _sha(out / "light" / "index.json")}}


# ---------------------------------------------------------------------------
# paper metadata: exact identifier resolution and a local cache
# ---------------------------------------------------------------------------

_sleep = time.sleep


def _http_json(url: str, timeout: int = 30) -> dict[str, Any]:
    """GET JSON, retrying 429/5xx and network errors a few times with backoff."""
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    failure: Exception | None = None
    for attempt in range(HTTP_TRIES):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                payload = json.loads(response.read().decode("utf-8"))
            if not isinstance(payload, dict):
                raise RuntimeError("literature provider returned an invalid payload shape")
            return payload
        except urllib.error.HTTPError as error:
            if error.code != 429 and error.code < 500:
                raise
            failure = error
        except (urllib.error.URLError, TimeoutError, ConnectionError) as error:
            failure = error
        if attempt + 1 < HTTP_TRIES:
            _sleep(HTTP_BACKOFF_SECONDS * 2 ** attempt)
    assert failure is not None
    raise failure


def _source_url(pmid: str, pmcid: str, doi: str) -> str:
    if doi:
        return "https://doi.org/" + urllib.parse.quote(doi, safe="/:();")
    if pmcid:
        return f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
    if pmid:
        return f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
    return ""


def _epmc_record(result: dict[str, Any]) -> dict[str, Any]:
    pmid = str(result.get("pmid", "") or "").strip()
    pmcid = str(result.get("pmcid", "") or "").strip()
    doi = str(result.get("doi", "") or "").strip()
    journal = str(result.get("journalTitle", "") or "").strip()
    if not journal and isinstance(result.get("journalInfo"), dict):
        journal = str((result["journalInfo"].get("journal") or {}).get("title", "") or "").strip()
    return {"pmid": pmid, "pmcid": pmcid, "doi": doi,
            "title": _plain(result.get("title")), "abstract": _plain(result.get("abstractText")),
            "journal": journal, "year": str(result.get("pubYear", "") or "").strip(),
            "authors": str(result.get("authorString", "") or "").strip(),
            "url": _source_url(pmid, pmcid, doi), "retrieved_at": _now()}


def _carries(record: dict[str, Any], kind: str, identifier: str) -> bool:
    """The returned record must itself name the id that was requested."""
    value = str(record.get(kind) or "").strip()
    return bool(value) and value.casefold() == identifier.strip().casefold()


def _pmcid_mapping(pmcid: str, timeout: int = 30) -> dict[str, str]:
    url = NCBI_IDCONV + "?" + urllib.parse.urlencode({"ids": pmcid, "format": "json", "tool": "candidate-qtl"})
    payload = _http_json(url, timeout=timeout)
    for record in payload.get("records") or []:
        if isinstance(record, dict) and str(record.get("pmcid", "")).upper() == pmcid.upper():
            return {"pmid": str(record.get("pmid", "") or ""), "doi": str(record.get("doi", "") or "")}
    return {}


def _europe_pmc(kind: str, identifier: str, timeout: int) -> dict[str, Any] | None:
    query = f"EXT_ID:{identifier} AND SRC:MED" if kind == "pmid" else (
        f'DOI:"{identifier}"' if kind == "doi" else f"PMCID:{identifier}")
    url = EUROPE_PMC_REST + "?" + urllib.parse.urlencode(
        {"query": query, "format": "json", "pageSize": 5, "resultType": "core"})
    payload = _http_json(url, timeout=timeout)
    result_list = payload.get("resultList")
    results = result_list.get("result") if isinstance(result_list, dict) else None
    if not isinstance(results, list) or any(not isinstance(result, dict) for result in results):
        raise RuntimeError("Europe PMC returned an invalid result list")
    for result in results:
        if isinstance(result, dict):
            record = _epmc_record(result)
            if _carries(record, kind, identifier):
                record["_source_provider"] = "europe_pmc_core"
                record["_source_url"] = url
                return record
    return None


def _crossref_doi(identifier: str, timeout: int) -> dict[str, Any] | None:
    """Exact DOI metadata only; never infer an abstract from title or snippets."""
    url = CROSSREF_WORKS + urllib.parse.quote(identifier, safe="")
    try:
        payload = _http_json(url, timeout=timeout)
    except urllib.error.HTTPError as error:
        if error.code == 404:
            return None
        raise
    message = payload.get("message")
    if payload.get("status") != "ok" or not isinstance(message, dict):
        raise RuntimeError("Crossref returned an invalid work response")
    if not _carries({"doi": message.get("DOI")}, "doi", identifier):
        return None
    titles = message.get("title")
    if not isinstance(titles, list) or any(not isinstance(title, str) for title in titles):
        raise RuntimeError("Crossref returned an invalid title list")
    abstract = message.get("abstract", "")
    if not isinstance(abstract, str):
        raise RuntimeError("Crossref returned an invalid abstract")
    published = message.get("published") or message.get("issued") or {}
    dates = published.get("date-parts", []) if isinstance(published, dict) else []
    year = str(dates[0][0]) if dates and isinstance(dates[0], list) and dates[0] else ""
    journal = message.get("container-title") or []
    authors = message.get("author") or []
    retrieved_at = _now()
    return {
        "pmid": "", "pmcid": "", "doi": str(message["DOI"]).strip(),
        "title": _plain(titles[0]) if titles else "", "abstract": _plain(abstract),
        "journal": _plain(journal[0]) if isinstance(journal, list) and journal else "",
        "year": year,
        "authors": "; ".join(" ".join(str(author.get(key, "") or "") for key in ("given", "family")).strip()
                             for author in authors if isinstance(author, dict)),
        "url": _source_url("", "", str(message["DOI"]).strip()), "retrieved_at": retrieved_at,
        "_source_provider": "crossref", "_source_url": url,
        "metadata_provenance": [{
            "provider": "crossref", "url": url, "retrieved_at": retrieved_at,
            "requested_doi": identifier,
            "payload_sha256": _sha256(json.dumps(payload, ensure_ascii=False, sort_keys=True).encode("utf-8")),
            "abstract_state": "retrieved" if _plain(abstract) else "not_in_metadata",
        }],
    }


def _fetch_paper(kind: str, identifier: str, timeout: int = 30) -> dict[str, Any] | None:
    """Europe PMC first; missing DOI metadata falls back to exact Crossref lookup."""
    primary_failure = None
    try:
        record = _europe_pmc(kind, identifier, timeout)
    except Exception as error:
        if kind != "doi":
            raise
        primary_failure = error
        record = None
    if record is None and kind == "doi":
        try:
            record = _crossref_doi(identifier, timeout)
        except Exception as error:
            if primary_failure is not None:
                raise RuntimeError(f"Europe PMC: {primary_failure}; Crossref: {error}") from error
            raise
        if record is not None:
            record["metadata_provenance"].insert(0, {
                "provider": "europe_pmc_core", "requested_doi": identifier,
                "url": EUROPE_PMC_REST + "?" + urllib.parse.urlencode({
                    "query": f'DOI:"{identifier}"', "format": "json", "pageSize": 5, "resultType": "core"}),
                "state": "provider_failure" if primary_failure is not None else "not_resolved",
                "detail": str(primary_failure) if primary_failure is not None else "no exact DOI record",
            })
        elif primary_failure is not None:
            raise primary_failure
    if record is None and kind == "pmcid":
        mapping = _pmcid_mapping(identifier, timeout=timeout)
        for alias in ("pmid", "doi"):
            if not mapping.get(alias):
                continue
            candidate = _europe_pmc(alias, mapping[alias], timeout)
            if candidate is None or (candidate["pmcid"] and candidate["pmcid"].upper() != identifier.upper()):
                continue
            candidate["pmcid"] = candidate["pmcid"] or identifier.upper()
            record = candidate
            break
    if record is None:
        return None
    record["aliases"] = _cache_aliases(f"{kind}:{identifier}", record)
    return record


def _load_metadata_cache(path_arg: str) -> dict[str, dict[str, Any]]:
    if not path_arg:
        return {}
    path = Path(path_arg)
    if not path.exists():
        raise ValueError(f"light metadata cache does not exist: {path}")
    payload = _read_json(path)
    if isinstance(payload, dict) and isinstance(payload.get("entries"), dict):
        payload = payload["entries"]
    if not isinstance(payload, dict):
        raise ValueError("light metadata cache must be a JSON object keyed by paper id")
    return {str(key): entry for key, entry in payload.items() if isinstance(entry, dict)}


def _cache_put(cache: dict[str, dict[str, Any]], record: dict[str, Any], requested: object) -> None:
    for key in dict.fromkeys([*_cache_aliases(requested, record), *record.get("aliases", [])]):
        cache.setdefault(key, record)


def _entry_aliases(entry: dict[str, Any]) -> list[str]:
    aliases = entry.get("aliases")
    return [str(alias) for alias in aliases] if isinstance(aliases, list) else []


def _cover_cache(cache: dict[str, dict[str, Any]], out_metadata: dict[str, Any],
                 out_source_verify: list[dict[str, Any]]) -> None:
    """Materialise both cache formats from the same resolved records."""
    records: dict[tuple[str, str, str], dict[str, Any]] = {}
    for key, entry in cache.items():
        if not str(entry.get("title", "")).strip():
            continue
        identity = (str(entry.get("pmid", "")), str(entry.get("pmcid", "")), str(entry.get("doi", "")))
        if not any(identity):
            # A hand-written cache may carry only title/abstract/url/retrieved_at.
            # The key already names the identifier, so resolve it without
            # inventing an alias the supplied entry does not contain.
            kind, _, identifier = str(key).partition(":")
            kind = kind.casefold()
            if kind == "pmid":
                identity = (identifier, "", "")
            elif kind == "pmcid":
                identity = ("", identifier.upper(), "")
            elif kind == "doi":
                identity = ("", "", identifier)
        if any(identity):
            records.setdefault(identity, entry)
        out_metadata[key] = {
            "pmid": identity[0], "pmcid": identity[1], "doi": identity[2],
            "title": str(entry.get("title", "")),
            "abstract": str(entry.get("abstract", "")),
            "url": str(entry.get("url", "")),
            "retrieved_at": str(entry.get("retrieved_at", "")),
            "journal": str(entry.get("journal", "")),
            "year": str(entry.get("year", "")),
            "authors": str(entry.get("authors", "")),
            "aliases": sorted({_cache_key(alias) for alias in _entry_aliases(entry)}),
            "_source_provider": str(entry.get("_source_provider", "local_record_cache")),
            "_source_url": str(entry.get("_source_url", entry.get("url", ""))),
            "metadata_provenance": entry.get("metadata_provenance", []),
            "abstract_provenance": entry.get("abstract_provenance", []),
        }
    from candidate_qtl.source_verify import _cache_key as source_cache_key

    for identity in sorted(records, key=lambda item: (item[0], item[2])):
        entry = records[identity]
        pmid, pmcid, doi = identity
        cache_key = (source_cache_key("pmid", pmid) if pmid else
                     source_cache_key("pmcid", pmcid) if pmcid else source_cache_key("doi", doi))
        out_source_verify.append({
            "cache_key": cache_key,
            "id": pmid or pmcid or doi,
            "pmid": pmid, "pmcid": pmcid, "doi": doi,
            "title": str(entry.get("title", "")),
            "abstractText": str(entry.get("abstract", "")),
            "journalTitle": str(entry.get("journal", "")),
            "pubYear": str(entry.get("year", "")),
            "authorString": str(entry.get("authors", "")),
            "retrieved_at": str(entry.get("retrieved_at", "")),
            "_source_provider": str(entry.get("_source_provider", "local_record_cache")),
            "_source_url": str(entry.get("_source_url", entry.get("url", ""))),
            "metadata_provenance": entry.get("metadata_provenance", []),
            "abstract_provenance": entry.get("abstract_provenance", []),
        })


# ---------------------------------------------------------------------------
# conversion to the canonical artifacts
# ---------------------------------------------------------------------------

def _copy_science_paper_fields(paper: dict[str, Any], assessment: dict[str, Any],
                              annotations: dict[str, Any]) -> None:
    """Transport recorded copy answers; never infer an unknown explanation or an attempted follow-up."""
    from candidate_qtl.observation_tiers import is_unknown

    fields = ("unknown_reason", "checked_unknown_reasons", "followup_state", "paper_checklist_reasons", "checklist_reasons")
    inputs = {f"assessment.{key}": copy.deepcopy(assessment[key]) for key in fields if key in assessment}
    inputs.update({f"observation.{key}": copy.deepcopy(paper[key]) for key in fields if key in paper})
    if "paper_checklist_reasons" in annotations:
        inputs["existing_annotations.paper_checklist_reasons"] = copy.deepcopy(annotations["paper_checklist_reasons"])
    if not inputs:
        return
    provenance: dict[str, Any] = {"input_fields": inputs, "applied_fields": {}, "conflicts": []}

    def explanation(value):
        return isinstance(value, str) and bool(value.strip()) and not is_unknown(value)

    def place(field, value, source):
        old = paper.get(field)
        if field == "paper_checklist_reasons" and field not in paper:
            old = annotations.get(field)
        unfilled = old in (None, "", {}) or field == "unknown_reason" and isinstance(old, str) and (not old.strip() or is_unknown(old))
        if not unfilled and old != value:
            provenance["conflicts"].append({"field": field, "existing_value": copy.deepcopy(old),
                                             "incoming_value": copy.deepcopy(value), "incoming_source": source})
            return
        paper[field] = copy.deepcopy(value)
        provenance["applied_fields"][field] = source

    canonical_reasons = assessment.get("paper_checklist_reasons")
    canonical_present = canonical_reasons not in (None, "", {})
    reasons = canonical_reasons if canonical_present else assessment.get("checklist_reasons")
    if isinstance(reasons, dict):
        place("paper_checklist_reasons", reasons,
              "assessment.paper_checklist_reasons" if canonical_present else "assessment.checklist_reasons")
    for field in ("unknown_reason", "followup_state"):
        value = assessment.get(field)
        valid = explanation(value) if field == "unknown_reason" else value in ("not_attempted", "completed")
        if valid:
            place(field, value, "assessment." + field)
        elif field in assessment and value not in (None, ""):
            provenance.setdefault("untransported_fields", []).append("assessment." + field)
    if not explanation(paper.get("unknown_reason")):
        aliases = {}
        for location, source in (("assessment", assessment), ("observation", paper)):
            if explanation(source.get("checked_unknown_reasons")):
                aliases[location + ".checked_unknown_reasons"] = source["checked_unknown_reasons"]
        for name in ("paper_checklist_reasons", "checklist_reasons"):
            source = assessment.get(name)
            if isinstance(source, dict) and explanation(source.get("checked_unknown_reasons")):
                aliases["assessment." + name + ".checked_unknown_reasons"] = source["checked_unknown_reasons"]
        if aliases:
            distinct = list(dict.fromkeys(aliases.values()))
            value = distinct[0] if len(distinct) == 1 else json.dumps(aliases, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
            place("unknown_reason", value, {"literal_checked_unknown_aliases": aliases})
        elif isinstance(reasons, dict):
            unknown = [key for key, answer in assessment["checklist"].items() if answer == "unknown"]
            # Do not let a complete checklist map conceal an unexplained
            # unknown reported identifier/species/correspondence.
            unknown.extend(key for key in ("title", "reported_gene", "reported_gene_id", "reported_species", "finding")
                           if is_unknown(paper.get(key)))
            if is_unknown(assessment.get("gene_correspondence")):
                unknown.append("gene_correspondence")
            if unknown and all(explanation(reasons.get(key)) for key in unknown):
                original = {key: reasons[key] for key in unknown}
                place("unknown_reason", json.dumps(original, ensure_ascii=False, sort_keys=True, separators=(",", ":")),
                      {"literal_unknown_answer_reasons": original})
            elif unknown:
                provenance["unexplained_unknown_fields"] = [key for key in unknown if not explanation(reasons.get(key))]
    provenance["state"] = "conflicting_original_values" if provenance["conflicts"] else "transported_original_fields"
    if "copy_science_field_transport" in annotations:
        provenance["previous_provenance"] = copy.deepcopy(annotations["copy_science_field_transport"])
    annotations["copy_science_field_transport"] = provenance


def _science_rows(path: Path) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]]]:
    if not path.exists():
        return {}, []
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or payload.get("contract") != SCIENCE_CONTRACT:
        raise ValueError("science-v2 results require the explicit science contract")
    observations = payload.get("paper_observations")
    decisions = payload.get("gene_decisions")
    if not isinstance(observations, list) or not isinstance(decisions, list):
        raise ValueError("science-v2 requires paper_observations and gene_decisions arrays")
    catalog: dict[str, dict[str, Any]] = {}
    conflicts = set()
    for observation in observations:
        if not isinstance(observation, dict) or not isinstance(observation.get("observation_id"), str):
            raise ValueError("science observation requires observation_id")
        key = observation["observation_id"]
        if not key.strip() or key in catalog:
            conflicts.add(key)
        catalog[key] = observation
    rows, issues = {}, []
    basis_tiers = {"paper_evidence": "U", "source_bound_function": "T4", "annotation": "T5", "identity_only": "T6", "no_bridge": "U"}
    for decision in decisions:
        if not isinstance(decision, dict):
            issues.append({"gene_id": "", "issue": "invalid_science_decision"})
            continue
        gene_id = str(decision.get("interval_gene_id", "")).strip()
        row = {**decision, "papers": [], "_science_state": decision.get("state"), "_science_basis": decision.get("basis")}
        try:
            if not gene_id or gene_id in rows:
                raise ValueError("missing or duplicate interval_gene_id")
            state, basis = decision.get("state"), decision.get("basis")
            if state not in {"assessed", "not_assessed", "retrieval_failed"}:
                raise ValueError("unknown science state")
            if state != "assessed":
                row["tier"] = "not_assessed"
            else:
                if basis not in basis_tiers:
                    raise ValueError("assessed decision requires explicit basis")
                assessments = decision.get("paper_assessments")
                if not isinstance(assessments, list):
                    raise ValueError("paper_assessments must be an array")
                if basis == "paper_evidence" and not assessments:
                    raise ValueError("paper_evidence requires a paper assessment")
                for assessment in assessments:
                    if not isinstance(assessment, dict):
                        raise ValueError("paper assessment must be an object")
                    key = assessment.get("observation_id")
                    if key not in catalog or key in conflicts:
                        raise ValueError("unknown or conflicting observation_id")
                    observation = catalog[key]
                    for field in ("paper_id", "title", "reported_gene", "reported_gene_id", "reported_species", "finding"):
                        if not isinstance(observation.get(field), str) or not observation[field].strip():
                            raise ValueError(f"observation requires {field}")
                    checklist = assessment.get("checklist")
                    if (not isinstance(checklist, dict) or set(checklist) != set(CHECKLIST_FIELDS)
                            or any(answer not in {"yes", "no", "unknown"} for answer in checklist.values())):
                        raise ValueError("all nine original checklist answers are required")
                    if assessment.get("gene_correspondence") not in {"same_gene", "homolog", "different_gene", "unknown"}:
                        raise ValueError("paper assessment requires correspondence")
                    if not isinstance(assessment.get("reason"), str) or not assessment["reason"].strip():
                        raise ValueError("paper assessment requires a copy-specific reason")
                    paper = {k: v for k, v in observation.items() if k != "observation_id"}
                    annotations = dict(paper.get("research_annotations", {}))
                    annotations["observation_id"] = key
                    # Unrecognized original context is retained as metadata, not
                    # silently dropped or interpreted as a scientific answer.
                    from candidate_qtl.paper_reports import LEGACY_RESEARCH_FIELDS
                    carried = {"paper_id", "title", "reported_gene", "reported_gene_id", "reported_species", "finding",
                               "reported_phenotype", "reported_process", "phenotype_comparison", "reported_organ", "reported_stage",
                               "source_excerpt", "source_location", "unknown_reason", "publication_type", "underlying_source_ids",
                               "corroboration_note", "followup_question", "followup_state", "paper_checklist_missing_fields_json",
                               "reported_function", "reported_annotation", "research_annotations", *LEGACY_RESEARCH_FIELDS}
                    extra = {k: v for k, v in paper.items() if k not in carried}
                    if extra:
                        annotations["original_observation_fields"] = extra
                    annotations["copy_assessment"] = {k: v for k, v in assessment.items() if k not in {"observation_id", "checklist", "gene_correspondence"}}
                    _copy_science_paper_fields(paper, assessment, annotations)
                    paper.update(gene_correspondence=assessment["gene_correspondence"], checklist=checklist,
                                 research_annotations=annotations)
                    row["papers"].append(paper)
                tiers = [_checklist_tier(paper["checklist"]) for paper in row["papers"]]
                # This is a transport declaration only; applicable-paper binding
                # and final strongest evidence remain the existing consumer's job.
                row["tier"] = min(tiers + ([basis_tiers[basis]] if basis in {"source_bound_function", "annotation"} else []),
                                  key=STATE_ORDER.index) if tiers else basis_tiers[basis]
        except (ValueError, TypeError) as error:
            issues.append({"gene_id": gene_id, "issue": "invalid_science_decision", "detail": str(error)})
            row.update(tier="not_assessed", _science_state="invalid", _conflict=True)
        if gene_id:
            rows[gene_id] = row
    return rows, issues


def _results_rows(path: Path) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]]]:
    """One row per gene. Identical duplicates collapse; conflicting ones leave the gene not_assessed."""
    issues: list[dict[str, Any]] = []
    rows: dict[str, dict[str, Any]] = {}
    if not path.exists():
        return rows, issues
    bibliography: dict[str, Any] = {}
    catalog_seen = False
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            row = None
        if not isinstance(row, dict):
            issues.append({"gene_id": "", "issue": "invalid_json_row", "line": number})
            continue
        if set(row) == {"bibliography"}:
            from candidate_qtl.paper_transport import expand_paper_observations
            if catalog_seen or rows:
                raise ValueError("bibliography must be a single header before gene rows")
            expand_paper_observations({"bibliography": row["bibliography"], "observations": []})
            bibliography, catalog_seen = row["bibliography"], True
            continue
        gene_id = str(row.get("interval_gene_id", "") or "").strip()
        if not gene_id:
            issues.append({"gene_id": "", "issue": "row_without_interval_gene_id", "line": number})
            continue
        raw_papers = row.get("papers", [])
        if bibliography or isinstance(raw_papers, dict) or (isinstance(raw_papers, list) and any(
                isinstance(paper, dict) and "paper_ref" in paper for paper in raw_papers)):
            from candidate_qtl.paper_transport import expand_paper_observations
            try:
                row["papers"] = expand_paper_observations(row.get("papers", []), bibliography,
                                                        id_field="paper_id", title_field="title")
            except ValueError as error:
                issues.append({"gene_id": gene_id, "issue": "invalid_paper_reference", "line": number, "detail": str(error)})
                row = {"interval_gene_id": gene_id, "tier": "not_assessed", "_conflict": True}
        if gene_id not in rows:
            rows[gene_id] = row
            continue
        if rows[gene_id].get("_conflict"):
            continue
        if json.dumps(row, sort_keys=True) == json.dumps(rows[gene_id], sort_keys=True):
            issues.append({"gene_id": gene_id, "issue": "duplicate_identical_row_collapsed", "line": number})
        else:
            issues.append({"gene_id": gene_id, "issue": "duplicate_conflicting_rows_not_assessed", "line": number})
            rows[gene_id] = {"interval_gene_id": gene_id, "tier": "not_assessed", "_conflict": True}
    return rows, issues


def _paper_state(cache: dict[str, dict[str, Any]], paper_id: object, title: object,
                 offline: bool, fetch: Any) -> tuple[dict[str, Any] | None, str, str]:
    """Return (resolved record, state, detail) for one declared paper."""
    kind, identifier = _citation(paper_id)
    if kind not in _LABELS:
        return None, "not_resolved", "paper_id is not a PMID, PMCID or DOI"
    key = _cache_key(f"{kind}:{identifier}")
    record = cache.get(key)
    if record is None:
        if offline:
            return None, "not_attempted", "offline single-cache run has no entry for this paper"
        try:
            record = fetch(kind, identifier)
        except Exception as error:  # a provider failure is not absence
            return None, "provider_failure", f"{type(error).__name__}: {error}"
        if record is None:
            return None, "not_resolved", "literature providers returned no exact identifier record"
        if not _carries(record, kind, identifier):
            return None, "not_resolved", "the returned record does not carry the requested identifier"
        record = {**record, "aliases": _cache_aliases(f"{kind}:{identifier}", record)}
        _cache_put(cache, record, f"{kind}:{identifier}")
    if record.get(kind) and not _carries(record, kind, identifier):
        return None, "not_resolved", "the cached record conflicts with the requested identifier"
    resolved_title = _plain(record.get("title"))
    if not resolved_title:
        return None, "not_resolved", "resolved record has no title"
    ratio = _title_similarity(title, resolved_title)
    if ratio < TITLE_THRESHOLD:
        return None, "title_mismatch", f"title similarity {ratio:.2f} below {TITLE_THRESHOLD}"
    return {**record, "title": resolved_title, "abstract": _plain(record.get("abstract"))}, "checked", \
        f"title similarity {ratio:.2f}"


def _strict_checklist(value: object) -> dict[str, str] | None:
    """The researcher's checklist, unchanged, or None when it is not all nine yes/no/unknown answers."""
    if (not isinstance(value, dict) or set(value) != set(CHECKLIST_FIELDS)
            or any(answer not in ("yes", "no", "unknown") for answer in value.values())):
        return None
    return dict(value)


def _checklist_tier(checks: dict[str, str]) -> str:
    from candidate_qtl.observation_tiers import checklist_tier

    return checklist_tier(checks)


def _check_paper(paper: object, cache: dict[str, dict[str, Any]], offline: bool, fetch: Any,
                 stable_sources: list[str], reference_species: str) -> tuple[dict[str, Any] | None, str, str, str]:
    """Hallucination guards plus 'can a discovery be formed'. Returns (kept, state, detail, paper_id)."""
    if not isinstance(paper, dict):
        return None, "invalid_record", "paper entry is not an object", ""
    paper_id = str(paper.get("paper_id", "") or "").strip()
    record, state, detail = _paper_state(cache, paper_id, paper.get("title", ""), offline, fetch)
    if record is None:
        return None, state, detail, paper_id
    checklist = _strict_checklist(paper.get("checklist"))
    correspondence = str(paper.get("gene_correspondence", "") or "").strip()
    reported_gene = str(paper.get("reported_gene", "") or "").strip()
    finding = str(paper.get("finding", "") or "").strip()
    if (checklist is None or correspondence not in {"same_gene", "homolog", "different_gene", "unknown"}
            or not reported_gene or not finding):
        return None, "paper_report_incomplete", (
            "reported_gene, gene_correspondence, finding and a checklist with all nine "
            "yes/no/unknown answers are required"), paper_id
    kind, identifier = _citation(paper_id)
    key = _cache_key(f"{kind}:{identifier}")
    alias_keys = {key, *_cache_aliases(paper_id, record), *(_cache_key(alias) for alias in _entry_aliases(record))}
    frozen = next((entry for entry in stable_sources if _cache_key(entry) in alias_keys), "")
    source_id = frozen or _canonical_id(kind, identifier, record)
    if not frozen:
        if ";" in source_id or "," in source_id:
            return None, "identifier_not_representable", "the identifier contains a separator the bundle splits on", paper_id
        reported_species = _species(paper.get("reported_species"))
        url = str(record.get("url", "") or "").strip()
        if not str(record.get("abstract", "")).strip():
            return None, "missing_abstract", "resolved record has no abstract for a new discovery", paper_id
        if not url.startswith("https://") or not str(record.get("retrieved_at", "") or "").strip():
            return None, "cache_entry_invalid", "a new source needs an https url and retrieved_at", paper_id
        if not reference_species or not reported_species:
            return None, "discovery_not_formable", (
                "a new source needs the frozen reference species and the reported species"), paper_id
    tokens = [token for token in _record_tokens(record)
              if token.casefold() != source_id.casefold() and ";" not in token and "," not in token]
    from candidate_qtl.paper_reports import LEGACY_RESEARCH_FIELDS
    observation_fields = {field: paper[field] for field in (
        "reported_phenotype", "reported_process", "phenotype_comparison", "reported_organ", "reported_stage",
        "source_excerpt", "source_location", "unknown_reason", "publication_type", "underlying_source_ids",
        "corroboration_note", "followup_question", "followup_state", "paper_checklist_missing_fields_json",
        "reported_function", "reported_annotation", "research_annotations", *sorted(LEGACY_RESEARCH_FIELDS)) if field in paper}
    if "checklist_reasons" in paper and isinstance(observation_fields.get("research_annotations", {}), dict):
        annotations = dict(observation_fields.get("research_annotations", {}))
        provenance = {"input_fields": {field: paper[field] for field in (
            "checklist_reasons", "paper_checklist_reasons") if field in paper}}
        if "checklist_reason_field_provenance" in annotations:
            provenance["previous_provenance"] = annotations["checklist_reason_field_provenance"]
        annotations["checklist_reason_field_provenance"] = provenance
        # Existing canonical placement retains precedence; both originals remain auditable.
        if "paper_checklist_reasons" not in observation_fields and "paper_checklist_reasons" not in annotations:
            observation_fields["paper_checklist_reasons"] = paper["checklist_reasons"]
        observation_fields["research_annotations"] = annotations
    location = observation_fields.get("source_location")
    if isinstance(location, (dict, list)) and isinstance(observation_fields.get("research_annotations", {}), dict):
        annotations = dict(observation_fields.get("research_annotations", {}))
        provenance = {"input_value": location}
        if "source_location_field_provenance" in annotations:
            provenance["previous_provenance"] = annotations["source_location_field_provenance"]
        annotations["source_location_field_provenance"] = provenance
        observation_fields["research_annotations"] = annotations
        try:
            observation_fields["source_location"] = (json.dumps(location, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False)
                                                       if location else "")
        except (ValueError, TypeError):
            return None, "invalid_record", "source_location contains values that cannot be serialized as JSON", paper_id
    return ({
        "source_id": source_id, "frozen": bool(frozen), "title": str(paper["title"]).strip(),
        "reported_gene": reported_gene, "reported_species": _species(paper.get("reported_species")),
        "reported_gene_id": str(paper.get("reported_gene_id", "") or "").strip(),
        "observation_fields": observation_fields,
        "gene_correspondence": correspondence, "finding": finding, "checklist": checklist,
        "record": record, "alias_keys": alias_keys, "alias_tokens": tokens,
    }, "checked", detail, paper_id)


def _paper_report(paper: dict[str, Any], reference_gene: str) -> dict[str, Any]:
    return {
        "source_id": paper["source_id"],
        "source_title": paper["title"],
        "reported_gene": paper["reported_gene"],
        "reported_gene_id": paper.get("reported_gene_id", ""),
        "reported_species": paper["reported_species"],
        "reference_gene_id": reference_gene,
        "gene_correspondence": paper["gene_correspondence"],
        "initial_tier_reason": paper["finding"],
        "paper_checklist": paper["checklist"],
        **paper.get("observation_fields", {}),
    }


def _frozen_database_source(unit: dict[str, Any], reference_gene: str) -> str:
    """A database accession frozen in this unit with provenance for exactly this reference gene."""
    from candidate_qtl.source_verify import classify_citation

    stable = {str(value).strip() for value in unit.get("stable_source_ids", [])}
    for record in unit.get("source_provenance", []) or []:
        if not isinstance(record, dict) or str(record.get("reference_gene", "")).strip() != reference_gene:
            continue
        if not any(str(record.get(field, "")).strip() for field in ("source_database", "provenance_detail", "source_url")):
            continue
        for source_id in record.get("source_ids", []) or []:
            source_id = str(source_id).strip()
            if source_id in stable and classify_citation(source_id)[0] == "database_accession":
                return source_id
    return ""


def _claim(gene_id: str, reference_gene: str, reference_species: str, declared: str, computed: str,
           function: str, reason: str, papers: list[dict[str, Any]], database_source: str) -> dict[str, Any]:
    claim: dict[str, Any] = {
        # A database-accession source is verified only for a curated membership claim type.
        "claim_type": "literature_bridge" if papers else "direct_functional_bridge",
        "member_interval_genes": [gene_id],
        "reference_gene": reference_gene,
        "reference_species": reference_species,
        "evidence_scope": "reference_gene",
        "biological_function": function,
        "phenotype_relation": reason,
        "mechanism": function,
        "uncertainty_note": "light Route B: focal-copy causality and cross-species transfer are not established",
        "initial_tier_reason": reason,
    }
    if declared == computed:
        claim["initial_tier"] = declared
    if papers:
        claim["source_ids"] = list(dict.fromkeys(paper["source_id"] for paper in papers))
        claim["paper_reports_json"] = [_paper_report(paper, reference_gene) for paper in papers]
        claim["source_title"] = papers[0]["title"]
    else:
        # No paper_reports_json at all: an empty list makes the importer blank the observation fields.
        claim["source_ids"] = [database_source]
        if declared == "T4":
            claim.update(evidence_types=["function"], reported_function=function)
        else:
            claim.update(evidence_types=["annotation"], reported_annotation=function)
    return claim


def _discovery(unit: dict[str, Any], reference_gene: str, reference_species: str,
               paper: dict[str, Any]) -> dict[str, Any]:
    record = paper["record"]
    return {
        "unit_id": unit["unit_id"],
        "reference_gene": reference_gene,
        "reference_species": reference_species,
        "reported_gene": paper["reported_gene"],
        "reported_gene_id": paper.get("reported_gene_id", ""),
        "reported_species": paper["reported_species"],
        "source_id": paper["source_id"],
        "source_title": paper["title"],
        "source_url": str(record["url"]).strip(),
        "retrieved_at": str(record["retrieved_at"]).strip(),
        "source_text": record["abstract"],
        "reference_function": paper["finding"],
        "source_identifier_aliases": ";".join(paper["alias_tokens"]),
    }


def _validated_index(bundle: dict[str, Any], bundle_sha256: str, index_path: Path,
                     light_index_sha256: str) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    if not light_index_sha256 or _sha(index_path) != light_index_sha256:
        raise ValueError("light index SHA-256 mismatch")
    index = _read_json(index_path)
    if not isinstance(index, dict) or index.get("contract") not in {CONTRACT, V2_CONTRACT}:
        raise ValueError("light index contract is unsupported")
    if index.get("case_id") != bundle["case_id"] or (index.get("bundle") or {}).get("sha256") != bundle_sha256:
        raise ValueError("light index does not bind this bundle")
    units = {unit["unit_id"]: unit for unit in bundle["units"]}
    batches = index.get("batches")
    if not isinstance(batches, list) or not batches:
        raise ValueError("light index requires batches")
    seen_units: set[str] = set()
    for batch in batches:
        if (not isinstance(batch, dict) or not isinstance(batch.get("bid"), str)
                or not re.fullmatch(r"b[0-9]{2,}", batch["bid"])):
            raise ValueError("light index contains an invalid batch")
        unit_ids = batch.get("unit_ids")
        if (not isinstance(unit_ids, list) or not unit_ids or any(uid not in units for uid in unit_ids)
                or len(set(unit_ids)) != len(unit_ids) or seen_units.intersection(unit_ids)):
            raise ValueError("light index batch unit roster does not match the bundle")
        seen_units.update(unit_ids)
        expected_genes = [gene for uid in unit_ids for gene in units[uid]["member_interval_genes"]]
        if set(batch.get("gene_ids") or []) != set(expected_genes):
            raise ValueError("light index batch gene roster does not match the bundle")
        directory = index_path.parent / batch["bid"]
        for name, key in (("material.md", "material_sha256"), ("task.md", "task_sha256"),
                          ("roster.json", "roster_sha256"), ("material.full.json", "material_full_sha256")):
            if key == "material_full_sha256" and key not in batch:
                continue  # Old frozen light batches have no full-material sidecar.
            target = directory / name
            if not target.is_file():
                raise ValueError(f"light batch file is missing: {target}")
            if str(batch.get(key, "")) != _sha(target):
                raise ValueError(f"light batch file hash mismatch: {target}")
    if not seen_units <= {unit["unit_id"] for unit in _requested_units(bundle)}:
        raise ValueError("light index names units that are not requested")
    return index, batches


def light_convert(bundle_arg: str, bundle_sha256: str, light_index_arg: str,
                  light_index_sha256: str, out_dir_arg: str, *,
                  cache_arg: str = "", offline: bool = False, fetch: Any = None) -> dict[str, Any]:
    """Check one light run and convert it into discoveries and an agent file."""
    fetch = fetch or _fetch_paper
    bundle_path, bundle = _load_bundle(bundle_arg, bundle_sha256)
    index_path = Path(light_index_arg).resolve()
    index, declared = _validated_index(bundle, bundle_sha256, index_path, light_index_sha256)
    science_v2 = index["contract"] == V2_CONTRACT
    units = {unit["unit_id"]: unit for unit in bundle["units"]}
    requested = _requested_units(bundle)
    seen_units = {unit_id for batch in declared for unit_id in batch["unit_ids"]}

    out = Path(out_dir_arg)
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError("light convert requires a new or empty output directory")
    out.mkdir(parents=True, exist_ok=True)
    cache = _load_metadata_cache(cache_arg)
    checks: dict[str, Any] = {"batches": {}, "gene_issues": [], "paper_issues": [],
                              "not_evaluated_units": [], "dropped_claims": [], "tier_comparison": []}
    totals = {"genes_seen": 0, "genes_requested": 0, "genes_not_assessed": 0, "genes_missing": 0,
              "invalid_gene_id": 0, "invalid_reference": 0, "invalid_tier": 0, "claims": 0,
              "claims_without_paper": 0, "tier_differs_from_declared": 0,
              "papers_seen": 0, "papers_checked": 0, "papers_not_attempted": 0,
              "papers_provider_failure": 0, "papers_not_resolved": 0, "papers_title_mismatch": 0,
              "papers_missing_abstract": 0, "papers_incomplete": 0, "papers_discovery_not_formable": 0,
              "papers_duplicate": 0, "discoveries": 0}
    state_counter = {"not_attempted": "papers_not_attempted", "provider_failure": "papers_provider_failure",
                     "not_resolved": "papers_not_resolved", "title_mismatch": "papers_title_mismatch",
                     "missing_abstract": "papers_missing_abstract", "invalid_record": "papers_incomplete",
                     "paper_report_incomplete": "papers_incomplete",
                     "cache_entry_invalid": "papers_discovery_not_formable",
                     "discovery_not_formable": "papers_discovery_not_formable",
                     "identifier_not_representable": "papers_discovery_not_formable"}
    claims_by_unit: dict[str, list[dict[str, Any]]] = {}
    judged_by_unit: dict[str, set[str]] = {}
    dropped_by_unit: dict[str, list[str]] = {}
    discoveries: list[dict[str, Any]] = []
    discovery_keys: set[tuple[str, str, str]] = set()
    member_coverage: dict[str, dict[str, Any]] = {
        gene: {"unit_id": unit["unit_id"], "gene_id": gene, "state": "not_assessed", "reason": "outside_light_run"}
        for unit in bundle["units"] for gene in unit["member_interval_genes"]}

    for batch in declared:
        bid = batch["bid"]
        results_path = index_path.parent / bid / "results.jsonl"
        batch_checks: dict[str, Any] = {
            "rows": 0, "missing_genes": [], "not_assessed_genes": [],
            "results_sha256": _sha(results_path) if results_path.is_file() else None}
        by_gene, row_issues = (_science_rows if science_v2 else _results_rows)(results_path)
        checks["gene_issues"].extend({"batch": bid, **issue} for issue in row_issues)
        batch_checks["rows"] = len(by_gene)
        batch_genes = {gene for unit_id in batch["unit_ids"] for gene in units[unit_id]["member_interval_genes"]}
        for gene_id in by_gene:
            if gene_id not in batch_genes:
                totals["invalid_gene_id"] += 1
                checks["gene_issues"].append({"batch": bid, "unit_id": "", "gene_id": gene_id,
                                              "issue": "gene_not_in_batch"})
        for unit_id in batch["unit_ids"]:
            unit = units[unit_id]
            roster = _roster(unit)
            judged_by_unit.setdefault(unit_id, set())
            claims_by_unit.setdefault(unit_id, [])
            stable_sources = [str(value).strip() for value in unit.get("stable_source_ids", [])
                              if str(value).strip()]
            for gene_id in unit["member_interval_genes"]:
                member_coverage[gene_id].update(state="not_assessed", reason="missing_gene_response", batch=bid)
                totals["genes_requested"] += 1
                row = by_gene.get(gene_id)
                if row is None:
                    totals["genes_missing"] += 1
                    batch_checks["missing_genes"].append(gene_id)
                    continue
                where = {"batch": bid, "unit_id": unit_id, "gene_id": gene_id}
                totals["genes_seen"] += 1
                member_coverage[gene_id].update(state="invalid", reason="machine_validation_pending")
                tier = str(row.get("tier", "") or "").strip()
                member_coverage[gene_id].update(raw_basis=row.get("_science_basis", "legacy_declared_tier"),
                                                derived_transport_tier=tier)
                if tier not in TIERS:
                    totals["invalid_tier"] += 1
                    checks["gene_issues"].append({**where, "issue": "invalid_tier", "value": tier})
                    continue
                if tier == "not_assessed":
                    state = row.get("_science_state", "invalid" if row.get("_conflict") else "not_assessed")
                    member_coverage[gene_id].update(state=state if state in {"not_assessed", "retrieval_failed", "invalid"} else "not_assessed",
                                                    reason=str(row.get("reason", "") or "not_assessed"))
                    totals["genes_not_assessed"] += 1
                    batch_checks["not_assessed_genes"].append(gene_id)
                    continue
                judged_by_unit[unit_id].add(gene_id)
                if tier in {"T6", "U"} and not (science_v2 and row.get("papers")):
                    if science_v2 and str(row.get("reference_gene", "")).casefold() not in {
                            endpoint.casefold() for endpoint in roster[gene_id]["reference_endpoints"]}:
                        member_coverage[gene_id].update(reason="reference_not_in_frozen_endpoints")
                        dropped_by_unit.setdefault(unit_id, []).append(gene_id)
                        continue
                    member_coverage[gene_id].update(state="assessed", reason="identity_only" if tier == "T6" else "no_bridge")
                    continue
                reference = str(row.get("reference_gene", "") or "").strip()
                frozen = next((value for value in roster[gene_id]["reference_endpoints"]
                               if reference and value.casefold() == reference.casefold()), "")
                if not frozen:
                    member_coverage[gene_id].update(reason="reference_not_in_frozen_endpoints")
                    totals["invalid_reference"] += 1
                    dropped_by_unit.setdefault(unit_id, []).append(gene_id)
                    checks["gene_issues"].append({
                        **where, "issue": ("reference_not_in_frozen_endpoints" if reference
                                           else "judged_gene_without_reference_endpoint"),
                        **({"value": reference} if reference else {})})
                    continue
                species = next((value for key, value in _endpoint_species(unit, _unit_gene(unit, gene_id)).items()
                                if key.casefold() == frozen.casefold()), "")
                declared_papers = row.get("papers")
                if declared_papers is None:
                    declared_papers = []
                if not isinstance(declared_papers, list):
                    checks["gene_issues"].append({**where, "issue": "papers_not_a_list"})
                    declared_papers = []
                kept: list[dict[str, Any]] = []
                dropped_papers: list[dict[str, str]] = []
                for paper in declared_papers:
                    totals["papers_seen"] += 1
                    accepted, state, detail, paper_id = _check_paper(
                        paper, cache, offline, fetch, stable_sources, species)
                    if accepted is not None and any(_paper_report(accepted, frozen) == _paper_report(other, frozen) for other in kept):
                        totals["papers_duplicate"] += 1
                        accepted, state, detail = None, "duplicate_paper", "the exact same observation is already kept for this gene"
                    if accepted is None:
                        if state != "duplicate_paper":
                            totals[state_counter.get(state, "papers_not_resolved")] += 1
                        checks["paper_issues"].append({**where, "paper_id": paper_id, "state": state,
                                                       "detail": detail})
                        dropped_papers.append({"paper_id": paper_id, "state": state})
                        continue
                    totals["papers_checked"] += 1
                    kept.append(accepted)
                if len({paper["source_id"] for paper in kept}) > MAX_PAPERS:
                    checks["gene_issues"].append({**where, "issue": "more_than_five_papers_preserved_for_canonical_selection"})
                if not kept and tier in {"T6", "U"}:
                    independent_basis = science_v2 and row.get("_science_basis") in {"identity_only", "no_bridge"}
                    member_coverage[gene_id].update(
                        state="assessed" if independent_basis or not dropped_papers else _paper_failure_state(dropped_papers),
                        reason=row["_science_basis"] if independent_basis else ("no_checked_paper_remained" if dropped_papers else "no_bridge"))
                    if independent_basis and dropped_papers:
                        # These papers were not the basis of the original T6/U
                        # judgment. Preserve retrieval diagnostics without
                        # turning that independent judgment into missing work.
                        member_coverage[gene_id]["ancillary_paper_issues"] = dropped_papers
                    elif dropped_papers:
                        dropped_by_unit.setdefault(unit_id, []).append(gene_id)
                    continue
                if tier in PAPER_TIERS and not kept:
                    member_coverage[gene_id].update(state=_paper_failure_state(dropped_papers),
                                                    reason="no_checked_paper_remained")
                    if science_v2 and row.get("_science_basis") in {"source_bound_function", "annotation"}:
                        tier = "T4" if row["_science_basis"] == "source_bound_function" else "T5"
                if tier in PAPER_TIERS and not kept:
                    dropped_by_unit.setdefault(unit_id, []).append(gene_id)
                    checks["dropped_claims"].append({
                        **where, "tier": tier, "dropped_papers": dropped_papers,
                        "reason": "tier_requires_a_checked_paper_but_none_remained"})
                    continue
                function = str(row.get("function", "") or "").strip()
                reason = str(row.get("reason", "") or "").strip()
                if not function or not reason:
                    member_coverage[gene_id].update(reason="function_or_reason_missing")
                    dropped_by_unit.setdefault(unit_id, []).append(gene_id)
                    checks["dropped_claims"].append({**where, "tier": tier, "dropped_papers": dropped_papers,
                                                     "reason": "function_or_reason_missing"})
                    continue
                database_source = ""
                if kept:
                    computed = min((_checklist_tier(paper["checklist"]) for paper in kept), key=STATE_ORDER.index)
                else:
                    computed = NO_PAPER_DOWNSTREAM_TIER[tier]
                    database_source = _frozen_database_source(unit, frozen)
                    if not database_source:
                        member_coverage[gene_id].update(reason="no_frozen_database_source_for_this_reference")
                        dropped_by_unit.setdefault(unit_id, []).append(gene_id)
                        checks["dropped_claims"].append({
                            **where, "tier": tier, "dropped_papers": dropped_papers,
                            "reason": "no_paper_and_no_frozen_database_source_for_this_reference"})
                        continue
                if computed != tier:
                    totals["tier_differs_from_declared"] += 1
                    checks["tier_comparison"].append({**where, "declared_tier": tier, "computed_tier": computed,
                                                      "dropped_papers": dropped_papers})
                claims_by_unit[unit_id].append(_claim(gene_id, frozen, species, tier, computed, function, reason,
                                                      kept, database_source))
                totals["claims"] += 1
                totals["claims_without_paper"] += int(not kept)
                member_coverage[gene_id].update(state=_paper_failure_state(dropped_papers) if dropped_papers else "assessed",
                                                reason="paper_check_incomplete" if dropped_papers else "original_judgment_transported")
                # The original explicit function judgment is independent of the
                # paper observations. Keep its own frozen DB source, rather
                # than copying an inapplicable paper's facts to that source.
                db_basis = row.get("_science_basis") if science_v2 else ("source_bound_function" if tier == "T4" else None)
                if kept and db_basis in {"source_bound_function", "annotation"}:
                    database_source = _frozen_database_source(unit, frozen)
                    if database_source:
                        db_tier = "T4" if db_basis == "source_bound_function" else "T5"
                        claims_by_unit[unit_id].append(_claim(
                            gene_id, frozen, species, db_tier, db_tier, function, reason, [], database_source))
                        totals["claims"] += 1
                        totals["claims_without_paper"] += 1
                for paper in kept:
                    if paper["frozen"]:
                        continue
                    key = (unit_id, frozen.casefold(), json.dumps(_discovery(unit, frozen, species, paper), sort_keys=True))
                    if key not in discovery_keys:
                        discovery_keys.add(key)
                        discoveries.append(_discovery(unit, frozen, species, paper))
                        totals["discoveries"] += 1
        checks["batches"][bid] = batch_checks

    responses = []
    for unit in requested:
        unit_id = unit["unit_id"]
        claims = claims_by_unit.get(unit_id, [])
        judged = judged_by_unit.get(unit_id, set())
        members = set(unit["member_interval_genes"])
        record: dict[str, Any] = {"prompt_id": f"route-b-semantic-{unit_id}", "route_b_unit_id": unit_id}
        if unit_id not in seen_units:
            responses.append({**record, "status": "not_evaluated", "reason": "outside light run"})
            checks["not_evaluated_units"].append({"unit_id": unit_id, "reason": "outside light run"})
            continue
        if claims:
            responses.append({**record, "raw_text": json.dumps({"claims": claims}, ensure_ascii=False),
                              "status": "completed"})
            continue
        dropped = sorted(set(dropped_by_unit.get(unit_id, [])))
        missing = sorted(members - judged)
        reason = ""
        if dropped:
            reason = "machine_check_dropped_claim:" + ";".join(dropped)
        elif missing:
            reason = "not_assessed_or_missing_gene:" + ";".join(missing)
        if reason:
            responses.append({**record, "status": "not_evaluated", "reason": reason})
            checks["not_evaluated_units"].append({"unit_id": unit_id, "reason": reason})
            continue
        responses.append({**record, "raw_text": json.dumps({"claims": []}), "status": "completed"})

    coverage_rows = list(member_coverage.values())
    coverage_counts = {state: sum(row["state"] == state for row in coverage_rows)
                       for state in ("assessed", "not_assessed", "invalid", "retrieval_failed")}
    partial_units = sorted({row["unit_id"] for row in coverage_rows if row["state"] != "assessed"})
    completion_gate = {"state": "partial" if partial_units else "complete", "may_claim_all_genes_assessed": not partial_units,
                       "partial_units": partial_units, "member_counts": coverage_counts,
                       "report_artifact": "member_coverage.json",
                       "rule": "Final case reports must use this gate; unit transport completed does not establish scientific completion."}
    coverage_path = out / "member_coverage.json"
    _write_json(coverage_path, {"contract": "agentic-route-b-member-coverage-1.0.0", "case_id": bundle["case_id"],
                              "bundle_sha256": bundle_sha256, "light_index_sha256": light_index_sha256,
                              "results": {bid: data["results_sha256"] for bid, data in checks["batches"].items()},
                              "members": coverage_rows, "counts": coverage_counts, "partial_units": partial_units,
                              "scientific_completion_gate": completion_gate,
                              "transport_completed_is_scientific_completion": False})
    out_metadata: dict[str, Any] = {}
    source_verify_records: list[dict[str, Any]] = []
    _cover_cache(cache, out_metadata, source_verify_records)
    _write_json(out / "metadata_cache.json", out_metadata)
    _write_json(out / "source_verify_cache.json", {"records": source_verify_records})
    discovery_path = None
    if discoveries:
        discovery_path = out / "discoveries.json"
        _write_json(discovery_path, {"case_id": bundle["case_id"],
                                     "parent_bundle_sha256": bundle_sha256,
                                     "discovered_sources": discoveries})
    _write_json(out / "agent_responses.json", {
        "case_id": bundle["case_id"], "executor": "light-route-b-convert", "model": "light",
        "responses": responses,
    })
    check_path = out / "light_check.json"
    _write_json(check_path, {"contract": index["contract"], "case_id": bundle["case_id"], **checks, "totals": totals})
    receipt = {
        "contract": index["contract"],
        "case_id": bundle["case_id"],
        "phenotype": bundle["phenotype"],
        "metadata_mode": "offline_cache_only" if offline else "cache_then_europe_pmc_crossref_doi",
        "input": {
            "bundle": {"path": str(bundle_path), "sha256": bundle_sha256},
            "light_index": {"path": str(index_path), "sha256": light_index_sha256},
            "metadata_cache": ({"path": str(Path(cache_arg).resolve()), "sha256": _sha(Path(cache_arg))}
                               if cache_arg else None),
            "results": {bid: {"file": f"{bid}/results.jsonl", "sha256": batch_checks["results_sha256"]}
                        for bid, batch_checks in checks["batches"].items()},
        },
        "outputs": {
            "member_coverage": {"file": coverage_path.name, "sha256": _sha(coverage_path)},
            "agent_responses": {"file": "agent_responses.json", "sha256": _sha(out / "agent_responses.json")},
            "discoveries": ({"file": discovery_path.name, "sha256": _sha(discovery_path)}
                            if discovery_path is not None else None),
            "metadata_cache": {"file": "metadata_cache.json", "sha256": _sha(out / "metadata_cache.json")},
            "source_verify_cache": {"file": "source_verify_cache.json",
                                    "sha256": _sha(out / "source_verify_cache.json")},
            "light_check": {"file": "light_check.json", "sha256": _sha(check_path)},
        },
        "discovery_state": ("written" if discovery_path is not None
                            else "none_required: every kept paper id is already frozen in its unit"),
        "counts": totals,
        "member_coverage": {"counts": coverage_counts, "partial_units": partial_units},
        "scientific_completion_gate": completion_gate,
        "checks": {"batches": checks["batches"], "dropped_claims": checks["dropped_claims"],
                   "tier_comparison": checks["tier_comparison"],
                   "not_evaluated_units": checks["not_evaluated_units"],
                   "gene_issue_count": len(checks["gene_issues"]),
                   "paper_issue_count": len(checks["paper_issues"])},
        "rerun": "A rerun can reuse the fetched metadata: pass this directory's metadata_cache.json as --cache.",
        "next_commands": _next_commands(bundle_path, out, discovery_path is not None),
        "limitations": ("Machine checks are hallucination guards only. A checked paper proves identity and "
                        "retrievability, not that the paper supports the stated bridge. No gene is selected, "
                        "excluded or scored here."),
    }
    _write_json(out / "light_convert_receipt.json", receipt)
    return {"contract": index["contract"], "case_id": bundle["case_id"], "out": str(out.resolve()),
            "member_coverage": {"path": str(coverage_path), "counts": coverage_counts, "partial_units": partial_units},
            "scientific_completion_gate": completion_gate,
            "counts": totals, "discoveries": str(discovery_path) if discovery_path else "",
            "agent_responses": str(out / "agent_responses.json"),
            "receipt_sha256": _sha(out / "light_convert_receipt.json"),
            "receipt": str(out / "light_convert_receipt.json")}


def verify_light_completion(receipt_arg: str, receipt_sha256: str) -> dict[str, Any]:
    """Final report gate: verify the sealed transport/coverage pair and inputs."""
    receipt_path = Path(receipt_arg).resolve()
    if not receipt_sha256 or _sha(receipt_path) != receipt_sha256:
        raise ValueError("light completion receipt SHA-256 mismatch")
    receipt = _read_json(receipt_path)
    for name in ("bundle", "light_index"):
        item = receipt["input"][name]
        if _sha(Path(item["path"])) != item["sha256"]:
            raise ValueError(f"light completion {name} hash mismatch")
    for name in ("member_coverage", "agent_responses"):
        item = receipt["outputs"][name]
        path = (receipt_path.parent / item["file"]).resolve()
        if path.parent != receipt_path.parent or _sha(path) != item["sha256"]:
            raise ValueError(f"light completion {name} hash mismatch")
    index_root = Path(receipt["input"]["light_index"]["path"]).parent
    for bid, item in receipt["input"]["results"].items():
        path = index_root / item["file"]
        if (_sha(path) if path.is_file() else None) != item["sha256"]:
            raise ValueError(f"light completion results hash mismatch: {bid}")
    coverage = _read_json(receipt_path.parent / receipt["outputs"]["member_coverage"]["file"])
    if coverage["scientific_completion_gate"] != receipt["scientific_completion_gate"]:
        raise ValueError("light completion gate differs from sealed coverage")
    return coverage["scientific_completion_gate"]


def _next_commands(bundle_path: Path, out: Path, has_discoveries: bool) -> list[dict[str, str]]:
    bundle_sha256 = _sha(bundle_path)
    commands = []
    bundle_dir = bundle_path.parent
    bundle_ref = "<BUNDLE_EXT_DIR>/route_b_semantic_bundle.json" if has_discoveries else str(bundle_path)
    if has_discoveries:
        commands.append({
            "step": "extend",
            "command": ("candidate-qtl route-b-semantic-bundle "
                        f"--extend-bundle {bundle_path} --parent-bundle-sha256 {bundle_sha256} "
                        f"--research-discoveries {out / 'discoveries.json'} --out-dir <BUNDLE_EXT_DIR>"),
        })
        prompts = "<BUNDLE_EXT_DIR>/route_b_semantic_prompts.jsonl"
    else:
        prompts = str(bundle_dir / "route_b_semantic_prompts.jsonl")
    commands.append({
        "step": "llm-run",
        "command": ("candidate-qtl llm-run --provider codex-agent --model light "
                    f"--agent-file {out / 'agent_responses.json'} --prompts {prompts} --out-dir <LLM_RUN_DIR>"),
    })
    commands.append({
        "step": "import",
        "command": (f"candidate-qtl route-b-semantic-import --bundle {bundle_ref} "
                    f"--prompts {prompts} --responses <LLM_RUN_DIR>/llm_candidate_responses.json "
                    "--llm-manifest <LLM_RUN_DIR>/llm_run_manifest.json --out-dir <IMPORT_DIR>"),
    })
    commands.append({
        "step": "source-verify",
        "command": ("candidate-qtl route-b-semantic-source-verify "
                    "--proposals <IMPORT_DIR>/semantic_claim_proposals.tsv "
                    f"--bundle {bundle_ref} "
                    f"--cache {out / 'source_verify_cache.json'} --out-dir <SOURCE_DIR>"),
    })
    commands.append({
        "step": "transfer-verify",
        "command": ("candidate-qtl route-b-semantic-transfer-verify "
                    "--proposals <IMPORT_DIR>/semantic_claim_proposals.tsv "
                    "--source-verification <SOURCE_DIR>/semantic_claim_source_verification.tsv "
                    f"--bundle {bundle_ref} --out-dir <TRANSFER_DIR>"),
    })
    commands.append({
        "step": "freeze",
        "command": ("candidate-qtl route-b-semantic-candidates "
                    "--claims <TRANSFER_DIR>/semantic_claims_verified.tsv --out-dir <FREEZE_DIR>"),
    })
    commands.append({"step": "final-report-completion-gate", "command": "",
                     "instruction": ("Before a final case report, call scripts.agentic_route_b_light.verify_light_completion "
                                     f"with {out / 'light_convert_receipt.json'} and the CLI-returned receipt_sha256. "
                                     "Consume its state/counts/partial_units; canonical unit completed alone is insufficient.")})
    return commands
