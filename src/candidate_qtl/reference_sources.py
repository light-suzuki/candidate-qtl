"""Phenotype-rich, source-backed reference records for Routes A/B/C (issue #232).

The corpus that ``reference-corpus-build`` could produce so far was a protein
*description* skeleton: a BLAST defline gives an accession, a species, a symbol
and a product name, and nothing else. Real-data testing showed exactly what that
costs -- a fasciation query ("fasciated stem with enlarged shoot apical
meristem") ranked NAC proteins at the top, because "No Apical Meristem" is a
domain *name*, while the CLAVATA machinery that actually produces the phenotype
scored zero. No retrieval backend fixes that: the biology was never in the
corpus.

What the routes need is not a better embedding over thin text but the record a
curator already wrote: what the gene does, **what its mutants look like**, which
pathway it sits in, and which experiment says so. UniProtKB carries all of that
in structured comment blocks (``FUNCTION``, ``DISRUPTION PHENOTYPE``,
``PATHWAY``, ``TISSUE SPECIFICITY``, ``DEVELOPMENTAL STAGE``), each tagged with
an ECO evidence code and the PubMed id it came from.

This module retrieves and parses those records:

``uniprot_fetch``
    Download one source-backed table per taxon over the documented REST API,
    recording the query, the UniProt release, the result count, the retrieval
    time and a SHA-256 for every file. Nothing is fetched without an explicit
    ``--allow-external-api``.
``read_uniprot_table``
    Parse a downloaded table into corpus rows, keeping the phenotype text apart
    from the function text, resolving each comment block's ECO codes into an
    evidence type, and carrying the **PubMed ids cited by the phenotype block
    itself** rather than the entry's whole bibliography.

Two rules hold throughout:

* **No phenotype argument exists here.** The fetch command takes taxa and a
  reviewed/unreviewed switch -- there is deliberately no free-text query, so a
  case phenotype cannot be smuggled into corpus construction and Route C cannot
  quietly become a second copy of Route A.
* **Missing is missing.** An entry with no curated disruption phenotype gets an
  empty ``mutant_phenotype``; that is ``knowledge_unavailable`` for the routes
  downstream, never evidence that the gene has no phenotype.
"""

from __future__ import annotations

import json
import re
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from . import __version__
from .core import runtime_path, safe_id, sha256

SCHEMA_VERSION = "1.0"
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
USER_AGENT = f"candidate-qtl/{__version__} reference-corpus-fetch (read-only reference lookup)"
PRESETS_RESOURCE = "refcorpus_presets.json"

# Requested in this order; the human-readable header UniProt returns for each is
# the key of ``UNIPROT_HEADERS`` below.
UNIPROT_FIELDS = [
    "accession", "id", "gene_primary", "gene_names", "protein_name", "organism_name", "organism_id",
    "protein_families", "cc_function", "cc_disruption_phenotype", "cc_developmental_stage",
    "cc_tissue_specificity", "cc_pathway", "go_id", "go_p", "lit_pubmed_id", "xref_tair",
]
UNIPROT_HEADERS = {
    "Entry": "accession", "Entry Name": "id", "Gene Names (primary)": "gene_primary",
    "Gene Names": "gene_names", "Protein names": "protein_name", "Organism": "organism_name",
    "Organism (ID)": "organism_id", "Protein families": "protein_families",
    "Function [CC]": "cc_function", "Disruption phenotype": "cc_disruption_phenotype",
    "Developmental stage": "cc_developmental_stage", "Tissue specificity": "cc_tissue_specificity",
    "Pathway": "cc_pathway", "Gene Ontology IDs": "go_id",
    "Gene Ontology (biological process)": "go_p", "PubMed ID": "lit_pubmed_id", "TAIR": "xref_tair",
}

FETCH_MANIFEST_NAME = "reference_fetch_manifest.json"

_ECO_BLOCK = re.compile(r"\{ECO:[^}]*\}")
_ECO_CODE = re.compile(r"ECO:(\d{7})")
_PUBMED = re.compile(r"PubMed:(\d+)", re.IGNORECASE)
_PREFIX = re.compile(
    r"\b(FUNCTION|DISRUPTION PHENOTYPE|PATHWAY|TISSUE SPECIFICITY|DEVELOPMENTAL STAGE|MISCELLANEOUS):\s*",
    re.IGNORECASE,
)
_LINK_NEXT = re.compile(r'<([^>]+)>;\s*rel="next"')

# ECO code -> evidence type, strongest first. A phenotype asserted from a
# published experiment is a different object from one inferred by similarity, and
# the routes downstream must be able to tell them apart.
EVIDENCE_TYPES = (
    ("experimental", {"0000269", "0000314", "0000315", "0000316", "0007005"}),
    ("curator_inference", {"0000305"}),
    ("sequence_similarity", {"0000250", "0000255", "0000312", "0007744"}),
    ("automatic_annotation", {"0000256", "0000259", "0000213", "0000501"}),
)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def clean_comment(value: object) -> str:
    """A UniProt comment block as readable text: no ECO braces, no block prefix."""
    text = _ECO_BLOCK.sub(" ", str(value or ""))
    text = _PREFIX.sub("", text)
    return re.sub(r"\s+", " ", text).strip(" ;.").strip()


def evidence_codes(value: object) -> list[str]:
    return sorted({f"ECO:{code}" for code in _ECO_CODE.findall(str(value or ""))})


def evidence_type(codes: list[str]) -> str:
    """Strongest evidence type present, or ``unspecified`` when nothing is tagged."""
    numbers = {code.split(":")[-1] for code in codes}
    for name, members in EVIDENCE_TYPES:
        if numbers & members:
            return name
    return "unspecified"


def pubmed_ids(value: object) -> list[str]:
    return sorted({identifier for identifier in _PUBMED.findall(str(value or ""))}, key=int)


def _split_ids(value: object, prefix: str = "") -> list[str]:
    items = [item.strip() for item in re.split(r"[;,|]+", str(value or "")) if item.strip()]
    return [f"{prefix}{item}" if prefix and not item.startswith(prefix) else item for item in items]


# --------------------------------------------------------------------------- #
# presets
# --------------------------------------------------------------------------- #

def _presets_path() -> Path:
    """Package resource path of ``refcorpus_presets.json`` (no repo-root copy)."""
    return Path(__file__).resolve().parent / "resources" / PRESETS_RESOURCE


def load_refcorpus_presets(path_arg: str | None = None) -> dict[str, dict[str, object]]:
    """Load the packaged, versioned corpus presets and validate their schema.

    Returns ``{preset_id: preset}``. A preset names the taxa a fetch should
    cover, whether each taxon needs unreviewed (TrEMBL) entries, and the gene
    symbols an acceptance gate must find in a corpus built from it (issue #299).
    """
    path = runtime_path(path_arg) if path_arg else _presets_path()
    if not path.exists():
        raise FileNotFoundError(f"reference corpus presets not found at {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"reference corpus presets {path} is not valid JSON: {error}") from error
    if not isinstance(payload, dict) or not isinstance(payload.get("presets"), list) or not payload["presets"]:
        raise ValueError(f"reference corpus presets {path} must be an object with a non-empty 'presets' list")
    presets: dict[str, dict[str, object]] = {}
    for preset in payload["presets"]:
        if not isinstance(preset, dict):
            raise ValueError(f"reference corpus presets {path}: every preset must be an object")
        preset_id = str(preset.get("preset_id", "")).strip()
        if not preset_id:
            raise ValueError(f"reference corpus presets {path}: a preset is missing 'preset_id'")
        if preset_id in presets:
            raise ValueError(f"reference corpus presets {path}: duplicate preset_id '{preset_id}'")
        if not str(preset.get("provider", "")).strip():
            raise ValueError(f"preset '{preset_id}': needs a non-empty 'provider'")
        taxa = preset.get("taxa")
        if not isinstance(taxa, list) or not taxa:
            raise ValueError(f"preset '{preset_id}': 'taxa' must be a non-empty list")
        for spec in taxa:
            if not isinstance(spec, dict) or not str(spec.get("species", "")).strip():
                raise ValueError(f"preset '{preset_id}': every taxon needs a non-empty 'species'")
            taxon_id = spec.get("taxon_id")
            if not isinstance(taxon_id, int) and not (isinstance(taxon_id, str) and str(taxon_id).strip().isdigit()):
                raise ValueError(f"preset '{preset_id}': taxon '{spec.get('species')}' needs a numeric 'taxon_id'")
            if not isinstance(spec.get("reviewed"), bool):
                raise ValueError(f"preset '{preset_id}': taxon '{spec.get('species')}' needs a boolean 'reviewed'")
        for field in ("required_terms", "optional_terms"):
            terms = preset.get(field, [])
            if not isinstance(terms, list) or any(not isinstance(term, str) or not str(term).strip() for term in terms):
                raise ValueError(f"preset '{preset_id}': '{field}' must be a list of non-empty symbols")
        presets[preset_id] = preset
    return presets


def expand_preset(preset_arg: str, path_arg: str | None = None) -> list[dict[str, object]]:
    """Resolve one preset id to its taxon list (``taxon_id`` as a string, like ``--taxon``)."""
    presets = load_refcorpus_presets(path_arg)
    preset = presets.get(str(preset_arg).strip())
    if preset is None:
        raise ValueError(f"unknown reference corpus preset '{preset_arg}'; available: {', '.join(sorted(presets))}")
    return [{
        "taxon_id": str(spec["taxon_id"]), "species": str(spec["species"]).strip(),
        "reviewed": bool(spec["reviewed"]), "note": str(spec.get("note", "")).strip(),
    } for spec in preset["taxa"]]


# --------------------------------------------------------------------------- #
# fetch
# --------------------------------------------------------------------------- #

def _request(url: str, timeout: int = 120) -> tuple[str, dict[str, str]]:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            payload = response.read().decode("utf-8")
            headers = {key.lower(): value for key, value in response.headers.items()}
    except urllib.error.HTTPError as error:
        raise RuntimeError(f"UniProt HTTP {error.code} for {url}") from error
    except urllib.error.URLError as error:
        raise RuntimeError(f"UniProt connection failed for {url}: {error.reason}") from error
    return payload, headers


def uniprot_fetch(
    taxon_args: list[str],
    out_dir_arg: str,
    allow_external_api: bool = False,
    include_unreviewed: bool = False,
    max_rows: int = 0,
    page_size: int = 500,
    preset_arg: str | None = None,
) -> dict[str, object]:
    """Download one reviewed reference table per taxon, with its release recorded.

    There is no free-text query argument on purpose: the only inputs are which
    organisms to retrieve and whether unreviewed entries are included. A corpus
    built from the case phenotype would make Route C a restatement of Route A.

    ``preset_arg`` (issue #299) expands a versioned preset from the packaged
    ``refcorpus_presets.json`` into taxa, each carrying its own reviewed flag
    (``reviewed: false`` means that taxon needs TrEMBL entries -- the
    nodulation vocabulary of Medicago/Lotus lives there). ``--taxon`` and
    ``--preset`` union; an explicit ``--taxon`` wins over the same id from the
    preset, and ``--include-unreviewed`` widens every taxon.
    """
    taxa = [str(value).strip() for value in (taxon_args or []) if str(value).strip()]
    preset_taxa: dict[str, dict[str, object]] = {}
    preset_id = str(preset_arg or "").strip() or None
    if preset_id:
        for spec in expand_preset(preset_id):
            # An explicit --taxon wins over the same id from the preset: the
            # caller's reviewed/unreviewed choice is the more specific one, so a
            # preset-listed taxon that is also on --taxon is fetched with the
            # explicit semantics and keeps no preset provenance.
            if str(spec["taxon_id"]) not in taxa:
                taxa.append(str(spec["taxon_id"]))
                preset_taxa[str(spec["taxon_id"])] = spec
    if not taxa:
        raise ValueError("reference-corpus-fetch requires at least one --taxon or --preset")
    if not allow_external_api:
        raise ValueError(
            "reference-corpus-fetch reaches UniProtKB over the network; pass --allow-external-api to allow it. "
            "Only the taxon id and the field list are transmitted: no phenotype, interval or candidate leaves this machine"
        )
    if page_size < 1 or page_size > 500:
        raise ValueError("--page-size must be between 1 and 500 (the UniProt REST maximum)")
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)

    sources: list[dict[str, object]] = []
    outputs: list[str] = []
    for taxon in taxa:
        preset_spec = preset_taxa.get(taxon)
        include_unreviewed_taxon = include_unreviewed or bool(preset_spec and not preset_spec["reviewed"])
        query = f"organism_id:{taxon}" if include_unreviewed_taxon else f"organism_id:{taxon} AND reviewed:true"
        url = f"{UNIPROT_SEARCH}?" + urllib.parse.urlencode(
            {"query": query, "format": "tsv", "fields": ",".join(UNIPROT_FIELDS), "size": page_size})
        header_line = ""
        lines: list[str] = []
        release = ""
        total = ""
        pages = 0
        while url:
            payload, headers = _request(url)
            pages += 1
            release = release or headers.get("x-uniprot-release", "")
            total = total or headers.get("x-total-results", "")
            page_lines = payload.splitlines()
            if page_lines:
                header_line = header_line or page_lines[0]
                lines.extend(page_lines[1:])
            if max_rows and len(lines) >= max_rows:
                lines = lines[:max_rows]
                break
            match = _LINK_NEXT.search(headers.get("link", ""))
            url = match.group(1) if match else ""
        if not header_line:
            raise RuntimeError(f"UniProt returned no table for taxon {taxon}")
        path = out / f"uniprot_{safe_id(taxon)}.tsv"
        path.write_text("\n".join([header_line, *lines]) + "\n", encoding="utf-8")
        outputs.append(str(path))
        sources.append({
            "taxon_id": taxon, "query": query, "url": url or "(last page)",
            "fields": list(UNIPROT_FIELDS), "source_database": "UniProtKB",
            "source_release": release, "reported_total_results": total,
            "rows_written": len(lines), "pages": pages, "reviewed_only": not include_unreviewed_taxon,
            # issue #299: preset provenance, additive so older manifest entries stay readable.
            "preset_id": preset_id if preset_spec else None,
            "reviewed": bool(preset_spec["reviewed"]) if preset_spec else not include_unreviewed,
            "max_rows": max_rows or None, "retrieved_at_utc": _now(),
            "path": str(path), "sha256": sha256(path),
        })

    manifest_path = out / FETCH_MANIFEST_NAME
    if manifest_path.exists():
        # Species are usually added one run at a time. Replacing the manifest
        # would strip the release and retrieval time off every table fetched
        # earlier, and those tables would silently lose their provenance.
        try:
            previous = json.loads(manifest_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            previous = {}
        fetched_now = {str(item["path"]) for item in sources}
        kept = [item for item in previous.get("sources", [])
                if isinstance(item, dict) and str(item.get("path", "")) not in fetched_now]
        sources = sorted([*kept, *sources], key=lambda item: str(item.get("taxon_id", "")))
        outputs = sorted({*outputs, *(str(item["path"]) for item in kept)})

    manifest = {
        "schema_version": SCHEMA_VERSION, "command": "reference-corpus-fetch",
        "candidate_qtl_version": __version__, "fetched_at_utc": _now(),
        "provider": "uniprotkb_rest", "endpoint": UNIPROT_SEARCH,
        "case_independent": True, "phenotype_used": False,
        "free_text_query_argument": False,
        "independence_rule": (
            "the only inputs are taxon ids and the reviewed/unreviewed switch; there is no phenotype, interval "
            "or free-text query argument, so a case cannot shape the corpus it is later evaluated against"
        ),
        "external_transmission": "taxon ids and the UniProt field list only",
        "preset": preset_id,
        "sources": sources,
        "row_total": sum(int(item.get("rows_written", 0)) for item in sources),
        "taxa": sorted({str(item.get("taxon_id", "")) for item in sources if str(item.get("taxon_id", ""))}),
        "next_step": "reference-corpus-build --from-uniprot <each table> --corpus-version <label>",
        "outputs": outputs,
    }
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), "tables": outputs, "row_total": manifest["row_total"],
            "taxa": manifest["taxa"],
            "releases": sorted({str(item.get("source_release", "")) for item in sources if item.get("source_release")})}


# --------------------------------------------------------------------------- #
# parse
# --------------------------------------------------------------------------- #

def _normalise_columns(frame: pd.DataFrame, path: Path) -> pd.DataFrame:
    renamed: dict[str, str] = {}
    for column in frame.columns:
        name = str(column).strip()
        field = UNIPROT_HEADERS.get(name, name if name in set(UNIPROT_FIELDS) else "")
        if field:
            renamed[column] = field
    frame = frame.rename(columns=renamed)
    if "accession" not in frame.columns:
        raise ValueError(
            f"UniProt table {path} has no Entry/accession column; export it with "
            f"fields={','.join(UNIPROT_FIELDS)}"
        )
    return frame


def _release_for(path: Path) -> tuple[str, str]:
    """(release, retrieved_at) recorded for this table by ``reference-corpus-fetch``."""
    manifest_path = path.parent / FETCH_MANIFEST_NAME
    if not manifest_path.exists():
        return "", ""
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return "", ""
    for source in manifest.get("sources", []):
        if isinstance(source, dict) and Path(str(source.get("path", ""))).name == path.name:
            return str(source.get("source_release", "")), str(source.get("retrieved_at_utc", ""))
    return "", ""


def read_uniprot_table(path: Path) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    """Parse one downloaded UniProt table into (rows, excluded, source metadata)."""
    frame = _normalise_columns(pd.read_csv(path, sep="\t", dtype=str).fillna(""), path)
    release, retrieved_at = _release_for(path)
    rows: list[dict[str, object]] = []
    excluded: list[dict[str, object]] = []
    for record in frame.to_dict("records"):
        accession = str(record.get("accession", "")).strip()
        if not accession:
            continue
        symbol = str(record.get("gene_primary", "")).strip() or (
            _split_ids(record.get("gene_names", ""))[0] if str(record.get("gene_names", "")).strip() else "")
        symbol = symbol.split()[0] if symbol else ""
        phenotype_raw = record.get("cc_disruption_phenotype", "")
        function_raw = record.get("cc_function", "")
        if not symbol:
            # No symbol means no family bucket, so the row could never be matched
            # to an interval later. Kept in the ledger with the phenotype it had,
            # so a curator can see what the naming gap cost.
            excluded.append({
                "reference_gene_id": accession, "provenance_tier": "retrieved", "source": str(path),
                "exclusion_reason": "no_gene_symbol_in_uniprot_record",
                "detail": clean_comment(phenotype_raw)[:200] or str(record.get("protein_name", ""))[:200],
            })
            continue
        codes = evidence_codes(phenotype_raw) or evidence_codes(function_raw)
        go_terms = _split_ids(record.get("go_id", ""))
        rows.append({
            "reference_gene_id": accession,
            "reference_species": str(record.get("organism_name", "")).strip(),
            "taxon_id": str(record.get("organism_id", "")).strip(),
            "gene_symbol": symbol,
            "description": str(record.get("protein_name", "")).strip(),
            "function_description": clean_comment(function_raw),
            "mutant_phenotype": clean_comment(phenotype_raw),
            "tissue_specificity": clean_comment(record.get("cc_tissue_specificity", "")),
            "developmental_stage": clean_comment(record.get("cc_developmental_stage", "")),
            "pathway": clean_comment(record.get("cc_pathway", "")),
            "protein_family": str(record.get("protein_families", "")).strip(),
            "go_terms": ";".join(go_terms),
            "keywords": clean_comment(record.get("go_p", "")),
            "source_id": f"UniProt:{accession}",
            "source_url": f"https://rest.uniprot.org/uniprotkb/{accession}",
            # The bibliography of the whole entry is not support for the phenotype
            # claim: only the PubMed ids the phenotype block itself cites are.
            "phenotype_source_ids": ";".join(f"PMID:{value}" for value in pubmed_ids(phenotype_raw)),
            "literature_source_ids": ";".join(f"PMID:{value}" for value in _split_ids(record.get("lit_pubmed_id", ""))),
            "evidence_codes": ";".join(codes),
            "evidence_type": evidence_type(codes),
            "source_database": "UniProtKB",
            "source_release": release,
            "retrieved_at": retrieved_at,
            "provenance_tier": "retrieved",
            "provenance_detail": f"uniprot:{path.name}",
        })
    meta = {
        "source": str(path), "reader": "uniprot tsv", "source_database": "UniProtKB",
        "source_release": release, "retrieved_at_utc": retrieved_at,
        "rows": len(rows), "excluded": len(excluded),
        "phenotype_rows": sum(1 for row in rows if row["mutant_phenotype"]),
        "pathway_rows": sum(1 for row in rows if row["pathway"]),
        "experimental_evidence_rows": sum(1 for row in rows if row["evidence_type"] == "experimental"),
    }
    return rows, excluded, meta


# --------------------------------------------------------------------------- #
# coverage
# --------------------------------------------------------------------------- #

COVERAGE_FIELDS = (
    "description", "function_description", "mutant_phenotype", "pathway", "go_terms",
    "tissue_specificity", "developmental_stage", "protein_family", "family_id", "orthogroup_id",
    "phenotype_source_ids", "literature_source_ids", "evidence_codes",
)
COVERAGE_COLUMNS = ["scope", "reference_species", "field", "rows", "non_empty", "coverage"]
ALL_SPECIES = "__all_species__"


def corpus_coverage(corpus: pd.DataFrame) -> pd.DataFrame:
    """Per-species and overall field coverage.

    Row count alone hides the failure this issue is about: a corpus can be large
    and still carry no mutant phenotype at all. Coverage is reported per field so
    "how much biology is in here" is answerable without opening the table.
    """
    rows: list[dict[str, object]] = []

    def add(scope: str, species: str, frame: pd.DataFrame) -> None:
        for field in COVERAGE_FIELDS:
            values = frame[field] if field in frame.columns else pd.Series(dtype=str)
            non_empty = int((values.astype(str).str.strip() != "").sum()) if len(frame) else 0
            rows.append({"scope": scope, "reference_species": species, "field": field,
                         "rows": int(len(frame)), "non_empty": non_empty,
                         "coverage": round(non_empty / len(frame), 6) if len(frame) else 0.0})

    add("all", ALL_SPECIES, corpus)
    if "reference_species" in corpus.columns:
        for species in sorted({str(value).strip() for value in corpus["reference_species"] if str(value).strip()}):
            add("species", species, corpus[corpus["reference_species"].astype(str).str.strip() == species])
    return pd.DataFrame(rows, columns=COVERAGE_COLUMNS)


def coverage_summary(coverage: pd.DataFrame) -> dict[str, object]:
    overall = coverage[coverage["scope"] == "all"]
    return {str(row["field"]): {"non_empty": int(row["non_empty"]), "rows": int(row["rows"]),
                                "coverage": float(row["coverage"])}
            for row in overall.to_dict("records")}
