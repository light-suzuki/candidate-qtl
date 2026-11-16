"""Build the reference resource that Routes B and C both need.

The monitored audit found Route B degraded to broad-function bins and Route C
unable to run at all. Neither failure was a logic failure: the reference side --
what a gene is, what its mutants look like, which family it belongs to, and who
says so -- simply did not exist on this workstation. Both routes consume that
same table, so one resource unblocks both.

**"Let a model build it" means retrieve and structure, never write from memory.**
Feeding model recollection into the routes whose job is to catch fabrication
would defeat the whole design, so every row carries a provenance tier and only
some tiers are usable:

``retrieved``
    Parsed from an authoritative source: a downloaded UniProtKB table
    (``--from-uniprot``, the phenotype-rich path added by issue #232), a local
    UniProt-derived protein database, or an exported dump. No model involvement.
    Usable.
``llm_verified``
    A model or agent proposed the row, and its accession was confirmed against a
    local authoritative record *and* its citations resolved. Usable, and marked
    as model-proposed for the rest of its life.
``llm_unverified``
    Proposed but not confirmed. Kept in the exclusion ledger with the reason and
    **never written to the corpus**, so it cannot reach an index or a Route B
    source-backed count.

**The corpus is case-independent by construction.** This command takes no
phenotype, no interval and no candidate argument. A corpus built from the case
phenotype would make Route C a second copy of Route A rather than an independent
route, and the convergence between them would be an artefact. Building it once,
phenotype-blind, and freezing it with a version is what keeps the routes
independent.

In benchmark mode a row that names a focal-lineage / answer-shaped marker is
barred, and the manifest records that the corpus contains model-proposed rows so
any number produced with it is labelled rather than quietly reported.
"""

from __future__ import annotations

import json
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from . import __version__
from .core import file_fingerprint, runtime_path, run_tool, safe_id, sha256
from .family_index import symbol_stem
from .identifiers import is_usable_symbol
from .query_eligibility import PISUM_MARKERS
from .reference_sources import corpus_coverage, coverage_summary, load_refcorpus_presets, read_uniprot_table
from .source_verify import classify_citation

SCHEMA_VERSION = "1.0"
CORPUS_FORMAT_VERSION = "1.0.0"

PROVENANCE_TIERS = ("retrieved", "llm_verified", "llm_unverified")
USABLE_TIERS = ("retrieved", "llm_verified")

CORPUS_COLUMNS = [
    "reference_gene_id", "reference_species", "taxon_id", "gene_symbol",
    "family_id", "orthogroup_id", "description", "function_description",
    "mutant_phenotype", "pathway", "pathway_id", "go_terms", "keywords",
    "source_id", "source_url", "provenance_tier", "provenance_detail",
    "verified_source_ids", "corpus_version",
    # issue #232: the curated biology a description-only corpus never carried.
    "tissue_specificity", "developmental_stage", "protein_family",
    "phenotype_source_ids", "literature_source_ids", "evidence_codes", "evidence_type",
    "source_database", "source_release", "retrieved_at",
]
# Fields whose presence decides whether a row carries biology or only a name.
# Used to pick the richer row when two sources describe the same accession.
_RICHNESS_COLUMNS = ("mutant_phenotype", "function_description", "pathway", "go_terms",
                     "tissue_specificity", "developmental_stage", "protein_family")
EXCLUSION_COLUMNS = ["reference_gene_id", "provenance_tier", "source", "exclusion_reason", "detail"]

# UniProt-style FASTA defline: "<description> OS=<species> OX=<taxon> GN=<symbol> PE=.. SV=.."
_DEFLINE = re.compile(
    r"^(?P<description>.*?)(?:\s+OS=(?P<species>.*?))?(?:\s+OX=(?P<taxon>\d+))?"
    r"(?:\s+GN=(?P<symbol>\S+))?(?:\s+PE=\d+)?(?:\s+SV=\d+)?\s*$"
)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _parse_defline(title: str) -> dict[str, str]:
    match = _DEFLINE.match(str(title).strip())
    if not match:
        return {"description": str(title).strip(), "species": "", "taxon": "", "symbol": ""}
    return {key: (match.group(key) or "").strip() for key in ("description", "species", "taxon", "symbol")}


def _row(base: dict[str, object] | None = None, **values: object) -> dict[str, object]:
    """One corpus row: every column present, ``base`` first, explicit values last.

    Keeping the override order explicit matters -- a model-proposed row may carry
    a species or family of its own, and the authoritative value must win.
    """
    row: dict[str, object] = {column: "" for column in CORPUS_COLUMNS}
    for source in (base or {}, values):
        row.update({key: value for key, value in source.items() if key in row})
    return row


def _from_blastdb(database: Path, binary: str) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    """Read every entry of a local protein BLAST database as retrieved rows.

    This is the honest skeleton: accession, description, species, taxon and gene
    symbol all come out of the database that is already on the machine, so the
    bulk of the corpus needs neither a network call nor a model.
    """
    command = [binary, "-db", str(database), "-entry", "all", "-outfmt", "%a\t%t"]
    completed = run_tool(command, text=True, capture_output=True)
    if completed.returncode:
        raise RuntimeError(f"blastdbcmd failed for {database}: {completed.stderr.strip()}")
    rows: list[dict[str, object]] = []
    excluded: list[dict[str, object]] = []
    for line in completed.stdout.splitlines():
        if not line.strip():
            continue
        accession, _, title = line.partition("\t")
        accession = accession.strip()
        parsed = _parse_defline(title)
        symbol = parsed["symbol"]
        if not accession:
            continue
        if not symbol:
            # No gene symbol means no family bucket and nothing to match an
            # interval against later. Kept in the ledger, never silently dropped.
            excluded.append({"reference_gene_id": accession, "provenance_tier": "retrieved",
                             "source": str(database), "exclusion_reason": "no_gene_symbol_in_defline",
                             "detail": parsed["description"][:200]})
            continue
        rows.append(_row(
            reference_gene_id=accession, reference_species=parsed["species"], taxon_id=parsed["taxon"],
            gene_symbol=symbol, family_id=symbol_stem(symbol), description=parsed["description"],
            source_id=f"UniProt:{accession}", source_url=f"https://rest.uniprot.org/uniprotkb/{accession}",
            provenance_tier="retrieved", provenance_detail=f"blastdb:{database.name}",
        ))
    meta = {
        "source": str(database), "reader": "blastdbcmd -entry all",
        "command": " ".join(command), "rows": len(rows), "excluded": len(excluded),
        "family_id_rule": "gene symbol stem (CLV1 -> CLV); a naming bucket, never an evolutionary family claim",
    }
    return rows, excluded, meta


def _from_tsv(path: Path) -> tuple[list[dict[str, object]], dict[str, object]]:
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if "reference_gene_id" not in frame.columns:
        raise ValueError(f"reference dump {path} requires reference_gene_id")
    rows = []
    for record in frame.to_dict("records"):
        symbol = str(record.get("gene_symbol", "")).strip()
        rows.append(_row(
            {key: str(record.get(key, "")).strip() for key in CORPUS_COLUMNS if key in record},
            provenance_tier="retrieved",
            provenance_detail=str(record.get("provenance_detail", "")) or f"dump:{path.name}",
            family_id=str(record.get("family_id", "")).strip() or (symbol_stem(symbol) if symbol else ""),
        ))
    return rows, {"source": str(path), "reader": "tsv dump", "rows": len(rows)}


def _payload_rows(path: Path) -> list[dict[str, object]]:
    if path.suffix.casefold() == ".tsv":
        return pd.read_csv(path, sep="\t", dtype=str).fillna("").to_dict("records")
    text = path.read_text(encoding="utf-8")
    if path.suffix.casefold() == ".jsonl":
        payloads = [json.loads(line) for line in text.splitlines() if line.strip()]
    else:
        parsed = json.loads(text)
        payloads = parsed if isinstance(parsed, list) else [parsed]
    records: list[dict[str, object]] = []
    for payload in payloads:
        if isinstance(payload, dict):
            for key in ("records", "rows", "entries", "genes"):
                value = payload.get(key)
                if isinstance(value, list):
                    records.extend(item for item in value if isinstance(item, dict))
                    break
            else:
                records.append(payload)
        elif isinstance(payload, list):
            records.extend(item for item in payload if isinstance(item, dict))
    return records


def preset_acceptance(corpus: pd.DataFrame, preset: dict[str, object]) -> dict[str, object]:
    """Check which preset terms exist as gene symbols in the corpus (issue #299).

    A term is ``found`` only when at least one row carries it as a usable gene
    symbol in the ``gene_symbol`` column -- never as a bare accession, AGI locus
    or locus tag (``is_usable_symbol``) -- compared with casefold so symbol case
    does not decide. Required terms failing this gate abort the build; optional
    terms are recorded only. This is the guard that would have caught the
    2026-08-09 failure, where a corpus was installed while the nodulation
    vocabulary silently stayed out of it.
    """
    symbols = corpus["gene_symbol"].astype(str).str.strip()
    folded = symbols[symbols.map(is_usable_symbol)].str.casefold()

    def check(term: str) -> dict[str, object]:
        matches = folded == str(term).strip().casefold()
        return {"found": bool(matches.any()), "count": int(matches.sum())}

    return {
        "preset_id": str(preset.get("preset_id", "")).strip(),
        "required_terms": {str(term): check(term) for term in preset.get("required_terms", [])},
        "optional_terms": {str(term): check(term) for term in preset.get("optional_terms", [])},
        "match_rule": (
            "casefold equality against the corpus gene_symbol column; only rows whose symbol "
            "passes is_usable_symbol count"
        ),
        "fail_closed_on_required": True,
    }


def _orthogroup_map(path_args: list[str] | None) -> dict[str, dict[str, str]]:
    """reference_gene_id -> {orthogroup_id, family_id} from a source-backed mapping.

    Kept separate from the record sources: a symbol stem is a naming bucket, and
    only an explicit mapping file may assert an orthogroup.
    """
    mapping: dict[str, dict[str, str]] = {}
    for value in path_args or []:
        path = runtime_path(value)
        frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        if "reference_gene_id" not in frame.columns:
            raise ValueError(f"orthogroup map {path} requires reference_gene_id")
        for record in frame.to_dict("records"):
            key = str(record["reference_gene_id"]).strip().casefold()
            if not key:
                continue
            entry = mapping.setdefault(key, {})
            for column in ("orthogroup_id", "family_id"):
                candidate = str(record.get(column, "")).strip()
                if candidate:
                    entry[column] = candidate
    return mapping


def reference_corpus_build(
    out_dir_arg: str,
    corpus_version: str,
    blastdb_args: list[str] | None = None,
    dump_args: list[str] | None = None,
    payload_args: list[str] | None = None,
    citation_cache_arg: str | None = None,
    benchmark_mode: bool = False,
    full_hash: bool = False,
    uniprot_args: list[str] | None = None,
    orthogroup_map_args: list[str] | None = None,
    preset_acceptance_arg: str | None = None,
) -> dict[str, object]:
    """Assemble a frozen, versioned, phenotype-blind reference corpus."""
    if not str(corpus_version).strip():
        raise ValueError("reference-corpus-build requires a non-empty --corpus-version")
    corpus_version = str(corpus_version).strip()
    if not (blastdb_args or dump_args or payload_args or uniprot_args):
        raise ValueError("reference-corpus-build requires at least one --from-uniprot / --from-blastdb / --from-dump / --from-payload input")
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    excluded: list[dict[str, object]] = []
    sources: list[dict[str, object]] = []

    binary = shutil.which("blastdbcmd") or ""
    for value in blastdb_args or []:
        if not binary:
            raise RuntimeError("blastdbcmd is not on PATH; cannot read a local protein database")
        database = runtime_path(value)
        retrieved, dropped, meta = _from_blastdb(database, binary)
        rows.extend(retrieved)
        excluded.extend(dropped)
        sources.append(meta)
    for value in uniprot_args or []:
        path = runtime_path(value)
        retrieved, dropped, meta = read_uniprot_table(path)
        rows.extend(_row(record, family_id=str(record.get("family_id", "")).strip()
                         or symbol_stem(str(record.get("gene_symbol", ""))))
                    for record in retrieved)
        excluded.extend(dropped)
        sources.append({**meta, "fingerprint": file_fingerprint(path, full_hash)})
    for value in dump_args or []:
        path = runtime_path(value)
        retrieved, meta = _from_tsv(path)
        rows.extend(retrieved)
        sources.append({**meta, "fingerprint": file_fingerprint(path, full_hash)})

    # Accessions proved by a retrieved source. A model-proposed row is only
    # promoted when its accession is one of these (or the row is dropped), so a
    # fabricated accession cannot enter the corpus.
    authoritative = {str(row["reference_gene_id"]).casefold(): row for row in rows}

    citation_cache: set[str] = set()
    if citation_cache_arg:
        cache_path = runtime_path(citation_cache_arg)
        files = sorted(cache_path.glob("*.json")) if cache_path.is_dir() else [cache_path]
        for item in files:
            payload = json.loads(item.read_text(encoding="utf-8"))
            entries = payload if isinstance(payload, list) else [payload]
            for entry in entries:
                if isinstance(entry, dict):
                    for key in ("pmid", "id", "doi", "pmcid"):
                        value = str(entry.get(key, "")).strip()
                        if value:
                            citation_cache.add(value.casefold())

    proposed = 0
    for value in payload_args or []:
        path = runtime_path(value)
        sources.append({"source": str(path), "reader": "model/agent payload",
                        "fingerprint": file_fingerprint(path, full_hash)})
        for record in _payload_rows(path):
            proposed += 1
            accession = str(record.get("reference_gene_id", "")).strip()
            citations = [item.strip() for item in re.split(r"[;,|]+", str(record.get("source_id", ""))) if item.strip()]
            blob = " ".join(str(record.get(key, "")) for key in
                            ("reference_gene_id", "reference_species", "gene_symbol", "description",
                             "function_description", "mutant_phenotype")).casefold()
            reason = ""
            if not accession:
                reason = "no_reference_gene_id"
            elif accession.casefold() not in authoritative:
                # Model memory is not an authority: the accession must appear in a
                # retrieved source before the row may be used.
                reason = "accession_not_in_any_retrieved_source"
            elif not citations:
                reason = "no_source_id_supplied"
            elif benchmark_mode and any(marker in blob for marker in PISUM_MARKERS):
                reason = "focal_lineage_or_answer_shaped_marker_in_benchmark_mode"
            verified: list[str] = []
            if not reason:
                for citation in citations:
                    citation_type, identifier = classify_citation(citation)
                    if citation_type == "database_accession" and identifier.casefold().endswith(accession.casefold()):
                        verified.append(citation)
                    elif citation_type in {"pmid", "pmcid", "doi"} and identifier.casefold() in citation_cache:
                        verified.append(citation)
                if not verified:
                    reason = "no_citation_resolved_against_the_local_cache"
            if reason:
                excluded.append({"reference_gene_id": accession or "(missing)", "provenance_tier": "llm_unverified",
                                 "source": str(path), "exclusion_reason": reason,
                                 "detail": str(record.get("description", ""))[:200]})
                continue
            anchor = authoritative[accession.casefold()]
            symbol = str(record.get("gene_symbol", "")).strip() or str(anchor["gene_symbol"])
            rows.append(_row(
                {key: str(record.get(key, "")).strip() for key in CORPUS_COLUMNS if key in record},
                reference_gene_id=accession,
                # Species and taxon come from the authoritative record, never from
                # the proposal: a model may not redefine what an accession is.
                reference_species=str(anchor["reference_species"]), taxon_id=str(anchor["taxon_id"]),
                gene_symbol=symbol, family_id=str(record.get("family_id", "")).strip() or symbol_stem(symbol),
                provenance_tier="llm_verified", provenance_detail=f"payload:{path.name}",
                verified_source_ids=";".join(verified),
            ))

    if not rows:
        raise ValueError("no usable reference record survived; see the exclusion ledger")

    corpus = pd.DataFrame(rows, columns=CORPUS_COLUMNS)
    corpus["corpus_version"] = corpus_version
    # A later model-verified row wins over the bare retrieved skeleton for the
    # same accession, because it carries the added function / phenotype text.
    corpus["_tier_rank"] = corpus["provenance_tier"].map({"llm_verified": 0, "retrieved": 1}).fillna(2)
    # Within a tier the richer record wins: the same accession read from a defline
    # and from a curated record is the same gene described twice, and the curated
    # phenotype is the reason this corpus exists.
    corpus["_richness"] = sum(
        (corpus[column].astype(str).str.strip() != "").astype(int) for column in _RICHNESS_COLUMNS)
    corpus = (corpus.sort_values(["reference_gene_id", "_tier_rank", "_richness"],
                                 ascending=[True, True, False], kind="stable")
              .drop_duplicates("reference_gene_id", keep="first")
              .drop(columns=["_tier_rank", "_richness"])
              .sort_values("reference_gene_id", kind="stable"))

    mapping = _orthogroup_map(orthogroup_map_args)
    if mapping:
        keys = corpus["reference_gene_id"].astype(str).str.casefold()
        for column in ("orthogroup_id", "family_id"):
            mapped = keys.map(lambda key: mapping.get(key, {}).get(column, ""))
            corpus[column] = mapped.where(mapped.astype(str) != "", corpus[column])

    # issue #299: the acceptance gate runs on the final, deduplicated corpus and
    # aborts before anything is written, so a corpus that lost its required
    # vocabulary can never silently replace a previously frozen one.
    acceptance: dict[str, object] | None = None
    if str(preset_acceptance_arg or "").strip():
        preset_id = str(preset_acceptance_arg).strip()
        presets = load_refcorpus_presets()
        if preset_id not in presets:
            raise ValueError(f"unknown reference corpus preset '{preset_id}'; available: {', '.join(sorted(presets))}")
        acceptance = preset_acceptance(corpus, presets[preset_id])
        missing = [term for term, record in acceptance["required_terms"].items() if not record["found"]]
        if missing:
            raise ValueError(
                f"preset acceptance failed for '{preset_id}': required gene symbols missing from the corpus "
                f"(as usable gene symbols): {', '.join(missing)}. Nothing was written, so the previous frozen "
                f"'{corpus_version}' corpus is untouched"
            )

    corpus_dir = out / safe_id(corpus_version)
    corpus_dir.mkdir(parents=True, exist_ok=True)
    corpus_path = corpus_dir / "reference_corpus.tsv"
    corpus.to_csv(corpus_path, sep="\t", index=False, lineterminator="\n")
    excluded_path = corpus_dir / "reference_corpus_excluded.tsv"
    pd.DataFrame(excluded, columns=EXCLUSION_COLUMNS).to_csv(excluded_path, sep="\t", index=False, lineterminator="\n")
    coverage = corpus_coverage(corpus)
    coverage_path = corpus_dir / "reference_corpus_coverage.tsv"
    coverage.to_csv(coverage_path, sep="\t", index=False, lineterminator="\n")

    tier_counts = {tier: int((corpus["provenance_tier"] == tier).sum()) for tier in USABLE_TIERS}
    manifest = {
        "schema_version": SCHEMA_VERSION, "corpus_format_version": CORPUS_FORMAT_VERSION,
        "candidate_qtl_version": __version__, "command": "reference-corpus-build",
        "corpus_version": corpus_version, "built_at_utc": _now(),
        "case_independent": True,
        "phenotype_used": False,
        "interval_arguments": [],
        "independence_rule": (
            "this command takes no phenotype, interval or candidate argument. A corpus built from a case "
            "phenotype would make Route C a second copy of Route A, and their agreement an artefact"
        ),
        "row_count": int(len(corpus)),
        "provenance_tier_counts": tier_counts,
        "model_proposed_rows": tier_counts.get("llm_verified", 0),
        "contains_model_proposed_rows": bool(tier_counts.get("llm_verified", 0)),
        "proposed_row_count": proposed,
        "excluded_row_count": len(excluded),
        "exclusion_reasons": {reason: sum(1 for row in excluded if row["exclusion_reason"] == reason)
                              for reason in sorted({str(row["exclusion_reason"]) for row in excluded})},
        "species": sorted({str(value) for value in corpus["reference_species"] if str(value)})[:20],
        "species_count": len({str(value) for value in corpus["reference_species"] if str(value).strip()}),
        "rows_per_species": {str(key): int(value) for key, value in
                             corpus["reference_species"].value_counts().sort_index().items() if str(key).strip()},
        # Row count says nothing about how much biology is here: a large corpus
        # with no mutant phenotype is exactly the failure this build addresses.
        "field_coverage": coverage_summary(coverage),
        "phenotype_row_count": int((corpus["mutant_phenotype"].astype(str).str.strip() != "").sum()),
        "pathway_row_count": int((corpus["pathway"].astype(str).str.strip() != "").sum()),
        "evidence_type_counts": {str(key): int(value) for key, value in
                                 corpus["evidence_type"].value_counts().sort_index().items() if str(key).strip()},
        "missing_phenotype_meaning": (
            "an empty mutant_phenotype is knowledge_unavailable for the routes downstream: the curated record "
            "carries no disruption phenotype. It is never evidence that the gene has no phenotype"
        ),
        "benchmark_mode": benchmark_mode,
        "sources": sources,
        "corpus_sha256": sha256(corpus_path),
        "excluded_sha256": sha256(excluded_path),
        "coverage_sha256": sha256(coverage_path),
        "promotion_gates": [
            "a retrieved row comes from a local authoritative source and involves no model",
            "a model-proposed row is promoted only when its accession appears in a retrieved source",
            "species and taxon of a promoted row are taken from the authoritative record, never from the proposal",
            "a model-proposed row needs at least one citation that resolves locally",
            "an unverified proposal is kept in the exclusion ledger and never written to the corpus",
            "in benchmark mode, focal-lineage / answer-shaped rows are barred",
        ],
        "downstream": {
            "route_c": f"retrieval-index-build --corpus {corpus_path} --database-version {corpus_version}",
            "route_b": f"route-b-prepare --reference-function {corpus_path}",
            "note": "one resource feeds both routes; it is versioned so a rebuild never silently changes a frozen run",
        },
        "outputs": {"corpus": str(corpus_path), "excluded": str(excluded_path),
                    "coverage": str(coverage_path),
                    "manifest": str(corpus_dir / "reference_corpus_manifest.json")},
    }
    if acceptance is not None:
        manifest["preset_acceptance"] = acceptance
    manifest_path = corpus_dir / "reference_corpus_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    result: dict[str, object] = {"manifest": str(manifest_path), "corpus": str(corpus_path), "excluded": str(excluded_path),
                                 "coverage": str(coverage_path),
                                 "corpus_version": corpus_version, "row_count": int(len(corpus)),
                                 "provenance_tier_counts": tier_counts, "excluded_row_count": len(excluded),
                                 "phenotype_row_count": manifest["phenotype_row_count"],
                                 "species_count": manifest["species_count"],
                                 "contains_model_proposed_rows": manifest["contains_model_proposed_rows"]}
    if acceptance is not None:
        result["preset_acceptance"] = acceptance
    return result
