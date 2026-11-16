"""Provider-curated bait sequences for frozen Direct-C records (review only).

A Direct-C record names a gene in a cited paper. BLAST projection needs a
protein sequence for that gene, and the frozen record usually carries only a
symbol. This module selects the sequence from one structured provider link:
a UniProtKB entry of the record's reference species that lists the record's
cited PubMed ID among its references and whose structured gene names include
the record's reported symbol. Free text, top hits and the interval are never
consulted. Zero or several qualifying entries fail closed.

The resulting edge says "the cited paper's gene symbol and this provider
sequence are the same curated entry". It authorizes review-only projection
provenance. It never changes support, Tier or rank.
"""
from __future__ import annotations

import json
import re
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from .blast_panel import PANEL_COLUMNS, _row_sha256
from .core import runtime_path, sha256

CONTRACT = "route-c-provider-bait-1.0.0"
EDGE_STATE = "verified_provider_sequence_reference"
EDGE_COLUMNS = (
    "record_id", "source_identity_field", "source_identity_value", "reference_species",
    "citation", "provider", "provider_accession", "provider_record_path",
    "provider_record_sha256", "matched_gene_name", "reference_positions",
    "query_sequence_sha256", "rank_contribution",
)
STATUS_COLUMNS = ("record_id", "state", "provider_accessions", "note")
_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search?query={query}&format=json&size=100"
_SPLIT = re.compile(r"[\s/;,()]+")


def _now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _species_prefix(species: str) -> str:
    parts = str(species or "").split()
    return (parts[0][:1] + parts[1][:1]).casefold() if len(parts) >= 2 else ""


def _name_keys(values: list[str], species: str) -> set[str]:
    """Casefolded name tokens, plus the genus/species-initial stripped form.

    "PsNIN" for Pisum sativum also yields "nin"; the prefix is derived from
    the record's own reference species, never guessed from the token.
    """
    prefix = _species_prefix(species)
    keys: set[str] = set()
    for value in values:
        for token in _SPLIT.split(str(value or "")):
            token = token.strip().casefold().replace("-", "")
            if len(token) < 2:
                continue
            keys.add(token)
            if prefix and token.startswith(prefix) and len(token) > len(prefix) + 1:
                keys.add(token[len(prefix):])
    return keys


def _entry_names(entry: dict) -> list[str]:
    names: list[str] = []
    for gene in entry.get("genes") or []:
        if not isinstance(gene, dict):
            continue
        if isinstance(gene.get("geneName"), dict):
            names.append(str(gene["geneName"].get("value", "")))
        for key in ("synonyms", "orderedLocusNames", "orfNames"):
            names.extend(str(item.get("value", "")) for item in gene.get(key) or []
                         if isinstance(item, dict))
    return [name for name in names if name.strip()]


def _record_pmids(record: dict) -> list[str]:
    from .citation_tokens import split_citations  # lazy: keep the core parser import-light

    pmids = []
    for token in split_citations(str(record.get("source_id", ""))):
        token = token.strip()
        if token.upper().startswith("PMID:") and token[5:].strip().isdigit():
            pmids.append(token[5:].strip())
    pmid = str(record.get("pmid", "")).strip()
    if pmid.isdigit():
        pmids.append(pmid)
    return sorted(set(pmids))


def entry_reference_positions(entry: dict, pmid: str) -> list[str]:
    positions: list[str] = []
    for reference in entry.get("references") or []:
        citation = reference.get("citation") if isinstance(reference, dict) else None
        if not isinstance(citation, dict):
            continue
        ids = {str(item.get("id", "")) for item in citation.get("citationCrossReferences") or []
               if isinstance(item, dict) and item.get("database") == "PubMed"}
        if pmid in ids:
            positions.extend(str(item) for item in reference.get("referencePositions") or [])
    return positions


def entry_matches(entry: dict, record: dict, pmid: str) -> str:
    """Return the matched provider gene name, or "" when the entry does not qualify."""
    species = str(record.get("reference_species", "")).strip()
    organism = entry.get("organism") if isinstance(entry.get("organism"), dict) else {}
    if not species or str(organism.get("scientificName", "")).casefold() != species.casefold():
        return ""
    if not entry_reference_positions(entry, pmid):
        return ""
    sequence = (entry.get("sequence") or {}).get("value") if isinstance(entry.get("sequence"), dict) else ""
    if not str(sequence or "").strip():
        return ""
    wanted = _name_keys([str(record.get("reported_gene_symbol", ""))], species)
    for name in _entry_names(entry):
        if _name_keys([name], species) & wanted:
            return name
    return ""


ALLELIC_MIN_RATIO = 0.95


def _molecule_types(entry: dict) -> set[str]:
    types: set[str] = set()
    for xref in entry.get("uniProtKBCrossReferences") or []:
        if isinstance(xref, dict) and xref.get("database") == "EMBL":
            types.update(str(item.get("value", "")) for item in xref.get("properties") or []
                         if isinstance(item, dict) and item.get("key") == "MoleculeType")
    return types


def allelic_representative(matches: dict[str, tuple[dict, str, str]]) -> str | None:
    """Pick one entry when every match is the same named gene with near-identical sequence.

    A paper that sequenced wild-type and mutant alleles yields several entries
    of one gene. Divergent sequences (different genes or paralogs) return None
    and the caller fails closed. The choice never looks at the interval:
    reviewed entry, then mRNA evidence, then the longest sequence, then accession.
    """
    from difflib import SequenceMatcher

    if len({name.casefold() for _, name, _ in matches.values()}) != 1:
        return None
    sequences = {acc: str(entry["sequence"]["value"]).upper() for acc, (entry, _, _) in matches.items()}
    items = sorted(sequences.items())
    for index, (_, left) in enumerate(items):
        for _, right in items[index + 1:]:
            if SequenceMatcher(None, left, right, autojunk=False).ratio() < ALLELIC_MIN_RATIO:
                return None

    def key(accession: str) -> tuple:
        entry = matches[accession][0]
        reviewed = "reviewed" in str(entry.get("entryType", "")).casefold() and "unreviewed" not in str(
            entry.get("entryType", "")).casefold()
        return (not reviewed, "mRNA" not in _molecule_types(entry), -len(sequences[accession]), accession)

    return min(matches, key=key)


def _search(pmid: str, cache: Path, allow_external_api: bool, timeout: int) -> list[dict]:
    path = cache / f"uniprot_lit_{pmid}.json"
    if not path.is_file():
        if not allow_external_api:
            raise FileNotFoundError(f"provider search cache missing for PMID {pmid}: {path}")
        url = _SEARCH_URL.format(query=urllib.parse.quote(f"lit_citation_id:{pmid}"))
        request = urllib.request.Request(url, headers={"User-Agent": "candidate-qtl route-c-provider-bait"})
        with urllib.request.urlopen(request, timeout=timeout) as response:  # noqa: S310 - fixed https host
            payload = json.loads(response.read().decode("utf-8"))
        cache.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    payload = json.loads(path.read_text(encoding="utf-8"))
    results = payload.get("results") if isinstance(payload, dict) else None
    if not isinstance(results, list):
        raise ValueError(f"provider search cache for PMID {pmid} lacks results")
    return [item for item in results if isinstance(item, dict)]


def build_provider_bait(frozen_arg: str, freeze_manifest_arg: str, out_dir_arg: str, *,
                        cache_arg: str | None = None, allow_external_api: bool = False,
                        timeout: int = 60) -> dict[str, object]:
    from .route_c_direct_rescue import validate_direct_rescue_frame

    frozen_path, freeze_path = runtime_path(frozen_arg), runtime_path(freeze_manifest_arg)
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    if freeze.get("command") != "route-c-direct-freeze":
        raise ValueError("provider bait requires a route-c-direct-freeze manifest")
    if str(freeze.get("frozen_sha256", "")).removeprefix("sha256:") != sha256(frozen_path):
        raise ValueError("Route C frozen artifact hash does not match freeze manifest")
    direct = validate_direct_rescue_frame(
        pd.read_csv(frozen_path, sep="\t", dtype=str, keep_default_na=False).fillna(""))
    out = runtime_path(out_dir_arg)
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"provider bait output directory is not empty: {out}")
    records_dir = out / "provider_records"
    records_dir.mkdir(parents=True, exist_ok=True)
    cache = runtime_path(cache_arg) if cache_arg else out / "provider_search_cache"
    retrieved_at = _now()
    edges, panel, statuses = [], [], []
    for record in sorted(direct.to_dict("records"), key=lambda row: str(row["record_id"])):
        record_id = str(record["record_id"]).strip()
        symbol = str(record.get("reported_gene_symbol", "")).strip()
        if record.get("direct_rescue_eligibility") != "eligible_direct_record":
            statuses.append({"record_id": record_id, "state": "not_evaluated_ineligible_record",
                             "provider_accessions": "", "note": str(record.get("direct_rescue_eligibility", ""))})
            continue
        pmids = _record_pmids(record)
        if not symbol or not pmids:
            statuses.append({"record_id": record_id, "state": "not_evaluated_missing_symbol_or_pmid",
                             "provider_accessions": "", "note": ""})
            continue
        matches: dict[str, tuple[dict, str, str]] = {}
        for pmid in pmids:
            for entry in _search(pmid, cache, allow_external_api, timeout):
                name = entry_matches(entry, record, pmid)
                accession = str(entry.get("primaryAccession", "")).strip()
                if name and accession:
                    matches.setdefault(accession, (entry, name, pmid))
        representative = allelic_representative(matches) if len(matches) > 1 else None
        if len(matches) != 1 and representative is None:
            statuses.append({"record_id": record_id,
                             "state": "no_provider_entry" if not matches else "ambiguous_provider_entries",
                             "provider_accessions": ";".join(sorted(matches)),
                             "note": "fail closed; no bait is selected"})
            continue
        accession = representative or next(iter(matches))
        entry, name, pmid = matches[accession]
        entry_path = records_dir / f"{accession}.json"
        entry_path.write_text(json.dumps(entry, ensure_ascii=False, indent=1, sort_keys=True) + "\n",
                              encoding="utf-8")
        sequence = str(entry["sequence"]["value"]).strip().upper()
        sequence_sha = _row_sha256(sequence)
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
        edges.append({
            "record_id": record_id, "source_identity_field": "reported_gene_symbol",
            "source_identity_value": symbol, "reference_species": str(record["reference_species"]),
            "citation": f"PMID:{pmid}", "provider": "UniProtKB", "provider_accession": accession,
            "provider_record_path": str(entry_path), "provider_record_sha256": sha256(entry_path),
            "matched_gene_name": name,
            "reference_positions": ";".join(entry_reference_positions(entry, pmid)),
            "query_sequence_sha256": sequence_sha, "rank_contribution": "none",
        })
        row = dict.fromkeys(PANEL_COLUMNS, "")
        row.update(panel_id=record_id, mechanism=f"record-level Direct C bait (frozen record {record_id})",
                   synonyms=";".join(sorted({symbol, *_entry_names(entry)})), reference_gene_name=symbol,
                   protein_accession=accession, taxon=str(record["reference_species"]), source="UniProtKB",
                   source_url=url, retrieved_at=retrieved_at, verification_status="verified",
                   verification_accession=accession, verification_sequence_sha256=sequence_sha,
                   verification_source_url=url, verification_retrieved_at=retrieved_at,
                   query_sequence=sequence, query_sequence_sha256=sequence_sha,
                   provider="route-c-provider-bait", model="none", rank_contribution="none")
        panel.append(row)
        statuses.append({"record_id": record_id,
                         "state": "provider_bait_selected_allelic_set" if representative else "provider_bait_selected",
                         "provider_accessions": ";".join(sorted(matches)),
                         "note": f"PMID:{pmid}; gene name {name}; bait {accession}"
                                 + ("; representative of near-identical entries" if representative else "")})
    edges_path, panel_path, status_path = (out / "provider_identity_edges.tsv", out / "panel.tsv",
                                           out / "provider_bait_status.tsv")
    pd.DataFrame(edges, columns=EDGE_COLUMNS).to_csv(edges_path, sep="\t", index=False, lineterminator="\n")
    pd.DataFrame(panel, columns=PANEL_COLUMNS).to_csv(panel_path, sep="\t", index=False, lineterminator="\n")
    pd.DataFrame(statuses, columns=STATUS_COLUMNS).to_csv(status_path, sep="\t", index=False, lineterminator="\n")
    manifest = {
        "contract": CONTRACT, "command": "route-c-provider-bait", "retrieved_at": retrieved_at,
        "interval_accessed": False, "rank_contribution": "none",
        "policy": [
            "one UniProtKB entry of the record's reference species that cites the record's PMID and lists its symbol",
            "zero or several qualifying entries fail closed; no top hit, alias guess or interval knowledge",
            "the edge is review-only projection provenance and never trait support",
        ],
        "inputs": {"frozen": {"path": str(frozen_path), "sha256": sha256(frozen_path)},
                   "freeze_manifest": {"path": str(freeze_path), "sha256": sha256(freeze_path)}},
        "outputs": {"edges": str(edges_path), "panel": str(panel_path), "status": str(status_path)},
        "edges_sha256": sha256(edges_path), "panel_sha256": sha256(panel_path),
        "selected_count": len(edges), "record_count": len(statuses),
    }
    manifest_path = out / "route_c_provider_bait_manifest.json"
    manifest["outputs"]["manifest"] = str(manifest_path)
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
                             encoding="utf-8")
    return {"manifest": str(manifest_path), "edges": str(edges_path), "panel": str(panel_path),
            "status": str(status_path), "selected_count": len(edges), "record_count": len(statuses)}


def _search_chain_rows(eligible_hits: Path) -> dict[str, dict]:
    """Return frozen bait-panel rows behind a blast-panel-search hits artifact."""
    search_path = eligible_hits.parent / "blast_panel_search_manifest.json"
    if not search_path.is_file():
        raise ValueError("provider identity edges require the blast-panel-search manifest beside the hits")
    search = json.loads(search_path.read_text(encoding="utf-8"))
    if (search.get("command") != "blast-panel-search"
            or search.get("eligible_interval_hits_sha256") != sha256(eligible_hits)):
        raise ValueError("provider identity edges: eligible hits are not the sealed search output")
    freeze_path = runtime_path(str(search.get("freeze", "")))
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    panel_path = runtime_path(str((freeze.get("outputs") or {}).get("panel", "")))
    if (freeze.get("command") != "blast-panel-freeze" or not panel_path.is_file()
            or sha256(panel_path) != freeze.get("panel_sha256")):
        raise ValueError("provider identity edges: bait panel freeze hash mismatch")
    rows = freeze.get("rows")
    if not isinstance(rows, dict):
        raise ValueError("provider identity edges: bait panel freeze lacks rows")
    return rows


def load_provider_edges(frame: pd.DataFrame, path: Path, records: dict[str, dict],
                        eligible_hits: Path | None) -> pd.DataFrame:
    """Recheck every edge against provider bytes, the frozen record and the BLAST chain."""
    frame = frame[list(EDGE_COLUMNS)].apply(lambda col: col.astype(str).str.strip())
    if frame.empty or frame.eq("").any().any():
        raise ValueError("provider identity edges require non-empty fields")
    if frame.record_id.duplicated().any():
        raise ValueError("provider identity edges allow one bait per record")
    if eligible_hits is None:
        raise ValueError("provider identity edges require the sealed eligible hits they explain")
    chain = _search_chain_rows(eligible_hits)
    for row in frame.to_dict("records"):
        record = records.get(row["record_id"])
        if record is None:
            raise ValueError("provider identity edge unknown record")
        if (row["source_identity_field"] != "reported_gene_symbol"
                or row["source_identity_value"] != str(record.get("reported_gene_symbol", "")).strip()
                or row["reference_species"] != str(record.get("reference_species", "")).strip()
                or row["provider"] != "UniProtKB" or row["rank_contribution"] != "none"):
            raise ValueError("provider identity edge does not match the frozen record")
        if row["citation"] not in {f"PMID:{pmid}" for pmid in _record_pmids(record)}:
            raise ValueError("provider identity edge citation is not the record's source")
        source = runtime_path(row["provider_record_path"])
        if not source.is_absolute():
            source = path.parent / source
        if not source.is_file() or sha256(source) != row["provider_record_sha256"]:
            raise ValueError("provider identity edge record hash mismatch")
        entry = json.loads(source.read_text(encoding="utf-8"))
        if str(entry.get("primaryAccession", "")) != row["provider_accession"]:
            raise ValueError("provider identity edge accession mismatch")
        if not entry_matches(entry, record, row["citation"].removeprefix("PMID:")):
            raise ValueError("provider identity edge no longer satisfies the provider rule")
        if _row_sha256(str(entry["sequence"]["value"]).strip()) != row["query_sequence_sha256"]:
            raise ValueError("provider identity edge sequence hash mismatch")
        bait = chain.get(row["record_id"])
        if (not isinstance(bait, dict) or bait.get("verification_accession") != row["provider_accession"]
                or bait.get("query_sequence_sha256") != row["query_sequence_sha256"]):
            raise ValueError("provider identity edge is not the frozen BLAST bait for this record")
    return frame.sort_values("record_id", kind="stable").reset_index(drop=True)


def provider_edge_links(frame: pd.DataFrame, record: dict, identity_type: str, identity: str,
                        observation: dict) -> bool:
    """Whether a source-bound identity reaches this BLAST hit through the provider edge.

    A symbol must be the edge's own endpoint. A locus must be the record's
    frozen locus and share a name with the provider entry's structured gene name.
    """
    if "provider_accession" not in frame:
        return False
    if (observation.get("mapping_method") != "blastp"
            or observation.get("source_identity_verification_state") != EDGE_STATE
            or observation.get("source_identity_field") != "reported_gene_symbol"):
        return False
    edges = frame[frame["record_id"] == str(record.get("record_id", ""))]
    if edges.empty:
        return False
    edge = edges.iloc[0]
    if observation.get("source_identity_value") != edge["source_identity_value"]:
        return False
    if identity_type == "gene_symbol":
        return identity == edge["source_identity_value"]
    if identity_type == "locus":
        species = str(edge["reference_species"])
        return (identity == str(record.get("reported_locus", "")).strip()
                and bool(_name_keys([identity], species) & _name_keys([str(edge["matched_gene_name"])], species)))
    return False
