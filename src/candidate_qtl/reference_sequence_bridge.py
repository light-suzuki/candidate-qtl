"""Post-freeze, source-bound sequence correspondence; never a Tier reviewer.

Only public accession identifiers can leave the workstation. Target proteins
are compared locally, all competitors remain visible, and partial annotation
models do not erase a useful locus correspondence.
"""
from __future__ import annotations

import io
import json
import math
import re
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from .blast_evidence import BLAST_OUTFMT, _resolve_runner, _tool_version, _wsl_path, build_protein_db
from .blast_panel import DEFAULT_THRESHOLDS
from .core import ExternalDependencyError, runtime_path, run_tool, sha256, write_fasta

BRIDGE_CONTRACT = "reference-sequence-bridge-1.0.0"
REQUEST_COLUMNS = ("request_id", "source_id", "reported_gene", "reference_accession", "reference_taxid", "focal_taxid")
CORRESPONDENCE_COLUMNS = (*REQUEST_COLUMNS, "cds_protein_id", "resolved_accession", "resolved_protein_id", "source_gene_names_json", "source_organism", "query_sequence_sha256",
    "sequence_state", "evaluation_state", "target_protein_id", "target_id", "pident", "query_coverage", "target_coverage",
    "query_length", "target_length", "query_start", "query_end", "target_start", "target_end", "query_segments", "target_segments",
    "evalue", "bitscore", "hsp_count", "alignment_scope", "eligibility", "locus_identity", "whole_annotation_model_identity", "rank_contribution", "note")
HIT_COLUMNS = ("panel_id", "target_id", "eligibility", "pident", "query_coverage", "target_coverage",
    "source_id", "reported_gene", "reference_accession", "cds_protein_id", "reference_taxid", "focal_taxid", "source_gene_names_json", "source_organism",
    "target_protein_id", "alignment_scope", "locus_identity", "whole_annotation_model_identity", "rank_contribution")
_ACCESSION = re.compile(r"^[A-Za-z][A-Za-z0-9_]*[0-9](?:\.[0-9]+)?$")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _fingerprint(path: Path) -> dict:
    return {"path": str(path.resolve()), "sha256": sha256(path)}


def _read_requests(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if set(frame.columns) - {*REQUEST_COLUMNS, "cds_protein_id"} or set(REQUEST_COLUMNS) - set(frame.columns):
        raise ValueError("sequence bridge request columns do not match the contract")
    if frame.empty or frame[list(REQUEST_COLUMNS)].eq("").any().any() or frame.request_id.duplicated().any():
        raise ValueError("sequence bridge requests require nonempty fields and unique request_id")
    if not frame.reference_accession.map(lambda item: bool(_ACCESSION.fullmatch(item))).all():
        raise ValueError("sequence bridge requires a public accession, never a sequence or URL")
    if not frame.reference_taxid.str.fullmatch(r"[1-9][0-9]*").all() or not frame.focal_taxid.str.fullmatch(r"[1-9][0-9]*").all():
        raise ValueError("sequence bridge requests require numeric reference/focal taxids")
    if "cds_protein_id" not in frame:
        frame["cds_protein_id"] = ""
    return frame


def _retrieve_record(accession: str, index: int, raw_dir: Path, records: Path | None,
                     allow_external_api: bool, timeout: int) -> tuple[object | None, dict]:
    """Cache-first official INSDC record; preserve failed bodies and provenance."""
    meta = {"reference_accession": accession, "retrieved_utc": _now(), "origin": "not_evaluated",
            "url": f"https://www.ebi.ac.uk/ena/browser/api/embl/{accession}?download=false", "state": "not_evaluated_record_unavailable"}
    text = None
    fmt = "embl"
    if records:
        for suffix, parser in ((".embl", "embl"), (".gb", "genbank"), (".gbk", "genbank")):
            cached = records / (accession + suffix)
            if cached.is_file():
                meta.update(origin="supplied_record_cache", input_record=_fingerprint(cached))
                try:
                    text, fmt = cached.read_text(encoding="utf-8"), parser
                except (OSError, UnicodeError) as error:
                    meta.update(state="not_evaluated_record_parse_failure", error=f"{type(error).__name__}: {error}")
                break
    if text is None and allow_external_api and "error" not in meta:
        meta["origin"] = "official_ena_live"
        try:
            request = urllib.request.Request(meta["url"], headers={"User-Agent": "candidate-qtl reference-sequence-bridge"})
            with urllib.request.urlopen(request, timeout=timeout) as response:
                text = response.read().decode("utf-8")
                meta["http_status"] = response.status
        except Exception as error:
            meta.update(state="not_evaluated_provider_failure", error=f"{type(error).__name__}: {error}")
            if hasattr(error, "read"):
                try:
                    text = error.read().decode("utf-8", errors="replace")
                except (OSError, UnicodeError):
                    pass
    saved = raw_dir / f"record_{index:06d}.txt"
    if text is not None:
        saved.write_text(text, encoding="utf-8")
        meta.update(raw_path=str(saved.resolve()), raw_sha256=sha256(saved))
    record = None
    if text is not None and "error" not in meta:
        try:
            record = SeqIO.read(io.StringIO(text), fmt)
            matches = record.id == accession if "." in accession else record.id.split(".")[0] == accession
            if not matches:
                meta.update(state="not_evaluated_accession_version_mismatch", actual_accession=record.id)
                record = None
            else:
                meta.update(state="record_parsed", resolved_accession=record.id)
        except Exception as error:
            meta.update(state="not_evaluated_record_parse_failure", error=f"{type(error).__name__}: {error}")
    (raw_dir / f"record_{index:06d}.meta.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return record, meta


def _extract_query(record, request: dict) -> tuple[str | None, dict]:
    taxids = {item.partition(":")[2] for feature in record.features if feature.type == "source"
              for item in feature.qualifiers.get("db_xref", []) if item.startswith("taxon:")}
    if taxids != {request["reference_taxid"]}:
        return None, {"sequence_state": "not_evaluated_reference_taxid_mismatch", "note": "official record taxid does not match request"}
    if request["reference_taxid"] != request["focal_taxid"]:
        return None, {"sequence_state": "not_evaluated_cross_species", "note": "same-species identity was not established"}
    features = [feature for feature in record.features if feature.type == "CDS"]
    selected = request.get("cds_protein_id", "")
    if selected:
        features = [feature for feature in features if selected in feature.qualifiers.get("protein_id", [])]
    if len(features) != 1:
        return None, {"sequence_state": "not_evaluated_ambiguous_cds", "note": "select an explicit source CDS protein_id; gene-name guesses are not accepted"}
    feature = features[0]
    translations = feature.qualifiers.get("translation", [])
    if len(translations) != 1:
        return None, {"sequence_state": "not_evaluated_missing_translation", "note": "no unique official CDS translation"}
    sequence = "".join(translations[0].split()).upper().removesuffix("*")
    if not sequence or "*" in sequence or not re.fullmatch("[A-Z]+", sequence):
        return None, {"sequence_state": "not_evaluated_invalid_translation", "note": "invalid official translation"}
    import hashlib
    source_tokens = sorted({token.strip() for field in ("gene", "locus_tag", "gene_synonym", "protein_id")
                            for value in feature.qualifiers.get(field, []) for token in value.split(";") if token.strip()})
    return sequence, {"sequence_state": "same_species_official_translation", "resolved_accession": record.id,
        "source_gene_names_json": json.dumps(source_tokens, ensure_ascii=False),
        "cds_protein_id": ";".join(feature.qualifiers.get("protein_id", [])),
        "resolved_protein_id": ";".join(feature.qualifiers.get("protein_id", [])),
        "query_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        "note": f"official CDS location={feature.location}; source genes={';'.join(feature.qualifiers.get('gene', []))}; no trait/Tier inference"}


def _base_row(request: dict) -> dict:
    return {**dict.fromkeys(CORRESPONDENCE_COLUMNS, ""), **request, "rank_contribution": "none",
        "locus_identity": "not_evaluated", "whole_annotation_model_identity": "not_evaluated", "eligibility": "not_evaluated", "source_gene_names_json": "[]"}


def _read_bridge_hsps(path: Path) -> list[dict]:
    """Keep every alignment separate, especially repeated domains in one model.

    A query may align both to its own locus and a homologous neighbouring gene
    inside a merged annotation. Combining the two HSPs' subject coordinates
    would manufacture full-model coverage and hide that annotation discrepancy.
    The original researcher can inspect all HSPs for genuinely split alignments.
    """
    hits = []
    with path.open(encoding="utf-8") as handle:
        for index, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            values = line.rstrip("\r\n").split("\t")
            if len(values) != 12:
                raise ValueError(f"sequence bridge raw HSP {index}: expected 12 columns")
            identity, length, qlen, slen = float(values[2]), int(values[3]), int(values[4]), int(values[5])
            qstart, qend = sorted((int(values[6]), int(values[7])))
            sstart, send = sorted((int(values[8]), int(values[9])))
            evalue, bitscore = float(values[10]), float(values[11])
            if (not all(math.isfinite(value) for value in (identity, evalue, bitscore))
                    or not 0 <= identity <= 100 or evalue < 0 or bitscore < 0 or length < 1
                    or not 1 <= qstart <= qend <= qlen or not 1 <= sstart <= send <= slen):
                raise ValueError(f"sequence bridge raw HSP {index}: invalid metrics or coordinates")
            hits.append({"qseqid": values[0], "sseqid": values[1], "pident": identity, "length": length,
                "qlen": qlen, "slen": slen, "qstart": qstart, "qend": qend, "sstart": sstart, "send": send,
                "evalue": evalue, "bitscore": bitscore, "hsp_count": 1, "alignment_scope": "single_hsp",
                "query_union_segments": f"{qstart}-{qend}", "target_union_segments": f"{sstart}-{send}",
                "query_coverage": (qend - qstart + 1) / qlen, "target_coverage": (send - sstart + 1) / slen})
    return sorted(hits, key=lambda row: (row["qseqid"], row["sseqid"], -row["bitscore"], row["qstart"], row["sstart"], row["qend"], row["send"]))


def _run_search(queries: Path, proteins: Path, out: Path, count: int, threads: int, runner: str) -> dict:
    selected, error = _resolve_runner(runner, ("blastp", "makeblastdb"))
    if error:
        return {"state": "not_evaluated_blast_unavailable", "error": error}
    db = out / "proteome_db"
    try:
        result = build_protein_db(str(proteins), str(db), threads, selected)
    except (OSError, ExternalDependencyError) as error:
        return {"state": "not_evaluated_blast_failure", "error": f"{type(error).__name__}: {error}"}
    if not result.get("available"):
        return {"state": "not_evaluated_blast_failure", "db": result}
    arg = lambda path: _wsl_path(path) if selected == "wsl" else str(path)
    raw = out / "all_hsps.tsv"
    command = (["wsl"] if selected == "wsl" else []) + ["blastp", "-query", arg(queries), "-db", arg(db),
        "-out", arg(raw), "-outfmt", BLAST_OUTFMT, "-evalue", str(DEFAULT_THRESHOLDS["max_evalue"]),
        "-num_threads", str(threads), "-max_target_seqs", str(max(1, count))]
    info = {"runner": selected, "command": command, "blastp_version": _tool_version("blastp", selected), "db": result}
    try:
        completed = run_tool(command, capture_output=True, text=True, check=False)
        info.update(returncode=completed.returncode, stdout=completed.stdout, stderr=completed.stderr)
        if completed.returncode:
            return {**info, "state": "not_evaluated_blast_failure"}
        return {**info, "state": "evaluated", "hits": _read_bridge_hsps(raw)}
    except (OSError, ExternalDependencyError, ValueError) as error:
        return {**info, "state": "not_evaluated_blast_failure", "error": f"{type(error).__name__}: {error}"}


def _model_conflicts(rows: list[dict]) -> None:
    """Distinct source genes matching disjoint model regions need model review."""
    for row in rows:
        if row["eligibility"] != "eligible":
            continue
        for other in rows:
            if other is row or other["eligibility"] != "eligible" or row["target_protein_id"] != other["target_protein_id"]:
                continue
            if row["reported_gene"] == other["reported_gene"] or row["resolved_protein_id"] == other["resolved_protein_id"]:
                continue
            if int(row["target_end"]) < int(other["target_start"]) or int(other["target_end"]) < int(row["target_start"]):
                row["whole_annotation_model_identity"] = "multiple_reported_genes_in_model"


def build_reference_sequence_bridge(requests_arg, proteome_arg, protein_gene_map_arg, out_dir_arg, *,
                                    records_arg=None, allow_external_api=False, threads=4, runner="auto", timeout=30) -> dict:
    """Resolve only source-selected CDSs, then compare the whole local proteome."""
    if threads < 1 or timeout < 1 or runner not in {"auto", "native", "wsl"}:
        raise ValueError("invalid sequence bridge runner, threads or timeout")
    requests_path = runtime_path(requests_arg)
    requests = _read_requests(requests_path)
    out = runtime_path(out_dir_arg).resolve()
    if out.exists():
        raise ValueError("sequence bridge requires a new output generation")
    out.mkdir(parents=True)
    raw_dir = out / "raw"
    raw_dir.mkdir()
    snapshot = out / "requests.tsv"
    snapshot.write_bytes(requests_path.read_bytes())
    freeze = {"contract": BRIDGE_CONTRACT, "created_utc": _now(), "requests_sha256": sha256(snapshot),
              "target_accessed": False, "rank_contribution": "none"}
    (out / "request_freeze.json").write_text(json.dumps(freeze, indent=2) + "\n", encoding="utf-8")
    cache = runtime_path(records_arg) if records_arg else None
    rows, queries, provenance = [], {}, []
    query_keys = {}
    for index, request in enumerate(requests.to_dict("records"), 1):
        row = _base_row(request)
        record, meta = _retrieve_record(request["reference_accession"], index, raw_dir, cache, allow_external_api, timeout)
        provenance.append(meta)
        row["sequence_state"] = meta["state"]
        if record is not None:
            row["source_organism"] = str(record.annotations.get("organism", ""))
            sequence, detail = _extract_query(record, request)
            row.update(detail)
            if sequence:
                query_id = f"q{index:06d}"
                queries[query_id] = sequence
                query_keys[query_id] = row
        rows.append(row)
    # Source-selected queries are recorded before any target sequence is opened.
    query_path = out / "queries.fa"
    write_fasta(query_path, list(queries.items()))
    source_freeze = {"request_freeze_sha256": sha256(out / "request_freeze.json"), "query_sha256": sha256(query_path),
                     "record_sha256": {meta["raw_path"]: meta["raw_sha256"] for meta in provenance if "raw_path" in meta},
                     "target_accessed": False}
    (out / "source_freeze.json").write_text(json.dumps(source_freeze, indent=2) + "\n", encoding="utf-8")
    proteins_path, map_path = runtime_path(proteome_arg), runtime_path(protein_gene_map_arg)
    mapping = pd.read_csv(map_path, sep="\t", dtype=str, keep_default_na=False)
    if not {"protein_id", "gene_id"} <= set(mapping) or mapping[["protein_id", "gene_id"]].eq("").any().any() or mapping.protein_id.duplicated().any():
        raise ValueError("sequence bridge protein-gene map requires unique protein_id and nonempty gene_id")
    gene_map = dict(zip(mapping.protein_id, mapping.gene_id))
    proteins = {}
    with proteins_path.open(encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in proteins or record.id not in gene_map:
                raise ValueError("sequence bridge proteome has duplicate or unmapped protein IDs")
            sequence = str(record.seq).upper().removesuffix("*")
            if not sequence or "*" in sequence or not re.fullmatch("[A-Z]+", sequence):
                raise ValueError("sequence bridge proteome requires valid proteins; terminal stop only")
            proteins[record.id] = sequence
    if not proteins:
        raise ValueError("sequence bridge proteome is empty")
    normalized = out / "target_proteins.fa"
    write_fasta(normalized, list(proteins.items()))
    search = _run_search(query_path, normalized, out, len(proteins), threads, runner) if queries else {"state": "not_evaluated_no_queries", "hits": []}
    raw_hits = search.pop("hits", [])
    all_hit_rows = []
    for hit in raw_hits:
        if hit["qseqid"] not in query_keys or hit["sseqid"] not in proteins:
            raise ValueError("sequence bridge BLAST returned unknown query or target IDs")
        base = query_keys[hit["qseqid"]]
        if hit["qlen"] != len(queries[hit["qseqid"]]) or hit["slen"] != len(proteins[hit["sseqid"]]):
            raise ValueError("sequence bridge BLAST sequence length mismatch")
        eligible = hit["qlen"] >= DEFAULT_THRESHOLDS["min_query_length"] and hit["pident"] >= DEFAULT_THRESHOLDS["min_identity"] and hit["query_coverage"] >= DEFAULT_THRESHOLDS["min_query_coverage"] and hit["evalue"] <= DEFAULT_THRESHOLDS["max_evalue"]
        all_hit_rows.append({**base, "evaluation_state": "evaluated_with_hit", "target_protein_id": hit["sseqid"], "target_id": gene_map[hit["sseqid"]],
            "pident": hit["pident"], "query_coverage": hit["query_coverage"], "target_coverage": hit["target_coverage"],
            "query_length": hit["qlen"], "target_length": hit["slen"], "query_start": hit["qstart"], "query_end": hit["qend"],
            "target_start": hit["sstart"], "target_end": hit["send"], "query_segments": hit["query_union_segments"], "target_segments": hit["target_union_segments"],
            "evalue": hit["evalue"], "bitscore": hit["bitscore"], "hsp_count": hit["hsp_count"], "alignment_scope": "single_hsp", "eligibility": "eligible" if eligible else "below_review_threshold",
            "locus_identity": "requires_researcher_assessment",
            "whole_annotation_model_identity": "model_coverage_passes_review_threshold" if hit["target_coverage"] >= DEFAULT_THRESHOLDS["min_target_coverage"] else "partial_model_correspondence"})
    seen = {row["request_id"] for row in all_hit_rows}
    for row in rows:
        if row["request_id"] in seen:
            continue
        if row["sequence_state"] == "same_species_official_translation":
            row["evaluation_state"] = "evaluated_no_hit" if search["state"] == "evaluated" else search["state"]
            row["eligibility"] = "no_hit" if search["state"] == "evaluated" else "not_evaluated"
        else:
            row["evaluation_state"] = "not_evaluated"
        all_hit_rows.append(row)
    _model_conflicts(all_hit_rows)
    correspondence = pd.DataFrame(all_hit_rows, columns=CORRESPONDENCE_COLUMNS).fillna("")
    correspondence = correspondence.sort_values(["request_id", "target_protein_id"], kind="stable")
    correspondence_path = out / "reference_sequence_correspondence.tsv"
    correspondence.to_csv(correspondence_path, sep="\t", index=False, lineterminator="\n")
    hit_rows = [{"panel_id": row["request_id"], **{key: row[key] for key in HIT_COLUMNS if key != "panel_id"}}
                for row in all_hit_rows if row["eligibility"] == "eligible"]
    hits_path = out / "eligible_hits.tsv"
    pd.DataFrame(hit_rows, columns=HIT_COLUMNS).sort_values(["panel_id", "target_id", "target_protein_id"], kind="stable").to_csv(hits_path, sep="\t", index=False, lineterminator="\n")
    inputs = {"requests": _fingerprint(requests_path), "proteome": _fingerprint(proteins_path), "protein_gene_map": _fingerprint(map_path)}
    for index, meta in enumerate(provenance):
        if "input_record" in meta:
            inputs[f"record_{index:06d}"] = meta["input_record"]
    manifest = {"contract": BRIDGE_CONTRACT, "created_utc": _now(), "inputs": inputs, "retrieval": provenance,
        "search": search, "thresholds": DEFAULT_THRESHOLDS, "source_selected_before_target_access": True,
        "external_transmission": "public_accession_identifiers_only; target proteins never uploaded",
        "tier_effect": "none", "support_effect": "none", "rank_contribution": "none",
        "outputs": {"requests_snapshot": str(snapshot), "correspondence": str(correspondence_path), "eligible_hits": str(hits_path)},
        "artifacts": {path.relative_to(out).as_posix(): sha256(path) for path in sorted(out.rglob("*")) if path.is_file()}}
    manifest_path = out / "reference_sequence_bridge_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), **manifest["outputs"], "request_count": len(requests), "eligible_hit_count": len(hit_rows),
            "evaluation_state": search["state"], "tier_effect": "none", "support_effect": "none"}


def load_reference_sequence_bridge(manifest_arg) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    """Verify all retained artifacts and original inputs before any consumption."""
    path = runtime_path(manifest_arg).resolve()
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("contract") != BRIDGE_CONTRACT or manifest.get("tier_effect") != "none" or manifest.get("support_effect") != "none":
        raise ValueError("unsupported sequence bridge contract or promotion")
    for name, digest in manifest.get("artifacts", {}).items():
        artifact = (path.parent / name).resolve()
        if not artifact.is_relative_to(path.parent) or not artifact.is_file() or sha256(artifact) != digest:
            raise ValueError("sequence bridge artifact hash mismatch")
    if not manifest.get("artifacts") or not manifest.get("inputs"):
        raise ValueError("sequence bridge requires complete artifacts and inputs")
    for descriptor in manifest["inputs"].values():
        original = runtime_path(descriptor["path"])
        if not original.is_file() or sha256(original) != descriptor["sha256"]:
            raise ValueError("sequence bridge input hash mismatch")
    outputs = {}
    for key in ("requests_snapshot", "correspondence", "eligible_hits"):
        output = runtime_path(manifest["outputs"][key]).resolve()
        if not output.is_relative_to(path.parent) or output.relative_to(path.parent).as_posix() not in manifest["artifacts"]:
            raise ValueError("sequence bridge output is not sealed")
        outputs[key] = output
    requests = _read_requests(outputs["requests_snapshot"])
    request_freeze = json.loads((path.parent / "request_freeze.json").read_text(encoding="utf-8"))
    source_freeze = json.loads((path.parent / "source_freeze.json").read_text(encoding="utf-8"))
    if (request_freeze.get("requests_sha256") != sha256(outputs["requests_snapshot"])
            or manifest["inputs"]["requests"]["sha256"] != sha256(outputs["requests_snapshot"])
            or source_freeze.get("request_freeze_sha256") != sha256(path.parent / "request_freeze.json")
            or source_freeze.get("query_sha256") != sha256(path.parent / "queries.fa")
            or request_freeze.get("target_accessed") is not False or source_freeze.get("target_accessed") is not False):
        raise ValueError("sequence bridge freeze chain mismatch")
    frame = pd.read_csv(outputs["correspondence"], sep="\t", dtype=str, keep_default_na=False)
    hits = pd.read_csv(outputs["eligible_hits"], sep="\t", dtype=str, keep_default_na=False)
    if tuple(frame.columns) != CORRESPONDENCE_COLUMNS or tuple(hits.columns) != HIT_COLUMNS:
        raise ValueError("sequence bridge output columns do not match contract")
    if set(frame.request_id) != set(requests.request_id) or not set(hits.panel_id) <= set(requests.request_id):
        raise ValueError("sequence bridge request roster mismatch")
    by_request = requests.set_index("request_id").to_dict("index")
    for row in frame.to_dict("records"):
        request = by_request[row["request_id"]]
        if any(row[key] != request[key] for key in REQUEST_COLUMNS if key != "request_id"):
            raise ValueError("sequence bridge correspondence source endpoint mismatch")
        if request.get("cds_protein_id") and row["cds_protein_id"] != request["cds_protein_id"]:
            raise ValueError("sequence bridge selected CDS endpoint mismatch")
    expected_hits = [{"panel_id": row["request_id"], **{key: row[key] for key in HIT_COLUMNS if key != "panel_id"}}
                     for row in frame.to_dict("records") if row["eligibility"] == "eligible"]
    canonical = lambda rows: sorted(json.dumps(row, sort_keys=True) for row in rows)
    if canonical(expected_hits) != canonical(hits.to_dict("records")):
        raise ValueError("sequence bridge eligible hits contradict correspondence")
    return manifest, frame, hits
