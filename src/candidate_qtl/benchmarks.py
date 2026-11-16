from __future__ import annotations

import json
import hashlib
import re
import shutil
from pathlib import Path
from time import perf_counter
from typing import Iterable

import pandas as pd
from Bio import SeqIO

from .benchmark_metrics import aggregate_resolution_metrics, resolution_metrics
from .benchmark_reporting import wilson_interval
from .core import attributes, file_fingerprint, open_text, runtime_path, run_tool, safe_id, subset_proteins
from .evidence_io import EVIDENCE_COLUMNS, candidate_frame as _candidate_frame, id_map as _id_map, evidence_row as _row
from .evidence_policy import evidence_build
from .hardening import prioritize_evidence
from .optional_deps import require_pysam
from .reference_extract import extract_candidates


DEFAULT_ANNOTATION_SOURCE = "doi:10.57745/5MDE99#CAMEOR_V2_ANNOTATION_1_prot_eggnog.emapper.annotations"


def _record_aliases(record) -> set[str]:
    aliases = {record.id, record.name}
    aliases.update(token for token in record.id.split("|") if token)
    match = re.search(r"(?:^|\s)GN=([^\s]+)", record.description)
    if match:
        aliases.add(match.group(1))
    return {alias for alias in aliases if alias}


def _homolog_sequence_index(paths: list[str]) -> dict[str, tuple[object, str]]:
    index: dict[str, tuple[object, str]] = {}
    for path_arg in paths:
        path = runtime_path(path_arg)
        with open_text(path) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                for alias in _record_aliases(record):
                    index.setdefault(alias, (record, str(path)))
    return index


def _write_homolog_queries(
    query_ids: list[str],
    sequence_index: dict[str, tuple[object, str]],
    output: Path,
    audited_fingerprints: dict[str, str] | None = None,
) -> dict[str, object]:
    records = []
    sources: list[dict[str, str]] = []
    missing: list[str] = []
    for number, query_id in enumerate(query_ids, 1):
        match = sequence_index.get(query_id)
        if match is None:
            missing.append(query_id)
            continue
        source_record, source_path = match
        record = source_record[:]
        record.id = f"query{number}_{safe_id(query_id)}"
        record.name = record.id
        record.description = f"requested_id={query_id} source_id={source_record.id}"
        seq_sha256 = hashlib.sha256(str(record.seq).upper().encode("ascii", "ignore")).hexdigest()
        if audited_fingerprints is not None:
            expected_sha256 = audited_fingerprints.get(query_id) or audited_fingerprints.get(source_record.id)
            if not expected_sha256:
                for alias in _record_aliases(source_record):
                    if alias in audited_fingerprints:
                        expected_sha256 = audited_fingerprints[alias]
                        break
            if expected_sha256 is None:
                raise ValueError(f"homolog query ID {query_id} (source: {source_record.id}) was not audited in query eligibility manifest")
            if seq_sha256 != expected_sha256:
                raise ValueError(f"homolog query {query_id} sequence SHA-256 fingerprint mismatch: computed {seq_sha256} != audited {expected_sha256}")
        records.append(record)
        sources.append(
            {
                "query_id": record.id,
                "requested_id": query_id,
                "source_record_id": source_record.id,
                "source_path": source_path,
                "length": str(len(record.seq)),
                "sha256": seq_sha256,
            }
        )
    if missing:
        raise ValueError("homolog query IDs were not found in --homolog-proteins: " + ", ".join(missing))
    output.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, str(output), "fasta")
    return {"output": str(output), "queries": sources, "query_count": len(records)}


def _wsl_path(path: Path) -> str:
    resolved = path.resolve()
    drive = resolved.drive.rstrip(":").casefold()
    if drive:
        suffix = resolved.as_posix().split(":", 1)[1]
        return f"/mnt/{drive}{suffix}"
    return resolved.as_posix()


def _run_interval_blastp(
    query_path: Path,
    subject_path: Path,
    output: Path,
    runner: str,
    num_threads: int,
) -> dict[str, object]:
    if runner not in {"auto", "native", "wsl"}:
        raise ValueError("blastp runner must be auto, native, or wsl")
    selected = runner
    if runner == "auto":
        selected = "native" if shutil.which("blastp") else "wsl"
    if selected == "native" and not shutil.which("blastp"):
        raise RuntimeError("blastp was not found on PATH")
    if selected == "wsl" and not shutil.which("wsl"):
        raise RuntimeError("WSL was not found; install BLAST+ natively or run with WSL")

    query_lengths = [len(record.seq) for record in SeqIO.parse(str(query_path), "fasta")]
    short_query = bool(query_lengths and min(query_lengths) < 100)
    blast_args = [
        "blastp",
        "-query",
        _wsl_path(query_path) if selected == "wsl" else str(query_path),
        "-subject",
        _wsl_path(subject_path) if selected == "wsl" else str(subject_path),
        "-outfmt",
        "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore",
        "-max_target_seqs",
        "5000",
        "-num_threads",
        str(num_threads),
        "-out",
        _wsl_path(output) if selected == "wsl" else str(output),
    ]
    if short_query:
        blast_args.extend(
            ["-word_size", "2", "-matrix", "PAM30", "-gapopen", "9", "-gapextend", "1", "-comp_based_stats", "0", "-seg", "no", "-evalue", "1000"]
        )
    command = (["wsl"] + blast_args) if selected == "wsl" else blast_args
    output.parent.mkdir(parents=True, exist_ok=True)
    completed = run_tool(command, capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"interval BLASTP failed ({selected}): {completed.stderr.strip()}")
    return {
        "output": str(output),
        "runner": selected,
        "short_query_mode": short_query,
        "command": command,
        "stderr": completed.stderr.strip(),
    }


def homolog_blast_evidence(
    candidates_arg: str,
    blast_arg: str,
    out_arg: str,
    source_id: str,
    min_relative_bitscore: float = 0.75,
    evidence_code: str = "curated_homolog_match",
    query_source_map: dict[str, str] | None = None,
    query_confidence_map: dict[str, str] | None = None,
    max_confidence: str = "high",
) -> dict[str, object]:
    """Convert interval-only homolog BLASTP results into tiered, auditable evidence."""
    candidates = _candidate_frame(candidates_arg)
    valid_genes = set(candidates["gene_id"])
    hits: list[dict[str, object]] = []
    with runtime_path(blast_arg).open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != 12:
                raise ValueError(f"{blast_arg}:{line_number}: expected 12 BLASTP columns")
            query_id, gene_id = fields[:2]
            if gene_id not in valid_genes:
                continue
            hits.append(
                {
                    "query_id": query_id,
                    "gene_id": gene_id,
                    "identity": float(fields[2]),
                    "alignment_length": int(fields[3]),
                    "query_length": max(1, int(fields[4])),
                    "subject_length": max(1, int(fields[5])),
                    "query_start": int(fields[6]),
                    "query_end": int(fields[7]),
                    "subject_start": int(fields[8]),
                    "subject_end": int(fields[9]),
                    "evalue": float(fields[10]),
                    "bitscore": float(fields[11]),
                }
            )
    frame = pd.DataFrame(hits)
    rows: list[dict[str, str]] = []
    audit_rows: list[dict[str, object]] = []
    if not frame.empty:
        frame = frame.sort_values("bitscore", ascending=False, kind="stable").drop_duplicates(["query_id", "gene_id"])
        best_by_query = frame.groupby("query_id")["bitscore"].max().to_dict()
        frame["relative_bitscore"] = [float(row.bitscore) / float(best_by_query[row.query_id]) for row in frame.itertuples()]
        frame["query_coverage"] = ((frame["query_end"] - frame["query_start"]).abs() + 1) / frame["query_length"]
        frame["subject_coverage"] = ((frame["subject_end"] - frame["subject_start"]).abs() + 1) / frame["subject_length"]
        short = frame["query_length"] < 100
        full_length_quality = (
            (frame["identity"] >= 20.0)
            & (frame["query_coverage"] >= 0.40)
            & (frame["evalue"] <= 1e-3)
        )
        conserved_domain_quality = (
            (~short)
            & (frame["identity"] >= 30.0)
            & (frame["query_coverage"] >= 0.14)
            & (frame["alignment_length"] >= 50)
            & (frame["evalue"] <= 1e-10)
        )
        short_peptide_quality = (
            short
            & (frame["identity"] >= 50.0)
            & (frame["alignment_length"] >= 10)
            & (frame["relative_bitscore"] >= 0.90)
        )
        frame["quality_class"] = "fail"
        frame.loc[conserved_domain_quality, "quality_class"] = "conserved_domain_only"
        frame.loc[short_peptide_quality, "quality_class"] = "short_peptide_motif"
        frame.loc[full_length_quality, "quality_class"] = "full_length_quality"
        frame["quality_pass"] = full_length_quality | conserved_domain_quality | short_peptide_quality
        best_gene_hits = frame.sort_values(["relative_bitscore", "bitscore"], ascending=False, kind="stable").drop_duplicates("gene_id")
        for hit in best_gene_hits.itertuples():
            relative = float(hit.relative_bitscore)
            confidence = "high" if relative >= 0.98 else "medium" if relative >= 0.90 else "low"
            if str(hit.quality_class) != "full_length_quality" and confidence == "high":
                confidence = "medium"
            if max_confidence == "medium" and confidence == "high":
                confidence = "medium"
            if max_confidence == "low":
                confidence = "low"
            source_confidence = (query_confidence_map or {}).get(str(hit.query_id), "high")
            confidence_order = {"low": 0, "medium": 1, "high": 2}
            if source_confidence not in confidence_order:
                raise ValueError(f"invalid query confidence for {hit.query_id}: {source_confidence}")
            if confidence_order[source_confidence] < confidence_order[confidence]:
                confidence = source_confidence
            selected = bool(hit.quality_pass) and relative >= min_relative_bitscore
            audit_rows.append(
                {
                    "gene_id": hit.gene_id,
                    "query_id": hit.query_id,
                    "identity": round(float(hit.identity), 3),
                    "alignment_length": int(hit.alignment_length),
                    "query_coverage": round(float(hit.query_coverage), 6),
                    "subject_coverage": round(float(hit.subject_coverage), 6),
                    "evalue": hit.evalue,
                    "bitscore": hit.bitscore,
                    "relative_bitscore": round(relative, 6),
                    "quality_gate": "pass" if hit.quality_pass else "fail",
                    "quality_class": hit.quality_class,
                    "evidence_confidence": confidence if selected else "",
                }
            )
            if selected:
                rows.append(
                    _row(
                        str(hit.gene_id),
                        evidence_code,
                        "support",
                        confidence,
                        f"{hit.query_id};relative_bitscore={relative:.6f}",
                        "interval_homolog_blastp",
                        (query_source_map or {}).get(str(hit.query_id), source_id),
                        (
                        "Trait-linked homolog panel was searched only against proteins in the QTL interval: "
                            f"identity={hit.identity:.1f}; qcov={hit.query_coverage:.3f}; "
                            f"scov={hit.subject_coverage:.3f}; evalue={hit.evalue:g}; relative_bitscore={relative:.3f}; "
                            f"quality_class={hit.quality_class}."
                        ),
                    )
                )
    output = runtime_path(out_arg)
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=EVIDENCE_COLUMNS).to_csv(output, sep="\t", index=False)
    audit_path = output.with_name(output.stem + "_audit.tsv")
    pd.DataFrame(audit_rows).to_csv(audit_path, sep="\t", index=False)
    return {
        "output": str(output),
        "audit": str(audit_path),
        "hit_gene_count": len(audit_rows),
        "supported_gene_count": len(rows),
        "source_id": source_id,
        "warning": "Similarity tiers prioritize assays but do not prove orthology or causality.",
    }


def _normalized(text: object) -> str:
    return re.sub(r"\s+", " ", str(text).replace("_", " ").replace("-", " ").casefold()).strip()


def _eggnog_records(path: Path):
    header: list[str] | None = None
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#query"):
                header = raw.lstrip("#").rstrip("\r\n").split("\t")
                continue
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if not header:
                raise ValueError(f"{path}: eggNOG annotation header beginning with #query is required")
            if len(fields) < len(header):
                fields.extend([""] * (len(header) - len(fields)))
            yield dict(zip(header, fields))


def eggnog_trait_evidence(
    candidates_arg: str,
    eggnog_arg: str,
    terms: list[str],
    out_arg: str,
    source_id: str = DEFAULT_ANNOTATION_SOURCE,
    max_evalue: float = 1e-5,
    record_index: dict[str, dict[str, str]] | None = None,
) -> dict[str, object]:
    """Create auditable trait-function evidence from official eggNOG annotations."""
    if max_evalue <= 0:
        raise ValueError("max_evalue must be greater than zero")
    candidates = _candidate_frame(candidates_arg)
    mapping = _id_map(candidates)
    valid_genes = set(candidates["gene_id"])
    normalized_terms = sorted({_normalized(term) for term in terms if _normalized(term)})
    if not normalized_terms:
        raise ValueError("at least one non-empty trait term is required")

    rows: list[dict[str, str]] = []
    audit: list[dict[str, object]] = []
    seen_genes: set[str] = set()
    annotation_path = runtime_path(eggnog_arg)
    records: Iterable[dict[str, str]]
    if record_index is None:
        records = _eggnog_records(annotation_path)
    else:
        records = (record_index[key] for key in mapping if key in record_index)
    for record in records:
        query = record.get("query") or record.get("query_name") or ""
        gene_id = mapping.get(query, query)
        if gene_id not in valid_genes:
            continue
        seen_genes.add(gene_id)
        evalue_text = record.get("seed_ortholog_evalue", "")
        try:
            evalue = float(evalue_text)
        except ValueError:
            evalue = float("inf")
        fields = {
            "preferred_name": record.get("Preferred_name", ""),
            "description": record.get("Description") or record.get("eggNOG free text desc.", ""),
            "gos": record.get("GOs", ""),
            "ec": record.get("EC", ""),
            "kegg_ko": record.get("KEGG_ko", ""),
            "kegg_pathway": record.get("KEGG_Pathway", ""),
            "cazy": record.get("CAZy", ""),
        }
        searchable = _normalized(" | ".join(fields.values()))
        matched = [term for term in normalized_terms if term in searchable] if evalue <= max_evalue else []
        audit.append(
            {
                "gene_id": gene_id,
                "transcript_id": query,
                "seed_ortholog_evalue": evalue_text,
                "quality_gate": "pass" if evalue <= max_evalue else "fail",
                "matched_terms": ";".join(matched),
                **fields,
            }
        )
        if matched:
            description = "; ".join(
                f"{key}={value}" for key, value in fields.items() if str(value).strip() not in {"", "-"}
            )
            rows.append(
                _row(
                    gene_id,
                    "trait_function_match",
                    "support",
                    "medium",
                    ";".join(matched),
                    "official_eggnog_trait_terms",
                    source_id,
                    f"Explicit predeclared trait terms matched the official V2 annotation: {description}",
                )
            )

    output = runtime_path(out_arg)
    output.parent.mkdir(parents=True, exist_ok=True)
    evidence = pd.DataFrame(rows, columns=EVIDENCE_COLUMNS).drop_duplicates()
    evidence.to_csv(output, sep="\t", index=False)
    audit_path = output.with_name(output.stem + "_audit.tsv")
    pd.DataFrame(audit).to_csv(audit_path, sep="\t", index=False)
    return {
        "output": str(output),
        "audit": str(audit_path),
        "candidate_count": len(candidates),
        "annotated_candidate_count": len(seen_genes),
        "matched_candidate_count": int(evidence["gene_id"].nunique()) if not evidence.empty else 0,
        "terms": normalized_terms,
        "source_id": source_id,
    }


def _answer_ids(value: object) -> list[str]:
    return [item.strip() for item in str(value).split("|") if item.strip()]


def _evaluate_rank(ranked_path: Path, answer_gene_id: str) -> dict[str, object]:
    ranked = pd.read_csv(ranked_path, sep="\t", dtype=str).fillna("")
    accepted = _answer_ids(answer_gene_id)
    # benchmark_metrics is the single source of truth for the tie-aware metrics
    # (Top-N conservative/optimistic, MRR, wrong-sole-top, safe-hold, abstention).
    metrics = resolution_metrics(ranked, accepted)
    if not metrics["answer_present"]:
        return {
            **metrics,
            "top1": metrics["top1_conservative"],
            "top5": metrics["top5_conservative"],
            "top10": metrics["top10_conservative"],
            "top20": metrics["top20_conservative"],
            "top5_possible": metrics["top5_optimistic"],
            "top10_possible": metrics["top10_optimistic"],
            "dry_group": "answer_absent",
            "score_breakdown": "",
            "top_eligibility": "",
            "failure_reason": "answer_not_in_interval_annotation",
        }
    answer_row = ranked[ranked["gene_id"].astype(str) == str(metrics["matched_answer_gene_id"])].iloc[0]
    breakdown = str(answer_row.get("score_breakdown", ""))
    rank_max = int(metrics["rank_max"])
    if "curated_homolog_match" not in breakdown and "trait_function_match" not in breakdown:
        failure_reason = "answer_has_no_profile_specific_support"
    elif "curated_homolog_match" not in breakdown and "homolog_informed" in ranked_path.as_posix():
        failure_reason = "answer_homolog_not_supported"
    elif rank_max > 20:
        failure_reason = "answer_below_top20"
    elif int(metrics["tie_size"]) > 10:
        failure_reason = "large_functional_tie"
    elif rank_max > 10:
        failure_reason = "answer_below_top10"
    else:
        failure_reason = "recovered"
    return {
        **metrics,
        "top1": metrics["top1_conservative"],
        "top5": metrics["top5_conservative"],
        "top10": metrics["top10_conservative"],
        "top20": metrics["top20_conservative"],
        "top5_possible": metrics["top5_optimistic"],
        "top10_possible": metrics["top10_optimistic"],
        "dry_group": str(answer_row.get("dry_group", "")),
        "score_breakdown": breakdown,
        "top_eligibility": str(answer_row.get("top_eligibility", "")),
        "failure_reason": failure_reason,
    }


def _expand_profiles(config: pd.DataFrame) -> pd.DataFrame:
    if "trait_terms" in config.columns:
        return config.copy()
    profile_columns = [
        ("phenotype_terms", "phenotype_only"),
        ("mechanism_terms", "mechanism_informed"),
    ]
    if not any(column in config.columns for column, _ in profile_columns):
        raise ValueError("benchmark config requires trait_terms or phenotype_terms/mechanism_terms")
    rows: list[dict[str, object]] = []
    for _, source in config.iterrows():
        for column, profile in profile_columns:
            terms = str(source.get(column, "")).strip()
            if not terms:
                continue
            row = source.to_dict()
            row["locus_id"] = str(source["case_id"])
            row["case_id"] = f"{source['case_id']}_{profile}"
            row["term_profile"] = profile
            row["trait_terms"] = terms
            rows.append(row)
        query_id = str(source.get("homolog_query_id", "")).strip()
        if query_id:
            row = source.to_dict()
            row["locus_id"] = str(source["case_id"])
            row["case_id"] = f"{source['case_id']}_homolog_informed"
            row["term_profile"] = "homolog_informed"
            row["trait_terms"] = str(source.get("mechanism_terms") or source.get("phenotype_terms") or "")
            rows.append(row)
    return pd.DataFrame(rows).fillna("")


def _gff_gene_coordinates(path: Path) -> dict[str, tuple[str, int, int]]:
    result: dict[str, tuple[str, int, int]] = {}
    with open_text(path) as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = attributes(fields[8])
            gene_id = attrs.get("ID") or attrs.get("Name")
            if gene_id:
                result[str(gene_id)] = (fields[0], int(fields[3]), int(fields[4]))
    return result


def _reference_lengths(fasta_path: Path) -> dict[str, int]:
    index = Path(str(fasta_path) + ".fai")
    if not index.exists():
        # Only building the index needs pysam; a pre-built .fai is read as text.
        require_pysam("benchmark-reference on an unindexed FASTA").faidx(str(fasta_path))
    result: dict[str, int] = {}
    with index.open(encoding="utf-8") as handle:
        for raw in handle:
            fields = raw.rstrip("\r\n").split("\t")
            result[fields[0]] = int(fields[1])
    return result


def _virtual_interval(
    answer_gene_id: str,
    interval_size: int,
    coordinates: dict[str, tuple[str, int, int]],
    reference_lengths: dict[str, int],
) -> str:
    hits = [coordinates[gene_id] for gene_id in _answer_ids(answer_gene_id) if gene_id in coordinates]
    if not hits:
        raise ValueError(f"answer genes are absent from GFF3: {answer_gene_id}")
    seqids = {item[0] for item in hits}
    if len(seqids) != 1:
        raise ValueError(f"acceptable answer genes span multiple sequence IDs: {answer_gene_id}")
    seqid = hits[0][0]
    center = (min(item[1] for item in hits) + max(item[2] for item in hits)) // 2
    half = interval_size // 2
    start = max(1, center - half)
    end = min(reference_lengths[seqid], start + interval_size - 1)
    start = max(1, end - interval_size + 1)
    return f"{seqid}:{start}-{end}"


def _bool_mean(frame: pd.DataFrame, column: str) -> float:
    if not len(frame) or column not in frame.columns:
        return 0.0
    return float(frame[column].astype(str).isin({"True", "1", "1.0"}).mean())


def _metric_block(frame: pd.DataFrame) -> dict[str, object]:
    numeric_rr_c = pd.to_numeric(frame["reciprocal_rank_conservative"], errors="coerce") if "reciprocal_rank_conservative" in frame else pd.Series(dtype=float)
    numeric_rr_o = pd.to_numeric(frame["reciprocal_rank_optimistic"], errors="coerce") if "reciprocal_rank_optimistic" in frame else pd.Series(dtype=float)
    return {
        "case_count": len(frame),
        # Top-1 (quality of the first gene to sequence) and Top-20 (does the truth
        # survive to the shortlist) are separate numbers, not one headline.
        "top1_recall": _bool_mean(frame, "top1"),
        "top3_recall": _bool_mean(frame, "top3_conservative"),
        "top5_recall": _bool_mean(frame, "top5"),
        "top10_recall": _bool_mean(frame, "top10"),
        "top20_recall": _bool_mean(frame, "top20"),
        "mrr_conservative": round(float(numeric_rr_c.mean()), 6) if len(numeric_rr_c.dropna()) else 0.0,
        "mrr_optimistic": round(float(numeric_rr_o.mean()), 6) if len(numeric_rr_o.dropna()) else 0.0,
        "wrong_sole_top_rate": _bool_mean(frame, "wrong_paralog_sole_top"),
        "safe_hold_rate": _bool_mean(frame, "safe_hold"),
        "median_candidate_count": float(frame["candidate_count"].median()) if len(frame) else 0.0,
        "median_rank_max": float(pd.to_numeric(frame["rank_max"], errors="coerce").median()) if len(frame) else None,
    }


def _wilson_suffix(resolution: dict, low_key: str, high_key: str) -> str:
    """Append '(95% CI low-high)' when the two Wilson bounds are defined."""
    low = resolution.get(f"{low_key}_wilson_95_low")
    high = resolution.get(f"{high_key}_wilson_95_high")
    if low is None or high is None:
        return ""
    return f" (95% CI {low:.1%}-{high:.1%})"


def _mrr_ci_suffix(resolution: dict) -> str:
    """Append the conservative/optimistic MRR bootstrap CI when defined."""
    values = [resolution.get(key) for key in (
        "mrr_conservative_ci_95_low", "mrr_conservative_ci_95_high",
        "mrr_optimistic_ci_95_low", "mrr_optimistic_ci_95_high")]
    if any(value is None for value in values):
        return ""
    low_c, high_c, low_o, high_o = values
    return f" (95% CI {low_c:.3f}-{high_c:.3f} / {low_o:.3f}-{high_o:.3f})"


def _write_markdown_report(summary: pd.DataFrame, out_dir: Path, manifest: dict[str, object]) -> Path:
    lines = [
        "# Single-reference virtual-QTL benchmark report",
        "",
        f"- reference: `{manifest['reference_fasta']['path']}`",
        f"- expanded cases: {len(summary)}",
        f"- unique loci: {summary['locus_id'].nunique() if 'locus_id' in summary else len(summary)}",
        f"- LLM used for scoring: `{manifest['llm_used_for_scoring']}`",
        "- Top-k rule: the full score-tie group must fit within k.",
        "",
        "## Ranking quality by profile",
        "",
        "Top-20 recall (does the truth survive to the shortlist) is shown next to",
        "Top-1 quality, MRR, and the wrong-sole-top / safe-hold split -- they are",
        "distinct properties and must not be collapsed into one score.",
        "",
        "| profile | n | Top-1 | Top-3 | Top-5 | Top-20 | MRR | wrong-sole-top | safe-hold | median cand |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for profile, group in summary.groupby("term_profile", dropna=False):
        metrics = _metric_block(group)
        lines.append(
            f"| {profile} | {metrics['case_count']} | {metrics['top1_recall']:.1%} | {metrics['top3_recall']:.1%} | "
            f"{metrics['top5_recall']:.1%} | {metrics['top20_recall']:.1%} | {metrics['mrr_conservative']:.3f} | "
            f"{metrics['wrong_sole_top_rate']:.1%} | {metrics['safe_hold_rate']:.1%} | {metrics['median_candidate_count']:.0f} |"
        )
    resolution = manifest.get("resolution_metrics", {})
    if resolution:
        lines.extend([
            "",
            "## Resolution quality (all cases)",
            "",
            f"- answer-present cases: {resolution.get('n_answer_present', 0)}; truth-missing cases: {resolution.get('n_truth_missing', 0)}",
            f"- Top-1 / Top-20 recall (conservative): {resolution.get('top1_recall_conservative', 0):.1%} / {resolution.get('top20_recall_conservative', 0):.1%}"
            + _wilson_suffix(resolution, "top1_recall_conservative", "top20_recall_conservative"),
            f"- MRR (conservative / optimistic): {resolution.get('mrr_conservative', 0):.3f} / {resolution.get('mrr_optimistic', 0):.3f}"
            + _mrr_ci_suffix(resolution),
            f"- wrong sole top (leads on score): {resolution.get('wrong_sole_top_count', 0)} case(s); of those confidently promoted to Top: {resolution.get('confident_wrong_sole_top_count', 0)}; safe hold: {resolution.get('safe_hold_rate', 0):.1%}",
            f"- truth-missing abstention / false-confident: {resolution.get('abstention_rate', 0):.1%} / {resolution.get('false_confident_rate', 0):.1%}",
        ])
    if "truth_tier" in summary.columns:
        lines.extend(["", "## Recall by truth tier and profile", "", "| truth tier | profile | n | Top-5 | Top-10 | Top-20 |", "|---|---|---:|---:|---:|---:|"])
        for (tier, profile), group in summary.groupby(["truth_tier", "term_profile"], dropna=False):
            metrics = _metric_block(group)
            lines.append(
                f"| {tier} | {profile} | {metrics['case_count']} | {metrics['top5_recall']:.1%} | "
                f"{metrics['top10_recall']:.1%} | {metrics['top20_recall']:.1%} |"
            )
    failures = summary[~summary["top10"]].copy()
    lines.extend(["", "## Top-10 failures", ""])
    if failures.empty:
        lines.append("None.")
    else:
        lines.extend(["| case | gene | profile | candidates | rank range | tie | reason |", "|---|---|---|---:|---:|---:|---|"])
        for row in failures.sort_values(["term_profile", "rank_max"], na_position="last").itertuples():
            lines.append(
                f"| {row.case_id} | {getattr(row, 'gene_name', row.answer_gene_id)} | {row.term_profile} | "
                f"{row.candidate_count} | {row.rank_min or '-'}-{row.rank_max or '-'} | {row.tie_size} | {row.failure_reason} |"
            )
    path = out_dir / "benchmark_report.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def reference_benchmark(
    config_arg: str,
    fasta_arg: str,
    gff_arg: str,
    protein_arg: str,
    eggnog_arg: str,
    out_dir_arg: str,
    seqid_map_arg: str | None = None,
    annotation_source_id: str = DEFAULT_ANNOTATION_SOURCE,
    interval_size: int = 0,
    homolog_protein_args: list[str] | None = None,
    blastp_runner: str = "auto",
    num_threads: int = 4,
    homolog_map_arg: str | None = None,
    resume: bool = False,
    query_eligibility_manifest_arg: str | None = None,
) -> dict[str, object]:
    """Run blind-with-respect-to-scoring virtual-QTL recovery cases on one reference."""
    if interval_size < 0:
        raise ValueError("interval_size must be zero or greater")
    if num_threads < 1:
        raise ValueError("num_threads must be at least one")
    config_path = runtime_path(config_arg)
    source_config = pd.read_csv(config_path, sep="\t", dtype=str).fillna("")
    homolog_map_path: Path | None = None
    if homolog_map_arg:
        homolog_map_path = runtime_path(homolog_map_arg)
        homolog_map = pd.read_csv(homolog_map_path, sep="\t", dtype=str).fillna("")
        required_map = {"case_id", "homolog_query_id", "homolog_source_id"}
        missing_map = sorted(required_map.difference(homolog_map.columns))
        if missing_map:
            raise ValueError(f"homolog map is missing columns: {', '.join(missing_map)}")
        if homolog_map["case_id"].duplicated().any():
            raise ValueError("homolog map case_id values must be unique")
        overlap = set(source_config.columns).intersection(set(homolog_map.columns).difference({"case_id"}))
        if overlap:
            raise ValueError("homolog map columns already exist in benchmark config: " + ", ".join(sorted(overlap)))
        source_config = source_config.merge(homolog_map, on="case_id", how="left", validate="one_to_one").fillna("")
    config = _expand_profiles(source_config)
    required = {"case_id", "trait", "answer_gene_id", "trait_terms", "answer_source_id"}
    missing = sorted(required.difference(config.columns))
    if missing:
        raise ValueError(f"benchmark config is missing columns: {', '.join(missing)}")
    if config["case_id"].duplicated().any() or config["case_id"].str.strip().eq("").any():
        raise ValueError("benchmark case_id values must be unique and non-empty")

    audited_fingerprints: dict[str, str] | None = None
    query_eligibility_manifest_path: Path | None = None
    if query_eligibility_manifest_arg:
        query_eligibility_manifest_path = runtime_path(query_eligibility_manifest_arg)
        if not query_eligibility_manifest_path.exists():
            raise ValueError(f"query eligibility manifest does not exist: {query_eligibility_manifest_path}")
        q_manifest = json.loads(query_eligibility_manifest_path.read_text(encoding="utf-8"))
        if not q_manifest.get("strict", False):
            raise ValueError("query eligibility manifest must have strict=True for benchmark execution")
        if q_manifest.get("excluded_count", 0) > 0:
            raise ValueError(f"strict query eligibility manifest contains {q_manifest.get('excluded_count')} excluded queries")
        if q_manifest.get("pisum_queries_present", 0) > 0:
            raise ValueError(f"strict query eligibility manifest contains Pisum queries (pisum_queries_present={q_manifest.get('pisum_queries_present')})")
        if q_manifest.get("answer_derived_queries", 0) > 0:
            raise ValueError(f"strict query eligibility manifest contains answer-derived queries (answer_derived_queries={q_manifest.get('answer_derived_queries')})")
        if q_manifest.get("unknown_provenance_queries", 0) > 0:
            raise ValueError(f"strict query eligibility manifest contains unknown-provenance queries (unknown_provenance_queries={q_manifest.get('unknown_provenance_queries')})")
        audited_fingerprints = {}
        for q in q_manifest.get("queries", []):
            qid = str(q.get("query_id") or "").strip()
            sha = str(q.get("sha256") or "").strip()
            sid = str(q.get("source_id") or "").strip()
            if qid and sha:
                audited_fingerprints[qid] = sha
            if sid and sha:
                audited_fingerprints[sid] = sha
        if not audited_fingerprints and (query_eligibility_manifest_path.parent / "query_provenance_audit.tsv").exists():
            audit_df = pd.read_csv(query_eligibility_manifest_path.parent / "query_provenance_audit.tsv", sep="\t", dtype=str).fillna("")
            for row in audit_df.to_dict("records"):
                qid = str(row.get("homolog_query_id") or "").strip()
                sha = str(row.get("query_sha256") or "").strip()
                sid = str(row.get("homolog_source_id") or "").strip()
                if qid and sha:
                    audited_fingerprints[qid] = sha
                if sid and sha:
                    audited_fingerprints[sid] = sha

    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    previous_case_state: dict[str, dict[str, object]] = {}
    previous_manifest_path = out_dir / "benchmark_manifest.json"
    if resume and previous_manifest_path.exists():
        previous_manifest = json.loads(previous_manifest_path.read_text(encoding="utf-8"))
        for previous_case in previous_manifest.get("cases", []):
            homolog = previous_case.get("homolog_evidence") or {}
            query_rows = (homolog.get("queries") or {}).get("queries", [])
            previous_case_state[str(previous_case.get("case_id", ""))] = {
                "interval": str(previous_case.get("interval", "")),
                "query_ids": [str(row.get("requested_id", "")) for row in query_rows],
            }
    fasta_path = runtime_path(fasta_arg)
    gff_path = runtime_path(gff_arg)
    coordinates: dict[str, tuple[str, int, int]] = {}
    reference_lengths: dict[str, int] = {}
    if "interval" not in config.columns:
        config["interval"] = ""
    if config["interval"].str.strip().eq("").any():
        if interval_size <= 0:
            raise ValueError("blank intervals require --interval-size")
        coordinates = _gff_gene_coordinates(gff_path)
        reference_lengths = _reference_lengths(fasta_path)
        config["interval"] = config.apply(
            lambda row: str(row["interval"]).strip()
            or _virtual_interval(str(row["answer_gene_id"]), interval_size, coordinates, reference_lengths),
            axis=1,
        )
    eggnog_index = {
        str(record.get("query") or record.get("query_name") or ""): record
        for record in _eggnog_records(runtime_path(eggnog_arg))
    }
    homolog_sequences: dict[str, tuple[object, str]] = {}
    if config.get("term_profile", pd.Series(dtype=str)).eq("homolog_informed").any():
        if not homolog_protein_args:
            raise ValueError("homolog_informed profiles require at least one --homolog-proteins FASTA")
        homolog_sequences = _homolog_sequence_index(homolog_protein_args)
    summary_rows: list[dict[str, object]] = []
    case_manifests: list[dict[str, object]] = []
    shared_cache: dict[str, dict[str, object]] = {}
    resumed_interval_count = 0
    reused_blast_count = 0
    for _, case in config.iterrows():
        started = perf_counter()
        case_id = str(case["case_id"])
        case_dir = out_dir / case_id
        interval = str(case["interval"])
        shared_key = hashlib.sha256(interval.encode("utf-8")).hexdigest()[:12]
        if shared_key not in shared_cache:
            shared_dir = out_dir / "00_shared" / f"{safe_id(interval)}_{shared_key}"
            extract_manifest = shared_dir / "01_extract" / "extract_manifest.json"
            official_path = shared_dir / "02_annotation" / "official_candidate_proteins.fa"
            evidence_manifest = shared_dir / "03_base_evidence" / "evidence_manifest.json"
            if resume and extract_manifest.exists() and official_path.exists() and evidence_manifest.exists():
                extract_result = json.loads(extract_manifest.read_text(encoding="utf-8"))
                if str(extract_result.get("interval")) != interval:
                    raise ValueError(f"resume interval mismatch in {extract_manifest}")
                candidate_path = Path(str(extract_result["outputs"]["candidates_tsv"]))
                official_subset = {
                    "input_candidates": str(candidate_path),
                    "input_proteins": str(runtime_path(protein_arg)),
                    "output": str(official_path),
                    "selected_proteins": sum(1 for _ in SeqIO.parse(str(official_path), "fasta")),
                    "unmatched_candidate_ids": [],
                    "resumed": True,
                }
                base_evidence = json.loads(evidence_manifest.read_text(encoding="utf-8")) | {
                    "manifest": str(evidence_manifest),
                    "resumed": True,
                }
                resumed_interval_count += 1
            else:
                extract_result = extract_candidates(
                    fasta_arg,
                    gff_arg,
                    interval,
                    str(shared_dir / "01_extract"),
                    seqid_map_arg,
                )
                candidate_path = Path(str(extract_result["outputs"]["candidates_tsv"]))
                official_subset = subset_proteins(
                    str(candidate_path),
                    protein_arg,
                    str(official_path),
                )
                base_evidence = evidence_build(
                    str(candidate_path),
                    str(shared_dir / "03_base_evidence"),
                    eggnog_args=[eggnog_arg],
                    interval=interval,
                )
            shared_cache[shared_key] = {
                "interval": interval,
                "extract": extract_result,
                "candidate_path": candidate_path,
                "official_subset": official_subset,
                "base_evidence": base_evidence,
            }
        shared = shared_cache[shared_key]
        extract_result = shared["extract"]
        candidate_path = Path(str(shared["candidate_path"]))
        official_subset = shared["official_subset"]
        term_result = eggnog_trait_evidence(
            str(candidate_path),
            eggnog_arg,
            [term for term in str(case["trait_terms"]).split("|") if term],
            str(case_dir / "02_annotation" / "trait_evidence.tsv"),
            annotation_source_id,
            record_index=eggnog_index,
        )
        manual_evidence_paths = [Path(str(term_result["output"]))]
        homolog_result: dict[str, object] | None = None
        if str(case.get("term_profile", "")) == "homolog_informed":
            query_ids = _answer_ids(case.get("homolog_query_id", ""))
            query_result = _write_homolog_queries(
                query_ids,
                homolog_sequences,
                case_dir / "02_annotation" / "homolog_queries.fa",
                audited_fingerprints=audited_fingerprints,
            )
            blast_path = case_dir / "02_annotation" / "homolog_vs_interval.tsv"
            previous = previous_case_state.get(case_id, {})
            reusable_blast = (
                resume
                and blast_path.exists()
                and previous.get("interval") == interval
                and previous.get("query_ids") == query_ids
            )
            if reusable_blast:
                blast_result = {
                    "output": str(blast_path),
                    "runner": "reused",
                    "short_query_mode": None,
                    "command": [],
                    "stderr": "",
                }
                reused_blast_count += 1
            else:
                blast_result = _run_interval_blastp(
                    Path(str(query_result["output"])),
                    Path(str(official_subset["output"])),
                    blast_path,
                    blastp_runner,
                    num_threads,
                )
            homolog_result = {
                "queries": query_result,
                "blastp": blast_result,
            }
            homolog_evidence = homolog_blast_evidence(
                str(candidate_path),
                str(blast_result["output"]),
                str(case_dir / "02_annotation" / "homolog_evidence.tsv"),
                str(case.get("homolog_source_id") or case.get("answer_source_id") or "curated_homolog"),
            )
            homolog_result["evidence"] = homolog_evidence
            manual_evidence_paths.append(Path(str(homolog_evidence["output"])))
        combined_manual = case_dir / "02_annotation" / "combined_manual_evidence.tsv"
        manual_frames = [pd.read_csv(path, sep="\t", dtype=str).fillna("") for path in manual_evidence_paths]
        pd.concat(manual_frames, ignore_index=True).drop_duplicates().to_csv(combined_manual, sep="\t", index=False)
        evidence_result = evidence_build(
            str(candidate_path),
            str(case_dir / "03_evidence"),
            manual_evidence_arg=str(combined_manual),
            existing_evidence_arg=str(shared["base_evidence"]["evidence"]),
            interval=interval,
        )
        rank_result = prioritize_evidence(
            str(candidate_path),
            str(evidence_result["evidence"]),
            str(case_dir / "04_rank"),
            top_n=5,
            reserve_n=5,
            rescue_n=3,
        )
        evaluation = _evaluate_rank(Path(str(rank_result["ranked_candidates"])), str(case["answer_gene_id"]))
        elapsed = round(perf_counter() - started, 3)
        summary_rows.append(
            {
                "case_id": case_id,
                "locus_id": case.get("locus_id", case_id),
                "term_profile": case.get("term_profile", "unspecified"),
                "truth_tier": case.get("truth_tier", "unspecified"),
                "gene_class": case.get("gene_class", "unspecified"),
                "gene_name": case.get("gene_name", ""),
                "homolog_source_type": case.get("homolog_source_type", ""),
                "trait": case["trait"],
                "interval": interval,
                "interval_size": int(interval.split("-")[-1]) - int(interval.split(":")[-1].split("-")[0]) + 1,
                "candidate_count": extract_result["candidate_count"],
                "official_protein_count": official_subset["selected_proteins"],
                "annotated_candidate_count": term_result["annotated_candidate_count"],
                "trait_matched_candidate_count": term_result["matched_candidate_count"],
                "homolog_supported_candidate_count": (
                    int(homolog_result["evidence"]["supported_gene_count"]) if homolog_result else 0
                ),
                "answer_gene_id": case["answer_gene_id"],
                "answer_source_id": case["answer_source_id"],
                "answer_mapping_source": case.get("answer_mapping_source", ""),
                **evaluation,
                "rank_fraction": (
                    round(float(evaluation["rank_max"]) / float(extract_result["candidate_count"]), 6)
                    if evaluation["rank_max"] is not None and extract_result["candidate_count"]
                    else None
                ),
                "elapsed_seconds": elapsed,
            }
        )
        case_manifests.append(
            {
                "case_id": case_id,
                "locus_id": case.get("locus_id", case_id),
                "truth_tier": case.get("truth_tier", "unspecified"),
                "gene_class": case.get("gene_class", "unspecified"),
                "homolog_source_type": case.get("homolog_source_type", ""),
                "term_profile": case.get("term_profile", "unspecified"),
                "interval": interval,
                "shared_interval_key": shared_key,
                "trait_terms": term_result["terms"],
                "extract": extract_result,
                "official_protein_subset": official_subset,
                "trait_evidence": term_result,
                "homolog_evidence": homolog_result,
                "evidence": evidence_result,
                "rank": rank_result,
                "evaluation": evaluation,
                "answer_gene_used_for_scoring": False,
            }
        )

    summary = pd.DataFrame(summary_rows)
    summary_path = out_dir / "benchmark_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    profile_metrics = {
        str(profile): _metric_block(group)
        for profile, group in summary.groupby("term_profile", dropna=False)
    }
    truth_tier_metrics = {
        str(tier): _metric_block(group)
        for tier, group in summary.groupby("truth_tier", dropna=False)
    }
    # Tie-aware resolution metrics: Top-1 quality, MRR, wrong-sole-top, and
    # abstention are reported alongside -- never folded into -- Top-20 recall.
    resolution_summary = aggregate_resolution_metrics(summary_rows)
    # Issue #367: Wilson 95% CIs for the headline Top-1 / Top-20 conservative
    # recall, over the same answer-present population as the recall itself.
    present_summary = [row for row in summary_rows if row.get("answer_present")]
    n_answer_present = len(present_summary)
    top1_hits = sum(1 for row in present_summary if row.get("top1_conservative"))
    top20_hits = sum(1 for row in present_summary if row.get("top20_conservative"))
    top1_low, top1_high = wilson_interval(top1_hits, n_answer_present)
    top20_low, top20_high = wilson_interval(top20_hits, n_answer_present)
    resolution_summary["top1_recall_conservative_wilson_95_low"] = top1_low
    resolution_summary["top1_recall_conservative_wilson_95_high"] = top1_high
    resolution_summary["top20_recall_conservative_wilson_95_low"] = top20_low
    resolution_summary["top20_recall_conservative_wilson_95_high"] = top20_high
    manifest = {
        "benchmark_type": "single-reference virtual-QTL known-gene recovery",
        "reference_fasta": file_fingerprint(runtime_path(fasta_arg)),
        "reference_gff3": file_fingerprint(runtime_path(gff_arg)),
        "official_protein_fasta": file_fingerprint(runtime_path(protein_arg)),
        "official_eggnog_annotations": file_fingerprint(runtime_path(eggnog_arg)),
        "config": str(config_path),
        "homolog_map": str(homolog_map_path) if homolog_map_path else None,
        "source_locus_count": len(source_config),
        "case_count": len(summary),
        "unique_interval_count": len(shared_cache),
        "virtual_interval_size": interval_size or None,
        "top1_recall": float(summary["top1"].mean()) if len(summary) else 0.0,
        "top5_recall": float(summary["top5"].mean()) if len(summary) else 0.0,
        "top10_recall": float(summary["top10"].mean()) if len(summary) else 0.0,
        "top20_recall": float(summary["top20"].mean()) if len(summary) else 0.0,
        "resolution_metrics": resolution_summary,
        "metric_policy": (
            "Top-20 recall (map-based) is reported but is not a standalone overall score. "
            "Metric changes are compared under a pre-defined protocol; weight fitting to benchmark answers is prohibited."
        ),
        "profile_metrics": profile_metrics,
        "truth_tier_metrics": truth_tier_metrics,
        "llm_used_for_scoring": False,
        "homolog_protein_fastas": [file_fingerprint(runtime_path(path)) for path in (homolog_protein_args or [])],
        "blastp_runner": blastp_runner,
        "resume_requested": resume,
        "resumed_interval_count": resumed_interval_count,
        "reused_blast_count": reused_blast_count,
        "pisum_queries_executed": False,
        "query_eligibility_manifest": file_fingerprint(query_eligibility_manifest_path) if query_eligibility_manifest_path else None,
        "answer_gene_policy": "answer IDs are read only after evidence construction and ranking",
        "cases": case_manifests,
        "summary": str(summary_path),
    }
    manifest_path = out_dir / "benchmark_manifest.json"
    report_path = _write_markdown_report(summary, out_dir, manifest)
    manifest["report"] = str(report_path)
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return {
        "summary": str(summary_path),
        "manifest": str(manifest_path),
        "report": str(report_path),
        "case_count": len(summary),
        "top1_recall": manifest["top1_recall"],
        "top5_recall": manifest["top5_recall"],
        "top10_recall": manifest["top10_recall"],
        "top20_recall": manifest["top20_recall"],
        "resolution_metrics": resolution_summary,
        "profile_metrics": profile_metrics,
    }
