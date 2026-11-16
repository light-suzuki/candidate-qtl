from __future__ import annotations

import json
import math
import shlex
import gzip
import re
from urllib.parse import unquote
from pathlib import Path
from statistics import mean

import pandas as pd
from Bio import SeqIO

from . import homology_qc, schemas
from .core import nonempty, runtime_path, safe_id, sha256
# Shared evidence I/O, vocabulary, and provenance helpers live in ``evidence_io``
# so ranking / benchmark / discovery / orthology code depends on one explicit
# boundary instead of importing private helpers from here. They keep their
# ``_``-prefixed local names below for the many internal call sites.
from .evidence_io import (
    CONFIDENCE_MULTIPLIERS,
    EVIDENCE_COLUMNS,
    annotation_aliases as _annotation_aliases,
    best_hits as _best_hits,
    candidate_frame as _candidate_frame,
    evidence_row as _row,
    id_map as _id_map,
)

VALID_DIRECTIONS = set(schemas.DIRECTIONS)


def annotation_keyword_evidence(
    candidates_arg: str,
    blast_arg: str,
    annotation_gff_arg: str | None,
    keyword_args: list[str],
    out_arg: str,
    subject_annotations_arg: str | None = None,
    use_stitle: bool = False,
) -> dict[str, object]:
    """Link best protein hits to GFF descriptions and score explicit trait keywords."""
    candidates = _candidate_frame(candidates_arg)
    mapping = _id_map(candidates)
    valid_genes = set(candidates["gene_id"])
    annotations: dict[str, tuple[str, str]] = {}
    annotation_path = runtime_path(annotation_gff_arg) if annotation_gff_arg else None
    if annotation_path:
        opener = gzip.open if annotation_path.suffix == ".gz" else open
        with opener(annotation_path, "rt", encoding="utf-8") as handle:
            for raw in handle:
                if not raw.strip() or raw.startswith("#"):
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) != 9 or fields[2] != "gene":
                    continue
                attrs = {}
                for item in fields[8].split(";"):
                    if "=" in item:
                        key, value = item.split("=", 1)
                        attrs[unquote(key)] = unquote(value)
                identifier = attrs.get("ID") or attrs.get("Name")
                description = attrs.get("Note") or attrs.get("description") or attrs.get("Name", "")
                if identifier:
                    for alias in _annotation_aliases(identifier):
                        annotations[alias] = (identifier, description)
    if subject_annotations_arg:
        subject_path = runtime_path(subject_annotations_arg)
        table = pd.read_csv(subject_path, sep="\t", dtype=str, header=None, comment="#").fillna("")
        if table.shape[1] < 2:
            raise ValueError("subject annotations must have at least two TSV columns: subject_id and description")
        for subject, description in table.iloc[:, :2].itertuples(index=False, name=None):
            for alias in _annotation_aliases(subject):
                annotations[alias] = (subject, description)
    best: dict[str, tuple[float, str]] = {}
    blast_path = runtime_path(blast_arg)
    with blast_path.open(encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 12:
                raise ValueError("annotation keyword BLAST must use the 12/13-column candidate-qtl outfmt")
            gene_id = mapping.get(fields[0], fields[0])
            if gene_id not in valid_genes:
                continue
            bitscore = float(fields[11])
            if gene_id not in best or bitscore > best[gene_id][0]:
                best[gene_id] = (bitscore, fields[1])
                if use_stitle:
                    if len(fields) < 13:
                        raise ValueError("--use-stitle requires BLAST outfmt with stitle as column 13")
                    annotations[fields[1]] = (fields[1], fields[12])
    keywords = sorted({term.strip().lower() for value in keyword_args for term in value.split(",") if term.strip()})
    if not keywords:
        raise ValueError("at least one non-empty --keyword is required")
    rows = []
    audit_rows = []
    for gene_id, (bitscore, subject) in best.items():
        annotation = next((annotations[alias] for alias in _annotation_aliases(subject) if alias in annotations), None)
        target_id, description = annotation if annotation else ("", "")
        lowered = description.lower()
        matched = [term for term in keywords if re.search(rf"(?<![0-9a-z]){re.escape(term)}(?![0-9a-z])", lowered)]
        audit_rows.append(
            {"gene_id": gene_id, "best_subject": subject, "target_gene_id": target_id, "bitscore": bitscore, "matched_keywords": ";".join(matched), "description": description}
        )
        if matched:
            rows.append(
                _row(
                    gene_id,
                    "trait_function_match",
                    "support",
                    "medium",
                    ";".join(matched),
                    "gff_annotation_keyword",
                    annotation_path or runtime_path(subject_annotations_arg or blast_arg),
                    f"Best-hit annotation matched explicit trait keywords ({', '.join(matched)}): {description}",
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
        "candidate_count": len(candidates),
        "annotated_best_hits": sum(bool(row["target_gene_id"]) for row in audit_rows),
        "matched_candidates": len(rows),
        "keywords": keywords,
        "warning": "Keyword matches are functional hypotheses, not causal evidence; keep the audit table and validate parental variants and cosegregation.",
    }


def _parse_blast(path: Path, mapping: dict[str, str], source_type: str) -> list[dict[str, str]]:
    grouped: dict[str, list[dict[str, object]]] = {}
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 7:
                raise ValueError(f"{path}: expected BLAST/DIAMOND outfmt 6 with at least 7 columns")
            query, subject = fields[0], fields[1]
            gene_id = mapping.get(query, query)
            if gene_id not in set(mapping.values()):
                continue
            identity, align_length = float(fields[2]), int(fields[3])
            if len(fields) >= 13:
                query_length = max(1, int(fields[4]))
                subject_length = max(1, int(fields[5]))
                evalue, bitscore, title = float(fields[10]), float(fields[11]), fields[12]
                query_coverage = align_length / query_length
                subject_coverage = align_length / subject_length
            else:
                evalue, bitscore, title = float(fields[4]), float(fields[5]), fields[6]
                query_coverage = math.nan
                subject_coverage = math.nan
            grouped.setdefault(gene_id, []).append(
                {
                    "subject": subject,
                    "identity": identity,
                    "query_coverage": query_coverage,
                    "subject_coverage": subject_coverage,
                    "evalue": evalue,
                    "bitscore": bitscore,
                    "title": title,
                }
            )
    rows: list[dict[str, str]] = []
    for gene_id, hits in grouped.items():
        hits.sort(key=lambda hit: (-float(hit["bitscore"]), float(hit["evalue"]), str(hit["subject"])))
        best = hits[0]
        coverage = float(best["query_coverage"])
        subject_coverage = float(best["subject_coverage"])
        confidence = homology_qc.confidence(best["identity"], coverage, subject_coverage, best["evalue"])
        quality = homology_qc.quality_class(best["identity"], coverage, subject_coverage)
        coverage_text = "unknown" if math.isnan(coverage) else f"{coverage:.3f}"
        subject_coverage_text = "unknown" if math.isnan(subject_coverage) else f"{subject_coverage:.3f}"
        rows.append(
            _row(
                gene_id,
                "homology_context",
                "context",
                confidence,
                best["subject"],
                source_type,
                path,
                f"generic best-hit conservation (context, not trait evidence; {quality}): best hit={best['subject']}; identity={best['identity']}; qcov={coverage_text}; scov={subject_coverage_text}; evalue={best['evalue']}; {best['title']}",
            )
        )
        if len(hits) > 1 and homology_qc.is_paralog_ambiguous(best["bitscore"], hits[1]["bitscore"]):
            rows.append(
                _row(
                    gene_id,
                    "paralog_ambiguous",
                    "against",
                    "medium",
                    f"top2_bitscore_ratio={homology_qc.bitscore_ratio(best['bitscore'], hits[1]['bitscore']):.3f}",
                    source_type,
                    path,
                    "Top two similarity hits are too close for a simple one-hit ortholog assumption.",
                )
            )
    return rows


def _parse_hmmer(path: Path, mapping: dict[str, str]) -> list[dict[str, str]]:
    best: dict[tuple[str, str], tuple[float, float]] = {}
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.split()
            if len(fields) < 8:
                continue
            domain, query = fields[0], fields[2]
            gene_id = mapping.get(query, query)
            if gene_id not in set(mapping.values()):
                continue
            evalue, score = float(fields[4]), float(fields[5])
            key = (gene_id, domain)
            if key not in best or score > best[key][1]:
                best[key] = (evalue, score)
    return [
        _row(
            gene_id,
            "domain_context",
            "context",
            "high" if evalue <= 1e-10 else "medium",
            domain,
            "hmmscan",
            path,
            f"generic domain detection (context, not trait evidence): domain={domain}; evalue={evalue}; score={score}",
        )
        for (gene_id, domain), (evalue, score) in sorted(best.items())
        if evalue <= 1e-3
    ]


def _parse_eggnog(path: Path, mapping: dict[str, str]) -> list[dict[str, str]]:
    header: list[str] | None = None
    rows: list[dict[str, str]] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#query"):
                header = raw.lstrip("#").rstrip("\r\n").split("\t")
                continue
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if header and len(fields) >= len(header):
                record = dict(zip(header, fields))
                query = record.get("query") or record.get("query_name") or fields[0]
                orthogroups = record.get("eggNOG_OGs") or record.get("eggNOG OGs", "")
                preferred_name = record.get("Preferred_name", "")
                description = (
                    record.get("Description")
                    or record.get("eggNOG free text desc.")
                    or preferred_name
                )
                gos = record.get("GOs", "")
            else:
                query = fields[0]
                orthogroups = fields[4] if len(fields) > 4 else ""
                description = fields[7] if len(fields) > 7 else ""
                gos = fields[9] if len(fields) > 9 else ""
            gene_id = mapping.get(query, query)
            if gene_id not in set(mapping.values()):
                continue
            if nonempty(orthogroups) and orthogroups != "-":
                rows.append(_row(gene_id, "orthogroup_supported", "support", "medium", orthogroups, "eggnog", path, description))
            if nonempty(gos) and gos != "-":
                rows.append(_row(gene_id, "go_annotation_context", "context", "low", gos, "eggnog_go", path, f"generic GO annotation presence (context, not trait evidence): {description}"))
    return rows


def _parse_expression(
    path: Path,
    mapping: dict[str, str],
    target_tissues: set[str],
    target_stages: set[str],
    min_expression: float,
    target_fold: float,
    target_conditions: set[str] | None = None,
) -> tuple[list[dict[str, str]], pd.DataFrame]:
    """Import long-form expression evidence without mixing normalization units.

    Optional columns (``stage``, ``condition``, ``unit``, ``replicate``,
    ``source_id``) refine the comparison. Target vs background enrichment is
    computed per ``unit`` so TPM, raw counts, and microarray intensities are
    never divided into one another. Expression is treated as functional
    evidence, never as causal proof.
    """

    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    required = {"gene_id", "tissue", "value"}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"expression table is missing required columns: {', '.join(missing)}")
    frame["gene_id"] = frame["gene_id"].map(lambda value: mapping.get(value, value))
    frame["_value"] = pd.to_numeric(frame["value"], errors="coerce")
    if frame["_value"].isna().any():
        raise ValueError("expression value must be numeric")
    frame["_unit"] = frame["unit"].str.strip() if "unit" in frame.columns else ""
    frame["_condition"] = frame["condition"].str.strip() if "condition" in frame.columns else ""

    if "is_target" in frame.columns and frame["is_target"].str.strip().ne("").any():
        frame["_target"] = frame["is_target"].str.strip().str.lower().isin({"1", "true", "yes", "target"})
    else:
        if not target_tissues:
            raise ValueError("expression import requires --target-tissue or an is_target column")
        tissue_match = frame["tissue"].str.strip().str.casefold().isin(target_tissues)
        if target_stages:
            if "stage" not in frame.columns:
                raise ValueError("expression table requires stage when --target-stage is used")
            stage_match = frame["stage"].str.strip().str.casefold().isin(target_stages)
            frame["_target"] = tissue_match & stage_match
        else:
            frame["_target"] = tissue_match
    if target_conditions:
        if "condition" not in frame.columns:
            raise ValueError("expression table requires a condition column when --target-condition is used")
        frame["_target"] = frame["_target"] & frame["_condition"].str.casefold().isin(
            {value.casefold() for value in target_conditions}
        )

    rows: list[dict[str, str]] = []
    summaries: list[dict[str, object]] = []
    for (gene_id, unit), group in frame.groupby(["gene_id", "_unit"], sort=True):
        targets = [float(value) for value in group.loc[group["_target"], "_value"]]
        background = [float(value) for value in group.loc[~group["_target"], "_value"]]
        if not targets:
            continue
        target_mean = mean(targets)
        background_mean = mean(background) if background else 0.0
        fold = target_mean / max(background_mean, 1e-9) if background else math.inf
        unit_label = unit or "unspecified_unit"
        conditions = sorted({value for value in group["_condition"].tolist() if value})
        source_ids = [value for value in group.get("source_id", pd.Series(dtype=str)).tolist() if nonempty(value)]
        source_id = source_ids[0] if source_ids else str(path)
        specificity = fold if math.isfinite(fold) else ""
        expressed = target_mean >= min_expression
        has_background = bool(background)
        # Detection is not specificity: without a same-unit background sample we
        # cannot claim the gene is target-specific, so the fold is not computed.
        if not has_background:
            specificity_class = "comparison_unavailable"
        elif math.isfinite(fold) and fold >= target_fold:
            specificity_class = "target_specific"
        else:
            specificity_class = "broadly_expressed"
        summaries.append(
            {
                "gene_id": gene_id,
                "unit": unit_label,
                "target_mean": target_mean,
                "background_mean": background_mean if has_background else "",
                "target_to_background_fold": specificity,
                "specificity_class": specificity_class,
                "expressed": "yes" if expressed else "no",
                "non_expression": "yes" if not expressed else "no",
                "comparison_available": "yes" if has_background else "no",
                "conditions_observed": ";".join(conditions),
                "target_observations": len(targets),
                "background_observations": len(background),
                "source_id": source_id,
            }
        )
        unit_note = f"unit={unit_label}; functional evidence only, not causal."
        if expressed:
            rows.append(_row(gene_id, "expression_detected", "support", "medium", target_mean, "expression", source_id, f"Expression detected in target tissue/stage. {unit_note}"))
            if not has_background:
                # No background measured: record the sampling gap as context, do
                # not convert it into positive target-specificity support.
                rows.append(_row(gene_id, "expression_comparison_unavailable", "context", "low", target_mean, "expression", source_id, f"Detected in the target, but no same-unit background was measured; specificity is unknown. {unit_note}"))
            elif fold >= target_fold:
                rows.append(_row(gene_id, "expression_target_specific", "support", "high", specificity, "expression", source_id, f"Target expression enriched over supplied background within the same unit. {unit_note}"))
        else:
            rows.append(_row(gene_id, "expression_absent", "against", "low", target_mean, "expression", source_id, f"Target expression is below the supplied threshold. {unit_note}"))
    return rows, pd.DataFrame(summaries)


def _position_evidence(
    candidates: pd.DataFrame,
    peak: str | None,
    interval: str | None,
    peak_window: int,
    boundary_window: int,
) -> tuple[list[dict[str, str]], pd.DataFrame]:
    # Positional interpretation lives in ranking_policy; the local import keeps
    # the dependency one-directional (ranking_policy imports shared helpers from
    # evidence_io, not from here).
    from .ranking_policy import position_evidence

    return position_evidence(candidates, peak, interval, peak_window, boundary_window)


def _validate_manual_evidence(path: Path, valid_genes: set[str]) -> list[dict[str, str]]:
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    missing = sorted(set(EVIDENCE_COLUMNS).difference(frame.columns))
    if missing:
        raise ValueError(f"manual evidence table is missing required columns: {', '.join(missing)}")
    if not set(frame["gene_id"]).issubset(valid_genes):
        raise ValueError("manual evidence contains gene_id values absent from the candidate table")
    if not set(frame["direction"]).issubset(VALID_DIRECTIONS):
        raise ValueError(f"manual evidence direction must be one of: {', '.join(sorted(VALID_DIRECTIONS))}")
    if not set(frame["confidence"]).issubset(CONFIDENCE_MULTIPLIERS):
        raise ValueError("manual evidence confidence must be high, medium, or low")
    if frame["source_id"].map(nonempty).eq(False).any():
        raise ValueError("manual evidence requires a source_id for every row")
    return frame[EVIDENCE_COLUMNS].to_dict("records")


def evidence_build(
    candidates_arg: str,
    out_dir_arg: str,
    blast_args: list[str] | None = None,
    diamond_args: list[str] | None = None,
    hmmer_args: list[str] | None = None,
    eggnog_args: list[str] | None = None,
    expression_arg: str | None = None,
    manual_evidence_arg: str | None = None,
    existing_evidence_arg: str | None = None,
    target_tissues: list[str] | None = None,
    target_stages: list[str] | None = None,
    min_expression: float = 1.0,
    target_fold: float = 2.0,
    peak: str | None = None,
    interval: str | None = None,
    peak_window: int = 0,
    boundary_window: int = 0,
    target_conditions: list[str] | None = None,
) -> dict[str, object]:
    candidates = _candidate_frame(candidates_arg)
    mapping = _id_map(candidates)
    valid_genes = set(candidates["gene_id"])
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    for arg in blast_args or []:
        rows.extend(_parse_blast(runtime_path(arg), mapping, "blastp"))
    for arg in diamond_args or []:
        rows.extend(_parse_blast(runtime_path(arg), mapping, "diamond"))
    for arg in hmmer_args or []:
        rows.extend(_parse_hmmer(runtime_path(arg), mapping))
    for arg in eggnog_args or []:
        rows.extend(_parse_eggnog(runtime_path(arg), mapping))
    expression_summary = pd.DataFrame()
    if expression_arg:
        target_set = {value.strip().casefold() for value in target_tissues or []}
        stage_set = {value.strip().casefold() for value in target_stages or []}
        condition_set = {value.strip() for value in target_conditions or [] if value.strip()}
        expression_rows, expression_summary = _parse_expression(runtime_path(expression_arg), mapping, target_set, stage_set, min_expression, target_fold, condition_set or None)
        rows.extend(expression_rows)
    position_rows, position_summary = _position_evidence(candidates, peak, interval, peak_window, boundary_window)
    rows.extend(position_rows)
    if manual_evidence_arg:
        rows.extend(_validate_manual_evidence(runtime_path(manual_evidence_arg), valid_genes))
    if existing_evidence_arg:
        rows.extend(_validate_manual_evidence(runtime_path(existing_evidence_arg), valid_genes))
    for _, candidate in candidates.iterrows():
        quality = str(candidate.get("annotation_quality", ""))
        cds_length = pd.to_numeric(pd.Series([candidate.get("cds_length", "0")]), errors="coerce").fillna(0).iloc[0]
        if quality in {"no_cds", "partial_or_unverified_orf"} or float(cds_length) == 0:
            rows.append(_row(str(candidate["gene_id"]), "annotation_incomplete", "rescue", "high", quality or "no_cds", "gff3", candidates_arg, "Incomplete coding annotation must not be silently discarded."))
    evidence = pd.DataFrame(rows, columns=EVIDENCE_COLUMNS).drop_duplicates()
    evidence = evidence[evidence["gene_id"].isin(valid_genes)].sort_values(["gene_id", "evidence_code", "source_id"], kind="stable")
    evidence_path = out_dir / "evidence.tsv"
    evidence.to_csv(evidence_path, sep="\t", index=False)
    if not expression_summary.empty:
        expression_summary.to_csv(out_dir / "expression_summary.tsv", sep="\t", index=False)
    if not position_summary.empty:
        position_summary.to_csv(out_dir / "position_summary.tsv", sep="\t", index=False)
    summary = candidates.copy()
    if not evidence.empty:
        codes = evidence.groupby("gene_id")["evidence_code"].apply(lambda values: ";".join(sorted(set(values)))).rename("evidence_codes")
        counts = evidence.groupby("gene_id").size().rename("evidence_count")
        summary = summary.merge(codes, on="gene_id", how="left").merge(counts, on="gene_id", how="left")
    if "evidence_codes" not in summary.columns:
        summary["evidence_codes"] = ""
    else:
        summary["evidence_codes"] = summary["evidence_codes"].fillna("")
    if "evidence_count" not in summary.columns:
        summary["evidence_count"] = 0
    else:
        summary["evidence_count"] = summary["evidence_count"].fillna(0).astype(int)
    summary_path = out_dir / "candidate_evidence_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    input_args = [candidates_arg, *(blast_args or []), *(diamond_args or []), *(hmmer_args or []), *(eggnog_args or [])]
    input_args.extend(value for value in [expression_arg, manual_evidence_arg, existing_evidence_arg] if value)
    fingerprints = []
    for value in input_args:
        path = runtime_path(value)
        fingerprints.append({"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    manifest = {
        "candidates": str(runtime_path(candidates_arg)),
        "candidate_count": len(candidates),
        "evidence_count": len(evidence),
        "evidence": str(evidence_path),
        "candidate_summary": str(summary_path),
        "input_fingerprints": fingerprints,
        "scoring_principle": "Evidence is structured and source-backed. Database hit counts do not directly increase score.",
    }
    manifest_path = out_dir / "evidence_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest | {"manifest": str(manifest_path)}


def orthology_plan(
    ranked_arg: str,
    protein_arg: str,
    target_proteome_arg: str,
    relative_proteome_args: list[str],
    out_dir_arg: str,
    top_n: int = 7,
    threads: int = 4,
) -> dict[str, object]:
    if top_n < 1 or threads < 1:
        raise ValueError("top_n and threads must be at least 1")
    ranked = _candidate_frame(ranked_arg)
    if "dry_group" in ranked.columns:
        selected = ranked[ranked["dry_group"].isin({"top", "rescue"})].head(top_n)
    else:
        selected = ranked.head(top_n)
    selected_ids = set(selected["gene_id"])
    with runtime_path(protein_arg).open(encoding="utf-8") as handle:
        proteins = [record for record in SeqIO.parse(handle, "fasta") if record.id in selected_ids]
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    selected_path = out_dir / "selected_candidate_proteins.fa"
    SeqIO.write(proteins, str(selected_path), "fasta")
    missing = sorted(selected_ids.difference({record.id for record in proteins}))
    if missing:
        raise ValueError(f"selected candidate proteins are missing: {', '.join(missing)}")
    target = runtime_path(target_proteome_arg)
    lines = ["#!/usr/bin/env bash", "set -euo pipefail", f"OUT={shlex.quote(str(out_dir))}", "mkdir -p \"$OUT/db\" \"$OUT/results\"", f"makeblastdb -in {shlex.quote(str(target))} -dbtype prot -parse_seqids -out \"$OUT/db/target\" >/dev/null"]
    result_sets = []
    for specification in relative_proteome_args:
        if "=" not in specification:
            raise ValueError("--relative-proteome must use LABEL=FASTA")
        label, path_text = specification.split("=", 1)
        label = safe_id(label)
        relative = runtime_path(path_text)
        forward = out_dir / "results" / f"{label}.forward.tsv"
        reverse = out_dir / "results" / f"{label}.reverse.tsv"
        lines.extend(
            [
                f"makeblastdb -in {shlex.quote(str(relative))} -dbtype prot -parse_seqids -out \"$OUT/db/{label}\" >/dev/null",
                f"blastp -query {shlex.quote(str(selected_path))} -db \"$OUT/db/{label}\" -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle' -max_target_seqs 10 -num_threads {threads} -out {shlex.quote(str(forward))}",
                f"awk '!seen[$1]++ {{print $2}}' {shlex.quote(str(forward))} > \"$OUT/results/{label}.best_ids.txt\"",
                f"blastdbcmd -db \"$OUT/db/{label}\" -entry_batch \"$OUT/results/{label}.best_ids.txt\" -out \"$OUT/results/{label}.best_hits.fa\"",
                f"blastp -query \"$OUT/results/{label}.best_hits.fa\" -db \"$OUT/db/target\" -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle' -max_target_seqs 10 -num_threads {threads} -out {shlex.quote(str(reverse))}",
            ]
        )
        result_sets.append({"label": label, "forward": str(forward), "reverse": str(reverse)})
    script_path = out_dir / "run_orthology_refinement.sh"
    script_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    script_path.chmod(0o755)
    manifest = {"selected_genes": sorted(selected_ids), "script": str(script_path), "result_sets": result_sets, "method": "top-candidate reciprocal best hit with top-two ambiguity check"}
    manifest_path = out_dir / "orthology_plan.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest | {"manifest": str(manifest_path)}


def orthology_import(manifest_arg: str, out_arg: str) -> dict[str, object]:
    manifest_path = runtime_path(manifest_arg)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    rows: list[dict[str, str]] = []
    for result_set in manifest.get("result_sets", []):
        label = result_set["label"]
        forward_path = runtime_path(result_set["forward"])
        reverse_path = runtime_path(result_set["reverse"])
        forward = _best_hits(forward_path)
        reverse = _best_hits(reverse_path)
        for gene_id, hits in forward.items():
            if not hits:
                continue
            best = hits[0]
            reverse_best = reverse.get(str(best["subject"]), [])
            is_rbh = bool(reverse_best and str(reverse_best[0]["subject"]) == gene_id)
            rows.append(
                _row(
                    gene_id,
                    "rbh_ortholog" if is_rbh else "paralog_ambiguous",
                    "support" if is_rbh else "against",
                    "high" if is_rbh else "medium",
                    best["subject"],
                    "reciprocal_blast",
                    f"{label}:{forward_path}",
                    f"relative={label}; reciprocal_best_hit={is_rbh}; reverse_best={reverse_best[0]['subject'] if reverse_best else 'none'}",
                )
            )
            forward_ambiguous = len(hits) > 1 and float(hits[1]["bitscore"]) >= 0.90 * float(best["bitscore"])
            reverse_ambiguous = len(reverse_best) > 1 and float(reverse_best[1]["bitscore"]) >= 0.90 * float(reverse_best[0]["bitscore"])
            if is_rbh and (forward_ambiguous or reverse_ambiguous):
                rows.append(
                    _row(
                        gene_id,
                        "paralog_ambiguous",
                        "against",
                        "medium",
                        f"forward_top2_close={forward_ambiguous};reverse_top2_close={reverse_ambiguous}",
                        "reciprocal_blast",
                        f"{label}:{forward_path}",
                        "RBH is present, but a near-tied forward or reverse hit indicates unresolved paralog redundancy.",
                    )
                )
    output = runtime_path(out_arg)
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=EVIDENCE_COLUMNS).to_csv(output, sep="\t", index=False)
    return {"output": str(output), "evidence_count": len(rows), "source_manifest": str(manifest_path)}


