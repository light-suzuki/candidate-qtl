"""Candidate-first automatic comparative triage for every interval gene.

This route never requires a curated reference-gene panel.  All candidate
proteins are searched against role-labelled local reference proteomes; mapping
robustness, functional annotations, phenotype matching, variants, ranking, and
missing-data states remain separate and auditable.
"""

from __future__ import annotations

import json
import shutil
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd
import yaml
from Bio import SeqIO

from .core import ExternalDependencyError, file_fingerprint, runtime_path, run_tool, safe_id
from .homolog_shared import _run_interval_blastp
from .optional_deps import require_pysam
from .hardening import classify_source, prioritize_evidence
from .homolog_panel import _llm_hypotheses, _mechanical_queries, _overlap
from .homology_qc import passes_quality, quality_class
from .shortlist import screening_shortlist
from .workflow_foundation import WorkflowState, write_foundation_outputs
from .workflow_extensions import extract_multi


ORTHOLOGY_COLUMNS = [
    "candidate_gene_id", "reference_label", "reference_role", "reference_species",
    "reference_gene_id", "relationship", "orthology_confidence", "mapping_providers",
    "query_coverage", "subject_coverage", "identity", "evalue", "bitscore",
    "top2_bitscore_ratio", "rbh_status", "paralog_ambiguity", "domain_concordance",
    "synteny_concordance", "mapping_status", "failure_reason", "mapping_rank",
]
FUNCTION_COLUMNS = [
    "candidate_gene_id", "reference_label", "reference_gene_id", "gene_description", "function_description",
    "go_terms", "mutant_phenotype", "pathway", "expression", "protein_family",
    "source_id", "retrieval_status", "review_status",
]
EVIDENCE_COLUMNS = [
    "gene_id", "evidence_code", "direction", "confidence", "value",
    "source_type", "source_id", "description",
]
RICH_FUNCTION_FIELDS = (
    "function_description", "go_terms", "mutant_phenotype", "pathway", "expression", "protein_family",
)


# Lanes a run can traverse. The canonical core is Route A / Route B / A/B
# comparison / Route C (docs/CORE-CONTRACT.ja.md §2); `auto-triage` is the
# panel-free screening entry and is none of them.
CANONICAL_LANES = ("route_a", "route_b", "route_c", "ab_comparison")

# Which lane carries the recovery duty for each use-case type
# (docs/CORE-CONTRACT.ja.md §2.5). Type 4 is outside the Top-20 guarantee.
USE_CASE_TYPE_REQUIRED_LANES = {
    "type1_reference_or_id_mismatch": ("route_a",),
    "type2_semantic_or_mechanism_mismatch": ("route_a", "route_c"),
    "type3_cross_species_transferable": ("route_a", "route_b"),
}


def _lane_coverage(llm_engaged: bool) -> dict[str, object]:
    """Declare, at the top of the manifest, which lanes this run did NOT traverse.

    A run of this command completes and emits a Stage-A budget whether or not any
    canonical route ran, and the only signal that none did used to be a nested
    ``llm.used: false``. That is not enough to stop a reader from quoting the
    result as a use-case-type measurement, and in 2026-08 one was quoted that way
    (issue #493). The absence is therefore stated in full, next to the result, in
    the vocabulary the contract uses.
    """
    seams_engaged = ["trait_hypotheses"] if llm_engaged else []
    return {
        "lanes_executed": ["deterministic_rank"],
        "lanes_absent": list(CANONICAL_LANES),
        "llm_seams_engaged": seams_engaged,
        "llm_seams_available_not_engaged": [] if llm_engaged else ["trait_hypotheses"],
        "use_case_types_evaluable": [],
        "use_case_type_required_lanes": {
            key: list(value) for key, value in USE_CASE_TYPE_REQUIRED_LANES.items()
        },
        "reporting_rule": (
            "no canonical route lane ran in this command; its output is a panel-free screening "
            "budget and must not be reported as a use-case-type retention result "
            "(docs/CORE-CONTRACT.ja.md section 2.5)"
        ),
    }


def _resolved(path_value: str, panel_path: Path) -> Path:
    path = runtime_path(path_value)
    return path if path.is_absolute() else panel_path.parent / path


def _source_is_backed(value: object) -> bool:
    return classify_source(str(value).strip()) in {"reference", "local_file"}


def _load_resource_panel(panel_arg: str) -> tuple[dict[str, object], Path]:
    path = runtime_path(panel_arg)
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not isinstance(payload.get("references"), list):
        raise ValueError("resource panel requires a references list")
    allowed_roles = {"close_anchor", "functional_reference", "outgroup"}
    labels: set[str] = set()
    for item in payload["references"]:
        if not isinstance(item, dict) or not all(item.get(key) for key in ("label", "role", "species", "proteome")):
            raise ValueError("each reference requires label, role, species, and proteome")
        if item["role"] not in allowed_roles:
            raise ValueError(f"invalid reference role {item['role']!r}")
        if str(item["label"]) in labels:
            raise ValueError(f"duplicate reference label: {item['label']}")
        labels.add(str(item["label"]))
    return payload, path


def _read_hits(path: Path) -> pd.DataFrame:
    columns = ["query_id", "subject_id", "identity", "alignment_length", "query_length", "subject_length", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bitscore"]
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=columns)
    frame = pd.read_csv(path, sep="\t", header=None, names=columns, dtype=str).fillna("")
    for column in ("identity", "alignment_length", "query_length", "subject_length", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bitscore"):
        frame[column] = pd.to_numeric(frame[column], errors="raise")
    frame["query_coverage"] = ((frame["query_end"] - frame["query_start"]).abs() + 1) / frame["query_length"].clip(lower=1)
    frame["subject_coverage"] = ((frame["subject_end"] - frame["subject_start"]).abs() + 1) / frame["subject_length"].clip(lower=1)
    frame["quality_class"] = [quality_class(row.identity, row.query_coverage, row.subject_coverage) for row in frame.itertuples()]
    frame["quality_pass"] = [passes_quality(row.identity, row.query_coverage, row.evalue) for row in frame.itertuples()]
    return frame.sort_values(["query_id", "bitscore", "evalue", "subject_id"], ascending=[True, False, True, True], kind="stable")


def _reference_annotations(reference: dict[str, object], panel_path: Path, descriptions: dict[str, str]) -> dict[str, dict[str, str]]:
    records = {
        key: {
            "reference_gene_id": key, "gene_description": value, "function_description": "", "go_terms": "",
            "mutant_phenotype": "", "pathway": "", "expression": "", "protein_family": "",
            "source_id": f"local_fasta:{reference['label']}", "retrieval_status": "fasta_description_only",
        }
        for key, value in descriptions.items()
    }
    annotation_value = str(reference.get("annotation", "")).strip()
    if not annotation_value:
        return records
    annotation_path = _resolved(annotation_value, panel_path)
    frame = pd.read_csv(annotation_path, sep="\t", dtype=str).fillna("")
    if "reference_gene_id" not in frame.columns:
        raise ValueError(f"{annotation_path}: annotation requires reference_gene_id")
    for row in frame.to_dict(orient="records"):
        gene_id = str(row["reference_gene_id"])
        aliases = [gene_id]
        aliases.extend(key for key in descriptions if key != gene_id and gene_id in key.split("|"))
        expression = str(row.get("expression", "")).strip()
        if not expression:
            expression_parts = []
            tissue = str(row.get("tissue_specificity", "")).strip()
            stage = str(row.get("developmental_stage", "")).strip()
            if tissue:
                expression_parts.append(f"tissue_specificity: {tissue}")
            if stage:
                expression_parts.append(f"developmental_stage: {stage}")
            expression = "; ".join(expression_parts)
        adapted = {
            "gene_description": str(row.get("gene_description", "") or row.get("description", "")),
            "function_description": str(row.get("function_description", "")),
            "go_terms": str(row.get("go_terms", "")),
            "mutant_phenotype": str(row.get("mutant_phenotype", "")),
            "pathway": str(row.get("pathway", "")),
            "expression": expression,
            "protein_family": str(row.get("protein_family", "")),
            "source_id": str(row.get("source_id", "")),
        }
        retrieval_status = "source_retrieved" if _source_is_backed(adapted["source_id"]) else "annotation_loaded_source_unverified"
        for alias in aliases:
            existing = records.setdefault(alias, {"reference_gene_id": alias})
            for column, value in adapted.items():
                if value or column not in existing:
                    existing[column] = value or str(existing.get(column, ""))
            existing["reference_gene_id"] = alias
            existing["retrieval_status"] = retrieval_status
    return records


def _functional_annotation_coverage(functions: pd.DataFrame) -> dict[str, object]:
    rows = int(len(functions))
    field_coverage: dict[str, dict[str, object]] = {}
    for field in ("gene_description", *RICH_FUNCTION_FIELDS):
        values = functions[field].astype(str).str.strip() if field in functions.columns else pd.Series(dtype=str)
        non_empty = int(values.ne("").sum()) if rows else 0
        field_coverage[field] = {
            "non_empty": non_empty,
            "rows": rows,
            "coverage": round(non_empty / rows, 6) if rows else 0.0,
        }
    retrieval_counts = (
        {str(key): int(value) for key, value in functions["retrieval_status"].value_counts().sort_index().items()}
        if rows and "retrieval_status" in functions.columns else {}
    )
    if rows:
        rich_mask = pd.Series(False, index=functions.index)
        for field in RICH_FUNCTION_FIELDS:
            if field in functions.columns:
                rich_mask |= functions[field].astype(str).str.strip().ne("")
        retrieval_mask = functions["retrieval_status"].eq("source_retrieved") if "retrieval_status" in functions.columns else pd.Series(False, index=functions.index)
        source_mask = functions["source_id"].map(_source_is_backed) if "source_id" in functions.columns else pd.Series(False, index=functions.index)
        source_backed_rich_rows = int((rich_mask & retrieval_mask & source_mask).sum())
    else:
        source_backed_rich_rows = 0
    return {
        "rows": rows,
        "field_coverage": field_coverage,
        "retrieval_status_counts": retrieval_counts,
        "source_backed_rich_rows": source_backed_rich_rows,
        "source_backed_rich_fraction": round(source_backed_rich_rows / rows, 6) if rows else 0.0,
        "status": "source_backed_function_available" if source_backed_rich_rows else "insufficient_source_backed_function",
    }


def _automatic_context_evidence(candidates: pd.DataFrame, orthology: pd.DataFrame, candidates_path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for gene_id in candidates["gene_id"].astype(str):
        rows.append({
            "gene_id": gene_id,
            "evidence_code": "interval_membership",
            "direction": "context",
            "confidence": "high",
            "value": "qtl_interval_candidate",
            "source_type": "candidate_universe",
            "source_id": str(candidates_path),
            "description": "Candidate is inside the supplied QTL interval; interval membership defines scope and does not score one gene over another.",
        })
    if orthology.empty:
        return rows
    for record in orthology[orthology["mapping_status"].eq("mapped")].to_dict(orient="records"):
        gene_id = str(record["candidate_gene_id"])
        label = str(record.get("reference_label", ""))
        reference_gene = str(record.get("reference_gene_id", ""))
        confidence = str(record.get("orthology_confidence", "low"))
        if confidence not in {"low", "medium", "high"}:
            confidence = "low"
        value = ";".join([
            f"reference={label}:{reference_gene}",
            f"relationship={record.get('relationship', '')}",
            f"rbh={record.get('rbh_status', '')}",
            f"identity={record.get('identity', '')}",
            f"query_coverage={record.get('query_coverage', '')}",
            f"subject_coverage={record.get('subject_coverage', '')}",
        ])
        rows.append({
            "gene_id": gene_id,
            "evidence_code": "homology_context",
            "direction": "context",
            "confidence": confidence,
            "value": value,
            "source_type": "candidate_first_reference_mapping",
            "source_id": f"reference:{label}:{reference_gene}",
            "description": "Quality-passing generic sequence/orthology mapping. It records conservation and mapping robustness but is not trait-specific support.",
        })
        if str(record.get("paralog_ambiguity", "")) == "top2_near_tie":
            rows.append({
                "gene_id": gene_id,
                "evidence_code": "paralog_ambiguous",
                "direction": "context",
                "confidence": "medium",
                "value": f"top2_bitscore_ratio={record.get('top2_bitscore_ratio', '')}",
                "source_type": "candidate_first_reference_mapping",
                "source_id": f"reference:{label}:{reference_gene}",
                "description": "Near-tied paralogs make copy-level interpretation ambiguous; this is epistemic context, not biological counterevidence.",
            })
    return rows


def _trait_match_evidence(row: object, matches: list[str]) -> dict[str, str]:
    """Convert a lexical trait/function match without inventing provenance.

    A product name or local FASTA description can be useful retrieval context,
    but it is not source-backed trait support. Only annotation rows that were
    resolved to a citable/local source are allowed to emit the scored
    ``automated_candidate_function_match`` code. All other lexical matches stay
    visible as unscored ``trait_context`` rather than disappearing.
    """
    source_id = str(getattr(row, "source_id", ""))
    retrieval_status = str(getattr(row, "retrieval_status", ""))
    source_backed = retrieval_status == "source_retrieved" and _source_is_backed(source_id)
    return {
        "gene_id": str(getattr(row, "candidate_gene_id")),
        "evidence_code": "automated_candidate_function_match" if source_backed else "trait_context",
        "direction": "support" if source_backed else "context",
        "confidence": "low",
        "value": ";".join(matches),
        "source_type": "candidate_first_reference_annotation",
        "source_id": source_id,
        "description": (
            f"Source-backed functional-reference text matched phenotype search terms: {', '.join(matches)}."
            if source_backed else
            f"Lexical match in FASTA-only or source-unverified annotation: {', '.join(matches)}; retained as context only."
        ),
    }


def _supporting_source_summary(evidence: pd.DataFrame) -> dict[str, str]:
    if evidence.empty or "direction" not in evidence.columns or "source_id" not in evidence.columns:
        return {}
    support = evidence[
        evidence["direction"].astype(str).eq("support")
        & evidence["source_id"].astype(str).str.strip().ne("")
    ]
    if support.empty:
        return {}
    return {
        str(gene_id): ";".join(sorted(set(group["source_id"].astype(str))))
        for gene_id, group in support.groupby("gene_id", sort=True)
    }


def _rank_semantics(evidence: pd.DataFrame) -> dict[str, object]:
    if evidence.empty:
        return {
            "support_evidence_count": 0,
            "against_evidence_count": 0,
            "context_evidence_count": 0,
            "rank_metric_eligible": False,
            "stage_a_is_budget_only": True,
            "rank_semantic_status": "not_evaluated_no_scoring_evidence",
            "rank_lane_warning": "no scoring evidence rows were available; Stage-A is a budget-only retrieval set",
        }
    directions = evidence["direction"].astype(str)
    support_count = int(directions.eq("support").sum())
    against_count = int(directions.eq("against").sum())
    context_count = int(directions.eq("context").sum())
    base = {
        "support_evidence_count": support_count,
        "against_evidence_count": against_count,
        "context_evidence_count": context_count,
    }
    if support_count == 0:
        return {
            **base,
            "rank_metric_eligible": False,
            "rank_is_penalty_only": True,
            "stage_a_is_budget_only": True,
            "rank_semantic_status": "not_evaluated_no_supporting_evidence",
            "rank_lane_warning": "evidence rows exist but none provide supporting evidence; Stage-A is budget-only and tied order must not be read as biological priority",
        }
    return {
        **base,
        "rank_metric_eligible": True,
        "stage_a_is_budget_only": False,
        "rank_semantic_status": "evaluable",
        "rank_lane_warning": "supporting evidence is present; rank remains an assay-priority heuristic, not causal probability",
    }


def _reciprocal_best_hits(
    top_subjects: set[str], reference_fasta: Path, native_fasta: Path | None,
    out_dir: Path, runner: str, threads: int, resume: bool = False,
    workflow_state: WorkflowState | None = None, cache_key: str = "reciprocal", provider_version: str = "",
) -> tuple[dict[str, str], str]:

    if not native_fasta:
        return {}, "skipped_native_proteome_unavailable"
    index = SeqIO.index(str(reference_fasta), "fasta")
    records = [index[value] for value in sorted(top_subjects) if value in index]
    index.close()
    if not records:
        return {}, "skipped_reference_ids_not_in_fasta"
    query = out_dir / "reciprocal_queries.faa"
    reciprocal_path = out_dir / "reciprocal_vs_native.tsv"
    signature = workflow_state.signature(
        [reference_fasta, native_fasta],
        {
            "step": "reciprocal_blastp_v1", "subjects": sorted(top_subjects), "runner": runner,
            "threads": threads, "provider_version": provider_version,
        },
    ) if workflow_state else ""
    if resume and workflow_state and workflow_state.can_reuse(cache_key, signature, [reciprocal_path]):
        result = {"output": str(reciprocal_path)}
        run_status = "resumed"
        workflow_state.record(cache_key, signature, "complete", [reciprocal_path], reused=True)
    else:
        SeqIO.write(records, str(query), "fasta")
        try:
            result = _run_interval_blastp(query, native_fasta, reciprocal_path, runner, threads)
            run_status = "completed"
            if workflow_state:
                workflow_state.record(cache_key, signature, "complete", [reciprocal_path])
        except Exception as error:
            if workflow_state:
                workflow_state.record(
                    cache_key, signature, "failed", failure_code="execution_failed", failure_message=str(error)
                )
            raise
    hits = _read_hits(runtime_path(str(result["output"])))
    if hits.empty:
        return {}, "completed_no_hits"
    top = hits.drop_duplicates("query_id")
    return dict(zip(top["query_id"].astype(str), top["subject_id"].astype(str))), run_status


def _variant_rows(parent_vcf: str | None, candidates: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    columns = ["gene_id", "chrom", "pos", "ref", "alt", "parent_genotypes", "genotype_differs", "effect_annotation", "review_status"]
    if not parent_vcf:
        return pd.DataFrame(columns=columns), "not_supplied"
    path = runtime_path(parent_vcf)
    intervals = [
        (str(row.gene_id), str(row.chrom), int(row.start), int(row.end))
        for row in candidates.itertuples()
        if all(hasattr(row, key) for key in ("chrom", "start", "end"))
    ]
    if not intervals:
        return pd.DataFrame(columns=columns), "skipped_candidate_coordinates_unavailable"
    rows: list[dict[str, object]] = []
    pysam = require_pysam('auto-triage with --parent-vcf')
    with pysam.VariantFile(str(path)) as handle:
        samples = list(handle.header.samples)
        for record in handle:
            for gene_id, chrom, start, end in intervals:
                if str(record.contig) != chrom or not start <= int(record.pos) <= end:
                    continue
                genotypes = {sample: "/".join("." if value is None else str(value) for value in (record.samples[sample].get("GT") or ())) for sample in samples}
                distinct = {value for value in genotypes.values() if value and value != "./."}
                effect = ""
                for key in ("ANN", "CSQ"):
                    if key in record.info:
                        effect = str(record.info[key])
                        break
                rows.append({
                    "gene_id": gene_id, "chrom": chrom, "pos": record.pos, "ref": record.ref,
                    "alt": ",".join(record.alts or ()), "parent_genotypes": ";".join(f"{key}={value}" for key, value in genotypes.items()),
                    "genotype_differs": len(distinct) > 1, "effect_annotation": effect,
                    "review_status": "unreviewed",
                })
    return pd.DataFrame(rows, columns=columns), "completed"


def _tool_version(runner: str) -> str:
    selected = runner
    if runner == "auto":
        selected = "native" if shutil.which("blastp") else "wsl"
    command = ["blastp", "-version"] if selected != "wsl" else ["wsl", "blastp", "-version"]
    try:
        completed = run_tool(command, capture_output=True, text=True, check=False)
    except (OSError, ExternalDependencyError):
        return "version_unavailable"
    return completed.stdout.splitlines()[0] if completed.returncode == 0 and completed.stdout else "version_unavailable"


def auto_triage(
    phenotype: str,
    resource_panel_arg: str,
    out_dir_arg: str,
    candidates_arg: str | None = None,
    candidate_proteins_arg: str | None = None,
    fasta_arg: str | None = None,
    gff3_arg: str | None = None,
    interval: str | None = None,
    base_evidence_arg: str | None = None,
    parent_vcf: str | None = None,
    trait_terms: list[str] | None = None,
    llm_provider: str | None = None,
    model: str | None = None,
    endpoint: str | None = None,
    allow_external_llm: bool = False,
    runner: str = "auto",
    threads: int = 4,
    target_size: int = 20,
    resume: bool = False,
    id_map_arg: str | None = None,
) -> dict[str, object]:
    """Run panel-free comparative screening for every QTL-interval candidate."""

    if not phenotype.strip():
        raise ValueError("phenotype is required")
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    extraction: dict[str, object] | None = None
    if candidates_arg or candidate_proteins_arg:
        if not candidates_arg or not candidate_proteins_arg:
            raise ValueError("--candidates and --candidate-proteins must be supplied together")
        candidates_path = runtime_path(candidates_arg)
        proteins_path = runtime_path(candidate_proteins_arg)
    else:
        if not all((fasta_arg, gff3_arg, interval)):
            raise ValueError("supply candidates/proteins, or fasta/gff3/interval")
        extraction = extract_multi(str(fasta_arg), str(gff3_arg), str(interval), str(out_dir / "01_extract"))
        candidates_path = runtime_path(str(extraction["outputs"]["candidates_tsv"]))
        proteins_path = runtime_path(str(extraction["outputs"]["candidate_proteins_fasta"]))
    inferred_id_map = runtime_path(id_map_arg) if id_map_arg else candidates_path.parent / "id_map.tsv"
    id_map = pd.read_csv(inferred_id_map, sep="\t", dtype=str).fillna("") if inferred_id_map.exists() else None
    candidates = pd.read_csv(candidates_path, sep="\t", dtype=str).fillna("")
    if "gene_id" not in candidates.columns or candidates["gene_id"].duplicated().any():
        raise ValueError("candidate table requires unique gene_id values")
    candidate_ids = set(candidates["gene_id"])
    with proteins_path.open(encoding="utf-8") as handle:
        protein_ids = {record.id for record in SeqIO.parse(handle, "fasta")}
    panel, panel_path = _load_resource_panel(resource_panel_arg)
    native = panel.get("native", {}) if isinstance(panel.get("native"), dict) else {}
    native_value = str(native.get("proteome", "")).strip()
    native_fasta = _resolved(native_value, panel_path) if native_value else None
    if native_fasta and not native_fasta.exists():
        native_fasta = None

    hypotheses, llm_meta = _llm_hypotheses(
        phenotype, str(native.get("species", "target crop")), "unknown", "unknown", out_dir / "trait_hypotheses",
        llm_provider, model, endpoint, allow_external_llm, resume,
    )
    search_rows = _mechanical_queries(phenotype)
    for phrase in trait_terms or []:
        search_rows.extend(_mechanical_queries(phrase))
    for row in hypotheses:
        search_rows.append({"query": row["query"], "origin": "llm", "mechanism": row["mechanism"], "known_gene": row["known_gene"], "rationale": row["rationale"]})
        if row["mechanism"]:
            search_rows.append({"query": row["mechanism"], "origin": "llm_mechanism", "mechanism": row["mechanism"], "known_gene": row["known_gene"], "rationale": row["rationale"]})
    search_terms = list(dict.fromkeys(str(row["query"]) for row in search_rows))
    pd.DataFrame(search_rows).to_csv(out_dir / "trait_search_hypotheses.tsv", sep="\t", index=False)

    orthology_rows: list[dict[str, object]] = []
    function_rows: list[dict[str, str]] = []
    provider_runs: list[dict[str, object]] = []
    workflow_state = WorkflowState(out_dir / "workflow_state.json")
    provider_version = _tool_version(runner)
    for reference in panel["references"]:
        label = str(reference["label"])
        role = str(reference["role"])
        species = str(reference["species"])
        reference_fasta = _resolved(str(reference["proteome"]), panel_path)
        reference["_resolved_proteome"] = str(reference_fasta)
        provider_dir = out_dir / "providers" / safe_id(label)
        provider_dir.mkdir(parents=True, exist_ok=True)
        started_at = datetime.now(UTC).isoformat()
        try:
            with reference_fasta.open(encoding="utf-8") as handle:
                descriptions = {record.id: record.description for record in SeqIO.parse(handle, "fasta")}
            forward_path = provider_dir / "candidate_vs_reference.tsv"
            forward_key = f"provider:{safe_id(label)}:forward_blastp"
            forward_signature = workflow_state.signature(
                [proteins_path, reference_fasta],
                {
                    "step": "candidate_vs_reference_blastp_v1", "runner": runner, "threads": threads,
                    "provider_version": provider_version,
                },
            )
            if resume and workflow_state.can_reuse(forward_key, forward_signature, [forward_path]):
                blast = {"output": str(forward_path), "runner": "resumed", "command": []}
                workflow_state.record(forward_key, forward_signature, "complete", [forward_path], reused=True)
            else:
                try:
                    blast = _run_interval_blastp(proteins_path, reference_fasta, forward_path, runner, threads)
                    workflow_state.record(forward_key, forward_signature, "complete", [forward_path])
                except Exception as error:
                    workflow_state.record(
                        forward_key, forward_signature, "failed", failure_code="execution_failed",
                        failure_message=str(error),
                    )
                    raise
            hits = _read_hits(runtime_path(str(blast["output"])))
            top_subjects = set(hits[hits["quality_pass"]].drop_duplicates("query_id")["subject_id"].astype(str))
            reciprocal, reciprocal_status = _reciprocal_best_hits(
                top_subjects, reference_fasta, native_fasta, provider_dir, runner, threads, resume,
                workflow_state, f"provider:{safe_id(label)}:reciprocal_blastp", provider_version,
            )
            annotation_status = "complete"
            annotation_error = ""
            try:
                annotations = _reference_annotations(reference, panel_path, descriptions)
            except Exception as error:
                annotation_status = "failed"
                annotation_error = str(error)
                annotations = _reference_annotations({**reference, "annotation": ""}, panel_path, descriptions)
            provider_runs.append({
                "label": label, "role": role, "status": "completed", "blast": blast,
                "reciprocal_status": reciprocal_status, "annotation_status": annotation_status,
                "annotation_error": annotation_error, "started_at": started_at,
                "finished_at": datetime.now(UTC).isoformat(),
            })
            for gene_id in sorted(candidate_ids):
                if gene_id not in protein_ids:
                    orthology_rows.append({"candidate_gene_id": gene_id, "reference_label": label, "reference_role": role, "reference_species": species, "mapping_status": "unmapped", "relationship": "unmapped", "failure_reason": "missing_candidate_protein", "mapping_rank": ""})
                    continue
                gene_hits = hits[(hits["query_id"] == gene_id) & hits["quality_pass"]].head(5)
                if gene_hits.empty:
                    orthology_rows.append({"candidate_gene_id": gene_id, "reference_label": label, "reference_role": role, "reference_species": species, "mapping_status": "unmapped", "relationship": "unmapped", "failure_reason": "no_quality_passing_hit", "mapping_rank": ""})
                    continue
                top = gene_hits.iloc[0]
                ratio = float(gene_hits.iloc[1]["bitscore"]) / float(top["bitscore"]) if len(gene_hits) > 1 and float(top["bitscore"]) else 0.0
                ambiguous = ratio >= 0.90
                selected_hits = gene_hits[gene_hits["bitscore"] >= 0.90 * float(top["bitscore"])] if ambiguous else gene_hits.head(1)
                for mapping_rank, (_, hit) in enumerate(selected_hits.iterrows(), start=1):
                    rbh = reciprocal.get(str(hit["subject_id"])) == gene_id
                    relationship = "one_to_many" if len(selected_hits) > 1 else "one_to_one_candidate" if rbh else "mapped_non_rbh"
                    confidence = "high" if rbh and not ambiguous and hit["query_coverage"] >= 0.7 and hit["subject_coverage"] >= 0.7 else "medium" if not ambiguous else "low"
                    orthology_rows.append({
                        "candidate_gene_id": gene_id, "reference_label": label, "reference_role": role,
                        "reference_species": species, "reference_gene_id": str(hit["subject_id"]),
                        "relationship": relationship, "orthology_confidence": confidence,
                        "mapping_providers": "blastp;reciprocal_blastp" if native_fasta else "blastp",
                        "query_coverage": round(float(hit["query_coverage"]), 6), "subject_coverage": round(float(hit["subject_coverage"]), 6),
                        "identity": round(float(hit["identity"]), 3), "evalue": hit["evalue"], "bitscore": hit["bitscore"],
                        "top2_bitscore_ratio": round(ratio, 6), "rbh_status": "yes" if rbh else reciprocal_status if not native_fasta else "no",
                        "paralog_ambiguity": "top2_near_tie" if ambiguous else "none_detected",
                        "domain_concordance": "not_run", "synteny_concordance": "not_run",
                        "mapping_status": "mapped", "failure_reason": "", "mapping_rank": mapping_rank,
                    })
                    annotation = annotations.get(str(hit["subject_id"]), {"reference_gene_id": str(hit["subject_id"]), "retrieval_status": "annotation_unavailable" if annotation_status == "complete" else "annotation_provider_failed"})
                    function_rows.append({
                        "candidate_gene_id": gene_id, "reference_label": label,
                        **{column: str(annotation.get(column, "")) for column in FUNCTION_COLUMNS if column not in {"candidate_gene_id", "reference_label", "review_status"}},
                        "review_status": "unreviewed",
                    })
        except Exception as error:
            provider_runs.append({"label": label, "role": role, "status": "failed", "error": str(error), "started_at": started_at, "finished_at": datetime.now(UTC).isoformat()})
            for gene_id in sorted(candidate_ids):
                orthology_rows.append({"candidate_gene_id": gene_id, "reference_label": label, "reference_role": role, "reference_species": species, "mapping_status": "provider_failure", "relationship": "unmapped", "failure_reason": str(error)[:300], "mapping_rank": ""})

    orthology = pd.DataFrame(orthology_rows, columns=ORTHOLOGY_COLUMNS).fillna("")
    orthology_path = out_dir / "02_automatic_orthologs.tsv"
    orthology.to_csv(orthology_path, sep="\t", index=False)
    functions = pd.DataFrame(function_rows, columns=FUNCTION_COLUMNS).fillna("")
    functions_path = out_dir / "03_reference_functions.tsv"
    functions.to_csv(functions_path, sep="\t", index=False)
    annotation_coverage = _functional_annotation_coverage(functions)

    trait_rows: list[dict[str, object]] = []
    evidence_rows: list[dict[str, str]] = _automatic_context_evidence(candidates, orthology, candidates_path)
    for row in functions.itertuples():
        functional_text = " ".join(str(getattr(row, column)) for column in ("gene_description", "function_description", "go_terms", "mutant_phenotype", "pathway", "expression", "protein_family"))
        matches = [term for term in search_terms if _overlap(term, functional_text)]
        if not matches:
            continue
        trait_rows.append({"gene_id": row.candidate_gene_id, "reference_label": row.reference_label, "reference_gene_id": row.reference_gene_id, "matched_terms": ";".join(matches), "functional_text": functional_text, "source_id": row.source_id, "review_status": "unreviewed"})
        evidence_rows.append(_trait_match_evidence(row, matches))
    trait_match = pd.DataFrame(trait_rows, columns=["gene_id", "reference_label", "reference_gene_id", "matched_terms", "functional_text", "source_id", "review_status"])
    trait_match_path = out_dir / "04_trait_function_matches.tsv"
    trait_match.to_csv(trait_match_path, sep="\t", index=False)

    variants, variant_status = _variant_rows(parent_vcf, candidates)
    variants_path = out_dir / "05_parental_variants.tsv"
    variants.to_csv(variants_path, sep="\t", index=False)
    automatic_evidence = pd.DataFrame(evidence_rows, columns=EVIDENCE_COLUMNS).drop_duplicates()
    if base_evidence_arg:
        base = pd.read_csv(runtime_path(base_evidence_arg), sep="\t", dtype=str).fillna("")
        combined = pd.concat([base, automatic_evidence], ignore_index=True).drop_duplicates()
    else:
        combined = automatic_evidence
    evidence_path = out_dir / "automatic_screening_evidence.tsv"
    combined.to_csv(evidence_path, sep="\t", index=False)
    supporting_source_summary = _supporting_source_summary(combined)
    automatic_source_backed_trait_support_gene_count = int(
        automatic_evidence[
            automatic_evidence["evidence_code"].eq("automated_candidate_function_match")
            & automatic_evidence["direction"].eq("support")
        ]["gene_id"].nunique()
    )
    rank_semantics = _rank_semantics(combined)
    rank = prioritize_evidence(str(candidates_path), str(evidence_path), str(out_dir / "rank"))
    ranked = pd.read_csv(runtime_path(str(rank["ranked_candidates"])), sep="\t", dtype=str).fillna("")
    ranked["_screening_order"] = pd.to_numeric(ranked["screening_rank"], errors="raise")
    ranked["_score_order"] = pd.to_numeric(ranked["dry_score"], errors="raise")
    ranked = ranked.sort_values(["_screening_order", "_score_order", "gene_id"], ascending=[True, False, True], kind="stable").drop(columns=["_screening_order", "_score_order"])
    mapped_orthology = orthology[orthology["mapping_status"].eq("mapped")].copy()
    mapped_orthology["_confidence_order"] = mapped_orthology["orthology_confidence"].map({"high": 0, "medium": 1, "low": 2}).fillna(3)
    mapped_orthology["_bitscore"] = pd.to_numeric(mapped_orthology["bitscore"], errors="coerce").fillna(0.0)
    best_orthology = mapped_orthology.sort_values(["_confidence_order", "_bitscore"], ascending=[True, False], kind="stable").drop_duplicates("candidate_gene_id").set_index("candidate_gene_id")
    function_summary = trait_match.groupby("gene_id")["matched_terms"].apply(lambda values: ";".join(sorted(set(";".join(values).split(";"))))).to_dict() if not trait_match.empty else {}
    variant_summary = variants.groupby("gene_id").size().to_dict() if not variants.empty else {}
    failure_summary = orthology[orthology["mapping_status"].ne("mapped")].groupby("candidate_gene_id")["failure_reason"].apply(lambda values: ";".join(sorted(set(filter(None, values))))).to_dict()
    ranked["automated_rank"] = range(1, len(ranked) + 1)
    ranked["automated_score"] = ranked["dry_score"]
    ranked["orthology_summary"] = ranked["gene_id"].map(lambda gene: str(best_orthology.loc[gene, "relationship"]) if gene in best_orthology.index else "unmapped")
    ranked["best_functional_reference"] = ranked["gene_id"].map(lambda gene: str(best_orthology.loc[gene, "reference_gene_id"]) if gene in best_orthology.index else "")
    ranked["trait_match_summary"] = ranked["gene_id"].map(function_summary).fillna("")
    ranked["parental_variant_summary"] = ranked["gene_id"].map(lambda gene: f"variant_count={variant_summary.get(gene, 0)}" if parent_vcf else "not_supplied")
    ranked["paralog_ambiguity"] = ranked["gene_id"].map(lambda gene: str(best_orthology.loc[gene, "paralog_ambiguity"]) if gene in best_orthology.index else "unknown")
    ranked["review_status"] = "unreviewed"
    ranked["supporting_source_ids"] = ranked["gene_id"].map(supporting_source_summary).fillna("")
    ranked["failure_or_missing_data_flags"] = ranked["gene_id"].map(failure_summary).fillna("")
    ranked_path = out_dir / "06_ranked_candidates.tsv"
    ranked.to_csv(ranked_path, sep="\t", index=False)
    shortlist = screening_shortlist(str(ranked_path), str(out_dir / "stage_a"), target_size=target_size)
    shortlist_manifest_path = runtime_path(str(shortlist["manifest"]))
    shortlist_manifest = json.loads(shortlist_manifest_path.read_text(encoding="utf-8"))
    shutil.copyfile(runtime_path(str(shortlist["shortlist"])), out_dir / "07_top_candidates_review.tsv")
    top = pd.read_csv(runtime_path(str(shortlist["shortlist"])), sep="\t", dtype=str).fillna("")
    assay = pd.DataFrame({
        "gene_id": top["gene_id"],
        "recommended_next_step": ["review_parental_variant_then_design_assay" if variant_summary.get(gene, 0) else "parental_sequence_then_recombinant_assay" for gene in top["gene_id"]],
        "review_status": "not_started",
    })
    assay_path = out_dir / "08_assay_plan.tsv"
    assay.to_csv(assay_path, sep="\t", index=False)

    mapped_counts = orthology.groupby("candidate_gene_id")["mapping_status"].apply(lambda values: any(value == "mapped" for value in values))
    foundation_outputs = write_foundation_outputs(
        out_dir=out_dir, candidates=candidates, protein_ids=protein_ids, panel=panel, panel_path=panel_path,
        orthology=orthology, functions=functions, trait_matches=trait_match, variants=variants,
        evidence=combined, ranked=ranked, provider_runs=provider_runs, provider_version=provider_version,
        inputs={
            "candidates": candidates_path, "candidate_proteins": proteins_path,
            "id_map": inferred_id_map if inferred_id_map.exists() else None,
            "base_evidence": runtime_path(base_evidence_arg) if base_evidence_arg else None,
            "parent_vcf": runtime_path(parent_vcf) if parent_vcf else None,
        },
        mode="llm-assisted" if llm_provider else "local",
        id_map=id_map,
    )
    manifest = {
        "schema_version": "1.0", "route": "candidate_first_auto_triage",
        "resume_requested": resume,
        "phenotype": phenotype, "candidate_count": len(candidates),
        "candidate_with_any_mapping_count": int(mapped_counts.sum()),
        "candidate_without_mapping_count": int((~mapped_counts).sum()),
        "resource_panel": file_fingerprint(panel_path, False),
        "reference_resources": [
            {
                "label": str(item["label"]), "role": str(item["role"]), "species": str(item["species"]),
                "proteome": (
                    file_fingerprint(_resolved(str(item["proteome"]), panel_path), False)
                    if _resolved(str(item["proteome"]), panel_path).exists()
                    else {"path": str(_resolved(str(item["proteome"]), panel_path)), "status": "missing"}
                ),
                "annotation": (
                    file_fingerprint(_resolved(str(item["annotation"]), panel_path), False)
                    if str(item.get("annotation", "")).strip() and _resolved(str(item["annotation"]), panel_path).exists()
                    else None
                ),
            }
            for item in panel["references"]
        ],
        "inputs": {"candidates": file_fingerprint(candidates_path, False), "candidate_proteins": file_fingerprint(proteins_path, False), "parent_vcf": file_fingerprint(runtime_path(parent_vcf), False) if parent_vcf else None},
        "provider": {"tool": "BLASTP", "version": provider_version, "runs": provider_runs},
        "llm": llm_meta, "llm_role": "trait_search_hypothesis_generation_only", "llm_writes_gene_score": False,
        "lane_coverage": _lane_coverage(bool(llm_meta.get("used"))),
        "parental_variant_status": variant_status,
        "functional_annotation_coverage": annotation_coverage,
        "trait_match_semantics": {
            "lexical_match_gene_count": int(trait_match["gene_id"].nunique()) if not trait_match.empty else 0,
            "source_backed_support_gene_count": automatic_source_backed_trait_support_gene_count,
            "unverified_matches_are_context_only": True,
        },
        "rank_semantics": {
            **rank_semantics,
            "boundary_tie_count": int(shortlist_manifest.get("boundary_tie_count", 0)),
            "boundary_tie_truncated": bool(shortlist_manifest.get("boundary_tie_truncated", False)),
            "non_biological_tie_break_disclosed": True,
        },
        "review_status": "unreviewed", "causal_probability": False,
        "independence_rule": "Mapping providers and annotation channels from the same reference chain form one comparative-function family.",
        "outputs": {"orthologs": str(orthology_path), "functions": str(functions_path), "trait_matches": str(trait_match_path), "variants": str(variants_path), "ranked": str(ranked_path), "top_review": str(out_dir / "07_top_candidates_review.tsv"), "assay_plan": str(assay_path), "workflow_state": str(workflow_state.path), **foundation_outputs},
        "created_at": datetime.now(UTC).isoformat(),
    }
    manifest_path = out_dir / "auto_triage_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    report = out_dir / "report.md"
    rank_note = (
        "Stage-A is a budget-only fallback because no supporting evidence was available; boundary ties are non-biological and must not be reported as a gene-priority rank."
        if rank_semantics["stage_a_is_budget_only"] else
        "Stage-A contains supporting evidence but remains an assay-priority heuristic, not a causal probability."
    )
    report.write_text(
        "# Candidate-first automatic triage\n\n"
        f"- Candidates: {len(candidates)}\n- Any reference mapping: {int(mapped_counts.sum())}\n"
        f"- Trait-function lexical matches: {trait_match['gene_id'].nunique() if not trait_match.empty else 0}\n"
        f"- Source-backed trait-support genes: {automatic_source_backed_trait_support_gene_count}\n"
        f"- Source-backed rich annotation rows: {annotation_coverage['source_backed_rich_rows']} / {annotation_coverage['rows']}\n"
        f"- Rank semantic status: {rank_semantics['rank_semantic_status']}\n"
        f"- Fixed Stage-A review set: {len(top)}\n- Parent VCF: {variant_status}\n\n"
        f"{rank_note}\n\nAll automatically retrieved functions remain unreviewed.\n",
        encoding="utf-8",
    )
    return {
        "ranked": str(ranked_path), "top_review": str(out_dir / "07_top_candidates_review.tsv"),
        "orthologs": str(orthology_path), "functions": str(functions_path), "manifest": str(manifest_path),
        "report": str(report), "workflow_state": str(workflow_state.path), "candidate_count": len(candidates),
        "top_count": len(top), **foundation_outputs,
    }
