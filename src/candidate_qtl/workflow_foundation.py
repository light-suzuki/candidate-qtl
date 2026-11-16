"""Canonical identity, provider-state, lineage, and cache primitives.

The workflow deliberately keeps these records separate from biological rank.
They describe what was mapped, what ran, what failed, and which observations a
derived result depends on; they are not extra evidence and are never scored.
"""

from __future__ import annotations

import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Iterable

import pandas as pd

from .core import file_fingerprint


FOUNDATION_SCHEMA_VERSION = "1.0"
COVERAGE_STATES = {
    "not_requested", "not_run", "running", "complete", "partial", "no_hit", "unmapped",
    "ambiguous", "resource_unavailable", "provider_missing", "input_invalid", "execution_failed",
    "parse_failed", "version_mismatch", "skipped_by_policy",
}

ENTITY_COLUMNS = [
    "canonical_candidate_id", "species", "assembly_id", "annotation_release", "namespace",
    "external_id", "entity_type", "relationship", "mapping_status", "mapping_method",
    "source_id", "resource_version", "confidence", "notes",
]
PROVIDER_RUN_COLUMNS = [
    "provider_name", "provider_version", "resource_name", "resource_version", "input_fingerprint",
    "output_fingerprint", "coverage_status", "run_status", "failure_code", "failure_message",
    "started_at", "finished_at", "records_in", "records_out",
]
COVERAGE_COLUMNS = [
    "canonical_candidate_id", "candidate_gene_id", "provider_name", "coverage_status", "record_count",
    "failure_code", "failure_message",
]
OBSERVATION_COLUMNS = [
    "observation_id", "canonical_candidate_id", "candidate_gene_id", "observation_type",
    "observation_category", "provider_name", "source_id", "source_family", "independence_family",
    "derived_from", "coverage_status", "review_status", "value_json",
]
LINEAGE_COLUMNS = [
    "observation_id", "derived_from", "relationship", "canonical_candidate_id", "source_family",
    "independence_family",
]


def _stable_id(prefix: str, *parts: object) -> str:
    payload = json.dumps(parts, ensure_ascii=False, sort_keys=True, separators=(",", ":"), default=str)
    return f"{prefix}_{hashlib.sha256(payload.encode('utf-8')).hexdigest()[:20]}"


def _fingerprint_identity(path: Path) -> dict[str, object]:
    if not path.exists():
        return {"path": str(path), "status": "missing"}
    fingerprint = file_fingerprint(path, False)
    identity = {
        key: fingerprint[key]
        for key in ("path", "size_bytes", "sha256", "method", "strength", "fai_path", "fai_sha256")
        if key in fingerprint
    }
    if fingerprint.get("strength") != "strong":
        identity["mtime_ns"] = fingerprint.get("mtime_ns")
    return identity


class WorkflowState:
    """Small deterministic step cache keyed by inputs, resources, and parameters."""

    def __init__(self, path: Path):
        self.path = path
        if path.exists():
            payload = json.loads(path.read_text(encoding="utf-8"))
            self.payload = payload if isinstance(payload, dict) else {}
        else:
            self.payload = {}
        self.payload.setdefault("schema_version", FOUNDATION_SCHEMA_VERSION)
        self.payload.setdefault("steps", {})

    @staticmethod
    def signature(inputs: Iterable[Path], parameters: dict[str, object]) -> str:
        payload = {
            "schema_version": FOUNDATION_SCHEMA_VERSION,
            "inputs": [_fingerprint_identity(path) for path in inputs],
            "parameters": parameters,
        }
        encoded = json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":"), default=str)
        return hashlib.sha256(encoded.encode("utf-8")).hexdigest()

    def can_reuse(self, step_name: str, signature: str, outputs: Iterable[Path]) -> bool:
        step = self.payload["steps"].get(step_name, {})
        output_list = list(outputs)
        if step.get("run_status") != "complete" or step.get("signature") != signature:
            return False
        if any(not path.exists() for path in output_list):
            return False
        recorded = step.get("output_fingerprints", [])
        return recorded == [_fingerprint_identity(path) for path in output_list]

    def record(
        self,
        step_name: str,
        signature: str,
        run_status: str,
        outputs: Iterable[Path] = (),
        dependencies: Iterable[str] = (),
        failure_code: str = "",
        failure_message: str = "",
        reused: bool = False,
    ) -> None:
        output_list = list(outputs)
        self.payload["steps"][step_name] = {
            "signature": signature,
            "dependencies": list(dependencies),
            "run_status": run_status,
            "failure_code": failure_code,
            "failure_message": failure_message,
            "reused": bool(reused),
            "output_fingerprints": [_fingerprint_identity(path) for path in output_list if path.exists()],
            "updated_at": datetime.now(UTC).isoformat(),
        }
        self.save()

    def save(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.path.with_suffix(self.path.suffix + ".tmp")
        temporary.write_text(
            json.dumps(self.payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        temporary.replace(self.path)


def _canonical_id(assembly_id: str, gene_id: str) -> str:
    return f"{assembly_id}:{gene_id}"


def _coverage_from_mapping(mapping_status: str, failure_reason: str, provider_failed: bool) -> tuple[str, str]:
    if provider_failed or mapping_status == "provider_failure":
        return "execution_failed", "provider_execution_failed"
    if mapping_status == "mapped":
        return "complete", ""
    if failure_reason == "missing_candidate_protein":
        return "input_invalid", "missing_candidate_protein"
    if failure_reason == "no_quality_passing_hit":
        return "no_hit", "no_quality_passing_hit"
    if mapping_status == "unmapped":
        return "unmapped", failure_reason or "unmapped"
    return "partial", failure_reason or "partial_result"


def write_foundation_outputs(
    *,
    out_dir: Path,
    candidates: pd.DataFrame,
    protein_ids: set[str],
    panel: dict[str, object],
    panel_path: Path,
    orthology: pd.DataFrame,
    functions: pd.DataFrame,
    trait_matches: pd.DataFrame,
    variants: pd.DataFrame,
    evidence: pd.DataFrame,
    ranked: pd.DataFrame,
    provider_runs: list[dict[str, object]],
    provider_version: str,
    inputs: dict[str, Path | None],
    mode: str,
    id_map: pd.DataFrame | None = None,
) -> dict[str, str]:
    """Write the canonical Issue-120 artifacts for an auto-triage run."""
    native = panel.get("native", {}) if isinstance(panel.get("native"), dict) else {}
    species = str(native.get("species", "unknown"))
    assembly_id = str(native.get("assembly_id") or native.get("assembly") or "assembly_unknown")
    annotation_release = str(native.get("annotation_release") or native.get("resource_version") or "release_unknown")
    references = {str(item["label"]): item for item in panel.get("references", []) if isinstance(item, dict)}
    canonical_by_gene = {str(gene): _canonical_id(assembly_id, str(gene)) for gene in candidates["gene_id"]}

    entity_rows: list[dict[str, object]] = []
    for row in candidates.to_dict(orient="records"):
        gene_id = str(row["gene_id"])
        canonical = canonical_by_gene[gene_id]
        source_id = f"native:{assembly_id}:{annotation_release}"
        entity_rows.append({
            "canonical_candidate_id": canonical, "species": species, "assembly_id": assembly_id,
            "annotation_release": annotation_release, "namespace": "native_gene", "external_id": gene_id,
            "entity_type": "gene", "relationship": "same_entity", "mapping_status": "exact",
            "mapping_method": "native_gff", "source_id": source_id, "resource_version": annotation_release,
            "confidence": "high", "notes": "canonical QTL-coordinate identity",
        })
        transcript = str(row.get("representative_transcript", "")).strip()
        if transcript:
            entity_rows.append({
                "canonical_candidate_id": canonical, "species": species, "assembly_id": assembly_id,
                "annotation_release": annotation_release, "namespace": "native_transcript", "external_id": transcript,
                "entity_type": "transcript", "relationship": "transcript_of", "mapping_status": "exact",
                "mapping_method": "native_gff", "source_id": source_id, "resource_version": annotation_release,
                "confidence": "high", "notes": "representative transcript; other isoforms remain in extraction tables",
            })
        if gene_id in protein_ids:
            entity_rows.append({
                "canonical_candidate_id": canonical, "species": species, "assembly_id": assembly_id,
                "annotation_release": annotation_release, "namespace": "native_protein", "external_id": gene_id,
                "entity_type": "protein", "relationship": "protein_of", "mapping_status": "exact",
                "mapping_method": "candidate_protein_fasta", "source_id": source_id,
                "resource_version": annotation_release, "confidence": "high", "notes": "",
            })

    if id_map is not None and not id_map.empty:
        required = {"gene_id", "transcript_id", "protein_id"}
        if not required.issubset(id_map.columns):
            raise ValueError("id_map requires gene_id, transcript_id, and protein_id")
        for row in id_map.fillna("").to_dict(orient="records"):
            gene_id = str(row["gene_id"])
            if gene_id not in canonical_by_gene:
                continue
            canonical = canonical_by_gene[gene_id]
            common = {
                "canonical_candidate_id": canonical, "species": species, "assembly_id": assembly_id,
                "annotation_release": annotation_release, "mapping_status": "exact",
                "mapping_method": "native_id_map", "source_id": f"native:{assembly_id}:{annotation_release}",
                "resource_version": annotation_release, "confidence": "high", "notes": "",
            }
            transcript_id = str(row.get("transcript_id", "")).strip()
            if transcript_id:
                entity_rows.append({
                    **common, "namespace": "native_transcript", "external_id": transcript_id,
                    "entity_type": "transcript", "relationship": "transcript_of",
                })
            protein_id = str(row.get("protein_id", "")).strip()
            if protein_id:
                entity_rows.append({
                    **common, "namespace": "native_protein", "external_id": protein_id,
                    "entity_type": "protein", "relationship": "protein_of",
                })
            sequence_record_id = str(row.get("sequence_record_id", "")).strip()
            if sequence_record_id:
                entity_rows.append({
                    **common, "namespace": "candidate_protein_fasta", "external_id": sequence_record_id,
                    "entity_type": "protein_sequence_record", "relationship": "sequence_record_for",
                })

    mapped = orthology[orthology["mapping_status"].eq("mapped")].copy()
    external_counts = mapped.groupby(["reference_label", "reference_gene_id"])["candidate_gene_id"].nunique().to_dict()
    for row in mapped.to_dict(orient="records"):
        gene_id = str(row["candidate_gene_id"])
        label = str(row["reference_label"])
        reference = references.get(label, {})
        external_id = str(row["reference_gene_id"])
        relationship = str(row.get("relationship", "ortholog"))
        count = external_counts.get((label, external_id), 1)
        mapping_status = "ambiguous" if str(row.get("paralog_ambiguity")) == "top2_near_tie" else "one_to_many" if count > 1 or relationship == "one_to_many" else "exact"
        entity_rows.append({
            "canonical_candidate_id": canonical_by_gene[gene_id], "species": str(row.get("reference_species", "")),
            "assembly_id": str(reference.get("assembly_id", "")),
            "annotation_release": str(reference.get("annotation_release", "")), "namespace": label,
            "external_id": external_id, "entity_type": "database_record", "relationship": relationship,
            "mapping_status": mapping_status, "mapping_method": str(row.get("mapping_providers", "")),
            "source_id": f"resource:{label}", "resource_version": str(reference.get("resource_version", "unknown")),
            "confidence": str(row.get("orthology_confidence", "")), "notes": str(row.get("paralog_ambiguity", "")),
        })
    entities = pd.DataFrame(entity_rows, columns=ENTITY_COLUMNS).drop_duplicates()
    entity_path = out_dir / "entity_aliases.tsv"
    entities.to_csv(entity_path, sep="\t", index=False)

    run_by_label = {str(item.get("label", "")): item for item in provider_runs}
    provider_table_rows: list[dict[str, object]] = []
    coverage_rows: list[dict[str, object]] = []
    for label, reference in references.items():
        run = run_by_label.get(label, {})
        failed = str(run.get("status", "")) == "failed"
        label_rows = orthology[orthology["reference_label"].eq(label)]
        coverage_values: list[str] = []
        for gene_id in sorted(canonical_by_gene):
            gene_rows = label_rows[label_rows["candidate_gene_id"].eq(gene_id)]
            states = [
                _coverage_from_mapping(str(row.get("mapping_status", "")), str(row.get("failure_reason", "")), failed)
                for row in gene_rows.to_dict(orient="records")
            ]
            mapped_count = int(gene_rows["mapping_status"].eq("mapped").sum()) if not gene_rows.empty else 0
            if failed:
                coverage_status, failure_code = "execution_failed", "provider_execution_failed"
            elif mapped_count > 1:
                coverage_status, failure_code = "ambiguous", "one_to_many_mapping"
            elif mapped_count == 1:
                coverage_status, failure_code = "complete", ""
            elif states:
                coverage_status, failure_code = states[0]
            else:
                coverage_status, failure_code = "not_run", "candidate_missing_from_provider_output"
            if coverage_status not in COVERAGE_STATES:
                raise ValueError(f"invalid coverage state: {coverage_status}")
            coverage_values.append(coverage_status)
            coverage_rows.append({
                "canonical_candidate_id": canonical_by_gene[gene_id], "candidate_gene_id": gene_id,
                "provider_name": f"blastp:{label}", "coverage_status": coverage_status,
                "record_count": mapped_count, "failure_code": failure_code,
                "failure_message": ";".join(sorted(set(filter(None, gene_rows.get("failure_reason", []))))),
            })
            function_gene_rows = functions[
                functions["candidate_gene_id"].eq(gene_id) & functions["reference_label"].eq(label)
            ] if not functions.empty else functions
            annotation_status = str(run.get("annotation_status", "not_run"))
            if annotation_status == "failed":
                annotation_coverage, annotation_failure = "execution_failed", "annotation_parse_failed"
            elif mapped_count == 0:
                annotation_coverage, annotation_failure = "not_run", "mapping_required"
            elif function_gene_rows.empty:
                annotation_coverage, annotation_failure = "partial", "annotation_missing_for_mapping"
            elif function_gene_rows["retrieval_status"].eq("source_retrieved").any():
                annotation_coverage, annotation_failure = "complete", ""
            else:
                annotation_coverage, annotation_failure = "partial", "fasta_description_or_annotation_unavailable"
            coverage_rows.append({
                "canonical_candidate_id": canonical_by_gene[gene_id], "candidate_gene_id": gene_id,
                "provider_name": f"annotation:{label}", "coverage_status": annotation_coverage,
                "record_count": len(function_gene_rows), "failure_code": annotation_failure,
                "failure_message": str(run.get("annotation_error", "")) if annotation_status == "failed" else "",
            })
        run_status = "failed" if failed else "complete"
        if failed:
            aggregate_coverage = "execution_failed"
        elif coverage_values and all(value == "no_hit" for value in coverage_values):
            aggregate_coverage = "no_hit"
        elif any(value in {"complete", "ambiguous"} for value in coverage_values):
            aggregate_coverage = "complete" if all(value == "complete" for value in coverage_values) else "partial"
        else:
            aggregate_coverage = "partial"
        blast = run.get("blast", {}) if isinstance(run.get("blast"), dict) else {}
        output = Path(str(blast.get("output", ""))) if blast.get("output") else None
        reference_path = Path(str(reference.get("_resolved_proteome", reference.get("proteome", ""))))
        provider_table_rows.append({
            "provider_name": f"blastp:{label}", "provider_version": provider_version,
            "resource_name": label, "resource_version": str(reference.get("resource_version", "unknown")),
            "input_fingerprint": json.dumps(
                [_fingerprint_identity(inputs["candidate_proteins"]), _fingerprint_identity(reference_path)],
                sort_keys=True,
            ),
            "output_fingerprint": json.dumps(_fingerprint_identity(output), sort_keys=True) if output else "",
            "coverage_status": aggregate_coverage, "run_status": run_status,
            "failure_code": "provider_execution_failed" if failed else "",
            "failure_message": str(run.get("error", "")), "started_at": str(run.get("started_at", "")),
            "finished_at": str(run.get("finished_at", "")), "records_in": len(candidates),
            "records_out": int(label_rows["mapping_status"].eq("mapped").sum()),
        })
        annotation_status = str(run.get("annotation_status", "not_run"))
        annotation_value = str(reference.get("annotation", "")).strip()
        annotation_path = panel_path.parent / annotation_value if annotation_value and not Path(annotation_value).is_absolute() else Path(annotation_value) if annotation_value else None
        provider_table_rows.append({
            "provider_name": f"annotation:{label}", "provider_version": str(reference.get("resource_version", "unknown")),
            "resource_name": label, "resource_version": str(reference.get("resource_version", "unknown")),
            "input_fingerprint": json.dumps(_fingerprint_identity(annotation_path), sort_keys=True) if annotation_path else json.dumps({"status": "not_requested"}),
            "output_fingerprint": "", "coverage_status": "execution_failed" if annotation_status == "failed" else "complete" if annotation_status == "complete" else "not_run",
            "run_status": "failed" if annotation_status == "failed" else "complete" if annotation_status == "complete" else "skipped",
            "failure_code": "annotation_parse_failed" if annotation_status == "failed" else "",
            "failure_message": str(run.get("annotation_error", "")), "started_at": str(run.get("started_at", "")),
            "finished_at": str(run.get("finished_at", "")), "records_in": int(label_rows["mapping_status"].eq("mapped").sum()),
            "records_out": int(functions["reference_label"].eq(label).sum()) if not functions.empty else 0,
        })
    provider_table = pd.DataFrame(provider_table_rows, columns=PROVIDER_RUN_COLUMNS)
    provider_path = out_dir / "provider_runs.tsv"
    provider_table.to_csv(provider_path, sep="\t", index=False)
    coverage = pd.DataFrame(coverage_rows, columns=COVERAGE_COLUMNS)
    coverage_path = out_dir / "candidate_provider_coverage.tsv"
    coverage.to_csv(coverage_path, sep="\t", index=False)

    observation_rows: list[dict[str, object]] = []
    mapping_ids: dict[tuple[str, str, str], str] = {}
    function_ids: dict[tuple[str, str, str], str] = {}
    evidence_ids: dict[str, list[str]] = {gene: [] for gene in canonical_by_gene}

    def add_observation(
        gene_id: str, observation_type: str, category: str, provider: str, source_id: str,
        source_family: str, independence_family: str, derived_from: list[str], coverage_status: str,
        review_status: str, value: dict[str, object], discriminator: object,
    ) -> str:
        observation_id = _stable_id("obs", canonical_by_gene[gene_id], observation_type, provider, discriminator)
        observation_rows.append({
            "observation_id": observation_id, "canonical_candidate_id": canonical_by_gene[gene_id],
            "candidate_gene_id": gene_id, "observation_type": observation_type,
            "observation_category": category, "provider_name": provider, "source_id": source_id,
            "source_family": source_family, "independence_family": independence_family,
            "derived_from": ";".join(sorted(set(filter(None, derived_from)))),
            "coverage_status": coverage_status, "review_status": review_status,
            "value_json": json.dumps(value, ensure_ascii=False, sort_keys=True, default=str),
        })
        return observation_id

    for row in orthology.to_dict(orient="records"):
        gene_id = str(row["candidate_gene_id"])
        label = str(row["reference_label"])
        failed = str(row.get("mapping_status")) == "provider_failure"
        coverage_status, _ = _coverage_from_mapping(str(row.get("mapping_status", "")), str(row.get("failure_reason", "")), failed)
        observation_id = add_observation(
            gene_id, "identity_mapping", "mapping_robustness", f"blastp:{label}", f"resource:{label}",
            f"mapping:{label}", f"comparative_function:{label}", [], coverage_status, "unreviewed", row,
            (label, row.get("reference_gene_id", ""), row.get("relationship", ""), row.get("failure_reason", "")),
        )
        if row.get("reference_gene_id"):
            mapping_ids[(gene_id, label, str(row["reference_gene_id"]))] = observation_id

    for row in functions.to_dict(orient="records"):
        gene_id = str(row["candidate_gene_id"])
        label = str(row["reference_label"])
        reference_gene = str(row["reference_gene_id"])
        parent = mapping_ids.get((gene_id, label, reference_gene), "")
        observation_id = add_observation(
            gene_id, "biological_annotation", "biological_annotation", f"annotation:{label}",
            str(row.get("source_id", "")), f"annotation:{label}", f"comparative_function:{label}",
            [parent], "complete" if row.get("retrieval_status") == "source_retrieved" else "partial",
            str(row.get("review_status", "unreviewed")), row, (label, reference_gene, row.get("source_id", "")),
        )
        function_ids[(gene_id, label, reference_gene)] = observation_id

    for row in variants.to_dict(orient="records"):
        gene_id = str(row["gene_id"])
        add_observation(
            gene_id, "parental_variant", "locus_specific_evidence", "parent_vcf", "local_file:parent_vcf",
            "parental_variant", "locus_specific_variant", [], "complete", str(row.get("review_status", "unreviewed")),
            row, (row.get("chrom"), row.get("pos"), row.get("ref"), row.get("alt")),
        )

    trait_parent: dict[tuple[str, str], str] = {}
    for row in trait_matches.to_dict(orient="records"):
        gene_id = str(row["gene_id"])
        label = str(row["reference_label"])
        reference_gene = str(row["reference_gene_id"])
        parent = function_ids.get((gene_id, label, reference_gene), "")
        observation_id = add_observation(
            gene_id, "trait_function_match", "derived_result", "trait_matcher", str(row.get("source_id", "")),
            "trait_match", "comparative_function", [parent], "complete", str(row.get("review_status", "unreviewed")),
            row, (reference_gene, row.get("matched_terms", "")),
        )
        trait_parent[(gene_id, str(row.get("source_id", "")))] = observation_id

    for index, row in enumerate(evidence.to_dict(orient="records")):
        gene_id = str(row["gene_id"])
        parent = trait_parent.get((gene_id, str(row.get("source_id", ""))), "")
        evidence_id = add_observation(
            gene_id, "screening_evidence", "derived_result", "evidence_normalizer", str(row.get("source_id", "")),
            str(row.get("source_type", "evidence")), "comparative_function", [parent], "complete",
            "unreviewed", row, (index, row.get("evidence_code", ""), row.get("source_id", "")),
        )
        evidence_ids.setdefault(gene_id, []).append(evidence_id)

    for row in ranked.to_dict(orient="records"):
        gene_id = str(row["gene_id"])
        add_observation(
            gene_id, "automated_rank", "derived_result", "deterministic_ranker", "local_file:06_ranked_candidates.tsv",
            "derived_rank", "derived_result", evidence_ids.get(gene_id, []), "complete", "unreviewed",
            {key: row.get(key, "") for key in ("automated_rank", "automated_score", "dry_group", "score_breakdown")},
            row.get("automated_rank", ""),
        )

    observations = pd.DataFrame(observation_rows, columns=OBSERVATION_COLUMNS)
    observation_path = out_dir / "normalized_observations.tsv"
    observations.to_csv(observation_path, sep="\t", index=False)
    lineage_rows: list[dict[str, str]] = []
    for row in observations.to_dict(orient="records"):
        for parent in filter(None, str(row["derived_from"]).split(";")):
            lineage_rows.append({
                "observation_id": str(row["observation_id"]), "derived_from": parent,
                "relationship": "derived_from", "canonical_candidate_id": str(row["canonical_candidate_id"]),
                "source_family": str(row["source_family"]), "independence_family": str(row["independence_family"]),
            })
    lineage = pd.DataFrame(lineage_rows, columns=LINEAGE_COLUMNS)
    lineage_path = out_dir / "evidence_lineage.tsv"
    lineage.to_csv(lineage_path, sep="\t", index=False)

    resource_manifest = {
        "schema_version": FOUNDATION_SCHEMA_VERSION,
        "mode": mode,
        "native": {"species": species, "assembly_id": assembly_id, "annotation_release": annotation_release},
        "resource_panel": _fingerprint_identity(panel_path),
        "inputs": {name: _fingerprint_identity(path) if path else {"status": "not_requested"} for name, path in inputs.items()},
        "references": [
            {
                "label": label, "role": str(reference.get("role", "")),
                "resource_version": str(reference.get("resource_version", "unknown")),
                "proteome": _fingerprint_identity(Path(str(reference.get("_resolved_proteome", reference.get("proteome", ""))))),
            }
            for label, reference in references.items()
        ],
    }
    resource_path = out_dir / "resource_manifest.json"
    resource_path.write_text(json.dumps(resource_manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "entity_aliases": str(entity_path), "provider_runs": str(provider_path),
        "candidate_provider_coverage": str(coverage_path), "normalized_observations": str(observation_path),
        "evidence_lineage": str(lineage_path), "resource_manifest": str(resource_path),
    }
