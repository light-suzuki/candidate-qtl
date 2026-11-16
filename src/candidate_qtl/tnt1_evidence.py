"""Normalize curated Medicago Tnt1 evidence without conflating insertion with causality.

Issue #524 requires Tnt1/FST observations to remain R108-background source
records until explicit identity evidence connects any reported A17-style gene ID
to the relevant R108 copy.  This module therefore produces a lossless evidence
ledger plus a deliberately narrow ``normalized_functional_evidence.tsv``
adapter: only genotype-supported or stronger causal evidence on a gene ID whose
accession matches the declared genetic background may enter that adapter.

No numeric score or candidate filtering is performed here.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import pandas as pd

from .core import file_fingerprint, runtime_path
from .reference_resources import FUNCTION_COLUMNS, load_reference_resource_registry, parse_lis_full_gene_id

SCHEMA_VERSION = "1.0"
EVIDENCE_STATES = {"supported", "not_supported", "not_evaluated"}

TNT1_COLUMNS = [
    "source_record_id",
    "reported_gene_id",
    "reported_gene_id_namespace",
    "reported_accession",
    "genetic_background",
    "identity_relation_state",
    "resource_match_state",
    "mutant_line",
    "fst_id",
    "insertion_site",
    "insertion_observed",
    "published_mutant_phenotype",
    "segregation_supported",
    "independent_alleles_supported",
    "complementation_or_rescue_supported",
    "evidence_class",
    "causal_gene_support_state",
    "route_b_function_transfer_ready",
    "route_b_block_reason",
    "phenotype_or_function",
    "gene_symbol",
    "ontology_ids",
    "source_ids",
    "source_url",
    "notes",
]

AUDIT_COLUMNS = ["record_index", "reported_gene_id", "status", "reason"]

_REQUIRED_COLUMNS = [
    "reported_gene_id",
    "genetic_background",
    "mutant_line",
    "insertion_observed",
    "published_mutant_phenotype",
    "segregation_supported",
    "independent_alleles_supported",
    "complementation_or_rescue_supported",
]

_TRUE_TOKENS = {"1", "true", "yes", "y", "supported", "present", "observed"}
_FALSE_TOKENS = {"0", "false", "no", "n", "not_supported", "absent", "refuted"}
_UNKNOWN_TOKENS = {"", "na", "n/a", "none", "unknown", "not_evaluated", "not evaluated", "unreported"}
_PUBLICATION_ID_PATTERN = re.compile(
    r"(?i)(?:^|[;,|\s])(?:pmid\s*:\s*\d+|pmcid\s*:\s*PMC\d+|(?:doi\s*:\s*)?10\.\d{4,9}/[^\s;,|]+)(?=$|[;,|\s])"
)


def _state(value: object, field: str) -> str:
    text = str(value or "").strip().casefold()
    if text in _TRUE_TOKENS:
        return "supported"
    if text in _FALSE_TOKENS:
        return "not_supported"
    if text in _UNKNOWN_TOKENS:
        return "not_evaluated"
    raise ValueError(f"{field} has invalid evidence state {value!r}; expected supported/not_supported/not_evaluated")


def _record_id(row: dict[str, object]) -> str:
    encoded = json.dumps(row, ensure_ascii=False, sort_keys=True, default=str, separators=(",", ":"))
    return "tnt1-" + hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:16]


def _has_publication_identifier(value: str) -> bool:
    """Require a publication ID, not an arbitrary ontology/database token."""
    return bool(_PUBLICATION_ID_PATTERN.search(value))


def _validate_state_consistency(states: dict[str, str]) -> None:
    """Reject explicit/incomplete evidence contradictions before precedence classification.

    ``not_evaluated`` remains unknown: it is never rewritten to negative.  But a
    higher-order support flag that logically asserts a published mutant
    phenotype cannot substitute for an unevaluated or explicitly unsupported
    phenotype field in the same curated record.
    """
    phenotype = states["published_mutant_phenotype"]
    insertion = states["insertion_observed"]
    contradictions: list[str] = []

    if states["segregation_supported"] == "supported" and phenotype != "supported":
        contradictions.append("segregation_supported_requires_published_mutant_phenotype")
    if states["independent_alleles_supported"] == "supported" and phenotype != "supported":
        contradictions.append("independent_alleles_supported_requires_published_mutant_phenotype")
    if states["complementation_or_rescue_supported"] == "supported" and phenotype != "supported":
        contradictions.append("complementation_or_rescue_supported_requires_published_mutant_phenotype")

    if insertion == "not_supported" and states["segregation_supported"] == "supported":
        contradictions.append("segregation_supported_conflicts_with_insertion_not_supported")
    if insertion == "not_supported" and (
        states["independent_alleles_supported"] == "supported"
        or states["complementation_or_rescue_supported"] == "supported"
    ):
        contradictions.append("strong_causal_flag_conflicts_with_insertion_not_supported")

    if contradictions:
        raise ValueError("contradictory Tnt1 evidence states: " + ";".join(contradictions))


def _evidence_class(states: dict[str, str]) -> tuple[str, str]:
    if states["complementation_or_rescue_supported"] == "supported" or states["independent_alleles_supported"] == "supported":
        return "strong_causal", "causal_supported"
    if states["published_mutant_phenotype"] == "supported" and states["segregation_supported"] == "supported":
        return "genotype_supported_mutant", "causal_supported"
    if states["published_mutant_phenotype"] == "supported":
        return "published_phenotype_only", "phenotype_associated"
    if states["insertion_observed"] == "supported":
        return "insertion_observed_only", "insertion_only"
    return "not_evaluated", "not_evaluated"


def normalize_tnt1_evidence(
    evidence_arg: str | Path,
    out_dir_arg: str | Path,
    *,
    registry_arg: str | Path | None = None,
    case_id: str = "",
) -> dict[str, object]:
    """Normalize a curated Tnt1 TSV into provenance-first evidence artifacts.

    The input is intentionally curated rather than scraped.  Every row must
    declare its genetic background and the gene ID exactly as reported/mapped by
    the source.  A background/ID-accession mismatch is preserved as
    ``requires_identity_bridge`` rather than silently rewritten.
    """
    evidence_path = runtime_path(evidence_arg)
    if not evidence_path.exists():
        raise ValueError(f"Tnt1 evidence does not exist: {evidence_path}")
    frame = pd.read_csv(evidence_path, sep="\t", dtype=str).fillna("")
    missing = [column for column in _REQUIRED_COLUMNS if column not in frame.columns]
    if missing:
        raise ValueError(f"Tnt1 evidence missing required columns: {', '.join(missing)}")

    registered_namespaces: set[str] = set()
    if registry_arg is not None:
        registry = load_reference_resource_registry(registry_arg)
        registered_namespaces = {row["native_gene_id_namespace"] for row in registry}

    rows: list[dict[str, str]] = []
    audit_rows: list[dict[str, object]] = []
    function_rows: list[dict[str, str]] = []

    for index, raw in enumerate(frame.to_dict("records"), 1):
        reported_gene_id = str(raw.get("reported_gene_id", "")).strip()
        parsed = parse_lis_full_gene_id(reported_gene_id)
        if parsed is None:
            audit_rows.append({"record_index": index, "reported_gene_id": reported_gene_id, "status": "invalid", "reason": "missing_or_invalid_full_lis_gene_id"})
            continue
        if parsed["species_code"].casefold() != "medtr":
            audit_rows.append({"record_index": index, "reported_gene_id": reported_gene_id, "status": "invalid", "reason": "reported_gene_id_is_not_medicago_truncatula"})
            continue
        background = str(raw.get("genetic_background", "")).strip()
        mutant_line = str(raw.get("mutant_line", "")).strip()
        if not background or not mutant_line:
            audit_rows.append({"record_index": index, "reported_gene_id": reported_gene_id, "status": "invalid", "reason": "genetic_background_and_mutant_line_are_required"})
            continue
        source_ids = str(raw.get("source_ids", "")).strip()
        source_url = str(raw.get("source_url", "")).strip()
        if not source_ids and not source_url:
            audit_rows.append({"record_index": index, "reported_gene_id": reported_gene_id, "status": "invalid", "reason": "source_ids_or_source_url_required"})
            continue

        try:
            states = {
                field: _state(raw.get(field, ""), field)
                for field in (
                    "insertion_observed",
                    "published_mutant_phenotype",
                    "segregation_supported",
                    "independent_alleles_supported",
                    "complementation_or_rescue_supported",
                )
            }
            _validate_state_consistency(states)
        except ValueError as exc:
            audit_rows.append({"record_index": index, "reported_gene_id": reported_gene_id, "status": "invalid", "reason": str(exc)})
            continue

        evidence_class, causal_state = _evidence_class(states)
        identity_state = (
            "native_background_match"
            if background.casefold() == parsed["accession"].casefold()
            else "requires_identity_bridge"
        )
        resource_state = "registered_namespace" if parsed["namespace"] in registered_namespaces else "unregistered_namespace"

        publication_id_present = _has_publication_identifier(source_ids)
        route_ready = causal_state == "causal_supported" and identity_state == "native_background_match" and publication_id_present
        if causal_state != "causal_supported":
            block_reason = "causal_gene_support_not_established"
        elif identity_state != "native_background_match":
            block_reason = "reported_gene_id_requires_explicit_background_identity_bridge"
        elif not publication_id_present:
            block_reason = "causal_evidence_missing_publication_identifier"
        else:
            block_reason = ""

        record = {
            "source_record_id": _record_id(raw),
            "reported_gene_id": parsed["full_id"],
            "reported_gene_id_namespace": parsed["namespace"],
            "reported_accession": parsed["accession"],
            "genetic_background": background,
            "identity_relation_state": identity_state,
            "resource_match_state": resource_state,
            "mutant_line": mutant_line,
            "fst_id": str(raw.get("fst_id", "")).strip(),
            "insertion_site": str(raw.get("insertion_site", "")).strip(),
            **states,
            "evidence_class": evidence_class,
            "causal_gene_support_state": causal_state,
            "route_b_function_transfer_ready": "true" if route_ready else "false",
            "route_b_block_reason": block_reason,
            "phenotype_or_function": str(raw.get("phenotype_or_function", "")).strip(),
            "gene_symbol": str(raw.get("gene_symbol", "")).strip(),
            "ontology_ids": str(raw.get("ontology_ids", "")).strip(),
            "source_ids": source_ids,
            "source_url": source_url,
            "notes": str(raw.get("notes", "")).strip(),
        }
        rows.append(record)
        audit_rows.append({
            "record_index": index,
            "reported_gene_id": parsed["full_id"],
            "status": "normalized",
            "reason": block_reason or "route_b_function_transfer_ready",
        })

        if route_ready:
            function_rows.append({
                "source_record_id": record["source_record_id"],
                "native_gene_id": parsed["full_id"],
                "native_species": str(raw.get("native_species", "Medicago truncatula")).strip() or "Medicago truncatula",
                "native_species_code": parsed["species_code"],
                "native_accession": parsed["accession"],
                "native_assembly": parsed["assembly_id"],
                "native_annotation": parsed["annotation_id"],
                "native_gene_id_namespace": parsed["namespace"],
                "resource_match_state": resource_state,
                "gene_symbol": record["gene_symbol"],
                "gene_model_pub_name": "",
                "classical_locus": "",
                "phenotype_or_function": record["phenotype_or_function"],
                "ontology_ids": record["ontology_ids"],
                "source_ids": source_ids,
                "source_confidence": evidence_class,
                "evidence_type": "tnt1_curated_mutant",
                "mutant_or_perturbation_background": background,
            })

    ledger = pd.DataFrame(rows, columns=TNT1_COLUMNS)
    functional = pd.DataFrame(function_rows, columns=FUNCTION_COLUMNS)
    audit = pd.DataFrame(audit_rows, columns=AUDIT_COLUMNS)

    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    ledger_path = out_dir / "tnt1_evidence.tsv"
    functional_path = out_dir / "normalized_functional_evidence.tsv"
    audit_path = out_dir / "tnt1_evidence_audit.tsv"
    manifest_path = out_dir / "tnt1_evidence_manifest.json"
    ledger.to_csv(ledger_path, sep="\t", index=False, lineterminator="\n")
    functional.to_csv(functional_path, sep="\t", index=False, lineterminator="\n")
    audit.to_csv(audit_path, sep="\t", index=False, lineterminator="\n")

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "producer_class": "tnt1_curated_evidence_normalizer",
        "case_id": case_id,
        "rank_contribution": "none",
        "cross_accession_identity_inferred": False,
        "input": file_fingerprint(evidence_path, full_hash=True),
        "registry": file_fingerprint(runtime_path(registry_arg), full_hash=True) if registry_arg is not None else None,
        "outputs": {
            "ledger": str(ledger_path),
            "normalized_functional_evidence": str(functional_path),
            "audit": str(audit_path),
        },
        "counts": {
            "input_records": int(len(frame)),
            "normalized_records": int(len(ledger)),
            "invalid_records": int((audit["status"] == "invalid").sum()) if not audit.empty else 0,
            "insertion_observed_only": int((ledger["evidence_class"] == "insertion_observed_only").sum()) if not ledger.empty else 0,
            "published_phenotype_only": int((ledger["evidence_class"] == "published_phenotype_only").sum()) if not ledger.empty else 0,
            "genotype_supported_mutant": int((ledger["evidence_class"] == "genotype_supported_mutant").sum()) if not ledger.empty else 0,
            "strong_causal": int((ledger["evidence_class"] == "strong_causal").sum()) if not ledger.empty else 0,
            "requires_identity_bridge": int((ledger["identity_relation_state"] == "requires_identity_bridge").sum()) if not ledger.empty else 0,
            "route_b_function_transfer_ready": int((ledger["route_b_function_transfer_ready"] == "true").sum()) if not ledger.empty else 0,
        },
        "policy": [
            "Tnt1 insertion/FST observation is distinct from causal mutant evidence",
            "genetic background is explicit on every normalized record",
            "a reported gene ID from another accession is preserved and requires an explicit identity bridge",
            "only genotype-supported or stronger causal evidence with native-background identity and a publication identifier enters normalized functional evidence",
            "explicitly contradictory curated evidence states are audited invalid before evidence-class precedence is applied; not_evaluated remains unknown rather than negative",
            "missing evidence remains not_evaluated rather than negative evidence",
            "no symbol-only identity mapping, numeric scoring, candidate filtering, or truth-tuned thresholding is performed",
        ],
    }
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "ledger": str(ledger_path),
        "normalized_functional_evidence": str(functional_path),
        "audit": str(audit_path),
        "manifest": str(manifest_path),
        "counts": manifest["counts"],
    }
