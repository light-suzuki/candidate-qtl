from __future__ import annotations

from pathlib import Path

import pandas as pd

from candidate_qtl.claim_binding import BOUND_CLAIM_COLUMNS


def bound_claim(
    *,
    route: str,
    claim_id: str,
    citation: str,
    subject_gene: str,
    organism: str,
    phenotype: str = "target phenotype",
    **overrides: str,
) -> dict[str, str]:
    row = {column: "" for column in BOUND_CLAIM_COLUMNS}
    row.update({
        "bound_claim_id": f"{route}:{claim_id}",
        "route": route,
        "claim_id": claim_id,
        "citation": citation,
        "source_record_id": citation,
        "subject_gene_identity": subject_gene,
        "subject_gene_verification_state": "verified_exact",
        "experiment_perturbation": "loss-of-function allele",
        "experiment_subject_gene_identity": subject_gene,
        "experiment_verification_state": "verified_same_subject",
        "predicate_relation": "causal_perturbation_supports_target_trait",
        "relation_verification_state": "verified",
        "object_phenotype_process": phenotype,
        "object_match_state": "target_trait",
        "organism": organism,
        "organism_verification_state": "verified_match",
        "accession": subject_gene,
        "accession_verification_state": "verified_match",
        "quoted_support": "The loss-of-function allele produced the target phenotype.",
        "locator": "abstract:sentence-3",
        "quote_binding_state": "verified_same_claim",
        "source_scope": "abstract",
        "binding_state": "verified_complete",
        "curation_method": "curated_human",
        "reviewed_by": "fixture-reviewer",
        "binding_reason": "all relation axes verified",
        "rank_contribution": "none",
    })
    row.update(overrides)
    return row


def write_route_a_bound_claims(
    root: Path,
    rows: list[dict[str, str]],
    *,
    name: str = "bound_claims.tsv",
) -> Path:
    claims = []
    for index, row in enumerate(rows, start=1):
        accession = str(row.get("proposed_stable_id", "")).strip()
        source = str(row.get("source_ids", "")).replace(",", ";").split(";")[0].strip()
        if not accession or not source:
            continue
        claims.append(bound_claim(
            route="a",
            claim_id=str(row.get("hypothesis_id", "")).strip() or f"mech_{index:03d}",
            citation=source,
            subject_gene=accession,
            organism=str(row.get("reference_species", "")).strip(),
            phenotype=str(row.get("pathway_or_process", "")).strip() or "target phenotype",
        ))
    path = root / name
    pd.DataFrame(claims, columns=BOUND_CLAIM_COLUMNS).to_csv(
        path, sep="\t", index=False, lineterminator="\n"
    )
    return path
