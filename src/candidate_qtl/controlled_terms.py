from __future__ import annotations

import hashlib
import json
import re

import pandas as pd

from .core import runtime_path


def promote_controlled_terms(
    evidence_arg: str,
    out_dir_arg: str,
    terms: list[str],
    demote_unapproved_existing: bool = False,
) -> dict[str, object]:
    """Promote audited annotation/homology rows using an approved vocabulary.

    This deterministic screening helper never promotes broad interval context,
    never changes genetic evidence codes, and never creates source provenance.
    """
    cleaned = list(dict.fromkeys(term.strip() for term in terms if term.strip()))
    if not cleaned:
        raise ValueError("at least one controlled term is required")
    too_broad = [term for term in cleaned if len(re.sub(r"\W", "", term)) < 3]
    if too_broad:
        raise ValueError("controlled terms must contain at least three alphanumeric characters")
    evidence_path = runtime_path(evidence_arg)
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(evidence_path, sep="\t", dtype=str).fillna("")
    required = {"gene_id", "evidence_code", "direction", "confidence", "value", "source_type", "source_id", "description"}
    if not required.issubset(frame.columns):
        raise ValueError("evidence table is missing canonical columns")
    def matches(text: str) -> list[str]:
        return [
            term for term in cleaned
            if re.search(rf"(?<![A-Za-z0-9]){re.escape(term)}(?![A-Za-z0-9])", text, flags=re.IGNORECASE)
        ]

    allowed_codes = {"homology_supported", "homology_context", "domain_context", "go_annotation_context"}
    promoted_rows = []
    audit_rows = []
    if demote_unapproved_existing:
        for index, row in frame[frame["evidence_code"].eq("trait_function_match")].iterrows():
            matched = matches(f"{row['value']} {row['description']}")
            if matched:
                audit_rows.append({
                    "gene_id": row["gene_id"], "source_evidence_code": "trait_function_match",
                    "matched_terms": ";".join(matched), "value": row["value"], "source_id": row["source_id"],
                    "action": "retained_approved_existing_match",
                })
            else:
                frame.at[index, "evidence_code"] = "trait_context"
                frame.at[index, "direction"] = "context"
                frame.at[index, "description"] = "Demoted: existing trait match did not match the approved mechanism vocabulary."
                audit_rows.append({
                    "gene_id": row["gene_id"], "source_evidence_code": "trait_function_match",
                    "matched_terms": "", "value": row["value"], "source_id": row["source_id"],
                    "action": "demoted_unapproved_existing_match",
                })
    existing = set(frame.loc[frame["evidence_code"].eq("trait_function_match"), "gene_id"])
    for row in frame[frame["evidence_code"].isin(allowed_codes)].itertuples(index=False):
        text = f"{getattr(row, 'value')} {getattr(row, 'description')}"
        matched = matches(text)
        if not matched:
            continue
        audit_rows.append({
            "gene_id": row.gene_id, "source_evidence_code": row.evidence_code,
            "matched_terms": ";".join(matched), "value": row.value,
            "source_id": row.source_id,
            "action": "already_had_trait_function_match" if row.gene_id in existing else "promoted_low_confidence",
        })
        if row.gene_id in existing:
            continue
        record = {column: str(getattr(row, column)) for column in frame.columns}
        record.update({
            "evidence_code": "trait_function_match", "direction": "support",
            "confidence": "low", "source_type": f"controlled_term:{row.source_type}",
            "description": f"Approved deterministic term match: {';'.join(matched)}. Functional screening hypothesis only.",
        })
        promoted_rows.append(record)
        existing.add(str(row.gene_id))
    combined = pd.concat([frame, pd.DataFrame(promoted_rows, columns=frame.columns)], ignore_index=True)
    evidence_out = out_dir / "evidence.tsv"
    audit_out = out_dir / "controlled_term_audit.tsv"
    combined.to_csv(evidence_out, sep="\t", index=False)
    pd.DataFrame(audit_rows, columns=[
        "gene_id", "source_evidence_code", "matched_terms", "value", "source_id", "action"
    ]).to_csv(audit_out, sep="\t", index=False)
    manifest = {
        "schema_version": "1.0", "input": str(evidence_path),
        "input_sha256": hashlib.sha256(evidence_path.read_bytes()).hexdigest(),
        "approved_terms": cleaned, "promoted_gene_count": len(promoted_rows),
        "demote_unapproved_existing": demote_unapproved_existing,
        "confidence": "low", "llm_used": False,
        "provenance_policy": "original source_id is retained; this step never creates a citable source",
        "interpretation": "Controlled term matches are screening hypotheses, not orthology or causal evidence.",
    }
    manifest_out = out_dir / "controlled_term_manifest.json"
    manifest_out.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return {"evidence": str(evidence_out), "audit": str(audit_out), "manifest": str(manifest_out), "promoted_gene_count": len(promoted_rows)}
