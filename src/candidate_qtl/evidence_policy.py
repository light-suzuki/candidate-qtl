from __future__ import annotations

import json

import pandas as pd

from .core import runtime_path
from .evidence import evidence_build as _base_evidence_build
from .evidence_io import EVIDENCE_COLUMNS, evidence_row as _row


def evidence_build(*args, **kwargs) -> dict[str, object]:
    """Apply canonical evidence-policy semantics after the base import.

    Two policy normalizations live here because this wrapper is the public
    ``evidence-build`` boundary used by the CLI and candidate pipeline:

    * invalid CDS phase is rescue/context for annotation quality, never silent
      exclusion;
    * ``paralog_ambiguous`` is epistemic copy/orthology uncertainty and is kept
      as unscored ``context`` rather than biological counterevidence (#421).

    The raw producer may still describe *why* the ambiguity was observed; this
    layer standardizes only its downstream evidence direction.
    """
    result = _base_evidence_build(*args, **kwargs)
    candidates_arg = args[0] if args else kwargs["candidates_arg"]
    candidates_path = runtime_path(candidates_arg)
    candidates = pd.read_csv(candidates_path, sep="\t", dtype=str).fillna("")

    evidence_path = runtime_path(str(result["evidence"]))
    evidence = pd.read_csv(evidence_path, sep="\t", dtype=str).fillna("")

    ambiguity_mask = evidence["evidence_code"].eq("paralog_ambiguous")
    ambiguity_rewritten = int((ambiguity_mask & ~evidence["direction"].eq("context")).sum())
    if ambiguity_mask.any():
        evidence.loc[ambiguity_mask, "direction"] = "context"

    invalid = pd.DataFrame()
    if "annotation_quality" in candidates.columns:
        invalid = candidates[candidates["annotation_quality"].eq("invalid_cds_phase")]
    additions = [
        _row(
            str(row["gene_id"]),
            "annotation_incomplete",
            "rescue",
            "high",
            "invalid_cds_phase",
            "gff3_phase_qc",
            candidates_path,
            "CDS phase is missing or inconsistent; retain the locus but do not trust its translated reference protein.",
        )
        for _, row in invalid.iterrows()
    ]
    if additions:
        evidence = pd.concat(
            [evidence, pd.DataFrame(additions, columns=EVIDENCE_COLUMNS)],
            ignore_index=True,
        )

    evidence = evidence.drop_duplicates().sort_values(
        ["gene_id", "evidence_code", "source_id"], kind="stable"
    )
    evidence.to_csv(evidence_path, sep="\t", index=False)

    # Candidate summary is cheap to rebuild and keeps the wrapper deterministic
    # after either policy normalization.
    summary_path = runtime_path(str(result["candidate_summary"]))
    summary = candidates.copy()
    if not evidence.empty:
        codes = evidence.groupby("gene_id")["evidence_code"].apply(
            lambda values: ";".join(sorted(set(values)))
        ).rename("evidence_codes")
        counts = evidence.groupby("gene_id").size().rename("evidence_count")
        summary = summary.merge(codes, on="gene_id", how="left").merge(
            counts,
            on="gene_id",
            how="left",
        )
    if "evidence_codes" not in summary.columns:
        summary["evidence_codes"] = ""
    else:
        summary["evidence_codes"] = summary["evidence_codes"].fillna("")
    if "evidence_count" not in summary.columns:
        summary["evidence_count"] = 0
    else:
        summary["evidence_count"] = summary["evidence_count"].fillna(0).astype(int)
    summary.to_csv(summary_path, sep="\t", index=False)

    manifest_path = runtime_path(str(result["manifest"]))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["evidence_count"] = len(evidence)
    manifest["invalid_cds_phase_rescue_count"] = len(invalid)
    manifest["paralog_ambiguous_context_count"] = int(ambiguity_mask.sum())
    manifest["paralog_ambiguous_direction_rewritten_count"] = ambiguity_rewritten
    manifest["paralog_ambiguity_policy"] = (
        "paralog_ambiguous is epistemic context and never native numerical counterevidence"
    )
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest | {"manifest": str(manifest_path)}
