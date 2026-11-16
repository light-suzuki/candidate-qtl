"""Leave-one-family-out corpus masking, generated from truth and frozen before prediction (issue #237).

A known-gene case can look easier than real discovery for a reason that has
nothing to do with the workflow's reasoning: the correct family may already sit
in the reference corpus by name, complete with a curated phenotype, an alias
list, and closely related examples one lexical hop away. Route B and Route C
both read that corpus, so a benchmark run against it partly measures "is the
answer already written down here" rather than "can this workflow discover
unfamiliar biology." Neither the honest Recall@K in ``docs/REFERENCE-CORPUS.ja.md``
nor any convergence number distinguishes the two unless the truth-shaped rows
are removed *before* prediction and the removal itself is audited.

This module is the removal step. ``corpus_mask`` takes a frozen reference
corpus and a mask specification -- one row of truth-derived facts per case,
built once during benchmark construction -- and produces three frozen
artifacts: the masked corpus (same schema as any other corpus, so it is a
drop-in input to ``retrieval-index-build`` / ``route-b-prepare``), an audit
table naming exactly which rows were removed and by which rule, and a manifest
that fingerprints the input corpus, the output corpus, and the mask rule
itself.

Three properties are load-bearing and each has a test:

1. **Masking happens strictly before prediction.** ``corpus-mask`` is a
   benchmark-construction command. Nothing downstream of it -- not
   ``retrieval-index-build``, not ``retrieval-run``, not ``route-b-prepare`` --
   accepts a mask spec or this manifest; they only ever see the masked corpus
   TSV, indistinguishable in shape from an unmasked one.
2. **The mask rule is truth-only and never reaches prediction.** The mask
   spec is read here and only here. The masked corpus carries no truth column,
   and this command refuses outright to build from a corpus that already has
   one, or to write a masked corpus that a fresh application of its own removal
   rule would still match (the self-check). A masking bug that leaves the
   answer behind is worse than no masking at all, because a downstream
   benchmark would score it as genuine novelty.
3. **Distant biology survives.** Matching is restricted to structured identity
   columns (``reference_gene_id``, ``gene_symbol``, ``family_id``,
   ``orthogroup_id``); free-text ``description`` / ``function_description`` /
   ``mutant_phenotype`` / ``pathway`` fields of *other* rows are never searched.
   ``leave_family_out`` removes the family's own rows and rows carrying a
   caller-supplied "clearly equivalent" alias -- not every row that happens to
   mention the same pathway or process.

The three tiers are nested by design: ``leave_gene_out`` removes the exact
reference gene; ``leave_family_out`` removes everything ``leave_gene_out``
would plus the family/subfamily; ``leave_close_orthogroup_out`` removes
everything ``leave_family_out`` would plus the whole orthogroup. Each stricter
tier is a superset of the one before it, so scoring the three tiers together
(``masking-benchmark``, in ``masking_benchmark.py``) reads as a monotone
novelty-stress gradient rather than three unrelated conditions.
"""

from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from . import __version__
from .core import file_fingerprint, runtime_path, safe_id, sha256
from .family_labels import families_match

SCHEMA_VERSION = "1.0"
MASK_FORMAT_VERSION = "1.0.0"

TIERS = ("leave_gene_out", "leave_family_out", "leave_close_orthogroup_out")

MASK_SPEC_REQUIRED_COLUMNS = ["case_id", "truth_gene_id"]
MASK_SPEC_OPTIONAL_COLUMNS = ("truth_family_id", "truth_orthogroup_id", "aliases")

AUDIT_COLUMNS = ["row_id", "reference_gene_id", "gene_symbol", "family_id", "orthogroup_id",
                 "tier", "matched_rule", "case_id"]

MASKED_CORPUS_NAME = "masked_reference_corpus.tsv"
AUDIT_NAME = "corpus_mask_audit.tsv"
MASK_MANIFEST_NAME = "corpus_mask_manifest.json"

_SPLIT = re.compile(r"[;,|]+")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _alias_set(raw: object) -> set[str]:
    return {item.strip().casefold() for item in _SPLIT.split(str(raw or "")) if item.strip()}


def _match_reason(
    row: dict[str, object],
    tier: str,
    truth_gene_key: str,
    truth_family: str,
    truth_orthogroup_key: str,
    alias_keys: set[str],
) -> str | None:
    """Which masking rule (if any) this corpus row trips, honouring tier nesting.

    Matching stays on structured identity columns only. Free text is never
    searched here, which is what lets more distant pathway/process rows
    survive ``leave_family_out`` intact.
    """
    gene_id = str(row.get("reference_gene_id", "")).strip()
    symbol = str(row.get("gene_symbol", "")).strip()
    family = str(row.get("family_id", "")).strip()
    orthogroup = str(row.get("orthogroup_id", "")).strip()

    if truth_gene_key and gene_id.casefold() == truth_gene_key:
        return "exact_reference_gene_id"
    if alias_keys and (gene_id.casefold() in alias_keys or symbol.casefold() in alias_keys):
        return "alias_match"

    if tier in ("leave_family_out", "leave_close_orthogroup_out") and truth_family:
        if family and families_match(family, truth_family):
            return "family_id_match"
        if orthogroup and families_match(orthogroup, truth_family):
            return "orthogroup_labelled_as_family_match"
        if symbol and families_match(symbol, truth_family):
            return "gene_symbol_family_match"

    if tier == "leave_close_orthogroup_out" and truth_orthogroup_key:
        if orthogroup and orthogroup.casefold() == truth_orthogroup_key:
            return "orthogroup_id_match"

    return None


def verify_masked_corpus(
    masked: pd.DataFrame,
    tier: str,
    truth_gene_id: str,
    truth_family_id: str = "",
    truth_orthogroup_id: str = "",
    aliases: set[str] | None = None,
) -> list[dict[str, str]]:
    """Reapply the masking predicate to a corpus that is supposed to already be masked.

    This is the self-check: ``corpus_mask`` calls it on the rows it is about to
    write, before anything reaches disk, and refuses the write if it returns
    anything. It is exported standalone so the same check can be run against a
    masked corpus produced any other way -- a leak here means the corpus still
    names the masked gene/family and must not be trusted as a novelty-stress
    input, whatever produced it.
    """
    if tier not in TIERS:
        raise ValueError(f"tier must be one of {', '.join(TIERS)}, got {tier!r}")
    truth_gene_key = str(truth_gene_id).strip().casefold()
    truth_orthogroup_key = str(truth_orthogroup_id).strip().casefold()
    alias_keys = {str(value).strip().casefold() for value in (aliases or set()) if str(value).strip()}
    leaks: list[dict[str, str]] = []
    for row in masked.to_dict("records"):
        reason = _match_reason(row, tier, truth_gene_key, str(truth_family_id).strip(), truth_orthogroup_key, alias_keys)
        if reason:
            leaks.append({"reference_gene_id": str(row.get("reference_gene_id", "")), "matched_rule": reason})
    return leaks


def _select_mask_spec_row(spec: pd.DataFrame, mask_spec_path: Path, case_id: str) -> dict[str, object]:
    rows = spec.to_dict("records")
    if case_id:
        matches = [row for row in rows if str(row["case_id"]).strip() == case_id]
        if not matches:
            raise ValueError(f"mask spec {mask_spec_path} has no row for --case-id {case_id!r}")
        if len(matches) > 1:
            raise ValueError(f"mask spec {mask_spec_path} has more than one row for --case-id {case_id!r}")
        return matches[0]
    if len(rows) == 1:
        return rows[0]
    raise ValueError(
        f"mask spec {mask_spec_path} has {len(rows)} rows; pass --case-id to select one (the spec may hold "
        "every case's truth-derived rule in one file, built once in the benchmark-construction phase)"
    )


def corpus_mask(
    corpus_arg: str,
    mask_spec_arg: str,
    tier: str,
    out_dir_arg: str,
    case_id: str = "",
    full_hash: bool = False,
) -> dict[str, object]:
    """Remove truth-shaped rows from a frozen corpus before prediction, and audit exactly what left."""
    if tier not in TIERS:
        raise ValueError(f"--tier must be one of {', '.join(TIERS)}, got {tier!r}")
    corpus_path = runtime_path(corpus_arg)
    mask_spec_path = runtime_path(mask_spec_arg)
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)

    corpus = pd.read_csv(corpus_path, sep="\t", dtype=str).fillna("")
    if "reference_gene_id" not in corpus.columns:
        raise ValueError(f"reference corpus {corpus_path} has no reference_gene_id column")
    truth_like = [column for column in corpus.columns if column.strip().casefold().startswith("truth")]
    if truth_like:
        raise ValueError(
            f"reference corpus {corpus_path} already carries truth-shaped column(s) {truth_like}; corpus-mask "
            "refuses to build a masked corpus from an input that names truth explicitly"
        )

    spec = pd.read_csv(mask_spec_path, sep="\t", dtype=str).fillna("")
    missing = [column for column in MASK_SPEC_REQUIRED_COLUMNS if column not in spec.columns]
    if missing:
        raise ValueError(f"mask spec {mask_spec_path} is missing columns: {', '.join(missing)}")
    record = _select_mask_spec_row(spec, mask_spec_path, case_id)
    case_id = str(record["case_id"]).strip()
    if not case_id:
        raise ValueError(f"mask spec {mask_spec_path} has a row with an empty case_id")

    truth_gene_id = str(record.get("truth_gene_id", "")).strip()
    if not truth_gene_id:
        raise ValueError(f"mask spec row for case {case_id!r} has no truth_gene_id")
    truth_family_id = str(record.get("truth_family_id", "")).strip()
    truth_orthogroup_id = str(record.get("truth_orthogroup_id", "")).strip()
    aliases = _alias_set(record.get("aliases", ""))

    if tier in ("leave_family_out", "leave_close_orthogroup_out") and not truth_family_id:
        raise ValueError(f"tier {tier!r} requires a non-empty truth_family_id in the mask spec for case {case_id!r}")
    if tier == "leave_close_orthogroup_out" and not truth_orthogroup_id:
        raise ValueError(f"tier {tier!r} requires a non-empty truth_orthogroup_id in the mask spec for case {case_id!r}")

    truth_gene_key = truth_gene_id.casefold()
    truth_orthogroup_key = truth_orthogroup_id.casefold()

    corpus_records = corpus.to_dict("records")
    gene_present = False
    own_symbol = ""
    for row in corpus_records:
        if str(row.get("reference_gene_id", "")).strip().casefold() == truth_gene_key:
            gene_present = True
            own_symbol = own_symbol or str(row.get("gene_symbol", "")).strip()
    # The truth gene's own symbol is a "clearly equivalent alias" of itself even
    # when the benchmark author forgot to list it explicitly in the spec.
    if own_symbol:
        aliases.add(own_symbol.casefold())

    kept_rows: list[dict[str, object]] = []
    audit_rows: list[dict[str, object]] = []
    rule_counts: dict[str, int] = {}
    for offset, row in enumerate(corpus_records, start=2):  # header is line 1
        reason = _match_reason(row, tier, truth_gene_key, truth_family_id, truth_orthogroup_key, aliases)
        if reason:
            rule_counts[reason] = rule_counts.get(reason, 0) + 1
            audit_rows.append({
                "row_id": offset, "reference_gene_id": str(row.get("reference_gene_id", "")),
                "gene_symbol": str(row.get("gene_symbol", "")), "family_id": str(row.get("family_id", "")),
                "orthogroup_id": str(row.get("orthogroup_id", "")), "tier": tier, "matched_rule": reason,
                "case_id": case_id,
            })
        else:
            kept_rows.append(row)

    masked = pd.DataFrame(kept_rows, columns=list(corpus.columns))
    if set(masked.columns) != set(corpus.columns):
        raise RuntimeError("internal error: masking must never change the corpus column schema")

    # Self-check: reapply the identical rule to what is about to be written.
    # Masking that silently fails is worse than no masking, because a
    # downstream benchmark would read the leftover row as genuine novelty.
    leaks = verify_masked_corpus(masked, tier, truth_gene_id, truth_family_id, truth_orthogroup_id, aliases)
    if leaks:
        examples = ", ".join(f"{leak['reference_gene_id']}({leak['matched_rule']})" for leak in leaks[:5])
        raise ValueError(
            f"corpus-mask self-check failed for case {case_id!r}, tier {tier!r}: {len(leaks)} row(s) still name "
            f"the masked gene/family after filtering (e.g. {examples}); refusing to write a leaking masked corpus"
        )

    case_dir = out / safe_id(case_id) / tier
    case_dir.mkdir(parents=True, exist_ok=True)
    masked_path = case_dir / MASKED_CORPUS_NAME
    masked.to_csv(masked_path, sep="\t", index=False, lineterminator="\n")
    audit_path = case_dir / AUDIT_NAME
    pd.DataFrame(audit_rows, columns=AUDIT_COLUMNS).to_csv(audit_path, sep="\t", index=False, lineterminator="\n")

    mask_rule = {
        "case_id": case_id, "tier": tier, "truth_gene_id": truth_gene_id,
        "truth_family_id": truth_family_id, "truth_orthogroup_id": truth_orthogroup_id,
        "aliases": sorted(aliases),
    }
    mask_fingerprint = hashlib.sha256(
        json.dumps(mask_rule, ensure_ascii=False, sort_keys=True).encode("utf-8")).hexdigest()

    manifest_path = case_dir / MASK_MANIFEST_NAME
    manifest = {
        "schema_version": SCHEMA_VERSION, "mask_format_version": MASK_FORMAT_VERSION,
        "candidate_qtl_version": __version__, "command": "corpus-mask",
        "built_at_utc": _now(),
        "case_id": case_id, "tier": tier,
        "mask_rule": mask_rule,
        "mask_fingerprint": mask_fingerprint,
        "input_row_count": int(len(corpus)),
        "removed_row_count": len(audit_rows),
        "masked_row_count": int(len(masked)),
        "removed_by_rule": rule_counts,
        "truth_gene_present_in_corpus": gene_present,
        "masking_had_no_effect": len(audit_rows) == 0,
        "self_check": {
            "performed": True, "leaks_found": 0,
            "rule": "the identical matching predicate is reapplied to the rows about to be written; any row that "
                    "still names the masked gene/family refuses the write instead of producing a silent leak",
        },
        "prediction_visibility": (
            "the mask spec and this manifest are read only by corpus-mask itself. retrieval-index-build, "
            "retrieval-run, route-b-prepare and every prediction command consume masked_reference_corpus.tsv like "
            "any other reference corpus and never see truth_gene_id, truth_family_id, aliases, or this manifest"
        ),
        "pathway_preservation_rule": (
            "matching is restricted to reference_gene_id / gene_symbol / family_id / orthogroup_id; the free-text "
            "description / function_description / mutant_phenotype / pathway fields of other rows are never "
            "searched, so rows that merely mention the same pathway or process survive leave_family_out"
        ),
        "tier_hierarchy": (
            "leave_gene_out rows are a subset of leave_family_out rows, which are a subset of "
            "leave_close_orthogroup_out rows: each stricter tier removes everything the one before it removes"
        ),
        "masking_timing": "masking runs strictly before prediction and is fingerprinted here",
        "no_truth_column_in_output": True,
        "corpus_sha256": sha256(corpus_path),
        "masked_corpus_sha256": sha256(masked_path),
        "audit_sha256": sha256(audit_path),
        "inputs": {"corpus": file_fingerprint(corpus_path, full_hash),
                   "mask_spec": file_fingerprint(mask_spec_path, full_hash)},
        "outputs": {"masked_corpus": str(masked_path), "audit": str(audit_path), "manifest": str(manifest_path)},
    }
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), "masked_corpus": str(masked_path), "audit": str(audit_path),
            "case_id": case_id, "tier": tier, "mask_fingerprint": mask_fingerprint,
            "input_row_count": int(len(corpus)), "removed_row_count": len(audit_rows),
            "masked_row_count": int(len(masked)), "removed_by_rule": rule_counts,
            "truth_gene_present_in_corpus": gene_present}


__all__ = ["corpus_mask", "verify_masked_corpus", "TIERS", "MASK_SPEC_REQUIRED_COLUMNS", "MASK_SPEC_OPTIONAL_COLUMNS"]
