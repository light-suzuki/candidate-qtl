"""Final machine competitor scan for semantic family/copy proposals (#505).

The semantic verifier establishes that a model-proposed candidate/reference pair
exists in the frozen machine facts.  That is necessary but not sufficient for a
*single* preferred copy: another unresolved family member may have the same
strong reference relation even if the model omitted it.

This final gate therefore treats the LLM/agent mapping as a reference hypothesis,
then enumerates *all* eligible family members with the same supplied reference
and one-to-one/RBH-like machine relation.  The machine, not the model, decides
whether the result is unique or a tie.
"""

from __future__ import annotations

import json
import re

import pandas as pd

from .core import file_fingerprint, runtime_path
from .family_semantic_resolution import ELIGIBLE_UNRESOLVED_STATES, SHORTLIST_COLUMNS

SCHEMA_VERSION = "1.0"
STRONG_ORTHOLOGY = frozenset({
    "1_1",
    "one_to_one",
    "one_one",
    "reciprocal_best_hit",
    "rbh",
})
FINAL_EXTRA_COLUMNS = [
    "semantic_match_basis",
    "semantic_verification_status",
    "semantic_reference_gene",
    "semantic_relation",
    "semantic_proposal_id",
]


def _split(value: object) -> set[str]:
    return {
        part.strip()
        for part in str(value or "").replace(",", ";").split(";")
        if part.strip()
    }


def _norm_relation(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value or "").strip().casefold()).strip("_")


def _strong(value: object) -> bool:
    return _norm_relation(value) in STRONG_ORTHOLOGY


def _atlas_references(row: dict[str, object]) -> set[str]:
    return _split(row.get("orthology_reference_gene", "") or row.get("reference_gene", ""))


def _atlas_relation(row: dict[str, object]) -> str:
    return str(
        row.get("orthology_relationship", "")
        or row.get("relationship", "")
        or row.get("orthology_class", "")
    ).strip()


def enforce_semantic_competitor_gate(
    preliminary_family_copy_arg: str,
    original_family_copy_arg: str,
    atlas_arg: str,
    out_dir_arg: str,
    *,
    base_shortlist_arg: str | None = None,
    case_id: str = "",
) -> dict[str, object]:
    """Recompute semantic copy promotion after scanning all machine-equivalent copies."""
    preliminary_path = runtime_path(preliminary_family_copy_arg)
    original_path = runtime_path(original_family_copy_arg)
    atlas_path = runtime_path(atlas_arg)
    preliminary = pd.read_csv(preliminary_path, sep="\t", dtype=str).fillna("")
    original = pd.read_csv(original_path, sep="\t", dtype=str).fillna("")
    atlas = pd.read_csv(atlas_path, sep="\t", dtype=str).fillna("")
    required = {"family_id", "candidate_gene_id", "disambiguation_status"}
    for label, frame in (("preliminary semantic copy", preliminary), ("original family copy", original)):
        missing = sorted(required - set(frame.columns))
        if missing:
            raise ValueError(f"{label} table is missing columns: {', '.join(missing)}")
    if "gene_id" not in atlas.columns:
        raise ValueError("semantic competitor gate atlas requires gene_id")
    if set(preliminary["candidate_gene_id"]) != set(original["candidate_gene_id"]):
        raise ValueError("preliminary semantic copy roster differs from original family-copy roster")

    atlas_by_gene = {str(row["gene_id"]): row for row in atlas.to_dict("records")}
    original_by_key = {
        (str(row["family_id"]), str(row["candidate_gene_id"])): row
        for row in original.to_dict("records")
    }
    preliminary_by_key = {
        (str(row["family_id"]), str(row["candidate_gene_id"])): row
        for row in preliminary.to_dict("records")
    }

    # One or more model-proposed mappings establish which supplied reference(s)
    # need to be checked.  The model does not establish uniqueness.
    reference_hypotheses: dict[str, dict[str, set[str]]] = {}
    for row in preliminary.to_dict("records"):
        if str(row.get("semantic_verification_status", "")) != "verified_copy_mapping":
            continue
        family_id = str(row["family_id"])
        reference = str(row.get("semantic_reference_gene", "")).strip()
        gene_id = str(row["candidate_gene_id"])
        if reference:
            reference_hypotheses.setdefault(family_id, {}).setdefault(reference, set()).add(gene_id)

    safe_promotions: dict[str, set[str]] = {}
    supporting_references: dict[tuple[str, str], set[str]] = {}
    model_mapped: set[tuple[str, str]] = set()
    for family_id, references in reference_hypotheses.items():
        family_original = original[original["family_id"].astype(str).eq(family_id)]
        eligible = {
            str(row["candidate_gene_id"])
            for row in family_original.to_dict("records")
            if str(row.get("disambiguation_status", "")) in ELIGIBLE_UNRESOLVED_STATES
        }
        for reference, proposed_genes in references.items():
            model_mapped.update((family_id, gene_id) for gene_id in proposed_genes)
            strong_competitors: set[str] = set()
            for gene_id in eligible:
                atlas_row = atlas_by_gene.get(gene_id, {})
                original_row = original_by_key.get((family_id, gene_id), {})
                relation_ok = _strong(_atlas_relation(atlas_row)) or _strong(original_row.get("orthology_class", ""))
                if reference in _atlas_references(atlas_row) and relation_ok:
                    strong_competitors.add(gene_id)
                    supporting_references.setdefault((family_id, gene_id), set()).add(reference)
            if proposed_genes & strong_competitors:
                safe_promotions.setdefault(family_id, set()).update(strong_competitors)

    output_columns = list(original.columns) + [
        column for column in FINAL_EXTRA_COLUMNS if column not in original.columns
    ]
    rows: list[dict[str, object]] = []
    competitor_expansion_count = 0
    for original_row in original.to_dict("records"):
        family_id = str(original_row["family_id"])
        gene_id = str(original_row["candidate_gene_id"])
        updated = dict(original_row)
        promoted = sorted(safe_promotions.get(family_id, set()))
        key = (family_id, gene_id)
        if gene_id in promoted:
            updated["disambiguation_status"] = "preferred_copy" if len(promoted) == 1 else "co_preferred_tie"
            references = sorted(supporting_references.get(key, set()))
            preliminary_info = preliminary_by_key.get(key, {})
            if key in model_mapped:
                updated["semantic_match_basis"] = "llm_reference_hypothesis+machine_competitor_scan"
                updated["semantic_verification_status"] = "verified_unique_or_tied_after_competitor_scan"
                updated["semantic_relation"] = str(preliminary_info.get("semantic_relation", ""))
                updated["semantic_proposal_id"] = str(preliminary_info.get("semantic_proposal_id", ""))
            else:
                competitor_expansion_count += 1
                updated["semantic_match_basis"] = "machine_equivalent_reference_competitor"
                updated["semantic_verification_status"] = "verified_co_preferred_competitor"
                updated["semantic_relation"] = "machine_equivalent_to_proposed_reference"
                updated["semantic_proposal_id"] = ""
            updated["semantic_reference_gene"] = ";".join(references)
        else:
            for column in FINAL_EXTRA_COLUMNS:
                updated[column] = ""
        rows.append(updated)
    final_frame = pd.DataFrame(rows, columns=output_columns)

    if base_shortlist_arg:
        base_shortlist = pd.read_csv(runtime_path(base_shortlist_arg), sep="\t", dtype=str).fillna("")
    else:
        base_shortlist = pd.DataFrame(columns=SHORTLIST_COLUMNS)
    for column in SHORTLIST_COLUMNS:
        if column not in base_shortlist.columns:
            base_shortlist[column] = ""
    shortlist = base_shortlist[SHORTLIST_COLUMNS].copy()
    promoted_families = set(safe_promotions)
    if promoted_families and len(shortlist):
        shortlist = shortlist[~shortlist["family_id"].astype(str).isin(promoted_families)].copy()

    additions: list[dict[str, object]] = []
    for family_id, genes in sorted(safe_promotions.items()):
        ordered = sorted(genes)
        for gene_id in ordered:
            base = original_by_key[(family_id, gene_id)]
            additions.append({
                "family_id": family_id,
                "candidate_gene_id": gene_id,
                "disambiguation_status": "preferred_copy" if len(ordered) == 1 else "co_preferred_tie",
                "leading_axes": "semantic_reference_hypothesis_machine_competitor_scan",
                "subfamily": str(base.get("subfamily", "")),
                "orthology_class": str(base.get("orthology_class", "")),
                "first_wet_priority": "first" if len(ordered) == 1 else "co_first",
                "rank_contribution": "none",
                "notes": (
                    "LLM/agent proposed a supplied-reference hypothesis; machine enumerated every "
                    "eligible family member with the same strong reference relation before setting copy priority"
                ),
            })
    if additions:
        shortlist = pd.concat([shortlist, pd.DataFrame(additions, columns=SHORTLIST_COLUMNS)], ignore_index=True)
    shortlist = shortlist.drop_duplicates(subset=["family_id", "candidate_gene_id"], keep="last")

    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    family_copy_path = out / "family_semantic_copy_final.tsv"
    shortlist_path = out / "family_semantic_copy_shortlist_final.tsv"
    final_frame.to_csv(family_copy_path, sep="\t", index=False, lineterminator="\n")
    shortlist.to_csv(shortlist_path, sep="\t", index=False, lineterminator="\n")

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id,
        "reference_hypothesis_count": int(sum(len(value) for value in reference_hypotheses.values())),
        "promoted_family_count": len(safe_promotions),
        "promoted_copy_count": int(sum(len(value) for value in safe_promotions.values())),
        "competitor_expansion_count": competitor_expansion_count,
        "rank_contribution": "none",
        "membership_effect": "none",
        "inputs": {
            "preliminary_family_copy": file_fingerprint(preliminary_path, full_hash=True),
            "original_family_copy": file_fingerprint(original_path, full_hash=True),
            "atlas": file_fingerprint(atlas_path, full_hash=True),
            "base_shortlist": file_fingerprint(runtime_path(base_shortlist_arg), full_hash=True) if base_shortlist_arg else None,
        },
        "outputs": {
            "family_copy": str(family_copy_path),
            "shortlist": str(shortlist_path),
            "manifest": str(out / "family_semantic_competitor_gate_manifest.json"),
        },
        "policy": [
            "LLM/agent proposes a supplied-reference hypothesis but never establishes uniqueness",
            "the machine scans every eligible unresolved family member for the same supplied reference and strong orthology",
            "one machine-equivalent copy becomes preferred_copy; several become co_preferred_tie even when the model omitted competitors",
            "all nonpromoted rows revert to their original deterministic disambiguation_status",
            "no candidate is deleted and rank_contribution/membership_effect remain none",
        ],
    }
    manifest_path = out / "family_semantic_competitor_gate_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return {
        "family_copy": str(family_copy_path),
        "shortlist": str(shortlist_path),
        "manifest": str(manifest_path),
        "promoted_family_count": len(safe_promotions),
        "promoted_copy_count": int(sum(len(value) for value in safe_promotions.values())),
        "competitor_expansion_count": competitor_expansion_count,
    }
