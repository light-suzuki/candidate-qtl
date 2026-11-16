"""Build a human-review table linking QTL candidates to trait-relevant reference genes.

This command deliberately does not infer a trait association from a BLAST hit.
The reference-gene panel is a source-backed, analyst-curated input; orthology
tables supply candidate-to-reference mappings from one or more providers.  The
output keeps every mapping and every unresolved candidate visible for review.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from .core import file_fingerprint, runtime_path


PANEL_REQUIRED = {"reference_gene_id", "source_id"}
ORTHOLOGY_REQUIRED = {"candidate_gene_id", "reference_gene_id", "relationship"}


def _provider_path(spec: str) -> tuple[str, Path]:
    if "=" not in spec:
        raise ValueError("--orthology must use PROVIDER=PATH")
    provider, raw_path = spec.split("=", 1)
    if not provider.strip() or not raw_path.strip():
        raise ValueError("--orthology must use non-empty PROVIDER=PATH")
    return provider.strip(), runtime_path(raw_path)


def _read_tsv(path: Path, required: set[str], label: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"{label} {path} is missing required columns: {', '.join(missing)}")
    return frame


def reference_candidate_table(
    candidates_arg: str,
    reference_genes_arg: str,
    orthology_args: list[str],
    out_dir_arg: str,
    annotation_arg: str | None = None,
) -> dict[str, object]:
    """Join a curated reference-gene panel to candidate orthology mappings.

    Multiple provider calls count as mapping robustness only.  They never
    become multiple independent biological evidence families.
    """

    candidates_path = runtime_path(candidates_arg)
    panel_path = runtime_path(reference_genes_arg)
    out_dir = runtime_path(out_dir_arg)
    candidates = _read_tsv(candidates_path, {"gene_id"}, "candidate table")
    panel = _read_tsv(panel_path, PANEL_REQUIRED, "reference-gene panel")
    if panel["reference_gene_id"].str.strip().eq("").any():
        raise ValueError("reference-gene panel contains an empty reference_gene_id")
    if panel["source_id"].str.strip().eq("").any():
        raise ValueError("reference-gene panel requires a source_id for every reference gene")

    panel = panel.drop_duplicates("reference_gene_id", keep="first").copy()
    panel["reference_gene_id"] = panel["reference_gene_id"].str.strip()
    candidate_ids = set(candidates["gene_id"].astype(str))

    mapping_frames: list[pd.DataFrame] = []
    input_fingerprints: dict[str, object] = {
        "candidates": file_fingerprint(candidates_path, False),
        "reference_genes": file_fingerprint(panel_path, False),
    }
    for spec in orthology_args:
        provider, path = _provider_path(spec)
        mapping = _read_tsv(path, ORTHOLOGY_REQUIRED, f"orthology provider {provider}")
        mapping["candidate_gene_id"] = mapping["candidate_gene_id"].str.strip()
        mapping["reference_gene_id"] = mapping["reference_gene_id"].str.strip()
        mapping = mapping[
            mapping["candidate_gene_id"].isin(candidate_ids)
            & mapping["reference_gene_id"].isin(set(panel["reference_gene_id"]))
        ].copy()
        mapping["provider"] = provider
        mapping["mapping_source_id"] = mapping.get("source_id", pd.Series("", index=mapping.index)).replace("", str(path))
        mapping_frames.append(mapping)
        input_fingerprints[f"orthology:{provider}"] = file_fingerprint(path, False)

    mappings = pd.concat(mapping_frames, ignore_index=True, sort=False) if mapping_frames else pd.DataFrame()
    annotation = pd.DataFrame(columns=["gene_id", "annotation"])
    if annotation_arg:
        annotation_path = runtime_path(annotation_arg)
        annotation = _read_tsv(annotation_path, {"gene_id"}, "annotation table")
        annotation = annotation.drop_duplicates("gene_id", keep="first")
        input_fingerprints["annotation"] = file_fingerprint(annotation_path, False)

    candidate_columns = [column for column in ("gene_id", "gene_name", "chrom", "start", "end", "dry_group", "screening_rank") if column in candidates.columns]
    candidate_view = candidates[candidate_columns].rename(columns={"gene_id": "candidate_gene_id", "gene_name": "candidate_gene_name"})
    if not mappings.empty:
        provider_counts = mappings.groupby(["candidate_gene_id", "reference_gene_id"])["provider"].agg(lambda values: ";".join(sorted(set(values)))).reset_index(name="mapping_providers")
        provider_counts["mapping_provider_count"] = provider_counts["mapping_providers"].str.split(";").str.len()
        relationships = mappings.groupby(["candidate_gene_id", "reference_gene_id"])["relationship"].agg(lambda values: ";".join(sorted(set(v for v in values if v)))).reset_index()
        mapping_sources = mappings.groupby(["candidate_gene_id", "reference_gene_id"])["mapping_source_id"].agg(lambda values: ";".join(sorted(set(v for v in values if v)))).reset_index()
        table = provider_counts.merge(relationships, on=["candidate_gene_id", "reference_gene_id"]).merge(mapping_sources, on=["candidate_gene_id", "reference_gene_id"])
        table = table.merge(panel, on="reference_gene_id", how="left").merge(candidate_view, on="candidate_gene_id", how="left")
    else:
        table = pd.DataFrame(columns=["candidate_gene_id", "reference_gene_id"])

    mapped_ids = set(table.get("candidate_gene_id", pd.Series(dtype=str)))
    unresolved = candidate_view[~candidate_view["candidate_gene_id"].isin(mapped_ids)].copy()
    unresolved["reference_gene_id"] = ""
    unresolved["mapping_providers"] = ""
    unresolved["mapping_provider_count"] = 0
    unresolved["relationship"] = ""
    unresolved["mapping_source_id"] = ""
    unresolved["mapping_status"] = "no_reference_panel_match"
    unresolved["review_required"] = True
    if not table.empty:
        relationship_text = table["relationship"].str.casefold()
        table["mapping_status"] = "candidate_reference_mapping"
        table["review_required"] = relationship_text.str.contains("one_to_many|many_to_many|ambiguous|paralog", regex=True, na=False) | table["mapping_provider_count"].eq(0)
    table = pd.concat([table, unresolved], ignore_index=True, sort=False)

    # An LLM-proposed panel row is a research lead until a person confirms that
    # the source really supports the stated trait relation.  Manual curation can
    # use the same table and set review_status=reviewed.
    if "review_status" in table.columns:
        pending_panel_review = table["review_status"].astype(str).str.casefold().ne("reviewed") & table["reference_gene_id"].astype(str).ne("")
        table.loc[pending_panel_review, "review_required"] = True

    if not annotation.empty:
        annotation = annotation.rename(columns={"gene_id": "candidate_gene_id"})
        table = table.merge(annotation, on="candidate_gene_id", how="left")
    if "annotation" not in table.columns:
        table["annotation"] = ""
    table["annotation"] = table["annotation"].fillna("")

    leading = [
        "candidate_gene_id", "candidate_gene_name", "chrom", "start", "end", "dry_group", "screening_rank",
        "reference_gene_id", "reference_gene_name", "reference_species", "trait_relation", "evidence_level",
        "source_id", "discovery_method", "review_status", "mapping_status", "relationship", "mapping_providers", "mapping_provider_count",
        "mapping_source_id", "annotation", "review_required",
    ]
    for column in leading:
        if column not in table.columns:
            table[column] = ""
    table = table[leading + [column for column in table.columns if column not in leading]]
    table = table.sort_values(["candidate_gene_id", "reference_gene_id", "mapping_providers"], kind="stable")

    out_dir.mkdir(parents=True, exist_ok=True)
    table_path = out_dir / "reference_candidate_table.tsv"
    review_path = out_dir / "reference_candidate_review.tsv"
    manifest_path = out_dir / "reference_candidate_manifest.json"
    table.to_csv(table_path, sep="\t", index=False)
    table[table["review_required"].astype(str).str.casefold().eq("true")].to_csv(review_path, sep="\t", index=False)
    manifest = {
        "reference_candidate_table": str(table_path),
        "review_table": str(review_path),
        "candidate_count": int(candidates["gene_id"].nunique()),
        "reference_gene_count": int(panel["reference_gene_id"].nunique()),
        "mapped_candidate_count": int(len(mapped_ids)),
        "unmapped_candidate_count": int(candidates[~candidates["gene_id"].isin(mapped_ids)]["gene_id"].nunique()),
        "review_count": int(table["review_required"].astype(str).str.casefold().eq("true").sum()),
        "input_fingerprints": input_fingerprints,
        "principle": "An LLM may propose the reference-gene panel, but every retained row needs a stable source ID and remains review_required until review_status=reviewed. Candidate-to-reference mapping and biological evidence are shown separately; multiple mapping providers indicate computational robustness only and do not create independent biological evidence.",
    }
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return {"table": str(table_path), "review": str(review_path), "manifest": str(manifest_path), "mapped_candidate_count": int(len(mapped_ids))}
