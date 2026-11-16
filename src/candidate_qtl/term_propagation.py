"""Deterministic trait-search-term propagation from reference genomes.

The dominant lever on candidate-triage quality is whether the true gene reaches
a source-backed function match, which depends on the trait search vocabulary.
This module raises that vocabulary *mechanically*, with no LLM: for each QTL
interval candidate it pulls the functional annotation (GO / description terms)
of that candidate's orthologs in the role-labelled reference genomes and
propagates those terms back to the candidate, with full provenance.

Design constraints (issue #130):
- The only assumed inputs are reference genomes: the native assembly plus, when
  available, model-plant / close-relative references (the reference-panel
  roles). A hand-curated per-trait table is optional, never required.
- Work stays scoped to the QTL interval candidates and their orthologs; no
  whole-genome recomputation, so it remains light on a 16 GB-class machine.
- Every propagated term keeps provenance (source ortholog gene, source genome,
  role, ontology id) and stays ``pending`` human review. The LLM is not part of
  this deterministic lane.

Term availability is provenance-carrying vocabulary, not causal evidence.
"""

from __future__ import annotations

import json
from collections import defaultdict

import pandas as pd

from . import schemas
from .core import file_fingerprint, runtime_path


PROPAGATION_VERSION = "1.1.0"
ORTHOLOGY_REQUIRED = {"candidate_gene_id", "reference_gene_id"}
ANNOTATION_REQUIRED = {"reference_gene_id", "term"}
ONTOLOGY_RESOLUTION_REQUIRED = {
    "trait_query_id",
    "resolution_status",
    "ontology_id",
    "ontology_name",
    "ontology_namespace",
    "rank_contribution",
}
RELATIONSHIP_CONFIDENCE = {
    "ortholog": 1.0, "one_to_one": 1.0, "one-to-one": 1.0, "rbh": 1.0, "co_ortholog": 0.8, "co-ortholog": 0.8,
    "one_to_many": 0.6, "many_to_one": 0.6, "many_to_many": 0.5, "paralog": 0.3, "homolog": 0.5,
}


def _normalize(term: str) -> str:
    return " ".join(str(term).strip().casefold().split())


def _load_trait_vocabulary(path_arg: str | None) -> tuple[set[str], set[str], str, int, int]:
    """Load legacy free-text terms or a frozen plant-ontology resolution artifact.

    The Phase-0 ontology resolver (#528/#529) is intentionally consumed through
    the existing ``trait_terms_arg`` surface rather than a new CLI option. Only
    resolved ontology IDs and namespace-unambiguous canonical ontology names
    become matching vocabulary. Parent IDs, raw user text, ambiguous rows and
    unresolved rows remain context/non-evidence and are not imported here.
    """
    if not path_arg:
        return set(), set(), "none", 0, 0
    path = runtime_path(path_arg)
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not lines:
        return set(), set(), "legacy_terms", 0, 0

    header = {column.strip() for column in lines[0].split("\t")}
    looks_like_resolution = bool(
        {"resolution_status", "ontology_id", "ontology_name", "ontology_namespace"} & header
    )
    if looks_like_resolution:
        missing = ONTOLOGY_RESOLUTION_REQUIRED - header
        if missing:
            raise ValueError(
                "ontology resolution table missing columns: " + ", ".join(sorted(missing))
            )
        frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        rank_values = {str(value).strip().casefold() for value in frame["rank_contribution"]}
        if rank_values - {"none"}:
            raise ValueError("ontology resolution rows must retain rank_contribution=none")
        resolved = frame[frame["resolution_status"].astype(str).str.strip().str.casefold().eq("resolved")]
        ontology_ids = {
            str(value).strip()
            for value in resolved["ontology_id"]
            if str(value).strip()
        }
        name_origins: dict[str, set[tuple[str, str]]] = defaultdict(set)
        for record in resolved.to_dict("records"):
            name = _normalize(record.get("ontology_name", ""))
            if not name:
                continue
            name_origins[name].add(
                (
                    str(record.get("ontology_namespace", "")).strip(),
                    str(record.get("ontology_id", "")).strip(),
                )
            )
        canonical_terms = {
            name for name, origins in name_origins.items() if len(origins) == 1
        }
        ambiguous_canonical_names = sum(len(origins) > 1 for origins in name_origins.values())
        return (
            canonical_terms,
            ontology_ids,
            "ontology_resolution",
            int(len(resolved)),
            ambiguous_canonical_names,
        )

    terms: set[str] = set()
    for line in lines:
        token = line.split("\t")[0].strip()
        if token and token.casefold() != "term":
            terms.add(_normalize(token))
    return terms, set(), "legacy_terms", len(terms), 0


def propagate_terms(
    candidates_arg: str,
    orthology_arg: str | None,
    reference_annotation_arg: str | None,
    out_arg: str,
    trait_terms_arg: str | None = None,
    full_hash: bool = False,
) -> dict[str, object]:
    """Propagate reference-ortholog function terms onto QTL interval candidates."""
    candidates_path = runtime_path(candidates_arg)
    out = runtime_path(out_arg)
    out.mkdir(parents=True, exist_ok=True)

    candidates = pd.read_csv(candidates_path, sep="\t", dtype=str).fillna("")
    if "gene_id" not in candidates.columns:
        raise ValueError("candidates table requires a gene_id column")
    if candidates["gene_id"].duplicated().any():
        raise ValueError("candidates table gene_id must be unique")
    candidate_ids = [str(g) for g in candidates["gene_id"]]
    candidate_set = set(candidate_ids)
    (
        trait_terms,
        trait_ontology_ids,
        trait_vocabulary_mode,
        resolved_trait_concepts,
        ambiguous_trait_canonical_names,
    ) = _load_trait_vocabulary(trait_terms_arg)

    skip_reasons: list[str] = []
    orthology = pd.DataFrame(columns=sorted(ORTHOLOGY_REQUIRED))
    if orthology_arg:
        orthology = pd.read_csv(runtime_path(orthology_arg), sep="\t", dtype=str).fillna("")
        missing = ORTHOLOGY_REQUIRED - set(orthology.columns)
        if missing:
            raise ValueError(f"orthology table missing columns: {', '.join(sorted(missing))}")
        orthology = orthology[orthology["candidate_gene_id"].isin(candidate_set)].copy()
    else:
        skip_reasons.append("no_orthology_mapping_native_only")

    annotation = pd.DataFrame(columns=sorted(ANNOTATION_REQUIRED))
    if reference_annotation_arg:
        annotation = pd.read_csv(runtime_path(reference_annotation_arg), sep="\t", dtype=str).fillna("")
        missing = ANNOTATION_REQUIRED - set(annotation.columns)
        if missing:
            raise ValueError(f"reference annotation table missing columns: {', '.join(sorted(missing))}")
    else:
        skip_reasons.append("no_reference_function_annotation_native_only")

    annotation_by_reference: dict[str, list[dict[str, str]]] = defaultdict(list)
    for record in annotation.to_dict("records"):
        reference_gene = str(record["reference_gene_id"]).strip()
        term = str(record["term"]).strip()
        if reference_gene and term:
            annotation_by_reference[reference_gene].append(record)

    rows: list[dict[str, object]] = []
    seen: set[tuple[str, str, str, str, str]] = set()
    candidates_with_ortholog: set[str] = set()
    orthologs_without_function: set[str] = set()
    source_genomes: set[str] = set()
    roles_used: set[str] = set()
    for record in orthology.to_dict("records"):
        candidate_gene = str(record["candidate_gene_id"]).strip()
        reference_gene = str(record["reference_gene_id"]).strip()
        if not candidate_gene or not reference_gene:
            continue
        candidates_with_ortholog.add(candidate_gene)
        relationship = str(record.get("relationship", "")).strip()
        role = str(record.get("reference_role", "")).strip() or "unspecified"
        if role not in schemas.TERM_REFERENCE_ROLES:
            role = "unspecified"
        roles_used.add(role)
        source_genome = str(record.get("source_genome", "")).strip()
        species = str(record.get("reference_species", "")).strip()
        if source_genome:
            source_genomes.add(source_genome)
        declared = str(record.get("confidence", "")).strip()
        base_weight = RELATIONSHIP_CONFIDENCE.get(relationship.casefold(), 0.5)
        try:
            confidence = min(base_weight, float(declared)) if declared else base_weight
        except ValueError:
            confidence = base_weight
        reference_terms = annotation_by_reference.get(reference_gene, [])
        if not reference_terms:
            orthologs_without_function.add(candidate_gene)
            continue
        for term_record in reference_terms:
            term = str(term_record["term"]).strip()
            normalized = _normalize(term)
            if not normalized:
                continue
            annotation_ontology_id = str(term_record.get("ontology_id", "")).strip()
            key = (candidate_gene, normalized, reference_gene, source_genome, annotation_ontology_id)
            if key in seen:
                continue
            seen.add(key)
            annotation_ontology_name = _normalize(term_record.get("ontology_name", ""))
            if not trait_terms and not trait_ontology_ids:
                relevance = "not_evaluated"
            elif trait_vocabulary_mode == "ontology_resolution":
                if annotation_ontology_id:
                    relevance = "trait_matched" if annotation_ontology_id in trait_ontology_ids else "unrelated"
                else:
                    relevance = (
                        "trait_matched"
                        if normalized in trait_terms or annotation_ontology_name in trait_terms
                        else "unrelated"
                    )
            elif normalized in trait_terms or any(t in normalized or normalized in t for t in trait_terms):
                relevance = "trait_matched"
            else:
                relevance = "unrelated"
            rows.append({
                "candidate_gene_id": candidate_gene,
                "term": term,
                "term_normalized": normalized,
                "ontology_id": annotation_ontology_id,
                "ontology_name": str(term_record.get("ontology_name", "")).strip(),
                "reference_gene_id": reference_gene,
                "source_genome": source_genome,
                "reference_role": role,
                "reference_species": species,
                "relationship": relationship,
                "term_source": str(term_record.get("source", "")).strip() or "reference_annotation",
                "match_type": "ortholog_propagated",
                "term_confidence": round(confidence, 6),
                "trait_relevance": relevance,
                "review_status": "pending",
                "notes": "" if reference_terms else "ortholog_without_reference_function",
            })

    frame = pd.DataFrame(rows, columns=schemas.PROPAGATED_TERM_COLUMNS).sort_values(
        ["candidate_gene_id", "term_normalized", "reference_gene_id", "source_genome"], kind="stable")
    schemas.validate(frame.fillna(""), schemas.PROPAGATED_TERM_SCHEMA)
    terms_path = out / "propagated_terms.tsv"
    frame.to_csv(terms_path, sep="\t", index=False, lineterminator="\n")

    expansion_rows: list[dict[str, object]] = []
    for (candidate_gene, normalized), group in frame.groupby(["candidate_gene_id", "term_normalized"], sort=True):
        expansion_rows.append({
            "candidate_gene_id": candidate_gene,
            "term_normalized": normalized,
            "term": group["term"].iloc[0],
            "ontology_ids": ";".join(sorted({o for o in group["ontology_id"] if o})),
            "source_genomes": ";".join(sorted({s for s in group["source_genome"] if s})),
            "reference_roles": ";".join(sorted(set(group["reference_role"]))),
            "reference_genes": ";".join(sorted(set(group["reference_gene_id"]))),
            "max_confidence": round(float(pd.to_numeric(group["term_confidence"], errors="coerce").max()), 6),
            "trait_relevance": "trait_matched" if (group["trait_relevance"] == "trait_matched").any() else group["trait_relevance"].iloc[0],
            "review_status": "pending",
        })
    expansion = pd.DataFrame(
        expansion_rows,
        columns=[
            "candidate_gene_id", "term_normalized", "term", "ontology_ids", "source_genomes",
            "reference_roles", "reference_genes", "max_confidence", "trait_relevance", "review_status",
        ],
    )
    expansion_path = out / "search_term_expansion.tsv"
    expansion.sort_values(["candidate_gene_id", "term_normalized"], kind="stable").to_csv(
        expansion_path, sep="\t", index=False, lineterminator="\n"
    )

    covered = sorted({str(g) for g in frame["candidate_gene_id"]}) if not frame.empty else []
    manifest = {
        "schema_version": "1.0",
        "propagation_version": PROPAGATION_VERSION,
        "interval_scoped": True,
        "candidate_count": len(candidate_ids),
        "candidates_with_ortholog": len(candidates_with_ortholog),
        "candidates_with_propagated_terms": len(covered),
        "term_coverage_fraction": round(len(covered) / len(candidate_ids), 6) if candidate_ids else 0.0,
        "orthologs_without_reference_function": sorted(orthologs_without_function - set(covered)),
        "candidates_without_ortholog": sorted(candidate_set - candidates_with_ortholog),
        "propagated_term_rows": len(frame),
        "distinct_terms": int(frame["term_normalized"].nunique()) if not frame.empty else 0,
        "trait_matched_terms": int((frame["trait_relevance"] == "trait_matched").sum()) if not frame.empty else 0,
        "trait_terms_provided": bool(trait_terms or trait_ontology_ids),
        "trait_vocabulary_mode": trait_vocabulary_mode,
        "resolved_trait_concept_count": resolved_trait_concepts,
        "trait_ontology_id_count": len(trait_ontology_ids),
        "ambiguous_trait_canonical_name_count": ambiguous_trait_canonical_names,
        "source_genomes": sorted(source_genomes),
        "reference_roles_used": sorted(roles_used),
        "skip_reasons": skip_reasons,
        "native_only": bool(skip_reasons and frame.empty),
        "fingerprint_mode": "strict_full_hash" if full_hash else "size_tiered",
        "inputs": {
            "candidates": file_fingerprint(candidates_path, full_hash),
            "orthology": file_fingerprint(runtime_path(orthology_arg), full_hash) if orthology_arg else {"status": "not_provided"},
            "reference_annotation": file_fingerprint(runtime_path(reference_annotation_arg), full_hash) if reference_annotation_arg else {"status": "not_provided"},
            "trait_terms": file_fingerprint(runtime_path(trait_terms_arg), full_hash) if trait_terms_arg else {"status": "not_provided"},
        },
        "scientific_limits": [
            "propagated terms are provenance-carrying search vocabulary, not causal evidence",
            "every term stays pending human review; this lane uses no LLM",
            "work is scoped to interval candidates and their orthologs; no whole-genome recomputation",
            "native-only input degrades to an empty propagation with a machine-readable skip reason",
            "frozen ontology resolution uses only resolved ontology IDs and namespace-unambiguous canonical names; parents, raw text, ambiguous and unresolved concepts are not promoted",
            "when a reference annotation carries an explicit ontology ID, only exact ID equality can trait-match; canonical-name fallback is allowed only when that ID is absent",
            "a canonical label resolved to multiple ontology namespaces/IDs is excluded from name-only fallback rather than flattening namespace identity",
            "same display term with distinct ontology IDs remains distinct provenance until the candidate-level expansion ledger aggregates it",
            "ontology matching changes trait_relevance only; term_confidence remains orthology-relationship-derived",
        ],
    }
    manifest_path = out / "term_propagation_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "propagated_terms": str(terms_path),
        "search_term_expansion": str(expansion_path),
        "manifest": str(manifest_path),
        "term_coverage_fraction": manifest["term_coverage_fraction"],
        "distinct_terms": manifest["distinct_terms"],
        "trait_matched_terms": manifest["trait_matched_terms"],
        "native_only": manifest["native_only"],
    }
