"""Phenotype-independent interval atlas (issue #168).

The workflow is strong at recovering interval candidates *once* a correct
pathway / gene-family / homolog bait is supplied. The largest remaining
uncertainty sits in front of that:

    phenotype -> which pathway? -> which gene family? -> which reference bait?

That step depends on a researcher (or an LLM) already knowing what to look for.
This module reduces the dependence by inverting the order: it first structures
**every gene in the QTL interval** into families, domain architectures,
orthogroups and broad functional bins, so hypotheses can be formed by *reading
the finite set of options actually present in the interval* rather than by
guessing the right bait up front.

Deliberate boundaries:

* The atlas is **phenotype-independent by construction**. No phenotype input
  reaches this module; overlaying one is a later, separate stage. The manifest
  records the atlas fingerprint so a later overlay cannot silently change how the
  interval looked before the hypothesis existed.
* Classification is **an exploration index, not causal evidence**. Nothing here
  writes a score, and family membership never promotes a gene.
* ``unknown`` / ``conflicting`` / ``dark`` genes are kept as first-class buckets.
  A gene is never dropped for lacking annotation.
* Broad-function bins come from a versioned, human-editable resource, not from a
  keyword list buried in code, and every assignment records which pattern and
  which source produced it.
* Several reference species agreeing through the same annotation chain is one
  piece of evidence, not several; the atlas records providers without voting.
* A gene with no family / orthology / domain signal is ``no_detectable_homology``,
  never ``lineage_restricted_supported``: a failed search is not evidence of
  lineage restriction (issue #187). Lineage claims require explicit
  multi-reference negative-search coverage, which the current inputs cannot
  supply; search coverage and claim basis are recorded in the manifest.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import pandas as pd

from .core import file_fingerprint, runtime_path
from .family_index import symbol_stem

SCHEMA_VERSION = "1.4"
BIN_COLUMNS = {"bin_id", "bin_label", "match_type", "pattern"}
MATCH_TYPES = ("go", "keyword")

# Issue #374 T3: gene-model completeness per interval gene, derived from the
# candidate table's own columns only. A phenotype never reaches this stage, so
# the status can back the manifest claim that nothing was filtered on phenotype
# before model completeness was assessed.
GENE_MODEL_STATUSES = ("full_model", "partial_model", "no_model")
# translation_status values that positively state there is no CDS at all
# (reference_extract.write_candidates emits "no_cds" when the gene has no CDS
# segments). Every other translation_status this repo produces
# (not_translated_invalid_phase / not_translated_length_not_multiple_of_3 /
# translated / translated_complete_orf / translated_partial_or_unverified_orf)
# is only ever set once a CDS sequence exists, so it positively indicates a
# model is present (issue #374 T3 F3): only "no_cds" is a negative claim.
GENE_MODEL_NO_CDS_TRANSLATION_STATUSES = frozenset({"no_cds"})

# Issue #184: family layers are kept as separate concepts. family_id /
# family_label / family_source remain as legacy aliases (their value is the
# family-index coarse value: a symbol-stem bucket or a curated family id
# depending on input); discovery_bucket_* / evolutionary_family_* /
# orthogroup / copy_set_* carry the separated semantics.
ATLAS_COLUMNS = [
    "gene_id", "chrom", "start", "end", "protein_length",
    "local_annotation_description", "local_annotation_go_terms",
    "family_id", "family_label", "family_source", "family_members_in_interval",
    "discovery_bucket_id", "discovery_bucket_label", "discovery_bucket_method",
    "evolutionary_family_id", "evolutionary_family_source", "evolutionary_family_status",
    "domain_ids", "domain_architecture", "domain_source",
    "family_index_best_subject", "family_index_subjects", "family_index_subject_endpoints",
    "family_index_assignment_status", "family_index_best_identity", "family_index_best_query_coverage",
    "family_index_best_evalue", "family_index_best_bitscore", "family_index_second_bitscore",
    "family_index_top2_bitscore_ratio", "family_index_symbol_source_rank",
    "family_index_symbol_source_bitscore_ratio",
    "orthogroup", "orthology_relationship", "orthology_reference_gene", "orthology_providers",
    "subfamily_id", "subfamily_source",
    "focal_paralog_count", "tandem_expansion",
    "broad_function_bins", "primary_broad_function", "broad_function_source", "broad_function_evidence",
    "conservation_class", "evidence_completeness", "annotation_rescue_needed",
    "copy_set_id", "copy_set_status",
    "classification_basis", "rank_contribution", "notes",
    "gene_model_completeness",
]
FAMILY_SUMMARY_COLUMNS = ["family_id", "family_label", "member_count", "member_gene_ids", "tandem_blocks",
                          "orthology_one_to_one", "domain_supported", "notes"]
DOMAIN_SUMMARY_COLUMNS = ["domain_architecture", "member_count", "member_gene_ids", "domain_ids", "notes"]
ORTHOGROUP_SUMMARY_COLUMNS = ["orthogroup", "member_count", "member_gene_ids", "relationships", "providers", "notes"]
FUNCTION_SUMMARY_COLUMNS = ["bin_id", "bin_label", "member_count", "member_gene_ids", "assignment_sources", "notes"]
UNCLASSIFIED_COLUMNS = ["gene_id", "chrom", "start", "end", "protein_length", "evidence_completeness",
                        "why_unclassified", "annotation_rescue_needed", "required_next_action"]
REVIEW_QUERY_COLUMNS = ["query_id", "scope", "scope_id", "scope_label", "member_count", "member_gene_ids",
                        "suggested_query", "query_type", "rationale", "phenotype_used"]

_CONSERVATION = ("conserved_cross_species", "relative_enriched", "no_detectable_homology",
                 "unresolved", "lineage_restricted_supported")


def default_bin_map_path() -> Path:
    """Broad-function bin map used when none is supplied.

    The authoritative, human-editable copy lives at the repo-root
    ``templates/broad_function_bins.tsv``; a byte-identical copy is vendored into
    the package so a clean wheel install can build an atlas without the source
    tree (a drift test keeps them in sync), mirroring the weight registry.
    """
    here = Path(__file__).resolve()
    source_tree = here.parents[2] / "templates" / "broad_function_bins.tsv"
    if source_tree.exists():
        return source_tree
    return here.parent / "resources" / "broad_function_bins.tsv"


def load_bin_map(path_arg: str | None = None) -> tuple[pd.DataFrame, dict[str, object]]:
    path = runtime_path(path_arg) if path_arg else default_bin_map_path()
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    missing = sorted(BIN_COLUMNS - set(frame.columns))
    if missing:
        raise ValueError(f"broad-function bin map is missing columns: {', '.join(missing)}")
    bad = sorted(set(frame["match_type"]) - set(MATCH_TYPES))
    if bad:
        raise ValueError(f"unknown match_type(s): {', '.join(bad)}; allowed {', '.join(MATCH_TYPES)}")
    meta = {"path": str(path), "rule_count": int(len(frame)),
            "bin_count": int(frame["bin_id"].nunique()),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    return frame, meta


# --------------------------------------------------------------------------- #
# optional input tables
# --------------------------------------------------------------------------- #

def _read_optional(path_arg: str | None, required: set[str], label: str) -> pd.DataFrame:
    if not path_arg:
        return pd.DataFrame(columns=sorted(required))
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{label} table is missing columns: {', '.join(missing)}")
    return frame


def _domain_records(frame: pd.DataFrame) -> dict[str, dict[str, object]]:
    records: dict[str, dict[str, object]] = {}
    if frame.empty:
        return records
    order_column = "start" if "start" in frame.columns else None
    for gene_id, group in frame.groupby("gene_id", sort=True):
        if order_column:
            group = group.assign(_order=pd.to_numeric(group[order_column], errors="coerce").fillna(0)).sort_values("_order", kind="stable")
        ids = [str(value) for value in group["domain_id"] if str(value)]
        names = [str(value) for value in group.get("domain_name", group["domain_id"]) if str(value)]
        records[str(gene_id)] = {
            "domain_ids": ";".join(ids),
            # Architecture keeps domain order, which is what discriminates
            # otherwise-similar family members later.
            "domain_architecture": "-".join(names or ids),
            "domain_source": ";".join(sorted({str(value) for value in group.get("source", pd.Series(dtype=str)) if str(value)})) or "supplied_table",
        }
    return records


def _orthology_records(frame: pd.DataFrame) -> dict[str, dict[str, object]]:
    records: dict[str, dict[str, object]] = {}
    if frame.empty:
        return records
    for gene_id, group in frame.groupby("gene_id", sort=True):
        relationships = sorted({str(value) for value in group.get("relationship", pd.Series(dtype=str)) if str(value)})
        providers = sorted({str(value) for value in group.get("provider", pd.Series(dtype=str)) if str(value)})
        references = sorted({str(value) for value in group.get("reference_gene_id", pd.Series(dtype=str)) if str(value)})
        orthogroups = sorted({str(value) for value in group.get("orthogroup", pd.Series(dtype=str)) if str(value)})
        records[str(gene_id)] = {
            "orthogroup": ";".join(orthogroups),
            "orthology_relationship": ";".join(relationships),
            "orthology_reference_gene": ";".join(references),
            # Providers are recorded, never counted as independent votes.
            "orthology_providers": ";".join(providers),
        }
    return records


def _family_subject_endpoint(subject: object) -> str:
    """Return only a syntactically declared UniProt accession wrapper payload.

    Other subject identifiers stay byte-for-byte exact. This is deliberately
    not a symbol, entry-name, family, or alias resolver.
    """
    raw = str(subject or "").strip()
    fields = raw.split("|")
    if len(fields) == 3 and fields[0] in {"sp", "tr"} and fields[1].strip() and fields[2].strip():
        return fields[1].strip()
    return raw


def _family_index_reference_endpoints(row: dict[str, object]) -> str:
    raw_subjects = [str(row.get("best_subject", "") or "")]
    raw_subjects.extend(str(row.get("family_subjects", "") or "").split(";"))
    endpoints = {_family_subject_endpoint(value) for value in raw_subjects if str(value).strip()}
    return ";".join(sorted(endpoints))


def _bin_assignments(description: str, go_terms: str, bins: pd.DataFrame) -> tuple[list[str], list[str], list[str]]:
    """Return (bin_ids, sources, evidence) for one gene."""
    lowered = description.casefold()
    terms = {value.strip() for value in re.split(r"[;,|\s]+", go_terms) if value.strip()}
    assigned: dict[str, tuple[str, str]] = {}
    for rule in bins.to_dict("records"):
        bin_id = str(rule["bin_id"])
        pattern = str(rule["pattern"])
        if str(rule["match_type"]) == "go":
            if pattern in terms and assigned.get(bin_id, ("", ""))[0] != "go_term":
                assigned[bin_id] = ("go_term", pattern)
        elif pattern and pattern.casefold() in lowered and bin_id not in assigned:
            assigned[bin_id] = ("keyword", pattern)
    if not assigned:
        return [], [], []
    # A GO-backed assignment outranks a keyword one when reporting the primary bin.
    ordered = sorted(assigned.items(), key=lambda item: (0 if item[1][0] == "go_term" else 1, item[0]))
    return ([key for key, _ in ordered],
            [value[0] for _, value in ordered],
            [f"{value[0]}:{value[1]}" for _, value in ordered])


def _conservation_class(orthology: dict[str, object], family_id: str, domain_ids: str) -> str:
    relationships = str(orthology.get("orthology_relationship", ""))
    references = str(orthology.get("orthology_reference_gene", ""))
    orthogroup = str(orthology.get("orthogroup", ""))
    if "one-to-one" in relationships or "one_to_one" in relationships:
        return "conserved_cross_species"
    if references:
        return "relative_enriched"
    if family_id or orthogroup or domain_ids:
        # Some signal exists (family / orthogroup / domain hit) but it does not
        # support an enrichment or lineage-restriction claim on its own.
        return "unresolved"
    # No family, orthology or domain signal. A failed search is not evidence of
    # lineage restriction, so the safe default is no_detectable_homology.
    return "no_detectable_homology"


def _family_index_search_coverage(family_index_arg: str | None) -> dict[str, object] | None:
    """Carry the family-index manifest's search provenance into the atlas manifest."""
    if not family_index_arg:
        return None
    manifest_path = runtime_path(family_index_arg).parent / "family_index_manifest.json"
    try:
        meta = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {"search_provenance": "family_index_manifest.json missing or unreadable; only the hit table was reused"}
    return {
        "search_provenance": "family_index_manifest.json",
        "index_scope": str(meta.get("index_scope") or "legacy_unspecified"),
        "reusable_as_genome_wide_index": str(meta.get("index_scope") or "reference_genome") == "reference_genome",
        "reference_species_searched": str(meta.get("reference_species") or ""),
        "reference_label": str(meta.get("reference_label") or ""),
        "search_method": meta.get("family_definition"),
        "thresholds": meta.get("quality_thresholds"),
        "searched_gene_count": int(meta.get("gene_count", 0)),
        "dark_gene_count": int(meta.get("dark_count", 0)),
    }


def _completeness(family_id: str, domains: dict[str, object], orthology: dict[str, object],
                  dark: bool, conflicting: bool) -> str:
    if dark:
        return "dark_or_short_orf"
    if conflicting:
        return "conflicting_annotation"
    if family_id and orthology.get("orthogroup"):
        return "high_confidence_family_assignment"
    if family_id:
        return "family_assignment"
    if domains.get("domain_ids"):
        return "domain_only_assignment"
    if orthology.get("orthogroup") or orthology.get("orthology_reference_gene"):
        return "orthogroup_only_assignment"
    return "no_informative_annotation"


def _gene_model_status(candidate: dict[str, object]) -> str:
    """full_model / partial_model / no_model from the candidate table alone.

    The status is derived only from columns the candidate table itself carries,
    never from the phenotype and never from the family / domain / orthology
    providers, so it stays valid for the B0 completeness claim even on a bare
    run. A missing column simply contributes nothing: with only
    ``protein_length`` present the status degrades to the full_model /
    no_model binary.

    Issue #374 T3 F3: ``cds_length`` / ``translation_status`` are interpreted,
    not merely tested for non-emptiness -- a non-empty string that positively
    states there is no model (``cds_length == 0``, or ``translation_status``
    on the ``GENE_MODEL_NO_CDS_TRANSLATION_STATUSES`` negative vocabulary,
    e.g. ``no_cds``) must read as ``no_model``, not ``partial_model``. A
    positive ``cds_length`` (> 0) or any other non-empty ``translation_status``
    upgrades an undeterminable protein length to ``partial_model``. Unknown /
    malformed values -- a non-numeric ``cds_length``, an unrecognized
    ``translation_status`` -- are never silently read as zero/absent: they
    keep ``partial_model`` so the ambiguity stays visible instead of being
    asserted as "no model". Both columns empty (or absent) stays ``no_model``,
    exactly as before.
    """
    protein_length = str(candidate.get("protein_length", "")).strip()
    if protein_length.isdigit() and int(protein_length) >= 1:
        return "full_model"
    cds_length = str(candidate.get("cds_length", "")).strip()
    translation_status = str(candidate.get("translation_status", "")).strip().casefold()
    if not cds_length and not translation_status:
        return "no_model"
    if translation_status in GENE_MODEL_NO_CDS_TRANSLATION_STATUSES:
        return "no_model"
    if cds_length.isdigit():
        return "no_model" if int(cds_length) == 0 else "partial_model"
    # cds_length is present but not a clean non-negative integer (malformed),
    # or translation_status is non-empty and outside the negative vocabulary
    # (unrecognized, or a real positive status such as translated /
    # translated_complete_orf / not_translated_invalid_phase, all of which
    # only ever accompany a CDS): a model likely exists, so the missingness
    # stays visible as partial_model rather than being collapsed to no_model.
    return "partial_model"


def interval_atlas(
    candidates_arg: str,
    out_dir_arg: str,
    family_index_arg: str | None = None,
    domains_arg: str | None = None,
    orthology_arg: str | None = None,
    annotations_arg: str | None = None,
    flags_arg: str | None = None,
    function_bins_arg: str | None = None,
    tandem_distance: int = 100_000,
    dark_max_aa: int = 200,
    case_id: str = "",
) -> dict[str, object]:
    """Structure every interval gene without using any phenotype information."""
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    candidates = pd.read_csv(runtime_path(candidates_arg), sep="\t", dtype=str).fillna("")
    if "gene_id" not in candidates.columns:
        raise ValueError("candidate table requires gene_id")
    if candidates["gene_id"].duplicated().any():
        raise ValueError("candidate table requires unique gene_id")

    family = pd.read_csv(runtime_path(family_index_arg), sep="\t", dtype=str).fillna("") if family_index_arg else pd.DataFrame()
    family_records = {str(row["gene_id"]): row for row in family.to_dict("records")} if not family.empty else {}
    domains = _domain_records(_read_optional(domains_arg, {"gene_id", "domain_id"}, "domain"))
    orthology = _orthology_records(_read_optional(orthology_arg, {"gene_id"}, "orthology"))
    annotations = _read_optional(annotations_arg, {"gene_id"}, "annotation")
    annotation_records = {str(row["gene_id"]): row for row in annotations.to_dict("records")} if not annotations.empty else {}
    flags = _read_optional(flags_arg, {"gene_id"}, "candidate flags")
    flag_records = {str(row["gene_id"]): row for row in flags.to_dict("records")} if not flags.empty else {}
    bins, bin_meta = load_bin_map(function_bins_arg)

    interval_gene_ids = set(candidates["gene_id"].astype(str))
    family_in_interval = family[family["gene_id"].astype(str).isin(interval_gene_ids)] if not family.empty else pd.DataFrame()
    family_index_meta = _family_index_search_coverage(family_index_arg)
    orthology_providers = sorted({value for record in orthology.values()
                                  for value in str(record.get("orthology_providers", "")).split(";") if value})
    no_hit_counts = {
        "family_index": int(family_in_interval["assignment_status"].ne("assigned").sum())
        if family_index_arg and "assignment_status" in family_in_interval.columns else None,
        "orthology": int(sum(1 for gene_id in interval_gene_ids if gene_id not in orthology)) if orthology_arg else None,
        "domains": int(sum(1 for gene_id in interval_gene_ids if gene_id not in domains)) if domains_arg else None,
    }

    rows: list[dict[str, object]] = []
    for candidate in candidates.to_dict("records"):
        gene_id = str(candidate["gene_id"])
        family_row = family_records.get(gene_id, {})
        domain = domains.get(gene_id, {})
        ortholog = orthology.get(gene_id, {})
        annotation = annotation_records.get(gene_id, {})
        flag = flag_records.get(gene_id, {})

        assigned = str(family_row.get("assignment_status", "")) == "assigned"
        family_id = str(family_row.get("best_subject_stem", "")) if assigned else ""
        family_label = str(family_row.get("best_subject_symbol", "")) if assigned else ""
        family_source = "family_index" if assigned else ""
        # Issue #184 separated family layers. Old family indexes (without the
        # new columns) fall back to the legacy coarse value for the bucket.
        discovery_bucket_id = str(family_row.get("discovery_bucket_id", "") or family_id)
        discovery_bucket_label = family_label
        discovery_bucket_method = str(family_row.get("discovery_bucket_method", "") or family_source)
        evolutionary_family_id = str(family_row.get("evolutionary_family_id", ""))
        evolutionary_family_source = str(family_row.get("evolutionary_family_source", ""))
        evolutionary_family_status = str(family_row.get("evolutionary_family_status", "")
                                        or ("unresolved" if assigned else ""))

        description = " ".join(str(value) for value in (
            annotation.get("description", ""), annotation.get("function_description", ""),
            annotation.get("protein_family", ""), candidate.get("description", "")) if value)
        go_terms = str(annotation.get("go_terms", "") or candidate.get("go_terms", ""))
        bin_ids, bin_sources, bin_evidence = _bin_assignments(description, go_terms, bins)

        # Two independent naming sources disagreeing is recorded, not silently resolved.
        reference_symbol = str(ortholog.get("orthology_reference_gene", "")).split(";")[0]
        reference_stem = symbol_stem(reference_symbol) if reference_symbol else ""
        conflicting = bool(family_id and reference_stem and family_id.casefold() != reference_stem.casefold())

        protein_length = str(candidate.get("protein_length", "") or family_row.get("protein_length", ""))
        dark_state = str(flag.get("dark_state", "") or candidate.get("dark_state", ""))
        short_orf = bool(protein_length.isdigit() and int(protein_length) < dark_max_aa)
        dark = bool(dark_state.startswith("dark")) or (short_orf and not assigned)

        completeness = _completeness(family_id, domain, ortholog, dark, conflicting)
        basis = [name for name, present in (
            ("family_index", bool(family_id)), ("domains", bool(domain.get("domain_ids"))),
            ("orthology", bool(ortholog.get("orthogroup") or ortholog.get("orthology_reference_gene"))),
            ("annotation_text", bool(description or go_terms)),
        ) if present]
        # Issue #184: the subfamily layer is the finest supported functional /
        # evolutionary branch. The atlas does not take a curated subfamily map,
        # so the orthogroup is the honest finest available proxy (the same rule
        # family-to-copy applies), otherwise it stays unresolved.
        subfamily_id = str(ortholog.get("orthogroup", ""))
        subfamily_source = "orthogroup" if subfamily_id else "unresolved"
        # Copy-set preview per gene: the copy-set id is the homology layer id
        # when one exists (curated evolutionary family / orthogroup / domain
        # architecture with provenance), else the discovery bucket as a
        # provisional placeholder. The definitive family-level copy-set decision
        # belongs to family-to-copy.
        layer_id = (evolutionary_family_id or subfamily_id
                    or (str(domain.get("domain_architecture", "")) if str(domain.get("domain_ids", "")) else ""))
        if dark or not assigned:
            copy_set_id, copy_set_status = "", "unresolved"
        elif not family_id:
            copy_set_id, copy_set_status = "", "not_applicable"
        else:
            copy_set_id = layer_id or family_id
            copy_set_status = "validated" if layer_id else "provisional"
        rows.append({
            "gene_id": gene_id,
            "chrom": str(candidate.get("chrom", "") or family_row.get("chrom", "")),
            "start": str(candidate.get("start", "") or family_row.get("start", "")),
            "end": str(candidate.get("end", "") or family_row.get("end", "")),
            "protein_length": protein_length,
            "local_annotation_description": "\n".join(dict.fromkeys(
                str(record.get(key, "")).strip() for record in (candidate, annotation)
                for key in ("description", "function_description", "protein_family")
                if str(record.get(key, "")).strip())),
            "local_annotation_go_terms": ";".join(dict.fromkeys(
                str(record.get("go_terms", "")).strip() for record in (candidate, annotation)
                if str(record.get("go_terms", "")).strip())),
            "family_id": family_id, "family_label": family_label, "family_source": family_source,
            "family_members_in_interval": 0,
            "discovery_bucket_id": discovery_bucket_id, "discovery_bucket_label": discovery_bucket_label,
            "discovery_bucket_method": discovery_bucket_method,
            "evolutionary_family_id": evolutionary_family_id,
            "evolutionary_family_source": evolutionary_family_source,
            "evolutionary_family_status": evolutionary_family_status,
            "domain_ids": str(domain.get("domain_ids", "")),
            "domain_architecture": str(domain.get("domain_architecture", "")),
            "domain_source": str(domain.get("domain_source", "")),
            # Keep family-index protein identities and their search/competition
            # context separate from orthology. Exact source endpoints are
            # parsed mechanically from the declared subject fields only.
            "family_index_best_subject": str(family_row.get("best_subject", "")),
            "family_index_subjects": str(family_row.get("family_subjects", "")),
            "family_index_subject_endpoints": _family_index_reference_endpoints(family_row),
            "family_index_assignment_status": str(family_row.get("assignment_status", "")),
            "family_index_best_identity": str(family_row.get("best_identity", "")),
            "family_index_best_query_coverage": str(family_row.get("best_query_coverage", "")),
            "family_index_best_evalue": str(family_row.get("best_evalue", "")),
            "family_index_best_bitscore": str(family_row.get("best_bitscore", "")),
            "family_index_second_bitscore": str(family_row.get("second_bitscore", "")),
            "family_index_top2_bitscore_ratio": str(family_row.get("top2_bitscore_ratio", "")),
            "family_index_symbol_source_rank": str(family_row.get("symbol_source_rank", "")),
            "family_index_symbol_source_bitscore_ratio": str(family_row.get("symbol_source_bitscore_ratio", "")),
            "orthogroup": str(ortholog.get("orthogroup", "")),
            "orthology_relationship": str(ortholog.get("orthology_relationship", "")),
            "orthology_reference_gene": str(ortholog.get("orthology_reference_gene", "")),
            "orthology_providers": str(ortholog.get("orthology_providers", "")),
            "subfamily_id": subfamily_id, "subfamily_source": subfamily_source,
            "focal_paralog_count": 0, "tandem_expansion": False,
            "broad_function_bins": ";".join(bin_ids),
            "primary_broad_function": bin_ids[0] if bin_ids else "unknown",
            "broad_function_source": ";".join(sorted(set(bin_sources))) if bin_sources else "none",
            "broad_function_evidence": ";".join(bin_evidence),
            "conservation_class": _conservation_class(ortholog, family_id, str(domain.get("domain_ids", ""))),
            "evidence_completeness": completeness,
            "annotation_rescue_needed": completeness in {"dark_or_short_orf", "no_informative_annotation"},
            "copy_set_id": copy_set_id, "copy_set_status": copy_set_status,
            "classification_basis": ";".join(basis) or "none",
            "rank_contribution": "none",
            "notes": "conflicting family naming between family index and orthology reference" if conflicting else "",
            "gene_model_completeness": _gene_model_status(candidate),
        })

    atlas = pd.DataFrame(rows, columns=ATLAS_COLUMNS) if rows else pd.DataFrame(columns=ATLAS_COLUMNS)
    if not atlas.empty:
        counts = atlas[atlas.family_id.ne("")]["family_id"].value_counts()
        atlas["family_members_in_interval"] = [int(counts.get(str(value), 0)) for value in atlas["family_id"]]
        atlas["focal_paralog_count"] = [max(0, int(value) - 1) for value in atlas["family_members_in_interval"]]
        atlas["tandem_expansion"] = _tandem_flags(atlas, tandem_distance)
        atlas = atlas.sort_values("gene_id", kind="stable")

    family_summary = _family_summary(atlas, tandem_distance)
    domain_summary = _domain_summary(atlas)
    orthogroup_summary = _orthogroup_summary(atlas)
    function_summary = _function_summary(atlas, bins)
    unclassified = _unclassified(atlas)
    queries = _review_queries(family_summary, function_summary, unclassified)

    paths = {
        "interval_gene_atlas": out / "interval_gene_atlas.tsv",
        "family_summary": out / "family_summary.tsv",
        "domain_architecture_summary": out / "domain_architecture_summary.tsv",
        "orthogroup_summary": out / "orthogroup_summary.tsv",
        "broad_function_summary": out / "broad_function_summary.tsv",
        "unclassified_dark_candidates": out / "unclassified_dark_candidates.tsv",
        "family_review_queries": out / "family_review_queries.tsv",
    }
    for frame, path in ((atlas, paths["interval_gene_atlas"]), (family_summary, paths["family_summary"]),
                        (domain_summary, paths["domain_architecture_summary"]), (orthogroup_summary, paths["orthogroup_summary"]),
                        (function_summary, paths["broad_function_summary"]), (unclassified, paths["unclassified_dark_candidates"]),
                        (queries, paths["family_review_queries"])):
        frame.to_csv(path, sep="\t", index=False, lineterminator="\n")

    manifest = {
        "schema_version": SCHEMA_VERSION, "case_id": case_id,
        "phenotype_used": False,
        "phenotype_independence": "no phenotype input reaches this stage; overlaying one is a separate, later stage",
        "atlas_sha256": hashlib.sha256(paths["interval_gene_atlas"].read_bytes()).hexdigest(),
        "gene_count": int(len(atlas)),
        "classified_gene_count": int((atlas["evidence_completeness"] != "no_informative_annotation").sum()) if not atlas.empty else 0,
        "unclassified_count": int(len(unclassified)),
        "family_count": int(len(family_summary)),
        "distinct_broad_functions": int(len(function_summary)),
        "annotation_rescue_needed_count": int(atlas["annotation_rescue_needed"].astype(str).eq("True").sum()) if not atlas.empty else 0,
        "evidence_completeness_counts": {str(k): int(v) for k, v in atlas["evidence_completeness"].value_counts().sort_index().items()} if not atlas.empty else {},
        "gene_model_completeness_counts": {str(k): int(v) for k, v in atlas["gene_model_completeness"].value_counts().sort_index().items()} if not atlas.empty else {},
        "gene_model_completeness_vocabulary": list(GENE_MODEL_STATUSES),
        "conservation_counts": {str(k): int(v) for k, v in atlas["conservation_class"].value_counts().sort_index().items()} if not atlas.empty else {},
        "lineage_claim": {
            "claim_basis": (
                "conserved_cross_species: one-to-one orthology row; "
                "relative_enriched: orthology reference-gene match; "
                "unresolved: family / orthogroup / domain signal without a reference-gene match; "
                "no_detectable_homology: no family, orthology or domain signal; "
                "lineage_restricted_supported: reserved for an explicit multi-reference / outgroup negative-search "
                "input and never emitted by the current atlas inputs"
            ),
            "vocabulary": list(_CONSERVATION),
            "reference_species_searched": (str(family_index_meta.get("reference_species_searched") or "not_supplied")
                                           if family_index_meta else "not_supplied"),
            "reference_roles": [role for role, present in (
                ("family_index", bool(family_index_arg)),
                ("orthology", bool(orthology_arg)),
                ("domains", bool(domains_arg)),
            ) if present],
            "search_method": family_index_meta.get("search_method") if family_index_meta else "not_supplied",
            "thresholds": family_index_meta.get("thresholds") if family_index_meta else "not_supplied",
            "provider_coverage": {
                "orthology_providers": orthology_providers or None,
                "family_index_reference_subjects": int(family_in_interval["best_subject"].ne("").nunique())
                if family_index_arg and not family_in_interval.empty and "best_subject" in family_in_interval.columns else None,
            },
            "no_hit_count": {
                **no_hit_counts,
                "no_detectable_homology_genes": int((atlas["conservation_class"] == "no_detectable_homology").sum())
                if not atlas.empty else 0,
            },
            "outgroup_coverage": ("not_supplied; a multi-species outgroup negative search is required before any "
                                  "lineage-restricted claim, so lineage_restricted_supported stays reserved"),
            "family_index_search_provenance": family_index_meta,
        },
        "broad_function_bin_map": bin_meta,
        "family_semantics": {
            "discovery_bucket": "discovery_bucket_id / discovery_bucket_label / discovery_bucket_method: naming / broad functional search bucket (symbol_stem / curated_keyword). Recall-first exploration index; alone it never defines a homologous family",
            "evolutionary_family": "evolutionary_family_id / evolutionary_family_source / evolutionary_family_status: the homologous-family layer (supported / provisional / unresolved). Currently supplied by the family-index curated family map; orthogroup and domain layers stay separate columns",
            "orthogroup": "the orthogroup column is the provider-defined evolutionary group (orthogroup_id); recorded with relationship and providers, never merged into the discovery bucket",
            "subfamily": "subfamily_id / subfamily_source: the finest supported functional / evolutionary branch. The atlas does not take a curated subfamily map, so the orthogroup is the honest finest available proxy (same rule as family-to-copy), otherwise unresolved",
            "copy_set": "copy_set_id / copy_set_status are a per-gene preview (validated / provisional / not_applicable / unresolved). copy_set_id is the homology layer id when one exists (curated evolutionary family / orthogroup / domain architecture with provenance), else the discovery bucket as a provisional placeholder; the definitive family-level copy-set decision is made by family-to-copy from the family-level homology gate",
            "gate": "a discovery bucket alone yields copy_set_status provisional; family-to-copy requires orthogroup/HOG, curated family, or domain+provenance support before copy discrimination and never emits a single preferred_copy from a bucket-only family (issue #184)",
            "dark": "genes without any family assignment stay unresolved and are never dropped; a missing evolutionary family assignment is not negative evidence",
        },
        "inputs": {
            "candidates": file_fingerprint(runtime_path(candidates_arg)),
            "family_index": file_fingerprint(runtime_path(family_index_arg)) if family_index_arg else None,
            "domains": file_fingerprint(runtime_path(domains_arg)) if domains_arg else None,
            "orthology": file_fingerprint(runtime_path(orthology_arg)) if orthology_arg else None,
            "annotations": file_fingerprint(runtime_path(annotations_arg)) if annotations_arg else None,
        },
        "safety": [
            "family / domain / orthogroup / broad-function assignment is an exploration index, never causal evidence",
            "no gene is promoted for belonging to a well-known family, and nothing here writes a score",
            "unknown, conflicting, orphan and dark genes are kept as first-class buckets and never dropped",
            "a failed homology search is not evidence of lineage restriction; no_detectable_homology / unresolved never lower a candidate's position (issue #187)",
            "the discovery bucket and the evolutionary homologous family are separate layers; a naming bucket alone never defines a valid copy set (issue #184)",
            "multiple reference species reached through one annotation chain are recorded, not counted as separate evidence",
            "the atlas fingerprint is recorded so a later phenotype overlay cannot silently change it",
        ],
        "outputs": {name: str(path) for name, path in sorted(paths.items())},
    }
    manifest_path = out / "interval_atlas_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), "gene_count": int(len(atlas)), "family_count": int(len(family_summary)),
            "unclassified_count": int(len(unclassified)), "review_query_count": int(len(queries)),
            "atlas_sha256": manifest["atlas_sha256"],
            **{name: str(path) for name, path in paths.items()}}


# --------------------------------------------------------------------------- #
# summaries
# --------------------------------------------------------------------------- #

def _tandem_flags(atlas: pd.DataFrame, distance: int) -> list[bool]:
    flags: dict[str, bool] = {}
    for family_id, group in atlas[atlas.family_id.ne("")].groupby("family_id", sort=True):
        positions: dict[str, list[tuple[int, str]]] = {}
        for row in group.to_dict("records"):
            start = str(row["start"])
            positions.setdefault(str(row["chrom"]), []).append((int(start) if start.isdigit() else 0, str(row["gene_id"])))
        clustered: set[str] = set()
        for entries in positions.values():
            entries.sort()
            for (previous, left), (current, right) in zip(entries, entries[1:]):
                if current - previous <= distance:
                    clustered.update({left, right})
        for gene_id in group["gene_id"].astype(str):
            flags[gene_id] = gene_id in clustered
    return [bool(flags.get(str(value), False)) for value in atlas["gene_id"]]


def _blocks(group: pd.DataFrame, distance: int) -> int:
    blocks = 0
    positions: dict[str, list[int]] = {}
    for row in group.to_dict("records"):
        start = str(row["start"])
        positions.setdefault(str(row["chrom"]), []).append(int(start) if start.isdigit() else 0)
    for starts in positions.values():
        starts.sort()
        blocks += 1
        for previous, current in zip(starts, starts[1:]):
            if current - previous > distance:
                blocks += 1
    return blocks


def _family_summary(atlas: pd.DataFrame, distance: int) -> pd.DataFrame:
    if atlas.empty:
        return pd.DataFrame(columns=FAMILY_SUMMARY_COLUMNS)
    rows: list[dict[str, object]] = []
    for family_id, group in atlas[atlas.family_id.ne("")].groupby("family_id", sort=True):
        rows.append({
            "family_id": family_id,
            "family_label": ";".join(sorted({str(value) for value in group["family_label"] if str(value)})),
            "member_count": int(len(group)),
            "member_gene_ids": ";".join(sorted(group["gene_id"].astype(str))),
            "tandem_blocks": _blocks(group, distance),
            "orthology_one_to_one": int(group["orthology_relationship"].astype(str).str.contains("one.to.one", regex=True, case=False).sum()),
            "domain_supported": int(group["domain_ids"].astype(str).ne("").sum()),
            "notes": "member count is an exploration index; it is not a causal score",
        })
    return pd.DataFrame(rows, columns=FAMILY_SUMMARY_COLUMNS).sort_values(["member_count", "family_id"], ascending=[False, True], kind="stable")


def _domain_summary(atlas: pd.DataFrame) -> pd.DataFrame:
    if atlas.empty:
        return pd.DataFrame(columns=DOMAIN_SUMMARY_COLUMNS)
    rows: list[dict[str, object]] = []
    for architecture, group in atlas[atlas.domain_architecture.ne("")].groupby("domain_architecture", sort=True):
        rows.append({
            "domain_architecture": architecture, "member_count": int(len(group)),
            "member_gene_ids": ";".join(sorted(group["gene_id"].astype(str))),
            "domain_ids": ";".join(sorted({value for row in group["domain_ids"] for value in str(row).split(";") if value})),
            "notes": "domain order is retained; it discriminates otherwise-similar family members later",
        })
    return pd.DataFrame(rows, columns=DOMAIN_SUMMARY_COLUMNS).sort_values(["member_count", "domain_architecture"], ascending=[False, True], kind="stable")


def _orthogroup_summary(atlas: pd.DataFrame) -> pd.DataFrame:
    if atlas.empty:
        return pd.DataFrame(columns=ORTHOGROUP_SUMMARY_COLUMNS)
    rows: list[dict[str, object]] = []
    for orthogroup, group in atlas[atlas.orthogroup.ne("")].groupby("orthogroup", sort=True):
        rows.append({
            "orthogroup": orthogroup, "member_count": int(len(group)),
            "member_gene_ids": ";".join(sorted(group["gene_id"].astype(str))),
            "relationships": ";".join(sorted({value for row in group["orthology_relationship"] for value in str(row).split(";") if value})),
            "providers": ";".join(sorted({value for row in group["orthology_providers"] for value in str(row).split(";") if value})),
            "notes": "providers are recorded, not counted as independent evidence",
        })
    return pd.DataFrame(rows, columns=ORTHOGROUP_SUMMARY_COLUMNS).sort_values(["member_count", "orthogroup"], ascending=[False, True], kind="stable")


def _function_summary(atlas: pd.DataFrame, bins: pd.DataFrame) -> pd.DataFrame:
    if atlas.empty:
        return pd.DataFrame(columns=FUNCTION_SUMMARY_COLUMNS)
    labels = {str(row["bin_id"]): str(row["bin_label"]) for row in bins.to_dict("records")}
    members: dict[str, list[str]] = {}
    sources: dict[str, set[str]] = {}
    for row in atlas.to_dict("records"):
        gene_id = str(row["gene_id"])
        assigned = [value for value in str(row["broad_function_bins"]).split(";") if value] or ["unknown"]
        for bin_id in assigned:
            members.setdefault(bin_id, []).append(gene_id)
            sources.setdefault(bin_id, set()).update(value for value in str(row["broad_function_source"]).split(";") if value)
    rows = [{
        "bin_id": bin_id, "bin_label": labels.get(bin_id, "unclassified" if bin_id == "unknown" else bin_id),
        "member_count": len(gene_ids), "member_gene_ids": ";".join(sorted(gene_ids)),
        "assignment_sources": ";".join(sorted(sources.get(bin_id, set()))) or "none",
        "notes": "exploration index only; bin membership never changes a gene score",
    } for bin_id, gene_ids in sorted(members.items())]
    return pd.DataFrame(rows, columns=FUNCTION_SUMMARY_COLUMNS).sort_values(["member_count", "bin_id"], ascending=[False, True], kind="stable")


def _unclassified(atlas: pd.DataFrame) -> pd.DataFrame:
    if atlas.empty:
        return pd.DataFrame(columns=UNCLASSIFIED_COLUMNS)
    scope = atlas[atlas.evidence_completeness.isin({"dark_or_short_orf", "no_informative_annotation", "conflicting_annotation"})]
    rows = [{
        "gene_id": str(row["gene_id"]), "chrom": str(row["chrom"]), "start": str(row["start"]), "end": str(row["end"]),
        "protein_length": str(row["protein_length"]), "evidence_completeness": str(row["evidence_completeness"]),
        "why_unclassified": str(row["notes"]) or f"classification basis: {row['classification_basis']}",
        "annotation_rescue_needed": row["annotation_rescue_needed"],
        "required_next_action": ("run annotation-rescue or supply a closer reference proteome"
                                 if str(row["evidence_completeness"]) != "conflicting_annotation"
                                 else "reconcile the disagreeing family assignments before using either"),
    } for row in scope.to_dict("records")]
    return pd.DataFrame(rows, columns=UNCLASSIFIED_COLUMNS)


def _review_queries(family_summary: pd.DataFrame, function_summary: pd.DataFrame, unclassified: pd.DataFrame) -> pd.DataFrame:
    """Family/process-level research queries. Deliberately phenotype-free."""
    rows: list[dict[str, object]] = []
    for row in family_summary.to_dict("records"):
        rows.append({
            "query_id": f"q_fam_{len(rows):03d}", "scope": "family", "scope_id": str(row["family_id"]),
            "scope_label": str(row["family_label"]), "member_count": row["member_count"],
            "member_gene_ids": str(row["member_gene_ids"]),
            "suggested_query": f"{row['family_label'] or row['family_id']} family: developmental role and mutant phenotypes in model plants",
            "query_type": "family_literature",
            "rationale": f"{row['member_count']} interval gene(s) belong to this family",
            "phenotype_used": False,
        })
    for row in function_summary.to_dict("records"):
        if str(row["bin_id"]) == "unknown":
            continue
        rows.append({
            "query_id": f"q_fun_{len(rows):03d}", "scope": "broad_function", "scope_id": str(row["bin_id"]),
            "scope_label": str(row["bin_label"]), "member_count": row["member_count"],
            "member_gene_ids": str(row["member_gene_ids"]),
            "suggested_query": f"{row['bin_label']}: processes and regulators reported for this class in model plants",
            "query_type": "process_literature",
            "rationale": f"{row['member_count']} interval gene(s) fall in this broad functional bin",
            "phenotype_used": False,
        })
    if len(unclassified):
        rows.append({
            "query_id": f"q_unk_{len(rows):03d}", "scope": "unclassified", "scope_id": "unclassified",
            "scope_label": "unknown / conflicting / dark", "member_count": int(len(unclassified)),
            "member_gene_ids": ";".join(sorted(unclassified["gene_id"].astype(str))),
            "suggested_query": "ortholog phenotype, expression and knockout literature for conserved unknown proteins in this interval",
            "query_type": "unknown_family_probe",
            "rationale": "these genes carry no informative family assignment and must not be dropped for that reason",
            "phenotype_used": False,
        })
    return pd.DataFrame(rows, columns=REVIEW_QUERY_COLUMNS) if rows else pd.DataFrame(columns=REVIEW_QUERY_COLUMNS)
