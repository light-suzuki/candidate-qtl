"""Build the local Route B provider stack, then run the reverse lookup (issue #216).

``interval-reverse-lookup`` consumes an interval atlas and an optional
reference-function table. Nothing built them. In the monitored audit the command
therefore completed for 6/6 cases while doing none of the intended work: the
atlas carried no family, orthogroup, domain or source-backed function, so the
route degraded to broad-function bins plus an "unknown" probe and still reported
completion.

This module closes that gap for **this workstation** rather than for a generic
deployment: it resolves the local resources that are actually configured here,
builds the provider tables from them, assembles the atlas, and runs B0.

    candidates (+ candidate proteins)
      -> local homology search against a configured protein database
      -> orthology/homology table  (gene_id -> reference gene, with provenance)
      -> interval-scoped family index when no reusable index was supplied
      -> locally cached, source-backed reference-function table
      -> interval-atlas
      -> interval-reverse-lookup (B0, phenotype-free)
      -> family/process hypotheses with source ids and relationship tiers
      -> family-to-copy or an explicit Wet handoff

Rules this module does not bend:

* **The interval atlas stays phenotype-independent.** ``interval-atlas`` never
  receives the phenotype, so the enumeration and classification of *every*
  interval gene -- including the dark and unknown ones -- is fixed before any
  phenotype filter can narrow attention. The reverse lookup may then use the
  phenotype: withholding it from B was protecting an independence claim between
  A and B that never held, because both traverse the same literature. What
  actually has to be protected is **procedural** independence: the two lanes must
  not read each other's intermediate output before the comparison, which the
  profile enforces (see ``guided_profiles._assert_lane_isolation``).
* **A local best hit is homology, not orthology.** The generated relationship is
  ``best_hit_homology``; nothing here claims one-to-one orthology, so a BLAST-only
  run can never reach the ``exact_ortholog`` tier.
* **The automatic family index is interval-scoped.** It makes the guided Route B
  connection usable from candidate proteins, but is never advertised as the
  reusable genome-wide index required by canonical benchmark production.
* **A missing provider is a recorded state, never negative evidence.** Every
  provider carries how it was resolved, and an unavailable one leaves the lane
  ``not_evaluated_for_reverse_lookup`` instead of degrading quietly into bins.
* **Nothing here writes a score.** Family membership, convergence and handoffs
  all carry ``rank_contribution=none``.
* **The atlas handoff is one canonical path.** ``<out-dir>/interval_gene_atlas.tsv``
  (root level, same convention as ``interval-atlas``) is what the next Route B
  stage consumes; ``atlas/`` is the interval-atlas working directory (issue #247).
"""

from __future__ import annotations

import json
import os
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from . import __version__
from .core import file_fingerprint, runtime_path, run_tool, sha256
from .discover import _auto_local_protein_db
from .dual_route import interval_reverse_lookup
from .family_candidate_bridge import dependency_cluster_id, orthogroup_dependency_token
from .family_index import family_index as build_family_index
from .interval_atlas import GENE_MODEL_STATUSES, interval_atlas
from .route_b_blast_evidence import run_local_blastp
from .route_b_identity import METADATA_COLUMNS, annotation_material, protein_query_ids, usable, write_identity_layer

SCHEMA_VERSION = "1.0"

# Resource -> environment variable consulted when no explicit argument is given.
# Each is optional; an unset one is reported, not guessed at.
RESOURCE_ENV = {
    "protein_db": "CANDIDATE_QTL_PROTEIN_DB",
    "reference_function": "CANDIDATE_QTL_REFERENCE_FUNCTION",
    "orthology": "CANDIDATE_QTL_ORTHOLOGY_TABLE",
    "domains": "CANDIDATE_QTL_DOMAIN_TABLE",
    "family_index": "CANDIDATE_QTL_FAMILY_INDEX",
    "annotations": "CANDIDATE_QTL_ANNOTATION_TABLE",
}

ORTHOLOGY_COLUMNS = ["gene_id", "reference_gene_id", "orthogroup", "relationship", "provider",
                     "percent_identity", "query_coverage", "evalue", "bitscore",
                     "raw_reference_gene_id", "max_hsp_identity", "covered_query_residues",
                     "identical_query_residues", "conflicting_query_residues", "identity_policy", *METADATA_COLUMNS]
# Keep the normalized Route-B reference table a lossless-enough projection of
# the reference corpus.  In particular, the three literature/source rosters
# are distinct evidence axes (and must not be collapsed into ``source_id``).
# The semantic stage uses these fields as provenance only; they never create a
# score or an extra evidence vote.
REFERENCE_FUNCTION_COLUMNS = [
    "reference_gene_id", "gene_id", "function_description", "pathway", "pathway_id",
    "mutant_phenotype", "source_id", "source_url", "reference_species",
    "literature_source_ids", "phenotype_source_ids", "verified_source_ids",
    "evidence_codes", "evidence_type", "provenance_tier", "provenance_detail",
    "source_database", "source_release", "retrieved_at", "corpus_version",
    "tissue_specificity", "developmental_stage", "protein_family", "go_terms", "keywords",
    "taxon_id", "gene_symbol", "family_id", "orthogroup_id",
    "mapping_state", "mapping_relationship", "function_transfer_scope",
    "protein_name", "domain_ids", "domain_description",
]
HANDOFF_COLUMNS = ["family_id", "scope", "relationship_tier", "member_count", "member_interval_genes",
                   "reference_genes", "source_ids", "next_action", "rank_contribution", "notes"]

# A local best hit is only usable as a homology link above these thresholds. They
# are recorded in the manifest so a reader can see what "homolog" meant here.
DEFAULT_MIN_IDENTITY = 30.0
DEFAULT_MIN_COVERAGE = 0.5
DEFAULT_MAX_EVALUE = 1e-5

READINESS = ("source_backed_reverse_lookup", "partial", "not_evaluated_for_reverse_lookup")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _resolve(name: str, explicit: str | None, *, must_exist: bool = True) -> dict[str, object]:
    """Resolve one provider input, recording *how* it was resolved."""
    if explicit:
        path = runtime_path(explicit)
        if must_exist and not path.exists():
            raise ValueError(f"{name} does not exist: {path}")
        return {"available": True, "path": str(path), "resolution": "explicit_argument"}
    env_name = RESOURCE_ENV[name]
    value = os.environ.get(env_name, "").strip()
    if value:
        path = runtime_path(value)
        if path.exists() or not must_exist:
            return {"available": True, "path": str(path), "resolution": f"environment:{env_name}"}
        return {"available": False, "path": str(path), "resolution": f"environment:{env_name}",
                "reason": "configured path does not exist"}
    return {"available": False, "path": None, "resolution": "not_configured",
            "reason": f"no --{name.replace('_', '-')} argument and {env_name} is unset"}


_MATERIAL_COLUMNS = {
    "orthology": ("reference_gene_id", "orthogroup", "orthology_reference_gene"),
    "family_index": ("best_subject", "best_subject_symbol", "best_subject_stem", "family_id"),
    "domains": ("domain_id", "domain_ids", "domain_architecture"),
    "annotations": ("description", "function_description", "protein_family", "go_terms"),
}
def _usable_material_value(value: object, field: str) -> bool:
    if field in {"description", "function_description", "protein_family", "domain_architecture"}:
        return annotation_material(value)
    return usable(value)


def _usable_material_by_gene(
    path: Path | None, gene_ids: set[str], provider: str,
) -> set[str]:
    """Return interval genes with usable, provider-specific material.

    Presence of a file or a row is not material: empty/off-interval rows must
    leave a gene eligible for the local reference fallback.  This is a
    preparation gate only and makes no biological claim about the material.
    """
    if path is None or not path.exists() or not gene_ids:
        return set()
    try:
        frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    except (OSError, pd.errors.ParserError, UnicodeDecodeError):
        return set()
    if "gene_id" not in frame.columns:
        return set()
    fields = _MATERIAL_COLUMNS.get(provider, tuple(column for column in frame.columns if column != "gene_id"))
    fields = tuple(column for column in fields if column in frame.columns)
    if not fields:
        return set()
    covered: set[str] = set()
    for record in frame.to_dict("records"):
        gene = str(record.get("gene_id", "")).strip()
        if gene in gene_ids and any(_usable_material_value(record.get(column), column) for column in fields):
            covered.add(gene)
    return covered


def _merge_orthology_table(supplied: Path, hits: list[dict[str, object]], out_dir: Path) -> Path:
    """Keep supplied rows and append all fallback hits in one lossless table."""
    existing = pd.read_csv(supplied, sep="\t", dtype=str).fillna("")
    retrieved = pd.DataFrame(hits, columns=ORTHOLOGY_COLUMNS).fillna("")
    columns = list(dict.fromkeys([*existing.columns, *retrieved.columns, *ORTHOLOGY_COLUMNS]))
    merged = pd.concat([
        existing.reindex(columns=columns, fill_value=""),
        retrieved.reindex(columns=columns, fill_value=""),
    ], ignore_index=True).drop_duplicates()
    path = out_dir / "route_b_orthology.tsv"
    merged.to_csv(path, sep="\t", index=False, lineterminator="\n")
    return path


def _snapshot_reference_query(
    proteins: Path | None, unresolved_genes: set[str], interval_gene_ids: set[str], out_dir: Path,
) -> tuple[Path | None, set[str], dict[str, object]]:
    """Freeze exact-ID unresolved queries without rewriting input records."""
    original_ids = protein_query_ids(proteins)
    query_ids: set[str] = set()
    query_path: Path | None = None
    records: list[bytes] = []
    if proteins is not None:
        query_path = out_dir / "route_b_reference_query.fa"
        if proteins.resolve() == query_path.resolve():
            raise ValueError("candidate protein FASTA must differ from the reference-query snapshot")
        for record in re.split(rb"(?<=[\r\n])(?=>)", proteins.read_bytes()):
            if not record.startswith(b">"):
                continue
            lines = record.splitlines()
            gene = lines[0][1:].decode("utf-8").strip().split()[0]
            if gene in unresolved_genes and any(line.strip() for line in lines[1:]):
                query_ids.add(gene)
                records.append(record)
        query_path.write_bytes(b"".join(records))
    meta: dict[str, object] = {
        "original": file_fingerprint(proteins, full_hash=True) if proteins else None,
        "query": file_fingerprint(query_path, full_hash=True) if query_path else None,
        "selection_rule": "exact protein IDs of interval genes lacking usable supplied material",
        "missing_sequence_is_not_negative_evidence": True,
    }
    for label, ids in {
        "original": original_ids,
        "query": query_ids,
        "skipped": original_ids - query_ids,
        "missing": unresolved_genes - query_ids,
        "out_of_interval": original_ids - interval_gene_ids,
    }.items():
        meta[f"{label}_gene_ids"] = sorted(ids)
        meta[f"{label}_gene_count"] = len(ids)
    return query_path, query_ids, meta


def _normalise_subject(subject: str) -> str:
    """UniProt-style ``sp|P12345|NAME_ARATH`` -> ``P12345``; anything else is kept."""
    parts = str(subject).split("|")
    if len(parts) >= 4 and parts[0].casefold() == "gi" and parts[2].casefold() in {"ref", "gb", "emb", "dbj"} and parts[3]:
        return parts[3]
    if len(parts) >= 2 and parts[0].casefold() in {"sp", "tr", "gi", "ref", "gb", "emb", "dbj"} and parts[1]:
        return parts[1]
    return str(subject)


def _run_blastp(proteins: Path, database: Path, out_dir: Path, threads: int,
                min_identity: float, min_coverage: float, max_evalue: float) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Measure identity and coverage on the same unique query residues."""
    return run_local_blastp(
        proteins, database, out_dir, threads, min_identity, min_coverage, max_evalue,
        run_tool=run_tool, find_binary=shutil.which, normalise_subject=_normalise_subject,
    )


def _reference_function_table(paths: list[Path], out_dir: Path) -> tuple[Path | None, dict[str, object]]:
    """Normalise locally cached reference records into the reverse-lookup contract.

    The Route C reference corpus (``retrieval-index-build``) and a hand-curated
    reference-function table are both accepted: the columns are a superset of one
    another, so one local resource can feed both routes.
    """
    if not paths:
        return None, {"available": False, "reason": "no reference-function resource resolved",
                      "source_backed_rows": 0}
    rows: list[dict[str, object]] = []
    inputs: dict[str, object] = {}
    for path in paths:
        if path.is_dir():
            from .uniprot_reference_material import read_reference_directory
            records, metadata = read_reference_directory(path)
            frame = pd.DataFrame(records, columns=REFERENCE_FUNCTION_COLUMNS).fillna("")
            inputs[str(path)] = metadata
        else:
            frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
            inputs[str(path)] = file_fingerprint(path)
        if not {"reference_gene_id", "gene_id"} & set(frame.columns):
            raise ValueError(f"reference-function resource {path} requires reference_gene_id and/or gene_id")
        for record in frame.to_dict("records"):
            rows.append({
                "reference_gene_id": str(record.get("reference_gene_id", "")).strip(),
                "gene_id": str(record.get("gene_id", "")).strip(),
                "function_description": str(record.get("function_description", "") or record.get("description", "")).strip(),
                "pathway": str(record.get("pathway", "")).strip(),
                "pathway_id": str(record.get("pathway_id", "")).strip(),
                "mutant_phenotype": str(record.get("mutant_phenotype", "")).strip(),
                "source_id": str(record.get("source_id", "")).strip(),
                "source_url": str(record.get("source_url", "")).strip(),
                "reference_species": str(record.get("reference_species", "")).strip(),
                **{
                    column: str(record.get(column, "")).strip()
                    for column in REFERENCE_FUNCTION_COLUMNS
                    if column not in {
                        "reference_gene_id", "gene_id", "function_description", "pathway", "pathway_id",
                        "mutant_phenotype", "source_id", "source_url", "reference_species",
                    }
                },
            })
    table = pd.DataFrame(rows, columns=REFERENCE_FUNCTION_COLUMNS).drop_duplicates()
    path = out_dir / "route_b_reference_function.tsv"
    table.to_csv(path, sep="\t", index=False, lineterminator="\n")
    source_fields = [column for column in (
        "source_id", "literature_source_ids", "phenotype_source_ids", "verified_source_ids",
    ) if column in table.columns]
    source_backed = int(table[source_fields].astype(str).apply(lambda row: any(value.strip() for value in row), axis=1).sum()) if source_fields and not table.empty else 0
    return path, {
        "available": True, "path": str(path), "row_count": int(len(table)),
        "source_backed_rows": source_backed,
        "unsourced_rows": int(len(table)) - source_backed,
        "inputs": inputs, "sha256": sha256(path),
        "rule": "a row without a source_id cannot make a function claim source-backed; it is kept but never counted as one",
    }


def _handoff(hypotheses_path: Path, out_dir: Path) -> tuple[Path, int]:
    """Verified families -> family-to-copy or an explicit Wet handoff (no score)."""
    frame = pd.read_csv(hypotheses_path, sep="\t", dtype=str).fillna("")
    rows: list[dict[str, object]] = []
    for record in frame.to_dict("records"):
        if str(record.get("verification_status", "")) != "source_backed":
            continue
        if str(record.get("scope", "")) not in {"family", "orthogroup"}:
            continue
        members = [value for value in str(record.get("member_interval_genes", "")).split(";") if value]
        rows.append({
            "family_id": str(record.get("gene_family", "")), "scope": str(record.get("scope", "")),
            "relationship_tier": str(record.get("relationship_tier", "")),
            "member_count": len(members), "member_interval_genes": ";".join(members),
            "reference_genes": str(record.get("reference_genes", "")),
            "source_ids": str(record.get("source_backing", "")),
            "next_action": "family-to-copy" if len(members) > 1 else "wet_handoff_single_interval_copy",
            "rank_contribution": "none",
            "notes": "a source-backed family reaching the interval is a handoff target, not causal evidence",
        })
    path = out_dir / "route_b_family_handoff.tsv"
    pd.DataFrame(rows, columns=HANDOFF_COLUMNS).to_csv(path, sep="\t", index=False, lineterminator="\n")
    return path, len(rows)


DEPENDENCY_CLUSTER_COLUMNS = ["dependency_cluster", "dependency_status", "root_source_id", "source_ids"]
DEPENDENCY_CLUSTER_FILENAME = "evidence_dependency_clusters.tsv"


def _dependency_clusters_table(orthology_path: Path | None, atlas_path: Path, out_dir: Path) -> tuple[Path, int]:
    """Derive the Route B evidence-dependency cluster table (issue #374).

    One cluster per orthogroup/HOG observed in the Route B orthology table and
    the frozen atlas, so rows derived from the same orthogroup/HOG share the
    same source/homology dependency. The table matches the
    ``family-candidate-bridge --dependency-clusters`` input schema
    (``dependency_cluster`` / ``dependency_status`` / ``root_source_id`` /
    ``source_ids``). Clusters label evidence dependency only: they never add
    score and ``rank_contribution`` stays ``none`` everywhere.

    When neither table carries an orthogroup/HOG value (for example the
    local-blastp-derived orthology table, whose orthogroup column is empty, and
    no supplied orthology table), no cluster row is produced: the file is still
    written header-only and the caller reports zero clusters, so bridge rows
    simply keep ``dependency_unknown`` -- nothing is fabricated.
    """
    clusters: dict[str, set[str]] = {}

    def register(orthogroup_cell: object, extra_source_ids: list[str]) -> None:
        for part in re.split(r"[;,\s]+", str(orthogroup_cell or "")):
            part = part.strip()
            if not part:
                continue
            token = orthogroup_dependency_token(part)
            members = clusters.setdefault(token, {token})
            members.update(value for value in extra_source_ids if value)

    if orthology_path is not None and orthology_path.exists():
        frame = pd.read_csv(orthology_path, sep="\t", dtype=str).fillna("")
        if "orthogroup" in frame.columns:
            source_column = "source_id" if "source_id" in frame.columns else None
            for record in frame.to_dict("records"):
                extras = [str(record[source_column]).strip()] if source_column else []
                register(record.get("orthogroup", ""), extras)
    atlas_frame = pd.read_csv(atlas_path, sep="\t", dtype=str).fillna("")
    if "orthogroup" in atlas_frame.columns:
        for record in atlas_frame.to_dict("records"):
            register(record.get("orthogroup", ""), [])

    rows: list[dict[str, str]] = []
    for token in sorted(clusters):
        rows.append({
            "dependency_cluster": dependency_cluster_id(token),
            "dependency_status": "root_source_known",
            "root_source_id": token,
            "source_ids": ";".join(sorted(clusters[token])),
        })
    path = out_dir / DEPENDENCY_CLUSTER_FILENAME
    pd.DataFrame(rows, columns=DEPENDENCY_CLUSTER_COLUMNS).to_csv(
        path, sep="\t", index=False, lineterminator="\n")
    return path, len(rows)


def route_b_prepare(
    candidates_arg: str,
    out_dir_arg: str,
    proteins_arg: str | None = None,
    protein_db_arg: str | None = None,
    reference_function_args: list[str] | None = None,
    orthology_arg: str | None = None,
    domains_arg: str | None = None,
    family_index_arg: str | None = None,
    annotations_arg: str | None = None,
    case_id: str = "",
    threads: int = 4,
    min_identity: float = DEFAULT_MIN_IDENTITY,
    min_coverage: float = DEFAULT_MIN_COVERAGE,
    max_evalue: float = DEFAULT_MAX_EVALUE,
    auto_local_db: bool = True,
    phenotype: str = "",
) -> dict[str, object]:
    """Assemble the local Route B provider stack and run the reverse lookup.

    The atlas is built without the phenotype either way; ``phenotype`` reaches
    only the reverse lookup, where it marks which entries of the fixed inventory
    look related. Enumeration first, filtering second.
    """
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    candidates_path = runtime_path(candidates_arg)
    if not candidates_path.exists():
        raise ValueError(f"candidate table does not exist: {candidates_path}")

    providers: dict[str, object] = {
        "orthology_supplied": _resolve("orthology", orthology_arg),
        "domains": _resolve("domains", domains_arg),
        "family_index": _resolve("family_index", family_index_arg),
        "annotations": _resolve("annotations", annotations_arg),
    }

    # --- reference functions -------------------------------------------------
    reference_paths: list[Path] = []
    for value in reference_function_args or []:
        path = runtime_path(value)
        if not path.exists():
            raise ValueError(f"reference-function resource does not exist: {path}")
        reference_paths.append(path)
    if not reference_paths:
        resolved = _resolve("reference_function", None)
        providers["reference_function_resolution"] = resolved
        if resolved["available"]:
            reference_paths.append(runtime_path(str(resolved["path"])))
    else:
        providers["reference_function_resolution"] = {"available": True, "resolution": "explicit_argument",
                                                      "path": [str(path) for path in reference_paths]}
    reference_function_path, reference_meta = _reference_function_table(reference_paths, out)
    providers["reference_function"] = reference_meta

    # --- orthology / homology ------------------------------------------------
    orthology_path: Path | None = None
    homology_rows: list[dict[str, object]] = []
    candidate_frame = pd.read_csv(candidates_path, sep="\t", dtype=str).fillna("")
    if "gene_id" not in candidate_frame.columns:
        raise ValueError("candidate table requires gene_id")
    candidate_gene_ids = candidate_frame["gene_id"].astype(str).str.strip().tolist()
    if any(not gene for gene in candidate_gene_ids) or len(candidate_gene_ids) != len(set(candidate_gene_ids)):
        raise ValueError("candidate table requires non-empty unique interval gene IDs")
    interval_gene_ids = set(candidate_gene_ids)
    material_by_gene: dict[str, set[str]] = {gene: set() for gene in interval_gene_ids}
    for provider_name in ("orthology", "family_index", "domains", "annotations"):
        provider_key = "orthology_supplied" if provider_name == "orthology" else provider_name
        provider_info = providers[provider_key]
        if provider_info["available"]:  # type: ignore[index]
            covered = _usable_material_by_gene(
                runtime_path(str(provider_info["path"])), interval_gene_ids, provider_name)
            for gene in covered:
                material_by_gene[gene].add(provider_name)
    # Candidate tables can carry inline functional annotation.  Treat only
    # informative text/GO fields as material; generic placeholders remain
    # unresolved and can trigger the reference fallback.
    candidate_fields = tuple(column for column in _MATERIAL_COLUMNS["annotations"] if column in candidate_frame.columns)
    for record in candidate_frame.to_dict("records"):
        gene = str(record.get("gene_id", "")).strip()
        if gene and any(_usable_material_value(record.get(column), column) for column in candidate_fields):
            material_by_gene[gene].add("candidate_annotations")
    unresolved_genes = {gene for gene, sources in material_by_gene.items() if not sources}
    proteins = runtime_path(proteins_arg) if proteins_arg else None
    if proteins is not None and not proteins.exists():
        raise ValueError(f"candidate protein FASTA does not exist: {proteins}")
    query_path, query_ids, query_meta = _snapshot_reference_query(
        proteins, unresolved_genes, interval_gene_ids, out)
    completed_query_ids: set[str] = set()
    supplied_orthology = providers["orthology_supplied"]["available"]  # type: ignore[index]
    supplied_orthology_path = (runtime_path(str(providers["orthology_supplied"]["path"]))
                               if supplied_orthology else None)  # type: ignore[index]
    if supplied_orthology:
        # Keep the supplied mapping as the baseline even when fallback search
        # is needed or unavailable; retrieval can append to it below.
        orthology_path = supplied_orthology_path
    if supplied_orthology:
        supplied_path = supplied_orthology_path
        # A complete, meaningful supplied table is authoritative and avoids
        # touching the local DB. Empty/off-interval tables still fall through.
        if not unresolved_genes:
            orthology_path = supplied_path
            providers["homology_search"] = {
                "available": False, "evaluation_state": "not_evaluated_supplied_orthology",
                "reason": "supplied interval reference material covers every gene; no local search was run",
                "material_coverage": {"gene_count": len(interval_gene_ids), "unresolved_gene_count": 0},
            }
    if unresolved_genes:
        database: Path | None = None
        if query_ids and protein_db_arg:
            database = runtime_path(protein_db_arg)
        elif query_ids and auto_local_db:
            database = _auto_local_protein_db()
        if not query_ids or database is None:
            providers["homology_search"] = {
                "available": False,
                "evaluation_state": "not_evaluated_no_protein_query" if not query_ids else "not_evaluated_provider_unavailable",
                "reason": ("no candidate protein FASTA supplied" if proteins is None
                           else "no protein sequences for materially unresolved interval genes" if not query_ids
                           else "no local protein database configured (--protein-db / CANDIDATE_QTL_PROTEIN_DB)"),
                "missing_provider_is_not_negative_evidence": True,
                "material_coverage": {"gene_count": len(interval_gene_ids),
                                      "unresolved_gene_count": len(unresolved_genes),
                                      "unresolved_genes": sorted(unresolved_genes)},
            }
        else:
            rows, meta = _run_blastp(query_path, database, out, threads, min_identity, min_coverage, max_evalue)
            if meta.get("available"):
                completed_query_ids = query_ids
            providers["homology_search"] = meta
            providers["homology_search"]["fallback_reason"] = (  # type: ignore[index]
                "unresolved_interval_gene_material")
            providers["homology_search"]["material_coverage"] = {  # type: ignore[index]
                "gene_count": len(interval_gene_ids),
                "unresolved_gene_count": len(unresolved_genes),
                "unresolved_genes": sorted(unresolved_genes),
            }
            homology_rows = rows
            if rows:
                if supplied_orthology:
                    orthology_path = _merge_orthology_table(
                        supplied_orthology_path, rows, out)  # type: ignore[arg-type]
                else:
                    orthology_path = out / "route_b_orthology.tsv"
                    pd.DataFrame(rows, columns=ORTHOLOGY_COLUMNS).to_csv(
                        orthology_path, sep="\t", index=False, lineterminator="\n")
                providers["homology_search"]["orthology_table"] = str(orthology_path)  # type: ignore[index]
                providers["homology_search"]["sha256"] = sha256(orthology_path)  # type: ignore[index]
                if not providers["family_index"]["available"]:  # type: ignore[index]
                    auto_dir = out / "auto_family_index"
                    generated = build_family_index(
                        blast_hits_arg=str(meta["hits_table"]), out_dir_arg=str(auto_dir),
                        genes_arg=str(candidates_path), reference_label=str(database),
                        reference_species="mixed_or_unspecified_protein_database",
                        min_identity=min_identity, min_query_coverage=min_coverage,
                        max_evalue=max_evalue, case_id=case_id,
                        index_scope="interval_candidates",
                    )
                    providers["family_index"] = {
                        "available": True,
                        "path": str(generated["family_index"]),
                        "manifest": str(generated["manifest"]),
                        "resolution": "generated_from_local_blastp",
                        "index_scope": "interval_candidates",
                        "reusable_as_genome_wide_index": False,
                        "reason": (
                            "no explicit family index was supplied; generated from this run's candidate proteins "
                            "so Route B can classify the interval without claiming genome-wide coverage"
                        ),
                    }
                    providers["homology_search"]["auto_family_index"] = {  # type: ignore[index]
                        "path": str(generated["family_index"]),
                        "manifest": str(generated["manifest"]),
                        "index_scope": "interval_candidates",
                    }
            elif supplied_orthology:
                # Preserve the supplied mapping even when the fallback search
                # produced no threshold-passing hits.
                orthology_path = supplied_orthology_path
    elif not supplied_orthology:
        providers["homology_search"] = {
            "available": False,
            "evaluation_state": "not_evaluated_sufficient_material",
            "reason": "supplied interval reference material covers every gene; no local search was run",
            "material_coverage": {"gene_count": len(interval_gene_ids), "unresolved_gene_count": 0},
        }
    query_meta["completed_query_gene_ids"] = sorted(completed_query_ids)
    query_meta["completed_query_gene_count"] = len(completed_query_ids)
    providers["homology_search"]["query_subset"] = query_meta  # type: ignore[index]

    # --- atlas + B0 ----------------------------------------------------------
    atlas_dir = out / "atlas"
    atlas = interval_atlas(
        candidates_arg=str(candidates_path), out_dir_arg=str(atlas_dir),
        family_index_arg=str(providers["family_index"]["path"]) if providers["family_index"]["available"] else None,  # type: ignore[index]
        domains_arg=str(providers["domains"]["path"]) if providers["domains"]["available"] else None,  # type: ignore[index]
        orthology_arg=str(orthology_path) if orthology_path else None,
        annotations_arg=str(providers["annotations"]["path"]) if providers["annotations"]["available"] else None,  # type: ignore[index]
        case_id=case_id,
    )
    # Canonical Route B handoff: the atlas the next Route B stage
    # consumes is exposed at <out-dir>/interval_gene_atlas.tsv -- the same
    # root-level convention interval-atlas uses on its own -- so a wrapper
    # following the documented path never has to infer the nested layout.
    # The full atlas working directory (all summaries) stays under atlas/.
    canonical_atlas = out / "interval_gene_atlas.tsv"
    shutil.copy2(Path(atlas["interval_gene_atlas"]), canonical_atlas)
    # Issue #374 T3: gene-model completeness aggregation for the manifest, from
    # the frozen canonical atlas (which was built without the phenotype).
    atlas_frame = pd.read_csv(canonical_atlas, sep="\t", dtype=str).fillna("")
    # Keep every interval gene in a separate phenotype-free identity inventory.
    # Do not mutate the already fingerprinted interval-atlas output in place.
    reference_records = (pd.read_csv(reference_function_path, sep="\t", dtype=str).fillna("").to_dict("records")
                         if reference_function_path else [])
    _, identity_meta = write_identity_layer(
        atlas_frame.to_dict("records"), homology_rows, providers["homology_search"],
        reference_records, out,
        query_ids=completed_query_ids,
    )
    providers["identity_layer"] = identity_meta
    gene_model_counts = {status: int((atlas_frame["gene_model_completeness"] == status).sum())
                         for status in GENE_MODEL_STATUSES}
    # Issue #374 review: evidence-dependency clusters for the family-candidate
    # bridge. One cluster per orthogroup/HOG seen in the Route B orthology
    # table and the frozen atlas; rows derived from the same orthogroup/HOG
    # share the same source/homology dependency. Clusters never add score.
    # Zero clusters (no orthogroup/HOG data anywhere) is a recorded state, and
    # the bridge rows keep dependency_unknown -- nothing is fabricated.
    dependency_path, dependency_count = _dependency_clusters_table(orthology_path, canonical_atlas, out)
    b0_dir = out / "route_b0"
    # The atlas above never saw the phenotype, so the inventory is already fixed;
    # the phenotype only marks entries within it.
    b0 = interval_reverse_lookup(
        atlas_arg=str(atlas["interval_gene_atlas"]), out_dir_arg=str(b0_dir),
        reference_function_arg=str(reference_function_path) if reference_function_path else None,
        phenotype=phenotype, case_id=case_id,
    )
    b0_status = json.loads(Path(b0["status"]).read_text(encoding="utf-8"))
    handoff_path, handoff_count = _handoff(Path(b0["family_process_hypotheses"]), out)

    completeness = str(b0_status.get("process_inference_completeness", "not_evaluated_for_reverse_lookup"))
    readiness = ("source_backed_reverse_lookup" if completeness == "complete"
                 else ("partial" if completeness == "partial" else "not_evaluated_for_reverse_lookup"))

    manifest = {
        "schema_version": SCHEMA_VERSION, "command": "route-b-prepare",
        "candidate_qtl_version": __version__, "case_id": case_id, "prepared_at_utc": _now(),
        "phenotype_used": bool(phenotype),
        "atlas_phenotype_independent": True,
        "gene_model_completeness": gene_model_counts,
        "gene_model_completeness_vocabulary": list(GENE_MODEL_STATUSES),
        "enumeration_before_filtering": (
            "interval-atlas enumerates and classifies every interval gene without the phenotype, so dark and "
            "unknown genes are on the list before any phenotype filter narrows attention; the phenotype reaches "
            "only the reverse lookup"
        ),
        "independence_claim": (
            "procedural, not informational: this lane never reads Route A's output, and Route A never reads this "
            "one, until the comparison. A and B traverse the same literature, so their agreement is corroboration "
            "over a shared evidence base -- never an independent confirmation"
        ),
        "providers": providers,
        "provider_coverage": b0_status.get("provider_coverage"),
        "identity_coverage": identity_meta,
        "reverse_lookup_status": b0_status.get("reverse_lookup_status"),
        "process_inference_completeness": completeness,
        "route_b_readiness": readiness,
        "readiness_vocabulary": list(READINESS),
        "readiness_rule": (
            "source_backed_reverse_lookup requires at least one family/orthogroup/domain provider AND at least one "
            "source-backed reference function. A run whose only signal is broad-function bins is "
            "not_evaluated_for_reverse_lookup and must never be reported as an equivalent Route B result"
        ),
        "unresolved_bucket": {
            "unclassified_interval_genes": b0_status.get("unclassified_genes"),
            "rule": "unknown / orphan / dark interval genes stay an explicit unresolved bucket and are never dropped",
        },
        "hypotheses": {
            "count": b0.get("hypothesis_count"),
            "by_scope": b0_status.get("hypotheses_by_scope"),
            "relationship_tiers": b0_status.get("relationship_tier_counts"),
        },
        "handoff": {
            "table": str(handoff_path), "source_backed_families": handoff_count,
            "rule": "a verified family is handed to family-to-copy or to a Wet decision; convergence and handoff add no score",
        },
        "inputs": {"candidates": file_fingerprint(candidates_path),
                   "proteins": query_meta["original"], "reference_query_proteins": query_meta["query"]},
        "outputs": {
            "interval_gene_atlas": str(canonical_atlas),
            "reference_query_proteins": str(query_path) if query_path else None,
            "identity_hits": identity_meta["identity_hits"]["path"],
            "identity_inventory": identity_meta["identity_inventory"]["path"],
            "interval_atlas_working_dir": str(atlas_dir),
            "interval_atlas_manifest": str(atlas["manifest"]),
            "route_b0_status": str(b0["status"]),
            "family_process_hypotheses": str(b0["family_process_hypotheses"]),
            "interval_reverse_lookup": str(b0["interval_reverse_lookup"]),
            "reference_function": str(reference_function_path) if reference_function_path else None,
            "orthology": str(orthology_path) if orthology_path else None,
            "dependency_clusters": str(dependency_path),
            "family_index": (str(providers["family_index"]["path"])
                             if providers["family_index"]["available"] else None),  # type: ignore[index]
            "family_index_manifest": (str(providers["family_index"].get("manifest") or "") or None),  # type: ignore[union-attr]
            "family_handoff": str(handoff_path),
            "manifest": str(out / "route_b_preparation_manifest.json"),
        },
        "dependency_clusters": {
            "table": str(dependency_path), "cluster_count": dependency_count,
            "consumed_by": "family-candidate-bridge --dependency-clusters (guided triple-route chains this automatically)",
            "derivation_rule": (
                "one dependency cluster per orthogroup/HOG observed in the Route B orthology table and the "
                "frozen atlas; rows derived from the same orthogroup/HOG share the same source/homology "
                "dependency. The cluster table labels dependency only: it never adds score and "
                "rank_contribution stays none everywhere"
            ),
            "no_data_rule": (
                "when neither the orthology table nor the atlas carries an orthogroup/HOG value (the "
                "local-blastp-derived orthology table has an empty orthogroup column and no orthology table "
                "was supplied), zero clusters are generated: the file is header-only, the bridge rows keep "
                "dependency_unknown, and nothing is fabricated"
            ),
        },
        "canonical_atlas_rule": (
            "the canonical Route B atlas handoff is <out-dir>/interval_gene_atlas.tsv, exposed by "
            "outputs.interval_gene_atlas and consumed as-is by the next Route B stage; the nested "
            "atlas/ directory is the interval-atlas working directory (summaries and review tables) "
            "and is not the handoff path. The canonical file is byte-identical to "
            "atlas/interval_gene_atlas.tsv from the same run"
        ),
        "frozen_output_fingerprints": {
            "interval_gene_atlas": file_fingerprint(canonical_atlas, full_hash=True),
            "family_process_hypotheses": file_fingerprint(Path(b0["family_process_hypotheses"]), full_hash=True),
        },
        "safety": [
            "a local best hit is recorded as best_hit_homology; no one-to-one orthology is claimed from it",
            "an auto-generated interval-scoped family index is never reported as a reusable genome-wide index",
            "a missing provider leaves the lane not_evaluated_for_reverse_lookup and is never negative evidence",
            "broad-function bins alone do not make Route B evaluated",
            "nothing here writes a score; every hypothesis and handoff carries rank_contribution=none",
            "dependency clusters label evidence dependency only; they never add score",
        ],
        "next_step": (
            "reverse literature: query the literature for each interval gene's reference homolog, then "
            "source-verify the citations before any claim counts"
        ),
    }
    manifest_path = out / "route_b_preparation_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "manifest": str(manifest_path), "route_b_readiness": readiness,
        "reverse_lookup_status": b0_status.get("reverse_lookup_status"),
        "process_inference_completeness": completeness,
        "interval_gene_atlas": str(canonical_atlas),
        "identity_hits": identity_meta["identity_hits"]["path"],
        "identity_inventory": identity_meta["identity_inventory"]["path"],
        "route_b0_status": str(b0["status"]),
        "interval_reverse_lookup": str(b0["interval_reverse_lookup"]),
        "orthology": str(orthology_path) if orthology_path else None,
        "reference_function": str(reference_function_path) if reference_function_path else None,
        "family_index": (str(providers["family_index"]["path"])
                         if providers["family_index"]["available"] else None),  # type: ignore[index]
        "family_index_manifest": (str(providers["family_index"].get("manifest") or "") or None),  # type: ignore[union-attr]
        "dependency_clusters": str(dependency_path) if dependency_count else None,
        "family_process_hypotheses": str(b0["family_process_hypotheses"]),
        "family_handoff": str(handoff_path), "source_backed_families": handoff_count,
        "hypothesis_count": b0.get("hypothesis_count"), "gene_count": atlas.get("gene_count"),
    }
