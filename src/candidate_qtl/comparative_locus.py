"""Comparative locus reconstruction for draft-reference QTL cases (issue #174).

This is an optional heavy-rescue lane. Expensive marker alignment, orthogroup
inference, synteny indexing and focal-sequence search are prepared outside the
case run and imported as versioned tables. The deterministic case-time step
reconstructs an expected locus without treating comparative expectation as
causal evidence or replacing a usable native interval.
"""

from __future__ import annotations

import json
from pathlib import Path
import statistics
import time
import tracemalloc

import pandas as pd

from .core import file_fingerprint, runtime_path

SCHEMA_VERSION = "1.0"
PROJECTION_COLUMNS = [
    "species", "left_marker_id", "right_marker_id", "seqid", "left_position", "right_position",
    "orientation", "marker_order", "projection_status", "hit_count_max", "resource_cache_version", "notes",
]
BLOCK_COLUMNS = [
    "species", "seqid", "gene_id", "orthogroup", "reference_order", "consensus_order",
    "normalized_order", "species_role", "phylogenetic_distance", "block_status",
]
CONSENSUS_COLUMNS = [
    "consensus_order", "orthogroup", "species_support_count", "eligible_species_count", "support_fraction",
    "copy_count_max", "expected_class", "order_median", "order_range", "reconstruction_confidence",
    "causal_confidence", "rank_contribution", "notes",
]
RECOVERY_COLUMNS = [
    "orthogroup", "expected_class", "focal_gene_id", "focal_source", "recovery_status", "source_id",
    "alignment_identity", "query_coverage", "candidate_class", "reconstruction_confidence",
    "causal_confidence", "rank_contribution", "notes",
]
UNIVERSE_COLUMNS = [
    "candidate_id", "gene_id", "orthogroup", "candidate_class", "focal_source", "recovery_status",
    "reconstruction_confidence", "causal_confidence", "rank_contribution", "native_coordinate_status",
    "source_id", "notes",
]
REVIEW_COLUMNS = ["reason", "target", "required_next_action", "allowed_choices"]


def _read(path_arg: str | None, required: set[str], label: str) -> pd.DataFrame:
    if not path_arg:
        return pd.DataFrame(columns=sorted(required))
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{label} table is missing columns: {', '.join(missing)}")
    return frame


def _frame(rows: list[dict[str, object]], columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame(rows).reindex(columns=columns) if rows else pd.DataFrame(columns=columns)


def _bool(value: object, default: bool = False) -> bool:
    text = str(value).strip().casefold()
    if not text:
        return default
    return text in {"1", "true", "yes", "available", "present"}


def _integer(value: object, default: int = 0) -> int:
    try:
        return int(float(str(value)))
    except ValueError:
        return default


def _fasta_ids(path_arg: str) -> list[str]:
    path = runtime_path(path_arg)
    ids = [line[1:].split()[0] for line in path.read_text(encoding="utf-8").splitlines() if line.startswith(">")]
    if len(ids) < 2:
        raise ValueError("marker FASTA requires at least two flanking marker sequences")
    if len(ids) != len(set(ids)):
        raise ValueError("marker FASTA IDs must be unique")
    return ids


def _project_species(
    species: str,
    hits: pd.DataFrame,
    marker_ids: list[str],
    cache_version: str,
) -> dict[str, object]:
    left_id, right_id = marker_ids[0], marker_ids[-1]
    scope = hits[hits["species"].astype(str) == species]
    left = scope[scope["marker_id"].astype(str) == left_id]
    right = scope[scope["marker_id"].astype(str) == right_id]
    base = {
        "species": species, "left_marker_id": left_id, "right_marker_id": right_id,
        "seqid": "", "left_position": "", "right_position": "", "orientation": "",
        "marker_order": "unknown", "projection_status": "insufficient_marker_support",
        "hit_count_max": 0, "resource_cache_version": cache_version, "notes": "",
    }
    if left.empty or right.empty:
        base["notes"] = "both flanking markers require a supplied projection hit"
        return base
    maximum = max(
        max(_integer(value, 1) for value in left.get("hit_count", pd.Series("1", index=left.index))),
        max(_integer(value, 1) for value in right.get("hit_count", pd.Series("1", index=right.index))),
    )
    base["hit_count_max"] = maximum
    if maximum != 1 or len(left) != 1 or len(right) != 1:
        base["projection_status"] = "multi_hit"
        base["notes"] = "non-unique marker placement is retained for review"
        return base
    lrow, rrow = left.iloc[0], right.iloc[0]
    lseq, rseq = str(lrow["seqid"]), str(rrow["seqid"])
    lpos, rpos = _integer(lrow["start"]), _integer(rrow["start"])
    orientation_values = {str(lrow.get("orientation", "+")), str(rrow.get("orientation", "+"))}
    base.update({
        "seqid": lseq if lseq == rseq else f"{lseq};{rseq}", "left_position": lpos,
        "right_position": rpos, "orientation": ";".join(sorted(orientation_values)),
    })
    if lseq != rseq:
        base["projection_status"] = "split_locus"
        base["notes"] = "flanking markers project to different reference sequences"
    elif len(orientation_values) != 1:
        base["projection_status"] = "orientation_discordant"
        base["notes"] = "flanking marker orientations disagree"
    else:
        base["marker_order"] = "forward" if lpos <= rpos else "reverse"
        base["projection_status"] = "collinear"
    return base


def _order_discordance(blocks: pd.DataFrame, consensus_order: dict[str, int]) -> dict[str, int]:
    discordance: dict[str, int] = {}
    for species, group in blocks.groupby("species", sort=True):
        observed = [str(value) for value in group.sort_values("normalized_order", kind="stable")["orthogroup"] if str(value)]
        values = [consensus_order[value] for value in observed if value in consensus_order]
        inversions = sum(1 for i, left in enumerate(values) for right in values[i + 1:] if left > right)
        discordance[str(species)] = inversions
    return discordance


def comparative_locus_reconstruct(
    markers_arg: str,
    reference_panel_arg: str,
    marker_projections_arg: str,
    reference_locus_genes_arg: str,
    out_dir_arg: str,
    focal_recovery_arg: str | None = None,
    native_genes_arg: str | None = None,
    absence_evidence_arg: str | None = None,
    min_reference_species: int = 2,
    core_support_fraction: float = 0.67,
    case_id: str = "",
    full_hash: bool = False,
) -> dict[str, object]:
    """Reconstruct an expected QTL candidate universe from comparative caches."""
    if min_reference_species < 2:
        raise ValueError("min_reference_species must be at least 2")
    if not 0 < core_support_fraction <= 1:
        raise ValueError("core_support_fraction must be in (0, 1]")
    started = time.perf_counter()
    tracemalloc.start()
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    marker_ids = _fasta_ids(markers_arg)
    panel = _read(reference_panel_arg, {"species", "role", "available", "cache_version"}, "reference panel")
    hits = _read(marker_projections_arg, {"marker_id", "species", "seqid", "start"}, "marker projection")
    locus_genes = _read(reference_locus_genes_arg, {"species", "gene_id", "orthogroup", "order_index"}, "reference locus gene")
    focal = _read(focal_recovery_arg, {"orthogroup", "focal_gene_id", "source_class"}, "focal recovery")
    native = _read(native_genes_arg, {"gene_id", "orthogroup"}, "native gene")
    absence = _read(absence_evidence_arg, {"orthogroup", "evidence_status"}, "absence evidence")
    if panel["species"].duplicated().any():
        raise ValueError("reference panel requires unique species rows")
    eligible_panel = panel[panel["available"].map(_bool)].copy()
    if eligible_panel.empty:
        raise ValueError("reference panel has no available precomputed resources")
    panel_by_species = {str(row["species"]): row for row in panel.to_dict("records")}

    projections = _frame([
        _project_species(
            str(row["species"]), hits, marker_ids, str(row.get("cache_version", "")),
        ) for row in eligible_panel.sort_values("species", kind="stable").to_dict("records")
    ], PROJECTION_COLUMNS)
    supported_species = sorted(projections.loc[projections["projection_status"].eq("collinear"), "species"].astype(str))
    projection_discordant_species = sorted(
        projections.loc[~projections["projection_status"].eq("collinear"), "species"].astype(str)
    )
    reviews: list[dict[str, str]] = []
    for row in projections[~projections["projection_status"].eq("collinear")].to_dict("records"):
        reviews.append({
            "reason": f"marker_projection_{row['projection_status']}", "target": str(row["species"]),
            "required_next_action": "review marker uniqueness/order or retain this species as discordant; do not force a locus",
            "allowed_choices": "replace marker cache | exclude species with reason | stop reconstruction",
        })
    if len(supported_species) < min_reference_species:
        reviews.append({
            "reason": "insufficient_reference_support", "target": case_id or "case",
            "required_next_action": f"obtain collinear projections in at least {min_reference_species} reference species",
            "allowed_choices": "add reference cache | keep unresolved | stop reconstruction",
        })

    block_rows: list[dict[str, object]] = []
    for species in supported_species:
        projection = projections[projections["species"].eq(species)].iloc[0]
        scope = locus_genes[locus_genes["species"].astype(str) == species].copy()
        if scope.empty:
            reviews.append({
                "reason": "reference_locus_gene_cache_missing", "target": species,
                "required_next_action": "populate the local gene-order/orthogroup cache for this projected block",
                "allowed_choices": "update cache | exclude species with reason",
            })
            continue
        scope["_order"] = pd.to_numeric(scope["order_index"], errors="coerce")
        if scope["_order"].isna().any():
            raise ValueError(f"reference locus gene order_index must be numeric for {species}")
        forward = str(projection["marker_order"]) == "forward"
        max_order = float(scope["_order"].max())
        for record in scope.sort_values("_order", kind="stable").to_dict("records"):
            normalized = float(record["_order"]) if forward else max_order - float(record["_order"]) + 1
            meta = panel_by_species[species]
            block_rows.append({
                "species": species, "seqid": str(projection["seqid"]), "gene_id": str(record["gene_id"]),
                "orthogroup": str(record["orthogroup"]), "reference_order": float(record["_order"]),
                "consensus_order": "", "normalized_order": normalized, "species_role": str(meta.get("role", "")),
                "phylogenetic_distance": str(meta.get("phylogenetic_distance", "")), "block_status": "included",
            })
    blocks = _frame(block_rows, BLOCK_COLUMNS)
    eligible_gene_species = sorted(set(blocks["species"].astype(str))) if not blocks.empty else []
    denominator = len(eligible_gene_species)

    consensus_rows: list[dict[str, object]] = []
    for orthogroup, group in blocks[blocks["orthogroup"].ne("")].groupby("orthogroup", sort=True):
        per_species = group.groupby("species").size()
        support = int(len(per_species))
        fraction = round(support / max(1, denominator), 6)
        copy_max = int(per_species.max())
        order_values = [float(value) for value in group["normalized_order"]]
        order_median = statistics.median(order_values)
        expected_class = "duplication" if copy_max > 1 else ("conserved_core" if fraction >= core_support_fraction else "lineage_variable")
        confidence = "high" if expected_class == "conserved_core" and support >= 3 else ("medium" if support >= 2 else "weak")
        consensus_rows.append({
            "consensus_order": "", "orthogroup": str(orthogroup), "species_support_count": support,
            "eligible_species_count": denominator, "support_fraction": fraction, "copy_count_max": copy_max,
            "expected_class": expected_class, "order_median": order_median,
            "order_range": f"{min(order_values):g}-{max(order_values):g}", "reconstruction_confidence": confidence,
            "causal_confidence": "not_assessed_here", "rank_contribution": "none",
            "notes": "comparative expectation defines candidate-universe membership, never causal priority",
        })
        if copy_max > 1:
            reviews.append({
                "reason": "reference_duplication_or_tandem_expansion", "target": str(orthogroup),
                "required_next_action": "retain copy multiplicity and inspect tandem/polyploid/synteny context before focal copy selection",
                "allowed_choices": "retain every copy | add synteny context | add ploidy-aware orthogroup cache",
            })
    unresolved = blocks[blocks["orthogroup"].eq("")] if not blocks.empty else blocks
    if not unresolved.empty:
        reviews.append({
            "reason": "unresolved_reference_family", "target": f"{len(unresolved)} reference gene(s)",
            "required_next_action": "resolve orthogroup/domain assignment or retain the unresolved family explicitly",
            "allowed_choices": "update orthogroup cache | retain unresolved",
        })
    consensus = _frame(consensus_rows, CONSENSUS_COLUMNS).sort_values(["order_median", "orthogroup"], kind="stable")
    if not consensus.empty:
        consensus["consensus_order"] = range(1, len(consensus) + 1)
    else:
        reviews.append({
            "reason": "no_resolved_orthogroup_consensus", "target": case_id or "case",
            "required_next_action": "resolve orthogroup/domain families in at least two projected reference blocks",
            "allowed_choices": "update orthogroup cache | retain unresolved | stop reconstruction",
        })
    consensus_order = {str(row["orthogroup"]): int(row["consensus_order"]) for row in consensus.to_dict("records")}
    if not blocks.empty:
        blocks["consensus_order"] = blocks["orthogroup"].map(consensus_order).fillna("")
    discordance = _order_discordance(blocks, consensus_order) if not blocks.empty else {}
    discordant_species = sorted(species for species, count in discordance.items() if count > 0)
    for species in discordant_species:
        reviews.append({
            "reason": "rearrangement_suspected", "target": species,
            "required_next_action": "review local synteny/inversion rather than forcing majority gene order",
            "allowed_choices": "retain discordant block | split block | exclude species with reason",
        })

    if consensus.empty or len(eligible_gene_species) < min_reference_species:
        locus_status = "insufficient_reference_support"
    elif discordant_species:
        locus_status = "rearrangement_suspected"
    elif projection_discordant_species:
        locus_status = "species_discordant"
    elif len(eligible_gene_species) >= 3:
        locus_status = "consensus_strong"
    else:
        locus_status = "consensus_partial"

    native_by_group: dict[str, list[dict[str, object]]] = {}
    for row in native.to_dict("records"):
        native_by_group.setdefault(str(row["orthogroup"]), []).append(row)
    focal_by_group: dict[str, list[dict[str, object]]] = {}
    for row in focal.to_dict("records"):
        focal_by_group.setdefault(str(row["orthogroup"]), []).append(row)
    absence_by_group = {str(row["orthogroup"]): row for row in absence.to_dict("records")}
    recovery_rows: list[dict[str, object]] = []
    universe_rows: list[dict[str, object]] = []
    accepted_sources = {"scaffold_recovered", "transcript_supported", "proteome_supported", "unplaced_recovered"}
    for expected in consensus.to_dict("records"):
        group_id = str(expected["orthogroup"])
        recoveries: list[dict[str, object]] = []
        for row in native_by_group.get(group_id, []):
            recoveries.append({**row, "focal_gene_id": str(row["gene_id"]), "source_class": "native_annotated", "source_id": str(row.get("source_id", "native_annotation"))})
        recoveries.extend(focal_by_group.get(group_id, []))
        if not recoveries:
            evidence_status = str(absence_by_group.get(group_id, {}).get("evidence_status", ""))
            possible_pav = evidence_status in {"callable_absent", "parental_presence_absence", "assembly_supported_absence"}
            status = "possible_pav" if possible_pav else "expected_missing"
            candidate_class = "presence_absence_uncertain" if possible_pav else "comparative_expected_missing"
            recoveries = [{"focal_gene_id": "", "source_class": candidate_class, "source_id": str(absence_by_group.get(group_id, {}).get("source_id", ""))}]
            reviews.append({
                "reason": status, "target": group_id,
                "required_next_action": "search focal assembly/unplaced scaffolds/transcriptome and inspect callable parental evidence before claiming gene loss",
                "allowed_choices": "add focal search hit | add callable absence evidence | retain unresolved",
            })
        if len(recoveries) > 1:
            reviews.append({
                "reason": "focal_duplication_or_ambiguous_recovery", "target": group_id,
                "required_next_action": "retain every recovered copy and resolve with locus-specific synteny/linkage evidence",
                "allowed_choices": "keep all copies | add placement evidence",
            })
        for row in recoveries:
            source_class = str(row.get("source_class", ""))
            if source_class not in accepted_sources | {"native_annotated", "comparative_expected_missing", "presence_absence_uncertain"}:
                source_class = "unplaced_recovered"
            if source_class == "native_annotated":
                recovery_status, candidate_class = "native_present", "native_annotated"
            elif source_class == "comparative_expected_missing":
                recovery_status, candidate_class = "expected_missing", source_class
            elif source_class == "presence_absence_uncertain":
                recovery_status, candidate_class = "possible_pav", source_class
            elif source_class == "transcript_supported":
                recovery_status, candidate_class = "transcript_only", source_class
            else:
                recovery_status, candidate_class = "scaffold_recovered", "scaffold_recovered"
            gene_id = str(row.get("focal_gene_id", ""))
            candidate_id = gene_id or f"expected:{group_id}"
            common = {
                "orthogroup": group_id, "expected_class": str(expected["expected_class"]),
                "focal_gene_id": gene_id, "focal_source": source_class, "recovery_status": recovery_status,
                "source_id": str(row.get("source_id", "")), "alignment_identity": str(row.get("identity", "")),
                "query_coverage": str(row.get("query_coverage", "")), "candidate_class": candidate_class,
                "reconstruction_confidence": str(expected["reconstruction_confidence"]),
                "causal_confidence": "not_assessed_here", "rank_contribution": "none",
                "notes": "membership/recovery evidence only; no causal score is assigned",
            }
            recovery_rows.append(common)
            universe_rows.append({
                "candidate_id": candidate_id, "gene_id": gene_id or candidate_id, "orthogroup": group_id,
                "candidate_class": candidate_class, "focal_source": source_class, "recovery_status": recovery_status,
                "reconstruction_confidence": str(expected["reconstruction_confidence"]),
                "causal_confidence": "not_assessed_here", "rank_contribution": "none",
                "native_coordinate_status": "native" if source_class == "native_annotated" else "not_native_or_unresolved",
                "source_id": str(row.get("source_id", "")),
                "notes": (
                    "transparent expected:<orthogroup> placeholder; no focal gene model was claimed"
                    if not gene_id else
                    "primary/native coordinates are never replaced by this heavy rescue lane"
                ),
            })

    recovery = _frame(recovery_rows, RECOVERY_COLUMNS).sort_values(["orthogroup", "focal_gene_id"], kind="stable")
    universe = _frame(universe_rows, UNIVERSE_COLUMNS).sort_values(["orthogroup", "candidate_id"], kind="stable")
    review = _frame(reviews, REVIEW_COLUMNS).drop_duplicates().sort_values(["reason", "target"], kind="stable")
    paths = {
        "comparative_marker_projection": out / "comparative_marker_projection.tsv",
        "reference_locus_blocks": out / "reference_locus_blocks.tsv",
        "consensus_gene_order": out / "consensus_gene_order.tsv",
        "expected_orthogroups": out / "expected_orthogroups.tsv",
        "focal_recovery_hits": out / "focal_recovery_hits.tsv",
        "reconstructed_candidate_universe": out / "reconstructed_candidate_universe.tsv",
        "comparative_reconstruction_review": out / "comparative_reconstruction_review.tsv",
    }
    for frame, path in (
        (projections, paths["comparative_marker_projection"]), (blocks, paths["reference_locus_blocks"]),
        (consensus, paths["consensus_gene_order"]), (consensus, paths["expected_orthogroups"]),
        (recovery, paths["focal_recovery_hits"]), (universe, paths["reconstructed_candidate_universe"]),
        (review, paths["comparative_reconstruction_review"]),
    ):
        frame.to_csv(path, sep="\t", index=False, lineterminator="\n")

    _, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    elapsed = round(time.perf_counter() - started, 6)
    input_args = {
        "markers": markers_arg, "reference_panel": reference_panel_arg, "marker_projections": marker_projections_arg,
        "reference_locus_genes": reference_locus_genes_arg, "focal_recovery": focal_recovery_arg,
        "native_genes": native_genes_arg, "absence_evidence": absence_evidence_arg,
    }
    manifest = {
        "schema_version": SCHEMA_VERSION, "case_id": case_id, "lane": "comparative_locus_reconstruction",
        "lane_class": "optional_heavy_rescue", "standard_lane_dependency": False,
        "marker_sequence_start": True, "native_chromosome_coordinate_required": False,
        "resource_strategy": "case-time consumes versioned precomputed marker/orthogroup/synteny/gene-order caches; no whole-genome OrthoFinder or alignment is launched",
        "reference_species_available": sorted(eligible_panel["species"].astype(str)),
        "reference_species_collinear": supported_species, "reference_species_with_gene_order": eligible_gene_species,
        "reference_species_projection_discordant": projection_discordant_species,
        "reference_cache_versions": {str(row["species"]): str(row["cache_version"]) for row in eligible_panel.to_dict("records")},
        "locus_reconstruction_status": locus_status, "gene_order_inversions_by_species": discordance,
        "expected_orthogroup_count": int(len(consensus)), "candidate_universe_count": int(len(universe)),
        "candidate_class_counts": {str(k): int(v) for k, v in universe["candidate_class"].value_counts().sort_index().items()} if not universe.empty else {},
        "review_required_count": int(len(review)), "min_reference_species": min_reference_species,
        "core_support_fraction": core_support_fraction,
        "confidence_contract": {
            "reconstruction_confidence": "confidence that a family belongs in the expected locus universe",
            "causal_confidence": "not assessed by this command",
            "rank_contribution": "always none",
        },
        "native_coordinate_contract": "a usable native interval remains primary and is never replaced; reconstructed candidates are a separate class",
        "absence_contract": "missing from the focal assembly is unresolved unless callable absence evidence supports possible_pav; it is never automatically called gene loss",
        "runtime": {"elapsed_seconds": elapsed, "peak_traced_memory_mb": round(peak_bytes / 1_048_576, 6), "measurement": "python_tracemalloc"},
        "inputs": {name: (file_fingerprint(runtime_path(value), full_hash=full_hash) if value else None) for name, value in input_args.items()},
        "outputs": {name: str(path) for name, path in paths.items()},
        "output_fingerprints": {name: file_fingerprint(path, full_hash=True) for name, path in paths.items()},
        "scientific_limits": [
            "cross-species majority does not define truth; discordant species and rearrangements remain review-required",
            "duplication, lineage variation and unresolved families are retained instead of collapsed",
            "comparative expectation changes candidate-universe membership only and never adds a causal score",
        ],
    }
    manifest_path = out / "comparative_reconstruction_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "manifest": str(manifest_path), "candidate_universe": str(paths["reconstructed_candidate_universe"]),
        "locus_status": locus_status, "expected_orthogroup_count": int(len(consensus)),
        "candidate_universe_count": int(len(universe)), "review_required_count": int(len(review)),
    }


BENCHMARK_COLUMNS = [
    "case_id", "degradation_type", "reference_species_count", "truth_orthogroup", "truth_absent_from_native",
    "causal_orthogroup_recovered", "missing_gene_recovered", "expected_orthogroup_count", "false_expected_count",
    "false_expected_inflation", "reconstructed_candidate_set_size", "expected_discordance", "discordance_detected",
    "runtime_seconds", "peak_traced_memory_mb", "manifest_path",
]


def _read_group_set(path_arg: object, required_column: str) -> set[str]:
    text = str(path_arg).strip()
    if not text:
        return set()
    frame = _read(text, {required_column}, "benchmark orthogroup")
    return {str(value) for value in frame[required_column] if str(value)}


def comparative_locus_benchmark(manifest_arg: str, out_dir_arg: str, full_hash: bool = False) -> dict[str, object]:
    """Score degraded-reference recovery without using truth during reconstruction."""
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    manifest_path = runtime_path(manifest_arg)
    cases_in = _read(manifest_arg, {"case_id", "degradation_type", "truth_orthogroup", "reconstructed_universe_tsv", "reconstruction_manifest_json"}, "comparative benchmark manifest")

    def resolve(value: object) -> Path:
        path = runtime_path(str(value))
        return path if path.is_absolute() else manifest_path.parent / path

    allowed_degradation = {
        "causal_annotation_deleted", "causal_contig_unplaced", "interval_models_deleted", "marker_only_coordinate_removed",
    }
    unknown = sorted(set(cases_in["degradation_type"].astype(str)) - allowed_degradation)
    if unknown:
        raise ValueError(f"unknown degradation_type(s): {', '.join(unknown)}")
    rows: list[dict[str, object]] = []
    fingerprint_paths: list[Path] = [manifest_path]
    for case in cases_in.sort_values("case_id", kind="stable").to_dict("records"):
        universe_path = resolve(case["reconstructed_universe_tsv"])
        reconstruction_path = resolve(case["reconstruction_manifest_json"])
        universe = _read(str(universe_path), {"orthogroup", "candidate_class", "candidate_id"}, "reconstructed universe")
        reconstruction = json.loads(reconstruction_path.read_text(encoding="utf-8"))
        truth = str(case["truth_orthogroup"])
        native_set = _read_group_set(resolve(case["native_orthogroups_tsv"]), "orthogroup") if str(case.get("native_orthogroups_tsv", "")).strip() else set()
        truth_set = _read_group_set(resolve(case["expected_truth_set_tsv"]), "orthogroup") if str(case.get("expected_truth_set_tsv", "")).strip() else set()
        truth_set = truth_set or {truth}
        reconstructed_set = {str(value) for value in universe["orthogroup"] if str(value)}
        truth_absent = truth not in native_set
        recovered = truth in reconstructed_set
        false_count = len(reconstructed_set - truth_set)
        expected_discordance = _bool(case.get("expected_discordance", ""))
        status = str(reconstruction.get("locus_reconstruction_status", ""))
        detected = status in {"species_discordant", "rearrangement_suspected", "insufficient_reference_support"}
        runtime = reconstruction.get("runtime", {})
        rows.append({
            "case_id": str(case["case_id"]), "degradation_type": str(case["degradation_type"]),
            "reference_species_count": len(reconstruction.get("reference_species_with_gene_order", [])),
            "truth_orthogroup": truth, "truth_absent_from_native": truth_absent,
            "causal_orthogroup_recovered": recovered, "missing_gene_recovered": recovered if truth_absent else "",
            "expected_orthogroup_count": len(reconstructed_set), "false_expected_count": false_count,
            "false_expected_inflation": round(false_count / max(1, len(reconstructed_set)), 6),
            "reconstructed_candidate_set_size": int(len(universe)), "expected_discordance": expected_discordance,
            "discordance_detected": detected if expected_discordance else "", "runtime_seconds": runtime.get("elapsed_seconds", ""),
            "peak_traced_memory_mb": runtime.get("peak_traced_memory_mb", ""), "manifest_path": str(reconstruction_path),
        })
        fingerprint_paths.extend([universe_path, reconstruction_path])
        for key in ("native_orthogroups_tsv", "expected_truth_set_tsv"):
            if str(case.get(key, "")).strip():
                fingerprint_paths.append(resolve(case[key]))
    cases = _frame(rows, BENCHMARK_COLUMNS)
    cases_path = out / "comparative_locus_benchmark_cases.tsv"
    cases.to_csv(cases_path, sep="\t", index=False, lineterminator="\n")
    missing_scope = cases[cases["truth_absent_from_native"].eq(True)]
    discordance_scope = cases[cases["expected_discordance"].eq(True)]
    by_refs: list[dict[str, object]] = []
    for count, group in cases.groupby("reference_species_count", sort=True):
        by_refs.append({
            "reference_species_count": int(count), "case_count": int(len(group)),
            "causal_orthogroup_recovery_rate": round(float(group["causal_orthogroup_recovered"].mean()), 6),
            "mean_false_expected_inflation": round(float(pd.to_numeric(group["false_expected_inflation"]).mean()), 6),
            "mean_candidate_set_size": round(float(pd.to_numeric(group["reconstructed_candidate_set_size"]).mean()), 3),
            "mean_runtime_seconds": round(float(pd.to_numeric(group["runtime_seconds"], errors="coerce").mean()), 6),
            "mean_peak_traced_memory_mb": round(float(pd.to_numeric(group["peak_traced_memory_mb"], errors="coerce").mean()), 6),
        })
    summary = pd.DataFrame(by_refs)
    summary_path = out / "comparative_locus_benchmark_by_reference_count.tsv"
    summary.to_csv(summary_path, sep="\t", index=False, lineterminator="\n")
    result = {
        "schema_version": SCHEMA_VERSION, "benchmark_type": "synthetic_or_degraded_reference_comparative_locus_recovery",
        "truth_separation": "truth orthogroups are read only by this evaluator after reconstruction artifacts already exist",
        "case_count": int(len(cases)),
        "causal_orthogroup_recovery_rate": round(float(cases["causal_orthogroup_recovered"].mean()), 6) if len(cases) else 0.0,
        "missing_gene_recovery_rate": round(float(missing_scope["missing_gene_recovered"].mean()), 6) if len(missing_scope) else None,
        "mean_false_expected_inflation": round(float(pd.to_numeric(cases["false_expected_inflation"]).mean()), 6) if len(cases) else 0.0,
        "mean_reconstructed_candidate_set_size": round(float(pd.to_numeric(cases["reconstructed_candidate_set_size"]).mean()), 3) if len(cases) else 0.0,
        "synteny_discordance_detection_rate": round(float(discordance_scope["discordance_detected"].mean()), 6) if len(discordance_scope) else None,
        "by_reference_species_count": by_refs,
        "inputs": [file_fingerprint(path, full_hash=full_hash) for path in dict.fromkeys(fingerprint_paths)],
        "outputs": {"cases": str(cases_path), "by_reference_count": str(summary_path)},
        "scientific_limits": [
            "a recovered causal orthogroup is candidate-universe recall, not causal-gene identification",
            "false expected inflation requires an explicit expected truth-set; otherwise the known causal orthogroup is the conservative sole truth",
            "runtime and memory describe deterministic case-time reconstruction, not offline resource preparation",
        ],
    }
    manifest_out = out / "comparative_locus_benchmark_manifest.json"
    manifest_out.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    report = out / "comparative_locus_benchmark_report.md"
    report.write_text("\n".join([
        "# Comparative locus degraded-reference benchmark", "",
        f"- Cases: {result['case_count']}",
        f"- Causal orthogroup recovery: {result['causal_orthogroup_recovery_rate']:.1%}",
        f"- Missing-gene recovery: {result['missing_gene_recovery_rate'] if result['missing_gene_recovery_rate'] is not None else 'not evaluated'}",
        f"- Mean false expected inflation: {result['mean_false_expected_inflation']:.1%}",
        f"- Mean candidate-set size: {result['mean_reconstructed_candidate_set_size']}", "",
        "Recovery measures candidate-universe recall only. Comparative reconstruction never establishes causality.", "",
    ]), encoding="utf-8")
    return {"manifest": str(manifest_out), "cases": str(cases_path), "summary": str(summary_path), "report": str(report), "case_count": int(len(cases))}
