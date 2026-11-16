"""Public Route-C direct pipeline facade.

Direct search uses the shared immutable literature cassette boundary from
issue #662. Source verification and freeze behavior remain in the v0.22
compatibility implementation. Projection and its benchmark consumer use the
provenance-complete, order-invariant contract introduced for issue #663.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from . import route_c_direct_pipeline_legacy as _legacy
from .core import file_fingerprint, runtime_path, sha256
from .route_c_direct_pipeline_legacy import (
    BENCHMARK_CONTAMINATION_STATES,
    DIRECT_BENCHMARK_CASE_COLUMNS,
    DIRECT_PROPOSAL_COLUMNS,
    SCHEMA_VERSION,
    SEARCH_COLUMNS,
    SPECIES_SCOPE_COLUMNS,
    _literature_query,
    direct_rescue_eligible_count,
    direct_rescue_freeze,
)
from .route_c_direct_projection import (
    PROJECTION_EVIDENCE_COLUMNS,
    PROJECTION_EVIDENCE_CONTRACT,
    PROJECTION_SELECTION_POLICY,
    direct_rescue_project,
    validate_projection_evidence_frame,
)
from .route_c_direct_search_cassette import direct_literature_search
from .route_c_direct_rescue import (
    SPECIES_DISTANCE,
    EVIDENCE_DIRECTNESS,
    validate_direct_rescue_frame,
    validate_projection_frame,
)


def _manifest_output_path(
    manifest_path: Path,
    manifest: dict[str, object],
    name: str,
) -> Path:
    outputs = manifest.get("outputs", {})
    if not isinstance(outputs, dict):
        raise ValueError("Route C projection manifest outputs must be an object")
    raw = str(outputs.get(name, "")).strip()
    if not raw:
        raise ValueError(f"Route C projection manifest lacks outputs.{name}")
    path = Path(raw)
    if not path.is_absolute():
        path = manifest_path.parent / path
    return path


def _verify_projection_evidence_generation(
    manifest_path: Path,
    manifest: dict[str, object],
) -> None:
    """Verify the optional v1 projection-evidence ledger before benchmark use."""
    contract = str(manifest.get("projection_evidence_contract", "")).strip()
    if not contract:
        return
    if contract != PROJECTION_EVIDENCE_CONTRACT:
        raise ValueError(f"unsupported Route C projection evidence contract: {contract}")
    evidence_path = _manifest_output_path(manifest_path, manifest, "projection_evidence")
    expected = str(manifest.get("projection_evidence_sha256", "")).strip().removeprefix("sha256:")
    if not expected or not evidence_path.exists() or sha256(evidence_path) != expected:
        raise ValueError("Route C benchmark projection evidence hash mismatch")
    validate_projection_evidence_frame(
        pd.read_csv(evidence_path, sep="\t", dtype=str).fillna("")
    )


def direct_rescue_benchmark(
    manifest_arg: str,
    out_dir_arg: str,
    full_hash: bool = False,
) -> dict[str, object]:
    """Score source recovery and projection using the projection v1 lineage.

    The case manifest continues to bind the raw atlas artifact SHA-256. Projection
    rows bind the canonical, order-independent interval-universe fingerprint
    declared by their own producer manifest. Keeping these identities separate
    prevents a harmless atlas-row permutation from invalidating the scientific
    projection generation while still authenticating the supplied input file.
    """
    manifest_path = runtime_path(manifest_arg)
    cases = _legacy._require_exact(
        pd.read_csv(manifest_path, sep="\t", dtype=str),
        DIRECT_BENCHMARK_CASE_COLUMNS,
        "Route C direct benchmark manifest",
    )
    rows: list[dict[str, object]] = []
    for case in cases.to_dict("records"):
        contamination = str(case["contamination_status"]).strip()
        if contamination not in BENCHMARK_CONTAMINATION_STATES:
            raise ValueError(
                "Route C direct benchmark contamination_status must be one of "
                f"{', '.join(sorted(BENCHMARK_CONTAMINATION_STATES))}; got {contamination or '<empty>'}"
            )
        direct_path = runtime_path(case["direct_rescue_tsv"])
        freeze_manifest_path = runtime_path(case["direct_freeze_manifest"])
        freeze_manifest = json.loads(freeze_manifest_path.read_text(encoding="utf-8"))
        freeze_valid = str(freeze_manifest.get("frozen_sha256", "")).removeprefix("sha256:") == sha256(direct_path)
        if not freeze_valid:
            raise ValueError(
                f"Route C benchmark frozen artifact hash mismatch for case {case['case_id']}"
            )
        direct = validate_direct_rescue_frame(
            pd.read_csv(direct_path, sep="\t", dtype=str).fillna("")
        )
        eligible = direct[
            direct["direct_rescue_eligibility"].eq("eligible_direct_record")
        ]
        truth_reference = str(case["truth_reference_gene_id"]).strip().casefold()
        clean = contamination == "clean"
        source_evaluated = bool(truth_reference and clean)
        source_hit = bool(
            source_evaluated
            and any(
                truth_reference
                in _legacy._split_identifiers(
                    ";".join(
                        [
                            str(row.get("reported_gene_id", "")),
                            str(row.get("reported_gene_symbol", "")),
                            str(row.get("reference_sequence_id", "")),
                        ]
                    )
                )
                for row in eligible.to_dict("records")
            )
        )
        truth_focal = str(case["truth_focal_gene_id"]).strip()
        projection_evaluated = bool(
            truth_focal and str(case["projection_tsv"]).strip() and clean
        )
        projection_hit = False
        if projection_evaluated:
            projection_path = runtime_path(case["projection_tsv"])
            projection_manifest_path = runtime_path(case["projection_manifest"])
            projection_manifest = json.loads(
                projection_manifest_path.read_text(encoding="utf-8")
            )
            if not isinstance(projection_manifest, dict):
                raise ValueError("Route C benchmark projection manifest must be a JSON object")
            if str(projection_manifest.get("frozen_sha256", "")).removeprefix("sha256:") != sha256(direct_path):
                raise ValueError(
                    f"Route C benchmark projection lineage mismatch for case {case['case_id']}"
                )
            if str(projection_manifest.get("projection_sha256", "")).removeprefix("sha256:") != sha256(projection_path):
                raise ValueError(
                    f"Route C benchmark projection artifact hash mismatch for case {case['case_id']}"
                )
            _verify_projection_evidence_generation(
                projection_manifest_path, projection_manifest
            )
            expected_qtl_id = str(case["expected_qtl_id"]).strip()
            if (
                not expected_qtl_id
                or str(projection_manifest.get("qtl_id", "")) != expected_qtl_id
            ):
                raise ValueError(
                    f"Route C benchmark projection qtl_id mismatch for case {case['case_id']}"
                )
            expected_atlas_sha = (
                str(case["expected_atlas_sha256"])
                .strip()
                .removeprefix("sha256:")
            )
            inputs = projection_manifest.get("inputs", {})
            atlas_input = inputs.get("atlas", {}) if isinstance(inputs, dict) else {}
            manifest_atlas_sha = (
                str(atlas_input.get("sha256", ""))
                if isinstance(atlas_input, dict)
                else ""
            ).removeprefix("sha256:")
            if not expected_atlas_sha or manifest_atlas_sha != expected_atlas_sha:
                raise ValueError(
                    f"Route C benchmark projection atlas mismatch for case {case['case_id']}"
                )
            projection = validate_projection_frame(
                pd.read_csv(projection_path, sep="\t", dtype=str).fillna("")
            )
            if set(projection["qtl_id"].astype(str)) - {"", expected_qtl_id}:
                raise ValueError(
                    f"Route C benchmark projection rows contain another qtl_id for case {case['case_id']}"
                )
            expected_fingerprint = str(
                projection_manifest.get("interval_universe_fingerprint", "")
            ).strip()
            if not expected_fingerprint:
                # Backward compatibility for projection manifests created before
                # the canonical fingerprint field existed.
                expected_fingerprint = f"sha256:{expected_atlas_sha}"
            if set(projection["interval_universe_fingerprint"].astype(str)) != {
                expected_fingerprint
            }:
                raise ValueError(
                    f"Route C benchmark projection rows contain another atlas for case {case['case_id']}"
                )
            projection_hit = bool(
                (
                    projection["projection_state"].eq("inside_current_qtl")
                    & projection["focal_gene_id"].eq(truth_focal)
                ).any()
            )
        rows.append(
            {
                "case_id": str(case["case_id"]),
                "contamination_status": contamination,
                "source_evaluation_status": (
                    "clean"
                    if source_evaluated
                    else (contamination if not clean else "not_evaluated")
                ),
                "source_recovery_hit": source_hit if source_evaluated else "",
                "eligible_direct_record_count": len(eligible),
                "projection_evaluation_status": (
                    "clean"
                    if projection_evaluated
                    else (contamination if not clean else "not_evaluated")
                ),
                "projection_hit": projection_hit if projection_evaluated else "",
            }
        )
    scored = pd.DataFrame(rows)
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    cases_path = out / "route_c_direct_benchmark_cases.tsv"
    scored.to_csv(cases_path, sep="\t", index=False, lineterminator="\n")

    def metric(status_column: str, hit_column: str) -> dict[str, object]:
        evaluated = scored[status_column].eq("clean")
        denominator = int(evaluated.sum())
        hits = int(
            scored.loc[evaluated, hit_column].astype(str).eq("True").sum()
        )
        return {
            "evaluated": denominator,
            "hits": hits,
            "recall": hits / denominator if denominator else None,
            "not_evaluated": int(scored[status_column].eq("not_evaluated").sum()),
            "contaminated": int(scored[status_column].eq("contaminated").sum()),
        }

    metrics = {
        "source_recovery": metric(
            "source_evaluation_status", "source_recovery_hit"
        ),
        "projection": metric("projection_evaluation_status", "projection_hit"),
    }
    summary = pd.DataFrame(
        [{"stage": key, **value} for key, value in metrics.items()]
    )
    summary_path = out / "route_c_direct_benchmark_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False, lineterminator="\n")
    output_manifest = {
        "schema_version": SCHEMA_VERSION,
        "command": "route-c-direct-benchmark",
        "case_count": len(scored),
        "metrics": metrics,
        "denominator_policy": (
            "not_evaluated and contaminated cases are not misses; source recovery and projection have separate denominators"
        ),
        "projection_fingerprint_policy": (
            "the case manifest authenticates the raw atlas input SHA-256; projection rows authenticate the canonical interval-universe fingerprint declared by the projection producer"
        ),
        "scoring_order": "prediction/freeze -> truth open -> score",
        "inputs": {"manifest": file_fingerprint(manifest_path, full_hash)},
        "outputs": {"cases": str(cases_path), "summary": str(summary_path)},
    }
    output_manifest_path = out / "route_c_direct_benchmark_manifest.json"
    output_manifest["outputs"]["manifest"] = str(output_manifest_path)
    output_manifest_path.write_text(
        json.dumps(output_manifest, ensure_ascii=False, indent=2, sort_keys=True)
        + "\n",
        encoding="utf-8",
    )
    return {
        "cases": str(cases_path),
        "summary": str(summary_path),
        "manifest": str(output_manifest_path),
        "metrics": metrics,
    }

__all__ = [
    "BENCHMARK_CONTAMINATION_STATES",
    "DIRECT_BENCHMARK_CASE_COLUMNS",
    "DIRECT_PROPOSAL_COLUMNS",
    "PROJECTION_EVIDENCE_COLUMNS",
    "PROJECTION_EVIDENCE_CONTRACT",
    "PROJECTION_SELECTION_POLICY",
    "SCHEMA_VERSION",
    "SEARCH_COLUMNS",
    "SPECIES_SCOPE_COLUMNS",
    "_literature_query",
    "direct_literature_search",
    "direct_rescue_benchmark",
    "direct_rescue_eligible_count",
    "direct_rescue_freeze",
    "direct_rescue_project",
    "validate_projection_evidence_frame",
]


def render_direct_proposal_prompt() -> str:
    """Render the interval-blind research contract from the strict schema enums.

    The harness copy is generated text, tested byte-for-byte against this source.
    Short display labels are deliberately not accepted or silently normalized.
    """
    return (
        "You own Route C: retrieve directly reported phenotype-gene/locus/variant relationships.\n"
        "Use only phenotype.json focal_species and phenotype; do not use QTL intervals, A/B outputs, truth or scores.\n"
        "Search public literature and sequence databases. Start with focal-species phenotype and cellular/developmental "
        "synonyms; use source-grounded gene synonyms, then broaden to related/model plants where useful.\n"
        "Generic function, homology and your proposed mechanism are not direct phenotype observations. "
        "Keep reported gene, experimental species/cultivar/assembly and reference sequence distinct. "
        "Never invent citations, identifiers or accessions.\n\n"
        f"Return TSV with exactly {len(DIRECT_PROPOSAL_COLUMNS)} columns, in this order:\n"
        + "\t".join(DIRECT_PROPOSAL_COLUMNS) + "\n"
        "Use the supplied route_c_direct_proposals.tsv header. record_id must be unique; "
        "query_phenotype preserves the input verbatim. Leave unknown fields empty.\n"
        "Use exactly one full literal value for each enum (no abbreviated S1/D0, translated labels or aliases):\n"
        "species_distance: " + " | ".join(sorted(SPECIES_DISTANCE)) + "\n"
        "evidence_directness: " + " | ".join(sorted(EVIDENCE_DIRECTNESS)) + "\n"
        "Example enum pair only: species_distance=S1_same_species; "
        "evidence_directness=D0_causal_perturbation_or_rescue. Choose from actual source facts.\n"
        "source_id is the PMID/PMCID/DOI actually reached; notes record inspected location and limitations. "
        "Multiple sources may be separate records; they are not independent focal genes.\n"
        "Save queries, URLs, tool names, raw responses, search coverage and access failures. "
        "Follow the execution budget; absence after a completed search differs from provider failure. "
        "For no findings return the header only and completed_empty; report failures separately.\n"
        "Do not project to the interval, rank, declare verified support or edit Wet membership. "
        "Source/identity/claim binding, freeze, projection and member validation happen downstream.\n"
    )


def render_harness_native_direct_prompt() -> str:
    """Ordinary C research using the assigned agent's native search/fetch tools.

    This prepends research scope only; the schema and verification instructions
    remain shared with the historical generated benchmark prompt. Neither the
    benchmark template nor the C freeze/identity/projection policy is modified.
    """
    return (
        "Independent Route C research. Ask the open question: Which genes, classical loci or variants have been reported in association with this phenotype in the focal species?\n"
        "Use only focal_species and phenotype from phenotype.json in a fresh research context. Do not read Route A/B results, QTL intervals, gene/orthogroup shortlists, benchmark truth, or parent research history.\n"
        "Use search/fetch tools actually available in your assigned harness. The dedicated route-c-direct-search command, a species-scope TSV and the Europe PMC provider are optional, not prerequisites. If search is unavailable, report not_evaluated rather than completed_empty.\n"
        "Search freely from the phenotype, synonyms, historical locus/mutant names, and reported gene/variant associations. Start with the focal species, then broaden to related/model plants where useful. Do not pre-filter by gene identity, family, uniqueness, or any interval.\n"
        "Genes, classical loci and variants can all be findings; a gene ID or genomic position need not be known. Preserve multiple source records for the same reported entity with different record_id values. Do not suppress overlap with other routes or remove records that might fall outside a future QTL.\n"
        "Record observed trait-gene/locus/variant relations and their actual evidence directness. Do not transform generic function or homology into a reported causal fact. Preserve source IDs, read scope, queries, retrieval URLs and failures; do not invent citations or accessions.\n"
        "Research first, then serialize to the existing C proposal contract below. Use the existing C source checks and freeze; project to the QTL only after freeze.\n\n"
        + render_direct_proposal_prompt()
    )
