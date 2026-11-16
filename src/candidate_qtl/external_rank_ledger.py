"""Coverage-aware external comparative-rank ledger.

Different OSS tools, external services, and network analyses produce their own
candidate rankings under different assumptions.  This ledger imports them into a
separate, auditable table and compares them against the native deterministic
``dry_score`` rank **without ever summing external scores into it**.

The consensus is a reciprocal-rank-fusion (RRF) screening view, not a causal
probability.  Rankers that share an information source are collapsed into one
``dependency_group`` so a single source cannot be double-counted, and an
uncovered candidate (never evaluated) is kept distinct from a covered candidate
that was evaluated and simply ranked low.
"""

from __future__ import annotations

import json
from collections import defaultdict

import pandas as pd

from . import schemas
from .core import file_fingerprint, runtime_path


LEDGER_VERSION = "1.0.0"
RRF_K = 60


def _numeric_rank(value: object) -> int | None:
    text = str(value).strip()
    if not text:
        return None
    try:
        rank = int(float(text))
    except ValueError:
        return None
    return rank if rank >= 1 else None


def external_rank_ledger(
    native_rank_arg: str,
    external_ranks_arg: str,
    out_arg: str,
    top_k: int = 20,
    rrf_k: int = RRF_K,
    discordance_spread: float = 0.5,
    full_hash: bool = False,
) -> dict[str, object]:
    """Build a coverage-aware external-rank ledger beside the native rank."""
    native_path = runtime_path(native_rank_arg)
    external_path = runtime_path(external_ranks_arg)
    out = runtime_path(out_arg)
    out.mkdir(parents=True, exist_ok=True)
    if top_k < 1:
        raise ValueError("top_k must be at least 1")
    if rrf_k < 1:
        raise ValueError("rrf_k must be at least 1")

    native = pd.read_csv(native_path, sep="\t", dtype=str).fillna("")
    if "gene_id" not in native.columns:
        raise ValueError("native rank table requires a gene_id column")
    if native["gene_id"].duplicated().any():
        raise ValueError("native rank table gene_id must be unique")
    if "screening_rank" in native.columns:
        native_rank = {str(row["gene_id"]): _numeric_rank(row["screening_rank"]) for _, row in native.iterrows()}
    else:
        # Fall back to dry_score order when an explicit screening_rank is absent.
        scored = native.copy()
        scored["_score"] = pd.to_numeric(scored.get("dry_score", 0), errors="coerce").fillna(0.0)
        scored = scored.sort_values(["_score", "gene_id"], ascending=[False, True], kind="stable").reset_index(drop=True)
        native_rank = {str(row["gene_id"]): index + 1 for index, row in scored.iterrows()}
    gene_universe = [str(g) for g in native["gene_id"]]

    external = pd.read_csv(external_path, sep="\t", dtype=str).fillna("")
    schemas.validate(external, schemas.EXTERNAL_RANK_INPUT_SCHEMA)
    for column in schemas.EXTERNAL_RANK_INPUT_COLUMNS:
        if column not in external.columns:
            external[column] = ""
    if external.duplicated(subset=["ranker_id", "gene_id"]).any():
        raise ValueError("external ranks must contain each (ranker_id, gene_id) pair at most once")

    # A dependency_group defaults to the declared source, else the ranker_id, so
    # rankers built on the same data are fused only once.
    def group_of(row: pd.Series) -> str:
        return str(row["dependency_group"]).strip() or str(row["source"]).strip() or str(row["ranker_id"]).strip()

    ranker_meta: dict[str, dict[str, str]] = {}
    covered_by_gene: dict[str, list[tuple[str, str, int]]] = defaultdict(list)  # gene -> [(ranker_id, family, rank)]
    group_best: dict[tuple[str, str], int] = {}  # (group, gene) -> best rank
    coverage_counts: dict[str, dict[str, int]] = defaultdict(lambda: {"covered": 0})
    for _, row in external.iterrows():
        ranker_id = str(row["ranker_id"]).strip()
        family = str(row["ranker_family"]).strip()
        gene = str(row["gene_id"]).strip()
        ranker_meta.setdefault(ranker_id, {
            "ranker_id": ranker_id, "ranker_family": family, "source": str(row["source"]).strip(),
            "version": str(row["version"]).strip(), "reference_species": str(row["reference_species"]).strip(),
            "dependency_group": group_of(row),
        })
        rank = _numeric_rank(row["rank"])
        covered = str(row["coverage"]).strip().casefold() == "covered" and rank is not None
        if not covered:
            continue
        if gene not in native_rank:
            # External rankers may only score candidates from the native set.
            raise ValueError(f"external rank references unknown candidate gene_id {gene!r}")
        coverage_counts[ranker_id]["covered"] += 1
        covered_by_gene[gene].append((ranker_id, family, rank))
        key = (group_of(row), gene)
        if key not in group_best or rank < group_best[key]:
            group_best[key] = rank

    # Reciprocal rank fusion across dependency groups (each group once).
    rrf_by_gene: dict[str, float] = defaultdict(float)
    contribution_rows: list[dict[str, object]] = []
    for (group, gene), rank in sorted(group_best.items()):
        contribution = 1.0 / (rrf_k + rank)
        rrf_by_gene[gene] += contribution
        contribution_rows.append({"gene_id": gene, "dependency_group": group, "group_best_rank": rank,
                                   "rrf_contribution": round(contribution, 8)})

    covered_genes = [g for g in gene_universe if g in rrf_by_gene]
    consensus_order = sorted(covered_genes, key=lambda g: (-rrf_by_gene[g], g))
    consensus_rank = {gene: index + 1 for index, gene in enumerate(consensus_order)}

    rows: list[dict[str, object]] = []
    for gene in gene_universe:
        native_r = native_rank.get(gene)
        covering = covered_by_gene.get(gene, [])
        families = sorted({family for _, family, _ in covering})
        ranks = [rank for _, _, rank in covering]
        native_top = native_r is not None and native_r <= top_k
        consensus_r = consensus_rank.get(gene)
        consensus_top = consensus_r is not None and consensus_r <= top_k
        if not covering:
            classification = "uncovered"
        elif native_top and consensus_top:
            classification = "consensus_top"
        elif consensus_top and not native_top:
            classification = "external_rescue"
        elif native_top and not consensus_top:
            classification = "native_only"
        else:
            classification = "discordant"
        spread = ""
        note = ""
        if ranks:
            rank_min, rank_max = min(ranks), max(ranks)
            denominator = max(len(gene_universe) - 1, 1)
            spread_value = round((rank_max - rank_min) / denominator, 6)
            spread = spread_value
            if len(families) > 1 and spread_value >= discordance_spread:
                note = "cross_family_rank_disagreement"
        else:
            rank_min = rank_max = ""
        rows.append({
            "gene_id": gene, "native_screening_rank": "" if native_r is None else native_r, "native_top": native_top,
            "covered_ranker_count": len(covering), "covering_families": ";".join(families),
            "consensus_rrf": "" if gene not in rrf_by_gene else round(rrf_by_gene[gene], 8),
            "consensus_rank": "" if consensus_r is None else consensus_r, "consensus_top": consensus_top,
            "classification": classification, "external_rank_min": rank_min, "external_rank_max": rank_max,
            "external_rank_spread": spread, "notes": note,
        })

    ledger = pd.DataFrame(rows, columns=schemas.EXTERNAL_RANK_LEDGER_COLUMNS).sort_values(["classification", "gene_id"], kind="stable")
    schemas.validate(ledger.fillna(""), schemas.EXTERNAL_RANK_LEDGER_SCHEMA)
    ledger_path = out / "external_rank_ledger.tsv"
    ledger.to_csv(ledger_path, sep="\t", index=False, lineterminator="\n")

    total_candidates = len(gene_universe)
    coverage_rows = []
    for ranker_id, meta in sorted(ranker_meta.items()):
        covered = coverage_counts[ranker_id]["covered"]
        coverage_rows.append({**meta, "candidates_total": total_candidates, "candidates_covered": covered,
                              "coverage_fraction": round(covered / total_candidates, 6) if total_candidates else 0.0})
    coverage_frame = pd.DataFrame(coverage_rows, columns=["ranker_id", "ranker_family", "dependency_group", "source", "version",
                                                          "reference_species", "candidates_total", "candidates_covered", "coverage_fraction"])
    coverage_path = out / "coverage_report.tsv"
    coverage_frame.sort_values(["ranker_family", "ranker_id"], kind="stable").to_csv(coverage_path, sep="\t", index=False, lineterminator="\n")

    discordance = ledger[ledger["classification"].isin(["discordant", "external_rescue", "native_only"]) | ledger["notes"].ne("")].copy()
    discordance_path = out / "rank_discordance.tsv"
    discordance.to_csv(discordance_path, sep="\t", index=False, lineterminator="\n")
    contribution_path = out / "consensus_contributions.tsv"
    pd.DataFrame(contribution_rows, columns=["gene_id", "dependency_group", "group_best_rank", "rrf_contribution"]).sort_values(
        ["gene_id", "dependency_group"], kind="stable").to_csv(contribution_path, sep="\t", index=False, lineterminator="\n")

    class_counts = {str(key): int(value) for key, value in ledger["classification"].value_counts().sort_index().items()}
    manifest = {
        "schema_version": "1.0", "ledger_version": LEDGER_VERSION, "rrf_k": rrf_k, "top_k": top_k,
        "consensus_method": "reciprocal_rank_fusion_over_dependency_groups",
        "candidate_count": total_candidates, "external_ranker_count": len(ranker_meta),
        "dependency_group_count": len({meta["dependency_group"] for meta in ranker_meta.values()}),
        "class_counts": class_counts, "families": sorted({meta["ranker_family"] for meta in ranker_meta.values()}),
        "native_rank": file_fingerprint(native_path, full_hash), "external_ranks": file_fingerprint(external_path, full_hash),
        "fingerprint_mode": "strict_full_hash" if full_hash else "size_tiered",
        "scientific_limits": [
            "external raw scores are never summed into the native dry_score",
            "the RRF consensus is a screening view, not a causal probability",
            "rankers sharing a dependency_group are fused once to avoid double counting",
            "uncovered means not evaluated; it is distinct from covered-but-low",
        ],
    }
    manifest_path = out / "external_rank_ledger_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"ledger": str(ledger_path), "coverage_report": str(coverage_path), "rank_discordance": str(discordance_path),
            "consensus_contributions": str(contribution_path), "manifest": str(manifest_path),
            "class_counts": class_counts, "external_ranker_count": len(ranker_meta)}
