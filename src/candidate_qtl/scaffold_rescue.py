"""Unplaced-scaffold rescue into the QTL candidate space (issue #173).

The candidate universe is normally the genes inside a physical interval on a
native chromosome. In minor crops and draft assemblies the causal gene - or part
of it - can sit on an unplaced scaffold, an unanchored contig, a haplotig, a
separated tandem-family scaffold, or repeat/gap-adjacent sequence. Restricting
the universe to the chromosome interval then discards the real candidate because
of assembly quality, not biology.

This optional route rescues scaffolds that are *plausibly linked* to the QTL locus
and keeps them as a clearly separate candidate class.

Deliberate boundaries:

* Being unplaced is not itself a reason to add a scaffold. Some positive
  placement evidence is required.
* Family or homology agreement alone never produces a high-confidence rescue: it
  says nothing about linkage to this locus.
* A single short marker hit is QC'd for repeats and multi-mapping before it counts,
  and a paired marker is preferred over a lone hit.
* Placement confidence is plausibility of belonging to the locus, not causal
  probability, and is kept in its own column.
* Primary-interval genes are never replaced or reordered by rescued ones.
* Placement ambiguity is preserved: several candidate scaffolds are kept rather
  than forced into one answer.
* No rescue hit is **not** evidence that the causal gene is absent from unplaced
  sequence.
"""

from __future__ import annotations

import json

import pandas as pd

from .core import file_fingerprint, runtime_path

SCHEMA_VERSION = "1.0"
RESCUE_CLASSES = ("primary_interval", "scaffold_rescue_high", "scaffold_rescue_medium",
                  "scaffold_rescue_weak", "unresolved_scaffold", "scaffold_rejected_as_unlinked")
# Support types that can independently argue placement, and the ones that cannot.
STRONG_SUPPORT = ("paired_marker", "synteny_bridge", "hic_or_linkage")
WEAK_SUPPORT = ("single_marker", "orthogroup_context")

CANDIDATE_COLUMNS = [
    "scaffold_id", "rescue_class", "placement_support_type", "marker_support", "synteny_support",
    "orthogroup_context_support", "hic_or_linkage_support", "conflicting_evidence",
    "candidate_gene_count", "placement_confidence", "causal_confidence", "review_status", "reason",
]
EVIDENCE_COLUMNS = ["scaffold_id", "support_type", "strength", "detail", "source", "qc_status", "notes"]
GENE_COLUMNS = ["gene_id", "scaffold_id", "rescue_class", "candidate_class", "placement_confidence",
                "rank_contribution", "notes"]
REVIEW_COLUMNS = ["reason", "target", "required_next_action", "allowed_choices"]


def _frame(rows: list[dict[str, object]], columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame(rows)[columns] if rows else pd.DataFrame(columns=columns)


def _read(path_arg: str | None, required: set[str], label: str) -> pd.DataFrame:
    if not path_arg:
        return pd.DataFrame(columns=sorted(required))
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{label} table is missing columns: {', '.join(missing)}")
    return frame


def _marker_support(hits: pd.DataFrame, scaffold_id: str, max_scaffolds: int) -> tuple[str, str, str, str]:
    """Return (support_type, detail, qc_status, conflict) for one scaffold's marker hits."""
    if hits.empty:
        return "", "", "no_marker_hit", ""
    scope = hits[hits.scaffold_id.astype(str) == scaffold_id]
    if scope.empty:
        return "", "", "no_marker_hit", ""
    # A marker hitting several scaffolds is a repeat, not a placement. Only the
    # markers that touch *this* scaffold can make it ambiguous.
    scaffolds_per_marker = hits.groupby(hits.marker_id.astype(str))["scaffold_id"].nunique().to_dict()
    markers = sorted({str(value) for value in scope["marker_id"]})
    multi = sorted(marker for marker in markers if int(scaffolds_per_marker.get(marker, 1)) > max_scaffolds)
    conflict = f"multi_mapping_markers={';'.join(multi)}" if multi else ""
    informative = [marker for marker in markers if marker not in multi]
    if not informative:
        return "", ";".join(markers), "rejected_multi_mapping", conflict
    # A pair only counts when both of its markers survived the repeat check.
    informative_scope = scope[scope.marker_id.astype(str).isin(informative)]
    pairs = informative_scope["pair_id"].astype(str) if "pair_id" in informative_scope.columns \
        else pd.Series("", index=informative_scope.index)
    complete_pairs = sorted({value for value in pairs
                             if value and informative_scope[pairs == value]["marker_id"].nunique() >= 2})
    if complete_pairs:
        return "paired_marker", f"pairs={';'.join(complete_pairs)}", "paired_hit", conflict
    return "single_marker", f"markers={';'.join(informative)}", "single_hit_only", conflict


def scaffold_rescue(
    scaffolds_arg: str,
    out_dir_arg: str,
    marker_hits_arg: str | None = None,
    synteny_arg: str | None = None,
    orthogroup_context_arg: str | None = None,
    placement_support_arg: str | None = None,
    scaffold_genes_arg: str | None = None,
    max_marker_scaffolds: int = 1,
    case_id: str = "",
) -> dict[str, object]:
    """Conservatively rescue unplaced scaffolds that are plausibly linked to the locus."""
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    scaffolds = _read(scaffolds_arg, {"scaffold_id"}, "scaffold")
    if scaffolds.empty:
        raise ValueError("scaffold table is empty")
    markers = _read(marker_hits_arg, {"scaffold_id", "marker_id"}, "marker hit")
    synteny = _read(synteny_arg, {"scaffold_id"}, "synteny bridge")
    orthogroup = _read(orthogroup_context_arg, {"scaffold_id"}, "orthogroup context")
    placement = _read(placement_support_arg, {"scaffold_id"}, "placement support")
    genes = _read(scaffold_genes_arg, {"gene_id", "scaffold_id"}, "scaffold gene")

    evidence_rows: list[dict[str, object]] = []
    candidate_rows: list[dict[str, object]] = []
    gene_rows: list[dict[str, object]] = []
    review_rows: list[dict[str, str]] = []

    for record in scaffolds.sort_values("scaffold_id", kind="stable").to_dict("records"):
        scaffold_id = str(record["scaffold_id"])
        conflicts: list[str] = []
        strong: list[str] = []
        weak: list[str] = []

        support_type, detail, qc_status, conflict = _marker_support(markers, scaffold_id, max_marker_scaffolds)
        if conflict:
            conflicts.append(conflict)
        if support_type:
            (strong if support_type in STRONG_SUPPORT else weak).append(support_type)
        if qc_status != "no_marker_hit":
            evidence_rows.append({
                "scaffold_id": scaffold_id, "support_type": support_type or "marker_rejected",
                "strength": "strong" if support_type in STRONG_SUPPORT else ("weak" if support_type else "none"),
                "detail": detail, "source": "marker_alignment", "qc_status": qc_status,
                "notes": "a paired marker is preferred; a marker hitting several scaffolds is a repeat, not a placement",
            })
        if qc_status == "rejected_multi_mapping":
            review_rows.append({"reason": "marker_multi_mapping", "target": scaffold_id,
                                "required_next_action": "repeat由来のmarker hitなのでplacementに使えない。unique markerかlonger flanking sequenceを用意する",
                                "allowed_choices": "unique markerを足す | このscaffoldをunresolvedのまま保持"})

        synteny_scope = synteny[synteny.scaffold_id.astype(str) == scaffold_id] if not synteny.empty else synteny
        if not synteny_scope.empty:
            blocks = sorted({str(row.get("block_id", "")) for row in synteny_scope.to_dict("records") if str(row.get("block_id", ""))})
            species = sorted({str(row.get("anchor_species", "")) for row in synteny_scope.to_dict("records") if str(row.get("anchor_species", ""))})
            strong.append("synteny_bridge")
            evidence_rows.append({
                "scaffold_id": scaffold_id, "support_type": "synteny_bridge", "strength": "strong",
                "detail": f"blocks={';'.join(blocks)}", "source": ";".join(species) or "synteny",
                "qc_status": "bridge_supported",
                "notes": "conserved block in a high-quality relative bridges the native gap",
            })

        placement_scope = placement[placement.scaffold_id.astype(str) == scaffold_id] if not placement.empty else placement
        if not placement_scope.empty:
            types = sorted({str(row.get("support_type", "")) for row in placement_scope.to_dict("records") if str(row.get("support_type", ""))})
            strong.append("hic_or_linkage")
            evidence_rows.append({
                "scaffold_id": scaffold_id, "support_type": "hic_or_linkage", "strength": "strong",
                "detail": ";".join(types), "source": ";".join(sorted({str(row.get("source", "")) for row in placement_scope.to_dict("records") if str(row.get("source", ""))})) or "assembly_metadata",
                "qc_status": "placement_metadata", "notes": "Hi-C / linkage / optical map / assembler placement",
            })

        orthogroup_scope = orthogroup[orthogroup.scaffold_id.astype(str) == scaffold_id] if not orthogroup.empty else orthogroup
        if not orthogroup_scope.empty:
            weak.append("orthogroup_context")
            evidence_rows.append({
                "scaffold_id": scaffold_id, "support_type": "orthogroup_context", "strength": "weak",
                "detail": ";".join(sorted({str(row.get("orthogroup", "")) for row in orthogroup_scope.to_dict("records") if str(row.get("orthogroup", ""))})),
                "source": "orthogroup", "qc_status": "context_only",
                # Gene content says what the scaffold contains, not where it belongs.
                "notes": "gene-content agreement alone never proves linkage to this locus",
            })

        gene_count = int((genes.scaffold_id.astype(str) == scaffold_id).sum()) if not genes.empty else 0
        distinct_strong = sorted(set(strong))
        if conflicts and not distinct_strong:
            rescue_class = "unresolved_scaffold"
            reason = "placement evidence is ambiguous or repeat-driven"
        elif len(distinct_strong) >= 2:
            rescue_class = "scaffold_rescue_high"
            reason = f"independent placement support: {', '.join(distinct_strong)}"
        elif len(distinct_strong) == 1:
            rescue_class = "scaffold_rescue_medium"
            reason = f"single strong placement support: {distinct_strong[0]}"
        elif weak:
            rescue_class = "scaffold_rescue_weak"
            reason = f"weak support only: {', '.join(sorted(set(weak)))}"
        else:
            # Not "no causal gene here" - only "not linked on the evidence supplied".
            rescue_class = "scaffold_rejected_as_unlinked"
            reason = "no placement evidence supplied for this scaffold"
        if rescue_class in {"scaffold_rescue_weak", "unresolved_scaffold"}:
            review_rows.append({"reason": f"scaffold_{rescue_class}", "target": scaffold_id,
                                "required_next_action": "unique marker / synteny bridge / Hi-C のいずれかを足すまで高信頼placementにしない",
                                "allowed_choices": "証拠を足す | weakのまま別classで保持 | 対象外にする"})

        confidence = {"scaffold_rescue_high": "high", "scaffold_rescue_medium": "medium",
                      "scaffold_rescue_weak": "weak"}.get(rescue_class, "unresolved")
        candidate_rows.append({
            "scaffold_id": scaffold_id, "rescue_class": rescue_class,
            "placement_support_type": ";".join(distinct_strong + sorted(set(weak))),
            "marker_support": support_type or "none", "synteny_support": "yes" if not synteny_scope.empty else "no",
            "orthogroup_context_support": "yes" if not orthogroup_scope.empty else "no",
            "hic_or_linkage_support": "yes" if not placement_scope.empty else "no",
            "conflicting_evidence": ";".join(conflicts), "candidate_gene_count": gene_count,
            "placement_confidence": confidence,
            # Kept separate on purpose: where a scaffold belongs says nothing about
            # whether a gene on it is causal.
            "causal_confidence": "not_assessed_here",
            "review_status": "review_required" if rescue_class in {"scaffold_rescue_weak", "unresolved_scaffold"} else "context",
            "reason": reason,
        })
        if gene_count and rescue_class not in {"scaffold_rejected_as_unlinked"}:
            for gene in genes[genes.scaffold_id.astype(str) == scaffold_id].to_dict("records"):
                gene_rows.append({
                    "gene_id": str(gene["gene_id"]), "scaffold_id": scaffold_id, "rescue_class": rescue_class,
                    "candidate_class": "scaffold_rescue", "placement_confidence": confidence,
                    "rank_contribution": "none",
                    "notes": "kept as a separate candidate class; primary-interval genes are never replaced",
                })

    if markers.empty and synteny.empty and placement.empty and orthogroup.empty:
        review_rows.append({"reason": "no_placement_evidence_supplied", "target": case_id or "case",
                            "required_next_action": "marker/flanking alignment、synteny bridge、Hi-C/linkage のいずれも無い。scaffoldは救済できない",
                            "allowed_choices": "証拠を供給 | primary intervalのみで進む"})

    candidates = _frame(candidate_rows, CANDIDATE_COLUMNS)
    evidence = _frame(evidence_rows, EVIDENCE_COLUMNS).sort_values(["scaffold_id", "support_type"], kind="stable")
    gene_candidates = _frame(gene_rows, GENE_COLUMNS).sort_values(["scaffold_id", "gene_id"], kind="stable")
    reviews = _frame(review_rows, REVIEW_COLUMNS).drop_duplicates().sort_values(["reason", "target"], kind="stable")

    paths = {"scaffold_rescue_candidates": out / "scaffold_rescue_candidates.tsv",
             "scaffold_placement_evidence": out / "scaffold_placement_evidence.tsv",
             "scaffold_gene_candidates": out / "scaffold_gene_candidates.tsv",
             "scaffold_rescue_review_required": out / "scaffold_rescue_review_required.tsv"}
    for frame, path in ((candidates, paths["scaffold_rescue_candidates"]), (evidence, paths["scaffold_placement_evidence"]),
                        (gene_candidates, paths["scaffold_gene_candidates"]), (reviews, paths["scaffold_rescue_review_required"])):
        frame.to_csv(path, sep="\t", index=False, lineterminator="\n")

    manifest = {
        "schema_version": SCHEMA_VERSION, "case_id": case_id,
        "scaffolds_examined": int(len(scaffolds)),
        "rescue_class_counts": {str(k): int(v) for k, v in candidates["rescue_class"].value_counts().sort_index().items()} if not candidates.empty else {},
        "rescued_gene_count": int(len(gene_candidates)),
        "max_marker_scaffolds": max_marker_scaffolds,
        "confidence_meaning": "placement_confidence is the plausibility that the scaffold belongs to this QTL locus; it is not a causal probability",
        "no_hit_interpretation": "a scaffold without placement evidence is 'not linked on the evidence supplied', never 'the causal gene is not on unplaced sequence'",
        "inputs": {name: (file_fingerprint(runtime_path(value)) if value else None) for name, value in (
            ("scaffolds", scaffolds_arg), ("marker_hits", marker_hits_arg), ("synteny", synteny_arg),
            ("orthogroup_context", orthogroup_context_arg), ("placement_support", placement_support_arg),
            ("scaffold_genes", scaffold_genes_arg))},
        "safety": [
            "being unplaced is not itself a reason to add a scaffold; positive placement evidence is required",
            "family or homology agreement alone never yields a high-confidence rescue",
            "a marker hitting several scaffolds is treated as a repeat and rejected for placement",
            "placement confidence and causal confidence are separate columns",
            "primary-interval genes are never replaced or reordered by rescued genes",
            "placement ambiguity is preserved instead of being collapsed into one scaffold",
            "rescued genes carry rank_contribution=none and form their own candidate class",
        ],
        "outputs": {name: str(path) for name, path in sorted(paths.items())},
    }
    manifest_path = out / "scaffold_rescue_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), "scaffolds_examined": int(len(scaffolds)),
            "rescued_gene_count": int(len(gene_candidates)),
            "rescue_class_counts": manifest["rescue_class_counts"],
            **{name: str(path) for name, path in paths.items()}}
