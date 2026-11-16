"""Laptop-safe staged orthology: execution profiles and out-of-interval paralogs.

Tier 1 keeps the query limited to QTL-interval candidate proteins while reusing
the whole focal proteome and a small reference panel as prebuilt search targets.
Searching the full focal proteome (not just the interval) is what lets us see
paralogs that live *outside* the QTL, which is exactly where a naive
interval-only reciprocal search would wrongly call a multi-copy gene a 1:1
ortholog. Execution profiles (accuracy64 / laptop16) only change thread count,
hit caps, and tree policy -- never the scientific labels or scoring.
"""

from __future__ import annotations

import json
import shlex
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from .core import runtime_path, safe_id, sha256
from .evidence_io import EVIDENCE_COLUMNS, best_hits as _best_hits, evidence_row as _row


# Defaults mirror templates/runtime_profiles.tsv; a --profiles-table overrides them.
DEFAULT_PROFILES: dict[str, dict[str, object]] = {
    "accuracy64": {
        "ram_target_gb": 64,
        "threads": 8,
        "reference_proteomes": "3-6",
        "max_hits_per_query": 50,
        "orthology_refinement": "Top+Rescue+ambiguous+selected_no_hit",
        "family_tree_policy": "ambiguous_top_families",
    },
    "laptop16": {
        "ram_target_gb": 16,
        "threads": 4,
        "reference_proteomes": "2-3",
        "max_hits_per_query": 25,
        "orthology_refinement": "Top+Rescue+ambiguous",
        "family_tree_policy": "one_family_at_a_time",
    },
}


def load_runtime_profile(profile: str, table_arg: str | None = None) -> dict[str, object]:
    """Resolve an execution profile from the built-in defaults or a TSV table."""
    if table_arg:
        frame = pd.read_csv(runtime_path(table_arg), sep="\t", dtype=str).fillna("")
        if "profile" not in frame.columns:
            raise ValueError("runtime profile table must contain a profile column")
        profiles = {str(row["profile"]): {k: row[k] for k in frame.columns if k != "profile"} for _, row in frame.iterrows()}
    else:
        profiles = DEFAULT_PROFILES
    if profile not in profiles:
        raise ValueError(f"unknown profile {profile!r}; available: {', '.join(sorted(profiles))}")
    resolved = {"profile": profile}
    for key, value in profiles[profile].items():
        resolved[key] = value
    return resolved


def _candidate_ids(candidates_arg: str) -> set[str]:
    frame = pd.read_csv(runtime_path(candidates_arg), sep="\t", dtype=str).fillna("")
    if "gene_id" not in frame.columns:
        raise ValueError("candidate table must contain gene_id")
    return {str(value) for value in frame["gene_id"] if str(value).strip()}


def orthology_tier1_plan(
    candidates_arg: str,
    protein_arg: str,
    focal_proteome_arg: str,
    reference_proteome_args: list[str],
    out_dir_arg: str,
    profile: str = "laptop16",
    profiles_table: str | None = None,
    aligner: str = "diamond",
    threads: int | None = None,
    max_hits: int | None = None,
) -> dict[str, object]:
    """Generate a Tier-1 search plan: candidate query vs focal self + reference panel."""
    if aligner not in {"diamond", "blastp"}:
        raise ValueError("aligner must be diamond or blastp")
    prof = load_runtime_profile(profile, profiles_table)
    n_threads = int(threads if threads is not None else prof["threads"])
    n_hits = int(max_hits if max_hits is not None else prof["max_hits_per_query"])

    candidate_ids = _candidate_ids(candidates_arg)
    with runtime_path(protein_arg).open(encoding="utf-8") as handle:
        query_proteins = [record for record in SeqIO.parse(handle, "fasta") if record.id in candidate_ids]
    if not query_proteins:
        raise ValueError("no candidate proteins matched the candidate table gene_id values")

    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    query_path = out_dir / "tier1_query_candidates.fa"
    SeqIO.write(query_proteins, str(query_path), "fasta")

    focal = runtime_path(focal_proteome_arg)
    lines = ["#!/usr/bin/env bash", "set -euo pipefail", f"OUT={shlex.quote(str(out_dir))}", 'mkdir -p "$OUT/db" "$OUT/results"']
    searches: list[dict[str, str]] = []

    def add_search(label: str, subject: Path, is_focal: bool) -> None:
        result = out_dir / "results" / f"{label}.tsv"
        if aligner == "diamond":
            lines.append(f"diamond makedb --in {shlex.quote(str(subject))} --db \"$OUT/db/{label}\" >/dev/null 2>&1")
            lines.append(
                f"diamond blastp --query {shlex.quote(str(query_path))} --db \"$OUT/db/{label}\" "
                f"--max-target-seqs {n_hits} --threads {n_threads} --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore "
                f"--out {shlex.quote(str(result))}"
            )
        else:
            lines.append(f"makeblastdb -in {shlex.quote(str(subject))} -dbtype prot -parse_seqids -out \"$OUT/db/{label}\" >/dev/null")
            lines.append(
                f"blastp -query {shlex.quote(str(query_path))} -db \"$OUT/db/{label}\" -max_target_seqs {n_hits} "
                f"-num_threads {n_threads} -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore' "
                f"-out {shlex.quote(str(result))}"
            )
        searches.append({"label": label, "subject": str(subject), "result": str(result), "role": "focal_self_search" if is_focal else "reference_panel"})

    add_search("focal_self", focal, True)
    for specification in reference_proteome_args:
        if "=" not in specification:
            raise ValueError("--reference-proteome must use LABEL=FASTA")
        label, path_text = specification.split("=", 1)
        add_search(safe_id(label), runtime_path(path_text), False)

    script_path = out_dir / "run_tier1_orthology.sh"
    script_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    script_path.chmod(0o755)

    manifest = {
        "tier": 1,
        "profile": prof,
        "threads": n_threads,
        "max_hits_per_query": n_hits,
        "query_candidate_count": len(query_proteins),
        "query_is_candidates_only": True,
        "runs_full_proteome_all_vs_all": False,
        "focal_self_search": True,
        "family_tree_policy": prof.get("family_tree_policy"),
        "script": str(script_path),
        "searches": searches,
        "principle": (
            "Query is limited to QTL-interval candidates; the whole focal proteome and a small "
            "reference panel are reused as prebuilt search targets so out-of-interval paralogs are "
            "still detected. Profiles change threads/hit caps/tree policy only, not the labels."
        ),
    }
    manifest_path = out_dir / "orthology_tier1_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest | {"manifest": str(manifest_path)}


def _classify_paralogs(self_bitscore: float, paralogs: list[tuple[str, float, bool]], close_ratio: float) -> dict[str, object]:
    """Classify a candidate's focal self-search paralogs.

    ``paralogs`` is a list of (subject_gene_id, bitscore, in_interval) excluding
    the self hit. ``close_ratio`` sets how near a paralog's bitscore must be to
    the self hit to count as a redundancy risk.
    """
    close = [p for p in paralogs if self_bitscore > 0 and p[1] >= close_ratio * self_bitscore]
    close_out = [p for p in close if not p[2]]
    close_in = [p for p in close if p[2]]
    if not close:
        relationship = "single_copy_in_focal"
        reason = ""
    elif close_out:
        relationship = "out_of_interval_paralog"
        reason = f"out_of_interval_paralog(n={len(close_out)})"
    else:
        relationship = "in_interval_paralog"
        reason = f"in_interval_paralog(n={len(close_in)})"
    closest = max(paralogs, key=lambda p: p[1])[0] if paralogs else ""
    return {
        "relationship_class": relationship,
        "ambiguity_reason": reason,
        "focal_paralog_count": len(paralogs),
        "close_paralog_count": len(close),
        "out_of_interval_paralog_count": len(close_out),
        "closest_paralog": closest,
    }


def _classify_orthogroup(focal_count: int, reference_counts: dict[str, int]) -> str:
    reference_total = sum(reference_counts.values())
    if reference_total == 0:
        return "lineage_specific_or_focal_only"
    max_reference = max(reference_counts.values()) if reference_counts else 0
    if focal_count <= 1 and max_reference <= 1:
        return "one_to_one"
    if focal_count <= 1 and max_reference > 1:
        return "one_to_many"
    if focal_count > 1 and max_reference <= 1:
        return "many_to_one"
    return "many_to_many"


def orthofinder_import(
    orthogroups_arg: str,
    candidates_arg: str,
    focal_species: str,
    out_arg: str,
) -> dict[str, object]:
    """Import an OrthoFinder Orthogroups.tsv and classify each candidate's orthogroup."""
    frame = pd.read_csv(runtime_path(orthogroups_arg), sep="\t", dtype=str).fillna("")
    id_column = frame.columns[0]
    species_columns = [column for column in frame.columns if column != id_column]
    if focal_species not in species_columns:
        raise ValueError(f"focal species {focal_species!r} not among Orthogroups columns: {', '.join(species_columns)}")
    reference_species = [column for column in species_columns if column != focal_species]
    candidate_ids = _candidate_ids(candidates_arg)

    def split_cell(cell: str) -> list[str]:
        return [gene.strip() for gene in str(cell).replace(",", " ").split() if gene.strip()]

    gene_to_row: dict[str, pd.Series] = {}
    for _, row in frame.iterrows():
        for gene in split_cell(row[focal_species]):
            gene_to_row[gene] = row

    evidence_rows: list[dict[str, str]] = []
    summary_rows: list[dict[str, object]] = []
    for gene_id in sorted(candidate_ids):
        row = gene_to_row.get(gene_id)
        if row is None:
            summary_rows.append({"gene_id": gene_id, "orthogroup": "", "relationship": "unassigned", "focal_paralog_count": 0, "reference_species_hit": 0})
            continue
        og_id = str(row[id_column])
        focal_members = split_cell(row[focal_species])
        focal_count = len(focal_members)
        reference_counts = {species: len(split_cell(row[species])) for species in reference_species}
        relationship = _classify_orthogroup(focal_count, reference_counts)
        paralog_count = max(0, focal_count - 1)
        summary_rows.append(
            {
                "gene_id": gene_id,
                "orthogroup": og_id,
                "relationship": relationship,
                "focal_paralog_count": paralog_count,
                "reference_species_hit": sum(1 for count in reference_counts.values() if count > 0),
            }
        )
        if relationship != "lineage_specific_or_focal_only":
            evidence_rows.append(
                _row(
                    gene_id,
                    "orthogroup_supported",
                    "support",
                    "high" if relationship == "one_to_one" else "medium",
                    f"{og_id};relationship={relationship}",
                    "orthofinder",
                    str(runtime_path(orthogroups_arg)),
                    f"OrthoFinder orthogroup {og_id} relationship {relationship} across {sum(1 for c in reference_counts.values() if c)} reference species.",
                )
            )
        if paralog_count > 0 or relationship == "many_to_many":
            evidence_rows.append(
                _row(
                    gene_id,
                    "paralog_ambiguous",
                    "against",
                    "medium",
                    f"paralog_count={paralog_count};relationship={relationship}",
                    "orthofinder",
                    str(runtime_path(orthogroups_arg)),
                    "Orthogroup contains in-species duplicates or a many-to-many relationship; resolve paralogs before a 1:1 claim.",
                )
            )

    return _write_orthology_evidence(out_arg, evidence_rows, summary_rows, ["gene_id", "orthogroup", "relationship", "focal_paralog_count", "reference_species_hit"], "orthofinder", runtime_path(orthogroups_arg))


def omamer_import(omamer_arg: str, candidates_arg: str, out_arg: str) -> dict[str, object]:
    """Import an OMAmer HOG placement table into orthogroup evidence."""
    frame = pd.read_csv(runtime_path(omamer_arg), sep="\t", dtype=str, comment="!").fillna("")
    columns = {column.lower(): column for column in frame.columns}
    query_col = next((columns[key] for key in ("qseqid", "query", "queryid", "id") if key in columns), frame.columns[0])
    hog_col = next((columns[key] for key in ("hogid", "hog", "hog_id", "hogid_str") if key in columns), None)
    if hog_col is None:
        raise ValueError("OMAmer table must contain a HOG id column (hogid/hog)")
    candidate_ids = _candidate_ids(candidates_arg)

    evidence_rows: list[dict[str, str]] = []
    summary_rows: list[dict[str, object]] = []
    for _, row in frame.iterrows():
        gene_id = str(row[query_col])
        if gene_id not in candidate_ids:
            continue
        hog = str(row[hog_col]).strip()
        assigned = hog and hog.lower() not in {"na", "n/a", "none", "-"}
        summary_rows.append({"gene_id": gene_id, "hog_id": hog if assigned else "", "assigned": "yes" if assigned else "no"})
        if assigned:
            evidence_rows.append(
                _row(
                    gene_id,
                    "orthogroup_supported",
                    "support",
                    "medium",
                    hog,
                    "omamer_hog",
                    str(runtime_path(omamer_arg)),
                    f"OMAmer placed the candidate in HOG {hog}.",
                )
            )

    return _write_orthology_evidence(out_arg, evidence_rows, summary_rows, ["gene_id", "hog_id", "assigned"], "omamer", runtime_path(omamer_arg))


def synteny_import(anchors_arg: str, candidates_arg: str, out_arg: str, anchor_format: str = "anchors") -> dict[str, object]:
    """Import MCScanX/jcvi synteny anchors into synteny-supported evidence."""
    if anchor_format not in {"anchors", "collinearity"}:
        raise ValueError("anchor_format must be anchors or collinearity")
    candidate_ids = _candidate_ids(candidates_arg)
    path = runtime_path(anchors_arg)
    partners: dict[str, set[str]] = {}
    import re

    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if anchor_format == "collinearity":
            match = re.match(r"^\s*\d+-\s*\d+:\s*(\S+)\s+(\S+)", raw)
            if not match:
                continue
            gene_a, gene_b = match.group(1), match.group(2)
        else:
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            gene_a, gene_b = fields[0], fields[1]
        for gene, partner in ((gene_a, gene_b), (gene_b, gene_a)):
            if gene in candidate_ids:
                partners.setdefault(gene, set()).add(partner)

    evidence_rows: list[dict[str, str]] = []
    summary_rows: list[dict[str, object]] = []
    for gene_id in sorted(candidate_ids):
        gene_partners = sorted(partners.get(gene_id, set()))
        summary_rows.append({"gene_id": gene_id, "syntenic_partner_count": len(gene_partners), "syntenic_partners": ";".join(gene_partners)})
        if gene_partners:
            evidence_rows.append(
                _row(
                    gene_id,
                    "synteny_supported",
                    "support",
                    "medium",
                    ";".join(gene_partners),
                    "synteny",
                    str(path),
                    f"Candidate lies in a syntenic block with {', '.join(gene_partners)}.",
                )
            )

    return _write_orthology_evidence(out_arg, evidence_rows, summary_rows, ["gene_id", "syntenic_partner_count", "syntenic_partners"], "synteny", path)


def _write_orthology_evidence(out_arg: str, evidence_rows: list[dict[str, str]], summary_rows: list[dict[str, object]], summary_columns: list[str], source: str, input_path: Path) -> dict[str, object]:
    out_path = runtime_path(out_arg)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(evidence_rows, columns=EVIDENCE_COLUMNS).to_csv(out_path, sep="\t", index=False)
    summary_path = out_path.with_name(out_path.stem + f"_{source}_summary.tsv")
    pd.DataFrame(summary_rows, columns=summary_columns).to_csv(summary_path, sep="\t", index=False)
    return {
        "evidence": str(out_path),
        "summary": str(summary_path),
        "source": source,
        "evidence_count": len(evidence_rows),
        "candidate_count": len(summary_rows),
        "input_sha256": sha256(input_path),
        "principle": "Third-party orthology/synteny is normalized into the same canonical evidence schema regardless of profile.",
    }


def focal_paralog_import(
    focal_blast_arg: str,
    candidates_arg: str,
    out_arg: str,
    close_ratio: float = 0.4,
) -> dict[str, object]:
    """Import focal self-search hits and flag out-of-interval paralog redundancy."""
    candidate_ids = _candidate_ids(candidates_arg)
    hits = _best_hits(runtime_path(focal_blast_arg))
    evidence_rows: list[dict[str, str]] = []
    summary_rows: list[dict[str, object]] = []
    for gene_id in sorted(candidate_ids):
        gene_hits = hits.get(gene_id, [])
        self_hits = [h for h in gene_hits if str(h["subject"]) == gene_id]
        self_bitscore = float(self_hits[0]["bitscore"]) if self_hits else (float(gene_hits[0]["bitscore"]) if gene_hits else 0.0)
        paralogs = [
            (str(h["subject"]), float(h["bitscore"]), str(h["subject"]) in candidate_ids)
            for h in gene_hits
            if str(h["subject"]) != gene_id
        ]
        classification = _classify_paralogs(self_bitscore, paralogs, close_ratio)
        summary_rows.append({"gene_id": gene_id, **classification})
        if classification["out_of_interval_paralog_count"] > 0 or classification["close_paralog_count"] > 0:
            evidence_rows.append(
                _row(
                    gene_id,
                    "paralog_ambiguous",
                    "against",
                    "medium" if classification["relationship_class"] == "out_of_interval_paralog" else "low",
                    classification["ambiguity_reason"],
                    "focal_self_search",
                    str(runtime_path(focal_blast_arg)),
                    (
                        "Focal-proteome self search found a near-scoring paralog"
                        + (" outside the QTL interval" if classification["relationship_class"] == "out_of_interval_paralog" else " inside the QTL interval")
                        + f"; closest={classification['closest_paralog']}. Resolve orthology before trusting a 1:1 assignment."
                    ),
                )
            )

    out_path = runtime_path(out_arg)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(evidence_rows, columns=EVIDENCE_COLUMNS).to_csv(out_path, sep="\t", index=False)
    summary_path = out_path.with_name(out_path.stem + "_paralog_summary.tsv")
    pd.DataFrame(
        summary_rows,
        columns=["gene_id", "relationship_class", "ambiguity_reason", "focal_paralog_count", "close_paralog_count", "out_of_interval_paralog_count", "closest_paralog"],
    ).to_csv(summary_path, sep="\t", index=False)

    return {
        "evidence": str(out_path),
        "paralog_summary": str(summary_path),
        "candidate_count": len(candidate_ids),
        "out_of_interval_paralog_genes": int(sum(1 for row in summary_rows if row["out_of_interval_paralog_count"] > 0)),
        "relationship_counts": pd.Series([row["relationship_class"] for row in summary_rows], dtype=str).value_counts().to_dict() if summary_rows else {},
        "input_sha256": sha256(runtime_path(focal_blast_arg)),
        "principle": "Out-of-interval paralogs are recorded as redundancy risk; a best hit alone is never called a 1:1 ortholog.",
    }
