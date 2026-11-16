from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
from pathlib import Path

import pandas as pd
import yaml
from Bio import SeqIO

from .core import file_fingerprint, fingerprint_is_weak, runtime_path, run_tool


from .optional_deps import require_pysam
from .hardening import prioritize_evidence
from .orchestration import primer_plan
from .primerblast_import import primerblast_import
from .next_experiment import select_next_experiment
from .elimination import candidate_elimination
from .narrowing import integrated_narrowing
from .shortlist import screening_shortlist
from .reference_extract import extract_candidates
from .wet_state import update_wet_state
from .hypothesis_ledger import update_hypothesis_ledger
from .annotation_rescue import annotation_rescue
from .causal_scenarios import (
    SCENARIO_COLUMNS,
    load_qtl_effect_contract as load_scenario_effect_contract,
    load_scenario_rows,
    render_scenario_report,
    sign_concordance_evidence,
    write_scenario_ledger,
)
from .discordance_scenarios import comparison_fingerprint, discordance_review_rows
from .multigenome_private import multigenome_private
from .promoter_review import promoter_review
from .structural_variants import sv_plan
from .dual_route import dual_route_compare
from .evidence_dependency import (
    dependency_cluster_summary,
    dependency_reporting_counts,
    normalize_evidence_kind_column,
)
from .qtl_effect_contract import annotate_sign_concordance, contract_frame, load_qtl_effect_contract
from . import schemas


def _primerblast_root_fingerprint(root: Path) -> dict[str, object]:
    """Hash bounded package code/config, never a root directory or its DBs."""
    if not root.is_dir():
        raise ValueError(f"primerblast_root_resource must be a directory: {root}")
    package = root / "primerblast_oss"
    files = [root / name for name in ("pyproject.toml", "setup.cfg", "setup.py") if (root / name).is_file()]
    entries = 0
    if package.is_dir():
        for directory, dirs, names in os.walk(package, followlinks=False):
            dirs[:] = sorted(name for name in dirs if name not in {"__pycache__", ".git"})
            entries += len(dirs) + len(names)
            if entries > 4096:
                raise ValueError("primerblast_oss package exceeds bounded code fingerprint scope (4096 entries)")
            files.extend(Path(directory) / name for name in sorted(names) if name.endswith(".py"))
    if sum(path.stat().st_size for path in files) > 16 * 1024 * 1024:
        raise ValueError("primerblast_oss code/config exceeds fingerprint limit (16 MiB)")
    fingerprints = {str(path.relative_to(root)): file_fingerprint(path, full_hash=True) for path in sorted(files)}
    return {"path": str(root), "resource_type": "directory", "method": "bounded_package_code_sha256",
            "strength": "strong" if (package / "__main__.py").is_file() else "weak",
            "scope": "primerblast_oss/**/*.py and root build configuration; excludes DB/data files",
            "entrypoint_present": (package / "__main__.py").is_file(), "files": fingerprints,
            "sha256": hashlib.sha256(json.dumps(fingerprints, sort_keys=True).encode()).hexdigest()}


def _promoter_observations_to_hypotheses(observations: pd.DataFrame, case_id: str) -> pd.DataFrame:
    """Typed ``promoter_or_utr`` hypotheses from ``promoter_review()``'s review_observations.tsv (issue #230).

    Only parental_variant / sv_pav overlaps that promoter_review itself classified
    as ``promoter_or_utr`` (region-based promoter/5'UTR/3'UTR overlap; see
    ``_REGION_HYPOTHESIS_TYPE`` in promoter_review.py) are promoted into the
    locus-hypothesis ledger here. Reference-only motif context (never causal) and
    other typed observations promoter_review can also emit (``te_regulatory``,
    ``sv_cnv_pav`` from SV/PAV/TE calls landing in a promoter/UTR window) stay
    visible in review_observations.tsv for human review without an automatic
    hypothesis, keeping this wiring to exactly what issue #230 asks for.
    """
    if observations.empty:
        return pd.DataFrame(columns=schemas.LOCUS_HYPOTHESIS_COLUMNS)
    eligible = observations[
        observations["feature_type"].isin({"parental_variant", "sv_pav"})
        & observations["links_to_hypothesis_type"].eq("promoter_or_utr")
    ]
    rows: list[dict[str, str]] = []
    for record in eligible.to_dict(orient="records"):
        gene_id = str(record.get("gene_id", "")).strip()
        observation_id = str(record.get("observation_id", "")).strip()
        if not gene_id or not observation_id:
            continue
        row = {column: "" for column in schemas.LOCUS_HYPOTHESIS_COLUMNS}
        row.update({
            "case_id": case_id, "hypothesis_type": "promoter_or_utr",
            "interval_or_element_id": str(record.get("region_id", "")) or observation_id,
            "member_gene_ids": gene_id,
            "proposed_mechanism": f"{record.get('observation_class', 'regulatory variant')} overlapping the "
                                  f"{record.get('region_class', 'promoter/UTR')} of {gene_id}",
            "status": "proposed", "status_reason": "workflow seed from promoter/UTR review route",
            "created_by": "workflow:promoter_review", "proposal_origin": "workflow", "review_status": "pending_review",
            "source_id": f"promoter_review:{observation_id}", "gene_disposition": "retained", "rejection_scope": "none",
        })
        rows.append(row)
    frame = pd.DataFrame(rows, columns=schemas.LOCUS_HYPOTHESIS_COLUMNS)
    if frame.empty:
        return frame
    return frame.drop_duplicates(
        subset=["case_id", "hypothesis_type", "interval_or_element_id", "member_gene_ids"]
    ).sort_values(["member_gene_ids", "interval_or_element_id"], kind="stable").reset_index(drop=True)


def liftover_markers(hits_arg: str | None, out_dir_arg: str, max_product: int = 5000,
                     marker_fasta_arg: str | None = None, blast_db_arg: str | None = None, threads: int = 4,
                     max_block_gap: int = 1_000_000) -> dict[str, object]:
    """Convert BLAST marker hits into pair-level loci and safe interval blocks."""
    out = runtime_path(out_dir_arg); out.mkdir(parents=True, exist_ok=True)
    if marker_fasta_arg:
        if not blast_db_arg: raise ValueError("--marker-fasta requires --db")
        hits_path = out / "marker_blast.tsv"
        command = ["blastn", "-task", "blastn-short", "-query", str(runtime_path(marker_fasta_arg)), "-db", str(runtime_path(blast_db_arg)),
                   "-outfmt", "6 qseqid sseqid pident length qlen qstart qend sstart send evalue", "-num_threads", str(threads), "-out", str(hits_path)]
        completed = run_tool(command, text=True, capture_output=True)
        if completed.returncode: raise RuntimeError(f"marker blastn failed: {completed.stderr.strip()}")
    elif hits_arg:
        hits_path = runtime_path(hits_arg)
    else:
        raise ValueError("provide marker BLAST --hits or --marker-fasta with --db")
    first = hits_path.read_text(encoding="utf-8").splitlines()[:1]
    has_header = bool(first and first[0].split("\t")[0] == "marker")
    if has_header:
        hits = pd.read_csv(hits_path, sep="\t", dtype=str).fillna("")
    else:
        raw = pd.read_csv(hits_path, sep="\t", dtype=str, header=None, comment="#").fillna("")
        if raw.shape[1] < 10: raise ValueError("headerless marker BLAST requires 10 columns: qseqid sseqid pident length qlen qstart qend sstart send evalue")
        raw = raw.iloc[:, :10]; raw.columns = ["query","chrom","identity","length","qlen","qstart","qend","sstart","send","evalue"]
        parsed = raw["query"].str.extract(r"^(?P<marker>.+?)[_:\-.](?P<primer>F|R)$", expand=True)
        if parsed.isna().any().any(): raise ValueError("marker query IDs must end in _F/_R (or :F/:R, -F/-R)")
        hits = pd.DataFrame({"marker":parsed.marker,"primer":parsed.primer,"chrom":raw.chrom,"start":raw.sstart,"end":raw.send,
                             "strand":raw.apply(lambda row: "+" if int(row.sstart) <= int(row.send) else "-", axis=1),"identity":raw.identity,
                             "query_coverage":pd.to_numeric(raw.length) / pd.to_numeric(raw.qlen),"evalue":raw.evalue})
    required = {"marker", "primer", "chrom", "start", "end", "strand", "identity"}
    if not required.issubset(hits.columns):
        raise ValueError(f"marker hits missing columns: {', '.join(sorted(required - set(hits.columns)))}")
    hits[["start", "end"]] = hits[["start", "end"]].astype(int)
    hits["identity"] = hits["identity"].astype(float)
    pairs, selected = [], []
    for marker, group in hits.groupby("marker", sort=True):
        f_hits = group[(group["primer"].str.upper() == "F") & (group["identity"] >= 90)]
        r_hits = group[(group["primer"].str.upper() == "R") & (group["identity"] >= 90)]
        products = []
        for _, f in f_hits.iterrows():
            for _, r in r_hits.iterrows():
                if f.chrom != r.chrom or f.strand == r.strand:
                    continue
                start, end = min(f.start, f.end, r.start, r.end), max(f.start, f.end, r.start, r.end)
                if end - start + 1 <= max_product:
                    products.append((str(f.chrom), start, end))
        products = sorted(set(products))
        status = "unique_pair" if len(products) == 1 else ("multi_pair" if products else "no_valid_pair")
        pairs.append({"marker": marker, "pair_product_count": len(products), "pair_status": status,
                      "products": ";".join(f"{c}:{s}-{e}" for c, s, e in products)})
        if len(products) == 1:
            selected.append((marker, *products[0]))
    blocks = []
    for chrom in sorted({row[1] for row in selected}):
        rows = sorted((r for r in selected if r[1] == chrom), key=lambda x: x[2])
        for marker, _, start, end in rows:
            if not blocks or blocks[-1]["chrom"] != chrom or start - blocks[-1]["end"] > max_block_gap:
                blocks.append({"block": len(blocks) + 1, "chrom": chrom, "start": start, "end": end, "markers": [marker]})
            else:
                blocks[-1]["end"] = max(blocks[-1]["end"], end); blocks[-1]["markers"].append(marker)
    ordered_markers = [row[0] for row in sorted(selected, key=lambda x: (x[1], x[2]))]
    numeric_order = [int(match.group(1)) for marker in ordered_markers if (match := re.search(r"(\d+)$", marker))]
    marker_order_state = "unknown"
    if len(numeric_order) == len(ordered_markers) and len(numeric_order) > 1:
        marker_order_state = "monotonic" if numeric_order in (sorted(numeric_order), sorted(numeric_order, reverse=True)) else "non_monotonic"
    split = len(blocks) != 1
    pd.DataFrame(pairs).to_csv(out / "marker_pairs.tsv", sep="\t", index=False)
    block_rows = [{**b, "markers": ";".join(b["markers"]), "interval": f"{b['chrom']}:{b['start']}-{b['end']}"} for b in blocks]
    pd.DataFrame(block_rows).to_csv(out / "interval_blocks.tsv", sep="\t", index=False)
    result = {"route": "marker_liftover", "synteny_state": "split" if split else ("collinear" if blocks and marker_order_state != "non_monotonic" else ("rearranged" if blocks else "failed")),
              "recommended_intervals": [r["interval"] for r in block_rows], "single_span_allowed": len(blocks) == 1,
              "marker_order_state": marker_order_state, "pair_table": str(out / "marker_pairs.tsv"), "blocks": str(out / "interval_blocks.tsv"), "blast_hits": str(hits_path), "max_block_gap": max_block_gap}
    (out / "marker_liftover.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def extract_multi(fasta: str, gff3: str, intervals: list[str], out_dir: str, seqid_map: str | None, full_hash: bool,
                  gff_qc: str | None, flank: int = 0, provenance: str | None = None) -> dict[str, object]:
    if provenance:
        pre = json.loads(runtime_path(provenance).read_text(encoding="utf-8"))
        if not pre.get("positional_ranking_allowed") or pre.get("position_status") == "unresolved":
            raise ValueError("coordinate provenance gate rejected positional extraction")
        if pre.get("synteny_state") in {"split", "reverse", "failed"} and len(intervals) == 1:
            raise ValueError("non-collinear provenance requires explicit multiple --interval values")
    if len(intervals) == 1 and flank == 0:
        return extract_candidates(fasta, gff3, intervals[0], out_dir, seqid_map, full_hash, gff_qc)
    out = runtime_path(out_dir); out.mkdir(parents=True, exist_ok=True)
    frames = []
    for index, interval in enumerate(intervals, 1):
        part = out / f"block_{index:02d}"
        extract_candidates(fasta, gff3, interval, str(part), seqid_map, full_hash, gff_qc)
        frame = pd.read_csv(part / "candidates.tsv", sep="\t", dtype=str).fillna("")
        # interval_block only disambiguates *multiple* blocks; for a single
        # interval it keeps discover's extract byte-identical to `extract`.
        if len(intervals) > 1:
            frame["interval_block"] = index
        frames.append(frame)
    merged = pd.concat(frames, ignore_index=True).drop_duplicates("gene_id").sort_values(["chrom", "start", "gene_id"])
    merged.to_csv(out / "candidates.tsv", sep="\t", index=False, lineterminator="\n")
    for fasta_name in ("candidate_cds.fa", "candidate_transcript_cds.fa", "candidate_proteins.fa"):
        seen: set[str] = set()
        with (out / fasta_name).open("w", encoding="utf-8", newline="\n") as handle:
            for index in range(1, len(intervals) + 1):
                part_fasta = out / f"block_{index:02d}" / fasta_name
                with part_fasta.open(encoding="utf-8") as source:
                    for record in SeqIO.parse(source, "fasta"):
                        if record.id in seen: continue
                        seen.add(record.id); SeqIO.write(record, handle, "fasta")
    pysam = require_pysam('extract')
    fasta_file = pysam.FastaFile(str(runtime_path(fasta)))
    with (out / "candidate_genomic.fa").open("w", encoding="utf-8", newline="\n") as handle:
        for _, row in merged.iterrows():
            start, end = max(1, int(row.start) - flank), min(fasta_file.get_reference_length(row.fasta_seqid), int(row.end) + flank)
            seq = fasta_file.fetch(row.fasta_seqid, start - 1, end).upper()
            handle.write(f">{row.gene_id} genomic={row.chrom}:{start}-{end} flank={flank}\n{seq}\n")
    fasta_file.close()
    result = {"intervals": intervals, "interval_count": len(intervals), "candidate_count": len(merged), "flank": flank,
              "synteny_state": "split" if len(intervals) > 1 else "collinear", "candidates_tsv": str(out / "candidates.tsv"),
              "candidate_genomic_fasta": str(out / "candidate_genomic.fa"), "candidate_proteins_fasta": str(out / "candidate_proteins.fa"),
              "candidate_cds_fasta": str(out / "candidate_cds.fa")}
    (out / "extract_manifest.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def direction_evidence(evidence_arg: str, expected: str, out_arg: str) -> dict[str, object]:
    frame = pd.read_csv(runtime_path(evidence_arg), sep="\t", dtype=str).fillna("")
    if "phenotype_direction" not in frame.columns:
        raise ValueError("evidence table requires phenotype_direction")
    opposite = {"increase": "decrease", "decrease": "increase", "gain": "loss", "loss": "gain"}
    frame["direction_match"] = frame["phenotype_direction"].map(lambda x: "match" if x == expected else ("opposite" if x == opposite.get(expected) else "unknown"))
    frame.loc[frame["direction_match"] == "opposite", "direction"] = "against"
    frame.to_csv(runtime_path(out_arg), sep="\t", index=False, lineterminator="\n")
    return {"output": str(runtime_path(out_arg)), "matches": int((frame.direction_match == "match").sum()), "opposites": int((frame.direction_match == "opposite").sum())}


def dark_candidates(candidates_arg: str, out_arg: str, max_aa: int = 200) -> dict[str, object]:
    frame = pd.read_csv(runtime_path(candidates_arg), sep="\t", dtype=str).fillna("")
    no_homology = frame.get("homology_hit", pd.Series([""] * len(frame))).str.strip().eq("")
    complete = frame.get("annotation_quality", pd.Series([""] * len(frame))).eq("complete_orf")
    protein_length = frame["protein_length"] if "protein_length" in frame.columns else pd.Series(0, index=frame.index)
    short = pd.to_numeric(protein_length, errors="coerce").fillna(0).between(1, max_aa)
    result = frame[no_homology & complete & short].copy(); result["candidate_scope"] = "dark"; result["dark_reason"] = "in_interval_complete_short_orf_no_homology"
    result.to_csv(runtime_path(out_arg), sep="\t", index=False, lineterminator="\n")
    return {"output": str(runtime_path(out_arg)), "dark_candidate_count": len(result)}


def candidate_flags(candidates_arg: str, out_arg: str, tandem_distance: int = 100000, max_aa: int = 200) -> dict[str, object]:
    frame = pd.read_csv(runtime_path(candidates_arg), sep="\t", dtype=str).fillna("")
    frame["tandem_redundancy_risk"] = ""
    family_col = "paralog_family" if "paralog_family" in frame.columns else ("gene_name" if "gene_name" in frame.columns else None)
    if family_col and {"chrom", "start"}.issubset(frame.columns):
        starts = pd.to_numeric(frame["start"], errors="coerce")
        for _, group in frame[frame[family_col].str.strip().ne("")].groupby(["chrom", family_col]):
            indices = list(group.index)
            for idx in indices:
                if any(other != idx and abs(starts[other] - starts[idx]) <= tandem_distance for other in indices):
                    frame.loc[idx, "tandem_redundancy_risk"] = "high_tandem_paralog"
    # Dark classification depends on whether homology was actually searched. When
    # no homology_hit column exists yet (pre-annotation), a short complete ORF is
    # only *provisionally* dark; a confirmed no-homology call needs a completed
    # search. This prevents processing order alone from deciding Dark.
    homology_searched = "homology_hit" in frame.columns
    no_homology = frame.get("homology_hit", pd.Series([""] * len(frame))).str.strip().eq("")
    complete = frame.get("annotation_quality", pd.Series([""] * len(frame))).eq("complete_orf")
    protein_length = frame["protein_length"] if "protein_length" in frame.columns else pd.Series(0, index=frame.index)
    short = pd.to_numeric(protein_length, errors="coerce").fillna(0).between(1, max_aa)
    short_complete = no_homology & complete & short
    frame["candidate_scope"] = "positional"
    frame["dark_state"] = ""
    if homology_searched:
        frame.loc[short_complete, "candidate_scope"] = "dark"
        frame.loc[short_complete, "dark_state"] = "confirmed_no_homology"
    else:
        frame.loc[short_complete, "candidate_scope"] = "dark_pending"
        frame.loc[short_complete, "dark_state"] = "homology_not_searched"
    frame["assay_warning"] = ""; frame.loc[frame.tandem_redundancy_risk.ne(""), "assay_warning"] = "use paralog-aware or long-read assay"
    frame.to_csv(runtime_path(out_arg), sep="\t", index=False, lineterminator="\n")
    return {"output": str(runtime_path(out_arg)), "dark_candidate_count": int((frame.candidate_scope == "dark").sum()),
            "dark_pending_count": int((frame.candidate_scope == "dark_pending").sum()), "homology_searched": homology_searched,
            "tandem_flag_count": int(frame.tandem_redundancy_risk.ne("").sum())}


def _read_route_manifest_status(manifest_path: Path, route: str) -> tuple[bool, str]:
    """Read a route's explicit semantic-evaluability artifact fail closed.

    Route B records ``semantic_status`` in ``interval_first_status.json``.
    Route A's ``phenotype_first_status.json`` predates that field, so its
    explicit ``verified_bait_count`` is its evaluability condition.  Legacy TSV
    status tables with ``semantic_status`` remain accepted.  A malformed,
    incomplete, or semantically not-evaluated artifact is *not* a biological
    absence: it returns ``False`` and is represented downstream by the
    ``provider_missing`` scenario row.
    """
    try:
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        try:
            frame = pd.read_csv(manifest_path, sep="\t", dtype=str, keep_default_na=False)
        except (ValueError, OSError):
            return False, "unknown"
        if "semantic_status" not in frame.columns:
            return False, "unknown"
        statuses = frame["semantic_status"].str.strip()
        if statuses.empty or statuses.eq("").any():
            return False, "unknown"
        if statuses.eq("evaluable").all():
            return True, "evaluated"
        if statuses.str.startswith("not_evaluated").any():
            return False, "not_evaluated"
        return False, "unknown"
    except (OSError, UnicodeDecodeError):
        return False, "unknown"
    if not isinstance(payload, dict):
        return False, "unknown"
    semantic_status = str(payload.get("semantic_status", "")).strip()
    if semantic_status.startswith("not_evaluated"):
        return False, "not_evaluated"
    if route == "a":
        if str(payload.get("route", "")).strip() != "phenotype_first":
            return False, "unknown"
        try:
            verified_baits = int(payload.get("verified_bait_count", 0) or 0)
        except (TypeError, ValueError):
            return False, "unknown"
        return verified_baits >= 1, ("evaluated" if verified_baits >= 1 else "not_evaluated")
    if str(payload.get("route", "")).strip() != "interval_first":
        return False, "unknown"
    evaluated = semantic_status == "evaluable"
    return evaluated, ("evaluated" if evaluated else "not_evaluated")


def _route_table_evaluable(path: Path | None) -> tuple[bool, str]:
    """Compatibility bridge for a direct dual-route hypothesis artifact.

    A directly registered route table is explicit evidence that the legacy
    caller ran that route.  It is evaluable only when it has at least one row;
    An empty or unreadable direct table fails closed; a precomputed comparison
    with neither a direct table nor a status manifest retains the historic
    compatibility inference in the caller.
    """
    if path is None:
        return False, "absent"
    try:
        frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except (ValueError, OSError, UnicodeDecodeError):
        return False, "unknown"
    return (not frame.empty), ("evaluated" if not frame.empty else "not_evaluated")


def case_run(case_arg: str, resources_arg: str, out_dir_arg: str, resume: bool = False, mode: str = "deterministic", strict_fingerprints: bool = False, *, narrowing_overrides: dict[str, str] | None = None) -> dict[str, object]:
    """Run a deterministic, source-backed case with step-level resumability."""
    case_path, resources_path, out = runtime_path(case_arg), runtime_path(resources_arg), runtime_path(out_dir_arg)
    case = yaml.safe_load(case_path.read_text(encoding="utf-8")) or {}
    registry = yaml.safe_load(resources_path.read_text(encoding="utf-8")) or {}
    if not isinstance(case, dict) or not case.get("case_id"):
        raise ValueError("case manifest requires case_id")
    if mode not in {"deterministic", "llm-assisted"}: raise ValueError("mode must be deterministic or llm-assisted")
    resource_map = registry.get("resources", registry)
    if not isinstance(resource_map, dict):
        raise ValueError("resource registry must be an object or contain resources")
    out.mkdir(parents=True, exist_ok=True)
    state_path = out / "case_state.json"
    state = json.loads(state_path.read_text(encoding="utf-8")) if resume and state_path.exists() else {"steps": {}}

    def resolve(value: object) -> Path | None:
        if not value:
            return None
        entry = resource_map.get(str(value), value)
        raw = entry.get("path") if isinstance(entry, dict) else entry
        path = runtime_path(str(raw))
        if not path.is_absolute(): path = resources_path.parent / path
        if not path.exists():
            raise ValueError(f"resource {value!r} does not exist: {path}")
        return path

    def resolve_narrowing(name: str, resource_key: str) -> Path | None:
        override = (narrowing_overrides or {}).get(name)
        if override:
            path = runtime_path(override)
            if not path.exists():
                raise ValueError(f"guided-profile narrowing override does not exist: {path}")
            return path
        return resolve(inputs.get(resource_key))

    def digest(paths: list[Path], config: object) -> str:
        hasher = hashlib.sha256(json.dumps(config, sort_keys=True, ensure_ascii=False).encode())
        for path in paths:
            stat = path.stat(); hasher.update(str(path).encode()); hasher.update(str(stat.st_size).encode()); hasher.update(str(stat.st_mtime_ns).encode())
            hasher.update(json.dumps(file_fingerprint(path, full_hash=strict_fingerprints), sort_keys=True, ensure_ascii=False).encode())
        return hasher.hexdigest()

    def output_fingerprints(outputs: list[Path]) -> dict[str, dict[str, object]]:
        # These are the stage's declared generated files, not its large input
        # references. Content must match even when size and mtime are unchanged.
        fingerprints = {}
        for path in outputs:
            fingerprint = file_fingerprint(path, full_hash=True)
            fingerprints[str(path)] = {key: fingerprint[key] for key in ("size_bytes", "sha256")}
        return fingerprints

    def cached(name: str, step_hash: str, outputs: list[Path]) -> bool:
        previous = state.get("steps", {}).get(name, {})
        if not resume or previous.get("input_hash") != step_hash:
            return False
        declared_outputs = previous.get("outputs")
        previous_fingerprints = previous.get("output_fingerprints")
        # Older manifests cannot prove output integrity. Recompute once rather
        # than trusting and recording whatever bytes happen to exist on resume.
        if (not isinstance(declared_outputs, list)
                or not all(isinstance(path, str) for path in declared_outputs)
                or not isinstance(previous_fingerprints, dict)
                or set(declared_outputs) != set(previous_fingerprints)
                or not {str(path) for path in outputs}.issubset(previous_fingerprints)):
            return False
        try:
            # Some stages (e.g. primer plans) declare additional outputs only
            # at completion, so validate the whole previous set, not a subset.
            return output_fingerprints([Path(path) for path in declared_outputs]) == previous_fingerprints
        except OSError:
            return False

    def complete(name: str, step_hash: str, outputs: list[Path], extra: dict | None = None) -> None:
        state.setdefault("steps", {})[name] = {
            "status": "complete", "input_hash": step_hash, "outputs": [str(p) for p in outputs],
            "output_fingerprints": output_fingerprints(outputs), **(extra or {}),
        }
        state_path.write_text(json.dumps(state, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    def write_stable(path: Path, text: str) -> None:
        if not path.exists() or path.read_text(encoding="utf-8") != text:
            path.parent.mkdir(parents=True, exist_ok=True); path.write_text(text, encoding="utf-8")

    source = case.get("source") or (case.get("sources") or [{}])[0]
    target = case.get("target") or {}
    inputs = case.get("inputs") or {}
    assay = case.get("assay") or {}
    interval_value = source.get("target_intervals") or source.get("target_interval") or target.get("interval")
    intervals = [str(x) for x in interval_value] if isinstance(interval_value, list) else ([str(interval_value)] if interval_value else [])
    marker_hits = resolve(source.get("marker_hits_resource") or inputs.get("marker_hits_resource"))
    qtlift_result_path = resolve(source.get("qtlift_result_resource"))
    route = "validated_target_interval" if intervals else ("marker_liftover" if marker_hits else "unresolved_position")
    review_rows: list[dict[str, str]] = []

    coordinate_path = out / "coordinate_provenance.json"
    coordinate_hash = digest([p for p in (marker_hits, qtlift_result_path) if p], {"source": source, "intervals": intervals})
    if not cached("coordinate", coordinate_hash, [coordinate_path]):
        marker_result = None
        qtlift_result = json.loads(qtlift_result_path.read_text(encoding="utf-8")) if qtlift_result_path else None
        if marker_hits:
            marker_result = liftover_markers(str(marker_hits), str(out / "marker_liftover"), int(source.get("max_marker_product", 5000)))
            intervals = marker_result["recommended_intervals"]
            if not intervals:
                route = "unresolved_position"
            elif marker_result["synteny_state"] == "split":
                route = "marker_liftover_split"
        allowed = bool(intervals)
        marker_provisional = bool(marker_result and marker_result.get("synteny_state") in {"split", "rearranged"})
        qtlift_final = (qtlift_result or {}).get("final") or {}
        qtlift_interval = f"{qtlift_final.get('contig')}:{qtlift_final.get('start')}-{qtlift_final.get('end')}" if qtlift_final.get("contig") and qtlift_final.get("start") and qtlift_final.get("end") else ""
        qtlift_confidence = str((qtlift_result or {}).get("confidence", ""))
        qtlift_safe = qtlift_confidence.casefold() in {"high", "validated"} and str((qtlift_result or {}).get("synteny_state", "")).casefold() not in {"split", "reverse", "failed"}
        if not marker_result and qtlift_interval and qtlift_safe:
            intervals = [qtlift_interval]; allowed = True; route = "qtlift_interval"
        coordinate = {"source_id": source.get("source_id", ""), "source_assembly": source.get("assembly", source.get("source_reference_name", "")),
                      "source_interval": source.get("interval", ""), "target_assembly": target.get("assembly", ""), "target_intervals": intervals,
                      "coordinate_method": route, "coordinate_basis": source.get("coordinate_basis", "unknown"),
                      "liftover_confidence": source.get("liftover_confidence", "medium_pair_supported" if marker_result else ("validated" if route == "validated_target_interval" else "")),
                      "marker_order_state": marker_result.get("marker_order_state", "") if marker_result else "",
                      "marker_pair_audit": marker_result.get("pair_table", "") if marker_result else "",
                      "qtlift_job_id": (qtlift_result or {}).get("job_id", ""), "qtlift_interval": qtlift_interval,
                      "qtlift_confidence": qtlift_confidence, "qtlift_synteny_state": (qtlift_result or {}).get("synteny_state", ""),
                      "qtlift_warnings": (qtlift_result or {}).get("warnings", []),
                      "literature_search_status": source.get("literature_search_status", "unknown"), "position_status": ("provisional_interval" if marker_provisional else "validated_interval") if allowed else "unresolved",
                      "positional_ranking_allowed": allowed}
        coordinate_path.write_text(json.dumps(coordinate, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        complete("coordinate", coordinate_hash, [coordinate_path], {"route": route})
    coordinate = json.loads(coordinate_path.read_text(encoding="utf-8")); intervals = coordinate["target_intervals"]; allowed = bool(coordinate["positional_ranking_allowed"])
    if coordinate.get("literature_search_status") != "complete":
        review_rows.append({"reason": "literature_search_incomplete", "target": str(source.get("source_id", "case source")), "required_next_action": "complete source interval/marker literature search", "allowed_choices": "complete search | document partial search | keep rescue-only"})
    if coordinate.get("marker_order_state") == "non_monotonic":
        review_rows.append({"reason":"marker_order_non_monotonic","target":str(source.get("source_id", "case source")),"required_next_action":"review split blocks independently; do not collapse to min-max span","allowed_choices":"accept per-block extraction | obtain synteny evidence | stop positional ranking"})
    if coordinate.get("qtlift_confidence", "").casefold() in {"low", "medium", "manual check"} or coordinate.get("qtlift_synteny_state", "").casefold() in {"split", "reverse", "failed"}:
        review_rows.append({"reason":"qtlift_low_or_rearranged","target":str(coordinate.get("qtlift_job_id", "QTLift")),"required_next_action":"retain per-block marker route or obtain additional synteny evidence; never collapse the QTLift span silently","allowed_choices":"use marker blocks | manual synteny review | stop positional ranking"})

    candidates_path = out / "candidates.tsv"
    existing_candidates = resolve(inputs.get("candidates_resource"))
    fasta = resolve(target.get("fasta_resource")); gff3 = resolve(target.get("gff3_resource"))
    extract_inputs = [p for p in (existing_candidates, fasta, gff3, coordinate_path) if p]
    extract_hash = digest(extract_inputs, {"intervals": intervals, "flank": assay.get("flank", 0)})
    if not cached("extract", extract_hash, [candidates_path]):
        if existing_candidates:
            shutil.copyfile(existing_candidates, candidates_path)
        elif allowed and fasta and gff3:
            extract_result = extract_multi(str(fasta), str(gff3), intervals, str(out / "extract"), None, False, None, int(assay.get("flank", 0)), str(coordinate_path))
            shutil.copyfile(extract_result["candidates_tsv"], candidates_path)
        elif not allowed:
            rescue_candidates = resolve(inputs.get("rescue_candidates_resource"))
            if rescue_candidates:
                shutil.copyfile(rescue_candidates, candidates_path)
            else:
                pd.DataFrame(columns=["gene_id"]).to_csv(candidates_path, sep="\t", index=False, lineterminator="\n")
        else:
            raise ValueError("candidate extraction requires target fasta_resource and gff3_resource, or inputs.candidates_resource")
        complete("extract", extract_hash, [candidates_path])
    candidates = pd.read_csv(candidates_path, sep="\t", dtype=str).fillna("")
    if "gene_id" not in candidates.columns or candidates["gene_id"].duplicated().any():
        raise ValueError("candidate table requires unique gene_id")

    annotation_rescue_dir = out / "annotation_rescue"
    annotation_rescue_table = None
    rescue_config = case.get("annotation_rescue") or {}
    rescue_queries = resolve(inputs.get("annotation_rescue_query_proteins_resource"))
    rescue_query_metadata = resolve(inputs.get("annotation_rescue_query_metadata_resource"))
    rescue_projection = resolve(inputs.get("annotation_rescue_projection_gff_resource"))
    if rescue_queries and not rescue_query_metadata:
        raise ValueError("annotation rescue queries require inputs.annotation_rescue_query_metadata_resource")
    if rescue_queries:
        if len(intervals) != 1 or not fasta:
            review_rows.append({"reason":"annotation_rescue_requires_single_interval_fasta","target":str(rescue_queries),"required_next_action":"run projection per validated interval block against the target FASTA","allowed_choices":"provide one validated interval | run each block separately | retain reference-missing review"})
        elif set(schemas.CANDIDATE_SCHEMA.required) - set(candidates.columns):
            review_rows.append({"reason":"annotation_rescue_native_coordinates_missing","target":str(candidates_path),"required_next_action":"provide canonical extract candidates with genomic coordinates","allowed_choices":"re-extract from target GFF3 | retain unmapped review"})
        else:
            rescue_outputs = [
                annotation_rescue_dir / "annotation_rescue.tsv", annotation_rescue_dir / "annotation_rescue_review.tsv",
                annotation_rescue_dir / "annotation_rescue_review_overflow.tsv", annotation_rescue_dir / "annotation_rescue_manifest.json",
                annotation_rescue_dir / "run_local_annotation_rescue.sh",
            ]
            rescue_inputs = [p for p in (fasta, candidates_path, rescue_queries, rescue_query_metadata, rescue_projection) if p]
            execute_local = bool(rescue_config.get("execute_local", False))
            if rescue_projection and execute_local:
                raise ValueError("annotation rescue must choose either a projection resource or execute_local")
            miniprot_name = str(rescue_config.get("miniprot", "miniprot"))
            # Only bind the miniprot binary into the cache key when we would run
            # it. When importing a precomputed projection GFF, the local miniprot
            # version is irrelevant and must not churn the resume cache.
            miniprot_path = shutil.which(miniprot_name) if execute_local else None
            miniprot_version = ""
            if execute_local and miniprot_path:
                version_result = run_tool([miniprot_path, "--version"], capture_output=True, text=True, check=False)
                miniprot_version = (version_result.stdout or version_result.stderr).strip().splitlines()[0] if (version_result.stdout or version_result.stderr).strip() else "unknown"
            rescue_hash = digest(rescue_inputs, {"interval":intervals[0], "config":rescue_config, "miniprot_path":miniprot_path, "miniprot_version":miniprot_version})
            if cached("annotation_rescue", rescue_hash, rescue_outputs):
                annotation_rescue_table = annotation_rescue_dir / "annotation_rescue.tsv"
            else:
                # The rescue lane is an optional review appendix. A failure here
                # (bad projection GFF, metadata gap, schema violation) must never
                # abort the deterministic primary pipeline; it degrades to a
                # surfaced review item instead.
                try:
                    annotation_rescue(
                        str(fasta), str(candidates_path), str(rescue_queries), str(rescue_query_metadata), intervals[0], str(annotation_rescue_dir),
                        projection_gff_arg=str(rescue_projection) if rescue_projection else None,
                        projection_complete=bool(rescue_config.get("projection_complete", False)),
                        execute=execute_local,
                        miniprot_executable=miniprot_name,
                        threads=int(rescue_config.get("threads", 4)), full_hash=strict_fingerprints,
                        review_limit=int(rescue_config.get("review_limit", 20)),
                    )
                except (ValueError, OSError) as error:
                    review_rows.append({"reason":"annotation_rescue_failed","target":str(rescue_queries),
                                        "required_next_action":f"rescue lane did not complete; inspect projection GFF / query metadata: {error}",
                                        "allowed_choices":"fix projection inputs and rerun | retain reference-missing review | proceed with deterministic rank only"})
                else:
                    complete("annotation_rescue", rescue_hash, rescue_outputs)
                    annotation_rescue_table = annotation_rescue_dir / "annotation_rescue.tsv"
    flagged_candidates_path = out / "candidate_flags.tsv"
    flags_hash = digest([candidates_path], {"tandem_distance": case.get("tandem_distance", 100000), "dark_max_aa": case.get("dark_max_aa", 200)})
    if not cached("candidate_flags", flags_hash, [flagged_candidates_path]):
        candidate_flags(str(candidates_path), str(flagged_candidates_path), int(case.get("tandem_distance", 100000)), int(case.get("dark_max_aa", 200)))
        complete("candidate_flags", flags_hash, [flagged_candidates_path])

    # Private-variant funnel: a formal stage between annotation and ranking.
    # It never writes evidence or a score; it produces the private/shared/uncallable
    # event ledger, the projection and callable-CDS denominators, and its own review
    # items. When it is not run, the manifest and report say `not_evaluated` rather
    # than leaving the absence of private variants to be read as a negative result.
    private_dir = out / "private_variants"
    private_model = case.get("private_variant_model") or {}
    private_status = "not_requested" if not private_model else "not_evaluated"
    private_focal_vcf = None
    private_outputs = [
        private_dir / name for name in (
            "focal_event_sharing_ledger.tsv", "candidate_control_matrix.tsv", "control_callability_sensitivity.tsv",
            "event_class_summary.tsv", "gene_private_variant_ledger.tsv", "projection_and_callable_denominators.tsv",
            "working_candidate_panel.tsv", "review_required.tsv", "private_variant_manifest.json")
    ]
    if private_model:
        private_control_vcfs: dict[str, Path] = {}
        private_callable_beds: dict[str, Path] = {}
        projection_ledger = reference_consequence = native_consequence = native_genes = None
        try:
            private_focal_vcf = resolve(inputs.get("private_focal_vcf_resource"))
            private_control_vcfs = {str(name): resolve(value) for name, value in (inputs.get("private_control_vcf_resources") or {}).items()}
            private_callable_beds = {str(name): resolve(value) for name, value in (inputs.get("private_callable_bed_resources") or {}).items()}
            for name, value in (private_model.get("callability_resources") or {}).items():
                # Resolve manifest-declared callability through the registry too, so a
                # resource id never reaches the funnel as a literal path.
                private_callable_beds.setdefault(str(name), resolve(value))
            projection_ledger = resolve(inputs.get("projection_ledger_resource"))
            reference_consequence = resolve(inputs.get("reference_frame_consequence_resource"))
            native_consequence = resolve(inputs.get("native_frame_consequence_resource"))
            native_genes = resolve(inputs.get("native_annotation_genes_resource"))
        except ValueError as error:
            private_focal_vcf = None
            review_rows.append({"reason":"private_variant_inputs_unresolved","target":str(case["case_id"]),
                                "required_next_action":f"register the declared private-variant resources: {error}",
                                "allowed_choices":"register resources | remove private_variant_model | keep the stage not_evaluated"})
        if private_focal_vcf and gff3:
            private_inputs = [private_focal_vcf, gff3] + sorted(private_control_vcfs.values()) + sorted(private_callable_beds.values())
            private_inputs += [path for path in (projection_ledger, reference_consequence, native_consequence, native_genes) if path]
            private_hash = digest(private_inputs, {"model": private_model, "intervals": intervals})
            if cached("private_variants", private_hash, private_outputs):
                private_status = "evaluated"
            else:
                # A defect in this lane must not abort the deterministic pipeline;
                # it degrades to a surfaced review item, like annotation rescue.
                try:
                    multigenome_private(
                        str(private_dir), str(private_focal_vcf), str(gff3), case_arg=str(case_path),
                        control_vcf_args=[f"{name}={path}" for name, path in sorted(private_control_vcfs.items())],
                        callable_bed_args=[f"{name}={path}" for name, path in sorted(private_callable_beds.items())],
                        projection_ledger_arg=str(projection_ledger) if projection_ledger else None,
                        native_genes_arg=str(native_genes) if native_genes else None,
                        reference_consequence_arg=str(reference_consequence) if reference_consequence else None,
                        native_consequence_arg=str(native_consequence) if native_consequence else None,
                        intervals=intervals or None, case_id=str(case["case_id"]), full_hash=strict_fingerprints,
                    )
                except (ValueError, OSError) as error:
                    review_rows.append({"reason":"private_variant_stage_failed","target":str(private_focal_vcf),
                                        "required_next_action":f"private-variant funnel did not complete; inspect the callsets and callability inputs: {error}",
                                        "allowed_choices":"fix inputs and rerun | keep the stage not_evaluated"})
                else:
                    complete("private_variants", private_hash, private_outputs)
                    private_status = "evaluated"
        elif not private_focal_vcf:
            review_rows.append({"reason":"private_variant_focal_callset_missing","target":str(case["case_id"]),
                                "required_next_action":"register inputs.private_focal_vcf_resource, or remove private_variant_model",
                                "allowed_choices":"register the focal VCF | remove private_variant_model"})
    if private_status != "evaluated":
        review_rows.append({"reason":"private_variant_not_evaluated","target":str(case["case_id"]),
                            "required_next_action":"private/shared status was not evaluated; do not read a low private-variant count as biological evidence",
                            "allowed_choices":"run multigenome-private | declare private_variant_model | keep not_evaluated"})

    evidence_path = out / "evidence.tsv"
    evidence_resources = [resolve(value) for value in inputs.get("evidence_resources", [])]
    evidence_resources = [p for p in evidence_resources if p]
    evidence_hash = digest(evidence_resources + [candidates_path], {"phenotype_direction": (case.get("trait") or {}).get("phenotype_direction", "")})
    excluded_path = out / "evidence_excluded.tsv"
    evidence_columns = ["gene_id","evidence_code","direction","confidence","value","source_type","source_id","description"]
    evidence_dependency_columns = ["evidence_kind", "root_source_id"]
    evidence_clusters_path = out / "evidence_dependency_clusters.tsv"
    evidence_dependency_manifest_path = out / "evidence_dependency_manifest.json"
    evidence_provider_failed = False
    if not cached("evidence", evidence_hash, [evidence_path, excluded_path, evidence_clusters_path, evidence_dependency_manifest_path]):
        # Issue #346: an evidence provider that cannot be read (missing entry,
        # unreadable file) must not abort the connection-diagnostic run. Like the
        # other optional lanes it degrades to a surfaced review item, and the
        # rank lane records an explicit not_evaluated_evidence_provider_failed
        # state instead of silently ranking on zero rows.
        # The provider-failed placeholder files keep the same column shape as
        # the attempted path: the 8 core columns plus any optional dependency
        # columns that were present in the input frames.
        present_dependency_columns: list[str] = []
        for resource_path in evidence_resources:
            try:
                header = tuple(pd.read_csv(resource_path, sep="\t", dtype=str, nrows=0).columns)
            except (ValueError, OSError):
                continue
            for column in evidence_dependency_columns:
                if column in header and column not in present_dependency_columns:
                    present_dependency_columns.append(column)
        try:
            frames = [pd.read_csv(path, sep="\t", dtype=str).fillna("") for path in evidence_resources]
            evidence = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=evidence_columns)
            for column in evidence_columns:
                if column not in evidence.columns: evidence[column] = ""
            evidence = evidence[evidence.gene_id.isin(set(candidates.gene_id))].copy()
            excluded = evidence[evidence.source_id.str.strip().eq("")].copy()
            evidence = evidence[evidence.source_id.str.strip().ne("")].copy()
            expected = str((case.get("trait") or {}).get("phenotype_direction", ""))
            if expected and "phenotype_direction" in evidence.columns:
                opposite = {"increase":"decrease","decrease":"increase","gain":"loss","loss":"gain"}
                mask = evidence.phenotype_direction.eq(opposite.get(expected, "__none__")); evidence.loc[mask, "direction"] = "against"
            # Issue #364: evidence.tsv feeds the rank step, which validates
            # evidence_kind against the strict EVIDENCE_DEPENDENCY_KINDS enum
            # (hardening.prioritize_evidence -> schemas.validate). Normalize
            # human-readable kinds (casefold, space -> underscore, blank ->
            # unknown) with the dependency module's own rule BEFORE the write
            # and fail closed on anything still out of vocabulary, so a
            # documented-acceptable input can never crash the rank step.
            # root_source_id stays free text (no enum).
            if "evidence_kind" in evidence.columns:
                evidence["evidence_kind"] = normalize_evidence_kind_column(evidence)
            # Optional dependency columns survive the canonical 8-column write
            # when any input resource carries them; column order keeps the 8
            # core columns first.
            evidence[evidence_columns + [c for c in evidence_dependency_columns if c in evidence.columns]].sort_values(["gene_id","evidence_code","source_id","direction"], kind="stable").to_csv(evidence_path, sep="\t", index=False, lineterminator="\n")
            excluded.to_csv(excluded_path, sep="\t", index=False, lineterminator="\n")
            dependency_counts = dependency_reporting_counts(evidence)
            dependency_cluster_summary(evidence).to_csv(evidence_clusters_path, sep="\t", index=False, lineterminator="\n")
            dependency_manifest = {
                "schema_version": schemas.SCHEMA_VERSION,
                "command": "case-run",
                "input_fingerprint": {str(path): file_fingerprint(path) for path in evidence_resources},
                "dependency_status": "evaluated",
                **dependency_counts,
                "cluster_summary": str(evidence_clusters_path.relative_to(out)),
            }
            write_stable(evidence_dependency_manifest_path, json.dumps(dependency_manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n")
            complete("evidence", evidence_hash, [evidence_path, excluded_path, evidence_clusters_path, evidence_dependency_manifest_path],
                     {"included": len(evidence), "excluded": len(excluded), **dependency_counts,
                      "evidence_dependency_manifest": str(evidence_dependency_manifest_path),
                      "evidence_dependency_clusters": str(evidence_clusters_path)})
        except (ValueError, OSError) as error:
            evidence_provider_failed = True
            if "unsupported evidence_kind" in str(error):
                review_rows.append({"reason": "evidence_kind_invalid",
                                    "target": str(evidence_resources) if evidence_resources else str(inputs.get("evidence_resources", [])),
                                    "required_next_action": f"fix evidence_kind values in the evidence resource to one of {sorted(schemas.EVIDENCE_DEPENDENCY_KINDS)}; ranking has no scoring evidence",
                                    "allowed_choices": "fix the evidence resource and rerun | keep the rank lane positional-fallback"})
            else:
                review_rows.append({"reason": "evidence_provider_failed",
                                    "target": str(evidence_resources) if evidence_resources else str(inputs.get("evidence_resources", [])),
                                    "required_next_action": f"evidence provider could not be read ({type(error).__name__}: {error}); ranking has no scoring evidence",
                                    "allowed_choices": "fix the evidence resource and rerun | keep the rank lane positional-fallback"})
            placeholder_columns = evidence_columns + present_dependency_columns
            pd.DataFrame(columns=placeholder_columns).to_csv(evidence_path, sep="\t", index=False, lineterminator="\n")
            pd.DataFrame(columns=placeholder_columns).to_csv(excluded_path, sep="\t", index=False, lineterminator="\n")
            dependency_cluster_summary(pd.DataFrame(columns=placeholder_columns)).to_csv(evidence_clusters_path, sep="\t", index=False, lineterminator="\n")
            dependency_manifest = {
                "schema_version": schemas.SCHEMA_VERSION,
                "command": "case-run",
                "input_fingerprint": {},
                "dependency_status": "not_evaluated",
                "evidence_count": 0,
                "independent_evidence_cluster_count": 0,
                "dependency_unknown_evidence_count": 0,
                "dependency_unknown_present": False,
                "rank_contribution": "none",
            }
            write_stable(evidence_dependency_manifest_path, json.dumps(dependency_manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n")

    # complete() flattens its ``extra`` dict into the step record, so the
    # evidence step's recorded counts/artifact paths live beside status/
    # input_hash/outputs/output_fingerprints.
    evidence_step_extra = {
        key: value for key, value in (state.get("steps", {}).get("evidence") or {}).items()
        if key not in ("status", "input_hash", "outputs", "output_fingerprints")
    }

    rank_dir = out / "rank"; ranked_path = rank_dir / "ranked_candidates.tsv"; contributions_path = rank_dir / "score_contributions.tsv"
    rank_hash = digest([flagged_candidates_path, evidence_path], {"ranking_policy": case.get("ranking_policy", "deterministic_v1"), "min_top_score": float(case.get("min_top_score", 2.0)), "top_n": int(assay.get("top_n", 3)), "reserve_n": int(assay.get("reserve_n", 3)), "rescue_n": int(assay.get("rescue_n", 3)), "assay": assay, "allowed": allowed})
    if not cached("rank", rank_hash, [ranked_path, contributions_path]):
        prioritize_evidence(str(flagged_candidates_path), str(evidence_path), str(rank_dir), None, int(assay.get("top_n", 3)), int(assay.get("reserve_n", 3)), float(case.get("min_top_score", 2.0)), int(assay.get("rescue_n", 3)))
        ranked = pd.read_csv(ranked_path, sep="\t", dtype=str).fillna("")
        if not allowed:
            ranked["dry_group"] = "rescue"; ranked["top_eligibility"] = "position_unresolved"
            ranked.to_csv(ranked_path, sep="\t", index=False, lineterminator="\n")
        complete("rank", rank_hash, [ranked_path, contributions_path])

    # Issue #346: an unwired evidence lane silently degenerates the rank into a
    # coordinate-ordered tie (dry_score=0.0, screening_rank=1 for every row) and
    # the run still reports REPORT_READY. The run completes exactly as before,
    # but the run-level state now says so explicitly instead of leaving the tie
    # to be read as biology. The status is derived from the *scored* contribution
    # table (score_contributions.tsv), so evidence rows that exist but contribute
    # no score (rescue-only or unweighted codes) stay distinguishable from a
    # genuinely informative rank lane.
    included_evidence = pd.read_csv(evidence_path, sep="\t", dtype=str).fillna("")
    scored_contributions = pd.read_csv(contributions_path, sep="\t", dtype=str).fillna("")
    evidence_present = len(included_evidence) > 0
    rank_is_positional_fallback = len(scored_contributions) == 0
    # Issue #397: a run whose evidence rows exist but carry no support-direction
    # row contributes only penalties and positional ties; it is declared
    # penalty-only instead of being read as evaluated biology. The test is the
    # evidence rows themselves (``direction == "support"`` present >= 1),
    # independent of scored_contributions: a support row is evaluable even if
    # the registry cannot score it, and against/context/rescue-only evidence
    # never is. support rows on the truth gene are irrelevant here.
    supporting_evidence_rows = int((included_evidence["direction"] == "support").sum())
    rank_is_penalty_only = evidence_present and not evidence_provider_failed and supporting_evidence_rows == 0
    if evidence_provider_failed:
        rank_semantic_status = "not_evaluated_evidence_provider_failed"
        rank_lane_warning = "rank lane is non-informative: evidence provider failed (positional fallback)"
    elif rank_is_penalty_only:
        rank_semantic_status = "not_evaluated_no_supporting_evidence"
        rank_lane_warning = "rank lane is non-informative: evidence present but no supporting evidence (penalty-only)"
    elif evidence_present:
        rank_semantic_status = "evaluable"
        rank_lane_warning = (
            "" if not rank_is_positional_fallback else
            "rank lane is non-informative: evidence present but no scored contribution (positional fallback)"
        )
    else:
        rank_semantic_status = "not_evaluated_no_scoring_evidence"
        rank_lane_warning = "rank lane is non-informative: no scoring evidence (positional fallback)"
    schemas.require_rank_semantic_status(rank_semantic_status)
    rank_metric_eligible = evidence_present and not evidence_provider_failed and supporting_evidence_rows > 0
    stage_a_is_budget_only = rank_is_positional_fallback or rank_is_penalty_only

    # Issue #330: the frozen Stage-A cut. The README describes
    # `ranking -> frozen Stage-A 20 -> review` as one chain, but screening-shortlist
    # was reachable only as a standalone CLI, so no orchestrated run ever produced
    # the frozen shortlist. This stage closes that arrow.
    #
    # It is deliberately a separate cut from rank's own top_n/reserve_n/rescue_n
    # grouping: `dry_group` answers "is the evidence sufficient to promote a claim",
    # screening-shortlist answers "which candidates are economical to review next"
    # (shortlist.py:20-30). The two are not interchangeable, so both are emitted.
    stage_a_dir = out / "stage_a"
    stage_a_outputs = [
        stage_a_dir / "deterministic_shortlist.tsv", stage_a_dir / "rescue_appendix.tsv",
        stage_a_dir / "boundary_ties.tsv", stage_a_dir / "shortlist_excluded.tsv",
        stage_a_dir / "shortlist_manifest.json",
    ]
    stage_a_target = 20
    stage_a_rescue = int(assay.get("shortlist_rescue_size", 5))
    stage_a_hash = digest([ranked_path], {"target_size": stage_a_target, "rescue_size": stage_a_rescue,
                                          "policy": "screening-shortlist-v1"})
    stage_a_ready = cached("stage_a", stage_a_hash, stage_a_outputs)
    if not stage_a_ready:
        try:
            screening_shortlist(str(ranked_path), str(stage_a_dir), target_size=stage_a_target,
                                rescue_size=stage_a_rescue)
        except Exception as error:  # noqa: BLE001 - one stage's failure must degrade to review, not kill the case
            review_rows.append({"reason": "stage_a_shortlist_failed", "target": str(case["case_id"]),
                                "required_next_action": f"the ranked table could not be cut into a frozen Stage-A shortlist ({type(error).__name__}: {error})",
                                "allowed_choices": "fix the ranked table and rerun | review the ranked table directly"})
        else:
            complete("stage_a", stage_a_hash, stage_a_outputs)
            stage_a_ready = True

    primer_dir = out / "primer"; primer_path = primer_dir / "primer_candidates.tsv"
    primer_results = [resolve(value) for value in inputs.get("primer_result_resources", [])]; primer_results = [p for p in primer_results if p]
    candidate_genomic = resolve(inputs.get("candidate_genomic_resource")); primerblast_root = resolve(assay.get("primerblast_root_resource")); nucleotide_db = resolve(target.get("nucleotide_db_resource"))
    primer_root_fingerprint = _primerblast_root_fingerprint(primerblast_root) if primerblast_root else None
    stage_a_shortlist_path = stage_a_dir / "deterministic_shortlist.tsv"
    primer_candidate_input = stage_a_shortlist_path if stage_a_ready else ranked_path
    primer_hash = digest(primer_results + [ranked_path, primer_candidate_input] + [p for p in (candidate_genomic, fasta, nucleotide_db) if p],
                         {"assay": assay, "candidate_scope": "frozen_stage_a_20", "primerblast_root": primer_root_fingerprint})
    primer_outputs = [primer_path]
    if not cached("primer", primer_hash, [primer_path]):
        primer_dir.mkdir(parents=True, exist_ok=True)
        if primer_results:
            primerblast_import([str(p) for p in primer_results], str(primer_dir), str(ranked_path))
        else:
            pd.DataFrame(columns=["gene_id","pair_specificity","product_size","coverage_status","validation_status","wet_validation_status","experimental_status","selection_reason"]).to_csv(primer_path, sep="\t", index=False, lineterminator="\n")
            if primerblast_root and fasta and nucleotide_db and candidate_genomic:
                if not primer_root_fingerprint["entrypoint_present"]:
                    raise ValueError("primerblast_root_resource lacks primerblast_oss/__main__.py for assay planning")
                plan = primer_plan("direct-sequencing", str(primer_dir / "plan"), str(primerblast_root), str(fasta), [str(nucleotide_db)],
                                   candidate_arg=str(primer_candidate_input), top_n=int(assay.get("top_n", 3)), gff_arg=str(gff3) if gff3 else None,
                                   candidate_genomic_fasta_arg=str(candidate_genomic), num_threads=int(assay.get("threads", 4)))
                primer_outputs.extend(Path(value) for key, value in plan.items() if key in {"script", "manifest"} and value)
        complete("primer", primer_hash, primer_outputs, {"primerblast_root_fingerprint": primer_root_fingerprint})

    wet_dir = out / "wet_state"
    wet_candidate_state = None
    wet_assay_registry = resolve(inputs.get("wet_assay_registry_resource"))
    wet_coverage = resolve(inputs.get("wet_coverage_resource"))
    wet_results = resolve(inputs.get("wet_results_resource"))
    wet_history = resolve(inputs.get("wet_history_resource"))
    exon_coverage = resolve(inputs.get("exon_coverage_resource"))
    wet_requested = bool(primer_results or exon_coverage or wet_assay_registry or wet_coverage or wet_results or wet_history)
    wet_outputs = [
        wet_dir / "assay_registry.tsv", wet_dir / "assay_target_coverage.tsv", wet_dir / "wet_results.tsv",
        wet_dir / "candidate_wet_state.tsv", wet_dir / "wet_state_history.tsv",
        wet_dir / "wet_review_required.tsv", wet_dir / "wet_state_manifest.json",
    ]
    if wet_requested:
        wet_inputs = [path for path in (primer_path, exon_coverage, wet_assay_registry, wet_coverage, wet_results, wet_history) if path]
        wet_hash = digest(wet_inputs, {"minimum_coverage_fraction": 1.0, "required_read_direction": "bidirectional"})
        if not cached("wet_state", wet_hash, wet_outputs):
            update_wet_state(
                str(wet_dir),
                assay_registry_arg=str(wet_assay_registry) if wet_assay_registry else None,
                coverage_arg=str(wet_coverage) if wet_coverage else None,
                results_arg=str(wet_results) if wet_results else None,
                previous_history_arg=str(wet_history) if wet_history else None,
                primer_results_arg=str(primer_path) if primer_results else None,
                exon_coverage_arg=str(exon_coverage) if exon_coverage else None,
            )
            complete("wet_state", wet_hash, wet_outputs)
        wet_candidate_state = wet_dir / "candidate_wet_state.tsv"

    hypothesis_dir = out / "hypotheses"
    hypothesis_ledger = None
    hypothesis_proposals = resolve(inputs.get("hypotheses_resource"))
    hypothesis_history = resolve(inputs.get("hypothesis_history_resource"))
    local_hypothesis_history = hypothesis_dir / "hypothesis_history.tsv"
    effective_hypothesis_history = hypothesis_history or (local_hypothesis_history if local_hypothesis_history.exists() else None)
    canonical_observations = resolve(inputs.get("canonical_observations_resource"))
    sv_events = resolve(inputs.get("sv_events_resource"))

    # --- QTL-effect context contract: optional, audit-only ---
    # The contract records what the mapped QTL says about allele effect so a
    # hypothesis can be checked for sign concordance. It is never forwarded to
    # prioritize_evidence or the rank/evidence digests (rank_contribution=none);
    # a declared-but-mismatched contract (top-level/nested case_id, header-only
    # TSV with a requested case, missing resource) fails the run closed.
    qtl_effect_contract_path = resolve(inputs.get("qtl_effect_contract_resource"))
    qtl_effect_contract = load_qtl_effect_contract(qtl_effect_contract_path, case_id=str(case["case_id"]))
    qtl_contract_fingerprint = file_fingerprint(qtl_effect_contract_path) if qtl_effect_contract_path else None
    qtl_contract_artifact = None
    if qtl_effect_contract_path is not None:
        qtl_contract_artifact = out / "audit" / "qtl_effect_contract.tsv"
        write_stable(qtl_contract_artifact, contract_frame(qtl_effect_contract).to_csv(sep="\t", index=False, lineterminator="\n"))
    hypothesis_audit_path = None
    qtl_contract_block: dict[str, object] = {
        "status": "used" if qtl_effect_contract_path is not None else "missing",
        "artifact": str(qtl_contract_artifact) if qtl_contract_artifact is not None else None,
        "fingerprint": (qtl_contract_fingerprint or {}).get("sha256"),
        "sign_audit": None,
        "sign_concordance_counts": {"concordant": 0, "incompatible": 0, "unknown": 0},
        "sign_concordance_review_required": False,
        "rank_contribution": "none",
    }

    # --- Promoter/UTR review route -> promoter_or_utr hypotheses (issue #230.2) ---
    # promoter_review() is a fully built, non-promoting review route (never touches
    # a gene score) that was previously only reachable as the standalone `promoter-
    # review` CLI subcommand. It runs here once its inputs are available and never
    # aborts the deterministic pipeline on failure; it degrades to a recorded review
    # item instead, the same way annotation-rescue and the private-variant funnel do.
    promoter_review_dir = out / "promoter_review"
    promoter_review_hypotheses_table = None
    promoter_config = case.get("promoter_review") or {}
    promoter_parental_variants = resolve(inputs.get("promoter_review_parental_variants_resource"))
    promoter_sv_calls = resolve(inputs.get("promoter_review_sv_calls_resource"))
    promoter_motif_hits = resolve(inputs.get("promoter_review_motif_hits_resource"))
    if not (promoter_parental_variants or promoter_sv_calls):
        review_rows.append({"reason": "promoter_review_not_requested", "target": str(case["case_id"]),
                            "required_next_action": "register inputs.promoter_review_parental_variants_resource and/or "
                                                    "inputs.promoter_review_sv_calls_resource to run the promoter/UTR review route",
                            "allowed_choices": "register review inputs | keep coding-only hypotheses"})
    elif not gff3:
        review_rows.append({"reason": "promoter_review_requires_gff3", "target": str(case["case_id"]),
                            "required_next_action": "register target.gff3_resource to run promoter/UTR review",
                            "allowed_choices": "register gff3 | skip promoter review"})
    else:
        promoter_manifest_path = promoter_review_dir / "promoter_review_manifest.json"
        promoter_inputs = [p for p in (flagged_candidates_path, gff3, promoter_parental_variants, promoter_sv_calls, promoter_motif_hits, exon_coverage) if p]
        promoter_hash = digest(promoter_inputs, {"windows": promoter_config.get("promoter_window"), "force": bool(promoter_config.get("force", False))})
        if cached("promoter_review", promoter_hash, [promoter_manifest_path]):
            candidate_table = promoter_review_dir / "promoter_or_utr_hypotheses.tsv"
            promoter_review_hypotheses_table = candidate_table if candidate_table.exists() else None
        else:
            try:
                promoter_result = promoter_review(
                    str(flagged_candidates_path), str(gff3), str(promoter_review_dir),
                    parental_variants_arg=str(promoter_parental_variants) if promoter_parental_variants else None,
                    sv_calls_arg=str(promoter_sv_calls) if promoter_sv_calls else None,
                    coverage_arg=str(exon_coverage) if exon_coverage else None,
                    motif_hits_arg=str(promoter_motif_hits) if promoter_motif_hits else None,
                    promoter_windows=[int(w) for w in (promoter_config.get("promoter_window") or [])] or None,
                    case_id=str(case["case_id"]), force=bool(promoter_config.get("force", False)),
                )
            except (ValueError, OSError) as error:
                review_rows.append({"reason": "promoter_review_failed", "target": str(promoter_review_dir),
                                    "required_next_action": f"promoter/UTR review did not complete; inspect inputs: {error}",
                                    "allowed_choices": "fix inputs and rerun | proceed with coding-only hypotheses"})
            else:
                complete("promoter_review", promoter_hash, [promoter_manifest_path])
                if not promoter_result.get("started"):
                    review_rows.append({"reason": "promoter_review_not_started", "target": str(promoter_review_dir),
                                        "required_next_action": str((promoter_result.get("route_gate") or {}).get(
                                            "note", "coding hypothesis not yet examined")),
                                        "allowed_choices": "provide coding coverage | pass promoter_review.force: true"})
                elif promoter_result.get("observation_count"):
                    observations = pd.read_csv(promoter_review_dir / "review_observations.tsv", sep="\t", dtype=str).fillna("")
                    promoter_hyp_rows = _promoter_observations_to_hypotheses(observations, str(case["case_id"]))
                    if len(promoter_hyp_rows):
                        promoter_hyp_path = promoter_review_dir / "promoter_or_utr_hypotheses.tsv"
                        promoter_hyp_rows.to_csv(promoter_hyp_path, sep="\t", index=False, lineterminator="\n")
                        promoter_review_hypotheses_table = promoter_hyp_path

    # --- Assembly SV/PAV plan -> sv_cnv_pav hypotheses via sv_events (issue #230.2) ---
    sv_plan_dir = out / "sv_plan"
    sv_plan_events_table = None
    sv_plan_config = case.get("sv_plan") or {}
    sv_plan_query_fasta = resolve(inputs.get("sv_plan_query_fasta_resource"))
    sv_plan_syri = resolve(inputs.get("sv_plan_syri_resource"))
    sv_plan_call_tables = {str(name): resolve(value) for name, value in (inputs.get("sv_plan_call_table_resources") or {}).items()}
    if sv_plan_query_fasta and fasta and len(intervals) == 1 and (sv_plan_syri or sv_plan_call_tables):
        sv_plan_outputs = [sv_plan_dir / "canonical_sv.tsv", sv_plan_dir / "sv_concordance.tsv", sv_plan_dir / "sv_evidence.tsv", sv_plan_dir / "sv_manifest.json"]
        sv_plan_inputs = [p for p in ([fasta, sv_plan_query_fasta, sv_plan_syri] + sorted(sv_plan_call_tables.values()) + [flagged_candidates_path]) if p]
        sv_plan_hash = digest(sv_plan_inputs, {"interval": intervals[0], "breakpoint_tolerance": int(sv_plan_config.get("breakpoint_tolerance", 100))})
        if cached("sv_plan", sv_plan_hash, sv_plan_outputs):
            sv_plan_events_table = sv_plan_dir / "canonical_sv.tsv"
        else:
            try:
                sv_plan(
                    str(fasta), str(sv_plan_query_fasta), intervals[0], str(sv_plan_dir),
                    candidates_arg=str(flagged_candidates_path),
                    syri_arg=str(sv_plan_syri) if sv_plan_syri else None,
                    call_specs=[f"{name}={path}" for name, path in sorted(sv_plan_call_tables.items())] or None,
                    breakpoint_tolerance=int(sv_plan_config.get("breakpoint_tolerance", 100)),
                )
            except (ValueError, OSError) as error:
                review_rows.append({"reason": "sv_plan_failed", "target": str(sv_plan_dir),
                                    "required_next_action": f"assembly SV/PAV plan did not complete; inspect reference/query/call-table inputs: {error}",
                                    "allowed_choices": "fix inputs and rerun | proceed without SV/PAV hypotheses"})
            else:
                complete("sv_plan", sv_plan_hash, sv_plan_outputs)
                sv_plan_events_table = sv_plan_dir / "canonical_sv.tsv"
    elif sv_plan_query_fasta or sv_plan_syri or sv_plan_call_tables:
        review_rows.append({"reason": "sv_plan_inputs_incomplete", "target": str(case["case_id"]),
                            "required_next_action": "sv-plan requires target.fasta_resource, inputs.sv_plan_query_fasta_resource, "
                                                    "a single validated interval, and a SyRI or call-table resource",
                            "allowed_choices": "register missing inputs | skip SV/PAV structural review"})

    # --- Multigenome-private events -> sv_events (issue #230.3) ---
    # focal_event_sharing_ledger.tsv rows are per-focal-event presence/absence
    # calls (a PAV signal), not the assembly-alignment SV calls sv-plan produces,
    # so the shapes differ (event_key vs call_id, gene_ids plural vs gene_id). This
    # adapter reprojects private_confirmed / private_possible rows onto sv-plan's
    # own call_id/gene_id contract instead of inventing a parallel schema, so both
    # feed the same sv_events_arg pathway below.
    multigenome_sv_events_table = None
    if private_status == "evaluated":
        sharing_path = private_dir / "focal_event_sharing_ledger.tsv"
        if sharing_path.exists():
            sharing = pd.read_csv(sharing_path, sep="\t", dtype=str).fillna("")
            candidate_events = sharing[sharing.get("event_class", pd.Series("", index=sharing.index)).isin({"private_confirmed", "private_possible"})]
            adapter_rows = []
            for record in candidate_events.to_dict(orient="records"):
                event_key = str(record.get("event_key", ""))
                genes = sorted({token.strip() for token in str(record.get("gene_ids", "")).replace("|", ";").split(";") if token.strip()})
                for gene_id in genes:
                    adapter_rows.append({"call_id": event_key, "gene_id": gene_id, "svtype": record.get("event_class", "")})
            if adapter_rows:
                multigenome_sv_events_table = private_dir / "private_events_as_sv_events.tsv"
                pd.DataFrame(adapter_rows, columns=["call_id", "gene_id", "svtype"]).drop_duplicates().sort_values(
                    ["gene_id", "call_id"], kind="stable"
                ).to_csv(multigenome_sv_events_table, sep="\t", index=False, lineterminator="\n")

    combined_sv_sources = [p for p in (sv_events, sv_plan_events_table, multigenome_sv_events_table) if p]
    combined_sv_events = None
    if combined_sv_sources:
        frames = []
        for path in combined_sv_sources:
            frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
            for column in ("call_id", "gene_id"):
                if column not in frame.columns:
                    frame[column] = ""
            frames.append(frame)
        merged_sv = pd.concat(frames, ignore_index=True, sort=False).fillna("")
        merged_sv = merged_sv[merged_sv["call_id"].str.strip().ne("")].drop_duplicates(subset=["call_id"]).sort_values(["gene_id", "call_id"], kind="stable")
        hypothesis_dir.mkdir(parents=True, exist_ok=True)
        combined_sv_events = hypothesis_dir / "sv_events_combined.tsv"
        merged_sv.to_csv(combined_sv_events, sep="\t", index=False, lineterminator="\n")

    # --- A/B dual-route comparison -> linked_qtl / unknown_locus_mechanism (issue #230.4) ---
    # dual_route_compare() is deterministic (no LLM/network call): given two
    # already-generated route hypothesis tables it only compares them. When a
    # precomputed comparison resource is supplied it is used as-is; otherwise both
    # route tables are required to compute one here.
    def resolve_bridge_override(name: str) -> Path | None:
        value = (narrowing_overrides or {}).get(name)
        if not value:
            return None
        path = runtime_path(str(value))
        if not path.exists():
            raise ValueError(f"guided-profile discordance override does not exist: {path}")
        return path

    dual_route_dir = out / "dual_route"
    dual_route_comparison_table = (
        resolve_bridge_override("discordance_comparison")
        or resolve(inputs.get("dual_route_comparison_resource"))
    )
    dual_route_a = resolve(inputs.get("dual_route_route_a_resource"))
    dual_route_b = resolve(inputs.get("dual_route_route_b_resource"))
    dual_route_a_manifest = (
        resolve_bridge_override("route_a_manifest")
        or resolve(inputs.get("dual_route_route_a_manifest_resource"))
    )
    dual_route_b_manifest = (
        resolve_bridge_override("route_b_manifest")
        or resolve(inputs.get("dual_route_route_b_manifest_resource"))
    )
    if not dual_route_comparison_table and dual_route_a and dual_route_b:
        dual_route_outputs = [
            dual_route_dir / "dual_route_hypothesis_comparison.tsv", dual_route_dir / "dual_route_review_required.tsv",
            dual_route_dir / "dual_route_convergence.json",
        ]
        dual_route_inputs = [p for p in (dual_route_a, dual_route_b, dual_route_a_manifest,
                                         dual_route_b_manifest) if p]
        dual_route_hash = digest(dual_route_inputs, {"case_id": case["case_id"]})
        if cached("dual_route_compare", dual_route_hash, dual_route_outputs):
            dual_route_comparison_table = dual_route_dir / "dual_route_hypothesis_comparison.tsv"
        else:
            try:
                dual_route_compare(
                    str(dual_route_a), str(dual_route_b), str(dual_route_dir), case_id=str(case["case_id"]),
                    route_b_manifest_arg=str(dual_route_b_manifest) if dual_route_b_manifest else None,
                )
            except (ValueError, OSError) as error:
                review_rows.append({"reason": "dual_route_compare_failed", "target": str(dual_route_dir),
                                    "required_next_action": f"A/B dual-route comparison did not complete; inspect route hypothesis tables: {error}",
                                    "allowed_choices": "fix route inputs and rerun | proceed without A/B escalation"})
            else:
                complete("dual_route_compare", dual_route_hash, dual_route_outputs)
                dual_route_comparison_table = dual_route_dir / "dual_route_hypothesis_comparison.tsv"
    elif not dual_route_comparison_table and (dual_route_a or dual_route_b):
        review_rows.append({"reason": "dual_route_compare_inputs_incomplete", "target": str(case["case_id"]),
                            "required_next_action": "A/B escalation needs both inputs.dual_route_route_a_resource and "
                                                    "inputs.dual_route_route_b_resource, or a precomputed inputs.dual_route_comparison_resource",
                            "allowed_choices": "register both route tables | register precomputed comparison | skip A/B escalation"})

    hypothesis_requested = bool(
        hypothesis_proposals or wet_candidate_state or hypothesis_history or annotation_rescue_table
        or combined_sv_events or dual_route_comparison_table or promoter_review_hypotheses_table
    )
    hypothesis_outputs = [
        hypothesis_dir / "locus_hypotheses.tsv", hypothesis_dir / "hypothesis_evidence_links.tsv",
        hypothesis_dir / "hypothesis_experiment_links.tsv", hypothesis_dir / "hypothesis_history.tsv",
        hypothesis_dir / "hypothesis_review_required.tsv", hypothesis_dir / "hypothesis_manifest.json",
    ]
    if hypothesis_requested:
        hypotheses_for_ledger = hypothesis_proposals
        if promoter_review_hypotheses_table is not None:
            promoter_rows = pd.read_csv(promoter_review_hypotheses_table, sep="\t", dtype=str).fillna("")
            if hypotheses_for_ledger:
                existing_rows = pd.read_csv(hypotheses_for_ledger, sep="\t", dtype=str).fillna("")
                for column in schemas.LOCUS_HYPOTHESIS_COLUMNS:
                    if column not in existing_rows.columns:
                        existing_rows[column] = ""
                combined_hypotheses = pd.concat([existing_rows[schemas.LOCUS_HYPOTHESIS_COLUMNS], promoter_rows], ignore_index=True)
            else:
                combined_hypotheses = promoter_rows
            hypothesis_dir.mkdir(parents=True, exist_ok=True)
            hypotheses_for_ledger = hypothesis_dir / "workflow_hypothesis_proposals.tsv"
            combined_hypotheses.to_csv(hypotheses_for_ledger, sep="\t", index=False, lineterminator="\n")
        hypothesis_inputs = [path for path in (
            hypotheses_for_ledger, hypothesis_history, canonical_observations, combined_sv_events,
            wet_candidate_state, wet_dir / "wet_results.tsv" if wet_candidate_state else None, ranked_path,
            annotation_rescue_table, dual_route_comparison_table, promoter_review_hypotheses_table,
        ) if path]
        hypothesis_config: dict[str, object] = {"case_id": case["case_id"], "policy": "locus-hypothesis-v2"}
        if qtl_effect_contract_path is not None and qtl_contract_fingerprint is not None:
            # A contract edit must invalidate the cached ledger/audit.
            hypothesis_config["qtl_effect_contract_sha256"] = qtl_contract_fingerprint.get("sha256")
        hypothesis_hash = digest(hypothesis_inputs, hypothesis_config)
        if not cached("hypothesis_ledger", hypothesis_hash, hypothesis_outputs):
            update_hypothesis_ledger(
                str(hypothesis_dir), hypotheses_arg=str(hypotheses_for_ledger) if hypotheses_for_ledger else None,
                case_id=str(case["case_id"]), evidence_arg=str(canonical_observations) if canonical_observations else None,
                wet_results_arg=str(wet_dir / "wet_results.tsv") if wet_candidate_state else None,
                sv_events_arg=str(combined_sv_events) if combined_sv_events else None,
                wet_candidate_state_arg=str(wet_candidate_state) if wet_candidate_state else None,
                annotation_rescue_arg=str(annotation_rescue_table) if annotation_rescue_table else None,
                dual_route_comparison_arg=str(dual_route_comparison_table) if dual_route_comparison_table else None,
                previous_history_arg=str(effective_hypothesis_history) if effective_hypothesis_history else None,
                ranked_arg=str(ranked_path),
            )
            complete("hypothesis_ledger", hypothesis_hash, hypothesis_outputs)
        hypothesis_ledger = hypothesis_dir / "locus_hypotheses.tsv"
        # Issue #360: sign-concordance audit on a copy of the ledger. The ledger
        # schema itself never carries expected_mechanism_sign, so it is merged
        # back from the proposal input; the audit never rewrites
        # locus_hypotheses.tsv and never touches the gene rank. `incompatible`
        # is review-visible only (no exclusion, no penalty, no reviewed marking).
        hypothesis_audit_path = hypothesis_dir / "locus_hypothesis_sign_audit.tsv"
        ledger_for_audit = pd.read_csv(hypothesis_ledger, sep="\t", dtype=str).fillna("")
        if hypotheses_for_ledger is not None:
            proposal_frame = pd.read_csv(hypotheses_for_ledger, sep="\t", dtype=str).fillna("")
            if "expected_mechanism_sign" in proposal_frame.columns:
                sign_rows = proposal_frame[["hypothesis_id", "expected_mechanism_sign"]].drop_duplicates(
                    subset=["hypothesis_id"], keep="first"
                )
                ledger_for_audit = ledger_for_audit.merge(sign_rows, on="hypothesis_id", how="left")
        annotated = annotate_sign_concordance(ledger_for_audit, qtl_effect_contract)
        annotated = annotated.sort_values("hypothesis_id", kind="stable")
        write_stable(hypothesis_audit_path, annotated.to_csv(sep="\t", index=False, lineterminator="\n"))
        counts = annotated["sign_concordance"].value_counts().to_dict()
        qtl_contract_block["sign_audit"] = str(hypothesis_audit_path)
        qtl_contract_block["sign_concordance_counts"] = {
            "concordant": int(counts.get("concordant", 0)),
            "incompatible": int(counts.get("incompatible", 0)),
            "unknown": int(counts.get("unknown", 0)),
        }
        qtl_contract_block["sign_concordance_review_required"] = (
            int(counts.get("incompatible", 0)) > 0
        )
        if state["steps"].get("hypothesis_ledger"):
            state["steps"]["hypothesis_ledger"]["qtl_effect_contract"] = qtl_contract_block

    # Issue #287: optional candidate-elimination stage between the hypothesis
    # ledger and next-experiment. It reuses the existing candidate_elimination()
    # and only runs when elimination evidence is registered; deterministic rank
    # is never changed by it. A failure degrades to review-required and a stale
    # remaining set is never forwarded.
    elimination_dir = out / "candidate_elimination"
    elimination_outputs = [
        elimination_dir / "candidate_survival.tsv", elimination_dir / "hypothesis_survival.tsv",
        elimination_dir / "elimination_events.tsv", elimination_dir / "remaining_candidate_set.tsv",
        elimination_dir / "elimination_review_required.tsv", elimination_dir / "elimination_manifest.json",
    ]
    elimination_evidence = resolve(inputs.get("elimination_evidence_resource"))
    elimination_history = resolve(inputs.get("elimination_history_resource"))
    elimination_rules = resolve(inputs.get("elimination_rules_resource"))
    remaining_candidates_path = None
    if elimination_evidence:
        elimination_inputs = [ranked_path, elimination_evidence] + [
            path for path in (hypothesis_ledger if hypothesis_ledger else None, wet_candidate_state,
                              elimination_history, elimination_rules) if path
        ]
        elimination_hash = digest(elimination_inputs, {"case_id": case["case_id"], "policy": "candidate-elimination-v1"})
        if not cached("candidate_elimination", elimination_hash, elimination_outputs):
            try:
                candidate_elimination(
                    candidates_arg=str(ranked_path), out_dir_arg=str(elimination_dir),
                    evidence_arg=str(elimination_evidence),
                    hypotheses_arg=str(hypothesis_ledger) if hypothesis_ledger else None,
                    wet_state_arg=str(wet_candidate_state) if wet_candidate_state else None,
                    rules_arg=str(elimination_rules) if elimination_rules else None,
                    previous_history_arg=str(elimination_history) if elimination_history else None,
                    case_id=str(case["case_id"]),
                )
            except Exception as error:  # noqa: BLE001 - one stage's failure must degrade to review, not kill the case
                review_rows.append({"reason": "candidate_elimination_failed", "target": str(case["case_id"]),
                                    "required_next_action": f"elimination evidence could not be applied ({type(error).__name__}: {error}); stale remaining sets are never reused",
                                    "allowed_choices": "fix elimination evidence and rerun | skip elimination"})
            else:
                complete("candidate_elimination", elimination_hash, elimination_outputs)
                remaining_candidates_path = elimination_dir / "remaining_candidate_set.tsv"
        else:
            remaining_candidates_path = elimination_dir / "remaining_candidate_set.tsv"

    # Issue #303: optional integrated-narrowing stage after candidate-elimination
    # and before next-experiment. It reuses integrated_narrowing() to compress
    # the ranked candidate set into an ordinal Wet-review order (N1..N4) from the
    # family-recovery / family-copy artifacts; it runs only when a family
    # artifact (family_to_copy_input.tsv or family_copy_discrimination.tsv) is
    # registered. The deterministic rank is never changed by it, every input is
    # part of the digest (stale-cache prevention), a failure degrades to
    # review-required, and stale narrowing outputs are never forwarded.
    # A triple-route guided profile can pass these artifacts directly through
    # narrowing_overrides when it also owns the case/resources manifests. The
    # rank is still produced here, once, under the normal case-run policy.
    narrowing_dir = out / "integrated_narrowing"
    narrowing_outputs = [
        narrowing_dir / "integrated_narrowing.tsv", narrowing_dir / "narrowing_manifest.json",
    ]
    narrowing_table = None
    narrowing_family_recovery = resolve_narrowing("family_recovery", "narrowing_family_recovery_resource")
    narrowing_convergence = resolve_narrowing("convergence", "narrowing_convergence_resource")
    narrowing_family_copy = resolve_narrowing("family_copy", "narrowing_family_copy_resource")
    narrowing_family_copy_shortlist = resolve_narrowing("family_copy_shortlist", "narrowing_family_copy_shortlist_resource")
    narrowing_family_copy_manifest = resolve_narrowing("family_copy_manifest", "narrowing_family_copy_manifest_resource")
    narrowing_remaining = resolve(inputs.get("narrowing_remaining_resource")) or remaining_candidates_path
    narrowing_provider_state = out / "narrowing_provider_state.json"
    bridge_identity_error = ""
    if narrowing_overrides:
        bridge_candidates_raw = str(narrowing_overrides.get("bridge_candidates", "")).strip()
        bridge_atlas_raw = str(narrowing_overrides.get("bridge_atlas", "")).strip()
        convergence_manifest_raw = str(narrowing_overrides.get("convergence_manifest", "")).strip()
        family_copy_manifest_raw = str(narrowing_overrides.get("family_copy_manifest", "")).strip()
        try:
            if not all((bridge_candidates_raw, bridge_atlas_raw, convergence_manifest_raw, family_copy_manifest_raw)):
                raise ValueError("triple-route bridge identity manifests/candidate universe are incomplete")
            bridge_candidates_path = runtime_path(bridge_candidates_raw)
            bridge_atlas_path = runtime_path(bridge_atlas_raw)
            convergence_manifest_path = runtime_path(convergence_manifest_raw)
            family_copy_manifest_path = runtime_path(family_copy_manifest_raw)
            narrowing_family_copy_manifest = family_copy_manifest_path
            convergence_manifest = json.loads(convergence_manifest_path.read_text(encoding="utf-8"))
            family_copy_manifest = json.loads(family_copy_manifest_path.read_text(encoding="utf-8"))
            expected_case_id = str(case["case_id"])
            for name, manifest in (("convergence", convergence_manifest), ("family-copy", family_copy_manifest)):
                if str(manifest.get("case_id", "")).strip() != expected_case_id:
                    raise ValueError(f"{name} manifest case_id does not match case-run {expected_case_id!r}")
            bridge_frame = pd.read_csv(bridge_candidates_path, sep="\t", dtype=str).fillna("")
            bridge_atlas = pd.read_csv(bridge_atlas_path, sep="\t", dtype=str).fillna("")
            if "gene_id" not in bridge_frame.columns or "gene_id" not in bridge_atlas.columns:
                raise ValueError("triple-route bridge candidates and atlas require gene_id")
            bridge_ids = {str(value).strip() for value in bridge_frame["gene_id"]}
            atlas_ids = {str(value).strip() for value in bridge_atlas["gene_id"]}
            ranked_ids = set(pd.read_csv(ranked_path, sep="\t", dtype=str).fillna("")["gene_id"].str.strip())
            if bridge_ids != ranked_ids or atlas_ids != ranked_ids:
                raise ValueError(
                    "triple-route candidate/atlas universe differs from case-run rank "
                    f"({len(bridge_ids)}/{len(atlas_ids)} != {len(ranked_ids)} genes)"
                )
            convergence_outputs = convergence_manifest.get("outputs") or {}
            family_outputs = family_copy_manifest.get("outputs") or {}
            bound_paths = (
                (convergence_outputs.get("family_to_copy_input"), narrowing_family_recovery),
                (convergence_outputs.get("triple_route_hypothesis_comparison"), narrowing_convergence),
                (family_outputs.get("family_copy_discrimination"), narrowing_family_copy),
                (family_outputs.get("family_copy_shortlist"), narrowing_family_copy_shortlist),
            )
            for recorded, actual in bound_paths:
                if not recorded or actual is None or runtime_path(str(recorded)).resolve() != actual.resolve():
                    raise ValueError("triple-route artifact is not the output recorded by its producing manifest")
            family_inputs = family_copy_manifest.get("inputs") or {}
            for name, actual in (("families", narrowing_family_recovery), ("atlas", bridge_atlas_path)):
                recorded = family_inputs.get(name) or {}
                actual_fp = file_fingerprint(actual)
                if not recorded.get("sha256") or recorded.get("sha256") != actual_fp.get("sha256"):
                    raise ValueError(f"family-copy manifest {name} fingerprint does not match the bridged artifact")
        except (OSError, ValueError, KeyError) as error:
            bridge_identity_error = str(error)
            review_rows.append({
                "reason": "triple_route_case_identity_mismatch",
                "target": str(case["case_id"]),
                "required_next_action": f"do not consume cross-case family artifacts: {error}",
                "allowed_choices": "rerun triple-route with matching case/candidate inputs | keep narrowing not_evaluated",
            })
            narrowing_family_recovery = None
            narrowing_convergence = None
            narrowing_family_copy = None
            narrowing_family_copy_shortlist = None
            narrowing_family_copy_manifest = None
    if narrowing_family_recovery or narrowing_family_copy:
        narrowing_inputs = [ranked_path] + [
            path for path in (narrowing_family_recovery, narrowing_convergence,
                              narrowing_family_copy, narrowing_family_copy_shortlist,
                              narrowing_remaining, narrowing_family_copy_manifest) if path
        ]
        narrowing_hash = digest(narrowing_inputs, {"case_id": case["case_id"],
                                                   "policy": "integrated-narrowing-v1"})
        if not cached("integrated_narrowing", narrowing_hash, narrowing_outputs):
            try:
                integrated_narrowing(
                    str(ranked_path), str(narrowing_dir),
                    family_recovery_arg=str(narrowing_family_recovery) if narrowing_family_recovery else None,
                    convergence_arg=str(narrowing_convergence) if narrowing_convergence else None,
                    family_copy_arg=str(narrowing_family_copy) if narrowing_family_copy else None,
                    family_copy_shortlist_arg=(str(narrowing_family_copy_shortlist)
                                               if narrowing_family_copy_shortlist else None),
                    remaining_arg=str(narrowing_remaining) if narrowing_remaining else None,
                    case_id=str(case["case_id"]),
                    family_copy_manifest_arg=(str(narrowing_family_copy_manifest)
                                              if narrowing_family_copy_manifest else None),
                )
            except Exception as error:  # noqa: BLE001 - one stage's failure must degrade to review, not kill the case
                review_rows.append({"reason": "integrated_narrowing_failed", "target": str(case["case_id"]),
                                    "required_next_action": f"family recovery / copy artifacts could not be narrowed ({type(error).__name__}: {error}); stale narrowing outputs are never reused",
                                    "allowed_choices": "fix family artifacts and rerun | skip narrowing"})
            else:
                complete("integrated_narrowing", narrowing_hash, narrowing_outputs)
                narrowing_table = narrowing_dir / "integrated_narrowing.tsv"
        else:
            narrowing_table = narrowing_dir / "integrated_narrowing.tsv"

    narrowing_state = {
        "schema_version": "1.0",
        "status": "evaluated" if narrowing_table else "not_evaluated",
        "reason": (
            "family_artifacts_available" if narrowing_table else
            "triple_route_case_identity_mismatch" if bridge_identity_error else
            "family_recovery_and_family_copy_inputs_missing_or_failed"
        ),
        "identity_error": bridge_identity_error or None,
        "rank_contribution": "none",
        "required_inputs": ["family_recovery", "family_copy"],
        "narrowing_table": str(narrowing_table) if narrowing_table else None,
    }
    write_stable(
        narrowing_provider_state,
        json.dumps(narrowing_state, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
    )

    stage_a_narrowed_outputs: list[Path] = []
    if narrowing_table:
        stage_a_narrowed_dir = out / "stage_a_narrowed"
        stage_a_narrowed_outputs = [
            stage_a_narrowed_dir / "deterministic_shortlist.tsv",
            stage_a_narrowed_dir / "rescue_appendix.tsv",
            stage_a_narrowed_dir / "boundary_ties.tsv",
            stage_a_narrowed_dir / "shortlist_excluded.tsv",
            stage_a_narrowed_dir / "shortlist_manifest.json",
        ]
        stage_a_narrowed_hash = digest(
            [narrowing_table],
            {"target_size": 20, "rescue_size": stage_a_rescue, "priority": "narrowing_priority"},
        )
        if not cached("stage_a_narrowed", stage_a_narrowed_hash, stage_a_narrowed_outputs):
            screening_shortlist(
                str(narrowing_table), str(stage_a_narrowed_dir), target_size=20,
                rescue_size=stage_a_rescue, priority_column="narrowing_priority",
                stage="A_narrowed_wet_priority",
            )
            complete("stage_a_narrowed", stage_a_narrowed_hash, stage_a_narrowed_outputs)

    next_dir = out / "next_experiment"
    next_options = next_dir / "next_experiment_options.tsv"
    next_decision_path = next_dir / "next_experiment_decision.json"
    next_dependencies = next_dir / "experiment_dependencies.tsv"
    next_review = next_dir / "experiment_review_required.tsv"
    mapping_gate = resolve(inputs.get("mapping_gate_resource"))
    candidate_decision = resolve(inputs.get("candidate_path_decision_resource"))
    supplied_flags = resolve(inputs.get("candidate_flags_resource"))
    assay_outcomes = resolve(inputs.get("assay_outcomes_resource"))
    exon_coverage = resolve(inputs.get("exon_coverage_resource"))
    informative_recombinants = resolve(inputs.get("informative_recombinants_resource"))
    next_inputs = [ranked_path, case_path, supplied_flags or flagged_candidates_path] + [
        path for path in (mapping_gate, candidate_decision, assay_outcomes, exon_coverage, informative_recombinants,
                          wet_candidate_state, hypothesis_ledger, remaining_candidates_path, narrowing_table) if path
    ]
    next_hash = digest(next_inputs, {"decision_support": case.get("decision_support", {})})
    if not cached("next_experiment", next_hash, [next_options, next_decision_path, next_dependencies, next_review]):
        select_next_experiment(
            str(ranked_path), str(case_path), str(next_dir),
            mapping_gate_arg=str(mapping_gate) if mapping_gate else None,
            candidate_path_decision_arg=str(candidate_decision) if candidate_decision else None,
            candidate_flags_arg=str(supplied_flags or flagged_candidates_path),
            assay_outcomes_arg=str(assay_outcomes) if assay_outcomes else None,
            exon_coverage_arg=str(exon_coverage) if exon_coverage else None,
            informative_recombinants_arg=str(informative_recombinants) if informative_recombinants else None,
            wet_candidate_state_arg=str(wet_candidate_state) if wet_candidate_state else None,
            hypotheses_arg=str(hypothesis_ledger) if hypothesis_ledger else None,
            remaining_candidates_arg=str(remaining_candidates_path) if remaining_candidates_path else None,
            narrowing_arg=str(narrowing_table) if narrowing_table else None,
        )
        complete("next_experiment", next_hash, [next_options, next_decision_path, next_dependencies, next_review])
    next_decision = json.loads(next_decision_path.read_text(encoding="utf-8"))

    # Issue #362 / PR #372: the locus causal scenario ledger is an independent
    # review artifact produced after next-experiment and before the report. It
    # never reads or writes rank/stage-A/score files (every row keeps
    # rank_contribution=none) and old cases without scenario inputs produce an
    # honest header-only ledger whose status stays unassessed. The #360
    # qtl_effect_contract module is NOT imported here (not merged yet); its
    # resource is read defensively and a broken contract only degrades the
    # recorded sign_contract_status, never the run.
    ledger_path = out / "locus_causal_scenarios.tsv"
    rows_input_status = "absent"
    scenario_rows_error = None
    scenario_rows_path = None
    try:
        scenario_rows_path = resolve(inputs.get("scenario_rows_resource"))
    except ValueError as error:
        scenario_rows_path = None
        rows_input_status = "invalid"
        scenario_rows_error = str(error)
    if scenario_rows_path:
        try:
            scenario_rows = load_scenario_rows(scenario_rows_path)
        except (ValueError, OSError) as error:
            scenario_rows = []
            rows_input_status = "invalid"
            scenario_rows_error = str(error)
        else:
            rows_input_status = "present"
    else:
        scenario_rows = []
    sign_contract_status = "absent"
    contract_path = None
    sign_concordance_value = None
    sign_evidence = None
    try:
        contract_path = resolve(inputs.get("qtl_effect_contract_resource"))
    except ValueError:
        contract_path = None
        sign_contract_status = "invalid"
    if contract_path:
        try:
            sign_concordance_value, sign_evidence = sign_concordance_evidence(
                load_scenario_effect_contract(contract_path), str(case["case_id"])
            )
        except (ValueError, OSError):
            sign_concordance_value = None
            sign_evidence = None
            sign_contract_status = "invalid"
        else:
            sign_contract_status = "present" if sign_concordance_value is not None else "unknown"
    if sign_evidence and sign_concordance_value in {"incompatible", "concordant"}:
        enriched_rows: list[dict[str, object]] = []
        for row in scenario_rows:
            enriched = dict(row)
            evidence_column = (
                "contradicting_evidence" if sign_concordance_value == "incompatible"
                else "supporting_evidence"
            )
            existing = str(enriched.get(evidence_column, "") or "")
            enriched[evidence_column] = f"{existing}; {sign_evidence}" if existing else sign_evidence
            sign_source = f"qtl_effect_contract:{case['case_id']}"
            source_ids = str(enriched.get("source_ids", "") or "")
            enriched["source_ids"] = f"{source_ids}; {sign_source}" if source_ids else sign_source
            enriched_rows.append(enriched)
        scenario_rows = enriched_rows
    # Issue #363: an A/B dual-route comparison artifact in scope is an entrance
    # to review-required causal scenarios, so its rows are mapped by
    # discordance_review_rows() and appended to the same ledger (the ledger is
    # rebuilt from all inputs every run, so resume never duplicates rows).
    # Route B's manifest (interval_first_status.json, JSON) can declare its
    # provider missing (a route status declaring not_evaluated): then every row
    # collapses to one provider_missing / not_evaluated row instead of mixing
    # into neither.  A legacy precomputed dual comparison without a status still
    # preserves its old behaviour, but every explicit Route A/B status artifact
    # is authoritative and read fail closed. A broken comparison/manifest only
    # records the status and warns; it never crashes the case and never touches
    # scenario_rows-derived rows.
    discordance_input_status = "absent"
    discordance_rows_count = 0
    route_a_status = "absent"
    route_b_status = "absent"
    discordance_error = None
    ledger_rows = list(scenario_rows)
    if dual_route_comparison_table:
        route_a_evaluated = True
        route_b_evaluated = True
        try:
            input_fingerprint = comparison_fingerprint(dual_route_comparison_table)
            if dual_route_a_manifest:
                route_a_evaluated, route_a_status = _read_route_manifest_status(dual_route_a_manifest, "a")
            elif dual_route_a:
                route_a_evaluated, route_a_status = _route_table_evaluable(dual_route_a)
            if dual_route_b_manifest:
                route_b_evaluated, route_b_status = _read_route_manifest_status(dual_route_b_manifest, "b")
            elif dual_route_b:
                route_b_evaluated, route_b_status = _route_table_evaluable(dual_route_b)
            ledger_rows.extend(discordance_review_rows(
                dual_route_comparison_table, case_id=str(case["case_id"]),
                input_fingerprint=input_fingerprint,
                route_a_evaluated=route_a_evaluated, route_b_evaluated=route_b_evaluated,
            ))
            discordance_rows_count = len(ledger_rows) - len(scenario_rows)
            discordance_input_status = "present"
        except (ValueError, OSError) as error:
            discordance_input_status = "invalid"
            discordance_error = str(error)
    scenario_hash = digest(
        [path for path in (
            scenario_rows_path, contract_path,
            dual_route_comparison_table if dual_route_comparison_table else None,
            dual_route_a_manifest if dual_route_a_manifest else None,
            dual_route_b_manifest if dual_route_b_manifest else None,
        ) if path],
        {"scenario_rows": scenario_rows_path is not None, "qtl_effect_contract": contract_path is not None,
         "discordance": dual_route_comparison_table is not None,
          "route_a_manifest": dual_route_a_manifest is not None,
          "route_b_manifest": dual_route_b_manifest is not None,
         "case_id": case["case_id"]},
    )
    if not cached("causal_scenarios", scenario_hash, [ledger_path]):
        # Invalid scenario rows (bad normalization, unknown scenario_type) must
        # degrade like the other optional lanes: write an honest header-only
        # ledger, record rows_input=invalid, and surface the error in the report
        # instead of crashing the case. The same holds for a broken discordance
        # input: any valid scenario rows still reach the ledger, the status is
        # recorded as invalid, and the report warns.
        try:
            write_scenario_ledger(ledger_rows, ledger_path, case_id=str(case["case_id"]))
        except (ValueError, OSError) as error:
            if discordance_input_status == "present":
                discordance_error = str(error)
            else:
                scenario_rows_error = str(error)
                rows_input_status = "invalid"
            pd.DataFrame(columns=SCENARIO_COLUMNS).to_csv(ledger_path, sep="\t", index=False, lineterminator="\n")
        written_ledger = pd.read_csv(ledger_path, sep="\t", dtype=str, keep_default_na=False)
        ledger_error = scenario_rows_error or discordance_error
        complete("causal_scenarios", scenario_hash, [ledger_path], {
            "scenario_count": int(len(written_ledger)),
            "sign_contract_status": sign_contract_status,
            "rows_input": rows_input_status,
            "discordance_rows": discordance_rows_count,
            "discordance_input": discordance_input_status,
            "route_a_status": route_a_status,
            "route_b_status": route_b_status,
            **({"error": ledger_error} if ledger_error else {}),
        })
    ledger_frame = pd.read_csv(ledger_path, sep="\t", dtype=str).fillna("")
    scenario_count = len(ledger_frame)

    llm_ranked_path = None
    if mode == "llm-assisted":
        llm_dir = out / "llm_assisted"; llm_dir.mkdir(parents=True, exist_ok=True)
        llm_resources = [resolve(value) for value in inputs.get("llm_evidence_resources", [])]; llm_resources = [p for p in llm_resources if p]
        reviewed_frames = []
        for path in llm_resources:
            frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
            if "reviewed" in frame.columns: frame = frame[frame.reviewed.str.casefold().isin({"true","yes","1"})]
            source_ids = frame["source_id"] if "source_id" in frame.columns else pd.Series("", index=frame.index)
            reviewed_frames.append(frame[source_ids.str.strip().ne("")])
        if reviewed_frames and sum(len(frame) for frame in reviewed_frames):
            deterministic = pd.read_csv(evidence_path, sep="\t", dtype=str).fillna("")
            assisted = pd.concat([deterministic, *reviewed_frames], ignore_index=True).drop_duplicates()
            assisted = assisted.sort_values(["gene_id","evidence_code","source_id","direction"], kind="stable")
            assisted_path = llm_dir / "evidence.tsv"; assisted.to_csv(assisted_path, sep="\t", index=False, lineterminator="\n")
            assisted_rank_dir = llm_dir / "rank"
            prioritize_evidence(str(flagged_candidates_path), str(assisted_path), str(assisted_rank_dir), None, int(assay.get("top_n", 3)), int(assay.get("reserve_n", 3)), float(case.get("min_top_score", 2.0)), int(assay.get("rescue_n", 3)))
            llm_ranked_path = assisted_rank_dir / "ranked_candidates.tsv"
        else:
            review_rows.append({"reason":"llm_assisted_evidence_missing","target":"llm-assisted rank","required_next_action":"provide human-reviewed source-backed canonical evidence","allowed_choices":"add reviewed evidence | use deterministic rank only"})

    review_path = out / "review_required.tsv"
    if excluded_path.exists() and len(pd.read_csv(excluded_path, sep="\t", dtype=str)):
        review_rows.append({"reason":"evidence_without_source_id","target":"excluded evidence rows","required_next_action":"add DOI/PMID/database identifier or keep excluded","allowed_choices":"add source ID | discard evidence"})
    if not primer_results:
        review_rows.append({"reason":"primer_results_missing","target":"Top candidates","required_next_action":"run generated PrimerBlast OSS plans and import JSON","allowed_choices":"run primer design | defer wet assay"})
    if not allowed: review_rows.append({"reason":"unresolved_position","target":str(source.get("source_id", "case")),"required_next_action":"provide validated interval or pair-level marker hits","allowed_choices":"provide interval | provide marker hits | retain rescue-only"})
    review_frame = pd.DataFrame(review_rows, columns=["reason","target","required_next_action","allowed_choices"]).drop_duplicates().sort_values(["reason","target"], kind="stable")
    write_stable(review_path, review_frame.to_csv(sep="\t", index=False, lineterminator="\n"))

    run_manifest_path = out / "run_manifest.json"
    resource_fingerprints = {}
    for resource_id, entry in sorted(resource_map.items()):
        try:
            path = resolve(resource_id)
            if path and path.is_file(): resource_fingerprints[resource_id] = file_fingerprint(path, full_hash=strict_fingerprints)
            elif path and resource_id == assay.get("primerblast_root_resource"):
                resource_fingerprints[resource_id] = primer_root_fingerprint
        except ValueError:
            resource_fingerprints[resource_id] = {"status":"missing", "entry":entry}
    case_fp = file_fingerprint(case_path, full_hash=strict_fingerprints)
    registry_fp = file_fingerprint(resources_path, full_hash=strict_fingerprints)
    all_fingerprints = {"case_manifest": case_fp, "resource_registry": registry_fp, **resource_fingerprints}
    weak_fingerprints = sorted(name for name, fp in all_fingerprints.items() if "path" in fp and fingerprint_is_weak(fp))
    fingerprint_warnings = []
    if strict_fingerprints:
        fingerprint_warnings.append("strict_fingerprints: all inputs hashed by full content")
    elif resume and weak_fingerprints:
        fingerprint_warnings.append("resume reused steps with metadata-only fingerprints for: " + ", ".join(weak_fingerprints))
    # Issue #364 adversarial fix: the top-level evidence_dependency block is
    # rebuilt from the dependency manifest ON DISK (both the evaluated and the
    # provider-failed paths write it), so a resumed re-run that failed after a
    # previously successful step record can never surface stale counts from
    # case_state.json's evidence step record.
    try:
        dependency_manifest_on_disk = json.loads(evidence_dependency_manifest_path.read_text(encoding="utf-8"))
        evidence_dependency_block = {
            key: dependency_manifest_on_disk[key]
            for key in ("dependency_status", "evidence_count",
                        "independent_evidence_cluster_count",
                        "dependency_unknown_evidence_count",
                        "dependency_unknown_present", "rank_contribution")
        }
        evidence_dependency_block.update({
            "input_fingerprint": dependency_manifest_on_disk.get("input_fingerprint", {}),
            "clusters": str(evidence_clusters_path.relative_to(out)),
            "manifest": str(evidence_dependency_manifest_path.relative_to(out)),
            "note": "reporting only; dependency clusters never contribute to the candidate rank",
        })
    except (OSError, ValueError, KeyError):
        evidence_dependency_block = {
            "dependency_status": "not_evaluated",
            "evidence_count": 0,
            "independent_evidence_cluster_count": 0,
            "dependency_unknown_evidence_count": 0,
            "dependency_unknown_present": False,
            "rank_contribution": "none",
            "note": "evidence dependency manifest missing; counts unavailable",
        }
    run_manifest = {"schema_version":"1.5", "case_manifest":case_fp, "resource_registry":registry_fp,
                    "resources":resource_fingerprints, "weak_fingerprints":weak_fingerprints, "fingerprint_warnings":fingerprint_warnings,
                    "fingerprint_mode":"strict_full_hash" if strict_fingerprints else "size_tiered",
                    "tools":{name:shutil.which(name) for name in ["blastn","blastp","primer3_core","minimap2","miniprot"]},
                    "steps":{key:value for key,value in state["steps"].items() if key != "report"}, "mode":mode, "llm_called":False,
                    "private_variant_status":private_status,
                    "private_variant_manifest":str(private_dir / "private_variant_manifest.json") if private_status == "evaluated" else None,
                    "private_variant_policy":"private status is a triage state only: it adds no score, promotes nothing, and a missing call is never a reference allele",
                    "deterministic_artifacts_preserved":True, "llm_assisted_rank":str(llm_ranked_path) if llm_ranked_path else None,
                    "rank_semantic_status":rank_semantic_status, "rank_is_positional_fallback":rank_is_positional_fallback,
                    "rank_metric_eligible":rank_metric_eligible, "stage_a_is_budget_only":stage_a_is_budget_only,
                    "rank_lane_warning":rank_lane_warning,
                    "scenario_ledger_status":{"rows_input":rows_input_status,
                                               "sign_contract_status":sign_contract_status,
                                               "discordance_input":discordance_input_status,
                                               "discordance_rows":discordance_rows_count,
                                               "route_a_status":route_a_status,
                                               "route_b_status":route_b_status,
                                              "scenario_count":scenario_count},
                    "evidence_dependency":evidence_dependency_block,
                    "qtl_effect_contract":qtl_contract_block}
    # Issue #397: the penalty-only declaration is recorded only when true so an
    # evidence-free lane keeps its pre-397 manifest byte-identical; the semantic
    # status alone still distinguishes every state for non-penalty runs.
    if rank_is_penalty_only:
        run_manifest["rank_is_penalty_only"] = True
    write_stable(run_manifest_path, json.dumps(run_manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n")
    report = out / "report.md"; ranked = pd.read_csv(ranked_path, sep="\t", dtype=str).fillna(""); primers = pd.read_csv(primer_path, sep="\t", dtype=str).fillna("")
    report_lines = [f"# {case['case_id']}", "", "## Provenance", "", f"- Route: `{coordinate['coordinate_method']}`", f"- Position status: `{coordinate['position_status']}`", f"- Source: `{coordinate.get('source_id','')}`", f"- Target intervals: `{'; '.join(intervals)}`", "- Deterministic rank is always preserved; LLM evidence never overwrites it.", "", "## Candidate groups", ""]
    counts = ranked["dry_group"].value_counts().to_dict() if "dry_group" in ranked else {}
    report_lines += [f"- {key}: {counts[key]}" for key in sorted(counts)]
    if rank_is_positional_fallback:
        report_lines += ["", f"> {rank_lane_warning}"]
    report_lines += ["", "## Ranked candidates", "", "| gene_id | group | score | evidence trace |", "|---|---:|---:|---|"]
    for row in ranked.head(20).itertuples():
        report_lines.append(f"| {row.gene_id} | {getattr(row, 'dry_group', '')} | {getattr(row, 'dry_score', '')} | {getattr(row, 'score_breakdown', '')} |")
    report_lines += ["", f"Full source-level trace: `{contributions_path.relative_to(out)}` and `{evidence_path.relative_to(out)}`."]
    report_lines += ["", "## Evidence dependency", ""]
    if evidence_provider_failed:
        report_lines += ["- Status: not_evaluated (evidence provider failed)"]
    else:
        dependency_counts = evidence_step_extra
        report_lines += [
            f"- Evidence rows: {dependency_counts.get('evidence_count', len(included_evidence))}",
            f"- Known independent dependency clusters: {dependency_counts.get('independent_evidence_cluster_count', 0)}",
            f"- Dependency-unknown rows: {dependency_counts.get('dependency_unknown_evidence_count', 0)} ({'present' if dependency_counts.get('dependency_unknown_present', False) else 'absent'})",
            "- The dependency-unknown bucket is NOT counted as an independent evidence.",
            "- `rank_contribution=none`: these counts never change the candidate rank.",
            "",
            "| dependency_cluster | dependency_status | root_source_id | evidence_count | evidence_kinds | routes | source_ids |",
            "|---|---|---|---|---|---|---|",
        ]
        dependency_clusters = pd.read_csv(evidence_clusters_path, sep="\t", dtype=str).fillna("")
        for record in dependency_clusters.itertuples():
            report_lines.append(
                f"| {record.dependency_cluster} | {record.dependency_status} | {record.root_source_id} | "
                f"{record.evidence_count} | {record.evidence_kinds} | {record.routes} | {record.source_ids} |"
            )
    if annotation_rescue_table:
        rescue_models = pd.read_csv(annotation_rescue_table, sep="\t", dtype=str).fillna("")
        rescue_counts = rescue_models["model_class"].value_counts().to_dict()
        report_lines += ["", "## Annotation rescue", "", f"- Projection audit rows: {len(rescue_models)} ({rescue_counts}).", "- This bounded appendix never changes the deterministic rank or fixed primary 20.", f"- Review: `{(annotation_rescue_dir / 'annotation_rescue_review.tsv').relative_to(out)}`."]
    report_lines += ["", "## Private variants (multigenome)", ""]
    if private_status == "evaluated":
        private_summary = pd.read_csv(private_dir / "event_class_summary.tsv", sep="\t", dtype=str).fillna("")
        private_counts = private_summary.groupby("event_class")["event_count"].apply(lambda values: sum(int(v) for v in values)).to_dict()
        denominators = pd.read_csv(private_dir / "projection_and_callable_denominators.tsv", sep="\t", dtype=str).fillna("")
        report_lines += [f"- Event classes: {private_counts}"]
        for row in denominators.itertuples():
            report_lines.append(f"- `{row.locus_id}`: {row.high_confidence_one_to_one}/{row.reference_genes_in_interval} high-confidence 1:1 projection, callable CDS fraction {row.callable_cds_fraction} (`{row.denominator_status}`).")
        report_lines += ["- `private_possible` means a required control was uncallable; it never removes a candidate.",
                         "- Private status adds no score and promotes nothing; see `private_variants/working_candidate_panel.tsv`."]
    else:
        report_lines += [f"- Status: `{private_status}`.",
                         "- Private/shared status was not evaluated, so a low private-variant count is NOT biological evidence here."]
    report_lines += ["", "## Assays", "", f"- Imported primer pairs: {len(primers)}", "- `in_silico` predictions and wet validation are separate columns."]
    if wet_candidate_state:
        wet_states = pd.read_csv(wet_candidate_state, sep="\t", dtype=str).fillna("")
        wet_reviews = pd.read_csv(wet_dir / "wet_review_required.tsv", sep="\t", dtype=str).fillna("")
        report_lines += [f"- Normalized Wet candidate states: {len(wet_states)}", f"- Wet review items: {len(wet_reviews)} (`{(wet_dir / 'wet_review_required.tsv').relative_to(out)}`).", "- Assay failure, partial coverage, and reviewed negative evidence remain distinct."]
    if hypothesis_ledger:
        hypotheses = pd.read_csv(hypothesis_ledger, sep="\t", dtype=str).fillna("")
        non_gene_types = {"annotation_missing", "linked_qtl", "intergenic_regulatory", "unknown_locus_mechanism"}
        non_gene_count = int(hypotheses["hypothesis_type"].isin(non_gene_types).sum()) if len(hypotheses) else 0
        report_lines += ["", "## Locus hypotheses", "", f"- Mechanism/locus hypotheses: {len(hypotheses)} (non-gene locus mechanisms: {non_gene_count}).",
                         "- Hypothesis status is separate from deterministic gene rank; LLM proposals remain pending review.",
                         "- Annotation-rescue, promoter/UTR review, SV-plan, multigenome-private, and A/B dual-route disagreement "
                         "seed non-gene hypotheses here deterministically; none of them change the gene score."]
    report_lines += ["", "## Next experiment decision support", "", f"- Primary action: `{next_decision['primary_recommendation']['experiment_id']}`", f"- Feasibility: `{next_decision['primary_recommendation']['feasibility_status']}`", "- This is rule-based ordinal decision support, not Bayesian information gain, causal probability, or automatic Wet execution.", "- Human approval: `pending`."]
    scenario_type_count = int(ledger_frame["scenario_type"].nunique()) if scenario_count else 0
    state_breakdown = (
        ", ".join(f"{state}: {int(count)}" for state, count in
                  ledger_frame["assessment_state"].value_counts().sort_index().items())
        if scenario_count else "none"
    )
    report_lines += ["", render_scenario_report(ledger_frame).rstrip(), "",
                     f"- Ledger summary: {scenario_count} scenario(s) across {scenario_type_count} type(s); "
                     f"states: {state_breakdown}; discordance-derived: {discordance_rows_count} "
                     f"({discordance_input_status}); see `{ledger_path.relative_to(out)}`."]
    if scenario_rows_error:
        report_lines += [f"- Warning: scenario_rows input invalid: {scenario_rows_error}; ledger left unassessed."]
    if discordance_error:
        report_lines += [f"- Warning: discordance input invalid: {discordance_error}; "
                         "discordance rows not added to the ledger."]
    report_lines += ["", "## Review checkpoints", "", f"- See `{review_path.name}` ({len(review_rows)} item(s)).", f"- See `{next_review.relative_to(out)}` for next-experiment missing inputs and constraints.", ""]
    write_stable(report, "\n".join(report_lines))
    final_artifacts = [coordinate_path,candidates_path,flagged_candidates_path,evidence_path,evidence_clusters_path,evidence_dependency_manifest_path,ranked_path,contributions_path,*(stage_a_outputs if stage_a_ready else []),primer_path,next_options,next_decision_path,next_dependencies,next_review,review_path,run_manifest_path,narrowing_provider_state,ledger_path]
    if private_status == "evaluated":
        final_artifacts.extend(private_outputs)
    if wet_candidate_state:
        final_artifacts.extend(wet_outputs)
    if hypothesis_ledger:
        final_artifacts.extend(hypothesis_outputs)
    if annotation_rescue_table:
        final_artifacts.extend([
            annotation_rescue_dir / "annotation_rescue.tsv", annotation_rescue_dir / "annotation_rescue_review.tsv",
            annotation_rescue_dir / "annotation_rescue_review_overflow.tsv", annotation_rescue_dir / "annotation_rescue_manifest.json",
        ])
    if promoter_review_hypotheses_table:
        final_artifacts.append(promoter_review_hypotheses_table)
    if sv_plan_events_table:
        final_artifacts.append(sv_plan_events_table)
    if multigenome_sv_events_table:
        final_artifacts.append(multigenome_sv_events_table)
    if dual_route_comparison_table:
        final_artifacts.append(dual_route_comparison_table)
    if narrowing_table:
        final_artifacts.extend(narrowing_outputs + stage_a_narrowed_outputs)
    if qtl_contract_artifact is not None:
        final_artifacts.append(qtl_contract_artifact)
    if hypothesis_audit_path is not None:
        final_artifacts.append(hypothesis_audit_path)
    final_hash = digest(final_artifacts, {"case_id":case["case_id"]})
    if not cached("report", final_hash, [report,run_manifest_path]): complete("report", final_hash, [report,run_manifest_path])
    final = {"case_id":case["case_id"],"state":"REPORT_READY","route":coordinate["coordinate_method"],"positional_ranking_allowed":allowed,
             "candidate_scope":"positional" if allowed else "rescue","candidate_count":len(candidates),"evidence":str(evidence_path),"ranked":str(ranked_path),
             "score_contributions":str(contributions_path),"stage_a_shortlist":str(stage_a_shortlist_path) if stage_a_ready else None,"stage_a_narrowed_shortlist":str(out / "stage_a_narrowed" / "deterministic_shortlist.tsv") if narrowing_table else None,"narrowing_provider_state":str(narrowing_provider_state),"integrated_narrowing":str(narrowing_table) if narrowing_table else None,"private_variant_status":private_status,"private_variant_manifest":str(private_dir / "private_variant_manifest.json") if private_status == "evaluated" else None,"private_variant_panel":str(private_dir / "working_candidate_panel.tsv") if private_status == "evaluated" else None,"annotation_rescue":str(annotation_rescue_table) if annotation_rescue_table else None,"annotation_rescue_review":str(annotation_rescue_dir / "annotation_rescue_review.tsv") if annotation_rescue_table else None,"primer_status":str(primer_path),"wet_candidate_state":str(wet_candidate_state) if wet_candidate_state else None,"wet_review_required":str(wet_dir / "wet_review_required.tsv") if wet_candidate_state else None,"locus_hypotheses":str(hypothesis_ledger) if hypothesis_ledger else None,"hypothesis_review_required":str(hypothesis_dir / "hypothesis_review_required.tsv") if hypothesis_ledger else None,"promoter_review_hypotheses":str(promoter_review_hypotheses_table) if promoter_review_hypotheses_table else None,"sv_plan_events":str(sv_plan_events_table) if sv_plan_events_table else None,             "dual_route_comparison":str(dual_route_comparison_table) if dual_route_comparison_table else None,"next_experiment_decision":str(next_decision_path),"next_experiment_options":str(next_options),"review_required":str(review_path),"report":str(report),"run_manifest":str(run_manifest_path),"llm_assisted_ranked":str(llm_ranked_path) if llm_ranked_path else None,"scenario_ledger":str(ledger_path),"scenario_count":scenario_count,"steps":state["steps"]}
    final.update({
        "remaining_candidate_set": str(remaining_candidates_path) if remaining_candidates_path else None,
        "narrowing_manifest": str(narrowing_dir / "narrowing_manifest.json") if narrowing_table else None,
        "rank_semantic_status": rank_semantic_status,
        "rank_is_positional_fallback": rank_is_positional_fallback,
        "rank_metric_eligible": rank_metric_eligible,
        "stage_a_is_budget_only": stage_a_is_budget_only,
        "rank_lane_warning": rank_lane_warning,
        "qtl_effect_contract": qtl_contract_block,
    })
    # Issue #397: same conditional penalty-only declaration as run_manifest.json.
    if rank_is_penalty_only:
        final["rank_is_penalty_only"] = True
    state.update(final); state_path.write_text(json.dumps(state, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return final | {"resumed": resume, "fingerprint_warnings": fingerprint_warnings, "weak_fingerprints": weak_fingerprints}
