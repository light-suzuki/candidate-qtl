"""Callability-aware multigenome private-variant funnel (issue #164).

Classifies every normalized focal event against a declared set of control
accessions, keeping *absence of a call* and *inability to call* strictly apart.

Deliberate safety boundaries, mirroring the rest of the review routes:

* A missing VCF record is never converted into a callable reference allele, and a
  failed projection is never converted into gene absence.
* ``private_possible`` (a required control is uncallable) never removes a
  candidate; it becomes an additional-callability action item.
* Private status is not folded into any score. Every emitted row carries
  ``rank_contribution=none``; promotion stays a human decision.
* Counts from different callsets are never summed. Every focal event is
  reconciled per control, and the control matrix keeps the per-control state.
* Projection and callable-CDS denominators are always reported, so a small
  private-variant count cannot be read as negative biological evidence when it is
  really an annotation or coverage gap.
"""

from __future__ import annotations

import bisect
import json
from pathlib import Path

import pandas as pd
import yaml

from .core import attributes, file_fingerprint, open_text, parse_interval, runtime_path

EVENT_CLASSES = ("private_confirmed", "private_possible", "shared_confirmed", "not_evaluated")
CANDIDATE_GROUPS = ("co_primary", "reserve", "annotation_rescue", "regulatory_sv", "evidence_insufficient")
PROJECTION_STATES = ("high_confidence_one_to_one", "any_hit", "partial", "unmapped", "not_evaluated")

SHARING_COLUMNS = [
    "event_key", "chrom", "pos", "ref", "alt", "focal_accession", "locus_id", "phenotype_controls",
    "event_class", "shared_with", "absent_callable_in", "absent_uncallable_in", "not_evaluated_in",
    "gene_ids", "reference_region_class", "reference_consequence", "native_consequence",
    "consequence_agreement", "protein_altering", "rank_contribution", "review_status", "notes",
]
CONTROL_MATRIX_COLUMNS = [
    "event_key", "control_accession", "control_role", "genotype_state", "event_observed",
    "callability_state", "callability_source", "resolution", "notes",
]
CALLABILITY_COLUMNS = [
    "control_accession", "control_role", "callability_resource", "callable_intervals",
    "callable_bp_in_intervals", "events_callable", "events_uncallable", "events_without_callability",
    "private_confirmed_blocked_by_callability", "notes",
]
EVENT_SUMMARY_COLUMNS = ["locus_id", "event_class", "event_count", "gene_count", "protein_altering_events", "notes"]
GENE_LEDGER_COLUMNS = [
    "gene_id", "locus_id", "chrom", "start", "end", "projection_state", "protein_effect_evaluated",
    "private_confirmed", "private_possible", "shared_confirmed", "not_evaluated",
    "protein_altering_private_confirmed", "regulatory_private_events", "consequence_disagreements",
    "rank_contribution", "review_status", "notes",
]
DENOMINATOR_COLUMNS = [
    "locus_id", "interval", "reference_genes_in_interval", "any_hit_projection", "high_confidence_one_to_one",
    "partial_projection", "unmapped", "projection_not_evaluated", "native_annotation_genes",
    "intended_cds_bp", "callable_cds_bp", "callable_cds_fraction", "protein_effect_evaluated_genes",
    "denominator_status", "required_next_action",
]
PANEL_COLUMNS = [
    "gene_id", "locus_id", "candidate_group", "private_confirmed", "private_possible",
    "protein_altering_private_confirmed", "projection_state", "rank_contribution", "gate",
    "required_next_action", "notes",
]
REVIEW_COLUMNS = ["reason", "target", "required_next_action", "allowed_choices"]

_MIN_CALLABLE_CDS_FRACTION = 0.8
_MIN_PROJECTION_FRACTION = 0.8


# --------------------------------------------------------------------------- #
# normalized event keys
# --------------------------------------------------------------------------- #

def normalized_event_key(chrom: str, pos: int, ref: str, alt: str) -> tuple[str, int, str, str]:
    """Minimal (parsimonious) representation of one biallelic REF/ALT pair.

    Trims the shared suffix and then the shared prefix so that equivalent indel
    spellings match across callsets. This does NOT left-shift over a repeat -
    that needs the reference FASTA, so inputs must already be
    ``bcftools norm -f <ref> -m -any`` normalized; the manifest records whether
    each input declared it.
    """
    ref, alt = ref.upper(), alt.upper()
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    trimmed = 0
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, trimmed = ref[1:], alt[1:], trimmed + 1
    return str(chrom), int(pos) + trimmed, ref, alt


def _event_id(key: tuple[str, int, str, str]) -> str:
    return f"{key[0]}:{key[1]}:{key[2]}>{key[3]}"


def _normalization_declared(path: Path) -> bool:
    with open_text(path) as handle:
        for line in handle:
            if not line.startswith("#"):
                return False
            if line.startswith("##bcftools_normCommand") or line.startswith("##candidate_qtl_normalized"):
                return True
    return False


def _alt_present(genotypes: list[str], alt_index: int) -> str:
    """Genotype state for one ALT allele across a record's samples."""
    if not genotypes:
        return "genotype_absent_in_vcf"
    seen_alt = False
    seen_call = False
    for genotype in genotypes:
        call = genotype.split(":", 1)[0].replace("|", "/")
        alleles = [value for value in call.split("/") if value]
        if not alleles or set(alleles) == {"."}:
            continue
        seen_call = True
        if str(alt_index) in alleles:
            seen_alt = True
    if seen_alt:
        return "alt_present"
    if seen_call:
        return "reference_or_other_allele"
    return "all_calls_missing"


def read_vcf_events(path: Path) -> dict[tuple[str, int, str, str], dict[str, object]]:
    """Read a VCF/gVCF into normalized biallelic events keyed by minimal representation."""
    events: dict[tuple[str, int, str, str], dict[str, object]] = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, _, ref, alt_field = fields[0], int(fields[1]), fields[2], fields[3], fields[4]
            genotypes = fields[9:] if len(fields) > 9 else []
            for index, alt in enumerate(alt_field.split(","), start=1):
                alt = alt.strip()
                if not alt or alt in {".", "<NON_REF>", "*"}:
                    continue
                key = normalized_event_key(chrom, pos, ref, alt)
                state = _alt_present(genotypes, index)
                previous = events.get(key)
                # A repeated key keeps the strongest observation; alt_present wins.
                if previous is None or (previous["genotype_state"] != "alt_present" and state == "alt_present"):
                    events[key] = {"chrom": key[0], "pos": key[1], "ref": key[2], "alt": key[3],
                                   "genotype_state": state, "raw_pos": pos, "raw_ref": ref, "raw_alt": alt}
    return events


# --------------------------------------------------------------------------- #
# callability
# --------------------------------------------------------------------------- #

class CallableRegions:
    """Per-chromosome sorted callable intervals from a BED file (0-based half-open)."""

    __slots__ = ("_by_chrom", "interval_count", "source")

    def __init__(self, path: Path | None) -> None:
        self._by_chrom: dict[str, tuple[list[int], list[int]]] = {}
        self.interval_count = 0
        self.source = str(path) if path else ""
        if path is None:
            return
        grouped: dict[str, list[tuple[int, int]]] = {}
        with open_text(path) as handle:
            for raw in handle:
                if not raw.strip() or raw.startswith(("#", "track", "browser")):
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) < 3:
                    continue
                grouped.setdefault(fields[0], []).append((int(fields[1]) + 1, int(fields[2])))
        for chrom, spans in grouped.items():
            spans.sort()
            merged: list[tuple[int, int]] = []
            for start, end in spans:
                if merged and start <= merged[-1][1] + 1:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end))
                else:
                    merged.append((start, end))
            self._by_chrom[chrom] = ([span[0] for span in merged], [span[1] for span in merged])
            self.interval_count += len(merged)

    @property
    def supplied(self) -> bool:
        return bool(self.source)

    def covers(self, chrom: str, start: int, end: int) -> bool:
        """True only when the whole 1-based inclusive span is callable."""
        bucket = self._by_chrom.get(str(chrom))
        if not bucket:
            return False
        starts, ends = bucket
        index = bisect.bisect_right(starts, start) - 1
        return index >= 0 and ends[index] >= end

    def callable_bp_span(self, chrom: str, start: int, end: int) -> int:
        """Callable base pairs inside a 1-based inclusive span."""
        bucket = self._by_chrom.get(str(chrom))
        if not bucket:
            return 0
        starts, ends = bucket
        total = 0
        index = bisect.bisect_left(starts, start)
        if index > 0 and ends[index - 1] >= start:
            index -= 1
        while index < len(starts) and starts[index] <= end:
            total += max(0, min(end, ends[index]) - max(start, starts[index]) + 1)
            index += 1
        return total


# --------------------------------------------------------------------------- #
# reference annotation
# --------------------------------------------------------------------------- #

def _read_reference_genes(path: Path) -> tuple[list[dict[str, object]], dict[str, list[tuple[str, int, int]]]]:
    """Return interval-agnostic gene records and per-gene CDS spans."""
    genes: list[dict[str, object]] = []
    cds: dict[str, list[tuple[str, int, int]]] = {}
    exons: dict[str, list[tuple[str, int, int]]] = {}
    transcript_gene: dict[str, str] = {}
    with open_text(path) as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                continue
            chrom, feature, start, end, strand = fields[0], fields[2], int(fields[3]), int(fields[4]), fields[6]
            attrs = attributes(fields[8])
            if feature == "gene":
                gene_id = attrs.get("ID") or attrs.get("Name")
                if gene_id:
                    genes.append({"gene_id": gene_id, "chrom": chrom, "start": start, "end": end, "strand": strand})
            elif feature in {"mRNA", "transcript"}:
                if attrs.get("ID") and attrs.get("Parent"):
                    transcript_gene[attrs["ID"]] = attrs["Parent"].split(",")[0]
            elif feature in {"CDS", "exon"}:
                target = cds if feature == "CDS" else exons
                for parent in attrs.get("Parent", "").split(","):
                    gene_id = transcript_gene.get(parent, parent)
                    if gene_id:
                        target.setdefault(gene_id, []).append((chrom, start, end))
    for gene_id, spans in exons.items():
        cds.setdefault(f"exon::{gene_id}", []).extend(spans)
    return genes, cds


def _merge_spans(spans: list[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    merged: list[tuple[str, int, int]] = []
    for chrom, start, end in sorted(spans):
        if merged and merged[-1][0] == chrom and start <= merged[-1][2] + 1:
            merged[-1] = (chrom, merged[-1][1], max(merged[-1][2], end))
        else:
            merged.append((chrom, start, end))
    return merged


class _GeneIndex:
    """Start-sorted per-chromosome gene index for event -> gene assignment."""

    def __init__(self, genes: list[dict[str, object]]) -> None:
        grouped: dict[str, list[tuple[int, int, dict[str, object]]]] = {}
        for gene in genes:
            grouped.setdefault(str(gene["chrom"]), []).append((int(gene["start"]), int(gene["end"]), gene))
        self._chroms: dict[str, tuple[list[int], list[tuple[int, int, dict[str, object]]], int]] = {}
        for chrom, rows in grouped.items():
            rows.sort(key=lambda row: (row[0], row[1], str(row[2]["gene_id"])))
            span = max((row[1] - row[0] for row in rows), default=0) + 1
            self._chroms[chrom] = ([row[0] for row in rows], rows, span)

    def overlapping(self, chrom: str, start: int, end: int) -> list[dict[str, object]]:
        bucket = self._chroms.get(str(chrom))
        if not bucket:
            return []
        starts, rows, span = bucket
        low = bisect.bisect_left(starts, start - span)
        high = bisect.bisect_right(starts, end)
        return [row[2] for row in rows[low:high] if row[1] >= start]


# --------------------------------------------------------------------------- #
# configuration
# --------------------------------------------------------------------------- #

def _load_model(case_arg: str | None, overrides: dict[str, object]) -> dict[str, object]:
    model: dict[str, object] = {}
    if case_arg:
        case = yaml.safe_load(runtime_path(case_arg).read_text(encoding="utf-8")) or {}
        if not isinstance(case, dict):
            raise ValueError("case manifest must be a mapping")
        model = dict(case.get("private_variant_model") or {})
        model.setdefault("case_id", case.get("case_id", ""))
        source = case.get("source") or {}
        intervals = source.get("target_intervals") or source.get("target_interval")
        if intervals and not model.get("loci"):
            values = intervals if isinstance(intervals, list) else [intervals]
            model["intervals"] = [str(value) for value in values]
    for key, value in overrides.items():
        if value:
            model[key] = value
    if not model.get("focal_accession"):
        raise ValueError("private_variant_model requires focal_accession")
    if not model.get("phenotype_controls"):
        raise ValueError("private_variant_model requires at least one phenotype control")
    return model


def _parse_named(values: list[str] | None, flag: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for item in values or []:
        if "=" not in item:
            raise ValueError(f"{flag} must be ACCESSION=PATH")
        name, path = item.split("=", 1)
        if not name or not path:
            raise ValueError(f"{flag} must be ACCESSION=PATH")
        parsed[name] = path
    return parsed


def _loci(model: dict[str, object]) -> list[dict[str, object]]:
    """Per-locus phenotype controls; a locus may override the case-level set."""
    default_controls = [str(name) for name in model.get("phenotype_controls") or []]
    declared = model.get("loci") or []
    loci: list[dict[str, object]] = []
    for entry in declared:
        if not isinstance(entry, dict) or not entry.get("interval"):
            raise ValueError("each locus requires locus_id and interval")
        chrom, start, end = parse_interval(str(entry["interval"]))
        loci.append({"locus_id": str(entry.get("locus_id") or entry["interval"]), "interval": str(entry["interval"]),
                     "chrom": chrom, "start": start, "end": end,
                     "phenotype_controls": [str(name) for name in (entry.get("phenotype_controls") or default_controls)]})
    if loci:
        return loci
    for interval in model.get("intervals") or []:
        chrom, start, end = parse_interval(str(interval))
        loci.append({"locus_id": str(interval), "interval": str(interval), "chrom": chrom, "start": start, "end": end,
                     "phenotype_controls": list(default_controls)})
    if not loci:
        loci.append({"locus_id": "genome_wide", "interval": "", "chrom": "", "start": 0, "end": 0,
                     "phenotype_controls": list(default_controls)})
    return loci


def _locus_for(loci: list[dict[str, object]], chrom: str, pos: int) -> dict[str, object]:
    for locus in loci:
        if not locus["interval"]:
            return locus
        if str(locus["chrom"]) == str(chrom) and int(locus["start"]) <= pos <= int(locus["end"]):
            return locus
    return {"locus_id": "outside_declared_loci", "interval": "", "chrom": "", "start": 0, "end": 0,
            "phenotype_controls": list(loci[0]["phenotype_controls"]) if loci else []}


# --------------------------------------------------------------------------- #
# supplied consequence frames
# --------------------------------------------------------------------------- #

def _load_consequences(path_arg: str | None, frame_label: str) -> dict[tuple[str, int, str, str], dict[str, set[str]]]:
    """Per-event, per-gene consequence terms from a supplied annotation table."""
    if not path_arg:
        return {}
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    required = {"chrom", "pos", "ref", "alt", "consequence"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{frame_label} consequence table is missing columns: {', '.join(missing)}")
    result: dict[tuple[str, int, str, str], dict[str, set[str]]] = {}
    for row in frame.to_dict("records"):
        key = normalized_event_key(str(row["chrom"]), int(row["pos"]), str(row["ref"]), str(row["alt"]))
        gene = str(row.get("gene_id", ""))
        terms = {term.strip() for term in str(row["consequence"]).replace(",", "&").split("&") if term.strip()}
        result.setdefault(key, {}).setdefault(gene, set()).update(terms)
    return result


_PROTEIN_ALTERING = {
    "missense_variant", "stop_gained", "stop_lost", "start_lost", "frameshift_variant",
    "inframe_insertion", "inframe_deletion", "splice_donor_variant", "splice_acceptor_variant",
    "protein_altering_variant", "coding_sequence_variant", "transcript_ablation",
}


def _is_protein_altering(terms: set[str]) -> bool:
    return bool({term.casefold() for term in terms} & _PROTEIN_ALTERING)


# --------------------------------------------------------------------------- #
# main entry point
# --------------------------------------------------------------------------- #

def multigenome_private(
    out_dir_arg: str,
    focal_vcf_arg: str,
    gff3_arg: str,
    case_arg: str | None = None,
    control_vcf_args: list[str] | None = None,
    callable_bed_args: list[str] | None = None,
    projection_ledger_arg: str | None = None,
    native_genes_arg: str | None = None,
    reference_consequence_arg: str | None = None,
    native_consequence_arg: str | None = None,
    focal_accession: str | None = None,
    phenotype_controls: list[str] | None = None,
    comparison_accessions: list[str] | None = None,
    intervals: list[str] | None = None,
    case_id: str = "",
    full_hash: bool = False,
) -> dict[str, object]:
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    model = _load_model(case_arg, {
        "focal_accession": focal_accession,
        "phenotype_controls": phenotype_controls,
        "comparison_accessions": comparison_accessions,
        # Resolved intervals (e.g. from a case-run marker liftover) win over the
        # ones written in the case manifest.
        "intervals": intervals,
        "case_id": case_id,
    })
    focal = str(model["focal_accession"])
    required_controls = [str(name) for name in model["phenotype_controls"]]
    comparison = [str(name) for name in (model.get("comparison_accessions") or []) if str(name) not in required_controls]
    require_callable = bool(model.get("require_all_controls_callable", True))
    loci = _loci(model)

    control_vcfs = _parse_named(control_vcf_args, "--control-vcf")
    declared_resources = {str(name): str(path) for name, path in (model.get("callability_resources") or {}).items()}
    callable_beds = {**declared_resources, **_parse_named(callable_bed_args, "--callable-bed")}

    # Per-locus overrides may name a control the case-level lists do not.
    locus_controls_declared = sorted({name for locus in loci for name in locus["phenotype_controls"]})
    all_declared = sorted({*required_controls, *comparison, *locus_controls_declared})
    required_declared = sorted({*required_controls, *locus_controls_declared})

    review_rows: list[dict[str, str]] = []
    for name in all_declared:
        if name not in control_vcfs:
            review_rows.append({"reason": "control_callset_missing", "target": name,
                                "required_next_action": f"supply --control-vcf {name}=<vcf>; without it every event stays not_evaluated for this accession",
                                "allowed_choices": "supply control callset | remove the accession from the manifest"})
    for name in required_declared:
        if name not in callable_beds:
            review_rows.append({"reason": "control_callability_missing", "target": name,
                                "required_next_action": f"supply --callable-bed {name}=<bed.gz>; absence of a VCF record is not evidence of a reference allele",
                                "allowed_choices": "supply callable regions | keep events as private_possible"})

    focal_path = runtime_path(focal_vcf_arg)
    focal_events = read_vcf_events(focal_path)
    normalization_declared = {focal: _normalization_declared(focal_path)}
    controls: dict[str, dict[tuple[str, int, str, str], dict[str, object]]] = {}
    for name, path_arg in sorted(control_vcfs.items()):
        path = runtime_path(path_arg)
        controls[name] = read_vcf_events(path)
        normalization_declared[name] = _normalization_declared(path)
    for name, declared in sorted(normalization_declared.items()):
        if not declared:
            review_rows.append({"reason": "normalization_not_declared", "target": name,
                                "required_next_action": "run `bcftools norm -f <reference> -m -any` and keep the header provenance; keys are trimmed to a minimal representation but not left-shifted",
                                "allowed_choices": "renormalize and rerun | accept minimal-representation matching only"})

    callability = {name: CallableRegions(runtime_path(callable_beds[name]) if name in callable_beds else None)
                   for name in all_declared}

    genes, cds_spans = _read_reference_genes(runtime_path(gff3_arg))
    gene_index = _GeneIndex(genes)
    reference_consequences = _load_consequences(reference_consequence_arg, "reference-frame")
    native_consequences = _load_consequences(native_consequence_arg, "native-frame")

    sharing_rows: list[dict[str, object]] = []
    matrix_rows: list[dict[str, object]] = []
    callability_stats = {name: {"callable": 0, "uncallable": 0, "missing": 0, "blocked": 0} for name in all_declared}

    for key in sorted(focal_events):
        event = focal_events[key]
        chrom, pos, ref, alt = key
        end = pos + max(len(ref), 1) - 1
        locus = _locus_for(loci, chrom, pos)
        locus_controls = [str(name) for name in locus["phenotype_controls"]]
        overlapping = gene_index.overlapping(chrom, pos, end)
        gene_ids = sorted({str(gene["gene_id"]) for gene in overlapping})
        region_class = _reference_region_class(gene_ids, cds_spans, chrom, pos, end)

        shared_with: list[str] = []
        absent_callable: list[str] = []
        absent_uncallable: list[str] = []
        not_evaluated_in: list[str] = []
        # Every declared accession is reconciled for every event. A case-level
        # control that this locus does not require is still recorded, as a
        # comparison accession rather than a phenotype control.
        for name in locus_controls + [item for item in all_declared if item not in locus_controls]:
            role = "phenotype_control" if name in locus_controls else "comparison"
            regions = callability.get(name) or CallableRegions(None)
            if name not in controls:
                state, resolution = "callset_absent", "not_evaluated"
                callable_state = "not_evaluated"
            else:
                control_event = controls[name].get(key)
                state = str(control_event["genotype_state"]) if control_event else "no_record"
                if control_event and state in {"alt_present", "genotype_absent_in_vcf"}:
                    resolution, callable_state = "shared", "not_required"
                else:
                    if not regions.supplied:
                        callable_state, resolution = "callability_not_supplied", "absent_uncallable"
                    elif regions.covers(chrom, pos, end):
                        callable_state, resolution = "callable", "absent_callable"
                    else:
                        callable_state, resolution = "uncallable", "absent_uncallable"
            if role == "phenotype_control":
                if resolution == "shared":
                    shared_with.append(name)
                elif resolution == "absent_callable":
                    absent_callable.append(name)
                elif resolution == "absent_uncallable":
                    absent_uncallable.append(name)
                else:
                    not_evaluated_in.append(name)
                stats = callability_stats.setdefault(name, {"callable": 0, "uncallable": 0, "missing": 0, "blocked": 0})
                if callable_state == "callable":
                    stats["callable"] += 1
                elif callable_state == "uncallable":
                    stats["uncallable"] += 1
                elif callable_state == "callability_not_supplied":
                    stats["missing"] += 1
            matrix_rows.append({
                "event_key": _event_id(key), "control_accession": name, "control_role": role,
                "genotype_state": state, "event_observed": resolution == "shared",
                "callability_state": callable_state, "callability_source": regions.source,
                "resolution": resolution,
                "notes": "absence of a VCF record is not a callable reference allele" if resolution == "absent_uncallable" else "",
            })

        if not locus_controls or not_evaluated_in:
            event_class = "not_evaluated"
        elif shared_with:
            event_class = "shared_confirmed"
        elif absent_uncallable:
            event_class = "private_possible"
        else:
            event_class = "private_confirmed"
        if event_class == "private_possible":
            for name in absent_uncallable:
                callability_stats.setdefault(name, {"callable": 0, "uncallable": 0, "missing": 0, "blocked": 0})["blocked"] += 1

        reference_terms = set().union(*reference_consequences.get(key, {}).values()) if key in reference_consequences else set()
        native_terms = set().union(*native_consequences.get(key, {}).values()) if key in native_consequences else set()
        if reference_terms and native_terms:
            agreement = "agree" if reference_terms == native_terms else ("partial" if reference_terms & native_terms else "disagree")
        elif reference_terms:
            agreement = "reference_frame_only"
        elif native_terms:
            agreement = "native_frame_only"
        else:
            agreement = "not_evaluated"
        protein_altering = _is_protein_altering(reference_terms) or _is_protein_altering(native_terms)
        review_status = "review_required" if agreement in {"partial", "disagree"} else "context"
        if agreement in {"partial", "disagree"}:
            review_rows.append({"reason": "dual_frame_consequence_disagreement", "target": _event_id(key),
                                "required_next_action": "reconcile reference-frame and native-frame CDS effect before using this event",
                                "allowed_choices": "recheck gene model | recheck projection | keep both frames as separate observations"})

        sharing_rows.append({
            "event_key": _event_id(key), "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
            "focal_accession": focal, "locus_id": locus["locus_id"], "phenotype_controls": ";".join(locus_controls),
            "event_class": event_class, "shared_with": ";".join(shared_with),
            "absent_callable_in": ";".join(absent_callable), "absent_uncallable_in": ";".join(absent_uncallable),
            "not_evaluated_in": ";".join(not_evaluated_in), "gene_ids": ";".join(gene_ids),
            "reference_region_class": region_class,
            "reference_consequence": "&".join(sorted(reference_terms)), "native_consequence": "&".join(sorted(native_terms)),
            "consequence_agreement": agreement, "protein_altering": protein_altering,
            "rank_contribution": "none", "review_status": review_status,
            "notes": "focal genotype: " + str(event["genotype_state"]),
        })

    sharing = _frame(sharing_rows, SHARING_COLUMNS)
    matrix = _frame(matrix_rows, CONTROL_MATRIX_COLUMNS)

    projection = _load_projection(projection_ledger_arg)
    gene_ledger, denominators, panel, extra_reviews = _gene_level_tables(
        sharing, genes, cds_spans, loci, projection, callability, required_declared,
        native_genes_arg, model)
    review_rows.extend(extra_reviews)

    summary_rows: list[dict[str, object]] = []
    for locus in loci + [{"locus_id": "outside_declared_loci"}]:
        locus_id = str(locus["locus_id"])
        scope = sharing[sharing["locus_id"] == locus_id] if not sharing.empty else sharing
        if scope.empty and locus_id == "outside_declared_loci":
            continue
        for event_class in EVENT_CLASSES:
            subset = scope[scope["event_class"] == event_class] if not scope.empty else scope
            genes_here = sorted({gene for value in subset["gene_ids"] for gene in str(value).split(";") if gene}) if not subset.empty else []
            summary_rows.append({
                "locus_id": locus_id, "event_class": event_class, "event_count": len(subset),
                "gene_count": len(genes_here),
                "protein_altering_events": int(subset["protein_altering"].astype(str).eq("True").sum()) if not subset.empty else 0,
                "notes": "counts are per callset-reconciled event; callsets are never summed",
            })
    summary = _frame(summary_rows, EVENT_SUMMARY_COLUMNS)

    reviews = _frame(review_rows, REVIEW_COLUMNS).drop_duplicates().sort_values(["reason", "target"], kind="stable")

    callability_rows = []
    for name in all_declared:
        regions = callability.get(name) or CallableRegions(None)
        stats = callability_stats.get(name, {"callable": 0, "uncallable": 0, "missing": 0, "blocked": 0})
        callable_bp = sum(regions.callable_bp_span(str(locus["chrom"]), int(locus["start"]), int(locus["end"]))
                          for locus in loci if locus["interval"])
        callability_rows.append({
            "control_accession": name, "control_role": "phenotype_control" if name in required_declared else "comparison",
            "callability_resource": regions.source, "callable_intervals": regions.interval_count,
            "callable_bp_in_intervals": callable_bp, "events_callable": stats["callable"],
            "events_uncallable": stats["uncallable"], "events_without_callability": stats["missing"],
            "private_confirmed_blocked_by_callability": stats["blocked"],
            "notes": "" if regions.supplied else "no callability resource: every absent event stays private_possible at best",
        })
    callability_frame = _frame(callability_rows, CALLABILITY_COLUMNS)

    paths = {
        "focal_event_sharing_ledger": out / "focal_event_sharing_ledger.tsv",
        "candidate_control_matrix": out / "candidate_control_matrix.tsv",
        "control_callability_sensitivity": out / "control_callability_sensitivity.tsv",
        "event_class_summary": out / "event_class_summary.tsv",
        "gene_private_variant_ledger": out / "gene_private_variant_ledger.tsv",
        "projection_and_callable_denominators": out / "projection_and_callable_denominators.tsv",
        "working_candidate_panel": out / "working_candidate_panel.tsv",
        "review_required": out / "review_required.tsv",
    }
    for frame, path in ((sharing, paths["focal_event_sharing_ledger"]), (matrix, paths["candidate_control_matrix"]),
                        (callability_frame, paths["control_callability_sensitivity"]), (summary, paths["event_class_summary"]),
                        (gene_ledger, paths["gene_private_variant_ledger"]), (denominators, paths["projection_and_callable_denominators"]),
                        (panel, paths["working_candidate_panel"]), (reviews, paths["review_required"])):
        frame.to_csv(path, sep="\t", index=False, lineterminator="\n")

    class_counts = {value: int((sharing["event_class"] == value).sum()) if not sharing.empty else 0 for value in EVENT_CLASSES}
    manifest = {
        "schema_version": "1.0",
        "case_id": str(model.get("case_id") or case_id),
        "focal_accession": focal,
        "phenotype_controls": required_controls,
        "comparison_accessions": comparison,
        "loci": [{"locus_id": locus["locus_id"], "interval": locus["interval"],
                  "phenotype_controls": locus["phenotype_controls"]} for locus in loci],
        "event_normalization": {
            "declared_rule": str(model.get("event_normalization", "left_normalized_ref_alt")),
            "applied_key": "minimal representation (shared suffix then prefix trimmed); multiallelic ALTs split per record",
            "left_shift_applied": False,
            "inputs_declaring_normalization": {name: bool(value) for name, value in sorted(normalization_declared.items())},
        },
        "require_all_controls_callable": require_callable,
        "callability_policy": (
            "a phenotype control that is uncallable, or has no callability resource, can never yield private_confirmed; "
            "a missing VCF record alone is never a callable reference allele"
        ),
        "callability_providers": {name: (callability[name].source or None) for name in sorted(callability)},
        "reference_assembly": str(model.get("reference_assembly", "")),
        "input_fingerprints": {
            "focal_vcf": file_fingerprint(focal_path, full_hash=full_hash),
            "gff3": file_fingerprint(runtime_path(gff3_arg), full_hash=full_hash),
            **{f"control_vcf:{name}": file_fingerprint(runtime_path(path), full_hash=full_hash) for name, path in sorted(control_vcfs.items())},
            **{f"callable_bed:{name}": file_fingerprint(runtime_path(path), full_hash=full_hash) for name, path in sorted(callable_beds.items())},
        },
        "focal_event_count": len(focal_events),
        "event_class_counts": class_counts,
        "gene_count": int(len(gene_ledger)),
        "review_required_count": int(len(reviews)),
        "outputs": {name: str(path) for name, path in sorted(paths.items())},
        "safety": [
            "a missing VCF record is never converted into a callable reference allele",
            "a failed or incomplete projection is never converted into gene absence",
            "private_possible never removes a candidate; it becomes an additional-callability action",
            "private status contributes no score (rank_contribution is always none) and promotes nothing automatically",
            "counts from different callsets are never summed; every event is reconciled per control",
        ],
    }
    manifest_path = out / "private_variant_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "manifest": str(manifest_path), "event_class_counts": class_counts,
        "focal_event_count": len(focal_events), "gene_count": int(len(gene_ledger)),
        "review_required_count": int(len(reviews)),
        **{name: str(path) for name, path in paths.items()},
    }


def _frame(rows: list[dict[str, object]], columns: list[str]) -> pd.DataFrame:
    if not rows:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(rows)[columns]


def _reference_region_class(gene_ids: list[str], cds_spans: dict[str, list[tuple[str, int, int]]],
                            chrom: str, start: int, end: int) -> str:
    if not gene_ids:
        return "intergenic"
    for gene_id in gene_ids:
        for span_chrom, span_start, span_end in cds_spans.get(gene_id, []):
            if span_chrom == chrom and span_end >= start and span_start <= end:
                return "cds"
    for gene_id in gene_ids:
        for span_chrom, span_start, span_end in cds_spans.get(f"exon::{gene_id}", []):
            if span_chrom == chrom and span_end >= start and span_start <= end:
                return "exon_noncoding"
    return "intragenic_noncoding"


def _load_projection(path_arg: str | None) -> dict[str, dict[str, str]]:
    """Per reference gene projection state from a projection ledger."""
    if not path_arg:
        return {}
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    if "gene_id" not in frame.columns:
        raise ValueError("projection ledger requires gene_id")
    if "projection_state" not in frame.columns:
        raise ValueError("projection ledger requires projection_state "
                         f"(one of: {', '.join(PROJECTION_STATES)})")
    result: dict[str, dict[str, str]] = {}
    for row in frame.to_dict("records"):
        state = str(row["projection_state"]).strip()
        complete = str(row.get("projection_run_complete", "true")).strip().casefold() not in {"false", "0", "no", ""}
        # Only a complete projection run may call a gene unmapped; an incomplete
        # run leaves it not_evaluated instead of implying gene absence.
        if state == "unmapped" and not complete:
            state = "not_evaluated"
        if state not in PROJECTION_STATES:
            state = "not_evaluated"
        result[str(row["gene_id"])] = {"projection_state": state, "projection_run_complete": str(complete)}
    return result


def _gene_level_tables(
    sharing: pd.DataFrame,
    genes: list[dict[str, object]],
    cds_spans: dict[str, list[tuple[str, int, int]]],
    loci: list[dict[str, object]],
    projection: dict[str, dict[str, str]],
    callability: dict[str, CallableRegions],
    required_controls: list[str],
    native_genes_arg: str | None,
    model: dict[str, object],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[dict[str, str]]]:
    reviews: list[dict[str, str]] = []
    gene_by_id = {str(gene["gene_id"]): gene for gene in genes}
    per_gene: dict[str, dict[str, object]] = {}
    if not sharing.empty:
        for row in sharing.to_dict("records"):
            for gene_id in str(row["gene_ids"]).split(";"):
                if not gene_id:
                    continue
                bucket = per_gene.setdefault(gene_id, {
                    "locus_id": str(row["locus_id"]), "private_confirmed": 0, "private_possible": 0,
                    "shared_confirmed": 0, "not_evaluated": 0, "protein_altering_private_confirmed": 0,
                    "regulatory_private_events": 0, "consequence_disagreements": 0, "protein_effect_evaluated": False,
                })
                bucket[str(row["event_class"])] = int(bucket[str(row["event_class"])]) + 1
                protein_altering = str(row["protein_altering"]) == "True"
                if str(row["consequence_agreement"]) != "not_evaluated":
                    bucket["protein_effect_evaluated"] = True
                if str(row["consequence_agreement"]) in {"partial", "disagree"}:
                    bucket["consequence_disagreements"] = int(bucket["consequence_disagreements"]) + 1
                if row["event_class"] == "private_confirmed" and protein_altering:
                    bucket["protein_altering_private_confirmed"] = int(bucket["protein_altering_private_confirmed"]) + 1
                if row["event_class"] in {"private_confirmed", "private_possible"} and str(row["reference_region_class"]) != "cds":
                    bucket["regulatory_private_events"] = int(bucket["regulatory_private_events"]) + 1

    gene_rows: list[dict[str, object]] = []
    panel_rows: list[dict[str, object]] = []
    for gene_id in sorted(per_gene):
        bucket = per_gene[gene_id]
        gene = gene_by_id.get(gene_id, {})
        state = projection.get(gene_id, {}).get("projection_state", "not_evaluated")
        gene_rows.append({
            "gene_id": gene_id, "locus_id": bucket["locus_id"], "chrom": str(gene.get("chrom", "")),
            "start": gene.get("start", ""), "end": gene.get("end", ""), "projection_state": state,
            "protein_effect_evaluated": bool(bucket["protein_effect_evaluated"]),
            "private_confirmed": bucket["private_confirmed"], "private_possible": bucket["private_possible"],
            "shared_confirmed": bucket["shared_confirmed"], "not_evaluated": bucket["not_evaluated"],
            "protein_altering_private_confirmed": bucket["protein_altering_private_confirmed"],
            "regulatory_private_events": bucket["regulatory_private_events"],
            "consequence_disagreements": bucket["consequence_disagreements"],
            "rank_contribution": "none",
            "review_status": "review_required" if bucket["private_possible"] or bucket["consequence_disagreements"] else "context",
            "notes": "private status is a triage state, not a score contribution",
        })
        group, gate, action = _candidate_group(bucket, state)
        panel_rows.append({
            "gene_id": gene_id, "locus_id": bucket["locus_id"], "candidate_group": group,
            "private_confirmed": bucket["private_confirmed"], "private_possible": bucket["private_possible"],
            "protein_altering_private_confirmed": bucket["protein_altering_private_confirmed"],
            "projection_state": state, "rank_contribution": "none", "gate": gate,
            "required_next_action": action,
            "notes": "multiple co_primary candidates are allowed; private_possible never removes a candidate",
        })

    native_gene_count = 0
    if native_genes_arg:
        native = pd.read_csv(runtime_path(native_genes_arg), sep="\t", dtype=str).fillna("")
        native_gene_count = int(native["gene_id"].nunique()) if "gene_id" in native.columns else len(native)

    denominator_rows: list[dict[str, object]] = []
    for locus in loci:
        if not locus["interval"]:
            continue
        chrom, start, end = str(locus["chrom"]), int(locus["start"]), int(locus["end"])
        in_interval = [gene for gene in genes if str(gene["chrom"]) == chrom and int(gene["end"]) >= start and int(gene["start"]) <= end]
        states = [projection.get(str(gene["gene_id"]), {}).get("projection_state", "not_evaluated") for gene in in_interval]
        intended_bp = 0
        callable_bp = 0
        for gene in in_interval:
            spans = _merge_spans([span for span in cds_spans.get(str(gene["gene_id"]), []) if span[0] == chrom])
            for span_chrom, span_start, span_end in spans:
                span_start, span_end = max(span_start, start), min(span_end, end)
                if span_end < span_start:
                    continue
                intended_bp += span_end - span_start + 1
                covered = [callability[name].callable_bp_span(span_chrom, span_start, span_end)
                           for name in required_controls if name in callability and callability[name].supplied]
                callable_bp += min(covered) if covered else 0
        high_confidence = states.count("high_confidence_one_to_one")
        any_hit = high_confidence + states.count("any_hit")
        evaluated = sum(1 for gene in in_interval if per_gene.get(str(gene["gene_id"]), {}).get("protein_effect_evaluated"))
        fraction = round(callable_bp / intended_bp, 4) if intended_bp else 0.0
        projection_fraction = high_confidence / len(in_interval) if in_interval else 0.0
        status = "sufficient"
        action = ""
        if in_interval and projection_fraction < _MIN_PROJECTION_FRACTION:
            status = "projection_incomplete"
            action = "run annotation-rescue or the SV/PAV lane before reading a low private-variant count as biological evidence"
        if intended_bp and fraction < _MIN_CALLABLE_CDS_FRACTION:
            status = "callability_incomplete" if status == "sufficient" else "projection_and_callability_incomplete"
            action = "extend control callability (coverage or assembly alignment) before treating absent events as reference alleles"
        if action:
            reviews.append({"reason": f"annotation_completeness_{status}", "target": str(locus["locus_id"]),
                            "required_next_action": action,
                            "allowed_choices": "run annotation-rescue | run SV/PAV lane | extend callability | keep denominators explicit"})
        denominator_rows.append({
            "locus_id": str(locus["locus_id"]), "interval": str(locus["interval"]),
            "reference_genes_in_interval": len(in_interval), "any_hit_projection": any_hit,
            "high_confidence_one_to_one": high_confidence, "partial_projection": states.count("partial"),
            "unmapped": states.count("unmapped"), "projection_not_evaluated": states.count("not_evaluated"),
            "native_annotation_genes": native_gene_count, "intended_cds_bp": intended_bp,
            "callable_cds_bp": callable_bp, "callable_cds_fraction": fraction,
            "protein_effect_evaluated_genes": evaluated,
            "denominator_status": status,
            "required_next_action": action,
        })

    if not projection and any(locus["interval"] for locus in loci):
        reviews.append({"reason": "projection_ledger_missing", "target": str(model.get("case_id") or "case"),
                        "required_next_action": "supply --projection-ledger; without it every reference gene stays projection not_evaluated",
                        "allowed_choices": "supply projection ledger | keep denominators as not_evaluated"})

    return (_frame(gene_rows, GENE_LEDGER_COLUMNS), _frame(denominator_rows, DENOMINATOR_COLUMNS),
            _frame(panel_rows, PANEL_COLUMNS), reviews)


def _candidate_group(bucket: dict[str, object], projection_state: str) -> tuple[str, str, str]:
    """Assign a working triage group. No group implies a score or a promotion."""
    if int(bucket["protein_altering_private_confirmed"]) > 0:
        return ("co_primary",
                "protein-altering private_confirmed event in every required phenotype control",
                "confirm the variant and its segregation; several co_primary candidates may coexist")
    if int(bucket["private_confirmed"]) > 0 and int(bucket["regulatory_private_events"]) > 0:
        return ("regulatory_sv",
                "private_confirmed event outside the reference CDS",
                "route through the promoter/UTR/SV review lane; this never adds score on its own")
    if int(bucket["private_confirmed"]) > 0:
        return ("reserve", "private_confirmed event without an evaluated protein effect",
                "supply reference-frame and native-frame consequences for this event")
    if int(bucket["private_possible"]) > 0:
        return ("reserve", "private status limited by control callability",
                "extend callability for the uncallable control, then re-evaluate; the candidate is not removed")
    if projection_state in {"not_evaluated", "partial", "unmapped"}:
        return ("annotation_rescue", f"projection state is {projection_state}",
                "run annotation rescue before reading the absence of private events as evidence")
    return ("evidence_insufficient", "no private event and no annotation gap identified",
            "no action; keep the gene visible with its denominators")
