"""Generic same-species cross-assembly evidence bridge (Issue #527 / #524).

This module is deliberately a small adapter between native LIS-like resource
artifacts and the existing Route B contracts. It never joins on a symbol or a
classical locus name. Every identity join is an exact full-ID join, and the
source functional card remains a source-native record until explicit mapping
evidence supplies a target ID.

LIS pangene membership may be supplied as a candidate-correspondence layer. It
only narrows which full-ID pairs need verification: pangene membership alone
never verifies one-copy identity.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Iterable

import pandas as pd

from .core import file_fingerprint, runtime_path

SCHEMA_VERSION = "1.2"
MAPPING_COLUMNS = [
    "source_full_id", "target_full_id", "source_namespace", "target_namespace",
    "mapping_scope", "mapping_basis", "relationship", "source_family_id",
    "target_family_id", "pangene_dataset", "pangene_set_id",
    "pangene_relation_present", "source_member_count_in_pangene",
    "target_member_count_in_pangene", "pangene_relationship",
    "reciprocal_sequence_state", "sequence_identity", "query_coverage",
    "subject_coverage", "synteny_state", "domain_concordance",
    "mapping_state", "mapping_reason",
]
ORTHOLOGY_COLUMNS = [
    "gene_id", "reference_gene_id", "orthogroup", "relationship", "provider",
    "percent_identity", "query_coverage", "evalue", "bitscore",
]
REFERENCE_FUNCTION_COLUMNS = [
    "reference_gene_id", "gene_id", "function_description", "pathway",
    "pathway_id", "mutant_phenotype", "source_id", "source_url",
    "reference_species", "gene_symbol", "evidence_type", "provenance_detail",
    "mapping_state", "mapping_relationship", "function_transfer_scope",
]
AUDIT_COLUMNS = [
    "source_full_id", "target_full_id", "audit_state", "relationship",
    "mapping_state", "reason",
]

_SOURCE_ID_COLUMNS = ("source_full_id", "source_gene_id", "source_id", "query_gene_id", "query_id", "source")
_TARGET_ID_COLUMNS = ("target_full_id", "target_gene_id", "target_id", "subject_gene_id", "subject_id", "target")
_GENE_COLUMNS = ("native_gene_id", "full_id", "gene_id", "gene", "source_full_id", "target_full_id")
_FAMILY_COLUMNS = ("family_id", "fam3_id", "fam3", "family", "family_stem")
_PANGENE_GENE_COLUMNS = ("gene", "gene_id", "native_gene_id", "member_gene_id")
_PANGENE_SET_COLUMNS = ("representative", "pangene_set_id", "representative_gene_id", "cluster_id")


def _column(frame: pd.DataFrame, candidates: Iterable[str], label: str) -> str:
    for name in candidates:
        if name in frame.columns:
            return name
    raise ValueError(f"{label} requires one of: {', '.join(candidates)}")


def _read(path_arg: str | Path, label: str) -> pd.DataFrame:
    path = runtime_path(path_arg)
    if not path.exists():
        raise ValueError(f"{label} does not exist: {path}")
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def _normalise_assignments(path_arg: str | Path, side: str) -> pd.DataFrame:
    frame = _read(path_arg, f"{side} fam3 assignment")
    gene_column = _column(frame, _GENE_COLUMNS, f"{side} fam3 assignment")
    family_column = _column(frame, _FAMILY_COLUMNS, f"{side} fam3 assignment")
    namespace_column = next((name for name in ("native_gene_id_namespace", "namespace") if name in frame.columns), "")
    rows: list[dict[str, str]] = []
    for row in frame.to_dict("records"):
        gene = str(row.get(gene_column, "")).strip()
        family = str(row.get(family_column, "")).strip()
        if not gene:
            continue
        rows.append({
            "gene": gene,
            "family": family,
            "namespace": str(row.get(namespace_column, "")).strip() if namespace_column else _native_namespace(gene, ""),
        })
    if not rows:
        raise ValueError(f"{side} fam3 assignment has no usable gene rows")
    return pd.DataFrame(rows).drop_duplicates().sort_values(["gene", "family"], kind="stable")


def _normalise_sequence(path_arg: str | Path) -> pd.DataFrame:
    frame = _read(path_arg, "pairwise sequence evidence")
    source_column = _column(frame, _SOURCE_ID_COLUMNS, "pairwise sequence evidence")
    target_column = _column(frame, _TARGET_ID_COLUMNS, "pairwise sequence evidence")
    state_column = next((name for name in ("reciprocal_sequence_state", "sequence_state", "reciprocal_state", "relationship") if name in frame.columns), "")
    identity_column = next((name for name in ("sequence_identity", "percent_identity", "pident", "identity") if name in frame.columns), "")
    query_column = next((name for name in ("query_coverage", "query_cov", "qcov") if name in frame.columns), "")
    subject_column = next((name for name in ("subject_coverage", "subject_cov", "scov") if name in frame.columns), "")
    if not state_column and not identity_column:
        raise ValueError("pairwise sequence evidence requires reciprocal_sequence_state or sequence_identity")
    synteny_column = next((name for name in ("synteny_state", "synteny_status") if name in frame.columns), "")
    domain_column = next((name for name in ("domain_concordance", "domain_state", "domain_status") if name in frame.columns), "")
    rows: list[dict[str, str]] = []
    for row in frame.to_dict("records"):
        source = str(row.get(source_column, "")).strip()
        target = str(row.get(target_column, "")).strip()
        if not source or not target:
            continue
        rows.append({
            "source": source, "target": target,
            "state": str(row.get(state_column, "")).strip() if state_column else "",
            "identity": str(row.get(identity_column, "")).strip() if identity_column else "",
            "query_coverage": str(row.get(query_column, "")).strip() if query_column else "",
            "subject_coverage": str(row.get(subject_column, "")).strip() if subject_column else "",
            "synteny": str(row.get(synteny_column, "")).strip() if synteny_column else "",
            "domain": str(row.get(domain_column, "")).strip() if domain_column else "",
        })
    return pd.DataFrame(rows, columns=["source", "target", "state", "identity", "query_coverage", "subject_coverage", "synteny", "domain"]).drop_duplicates()


def _normalise_pangene_hash(path_arg: str | Path) -> pd.DataFrame:
    """Load the LIS pangene hash table into gene -> set-ID pairs.

    LIS datastore-specifications documents the canonical hsh table as two
    tab-delimited columns with the set ID in the first column and a gene in the
    second. Some frozen adapters/fixtures expose a named header (``gene`` and
    ``representative``), so both the native headerless order and named variants
    are accepted. All joins remain exact full-ID joins. Pandas compression
    inference keeps the native ``*.hsh.tsv.gz`` artifact directly usable.
    """
    path = runtime_path(path_arg)
    if not path.exists():
        raise ValueError(f"pangene hash does not exist: {path}")
    # Native LIS *.hsh.tsv[.gz] files are headerless and ordered as
    # pangene_set_id<TAB>gene. Read the first two fields without interpreting
    # the first biological row as a header, then detect a named fixture header
    # when present.
    frame = pd.read_csv(path, sep="\t", dtype=str, compression="infer", header=None,
                        usecols=[0, 1]).fillna("")
    if frame.empty:
        return pd.DataFrame(columns=["gene", "pangene_set_id"])
    first = [str(value).strip().casefold() for value in frame.iloc[0].tolist()]
    gene_aliases = {value.casefold() for value in _PANGENE_GENE_COLUMNS}
    set_aliases = {value.casefold() for value in _PANGENE_SET_COLUMNS}
    if first[0] in gene_aliases and first[1] in set_aliases:
        gene_index, set_index = 0, 1
        frame = frame.iloc[1:]
    elif first[0] in set_aliases and first[1] in gene_aliases:
        gene_index, set_index = 1, 0
        frame = frame.iloc[1:]
    else:
        gene_index, set_index = 1, 0
    rows = []
    for row in frame.to_dict("records"):
        gene = str(row.get(gene_index, "")).strip()
        representative = str(row.get(set_index, "")).strip()
        if gene and representative:
            rows.append({"gene": gene, "pangene_set_id": representative})
    return pd.DataFrame(rows, columns=["gene", "pangene_set_id"]).drop_duplicates()


def _pangene_pairs(pangene: pd.DataFrame, source_genes: set[str], target_genes: set[str]) -> tuple[dict[tuple[str, str], dict[str, str]], dict[str, object]]:
    """Return exact source/target candidate pairs and namespace coverage.

    A representative can contain duplicated/split members from either side. The
    full Cartesian product is retained; multiplicity is evidence, not a reason
    to pick one copy.
    """
    pairs: dict[tuple[str, str], dict[str, str]] = {}
    source_seen: set[str] = set()
    target_seen: set[str] = set()
    for pangene_set_id, group in pangene.groupby("pangene_set_id", sort=True):
        members = set(str(value) for value in group["gene"])
        sources = sorted(members & source_genes)
        targets = sorted(members & target_genes)
        source_seen.update(sources)
        target_seen.update(targets)
        if not sources or not targets:
            continue
        source_count, target_count = len(sources), len(targets)
        relation = ("one_to_one_candidate" if source_count == 1 and target_count == 1
                    else "one_to_many_candidate" if source_count == 1 and target_count > 1
                    else "many_to_one_candidate" if source_count > 1 and target_count == 1
                    else "many_to_many_candidate")
        for source in sources:
            for target in targets:
                pairs[(source, target)] = {
                    "pangene_set_id": str(pangene_set_id),
                    "source_member_count_in_pangene": str(source_count),
                    "target_member_count_in_pangene": str(target_count),
                    "pangene_relationship": relation,
                }
    coverage = {
        "source_gene_count": len(source_genes),
        "target_gene_count": len(target_genes),
        "source_gene_in_pangene_count": len(source_seen),
        "target_gene_in_pangene_count": len(target_seen),
        "source_namespace_represented": bool(source_seen),
        "target_namespace_represented": bool(target_seen),
    }
    return pairs, coverage


def _native_namespace(value: str, supplied: str) -> str:
    if supplied:
        return supplied
    parts = value.split(".", 4)
    return ".".join(parts[:4]) if len(parts) == 5 else value.rsplit(".", 1)[0]


def _first_value(row: dict[str, str], *keys: str) -> str:
    for key in keys:
        value = str(row.get(key, "")).strip()
        if value:
            return value
    return ""


def _bridge_provenance(function_row: dict[str, str], map_row: dict[str, str],
                       source_id: str, target_id: str) -> str:
    """Preserve root evidence identity and an append-only, validated bridge path.

    Rich native metadata is preserved only when the source record supplied it;
    a parseable full gene ID is never decomposed to fabricate accession,
    assembly, annotation, or namespace metadata. Once bridge-owned provenance
    exists, it is immutable and must be internally coherent before another hop
    may be appended.
    """
    prior_text = str(function_row.get("provenance_detail", "")).strip()
    existing_bridge: dict[str, object] | None = None
    upstream_detail: object = ""
    if prior_text:
        try:
            parsed = json.loads(prior_text)
        except json.JSONDecodeError:
            upstream_detail = prior_text
        else:
            if isinstance(parsed, dict) and "assembly_bridge" in parsed:
                bridge = parsed.get("assembly_bridge")
                if not isinstance(bridge, dict):
                    raise ValueError("assembly_bridge provenance_detail must contain an object")
                existing_bridge = bridge
                upstream_detail = parsed.get("upstream_provenance_detail", "")
            else:
                upstream_detail = prior_text

    if existing_bridge is None:
        # First bridge: root identity comes from an explicit native gene when
        # present, otherwise from the exact source join ID. Do not mine the ID
        # string for unsupplied accession/assembly/annotation metadata.
        root_gene = _first_value(function_row, "native_gene_id") or source_id
        if not root_gene:
            raise ValueError("assembly_bridge root evidence requires a non-empty gene_id")
        evidence_native: dict[str, object] = {
            "gene_id": root_gene,
            "species": _first_value(function_row, "native_species", "reference_species"),
            "species_code": _first_value(function_row, "native_species_code"),
            "accession": _first_value(function_row, "native_accession"),
            "assembly": _first_value(function_row, "native_assembly"),
            "annotation": _first_value(function_row, "native_annotation"),
            "gene_id_namespace": _first_value(function_row, "native_gene_id_namespace"),
            "source_record_id": _first_value(function_row, "source_record_id"),
            "ontology_ids": _first_value(function_row, "ontology_ids"),
            "source_confidence": _first_value(function_row, "source_confidence"),
            "gene_symbol": _first_value(function_row, "gene_symbol"),
            "mutant_background": _first_value(function_row, "mutant_or_perturbation_background"),
        }
        mapping_path = [source_id]
        mapping_hops: list[dict[str, object]] = []
    else:
        existing_native = existing_bridge.get("evidence_native")
        if not isinstance(existing_native, dict):
            raise ValueError("assembly_bridge evidence_native must contain an object")
        root_gene = str(existing_native.get("gene_id", "")).strip()
        if not root_gene:
            raise ValueError("assembly_bridge evidence_native requires a non-empty gene_id")
        evidence_native = dict(existing_native)

        existing_path = existing_bridge.get("mapping_path")
        if not isinstance(existing_path, list):
            raise ValueError("assembly_bridge mapping_path must be a list")
        if not existing_path or any(not isinstance(value, str) or not value.strip() for value in existing_path):
            raise ValueError("assembly_bridge mapping_path requires non-empty string gene IDs")
        mapping_path = [value.strip() for value in existing_path]
        if mapping_path[0] != root_gene:
            raise ValueError("assembly_bridge mapping_path must start at the root evidence gene")
        if mapping_path[-1] != source_id:
            raise ValueError("assembly_bridge mapping_path does not terminate at the current source gene")

        existing_hops = existing_bridge.get("mapping_hops")
        if not isinstance(existing_hops, list):
            raise ValueError("assembly_bridge mapping_hops must be a list")
        mapping_hops = []
        for index, hop in enumerate(existing_hops):
            if not isinstance(hop, dict):
                raise ValueError("assembly_bridge mapping_hops entries must be objects")
            hop_source = str(hop.get("source_full_id", "")).strip()
            hop_target = str(hop.get("target_full_id", "")).strip()
            if not hop_source or not hop_target:
                raise ValueError("assembly_bridge mapping_hops entries require source_full_id and target_full_id")
            mapping_hops.append(dict(hop))
            if index + 1 >= len(mapping_path):
                raise ValueError("assembly_bridge mapping_hops is longer than mapping_path")
            if hop_source != mapping_path[index] or hop_target != mapping_path[index + 1]:
                raise ValueError("assembly_bridge mapping_hops does not match mapping_path")
        if len(mapping_hops) != len(mapping_path) - 1:
            raise ValueError("assembly_bridge mapping_hops must cover every mapping_path edge")

    if mapping_path[-1] != source_id:
        raise ValueError("assembly_bridge current source does not match provenance path")
    if target_id == source_id:
        raise ValueError("assembly_bridge target must differ from source")

    mapping_hops.append({
        "source_full_id": source_id,
        "target_full_id": target_id,
        "mapping_scope": str(map_row.get("mapping_scope", "")),
        "mapping_basis": str(map_row.get("mapping_basis", "")),
        "relationship": str(map_row.get("relationship", "")),
        "mapping_state": str(map_row.get("mapping_state", "")),
        "mapping_reason": str(map_row.get("mapping_reason", "")),
        "pangene_dataset": str(map_row.get("pangene_dataset", "")),
    })
    mapping_path.append(target_id)

    payload: dict[str, object] = {}
    if upstream_detail not in ("", None):
        payload["upstream_provenance_detail"] = upstream_detail
    payload["assembly_bridge"] = {
        "evidence_native": evidence_native,
        "mapping_path": mapping_path,
        "mapping_hops": mapping_hops,
    }
    return json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _state_is_reciprocal(value: str) -> bool:
    normalized = re.sub(r"[^a-z0-9]+", "_", value.casefold()).strip("_")
    return normalized in {
        "reciprocal", "unique_reciprocal", "reciprocal_best_hit", "rbh",
        "verified", "verified_reciprocal", "reciprocal_sequence_verified",
        "one_to_one", "true", "yes", "pass",
    }


def _state_is_conflict(value: str) -> bool:
    normalized = re.sub(r"[^a-z0-9]+", "_", value.casefold())
    return any(token in normalized for token in ("conflict", "incompatible", "contradict"))


def _mapping_reason(*, family_match: bool, relation: str, sequence: str, synteny: str, domain: str, has_sequence: bool, pangene_present: bool) -> str:
    if not family_match:
        if pangene_present and not has_sequence:
            return "pangene_candidate_without_fam3_concordance; identity_not_verified"
        return "family_mismatch; pairwise_sequence_is_not_silent_identity" if has_sequence else "source_and_target_fam3_mismatch"
    if _state_is_conflict(synteny):
        return "synteny_conflict_prevents_forced_assignment"
    if relation != "one_to_one":
        return f"family_concordance_but_{relation}"
    if _state_is_reciprocal(sequence):
        return "unique_reciprocal_sequence_and_family_concordance"
    if has_sequence:
        return "pairwise_sequence_present_but_not_explicitly_reciprocal"
    if pangene_present:
        return "pangene_candidate_and_fam3_concordance_without_reciprocal_sequence"
    return "same_fam3_without_reciprocal_sequence_evidence"


def _adapter_eligible(row: dict[str, str]) -> tuple[bool, str]:
    """Return whether a mapping may enter the existing Route B adapters.

    The full mapping table remains lossless. Adapters are intentionally
    narrower: a verified 1:1 is safe, and an ambiguous 1:many is safe only as
    a family-level set with concordant fam3 evidence. Provisional, family
    mismatch, synteny-conflict, and not-evaluated rows stay audit-only.
    """
    state = str(row.get("mapping_state", ""))
    relationship = str(row.get("relationship", ""))
    if state == "verified" and relationship == "one_to_one":
        return True, "verified_one_to_one"
    if (state == "ambiguous" and relationship == "one_to_many"
            and str(row.get("source_family_id", "")).strip()
            and str(row.get("target_family_id", "")).strip()
            and not _state_is_conflict(str(row.get("synteny_state", "")))
            and "fam3_mismatch" not in str(row.get("mapping_basis", ""))):
        return True, "ambiguous_one_to_many_family_set"
    return False, "mapping_state_not_safe_for_route_b_adapter"


def build_assembly_bridge(
    functional_evidence_arg: str | Path,
    source_fam3_arg: str | Path,
    target_fam3_arg: str | Path,
    sequence_evidence_arg: str | Path,
    out_dir_arg: str | Path,
    *,
    source_namespace: str = "",
    target_namespace: str = "",
    mapping_scope: str = "same_species_cross_assembly",
    case_id: str = "",
    pangene_hash_arg: str | Path | None = None,
    pangene_dataset: str = "",
) -> dict[str, object]:
    """Build mapping and Route B adapters from explicit native-ID evidence."""
    if not mapping_scope.strip():
        raise ValueError("mapping_scope must not be empty")
    if pangene_hash_arg and not pangene_dataset.strip():
        raise ValueError("pangene_hash requires an explicit pangene_dataset generation")
    functions = _read(functional_evidence_arg, "normalized functional evidence")
    function_gene_column = _column(functions, ("native_gene_id", "source_full_id", "gene_id"), "normalized functional evidence")
    source = _normalise_assignments(source_fam3_arg, "source")
    target = _normalise_assignments(target_fam3_arg, "target")
    sequence = _normalise_sequence(sequence_evidence_arg)
    source_gene_set = set(str(value) for value in source["gene"])
    target_gene_set = set(str(value) for value in target["gene"])
    pangene = _normalise_pangene_hash(pangene_hash_arg) if pangene_hash_arg else pd.DataFrame(columns=["gene", "pangene_set_id"])
    pangene_pairs, pangene_coverage = _pangene_pairs(pangene, source_gene_set, target_gene_set) if pangene_hash_arg else ({}, {
        "source_gene_count": len(source_gene_set), "target_gene_count": len(target_gene_set),
        "source_gene_in_pangene_count": 0, "target_gene_in_pangene_count": 0,
        "source_namespace_represented": False, "target_namespace_represented": False,
    })
    pangene_targets_by_source: dict[str, set[str]] = {}
    for source_id, target_id in pangene_pairs:
        pangene_targets_by_source.setdefault(source_id, set()).add(target_id)

    source_by_gene: dict[str, list[dict[str, str]]] = {}
    target_by_family: dict[str, list[dict[str, str]]] = {}
    target_by_gene: dict[str, list[dict[str, str]]] = {}
    for row in source.to_dict("records"):
        source_by_gene.setdefault(row["gene"], []).append(row)
    for row in target.to_dict("records"):
        target_by_gene.setdefault(row["gene"], []).append(row)
        if row["family"]:
            target_by_family.setdefault(row["family"], []).append(row)

    function_rows = functions.to_dict("records")
    function_ids = {str(row.get(function_gene_column, "")).strip() for row in function_rows if str(row.get(function_gene_column, "")).strip()}
    sequence_sources = set(sequence["source"]) if not sequence.empty else set()
    source_ids = sorted(function_ids | sequence_sources | set(source_by_gene))
    mapping_rows: list[dict[str, str]] = []
    audit_rows: list[dict[str, str]] = []

    # Multiplicity is calculated globally from all candidates that can actually
    # be considered for each source. If a cross-side pangene relation exists,
    # it is the candidate-generation layer and fam3 becomes a concordance check;
    # sources without a pangene candidate retain the old fam3 fallback. Explicit
    # sequence pairs are always retained so independent verification evidence is
    # never hidden by pangene coverage.
    source_to_targets: dict[str, set[str]] = {}
    target_to_sources: dict[str, set[str]] = {}
    for row in source.to_dict("records"):
        if pangene_hash_arg and row["gene"] in pangene_targets_by_source:
            continue
        for target_row in target_by_family.get(row["family"], []):
            source_to_targets.setdefault(row["gene"], set()).add(target_row["gene"])
            target_to_sources.setdefault(target_row["gene"], set()).add(row["gene"])
    for source_id, target_id in pangene_pairs:
        source_to_targets.setdefault(source_id, set()).add(target_id)
        target_to_sources.setdefault(target_id, set()).add(source_id)
    for row in sequence[["source", "target"]].drop_duplicates().to_dict("records") if not sequence.empty else []:
        source_to_targets.setdefault(str(row["source"]), set()).add(str(row["target"]))
        target_to_sources.setdefault(str(row["target"]), set()).add(str(row["source"]))

    for source_id in source_ids:
        source_assignments = source_by_gene.get(source_id, [])
        source_families = {row["family"] for row in source_assignments if row["family"]}
        family_targets = [row for family in sorted(source_families) for row in target_by_family.get(family, [])]
        explicit_pairs = sequence[sequence["source"].eq(source_id)] if not sequence.empty else pd.DataFrame()
        candidate_targets = set(str(value) for value in explicit_pairs["target"]) if not explicit_pairs.empty else set()
        pangene_targets = set(pangene_targets_by_source.get(source_id, set()))
        candidate_targets.update(pangene_targets)
        if not pangene_targets:
            candidate_targets.update(row["gene"] for row in family_targets)
        if not candidate_targets:
            source_namespace_value = source_namespace or (source_assignments[0]["namespace"] if source_assignments else _native_namespace(source_id, ""))
            mapping_rows.append({
                "source_full_id": source_id, "target_full_id": "", "source_namespace": source_namespace_value,
                "target_namespace": target_namespace, "mapping_scope": mapping_scope,
                "mapping_basis": "source_functional_evidence+fam3_assignments", "relationship": "unresolved",
                "source_family_id": ";".join(sorted(source_families)), "target_family_id": "",
                "pangene_dataset": pangene_dataset if pangene_hash_arg else "", "pangene_set_id": "",
                "pangene_relation_present": "false", "source_member_count_in_pangene": "0",
                "target_member_count_in_pangene": "0", "pangene_relationship": "",
                "reciprocal_sequence_state": "", "sequence_identity": "", "query_coverage": "",
                "subject_coverage": "", "synteny_state": "", "domain_concordance": "",
                "mapping_state": "not_evaluated",
                "mapping_reason": "source_or_target_fam3_assignment_missing; absence_is_not_negative_evidence",
            })
            audit_rows.append({"source_full_id": source_id, "target_full_id": "", "audit_state": "unresolved",
                               "relationship": "unresolved", "mapping_state": "not_evaluated",
                               "reason": "source_or_target_fam3_assignment_missing"})
            continue

        for target_id in sorted(candidate_targets):
            source_family_ids = sorted(source_families)
            target_assignments = target_by_gene.get(target_id, [])
            target_families = {row["family"] for row in target_assignments if row["family"]}
            family_match = bool(source_families & target_families)
            sequence_rows = explicit_pairs[explicit_pairs["target"].eq(target_id)] if not explicit_pairs.empty else pd.DataFrame()
            seq = sequence_rows.iloc[0].to_dict() if not sequence_rows.empty else {}
            relation = ("one_to_one" if len(source_to_targets.get(source_id, set())) == 1 and len(target_to_sources.get(target_id, set())) == 1
                        else "one_to_many" if len(source_to_targets.get(source_id, set())) > 1 and len(target_to_sources.get(target_id, set())) == 1
                        else "many_to_one" if len(source_to_targets.get(source_id, set())) == 1 and len(target_to_sources.get(target_id, set())) > 1
                        else "many_to_many")
            has_sequence = bool(seq)
            synteny = str(seq.get("synteny", ""))
            domain = str(seq.get("domain", ""))
            pangene_info = pangene_pairs.get((source_id, target_id), {})
            pangene_present = bool(pangene_info)
            state = ("verified" if family_match and relation == "one_to_one" and _state_is_reciprocal(str(seq.get("state", ""))) and not _state_is_conflict(synteny)
                     else "ambiguous" if relation != "one_to_one" or _state_is_conflict(synteny)
                     else "provisional")
            mapping_rows.append({
                "source_full_id": source_id,
                "target_full_id": target_id,
                "source_namespace": source_namespace or (source_assignments[0]["namespace"] if source_assignments else _native_namespace(source_id, "")),
                "target_namespace": target_namespace or (target_assignments[0]["namespace"] if target_assignments else _native_namespace(target_id, "")),
                "mapping_scope": mapping_scope,
                "mapping_basis": "+".join(value for value in (
                    "pangene_correspondence" if pangene_present else "",
                    "fam3_concordance" if family_match else "fam3_mismatch",
                    "pairwise_sequence" if has_sequence else "",
                    "synteny" if synteny else "",
                    "domain" if domain else "",
                ) if value),
                "relationship": relation,
                "source_family_id": ";".join(source_family_ids),
                "target_family_id": ";".join(sorted(target_families)),
                "pangene_dataset": pangene_dataset if pangene_present else "",
                "pangene_set_id": str(pangene_info.get("pangene_set_id", "")),
                "pangene_relation_present": "true" if pangene_present else "false",
                "source_member_count_in_pangene": str(pangene_info.get("source_member_count_in_pangene", "0")),
                "target_member_count_in_pangene": str(pangene_info.get("target_member_count_in_pangene", "0")),
                "pangene_relationship": str(pangene_info.get("pangene_relationship", "")),
                "reciprocal_sequence_state": str(seq.get("state", "")),
                "sequence_identity": str(seq.get("identity", "")),
                "query_coverage": str(seq.get("query_coverage", "")),
                "subject_coverage": str(seq.get("subject_coverage", "")),
                "synteny_state": synteny,
                "domain_concordance": domain,
                "mapping_state": state,
                "mapping_reason": _mapping_reason(family_match=family_match, relation=relation,
                                                   sequence=str(seq.get("state", "")), synteny=synteny,
                                                   domain=domain, has_sequence=has_sequence,
                                                   pangene_present=pangene_present),
            })
            audit_rows.append({"source_full_id": source_id, "target_full_id": target_id,
                               "audit_state": "mapped", "relationship": relation,
                               "mapping_state": state,
                               "reason": mapping_rows[-1]["mapping_reason"]})

    mapping = pd.DataFrame(mapping_rows, columns=MAPPING_COLUMNS).sort_values(["source_full_id", "target_full_id"], kind="stable")
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    mapping_path = out / "assembly_bridge_mapping.tsv"
    mapping.to_csv(mapping_path, sep="\t", index=False, lineterminator="\n")

    # Make the fail-closed adapter boundary visible in the audit rather than
    # silently dropping rows from the legacy Route B-shaped outputs.
    mapping_by_pair = {
        (str(row["source_full_id"]), str(row["target_full_id"])): row
        for row in mapping.to_dict("records")
    }
    for audit in audit_rows:
        pair = (str(audit["source_full_id"]), str(audit["target_full_id"]))
        mapped = mapping_by_pair.get(pair)
        if mapped and mapped["target_full_id"] and not _adapter_eligible(mapped)[0]:
            audit["audit_state"] = "adapter_excluded"
            audit["reason"] = f"{audit['reason']};route_b_adapter_excluded"

    function_by_gene: dict[str, list[dict[str, str]]] = {}
    for row in function_rows:
        key = str(row.get(function_gene_column, "")).strip()
        if key:
            function_by_gene.setdefault(key, []).append(row)
    orthology_rows: list[dict[str, str]] = []
    reference_rows: list[dict[str, str]] = []
    for map_row in mapping.to_dict("records"):
        source_id = str(map_row["source_full_id"])
        target_id = str(map_row["target_full_id"])
        if not target_id or source_id not in function_by_gene:
            continue
        relation = str(map_row["relationship"])
        state = str(map_row["mapping_state"])
        adapter_ok, _adapter_reason = _adapter_eligible(map_row)
        if not adapter_ok:
            continue
        route_relation = "one_to_one" if state == "verified" else relation
        if state == "provisional":
            route_relation = "family_level_only"
        orthology_rows.append({
            "gene_id": target_id, "reference_gene_id": source_id,
            "orthogroup": str(map_row["source_family_id"]), "relationship": route_relation,
            "provider": "pisum_assembly_bridge", "percent_identity": str(map_row["sequence_identity"]),
            "query_coverage": str(map_row["query_coverage"]), "evalue": "", "bitscore": "",
        })
        for function_row in function_by_gene[source_id]:
            functional_text = str(function_row.get("phenotype_or_function", "") or function_row.get("function_description", "")).strip()
            reference_rows.append({
                "reference_gene_id": source_id,
                "gene_id": target_id if state == "verified" and relation == "one_to_one" else "",
                "function_description": functional_text,
                "pathway": str(function_row.get("pathway", "")).strip(),
                "pathway_id": str(function_row.get("pathway_id", "")).strip(),
                "mutant_phenotype": str(function_row.get("mutant_phenotype", "")).strip(),
                "source_id": str(function_row.get("source_ids", "") or function_row.get("source_id", "")).strip(),
                "source_url": str(function_row.get("source_url", "")).strip(),
                "reference_species": _first_value(function_row, "reference_species", "native_species"),
                "gene_symbol": _first_value(function_row, "gene_symbol"),
                "evidence_type": _first_value(function_row, "evidence_type"),
                "provenance_detail": _bridge_provenance(function_row, map_row, source_id, target_id),
                "mapping_state": state,
                "mapping_relationship": relation,
                "function_transfer_scope": "exact_gene" if state == "verified" and relation == "one_to_one" else "family_set",
            })
    orthology_path = out / "route_b_orthology.tsv"
    reference_path = out / "route_b_reference_function.tsv"
    pd.DataFrame(orthology_rows, columns=ORTHOLOGY_COLUMNS).drop_duplicates().to_csv(orthology_path, sep="\t", index=False, lineterminator="\n")
    pd.DataFrame(reference_rows, columns=REFERENCE_FUNCTION_COLUMNS).drop_duplicates().to_csv(reference_path, sep="\t", index=False, lineterminator="\n")
    audit_path = out / "assembly_bridge_audit.tsv"
    pd.DataFrame(audit_rows, columns=AUDIT_COLUMNS).drop_duplicates().to_csv(audit_path, sep="\t", index=False, lineterminator="\n")

    counts = {
        "native_function_card_count": int(len(function_rows)),
        "mapping_attempted": int(mapping["target_full_id"].ne("").sum()),
        "mapping_verified_one_to_one": int((mapping["mapping_state"].eq("verified") & mapping["relationship"].eq("one_to_one")).sum()),
        "mapping_one_to_many": int(mapping["relationship"].eq("one_to_many").sum()),
        "mapping_unresolved": int(mapping["mapping_state"].eq("not_evaluated").sum()),
        "source_card_without_target_mapping": int(sum(1 for row in mapping_rows if not row["target_full_id"])),
        "family_concordant_count": int(sum(bool(str(row["source_family_id"]) and str(row["target_family_id"]) and set(str(row["source_family_id"]).split(";")) & set(str(row["target_family_id"]).split(";"))) for row in mapping_rows)),
        "reciprocal_sequence_verified_count": int(sum(_state_is_reciprocal(row["reciprocal_sequence_state"]) for row in mapping_rows)),
        "synteny_verified_count": int(sum(str(row["synteny_state"]).casefold() in {"verified", "concordant", "consistent"} for row in mapping_rows)),
        "pangene_candidate_mapping_count": int(mapping["pangene_relation_present"].eq("true").sum()),
        "pangene_one_to_one_candidate_count": int(mapping["pangene_relationship"].eq("one_to_one_candidate").sum()),
        "pangene_ambiguous_candidate_count": int(mapping["pangene_relationship"].isin({"one_to_many_candidate", "many_to_one_candidate", "many_to_many_candidate"}).sum()),
        "final_target_reference_count": int(len(reference_rows)),
        "route_b_adapter_mapping_count": int(sum(_adapter_eligible(row)[0] for row in mapping.to_dict("records") if str(row["target_full_id"]))),
        "route_b_adapter_excluded_mapping_count": int(sum((not _adapter_eligible(row)[0]) for row in mapping.to_dict("records") if str(row["target_full_id"]))),
    }
    input_pairs: list[tuple[str, str | Path]] = [
        ("functional_evidence", functional_evidence_arg), ("source_fam3", source_fam3_arg),
        ("target_fam3", target_fam3_arg), ("sequence_evidence", sequence_evidence_arg),
    ]
    if pangene_hash_arg:
        input_pairs.append(("pangene_hash", pangene_hash_arg))
    manifest = {
        "schema_version": SCHEMA_VERSION, "command": "assembly-bridge", "case_id": case_id,
        "mapping_scope": mapping_scope, "source_namespace": source_namespace,
        "target_namespace": target_namespace, "rank_contribution": "none",
        "pangene_dataset": pangene_dataset if pangene_hash_arg else None,
        "pangene_namespace_coverage": pangene_coverage if pangene_hash_arg else None,
        "inputs": {name: file_fingerprint(runtime_path(path), full_hash=True) for name, path in input_pairs},
        "outputs": {"mapping": str(mapping_path), "route_b_orthology": str(orthology_path),
                    "route_b_reference_function": str(reference_path), "audit": str(audit_path)},
        "counts": counts,
        "policy": [
            "full native IDs and namespaces are preserved; no symbol-only join is performed",
            "root functional-evidence identity and ordered mapping hops are retained in provenance_detail across chained bridges",
            "pangene membership is candidate correspondence only and never verifies one-copy identity",
            "when a cross-side pangene relation exists it narrows candidate generation; fam3 becomes a concordance check",
            "sources without a cross-side pangene candidate retain the previous fam3 fallback; absence is not negative evidence",
            "pangene split/merge multiplicity is retained in the full mapping table",
            "same fam3 alone is provisional and never verified",
            "only unique reciprocal sequence plus family concordance can verify a one-to-one relation",
            "family mismatch and synteny conflict remain provisional/ambiguous, never forced",
            "mapping is non-scoring and direct same-species provenance stays distinct from cross-species evidence",
        ],
        "fingerprints": {"mapping": file_fingerprint(mapping_path, full_hash=True),
                         "route_b_orthology": file_fingerprint(orthology_path, full_hash=True),
                         "route_b_reference_function": file_fingerprint(reference_path, full_hash=True),
                         "audit": file_fingerprint(audit_path, full_hash=True)},
    }
    manifest_path = out / "assembly_bridge_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"mapping": str(mapping_path), "route_b_orthology": str(orthology_path),
            "route_b_reference_function": str(reference_path), "audit": str(audit_path),
            "manifest": str(manifest_path), "counts": counts}


build_same_species_assembly_bridge = build_assembly_bridge
