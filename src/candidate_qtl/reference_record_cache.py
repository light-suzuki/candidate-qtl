"""Fetch a frozen UniProt record cache for Route A's named references.

The cache is an input to the existing offline reference-gene-resolve stage.
Neither a name hit nor a saved record verifies a biological claim.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import shutil
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .core import runtime_path
from .reference_gene_resolution import (
    _fetch_uniprot,
)
from .reference_gene_resolution_legacy import _resolve_reference_candidates
from .uniprot_error_policy import UNIPROT_ENTRY_ENDPOINT
from .reference_name_pairs import (IDENTITY_DECLARED_ACCESSION, IDENTITY_MISMATCH, atomic_reference_pairs,
                                   declared_record_identity, declared_uniprot_accessions,
                                   reference_lookup_handles, reference_identity_names, reference_name_for_input)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def validate_reference_record_cache(manifest_arg: str, *, mechanisms_arg: str | None = None,
                                    records_dir_arg: str | None = None) -> tuple[dict, list[dict]]:
    """Validate frozen input, ledger and exact record roster before reuse/retry."""
    path = runtime_path(manifest_arg).resolve()
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict):
        raise ValueError("reference record cache manifest must be an object")
    if manifest.get("contract") != "route-a-reference-record-cache-1.0.0":
        raise ValueError("unsupported reference record cache contract")
    mechanisms = runtime_path(mechanisms_arg or manifest["mechanisms"]).resolve()
    if _sha(mechanisms) != manifest.get("mechanisms_sha256"):
        raise ValueError("reference cache mechanisms hash mismatch")
    ledger_path = runtime_path(manifest["query_ledger"]).resolve()
    if _sha(ledger_path) != manifest.get("query_ledger_sha256"):
        raise ValueError("reference cache query ledger hash mismatch")
    ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
    expected = manifest.get("records_sha256", {})
    if not isinstance(ledger, list) or any(not isinstance(item, dict) for item in ledger) or not isinstance(expected, dict):
        raise ValueError("reference cache malformed ledger/record roster")
    if any(not re.fullmatch(r"[A-Z0-9]+", accession) or not re.fullmatch(r"[0-9a-f]{64}", str(digest))
           for accession, digest in expected.items()):
        raise ValueError("reference cache malformed accession/hash")
    records_dir = runtime_path(records_dir_arg or manifest.get("records_dir") or manifest.get("partial_records_dir")).resolve()
    if not records_dir.is_dir() or (records_dir_arg and manifest.get("records_dir")
                                   and records_dir != runtime_path(manifest["records_dir"]).resolve()):
        raise ValueError("reference cache record directory identity mismatch")
    actual = {p.stem: _sha(p) for p in records_dir.glob("*.json")}
    if actual != expected:
        raise ValueError("reference cache record roster/hash mismatch")
    if any(not isinstance(item.get("reference_gene_name"), str) or not isinstance(item.get("reference_species"), str)
           for item in ledger):
        raise ValueError("reference cache missing query identity")
    keys = [(item["reference_gene_name"], item["reference_species"]) for item in ledger]
    if len(keys) != len(set(keys)) or len(keys) != manifest.get("named_query_count"):
        raise ValueError("reference cache duplicate/missing query identity")
    with mechanisms.open(encoding="utf-8", newline="") as stream:
        input_rows = list(csv.DictReader(stream, delimiter="\t"))
    input_ids = [(row.get("case_id", ""), row.get("hypothesis_id", "")) for row in input_rows]
    if (any(not all(pair) or any(value != value.strip() for value in pair) for pair in input_ids)
            or len(input_ids) != len({tuple(value.casefold() for value in pair) for pair in input_ids})):
        raise ValueError("reference cache input duplicate/malformed case/hypothesis identity")
    grouped = {}
    raw_inputs_by_pair: dict[tuple[str, str], list[dict[str, str]]] = {}
    unresolved_keys = []
    for row in input_rows:
        if any(not row.get(key) for key in ("case_id", "hypothesis_id")):
            raise ValueError("reference cache input missing case/hypothesis identity")
        pairs, _ = atomic_reference_pairs(reference_name_for_input(row), row.get("reference_species", ""))
        if not pairs:
            unresolved_keys.append((row["case_id"], row["hypothesis_id"]))
        for pair in pairs:
            grouped.setdefault(pair, []).append(f"{row['case_id']}/{row['hypothesis_id']}")
            raw_inputs_by_pair.setdefault(pair, []).append(row)
    if set(keys) != set(grouped):
        raise ValueError("reference cache query/input identity roster mismatch")
    if "cache_complete" in manifest:
        unresolved = manifest.get("unresolved_inputs", [])
        if (not isinstance(unresolved, list) or any(not isinstance(row, dict) for row in unresolved)
                or sorted((row.get("case_id"), row.get("hypothesis_id")) for row in unresolved) != sorted(unresolved_keys)):
            raise ValueError("reference cache unresolved input identity mismatch")
        complete = not unresolved_keys and all(item.get("status") in {"complete", "no_provider_match"} for item in ledger)
        if manifest["cache_complete"] is not complete:
            raise ValueError("reference cache completion state contradicts query/input roster")
    for item in ledger:
        if item.get("status") not in {"complete", "no_provider_match", "not_evaluated"}:
            raise ValueError("reference cache invalid query status")
        if (not isinstance(item.get("query_urls", []), list)
                or not isinstance(item.get("lookup_failures", []), list)
                or not isinstance(item.get("record_fetch_failures", {}), dict)
                or not isinstance(item.get("record_identity_failures", {}), dict)):
            raise ValueError("reference cache malformed lookup provenance")
        if any(not isinstance(item.get(key, []), list) for key in ("cached_accessions", "accessions", "hypothesis_keys")):
            raise ValueError("reference cache malformed accession/identity roster")
        if not isinstance(item.get("provisional_accessions", []), list) or not set(
                item.get("provisional_accessions", [])).issubset(set(item.get("cached_accessions", []))):
            raise ValueError("reference cache provisional accessions are not cached records")
        if item.get("provisional_accessions") and set(item["cached_accessions"]) - set(item["provisional_accessions"]):
            raise ValueError("reference cache keeps a provisional accession beside an exact record")
        if not set(item.get("cached_accessions", [])).issubset(expected):
            raise ValueError("reference cache ledger record identity mismatch")
        if item.get("hypothesis_keys") != sorted(grouped[(item["reference_gene_name"], item["reference_species"])]):
            raise ValueError("reference cache case/hypothesis identity mismatch")
        if item["status"] == "complete" and (sorted(item.get("cached_accessions", [])) != item.get("accessions")
                                              or len(item.get("accessions", [])) != item.get("provider_hit_count")):
            raise ValueError("reference cache incomplete query declared complete")
        for accession in item.get("cached_accessions", []):
            record = json.loads((records_dir / f"{accession}.json").read_text(encoding="utf-8"))
            item_inputs = item.get("raw_inputs") or raw_inputs_by_pair[(item["reference_gene_name"], item["reference_species"])]
            identity_names = reference_identity_names(item_inputs, item["reference_gene_name"], item["reference_species"])
            identity = declared_record_identity(record, accession, identity_names, item["reference_species"],
                                                declared_uniprot_accessions(item_inputs))
            # A provisional (declared accession + species, name unconfirmed)
            # record is cached only when the ledger lists it as provisional.
            if identity == IDENTITY_MISMATCH or (
                    identity == IDENTITY_DECLARED_ACCESSION
                    and accession not in item.get("provisional_accessions", [])):
                raise ValueError("reference cache record name/species/primary identity mismatch")
    return manifest, ledger


def build_reference_record_cache(
    mechanisms_arg: str,
    out_dir_arg: str,
    *,
    allow_external_api: bool = False,
    timeout: int = 60,
    max_queries: int = 32,
    max_records: int = 128,
    max_elapsed_seconds: int = 1800,
    resume_manifest_arg: str | None = None,
) -> dict[str, Any]:
    """Retrieve every exact provider alternative, with explicit finite limits.

    Healthy records remain available during a failed query. Missing endpoints
    stay not_evaluated; retry creates a new checkpoint and reissues failed queries.
    """
    if not allow_external_api:
        raise ValueError("record retrieval requires explicit --allow-external-api")
    if min(timeout, max_queries, max_records, max_elapsed_seconds) < 1:
        raise ValueError("timeout and retrieval limits must be positive")
    mechanisms = runtime_path(mechanisms_arg).resolve()
    output = runtime_path(out_dir_arg).resolve()
    if output.exists():
        raise ValueError("output directory must not already exist")
    with mechanisms.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"hypothesis_id", "case_id", "reference_gene_name", "reference_species"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("mechanisms lack reference-name identity columns")
        rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError("malformed mechanism row")
    ids = [(row["case_id"], row["hypothesis_id"]) for row in rows]
    if (any(not all(pair) or any(value != value.strip() for value in pair) for pair in ids)
            or len(ids) != len({tuple(value.casefold() for value in pair) for pair in ids})):
        raise ValueError("mechanism IDs must be unique and nonempty per case")
    grouped: dict[tuple[str, str], list[str]] = {}
    raw_inputs: dict[tuple[str, str], list[dict[str, str]]] = {}
    unresolved_inputs: list[dict[str, str]] = []
    for row in rows:
        name, species = reference_name_for_input(row), row["reference_species"].strip()
        pairs, reason = atomic_reference_pairs(name, species)
        original = dict(row)
        if not pairs:
            unresolved_inputs.append({**original, "status": "not_evaluated", "reason": reason})
        for pair in pairs:
            grouped.setdefault(pair, []).append(f"{row['case_id']}/{row['hypothesis_id']}")
            raw_inputs.setdefault(pair, []).append(original)
    if len(grouped) > max_queries:
        raise ValueError(f"{len(grouped)} unique named references exceed max_queries={max_queries}")

    previous, previous_ledger = ({}, [])
    if resume_manifest_arg:
        previous, previous_ledger = validate_reference_record_cache(resume_manifest_arg, mechanisms_arg=str(mechanisms))
        if {(item["reference_gene_name"], item["reference_species"]) for item in previous_ledger} != set(grouped):
            raise ValueError("resume query identity roster mismatch")
    output.mkdir(parents=True)
    records_dir = output / "records"
    records_dir.mkdir()
    deadline = time.monotonic() + max_elapsed_seconds
    ledger: list[dict[str, Any]] = []
    records: dict[str, str] = {}
    if previous:
        old_records = Path(previous.get("records_dir") or previous["partial_records_dir"])
        for accession, digest in previous["records_sha256"].items():
            shutil.copyfile(old_records / f"{accession}.json", records_dir / f"{accession}.json")
            records[accession] = digest
    previous_by_pair = {(item["reference_gene_name"], item["reference_species"]): item for item in previous_ledger}
    for name, species in sorted(grouped, key=lambda pair: (pair[1].casefold(), pair[0].casefold())):
        entry: dict[str, Any] = {
            "reference_gene_name": name, "reference_species": species,
            "hypothesis_keys": sorted(grouped[(name, species)]),
            "query": "", "provider_hit_count": 0, "accessions": [],
            "status": "not_evaluated", "reason": "", "logical_query_attempts": 0,
            "retrieved_at": "", "search_url": "", "query_urls": [], "record_urls": {},
            "raw_inputs": raw_inputs[(name, species)], "cached_accessions": [],
            "lookup_handles": reference_lookup_handles(raw_inputs[(name, species)]),
            "lookup_failures": [], "record_fetch_failures": {}, "record_identity_failures": {},
        }
        prior = previous_by_pair.get((name, species))
        if prior and prior["status"] in {"complete", "no_provider_match"}:
            ledger.append({**prior, "logical_query_attempts": 0, "reused_from_manifest": str(resume_manifest_arg)})
            continue
        if time.monotonic() >= deadline:
            entry["reason"] = "elapsed_budget_exhausted"
            ledger.append(entry)
            continue
        try:
            entry["logical_query_attempts"] = 1
            entry["retrieved_at"] = datetime.now(timezone.utc).isoformat()
            accessions, status, query, hits, lookup_failures, query_urls = _resolve_reference_candidates(
                name, species, raw_inputs[(name, species)], timeout)
            entry.update(query=query, query_urls=query_urls, lookup_failures=lookup_failures, provider_hit_count=hits,
                         accessions=sorted(set(accessions)))
            if len(query_urls) == 1:
                entry["search_url"] = query_urls[0]
            if hits != len(entry["accessions"]):
                raise ValueError("provider hit count does not match the complete accession roster")
            if not accessions:
                entry.update(status="not_evaluated" if lookup_failures else
                             "no_provider_match" if status == "unresolved_gene_name" else "not_evaluated",
                             reason=json.dumps(lookup_failures, sort_keys=True) if lookup_failures else status)
                ledger.append(entry)
                continue
            if len(records | {accession: "" for accession in accessions}) > max_records:
                entry["reason"] = "record_budget_exceeded"
                ledger.append(entry)
                continue
            fetched: list[str] = []
            entry["cached_accessions"] = fetched
            failures = {}
            identity_failures = {}
            provisional_records: dict[str, dict] = {}
            for accession in entry["accessions"]:
                try:
                    if time.monotonic() >= deadline:
                        raise TimeoutError("elapsed_budget_exhausted")
                    if not re.fullmatch(r"[A-Z0-9]+", accession):
                        raise ValueError("provider returned an unsafe accession")
                    entry["record_urls"][accession] = UNIPROT_ENTRY_ENDPOINT + accession + ".json"
                    if accession in records:
                        record = json.loads((records_dir / f"{accession}.json").read_text(encoding="utf-8"))
                    else:
                        record = _fetch_uniprot(accession, timeout)
                    if not isinstance(record, dict):
                        raise ValueError(f"provider returned an invalid record: {accession}")
                    if record.get("primaryAccession") != accession:
                        identity_failures[accession] = f"provider record identity mismatch: {accession}"
                        continue
                    identity_names = reference_identity_names(raw_inputs[(name, species)], name, species)
                    identity = declared_record_identity(record, accession, identity_names, species,
                                                        declared_uniprot_accessions(raw_inputs[(name, species)]))
                    if identity == IDENTITY_MISMATCH:
                        identity_failures[accession] = f"provider record name/species mismatch: {accession}"
                        continue
                    if identity == IDENTITY_DECLARED_ACCESSION:
                        # The worker's own accession exists with the declared
                        # species; only the name is unconfirmed. Hold it until
                        # every alternative is checked (see below).
                        provisional_records[accession] = record
                        continue
                    if accession not in records:
                        path = records_dir / f"{accession}.json"
                        _write_json(path, record)
                        records[accession] = _sha(path)
                    fetched.append(accession)
                except (RuntimeError, TimeoutError, OSError, ValueError) as error:
                    failures[accession] = f"{type(error).__name__}: {error}"
            if provisional_records and fetched:
                # An exact name/species record exists: the declared accession
                # names another gene and is not adopted, even provisionally.
                for accession in provisional_records:
                    identity_failures[accession] = f"provider record name/species mismatch: {accession}"
            elif provisional_records:
                # Keep the record as provisional review material, never support.
                for accession, record in sorted(provisional_records.items()):
                    if accession not in records:
                        path = records_dir / f"{accession}.json"
                        _write_json(path, record)
                        records[accession] = _sha(path)
                    fetched.append(accession)
                entry["provisional_accessions"] = sorted(provisional_records)
            all_failures = {f"lookup:{index}": failure for index, failure in enumerate(lookup_failures)}
            all_failures.update(failures)
            all_failures.update(identity_failures)
            entry["failed_accessions"] = {**failures, **identity_failures}
            entry["record_fetch_failures"] = failures
            entry["record_identity_failures"] = identity_failures
            entry.update(status="not_evaluated" if all_failures else "complete",
                         reason=json.dumps(all_failures, sort_keys=True) if all_failures else status)
        except (RuntimeError, TimeoutError, OSError, ValueError) as error:
            entry.update(status="not_evaluated", reason=f"{type(error).__name__}: {error}")
        ledger.append(entry)

    ledger_path = output / "query_ledger.json"
    _write_json(ledger_path, ledger)
    complete = not unresolved_inputs and all(item["status"] in {"complete", "no_provider_match"} for item in ledger)
    manifest = {
        "contract": "route-a-reference-record-cache-1.0.0",
        "mechanisms": str(mechanisms), "mechanisms_sha256": _sha(mechanisms),
        "records_dir": str(records_dir), "partial_records_dir": "",
        "records_sha256": dict(sorted(records.items())),
        "query_ledger": str(ledger_path), "query_ledger_sha256": _sha(ledger_path),
        "named_query_count": len(grouped), "record_count": len(records),
        "logical_query_attempt_count": sum(item["logical_query_attempts"] for item in ledger),
        "status_counts": {status: sum(item["status"] == status for item in ledger)
                          for status in ("complete", "no_provider_match", "not_evaluated")},
        "cache_usable_for_prepare": True, "cache_complete": complete,
        "unresolved_inputs": unresolved_inputs,
        "failed_query_count": sum(item["status"] == "not_evaluated" for item in ledger),
        "resume_manifest": str(runtime_path(resume_manifest_arg).resolve()) if resume_manifest_arg else None,
        "resume_manifest_sha256": _sha(runtime_path(resume_manifest_arg)) if resume_manifest_arg else None,
        "limits": {"timeout": timeout, "max_queries": max_queries,
                   "max_records": max_records, "max_elapsed_seconds": max_elapsed_seconds},
        "interpretation": "Record identity context only; source/claim and sequence gates remain separate",
    }
    manifest_path = output / "manifest.json"
    _write_json(manifest_path, manifest)
    return {**manifest, "manifest": str(manifest_path), "manifest_sha256": _sha(manifest_path)}
