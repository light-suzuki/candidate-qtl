"""Fail-closed import of Route-B semantic LLM proposals (#403).

The common ``llm-run`` command is the provider adapter. This module maps its
normalized responses back to the frozen Route-B semantic units and converts
model JSON into the shared semantic-claim schema.

This is still a *proposal* stage. It can mechanically reject invented source
identifiers or reference genes, but it cannot establish that a paper actually
supports the asserted biological relation. Therefore well-formed, in-bundle
claims remain ``proposed`` until later source-content and transfer verification.
"""

from __future__ import annotations

import hashlib
import json

import pandas as pd

from .core import runtime_path
from .route_b_response_validation import route_b_claim_schema_error
from .semantic_claims import CLAIM_COLUMNS, EVIDENCE_SCOPES, normalize_semantic_claims

SCHEMA_VERSION = "1.1"
STATUS_COLUMNS = [
    "unit_id", "prompt_id", "parse_status", "model_claim_count",
    "imported_claim_row_count", "unsupported_claim_row_count", "rejected_claim_row_count", "status", "note",
]


def _sha256_json(value: object) -> str:
    encoded = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


_ARTIFACT_HASH_KEYS = ("bundle_sha256", "prompt_sha256", "response_sha256")


def _sha256_file(path_arg: str) -> str:
    path = runtime_path(path_arg)
    if not path.exists():
        raise ValueError(f"semantic provenance artifact does not exist: {path}")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_route_b_semantic_artifacts(
    bundle_arg: str,
    prompts_arg: str,
    responses_arg: str,
    llm_manifest_arg: str | None,
    *,
    require_manifest: bool = False,
) -> dict[str, str] | None:
    """Validate frozen bundle/prompt/response bytes against an LLM manifest.

    A direct proposal import may remain usable without a manifest because it
    only emits non-verified proposals.  The canonical extension passes
    ``require_manifest=True`` before it can consume any response as a frozen
    production input.  If a manifest is supplied, all three artifact hashes
    are mandatory; provider/model-only legacy metadata is never accepted as
    frozen provenance.
    """
    if not llm_manifest_arg:
        if require_manifest:
            raise ValueError(
                "canonical Route B semantic responses require llm_run_manifest.json "
                "with bundle_sha256, prompt_sha256, and response_sha256"
            )
        return None

    manifest_path = runtime_path(llm_manifest_arg)
    if not manifest_path.exists():
        raise ValueError(f"LLM run manifest does not exist: {manifest_path}")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"invalid LLM run manifest: {manifest_path}") from error
    if not isinstance(manifest, dict):
        raise ValueError(f"LLM run manifest must be a JSON object: {manifest_path}")

    missing = [
        key for key in _ARTIFACT_HASH_KEYS
        if not isinstance(manifest.get(key), str) or not manifest[key].strip()
    ]
    if missing:
        raise ValueError(
            "LLM run manifest is missing frozen artifact hashes: "
            + ", ".join(missing)
        )

    current = {
        "bundle_sha256": _sha256_file(bundle_arg),
        "prompt_sha256": _sha256_file(prompts_arg),
        "response_sha256": _sha256_file(responses_arg),
    }
    mismatches = [
        f"{key}: manifest={manifest[key]} current={current[key]}"
        for key in _ARTIFACT_HASH_KEYS
        if manifest[key].strip().casefold() != current[key]
    ]
    if mismatches:
        raise ValueError(
            "LLM run manifest frozen artifact hash mismatch: "
            + "; ".join(mismatches)
        )
    return current


_EXECUTION_EVIDENCE_KEYS = (
    "distinct_response_count",
    "empty_claim_response_count",
    "execution_evidence",
    "counts_as_per_prompt_evaluation",
)


def _llm_execution_evidence(llm_manifest_arg: str | None) -> dict[str, object] | None:
    """Carry llm-run's per-prompt-variation verdict onto the proposal import.

    The import manifest is the first artifact that can be read as "Route B ran
    for this case", so the evidence has to travel with it. Without this the
    verdict lives only in the llm-run manifest, which every later consumer --
    the Wet handoff, a scorer, a report author -- would have to go find.
    """
    if not llm_manifest_arg:
        return None
    manifest_path = runtime_path(llm_manifest_arg)
    if not manifest_path.exists():
        return None
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(manifest, dict):
        return None
    evidence = {key: manifest[key] for key in _EXECUTION_EVIDENCE_KEYS if key in manifest}
    return evidence or None


def _reference_tokens(value: object) -> list[str]:
    """Split a reference-gene field into the identifiers it actually names.

    The atlas joins multiple orthology references with ``;``, so the value a
    model echoes back can name several genes. Order is preserved and blanks
    dropped; nothing is normalised beyond whitespace, because grounding is an
    exact-membership test against the frozen bundle.
    """
    text = str(value or "")
    return [token.strip() for token in text.replace(",", ";").split(";") if token.strip()]


def _source_ids(value: object) -> list[str]:
    if isinstance(value, list):
        raw = value
    else:
        raw = str(value or "").replace(",", ";").split(";")
    return sorted({str(item).strip() for item in raw if str(item).strip()})


def _prompt_records(path_arg: str) -> dict[str, dict[str, object]]:
    path = runtime_path(path_arg)
    if not path.exists():
        raise ValueError(f"Route B semantic prompts do not exist: {path}")
    records: dict[str, dict[str, object]] = {}
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        record = json.loads(line)
        if not isinstance(record, dict) or not record.get("prompt_id") or not record.get("route_b_unit_id"):
            raise ValueError(f"invalid Route B semantic prompt at line {line_number}")
        prompt_id = str(record["prompt_id"])
        if prompt_id in records:
            raise ValueError(f"duplicate Route B semantic prompt_id: {prompt_id}")
        records[prompt_id] = record
    return records


def _semantic_claim_row_errors(rows: list[dict[str, str]]) -> dict[int, str]:
    """Retain single-row quarantine errors while validating valid units once.

    The shared normalizer remains the only validation policy. If a unit batch
    fails, validate its rows individually to preserve each original error and
    avoid rejecting valid neighbors. Cross-row checks still run on the final
    accepted table in ``import_route_b_semantic_proposals``.
    """
    if not rows:
        return {}
    try:
        normalize_semantic_claims(pd.DataFrame(rows, columns=CLAIM_COLUMNS))
    except ValueError as error:
        if len(rows) == 1:
            return {0: str(error)}
        errors = {}
        for index, row in enumerate(rows):
            try:
                normalize_semantic_claims(pd.DataFrame([row], columns=CLAIM_COLUMNS))
            except ValueError as error:
                errors[index] = str(error)
        return errors
    return {}


def import_route_b_semantic_proposals(
    bundle_arg: str,
    prompts_arg: str,
    responses_arg: str,
    out_dir_arg: str,
    *,
    llm_manifest_arg: str | None = None,
    case_id: str = "",
) -> dict[str, object]:
    """Import normalized ``llm-run`` responses as non-scoring claim proposals."""
    bundle_path = runtime_path(bundle_arg)
    responses_path = runtime_path(responses_arg)
    if not bundle_path.exists():
        raise ValueError(f"Route B semantic bundle does not exist: {bundle_path}")
    if not responses_path.exists():
        raise ValueError(f"LLM normalized responses do not exist: {responses_path}")

    bundle = json.loads(bundle_path.read_text(encoding="utf-8"))
    if not isinstance(bundle, dict) or bundle.get("route") != "interval_first" or not isinstance(bundle.get("units"), list):
        raise ValueError("invalid Route B semantic bundle")
    units = {str(unit.get("unit_id", "")): unit for unit in bundle["units"] if isinstance(unit, dict)}
    if units and "" in units:
        raise ValueError("Route B semantic bundle requires non-empty unique unit_id values")
    if len(units) != len(bundle["units"]):
        raise ValueError("Route B semantic bundle contains duplicate unit_id values")

    # A unit that carries no reviewable material is deliberately not sent to a
    # provider (issue #496 Phase 2). It stays in the bundle and is imported as
    # an explicit not_evaluated row, so the roster check runs against the units
    # that were actually requested rather than silently accepting a gap.
    not_requested = {
        unit_id for unit_id, unit in units.items()
        if str(unit.get("request_state", "requested")) != "requested"
    }
    requested_units = {unit_id for unit_id in units if unit_id not in not_requested}
    prompts = _prompt_records(prompts_arg)
    prompt_to_unit = {prompt_id: str(record["route_b_unit_id"]) for prompt_id, record in prompts.items()}
    if set(prompt_to_unit.values()) != requested_units:
        missing = sorted(requested_units - set(prompt_to_unit.values()))
        extra = sorted(set(prompt_to_unit.values()) - requested_units)
        raise ValueError(f"prompt/bundle unit mismatch: missing={missing}, extra={extra}")

    payload = json.loads(responses_path.read_text(encoding="utf-8"))
    responses = payload.get("candidates", []) if isinstance(payload, dict) else []
    if not isinstance(responses, list):
        raise ValueError("llm-run normalized responses require candidates list")
    by_prompt: dict[str, dict[str, object]] = {}
    for response in responses:
        if not isinstance(response, dict) or not response.get("prompt_id"):
            raise ValueError("every semantic response requires prompt_id")
        prompt_id = str(response["prompt_id"])
        if prompt_id in by_prompt:
            raise ValueError(f"duplicate semantic response prompt_id: {prompt_id}")
        by_prompt[prompt_id] = response
    if set(by_prompt) != set(prompts):
        missing = sorted(set(prompts) - set(by_prompt))
        extra = sorted(set(by_prompt) - set(prompts))
        raise ValueError(f"semantic response roster mismatch: missing={missing}, extra={extra}")

    artifact_hashes = validate_route_b_semantic_artifacts(
        str(bundle_path), prompts_arg, str(responses_path), llm_manifest_arg,
    )
    provider = ""
    model = ""
    llm_manifest: dict[str, object] | None = None
    if llm_manifest_arg:
        llm_manifest_path = runtime_path(llm_manifest_arg)
        llm_manifest = json.loads(llm_manifest_path.read_text(encoding="utf-8"))
        if isinstance(llm_manifest, dict):
            provider = str(llm_manifest.get("provider", ""))
            model = str(llm_manifest.get("model", ""))

    rows: list[dict[str, str]] = []
    rejected_rows: list[dict[str, object]] = []
    status_rows: list[dict[str, object]] = []
    for prompt_id in sorted(prompts):
        prompt = prompts[prompt_id]
        response = by_prompt[prompt_id]
        unit_id = prompt_to_unit[prompt_id]
        unit = units[unit_id]
        identity = unit.get("identity", {}) if isinstance(unit.get("identity"), dict) else {}
        family_ids = [str(item) for item in identity.get("evolutionary_family_ids", []) if str(item).strip()]
        if not family_ids:
            family_ids = [str(item) for item in identity.get("family_ids", []) if str(item).strip()]
        orthogroups = [str(item) for item in identity.get("orthogroup_ids", []) if str(item).strip()]
        subfamilies = [str(item) for item in identity.get("subfamily_ids", []) if str(item).strip()]
        members = sorted({str(item).strip() for item in unit.get("member_interval_genes", []) if str(item).strip()})
        known_sources = {str(item).strip() for item in unit.get("stable_source_ids", []) if str(item).strip()}
        # The roster is tokenised for the same reason the claim side is: the atlas
        # stores several orthology references in one semicolon-joined field, and a
        # claim naming one of them must be groundable against it. The
        # joined form is kept as well, so a claim echoing the field verbatim still
        # matches and no previously accepted input becomes ungrounded.
        known_reference_genes = set()
        for record in unit.get("genes", []):
            if not isinstance(record, dict):
                continue
            joined = str(record.get("orthology_reference_gene", "")).strip()
            if joined:
                known_reference_genes.add(joined)
                known_reference_genes.update(_reference_tokens(joined))
            # A family-index subject endpoint is an exact frozen retrieval
            # handle. It may ground a proposal reference, but says nothing
            # about orthology or trait support by itself.
            family_endpoints = str(record.get("family_index_subject_endpoints", "")).strip()
            if family_endpoints:
                known_reference_genes.add(family_endpoints)
                known_reference_genes.update(_reference_tokens(family_endpoints))
            # All acquired BLAST alternatives are frozen retrieval handles,
            # including hits absent from the atlas's family/orthology fields.
            # Grounding one never establishes orthology or source support.
            for hit in record.get("sequence_similarity_hits", []):
                if isinstance(hit, dict) and str(hit.get("hit_accession", "")).strip():
                    known_reference_genes.add(str(hit["hit_accession"]).strip())
        known_reference_genes.update(
            str(record.get("reference_gene", "")).strip()
            for record in unit.get("literature_link_candidates", [])
            if isinstance(record, dict) and str(record.get("reference_gene", "")).strip()
        )
        known_reference_genes.update(
            str(record.get("reference_gene", "")).strip()
            for record in unit.get("reference_material_cards", [])
            if isinstance(record, dict) and str(record.get("reference_gene", "")).strip()
        )
        from .observation_tiers import gene_endpoint_key
        from .paper_reports import source_tokens
        frozen_references = {}
        for reference in sorted(known_reference_genes):
            if ";" not in reference and gene_endpoint_key(reference):
                frozen_references.setdefault(gene_endpoint_key(reference), reference)
        frozen_sources = {}
        for source in sorted(known_sources):
            for key in source_tokens([source]):
                frozen_sources.setdefault(key, source)

        parse_status = str(response.get("parse_status", ""))
        model_json = response.get("model_json")
        if parse_status != "json" and (
            (isinstance(model_json, dict) and model_json.get("claims"))
            or response.get("claims")
        ):
            raise ValueError(
                f"Route B response {prompt_id} has claims with non-json parse_status={parse_status!r}"
            )
        claims = model_json.get("claims", []) if parse_status == "json" and isinstance(model_json, dict) else []
        # Merge quarantines an all-invalid response as not_evaluated. Preserve
        # its invalid claims in the importer ledger, without making them eligible.
        pre_rejected = response.get("rejected_claims", []) if parse_status == "not_evaluated" else []
        if not isinstance(claims, list):
            claims = []
        imported = 0
        unsupported = 0
        rejected = 0
        for rejected_claim in pre_rejected:
            rejected_rows.append({"unit_id": unit_id, "prompt_id": prompt_id,
                                  "claim_index": rejected_claim["claim_index"], "source_index": 0,
                                  "reason": rejected_claim["reason"], "raw_claim": rejected_claim["raw_claim"]})
            rejected += 1
        claim_rejected_start = len(rejected_rows)
        mapped_rows: list[tuple[dict[str, str], int, int, dict[str, object]]] = []
        fingerprints: dict[str, str] | None = None
        for claim_index, claim in enumerate(claims, 1):
            schema_error = route_b_claim_schema_error(claim)
            if schema_error:
                rejected_rows.append({"unit_id": unit_id, "prompt_id": prompt_id,
                                      "claim_index": claim_index, "source_index": 0,
                                      "reason": schema_error, "raw_claim": claim})
                rejected += 1
                continue
            claim_type = str(claim.get("claim_type", "")).strip()
            claim_members = members
            if "member_interval_genes" in claim:
                declared = claim["member_interval_genes"]
                if isinstance(declared, str):
                    declared = declared.replace(",", ";").split(";")
                valid = isinstance(declared, list) and all(isinstance(value, str) and value.strip() for value in declared)
                claim_members = sorted({value.strip() for value in declared}) if valid else []
                if not claim_members or not set(claim_members) <= set(members):
                    rejected_rows.append({"unit_id": unit_id, "prompt_id": prompt_id,
                                          "claim_index": claim_index, "source_index": 0,
                                          "reason": "member_interval_genes must be a non-empty subset of the frozen semantic unit",
                                          "raw_claim": claim})
                    rejected += 1
                    continue
            evidence_scope = str(claim.get("evidence_scope", "unknown")).strip() or "unknown"
            if evidence_scope not in EVIDENCE_SCOPES:
                evidence_scope = "unknown"
            reference_gene = str(claim.get("reference_gene", "")).strip()
            cited_sources = _source_ids(claim.get("source_ids", []))
            # Different spellings of the same frozen ID reuse its exact material
            # key. Raw responses/hashes remain unchanged; new IDs stay unsupported.
            cited_sources = list(dict.fromkeys(frozen_sources.get(next(iter(source_tokens([source]))), source) for source in cited_sources))
            source_rows = cited_sources or [""]
            # The bundle hands the model ``orthology_reference_gene`` as a
            # semicolon-joined field, so a one-to-many family arrives as
            # "sp|P93003|TFL1_ARATH;sp|Q9ZNV5|CEN_ARATH". A model that obeys
            # "use only the supplied unit" copies it back whole, and matching the
            # joined string against the frozen roster found nothing -- which
            # rejected precisely the paralog-rich families this lane exists to
            # handle. Ground each token instead. The closed world is
            # unchanged: a token absent from the bundle still grounds nothing, and
            # a claim whose tokens are all absent is still unsupported.
            reference_tokens = _reference_tokens(reference_gene)
            reference_tokens = [token if token in known_reference_genes else frozen_references.get(gene_endpoint_key(token), token)
                                for token in reference_tokens]
            if reference_tokens:
                reference_gene = ";".join(reference_tokens)
            grounded_references = [token for token in reference_tokens if token in known_reference_genes]
            reference_unknown = bool(reference_tokens) and not grounded_references

            for source_index, source_id in enumerate(source_rows, 1):
                problems: list[str] = []
                if source_id and source_id not in known_sources:
                    problems.append(f"source_id_not_in_frozen_bundle:{source_id}")
                if reference_unknown:
                    problems.append(f"reference_gene_not_in_frozen_bundle:{reference_gene}")
                grounded_note = (
                    f"reference_gene_grounded:{';'.join(grounded_references)}"
                    if grounded_references else ""
                )
                verification_status = "unsupported" if problems else "proposed"
                verification_reason = (
                    ";".join(problems)
                    if problems
                    else "; ".join(part for part in (
                        grounded_note,
                        "source_id_present_in_frozen_bundle; semantic_relation_verification_pending"
                        if source_id
                        else "no_source_id_supplied; semantic_relation_verification_pending",
                    ) if part)
                )

                if fingerprints is None:
                    fingerprints = {
                        "prompt_fingerprint": _sha256_json(prompt),
                        "input_fingerprint": _sha256_json(unit),
                        "output_fingerprint": _sha256_json(model_json if isinstance(model_json, dict) else response.get("raw_text", "")),
                    }
                row = {column: "" for column in CLAIM_COLUMNS}
                from .relevance_tiers import INITIAL_FIELDS
                row.update({key: str(claim.get(key, "")).strip() for key in INITIAL_FIELDS})
                if bundle.get("b_review_membership_policy") == "initial_bibliography_v1":
                    row.update(candidate_admission_policy="initial_bibliography_v1", focal_species=str(bundle.get("focal_species", "")))
                row.update({
                    "claim_id": f"{unit_id}-c{claim_index:03d}-s{source_index:02d}",
                    "route": "interval_first",
                    "semantic_unit_id": unit_id,
                    "family_id": family_ids[0] if family_ids else "",
                    "orthogroup_id": orthogroups[0] if orthogroups else "",
                    "subfamily_or_clade": subfamilies[0] if subfamilies else "",
                    "member_interval_genes": ";".join(claim_members),
                    "claim_type": claim_type,
                    "biological_function": str(claim.get("biological_function", "")).strip(),
                    "source_title": str(claim.get("source_title", "")).strip(),
                    "source_identifier_aliases": str(claim.get("source_identifier_aliases", "")).strip(),
                    "mechanism": str(claim.get("mechanism", "")).strip(),
                    "phenotype_relation": str(claim.get("phenotype_relation", "")).strip(),
                    "reference_gene": reference_gene,
                    "reference_species": str(claim.get("reference_species", "")).strip(),
                    "evidence_scope": evidence_scope,
                    "source_id": source_id,
                    "source_url": "",
                    "alternative_hypothesis": str(claim.get("alternative_hypothesis", "")).strip(),
                    "uncertainty_note": str(claim.get("uncertainty_note", "")).strip(),
                    "provider": provider,
                    "model": model,
                    **fingerprints,
                    "verification_status": verification_status,
                    "verification_reason": verification_reason,
                    "rank_contribution": "none",
                })
                from .paper_reports import imported_paper_fields
                row.update(imported_paper_fields(claim, source_id))
                mapped_rows.append((row, claim_index, source_index, claim))

        row_errors = _semantic_claim_row_errors([item[0] for item in mapped_rows])
        for row_index, (row, claim_index, source_index, claim) in enumerate(mapped_rows):
            if row_index in row_errors:
                rejected_rows.append({"unit_id": unit_id, "prompt_id": prompt_id,
                                      "claim_index": claim_index, "source_index": source_index,
                                      "reason": row_errors[row_index], "raw_claim": claim, "mapped_row": row})
                rejected += 1
                continue
            rows.append(row)
            imported += 1
            unsupported += int(row["verification_status"] == "unsupported")
        # Schema errors were collected before batched semantic validation.
        # Restore the original claim/source order without reordering any
        # pre-rejected response records supplied by merge.
        rejected_rows[claim_rejected_start:] = sorted(
            rejected_rows[claim_rejected_start:],
            key=lambda record: (record["claim_index"], record["source_index"]),
        )

        if rejected:
            status = "partial_import" if imported else "not_evaluated_all_claim_rows_rejected"
            note = "invalid claim rows quarantined without changing their scope or anchors"
        elif parse_status != "json":
            status = "unparsed"
            note = "model output was not parsed as JSON; no claim is eligible for verification"
        elif not claims:
            status = "no_claims"
            note = "model returned no semantic bridge for this frozen unit"
        elif unsupported:
            status = "contains_unsupported_claims"
            note = "one or more model citations/reference genes were not present in the frozen input bundle"
        else:
            status = "proposals_imported"
            note = "claims remain proposed until source-content/relation and transfer verification"
        status_rows.append({
            "unit_id": unit_id,
            "prompt_id": prompt_id,
            "parse_status": parse_status,
            "model_claim_count": len(claims) + len(pre_rejected),
            "imported_claim_row_count": imported,
            "unsupported_claim_row_count": unsupported,
            "rejected_claim_row_count": rejected,
            "status": status,
            "note": note,
        })

    for unit_id in sorted(not_requested):
        status_rows.append({
            "unit_id": unit_id,
            "prompt_id": "",
            "parse_status": "not_requested",
            "model_claim_count": 0,
            "imported_claim_row_count": 0,
            "unsupported_claim_row_count": 0,
            "rejected_claim_row_count": 0,
            "status": "not_evaluated_no_reviewable_material",
            "note": ("no family, domain, reference function or literature was supplied for this unit, "
                     "so no request was issued; this is not evaluated, not absence"),
        })

    if bundle.get("b_review_membership_policy") == "initial_bibliography_v1":
        from .paper_reports import stamp_paper_handoff
        rows = stamp_paper_handoff(rows, route="interval_first")
    claims_frame = pd.DataFrame(rows, columns=CLAIM_COLUMNS)
    if not claims_frame.empty:
        claims_frame = normalize_semantic_claims(claims_frame)
    status_frame = pd.DataFrame(status_rows, columns=STATUS_COLUMNS).sort_values("unit_id", kind="stable")

    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    claims_path = out / "semantic_claim_proposals.tsv"
    status_path = out / "semantic_claim_proposal_status.tsv"
    rejected_path = out / "semantic_claim_rejected_rows.json"
    rejected_path.write_text(json.dumps(rejected_rows, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    claims_frame.to_csv(claims_path, sep="\t", index=False, lineterminator="\n")
    status_frame.to_csv(status_path, sep="\t", index=False, lineterminator="\n")

    proposed_count = int((claims_frame["verification_status"] == "proposed").sum()) if len(claims_frame) else 0
    unsupported_count = int((claims_frame["verification_status"] == "unsupported").sum()) if len(claims_frame) else 0
    import_state = "not_evaluated_all_claim_rows_rejected" if rejected_rows and claims_frame.empty else (
        "partial_import" if rejected_rows else "imported")
    from .paper_evidence_table import write_paper_evidence_table
    paper_evidence = write_paper_evidence_table(claims_frame.to_dict("records"), out,
                                              route="interval_first", source_artifact=str(responses_path))
    manifest = {
        "paper_evidence": paper_evidence,
        "schema_version": SCHEMA_VERSION,
        "case_id": case_id or str(bundle.get("case_id", "")),
        "route": "interval_first",
        "stage": "semantic_proposal_import",
        "import_state": import_state,
        "semantic_unit_count": len(units),
        "response_count": len(by_prompt),
        "claim_row_count": len(claims_frame),
        "proposed_claim_row_count": proposed_count,
        "unsupported_claim_row_count": unsupported_count,
        "verified_claim_row_count": 0,
        "rejected_claim_row_count": len(rejected_rows),
        "rejected_rows_sha256": _sha256_file(str(rejected_path)),
        "provider": provider,
        "model": model,
        "llm_manifest": str(runtime_path(llm_manifest_arg)) if llm_manifest_arg else None,
        "llm_execution_evidence": _llm_execution_evidence(llm_manifest_arg),
        "artifact_hashes": artifact_hashes,
        "bundle_sha256": artifact_hashes.get("bundle_sha256") if artifact_hashes else None,
        "prompt_sha256": artifact_hashes.get("prompt_sha256") if artifact_hashes else None,
        "response_sha256": artifact_hashes.get("response_sha256") if artifact_hashes else None,
        "rank_contribution": "none",
        "policy": [
            "every Route B semantic claim retains semantic_unit_id from the exact frozen unit shown to the model",
            "the response roster must exactly match the frozen semantic prompt roster",
            "source ids and reference genes invented outside the frozen unit are unsupported",
            "a source id present in the retrieval bundle proves only input provenance, not that the source supports the model claim",
            "well-formed claims remain proposed until deterministic/source-content relation and transfer verification",
            "no proposal is narrowing-eligible and no numerical gene score is written at this stage",
            "llm_execution_evidence reports whether the response roster can evidence per-prompt evaluation; a uniform roster is not countable as Route B executed",
        ],
        "outputs": {
            **paper_evidence["outputs"],
            "proposals": str(claims_path),
            "status": str(status_path),
            "rejected_rows": str(rejected_path),
            "manifest": str(out / "semantic_claim_proposal_manifest.json"),
        },
    }
    manifest_path = out / "semantic_claim_proposal_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "paper_evidence": paper_evidence, **paper_evidence["outputs"],
        "proposals": str(claims_path),
        "import_state": import_state,
        "status": str(status_path),
        "manifest": str(manifest_path),
        "claim_row_count": len(claims_frame),
        "proposed_claim_row_count": proposed_count,
        "unsupported_claim_row_count": unsupported_count,
        "verified_claim_row_count": 0,
        "rejected_claim_row_count": len(rejected_rows),
        "rejected_rows": str(rejected_path),
    }
