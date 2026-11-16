"""Offline building blocks for practical packets and canonical verification/Wet.

Change class: CORE_INTEGRATION (issue #744). This module wires existing
canonical functions for ``prepare_practical_run.py`` packet outputs;
it adds no Route, ranker, score, threshold, or evidence meaning, and it never
calls a model or the network. Every external lookup stays offline
(``allow_external_api=False``); missing data stays explicit and is never
treated as negative evidence.

Packet layout consumed here is produced by
``prepare_practical_run.packet`` (see ``materialize_expanded_input``):
``phenotype.json`` (+ ``case.json``) for A/C lanes, and ``candidates.tsv``,
``candidate_proteins.fa``, ``reference_corpus.tsv``, ``case.json`` for B lanes.
LLM/atom responses remain external inputs: this adapter imports and verifies
frozen response artifacts, exactly like the existing native adapter path.
It does not execute A projection, author reviews, or orchestrate a packet-to-Wet
run. Those steps remain caller responsibilities; see the harness README.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import tempfile


_A_ENVELOPE_KEYS = {
    "hypotheses", "research_log", "direct_answer_exposure", "limitations",
}


def _require_file(path: Path, label: str) -> Path:
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"missing {label}: {path}")
    return path


def _unique_gene_ids(candidates_path: Path) -> list[str]:
    with candidates_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "gene_id" not in reader.fieldnames:
            raise ValueError("candidate inventory requires gene_id")
        ids = [(row["gene_id"] or "").strip() for row in reader]
    if not ids or "" in ids or len(ids) != len(set(ids)):
        raise ValueError("candidate inventory requires unique non-empty gene IDs")
    return ids


def read_packet_b(packet_dir: str | Path) -> dict[str, object]:
    """Validate a materialized B packet directory and return its file map."""
    root = Path(packet_dir)
    files = {
        "candidates": _require_file(root / "candidates.tsv", "B candidates"),
        "case": _require_file(root / "case.json", "B case metadata"),
        "phenotype": _require_file(root / "phenotype.json", "B phenotype"),
    }
    gene_ids = _unique_gene_ids(files["candidates"])
    try:
        case_metadata = json.loads(files["case"].read_text(encoding="utf-8"))
        interval = case_metadata["interval"]
        case_id = case_metadata["case_id"]
    except (json.JSONDecodeError, KeyError, TypeError) as error:
        raise ValueError(f"B case metadata requires case_id and interval: {error}")
    if not isinstance(interval, str) or not interval.strip():
        raise ValueError("B case metadata requires a non-empty interval")
    if not isinstance(case_id, str) or not case_id.strip():
        raise ValueError("B case metadata requires a non-empty case_id")
    return {"dir": str(root), "files": {k: str(v) for k, v in files.items()},
            "gene_ids": gene_ids, "interval": interval, "case_id": case_id}


def build_atlas(packet: dict[str, object], out_dir: str | Path, *,
                case_id: str, family_index: str | None = None,
                domains: str | None = None, orthology: str | None = None,
                annotations: str | None = None) -> dict[str, object]:
    """Build the phenotype-free interval atlas from a validated B packet.

    ``reference_corpus.tsv`` is deliberately not consumed here: it is B
    semantic reference material for the bundle stage, not an atlas input, and
    this adapter never invents a family mapping from it.
    """
    from candidate_qtl.interval_atlas import interval_atlas

    if not case_id.strip():
        raise ValueError("build_atlas requires a non-empty case_id")
    result = interval_atlas(str(packet["files"]["candidates"]), str(out_dir),
                            family_index_arg=family_index, domains_arg=domains,
                            orthology_arg=orthology, annotations_arg=annotations,
                            case_id=case_id)
    result["atlas"] = result["interval_gene_atlas"]
    return result


def import_a_payload(payload_path: str | Path, out_dir: str | Path, *,
                     case_id: str, model: str,
                     payload_source: str = "llm") -> dict[str, object]:
    """Import frozen A hypotheses (LLM/agent/human JSON) to canonical TSVs."""
    from candidate_qtl.mechanism_payload_adapter import forbidden_columns
    from candidate_qtl.phenotype_investigation import mechanism_import

    if not case_id.strip():
        raise ValueError("import_a_payload requires a non-empty case_id")
    if not model.strip():
        raise ValueError("import_a_payload requires the real proposer model ref")
    source = Path(payload_path)
    try:
        raw = source.read_bytes()
        envelope = json.loads(raw.decode("utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"A payload must be a readable JSON envelope: {error}") from error
    if not isinstance(envelope, dict):
        raise ValueError("A payload envelope must be a JSON object")
    unknown = [key for key in envelope if key not in _A_ENVELOPE_KEYS]
    forbidden = forbidden_columns([str(key) for key in unknown])
    if forbidden:
        raise ValueError(f"A payload envelope contains forbidden key(s): {forbidden}")
    if unknown:
        raise ValueError(f"A payload envelope contains unknown top-level key(s): {unknown}")
    hypotheses = envelope.get("hypotheses")
    if not isinstance(hypotheses, list) or any(not isinstance(row, dict) for row in hypotheses):
        raise ValueError("A payload envelope hypotheses must be a list of objects")
    if "research_log" in envelope and (
        not isinstance(envelope["research_log"], list)
        or any(not isinstance(row, dict) for row in envelope["research_log"])
    ):
        raise ValueError("A payload envelope research_log must be a list of objects")
    if "direct_answer_exposure" in envelope and not isinstance(
        envelope["direct_answer_exposure"], str
    ):
        raise ValueError("A payload envelope direct_answer_exposure must be a string")
    if "limitations" in envelope and not (
        isinstance(envelope["limitations"], (str, dict))
        or (isinstance(envelope["limitations"], list)
            and all(isinstance(item, str) for item in envelope["limitations"]))
    ):
        raise ValueError("A payload envelope limitations must be a string, list of strings or JSON object")

    def reject_nested_leaks(value: object, location: str) -> None:
        if isinstance(value, dict):
            leaked = forbidden_columns([
                str(key) for key in value if key != "direct_answer_exposure"
            ])
            if leaked:
                raise ValueError(
                    f"A payload envelope contains forbidden nested key(s) at {location}: {leaked}"
                )
            for key, nested in value.items():
                if key == "direct_answer_exposure" and not isinstance(nested, str):
                    raise ValueError(
                        f"A payload envelope {location}.direct_answer_exposure must be a string"
                    )
                reject_nested_leaks(nested, f"{location}.{key}")
        elif isinstance(value, list):
            for index, nested in enumerate(value):
                reject_nested_leaks(nested, f"{location}[{index}]")

    for key in ("research_log", "direct_answer_exposure", "limitations"):
        if key in envelope:
            reject_nested_leaks(envelope[key], key)

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    original_copy = out / "raw_a_payload_original.json"
    original_copy.write_bytes(raw)
    normalized_bytes = (json.dumps({"hypotheses": hypotheses}, ensure_ascii=False,
                                   indent=2) + "\n").encode("utf-8")
    with tempfile.TemporaryDirectory(prefix="a-payload-") as temporary:
        normalized_path = Path(temporary) / "hypotheses.json"
        normalized_path.write_bytes(normalized_bytes)
        result = mechanism_import(
            str(normalized_path), str(out), payload_source=payload_source,
            model=model, case_id=case_id)
    provenance_path = out / "a_payload_envelope_provenance.json"
    provenance = {
        "schema": "type12-a-payload-envelope-v1",
        "original_payload": str(original_copy),
        "original_payload_source": str(source),
        "original_payload_sha256": hashlib.sha256(raw).hexdigest(),
        "normalized_hypotheses_sha256": hashlib.sha256(normalized_bytes).hexdigest(),
        "retained_metadata": {
            key: envelope[key] for key in _A_ENVELOPE_KEYS - {"hypotheses"}
            if key in envelope
        },
        "canonical_hypotheses": str(result["hypotheses"]),
    }
    provenance_path.write_text(json.dumps(provenance, ensure_ascii=False,
                                          indent=2, sort_keys=True) + "\n",
                               encoding="utf-8")
    result["raw_a_payload"] = str(original_copy)
    result["a_payload_envelope_provenance"] = str(provenance_path)
    result["raw_a_payload_sha256"] = provenance["original_payload_sha256"]
    return result


def verify_and_freeze_a(hypotheses_tsv: str | Path, out_dir: str | Path, *,
                        focal_species: str, record_cache: str | None = None,
                        source_records: str | None = None,
                        bound_claims: str | None = None,
                        agent_claim_review: str | None = None,
                        source_review_records: str | None = None,
                        case_id: str = "") -> dict[str, object]:
    """Verify A hypotheses fully offline and freeze bait panels.

    Without bound claims (or agent review), verified sequences land in the
    review-only projection panel, never in the supported panel: review status
    is preserved, support is never promoted.
    """
    from candidate_qtl.llm_pathway import bait_panel_freeze, hypothesis_verify

    verify_dir = Path(out_dir) / "verify"
    verify = hypothesis_verify(str(hypotheses_tsv), str(verify_dir), focal_species,
                               record_cache_arg=record_cache,
                               allow_external_api=False,
                               source_records_arg=source_records,
                               bound_claims_arg=bound_claims,
                               agent_claim_review_arg=agent_claim_review,
                               source_review_records_arg=source_review_records)
    panel_dir = Path(out_dir) / "panel"
    panel = bait_panel_freeze(
        str(hypotheses_tsv), str(verify_dir / "proposal_validation.tsv"),
        str(verify_dir / "verified_sequences.json"), str(panel_dir),
        case_id=case_id,
        reviewable_sequences_arg=str(verify_dir / "reviewable_sequences.json"))
    return {"verify": verify, "panel": panel,
            "validation": str(verify_dir / "proposal_validation.tsv"),
            "verified_sequences": str(verify_dir / "verified_sequences.json"),
            "reviewable_sequences": str(verify_dir / "reviewable_sequences.json"),
            "verified_panel": str(panel_dir / "verified_bait_panel.tsv"),
            "projection_panel": str(panel_dir / "projection_bait_panel.tsv")}


def build_bundle(atlas_path: str | Path, out_dir: str | Path, *,
                 phenotype: str, literature: str | None = None,
                 case_id: str = "") -> dict[str, object]:
    """Build the phenotype-aware semantic bundle from a frozen atlas."""
    from candidate_qtl.route_b_semantic_bundle import build_route_b_semantic_bundle

    return build_route_b_semantic_bundle(str(atlas_path), str(out_dir),
                                         phenotype=phenotype,
                                         literature_arg=literature,
                                         case_id=case_id)


def import_verify_freeze_b(bundle_path: str | Path, prompts_path: str | Path,
                           responses_path: str | Path, out_dir: str | Path, *,
                           case_id: str = "",
                           source_cache: str | None = None,
                           bound_claims: str | None = None,
                           agent_claim_review: str | None = None,
                           source_review_records: str | None = None,
                           record_cache: str | None = None) -> dict[str, object]:
    """Import frozen B responses and run source/transfer verification offline.

    Empty/unread units stay ``not_evaluated``; they are never auto-filled and
    never converted to negatives.
    """
    from candidate_qtl.route_b_semantic_candidates import (
        freeze_route_b_semantic_candidates,
    )
    from candidate_qtl.route_b_semantic_import import (
        import_route_b_semantic_proposals,
    )
    from candidate_qtl.route_b_semantic_source_verify import (
        verify_route_b_semantic_sources,
    )
    from candidate_qtl.route_b_semantic_transfer_verify import (
        verify_route_b_semantic_transfer,
    )

    out = Path(out_dir)
    imported = import_route_b_semantic_proposals(
        str(bundle_path), str(prompts_path), str(responses_path),
        str(out / "import"), case_id=case_id)
    sourced = verify_route_b_semantic_sources(
        str(imported["proposals"]), str(out / "source"),
        bundle_arg=str(bundle_path), cache_arg=source_cache, allow_external_api=False,
        case_id=case_id)
    transfer = verify_route_b_semantic_transfer(
        str(imported["proposals"]), str(sourced["source_verification"]),
        str(bundle_path), str(out / "transfer"), case_id=case_id,
        bound_claims_arg=bound_claims,
        agent_claim_review_arg=agent_claim_review,
        source_review_records_arg=source_review_records,
        record_cache_arg=record_cache)
    frozen = freeze_route_b_semantic_candidates(
        str(transfer["claims"]), str(out / "candidates"), case_id=case_id)
    return {"import": imported, "source": sourced, "transfer": transfer,
            "freeze": frozen,
            "proposals": str(imported["proposals"]),
            "source_verification": str(sourced["source_verification"]),
            "transfer_claims": str(transfer["claims"]),
            "transfer_verification": str(transfer["transfer_verification"]),
            "transfer_manifest": str(transfer["manifest"]),
            "candidates": str(frozen["candidates"]),
            "candidates_manifest": str(frozen["manifest"]),
            "review_candidates": str(frozen["review_candidates"]),
            "review_manifest": str(frozen["review_manifest"])}


def handoff_to_wet(out_dir: str | Path, *, case_id: str, wet_budget: int = 20,
                   route_b_claims: str | None = None,
                   route_b_candidates: str | None = None,
                   route_b_claims_manifest: str | None = None,
                   route_b_candidates_manifest: str | None = None,
                   route_a_bridge: str | None = None,
                   route_a_manifest: str | None = None,
                   interval_atlas: str | None = None,
                   review_lanes: list[str] | None = None,
                   review_manifests: list[str] | None = None) -> dict[str, object]:
    """Run the canonical Wet handoff on frozen canonical inputs.

    Supported members are always kept; review candidates enter ahead of
    generic fill without support promotion; fill carries no support flags and
    no score. ``wet_budget`` is an operational target, not a truncation cap.
    """
    from candidate_qtl.wet_candidate_handoff import wet_candidate_handoff

    if not case_id.strip():
        raise ValueError("handoff_to_wet requires a non-empty case_id")
    return wet_candidate_handoff(
        str(out_dir), case_id=case_id, wet_budget=wet_budget,
        route_b_claims_arg=route_b_claims,
        route_b_candidates_arg=route_b_candidates,
        route_b_claims_manifest_arg=route_b_claims_manifest,
        route_b_candidates_manifest_arg=route_b_candidates_manifest,
        route_a_bridge_arg=route_a_bridge,
        route_a_manifest_arg=route_a_manifest,
        interval_atlas_arg=interval_atlas,
        review_candidates_arg=review_lanes,
        review_candidates_manifest_arg=review_manifests)
