"""Family candidate bridge: normalize Route A / B / C family evidence into one table.

T1 scope is the deterministic normalization core only. The bridge merges three
independently produced route artifacts -- Route A (phenotype-first hypotheses
with verified bait status and interval family match), Route B (orthology table
and interval gene atlas), and Route C (frozen retrieval interval match) -- into
one reporting table with a strict mapping tier and an explicit
``not_evaluated`` vocabulary. It never scores: ``rank_contribution`` is always
``none``, there are no score columns, and the bridge is a parallel reporting
artifact -- ``family_to_copy.py`` / ``narrowing.py`` wiring is unchanged.

Deliberate boundaries (mirroring the honesty rules of the rest of the
workflow):

* A falsified family (``falsified_no_member_after_complete_sequence_search`` /
  ``falsified_family_not_detected_genome_wide``) is not family evidence.  A
  separately verified Route A bait may nevertheless project its complete
  threshold-passing BLAST target set into the interval; multiple targets stay
  explicitly paralog-ambiguous and the family-index disagreement is recorded.
* A hypothesis without a verified bait is never a match: one
  ``tier6_not_evaluated`` row per proposed member gene with
  ``not_evaluated_no_verified_bait``, and ``not_run`` provider -- never a miss,
  never a blastp claim.
* ``name_only`` (naming-bucket) rows are strictly the weakest evaluated tier
  (``tier5_name_only``) and are never promoted by any other column.
* Coverage / identity below threshold never promotes: source-backed rows are
  preserved at ``tier2_source_backed_ambiguous`` (evidence kept, ambiguity
  ``coverage_insufficient``, review required); sequence-only rows become
  ``tier6_not_evaluated`` with ``not_evaluated_coverage_below_threshold``.
* Route B ``best_hit_homology`` is a local best hit, not a reciprocal one
  (route_b_prepare writes exactly this label and documents that it claims no
  one-to-one orthology), so it maps to ``blastp``.
* Route C rows never claim source-backed orthology: an orthogroup-key match
  from the frozen retrieval corpus is a sequence-match claim, so Route C
  ``orthogroup`` rows are capped at ``tier3_sequence_match_unique`` /
  ``tier4_sequence_match_ambiguous``.
* A row's ``source_id`` is never promoted to a dependency root on its own:
  roots come only from the supplied dependency-clusters table, from the
  verified-bait ``function_source_ids`` (exactly one value), or from the
  Route C ``database_version``.
* Nothing absolute leaks into the artifacts: fingerprints are content hashes
  only, and the manifest ``outputs`` are relative filenames.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re

import pandas as pd

from .core import file_fingerprint, runtime_path, sha256
from .family_labels import families_match_expanded
# The schema enums are re-exported here (public bridge vocabulary).
from .schemas import (  # noqa: F401
    BRIDGE_AMBIGUITY_STATUSES,
    BRIDGE_COLUMNS,
    BRIDGE_MAPPING_TIERS,
    BRIDGE_MATCH_TYPES,
    BRIDGE_NOT_EVALUATED_REASONS,
    BRIDGE_PROVIDER_STATUSES,
    BRIDGE_REVIEW_STATUSES,
    FAMILY_CANDIDATE_BRIDGE_SCHEMA,
)

BRIDGE_SCHEMA_VERSION = "1.0"
ROUTES = frozenset({"a", "b", "c"})
MIN_IDENTITY_FRACTION = 0.3
MIN_COVERAGE_FRACTION = 0.5

# Orthology-class relationships that claim reciprocal / table-backed
# orthology. ``best_hit_homology`` (the route_b_prepare label for a local best
# hit) is deliberately absent: it explicitly claims no one-to-one orthology,
# so it maps to blastp.
_RBH_RELATIONSHIPS = frozenset({
    "one_to_one", "one_to_many", "many_to_many", "close_paralog",
    "family_level_only", "unmapped_ambiguous",
})
_AMBIGUOUS_RELATIONSHIPS = frozenset({"one_to_many", "many_to_many", "close_paralog", "unmapped_ambiguous"})
_NOT_EVALUATED_PROVIDER_STATUSES = frozenset({"failed", "unavailable", "not_run"})
_COMPLETE_PROVIDER_STATUSES = frozenset({"complete", "imported_complete"})
_FALSIFIED_FAMILY_STATUSES = frozenset({
    "falsified_no_member_after_complete_sequence_search",
    "falsified_family_not_detected_genome_wide",
})
# Route A states where the family question could not be answered, so the
# verified bait's own sequence evidence is the only route into the interval
#. The bait was searched against the interval proteins regardless,
# and that evidence is independent of any family-name space: an
# Arabidopsis-named bait matches a Medicago-named interval protein without
# either side sharing a symbol. The projection is positive-only -- absence
# still needs the genome-wide index -- and the paralog guard stays explicitly
# unevaluated.
#
# Both members describe the same epistemic state. An index that does not exist
# and an index that does not carry this family differ in why the question went
# unanswered, not in whether it was answered. Admitting only the first silently
# discarded the bait hits in every run where an index existed but was
# incomplete, which is the ordinary case rather than the exception.
_SEQUENCE_PROJECTION_FAMILY_STATUSES = frozenset({
    "not_evaluated_missing_reference_genome_family_index",
    "not_evaluated_family_unresolved_in_index",
})
_UNCLASSIFIED_FAMILY_STATUSES = frozenset({"surviving_sequence_only_annotation_gap", "hold_absent_annotation_only"})
_MATCHED_FAMILY_STATUSES = frozenset({
    "concrete_interval_gene_matched", "family_present_copy_unresolved",
    "surviving_family_present_in_interval",
})
_TIER6 = "tier6_not_evaluated"
_REVIEW_NOT_EVALUATED = "not_evaluated"
_DEPENDENCY_UNKNOWN = "dependency_unknown"
ROUTE_A_REVIEW_CANDIDATE_COLUMNS = [
    "case_id", "candidate_gene_id", "route", "support_status", "membership_basis", "review_status", "eligibility",
    "hypothesis_id", "reference_accession", "source_id", "identity", "query_coverage",
    "target_coverage", "evalue", "ambiguity_status", "provider_status", "rank_contribution",
    "reason", "input_fingerprint", "relevance_evidence_json",
]


def _read(path_arg: str | None, required: set[str], label: str) -> pd.DataFrame:
    if not path_arg:
        return pd.DataFrame(columns=sorted(required))
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{label} table is missing columns: {', '.join(missing)}")
    return frame


def _tokens(value: object, separator: str = ";") -> list[str]:
    return [part.strip() for part in str(value or "").replace(",", ";").split(separator) if part.strip()]


_NEGATED_PI_FAMILY = re.compile(
    r"\b(?:non|not|without|excluding|except)\b[^;\n]{0,24}"
    r"\b(?:pi(?:[- ]class)?|pistillata)\b",
    re.I,
)


def _positive_pi_family_label(value: object) -> bool:
    """Reject a negated PI/PISTILLATA claim before symbol token matching."""
    return not _NEGATED_PI_FAMILY.search(str(value or ""))


def _row_input_fingerprint(source_row: dict[str, object], source_columns: list[str]) -> str:
    text = "\t".join(str(source_row.get(column, "")) for column in sorted(source_columns))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _row_output_fingerprint(bridge_row: dict[str, object]) -> str:
    columns = sorted(column for column in BRIDGE_COLUMNS if column != "output_fingerprint")
    text = "\t".join(str(bridge_row.get(column, "")) for column in columns)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _finalize_row(row: dict[str, str]) -> dict[str, str]:
    """Stamp the output fingerprint over the fully-resolved row.

    Called once, after dependency clusters are assigned, so the fingerprint
    covers every final column value.
    """
    row["output_fingerprint"] = _row_output_fingerprint(row)
    return row


def _cluster_id(root_source_id: str) -> str:
    digest = hashlib.sha256(root_source_id.encode("utf-8")).hexdigest()[:20]
    return f"dep_{digest}"


_ORTHOGROUP_TOKEN_PREFIX = "orthogroup:"


def orthogroup_dependency_token(value: object) -> str:
    """The source token for a row whose homology claim is an orthogroup/HOG id.

    Route B rows without an explicit orthology ``source_id`` carry this token
    as their source identity, and the route-b-prepared dependency-clusters
    table (``evidence_dependency_clusters.tsv``) lists the same token in its
    ``source_ids``, so rows derived from the same orthogroup/HOG resolve to the
    same dependency cluster. Resolution still happens only through the supplied
    clusters table: the token alone is never promoted to a root.
    """
    return f"{_ORTHOGROUP_TOKEN_PREFIX}{str(value or '').strip()}"


def dependency_cluster_id(token: str) -> str:
    """Cluster id of one dependency token (same derivation as _cluster_id)."""
    return _cluster_id(token)


def _coverage_pass(identity: object, query_coverage: object, target_coverage: object) -> tuple[bool, bool]:
    """Return (pass, applicable) for the coverage/identity thresholds.

    Empty values are not applicable (pass). ``identity`` values below 1.0 are
    fractions, values >= 1.0 are percents divided by 100. Unparseable values
    fail closed: a threshold that cannot be confirmed never promotes.
    """
    applicable = False
    for value, minimum, is_identity in (
        (identity, MIN_IDENTITY_FRACTION, True),
        (query_coverage, MIN_COVERAGE_FRACTION, False),
        (target_coverage, MIN_COVERAGE_FRACTION, False),
    ):
        text = str(value).strip()
        if not text:
            continue
        applicable = True
        try:
            number = float(text)
        except ValueError:
            return False, True
        if is_identity and number >= 1.0:
            number = number / 100.0
        if number < minimum:
            return False, True
    return True, applicable


def _route_a_blast_projection_eligible(blast: dict[str, object]) -> bool:
    """Return whether one frozen Route A BLAST row can project its top hit.

    This is deliberately the same defensive gate already used when a BLAST
    layer reaches ``family_to_copy``: an evaluable, verified, unique hit with
    a concrete target and the existing bridge identity/coverage floors.  The
    bridge does not infer a family member from a falsified family; it only
    accepts the actual top-hit target recorded by ``blast_evidence.tsv``.
    """
    if str(blast.get("provider_status", "")).strip() != "evaluable":
        return False
    if str(blast.get("review_status", "")).strip() != "verified":
        return False
    if str(blast.get("ambiguity_status", "")).strip() != "unique":
        return False
    if not str(blast.get("top_hit_target_id", "")).strip():
        return False
    values = (blast.get("top_hit_identity", ""), blast.get("top_hit_query_coverage", ""),
              blast.get("top_hit_target_coverage", ""))
    if any(not str(value).strip() for value in values):
        return False
    passes, applicable = _coverage_pass(*values)
    return bool(passes and applicable)


_LONG_FORM_ELIGIBLE_NAMES = (
    "eligible_interval_hits.tsv",
    "blast_eligible_hits.tsv",
    "eligible_blast_hits.tsv",
    "blast_panel_eligible_hits.tsv",
    "blast_evidence_long.tsv",
)


def _discover_long_form_blast_artifact(short_path: str | None) -> Path | None:
    """Find the immutable long-form hit artifact next to blast evidence.

    ``blast_evidence.tsv`` is intentionally a one-row-per-bait summary.  The
    search stage also emits a long-form, threshold-filtered table for copy
    review.  Keep discovery permissive across the v1/v2 artifact names, but
    never treat an arbitrary sibling TSV as evidence: only a manifest pointer
    or the names above are accepted.
    """
    if not short_path:
        return None
    source = runtime_path(short_path)
    if not source.exists():
        return None
    manifest = source.parent / "blast_panel_search_manifest.json"
    if manifest.exists():
        try:
            payload = json.loads(manifest.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            payload = {}
        def pointers(value: object):
            if isinstance(value, dict):
                for key, child in value.items():
                    key_text = str(key).casefold()
                    if isinstance(child, str) and child.casefold().endswith(".tsv") \
                            and ("eligible" in key_text or "long" in key_text):
                        yield child
                    yield from pointers(child)
            elif isinstance(value, list):
                for child in value:
                    yield from pointers(child)
        for pointer in pointers(payload):
            candidate = Path(pointer)
            if not candidate.is_absolute():
                candidate = manifest.parent / candidate
            if candidate.exists():
                expected = str(payload.get("eligible_interval_hits_sha256", "")).strip()
                if expected and sha256(candidate) != expected:
                    raise ValueError(
                        "eligible interval hits fingerprint does not match blast-panel-search manifest"
                    )
                return candidate
    if source.name.casefold() in {name.casefold() for name in _LONG_FORM_ELIGIBLE_NAMES}:
        return source
    for name in _LONG_FORM_ELIGIBLE_NAMES:
        candidate = source.parent / name
        if candidate.exists():
            return candidate
    return None


def _long_form_blast_hits(
    short_path: str | None,
    short_rows: dict[str, dict[str, object]],
) -> dict[str, list[dict[str, str]]]:
    """Load all eligible interval targets, keyed by verified bait accession.

    The helper accepts the stable names used by the search stage and a small
    alias set so the bridge remains compatible with an artifact generated by a
    newer search worker.  Eligibility is never inferred from rank: rows must
    be explicitly eligible or carry complete threshold metrics that pass the
    existing bridge floors.
    """
    path = _discover_long_form_blast_artifact(short_path)
    if path is None:
        return {}
    try:
        frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False).fillna("")
    except (OSError, ValueError, pd.errors.ParserError):
        return {}
    def pick(record: dict[str, object], *names: str) -> str:
        for name in names:
            value = str(record.get(name, "")).strip()
            if value:
                return value
        return ""
    out: dict[str, list[dict[str, str]]] = {}
    for record in frame.to_dict("records"):
        panel = pick(record, "panel_id", "query_id", "bait_id", "homolog_query_id")
        target = pick(record, "target_id", "target_gene_id", "hit_target_id", "sseqid", "gene_id")
        if not panel or not target or panel not in short_rows:
            continue
        if pick(record, "review_status").casefold() not in {"verified", "review_required"}:
            continue
        if pick(record, "provider_status").casefold() != "evaluable":
            continue
        if pick(record, "rank_contribution").casefold() != "none":
            continue
        eligibility = pick(record, "eligibility", "hit_eligibility", "threshold_status", "eligible")
        if eligibility and eligibility.casefold() not in {"eligible", "threshold_pass", "threshold_passing", "passed", "true", "1"}:
            continue
        identity = pick(record, "identity", "pident", "hit_identity", "top_hit_identity")
        query_coverage = pick(record, "query_coverage", "qcov", "hit_query_coverage", "top_hit_query_coverage")
        target_coverage = pick(record, "target_coverage", "tcov", "hit_target_coverage", "top_hit_target_coverage")
        evalue = pick(record, "evalue", "hit_evalue", "top_hit_evalue")
        passes, applicable = _coverage_pass(identity, query_coverage, target_coverage)
        # The long-form producer is authoritative about threshold eligibility,
        # but incomplete/failed rows are never allowed through by discovery.
        if not applicable or not passes:
            continue
        out.setdefault(panel, []).append({
            "panel_id": panel,
            "target_id": target,
            "identity": identity,
            "query_coverage": query_coverage,
            "target_coverage": target_coverage,
            "evalue": evalue,
            "bitscore": pick(record, "bitscore", "top_hit_bitscore"),
            "source_row": str(record.get("row_id", "")),
        })
    for panel, rows in out.items():
        dedup: dict[str, dict[str, str]] = {}
        for row in rows:
            dedup.setdefault(row["target_id"], row)
        out[panel] = sorted(dedup.values(), key=lambda row: (row["target_id"], row["bitscore"]))
    return out


def route_a_homolog_projection(row: dict, members: list[str]) -> str:
    """Describe an already frozen, coverage-qualified A projection for Tier use."""
    from math import isfinite
    values = [row.get(key, "") for key in ("identity", "query_coverage", "target_coverage")]
    try:
        usable = all(str(value).strip() and isfinite(float(value)) for value in values)
    except (ValueError, TypeError):
        usable = False
    if (not usable or not _coverage_pass(*values)[0]
            or row.get("ambiguity_status") == "partial_alignment_unresolved"):
        return ""
    return json.dumps({"basis": "qualified_sequence_hit", "route": "a",
                       "candidate_gene_id": row["candidate_gene_id"],
                       "reference_gene_id": row["reference_accession"],
                       "member_genes": sorted(set(members)), "copy_count": len(set(members)),
                       "ambiguity": row.get("ambiguity_status", ""),
                       "evidence": {key: row.get(key, "") for key in
                                    ("identity", "query_coverage", "target_coverage", "evalue", "input_fingerprint")}},
                      ensure_ascii=False, sort_keys=True)


def _route_a_review_projection_rows(
    projection_panel: pd.DataFrame,
    blast_evidence: pd.DataFrame,
    long_form_by_accession: dict[str, list[dict[str, str]]],
    interval_gene_ids: set[str] | None,
    case_id: str,
    partial_by_accession: dict[str, list[dict[str, str]]] | None = None,
    hypotheses: list[dict[str, object]] | None = None,
) -> list[dict[str, str]]:
    """Emit review-only interval candidates from identity-verified BLAST baits.

    This deliberately does not use family labels. Frozen identity/species/
    sequence is mandatory. Incomplete claims or significant partial alignments
    remain review-only within the supplied interval universe.
    """
    if projection_panel.empty or blast_evidence.empty or interval_gene_ids is None:
        return []
    blast_by_id = {str(row.get("panel_id", "")).strip(): row for row in blast_evidence.to_dict("records")}
    hypotheses_by_id: dict[str, list[dict[str, object]]] = {}
    for hypothesis in hypotheses or []:
        hypotheses_by_id.setdefault(str(hypothesis.get("hypothesis_id", "")).strip(), []).append(hypothesis)

    def evidence_for(panel: dict[str, object], target: str) -> str:
        hypothesis_id = str(panel.get("hypothesis_id", "")).strip()
        accession = str(panel.get("bait_id", "")).strip()
        evidence: list[dict[str, object]] = []
        for hypothesis in hypotheses_by_id.get(hypothesis_id, []):
            references = {
                str(hypothesis.get(key, "")).strip()
                for key in ("proposed_stable_id", "reference_accession", "reference_gene_id")
                if str(hypothesis.get(key, "")).strip()
            }
            if references and accession not in references:
                continue
            active_hypothesis = dict(hypothesis)
            if hypothesis.get("candidate_admission_policy") == "initial_bibliography_v1":
                # The sealed hypothesis input keeps this provider audit. Do not
                # promote its bibliography into the active evidence view.
                active_hypothesis.pop("reference_record_source_ids", None)
                provenance = active_hypothesis.pop("retrieval_provenance", "")
                if provenance:
                    try:
                        provenance = json.loads(str(provenance))
                    except (ValueError, TypeError):
                        provenance = None
                    if isinstance(provenance, dict):
                        provenance.pop("provider_linked_source_ids", None)
                        active_hypothesis["retrieval_provenance"] = json.dumps(
                            provenance, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
                    else:
                        active_hypothesis["retrieval_provenance_state"] = "audit_only_unparsed"
            # Lossless vocabulary aliases consumed by relevance_tiers; source
            # fields remain in the original record as well.
            evidence.append({
                **active_hypothesis,
                "candidate_gene_id": target,
                "family_id": str(hypothesis.get("family_id", "")).strip()
                or str(hypothesis.get("gene_family", "")).strip(),
                "reference_gene": str(hypothesis.get("reference_gene", "")).strip()
                or str(hypothesis.get("reference_gene_name", "")).strip()
                or str(hypothesis.get("proposed_stable_id", "")).strip(),
                "phenotype_relation": str(hypothesis.get("phenotype_relation", "")).strip()
                or str(hypothesis.get("expected_role", "")).strip(),
                "evidence_route": "a",
                "evidence_origin": "route_a_hypothesis",
            })
        return json.dumps(
            sorted(evidence, key=lambda row: json.dumps(row, ensure_ascii=False, sort_keys=True)),
            ensure_ascii=False, sort_keys=True, separators=(",", ":"),
        )
    out: list[dict[str, str]] = []
    for panel in projection_panel.to_dict("records"):
        review_only = (str(panel.get("projection_eligibility", "")).casefold() == "review_only"
                       and str(panel.get("review_status", "")).casefold() == "review_required")
        verified = (str(panel.get("projection_eligibility", "")).casefold() == "verified_support"
                    and str(panel.get("review_status", "")).casefold() == "verified")
        if not (review_only or verified):
            continue
        accession = str(panel.get("bait_id", "")).strip()
        blast = blast_by_id.get(accession)
        if blast is None or str(blast.get("provider_status", "")).strip() != "evaluable":
            continue
        from .observation_tiers import checklist_tier
        from .paper_reports import PAPER_FIELD, normalize_paper_reports
        focal_report = any(checklist_tier(report.get("paper_checklist")) == "T1"
                           for hypothesis in hypotheses_by_id.get(str(panel.get("hypothesis_id", "")), [])
                           for report in json.loads(normalize_paper_reports(hypothesis.get(PAPER_FIELD)) or "[]"))
        hits = list(long_form_by_accession.get(accession, [])) if review_only or focal_report else []
        if not hits and (review_only or focal_report):
            target = str(blast.get("top_hit_target_id", "")).strip()
            passes, applicable = _coverage_pass(
                blast.get("top_hit_identity", ""), blast.get("top_hit_query_coverage", ""),
                blast.get("top_hit_target_coverage", ""),
            )
            if target and applicable and passes:
                hits = [{
                    "target_id": target,
                    "identity": str(blast.get("top_hit_identity", "")).strip(),
                    "query_coverage": str(blast.get("top_hit_query_coverage", "")).strip(),
                    "target_coverage": str(blast.get("top_hit_target_coverage", "")).strip(),
                    "evalue": str(blast.get("top_hit_evalue", "")).strip(),
                }]
        hits.extend((partial_by_accession or {}).get(accession, []))
        for hit in hits:
            target = str(hit.get("target_id", "")).strip()
            if target not in interval_gene_ids:
                continue
            identity = str(hit.get("identity", "")).strip()
            query_coverage = str(hit.get("query_coverage", "")).strip()
            target_coverage = str(hit.get("target_coverage", "")).strip()
            passes, applicable = _coverage_pass(identity, query_coverage, target_coverage)
            partial = hit.get("projection_basis") == "partial_alignment_review_only"
            if not partial and (not applicable or not passes):
                continue
            fingerprint = hashlib.sha256((
                _row_input_fingerprint(panel, sorted(panel)) + "\n" +
                _row_input_fingerprint(hit, sorted(hit))
            ).encode("utf-8")).hexdigest()
            out.append({
                "case_id": case_id, "candidate_gene_id": target, "route": "a",
                "support_status": "review_required", "membership_basis": "route_a_review_projection", "review_status": "review_required",
                "eligibility": "review_only", "hypothesis_id": str(panel.get("hypothesis_id", "")),
                "reference_accession": accession, "source_id": str(panel.get("function_source_ids", "")),
                "identity": identity, "query_coverage": query_coverage, "target_coverage": target_coverage,
                "evalue": str(hit.get("evalue", "")).strip(),
                "ambiguity_status": "partial_alignment_unresolved" if partial else str(blast.get("ambiguity_status", "")).strip(),
                "provider_status": "complete", "rank_contribution": "none",
                "reason": "partial_alignment_requires_domain_or_copy_review" if partial else "source_claim_binding_not_verified", "input_fingerprint": fingerprint,
                "relevance_evidence_json": evidence_for(panel, target),
            })
    for row in out:
        # Keep the projection/partial-alignment record as well as its original
        # hypothesis. Nested evidence must not hide mapping uncertainty.
        original = {key: value for key, value in row.items() if key != "relevance_evidence_json"}
        members = [item["candidate_gene_id"] for item in out
                   if (item["hypothesis_id"], item["reference_accession"]) == (row["hypothesis_id"], row["reference_accession"])]
        projection = route_a_homolog_projection(row, members)
        evidence = [{**item, "homolog_projection_json": projection}
                    for item in json.loads(row["relevance_evidence_json"])]
        row["relevance_evidence_json"] = json.dumps(
            [*evidence, original],
            ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return sorted(out, key=lambda row: (row["candidate_gene_id"], row["reference_accession"], row["hypothesis_id"]))


def _derive_ambiguity(
    *,
    provider_status: str,
    reason: str,
    identity: str,
    query_coverage: str,
    target_coverage: str,
    matched_gene_count: int,
    paralog_count: str,
    family_status: str,
    relationship: str,
    tandem_ambiguous: bool,
) -> str:
    """Independent of the mapping tier (design order (1)-(4))."""
    if provider_status in _NOT_EVALUATED_PROVIDER_STATUSES or reason:
        return "not_evaluated"
    coverage_ok, _ = _coverage_pass(identity, query_coverage, target_coverage)
    if not coverage_ok:
        return "coverage_insufficient"
    if (matched_gene_count >= 2 or _paralog_number(paralog_count) >= 2
            or family_status == "family_present_copy_unresolved"
            or relationship in _AMBIGUOUS_RELATIONSHIPS):
        return "tandem_ambiguous" if tandem_ambiguous else "paralog_ambiguous"
    return "resolved_single_copy"


def _paralog_number(value: object) -> int:
    try:
        return int(str(value).strip())
    except ValueError:
        return 0


def _derive_review_status(mapping_tier: str, ambiguity_status: str, match_type: str) -> str:
    if mapping_tier == _TIER6:
        return _REVIEW_NOT_EVALUATED
    if match_type == "name_only":
        return "review_required"
    if ambiguity_status in {"paralog_ambiguous", "tandem_ambiguous", "coverage_insufficient"}:
        return "review_required"
    return "not_required"


def _derive_tier(
    *,
    route: str,
    match_type: str,
    ambiguity_status: str,
    provider_status: str,
    reason: str,
    identity: str,
    query_coverage: str,
    target_coverage: str,
) -> tuple[str, str]:
    """Exact decision order; returns (mapping_tier, reason)."""
    if provider_status in _NOT_EVALUATED_PROVIDER_STATUSES or reason:
        return _TIER6, reason
    if match_type == "name_only":
        return "tier5_name_only", reason
    source_backed = match_type in {"curated", "rbh", "orthogroup"} and not (route == "c" and match_type == "orthogroup")
    sequence_match = match_type in {"domain", "blastp"} or (route == "c" and match_type == "orthogroup")
    if source_backed:
        coverage_ok, _ = _coverage_pass(identity, query_coverage, target_coverage)
        if (ambiguity_status == "resolved_single_copy" and provider_status in _COMPLETE_PROVIDER_STATUSES
                and coverage_ok and not reason):
            return "tier1_source_backed_unique", ""
        # Below threshold the evidence is preserved, never erased: tier2 with
        # coverage_insufficient and no reason (review required).
        return "tier2_source_backed_ambiguous", ""
    if sequence_match:
        coverage_ok, _ = _coverage_pass(identity, query_coverage, target_coverage)
        if not coverage_ok:
            return _TIER6, "not_evaluated_coverage_below_threshold"
        if ambiguity_status == "resolved_single_copy" and provider_status in _COMPLETE_PROVIDER_STATUSES:
            return "tier3_sequence_match_unique", reason
        return "tier4_sequence_match_ambiguous", reason
    return _TIER6, reason


def _route_b_match_type(record: dict[str, object]) -> str:
    relationship = str(record.get("relationship", "")).strip()
    orthogroup = str(record.get("orthogroup", "")).strip()
    if relationship:
        if relationship in _RBH_RELATIONSHIPS:
            return "rbh"
        # best_hit_homology (or any unlisted label) is a best-hit claim, never
        # reciprocal orthology: blastp, and it downgrades the whole row even
        # when an orthogroup id is present.
        return "blastp"
    if orthogroup:
        return "orthogroup"
    return "blastp"


def _build_cluster_map(clusters: pd.DataFrame) -> tuple[dict[str, str], dict[str, str]]:
    """Return (source_token -> cluster_id, cluster_id -> root_source_id).

    Tokens listed in more than one cluster are excluded (ambiguous); rows that
    resolve through them fall back to ``dependency_unknown``.
    """
    token_to_clusters: dict[str, set[str]] = {}
    cluster_roots: dict[str, str] = {}
    for row in clusters.to_dict("records"):
        cluster = str(row["dependency_cluster"]).strip()
        if not cluster:
            continue
        members = [part for part in re.split(r"[;,\s]+", str(row.get("source_ids", ""))) if part]
        root = str(row.get("root_source_id", "")).strip()
        if root:
            members.append(root)
            cluster_roots.setdefault(cluster, root)
        for member in members:
            token_to_clusters.setdefault(member, set()).add(cluster)
    mapping = {token: next(iter(value)) for token, value in token_to_clusters.items() if len(value) == 1}
    return mapping, cluster_roots


def _resolve_cluster(
    source_id: str,
    cluster_map: dict[str, str] | None,
    cluster_roots: dict[str, str] | None,
) -> tuple[str, str]:
    """Return (dependency_cluster_id, root_source_id) for a source_id."""
    if not cluster_map:
        return _DEPENDENCY_UNKNOWN, ""
    tokens = _tokens(source_id, ";")
    hits = {cluster_map[token] for token in tokens if token in cluster_map}
    if len(hits) == 1:
        cluster = next(iter(hits))
        return cluster, (cluster_roots or {}).get(cluster, "")
    return _DEPENDENCY_UNKNOWN, ""


def _route_a_rows(
    hypotheses: pd.DataFrame,
    family_status: pd.DataFrame,
    bait_panel: pd.DataFrame,
    case_id: str,
    blast_evidence: pd.DataFrame | None = None,
    interval_gene_ids: set[str] | None = None,
    long_form_arg: str | None = None,
) -> list[dict[str, str]]:
    """Route A rows (dependency cluster assigned by the caller).

    ``blast_evidence`` (blast_evidence.tsv from ``blast-panel-search``) is
    joined onto matched-family rows by ``reference_accession`` (the verified
    bait accession, which is the blast ``panel_id``). When a bait was
    evaluated with a threshold-passing top hit, its identity / query coverage /
    target coverage / E-value and paralog ambiguity are transcribed into the
    bridge row. A blast hit is evidence of similarity only: it never promotes
    a matched family-status row to a higher tier. The sole exception is the
    explicit falsified-family rescue above, where an eligible hit creates one
    sequence-match row for its actual target; a missing / ``not_evaluated``
    blast row leaves the family-status evidence unchanged (#374 Criterion 3).
    """
    status_by_hypothesis: dict[str, dict[str, object]] = {}
    for row in family_status.to_dict("records"):
        status_by_hypothesis.setdefault(str(row["hypothesis_id"]), row)
    panel_by_hypothesis: dict[str, dict[str, object]] = {}
    for row in bait_panel.to_dict("records"):
        panel_by_hypothesis.setdefault(str(row["hypothesis_id"]), row)

    blast_by_accession: dict[str, dict[str, object]] = {}
    if blast_evidence is not None and not blast_evidence.empty:
        for row in blast_evidence.to_dict("records"):
            blast_by_accession.setdefault(str(row["panel_id"]).strip(), row)
    long_form_path = _discover_long_form_blast_artifact(long_form_arg or "")
    long_form_by_accession = _long_form_blast_hits(long_form_arg or "", blast_by_accession)

    hypothesis_columns = sorted(hypotheses.columns)
    status_columns = sorted(family_status.columns) if not family_status.empty else hypothesis_columns

    rows: list[dict[str, str]] = []
    for hypothesis in hypotheses.to_dict("records"):
        hypothesis_id = str(hypothesis["hypothesis_id"])
        status = status_by_hypothesis.get(hypothesis_id)
        panel = panel_by_hypothesis.get(hypothesis_id)
        verified = panel is not None
        function_source_ids = str(panel.get("function_source_ids", "")).strip() if panel else ""
        root_tokens = _tokens(function_source_ids, ";")
        verified_root = function_source_ids if len(root_tokens) == 1 else ""
        verified_root_cluster = _cluster_id(verified_root) if verified_root else _DEPENDENCY_UNKNOWN
        source_id = str(status.get("source_ids", "")).strip() if status is not None else ""
        if not source_id:
            source_id = function_source_ids
        provider_version = str(panel.get("source_provider", "")).strip() if panel else ""
        # The verified bait accession is the blast join key. Route A hypotheses
        # tables vary in what they carry: prefer the proposed accession, then
        # the verified bait panel's bait_id (the frozen accession), then the
        # gene name. The panel fallback is what makes the blast evidence join
        # work end-to-end when the hypothesis table only carries a gene name
        # (#374 Criterion 3).
        reference_accession = (str(hypothesis.get("proposed_stable_id", "")).strip()
                               or (str(panel.get("bait_id", "")).strip() if panel else "")
                               or str(hypothesis.get("reference_gene_name", "")).strip()
                               or str(hypothesis.get("reference_genes", "")).strip())
        family_id = str(hypothesis.get("gene_family", "")).strip() or hypothesis_id
        member_tokens = _tokens(hypothesis.get("member_interval_genes", ""))

        def unverified_row(gene_id: str, paralog_count: str, source: dict[str, object],
                           source_columns: list[str]) -> dict[str, str]:
            return {
                "case_id": case_id, "gene_id": gene_id, "family_id": family_id, "route": "a",
                "reference_accession": reference_accession, "source_id": "",
                "match_type": "", "identity": "", "query_coverage": "", "target_coverage": "",
                "evalue": "", "paralog_count": paralog_count, "ambiguity_status": "not_evaluated",
                "provider_status": "not_run", "provider_version": "",
                "input_fingerprint": _row_input_fingerprint(source, source_columns),
                "review_status": _REVIEW_NOT_EVALUATED, "dependency_cluster_id": _DEPENDENCY_UNKNOWN,
                "mapping_tier": _TIER6, "rank_contribution": "none",
                "reason": "not_evaluated_no_verified_bait",
            }

        if status is None:
            for gene_id in member_tokens:
                row = unverified_row(gene_id, str(len(member_tokens)), hypothesis, hypothesis_columns)
                rows.append(_finalize_row(row))
            continue

        family_status_value = str(status.get("family_status", "")).strip()
        interval_tokens = _tokens(status.get("interval_gene_ids", ""))
        gene_tokens = interval_tokens or member_tokens

        if (family_status_value in _FALSIFIED_FAMILY_STATUSES
                or family_status_value in _SEQUENCE_PROJECTION_FAMILY_STATUSES):
            # Family-index non-detection is not family evidence, but it must
            # not erase an independently verified direct projection.  The
            # projection is deliberately one row for the actual BLAST top hit
            # -- never one row per proposed/member family token.  The same
            # sequence evidence is the only way in when no reference-genome
            # index existed at all; the two cases differ in what
            # may be concluded, not in how the interval members are found, so
            # they are recorded under separate provenance keys.
            projection_is_falsified = family_status_value in _FALSIFIED_FAMILY_STATUSES
            provenance_key = ("_route_a_falsified_projection" if projection_is_falsified
                              else "_route_a_sequence_projection")
            projected_any = False
            blast = blast_by_accession.get(reference_accession)
            long_hits = long_form_by_accession.get(reference_accession, [])
            # A direct projection is safe only when the caller supplied the
            # interval-gene universe used by the canonical route-B prepare
            # stage.  Without that universe, a hand-written blast row could
            # smuggle an arbitrary identifier into Route A, so fail closed.
            target_is_interval_gene = (
                interval_gene_ids is not None
                and str(blast.get("top_hit_target_id", "")).strip() in interval_gene_ids
            ) if blast is not None else False
            # The summary row intentionally retains only top_hit_target_id for
            # compact auditing.  When the search stage supplies its separate
            # long-form threshold-passing artifact, preserve every eligible
            # interval target through the falsified-family bridge.  Sequence
            # order is not copy evidence: two or more eligible targets are all
            # emitted as tier4/paralog_ambiguous for downstream copy review.
            if (verified and blast is not None and long_form_path is not None
                    and str(blast.get("provider_status", "")).strip() == "evaluable"
                    and str(blast.get("review_status", "")).strip() == "verified"
                    and interval_gene_ids is not None):
                eligible_hits = [hit for hit in long_hits if hit["target_id"] in interval_gene_ids]
                for hit in eligible_hits:
                    identity = hit["identity"]
                    query_coverage = hit["query_coverage"]
                    target_coverage = hit["target_coverage"]
                    evalue = hit["evalue"]
                    # The long-form table preserves all threshold-passing
                    # targets, but it must not erase a conservative compact
                    # verdict.  A single retained row can still be one of
                    # several hits when the compact bait summary classified
                    # the query as ambiguous (or otherwise non-unique).
                    compact_ambiguity = str(blast.get("ambiguity_status", "")).strip()
                    ambiguity = (
                        "paralog_ambiguous"
                        if len(eligible_hits) > 1 or compact_ambiguity != "unique"
                        else "resolved_single_copy"
                    )
                    mapping_tier, tier_reason = _derive_tier(
                        route="a", match_type="blastp", ambiguity_status=ambiguity,
                        provider_status="complete", reason="", identity=identity,
                        query_coverage=query_coverage, target_coverage=target_coverage,
                    )
                    blast_fingerprint = _row_input_fingerprint(hit, sorted(hit))
                    family_status_fingerprint = _row_input_fingerprint(status, status_columns)
                    row = {
                        "case_id": case_id, "gene_id": hit["target_id"], "family_id": family_id,
                        "route": "a", "reference_accession": reference_accession,
                        "source_id": source_id, "match_type": "blastp", "identity": identity,
                        "query_coverage": query_coverage, "target_coverage": target_coverage,
                        "evalue": evalue, "paralog_count": str(max(0, len(eligible_hits) - 1)),
                        "ambiguity_status": ambiguity, "provider_status": "complete",
                        "provider_version": provider_version,
                        "input_fingerprint": hashlib.sha256(
                            f"family_status:{family_status_fingerprint}\nblast:{blast_fingerprint}".encode("utf-8")
                        ).hexdigest(),
                        "review_status": _derive_review_status(mapping_tier, ambiguity, "blastp"),
                        "dependency_cluster_id": verified_root_cluster,
                        "mapping_tier": mapping_tier, "rank_contribution": "none",
                        "reason": tier_reason,
                        provenance_key: {
                            "hypothesis_id": hypothesis_id, "family_id": family_id,
                            "family_status": family_status_value,
                            "paralog_guard_evaluated": projection_is_falsified,
                            "family_status_input_fingerprint": family_status_fingerprint,
                            "reference_accession": reference_accession,
                            "blast_panel_id": str(blast.get("panel_id", "")).strip(),
                            "blast_provider_status": str(blast.get("provider_status", "")).strip(),
                            "blast_review_status": str(blast.get("review_status", "")).strip(),
                            "blast_ambiguity_status": "ambiguous_paralog" if len(eligible_hits) > 1 else str(blast.get("ambiguity_status", "")).strip(),
                            "blast_top_hit_target_id": str(blast.get("top_hit_target_id", "")).strip(),
                            "blast_projected_target_id": hit["target_id"],
                            "blast_eligible_target_count": len(eligible_hits),
                            "blast_db_fingerprint": str(blast.get("db_fingerprint", "")).strip(),
                            "blast_interval_proteins_sha256": str(blast.get("interval_proteins_sha256", "")).strip(),
                            "blast_blastp_version": str(blast.get("blastp_version", "")).strip(),
                            "blast_input_fingerprint": blast_fingerprint,
                            "long_form_artifact": str(long_form_path.name),
                            "interval_gene_universe_supplied": True,
                            "target_in_interval_gene_universe": True,
                        },
                    }
                    rows.append(_finalize_row(row))
                    projected_any = True
                # A present long-form artifact is authoritative.  If none of
                # its targets belongs to the supplied interval universe, do
                # not fall back to top-hit-only projection.
                if projection_is_falsified or projected_any:
                    continue
            if (verified and blast is not None and target_is_interval_gene
                    and _route_a_blast_projection_eligible(blast)):
                target_gene_id = str(blast.get("top_hit_target_id", "")).strip()
                identity = str(blast.get("top_hit_identity", "")).strip()
                query_coverage = str(blast.get("top_hit_query_coverage", "")).strip()
                target_coverage = str(blast.get("top_hit_target_coverage", "")).strip()
                evalue = str(blast.get("top_hit_evalue", "")).strip()
                blast_paralog_count = str(blast.get("paralog_count", "")).strip()
                paralog_count = blast_paralog_count or "0"
                ambiguity = _derive_ambiguity(
                    provider_status="complete", reason="", identity=identity,
                    query_coverage=query_coverage, target_coverage=target_coverage,
                    matched_gene_count=1, paralog_count=paralog_count,
                    family_status=family_status_value, relationship="",
                    tandem_ambiguous=False,
                )
                # The source BLAST verdict is authoritative for uniqueness;
                # retain the existing bridge/paralog interpretation as a
                # defensive additional guard rather than relaxing it.
                if str(blast.get("ambiguity_status", "")).strip() == "ambiguous_paralog":
                    ambiguity = "paralog_ambiguous"
                mapping_tier, tier_reason = _derive_tier(
                    route="a", match_type="blastp", ambiguity_status=ambiguity,
                    provider_status="complete", reason="", identity=identity,
                    query_coverage=query_coverage, target_coverage=target_coverage,
                )
                blast_columns = sorted(blast)
                blast_fingerprint = _row_input_fingerprint(blast, blast_columns)
                family_status_fingerprint = _row_input_fingerprint(status, status_columns)
                row = {
                    "case_id": case_id, "gene_id": target_gene_id, "family_id": family_id,
                    "route": "a", "reference_accession": reference_accession,
                    "source_id": source_id, "match_type": "blastp", "identity": identity,
                    "query_coverage": query_coverage, "target_coverage": target_coverage,
                    "evalue": evalue, "paralog_count": paralog_count,
                    "ambiguity_status": ambiguity, "provider_status": "complete",
                    "provider_version": provider_version,
                    # Combine both input records so the bridge row is tied to
                    # the family-index disagreement and the BLAST evidence.
                    "input_fingerprint": hashlib.sha256(
                        f"family_status:{family_status_fingerprint}\nblast:{blast_fingerprint}".encode("utf-8")
                    ).hexdigest(),
                    "review_status": _derive_review_status(mapping_tier, ambiguity, "blastp"),
                    "dependency_cluster_id": verified_root_cluster,
                    "mapping_tier": mapping_tier, "rank_contribution": "none",
                    "reason": tier_reason,
                    provenance_key: {
                        "hypothesis_id": hypothesis_id, "family_id": family_id,
                        "family_status": family_status_value,
                        "paralog_guard_evaluated": projection_is_falsified,
                        "family_status_input_fingerprint": family_status_fingerprint,
                        "reference_accession": reference_accession,
                        "blast_panel_id": str(blast.get("panel_id", "")).strip(),
                        "blast_provider_status": str(blast.get("provider_status", "")).strip(),
                        "blast_review_status": str(blast.get("review_status", "")).strip(),
                        "blast_ambiguity_status": str(blast.get("ambiguity_status", "")).strip(),
                        "blast_top_hit_target_id": target_gene_id,
                        "blast_db_fingerprint": str(blast.get("db_fingerprint", "")).strip(),
                        "blast_interval_proteins_sha256": str(blast.get("interval_proteins_sha256", "")).strip(),
                        "blast_blastp_version": str(blast.get("blastp_version", "")).strip(),
                        "blast_input_fingerprint": blast_fingerprint,
                        "interval_gene_universe_supplied": interval_gene_ids is not None,
                        "target_in_interval_gene_universe": target_is_interval_gene,
                    },
                }
                rows.append(_finalize_row(row))
                projected_any = True
            # A falsified family with no independently eligible BLAST row stays
            # absent.  Missing/failed/ambiguous/below-threshold evidence is not
            # converted into a negative bridge row.
            if projection_is_falsified or projected_any:
                continue
            # No sequence projection was possible for a family whose index was
            # simply unavailable: the question is still unanswered, so fall
            # through to the existing not_evaluated handling rather than
            # silently dropping the hypothesis.

        if family_status_value in _UNCLASSIFIED_FAMILY_STATUSES:
            if not interval_tokens:
                continue
            for gene_id in interval_tokens:
                row = {
                    "case_id": case_id, "gene_id": gene_id, "family_id": family_id, "route": "a",
                    "reference_accession": reference_accession, "source_id": source_id,
                    "match_type": "", "identity": "", "query_coverage": "", "target_coverage": "",
                    "evalue": "", "paralog_count": str(len(interval_tokens)),
                    "ambiguity_status": "unclassified_gene",
                    "provider_status": "complete", "provider_version": provider_version,
                    "input_fingerprint": _row_input_fingerprint(status, status_columns),
                    "review_status": _REVIEW_NOT_EVALUATED,
                    "dependency_cluster_id": verified_root_cluster,
                    "mapping_tier": _TIER6, "rank_contribution": "none",
                    "reason": "not_evaluated_unclassified_gene",
                }
                rows.append(_finalize_row(row))
            continue

        if family_status_value.startswith("not_evaluated") or family_status_value == "not_evaluated":
            if family_status_value == "not_evaluated_interval_unresolved":
                reason = "not_evaluated_interval_unresolved"
            elif family_status_value == "not_evaluated_no_verified_bait":
                reason = "not_evaluated_no_verified_bait"
            elif family_status_value == "not_evaluated_family_unresolved_in_index":
                # Zero resolution of the declared identifiers in the index: an
                # index-vocabulary gap held as not-evaluated, never a negative.
                reason = "not_evaluated_family_unresolved_in_index"
            else:
                reason = "not_evaluated_provider_missing"
            for gene_id in gene_tokens:
                row = {
                    "case_id": case_id, "gene_id": gene_id, "family_id": family_id, "route": "a",
                    "reference_accession": reference_accession, "source_id": source_id,
                    "match_type": "", "identity": "", "query_coverage": "", "target_coverage": "",
                    "evalue": "", "paralog_count": str(len(gene_tokens)),
                    "ambiguity_status": "not_evaluated",
                    "provider_status": "not_run", "provider_version": provider_version,
                    "input_fingerprint": _row_input_fingerprint(status, status_columns),
                    "review_status": _REVIEW_NOT_EVALUATED,
                    "dependency_cluster_id": verified_root_cluster,
                    "mapping_tier": _TIER6, "rank_contribution": "none", "reason": reason,
                }
                rows.append(_finalize_row(row))
            continue

        if family_status_value in _MATCHED_FAMILY_STATUSES:
            if not verified:
                for gene_id in member_tokens:
                    row = unverified_row(gene_id, str(len(member_tokens)), status, status_columns)
                    rows.append(_finalize_row(row))
                continue
            match_type = ("curated" if family_status_value in {"concrete_interval_gene_matched",
                                                               "family_present_copy_unresolved"}
                          else "blastp")
            # Blast evidence (blast-panel-search) joined by reference_accession
            # == blast panel_id. Only an evaluable row with a threshold-passing
            # top hit contributes identity/coverage/E-value; the blast paralog
            # count and ambiguity are additional copy-level signals. A missing
            # or not_evaluated blast row leaves the family-status evidence
            # untouched (a blast provider gap is not a miss).
            blast = blast_by_accession.get(reference_accession)
            blast_applicable = (
                blast is not None
                and str(blast.get("provider_status", "")).strip() == "evaluable"
                and str(blast.get("top_hit_target_id", "")).strip() != ""
            )
            identity = str(blast.get("top_hit_identity", "")).strip() if blast_applicable else ""
            query_coverage = str(blast.get("top_hit_query_coverage", "")).strip() if blast_applicable else ""
            target_coverage = str(blast.get("top_hit_target_coverage", "")).strip() if blast_applicable else ""
            evalue = str(blast.get("top_hit_evalue", "")).strip() if blast_applicable else ""
            blast_paralog_count = str(blast.get("paralog_count", "")).strip() if blast is not None else ""
            paralog_count = blast_paralog_count or str(len(gene_tokens))
            for gene_id in gene_tokens:
                ambiguity = _derive_ambiguity(
                    provider_status="complete", reason="", identity=identity, query_coverage=query_coverage,
                    target_coverage=target_coverage, matched_gene_count=len(gene_tokens),
                    paralog_count=paralog_count, family_status=family_status_value, relationship="",
                    tandem_ambiguous=False,
                )
                # A blast ambiguous_paralog verdict is copy-level ambiguity even
                # when the family-status gene count is one.
                if (blast_applicable and ambiguity == "resolved_single_copy"
                        and str(blast.get("ambiguity_status", "")).strip() == "ambiguous_paralog"):
                    ambiguity = "paralog_ambiguous"
                mapping_tier, tier_reason = _derive_tier(
                    route="a", match_type=match_type, ambiguity_status=ambiguity,
                    provider_status="complete", reason="", identity=identity, query_coverage=query_coverage,
                    target_coverage=target_coverage,
                )
                row = {
                    "case_id": case_id, "gene_id": gene_id, "family_id": family_id, "route": "a",
                    "reference_accession": reference_accession, "source_id": source_id,
                    "match_type": match_type, "identity": identity, "query_coverage": query_coverage,
                    "target_coverage": target_coverage,
                    "evalue": evalue, "paralog_count": paralog_count, "ambiguity_status": ambiguity,
                    "provider_status": "complete", "provider_version": provider_version,
                    "input_fingerprint": _row_input_fingerprint(status, status_columns),
                    "review_status": _derive_review_status(mapping_tier, ambiguity, match_type),
                    "dependency_cluster_id": verified_root_cluster, "mapping_tier": mapping_tier,
                    "rank_contribution": "none", "reason": tier_reason,
                }
                rows.append(_finalize_row(row))
            continue

        # Unknown family status: fail closed to the unverified-bait state.
        for gene_id in member_tokens:
            row = unverified_row(gene_id, str(len(member_tokens)), status, status_columns)
            rows.append(_finalize_row(row))

    return rows


def _route_a_interval_family_projection_rows(
    hypotheses: pd.DataFrame,
    family_status: pd.DataFrame,
    bait_panel: pd.DataFrame,
    family_index: pd.DataFrame,
    case_id: str,
) -> list[dict[str, str]]:
    """Project verified Route A families onto a frozen interval family index.

    This is deliberately positive-only.  An ``interval_candidates`` index can
    show that an interval protein resembles a named family, but it cannot show
    that the family is absent from the genome.  Therefore only the Route A
    state that explicitly says the reusable reference-genome index was missing
    is eligible here; no match produces no bridge row and changes no status.
    """
    if family_index.empty:
        return []

    status_by_hypothesis = {
        str(row["hypothesis_id"]): row for row in family_status.to_dict("records")
    }
    panel_by_hypothesis = {
        str(row["hypothesis_id"]): row for row in bait_panel.to_dict("records")
    }
    index_columns = sorted(family_index.columns)
    hypothesis_columns = sorted(hypotheses.columns)
    status_columns = sorted(family_status.columns)

    def labels(record: dict[str, object]) -> list[str]:
        values: list[str] = []
        for column in (
            "best_subject_symbol", "best_subject_stem", "family_symbols", "family_stems",
            "discovery_bucket_id", "evolutionary_family_id",
        ):
            values.extend(_tokens(record.get(column, "")))
        return list(dict.fromkeys(value for value in values if value))

    output: list[dict[str, str]] = []
    for hypothesis in hypotheses.to_dict("records"):
        hypothesis_id = str(hypothesis["hypothesis_id"])
        status = status_by_hypothesis.get(hypothesis_id)
        panel = panel_by_hypothesis.get(hypothesis_id)
        if status is None or panel is None:
            continue
        if str(status.get("family_status", "")).strip() != \
                "not_evaluated_missing_reference_genome_family_index":
            continue

        family_id = str(hypothesis.get("gene_family", "")).strip() or hypothesis_id
        query_labels = list(dict.fromkeys(value for value in (
            family_id,
            str(hypothesis.get("reference_gene_name", "")).strip(),
            str(hypothesis.get("reference_genes", "")).strip(),
        ) if value and _positive_pi_family_label(value)))
        if not query_labels:
            continue

        matches: list[tuple[dict[str, object], str, str]] = []
        for index_row in family_index.to_dict("records"):
            if str(index_row.get("index_scope", "")).strip() != "interval_candidates":
                continue
            if str(index_row.get("assignment_status", "")).strip() != "assigned":
                continue
            gene_id = str(index_row.get("gene_id", "")).strip()
            if not gene_id:
                continue
            identity = str(index_row.get("best_identity", "")).strip()
            query_coverage = str(index_row.get("best_query_coverage", "")).strip()
            coverage_ok, applicable = _coverage_pass(identity, query_coverage, "")
            if not coverage_ok or not applicable:
                continue
            matched_pair = next((
                (query_label, index_label)
                for query_label in query_labels
                for index_label in labels(index_row)
                if families_match_expanded(query_label, index_label)
            ), None)
            if matched_pair is not None:
                matches.append((index_row, matched_pair[0], matched_pair[1]))

        if not matches:
            continue
        # One family may have several interval copies.  Sequence order or best
        # bitscore is not copy evidence, so retain the complete matching set.
        unique_matches: dict[str, tuple[dict[str, object], str, str]] = {}
        for match in matches:
            unique_matches.setdefault(str(match[0]["gene_id"]).strip(), match)
        matched_count = len(unique_matches)
        ambiguity = "resolved_single_copy" if matched_count == 1 else "paralog_ambiguous"
        paralog_count = str(max(0, matched_count - 1))
        function_sources = str(panel.get("function_source_ids", "")).strip()
        source_id = str(status.get("source_ids", "")).strip() or function_sources
        root_tokens = _tokens(function_sources)
        dependency_cluster = _cluster_id(function_sources) if len(root_tokens) == 1 else _DEPENDENCY_UNKNOWN
        reference_accession = (
            str(hypothesis.get("proposed_stable_id", "")).strip()
            or str(panel.get("bait_id", "")).strip()
            or str(hypothesis.get("reference_gene_name", "")).strip()
        )
        hypothesis_fingerprint = _row_input_fingerprint(hypothesis, hypothesis_columns)
        status_fingerprint = _row_input_fingerprint(status, status_columns)
        for gene_id, (index_row, query_label, index_label) in unique_matches.items():
            identity = str(index_row.get("best_identity", "")).strip()
            query_coverage = str(index_row.get("best_query_coverage", "")).strip()
            evalue = str(index_row.get("best_evalue", "")).strip()
            mapping_tier, tier_reason = _derive_tier(
                route="a", match_type="blastp", ambiguity_status=ambiguity,
                provider_status="complete", reason="", identity=identity,
                query_coverage=query_coverage, target_coverage="",
            )
            index_fingerprint = _row_input_fingerprint(index_row, index_columns)
            row = {
                "case_id": case_id, "gene_id": gene_id, "family_id": family_id, "route": "a",
                "reference_accession": reference_accession, "source_id": source_id,
                "match_type": "blastp", "identity": identity,
                "query_coverage": query_coverage, "target_coverage": "", "evalue": evalue,
                "paralog_count": paralog_count, "ambiguity_status": ambiguity,
                "provider_status": "complete", "provider_version": "route_b0_interval_family_index",
                "input_fingerprint": hashlib.sha256(
                    f"hypothesis:{hypothesis_fingerprint}\nstatus:{status_fingerprint}\n"
                    f"family_index:{index_fingerprint}".encode("utf-8")
                ).hexdigest(),
                "review_status": _derive_review_status(mapping_tier, ambiguity, "blastp"),
                "dependency_cluster_id": dependency_cluster, "mapping_tier": mapping_tier,
                "rank_contribution": "none", "reason": tier_reason,
                "_route_a_interval_family_projection": {
                    "hypothesis_id": hypothesis_id, "family_id": family_id, "gene_id": gene_id,
                    "query_label": query_label, "interval_family_label": index_label,
                    "index_scope": "interval_candidates", "assignment_status": "assigned",
                    "matched_interval_copy_count": matched_count,
                    "best_subject": str(index_row.get("best_subject", "")).strip(),
                    "best_identity": identity, "best_query_coverage": query_coverage,
                    "best_evalue": evalue, "family_index_input_fingerprint": index_fingerprint,
                    "positive_only": True, "absence_inference_allowed": False,
                },
            }
            output.append(_finalize_row(row))
    return output


def _route_b_rows(
    orthology: pd.DataFrame,
    atlas: pd.DataFrame,
    case_id: str,
) -> tuple[list[dict[str, str]], int]:
    rows: list[dict[str, str]] = []
    skipped = 0
    atlas_rows: dict[str, dict[str, object]] = {}
    if not atlas.empty:
        atlas_rows = {str(row["gene_id"]): row for row in atlas.to_dict("records")}
    orthology_columns = sorted(orthology.columns)
    atlas_columns = sorted(atlas.columns) if not atlas.empty else []

    for record in orthology.to_dict("records"):
        gene_id = str(record["gene_id"])
        atlas_record = atlas_rows.get(gene_id, {})
        family_id = (str(atlas_record.get("evolutionary_family_id", "")).strip()
                     or str(atlas_record.get("family_id", "")).strip()
                     or str(record.get("orthogroup", "")).strip())
        if not family_id:
            skipped += 1
            continue
        match_type = _route_b_match_type(record)
        paralog_count = str(atlas_record.get("focal_paralog_count", "")) if not atlas.empty else ""
        source_id = str(record.get("source_id", "")).strip() if "source_id" in orthology.columns else ""
        # Without an explicit source_id, the orthogroup/HOG id is the identity
        # of the homology claim: rows derived from the same orthogroup then
        # resolve to the same dependency cluster through the supplied clusters
        # table (never auto-promoted without it).
        if not source_id:
            source_id = orthogroup_dependency_token(record.get("orthogroup", "")) if str(record.get("orthogroup", "")).strip() else ""
        provider_status = "complete"
        ambiguity = _derive_ambiguity(
            provider_status=provider_status, reason="",
            identity=str(record.get("percent_identity", "")),
            query_coverage=str(record.get("query_coverage", "")),
            target_coverage=str(record.get("target_coverage", "")),
            matched_gene_count=1, paralog_count=paralog_count, family_status="",
            relationship=str(record.get("relationship", "")).strip(),
            tandem_ambiguous=_tandem_marker(atlas_record),
        )
        mapping_tier, tier_reason = _derive_tier(
            route="b", match_type=match_type, ambiguity_status=ambiguity,
            provider_status=provider_status, reason="",
            identity=str(record.get("percent_identity", "")),
            query_coverage=str(record.get("query_coverage", "")),
            target_coverage=str(record.get("target_coverage", "")),
        )
        row = {
            "case_id": case_id, "gene_id": gene_id, "family_id": family_id, "route": "b",
            "reference_accession": str(record.get("reference_gene_id", "")).strip(), "source_id": source_id,
            "match_type": match_type, "identity": str(record.get("percent_identity", "")),
            "query_coverage": str(record.get("query_coverage", "")),
            "target_coverage": str(record.get("target_coverage", "")),
            "evalue": str(record.get("evalue", "")),
            "paralog_count": paralog_count, "ambiguity_status": ambiguity,
            "provider_status": provider_status, "provider_version": str(record.get("provider", "")).strip(),
            "input_fingerprint": _row_input_fingerprint(record, orthology_columns),
            "review_status": _derive_review_status(mapping_tier, ambiguity, match_type),
            "dependency_cluster_id": _DEPENDENCY_UNKNOWN, "mapping_tier": mapping_tier,
            "rank_contribution": "none", "reason": tier_reason,
        }
        rows.append(_finalize_row(row))

    if not atlas.empty:
        for gene_id in sorted(atlas_rows):
            atlas_record = atlas_rows[gene_id]
            evolutionary = str(atlas_record.get("evolutionary_family_id", "")).strip()
            domain_ok = bool(str(atlas_record.get("domain_architecture", "")).strip()
                             and str(atlas_record.get("domain_ids", "")).strip()
                             and str(atlas_record.get("domain_source", "")).strip())
            orthogroup = str(atlas_record.get("orthogroup", "")).strip()
            bucket = str(atlas_record.get("family_id", "")).strip() or str(atlas_record.get("family_label", "")).strip()
            axes: list[tuple[str, str]] = []
            if evolutionary:
                axes.append(("curated", evolutionary))
            if domain_ok:
                axes.append(("domain", str(atlas_record.get("domain_architecture", "")).strip()))
            if orthogroup:
                axes.append(("orthogroup", orthogroup))
            if not (evolutionary or domain_ok or orthogroup) and bucket:
                axes.append(("name_only", bucket))
            for match_type, family_id in axes:
                paralog_count = str(atlas_record.get("focal_paralog_count", "")).strip()
                ambiguity = _derive_ambiguity(
                    provider_status="complete", reason="", identity="", query_coverage="",
                    target_coverage="", matched_gene_count=1, paralog_count=paralog_count,
                    family_status="", relationship="", tandem_ambiguous=_tandem_marker(atlas_record),
                )
                mapping_tier, tier_reason = _derive_tier(
                    route="b", match_type=match_type, ambiguity_status=ambiguity,
                    provider_status="complete", reason="", identity="", query_coverage="",
                    target_coverage="",
                )
                row = {
                    "case_id": case_id, "gene_id": gene_id, "family_id": family_id, "route": "b",
                    "reference_accession": str(atlas_record.get("orthology_reference_gene", "")).split(";")[0].strip(),
                    "source_id": orthogroup_dependency_token(family_id) if match_type == "orthogroup" else "",
                    "match_type": match_type, "identity": "", "query_coverage": "",
                    "target_coverage": "", "evalue": "", "paralog_count": paralog_count,
                    "ambiguity_status": ambiguity, "provider_status": "complete", "provider_version": "",
                    "input_fingerprint": _row_input_fingerprint(atlas_record, atlas_columns),
                    "review_status": _derive_review_status(mapping_tier, ambiguity, match_type),
                    "dependency_cluster_id": _DEPENDENCY_UNKNOWN, "mapping_tier": mapping_tier,
                    "rank_contribution": "none", "reason": tier_reason,
                }
                rows.append(_finalize_row(row))

    return rows, skipped


def _tandem_marker(atlas_record: dict[str, object]) -> bool:
    value = str(atlas_record.get("tandem_expansion", "")).strip().casefold()
    return value in {"true", "1", "yes"}


def _route_c_rows(
    match: pd.DataFrame,
    merged: pd.DataFrame,
    case_id: str,
) -> tuple[list[dict[str, str]], int, int]:
    merged_lookup: dict[tuple[str, str], dict[str, object]] = {}
    for row in merged.to_dict("records"):
        merged_lookup[(str(row["reference_gene_id"]), str(row["rank"]))] = row
    match_columns = sorted(match.columns)

    rows: list[dict[str, str]] = []
    skipped_no_family = 0
    skipped_empty_genes = 0
    for record in match.to_dict("records"):
        reference = str(record["reference_gene_id"])
        rank = str(record["rank"])
        family_id = (str(record.get("family_id", "")).strip()
                     or str(record.get("orthogroup_id", "")).strip())
        if not family_id:
            skipped_no_family += 1
            continue
        merged_record = merged_lookup.get((reference, rank), {})
        database_version = str(merged_record.get("database_version", "")).strip()
        model_version = str(merged_record.get("model_version", "")).strip()
        basis = {part.strip() for part in str(record.get("match_basis", "")).split(";") if part.strip()}
        match_type = "orthogroup" if "same_orthogroup" in basis else "name_only"
        gene_tokens = _tokens(record.get("matched_interval_gene", ""))
        if not gene_tokens:
            # A match row that maps to no interval gene is not bridge evidence
            # for any gene; record the skip instead of emitting a row.
            skipped_empty_genes += 1
            continue
        paralog_count = str(len(gene_tokens))
        dependency_cluster_id = _dependency_root(database_version)
        for gene_id in gene_tokens:
            ambiguity = _derive_ambiguity(
                provider_status="complete", reason="", identity="", query_coverage="",
                target_coverage="", matched_gene_count=len(gene_tokens), paralog_count=paralog_count,
                family_status="", relationship="", tandem_ambiguous=False,
            )
            mapping_tier, tier_reason = _derive_tier(
                route="c", match_type=match_type, ambiguity_status=ambiguity,
                provider_status="complete", reason="", identity="", query_coverage="",
                target_coverage="",
            )
            row = {
                "case_id": case_id, "gene_id": gene_id, "family_id": family_id, "route": "c",
                "reference_accession": reference, "source_id": database_version,
                "match_type": match_type, "identity": "", "query_coverage": "", "target_coverage": "",
                "evalue": "", "paralog_count": paralog_count, "ambiguity_status": ambiguity,
                "provider_status": "complete", "provider_version": model_version,
                "input_fingerprint": _row_input_fingerprint(record, match_columns),
                "review_status": _derive_review_status(mapping_tier, ambiguity, match_type),
                "dependency_cluster_id": dependency_cluster_id, "mapping_tier": mapping_tier,
                "rank_contribution": "none", "reason": tier_reason,
            }
            rows.append(_finalize_row(row))

    return rows, skipped_no_family, skipped_empty_genes


def _dependency_root(database_version: str) -> str:
    if database_version:
        return _cluster_id(database_version)
    return _DEPENDENCY_UNKNOWN


def family_candidate_bridge(
    out_dir_arg: str,
    route_a_hypotheses_arg: str | None = None,
    route_a_family_status_arg: str | None = None,
    route_a_bait_panel_arg: str | None = None,
    route_b_orthology_arg: str | None = None,
    route_b_atlas_arg: str | None = None,
    route_c_retrieval_arg: str | None = None,
    route_c_match_arg: str | None = None,
    dependency_clusters_arg: str | None = None,
    genes_arg: str | None = None,
    case_id: str = "",
    route_a_blast_arg: str | None = None,
    route_a_blast_eligible_arg: str | None = None,
    route_b_family_index_arg: str | None = None,
    route_a_projection_panel_arg: str | None = None,
) -> dict[str, object]:
    """Normalize Route A / B / C family evidence into one reporting table.

    ``route_a_blast_arg`` is the compact ``blast_evidence.tsv`` written by
    ``blast-panel-search`` for the Route A verified baits; its identity /
    coverage / E-value / paralog ambiguity are transcribed onto matched Route A
    rows by reference_accession (issue #374 Criterion 3). When the sibling
    long-form eligible-hit artifact is present, every threshold-passing target
    is retained for a falsified family; the compact top hit is never used as a
    copy selector (issue #532). A blast hit is similarity evidence only for an
    existing matched family row; falsified-family rows remain reviewable and
    carry no rank contribution.
    """
    case_id = case_id.strip()
    if not case_id:
        raise ValueError("family-candidate-bridge requires a non-empty --case-id")
    if not (route_a_hypotheses_arg or route_b_orthology_arg or route_b_atlas_arg or route_c_match_arg):
        raise ValueError("family-candidate-bridge requires at least one route input "
                         "(--route-a-hypotheses / --route-b-orthology / --route-b-atlas / --route-c-match)")

    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)

    hypotheses = _read(route_a_hypotheses_arg, {"hypothesis_id", "member_interval_genes"}, "route A hypotheses")
    family_status = _read(route_a_family_status_arg, {"hypothesis_id"}, "route A family status")
    bait_panel = _read(route_a_bait_panel_arg, {"hypothesis_id"}, "route A bait panel")
    projection_panel = _read(route_a_projection_panel_arg, {"hypothesis_id"}, "Route A projection bait panel")
    blast_evidence = _read(route_a_blast_arg, {"panel_id"}, "route A blast evidence")
    long_form_arg = route_a_blast_eligible_arg or route_a_blast_arg
    long_form_path = _discover_long_form_blast_artifact(long_form_arg)
    long_form_by_accession = _long_form_blast_hits(long_form_arg, {
        str(row["panel_id"]).strip(): row for row in blast_evidence.to_dict("records")
    }) if long_form_path is not None else {}
    orthology = _read(route_b_orthology_arg, {"gene_id", "reference_gene_id"}, "route B orthology")
    atlas = _read(route_b_atlas_arg, {"gene_id"}, "route B atlas")
    route_b_family_index = _read(
        route_b_family_index_arg, {"gene_id", "assignment_status", "index_scope"},
        "Route B family index",
    )
    merged = _read(route_c_retrieval_arg, {"reference_gene_id", "rank"}, "route C retrieval")
    match = _read(route_c_match_arg, {"reference_gene_id", "rank", "matched_interval_gene"}, "route C interval match")
    clusters = _read(dependency_clusters_arg, {"dependency_cluster", "source_ids"}, "dependency clusters")
    genes_frame = _read(genes_arg, {"gene_id"}, "genes")

    cluster_map, cluster_roots = _build_cluster_map(clusters) if not clusters.empty else ({}, {})

    rows: list[dict[str, str]] = []
    skipped = {"route_b_rows_skipped_no_family": 0, "route_c_rows_skipped_no_family": 0}
    cluster_roots_by_id: dict[str, str] = dict(cluster_roots)
    for panel_record in bait_panel.to_dict("records"):
        root_tokens = _tokens(panel_record.get("function_source_ids", ""), ";")
        if len(root_tokens) == 1:
            cluster_roots_by_id.setdefault(_cluster_id(root_tokens[0]), root_tokens[0])

    interval_gene_ids = (
        {str(value).strip() for value in genes_frame["gene_id"] if str(value).strip()}
        if genes_arg is not None else None
    )
    route_a_rows = _route_a_rows(
        hypotheses, family_status, bait_panel, case_id, blast_evidence, interval_gene_ids,
        long_form_arg,
    )
    interval_projection_rows = _route_a_interval_family_projection_rows(
        hypotheses, family_status, bait_panel, route_b_family_index, case_id,
    )
    # The same hypothesis/reference may reach an interval copy both by name
    # and by sequence. Prefer its sequence row, while retaining
    # independent references and hypotheses about the same interval gene.
    def projection_key(row: dict[str, object], provenance_field: str) -> tuple[str, ...] | None:
        provenance = row.get(provenance_field)
        if not isinstance(provenance, dict):
            return None
        key = (str(row.get("case_id", "")), str(row.get("gene_id", "")),
               str(row.get("reference_accession", "")), str(provenance.get("hypothesis_id", "")))
        return key if all(key) else None

    sequence_projected_keys = {
        key for row in route_a_rows
        if (key := projection_key(row, "_route_a_sequence_projection")) is not None
    }
    deduplicated_interval_projections = sorted(
        str(row.get("gene_id", "")) for row in interval_projection_rows
        if projection_key(row, "_route_a_interval_family_projection") in sequence_projected_keys
    )
    interval_projection_rows = [
        row for row in interval_projection_rows
        if projection_key(row, "_route_a_interval_family_projection") not in sequence_projected_keys
    ]
    route_a_rows.extend(interval_projection_rows)
    from .partial_alignment_review import partial_review_hits
    partial_hits = partial_review_hits(route_a_blast_arg, projection_panel.to_dict("records"), out) if not projection_panel.empty else {}
    partial_path = out / "partial_alignment_review.json"
    review_projection_rows = _route_a_review_projection_rows(
        projection_panel, blast_evidence, long_form_by_accession, interval_gene_ids, case_id,
        partial_hits, hypotheses.to_dict("records"),
    )
    route_b_rows, route_b_skipped = _route_b_rows(orthology, atlas, case_id)
    route_c_rows, route_c_skipped_no_family, route_c_skipped_empty_genes = _route_c_rows(match, merged, case_id)
    skipped["route_b_rows_skipped_no_family"] = route_b_skipped
    skipped["route_c_rows_skipped_no_family"] = route_c_skipped_no_family
    skipped["route_c_rows_skipped_empty_interval_genes"] = route_c_skipped_empty_genes

    # Dependency resolution. Route A roots come from the verified-bait
    # function_source_ids (exactly one value) or the clusters table; Route B
    # resolves only through the clusters table (never auto-promoted); Route C
    # always roots on the declared database_version.
    for row in route_a_rows:
        if cluster_map:
            resolved_cluster, root = _resolve_cluster(row["source_id"], cluster_map, cluster_roots_by_id)
            if resolved_cluster != _DEPENDENCY_UNKNOWN:
                row["dependency_cluster_id"] = resolved_cluster
            else:
                root = cluster_roots_by_id.get(row["dependency_cluster_id"], "")
        else:
            # Rows already carry the cluster derived from the verified-bait
            # function_source_ids (exactly one value) in _route_a_rows; the
            # family-status source_ids column is display-only and never
            # promoted to a root on its own.
            root = cluster_roots_by_id.get(row["dependency_cluster_id"], "")
        if row["dependency_cluster_id"] != _DEPENDENCY_UNKNOWN:
            cluster_roots_by_id[row["dependency_cluster_id"]] = root
    for row in route_b_rows:
        row["dependency_cluster_id"], root = _resolve_cluster(row["source_id"], cluster_map, cluster_roots_by_id)
        if row["dependency_cluster_id"] != _DEPENDENCY_UNKNOWN:
            cluster_roots_by_id[row["dependency_cluster_id"]] = root
    for row in route_c_rows:
        row["dependency_cluster_id"] = _dependency_root(str(row["source_id"]))
        if row["dependency_cluster_id"] != _DEPENDENCY_UNKNOWN:
            cluster_roots_by_id[row["dependency_cluster_id"]] = str(row["source_id"])

    rows = route_a_rows + route_b_rows + route_c_rows
    rows = [_finalize_row(row) for row in rows]

    if rows:
        frame = pd.DataFrame(rows)[BRIDGE_COLUMNS]
    else:
        frame = pd.DataFrame(columns=BRIDGE_COLUMNS)
    frame = frame.drop_duplicates().sort_values(
        ["case_id", "gene_id", "family_id", "route", "match_type", "dependency_cluster_id", "source_id"],
        kind="stable",
    )
    bridge_path = out / "family_candidate_bridge.tsv"
    frame.to_csv(bridge_path, sep="\t", index=False, lineterminator="\n")
    review_projection_path = out / "route_a_review_candidates.tsv"
    pd.DataFrame(review_projection_rows, columns=ROUTE_A_REVIEW_CANDIDATE_COLUMNS).to_csv(
        review_projection_path, sep="\t", index=False, lineterminator="\n"
    )
    review_manifest_path = out / "route_a_review_candidates_manifest.json"

    input_records: dict[str, object] = {}
    for name, arg in (
        ("route_a_hypotheses", route_a_hypotheses_arg), ("route_a_family_status", route_a_family_status_arg),
        ("route_a_bait_panel", route_a_bait_panel_arg), ("route_a_blast", route_a_blast_arg),
        ("route_a_projection_panel", route_a_projection_panel_arg),
        ("route_a_partial_alignment_review", str(partial_path) if partial_hits else None),
        ("route_a_blast_hits", route_a_blast_eligible_arg),
        ("route_b_orthology", route_b_orthology_arg),
        ("route_b_family_index", route_b_family_index_arg),
        ("route_b_atlas", route_b_atlas_arg), ("route_c_retrieval", route_c_retrieval_arg),
        ("route_c_match", route_c_match_arg), ("dependency_clusters", dependency_clusters_arg),
        ("genes", genes_arg),
    ):
        if arg:
            fingerprint = file_fingerprint(runtime_path(arg), full_hash=True)
            input_records[name] = {
                "present": True, "path": str(runtime_path(arg)),
                "sha256": fingerprint.get("sha256"), "size_bytes": fingerprint.get("size_bytes"),
                "mtime_ns": fingerprint.get("mtime_ns"), "method": fingerprint.get("method"),
                "strength": fingerprint.get("strength"),
            }
        else:
            input_records[name] = {"present": False, "sha256": None, "size_bytes": None,
                                   "mtime_ns": None, "method": None, "strength": None}

    if partial_hits:
        partial_inputs = json.loads(partial_path.read_text(encoding="utf-8"))["inputs"]
        for name, entry in partial_inputs.items():
            fingerprint = file_fingerprint(runtime_path(entry["path"]), full_hash=True)
            if fingerprint["sha256"] != entry["sha256"]:
                raise ValueError("partial alignment input changed during projection")
            input_records["route_a_partial_" + name] = {"present": True, **fingerprint}

    def _counts(column: str) -> dict[str, int]:
        if frame.empty:
            return {}
        return {str(key): int(value) for key, value in frame[column].value_counts().sort_index().items()}

    falsified_family_projections = []
    for row in route_a_rows:
        provenance = row.get("_route_a_falsified_projection")
        if isinstance(provenance, dict):
            falsified_family_projections.append({
                **provenance,
                "bridge_output_fingerprint": str(row.get("output_fingerprint", "")),
            })
    falsified_family_projections.sort(
        key=lambda value: (str(value.get("hypothesis_id", "")), str(value.get("blast_panel_id", "")),
                           str(value.get("blast_top_hit_target_id", "")))
    )
    interval_family_projections = []
    for row in route_a_rows:
        provenance = row.get("_route_a_interval_family_projection")
        if isinstance(provenance, dict):
            interval_family_projections.append({
                **provenance,
                "bridge_output_fingerprint": str(row.get("output_fingerprint", "")),
            })
    interval_family_projections.sort(
        key=lambda value: (str(value.get("hypothesis_id", "")), str(value.get("gene_id", "")))
    )
    sequence_family_projections = []
    for row in route_a_rows:
        provenance = row.get("_route_a_sequence_projection")
        if isinstance(provenance, dict):
            sequence_family_projections.append({
                **provenance,
                "projected_gene_id": str(row.get("gene_id", "")),
                "mapping_tier": str(row.get("mapping_tier", "")),
                "ambiguity_status": str(row.get("ambiguity_status", "")),
                "bridge_output_fingerprint": str(row.get("output_fingerprint", "")),
            })
    sequence_family_projections.sort(
        key=lambda value: (str(value.get("hypothesis_id", "")), str(value.get("projected_gene_id", "")))
    )

    cluster_entries: dict[str, dict[str, object]] = {}
    for record in frame.to_dict("records"):
        cluster_id = str(record["dependency_cluster_id"])
        entry = cluster_entries.setdefault(cluster_id, {
            "dependency_cluster_id": cluster_id, "dependency_status": "",
            "root_source_id": "", "routes": set(), "gene_ids": set(), "row_count": 0,
        })
        root = cluster_roots_by_id.get(cluster_id, "")
        if root:
            entry["root_source_id"] = root
        entry["routes"].add(str(record["route"]))
        entry["gene_ids"].add(str(record["gene_id"]))
        entry["row_count"] += 1

    dependency_clusters = []
    for cluster_id in sorted(cluster_entries):
        entry = cluster_entries[cluster_id]
        root = str(entry["root_source_id"])
        dependency_clusters.append({
            "dependency_cluster_id": cluster_id,
            "dependency_status": "root_source_known" if root else "dependency_unknown",
            "root_source_id": root,
            "routes": ";".join(sorted(entry["routes"])),
            "gene_count": len(entry["gene_ids"]),
            "row_count": entry["row_count"],
        })

    manifest: dict[str, object] = {
        "schema_version": BRIDGE_SCHEMA_VERSION,
        "contract": "family-candidate-bridge-1.0.0",
        "command": "family-candidate-bridge",
        "case_id": case_id,
        "inputs": input_records,
        "row_counts": {
            "total": int(len(frame)),
            "by_route": _counts("route"),
            "by_mapping_tier": _counts("mapping_tier"),
            "by_ambiguity_status": _counts("ambiguity_status"),
            "by_provider_status": _counts("provider_status"),
            "by_review_status": _counts("review_status"),
        },
        "dependency_clusters": dependency_clusters,
        "skipped": skipped,
        "route_a_blast": {
            "input_present": bool(route_a_blast_arg),
            "evidence_rows": int(len(blast_evidence)) if not blast_evidence.empty else 0,
            "long_form_input_present": long_form_path is not None,
            "long_form_input": str(long_form_path.name) if long_form_path is not None else None,
            "long_form_input_sha256": sha256(long_form_path) if long_form_path is not None else None,
            "long_form_evidence_rows": int(sum(len(rows) for rows in long_form_by_accession.values())),
            # Route A bridge rows that carried an identity / coverage / evalue
            # transcribed from an evaluable blast-panel-search top hit.
            "rows_joined": int(sum(1 for row in route_a_rows
                                   if str(row.get("identity", "")).strip() != "")),
            "falsified_family_projection_count": len(falsified_family_projections),
            "falsified_family_projections": falsified_family_projections,
            "falsified_family_projection_policy": (
                "family-index falsification remains auditable disagreement; every independently "
                "verified, evaluable, threshold-passing long-form interval target is emitted, "
                "and multiple targets remain tier4_sequence_match_ambiguous/paralog_ambiguous "
                "with no top-hit copy selection"
            ),
        },
        "route_a_interval_family_projection": {
            "input_present": bool(route_b_family_index_arg),
            "projection_count": len(interval_family_projections),
            "projections": interval_family_projections,
            "deduplicated_into_sequence_projection": deduplicated_interval_projections,
            "policy": (
                "after the independent Route A and Route B freezes, a verified Route A family may use "
                "Route B0's interval_candidates family index as positive similarity evidence. Every "
                "matching interval copy is retained; no-hit is not absence, does not falsify the family, "
                "and adds no negative row or score"
            ),
        },
        "route_a_sequence_family_projection": {
            "projection_count": len(sequence_family_projections),
            "projections": sequence_family_projections,
            "paralog_guard_evaluated": False if sequence_family_projections else None,
            "policy": (
                "when no reference-genome family index existed, a verified Route A bait still reaches "
                "its interval members through the frozen blast-panel-search evidence, which does not "
                "depend on either side using the same family or gene symbol. Positive-only: no hit is "
                "not absence and never falsifies the family. Every threshold-passing interval copy is "
                "retained, the genome-wide paralog guard stays unevaluated, and no row carries rank"
            ),
        },
        "route_a_review_projection": {
            "input_present": bool(route_a_projection_panel_arg),
            "candidate_count": len(review_projection_rows),
            "policy": (
                "identity/species/sequence-verified but source-claim-unbound Route A baits may project "
                "through evaluable BLAST hits inside the supplied interval universe; significant partial "
                "alignments retain frozen identity/E-value/minimum-span gates and remain unresolved review; "
                "these rows are review_only and never enter supported_member_union"
            ),
        },
        "rank_contribution": "none",
        "thresholds": {"min_identity_fraction": MIN_IDENTITY_FRACTION,
                       "min_coverage_fraction": MIN_COVERAGE_FRACTION},
        "outputs": {
            "family_candidate_bridge": "family_candidate_bridge.tsv",
            "route_a_review_candidates": "route_a_review_candidates.tsv",
            "route_a_review_candidates_manifest": "route_a_review_candidates_manifest.json",
            "manifest": "family_candidate_bridge_manifest.json",
            "family_candidate_bridge_tsv_sha256": sha256(bridge_path),
            "route_a_review_candidates_sha256": sha256(review_projection_path),
        },
    }
    review_manifest = {
        "contract": "route-review-candidates-1.0.0",
        "stage": "route_a_review_candidate_freeze",
        "case_id": case_id,
        "output": review_projection_path.name,
        "artifact_sha256": sha256(review_projection_path),
        "inputs": {name: entry for name, entry in input_records.items()
                   if (name in {"route_a_hypotheses", "route_a_projection_panel", "route_a_blast", "route_a_blast_hits", "genes"}
                       or name.startswith("route_a_partial_"))
                   and bool(entry.get("present"))},
        "rank_contribution": "none",
        "support_policy": "review_required only; canonical support and supported_member_union are unchanged",
    }
    review_manifest_path.write_text(json.dumps(review_manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if genes_arg is not None:
        genes = sorted({str(value).strip() for value in genes_frame["gene_id"] if str(value).strip()})
        bridge_gene_ids = set(frame["gene_id"].astype(str)) if not frame.empty else set()
        manifest["genes_without_bridge_rows"] = sorted(set(genes) - bridge_gene_ids)

    manifest_path = out / "family_candidate_bridge_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "manifest": str(manifest_path),
        "bridge_tsv": str(bridge_path),
        "route_a_review_candidates": str(review_projection_path),
        "review_candidates": str(review_projection_path),
        "review_manifest": str(review_manifest_path),
        "row_count": int(len(frame)),
        "by_route": _counts("route"),
        "rank_contribution": "none",
    }
