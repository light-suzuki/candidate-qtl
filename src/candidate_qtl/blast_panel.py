"""Deterministic blast-panel freeze layer (issue #374 / #355).

Issue #355 introduced ``route-a-interval-match``: only *verified* baits frozen
by ``bait-panel-freeze`` may be forward-matched to concrete genes inside the
QTL interval, and a blast hit is evidence of similarity only -- never causal
proof, never a gene score. Issue #374 ("recall-first Top-20 benchmark
objective") asks for the blast evidence path to be a closed, auditable layer:
the panel that becomes BLAST queries is frozen byte-for-byte *before* any
interval protein is opened, and the freeze manifest pins every row fingerprint
so a later stage (``blast_evidence.py`` / the ``blast-panel-search`` command,
written by a different owner) can verify immutability before use.

This module is the deterministic core of that split:

* ``freeze_blast_panel`` takes exactly one input -- a panel TSV -- and writes
  ``frozen_bait_panel.tsv`` plus ``blast_panel_freeze_manifest.json``. No
  interval, candidate, truth or score input is accepted structurally; a
  forbidden column token is a hard error.
* ``verify_blast_panel_freeze`` re-hashes rows and files and refuses any drift.
* ``validate_panel_row`` / ``_generic_motif`` implement the row-level verdicts
  (accession format, sequence sanity, generic-motif rejection, fingerprint,
  provenance) without touching anything outside stdlib + pandas + core. A row
  may be ``projection_eligibility=review_only`` only when its provider
  identity/species/sequence evidence is present; it stays review-only through
  projection and never becomes scientific support.

Import policy (pinned by ``tests/test_blast_panel_rank_invariance.py``):
``llm_pathway`` / ``family_index`` / ``hardening`` / ``shortlist`` /
``ranking_policy`` must never be imported here, so freezing a panel can never
mutate ranking state or pull in a model pathway.
"""

from __future__ import annotations

import hashlib
import json
import re
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd

from . import __version__
from .core import file_fingerprint, runtime_path, safe_id, sha256
from .identifiers import is_uniprot_accession

COMMAND = "blast-panel-freeze"
FREEZE_VERSION = "blast-panel-freeze-1.2.0"
SCHEMA_VERSION = "1.1"

PANEL_COLUMNS = [
    "panel_id", "mechanism", "family", "synonyms", "reference_gene_name",
    "protein_accession", "taxon", "source", "source_url", "retrieved_at",
    "verification_status", "verification_accession", "verification_sequence_sha256",
    "verification_source_url", "verification_retrieved_at",
    "trait_direction", "query_sequence", "query_sequence_sha256", "provider",
    "model", "prompt_fingerprint", "review_status", "review_reason",
    "rank_contribution", "projection_eligibility",
]

REVIEW_STATUS_VOCAB = ("verified", "review_required", "rejected", "not_evaluated")

# Issue #409: a canonical stable ID carries the database it came from
# ("UniProt:Q9TES1", "TAIR:AT1G69120"). Normalizing the prefix away would drop
# the provenance that source verification and accession checks read, so the
# colon form is accepted for a fixed allowlist of stable-ID databases instead.
# Everything outside the allowlist stays fail-closed exactly as before, so a
# free-form "bait:1" is still refused.
PANEL_ID_NAMESPACES = frozenset({
    "uniprot", "tair", "refseq", "ncbi", "genbank", "embl", "ddbj",
    "ensembl", "ensemblplants", "phytozome", "gramene", "medtr", "medicago",
    "lotus", "panther", "interpro", "pfam",
})
_PANEL_ID_BARE = re.compile(r"^[A-Za-z0-9_.-]+$")
# Two characters minimum, leading alphanumeric: "TAIR:" (empty local part) and
# "TAIR:-" must not pass.
_PANEL_ID_LOCAL = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]+$")


def _panel_id_is_safe(panel_id: str) -> bool:
    """Whether ``panel_id`` may be frozen as-is.

    Two accepted forms: the historical bare token ``[A-Za-z0-9_.-]+``, and
    ``namespace:local`` where the namespace is a known stable-ID database
    (case-insensitive) and the local part is itself a safe token. A second
    colon, whitespace, or an unknown namespace fails.

    The panel_id keeps its canonical colon form in the manifest keys and the
    frozen TSV. Downstream ``blast_evidence`` normalizes with ``safe_id`` on
    both sides -- FASTA id generation and hit matching -- so the query path
    cannot drift apart from the frozen record.
    """
    if _PANEL_ID_BARE.match(panel_id):
        return True
    namespace, separator, local = panel_id.partition(":")
    if not separator or namespace.casefold() not in PANEL_ID_NAMESPACES:
        return False
    return bool(_PANEL_ID_LOCAL.match(local))

# Any column name that tokenizes onto one of these is refused: an interval /
# candidate / truth / score shaped column cannot ride into a freeze.
FORBIDDEN_COLUMN_TOKENS = frozenset({"interval", "candidate", "truth", "answer", "causal", "score", "gene_id"})

TRAIT_DIRECTION_VOCAB = {"increase", "decrease", "gain", "loss", "unknown"}

# UniProt uses the shared canonical grammar above. These non-UniProt protein
# identifier forms are also accepted for local/reference-backed panels.
NON_UNIPROT_ACCESSION_PATTERNS = (
    re.compile(r"^(XP|NP|YP|WP)_[0-9]+(\.[0-9]+)?$"),
    re.compile(r"^ENSP[0-9]+(\.[0-9]+)?$"),
)

DEFAULT_THRESHOLDS: dict[str, float | int] = {
    "min_query_length": 30,
    "min_identity": 30.0,
    "min_query_coverage": 0.6,
    "min_target_coverage": 0.6,
    "max_evalue": 1e-5,
    "paralog_tie_bitscore_ratio": 0.9,
    "max_single_aa_fraction": 0.25,
    "max_homopolymer_run": 10,
    "max_unknown_aa_fraction": 0.5,
}

# 20 canonical amino acids; X is a tolerated "unknown" placeholder handled
# separately by the unknown-fraction threshold.
AA_ALPHABET = frozenset("ACDEFGHIKLMNPQRSTVWY")

OUTPUT_PANEL_NAME = "frozen_bait_panel.tsv"
OUTPUT_MANIFEST_NAME = "blast_panel_freeze_manifest.json"

TRUTH_ABSENCE_ASSERTION = (
    "no interval, candidate, truth or score input was supplied to this freeze; "
    "the panel TSV is the only input, and it is frozen before any interval "
    "protein is opened"
)

POLICY = [
    "the panel TSV is the only input; interval, candidate, truth and score inputs are structurally absent",
    "interval_accessed is false at freeze time; only blast-panel-search may open interval proteins, after verifying this manifest",
    "every row records rank_contribution=none; a panel row carrying any other value is rejected",
    "verified-support rows and identity-verified review-only rows may become BLAST queries; review-only projections remain review-only and never support",
    "frozen artifacts are immutable: verify re-hashes rows and files and refuses any drift",
    "blast hits are evidence of similarity only; they are never causal proof and never add to a gene score",
]

VERIFIED = "verified"
REVIEW_REQUIRED = "review_required"
REJECTED = "rejected"


# --------------------------------------------------------------------------- #
# row-level validation
# --------------------------------------------------------------------------- #

def _accession_format_valid(accession: str) -> bool:
    return is_uniprot_accession(accession) or any(pattern.fullmatch(accession) for pattern in NON_UNIPROT_ACCESSION_PATTERNS)


def _row_sha256(sequence: str) -> str:
    """SHA-256 of the normalized (upper-cased) amino-acid sequence bytes."""
    return hashlib.sha256(sequence.upper().encode("utf-8")).hexdigest()


def _single_aa_fraction(sequence: str) -> tuple[str, float]:
    counts: dict[str, int] = {}
    for residue in sequence:
        counts[residue] = counts.get(residue, 0) + 1
    if not counts:
        return "", 0.0
    most_common = max(counts, key=counts.get)
    return most_common, counts[most_common] / len(sequence)


def _longest_homopolymer_run(sequence: str) -> int:
    longest = 0
    current = 0
    previous = ""
    for residue in sequence:
        if residue == previous:
            current += 1
        else:
            current = 1
            previous = residue
        if current > longest:
            longest = current
    return longest


LOCAL_TANDEM_MIN_LEN = 12
LOCAL_TANDEM_MAX_PERIOD = 6


def _has_local_tandem_array(
    sequence: str,
    min_len: int = LOCAL_TANDEM_MIN_LEN,
    max_period: int = LOCAL_TANDEM_MAX_PERIOD,
) -> bool:
    """Return whether the sequence holds a perfect local tandem repeat span.

    For every start position and every period p in {2..max_period}, checks
    whether the contiguous span of n = ceil(min_len / p) * p residues equals
    the p-mer unit repeated n // p times, i.e. a tandem repeat region of at
    least min_len residues. min_len=12 is the local window length used by SEG
    (Wootton & Federhen 1993) and matches the local-segment definitions of
    pfilt / CAST; these constants are definitional, not tuned to any benchmark.
    """
    for start in range(len(sequence) - min_len + 1):
        for period in range(2, max_period + 1):
            span = ((min_len + period - 1) // period) * period
            unit = sequence[start : start + period]
            if sequence[start : start + span] == unit * (span // period):
                return True
    return False


def _generic_motif(sequence: str) -> tuple[bool, str]:
    """Report whether a sequence looks like a hand-written generic motif.

    Triggers, in order: (a) one amino acid is >= 25% of the residues,
    (b) a homopolymer run of >= 10, (c) a local tandem repeat span of >= 12
    residues (perfect repetition of a 2..6-residue unit). Returns
    ``(triggered, detail)`` where detail names the first rule that fired.
    """
    _, fraction = _single_aa_fraction(sequence)
    if fraction >= DEFAULT_THRESHOLDS["max_single_aa_fraction"]:
        return True, "single_amino_acid_fraction"
    if _longest_homopolymer_run(sequence) >= DEFAULT_THRESHOLDS["max_homopolymer_run"]:
        return True, "homopolymer_run"
    if _has_local_tandem_array(sequence):
        return True, "tandem_repeat"
    return False, ""


def _normalize_trait_direction(value: str) -> str:
    folded = value.strip().casefold()
    return folded if folded in TRAIT_DIRECTION_VOCAB else "unknown"


def validate_panel_row(row: dict[str, str]) -> dict[str, object]:
    """Run the ordered row checks and return a verdict (never raises).

    Check order: V-ACC -> V-SEQ -> V-AA -> V-LEN -> V-GENERIC -> V-SHA ->
    V-PROV. Input ``review_status`` is sticky: ``rejected``,
    ``review_required`` and ``not_evaluated`` are kept regardless of the
    mechanical checks (with the input review_reason); ``verified`` or an empty
    value lets the mechanical checks decide -- all checks passing means
    ``verified``, any failure means ``review_required(reason=...)``. Any input
    status outside ``REVIEW_STATUS_VOCAB`` is a fail-safe
    ``review_required(reason=invalid_input_review_status)``.
    """
    values = {column: str(row.get(column, "")).strip() for column in PANEL_COLUMNS}
    sequence = values["query_sequence"]
    input_status = values["review_status"]
    if input_status == REJECTED:
        values["review_status"] = REJECTED
        values["review_reason"] = values["review_reason"] or "rejected_by_input"
    elif input_status in (REVIEW_REQUIRED, "not_evaluated"):
        values["review_status"] = input_status
        values["review_reason"] = values["review_reason"]
    elif input_status in (VERIFIED, ""):
        values["review_status"] = VERIFIED
        values["review_reason"] = ""
    else:
        values["review_status"] = REVIEW_REQUIRED
        values["review_reason"] = "invalid_input_review_status"

    verdict = _validate_sequence_checks(values)
    if verdict is not None and input_status in (VERIFIED, ""):
        values["review_status"], values["review_reason"] = verdict
    values["trait_direction"] = _normalize_trait_direction(values["trait_direction"])
    eligibility = values["projection_eligibility"].casefold()
    if eligibility not in {"", "verified_support", "review_only"}:
        values["review_status"] = REVIEW_REQUIRED
        values["review_reason"] = "invalid_projection_eligibility"
    elif eligibility == "review_only" and values["review_status"] != REVIEW_REQUIRED:
        values["review_status"] = REVIEW_REQUIRED
        values["review_reason"] = values["review_reason"] or "review_only_projection"
    normalized = sequence.upper()
    values["query_sequence"] = normalized
    values["query_sequence_sha256"] = _row_sha256(normalized)
    return {
        **values,
        "query_length": len(normalized),
        "is_generic_motif": bool(_generic_motif(normalized)[0]),
    }


def _validate_sequence_checks(values: dict[str, str]) -> tuple[str, str] | None:
    """The ordered checks below accession/provenance; None means all passed."""
    accession = values["protein_accession"]
    if not _accession_format_valid(accession):
        return REVIEW_REQUIRED, "bad_accession_format"
    sequence = values["query_sequence"]
    if not sequence:
        return REVIEW_REQUIRED, "empty_sequence"
    unknown_fraction = sequence.upper().count("X") / len(sequence)
    invalid_residues = [residue for residue in sequence.upper() if residue not in AA_ALPHABET and residue != "X"]
    if invalid_residues or unknown_fraction > DEFAULT_THRESHOLDS["max_unknown_aa_fraction"]:
        return REVIEW_REQUIRED, "invalid_aa_or_high_unknown"
    if len(sequence) < DEFAULT_THRESHOLDS["min_query_length"]:
        return REVIEW_REQUIRED, "short_query"
    if _generic_motif(sequence.upper())[0]:
        return REVIEW_REQUIRED, "generic_motif"
    expected = values["query_sequence_sha256"]
    if not expected or _row_sha256(sequence.upper()) != expected.strip():
        return REVIEW_REQUIRED, "sequence_fingerprint_mismatch"
    verification_fields = (
        "verification_status", "verification_accession", "verification_sequence_sha256",
        "verification_source_url", "verification_retrieved_at",
    )
    if not values["source"] or not values["source_url"] or not values["retrieved_at"]:
        return REVIEW_REQUIRED, "missing_provenance"
    if any(not values[field] for field in verification_fields):
        return REVIEW_REQUIRED, "missing_verification_evidence"
    if values["verification_status"] != VERIFIED:
        return REVIEW_REQUIRED, "verification_not_verified"
    if values["verification_accession"] != accession:
        return REVIEW_REQUIRED, "verification_accession_mismatch"
    if values["verification_sequence_sha256"] != _row_sha256(sequence.upper()):
        return REVIEW_REQUIRED, "verification_sequence_fingerprint_mismatch"
    return None


# --------------------------------------------------------------------------- #
# freeze
# --------------------------------------------------------------------------- #

def _read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False).fillna("")


def _frame(rows: list[dict[str, str]], columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=columns)


def _write(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n", encoding="utf-8")


def _forbidden_columns(frame: pd.DataFrame) -> list[str]:
    """Refuse columns whose normalized form contains a forbidden token.

    A column is lower-cased and stripped of every non-alphanumeric character
    (``answer-gene_id`` -> ``answergeneid``), and each FORBIDDEN_COLUMN_TOKENS
    entry is normalized the same way; any substring match rejects the column.
    ``reference_gene_name`` (-> ``referencegenename``) and
    ``rank_contribution`` (-> ``rankcontribution``) cannot be false positives.
    """
    normalized_tokens = frozenset(re.sub(r"[^a-z0-9]+", "", token.casefold()) for token in FORBIDDEN_COLUMN_TOKENS)
    forbidden: list[str] = []
    for column in frame.columns:
        normalized = re.sub(r"[^a-z0-9]+", "", str(column).casefold())
        if any(token in normalized for token in normalized_tokens):
            forbidden.append(str(column))
    return forbidden


def _now() -> str:
    return datetime.now(UTC).isoformat(timespec="seconds")


def freeze_blast_panel(
    panel_arg: str,
    out_dir_arg: str,
    case_id: str = "",
    leak_check_candidates_arg: str | None = None,
    min_panel_size: int = 2,
) -> dict[str, object]:
    """Freeze a panel TSV into an immutable bait panel plus a manifest.

    Only the panel TSV is read. Column / duplicate / rank_contribution rules
    fail closed with ValueError; every row is validated and normalized
    (sequence upper-cased, SHA-256 recomputed); the freeze always completes,
    reporting ``panel_status=review_required_small_panel`` when fewer than
    ``min_panel_size`` rows are verified.
    """
    panel_path = runtime_path(panel_arg)
    out_dir = runtime_path(out_dir_arg)
    if leak_check_candidates_arg is not None:
        raise ValueError(
            "leak_check_candidates_arg is no longer accepted: candidate input before the "
            "freeze would violate the pre-interval firewall"
        )
    if not panel_path.exists():
        raise ValueError(f"panel TSV does not exist: {panel_path}")
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        frame = _read_tsv(panel_path)
    except (pd.errors.EmptyDataError, pd.errors.ParserError) as error:
        raise ValueError(f"panel TSV is not a readable table: {panel_path}") from error

    # ``projection_eligibility`` was added after the first frozen-panel
    # contract.  Legacy panels default to canonical support, while new panels
    # state the distinction explicitly.
    missing = [column for column in PANEL_COLUMNS
               if column not in frame.columns and column != "projection_eligibility"]
    if missing:
        raise ValueError(f"panel TSV is missing required columns: {', '.join(sorted(missing))}")

    forbidden = _forbidden_columns(frame)
    if forbidden:
        raise ValueError(f"panel TSV contains forbidden columns: {', '.join(sorted(forbidden))}")

    panel_ids = [str(value).strip() for value in frame["panel_id"].fillna("")]
    for panel_id in panel_ids:
        if not panel_id:
            raise ValueError("panel TSV has an empty panel_id value")
        if not _panel_id_is_safe(panel_id):
            raise ValueError(
                f"panel TSV has an unsafe panel_id value: {panel_id!r} "
                "(allowed: [A-Za-z0-9_.-]+, or <namespace>:<local> for a known stable-ID database "
                f"[{', '.join(sorted(PANEL_ID_NAMESPACES))}])"
            )
    duplicate_ids = sorted({panel_id for panel_id in panel_ids if panel_ids.count(panel_id) > 1})
    if duplicate_ids:
        raise ValueError(f"panel TSV has duplicate panel_id values: {', '.join(duplicate_ids)}")
    # Issue #409: admitting the colon form creates a collision class the bare
    # form did not have -- "UniProt:Q9TES1" and "UniProt_Q9TES1" are distinct
    # panel_ids that become one FASTA id downstream (blast_evidence normalizes
    # with safe_id on both the write and the lookup side). Refuse them here
    # rather than let two baits silently share one query record.
    normalized_ids = [safe_id(panel_id) for panel_id in panel_ids]
    collisions = sorted({normalized for normalized in normalized_ids if normalized_ids.count(normalized) > 1})
    if collisions:
        raise ValueError(
            "panel TSV has panel_id values that collide once normalized for the query FASTA: "
            f"{', '.join(collisions)}"
        )

    if not (frame["rank_contribution"].fillna("").astype(str).eq("none")).all():
        raise ValueError("panel TSV must record rank_contribution=none on every row")

    rows: list[dict[str, str]] = []
    row_verdicts: dict[str, dict[str, object]] = {}
    for record in frame.to_dict("records"):
        values = {column: str(record.get(column, "")).strip() for column in PANEL_COLUMNS}
        values["projection_eligibility"] = values["projection_eligibility"] or "verified_support"
        verdict = validate_panel_row(values)
        row_verdicts[str(verdict["panel_id"])] = verdict
        rows.append({column: str(verdict[column]) for column in PANEL_COLUMNS})

    validated = _frame(rows, PANEL_COLUMNS).sort_values("panel_id", kind="stable")
    panel_out = out_dir / OUTPUT_PANEL_NAME
    _write(validated, panel_out)

    verified_count = int(sum(1 for verdict in row_verdicts.values() if verdict["review_status"] == VERIFIED))
    review_required_count = int(sum(1 for verdict in row_verdicts.values() if verdict["review_status"] != VERIFIED))
    panel_status = VERIFIED if verified_count >= min_panel_size else "review_required_small_panel"

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "freeze_version": FREEZE_VERSION,
        "command": COMMAND,
        "candidate_qtl_version": __version__,
        "case_id": case_id,
        "frozen_at_utc": _now(),
        "panel_status": panel_status,
        "row_count": int(len(validated)),
        "verified_count": verified_count,
        "review_required_count": review_required_count,
        "min_panel_size": int(min_panel_size),
        "interval_accessed": False,
        "truth_inputs_present": False,
        "truth_absence_assertion": TRUTH_ABSENCE_ASSERTION,
        "input_allowlist": ["panel"],
        "inputs": {
            "panel": file_fingerprint(panel_path, full_hash=True),
            OUTPUT_PANEL_NAME: file_fingerprint(panel_out, full_hash=True),
        },
        "panel_sha256": sha256(panel_out),
        "rows": {
            str(verdict["panel_id"]): {
                "query_sequence_sha256": str(verdict["query_sequence_sha256"]),
                "query_length": int(verdict["query_length"]),
                "verification_status": str(verdict["verification_status"]),
                "verification_accession": str(verdict["verification_accession"]),
                "verification_sequence_sha256": str(verdict["verification_sequence_sha256"]),
                "review_status": str(verdict["review_status"]),
                "review_reason": str(verdict["review_reason"]),
                "rank_contribution": str(verdict["rank_contribution"]),
            }
            for verdict in row_verdicts.values()
        },
        "thresholds": dict(DEFAULT_THRESHOLDS),
        "policy": list(POLICY),
        "outputs": {
            "panel": str(panel_out),
            "freeze_manifest": str(out_dir / OUTPUT_MANIFEST_NAME),
        },
    }
    manifest_path = out_dir / OUTPUT_MANIFEST_NAME
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    return {
        "manifest": str(manifest_path),
        "panel": str(panel_out),
        "row_count": int(len(validated)),
        "verified_count": verified_count,
        "review_required_count": review_required_count,
        "panel_status": panel_status,
        "interval_accessed": False,
        "truth_inputs_present": False,
        "frozen_at_utc": str(manifest["frozen_at_utc"]),
    }


# --------------------------------------------------------------------------- #
# convert the bait-panel-freeze output (issue #374 / Criterion 3)
# --------------------------------------------------------------------------- #

def panel_from_verified_bait_panel(panel_arg: str, sequences_arg: str) -> pd.DataFrame:
    """Convert a ``bait-panel-freeze`` (llm_pathway) verified bait panel into
    the blast-panel freeze input schema.

    The guided triple-route profile freezes the Route A panel with
    ``bait-panel-freeze`` and re-freezes it here so ``blast-panel-search`` can
    verify it (fail closed) before opening the interval proteins. The upstream
    per-row ``review_status`` is inherited, never upgraded: only rows the bait
    freeze already marked ``verified`` convert to verified blast-panel rows; a
    row the upstream freeze marked ``review_required`` for a bait-level reason
    converts to a ``review_required`` blast-panel row, so the search skips the
    row and neither is ever a failure or a silent promotion. Panel-level
    diversity review is carried separately by the freeze manifest (#538).

    The query/verification SHA-256 fingerprints are recomputed with the
    blast-panel normalization (upper-cased bytes) so every converted row passes
    ``validate_panel_row``; the freeze is the single source of truth after that.
    """
    panel = _read_tsv(runtime_path(panel_arg))
    missing = sorted({"hypothesis_id", "bait_id", "sequence_sha256"} - set(panel.columns))
    if missing:
        raise ValueError(f"verified bait panel is missing columns: {', '.join(missing)}")
    sequences = json.loads(runtime_path(sequences_arg).read_text(encoding="utf-8"))

    rows: list[dict[str, str]] = []
    seen_baits: dict[str, tuple[str, str]] = {}
    for record in panel.sort_values("hypothesis_id", kind="stable").to_dict("records"):
        hypothesis_id = str(record.get("hypothesis_id", "")).strip()
        bait_id = str(record.get("bait_id", "")).strip()
        sequence = str(sequences.get(hypothesis_id, "")).upper()
        if not bait_id or not sequence:
            continue
        bait_key = bait_id.casefold()
        previous = seen_baits.get(bait_key)
        if previous is not None:
            previous_sequence, previous_hypothesis_id = previous
            if previous_sequence != sequence:
                raise ValueError(
                    "verified bait panel maps one bait_id to different sequences: "
                    f"{bait_id} ({previous_hypothesis_id}, {hypothesis_id})"
                )
            # Upstream intentionally preserves hypothesis-level provenance, so
            # several hypotheses may resolve to the same accession. BLAST panel
            # IDs are accession-level and must be unique; keep the first row in
            # stable hypothesis order when the verified sequence is identical.
            continue
        seen_baits[bait_key] = (sequence, hypothesis_id)
        # The upstream bait-panel freeze derives a per-row review_status from
        # the panel status (verified for a frozen panel, review_required for a
        # small one). Inherit it: a non-verified row is never promoted to
        # VERIFIED by this converter, so the downstream freeze and search keep
        # the fail-closed policy.
        upstream_status = str(record.get("review_status", "")).strip().casefold()
        if upstream_status == VERIFIED:
            row_review_status = VERIFIED
            row_review_reason = ""
        else:
            row_review_status = REVIEW_REQUIRED
            row_review_reason = f"upstream_bait_review_status:{upstream_status or 'missing'}"
        source = str(record.get("source_provider", "")).strip()
        source_url = str(record.get("source_url", "")).strip()
        retrieved_at = str(record.get("retrieved_at", "")).strip()
        digest = _row_sha256(sequence)
        rows.append({
            "panel_id": bait_id,
            "mechanism": str(record.get("pathway_or_process", "")).strip(),
            "family": str(record.get("gene_family", "")).strip(),
            "synonyms": "",
            "reference_gene_name": str(record.get("reference_gene_name", "")).strip(),
            "protein_accession": bait_id,
            "taxon": str(record.get("reference_species", "")).strip(),
            "source": source,
            "source_url": source_url,
            "retrieved_at": retrieved_at,
            "verification_status": VERIFIED,
            "verification_accession": bait_id,
            "verification_sequence_sha256": digest,
            "verification_source_url": source_url,
            "verification_retrieved_at": retrieved_at,
            "trait_direction": "unknown",
            "query_sequence": sequence,
            "query_sequence_sha256": digest,
            "provider": source,
            "model": "",
            "prompt_fingerprint": "",
            "review_status": row_review_status,
            "review_reason": row_review_reason,
            "rank_contribution": "none",
            "projection_eligibility": str(record.get("projection_eligibility", "verified_support")).strip() or "verified_support",
        })
    return _frame(rows, PANEL_COLUMNS)


# --------------------------------------------------------------------------- #
# verify
# --------------------------------------------------------------------------- #

def _verify_fingerprint(fingerprint: dict[str, object], actual_path: Path, item: str) -> None:
    if not actual_path.exists():
        raise ValueError(f"{item}: frozen input file is missing: {actual_path}")
    recorded = str(fingerprint.get("sha256") or "")
    if not recorded:
        raise ValueError(f"{item}: fingerprint has no sha256")
    try:
        digest = sha256(actual_path)
    except OSError as error:
        raise ValueError(f"{item}: cannot hash frozen input file: {actual_path}: {error}") from error
    if digest != recorded:
        raise ValueError(f"{item}: file content drifted from the frozen fingerprint: {actual_path}")


def _verify_thresholds(thresholds: object) -> None:
    """Validate the manifest thresholds: keys and value ranges (fail closed).

    The freeze takes no threshold input, so every key and value must exactly
    match DEFAULT_THRESHOLDS. This prevents a modified manifest from changing
    the later search contract.
    """
    if not isinstance(thresholds, dict):
        raise ValueError("thresholds: manifest thresholds must be an object")
    if set(thresholds) != set(DEFAULT_THRESHOLDS):
        raise ValueError(f"thresholds: manifest thresholds keys must match {sorted(DEFAULT_THRESHOLDS)}")

    if thresholds != DEFAULT_THRESHOLDS:
        raise ValueError("thresholds: manifest thresholds must exactly match the frozen defaults")


def verify_blast_panel_freeze(freeze_arg: str, panel_arg: str | None = None) -> dict[str, object]:
    """Verify a freeze manifest against the actual artifacts (first violation wins).

    Checks, in order: (1) manifest JSON parses, (2) command and freeze_version
    match, (3) interval_accessed is the strict bool False, (4)
    truth_inputs_present is the strict bool False, (5) thresholds match
    DEFAULT_THRESHOLDS keys and documented value ranges, (6) every
    input_allowlist member has a fingerprint in inputs, (7) every input file
    re-hashes to its fingerprint, (8) panel_sha256 matches the frozen panel
    file, (9) every manifest row fingerprint matches the frozen panel row.
    """
    freeze_path = runtime_path(freeze_arg)
    if not freeze_path.exists():
        raise ValueError(f"manifest_json: freeze manifest does not exist: {freeze_path}")
    try:
        freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"manifest_json: freeze manifest is not valid JSON: {error}") from error

    if str(freeze.get("command")) != COMMAND or str(freeze.get("freeze_version")) != FREEZE_VERSION:
        raise ValueError(f"command_and_freeze_version: expected command={COMMAND} freeze_version={FREEZE_VERSION}")

    if freeze.get("interval_accessed") is not False:
        raise ValueError("interval_accessed: must be the strict boolean false at freeze time")
    if freeze.get("truth_inputs_present") is not False:
        raise ValueError("truth_inputs_present: must be the strict boolean false at freeze time")
    _verify_thresholds(freeze.get("thresholds"))

    inputs = freeze.get("inputs", {})
    if not isinstance(inputs, dict):
        raise ValueError("inputs: manifest inputs must be an object")
    allowlist = freeze.get("input_allowlist", [])
    if allowlist != ["panel"] or set(inputs) != {"panel", OUTPUT_PANEL_NAME}:
        raise ValueError("input_allowlist: only the source panel and frozen panel fingerprints are allowed")

    for name, fingerprint in inputs.items():
        if not isinstance(fingerprint, dict):
            raise ValueError(f"inputs_fingerprints: input {name} has no fingerprint object")
        if name == OUTPUT_PANEL_NAME and panel_arg:
            path = runtime_path(str(panel_arg))
        else:
            recorded_path = str(fingerprint.get("path") or "")
            if not recorded_path:
                raise ValueError(f"inputs_fingerprints: input {name} has no recorded path")
            path = runtime_path(recorded_path)
        _verify_fingerprint(fingerprint, path, "inputs_fingerprints")

    panel_path = runtime_path(str(panel_arg)) if panel_arg else out_panel_from_freeze(freeze, freeze_path)
    if not panel_path.exists():
        raise ValueError(f"panel_sha256: frozen panel file is missing: {panel_path}")
    if sha256(panel_path) != str(freeze.get("panel_sha256") or ""):
        raise ValueError(f"panel_sha256: frozen panel file drifted: {panel_path}")

    panel_rows = _read_tsv(panel_path)
    panel_by_id = {str(record["panel_id"]): record for record in panel_rows.to_dict("records")}
    rows = freeze.get("rows", {})
    if not isinstance(rows, dict):
        raise ValueError("row_sequence_fingerprints: manifest rows must be an object")
    for panel_id, record in rows.items():
        if panel_id not in panel_by_id:
            raise ValueError(f"row_sequence_fingerprints: panel_id missing from frozen panel: {panel_id}")
        expected = str(record.get("query_sequence_sha256") or "")
        actual = str(panel_by_id[panel_id]["query_sequence_sha256"])
        if expected != actual:
            raise ValueError(f"row_sequence_fingerprints: row fingerprint drifted for {panel_id}")

    return {
        "verified": True,
        "freeze": str(freeze_path),
        "panel": str(panel_path),
        "checks": {
            "manifest_json": "pass",
            "command_and_freeze_version": "pass",
            "interval_accessed": "pass",
            "truth_inputs_present": "pass",
            "thresholds": "pass",
            "input_allowlist": "pass",
            "inputs_fingerprints": "pass",
            "panel_sha256": "pass",
            "row_sequence_fingerprints": "pass",
        },
        "interval_accessed": False,
    }


def out_panel_from_freeze(freeze: dict[str, object], freeze_path: Path) -> Path:
    """The frozen panel path recorded in the manifest, with a fallback next to the manifest."""
    outputs = freeze.get("outputs", {})
    if isinstance(outputs, dict):
        recorded = str(outputs.get("panel") or "")
        if recorded:
            return runtime_path(recorded)
    return freeze_path.parent / OUTPUT_PANEL_NAME
