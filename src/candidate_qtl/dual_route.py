"""Forward-lookup and reverse-lookup discovery routes (issue #169).

Two independent ways to reach a pathway hypothesis, kept separate until the end:

* **Route A - phenotype-first (順引き).** Pin down the gene or family from the
  phenotype and known biology, then look for it in the interval. Strong when the
  trait or its pathway is well studied in a model or relative species. Implemented
  by the existing ``llm-pathway-*`` / ``hypothesis-verify`` / ``bait-panel-freeze``
  / ``pathway-interval-falsify`` chain; this module only reads their artifacts and
  detects the failure signals that should trigger Route B.

* **Route B - interval-first reverse lookup (逆引き).** Start from the genes that
  are actually in the interval and ask what they are known to be. The functional
  description may come from any annotated reference: a model-species homolog, a
  relative, or the focal reference genome's own annotation. Hypotheses come back
  from the interval rather than from a guess made before looking.

The two routes fail in **mirror-image** ways, which is the reason for running both:

* Route A can hallucinate a gene that is not there. When the pinned-down family
  turns out to be absent from the interval, the natural drift is toward ever more
  minor candidates, and a model asked to keep answering will invent one. The
  guard is anchoring: a hypothesis with no interval member is marked
  ``unanchored``, never quietly carried.
* Route B can rationalise a gene that *is* there. Everything it proposes exists,
  but a plausible post-hoc story can be told about almost any gene once you are
  looking at it. The guard is that Route B never sees the phenotype-first answer
  and is scored in its own lane.

The two are deliberately **not** merged while they are being generated. Mixing
them would let a hypothesis formed after seeing the interval be scored as
phenotype-only discovery. They are compared only at the end, which keeps
``mechanistic_convergence`` / ``family_or_subfamily_convergence`` /
``gene_overlap_only`` / ``process_related`` / ``text_only_match`` /
``phenotype_only_supported`` / ``interval_only_supported`` / ``discordant`` /
``both_weak`` apart and reports each route's own risk. A shared interval gene
alone is ``gene_overlap_only``: it never upgrades to a shared mechanism by
itself (issue #186).

Boundaries:

* Route B is **recall-first**. An adjacent family or process (PXY/TDR when the
  first guess was CLV-WUS) is kept, but not at the same strength: every hypothesis
  carries a ``relationship_tier`` rather than a single merged score.
* Family or process membership alone never changes a gene score. Everything this
  module writes carries ``rank_contribution=none``.
* Unknown / orphan / dark interval genes are carried into Route B as their own
  probe instead of being dropped for lacking annotation.
* Route B may legitimately see the interval, so its artifacts are labelled
  ``interval_informed``; they must never be counted as blind phenotype-only
  discovery performance.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

from .core import file_fingerprint, runtime_path, safe_id
from .family_labels import alias_canonical, family_identity_keys
from .hardening import classify_source
from .identifiers import is_usable_symbol
from .phenotype_context import structured_context_block

SCHEMA_VERSION = "1.4"
ROUTE_A = "phenotype_first"
ROUTE_B = "interval_first"

FAMILY_TO_COPY_INPUT_COLUMNS = ("hypothesis_id", "family_id", "gene_family", "member_interval_genes")

# Ordered strongest -> weakest. Kept as a tier, never collapsed into one number.
RELATIONSHIP_TIERS = (
    "exact_ortholog", "same_subfamily", "same_family",
    "same_signaling_architecture", "same_broad_process", "unclassified",
)
# Per-pair convergence decomposition. "discordant" / "both_weak"
# never appear on a pair row: they are overall-level states kept here so the
# full vocabulary is documented in one place.
CONVERGENCE_STATES = (
    "mechanistic_convergence", "family_or_subfamily_convergence", "gene_overlap_only",
    "process_related", "text_only_match", "phenotype_only_supported",
    "interval_only_supported", "discordant", "both_weak",
)
# Overall strictness, strongest -> weakest.
OVERALL_STATES = (
    "strict_mechanistic_convergent", "family_convergent_only", "gene_overlap_only",
    "process_related", "text_only_match", "discordant", "single_route_supported",
    "both_weak",
)
MATCH_TIERS = (
    "tier1_mechanistic", "tier2_family", "tier3_gene_overlap",
    "tier4_process_related", "tier5_text_only",
)
SOURCE_INDEPENDENCE_STATUSES = (
    "shared_source_chain", "independent_source_chains", "source_ids_incomplete",
    "not_applicable",
)
ROUTE_A_SIGNALS = ("no_supported_pathway", "no_verified_bait",
                   "expected_family_absent_from_interval", "pathway_supported_but_interval_mismatch")

# Curated process vocabulary: canonical concept id -> (normalized aliases, parent
# concept id, known stable database IDs). Free-text labels that normalize to the
# same concept id are treated as the same concept only through this registry, and
# a stable database ID listed here resolves to the same concept as its aliases
# (e.g. "GO:0001234" and "CLV-WUS pathway" both mean process:clv-wus). Entries are
# exemplars to be extended per species / database; embedding similarity alone
# never upgrades a match.
PROCESS_VOCABULARY: dict[str, tuple[frozenset[str], str, frozenset[str]]] = {
    "process:clv-wus": (frozenset({
        "clavata-wuschel", "clv-wus pathway", "clv-wus", "clv wus",
        "stem cell maintenance via clv signaling",
        "stem-cell maintenance via clv signaling",
    }), "process:receptor-kinase-signaling", frozenset({"go:0001234"})),
    "process:cle-receptor-vascular": (frozenset({
        "cle-receptor vascular signalling", "cle-receptor vascular signaling",
        "cle receptor vascular signalling", "cle receptor vascular signaling",
        "vascular patterning via cle-receptor signaling",
        "vascular patterning via cle-receptor signalling",
    }), "process:receptor-kinase-signaling", frozenset({"go:0005678"})),
}

# Stable concept IDs embedded in free text (GO/PO/TO/pathway databases/UniProt/PMID).
_STABLE_ID_PATTERN = re.compile(
    r"(?:GO|PO|TO|EC|ko|path|map|UniProt|KEGG|MetaCyc|PlantCyc|AraCyc|SMPDB|PMID):[A-Za-z0-9._-]+"
    r"|R-[A-Za-z]{2,}-[A-Za-z0-9._-]+",
    re.IGNORECASE,
)

HYPOTHESIS_COLUMNS = [
    "hypothesis_id", "route", "scope", "pathway_or_process", "process_id", "gene_family", "reference_genes",
    "member_interval_genes", "member_count", "anchor_status", "relationship_tier", "source_backing",
    "proposal_source", "counterevidence", "verification_status", "phenotype_overlay_status",
    "phenotype_overlap_terms", "shared_route_overlay_status", "shared_route_overlap_terms",
    "rank_contribution", "notes",
]
LOOKUP_COLUMNS = [
    "gene_id", "family_id", "family_label", "reference_gene", "annotation_origin", "orthology_relationship",
    "orthogroup", "domain_architecture", "broad_function_bins", "relationship_tier", "reference_function",
    "reference_pathway", "reference_pathway_id", "reference_mutant_phenotype", "reference_source_id",
    "reference_source_url", "source_url", "literature_source_ids", "phenotype_source_ids", "verified_source_ids",
    "evidence_codes", "evidence_type", "provenance_tier", "provenance_detail", "source_database",
    "source_release", "retrieved_at", "corpus_version", "tissue_specificity", "developmental_stage",
    "protein_family", "go_terms", "keywords", "reference_taxon_id", "reference_gene_symbol",
    "reference_family_id", "reference_orthogroup_id",
    "reference_species", "protein_name", "reference_domain_ids", "reference_domain_description",
    "mapping_state", "mapping_relationship", "function_transfer_scope",
    "lookup_status", "phenotype_overlay_status", "phenotype_overlap_terms",
    "shared_route_overlay_status", "shared_route_overlap_terms", "notes",
]
OVERLAY_VOCABULARY_COLUMNS = ["vocabulary_source", "term", "matched_phenotype_token"]
COMPARISON_COLUMNS = [
    "comparison_id", "route_a_hypothesis_id", "route_b_hypothesis_id", "pathway_or_process",
    "gene_family", "shared_interval_genes", "match_basis", "match_tier", "relationship_tier",
    "route_a_source_backing", "route_b_source_backing", "convergence_status",
    "route_a_source_ids", "route_b_source_ids", "shared_source_ids",
    "source_independence_status", "interval_overlap_status", "route_a_anchor_status",
    "rank_contribution", "notes",
]
REVIEW_COLUMNS = ["reason", "target", "required_next_action", "allowed_choices"]

REVERSE_SYSTEM_PROMPT = """You are shown the annotated identity of every gene that is actually present in a
QTL interval. Propose which biological processes or mechanisms those genes could implicate.

You are NOT selecting a causal gene and NOT ranking the interval genes.

Rules:
- Work ONLY from the entries you are shown. Every proposal must name at least one
  entry from the list. Do not introduce a gene, family or accession that is absent
  from the list: an unanchored proposal is rejected downstream, so guessing costs
  you and gains nothing.
- If nothing in the list plausibly explains the phenotype, say so. "No supported
  hypothesis from this interval" is a valid and useful answer.
- Propose several competing processes, not one.
- For each proposal, name the gene family it rests on, the mechanistic link, at
  least one counterevidence or alternative explanation, and your uncertainty.
- Adjacent mechanisms are welcome (for example a different receptor-ligand pair in
  the same signalling architecture); say so explicitly rather than overstating the match.
- Every proposal is re-verified against an authoritative database afterwards."""


def _frame(rows: list[dict[str, object]], columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame(rows).reindex(columns=columns) if rows else pd.DataFrame(columns=columns)


def _read(path_arg: str | None, required: set[str], label: str) -> pd.DataFrame:
    if not path_arg:
        return pd.DataFrame(columns=sorted(required))
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{label} table is missing columns: {', '.join(missing)}")
    return frame


def _genes(value: object) -> set[str]:
    return {item.strip() for item in str(value).replace(",", ";").split(";") if item.strip()}


def _function_lookup(path_arg: str | None) -> tuple[dict[str, dict[str, object]], dict[str, dict[str, object]]]:
    """Functional descriptions keyed by reference gene and/or by focal gene.

    The reference does not have to be Arabidopsis, or any other species: a table
    keyed on ``gene_id`` lets the focal reference genome's own annotation drive the
    reverse lookup when no cross-species homolog is available.
    """
    if not path_arg:
        return {}, {}
    frame = pd.read_csv(runtime_path(path_arg), sep="\t", dtype=str).fillna("")
    if not {"reference_gene_id", "gene_id"} & set(frame.columns):
        raise ValueError("reference function table requires reference_gene_id and/or gene_id")
    def merge(rows: list[dict[str, object]]) -> dict[str, object]:
        """Merge duplicate corpus records without discarding provenance.

        A corpus can contain one row for the UniProt record and additional rows
        for phenotype/literature evidence.  Selecting the last row silently
        lost those source rosters before B0.  Text fields use the first stable
        value, while semicolon-delimited evidence fields are unioned.
        """
        merged: dict[str, object] = {}
        roster_fields = {
            "source_id", "literature_source_ids", "phenotype_source_ids", "verified_source_ids",
            "evidence_codes",
        }
        for row in rows:
            for key, value in row.items():
                text = str(value or "").strip()
                if not text:
                    continue
                if key in roster_fields:
                    values = {part.strip() for part in re.split(r"[;,|]+", str(merged.get(key, ""))) if part.strip()}
                    values.update(part.strip() for part in re.split(r"[;,|]+", text) if part.strip())
                    merged[key] = ";".join(sorted(values, key=lambda item: (item.casefold(), item)))
                elif not str(merged.get(key, "")).strip():
                    merged[key] = text
        return merged

    reference_rows: dict[str, list[dict[str, object]]] = {}
    focal_rows: dict[str, list[dict[str, object]]] = {}
    for row in frame.to_dict("records"):
        reference = str(row.get("reference_gene_id", "")).strip()
        focal = str(row.get("gene_id", "")).strip()
        if reference:
            reference_rows.setdefault(reference.casefold(), []).append(row)
        if focal:
            focal_rows.setdefault(focal.casefold(), []).append(row)
    by_reference = {key: merge(rows) for key, rows in reference_rows.items()}
    by_focal: dict[str, dict[str, object]] = {}
    for key, rows in focal_rows.items():
        references = {
            str(row.get("reference_gene_id", "")).strip().casefold()
            for row in rows if str(row.get("reference_gene_id", "")).strip()
        }
        # Never attach the function text from one reference to the unioned
        # source roster of another.  When the atlas supplies an explicit
        # reference, by_reference above remains authoritative; this focal-only
        # fallback is allowed only when all rows describe one reference.
        if len(references) <= 1:
            by_focal[key] = merge(rows)
    return by_reference, by_focal


def _tokens_from_reference(record: dict[str, object], fields: tuple[str, ...]) -> list[str]:
    """Return one stable, case-insensitive source roster from corpus fields."""
    found: dict[str, str] = {}
    for field in fields:
        for value in re.split(r"[;,|]+", str(record.get(field, "") or "")):
            token = value.strip()
            if token:
                found.setdefault(token.casefold(), token)
    return sorted(found.values(), key=lambda item: (item.casefold(), item))


# --------------------------------------------------------------------------- #
# Route A: phenotype-first status
# --------------------------------------------------------------------------- #

def phenotype_first_status(
    out_dir_arg: str,
    hypotheses_arg: str | None = None,
    validation_arg: str | None = None,
    panel_arg: str | None = None,
    falsification_arg: str | None = None,
    case_id: str = "",
    phenotype_context: str | None = None,
    phenotype_text: str = "",
) -> dict[str, object]:
    """Normalize Route A artifacts and detect the signals that should trigger Route B.

    ``phenotype_context`` (optional ``phenotype_context.yaml``) is recorded in
    the status for audit only; Route A hypothesis normalization reads the
    mechanism/verification tables, never the structured context fields.
    """
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    hypotheses = _read(hypotheses_arg, {"hypothesis_id"}, "phenotype-first hypothesis")
    validation = _read(validation_arg, {"hypothesis_id", "verification_status"}, "proposal validation")
    panel = _read(panel_arg, {"bait_id"}, "verified bait panel")
    falsification = _read(falsification_arg, {"hypothesis_id", "verdict"}, "interval falsification")

    verified = set(validation[validation.verification_status.eq("verified")]["hypothesis_id"]) if not validation.empty else set()
    # The fine level is the resolved one when it was computed; fall back to coarse.
    # A computed fine verdict may still have no resolved member (for example
    # absent_annotation_only). In that case the coarse interval member is retained:
    # fine refines copy resolution, it must not erase evidence that the family is
    # physically present in the interval.
    verdicts: dict[str, str] = {}
    interval_genes: dict[str, str] = {}
    if not falsification.empty:
        for level in ("coarse", "fine"):
            scope = falsification[falsification.get("level", pd.Series("", index=falsification.index)).eq(level)] if "level" in falsification.columns else falsification
            for row in scope.to_dict("records"):
                verdict = str(row["verdict"])
                if verdict.startswith("not_computed_"):
                    continue
                hypothesis_id = str(row["hypothesis_id"])
                verdicts[hypothesis_id] = verdict
                members = str(row.get("interval_gene_ids", ""))
                if members or not interval_genes.get(hypothesis_id):
                    interval_genes[hypothesis_id] = members

    rows: list[dict[str, object]] = []
    # Keep source-only hypothesis fields available to downstream Route A
    # consumers.  The canonical columns below are still normalized exactly as
    # before; fields such as mechanism_statement, expected_role, and a stable
    # reference identifier are evidence, not alternate normalized columns.
    source_columns = [
        str(column) for column in hypotheses.columns
        if str(column) not in HYPOTHESIS_COLUMNS
    ]
    for record in hypotheses.to_dict("records"):
        hypothesis_id = str(record["hypothesis_id"])
        verdict = verdicts.get(hypothesis_id, "not_evaluated")
        members = interval_genes.get(hypothesis_id, "")
        rows.append({
            **record,
            "hypothesis_id": hypothesis_id, "route": ROUTE_A, "scope": "pathway",
            "pathway_or_process": str(record.get("pathway_or_process", "")),
            "process_id": str(record.get("process_id", "") or record.get("pathway_id", "")),
            "gene_family": str(record.get("gene_family", "")),
            "reference_genes": str(record.get("reference_gene_name", "") or record.get("proposed_stable_id", "")),
            "member_interval_genes": members, "member_count": len(_genes(members)),
            # A phenotype-first hypothesis with no interval member is where the drift
            # toward ever more minor candidates - and invention - starts.
            "anchor_status": "interval_anchored" if _genes(members) else "unanchored",
            "relationship_tier": "exact_ortholog" if hypothesis_id in verified else "same_family",
            "source_backing": "verified_bait" if hypothesis_id in verified else "",
            "proposal_source": str(record.get("proposal_source", "")) or "phenotype_first",
            "counterevidence": str(record.get("counterevidence", "")),
            "verification_status": "verified" if hypothesis_id in verified else str(record.get("review_status", "pending_verification")),
            "rank_contribution": "none",
            "notes": f"interval verdict: {verdict}",
        })
    normalized = _frame(rows, [*HYPOTHESIS_COLUMNS, *source_columns])
    normalized.to_csv(out / "phenotype_first_hypotheses.tsv", sep="\t", index=False, lineterminator="\n")
    panel.to_csv(out / "phenotype_first_verified_baits.tsv", sep="\t", index=False, lineterminator="\n")

    falsified = sorted(key for key, value in verdicts.items() if value.startswith("falsified"))
    surviving = sorted(key for key, value in verdicts.items() if value.startswith("surviving"))
    signals: list[dict[str, object]] = []
    if hypotheses.empty:
        signals.append({"signal": "no_supported_pathway", "detail": "no phenotype-first hypothesis was produced"})
    elif not any(str(row["source_backing"]) for row in rows) and not len(panel):
        signals.append({"signal": "no_supported_pathway", "detail": "no hypothesis reached source-backed status"})
    if not len(panel) and not verified:
        signals.append({"signal": "no_verified_bait", "detail": "no proposal survived accession/species/sequence verification"})
    if falsified:
        signals.append({"signal": "expected_family_absent_from_interval",
                        "detail": f"{len(falsified)} hypothesis family(ies) have no member in the interval: {', '.join(falsified[:10])}"})
    if (len(panel) or verified) and verdicts and not surviving:
        signals.append({"signal": "pathway_supported_but_interval_mismatch",
                        "detail": "verified baits exist but no hypothesis family is present in the interval"})

    unanchored = [str(row["hypothesis_id"]) for row in rows if row["anchor_status"] == "unanchored"]
    anchored_rate = round(1 - len(unanchored) / len(rows), 4) if rows else 0.0
    status = {
        "schema_version": SCHEMA_VERSION, "route": ROUTE_A, "case_id": case_id,
        "hypothesis_count": int(len(hypotheses)), "verified_bait_count": int(len(panel) or len(verified)),
        "surviving_hypotheses": surviving, "falsified_hypotheses": falsified,
        "unanchored_hypotheses": unanchored, "interval_anchored_rate": anchored_rate,
        "risk": {
            "failure_mode": "hallucination: when the pinned-down family is absent from the interval, the route drifts toward ever more minor candidates and can invent one",
            "unanchored_hypothesis_count": len(unanchored),
            "hallucination_risk": ("high" if rows and anchored_rate < 0.5 else ("medium" if unanchored else "low")),
            "mitigation": "run interval-reverse-lookup, which can only propose what is present in the interval",
        },
        "signals": signals, "signal_names": sorted({str(item["signal"]) for item in signals}),
        "route_b_recommended": bool(signals) or bool(unanchored),
        "route_b_reason": ("Route A produced a failure signal; run interval-reverse-lookup as an independent route"
                           if signals else "Route A is supported; Route B may still be run independently"),
        "independence": "Route A artifacts stay separate from Route B until dual-route-compare",
        "structured_context": structured_context_block(
            phenotype_context, case_id=case_id, phenotype_text=phenotype_text,
            query_usage={
                "usage": "stored_only",
                "reason": ("Route A hypothesis normalization reads the mechanism/verification/panel tables "
                           "only; no structured field enters it"),
            },
        ),
        "inputs": {name: (file_fingerprint(runtime_path(value)) if value else None) for name, value in (
            ("hypotheses", hypotheses_arg), ("validation", validation_arg),
            ("panel", panel_arg), ("falsification", falsification_arg))},
    }
    status_path = out / "phenotype_first_status.json"
    status_path.write_text(json.dumps(status, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"status": str(status_path), "hypotheses": str(out / "phenotype_first_hypotheses.tsv"),
            "verified_baits": str(out / "phenotype_first_verified_baits.tsv"),
            "hypothesis_count": int(len(hypotheses)), "signal_names": status["signal_names"],
            "route_b_recommended": status["route_b_recommended"]}


# --------------------------------------------------------------------------- #
# Route B: interval-first reverse lookup
# --------------------------------------------------------------------------- #

def _tier(row: dict[str, object]) -> str:
    relationship = str(row.get("orthology_relationship", "")).casefold()
    if "one-to-one" in relationship or "one_to_one" in relationship:
        return "exact_ortholog"
    if str(row.get("orthogroup", "")):
        return "same_subfamily"
    if str(row.get("family_id", "")):
        return "same_family"
    if str(row.get("domain_architecture", "")):
        return "same_signaling_architecture"
    if str(row.get("broad_function_bins", "")):
        return "same_broad_process"
    return "unclassified"


def _best_tier(tiers: list[str]) -> str:
    for tier in RELATIONSHIP_TIERS:
        if tier in tiers:
            return tier
    return "unclassified"


_PHENOTYPE_STOPWORDS = frozenset({
    "about", "after", "against", "among", "before", "between", "causes", "drives",
    "from", "gene", "genes", "loss", "function", "plant", "plants", "protein",
    "result", "results", "that", "their", "this", "trait", "with", "without",
})


def _phenotype_terms(value: object) -> set[str]:
    return {
        token.casefold() for token in re.findall(r"[A-Za-z][A-Za-z0-9_-]{3,}", str(value or ""))
        if token.casefold() not in _PHENOTYPE_STOPWORDS
    }


def _phenotype_overlay(phenotype_terms: set[str], evidence: object) -> tuple[str, str]:
    if not phenotype_terms:
        return "not_applied", ""
    overlap = sorted(phenotype_terms & _phenotype_terms(evidence))
    return ("phenotype_consistent" if overlap else "phenotype_not_matched", ";".join(overlap))


def _shared_route_vocabulary(phenotype: str, shared_route_arg: str | None) -> tuple[set[str], list[dict[str, str]], dict[str, object]]:
    """Derive an auditable B1 overlay vocabulary from the pre-fixed shared route.

    This deliberately reuses the blind-route token and route derivation contract.
    The raw phenotype overlay remains separate: route-derived pathway/keyword
    terms are diagnostic context, never a score or rank input.
    """
    if not shared_route_arg:
        return set(), [], {
            "input_provided": False, "matched_phenotype_tokens": [],
            "pathway_families": [], "annotation_keywords": [],
            "status": "not_applied_no_shared_route",
        }
    route = _read(shared_route_arg, set(ROUTE_COLUMNS), "shared phenotype route")
    tokens = phenotype_tokens(phenotype)
    pathways, keywords = derive_route(tokens, route)
    matched_tokens = sorted({
        str(row["phenotype_token"]).strip().casefold()
        for row in route.to_dict("records")
        if str(row["phenotype_token"]).strip().casefold() in tokens
    })
    rows: list[dict[str, str]] = []
    for row in route.to_dict("records"):
        token = str(row["phenotype_token"]).strip().casefold()
        if token not in matched_tokens:
            continue
        family = str(row["pathway_family"]).strip()
        if family:
            rows.append({"vocabulary_source": "shared_route_pathway_family", "term": family,
                         "matched_phenotype_token": token})
        for keyword in re.split(r"[|;,]+", str(row["annotation_keywords"])):
            keyword = keyword.strip()
            if keyword:
                rows.append({"vocabulary_source": "shared_route_annotation_keyword", "term": keyword,
                             "matched_phenotype_token": token})
    vocabulary = _phenotype_terms(" ".join(sorted(pathways | keywords)))
    return vocabulary, rows, {
        "input_provided": True,
        "matched_phenotype_tokens": matched_tokens,
        "pathway_families": sorted(pathways),
        "annotation_keywords": sorted(keywords),
        "status": "applied" if vocabulary else "not_applied_no_shared_route_match",
    }


def interval_reverse_lookup(
    atlas_arg: str,
    out_dir_arg: str,
    reference_function_arg: str | None = None,
    llm_proposals_arg: str | None = None,
    phenotype: str = "",
    case_id: str = "",
    min_members: int = 1,
    min_semantic_coverage: float = 0.0,
    shared_route_arg: str | None = None,
) -> dict[str, object]:
    """Derive family/process hypotheses from what the interval actually contains."""
    if not 0.0 <= min_semantic_coverage <= 1.0:
        raise ValueError("min_semantic_coverage must be between 0 and 1")
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    atlas = _read(atlas_arg, {"gene_id"}, "interval atlas")
    if atlas.empty:
        raise ValueError("interval atlas is empty; run interval-atlas first")
    by_reference, by_focal = _function_lookup(reference_function_arg)
    proposals = _read(llm_proposals_arg, {"hypothesis_id"}, "LLM proposal")

    review_rows: list[dict[str, str]] = []
    lookup_rows: list[dict[str, object]] = []
    phenotype_terms = _phenotype_terms(phenotype)
    shared_route_terms, vocabulary_rows, shared_route_summary = _shared_route_vocabulary(phenotype, shared_route_arg)
    vocabulary_rows = [
        {"vocabulary_source": "raw_phenotype_token", "term": term, "matched_phenotype_token": term}
        for term in sorted(phenotype_terms)
    ] + vocabulary_rows
    for record in atlas.to_dict("records"):
        gene_id = str(record["gene_id"])
        reference = str(record.get("orthology_reference_gene", "")).split(";")[0] or str(record.get("family_label", ""))
        # Any annotated reference works: a model-species homolog, a relative, or
        # the focal reference genome's own annotation keyed on the gene itself.
        function = by_reference.get(reference.casefold(), {})
        origin = "reference_homolog" if function else ""
        if not function:
            function = by_focal.get(gene_id.casefold(), {})
            origin = "focal_reference_annotation" if function else ("homolog_without_function" if reference else "none")
        tier = _tier(record)
        classified = bool(str(record.get("family_id", "")) or str(record.get("domain_architecture", ""))
                          or str(record.get("orthogroup", "")))
        overlay_status, overlap_terms = _phenotype_overlay(phenotype_terms, " ".join(str(value) for value in (
            record.get("family_label", ""), record.get("domain_architecture", ""),
            record.get("broad_function_bins", ""), function.get("function_description", ""),
            function.get("pathway", ""), function.get("mutant_phenotype", ""),
        )))
        shared_overlay_status, shared_overlap_terms = _phenotype_overlay(shared_route_terms, " ".join(str(value) for value in (
            record.get("family_label", ""), record.get("domain_architecture", ""),
            record.get("broad_function_bins", ""), function.get("function_description", ""),
            function.get("pathway", ""), function.get("mutant_phenotype", ""),
        )))
        if shared_route_summary["status"] != "applied":
            shared_overlay_status, shared_overlap_terms = str(shared_route_summary["status"]), ""
        # Keep the corpus's independent provenance axes separate in the B0
        # lookup.  ``reference_source_id`` is a compatibility field for the
        # existing coverage gate; it is a de-duplicated roster, not a vote
        # count.  The semantic bundle also receives each component field.
        source_roster = _tokens_from_reference(function, (
            "source_id", "literature_source_ids", "phenotype_source_ids", "verified_source_ids",
        ))
        lookup_rows.append({
            "gene_id": gene_id, "family_id": str(record.get("family_id", "")),
            "family_label": str(record.get("family_label", "")), "reference_gene": reference,
            "annotation_origin": origin,
            "orthology_relationship": str(record.get("orthology_relationship", "")),
            "orthogroup": str(record.get("orthogroup", "")),
            "domain_architecture": str(record.get("domain_architecture", "")),
            "broad_function_bins": str(record.get("broad_function_bins", "")),
            "relationship_tier": tier,
            "reference_function": str(function.get("function_description", "")),
            "reference_pathway": str(function.get("pathway", "")),
            "reference_pathway_id": str(function.get("pathway_id", "")),
            "reference_mutant_phenotype": str(function.get("mutant_phenotype", "")),
            "reference_source_id": ";".join(source_roster),
            "reference_source_url": str(function.get("source_url", "")),
            "source_url": str(function.get("source_url", "")),
            "literature_source_ids": str(function.get("literature_source_ids", "")),
            "phenotype_source_ids": str(function.get("phenotype_source_ids", "")),
            "verified_source_ids": str(function.get("verified_source_ids", "")),
            "evidence_codes": str(function.get("evidence_codes", "")),
            "evidence_type": str(function.get("evidence_type", "")),
            "provenance_tier": str(function.get("provenance_tier", "")),
            "provenance_detail": str(function.get("provenance_detail", "")),
            "source_database": str(function.get("source_database", "")),
            "source_release": str(function.get("source_release", "")),
            "retrieved_at": str(function.get("retrieved_at", "")),
            "corpus_version": str(function.get("corpus_version", "")),
            "tissue_specificity": str(function.get("tissue_specificity", "")),
            "developmental_stage": str(function.get("developmental_stage", "")),
            "protein_family": str(function.get("protein_family", "")),
            "go_terms": str(function.get("go_terms", "")),
            "keywords": str(function.get("keywords", "")),
            "reference_taxon_id": str(function.get("taxon_id", "")),
            "reference_gene_symbol": str(function.get("gene_symbol", "")),
            "reference_family_id": str(function.get("family_id", "")),
            "reference_orthogroup_id": str(function.get("orthogroup_id", "")),
            "reference_species": str(function.get("reference_species", "")),
            "protein_name": str(function.get("protein_name", "")),
            "reference_domain_ids": str(function.get("domain_ids", "")),
            "reference_domain_description": str(function.get("domain_description", "")),
            "mapping_state": str(function.get("mapping_state", "")),
            "mapping_relationship": str(function.get("mapping_relationship", "")),
            "function_transfer_scope": str(function.get("function_transfer_scope", "")),
            "lookup_status": ("source_backed" if source_roster else
                              ("homolog_only" if classified else "unclassified")),
            "phenotype_overlay_status": overlay_status,
            "phenotype_overlap_terms": overlap_terms,
            "shared_route_overlay_status": shared_overlay_status,
            "shared_route_overlap_terms": shared_overlap_terms,
            "notes": "" if classified else "no informative classification; carried forward as an unknown probe, not dropped",
        })
    lookup = _frame(lookup_rows, LOOKUP_COLUMNS).sort_values("gene_id", kind="stable")

    coverage_masks = {
        "family": lookup["family_id"].astype(str).str.strip().ne(""),
        "orthogroup": lookup["orthogroup"].astype(str).str.strip().ne(""),
        "domain": lookup["domain_architecture"].astype(str).str.strip().ne(""),
        "broad_function": lookup["broad_function_bins"].astype(str).str.strip().ne(""),
        # A source id with nothing behind it is not a reference function (issue
        # #476). The old mask asked only whether ``reference_source_id`` was
        # non-empty, so an atlas whose reference rows carry a product name and a
        # pipeline label (``local_fasta:<panel>``) -- no FUNCTION, no PATHWAY, no
        # DISRUPTION PHENOTYPE -- still reported ``source_backed_reverse_lookup``,
        # the strongest state, and was indistinguishable from a run backed by
        # curated UniProt text. Membership now requires a declared source *and*
        # at least one piece of the functional content the reverse lookup exists
        # to read. Whether that source is independently citable is a separate
        # axis, reported below and never merged into this one.
        "source_backed_reference_function": (
            lookup["reference_source_id"].astype(str).str.strip().ne("")
            & (
                lookup["reference_function"].astype(str).str.strip().ne("")
                | lookup["reference_pathway"].astype(str).str.strip().ne("")
                | lookup["reference_mutant_phenotype"].astype(str).str.strip().ne("")
            )
        ),
    }
    semantic_mask = pd.Series(False, index=lookup.index)
    for mask in coverage_masks.values():
        semantic_mask = semantic_mask | mask.reset_index(drop=True)
    semantically_annotated_genes = int(semantic_mask.sum())
    semantic_coverage_fraction = round(semantically_annotated_genes / max(1, len(lookup)), 6)
    semantically_evaluable = bool(
        semantically_annotated_genes
        and semantic_coverage_fraction >= min_semantic_coverage
    )

    # Issue #216: `semantically_evaluable` counts a broad-function bin as usable
    # annotation, which is right for "did this command have anything to group by"
    # and wrong for "did the intended reverse lookup happen". A run whose only
    # signal is broad bins never performed
    # interval gene -> homolog -> source-backed function, so it is reported
    # separately and must not read as an equivalent Route B result.
    structural_mask = (coverage_masks["family"] | coverage_masks["orthogroup"] | coverage_masks["domain"]).reset_index(drop=True)
    structural_genes = int(structural_mask.sum())
    source_backed_genes_count = int(coverage_masks["source_backed_reference_function"].sum())
    broad_only_genes = int(coverage_masks["broad_function"].sum())
    if source_backed_genes_count and structural_genes:
        reverse_lookup_status, process_inference_completeness = "source_backed_reverse_lookup", "complete"
    elif source_backed_genes_count or structural_genes:
        reverse_lookup_status = "source_backed_function_without_homology_structure" if source_backed_genes_count else "structural_providers_without_source_backed_function"
        process_inference_completeness = "partial"
    elif broad_only_genes:
        reverse_lookup_status, process_inference_completeness = "broad_bin_only", "not_evaluated_for_reverse_lookup"
    else:
        reverse_lookup_status, process_inference_completeness = "no_provider", "not_evaluated_for_reverse_lookup"

    hypotheses: list[dict[str, object]] = []
    for family_id, group in lookup[lookup.family_id.ne("")].groupby("family_id", sort=True):
        records = group.to_dict("records")
        if len(records) < min_members:
            continue
        pathways = sorted({str(row["reference_pathway"]) for row in records if str(row["reference_pathway"])})
        process_ids = sorted({str(row["reference_pathway_id"]) for row in records
                              if str(row.get("reference_pathway_id", ""))})
        bins = sorted({value for row in records for value in str(row["broad_function_bins"]).split(";") if value})
        sources = sorted({str(row["reference_source_id"]) for row in records if str(row["reference_source_id"])})
        hypotheses.append({
            "hypothesis_id": f"rev_fam_{safe_id(str(family_id))}", "route": ROUTE_B, "scope": "family",
            "pathway_or_process": "; ".join(pathways) or (bins[0] if bins else "unknown_process"),
            "process_id": "; ".join(process_ids),
            "gene_family": str(family_id),
            "reference_genes": ";".join(sorted({str(row["reference_gene"]) for row in records if str(row["reference_gene"])})),
            "member_interval_genes": ";".join(sorted(str(row["gene_id"]) for row in records)),
            "member_count": len(records), "anchor_status": "interval_anchored",
            "relationship_tier": _best_tier([str(row["relationship_tier"]) for row in records]),
            "source_backing": ";".join(sources),
            "proposal_source": "interval_reverse_lookup",
            "counterevidence": "",
            "verification_status": "source_backed" if sources else "pending_verification",
            "rank_contribution": "none",
            "notes": "derived from the interval, not from a phenotype-only guess",
        })
    for orthogroup, group in lookup[lookup.orthogroup.ne("")].groupby("orthogroup", sort=True):
        records = group.to_dict("records")
        if len(records) < min_members:
            continue
        pathways = sorted({str(row["reference_pathway"]) for row in records if str(row["reference_pathway"])})
        process_ids = sorted({str(row["reference_pathway_id"]) for row in records if str(row["reference_pathway_id"])})
        sources = sorted({str(row["reference_source_id"]) for row in records if str(row["reference_source_id"])})
        hypotheses.append({
            "hypothesis_id": f"rev_og_{safe_id(str(orthogroup))}", "route": ROUTE_B, "scope": "orthogroup",
            "pathway_or_process": "; ".join(pathways) or "unknown_process", "process_id": "; ".join(process_ids),
            "gene_family": str(orthogroup),
            "reference_genes": ";".join(sorted({str(row["reference_gene"]) for row in records if str(row["reference_gene"])})),
            "member_interval_genes": ";".join(sorted(str(row["gene_id"]) for row in records)),
            "member_count": len(records), "anchor_status": "interval_anchored", "relationship_tier": "same_subfamily",
            "source_backing": ";".join(sources), "proposal_source": "interval_reverse_lookup",
            "counterevidence": "", "verification_status": "source_backed" if sources else "pending_verification",
            "rank_contribution": "none",
            "notes": "orthogroup is a supported evolutionary grouping from the interval; it is not causal evidence",
        })
    bin_members: dict[str, list[dict[str, object]]] = {}
    for row in lookup.to_dict("records"):
        for bin_id in (value for value in str(row["broad_function_bins"]).split(";") if value):
            bin_members.setdefault(bin_id, []).append(row)
    for bin_id, records in sorted(bin_members.items()):
        if len(records) < max(2, min_members):
            continue
        hypotheses.append({
            "hypothesis_id": f"rev_bin_{safe_id(bin_id)}", "route": ROUTE_B, "scope": "broad_process",
            "pathway_or_process": bin_id, "gene_family": "",
            "reference_genes": ";".join(sorted({str(row["reference_gene"]) for row in records if str(row["reference_gene"])})),
            "member_interval_genes": ";".join(sorted(str(row["gene_id"]) for row in records)),
            "member_count": len(records), "anchor_status": "interval_anchored",
            "relationship_tier": "same_broad_process",
            "source_backing": ";".join(sorted({str(row["reference_source_id"]) for row in records if str(row["reference_source_id"])})),
            "proposal_source": "interval_reverse_lookup", "counterevidence": "",
            "verification_status": "pending_verification", "rank_contribution": "none",
            "notes": "broad-process grouping is the weakest tier and is kept for recall, not precision",
        })
    unclassified = lookup[lookup.lookup_status.eq("unclassified")]
    if len(unclassified):
        hypotheses.append({
            "hypothesis_id": "rev_unknown_probe", "route": ROUTE_B, "scope": "unknown_probe",
            "pathway_or_process": "unknown / conserved-unknown proteins", "gene_family": "",
            "reference_genes": "", "member_interval_genes": ";".join(sorted(unclassified["gene_id"].astype(str))),
            "member_count": int(len(unclassified)), "anchor_status": "interval_anchored",
            "relationship_tier": "unclassified", "source_backing": "",
            "proposal_source": "interval_reverse_lookup", "counterevidence": "",
            "verification_status": "pending_verification", "rank_contribution": "none",
            "notes": "unknown / orphan / dark interval genes are never removed by this route",
        })
        review_rows.append({"reason": "unclassified_interval_genes_carried_forward", "target": str(len(unclassified)),
                            "required_next_action": "run annotation-rescue or add a closer reference proteome; do not treat these as excluded",
                            "allowed_choices": "annotation rescue | add proteome | keep as an explicit unknown probe"})
    reference_symbols = {str(value).casefold() for value in lookup["reference_gene"] if str(value)}
    family_labels = {str(value).casefold() for value in lookup["family_label"] if str(value)}
    for record in proposals.to_dict("records"):
        family = str(record.get("gene_family", ""))
        reference_name = str(record.get("reference_gene_name", "") or record.get("proposed_stable_id", ""))
        members = sorted(lookup[lookup.family_id.astype(str).str.casefold() == family.casefold()]["gene_id"].astype(str)) if family else []
        if not members and reference_name:
            members = sorted(lookup[lookup.reference_gene.astype(str).str.casefold() == reference_name.casefold()]["gene_id"].astype(str))
        # Anti-hallucination anchor: a proposal naming something that is not in the
        # interval listing is retained but marked, never counted as a finding.
        anchored = bool(members) or family.casefold() in family_labels or reference_name.casefold() in reference_symbols
        hypotheses.append({
            "hypothesis_id": f"rev_llm_{safe_id(str(record['hypothesis_id']))}", "route": ROUTE_B, "scope": "llm_proposal",
            "pathway_or_process": str(record.get("pathway_or_process", "")), "gene_family": family,
            "reference_genes": reference_name,
            "member_interval_genes": ";".join(members), "member_count": len(members),
            "anchor_status": "interval_anchored" if anchored else "unanchored",
            "relationship_tier": "same_family" if members else "unclassified",
            "source_backing": "", "proposal_source": str(record.get("proposal_source", "")) or "llm",
            "counterevidence": str(record.get("counterevidence", "")),
            "verification_status": "pending_verification" if anchored else "rejected_unanchored",
            "rank_contribution": "none",
            "notes": ("LLM proposal from the interval-informed lane; never blind discovery performance" if anchored
                      else "names a gene/family that is not in the interval listing: the hallucination guard"),
        })
        if anchored:
            review_rows.append({"reason": "llm_proposal_requires_verification", "target": str(record["hypothesis_id"]),
                                "required_next_action": "run hypothesis-verify before this proposal may become a bait",
                                "allowed_choices": "verify | discard | keep as pending"})
        else:
            review_rows.append({"reason": "llm_proposal_not_anchored_in_interval", "target": str(record["hypothesis_id"]),
                                "required_next_action": "区間の一覧に無い遺伝子/familyを挙げている。捏造の可能性があるので、baitへは進めずproposalとしてのみ保持する",
                                "allowed_choices": "破棄 | 区間外仮説として別扱い | 区間を見直す"})

    hypothesis_frame = _frame(hypotheses, HYPOTHESIS_COLUMNS).fillna("").sort_values("hypothesis_id", kind="stable")
    if not hypothesis_frame.empty:
        overlay_by_gene = {
            str(row["gene_id"]): {
                "raw": (str(row["phenotype_overlay_status"]), str(row["phenotype_overlap_terms"])),
                "shared": (str(row["shared_route_overlay_status"]), str(row["shared_route_overlap_terms"])),
            }
            for row in lookup.to_dict("records")
        }
        for index, row in hypothesis_frame.iterrows():
            members = _genes(row.get("member_interval_genes", ""))
            raw_statuses = {overlay_by_gene.get(gene, {"raw": ("not_applied", "")})["raw"][0] for gene in members}
            raw_terms = sorted({term for gene in members for term in overlay_by_gene.get(gene, {"raw": ("", "")})["raw"][1].split(";") if term})
            shared_statuses = {overlay_by_gene.get(gene, {"shared": ("not_applied_no_shared_route", "")})["shared"][0] for gene in members}
            shared_terms = sorted({term for gene in members for term in overlay_by_gene.get(gene, {"shared": ("", "")})["shared"][1].split(";") if term})
            if not phenotype_terms:
                raw_status = "not_applied"
            elif "phenotype_consistent" in raw_statuses:
                raw_status = "phenotype_consistent"
            else:
                raw_status = "phenotype_not_matched"
            if shared_route_summary["status"] != "applied":
                shared_status = str(shared_route_summary["status"])
            elif "phenotype_consistent" in shared_statuses:
                shared_status = "phenotype_consistent"
            else:
                shared_status = "phenotype_not_matched"
            hypothesis_frame.at[index, "phenotype_overlay_status"] = raw_status
            hypothesis_frame.at[index, "phenotype_overlap_terms"] = ";".join(raw_terms)
            hypothesis_frame.at[index, "shared_route_overlay_status"] = shared_status
            hypothesis_frame.at[index, "shared_route_overlap_terms"] = ";".join(shared_terms)
    if hypothesis_frame.empty:
        review_rows.append({"reason": "reverse_lookup_produced_no_hypothesis", "target": case_id or "case",
                            "required_next_action": "the interval carries no usable classification; run interval-atlas with orthology/domain inputs or annotation-rescue",
                            "allowed_choices": "add classification inputs | annotation rescue | abstain"})
    if process_inference_completeness != "complete":
        review_rows.append({
            "reason": f"route_b_reverse_lookup_not_complete:{reverse_lookup_status}",
            "target": case_id or "case",
            "required_next_action": (
                "run `route-b-prepare` to build an orthology / domain / source-backed reference-function "
                "atlas from the local resources; broad-function bins alone are not a reverse lookup"
            ),
            "allowed_choices": "route-b-prepare | add providers | annotation rescue | mark not_evaluated_for_reverse_lookup",
        })
    if not semantically_evaluable:
        review_rows.append({
            "reason": "route_b_semantic_inputs_insufficient",
            "target": case_id or "case",
            "required_next_action": (
                "build the interval atlas with family/orthogroup/domain/function providers; "
                "command completion alone is not a Route B benchmark"
            ),
            "allowed_choices": "add providers | annotation rescue | mark not_evaluated_provider_missing",
        })

    prompt_path, prompts_path = _write_reverse_prompt(out, lookup, phenotype, case_id)
    reviews = _frame(review_rows, REVIEW_COLUMNS).drop_duplicates().sort_values(["reason", "target"], kind="stable") if review_rows else pd.DataFrame(columns=REVIEW_COLUMNS)

    paths = {"interval_reverse_lookup": out / "interval_reverse_lookup.tsv",
             "family_process_hypotheses": out / "family_process_hypotheses.tsv",
             "reverse_lookup_review_required": out / "reverse_lookup_review_required.tsv",
             "phenotype_overlay_vocabulary": out / "phenotype_overlay_vocabulary.tsv"}
    lookup.to_csv(paths["interval_reverse_lookup"], sep="\t", index=False, lineterminator="\n")
    hypothesis_frame.to_csv(paths["family_process_hypotheses"], sep="\t", index=False, lineterminator="\n")
    reviews.to_csv(paths["reverse_lookup_review_required"], sep="\t", index=False, lineterminator="\n")
    _frame(vocabulary_rows, OVERLAY_VOCABULARY_COLUMNS).sort_values(
        ["vocabulary_source", "term", "matched_phenotype_token"], kind="stable"
    ).to_csv(paths["phenotype_overlay_vocabulary"], sep="\t", index=False, lineterminator="\n")

    lane_variant = "b1_phenotype_overlay" if phenotype else "b0_phenotype_free"
    semantic_status = "evaluable" if semantically_evaluable else "not_evaluated_provider_missing"
    coverage = {
        f"{name}_gene_count": int(mask.sum())
        for name, mask in coverage_masks.items()
    }
    # The status words say whether the lookup could run. These counts say what it
    # actually had to read, so "source-backed" can be checked rather than trusted
    #. ``citable`` applies the same source grading the auto-triage
    # lane uses; a focal-genome annotation marker such as
    # ``local:annotation`` is a legitimate reference here but is not citable, so
    # the two are counted separately and neither is scored.
    declared_source_gene_count = int(lookup["reference_source_id"].astype(str).str.strip().ne("").sum())
    citable_source_gene_count = int(
        lookup["reference_source_id"].map(lambda value: classify_source(str(value or "").strip()) in {"reference", "local_file"}).sum()
    )
    coverage.update({
        "semantically_annotated_gene_count": semantically_annotated_genes,
        "semantic_coverage_fraction": semantic_coverage_fraction,
        "min_semantic_coverage": min_semantic_coverage,
        "declared_reference_source_gene_count": declared_source_gene_count,
        "declared_source_without_function_text_gene_count": (
            declared_source_gene_count - int(coverage_masks["source_backed_reference_function"].sum())
        ),
        "citable_reference_source_gene_count": citable_source_gene_count,
        "reference_function_text_gene_count": int(lookup["reference_function"].astype(str).str.strip().ne("").sum()),
        "reference_pathway_text_gene_count": int(lookup["reference_pathway"].astype(str).str.strip().ne("").sum()),
        "reference_mutant_phenotype_text_gene_count": int(
            lookup["reference_mutant_phenotype"].astype(str).str.strip().ne("").sum()
        ),
    })
    overlay_counts = {
        str(key): int(value)
        for key, value in lookup["phenotype_overlay_status"].value_counts().sort_index().items()
    }
    shared_overlay_counts = {
        str(key): int(value)
        for key, value in lookup["shared_route_overlay_status"].value_counts().sort_index().items()
    }
    # Same definition as the coverage mask, so process_inference_status and
    # reverse_lookup_status cannot disagree about what "source backed" means.
    source_backed_count = int(coverage_masks["source_backed_reference_function"].sum())
    if not reference_function_arg:
        process_inference_status = "unavailable_no_reference_function"
    elif source_backed_count:
        process_inference_status = "available_source_backed"
    else:
        process_inference_status = "attempted_no_source_backed_function"

    # Issue #305: family-label quality. `semantically_evaluable` only asks
    # whether *any* provider annotated the interval; it never asks whether the
    # family labels are usable gene symbols. A TrEMBL-derived atlas whose
    # family_label values are all UniProt accessions still passes the coverage
    # gate yet is not a symbol-resolvable Route B reference. Quality is judged
    # on family_label (the symbol side, best_subject_symbol), not on family_id
    # (the stem side), because identifier-shaped values land in the symbol
    # column when no curated symbol exists; the shared identifier grammar
    # (identifiers.is_usable_symbol) is the single source of truth and no
    # regex is duplicated here. Thresholds are fixed scientific-sanity values
    # and are never tuned against benchmark truth.
    family_label_values = lookup["family_label"].astype(str).str.strip()
    family_annotated_genes = int(family_label_values.ne("").sum())
    usable_symbol_genes = int(sum(is_usable_symbol(value) for value in family_label_values if value))
    usable_symbol_fraction = round(usable_symbol_genes / max(1, family_annotated_genes), 6)
    accession_or_locus_tag_fraction = round(1.0 - usable_symbol_fraction, 6)
    unresolved_family_fraction = round((int(len(lookup)) - usable_symbol_genes) / max(1, len(lookup)), 6)
    if family_annotated_genes == 0:
        family_label_quality = "not_evaluable"
    elif usable_symbol_fraction < 0.5:
        family_label_quality = "degraded"
    else:
        family_label_quality = "healthy"
    semantic_quality = {
        "usable_symbol_fraction": usable_symbol_fraction,
        "accession_or_locus_tag_fraction": accession_or_locus_tag_fraction,
        "unresolved_family_fraction": unresolved_family_fraction,
        "family_label_quality": family_label_quality,
        "rule": (
            "family_label_quality is judged on the family_label column (the symbol side, "
            "best_subject_symbol); a label is usable when is_usable_symbol() matches no "
            "identifier grammar (UniProt accession, AGI locus, locus tag). "
            "usable_symbol_fraction = usable family labels / genes with a non-empty family_label "
            "(family-annotated denominator); accession_or_locus_tag_fraction = (family-annotated "
            "- usable) / family-annotated; unresolved_family_fraction = (interval genes - usable) "
            "/ interval genes, so a gene with no family label counts as unresolved. "
            "family_label_quality = not_evaluable when no gene carries a family_label (consistent "
            "with provider missing); degraded when usable_symbol_fraction < 0.5 (including 0: all "
            "labels are accessions or locus tags); healthy when usable_symbol_fraction >= 0.5. "
            "Thresholds are fixed scientific-sanity values and are never tuned against benchmark truth."
        ),
    }

    status = {
        "schema_version": SCHEMA_VERSION, "route": ROUTE_B, "case_id": case_id,
        "command_completed": True,
        "lane_variant": lane_variant,
        "lane": "interval_informed",
        "interval_information_in_model_input": True,
        "phenotype_used": bool(phenotype),
        "blind_discovery_claim": False,
        "semantically_evaluable": semantically_evaluable,
        "semantic_status": semantic_status,
        "semantic_evaluability_reason": (
            "the interval atlas contains usable family/orthogroup/domain/broad-function/source-backed function annotations"
            if semantically_evaluable else
            "the command completed, but no interval gene has a usable Route B semantic provider annotation"
        ),
        "provider_coverage": coverage,
        "semantic_quality": semantic_quality,
        "process_inference_status": process_inference_status,
        # Strict Route B evaluability. `semantic_status` answers
        # "was there anything to group by"; these two answer "did the intended
        # interval-gene -> homolog -> source-backed-function lookup happen".
        "reverse_lookup_status": reverse_lookup_status,
        "process_inference_completeness": process_inference_completeness,
        "reverse_lookup_evaluability": {
            "structural_provider_gene_count": structural_genes,
            "source_backed_function_gene_count": source_backed_genes_count,
            "broad_function_gene_count": broad_only_genes,
            "rule": (
                "complete requires at least one family/orthogroup/domain provider AND at least one "
                "source-backed reference function; broad-function bins alone are "
                "not_evaluated_for_reverse_lookup and are never equivalent to a source-backed Route B result"
            ),
            "missing_provider_is_not_negative_evidence": True,
        },
        "phenotype_overlay": {
            "status_counts": overlay_counts,
            "raw_phenotype_terms": sorted(phenotype_terms),
            "raw_phenotype_status_counts": overlay_counts,
            "shared_route": {
                **shared_route_summary,
                "status_counts": shared_overlay_counts,
                "rule": "pathway_family and annotation_keywords are derived only through phenotype-token matches in the supplied shared route",
            },
            "non_scoring": True,
            "rule": "raw-phenotype and shared-route overlays are reported separately; neither changes tier, rank, or score",
        },
        "relationship_tier_resolution": {
            "one_to_one_gene_count": int(lookup["relationship_tier"].eq("exact_ortholog").sum()),
            "orthogroup_gene_count": int(coverage_masks["orthogroup"].sum()),
            "family_gene_count": int(coverage_masks["family"].sum()),
            "domain_gene_count": int(coverage_masks["domain"].sum()),
            "broad_function_gene_count": int(coverage_masks["broad_function"].sum()),
            "limiting_inputs": [
                name for name, count in (
                    ("orthology_relationship", int(lookup["orthology_relationship"].astype(str).str.strip().ne("").sum())),
                    ("orthogroup", int(coverage_masks["orthogroup"].sum())),
                    ("family", int(coverage_masks["family"].sum())),
                    ("domain", int(coverage_masks["domain"].sum())),
                ) if count == 0
            ],
        },
        "separation_rule": "this lane sees the interval; its results must never be reported as blind phenotype-only discovery",
        "interval_gene_count": int(len(lookup)),
        "source_backed_genes": int(lookup.lookup_status.eq("source_backed").sum()),
        "unclassified_genes": int(lookup.lookup_status.eq("unclassified").sum()),
        "hypothesis_count": int(len(hypothesis_frame)),
        "hypotheses_by_scope": {str(k): int(v) for k, v in hypothesis_frame["scope"].value_counts().sort_index().items()} if not hypothesis_frame.empty else {},
        "relationship_tier_counts": {str(k): int(v) for k, v in hypothesis_frame["relationship_tier"].value_counts().sort_index().items()} if not hypothesis_frame.empty else {},
        "llm_proposals_imported": int(len(proposals)),
        "unanchored_proposals": int(hypothesis_frame["anchor_status"].eq("unanchored").sum()) if not hypothesis_frame.empty else 0,
        "annotation_origin_counts": {str(k): int(v) for k, v in lookup["annotation_origin"].value_counts().sort_index().items()},
        "risk": {
            "failure_mode": "post-hoc rationalisation: everything proposed exists, but a plausible story can be told about almost any gene once you are looking at it",
            "unclassified_gene_fraction": round(int(lookup.lookup_status.eq("unclassified").sum()) / max(1, len(lookup)), 4),
            "post_hoc_rationalization_risk": "structural: mitigated by never showing this route the phenotype-first answer and by scoring it in its own lane",
            "hallucination_guard": "a proposal naming a gene or family absent from the interval listing is marked unanchored and rejected as a bait",
            "mitigation": "compare with the phenotype-first route; only strict mechanistic convergence has both mirror-image risks reduced at once (issue #186)",
        },
        "prompt": str(prompt_path), "prompts_jsonl": str(prompts_path),
        "inputs": {"atlas": file_fingerprint(runtime_path(atlas_arg)),
                   "reference_function": file_fingerprint(runtime_path(reference_function_arg)) if reference_function_arg else None,
                   "llm_proposals": file_fingerprint(runtime_path(llm_proposals_arg)) if llm_proposals_arg else None,
                   "shared_route": file_fingerprint(runtime_path(shared_route_arg)) if shared_route_arg else None},
        "safety": [
            "adjacent families and processes are retained but carry a weaker relationship_tier, never the same strength as an exact ortholog",
            "family or process membership alone never changes a gene score",
            "unknown / orphan / dark interval genes are carried forward as an explicit probe, never removed",
            "LLM proposals stay pending_verification and are not canonical evidence",
        ],
        "outputs": {name: str(path) for name, path in sorted(paths.items())},
        "frozen_output_fingerprints": {
            name: file_fingerprint(path, full_hash=True)
            for name, path in sorted(paths.items())
        },
    }
    status_path = out / "interval_first_status.json"
    status_path.write_text(json.dumps(status, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"status": str(status_path), "hypothesis_count": int(len(hypothesis_frame)),
            "interval_gene_count": int(len(lookup)), "prompt": str(prompt_path),
            "prompts_jsonl": str(prompts_path),
            **{name: str(path) for name, path in paths.items()}}


def _write_reverse_prompt(out: Path, lookup: pd.DataFrame, phenotype: str, case_id: str) -> tuple[Path, Path]:
    """Prompt package for the interval-informed LLM lane.

    The model is shown the model-species homologs present in the interval - not
    the focal gene ids, and never asked to rank candidates.
    """
    lines = [f"Case: {case_id or 'unnamed'}",
             f"Phenotype: {phenotype or 'not supplied (pure reverse lookup)'}",
             "", f"The {len(lookup)} gene(s) present in this QTL interval, with whatever annotated",
             "identity is available (model-species homolog, relative, or the reference's own annotation):"]
    for row in lookup.to_dict("records"):
        reference = str(row["reference_gene"]) or "(no annotated identity)"
        detail = "; ".join(value for value in (
            str(row["family_label"]) and f"family {row['family_label']}",
            str(row["domain_architecture"]) and f"domains {row['domain_architecture']}",
            str(row["orthology_relationship"]) and f"orthology {row['orthology_relationship']}",
            str(row["reference_function"]) and f"function {row['reference_function']}",
        ) if value)
        lines.append(f"- {reference}{': ' + detail if detail else ''}")
    lines += ["", 'Return JSON only:',
              '{"hypotheses": [{"pathway_or_process": "", "gene_family": "", "reference_gene_name": "",'
              ' "reference_species": "", "proposed_stable_id": "", "expected_role": "", "mechanistic_link": "",'
              ' "counterevidence": "", "confidence": "low|medium|high", "uncertainty": ""}]}']
    user = "\n".join(lines)
    prompt_path = out / "reverse_lookup_prompt.txt"
    prompt_path.write_text(f"{REVERSE_SYSTEM_PROMPT}\n\n{user}\n", encoding="utf-8")
    prompts_path = out / "reverse_lookup_prompts.jsonl"
    prompts_path.write_text(json.dumps({
        "prompt_id": f"reverse::{safe_id(case_id or 'case')}",
        "gene_id": f"interval_reverse_lookup:{case_id or 'case'}",
        "system": REVERSE_SYSTEM_PROMPT, "user": user,
    }, ensure_ascii=False) + "\n", encoding="utf-8")
    return prompt_path, prompts_path


# --------------------------------------------------------------------------- #
# comparison
# --------------------------------------------------------------------------- #

def _normalized_text(value: object) -> str:
    return " ".join(str(value or "").split()).strip().strip(".,;:!?()[]\"'").casefold()


def _vocabulary_concept(text: str) -> str:
    for concept, (aliases, _parent, _ids) in PROCESS_VOCABULARY.items():
        if text in aliases:
            return concept
    return ""


def _vocabulary_concept_for_id(stable_id: str) -> str:
    """Map a stable database ID to a curated concept when the registry knows it."""
    folded = stable_id.casefold()
    for concept, (_aliases, _parent, known_ids) in PROCESS_VOCABULARY.items():
        if folded in known_ids:
            return concept
    return ""


def _process_identity(record: dict[str, object]) -> tuple[str, str, str]:
    """(identity key, normalized free text, vocabulary concept) for a process mention.

    Identity order: explicit ``process_id``/``pathway_id`` column, then a stable
    database ID embedded in the free text, then the curated vocabulary, then the
    normalized free text itself as a weak alias key. A stable ID that the curated
    vocabulary knows resolves to the same concept as its aliases; unknown IDs stay
    ``id:...``. Free text is never promoted to a stable identity on its own
    (issue #186).
    """
    text = _normalized_text(record.get("pathway_or_process"))
    explicit = _normalized_text(record.get("process_id") or record.get("pathway_id"))
    if explicit:
        concept = _vocabulary_concept_for_id(explicit)
        if concept:
            return f"concept:{concept}", text, concept
        return f"id:{explicit}", text, ""
    embedded = _STABLE_ID_PATTERN.search(text)
    if embedded:
        concept = _vocabulary_concept_for_id(embedded.group(0))
        if concept:
            return f"concept:{concept}", text, concept
        return f"id:{embedded.group(0).upper()}", text, ""
    concept = _vocabulary_concept(text)
    if concept:
        return f"concept:{concept}", text, concept
    return f"text:{text}", text, ""


def _family_keys(record: dict[str, object]) -> tuple[set[str], set[str]]:
    """Auditable identity-key sets for a hypothesis family and orthogroup."""
    return (family_identity_keys(record.get("gene_family")),
            family_identity_keys(record.get("orthogroup")))


def _source_ids(record: dict[str, object]) -> set[str]:
    """Reference-record ids behind a route (e.g. UniProt:Q8GRU6, PMID:1001).

    Status strings such as ``verified_bait`` carry no record id and are dropped,
    so a bait status alone never counts as a shared source chain.
    """
    ids: set[str] = set()
    for column in ("source_backing", "source_id"):
        for token in str(record.get(column, "")).replace(",", ";").split(";"):
            token = token.strip()
            if token and ":" in token:
                ids.add(token)
    return ids


def _shared_source_ids(source_a: set[str], source_b: set[str]) -> list[str]:
    folded_b = {value.casefold() for value in source_b}
    return sorted(token for token in source_a if token.casefold() in folded_b)


def _independence_status(source_a: set[str], source_b: set[str]) -> str:
    if _shared_source_ids(source_a, source_b):
        return "shared_source_chain"
    if source_a and source_b:
        return "independent_source_chains"
    return "source_ids_incomplete"


def _convergence_status(
    *,
    process_match: bool,
    stable_process: bool,
    family_match: bool,
    shared_genes: bool,
    backed_a: bool,
    backed_b: bool,
    process_related: bool,
    text_equal: bool,
) -> tuple[str, str]:
    """Pair-level convergence tier (issue #186).

    Tier 1 requires a stable same process concept AND family/subfamily or interval
    member consistency AND source backing on both routes. Without backing, a
    same-concept pair drops to the family / gene-overlap tier below it.
    """
    if process_match and stable_process and (family_match or shared_genes) and backed_a and backed_b:
        return "mechanistic_convergence", "tier1_mechanistic"
    if family_match:
        return "family_or_subfamily_convergence", "tier2_family"
    if shared_genes:
        return "gene_overlap_only", "tier3_gene_overlap"
    if process_related or (process_match and stable_process):
        return "process_related", "tier4_process_related"
    if text_equal:
        return "text_only_match", "tier5_text_only"
    return "", ""


def _highest_tier(tiers: list[str]) -> str:
    for tier in MATCH_TIERS:
        if tier in tiers:
            return tier
    return ""


def _read_route_b_phenotype_seen(manifest_path_arg: str | None) -> bool | None:
    """Record whether Route B saw the phenotype (B0 = false, B1 = true).

    Reads ``phenotype_used`` from the interval_reverse_lookup status manifest;
    ``None`` means the caller did not supply the manifest.
    """
    if not manifest_path_arg:
        return None
    manifest_path = runtime_path(manifest_path_arg)
    if not manifest_path.exists():
        raise ValueError(f"route B manifest not found: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if "phenotype_used" not in manifest:
        return None
    return bool(manifest["phenotype_used"])


def dual_route_compare(
    route_a_arg: str,
    route_b_arg: str,
    out_dir_arg: str,
    case_id: str = "",
    route_b_manifest_arg: str | None = None,
    benchmark_mode: bool = False,
) -> dict[str, object]:
    """Compare the two routes only after both have been generated independently.

    Convergence is decomposed into tiers (strict mechanistic > family/subfamily >
    shared interval gene > process-related > free-text-only) so a shared interval
    gene alone is never presented as a shared mechanism. The artifact records the
    source chains behind each route and whether Route B saw the phenotype
    (B0 phenotype-blind vs B1 phenotype overlay, issue #186).
    """
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    route_a = _read(route_a_arg, {"hypothesis_id"}, "route A hypothesis")
    route_b = _read(route_b_arg, {"hypothesis_id"}, "route B hypothesis")
    phenotype_seen_by_route_b = _read_route_b_phenotype_seen(route_b_manifest_arg)

    rows: list[dict[str, object]] = []
    matched_a: set[str] = set()
    matched_b: set[str] = set()
    for record_a in route_a.to_dict("records"):
        genes_a = _genes(record_a.get("member_interval_genes", ""))
        family_a, orthogroup_a = _family_keys(record_a)
        process_key_a, text_a, concept_a = _process_identity(record_a)
        source_a = _source_ids(record_a)
        backed_a = bool(str(record_a.get("source_backing", "")).strip())
        for record_b in route_b.to_dict("records"):
            genes_b = _genes(record_b.get("member_interval_genes", ""))
            family_b, orthogroup_b = _family_keys(record_b)
            process_key_b, text_b, concept_b = _process_identity(record_b)
            source_b = _source_ids(record_b)
            backed_b = bool(str(record_b.get("source_backing", "")).strip())
            shared = sorted(genes_a & genes_b)
            family_match = bool((family_a and family_a & family_b)
                                or (orthogroup_a and orthogroup_a & orthogroup_b))
            process_match = bool(process_key_a) and process_key_a == process_key_b
            stable_process = process_match and not process_key_a.startswith("text:")
            text_equal = bool(text_a) and text_a == text_b
            parent_a = PROCESS_VOCABULARY[concept_a][1] if concept_a in PROCESS_VOCABULARY else ""
            parent_b = PROCESS_VOCABULARY[concept_b][1] if concept_b in PROCESS_VOCABULARY else ""
            process_related = (not process_match) and bool(parent_a) and parent_a == parent_b

            basis: list[str] = []
            if shared:
                basis.append("shared_interval_gene")
            if family_match:
                basis.append("same_gene_family")
            if process_match:
                basis.append("same_process_concept" if stable_process else "same_pathway_text")
            elif process_related:
                basis.append("process_related_concept")
            elif text_equal:
                basis.append("same_pathway_text")
            if not basis:
                continue
            status, tier = _convergence_status(
                process_match=process_match, stable_process=stable_process,
                family_match=family_match, shared_genes=bool(shared),
                backed_a=backed_a, backed_b=backed_b,
                process_related=process_related, text_equal=text_equal,
            )
            shared_source = _shared_source_ids(source_a, source_b)
            matched_a.add(str(record_a["hypothesis_id"]))
            matched_b.add(str(record_b["hypothesis_id"]))
            rows.append({
                "comparison_id": f"cmp_{len(rows):03d}",
                "route_a_hypothesis_id": str(record_a["hypothesis_id"]),
                "route_b_hypothesis_id": str(record_b["hypothesis_id"]),
                "pathway_or_process": str(record_a.get("pathway_or_process", "")) or str(record_b.get("pathway_or_process", "")),
                "gene_family": str(record_a.get("gene_family", "")) or str(record_b.get("gene_family", "")),
                "shared_interval_genes": ";".join(shared), "match_basis": ";".join(basis),
                "match_tier": tier,
                "relationship_tier": _best_tier([str(record_a.get("relationship_tier", "unclassified")),
                                                 str(record_b.get("relationship_tier", "unclassified"))]),
                "route_a_source_backing": str(record_a.get("source_backing", "")),
                "route_b_source_backing": str(record_b.get("source_backing", "")),
                "convergence_status": status, "rank_contribution": "none",
                "route_a_source_ids": ";".join(sorted(source_a)),
                "route_b_source_ids": ";".join(sorted(source_b)),
                "shared_source_ids": ";".join(shared_source),
                "source_independence_status": _independence_status(source_a, source_b),
                "interval_overlap_status": "shared_interval_genes" if shared else "no_shared_interval_genes",
                "route_a_anchor_status": str(record_a.get("anchor_status", "")),
                "notes": "convergence is agreement between two independent routes, not additional causal evidence",
            })
    for record_a in route_a.to_dict("records"):
        if str(record_a["hypothesis_id"]) in matched_a:
            continue
        rows.append({
            "comparison_id": f"cmp_{len(rows):03d}", "route_a_hypothesis_id": str(record_a["hypothesis_id"]),
            "route_b_hypothesis_id": "", "pathway_or_process": str(record_a.get("pathway_or_process", "")),
            "gene_family": str(record_a.get("gene_family", "")), "shared_interval_genes": "", "match_basis": "",
            "match_tier": "", "relationship_tier": str(record_a.get("relationship_tier", "unclassified")),
            "route_a_source_backing": str(record_a.get("source_backing", "")), "route_b_source_backing": "",
            "convergence_status": "phenotype_only_supported", "rank_contribution": "none",
            "route_a_source_ids": "", "route_b_source_ids": "", "shared_source_ids": "",
            "source_independence_status": "not_applicable", "interval_overlap_status": "no_shared_interval_genes",
            # Carried so a downstream reader can tell "Route A named something that is
            # not in the interval listing at all" (unanchored) apart from "Route A is
            # anchored in the interval, Route B simply never reached it". Both are
            # phenotype_only_supported, but only the first is evidence about interval
            # membership; collapsing them would overstate what the disagreement shows.
            "route_a_anchor_status": str(record_a.get("anchor_status", "")),
            "notes": "the interval-first route did not reach this hypothesis",
        })
    for record_b in route_b.to_dict("records"):
        if str(record_b["hypothesis_id"]) in matched_b:
            continue
        rows.append({
            "comparison_id": f"cmp_{len(rows):03d}", "route_a_hypothesis_id": "",
            "route_b_hypothesis_id": str(record_b["hypothesis_id"]),
            "pathway_or_process": str(record_b.get("pathway_or_process", "")),
            "gene_family": str(record_b.get("gene_family", "")),
            "shared_interval_genes": ";".join(sorted(_genes(record_b.get("member_interval_genes", "")))),
            "match_basis": "", "match_tier": "",
            "relationship_tier": str(record_b.get("relationship_tier", "unclassified")),
            "route_a_source_backing": "", "route_b_source_backing": str(record_b.get("source_backing", "")),
            "convergence_status": "interval_only_supported", "rank_contribution": "none",
            "route_a_source_ids": "", "route_b_source_ids": "", "shared_source_ids": "",
            "source_independence_status": "not_applicable", "interval_overlap_status": "no_shared_interval_genes",
            "route_a_anchor_status": "",
            "notes": "new hypothesis reached only by reading the interval",
        })

    comparison = _frame(rows, COMPARISON_COLUMNS)
    backed_a = int(route_a.get("source_backing", pd.Series(dtype=str)).astype(str).ne("").sum()) if not route_a.empty else 0
    backed_b = int(route_b.get("source_backing", pd.Series(dtype=str)).astype(str).ne("").sum()) if not route_b.empty else 0
    matched = comparison[comparison.match_tier.ne("")] if not comparison.empty else comparison
    best_tier = _highest_tier(list(matched["match_tier"])) if not matched.empty else ""

    # Overall strictness: the strongest matched tier names the overall
    # state, and a gene-overlap-only match never becomes a strong convergence.
    if not len(route_a) and not len(route_b):
        overall = "both_weak"
    elif not backed_a and not backed_b:
        overall = "both_weak"
    elif best_tier == "tier1_mechanistic":
        overall = "strict_mechanistic_convergent"
    elif best_tier == "tier2_family":
        overall = "family_convergent_only"
    elif best_tier == "tier3_gene_overlap":
        overall = "gene_overlap_only"
    elif best_tier == "tier4_process_related":
        overall = "process_related"
    elif best_tier == "tier5_text_only":
        overall = "text_only_match"
    elif len(route_a) and len(route_b):
        overall = "discordant"
    elif len(route_a) or len(route_b):
        overall = "single_route_supported"
    else:
        overall = "both_weak"

    unanchored_a = int(route_a.get("anchor_status", pd.Series(dtype=str)).astype(str).eq("unanchored").sum()) if not route_a.empty else 0
    unanchored_b = int(route_b.get("anchor_status", pd.Series(dtype=str)).astype(str).eq("unanchored").sum()) if not route_b.empty else 0
    strong_statuses = ("mechanistic_convergence", "family_or_subfamily_convergence")
    convergent_backed = int(comparison[comparison.convergence_status.isin(strong_statuses)
                                       & comparison.route_a_source_backing.astype(str).ne("")].shape[0]) if not comparison.empty else 0
    strict_rows = int((comparison["convergence_status"] == "mechanistic_convergence").sum()) if not comparison.empty else 0
    matched_rows = int(len(matched))

    review_rows: list[dict[str, str]] = []
    if overall == "discordant":
        review_rows.append({"reason": "routes_discordant", "target": case_id or "case",
                            "required_next_action": "両routeが別のmechanismを指している。どちらも捨てず、識別できる実験を選ぶ",
                            "allowed_choices": "両方を保持して次実験で識別 | 追加文献検証 | 追加mapping"})
    if overall in ("gene_overlap_only", "process_related", "text_only_match"):
        review_rows.append({"reason": "weak_convergence_only", "target": case_id or "case",
                            "required_next_action": "一致がmechanism同一を保証しないtierで止まっている。stable process ID / family検証 / 識別実験でtierを上げるか、弱い一致として保持して次実験で判断する",
                            "allowed_choices": "process ID正規化の追加 | family/subfamily検証 | 識別実験 | 弱い一致として保持"})
    if overall == "both_weak":
        review_rows.append({"reason": "both_routes_weak", "target": case_id or "case",
                            "required_next_action": "source-backedな仮説がどちらのrouteからも出ていない。追加mappingやannotation rescueへ戻る",
                            "allowed_choices": "追加mapping | annotation rescue | regulatory/SV lane | abstain"})
    if unanchored_a or unanchored_b:
        review_rows.append({"reason": "unanchored_hypotheses_present", "target": f"A={unanchored_a} B={unanchored_b}",
                            "required_next_action": "区間内にメンバーを持たない仮説がある。捏造リスクを負っているので、区間外仮説として別扱いにするか区間の妥当性を見直す",
                            "allowed_choices": "区間外仮説として保持 | 区間・座標を見直す | 破棄"})
    reviews = _frame(review_rows, REVIEW_COLUMNS) if review_rows else pd.DataFrame(columns=REVIEW_COLUMNS)

    paths = {"dual_route_hypothesis_comparison": out / "dual_route_hypothesis_comparison.tsv",
             "dual_route_review_required": out / "dual_route_review_required.tsv"}
    comparison.to_csv(paths["dual_route_hypothesis_comparison"], sep="\t", index=False, lineterminator="\n")
    reviews.to_csv(paths["dual_route_review_required"], sep="\t", index=False, lineterminator="\n")

    handoff_rows = []
    for row in comparison.to_dict("records"):
        shared = str(row.get("shared_interval_genes", "")).strip()
        if not shared:
            continue
        cmp_id = str(row["comparison_id"])
        family_name = str(row.get("route_a_family", "") or row.get("route_b_family", "") or row.get("family_id", "")).strip()
        alias = alias_canonical(family_name) or ""
        mech_id = alias or (family_name if " " not in family_name and "(" not in family_name and ";" not in family_name else "")
        handoff_rows.append({
            "hypothesis_id": cmp_id,
            "family_id": mech_id,
            "gene_family": family_name,
            "member_interval_genes": shared,
        })
    if handoff_rows:
        handoff = pd.DataFrame(handoff_rows, columns=FAMILY_TO_COPY_INPUT_COLUMNS)
        paths["family_to_copy_input"] = out / "family_to_copy_input.tsv"
        handoff.to_csv(paths["family_to_copy_input"], sep="\t", index=False, lineterminator="\n")

    matched_rows_data = matched.to_dict("records")
    if "shared_source_chain" in {str(row["source_independence_status"]) for row in matched_rows_data}:
        source_independence = "shared_source_chain"
    elif matched_rows_data and all(str(row["source_independence_status"]) == "independent_source_chains"
                                   for row in matched_rows_data):
        source_independence = "independent_source_chains"
    elif matched_rows_data:
        source_independence = "source_ids_incomplete"
    else:
        source_independence = "not_applicable"
    matched_a_ids = sorted({token for row in matched_rows_data for token in str(row["route_a_source_ids"]).split(";") if token})
    matched_b_ids = sorted({token for row in matched_rows_data for token in str(row["route_b_source_ids"]).split(";") if token})
    matched_shared_ids = sorted({token for row in matched_rows_data for token in str(row["shared_source_ids"]).split(";") if token})
    overlap_status = "shared_interval_genes" if any(str(row["shared_interval_genes"]) for row in matched_rows_data) else "no_shared_interval_genes"

    convergence = {
        "schema_version": SCHEMA_VERSION, "case_id": case_id, "overall_status": overall,
        "route_a_hypotheses": int(len(route_a)), "route_b_hypotheses": int(len(route_b)),
        "route_a_source_backed": backed_a, "route_b_source_backed": backed_b,
        "benchmark_mode": benchmark_mode,
        "route_risk": {
            "route_a_failure_mode": "hallucinating a gene that is not in the interval",
            "route_b_failure_mode": "rationalising a gene that is in the interval",
            "route_a_unanchored": unanchored_a, "route_b_unanchored": unanchored_b,
            "strict_mechanistic_rows": strict_rows,
            "convergent_and_source_backed": convergent_backed,
            "why_two_routes": (
                "the routes fail in mirror-image ways, so a hypothesis reached independently by both has "
                "its existence checked by the interval and its relevance checked by the phenotype"
            ),
            "residual_risk": (
                "hypotheses that only one route reached still carry that route's failure mode; "
                "only strict mechanistic convergence lowers both mirror-image risks at once, and it is "
                "not proof and adds no score"
            ),
        },
        "status_counts": {str(k): int(v) for k, v in comparison["convergence_status"].value_counts().sort_index().items()} if not comparison.empty else {},
        "match_tier_counts": {str(k): int(v) for k, v in comparison["match_tier"].value_counts().sort_index().items()} if not comparison.empty else {},
        "relationship_tier_counts": {str(k): int(v) for k, v in comparison["relationship_tier"].value_counts().sort_index().items()} if not comparison.empty else {},
        "independence": {
            "design": "the two routes were generated separately; this comparison is the only place they meet",
            "route_a_source_ids": matched_a_ids,
            "route_b_source_ids": matched_b_ids,
            "shared_source_ids": matched_shared_ids,
            "source_independence_status": source_independence,
            "interval_overlap_status": overlap_status,
            "phenotype_seen_by_route_b": phenotype_seen_by_route_b,
        },
        "benchmark_rule": "interval-first results are interval-informed and must be scored in their own lane, never pooled with blind phenotype-only recall; B0 (phenotype-blind) and B1 (phenotype overlay) runs are recorded via phenotype_seen_by_route_b and aggregated separately",
        "safety": [
            "convergence is agreement between routes, not extra causal evidence, and adds no score",
            "a shared interval gene alone is gene_overlap_only, never strict mechanistic convergence",
            "free-text matches without a stable concept ID stay text_only_match and never upgrade to strict convergence",
            "a shared source chain is recorded so agreement on the same reference record is not presented as independent",
            "a discordant result keeps both mechanisms alive instead of picking one",
            "both_weak routes back to additional mapping / annotation rescue rather than forcing a candidate",
        ],
        "inputs": {
            "route_a": file_fingerprint(runtime_path(route_a_arg)),
            "route_b": file_fingerprint(runtime_path(route_b_arg)),
            "route_b_manifest": file_fingerprint(runtime_path(route_b_manifest_arg)) if route_b_manifest_arg else None,
        },
        "outputs": {name: str(path) for name, path in sorted(paths.items())},
    }
    convergence_path = out / "dual_route_convergence.json"
    convergence_path.write_text(json.dumps(convergence, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"convergence": str(convergence_path), "overall_status": overall,
            "matched_pairs": matched_rows, "comparison_rows": int(len(comparison)),
            "phenotype_seen_by_route_b": phenotype_seen_by_route_b,
            **{name: str(path) for name, path in paths.items()}}
from .route_vocabulary import ROUTE_COLUMNS, derive_route, phenotype_tokens
