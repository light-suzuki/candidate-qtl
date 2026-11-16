"""What a retrieved source actually supports, not merely that it exists (issue #226).

``source-verify`` closed the fabrication hole: a citation that does not resolve
in the literature database is rejected, because model memory is not a citable
source. But a resolvable citation is a weaker guarantee than it looks. A real
paper can mention a gene without supporting the causal statement it was cited
for, and these are not the same object:

* a knockout that reproduces the phenotype and a complementation that rescues it;
* an overexpression / RNAi phenotype;
* an expression correlation;
* a GWAS / QTL association;
* a review that repeats someone else's conclusion;
* a paper that names the gene and says nothing about the mechanism.

This module classifies which of those the retrieved record looks like, from the
title/abstract text that came back, and keeps two axes apart that the previous
single ``verification_status`` mixed together:

``citation_verified``
    Did the identifier resolve to a real record? A citation axis only.
``claim_support_status``
    Does the retrieved text support the claim it was cited for? An evidence axis.

Deliberate limits:

* **Cue matching is a triage, not a reading.** A cue lexicon can tell "this
  abstract describes a complementation experiment" from "this abstract reports a
  correlation"; it cannot tell whether the experiment supports *this* mechanism.
  Anything short of a clear cue lands in ``review_required``, never in a silent
  pass or a rejection.
* **Nothing here writes a score.** Classification decides what a reader must
  check, not what a candidate is worth.
* **An agent may extract, but not assert.** An extraction payload passes a
  deterministic gate: its quoted span must occur verbatim in the retrieved text,
  or the row is rejected with a reason. A quote that cannot be located is not
  evidence, whoever produced it.
"""

from __future__ import annotations

import re

# The evidence relation between the source and the claim it was cited for.
CLAIM_RELATIONS = (
    "direct_genetic_test",
    "perturbation_phenotype",
    "expression_association",
    "genetic_association",
    "review_or_background",
    "mention_without_mechanism",
    "no_claim_term_in_source",
    "not_evaluated_no_text",
)
# What that relation means for the claim. Only the first two are mechanistic
# support; the rest are explicitly weaker states, kept distinct so a report
# cannot flatten them into "sourced".
CLAIM_SUPPORT_STATUSES = (
    "supported_direct_genetic_test",
    "supported_perturbation_phenotype",
    "associative_support_only",
    "background_or_review_only",
    "unsupported_mention",
    "review_required",
    "not_evaluated",
)
RETRIEVAL_SCOPES = ("title_abstract", "title_only", "full_text_cached", "none")

RELATION_TO_SUPPORT = {
    "direct_genetic_test": "supported_direct_genetic_test",
    "perturbation_phenotype": "supported_perturbation_phenotype",
    "expression_association": "associative_support_only",
    "genetic_association": "associative_support_only",
    "review_or_background": "background_or_review_only",
    "mention_without_mechanism": "review_required",
    "no_claim_term_in_source": "unsupported_mention",
    "not_evaluated_no_text": "not_evaluated",
}

# Ordered: the strongest experiment described wins. Each entry is
# (relation, compiled pattern). Patterns stay conservative -- a cue that fires on
# ordinary prose would turn "review_required" into a false pass.
_CUES: tuple[tuple[str, re.Pattern[str]], ...] = (
    ("direct_genetic_test", re.compile(
        r"\b(knock[- ]?outs?|knockout|loss[- ]of[- ]function|null (?:allele|mutants?)|"
        r"complementation|complemented|rescued? the (?:mutant )?phenotype|"
        r"crispr|cas9|t-dna insertion|ems (?:mutagenesis|mutants?)|allelic series|"
        r"transposon insertion|gene editing|targeted mutagenesis)\b", re.IGNORECASE)),
    ("perturbation_phenotype", re.compile(
        r"\b(over[- ]?express(?:ion|ing|ed)|rna ?i|knock[- ]?down|silencing|vigs|"
        r"antisense|ectopic expression|dominant[- ]negative|amirna)\b", re.IGNORECASE)),
    ("genetic_association", re.compile(
        r"\b(genome[- ]wide association|gwas|quantitative trait loc(?:us|i)|\bqtl\b|"
        r"linkage (?:analysis|mapping)|marker[- ]trait association|associated with the (?:trait|phenotype))\b",
        re.IGNORECASE)),
    ("expression_association", re.compile(
        r"\b(expression (?:level|profil|pattern|analysis)|transcript (?:level|abundance)|"
        r"qrt[- ]?pcr|rt[- ]?qpcr|rna[- ]seq|transcriptome|up[- ]?regulated|down[- ]?regulated|"
        r"co[- ]?expressed)\b", re.IGNORECASE)),
    ("review_or_background", re.compile(
        r"\b(this review|we review|a review of|review(?:s)? (?:the|recent)|perspectives?|"
        r"overview of|advances in)\b", re.IGNORECASE)),
)


def classify_claim(
    text: str,
    claim_terms_found: list[str],
    claim_terms: list[str],
    retrieval_scope: str,
    *,
    cue_cache: dict[str, tuple[str, str, tuple[str, ...]]] | None = None,
) -> tuple[str, str, list[str]]:
    """Return (relation, claim_support_status, cues) for one retrieved record.

    ``claim_terms_found`` is what the existing term check already located; this
    function only decides what *kind* of statement the source makes about them.
    An optional invocation-local cache reuses only the source cue scan. The
    claim-specific and retrieval gates still run for every claim.
    """
    if retrieval_scope == "none" or not str(text).strip():
        return "not_evaluated_no_text", "not_evaluated", []
    if claim_terms and not claim_terms_found:
        # It resolves and it is a real paper -- but it never mentions the subject
        # it was cited for. That is a citation problem, not counterevidence.
        return "no_claim_term_in_source", "unsupported_mention", []
    source_text = str(text)
    if cue_cache is not None and source_text in cue_cache:
        relation, status, cached_cues = cue_cache[source_text]
        return relation, status, list(cached_cues)
    cues: list[str] = []
    relation = ""
    for name, pattern in _CUES:
        found = sorted({match.group(0).casefold() for match in pattern.finditer(source_text)})
        if found:
            cues.extend(found)
            relation = relation or name
    if not relation:
        # Terms present, no experiment described in the retrievable text. A human
        # has to read it; the workflow will not guess in either direction.
        relation, status = "mention_without_mechanism", "review_required"
    else:
        status = RELATION_TO_SUPPORT[relation]
    if cue_cache is not None:
        cue_cache[source_text] = relation, status, tuple(cues)
    return relation, status, cues


def scope_for(record: dict[str, object], full_text_cached: bool = False) -> str:
    """How much of the source the classification actually saw."""
    if full_text_cached:
        return "full_text_cached"
    if str(record.get("abstractText", "")).strip():
        return "title_abstract"
    if str(record.get("title", "")).strip():
        return "title_only"
    return "none"


def classifiable_text(record: dict[str, object]) -> str:
    return " ".join(str(record.get(key, "")) for key in
                    ("title", "abstractText", "journalTitle", "keywordList", "fullText"))


# --------------------------------------------------------------------------- #
# agent extraction gate
# --------------------------------------------------------------------------- #

EXTRACTION_COLUMNS = ["citation", "claim_type", "evidence_relation", "quoted_support", "locator"]
EXTRACTION_STATUSES = ("accepted", "rejected_unknown_relation", "rejected_quote_not_in_retrieved_text",
                       "rejected_unknown_citation", "rejected_missing_field")


def _normalise(value: object) -> str:
    return re.sub(r"\s+", " ", str(value or "")).strip().casefold()


def gate_extraction(
    rows: list[dict[str, object]],
    retrieved_text: dict[str, str],
    hypothesis_citations: dict[str, set[str]] | None = None,
) -> list[dict[str, object]]:
    """Validate agent-extracted claim rows against the text that was retrieved.

    The gate is deterministic and unforgiving in one specific way: a quoted span
    that does not occur in the retrieved record is rejected. An extraction that
    passes annotates the row -- it never creates evidence of its own, and it can
    never promote a citation that did not resolve.

    When a hypothesis/citation roster is supplied, explicit targets must belong
    to it. Conflicting duplicates and scoped/unscoped overlaps fail closed;
    identical rows are deduplicated. No target is inferred from quoted text.
    """
    gated: list[dict[str, object]] = []
    normalised_sources: dict[str, str] = {}
    if hypothesis_citations is not None:
        unique: dict[tuple[str, str], dict[str, object]] = {}
        citation_scopes: dict[str, set[bool]] = {}
        for row in rows:
            citation = str(row.get("citation", "")).strip()
            target = str(row.get("hypothesis_id", "")).strip()
            if target:
                if target not in hypothesis_citations:
                    raise ValueError(f"claim extraction has unknown hypothesis_id: {target}")
                if citation not in hypothesis_citations[target]:
                    raise ValueError(f"claim extraction target/citation mismatch: {target}: {citation}")
            canonical = {column: str(row.get(column, "")).strip() for column in EXTRACTION_COLUMNS}
            if target:
                canonical["hypothesis_id"] = target
            key = (target, citation)
            if key in unique and unique[key] != canonical:
                raise ValueError(f"conflicting claim extractions for hypothesis/citation: {key}")
            unique[key] = canonical
            citation_scopes.setdefault(citation, set()).add(bool(target))
        for citation, scopes in citation_scopes.items():
            if len(scopes) > 1:
                raise ValueError(f"ambiguous scoped/unscoped claim extractions: {citation}")
        rows = list(unique.values())
    for row in rows:
        citation = str(row.get("citation", "")).strip()
        relation = str(row.get("evidence_relation", "")).strip()
        quote = str(row.get("quoted_support", "")).strip()
        missing = [column for column in ("citation", "evidence_relation", "quoted_support")
                   if not str(row.get(column, "")).strip()]
        if missing:
            status = "rejected_missing_field"
        elif citation not in retrieved_text:
            status = "rejected_unknown_citation"
        elif relation not in CLAIM_RELATIONS:
            status = "rejected_unknown_relation"
        else:
            if citation not in normalised_sources:
                normalised_sources[citation] = _normalise(retrieved_text[citation])
            status = ("accepted" if _normalise(quote) in normalised_sources[citation]
                      else "rejected_quote_not_in_retrieved_text")
        gated.append({
            "citation": citation, "claim_type": str(row.get("claim_type", "")).strip(),
            "evidence_relation": relation, "quoted_support": quote,
            "locator": str(row.get("locator", "")).strip(),
            "extraction_status": status,
            "claim_support_status": RELATION_TO_SUPPORT.get(relation, "review_required")
            if status == "accepted" else "review_required",
        })
        target = str(row.get("hypothesis_id", "")).strip()
        if target:
            gated[-1]["hypothesis_id"] = target
    return gated
