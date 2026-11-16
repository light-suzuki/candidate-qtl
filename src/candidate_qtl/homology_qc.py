"""Shared homology quality policy for discover and the evidence pipeline.

Both the simple `discover` entry point and the advanced evidence-build path must
grade a protein hit the same way: same identity/coverage/e-value thresholds,
same near-tie paralog rule, same quality classes. Centralizing the policy here
keeps a hit from meaning one thing in `discover` and another in `evidence-build`.
"""

from __future__ import annotations

import math
import re

# Quality gate (a hit worth reporting at all).
MIN_IDENTITY = 30.0
MIN_QUERY_COVERAGE = 0.50
MAX_EVALUE = 1e-5

# High-confidence, full-length thresholds.
FULL_QUERY_COVERAGE = 0.70
FULL_SUBJECT_COVERAGE = 0.70
FULL_IDENTITY = 50.0
FULL_EVALUE = 1e-10

# Medium-confidence subject-coverage floor.
MED_SUBJECT_COVERAGE = 0.30

# Two hits whose bitscores are this close cannot be called one clean ortholog.
PARALOG_BITSCORE_RATIO = 0.90


def _num(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def passes_quality(identity: object, query_coverage: object, evalue: object) -> bool:
    """The minimum bar for a hit to count as quality-gated homology."""
    qcov = _num(query_coverage)
    return _num(identity) >= MIN_IDENTITY and (math.isnan(qcov) or qcov >= MIN_QUERY_COVERAGE) and _num(evalue) <= MAX_EVALUE


def confidence(identity: object, query_coverage: object, subject_coverage: object, evalue: object) -> str:
    """Grade a best hit high / medium / low (shared with _parse_blast)."""
    ident, qcov, scov, ev = _num(identity), _num(query_coverage), _num(subject_coverage), _num(evalue)
    if not math.isnan(qcov) and qcov >= FULL_QUERY_COVERAGE and scov >= FULL_SUBJECT_COVERAGE and ident >= FULL_IDENTITY and ev <= FULL_EVALUE:
        return "high"
    if (math.isnan(qcov) or qcov >= MIN_QUERY_COVERAGE) and (math.isnan(scov) or scov >= MED_SUBJECT_COVERAGE) and ident >= MIN_IDENTITY and ev <= MAX_EVALUE:
        return "medium"
    return "low"


def quality_class(identity: object, query_coverage: object, subject_coverage: object) -> str:
    """Distinguish a full-length match from a domain-only or partial hit."""
    ident, qcov, scov = _num(identity), _num(query_coverage), _num(subject_coverage)
    if math.isnan(qcov) or math.isnan(scov):
        return "coverage_unknown"
    if qcov >= FULL_QUERY_COVERAGE and scov >= FULL_SUBJECT_COVERAGE:
        return "full_length_match"
    if ident >= FULL_IDENTITY and (qcov < FULL_QUERY_COVERAGE or scov < FULL_SUBJECT_COVERAGE):
        return "conserved_domain_only"
    return "short_or_partial_match"


def is_paralog_ambiguous(best_bitscore: object, second_bitscore: object | None) -> bool:
    """True when the second hit is within PARALOG_BITSCORE_RATIO of the best."""
    best = _num(best_bitscore)
    if second_bitscore is None or best <= 0:
        return False
    return _num(second_bitscore) >= PARALOG_BITSCORE_RATIO * best


def bitscore_ratio(best_bitscore: object, second_bitscore: object) -> float:
    best = _num(best_bitscore)
    return _num(second_bitscore) / best if best > 0 else 0.0


# --------------------------------------------------------------------------- #
# Multi-species annotation parsing and concordance handling
# --------------------------------------------------------------------------- #

# Mandatory disclaimer for homology-only candidate rationales.
HOMOLOGY_DISCLAIMER = "※ 全長相同性のみで判定。活性部位残基の保存は未確認。"

# Headline-hit selection disclosure. The closest hit is the
# best-scoring qualifying hit; the most informative hit is the first
# qualifying hit whose annotation states a function. When the two differ,
# both are shown with these labels so the reader can tell which rule chose
# each line.
HEADLINE_LABEL_CLOSEST = "最近縁"
HEADLINE_LABEL_INFORMATIVE = "機能記載"
NO_FUNCTION_NOTE = "（機能記載なし）"

# Concordance evaluation states.
CONCORDANCE_CONCORDANT = "concordant"
CONCORDANCE_DISCORDANT = "discordant"
CONCORDANCE_SINGLE_INFORMATIVE = "single_informative"
CONCORDANCE_ALL_UNINFORMATIVE = "all_uninformative"
CONCORDANCE_STATES = (
    CONCORDANCE_CONCORDANT,
    CONCORDANCE_DISCORDANT,
    CONCORDANCE_SINGLE_INFORMATIVE,
    CONCORDANCE_ALL_UNINFORMATIVE,
)

# Uninformative annotation patterns (case-insensitive).
UNINFORMATIVE_TERMS = (
    "uncharacterized",
    "hypothetical",
    "unknown protein",
    "protein of unknown function",
    "unnamed protein",
    "predicted protein",
    "expressed protein",
    "unidentified",
    "unclassified",
    "cDNA clone",
)

# Functional class words: a protein name containing one of these is treated as
# an informative annotation. A name with none of these (e.g. "WD40-1 protein")
# does not state a function and is treated as uninformative.
FUNCTIONAL_CLASS_WORDS = (
    "synthase",
    "synthetase",
    "reductase",
    "oxidase",
    "oxygenase",
    "hydroxylase",
    "monooxygenase",
    "dioxygenase",
    "isomerase",
    "kinase",
    "phosphatase",
    "hydrolase",
    "ligase",
    "lyase",
    "transferase",
    "polymerase",
    "dehydrogenase",
    "decarboxylase",
    "desaturase",
    "cyclase",
    "mutase",
    "epimerase",
    "dehydratase",
    "transporter",
    "translocase",
    "channel",
    "receptor",
    "factor",
    "pump",
    "peptidase",
    "protease",
    "peroxidase",
    "elongase",
    "methyltransferase",
    "acyltransferase",
)


def extract_species(text: object) -> str:
    """Extract species name from a UniProt / BLAST description or defline."""
    if not text:
        return ""
    s = str(text).strip()
    match = re.search(r"\bOS=([^=]+?)(?=\s+[A-Z]{2}=|$)", s)
    if match:
        return match.group(1).strip()
    return ""


def is_reviewed_entry(subject: object) -> bool:
    """True for Swiss-Prot (reviewed) entries (sp|...), False for TrEMBL/other."""
    if not subject:
        return False
    s = str(subject).strip()
    return s.startswith("sp|") or "|sp|" in s or bool(re.match(r"^sp\|", s, re.IGNORECASE))


def extract_accession(subject: object) -> str:
    """Extract bare protein accession from pipe-delimited or prefixed string."""
    if not subject:
        return ""
    s = str(subject).strip()
    # Match pipe forms: sp|P13114|CHSY_ARATH, tr|A0A123|..., or generic |ACC|
    pipe_match = re.match(r"^[a-z]{2}\|([^|\s]+)(?:\|.*)?$", s, re.IGNORECASE)
    if pipe_match:
        return pipe_match.group(1).strip()
    # Match database prefixes: UniProt:P13114, TAIR:AT1G69120, etc.
    db_match = re.match(r"^(?:uniprot|refseq|ncbi|tair|embl|ddbj|ensembl):([^|\s]+)$", s, re.IGNORECASE)
    if db_match:
        return db_match.group(1).strip()
    # Fallback to bare token before whitespace
    return s.split()[0] if s else ""


def extract_protein_name(description: object, subject: object = "") -> str:
    """Extract clean protein functional name from BLAST description / defline."""
    text = str(description or "").strip()
    subj = str(subject or "").strip()
    if not text and subj:
        text = subj
    if subj and text.startswith(subj):
        text = text[len(subj):].strip()
    # Remove leading db|acc|name prefix if still present
    text = re.sub(r"^[a-z]{2}\|[^|\s]+\|[^|\s]*\s*", "", text, flags=re.IGNORECASE).strip()
    # Cut before UniProt annotation tags: OS=, OX=, GN=, PE=, SV=
    cut = re.search(r"\s+(?:OS|OX|GN|PE|SV)=", text)
    if cut:
        text = text[:cut.start()].strip()
    return text


def is_informative_annotation(text: object, subject: object = "") -> bool:
    """True if annotation states a specific function rather than a placeholder.

    A protein name is informative only when it contains a functional class word
    (e.g. synthase/reductase/isomerase) or comes from a reviewed Swiss-Prot
    entry. Names like "WD40-1 protein" or "CASTOR" that name a protein without
    stating a function are treated as uninformative (#422).
    """
    if not text:
        return False
    s = str(text).strip().casefold()
    if not s or s == "nan" or s == "none":
        return False
    for term in UNINFORMATIVE_TERMS:
        if term.casefold() in s:
            return False
    if any(word in s for word in FUNCTIONAL_CLASS_WORDS):
        return True
    # Reviewed Swiss-Prot entries carry a curated name even when the plain
    # protein name alone has no functional class word (e.g. TTG1_ARATH).
    if subject and is_reviewed_entry(subject):
        return True
    # If it's just a locus tag or accession
    if re.fullmatch(r"^[a-z0-9_.-]+$", s) and any(c.isdigit() for c in s) and len(s) <= 15:
        # e.g. At1g01010 or P12345 alone without descriptive words
        if not re.search(r"(?:synthase|kinase|transporter|protein|factor|reductase|oxidase|isomerase|ligase|hydrolase|domain)", s):
            return False
    return False


def _normalize_function_signature(protein_name: str) -> set[str]:
    """Normalize protein name to functional token set for concordance comparison."""
    lowered = protein_name.casefold()
    # Remove common filler tokens
    stopwords = {"protein", "putative", "probable", "like", "family", "homolog",
                 "chloroplastic", "cytosolic", "mitochondrial", "nuclear", "membrane",
                 "isoform", "subunit", "chain", "type", "at", "os", "the", "and", "of", "in"}
    tokens = {t for t in re.findall(r"[a-z0-9]+", lowered) if len(t) >= 3 and t not in stopwords}
    return tokens


def _functional_class_set(protein_name: str) -> set[str]:
    """Functional class words found in a protein name (e.g. {'hydroxylase'})."""
    lowered = protein_name.casefold()
    return {word for word in FUNCTIONAL_CLASS_WORDS if word in lowered}


def assess_annotation_concordance(hits: list[dict[str, object]]) -> tuple[str, str]:
    """Assess functional concordance across multiple homology hits.

    Returns:
        (concordance_state, summary_details)
        where concordance_state is one of:
        - 'concordant'
        - 'discordant'
        - 'single_informative'
        - 'all_uninformative'
    """
    if not hits:
        return CONCORDANCE_ALL_UNINFORMATIVE, "no homology hits"

    # Filter to quality passing hits
    qualifying_hits = [
        hit for hit in hits
        if passes_quality(hit.get("identity"), hit.get("query_coverage"), hit.get("evalue"))
    ]
    if not qualifying_hits:
        return CONCORDANCE_ALL_UNINFORMATIVE, "no quality-passing homology hits"

    # Extract informative annotations
    informative_entries: list[tuple[str, str, dict[str, object]]] = []
    for hit in qualifying_hits:
        pname = extract_protein_name(hit.get("description", ""), hit.get("subject", ""))
        species = extract_species(hit.get("description", "") or hit.get("subject", ""))
        if is_informative_annotation(pname, hit.get("subject", "")):
            informative_entries.append((pname, species, hit))

    if not informative_entries:
        return CONCORDANCE_ALL_UNINFORMATIVE, f"{len(qualifying_hits)} qualifying hits, all uninformative/uncharacterized"

    if len(informative_entries) == 1:
        pname, species, _ = informative_entries[0]
        return CONCORDANCE_SINGLE_INFORMATIVE, f"1 informative hit: '{pname}' ({species or 'unknown species'})"

    # Multiple informative hits: compare functional class words first. Two
    # different enzyme classes (e.g. hydroxylase vs monooxygenase) are a
    # discordance even when the names share a generic token such as
    # "flavonoid".
    class_sets = [_functional_class_set(entry[0]) for entry in informative_entries]
    nonempty_class_sets = [cs for cs in class_sets if cs]
    if len(nonempty_class_sets) >= 2:
        base = nonempty_class_sets[0]
        if any(base.isdisjoint(cs) for cs in nonempty_class_sets[1:]):
            species_summary = ", ".join(f"{sp or 'unknown'}:{nm}" for nm, sp, _ in informative_entries[:3])
            return CONCORDANCE_DISCORDANT, f"discordant annotations across {len(informative_entries)} hits ({species_summary})"

    # Fall back to shared-token comparison when class words are absent or overlap.
    names = [entry[0] for entry in informative_entries]
    sig0 = _normalize_function_signature(names[0])

    is_concordant = True
    for name in names[1:]:
        sig_i = _normalize_function_signature(name)
        # Check token intersection or substring containment
        if sig0 and sig_i:
            overlap = sig0 & sig_i
            # If no functional tokens overlap and neither name is a substring of the other
            if not overlap and (name.casefold() not in names[0].casefold() and names[0].casefold() not in name.casefold()):
                is_concordant = False
                break
        elif name.casefold() != names[0].casefold():
            is_concordant = False
            break

    species_summary = ", ".join(f"{sp or 'unknown'}:{nm}" for nm, sp, _ in informative_entries[:3])
    if is_concordant:
        return CONCORDANCE_CONCORDANT, f"{len(informative_entries)} concordant informative hits ({species_summary})"
    else:
        # Note: Invariant 3: Do NOT resolve discordance by taxonomy proximity.
        return CONCORDANCE_DISCORDANT, f"discordant annotations across {len(informative_entries)} hits ({species_summary})"


def _hit_rationale_line(
    hit: dict[str, object],
    *,
    label: str = "",
    mark_uninformative: bool = False,
) -> str:
    """One rationale line for a single hit row, grounded in that row only.

    The line carries species, accession, identity and both coverages from the
    hit row plus its annotation name. When `mark_uninformative` is set, a note
    states that the annotation does not describe a function (issue #424).
    """
    pname = extract_protein_name(hit.get("description", ""), hit.get("subject", ""))
    if not pname:
        pname = "(no protein name in hit row)"
    species = extract_species(hit.get("description", "") or hit.get("subject", ""))
    acc = extract_accession(hit.get("subject", ""))
    target_parts = []
    if species:
        target_parts.append(species)
    if acc:
        target_parts.append(f"({acc})")
    target = " ".join(target_parts) if target_parts else str(hit.get("subject", "homolog"))
    identity = _num(hit.get("identity"))
    qcov = _num(hit.get("query_coverage"))
    scov = _num(hit.get("subject_coverage"))
    ident_text = f"{identity:.1f}%" if not math.isnan(identity) else "identity n/a"
    qcov_text = f"{qcov:.2f}" if not math.isnan(qcov) else "n/a"
    scov_text = f"{scov:.2f}" if not math.isnan(scov) else "n/a"
    line = (
        f"{target}, identity {ident_text} / "
        f"query coverage {qcov_text} / subject coverage {scov_text}, "
        f"annotated as {pname}"
    )
    if mark_uninformative:
        line = f"{line}{NO_FUNCTION_NOTE}"
    return f"{label}: {line}" if label else line


def format_grounded_rationale(
    subject_or_species: object = "",
    protein_name: object = "",
    accession: object = "",
    *,
    hits: list[dict[str, object]] | None = None,
    include_disclaimer: bool = True,
) -> str:
    """Format a candidate homology rationale grounded strictly in actual hit rows.

    Enforces strict phrasing: 'homologous to X, annotated as Y' with identity and
    coverage from the actual hit row, plus a mandatory disclaimer:
    '※ 全長相同性のみで判定。活性部位残基の保存は未確認。'
    A discordant hit set appends the best other informative hit as a note; the
    note never resolves the discordance by taxonomy (#422).

    The headline-hit selection rule is disclosed (issue #424): the closest hit
    is the best-scoring qualifying hit, the most informative hit is the first
    qualifying hit whose annotation states a function. When the two differ,
    both are shown with their selection-rule labels (最近縁 / 機能記載).
    `hits` is expected in bitscore-descending order (best first); callers sort
    with `_sorted_top` (discover.py).
    """
    if hits is not None:
        qualifying = [
            h for h in hits
            if passes_quality(h.get("identity"), h.get("query_coverage"), h.get("evalue"))
        ]
        if not qualifying:
            rationale = "no qualifying homology hits"
            if include_disclaimer:
                rationale = f"{rationale}. {HOMOLOGY_DISCLAIMER}"
            return rationale

        closest = qualifying[0]
        informative = []
        for h in qualifying:
            pname = extract_protein_name(h.get("description", ""), h.get("subject", ""))
            if is_informative_annotation(pname, h.get("subject", "")):
                informative.append((pname, h))
        if informative and closest is not informative[0][1]:
            # The closest hit carries no functional annotation and a different
            # hit does: disclose both with their selection-rule labels.
            rationale = (
                f"{_hit_rationale_line(closest, label=HEADLINE_LABEL_CLOSEST, mark_uninformative=True)}"
                f"<br>{_hit_rationale_line(informative[0][1], label=HEADLINE_LABEL_INFORMATIVE)}"
            )
            chosen = informative[0][1]
            pname = informative[0][0]
        else:
            if informative:
                chosen_pname, chosen = informative[0]
            else:
                chosen = closest
                chosen_pname = extract_protein_name(chosen.get("description", ""), chosen.get("subject", ""))
            pname = chosen_pname
            rationale = "homologous to " + _hit_rationale_line(chosen)

        # When multiple informative hits disagree, note a hit whose functional
        # class conflicts with the chosen one (falling back to the best other
        # informative hit), without resolving the disagreement.
        state, _ = assess_annotation_concordance(qualifying)
        if state == CONCORDANCE_DISCORDANT and len(informative) > 1:
            chosen_classes = _functional_class_set(pname)
            other = None
            if chosen_classes:
                other = next((entry for entry in informative if entry[1] is not chosen and chosen_classes.isdisjoint(_functional_class_set(entry[0]))), None)
            if other is None:
                other = informative[1] if informative[1][1] is not chosen else informative[0]
            other_pname, other_hit = other
            other_species = extract_species(other_hit.get("description", "") or other_hit.get("subject", ""))
            other_acc = extract_accession(other_hit.get("subject", ""))
            other_target = " ".join(part for part in [other_species, f"({other_acc})" if other_acc else ""] if part) or "another hit"
            rationale = f"{rationale}. ※ 他のヒットの注釈は {other_target} の {other_pname} で不一致"
    else:
        target_str = str(subject_or_species or "").strip()
        acc_str = str(accession or "").strip()
        if acc_str and acc_str not in target_str:
            target = f"{target_str} ({acc_str})" if target_str else acc_str
        else:
            target = target_str or "reference protein"

        pname = str(protein_name or "").strip() or "uncharacterized protein"
        rationale = f"homologous to {target}, annotated as {pname}"

    if include_disclaimer:
        rationale = f"{rationale}. {HOMOLOGY_DISCLAIMER}"
    return rationale
