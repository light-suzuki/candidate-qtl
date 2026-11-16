"""Shared identifier grammar for accession / locus-tag / symbol decisions.

The workflow has to decide in several places whether a token is a real gene
symbol or merely an identifier: family-index symbol extraction (``GN=`` fields
and UniProt deflines), the citable-source grammar in ``hardening``, the Route B
semantic gate, and benchmark record parsing.  Every consumer previously kept
its own copy of the grammar, which let it drift (issue #297): the legume
family-symbol fixes (2addf05 / afa33de) implemented a 6/10-character UniProt
pattern that missed ``F8T7T1`` and misjudged ``A0ACM8R380`` as a usable symbol.

This module is the single source of truth.  Consumers import from here instead
of re-implementing accession grammars.

Decision rules are deliberately conservative:

- a token that matches an identifier grammar is NEVER a usable gene symbol
- locus-tag detection is prefix-table driven and case-sensitive, so real
  symbols such as ``SBE1``, ``GAI``, ``CHI`` or ``FTa2`` are never mistaken
  for tags
- free natural-language protein names are never auto-symbolized here; a
  protein name only becomes a family candidate when it is a symbol-like token
  itself (e.g. ``FTa2``) or when a curated alias registry says so (see
  ``symbol_candidates`` with ``alias_lookup``)
"""

from __future__ import annotations

import re

# Canonical UniProt accession grammar (6 or 10 characters).  Mirrors
# hardening's citable-source check; a consistency test pins the two copies
# together so they cannot drift unnoticed.
#
# The 6-character form is covered by ``[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}``
# with one block (e.g. the real Medicago entry ``F8T7T1``), so accessions are
# not restricted to the historical ``[OPQ]`` prefix.
UNIPROT_ACCESSION_PATTERN = (
    r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]"
    r"|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"
)
UNIPROT_ACCESSION = re.compile(rf"^(?:{UNIPROT_ACCESSION_PATTERN})$")

# Arabidopsis AGI locus (TAIR-style).
AGI_PATTERN = r"AT[1-5CM]G[0-9]{5}(?:\.[0-9]+)?"
AGI = re.compile(rf"^(?:{AGI_PATTERN})$")

# UniProt defline id, e.g. ``tr|A0ACM8R380|A0ACM8R380_MEDTR`` or
# ``sp|Q9LZS3|Q9LZS3_ARATH``.  Group 1 = accession, group 2 = entry name.
UNIPROT_ID = re.compile(r"^[a-z]{2}\|([^|]+)\|([^|]+)$", re.IGNORECASE)

# Unified ``GN=`` field extraction.  A gene-name field in a UniProt defline is
# a single non-whitespace, non-``;`` token.
GN_FIELD = re.compile(r"(?:^|\s)GN=([^\s;]+)", re.IGNORECASE)

# Curated, case-sensitive locus-tag prefixes followed by a digit-bearing id
# (e.g. ``MTR_7g085020``, ``LotjaGi1g1v0638100.1``, ``Os01g01010``).  The
# case-sensitive table is what protects symbols such as ``SBE1`` (no prefix)
# or ``GAI`` (``GA`` is not a prefix) from being misclassified.  Additions are
# curated; do not loosen this table ad hoc.
LOCUS_TAG_PREFIXES: tuple[str, ...] = (
    "MTR_", "MTR.", "Medtr",          # Medicago truncatula
    "LotjaGi", "Lj",                  # Lotus japonicus
    "GLYMA_",                         # Glycine max
    "LOC_",                           # generic LOC_Os... style
    "Os", "Sb", "Zm", "Hv", "Ta",     # rice / sorghum / maize / barley / wheat
    "Solyc", "Sotub", "Potri", "Prupe", "Phvul", "Vigun", "Vigan",
    "Manes", "Bradi", "Seita", "HORVU", "Traes", "Aqcoe", "FvH4",
    "Ca_", "Csa", "Cucsa", "Bvu", "Tu", "Brara", "Bol", "Bna", "Thalv",
)
# Remainder after a known prefix: at least one digit, optionally preceded by
# letters (``LOC_Os01g01010``) and optionally followed by a version suffix.
_LOCUS_TAG_REMAINDER = re.compile(r"^[A-Za-z]*\d[A-Za-z0-9]*(?:\.\d+)?$")

# UniProt defline trailing annotations that delimit the protein name field.
_ANNOTATION_FIELD = re.compile(r"\s(?:OS|GN|PE|SV)=")

# Remainder rules after a known locus-tag prefix.  The generic rule requires a
# leading digit (``MTR_7g085020``, ``Os01g01010``, ``LotjaGi1g1v0638100.1``),
# so a real symbol with a species prefix (``OsSWEET1``, ``LjNFP2``, ``TaGW2``)
# is never misread as a locus tag.  A few prefix families are exceptions whose
# own convention puts letters before the digits (``LOC_Os01g01010``,
# ``TraesCS5A02G401300``, ``BraraA01001``, ``BnaA01g12340D``, ``BvuBv9_000010``,
# ``TuG000010``); for those, a letter-leading remainder is a locus tag.
_LOCUS_TAG_DIGIT_REMAINDER = re.compile(r"^\d[A-Za-z0-9]*(?:\.\d+)?$")
_LOCUS_TAG_LETTER_REMAINDER = re.compile(r"^[A-Za-z]*\d[A-Za-z0-9_]*(?:\.\d+)?$")
LETTER_FIRST_LOCUS_PREFIXES: frozenset[str] = frozenset({
    "LOC_", "Traes", "Brara", "Bna", "Bvu", "Tu",
})

# Issue #344: the ``LotjaGi`` (JGI/LjG) family has its own remainder rule
# because the pinned legume proteome (frozen DB ``legume_plant_pep3``) contains
# identifier-only deflines in two shapes that the digit-leading rule above
# rejects (45,662 rows total, 4,574 unique ids, both observed in real data):
#
# - an ``_LC.<version>`` low-confidence gene-model suffix on the regular
#   chromosome/gene form, e.g. ``LotjaGi2g1v0260700_LC.1`` (4,520 unique ids)
# - letter-leading assembly-variant ids ``LotjaGiC1v...`` / ``LotjaGiM1v...``
#   (``C`` = chromosome, ``M`` = mitochondrion), e.g.
#   ``LotjaGiC1v0000600.1`` / ``LotjaGiM1v0000600.1`` (54 unique ids)
#
# Both branches are scoped deliberately: the ``_LC.`` suffix requires its
# version and only follows the digit-leading form, and the letter-leading
# branch is limited to the observed ``C1v`` / ``M1v`` markers, so the grammar
# is not widened into an unbounded identifier pattern.
_LOTJA_REMAINDER = re.compile(
    r"^(?:[CM]1v\d+(?:\.\d+)?|\d[A-Za-z0-9]*(?:_LC\.\d+)?(?:\.\d+)?)$"
)

_PIPE_ACCESSION = re.compile(r"^[a-z]{2}\|([^|]+)(?:\|.*)?$", re.IGNORECASE)
_DB_PREFIX = re.compile(
    r"^(?:uniprot|uniprotkb|tair|refseq|ncbi|genbank|ensembl|ensemblplants|embl|ddbj|phytozome|gramene|medtr|medicago|lotus):",
    re.IGNORECASE,
)


def canonical_accession(token: object) -> str:
    """Extract canonical bare accession or identifier from common database forms.

    Normalizes:
    - ``sp|P13114|CHSY_ARATH`` -> ``P13114``
    - ``sp|P13114|`` -> ``P13114``
    - ``tr|A0A1I9LSP1|A0A1I9LSP1_ARATH`` -> ``A0A1I9LSP1``
    - ``UniProt:P13114`` / ``UniProtKB:P13114`` -> ``P13114``
    - ``TAIR:AT1G69120`` -> ``AT1G69120``
    - ``P13114`` -> ``P13114``
    - ``AT1G69120`` -> ``AT1G69120``
    """
    if token is None:
        return ""
    text = str(token).strip()
    if not text:
        return ""
    pipe_match = _PIPE_ACCESSION.match(text)
    if pipe_match:
        return pipe_match.group(1).strip()
    stripped = _DB_PREFIX.sub("", text)
    return stripped.strip()


def is_uniprot_accession(token: object) -> bool:
    """True when ``token`` is a canonical 6- or 10-character UniProt accession."""
    return bool(isinstance(token, str) and token and UNIPROT_ACCESSION.fullmatch(token))


def is_agi_locus(token: object) -> bool:
    """True when ``token`` is an Arabidopsis AGI locus id (e.g. AT1G01010)."""
    return bool(isinstance(token, str) and token and AGI.fullmatch(token))


def is_locus_tag(token: object) -> bool:
    """True when ``token`` is a known plant locus-tag id.

    The remainder after a known case-sensitive prefix must start with a digit,
    except for the letter-first prefix families (``LOC_Os01g01010``,
    ``TraesCS5A02G401300``, ``BraraA01001``, ``BnaA01g12340D``, ``BvuBv9_000010``,
    ``TuG000010``), so ``MTR_7g085020`` / ``LotjaGi1g1v0638100.1`` /
    ``Os01g01010`` qualify while symbols such as ``SBE1``, ``FTa2`` or
    species-prefixed real symbols (``OsSWEET1``, ``LjNFP2``) never do.

    The ``LotjaGi`` family additionally accepts the two identifier shapes
    observed in real data (issue #344): the ``_LC.<version>`` low-confidence
    gene-model suffix (``LotjaGi2g1v0260700_LC.1``) and the letter-leading
    assembly-variant ids ``LotjaGiC1v...`` / ``LotjaGiM1v...``
    (``LotjaGiC1v0000600.1``).  Unobserved shapes stay rejected.
    """
    if not isinstance(token, str) or not token:
        return False
    for prefix in LOCUS_TAG_PREFIXES:
        if token.startswith(prefix):
            remainder = token[len(prefix):]
            if prefix == "LotjaGi":
                return bool(_LOTJA_REMAINDER.fullmatch(remainder))
            pattern = (_LOCUS_TAG_LETTER_REMAINDER if prefix in LETTER_FIRST_LOCUS_PREFIXES
                       else _LOCUS_TAG_DIGIT_REMAINDER)
            return bool(pattern.fullmatch(remainder))
    return False


def gene_name_field(text: object) -> str | None:
    """Return the ``GN=`` value from a UniProt-style defline, or None.

    Old-style deflines carry ``GN=Name=FT`` / ``GN=OrderedLocusNames=...``; a
    value that itself contains ``=`` is an attribute wrapper, not a symbol, so
    it is rejected and the caller falls back to the protein-name path.
    """
    if not text:
        return None
    match = GN_FIELD.search(str(text))
    if not match:
        return None
    value = match.group(1)
    return value if "=" not in value else None


def uniprot_name_field(subject_id: str, subject_title: str) -> str | None:
    """Return the protein-name field of a UniProt defline, or None.

    Only UniProt-shaped subject ids are decoded (``xx|ACC|ENTRY_NAME``).  The
    returned name is the text before the first ``OS=`` / ``GN=`` / ``PE=`` /
    ``SV=`` annotation, e.g. ``Flowering locus protein T`` or ``FTa2``.
    Free-text titles of non-UniProt subjects are not interpreted.
    """
    match = UNIPROT_ID.match(subject_id)
    if not subject_title or not match:
        return None
    title = str(subject_title).strip()
    prefix = match.group(0)
    if title.startswith(prefix):
        title = title[len(prefix):].lstrip()
    cut = _ANNOTATION_FIELD.search(title)
    if cut:
        title = title[: cut.start()]
    name = title.strip()
    return name or None


def classify_gene_name_field(value: object) -> str:
    """Classify a single token from a ``GN=``-style field.

    Returns one of: ``gene_symbol`` / ``locus_tag`` / ``uniprot_accession`` /
    ``agi_locus`` / ``empty``.  ``gene_symbol`` is the default for any
    non-empty token that matches no identifier grammar (conservative: an
    identifier that slips past the grammar still looks like a symbol, so the
    caller may additionally carry provenance columns for audit).
    """
    if value is None or not str(value).strip():
        return "empty"
    token = str(value).strip()
    if is_uniprot_accession(token):
        return "uniprot_accession"
    if is_agi_locus(token):
        return "agi_locus"
    if is_locus_tag(token):
        return "locus_tag"
    return "gene_symbol"


def is_usable_symbol(value: object) -> bool:
    """Whether a token is usable as a family symbol (not an identifier).

    Empty values, UniProt accessions, AGI loci and locus tags are never usable,
    and neither is a token that starts with a digit (locus-tag fragments such
    as ``7g085020``).  Anything else -- including descriptive protein names --
    counts as usable at the identifier level; semantic quality is judged
    separately by the caller (e.g. the Route B semantic gate).
    """
    if value is None or not str(value).strip():
        return False
    token = str(value).strip()
    return not (
        token[:1].isdigit()
        or is_uniprot_accession(token)
        or is_agi_locus(token)
        or is_locus_tag(token)
    )


def symbol_candidates(
    subject_id: str,
    subject_title: str,
    alias_lookup=None,
) -> list[dict]:
    """Ordered, provenance-annotated family-symbol candidates for a hit.

    ``alias_lookup`` is an optional callable ``str -> str | None`` that maps a
    natural-language protein name to a curated canonical symbol; when provided
    it is applied to the protein-name field.  Candidate sources are:

    - ``gn`` (highest confidence): the ``GN=`` value, when it is a usable symbol
    - ``curated_alias``: protein name resolved through the curated alias map
    - ``protein_name``: a symbol-like single-token protein name (e.g. ``FTa2``)

    Multi-word free-text protein names are never auto-symbolized; without a
    curated alias they produce no candidate.  Callers take the first candidate
    or use the full list for audit.
    """
    candidates: list[dict] = []
    gn = gene_name_field(subject_title)
    if gn and is_usable_symbol(gn):
        candidates.append({"value": gn, "source": "gn", "confidence": "high"})
    name = uniprot_name_field(subject_id, subject_title)
    if name:
        canonical = alias_lookup(name) if alias_lookup else None
        if canonical:
            candidates.append(
                {"value": canonical, "source": "curated_alias", "confidence": "medium"}
            )
        elif is_usable_symbol(name) and " " not in name and len(name) <= 24:
            candidates.append(
                {"value": name, "source": "protein_name", "confidence": "medium"}
            )
    return candidates
