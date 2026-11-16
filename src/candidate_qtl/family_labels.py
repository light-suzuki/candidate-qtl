"""Shared, auditable normalization for biological family labels.

Route A often emits descriptive labels while homology-backed routes emit gene
symbols.  Exact string equality therefore undercounts genuine family recovery.
This module provides a deliberately small many-to-one alias registry and keeps
the raw normalized label as a key.  It does not infer pathway membership and it
never changes a score.
"""

from __future__ import annotations

import re


_SPACE = re.compile(r"\s+")
_NON_ALNUM = re.compile(r"[^a-z0-9]+")
_TRAILING_COPY = re.compile(r"^([a-z][a-z0-9_-]*?[a-z])[-_.]?\d+$", re.IGNORECASE)
# Uppercase symbol-like tokens in a descriptive Route A label ("CIK", "D14",
# "CLV1", "MAX2").  Tokens must be fully uppercase alphanumeric and at least
# two characters so generic words ("protein", "factor") never become keys.
_SYMBOL_TOKEN = re.compile(r"[A-Z][A-Z0-9]+")
# Trailing Roman numerals converted to Arabic digits in the expanded tokenizer
#: SBEI -> SBE1, SBEII -> SBE2, GBSSI -> GBSS1.  Longest suffix
# first so "SBEIII" converts as 3, not as 1.  X/V/L/C/D/M are not numerals
# here, so X-ending names (KNOX, POLLUX, GA20OX, LYKX, nodX) are untouched.
_ROMAN_MAP = {"III": "3", "II": "2", "I": "1"}
_ROMAN_SUFFIX = re.compile(r"(III|II|I)$", re.I)
# Real words whose trailing I/II/III is a letter, not a numeral; converting
# them would fabricate a false synonym with a real gene (COPI -> COP1 pairs
# with the photomorphogenesis E3 ligase).  Curated; do not grow ad hoc.
_ROMAN_EXCLUDED = frozenset({"copi", "taxi", "kari", "psii", "coii", "nodi"})
# Species prefix stripped before the Roman-suffix check ("PsSBEI" -> "SBEI").
# Only legume/plant prefixes that commonly prefix symbols as PsSBEI / MtSBEI.
# Case-sensitive on purpose: "Sb" must not match "SB" inside "SBEI".
_SPECIES_PREFIX = re.compile(r"^(?:Ps|Mt|Lj|Os|At|Ta|Gm|Zm|Sl|St|Ca|Sb|Hv)")

# Curated equivalence classes needed at the discovery-family boundary.  The
# patterns are intentionally specific: broad words such as "protein" or
# "transcription factor" are never keys by themselves.
_ALIASES: tuple[tuple[str, tuple[re.Pattern[str], ...]], ...] = (
    ("phy", (
        re.compile(r"\bphytochrome\b", re.I),
        re.compile(r"\bphy(?:a|b|c|d|e)?\b", re.I),
    )),
    ("mlo", (
        re.compile(r"\bmlo\b", re.I),
        re.compile(r"\bmildew\s+resistance\s+locus\s+o\b", re.I),
    )),
    ("ttg", (
        re.compile(r"\bttg\d*\b", re.I),
        re.compile(r"\btransparent\s+testa\s+glabra\s*\d*\b", re.I),
    )),
    ("nlp", (
        re.compile(r"\bnlp\d*\b", re.I),
        re.compile(r"\bnin[- ]like\b", re.I),
        re.compile(r"\brwp[- ]rk\b", re.I),
    )),
    ("lfy", (
        re.compile(r"\blfy\b", re.I),
        re.compile(r"\bleafy\b", re.I),
        re.compile(r"\bfloricaula\b", re.I),
    )),
    ("ft", (
        re.compile(r"\bft\b", re.I),
        re.compile(r"\bflowering\s+locus\s+t\b", re.I),
        re.compile(r"\bpebp\b.*\bft(?:[- ]clade)?\b", re.I),
        # Issue #297: TrEMBL protein names spell the FT family "Flowering
        # locus protein T" (the protein head noun sits between locus and T)
        # or "Flowering locus T-like"; the FT-clade is also written without
        # a preceding PEBP word ("FT clade"), matching the TFL1/CEN/SP clade.
        re.compile(r"\bflowering\s+locus\s+(?:protein\s+)?t(?:-like)?\b", re.I),
        re.compile(r"\bft\s+clade\b", re.I),
    )),
    ("if4e", (
        re.compile(r"\be?if4e\d*\b", re.I),
        re.compile(r"\beukaryotic\s+translation\s+initiation\s+factor\s+4e\b", re.I),
    )),
    ("wus", (
        re.compile(r"\bwus\b", re.I),
        re.compile(r"\bwuschel\b", re.I),
    )),
    # Issue #297: NSP2 (pea Sym7) is the nodulation-signaling pathway 2 GRAS
    # transcription factor; the Medicago TrEMBL defline spells the protein
    # "Nodulation-signaling pathway 2".  Both the "nodulation" and "nodule"
    # spellings of the head noun are established.
    ("nsp2", (
        re.compile(r"\bnsp2\b", re.I),
        re.compile(r"\bnodul(?:e|ation)[- ]?sign(?:al)?ing\s+pathway\s+2\b", re.I),
    )),
    # Issue #301: SBEI and SBE1 are the same starch-branching-enzyme I locus
    # (pea R/r, Mendel trait 2); the Roman and Arabic spellings are both
    # established in the literature and in UniProt (SBEI_PEA).
    ("sbe1", (
        re.compile(r"\bsbei\b", re.I),
    )),
    # Issue #301: GBSSI and GBSS1 are the same granule-bound starch synthase I
    # locus (pea lam); both spellings are established (Q43092 GBSS1_PEA).
    ("gbss1", (
        re.compile(r"\bgbssi\b", re.I),
    )),
    # Issue #301: storage-protein family head nouns (pea vic/leg/psl/sca/cvc
    # cases).  Only the protein-family head noun itself is a pattern; function
    # words (kinase, receptor, protein, enzyme, factor, ...) are never keys.
    # "\bvicilin\b" does not match inside "convicilin" (no word boundary), so
    # the two families stay separate.
    ("vicilin", (
        re.compile(r"\bvicilin\b", re.I),
    )),
    ("convicilin", (
        re.compile(r"\bconvicilin\b", re.I),
    )),
    ("legumin", (
        re.compile(r"\blegumin\b", re.I),
    )),
    ("lectin", (
        re.compile(r"\blectin\b", re.I),
    )),
    ("albumin", (
        re.compile(r"\balbumin\b", re.I),
    )),
    # Issue #411: BAN (BANYULS) is anthocyanidin reductase (ANR) in Arabidopsis
    # (UniProtKB:Q9SEV0, GN=BAN, Synonyms=ANR, Description="Anthocyanidin reductase").
    ("anr", (
        re.compile(r"\banr\b", re.I),
        re.compile(r"\bban\b", re.I),
        re.compile(r"\bbanyuls\b", re.I),
        re.compile(r"\banthocyanidin\s+reductase\b", re.I),
    )),
    # Issue #411: LUP (lupeol synthase) family in Arabidopsis
    # (UniProtKB:P93033, GN=LUP1, Description="Lupeol synthase 1").
    ("lup", (
        re.compile(r"\blup\d*\b", re.I),
        re.compile(r"\blupeol\s+synthase\b", re.I),
    )),
    # Issue #419: BRC1 / TB1 / TCP18 are the same shoot-branching TCP
    # transcription factor. Arabidopsis BRANCHED1 has locus name TCP18
    # (UniProtKB:A1YKT1, GN=TCP18, Synonyms=BRC1/TBL1, Description="Protein
    # BRANCHED 1"), maize TEOSINTE BRANCHED 1 is the ortholog (UniProtKB:Q93WI2,
    # GN=TB1), and pea PsBRC1 is the pea ortholog. The species-prefixed
    # spelling (PsBRC1) is covered by the optional ps- prefix on brc1. The
    # branched pattern is spelled as explicit alternatives ("BRANCHED1" and
    # "BRANCHED 1") rather than a zero-width "\s*" so descriptions such as
    # "branched 1,3-beta-glucanase" cannot ride the alias.
    ("brc1", (
        re.compile(r"\b(?:ps)?brc1\b", re.I),
        re.compile(r"\bbranched1\b", re.I),
        re.compile(r"\bbranched\s+1\b", re.I),
        re.compile(r"\btb1\b", re.I),
        re.compile(r"\bteosinte\s+branched\s*1\b", re.I),
        re.compile(r"\btcp18\b", re.I),
    )),
    # Arabidopsis floral-homeotic PISTILLATA is named PI in the curated
    # UniProtKB gene field. Route A commonly emits the compact PI symbol
    # inside an ABC-model/MADS family label, while homology databases expose
    # the full PISTILLATA protein name. Bare ``PI`` is deliberately excluded:
    # it is ambiguous with phosphoinositide (PI) kinase families.
    ("pistillata", (
        re.compile(
            r"^(?![^\n]*\b(?:non|not|without|excluding|except)\b[- ]*pistillata\b).*"
            r"\bpistillata\b",
            re.I,
        ),
        # Parenthesized/listed symbol roster, e.g. the ABC-model label
        # "MADS-box floral homeotic factors (AP1/AP3/PI/AG)". Requiring a
        # delimiter immediately after PI excludes "PI kinase" even when a
        # mixed label also happens to mention MADS.
        re.compile(
            r"^(?![^\n]*(?:\b(?:non|not|without|excluding|except)\b[- ]*pi\b|"
            r"\bpi[- ](?:independent|negative)\b)).*"
            r"\bmads(?:[- ]box)?\b[^()\n]{0,120}"
            r"\([^)\n]{0,120}\bpi\s*(?=[/,;)])",
            re.I,
        ),
        re.compile(
            r"^(?![^\n]*(?:\b(?:non|not|without|excluding|except)\b[- ]*pi\b|"
            r"\bpi[- ](?:independent|negative)\b)).*"
            r"\bpi[- ]class\b[^\n]{0,80}\bmads(?:[- ]box)?\b",
            re.I,
        ),
        re.compile(
            r"^(?![^\n]*(?:\b(?:non|not|without|excluding|except)\b[- ]*pi\b|"
            r"\bpi[- ](?:independent|negative)\b)).*"
            r"\bmads(?:[- ]box)?\b[^\n]{0,80}\bpi[- ]class\b",
            re.I,
        ),
    )),
)

# Auditable public database citations for every curated alias in _ALIASES.
ALIAS_PROVENANCE: dict[str, dict[str, str]] = {
    "phy": {"database": "UniProtKB", "accession": "P14712", "field": "GN/Description", "note": "Phytochrome family"},
    "mlo": {"database": "UniProtKB", "accession": "Q9S7U8", "field": "GN/Description", "note": "MLO family"},
    "ttg": {"database": "UniProtKB", "accession": "Q9FN69", "field": "GN/Description", "note": "TTG1 / Transparent Testa Glabra"},
    "nlp": {"database": "UniProtKB", "accession": "Q9M341", "field": "GN/Description", "note": "NIN-like protein"},
    "lfy": {"database": "UniProtKB", "accession": "P29382", "field": "GN/Description", "note": "LEAFY / FLORICAULA"},
    "ft": {"database": "UniProtKB", "accession": "Q9SVR7", "field": "GN/Description", "note": "Flowering Locus T / PEBP"},
    "if4e": {"database": "UniProtKB", "accession": "P48593", "field": "GN/Description", "note": "eIF4E translation initiation factor"},
    "wus": {"database": "UniProtKB", "accession": "Q9SBG4", "field": "GN/Description", "note": "WUSCHEL"},
    "nsp2": {"database": "UniProtKB", "accession": "Q5NE24", "field": "GN/Description", "note": "Nodulation-signaling pathway 2"},
    "sbe1": {"database": "UniProtKB", "accession": "Q43092", "field": "GN/Description", "note": "Starch branching enzyme 1 / SBEI"},
    "gbss1": {"database": "UniProtKB", "accession": "Q43093", "field": "GN/Description", "note": "Granule-bound starch synthase 1 / GBSSI"},
    "vicilin": {"database": "UniProtKB", "accession": "P08438", "field": "Description", "note": "Vicilin storage protein"},
    "convicilin": {"database": "UniProtKB", "accession": "P13918", "field": "Description", "note": "Convicilin storage protein"},
    "legumin": {"database": "UniProtKB", "accession": "P02857", "field": "Description", "note": "Legumin storage protein"},
    "lectin": {"database": "UniProtKB", "accession": "P02867", "field": "Description", "note": "Lectin storage protein"},
    "albumin": {"database": "UniProtKB", "accession": "P13917", "field": "Description", "note": "Albumin storage protein"},
    "anr": {"database": "UniProtKB", "accession": "Q9SEV0", "field": "GN=BAN, Synonyms=ANR, Description=Anthocyanidin reductase", "note": "BANYULS / anthocyanidin reductase"},
    "lup": {"database": "UniProtKB", "accession": "P93033", "field": "GN=LUP1, Description=Lupeol synthase 1", "note": "Lupeol synthase family"},
    "brc1": {"database": "UniProtKB", "accession": "A1YKT1; Q93WI2", "field": "A1YKT1: GN=TCP18, Synonyms=BRC1/TBL1, Description='Protein BRANCHED 1'; Q93WI2: GN=TB1, Description='Transcription factor TEOSINTE BRANCHED 1'", "note": "BRANCHED1 / TEOSINTE BRANCHED1 / TCP18 (pea PsBRC1 ortholog)"},
    "pistillata": {"database": "UniProtKB", "accession": "P48007", "field": "GN=PI, Description=Floral homeotic protein PISTILLATA", "note": "PI in MADS/floral-homeotic context / PISTILLATA"},
}

# Expanded-only curated aliases.  FTa1/FTa2 (GIGAS) and
# FTb1/FTb2 are FT-clade PEBP paralogs whose protein names ("FTa2") must meet
# a strict truth symbol "FT" in the Route A (expanded) matcher, but the
# paralog patterns deliberately stay out of the strict layer: leave-family-out
# masking and decoy checks must never conflate "FT" with its FTa/FTb names
# without a shared explicit symbol, and families_match("FTa2 protein", "FT")
# stays False (the tokenizer bridge is expanded-only by design).
_EXPANDED_ALIASES: tuple[tuple[str, tuple[re.Pattern[str], ...]], ...] = (
    ("ft", (
        re.compile(r"\bfta\d*\b", re.I),
        re.compile(r"\bftb\d*\b", re.I),
    )),
)


def normalized_family_label(value: object) -> str:
    """Normalize punctuation and whitespace without discarding label content."""
    text = _NON_ALNUM.sub(" ", str(value or "").casefold())
    return _SPACE.sub(" ", text).strip()


def _camel_symbol_key(token: str) -> str | None:
    """Symbol key from a lowercase-mixed symbol token ("CCaMK", "Nup85", "GA3ox").

    Conservative structure guard (issue #301): the token must start uppercase
    and contain either a second uppercase letter or a digit, so generic words
    (Protein, Kinase, Receptor, Factor, Family, Complex, Domain, Enzyme,
    Locus, Le, Nup, Pea, ...) never become keys.  Only the full token becomes
    a key -- no leading-run sub-key -- so "FTa2" cannot leak a bare "ft"
    symbol key and "Nup85" keeps "nup85" without a bare "nup".
    """
    if not token[:1].isupper():
        return None
    if sum(1 for ch in token if ch.isupper()) < 2 and not any(ch.isdigit() for ch in token):
        return None
    return f"symbol:{token.casefold()}"


def _roman_symbol_key(token: str) -> str | None:
    """Arabic-digit key for a trailing Roman I/II/III suffix ("SBEI" -> SBE1).

    Expanded-tokenizer only (issue #301).  The species prefix (PsSBEI -> SBEI)
    is stripped first, then the suffix must sit at the absolute end of a token
    of at least four characters, so short real symbols that end in "I" by
    coincidence (CHI, LOX, LUX, WOX, UNI, GAI, CLV, BBI, CRI, DOI, NFP) and
    X-ending names (KNOX, POLLUX, GA20OX, LYKX, nodX) are never converted;
    "SBEIIa" keeps its suffix because "a" follows.  Longer real words whose
    trailing I/II/III is a letter, not a numeral, are excluded explicitly
    (COPI -> COP1 is a real false-pairing with the E3 ligase).  The converted
    key is added next to the original key, which is always kept.
    """
    stripped = _SPECIES_PREFIX.sub("", token, count=1)
    if len(stripped) < 4:
        return None
    compact = stripped.casefold()
    if compact in _ROMAN_EXCLUDED:
        return None
    match = _ROMAN_SUFFIX.search(stripped)
    if match is None:
        return None
    return f"symbol:{stripped[: match.start()].casefold()}{_ROMAN_MAP[match.group(1).upper()]}"


def _symbol_tokens_split(value: str) -> tuple[set[str], set[str]]:
    """Split ``_symbol_tokens`` keys into exact-token and stem-only keys.

    See ``_identity_keys_split`` for what "exact" vs "stem" means; this is the
    same split applied to the uppercase symbol tokens pulled out of a
    descriptive Route A label rather than to the label/value itself.  The
    exact bucket also carries camelCase keys (CCaMK, Nup85, GA3ox) and
    Roman-converted keys (SBEI -> SBE1); all three sources stay expanded-only.
    """
    exact: set[str] = set()
    stem_keys: set[str] = set()
    for token in re.split(r"[^A-Za-z0-9]+", value):
        if not token:
            continue
        if _SYMBOL_TOKEN.fullmatch(token):
            compact = token.casefold()
            exact.add(f"symbol:{compact}")
            stem = _TRAILING_COPY.sub(r"\1", compact)
            if len(stem) >= 3 and stem != compact:
                stem_keys.add(f"symbol:{stem}")
            roman = _roman_symbol_key(token)
            if roman is not None:
                exact.add(roman)
        else:
            camel = _camel_symbol_key(token)
            if camel is not None:
                exact.add(camel)
                roman = _roman_symbol_key(token)
                if roman is not None:
                    exact.add(roman)
    return exact, stem_keys


def _symbol_tokens(value: str) -> set[str]:
    """Symbol keys from uppercase tokens inside a descriptive label.

    Route A emits labels such as "CIK LRR receptor kinases" while homology-
    backed lanes and the benchmark truth use gene symbols ("CIK1").  Adding the
    stemmed uppercase tokens as symbol keys lets the two namespaces meet
    (issue #248).  Conservative: a token must be fully uppercase alphanumeric
    and at least two characters, so lowercase generic words are never keys.
    """
    exact, stem_keys = _symbol_tokens_split(value)
    return exact | stem_keys


def _identity_keys_split(value: object) -> tuple[set[str], set[str]]:
    """Split a value's identity keys into exact-match keys and stem-only keys.

    "Exact" keys (normalized label, full symbol, curated alias) express a
    specific identity claim. "Stem" keys additionally fold off a trailing
    paralog numeral (``CIK1`` -> ``symbol:cik``), so ``CIK1`` and ``CIK6``
    share a stem key without sharing an exact one (issue #268: this is why
    stem-level family matches can conflate distinct paralogs). When the
    stemmed symbol equals the full symbol, ``family_identity_keys`` already
    collapses them into a single ``symbol:`` key, so nothing is duplicated
    into ``stem_keys`` in that case -- it stays an exact key.
    """
    raw = str(value or "").strip()
    if not raw:
        return set(), set()
    normalized = normalized_family_label(raw)
    exact: set[str] = {f"label:{normalized}"} if normalized else set()
    stem_keys: set[str] = set()
    compact = _NON_ALNUM.sub("", raw.casefold())
    if compact:
        exact.add(f"symbol:{compact}")
        stem = _TRAILING_COPY.sub(r"\1", compact)
        if len(stem) >= 3 and stem != compact:
            stem_keys.add(f"symbol:{stem}")
    for canonical, patterns in _ALIASES:
        if any(pattern.search(raw) for pattern in patterns):
            exact.add(f"alias:{canonical}")
    return exact, stem_keys


def family_identity_keys(value: object) -> set[str]:
    """Return conservative identity keys for a family label or gene symbol."""
    exact, stem_keys = _identity_keys_split(value)
    return exact | stem_keys


def _expanded_alias_keys(value: object) -> set[str]:
    """Curated alias keys that only exist in the expanded (Route A) matcher.

    ``_EXPANDED_ALIASES`` bridges Route A protein names to strict truth
    symbols without loosening strict matching anywhere else: masking,
    decoys and negative checks keep using ``family_identity_keys``.
    """
    raw = str(value or "")
    return {
        f"alias:{canonical}"
        for canonical, patterns in _EXPANDED_ALIASES
        if any(pattern.search(raw) for pattern in patterns)
    }


def family_identity_keys_expanded(value: object) -> set[str]:
    """Strict keys plus uppercase symbol tokens from a descriptive label.

    Route A emits labels such as "CIK LRR receptor kinases" while homology-
    backed lanes and the benchmark truth use gene symbols ("CIK1").  The
    expanded keys let the two namespaces meet for Route A evaluation
    (issue #248) without loosening strict matching elsewhere: masking
    (leave-family-out) and negative-control checks keep using
    ``family_identity_keys`` so they never over-remove or over-flag.
    """
    keys = family_identity_keys(value)
    raw = str(value or "").strip()
    keys |= _symbol_tokens(raw)
    keys |= _expanded_alias_keys(raw)
    return keys


def families_match_expanded(left: object, right: object) -> bool:
    """Whether a descriptive label and a symbol share an explicit identity key.

    Only the left-hand value is expanded, so a Route A label can meet a truth
    symbol while the truth side stays strict.
    """
    left_keys = family_identity_keys_expanded(left)
    return bool(left_keys and left_keys & family_identity_keys(right))


def family_set_contains_expanded(
    values: set[str] | list[str] | tuple[str, ...], truth: object
) -> bool:
    """Like family_set_contains, but descriptive values get expanded keys."""
    return any(families_match_expanded(value, truth) for value in values)


def families_match(left: object, right: object) -> bool:
    """Whether two non-empty labels share an explicit identity key."""
    left_keys = family_identity_keys(left)
    return bool(left_keys and left_keys & family_identity_keys(right))


def canonical_family_key(value: object) -> str:
    """Stable grouping key; curated aliases win, then the shortest symbol key."""
    keys = family_identity_keys(value)
    for prefix in ("alias:", "symbol:", "label:"):
        candidates = sorted((key for key in keys if key.startswith(prefix)), key=lambda key: (len(key), key))
        if candidates:
            return candidates[0]
    return ""


def family_set_contains(values: set[str] | list[str] | tuple[str, ...], truth: object) -> bool:
    """Match one truth family against a collection using the shared registry."""
    return any(families_match(value, truth) for value in values)


def alias_canonical(value: object) -> str | None:
    """First curated alias canonical matched by the value, or None.

    ``_ALIASES`` is checked in registry order, then the expanded-only set,
    so "Flowering locus protein T" resolves to "ft" and an unmatched protein
    name (e.g. "CCaMK kinase family") resolves to None.  This is a read-only
    lookup over the curated patterns; it never modifies matching behaviour.
    """
    raw = str(value or "")
    for canonical, patterns in _ALIASES:
        if any(pattern.search(raw) for pattern in patterns):
            return canonical
    for canonical, patterns in _EXPANDED_ALIASES:
        if any(pattern.search(raw) for pattern in patterns):
            return canonical
    return None


# --------------------------------------------------------------------- #
# Issue #268: report *how* a match was made, not just whether it was.   #
# A stem key alone (CIK1 <-> CIK6, BAM1 <-> BAM3, CLE12 <-> CLE25, ...)  #
# is a weaker claim than a shared exact symbol, normalized label, or    #
# curated alias -- paralogs collapse into the same family under stem    #
# matching. These functions never change families_match / family_set_  #
# contains (or their *_expanded siblings): every current hit stays a    #
# hit. They only add a read of which tier produced it.                  #
# --------------------------------------------------------------------- #

MATCH_LEVEL_EXACT = "exact_symbol"
MATCH_LEVEL_STEM = "stem"


def family_match_level(left: object, right: object) -> str:
    """How two labels matched: "exact_symbol", "stem", or "" for no match."""
    left_exact, left_stem = _identity_keys_split(left)
    right_exact, right_stem = _identity_keys_split(right)
    if left_exact & right_exact:
        return MATCH_LEVEL_EXACT
    if (left_exact | left_stem) & (right_exact | right_stem):
        return MATCH_LEVEL_STEM
    return ""


def family_match_level_expanded(left: object, right: object) -> str:
    """Expanded (Route A) variant of ``family_match_level``.

    Only the left-hand value is expanded with uppercase symbol tokens,
    mirroring ``families_match_expanded``.
    """
    left_exact, left_stem = _identity_keys_split(left)
    raw = str(left or "").strip()
    token_exact, token_stem = _symbol_tokens_split(raw)
    left_exact = left_exact | token_exact | _expanded_alias_keys(raw)
    left_stem = left_stem | token_stem
    right_exact, right_stem = _identity_keys_split(right)
    if left_exact & right_exact:
        return MATCH_LEVEL_EXACT
    if (left_exact | left_stem) & (right_exact | right_stem):
        return MATCH_LEVEL_STEM
    return ""


def family_set_match_level(values: set[str] | list[str] | tuple[str, ...], truth: object) -> str:
    """Best match level across a collection: "exact_symbol" > "stem" > ""."""
    best = ""
    for value in values:
        level = family_match_level(value, truth)
        if level == MATCH_LEVEL_EXACT:
            return level
        if level == MATCH_LEVEL_STEM:
            best = level
    return best


def family_set_match_level_expanded(values: set[str] | list[str] | tuple[str, ...], truth: object) -> str:
    """Like ``family_set_match_level``, but values get expanded keys (Route A)."""
    best = ""
    for value in values:
        level = family_match_level_expanded(value, truth)
        if level == MATCH_LEVEL_EXACT:
            return level
        if level == MATCH_LEVEL_STEM:
            best = level
    return best


def split_respecting_brackets(text: object, separator: str = ";") -> list[str]:
    """Split text by separator only when outside parentheses, brackets, and quotes.

    Ensures that descriptive strings such as 'MBW complex (TT2/TT8/TTG1 in Arabidopsis; AN2/AN1 in Petunia)'
    are not broken in the middle of parentheses or quotes.
    """
    raw = str(text or "").strip()
    if not raw:
        return []
    parts: list[str] = []
    current: list[str] = []
    stack: list[str] = []
    in_quote: str | None = None
    bracket_pairs = {"(": ")", "[": "]", "{": "}"}
    for ch in raw:
        if in_quote:
            if ch == in_quote:
                in_quote = None
            current.append(ch)
        elif ch in {'"', "'"}:
            in_quote = ch
            current.append(ch)
        elif ch in bracket_pairs:
            stack.append(bracket_pairs[ch])
            current.append(ch)
        elif stack and ch == stack[-1]:
            stack.pop()
            current.append(ch)
        elif not stack and ch == separator:
            token = "".join(current).strip()
            if token:
                parts.append(token)
            current = []
        else:
            current.append(ch)
    token = "".join(current).strip()
    if token:
        parts.append(token)
    return parts
