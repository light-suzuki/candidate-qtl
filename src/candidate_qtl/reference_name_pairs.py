"""Bounded, explicit Route A name/species alternatives; never a cross-product."""

from __future__ import annotations

import re

from .identifiers import is_agi_locus, is_uniprot_accession


def reference_handle(value: object, *, allow_untyped: bool = False) -> tuple[str, str] | None:
    """Parse database namespaces, never biological aliases or gene-name guesses."""
    match = re.fullmatch(r"(?:(uniprot(?:kb)?|tair|genbank|embl|refseq|ensembl(?:plants)?|gramene|phytozome|geneid|ncbi\s*gene)\s*[: ]\s*)?([^\s:]+)", str(value).strip(), re.I)
    if not match or (not match[1] and not allow_untyped):
        return None
    namespace, token = (match[1] or "").casefold().replace(" ", ""), match[2]
    if namespace not in {"gramene", "phytozome"}:
        token = token.upper()
    if namespace in {"", "uniprot", "uniprotkb"} and is_uniprot_accession(token):
        return "uniprot", token
    if namespace in {"", "tair"} and is_agi_locus(token):
        return "tair", token
    if namespace in {"geneid", "ncbigene"} and token.isdigit():
        return "geneid", token
    kind = {"genbank": "embl", "embl": "embl", "refseq": "refseq", "ensembl": "ensemblplants", "ensemblplants": "ensemblplants", "gramene": "gramene", "phytozome": "phytozome"}.get(namespace)
    if kind and re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]*[0-9][A-Za-z0-9_.-]*", token):
        return kind, token
    return None


def record_reference_handles(record: dict) -> set[tuple[str, str]]:
    """Use explicit authoritative cross-references; preserve accession versions."""
    handles = set()
    primary = reference_handle(record.get("primaryAccession", ""), allow_untyped=True)
    if primary:
        handles.add(primary)
    for accession in record.get("secondaryAccessions", []):
        secondary = reference_handle(f"UniProt:{accession}")
        if secondary:
            handles.add(secondary)
    for entry in record.get("uniProtKBCrossReferences", []):
        if not isinstance(entry, dict):
            continue
        handle = reference_handle(f"{entry.get('database', '')}:{entry.get('id', '')}")
        if handle:
            handles.add(handle)
        if str(entry.get("database", "")).casefold() == "embl":
            for prop in entry.get("properties", []):
                if isinstance(prop, dict) and prop.get("key") == "ProteinId":
                    protein = reference_handle(f"GenBank:{prop.get('value', '')}")
                    if protein:
                        handles.add(protein)
    return handles


def reference_name_for_input(row: dict) -> str:
    """An explicit external ID is also a lookup handle when no name is known."""
    name = str(row.get("reference_gene_name", "")).strip()
    handle = reference_handle(row.get("proposed_stable_id", ""), allow_untyped=True)
    if not name and handle and handle[0] not in {"uniprot", "tair"}:
        return f"{handle[0]}:{handle[1]}"
    return name


def atomic_reference_pairs(name: str, species: str, *, max_pairs: int = 16) -> tuple[list[tuple[str, str]], str]:
    """Pair parallel semicolon clauses or share one explicitly supplied species.

    Slash/pipe-separated names in a clause are alternatives for that clause's
    species. Parentheses are annotations, not a species alias dictionary.
    Different numbers of name/species clauses cannot establish a pairing.
    """
    if not name.strip() or not species.strip():
        return [], "incomplete_name_species"
    names = [value.strip() for value in name.split(";")]
    species_parts = [" ".join(value.split()) for value in species.split(";")]
    if any(not value for value in names + species_parts):
        return [], "ambiguous_name_species_pairing"
    if len(species_parts) != 1 and len(species_parts) != len(names):
        return [], "ambiguous_name_species_pairing"
    pairs = []
    for index, clause in enumerate(names):
        target_species = species_parts[0] if len(species_parts) == 1 else species_parts[index]
        if any(char in target_species for char in "/|(),"):
            return [], "ambiguous_name_species_pairing"
        # Split outside annotations, then remove one trailing gloss per endpoint.
        for alternative in re.split(r"[/|](?![^()]*\))", clause):
            alternative = re.sub(r"\s*\([^()]*\)\s*$", "", alternative).strip()
            if any(char in alternative for char in "()"):
                return [], "ambiguous_name_species_pairing"
            alternative = " ".join(alternative.split())
            if not alternative or (not reference_handle(alternative) and not re.fullmatch(r"[A-Za-z][A-Za-z0-9_. -]*", alternative)):
                return [], "ambiguous_name_species_pairing"
            pair = (alternative, target_species)
            if pair not in pairs:
                pairs.append(pair)
    if len(pairs) > max_pairs:
        return [], "atomic_reference_budget_exceeded"
    return pairs, "atomic_name_species"


# Rank markers that make a record name an infraspecific taxon of a binomial.
_INFRASPECIFIC_MARKERS = frozenset({"subsp.", "ssp.", "var.", "subvar.", "f.", "cv."})


def record_species_matches(record: dict, species: str) -> bool:
    """Compare the retrieved scientific name with the declared species.

    Whitespace and case are presentation differences. A record that names an
    infraspecific taxon of the declared binomial (``Oryza sativa subsp.
    japonica`` for ``Oryza sativa``) belongs to that species. The reverse is not
    true: a declared subspecies is not confirmed by a species-level record.
    A shared genus, an unconfirmed alias or a different species is never a
    species receipt. This is reference-sequence identity, not a paper audit.
    """
    organism = record.get("organism", {})
    if not isinstance(organism, dict):
        return False
    normalized = lambda value: " ".join(str(value).split()).casefold()
    retrieved, declared = normalized(organism.get("scientificName", "")), normalized(species)
    if not declared:
        return False
    if retrieved == declared:
        return True
    if len(declared.split()) != 2 or not retrieved.startswith(declared + " "):
        return False
    rest = retrieved[len(declared) + 1:].split()
    return len(rest) >= 2 and rest[0] in _INFRASPECIFIC_MARKERS


def species_prefixed_name(name: str, species: str) -> str:
    """Return ``name`` without a leading species abbreviation, or ``""``.

    ``SlGLK2`` for Solanum lycopersicum becomes ``GLK2``. The prefix must be the
    genus initial plus the epithet initial of the declared species, followed by
    an uppercase letter or digit, so ordinary symbols are left alone.
    """
    words = " ".join(str(species).split()).split(" ")
    value = str(name).strip()
    if len(words) < 2 or not words[0] or not words[1] or len(value) < 4:
        return ""
    prefix = words[0][0].upper() + words[1][0].lower()
    rest = value[2:]
    if value[:2] != prefix or not (rest[0].isupper() or rest[0].isdigit()):
        return ""
    return rest


def exact_record_candidates(records: dict, name: str | list[str] | tuple[str, ...] | set[str], species: str) -> list[str]:
    """Match any explicitly declared name for one atomic pair and exact species."""
    normalized = lambda value: " ".join(str(value).split()).casefold()
    raw_names = name if isinstance(name, (list, tuple, set)) else [name]
    wanted = {normalized(value) for value in raw_names if normalized(value)}
    # A species-abbreviation prefix (SlGLK2) names the same symbol (GLK2) in
    # that species; the species check below still applies.
    wanted |= {normalized(stripped) for value in raw_names
               if (stripped := species_prefixed_name(str(value), species))}
    typed = {reference_handle(value) for value in raw_names} - {None}
    if not wanted:
        return []
    found = []
    for accession, record in records.items():
        if not isinstance(record, dict) or record.get("primaryAccession") != accession:
            continue
        if not record_species_matches(record, species):
            continue
        record_handles = {(kind, value.casefold()) for kind, value in record_reference_handles(record)}
        if any((kind, value.casefold()) in record_handles for kind, value in typed):
            found.append(accession)
            continue
        for gene in record.get("genes", []) if isinstance(record.get("genes"), list) else []:
            if not isinstance(gene, dict):
                continue
            main = gene.get("geneName", {})
            names = [main.get("value", "")] if isinstance(main, dict) else []
            for field in ("synonyms", "orderedLocusNames", "orfNames"):
                values = gene.get(field, [])
                names += [entry.get("value", "") for entry in values if isinstance(entry, dict)] if isinstance(values, list) else []
            if any(value and normalized(value) in wanted for value in names):
                found.append(accession)
                break
    return sorted(set(found))


def reference_lookup_names(raw_name: str, atomic_name: str = "", species: str = "",
                           input_species: str = "") -> list[str]:
    """Return bounded complete name forms worth querying for one atomic pair.

    Keep formal names intact (for example ``TERMINAL FLOWER 1``); the
    authoritative record check remains responsible for deciding whether a
    result is actually the declared gene.
    """
    source = str(raw_name or "")
    names = [str(atomic_name or "").strip()]
    clauses = [source]
    if species:
        raw_clauses = [value.strip() for value in source.split(";")]
        species_clauses = [" ".join(value.split()) for value in str(input_species or species).split(";")]
        clauses = []
        if len(species_clauses) == 1 or len(species_clauses) == len(raw_clauses):
            target = " ".join(str(species).split()).casefold()
            for index, clause in enumerate(raw_clauses):
                clause_species = species_clauses[0] if len(species_clauses) == 1 else species_clauses[index]
                if clause_species.casefold() != target:
                    continue
                for branch in re.split(r"[/|](?![^()]*\))", clause):
                    alternative = re.sub(r"\s*\([^()]*\)\s*$", "", branch).strip()
                    if " ".join(alternative.split()).casefold() == " ".join(str(atomic_name or "").split()).casefold():
                        clauses.append(branch)
    # Parenthetical text is often a formal name or alternate symbol. Do not
    # interpret it as a species mapping; atomic_reference_pairs owns that.
    for clause in clauses:
        names.extend(value for value in re.findall(r"\(([^()]*)\)", clause) if value.strip())
        for value in re.split(r"[;/|]", clause):
            value = re.sub(r"\s*\([^()]*\)\s*$", "", value).strip()
            if value:
                names.append(value)
    output = []
    for value in names:
        value = " ".join(value.split())
        if value and value not in output:
            output.append(value)
    return output


def reference_lookup_names_for_inputs(raw_inputs: list[dict], atomic_name: str, species: str) -> list[str]:
    """Collect only names explicitly attached to this atomic name/species pair."""
    names = [str(atomic_name or "").strip()]
    for row in raw_inputs:
        raw_name = str(row.get("input_reference_gene_name") or row.get("reference_gene_name", ""))
        raw_species = str(row.get("input_reference_species") or row.get("reference_species", ""))
        names.extend(reference_lookup_names(raw_name, atomic_name, species or raw_species, raw_species))
    output = []
    for name in names:
        name = " ".join(str(name).split())
        if name and name.casefold() not in {value.casefold() for value in output}:
            output.append(name)
    return output


def reference_lookup_handles(raw_inputs: list[dict]) -> list[tuple[str, str]]:
    """Extract only explicitly typed, syntactically valid ID lookup handles."""
    handles: list[tuple[str, str]] = []
    for row in raw_inputs:
        proposed = str(row.get("proposed_stable_id", "")).strip()
        values = [(proposed, True)]
        values.append((str(row.get("reference_gene_name", "")), False))
        source_ids = str(row.get("source_ids", ""))
        values.extend((value.strip(), False) for value in re.split(r"[;,|]+", source_ids) if value.strip())
        for value, allow_untyped in values:
            handle = reference_handle(value, allow_untyped=allow_untyped)
            if handle and handle not in handles:
                handles.append(handle)
    return handles


def reference_identity_names(raw_inputs: list[dict], atomic_name: str, species: str) -> list[str]:
    """Declared names plus typed external handles, without guessing aliases.

    A supplied UniProt accession alone cannot override a conflicting gene name.
    Other database IDs are checked against that record's explicit cross-references.
    """
    names = reference_lookup_names_for_inputs(raw_inputs, atomic_name, species)
    names.extend(f"{kind}:{value}" for kind, value in reference_lookup_handles(raw_inputs)
                 if kind != "uniprot")
    return list(dict.fromkeys(names))


# Identity outcomes for one retrieved record against a named, declared target.
IDENTITY_EXACT = "exact_name_species"
IDENTITY_REVIEWED_ALIAS = "reviewed_name_alias"
IDENTITY_DECLARED_ACCESSION = "declared_accession_species_name_unconfirmed"
IDENTITY_MISMATCH = "mismatch"
CONFIRMED_IDENTITIES = frozenset({IDENTITY_EXACT, IDENTITY_REVIEWED_ALIAS})
# Identities that may carry the declared target's papers to projection.
# The declared-accession state is provisional: review-only, never support.
PROJECTABLE_IDENTITIES = CONFIRMED_IDENTITIES | {IDENTITY_DECLARED_ACCESSION}
IDENTITY_REVIEW_COLUMNS = ("reference_gene_name", "reference_species", "accession",
                           "record_gene_name", "verdict", "evidence_source_ids", "reviewer", "rationale")


def declared_uniprot_accessions(raw_inputs: list[dict]) -> set[str]:
    """UniProt accessions the worker itself declared as ``proposed_stable_id``.

    Accessions merely found by a name or citation search are not declarations.
    """
    accessions = set()
    for row in raw_inputs:
        handle = reference_handle(str(row.get("proposed_stable_id", "")).strip(), allow_untyped=True)
        if handle and handle[0] == "uniprot":
            accessions.add(handle[1].upper())
    return accessions


def record_gene_names(record: dict) -> list[str]:
    """Gene names, synonyms and locus names written on a provider record."""
    names: list[str] = []
    for gene in record.get("genes", []) if isinstance(record.get("genes"), list) else []:
        if not isinstance(gene, dict):
            continue
        main = gene.get("geneName", {})
        if isinstance(main, dict) and main.get("value"):
            names.append(str(main["value"]))
        for field in ("synonyms", "orderedLocusNames", "orfNames"):
            values = gene.get(field, [])
            if isinstance(values, list):
                names.extend(str(entry.get("value", "")) for entry in values
                             if isinstance(entry, dict) and entry.get("value"))
    return list(dict.fromkeys(names))


def declared_record_identity(record: dict, accession: str, names, species: str,
                             declared_accessions: set[str] | frozenset[str] = frozenset(),
                             identity_reviews: dict | None = None) -> str:
    """Decide how far one record confirms a named, species-qualified target.

    * ``exact_name_species``: the record carries a declared name and species.
    * ``reviewed_name_alias``: the names differ, but an independent identity
      review approved this exact name/species/accession and named a gene name
      that the record really carries (the record check stays mechanical).
    * ``declared_accession_species_name_unconfirmed``: the worker declared this
      accession, the record exists with that primary accession and the species
      matches, but no declared name is on the record. Provisional only.
    * ``mismatch``: anything else (another species, an accession found only by
      search, or a record whose primary accession differs).
    """
    if not isinstance(record, dict) or record.get("primaryAccession") != accession:
        return IDENTITY_MISMATCH
    if accession in exact_record_candidates({accession: record}, names, species):
        return IDENTITY_EXACT
    if accession.upper() not in {value.upper() for value in declared_accessions}:
        return IDENTITY_MISMATCH
    if not record_species_matches(record, species):
        return IDENTITY_MISMATCH
    raw_names = names if isinstance(names, (list, tuple, set)) else [names]
    normalized = lambda value: " ".join(str(value).split()).casefold()
    on_record = {normalized(value) for value in record_gene_names(record)}
    for name in raw_names:
        review = (identity_reviews or {}).get((normalized(name), normalized(species), accession.upper()))
        if review and normalized(review.get("record_gene_name", "")) in on_record:
            return IDENTITY_REVIEWED_ALIAS
    return IDENTITY_DECLARED_ACCESSION


def load_identity_reviews(path) -> dict:
    """Read an identity review ledger; only ``same_gene`` rows with evidence count.

    Keyed by normalized (worker's reference name, species, accession). A row
    that is not ``same_gene`` or lacks evidence or a reviewer is ignored, so a
    missing or negative review leaves the target provisional.
    """
    import csv
    from pathlib import Path as _Path
    if not path:
        return {}
    file = _Path(path)
    with file.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if rows and set(IDENTITY_REVIEW_COLUMNS) - set(rows[0]):
        raise ValueError("identity review ledger is missing columns: "
                         + ", ".join(sorted(set(IDENTITY_REVIEW_COLUMNS) - set(rows[0]))))
    normalized = lambda value: " ".join(str(value).split()).casefold()
    reviews = {}
    for row in rows:
        if (str(row.get("verdict", "")).strip() != "same_gene"
                or not str(row.get("evidence_source_ids", "")).strip()
                or not str(row.get("reviewer", "")).strip()
                or not str(row.get("record_gene_name", "")).strip()):
            continue
        accession = str(row.get("accession", "")).strip().split(":")[-1].upper()
        key = (normalized(row.get("reference_gene_name", "")), normalized(row.get("reference_species", "")), accession)
        if key in reviews and reviews[key] != row:
            raise ValueError(f"conflicting identity reviews for {key}")
        reviews[key] = dict(row)
    return reviews
