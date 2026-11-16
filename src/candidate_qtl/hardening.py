from __future__ import annotations

import gzip
import hashlib
import json
import re
from urllib.parse import urlencode

import pandas as pd

from . import schemas
from .core import europe_pmc_search, file_fingerprint, open_text, runtime_path, safe_id
from .evidence import orthology_plan
from .evidence_io import (
    CONFIDENCE_MULTIPLIERS,
    EVIDENCE_COLUMNS,
    STRONG_SINGLE_CODES,
    annotation_aliases as _annotation_aliases,
    candidate_frame as _candidate_frame,
    id_map as _id_map,
    evidence_row as _row,
)
from .evidence_policy import evidence_build
from .ranking_policy import INDEPENDENCE_FAMILIES, load_weight_registry, load_weights


LLM_IMPORT_CODES = {
    "direct_model_genetic",
    "ortholog_model_genetic",
    "known_qtg",
    "trait_function_match",
    "domain_supported",
    "synteny_supported",
    "counterevidence",
}

# A citable source is a stable external identifier (publication, database
# accession, ontology term) or a concrete local analysis file -- something that
# can anchor an independence claim. A bare free-text label or empty value cannot.
#
# Recognition is allowlist-driven. A namespaced id is citable only
# when its local part matches that namespace's accession grammar, and a bare
# token is citable only when it matches a specific, well-known accession grammar
# as a whole uppercase token. Unknown namespaces and unknown compact tokens
# (``gene1234``, ``trait2026``, ``abc_12345``, ``GENE0001``) are NOT auto-citable;
# they fall through to ``label`` so a free label can never inflate
# independent_source_count or promote a candidate to Top.
_UNIPROT = r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"  # UniProt 6/10-char
_AGI = r"AT[1-5CM]G[0-9]{5}(?:\.[0-9]+)?"                                                  # Arabidopsis AGI locus
_REFSEQ = r"(?:AC|AP|NC|NG|NM|NP|NR|NT|NW|NZ|WP|XM|XP|XR|YP)_[0-9]{5,}(?:\.[0-9]+)?"       # RefSeq (underscore prefix)
_INSDC = r"(?:[A-Z][0-9]{5}|[A-Z]{2}[0-9]{6}|[A-Z]{2}[0-9]{8}|[A-Z]{3}[0-9]{5}|[A-Z]{3}[0-9]{7})(?:\.[0-9]+)?"  # GenBank/EMBL/DDBJ nt/protein
_PMC = r"PMC[0-9]{4,}"                                                                     # PubMed Central id
_ACCESSION_TOKEN = r"[A-Za-z0-9][A-Za-z0-9._-]*"  # a single non-space id for namespaces whose grammars are diverse

# Per-namespace grammar for ``namespace:local`` ids. The namespace prefix is
# matched case-insensitively; the local part is validated against the namespace's
# own grammar, so ``UniProt:Q9LZS3`` is citable but ``PMID:fixture`` and
# ``UniProt:candidate`` are not. Namespaces with heterogeneous id shapes
# (Ensembl, Phytozome) accept a permissive single token,
# which still rejects free text with whitespace.
_NAMESPACE_GRAMMARS = {
    "pmid": r"[0-9]+", "pmc": r"(?:PMC)?[0-9]+", "pmcid": r"(?:PMC)?[0-9]+",
    "doi": r"10\.[0-9]{4,}/\S+",
    "uniprot": _UNIPROT, "uniprotkb": _UNIPROT, "swissprot": _UNIPROT, "sprot": _UNIPROT, "trembl": _UNIPROT,
    "tair": _AGI, "araport": _AGI,
    "refseq": _REFSEQ, "ncbiprotein": rf"(?:{_REFSEQ}|{_INSDC})", "ncbi": rf"(?:{_REFSEQ}|{_INSDC})",
    "genbank": _INSDC, "embl": _INSDC, "ddbj": _INSDC,
    "ensembl": _ACCESSION_TOKEN, "ensemblplants": _ACCESSION_TOKEN, "phytozome": _ACCESSION_TOKEN,
    "gramene": _ACCESSION_TOKEN, "medtr": _ACCESSION_TOKEN, "medicago": _ACCESSION_TOKEN, "lotus": _ACCESSION_TOKEN,
    "panther": _ACCESSION_TOKEN, "interpro": r"IPR[0-9]{6}", "pfam": r"PF[0-9]{5}",
    "go": r"[0-9]{7}", "po": r"[0-9]{7}", "to": r"[0-9]{7}", "co": _ACCESSION_TOKEN,
}
# Ontology / literature tokens that are citable wherever they appear in the id.
_ONTOLOGY_OR_DOI = re.compile(r"(?i)(pmid:\s*\d|10\.\d{4,}/|\bGO:\d|\bTO:\d|\bPO:\d|\bCO_\d)")
_NAMESPACED = re.compile(r"^\s*([A-Za-z][A-Za-z0-9]*)\s*[:|]\s*(\S.*?)\s*$")
_BARE_ACCESSION = re.compile(rf"^(?:{_UNIPROT}|{_AGI}|{_REFSEQ}|{_INSDC}|{_PMC})$")  # whole uppercase token, case-sensitive
_FILE_SUFFIXES = (".tsv", ".csv", ".fa", ".fasta", ".fna", ".json", ".gff3", ".gff", ".vcf", ".paf", ".txt", ".tbl", ".emapper.annotations")


def _namespaced_reference(text: str) -> bool:
    """True when ``text`` is ``namespace:local`` with a known namespace whose grammar accepts ``local``."""
    match = _NAMESPACED.match(text)
    if not match:
        return False
    grammar = _NAMESPACE_GRAMMARS.get(match.group(1).lower())
    return bool(grammar and re.fullmatch(grammar, match.group(2), re.IGNORECASE))


def classify_source(source_id: object) -> str:
    """Grade a source_id: reference (citable id) / local_file / label / none.

    ``reference`` is allowlist-driven (issue #153):

    * ontology / literature tokens (``PMID:12345``, ``GO:0006355``, a DOI);
    * a namespaced database id whose local part matches that namespace's
      accession grammar (``UniProt:Q9LZS3``, ``TAIR:AT1G69120``,
      ``RefSeq:NP_001078807``) -- a known namespace with a malformed local part
      (``PMID:fixture``) is not citable;
    * a bare accession matched as a whole uppercase token against a specific,
      well-known grammar (UniProt, AGI, RefSeq, GenBank/EMBL/DDBJ, PMC).

    Unknown namespaces and unknown compact tokens (``gene1234``, ``trait2026``,
    ``abc_12345``, ``GENE0001``) grade ``label``: they cannot anchor an
    independence claim, so a free label never inflates independent_source_count
    or promotes a candidate to Top. Paths / known data suffixes grade
    ``local_file``; empty or placeholder values grade ``none``.
    """
    text = str(source_id).strip()
    if not text or text.lower() in {"na", "none", "required", "pending", "-"}:
        return "none"
    if _ONTOLOGY_OR_DOI.search(text) or _namespaced_reference(text):
        return "reference"
    if " " not in text and _BARE_ACCESSION.match(text):
        return "reference"
    lowered = text.lower()
    if "/" in text or "\\" in text or lowered.endswith(_FILE_SUFFIXES):
        return "local_file"
    return "label"


def independence_group(source_id: object, explicit_group: object = "") -> str:
    """Resolve the independence group for a support row.

    An explicit ``independence_group`` value wins (so different papers reusing one
    dataset can be tied together). Otherwise a citable source anchors its own
    group; labels and empty sources cannot prove independence and return ""."""
    explicit = str(explicit_group).strip()
    if explicit:
        return explicit
    kind = classify_source(source_id)
    return str(source_id).strip() if kind in {"reference", "local_file"} else ""


# Generic conservation (BLAST best hit, HMMER domain, GO presence) is imported
# as context, not support, so it never adds trait score; it is still tracked
# to retain conserved-but-uncharacterized candidates as Rescue. Rows whose only
# backing is such a pipeline artifact collapse to one shared group so several
# accessions from one analysis chain cannot masquerade as independent studies.
CONSERVATION_CONTEXT = {"homology_context", "domain_context", "go_annotation_context", "homology_supported", "curated_homolog_match", "automated_homolog_match"}


def support_group(source_id: object, explicit_group: object, code: str) -> str:
    """Resolve the independence group of an evidence row (public).

    The same resolver the scorer uses, exposed so ranking-qa source-dropout
    analysis drops rows by the exact group the scoring layer recognises:
    an explicit ``independence_group`` wins; uncited comparative-function rows
    collapse to ``comparative_annotation_chain``; otherwise a citable source
    anchors its own group and labels/empties resolve to "".
    """
    explicit = str(explicit_group).strip()
    if explicit:
        return explicit
    if code in CONSERVATION_CONTEXT:
        return "comparative_annotation_chain"
    if (
        INDEPENDENCE_FAMILIES.get(code) == "comparative_function"
        and classify_source(source_id) == "local_file"
    ):
        return "comparative_annotation_chain"
    return independence_group(source_id, "")


def annotation_keyword_evidence(
    candidates_arg: str,
    blast_arg: str,
    annotation_gff_arg: str | None,
    keyword_args: list[str],
    out_arg: str,
    min_identity: float = 30.0,
    min_query_coverage: float = 0.50,
    min_subject_coverage: float = 0.30,
    max_evalue: float = 1e-5,
    subject_annotations_arg: str | None = None,
    use_stitle: bool = False,
) -> dict[str, object]:
    """Link quality-gated protein hits to GFF descriptions and trait keywords."""
    if not 0 <= min_identity <= 100:
        raise ValueError("min_identity must be between 0 and 100")
    if not 0 <= min_query_coverage <= 1 or not 0 <= min_subject_coverage <= 1:
        raise ValueError("coverage thresholds must be between 0 and 1")
    if max_evalue <= 0:
        raise ValueError("max_evalue must be greater than 0")

    candidates = _candidate_frame(candidates_arg)
    mapping = _id_map(candidates)
    valid_genes = set(candidates["gene_id"])
    annotations: dict[str, tuple[str, str]] = {}
    annotation_path = runtime_path(annotation_gff_arg) if annotation_gff_arg else None
    from urllib.parse import unquote

    if annotation_path:
        opener = gzip.open if annotation_path.suffix == ".gz" else open
        with opener(annotation_path, "rt", encoding="utf-8") as handle:
            for raw in handle:
                if not raw.strip() or raw.startswith("#"):
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) != 9 or fields[2] != "gene":
                    continue
                attrs: dict[str, str] = {}
                for item in fields[8].split(";"):
                    if "=" in item:
                        key, value = item.split("=", 1)
                        attrs[unquote(key)] = unquote(value)
                identifier = attrs.get("ID") or attrs.get("Name")
                description = attrs.get("Note") or attrs.get("description") or attrs.get("Name", "")
                if identifier:
                    for alias in _annotation_aliases(identifier):
                        annotations[alias] = (identifier, description)
    if subject_annotations_arg:
        table = pd.read_csv(runtime_path(subject_annotations_arg), sep="\t", dtype=str, header=None, comment="#").fillna("")
        if table.shape[1] < 2:
            raise ValueError("subject annotations must have at least two TSV columns: subject_id and description")
        for subject, description in table.iloc[:, :2].itertuples(index=False, name=None):
            for alias in _annotation_aliases(subject):
                annotations[alias] = (subject, description)

    best: dict[str, dict[str, object]] = {}
    blast_path = runtime_path(blast_arg)
    with blast_path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 13:
                raise ValueError(
                    f"{blast_path}:{line_number}: annotation-keywords requires the 13-column candidate-qtl BLAST/DIAMOND outfmt"
                )
            gene_id = mapping.get(fields[0], fields[0])
            if gene_id not in valid_genes:
                continue
            identity = float(fields[2])
            align_length = int(fields[3])
            query_length = max(1, int(fields[4]))
            subject_length = max(1, int(fields[5]))
            evalue = float(fields[10])
            bitscore = float(fields[11])
            record = {
                "subject": fields[1],
                "identity": identity,
                "query_coverage": align_length / query_length,
                "subject_coverage": align_length / subject_length,
                "evalue": evalue,
                "bitscore": bitscore,
            }
            if use_stitle:
                if len(fields) < 14:
                    raise ValueError("--use-stitle requires stitle as column 14 after the 13-column candidate-qtl fields")
                annotations[fields[1]] = (fields[1], fields[13])
            if gene_id not in best or bitscore > float(best[gene_id]["bitscore"]):
                best[gene_id] = record

    keywords = sorted({term.strip().casefold() for value in keyword_args for term in value.split(",") if term.strip()})
    if not keywords:
        raise ValueError("at least one non-empty --keyword is required")

    rows: list[dict[str, str]] = []
    audit_rows: list[dict[str, object]] = []
    for gene_id, hit in best.items():
        annotation = next((annotations[alias] for alias in _annotation_aliases(str(hit["subject"])) if alias in annotations), None)
        target_id, description = annotation if annotation else ("", "")
        quality_pass = (
            float(hit["identity"]) >= min_identity
            and float(hit["query_coverage"]) >= min_query_coverage
            and float(hit["subject_coverage"]) >= min_subject_coverage
            and float(hit["evalue"]) <= max_evalue
        )
        lowered = description.casefold()
        matched = [term for term in keywords if re.search(rf"(?<![0-9a-z]){re.escape(term)}(?![0-9a-z])", lowered)] if quality_pass else []
        audit_rows.append(
            {
                "gene_id": gene_id,
                "best_subject": hit["subject"],
                "target_gene_id": target_id,
                "bitscore": hit["bitscore"],
                "identity": hit["identity"],
                "query_coverage": round(float(hit["query_coverage"]), 6),
                "subject_coverage": round(float(hit["subject_coverage"]), 6),
                "evalue": hit["evalue"],
                "quality_gate": "pass" if quality_pass else "fail",
                "matched_keywords": ";".join(matched),
                "description": description,
            }
        )
        if matched:
            rows.append(
                _row(
                    gene_id,
                    "trait_function_match",
                    "support",
                    "medium",
                    ";".join(matched),
                    "gff_annotation_keyword",
                    annotation_path or runtime_path(subject_annotations_arg or blast_arg),
                    (
                        f"Quality-gated best-hit annotation matched explicit trait keywords ({', '.join(matched)}): "
                        f"identity={hit['identity']}; qcov={float(hit['query_coverage']):.3f}; "
                        f"scov={float(hit['subject_coverage']):.3f}; evalue={hit['evalue']}; {description}"
                    ),
                )
            )

    output = runtime_path(out_arg)
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=EVIDENCE_COLUMNS).to_csv(output, sep="\t", index=False)
    audit_path = output.with_name(output.stem + "_audit.tsv")
    pd.DataFrame(audit_rows).to_csv(audit_path, sep="\t", index=False)
    return {
        "output": str(output),
        "audit": str(audit_path),
        "candidate_count": len(candidates),
        "annotated_best_hits": sum(bool(row["target_gene_id"]) for row in audit_rows),
        "quality_pass_count": sum(row["quality_gate"] == "pass" for row in audit_rows),
        "matched_candidates": len(rows),
        "keywords": keywords,
        "thresholds": {
            "min_identity": min_identity,
            "min_query_coverage": min_query_coverage,
            "min_subject_coverage": min_subject_coverage,
            "max_evalue": max_evalue,
        },
        "warning": "Keyword matches are functional hypotheses, not causal evidence; failed homology gates never contribute trait-function score.",
    }


def _strongest_by_direction(code_rows: pd.DataFrame, weight: float) -> list[dict[str, object]]:
    chosen: list[dict[str, object]] = []
    for direction, sign in (("support", 1.0), ("against", -1.0)):
        subset = code_rows[code_rows["direction"] == direction]
        if subset.empty:
            continue
        candidates: list[tuple[float, pd.Series]] = []
        for _, item in subset.iterrows():
            multiplier = CONFIDENCE_MULTIPLIERS.get(str(item["confidence"]), 0.0)
            candidates.append((sign * weight * multiplier, item))
        score, item = max(candidates, key=lambda pair: abs(pair[0]))
        chosen.append(
            {
                "contribution": round(score, 3),
                "direction": direction,
                "confidence": item["confidence"],
                "source_id": item["source_id"],
            }
        )
    return chosen


def prioritize_evidence(
    candidates_arg: str,
    evidence_arg: str,
    out_dir_arg: str,
    weights_arg: str | None = None,
    top_n: int = 7,
    reserve_n: int = 3,
    min_top_score: float = 2.0,
    rescue_n: int = 3,
    case_id: str = "",
) -> dict[str, object]:
    """Rank source-backed Dry evidence with conflict flags and a bounded Rescue batch."""
    if min(top_n, reserve_n, rescue_n) < 0 or min_top_score < 0:
        raise ValueError("top_n, reserve_n, rescue_n, and min_top_score must be zero or greater")
    candidates = _candidate_frame(candidates_arg)
    evidence = pd.read_csv(runtime_path(evidence_arg), sep="\t", dtype=str).fillna("")
    schemas.validate(evidence, schemas.EVIDENCE_SCHEMA)
    weights = load_weights(weights_arg)
    if weights_arg:
        supplied = pd.read_csv(runtime_path(weights_arg), sep="\t", dtype=str).fillna("")
        if "base_weight_or_tier" in supplied.columns:
            registry, registry_meta = load_weight_registry(weights_arg)
        else:
            registry = pd.DataFrame({"evidence_code": list(weights), "base_weight_or_tier": list(weights.values())})
            registry_meta = {"path": str(runtime_path(weights_arg)), "version": "custom-unversioned", "sha256": hashlib.sha256(runtime_path(weights_arg).read_bytes()).hexdigest()}
    else:
        registry, registry_meta = load_weight_registry()
    unscored_codes = sorted(set(evidence.loc[evidence["direction"].isin({"support", "against"}), "evidence_code"]).difference(weights))

    scored_rows: list[dict[str, object]] = []
    conflict_by_gene: dict[str, set[str]] = {}
    for gene_id, group in evidence.groupby("gene_id", sort=False):
        support_codes_in_gene: set[str] = set()
        counterevidence_against = False
        for code, code_rows in group.groupby("evidence_code", sort=False):
            if code not in weights:
                continue
            directions = set(code_rows["direction"]).intersection({"support", "against"})
            if directions == {"support", "against"}:
                conflict_by_gene.setdefault(gene_id, set()).add(code)
            if "support" in directions:
                support_codes_in_gene.add(code)
            if code == "counterevidence" and "against" in directions:
                counterevidence_against = True
            for chosen in _strongest_by_direction(code_rows, weights[code]):
                scored_rows.append({"gene_id": gene_id, "evidence_code": code, **chosen})
        if counterevidence_against:
            cross_code_support = sorted(support_codes_in_gene - {"counterevidence"})
            if cross_code_support:
                conflict_by_gene.setdefault(gene_id, set()).update(
                    f"{support_code}+counterevidence" for support_code in cross_code_support
                )

    contribution_columns = ["gene_id", "evidence_code", "direction", "contribution", "confidence", "source_id"]
    contributions = pd.DataFrame(scored_rows, columns=contribution_columns)
    if contributions.empty:
        contributions["independence_family"] = pd.Series(dtype=str)
        contributions["used_for_score"] = pd.Series(dtype=bool)
        totals = pd.Series(dtype=float)
    else:
        contributions["independence_family"] = contributions["evidence_code"].map(INDEPENDENCE_FAMILIES).fillna(
            contributions["evidence_code"]
        )
        contributions["used_for_score"] = False
        for _, family_rows in contributions.groupby(
            ["gene_id", "independence_family", "direction"], sort=False
        ):
            chosen_index = family_rows["contribution"].abs().idxmax()
            contributions.at[chosen_index, "used_for_score"] = True
        totals = contributions[contributions["used_for_score"]].groupby("gene_id")["contribution"].sum()
    trace_columns = ["evidence_code", "directness", "trait_specificity", "independence_rule", "rationale", "reference_type", "reference_id", "review_status", "registry_version"]
    available_trace = [column for column in trace_columns if column in registry.columns]
    if available_trace:
        contributions = contributions.merge(registry[available_trace], on="evidence_code", how="left", validate="many_to_one")
    ranked = candidates.copy()
    ranked["dry_score"] = ranked["gene_id"].map(totals).fillna(0.0).round(3)
    ranked["screening_rank"] = ranked["dry_score"].rank(method="min", ascending=False).astype(int)
    ranked["screening_rank_max"] = ranked["dry_score"].rank(method="max", ascending=False).astype(int)
    ranked["screening_tie_size"] = ranked.groupby("dry_score")["gene_id"].transform("size").astype(int)
    if contributions.empty:
        breakdowns = pd.Series(dtype=str)
    else:
        used_contributions = contributions[contributions["used_for_score"]]
        breakdowns = pd.Series(
            {
                gene_id: ";".join(f"{row.evidence_code}/{row.direction}:{row.contribution:+g}" for row in group.itertuples())
                for gene_id, group in used_contributions.groupby("gene_id")
            },
            dtype=str,
        )
    ranked["score_breakdown"] = ranked["gene_id"].map(breakdowns).fillna("no_scored_evidence:0")
    ranked["evidence_conflicts"] = ranked["gene_id"].map(
        lambda gene_id: ";".join(sorted(conflict_by_gene.get(gene_id, set())))
    )

    rescue_codes = evidence[evidence["direction"] == "rescue"].groupby("gene_id")["evidence_code"].apply(lambda values: ";".join(sorted(set(values))))
    support_codes = evidence[evidence["direction"] == "support"].groupby("gene_id")["evidence_code"].apply(set)
    all_codes = evidence.groupby("gene_id")["evidence_code"].apply(set)

    # Distinct *citable* provenance per gene. support_family_count only measures
    # non-redundant scoring channels; it does not prove that two codes came from
    # independent studies/analyses. Independence needs distinct citable sources
    # (or an explicit independence_group), so a single paper's annotation +
    # expression is one group, not two.
    support_rows = evidence[evidence["direction"] == "support"].copy()
    explicit_groups = support_rows["independence_group"] if "independence_group" in support_rows.columns else pd.Series("", index=support_rows.index)

    # Comparative-function evidence (homology / domain / orthology / GO context)
    # for one gene comes from a single analysis chain, so several accessions from
    # that chain must collapse to one independent source -- otherwise recognising
    # bare UniProt/AGI accessions as citable (below) would let multiple hits from
    # one BLAST run masquerade as independent studies. An explicit
    # independence_group still wins, and genetic / position / expression evidence
    # with distinct citable sources keeps anchoring its own group.
    #
    # The same holds for any other comparative-function row whose only backing is
    # a pipeline artifact. ``annotation-keywords`` emits trait_function_match with
    # the *annotation file* as source_id, and that file holds the annotation of
    # the very best hit that produced the homology row -- one blastp run, one
    # chain. Grading that file ``local_file`` let it anchor a second independence
    # group, so a single BLAST run reached ``independent_experimental_support``
    # and, with any positional evidence on top, Top. A file path is not an
    # independent study, so an uncited comparative row joins the comparative
    # chain. A comparative row backed by a real citable reference (a curated
    # PMID/DOI or a resolved accession) still anchors its own group, which is what
    # keeps separately sourced trait-function and domain evidence independent.
    # ``support_group`` is the module-level resolver shared with ranking-qa.
    support_rows["_indep_group"] = [
        support_group(sid, grp, code)
        for sid, grp, code in zip(support_rows["source_id"], explicit_groups, support_rows["evidence_code"])
    ]
    independent_groups = support_rows[support_rows["_indep_group"] != ""].groupby("gene_id")["_indep_group"].apply(set)
    citable_support_codes = support_rows[support_rows["_indep_group"] != ""].groupby("gene_id")["evidence_code"].apply(set)

    def _channels(gene_id: str) -> set[str]:
        return {INDEPENDENCE_FAMILIES[code] for code in support_codes.get(gene_id, set()) if code in INDEPENDENCE_FAMILIES}

    def _independent_count(gene_id: str) -> int:
        return len(independent_groups.get(gene_id, set()))

    ranked["rescue_reason"] = ranked["gene_id"].map(rescue_codes).fillna("")
    ranked["support_family_count"] = ranked["gene_id"].map(lambda gene_id: len(_channels(gene_id)))
    ranked["independent_source_count"] = ranked["gene_id"].map(_independent_count)
    ranked["top_eligibility"] = ranked["gene_id"].map(
        lambda gene_id: (
            "conflicting_evidence"
            if conflict_by_gene.get(gene_id)
            else "strong_genetic_evidence"
            if citable_support_codes.get(gene_id, set()).intersection(STRONG_SINGLE_CODES)
            # >=2 distinct citable sources = genuinely independent experimental support.
            else "independent_experimental_support"
            if _independent_count(gene_id) >= 2
            # >=2 scoring channels but not independently sourced: non-redundant, not independent.
            else "multiple_non_redundant_channels"
            if len(_channels(gene_id)) >= 2
            else "insufficient_independent_evidence"
        )
    )

    for index, row in ranked.iterrows():
        codes = support_codes.get(row["gene_id"], set())
        gene_all_codes = all_codes.get(row["gene_id"], set())
        if gene_all_codes.intersection(CONSERVATION_CONTEXT) and not codes.intersection({"trait_function_match", "automated_candidate_function_match", "known_qtg", "direct_model_genetic", "ortholog_model_genetic", "orthogroup_supported"}):
            ranked.at[index, "rescue_reason"] = ";".join(filter(None, [ranked.at[index, "rescue_reason"], "conserved_uncharacterized"]))
        if "expression_target_specific" in codes and float(row["dry_score"]) < min_top_score:
            ranked.at[index, "rescue_reason"] = ";".join(filter(None, [ranked.at[index, "rescue_reason"], "target_expression_rescue"]))
        if "peak_proximal" in codes and float(row["dry_score"]) < min_top_score:
            ranked.at[index, "rescue_reason"] = ";".join(filter(None, [ranked.at[index, "rescue_reason"], "peak_proximity_rescue"]))
        if "recombinant_core" in codes and float(row["dry_score"]) < min_top_score:
            ranked.at[index, "rescue_reason"] = ";".join(filter(None, [ranked.at[index, "rescue_reason"], "recombinant_core_rescue"]))

    def rescue_priority(gene_id: str, reasons: str) -> float:
        codes = support_codes.get(gene_id, set())
        score = 0.0
        if "annotation_incomplete" in reasons:
            score += 3
        if "target_expression_rescue" in reasons:
            score += 2
        if "peak_proximity_rescue" in reasons:
            score += 2
        if "recombinant_core_rescue" in reasons:
            score += 2
        if "conserved_uncharacterized" in reasons:
            score += 1
        if "synteny_supported" in codes:
            score += 1
        return score

    ranked["rescue_priority"] = [rescue_priority(str(row.gene_id), str(row.rescue_reason)) for row in ranked.itertuples()]
    ranked = ranked.sort_values(["dry_score", "gene_id"], ascending=[False, True], kind="stable")
    ranked["dry_group"] = "evidence_insufficient"
    scorable = (ranked["dry_score"] >= min_top_score) & ranked["rescue_reason"].eq("")
    # Top requires genuine independence (strong genetic, or >=2 citable sources).
    top_eligible = ranked.index[scorable & ranked["top_eligibility"].isin({"strong_genetic_evidence", "independent_experimental_support"})]
    top_selected = set(top_eligible[:top_n])
    ranked.loc[top_eligible[:top_n], "dry_group"] = "top"
    # Reserve may include non-redundant multi-channel genes that are not yet
    # independently sourced -- worth reading if Top fails, but not Top.
    reserve_eligible = ranked.index[
        scorable
        & ranked["top_eligibility"].isin({"strong_genetic_evidence", "independent_experimental_support", "multiple_non_redundant_channels"})
        & ~ranked.index.isin(top_selected)
    ]
    ranked.loc[reserve_eligible[:reserve_n], "dry_group"] = "reserve"

    rescue_candidates = ranked[ranked["rescue_reason"].ne("")].sort_values(
        ["rescue_priority", "dry_score", "gene_id"], ascending=[False, False, True], kind="stable"
    )
    selected_rescue = rescue_candidates.head(rescue_n).index
    ranked.loc[selected_rescue, "dry_group"] = "rescue"
    ranked.loc[ranked["rescue_reason"].ne("") & ~ranked.index.isin(selected_rescue), "dry_group"] = "evidence_insufficient"

    order = {"top": 0, "reserve": 1, "rescue": 2, "evidence_insufficient": 3}
    ranked["_order"] = ranked["dry_group"].map(order)
    ranked = ranked.sort_values(["_order", "dry_score", "rescue_priority", "gene_id"], ascending=[True, False, False, True], kind="stable")
    ranked["dry_rank"] = ranked.groupby("dry_group", sort=False).cumcount() + 1
    ranked = ranked.drop(columns="_order")

    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    ranked_path = out_dir / "ranked_candidates.tsv"
    contribution_path = out_dir / "score_contributions.tsv"
    weights_path = out_dir / "weights_used.tsv"
    ranked.to_csv(ranked_path, sep="\t", index=False)
    contributions.to_csv(contribution_path, sep="\t", index=False)
    weights_used = registry[registry["evidence_code"].isin(weights)].copy() if "review_status" in registry.columns else registry.copy()
    weights_used["weight"] = weights_used["evidence_code"].map(weights)
    weights_used.to_csv(weights_path, sep="\t", index=False)

    top = ranked[ranked["dry_group"] == "top"]
    functional_codes = {"direct_model_genetic", "ortholog_model_genetic", "known_qtg", "trait_function_match"}
    top_with_function = sum(bool(support_codes.get(gene_id, set()).intersection(functional_codes)) for gene_id in top["gene_id"])
    fallback_reasons: list[str] = []
    alerts: list[str] = []
    if unscored_codes:
        alerts.append("unweighted_evidence_codes:" + ",".join(unscored_codes))
    if conflict_by_gene:
        alerts.append("evidence_conflicts_require_review:" + ",".join(sorted(conflict_by_gene)))
    if len(top) < top_n:
        fallback_reasons.append("fewer_supported_top_candidates_than_requested")
    if top_with_function == 0:
        fallback_reasons.append("no_top_candidate_has_functional_or_genetic_support")
    if (ranked.get("annotation_quality", pd.Series("", index=ranked.index)) == "no_cds").any():
        alerts.append("interval_contains_gene_models_without_cds")
    decision = {
        "case_id": case_id,
        "mode": "broaden_to_regulatory_sv_te_cnv" if fallback_reasons else "test_coding_candidates_first",
        "fallback_reasons": fallback_reasons,
        "alerts": alerts,
        "top_candidate_count": len(top),
        "rescue_candidate_count": int((ranked["dry_group"] == "rescue").sum()),
        "rescue_limit": rescue_n,
        "top_with_functional_or_genetic_support": top_with_function,
        "top_with_independent_experimental_support": int((top["top_eligibility"] == "independent_experimental_support").sum()) if not top.empty else 0,
        "reserve_multiple_non_redundant_channels": int((ranked[ranked["dry_group"] == "reserve"]["top_eligibility"] == "multiple_non_redundant_channels").sum()),
        "independence_note": (
            "support_family_count counts non-redundant scoring channels (a dedup guard), NOT independent studies. "
            "independent_source_count counts distinct citable sources; Top requires strong genetic evidence or "
            "independent_experimental_support (>=2 citable sources). multiple_non_redundant_channels from one source is Reserve, not Top."
        ),
        "interpretation": "This is an early branch decision, not a claim that regulatory or structural variation is causal.",
        "screening_rank_interpretation": (
            "screening_rank orders functional suspicion before Wet work; dry_group separately records whether independent evidence is sufficient."
        ),
        "weight_registry": registry_meta,
        "weight_interpretation": "Ordinal pre-Wet assay-prioritization heuristic; not a causal probability, posterior, or statistical significance.",
    }
    decision_path = out_dir / "candidate_path_decision.json"
    decision_path.write_text(json.dumps(decision, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return {
        "ranked_candidates": str(ranked_path),
        "score_contributions": str(contribution_path),
        "weights_used": str(weights_path),
        "decision": str(decision_path),
        "groups": {key: int(value) for key, value in ranked["dry_group"].value_counts().items()},
        "weight_registry": registry_meta,
    }


def candidate_pipeline(
    candidates_arg: str,
    out_dir_arg: str,
    blast_args: list[str] | None = None,
    diamond_args: list[str] | None = None,
    hmmer_args: list[str] | None = None,
    eggnog_args: list[str] | None = None,
    expression_arg: str | None = None,
    manual_evidence_arg: str | None = None,
    refined_evidence_arg: str | None = None,
    target_tissues: list[str] | None = None,
    target_stages: list[str] | None = None,
    min_expression: float = 1.0,
    target_fold: float = 2.0,
    peak: str | None = None,
    interval: str | None = None,
    peak_window: int = 0,
    boundary_window: int = 0,
    weights_arg: str | None = None,
    top_n: int = 7,
    reserve_n: int = 3,
    min_top_score: float = 2.0,
    rescue_n: int = 3,
    protein_arg: str | None = None,
    target_proteome_arg: str | None = None,
    relative_proteome_args: list[str] | None = None,
    threads: int = 4,
) -> dict[str, object]:
    out_dir = runtime_path(out_dir_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    evidence_result = evidence_build(
        candidates_arg,
        str(out_dir / "01_evidence"),
        blast_args,
        diamond_args,
        hmmer_args,
        eggnog_args,
        expression_arg,
        manual_evidence_arg,
        None,
        target_tissues,
        target_stages,
        min_expression,
        target_fold,
        peak,
        interval,
        peak_window,
        boundary_window,
    )
    primary_rank = prioritize_evidence(
        candidates_arg,
        str(evidence_result["evidence"]),
        str(out_dir / "02_primary_rank"),
        weights_arg,
        top_n,
        reserve_n,
        min_top_score,
        rescue_n,
    )
    stages: list[dict[str, object]] = [
        {"stage": 1, "name": "structured evidence", "status": "complete", "result": evidence_result},
        {"stage": 2, "name": "primary rank", "status": "complete", "result": primary_rank},
    ]
    orthology_result = None
    if protein_arg and target_proteome_arg and relative_proteome_args:
        orthology_result = orthology_plan(
            str(primary_rank["ranked_candidates"]), protein_arg, target_proteome_arg, relative_proteome_args,
            str(out_dir / "03_orthology_refinement"), top_n, threads,
        )
        stages.append({"stage": 3, "name": "top-candidate orthology refinement", "status": "planned; run script, then orthology-import", "result": orthology_result})
    else:
        stages.append({"stage": 3, "name": "top-candidate orthology refinement", "status": "not requested", "requires": ["--protein-fasta", "--target-proteome", "--relative-proteome"]})

    final_rank = primary_rank
    if refined_evidence_arg:
        refined_result = evidence_build(
            candidates_arg,
            str(out_dir / "04_refined_evidence"),
            existing_evidence_arg=str(evidence_result["evidence"]),
            manual_evidence_arg=refined_evidence_arg,
        )
        final_rank = prioritize_evidence(
            candidates_arg,
            str(refined_result["evidence"]),
            str(out_dir / "05_final_rank"),
            weights_arg,
            top_n,
            reserve_n,
            min_top_score,
            rescue_n,
        )
        stages.extend([
            {"stage": 4, "name": "refined evidence merge", "status": "complete", "result": refined_result},
            {"stage": 5, "name": "final rank", "status": "complete", "result": final_rank},
        ])

    manifest = {
        "candidates": str(runtime_path(candidates_arg)),
        "stages": stages,
        "primary_ranked_candidates": primary_rank["ranked_candidates"],
        "final_ranked_candidates": final_rank["ranked_candidates"],
        "decision": final_rank["decision"],
        "next_action": (
            "Run the orthology refinement script and orthology-import, then rerun candidate-pipeline with --refined-evidence."
            if orthology_result and not refined_evidence_arg
            else "Review Top/Rescue candidates and the coding-vs-regulatory/SV decision before parental sequencing."
        ),
    }
    manifest_path = out_dir / "candidate_pipeline.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return {
        "manifest": str(manifest_path),
        "primary_ranked_candidates": primary_rank["ranked_candidates"],
        "final_ranked_candidates": final_rank["ranked_candidates"],
        "decision": final_rank["decision"],
        "orthology_plan": orthology_result["manifest"] if orthology_result else None,
        "next_action": manifest["next_action"],
    }


def _select_ranked_candidates(frame: pd.DataFrame, top_n: int) -> pd.DataFrame:
    selected = frame.copy()
    if "deterministic_review_order" in selected.columns:
        # A frozen Stage-A shortlist must be passed to humans and models in the
        # same order.  Eligibility groups are interpretation, not a second
        # candidate filter at this boundary.
        selected["_group"] = 0
        selected["_rank"] = pd.to_numeric(
            selected["deterministic_review_order"], errors="coerce"
        ).fillna(999999)
        return selected.sort_values(["_rank", "gene_id"], kind="stable").head(top_n)
    if "dry_group" in selected.columns:
        order = {"top": 0, "reserve": 1, "rescue": 2, "evidence_insufficient": 3, "not_selected_yet": 3}
        selected["_group"] = selected["dry_group"].str.casefold().map(order).fillna(9)
    else:
        selected["_group"] = 0
    for column in ("dry_rank", "program_rank", "rank"):
        if column in selected.columns:
            selected["_rank"] = pd.to_numeric(selected[column], errors="coerce").fillna(999999)
            break
    else:
        selected["_rank"] = range(1, len(selected) + 1)
    return selected.sort_values(["_group", "_rank", "gene_id"], kind="stable").head(top_n)


def _candidate_literature_query(evidence_rows: list[dict[str, str]], trait: str) -> str:
    """Build a literature query from public homolog/function terms, not local IDs."""
    aliases: list[str] = []
    descriptions: list[str] = []
    for item in evidence_rows:
        if item.get("direction") not in {"support", "against"}:
            continue
        value = str(item.get("value", "")).strip()
        if not value:
            continue
        match = re.search(r"\bGN=([A-Za-z0-9_-]+)", value)
        if match:
            aliases.append(match.group(1))
        description = re.split(r"\s+OS=", value, maxsplit=1)[0].strip()
        if description and not description.casefold().startswith(("chr", "scaffold")):
            descriptions.append(description[:100])
        for token in re.findall(r"\b[A-Z][A-Z0-9-]{1,11}\b", value):
            if token not in {"OS", "OX", "GN", "PE", "SV", "LOF"}:
                aliases.append(token)
    aliases = list(dict.fromkeys(aliases))[:4]
    descriptions = list(dict.fromkeys(descriptions))[:2]
    trait_folded = trait.casefold()
    if "powdery mildew" in trait_folded:
        trait_terms = ['"powdery mildew"', "resistance"]
    elif "fasciat" in trait_folded:
        trait_terms = ["fasciation", "fasciated", '"shoot apical meristem"']
    else:
        cleaned = re.sub(r"\([^)]*\)", "", trait).strip()
        trait_terms = [f'"{cleaned}"'] if cleaned else []
    function_terms = aliases + [f'"{value}"' for value in descriptions]
    left = " OR ".join(trait_terms)
    right = " OR ".join(function_terms)
    if left and right:
        return f"({left}) AND ({right})"
    return left or right or f'"{trait}"'


def research_packet(
    candidate_arg: str,
    out_arg: str,
    species: str,
    trait: str,
    top_n: int,
    fetch_literature: bool,
    literature_limit: int,
    evidence_arg: str | None = None,
    contributions_arg: str | None = None,
) -> dict[str, object]:
    """Create a machine-readable OSS-to-LLM packet with source-backed evidence."""
    candidate_path, out_dir = runtime_path(candidate_arg), runtime_path(out_arg)
    out_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(candidate_path, sep="\t", dtype=str).fillna("")
    if "gene_id" not in frame.columns:
        raise ValueError("candidate table must contain gene_id")
    selected = _select_ranked_candidates(frame, top_n)
    selected_ids = set(selected["gene_id"])

    evidence_by_gene: dict[str, list[dict[str, str]]] = {gene_id: [] for gene_id in selected_ids}
    if evidence_arg:
        evidence = pd.read_csv(runtime_path(evidence_arg), sep="\t", dtype=str).fillna("")
        missing = sorted(set(EVIDENCE_COLUMNS).difference(evidence.columns))
        if missing:
            raise ValueError(f"evidence table is missing required columns: {', '.join(missing)}")
        for gene_id, group in evidence[evidence["gene_id"].isin(selected_ids)].groupby("gene_id"):
            evidence_by_gene[gene_id] = group[EVIDENCE_COLUMNS].to_dict("records")

    contributions_by_gene: dict[str, list[dict[str, str]]] = {gene_id: [] for gene_id in selected_ids}
    if contributions_arg:
        contributions = pd.read_csv(runtime_path(contributions_arg), sep="\t", dtype=str).fillna("")
        required = {"gene_id", "evidence_code", "contribution", "source_id"}
        if not required.issubset(contributions.columns):
            raise ValueError("score contributions require gene_id, evidence_code, contribution, and source_id")
        for gene_id, group in contributions[contributions["gene_id"].isin(selected_ids)].groupby("gene_id"):
            contributions_by_gene[gene_id] = group.to_dict("records")

    packets: list[dict[str, object]] = []
    prompts = [
        "# Candidate-gene evidence review",
        "",
        "Use only cited sources. Separate support, counterevidence, uncertainty, and recommended Wet tests. Do not infer causality from functional plausibility.",
    ]
    literature_files: list[str] = []
    for _, row in selected.iterrows():
        gene_id = str(row["gene_id"])
        evidence_rows = evidence_by_gene.get(gene_id, [])
        query_terms = _candidate_literature_query(evidence_rows, trait)
        literature: list[dict[str, object]] = []
        if fetch_literature:
            literature_path = out_dir / f"{safe_id(gene_id)}_europe_pmc.json"
            europe_pmc_search(query_terms, str(literature_path), literature_limit)
            literature = json.loads(literature_path.read_text(encoding="utf-8")).get("records", [])
            literature_files.append(str(literature_path))
        packet = {
            "gene_id": gene_id,
            "species": species,
            "trait": trait,
            "candidate": {key: str(value) for key, value in row.items() if not str(key).startswith("_")},
            "evidence": evidence_rows,
            "score_contributions": contributions_by_gene.get(gene_id, []),
            "literature_query": query_terms,
            "literature": literature,
            "uniprot_search_url": "https://rest.uniprot.org/uniprotkb/search?" + urlencode({"query": query_terms, "format": "json", "size": 10}),
            "required_output_schema": {
                "gene_id": gene_id,
                "claims": [
                    {
                        "evidence_code": "one supported canonical code",
                        "direction": "support|against|context",
                        "confidence": "high|medium|low",
                        "value": "short value",
                        "source_id": "PMID, DOI, accession, or stable URL",
                        "description": "claim tied to the source",
                        "reviewed": True,
                    }
                ],
            },
        }
        packets.append(packet)
        prompts.extend([
            "",
            f"## {gene_id}",
            f"Species: {species}",
            f"Trait: {trait}",
            f"Packet: {safe_id(gene_id)}.json",
            "Return JSON matching required_output_schema. Claims without a stable source_id or requiring speculation must be omitted.",
        ])
        (out_dir / f"{safe_id(gene_id)}.json").write_text(json.dumps(packet, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    packet_index = out_dir / "research_packet.json"
    packet_index.write_text(json.dumps({"species": species, "trait": trait, "candidates": packets}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    selected_prompts = out_dir / "ai_research_prompts.md"
    selected_prompts.write_text("\n".join(prompts) + "\n", encoding="utf-8")
    return {
        "packet": str(packet_index),
        "prompts": str(selected_prompts),
        "candidate_count": len(packets),
        "literature_files": literature_files,
        "ai_optional": True,
        "warning": "LLM output is not evidence until llm-evidence-import validates reviewed=true and a stable source_id.",
    }


# Every LLM-imported row whose source_id was not resolved against a real corpus
# shares this independence group. `support_group` honours an explicit group over
# the source_id, so N unverified accessions collapse to one group and can never
# reach `independent_source_count >= 2` on their own. A well-formed accession is
# not evidence that the accession exists.
UNVERIFIED_LLM_SOURCE_GROUP = "llm_curation_source_unverified"

# Token characters in the accession grammars classify_source recognises. Used to
# harvest candidate identifiers out of whatever corpus files the reviewer points
# at (FASTA headers, TSV columns, plain accession lists) without needing to know
# their layout.
_CORPUS_TOKEN = re.compile(r"[A-Za-z0-9_.:\-]+")


def _load_source_corpus(paths: list[str] | None) -> tuple[set[str], list[dict[str, object]]]:
    """Collect the identifiers that a source_id may resolve against.

    Deliberately format-agnostic: a reviewer should be able to point at the
    reference proteome, an annotation table, and a citation list without
    converting them first. Tokens are harvested, plus the local part of any
    ``namespace:local`` token so ``UniProt:Q9LZS3`` in a claim matches a bare
    ``Q9LZS3`` in a FASTA header.
    """
    corpus: set[str] = set()
    manifest: list[dict[str, object]] = []
    for raw in paths or []:
        path = runtime_path(raw)
        if not path.exists():
            raise ValueError(f"--resolve-against file does not exist: {path}")
        before = len(corpus)
        with open_text(path) as handle:
            for line in handle:
                for token in _CORPUS_TOKEN.findall(line):
                    corpus.add(token)
                    corpus.add(token.casefold())
                    if ":" in token:
                        local = token.split(":", 1)[1]
                        if local:
                            corpus.add(local)
                            corpus.add(local.casefold())
        manifest.append({"path": str(path), "fingerprint": file_fingerprint(path, False), "new_tokens": len(corpus) - before})
    return corpus, manifest


def _resolves(source_id: str, corpus: set[str]) -> bool:
    candidates = {source_id, source_id.casefold()}
    if ":" in source_id:
        local = source_id.split(":", 1)[1]
        candidates.update({local, local.casefold()})
    return bool(candidates & corpus)


def llm_evidence_import(
    candidate_arg: str,
    input_arg: str,
    out_arg: str,
    resolve_against: list[str] | None = None,
) -> dict[str, object]:
    """Import only human-reviewed, source-backed LLM research claims.

    Two gates beyond the original format checks (#189):

    * ``reviewed: true`` must be attributed. It is a self-declaration inside a
      file a model can write, so it proves nothing on its own; requiring
      ``reviewed_by`` at least makes the claim attributable and auditable.
    * A source_id that merely *looks* like an accession no longer anchors an
      independence group. With ``--resolve-against`` it must be found in a real
      local corpus; without one, accepted rows are marked unverified and share a
      single independence group, so fabricated-but-well-formed accessions cannot
      manufacture ``independent_experimental_support``.
    """
    candidates = _candidate_frame(candidate_arg)
    valid_genes = set(candidates["gene_id"])
    payload = json.loads(runtime_path(input_arg).read_text(encoding="utf-8"))
    records = payload.get("candidates", payload) if isinstance(payload, dict) else payload
    if not isinstance(records, list):
        raise ValueError("LLM evidence JSON must be a list or an object with a candidates list")
    corpus, corpus_manifest = _load_source_corpus(resolve_against)
    verification_performed = bool(corpus_manifest)
    rows: list[dict[str, str]] = []
    rejected: list[dict[str, str]] = []
    audit: list[dict[str, object]] = []
    for record in records:
        if not isinstance(record, dict):
            rejected.append({"gene_id": "", "reason": "record_not_object"})
            continue
        gene_id = str(record.get("gene_id", ""))
        claims = record.get("claims", [])
        if gene_id not in valid_genes or not isinstance(claims, list):
            rejected.append({"gene_id": gene_id, "reason": "unknown_gene_or_claims_not_list"})
            continue
        for claim in claims:
            if not isinstance(claim, dict):
                rejected.append({"gene_id": gene_id, "reason": "claim_not_object"})
                continue
            code = str(claim.get("evidence_code", ""))
            direction = str(claim.get("direction", ""))
            confidence = str(claim.get("confidence", ""))
            source_id = str(claim.get("source_id", "")).strip()
            reviewed = claim.get("reviewed") is True
            reviewed_by = str(claim.get("reviewed_by", "")).strip()
            reviewed_at = str(claim.get("reviewed_at", "")).strip()
            resolved = _resolves(source_id, corpus) if source_id else False
            reason = ""
            if not reviewed:
                reason = "not_human_reviewed"
            elif not reviewed_by:
                # `reviewed: true` alone is a bare assertion in a file the model
                # may have authored end to end. Attribution does not prove a human
                # reviewed it, but it makes the claim answerable to someone.
                reason = "missing_reviewer_attribution"
            elif code not in LLM_IMPORT_CODES:
                reason = "unsupported_evidence_code"
            elif direction not in {"support", "against", "context"}:
                reason = "unsupported_direction"
            elif confidence not in CONFIDENCE_MULTIPLIERS:
                reason = "unsupported_confidence"
            elif not source_id or source_id.casefold() in {"none", "unknown", "n/a"}:
                reason = "missing_source_id"
            elif verification_performed and not resolved:
                reason = "unresolvable_source_id"
            if reason:
                rejected.append({"gene_id": gene_id, "evidence_code": code, "source_id": source_id, "reason": reason})
                continue
            row = _row(
                gene_id,
                code,
                direction,
                confidence,
                claim.get("value", ""),
                "llm_reviewed_curation",
                source_id,
                str(claim.get("description", "")),
            )
            row["independence_group"] = "" if resolved else UNVERIFIED_LLM_SOURCE_GROUP
            rows.append(row)
            audit.append(
                {
                    "gene_id": gene_id,
                    "evidence_code": code,
                    "source_id": source_id,
                    "reviewed_by": reviewed_by,
                    "reviewed_at": reviewed_at or "not_recorded",
                    "review_attestation": "self_declared_not_verified",
                    "source_verification": "resolved_in_local_corpus" if resolved else ("unresolved" if verification_performed else "not_performed"),
                    "independence_group": row["independence_group"],
                }
            )
    output = runtime_path(out_arg)
    output.parent.mkdir(parents=True, exist_ok=True)
    columns = [*EVIDENCE_COLUMNS, "independence_group"]
    pd.DataFrame(rows, columns=columns).drop_duplicates().to_csv(output, sep="\t", index=False)
    rejected_path = output.with_name(output.stem + "_rejected.tsv")
    pd.DataFrame(rejected).to_csv(rejected_path, sep="\t", index=False)
    audit_path = output.with_name(output.stem + "_review_audit.tsv")
    pd.DataFrame(audit).to_csv(audit_path, sep="\t", index=False)
    unverified = sum(1 for row in rows if row["independence_group"] == UNVERIFIED_LLM_SOURCE_GROUP)
    return {
        "output": str(output),
        "accepted_count": len(rows),
        "rejected_count": len(rejected),
        "rejected": str(rejected_path),
        "review_audit": str(audit_path),
        "source_verification": "performed" if verification_performed else "not_performed",
        "resolve_against": corpus_manifest,
        "resolved_source_count": len(rows) - unverified,
        "unverified_source_count": unverified,
        # Recorded, not silently assumed: a run with no corpus is not the same as
        # a run whose sources were checked, and the manifest must show which.
        "human_review_attestation": "self_declared_not_verified",
        "warning": (
            "Accepted rows remain curated hypotheses; merge them into canonical evidence and inspect conflicts before ranking. "
            "reviewed=true is a self-declaration, not proof of human review. "
            + (
                f"{unverified} row(s) have unverified sources and share independence group '{UNVERIFIED_LLM_SOURCE_GROUP}', so they cannot reach independent_experimental_support on their own; pass --resolve-against to verify them."
                if unverified
                else "All accepted sources resolved against the supplied corpus."
            )
        ),
    }
