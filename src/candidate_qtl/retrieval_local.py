"""Workstation-local Route C provider: frozen plant reference indexes (issue #217).

``retrieval-import`` / ``retrieval-merge`` / ``retrieval-interval-match`` define a
strict Route C contract, but nothing in the repository could *produce* a
retrieval on this workstation: the monitored audit found no configured executor,
no model, and no indexed reference database, so Route C was an import interface
rather than an executable route and every case came back
``not_evaluated_provider_unavailable``.

This module makes Route C runnable locally, with two commands:

``retrieval-index-build``
    Turn one or more locally cached, source-backed reference corpora into an
    immutable, versioned index: a frozen document table with row-level
    provenance, a token->posting index, an exclusion ledger for rows that carry
    no family/orthogroup or no usable text, and a manifest with SHA-256 for
    every artifact plus the corpus fingerprints, preprocessing description,
    model name/version and build timestamp.

``retrieval-run``
    Phenotype text in, Route C retrieval contract out. **The command has no
    interval argument of any kind** -- not an interval, not an atlas, not a
    candidate table -- so Phase 1 cannot receive interval information even by
    mistake. It re-verifies the index fingerprints before retrieving, records
    model/database versions and hashes, and writes a frozen TSV that
    ``retrieval-import`` accepts unchanged.

**Provider choice.** The default provider is a deterministic local lexical
TF-IDF scorer with no third-party dependency. That is a deliberate first
version, not a stand-in for an embedding model: it runs on this workstation
today, it is reproducible bit-for-bit, and it never needs to send a private
phenotype anywhere. Whether a neural embedding backend is installed is audited
and recorded in every manifest, so upgrading to one is a provider swap with a
new ``model_version`` and a new index, not a rewrite. Nothing about the honesty
rules changes with the provider: retrieval similarity stays a discovery route,
never causal evidence, and ``rank_contribution`` stays ``none``.

**Failure is not a miss.** A missing index, an unverifiable fingerprint, or a
query that matches nothing is reported as an explicit status, never as a
biological result.
"""

from __future__ import annotations

import importlib.util
import json
import math
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from . import __version__
from .core import file_fingerprint, runtime_path, safe_id, sha256
from .phenotype_context import structured_context_block
from .retrieval_route import RETRIEVAL_COLUMNS

SCHEMA_VERSION = "1.0"
INDEX_FORMAT_VERSION = "1.0.0"
# The scorer is versioned independently of the package: a change to tokenisation
# or weighting must invalidate an existing index, which a reader can only see if
# the model version moves.
MODEL_NAME = "candidate-qtl-local-lexical"
# 1.1.0 indexes the curated phenotype/tissue/stage/family text that issue #232
# added to the corpus contract. The indexed field set is part of the scorer, so
# an index built by 1.0.0 is not comparable and must be rebuilt.
MODEL_VERSION = "tfidf-1.1.0"
MODALITY = "text_lexical"

CORPUS_REQUIRED_COLUMNS = ["reference_gene_id", "reference_species", "source_id"]
CORPUS_TEXT_COLUMNS = ("description", "function_description", "mutant_phenotype", "pathway", "go_terms", "keywords",
                       "tissue_specificity", "developmental_stage", "protein_family")
CORPUS_FAMILY_COLUMNS = ("family_id", "orthogroup_id")
DOCUMENT_COLUMNS = [
    "doc_id", "reference_gene_id", "reference_species", "family_id", "orthogroup_id",
    "document_text", "source_id", "source_url", "corpus_file", "corpus_row",
]
EXCLUSION_COLUMNS = ["reference_gene_id", "corpus_file", "corpus_row", "exclusion_reason"]

DOCUMENTS_NAME = "documents.tsv"
POSTINGS_NAME = "postings.json"
EXCLUDED_NAME = "excluded_rows.tsv"
INDEX_MANIFEST_NAME = "index_manifest.json"
RETRIEVAL_NAME = "route_c_local_retrieval.tsv"
RUN_STATUS_NAME = "route_c_local_retrieval_status.json"

# Deliberately small and shared: words that carry no biology on their own must
# not drive a retrieval hit. Kept beside the blind-benchmark stopword list in
# spirit -- generic phenotype vocabulary is not a pathway.
_STOPWORDS = frozenset(
    "the a an and or of to in on with without at is are be for from by as that this these those "
    "increased decreased reduced increase decrease reduction higher lower more less altered abnormal "
    "plant plants phenotype trait gene genes protein proteins mutant mutants type wild".split()
)
_TOKEN = re.compile(r"[a-z][a-z0-9-]{2,}")

# Neural backends are audited, never required. A false here is a provider fact
# recorded in the manifest, not a failure.
_OPTIONAL_BACKENDS = ("sentence_transformers", "transformers", "faiss", "chromadb")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def provider_audit() -> dict[str, object]:
    """Which retrieval backends this workstation actually has."""
    return {
        "active_provider": "local_lexical_tfidf",
        "active_provider_requires_third_party": False,
        "external_transmission": "none: the phenotype text and the corpora never leave this machine",
        "optional_backends": {name: importlib.util.find_spec(name) is not None for name in _OPTIONAL_BACKENDS},
        "upgrade_path": (
            "an embedding backend is a provider swap: build a new index with its own database_version and "
            "model_version; existing frozen retrievals stay valid and comparable only within their own database"
        ),
    }


def _tokens(text: str) -> list[str]:
    return [token for token in _TOKEN.findall(str(text).casefold()) if token not in _STOPWORDS]


def _document_text(row: dict[str, object]) -> str:
    """One flat text field per reference record, with tabs/newlines removed so the
    frozen table stays a well-formed TSV."""
    parts = [str(row.get(column, "")).strip() for column in CORPUS_TEXT_COLUMNS]
    return re.sub(r"\s+", " ", " ".join(part for part in parts if part)).strip()


def _vector(tokens: list[str], idf: dict[str, float]) -> dict[str, float]:
    """L2-normalised TF-IDF vector. Unknown tokens are dropped, not smoothed in."""
    counts = Counter(token for token in tokens if token in idf)
    weights = {token: (1.0 + math.log(count)) * idf[token] for token, count in counts.items()}
    norm = math.sqrt(sum(value * value for value in weights.values()))
    return {token: value / norm for token, value in weights.items()} if norm else {}


# --------------------------------------------------------------------------- #
# index build
# --------------------------------------------------------------------------- #

def retrieval_index_build(
    corpus_args: list[str],
    out_dir_arg: str,
    database_version: str,
    corpus_label: str = "",
    full_hash: bool = False,
) -> dict[str, object]:
    """Freeze local source-backed reference records into a versioned index."""
    if not corpus_args:
        raise ValueError("retrieval-index-build requires at least one --corpus input")
    if not str(database_version).strip():
        raise ValueError("retrieval-index-build requires a non-empty --database-version")
    database_version = str(database_version).strip()
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)

    documents: list[dict[str, object]] = []
    excluded: list[dict[str, object]] = []
    corpus_inputs: dict[str, object] = {}
    seen_genes: set[str] = set()
    for corpus_arg in corpus_args:
        corpus_path = runtime_path(corpus_arg)
        frame = pd.read_csv(corpus_path, sep="\t", dtype=str).fillna("")
        missing = [column for column in CORPUS_REQUIRED_COLUMNS if column not in frame.columns]
        if missing:
            raise ValueError(f"reference corpus {corpus_path} is missing columns: {', '.join(missing)}")
        if not any(column in frame.columns for column in CORPUS_FAMILY_COLUMNS):
            raise ValueError(f"reference corpus {corpus_path} needs family_id and/or orthogroup_id for later interval matching")
        corpus_inputs[str(corpus_path)] = file_fingerprint(corpus_path, full_hash)
        for offset, row in enumerate(frame.to_dict("records"), start=2):  # header is line 1
            gene = str(row["reference_gene_id"]).strip()
            family = str(row.get("family_id", "")).strip()
            orthogroup = str(row.get("orthogroup_id", "")).strip()
            text = _document_text(row)
            reason = ""
            if not gene:
                reason = "empty_reference_gene_id"
            elif not (family or orthogroup):
                # Unmappable to an interval later: excluded before Top-K so it can
                # never consume a retrieval slot, but kept in the ledger.
                reason = "missing_family_and_orthogroup"
            elif not str(row.get("source_id", "")).strip():
                reason = "function_not_source_backed"
            elif not text:
                reason = "no_usable_description_text"
            elif not _tokens(text):
                reason = "text_carries_no_indexable_token"
            elif gene.casefold() in seen_genes:
                reason = "duplicate_reference_gene_id"
            if reason:
                excluded.append({"reference_gene_id": gene, "corpus_file": str(corpus_path),
                                 "corpus_row": offset, "exclusion_reason": reason})
                continue
            seen_genes.add(gene.casefold())
            documents.append({
                "doc_id": f"doc_{len(documents):06d}", "reference_gene_id": gene,
                "reference_species": str(row["reference_species"]).strip(),
                "family_id": family, "orthogroup_id": orthogroup,
                "document_text": text, "source_id": str(row["source_id"]).strip(),
                "source_url": str(row.get("source_url", "")).strip(),
                "corpus_file": str(corpus_path), "corpus_row": offset,
            })

    if not documents:
        raise ValueError("no indexable reference record survived the exclusion rules; see the exclusion ledger")

    tokenised = [_tokens(str(record["document_text"])) for record in documents]
    document_frequency = Counter(token for tokens in tokenised for token in set(tokens))
    total = len(documents)
    idf = {token: math.log((total + 1) / (count + 1)) + 1.0 for token, count in document_frequency.items()}
    postings: dict[str, dict[str, float]] = {}
    for index, tokens in enumerate(tokenised):
        for token, weight in _vector(tokens, idf).items():
            postings.setdefault(token, {})[str(index)] = round(weight, 8)

    index_dir = out / safe_id(database_version)
    index_dir.mkdir(parents=True, exist_ok=True)
    documents_path = index_dir / DOCUMENTS_NAME
    pd.DataFrame(documents, columns=DOCUMENT_COLUMNS).to_csv(documents_path, sep="\t", index=False, lineterminator="\n")
    postings_path = index_dir / POSTINGS_NAME
    postings_path.write_text(
        json.dumps({"idf": {token: round(value, 8) for token, value in sorted(idf.items())},
                    "postings": {token: postings[token] for token in sorted(postings)}},
                   ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    excluded_path = index_dir / EXCLUDED_NAME
    pd.DataFrame(excluded, columns=EXCLUSION_COLUMNS).to_csv(excluded_path, sep="\t", index=False, lineterminator="\n")

    manifest = {
        "schema_version": SCHEMA_VERSION, "index_format_version": INDEX_FORMAT_VERSION,
        "candidate_qtl_version": __version__,
        "database_version": database_version,
        "corpus_label": corpus_label or database_version,
        "built_at_utc": _now(),
        "command": "retrieval-index-build",
        "model_name": MODEL_NAME, "model_version": MODEL_VERSION,
        "provider": provider_audit(),
        "preprocessing": {
            "tokenizer": "casefold, [a-z][a-z0-9-]{2,}, shared stopword list, no stemming",
            "weighting": "tf=1+log(count), idf=log((N+1)/(df+1))+1, L2-normalised; score is cosine similarity in [0,1]",
            "text_fields": list(CORPUS_TEXT_COLUMNS),
            "determinism": "no randomness and no learned parameters; the same corpus rebuilds byte-identically",
        },
        "document_count": len(documents),
        "token_count": len(idf),
        "excluded_row_count": len(excluded),
        "exclusion_reasons": {reason: sum(1 for row in excluded if row["exclusion_reason"] == reason)
                              for reason in sorted({str(row["exclusion_reason"]) for row in excluded})},
        "species": sorted({str(record["reference_species"]) for record in documents if record["reference_species"]}),
        "source_backed": "every indexed row carries a source_id; rows without one are excluded, never indexed",
        "interval_information_in_index": False,
        "corpus_inputs": corpus_inputs,
        "artifact_sha256": {
            DOCUMENTS_NAME: sha256(documents_path),
            POSTINGS_NAME: sha256(postings_path),
            EXCLUDED_NAME: sha256(excluded_path),
        },
        "outputs": {"documents": str(documents_path), "postings": str(postings_path),
                    "excluded": str(excluded_path), "manifest": str(index_dir / INDEX_MANIFEST_NAME)},
        "safety": [
            "the index holds reference records only; no QTL interval, candidate id, or private sequence enters it",
            "rows without a family/orthogroup or without a source id are excluded before Top-K and kept in a ledger",
            "retrieval similarity computed from this index is a discovery route, never causal evidence",
        ],
    }
    manifest_path = index_dir / INDEX_MANIFEST_NAME
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), "documents": str(documents_path), "postings": str(postings_path),
            "excluded": str(excluded_path), "database_version": database_version,
            "document_count": len(documents), "excluded_row_count": len(excluded), "token_count": len(idf)}


# --------------------------------------------------------------------------- #
# phenotype-only retrieval
# --------------------------------------------------------------------------- #

def _load_index(index_arg: str) -> tuple[dict[str, object], pd.DataFrame, dict[str, object], Path]:
    manifest_path = runtime_path(index_arg)
    if manifest_path.is_dir():
        manifest_path = manifest_path / INDEX_MANIFEST_NAME
    if not manifest_path.exists():
        raise ValueError(f"index manifest not found: {manifest_path}; run retrieval-index-build first")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    index_dir = manifest_path.parent
    documents_path = index_dir / DOCUMENTS_NAME
    postings_path = index_dir / POSTINGS_NAME
    mismatched = [
        name for name, path in ((DOCUMENTS_NAME, documents_path), (POSTINGS_NAME, postings_path))
        if not path.exists() or sha256(path) != str(manifest.get("artifact_sha256", {}).get(name, ""))
    ]
    if mismatched:
        raise ValueError(
            f"frozen index artifact(s) do not match {manifest_path}: {', '.join(mismatched)}; "
            "rebuild the index rather than retrieving against a changed one"
        )
    if str(manifest.get("model_version")) != MODEL_VERSION:
        raise ValueError(
            f"index {manifest.get('database_version')} was built by model_version {manifest.get('model_version')} "
            f"but this build is {MODEL_VERSION}; rebuild the index (scores are not comparable across scorers)"
        )
    documents = pd.read_csv(documents_path, sep="\t", dtype=str).fillna("")
    postings = json.loads(postings_path.read_text(encoding="utf-8"))
    return manifest, documents, postings, manifest_path


def retrieval_run(
    index_arg: str,
    phenotype: str,
    out_dir_arg: str,
    top_k: int = 50,
    case_id: str = "",
    phenotype_context: str | None = None,
) -> dict[str, object]:
    """Retrieve reference records from the phenotype text alone.

    There is deliberately no interval / atlas / candidate argument on this
    command: Phase 1 must be unable to see the QTL interval. Interval matching
    happens later, in ``retrieval-interval-match``, against the frozen output.

    ``phenotype_context`` (optional ``phenotype_context.yaml``) is recorded in
    the status for audit only; the TF-IDF query stays the phenotype free text
    and no structured field is appended to it.
    """
    if not str(phenotype).strip():
        raise ValueError("retrieval-run requires a non-empty --phenotype")
    if top_k < 1:
        raise ValueError("top-k must be at least 1")
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    manifest, documents, postings, manifest_path = _load_index(index_arg)

    idf = {token: float(value) for token, value in postings["idf"].items()}
    query_tokens = _tokens(phenotype)
    query = _vector(query_tokens, idf)
    scores: dict[int, float] = {}
    for token, weight in query.items():
        for document_index, document_weight in postings["postings"].get(token, {}).items():
            scores[int(document_index)] = scores.get(int(document_index), 0.0) + weight * float(document_weight)

    records = documents.to_dict("records")
    # Deterministic order: score desc, then reference_gene_id. A zero-scoring
    # document is not retrieved at all -- an empty result is an honest state.
    ranked = sorted(((round(score, 6), str(records[index]["reference_gene_id"]), index)
                     for index, score in scores.items() if score > 0),
                    key=lambda item: (-item[0], item[1]))[:top_k]
    rows: list[dict[str, object]] = []
    for rank, (score, _gene, index) in enumerate(ranked, start=1):
        record = records[index]
        rows.append({
            "reference_gene_id": str(record["reference_gene_id"]),
            "reference_species": str(record["reference_species"]),
            "family_id": str(record["family_id"]), "orthogroup_id": str(record["orthogroup_id"]),
            "retrieval_score": score, "matched_modality": MODALITY,
            "matched_description_or_phenotype": str(record["document_text"])[:400],
            "model_name": MODEL_NAME, "model_version": MODEL_VERSION,
            "database_version": str(manifest["database_version"]), "rank": rank,
        })
    retrieval = pd.DataFrame(rows, columns=RETRIEVAL_COLUMNS)
    retrieval_path = out / RETRIEVAL_NAME
    retrieval.to_csv(retrieval_path, sep="\t", index=False, lineterminator="\n")

    matched_tokens = sorted(set(query) & set(postings["postings"]))
    status = {
        "schema_version": SCHEMA_VERSION, "route": "embedding_retrieval", "phase": "phase_1_phenotype_only",
        "command": "retrieval-run", "candidate_qtl_version": __version__,
        "case_id": case_id,
        "retrieved_at_utc": _now(),
        "lane": "phenotype_only",
        "interval_information_in_model_input": False,
        "interval_arguments": [],
        "separation_rule": (
            "this command accepts no interval, atlas, or candidate argument; interval matching runs only "
            "against the frozen output in a separate Phase 2 command (retrieval-interval-match)"
        ),
        "provider": provider_audit(),
        "model_name": MODEL_NAME, "model_version": MODEL_VERSION,
        "database_version": str(manifest["database_version"]),
        "database_built_at_utc": manifest.get("built_at_utc"),
        "index_manifest": str(manifest_path),
        "index_artifact_sha256": manifest.get("artifact_sha256"),
        "index_verified_before_retrieval": True,
        "corpus_sources": sorted(str(value) for value in manifest.get("corpus_inputs", {})),
        "preprocessing": manifest.get("preprocessing"),
        "indexed_document_count": int(len(documents)),
        "excluded_row_count": manifest.get("excluded_row_count"),
        "query_token_count": len(query_tokens),
        "query_tokens_in_index": matched_tokens,
        "top_k_requested": top_k,
        "retrieval_count": int(len(retrieval)),
        "retrieval_status": "retrieved" if len(retrieval) else "not_evaluated_no_indexed_token_matched_the_phenotype",
        "retrieval_sha256": sha256(retrieval_path),
        "rank_contribution": "none",
        "causal_claim": "similarity is an independent discovery route; it is never causal evidence and never enters a gene score",
        "structured_context": structured_context_block(
            phenotype_context, case_id=case_id, phenotype_text=phenotype,
            query_usage={
                "usage": "stored_only",
                "reason": ("the TF-IDF retrieval query is the phenotype free text exactly as before; "
                           "no structured field is appended to the query"),
            },
        ),
        "external_transmission": "none",
        "outputs": {"retrieval": str(retrieval_path), "status": str(out / RUN_STATUS_NAME)},
        "next_step": (
            "candidate-qtl retrieval-import --retrieval "
            f"{retrieval_path} --out-dir <phase1 freeze dir>; interval matching only afterwards"
        ),
    }
    (out / RUN_STATUS_NAME).write_text(json.dumps(status, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"retrieval": str(retrieval_path), "status": str(out / RUN_STATUS_NAME),
            "retrieval_count": int(len(retrieval)), "retrieval_status": status["retrieval_status"],
            "database_version": str(manifest["database_version"]), "model_version": MODEL_VERSION,
            "rank_contribution": "none"}
