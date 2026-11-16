"""Local semantic and provider-fusion retrieval for Route C (issue #225).

Issue #217 delivered exactly one Route C provider: a deterministic local
lexical TF-IDF scorer (``retrieval_local.py``). It is a genuinely useful
reproducible baseline, but exact-token TF-IDF cannot see synonymy or
paraphrase -- "enlarged meristem" and "fasciated shoot apex" share no token,
even though they can describe the same biology.

This module adds three things, all inside Route C's existing honesty
contract (phenotype-only query -> frozen retrieval -> interval match only
after freeze; ``rank_contribution`` always ``none``):

``retrieval_semantic_index_build`` / ``retrieval_semantic_run``
    A second, versioned local provider: TF-IDF followed by truncated SVD
    (latent semantic analysis / LSI), cosine similarity in the reduced
    space. This is honestly named as **latent-semantic, not neural** -- it
    captures co-occurrence structure across the whole corpus, which lets it
    bridge some paraphrase gaps that raw term-matching cannot, but it is
    still a linear, corpus-scale model, not a learned language embedding.
    It needs no new dependency: ``numpy`` is already installed as a
    transitive dependency of ``pandas``. SVD sign is not unique across
    LAPACK implementations, so every retained component is canonicalised
    (its largest-magnitude term loading is forced positive) before the
    index is frozen, and the rank ``k`` is a pinned constant, never a CLI
    knob -- so the same corpus always rebuilds byte-identically.

``retrieval_neural_index_build`` / ``retrieval_neural_run``
    A third, optional-dependency provider id backed by ``sentence_transformers``,
    gated the same way ``retrieval_local.provider_audit()`` gates neural
    backends: audited, never assumed. This workstation's venv does not have
    ``sentence_transformers`` installed (confirmed by the same
    ``importlib.util.find_spec`` probe used elsewhere in this package), so
    the live encode/index path has never been exercised here and is not
    implemented in this build -- a stub that pretended to embed without ever
    having been run would be an untested guess dressed up as a result. The
    only honest behaviour when the backend is missing is to report
    ``not_evaluated_provider_unavailable`` and stop: never a crash, never a
    silent fallback to another provider, never a fabricated retrieval row.

``retrieval_hybrid_merge``
    A hybrid mode that fuses two or more provider retrievals by
    reciprocal-rank fusion (RRF) over each provider's *within-provider*
    rank. This is the provider-level counterpart of
    ``retrieval_route.retrieval_merge``'s within-database round-robin merge:
    that function never compares raw scores across *databases* (it says so
    in its own status JSON, ``cross_database_raw_scores_compared: false``);
    this one never compares raw scores across *providers*
    (``cross_provider_raw_scores_compared: false``). An LSA cosine can be
    negative and a lexical TF-IDF cosine cannot -- two incomparable scales --
    so the fused score is computed from rank alone.

Every provider here reuses ``RETRIEVAL_COLUMNS`` from ``retrieval_route.py``
by import (never redefines it), keeps the phenotype-only run signature with
no interval / atlas / candidate argument, and pins ``rank_contribution``
to ``"none"``.
"""

from __future__ import annotations

import importlib.util
import json
import math
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__
from .core import file_fingerprint, runtime_path, safe_id, sha256
from .family_labels import canonical_family_key
from .retrieval_local import (
    CORPUS_FAMILY_COLUMNS,
    CORPUS_REQUIRED_COLUMNS,
    CORPUS_TEXT_COLUMNS,
    DOCUMENT_COLUMNS,
    EXCLUSION_COLUMNS,
)
from .retrieval_route import FAMILY_COLUMNS, RETRIEVAL_COLUMNS, validate_retrieval_frame

# --------------------------------------------------------------------------- #
# provider ids shared with commands/retrieval.py's --provider flag
# --------------------------------------------------------------------------- #
PROVIDER_LEXICAL = "lexical"
PROVIDER_LSA = "lsa"
PROVIDER_NEURAL = "neural"
PROVIDERS = (PROVIDER_LEXICAL, PROVIDER_LSA, PROVIDER_NEURAL)

# --------------------------------------------------------------------------- #
# LSA (latent semantic analysis: TF-IDF -> truncated SVD -> cosine)
# --------------------------------------------------------------------------- #
SCHEMA_VERSION = "1.0"
INDEX_FORMAT_VERSION = "1.0.0"
# Independently versioned from retrieval_local.MODEL_VERSION on purpose: a
# tokenizer or weighting change in the lexical provider must not silently
# change this provider's behaviour, and vice versa -- each provider swap is
# its own model_version bump against its own index.
MODEL_NAME = "candidate-qtl-local-lsa"
MODEL_VERSION = "lsa-svd-1.0.0"
MODALITY = "text_latent_semantic"

# Pinned, not a CLI flag: making the SVD rank tunable would turn it into a
# knob that could be nudged to make a benchmark look better, which is exactly
# what this project's scoring/weighting rule forbids. A corpus with fewer
# available singular values than this cap simply uses all of them (recorded
# as lsa_rank_k_effective in the manifest).
LSA_RANK_K = 64

CORPUS_TEXT_COLUMN_NAMES = CORPUS_TEXT_COLUMNS

DOCUMENTS_NAME = "lsa_documents.tsv"
EXCLUDED_NAME = "lsa_excluded_rows.tsv"
LATENT_NAME = "lsa_latent_index.json"
INDEX_MANIFEST_NAME = "lsa_index_manifest.json"
RETRIEVAL_NAME = "route_c_lsa_retrieval.tsv"
RUN_STATUS_NAME = "route_c_lsa_retrieval_status.json"

# Kept identical in content to retrieval_local's stopword/token rules
# (same tokenizer philosophy) but owned locally so a bump to one scorer's
# tokenizer never silently changes the other's.
_STOPWORDS = frozenset(
    "the a an and or of to in on with without at is are be for from by as that this these those "
    "increased decreased reduced increase decrease reduction higher lower more less altered abnormal "
    "plant plants phenotype trait gene genes protein proteins mutant mutants type wild".split()
)
_TOKEN = re.compile(r"[a-z][a-z0-9-]{2,}")

_OPTIONAL_NEURAL_BACKENDS = ("sentence_transformers", "transformers", "faiss", "chromadb")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _tokens(text: str) -> list[str]:
    return [token for token in _TOKEN.findall(str(text).casefold()) if token not in _STOPWORDS]


def _document_text(row: dict[str, object]) -> str:
    parts = [str(row.get(column, "")).strip() for column in CORPUS_TEXT_COLUMN_NAMES]
    return re.sub(r"\s+", " ", " ".join(part for part in parts if part)).strip()


def _vector(tokens: list[str], idf: dict[str, float]) -> dict[str, float]:
    """L2-normalised TF-IDF vector, identical formula to retrieval_local's --
    the two providers apply the same term weighting before diverging: the
    lexical provider stops there, this one feeds the vector into SVD."""
    counts = Counter(token for token in tokens if token in idf)
    weights = {token: (1.0 + math.log(count)) * idf[token] for token, count in counts.items()}
    norm = math.sqrt(sum(value * value for value in weights.values()))
    return {token: value / norm for token, value in weights.items()} if norm else {}


def provider_audit() -> dict[str, object]:
    """What this workstation's LSA provider actually depends on, and what
    optional neural backends are also visible (upgrade-path transparency,
    same probe retrieval_local.provider_audit() uses)."""
    return {
        "active_provider": "local_lsa_svd",
        "active_provider_requires_third_party": True,
        "third_party_dependency": "numpy (already required transitively by pandas; no new install for this provider)",
        "external_transmission": "none: the phenotype text and the corpora never leave this machine",
        "optional_neural_backends": {name: importlib.util.find_spec(name) is not None for name in _OPTIONAL_NEURAL_BACKENDS},
        "upgrade_path": (
            "the neural adapter (--provider neural) is a separate provider id with its own model_name/model_version; "
            "it never silently replaces this LSA provider, and this LSA provider never silently claims to be neural"
        ),
    }


def retrieval_semantic_index_build(
    corpus_args: list[str],
    out_dir_arg: str,
    database_version: str,
    corpus_label: str = "",
    full_hash: bool = False,
) -> dict[str, object]:
    """Freeze local source-backed reference records into a versioned LSA index.

    Mirrors ``retrieval_local.retrieval_index_build``'s index-build/run split,
    immutable-artifact + SHA-256-manifest discipline, and corpus exclusion
    rules exactly; only the scorer differs. TF-IDF term-document matrix ->
    truncated SVD (``numpy.linalg.svd``) -> cosine similarity in the reduced
    latent space. This is latent-semantic analysis, not a neural embedding.
    """
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
    vocabulary = sorted(idf)
    vocab_index = {token: position for position, token in enumerate(vocabulary)}

    matrix = np.zeros((total, len(vocabulary)), dtype=float)
    for row_index, tokens in enumerate(tokenised):
        for token, weight in _vector(tokens, idf).items():
            matrix[row_index, vocab_index[token]] = weight

    n_docs, n_terms = matrix.shape
    rank_available = max(1, min(n_docs, n_terms))
    k_effective = min(LSA_RANK_K, rank_available)
    u, s, vt = np.linalg.svd(matrix, full_matrices=False)
    u_k = np.array(u[:, :k_effective])
    s_k = np.array(s[:k_effective])
    vt_k = np.array(vt[:k_effective, :])

    # Determinism guard: SVD sign is not unique across LAPACK
    # implementations/versions. Canonicalise every retained component so its
    # largest-magnitude term loading is positive, flipping the matching
    # document loadings in step so U @ diag(S) @ Vt stays consistent. Without
    # this, the *same corpus* could rebuild with flipped-sign components on a
    # different machine/numpy build even though the similarity geometry is
    # identical -- byte-identical rebuilds require removing that ambiguity.
    for component in range(k_effective):
        peak = int(np.argmax(np.abs(vt_k[component])))
        if vt_k[component, peak] < 0:
            vt_k[component, :] *= -1.0
            u_k[:, component] *= -1.0

    doc_latent = u_k * s_k[np.newaxis, :]
    doc_norms = np.linalg.norm(doc_latent, axis=1)
    doc_latent_unit = np.divide(
        doc_latent, doc_norms[:, np.newaxis],
        out=np.zeros_like(doc_latent), where=doc_norms[:, np.newaxis] > 0,
    )

    index_dir = out / safe_id(database_version)
    index_dir.mkdir(parents=True, exist_ok=True)
    documents_path = index_dir / DOCUMENTS_NAME
    pd.DataFrame(documents, columns=DOCUMENT_COLUMNS).to_csv(documents_path, sep="\t", index=False, lineterminator="\n")
    excluded_path = index_dir / EXCLUDED_NAME
    pd.DataFrame(excluded, columns=EXCLUSION_COLUMNS).to_csv(excluded_path, sep="\t", index=False, lineterminator="\n")

    latent_path = index_dir / LATENT_NAME
    latent_payload = {
        "idf": {token: round(float(value), 10) for token, value in sorted(idf.items())},
        "vocabulary": vocabulary,
        "singular_values": [round(float(value), 10) for value in s_k.tolist()],
        "term_loadings": [[round(float(value), 10) for value in row] for row in vt_k.tolist()],
        "document_latent_unit": [[round(float(value), 10) for value in row] for row in doc_latent_unit.tolist()],
    }
    latent_path.write_text(json.dumps(latent_payload, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")

    manifest = {
        "schema_version": SCHEMA_VERSION, "index_format_version": INDEX_FORMAT_VERSION,
        "candidate_qtl_version": __version__,
        "database_version": database_version,
        "corpus_label": corpus_label or database_version,
        "built_at_utc": _now(),
        "command": "retrieval-index-build",
        "provider_flag": PROVIDER_LSA,
        "model_name": MODEL_NAME, "model_version": MODEL_VERSION, "modality": MODALITY,
        "provider": provider_audit(),
        "preprocessing": {
            "tokenizer": "casefold, [a-z][a-z0-9-]{2,}, shared stopword list, no stemming",
            "weighting": "tf=1+log(count), idf=log((N+1)/(df+1))+1, L2-normalised TF-IDF row vector (same formula as the lexical provider)",
            "text_fields": list(CORPUS_TEXT_COLUMN_NAMES),
            "dimensionality_reduction": "truncated SVD (numpy.linalg.svd) of the dense TF-IDF term-document matrix -- latent semantic analysis / LSI, not a neural embedding",
            "similarity": "cosine similarity between the query and each document, both projected into the retained singular-vector subspace",
            "lsa_rank_k_configured": LSA_RANK_K,
            "lsa_rank_k_effective": k_effective,
            "rank_effective_reason": (
                "equal to the configured rank" if k_effective == LSA_RANK_K
                else f"clamped to min(n_docs={n_docs}, n_terms={n_terms}) because the corpus has fewer singular values than the configured rank"
            ),
            "sign_canonicalisation_rule": (
                "for every retained component, the term-loading entry with the largest absolute value is forced "
                "positive; the matching document-loading column is flipped in the same step, so the SVD sign "
                "ambiguity never changes the stored artifacts or a later retrieval"
            ),
            "determinism": (
                "no randomness; numpy.linalg.svd is a deterministic LAPACK call and sign canonicalisation removes "
                "the one remaining ambiguity, so the same corpus rebuilds byte-identically"
            ),
        },
        "document_count": n_docs,
        "term_count": n_terms,
        "excluded_row_count": len(excluded),
        "exclusion_reasons": {reason: sum(1 for row in excluded if row["exclusion_reason"] == reason)
                              for reason in sorted({str(row["exclusion_reason"]) for row in excluded})},
        "species": sorted({str(record["reference_species"]) for record in documents if record["reference_species"]}),
        "source_backed": "every indexed row carries a source_id; rows without one are excluded, never indexed",
        "interval_information_in_index": False,
        "corpus_inputs": corpus_inputs,
        "artifact_sha256": {
            DOCUMENTS_NAME: sha256(documents_path),
            LATENT_NAME: sha256(latent_path),
            EXCLUDED_NAME: sha256(excluded_path),
        },
        "outputs": {"documents": str(documents_path), "latent": str(latent_path),
                    "excluded": str(excluded_path), "manifest": str(index_dir / INDEX_MANIFEST_NAME)},
        "safety": [
            "the index holds reference records only; no QTL interval, candidate id, or private sequence enters it",
            "rows without a family/orthogroup or without a source id are excluded before Top-K and kept in a ledger",
            "retrieval similarity computed from this index is a discovery route, never causal evidence",
            "this is a latent-semantic (TF-IDF + truncated SVD) provider, not a neural embedding -- see --provider neural for the (currently unavailable) neural adapter",
        ],
    }
    manifest_path = index_dir / INDEX_MANIFEST_NAME
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), "documents": str(documents_path), "latent": str(latent_path),
            "excluded": str(excluded_path), "database_version": database_version,
            "document_count": n_docs, "excluded_row_count": len(excluded), "term_count": n_terms}


def _load_semantic_index(index_arg: str) -> tuple[dict[str, object], pd.DataFrame, dict[str, object], Path]:
    manifest_path = runtime_path(index_arg)
    if manifest_path.is_dir():
        manifest_path = manifest_path / INDEX_MANIFEST_NAME
    if not manifest_path.exists():
        raise ValueError(f"LSA index manifest not found: {manifest_path}; run retrieval-index-build --provider lsa first")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    index_dir = manifest_path.parent
    documents_path = index_dir / DOCUMENTS_NAME
    latent_path = index_dir / LATENT_NAME
    mismatched = [
        name for name, path in ((DOCUMENTS_NAME, documents_path), (LATENT_NAME, latent_path))
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
    latent = json.loads(latent_path.read_text(encoding="utf-8"))
    return manifest, documents, latent, manifest_path


def retrieval_semantic_run(
    index_arg: str,
    phenotype: str,
    out_dir_arg: str,
    top_k: int = 50,
    case_id: str = "",
) -> dict[str, object]:
    """Retrieve reference records from the phenotype text alone, via LSA cosine
    similarity in the frozen index's latent space.

    There is deliberately no interval / atlas / candidate argument on this
    command, matching every other Phase 1 retrieval command exactly.
    """
    if not str(phenotype).strip():
        raise ValueError("retrieval-run requires a non-empty --phenotype")
    if top_k < 1:
        raise ValueError("top-k must be at least 1")
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    manifest, documents, latent, manifest_path = _load_semantic_index(index_arg)

    idf = {token: float(value) for token, value in latent["idf"].items()}
    vocabulary = list(latent["vocabulary"])
    vocab_index = {token: position for position, token in enumerate(vocabulary)}
    term_loadings = np.array(latent["term_loadings"], dtype=float)
    doc_latent_unit = np.array(latent["document_latent_unit"], dtype=float)

    query_tokens = _tokens(phenotype)
    query_vector = _vector(query_tokens, idf)
    matched_tokens = sorted(token for token in query_vector if token in vocab_index)
    query_latent = np.zeros(term_loadings.shape[0] if term_loadings.size else 0, dtype=float)
    for token in matched_tokens:
        query_latent += query_vector[token] * term_loadings[:, vocab_index[token]]
    query_norm = float(np.linalg.norm(query_latent)) if query_latent.size else 0.0
    if query_norm > 0 and doc_latent_unit.size:
        scores = doc_latent_unit @ (query_latent / query_norm)
    else:
        scores = np.zeros(len(documents))

    records = documents.to_dict("records")
    ranked = sorted(
        ((round(float(scores[index]), 6), str(records[index]["reference_gene_id"]), index)
         for index in range(len(records)) if scores[index] > 0),
        key=lambda item: (-item[0], item[1]),
    )[:top_k]
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

    status = {
        "schema_version": SCHEMA_VERSION, "route": "embedding_retrieval", "phase": "phase_1_phenotype_only",
        "command": "retrieval-run", "provider_flag": PROVIDER_LSA, "candidate_qtl_version": __version__,
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
        "model_name": MODEL_NAME, "model_version": MODEL_VERSION, "modality": MODALITY,
        "provider_type_declaration": "latent-semantic (TF-IDF + truncated SVD), not a neural embedding",
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


# --------------------------------------------------------------------------- #
# neural adapter (sentence_transformers), correctly gated (issue #225 #2)
# --------------------------------------------------------------------------- #
NEURAL_BACKEND_MODULE = "sentence_transformers"
NEURAL_SCHEMA_VERSION = "1.0"
MODEL_NAME_NEURAL = "candidate-qtl-local-neural-sentence-transformers"
MODEL_VERSION_NEURAL = "st-adapter-1.0.0"
MODALITY_NEURAL = "text_dense_neural_embedding"
NEURAL_INDEX_MANIFEST_NAME = "neural_index_manifest.json"
NEURAL_RETRIEVAL_NAME = "route_c_neural_retrieval.tsv"
NEURAL_RUN_STATUS_NAME = "route_c_neural_retrieval_status.json"

# Matches the "not_evaluated_<reason>" vocabulary already used throughout the
# package (e.g. retrieval_local's not_evaluated_no_indexed_token_matched_the_phenotype,
# and the phrase named in retrieval_local.py's own module docstring for the
# pre-#217 state of Route C).
NOT_EVALUATED_PROVIDER_UNAVAILABLE = "not_evaluated_provider_unavailable"
NOT_EVALUATED_ADAPTER_NOT_IMPLEMENTED = "not_evaluated_adapter_not_implemented"


def neural_provider_available() -> bool:
    """Whether sentence_transformers is importable here. A fact, not a wish:
    this workstation's venv has only numpy among the optional retrieval
    backends (biopython/pandas/pysam/PyYAML are the base dependencies)."""
    return importlib.util.find_spec(NEURAL_BACKEND_MODULE) is not None


def neural_provider_audit() -> dict[str, object]:
    available = neural_provider_available()
    return {
        "active_provider": "local_neural_sentence_transformers",
        "active_provider_requires_third_party": True,
        "backend_module": NEURAL_BACKEND_MODULE,
        "backend_importable": available,
        "external_transmission": (
            "none: local sentence-transformers inference only, no network call at query time"
            if available else "not_applicable: provider unavailable, no inference attempted"
        ),
        "reason_if_unavailable": "" if available else f"{NEURAL_BACKEND_MODULE} is not importable in this Python environment",
    }


def retrieval_neural_index_build(
    corpus_args: list[str],
    out_dir_arg: str,
    database_version: str,
    corpus_label: str = "",
    full_hash: bool = False,
) -> dict[str, object]:
    """Adapter entry point for a sentence-transformers-backed index.

    ``sentence_transformers`` is not importable on this workstation, so the
    live encode-and-index path has never been run here and is deliberately
    NOT implemented in this build: writing an unexercised embedding pipeline
    would be an untested guess dressed up as a result, which Route C's
    honesty contract forbids (``retrieval_local.py``: "failure is not a
    miss"). When the backend is missing this always reports
    ``not_evaluated_provider_unavailable`` and writes no index artifact --
    never a crash, never a silent fallback to another provider.
    """
    del corpus_args  # declared input surface only; nothing is read while unavailable/unimplemented
    if not str(database_version).strip():
        raise ValueError("retrieval-index-build requires a non-empty --database-version")
    database_version = str(database_version).strip()
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    audit = neural_provider_audit()
    index_dir = out / safe_id(database_version)
    index_dir.mkdir(parents=True, exist_ok=True)

    if audit["backend_importable"]:
        index_build_status = NOT_EVALUATED_ADAPTER_NOT_IMPLEMENTED
        reason = (
            f"{NEURAL_BACKEND_MODULE} is importable here, but the live encode/index path is not implemented in "
            "this build; add it together with an integration test run in an environment that actually has "
            f"{NEURAL_BACKEND_MODULE} installed, rather than shipping an unverified path"
        )
    else:
        index_build_status = NOT_EVALUATED_PROVIDER_UNAVAILABLE
        reason = audit["reason_if_unavailable"]

    manifest = {
        "schema_version": NEURAL_SCHEMA_VERSION,
        "candidate_qtl_version": __version__,
        "database_version": database_version,
        "corpus_label": corpus_label or database_version,
        "built_at_utc": _now(),
        "command": "retrieval-index-build",
        "provider_flag": PROVIDER_NEURAL,
        "model_name": MODEL_NAME_NEURAL, "model_version": MODEL_VERSION_NEURAL, "modality": MODALITY_NEURAL,
        "provider": audit,
        "index_build_status": index_build_status,
        "reason": reason,
        "document_count": 0,
        "external_transmission": audit["external_transmission"],
        "safety": [
            "no embedding is computed and no real index artifact is written while the backend is unavailable or unimplemented",
            "a missing/unimplemented backend is reported as an explicit status, never as a fabricated retrieval result",
        ],
        "outputs": {"manifest": str(index_dir / NEURAL_INDEX_MANIFEST_NAME)},
    }
    manifest_path = index_dir / NEURAL_INDEX_MANIFEST_NAME
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"manifest": str(manifest_path), "database_version": database_version,
            "document_count": 0, "index_build_status": index_build_status}


def retrieval_neural_run(
    index_arg: str,
    phenotype: str,
    out_dir_arg: str,
    top_k: int = 50,
    case_id: str = "",
) -> dict[str, object]:
    """Phenotype-only retrieval against the neural adapter index.

    No interval / atlas / candidate argument, matching every other Phase 1
    retrieval command. On this workstation this always resolves to
    ``not_evaluated_provider_unavailable``: it never fabricates a retrieval
    row and never falls back to another provider.
    """
    if not str(phenotype).strip():
        raise ValueError("retrieval-run requires a non-empty --phenotype")
    if top_k < 1:
        raise ValueError("top-k must be at least 1")
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)
    manifest_path = runtime_path(index_arg)
    if manifest_path.is_dir():
        manifest_path = manifest_path / NEURAL_INDEX_MANIFEST_NAME
    if not manifest_path.exists():
        raise ValueError(f"neural index manifest not found: {manifest_path}; run retrieval-index-build --provider neural first")
    index_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    audit = neural_provider_audit()

    retrieval = pd.DataFrame(columns=RETRIEVAL_COLUMNS)
    retrieval_path = out / NEURAL_RETRIEVAL_NAME
    retrieval.to_csv(retrieval_path, sep="\t", index=False, lineterminator="\n")

    retrieval_status = str(index_manifest.get("index_build_status", NOT_EVALUATED_PROVIDER_UNAVAILABLE))
    reason = index_manifest.get("reason", audit["reason_if_unavailable"])
    status = {
        "schema_version": NEURAL_SCHEMA_VERSION, "route": "embedding_retrieval", "phase": "phase_1_phenotype_only",
        "command": "retrieval-run", "provider_flag": PROVIDER_NEURAL, "candidate_qtl_version": __version__,
        "case_id": case_id, "retrieved_at_utc": _now(), "lane": "phenotype_only",
        "interval_information_in_model_input": False, "interval_arguments": [],
        "separation_rule": (
            "this command accepts no interval, atlas, or candidate argument; interval matching runs only "
            "against the frozen output in a separate Phase 2 command (retrieval-interval-match)"
        ),
        "provider": audit,
        "model_name": MODEL_NAME_NEURAL, "model_version": MODEL_VERSION_NEURAL, "modality": MODALITY_NEURAL,
        "database_version": str(index_manifest.get("database_version", "")),
        "index_manifest": str(manifest_path),
        "index_build_status": index_manifest.get("index_build_status"),
        "top_k_requested": top_k,
        "retrieval_count": 0,
        "retrieval_status": retrieval_status,
        "reason": reason,
        "retrieval_sha256": sha256(retrieval_path),
        "rank_contribution": "none",
        "causal_claim": "similarity is an independent discovery route; it is never causal evidence and never enters a gene score",
        "external_transmission": audit["external_transmission"],
        "fabrication_guard": "no retrieval row is ever emitted when the backend is unavailable or unimplemented; the empty contract-shaped table above is the honest output",
        "outputs": {"retrieval": str(retrieval_path), "status": str(out / NEURAL_RUN_STATUS_NAME)},
    }
    (out / NEURAL_RUN_STATUS_NAME).write_text(json.dumps(status, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"retrieval": str(retrieval_path), "status": str(out / NEURAL_RUN_STATUS_NAME),
            "retrieval_count": 0, "retrieval_status": retrieval_status,
            "database_version": status["database_version"], "model_version": MODEL_VERSION_NEURAL,
            "rank_contribution": "none"}


# --------------------------------------------------------------------------- #
# hybrid fusion: reciprocal-rank fusion across providers (issue #225 #4)
# --------------------------------------------------------------------------- #
# Cormack, Clarke & Buettcher 2009's constant. Pinned as code, not exposed as
# a CLI flag -- for the same reason LSA_RANK_K is pinned: it is a fusion
# mechanism parameter, not a score weight to tune towards a benchmark result.
RRF_K = 60
HYBRID_MODEL_NAME = "candidate-qtl-hybrid-rrf"
HYBRID_MODEL_VERSION = "rrf-k60-1.0.0"
HYBRID_MODALITY = "hybrid_reciprocal_rank_fusion"
HYBRID_SCHEMA_VERSION = "1.0"
HYBRID_MERGED_NAME = "route_c_retrieval_hybrid.tsv"
HYBRID_EXCLUDED_NAME = "route_c_retrieval_hybrid_excluded.tsv"
HYBRID_STATUS_NAME = "route_c_retrieval_hybrid_status.json"


def retrieval_hybrid_merge(
    retrieval_args: list[str],
    out_dir_arg: str,
    top_k: int = 50,
    case_id: str = "",
) -> dict[str, object]:
    """Fuse two or more provider retrievals by reciprocal-rank fusion (RRF)
    over each provider's within-provider rank.

    This is the provider-level counterpart of
    ``retrieval_route.retrieval_merge``'s within-database round-robin: that
    function never compares raw scores across *databases*
    (``cross_database_raw_scores_compared: false``); this one never compares
    raw scores across *providers* (``cross_provider_raw_scores_compared:
    false``). A row's fused ``retrieval_score`` is
    ``sum(1 / (RRF_K + within_provider_rank))`` over the providers that
    retrieved it -- the source providers' own ``retrieval_score`` column is
    never read by this function, so an LSA cosine (which can be negative)
    and a lexical TF-IDF cosine (which cannot) are never put on the same
    axis. Rows lacking family/orthogroup annotation are excluded before
    Top-K, matching ``retrieval_merge``.
    """
    if len(retrieval_args) < 2:
        raise ValueError("retrieval-hybrid-merge requires at least two --retrieval inputs")
    if top_k < 1:
        raise ValueError("top-k must be at least 1")
    out = runtime_path(out_dir_arg)
    out.mkdir(parents=True, exist_ok=True)

    providers: list[tuple[str, Path, pd.DataFrame]] = []
    seen_providers: set[str] = set()
    excluded_rows: list[dict[str, object]] = []
    provider_counts: dict[str, dict[str, object]] = {}
    for value in retrieval_args:
        path = runtime_path(value)
        frame = validate_retrieval_frame(
            pd.read_csv(path, sep="\t", dtype=str).fillna(""),
            f"retrieval input {path}", require_family=False,
        ).copy()
        provider_keys = sorted({f"{row['model_name']}={row['model_version']}" for row in frame.to_dict("records")})
        if len(provider_keys) != 1:
            raise ValueError(f"retrieval input {path} must contain exactly one model_name/model_version pair (one provider)")
        provider = provider_keys[0]
        if provider in seen_providers:
            raise ValueError(f"retrieval-hybrid-merge provider must be unique per input: {provider}")
        seen_providers.add(provider)
        frame["rank"] = pd.to_numeric(frame["rank"], errors="coerce").astype(int)
        annotated_mask = frame[list(FAMILY_COLUMNS)].astype(str).apply(
            lambda series: series.str.strip().ne("").any(), axis=1,
        )
        for row in frame[~annotated_mask].to_dict("records"):
            excluded_rows.append({
                **{column: row.get(column, "") for column in RETRIEVAL_COLUMNS},
                "source_provider": provider, "source_input": str(path),
                "exclusion_reason": "missing_family_or_orthogroup",
            })
        annotated_frame = frame[annotated_mask]
        providers.append((provider, path, annotated_frame))
        provider_counts[provider] = {
            "input_rows": int(len(frame)), "annotated_rows": int(annotated_mask.sum()),
            "excluded_unannotated_rows": int((~annotated_mask).sum()),
        }

    # RRF fusion: only `rank` is ever read from the source rows below. Each
    # unique (gene, family, orthogroup) identity accumulates
    # 1 / (RRF_K + rank) from every provider that retrieved it.
    fused: dict[tuple[str, str, str], dict[str, object]] = {}
    for provider, _path, frame in providers:
        for row in frame.to_dict("records"):
            key = (
                str(row.get("reference_gene_id", "")).casefold(),
                canonical_family_key(row.get("family_id", "")),
                canonical_family_key(row.get("orthogroup_id", "")),
            )
            contribution = 1.0 / (RRF_K + int(row["rank"]))
            entry = fused.setdefault(key, {"row": row, "rrf_score": 0.0, "database_versions": set(), "_best_rank": None})
            entry["rrf_score"] = float(entry["rrf_score"]) + contribution
            entry["database_versions"].add(str(row.get("database_version", "")).strip())
            # The representative row is whichever provider retrieved this
            # identity at the best (numerically lowest) within-provider rank
            # -- a deterministic, rank-only tie-break, never a score comparison.
            if entry["_best_rank"] is None or int(row["rank"]) < entry["_best_rank"]:
                entry["row"] = row
                entry["_best_rank"] = int(row["rank"])

    ranked_keys = sorted(
        fused.items(),
        key=lambda item: (-round(float(item[1]["rrf_score"]), 10), str(item[1]["row"].get("reference_gene_id", ""))),
    )[:top_k]

    rows: list[dict[str, object]] = []
    for rank, (_key, entry) in enumerate(ranked_keys, start=1):
        source_row = entry["row"]
        database_versions = sorted(value for value in entry["database_versions"] if value)
        rows.append({
            "reference_gene_id": str(source_row.get("reference_gene_id", "")),
            "reference_species": str(source_row.get("reference_species", "")),
            "family_id": str(source_row.get("family_id", "")),
            "orthogroup_id": str(source_row.get("orthogroup_id", "")),
            "retrieval_score": round(float(entry["rrf_score"]), 10),
            "matched_modality": HYBRID_MODALITY,
            "matched_description_or_phenotype": str(source_row.get("matched_description_or_phenotype", ""))[:400],
            "model_name": HYBRID_MODEL_NAME, "model_version": HYBRID_MODEL_VERSION,
            "database_version": "+".join(database_versions) if database_versions else "unknown",
            "rank": rank,
        })
    merged = pd.DataFrame(rows, columns=RETRIEVAL_COLUMNS)
    if merged.empty:
        raise ValueError("retrieval-hybrid-merge produced no annotated rows")
    validate_retrieval_frame(merged, "hybrid-fused retrieval")
    merged_path = out / HYBRID_MERGED_NAME
    merged.to_csv(merged_path, sep="\t", index=False, lineterminator="\n")
    excluded_path = out / HYBRID_EXCLUDED_NAME
    pd.DataFrame(
        excluded_rows, columns=RETRIEVAL_COLUMNS + ["source_provider", "source_input", "exclusion_reason"],
    ).to_csv(excluded_path, sep="\t", index=False, lineterminator="\n")

    status_path = out / HYBRID_STATUS_NAME
    status = {
        "schema_version": HYBRID_SCHEMA_VERSION, "route": "embedding_retrieval_hybrid", "case_id": case_id,
        "lane": "phenotype_only", "interval_information_in_model_input": False, "interval_arguments": [],
        "fusion_policy": "reciprocal_rank_fusion_by_within_provider_rank",
        "rrf_constant_k": RRF_K,
        "cross_provider_raw_scores_compared": False,
        "cross_provider_raw_scores_summed": False,
        "fusion_rule": (
            "fused_score = sum over contributing providers of 1 / (RRF_K + within_provider_rank); "
            "the source providers' own retrieval_score column is never read by this function"
        ),
        "providers": [provider for provider, _path, _frame in providers],
        "top_k": top_k,
        "selected_rows": int(len(merged)),
        "excluded_rows": int(len(excluded_rows)),
        "per_provider": provider_counts,
        "rank_contribution": "none",
        "inputs": {str(path): file_fingerprint(path) for _provider, path, _frame in providers},
        "output_fingerprints": {
            "merged": file_fingerprint(merged_path, full_hash=True),
            "excluded": file_fingerprint(excluded_path, full_hash=True),
        },
        "outputs": {"merged": str(merged_path), "excluded": str(excluded_path), "status": str(status_path)},
    }
    status_path.write_text(json.dumps(status, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"merged": str(merged_path), "excluded": str(excluded_path), "status": str(status_path),
            "selected_rows": int(len(merged)), "excluded_rows": int(len(excluded_rows))}
