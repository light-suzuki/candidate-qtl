"""Neutral import boundary for homolog screening primitives (#669)."""
from __future__ import annotations
from pathlib import Path
from typing import Any

def _run_interval_blastp(query_path: Path, subject_path: Path, output: Path, runner: str, num_threads: int) -> dict[str, object]:
    from .benchmarks import _run_interval_blastp as implementation
    return implementation(query_path, subject_path, output, runner, num_threads)

def homolog_blast_evidence(*args: Any, **kwargs: Any) -> dict[str, object]:
    from .benchmarks import homolog_blast_evidence as implementation
    return implementation(*args, **kwargs)

def evaluate_rank(ranked_path: Path, answer_gene_id: str) -> dict[str, object]:
    from .benchmarks import _evaluate_rank as implementation
    return implementation(ranked_path, answer_gene_id)
