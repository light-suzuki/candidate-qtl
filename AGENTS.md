# AGENTS.md

Repository: `light-suzuki/candidate-qtl`.

This file is the repository-level execution contract for AI coding/research agents and for humans doing the same work. It applies to the whole repository unless a more specific nested `AGENTS.md` applies (for example `src/candidate_qtl/AGENTS.md`). All agents use the same artifact and evidence contracts, not product-specific prompt conventions. The repository is Japanese-first; most detailed contracts are in Japanese.

Scientific policy has one owner: [canonical Wet membership](docs/CANONICAL-WET-MEMBERSHIP.ja.md) (`docs/CANONICAL-WET-MEMBERSHIP.ja.md`). Read it before making Tier/adoption decisions. Use [canonical Route B](docs/CANONICAL-ROUTE-B.ja.md) and the [ordinary-case runbook](docs/AGENTIC-QTL-CASE.ja.md) for execution. Do not create another policy document or candidate selector to explain the same seam.

## 1. Read order

Before non-trivial work, read:

1. `CORE-README.ja.md` — goal and architecture
2. `docs/CORE-CONTRACT.ja.md` — scientific invariants and change boundaries
3. the nearest nested `AGENTS.md`
4. only the implementation/runbook needed for the task

## 2. Canonical idea

The tool is a flow: from a QTL interval and a phenotype to an experimentally manageable Wet-validation candidate set of about 20 genes, retaining the true causal candidate(s). Support verification is a separate evidence attribute, not a prerequisite for Wet inclusion. The goal is not a causal Top-1 gene or an arbitrary 20-row file.

For an ordinary request such as "Are there genes in this QTL interval that could cause this phenotype?", a parent agent owns the end-to-end investigation. It carries the supplied FASTA/GFF3, QTL interval and phenotype through Route A and Direct C to an adopted Wet set and a readable report. **Route A and Direct C are required; interval-wide Route B is an explicit option** that the user or a specific question requests (for example, to examine every interval gene). The parent splits work and Wet boundary candidates into bounded tasks for workers, merges every assigned ID, and resumes frozen checkpoints. Do not ask the researcher to prepare per-gene TSVs, divide 1,000 genes manually, or stitch A/B/C results. Report missing execution as `not_evaluated`. Direct C may return no findings, but attempt and record it.

Without B, the Wet set holds only A/C supported and review candidates; it is never padded from interval annotation. When `wet-candidate-handoff` (or `practical_packet_run.py finish`) reports `route_a_rethink.state` = `required` (no Wet membership) or `recommended` (fewer than `--min-route-candidates`, default 3), rethink Route A once in a fresh generation: a new A worker with the same phenotype-side inputs, a new output directory, the earlier A freeze preserved, and no interval gene list or earlier candidate names. Then report any shortfall as is, or offer optional B. The threshold is a workload signal, not a scientific cutoff; it never pads or trims membership.

The ordinary wrapper `workflow/practical_packet_run.py` (`prepare` → `review` → `finish`) transports and checks saved responses; it does not launch researchers or provide search tools. Branch on the sealed stage policy: `ready_for_initial_candidates` finishes without `--reviews`; `awaiting_independent_reviews` (legacy content policy) keeps its original review requirements. A generated request or a complete response roster alone does not finish the case.

**Twenty is a guide, not a minimum count, truncation cap, or completion condition; never truncate retained candidates merely to reach 20.** Preserve existing eligible weaker review candidates, alternative mechanisms and ambiguous copies. A completed investigation may return 8 or 12 genes without count-driven exploration. Do not require every review candidate to have direct trait evidence, a verified claim or a long mechanistic explanation. Do not invent IDs, duplicate genes, upgrade evidence or claim causal-gene retention from the row count.

```text
Route A + Direct C, each independently frozen
  + optional Route B (interval-wide research)
  -> claim review -> per-paper checklist (nine answers) -> Tier T1-T6/U
  -> supported_member_union (verified subset only) + explicit review candidates
  -> family/copy context
  -> Wet: full gene ledger; T1 + existing eligible support/explicit review
  -> freeze wet_candidates.tsv
  -> only the first overflowing Tier -> boundary review; later Tiers stay ledger-only
```

## 3. Fresh context and lane isolation

Separate worker IDs do not establish route independence. Start scientific workers without the parent's conversation history, and record their context mode, model and allowed input artifacts in the run plan/receipt. If a harness forks full history by default, disable it explicitly. Route A receives phenotype-side inputs and schema only; Route B receives its own frozen material and the phenotype, not A outputs or parent notes; Direct C receives species and phenotype only. A parent that has read interval data may orchestrate but must not author the formal blind A response. Preserve any leaked generation as diagnostic and rerun the affected lanes without supplying its discoveries.

Check that each worker actually has the search/fetch tools it needs; a prompt saying "Web-enabled" does not supply tools. Record queries, retrieved sources and failures.

## 4. Repository hierarchy

| Layer | Responsibility |
|---|---|
| CORE | Route A, Direct C freeze/projection/members, freeze, supported union, family/copy, Wet handoff |
| OPTIONAL | interval-wide Route B research, Direct-C dedicated search, second-stage review, expression/synteny/SV/annotation evidence |
| SCREENING | `discover`, `auto-triage`, preliminary shortlists |
| DIAGNOSTIC | budget portfolios, ablations, legacy lexical/LSA/neural retrieval |
| HISTORICAL | frozen older generations; preserve provenance, not current authority |

Module existence is not canonical reachability. Benchmark datasets and evaluation tooling are not part of the public release.

## 5. Route A and Route B

Route A is **phenotype-first**: phenotype -> mechanism/family/reference hypothesis -> source/accession/sequence verification -> deterministic interval projection. Initial generation is interval-blind. Multiple copies may remain. Unverified claims never become trusted support.

Route B is **interval-first**: interval -> all genes -> family/orthogroup/domain/reference **compression** -> semantic interpretation -> source/transfer verification -> interval members. B0 is phenotype-free and cannot see A output before freeze. Family/subfamily/orthogroup is an **intermediate representation**, not the final answer. Process large representations in bounded shards that preserve all requested units, IDs and evidence; no per-shard Top-K quotas. Distinguish failed/unprocessed shards from reviewed units with no claims.

In ordinary Route B work the assigned researcher may use general Web and literature search throughout. Internal knowledge and search snippets are not verified evidence; freeze new material in a new bundle before source-bound claims. The parent owns any run-specific search restriction, round budget or early stop, and records it; a bounded partial result is valid but is not completed literature investigation.

New ordinary A/B runs use `initial_bibliography_v1`. Frozen legacy generations (for example `source_function_v1`) keep their original policy.

## 6. Direct Route C

C independently collects reported trait–gene/locus/variant facts without reading A/B outputs, QTL intervals, preselected genes/families, or the parent's research context. A duplicate of an A finding is legitimate; a source outside the eventual QTL is recorded, not pre-filtered. Do not convert generic function/homology guesses into direct empirical reports.

`independent C worker research (or optional route-c-direct-search) -> 17-column proposals -> route-c-direct-freeze -> route-c-direct-project -> route-c-direct-members -> Wet`

Support promotion needs a sealed source binding produced by an attested `curated_human` or `independent_agent_review` that passes the machine checks (verbatim quote of the frozen source, verified citation, frozen identity, retrieved scope) plus an authorized projection. Never fabricate either label. Cross-species, orthology/BLAST-only and conflicting projections stay review-only; outside/no-hit/unresolved stay in the ledger. `no_interval_match != outside_current_qtl`.

## 7. Claim review and Tier

The original researcher records paper ID/title, reported gene name/ID/species, all nine checklist answers and reasons. Machines transport these explicit values and audit database metadata separately; never fill or overwrite `reported_gene`/`reported_gene_id` from database lookups and call that a hallucination check. Unknowns keep checked reasons. A decision table derives T1–T6/U; the strongest applicable paper sets the gene Tier, with no paper-count bonus. Tier is relevance to the phenotype, separate from support, mapping quality and copy state. See the [LLM boundary](docs/LLM-BOUNDARY.ja.md).

## 8. Family is an intermediate representation

Use family/orthogroup/subfamily to compress, bridge species, ground identifiers and preserve paralog ambiguity. Do not optimize around `family -> one gene`. Copy discrimination needs copy-specific evidence; indistinguishable copies are kept together.

## 9. Wet review budget

Follow `docs/CANONICAL-WET-MEMBERSHIP.ja.md`, not duplicate numeric policies.

The ordinary handoff saves the full relevance inventory and admits T1 plus existing eligible supported/explicit-review genes and copies. It uses `upper_evidence_scoped_v3`: when eligible reported T1–T3 exist, that upper pool plus protected Direct C is the active scope; lower candidates are kept in the ledger as `not_selected_upper_tier_priority`, not as biological exclusions. Do not pad a short upper set.

Only the first Tier whose cumulative membership exceeds 25 enters boundary comparison. Higher Tiers stay frozen; later Tiers stay ledger-only without another comparison. Twenty-five is not a truncation cap, and within-Tier ranking or forced reduction is not required. `<=25` can complete without a boundary request, including 8 or 12 eligible genes. Boundary decisions use the original actor's reasoned assessments (`original_assessment_v1`) with hash/roster/identity checks. Deferral is not proof of non-causality; do not silently re-add explicitly deferred genes through generic fill. Use `--budget-policy strict` only for explicitly fixed capacity.

## 10. Scientific invariants

- not_evaluated / provider unavailable / failed search / not measured is not biological negative evidence;
- annotation absence is not sequence absence; paralog ambiguity is not counterevidence;
- generic homology, cross-species reports and route convergence do not by themselves establish focal trait support;
- LLM/model memory is not a citable source; unverified source/accession/function claims are not trusted evidence;
- a keyword, symbol or real quote does not prove gene/species/trait/experiment binding;
- LLM prose/confidence does not create numerical causal scores;
- known answers used for evaluation must not shape production aliases, prompts, weights or thresholds.

LLM output alone cannot silently add/remove genes from a frozen set or change support flags. Adoption changes go through the explicit, recorded, hash-validated consumer; the original freeze stays immutable.

## 11. Changes, validation and hygiene

- State one change class: `CORE_FIX`, `CORE_INTEGRATION`, `EXPERIMENTAL` or `GOVERNANCE` (documentation/contract sync without behavior change). Do not mix scientific selection changes with repository cleanup.
- Repair the existing A/B/C -> Wet chain before expansion. Do not add a Route, ranker, fixed family quota or provider merely to use the budget. Do not hide wiring defects with weights, aliases or fixture exceptions.
- Inspect the diff for unrelated changes and run targeted tests: `python scripts/run_test_partition.py --partition core` or `scripts/local_ci.sh test lint`. A unit/mock test is not a real-model end-to-end result. Report partial testing as partial.
- Never rewrite frozen runs to show improvement. Do not commit private sequences, VCFs, chromatograms, phenotype data or other private research data.
