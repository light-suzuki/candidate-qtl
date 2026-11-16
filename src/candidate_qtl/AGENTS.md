# AGENTS.md — production package scope

Issue #1002 makes first-stage T1 a focal-species same/very-close reported-phenotype candidate with a grounded interval homolog projection. Preserve every eligible copy and reported-gene/proxy/member roles. same_gene_as_focal may remain unknown/no; strict post-freeze identity is optional later work. Support and Direct C stay unchanged. The user requested the new rule as the standard without an additional Tier-policy version; preserve old artifacts for comparison.

Ordinary A/B research screens about 20 distinct paper candidates per original reference gene or shared mechanism. Record each paper's original checklist and reason, derive its initial Tier, then select five useful distinct papers in strongest-initial-Tier order for the next-stage evidence handoff. Five is the primary evidence set to carry forward; cache-first bibliography ID/title checks consume that same set. DOI/PMID aliases of one paper share one slot. Include originals and reviews; reuse shared research across copies. Retain every selected paper's facts and reasons so the next stage can consult multiple reports. Preserve all original observations, counterevidence, unselected papers and selection reasons; the handoff selection does not restrict the existing strongest-applicable-paper gene Tier rule. If a strong paper is inaccessible or inapplicable, try another access path or paper. Record screened IDs, queries, retrieval failures and any shortfall; one found paper alone is not investigation completion. These are research targets, never candidate admission quotas or evidence-count scores. Provider bibliographies remain retrieval context. The production source-verify default is ID/title metadata only; explicit content_v1 retains strict checks. New B Web sources must return structured material/bindings and pass route-b-semantic-bundle --extend-bundle into a new generation before claim merge. Never silently discard discoveries or admit unfrozen sources.

New A/B paper reports use paper_checklist_v1: the original research LLM records the paper title/ID, reported gene name/ID/species, yes/no/unknown checklist answers and one short reason during the existing literature investigation. A short decision table derives T1-T6/U; the strongest applicable paper sets the gene Tier, with no paper-count bonus. Database metadata checks are cache-first and audit-only: name/display differences and missing aliases never veto candidates. Match IDs and provider gene names/aliases; species labels are context, not an extra metadata match gate. Reported genes, proxy references and focal loci are separate roles. No extra paper-content reviewer, Tier-only search or strict name-string matching is required. Legacy observation/declared-Tier formats remain readable; missing facts are not negative evidence.

New ordinary A/B runs use `initial_bibliography_v1`: adopt actual first-stage candidates from the original biological judgment, cache-first paper ID/title checks and frozen identity/member correspondence. Assign initial T1–T6/U in the original proposal. Do not require independent content review or exact Tier confirmation in stage one; selective content checks belong to stage-two narrowing. Candidate adoption is separate from verified support; strict Direct C and identity/sequence/hash gates stay unchanged. Old `source_function_v1` bundles and unversioned generations replay their original policies. See the [canonical Wet policy](../../docs/CANONICAL-WET-MEMBERSHIP.ja.md).

Paper evidence tables are views of the same original paper_reports_json, never a second scientific reviewer. Record all nine checklist answers and checked unknown reasons; preserve source_excerpt/source_location when available and keep unfilled answers in paper_checklist_missing_fields_json. For initial checklist/bibliography recording, missing optional publication metadata, excerpts or locations alone do not request another research turn; keep their missingness visible. Explicit strict content and legacy recording requirements remain unchanged. mechanism-import and route-b-semantic-import export paper_evidence_table.tsv and one bounded paper_followup_requests.json batch for the original researcher. Complete required scientific blanks and concrete upper-Tier leads, not every hypothetical unknown; record followup_state=completed after the attempt. Use a new payload/output generation and preserve route independence/freezes. Requests are not model execution; recording_state=needs_completion is separate from candidate membership. Existing cached DOI/title and gene-name/ID audits feed the table without extra retrieval. Missing gene IDs with known names/species remain explicit unknowns, not an automatic Tier penalty.

Opt-in Light science-v2 preserves the original LLM-written paper ID/title and reported gene name/ID/species. Python derives Tier from the original checklist and existing applicability/projection gates; the researcher need not transcribe a Tier number. Shared observations use explicit references, never inferred biological equivalence. Keep copy-specific checklists/reasons and independent source-bound database function hypotheses. Member coverage distinguishes assessed, missing, invalid and failed retrieval even when a unit transports valid claims as completed. Legacy Light input and frozen generations remain readable.

This nested contract applies to `src/candidate_qtl/**` in addition to the repository-root `AGENTS.md`.

Before editing production code, read:

1. `CORE-README.ja.md`
2. `docs/CORE-CONTRACT.ja.md`
3. this file
4. the active GitHub Issue and dependency state
5. only the implementation/runbook docs required for that Issue

Issue #469 is a **closed historical stabilization roadmap**, not the current work controller. Use current `main`, current contracts, and active Issues.

## 1. Module classes

- **CORE** — Route A, Direct Route C freeze/projection/members, freeze/merge, `supported_member_union`, family/copy, Wet handoff
- **OPTIONAL** — interval-wide Route B research, Direct-C dedicated search, bounded refinement, extra evidence providers
- **SCREENING** — `discover`, `auto-triage`, preliminary shortlist producers
- **DIAGNOSTIC / EXPERIMENTAL** — budget portfolios, benchmarks, standalone retrieval/review seams
- **HISTORICAL** — superseded or frozen generations retained for provenance

A module existing under `src/candidate_qtl/` does not automatically make it production core.

## 2. Required change class

State exactly one:

- `CORE_FIX`
- `CORE_INTEGRATION`
- `BENCH_ONLY`
- `EXPERIMENTAL`
- `GOVERNANCE`

Documentation-only synchronization is `GOVERNANCE`.

## 3. Canonical production target

Base canonical correctness is Route A + Direct C (A/C-first, 2026-11-08). Route B is optional.

```text
Route A freeze + Direct-C verified-member input
        │
optional Route B freeze
        ↓
supported member union
        ↓
family/copy discrimination
        ↓
Wet budget gate
   ordinary target handoff -> full gene relevance ledger
           -> retain T1 + existing eligible supported/review candidates
           -> no count-driven T4/U investigation
           -> only first overflowing tier becomes boundary scope; later tiers stay ledger-only
           -> finish adopts/defers that sealed scope, without count-driven continuation
   legacy strict/supported_only/refinement -> tier1_verified_support + tier2_reference_phenotype -> tier3_family_or_domain; exclude tier4_dark
```

Review fill rows carry `membership_basis=interval_review_coverage`; they are not verified support and carry no route flag or score. `--wet-fill-mode supported_only`, strict capacity, and explicit refinement are diagnostic/compatibility behavior.

Initial 12 and final 20 are operational targets. Around 15-30 final candidates
requires no count-only warning, justification, or stop; this range is not a hard limit. Boundary review may complete with more than 25 candidates; the threshold selects the review scope and never truncates the adopted set.
Scientific verification gates remain mandatory for verified support. Grounded
A/B proposals with incomplete source/claim verification may enter Wet as explicit
review candidates ahead of generic interval fill, without entering the supported
union. Identity, sequence and interval-membership integrity remain mandatory.
Strict capacity is opt-in only.

Legacy B bundles use `b_review_membership_policy=source_function_v1`: checked reference-function content, interval correspondence and a specific conditional phenotype bridge are required for B-only review eligibility, without requiring focal causal papers or identical reported phenotypes. Source IDs or generic class/GO labels alone do not suffice. Preserve the ledger and independent A/C eligibility; old bundles replay their legacy policy. Exact shared function judgments retain separate claim/copy application states. The canonical Wet contract owns these rules.

Default relevance selection uses #817's `upper_evidence_scoped_v3`. Eligible reported
new T1-T3 (or legacy T1/T2) activates the upper pool plus protected C; legacy T3
is not a reported T3. Preserve original eligibility separately from active IDs.
Seal the policy version, original eligible IDs and active IDs with the inventory hash.
Lower candidates stay in the ledger as `not_selected_upper_tier_priority`, never
biological negatives. No upper candidates means the original eligible scope.
Twenty is not a minimum or completion condition. At active overflow above 25,
compare only the first overflowing tier using existing evidence; higher tiers stay
frozen, later tiers stay ledger-only. No within-tier ranking or forced reduction;
do not add a Tier/content reviewer. Only sealed v3 boundaries consume original actor
decisions as `original_assessment_v1` with hash/roster/identity checks and empty reviewer fields.
Missing/unknown decisions retain candidates and incomplete history. Old policies and
generic selection still require independent review. A/B candidates may be deferred by reasoned comparison;
Direct-C protection stays separate. Retain pending membership when review
cannot complete; request generation is not model execution. Old unversioned
scopes and `evidence_scoped_v2` keep their original policies. Completion/partial
review cannot restore lower candidates or start a lower round. `supported_only`,
strict capacity, and explicit refinement remain legacy diagnostic paths.

`wet_candidate_trace.tsv` records union/review presence, bound claim count,
relevance and first-stage disposition. Its
`upstream_llm_assessment=not_evaluated_by_handoff` value is an observability
contract, not a model-performance or causal-truth judgment.

## 4. Route A / Route B

### Route A — phenotype-first

`phenotype -> mechanism/family/reference -> verification -> deterministic interval projection`

- initial hypothesis is interval-blind;
- source/accession/species/sequence must be verified for supported membership;
- identity/species/sequence-confirmed proposals may project as review-only candidates while source/claim verification remains incomplete;
- multiple interval copies may remain;
- Route B and truth remain hidden before freeze.

### Route B — interval-first

`interval -> all genes -> family/orthogroup/domain compression -> semantic interpretation -> source/transfer verification`

- B0 is phenotype-free;
- family/orthogroup is an intermediate representation, not a final answer;
- raw hundreds of genes are not dumped directly into an LLM;
- Route A remains hidden before freeze;
- group support expands mechanically to all plausible copies.

Family/copy resolution uses the canonical semantic resolution and competitor-gate modules.

## 5. Direct Route C — current implementation boundary

Current main implements:

```text
independent harness-native research (optional route-c-direct-search)
-> existing 17-column proposal -> route-c-direct-freeze
-> route-c-direct-project
-> route-c-direct-members
-> optional supported-member-union / family-copy / Wet handoff
```

The guided path consumes a **pre-frozen** Direct-C artifact after B0 interval-atlas construction. Before freeze the independent C researcher receives only phenotype and focal species; no A/B findings, gene/family shortlist, QTL interval, truth, or parent research history. Use available harness-native Web tools rather than require a dedicated retrieval provider. Preserve duplicate reported genes/loci/variants and eventual no-interval-match records; proposal extraction and BLAST execution remain external to the profile.

The current member adapter can emit `verified_direct_support` only for a narrow S0/S1 D0/D1 same-cultivar/same-species exact/assembly-mapped subset. S0/S1 D2, cross-species, orthology/BLAST-only projections remain `reference_phenotype_review`; outside/no-hit/unresolved rows remain ledger-only.

The Route-C-specific structured source claim binding from #661 is implemented and the issue is closed. Support requires a sealed, authorizing binding sidecar. The producer accepts an attested `curated_human` or `independent_agent_review`; both pass the same machine checks of quote, citation, frozen identity and retrieved scope, and the method is recorded. Remaining open work includes:

- #663: order-invariant, provenance-complete projection reduction
- #662 (closed): retrieval failure isolation and saved-material reuse; strict search lineage is optional

A machine-checked complete binding plus projection authorization promotes Direct-C support without a separate human review. A conflicting, unauthorized or incomplete row falls back to review-only.

Additional rules:

- Direct C is part of base correctness; Route B is optional;
- freeze occurs before QTL projection;
- `no_interval_match != outside_current_qtl`;
- no-hit / provider failure / not evaluated is not negative evidence;
- generic homology or cross-species direct evidence does not automatically transfer focal-gene trait support;
- no separate C score or Top-20 producer;
- legacy lexical/LSA/neural retrieval is diagnostic/historical, not current Direct C.

## 6. Stabilization guardrails

Without an explicit active architecture Issue, do not introduce:

- a new Route;
- a new top-level CLI command family;
- a parallel Top-20 producer presented as canonical;
- a new adaptive stage;
- ranking weight/threshold tuning;
- benchmark-truth-derived aliases or exceptions;
- a hidden requirement that makes Route B mandatory for the A/C workflow;
- a policy that always collapses a family to one gene;
- an evidence meaning that bypasses existing verification contracts.

## 7. Scientific invariants

- `not_evaluated` is not negative evidence.
- Provider/search/measurement failure is not biological absence.
- Annotation absence is not sequence absence.
- Generic homology/conservation is context, not trait-specific support.
- Unverified/source-unknown text is not trusted support.
- `paralog_ambiguous` is uncertainty, not automatic counterevidence.
- LLM prose does not write numerical gene scores.
- Convergence is provenance, not a score bonus.
- Benchmark truth does not shape production behavior inside the same run.
- A symbol/cue occurrence cannot authorize a different unverified stable ID.
- A real quote cannot bypass gene/species/trait/experiment binding.
- Interval no-hit is not outside-QTL mapping.

## 8. LLM boundary

The original research LLM is responsible for interpreting other-plant reports and judging their correspondence to a sequence proxy. Ordinary reported-gene → proxy → focal-gene hypotheses may proceed without the reported gene's own sequence or a new proxy phenotype paper; no extra biological reviewer is required solely for this bridge. Preserve the original reported facts and uncertainty, and keep actual proxy identity/sequence and frozen member/hash checks. Full characterization of all interval genes is not a completion requirement. Use the [canonical proxy contract](../../docs/CANONICAL-WET-MEMBERSHIP.ja.md#報告遺伝子代理参照対象座位を分ける); the optional same-species identity assessment is a separate path.

Allowed:

- Route A proposes mechanisms/families/reference targets; machine verifies them.
- Route B interprets bounded family/orthogroup units; machine verifies source and transfer.
- over-budget review proposes which ambiguity/evidence to resolve.

Not allowed:

- prose/confidence directly scoring genes;
- unverified text selecting a sole preferred copy;
- LLM review adding/removing genes from a frozen set;
- provider failure becoming negative evidence;
- an agent-supplied Direct-C relation label creating support without structured binding.

## 9. Minimal-diff rule

For `CORE_FIX`:

1. identify one root cause;
2. touch only required production files;
3. add a regression test;
4. do not mix refactors, providers, metrics, or benchmark artifacts;
5. do not update performance headlines in the same PR;
6. verify the real canonical consumer path, not only a helper function.

For `GOVERNANCE`:

- do not change production behavior;
- synchronize every source-of-truth document that would otherwise contradict the change;
- add lightweight contract-consistency checks when a stale claim has already recurred.

## 10. Benchmark-discovered defects

Use:

`diagnostic frozen -> CORE_FIX Issue -> production PR -> realistic validation -> new benchmark generation`

Do not tune production from inside a truth-aware scoring loop.

## 11. Done means reachable and auditable

When claiming integration, verify:

- producer artifact reaches the intended consumer;
- the route is exposed by the claimed CLI/profile;
- final Wet handoff receives the canonical supported-member representation;
- missing providers remain explicit;
- source/projection/freeze hashes are checked;
- documented Direct-C support policy matches `route_c_direct_members.py` and `supported_member_union.py`;
- known blockers are not described as implemented guarantees.
