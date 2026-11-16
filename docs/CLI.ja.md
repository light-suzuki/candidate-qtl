# candidate-qtl CLI

candidate-qtl はQTLの統計検出器ではない。既に得たQTLについて、根拠を保全した候補遺伝子トリアージ（区間内candidate抽出、既存OSSを使うDry解析の準備、親間direct sequencing primer、確認済み親間変異からの選抜marker、F2共分離集計）と複数経路の仮説探索・次実験の意思決定支援を接続するワークフローCLIである（主対象: Windows 11 + WSL2。Linuxでも実行可能）。

## 座標 provenance preflight

候補抽出より前に、source 座標・marker・文献探索状態から最も強い経路を非対話で決定する。

~~~bash
candidate-qtl qtl-source-init \
  --manifest templates/qtl_source_manifest.yaml \
  --out-dir out/case/preflight
~~~

`qtl_source_preflight.json` は `qtlift_interval`、`marker_liftover`、
`marker_recovery`、`unresolved_position` のいずれかを返す。位置未解決時は
`positional_ranking_allowed=false` と `candidate_scope=rescue` を固定し、必要な介入は
`review_required.tsv` のみに集約する。

## 非対話 case runner

~~~bash
candidate-qtl case-run --case templates/case_manifest.yaml \
  --resources templates/local_resources.yaml --mode deterministic \
  --out-dir out/example --resume
~~~

同一manifest/resource hashの完了済みrunは再計算しない。位置未解決caseはRescue-only reportになる。
marker入口は `liftover-markers`、split区間は複数の `extract --interval`、primer用配列は
`extract --flank N` を使う。`direction-evidence` は逆方向をcounterevidenceへ変換し、
`candidate-flags` はtandem paralogと短鎖no-homology ORFを明示する。

`case-run` はcoordinate、extract、candidate flags、evidence normalization、deterministic rank、
primer import/plan、reportをstep別hashで再開する。source IDのないevidenceは
`evidence_excluded.tsv`へ隔離され、Topへ影響しない。Topはgene ID順ではなく、強い遺伝証拠、
または独立した複数evidence familyを必要とする。

`--mode llm-assisted`はhuman-reviewedかつsource-backedなcanonical evidenceだけを
`llm_assisted/`へ取り込み、deterministic rankを上書きしない。`run_manifest.json`には
resource fingerprint、利用tool、step hash、LLM呼出有無を保存する。

~~~bash
candidate-qtl liftover-markers --marker-fasta markers.fa --db target_genome_db \
  --out-dir out/marker_liftover --max-product 5000 --max-block-gap 1000000
~~~

marker IDは `_F` / `_R` で終える。pairとして向き・距離が成立する産物だけを採用し、
marker順序が非単調、複数contig、または離れたblockの場合は単一min-max spanを許可しない。

## 最小入力からの候補発見

~~~powershell
candidate-qtl discover `
  --fasta genome.fa --gff3 genes.gff3 `
  --interval chr5:71000000-112000000 `
  --phenotype "fasciated stem" --direction increase --inheritance recessive `
  --out-dir results
~~~

注釈providerは、明示した`--protein-db`、設定済みローカルDB、NCBI BLAST APIの順に選ぶ。
ローカルDBがなくても既定では外部APIへfallbackする。外部への配列送信を禁止する場合だけ
`--offline`を付ける。API候補上限を超えた場合は一部候補だけを恣意的に送らず、
`review_required.tsv`へローカルDB利用または区間縮小を要求する。

主要出力は`candidate_shortlist.tsv`、全候補順位、canonical evidence、source trace、
Dark/Rescue分類、Markdown reportである。候補の最終判断は人間またはLLMが行い、
このコマンド自体は原因遺伝子を断定しない。

## Windowsからの実行

PowerShellではWindowsパスをそのまま渡せる。

~~~powershell
Set-Location /path/to/candidate-qtl
.\scripts\candidate-qtl.ps1 doctor --json
.\scripts\candidate-qtl.ps1 extract --fasta C:\private\assembly.fa --gff3 C:\private\assembly.gff3 --interval chr1:100000-500000 --out-dir C:\private\qtl_01\extract --json
~~~

内部ではWSLのPythonで実行し、Windowsパスを `/mnt/...` に変換する。WSLから直接実行する場合は以下を使う。

~~~bash
PYTHONPATH=src python3 -m candidate_qtl doctor --json
PYTHONPATH=src python3 -m candidate_qtl --help
~~~

### ネイティブ Windows で動く範囲

pysam(Linux/macOS専用)は使用箇所で import するため、package の import と**表計算だけのcommandはネイティブ Windows の Python でも動く**。順位の確認、Wet状態の更新、次実験の比較といった日常操作でWSLを往復しなくてよい。

~~~powershell
py -m pip install pandas PyYAML biopython
py -m pip install --no-deps -e .
py -m candidate_qtl doctor --json    # pysam_available / pysam_dependent_commands を確認
py -m candidate_qtl evidence-rank --candidates cand.tsv --evidence evidence.tsv --out-dir rank
~~~

配列 I/O を行うcommand(`extract`、`discover`、`case-run`、`auto-triage`、`annotation-rescue`、
`exon-sequencing-plan`、`gff-preflight`、`guided-profile`、`benchmark-reference`、
`candidate-pipeline`)は pysam が必要で、ネイティブ Windows では実行方法を示して終了コード1で止まる。

~~~text
error: extract needs pysam, which has no native Windows build.
  - In WSL or Linux: python -m pip install -c constraints.txt -e '.[alignment]'
  - From Windows PowerShell: use scripts\candidate-qtl.ps1, which forwards to WSL.
~~~

正確な一覧は `doctor` の `pysam_dependent_commands` が返す。ゲノムを扱う本番runは従来どおり
WSL(`scripts\candidate-qtl.ps1`)経由で実行する。

## コマンド

コマンドの完全な一覧と引数は生成物 [`CLI-REFERENCE.md`](CLI-REFERENCE.md) が正本です（`python scripts/generate_cli_reference.py` で再生成。CLIを意図的に変えたら再生成してコミットしてください）。command数も生成物から取得し、手書き文書へ固定値を転記しません。以下は全コマンドの一覧と役割です。役割文は生成物の各コマンド節のhelp冒頭文（英語）をそのまま記載しています。

| コマンド | 役割 |
|---|---|
| exon-coverage-check | Detect CDS windows PrimerBlast OSS results do not span and plan relaxed retries. |
| benchmark-fixture-validate | Validate a committed benchmark fixture against its own manifests and report (T8 recurrence prevention, issues #348/#351/#352/#353/#354): freeze manifest files[] inventory (size+SHA-256), catalog ledger rows vs catalog_total, rank_metrics row-position/tie-aware/Stage-A columns, firewall case-level review_required propagation into the waterfall, started<=frozen<=completed timestamp order, and report headline counts vs the score tables. Exit code 1 on errors, 0 on warnings only; --json prints the full check manifest. |
| command-log-validate | Issue #381: audit command_log.tsv files mechanically. Legacy mode (default) classifies the 12 measured T11 header variants, reports argv completeness (full/medium/short/none), exit-code / start+end UTC / wall-time coverage, timestamp notation (iso_z / run_id_style / datetime_no_zone), missingness (empty or '-' cells) per canonical field, stage coverage with measured name normalization, duplicate stage rows, end-before-start rows and naive timestamps; a directory scans cases/*/command_log.tsv recursively, and aggregated multi-header files are split into per-header blocks. Strict mode requires the canonical 13-column schema (docs/COMMAND-LOG.ja.md) and fails on missing required fields, invalid argv_json/status/exit_code, naive or unparseable timestamps and completed<started. Expected-stage assertions run only with --manifest ({"schema_version": "1.0", "expected_stages": [...]}); without it the report says not_evaluated_expected_stage_set_missing. Exit code 1 on structural/parse errors, 0 on findings only; never rewrites the log. |
| benchmark-report-compare | Issue #367: paired two-run benchmark comparison. Compares two per-case evaluation tables (supervised_score_cases.tsv-style case_state or a generic evaluation_status column) under the #367 reporting contract: the roster/truth/metric-contract fingerprints must match (any missing or mismatched fingerprint is not_comparable, fail-closed), then every case's transition and the paired hit-rate delta with its seed-fixed bootstrap 95% CI are written out. Never changes scores, ranks, candidate sets, or frozen fixtures; rank_contribution=none. |
| benchmark-reference | Run known-gene virtual-QTL recovery cases on one reference without using answer IDs for scoring. |
| primerblast-import | Normalize PrimerBlast OSS JSON reports into an orderable marker ledger. |
| assay-plan | Legacy top-candidate primerblast-oss assay plan. |
| primer-plan | Generate primerblast-oss plans for direct sequencing, parental-variant markers, or a coarse QTL marker grid. |
| exon-sequencing-plan | Tile every Top-candidate representative CDS for parental PCR/Sanger. |
| collect-sanger-markers | Collect PrimerBlast OSS tile results into an ordering table. |
| cross-crop-benchmark | Run a reproducible cross-crop, coverage-aware meta-ranking benchmark across resource conditions. |
| robustness-benchmark | Drive the rank gate under annotation-degradation and hard-negative variants; verify no plausible-but-wrong lone Top and a safe hold. |
| discovery-benchmark | Measure honest recall with the REAL answer left out (cross-species pathway homology kept). Re-ranks with real-answer-derived evidence removed. No weight is changed. |
| blind-benchmark | Phenotype-first blind (autonomous-floor) benchmark: derive pathway/family from the phenotype alone, no human-chosen homolog. Records results as-is, including zero recall. Distinct from the map-based benchmark. |
| budget-run | Deterministic budget engine: tier0 inventory, tier1 clusters, tier2 budgeted sequence-search plan, per-gene review handoff, the resource-usage ledger under a declared budget profile, and the deterministic / LLM-assisted Top-N portfolio. Never ranks or scores. |
| budget-llm-review | Batched, budget-gated LLM cluster review + critic + Top-N portfolio reducer over an existing budget-run output directory (issue #375 T5). Reviews only Tier1/2-selected clusters in small batches, rejects unverified mechanism/source ids, and never converts model prose into a number. |
| build-executor-package | Build a physically separated Executor package (code + executor manifests + declared reference metadata only) with an allowlist package manifest and a run_capsule.json. Fails closed when any truth/scoring artifact or link would enter the package. |
| audit-predict | Answerless prediction phase: derive the blind route from the phenotype alone, emit one deterministic ranking per case plus a route-status table, and freeze them with SHA-256/UTC/input fingerprints. Rejects any answer/truth column; truth-dependent lanes are recorded not_evaluated. |
| audit-score | Post-freeze scoring phase: re-hash every frozen ranking, refuse to score a mutated freeze, then score it against the truth. Writes only into its own out-dir, so scoring can never change a frozen prediction. |
| knowledge-coverage | Record per-candidate interpretation coverage (annotation / family / homolog / reference function / literature / ontology / expression / variant) in three states: available, knowledge_unavailable (the resource was consulted and holds nothing) and not_evaluated (no provider was supplied). Neither absence is negative evidence and none of it changes a rank. |
| mechanism-coverage-benchmark | Score Route A phenotype-to-mechanism coverage separately from candidate recall, always next to what the coverage cost: distinct mechanisms after collapsing rephrasings, duplicate count, mechanisms-per-hit, and failure attribution (mechanism_not_proposed / mechanism_proposed_but_unverified / verified_mechanism_no_interval_support). Reads frozen artifacts and never re-runs a route. |
| benchmark-stratify | Recompute per-case benchmark metrics inside each stratum (knowledge-coverage tier, lineage specificity, interval-size band ...). Strata are reported side by side and never averaged: recall over well-studied cases and over knowledge-limited cases are two measurements. |
| masking-benchmark | Issue #237: score frozen Route A/B/C artifacts produced under corpus-mask's tiers (and the unmasked baseline) as a novelty-stress lane kept separate from ordinary benchmarking. Reports exact_family_recovery and neighbouring_process_recovery as separate columns, candidate_inflation versus the same case's unmasked candidate/hypothesis set, and per-metric not_evaluated denominators. Never re-runs a route. |
| benchmark-knowledge-mask | Type 1–3 benchmark construction only: Type 2 (`leave_truth_semantic_edges_out`) removes only truth-linked function/phenotype/pathway/source edges from Route B reference-function, while neutral assembly mapping and orthology/family projection rows remain intact; Type 1/3 tiers additionally remove nested exact-gene, family/orthogroup, and declared sequence-homolog knowledge edges. Candidate-universe/presence rows are preserved; raw truth IDs and mask paths/content are not emitted. A Type 2 no-op is `not_evaluated`. |
| reference-panel-benchmark | Issue #238: fingerprint reference-species panels and/or score frozen Route A/B artifacts produced under them. Reports coverage by species/panel, A/B family and exact-gene recall where evaluable, distinct-primary-source-id counts (a mirrored experiment is one source, not several), and stratifies failure into missing_ortholog / paralog_ambiguity / absent_phenotype_function / annotation_gap via knowledge-coverage's benchmark-stratify. Declares no universal best panel; exposes per-case winners instead. |
| next-experiment-benchmark | Score rule-based next-experiment escalation on frozen per-case artifacts: exact-action and appropriate-escalation accuracy, blocked-technology detection kept separate from inappropriate choices, and missing-provider cases reported not_evaluated. |
| query-eligibility | Gate a homolog-query manifest to non-Pisum functional homologs only; bar any Pisum-derived query before the benchmark runs. |
| discovery-chain-benchmark | A/B discovery chain end-to-end benchmark: aggregate externally-run lane artifacts (A-only / B0-only / union / convergence / full chain) into lane metrics and first-failure taxonomy. The harness never runs the lanes. |
| supervised-benchmark-validate | Validate the supervised-benchmark-v1 Executor-visible case manifest (docs/BENCHMARK-MANIFEST-SCHEMA.ja.md #1) and unconditionally reject any truth-shaped column (answer_*/truth_*/causal_*/known_gene*/*solution*/*label*, prefixed or suffixed, case-insensitive) before any case is processed. Row-level issues (missing field, unknown benchmark_profile/allowed_external_search) are reported per case and do not by themselves fail the command. |
| supervised-benchmark-truth-validate | Validate the Supervisor-only truth package (docs/BENCHMARK-MANIFEST-SCHEMA.ja.md #2) against its own schema, separate from the Executor-visible manifest: truth_unit_type (gene/regulatory_region/sv/pav/multi_gene/locus), truth_provenance_tier (tier1/tier2/tier3/unresolved), profile_disposition (primary/secondary/ceiling_only/stress_only/case_study_only/excluded). Never accepted by an Executor command. |
| supervised-benchmark-predict-prepare | Bridge a validated Executor case manifest (supervised-benchmark-validate schema) into the exact audit-predict manifest (case_id/phenotype_text/interval/interval_genes): derives the interval from interval_id when it parses as CHROM:START-END (or reads an explicit `interval` column), resolves candidate_input relative to the manifest, rejects truth-shaped columns unconditionally, verifies the 1-to-1 case_id join, fingerprints every input, and emits the exact post-freeze audit/score commands. |
| supervised-benchmark-source-log-freeze | Freeze the Executor source log (and optional case roster) with SHA-256 before any truth file is opened, and emit the exact post-freeze supervised-benchmark-audit / supervised-benchmark-score commands. supervised-benchmark-score re-verifies this freeze with --source-log-freeze before opening the truth package, and refuses to score a source log that changed after the freeze. |
| supervised-benchmark-audit | Grade a per-case Executor source/query log against the Supervisor-only truth package's identity and known-answer-derived-source lists, in audit-after-run Mode B (no live interception required). Classifies by provenance, not proximity: truth_family_id is never consulted, so a phenotype-derived non-focal homolog is GREEN/allow even when it shares the truth family; only an exact match against the truth's own identity or a known answer-derived source turns a row RED, and RED becomes contaminated only if the row records it was actually used. Rolls up to a per-case clean/review_required/contaminated state. |
| supervised-benchmark-score | Classify every case frozen by audit-predict as hit/miss/not_evaluated/contaminated. Refuses to score if the frozen prediction's SHA-256 no longer matches the freeze manifest (reusing audit_freeze's own verification helper). Gene-truth clean cases are scored by calling audit-score unchanged (no parallel ranking system); when --discovery-chain-cases is supplied, clean cases without a scoreable single-gene truth instead consume discovery-chain-benchmark's own union_hit/full_chain_hit columns. Contaminated and not_evaluated cases are excluded from both the hit count and the recall denominator, and are never converted into a miss. |
| supervised-benchmark-route-family-metrics | Issue #350: aggregate the route family recovery lanes (Route A / B0 / B1 / union / A u B0 / C retrieval / gene rank) with explicit per-metric denominators. Every denominator is built from per-case requirements, never hard-coded: family recall = family_scoreable (truth_family_id present) + route evaluable + clean source audit; A u B0 additionally requires the union to have stayed on B0; C = family_scoreable + C evaluable (route_c_benchmark_cases.tsv) + clean; gene rank = gene truth + in universe + rank_metric_eligible + clean. Per (case, metric) the table carries eligible/hit/miss/not_evaluated/not_scoreable/review_required/contaminated; a not_evaluated case is never converted into a miss. Consumes frozen artifacts only. |
| reference-panel | Validate and availability-audit a comparative reference-resource panel (native required, others optional with machine-readable skip reasons). |
| external-rank | Import external comparative ranks into a coverage-aware ledger (RRF consensus) without touching the native dry_score. |
| propagate-terms | Deterministically propagate reference-ortholog function terms onto QTL interval candidates (no LLM); native-only input degrades gracefully. |
| comparative-locus-reconstruct | Reconstruct an optional candidate universe from flanking markers and precomputed multi-species locus caches; never adds causal score. |
| comparative-locus-benchmark | Evaluate frozen comparative reconstructions on synthetic/degraded-reference cases. |
| doctor | Inspect required and optional WSL tools. |
| provenance-ledger | Build and audit the normalized provenance DAG: nodes, edges, and shared derivation paths. |
| discover | Run preliminary candidate screening from genome, interval, and phenotype inputs; promote reviewed outputs to case-run for the full workflow. |
| discover-promote | Create a reviewable case-run manifest from preliminary discover outputs. |
| homolog-panel-auto | Generate a sequence-backed, pending-review homolog panel from phenotype text without inspecting the QTL interval. |
| homolog-panel-match | Match an automated homolog panel only against interval proteins and optionally freeze the deterministic Stage-A 20. |
| homolog-panel-benchmark | Blindly compare automated phenotype-to-homolog retrieval on frozen known-gene cases. |
| auto-triage | Process every interval candidate against role-labelled reference proteomes without a curated reference-gene panel. |
| candidate-elimination | Manage candidate and mechanism survival state. Only hard positional/genetic evidence may exclude a gene; mechanism-specific negatives reject one mechanism; assay failure, partial coverage and unmeasured data never eliminate anything. |
| direction-evidence | Mark phenotype-direction matches and counterevidence deterministically. |
| evidence-build | Normalize annotation, expression, position, and manual evidence into one source-backed ledger. |
| evidence-term-promote | Promote approved annotation terms into low-confidence deterministic screening evidence. |
| peak-context | Import QTL peak, LOD support interval, and recombinant core into auditable position evidence. |
| trait-terms | Build a controlled, approval-gated trait-term set for keyword evidence. |
| trait-ontology-map | Map free trait terms to ontology IDs using a local synonym table only. |
| qtl-source-init | Validate coordinate provenance and choose a safe liftover route. |
| extract | Extract genes and sequences from one QTL interval. |
| gff-preflight | Cross-check candidate-qtl GFF interpretation against optional AGAT/gffread outputs. |
| annotation-rescue | Project known proteins only into a QTL interval and audit missing/broken native models without changing rank. |
| liftover-markers | Resolve pair-level marker hits into safe target interval blocks. |
| dark-candidates | Export short complete no-homology ORFs as an explicit dark bucket. |
| candidate-flags | Flag tandem redundancy and dark short-ORF candidates. |
| annotation-keywords | Link best protein hits to GFF descriptions and emit auditable trait-keyword evidence. |
| annotate-plan | Create candidate-only protein, GO evidence template, and local homology/domain commands. |
| research-packet | Create AI-optional candidate research packets and prompts. |
| llm-prompt-export | Export provider-neutral, reproducible candidate research prompts as JSONL. |
| homolog-linking-export | Export provider-neutral prompt package for agent-in-the-loop homolog linking. |
| homolog-linking-import | Mechanically validate and import agent linking proposals with strict closed-world and grounding substring gates. |
| homolog-linking-counterfactual | Evaluate what truth Stage-A rank would be if linking proposals acted as support, without altering actual rank. |
| homolog-linking-probe-export | Export control memorization probe prompt with phenotype ONLY (no interval, no candidates). |
| homolog-linking-score | Score linking results stratified by memorization probe into recall_positive / recall_negative. |
| llm-shortlist-prompt | Export one optional Stage-B prompt that ranks the exact frozen Stage-A set. |
| llm-shortlist-import | Validate an LLM batch ranking without changing Stage A or scoring evidence. |
| llm-review-benchmark | Compare deterministic and optional LLM review ranks on identical frozen shortlists. |
| llm-run | Run prompt JSONL through OpenAI Responses or a local OpenAI-compatible API. |
| llm-evidence-import | Import only human-reviewed, source-backed model claims. Each claim needs reviewed=true plus reviewed_by; unresolved source IDs are marked and cannot count as independent. |
| llm-pathway-prompt | Export a versioned phenotype-only prompt package asking for competing pathway/gene-family hypotheses and non-focal reference baits. The model never sees or ranks QTL interval candidates. |
| llm-pathway-import | Normalize raw model output into the strict pathway-proposal schema, keeping the raw payload verbatim. Every proposal stays pending_verification. |
| hypothesis-verify | Verify proposed accessions against an authoritative database: existence, species match, retrievable sequence, source-backed function. Unverified proposals never become baits. |
| bait-panel-freeze | Freeze verified proposals into an immutable bait panel, FASTA, and a query manifest consumable by query-eligibility. |
| blast-panel-freeze | Freeze a panel TSV into an immutable bait panel plus a manifest. Only the panel is read; interval, candidate, truth and score inputs are structurally absent. |
| blast-panel-verify | Re-hash a blast-panel freeze manifest against its frozen artifacts and refuse any drift. |
| blast-panel-search | BLAST the verified frozen baits against the interval proteins and emit auditable blast evidence. Missing tools or a provider failure degrade to not_evaluated_*, never a miss. |
| phenotype-pathway-benchmark | Score the phenotype-only proposal lane stage by stage (pathway/family recall, accession existence, verification pass, downstream causal recovery). Never pooled with the map-based lane. |
| scaffold-rescue | Conservatively rescue unplaced / unanchored scaffolds that are plausibly linked to the QTL locus, as a separate candidate class. Being unplaced is not itself a reason; family agreement alone is never high confidence. |
| cosegregation | Summarize phenotype distribution by a candidate marker genotype. |
| qtl2-plan | Validate F2 tables and generate or run an R/qtl2 interval scan. |
| guided-profile | Start from a screening, full-case, answerless audit, or triple-route profile without hiding the stable commands or review gates. |
| case-run | Run/resume a deterministic manifest-driven case and emit a provenance report. |
| challenge-benchmark-aggregate | Issue #366: aggregate the pre-registered discriminative challenge fixture against a predictions directory (docs/BENCHMARK-CHALLENGE.ja.md §8): verifies the fixture freeze, reads run_manifest / truth / difficulty_tags, consumes per-case rankings and route family lists only, and emits seven deterministic artifacts. Never re-runs a route; B0 unavailable is not_evaluated, never a miss; sv/multi_gene truth units stay out of the gene denominator. |
| paired-benchmark | Paired benchmark on one case roster: (a) deterministic ABC vs (b) ABC + LLM multi-bait. Reads frozen per-arm artifacts only (route A/B/C families, ranked, narrowing, family-to-copy, multi-bait) and reports route coverage / family recall / preferred-copy evaluability / Stage-A Top-20 / narrowed Top-20 / reduction ratio / wrong-confident-promotion separately per arm. An arm without its own artifact is not_evaluated for that metric, never a miss; an unconfigured or failed LLM provider makes the multi-bait arm not_evaluated while the deterministic arm keeps running. Weights, thresholds and deterministic scores are never changed. |
| orchestrate | Generate a laptop-safe OSS workflow for candidate annotation and primer/marker design. |
| family-index | Build the cached genome-wide family index for one reference genome from a BLAST/DIAMOND hit table. Run once per reference; per-case falsification is then table lookup only. |
| interval-atlas | Structure every QTL-interval gene by family / domain / orthogroup / broad function without using any phenotype. An exploration index for forming hypotheses from what the interval actually contains; it never scores or promotes a gene. |
| family-to-copy | Compare the members of a recovered gene family against each other and record what does and does not discriminate between the copies. A first-Wet priority layer, never a causal claim. |
| family-candidate-bridge | Normalize Route A / B / C family evidence (phenotype-first hypotheses, orthology/atlas, frozen retrieval match) into one reporting table with a strict mapping tier (tier1_source_backed_unique .. tier6_not_evaluated). Reporting-only: rank_contribution is always none, no score columns, and the family_to_copy / narrowing wiring is untouched. At least one route input is required. |
| assembly-bridge | Normalize native same-species cross-assembly evidence into a full-ID-preserving mapping plus the existing Route B orthology/reference-function TSVs. Same fam3 alone stays provisional; only unique reciprocal sequence with family concordance can verify 1:1, and missing mappings stay not_evaluated. |
| multi-reference-locus-union | Build one candidate universe from a primary interval and repeatable alternate-reference intervals. Only verified 1:1 equivalents collapse to one locus entity; one-to-many copies, alternate-only genes and PAV candidates remain explicit. Emits a separate Route-B compression input so universe recall is not confused with conditional compression. |
| phenotype-first-status | Route A (forward lookup): normalize the phenotype-first artifacts and detect the failure signals that should trigger the interval-first route. |
| interval-reverse-lookup | Route B (reverse lookup): start from the interval's genes, go out to their model-species homologs, and derive family/process hypotheses from what is actually there. Interval-informed, so never scored as blind discovery. |
| integrated-narrowing | Combine base biological ranking (ranked_candidates.tsv) with family recovery and family-to-copy state into an experimental review order: ordinal narrowing_tier / narrowing_priority / first-Wet priority. A narrowing / Wet-routing layer, never a causal claim: dry_score / dry_rank / dry_group / top_eligibility are passed through untouched and rank_contribution stays none. |
| route-b-prepare | Build the local Route B provider stack (homology search -> orthology table, cached source-backed reference functions, interval atlas) and run the phenotype-free B0 reverse lookup. Reports readiness: a broad-bin-only run is not_evaluated_for_reverse_lookup, never an equivalent Route B result. |
| reverse-literature | Route B stage 2: for every interval gene's reference homolog, query the literature focused by the phenotype. A returned article is a link candidate, not evidence; an offline query is not_evaluated, never 'no literature exists'. |
| dual-route-compare | Compare the independently generated forward and reverse routes and record the convergence tier: strict mechanistic / family / gene overlap / process-related / text-only, plus phenotype-only / interval-only / discordant / both-weak. Convergence adds no score. |
| pathway-interval-falsify | Falsify pathway/gene-family hypotheses against the QTL interval, coarse bucket first and the fine ortholog level only for survivors. Enrichment filters hypotheses; it never scores a gene. |
| route-a-interval-match | Route A stage: forward-match only the verified/source-backed baits frozen by bait-panel-freeze to concrete genes inside the QTL interval via a reference_genome-scoped family index. Wires pathway-interval-falsify with verified targets only; zero verified baits, a missing index or an unresolved interval degrade to not_evaluated, never to a miss. |
| homology-concordance | Compare BLAST/DIAMOND/MMseqs2, domain, and orthology providers without vote inflation. |
| reference-candidate-table | Create a review table linking QTL candidates to source-backed reference genes through orthology mappings. |
| orthology-plan | Generate reciprocal-BLAST refinement only for Top/Rescue candidates. |
| orthology-import | Convert reciprocal-BLAST refinement results into canonical evidence rows. |
| orthology-tier1-plan | Profile-aware Tier-1 plan: candidate query vs focal self search and a reference panel. |
| focal-paralog-import | Flag out-of-interval paralog redundancy from a focal self-search. |
| orthofinder-import | Import an OrthoFinder Orthogroups.tsv and classify 1:1/1:many/many:many. |
| omamer-import | Import an OMAmer HOG placement table into orthogroup evidence. |
| synteny-import | Import MCScanX/jcvi synteny anchors into synteny evidence. |
| rank | Rank a candidate ledger with deterministic, inspectable rules. |
| dry-prioritize | Rank parental-sequencing candidates using pre-Wet evidence only. |
| evidence-rank | Rank candidates from the canonical evidence ledger with capped, inspectable contributions. |
| ranking-qa | Audit ranking stability across eligibility, family dropout, weight perturbation, confidence, provider, source-group dropout, and provider-dropout paths. |
| ranking-stability-benchmark | Report whether frozen-case truth genes stay within Top-K under weight perturbation and source/provider dropout. |
| screening-shortlist | Freeze the non-LLM Stage-A screening shortlist before human or LLM review. |
| evidence-tiers | Issue #229: report target_size-independent biological evidence tiers (strong/plausible/unresolved/contradictory/knowledge_limited) side by side with the configurable review_budget/assay_capacity operational cut. The operational cut never implies a biological exclusion. |
| evidence-tier-benchmark | Issue #229: Recall@budget alongside evidence-tier recall, stratified by interval-size band via knowledge_coverage.benchmark_stratify. Reads frozen evidence-rank outputs only; never re-ranks. |
| weight-policy-benchmark | Compare numeric weight baselines and structural ranking gates on a frozen known-gene benchmark. |
| mapping-gate | Assess whether QTL mapping is adequate before candidate triage. |
| next-experiment | Rank next experimental actions separately from biological candidate rank. |
| candidate-pipeline | Run consistent primary evidence/ranking and optionally plan or integrate Top-candidate orthology refinement. |
| phenotype-investigate | Route A stage 1: build the interval-blind phenotype research packet (direction/inheritance/tissue/stage/onset/penetrance/synonyms/ontology/uncertainty) plus a prompt package and a blank payload template an LLM or a human fills identically. |
| phenotype-localization-packet | Route A coverage: phenotype-only localization promptを作る。 |
| phenotype-localization-freeze | Route A coverage: localization responseを厳密検証してhash-freezeする。 |
| mechanism-import | Route A stage 2: import an LLM or human mechanism payload through one strict schema. Fewer than the required competing mechanisms, or a hypothesis without a counter-hypothesis, discriminating observation, source id or model-species evidence, is review-required. |
| mechanism-coverage-audit | Route A coverage: 初期12件を監査し、必要時だけ最大3件のblind補充promptを出す。 |
| reference-gene-resolve | Route A stage 2.5: provider由来の参照accession候補を、親仮説との対応・個別の出典とともに保持する。名前検索はunreviewedも含み、件数で絞らない。詳細は 参照候補の保持。 |
| reference-gene-cache | 凍結したRoute A仮説の参照名・種から、明示的なAPI許可でUniProt生レコードを取得しSHA付きcacheを作る。canonical prepareはmanifestを照合してofflineで解決する。 |
| source-verify | Verify literature IDs/titles with Europe PMC and an official NCBI PubMed/PMC E-utilities fallback. The default initial_bibliography_v1 checks metadata only; explicit content_v1 also classifies what the source supports. An unresolvable PMID/DOI is rejected only after a completed lookup; citation_verified and claim_support_status remain separate, and provider failure is not_evaluated rather than negative evidence. |
| claim-verification-benchmark | Score the claim-support classifier against labelled examples (mechanistic support, associative support, irrelevant mention, fabricated citation): per-class accuracy and a confusion table, with cases whose record could not be retrieved reported not_evaluated instead of wrong. |
| investigation-attribution | Attribute a Route A outcome to the stage where it stopped: phenotype research, verification, bait freeze, interval match, or candidate recovery. Truth-free, so it runs inside a prediction-only audit. |
| literature | Run one read-only Europe PMC literature query and save normalized JSON. |
| research | Research selected candidates through read-only APIs without using an AI. |
| reference-corpus-fetch | Download source-backed reference records (UniProtKB: function, disruption phenotype, pathway, tissue, GO, evidence codes and the PubMed ids each block cites) for one or more taxa. Requires --allow-external-api; there is no free-text query argument, so a case phenotype cannot shape the corpus. |
| reference-corpus-build | Build the versioned, phenotype-blind reference corpus that Routes B and C both consume. Retrieved rows come from local authoritative sources; a model-proposed row is promoted only when its accession appears in a retrieved source and a citation resolves. Unverified proposals stay in an exclusion ledger. |
| reference-corpus-compare | Compare two or more frozen corpora: per-field and per-species coverage (mutant phenotype, pathway, GO, evidence) against a baseline label, plus optional Route C Recall@K recomputed from frozen retrievals per corpus. Ranking weights are untouched and scores are never compared across corpora. |
| corpus-mask | Issue #237: remove the truth gene / its family / its close orthogroup from a frozen reference corpus BEFORE prediction, so a benchmark can measure recovery of genuinely unfamiliar biology. Mask rules come only from a benchmark-construction mask spec and are fingerprinted in the manifest; no prediction command ever sees them, only the masked corpus TSV. Refuses to write a masked corpus that still names the masked gene/family (self-check). |
| retrieval-index-build | Route C: build an immutable, versioned local reference index from source-backed plant records. Rows without a family/orthogroup or without a source id are excluded before Top-K and kept in a ledger; every artifact is SHA-256 fingerprinted. |
| retrieval-run | Route C: run a phenotype-only local retrieval against a frozen index and emit the retrieval contract. The command has no interval/atlas/candidate argument, so Phase 1 cannot see the QTL interval. |
| retrieval-merge | Route C: merge two or more per-database retrieval TSVs by within-database rank. Raw scores are never compared across databases; unannotated rows are excluded before Top-K selection. |
| retrieval-hybrid-merge | Route C: fuse two or more per-provider retrieval TSVs (e.g. lexical + lsa) by reciprocal-rank fusion over each provider's within-provider rank. Raw retrieval_score values are never compared or summed across providers -- only rank is read. |
| retrieval-import | Route C: import an external phenotype/protein retrieval TSV, validate it against the documented contract, and freeze it. No interval input exists on this command: the retrieval stays phenotype-only and interval-blind until the separate interval-match command. |
| retrieval-interval-match | Route C: match the frozen retrieval artifact against the QTL interval (interval atlas). Runs only after retrieval-import; the frozen retrieval is never modified here. |
| triple-route-compare | Compare all three independently generated routes (A: phenotype-first, B: interval-first, C: embedding retrieval). Convergence is recorded as state and adds no score. |
| scenario-create | Create a locus-level causal-scenario ledger from frozen route and evidence artifacts without changing rank or score. |
| scenario-report | Render the causal-scenario ledger as a reporting-only summary without changing rank or score. |
| tnt1-evidence-normalize | Normalize curated Medicago Tnt1 evidence while keeping insertion observation, causal support, genetic background, publication provenance, and cross-assembly identity separate; only native-background causal evidence with a complete publication ID enters the Route B adapter. |
| route-c-benchmark | Score frozen Route C retrieval artifacts on known-gene cases: Family Recall@5/@10/@20, Gene Recall@20, and A-only / C-only / A+C union. Interval-blind by construction. |
| route-c-direct-search | Route C Direct Rescue: phenotypeと明示species scopeからEurope PMCをinterval-blind検索する。retrieval hitは候補リンクでありsupportではない。 |
| route-c-direct-freeze | Route C Direct Rescue: agent/researcher proposalをgeneric source verificationへ通し、citation・direct claim・gene/locus mentionを確認してSHA-256 freezeする。 |
| route-c-direct-project | Route C Direct Rescue: freeze hashを再検証後、interval atlasと既存eligible BLAST hitsへ投影する。homologyだけでtrait supportは転移しない。 |
| route-c-provider-bait | Route C: 記号だけの凍結recordに、引用PMIDを参考文献に持ち同じgene名のUniProtKB entryを1件だけbaitとして選ぶ。0件・別配列の複数件はfail closed。区間は見ない。出力のprovider_identity_edges.tsvはprojectの `--identity-mapping` へ渡し、review-onlyのsource identityになる。 |
| route-c-direct-members | Direct Cのfreeze/projection lineageを検証し、同種の直接支持、異種の類似表現型review、ledger-onlyを分離する。 |
| route-c-direct-benchmark | Direct source recoveryとpost-freeze projectionを別分母で採点する。expected QTL ID / atlas SHA-256を照合し、別interval artifactの混入を拒否する。contaminated/not_evaluatedはmissにしない。 |
| route-c-provider-comparison | Route C: run route-c-benchmark once per provider (e.g. lexical/lsa/hybrid) over the same frozen known-gene cases and report Family/Gene Recall@K side by side. Never restructures route-c-benchmark's own scoring and never pools a provider's raw retrieval_score. |
| phenotype-variant-benchmark | Score frozen per-formulation lane artifacts for phenotype wording / description-granularity robustness: Route A family overlap, Route C Recall@K and retrieval overlap, Route B overlay stability with an atlas-fingerprint check, A/B four-states, Stage-A stability, and stable / rank-movement / catastrophic-family-dropout / stage-A-dropout classification. Never re-runs a lane and never changes a score. |
| variant-evidence | Normalize a parental VCF plan and compare SnpEff/VEP transcript consequences without treating predictions as Wet evidence. |
| multigenome-private | Classify focal events against multiple control accessions as private_confirmed / private_possible / shared_confirmed / not_evaluated. Keeps control callability and projection denominators explicit; never adds score and never promotes a candidate. |
| sv-plan | Plan assembly SV/PAV analysis and reconcile SyRI/read-based calls as unscored context. |
| promoter-review | Limited non-promoting promoter/5'UTR/3'UTR/SV review route after the coding hypothesis is examined. Emits strand-aware regions and intersects parental variants / SV-PAV-TE as unscored, review-required context; never adds to a gene score. |
| allele-compare | Call candidate variants from a pair of homologous FASTA sequences with minimap2. |
| wet-state | Normalize assay coverage/results and derive append-only candidate Wet states. |
| hypothesis-ledger | Track locus/mechanism hypotheses separately from gene rank. |
| candidate-budget-gate | Freeze a supported Wet-candidate set and request bounded adaptive refinement only when it exceeds the Wet budget. |
| llm-usage | Aggregate provider-reported token usage/call counts from llm-run raw responses without estimating missing tokens or prices. |
| route-b-semantic-bundle | After B0 freeze, compress the complete interval atlas into identity-preserving family/orthogroup/gene units and emit provider-neutral semantic prompts. |
| route-b-semantic-import | Import normalized llm-run responses for frozen Route-B semantic prompts. Unknown sources/reference genes are unsupported; all other rows remain proposals pending relation verification. |
| route-b-semantic-source-verify | Verify proposal citations/source relations with the existing source-verify engine while keeping focal-family transfer as a separate unresolved axis. |
| route-b-semantic-transfer-verify | Combine source-level mechanistic support with frozen family/orthogroup/copy identity. Group-level claims may verify without choosing one copy; interval_gene claims require one-to-one/RBH-like transfer. |
| route-b-semantic-candidates | Freeze Route-B semantic candidate members from verified claims before comparison. |
| semantic-claims-freeze | Validate and freeze structured Route-B semantic claims; only verified source-backed rows are narrowing-eligible. |
| supported-member-union | Build the non-scoring union of supported Route-A bridge members and final verified Route-B semantic members. |
| wet-candidate-handoff | Canonical path from frozen Route A/B artifacts to a Wet-reviewable candidate set and stage_a_narrowed. |
| adaptive-synthesis-request | Emit the bounded ambiguity-resolution request for an over-budget candidate set; refuses to run within budget. |
| adaptive-synthesis-import | Normalize reviewer ambiguity-resolution proposals; rank_contribution stays none and no proposal narrows the set. |
| adaptive-refinement-apply | Narrow a candidate set using machine-verified verdicts only; unsourced or not_evaluated verdicts never remove a candidate. |

以下は代表的なコマンドの詳細（入力・出力・役割）です。全コマンドの網羅的な一覧は上の表が担い、引数レベルの正本は [`CLI-REFERENCE.md`](CLI-REFERENCE.md) です。

| コマンド | 入力 | 出力 | 役割 |
|---|---|---|---|
| doctor | なし | JSON/text | WSLと必須・任意OSSの可用性 |
| extract | FASTA, GFF3, interval | 候補TSV、genomic/CDS FASTA、manifest | QTL区間だけを物理的に抽出 |
| peak-context | candidate + LODプロファイル/peak metadata/`--peak`/`--interval` | peak_context TSV、canonical position evidence、manifest | QTL peak・support interval・組換えcoreを段階付き位置証拠へ取込 |
| trait-ontology-map | trait terms TSV + ローカルontology synonym TSV | 正規化trait-ontology TSV | 外部API不要でPTO/PO/CO IDを付与、match_typeと承認状態を保持 |
| trait-terms build | raw trait terms + 任意のontology/manual synonym/LLM proposal | 承認済みterm TSV、audit、direct/context語リスト | 語を分類・承認gateし、benchmarkでanswer-leakageを拒否 |
| reference-candidate-table | candidate + source-backed reference gene panel + provider別orthology TSV | candidate-reference対応表、review表、manifest | 自由形質語から別途作った既知gene panelをアンカーにし、候補とのorthology・annotation・paralog曖昧性を人が読む。provider複数一致は頑健性であり独立票ではない |
| auto-triage | candidate/protein/任意id-mapまたはFASTA/GFF3/interval + 自然文形質 + role付きreference resource panel + 任意VCF | 全候補ortholog、reference機能、形質match、variant、rank、固定20件、canonical identity/provider/coverage/lineage/state | curated gene panelなしのcandidate-first経路。fingerprint一致時だけresumeし、任意LLMは形質機構語だけを追加 |
| homolog-panel-auto | 自由形質文 + 種 + 任意LLM provider | UniProtで検証したpending-review panel、FASTA、raw API、manifest | LLMなしは固定語形展開、LLMありは検索仮説を追加。候補ID・座標・配列は送信しない |
| homolog-panel-match | 自動panel + 区間candidate/protein + 任意base evidence | interval-only BLAST evidence、rank、固定20件 | 自動panel証拠をlow/mediumの1 comparative familyとして統合。orthology・因果を断定しない |
| homolog-panel-benchmark | frozen known-gene benchmark manifest | case別fixed-20 recall、truthへの直接homolog support、failure、manifest | primaryはraw trait textのみ。curated termsは`--benchmark-curated-terms`を明示したablationに限定。provider/model/prompt fingerprint一致時だけresume |
| primerblast-import | PrimerBlast OSS JSON（design/CAPS/InDel）+ 任意candidate | primer_candidates TSV、KASP/ASP/HRM scaffold、manifest | primer pair・specificity rank・off-target・CAPS enzyme・riskを正規化し親/F1 validation gateを保持 |
| orthology-tier1-plan | candidate + protein + focal proteome + reference panel | Tier1検索script、manifest | profile(accuracy64/laptop16)でthread・hit上限・tree方針を制御、queryは候補のみ |
| focal-paralog-import | focal自己検索BLAST + candidate | paralog evidence、paralog_summary | QTL区間外paralogを検出しambiguity_reasonを出力、best hitを1:1と断定しない |
| orthofinder-import | OrthoFinder Orthogroups.tsv + candidate | orthogroup evidence、summary | 1:1/1:many/many:many/lineage-specificを分類しparalog数を記録 |
| omamer-import | OMAmer HOG table + candidate | orthogroup evidence | HOG配置をcanonical evidenceへ変換 |
| synteny-import | jcvi .anchors / MCScanX .collinearity + candidate | synteny evidence | シンテニーブロックの相手遺伝子をsynteny_supportedへ |
| exon-coverage-check | exon_sequencing_plan + PrimerBlast tile JSON | coverage report、gap一覧、retry plan | 各CDS windowがamplicon被覆されたか検証しgapに緩和retryを計画 |
| annotation-keywords | candidate、BLASTP、相同遺伝子GFF、trait keyword | canonical evidence、全best-hit監査表 | AIなしで相同注釈と形質語を照合 |
| allele-compare | homologous candidate FASTA 2本 | minimap2 PAF、variant TSV、任意のmarker候補TSV | 上位候補だけassembly/親候補間で比較 |
| benchmark-reference | 1 referenceのFASTA/GFF/protein/eggNOG + known-gene/homolog config | profile別rank、Top-k recall、manifest | answer IDを採点後にだけ使う約10 Mb virtual-QTL benchmark。`--resume`対応 |
| screening-shortlist | ranked candidate TSV | primary 20件、Rescue appendix、境界同点監査、manifest | LLM・人手評価の前に同じStage A集合を固定。`top_eligibility`は検索gateにしない |
| evidence-tiers | evidence-rank/candidate-pipelineのranked_candidates.tsv(+任意knowledge_coverage) | target_size非依存のbiological evidence tier(strong/plausible/unresolved/contradictory/knowledge_limited)、review_budget/assay_capacityのoperational cut、manifest | operational cutはbiologicalな除外を意味しない。tierとcutは別列 |
| evidence-tier-benchmark | case別manifest(evidence-tiers/evidence-rankの出力) | Recall@budgetとevidence-tier recallの並置、interval-sizeバンド別layer化、manifest | 凍結済み出力を読むだけで再ランクしない |
| weight-policy-benchmark | frozen known-gene benchmark manifest | current/uniform/tier/eligibility比較、構造ablation、Wet-budget、bootstrap CI | ordinalな実験順序の検証であり因果確率ではない |
| mapping-gate | mapping metrics TSV | mapping status、review不足、machine-readable gate | 候補順位より前にlinked peak、marker密度、recombinant、低組換えを判定 |
| next-experiment | deterministic rank + case manifest + 任意mapping/assay/SV/recombinant状態 | action options/decision/dependency/review artifact | 生物学的順位を変更せず、marker、Sanger、targeted sequencing、SV、expressionの次実験だけをordinalに比較。人承認pending |
| wet-state | assay registry + 実測coverage + parent/F1/F2/recombinant result + 任意previous history | assay/coverage/result/candidate/history/reviewの6表 | assay失敗・partial coverage・範囲限定陰性を分離。coding negativeでもregulatory/SVを保持し、禁止遷移を履歴へ残す |
| hypothesis-ledger | locus/mechanism proposal + canonical evidence/Wet/SV ID + 任意history/rank | hypothesis/evidence-link/experiment-link/history/reviewの5表 | gene rankと独立。coding negativeはmechanismだけweakenしregulatory/SVへ分岐。LLM proposalはpending review |
| annotation-rescue | target FASTA + native candidates + source-backed query protein/metadata + 任意miniprot GFF | 全projection ledger + bounded review + manifest + local script | QTL区間限定。provider failureとunmappedを分離し、rank/fixed20を変更しない |
| reference-panel | native/close_anchor/functional_reference/outgroup roleのresource宣言TSV | 検証済みpanel + availability/skip理由 + provenance manifest | nativeのみ必須。optional欠落と外部APIをmachine-readableに分離。known-gene panelと別物 |
| external-rank | native ranked_candidates + 外部rank long TSV | coverage-aware ledger + coverage/discordance/contribution表 + manifest | 外部生スコアをdry_scoreへ加算しない。RRF consensusはscreening view。未評価と低順位を分離 |
| cross-crop-benchmark | case別config(native rank + 任意external/panel/rescue + known gene) | case/summary/examples/coverage表 + manifest + report | 作物群×resource条件で過適合を検証。uncoveredを最下位化せず、meta非改善も報告 |
| propagate-terms | 区間候補 + 任意orthology + 任意reference機能注釈(+任意trait語) | 伝播語ledger + 候補別検索語彙 + manifest | 参照ゲノムのオーソログ機能語を決定論的に伝播。LLM非依存・区間限定・native単独でも完走 |
| robustness-benchmark | base candidates + evidence + causal gene | 条件別safety表 + summary + manifest + report | 注釈劣化・paralog/tandem/decoy/truth欠損下でrankゲートを駆動し、誤sole-Top昇格の有無と安全保留を検証 |
| discovery-benchmark | case別config(candidates + evidence + answer gene) | leave-answer-out honest recall + leakage表 + manifest | 答えを記述する証拠を除いて再ランクし、真の識別力とleakage寄与を測る。weight不変・honest recallが主指標 |
| query-eligibility | homolog query manifest(+任意query FASTA) | eligible/excluded query表 + provenance監査 + manifest | Pisum由来query(別assembly/projection/renamed/answer-backmapped)を実行前に禁止し、非Pisum機能ホモログのみ許可。strictでfail |
| llm-shortlist-prompt / llm-shortlist-import | Stage A 20件のresearch packet / model response | 同じ20件だけの任意Stage B順位、candidate-set監査 | LLMは集合を追加・削除せず、未reviewの読み順だけを返す |
| llm-review-benchmark | 同じStage A集合とLLMなし/あり順位、既知truth | profile別Top-k・truth rank | 小規模Stage B pilotを本体の非LLM benchmarkと分離 |
| evidence-term-promote | canonical evidence + 承認済み機構語 | low-confidence trait match + audit | LLM不要。sourceを新造せず、広い位置情報を機能証拠へ変換しない |
| annotate-plan | candidate TSV + protein FASTA | candidate protein FASTA、GO台帳、局所annotation script | BLAST/DIAMOND/HMMER/eggNOGをcandidateだけへ適用 |
| evidence-build | candidate TSV + 注釈/発現/位置/手動証拠 | canonical evidence TSV | 異種OSS出力を同じ証拠形式へ統合 |
| evidence-rank | candidate TSV + evidence TSV | rank、加減点、重み、分岐JSON | source-backed evidenceだけで一次/最終順位を計算 |
| orthology-plan | 一次rank + proteome | 上位候補だけのreciprocal BLAST script | 全候補への重い計算を避ける |
| orthology-import | reciprocal BLAST結果 | canonical evidence TSV | RBHとparalog ambiguityを再順位へ戻す |
| candidate-pipeline | 上記入力一式 | stage別出力とmanifest | 一次順位から上位精査・最終順位までを一貫管理 |
| orchestrate | candidate TSV + 利用可能なDB・ゲノム | 工程別script、manifest、run script | OSS annotationとPrimerBlast OSS工程をまとめて準備 |
| primer-plan | stage別入力 | PrimerBlast OSS実行script、plan TSV、manifest | sequencing primer、variant marker、QTL grid marker |
| literature | 検索語 | Europe PMC JSON | 認証不要のread-only文献API |
| research | candidate TSV、種、形質 | API文献summary/evidence TSV | AIなしの候補リサーチ |
| research-packet | candidate TSV、種、形質 | 検索表、任意AI用prompt | LLMを使う場合だけの補助 |
| retrieval-merge | DBごとのretrieval結果TSVを2つ以上 | DB内rankのround-robin統合 + excluded ledger + status | DB間raw scoreを比較せず、family/orthogroup未注釈行にTop-K slotを消費させない |
| retrieval-import | 外部retrieval結果TSV（contract11列）+ 任意phenotype | freeze済みretrieval TSV + status | Route C: phenotype/protein embedding retrieval をvalidateしてfreeze。**interval引数は存在せず**、contract外の列と未監査のmulti-DB表は拒否 |
| retrieval-interval-match | freeze済みretrieval + interval atlas | interval照合TSV + status | freeze後にのみinterval照合。freeze済みartifactは書き換えない |
| triple-route-compare | route A/B仮説 + route C freeze | 3経路比較TSV + convergence JSON | A/B/Cのconvergenceを状態として記録。scoreは一切動かさない |
| route-c-benchmark | case別manifest（retrieval + truth family/gene + 任意route A） | Family Recall@5/10/20・Gene Recall@20・A/C/union表 + report | phenotype-only / interval-blind のretrieval性能を計測。causal evidenceにしない |
| discovery-chain-benchmark | case別manifest（phenotype + truth family/gene + 各laneの結果artifact TSV + 任意atlas/rescue/orthology/reference-function/literature/interval candidate） | lane別metrics + failure taxonomy + **joint A/B 4-state** + **Route B funnel** + **narrowing-efficiency** + **negative controls** + manifest + report | A/B discovery chain全体をend-to-endで集計。lane自体は実行せず、外部実行済みartifactを入力としてスコアする。`control_type` で shuffled / decoy interval / lexical regression のnegative laneを同レポートに混載できる |
| ranking-qa | candidate + evidence（+任意weights/top-k/perturbation/source-dropout/provider-dropout） | rank stability / source dropout stability / weight perturbation summary + manifest | canonical rankを変更せず、weight・source-group・providerのdropout/perturbation下でのrank変動を診断。score差はprobabilityではない |
| ranking-stability-benchmark | frozen case manifest（candidates/evidence/ranked_tsv/answer_gene_id） | stability benchmark cases/summary/manifest/report | frozen prediction artifactのtruth geneがperturbation・dropout下でTop-K内に残るかを報告。入力欠落・truth不在はnot_evaluated |
| next-experiment-benchmark | case別manifest（expected action class + ranked/case manifest + 任意mapping/assay/SV/recombinant artifact） | exact-action / appropriate-escalation accuracy + blocked検出 + report | next-experimentの決定層を評価。blocked≠不適切、missing provider≠反証。ordinalのまま |
| phenotype-variant-benchmark | case別manifest（truth + variants.tsv + variant別lane artifact） | Route A/C overlap・Route B overlay stability・4-state・Stage-A stability + manifest + report | phenotype言い換え・粒度へのロバストネスを計測。variantはtruth前にfreeze、答え用語の注入なし、atlas fingerprint不変を検証 |
| dry-prioritize | Dry evidence TSV | score内訳、Top/Reserve/Rescue群付きTSV | Wet情報を拒否し、親を読む順番だけを決める |

ZW6のようにGFFとFASTAの染色体IDが異なる場合は `extract --seqid-map GFF_ID=FASTA_ID` を使う。multi-Gb FASTAは既定でmetadataとFAI hashを記録し、全配列SHA-256が必要な場合だけ `--full-hash` を加える。

`allele-compare --reference-candidates candidates.tsv` を指定すると、gene-oriented alignment座標をゲノム座標へ戻した `marker_candidates.tsv` も作る。assembly比較だけで得た変異は常に `unvalidated` であり、実親のSanger/reads確認前には `primer-plan --stage variant-marker` を通過しない。

複数アラインメントは構造として分類する。混在strand（`inversion_suspected`）、query断片が分離する`split_alignment`、同一query領域が複数locusへ載る`duplication_suspected`（`ambiguous_primary_alignment`）、片側だけ整列する`presence_absence_or_large_indel`を検出し、これらやpartial・reverse best alignmentではmarkerを出さず、理由を`marker_suppression_reasons`とmanifestへ残す。complex eventを単純SNV/InDelへ誤分類しない。InDelはreference配列に対しleft-normalizeし（`left_normalized`列）、`--min-coverage`でcoverage閾値、`--sequence-type genomic|cds|transcript`で配列種別（frame-aware判定）を明示できる。

候補順位のbenchmarkはaccession pair比較ではなく、`benchmark-reference`で1 referenceずつ行う。`screening_rank`は機能的な怪しさ、`dry_group`は独立証拠の十分性であり、同じ意味ではない。
| rank | candidate ledger TSV | ranked TSV | 既存の全工程台帳用legacy ranking。新規のWet前順位付けには使わない |
| assay-plan | ranked TSV | assay plan TSV + script | 旧PrimerBlast OSS接続。新規利用はprimer-planを推奨 |
| cosegregation | F2 TSV | genotype別表現型要約TSV/JSON | Wet結果を仮説更新に使う |

## 一括オーケストレーション

~~~powershell
.\scripts\candidate-qtl.ps1 orchestrate `
  --candidates C:\private\qtl\dry_candidates.tsv `
  --out-dir C:\private\qtl\orchestration `
  --protein-fasta C:\private\reference\proteins.fa `
  --diamond-db C:\private\db\arabidopsis.dmnd `
  --hmmer-db C:\private\db\Pfam-A.hmm `
  --primerblast-root /path/to/primerblast-oss `
  --genome C:\private\reference\genome.fa `
  --gff3 C:\private\reference\genes.gff3 `
  --nucleotide-db C:\private\db\genome `
  --candidate-genomic-fasta C:\private\qtl\extract\candidate_genomic.fa `
  --interval Chr1:100000-500000 `
  --top-n 8 --num-threads 4 --json
~~~

入力が存在する工程だけscriptを生成する。親間変異TSVをまだ渡していない場合、選抜marker工程はmanifest上でdeferredとなる。

## 旧自由記述表からのWet前Dry順位付け（互換用）

新規解析は`candidate-pipeline`とcanonical `evidence.tsv`を使う。既存の`templates/dry_candidate_prioritization.tsv`を再利用する場合だけ、以下の`dry-prioritize`を使う。

~~~powershell
.\scripts\candidate-qtl.ps1 dry-prioritize `
  --candidates C:\private\qtl\dry_candidate_evidence.tsv `
  --out C:\private\qtl\dry_candidates.tsv `
  --top-n 7 --reserve-n 3 --json
~~~

出力の `score_breakdown` が各加減点を示す。`rescue_reason` がある候補は知識スコアと別の `rescue` 群へ残る。親間変異、marker validation、F2共分離などのWet列に値がある入力はエラーにし、Dry順位への情報漏洩を防ぐ。

全候補は `top`／`reserve`／`rescue`／`not_selected_yet` のいずれかに分類され、`selection_reason` 列に選択・未選択の理由を明示する（`not_selected_yet` は原因候補からの除外ではない）。`top`・`reserve` には親間direct sequencing用の `recommended_regions`（全CDS・exon-intron junction・UTR・機構上妥当ならpromoter）を付す。既存のWet台帳・`cosegregation` は独立した後段として維持する。

## PrimerBlast OSSの段階別利用

### 親間direct sequencing primer

~~~powershell
.\scripts\candidate-qtl.ps1 primer-plan `
  --stage direct-sequencing `
  --candidates C:\private\qtl\dry_candidates.tsv --top-n 8 `
  --out-dir C:\private\qtl\primers\sequencing `
  --genome C:\private\reference\genome.fa `
  --gff3 C:\private\reference\genes.gff3 `
  --candidate-genomic-fasta C:\private\qtl\extract\candidate_genomic.fa `
  --db C:\private\db\genome `
  --primerblast-root /path/to/primerblast-oss `
  --json
~~~

GFF3-aware gene assayに加え、genomic FASTAを渡すと長い遺伝子を覆うtiling primerも生成する。

### 確認済み親間変異からCAPS/InDel選抜marker

~~~powershell
.\scripts\candidate-qtl.ps1 primer-plan `
  --stage variant-marker `
  --variants C:\private\qtl\parental_variants.tsv `
  --candidates C:\private\qtl\dry_candidates.tsv --top-n 8 `
  --out-dir C:\private\qtl\primers\selection_markers `
  --genome C:\private\reference\genome.fa `
  --vcf C:\private\qtl\parents.vcf.gz `
  --db C:\private\db\genome `
  --primerblast-root /path/to/primerblast-oss `
  --json
~~~

入力は `templates/parental_variants.tsv` に合わせる。SNPはCAPS候補、InDelはサイズ差marker候補として処理し、JSONとprimer発注用order tableを出す。

### PrimerBlast OSS結果の取込（primerblast-import）

`primer-plan` が生成・実行したPrimerBlast OSSのJSON（design report・CAPS/InDel/assay report）を `primerblast-import` で `primer_candidates.tsv` へ正規化する。primer pair、`specificity_rank`（off-target少・gel区別可・低risk順）、`n_off_target_total`、CAPS enzyme・`caps_min_gel_gap`・`caps_distinguishable`、複数DBの `offtarget_<db>` 列を保持する。marker typeは `caps`／`dcaps_needed`／`indel_size`／`sequencing_or_pcr` を区別し、natural CAPSが無いSNPは `dcaps_needed` として dCAPS/代替markerへ回す。

予測risk（`predicted_risk`）と実験結果（`parent_a_result`/`parent_b_result`/`f1_result`/`experimental_status`）を別列に分け、`validation_status` は既定 `unvalidated`。全F2展開の前に親・F1検証を要求するgateを保つ。`--candidates` でcandidate ledgerへgene IDをjoinし、`--export kasp|asp|hrm` でallele-specific発注scaffoldを出せる。

~~~powershell
.\scripts\candidate-qtl.ps1 primerblast-import `
  --results C:\private\qtl\primers\selection_markers\results\GeneA.marker.json `
  --candidates C:\private\qtl\dry_candidates.tsv `
  --out-dir C:\private\qtl\primers\imported --export kasp --json
~~~

### QTL区間を粗く狭めるgrid marker

~~~powershell
.\scripts\candidate-qtl.ps1 primer-plan `
  --stage qtl-grid `
  --interval Chr1:100000-500000 --n-markers 8 `
  --out-dir C:\private\qtl\primers\grid `
  --genome C:\private\reference\genome.fa `
  --db C:\private\db\genome `
  --primerblast-root /path/to/primerblast-oss `
  --json
~~~

## 実験上の制約

生成scriptは計算案であり、Wet検証を置き換えない。

- direct sequencing primerは両親とF1で単一増幅を確認する。
- marker化する変異はSanger traceまたは十分な親VCF evidenceで確認する。
- 選抜markerは両親・F1・代表F2で検証してから全F2へ展開する。
- 複数referenceがある場合はDBを反復指定し、multi-reference specificityを確認する。

詳細は ORCHESTRATION.ja.md と [LOCAL-ANNOTATION.ja.md](LOCAL-ANNOTATION.ja.md) を参照する。

## JSON契約

成功時は結果オブジェクトをstdoutへ出す。失敗時はJSONを指定した場合でもエラーをstderrへ出す。FASTAやVCFの内容、秘密情報、未公開表現型はstdoutへ出さない。

`--json` 指定時のエラーは stderr の1行目に機械可読なオブジェクトで出る。

~~~json
{"error": "interval must be CHROM:START-END with 1-based inclusive coordinates", "command": "extract"}
~~~

## 終了コードとエラー表示

| code | 意味 |
|---:|---|
| 0 | 成功 |
| 1 | 想定内の失敗。入力が拒否された、ファイルが無い、gateが止めた、外部toolが非ゼロで終了した |
| 2 | candidate-qtl 側の不具合。入力の問題ではない |

想定内の失敗(1)は**tracebackを出さず**、修正方法を含む1行だけを stderr に出す。

~~~text
$ candidate-qtl extract --fasta genome.fa --gff3 genes.gff3 --interval chrZ:1-900 --out-dir out
error: GFF sequence ID 'chrZ' maps to FASTA sequence ID 'chrZ', which is absent;
available examples: chr1. Supply --seqid-map GFF_SEQID=FASTA_SEQID.
~~~

不具合(2)の場合だけ traceback を残す。想定内の失敗でも traceback が要るときは `--traceback`
または環境変数 `CANDIDATE_QTL_TRACEBACK=1` を使う。

~~~bash
candidate-qtl --traceback extract ...    # flag は subcommand より前
CANDIDATE_QTL_TRACEBACK=1 candidate-qtl extract ...
~~~

`guided-profile` のように失敗を例外にせず audit manifest へ記録するcommandは、manifest を
書き出したうえで結果の `exit_code` を終了コードにする(該当節を参照)。

Python APIとしての `candidate_qtl.cli.main()` は従来どおり例外をそのまま送出する。終了コードへ
変換するのはプロセス入口の `run()` で、console script `candidate-qtl` と
`python -m candidate_qtl` の双方がこれを使う。

## Guided profile

`guided-profile` は既存commandを削除せず、通常運用の入口をまとめます。profileは `screening` / `full-case` / `audit-benchmark` / `audit-predict` / `triple-route`。
`audit-predict` / `triple-route` のRoute Aは `--phenotype-localization-payload` を先にfreezeし、目標12件のwell-formedな仮説をimportした後にcausal-layer coverageを通す（行数だけでは拒否しない）。補充再開は `--gap-fill-mechanisms`、`--gap-fill-mechanism-import-manifest`、`--gap-request-audit` の3点を必須とし、coverage不足のportfolioはreference verificationへ進まない。

- `screening`: FASTA/GFF3/QTL/phenotypeから予備shortlistまで。
- `full-case`: screening昇格またはreview済みcase-run。
- `audit-benchmark`: `audit-predict` の後方互換alias。truth/answer列を拒否し、answerless predictionをfreezeする。採点は別processでのみ行う。
- `audit-predict`: truth入力なしで prediction だけを freeze するlane（`audit-score` で別プロセスから採点）。
- `triple-route`: A（phenotype-first）/ B（interval-first）/ C（retrieval）の3経路を独立laneで生成し、最後の `triple-route-compare` だけで合流。Bは **B0（phenotype-free atlas freeze → reverse lookup）→ B1（frozen atlasへのphenotype overlay）→ reverse-literature** の順で走り、atlasはoverlayより前にfingerprint（`atlas_sha256`）されます。reverse-literatureはB0のlookup出力から自動chainingされ、offlineでcacheが無ければ `not_evaluated_offline` を記録して完了します（`--reverse-lit-cache` でquery-keyed cacheを指定可能）。Route Aのsource/accession検証を許可する`--allow-external-api`だけでは、診断専用のlegacy B1外部検索は有効になりません。B1からEurope PMCを検索する場合だけ`--allow-reverse-literature-api`を追加します。compare出力にはA/Bの4状態（`a_and_b_supported` / `a_only_supported` / `b_only_supported` / `neither_supported`）と `literature_link_candidate_count` が載ります。B1は `--phenotype` が無いと `not_runnable_missing_input` になり、失敗扱いしません。

既定はdry-runです。`--execute`時もprofile固有のrank logicは持たず、missing inputはmanifestへ残し、review-required項目は停止条件として扱います。

provider選択は `discover` と同一です。`screening` はローカルprotein BLAST DBを既定で自動探索し、`--protein-db` で明示指定、`--no-auto-local-db` で無効化できます。選択方針(`provider_policy`)と実際に使われたprovider(`annotation_provider`)は両方 `guided_profile_manifest.json` に残ります。

`--execute` 中にstageが失敗しても manifest は必ず残ります。`audit-benchmark` / `audit-predict` の各laneは独立で、先頭laneが失敗しても入力の揃った後続laneは実行されます。どちらもtruth入力を受け付けず、予測をfreezeした後に別processでのみ採点します。`screening` / `full-case` はfail-fastですが、停止前にplanned commandとresolved argumentsを含むmanifestを書き出します。stage別に `completed` / `failed`(例外型・message・開始終了時刻)/ `not_runnable_missing_input` / `skipped_after_failure` が記録され、全体は `overall_status`(`completed` / `partial` / `failed` / `interrupted` / `dry_run`)になります。失敗があった場合の終了コードは1です。

## Interval atlas（phenotype非依存）

`interval-atlas` は、**phenotypeを一切使わずに**QTL区間内の全遺伝子を family / domain architecture / orthogroup / broad function / conservation / evidence completeness へ整理します。正しいbaitを先に思いつく必要を減らし、区間内に実在する有限の選択肢から仮説を立てるためのroute（`family_review_queries.tsv` はfamily・process単位のクエリを出します）。

分類はあくまで探索indexで、causal evidenceではありません。全出力の `rank_contribution` は `none` で、unknown / conflicting / dark な遺伝子も独立bucketとして残り、1件も落ちません。broad function binは版管理された `templates/broad_function_bins.tsv`（wheelにも同梱）で、どのパターンが当たったかが遺伝子ごとに残ります。manifestは `phenotype_used: false` と `atlas_sha256` を記録するので、後段のphenotype overlayがatlasを書き換えていないことを監査できます。詳細は [INTERVAL-ATLAS.ja.md](INTERVAL-ATLAS.ja.md)。

`knowledge-coverage` は「この候補が良いか」ではなく「**そもそも何か言えたのか**」を候補ごとに記録します（annotation / family / homolog / reference function / literature / ontology / expression / variant）。`available` / `knowledge_unavailable`（資源を見たが何も無い＝unknown、unlikelyではない）/ `not_evaluated`（providerを渡していない）の3状態を潰さず、知識不足が制約の候補は**降格ではなくエスカレーション**（`escalate_knowledge_limited_candidate`）にします。`benchmark-stratify` は任意のper-caseベンチマーク表を層別（知識tier・系統特異性・区間サイズ帯など）に再集計し、層を平均して1つの数値にしません。詳細は KNOWLEDGE-COVERAGE.ja.md。

## 順引き / 逆引きの2経路

仮説へ至る道は2つあり、**壊れ方が鏡像**です。順引き（`phenotype-first-status`）は表現型から遺伝子/familyを決め打ちして区間で探すので、**区間に無い遺伝子を捏造する**リスクを負います。逆引き（`interval-reverse-lookup`）は区間に実在する遺伝子から出発するので捏造はしませんが、**区間にある遺伝子を後付けで正当化する**リスクを負います。両方走らせて `dual-route-compare` で突き合わせると、存在は区間に、関連性は表現型に検証されます。

逆引きの参照はアラビ限定ではありません。`--reference-function` は `reference_gene_id`（他種ホモログ）でも `gene_id`（**焦点リファレンス自身の注釈**）でも引けます。LLMを使う場合、promptには区間に実在する遺伝子の注釈済み正体だけを並べ、一覧に無いものを名指しした提案は `rejected_unanchored` としてbaitへ進めません。「この区間からは支持できる仮説が無い」も有効な答えとして明示しています。

生成中は2経路を混ぜません（混ぜるとinterval-informedな仮説をblind発見性能として集計してしまう）。合流は `dual-route-compare` だけで、`mechanistic_convergence` / `family_or_subfamily_convergence` / `gene_overlap_only` / `process_related` / `text_only_match` / `phenotype_only_supported` / `interval_only_supported` / `discordant` / `both_weak` を区別し、convergenceはscoreを一切動かしません。同じinterval geneを共有するだけの一致は `gene_overlap_only` であり、strict convergenceにはなりません。詳細は [DUAL-ROUTE-DISCOVERY.ja.md](DUAL-ROUTE-DISCOVERY.ja.md)。

### 複数assemblyの候補集合とRoute B圧縮

精度の低いassemblyへ全候補を強制投影しません。`multi-reference-locus-union` は、QTL座標・marker・syntenyなどで各assembly側の対応locusを先に定義し、各区間の遺伝子を和集合として保持します。verified 1:1だけを同一entityへ統合し、one-to-manyは全copy、alternate-only/PAVは元assembly IDのまま残します。したがって、V2にしか存在しない遺伝子も、別assemblyに写らないことを理由にcandidate universeから消えません。

`candidate_universe.tsv` はrecall評価用の全候補、`compression_candidates.tsv` はverified 1:1の重複だけを除いたRoute B入力です。この2つの件数・truth retentionは別々に測定します。`assembly-bridge` のfam3はfamily-levelの圧縮・照合には使えますが、fam3一致だけでsingle-copy identityを確定しません。family不一致、synteny conflict、unresolved mappingはreview/not_evaluatedとして保持します。

Type 1–3のmaskは候補の存在を隠す操作ではなく、prediction前の知識edgeを隠す操作です。`benchmark-knowledge-mask` のType 2（`leave_truth_semantic_edges_out`）は、truth固有のfunction / phenotype / pathway / source edgeだけをRoute Bのreference-functionから除き、assembly mappingとorthology/familyの中立な投影は保持します。Type 1/3はexact gene → family/orthogroup → declared sequence homologの順に入れ子で機能・対応関係を除きます。いずれもcandidate-universe/presence行は削除せず、Type 2のsemantic除去が0行なら`not_evaluated`です。これによりType 3でも「別assemblyに正解がある」ことと「既知homolog知識を使った」ことを混同しません。

### Route B semantic investigation & candidate freeze

Route B の詳細調査は、B0（phenotype-free interval atlas）を凍結した後に開始します。

1. `route-b-semantic-bundle`: 凍結 atlas を family/orthogroup/gene 単位の bounded unit に集約し、provider-neutral な prompt を生成（生 gene 列挙の dump を禁止）。`--lookup` に同じ B0 run の `interval_reverse_lookup.tsv` を渡した場合は、source-backed function/pathway/mutant/source ID も凍結 unit へ保持する。B0 coverage manifest（`assigned` / `unassigned` / `unresolved` / `failed` の厳密合計一致）を出力。
2. `route-b-semantic-import`: LLM 応答から family/subfamily/species/source/evidence-scope 付き structured claims を取り込み（初期状態は `proposed`、未認識 source/gene は `unsupported`）。
3. `route-b-semantic-source-verify`: 文献 source/citation を local cache / Europe PMC で機械検証し、Europe PMCで解決できない場合はNCBI PubMed/PMC E-utilitiesへcitation IDだけを送ってfallbackする。provider failureはmiss/unsupportedに変換しない。
4. `route-b-semantic-transfer-verify`: source 支持と区間ホモロジー（group-level vs 1:1 RBH copy-level）を突き合わせ、4値（`verified` / `unsupported` / `contradicted` / `not_evaluated`）を判定。
5. `route-b-semantic-candidates`: 検証済み（`verification_status=verified`）claim の支持 member を `route_b_verified_candidates.tsv` として比較前に凍結（複数 copy は単一遺伝子へ強制縮退させず保持）。
6. `semantic-claims-freeze`: structured claims の fingerprint を固定し、`semantic_claims_frozen.tsv` を出力。

### Wet候補への接続と bounded adaptive synthesis

A/B を凍結したあと、Wet reviewer が扱える候補集合まで運ぶ経路です。

1. `wet-candidate-handoff`: 凍結済み A/B artifact から supported member union を作り、Wet budget（既定20）と突き合わせ、必要な場合だけ次段を要求し、budget 内なら `integrated-narrowing` まで進めて `stage_a_narrowed` を出す。`ranked_candidates.tsv` と `stage_a/` は入力であり、書き換えません。
2. `adaptive-synthesis-request`: 候補数が budget を超えたときだけ発動。**どの遺伝子を残すかは尋ねず**、どの family/subfamily ambiguity を解けば安全に縮むか、どの追加資料が判別に効くか、どの competing explanation を残すべきかだけを問う。budget 内の集合に対して実行しようとすると fail-closed で拒否します。
3. `adaptive-synthesis-import`: 回答を構造化して取り込む。`rank_contribution=none` / `review_required=yes` 固定で、score 列を持つ提案と、要求に含まれない遺伝子を挙げた提案は拒否します。provider が無ければ `not_evaluated_no_provider` として空表を出し、workflow は完走します。
4. `adaptive-refinement-apply`: **machine 検証済みの verdict だけ**で候補を削る。`support_falsified` は `verification_source` が無ければ拒否、`not_evaluated` は候補を残します。

1 round で budget 内に収まらない場合は、無理に Top-N へ潰さず `review_required_over_budget` で停止します。完了判定は exit code ではなく `handoff_state` です。

## 第3の経路: Route C（embedding retrieval）

phenotype を embedding した類似度検索を、A/Bと独立した第3 discovery routeとして扱えます。複数DBは `retrieval-merge` でDB内rankを公平に統合し、その後 `retrieval-import` でfreezeします。**retrieval生成時はQTL intervalを一切入力しません**（interval引数が存在せず、contract外の列は拒否）。`retrieval-interval-match` はfreeze済みartifactにだけintervalを照合し、最後に `triple-route-compare` でA/B/Cのconvergence状態を記録します。一致はscoreに加算せず、causal evidenceにもしません。詳細は [ROUTE-C-RETRIEVAL.ja.md](ROUTE-C-RETRIEVAL.ja.md)。

A/B/C が共有する参照コーパスは `reference-corpus-fetch`（UniProtKB から curated function / disruption phenotype / pathway / evidence code を taxon 単位で取得。自由文クエリ引数は存在しません）→ `reference-corpus-build --from-uniprot` → `reference-corpus-compare`（フィールド別・種別カバレッジと Recall@K の差分）で作ります。**行数はカバレッジではありません** — description だけの39,254行コーパスは fasciation クエリで CLAVATA を1件も返しませんでした。詳細と実測値は REFERENCE-CORPUS.ja.md。

## unplaced scaffold の救済

minor crop や draft assembly では、原因遺伝子（またはその一部）が unplaced / unanchored scaffold、haplotig、tandem familyの分離scaffold、repeat/gap隣接配列に載っていることがあります。chromosome上の区間だけをcandidate universeにすると、**assembly品質の問題で真の候補を最初から落とします**。

`scaffold-rescue` は、QTL locusへ**plausibly linked**なscaffoldを保守的に救済し、primary interval とは**別のcandidate class**として保持します。

| 支持の種類 | 強さ |
| --- | --- |
| paired marker（両端が同一scaffoldへunique） | strong |
| synteny bridge（近縁種の保存blockが native gap を橋渡し） | strong |
| Hi-C / linkage map / optical map / assembler placement | strong |
| 単独marker hit | weak |
| orthogroup / gene-content context | weak |

独立strong 2つ以上 → `scaffold_rescue_high`、strong 1つ → `medium`、weakのみ → `weak`、repeat由来などambiguous → `unresolved_scaffold`、支持ゼロ → `scaffold_rejected_as_unlinked`。

安全原則:

- **unplacedであるだけでは候補に加えません。** 何らかの積極的なplacement証拠が要ります
- **family/homology一致だけでは high-confidence にしません**（そのscaffoldが何を含むかは言えても、どこに属するかは言えない）
- 複数scaffoldへ当たるmarkerは**repeatとして棄却**します
- **placement confidence と causal confidence は別列**です
- primary interval の遺伝子を rescue 遺伝子で置換・並べ替えしません
- placement ambiguity は1本に決めず保持します
- **救済ヒットが無いことは「unplaced配列上に原因が無い」という negative evidence ではありません**

## native reference が不十分な場合の comparative locus reconstruction

`comparative-locus-reconstruct` は flanking marker sequence と、事前計算済みの複数近縁種 marker projection / gene-order / orthogroup cache から expected locus を再構成します。focal assembly・unplaced scaffold・transcriptome の事前検索結果を取り込み、`native_annotated / scaffold_recovered / transcript_supported / comparative_expected_missing / presence_absence_uncertain` を別 class で出します。

これは optional heavy rescue lane であり、標準経路へ whole-genome alignment や OrthoFinder を要求しません。多数決を真実とせず、species disagreement・rearrangement・duplication・possible PAV を review-required に保持し、`rank_contribution=none` のため causal score を動かしません。`comparative-locus-benchmark` は frozen reconstruction にだけ degraded-reference truth を照合し、回復率・false expected inflation・candidate set size・runtime/memory を評価します。詳細は COMPARATIVE-LOCUS-RECONSTRUCTION.ja.md。

## 候補を安全に落とす（elimination layer）

`candidate-elimination` は、候補と機構仮説の**生存状態**を管理します。scoreを引くのではなく「何がどの証拠で落ちたか」を残すことで、真の原因候補を捨てずにcandidate setを圧縮します。

中核はstrengthと許される遷移の**構造的な束縛**です。registry読込時に検証され、違反する規則はエラーになります。

| strength | 許される遷移 | 例 |
| --- | --- | --- |
| `hard` | **gene-levelの除外のみ**（位置・遺伝証拠） | 情報のある組換え体で causal block 外、座標で区間外、明確な非共分離 |
| `strong` | **機構1つの棄却/弱化のみ**（gene全体は落とさない） | 規定coverageでfull CDS+splice境界を読んで親間差なし → `coding_rejected` |
| `moderate` | **weakenのみ** | 対象組織で他copyより明確に不整合なexpression、alternate paralogが強く支持 |
| `none` | `unresolved` / `technical_failure` / `coverage_incomplete` のみ | assay failure、PCR非増幅、部分sequencing、expression未測定、annotation無し、provider failure |

**absence of evidence を evidence of absence にしません。** coding negative は coding 機構だけを棄却し、gene は regulatory / SV 分岐を保ったまま生存します。quality/coverage gateを満たさない証拠は、registryの部分coverage版ruleへ降格して**何も除外しません**。

state変更はすべて rule id と evidence id 付きで `elimination_events.tsv` へ**追記のみ**で残ります。manifestに `candidate_reduction_ratio` / `hard_excluded_count` / `mechanism_only_rejected_count` / `unresolved_count` / `false_exclusion_risk_flags` と、`next-experiment` 向けの「残候補を最も多く識別できる実験を選ぶ」ヒントを出します。

rule registry は版管理された `templates/elimination_rules.tsv`（wheel同梱・drift test付き）で、`--rules` で差し替えられます。

## family から copy へ

`family-to-copy` は、familyまで絞れた後の「**どのcopyか**」を扱います。同一family内のmemberを orthology class / subfamily / synteny / domain architecture / expression / position の各軸で互いに比較し、軸ごとに独立列として保持します（合算しません）。

**homology gate:** family仮説の `gene_family` は探索用naming bucket（symbol stem等）のことがあり、**bucket単独ではcopy集合を定義しません**。memberが homology層（同一orthogroup/HOG・同一curated family・同一domain architecture+provenance）を共有する場合だけ通常のcopy discriminationを実行し（`copy_set_status=validated`）、bucketのみの場合は `family_resolution_required` として単独 `preferred_copy` を作りません。bucket内で別のhomologous familyに属するmemberは `not_in_validated_copy_set`（`not_applicable`）としてparalog比較に混ぜません。`family_copy_discrimination.tsv` に `copy_set_id` / `copy_set_status` 列を追加し、manifestの `homology_gate` に根拠と該当familyを記録します。

**target-relative copy selection:** copy 選択は **verified target specification**（`--family-target` TSV）への一致で行います。`subfamily_vs_target` / `domain_vs_target` は target の expected 値と copy を比較し（`target_match` / `target_mismatch`）、family 内 modal 比較（`domain_architecture_match` / `domain_vs_family_modal` / `subfamily_family_modal`）は**記述統計と outlier / truncation QC 専用**で選択には使いません。target が無い場合は modal だけで単独 preferred を作らず（`copy_target_missing` で target 解決 / orthology / synteny / position を要求）、recombinant core 等の locus-specific 独立証拠がある場合だけ例外として preferred を作ります。`pending_review`（未検証 / LLM proposal）target は単独 preferred の根拠にせず（`copy_target_pending_review`）、複数 copy が target に一致する duplicated target も target match だけで sole preferred にしません（`copy_target_duplicated`）。出力に `target_id` / `target_verification_status` / `target_reference_gene` / `target_subfamily` / `target_domain_architecture` / `subfamily_vs_target` / `domain_vs_target` / `domain_vs_family_modal` / `target_specific_axes` / `selection_basis`（`target_relative` / `locus_specific_only` / `unresolved_no_target`）を追加しました。

`preferred_copy` は **first Wet / first sequencing の優先順位**であって causal claim ではありません。BLAST best hit だけでは one-to-one ortholog とせず、複数orthology providerの一致は robustness であって独立票ではなく、同一evidence chainから派生した軸は二重計上しません。expression は copy間で差があるときだけ tie-break に使い、`unmeasured` を非発現へ変換しません。broad interval membership は copy間を区別しないので加点しません。**区別できる軸が1つも無ければ単独1位を作らず**、`tandem_ambiguous` / `paralog_ambiguous` のまま保持して追加の識別証拠を要求します。subfamily は統制map か orthogroup から導出し、どちらも無ければ `subfamily_unresolved` と正直に出します（symbol stem からは作りません）。

## Phenotype → pathway proposal → 検証済み bait panel

`llm-pathway-prompt` → `llm-pathway-import` → `source-verify` → `hypothesis-verify` → `bait-panel-freeze` → `query-eligibility` の順で、phenotypeからの pathway / gene-family 仮説を**検証済みbaitだけ**に絞ってから既存のmapping・rankへ渡します。canonical Route Aでは、機構を支持する検証済みPMID/PMCID/DOIと、解決済みaccession・種・配列を別々に検証し、両方のprovider/retrieval/content fingerprintを同じbaitにfreezeします。accession単独はpublication supportにならず、provider障害はnegativeではなく`review_required`です。LLM proposalはcanonical evidenceにならず、gene scoreにもTop判定にも一切影響しません。未検証accession・種不一致・取得不能配列はpanelに入らず、理由付きで rejected / review ledgerに残ります。blind modeでは `--candidates` を渡すとQTL区間候補IDのprompt混入を監査し、混入時はhard errorで停止します。

提案が揃ったら `family-index`（リファレンスごとに1回・キャッシュ）と `pathway-interval-falsify`（case毎・表引きのみ）で、**区間そのものを使って仮説を落とします**。粗い家系バケットで区間内メンバー0の仮説を先に切り、生き残ったものだけ精密なortholog階層へ進むカスケードなので、普通のラップトップで完結します。反証されるのは区間内メンバーが0のときだけで、存在するが濃縮していない家系は生存扱いです（recall優先）。濃縮の超幾何は**仮説フィルタであって遺伝子スコアではありません**。

triple-route では `route-a-interval-match` がこの反証を**検証済みbaitだけ**に適用します。`bait-panel-freeze` で凍結されたverified/source-backed targetのみを `pathway-interval-falsify` へ渡し、`reference_genome` scopeのfamily indexで区間内の**具体的なgene ID**へforward matchします。verified bait 0件・index未指定・区間未解決は `not_evaluated_*` として完走し、annotationにmemberが無いだけの不在は rescue 完了まで `absent_annotation_only` のholdです（rescue providerが失敗/部分のままのno-hitも同様）。出力の `route_a_interval_family_match.tsv` は hypothesis / family / verdict / family_status / match_level / interval_gene_ids / rescue status / source_ids / `rank_contribution=none` を持ち、verdictは `phenotype-first-status` へ渡って member_interval_genes・anchor_status となり、compare → family-to-copy のhandoffへ届きます。このstageはrankを変更しません。

段階別の評価は `phenotype-pathway-benchmark` で行い、既存58-caseのmap-based laneとは別laneとして集計します。詳細は [PHENOTYPE-PATHWAY-PROPOSAL.ja.md](PHENOTYPE-PATHWAY-PROPOSAL.ja.md) を参照してください。

`source-verify` のCLI既定は `initial_bibliography_v1` です。cache-firstで文献ID・タイトルの書誌情報を照合し、`--allow-external-api` を付けた場合は Europe PMC と公式NCBI PubMed/PMC E-utilities fallbackへ照会できます。既定では抄録・本文を読み直して支持内容を分類しません。`citation_verified` と `claim_support_status` は別列で、provider未取得・障害は `not_evaluated` として扱います。

明示的な `--verification-policy content_v1` では従来の内容確認を行い、その論文が何によって主張を支持しているかを分類します。このstrict経路では、実在する論文が遺伝子名を挙げているだけなら `unsupported_mention`、cueが立たない場合は `review_required` とし、agent抽出の引用文が取得済み原文に実在することも検査します。分類器自体は `claim-verification-benchmark` で採点します。詳細は [CLAIM-VERIFICATION.ja.md](CLAIM-VERIFICATION.ja.md)。

## Promoter / UTR / SV review

`promoter-review` はcoding scopeを検査した後だけ開始する、非昇格型のreview routeです。strand/isoform別promoter・5'UTR・3'UTR、親間SNV/indel、SV/PAV/TE、repeat/gap/low-quality flagを共通schemaに正規化します。`--hypotheses locus_hypotheses.tsv` を与えるとobservation IDをtyped contextとしてlinkしますが、gene scoreやTop判定には加点しません。reference-only motif hitも常にcontextです。

feature表(親間variant / SV-PAV / motif)はchunk単位で読み、対象region の chrom・範囲外の行を読み込み時点で捨てます。交差判定はchrom別のstart-sorted区間索引で、全region×全featureのPythonループは行いません。数十万〜数百万行の genome-wide table でも常駐メモリは対象region周辺だけで、判定結果は全走査実装と同一です。読み取り行数・保持行数・region境界・索引方式は `promoter_review_manifest.json` の `feature_intersection` に残ります。

## Multigenome private variants

`multigenome-private` は焦点系統の正規化eventを、宣言した対照系統群に対して決定論的に分類します。

- `private_confirmed`: 全ての必須表現型対照でcallableかつ同一eventが不在。
- `private_possible`: 同一eventは観測されないが、1系統以上の必須対照がuncallable(または callability 入力が無い)。
- `shared_confirmed`: 1系統以上の必須対照に同一eventが存在(ALT alleleがgenotypeに現れる)。
- `not_evaluated`: 必須比較またはcallability evidenceが未実施・不完全。

安全境界:

- VCFにレコードが無いだけではcallable referenceと判定しません。genotype `0/0` かつcallableのときだけ「対照はreference」です。
- 投影失敗・部分投影をgene absenceに変換しません。完全なprojection runのno-hitだけが `unmapped`、それ以外は `not_evaluated` です。
- `private_possible` は候補除外に使わず、追加callability確認の `required_next_action` になります。
- private statusはscoreへ一切加算されません(全出力の `rank_contribution` は `none`)。`working_candidate_panel.tsv` の `co_primary` は複数保持できます。
- 異なるcallsetの件数は加算せず、`candidate_control_matrix.tsv` にevent×accessionのreconciliationを出します。
- 区間内reference gene数・any-hit/high-confidence 1:1 投影数・callable CDS bp / intended CDS bp は `projection_and_callable_denominators.tsv` に必ず出力され、閾値未満なら `annotation-rescue` またはSV/PAV laneが `required_next_action` になります。

event keyは多アレルをALTごとに分割し、共通suffix→共通prefixの順にtrimした最小表現です。repeat上のleft-shiftは行わないので、入力は事前に `bcftools norm -f <reference> -m -any` で正規化してください。header由来の宣言有無は manifest の `event_normalization.inputs_declaring_normalization` に残り、未宣言はreview項目になります。

座ごとに異なる表現型対照は case manifest の `private_variant_model.loci[].phenotype_controls` で指定できます。`case-run` では annotation の後・ranking の前の正式stageとして実行され、未実施時は `run_manifest.json` の `private_variant_status` と report に `not_evaluated` が明示されます。
