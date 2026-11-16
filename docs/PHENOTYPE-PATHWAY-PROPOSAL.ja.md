# Phenotype → pathway 仮説 → 検証済み bait panel

本workflowは、**正当な非焦点種のpathwayホモログをbaitとして与えた後**のマップベース優先順位付けには強く、phenotypeだけから「どのpathway / gene familyを調べるか」を自律導出する力は弱い（ベンチマーク所見 の自律床 0/58）。

このrouteは、その前段に**検証可能なproposal層**を足す。既存の `homolog-panel-auto`、`query-eligibility`、`blind-benchmark`、deterministic rank は置き換えない。

```text
phenotype + direction + inheritance + organ/stage + focal species
  -> llm-pathway-prompt          versioned prompt/input package
  -> (llm-run など任意のprovider)  raw output
  -> llm-pathway-import          strict schema へ正規化(raw は verbatim 保存)
  -> hypothesis-verify           accession / species / source / sequence 検証
  -> bait-panel-freeze           immutable panel + FASTA + query manifest
  -> pathway-interval-falsify    区間で仮説を反証・濃縮判定(coarse -> fine カスケード)
  -> query-eligibility -> 区間内homology/orthology mapping -> 既存 deterministic rank
```

`family-index` はこの手前に1回だけ置く前処理で、リファレンスごとにキャッシュされます。

## 安全境界（実装とtestで固定）

- LLM proposalは canonical evidence にならない。このroute はどのファイルにも `evidence_code` / `dry_score` / `top_eligibility` を書かない。
- model confidence は proposal 属性として保存するだけで、gene score には一切加算しない。frozen panel には `confidence` 列自体を持たせない。
- 未検証accession・model が発明したaccessionはbait panelに入らない。
- model memory は citable source ではない。accession・種・配列を authoritative database から**再取得**し、fingerprintを保存できた提案だけがbaitになる。
- blind modeはQTL区間候補を見ない。`--candidates` を渡すと、その gene_id が prompt に混入していないかを監査し、混入時はhard errorで停止する（候補表はprompt本文には決して入らない）。監査結果は `llm_input_manifest.json` の `blind_leak_audit` に残る。
- 検証失敗した提案は削除せず、理由付きで `proposal_rejected.tsv` / `review_required.tsv` に残す。
- 単一pathway・単一baitは `review_required`。自動承認しない。

## コマンド

### `llm-pathway-prompt`

`prompt.txt`、`llm_pathway_prompts.jsonl`（既存 `llm-run` と同じJSONL契約）、`llm_input_manifest.json` を出す。manifestには case、focal species、phenotype、direction、inheritance、organ、stage、researcher constraints、known exclusions、execution mode、prompt template version と sha256、blind leak監査結果、`qtl_candidate_information_in_model_input` が入る。

- `--mode blind`: phenotype・direction/inheritance・organ/stage・focal species のみ。researcher contextは拒否する。
- `--mode assisted`: `--assisted-context` で研究者仮説・検証済み文献・ontology termを追加できる。manifest上でblindと明確に分離される。

### `llm-pathway-import`

raw JSON / JSONL / テキスト（コードフェンス・前置き付きも可）を strict schema へ正規化する。rawはverbatimコピーとsha256で残る。LLMを使わず手書きした提案も同じschemaでimportできる。

出力 `llm_pathway_hypotheses.tsv` の各行は `hypothesis_id` / `pathway_or_process` / `gene_family` / `reference_gene_name` / `reference_species` / `proposed_stable_id` / `expected_role` / `mechanistic_link` / `counterevidence` / `confidence` / `uncertainty` / `proposal_source` / `model` / `prompt_version` / `parser_version` / `execution_mode` / `review_status=pending_verification` を持つ。

parseできなければ提案を創作せず、`no_parsable_proposal` としてabstainする。distinct pathwayが2未満なら `single_pathway_overconfidence` を立てる。

### `hypothesis-verify`

提案ごとに次を検証する。

1. stable accessionが authoritative database に実在する
2. 取得したrecordの species が申告 species と一致する（属レベル）
3. 配列を再取得でき、sha256 fingerprintを保存できる
4. reference gene が focal species 由来でない
5. `--benchmark-mode` では focal-lineage / answer-backmapped marker を追加で除外する（`query_eligibility` と同じmarker集合）
6. 機能が source-backed（PMID または reviewed entry）である
7. 同一familyのparalog / 代替baitを `paralog_alternatives` に保持する

`--record-cache` にUniProt形状のJSONを渡せばオフラインで検証できる。`--allow-external-api` を付けたときだけ UniProt REST を叩く。**オフラインでcacheに無いaccessionは `verified` にならず `review_required`**（`verification_not_performed_offline`）になり、panelへは入らない。

### `bait-panel-freeze`

`verified` の提案だけを凍結し、`verified_bait_panel.tsv` / `verified_bait_panel.fasta` / `bait_panel_manifest.json` と、**`query-eligibility` がそのまま読める `verified_bait_panel_query_manifest.tsv`** を出す。selection_basis / sequence_origin は `schemas` の統制語彙（`pathway_member` / `non_pisum_reference_proteome`）を再利用する。保存した配列がfingerprintと一致しなければ凍結せずreviewへ回す。

検証済みbaitが1件だけの場合、panel全体はmechanism coverage不足として
`panel_status=review_required_small_panel` のままにする。一方、accession・source・species・sequenceを
通過したbait行の `review_status=verified` は保持し、非加点の決定論的interval projectionには使用する。
panelの網羅性不足とbait identityの未検証を同一状態に潰さず、copy/paralog曖昧性は後段で保持する。

### `phenotype-pathway-benchmark`

最終gene rankだけでなく段階別に測る。config TSV に case ごとの `expected_pathway` / `expected_gene_family` / `causal_gene_id` と各段階の成果物パスを与える。

- pathway Top-1 / Top-3 / Top-5 recall
- causal gene-family recall
- valid な非焦点種baitを1件以上作れたcase率
- stable accession実在率 / accession-species整合率 / verification pass率
- verified baitからcausal homologをTop-20へ回収した率
- wrong single-pathway overconfidence率、abstention率
- failure reason内訳: `pathway_miss` / `family_miss` / `retrieval_miss` / `accession_hallucination` / `verification_miss` / `panel_miss` / `mapping_or_ranking_miss` / `abstention_or_unparsable_output`

manifestの `lane` は `phenotype_only_proposal_lane` で、既存58-caseのmap-based homolog-informed laneとは**混ぜて集計しない**ことを明記する。

## 区間による仮説の反証（`family-index` / `pathway-interval-falsify`）

提案を集めた後、**区間そのものを使って仮説を落とす**段です。blindの線引きは次のとおり:

- **提案生成はblind**（区間を見せたら、見た遺伝子を選んでいるだけになる）
- **提案の反証・順位付けに区間を使うのは正当**なマップベースの手順

### 計算予算

このworkflowの核はリファレンス + Arabidopsis/近縁種 + QTL区間であり、**普通のラップトップで完結すること**が要件です。重い処理は1箇所に押し込んであります。

| 段 | 実行頻度 | 中身 |
| --- | --- | --- |
| `family-index` | **リファレンスごとに1回** | 全遺伝子 → モデル/近縁種プロテオームのhit表をfamily割り当てへ。fingerprint付きでキャッシュ |
| `pathway-interval-falsify` coarse | case毎 | 広い家系バケットの**表引きだけ**。配列検索なし |
| `pathway-interval-falsify` fine | coarse生存分のみ | best-hit reference gene 単位に精密化 |

両方とも manifest の `compute` に **実時間と peak memory** を記録するので、「軽い」は主張ではなく検証可能な数字です。

`family-index` は BLAST/DIAMOND の outfmt 6 表を入力に取ります。検索ツールは運用側の選択に委ねる設計で、家系ラベルは `--family-map`（Pfam・orthogroup等）で差し替えられます。既定はgene symbolのstem（`CLV1`→`CLV`、`PIN8`→`PIN`。`AT1G75820` のようなlocus tagはstemしない）。

### 判定

| verdict | 意味 |
| --- | --- |
| `surviving_enriched` | 区間内にメンバーがあり、BH補正後も濃縮 |
| `surviving_uninformative` | 区間内にメンバーはあるが濃縮なし。**生存扱い**（recall優先） |
| `falsified_no_member_in_interval` | index内の解決は成立したうえで、配列検索でも区間内に無い。強い反証 |
| `not_evaluated_family_unresolved_in_index` | 宣言identifier（`gene_family` / `reference_gene_name` / `proposed_stable_id`）がindex全体でどこにも解決しない。index語彙の不足（検出失敗）の可能性があり、**反証ではなく評価不能（hold）**としてreview項目になる |
| `not_computed_coarse_family_unresolved_in_index` | coarseが上記holdになったためfine計算を省略。反証には数えない |
| `not_computed_coarse_falsified` | 粗い階層で落ちたので精密計算を省略（カスケード） |

> **family indexの語彙依存性**: 家系ラベルはbest-BLAST-symbol stem由来のため、文献標準名がそのままindex語彙に一致するとは限らない（例: 真のNFR1オーソログが`CERK/LYK`ステムで登録される）。payloadの宣言identifierがindex全体で0解決のとき、zero-hitは真の不在ではなくindex語彙ギャップの可能性が残るため、このworkflowは`not_evaluated_family_unresolved_in_index`として評価不能のまま保持し、生物学的不在主張（`falsified_*`）には変換しない。

### 安全境界

- **反証は `interval_members == 0` のときだけ**。存在するが濃縮していない家系は生き残ります。濃縮は加点であってゲートではない。
- **enrichmentは仮説フィルタであって遺伝子ランカーではない。** メンバー1個の家系では超幾何が退化し（m=1, k=1 で p = n/N、区間内の全遺伝子で同じ値）、遺伝子の判別に使えません。遺伝子の順位付けはorthology信頼度（`top2_bitscore_ratio` 等）が担い、全出力の `rank_contribution` は `none` です。
- **tandem arrayは保守的p値の前に畳みます。** 遺伝子家系は物理的に固まるので、連鎖が濃縮に見えるのを防ぐため、区間内メンバーをblockへ畳んだ数を観測値として使います（`conservative_p` ≧ `enrichment_p`）。
- **家系未割り当ての遺伝子は分母に残したまま盲点として報告します。** 家系の真のサイズは測定値以上にしかなりません。
- 候補が30件を超えると「人が全件見られる規模」を外れたとしてreview項目になります。到達目標は5件固定ではなく**人手でレビュー可能な件数**です。

各候補には `falsifier`（何が観測されればこの候補が落ちるか）が付きます。長い実験期間中に候補集合を見返すための欄です。

## 再現情報

`llm_input_manifest.json` に prompt template version / hash・外部送信範囲・blind/assisted・QTL候補混入の有無、`llm_pathway_import_manifest.json` に model・provider・parser/schema version・raw と parsed 双方、`hypothesis_verification_manifest.json` に verification provider・取得日時・sequence fingerprint・rejection内訳、`bait_panel_manifest.json` に panel/FASTA の sha256 と bait ごとの provenance が残る。同一入力を別 model / provider で走らせて段階別に比較できる。
