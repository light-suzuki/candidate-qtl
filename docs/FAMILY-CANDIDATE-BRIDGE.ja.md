# Family Candidate Bridge（家族候補ブリッジ）

## 目的とスコープ

`family-candidate-bridge` は、Route A（phenotype-first 仮説 + 検証済み bait + 区間 family match）、
Route B（orthology 表 + interval gene atlas）、Route C（frozen retrieval の区間 match）という
独立に生成された 3 経路の家族エビデンスを、**1 枚の報告用テーブル**（`family_candidate_bridge.tsv`）に
決定的に正規化するコマンドです。

- **T1 スコープは決定的正規化コアのみ**。BLAST の実行、canonical wiring への組み込み、
  benchmark フィクスチャは後続スライスであり、本モジュールには含まれない。
- ブリッジは**並列の報告アーティファクト**であり、`family_to_copy.py` / `narrowing.py` の配線は
  変更しない。ランキングには一切供給されない（全行 `rank_contribution=none`、スコア列なし）。
- 出力は `family_candidate_bridge.tsv` と `family_candidate_bridge_manifest.json`。

## 入力（Route A / B / C）

| フラグ | アーティファクト | 必須列 | 備考 |
|---|---|---|---|
| `--route-a-hypotheses` | `phenotype_first_hypotheses.tsv` | `hypothesis_id`, `member_interval_genes` | `gene_family` / `proposed_stable_id` は任意 |
| `--route-a-family-status` | `route_a_interval_family_match.tsv` | `hypothesis_id` | `family_status` / `interval_gene_ids` / `source_ids` |
| `--route-a-bait-panel` | `verified_bait_panel.tsv` | `hypothesis_id` | `function_source_ids` / `source_provider` |
| `--route-b-orthology` | orthology 表 | `gene_id`, `reference_gene_id` | `relationship` / `orthogroup` / `percent_identity` / `query_coverage` / `evalue` / `provider` / `source_id`(pass-through) |
| `--route-b-atlas` | `interval_gene_atlas.tsv` | `gene_id` | `evolutionary_family_id` / `orthogroup` / `domain_*` / `family_id` / `family_label` / `focal_paralog_count` / `tandem_expansion` |
| `--route-c-retrieval` | `route_c_retrieval_merged.tsv` | `reference_gene_id`, `rank` | `database_version` / `model_version` |
| `--route-c-match` | `route_c_interval_match.tsv` | `reference_gene_id`, `rank`, `matched_interval_gene` | `match_basis` / `family_id` / `orthogroup_id` |
| `--dependency-clusters` | `evidence_dependency_clusters.tsv` | `dependency_cluster`, `source_ids` | `dependency_status` / `root_source_id` |
| `--genes` | 1 列（`gene_id`）TSV | `gene_id` | ブリッジ行を持たない遺伝子をマニフェストに報告 |

少なくとも 1 つの Route 入力が必須（無ければ非ゼロ終了で明確なエラー）。

### Route A の行生成

- `concrete_interval_gene_matched`（bait 検証済み）→ `curated`。マッチ遺伝子が 1 つなら
  `resolved_single_copy` / tier1、2 つ以上なら `paralog_ambiguous` / tier2。
- `family_present_copy_unresolved`（bait 検証済み）→ `curated` / `paralog_ambiguous` / tier2。
- `surviving_family_present_in_interval`：bait 検証済みなら `blastp`（検証済み ortholog 主張なしの
  家族存在）。bait 未検証なら tier6 + `not_evaluated_no_verified_bait`（`blastp` にはしない）。
- `surviving_sequence_only_annotation_gap` / `hold_absent_annotation_only` →
  `unclassified_gene` / tier6 / `not_evaluated_unclassified_gene`（annotation gap は provider 障害ではない
  ので provider は `complete`）。`interval_gene_ids` が無い場合は行を出さない。
- `falsified_*` の 2 状態 → **行を出さない**（偽証された家族はブリッジエビデンスではない）。
- `not_evaluated*` → tier6。reason は `not_evaluated_no_verified_bait` /
  `not_evaluated_missing_reference_genome_family_index` → `not_evaluated_provider_missing`、
  `not_evaluated_interval_unresolved` → `not_evaluated_interval_unresolved`。
- bait パネルに無い仮説（またはパネル未供給）→ `member_interval_genes` の各トークンに
  1 行 tier6（`not_evaluated_no_verified_bait`、provider `not_run`、ambiguity `not_evaluated`）。
  family status 行が完全に無い場合も同じ扱い。

### Route B の行生成

- orthology 各行（`gene_id` + `reference_gene_id` 必須）：
  `family_id` = `atlas.evolutionary_family_id` → `atlas.family_id` → `orthology.orthogroup` の
  最初の非空。全部空なら行をスキップし `skipped.route_b_rows_skipped_no_family` に数える。
- `match_type`：reciprocal クラス（`one_to_one` / `one_to_many` / `many_to_many` /
  `close_paralog` / `family_level_only` / `unmapped_ambiguous`）を持つ `relationship` は `rbh`。
  **`best_hit_homology` は `blastp`**（route_b_prepare が「local best hit は 1:1 orthology の証拠ではない」
  と明記するラベル。reciprocal でない best-hit 主張であり、orthogroup 列があっても行全体を
  `blastp` に格下げする）。`relationship` が無く `orthogroup` があれば `orthogroup`、どちらも無ければ
  `blastp`。
- `identity` ← `percent_identity`、`query_coverage` ← `query_coverage`、
  `target_coverage` ← orthology の `target_coverage` 列を pass-through（無ければ空）、
  `evalue` ← `evalue`、`paralog_count` ← `atlas.focal_paralog_count`（atlas 未供給なら空）。
- `provider_status` は常に `complete`、`provider_version` ← `provider`。
- `source_id`：orthology 入力に `source_id` 列があれば pass-through（無ければ空）。
- atlas 由来行：`--route-b-atlas` 供給時、遺伝子ごとにエビデンス軸ごと 1 行を追加で出力。
  `evolutionary_family_id` → `curated`、`domain_architecture`+`domain_ids`+`domain_source` → `domain`、
  `orthogroup` → `orthogroup`、どれも無く `family_id`/`family_label` のみ → `name_only`。
  atlas 行は provider `complete`、identity/coverage は空（pass-through のみ）。

### Route C の行生成

- match 行ごと：`gene_id` = `matched_interval_gene` の `;` トークン各々。
  `family_id` = `family_id` 無ければ `orthogroup_id`。両方無い行はスキップし
  `skipped.route_c_rows_skipped_no_family` に数える。
- `match_type`：`match_basis` に `same_orthogroup` を含めば `orthogroup`、それ以外 `name_only`。
- `source_id` ← `merged.database_version`（`(reference_gene_id, rank)` で結合）、
  `provider_version` ← `merged.model_version`。
- `paralog_count` = トークン数。複数遺伝子 → `paralog_ambiguous`、1 遺伝子 → `resolved_single_copy`。
- identity / coverage / evalue は空（`retrieval_score` は identity ではない）。

## 出力列（22 列、順序固定）

| # | 列 | 説明 |
|---|---|---|
| 1 | `case_id` | ケース識別子（既定 `""`） |
| 2 | `gene_id` | 区間遺伝子 |
| 3 | `family_id` | 家族（A: gene_family/仮説ID、B: 上記優先順、C: family_id/orthogroup_id） |
| 4 | `route` | `a` / `b` / `c` |
| 5 | `reference_accession` | 参照アクセッション（A: proposed_stable_id、B: reference_gene_id、C: reference_gene_id） |
| 6 | `source_id` | 行のソース ID（B は pass-through、C は database_version） |
| 7 | `match_type` | 空 / `curated` / `rbh` / `orthogroup` / `domain` / `blastp` / `name_only` |
| 8 | `identity` | 空または文字列（`0` と空は区別される） |
| 9 | `query_coverage` | 空または文字列 |
| 10 | `target_coverage` | 空または文字列 |
| 11 | `evalue` | 空または文字列 |
| 12 | `paralog_count` | 空または文字列 |
| 13 | `ambiguity_status` | 下記 enum |
| 14 | `provider_status` | 下記 enum |
| 15 | `provider_version` | provider バージョン |
| 16 | `input_fingerprint` | ソース行の sha256（ソース列名ソート順でタブ結合） |
| 17 | `output_fingerprint` | ブリッジ行自身の sha256（`output_fingerprint` 以外 21 列を列名ソートでタブ結合） |
| 18 | `review_status` | 下記 enum |
| 19 | `dependency_cluster_id` | `dep_<sha256(root)[:20]>` または `dependency_unknown` |
| 20 | `mapping_tier` | 下記 enum |
| 21 | `rank_contribution` | 常に `none` |
| 22 | `reason` | 空または下記 not_evaluated reason |

## Enum 語彙（`schemas.py` に定義、ブリッジから re-export）

- `match_type`：`curated` / `rbh` / `orthogroup` / `domain` / `blastp` / `name_only`（空可）
- `mapping_tier`（強い→弱い、タプル）：
  1. `tier1_source_backed_unique`
  2. `tier2_source_backed_ambiguous`
  3. `tier3_sequence_match_unique`
  4. `tier4_sequence_match_ambiguous`
  5. `tier5_name_only`
  6. `tier6_not_evaluated`
- `ambiguity_status`：`resolved_single_copy` / `paralog_ambiguous` / `tandem_ambiguous` /
  `coverage_insufficient` / `unclassified_gene` / `not_evaluated`
- `provider_status`：`complete` / `failed` / `unavailable` / `imported_partial` /
  `imported_complete` / `not_run`
- `review_status`：`not_required` / `review_required` / `review_complete` / `not_evaluated`
- `reason`（not_evaluated reason、空可）：
  `not_evaluated_provider_missing` / `not_evaluated_provider_failed` /
  `not_evaluated_provider_unavailable` / `not_evaluated_no_verified_bait` /
  `not_evaluated_no_reference_protein` / `not_evaluated_coverage_below_threshold` /
  `not_evaluated_interval_unresolved` / `not_evaluated_unclassified_gene`

スキーマは `FAMILY_CANDIDATE_BRIDGE_SCHEMA`（`schemas.validate` で検証可）。
`mapping_tier` はタプルでも `Series.isin` / pandera の `sorted()` の両方で動作する。
`case_id` / `gene_id` / `family_id` / `route` / `ambiguity_status` / `provider_status` /
`mapping_tier` / `rank_contribution` は非空必須。

## mapping tier の決定順（行ごとに厳密にこの順）

1. `provider_status` が `failed` / `unavailable` / `not_run`、または `reason` 非空 →
   `tier6_not_evaluated`。
2. `match_type == "name_only"` → `tier5_name_only`（他の列が何であってもこれ以上には昇格しない）。
3. `match_type` が `curated` / `rbh` / `orthogroup`（**Route C の orthogroup を除く**）：
   - `ambiguity_status == resolved_single_copy` かつ `provider_status` が `complete` /
     `imported_complete` かつ coverage/identity パス かつ `reason == ""` →
     `tier1_source_backed_unique`
   - それ以外 → `tier2_source_backed_ambiguous`
4. `match_type` が `domain` / `blastp`（**および Route C の `orthogroup`**）：
   - `resolved_single_copy` + provider `complete` / `imported_complete` + coverage/identity パス +
     `reason == ""` → `tier3_sequence_match_unique`
   - それ以外 → `tier4_sequence_match_ambiguous`
5. それ以外（`match_type` 空）→ `tier6_not_evaluated`。

> **Route C のキャップ**：frozen retrieval corpus の orthogroup キー一致は source-backed な
> orthology 主張ではないため、Route C の `orthogroup` 行は sequence-match 階層
> （tier3 / tier4）に留まる。`name_only` 行は常に tier5。

### Coverage / identity のパス判定

- 値が空 → パス（該当なし）。
- `identity`：値 < 1.0 は小数、値 >= 1.0 はパーセント（100 で除算）。`MIN_IDENTITY_FRACTION = 0.3`
  未満で失敗。
- `query_coverage` / `target_coverage`：小数（0–1）。`MIN_COVERAGE_FRACTION = 0.5` 未満で失敗。
  T1 では percent 表記の coverage は正規化しない（percent 正規化は `identity` のみ）。
- 数値に変換できない値は fail-closed（確認できない閾値は昇格させない）。

### 閾値未満の行の扱い（state table (e) の分割を実装）

- **source-backed 行**（`curated` / `rbh` / `orthogroup`）：エビデンスは保持し、
  `tier2_source_backed_ambiguous` に留める。`ambiguity_status = coverage_insufficient`、
  `reason = ""`（エビデンスは消さない）、`review_status = review_required`。
- **sequence-only 行**（`domain` / `blastp`）：`tier6_not_evaluated`、
  `ambiguity_status = coverage_insufficient`、`reason = not_evaluated_coverage_below_threshold`、
  `review_status = not_evaluated`。

## ambiguity の導出（tier とは独立）

1. provider が `not_run` / `failed` / `unavailable`、または `reason` 非空 → `not_evaluated`
2. coverage / identity 閾値未満 → `coverage_insufficient`
3. マッチ遺伝子 >= 2、または `paralog_count` >= 2、または family_status が
   `family_present_copy_unresolved`、または Route B の relationship が
   `one_to_many` / `many_to_many` / `close_paralog` / `unmapped_ambiguous` →
   `paralog_ambiguous`（atlas の `tandem_expansion` が truthy なら `tandem_ambiguous`）
4. それ以外 → `resolved_single_copy`

`unclassified_gene` は Route A の family_status が `surviving_sequence_only_annotation_gap` /
`hold_absent_annotation_only` かつ `interval_gene_ids` が存在する行にのみ設定される。

## dependency clustering

- `--dependency-clusters` 供給時：`source_ids`（`;` / `,` 区切り）と `root_source_id` から
  ソース → クラスタの写像を作る。行の `source_id` が写像に載っていればそのクラスタ。
  複数クラスタに載るソースは除外（`dependency_unknown`）、載らない場合も `dependency_unknown`。
- クラスタ表なし：
  - Route A：検証済み bait の `function_source_ids` が非空 1 値のみ → それを root に
    `dep_sha256(root)`。0 個 / 2 個以上 → `dependency_unknown`。
  - Route B：常に `dependency_unknown`。
  - Route C：常に `database_version` を root に `dep_sha256(database_version)`（宣言された root）。
- クラスタ ID は `dep_` + sha256(root).hexdigest()[:20]。
- **行の `source_id` は決して自動で root に昇格しない**（クラスタ表か検証済み bait の
  `function_source_ids` 経由のみ）。
- マニフェストの `dependency_clusters` はクラスタごとに 1 エントリ：
  `dependency_cluster_id` / `dependency_status`（`root_source_known` | `dependency_unknown`）/
  `root_source_id` / `routes`（`;` 結合ソート）/ `gene_count`（異なる遺伝子数）/
  `row_count`。`dependency_cluster_id` でソート。

## fingerprint と決定性の契約

- `input_fingerprint`：ソース行の列（ソースフレームの全列）を列名ソートでタブ結合した
  sha256。ソースファイルが同一なら同一。
- `output_fingerprint`：`output_fingerprint` 以外の 21 列を列名ソートでタブ結合した sha256。
  依存クラスタ解決後の最終行に一度だけ打刻される。
- 完全重複行は `drop_duplicates` で除去。ソートキーは
  `["case_id","gene_id","family_id","route","match_type","dependency_cluster_id","source_id"]`
  （stable sort）。TSV は `sep="\t"`、`lineterminator="\n"`。
- **絶対パスは TSV / マニフェストに一切現れない**：入力レコードは `path` を持たず、
  マニフェストの `outputs` は相対ファイル名のみ。
- 同一入力 + 異なる out-dir で TSV もマニフェストもバイト同一。

## マニフェスト

`schema_version: "1.0"`、`command: "family-candidate-bridge"`、`case_id`、
`inputs`（各入力名ごとに `present` / `sha256` / `size_bytes` / `mtime_ns` / `method` / `strength`）、
`row_counts`（`total` と `by_route` / `by_mapping_tier` / `by_ambiguity_status` /
`by_provider_status` / `by_review_status`）、`dependency_clusters`、
`skipped`（`route_b_rows_skipped_no_family` / `route_c_rows_skipped_no_family`）、
`rank_contribution: "none"`、`thresholds`、`outputs`、`genes_without_bridge_rows`（`--genes` 供給時のみ）。

## 現状の制約と後続スライス

- T1 は決定的正規化のみ。BLAST 実行、`family_to_copy` / `narrowing` への canonical 配線、
  benchmark 参加は未実装（後続スライス）。
- `family_to_copy.py` / `narrowing.py` / `workflow_extensions.py` の配線は変更していない。
- ブリッジは報告専用：スコア列・ランク列は存在せず、`rank_contribution` は全行 `none`。
