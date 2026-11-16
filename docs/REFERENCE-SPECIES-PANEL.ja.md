# reference-species panel robustness

## 問題

Route A / B / C はどれも同じ参照コーパスを読む。そのコーパスは生物学の中立な
サンプルではない — Arabidopsis は深く、豆科・単子葉はそれより薄く、系統特異的
あるいは最近重複した遺伝子には Arabidopsis のカウンターパートが無いことすらある。
強い結果に見えても、それが「原因遺伝子の生物学」ではなく「UniProt curation史に
おける Arabidopsis の厚さ」を利用しているだけかもしれない。この依存は、
reference-species panel を実際に振って何が変わるかを見ない限り見えない。

## `reference-panel-benchmark` の2つの役割

### 1. `--from-corpus` — panelの構築とfingerprint

```bash
candidate-qtl-bench reference-panel-benchmark \
  --from-corpus reference_corpus.tsv --corpus-version v1 \
  --panel arabidopsis_only=3702 \
  --panel legume_only=3880,3847 \
  --panel multi_species=3702,3880,3847,4530 \
  --out-dir <dir>
```

各 `--panel LABEL=taxon_id_or_species,...` は `taxon_id` / `reference_species`
でコーパスをdeterministicにフィルタし、`panel_<label>.tsv` として書き出して
SHA-256をとる。`panel_membership.tsv` に membership token・行数・種一覧・
`corpus_version`・sha256 が並ぶ。panel構成は他のこのリポジトリの入力と同じ厳格
さでfingerprintされる。

`panel_species_coverage.tsv` は各panelのコーパス部分集合に
`reference_sources.corpus_coverage` をそのままかけたもの（フィールド×種の充足
率）。Route B が「本当にどの種のannotationに依存しているか」が panel 単位で見え
る。

### 2. `--manifest` — 凍結artifactの採点

```bash
candidate-qtl-bench reference-panel-benchmark --manifest manifest.tsv --out-dir <dir>
```

manifestは1行 = 1(case_id, panel_label)。`route_a_tsv` / `route_b_tsv`
（family列）、`orthology_tsv`（gene_id/reference_gene_id/relationship/source_id）、
`reference_function_tsv`（reference_gene_id/function_description/mutant_phenotype/
pathway/source_id）を任意で渡す。何も再実行しない — 凍結artifactを読むだけ。

## 「行が増えた」は「証拠が増えた」ではない

`distinct_primary_source_count` は `source_id` /
`phenotype_source_ids` / `literature_source_ids` を casefoldで重複除去して数える。
同じ一次実験が複数DB・複数種annotationにミラーされても、それは1つの証拠。
`redundant_source_ratio = reference_function_row_count / distinct_primary_source_count`
が1を大きく超えるとき、そのpanelの追加行はほぼ重複コピーであり、独立証拠の増加
ではない。

実測例（テストfixture）: legume_only panelは1行1source（ratio 1.0）。
multi_species panelはそこへ「同じ UniProt:Medtr_X をlegume DB向けに再輸出した行」
を1つ足しただけで3行になった — 行数は3倍でも distinct source は1つしか増えて
いない（ratio 1.5）。「種を足したら証拠が3倍になった」と誤読させない数字がこれ。

## 失敗の層別化

`failure_stratum` は4状態 + 「失敗なし」+ `not_evaluated`:

| stratum | 意味 |
|---|---|
| `missing_ortholog` | このpanelでは候補となるreference geneが見つからない |
| `paralog_ambiguity` | reference geneが複数（またはrelationship列がambiguousと明記） |
| `absent_phenotype_function` | orthologは見つかったが、そのcurated recordに function/phenotype/pathway が無い（`knowledge_unavailable`） |
| `annotation_gap` | 候補側（focal gene）自身の `focal_gene_annotation` が空 |
| （空文字） | このpanelで十分に被覆された |
| `not_evaluated` | orthology_tsv も reference_function_tsv も渡されていない |

`exact_gene_recall` と `failure_stratum` は別軸: テストfixtureの `thin_case` は
正しい遺伝子を見つけている（`exact_gene_recall=True`）が、そのreference record
に生物学的記述が無い（`absent_phenotype_function`）。「正しい遺伝子を見つけた」
と「その遺伝子について十分に語れる」は別の主張。

per-panel・per-stratum の recall 再計算は独自実装せず、
`knowledge_coverage.benchmark_stratify` をそのまま呼ぶ
（`stratified/benchmark_strata_summary.tsv`）。

## 唯一の最良panelは宣言しない

`panel_case_winners.tsv` はcaseごとに「実際に truth を回収できたpanel」を記録す
る。legume caseではlegume_onlyとmulti_speciesが勝ち、無関係なcaseでは同じpanel
が負けることがある — その依存が測りたいものであり、平均して消してよいノイズで
はない。manifestにも
「この benchmark は単一の最良 reference-species panel を宣言しない」と明文で
書く。

## 非目標

- panel間のスコアを比較しない（`reference-corpus-compare` と同じ原則）。
- 行数を独立証拠の量として扱わない。
- 「多様な種を足せば必ず良くなる」という結論を先取りしない — 減らないことも、
  冗長なだけのこともある。

関連: REFERENCE-CORPUS.ja.md /
NOVELTY-STRESS-BENCHMARK.ja.md
